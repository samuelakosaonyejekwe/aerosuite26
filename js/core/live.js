// Live-data connectors. Three sources, in this order of independence:
//  1. data bundled with the app (the world airport and runway database) — no network at all;
//  2. requests straight from the user's own browser to public, CORS-enabled providers;
//  3. the cloud snapshot, data/snapshot.json, which a scheduled job (tools/snapshot.mjs) rebuilds every few
//     hours on the hosting side: it gives every location-independent feed an instant first value at start-up
//     and is the only route to providers that browsers cannot reach (no CORS): carbon price, jet-fuel spot
//     price and policy interest rates.
// Responses are cached in IndexedDB with a time-to-live and served stale-while-revalidate, so the app keeps
// working offline with the last known values.
//
// Which provider a feed may use is decided by the deployment configuration (config.json at the site root) and
// the audited licence registry (js/data/licences.json). With "commercial": true a source is contacted only
// when the registry classes it "commercial-ok", or the operator supplied the key or address it needs, or the
// operator listed its id under "accept". Every gate goes through allowed() below.

import { idb, state, patchCaseMany, emit, ls } from './store.js';
import { isa } from './atmosphere.js';
import * as airports from './airports.js';
import { gridWeather, gridMarine, gridAir } from './gridwx.js';
import { searchPlaces, nearestPlace, terrain, climateAt } from './places.js';

// ---- unit conversions shared with tools/snapshot.mjs --------------------------------------------
export const JET_DENSITY_KG_L = 0.804, L_PER_BBL = 158.987, L_PER_US_GAL = 3.785411784; // Jet A-1 at 15 °C, mid-specification
export const usdKgFromUsdGal = (usdGal) => usdGal / (L_PER_US_GAL * JET_DENSITY_KG_L);
export const usdKgFromUsdBbl = (usdBbl) => usdBbl / (L_PER_BBL * JET_DENSITY_KG_L);
/** EUR amount → USD, given the ECB reference rate expressed as EUR per 1 USD (the form the FX feed uses). */
export const usdFromEur = (eur, eurPerUsd) => eur / eurPerUsd;

const J = (r) => { if (!r.ok) throw new Error(`HTTP ${r.status}`); return r.json(); };
const T = (r) => { if (!r.ok) throw new Error(`HTTP ${r.status}`); return r.text(); };
const get = (url, opt) => fetch(url, { ...opt, cache: 'no-store', signal: AbortSignal.timeout ? AbortSignal.timeout(20000) : undefined });
/** Per-request network handle that notices when the service worker answered from its offline copy. */
function makeNet() {
  const n = { cachedAt: null, get: async (url, opt) => { const r = await get(url, opt), c = r.headers.get('x-sw-cached-at'); if (c) n.cachedAt = Math.min(n.cachedAt ?? Infinity, Number(c)); return r; } };
  return n;
}
const MIN = 60e3, HOUR = 3600e3, DAY = 86400e3;
const r3 = (v) => Math.round(v * 1000) / 1000;

const LEVELS = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150];
const pageBase = () => { try { return typeof document !== 'undefined' && /^https?:/.test(document.baseURI) ? document.baseURI : null; } catch { return null; } };
const isNode = typeof process !== 'undefined' && !!process.versions?.node && typeof window === 'undefined';

// ---- deployment configuration (config.json) -----------------------------------------------------
/** Defaults used when config.json is missing. Every key is documented in README.md ("Deployment configuration"). */
export const CONFIG_DEFAULTS = Object.freeze({ commercial: false, openMeteoApiKey: '', weatherProvider: 'auto', carbonPriceUrl: '', attribution: true, overpassUrl: '', openAlexApiKey: '', accept: [] });
export const config = { ...CONFIG_DEFAULTS, accept: [] };
const httpsUrl = (v) => { try { const u = new URL(String(v)); return u.protocol === 'https:' || u.hostname === 'localhost' ? u.href : ''; } catch { return ''; } };
/** Keep only known keys with the right types; anything else falls back to the default. */
export function cleanConfig(o) {
  const c = { ...CONFIG_DEFAULTS, accept: [] }; if (!o || typeof o !== 'object') return c;
  if (typeof o.commercial === 'boolean') c.commercial = o.commercial;
  for (const k of ['openMeteoApiKey', 'openAlexApiKey']) if (typeof o[k] === 'string') c[k] = o[k].trim();
  if (['auto', 'open-meteo', 'met-norway', 'nws', 'noaa-grid', 'none'].includes(o.weatherProvider)) c.weatherProvider = o.weatherProvider;
  if (typeof o.carbonPriceUrl === 'string' && o.carbonPriceUrl.trim()) { // an https address, or a path on this site
    const b = pageBase(); try { const u = new URL(o.carbonPriceUrl.trim(), b || undefined); c.carbonPriceUrl = httpsUrl(u.href) || (b && u.origin === new URL(b).origin ? u.href : ''); } catch { /* not an address */ }
  }
  if (typeof o.overpassUrl === 'string') c.overpassUrl = httpsUrl(o.overpassUrl.trim());
  if (Array.isArray(o.accept)) c.accept = o.accept.filter((x) => typeof x === 'string');
  c.attribution = c.commercial ? true : o.attribution !== false; // required notices cannot be switched off in a commercial deployment
  return c;
}
let configP = null;
/** Set the configuration directly (tests, the snapshot tool). */
export function setConfig(o) { Object.assign(config, cleanConfig(o)); configP = Promise.resolve(config); return config; }
/**
 * Read config.json once: this origin's copy, else the copy saved on this device, else the defaults. The saved
 * copy means a device that has seen "commercial": true stays in commercial mode while offline.
 */
export function loadConfig({ force = false } = {}) {
  if (configP && !force) return configP;
  configP = (async () => {
    let raw = globalThis.__AEROSUITE_CONFIG__ || null; const base = pageBase();
    if (!raw && base) { try { const r = await get(new URL('config.json', base).href); if (r.status === 404) raw = {}; else raw = await J(r); if (raw && typeof raw === 'object') ls.set('config', raw); } catch { raw = ls.get('config', null); } }
    Object.assign(config, cleanConfig(raw || {}));
    return config;
  })();
  return configP;
}

// ---- licence registry (js/data/licences.json) ---------------------------------------------------
let regP = null;
/** The audited registry of sources and components. Bundled with the app; never fetched from elsewhere. */
export function loadLicences() {
  regP ||= (async () => {
    if (globalThis.__AEROSUITE_LICENCES__) return globalThis.__AEROSUITE_LICENCES__;
    let url = null; try { if (import.meta.url) url = new URL('../data/licences.json', import.meta.url); } catch { /* bundled copy */ }
    if (!url && pageBase()) url = new URL('js/data/licences.json', pageBase());
    if (!url) throw new Error('The licence registry is not available in this copy of the app');
    if (isNode && url.protocol === 'file:') { const fs = 'node:fs/promises', { readFile } = await import(fs); return JSON.parse(await readFile(url, 'utf8')); }
    return fetch(url).then(J);
  })().catch((e) => { regP = null; throw e; });
  return regP;
}
/** An error raised because the deployment's licence policy rules a source out (not a network failure). */
const policy = (msg) => Object.assign(new Error(msg), { policy: true });
/**
 * May this deployment contact the registry source `id`? Always in a non-commercial deployment. In a commercial
 * one: only "commercial-ok" sources, sources whose key or address the operator configured, and ids the operator
 * accepted. An unknown id, or a registry that cannot be read, is refused.
 */
export async function allowed(id, cfg = config) {
  if (!cfg.commercial) return true;
  if (cfg.accept.includes(id)) return true;
  let reg; try { reg = await loadLicences(); } catch { return false; }
  const s = reg.sources.find((x) => x.id === id); if (!s) return false;
  if (s.usedInCommercial === false) return false; // cleared licence, but the deployment policy keeps it out of commercial mode (see the registry entry)
  if (s.class === 'commercial-ok') return true;
  if (id === 'open-meteo-customer') return !!cfg.openMeteoApiKey;
  if (id === 'operator-carbon') return !!cfg.carbonPriceUrl;
  return false;
}
/** Same question without the registry, for labels and cache keys: is a source that is not commercial-ok open here? */
const open = (id) => !config.commercial || config.accept.includes(id);
/** Cache-key prefix: values fetched under non-commercial terms are never reused by a commercial deployment. */
const modeKey = () => (config.commercial ? 'c:' : '');

// ---- Open-Meteo: free hosts (non-commercial only) or customer hosts with the operator's key -----
const omId = () => (config.openMeteoApiKey ? 'open-meteo-customer' : open('open-meteo-free') ? 'open-meteo-free' : null);
/** Address on the right Open-Meteo host. `sub` is '', 'archive-', 'air-quality-', 'marine-' or 'geocoding-'. */
async function omUrl(sub, path, query) {
  const key = config.openMeteoApiKey;
  if (key && (await allowed('open-meteo-customer'))) return `https://customer-${sub}api.open-meteo.com${path}?${query}&apikey=${encodeURIComponent(key)}`;
  if (await allowed('open-meteo-free')) return `https://${sub}api.open-meteo.com${path}?${query}`;
  throw policy('Open-Meteo allows commercial use only with a subscription key (openMeteoApiKey in config.json)');
}

// ---- weather providers: each returns the same normalised object ---------------------------------
// Fields: time (UTC, yyyy-mm-ddThh:mm), elev_m, T_C, rh (0–1), p_hPa (station), qnh_hPa, wind_ms, wind_dir_deg,
// gust_ms, precip_mm_h, cloud_pct, code, visibility_m, freezing_level_m, aloft[{ hPa, alt_m, T_C, speed_ms, dir_deg }].
// A provider that cannot supply a field leaves it null and names it in `missing`; finishWeather() then fills the
// documented stand-ins (station pressure from sea-level pressure, freezing level from the standard lapse rate)
// and writes a note that the pages show. Winds aloft are never invented: an empty list makes the mission
// analyses fall back to still air and the standard atmosphere.
//   Open-Meteo   — everything.
//   MET Norway   — no station pressure (derived), no visibility, no freezing level, no winds aloft; gusts only in the Nordic area.
//   US NWS       — latest station observation, United States only: no freezing level, no winds aloft; gusts only when reported.
const stationPressure = (qnh_hPa, elev_m, T_C) => qnh_hPa * Math.exp((-9.80665 * elev_m) / (287.05 * (T_C + 273.15 + 0.00325 * elev_m)));
const utcMinute = (t) => new Date(t).toISOString().slice(0, 16);
const lat4 = (v) => String(Math.round(v * 1e4) / 1e4);
export const WEATHER = {
  'open-meteo': {
    label: 'Open-Meteo', reg: () => omId(),
    async load(p, get) {
      const hv = ['freezing_level_height', 'visibility', ...LEVELS.flatMap((l) => [`temperature_${l}hPa`, `wind_speed_${l}hPa`, `wind_direction_${l}hPa`, `geopotential_height_${l}hPa`])].join(',');
      const d = await get(await omUrl('', '/v1/forecast', `latitude=${p.lat}&longitude=${p.lon}&current=temperature_2m,relative_humidity_2m,surface_pressure,pressure_msl,wind_speed_10m,wind_direction_10m,wind_gusts_10m,precipitation,cloud_cover,weather_code&hourly=${hv}&forecast_hours=1&wind_speed_unit=ms&timezone=UTC`)).then(J);
      const c = d.current, hr = d.hourly || {}, at = (k) => (hr[k] ? hr[k][0] : null);
      return {
        time: c.time, elev_m: d.elevation, T_C: c.temperature_2m, rh: c.relative_humidity_2m / 100, p_hPa: c.surface_pressure, qnh_hPa: c.pressure_msl, wind_ms: c.wind_speed_10m, wind_dir_deg: c.wind_direction_10m, gust_ms: c.wind_gusts_10m,
        precip_mm_h: c.precipitation, cloud_pct: c.cloud_cover, code: c.weather_code, visibility_m: at('visibility'), freezing_level_m: at('freezing_level_height'),
        aloft: LEVELS.map((l) => ({ hPa: l, alt_m: at(`geopotential_height_${l}hPa`), T_C: at(`temperature_${l}hPa`), speed_ms: at(`wind_speed_${l}hPa`), dir_deg: at(`wind_direction_${l}hPa`) })).filter((a) => a.alt_m != null && a.T_C != null),
      };
    },
  },
  'met-norway': {
    label: 'MET Norway', reg: () => 'met-norway',
    async load(p, get) {
      // MET Norway asks every client to identify itself; a page identifies itself by its Origin header, which a copy opened from disk does not have.
      if (typeof location !== 'undefined' && location.protocol === 'file:') throw policy('MET Norway needs the app to be served from a web address');
      const d = await get(`https://api.met.no/weatherapi/locationforecast/2.0/complete?lat=${lat4(p.lat)}&lon=${lat4(p.lon)}`).then(J); // at most 4 decimals (terms of service)
      const ts = d.properties?.timeseries || [], now = Date.now(), k = Math.max(0, ts.findIndex((x) => Date.parse(x.time) > now) - 1), e = ts[k < ts.length ? k : 0]; if (!e) throw new Error('empty forecast');
      const i = e.data.instant.details, n1 = e.data.next_1_hours?.details, n6 = e.data.next_6_hours?.details, num = (v) => (Number.isFinite(v) ? v : null);
      return {
        time: utcMinute(e.time), elev_m: num(d.geometry?.coordinates?.[2]), T_C: num(i.air_temperature), rh: Number.isFinite(i.relative_humidity) ? i.relative_humidity / 100 : null, p_hPa: null, qnh_hPa: num(i.air_pressure_at_sea_level),
        wind_ms: num(i.wind_speed), wind_dir_deg: num(i.wind_from_direction), gust_ms: num(i.wind_speed_of_gust), precip_mm_h: n1 ? num(n1.precipitation_amount) : n6 && Number.isFinite(n6.precipitation_amount) ? n6.precipitation_amount / 6 : 0,
        cloud_pct: num(i.cloud_area_fraction), code: null, visibility_m: null, freezing_level_m: null, aloft: [],
      };
    },
  },
  nws: {
    label: 'US National Weather Service', reg: () => 'nws',
    async load(p, get) {
      const pr = await get(`https://api.weather.gov/points/${lat4(p.lat)},${lat4(p.lon)}`); if (pr.status === 404) throw new Error('outside the area covered (United States only)');
      const pt = await J(pr), su = pt.properties?.observationStations; if (!su || !/^https:\/\/api\.weather\.gov\//.test(su)) throw new Error('no observing station listed');
      const st = (await get(`${su}?limit=1`).then(J)).features?.[0]?.properties?.stationIdentifier; if (!/^[A-Z0-9]{3,8}$/.test(st || '')) throw new Error('no observing station nearby');
      const o = (await get(`https://api.weather.gov/stations/${st}/observations/latest`).then(J)).properties, v = (x) => (Number.isFinite(x?.value) ? x.value : null);
      const ms = (x) => (v(x) == null ? null : /km_h/.test(x.unitCode || '') ? v(x) / 3.6 : v(x)), T = v(o.temperature); if (T == null) throw new Error(`station ${st} reports no temperature`);
      const okta = { SKC: 0, CLR: 0, FEW: 19, SCT: 44, BKN: 75, OVC: 100, VV: 100 }, cl = (o.cloudLayers || []).map((c) => okta[c.amount]).filter(Number.isFinite);
      const alt = v(o.barometricPressure), slp = v(o.seaLevelPressure);
      return {
        time: utcMinute(o.timestamp), elev_m: v(o.elevation), T_C: T, rh: v(o.relativeHumidity) == null ? null : v(o.relativeHumidity) / 100, p_hPa: null, qnh_hPa: alt != null ? alt / 100 : slp != null ? slp / 100 : null,
        wind_ms: ms(o.windSpeed), wind_dir_deg: v(o.windDirection), gust_ms: ms(o.windGust), precip_mm_h: v(o.precipitationLastHour) ?? 0, cloud_pct: cl.length ? Math.max(...cl) : null, code: null,
        visibility_m: v(o.visibility), freezing_level_m: null, aloft: [], station: st,
      };
    },
  },
};
WEATHER['noaa-grid'] = {
  label: 'NOAA GFS forecast grid', reg: () => 'noaa-gfs',
  note: 'Interpolated from the NOAA GFS forecast grid (1° at the surface, 2.5° aloft): regional model values, not a station observation.',
  async load(p) { return gridWeather(p.lat, p.lon, { elev_m: Number.isFinite(p.elev_m) ? p.elev_m : await bundledElevation(p) }); },
};
/** Site elevation from the data bundled with the app: an airport within 10 km, else a city within 15 km, else the 0.5° terrain grid. */
export async function bundledElevation(p) {
  const a = (await airports.nearest(p.lat, p.lon, 10, 1).catch(() => []))[0]; if (Number.isFinite(a?.elev_m)) return a.elev_m;
  const c = await nearestPlace(p.lat, p.lon, 15).catch(() => null); if (Number.isFinite(c?.elev_m)) return c.elev_m;
  return (await terrain([p]).catch(() => [null]))[0];
}
/**
 * Complete a provider's reply from the NOAA GFS grid: whatever it left empty (gusts, visibility, freezing level,
 * winds and temperatures aloft, …) is taken from the grid, so no deployment loses a field because of its provider.
 */
async function gridFill(w, p) {
  const gaps = ['gust_ms', 'visibility_m', 'freezing_level_m', 'rh', 'wind_ms', 'cloud_pct', 'qnh_hPa'].some((k) => !Number.isFinite(w[k])) || !w.aloft?.length;
  if (!gaps || !(await allowed('noaa-gfs'))) return w;
  let g; try { g = await gridWeather(p.lat, p.lon, { elev_m: Number.isFinite(p.elev_m) ? p.elev_m : Number.isFinite(w.elev_m) ? w.elev_m : null }); } catch { return w; }
  const out = { ...w }, filled = [];
  for (const [k, label] of [['gust_ms', 'gusts'], ['visibility_m', 'visibility'], ['freezing_level_m', 'freezing level'], ['rh', 'humidity'], ['wind_ms', 'wind'], ['cloud_pct', 'cloud cover'], ['qnh_hPa', 'pressure']]) if (!Number.isFinite(out[k]) && Number.isFinite(g[k])) { out[k] = g[k]; if (k === 'wind_ms') out.wind_dir_deg = g.wind_dir_deg; filled.push(label); }
  if (filled.includes('gusts') && Number.isFinite(out.wind_ms)) out.gust_ms = Math.max(out.gust_ms, out.wind_ms);
  if (!out.aloft?.length && g.aloft.length) { out.aloft = g.aloft; filled.push('winds and temperatures aloft'); }
  if (out.code == null) out.code = g.code;
  return filled.length ? { ...out, filled, fill_cycle: g.cycle } : w;
}
/** Providers to try, in order, under the current configuration. */
export function weatherPlan(cfg = config) {
  const w = cfg.weatherProvider, om = !cfg.commercial || !!cfg.openMeteoApiKey || cfg.accept.includes('open-meteo-free');
  if (w === 'none') return [];
  // Open-Meteo (finest) where its terms allow; then a point forecast or observation completed from the NOAA grid; the grid alone
  // is the last step and also works from the offline copy.
  const auto = cfg.commercial ? [...(om ? ['open-meteo'] : []), 'met-norway', 'nws', 'noaa-grid'] : ['open-meteo', 'met-norway', 'noaa-grid'];
  if (w === 'auto' || !WEATHER[w] || (w === 'open-meteo' && !om)) return auto;
  return [w, ...auto.filter((x) => x !== w)];
}
/** Fill the documented stand-ins for fields a provider cannot supply and say so. */
export function finishWeather(w, id) {
  const P = WEATHER[id], missing = [], notes = [], out = { ...w, provider: id, provider_label: P.label + (w.station ? ` (station ${w.station})` : '') + (w.filled?.length ? ' + NOAA GFS grid' : '') };
  if (!Number.isFinite(out.T_C)) throw new Error('no temperature in the reply');
  if (!Number.isFinite(out.elev_m)) out.elev_m = null;
  if (!Number.isFinite(out.qnh_hPa) && !Number.isFinite(out.p_hPa)) { missing.push('pressure'); out.qnh_hPa = 1013.25; notes.push('no pressure reported: standard sea-level pressure (1013.25 hPa) assumed'); }
  if (!Number.isFinite(out.p_hPa)) { out.p_hPa = Math.round(stationPressure(out.qnh_hPa, out.elev_m || 0, out.T_C) * 10) / 10; if (!missing.includes('pressure')) notes.push('station pressure derived from sea-level pressure and elevation'); }
  if (!Number.isFinite(out.rh)) { missing.push('humidity'); out.rh = null; }
  if (!Number.isFinite(out.wind_ms)) { missing.push('wind'); out.wind_ms = null; out.wind_dir_deg = null; }
  if (!Number.isFinite(out.gust_ms)) { missing.push('gusts'); out.gust_ms = null; notes.push('no gust data: gust taken equal to the mean wind'); }
  if (!Number.isFinite(out.visibility_m)) { missing.push('visibility'); out.visibility_m = null; }
  if (!Number.isFinite(out.freezing_level_m)) { missing.push('freezing level'); out.freezing_level_m = Math.round((out.elev_m || 0) + Math.max(0, out.T_C) / 0.0065); out.freezing_level_estimated = true; notes.push('freezing level estimated from the surface temperature with the standard lapse rate (6.5 K/km)'); }
  if (!out.aloft?.length) { missing.push('winds aloft'); out.aloft = []; notes.push('no winds or temperatures aloft: mission analyses use still air and the standard atmosphere (ISA)'); }
  out.missing = missing;
  out.note = [w.filled?.length ? `${w.filled.join(', ').replace(/^./, (c) => c.toUpperCase())}: NOAA GFS forecast grid (1° at the surface, 2.5° aloft), which ${P.label} does not supply.` : '', P.note || '',
    notes.length ? `${missing.length ? `${P.label} does not supply ${missing.join(', ')}. ` : ''}${notes.map((n) => n[0].toUpperCase() + n.slice(1)).join('. ')}.` : ''].filter(Boolean).join(' ');
  return out;
}

// ---- carbon price sources usable in a commercial deployment ------------------------------------
const GOVUK = 'https://www.gov.uk', UK_ETS_PAGE = '/government/publications/determinations-of-the-uk-ets-carbon-price';
/**
 * UK Emissions Trading Scheme carbon price determined by the UK ETS Authority for the current scheme year
 * (GOV.UK content API, Open Government Licence v3.0). A statutory reference figure fixed once a year in GBP
 * per tonne CO2e — it is neither an EU ETS price nor a daily market price, and is labelled accordingly.
 */
export async function ukEtsCarbon(get) {
  const idx = await get(`${GOVUK}/api/content${UK_ETS_PAGE}`).then(J);
  const att = (idx.details?.attachments || []).map((a) => ({ url: a.url, year: Math.max(...(String(a.title).match(/20\d\d/g) || []).map(Number)) })).filter((a) => /^\/government\/publications\//.test(a.url || '') && Number.isFinite(a.year)).sort((a, b) => b.year - a.year)[0];
  if (!att) throw new Error('no determination listed');
  const doc = await get(`${GOVUK}/api/content${att.url}`).then(J), text = String(doc.details?.body || '').replace(/<[^>]+>/g, ' ').replace(/&pound;|&#163;/g, '£').replace(/&nbsp;|&#160;/g, ' ').replace(/\s+/g, ' ');
  const m = /scheme year beginning on 1 January (20\d\d) is £\s?(\d+(?:\.\d+)?)/.exec(text); if (!m) throw new Error('determination text not recognised');
  const year = Number(m[1]), price = Number(m[2]); if (!(price > 1 && price < 1000)) throw new Error('no plausible price');
  const pub = Date.parse(doc.first_published_at || doc.details?.first_public_at || doc.public_updated_at);
  return { date: Number.isFinite(pub) ? new Date(pub).toISOString().slice(0, 10) : `${year - 1}-12-01`, price, currency: 'GBP', gbp_t: price, scheme_year: year, valid_from: `${year}-01-01`, valid_to: `${year}-12-31`, market: 'UK ETS',
    instrument: `UK Emissions Trading Scheme carbon price determined by the UK ETS Authority for scheme year ${year}, GBP per tonne CO₂e (annual statutory figure; not an EU ETS price and not a daily market price)`,
    source: 'UK ETS Authority, GOV.UK', attribution: 'Contains public sector information licensed under the Open Government Licence v3.0.', url: GOVUK + att.url, registry: ['uk-ets-price'] };
}
/** Carbon price from the operator's own JSON: { "date": "2026-10-08", "price": 85.07, "currency": "EUR", "market": "…", "source": "…", "attribution": "…" }. */
export async function operatorCarbon(get, url = config.carbonPriceUrl) {
  const d = await get(url).then(J), price = Number(d?.price), currency = String(d?.currency || '').toUpperCase(), t = Date.parse(d?.date);
  if (!(price > 0 && price < 1e5) || !/^[A-Z]{3}$/.test(currency) || !Number.isFinite(t)) throw new Error('operator carbon price file needs "date", "price" and a three-letter "currency"');
  const text = (v, n) => String(v || '').slice(0, n);
  return { date: new Date(t).toISOString().slice(0, 10), price, currency, market: text(d.market, 80) || 'Carbon price', instrument: text(d.market, 80) || 'Carbon price supplied by the operator of this site', source: text(d.source, 160) || 'supplied by the operator of this site', attribution: text(d.attribution, 300), operator: true, registry: ['operator-carbon'] };
}
/** Brent rows [{ t, v }] (oldest first) → the object the oil feed publishes. */
export function brentSeries(rows) {
  const tail = rows.slice(-260), last = tail[tail.length - 1], px = tail.map((r) => r.v), rets = px.slice(1).map((v, i) => Math.log(v / px[i]));
  if (!last || !(last.v > 5 && last.v < 500)) throw new Error('no plausible Brent observation');
  const m = rets.reduce((a, b) => a + b, 0) / rets.length, sd = Math.sqrt(rets.reduce((a, b) => a + (b - m) ** 2, 0) / (rets.length - 1));
  return { date: last.t, brent_usd_bbl: last.v, vol_annual: sd * Math.sqrt(252), history: tail.filter((_, i) => i % 5 === 0 || i === tail.length - 1) };
}
/** Brent from the open datasets/oil-prices package (EIA series repackaged, PDDL). */
export async function brentFromDataset(get) {
  const txt = await get('https://raw.githubusercontent.com/datasets/oil-prices/main/data/brent-daily.csv').then(T);
  return brentSeries(txt.trim().split('\n').slice(1).map((l) => l.split(',')).filter((r) => r.length >= 2 && Number.isFinite(Number(r[1]))).map((r) => ({ t: r[0], v: Number(r[1]) })));
}
/** Euro reference rates from the ECB Data Portal → rates per US dollar (own calculation), ECB values kept unchanged in eur_rates. */
export async function fxFromEcb(get) {
  const lines = (await get('https://data-api.ecb.europa.eu/service/data/EXR/D..EUR.SP00.A?lastNObservations=1&format=csvdata&detail=dataonly').then(T)).trim().split(/\r?\n/), hd = lines[0].split(','), iC = hd.indexOf('CURRENCY'), iT = hd.indexOf('TIME_PERIOD'), iV = hd.indexOf('OBS_VALUE');
  if (iC < 0 || iT < 0 || iV < 0) throw new Error('unexpected reply');
  const rows = lines.slice(1).map((l) => l.split(',')).map((r) => ({ c: r[iC], t: r[iT], v: Number(r[iV]) })).filter((r) => /^[A-Z]{3}$/.test(r.c) && r.v > 0), date = rows.map((r) => r.t).sort().pop();
  const eur = Object.fromEntries(rows.filter((r) => r.t === date).map((r) => [r.c, r.v])); if (!(eur.USD > 0)) throw new Error('no US dollar reference rate');
  const rates = { USD: 1, EUR: 1 / eur.USD }; for (const [c, v] of Object.entries(eur)) if (c !== 'USD') rates[c] = v / eur.USD;
  return { date, base: 'USD', rates, eur_rates: eur, source: 'ECB', derived: 'Rates per USD are calculated from the ECB euro reference rates.' };
}

// ---- aerodromes: bundled OurAirports database, optional OpenStreetMap enrichment ----------------
const RADIUS_KM = 40;
/** Bundled database → the shape the case page uses. Closed runways are counted, not offered. */
async function bundledAerodromes(p) {
  const [list, m] = await Promise.all([airports.nearest(p.lat, p.lon, RADIUS_KM, 25), airports.meta()]);
  const fields = list.map((a) => { const rws = a.runways.filter((r) => !r.closed).sort((x, y) => y.len_m - x.len_m); return { name: a.name, ident: a.ident, icao: a.icao || (/^[A-Z]{4}$/.test(a.ident) ? a.ident : ''), iata: a.iata, kind: a.type, ele_m: a.elev_m, dist_km: a.dist_km, lat: a.lat, lon: a.lon, country: a.country, runways: rws, closed_runways: a.runways.length - rws.length, longest_m: rws[0]?.len_m || 0, source: 'OurAirports' }; });
  return { fields, nRunways: fields.reduce((n, f) => n + f.runways.length, 0), source: 'OurAirports', dataset_date: m.dataset_date };
}
/** Can extra detail from OpenStreetMap be requested in this deployment? (Commercial: only through the operator's own Overpass server.) */
export const osmAvailable = () => !!config.overpassUrl || open('overpass-public');
const OVERPASS_PUBLIC = ['https://overpass-api.de/api/interpreter', 'https://overpass.private.coffee/api/interpreter', 'https://maps.mail.ru/osm/tools/overpass/api/interpreter'];
async function overpassAerodromes(p, get) {
  const endpoints = config.overpassUrl ? [config.overpassUrl] : (await allowed('overpass-public')) ? OVERPASS_PUBLIC : [];
  if (!endpoints.length) throw policy('The public OpenStreetMap query servers are not for commercial use; the operator can configure its own (overpassUrl)');
  const q = `[out:json][timeout:25];(nwr["aeroway"~"^(aerodrome|heliport)$"](around:${RADIUS_KM * 1000},${p.lat},${p.lon});way["aeroway"="runway"](around:${RADIUS_KM * 1000},${p.lat},${p.lon}););out tags geom 250;`;
  // public Overpass instances are volunteer-run and sometimes busy: try each in turn
  let d = null, lastErr = null;
  for (const ep of endpoints) {
    try { d = await get(`${ep}?data=${encodeURIComponent(q)}`).then(J); break; } catch (e) { lastErr = e; }
  }
  if (!d) throw new Error(lastErr?.message === 'Failed to fetch' ? 'Map service busy' : lastErr?.message || 'Map service busy');
  const hav = (a, b) => airports.distanceKm(a, b) * 1000;
  const centre = (e) => (e.lat != null ? { lat: e.lat, lon: e.lon } : e.bounds ? { lat: (e.bounds.minlat + e.bounds.maxlat) / 2, lon: (e.bounds.minlon + e.bounds.maxlon) / 2 } : null);
  const fields = [], runways = [];
  for (const e of d.elements || []) {
    const t = e.tags || {}, c = centre(e); if (!c) continue;
    if (t.aeroway === 'runway') {
      const g = e.geometry || []; let len = Number(String(t.length || '').replace(/[^0-9.]/g, '')) || 0, hdg = null;
      if (g.length >= 2) { const a = g[0], b = g[g.length - 1]; if (!len) len = hav(a, b); hdg = airports.bearingDeg(a, b) % 180; }
      runways.push({ ref: t.ref || '', len_m: Math.round(len), surface: t.surface || 'unknown', width_m: Number(t.width) || null, heading_deg: hdg == null ? null : Math.round(hdg), heading_source: 'map geometry', ...c });
    } else fields.push({ name: t.name || t['name:en'] || t.icao || 'Unnamed aerodrome', icao: t.icao || '', iata: t.iata || '', kind: t.aeroway, ele_m: t.ele != null && t.ele !== '' && Number.isFinite(Number(t.ele)) ? Number(t.ele) : null, dist_km: hav(p, c) / 1000, ...c, source: 'OpenStreetMap' });
  }
  for (const f of fields) { f.runways = runways.filter((r) => hav(f, r) < 4500).sort((a, b) => b.len_m - a.len_m); f.longest_m = f.runways[0]?.len_m || 0; }
  fields.sort((a, b) => a.dist_km - b.dist_km);
  return { fields, nRunways: runways.length, source: 'OpenStreetMap' };
}
/** Keep the bundled record as the reference; let the map fill gaps (missing runways, surface, width) and add fields it alone knows. */
function mergeAerodromes(base, osm) {
  const used = new Set(), same = (a, b) => (a.icao && a.icao === b.icao) || (a.iata && a.iata === b.iata) || airports.distanceKm(a, b) < 1.5;
  const fields = base.fields.map((f) => {
    const o = osm.fields.find((x) => !used.has(x) && same(f, x)); if (!o) return f; used.add(o);
    if (!f.runways.length && o.runways.length) return { ...f, runways: o.runways, longest_m: o.longest_m, source: 'OurAirports + OpenStreetMap runways' };
    let filled = false;
    const runways = f.runways.map((r) => { const m = o.runways.find((x) => (x.ref && x.ref.replace(/\s/g, '') === r.ref) || Math.abs(x.len_m - r.len_m) < 0.06 * r.len_m); if (!m) return r; const surface = r.surface === 'unknown' && m.surface !== 'unknown' ? m.surface : r.surface, width_m = r.width_m ?? m.width_m; if (surface !== r.surface || width_m !== r.width_m) filled = true; return { ...r, surface, width_m }; });
    return filled ? { ...f, runways, source: 'OurAirports + OpenStreetMap details' } : f;
  });
  for (const o of osm.fields) if (!used.has(o) && !fields.some((f) => same(f, o))) fields.push(o);
  fields.sort((a, b) => a.dist_km - b.dist_km);
  return { fields, nRunways: fields.reduce((n, f) => n + f.runways.length, 0), source: 'OurAirports + OpenStreetMap', dataset_date: base.dataset_date };
}
/** Airports matching a code or name, as place-search results (same shape as geocode()). Works offline. */
export async function searchAirports(text, limit = 6) {
  return (await airports.search(text, limit)).map((a) => ({ name: a.name, admin: [...new Set([a.icao || a.ident, a.iata, a.country].filter(Boolean))].join(' · '), lat: a.lat, lon: a.lon, elev_m: a.elev_m, country: a.country, airport: true, kind: a.type, longest_m: a.longest_m }));
}

// ---- cloud snapshot -----------------------------------------------------------------------------
const SNAPSHOT_PATH = 'data/snapshot.json';
let snap = null, snapAt = 0, snapBusy = null;
/**
 * The newest reachable cloud snapshot: this origin's copy and the copies on the hosts listed in mirrors.json
 * (a mirror that deploys from the repository carries an older file than the host running the schedule).
 * Falls back to the last snapshot saved on this device. Returns null when there has never been one.
 */
export async function loadSnapshot({ force = false } = {}) {
  if (snap !== null && Date.now() - snapAt < (force ? MIN : 10 * MIN)) return snap || null;
  if (snapBusy) return snapBusy;
  snapBusy = (async () => {
    const base = pageBase(), urls = base ? [new URL(SNAPSHOT_PATH, base).href] : [];
    let mirrors = globalThis.__AEROSUITE_MIRRORS__ || ls.get('mirrors', []);
    if (base) try { mirrors = ((await get(new URL('mirrors.json', base).href).then(J)).mirrors || []).map((m) => m.url).filter(Boolean); ls.set('mirrors', mirrors); } catch { /* offline: use the saved list */ }
    for (const m of mirrors) try { const u = new URL(SNAPSHOT_PATH, String(m).replace(/\/?$/, '/')).href; if (!urls.includes(u)) urls.push(u); } catch { /* malformed address */ }
    const when = (x) => Date.parse(x?.generated) || 0;
    const got = (await Promise.all(urls.map((u) => get(u).then(J).then((d) => (d && d.feeds && when(d) ? { ...d, from: u } : null)).catch(() => null)))).filter(Boolean).sort((a, b) => when(b) - when(a));
    let saved = await idb.get('snapshot'); const carried = globalThis.__AEROSUITE_ESSENTIAL__?.snapshot; // the single-file copy carries the snapshot of its build
    if (carried?.feeds && when(carried) > when(saved)) saved = { ...carried, from: null, carried: true };
    const best = got[0] && when(got[0]) >= when(saved) ? got[0] : saved || null;
    if (best && best === got[0] && when(best) > when(saved)) idb.set('snapshot', best);
    snap = best || false; snapAt = Date.now();
    emit('live', { id: 'snapshot' });
    return best;
  })().finally(() => { snapBusy = null; });
  return snapBusy;
}
/** Loader for feeds that exist only in the snapshot (their providers do not allow browser requests). */
/** A snapshot feed may be used if every registry source it was built from is allowed here (a commercial deployment refuses feeds that do not say). */
export async function feedAllowed(f) {
  if (!config.commercial) return true;
  if (!Array.isArray(f?.registry) || !f.registry.length) return false;
  for (const id of f.registry) if (!(await allowed(id))) return false;
  return true;
}
function fromSnapshot(id) {
  return async () => {
    const s = await loadSnapshot({ force: true }), f = s?.feeds?.[id];
    if (!f?.data) throw new Error(!s ? 'Cloud snapshot not reachable yet' : f?.error ? `Provider unavailable: ${f.error}` : 'Not in the cloud snapshot');
    if (!(await feedAllowed(f))) throw policy('The cloud snapshot was built from a source this deployment may not use');
    return { ...f.data, _ts: f.ts, _source: f.source, _url: f.url, _error: f.ok ? null : f.error || 'last fetch failed' };
  };
}
const SNAPSHOT_PARAMS = { macro: { country: 'WLD' } };
/** Put snapshot values into the on-device cache wherever they are newer than what is there. Returns how many were used. */
export async function seedFromSnapshot(opt) {
  await loadConfig();
  const s = await loadSnapshot(opt); if (!s) return 0;
  let n = 0;
  for (const [id, f] of Object.entries(s.feeds)) {
    const c = CONNECTORS[id]; if (!c || !f?.data || !Number.isFinite(f.ts) || !(await feedAllowed(f))) continue;
    if (id === 'carbon') continue; // the carbon panel is filtered market by market in its own loader
    const params = SNAPSHOT_PARAMS[id] || {}, key = c.key(params), ck = `live.${id}.${key}`, cur = await idb.get(ck);
    if (cur && cur.ts >= f.ts) continue;
    const data = c.cloud || id === 'oil' ? { ...f.data, _ts: f.ts, _source: f.source, _url: f.url, _error: f.ok ? null : f.error || 'last fetch failed' } : f.data;
    await idb.set(ck, { data, ts: f.ts, via: 'snapshot' }); setStatus(id, { ts: f.ts, ok: true, error: null, key, params, via: 'snapshot' }); n++;
  }
  return n;
}

/** Connector catalogue. `key(p)` identifies a cache entry; `load(p)` fetches and normalises. */
export const CONNECTORS = {
  weather: {
    title: 'Surface weather and winds aloft',
    get provider() { const pl = weatherPlan(); return pl.length ? pl.map((id) => WEATHER[id].label + (id === 'met-norway' || id === 'nws' ? ' (completed from the NOAA GFS grid)' : '')).join(', then ') : 'Switched off (weatherProvider: none)'; },
    get home() { return { 'open-meteo': 'https://open-meteo.com', 'met-norway': 'https://api.met.no', nws: 'https://www.weather.gov/documentation/services-web-api', 'noaa-grid': 'https://nomads.ncep.noaa.gov/' }[weatherPlan()[0]] || 'https://nomads.ncep.noaa.gov/'; },
    get ttl() { return weatherPlan()[0] === 'met-norway' ? 30 * MIN : 15 * MIN; }, // MET Norway asks clients not to re-request before the reply expires (about half an hour)
    needs: 'site', reg: () => { const ids = weatherPlan().map((id) => WEATHER[id].reg()).filter(Boolean); return ids.length && !ids.includes('noaa-gfs') ? [...ids, 'noaa-gfs'] : ids; },
    use: 'Site temperature, pressure, humidity, wind and gusts for field performance, icing, thermal and ECS loads; winds and temperatures aloft for mission fuel.',
    key: (p) => `${modeKey()}${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) {
      const errs = [];
      for (const id of weatherPlan()) {
        const P = WEATHER[id], reg = P.reg();
        if (!reg || !(await allowed(reg))) { errs.push(`${P.label}: not licensed for this deployment`); continue; }
        try { const w = await P.load(p, get); return finishWeather(id === 'open-meteo' || id === 'noaa-grid' ? w : await gridFill(w, p), id); } catch (e) { errs.push(`${P.label}: ${e.message || 'request failed'}`); }
      }
      throw errs.length ? new Error(errs.length === 1 ? errs[0].replace(/^[^:]+: /, '') : errs.join('; ')) : policy('The weather feed is switched off in this deployment');
    },
  },
  climate: {
    title: 'Site design temperatures',
    get provider() { return omId() ? 'Open-Meteo historical reanalysis (ERA5, last 12 months); NCEP-DOE Reanalysis 2 statistics bundled with the app as fall-back' : 'NCEP-DOE Reanalysis 2 statistics (NOAA PSL), bundled with the app'; },
    get home() { return omId() ? 'https://open-meteo.com/en/docs/historical-weather-api' : 'https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html'; }, ttl: 14 * DAY, needs: 'site', reg: () => [omId(), 'noaa-reanalysis2'].filter(Boolean),
    use: 'Hot-day and cold-day temperatures for hot-and-high take-off limits, cooling margins and cold-soak. Works without a connection.',
    key: (p) => `${modeKey()}${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) {
      // refinement where Open-Meteo may be used: the site's own last twelve months from the ERA5 archive
      if (omId()) try {
        const end = new Date(Date.now() - 6 * DAY), start = new Date(end.getTime() - 365 * DAY), f = (x) => x.toISOString().slice(0, 10);
        const q = (a, f2) => { const s = a.filter((v) => v != null).sort((x, y) => x - y); return s.length ? s[Math.min(s.length - 1, Math.floor(f2 * s.length))] : null; };
        const d = await get(await omUrl('archive-', '/v1/archive', `latitude=${p.lat}&longitude=${p.lon}&start_date=${f(start)}&end_date=${f(end)}&daily=temperature_2m_max,temperature_2m_min,wind_speed_10m_max&wind_speed_unit=ms&timezone=UTC`)).then(J);
        const out = { from: f(start), to: f(end), hot99_C: q(d.daily.temperature_2m_max, 0.99), hot_mean_C: q(d.daily.temperature_2m_max, 0.5), cold01_C: q(d.daily.temperature_2m_min, 0.01), wind99_ms: q(d.daily.wind_speed_10m_max, 0.99), provider: 'open-meteo', basis: 'ERA5 reanalysis at the site, last 12 months' };
        if (out.hot99_C != null) return out;
      } catch { /* fall back to the bundled statistics */ }
      // foundation in every mode: multi-year reanalysis statistics bundled with the app (no network)
      if (!(await allowed('noaa-reanalysis2'))) throw policy('No climate source is licensed for this deployment');
      const c = await climateAt(p.lat, p.lon, Number.isFinite(p.elev_m) ? p.elev_m : await bundledElevation(p));
      return { ...c, provider: 'noaa-reanalysis2', basis: `NCEP-DOE Reanalysis 2, ${c.from}–${c.to} (wind ${c.wind_period}), regional values on a 2.5° grid (the strong-wind figure is conservative: the reanalysis is windier than finer ones over mid-latitude land)${c.elevation_correction_K ? `, moved ${c.elevation_correction_K > 0 ? '+' : ''}${c.elevation_correction_K} K to the site elevation` : ''}` };
    },
  },
  airquality: {
    title: 'Air quality and dust', get provider() { return omId() ? 'Open-Meteo / Copernicus CAMS; NOAA GEFS-Aerosols grid as fall-back' : 'NOAA GEFS-Aerosols forecast grid (via the cloud snapshot)'; },
    get home() { return omId() ? 'https://open-meteo.com/en/docs/air-quality-api' : 'https://registry.opendata.aws/noaa-gefs/'; }, ttl: HOUR, needs: 'site', reg: () => [omId(), 'noaa-gefs-aerosols', ...(config.accept.includes('nasa-geos-cf') ? ['nasa-geos-cf'] : [])].filter(Boolean),
    use: 'Dust and particulates for engine erosion, filter loading and visibility; a local air-quality baseline for emissions assessment.',
    key: (p) => `${modeKey()}${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) {
      if (omId()) try { const d = await get(await omUrl('air-quality-', '/v1/air-quality', `latitude=${p.lat}&longitude=${p.lon}&current=pm10,pm2_5,dust,aerosol_optical_depth,nitrogen_dioxide,ozone&timezone=UTC`)).then(J); return { time: d.current.time, pm10: d.current.pm10, pm2_5: d.current.pm2_5, dust: d.current.dust, aod: d.current.aerosol_optical_depth, no2: d.current.nitrogen_dioxide, o3: d.current.ozone, provider: 'open-meteo' }; } catch { /* fall back to the grid */ }
      if (!(await allowed('noaa-gefs-aerosols'))) throw policy('No air-quality source is licensed for this deployment');
      const g = await gridAir(p.lat, p.lon), gasOk = g.gases === 'nasa-geos-cf' && config.accept.includes('nasa-geos-cf'); if (!gasOk) { g.no2 = null; g.o3 = null; }
      return { ...g, provider: 'noaa-gefs-aerosols', note: `Particulates, dust and aerosol optical depth: NOAA GEFS-Aerosols forecast (2.5° grid). ${gasOk ? 'Nitrogen dioxide and ozone: NASA GEOS-CF forecast (2.5° grid), accepted by the operator of this site.' : 'Nitrogen dioxide and ozone are not shown: no keyless source with an explicit licence for commercial reuse publishes them worldwide. They appear with an Open-Meteo key, or once the operator accepts the NASA GEOS-CF source ("accept": ["nasa-geos-cf"]).'}` };
    },
  },
  marine: {
    title: 'Sea state', get provider() { return omId() ? 'Open-Meteo marine models; DWD wave model grid as fall-back' : 'DWD global wave model GWAM grid (via the cloud snapshot)'; },
    get home() { return omId() ? 'https://open-meteo.com/en/docs/marine-weather-api' : 'https://opendata.dwd.de/weather/maritime/wave_models/gwam/'; }, ttl: HOUR, needs: 'site', reg: () => [omId(), 'dwd-gwam'].filter(Boolean),
    use: 'Wave height for offshore helicopter operations, ditching and flotation assessments.',
    key: (p) => `${modeKey()}${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) {
      if (omId()) try { const d = await get(await omUrl('marine-', '/v1/marine', `latitude=${p.lat}&longitude=${p.lon}&current=wave_height,wave_period,wave_direction&timezone=UTC`)).then(J); return { time: d.current.time, wave_height_m: d.current.wave_height, wave_period_s: d.current.wave_period, wave_dir_deg: d.current.wave_direction, provider: 'open-meteo' }; } catch { /* fall back to the grid */ }
      if (!(await allowed('dwd-gwam'))) throw policy('No sea-state source is licensed for this deployment');
      return { ...(await gridMarine(p.lat, p.lon)), provider: 'dwd-gwam' };
    },
  },
  aerodromes: {
    title: 'Nearby aerodromes and runways', get provider() { return `OurAirports (public domain), bundled with the app${osmAvailable() ? '; optional OpenStreetMap enrichment' : ''}`; }, home: 'https://ourairports.com/data/', ttl: 30 * DAY, needs: 'site', bundled: true,
    reg: () => ['ourairports', ...(osmAvailable() ? ['osm-data', ...(config.overpassUrl ? [] : ['overpass-public'])] : [])],
    use: 'Runway length, width, surface, true heading, threshold elevation and displaced thresholds for take-off and landing analyses and alternates. Works without a connection.',
    key: (p) => `oa:${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) {
      // primary: the bundled database (always there, also offline); OpenStreetMap only on request or if the bundle is missing
      let base = null, baseErr = null, osm = null, osmErr = null;
      try { base = await bundledAerodromes(p); } catch (e) { baseErr = e; }
      if ((p.osm && osmAvailable()) || !base) { try { osm = await overpassAerodromes(p, get); } catch (e) { osmErr = e.message || 'Map service busy'; if (!base) throw new Error(baseErr?.message || osmErr); } }
      const out = base && osm ? mergeAerodromes(base, osm) : base || osm;
      return { ...out, fields: out.fields.slice(0, 25), osm: !!osm, osm_error: osm ? null : p.osm ? osmErr : null };
    },
  },
  spaceweather: {
    title: 'Space weather (planetary Kp index)', provider: 'NOAA Space Weather Prediction Center', home: 'https://www.swpc.noaa.gov', ttl: 30 * MIN, reg: () => ['noaa-swpc'],
    use: 'Geomagnetic activity affecting GNSS accuracy, HF communication and high-latitude radiation exposure.',
    key: () => 'kp',
    async load(p, get) { const d = await get('https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json').then(J); const rows = Array.isArray(d[0]) ? d.slice(1).map((r) => ({ time_tag: r[0], Kp: Number(r[1]) })) : d; const last = rows[rows.length - 1]; return { time: last.time_tag, kp: Number(last.Kp), history: rows.slice(-24).map((r) => ({ t: r.time_tag, kp: Number(r.Kp) })) }; },
  },
  fx: {
    title: 'Currency exchange rates', provider: 'European Central Bank euro reference rates (ECB Data Portal; Frankfurter as fall-back)', home: 'https://data.ecb.europa.eu', ttl: 12 * HOUR, reg: () => ['ecb-statistics', 'frankfurter'],
    use: 'Convert costs and revenues between currencies in the economics suite.',
    key: () => 'usd',
    async load(p, get) {
      try { return await fxFromEcb(get); }
      catch { const d = await get('https://api.frankfurter.dev/v1/latest?base=USD').then(J); return { date: d.date, base: 'USD', rates: { USD: 1, ...d.rates }, source: 'ECB via Frankfurter' }; }
    },
  },
  oil: {
    title: 'Crude oil benchmark (Brent) and derived jet-fuel price', get provider() { return open('oil-prices-dataset') ? 'US EIA daily spot series via the open "oil-prices" dataset' : 'US Energy Information Administration daily spot series (via the cloud snapshot)'; },
    get home() { return open('oil-prices-dataset') ? 'https://github.com/datasets/oil-prices' : 'https://www.eia.gov/dnav/pet/hist/RBRTED.htm'; }, ttl: 12 * HOUR, reg: () => (open('oil-prices-dataset') ? ['oil-prices-dataset', 'eia'] : ['eia']),
    use: 'Fuel price for operating cost, mission cost index and fuel-price risk. Jet fuel is derived from Brent with an adjustable refining (crack) margin, so treat it as an indicator, not a quotation.',
    key: () => 'brent',
    // non-commercial: straight from the open dataset; commercial: the EIA series the snapshot job reads from eia.gov
    async load(p, get) { return (await allowed('oil-prices-dataset')) ? brentFromDataset(get) : fromSnapshot('oil')(); },
  },
  jetfuel: {
    title: 'Jet-fuel spot price', provider: 'US Energy Information Administration daily spot series (via the cloud snapshot)', home: 'https://www.eia.gov/dnav/pet/pet_pri_spt_s1_d.htm', ttl: 6 * HOUR, cloud: true, reg: () => ['eia'],
    use: 'Fuel price for operating cost and fuel-price risk: the U.S. Gulf Coast kerosene-type jet fuel spot price, converted to USD/kg. A regional wholesale benchmark, not an into-plane price at your airport.',
    key: () => 'usgc', load: fromSnapshot('jetfuel'),
  },
  carbon: {
    title: 'Carbon prices (official, by market)',
    get provider() { return `${config.carbonPriceUrl ? 'Operator of this site; ' : ''}European Commission (EU ETS), UK ETS Authority, Gouvernement du Québec (California–Québec market) — via the cloud snapshot`; },
    home: 'https://taxation-customs.ec.europa.eu/carbon-border-adjustment-mechanism/price-cbam-certificates_en', ttl: 6 * HOUR, cloud: true,
    reg: () => [...(config.carbonPriceUrl ? ['operator-carbon'] : []), 'ec-cbam-price', 'uk-ets-price', 'quebec-wci-auction', ...(config.accept.includes('eex-auction') ? ['eex-auction'] : [])],
    use: 'Carbon cost in the economics suite. Official allowance prices published under open licences, one per market: the case takes the market of the site\'s country (the EU ETS price where the country has no listed market) unless you choose another; converted to USD with the exchange-rate feed. Each is a published reference price with its own date — not a live quotation.',
    key: () => carbonMode(),
    async load(p, get) {
      const markets = {}, errors = {}; let ts = 0, src = null, url = null, err = null, direct = false;
      const s = await loadSnapshot({ force: true }).catch(() => null), f = s?.feeds?.carbon;
      if (f?.data?.markets) for (const [id, m] of Object.entries(f.data.markets)) { if (!Array.isArray(m?.registry) || !m.registry.length) continue; let ok = true; for (const r of m.registry) if (!(await allowed(r)) || (r === 'eex-auction' && !config.accept.includes('eex-auction'))) ok = false; if (ok) markets[id] = m; }
      if (Object.keys(markets).length) { ts = f.ts; src = f.source; url = f.url; err = f.ok ? null : f.error || 'last fetch failed'; }
      else if (await allowed('uk-ets-price')) try { const d = await ukEtsCarbon(get); markets['uk-ets'] = d; ts = Date.now(); src = d.source; url = d.url; direct = true; } catch (e) { errors['uk-ets'] = e.message; } // no snapshot: the one market a browser can read directly
      if (config.carbonPriceUrl && (await allowed('operator-carbon'))) try { markets.operator = { ...(await operatorCarbon(get)), id: 'operator', scheme: 'Operator', countries: [] }; ts = Date.now(); direct = true; } catch (e) { errors.operator = e.message; }
      if (!Object.keys(markets).length) throw new Error(!s ? 'Cloud snapshot not reachable yet' : Object.values(errors)[0] || 'No carbon price in the cloud snapshot');
      return { markets, order: Object.keys(markets), errors, _ts: ts, _source: src, _url: url, _error: err, ...(direct ? { _direct: true } : {}) };
    },
  },
  rates: {
    title: 'Policy interest rates (US and euro area)', provider: 'Federal Reserve Board H.15 and European Central Bank Data Portal (via the cloud snapshot)', home: 'https://data.ecb.europa.eu', ttl: 6 * HOUR, cloud: true, reg: () => ['frb-h15', 'nyfed-markets', 'ecb-statistics'],
    use: 'Reference points for discount-rate and financing assumptions. Shown for information; not written into the case.',
    key: () => 'policy', load: fromSnapshot('rates'),
  },
  macro: {
    title: 'Inflation and lending rates', provider: 'World Bank Open Data', home: 'https://data.worldbank.org', ttl: 7 * DAY, reg: () => ['world-bank'],
    use: 'Inflation and interest-rate assumptions for discounted cash flow, financing and life-cycle cost.',
    key: (p) => p.country || 'WLD',
    async load(p, get) {
      const cc = p.country || 'WLD', one = async (ind) => { try { const d = await get(`https://api.worldbank.org/v2/country/${cc}/indicator/${ind}?format=json&mrnev=1`).then(J); const r = d[1]?.[0]; return r ? { value: r.value, year: r.date, country: r.country?.value } : null; } catch { return null; } };
      const [infl, lend, gdp] = await Promise.all([one('FP.CPI.TOTL.ZG'), one('FR.INR.LEND'), one('NY.GDP.MKTP.KD.ZG')]);
      if (!infl && !lend && !gdp) throw new Error('No data for this country');
      return { country: infl?.country || lend?.country || cc, inflation_pct: infl?.value ?? null, inflation_year: infl?.year, lending_pct: lend?.value ?? null, lending_year: lend?.year, gdp_growth_pct: gdp?.value ?? null, gdp_year: gdp?.year };
    },
  },
  grid: {
    title: 'Electricity grid carbon intensity', get provider() { return `Ember yearly electricity data (national averages, bundled with the app)${gbLive() ? '; NESO Carbon Intensity API for the current value in Great Britain' : ''}`; },
    home: 'https://ember-energy.org/data/yearly-electricity-data/', ttl: 30 * MIN, bundledNote: true, reg: () => ['grid-factors', ...(gbLive() ? ['gb-carbon-intensity'] : [])],
    use: 'Well-to-wake CO₂ of battery charging for electric aircraft: the national average of the site\'s country for the latest year published (world average where the country is not listed). Enter your own factor for a specific supplier or hour.',
    key: (p) => `cc:${p?.country || 'WLD'}`,
    async load(p, get) {
      const t = await gridFactors(), cc = String(p?.country || '').toUpperCase(), rec = t.countries?.[cc] || null, base = rec || t.world;
      const out = { country: rec ? cc : 'World', gCO2_kWh: base.g, year: base.year, basis: `${rec ? 'national' : 'world'} average ${base.year}`, national_g: base.g, source: t.source, method: t.method };
      if (cc === 'GB' && gbLive() && (await allowed('gb-carbon-intensity'))) try { const d = await get('https://api.carbonintensity.org.uk/intensity').then(J), x = d.data[0], v = x.intensity.actual ?? x.intensity.forecast; if (Number.isFinite(v)) Object.assign(out, { gCO2_kWh: v, index: x.intensity.index, from: x.from, basis: 'current half-hour, Great Britain (NESO)' }); } catch { /* keep the national average */ }
      return out;
    },
  },
  literature: {
    title: 'Latest research literature', provider: 'OpenAlex scholarly index', home: 'https://openalex.org', ttl: DAY, reg: () => ['openalex'],
    use: 'Recent and most-cited publications for each suite, to check methods and find validation data.',
    key: (p) => `${p.q}|${p.sort}`,
    async load(p, get) {
      const sort = p.sort === 'cited' ? 'cited_by_count:desc' : 'publication_date:desc', from = p.sort === 'cited' ? '' : `,from_publication_date:${new Date(Date.now() - 540 * DAY).toISOString().slice(0, 10)}`;
      const d = await get(`https://api.openalex.org/works?search=${encodeURIComponent(p.q)}&filter=type:article${from}&sort=${sort}&per-page=12&select=id,title,publication_date,doi,cited_by_count,primary_location,authorships,open_access${config.openAlexApiKey ? `&api_key=${encodeURIComponent(config.openAlexApiKey)}` : ''}`).then(J);
      return { total: d.meta?.count, items: (d.results || []).filter((w) => w.title).map((w) => ({ title: w.title, date: w.publication_date, cited: w.cited_by_count, url: w.doi || w.id, venue: w.primary_location?.source?.display_name || '', authors: (w.authorships || []).slice(0, 3).map((a) => a.author?.display_name).filter(Boolean).join(', ') + ((w.authorships || []).length > 3 ? ' et al.' : ''), oa: !!w.open_access?.is_oa })) };
    },
  },
  opensource: {
    title: 'Open-source solvers and tools', get provider() { return `Curated list bundled with the app${open('github-api') ? '; GitHub repository search for further projects' : ''}`; }, home: 'https://github.com', ttl: DAY, reg: () => (open('github-api') ? ['tools-list', 'github-api'] : ['tools-list']),
    use: 'Maintained open-source codes for the higher-fidelity models each suite hands off, and reference implementations to cross-check against. Works without a connection.',
    key: (p) => modeKey() + (p.suite || '') + '|' + p.q,
    async load(p, get) {
      // foundation in every mode: the curated list shipped with the app (js/data/tools.json)
      const t = await toolsList().catch(() => null), curated = (t?.suites?.[p.suite] || []).map((x) => ({ name: x.name, desc: x.what, stars: null, pushed: '', url: x.url, lang: '', license: x.licence, curated: true }));
      // extra, non-commercial deployments only: a live repository search (GitHub's terms carry no explicit grant for embedding it in a commercial product)
      let extra = [], total = null;
      if (p.q && (await allowed('github-api')) && open('github-api')) try { const d = await get(`https://api.github.com/search/repositories?q=${encodeURIComponent(p.q)}&sort=stars&order=desc&per_page=10`, { headers: { Accept: 'application/vnd.github+json' } }).then(J); total = d.total_count; extra = (d.items || []).map((r) => ({ name: r.full_name, desc: r.description || '', stars: r.stargazers_count, pushed: r.pushed_at, url: r.html_url, lang: r.language || '', license: r.license?.spdx_id || '' })); } catch (e) { if (!curated.length) throw e; }
      if (!curated.length && !extra.length) throw new Error('No tool list for this suite');
      const seen = new Set(curated.map((x) => x.url.replace(/\/$/, '').toLowerCase()));
      return { total, curated: curated.length, items: [...curated, ...extra.filter((x) => !seen.has(x.url.toLowerCase()))] };
    },
  },
};

/** Bundled JSON next to the app (licence registry style): read once, from the module's own location, the page or the single-file copy. */
const bundledJson = (file, globalName) => { let pr = null; return () => (pr ||= (async () => {
  if (globalThis[globalName]) return globalThis[globalName];
  let url = null; try { if (import.meta.url) url = new URL(`../data/${file}`, import.meta.url); } catch { /* bundled copy */ }
  if (!url && pageBase()) url = new URL(`js/data/${file}`, pageBase());
  if (!url) throw new Error(`${file} is not available in this copy of the app`);
  if (isNode && url.protocol === 'file:') { const fs = 'node:fs/promises', { readFile } = await import(fs); return JSON.parse(await readFile(url, 'utf8')); }
  return fetch(url).then(J);
})().catch((e) => { pr = null; throw e; })); };
/** National grid emission factors (js/data/grid-factors.json). */
export const gridFactors = bundledJson('grid-factors.json', '__AEROSUITE_GRID_FACTORS__');
/** Curated open-source tools per suite (js/data/tools.json). */
export const toolsList = bundledJson('tools.json', '__AEROSUITE_TOOLS__');
/** The live GB grid value is a refinement for non-commercial deployments (its API terms carry an indemnity an operator must accept first). */
const gbLive = () => open('gb-carbon-intensity');

/** Cache key of the carbon feed: the official panel, with or without the operator's own price. */
export const carbonMode = () => (config.carbonPriceUrl ? 'panel+operator' : 'panel');
/**
 * The carbon market the case uses. `choice` is a market id or 'auto': the operator's price if there is one, else the
 * market whose scheme covers the whole of the site's country, else the EU ETS price (the largest market) as reference.
 * Accepts the panel ({ markets }) or a single record (returned as it is).
 */
export function carbonChoice(data, country = '', choice = 'auto') {
  if (!data) return null; if (!data.markets) return data;
  const ms = data.markets, ids = data.order?.filter((id) => ms[id]) || Object.keys(ms), cc = String(country || '').toUpperCase();
  if (choice && choice !== 'auto' && ms[choice]) return ms[choice];
  if (ms.operator) return ms.operator;
  const home = ids.find((id) => !ms[id].subnational && (ms[id].countries || []).includes(cc));
  return ms[home] || ms[ids.find((id) => ms[id].scheme === 'EU ETS')] || ms[ids[0]] || null;
}
/** Place-name search is always available: cities and airports are bundled with the app. */
export const placeSearchAvailable = () => true;
/**
 * Place search (interactive, not cached as a connector row). Cities come from the GeoNames database bundled with the
 * app, so the search works offline in every deployment; where Open-Meteo may be used its geocoder adds smaller places.
 */
export async function geocode(name) {
  await loadConfig();
  const local = await searchPlaces(name, 8).catch(() => []); let online = [];
  if (omId() && globalThis.navigator?.onLine !== false) try { const d = await get(await omUrl('geocoding-', '/v1/search', `name=${encodeURIComponent(name)}&count=8&language=en&format=json`)).then(J); online = (d.results || []).map((r) => ({ name: r.name, admin: [r.admin1, r.country].filter(Boolean).join(', '), lat: r.latitude, lon: r.longitude, elev_m: r.elevation, country: r.country_code, tz: r.timezone, source: 'Open-Meteo' })); } catch { /* the bundled list answers */ }
  const out = [...local]; for (const o of online) if (!out.some((x) => x.country === o.country && Math.abs(x.lat - o.lat) < 0.08 && Math.abs(x.lon - o.lon) < 0.08)) out.push(o);
  return out.slice(0, 12);
}
/** Terrain elevation for points [{ lat, lon }]: Open-Meteo where it may be used; otherwise the bundled data (airport within 10 km, city within 15 km, else the 0.5° ETOPO terrain grid). */
export async function elevations(points) {
  await loadConfig();
  if (omId()) try { const out = []; for (let i = 0; i < points.length; i += 100) { const ch = points.slice(i, i + 100); const d = await get(await omUrl('', '/v1/elevation', `latitude=${ch.map((p) => p.lat.toFixed(4)).join(',')}&longitude=${ch.map((p) => p.lon.toFixed(4)).join(',')}`)).then(J); out.push(...d.elevation); } return out; } catch { /* bundled data below */ }
  return Promise.all(points.map((p) => bundledElevation(p)));
}
/** Coarse terrain profile along points, from the bundled ETOPO grid only (no network): mean height of 0.5° cells. */
export const terrainProfile = (points) => terrain(points);
/** Attribution lines (text + link) for registry ids, from the licence registry. */
export async function attributions(ids) {
  const reg = await loadLicences().catch(() => null); if (!reg) return [];
  const seen = new Set(), out = [];
  for (const id of ids || []) { const s = reg.sources.find((x) => x.id === id); if (s?.attribution && !seen.has(s.attribution)) { seen.add(s.attribution); out.push({ id, text: s.attribution, url: s.attributionUrl || s.url || null, licence: s.licence, class: s.class }); } }
  return out;
}

// ---- cache + status ---------------------------------------------------------------------
export const status = ls.get('liveStatus', {}); // { id: { ts, ok, error, key } }
const inflight = new Map();
const setStatus = (id, s) => { status[id] = { ...(status[id] || {}), ...s }; ls.set('liveStatus', status); emit('live', { id }); };

/**
 * Fetch through the cache. Returns { data, ts, stale, error, fromCache }.
 * Fresh cache → returned immediately. Stale cache → returned immediately and refreshed in the background
 * (unless force, which waits). Offline or provider error → last cached value with `error` set.
 */
export async function live(id, params = {}, { force = false } = {}) {
  await loadConfig(); // the configuration decides the provider and the cache key
  const c = CONNECTORS[id], ck = `live.${id}.${c.key(params)}`, cached = await idb.get(ck), now = Date.now();
  const fresh = cached && now - cached.ts < c.ttl;
  if (fresh && !force) return { ...cached, stale: false, fromCache: true };
  const go = () => {
    if (inflight.has(ck)) return inflight.get(ck);
    const p = (async () => {
      try {
        const net = makeNet(), data = await c.load(params, net.get), offlineCopy = net.cachedAt != null, rec = { data, ts: net.cachedAt ?? (Number.isFinite(data?._ts) ? data._ts : Date.now()) };
        if (offlineCopy && cached && cached.ts >= rec.ts) return { ...cached, stale: true, fromCache: true, error: 'Offline' };
        if (data?._direct) rec.via = 'direct'; await idb.set(ck, rec); setStatus(id, { ts: rec.ts, ok: !offlineCopy, error: offlineCopy ? 'Offline' : null, key: c.key(params), params, via: data?._direct ? 'direct' : c.cloud || Number.isFinite(data?._ts) ? 'snapshot' : c.bundled && !data?.osm ? 'bundled' : 'direct' });
        return { ...rec, stale: offlineCopy, fromCache: offlineCopy, error: offlineCopy ? 'Offline' : undefined };
      }
      catch (e) { const msg = e.policy ? e.message : globalThis.navigator?.onLine === false ? 'Offline' : e.name === 'TimeoutError' ? 'Timed out' : e.message || 'Request failed'; setStatus(id, { ok: false, error: msg, policy: !!e.policy, tried: Date.now(), params }); return cached ? { ...cached, stale: true, fromCache: true, error: msg } : { data: null, ts: 0, stale: true, error: msg }; }
      finally { inflight.delete(ck); }
    })();
    inflight.set(ck, p); return p;
  };
  if (cached && !force) { go(); return { ...cached, stale: true, fromCache: true, refreshing: true }; }
  return go();
}

// ---- applying live data to the shared case ----------------------------------------------------
export const jetFromBrent = (usdBbl, crack) => usdKgFromUsdBbl(usdBbl * (1 + crack));
/** A quoted price older than this is not written into the case (the derived estimate or the manual value stays). */
const MAX_QUOTE_AGE = 35 * DAY;
const recent = (date) => Date.now() - Date.parse(date) < MAX_QUOTE_AGE;
/** Jet-fuel price for the case: the quoted spot series when it is recent, else the Brent-derived estimate. */
export function fuelPrice(jet, oil, crack) {
  if (Number.isFinite(jet?.usd_kg) && recent(jet.date)) return { usd_kg: jet.usd_kg, quoted: true, source: `US Gulf Coast jet fuel spot ${jet.usd_gal} USD/gal on ${jet.date} (EIA), ${JET_DENSITY_KG_L} kg/L` };
  if (Number.isFinite(oil?.brent_usd_bbl)) return { usd_kg: jetFromBrent(oil.brent_usd_bbl, crack), quoted: false, source: `Estimated from Brent ${oil.brent_usd_bbl} USD/bbl on ${oil.date} plus ${Math.round(crack * 100)}% refining margin` };
  return null;
}
/**
 * Carbon price for the case in USD/t, converted with the newest exchange rate available. Accepts the EU allowance
 * auction record ({ eur_t }), the UK ETS determination and an operator-supplied price ({ price, currency }).
 * A price fixed for a period (valid_to) is used until the period ends; any other price only while it is recent.
 */
export function carbonPrice(carbon, fx) {
  const cur = carbon?.currency || (Number.isFinite(carbon?.eur_t) ? 'EUR' : null), price = Number.isFinite(carbon?.price) ? carbon.price : carbon?.eur_t;
  if (!cur || !Number.isFinite(price)) return null;
  if (carbon.valid_to ? Date.now() > Date.parse(carbon.valid_to) + MAX_QUOTE_AGE : !recent(carbon.date)) return null;
  const live = Number.isFinite(fx?.rates?.[cur]), rate = cur === 'USD' ? 1 : live ? fx.rates[cur] : cur === 'EUR' ? carbon.eur_per_usd : carbon.fx_per_usd, fxDate = live ? fx.date : carbon.fx_date;
  if (!Number.isFinite(rate) || rate <= 0) return null;
  const what = carbon.scheme_year ? `UK ETS carbon price for scheme year ${carbon.scheme_year}: ${price} GBP/t (UK ETS Authority determination, Open Government Licence v3.0)` : carbon.operator || carbon.market ? `${carbon.market} ${carbon.period || carbon.date}: ${price} ${cur}/t (${carbon.source})` : `EU ETS allowance (EUA) auction ${carbon.date}: ${price} EUR/t (EEX)`;
  return { usd_t: usdFromEur(price, rate), source: cur === 'USD' ? what : `${what}, ${Math.round(rate * 1e4) / 1e4} ${cur} per USD on ${fxDate}` };
}

/** Pull everything relevant for the current case site and write it into the case (with provenance). */
export async function refreshSite({ force = false } = {}) {
  await loadConfig();
  const s = state.case.site; if (s.lat == null || s.lon == null) return null;
  const p = { lat: s.lat, lon: s.lon }, [w, cl, mar] = await Promise.all([live('weather', p, { force }), live('climate', p, { force }), live('marine', p, { force })]);
  live('airquality', p, { force }).catch(() => {}); // shown on the Live data page; not written into the case
  const patch = {};
  // a field the provider does not supply keeps the value already in the case (see finishWeather for the stand-ins)
  if (w.data) Object.assign(patch, { T_C: w.data.T_C, p_hPa: w.data.p_hPa, rh: w.data.rh ?? s.rh, wind_ms: w.data.wind_ms ?? s.wind_ms, wind_dir_deg: w.data.wind_dir_deg ?? s.wind_dir_deg, gust_ms: w.data.gust_ms ?? w.data.wind_ms ?? s.gust_ms, precip_mm_h: w.data.precip_mm_h ?? 0, cloud_pct: w.data.cloud_pct ?? s.cloud_pct,
    visibility_m: w.data.visibility_m ?? s.visibility_m, freezing_level_m: w.data.freezing_level_m ?? s.freezing_level_m, winds_aloft: w.data.aloft || [], source: `${w.data.provider_label || 'Open-Meteo'}, ${w.data.time}Z${w.stale ? ' (cached)' : ''}`, wx_provider: w.data.provider || 'open-meteo', wx_note: w.data.note || '', updated: w.ts });
  if (w.data && patch.precip_mm_h > 0.1) { patch.runway_mu_brake = patch.T_C <= 0 ? 0.1 : 0.25; patch.runway_state = patch.T_C <= 0 ? 'contaminated (freezing precipitation)' : 'wet'; } else if (w.data) { patch.runway_mu_brake = 0.4; patch.runway_state = 'dry'; }
  if (cl.data?.hot99_C != null) patch.design_hot_C = cl.data.hot99_C;
  if (mar.data?.wave_height_m != null) patch.wave_height_m = mar.data.wave_height_m;
  if (state.settings.autoApplyLive && Object.keys(patch).length) patchCaseMany('site', patch);
  return { weather: w, climate: cl, marine: mar };
}
/** Global (location-independent) feeds: FX, oil, jet fuel, carbon, macro-economics, interest rates, space weather. */
export async function refreshGlobal({ force = false } = {}) {
  await loadConfig();
  const [fx, oil, kp, macro, jet, carbon, rates] = await Promise.all([live('fx', {}, { force }), live('oil', {}, { force }), live('spaceweather', {}, { force }), live('macro', { country: state.case.site.country || 'WLD' }, { force }), live('jetfuel', {}, { force }), live('carbon', {}, { force }), live('rates', {}, { force })]);
  live('grid', { country: state.case.site.country || '' }, { force }).catch(() => {}); // national grid factor for the site's country (shown on the Live data page)
  if (state.settings.autoApplyLive) {
    if (kp.data) patchCaseMany('site', { kp_index: kp.data.kp }, true);
    const e = {};
    const fp = fuelPrice(jet.data, oil.data, state.settings.fuelCrack);
    if (fp && state.case.prop.fuel !== 'Liquid hydrogen' && state.case.prop.type !== 'electric') { const k = state.case.prop.fuel === 'Avgas 100LL' ? 2.2 : state.case.prop.fuel?.startsWith('SAF') ? 2.84 : 1; e.fuel_usd_kg = Math.round(fp.usd_kg * k * 1000) / 1000; e.fuel_source = fp.source + (k !== 1 ? ` × ${k} for ${state.case.prop.fuel} (indicative ratio)` : ''); const vol = fp.quoted && Number.isFinite(jet.data.vol_annual) ? jet.data.vol_annual : oil.data?.vol_annual; if (Number.isFinite(vol)) e.fuel_vol = Math.round(vol * 1000) / 1000; }
    const cp = carbonPrice(carbonChoice(carbon.data, state.case.site.country, state.settings.carbonMarket), fx.data);
    if (cp) { e.carbon_usd_t = Math.round(cp.usd_t * 100) / 100; e.carbon_source = cp.source; }
    if (macro.data?.inflation_pct != null && macro.data.inflation_pct > -5 && macro.data.inflation_pct < 60) { e.inflation = Math.round(macro.data.inflation_pct * 10) / 1000; e.inflation_source = `${macro.data.country} ${macro.data.inflation_year}`; }
    if (fx.data) { e.fx_rates = fx.data.rates; e.fx_date = fx.data.date; }
    if (Object.keys(e).length) patchCaseMany('econ', e);
  }
  return { fx, oil, kp, macro, jet, carbon, rates };
}

// ---- scheduling ---------------------------------------------------------------------------
let timer = null;
export function startLive() {
  const tick = (force = false) => { if (!state.settings.autoRefreshLive || navigator.onLine === false) return; refreshSite({ force }).catch(() => {}); refreshGlobal({ force }).catch(() => {}); };
  // first paint: values from the cloud snapshot (at most a few hours old) go into the cache before the per-feed
  // requests start; a slow or missing snapshot never delays them by more than a moment
  loadConfig().then(() => (state.settings.autoRefreshLive ? Promise.race([seedFromSnapshot().catch(() => 0), new Promise((r) => setTimeout(r, 2500))]) : null)).then(() => tick());
  clearInterval(timer); timer = setInterval(tick, 5 * MIN);                        // TTLs decide what is actually re-fetched
  window.addEventListener('online', () => tick());
  document.addEventListener('visibilitychange', () => { if (document.visibilityState === 'visible') tick(); });
  navigator.serviceWorker?.addEventListener('message', (ev) => { if (ev.data?.type === 'periodic-refresh') tick(); });
}

/** Per-suite search phrases for the literature and open-source connectors. */
export const SUITE_QUERIES = {
  cfd: ['aircraft computational fluid dynamics turbulence modelling', 'topic:cfd'], fea: ['aircraft structure finite element analysis', 'topic:finite-element-analysis'],
  aeroelastic: ['aircraft flutter aeroelasticity', 'aeroelasticity flutter'], flightdyn: ['aircraft flight dynamics stability handling qualities', 'topic:flight-dynamics'],
  performance: ['aircraft performance takeoff climb range', 'aircraft performance analysis'], rotorcraft: ['helicopter rotor aeromechanics', 'rotorcraft helicopter simulation'],
  propulsion: ['aircraft gas turbine engine performance cycle', 'gas turbine engine cycle simulation'], propeller: ['propeller blade element momentum performance', 'propeller blade element momentum'],
  fatigue: ['aircraft fatigue crack growth damage tolerance', 'fatigue crack growth'], vibration: ['rotor dynamics modal analysis aircraft vibration', 'topic:rotordynamics'],
  acoustics: ['aircraft aeroacoustics noise prediction', 'aeroacoustics noise'], thermal: ['aircraft thermal management heat transfer', 'heat transfer simulation'],
  icing: ['aircraft icing ice accretion simulation', 'aircraft icing'], gear: ['aircraft landing gear dynamics shimmy', 'landing gear simulation'],
  crash: ['aircraft crashworthiness bird strike impact', 'crashworthiness impact simulation'], control: ['flight control system fault tolerant autopilot', 'topic:flight-controller'],
  avionics: ['GNSS inertial navigation sensor fusion aircraft', 'topic:sensor-fusion navigation'], hydmech: ['aircraft hydraulic actuator system gearbox', 'hydraulic system simulation'],
  electrical: ['hybrid electric aircraft propulsion battery', 'electric aircraft battery model'], fuelecs: ['aircraft fuel system environmental control system', 'fuel sloshing simulation'],
  composites: ['composite laminate failure delamination aerospace', 'topic:composite-materials laminate'], safety: ['aircraft system safety assessment reliability fault tree', 'topic:fault-tree-analysis'],
  mdao: ['multidisciplinary design optimization aircraft', 'topic:mdao'], mission: ['aircraft trajectory optimization fuel mission', 'aircraft trajectory optimization'],
  vvuq: ['verification validation uncertainty quantification simulation', 'topic:uncertainty-quantification'], economics: ['aircraft direct operating cost life cycle cost airline economics', 'aircraft cost model'],
};
export const stdAtSite = (s) => isa(s.elev_m || 0).T - 273.15;
