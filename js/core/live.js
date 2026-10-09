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
  if (['auto', 'open-meteo', 'met-norway', 'nws', 'none'].includes(o.weatherProvider)) c.weatherProvider = o.weatherProvider;
  if (typeof o.carbonPriceUrl === 'string' && o.carbonPriceUrl.trim()) { const b = pageBase(); try { c.carbonPriceUrl = httpsUrl(new URL(o.carbonPriceUrl.trim(), b || undefined).href) || (b && new URL(o.carbonPriceUrl.trim(), b).origin === new URL(b).origin ? new URL(o.carbonPriceUrl.trim(), b).href : ''); } catch { /* not an address */ } }
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
/** Providers to try, in order, under the current configuration. */
export function weatherPlan(cfg = config) {
  const w = cfg.weatherProvider, om = !cfg.commercial || !!cfg.openMeteoApiKey || cfg.accept.includes('open-meteo-free');
  if (w === 'none') return [];
  const auto = cfg.commercial ? [...(om ? ['open-meteo'] : []), 'met-norway', 'nws'] : ['open-meteo', 'met-norway'];
  if (w === 'auto' || !WEATHER[w] || (w === 'open-meteo' && !om)) return auto;
  return [w, ...auto.filter((x) => x !== w)];
}
/** Fill the documented stand-ins for fields a provider cannot supply and say so. */
export function finishWeather(w, id) {
  const P = WEATHER[id], missing = [], notes = [], out = { ...w, provider: id, provider_label: P.label + (w.station ? ` (station ${w.station})` : '') };
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
  out.missing = missing; out.note = notes.length ? `${P.label} does not supply ${missing.join(', ')}. ${notes.map((n) => n[0].toUpperCase() + n.slice(1)).join('. ')}.` : '';
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
async function overpassAerodromes(p, get) {
  const q = `[out:json][timeout:25];(nwr["aeroway"~"^(aerodrome|heliport)$"](around:${RADIUS_KM * 1000},${p.lat},${p.lon});way["aeroway"="runway"](around:${RADIUS_KM * 1000},${p.lat},${p.lon}););out tags geom 250;`;
  // public Overpass instances are volunteer-run and sometimes busy: try each in turn
  let d = null, lastErr = null;
  for (const ep of ['https://overpass-api.de/api/interpreter', 'https://overpass.kumi.systems/api/interpreter', 'https://overpass.private.coffee/api/interpreter', 'https://maps.mail.ru/osm/tools/overpass/api/interpreter']) {
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
const pageBase = () => { try { return typeof document !== 'undefined' && /^https?:/.test(document.baseURI) ? document.baseURI : null; } catch { return null; } };
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
    const saved = await idb.get('snapshot'), best = got[0] && when(got[0]) >= when(saved) ? got[0] : saved || null;
    if (best && best === got[0] && when(best) > when(saved)) idb.set('snapshot', best);
    snap = best || false; snapAt = Date.now();
    emit('live', { id: 'snapshot' });
    return best;
  })().finally(() => { snapBusy = null; });
  return snapBusy;
}
/** Loader for feeds that exist only in the snapshot (their providers do not allow browser requests). */
function fromSnapshot(id) {
  return async () => {
    const s = await loadSnapshot({ force: true }), f = s?.feeds?.[id];
    if (!f?.data) throw new Error(!s ? 'Cloud snapshot not reachable yet' : f?.error ? `Provider unavailable: ${f.error}` : 'Not in the cloud snapshot');
    return { ...f.data, _ts: f.ts, _source: f.source, _url: f.url, _error: f.ok ? null : f.error || 'last fetch failed' };
  };
}
const SNAPSHOT_PARAMS = { macro: { country: 'WLD' } };
/** Put snapshot values into the on-device cache wherever they are newer than what is there. Returns how many were used. */
export async function seedFromSnapshot(opt) {
  const s = await loadSnapshot(opt); if (!s) return 0;
  let n = 0;
  for (const [id, f] of Object.entries(s.feeds)) {
    const c = CONNECTORS[id]; if (!c || !f?.data || !Number.isFinite(f.ts)) continue;
    const params = SNAPSHOT_PARAMS[id] || {}, key = c.key(params), ck = `live.${id}.${key}`, cur = await idb.get(ck);
    if (cur && cur.ts >= f.ts) continue;
    const data = c.cloud ? { ...f.data, _ts: f.ts, _source: f.source, _url: f.url, _error: f.ok ? null : f.error || 'last fetch failed' } : f.data;
    await idb.set(ck, { data, ts: f.ts, via: 'snapshot' }); setStatus(id, { ts: f.ts, ok: true, error: null, key, params, via: 'snapshot' }); n++;
  }
  return n;
}

/** Connector catalogue. `key(p)` identifies a cache entry; `load(p)` fetches and normalises. */
export const CONNECTORS = {
  weather: {
    title: 'Surface weather and winds aloft', provider: 'Open-Meteo (national weather-service models)', home: 'https://open-meteo.com', ttl: 15 * MIN, needs: 'site',
    use: 'Site temperature, pressure, humidity, wind and gusts for field performance, icing, thermal and ECS loads; winds and temperatures aloft for mission fuel.',
    key: (p) => `${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) {
      const hv = ['freezing_level_height', 'visibility', ...LEVELS.flatMap((l) => [`temperature_${l}hPa`, `wind_speed_${l}hPa`, `wind_direction_${l}hPa`, `geopotential_height_${l}hPa`])].join(',');
      const d = await get(`https://api.open-meteo.com/v1/forecast?latitude=${p.lat}&longitude=${p.lon}&current=temperature_2m,relative_humidity_2m,surface_pressure,pressure_msl,wind_speed_10m,wind_direction_10m,wind_gusts_10m,precipitation,cloud_cover,weather_code&hourly=${hv}&forecast_hours=1&wind_speed_unit=ms&timezone=UTC`).then(J);
      const c = d.current, hr = d.hourly || {}, at = (k) => (hr[k] ? hr[k][0] : null);
      return {
        time: c.time, elev_m: d.elevation, T_C: c.temperature_2m, rh: c.relative_humidity_2m / 100, p_hPa: c.surface_pressure, qnh_hPa: c.pressure_msl, wind_ms: c.wind_speed_10m, wind_dir_deg: c.wind_direction_10m, gust_ms: c.wind_gusts_10m,
        precip_mm_h: c.precipitation, cloud_pct: c.cloud_cover, code: c.weather_code, visibility_m: at('visibility'), freezing_level_m: at('freezing_level_height'),
        aloft: LEVELS.map((l) => ({ hPa: l, alt_m: at(`geopotential_height_${l}hPa`), T_C: at(`temperature_${l}hPa`), speed_ms: at(`wind_speed_${l}hPa`), dir_deg: at(`wind_direction_${l}hPa`) })).filter((a) => a.alt_m != null && a.T_C != null),
      };
    },
  },
  climate: {
    title: 'Site design temperatures (last 12 months)', provider: 'Open-Meteo historical reanalysis (ERA5)', home: 'https://open-meteo.com/en/docs/historical-weather-api', ttl: 14 * DAY, needs: 'site',
    use: 'Hot-day and cold-day temperatures for hot-and-high take-off limits, cooling margins and cold-soak.',
    key: (p) => `${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) {
      const end = new Date(Date.now() - 6 * DAY), start = new Date(end.getTime() - 365 * DAY), f = (x) => x.toISOString().slice(0, 10);
      const d = await get(`https://archive-api.open-meteo.com/v1/archive?latitude=${p.lat}&longitude=${p.lon}&start_date=${f(start)}&end_date=${f(end)}&daily=temperature_2m_max,temperature_2m_min,wind_speed_10m_max&wind_speed_unit=ms&timezone=UTC`).then(J);
      const q = (a, f2) => { const s = a.filter((v) => v != null).sort((x, y) => x - y); return s.length ? s[Math.min(s.length - 1, Math.floor(f2 * s.length))] : null; };
      return { from: f(start), to: f(end), hot99_C: q(d.daily.temperature_2m_max, 0.99), hot_mean_C: q(d.daily.temperature_2m_max, 0.5), cold01_C: q(d.daily.temperature_2m_min, 0.01), wind99_ms: q(d.daily.wind_speed_10m_max, 0.99) };
    },
  },
  airquality: {
    title: 'Air quality and dust', provider: 'Open-Meteo / Copernicus CAMS', home: 'https://open-meteo.com/en/docs/air-quality-api', ttl: HOUR, needs: 'site',
    use: 'Dust and particulates for engine erosion, filter loading and visibility; a local air-quality baseline for emissions assessment.',
    key: (p) => `${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) { const d = await get(`https://air-quality-api.open-meteo.com/v1/air-quality?latitude=${p.lat}&longitude=${p.lon}&current=pm10,pm2_5,dust,aerosol_optical_depth,nitrogen_dioxide,ozone&timezone=UTC`).then(J); return { time: d.current.time, pm10: d.current.pm10, pm2_5: d.current.pm2_5, dust: d.current.dust, aod: d.current.aerosol_optical_depth, no2: d.current.nitrogen_dioxide, o3: d.current.ozone }; },
  },
  marine: {
    title: 'Sea state', provider: 'Open-Meteo marine models', home: 'https://open-meteo.com/en/docs/marine-weather-api', ttl: HOUR, needs: 'site',
    use: 'Wave height for offshore helicopter operations, ditching and flotation assessments.',
    key: (p) => `${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) { const d = await get(`https://marine-api.open-meteo.com/v1/marine?latitude=${p.lat}&longitude=${p.lon}&current=wave_height,wave_period,wave_direction&timezone=UTC`).then(J); return { time: d.current.time, wave_height_m: d.current.wave_height, wave_period_s: d.current.wave_period, wave_dir_deg: d.current.wave_direction }; },
  },
  aerodromes: {
    title: 'Nearby aerodromes and runways', provider: 'OurAirports (public domain), bundled with the app; optional OpenStreetMap enrichment', home: 'https://ourairports.com/data/', ttl: 30 * DAY, needs: 'site', bundled: true,
    use: 'Runway length, width, surface, true heading, threshold elevation and displaced thresholds for take-off and landing analyses and alternates. Works without a connection.',
    key: (p) => `oa:${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) {
      // primary: the bundled database (always there, also offline); OpenStreetMap only on request or if the bundle is missing
      let base = null, baseErr = null, osm = null, osmErr = null;
      try { base = await bundledAerodromes(p); } catch (e) { baseErr = e; }
      if (p.osm || !base) { try { osm = await overpassAerodromes(p, get); } catch (e) { osmErr = e.message || 'Map service busy'; if (!base) throw new Error(baseErr?.message || osmErr); } }
      const out = base && osm ? mergeAerodromes(base, osm) : base || osm;
      return { ...out, fields: out.fields.slice(0, 25), osm: !!osm, osm_error: osm ? null : p.osm ? osmErr : null };
    },
  },
  spaceweather: {
    title: 'Space weather (planetary Kp index)', provider: 'NOAA Space Weather Prediction Center', home: 'https://www.swpc.noaa.gov', ttl: 30 * MIN,
    use: 'Geomagnetic activity affecting GNSS accuracy, HF communication and high-latitude radiation exposure.',
    key: () => 'kp',
    async load(p, get) { const d = await get('https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json').then(J); const rows = Array.isArray(d[0]) ? d.slice(1).map((r) => ({ time_tag: r[0], Kp: Number(r[1]) })) : d; const last = rows[rows.length - 1]; return { time: last.time_tag, kp: Number(last.Kp), history: rows.slice(-24).map((r) => ({ t: r.time_tag, kp: Number(r.Kp) })) }; },
  },
  fx: {
    title: 'Currency exchange rates', provider: 'Frankfurter (European Central Bank reference rates)', home: 'https://frankfurter.dev', ttl: 12 * HOUR,
    use: 'Convert costs and revenues between currencies in the economics suite.',
    key: () => 'usd',
    async load(p, get) { const d = await get('https://api.frankfurter.dev/v1/latest?base=USD').then(J); return { date: d.date, base: 'USD', rates: { USD: 1, ...d.rates } }; },
  },
  oil: {
    title: 'Crude oil benchmark (Brent) and derived jet-fuel price', provider: 'US EIA daily spot series via the open "oil-prices" dataset', home: 'https://github.com/datasets/oil-prices', ttl: 12 * HOUR,
    use: 'Fuel price for operating cost, mission cost index and fuel-price risk. Jet fuel is derived from Brent with an adjustable refining (crack) margin, so treat it as an indicator, not a quotation.',
    key: () => 'brent',
    async load(p, get) {
      const txt = await get('https://raw.githubusercontent.com/datasets/oil-prices/main/data/brent-daily.csv').then(T), rows = txt.trim().split('\n').slice(1).map((l) => l.split(',')).filter((r) => r.length >= 2 && Number.isFinite(Number(r[1])));
      const tail = rows.slice(-260), last = tail[tail.length - 1], px = tail.map((r) => Number(r[1])), rets = px.slice(1).map((v, i) => Math.log(v / px[i]));
      const m = rets.reduce((a, b) => a + b, 0) / rets.length, sd = Math.sqrt(rets.reduce((a, b) => a + (b - m) ** 2, 0) / (rets.length - 1));
      return { date: last[0], brent_usd_bbl: Number(last[1]), vol_annual: sd * Math.sqrt(252), history: tail.filter((_, i) => i % 5 === 0 || i === tail.length - 1).map((r) => ({ t: r[0], v: Number(r[1]) })) };
    },
  },
  jetfuel: {
    title: 'Jet-fuel spot price', provider: 'US Energy Information Administration daily spot series (via the cloud snapshot)', home: 'https://www.eia.gov/dnav/pet/pet_pri_spt_s1_d.htm', ttl: 6 * HOUR, cloud: true,
    use: 'Fuel price for operating cost and fuel-price risk: the U.S. Gulf Coast kerosene-type jet fuel spot price, converted to USD/kg. A regional wholesale benchmark, not an into-plane price at your airport.',
    key: () => 'usgc', load: fromSnapshot('jetfuel'),
  },
  carbon: {
    title: 'Carbon price (EU ETS allowance)', provider: 'EEX primary auctions of EU allowances, the EU common auction platform (via the cloud snapshot)', home: 'https://www.eex.com/en/markets/environmentals/eu-ets1-eu-ets2-auctions/eu-ets1-auctions', ttl: 6 * HOUR, cloud: true,
    use: 'Carbon cost in the economics suite: the clearing price of the latest EU allowance (EUA) auction in EUR per tonne of CO₂, converted to USD with the exchange-rate feed. Applies to flights covered by the EU ETS; other schemes price carbon differently.',
    key: () => 'eua', load: fromSnapshot('carbon'),
  },
  rates: {
    title: 'Policy interest rates (US and euro area)', provider: 'Federal Reserve (via FRED) and European Central Bank Data Portal (via the cloud snapshot)', home: 'https://data.ecb.europa.eu', ttl: 6 * HOUR, cloud: true,
    use: 'Reference points for discount-rate and financing assumptions. Shown for information; not written into the case.',
    key: () => 'policy', load: fromSnapshot('rates'),
  },
  macro: {
    title: 'Inflation and lending rates', provider: 'World Bank Open Data', home: 'https://data.worldbank.org', ttl: 7 * DAY,
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
    title: 'Electricity grid carbon intensity (Great Britain)', provider: 'National Energy System Operator Carbon Intensity API', home: 'https://carbonintensity.org.uk', ttl: 30 * MIN,
    use: 'Well-to-wake CO₂ of battery charging for electric aircraft. Regional indicator: enter your own grid factor for other regions.',
    key: () => 'gb',
    async load(p, get) { const d = await get('https://api.carbonintensity.org.uk/intensity').then(J); const x = d.data[0]; return { from: x.from, gCO2_kWh: x.intensity.actual ?? x.intensity.forecast, index: x.intensity.index }; },
  },
  literature: {
    title: 'Latest research literature', provider: 'OpenAlex scholarly index', home: 'https://openalex.org', ttl: DAY,
    use: 'Recent and most-cited publications for each suite, to check methods and find validation data.',
    key: (p) => `${p.q}|${p.sort}`,
    async load(p, get) {
      const sort = p.sort === 'cited' ? 'cited_by_count:desc' : 'publication_date:desc', from = p.sort === 'cited' ? '' : `,from_publication_date:${new Date(Date.now() - 540 * DAY).toISOString().slice(0, 10)}`;
      const d = await get(`https://api.openalex.org/works?search=${encodeURIComponent(p.q)}&filter=type:article${from}&sort=${sort}&per-page=12&select=id,title,publication_date,doi,cited_by_count,primary_location,authorships,open_access`).then(J);
      return { total: d.meta?.count, items: (d.results || []).filter((w) => w.title).map((w) => ({ title: w.title, date: w.publication_date, cited: w.cited_by_count, url: w.doi || w.id, venue: w.primary_location?.source?.display_name || '', authors: (w.authorships || []).slice(0, 3).map((a) => a.author?.display_name).filter(Boolean).join(', ') + ((w.authorships || []).length > 3 ? ' et al.' : ''), oa: !!w.open_access?.is_oa })) };
    },
  },
  opensource: {
    title: 'Open-source solvers and tools', provider: 'GitHub repository search', home: 'https://github.com', ttl: DAY,
    use: 'Actively maintained open-source codes for the higher-fidelity models each suite hands off, and reference implementations to cross-check against.',
    key: (p) => p.q,
    async load(p, get) { const d = await get(`https://api.github.com/search/repositories?q=${encodeURIComponent(p.q)}&sort=stars&order=desc&per_page=10`, { headers: { Accept: 'application/vnd.github+json' } }).then(J); return { total: d.total_count, items: (d.items || []).map((r) => ({ name: r.full_name, desc: r.description || '', stars: r.stargazers_count, pushed: r.pushed_at, url: r.html_url, lang: r.language || '', license: r.license?.spdx_id || '' })) }; },
  },
};

/** Geocoding (not cached as a connector row: it is an interactive search). */
export async function geocode(name) {
  const d = await get(`https://geocoding-api.open-meteo.com/v1/search?name=${encodeURIComponent(name)}&count=8&language=en&format=json`).then(J);
  return (d.results || []).map((r) => ({ name: r.name, admin: [r.admin1, r.country].filter(Boolean).join(', '), lat: r.latitude, lon: r.longitude, elev_m: r.elevation, country: r.country_code, tz: r.timezone }));
}
/** Terrain elevation for up to 100 points per call. */
export async function elevations(points) {
  const out = [];
  for (let i = 0; i < points.length; i += 100) { const ch = points.slice(i, i + 100); const d = await get(`https://api.open-meteo.com/v1/elevation?latitude=${ch.map((p) => p.lat.toFixed(4)).join(',')}&longitude=${ch.map((p) => p.lon.toFixed(4)).join(',')}`).then(J); out.push(...d.elevation); }
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
  const c = CONNECTORS[id], ck = `live.${id}.${c.key(params)}`, cached = await idb.get(ck), now = Date.now();
  const fresh = cached && now - cached.ts < c.ttl;
  if (fresh && !force) return { ...cached, stale: false, fromCache: true };
  const go = () => {
    if (inflight.has(ck)) return inflight.get(ck);
    const p = (async () => {
      try {
        const net = makeNet(), data = await c.load(params, net.get), offlineCopy = net.cachedAt != null, rec = { data, ts: net.cachedAt ?? (c.cloud && Number.isFinite(data?._ts) ? data._ts : Date.now()) };
        if (offlineCopy && cached && cached.ts >= rec.ts) return { ...cached, stale: true, fromCache: true, error: 'Offline' };
        await idb.set(ck, rec); setStatus(id, { ts: rec.ts, ok: !offlineCopy, error: offlineCopy ? 'Offline' : null, key: c.key(params), params, via: c.cloud ? 'snapshot' : c.bundled && !data?.osm ? 'bundled' : 'direct' });
        return { ...rec, stale: offlineCopy, fromCache: offlineCopy, error: offlineCopy ? 'Offline' : undefined };
      }
      catch (e) { const msg = navigator.onLine === false ? 'Offline' : e.name === 'TimeoutError' ? 'Timed out' : e.message || 'Request failed'; setStatus(id, { ok: false, error: msg, tried: Date.now(), params }); return cached ? { ...cached, stale: true, fromCache: true, error: msg } : { data: null, ts: 0, stale: true, error: msg }; }
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
/** Carbon price for the case in USD/t: the EUA auction price converted with the newest exchange rate available. */
export function carbonPrice(carbon, fx) {
  if (!Number.isFinite(carbon?.eur_t) || !recent(carbon.date)) return null;
  const rate = Number.isFinite(fx?.rates?.EUR) ? fx.rates.EUR : carbon.eur_per_usd, fxDate = Number.isFinite(fx?.rates?.EUR) ? fx.date : carbon.fx_date;
  if (!Number.isFinite(rate) || rate <= 0) return null;
  return { usd_t: usdFromEur(carbon.eur_t, rate), source: `EU ETS allowance (EUA) auction ${carbon.date}: ${carbon.eur_t} EUR/t (EEX), ${Math.round(rate * 1e4) / 1e4} EUR per USD on ${fxDate}` };
}

/** Pull everything relevant for the current case site and write it into the case (with provenance). */
export async function refreshSite({ force = false } = {}) {
  const s = state.case.site; if (s.lat == null || s.lon == null) return null;
  const p = { lat: s.lat, lon: s.lon }, [w, cl, mar] = await Promise.all([live('weather', p, { force }), live('climate', p, { force }), live('marine', p, { force })]);
  const patch = {};
  if (w.data) Object.assign(patch, { T_C: w.data.T_C, p_hPa: w.data.p_hPa, rh: w.data.rh, wind_ms: w.data.wind_ms, wind_dir_deg: w.data.wind_dir_deg, gust_ms: w.data.gust_ms, precip_mm_h: w.data.precip_mm_h, cloud_pct: w.data.cloud_pct,
    visibility_m: w.data.visibility_m ?? s.visibility_m, freezing_level_m: w.data.freezing_level_m ?? s.freezing_level_m, winds_aloft: w.data.aloft, source: `Open-Meteo, ${w.data.time}Z${w.stale ? ' (cached)' : ''}`, updated: w.ts });
  if (w.data && patch.precip_mm_h > 0.1) { patch.runway_mu_brake = patch.T_C <= 0 ? 0.1 : 0.25; patch.runway_state = patch.T_C <= 0 ? 'contaminated (freezing precipitation)' : 'wet'; } else if (w.data) { patch.runway_mu_brake = 0.4; patch.runway_state = 'dry'; }
  if (cl.data?.hot99_C != null) patch.design_hot_C = cl.data.hot99_C;
  if (mar.data?.wave_height_m != null) patch.wave_height_m = mar.data.wave_height_m;
  if (state.settings.autoApplyLive && Object.keys(patch).length) patchCaseMany('site', patch);
  return { weather: w, climate: cl, marine: mar };
}
/** Global (location-independent) feeds: FX, oil, jet fuel, carbon, macro-economics, interest rates, space weather. */
export async function refreshGlobal({ force = false } = {}) {
  const [fx, oil, kp, macro, jet, carbon, rates] = await Promise.all([live('fx', {}, { force }), live('oil', {}, { force }), live('spaceweather', {}, { force }), live('macro', { country: state.case.site.country || 'WLD' }, { force }), live('jetfuel', {}, { force }), live('carbon', {}, { force }), live('rates', {}, { force })]);
  if (state.settings.autoApplyLive) {
    if (kp.data) patchCaseMany('site', { kp_index: kp.data.kp }, true);
    const e = {};
    const fp = fuelPrice(jet.data, oil.data, state.settings.fuelCrack);
    if (fp && state.case.prop.fuel !== 'Liquid hydrogen' && state.case.prop.type !== 'electric') { const k = state.case.prop.fuel === 'Avgas 100LL' ? 2.2 : state.case.prop.fuel?.startsWith('SAF') ? 2.84 : 1; e.fuel_usd_kg = Math.round(fp.usd_kg * k * 1000) / 1000; e.fuel_source = fp.source + (k !== 1 ? ` × ${k} for ${state.case.prop.fuel} (indicative ratio)` : ''); const vol = fp.quoted && Number.isFinite(jet.data.vol_annual) ? jet.data.vol_annual : oil.data?.vol_annual; if (Number.isFinite(vol)) e.fuel_vol = Math.round(vol * 1000) / 1000; }
    const cp = carbonPrice(carbon.data, fx.data);
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
  (state.settings.autoRefreshLive ? Promise.race([seedFromSnapshot().catch(() => 0), new Promise((r) => setTimeout(r, 2500))]) : Promise.resolve()).then(() => tick());
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
