// Live-data connectors. Every request goes straight from the user's own browser to a public,
// CORS-enabled provider — there is no server of ours in the middle, so data stays fresh for every user
// regardless of where the app was built or hosted. Responses are cached in IndexedDB with a time-to-live
// and served stale-while-revalidate, so the app keeps working offline with the last known values.

import { idb, state, patchCaseMany, emit, ls } from './store.js';
import { isa } from './atmosphere.js';

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
    title: 'Nearby aerodromes and runways', provider: 'OpenStreetMap via Overpass API', home: 'https://www.openstreetmap.org', ttl: 30 * DAY, needs: 'site',
    use: 'Runway length, surface, heading and elevation for take-off and landing analyses and alternates.',
    key: (p) => `${r3(p.lat)},${r3(p.lon)}`,
    async load(p, get) {
      const q = `[out:json][timeout:25];(nwr["aeroway"~"^(aerodrome|heliport)$"](around:40000,${p.lat},${p.lon});way["aeroway"="runway"](around:40000,${p.lat},${p.lon}););out tags geom 250;`;
      // public Overpass instances are volunteer-run and sometimes busy: try each in turn
      let d = null, lastErr = null;
      for (const ep of ['https://overpass-api.de/api/interpreter', 'https://overpass.kumi.systems/api/interpreter', 'https://overpass.private.coffee/api/interpreter', 'https://maps.mail.ru/osm/tools/overpass/api/interpreter']) {
        try { d = await get(`${ep}?data=${encodeURIComponent(q)}`).then(J); break; } catch (e) { lastErr = e; }
      }
      if (!d) throw new Error(lastErr?.message === 'Failed to fetch' ? 'Map service busy' : lastErr?.message || 'Map service busy');
      const hav = (a, b) => { const R = 6371008.8, f1 = (a.lat * Math.PI) / 180, f2 = (b.lat * Math.PI) / 180, df = f2 - f1, dl = ((b.lon - a.lon) * Math.PI) / 180, s = Math.sin(df / 2) ** 2 + Math.cos(f1) * Math.cos(f2) * Math.sin(dl / 2) ** 2; return 2 * R * Math.asin(Math.sqrt(s)); };
      const centre = (e) => (e.lat != null ? { lat: e.lat, lon: e.lon } : e.bounds ? { lat: (e.bounds.minlat + e.bounds.maxlat) / 2, lon: (e.bounds.minlon + e.bounds.maxlon) / 2 } : null);
      const fields = [], runways = [];
      for (const e of d.elements || []) {
        const t = e.tags || {}, c = centre(e); if (!c) continue;
        if (t.aeroway === 'runway') {
          const g = e.geometry || []; let len = Number(String(t.length || '').replace(/[^0-9.]/g, '')) || 0, hdg = null;
          if (g.length >= 2) { const a = g[0], b = g[g.length - 1]; if (!len) len = hav(a, b); hdg = (Math.atan2((b.lon - a.lon) * Math.cos((a.lat * Math.PI) / 180), b.lat - a.lat) * 180) / Math.PI; hdg = (hdg + 360) % 180; }
          runways.push({ ref: t.ref || '', len_m: Math.round(len), surface: t.surface || 'unknown', width_m: Number(t.width) || null, heading_deg: hdg == null ? null : Math.round(hdg), ...c });
        } else fields.push({ name: t.name || t['name:en'] || t.icao || 'Unnamed aerodrome', icao: t.icao || '', iata: t.iata || '', kind: t.aeroway, ele_m: t.ele != null && t.ele !== '' ? Number(t.ele) : null, dist_km: hav(p, c) / 1000, ...c });
      }
      for (const f of fields) { f.runways = runways.filter((r) => hav(f, r) < 4500).sort((a, b) => b.len_m - a.len_m); f.longest_m = f.runways[0]?.len_m || 0; }
      fields.sort((a, b) => a.dist_km - b.dist_km);
      return { fields: fields.slice(0, 25), nRunways: runways.length };
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
        const net = makeNet(), data = await c.load(params, net.get), offlineCopy = net.cachedAt != null, rec = { data, ts: net.cachedAt ?? Date.now() };
        if (offlineCopy && cached && cached.ts >= rec.ts) return { ...cached, stale: true, fromCache: true, error: 'Offline' };
        await idb.set(ck, rec); setStatus(id, { ts: rec.ts, ok: !offlineCopy, error: offlineCopy ? 'Offline' : null, key: c.key(params), params });
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
const KG_PER_BBL = 158.987 * 0.804; // litres per barrel × Jet A-1 density
export const jetFromBrent = (usdBbl, crack) => (usdBbl * (1 + crack)) / KG_PER_BBL;

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
/** Global (location-independent) feeds: FX, oil, macro-economics, space weather. */
export async function refreshGlobal({ force = false } = {}) {
  const [fx, oil, kp, macro] = await Promise.all([live('fx', {}, { force }), live('oil', {}, { force }), live('spaceweather', {}, { force }), live('macro', { country: state.case.site.country || 'WLD' }, { force })]);
  if (state.settings.autoApplyLive) {
    if (kp.data) patchCaseMany('site', { kp_index: kp.data.kp }, true);
    const e = {};
    if (oil.data && state.case.prop.fuel !== 'Liquid hydrogen' && state.case.prop.type !== 'electric') { const base = jetFromBrent(oil.data.brent_usd_bbl, state.settings.fuelCrack); e.fuel_usd_kg = Math.round(base * (state.case.prop.fuel === 'Avgas 100LL' ? 2.2 : state.case.prop.fuel?.startsWith('SAF') ? 2.5 : 1) * 1000) / 1000; e.fuel_source = `Brent ${oil.data.brent_usd_bbl} USD/bbl on ${oil.data.date}`; if (Number.isFinite(oil.data.vol_annual)) e.fuel_vol = Math.round(oil.data.vol_annual * 1000) / 1000; }
    if (macro.data?.inflation_pct != null && macro.data.inflation_pct > -5 && macro.data.inflation_pct < 60) { e.inflation = Math.round(macro.data.inflation_pct * 10) / 1000; e.inflation_source = `${macro.data.country} ${macro.data.inflation_year}`; }
    if (fx.data) { e.fx_rates = fx.data.rates; e.fx_date = fx.data.date; }
    if (Object.keys(e).length) patchCaseMany('econ', e);
  }
  return { fx, oil, kp, macro };
}

// ---- scheduling ---------------------------------------------------------------------------
let timer = null;
export function startLive() {
  const tick = (force = false) => { if (!state.settings.autoRefreshLive || navigator.onLine === false) return; refreshSite({ force }).catch(() => {}); refreshGlobal({ force }).catch(() => {}); };
  tick();
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
