// Licence and data-source compliance checks: `node tests/licences.mjs`
//  - the audited registry (js/data/licences.json) is complete and well formed;
//  - every live-data connector, snapshot feed, service-worker host and content-security-policy host has a registry entry;
//  - with "commercial": true no request goes to a host whose source is not classed commercial-ok, unless the operator
//    supplied the key or address for it — checked by running every connector and every snapshot feed against a fetch
//    that records the hosts asked for;
//  - the weather providers return the same normalised object and say what they could not supply;
//  - the required attribution statements and notices are present and shipped.
// No network is used.

import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { CONNECTORS, WEATHER, CONFIG_DEFAULTS, config, setConfig, cleanConfig, allowed, feedAllowed, attributions, weatherPlan, finishWeather, carbonMode, carbonPrice, ukEtsCarbon, operatorCarbon, fxFromEcb, geocode, elevations, osmAvailable, placeSearchAvailable, loadLicences } from '../js/core/live.js';
import { FEEDS, buildSnapshot } from '../tools/snapshot.mjs';

const root = fileURLToPath(new URL('..', import.meta.url)), text = (f) => readFileSync(root + f, 'utf8');
let fail = 0, n = 0;
const ok = (cond, msg) => { n++; if (!cond) { fail++; console.log('  FAIL ' + msg); } };
const section = (t) => console.log(t);
const near = (a, b, tol) => Math.abs(a - b) <= tol;

const reg = await loadLicences(), byId = new Map(reg.sources.map((s) => [s.id, s]));
const CLASSES = ['commercial-ok', 'commercial-with-key', 'unclear', 'not-allowed'];
const hostOwner = new Map(); for (const s of reg.sources) for (const h of s.hosts || []) hostOwner.set(h, s);

// ---- registry ------------------------------------------------------------------------------------
section('Licence registry');
ok(reg.retrieved === '2026-10-09' && reg.sources.length >= 25 && reg.components.length >= 6, 'registry header, sources and components');
ok(new Set(reg.sources.map((s) => s.id)).size === reg.sources.length, 'source ids are unique');
for (const s of reg.sources) {
  ok(typeof s.name === 'string' && typeof s.usedFor === 'string' && s.usedFor.length > 20 && typeof s.licence === 'string' && CLASSES.includes(s.class) && s.retrieved === '2026-10-09' && Array.isArray(s.obligations) && s.obligations.length > 0 && Array.isArray(s.hosts) && 'usedInCommercial' in s && 'replacement' in s && 'attribution' in s, `${s.id}: id, name, use, licence, class, retrieval date, obligations, hosts, commercial use, replacement, attribution`);
  if (!['operator-carbon', 'third-party-citations'].includes(s.id)) ok(/^https:\/\//.test(s.url || '') && typeof s.quote === 'string' && s.quote.length >= 30, `${s.id}: terms address and a quoted sentence`);
  if (s.class !== 'commercial-ok') ok(s.usedInCommercial !== true || s.needsLegalDecision === true, `${s.id}: a source that is not commercial-ok is not used unconditionally in commercial mode`);
  if (s.class === 'not-allowed' || s.class === 'unclear') ok(typeof s.replacement === 'string' && s.replacement.length > 10, `${s.id}: names its replacement`);
  if (s.class === 'commercial-ok' && s.kind !== 'citations') ok(typeof s.attribution === 'string' && s.attribution.length > 10, `${s.id}: carries an attribution text`);
}
for (const c of reg.components) ok(typeof c.name === 'string' && typeof c.version === 'string' && typeof c.licence === 'string' && typeof c.modifications === 'string' && Array.isArray(c.obligations) && typeof c.shipped === 'boolean' && (!c.shipped || existsSync(root + c.licenceFile)) && c.files.every((f) => existsSync(root + f)), `component ${c.id}: version, licence, modifications, obligations, and its files and licence text exist`);
ok(byId.get('open-meteo-free').class === 'not-allowed' && byId.get('open-meteo-customer').class === 'commercial-with-key' && byId.get('fred').class === 'not-allowed' && byId.get('fred').used === false && byId.get('eex-auction').class === 'not-allowed' && byId.get('github-api').class === 'unclear' && byId.get('overpass-public').class === 'commercial-with-key', 'restricted sources are classed as such');

// ---- configuration -------------------------------------------------------------------------------
section('Deployment configuration');
const fileCfg = JSON.parse(text('config.json'));
ok(['commercial', 'openMeteoApiKey', 'weatherProvider', 'carbonPriceUrl', 'attribution'].every((k) => k in fileCfg) && fileCfg.commercial === false && fileCfg.openMeteoApiKey === '' && fileCfg.weatherProvider === 'auto' && fileCfg.carbonPriceUrl === '' && fileCfg.attribution === true, 'config.json holds the documented keys with their defaults');
ok(Object.keys(fileCfg).every((k) => k in CONFIG_DEFAULTS) && JSON.stringify(cleanConfig(fileCfg)) === JSON.stringify(cleanConfig({})), 'config.json has no unknown keys and equals the built-in defaults');
ok(cleanConfig({ commercial: true, attribution: false }).attribution === true && cleanConfig({ attribution: false }).attribution === false, 'attribution cannot be switched off in a commercial deployment');
ok(cleanConfig({ commercial: 'yes', weatherProvider: 'x', carbonPriceUrl: 'javascript:alert(1)', overpassUrl: 'http://insecure.example/api', accept: 'all' }).commercial === false && cleanConfig({ weatherProvider: 'x' }).weatherProvider === 'auto' && cleanConfig({ carbonPriceUrl: 'javascript:alert(1)' }).carbonPriceUrl === '' && cleanConfig({ overpassUrl: 'http://insecure.example/api' }).overpassUrl === '' && cleanConfig({ accept: 'all' }).accept.length === 0, 'malformed configuration values fall back to the defaults');
for (const k of Object.keys(CONFIG_DEFAULTS)) ok(new RegExp('`' + k + '`').test(text('README.md')), `README documents the "${k}" key`);

// ---- connectors and feeds have registry entries --------------------------------------------------
section('Every connector and snapshot feed is in the registry');
const NC = {}, COM = { commercial: true }, KEYED = { commercial: true, openMeteoApiKey: 'KEY123', carbonPriceUrl: 'https://prices.operator.example/carbon.json', overpassUrl: 'https://overpass.operator.example/api/interpreter', openAlexApiKey: 'OA1' };
for (const [name, cfg] of [['non-commercial', NC], ['commercial', COM], ['commercial with keys', KEYED]]) {
  setConfig(cfg);
  for (const [id, c] of Object.entries(CONNECTORS)) {
    const ids = typeof c.reg === 'function' ? c.reg() : null;
    ok(Array.isArray(ids) && ids.every((x) => byId.has(x)), `${name}: connector ${id} names registry sources (${ids})`);
    if (cfg === NC) ok(ids?.length > 0, `non-commercial: connector ${id} has at least one source`);
    if (cfg !== NC) for (const x of ids || []) ok(await allowed(x), `${name}: connector ${id} lists only sources it may use (${x})`);
    ok(typeof c.title === 'string' && typeof c.provider === 'string' && /^https:\/\//.test(c.home) && c.ttl > 0 && typeof c.use === 'string', `${name}: connector ${id} title, provider, home, ttl, use`);
  }
}
for (const [id, f] of Object.entries(FEEDS)) ok(Array.isArray(f.registry) && f.registry.length > 0 && f.registry.every((x) => byId.has(x)), `snapshot feed ${id} names registry sources`);
const snapFile = JSON.parse(text('data/snapshot.json'));
for (const [id, f] of Object.entries(snapFile.feeds)) ok(Array.isArray(f.registry) && f.registry.length > 0 && f.registry.every((x) => byId.has(x)), `data/snapshot.json: feed ${id} records its registry sources`);
ok(typeof snapFile.commercial === 'boolean', 'data/snapshot.json records the mode it was built in');
for (const f of Object.values(snapFile.feeds)) for (const u of f.urls || []) ok(hostOwner.has(new URL(u).hostname) && f.registry.includes(hostOwner.get(new URL(u).hostname).id), `data/snapshot.json: ${new URL(u).hostname} belongs to a source the feed declares`);

// ---- hosts in the service worker and the content-security policy ---------------------------------
section('Hosts');
const swHosts = [...text('sw.js').match(/const LIVE_HOSTS = \[([\s\S]*?)\];/)[1].matchAll(/'([^']+)'/g)].map((m) => m[1]);
const csp = text('index.html').match(/connect-src 'self'([^;"]*)/)[1].trim().split(/\s+/).map((u) => u.replace(/^https:\/\//, ''));
const mirrorHosts = (JSON.parse(text('mirrors.json')).mirrors || []).map((m) => new URL(m.url).hostname);
for (const h of swHosts) ok(hostOwner.has(h), `service worker host ${h} has a registry entry`);
for (const h of csp) ok(hostOwner.has(h) || mirrorHosts.includes(h) || h === '*.blob.core.windows.net', `content-security-policy host ${h} has a registry entry`);
const usedHosts = new Set(reg.sources.filter((s) => s.used !== false && s.kind !== 'web-data').flatMap((s) => s.hosts)); // "web-data" sources are read only by the snapshot job …
usedHosts.add('raw.githubusercontent.com'); // … except the oil dataset, which a non-commercial browser reads directly
for (const h of usedHosts) ok(csp.includes(h) || fileCfg.commercial, `host ${h} used by a connector is allowed by the content-security policy`);
ok(!csp.some((h) => /fred\.stlouisfed|overpass\.kumi|eex/.test(h)) && !swHosts.some((h) => /fred\.stlouisfed|overpass\.kumi|eex/.test(h)), 'hosts no longer used are gone from the policy and the service worker');
ok(!/stlouisfed\.org/.test(text('js/core/live.js')) && !/https:\/\/fred\.stlouisfed\.org/.test(text('tools/snapshot.mjs')), 'no code path contacts FRED');

// ---- commercial mode: which hosts are asked ------------------------------------------------------
section('Commercial mode contacts only permitted hosts');
const realFetch = globalThis.fetch;
/** Run every connector (and the interactive look-ups) against a fetch that records the host and answers "unavailable", so every fall-back is tried too. */
async function hostsAsked(cfg, { snapshot = false } = {}) {
  setConfig(cfg);
  const seen = new Set(), rec = async (url) => { seen.add(new URL(String(url)).hostname); return new Response('', { status: 503 }); };
  globalThis.fetch = rec;
  const site = { lat: 51.4706, lon: -0.4619 }, params = { weather: site, climate: site, airquality: site, marine: site, aerodromes: { ...site, osm: true }, macro: { country: 'GB' }, literature: { q: 'aircraft icing', sort: 'recent' }, opensource: { q: 'aircraft icing' } };
  try {
    for (const [id, c] of Object.entries(CONNECTORS)) await c.load(params[id] || {}, rec).catch(() => null);
    await geocode('London').catch(() => null); await elevations([site]).catch(() => null);
    if (snapshot) { const s = await buildSnapshot({ previous: { feeds: { carbon: { ok: true, ts: 1, data: { eur_t: 80, date: '2026-10-08' }, registry: ['eex-auction'] } } } }); seen.snapshot = s; }
  } finally { globalThis.fetch = realFetch; }
  return seen;
}
const classOf = (h) => hostOwner.get(h)?.class;
const com = await hostsAsked(COM, { snapshot: true });
ok(com.size >= 8, `commercial mode still reaches its providers (${[...com].join(', ')})`);
for (const h of com) ok(classOf(h) === 'commercial-ok', `commercial mode asked ${h} (${hostOwner.get(h)?.id || 'not in the registry'}: ${classOf(h)})`);
for (const h of ['api.open-meteo.com', 'archive-api.open-meteo.com', 'air-quality-api.open-meteo.com', 'marine-api.open-meteo.com', 'geocoding-api.open-meteo.com', 'customer-api.open-meteo.com', 'overpass-api.de', 'overpass.private.coffee', 'maps.mail.ru', 'api.github.com', 'raw.githubusercontent.com', 'public.eex-group.com', 'fred.stlouisfed.org', 'power.larc.nasa.gov']) ok(!com.has(h), `commercial mode does not contact ${h}`);
for (const h of ['api.met.no', 'api.weather.gov', 'www.gov.uk', 'www.eia.gov', 'www.federalreserve.gov', 'markets.newyorkfed.org', 'data-api.ecb.europa.eu', 'api.frankfurter.dev']) ok(com.has(h), `commercial mode uses ${h}`);
ok(com.snapshot.commercial === true && com.snapshot.feeds.carbon.data === null && com.snapshot.feeds.carbon.registry.join() === 'uk-ets-price', 'commercial snapshot does not carry over an EEX value left by an earlier non-commercial run');

const keyed = await hostsAsked(KEYED, { snapshot: true });
for (const h of keyed) ok(classOf(h) === 'commercial-ok' || /^customer-[a-z-]*api\.open-meteo\.com$/.test(h) || h === 'prices.operator.example' || h === 'overpass.operator.example', `commercial mode with keys asked ${h}`);
for (const h of ['customer-api.open-meteo.com', 'customer-archive-api.open-meteo.com', 'customer-air-quality-api.open-meteo.com', 'customer-marine-api.open-meteo.com', 'customer-geocoding-api.open-meteo.com', 'prices.operator.example', 'overpass.operator.example']) ok(keyed.has(h), `with the operator's key or address, ${h} is used`);
ok(!keyed.has('api.open-meteo.com') && !keyed.has('overpass-api.de') && !keyed.has('public.eex-group.com') && !keyed.has('api.github.com'), 'keys do not unlock the free or unlicensed hosts');

const accepted = await hostsAsked({ commercial: true, accept: ['nasa-power'] });
ok(accepted.has('power.larc.nasa.gov') && !accepted.has('api.open-meteo.com'), 'a source the operator accepts ("accept") is used; nothing else changes');

const non = await hostsAsked(NC, { snapshot: true });
for (const h of ['api.open-meteo.com', 'archive-api.open-meteo.com', 'geocoding-api.open-meteo.com', 'overpass-api.de', 'api.github.com', 'raw.githubusercontent.com']) ok(non.has(h), `control: the non-commercial configuration does contact ${h}`);
ok(!non.has('public.eex-group.com'), 'the EEX auction report is not contacted in any default configuration (its terms require written approval)');
ok(!non.has('fred.stlouisfed.org'), 'FRED is not contacted in any mode');
for (const h of non) ok(hostOwner.has(h), `non-commercial mode: ${h} has a registry entry`);

setConfig(COM);
ok(carbonMode() === 'ukets' && !osmAvailable() && !placeSearchAvailable() && weatherPlan().join() === 'met-norway,nws', 'commercial defaults: UK ETS carbon price, no public Overpass, no online geocoder, MET Norway then NWS');
ok((await geocode('London').then(() => 'ok', (e) => (e.policy ? 'policy' : 'other'))) === 'policy', 'place-name search is refused by policy, not by a network error');
{ const e = await elevations([{ lat: 51.4706, lon: -0.4619 }, { lat: 0, lon: -30 }]); ok(near(e[0], 25, 15) && e[1] === null, `terrain elevation falls back to the bundled airport database (${e})`); }
{ const o = await CONNECTORS.opensource.load({ q: 'aircraft icing' }, async () => { throw new Error('no network'); }); ok(o.linkOnly && /^https:\/\/github\.com\/search\?q=/.test(o.items[0].url), 'open-source tools: a link to GitHub search instead of an API call'); }
{ const a = await CONNECTORS.aerodromes.load({ lat: 6.5774, lon: 3.3212, osm: true }, async () => { throw new Error('no network'); }); ok(a.source === 'OurAirports' && a.osm === false && a.fields[0].icao === 'DNMM', 'aerodromes: bundled database only'); }
ok(!(await feedAllowed({ data: {}, registry: ['eex-auction'] })) && !(await feedAllowed({ data: {} })) && (await feedAllowed({ data: {}, registry: ['eia'] })) && !(await allowed('no-such-source')), 'snapshot feeds from unlicensed or undeclared sources are refused');
ok(weatherPlan({ ...CONFIG_DEFAULTS, commercial: true, weatherProvider: 'open-meteo' }).join() === 'met-norway,nws' && weatherPlan({ ...CONFIG_DEFAULTS, commercial: true, openMeteoApiKey: 'k' })[0] === 'open-meteo' && weatherPlan({ ...CONFIG_DEFAULTS, weatherProvider: 'none' }).length === 0 && weatherPlan({ ...CONFIG_DEFAULTS })[0] === 'open-meteo' && weatherPlan({ ...CONFIG_DEFAULTS, commercial: true, weatherProvider: 'nws' })[0] === 'nws', 'weather provider order follows the configuration');
ok(CONNECTORS.weather.key({ lat: 1, lon: 2 }) !== (setConfig(NC), CONNECTORS.weather.key({ lat: 1, lon: 2 })), 'values cached under non-commercial terms are not reused by a commercial deployment');

// ---- provider normalisation ----------------------------------------------------------------------
section('Weather providers return the same object');
const json = (o) => async () => new Response(JSON.stringify(o), { status: 200 });
const FIELDS = ['time', 'elev_m', 'T_C', 'rh', 'p_hPa', 'qnh_hPa', 'wind_ms', 'wind_dir_deg', 'gust_ms', 'precip_mm_h', 'cloud_pct', 'code', 'visibility_m', 'freezing_level_m', 'aloft', 'provider', 'provider_label', 'missing', 'note'];
setConfig(COM);
const soon = new Date(Date.now() + 20 * 60e3).toISOString(), before = new Date(Date.now() - 40 * 60e3).toISOString();
const met = finishWeather(await WEATHER['met-norway'].load({ lat: 6.57741234, lon: 3.3212 }, async (u) => { ok(/lat=6\.5774&lon=3\.3212$/.test(u), 'MET Norway request truncates coordinates to four decimals'); return json({ geometry: { coordinates: [3.3212, 6.5774, 500] }, properties: { timeseries: [{ time: before, data: { instant: { details: { air_pressure_at_sea_level: 1013.2, air_temperature: 24, cloud_area_fraction: 79.7, relative_humidity: 95, wind_from_direction: 299.1, wind_speed: 0.9 } }, next_1_hours: { details: { precipitation_amount: 0.4 } } } }, { time: soon, data: { instant: { details: { air_temperature: 30 } } } }] } })(); }), 'met-norway');
ok(FIELDS.every((k) => k in met), 'MET Norway: every normalised field is present');
ok(met.T_C === 24 && met.rh === 0.95 && met.qnh_hPa === 1013.2 && met.wind_ms === 0.9 && met.precip_mm_h === 0.4 && met.elev_m === 500 && near(met.p_hPa, 957, 4), `MET Norway: values, and station pressure derived from sea-level pressure at 500 m (${met.p_hPa} hPa)`);
ok(met.gust_ms === null && met.visibility_m === null && met.aloft.length === 0 && met.freezing_level_estimated && near(met.freezing_level_m, 500 + 24 / 0.0065, 2) && ['gusts', 'visibility', 'freezing level', 'winds aloft'].every((m) => met.missing.includes(m)) && /MET Norway does not supply/.test(met.note) && /ISA/.test(met.note), 'MET Norway: missing fields are named, with the ISA stand-ins and a note');
const nwsReplies = { points: { properties: { observationStations: 'https://api.weather.gov/gridpoints/OKX/42,40/stations' } }, stations: { features: [{ properties: { stationIdentifier: 'KJFK' } }] }, latest: { properties: { timestamp: '2026-10-09T01:25:00+00:00', elevation: { value: 7 }, temperature: { value: 20 }, relativeHumidity: { value: 37.35 }, windDirection: { value: 270 }, windSpeed: { unitCode: 'wmoUnit:km_h-1', value: 14.832 }, windGust: { unitCode: 'wmoUnit:km_h-1', value: null }, barometricPressure: { value: 101625.52 }, seaLevelPressure: { value: null }, visibility: { value: 16093.44 }, cloudLayers: [{ amount: 'FEW' }, { amount: 'BKN' }] } } };
const nws = finishWeather(await WEATHER.nws.load({ lat: 40.6413, lon: -73.7781 }, async (u) => json(/\/points\//.test(u) ? nwsReplies.points : /observations\/latest$/.test(u) ? nwsReplies.latest : nwsReplies.stations)()), 'nws');
ok(FIELDS.every((k) => k in nws) && nws.time === '2026-10-09T01:25' && nws.T_C === 20 && near(nws.wind_ms, 4.12, 0.01) && near(nws.qnh_hPa, 1016.26, 0.01) && near(nws.p_hPa, 1015.4, 0.3) && nws.cloud_pct === 75 && near(nws.visibility_m, 16093, 1) && nws.gust_ms === null && /station KJFK/.test(nws.provider_label) && nws.missing.includes('winds aloft') && !nws.missing.includes('visibility'), 'US NWS: latest observation normalised (km/h to m/s, Pa to hPa, cloud layers)');
ok((await WEATHER.nws.load({ lat: 51, lon: 0 }, async () => new Response('', { status: 404 })).then(() => '', (e) => e.message)).includes('United States'), 'US NWS: a place outside its coverage is reported as such');
setConfig(NC);
const omReply = { elevation: 38, current: { time: '2026-10-09T01:15', temperature_2m: 25, relative_humidity_2m: 80, surface_pressure: 1008, pressure_msl: 1012, wind_speed_10m: 3, wind_direction_10m: 200, wind_gusts_10m: 6, precipitation: 0, cloud_cover: 50, weather_code: 2 }, hourly: { visibility: [9000], freezing_level_height: [4500], temperature_500hPa: [-6], wind_speed_500hPa: [12], wind_direction_500hPa: [250], geopotential_height_500hPa: [5800] } };
const om = await CONNECTORS.weather.load({ lat: 6.5, lon: 3.3 }, async (u) => { ok(/^https:\/\/api\.open-meteo\.com\/v1\/forecast\?/.test(u) && !/apikey/.test(u), 'non-commercial: free Open-Meteo host, no key'); return json(omReply)(); });
ok(FIELDS.every((k) => k in om) && om.provider === 'open-meteo' && om.missing.length === 0 && om.note === '' && om.aloft.length === 1 && om.gust_ms === 6 && om.p_hPa === 1008, 'Open-Meteo: complete object, nothing missing');
setConfig(KEYED);
await CONNECTORS.weather.load({ lat: 6.5, lon: 3.3 }, async (u) => { ok(/^https:\/\/customer-api\.open-meteo\.com\/v1\/forecast\?/.test(u) && /&apikey=KEY123$/.test(u), 'commercial with a key: customer host with the key'); return json(omReply)(); });
setConfig(COM);
const wx = await CONNECTORS.weather.load({ lat: 40.64, lon: -73.78 }, async (u) => (/api\.met\.no/.test(u) ? new Response('', { status: 503 }) : json(/\/points\//.test(u) ? nwsReplies.points : /observations\/latest$/.test(u) ? nwsReplies.latest : nwsReplies.stations)()));
ok(wx.provider === 'nws', 'commercial: when MET Norway does not answer, the next licensed provider is used');

section('Other feeds');
const fx = await fxFromEcb(async () => new Response('KEY,FREQ,CURRENCY,CURRENCY_DENOM,EXR_TYPE,EXR_SUFFIX,TIME_PERIOD,OBS_VALUE\nEXR.D.ARS.EUR.SP00.A,D,ARS,EUR,SP00,A,2020-10-30,91.5953\nEXR.D.GBP.EUR.SP00.A,D,GBP,EUR,SP00,A,2026-10-08,0.84698\nEXR.D.USD.EUR.SP00.A,D,USD,EUR,SP00,A,2026-10-08,1.1186\n'));
ok(fx.date === '2026-10-08' && fx.rates.USD === 1 && near(fx.rates.EUR, 1 / 1.1186, 1e-9) && near(fx.rates.GBP, 0.84698 / 1.1186, 1e-9) && !('ARS' in fx.rates) && fx.eur_rates.USD === 1.1186 && /calculated/.test(fx.derived), 'ECB reference rates: unchanged ECB values kept, per-USD rates derived, discontinued currencies dropped');
const uk = await ukEtsCarbon(async (u) => json(/determinations-of-the-uk-ets-carbon-price$/.test(u) ? { details: { attachments: [{ title: 'UK ETS: Carbon price for use in civil penalties, 2025', url: '/government/publications/determinations-of-the-uk-ets-carbon-price/x-2025' }, { title: 'UK ETS: Carbon price for use in civil penalties, 2026', url: '/government/publications/determinations-of-the-uk-ets-carbon-price/x-2026' }, { title: 'UK ETS: Carbon prices for use in civil penalties, 2021 and 2022', url: '/government/publications/determinations-of-the-uk-ets-carbon-price/x-2022' }] } } : (ok(/x-2026$/.test(u), 'UK ETS: the newest determination is read'), { first_published_at: '2025-11-28T17:00:00Z', details: { body: '<p>the carbon price per tonne of carbon dioxide equivalent for the scheme year beginning on 1 January 2026 is &pound;49.41.</p>' } }))());
ok(uk.price === 49.41 && uk.currency === 'GBP' && uk.scheme_year === 2026 && uk.valid_to === '2026-12-31' && uk.date === '2025-11-28' && /not an EU ETS price/.test(uk.instrument) && uk.registry[0] === 'uk-ets-price' && /Open Government Licence v3\.0/.test(uk.attribution), 'UK ETS determination: price, currency, scheme year, validity, label and licence statement');
{ const day = 86400e3, iso = (t) => new Date(t).toISOString().slice(0, 10), cur = { ...uk, valid_to: iso(Date.now() + 30 * day) }, cp = carbonPrice(cur, { rates: { GBP: 0.75 }, date: '2026-10-08' });
  ok(near(cp.usd_t, 49.41 / 0.75, 1e-9) && /UK ETS carbon price for scheme year 2026/.test(cp.source) && /GBP per USD/.test(cp.source) && !/EUA|EEX/.test(cp.source), 'UK ETS price converts with the GBP rate and is labelled as UK ETS');
  ok(carbonPrice({ ...uk, valid_to: iso(Date.now() - 60 * day) }, { rates: { GBP: 0.75 } }) === null && carbonPrice(cur, null) === null, 'an expired determination, or one without an exchange rate, is not written into the case'); }
const op = await operatorCarbon(json({ date: '2026-10-08', price: 85.07, currency: 'eur', market: 'EU ETS (EUA Dec-26)', source: 'licensed from Vendor X', attribution: '© Vendor X' }), 'https://prices.operator.example/carbon.json');
ok(op.price === 85.07 && op.currency === 'EUR' && op.operator === true && op.registry[0] === 'operator-carbon', 'operator-supplied carbon price is read');
{ const cp = carbonPrice({ ...op, date: new Date().toISOString().slice(0, 10) }, { rates: { EUR: 0.9 }, date: '2026-10-08' }); ok(near(cp.usd_t, 85.07 / 0.9, 1e-9) && /EU ETS \(EUA Dec-26\)/.test(cp.source) && /licensed from Vendor X/.test(cp.source), 'operator price converts and keeps its market and source labels'); }
ok((await operatorCarbon(json({ price: 'x' }), 'https://prices.operator.example/c.json').then(() => 'ok', () => 'refused')) === 'refused', 'a malformed operator file is refused');

// ---- attribution and notices ---------------------------------------------------------------------
section('Attribution and notices');
const attr = (id) => byId.get(id).attribution || '';
ok(/© OpenStreetMap contributors/.test(attr('osm-data')) && /openstreetmap\.org\/copyright/.test(byId.get('osm-data').attributionUrl), 'OpenStreetMap attribution with the copyright link');
ok(/World Bank/.test(attr('world-bank')) && /CC BY 4\.0/.test(attr('world-bank')), 'World Bank CC BY 4.0 attribution');
ok(/MET Norway/.test(attr('met-norway')) && /CC BY 4\.0/.test(attr('met-norway')), 'MET Norway attribution');
ok(/Source: ECB statistics/.test(attr('ecb-statistics')) && /own calculations/.test(attr('ecb-statistics')), 'ECB source statement, with derived rates labelled');
ok(/Reproduction is authorised, provided the source is acknowledged/.test(attr('icao-edb')) && /Reproduction is authorised, provided the source is acknowledged/.test(attr('easa-noise')), 'EASA acknowledgement');
ok(/Contains public sector information licensed under the Open Government Licence v3\.0\./.test(attr('uk-ets-price')), 'Open Government Licence statement');
ok(/Open-Meteo\.com/.test(attr('open-meteo-customer')) && /U\.S\. Energy Information Administration/.test(attr('eia')) && /Federal Reserve/.test(attr('frb-h15')) && /OpenAlex/.test(attr('openalex')) && /OurAirports/.test(attr('ourairports')), 'Open-Meteo, EIA, Federal Reserve, OpenAlex and OurAirports credits');
setConfig(COM);
{ const a = await attributions(CONNECTORS.weather.reg()); ok(a.length === 2 && a.some((x) => /MET Norway/.test(x.text) && /^https:/.test(x.url)), 'attribution lines resolve for the connectors in use'); }
const notices = text('THIRD_PARTY_NOTICES.md'), licence = text('LICENSE'), views = text('js/ui/views/livehub.js') + text('js/ui/views/about.js'), caseView = text('js/ui/views/case.js');
ok(/Copyright \(c\) 2026 Samuel Akosa Onyejekwe\. All rights reserved\./.test(licence) && /proprietary/.test(licence) && /THIRD_PARTY_NOTICES\.md/.test(licence) && /Lesser General Public License/.test(licence) && !/Permission is hereby granted/.test(licence), 'LICENSE: proprietary, no open-source grant, third-party components under their own licences');
for (const c of reg.components) ok(notices.includes(c.name.split(' (')[0].split(',')[0]) && (!c.shipped || notices.includes(c.licenceFile)), `THIRD_PARTY_NOTICES.md lists ${c.id}${c.shipped ? ' and where its licence text is' : ''}`);
ok(/0\.0\.23/.test(notices) && /0\.0\.7/.test(notices) && /var process;/.test(notices) && /export default occtimportjs;/.test(notices) && /export default createLazPerf;/.test(notices) && /Replacing the library/.test(notices) && /Open CASCADE exception/.test(notices) && /written offer/.test(notices), 'THIRD_PARTY_NOTICES.md: versions, modifications, the Open CASCADE exception, source offer and how to replace the LGPL kernel');
ok(/^\/\/ occt-import-js 0\.0\.23 \(LGPL-2\.1/.test(text('js/vendor/occt/occt-import-js.js')) && /Modified from/.test(text('js/vendor/occt/occt-import-js.js').slice(0, 400)) && /Modified from/.test(text('js/vendor/lazperf/laz-perf.js').slice(0, 300)), 'modified vendored files state the change in their header');
for (const t of ['© OpenStreetMap contributors', 'World Bank, World Development Indicators (CC BY 4.0)', 'MET Norway', 'ECB statistics', 'Reproduction is authorised, provided the source is acknowledged', 'Open Government Licence v3.0']) ok(notices.includes(t), `THIRD_PARTY_NOTICES.md carries “${t}”`);
ok(/licenceSection/.test(text('js/ui/views/about.js')) && /export async function licenceSection/.test(views) && /Data sources and licences/.test(views) && /dataTable\(/.test(views) && /attributionLine\(/.test(views), 'the Live data and About pages show the registry tables and attribution lines');
ok(/© OpenStreetMap contributors/.test(caseView) && /openstreetmap\.org\/copyright/.test(caseView) && /OurAirports \(public domain\)/.test(caseView) && /attributions\(/.test(caseView) && /wx_note/.test(caseView), 'the location page credits its data where it is shown and notes what the weather provider lacks');
const assets = JSON.parse(text('sw.js').match(/const ASSETS = (\[[^\n]*\]);/)[1]);
for (const f of ['LICENSE', 'THIRD_PARTY_NOTICES.md', 'config.json', 'js/data/licences.json', ...reg.components.filter((c) => c.shipped).map((c) => c.licenceFile)]) ok(assets.includes(f), `${f} is in the deployed asset list`);
const pages = text('.github/workflows/pages.yml');
ok(/LICENSE/.test(pages) && /THIRD_PARTY_NOTICES\.md/.test(pages) && /config\.json/.test(pages), 'the publish workflow copies LICENSE, THIRD_PARTY_NOTICES.md and config.json into the site');
ok(/^UNLICENSED$/.test(JSON.parse(text('package.json')).license), 'package.json does not declare an open-source licence');

console.log(`${n} checks, ${fail} failed`);
process.exit(fail ? 1 : 0);
