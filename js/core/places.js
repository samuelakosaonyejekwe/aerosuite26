// Place names, design-temperature climatology and coarse terrain bundled with the app — no network at run time:
//   js/data/places/   GeoNames cities with more than 15 000 inhabitants (CC BY 4.0), built by `tools/fetch-data.mjs --places`;
//   js/data/climate/  hot-day, cold-day and strong-wind statistics from the NCEP-DOE Reanalysis 2 (US Government, public
//                     domain), built by tools/fetch-climate.mjs;
//   js/data/terrain/  0.5° mean terrain height from NOAA ETOPO 2022 (public domain), built by `tools/fetch-data.mjs --terrain`.
// The files are plain files of this app, stored on the device by the service worker. Works in the browser and in Node.

import { fold, distanceKm } from './airports.js';

const isNode = typeof process !== 'undefined' && !!process.versions?.node && typeof window === 'undefined';
function base(sub) {
  try { if (import.meta.url) return new URL(`../data/${sub}/`, import.meta.url); } catch { /* bundled copy: no module URL */ }
  if (typeof document !== 'undefined' && /^https?:/.test(document.baseURI)) return new URL(`js/data/${sub}/`, document.baseURI);
  throw new Error('The bundled place, climate and terrain data are not available in this copy of the app');
}
async function read(sub, name, as) {
  // the single-file copy carries a subset inline (cities of 100 000 or more, the climate grid; tools/build.mjs)
  const inline = globalThis.__AEROSUITE_STANDALONE__ ? globalThis.__AEROSUITE_ESSENTIAL__?.files : null;
  if (inline) {
    const v = inline[`${sub}/${name}`]; if (v == null) throw new Error(`${sub}/${name} is not part of the single-file copy`);
    if (as === 'bin') { const b = atob(v.b64), u = new Uint8Array(b.length); for (let i = 0; i < b.length; i++) u[i] = b.charCodeAt(i); return u.buffer; }
    return as === 'json' ? JSON.parse(v) : v;
  }
  const url = new URL(name, base(sub));
  if (isNode && url.protocol === 'file:') { const fs = 'node:fs/promises', { readFile } = await import(fs), b = await readFile(url); return as === 'bin' ? b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) : as === 'json' ? JSON.parse(b.toString('utf8')) : b.toString('utf8'); }
  const r = await fetch(url); if (!r.ok) throw new Error(`data file ${sub}/${name}: HTTP ${r.status}`);
  return as === 'bin' ? r.arrayBuffer() : as === 'json' ? r.json() : r.text();
}
const memo = new Map();
const once = (key, make) => { if (!memo.has(key)) memo.set(key, make().catch((e) => { memo.delete(key); throw e; })); return memo.get(key); };

// ---- places -------------------------------------------------------------------------------------
export const placesMeta = () => once('places.index', () => read('places', 'index.json', 'json'));
function parseCities(text) {
  const out = []; let cc = '', admin = '';
  for (const l of text.split('\n')) {
    if (!l) continue; if (l[0] === '#') { [cc, admin = ''] = l.slice(1).split('|'); continue; }
    const f = l.split('|'), name = f[0], ascii = f[1] || name;
    out.push({ name, country: cc, admin, lat: Number(f[2]) / 1000, lon: Number(f[3]) / 1000, elev_m: f[4] === '' ? null : Number(f[4]), pop: Number(f[5]) * 1000, key: fold(ascii), key2: ascii === name ? '' : fold(name) });
  }
  return out;
}
const cities = () => once('places', async () => parseCities(await read('places', 'cities.txt', 'text')));
/** Second tier (5 000–15 000 inhabitants): one file per first letter, read only when a search starts with that letter. */
const smallCities = (letter) => once(`places.small.${letter}`, async () => { const m = await placesMeta().catch(() => null); if (!m?.small?.buckets?.[letter]) return []; return parseCities(await read('places', `small/${letter}.txt`, 'text')); });
let regionNames = null;
const countryName = (cc) => { try { regionNames ||= new Intl.DisplayNames(['en'], { type: 'region' }); return regionNames.of(cc) || cc; } catch { return cc; } };
/**
 * Cities matching a name, best first: an exact name, then names that start with the text, then names containing it
 * as a word; larger places first within each class. Returns place-search rows: { name, admin, lat, lon, elev_m, country, population }.
 */
export async function searchPlaces(text, limit = 8) {
  const q = fold(text); if (q.length < 2) return [];
  const hits = [];
  for (const c of await cities()) {
    const cls = (k) => (!k ? 9 : k === q ? 0 : k.startsWith(q + ' ') ? 1 : k.startsWith(q) ? 2 : k.includes(' ' + q) ? 3 : 9), rank = Math.min(cls(c.key), cls(c.key2));
    if (rank < 9) hits.push({ c, rank });
  }
  // smaller places: only names that start with the text, from the one file for its first letter (not available in the single-file copy)
  if (q.length >= 3 && hits.filter((x) => x.rank <= 2).length < limit) { const letter = /^[a-z]/.test(q) ? q[0] : '0'; for (const c of await smallCities(letter).catch(() => [])) { const k = (x) => (!x ? 9 : x === q ? 0 : x.startsWith(q + ' ') ? 1 : x.startsWith(q) ? 2 : 9), rank = Math.min(k(c.key), k(c.key2)); if (rank < 9) hits.push({ c, rank: rank + 0.5 }); } }
  hits.sort((a, b) => a.rank - b.rank || b.c.pop - a.c.pop);
  return hits.slice(0, limit).map(({ c }) => ({ name: c.name, admin: [c.admin && c.admin !== c.name ? c.admin : '', countryName(c.country)].filter(Boolean).join(', '), lat: c.lat, lon: c.lon, elev_m: c.elev_m, country: c.country, population: c.pop, source: 'GeoNames' }));
}
/** The nearest bundled city within radiusKm: { name, country, admin, lat, lon, elev_m, dist_km } or null. */
export async function nearestPlace(lat, lon, radiusKm = 25) {
  if (!Number.isFinite(lat) || !Number.isFinite(lon)) return null;
  const p = { lat, lon }, dLat = radiusKm / 111.19, dLon = dLat / Math.max(0.05, Math.cos((lat * Math.PI) / 180)); let best = null;
  for (const c of await cities()) { if (Math.abs(c.lat - lat) > dLat) continue; const dl = Math.abs(c.lon - lon); if (Math.min(dl, 360 - dl) > dLon) continue; const d = distanceKm(p, c); if (d <= radiusKm && (!best || d < best.dist_km)) best = { name: c.name, country: c.country, admin: c.admin, lat: c.lat, lon: c.lon, elev_m: c.elev_m, dist_km: d }; }
  return best;
}

// ---- terrain ------------------------------------------------------------------------------------
export const terrainMeta = () => once('terrain.index', () => read('terrain', 'index.json', 'json'));
const terrainGrid = () => once('terrain', async () => { const m = await terrainMeta(), buf = await read('terrain', m.file, 'bin'); if (buf.byteLength !== m.grid.nj * m.grid.ni * 2) throw new Error('terrain file has an unexpected size'); return { m, v: new DataView(buf) }; });
/** Mean terrain height [m] of the 0.5° cells around each point, interpolated (sea = 0). Coarse: peaks are higher. */
export async function terrain(points) {
  const { m, v } = await terrainGrid(), g = m.grid, at = (j, i) => v.getInt16(2 * (Math.min(g.nj - 1, Math.max(0, j)) * g.ni + (((i % g.ni) + g.ni) % g.ni)), true);
  return points.map((p) => { const y = (g.lat_north - g.step_deg / 2 - p.lat) / g.step_deg, x = (p.lon - g.lon_west - g.step_deg / 2) / g.step_deg, j0 = Math.floor(y), i0 = Math.floor(x), fy = y - j0, fx = x - i0; return Math.round((at(j0, i0) * (1 - fx) + at(j0, i0 + 1) * fx) * (1 - fy) + (at(j0 + 1, i0) * (1 - fx) + at(j0 + 1, i0 + 1) * fx) * fy); });
}

// ---- climate ------------------------------------------------------------------------------------
export const climateMeta = () => once('climate.index', () => read('climate', 'index.json', 'json'));
const climateGrid = () => once('climate', async () => { const m = await climateMeta(), buf = await read('climate', m.file, 'bin'); if (buf.byteLength !== m.fields.length * m.grid.nj * m.grid.ni * 2) throw new Error('climate file has an unexpected size'); return { m, v: new DataView(buf) }; });
/**
 * Design temperatures and strong-wind value at a point from the bundled reanalysis statistics. The strong-wind value
 * is the land one (coastal sea cells replaced by nearby land values) unless the point is at sea. Temperatures are moved
 * from the reanalysis terrain height to `elev_m` (when given) with the standard lapse rate, 6.5 K/km.
 */
export async function climateAt(lat, lon, elev_m = null) {
  const { m, v } = await climateGrid(), g = m.grid, n = g.nj * g.ni, y = Math.min(g.nj - 1, Math.max(0, (g.lat0 - lat) / g.step_deg)), x = (((lon - g.lon0) % 360) + 360) % 360 / g.step_deg, j0 = Math.min(g.nj - 2, Math.floor(y)), i0 = Math.floor(x) % g.ni, i1 = (i0 + 1) % g.ni, fy = y - j0, fx = x - Math.floor(x);
  const val = (k, j, i) => v.getInt16(2 * (k * n + j * g.ni + i), true) * m.fields[k].scale, out = {};
  m.fields.forEach((f, k) => { out[f.id] = (val(k, j0, i0) * (1 - fx) + val(k, j0, i1) * fx) * (1 - fy) + (val(k, j0 + 1, i0) * (1 - fx) + val(k, j0 + 1, i1) * fx) * fy; });
  // a point with (almost) no land in the reanalysis around it is a site at sea: use the unmodified fields there
  const atSea = out.land_pct != null && out.land_pct < 10; if (atSea && out.wind99_sea_ms != null) out.wind99_ms = out.wind99_sea_ms;
  const dz = Number.isFinite(elev_m) && Math.abs(elev_m - out.orog_m) < 4000 ? elev_m - out.orog_m : 0, dT = -0.0065 * dz, r1 = (t) => Math.round(t * 10) / 10;
  return { from: String(m.period.temperature[0]), to: String(m.period.temperature[1]), hot99_C: r1(out.hot99_C + dT), hot_mean_C: r1(out.hot50_C + dT), cold01_C: r1(out.cold01_C + dT), wind99_ms: r1(out.wind99_ms), model_elev_m: Math.round(out.orog_m), elevation_correction_K: r1(dT), wind_period: m.period.wind.join('–'), source: m.source, surface: atSea ? 'sea' : 'land' };
}
