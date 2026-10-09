// World airports and runways, bundled with the app (js/data/airports/, built by tools/fetch-data.mjs from
// the public-domain OurAirports datasets). Nothing here needs the network: tiles are plain files of this
// app, stored on the device by the service worker, and only the 10° × 10° tiles a query touches are read.
// Works in the browser (fetch) and in Node (file system).

export const TILE_DEG = 10;
export const TYPES = { L: 'large airport', M: 'medium airport', S: 'small airport', H: 'heliport', W: 'seaplane base', B: 'balloonport' };
/** Plain-word runway surfaces; the tile files store the index. */
export const SURFACES = ['unknown', 'asphalt', 'concrete', 'paved', 'grass', 'gravel', 'dirt', 'unpaved', 'sand', 'compacted', 'snow', 'ice', 'water', 'metal', 'wood'];
/** Common name endings, stored as "~<index in base 36>" to keep the files small. */
export const SUFFIXES = ['', ' Airport', ' Heliport', ' International Airport', ' Regional Airport', ' Municipal Airport', ' Airstrip', ' Seaplane Base', ' County Airport', ' Airfield', ' Field', ' Ranch Airport', ' Hospital Heliport', ' Air Base', ' Aerodrome', ' Landing Strip', ' Field Airport', ' Farm Airport', ' Memorial Airport', ' Airpark', ' Medical Center Heliport'];

const R_EARTH = 6371008.8, RAD = Math.PI / 180;
export const fold = (s) => String(s || '').normalize('NFD').replace(/[̀-ͯ]/g, '').toLowerCase().replace(/[^a-z0-9]+/g, ' ').trim();
export function tileOf(lat, lon) {
  const n = 180 / TILE_DEG, m = 360 / TILE_DEG, i = Math.min(n - 1, Math.max(0, Math.floor((lat + 90) / TILE_DEG))), j = Math.min(m - 1, Math.max(0, Math.floor((lon + 180) / TILE_DEG)));
  return { id: `${i}_${j}`, lat0: i * TILE_DEG - 90, lon0: j * TILE_DEG - 180 };
}
export function distanceKm(a, b) {
  const f1 = a.lat * RAD, f2 = b.lat * RAD, df = f2 - f1, dl = (b.lon - a.lon) * RAD, s = Math.sin(df / 2) ** 2 + Math.cos(f1) * Math.cos(f2) * Math.sin(dl / 2) ** 2;
  return (2 * R_EARTH * Math.asin(Math.min(1, Math.sqrt(s)))) / 1000;
}
/** Initial true bearing from a to b, degrees 0–360. */
export function bearingDeg(a, b) {
  const f1 = a.lat * RAD, f2 = b.lat * RAD, dl = (b.lon - a.lon) * RAD;
  return ((Math.atan2(Math.sin(dl) * Math.cos(f2), Math.cos(f1) * Math.sin(f2) - Math.sin(f1) * Math.cos(f2) * Math.cos(dl)) / RAD) + 360) % 360;
}

// ---- file access (browser: fetch relative to this module or the page; Node: file system) -------------
const isNode = typeof process !== 'undefined' && !!process.versions?.node && typeof window === 'undefined';
let baseUrl = null;
function base() {
  if (baseUrl) return baseUrl;
  try { if (import.meta.url) baseUrl = new URL('../data/airports/', import.meta.url); } catch { /* bundled copy: no module URL */ }
  if (!baseUrl && typeof document !== 'undefined' && /^https?:/.test(document.baseURI)) baseUrl = new URL('js/data/airports/', document.baseURI);
  if (!baseUrl) throw new Error('The bundled airport database is not available in this copy of the app');
  return baseUrl;
}
async function readText(name) {
  // the single-file copy has no data files beside it: it carries a subset inline (large and medium airports; tools/build.mjs)
  const inline = globalThis.__AEROSUITE_ESSENTIAL__?.files;
  if (inline && globalThis.__AEROSUITE_STANDALONE__) { const t = inline[`airports/${name}`]; if (typeof t !== 'string') throw new Error(`Airport database file ${name} is not part of the single-file copy`); return t; }
  const url = new URL(name, base());
  if (isNode && url.protocol === 'file:') { const fs = 'node:fs/promises', { readFile } = await import(fs); return readFile(url, 'utf8'); }
  const r = await fetch(url); if (!r.ok) throw new Error(`Airport database file ${name}: HTTP ${r.status}`);
  return r.text();
}
const memo = new Map();
const once = (key, make) => { if (!memo.has(key)) memo.set(key, make().catch((e) => { memo.delete(key); throw e; })); return memo.get(key); };

/** Dataset description: source, licence, dataset date, counts and the list of tiles. */
export const meta = () => once('index', async () => JSON.parse(await readText('index.json')));

const r4 = (v) => Math.round(v * 1e4) / 1e4, r5 = (v) => Math.round(v * 1e5) / 1e5, numOr = (s) => (s === '' || s == null ? null : Number(s));
/** Designator of the opposite runway end: 09L → 27R, 18 → 36, 04C → 22C; '' when there is none (helipads, water lanes). */
export function reciprocal(id) { const m = /^(\d{1,2})([LRC]?)$/.exec(String(id || '')); if (!m || Number(m[1]) < 1 || Number(m[1]) > 36) return ''; return String(((Number(m[1]) + 17) % 36) + 1).padStart(2, '0') + ({ L: 'R', R: 'L', C: 'C' }[m[2]] || ''); }
function parseRunway(txt, ap) {
  const f = txt.split(','), [leId, heId = reciprocal(leId)] = f[0].split('/'), flags = Number(f[4] || 0);
  const end = (ident, k, hdg) => ({ ident: ident || '', heading_degT: hdg, lat: f[k] !== '' && f[k] != null ? r5(ap.lat + Number(f[k]) / 1e5) : null, lon: f[k + 1] !== '' && f[k + 1] != null ? r5(ap.lon + Number(f[k + 1]) / 1e5) : null, elev_m: numOr(f[k + 2]), displaced_m: Number(f[k + 3] || 0) });
  let leH = numOr(f[5]), heH = numOr(f[6]);
  const le = end(leId, 7, leH), he = end(heId, 11, heH);
  let src = leH != null ? 'surveyed' : null;
  if (leH == null && heH != null) { leH = (heH + 180) % 360; src = 'surveyed'; }
  if (leH == null && le.lat != null && he.lat != null) { leH = Math.round(bearingDeg(le, he) * 10) / 10; src = 'thresholds'; }
  if (leH == null) { const n = Number(String(leId).match(/^\d{1,2}/)?.[0]); if (n >= 1 && n <= 36) { leH = n * 10; src = 'designator'; } } // magnetic, to the nearest 10°
  if (heH == null && leH != null) heH = (leH + 180) % 360;
  le.heading_degT = leH; he.heading_degT = heH;
  return { ref: [leId, heId].filter(Boolean).join('/'), len_m: Number(f[1]) || 0, width_m: numOr(f[2]), surface: SURFACES[Number(f[3] || 0)] || 'unknown', lighted: !!(flags & 1), closed: !!(flags & 2), heading_deg: leH, heading_source: src, le, he,
    lat: le.lat != null && he.lat != null ? r5((le.lat + he.lat) / 2) : ap.lat, lon: le.lon != null && he.lon != null ? r5((le.lon + he.lon) / 2) : ap.lon };
}
function parseLine(line, lat0, lon0, country, type) {
  const f = line.split('|'), ap = { ident: f[0], icao: f[1] === '1' ? f[0] : f[1] || '', iata: f[2] || '', name: (f[3] || f[0]).replace(/~([0-9a-z])$/, (_, c) => SUFFIXES[parseInt(c, 36)] || ''), type: TYPES[type] || 'aerodrome', lat: r4(lat0 + Number(f[4]) / 1e4), lon: r4(lon0 + Number(f[5]) / 1e4), elev_m: numOr(f[6]), country };
  ap.runways = f[7] ? f[7].split(';').map((t) => parseRunway(t, ap)) : [];
  ap.longest_m = ap.runways.filter((r) => !r.closed).reduce((m, r) => Math.max(m, r.len_m), 0);
  return ap;
}
const tile = (id) => once(`t${id}`, async () => {
  const idx = await meta(); if (!idx.tiles[id]) return [];
  const [i, j] = id.split('_').map(Number), lat0 = i * TILE_DEG - 90, lon0 = j * TILE_DEG - 180;
  const out = []; let cc = '', ty = '';
  for (const l of (await readText(`t${id}.txt`)).split('\n')) { if (!l) continue; if (l[0] === '#') { [cc, ty] = l.slice(1).split('|'); continue; } out.push(parseLine(l, lat0, lon0, cc, ty)); }
  return out;
});
const codes = () => once('codes', async () => { const m = new Map(); for (const l of (await readText('codes.txt')).split('\n')) { const k = l.indexOf(':'); if (k < 0) continue; const id = l.slice(0, k); for (const c of l.slice(k + 1).split(' ')) { if (!m.has(c)) m.set(c, []); m.get(c).push(id); } } return m; });
const names = () => once('names', async () => { const out = []; let ty = ''; for (const l of (await readText('names.txt')).split('\n')) { if (!l) continue; if (l[0] === '#') { ty = l.slice(1); continue; } const f = l.split('|'), name = (f[2] || '').replace(/~([0-9a-z])$/, (_, c) => SUFFIXES[parseInt(c, 36)] || ''); out.push({ ident: f[0], type: ty, tile: f[1], name, municipality: f[3] || '', text: ` ${fold(name)} ${fold(f[3])} ` }); } return out; });

/**
 * Aerodromes within radiusKm of a point, nearest first. Each: { ident, icao, iata, name, type, lat, lon, elev_m,
 * country, dist_km, longest_m, runways: [{ ref, len_m, width_m, surface, lighted, closed, heading_deg, le, he, lat, lon }] }.
 */
export async function nearest(lat, lon, radiusKm = 40, limit = 25) {
  if (!Number.isFinite(lat) || !Number.isFinite(lon)) return [];
  const p = { lat, lon }, dLat = radiusKm / 111.19, cos = Math.cos(Math.min(89.9, Math.abs(lat) + dLat) * RAD), dLon = Math.min(180, dLat / Math.max(cos, 1e-6));
  const ids = new Set();
  for (let la = Math.max(-90, lat - dLat); ; la += TILE_DEG) {
    const laC = Math.min(la, Math.min(90, lat + dLat));
    for (let lo = lon - dLon; ; lo += TILE_DEG) { const loC = Math.min(lo, lon + dLon); ids.add(tileOf(laC, ((((loC + 180) % 360) + 360) % 360) - 180).id); if (loC >= lon + dLon) break; }
    if (laC >= Math.min(90, lat + dLat)) break;
  }
  const out = [];
  for (const list of await Promise.all([...ids].map(tile))) for (const a of list) { if (Math.abs(a.lat - lat) > dLat) continue; const d = distanceKm(p, a); if (d <= radiusKm) out.push({ ...a, dist_km: d }); }
  out.sort((a, b) => a.dist_km - b.dist_km);
  return out.slice(0, limit);
}

/** Exact look-up by ICAO, IATA or local identifier (case-insensitive). Returns one aerodrome or null. */
export async function findByCode(code) {
  const c = String(code || '').trim().toUpperCase(); if (!/^[A-Z0-9-]{2,8}$/.test(c)) return null;
  const ids = (await codes()).get(c); if (!ids) return null;
  const hits = (await Promise.all(ids.map(tile))).flat().filter((a) => a.ident === c || a.icao === c || a.iata === c);
  const rank = (a) => (a.icao === c ? 0 : a.iata === c ? 1 : 2) * 10 + Object.values(TYPES).indexOf(a.type);
  return hits.sort((a, b) => rank(a) - rank(b))[0] || null;
}

/**
 * Text search over codes, airport names and the towns they serve. An exact code match comes first; names are
 * searched among airports with airline service or an IATA code, all large and medium airports, and any
 * tile already in memory. Returns full aerodrome records, best first.
 */
export async function search(text, limit = 8) {
  const q = fold(text); if (q.length < 2) return [];
  const out = [], seen = new Set(), push = (a) => { if (a && !seen.has(a.ident)) { seen.add(a.ident); out.push(a); } };
  if (/^[a-z0-9-]{2,8}$/i.test(String(text).trim())) push(await findByCode(text));
  const words = q.split(' '), match = (t) => words.every((w) => t.includes(' ' + w)), order = 'LMSWHB';
  const hits = (await names()).filter((n) => match(n.text)).map((n) => ({ n, score: order.indexOf(n.type) * 10 + (n.text.includes(` ${q} `) ? 0 : 3) })).sort((a, b) => a.score - b.score).slice(0, limit * 2);
  const tiles = await Promise.all(hits.map((h) => tile(h.n.tile)));
  hits.forEach((h, i) => push(tiles[i].find((a) => a.ident === h.n.ident)));
  if (out.length < limit) for (const [k, pr] of memo) { if (!/^t\d/.test(k)) continue; for (const a of await pr) if (out.length < limit * 2 && match(` ${fold(a.name)} `)) push(a); } // local strips near places already looked at
  return out.slice(0, limit);
}
