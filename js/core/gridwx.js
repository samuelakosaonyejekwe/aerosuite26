// Reader for the global forecast grids the snapshot job publishes under data/grid/ (tools/grid.mjs): NOAA GFS
// weather (surface and eleven pressure levels), DWD wave model sea state, NOAA GEFS-Aerosols air quality. The
// files are plain static files of the site, kept by the service worker for offline use; this module interpolates
// them bilinearly in space and linearly in time to a point. The files are difference-coded and gzipped; the browser's
// own gzip decompression (DecompressionStream) unpacks them. Works in the browser (fetch) and in Node (file system).

const isNode = typeof process !== 'undefined' && !!process.versions?.node && typeof window === 'undefined';
const H = 3600e3, RAD = Math.PI / 180;
/** A slice may be used this long after its valid time when no later slice exists (then the grid counts as out of date). */
export const MAX_AGE_H = 6;

let dirOverride = null, clock = () => Date.now();
/** Tests: a fixed clock, so stored grids can be read whatever the date. */
export function setGridClock(fn) { clock = fn || (() => Date.now()); }
/** Tests: read the grid files from another directory (a file: URL ending in "/"). */
export function setGridDir(url) { dirOverride = url; memo.clear(); }
function dir() {
  if (dirOverride) return dirOverride;
  if (typeof document !== 'undefined' && /^https?:/.test(document.baseURI)) return new URL('data/grid/', document.baseURI);
  if (isNode) { try { if (import.meta.url) return new URL('../../data/grid/', import.meta.url); } catch { /* bundled copy */ } }
  throw new Error('The forecast grids need the app to be served from a web address');
}
async function read(name, as) {
  const url = new URL(name, dir());
  if (isNode && url.protocol === 'file:') { const fs = 'node:fs/promises', { readFile } = await import(fs), b = await readFile(url); return as === 'json' ? JSON.parse(b.toString('utf8')) : b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength); }
  const r = await fetch(url, { cache: 'no-store' }); if (!r.ok) throw new Error(`forecast grid file ${name}: HTTP ${r.status}`);
  return as === 'json' ? r.json() : r.arrayBuffer();
}
const memo = new Map();
const once = (key, make) => { if (!memo.has(key)) memo.set(key, make().catch((e) => { memo.delete(key); throw e; })); return memo.get(key); };
let indexAt = 0;
/** data/grid/index.json (re-read at most every ten minutes). */
export function gridIndex({ force = false } = {}) {
  if (force || Date.now() - indexAt > 10 * 60e3) { memo.delete('index'); indexAt = Date.now(); }
  return once('index', async () => { const idx = await read('index.json', 'json'); if (![1, 2].includes(idx?.schema) || !idx.products) throw new Error('forecast grid index not recognised'); for (const k of [...memo.keys()]) if (k !== 'index' && !k.endsWith(`@${idx.generated}`)) memo.delete(k); return idx; });
}

const SIZE = { i8: 1, u8: 1, i16: 2 };
async function gunzip(buf) {
  if (isNode) { const z = 'node:zlib', { gunzipSync } = await import(z), b = gunzipSync(new Uint8Array(buf)); return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength); }
  if (typeof DecompressionStream === 'undefined') throw new Error('this browser cannot unpack the forecast grids (no gzip decompression)');
  return new Response(new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
}
/** Undo the "delta-gzip" file encoding (tools/grid.mjs, encodeSlice): gunzip, then add up the differences field by field. */
async function decodeSlice(buf, fields, grids) {
  const head = new Uint8Array(buf, 0, 2), d = new Uint8Array(head[0] === 0x1f && head[1] === 0x8b ? await gunzip(buf) : buf), out = new Uint8Array(d.length), view = new DataView(out.buffer); let o = 0;
  for (const f of fields) {
    const g = grids[f.grid], n = g.nj * g.ni, wide = f.type === 'i16', prev = new Int32Array(g.ni); let left = 0;
    for (let j = 0; j < g.nj; j++) for (let i = 0; i < g.ni; i++) { const k = j * g.ni + i, p = i ? left : j ? prev[0] : 0, v = wide ? ((d[o + k] | (d[o + n + k] << 8)) + p) & 0xffff : (d[o + k] + p) & 255; if (wide) view.setUint16(o + 2 * k, v, true); else out[o + k] = v; left = v; if (!i) prev[0] = v; }
    o += (wide ? 2 : 1) * n;
  }
  if (o !== d.length) throw new Error('forecast grid file has an unexpected size');
  return out.buffer;
}
/** Read one grid file as stored (plain, or "delta-gzip"). */
const slice = async (t, fields, grids) => { const g = t.grids ? { ...grids, ...t.grids } : grids, raw = await read(t.file, 'bin'); return openSlice(t.enc === 'delta-gzip' ? await decodeSlice(raw, fields, g) : raw, fields, g); };
/** One slice file → { fieldId: { f, g, view, at } } with byte offsets worked out from the field list. */
function openSlice(buf, fields, grids) {
  const view = new DataView(buf), out = {}; let o = 0;
  for (const f of fields) { const g = grids[f.grid]; out[f.id] = { f, g, at: o }; o += SIZE[f.type] * g.nj * g.ni; }
  if (o !== buf.byteLength) throw new Error('forecast grid file has an unexpected size');
  const raw = (e, k) => (e.f.type === 'i16' ? view.getInt16(e.at + 2 * k, true) : e.f.type === 'i8' ? view.getInt8(e.at + k) : view.getUint8(e.at + k));
  const val = (e, j, i) => { const r = raw(e, j * e.g.ni + i); if (e.f.miss != null && r === e.f.miss) return NaN; let v = r * e.f.scale; if (e.f.sq) v *= v; return v + (e.f.offset || 0); };
  /** Bilinear value at a point; missing neighbours are left out (NaN when all four are missing). */
  const sample = (id, lat, lon) => {
    const e = out[id]; if (!e) return NaN;
    const g = e.g, y = Math.min(g.nj - 1, Math.max(0, (90 - lat) / g.step)), x = ((((lon % 360) + 360) % 360) / g.step), j0 = Math.min(g.nj - 2, Math.floor(y)), i0 = Math.floor(x) % g.ni, i1 = (i0 + 1) % g.ni, fy = y - j0, fx = x - Math.floor(x);
    let s = 0, w = 0; for (const [v, k] of [[val(e, j0, i0), (1 - fy) * (1 - fx)], [val(e, j0, i1), (1 - fy) * fx], [val(e, j0 + 1, i0), fy * (1 - fx)], [val(e, j0 + 1, i1), fy * fx]]) if (v === v && k > 0) { s += v * k; w += k; }
    return w > 0 ? s / w : NaN;
  };
  /** Nearest non-missing value within `reach` cells (for fields with gaps, e.g. sea state next to a coast). */
  const nearest = (id, lat, lon, reach = 2) => {
    const e = out[id]; if (!e) return NaN; const g = e.g, j = Math.round((90 - lat) / g.step), i = Math.round((((lon % 360) + 360) % 360) / g.step);
    let best = NaN, bd = Infinity;
    for (let dj = -reach; dj <= reach; dj++) for (let di = -reach; di <= reach; di++) { const jj = j + dj; if (jj < 0 || jj >= g.nj) continue; const v = val(e, jj, (((i + di) % g.ni) + g.ni) % g.ni); if (v !== v) continue; const d = dj * dj + (di * Math.cos(lat * RAD)) ** 2; if (d < bd) { bd = d; best = v; } }
    return best;
  };
  return { sample, nearest };
}
/** The product's slices around `now`: { a, b, w } (value = a·(1−w) + b·w), the valid time used, and whether it is extrapolated. */
async function bracket(product, now) {
  const idx = await gridIndex(), p = idx.products?.[product];
  if (!p?.times?.length) throw new Error(p?.error ? `forecast grid not available (${p.error})` : 'forecast grid not available');
  const ts = p.times.map((t) => ({ ...t, ms: Date.parse(t.valid) })).sort((x, y) => x.ms - y.ms), last = ts[ts.length - 1];
  if (now > last.ms + MAX_AGE_H * H) throw new Error(`forecast grid is out of date (last valid time ${last.valid.slice(0, 16)}Z)`);
  let i = ts.findIndex((t) => t.ms > now); if (i < 0) i = ts.length; const lo = ts[Math.max(0, i - 1)], hi = ts[Math.min(ts.length - 1, i)];
  const open = (t) => once(`${product}/${t.file}@${idx.generated}`, () => slice(t, p.fields, p.grids));
  const w = hi.ms > lo.ms ? Math.min(1, Math.max(0, (now - lo.ms) / (hi.ms - lo.ms))) : 0;
  return { p, idx, a: await open(lo), b: await open(hi), w, res: lo.resolution || null, time: new Date(Math.min(Math.max(now, ts[0].ms), last.ms)).toISOString().slice(0, 16), cycle: p.cycle };
}
const mix = (x, y, w) => (x !== x ? y : y !== y ? x : x * (1 - w) + y * w);
const r1 = (v, d = 1) => (v === v ? Math.round(v * 10 ** d) / 10 ** d : null);
const windOf = (u, v) => ({ speed: Math.hypot(u, v), dir: (Math.atan2(-u, -v) / RAD + 360) % 360 });
const stationPressure = (qnh_hPa, elev_m, T_C) => qnh_hPa * Math.exp((-9.80665 * elev_m) / (287.05 * (T_C + 273.15 + 0.00325 * elev_m)));

/**
 * Weather at a point from the NOAA GFS grid, as the normalised object the weather providers return. `elev_m` (the
 * site elevation, when known) is used to correct the 2 m temperature from the model terrain height with the
 * standard lapse rate and to derive station pressure; without it the model terrain height stands in.
 */
export async function gridWeather(lat, lon, { elev_m = null, now = clock() } = {}) {
  const { p, idx, a, b, w, time, cycle, res } = await bracket('wx', now), s = (id) => mix(a.sample(id, lat, lon), b.sample(id, lat, lon), w);
  const orogSlice = await once(`orog@${idx.generated}`, () => slice(p.static.orog, [{ id: 'orog', grid: p.static.orog.grid, type: p.static.orog.type, scale: p.static.orog.scale }], p.grids)), orog = orogSlice.sample('orog', lat, lon);
  const site = Number.isFinite(elev_m) ? elev_m : orog, dz = site - orog, T = s('t2m_C') - (Math.abs(dz) < 3000 ? 0.0065 * dz : 0), qnh = s('qnh_hPa'), sfc = windOf(s('u10_ms'), s('v10_ms'));
  const precip = s('precip_mm_h'), cloud = s('cloud_pct'), vis = s('vis_m'), fz = s('fz_m');
  const aloft = p.levels_hPa.map((l) => { const wd = windOf(s(`u_${l}`), s(`v_${l}`)); return { hPa: l, alt_m: Math.round(s(`z_${l}`)), T_C: r1(s(`T_${l}`)), speed_ms: r1(wd.speed), dir_deg: Math.round(wd.dir) % 360 }; }).filter((x) => Number.isFinite(x.alt_m) && x.alt_m > site + 30); // levels below the ground are left out
  const code = precip > 0.1 ? (T <= 0 ? 71 : precip > 4 ? 63 : 61) : vis < 1000 ? 45 : cloud > 87 ? 3 : cloud > 50 ? 2 : cloud > 12 ? 1 : 0; // WMO-style code derived from the fields above
  return { time, elev_m: Number.isFinite(elev_m) ? elev_m : Math.round(orog), T_C: r1(T), rh: r1(Math.min(100, Math.max(0, s('rh_pct'))) / 100, 2), p_hPa: r1(stationPressure(qnh, site, T)), qnh_hPa: r1(qnh), wind_ms: r1(sfc.speed), wind_dir_deg: Math.round(sfc.dir) % 360,
    gust_ms: r1(Math.max(s('gust_ms'), sfc.speed)), precip_mm_h: r1(precip, 2), cloud_pct: Math.round(cloud), code, visibility_m: Math.round(vis / 100) * 100, freezing_level_m: Math.round(Math.max(0, fz)), aloft,
    model_elev_m: Math.round(orog), cycle, grid: res || p.resolution };
}
/** Sea state at a point from the wave grid (null wave height when none of the surrounding 2° cells holds sea: inland). */
export async function gridMarine(lat, lon, { now = clock() } = {}) {
  const { a, b, w, time, cycle, p } = await bracket('sea', now), s = (id) => mix(a.sample(id, lat, lon), b.sample(id, lat, lon), w); // a value needs sea in one of the four 2° cells around the point
  const hs = s('hs_m'), dA = a.nearest('dir_deg', lat, lon, 1), dB = b.nearest('dir_deg', lat, lon, 1), dir = dA !== dA ? dB : dB !== dB ? dA : (Math.atan2(Math.sin(dA * RAD) * (1 - w) + Math.sin(dB * RAD) * w, Math.cos(dA * RAD) * (1 - w) + Math.cos(dB * RAD) * w) / RAD + 360) % 360;
  return { time, wave_height_m: hs === hs ? r1(hs) : null, wave_period_s: hs === hs ? r1(s('period_s')) : null, wave_dir_deg: hs === hs && dir === dir ? Math.round(dir) % 360 : null, cycle, grid: p.resolution };
}
/** Particulates, dust and aerosol optical depth at a point from the aerosol grid. Gases (NO₂, O₃) are not part of it. */
export async function gridAir(lat, lon, { now = clock() } = {}) {
  const { a, b, w, time, cycle, p } = await bracket('air', now), s = (id) => mix(a.sample(id, lat, lon), b.sample(id, lat, lon), w);
  // gases are present only where the operator has accepted the NASA GEOS-CF source (tools/grid.mjs builds them only then)
  let no2 = null, o3 = null, gases = null; try { const g = await bracket('gas', now), q = (id) => mix(g.a.sample(id, lat, lon), g.b.sample(id, lat, lon), g.w); no2 = r1(q('no2')); o3 = r1(q('o3')); gases = g.p.registry; } catch { /* not published in this deployment */ }
  return { time, pm10: r1(s('pm10')), pm2_5: r1(s('pm2_5')), dust: r1(s('dust')), aod: r1(s('aod'), 2), no2, o3, gases, cycle, grid: p.resolution };
}
