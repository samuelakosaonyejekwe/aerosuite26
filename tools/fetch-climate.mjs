// Site design-temperature climatology, built once: `node tools/fetch-climate.mjs [--from 2015 --to 2024] [--wind-from 2022]`
// Reads the NCEP-DOE Reanalysis 2 (a US Government product, distributed by the NOAA Physical Sciences Laboratory)
// through the PSL THREDDS OPeNDAP service, computes per grid cell
//   hot99   99th percentile of the daily maximum 2 m temperature,
//   hot50   median of the daily maximum 2 m temperature,
//   cold01  1st percentile of the daily minimum 2 m temperature,
//   wind99  99th percentile of the daily maximum 10 m wind speed (the largest of the four 6-hourly values of each day),
// over the stated years, and writes them with the model terrain height on a 2.5° × 2.5° global grid to
// js/data/climate/ (index.json + grid.bin, int16). The app reads them through js/core/climate.js, corrects the
// temperatures to the site elevation with the standard lapse rate, and needs no network for them.
//
// The reanalysis grid is coarse (T62, about 1.9°): the values are regional design figures, not station extremes.

import { writeFile, mkdir } from 'node:fs/promises';
import { join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = resolve(fileURLToPath(new URL('..', import.meta.url))), outDir = join(root, 'js/data/climate');
const args = process.argv.slice(2), arg = (k, d) => (args.includes(k) ? Number(args[args.indexOf(k) + 1]) : d);
const Y0 = arg('--from', 2015), Y1 = arg('--to', 2024), W0 = arg('--wind-from', 2022);
const BASE = 'https://psl.noaa.gov/thredds/dodsC/Datasets', UA = 'AeroSuite26-data/1.0 (one-off climatology build)';
const NLAT = 94, NLON = 192;

/** One OPeNDAP request → Float32Array of the first array in the binary reply (big-endian floats after "Data:"). */
async function dods(path, query) {
  for (let k = 0; ; k++) {
    try {
      const r = await fetch(`${BASE}/${path}.dods?${query}`, { headers: { 'user-agent': UA }, signal: AbortSignal.timeout(600000) }); if (!r.ok) throw new Error(`HTTP ${r.status}`);
      const buf = Buffer.from(await r.arrayBuffer()), at = buf.indexOf('\nData:\n'); if (at < 0) throw new Error('not an OPeNDAP data reply');
      const p = at + 7, n = buf.readUInt32BE(p); if (buf.readUInt32BE(p + 4) !== n || p + 8 + 4 * n > buf.length) throw new Error('truncated reply');
      const out = new Float32Array(n); for (let i = 0; i < n; i++) out[i] = buf.readFloatBE(p + 8 + 4 * i);
      return out;
    } catch (e) { if (k >= 3) throw new Error(`${path}: ${e.message}`); await new Promise((res) => setTimeout(res, 5000 * (k + 1))); }
  }
}
/** A whole year of one variable, requested in pieces of `step` time steps (the service drops very large replies). */
async function series(path, v, nT, step = 61) {
  const cells = NLAT * NLON, out = new Float32Array(nT * cells);
  for (let t0 = 0; t0 < nT; t0 += step) { const t1 = Math.min(nT, t0 + step) - 1, a = await dods(path, `${v}[${t0}:1:${t1}][0:1:0][0:1:${NLAT - 1}][0:1:${NLON - 1}]`); if (a.length !== (t1 - t0 + 1) * cells) throw new Error(`${path}: ${a.length} values for steps ${t0}–${t1}`); out.set(a, t0 * cells); }
  return out;
}
const days = (y) => ((y % 4 === 0 && y % 100 !== 0) || y % 400 === 0 ? 366 : 365);
const ok = (v) => v > -1e30 && v < 1e30 && v === v;

/** Percentile per cell from a list of per-year arrays laid out [time][cell]. */
function percentiles(years, qs) {
  const cells = NLAT * NLON, total = years.reduce((a, y) => a + y.length / cells, 0), col = new Float32Array(total), out = qs.map(() => new Float32Array(cells));
  for (let c = 0; c < cells; c++) {
    let n = 0; for (const y of years) for (let t = 0, T = y.length / cells; t < T; t++) { const v = y[t * cells + c]; if (ok(v)) col[n++] = v; }
    const s = col.subarray(0, n).sort();
    qs.forEach((q, k) => { out[k][c] = n ? s[Math.min(n - 1, Math.floor(q * n))] : NaN; });
  }
  return out;
}

// The per-cell statistics on the reanalysis grid can be kept (--save-native <file>) and reused (--native <file>), so the
// land/sea treatment and the regridding can be changed without downloading the reanalysis again.
const sArg = (k) => (args.includes(k) ? args[args.indexOf(k) + 1] : null), nativeIn = sArg('--native'), nativeOut = sArg('--save-native');
let lat, hot99, hot50, cold01, wind99, orog, land;
if (nativeIn) {
  const { readFile } = await import('node:fs/promises'), b = await readFile(nativeIn), f = (k) => new Float32Array(b.buffer.slice(b.byteOffset + k * 4 * NLAT * NLON + 4 * NLAT, b.byteOffset + (k + 1) * 4 * NLAT * NLON + 4 * NLAT));
  lat = new Float32Array(b.buffer.slice(b.byteOffset, b.byteOffset + 4 * NLAT)); [hot99, hot50, cold01, wind99, orog, land] = [0, 1, 2, 3, 4, 5].map(f);
  console.log(`per-cell statistics read from ${nativeIn}`);
} else {
console.log(`NCEP-DOE Reanalysis 2: temperatures ${Y0}–${Y1}, wind ${W0}–${Y1}`);
lat = await dods('ncep.reanalysis2.dailyavgs/gaussian_grid/tmax.2m.gauss.' + Y1 + '.nc', 'lat'); const lon = await dods('ncep.reanalysis2.dailyavgs/gaussian_grid/tmax.2m.gauss.' + Y1 + '.nc', 'lon');
if (lat.length !== NLAT || lon.length !== NLON || !(lat[0] > lat[1])) throw new Error('unexpected reanalysis grid');
const daily = async (v) => { const ys = []; for (let y = Y0; y <= Y1; y++) { const a = await series(`ncep.reanalysis2.dailyavgs/gaussian_grid/${v}.2m.gauss.${y}.nc`, v, days(y)); ys.push(a); process.stdout.write(` ${v} ${y}`); } console.log(); return ys; };
[hot99, hot50] = percentiles(await daily('tmax'), [0.99, 0.5]);
[cold01] = percentiles(await daily('tmin'), [0.01]);
const windYears = [];
for (let y = W0; y <= Y1; y++) {
  const n = days(y) * 4, cells = NLAT * NLON;
  const u = await series(`ncep.reanalysis2/gaussian_grid/uwnd.10m.gauss.${y}.nc`, 'uwnd', n), v = await series(`ncep.reanalysis2/gaussian_grid/vwnd.10m.gauss.${y}.nc`, 'vwnd', n);
  const mx = new Float32Array(days(y) * cells);
  for (let t = 0; t < n; t++) { const d = (t >> 2) * cells, o = t * cells; for (let c = 0; c < cells; c++) { const a = u[o + c], b = v[o + c]; if (ok(a) && ok(b)) { const s = Math.hypot(a, b); if (s > mx[d + c]) mx[d + c] = s; } } }
  windYears.push(mx); process.stdout.write(` wind ${y}`);
}
console.log();
[wind99] = percentiles(windYears, [0.99]);
orog = await dods('ncep.reanalysis2/gaussian_grid/hgt.sfc.gauss.nc', `hgt[0:1:0][0:1:${NLAT - 1}][0:1:${NLON - 1}]`);
if (orog.length !== NLAT * NLON) throw new Error('terrain height: unexpected size');

land = await dods('ncep.reanalysis/surface_gauss/land.sfc.gauss.nc', `land[0:1:0][0:1:${NLAT - 1}][0:1:${NLON - 1}]`);
if (land.length !== NLAT * NLON) throw new Error('land mask: unexpected size');
if (nativeOut) { const b = Buffer.alloc(4 * NLAT + 6 * 4 * NLAT * NLON); Buffer.from(lat.buffer).copy(b, 0); [hot99, hot50, cold01, wind99, orog, land].forEach((f, k) => Buffer.from(f.buffer, f.byteOffset, f.byteLength).copy(b, 4 * NLAT + k * 4 * NLAT * NLON)); await writeFile(nativeOut, b); }
}

// ---- land values for land sites -------------------------------------------------------------------
// A reanalysis cell on a coast is partly sea: its wind is far stronger than over the land an aerodrome stands on.
// For the land wind field every sea cell takes the distance-weighted mean of the
// land cells within three cells of it (left unchanged where there is none: open ocean), so that interpolation to a
// coastal site mixes land values only. The unmodified fields are kept for sites at sea.
function landFilled(f) {
  const out = Float32Array.from(f);
  for (let j = 0; j < NLAT; j++) for (let i = 0; i < NLON; i++) {
    if (land[j * NLON + i] >= 0.5) continue;
    for (let r = 1; r <= 3; r++) {
      let s = 0, w = 0;
      for (let dj = -r; dj <= r; dj++) for (let di = -r; di <= r; di++) { const jj = j + dj; if (jj < 0 || jj >= NLAT || Math.max(Math.abs(dj), Math.abs(di)) !== r) continue; const c = jj * NLON + ((i + di + NLON) % NLON); if (land[c] >= 0.5) { const k = 1 / (dj * dj + di * di); s += f[c] * k; w += k; } }
      if (w > 0) { out[j * NLON + i] = s / w; break; }
    }
  }
  return out;
}

// ---- Gaussian grid → 2.5° regular grid (bilinear; latitude rows are not evenly spaced) ------------
const STEP = 2.5, NJ = 73, NI = 144;
function regrid(f) {
  const out = new Float32Array(NJ * NI);
  for (let j = 0; j < NJ; j++) {
    const la = 90 - j * STEP; let j0 = 0; while (j0 < NLAT - 2 && lat[j0 + 1] > la) j0++;
    const fy = Math.min(1, Math.max(0, (lat[j0] - la) / (lat[j0] - lat[j0 + 1])));
    for (let i = 0; i < NI; i++) { const x = (i * STEP) / (360 / NLON), i0 = Math.floor(x) % NLON, i1 = (i0 + 1) % NLON, fx = x - Math.floor(x); out[j * NI + i] = (f[j0 * NLON + i0] * (1 - fx) + f[j0 * NLON + i1] * fx) * (1 - fy) + (f[(j0 + 1) * NLON + i0] * (1 - fx) + f[(j0 + 1) * NLON + i1] * fx) * fy; }
  }
  return out;
}
const L = landFilled, landPct = Float32Array.from(land, (v) => 100 * v);
// Temperatures stay as the reanalysis gives them (checked against ERA5 site values at coastal airports, the unmodified
// cells are closer than land-only values, which run hot). The strong-wind value is the one a partly-sea cell distorts
// most: for land sites it comes from the land-filled field; the unmodified field serves sites at sea.
const FIELDS = [['hot99_C', hot99, 0.01, -273.15], ['hot50_C', hot50, 0.01, -273.15], ['cold01_C', cold01, 0.01, -273.15], ['wind99_ms', L(wind99), 0.01, 0], ['orog_m', orog, 1, 0], ['wind99_sea_ms', wind99, 0.01, 0], ['land_pct', landPct, 1, 0]];
const bin = Buffer.alloc(FIELDS.length * NJ * NI * 2); let o = 0; const stats = {};
for (const [id, f, scale, add] of FIELDS) { const g = regrid(f); let mn = Infinity, mx = -Infinity; for (const v of g) { const x = v + add; if (!(x === x)) throw new Error(`${id}: missing value in the grid`); mn = Math.min(mn, x); mx = Math.max(mx, x); bin.writeInt16LE(Math.max(-32768, Math.min(32767, Math.round(x / scale))), o); o += 2; } stats[id] = [Math.round(mn * 100) / 100, Math.round(mx * 100) / 100]; }
if (!(stats.hot99_C[1] > 35 && stats.hot99_C[1] < 60 && stats.cold01_C[0] < -40 && stats.cold01_C[0] > -95 && stats.wind99_sea_ms[1] > 10 && stats.wind99_sea_ms[1] < 60)) throw new Error('implausible statistics: ' + JSON.stringify(stats));
await mkdir(outDir, { recursive: true });
await writeFile(join(outDir, 'grid.bin'), bin);
await writeFile(join(outDir, 'index.json'), JSON.stringify({
  schema: 1, source: 'NCEP-DOE Reanalysis 2', publisher: 'NOAA Physical Sciences Laboratory (data: NOAA/NCEP and US Department of Energy)', home: 'https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html', service: BASE, license: 'US Government product, public domain; acknowledgement requested',
  attribution: 'NCEP-DOE Reanalysis 2 data provided by the NOAA PSL, Boulder, Colorado, USA, from their website at https://psl.noaa.gov', generated: new Date().toISOString(), registry: 'noaa-reanalysis2',
  period: { temperature: [Y0, Y1], wind: [W0, Y1] }, method: 'Per reanalysis cell (T62 Gaussian grid, about 1.9°): 99th percentile and median of the daily maximum 2 m temperature, 1st percentile of the daily minimum 2 m temperature, 99th percentile of the daily maximum of the four 6-hourly 10 m wind speeds; then bilinear interpolation to 2.5°. Strong wind is stored twice: wind99_ms for land sites, in which sea cells of the reanalysis first take the distance-weighted mean of the land cells within three cells (a coastal cell is partly sea and far windier than the land beside it), and wind99_sea_ms, unmodified, for sites at sea; land_pct is the reanalysis land share. The reanalysis 10 m wind runs high over mid-latitude land compared with finer reanalyses: treat the strong-wind value as a conservative regional figure. Regional values: the app corrects temperature to the site elevation with 6.5 K/km from the model terrain height.',
  grid: { step_deg: STEP, nj: NJ, ni: NI, lat0: 90, lon0: 0 }, type: 'int16le', fields: FIELDS.map(([id, , scale]) => ({ id, scale })), range: stats, file: 'grid.bin', bytes: bin.length,
}, null, 1) + '\n');
console.log(`js/data/climate: ${bin.length} bytes`, JSON.stringify(stats));
