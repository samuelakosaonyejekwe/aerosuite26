// Global forecast grids for the app, built by the scheduled snapshot job: `node tools/grid.mjs [--out data/grid]`
// (also called from tools/snapshot.mjs). Everything comes from open government model output that may be reused
// commercially, is fetched server-side without keys, and is written as small static files the app interpolates
// to any site — so weather aloft, freezing level, gusts, visibility, sea state and dust are available to every
// deployment, online or from the offline copy, without a per-user weather API:
//   wx   NOAA/NCEP GFS 1° (US Government, public domain) through the NOMADS grib filter: surface fields and wind,
//        temperature and geopotential height at eleven pressure levels — on a 1° / 2.5° grid for the three valid
//        times nearest to now, on a 2.5° / 5° grid for the two later ones;
//   sea  DWD global wave model GWAM 0.25° (Deutscher Wetterdienst open data, CC BY 4.0): significant wave height,
//        mean period and mean direction, averaged over the sea points of each 2° cell;
//   air  NOAA GEFS-Aerosols 0.25° (public domain) from the NOAA Open Data bucket: PM2.5, PM10, dust and aerosol
//        optical depth on a 2.5° grid.
// Each product holds a few valid times covering about the next 24 hours, so a job that runs every three hours
// (or misses a few runs) always leaves a valid pair of slices to interpolate between. Values are quantised to
// one or two bytes (see FIELDS), difference-coded and gzipped (encodeSlice); data/grid/index.json describes the
// layout and is all the app needs to read them.

import { writeFile, mkdir, readFile, readdir, rm } from 'node:fs/promises';
import { join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { gribMessages, gribValues, gridReader } from './grib2.mjs';
import { bunzip2 } from './bunzip2.mjs';
import { gzipSync } from 'node:zlib';
import net from 'node:net';

// Some servers answer the first connection attempt slowly; Node's default of 250 ms per address is too short for them.
net.setDefaultAutoSelectFamilyAttemptTimeout?.(5000);

const root = resolve(fileURLToPath(new URL('..', import.meta.url)));
const UA = 'AeroSuite26-snapshot/1.0 (scheduled public-data snapshot; static web app)', H = 3600e3;
export const LEVELS = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150];
const grid = (step) => ({ step, nj: Math.round(180 / step) + 1, ni: Math.round(360 / step) }); // rows from 90°N, columns from 0°E
/** Weather: the slices nearest in time on the fine grids (1° surface, 2.5° aloft), the later ones on the coarse grids. Sea state and aerosols use `w` and `s` of the coarse set. */
export const FINE = { s: grid(1), a: grid(2.5) }, GRIDS = { s: grid(2.5), a: grid(5), w: grid(2) }, FINE_SLICES = 3;
/** Field layout of each slice file, in file order. value = raw × scale + offset (squared first when `sq`); `miss` marks "no value". */
export const FIELDS = {
  wx: [
    { id: 't2m_C', grid: 's', type: 'i16', scale: 0.01 }, { id: 'rh_pct', grid: 's', type: 'u8', scale: 1 }, { id: 'qnh_hPa', grid: 's', type: 'i16', scale: 0.1, offset: 1000 },
    { id: 'u10_ms', grid: 's', type: 'i8', scale: 0.5 }, { id: 'v10_ms', grid: 's', type: 'i8', scale: 0.5 }, { id: 'gust_ms', grid: 's', type: 'u8', scale: 0.5 }, { id: 'vis_m', grid: 's', type: 'u8', scale: 100 },
    { id: 'cloud_pct', grid: 's', type: 'u8', scale: 1 }, { id: 'precip_mm_h', grid: 's', type: 'u8', scale: 0.05, sq: true }, { id: 'fz_m', grid: 's', type: 'i16', scale: 1 },
    ...LEVELS.flatMap((l) => [{ id: `u_${l}`, grid: 'a', type: 'i8', scale: 1 }, { id: `v_${l}`, grid: 'a', type: 'i8', scale: 1 }, { id: `T_${l}`, grid: 'a', type: 'u8', scale: 0.5, offset: -93.15 }, { id: `z_${l}`, grid: 'a', type: 'i16', scale: 1 }]),
  ],
  sea: [{ id: 'hs_m', grid: 'w', type: 'u8', scale: 0.1, miss: 255 }, { id: 'period_s', grid: 'w', type: 'u8', scale: 0.1, miss: 255 }, { id: 'dir_deg', grid: 'w', type: 'u8', scale: 2, miss: 255 }],
  gas: [{ id: 'no2', grid: 's', type: 'i16', scale: 0.1 }, { id: 'o3', grid: 's', type: 'i16', scale: 0.1 }],
  air: [{ id: 'pm2_5', grid: 's', type: 'i16', scale: 0.1 }, { id: 'pm10', grid: 's', type: 'i16', scale: 0.1 }, { id: 'dust', grid: 's', type: 'i16', scale: 0.1 }, { id: 'aod', grid: 's', type: 'u8', scale: 0.02 }],
};
const BYTES = { i8: 1, u8: 1, i16: 2 }, RANGE = { i8: [-128, 127], u8: [0, 255], i16: [-32768, 32767] };
export const sliceBytes = (product, grids = GRIDS) => FIELDS[product].reduce((n, f) => n + BYTES[f.type] * grids[f.grid].nj * grids[f.grid].ni, 0);

/** Pack { fieldId: Float32Array on the field's grid } into one slice buffer. */
export function packSlice(product, data, grids = GRIDS) {
  const buf = Buffer.alloc(sliceBytes(product, grids)); let o = 0;
  for (const f of FIELDS[product]) {
    const g = grids[f.grid], src = data[f.id], [lo, hi] = RANGE[f.type]; if (!src || src.length !== g.nj * g.ni) throw new Error(`${product}: field ${f.id} is missing`);
    for (let k = 0; k < src.length; k++) {
      let v = src[k], raw;
      if (v !== v) { if (f.miss == null) throw new Error(`${product}: ${f.id} has a missing value`); raw = f.miss; }
      else { v -= f.offset || 0; if (f.sq) v = Math.sqrt(Math.max(0, v)); raw = Math.max(lo, Math.min(f.miss != null ? f.miss - 1 : hi, Math.round(v / f.scale))); }
      if (f.type === 'i16') { buf.writeInt16LE(raw, o); o += 2; } else if (f.type === 'i8') buf.writeInt8(raw, o++); else buf.writeUInt8(raw, o++);
    }
  }
  return buf;
}

/**
 * File encoding "delta-gzip": every field is replaced by the difference from its left neighbour (from the value above at the
 * start of a row), 16-bit differences are split into a low-byte plane and a high-byte plane, and the result is gzipped.
 * Smooth fields shrink to less than half; js/core/gridwx.js reverses it.
 */
export function encodeSlice(raw, fields, grids = GRIDS) {
  const out = Buffer.alloc(raw.length); let o = 0;
  for (const f of fields) {
    const g = grids[f.grid], n = g.nj * g.ni, wide = f.type === 'i16', at = (k) => (wide ? raw.readInt16LE(o + 2 * k) : raw[o + k]);
    for (let j = 0; j < g.nj; j++) for (let i = 0; i < g.ni; i++) { const k = j * g.ni + i, d = at(k) - (i ? at(k - 1) : j ? at(k - g.ni) : 0); if (wide) { out[o + k] = d & 255; out[o + n + k] = (d >> 8) & 255; } else out[o + k] = d & 255; }
    o += (wide ? 2 : 1) * n;
  }
  return gzipSync(out, { level: 9 });
}
const ENC = 'delta-gzip';

function makeGet() {
  return async (url, opt = {}) => {
    let last;
    for (let k = 0; k < 3; k++) {
      try { const r = await fetch(url, { ...opt, headers: { 'user-agent': UA, ...(opt.headers || {}) }, signal: AbortSignal.timeout(120000) }); if (r.status >= 500 && k < 2) { last = new Error(`HTTP ${r.status}`); await new Promise((res) => setTimeout(res, 3000)); continue; } return r; }
      catch (e) { last = e; await new Promise((res) => setTimeout(res, 3000)); }
    }
    throw new Error(last?.cause?.code || last?.message || 'request failed');
  };
}
const bytesOf = async (r, what) => { if (!r.ok) throw new Error(`${what}: HTTP ${r.status}`); return Buffer.from(await r.arrayBuffer()); };
const ymd = (t) => new Date(t).toISOString().slice(0, 10).replace(/-/g, ''), hh = (t) => String(new Date(t).getUTCHours()).padStart(2, '0'), f3 = (f) => String(f).padStart(3, '0');
/** Sample a decoded 1°/0.25° field at every point of a target grid (bilinear). */
function toGrid(reader, g, fn = (v) => v) { const out = new Float32Array(g.nj * g.ni); for (let j = 0; j < g.nj; j++) for (let i = 0; i < g.ni; i++) out[j * g.ni + i] = fn(reader.sample(90 - j * g.step, i * g.step)); return out; }
/** Valid times for a model run that started at `cycle`: the newest 3-hourly step not after `now` (at least +3 h, so forecast-only fields exist), then every 6 h. */
const validSteps = (cycle, now, n) => { const f0 = Math.max(3, Math.floor((now - cycle) / (3 * H)) * 3); return Array.from({ length: n }, (_, k) => f0 + 6 * k); };

// ---- wx: NOAA GFS -------------------------------------------------------------------------------
const GFS_VARS = ['TMP', 'RH', 'UGRD', 'VGRD', 'HGT', 'PRMSL', 'GUST', 'VIS', 'TCDC', 'PRATE'], GFS_LEVS = [...LEVELS.map((l) => `${l}_mb`), '2_m_above_ground', '10_m_above_ground', 'surface', 'mean_sea_level', '0C_isotherm', 'entire_atmosphere'];
const gfsUrl = (cycle, f) => `https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_1p00.pl?dir=%2Fgfs.${ymd(cycle)}%2F${hh(cycle)}%2Fatmos&file=gfs.t${hh(cycle)}z.pgrb2.1p00.f${f3(f)}&${GFS_VARS.map((v) => `var_${v}=on`).join('&')}&${GFS_LEVS.map((l) => `lev_${l}=on`).join('&')}`;
/** One GFS file → { fields on the app grids, orog }. Refuses a file that lacks any field rather than writing a partial slice. */
export function gfsSlice(buf, grids = GRIDS) {
  const ms = gribMessages(buf), pick = (cat, num, lt, lv, tpl) => { const c = ms.filter((m) => m.discipline === 0 && m.category === cat && m.number === num && m.levelType === lt && (lv == null || Math.abs(m.level - lv) < 0.5)); const m = (tpl != null && c.find((x) => x.template === tpl)) || c[0]; if (!m) throw new Error(`GFS file lacks field ${cat}/${num} at level type ${lt}${lv != null ? ` ${lv}` : ''}`); return gridReader(m, gribValues(m)); };
  const S = grids.s, A = grids.a, out = {};
  out.t2m_C = toGrid(pick(0, 0, 103, 2), S, (v) => v - 273.15); out.rh_pct = toGrid(pick(1, 1, 103, 2), S); out.qnh_hPa = toGrid(pick(3, 1, 101), S, (v) => v / 100);
  out.u10_ms = toGrid(pick(2, 2, 103, 10), S); out.v10_ms = toGrid(pick(2, 3, 103, 10), S); out.gust_ms = toGrid(pick(2, 22, 1), S); out.vis_m = toGrid(pick(19, 0, 1), S);
  out.cloud_pct = toGrid(pick(6, 1, 10), S); out.precip_mm_h = toGrid(pick(1, 7, 1, null, 0), S, (v) => v * 3600); out.fz_m = toGrid(pick(3, 5, 4), S);
  for (const l of LEVELS) { out[`u_${l}`] = toGrid(pick(2, 2, 100, l * 100), A); out[`v_${l}`] = toGrid(pick(2, 3, 100, l * 100), A); out[`T_${l}`] = toGrid(pick(0, 0, 100, l * 100), A, (v) => v - 273.15); out[`z_${l}`] = toGrid(pick(3, 5, 100, l * 100), A); }
  const orog = toGrid(pick(3, 5, 1), S), t = out.t2m_C, z5 = out.z_500;
  let mn = Infinity, mx = -Infinity; for (const v of t) { mn = Math.min(mn, v); mx = Math.max(mx, v); }
  if (!(mn > -95 && mn < -10 && mx > 25 && mx < 62) || !(z5[0] > 4600 && z5[0] < 6100)) throw new Error(`implausible GFS field values (2 m temperature ${mn.toFixed(1)}…${mx.toFixed(1)} °C)`);
  return { fields: out, orog };
}
async function buildWx(get, now, outDir, nSlices = 5) {
  let lastErr = null;
  for (let back = 0; back < 5; back++) { // newest run first; a run is usable once every file we need has been published
    const cycle = Math.floor((now - 3 * H) / (6 * H)) * 6 * H - back * 6 * H, steps = validSteps(cycle, now, nSlices), slices = [];
    try {
      let orog = null;
      for (const [k, f] of steps.entries()) {
        const buf = await bytesOf(await get(gfsUrl(cycle, f)), `GFS ${ymd(cycle)} ${hh(cycle)}Z f${f3(f)}`); if (buf.toString('latin1', 0, 4) !== 'GRIB') throw new Error(`GFS ${ymd(cycle)} ${hh(cycle)}Z f${f3(f)} is not published yet`);
        const grids = k < FINE_SLICES ? FINE : GRIDS, sl = gfsSlice(buf, grids); orog ||= sl.orog; // the first slice is on the fine grid, and so is the terrain height
        slices.push({ f, fine: k < FINE_SLICES, grids: { s: grids.s, a: grids.a }, buf: encodeSlice(packSlice('wx', sl.fields, grids), FIELDS.wx, grids) });
      }
      const times = []; for (const [k, sl] of slices.entries()) { await writeFile(join(outDir, `wx_${k}.bin.gz`), sl.buf); times.push({ valid: new Date(cycle + sl.f * H).toISOString(), forecast_h: sl.f, file: `wx_${k}.bin.gz`, enc: ENC, grids: sl.grids, resolution: sl.fine ? 'surface 1°, pressure levels 2.5°' : 'surface 2.5°, pressure levels 5°' }); }
      const of = [{ id: 'orog', grid: 's', type: 'i16', scale: 1 }], oraw = Buffer.alloc(orog.length * 2); orog.forEach((v, i) => oraw.writeInt16LE(Math.round(v), 2 * i)); const ob = encodeSlice(oraw, of, FINE); await writeFile(join(outDir, 'orog.bin.gz'), ob);
      return { ok: true, source: 'NOAA/NCEP Global Forecast System (GFS), 1° output', publisher: 'US National Weather Service, National Centers for Environmental Prediction', home: 'https://nomads.ncep.noaa.gov/', registry: 'noaa-gfs', attribution: 'Forecast grid: NOAA/NWS/NCEP Global Forecast System (public domain)',
        cycle: new Date(cycle).toISOString(), resolution: `surface fields 1°, pressure levels 2.5° for the first ${FINE_SLICES} valid times; 2.5° and 5° for the later ones`, levels_hPa: LEVELS, grids: GRIDS, fields: FIELDS.wx, times, static: { orog: { file: 'orog.bin.gz', enc: ENC, grids: { s: FINE.s }, grid: 's', type: 'i16', scale: 1 } }, bytes: slices.reduce((n, x) => n + x.buf.length, 0) + ob.length };
    } catch (e) { lastErr = e; }
  }
  throw lastErr;
}

// ---- sea: DWD GWAM ------------------------------------------------------------------------------
const gwamUrl = (cycle, v, f) => `https://opendata.dwd.de/weather/maritime/wave_models/gwam/grib/${hh(cycle)}/${v}/GWAM_${v.toUpperCase()}_${ymd(cycle)}${hh(cycle)}_${f3(f)}.grib2.bz2`;
/** Mean over the sea points inside each 2° cell (NaN where the cell holds no sea); directions are averaged as vectors. */
function seaCells(reader, g, mode) {
  const out = new Float32Array(g.nj * g.ni), R = Math.PI / 180;
  for (let j = 0; j < g.nj; j++) for (let i = 0; i < g.ni; i++) {
    let s = 0, c = 0, n = 0;
    for (let dy = -0.75; dy <= 0.76; dy += 0.5) for (let dx = -0.75; dx <= 0.76; dx += 0.5) { const la = 90 - j * g.step + dy; if (la > 90 || la < -90) continue; const v = reader.sample(la, i * g.step + dx); if (v !== v) continue; if (mode === 'dir') { s += Math.sin(v * R); c += Math.cos(v * R); } else s += v; n++; }
    out[j * g.ni + i] = !n ? NaN : mode === 'dir' ? ((Math.atan2(s, c) / R) + 360) % 360 : s / n;
  }
  return out;
}
async function buildSea(get, now, outDir, wantValid) {
  let lastErr = null;
  for (let back = 0; back < 4; back++) {
    const cycle = Math.floor((now - 4 * H) / (12 * H)) * 12 * H - back * 12 * H, steps = wantValid.map((t) => (t - cycle) / H).filter((f) => f >= 0 && f <= 174 && f % 3 === 0);
    if (steps.length < 2) continue;
    try {
      const slices = [];
      for (const f of steps) {
        const one = async (v) => { const m = gribMessages(bunzip2(await bytesOf(await get(gwamUrl(cycle, v, f)), `GWAM ${v} ${ymd(cycle)}${hh(cycle)} +${f}`)))[0]; if (!m) throw new Error(`GWAM ${v}: empty file`); return gridReader(m, gribValues(m)); };
        const hs = seaCells(await one('swh'), GRIDS.w), per = seaCells(await one('tm10'), GRIDS.w), dir = seaCells(await one('mwd'), GRIDS.w, 'dir');
        let n = 0, mx = 0; for (const v of hs) if (v === v) { n++; mx = Math.max(mx, v); } if (n < 5000 || !(mx > 2 && mx < 25)) throw new Error(`implausible wave field (${n} sea cells, highest ${mx.toFixed(1)} m)`);
        slices.push({ f, buf: encodeSlice(packSlice('sea', { hs_m: hs, period_s: per, dir_deg: dir }), FIELDS.sea) });
      }
      const times = []; for (const [k, s] of slices.entries()) { await writeFile(join(outDir, `sea_${k}.bin.gz`), s.buf); times.push({ valid: new Date(cycle + s.f * H).toISOString(), forecast_h: s.f, file: `sea_${k}.bin.gz`, enc: ENC }); }
      return { ok: true, source: 'DWD global wave model GWAM, 0.25° output', publisher: 'Deutscher Wetterdienst (DWD), open data server', home: 'https://opendata.dwd.de/weather/maritime/wave_models/gwam/', registry: 'dwd-gwam', attribution: 'Sea state: Deutscher Wetterdienst (DWD), global wave model GWAM, CC BY 4.0 — averaged to 2° cells',
        cycle: new Date(cycle).toISOString(), resolution: '2° cells: mean over the sea points of the 0.25° model grid', grids: GRIDS, fields: FIELDS.sea, times, bytes: slices.reduce((n, s) => n + s.buf.length, 0) };
    } catch (e) { lastErr = e; }
  }
  throw lastErr || new Error('no wave model run covers the wanted times');
}

// ---- air: NOAA GEFS-Aerosols --------------------------------------------------------------------
const chemUrl = (cycle, f) => `https://noaa-gefs-pds.s3.amazonaws.com/gefs.${ymd(cycle)}/${hh(cycle)}/chem/pgrb2ap25/gefs.chem.t${hh(cycle)}z.a2d_0p25.f${f3(f)}.grib2`;
const AIR_PICK = { pm2_5: (l) => /:PMTF:surface:/.test(l) && /aerosol=Total aerosol/.test(l), pm10: (l) => /:PMTC:surface:/.test(l) && /aerosol=Total aerosol/.test(l), dust: (l) => /:PMTC:surface:/.test(l) && /aerosol=Dust dry/.test(l), aod: (l) => /:AOTK:entire atmosphere:/.test(l) && /aerosol=Total aerosol/.test(l) && /wavelength >=5\.45e-07,<=5\.65e-07/.test(l) };
async function buildAir(get, now, outDir, wantValid) {
  let lastErr = null;
  for (let back = 0; back < 8; back++) {
    const cycle = Math.floor((now - 6 * H) / (6 * H)) * 6 * H - back * 6 * H, steps = wantValid.map((t) => (t - cycle) / H).filter((f) => f >= 0 && f <= 120 && f % 3 === 0);
    if (steps.length < 2) continue;
    try {
      const slices = [];
      for (const f of steps) {
        const url = chemUrl(cycle, f), ir = await get(`${url}.idx`); if (!ir.ok) throw new Error(`GEFS-Aerosols ${ymd(cycle)} ${hh(cycle)}Z f${f3(f)}: HTTP ${ir.status}`);
        const idx = (await ir.text()).trim().split('\n').map((l) => ({ l, at: Number(l.split(':')[1]) })), data = {};
        for (const [id, test] of Object.entries(AIR_PICK)) {
          const k = idx.findIndex((x) => test(x.l)); if (k < 0 || !Number.isFinite(idx[k].at)) throw new Error(`GEFS-Aerosols index lacks ${id}`);
          const range = `bytes=${idx[k].at}-${k + 1 < idx.length ? idx[k + 1].at - 1 : ''}`, m = gribMessages(await bytesOf(await get(url, { headers: { range } }), `GEFS-Aerosols ${id}`))[0]; if (!m) throw new Error(`GEFS-Aerosols ${id}: no message in the byte range`);
          data[id] = toGrid(gridReader(m, gribValues(m)), GRIDS.s, (v) => Math.max(0, v));
        }
        let mean = 0; for (const v of data.pm2_5) mean += v / data.pm2_5.length; if (!(mean > 0.5 && mean < 100) || data.aod.some((v) => v > 12)) throw new Error(`implausible aerosol field (mean PM2.5 ${mean.toFixed(2)} µg/m³)`);
        slices.push({ f, buf: encodeSlice(packSlice('air', data), FIELDS.air) });
      }
      const times = []; for (const [k, s] of slices.entries()) { await writeFile(join(outDir, `air_${k}.bin.gz`), s.buf); times.push({ valid: new Date(cycle + s.f * H).toISOString(), forecast_h: s.f, file: `air_${k}.bin.gz`, enc: ENC }); }
      return { ok: true, source: 'NOAA GEFS-Aerosols (GEFS chemistry member), 0.25° output', publisher: 'US National Weather Service, National Centers for Environmental Prediction; NOAA Open Data Dissemination', home: 'https://registry.opendata.aws/noaa-gefs/', registry: 'noaa-gefs-aerosols', attribution: 'Aerosol forecast: NOAA/NWS/NCEP GEFS-Aerosols (public domain)',
        cycle: new Date(cycle).toISOString(), resolution: '2.5° (sampled from the 0.25° output)', units: { pm2_5: 'µg/m³', pm10: 'µg/m³', dust: 'µg/m³ (dust up to 10 µm)', aod: 'aerosol optical depth at 550 nm' }, grids: GRIDS, fields: FIELDS.air, times, bytes: slices.reduce((n, s) => n + s.buf.length, 0) };
    } catch (e) { lastErr = e; }
  }
  throw lastErr || new Error('no aerosol run covers the wanted times');
}

// ---- gas: NASA GEOS-CF (only where the operator has accepted the source) --------------------------
// Surface nitrogen dioxide and ozone from NASA GMAO's GEOS Composition Forecast through its OPeNDAP server (plain text,
// every tenth point of the 0.25° grid = 2.5°). NASA publishes no explicit reuse grant for this product and describes
// its forecasts as experimental, so the licence registry classes it "unclear": it is fetched only when the deployment
// lists "nasa-geos-cf" under "accept" in config.json.
const GEOSCF = 'https://opendap.nccs.nasa.gov/dods/gmao/geos-cf/v2/fcst/aqc_tavg_1hr_glo_L1440x721_slv.latest';
const UG = { no2: 46.01 / 24.45, o3: 48.0 / 24.45 }; // µg/m³ per ppb at 25 °C and 1013.25 hPa
/** GrADS time stamp of the first step, e.g. "09:30z08oct2026" → ms. */
export const gradsTime = (t) => { const m = /^(\d\d):(\d\d)z(\d\d)([a-z]{3})(\d{4})$/i.exec(t || ''), mo = m ? ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'].indexOf(m[4].toLowerCase()) : -1; return mo < 0 ? NaN : Date.UTC(Number(m[5]), mo, Number(m[3]), Number(m[1]), Number(m[2])); };
/** One variable block of an OPeNDAP text reply ("name, [1][1][73][144]" then one line per latitude, south to north, from 180°W) → Float32Array on grid s. */
export function geosCfField(text, name) {
  const g = GRIDS.s, at = text.indexOf(`${name}.${name}, [`) >= 0 ? text.indexOf(`${name}.${name}, [`) : text.indexOf(`${name}, [`); if (at < 0) throw new Error(`GEOS-CF reply lacks ${name}`);
  const lines = text.slice(at).split('\n').slice(1, 1 + g.nj), out = new Float32Array(g.nj * g.ni); if (lines.length !== g.nj) throw new Error(`GEOS-CF ${name}: truncated reply`);
  lines.forEach((l, r) => { const v = l.split(', ').slice(1).map(Number); if (v.length !== g.ni || v.some((x) => !(x >= 0 && x < 1e-3))) throw new Error(`GEOS-CF ${name}: unexpected row`); const j = g.nj - 1 - r; for (let i = 0; i < g.ni; i++) out[j * g.ni + ((i + g.ni / 2) % g.ni)] = v[i] * 1e9 * UG[name]; });
  return out;
}
async function buildGas(get, now, outDir, wantValid) {
  const das = await (await get(`${GEOSCF}.das`)).text(), t0 = gradsTime(/time \{[\s\S]*?minimum "([^"]+)"/.exec(das)?.[1]), nT = 120; if (!Number.isFinite(t0)) throw new Error('GEOS-CF: start time not found');
  const steps = wantValid.map((t) => Math.round((t - t0) / H)).filter((k) => k >= 0 && k < nT); if (steps.length < 2) throw new Error('the current GEOS-CF forecast does not cover the wanted times');
  const slices = [];
  for (const k of steps) {
    const q = (v) => `${v}%5B${k}:1:${k}%5D%5B0:1:0%5D%5B0:10:720%5D%5B0:10:1439%5D`, r = await get(`${GEOSCF}.ascii?${q('no2')},${q('o3')}`); if (!r.ok) throw new Error(`GEOS-CF: HTTP ${r.status}`);
    const text = await r.text(), data = { no2: geosCfField(text, 'no2'), o3: geosCfField(text, 'o3') }; let mean = 0; for (const v of data.o3) mean += v / data.o3.length; if (!(mean > 10 && mean < 150)) throw new Error(`implausible ozone field (mean ${mean.toFixed(1)} µg/m³)`);
    slices.push({ t: t0 + k * H, buf: encodeSlice(packSlice('gas', data), FIELDS.gas) });
  }
  const times = []; for (const [k, sl] of slices.entries()) { await writeFile(join(outDir, `gas_${k}.bin.gz`), sl.buf); times.push({ valid: new Date(sl.t).toISOString(), file: `gas_${k}.bin.gz`, enc: ENC }); }
  return { ok: true, source: 'NASA GMAO GEOS Composition Forecast (GEOS-CF v2), 0.25° output', publisher: 'NASA Global Modeling and Assimilation Office', home: 'https://gmao.gsfc.nasa.gov/gmao-products/geos-cf/', registry: 'nasa-geos-cf', attribution: 'NO₂ and ozone: NASA GMAO GEOS Composition Forecast (GEOS-CF)',
    cycle: new Date(t0).toISOString(), resolution: '2.5° (every tenth point of the 0.25° output)', units: { no2: 'µg/m³ (from mole fraction at 25 °C, 1013.25 hPa)', o3: 'µg/m³ (from mole fraction at 25 °C, 1013.25 hPa)' }, grids: GRIDS, fields: FIELDS.gas, times, bytes: slices.reduce((n, x) => n + x.buf.length, 0) };
}

/** Build every product into `outDir`. A product that fails keeps its previous entry (and files) and is marked ok: false. */
export async function buildGrids({ outDir = join(root, 'data/grid'), now = Date.now(), only = null, accept = [] } = {}) {
  await mkdir(outDir, { recursive: true });
  let prev = null; try { prev = JSON.parse(await readFile(join(outDir, 'index.json'), 'utf8')); } catch { /* first run */ }
  const get = makeGet(), products = {}, want = (id) => !only || only.includes(id);
  const run = async (id, fn) => { if (!want(id)) { if (prev?.products?.[id]) products[id] = prev.products[id]; return; } try { products[id] = await fn(); } catch (e) { products[id] = { ...(prev?.products?.[id] || {}), ok: false, error: String(e.message || e).slice(0, 200), tried: new Date(now).toISOString() }; } };
  await run('wx', () => buildWx(get, now, outDir));
  const base = (products.wx?.times || []).map((t) => Date.parse(t.valid)), valid = base.length ? base : validSteps(Math.floor(now / (6 * H)) * 6 * H - 6 * H, now, 5).map((f) => Math.floor(now / (6 * H)) * 6 * H - 6 * H + f * H);
  await run('sea', () => buildSea(get, now, outDir, valid));
  await run('air', () => buildAir(get, now, outDir, valid.filter((_, k) => k % 2 === 0)));
  if (accept.includes('nasa-geos-cf')) await run('gas', () => buildGas(get, now, outDir, valid.filter((_, k) => k % 2 === 0))); // never without the operator's acceptance
  const index = { schema: 2, generated: new Date(now).toISOString(), generator: 'tools/grid.mjs', about: 'Global forecast grids interpolated by the app (js/core/gridwx.js). Rows run from 90°N southwards, columns from 0°E eastwards; int16 values are little-endian. Files marked enc "delta-gzip" hold, gzipped, each field as differences from the left neighbour (from the value above at the start of a row), 16-bit differences as a low-byte plane followed by a high-byte plane. A valid time may carry its own grids.', products };
  const used = new Set(['index.json', ...Object.values(products).flatMap((p) => [...(p.times || []).map((t) => t.file), ...Object.values(p.static || {}).map((x) => x.file)])]);
  for (const f of await readdir(outDir)) if (!used.has(f) && /^(wx|sea|air|gas)_\d+\.bin(\.gz)?$|^orog\.bin(\.gz)?$/.test(f)) await rm(join(outDir, f)); // files of an earlier layout or run
  await writeFile(join(outDir, 'index.json'), JSON.stringify(index, null, 1) + '\n');
  return index;
}

const isMain = process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url);
if (isMain) {
  const args = process.argv.slice(2), arg = (k) => (args.includes(k) ? args[args.indexOf(k) + 1] : null);
  let accept = arg('--accept')?.split(',') || []; if (!arg('--accept')) try { accept = JSON.parse(await readFile(join(root, 'config.json'), 'utf8')).accept || []; } catch { /* defaults */ }
  const idx = await buildGrids({ outDir: arg('--out') ? resolve(arg('--out')) : undefined, only: arg('--only')?.split(',') || null, accept });
  for (const [id, p] of Object.entries(idx.products)) console.log(`${p.ok ? 'ok  ' : 'FAIL'} ${id.padEnd(4)} ${p.ok ? `${p.cycle} · ${p.times.map((t) => (t.forecast_h != null ? `+${t.forecast_h}h` : t.valid.slice(5, 13))).join(' ')} · ${(p.bytes / 1e3).toFixed(0)} kB` : p.error}`);
}
