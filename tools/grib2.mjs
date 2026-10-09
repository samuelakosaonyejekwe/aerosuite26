// Minimal GRIB edition 2 reader for the fields the snapshot job takes from NOAA model output (GFS, GFS-Wave,
// GEFS-Aerosols): regular latitude/longitude grids (template 3.0) packed with simple packing (5.0), complex
// packing (5.2) or complex packing with spatial differencing (5.3), with or without a bitmap. Anything else
// (JPEG 2000, PNG, spectral data, other grids) is refused with a clear error rather than decoded wrongly.
// Reference: WMO Manual on Codes, FM 92 GRIB edition 2; NCEP g2 library (comunpack) for the group layout.

const sm16 = (b, o) => { const v = b.readUInt16BE(o); return v & 0x8000 ? -(v & 0x7fff) : v; };          // sign-and-magnitude integers
const sm32 = (b, o) => { const v = b.readUInt32BE(o); return v & 0x80000000 ? -(v & 0x7fffffff) : v; };
const smN = (b, o, n) => { let v = 0; for (let i = 0; i < n; i++) v = v * 256 + b[o + i]; const top = 2 ** (8 * n - 1); return v >= top ? -(v - top) : v; };

/** Sequential big-endian bit reader. */
class Bits {
  constructor(buf, pos = 0) { this.b = buf; this.p = pos * 8; }
  read(n) { // n ≤ 32
    if (n === 0) return 0;
    let v = 0, p = this.p; const b = this.b;
    for (let left = n; left > 0;) { const byte = b[p >> 3], avail = 8 - (p & 7), take = Math.min(avail, left); v = v * 2 ** take + ((byte >> (avail - take)) & ((1 << take) - 1)); p += take; left -= take; }
    this.p = p; return v;
  }
  align() { this.p = (this.p + 7) & ~7; }
}

/** Split a buffer holding one or more GRIB2 messages into message descriptors (no data decoded yet). */
export function gribMessages(buf) {
  const out = []; let p = 0;
  while (p + 16 <= buf.length) {
    if (buf.toString('latin1', p, p + 4) !== 'GRIB') { p++; continue; }
    if (buf[p + 7] !== 2) throw new Error('not GRIB edition 2');
    const total = Number(buf.readBigUInt64BE(p + 8)), end = p + total, discipline = buf[p + 6];
    let q = p + 16, grid = null, prod = null, drs = null, bitmap = null, refTime = null;
    while (q < end - 4) {
      const len = buf.readUInt32BE(q), sec = buf[q + 4], s = buf.subarray(q, q + len);
      if (sec === 1) refTime = Date.UTC(s.readUInt16BE(12), s[14] - 1, s[15], s[16], s[17], s[18]);
      else if (sec === 3) {
        const tpl = s.readUInt16BE(12); if (tpl !== 0) throw new Error(`grid template 3.${tpl} is not supported`);
        grid = { ni: s.readUInt32BE(30), nj: s.readUInt32BE(34), la1: sm32(s, 46) / 1e6, lo1: sm32(s, 50) / 1e6, la2: sm32(s, 55) / 1e6, lo2: sm32(s, 59) / 1e6, di: s.readUInt32BE(63) / 1e6, dj: s.readUInt32BE(67) / 1e6, scan: s[71] };
        if (grid.scan & 0x30) throw new Error('grid scanning mode is not supported');
      } else if (sec === 4) {
        const tpl = s.readUInt16BE(7), unit = s[17], ft = s.readUInt32BE(18), hours = { 0: 1 / 60, 1: 1, 2: 24, 10: 3, 11: 6, 12: 12, 13: 1 / 3600 }[unit];
        prod = { template: tpl, discipline, category: s[9], number: s[10], forecast_h: hours == null ? null : ft * hours, levelType: s[22], level: s[23] === 255 ? null : sm32(s, 24) / 10 ** (s[23] > 127 ? -(s[23] & 127) : s[23]) };
        if (tpl === 48 || tpl === 44) { // aerosol products: the time and level fields sit after the aerosol description
          const o = tpl === 48 ? 41 : 30, h2 = { 0: 1 / 60, 1: 1, 2: 24, 10: 3, 11: 6, 12: 12, 13: 1 / 3600 }[s[o]];
          Object.assign(prod, { aerosol: s.readUInt16BE(11), forecast_h: h2 == null ? null : s.readUInt32BE(o + 1) * h2, levelType: s[o + 5], level: s[o + 6] === 255 ? null : sm32(s, o + 7) / 10 ** (s[o + 6] > 127 ? -(s[o + 6] & 127) : s[o + 6]), sizeType: s[13], sizeLo: sm32(s, 15) * 10 ** -s[14], sizeHi: sm32(s, 20) * 10 ** -s[19] });
          if (tpl === 48) prod.wavelength = sm32(s, 26) * 10 ** -s[25];
        }
        if (![0, 1, 8, 11, 44, 48].includes(tpl)) prod.unsupported = true;
      } else if (sec === 5) {
        const tpl = s.readUInt16BE(9); drs = { n: s.readUInt32BE(5), template: tpl, R: s.readFloatBE(11), E: sm16(s, 15), D: sm16(s, 17), nbits: s[19] };
        if (tpl === 2 || tpl === 3) Object.assign(drs, { missMgmt: s[22], ng: s.readUInt32BE(31), refW: s[35], bitsW: s[36], refL: s.readUInt32BE(37), incL: s[41], lastL: s.readUInt32BE(42), bitsL: s[46] });
        if (tpl === 3) Object.assign(drs, { order: s[47], extra: s[48] });
      } else if (sec === 6) { const ind = s[5]; if (ind === 0) bitmap = s.subarray(6); else if (ind !== 255) throw new Error('predefined bitmaps are not supported'); else bitmap = null; }
      else if (sec === 7) { out.push({ ...prod, refTime, grid, drs, bitmap, data: s.subarray(5) }); }
      q += len;
    }
    p = end;
  }
  return out;
}

/** Decode one message → Float32Array of ni × nj values in the file's scan order (NaN where the bitmap or the packing marks a value missing). */
export function gribValues(m) {
  const { drs, data, grid } = m, n = drs.n, total = grid.ni * grid.nj, vals = new Float64Array(n), miss = new Uint8Array(n);
  if (drs.template === 0) { const br = new Bits(data); if (drs.nbits) for (let i = 0; i < n; i++) vals[i] = br.read(drs.nbits); }
  else if (drs.template === 2 || drs.template === 3) {
    let pos = 0, v1 = 0, v2 = 0, minsd = 0;
    if (drs.template === 3) { if (drs.order !== 1 && drs.order !== 2) throw new Error('spatial differencing order is not supported'); const k = drs.extra; if (k) { v1 = smN(data, 0, k); pos = k; if (drs.order === 2) { v2 = smN(data, pos, k); pos += k; } minsd = smN(data, pos, k); pos += k; } }
    const br = new Bits(data, pos), ng = drs.ng, gref = new Float64Array(ng), gw = new Uint8Array(ng), gl = new Uint32Array(ng);
    for (let j = 0; j < ng; j++) gref[j] = br.read(drs.nbits); br.align();
    for (let j = 0; j < ng; j++) gw[j] = br.read(drs.bitsW) + drs.refW; br.align();
    for (let j = 0; j < ng; j++) gl[j] = br.read(drs.bitsL) * drs.incL + drs.refL; br.align();
    gl[ng - 1] = drs.lastL;
    const mm = drs.missMgmt, refMiss = 2 ** drs.nbits - 1;
    if (mm > 2) throw new Error('missing-value management is not supported');
    let i = 0;
    for (let j = 0; j < ng; j++) {
      const w = gw[j], len = gl[j], top = w ? 2 ** w - 1 : 0; if (i + len > n) throw new Error('group lengths exceed the number of values');
      for (let k = 0; k < len; k++, i++) {
        if (w === 0) { if (mm && gref[j] >= refMiss - (mm - 1)) miss[i] = 1; else vals[i] = gref[j]; }
        else { const x = br.read(w); if (mm && x >= top - (mm - 1)) miss[i] = 1; else vals[i] = gref[j] + x; }
      }
    }
    if (i !== n) throw new Error('group lengths do not add up to the number of values');
    if (drs.template === 3) { // undo the differencing over the non-missing values
      let a = 0, b = 0, seen = 0;
      for (let k = 0; k < n; k++) {
        if (miss[k]) continue;
        if (seen === 0) vals[k] = v1; else if (seen === 1 && drs.order === 2) vals[k] = v2;
        else vals[k] = drs.order === 1 ? vals[k] + minsd + a : vals[k] + minsd + 2 * a - b;
        b = a; a = vals[k]; seen++;
      }
    }
  } else throw new Error(`data packing 5.${drs.template} is not supported`);
  const f = 2 ** drs.E, d = 10 ** -drs.D, out = new Float32Array(total);
  if (m.bitmap) { const bm = m.bitmap; let k = 0; for (let i = 0; i < total; i++) { if ((bm[i >> 3] >> (7 - (i & 7))) & 1) { out[i] = miss[k] ? NaN : (drs.R + vals[k] * f) * d; k++; } else out[i] = NaN; } if (k !== n) throw new Error('bitmap does not match the number of packed values'); }
  else { if (n !== total) throw new Error('number of values does not match the grid'); for (let i = 0; i < total; i++) out[i] = miss[i] ? NaN : (drs.R + vals[i] * f) * d; }
  return out;
}

/** Value access on a decoded global lat/lon field: (lat index from the north, lon index from 0°E), whatever the file's row order. */
export function gridReader(m, values) {
  const g = m.grid, north = !(g.scan & 0x40), di = g.di || 360 / g.ni, dj = g.dj || 180 / (g.nj - 1), lo1 = ((g.lo1 % 360) + 360) % 360, latN = north ? g.la1 : g.la2;
  if (g.scan & 0x80) throw new Error('west-going rows are not supported');
  const at = (j, i) => values[(north ? j : g.nj - 1 - j) * g.ni + i];
  /** Bilinear value at a latitude/longitude; NaN neighbours are left out (NaN if all four are missing). */
  const sample = (lat, lon) => {
    const y = Math.min(g.nj - 1, Math.max(0, (latN - lat) / dj)), x = ((((lon - lo1) % 360) + 360) % 360) / di, j0 = Math.min(g.nj - 2, Math.floor(y)), i0 = Math.floor(x) % g.ni, fy = y - j0, fx = x - Math.floor(x), i1 = (i0 + 1) % g.ni;
    if (i0 + 1 >= g.ni && g.ni * di < 359.5) return at(j0, i0);
    let s = 0, w = 0; for (const [v, k] of [[at(j0, i0), (1 - fy) * (1 - fx)], [at(j0, i1), (1 - fy) * fx], [at(j0 + 1, i0), fy * (1 - fx)], [at(j0 + 1, i1), fy * fx]]) if (v === v && k > 0) { s += v * k; w += k; }
    return w > 0 ? s / w : NaN;
  };
  return { at, sample, latN, di, dj };
}
