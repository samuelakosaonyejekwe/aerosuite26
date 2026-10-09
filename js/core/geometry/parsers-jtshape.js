// JT tessellation decoding: Shape LOD segments → triangles.
// Written from the published JT File Format Reference (data-compression chapter and the decoding appendices):
// the Int32 compressed data packets (CDP Mk. 1 with Bitlength / Huffman / Arithmetic codecs and probability
// contexts, CDP Mk. 2 with Bitlength, Arithmetic and the Chopper), predictor unpacking, quantised and lossless
// vertex coordinate arrays, the JT 8 tri-strip representation and the JT 9 topologically compressed mesh
// (the "dual VFMesh" polygon-mesh coder).

import { view, inflate } from './parsers-util.js';

/** Sequential little/big-endian reader with bounds checks. */
export class JtReader {
  constructor(b, p, le) { this.b = b; this.dv = view(b); this.p = p; this.le = le; }
  need(n) { if (this.p + n > this.b.length || n < 0) throw new Error('JT shape data ends unexpectedly (truncated or mis-decoded segment)'); }
  u8() { this.need(1); return this.b[this.p++]; }
  i16() { this.need(2); const v = this.dv.getInt16(this.p, this.le); this.p += 2; return v; }
  i32() { this.need(4); const v = this.dv.getInt32(this.p, this.le); this.p += 4; return v; }
  u32() { this.need(4); const v = this.dv.getUint32(this.p, this.le); this.p += 4; return v; }
  f32() { this.need(4); const v = this.dv.getFloat32(this.p, this.le); this.p += 4; return v; }
  skip(n) { this.need(n); this.p += n; }
  words(n) { if (!(n >= 0 && n <= 1e8)) throw new Error('implausible code-text length'); this.need(4 * n); const w = new Uint32Array(n); for (let i = 0; i < n; i++) { w[i] = this.dv.getUint32(this.p, this.le); this.p += 4; } return w; }
}
/** MSB-first bit reader over 32-bit code-text words; reads past the end yield zeros. */
class Bits {
  constructor(w) { this.w = w; this.i = 0; this.n = w.length * 32; }
  bit() { const k = this.i++; return k < this.n ? (this.w[k >>> 5] >>> (31 - (k & 31))) & 1 : 0; }
  u(n) { let v = 0; for (let k = 0; k < n; k++) v = v * 2 + this.bit(); return v; }
  s(n) { if (n === 0) return 0; const v = this.u(n); return n < 32 && v >= 2 ** (n - 1) ? v - 2 ** n : v | 0; }
}
/** MSB-first bit reader over raw bytes (probability context tables); align() skips to the next byte. */
class RawBits {
  constructor(R) { this.R = R; this.cur = 0; this.left = 0; }
  u(n) { let v = 0; for (let k = 0; k < n; k++) { if (!this.left) { this.cur = this.R.u8(); this.left = 8; } v = v * 2 + ((this.cur >> --this.left) & 1); } return v; }
  align() { this.left = 0; }
}

/** How often each codec / representation was used while decoding (reset by the caller); reported in the model metadata. */
export const jtCodecUse = {};
/**
 * Paths that no real file available during development exercised. They follow the published description but are
 * switched off unless the caller opts in (opts.jtUnverified), so that unverified decoding is never returned silently.
 */
export const jtOptions = { allowUnverified: false };
const unverified = (what) => { if (!jtOptions.allowUnverified) throw new Error(`${what} could not be verified against real data and is disabled (pass opts.jtUnverified = true to try it)`); };
const used = (k) => { jtCodecUse[k] = (jtCodecUse[k] || 0) + 1; };

// ---- predictors ----
export const PRED = { Lag1: 0, Lag2: 1, Stride1: 2, Stride2: 3, StripIndex: 4, Ramp: 5, Xor1: 6, Xor2: 7, None: 8 };
export function unpackResiduals(r, pred) {
  if (pred === PRED.None) return r;
  const n = r.length, v = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    if (i < 4) { v[i] = r[i]; continue; }
    const v1 = v[i - 1], v2 = v[i - 2], v4 = v[i - 4]; let p;
    switch (pred) {
      case PRED.Lag2: case PRED.Xor2: p = v2; break; case PRED.Stride1: p = v1 + (v1 - v2); break; case PRED.Stride2: p = v2 + (v2 - v4); break;
      case PRED.StripIndex: p = v2 - v4 < 8 && v2 - v4 > -8 ? v2 + (v2 - v4) : v2 + 2; break; case PRED.Ramp: p = i; break; default: p = v1;
    }
    v[i] = pred === PRED.Xor1 || pred === PRED.Xor2 ? r[i] ^ p : (r[i] + p) | 0;
  }
  return v;
}

// ---- codecs ----
/** Arithmetic decoding of `count` symbols with probability contexts; escape symbols (−2) take the next out-of-band value. */
function arithmetic(bits, contexts, count, oob, symbolCount = count) {
  const out = new Int32Array(count); let low = 0, high = 0xffff, code = bits.u(16), ctx = 0, k = 0, o = 0;
  for (let s = 0; s < symbolCount && k < count; s++) {
    const C = contexts[ctx] || contexts[0], total = C.total; if (!(total > 0)) throw new Error('empty probability context');
    const range = high - low + 1, rescaled = Math.floor(((code - low + 1) * total - 1) / range);
    let e = 0, cum = 0; while (e < C.entries.length - 1 && cum + C.entries[e].count <= rescaled) { cum += C.entries[e].count; e++; }
    const en = C.entries[e];
    if (en.symbol !== -2) out[k++] = en.value; else if (ctx === 0) { if (o >= oob.length) throw new Error('out-of-band values exhausted'); out[k++] = oob[o++]; }
    ctx = en.next;
    high = low + Math.floor((range * (cum + en.count)) / total) - 1; low = low + Math.floor((range * cum) / total);
    for (;;) {
      if (((~(high ^ low)) & 0x8000) !== 0) { /* most significant bits agree: shift them out */ }
      else if ((low & 0x4000) && !(high & 0x4000)) { code ^= 0x4000; low &= 0x3fff; high |= 0x4000; }
      else break;
      low = (low << 1) & 0xffff; high = ((high << 1) | 1) & 0xffff; code = ((code << 1) | bits.bit()) & 0xffff;
    }
  }
  if (k < count) throw new Error('arithmetic code text ended before all values were decoded');
  return out;
}
/** Bitlength codec of CDP Mk. 1: adaptive field width, changed in steps of two bits by a prefix code. */
function bitlength1(bits, totalBits, count) {
  const out = new Int32Array(count); let k = 0, width = 0;
  while (bits.i < totalBits && k < count) {
    if (bits.bit() === 1) { const dir = bits.bit(); for (;;) { width += dir ? 2 : -2; if (width < 0 || width > 32) throw new Error('bitlength field width out of range'); if (bits.bit() !== dir) break; } }
    out[k++] = bits.s(width);
  }
  if (k < count) throw new Error('bitlength code text ended early');
  return out;
}
/** Bitlength codec of CDP Mk. 2: fixed-width, or variable-width runs around a mean. */
function bitlength2(bits, count) {
  const out = new Int32Array(count);
  if (bits.bit() === 0) {
    const nMin = bits.u(6), nMax = bits.u(6), min = bits.s(nMin), max = bits.s(nMax); let range = max - min, width = 0;
    if (range > 0) { width = 1; for (range >>>= 1; range; range >>>= 1) width++; }
    for (let i = 0; i < count; i++) out[i] = (width ? bits.u(width) : 0) + min;
    return out;
  }
  const mean = bits.s(32), chg = bits.u(3), run = bits.u(3), maxDecr = -(2 ** (chg - 1)), maxIncr = 2 ** (chg - 1) - 1; let width = 0, k = 0;
  while (k < count) {
    let d; do { d = bits.s(chg); width += d; if (width < 0 || width > 32) throw new Error('bitlength field width out of range'); } while ((d === maxDecr || d === maxIncr) && chg > 0);
    const len = bits.u(run); if (len === 0 && bits.i >= bits.n) throw new Error('bitlength code text ended early');
    for (let j = 0; j < len && k < count; j++) out[k++] = (bits.s(width) + mean) | 0;
  }
  return out;
}
/** JT 10 "nibbler" integer: 4-bit groups, least significant first, each followed by a continuation bit; sign-extended. */
function nibble(bits) {
  let v = 0, n = 0;
  for (;;) { v += bits.u(4) * 2 ** n; n += 4; if (bits.bit() === 0) break; if (n >= 32) break; }
  if (n >= 32) return v | 0;
  return v >= 2 ** (n - 1) ? v - 2 ** n : v;
}
/** Bitlength codec of the JT 10 Int32 CDP: as Mk. 2 but with nibbler-coded minimum / maximum / mean and fixed 4-bit block fields. */
function bitlength3(bits, count) {
  const out = new Int32Array(count);
  if (bits.bit() === 0) {
    const min = nibble(bits), max = nibble(bits); let range = (max - min) >>> 0, width = 0;
    for (; range; range >>>= 1) width++;
    for (let i = 0; i < count; i++) out[i] = ((width ? bits.u(width) : 0) + min) | 0;
    return out;
  }
  used('JT10 variable-width bitlength');
  const mean = nibble(bits), chg = 4, run = 4, maxDecr = -8, maxIncr = 7; let width = 0, k = 0;
  while (k < count) {
    let d; do { d = bits.s(chg); width += d; if (width < 0 || width > 32) throw new Error('bitlength field width out of range'); } while (d === maxDecr || d === maxIncr);
    const len = bits.u(run); if (len === 0 && bits.i >= bits.n) throw new Error('bitlength code text ended early');
    for (let j = 0; j < len && k < count; j++) out[k++] = (bits.s(width) + mean) | 0;
  }
  return out;
}
/** Probability context of the JT 10 Int32 CDP: entry count, 6-bit count width, 7-bit value width, minimum value; entries carry an escape flag. */
function readContext10(R) {
  const rb = new RawBits(R), n = rb.u(16), occBits = rb.u(6), valBits = rb.u(7), minValue = rb.u(32) | 0, entries = []; let total = 0;
  if (valBits > 32) throw new Error('implausible probability context');
  for (let i = 0; i < n; i++) { const esc = rb.u(1), count = rb.u(occBits), value = (rb.u(valBits) + minValue) | 0; entries.push({ symbol: esc ? -2 : i, count, value, next: 0 }); total += count; }
  rb.align();
  return [{ entries, total }];
}
/**
 * Huffman codec of CDP Mk. 1. The code table is not stored: it is rebuilt from the occurrence counts with the
 * reference min-heap (insertion sifts up while the parent is strictly greater; removal moves the last element
 * down past children that are not greater), the first node taken becoming the "1" branch.
 */
function huffman(bits, contexts, count, oob) {
  const C = contexts[0], heap = [];
  const add = (nd) => { heap.push(nd); let i = heap.length; while (i !== 1 && heap[(i >> 1) - 1].w > nd.w) { heap[i - 1] = heap[(i >> 1) - 1]; i >>= 1; } heap[i - 1] = nd; };
  const top = () => { const first = heap[0], y = heap[heap.length - 1]; let size = heap.length - 1, i = 1, ci = 2; while (ci <= size) { if (ci < size && heap[ci - 1].w > heap[ci].w) ci++; if (y.w < heap[ci - 1].w) break; heap[i - 1] = heap[ci - 1]; i = ci; ci *= 2; } heap[i - 1] = y; heap.pop(); return first; };
  for (const en of C.entries) add({ w: en.count, en });
  if (!heap.length) throw new Error('empty Huffman table');
  while (heap.length > 1) { const one = top(), zero = top(); add({ w: one.w + zero.w, one, zero }); }
  const root = heap[0], out = new Int32Array(count); let o = 0;
  for (let k = 0; k < count; k++) {
    let nd = root; while (nd.one) { nd = bits.bit() ? nd.one : nd.zero; }
    if (bits.i > bits.n + 64) throw new Error('Huffman code text ended early');
    if (nd.en.symbol !== -2) out[k] = nd.en.value; else { if (o >= oob.length) throw new Error('out-of-band values exhausted'); out[k] = oob[o++]; }
  }
  return out;
}

function readContexts(R, mk2, tableCount) {
  const rb = new RawBits(R), ctxs = []; let valBits = 0, minValue = 0;
  for (let t = 0; t < tableCount; t++) {
    const n = rb.u(mk2 ? 16 : 32), symBits = rb.u(6), occBits = rb.u(6); if (t === 0) valBits = rb.u(6); const nextBits = mk2 ? 0 : rb.u(6); if (t === 0) minValue = rb.u(32) | 0;
    if (!(n >= 0 && n <= 1e6)) throw new Error('implausible probability context size');
    const entries = []; let total = 0;
    for (let i = 0; i < n; i++) { const symbol = rb.u(symBits) - 2, count = rb.u(occBits), value = t === 0 ? (rb.u(valBits) + minValue) | 0 : 0, next = nextBits ? rb.u(nextBits) : 0; entries.push({ symbol, count, value, next }); total += count; }
    if (t > 0) { const map = new Map(ctxs[0].entries.map((e) => [e.symbol, e.value])); for (const e of entries) e.value = map.get(e.symbol) ?? 0; }
    ctxs.push({ entries, total });
  }
  rb.align();
  return ctxs;
}

/**
 * Move-to-Front pseudo-codec of the JT 10 Int32 CDP: "window offsets" index a short list of recently used values, an
 * offset of -1 takes the next of the separately coded "window values"; the value used moves to the front of the list.
 */
const jtMtf = { window: 16, move: true };
function moveToFront(values, offsets, count) {
  if (offsets.length !== count) throw new Error('Move-to-Front offset array length mismatch');
  const out = new Int32Array(count), win = []; let nv = 0;
  for (let i = 0; i < count; i++) {
    const k = offsets[i]; let v;
    if (k === -1) { if (nv >= values.length) throw new Error('Move-to-Front window values exhausted'); v = values[nv++]; }
    else { if (!(k >= 0 && k < win.length)) throw new Error('Move-to-Front offset outside the window'); v = win[k]; if (jtMtf.move) win.splice(k, 1); }
    if (k === -1 || jtMtf.move) { win.unshift(v); if (win.length > jtMtf.window) win.pop(); }
    out[i] = v;
  }
  if (nv !== values.length) throw new Error('Move-to-Front window values left over');
  return out;
}
/** Int32 Compressed Data Packet Mk. 2 → residual values. */
export function cdp2(R, depth = 0) {
  if (depth > 8) throw new Error('CDP nesting too deep');
  const count = R.i32(); if (count === 0) return new Int32Array(0);
  if (!(count > 0 && count <= 2e8)) throw new Error('implausible CDP value count');
  const codec = R.u8(); used(`${R.v10 ? 'CDP3' : 'CDP2'} ${['null', 'bitlength', '?', 'arithmetic', 'chopper', 'move-to-front'][codec] ?? codec}`);
  if (codec === 5 && R.v10) {
    const values = cdp2(R, depth + 1), offsets = cdp2(R, depth + 1);
    return moveToFront(values, offsets, count);
  }
  if (codec === 4) {
    const chop = R.u8(); if (chop === 0) return cdp2(R, depth + 1);
    const bias = R.i32(), span = R.u8(), msb = cdp2(R, depth + 1), lsb = cdp2(R, depth + 1), out = new Int32Array(count), sh = span - chop;
    if (msb.length !== count || lsb.length !== count) throw new Error('chopped CDP fields disagree in length');
    for (let i = 0; i < count; i++) out[i] = ((lsb[i] | (msb[i] << sh)) + bias) | 0;
    return out;
  }
  const bitLen = R.i32(); if (!(bitLen >= 0)) throw new Error('negative code-text length');
  const words = R.words((bitLen + 31) >>> 5);
  if (codec === 0) { unverified('the Null codec of the Int32 compressed data packet Mk. 2'); const out = new Int32Array(count); for (let i = 0; i < count && i < words.length; i++) out[i] = words[i] | 0; return out; }
  if (codec === 1) return R.v10 ? bitlength3(new Bits(words), count) : bitlength2(new Bits(words), count);
  if (codec === 3) { const ctx = R.v10 ? readContext10(R) : readContexts(R, true, 1), oob = R.v10 && !ctx[0].entries.some((e) => e.symbol === -2) ? new Int32Array(0) : cdp2(R, depth + 1); return arithmetic(new Bits(words), ctx, count, oob); }   // JT 10 stores out-of-band values only when the table has an escape entry
  throw new Error(`unknown CDP Mk. 2 codec ${codec}`);
}
/** Int32 Compressed Data Packet (Mk. 1, JT 8) → residual values. */
export function cdp1(R, depth = 0) {
  if (depth > 8) throw new Error('CDP nesting too deep');
  const codec = R.u8(); used(`CDP1 ${['null', 'bitlength', 'huffman', 'arithmetic'][codec] ?? codec}`);
  if (codec === 0) { unverified('the Null codec of the Int32 compressed data packet'); const n = R.i32(); const w = R.words(n); return Int32Array.from(w, (x) => x | 0); }
  if (codec !== 1 && codec !== 2 && codec !== 3) throw new Error(`unknown CDP codec ${codec}`);
  let ctxs = null, oob = new Int32Array(0), tables = 0;
  if (codec >= 2) { tables = R.u8(); if (!(tables >= 1 && tables <= 8)) throw new Error('implausible probability context count'); ctxs = readContexts(R, false, tables); const nOob = R.i32(); if (nOob > 0) oob = cdp1(R, depth + 1); }
  const bitLen = R.i32(), count = R.i32(); if (!(count >= 0 && count <= 2e8 && bitLen >= 0)) throw new Error('implausible CDP value count');
  const symbols = tables > 1 ? R.i32() : count, n = R.i32(), words = R.words(n), bits = new Bits(words);
  if (codec === 1) return bitlength1(bits, bitLen, count);
  if (codec === 2) return huffman(bits, ctxs, count, oob);
  return arithmetic(bits, ctxs, count, oob, symbols);
}

// ---- vertex coordinates ----
function quantizers(R) { const q = []; for (let k = 0; k < 3; k++) q.push({ min: R.f32(), max: R.f32(), bits: R.u8() }); return q; }
const dequant = (code, q) => { const maxCode = q.bits < 32 ? 2 ** q.bits - 1 : 0xffffffff; return q.min + ((code >>> 0) * (q.max - q.min)) / (maxCode || 1); };
const F32 = new Float32Array(1), U32 = new Uint32Array(F32.buffer);
/** JT 9 "Compressed Vertex Coordinate Array": quantised codes or exponent/mantissa pairs, Lag1-predicted. */
function coords9(R) {
  const n = R.i32(), nc = R.u8(), q = quantizers(R); if (!(n >= 0 && n <= 5e7 && nc >= 2 && nc <= 4)) throw new Error('implausible vertex coordinate array');
  const out = new Float64Array(3 * n);
  for (let c = 0; c < nc; c++) {
    if (R.v10) used(q[0].bits > 0 ? 'JT10 quantised coordinates' : 'JT10 lossless coordinates'); else used(q[0].bits > 0 ? 'JT9 quantised coordinates' : 'JT9 lossless coordinates');
    if (R.v10 && q[0].bits === 0) {                          // JT 10: the IEEE-754 bit pattern of each coordinate as one Lag1-predicted integer
      const raw = unpackResiduals(cdp2(R), PRED.Lag1); if (raw.length !== n) throw new Error('vertex coordinate array length mismatch');
      if (c < 3) for (let i = 0; i < n; i++) { U32[0] = raw[i] >>> 0; out[3 * i + c] = F32[0]; }
      continue;
    }
    if (q[0].bits > 0) { const codes = unpackResiduals(cdp2(R), PRED.Lag1); if (codes.length !== n) throw new Error('vertex code array length mismatch'); if (c < 3) for (let i = 0; i < n; i++) out[3 * i + c] = dequant(codes[i], q[c]); }
    else { const ex = unpackResiduals(cdp2(R), PRED.Lag1), ma = unpackResiduals(cdp2(R), PRED.Lag1); if (ex.length !== n || ma.length !== n) throw new Error('vertex exponent/mantissa length mismatch'); if (c < 3) for (let i = 0; i < n; i++) { U32[0] = ((ex[i] << 23) | ma[i]) >>> 0; out[3 * i + c] = F32[0]; } }
  }
  R.i32();                                                 // hash
  return out;
}

// ---- JT 9 topologically compressed mesh (dual VFMesh decoder) ----
function decodeDualMesh(deg, val, grp, flags, masks, large, splitFace, splitPos) {
  // dual vertices = triangles of the model (valence 3); dual faces = model vertices (degree = incident triangles)
  const vVal = [], vGrp = [], vFlag = [], vFI = [], vf = [], fDeg = [], fEmpty = [], fVI = [], fv = [];
  const pos = { deg: new Int32Array(8), val: 0, grp: 0, flag: 0, mask: new Int32Array(8), large: 0, sf: 0, sp: 0 };
  const active = [], removed = new Set();
  const valence = (v) => vVal[v], degree = (f) => fDeg[f], face = (v, s) => vf[vFI[v] + s], vtx = (f, s) => fv[fVI[f] + s];
  const setVtxFace = (v, s, f) => { vf[vFI[v] + s] = f; };
  const setFaceVtx = (f, s, v) => { if (fv[fVI[f] + s] !== v) { if (fv[fVI[f] + s] === -1) fEmpty[f]--; fv[fVI[f] + s] = v; } };
  const findVtxSlot = (f, t) => { for (let s = 0; s < fDeg[f]; s++) if (fv[fVI[f] + s] === t) return s; return -1; };
  const findFaceSlot = (v, t) => { for (let s = 0; s < vVal[v]; s++) if (vf[vFI[v] + s] === t) return s; return -1; };
  const inc = (i, n) => (i + 1) % n, dec = (i, n) => (i + n - 1) % n;
  const ioVtx = () => {
    if (pos.val >= val.length) return -1;
    const v = vVal.length, cv = val[pos.val++]; if (!(cv >= 0 && cv <= 65535)) throw new Error('implausible valence');
    vVal.push(cv); vFI.push(vf.length); for (let i = 0; i < cv; i++) vf.push(-1);
    vGrp.push(pos.grp < grp.length ? grp[pos.grp++] : -1); vFlag.push(pos.flag < flags.length ? flags[pos.flag++] : 0);
    return v;
  };
  const faceContext = (v) => {
    const cv = vVal[v]; let known = 0, tot = 0; for (let i = 0; i < cv; i++) { const f = face(v, i); if (f >= 0) { known++; tot += fDeg[f]; } }
    if (cv === 3) return tot < known * 6 ? 0 : tot === known * 6 ? 1 : 2;
    if (cv === 4) return tot < known * 4 ? 3 : tot === known * 4 ? 4 : 5;
    return cv === 5 ? 6 : 7;
  };
  const ioFace = (v) => {
    const c = faceContext(v), sym = pos.deg[c] < deg[c].length ? deg[c][pos.deg[c]++] : -1;
    if (sym === 0) return -1;                               // split: the face already exists
    if (sym < 0) throw new Error('face-degree symbols exhausted');
    const f = fDeg.length; if (sym > 1e6) throw new Error('implausible face degree');
    fDeg.push(sym); fEmpty.push(sym); fVI.push(fv.length); for (let i = 0; i < sym; i++) fv.push(-1);
    if (sym <= 64) { const g = Math.min(7, Math.max(0, sym - 2)); if (pos.mask[g] < masks[g].length) pos.mask[g]++; } else pos.large += (sym + 31) >>> 5;   // attribute masks are consumed but not needed for positions
    return f;
  };
  const addVtxToFace = (v, jFSlot, f, iVSlot) => {
    const cw = dec(iVSlot, fDeg[f]), ccw = inc(iVSlot, fDeg[f]);
    setFaceVtx(f, iVSlot, v);
    const fp = vtx(f, cw);
    if (fp !== -1) { let ip = findFaceSlot(fp, f); const s = inc(jFSlot, vVal[v]); if (ip >= 0 && face(v, s) === -1) { ip = dec(ip, vVal[fp]); setVtxFace(v, s, face(fp, ip)); } }
    const fn = vtx(f, ccw);
    if (fn !== -1) { let ix = findFaceSlot(fn, f); const s = dec(jFSlot, vVal[v]); if (ix >= 0 && face(v, s) === -1) { ix = inc(ix, vVal[fn]); setVtxFace(v, s, face(fn, ix)); } }
  };
  const activateF = (v, slot) => {
    let f = ioFace(v);
    if (f >= 0) { setVtxFace(v, slot, f); setFaceVtx(f, 0, v); active.push(f); }
    else {
      const off = pos.sf < splitFace.length ? splitFace[pos.sf++] : -1, sp = pos.sp < splitPos.length ? splitPos[pos.sp++] : -1;
      if (!(off >= 1 && off <= active.length) || sp < 0) throw new Error('invalid split symbol');
      f = active[active.length - off]; if (sp >= fDeg[f]) throw new Error('invalid split position');
      setVtxFace(v, slot, f); addVtxToFace(v, slot, f, sp);
    }
    return f;
  };
  const activateV = (f, slot) => { const v = ioVtx(); if (v < 0) throw new Error('valence symbols exhausted'); if (vVal[v] < 1) throw new Error('vertex without faces'); setVtxFace(v, 0, f); addVtxToFace(v, 0, f, slot); return v; };
  const completeV = (v, slot) => {
    const cv = vVal[v]; let vp = face(v, 0), jp = slot, i = 1, vn;
    while (i < cv && (vn = face(v, i)) !== -1) {             // faces already known in the "next" direction
      jp = dec(jp, fDeg[vp]); const v2 = vtx(vp, jp); if (v2 === -1) break;
      let jn = findVtxSlot(vn, v2); if (jn < 0) throw new Error('inconsistent mesh topology'); jn = dec(jn, fDeg[vn]);
      addVtxToFace(v, i, vn, jn); vp = vn; jp = jn; i++;
    }
    if (i >= cv) return;
    const ilast = i; vp = face(v, 0); jp = slot; i = cv - 1;
    while (i >= ilast && (vn = face(v, i)) !== -1) {         // … and in the "previous" direction
      jp = inc(jp, fDeg[vp]); const v2 = vtx(vp, jp); if (v2 === -1) break;
      let jn = findVtxSlot(vn, v2); if (jn < 0) throw new Error('inconsistent mesh topology'); jn = inc(jn, fDeg[vn]);
      addVtxToFace(v, i, vn, jn); vp = vn; jp = jn; i--;
    }
    for (let k = ilast; k <= i; k++) activateF(v, k);
  };
  const nextActive = () => {
    while (active.length && removed.has(active[active.length - 1])) active.pop();
    let best = -1, lowest = Infinity;
    for (let i = active.length - 1; i >= Math.max(0, active.length - 16); i--) { const f = active[i]; if (removed.has(f)) { active.splice(i, 1); continue; } if (fEmpty[f] < lowest) { lowest = fEmpty[f]; best = f; } }
    return best;
  };
  for (let guardN = 0; guardN < 1e7; guardN++) {
    const v0 = ioVtx(); if (v0 < 0) break;
    for (let i = 0; i < vVal[v0]; i++) activateF(v0, i);
    for (let f; (f = nextActive()) !== -1;) {
      for (let s, g = 0; (s = findVtxSlot(f, -1)) !== -1; g++) { if (g > 1e6) throw new Error('mesh decoder did not terminate'); const v = activateV(f, s); completeV(v, s); }
      removed.add(f);
    }
  }
  if (pos.val !== val.length || pos.sf !== splitFace.length || pos.sp !== splitPos.length || deg.some((d, c) => pos.deg[c] !== d.length)) throw new Error('topology symbols were not all consumed (mesh decoded inconsistently)');
  const tri = [], triGroup = [];
  for (let v = 0; v < vVal.length; v++) { if (vGrp[v] < 0) continue; const n = vVal[v]; for (let k = 2; k < n; k++) { const a = face(v, 0), b = face(v, k - 1), c = face(v, k); if (a >= 0 && b >= 0 && c >= 0) { tri.push(a, b, c); triGroup.push(vGrp[v]); } } }
  return { tri, vertices: fDeg.length, groups: triGroup };
}

const TRISTRIP_LOD = '10dd10ab-2ac8-11d1-9b6b0080c7bb5997';
/**
 * Decode one Shape LOD segment payload (the bytes after the segment header) of a tri-strip set.
 * Returns { positions: Float64Array, triangles: number[] } in the part's local coordinates.
 */
const TOPO_REP_V10 = 'f830a5ad-be4c-4fbc-9b5fb9269278d2e1';
const guidAt = (b, p, le) => { const dv = new DataView(b.buffer, b.byteOffset, b.byteLength), h = (v, n) => v.toString(16).padStart(n, '0'); let s = h(dv.getUint32(p, le), 8) + '-' + h(dv.getUint16(p + 4, le), 4) + '-' + h(dv.getUint16(p + 6, le), 4) + '-'; for (let i = 8; i < 16; i++) s += h(b[p + i], 2); return s; };
/** Triangles as a typed array; any triangle that points outside the vertex array is dropped. */
function packLod(positions, tri, groups) {
  const nv = positions.length / 3, T = new Uint32Array(tri.length), G = groups ? new Int32Array(groups.length) : null; let o = 0;
  for (let i = 0; i + 2 < tri.length; i += 3) { const a = tri[i], b = tri[i + 1], c = tri[i + 2]; if (!(a >= 0 && a < nv && b >= 0 && b < nv && c >= 0 && c < nv)) continue; T[o] = a; T[o + 1] = b; T[o + 2] = c; if (G) G[o / 3] = groups[i / 3]; o += 3; }
  return { positions, triangles: o === T.length ? T : T.slice(0, o), faceGroups: G ? (o === T.length ? G : G.slice(0, o / 3)) : null };
}
export async function decodeShapeLOD(payload, le, major, typeId) {
  if (typeId !== TRISTRIP_LOD) throw new Error('not a tri-strip set shape (polyline, point and primitive sets carry no triangles)');
  const R = new JtReader(payload, 4 + 16 + 1 + 4, le);       // element length, object type, base type, object id
  if (major >= 9) {
    if (major >= 10) {
      R.v10 = true;
      // JT 10: U8 version fields, and the topologically compressed representation is a nested element of its own
      R.u8(); R.u8(); R.skip(8);                            // base shape LOD and vertex shape LOD versions, vertex bindings
      R.i32(); const rep = guidAt(payload, R.p, le); R.skip(16);
      if (rep !== TOPO_REP_V10) throw new Error(`unsupported JT 10 mesh representation (${rep})`);
      R.u8(); R.i32();                                      // base type, object id
      R.u8(); R.skip(4); R.u8();                            // version, vertex records reference, version
      used('JT10 TopoMesh element');
    } else {
      R.i16();                                              // base shape LOD data version
      R.i16(); R.skip(8);                                   // vertex shape LOD data: version, vertex bindings
      R.i16(); R.i32();                                     // TopoMesh LOD data: version, vertex records object id
      R.i16();                                              // TopoMesh topologically compressed LOD data version
    }
    const deg = []; for (let i = 0; i < 8; i++) deg.push(cdp2(R));
    const val = cdp2(R), grp = cdp2(R), flags = unpackResiduals(cdp2(R), PRED.Lag1), masks = []; for (let i = 0; i < 8; i++) masks.push(cdp2(R));
    cdp2(R); if (major < 10) cdp2(R);                       // attribute mask words above bit 31 of context 7 (two packets in JT 9, one in JT 10)
    const large = cdp2(R), splitFace = unpackResiduals(cdp2(R), PRED.Lag1), splitPos = cdp2(R);
    R.u32();                                                // composite hash
    const mesh = decodeDualMesh(deg, val, grp, flags, masks, large, splitFace, splitPos);
    R.skip(8); R.skip(4);                                   // vertex bindings, quantisation parameters
    const nv = R.i32(); if (nv === 0) return { positions: new Float64Array(0), triangles: [] };
    R.i32();                                                // number of vertex attribute records
    const positions = coords9(R);
    if (positions.length / 3 < mesh.vertices) throw new Error('fewer vertex coordinates than mesh vertices');
    return packLod(positions, mesh.tri, mesh.groups);
  }
  // JT 8: vertex-based shape, tri-strips as primitive index lists
  R.skip(8);                                                // vertex shape LOD data: bindings and quantisation parameters (repeated below)
  R.i16();                                                  // vertex-based shape compressed rep data version
  const normalBinding = R.u8(); R.u8(); R.u8(); const bitsVertex = R.u8(); R.skip(3);
  const prim = unpackResiduals(cdp1(R), PRED.Stride1); if (prim.length < 2) return { positions: new Float64Array(0), triangles: [] };
  const nVerts = prim[prim.length - 1]; if (!(nVerts >= 0 && nVerts <= 5e7)) throw new Error('implausible vertex count');
  let positions, index = null;
  used(bitsVertex === 0 ? 'JT8 lossless vertex data' : 'JT8 quantised vertex data');
  if (bitsVertex === 0) {
    const usize = R.i32(), csize = R.i32(); let D = R, raw = null;
    if (csize > 0) { R.need(csize); raw = await inflate(payload.subarray(R.p, R.p + csize), Math.max(usize, 16) + 1024, 'deflate'); R.skip(csize); D = new JtReader(raw, 0, le); }
    positions = new Float64Array(3 * nVerts);
    for (let i = 0; i < nVerts; i++) { if (normalBinding) D.skip(12); positions[3 * i] = D.f32(); positions[3 * i + 1] = D.f32(); positions[3 * i + 2] = D.f32(); }
  } else {
    const q = quantizers(R), n = R.i32(); if (!(n >= 0 && n <= 5e7)) throw new Error('implausible vertex count');
    positions = new Float64Array(3 * n);
    for (let c = 0; c < 3; c++) { const codes = unpackResiduals(cdp1(R), PRED.Lag1); if (codes.length !== n) throw new Error('vertex code array length mismatch'); for (let i = 0; i < n; i++) positions[3 * i + c] = dequant(codes[i], q[c]); }
    if (normalBinding) { R.u8(); const nn = R.i32(); for (let c = 0; c < 4; c++) { const a = cdp1(R); if (a.length !== nn) throw new Error('normal code array length mismatch'); } }
    index = unpackResiduals(cdp1(R), PRED.StripIndex);
  }
  const tri = [];
  for (let s = 0; s + 1 < prim.length; s++) for (let o = prim[s], k = 0; o + 2 < prim[s + 1]; o++, k++) { const a = o, b = o + (k % 2 ? 2 : 1), c = o + (k % 2 ? 1 : 2); const i0 = index ? index[a] : a, i1 = index ? index[b] : b, i2 = index ? index[c] : c; if (i0 === undefined || i1 === undefined || i2 === undefined || i0 === i1 || i1 === i2 || i0 === i2) continue; tri.push(i0, i1, i2); }   // repeated indices are strip stitching, not triangles
  return packLod(positions, tri, null);
}
