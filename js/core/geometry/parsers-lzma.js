// XZ container / LZMA2 decoder (JT 10 compresses its segments this way). Written from the public .xz file format
// and LZMA specifications: stream header, block headers with a single LZMA2 filter, LZMA2 chunks, and the LZMA
// range decoder with literal / match / repeat-match models. Integrity checks (CRC) are skipped, not verified.

const PROB_INIT = 1024;
class Lzma {
  constructor() { this.out = []; this.setProps(3, 0, 2); }
  setProps(lc, lp, pb) { this.lc = lc; this.lp = lp; this.pb = pb; this.resetState(); }
  resetState() {
    const mk = (n) => new Uint16Array(n).fill(PROB_INIT);
    this.isMatch = mk(12 << 4); this.isRep = mk(12); this.isRepG0 = mk(12); this.isRepG1 = mk(12); this.isRepG2 = mk(12); this.isRep0Long = mk(12 << 4);
    this.lit = mk(0x300 << (this.lc + this.lp)); this.posSlot = mk(4 << 6); this.posDec = mk(115); this.align = mk(16);
    this.len = { choice: mk(2), low: mk(16 << 3), mid: mk(16 << 3), high: mk(256) }; this.repLen = { choice: mk(2), low: mk(16 << 3), mid: mk(16 << 3), high: mk(256) };
    this.state = 0; this.rep0 = this.rep1 = this.rep2 = this.rep3 = 0;
  }
  initRange(src, p, end) { if (p + 5 > end) throw new Error('LZMA chunk is truncated'); this.src = src; this.end = end; this.range = 0xffffffff; this.code = ((src[p + 1] << 24) | (src[p + 2] << 16) | (src[p + 3] << 8) | src[p + 4]) >>> 0; this.p = p + 5; }
  norm() { if (this.range < 0x1000000) { this.range = (this.range << 8) >>> 0; this.code = ((this.code << 8) | (this.p < this.end ? this.src[this.p] : 0)) >>> 0; this.p++; } }
  bit(probs, i) {
    const pr = probs[i], bound = (this.range >>> 11) * pr;
    if (this.code < bound) { this.range = bound; probs[i] = pr + ((2048 - pr) >> 5); this.norm(); return 0; }
    this.range = (this.range - bound) >>> 0; this.code = (this.code - bound) >>> 0; probs[i] = pr - (pr >> 5); this.norm(); return 1;
  }
  tree(probs, off, n) { let m = 1; for (let i = 0; i < n; i++) m = (m << 1) | this.bit(probs, off + m); return m - (1 << n); }
  rtree(probs, off, n) { let m = 1, s = 0; for (let i = 0; i < n; i++) { const b = this.bit(probs, off + m); m = (m << 1) | b; s |= b << i; } return s; }
  direct(n) { let r = 0; for (; n > 0; n--) { this.range >>>= 1; let t = 0; if (this.code >= this.range) { this.code = (this.code - this.range) >>> 0; t = 1; } r = ((r << 1) | t) >>> 0; this.norm(); } return r; }
  length(L, ps) { if (!this.bit(L.choice, 0)) return this.tree(L.low, ps << 3, 3); if (!this.bit(L.choice, 1)) return 8 + this.tree(L.mid, ps << 3, 3); return 16 + this.tree(L.high, 0, 8); }
  /** Decode `count` bytes into this.out. */
  run(count, cap) {
    const out = this.out, target = out.length + count; if (target > cap) throw new Error('LZMA data expands beyond the import limit');
    const pbMask = (1 << this.pb) - 1, lpMask = (1 << this.lp) - 1;
    while (out.length < target) {
      const pos = out.length, ps = pos & pbMask;
      if (!this.bit(this.isMatch, (this.state << 4) + ps)) {
        const prev = pos ? out[pos - 1] : 0, base = 0x300 * (((pos & lpMask) << this.lc) + (prev >> (8 - this.lc))); let sym = 1;
        if (this.state >= 7) { let mb = out[pos - this.rep0 - 1] ?? 0; do { const m = (mb >> 7) & 1; mb <<= 1; const bt = this.bit(this.lit, base + ((1 + m) << 8) + sym); sym = (sym << 1) | bt; if (m !== bt) break; } while (sym < 0x100); }
        while (sym < 0x100) sym = (sym << 1) | this.bit(this.lit, base + sym);
        out.push(sym & 0xff); this.state = this.state < 4 ? 0 : this.state < 10 ? this.state - 3 : this.state - 6; continue;
      }
      let len;
      if (!this.bit(this.isRep, this.state)) {
        this.rep3 = this.rep2; this.rep2 = this.rep1; this.rep1 = this.rep0; len = this.length(this.len, ps); this.state = this.state < 7 ? 7 : 10;
        const slot = this.tree(this.posSlot, Math.min(len, 3) << 6, 6);
        if (slot < 4) this.rep0 = slot;
        else { const nb = (slot >> 1) - 1; let dist = ((2 | (slot & 1)) << nb) >>> 0; if (slot < 14) dist += this.rtree(this.posDec, dist - slot, nb); else { dist = (dist + ((this.direct(nb - 4) << 4) >>> 0)) >>> 0; dist += this.rtree(this.align, 0, 4); } this.rep0 = dist >>> 0; }
        if (this.rep0 === 0xffffffff) return;
      } else {
        if (!this.bit(this.isRepG0, this.state)) { if (!this.bit(this.isRep0Long, (this.state << 4) + ps)) { this.state = this.state < 7 ? 9 : 11; if (this.rep0 >= pos) throw new Error('LZMA data is corrupt (distance before the start)'); out.push(out[pos - this.rep0 - 1]); continue; } }
        else { let dist; if (!this.bit(this.isRepG1, this.state)) dist = this.rep1; else { if (!this.bit(this.isRepG2, this.state)) dist = this.rep2; else { dist = this.rep3; this.rep3 = this.rep2; } this.rep2 = this.rep1; } this.rep1 = this.rep0; this.rep0 = dist; }
        len = this.length(this.repLen, ps); this.state = this.state < 7 ? 8 : 11;
      }
      len += 2; if (this.rep0 >= out.length) throw new Error('LZMA data is corrupt (distance before the start)');
      for (let i = 0, from = out.length - this.rep0 - 1; i < len && out.length < target; i++) out.push(out[from + i]);
    }
  }
}

/** Decode the LZMA2 chunk sequence starting at src[p]; returns { bytes, end }. */
function lzma2(src, p, cap) {
  const d = new Lzma(); let haveProps = false;
  for (let g = 0; g < 1e7; g++) {
    if (p >= src.length) throw new Error('LZMA2 stream is truncated');
    const ctl = src[p++]; if (ctl === 0) break;
    if (ctl === 1 || ctl === 2) { const n = ((src[p] << 8) | src[p + 1]) + 1; p += 2; if (p + n > src.length || d.out.length + n > cap) throw new Error('LZMA2 uncompressed chunk is truncated'); for (let i = 0; i < n; i++) d.out.push(src[p + i]); p += n; continue; }
    if (ctl < 0x80) throw new Error('invalid LZMA2 control byte');
    const unpacked = (((ctl & 0x1f) << 16) | (src[p] << 8) | src[p + 1]) + 1, packed = ((src[p + 2] << 8) | src[p + 3]) + 1, reset = (ctl >> 5) & 3; p += 4;
    if (reset >= 2) { const pr = src[p++]; if (pr >= 225) throw new Error('invalid LZMA properties'); const lc = pr % 9, r = Math.floor(pr / 9); d.setProps(lc, r % 5, Math.floor(r / 5)); haveProps = true; }
    else if (reset === 1) d.resetState();
    if (!haveProps) throw new Error('LZMA2 chunk without properties');
    if (p + packed > src.length) throw new Error('LZMA2 chunk is truncated');
    d.initRange(src, p, p + packed); d.run(unpacked, cap); p += packed;
  }
  return { bytes: d.out, end: p };
}

/** Decompress an .xz stream (LZMA2 filter only). */
export function unxz(src, cap = 512e6) {
  if (src.length < 12 || src[0] !== 0xfd || src[1] !== 0x37 || src[2] !== 0x7a || src[3] !== 0x58 || src[4] !== 0x5a || src[5] !== 0) throw new Error('not an XZ stream');
  const checkSize = [0, 4, 4, 4, 8, 8, 8, 16, 16, 16, 32, 32, 32, 64, 64, 64][src[7] & 15]; let p = 12; const parts = []; let total = 0;
  for (let g = 0; g < 1e5 && p < src.length && src[p] !== 0; g++) {            // a zero byte here starts the index
    const start = p, hsize = (src[p] + 1) * 4, flags = src[p + 1]; if ((flags & 3) !== 0) throw new Error('XZ blocks with filter chains other than plain LZMA2 are not supported');
    let q = p + 2; const varint = () => { let v = 0, s = 0; for (;;) { const c = src[q++]; v += (c & 0x7f) * 2 ** s; if (!(c & 0x80)) return v; s += 7; if (s > 63) throw new Error('bad XZ integer'); } };
    if (flags & 0x40) varint(); if (flags & 0x80) varint();
    if (varint() !== 0x21) throw new Error('XZ filter is not LZMA2');
    p = start + hsize; const r = lzma2(src, p, cap - total); parts.push(r.bytes); total += r.bytes.length; p = r.end;
    while ((p - start) % 4) p++; p += checkSize;
  }
  if (p >= src.length) throw new Error('XZ stream is truncated (no index after the last block)');
  const out = new Uint8Array(total); let o = 0; for (const a of parts) { out.set(a, o); o += a.length; }
  return out;
}
