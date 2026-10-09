// bzip2 decompression (single stream), for the .grib2.bz2 files of the DWD open-data server. Node's zlib has no
// bzip2, and the snapshot job has no dependencies. Block check sums are not verified; a damaged stream is caught
// by the structural checks here and by the GRIB reader's own length checks.

export function bunzip2(buf) {
  let pos = 0, cur = 0, left = 0;
  const bits = (n) => { let v = 0; while (n > 0) { if (left === 0) { if (pos >= buf.length) throw new Error('bzip2: unexpected end of data'); cur = buf[pos++]; left = 8; } const take = Math.min(n, left); v = v * (1 << take) + ((cur >> (left - take)) & ((1 << take) - 1)); left -= take; n -= take; } return v; };
  if (bits(8) !== 0x42 || bits(8) !== 0x5a || bits(8) !== 0x68) throw new Error('not a bzip2 stream');
  const level = bits(8) - 0x30; if (level < 1 || level > 9) throw new Error('bzip2: bad block size');
  const blockSize = level * 100000, tt = new Uint32Array(blockSize), chunks = [];
  for (;;) {
    const m1 = bits(24), m2 = bits(24);
    if (m1 === 0x177245 && m2 === 0x385090) break;
    if (m1 !== 0x314159 || m2 !== 0x265359) throw new Error('bzip2: bad block header');
    bits(32); if (bits(1)) throw new Error('bzip2: randomised blocks are not supported');
    const origPtr = bits(24), used = bits(16), symToByte = [];
    for (let i = 0; i < 16; i++) if (used & (0x8000 >> i)) { const b = bits(16); for (let j = 0; j < 16; j++) if (b & (0x8000 >> j)) symToByte.push(i * 16 + j); }
    const symTotal = symToByte.length + 2, groups = bits(3), nSel = bits(15); if (groups < 2 || groups > 6 || !nSel) throw new Error('bzip2: bad table count');
    const order = [0, 1, 2, 3, 4, 5].slice(0, groups), selectors = new Uint8Array(nSel);
    for (let i = 0; i < nSel; i++) { let j = 0; while (bits(1)) if (++j >= groups) throw new Error('bzip2: bad selector'); const t = order[j]; order.splice(j, 1); order.unshift(t); selectors[i] = t; }
    const tables = [];
    for (let g = 0; g < groups; g++) {
      const len = new Uint8Array(symTotal); let t = bits(5), minLen = 32, maxLen = 0;
      for (let i = 0; i < symTotal; i++) { for (;;) { if (t < 1 || t > 20) throw new Error('bzip2: bad code length'); if (!bits(1)) break; t += bits(1) ? -1 : 1; } len[i] = t; if (t < minLen) minLen = t; if (t > maxLen) maxLen = t; }
      const perm = [], count = new Int32Array(maxLen + 2), limit = new Int32Array(maxLen + 2), base = new Int32Array(maxLen + 2);
      for (let l = minLen; l <= maxLen; l++) for (let i = 0; i < symTotal; i++) if (len[i] === l) perm.push(i);
      for (let i = 0; i < symTotal; i++) count[len[i]]++;
      for (let l = minLen, code = 0, idx = 0; l <= maxLen; l++) { base[l] = idx - code; code += count[l]; idx += count[l]; limit[l] = code - 1; code <<= 1; }
      tables.push({ minLen, maxLen, limit, base, perm });
    }
    const byteCount = new Int32Array(256), mtf = symToByte.slice();
    let n = 0, runPos = 0, runLen = 0, sel = 0, symLeft = 0, tab = null;
    for (;;) {
      if (symLeft-- === 0) { symLeft = 49; if (sel >= nSel) throw new Error('bzip2: ran out of selectors'); tab = tables[selectors[sel++]]; }
      let l = tab.minLen, code = bits(l);
      while (code > tab.limit[l]) { if (++l > tab.maxLen) throw new Error('bzip2: bad code'); code = code * 2 + bits(1); }
      const sym = tab.perm[code + tab.base[l]];
      if (sym <= 1) { if (!runPos) { runPos = 1; runLen = 0; } runLen += runPos << sym; runPos <<= 1; continue; }
      if (runPos) { runPos = 0; if (n + runLen > blockSize) throw new Error('bzip2: block overflow'); const b = mtf[0]; byteCount[b] += runLen; for (; runLen > 0; runLen--) tt[n++] = b; }
      if (sym === symTotal - 1) break;
      if (n >= blockSize) throw new Error('bzip2: block overflow');
      const b = mtf[sym - 1]; mtf.splice(sym - 1, 1); mtf.unshift(b); byteCount[b]++; tt[n++] = b;
    }
    if (origPtr >= n && n) throw new Error('bzip2: bad block pointer');
    for (let i = 0, sum = 0; i < 256; i++) { const c = byteCount[i]; byteCount[i] = sum; sum += c; }
    for (let i = 0; i < n; i++) { const b = tt[i] & 0xff; tt[byteCount[b]++] |= i << 8; }
    let out = Buffer.allocUnsafe(n + 1024), o = 0, p = n ? tt[origPtr] : 0, current = p & 0xff, run = -1; p >>>= 8;
    for (let k = 0; k < n; k++) {
      const previous = current; p = tt[p]; current = p & 0xff; p >>>= 8;
      let copies = 1, byte = current;
      if (run++ === 3) { copies = current; byte = previous; current = -1; }
      if (o + copies > out.length) { const bigger = Buffer.allocUnsafe(Math.max(out.length * 2, o + copies + 1024)); out.copy(bigger, 0, 0, o); out = bigger; }
      for (; copies > 0; copies--) out[o++] = byte;
      if (current !== previous) run = 0;
    }
    chunks.push(out.subarray(0, o));
  }
  return Buffer.concat(chunks);
}
