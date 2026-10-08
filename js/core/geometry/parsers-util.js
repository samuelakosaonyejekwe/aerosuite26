// Shared low-level helpers for the geometry/mesh readers: bounded byte scanning, a tiny XML tree parser,
// a ZIP central-directory reader, the element-type tables and the boundary-face extraction used for display.
// Everything here treats its input as untrusted: loops are bounded by the buffer length and counts are
// validated before anything is allocated.

export const LIMITS = { verts: 40e6, elems: 60e6, inflate: 600e6 };

/** Validate a count read from a file before it is used for a loop or an allocation. */
export function guard(n, what, max = LIMITS.elems) {
  if (!Number.isFinite(n) || n < 0 || Math.floor(n) !== n) throw new Error(`invalid ${what} count (${n})`);
  if (n > max) throw new Error(`${what} count ${n} exceeds the import limit of ${max}`);
  return n;
}

const UTF8 = new TextDecoder('utf-8');
/** Decode bytes[s,e) to a string (fast path for short ASCII runs). */
export function str(b, s, e) {
  if (e > b.length) e = b.length;
  if (e <= s) return '';
  if (e - s > 64) return UTF8.decode(b.subarray(s, e));
  let o = '';
  for (let i = s; i < e; i++) { const c = b[i]; if (c > 127) return UTF8.decode(b.subarray(s, e)); o += String.fromCharCode(c); }
  return o;
}
/** Whole-buffer text (UTF-8, BOM stripped). */
export function text(b) { return UTF8.decode(b[0] === 0xef && b[1] === 0xbb && b[2] === 0xbf ? b.subarray(3) : b); }
/** Number with Fortran 'D' exponents tolerated. */
export function toNum(t) { const v = +t; return v === v || !t ? v : +t.replace(/[dD]/, 'e'); }
/** All numbers in a whitespace/comma separated string. */
export function nums(s) { const o = []; for (const t of s.split(/[\s,;]+/)) if (t) o.push(toNum(t)); return o; }
export const view = (b) => new DataView(b.buffer, b.byteOffset, b.byteLength);
/** True when the head of the buffer looks like text rather than binary data. */
export function looksText(b, n = 4096) {
  const m = Math.min(b.length, n); let bad = 0;
  for (let i = 0; i < m; i++) { const c = b[i]; if (c === 0 || (c < 9) || (c > 13 && c < 32 && c !== 27)) bad++; }
  return m > 0 && bad <= m * 0.01;
}

/** Bounded token/line scanner over a byte buffer (no giant strings, no regex backtracking). */
export class Scanner {
  constructor(b, { pos = 0, comment = '', sep = '' } = {}) {
    this.b = b; this.p = pos; this.n = b.length; this.cc = comment ? comment.charCodeAt(0) : -1;
    this.ws = new Uint8Array(256);
    for (const c of [9, 10, 11, 12, 13, 32]) this.ws[c] = 1;
    for (const ch of sep) this.ws[ch.charCodeAt(0)] = 1;
  }
  eof() { return this.p >= this.n; }
  /** Rest of the current line (without the line break); null at end of input. */
  line() {
    if (this.p >= this.n) return null;
    const b = this.b; let e = this.p;
    while (e < this.n && b[e] !== 10) e++;
    const s = str(b, this.p, e > this.p && b[e - 1] === 13 ? e - 1 : e);
    this.p = e + 1;
    return s;
  }
  /** Next whitespace-delimited token (comments skipped); null at end of input. */
  token() {
    const { b, n, ws, cc } = this; let p = this.p;
    for (;;) {
      while (p < n && ws[b[p]]) p++;
      if (p < n && b[p] === cc) { while (p < n && b[p] !== 10) p++; continue; }
      break;
    }
    if (p >= n) { this.p = p; return null; }
    const s = p;
    while (p < n && !ws[b[p]]) p++;
    this.p = p;
    return str(b, s, p);
  }
  peek() { const p = this.p, t = this.token(); this.p = p; return t; }
  num() { const t = this.token(); return t === null ? NaN : toNum(t); }
  /** Integer token; throws a clear error when the input ends or is not a number. */
  int(what = 'integer') { const t = this.token(), v = t === null ? NaN : parseInt(t, 10); if (v !== v) throw new Error(`unexpected ${t === null ? 'end of file' : `token "${t.slice(0, 20)}"`} while reading ${what}`); return v; }
  /** Float token with the same end-of-input check. */
  float(what = 'number') { const t = this.token(), v = t === null ? NaN : toNum(t); if (v !== v) throw new Error(`unexpected ${t === null ? 'end of file' : `token "${t.slice(0, 20)}"`} while reading ${what}`); return v; }
}

// ---------- element tables ----------
// n = nodes per element, c = corner nodes (always listed first), dim = topological dimension.
export const ET = {
  point1: { n: 1, c: 1, dim: 0 }, line2: { n: 2, c: 2, dim: 1 }, line3: { n: 3, c: 2, dim: 1 },
  tri3: { n: 3, c: 3, dim: 2 }, tri6: { n: 6, c: 3, dim: 2 }, quad4: { n: 4, c: 4, dim: 2 }, quad8: { n: 8, c: 4, dim: 2 }, quad9: { n: 9, c: 4, dim: 2 },
  tet4: { n: 4, c: 4, dim: 3 }, tet10: { n: 10, c: 4, dim: 3 }, hex8: { n: 8, c: 8, dim: 3 }, hex20: { n: 20, c: 8, dim: 3 }, hex27: { n: 27, c: 8, dim: 3 },
  wedge6: { n: 6, c: 6, dim: 3 }, wedge15: { n: 15, c: 6, dim: 3 }, wedge18: { n: 18, c: 6, dim: 3 }, pyr5: { n: 5, c: 5, dim: 3 }, pyr13: { n: 13, c: 5, dim: 3 }, pyr14: { n: 14, c: 5, dim: 3 },
};
export const shapeOf = (type) => type.replace(/\d+$/, '');
// Corner-node faces, wound so the normal points out of a positively oriented element
// (hex: 0-3 bottom ring, 4-7 top ring; wedge: 0-2 bottom, 3-5 top; pyramid: 0-3 base, 4 apex).
export const FACES = {
  tet: [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
  hex: [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]],
  wedge: [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]],
  pyr: [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
};
export const EDGES = {
  tri: [[0, 1], [1, 2], [2, 0]], quad: [[0, 1], [1, 2], [2, 3], [3, 0]],
  tet: [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]],
  hex: [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]],
  wedge: [[0, 1], [1, 2], [2, 0], [3, 4], [4, 5], [5, 3], [0, 3], [1, 4], [2, 5]],
  pyr: [[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [1, 4], [2, 4], [3, 4]],
};

/** Accumulates nodes (optionally with file ids), element blocks and groups while a reader runs. */
export class Mesh {
  constructor() { this.xyz = []; this.ids = null; this.blocks = new Map(); this.groups = []; this.gkey = new Map(); this.badRefs = 0; }
  get nNodes() { return this.xyz.length / 3; }
  node(x, y, z, id) { const i = this.xyz.length / 3; this.xyz.push(x, y, z); if (id !== undefined) (this.ids ??= new Map()).set(id, i); return i; }
  /** File node id → 0-based index (-1 when the id is unknown). */
  idx(id) { const i = this.ids ? this.ids.get(id) : id; return i === undefined ? -1 : i; }
  /** Group index for (kind, tag); created on first use. Group `id` always equals its index in model.groups. */
  group(kind, tag, name) {
    const k = kind + '\u0001' + tag; let g = this.gkey.get(k);
    if (g === undefined) { g = this.groups.length; this.groups.push({ id: g, name: name ?? String(tag), kind, count: 0, tag }); this.gkey.set(k, g); }
    return g;
  }
  /** Add one element given 0-based node indices; elements with unknown nodes are counted and dropped. */
  elem(type, nodes, g = -1) {
    const n = ET[type].n;
    for (let i = 0; i < n; i++) { const v = nodes[i]; if (!(v >= 0)) { this.badRefs++; return; } }
    let blk = this.blocks.get(type);
    if (!blk) this.blocks.set(type, (blk = { conn: [], group: [], any: false }));
    for (let i = 0; i < n; i++) blk.conn.push(nodes[i]);
    blk.group.push(g);
    if (g >= 0) { blk.any = true; this.groups[g].count++; }
  }
  /** Add one element given file node ids. */
  elemIds(type, ids, g = -1) { const n = ET[type].n, a = new Array(n); for (let i = 0; i < n; i++) a[i] = this.idx(ids[i]); this.elem(type, a, g); }
  /** Queue an element whose nodes may be defined later in the file; resolved in result(). */
  later(type, ids, g = -1) { (this.pending ??= []).push(type, ids, g); }
  result(extra = {}) {
    if (this.pending) { const q = this.pending; this.pending = null; for (let i = 0; i < q.length; i += 3) this.elemIds(q[i], q[i + 1], q[i + 2]); }
    const elements = [];
    for (const [type, blk] of this.blocks) {
      const e = { type, nodesPer: ET[type].n, count: blk.group.length, conn: Uint32Array.from(blk.conn) };
      if (blk.any) e.group = Int32Array.from(blk.group);
      elements.push(e);
    }
    const out = { positions: Float64Array.from(this.xyz), elements, groups: this.groups, ...extra };
    if (this.badRefs) (out.warnings ??= []).push(`${this.badRefs} element(s) referenced undefined nodes and were dropped.`);
    return out;
  }
}

/** Append the cells of a structured ni×nj×nk block (node index = base + i + ni·(j + nj·k)) as hex8 / quad4 / line2. */
export function structuredBlock(M, ni, nj, nk, base, g) {
  const dims = [ni, nj, nk].map((d, a) => ({ d, a })).filter((o) => o.d > 1);
  const stride = [1, ni, ni * nj], id = (i, j, k) => base + i + ni * (j + nj * k);
  if (dims.length === 3) {
    for (let k = 0; k < nk - 1; k++) for (let j = 0; j < nj - 1; j++) for (let i = 0; i < ni - 1; i++) {
      const a = id(i, j, k);
      M.elem('hex8', [a, a + 1, a + 1 + ni, a + ni, a + stride[2], a + 1 + stride[2], a + 1 + ni + stride[2], a + ni + stride[2]], g);
    }
  } else if (dims.length === 2) {
    const [p, q] = dims, sp = stride[p.a], sq = stride[q.a];
    for (let v = 0; v < q.d - 1; v++) for (let u = 0; u < p.d - 1; u++) { const a = base + u * sp + v * sq; M.elem('quad4', [a, a + sp, a + sp + sq, a + sq], g); }
  } else if (dims.length === 1) {
    const s = stride[dims[0].a];
    for (let u = 0; u < dims[0].d - 1; u++) M.elem('line2', [base + u * s, base + (u + 1) * s], g);
  }
}

/**
 * Display triangles for a set of element blocks: surface elements are triangulated directly; for volume
 * elements only the boundary faces (faces referenced by exactly one element) are kept. Surface elements that
 * coincide with a volume boundary face are not drawn twice.
 */
export function displayTriangles(elements, nVerts) {
  const big = nVerts > 200000, N = nVerts;
  const key = (a, b, c) => {
    if (a > b) { const t = a; a = b; b = t; } if (b > c) { const t = b; b = c; c = t; } if (a > b) { const t = a; a = b; b = t; }
    return big ? a + '_' + b + '_' + c : (a * N + b) * N + c;
  };
  const key4 = (f) => { const s = [f[0], f[1], f[2], f[3]].sort((x, y) => x - y); return key(s[0], s[1], s[2]); };
  const faces = new Map();
  let hasVol = false;
  for (const el of elements) {
    const info = ET[el.type]; if (!info || info.dim !== 3) continue;
    hasVol = true;
    const F = FACES[shapeOf(el.type)], np = el.nodesPer, c = el.conn;
    for (let e = 0; e < el.count; e++) {
      const o = e * np;
      for (const f of F) {
        const n = f.map((i) => c[o + i]), k = n.length === 3 ? key(n[0], n[1], n[2]) : key4(n);
        if (faces.has(k)) faces.set(k, null); else faces.set(k, n);
      }
    }
  }
  const tri = [];
  const push = (n) => { tri.push(n[0], n[1], n[2]); if (n.length === 4) tri.push(n[0], n[2], n[3]); };
  for (const n of faces.values()) if (n) push(n);
  for (const el of elements) {
    const info = ET[el.type]; if (!info || info.dim !== 2) continue;
    const np = el.nodesPer, c = el.conn, nc = info.c;
    for (let e = 0; e < el.count; e++) {
      const o = e * np, n = nc === 3 ? [c[o], c[o + 1], c[o + 2]] : [c[o], c[o + 1], c[o + 2], c[o + 3]];
      if (hasVol && faces.get(nc === 3 ? key(n[0], n[1], n[2]) : key4(n))) continue;
      push(n);
    }
  }
  return Uint32Array.from(tri);
}

/** Line segments (pairs) from 1-D element blocks. */
export function displayLines(elements) {
  const out = [];
  for (const el of elements) {
    if (ET[el.type]?.dim !== 1) continue;
    for (let e = 0; e < el.count; e++) out.push(el.conn[e * el.nodesPer], el.conn[e * el.nodesPer + 1]);
  }
  return Uint32Array.from(out);
}

// ---------- XML (small, forgiving tree parser; namespaces are stripped) ----------
const ENT = { amp: '&', lt: '<', gt: '>', quot: '"', apos: "'" };
const unent = (s) => (s.includes('&') ? s.replace(/&(#x[0-9a-f]+|#\d+|\w+);/gi, (m, e) => (e[0] === '#' ? String.fromCodePoint(e[1] === 'x' || e[1] === 'X' ? parseInt(e.slice(2), 16) : +e.slice(1)) : ENT[e] ?? m)) : s);
export function xmlParse(src, maxNodes = 30e6) {
  const root = { name: '#root', attrs: {}, children: [], text: '' }, stack = [root], n = src.length;
  let i = 0, count = 0;
  while (i < n) {
    const lt = src.indexOf('<', i);
    const top = stack[stack.length - 1];
    if (lt < 0) break;
    if (lt > i) { const t = src.slice(i, lt); if (/\S/.test(t)) top.text += t; }
    if (src.startsWith('<!--', lt)) { const e = src.indexOf('-->', lt); i = e < 0 ? n : e + 3; continue; }
    if (src.startsWith('<![CDATA[', lt)) { const e = src.indexOf(']]>', lt); top.text += src.slice(lt + 9, e < 0 ? n : e); i = e < 0 ? n : e + 3; continue; }
    const gt = src.indexOf('>', lt);
    if (gt < 0) break;
    const c1 = src[lt + 1];
    if (c1 === '?' || c1 === '!') { i = gt + 1; continue; }
    if (c1 === '/') { if (stack.length > 1) stack.pop(); i = gt + 1; continue; }
    const self = src[gt - 1] === '/', tag = src.slice(lt + 1, self ? gt - 1 : gt), sp = tag.search(/\s/);
    let name = sp < 0 ? tag : tag.slice(0, sp);
    const colon = name.indexOf(':'); if (colon >= 0) name = name.slice(colon + 1);
    const node = { name, attrs: {}, children: [], text: '' };
    if (sp >= 0) { const re = /([\w:.-]+)\s*=\s*(?:"([^"]*)"|'([^']*)')/g; let m; while ((m = re.exec(tag))) node.attrs[m[1].includes(':') ? m[1].slice(m[1].indexOf(':') + 1) : m[1]] = unent(m[2] ?? m[3]); }
    if (++count > maxNodes) throw new Error('XML document has too many elements');
    top.children.push(node);
    if (!self) { if (stack.length > 2000) throw new Error('XML nesting too deep'); stack.push(node); }
    i = gt + 1;
  }
  return root;
}
/** All descendants with the given tag name (iterative, document order). */
export function xmlAll(node, name, out = []) {
  const st = [node];
  while (st.length) { const nd = st.pop(); if (nd.name === name && nd !== node) out.push(nd); for (let i = nd.children.length - 1; i >= 0; i--) st.push(nd.children[i]); }
  return out;
}
export const xmlFirst = (node, name) => xmlAll(node, name)[0] ?? null;
export const xmlChild = (node, name) => node.children.find((c) => c.name === name) ?? null;

// ---------- ZIP (stored + deflate entries via the central directory) ----------
async function inflateRaw(data, cap) {
  const ds = new DecompressionStream('deflate-raw'), w = ds.writable.getWriter();
  w.write(data).catch(() => {}); w.close().catch(() => {});
  const r = ds.readable.getReader(), chunks = []; let total = 0;
  for (;;) {
    const { done, value } = await r.read();
    if (done) break;
    total += value.length;
    if (total > cap) { r.cancel().catch(() => {}); throw new Error('compressed entry expands beyond the import limit'); }
    chunks.push(value);
  }
  const out = new Uint8Array(total); let o = 0;
  for (const c of chunks) { out.set(c, o); o += c.length; }
  return out;
}
/** List the entries of a ZIP archive: [{ name, method, size, read() }]. Throws on anything that is not a plain ZIP. */
export function zipEntries(b) {
  const dv = view(b), n = b.length; let e = -1;
  for (let i = n - 22; i >= Math.max(0, n - 66000); i--) if (dv.getUint32(i, true) === 0x06054b50) { e = i; break; }
  if (e < 0) throw new Error('ZIP end-of-central-directory record not found (truncated archive?)');
  const count = dv.getUint16(e + 10, true); let p = dv.getUint32(e + 16, true);
  if (p === 0xffffffff) throw new Error('ZIP64 archives are not supported');
  const out = [];
  for (let k = 0; k < count && p + 46 <= n; k++) {
    if (dv.getUint32(p, true) !== 0x02014b50) break;
    const method = dv.getUint16(p + 10, true), csize = dv.getUint32(p + 20, true), size = dv.getUint32(p + 24, true);
    const nl = dv.getUint16(p + 28, true), xl = dv.getUint16(p + 30, true), cl = dv.getUint16(p + 32, true), lho = dv.getUint32(p + 42, true);
    const name = str(b, p + 46, p + 46 + nl);
    out.push({
      name, method, size,
      async read() {
        if (lho + 30 > n || dv.getUint32(lho, true) !== 0x04034b50) throw new Error(`ZIP entry "${name}" has a bad local header`);
        const start = lho + 30 + dv.getUint16(lho + 26, true) + dv.getUint16(lho + 28, true);
        if (start + csize > n) throw new Error(`ZIP entry "${name}" is truncated`);
        const data = b.subarray(start, start + csize);
        if (method === 0) return data;
        if (method !== 8) throw new Error(`ZIP entry "${name}" uses unsupported compression method ${method}`);
        if (typeof DecompressionStream === 'undefined') throw new Error('this runtime has no DecompressionStream; cannot inflate ZIP content');
        return inflateRaw(data, LIMITS.inflate);
      },
    });
    p += 46 + nl + xl + cl;
  }
  return out;
}

/** base64 → bytes (atob exists in browsers, workers and Node). */
export function base64(s) {
  const bin = atob(s.replace(/\s+/g, '')), out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

// ---------- 4×4 transforms (column-major, as in glTF) ----------
export const M4 = {
  I: () => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
  mul(a, b) { const o = new Array(16); for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) { let s = 0; for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k]; o[c * 4 + r] = s; } return o; },
  trs(t = [0, 0, 0], q = [0, 0, 0, 1], s = [1, 1, 1]) {
    const [x, y, z, w] = q;
    return [(1 - 2 * (y * y + z * z)) * s[0], 2 * (x * y + z * w) * s[0], 2 * (x * z - y * w) * s[0], 0,
      2 * (x * y - z * w) * s[1], (1 - 2 * (x * x + z * z)) * s[1], 2 * (y * z + x * w) * s[1], 0,
      2 * (x * z + y * w) * s[2], 2 * (y * z - x * w) * s[2], (1 - 2 * (x * x + y * y)) * s[2], 0, t[0], t[1], t[2], 1];
  },
  det3: (m) => m[0] * (m[5] * m[10] - m[6] * m[9]) - m[4] * (m[1] * m[10] - m[2] * m[9]) + m[8] * (m[1] * m[6] - m[2] * m[5]),
  apply: (m, x, y, z) => [m[0] * x + m[4] * y + m[8] * z + m[12], m[1] * x + m[5] * y + m[9] * z + m[13], m[2] * x + m[6] * y + m[10] * z + m[14]],
};

export const UNIT_WORDS = { m: 'm', meter: 'm', metre: 'm', meters: 'm', metres: 'm', mm: 'mm', millimeter: 'mm', millimetre: 'mm', millimeters: 'mm', cm: 'cm', centimeter: 'cm', centimetre: 'cm', in: 'in', inch: 'in', inches: 'in', ft: 'ft', foot: 'ft', feet: 'ft' };
/** Map a unit word from a file to the platform's unit ids; unknown words give null (never guessed). */
export const unitFromWord = (w) => (w ? UNIT_WORDS[String(w).trim().toLowerCase()] ?? null : null);
