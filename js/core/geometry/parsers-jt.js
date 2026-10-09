// JT (ISO 14306) reader - container and logical scene graph, written from the published "JT File Format Reference".
// What is decoded here: the file header, the table of contents (32-bit offsets, or 64-bit in JT 10), segment
// headers, ZLIB- and (JT 10) LZMA2/XZ-compressed segments, the LSG in its JT 8 / JT 9 / JT 10 layouts (assembly /
// instance / part / LOD / shape nodes, geometric transforms, property atoms and the property table: names, units,
// late-loaded segment references) and XT B-Rep segments, whose embedded Parasolid stream is handed to the
// Parasolid translator. The compressed tessellation (Shape LOD segments) is decoded in parsers-jtshape.js.

import { str, view, inflate } from './parsers-util.js';
import { unxz } from './parsers-lzma.js';

const guidAt = (b, p, le) => { const dv = view(b), h = (v, n) => v.toString(16).padStart(n, '0'); let s = h(dv.getUint32(p, le), 8) + '-' + h(dv.getUint16(p + 4, le), 4) + '-' + h(dv.getUint16(p + 6, le), 4) + '-'; for (let i = 8; i < 16; i++) s += h(b[p + i], 2); return s; };
const EOE = 'ffffffff-ffff-ffff-ffffffffffffffff';
const T = {
  '10dd1035-2ac8-11d1-9b6b0080c7bb5997': 'base', '10dd103e-2ac8-11d1-9b6b0080c7bb5997': 'partition', '10dd101b-2ac8-11d1-9b6b0080c7bb5997': 'group', '10dd102a-2ac8-11d1-9b6b0080c7bb5997': 'instance',
  'ce357244-38fb-11d1-a506006097bdc6e1': 'part', 'ce357245-38fb-11d1-a506006097bdc6e1': 'metadata', '10dd102c-2ac8-11d1-9b6b0080c7bb5997': 'lod', '10dd104c-2ac8-11d1-9b6b0080c7bb5997': 'rangeLod', '10dd10f3-2ac8-11d1-9b6b0080c7bb5997': 'switch',
  '10dd1077-2ac8-11d1-9b6b0080c7bb5997': 'triStripShape', '10dd1046-2ac8-11d1-9b6b0080c7bb5997': 'polylineShape', 'e40373c1-1ad9-11d3-9daf00a0c9c7dd42': 'primitiveShape',
  '10dd1083-2ac8-11d1-9b6b0080c7bb5997': 'transform', '10dd1030-2ac8-11d1-9b6b0080c7bb5997': 'material',
  '10dd106e-2ac8-11d1-9b6b0080c7bb5997': 'stringProp', '10dd102b-2ac8-11d1-9b6b0080c7bb5997': 'intProp', '10dd1019-2ac8-11d1-9b6b0080c7bb5997': 'floatProp', '10dd1004-2ac8-11d1-9b6b0080c7bb5997': 'refProp', 'ce357246-38fb-11d1-a506006097bdc6e1': 'dateProp', 'e0b05be5-fbbd-11d1-a3a700aa00d10954': 'lateLoadedProp',
};
const SEGMENT = { 1: 'Logical Scene Graph', 2: 'JT B-Rep', 3: 'PMI Data', 4: 'Meta Data', 6: 'Shape', 7: 'Shape LOD0', 8: 'Shape LOD1', 9: 'Shape LOD2', 10: 'Shape LOD3', 11: 'Shape LOD4', 12: 'Shape LOD5', 13: 'Shape LOD6', 14: 'Shape LOD7', 15: 'Shape LOD8', 16: 'Shape LOD9', 17: 'XT B-Rep', 18: 'Wireframe', 20: 'ULP', 24: 'LWPA' };
const ZLIB_TYPES = new Set([1, 2, 3, 4, 17, 18, 20, 24]);
const I4 = () => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
const mul4 = (a, b) => { const o = new Array(16); for (let r = 0; r < 4; r++) for (let c = 0; c < 4; c++) { let s = 0; for (let k = 0; k < 4; k++) s += a[r * 4 + k] * b[k * 4 + c]; o[r * 4 + c] = s; } return o; };

/**
 * Parse the LSG element stream with one of the node-data layouts: 'v9' has I16 version fields (JT 9.x),
 * 'v8' has none (JT 8.x), 'v10' uses U8 version fields. Returns null when the stream does not fit the layout.
 */
function parseLSG(d, le, layout) {
  const dv = view(d), n = d.length; let p = 0;
  const nodes = new Map(), attrs = new Map(), props = new Map(); let order = [];
  const ver = (q) => (layout === 'v9' ? q + 2 : layout === 'v10' ? q + 1 : q);
  const mb = (q) => { const k = dv.getInt32(q, le); if (!(k >= 0 && q + 4 + 2 * k <= n)) throw new Error('bad string'); let s = ''; for (let i = 0; i < k; i++) s += String.fromCharCode(dv.getUint16(q + 4 + 2 * i, le)); return { s, end: q + 4 + 2 * k }; };
  const section = (handler) => {
    for (let g = 0; g < 5e6; g++) {
      if (p + 20 > n) throw new Error('element stream ends early');
      const len = dv.getInt32(p, le), id = guidAt(d, p + 4, le), end = p + 4 + len; if (!(len >= 16 && end <= n)) throw new Error('bad element length');
      if (id === EOE) { p = end; return; }
      const baseType = d[p + 20], objId = dv.getInt32(p + 21, le); handler(T[id] || id, baseType, objId, p + 25, end); p = end;
    }
  };
  try {
    section((kind, baseType, id, q, end) => {
      if (baseType === 3) {                                  // attribute: base attribute data, then type-specific data
        let r = layout === 'v10' ? q + 9 + (d[q] >= 2 ? 1 : 0) : ver(q) + 5;   // JT 10.5 (base attribute version 2) adds one byte         // version, state flags, field-inhibit flags (U8 + U32 + U32 in JT 10; 5 bytes after the version before)
        const a = { kind, id }; attrs.set(id, a);
        if (kind === 'transform') {
          r = ver(r); const mask = dv.getUint16(r, le); r += 2; let nvals = 0; for (let i = 0; i < 16; i++) if (mask & (1 << i)) nvals++;
          // the stored elements are F32 in JT 8 and F64 in the JT 9 files examined: take the width from the element length
          const room = end - r, width = !nvals ? 4 : room === 8 * nvals || room === 8 * nvals + 4 ? 8 : room === 4 * nvals || room === 4 * nvals + 4 ? 4 : 0; if (!width) throw new Error('transform element has an unexpected size');   // JT 10.5 appends a 4-byte field
          const m = I4(); for (let i = 0; i < 16; i++) if (mask & (0x8000 >> i)) { m[i] = width === 8 ? dv.getFloat64(r, le) : dv.getFloat32(r, le); r += width; } a.m = m;
        }
        return;
      }
      if (baseType > 2 && baseType !== 255) return;
      let r = ver(q); r += 4;                                // node flags
      const na = dv.getInt32(r, le); r += 4; if (!(na >= 0 && r + 4 * na <= end)) throw new Error('bad attribute count');
      const node = { kind, id, attrs: [], children: [] }; for (let i = 0; i < na; i++) { node.attrs.push(dv.getInt32(r, le)); r += 4; }
      if (baseType === 1) { r = ver(r); const nc = dv.getInt32(r, le); r += 4; if (!(nc >= 0 && r + 4 * nc <= end)) throw new Error('bad child count'); for (let i = 0; i < nc; i++) { node.children.push(dv.getInt32(r, le)); r += 4; } if (kind === 'partition') { try { const q2 = (layout === 'v10' ? r + 1 : r) + 4, nm = mb(q2); if (nm.s && nm.end <= end && /^[\x20-\uffff]+$/.test(nm.s)) node.fileName = nm.s; } catch { /* no file name */ } } }
      else if (kind === 'instance') { r = ver(r); node.children.push(dv.getInt32(r, le)); }
      nodes.set(id, node); order.push(id);
    });
    section((kind, baseType, id, q, end) => {
      let r = ver(q) + 4;                                    // base property atom data: version, state flags
      if (kind === 'stringProp') { r = ver(r); props.set(id, mb(r).s); }
      else if (kind === 'intProp') { r = ver(r); props.set(id, dv.getInt32(r, le)); }
      else if (kind === 'floatProp') { r = ver(r); props.set(id, dv.getFloat32(r, le)); }
      else if (kind === 'lateLoadedProp') { r = ver(r); if (r + 20 > end) throw new Error('late-loaded atom overruns'); props.set(id, { segment: guidAt(d, r, le), segmentType: dv.getInt32(r + 16, le) }); }
      else props.set(id, null);
    });
    // property table
    p += 2; const count = dv.getInt32(p, le); p += 4; const table = new Map();     // the table carries a version number in every layout
    if (!(count >= 0 && count < 5e6)) throw new Error('bad property table');
    for (let i = 0; i < count; i++) { const el = dv.getInt32(p, le); p += 4; const pairs = []; for (let g = 0; g < 1e5; g++) { const k = dv.getInt32(p, le); p += 4; if (k === 0) break; const v = dv.getInt32(p, le); p += 4; pairs.push([props.get(k), props.get(v)]); } table.set(el, pairs); }
    // plausibility: every child reference resolves
    let refs = 0, hit = 0; for (const nd of nodes.values()) for (const c of nd.children) { refs++; if (nodes.has(c)) hit++; }
    if (!nodes.size || hit < refs) return null;
    return { nodes, attrs, props, table, root: order[0], layout };
  } catch (e) { if (typeof process === 'object' && process.env && process.env.JT_DEBUG) console.error('LSG layout', layout, 'rejected:', e.message); return null; }
}

/** Read the JT container and scene graph. Returns { meta, lsg, segments, xtStreams: [{ segment, bytes }] }. */
export async function parseJT(b) {
  const head = str(b, 0, 80), vm = /^Version\s+(\d+)\.(\d+)/.exec(head); if (!vm) throw new Error('not a JT file (no version header)');
  if (b.length < 125) throw Object.assign(new Error('the JT file ends inside its header (truncated)'), { meta: { version: `${vm[1]}.${vm[2]}` } });
  const major = +vm[1], minor = +vm[2], le = b[80] === 0, dv = view(b);
  const meta = { version: `${major}.${minor}`, byteOrder: le ? 'little-endian' : 'big-endian', writer: head.replace(/^Version\s+[\d.]+\s+JT/, '').trim() || null, segments: [], lsgLayout: null, nodeCensus: {}, partNames: [], units: null, shapeSegments: 0, xtSegments: 0, jtBrepSegments: 0 };
  // JT 10 widened the TOC offsets to 64 bits
  const wide = major >= 10, tocAt = wide ? Number(dv.getBigUint64(85, le)) : dv.getInt32(85, le), lsgId = guidAt(b, wide ? 93 : 89, le);
  if (!(tocAt > 0 && tocAt + 4 <= b.length)) throw Object.assign(new Error('the table of contents lies outside the file (truncated)'), { meta });
  const count = dv.getInt32(tocAt, le), esz = wide ? 32 : 28; if (!(count >= 0 && tocAt + 4 + count * esz <= b.length)) throw Object.assign(new Error('the table of contents is truncated'), { meta });
  const segments = [];
  for (let i = 0; i < count; i++) { const q = tocAt + 4 + i * esz, id = guidAt(b, q, le), offset = wide ? Number(dv.getBigUint64(q + 16, le)) : dv.getInt32(q + 16, le), length = dv.getInt32(q + (wide ? 24 : 20), le), type = dv.getUint32(q + (wide ? 28 : 24), le) >>> 24; segments.push({ id, offset, length, type }); }
  const census = {}; for (const s of segments) { const k = SEGMENT[s.type] ?? `type ${s.type}`; census[k] = (census[k] || 0) + 1; if (s.type >= 6 && s.type <= 16) meta.shapeSegments++; if (s.type === 17) meta.xtSegments++; if (s.type === 2) meta.jtBrepSegments++; }
  meta.segments = census;
  /** Element bytes of a segment (inflated when the segment says so). */
  const payload = async (s) => {
    const q = s.offset + 24; if (!(s.offset >= 0 && q + 9 <= b.length) || s.offset + s.length > b.length) throw new Error('segment lies outside the file');
    if (!ZLIB_TYPES.has(s.type)) return b.subarray(q, s.offset + s.length);
    const flag = dv.getInt32(q, le), clen = dv.getInt32(q + 4, le), alg = b[q + 8];
    if (flag === 2 && alg === 2) { if (!(clen > 1 && q + 8 + clen <= b.length)) throw new Error('compressed segment is truncated'); return inflate(b.subarray(q + 9, q + 8 + clen), 512e6, 'deflate'); }
    if (flag === 3 && alg === 3) { if (!(clen > 1 && q + 8 + clen <= b.length)) throw new Error('compressed segment is truncated'); return unxz(b.subarray(q + 9, q + 8 + clen)); }   // JT 10: XZ / LZMA2
    return b.subarray(q + 9, s.offset + s.length);
  };
  const lsgSeg = segments.find((s) => s.id === lsgId) || segments.find((s) => s.type === 1); let lsg = null;
  if (lsgSeg) {
    const d = await payload(lsgSeg);
    for (const layout of major >= 10 ? ['v10', 'v9', 'v8'] : major >= 9 ? ['v9', 'v8', 'v10'] : ['v8', 'v9', 'v10']) { lsg = parseLSG(d, le, layout); if (lsg) break; }
  }
  if (lsg) {
    meta.lsgLayout = lsg.layout;
    for (const nd of lsg.nodes.values()) { meta.nodeCensus[nd.kind] = (meta.nodeCensus[nd.kind] || 0) + 1; nd.props = {}; for (const [k, v] of lsg.table.get(nd.id) || []) if (typeof k === 'string') nd.props[k] = v; const nm = nd.props.JT_PROP_NAME; if (typeof nm === 'string') nd.name = nm.replace(/\.(part|asm)?;\d+;\d+:?$/i, '').replace(/;\d+;\d+:?$/, ''); }
    for (const nd of lsg.nodes.values()) { if (nd.kind === 'part' && nd.name && meta.partNames.length < 500) meta.partNames.push(nd.name); const u = nd.props.JT_PROP_MEASUREMENT_UNITS; if (typeof u === 'string' && !meta.units) meta.units = u; }
  }
  // XT B-Rep segments: the Parasolid stream starts at its "PS" / "B" flag sequence inside the element
  const xtStreams = [];
  for (const s of segments) {
    if (s.type !== 17) continue;
    try { const d = await payload(s); let at = -1; for (let i = 0; i + 30 < d.length && i < 400; i++) if (((d[i] === 0x50 && d[i + 1] === 0x53 && d[i + 2] === 0 && d[i + 3] === 0) || d[i] === 0x42 || d[i] === 0x54) && /^.{1,12}: TRANSMIT FILE/.test(str(d, i, i + 40))) { at = i; break; } if (at >= 0) xtStreams.push({ segment: s.id, bytes: d.subarray(at) }); } catch { /* unreadable segment: counted in the census only */ }
  }
  return { meta, lsg, segments, xtStreams, payload, le, major, minor, guidAt: (d, q) => guidAt(d, q, le) };
}
/** Product of two scene-graph matrices in the row-vector convention: first `a`, then `b`. */
export function jtMul(a, b) { return mul4(a, b);
}

/**
 * Tessellated shapes to show: walks the scene graph accumulating transforms, takes only the first (most detailed)
 * child of every LOD node, and returns [{ segment, name, matrix, underXT }] for each tri-strip shape node, where
 * underXT says that the owning part also carries an XT B-Rep segment.
 */
export function jtShapes(lsg, lod = 0) {
  const out = [];
  const walk = (id, M, depth, path, underXT) => {
    const nd = lsg.nodes.get(id); if (!nd || depth > 64) return;
    let m = M; for (const a of nd.attrs) { const at = lsg.attrs.get(a); if (at && at.m) m = mul4(at.m, m); }
    const here = nd.name || path; let xt = underXT;
    if (nd.kind === 'part') xt = Object.values(nd.props || {}).some((v) => v && typeof v === 'object' && v.segmentType === 17);
    if (nd.kind === 'triStripShape') { for (const v of Object.values(nd.props || {})) if (v && typeof v === 'object' && v.segmentType >= 6 && v.segmentType <= 16) out.push({ segment: v.segment, name: here || null, matrix: m, underXT: xt }); return; }
    if (nd.kind === 'partition' && id !== lsg.root && nd.fileName && !nd.children.length) { out.push({ external: nd.fileName, name: here || null, matrix: m }); return; }   // geometry in a separate file
    const pickLod = Math.min(Math.max(0, lod | 0), Math.max(0, nd.children.length - 1)), kids = nd.kind === 'rangeLod' || nd.kind === 'lod' ? nd.children.slice(pickLod, pickLod + 1) : nd.children;
    for (const c of kids) walk(c, m, depth + 1, here, xt);
  };
  if (lsg) walk(lsg.root, I4(), 0, null, false);
  return out;
}

/**
 * Placement of every part that owns one of the given segments: walks the scene graph from the root, accumulating
 * geometric transforms (row-vector convention p' = p·A·M), and returns [{ segment, name, matrix }].
 */
export function jtPlacements(lsg, segmentIds) {
  const out = [], want = new Set(segmentIds);
  const walk = (id, M, depth, path) => {
    const nd = lsg.nodes.get(id); if (!nd || depth > 64) return;
    let m = M; for (const a of nd.attrs) { const at = lsg.attrs.get(a); if (at && at.m) m = mul4(at.m, m); }
    const here = nd.name || path;
    for (const v of Object.values(nd.props || {})) if (v && typeof v === 'object' && want.has(v.segment)) out.push({ segment: v.segment, name: here || null, matrix: m });
    for (const c of nd.children) walk(c, m, depth + 1, here);
  };
  if (lsg) walk(lsg.root, I4(), 0, null);
  return out;
}
