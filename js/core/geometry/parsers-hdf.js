// Readers for mesh files stored in scientific containers:
//   classic NetCDF (CDF-1/2/5, parsed natively) and NetCDF-4/HDF5 → Exodus II
//   HDF5 (through the vendored h5wasm kernel)                       → CGNS, MED (Salome), Fluent .msh.h5
// Legacy ADF-format CGNS is not HDF5 and stays metadata-only.

import { Mesh, ET, guard, str, view, structuredBlock, unitFromWord, LIMITS } from './parsers-util.js';
import { loadH5 } from './parsers-wasm.js';

const num = (v) => (typeof v === 'bigint' ? Number(v) : v);
/** Any numeric dataset value → plain typed array of numbers (64-bit integers are narrowed). */
const numeric = (v) => (v instanceof BigInt64Array || v instanceof BigUint64Array ? Float64Array.from(v, Number) : ArrayBuffer.isView(v) ? v : Array.isArray(v) ? Float64Array.from(v, (x) => Number(x)) : v === null || v === undefined ? new Float64Array(0) : Float64Array.of(Number(v)));
const cstr = (v) => (typeof v === 'string' ? v : ArrayBuffer.isView(v) ? String.fromCharCode(...Array.from(v.subarray(0, 4096), (c) => c & 255)) : Array.isArray(v) ? v.join('') : String(v ?? '')).replace(/\0[\s\S]*$/, '').trim();

// ---------- classic NetCDF ----------
const NC_SIZE = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8, 7: 1, 8: 2, 9: 4, 10: 8, 11: 8 };
/** Parse a classic NetCDF header. Returns the same accessor interface as the HDF5 wrapper below. */
export function netcdfClassic(b) {
  const dv = view(b), ver = b[3], wide = ver === 5, n = b.length; let p = 4;
  const need = (k) => { if (p + k > n) throw new Error('NetCDF header is truncated'); };
  const u32 = () => { need(4); const v = dv.getUint32(p); p += 4; return v; };
  const cnt = () => { if (!wide) return u32(); need(8); const v = Number(dv.getBigUint64(p)); p += 8; return v; };
  const name = () => { const k = cnt(); if (k > 4096) throw new Error('NetCDF name length is not plausible'); need(k); const s = str(b, p, p + k); p += k + ((4 - (k % 4)) % 4); return s; };
  const values = (type, k, at) => {
    const s = NC_SIZE[type]; if (!s) throw new Error(`NetCDF type ${type} is not recognised`);
    if (at + k * s > n) throw new Error('NetCDF data extends past the end of the file (truncated)');
    if (type === 2) return b.subarray(at, at + k);
    const out = type === 5 || type === 6 || type >= 10 ? new Float64Array(k) : new Int32Array(k);
    for (let i = 0, o = at; i < k; i++, o += s) out[i] = type === 1 ? dv.getInt8(o) : type === 3 ? dv.getInt16(o) : type === 4 ? dv.getInt32(o) : type === 5 ? dv.getFloat32(o) : type === 6 ? dv.getFloat64(o) : type === 7 ? dv.getUint8(o) : type === 8 ? dv.getUint16(o) : type === 9 ? dv.getUint32(o) : Number(type === 10 ? dv.getBigInt64(o) : dv.getBigUint64(o));
    return type === 9 ? Float64Array.from(out, (x) => x >>> 0) : out;
  };
  const atts = () => {
    const tag = u32(), k = cnt(), o = {}; if (tag !== 0 && tag !== 0x0c) throw new Error('NetCDF attribute list is malformed');
    for (let i = 0; i < guard(k, 'NetCDF attribute', 1e5); i++) { const nm = name(), type = u32(), ne = guard(cnt(), 'NetCDF attribute value', 1e8), s = NC_SIZE[type] || 1, v = values(type, ne, p); o[nm] = type === 2 ? cstr(v) : ne === 1 ? v[0] : Array.from(v); p += ne * s + ((4 - ((ne * s) % 4)) % 4); }
    return o;
  };
  const numrecs = cnt(), dims = [], vars = new Map();
  { const tag = u32(), k = cnt(); if (tag !== 0 && tag !== 0x0a) throw new Error('NetCDF dimension list is malformed'); for (let i = 0; i < guard(k, 'NetCDF dimension', 1e5); i++) dims.push({ name: name(), size: cnt() }); }
  const gatt = atts();
  { const tag = u32(), k = cnt(); if (tag !== 0 && tag !== 0x0b) throw new Error('NetCDF variable list is malformed');
    for (let i = 0; i < guard(k, 'NetCDF variable', 1e5); i++) {
      const nm = name(), nd = guard(cnt(), 'NetCDF variable dimension', 64), ids = []; for (let j = 0; j < nd; j++) ids.push(cnt());
      const at = atts(), type = u32(), vsize = cnt(); let begin; if (ver === 1) begin = u32(); else { need(8); begin = Number(dv.getBigUint64(p)); p += 8; }
      vars.set(nm, { ids, at, type, vsize, begin, record: nd > 0 && dims[ids[0]]?.size === 0 });
    } }
  let recsize = 0; for (const v of vars.values()) if (v.record) recsize += v.vsize;
  const dim = (nm) => { const d = dims.find((x) => x.name === nm); return d ? (d.size === 0 ? numrecs : d.size) : undefined; };
  return {
    container: `NetCDF classic (CDF-${ver})`, has: (nm) => vars.has(nm), dim, names: () => [...vars.keys()], dims: () => Object.fromEntries(dims.map((d) => [d.name, d.size === 0 ? numrecs : d.size])),
    attr: (nm, a) => (nm === null ? gatt[a] : vars.get(nm)?.at[a]),
    read(nm) {
      const v = vars.get(nm); if (!v) return null;
      const shape = v.ids.map((id, k) => { const d = dims[id]; if (!d) throw new Error('NetCDF variable refers to an unknown dimension'); return k === 0 && v.record ? numrecs : d.size; }), total = guard(shape.reduce((s, x) => s * x, 1), `NetCDF ${nm} value`, LIMITS.elems * 27);
      if (!v.record) return { data: values(v.type, total, v.begin), shape };
      const per = total / (numrecs || 1), s = NC_SIZE[v.type], out = v.type === 2 ? new Uint8Array(total) : new Float64Array(total);
      for (let r = 0; r < numrecs; r++) out.set(values(v.type, per, v.begin + r * recsize), r * per);
      return { data: out, shape, s };
    },
  };
}

// ---------- HDF5 access ----------
let h5seq = 0;
/** Open `bytes` as an HDF5 file in the kernel's in-memory file system, run fn(file, h5), always close and remove it. */
async function withH5(b, ctx, fn) {
  if (ctx.opts.wasm === false) throw new Error('HDF5 decoding was switched off (opts.wasm = false)');
  const h5 = await loadH5(), path = `/import_${++h5seq}.h5`;
  let f = null;
  try {
    try { h5.FS.writeFile(path, b); f = new h5.File(path, 'r'); }
    catch (e) { throw new Error(`the HDF5 library could not open the file (${String(e?.message || e).split('\n').find((l) => /\(\)|error/i.test(l))?.replace(/^.*?in /, '').trim() || 'corrupt or truncated container'})`); }
    return await fn(f, h5);
  } finally { try { if (f) f.close(); } catch { /* already closed */ } try { h5.FS.unlink(path); } catch { /* never written */ } }
}
const isGroup = (e) => !!e && (e.type === 'Group' || typeof e.keys === 'function');
const isData = (e) => !!e && e.type === 'Dataset';
const child = (g, k) => { try { return g.get(k); } catch { return null; } };
const attrOf = (e, a) => { try { const x = e.attrs?.[a]; return x ? x.value : undefined; } catch { return undefined; } };
const valueOf = (e) => { try { return isData(e) ? e.value : null; } catch { return null; } };
/** Fixed-width character table (rows × width, as bytes or as an array of 1-character strings) → row strings. */
function charRows(v, shape) {
  if (!v || !shape || shape.length < 2) return typeof v === 'string' ? [v] : Array.isArray(v) && shape?.length === 1 ? v.map(cstr) : [];
  const [rows, w] = shape, out = [];
  for (let r = 0; r < rows; r++) out.push(Array.isArray(v) ? cstr(v.slice(r * w, (r + 1) * w).map((c) => (c === '' ? '\0' : c)).join('')) : cstr(v.subarray(r * w, (r + 1) * w)));
  return out;
}

// ---------- Exodus II ----------
const EXO_SIDES = {
  hex: [[0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [0, 4, 7, 3], [0, 3, 2, 1], [4, 5, 6, 7]], tet: [[0, 1, 3], [1, 2, 3], [0, 3, 2], [0, 2, 1]],
  wedge: [[0, 1, 4, 3], [1, 2, 5, 4], [0, 3, 5, 2], [0, 2, 1], [3, 4, 5]], pyr: [[0, 1, 4], [1, 2, 4], [2, 3, 4], [0, 4, 3], [0, 3, 2, 1]],
  quad: [[0, 1], [1, 2], [2, 3], [3, 0]], tri: [[0, 1], [1, 2], [2, 0]],
};
function exoType(name, npe) {
  const t = String(name || '').toUpperCase();
  if (/^HEX/.test(t)) return { 8: 'hex8', 20: 'hex20', 27: 'hex27' }[npe];
  if (/^TET/.test(t)) return { 4: 'tet4', 10: 'tet10' }[npe];
  if (/^(WEDGE|PENTA|PRISM)/.test(t)) return { 6: 'wedge6', 15: 'wedge15', 18: 'wedge18' }[npe];
  if (/^PYR/.test(t)) return { 5: 'pyr5', 13: 'pyr13', 14: 'pyr14' }[npe];
  if (/^(QUAD|SHELL)/.test(t)) return { 4: 'quad4', 8: 'quad8', 9: 'quad9' }[npe];
  if (/^TRI/.test(t)) return { 3: 'tri3', 6: 'tri6' }[npe];
  if (/^(BAR|BEAM|TRUSS|EDGE|ROD)/.test(t)) return { 2: 'line2', 3: 'line3' }[npe];
  return undefined;
}
/** Build the model from an Exodus variable accessor (classic NetCDF or HDF5). */
function exodus(src, ctx) {
  const nn = src.dim('num_nodes'), nd = src.dim('num_dim') ?? 3, nblk = src.dim('num_el_blk') ?? 0;
  if (!(nn >= 0) || !(src.has('coordx') || src.has('coord'))) throw new Error('not an Exodus II database (no num_nodes / coordinate variables)');
  guard(nn, 'Exodus node', LIMITS.verts);
  const M = new Mesh(), meta = { container: src.container, title: src.attr(null, 'title') ?? null, apiVersion: src.attr(null, 'api_version') ?? null, version: src.attr(null, 'version') ?? null, dimension: nd, nodes: nn, elements: src.dim('num_elem') ?? 0, blocks: [], sideSets: [], nodeSets: [], timeSteps: src.dim('time_step') ?? 0, skippedBlocks: 0 };
  let X, Y, Z;
  if (src.has('coordx')) { X = src.read('coordx').data; Y = nd > 1 && src.has('coordy') ? src.read('coordy').data : null; Z = nd > 2 && src.has('coordz') ? src.read('coordz').data : null; }
  else { const c = src.read('coord').data; X = c.subarray(0, nn); Y = nd > 1 ? c.subarray(nn, 2 * nn) : null; Z = nd > 2 ? c.subarray(2 * nn, 3 * nn) : null; }
  if (X.length < nn) throw new Error('Exodus coordinate arrays are shorter than num_nodes');
  for (let i = 0; i < nn; i++) M.node(X[i], Y ? Y[i] : 0, Z ? Z[i] : 0);
  const ids = src.has('eb_prop1') ? src.read('eb_prop1').data : [], names = src.has('eb_names') ? (() => { const r = src.read('eb_names'); return charRows(r.data, r.shape); })() : [];
  const blocks = []; let e0 = 0;
  for (let k = 1; k <= guard(nblk, 'Exodus block', 1e6); k++) {
    const vn = `connect${k}`; if (!src.has(vn)) continue;
    const r = src.read(vn), ne = r.shape[0], npe = r.shape[1] ?? 1, tn = cstr(src.attr(vn, 'elem_type')), type = exoType(tn, npe), id = ids[k - 1] ?? k, nm = names[k - 1] || `block ${id}`;
    blocks.push({ start: e0, count: ne, npe, data: r.data, type }); e0 += ne;
    meta.blocks.push({ id, name: nm, elemType: tn, nodesPerElement: npe, elements: ne, read: !!type });
    if (!type) { meta.skippedBlocks++; continue; }
    const g = M.group('zone', k, nm); M.groups[g].tag = id;
    const a = new Array(npe); for (let e = 0; e < ne; e++) { for (let j = 0; j < npe; j++) { const v = r.data[e * npe + j] - 1; a[j] = v >= 0 && v < nn ? v : -1; } M.elem(type, a, g); }
  }
  // side sets → boundary faces (edges for planar meshes); node sets → census
  const ssIds = src.has('ss_prop1') ? src.read('ss_prop1').data : [], ssNames = src.has('ss_names') ? (() => { const r = src.read('ss_names'); return charRows(r.data, r.shape); })() : [];
  for (let k = 1; k <= (src.dim('num_side_sets') ?? 0) && k < 1e6; k++) {
    if (!src.has(`elem_ss${k}`) || !src.has(`side_ss${k}`)) continue;
    const el = src.read(`elem_ss${k}`).data, sd = src.read(`side_ss${k}`).data, nm = ssNames[k - 1] || `side set ${ssIds[k - 1] ?? k}`, g = M.group('boundary', 'ss' + k, nm); let miss = 0;
    for (let i = 0; i < el.length; i++) {
      const e = el[i] - 1, blk = blocks.find((q) => e >= q.start && e < q.start + q.count), shape = blk?.type?.replace(/\d+$/, ''), tab = shape && (ET[blk.type].dim === 3 || nd === 2) ? EXO_SIDES[shape] : null, f = tab?.[sd[i] - 1];
      if (!f) { miss++; continue; }
      const o = (e - blk.start) * blk.npe, nodes = f.map((j) => blk.data[o + j] - 1);
      M.elem(nodes.length === 2 ? 'line2' : nodes.length === 3 ? 'tri3' : 'quad4', nodes.map((v) => (v >= 0 && v < nn ? v : -1)), g);
    }
    meta.sideSets.push({ id: ssIds[k - 1] ?? k, name: nm, sides: el.length, unresolved: miss });
  }
  const nsIds = src.has('ns_prop1') ? src.read('ns_prop1').data : [], nsNames = src.has('ns_names') ? (() => { const r = src.read('ns_names'); return charRows(r.data, r.shape); })() : [];
  for (let k = 1; k <= (src.dim('num_node_sets') ?? 0) && k < 1e6; k++) { const c = src.dim(`num_nod_ns${k}`) ?? 0, nm = nsNames[k - 1] || `node set ${nsIds[k - 1] ?? k}`; meta.nodeSets.push({ id: nsIds[k - 1] ?? k, name: nm, nodes: c }); M.groups.push({ id: M.groups.length, name: nm, kind: 'boundary', count: 0, nodes: c, tag: nsIds[k - 1] ?? k }); }
  if (meta.skippedBlocks) ctx.warn(`${meta.skippedBlocks} element block(s) of an unmapped type (${meta.blocks.filter((q) => !q.read).map((q) => q.elemType).join(', ')}) were skipped.`);
  if (meta.timeSteps) ctx.warn(`The database holds ${meta.timeSteps} time step(s) of results; only the mesh is read.`);
  return M.result({ meta, kind: nd === 2 ? 'surface-mesh' : undefined });
}
/** HDF5 (NetCDF-4) variables behind the same accessor interface. */
function h5Vars(f) {
  const ds = (nm) => { const e = child(f, nm); return isData(e) ? e : null; };
  return {
    container: 'NetCDF-4 / HDF5', has: (nm) => !!ds(nm), names: () => f.keys(),
    dim: (nm) => { const d = ds(nm); return d ? d.shape?.[0] ?? undefined : undefined; },
    attr: (nm, a) => { const v = attrOf(nm === null ? f : ds(nm), a); return ArrayBuffer.isView(v) && v.length === 1 ? num(v[0]) : typeof v === 'bigint' ? Number(v) : Array.isArray(v) && v.length === 1 ? v[0] : v; },
    read(nm) { const d = ds(nm); if (!d) return null; const v = d.value, shape = d.shape || []; return { data: typeof v === 'string' || (Array.isArray(v) && typeof v[0] === 'string') ? v : numeric(v), shape }; },
  };
}
export async function readExodus(b, ctx) {
  if (b[0] === 0x43 && b[1] === 0x44 && b[2] === 0x46) return exodus(netcdfClassic(b), ctx);
  return withH5(b, ctx, (f) => exodus(h5Vars(f), ctx));
}

// ---------- CGNS (HDF5) ----------
const CGNS_ET = { 3: 'line2', 4: 'line3', 5: 'tri3', 6: 'tri6', 7: 'quad4', 8: 'quad8', 9: 'quad9', 10: 'tet4', 11: 'tet10', 12: 'pyr5', 13: 'pyr14', 14: 'wedge6', 15: 'wedge15', 16: 'wedge18', 17: 'hex8', 18: 'hex20', 19: 'hex27', 21: 'pyr13' };
const cgLabel = (g) => cstr(attrOf(g, 'label'));
const cgName = (g, key) => cstr(attrOf(g, 'name')) || key;
const cgKids = (g, label) => { const out = []; for (const k of g.keys()) { if (k[0] === ' ') continue; const c = child(g, k); if (isGroup(c) && (!label || cgLabel(c) === label)) out.push({ key: k, name: cgName(c, k), node: c }); } return out; };
const cgData = (g) => valueOf(child(g, ' data'));
const cgText = (g) => cstr(cgData(g));
export async function readCGNS(b, ctx) {
  if (b[0] !== 0x89) throw new Error('this CGNS file uses the legacy ADF container, which is not decoded; convert it with cgnsconvert -h (ADF → HDF5)');
  return withH5(b, ctx, (f) => {
    const M = new Mesh(), meta = { container: 'HDF5', cgnsVersion: null, bases: [], zones: [], families: [], boundaryConditions: [], polyhedral: false };
    let units = { length: null, source: null }, unsupported = 0, vertexBCs = 0, structuredBCs = 0;
    const ver = cgKids(f, 'CGNSLibraryVersion_t')[0]; if (ver) { const v = cgData(ver.node); meta.cgnsVersion = v ? +Number(numeric(v)[0]).toFixed(2) : null; }
    const bases = cgKids(f, 'CGNSBase_t'); if (!bases.length) throw new Error('no CGNSBase_t node was found (not a CGNS/HDF5 tree)');
    const unitFrom = (g) => { const du = cgKids(g, 'DimensionalUnits_t')[0]; if (!du) return; const v = cgData(du.node); if (!ArrayBuffer.isView(v) || v.length < 64) return; const w = cstr(v.subarray(32, 64)); meta.lengthUnit = w; const u = unitFromWord(w); if (u) units = { length: u, source: 'file' }; };
    for (const base of bases) {
      const bd = numeric(cgData(base.node) ?? []), cellDim = bd[0] ?? 3; meta.bases.push({ name: base.name, cellDimension: cellDim, physicalDimension: bd[1] ?? 3 }); unitFrom(base.node);
      for (const fam of cgKids(base.node, 'Family_t')) meta.families.push(fam.name);
      const zones = cgKids(base.node, 'Zone_t');
      for (const z of zones) {
        const zt = cgKids(z.node, 'ZoneType_t')[0], ztype = zt ? cgText(zt.node) : 'Structured', size = numeric(cgData(z.node) ?? []), prefix = zones.length > 1 || bases.length > 1 ? `${z.name} / ` : '';
        const gc = cgKids(z.node, 'GridCoordinates_t')[0]; if (!gc) { ctx.warn(`Zone "${z.name}" has no GridCoordinates and was skipped.`); continue; }
        unitFrom(z.node); unitFrom(gc.node);
        const arr = (re) => { const a = cgKids(gc.node, 'DataArray_t').find((c) => re.test(c.name)); return a ? numeric(cgData(a.node) ?? []) : null; };
        const X = arr(/^CoordinateX$/i) || arr(/^CoordinateR$/i), Y = arr(/^CoordinateY$/i) || arr(/^CoordinateTheta$/i), Z = arr(/^CoordinateZ$/i);
        if (!X) { ctx.warn(`Zone "${z.name}" has no CoordinateX array and was skipped.`); continue; }
        const base0 = M.nNodes, nv = X.length; guard(base0 + nv, 'CGNS vertex', LIMITS.verts);
        for (let i = 0; i < nv; i++) M.node(X[i], Y ? Y[i] : 0, Z ? Z[i] : 0);
        const zrec = { name: z.name, base: base.name, type: ztype, vertices: nv, sections: [], bcs: [] }; meta.zones.push(zrec);
        // boundary conditions given on elements/faces: [{ g, lo, hi } | { g, set }]
        const bcs = [];
        for (const zbc of cgKids(z.node, 'ZoneBC_t')) for (const bc of cgKids(zbc.node, 'BC_t')) {
          const loc = cgKids(bc.node, 'GridLocation_t')[0], where = loc ? cgText(loc.node) : 'Vertex', fam = cgKids(bc.node, 'FamilyName_t')[0], family = fam ? cgText(fam.node) : null, bctype = cgText(bc.node);
          const pr = cgKids(bc.node, 'IndexRange_t').find((c) => /^(PointRange|ElementRange)$/.test(c.name)), pl = cgKids(bc.node, 'IndexArray_t').find((c) => /^(PointList|ElementList)$/.test(c.name));
          zrec.bcs.push({ name: bc.name, type: bctype, family, location: where }); meta.boundaryConditions.push({ zone: z.name, name: bc.name, type: bctype, family, location: where });
          if (ztype !== 'Unstructured') { structuredBCs++; continue; }
          const onElems = /Center$/.test(where) || /^Element/.test(pr?.name || pl?.name || ''); if (!onElems) { vertexBCs++; continue; }
          const g = M.group('boundary', `${z.name}/${bc.name}`, `${prefix}${family || bc.name}`); M.groups[g].bcType = bctype; M.groups[g].family = family; if (family) M.groups[g].bcName = bc.name;
          if (pr) { const r = numeric(cgData(pr.node) ?? []); if (r.length >= 2) bcs.push({ g, lo: r[0], hi: r[1] }); } else if (pl) bcs.push({ g, set: new Set(numeric(cgData(pl.node) ?? [])) });
        }
        const bcOf = (no) => { for (const q of bcs) if (q.set ? q.set.has(no) : no >= q.lo && no <= q.hi) return q.g; return -1; };
        if (ztype === 'Structured') {
          const idim = Math.max(1, Math.round(size.length / 3)), d = [size[0] || 1, idim > 1 ? size[1] : 1, idim > 2 ? size[2] : 1];
          if (d[0] * d[1] * d[2] !== nv) { ctx.warn(`Structured zone "${z.name}": the zone size does not match the coordinate count; cells were not built.`); continue; }
          zrec.dimensions = d; structuredBlock(M, d[0], d[1], d[2], base0, M.group('zone', `${base.name}/${z.name}`, z.name));
          continue;
        }
        // unstructured: element sections
        const ngon = [], nface = [], at = (v) => (v >= 1 && v <= nv ? base0 + v - 1 : -1);
        for (const s of cgKids(z.node, 'Elements_t')) {
          const er = cgKids(s.node, 'IndexRange_t').find((c) => c.name === 'ElementRange'), sd = numeric(cgData(s.node) ?? []), et = sd[0], rng = er ? numeric(cgData(er.node) ?? []) : [], first = rng[0] ?? 1, count = (rng[1] ?? 0) - first + 1;
          const das = cgKids(s.node, 'DataArray_t'), cn = das.find((c) => c.name === 'ElementConnectivity'), so = das.find((c) => c.name === 'ElementStartOffset'), conn = cn ? numeric(cgData(cn.node) ?? []) : null, offs = so ? numeric(cgData(so.node) ?? []) : null;
          zrec.sections.push({ name: s.name, elementType: et, first, count });
          if (!conn || !(count > 0)) continue;
          guard(count, 'CGNS element');
          if (CGNS_ET[et]) {
            const t = CGNS_ET[et], npe = ET[t].n, sg = M.group(ET[t].dim === cellDim ? 'zone' : 'boundary', `${z.name}/${s.name}`, `${prefix}${s.name}`), a = new Array(npe);
            for (let e = 0; e < count && (e + 1) * npe <= conn.length; e++) { for (let j = 0; j < npe; j++) a[j] = at(conn[e * npe + j]); const bg = bcOf(first + e); M.elem(t, a, bg >= 0 ? bg : sg); }
          } else if (et === 20) {                                 // MIXED: [type, nodes…] per element
            const sg = M.group('zone', `${z.name}/${s.name}`, `${prefix}${s.name}`);
            for (let e = 0, p = 0; e < count && p < conn.length; e++) {
              const t = CGNS_ET[conn[offs ? offs[e] : p]], o = (offs ? offs[e] : p) + 1; if (!t) { unsupported += count - e; break; }
              const npe = ET[t].n, a = new Array(npe); for (let j = 0; j < npe; j++) a[j] = at(conn[o + j]);
              const bg = bcOf(first + e); M.elem(t, a, bg >= 0 ? bg : ET[t].dim === cellDim ? sg : M.group('boundary', `${z.name}/${s.name}/b`, `${prefix}${s.name} (faces)`)); p = o + npe;
            }
          } else if (et === 22 || et === 23) {                     // NGON_n / NFACE_n: [n, ids…] (CGNS 3) or offsets (CGNS 4)
            const lists = []; if (offs) for (let e = 0; e + 1 < offs.length && e < count; e++) lists.push(conn.subarray(offs[e], offs[e + 1])); else for (let e = 0, p = 0; e < count && p < conn.length; e++) { const k = conn[p]; lists.push(conn.subarray(p + 1, p + 1 + k)); p += 1 + k; }
            (et === 22 ? ngon : nface).push({ first, lists, name: s.name });
          } else unsupported += count;
        }
        if (ngon.length) {                                         // polyhedral zone: show the faces that bound exactly one cell
          meta.polyhedral = true; const use = new Map(); let cells = 0;
          for (const nf of nface) for (const l of nf.lists) { cells++; for (const v of l) { const id = Math.abs(v); use.set(id, (use.get(id) || 0) + 1); } }
          const sg = M.group('boundary', `${z.name}/ngon`, `${prefix}polyhedral boundary`);
          for (const ng of ngon) ng.lists.forEach((l, i) => { const no = ng.first + i; if (nface.length && use.get(no) !== 1) return; const bg = bcOf(no), g = bg >= 0 ? bg : sg, n = Array.from(l, at); if (n.length === 3) M.elem('tri3', n, g); else if (n.length === 4) M.elem('quad4', n, g); else for (let j = 2; j < n.length; j++) M.elem('tri3', [n[0], n[j - 1], n[j]], g); });
          zrec.polyhedralCells = cells;
          ctx.warn(`Zone "${z.name}" is polyhedral (NGON_n${nface.length ? ' / NFACE_n' : ''}): ${nface.length ? 'its boundary faces are shown as a surface' : 'ALL faces are shown because no NFACE_n section defines the cells'}; polyhedral cells are not converted to volume elements.`);
        }
      }
    }
    if (unsupported) ctx.warn(`${unsupported} element(s) of unsupported CGNS types (e.g. cubic or user-defined) were skipped.`);
    if (vertexBCs) ctx.warn(`${vertexBCs} boundary condition(s) are defined on vertices; they are listed in the metadata but not mapped to elements.`);
    if (structuredBCs) ctx.warn(`${structuredBCs} boundary condition(s) on structured zones are listed in the metadata but not mapped to cell faces.`);
    if (!M.nNodes) throw Object.assign(new Error('the CGNS tree holds no zone with grid coordinates'), { meta });
    const hasVol = [...M.blocks.keys()].some((t) => ET[t].dim === 3);
    return M.result({ meta, units, kind: meta.polyhedral && !hasVol ? 'surface-mesh' : undefined });
  });
}

// ---------- MED (Salome) ----------
const MED_ET = { SE2: 'line2', SE3: 'line3', TR3: 'tri3', TR6: 'tri6', QU4: 'quad4', QU8: 'quad8', QU9: 'quad9', TE4: 'tet4', T10: 'tet10', HE8: 'hex8', H20: 'hex20', H27: 'hex27', PE6: 'wedge6', P15: 'wedge15', PY5: 'pyr5', P13: 'pyr13' };
const ORIENT = { tet: [[1, 0, 2, 0, 3, 0], [[0, 1]]], hex: [[1, 0, 3, 0, 4, 0], [[0, 4], [1, 5], [2, 6], [3, 7]]], wedge: [[1, 0, 2, 0, 3, 0], [[0, 3], [1, 4], [2, 5]]], pyr: [[1, 0, 3, 0, 4, 0], [[1, 3]]] };
/** MED lists volume cells with the opposite handedness to the platform's convention: flip blocks that are uniformly negative. */
function toPositiveOrientation(M) {
  const X = M.xyz, flipped = [];
  for (const [type, blk] of M.blocks) {
    const o = ORIENT[type.replace(/\d+$/, '')]; if (!o || ET[type].dim !== 3) continue;
    const np = ET[type].n, ne = blk.group.length, c = blk.conn, [q, sw] = o; let neg = 0;
    for (let e = 0; e < ne; e++) { const v = (k) => { const a = c[e * np + q[2 * k]], b0 = c[e * np + q[2 * k + 1]]; return [X[3 * a] - X[3 * b0], X[3 * a + 1] - X[3 * b0 + 1], X[3 * a + 2] - X[3 * b0 + 2]]; }, u = v(0), w = v(1), h = v(2); if (u[0] * (w[1] * h[2] - w[2] * h[1]) - u[1] * (w[0] * h[2] - w[2] * h[0]) + u[2] * (w[0] * h[1] - w[1] * h[0]) < 0) neg++; }
    if (ne && neg === ne) { for (let e = 0; e < ne; e++) for (const [a, b0] of sw) { const t = c[e * np + a]; c[e * np + a] = c[e * np + b0]; c[e * np + b0] = t; } flipped.push(`${ne} ${type}`); }
  }
  return flipped;
}
export async function readMED(b, ctx) {
  return withH5(b, ctx, (f) => {
    const info = child(f, 'INFOS_GENERALES'), meta = { container: 'HDF5', medVersion: info ? ['MAJ', 'MIN', 'REL'].map((a) => num(numeric(attrOf(info, a))[0])).join('.') : null, meshes: [], families: 0, skipped: {} };
    const ens = child(f, 'ENS_MAA'); if (!isGroup(ens) || !ens.keys().length) throw Object.assign(new Error('no mesh (ENS_MAA) was found in the MED file'), { meta });
    const M = new Mesh(); let units = { length: null, source: null };
    for (const mname of ens.keys()) {
      const mg = child(ens, mname); if (!isGroup(mg)) continue;
      const esp = num(numeric(attrOf(mg, 'ESP'))[0]) || 3, step = mg.keys().map((k) => child(mg, k)).find(isGroup);
      const uni = cstr(attrOf(mg, 'UNI')).split(/\s+/).filter(Boolean), u = uni.length && uni.every((w) => w === uni[0]) ? unitFromWord(uni[0]) : null; if (u) units = { length: u, source: 'file' };
      const rec = { name: mname, dimension: num(numeric(attrOf(mg, 'DIM'))[0]) || null, spaceDimension: esp, description: cstr(attrOf(mg, 'DES')) || null, unit: uni[0] ?? null, nodes: 0, elements: {} }; meta.meshes.push(rec);
      const coo = step && child(step, 'NOE') && valueOf(child(child(step, 'NOE'), 'COO'));
      if (!coo) { ctx.warn(`MED mesh "${mname}" has no node coordinates (structured MED grids are not read) and was skipped.`); continue; }
      const C = numeric(coo), nn = Math.floor(C.length / esp), base0 = M.nNodes; guard(base0 + nn, 'MED node', LIMITS.verts); rec.nodes = nn;
      for (let i = 0; i < nn; i++) M.node(C[i], esp > 1 ? C[nn + i] : 0, esp > 2 ? C[2 * nn + i] : 0);          // coordinates are stored axis by axis
      // families: number → first group name
      const famName = new Map(), fas = child(child(f, 'FAS') || f, mname), el = isGroup(fas) ? child(fas, 'ELEME') : null;
      if (isGroup(el)) for (const k of el.keys()) { const fg = child(el, k); if (!isGroup(fg)) continue; const no = num(numeric(attrOf(fg, 'NUM'))[0]), gro = child(fg, 'GRO'), nom = isGroup(gro) ? child(gro, 'NOM') : null, v = valueOf(nom), names = ArrayBuffer.isView(v) ? charRows(v, [Math.floor(v.length / 80), 80]) : Array.isArray(v) ? v.map(cstr) : []; famName.set(no, names.filter(Boolean).join(' + ') || k); meta.families++; }
      const mai = step ? child(step, 'MAI') : null;
      if (isGroup(mai)) for (const tk of mai.keys()) {
        const tg = child(mai, tk), nod = isGroup(tg) ? valueOf(child(tg, 'NOD')) : null, t = MED_ET[tk];
        if (!nod) continue;
        if (!t) { meta.skipped[tk] = (meta.skipped[tk] || 0) + 1; continue; }
        const N = numeric(nod), np = ET[t].n, ne = Math.floor(N.length / np), fam = numeric(valueOf(child(tg, 'FAM')) ?? []), a = new Array(np); guard(ne, 'MED element'); rec.elements[tk] = ne;
        for (let e = 0; e < ne; e++) { for (let j = 0; j < np; j++) { const v = N[j * ne + e]; a[j] = v >= 1 && v <= nn ? base0 + v - 1 : -1; } const fn = fam.length > e ? fam[e] : 0; M.elem(t, a, fn !== 0 ? M.group('component', `${mname}/${fn}`, famName.get(fn) ?? `family ${fn}`) : -1); }
      }
    }
    if (!M.nNodes) throw Object.assign(new Error('the MED file holds no unstructured mesh with node coordinates'), { meta });
    const sk = Object.keys(meta.skipped); if (sk.length) ctx.warn(`MED element types not mapped and skipped: ${sk.join(', ')} (points, polygons and polyhedra are not converted).`);
    const flipped = toPositiveOrientation(M); meta.reorderedBlocks = flipped;
    if (flipped.length) ctx.warn(`MED lists volume cells with the opposite handedness: ${flipped.join(', ')} were renumbered to the platform's node-ordering convention (geometry unchanged).`);
    return M.result({ meta, units, kind: [...M.blocks.keys()].some((t) => ET[t].dim === 3) ? 'structural-mesh' : undefined });
  });
}

// ---------- Fluent CFF mesh (.msh.h5 / .cas.h5) ----------
export async function readFluentH5(b, ctx) {
  return withH5(b, ctx, (f) => {
    const meshes = child(f, 'meshes'), mk = isGroup(meshes) ? meshes.keys()[0] : null, mesh = mk ? child(meshes, mk) : null, meta = { container: 'HDF5 (Fluent CFF)', dimension: null, nodes: 0, faces: 0, cells: 0, zones: [], volumeCellsReconstructed: false };
    if (!isGroup(mesh)) throw Object.assign(new Error('no /meshes group was found — this HDF5 file is not a Fluent CFF mesh/case'), { meta });
    for (const [k, a] of [['dimension', 'dimension'], ['nodes', 'nodeCount'], ['faces', 'faceCount'], ['cells', 'cellCount']]) { const v = attrOf(mesh, a); if (v !== undefined) meta[k] = num(numeric(v)[0]); }
    const sections = (g) => (isGroup(g) ? g.keys().map((k) => ({ k, e: child(g, k) })).sort((x, y) => (+x.k || 0) - (+y.k || 0)) : []);
    const M = new Mesh(), nodesG = child(mesh, 'nodes'), facesG = child(mesh, 'faces');
    for (const s of sections(isGroup(nodesG) ? child(nodesG, 'coords') : null)) {
      if (!isData(s.e)) continue; const v = numeric(s.e.value), d = s.e.shape?.[1] ?? meta.dimension ?? 3; meta.dimension ??= d;
      guard(M.nNodes + v.length / d, 'Fluent node', LIMITS.verts); for (let i = 0; i + d <= v.length; i += d) M.node(v[i], v[i + 1], d > 2 ? v[i + 2] : 0);
    }
    if (!M.nNodes) throw Object.assign(new Error('no node coordinates (/meshes/*/nodes/coords) were found'), { meta });
    const nv = M.nNodes;
    // face zones: id / name / minId / maxId / zoneType
    const zt = isGroup(facesG) ? child(facesG, 'zoneTopology') : null, zones = [];
    if (isGroup(zt)) {
      const col = (k) => numeric(valueOf(child(zt, k)) ?? []), id = col('id'), lo = col('minId'), hi = col('maxId'), ty = col('zoneType'), nmv = valueOf(child(zt, 'name')), names = (Array.isArray(nmv) ? nmv.map(cstr).join(';') : cstr(nmv)).split(';').map((s) => s.trim());
      for (let i = 0; i < id.length; i++) zones.push({ id: id[i], lo: lo[i], hi: hi[i], type: ty[i] ?? null, name: names[i] || `zone ${id[i]}` });
    }
    const cz = isGroup(child(mesh, 'cells')) ? child(child(mesh, 'cells'), 'zoneTopology') : null;
    if (isGroup(cz)) { const id = numeric(valueOf(child(cz, 'id')) ?? []), lo = numeric(valueOf(child(cz, 'minId')) ?? []), hi = numeric(valueOf(child(cz, 'maxId')) ?? []), nmv = valueOf(child(cz, 'name')), names = (Array.isArray(nmv) ? nmv.map(cstr).join(';') : cstr(nmv)).split(';'); for (let i = 0; i < id.length; i++) { const c = (hi[i] ?? 0) - (lo[i] ?? 1) + 1; meta.zones.push({ id: id[i], kind: 'cells', name: names[i]?.trim() || null, count: c }); M.groups.push({ id: M.groups.length, name: names[i]?.trim() || `cell zone ${id[i]}`, kind: 'zone', count: c, tag: id[i] }); } }
    // faces: per section a node-count array and a flat 1-based node list; c1 = 0 marks a boundary face
    const cat = (g) => { const parts = sections(g).filter((s) => isData(s.e)).map((s) => numeric(s.e.value)); if (parts.length === 1) return parts[0]; const out = new Float64Array(parts.reduce((s, p) => s + p.length, 0)); let o = 0; for (const p of parts) { out.set(p, o); o += p.length; } return out; };
    const c1 = isGroup(facesG) ? cat(child(facesG, 'c1')) : new Float64Array(0); let fid = 0, internal = 0;
    for (const s of sections(isGroup(facesG) ? child(facesG, 'nodes') : null)) {
      if (!isGroup(s.e)) continue; const nn = numeric(valueOf(child(s.e, 'nnodes')) ?? []), nd = numeric(valueOf(child(s.e, 'nodes')) ?? []);
      for (let i = 0, p = 0; i < nn.length && p + nn[i] <= nd.length; i++, p += nn[i]) {
        const no = ++fid, zone = zones.find((q) => no >= q.lo && no <= q.hi), boundary = c1.length >= no ? c1[no - 1] === 0 : zone ? zone.type !== 2 : true;
        if (!boundary) { internal++; continue; }
        const k = nn[i], n = new Array(k); for (let j = 0; j < k; j++) { const v = nd[p + j]; n[j] = v >= 1 && v <= nv ? v - 1 : -1; }
        const g = zone ? M.group('boundary', zone.id, zone.name) : -1;
        if (k === 2) M.elem('line2', n, g); else if (k === 3) M.elem('tri3', n, g); else if (k === 4) M.elem('quad4', n, g); else for (let j = 2; j < k; j++) M.elem('tri3', [n[0], n[j - 1], n[j]], g);
      }
    }
    for (const q of zones) meta.zones.push({ id: q.id, kind: 'faces', name: q.name, type: q.type, count: q.hi - q.lo + 1 });
    meta.faces ||= fid; meta.internalFaces = internal;
    ctx.warn(`Fluent CFF mesh: ${meta.nodes || nv} nodes, ${meta.faces} faces, ${meta.cells} cells. Boundary faces are shown as a surface; cell connectivity is NOT reconstructed, so this model carries no volume elements.`);
    return M.result({ meta, kind: M.blocks.size ? 'surface-mesh' : 'pointcloud' });
  });
}

/** Generic HDF5 container: list the top of the object tree so the user can tell what the file is. */
export async function readHDF5Tree(b, ctx) {
  const meta = { container: 'HDF5', superblockVersion: b[8], tree: [] };
  try {
    await withH5(b, ctx, (f) => {
      const walk = (g, path, depth) => { for (const k of g.keys()) { if (meta.tree.length >= 200) return; const e = child(g, k), p = `${path}/${k}`; if (isData(e)) meta.tree.push(`${p}  [${(e.shape || []).join('×')}] ${typeof e.dtype === 'string' ? e.dtype : 'compound'}`); else if (isGroup(e)) { meta.tree.push(p + '/'); if (depth < 3) walk(e, p, depth + 1); } } };
      walk(f, '', 0);
      meta.rootAttributes = Object.keys(f.attrs || {});
    });
    ctx.warn('Generic HDF5 container: the object tree is listed in the metadata, but its layout is not one of the mesh schemas read here (CGNS, MED, Exodus/NetCDF-4, Fluent CFF).');
  } catch (e) { ctx.warn(`The HDF5 object tree could not be listed — ${e.message}.`); }
  return { kind: 'metadata-only', meta };
}

/** Entry point for every HDF5 / NetCDF mesh container: the schema is taken from the content, not the extension. */
export async function readHDF5Any(b, ctx) {
  if (b[0] === 0x43 && b[1] === 0x44 && b[2] === 0x46) return { ...exodus(netcdfClassic(b), ctx), formatId: 'exodus' };
  if (b[0] !== 0x89 && str(b, 0, 16).startsWith('@(#)ADF')) throw Object.assign(new Error('this CGNS file uses the legacy ADF container, which is not decoded; convert it with `cgnsconvert -h` (ADF → HDF5)'), { meta: { container: 'ADF', adfVersion: /ADF Database Version\s+(\S+)/.exec(str(b, 0, 64))?.[1] ?? null } });
  if (b[0] !== 0x89) throw Object.assign(new Error('the file has neither an HDF5 nor a NetCDF signature'), { meta: { container: 'unconfirmed (extension only)' } });
  let schema;
  try {
    schema = await withH5(b, ctx, (f) => {
      const keys = f.keys();
      if (keys.includes('ENS_MAA')) return 'med';
      if (keys.includes('meshes') && isGroup(child(f, 'meshes'))) return 'fluent-h5';
      if (keys.includes(' format') || keys.includes('CGNSLibraryVersion') || keys.some((k) => k[0] !== ' ' && cgLabel(child(f, k)) === 'CGNSBase_t')) return 'cgns';
      if (keys.includes('coordx') || keys.includes('coord') || keys.includes('num_nodes')) return 'exodus';
      return 'hdf5';
    });
  } catch (e) { throw Object.assign(e, { meta: { container: 'HDF5', superblockVersion: b[8] } }); }
  const part = schema === 'med' ? await readMED(b, ctx) : schema === 'fluent-h5' ? await readFluentH5(b, ctx) : schema === 'cgns' ? await readCGNS(b, ctx) : schema === 'exodus' ? await readExodus(b, ctx) : await readHDF5Tree(b, ctx);
  return { ...part, formatId: schema };
}
