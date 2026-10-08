// Readers for tessellated surfaces and point clouds:
// STL, OBJ, PLY, OFF, glTF/GLB, 3MF, AMF, COLLADA, VRML/X3D, XYZ, PCD, LAS and the platform's own JSON.
// Each reader takes (bytes, ctx) and returns a partial model; parsers.js validates and completes it.

import { Scanner, Mesh, guard, str, text, nums, toNum, view, xmlParse, xmlAll, xmlFirst, xmlChild, zipEntries, base64, M4, unitFromWord, LIMITS } from './parsers-util.js';

// ---------- STL ----------
export function readSTL(b, ctx) {
  const n = b.length >= 84 ? view(b).getUint32(80, true) : -1;
  const low = str(b, 0, Math.min(b.length, 4096)).trimStart().toLowerCase();
  const asciiLike = low.startsWith('solid') && /\bfacet\b|\bendsolid\b/.test(low);
  // A binary STL is recognised by its size; many exporters write a header that begins with "solid".
  if (n >= 0 && (84 + 50 * n === b.length && !(asciiLike && n === 0) || (!asciiLike && n > 0 && 84 + 50 * n <= b.length))) {
    guard(n, 'STL triangle', LIMITS.elems / 3);
    if (84 + 50 * n < b.length) ctx.warn(`${b.length - 84 - 50 * n} trailing bytes after the last STL triangle were ignored.`);
    const dv = view(b), pos = new Float64Array(9 * n), tri = new Uint32Array(3 * n);
    let attr = 0;
    for (let t = 0, o = 84; t < n; t++, o += 50) {
      for (let k = 0; k < 9; k++) pos[9 * t + k] = dv.getFloat32(o + 12 + 4 * k, true);
      if (dv.getUint16(o + 48, true)) attr++;
    }
    for (let i = 0; i < 3 * n; i++) tri[i] = i;
    const header = str(b, 0, 80).replace(/[^\x20-\x7e]+/g, ' ').trim();
    return { positions: pos, triangles: tri, kind: 'surface', meta: { encoding: 'binary', header, solidPrefixedBinary: header.toLowerCase().startsWith('solid'), trianglesWithAttributeBytes: attr, sharedVertices: false } };
  }
  if (!asciiLike) throw new Error('neither a consistent binary STL (84 + 50·n bytes) nor an ASCII STL with facets');
  const sc = new Scanner(b), xyz = [], solids = [];
  for (let t; (t = sc.token()) !== null;) {
    if (t.length === 6 && (t === 'vertex' || t === 'VERTEX' || t.toLowerCase() === 'vertex')) xyz.push(sc.num(), sc.num(), sc.num());
    else if (t.length === 5 && t.toLowerCase() === 'solid') solids.push({ name: (sc.line() ?? '').trim(), start: xyz.length / 9 });
  }
  const extra = xyz.length % 9;
  if (extra) { ctx.warn('ASCII STL ends inside a facet (truncated file?); the incomplete facet was dropped.'); xyz.length -= extra; }
  const nt = xyz.length / 9, tri = new Uint32Array(3 * nt);
  for (let i = 0; i < 3 * nt; i++) tri[i] = i;
  const out = { positions: Float64Array.from(xyz), triangles: tri, kind: 'surface', meta: { encoding: 'ascii', sharedVertices: false } };
  out.groups = solids.map((s, i) => ({ id: i, name: s.name || `solid ${i + 1}`, kind: 'solid', count: (i + 1 < solids.length ? solids[i + 1].start : nt) - s.start }));
  if (solids.length > 1) {
    const g = new Int32Array(nt); let k = 0;
    for (let t = 0; t < nt; t++) { while (k + 1 < solids.length && solids[k + 1].start <= t) k++; g[t] = k; }
    out.elements = [{ type: 'tri3', nodesPer: 3, count: nt, conn: tri.slice(), group: g }];
  }
  return out;
}

// ---------- OBJ ----------
export function readOBJ(b, ctx) {
  const sc = new Scanner(b), M = new Mesh(), lines = [], mtl = new Set(), libs = [];
  let g = -1, quads = 0, ngons = 0, bad = 0;
  for (let ln; (ln = sc.line()) !== null;) {
    const c0 = ln.charCodeAt(0), c1 = ln.charCodeAt(1);
    if (c0 === 118 && (c1 === 32 || c1 === 9)) { const p = ln.split(/\s+/); M.node(toNum(p[1] ?? ''), toNum(p[2] ?? ''), toNum(p[3] ?? '')); }
    else if ((c0 === 102 || c0 === 108) && (c1 === 32 || c1 === 9)) {
      const p = ln.trim().split(/\s+/), nv = M.nNodes, idx = [];
      for (let i = 1; i < p.length; i++) { const v = parseInt(p[i], 10); idx.push(v < 0 ? nv + v : v - 1); }
      if (idx.some((v) => !(v >= 0 && v < nv))) { bad++; continue; }
      if (c0 === 108) { for (let i = 1; i < idx.length; i++) lines.push(idx[i - 1], idx[i]); continue; }
      if (idx.length === 4) quads++; else if (idx.length > 4) ngons++;
      for (let i = 2; i < idx.length; i++) M.elem('tri3', [idx[0], idx[i - 1], idx[i]], g);
    } else if ((c0 === 103 || c0 === 111) && (c1 === 32 || c1 === 9 || ln.length === 1)) { const name = ln.slice(1).trim() || 'default'; g = M.group('component', name, name); }
    else if (ln.startsWith('usemtl')) mtl.add(ln.slice(6).trim());
    else if (ln.startsWith('mtllib')) libs.push(ln.slice(6).trim());
  }
  if (bad) ctx.warn(`${bad} face/line statement(s) referenced vertices that do not exist and were skipped.`);
  return M.result({ lines: Uint32Array.from(lines), kind: 'surface', meta: { quadsTriangulated: quads, polygonsTriangulated: ngons, materials: [...mtl], materialLibraries: libs } });
}

// ---------- PLY ----------
const PLY_SIZE = { char: 1, int8: 1, uchar: 1, uint8: 1, short: 2, int16: 2, ushort: 2, uint16: 2, int: 4, int32: 4, uint: 4, uint32: 4, float: 4, float32: 4, double: 8, float64: 8 };
export function readPLY(b, ctx) {
  const sc = new Scanner(b), elems = [], comments = [];
  let fmt = null, ended = false;
  for (let ln, k = 0; k < 10000 && (ln = sc.line()) !== null; k++) {
    const p = ln.trim().split(/\s+/);
    if (p[0] === 'format') fmt = p[1];
    else if (p[0] === 'comment' || p[0] === 'obj_info') comments.push(ln.slice(p[0].length).trim());
    else if (p[0] === 'element') elems.push({ name: p[1], count: guard(parseInt(p[2], 10), `PLY ${p[1]}`), props: [] });
    else if (p[0] === 'property' && elems.length) {
      if (p[1] === 'list') elems[elems.length - 1].props.push({ list: true, ct: p[2], it: p[3], name: p[4] });
      else elems[elems.length - 1].props.push({ list: false, type: p[1], name: p[2] });
    } else if (p[0] === 'end_header') { ended = true; break; }
  }
  if (!ended || !fmt) throw new Error('PLY header is incomplete (no format / end_header)');
  for (const e of elems) for (const q of e.props) for (const t of q.list ? [q.ct, q.it] : [q.type]) if (!PLY_SIZE[t]) throw new Error(`PLY property type "${t}" is not recognised`);
  const ascii = fmt === 'ascii', le = fmt === 'binary_little_endian';
  if (!ascii && !le && fmt !== 'binary_big_endian') throw new Error(`unknown PLY format "${fmt}"`);
  const dv = view(b); let p = sc.p;
  const rd = (t) => {
    if (ascii) return sc.float('PLY value');
    const s = PLY_SIZE[t];
    if (p + s > b.length) throw new Error('PLY data ends before the counts in the header are satisfied (truncated file)');
    let v;
    switch (t) {
      case 'char': case 'int8': v = dv.getInt8(p); break; case 'uchar': case 'uint8': v = dv.getUint8(p); break;
      case 'short': case 'int16': v = dv.getInt16(p, le); break; case 'ushort': case 'uint16': v = dv.getUint16(p, le); break;
      case 'int': case 'int32': v = dv.getInt32(p, le); break; case 'uint': case 'uint32': v = dv.getUint32(p, le); break;
      case 'float': case 'float32': v = dv.getFloat32(p, le); break; default: v = dv.getFloat64(p, le);
    }
    p += s; return v;
  };
  let pos = new Float64Array(0); const tri = []; let nv = 0, polys = 0, bad = 0;
  for (const e of elems) {
    const fixed = e.props.every((q) => !q.list) ? e.props.reduce((s, q) => s + PLY_SIZE[q.type], 0) : 0;
    if (!ascii && fixed && p + fixed * e.count > b.length) throw new Error(`PLY ${e.name} data is truncated (${e.count} records declared)`);
    if (e.name === 'vertex') {
      guard(e.count, 'PLY vertex', LIMITS.verts);
      const ix = e.props.findIndex((q) => q.name === 'x'), iy = e.props.findIndex((q) => q.name === 'y'), iz = e.props.findIndex((q) => q.name === 'z');
      if (ix < 0 || iy < 0) throw new Error('PLY vertex element has no x/y properties');
      pos = new Float64Array(3 * e.count); nv = e.count;
      for (let i = 0; i < e.count; i++) for (let k = 0; k < e.props.length; k++) {
        const q = e.props[k];
        if (q.list) { const c = rd(q.ct); for (let j = 0; j < c; j++) rd(q.it); continue; }
        const v = rd(q.type);
        if (k === ix) pos[3 * i] = v; else if (k === iy) pos[3 * i + 1] = v; else if (k === iz) pos[3 * i + 2] = v;
      }
    } else {
      const isFace = e.name === 'face';
      for (let i = 0; i < e.count; i++) for (const q of e.props) {
        if (!q.list) { rd(q.type); continue; }
        const c = rd(q.ct);
        if (!(c >= 0 && c <= 1e6)) throw new Error('PLY list length is not plausible (corrupt file)');
        if (isFace && (q.name === 'vertex_indices' || q.name === 'vertex_index')) {
          let a = -1, prev = -1;
          if (c > 3) polys++;
          for (let j = 0; j < c; j++) {
            const v = rd(q.it);
            if (j === 0) a = v; else if (j > 1) { if (a >= 0 && a < nv && prev >= 0 && prev < nv && v >= 0 && v < nv) tri.push(a, prev, v); else bad++; }
            prev = v;
          }
        } else for (let j = 0; j < c; j++) rd(q.it);
      }
    }
  }
  if (bad) ctx.warn(`${bad} PLY triangle(s) referenced vertices outside the vertex list and were skipped.`);
  return { positions: pos, triangles: Uint32Array.from(tri), kind: tri.length ? 'surface' : 'pointcloud', meta: { encoding: fmt, comments, elements: elems.map((e) => ({ name: e.name, count: e.count, properties: e.props.map((q) => q.name) })), polygonsTriangulated: polys } };
}

// ---------- OFF ----------
export function readOFF(b, ctx) {
  const sc = new Scanner(b);
  const next = () => { for (let ln; (ln = sc.line()) !== null;) { const h = ln.indexOf('#'); if (h >= 0) ln = ln.slice(0, h); ln = ln.trim(); if (ln) return ln; } return null; };
  let ln = next();
  if (ln === null) throw new Error('empty OFF file');
  const m = /^((?:ST)?C?N?4?n?OFF)\b\s*(.*)$/.exec(ln);
  let variant = 'OFF';
  if (m) { variant = m[1]; ln = m[2] || next(); }
  if (/4|n/.test(variant)) throw new Error(`${variant} (higher-dimensional OFF) is not supported`);
  const c = nums(ln ?? ''), nv = guard(c[0], 'OFF vertex', LIMITS.verts), nf = guard(c[1], 'OFF face');
  const pos = new Float64Array(3 * nv), tri = []; let polys = 0, bad = 0;
  for (let i = 0; i < nv; i++) { ln = next(); if (ln === null) throw new Error('OFF file ends inside the vertex list'); const v = nums(ln); pos[3 * i] = v[0]; pos[3 * i + 1] = v[1]; pos[3 * i + 2] = v[2] ?? 0; }
  for (let i = 0; i < nf; i++) {
    ln = next(); if (ln === null) { ctx.warn(`OFF file ends after ${i} of ${nf} faces.`); break; }
    const v = nums(ln), k = v[0];
    if (!(k >= 3 && k < v.length)) { bad++; continue; }
    if (k > 3) polys++;
    for (let j = 3; j <= k; j++) { const a = v[1], p = v[j - 1], q = v[j]; if (a >= 0 && a < nv && p >= 0 && p < nv && q >= 0 && q < nv) tri.push(a, p, q); else bad++; }
  }
  if (bad) ctx.warn(`${bad} OFF face record(s) were malformed or out of range and were skipped.`);
  return { positions: pos, triangles: Uint32Array.from(tri), kind: 'surface', meta: { variant, polygonsTriangulated: polys } };
}

// ---------- glTF 2.0 / GLB ----------
const GL_SIZE = { 5120: 1, 5121: 1, 5122: 2, 5123: 2, 5125: 4, 5126: 4 }, GL_COMP = { SCALAR: 1, VEC2: 2, VEC3: 3, VEC4: 4, MAT4: 16 };
export function readGLTF(b, ctx) {
  let json, bin = null;
  if (b[0] === 0x67 && b[1] === 0x6c && b[2] === 0x54 && b[3] === 0x46) {
    const dv = view(b), ver = dv.getUint32(4, true);
    if (ver !== 2) throw new Error(`GLB container version ${ver} is not supported (only glTF 2.0)`);
    for (let p = 12; p + 8 <= b.length;) {
      const len = dv.getUint32(p, true), type = dv.getUint32(p + 4, true);
      if (p + 8 + len > b.length) throw new Error('GLB chunk extends past the end of the file (truncated)');
      if (type === 0x4e4f534a) json = JSON.parse(text(b.subarray(p + 8, p + 8 + len)));
      else if (type === 0x004e4942 && !bin) bin = b.subarray(p + 8, p + 8 + len);
      p += 8 + len + ((4 - (len % 4)) % 4);
    }
    if (!json) throw new Error('GLB file has no JSON chunk');
  } else json = JSON.parse(text(b));
  if (!json || typeof json !== 'object' || !json.asset) throw new Error('not a glTF document (no asset object)');
  if (json.asset.version && !String(json.asset.version).startsWith('2')) throw new Error(`glTF version ${json.asset.version} is not supported (only 2.0)`);
  const comp = ctx.opts.companions || [];
  const buffers = (json.buffers || []).map((bf, i) => {
    if (!bf.uri) return i === 0 ? bin : null;
    if (bf.uri.startsWith('data:')) { const c = bf.uri.indexOf(','); return c < 0 ? null : base64(bf.uri.slice(c + 1)); }
    const want = decodeURIComponent(bf.uri).split('/').pop().toLowerCase(), hit = comp.find((f) => String(f.name).split(/[\\/]/).pop().toLowerCase() === want);
    if (!hit) ctx.warn(`External buffer "${bf.uri}" was not supplied; pass it in opts.companions or export a single-file .glb.`);
    return hit ? hit.bytes : null;
  });
  const accessor = (ai) => {
    const a = json.accessors?.[ai]; if (!a) return null;
    if (a.sparse) { ctx.warn('A sparse accessor was skipped (not supported).'); return null; }
    const bv = json.bufferViews?.[a.bufferView], buf = bv ? buffers[bv.buffer] : null, nc = GL_COMP[a.type], cs = GL_SIZE[a.componentType];
    if (!buf || !nc || !cs) return null;
    const count = guard(a.count, 'glTF accessor', LIMITS.elems), stride = bv.byteStride || nc * cs, off = (bv.byteOffset || 0) + (a.byteOffset || 0);
    if (count && off + stride * (count - 1) + nc * cs > buf.length) throw new Error('glTF accessor reaches past the end of its buffer');
    const dv = view(buf), out = a.componentType === 5126 ? new Float64Array(count * nc) : new Float64Array(count * nc);
    const norm = a.normalized ? { 5120: 127, 5121: 255, 5122: 32767, 5123: 65535 }[a.componentType] || 1 : 1;
    for (let i = 0; i < count; i++) for (let k = 0; k < nc; k++) {
      const p = off + i * stride + k * cs; let v;
      switch (a.componentType) { case 5120: v = dv.getInt8(p); break; case 5121: v = dv.getUint8(p); break; case 5122: v = dv.getInt16(p, true); break; case 5123: v = dv.getUint16(p, true); break; case 5125: v = dv.getUint32(p, true); break; default: v = dv.getFloat32(p, true); }
      out[i * nc + k] = norm === 1 ? v : Math.max(v / norm, -1);
    }
    return out;
  };
  const M = new Mesh(); const lines = []; let skipped = 0, nonTri = 0;
  const addMesh = (mi, mat, label) => {
    const mesh = json.meshes?.[mi]; if (!mesh) return;
    const g = M.group('component', M.groups.length, label || mesh.name || `mesh ${mi}`), flip = M4.det3(mat) < 0;
    for (const pr of mesh.primitives || []) {
      if (pr.extensions && (pr.extensions.KHR_draco_mesh_compression || pr.extensions.EXT_meshopt_compression)) { skipped++; continue; }
      const P = accessor(pr.attributes?.POSITION); if (!P) { skipped++; continue; }
      const base = M.nNodes, nv = P.length / 3, mode = pr.mode ?? 4;
      for (let i = 0; i < nv; i++) { const q = M4.apply(mat, P[3 * i], P[3 * i + 1], P[3 * i + 2]); M.node(q[0], q[1], q[2]); }
      let I = pr.indices !== undefined ? accessor(pr.indices) : null;
      if (pr.indices !== undefined && !I) { skipped++; continue; }
      if (!I) { I = new Float64Array(nv); for (let i = 0; i < nv; i++) I[i] = i; }
      const T = (a, c, d) => { if (a < nv && c < nv && d < nv) M.elem('tri3', flip ? [base + a, base + d, base + c] : [base + a, base + c, base + d], g); };
      if (mode === 4) for (let i = 0; i + 2 < I.length; i += 3) T(I[i], I[i + 1], I[i + 2]);
      else if (mode === 5) for (let i = 0; i + 2 < I.length; i++) (i % 2 ? T(I[i + 1], I[i], I[i + 2]) : T(I[i], I[i + 1], I[i + 2]));
      else if (mode === 6) for (let i = 1; i + 1 < I.length; i++) T(I[0], I[i], I[i + 1]);
      else if (mode === 1) { for (let i = 0; i + 1 < I.length; i += 2) if (I[i] < nv && I[i + 1] < nv) lines.push(base + I[i], base + I[i + 1]); nonTri++; }
      else if (mode === 3 || mode === 2) { for (let i = 0; i + 1 < I.length; i++) if (I[i] < nv && I[i + 1] < nv) lines.push(base + I[i], base + I[i + 1]); nonTri++; }
      else nonTri++;
    }
  };
  const nodes = json.nodes || [], seen = new Set();
  const walk = (ni, parent, depth) => {
    const nd = nodes[ni]; if (!nd || depth > 64 || seen.has(ni)) return;
    seen.add(ni);
    const local = Array.isArray(nd.matrix) && nd.matrix.length === 16 ? nd.matrix : M4.trs(nd.translation, nd.rotation, nd.scale), mat = M4.mul(parent, local);
    if (nd.mesh !== undefined) addMesh(nd.mesh, mat, nd.name);
    for (const c of nd.children || []) walk(c, mat, depth + 1);
    seen.delete(ni);
  };
  const scene = json.scenes?.[json.scene ?? 0];
  if (scene?.nodes?.length) for (const r of scene.nodes) walk(r, M4.I(), 0);
  else if (nodes.length) { const kids = new Set(nodes.flatMap((nd) => nd.children || [])); nodes.forEach((_, i) => { if (!kids.has(i)) walk(i, M4.I(), 0); }); }
  else (json.meshes || []).forEach((_, i) => addMesh(i, M4.I(), null));
  if (skipped) ctx.warn(`${skipped} glTF primitive(s) were skipped (compressed with Draco/meshopt, sparse, or missing buffer data).`);
  return M.result({
    lines: Uint32Array.from(lines), kind: 'surface', units: { length: 'm', source: 'file' },
    meta: { container: bin || b[0] === 0x67 ? 'glb' : 'gltf', generator: json.asset.generator ?? null, version: json.asset.version ?? null, upAxis: 'Y', unitNote: 'glTF 2.0 defines all lengths in metres and +Y as up; use transform({ swapYZ: true }) for a Z-up frame.', nodes: nodes.length, meshes: (json.meshes || []).length, nonTrianglePrimitives: nonTri, extensionsUsed: json.extensionsUsed || [] },
  });
}

// ---------- 3MF ----------
// 3MF transform: 12 numbers m00 m01 m02 m10 … m32 of a row-vector affine matrix → column-major 4×4.
const m3mf = (s) => { const v = s ? nums(s) : []; return v.length === 12 && v.every(Number.isFinite) ? [v[0], v[1], v[2], 0, v[3], v[4], v[5], 0, v[6], v[7], v[8], 0, v[9], v[10], v[11], 1] : M4.I(); };
export async function read3MF(b, ctx) {
  let xml;
  if (b[0] === 0x50 && b[1] === 0x4b) {
    const entries = zipEntries(b), e = entries.find((x) => /^\/?3D\/[^/]*\.model$/i.test(x.name)) || entries.find((x) => /\.model$/i.test(x.name));
    if (!e) throw new Error('the ZIP container holds no 3D/*.model part');
    if (entries.filter((x) => /\.model$/i.test(x.name)).length > 1) ctx.warn('The package holds several model parts (production extension); only the root model part was read.');
    xml = text(await e.read());
  } else xml = text(b);
  const root = xmlFirst(xmlParse(xml), 'model');
  if (!root) throw new Error('no <model> element in the 3MF model part');
  const objects = new Map();
  for (const o of xmlAll(root, 'object')) {
    const mesh = xmlChild(o, 'mesh'), rec = { name: o.attrs.name || `object ${o.attrs.id}`, v: [], t: [], comps: [] };
    if (mesh) {
      for (const v of xmlChild(mesh, 'vertices')?.children || []) rec.v.push(+v.attrs.x, +v.attrs.y, +v.attrs.z);
      for (const t of xmlChild(mesh, 'triangles')?.children || []) rec.t.push(+t.attrs.v1, +t.attrs.v2, +t.attrs.v3);
    }
    for (const c of xmlChild(o, 'components')?.children || []) rec.comps.push({ id: c.attrs.objectid, m: m3mf(c.attrs.transform) });
    objects.set(o.attrs.id, rec);
  }
  const M = new Mesh(); let bad = 0;
  const emit = (id, mat, g, depth) => {
    const o = objects.get(id); if (!o || depth > 32) return;
    const base = M.nNodes, nv = o.v.length / 3, flip = M4.det3(mat) < 0;
    guard(base + nv, '3MF vertex', LIMITS.verts);
    for (let i = 0; i < nv; i++) { const q = M4.apply(mat, o.v[3 * i], o.v[3 * i + 1], o.v[3 * i + 2]); M.node(q[0], q[1], q[2]); }
    for (let i = 0; i < o.t.length; i += 3) { const a = o.t[i], c = o.t[i + 1], d = o.t[i + 2]; if (a < nv && c < nv && d < nv && a >= 0 && c >= 0 && d >= 0) M.elem('tri3', flip ? [base + a, base + d, base + c] : [base + a, base + c, base + d], g); else bad++; }
    for (const c of o.comps) emit(c.id, M4.mul(mat, c.m), g, depth + 1);
  };
  const items = xmlAll(root, 'item');
  if (items.length) for (const it of items) { const o = objects.get(it.attrs.objectid); if (o) emit(it.attrs.objectid, m3mf(it.attrs.transform), M.group('component', M.groups.length, o.name), 0); }
  else for (const [id, o] of objects) if (o.v.length) emit(id, M4.I(), M.group('component', M.groups.length, o.name), 0);
  if (bad) ctx.warn(`${bad} 3MF triangle(s) referenced vertices outside their object and were skipped.`);
  const unitAttr = root.attrs.unit ?? null, length = unitFromWord(unitAttr);
  if (unitAttr && !length) ctx.warn(`The 3MF unit "${unitAttr}" has no platform equivalent; assign the units manually.`);
  const metadata = {}; for (const md of xmlAll(root, 'metadata')) if (md.attrs.name) metadata[md.attrs.name] = md.text.trim();
  return M.result({ kind: 'surface', units: { length, source: length ? 'file' : null }, meta: { unit: unitAttr, unitNote: unitAttr ? null : 'No unit attribute; the 3MF specification default is millimetre, but the units are left for the user to confirm.', objects: objects.size, buildItems: items.length, metadata } });
}

// ---------- AMF ----------
export async function readAMF(b, ctx) {
  let xml;
  if (b[0] === 0x50 && b[1] === 0x4b) {
    const entries = zipEntries(b), e = entries.find((x) => /\.amf$/i.test(x.name)) || entries[0];
    if (!e) throw new Error('the zipped AMF holds no entries');
    xml = text(await e.read());
  } else xml = text(b);
  const root = xmlFirst(xmlParse(xml), 'amf');
  if (!root) throw new Error('no <amf> root element');
  const M = new Mesh(); let bad = 0;
  const val = (nd, k) => +(xmlChild(nd, k)?.text ?? NaN);
  for (const o of xmlAll(root, 'object')) for (const mesh of xmlAll(o, 'mesh')) {
    const base = M.nNodes;
    for (const v of xmlAll(mesh, 'coordinates')) M.node(val(v, 'x'), val(v, 'y'), val(v, 'z'));
    const nv = M.nNodes - base;
    xmlAll(mesh, 'volume').forEach((vol, k) => {
      const name = xmlAll(vol, 'metadata').find((m) => m.attrs.type === 'name')?.text.trim() || `object ${o.attrs.id ?? ''} volume ${k + 1}`;
      const g = M.group('solid', M.groups.length, name);
      for (const t of xmlAll(vol, 'triangle')) { const a = val(t, 'v1'), c = val(t, 'v2'), d = val(t, 'v3'); if (a >= 0 && a < nv && c >= 0 && c < nv && d >= 0 && d < nv) M.elem('tri3', [base + a, base + c, base + d], g); else bad++; }
    });
  }
  if (bad) ctx.warn(`${bad} AMF triangle(s) had missing or out-of-range vertex references and were skipped.`);
  if (xmlFirst(root, 'constellation')) ctx.warn('The AMF constellation (object instances/placements) is not applied; objects are shown in their own frames.');
  const unitAttr = root.attrs.unit ?? null, length = unitFromWord(unitAttr);
  if (unitAttr && !length) ctx.warn(`The AMF unit "${unitAttr}" has no platform equivalent; assign the units manually.`);
  return M.result({ kind: 'surface', units: { length, source: length ? 'file' : null }, meta: { unit: unitAttr, version: root.attrs.version ?? null } });
}

// ---------- COLLADA ----------
export function readDAE(b, ctx) {
  const root = xmlFirst(xmlParse(text(b)), 'COLLADA');
  if (!root) throw new Error('no <COLLADA> root element');
  const M = new Mesh(); let bad = 0;
  for (const geo of xmlAll(root, 'geometry')) {
    const mesh = xmlChild(geo, 'mesh'); if (!mesh) continue;
    const sources = new Map();
    for (const s of mesh.children.filter((c) => c.name === 'source')) { const fa = xmlChild(s, 'float_array'); if (fa) sources.set('#' + s.attrs.id, { data: nums(fa.text), stride: +(xmlFirst(s, 'accessor')?.attrs.stride ?? 3) }); }
    const vertsEl = xmlChild(mesh, 'vertices'), posRef = vertsEl?.children.find((c) => c.name === 'input' && c.attrs.semantic === 'POSITION')?.attrs.source;
    const src = sources.get(posRef); if (!src) continue;
    const base = M.nNodes, nv = Math.floor(src.data.length / src.stride), g = M.group('component', M.groups.length, geo.attrs.name || geo.attrs.id || 'geometry');
    for (let i = 0; i < nv; i++) M.node(src.data[i * src.stride], src.data[i * src.stride + 1], src.stride > 2 ? src.data[i * src.stride + 2] : 0);
    const T = (a, c, d) => { if (a >= 0 && a < nv && c >= 0 && c < nv && d >= 0 && d < nv) M.elem('tri3', [base + a, base + c, base + d], g); else bad++; };
    for (const prim of mesh.children) {
      if (!['triangles', 'polylist', 'polygons'].includes(prim.name)) continue;
      const inputs = prim.children.filter((c) => c.name === 'input'), vin = inputs.find((c) => c.attrs.semantic === 'VERTEX'); if (!vin) continue;
      const step = Math.max(...inputs.map((c) => +(c.attrs.offset ?? 0))) + 1, off = +(vin.attrs.offset ?? 0);
      const poly = (p) => { const k = Math.floor(p.length / step); for (let j = 2; j < k; j++) T(p[off], p[(j - 1) * step + off], p[j * step + off]); };
      if (prim.name === 'polygons') { for (const p of prim.children.filter((c) => c.name === 'p')) poly(nums(p.text)); continue; }
      const p = nums(xmlChild(prim, 'p')?.text ?? '');
      if (prim.name === 'triangles') for (let i = 0; i + 3 * step <= p.length; i += 3 * step) T(p[i + off], p[i + step + off], p[i + 2 * step + off]);
      else { let o = 0; for (const vc of nums(xmlChild(prim, 'vcount')?.text ?? '')) { if (!(vc >= 0) || o + vc * step > p.length) break; poly(p.slice(o, o + vc * step)); o += vc * step; } }
    }
  }
  if (bad) ctx.warn(`${bad} COLLADA triangle(s) had out-of-range indices and were skipped.`);
  const unit = xmlFirst(root, 'unit'), meter = unit ? +unit.attrs.meter : NaN, up = xmlFirst(root, 'up_axis')?.text.trim() ?? null;
  const length = { 1: 'm', 0.001: 'mm', 0.01: 'cm', 0.0254: 'in', 0.3048: 'ft' }[meter] ?? null;
  if (xmlAll(root, 'node').some((nd) => nd.children.some((c) => ['matrix', 'translate', 'rotate', 'scale'].includes(c.name)))) ctx.warn('The COLLADA scene uses node transforms; they are NOT applied — geometry is shown in its local frame.');
  return M.result({ kind: 'surface', units: { length, source: length ? 'file' : null }, meta: { unitMeter: Number.isFinite(meter) ? meter : null, unitName: unit?.attrs.name ?? null, upAxis: up, authoringTool: xmlFirst(root, 'authoring_tool')?.text.trim() ?? null } });
}

// ---------- VRML 2.0 / X3D ----------
function faceSet(M, pts, idx, g) {
  const base = M.nNodes, nv = Math.floor(pts.length / 3); let bad = 0, poly = [];
  for (let i = 0; i < nv; i++) M.node(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]);
  const flush = () => { for (let j = 2; j < poly.length; j++) { const a = poly[0], c = poly[j - 1], d = poly[j]; if (a < nv && c < nv && d < nv) M.elem('tri3', [base + a, base + c, base + d], g); else bad++; } poly = []; };
  for (const v of idx) { if (v < 0) flush(); else poly.push(v); }
  flush();
  return bad;
}
export function readVRML(b, ctx) {
  const src = text(b).replace(/#[^\n]*/g, ''), M = new Mesh(); let bad = 0, p = 0, sets = 0;
  for (let guardN = 0; guardN < 1e6; guardN++) {
    const at = src.indexOf('IndexedFaceSet', p); if (at < 0) break;
    const open = src.indexOf('{', at); if (open < 0) break;
    let depth = 0, end = open;
    for (; end < src.length; end++) { const c = src.charCodeAt(end); if (c === 123) depth++; else if (c === 125 && --depth === 0) break; }
    const body = src.slice(open, end), pm = /\bpoint\s*\[([^\]]*)\]/.exec(body), im = /\bcoordIndex\s*\[([^\]]*)\]/.exec(body);
    if (pm && im) { bad += faceSet(M, nums(pm[1]), nums(im[1]), M.group('component', M.groups.length, `IndexedFaceSet ${++sets}`)); }
    p = end + 1;
  }
  if (/\b(Transform|PROTO|USE)\b/.test(src)) ctx.warn('The VRML scene uses Transform / PROTO / USE nodes; they are NOT evaluated — face sets are shown in their local frames.');
  if (!sets) ctx.warn('No IndexedFaceSet with inline coordinates was found.');
  if (bad) ctx.warn(`${bad} face(s) had out-of-range indices and were skipped.`);
  return M.result({ kind: 'surface', meta: { faceSets: sets } });
}
export function readX3D(b, ctx) {
  const root = xmlParse(text(b)), M = new Mesh(); let bad = 0, sets = 0;
  for (const tag of ['IndexedFaceSet', 'IndexedTriangleSet']) for (const fs of xmlAll(root, tag)) {
    const co = xmlChild(fs, 'Coordinate'); if (!co?.attrs.point) continue;
    let idx = nums(fs.attrs.coordIndex ?? fs.attrs.index ?? '');
    if (tag === 'IndexedTriangleSet') { const t = []; for (let i = 0; i + 2 < idx.length; i += 3) t.push(idx[i], idx[i + 1], idx[i + 2], -1); idx = t; }
    bad += faceSet(M, nums(co.attrs.point), idx, M.group('component', M.groups.length, fs.attrs.DEF || `${tag} ${++sets}`)); sets += fs.attrs.DEF ? 1 : 0;
  }
  if (xmlFirst(root, 'Transform') || xmlAll(root, 'Coordinate').some((c) => c.attrs.USE)) ctx.warn('The X3D scene uses Transform nodes or USE references; they are NOT evaluated — face sets are shown in their local frames.');
  if (!M.nNodes) ctx.warn('No IndexedFaceSet / IndexedTriangleSet with inline coordinates was found.');
  if (bad) ctx.warn(`${bad} face(s) had out-of-range indices and were skipped.`);
  return M.result({ kind: 'surface', meta: { faceSets: M.groups.length } });
}

// ---------- XYZ / delimited point tables ----------
export function readXYZ(b, ctx) {
  const sc = new Scanner(b), xyz = [], skipped = []; let cols = 0, two = 0, three = 0;
  for (let ln; (ln = sc.line()) !== null;) {
    const t = ln.trim(); if (!t) continue;
    const p = t.split(/[\s,;]+/), x = toNum(p[0]), y = toNum(p[1] ?? '');
    if (!(x === x && y === y) || p[1] === undefined || p[0] === '') { if (skipped.length < 5) skipped.push(t.slice(0, 80)); else skipped.count = (skipped.count || 5) + 1; continue; }
    const z = p.length > 2 ? toNum(p[2]) : 0;
    if (p.length > 2 && z === z) three++; else two++;
    xyz.push(x, y, z === z ? z : 0); cols = Math.max(cols, p.length);
  }
  if (!xyz.length) throw new Error('no numeric coordinate rows were found');
  guard(xyz.length / 3, 'point', LIMITS.verts);
  const planar = two > 0 && three === 0, n = xyz.length / 3, lines = [];
  if (planar) for (let i = 1; i < n; i++) lines.push(i - 1, i);
  if (two && three) ctx.warn('Rows have a mix of two and three coordinates; missing z values were read as 0.');
  if (planar) ctx.warn('Two-column table: read as planar (x, y) coordinates with z = 0 and joined in file order into a polyline.');
  if (cols > 3) ctx.warn(`Rows carry ${cols} columns; only the first three are read as coordinates.`);
  return { positions: Float64Array.from(xyz), lines: Uint32Array.from(lines), kind: 'pointcloud', meta: { columns: cols, planar, skippedHeaderLines: skipped.slice(0, 5), skippedLineCount: skipped.count || skipped.length } };
}

// ---------- PCD ----------
export function readPCD(b, ctx) {
  const sc = new Scanner(b), h = {}; let data = null;
  for (let ln, k = 0; k < 200 && (ln = sc.line()) !== null; k++) {
    const t = ln.trim(); if (!t || t[0] === '#') continue;
    const p = t.split(/\s+/); h[p[0].toUpperCase()] = p.slice(1);
    if (p[0].toUpperCase() === 'DATA') { data = (p[1] || '').toLowerCase(); break; }
  }
  if (!data || !h.FIELDS) throw new Error('PCD header is incomplete (FIELDS / DATA missing)');
  const fields = h.FIELDS, size = (h.SIZE || fields.map(() => 4)).map(Number), type = h.TYPE || fields.map(() => 'F'), cnt = (h.COUNT || fields.map(() => 1)).map(Number);
  const ix = fields.indexOf('x'), iy = fields.indexOf('y'), iz = fields.indexOf('z');
  if (ix < 0 || iy < 0 || iz < 0) throw new Error('PCD file has no x/y/z fields');
  const n = guard(parseInt(h.POINTS?.[0] ?? (h.WIDTH?.[0] * (h.HEIGHT?.[0] ?? 1)), 10), 'PCD point', LIMITS.verts);
  const meta = { version: h.VERSION?.[0] ?? null, fields, width: +h.WIDTH?.[0] || null, height: +h.HEIGHT?.[0] || null, viewpoint: h.VIEWPOINT?.map(Number) ?? null, encoding: data, declaredPoints: n };
  if (data === 'binary_compressed') throw Object.assign(new Error('PCD binary_compressed (LZF) data is not decoded; re-save as ascii or binary'), { meta });
  const xyz = []; let dropped = 0;
  const add = (x, y, z) => { if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) xyz.push(x, y, z); else dropped++; };
  if (data === 'ascii') {
    const col = []; let c = 0; for (let f = 0; f < fields.length; f++) { col.push(c); c += cnt[f]; }
    for (let i = 0, ln; i < n && (ln = sc.line()) !== null;) { const p = ln.trim().split(/\s+/); if (p.length < c) continue; add(toNum(p[col[ix]]), toNum(p[col[iy]]), toNum(p[col[iz]])); i++; }
  } else if (data === 'binary') {
    const off = []; let rec = 0; for (let f = 0; f < fields.length; f++) { off.push(rec); rec += size[f] * cnt[f]; }
    if (!(rec > 0) || sc.p + rec * n > b.length) throw new Error('PCD binary data is shorter than the declared point count (truncated)');
    const dv = view(b), get = (f, p) => (type[f] === 'F' ? (size[f] === 8 ? dv.getFloat64(p, true) : dv.getFloat32(p, true)) : type[f] === 'U' ? (size[f] === 1 ? dv.getUint8(p) : size[f] === 2 ? dv.getUint16(p, true) : dv.getUint32(p, true)) : size[f] === 1 ? dv.getInt8(p) : size[f] === 2 ? dv.getInt16(p, true) : dv.getInt32(p, true));
    for (let i = 0, p = sc.p; i < n; i++, p += rec) add(get(ix, p + off[ix]), get(iy, p + off[iy]), get(iz, p + off[iz]));
  } else throw new Error(`unknown PCD DATA encoding "${data}"`);
  if (dropped) ctx.warn(`${dropped} non-finite point(s) (invalid returns) were dropped.`);
  if (xyz.length / 3 + dropped < n) ctx.warn(`The file declares ${n} points but only ${xyz.length / 3 + dropped} were present.`);
  return { positions: Float64Array.from(xyz), kind: 'pointcloud', meta };
}

// ---------- LAS 1.x ----------
/** Public header block shared by LAS and LAZ. */
export function lasHeader(b) {
  if (b.length < 227) throw new Error('LAS header is truncated');
  const dv = view(b), major = b[24], minor = b[25], fmt = b[104];
  const h = {
    version: `${major}.${minor}`, systemIdentifier: str(b, 26, 58).replace(/\0+$/, '').trim(), generatingSoftware: str(b, 58, 90).replace(/\0+$/, '').trim(),
    headerSize: dv.getUint16(94, true), offsetToPoints: dv.getUint32(96, true), vlrCount: dv.getUint32(100, true), pointFormat: fmt & 0x3f, compressed: (fmt & 0xc0) !== 0,
    recordLength: dv.getUint16(105, true), pointCount: dv.getUint32(107, true),
    scale: [dv.getFloat64(131, true), dv.getFloat64(139, true), dv.getFloat64(147, true)], offset: [dv.getFloat64(155, true), dv.getFloat64(163, true), dv.getFloat64(171, true)],
    bbox: { min: [dv.getFloat64(187, true), dv.getFloat64(203, true), dv.getFloat64(219, true)], max: [dv.getFloat64(179, true), dv.getFloat64(195, true), dv.getFloat64(211, true)] },
  };
  if (major === 1 && minor >= 4 && b.length >= 255) { const c = Number(dv.getBigUint64(247, true)); if (c > 0) h.pointCount = c; }
  return h;
}
export function readLAS(b, ctx) {
  const h = lasHeader(b);
  if (h.compressed) throw new Error('point records are LASzip-compressed (LAZ); decompress to .las first');
  const rec = h.recordLength, avail = rec >= 12 && h.offsetToPoints <= b.length ? Math.floor((b.length - h.offsetToPoints) / rec) : 0;
  if (!avail) throw new Error('LAS header points at no readable point records');
  let n = h.pointCount;
  if (n > avail) { ctx.warn(`The header declares ${n} points but the file holds only ${avail}; the file is truncated.`); n = avail; }
  const maxPoints = Math.max(1, Math.floor(ctx.opts.maxPoints ?? 200000)), stride = Math.max(1, Math.ceil(n / maxPoints)), m = Math.ceil(n / stride);
  const dv = view(b), pos = new Float64Array(3 * m);
  for (let i = 0, k = 0; i < n; i += stride, k++) {
    const p = h.offsetToPoints + i * rec;
    pos[3 * k] = dv.getInt32(p, true) * h.scale[0] + h.offset[0]; pos[3 * k + 1] = dv.getInt32(p + 4, true) * h.scale[1] + h.offset[1]; pos[3 * k + 2] = dv.getInt32(p + 8, true) * h.scale[2] + h.offset[2];
  }
  if (stride > 1) ctx.warn(`Point cloud subsampled: every ${stride}th of ${n} points was read (${m} points; raise opts.maxPoints to read more).`);
  return { positions: pos, kind: 'pointcloud', meta: { ...h, pointsRead: m, stride, unitNote: 'LAS units are defined by the coordinate reference system records, which are not interpreted; confirm the units.' } };
}

// ---------- platform JSON (written by exportModel(model, 'json')) ----------
export function readJSON(b) {
  const j = JSON.parse(text(b));
  if (!j || j.aerosuiteGeometry !== 1 || !Array.isArray(j.positions)) throw new Error('not a platform geometry JSON document');
  return {
    positions: Float64Array.from(j.positions), ...(j.triangles?.length ? { triangles: Uint32Array.from(j.triangles) } : {}), lines: Uint32Array.from(j.lines || []),
    elements: (j.elements || []).map((e) => ({ type: e.type, nodesPer: e.nodesPer, count: e.count, conn: Uint32Array.from(e.conn || []), ...(e.group ? { group: Int32Array.from(e.group) } : {}) })),
    groups: j.groups || [], kind: j.kind, units: j.units?.length ? j.units : undefined, meta: { ...(j.meta || {}), sourceFormat: j.format ?? null, sourceName: j.name ?? null, priorLog: j.log || [] }, warnings: (j.warnings || []).map((w) => `(from source) ${w}`),
  };
}
