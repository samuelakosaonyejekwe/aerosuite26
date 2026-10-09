// Geometry / mesh interoperability tests: `node tests/geometry.mjs`
// Every fixture is built in memory. Exits non-zero when any check fails.

import zlib from 'node:zlib';
import * as G from '../js/core/geometry/index.js';
import { loadH5 } from '../js/core/geometry/parsers-wasm.js';

let pass = 0, fail = 0, section = '';
const failures = [];
const ok = (name, cond, detail = '') => { if (cond) pass++; else { fail++; failures.push(`[${section}] ${name}${detail ? ' — ' + detail : ''}`); } };
const near = (name, a, b, tol = 1e-9) => ok(name, Number.isFinite(a) && Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `got ${a}, expected ${b} (tol ${tol})`);
const sec = (s) => { section = s; };
const enc = (s) => new TextEncoder().encode(s);
const cat = (...parts) => { const n = parts.reduce((s, p) => s + p.length, 0), o = new Uint8Array(n); let k = 0; for (const p of parts) { o.set(p, k); k += p.length; } return o; };
const imp = (name, data, opts) => G.importFile(name, typeof data === 'string' ? enc(data) : data, opts);

// ---------- fixtures ----------
const CV = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
const CT = [[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]];
const CQ = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]];
// Kuhn split of the unit cube into 6 positively oriented tetrahedra
const TETS = (() => { const id = (p) => CV.findIndex((v) => v[0] === p[0] && v[1] === p[1] && v[2] === p[2]), out = []; for (const [a, b, c] of [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]) { const p1 = [0, 0, 0]; p1[a] = 1; const p2 = p1.slice(); p2[b] = 1; const t = [0, id(p1), id(p2), 6], even = [[0, 1, 2], [1, 2, 0], [2, 0, 1]].some((q) => q[0] === a && q[1] === b && q[2] === c); out.push(even ? t : [t[0], t[2], t[1], t[3]]); } return out; })();
// boundary triangles of the tetrahedralised cube (faces used by exactly one tet, outward)
const BT = (() => { const m = new Map(); for (const t of TETS) for (const f of [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]) { const n = f.map((i) => t[i]), k = n.slice().sort().join(); m.set(k, m.has(k) ? null : n); } return [...m.values()].filter(Boolean); })();
const mk = (verts, tris, extra = {}) => ({ name: 'fixture', format: 'test', formatName: 'test', support: 'native', kind: 'surface', positions: Float64Array.from(verts.flat()), triangles: Uint32Array.from(tris.flat()), lines: new Uint32Array(0), elements: [], groups: [], units: { length: null, source: null }, meta: {}, log: [], warnings: [], ...extra });
const cube = () => mk(CV, CT);
const tetCube = () => { const m = mk(CV, BT, { kind: 'volume-mesh', elements: [{ type: 'tet4', nodesPer: 4, count: 6, conn: Uint32Array.from(TETS.flat()), group: new Int32Array(6) }], groups: [{ id: 0, name: 'solid body', kind: 'zone', count: 6 }] }); return m; };
function uvSphere(r, nu, nv) {
  const V = [[0, 0, r]], T = [];
  for (let i = 1; i < nv; i++) for (let j = 0; j < nu; j++) { const th = (Math.PI * i) / nv, ph = (2 * Math.PI * j) / nu; V.push([r * Math.sin(th) * Math.cos(ph), r * Math.sin(th) * Math.sin(ph), r * Math.cos(th)]); }
  V.push([0, 0, -r]); const S = V.length - 1, id = (i, j) => 1 + (i - 1) * nu + (j % nu);
  for (let j = 0; j < nu; j++) { T.push([0, id(1, j), id(1, j + 1)]); T.push([S, id(nv - 1, j + 1), id(nv - 1, j)]); }
  for (let i = 1; i < nv - 1; i++) for (let j = 0; j < nu; j++) { T.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]); T.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]); }
  return mk(V, T);
}
function cylinder(r, L, n) {
  const V = [[0, 0, 0], [0, 0, L]], T = [];
  for (let j = 0; j < n; j++) { const a = (2 * Math.PI * j) / n; V.push([r * Math.cos(a), r * Math.sin(a), 0], [r * Math.cos(a), r * Math.sin(a), L]); }
  for (let j = 0; j < n; j++) { const b0 = 2 + 2 * j, t0 = b0 + 1, b1 = 2 + 2 * ((j + 1) % n), t1 = b1 + 1; T.push([0, b1, b0], [1, t0, t1], [b0, b1, t1], [b0, t1, t0]); }
  return mk(V, T);
}
function naca(t, n = 120, camber = 0) {
  const up = [], lo = [];
  for (let i = 0; i <= n; i++) { const x = 0.5 * (1 - Math.cos((Math.PI * i) / n)), y = 5 * t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1036 * x ** 4), yc = camber * 4 * x * (1 - x); up.push([x, yc + y]); lo.push([x, yc - y]); }
  return [...up.reverse(), ...lo.slice(1, -1)];        // TE → upper → LE → lower → (TE)
}
const cubeCheck = async (label, model, { nv = null, vol = true, area = 6 } = {}) => {
  const a = G.analyse(model);
  if (nv !== null) ok(`${label}: vertex count`, a.nVerts === nv, `got ${a.nVerts}`);
  ok(`${label}: bbox`, a.bbox.size.every((s) => Math.abs(s - 1) < 1e-6), JSON.stringify(a.bbox.size));
  near(`${label}: area`, a.area, area, 1e-6);
  ok(`${label}: watertight`, a.watertight === true, `boundary ${a.boundaryEdges}, non-manifold ${a.nonManifoldEdges}`);
  if (vol) near(`${label}: volume`, a.volume, 1, 1e-6);
  ok(`${label}: log has provenance`, model.log.length >= 3 && model.log[0].step === 'source', JSON.stringify(model.log[0]));
  return a;
};
const fmtId = (label, m, id, support) => { ok(`${label}: format ${id}`, m.format === id, `got ${m.format}`); if (support) ok(`${label}: support ${support}`, m.support === support, m.support); };

// ---------- ZIP builder (test fixture only) ----------
const CRC = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return (b) => { let c = 0xffffffff; for (const v of b) c = t[(c ^ v) & 255] ^ (c >>> 8); return (c ^ 0xffffffff) >>> 0; }; })();
function zip(files) {
  const locals = [], centrals = []; let off = 0;
  for (const [name, data, store] of files) {
    const nm = enc(name), comp = store ? data : new Uint8Array(zlib.deflateRawSync(data)), h = new DataView(new ArrayBuffer(30));
    h.setUint32(0, 0x04034b50, true); h.setUint16(4, 20, true); h.setUint16(8, store ? 0 : 8, true); h.setUint32(14, CRC(data), true); h.setUint32(18, comp.length, true); h.setUint32(22, data.length, true); h.setUint16(26, nm.length, true);
    const c = new DataView(new ArrayBuffer(46));
    c.setUint32(0, 0x02014b50, true); c.setUint16(4, 20, true); c.setUint16(6, 20, true); c.setUint16(10, store ? 0 : 8, true); c.setUint32(16, CRC(data), true); c.setUint32(20, comp.length, true); c.setUint32(24, data.length, true); c.setUint16(28, nm.length, true); c.setUint32(42, off, true);
    locals.push(new Uint8Array(h.buffer), nm, comp); centrals.push(new Uint8Array(c.buffer), nm); off += 30 + nm.length + comp.length;
  }
  const cd = cat(...centrals), e = new DataView(new ArrayBuffer(22));
  e.setUint32(0, 0x06054b50, true); e.setUint16(8, files.length, true); e.setUint16(10, files.length, true); e.setUint32(12, cd.length, true); e.setUint32(16, off, true);
  return cat(...locals, cd, new Uint8Array(e.buffer));
}

const truncatable = [];      // [name, bytes] of every fixture, reused by the malformed-input sweep
const keep = (name, data) => { const b = typeof data === 'string' ? enc(data) : data; truncatable.push([name, b]); return b; };

// =====================================================================================
sec('formats');
{
  ok('FORMATS is a non-empty array', Array.isArray(G.FORMATS) && G.FORMATS.length > 30);
  const ids = new Set();
  for (const f of G.FORMATS) {
    ok(`profile ${f.id} complete`, !!(f.id && f.name && f.ext?.length && f.category && f.support && f.reads && f.preserves && f.limits), JSON.stringify(f).slice(0, 120));
    ok(`profile ${f.id} enums`, ['cad', 'tessellated', 'mesh', 'pointcloud', 'geospatial'].includes(f.category) && ['native', 'partial', 'metadata'].includes(f.support));
    if (f.support !== 'native') ok(`profile ${f.id} has a conversion pathway`, f.pathway.length > 20);
    ok(`profile ${f.id} unique id`, !ids.has(f.id)); ids.add(f.id);
  }
  for (const e of ['step', 'stp', 'iges', 'igs', 'x_t', 'x_b', 'sat', 'sab', 'jt', 'stl', 'obj', 'ply', '3mf', 'amf', 'gltf', 'glb', 'dae', 'wrl', 'x3d', 'off', 'vtk', 'vtp', 'las', 'laz', 'e57', 'pcd', 'xyz', 'cgns', 'msh', 'vtu', 'vts', 'vtr', 'vti', 'vtm', 'e', 'exo', 'med', 'unv', 'bdf', 'nas', 'dat', 'inp', 'cdb', 'h5', 'su2', 'plt', 'szplt', 'p3d', 'q', 'tif', 'catpart', 'sldprt'])
    ok(`extension .${e} has a profile`, G.FORMATS.some((f) => f.ext.includes(e)));
}

// ---------- native readers ----------
sec('STL');
{
  const ascii = keep('cube.stl', G.exportModel(cube(), 'stl'));
  let m = await imp('cube.stl', ascii); fmtId('ascii', m, 'stl', 'native');
  await cubeCheck('ascii', m, { nv: 36 }); ok('ascii: units not guessed', m.units.length === null); ok('ascii: 12 triangles', m.triangles.length === 36); ok('ascii: kind surface', m.kind === 'surface');
  ok('ascii: duplicate vertices reported, not repaired', G.analyse(m).duplicateVerts === 28);
  const bin = (header) => { const b = new Uint8Array(84 + 50 * 12), dv = new DataView(b.buffer); b.set(enc(header).subarray(0, 80)); dv.setUint32(80, 12, true); CT.forEach((t, i) => t.forEach((v, k) => CV[v].forEach((c, j) => dv.setFloat32(84 + 50 * i + 12 + 12 * k + 4 * j, c, true)))); return b; };
  m = await imp('cube.stl', keep('cubeb.stl', bin('binary cube'))); await cubeCheck('binary', m, { nv: 36 }); ok('binary: encoding', m.meta.encoding === 'binary');
  const trap = bin('solid cube exported by some CAD system');
  ok('solid-prefixed binary detected as stl', G.detectFormat('x.stl', trap).id === 'stl');
  m = await imp('trap.stl', trap); await cubeCheck('solid-prefixed binary', m, { nv: 36 }); ok('solid-prefixed binary: flagged', m.meta.solidPrefixedBinary === true && m.meta.encoding === 'binary');
  m = await imp('two.stl', 'solid a\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid a\nsolid b\nfacet normal 0 0 1\nouter loop\nvertex 0 0 1\nvertex 1 0 1\nvertex 0 1 1\nendloop\nendfacet\nendsolid b\n');
  ok('two solids → two groups', m.groups.length === 2 && m.groups[0].name === 'a' && m.groups[1].kind === 'solid' && m.elements[0]?.group?.[1] === 1, JSON.stringify(m.groups));
}

sec('OBJ');
{
  const src = '# cube\no Box\n' + CV.map((v) => `v ${v.join(' ')}`).join('\n') + '\ng sides\n' + CQ.slice(0, 3).map((q) => `f ${q.map((v) => `${v + 1}/1/1`).join(' ')}`).join('\n') + '\ng rest\nusemtl steel\n' + CQ.slice(3).map((q) => `f ${q.map((v) => v - 8).join(' ')}`).join('\n') + '\n';
  const m = await imp('cube.obj', keep('cube.obj', src)); fmtId('obj', m, 'obj', 'native');
  await cubeCheck('quads + negative indices', m, { nv: 8 });
  ok('groups from o/g', m.groups.map((g) => g.name).join() === 'Box,sides,rest' && m.groups[1].count === 6, JSON.stringify(m.groups));
  ok('quads triangulated and counted', m.triangles.length === 36 && m.meta.quadsTriangulated === 6); ok('material recorded', m.meta.materials[0] === 'steel');
  const rt = await imp('rt.obj', G.exportModel(cube(), 'obj')); await cubeCheck('export round trip', rt, { nv: 8 });
}

sec('PLY');
{
  const head = (fmt) => `ply\nformat ${fmt} 1.0\ncomment made by test\nelement vertex 8\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nelement face 6\nproperty list uchar int vertex_indices\nend_header\n`;
  let m = await imp('cube.ply', keep('cube.ply', head('ascii') + CV.map((v) => v.join(' ') + ' 255').join('\n') + '\n' + CQ.map((q) => '4 ' + q.join(' ')).join('\n') + '\n'));
  fmtId('ply', m, 'ply', 'native'); await cubeCheck('ascii', m, { nv: 8 }); ok('ascii: comment kept', m.meta.comments[0] === 'made by test');
  for (const le of [true, false]) {
    const body = new Uint8Array(8 * 13 + 6 * 17), dv = new DataView(body.buffer); let p = 0;
    for (const v of CV) { for (const c of v) { dv.setFloat32(p, c, le); p += 4; } dv.setUint8(p++, 200); }
    for (const q of CQ) { dv.setUint8(p++, 4); for (const v of q) { dv.setInt32(p, v, le); p += 4; } }
    m = await imp('cube.ply', keep(`cube_${le ? 'le' : 'be'}.ply`, cat(enc(head(le ? 'binary_little_endian' : 'binary_big_endian')), body)));
    await cubeCheck(le ? 'binary LE' : 'binary BE', m, { nv: 8 });
  }
  m = await imp('pts.ply', 'ply\nformat ascii 1.0\nelement vertex 3\nproperty double x\nproperty double y\nproperty double z\nend_header\n0 0 0\n1 0 0\n0 2 0\n');
  ok('point-cloud PLY', m.kind === 'pointcloud' && m.positions.length === 9 && G.analyse(m).classification.label === 'point cloud / scan');
}

sec('OFF');
{
  const m = await imp('cube.off', keep('cube.off', 'OFF\n# comment\n8 6 12\n' + CV.map((v) => v.join(' ')).join('\n') + '\n' + CQ.map((q) => '4 ' + q.join(' ') + ' 255 0 0').join('\n') + '\n'));
  fmtId('off', m, 'off', 'native'); await cubeCheck('off', m, { nv: 8 });
}

sec('glTF');
{
  const bin = new Uint8Array(8 * 12 + 36 * 2), dv = new DataView(bin.buffer);
  CV.forEach((v, i) => v.forEach((c, k) => dv.setFloat32(12 * i + 4 * k, c, true))); CT.flat().forEach((v, i) => dv.setUint16(96 + 2 * i, v, true));
  const doc = (uri) => ({ asset: { version: '2.0', generator: 'test' }, scene: 0, scenes: [{ nodes: [0] }], nodes: [{ name: 'root', children: [1] }, { name: 'Wing box', mesh: 0, translation: [1, 2, 3] }], meshes: [{ name: 'cube', primitives: [{ attributes: { POSITION: 0 }, indices: 1, mode: 4 }] }], buffers: [uri ? { uri, byteLength: bin.length } : { byteLength: bin.length }], bufferViews: [{ buffer: 0, byteOffset: 0, byteLength: 96 }, { buffer: 0, byteOffset: 96, byteLength: 72 }], accessors: [{ bufferView: 0, componentType: 5126, count: 8, type: 'VEC3' }, { bufferView: 1, componentType: 5123, count: 36, type: 'SCALAR' }] });
  let json = enc(JSON.stringify(doc(null))); while (json.length % 4) json = cat(json, enc(' '));
  const hdr = new DataView(new ArrayBuffer(12)), c1 = new DataView(new ArrayBuffer(8)), c2 = new DataView(new ArrayBuffer(8));
  hdr.setUint32(0, 0x46546c67, true); hdr.setUint32(4, 2, true); hdr.setUint32(8, 12 + 8 + json.length + 8 + bin.length, true);
  c1.setUint32(0, json.length, true); c1.setUint32(4, 0x4e4f534a, true); c2.setUint32(0, bin.length, true); c2.setUint32(4, 0x004e4942, true);
  let m = await imp('cube.glb', keep('cube.glb', cat(new Uint8Array(hdr.buffer), new Uint8Array(c1.buffer), json, new Uint8Array(c2.buffer), bin)));
  fmtId('glb', m, 'gltf', 'native'); const a = await cubeCheck('glb', m, { nv: 8 });
  ok('glb: node transform applied', Math.abs(a.bbox.min[0] - 1) < 1e-6 && Math.abs(a.bbox.min[1] - 2) < 1e-6 && Math.abs(a.bbox.min[2] - 3) < 1e-6, JSON.stringify(a.bbox.min));
  ok('glb: node name → group', m.groups[0]?.name === 'Wing box' && m.groups[0].count === 12); ok('glb: metres by specification', m.units.length === 'm' && m.units.source === 'file');
  m = await imp('cube.gltf', keep('cube.gltf', JSON.stringify(doc('data:application/octet-stream;base64,' + Buffer.from(bin).toString('base64')))));
  await cubeCheck('gltf embedded base64', m, { nv: 8 });
  m = await imp('cube.gltf', JSON.stringify(doc('cube.bin')));
  ok('gltf external buffer missing → warning, no throw', m.warnings.some((w) => /cube\.bin/.test(w)));
  m = await imp('cube.gltf', JSON.stringify(doc('cube.bin')), { companions: [{ name: 'cube.bin', bytes: bin }] }); await cubeCheck('gltf companion buffer', m, { nv: 8 });
}

sec('3MF / AMF');
{
  const model = `<?xml version="1.0" encoding="UTF-8"?>\n<model unit="inch" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02"><metadata name="Title">Cube</metadata><resources><object id="1" name="Bracket" type="model"><mesh><vertices>${CV.map((v) => `<vertex x="${v[0]}" y="${v[1]}" z="${v[2]}"/>`).join('')}</vertices><triangles>${CT.map((t) => `<triangle v1="${t[0]}" v2="${t[1]}" v3="${t[2]}"/>`).join('')}</triangles></mesh></object></resources><build><item objectid="1" transform="1 0 0 0 1 0 0 0 1 5 0 0"/></build></model>`;
  const z = zip([['[Content_Types].xml', enc('<Types/>'), true], ['3D/3dmodel.model', enc(model), false]]);
  let m = await imp('cube.3mf', keep('cube.3mf', z)); fmtId('3mf', m, '3mf', 'native');
  const a = await cubeCheck('3mf (deflated)', m, { nv: 8 }); ok('3mf: unit attribute → units', m.units.length === 'in' && m.units.source === 'file'); ok('3mf: build transform applied', Math.abs(a.bbox.min[0] - 5) < 1e-9);
  ok('3mf: object name → group', m.groups[0]?.name === 'Bracket'); ok('3mf: metadata', m.meta.metadata.Title === 'Cube');
  ok('3mf detected without the extension', G.detectFormat('archive.zip', z).id === '3mf');
  m = await imp('nounit.3mf', zip([['3D/3dmodel.model', enc(model.replace(' unit="inch"', '')), true]])); ok('3mf without unit → null (not guessed)', m.units.length === null && m.triangles.length === 36);
  const amf = `<?xml version="1.0"?><amf unit="millimeter" version="1.1"><object id="0"><mesh><vertices>${CV.map((v) => `<vertex><coordinates><x>${v[0]}</x><y>${v[1]}</y><z>${v[2]}</z></coordinates></vertex>`).join('')}</vertices><volume materialid="1"><metadata type="name">Core</metadata>${CT.map((t) => `<triangle><v1>${t[0]}</v1><v2>${t[1]}</v2><v3>${t[2]}</v3></triangle>`).join('')}</volume></mesh></object></amf>`;
  m = await imp('cube.amf', keep('cube.amf', amf)); fmtId('amf', m, 'amf', 'native'); await cubeCheck('amf xml', m, { nv: 8 }); ok('amf: units', m.units.length === 'mm'); ok('amf: volume name', m.groups[0]?.name === 'Core');
  m = await imp('cube.amf', zip([['cube.amf', enc(amf), false]])); await cubeCheck('amf zipped', m, { nv: 8 });
}

sec('VTK legacy');
{
  let m = await imp('cube.vtk', keep('cube.vtk', G.exportModel(cube(), 'vtk'))); fmtId('polydata', m, 'vtk', 'native'); await cubeCheck('polydata', m, { nv: 8 }); ok('polydata dataset', m.meta.dataset === 'POLYDATA');
  m = await imp('tets.vtk', keep('tets.vtk', G.exportModel(tetCube(), 'vtk'))); await cubeCheck('unstructured tets', m, { nv: 8 });
  ok('tets: element block', m.elements.length === 1 && m.elements[0].type === 'tet4' && m.elements[0].count === 6 && m.kind === 'volume-mesh'); ok('tets: boundary faces only', m.triangles.length === 36, `${m.triangles.length / 3} triangles`);
  const hexSrc = '# vtk DataFile Version 3.0\nhex\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 8 float\n' + CV.map((v) => v.join(' ')).join('\n') + '\nCELLS 1 9\n8 0 1 2 3 4 5 6 7\nCELL_TYPES 1\n12\nCELL_DATA 1\nSCALARS pressure float 1\nLOOKUP_TABLE default\n1.5\n';
  m = await imp('hex.vtk', keep('hex.vtk', hexSrc)); await cubeCheck('hex cell', m, { nv: 8 }); ok('hex: type + cell data listed', m.elements[0].type === 'hex8' && m.meta.cellData[0] === 'pressure');
  m = await imp('v5.vtk', '# vtk DataFile Version 5.1\nnew\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 8 float\n' + CV.map((v) => v.join(' ')).join('\n') + '\nCELLS 2 8\nOFFSETS vtktypeint64\n0 8\nCONNECTIVITY vtktypeint64\n0 1 2 3 4 5 6 7\nCELL_TYPES 1\n12\n');
  await cubeCheck('5.1 offsets/connectivity layout', m, { nv: 8 });
  m = await imp('sg.vtk', '# vtk DataFile Version 3.0\nsg\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS 2 2 2\nPOINTS 8 float\n0 0 0 1 0 0 0 1 0 1 1 0 0 0 1 1 0 1 0 1 1 1 1 1\n'); await cubeCheck('structured grid', m, { nv: 8 });
  m = await imp('bin.vtk', '# vtk DataFile Version 3.0\nb\nBINARY\nDATASET POLYDATA\nPOINTS 1 float\n\0\0\0\0'); ok('binary legacy → metadata-only + warning', m.kind === 'metadata-only' && m.warnings.some((w) => /BINARY/.test(w)) && m.meta.dataset === 'POLYDATA');
}

sec('Gmsh');
{
  const tc = tetCube(); tc.elements.push({ type: 'tri3', nodesPer: 3, count: 12, conn: Uint32Array.from(BT.flat()), group: new Int32Array(12).fill(1) }); tc.groups.push({ id: 1, name: 'wall', kind: 'boundary', count: 12 });
  let m = await imp('cube.msh', keep('cube22.msh', G.exportModel(tc, 'msh'))); fmtId('2.2', m, 'gmsh', 'native'); await cubeCheck('2.2', m, { nv: 8 });
  ok('2.2: version', m.meta.version === '2.2'); ok('2.2: elements', m.elements.find((e) => e.type === 'tet4')?.count === 6 && m.elements.find((e) => e.type === 'tri3')?.count === 12);
  const gz = m.groups.find((g) => g.name === 'solid body'), gw = m.groups.find((g) => g.name === 'wall');
  ok('2.2: physical names → groups', gz?.kind === 'zone' && gz.count === 6 && gw?.kind === 'boundary' && gw.count === 12, JSON.stringify(m.groups));
  ok('2.2: boundary faces not duplicated by surface elements', m.triangles.length === 36);
  const v41 = `$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$PhysicalNames\n2\n2 7 "skin"\n3 9 "fluid volume"\n$EndPhysicalNames\n$Entities\n0 0 1 1\n1 0 0 0 1 1 1 1 7 0\n1 0 0 0 1 1 1 1 9 1 1\n$EndEntities\n$Nodes\n1 8 1 8\n3 1 0 8\n${[1, 2, 3, 4, 5, 6, 7, 8].join('\n')}\n${CV.map((v) => v.join(' ')).join('\n')}\n$EndNodes\n$Elements\n2 18 1 18\n2 1 2 12\n${BT.map((t, i) => `${i + 1} ${t.map((v) => v + 1).join(' ')}`).join('\n')}\n3 1 4 6\n${TETS.map((t, i) => `${i + 13} ${t.map((v) => v + 1).join(' ')}`).join('\n')}\n$EndElements\n`;
  m = await imp('cube41.msh', keep('cube41.msh', v41)); await cubeCheck('4.1', m, { nv: 8 }); ok('4.1: version', m.meta.version === '4.1');
  ok('4.1: entity → physical groups', m.groups.find((g) => g.name === 'fluid volume')?.count === 6 && m.groups.find((g) => g.name === 'skin')?.kind === 'boundary', JSON.stringify(m.groups));
  m = await imp('b.msh', '$MeshFormat\n2.2 1 8\n\x01\0\0\0\n$EndMeshFormat\n'); ok('binary msh → metadata-only', m.kind === 'metadata-only' && m.warnings.some((w) => /binary/i.test(w)));
}

sec('SU2');
{
  let m = await imp('cube.su2', keep('cube.su2', G.exportModel(tetCube(), 'su2'))); fmtId('su2', m, 'su2', 'native'); await cubeCheck('export round trip', m, { nv: 8 });
  ok('marker from export', m.groups.some((g) => g.kind === 'boundary' && g.name === 'boundary' && g.count === 12), JSON.stringify(m.groups));
  const src = `% hand-written\nNDIME= 3\nNELEM= 6\n${TETS.map((t, i) => `10 ${t.join(' ')} ${i}`).join('\n')}\nNPOIN= 8\n${CV.map((v, i) => `${v.join(' ')} ${i}`).join('\n')}\nNMARK= 2\nMARKER_TAG= wall\nMARKER_ELEMS= 10\n${BT.slice(0, 10).map((t) => `5 ${t.join(' ')}`).join('\n')}\nMARKER_TAG= farfield\nMARKER_ELEMS= 2\n${BT.slice(10).map((t) => `5 ${t.join(' ')}`).join('\n')}\n`;
  m = await imp('hand.su2', src); await cubeCheck('handcrafted', m, { nv: 8 });
  ok('MARKER tags → boundary groups', m.groups.length === 2 && m.groups[0].name === 'wall' && m.groups[0].count === 10 && m.groups[1].name === 'farfield' && m.groups[1].count === 2);
  ok('classified as CFD volume mesh', G.analyse(m).classification.label === 'CFD volume mesh');
}

sec('Nastran');
{
  const free = `SOL 101\nCEND\nTITLE = cube\nBEGIN BULK\n$ comment\n${CV.map((v, i) => `GRID,${i + 1},,${v[0]}.,${v[1]}.,${v[2]}.`).join('\n')}\n${TETS.map((t, i) => `CTETRA,${i + 1},5,${t.map((v) => v + 1).join(',')}`).join('\n')}\nPSOLID,5,3\nMAT1,3,7.0+10,,0.33,2.7-9\nENDDATA\n`;
  let m = await imp('cube.bdf', keep('cube_free.bdf', free)); fmtId('free', m, 'nastran', 'native'); await cubeCheck('free field', m, { nv: 8 });
  ok('free: kind structural', m.kind === 'structural-mesh' && m.elements[0].type === 'tet4' && m.elements[0].count === 6);
  ok('free: property group', m.groups.find((g) => g.kind === 'property')?.count === 6 && /PSOLID 5/.test(m.groups.find((g) => g.kind === 'property').name), JSON.stringify(m.groups));
  ok('free: material census', m.meta.materials[0]?.E === 7e10 && Math.abs(m.meta.materials[0].rho - 2.7e-9) < 1e-20 && m.groups.find((g) => g.kind === 'material')?.count === 6, JSON.stringify(m.meta.materials));
  ok('free: card census', m.meta.cards.GRID === 8 && m.meta.cards.CTETRA === 6); ok('free: units null', m.units.length === null);
  ok('free: classified structural', G.analyse(m).classification.label === 'structural mesh');
  const F = (...f) => f.map((s) => String(s).padEnd(8).slice(0, 8)).join('').trimEnd();
  const fixed = `BEGIN BULK\n${CV.map((v, i) => F('GRID', i + 1, '', v[0] + '.0', v[1] + '.0', v[2] + '.0')).join('\n')}\n${F('CHEXA', 1, 2, 1, 2, 3, 4, 5, 6, '+C1')}\n${F('+C1', 7, 8)}\n${F('CQUAD4', 2, 9, 1, 2, 3, 4)}\n${F('CBAR', 3, 11, 1, 7, '0.', '0.', '1.')}\n${F('PSHELL', 9, 3, '0.002')}\n${F('PSOLID', 2, 3)}\n${F('MAT1', 3, '7.+10', '', '.33')}\nENDDATA\n`;
  m = await imp('cube.nas', keep('cube_fixed.nas', fixed)); await cubeCheck('fixed field + continuation', m, { nv: 8 });
  ok('fixed: hex + quad + bar', m.elements.map((e) => `${e.type}:${e.count}`).sort().join() === 'hex8:1,line2:1,quad4:1', m.elements.map((e) => `${e.type}:${e.count}`).join());
  ok('fixed: PSHELL thickness in group name', m.groups.some((g) => /PSHELL 9 \(t = 0.002\)/.test(g.name)), JSON.stringify(m.groups.map((g) => g.name))); ok('fixed: bar as line', m.lines.length === 2);
  const large = `GRID*                  1               0             0.5            0.25*\n*                   0.75\nGRID,2,,1.,0.,0.\nGRID,3,7,0.,1.,0.\nCTRIA3,1,1,1,2,3\n`;
  m = await imp('l.bdf', large); ok('large field GRID*', Math.abs(m.positions[0] - 0.5) < 1e-12 && Math.abs(m.positions[1] - 0.25) < 1e-12 && Math.abs(m.positions[2] - 0.75) < 1e-12, Array.from(m.positions.slice(0, 3)).join());
  ok('local coordinate system reported, not applied', m.meta.gridsInLocalSystems === 1 && m.warnings.some((w) => /CP/.test(w)));
}

sec('Abaqus');
{
  const src = `*HEADING\ncube model\n** comment\n*PART, NAME=Block\n*NODE\n${CV.map((v, i) => `${i + 1}, ${v[0]}., ${v[1]}., ${v[2]}.`).join('\n')}\n*ELEMENT, TYPE=C3D4, ELSET=BODY\n${TETS.map((t, i) => `${i + 1}, ${t.map((v) => v + 1).join(', ')}`).join('\n')}\n*NSET, NSET=FIXED, GENERATE\n1, 4, 1\n*ELSET, ELSET=ALL\n1, 2, 3, 4, 5, 6\n*SOLID SECTION, ELSET=BODY, MATERIAL=AL7075\n*END PART\n*ASSEMBLY, NAME=A\n*INSTANCE, NAME=Block-1, PART=Block\n10., 0., 0.\n*END INSTANCE\n*END ASSEMBLY\n*MATERIAL, NAME=AL7075\n*ELASTIC\n71000., 0.33\n`;
  const m = await imp('cube.inp', keep('cube.inp', src)); fmtId('abaqus', m, 'abaqus', 'native'); const a = await cubeCheck('C3D4 part', m, { nv: 8 });
  ok('instance translation applied', Math.abs(a.bbox.min[0] - 10) < 1e-9, JSON.stringify(a.bbox.min));
  ok('element type mapped', m.elements[0].type === 'tet4' && m.meta.elementTypes.C3D4 === 6); ok('elset → group', m.groups[0].name === 'Block.BODY' && m.groups[0].count === 6, JSON.stringify(m.groups));
  ok('nset/elset census', m.meta.nsets[0].name === 'FIXED' && m.meta.nsets[0].count === 4 && m.meta.elsets[0].count === 6); ok('material name', m.meta.materials[0] === 'AL7075' && m.groups.some((g) => g.kind === 'material' && g.name === 'AL7075'));
  ok('heading + section', m.meta.heading === 'cube model' && m.meta.sections[0].material === 'AL7075'); ok('kind structural', m.kind === 'structural-mesh');
  const m2 = await imp('h.inp', `*NODE\n${CV.map((v, i) => `${i + 1}, ${v.join(', ')}`).join('\n')}\n*ELEMENT, TYPE=C3D8R\n1, 1, 2, 3, 4,\n5, 6, 7, 8\n*ELEMENT, TYPE=S4R, ELSET=SKIN\n2, 1, 2, 3, 4\n*ELEMENT, TYPE=SPRINGA\n3, 1, 2\n`);
  ok('flat deck: C3D8R over two lines + S4R; unknown type warned', m2.elements.find((e) => e.type === 'hex8')?.count === 1 && m2.elements.find((e) => e.type === 'quad4')?.count === 1 && m2.warnings.some((w) => /SPRINGA/.test(w)), JSON.stringify(m2.warnings));
}

sec('UNV');
{
  const I = (...v) => v.map((x) => String(x).padStart(10)).join('');
  const src = `    -1\n   164\n         5mm (milli-newton)  2\n  1.00000000000000000D+03  1.00000000000000000D+03  1.00000000000000000D+00\n  2.73149999999999977D+02\n    -1\n    -1\n  2411\n${CV.map((v, i) => `${I(i + 1, 1, 1, 11)}\n   ${v.map((c) => c.toExponential(16).replace('e', 'D')).join('   ')}`).join('\n')}\n    -1\n    -1\n  2412\n${TETS.map((t, i) => `${I(i + 1, 111, 2, 1, 7, 4)}\n${I(...t.map((v) => v + 1))}`).join('\n')}\n${I(7, 21, 2, 1, 7, 2)}\n${I(0, 0, 0)}\n${I(1, 7)}\n    -1\n    -1\n  2467\n${I(1, 0, 0, 0, 0, 0, 0, 3)}\nroot rib\n${I(8, 1, 0, 0, 8, 2, 0, 0)}\n${I(8, 3, 0, 0)}\n    -1\n`;
  const m = await imp('cube.unv', keep('cube.unv', src)); fmtId('unv', m, 'unv', 'native'); await cubeCheck('2411/2412', m, { nv: 8 });
  ok('dataset 164 → units mm', m.units.length === 'mm' && m.units.source === 'file' && /MM/.test(m.meta.unitSystem), JSON.stringify(m.units)); ok('D exponents parsed', m.positions[3] === 1);
  ok('tets + beam with orientation record', m.elements.find((e) => e.type === 'tet4')?.count === 6 && m.elements.find((e) => e.type === 'line2')?.count === 1, m.elements.map((e) => e.type + e.count).join());
  ok('permanent group', m.groups.find((g) => g.name === 'root rib')?.count === 3, JSON.stringify(m.groups)); ok('dataset census', m.meta.datasets[2411] === 1 && m.meta.datasets[2412] === 1);
}

sec('point clouds');
{
  let m = await imp('scan.xyz', keep('scan.xyz', '# x y z i\n' + CV.map((v) => v.join(' ') + ' 0.5').join('\n') + '\n')); fmtId('xyz', m, 'xyz', 'native');
  ok('xyz: 8 points, kind, units', m.positions.length === 24 && m.kind === 'pointcloud' && m.units.length === null && m.meta.columns === 4);
  const a = G.analyse(m); ok('xyz: analyse without faces', a.bbox.size.every((s) => s === 1) && a.nTris === 0 && a.volume === null && a.watertight === false && a.classification.label === 'point cloud / scan');
  m = await imp('pts.csv', 'x,y,z\n0,0,0\n1,2,3\n'); ok('csv with header', m.positions.length === 6 && m.positions[5] === 3 && m.format === 'xyz');
  m = await imp('naca.dat', 'NACA 0012\n' + naca(0.12, 30).map((p) => p.map((v) => v.toFixed(6)).join('  ')).join('\n') + '\n'); ok('two-column aerofoil table → planar polyline', m.format === 'xyz' && m.meta.planar && m.lines.length === 2 * (m.positions.length / 3 - 1) && m.warnings.some((w) => /planar/.test(w)), m.format);
  const pcdHead = (data) => `# .PCD v0.7 - Point Cloud Data file format\nVERSION 0.7\nFIELDS x y z intensity\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH 8\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS 8\nDATA ${data}\n`;
  m = await imp('c.pcd', keep('c.pcd', pcdHead('ascii') + CV.map((v) => v.join(' ') + ' 1').join('\n') + '\n')); fmtId('pcd', m, 'pcd', 'native'); ok('pcd ascii', m.positions.length === 24 && m.kind === 'pointcloud' && m.meta.width === 8);
  const pb = new Uint8Array(8 * 16), dv = new DataView(pb.buffer); CV.forEach((v, i) => v.forEach((c, k) => dv.setFloat32(16 * i + 4 * k, c, true))); dv.setFloat32(0, NaN, true);
  m = await imp('b.pcd', keep('b.pcd', cat(enc(pcdHead('binary')), pb))); ok('pcd binary, NaN point dropped + warned', m.positions.length === 21 && m.warnings.some((w) => /non-finite/.test(w)));
  m = await imp('z.pcd', pcdHead('binary_compressed') + 'xxxx'); ok('pcd compressed → metadata-only', m.kind === 'metadata-only' && m.meta.encoding === 'binary_compressed');
  const las = (n, fmtByte = 0) => { const b = new Uint8Array(227 + 20 * n), d = new DataView(b.buffer); b.set(enc('LASF')); b[24] = 1; b[25] = 2; b.set(enc('TEST SYSTEM'), 26); d.setUint16(94, 227, true); d.setUint32(96, 227, true); b[104] = fmtByte; d.setUint16(105, 20, true); d.setUint32(107, n, true); [0.01, 0.01, 0.01, 100, 200, 300].forEach((v, i) => d.setFloat64(131 + 8 * i, v, true)); for (let i = 0; i < n; i++) { d.setInt32(227 + 20 * i, i * 100, true); d.setInt32(231 + 20 * i, i * 50, true); d.setInt32(235 + 20 * i, 7, true); } return b; };
  m = await imp('scan.las', keep('scan.las', las(1000))); fmtId('las', m, 'las', 'partial'); ok('las: all points, scale/offset applied', m.positions.length === 3000 && Math.abs(m.positions[3] - 101) < 1e-9 && Math.abs(m.positions[4] - 200.5) < 1e-9 && Math.abs(m.positions[5] - 300.07) < 1e-9, Array.from(m.positions.slice(3, 6)).join());
  ok('las: header metadata', m.meta.version === '1.2' && m.meta.pointCount === 1000 && m.meta.systemIdentifier === 'TEST SYSTEM' && m.units.length === null);
  m = await imp('scan.las', las(1000), { maxPoints: 100 }); ok('las: stride subsampling to maxPoints', m.positions.length === 300 && m.meta.stride === 10 && m.warnings.some((w) => /subsampled/.test(w)));
  m = await imp('scan.laz', las(10, 0x80)); ok('laz → metadata-only with header', m.format === 'laz' && m.kind === 'metadata-only' && m.meta.pointCount === 10 && m.warnings.some((w) => /LASzip|pathway/i.test(w)));
}

// ---------- partial readers ----------
sec('VTK XML');
{
  const vtu = (fmt) => `<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n<UnstructuredGrid><Piece NumberOfPoints="8" NumberOfCells="6">\n<Points><DataArray type="Float64" NumberOfComponents="3" format="${fmt}">${CV.flat().join(' ')}</DataArray></Points>\n<Cells><DataArray type="Int32" Name="connectivity" format="${fmt}">${TETS.flat().join(' ')}</DataArray><DataArray type="Int32" Name="offsets" format="${fmt}">4 8 12 16 20 24</DataArray><DataArray type="UInt8" Name="types" format="${fmt}">10 10 10 10 10 10</DataArray></Cells>\n<CellData><DataArray type="Float64" Name="stress" format="${fmt}">1 2 3 4 5 6</DataArray></CellData>\n</Piece></UnstructuredGrid></VTKFile>\n`;
  let m = await imp('cube.vtu', keep('cube.vtu', vtu('ascii'))); fmtId('vtu', m, 'vtkxml', 'partial'); await cubeCheck('vtu ascii', m, { nv: 8 }); ok('vtu: tets + array names', m.elements[0].type === 'tet4' && m.meta.dataArrays.includes('stress'));
  m = await imp('cube.vtu', vtu('binary')); ok('vtu binary → metadata-only + warning + pathway', m.kind === 'metadata-only' && m.warnings.some((w) => /inline ASCII/.test(w)) && m.meta.pieces[0].NumberOfPoints === '8');
  m = await imp('cube.vtp', `<VTKFile type="PolyData"><PolyData><Piece NumberOfPoints="8" NumberOfPolys="6"><Points><DataArray format="ascii">${CV.flat().join(' ')}</DataArray></Points><Polys><DataArray Name="connectivity" format="ascii">${CQ.flat().join(' ')}</DataArray><DataArray Name="offsets" format="ascii">4 8 12 16 20 24</DataArray></Polys></Piece></PolyData></VTKFile>`);
  await cubeCheck('vtp polys', m, { nv: 8 });
  m = await imp('grid.vts', '<VTKFile type="StructuredGrid"><StructuredGrid WholeExtent="0 1 0 1 0 1"/></VTKFile>'); ok('vts recognised, not read', m.format === 'vtkxml' && m.kind === 'metadata-only' && m.meta.type === 'StructuredGrid');
}

sec('ANSYS CDB');
{
  const i9 = (...v) => v.map((x) => String(x).padStart(9)).join(''), e21 = (x) => x.toExponential(13).replace(/e([+-])(\d+)$/, (s, sg, d) => `E${sg}${d.padStart(3, '0')}`).padStart(21);
  const src = `/COM,ANSYS RELEASE 2023\n/PREP7\n/UNITS,MPA\nET,1,185\nET,2,181\nNBLOCK,6,SOLID,8,8\n(3i9,6e21.13e3)\n${CV.map((v, i) => i9(i + 1, 0, 0) + v.filter((c, k) => k < 2 || c !== 0 || true).map(e21).join('')).join('\n')}\nN,R5.3,LOC,       -1,\nEBLOCK,19,SOLID,3,3\n(19i9)\n${i9(1, 1, 1, 1, 0, 0, 0, 0, 8, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8)}\n${i9(1, 1, 1, 1, 0, 0, 0, 0, 8, 0, 2, 1, 2, 3, 3, 5, 5, 5, 5)}\n${i9(2, 2, 1, 1, 0, 0, 0, 0, 4, 0, 3, 1, 2, 3, 3)}\n       -1\nCMBLOCK,SKIN,ELEM,2\n(8i10)\n         1         3\nFINISH\n`;
  const m = await imp('cube.cdb', keep('cube.cdb', src)); fmtId('cdb', m, 'cdb', 'partial');
  const a = G.analyse(m); ok('nblock: 8 nodes, unit bbox', a.nVerts === 8 && a.bbox.size.every((s) => Math.abs(s - 1) < 1e-12), JSON.stringify(a.bbox.size));
  ok('eblock: hex + degenerate tet + degenerate shell tri', m.elements.map((e) => `${e.type}:${e.count}`).sort().join() === 'hex8:1,tet4:1,tri3:1', m.elements.map((e) => `${e.type}:${e.count}`).join());
  ok('ET table + units + components', m.meta.elementTypes[1] === 185 && m.units.length === 'mm' && m.meta.components[0].name === 'SKIN'); ok('material groups', m.groups.length === 2 && m.groups[0].kind === 'material' && m.groups[0].count === 2);
  ok('kind structural', m.kind === 'structural-mesh');
}

sec('Fluent');
{
  const hx = (v) => v.toString(16);
  const src = `(0 "written by test (with parens)")\n(2 3)\n(10 (0 1 8 0 3))\n(10 (1 1 8 1 3)(\n${CV.map((v) => v.map((c) => c.toExponential(8)).join(' ')).join('\n')}\n))\n(12 (0 1 1 0))\n(12 (2 1 1 1 4))\n(13 (0 1 ${hx(6)} 0))\n(13 (3 1 4 3 4)(\n${CQ.slice(0, 4).map((q) => q.map((v) => hx(v + 1)).join(' ') + ' 1 0').join('\n')}\n))\n(13 (4 5 6 a 0)(\n${CQ.slice(4).map((q) => '4 ' + q.map((v) => hx(v + 1)).join(' ') + ' 1 0').join('\n')}\n))\n(45 (2 fluid air)())\n(45 (3 wall wing-surface)())\n(45 (4 velocity-inlet inlet)())\n`;
  const m = await imp('case.msh', keep('fluent.msh', src)); fmtId('fluent', m, 'fluent', 'partial');
  await cubeCheck('boundary face zones', m, { nv: 8 }); ok('zone names → groups', m.groups.find((g) => g.name === 'wing-surface')?.count === 4 && m.groups.find((g) => g.name === 'inlet')?.count === 2 && m.groups.find((g) => g.kind === 'zone')?.name === 'air', JSON.stringify(m.groups));
  ok('census', m.meta.nodes === 8 && m.meta.faces === 6 && m.meta.cells === 1 && m.meta.dimension === 3); ok('honest: no volume elements', m.kind === 'surface-mesh' && !m.elements.some((e) => e.type === 'hex8') && m.warnings.some((w) => /NOT reconstructed/.test(w)));
}

sec('Tecplot');
{
  const src = `TITLE = "cube"\nVARIABLES = "X", "Y", "Z", "P"\nZONE T="block", N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK\n${CV.map((v) => v.join(' ') + ' 101325').join('\n')}\n1 2 3 4 5 6 7 8\nZONE T="tets" N=8 E=6 F=FEPOINT ET=TETRAHEDRON\n${CV.map((v) => `${v[0] + 2} ${v[1]} ${v[2]} 0`).join('\n')}\n${TETS.map((t) => t.map((v) => v + 1).join(' ')).join('\n')}\n`;
  let m = await imp('cube.dat', keep('tec.dat', src)); fmtId('tecplot', m, 'tecplot', 'partial');
  ok('two FE zones', m.elements.find((e) => e.type === 'hex8')?.count === 1 && m.elements.find((e) => e.type === 'tet4')?.count === 6 && m.groups.length === 2 && m.groups[1].name === 'tets', JSON.stringify(m.groups));
  const a = G.analyse(m); near('area of two cubes', a.area, 12, 1e-9); near('volume of two cubes', a.volume, 2, 1e-9); ok('variables listed', m.meta.variables.join() === 'X,Y,Z,P' && m.meta.title === 'cube');
  m = await imp('tri.dat', 'VARIABLES = x y\nZONE N=4, E=2, F=FEPOINT, ET=TRIANGLE\n0 0\n1 0\n1 1\n0 1\n1 2 3\n1 3 4\n'); ok('2-D FETRIANGLE zone', m.elements[0].type === 'tri3' && m.elements[0].count === 2 && Math.abs(G.analyse(m).area - 1) < 1e-12);
  m = await imp('ord.dat', 'VARIABLES = "X" "Y" "Z"\nZONE I=2, J=2, K=2, F=POINT\n0 0 0\n1 0 0\n0 1 0\n1 1 0\n0 0 1\n1 0 1\n0 1 1\n1 1 1\n'); await cubeCheck('ordered IJK zone', m, { nv: 8 });
  m = await imp('blk.dat', 'VARIABLES = "X" "Y" "Z"\nZONE N=8, E=1, DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n' + [0, 1, 2].map((k) => CV.map((v) => v[k]).join(' ')).join('\n') + '\n1 2 3 4 5 6 7 8\n'); await cubeCheck('BLOCK packing', m, { nv: 8 });
}

sec('Plot3D');
{
  const blk = (ox) => [0, 1, 2].map((k) => [0, 1, 2, 3, 4, 5, 6, 7].map((i) => [(i & 1) + ox, (i >> 1) & 1, (i >> 2) & 1][k].toFixed(6)).join(' ')).join('\n');
  let m = await imp('grid.xyz', keep('grid.xyz', `2 2 2\n${blk(0)}\n`)); fmtId('single block', m, 'plot3d', 'partial'); await cubeCheck('single block', m, { nv: 8 }); ok('hex cell', m.elements[0].type === 'hex8' && m.kind === 'volume-mesh');
  m = await imp('grid.p3d', keep('grid.p3d', `2\n2 2 2\n2 2 2\n${blk(0)}\n${blk(3)}\n`)); ok('multi-block', m.elements[0].count === 2 && m.groups.length === 2 && m.meta.blocks.length === 2 && Math.abs(G.analyse(m).volume - 2) < 1e-9);
  m = await imp('surf.x', '3 2 1\n0 1 2 0 1 2\n0 0 0 1 1 1\n0 0 0 0 0 0\n'); ok('k = 1 block → quads', m.elements[0].type === 'quad4' && m.elements[0].count === 2 && Math.abs(G.analyse(m).area - 2) < 1e-12);
}

sec('OpenFOAM');
{
  const hd = (cls, obj, note = '') => `/*--------------------------------*- C++ -*----------------------------------*\\\n| =========                 |\n\\*---------------------------------------------------------------------------*/\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       ${cls};\n${note}    location    "constant/polyMesh";\n    object      ${obj};\n}\n// * * * * * * * * * * * * * //\n\n`;
  const points = hd('vectorField', 'points') + `8\n(\n${CV.map((v) => `(${v.join(' ')})`).join('\n')}\n)\n`, faces = hd('faceList', 'faces') + `6\n(\n${CQ.map((q) => `4(${q.join(' ')})`).join('\n')}\n)\n`;
  const boundary = hd('polyBoundaryMesh', 'boundary') + '2\n(\n    walls\n    {\n        type            wall;\n        nFaces          4;\n        startFace       0;\n    }\n    inlet\n    {\n        type            patch;\n        nFaces          2;\n        startFace       4;\n    }\n)\n', owner = hd('labelList', 'owner', '    note        "nPoints:8  nCells:1  nFaces:6  nInternalFaces:0";\n') + '6\n(\n0\n0\n0\n0\n0\n0\n)\n';
  keep('points', points);
  let m = await imp('points', points, { companions: [{ name: 'faces', bytes: enc(faces) }, { name: 'boundary', bytes: enc(boundary) }, { name: 'owner', bytes: enc(owner) }] });
  fmtId('polyMesh', m, 'openfoam', 'partial'); await cubeCheck('points + faces + boundary', m, { nv: 8 });
  ok('patches → boundary groups', m.groups.length === 2 && m.groups[0].name === 'walls (wall)' && m.groups[0].count === 4 && m.groups[1].count === 2, JSON.stringify(m.groups)); ok('cell count from owner', m.meta.cells === 1);
  m = await imp('faces', faces, { companions: [{ name: 'points', bytes: enc(points) }] }); ok('main file may be faces; no boundary → warning', m.triangles.length === 36 && m.warnings.some((w) => /No "boundary" file/.test(w)));
  m = await imp('points', points); ok('points alone → metadata-only, asks for companions', m.kind === 'metadata-only' && m.warnings.some((w) => /opts\.companions/.test(w)));
}

sec('COLLADA / VRML / X3D');
{
  const dae = `<?xml version="1.0"?><COLLADA xmlns="http://www.collada.org/2005/11/COLLADASchema" version="1.4.1"><asset><contributor><authoring_tool>TestTool</authoring_tool></contributor><unit name="centimeter" meter="0.01"/><up_axis>Z_UP</up_axis></asset><library_geometries><geometry id="g1" name="Box"><mesh><source id="pos"><float_array id="pa" count="24">${CV.flat().join(' ')}</float_array><technique_common><accessor source="#pa" count="8" stride="3"/></technique_common></source><source id="nrm"><float_array id="na" count="3">0 0 1</float_array></source><vertices id="vtx"><input semantic="POSITION" source="#pos"/></vertices><polylist count="6"><input semantic="VERTEX" source="#vtx" offset="0"/><input semantic="NORMAL" source="#nrm" offset="1"/><vcount>4 4 4 4 4 4</vcount><p>${CQ.flat().map((v) => `${v} 0`).join(' ')}</p></polylist></mesh></geometry></library_geometries></COLLADA>`;
  let m = await imp('cube.dae', keep('cube.dae', dae)); fmtId('dae', m, 'dae', 'partial'); await cubeCheck('polylist with two inputs', m, { nv: 8 }); ok('dae: unit + up axis + name', m.units.length === 'cm' && m.meta.upAxis === 'Z_UP' && m.groups[0].name === 'Box' && m.meta.authoringTool === 'TestTool');
  const wrl = `#VRML V2.0 utf8\nShape { geometry IndexedFaceSet { coord Coordinate { point [ ${CV.map((v) => v.join(' ')).join(', ')} ] } coordIndex [ ${CQ.map((q) => q.join(', ') + ', -1').join(', ')} ] } }\n`;
  m = await imp('cube.wrl', keep('cube.wrl', wrl)); fmtId('vrml', m, 'vrml', 'partial'); await cubeCheck('IndexedFaceSet', m, { nv: 8 });
  const x3d = `<?xml version="1.0"?><X3D profile="Interchange" version="3.3"><Scene><Shape><IndexedFaceSet DEF="Hull" coordIndex="${CQ.map((q) => q.join(' ') + ' -1').join(' ')}"><Coordinate point="${CV.map((v) => v.join(' ')).join(', ')}"/></IndexedFaceSet></Shape></Scene></X3D>`;
  m = await imp('cube.x3d', keep('cube.x3d', x3d)); fmtId('x3d', m, 'x3d', 'partial'); await cubeCheck('X3D IndexedFaceSet', m, { nv: 8 }); ok('x3d DEF name', m.groups[0].name === 'Hull');
}

sec('STEP');
{
  const head = (schema) => `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('test part'),'2;1');\nFILE_NAME('cube.stp','2024-05-01T10:00:00',('A. Engineer'),('Test Org'),'TestPre 1.0','TestCAD 2024','');\nFILE_SCHEMA(('${schema}'));\nENDSEC;\nDATA;\n`;
  const brep = head('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }') + CV.map((v, i) => `#${i + 1}=CARTESIAN_POINT('',(${v.map((c) => (c * 25).toFixed(1)).join(',')}));`).join('\n') + `\n#20=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\n#21=( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );\n#22=( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#25)) GLOBAL_UNIT_ASSIGNED_CONTEXT((#20,#21)) REPRESENTATION_CONTEXT('','') );\n#30=PRODUCT('P-001','Wing rib; left','',(#31));\n#31=PRODUCT_CONTEXT('',#32,'mechanical');\n#40=ADVANCED_FACE('',(#41),#42,.T.);\n#43=ADVANCED_FACE('',(#41),#42,.F.);\n#44=MANIFOLD_SOLID_BREP('',#45);\n#46=B_SPLINE_SURFACE_WITH_KNOTS('',3,3,((#1,#2),(#3,#4)),.UNSPECIFIED.,.F.,.F.,.F.,(4,4),(4,4),(0.,1.),(0.,1.),.UNSPECIFIED.);\n#50=POLYLINE('',(#1,#2,#3,#4));\n#51=NEXT_ASSEMBLY_USAGE_OCCURRENCE('1','','',#60,#61,$);\nENDSEC;\nEND-ISO-10303-21;\n`;
  let m = await imp('cube.stp', keep('cube.stp', brep)); fmtId('B-rep', m, 'step', 'partial');
  ok('header: schema / AP / system', m.meta.ap === 'AP214' && /AUTOMOTIVE_DESIGN/.test(m.meta.schema[0]) && m.meta.originatingSystem === 'TestCAD 2024' && m.meta.author === 'A. Engineer' && m.meta.timeStamp === '2024-05-01T10:00:00', JSON.stringify(m.meta).slice(0, 300));
  ok('units from SI_UNIT via the global unit context', m.units.length === 'mm' && m.units.source === 'file' && m.meta.lengthUnit === 'MILLIMETRE', JSON.stringify(m.units));
  ok('product name (with semicolon inside the string)', m.meta.products[0] === 'Wing rib; left' && m.meta.assemblyUsages === 1, JSON.stringify(m.meta.products));
  ok('entity census', m.meta.census.CARTESIAN_POINT === 8 && m.meta.census.ADVANCED_FACE === 2 && m.meta.census.MANIFOLD_SOLID_BREP === 1 && m.meta.census.B_SPLINE_SURFACE_WITH_KNOTS === 1 && m.meta.census.SI_UNIT === 2, JSON.stringify(m.meta.census));
  const a = G.analyse(m); ok('cartesian points as cloud for extents', m.kind === 'cad-brep' && a.nVerts === 8 && a.bbox.size.every((s) => Math.abs(s - 25) < 1e-9) && m.triangles.length === 0, JSON.stringify(a.bbox.size));
  ok('polyline read', m.lines.length === 6); ok('says clearly that B-rep is not tessellated', m.warnings.some((w) => /NOT tessellated/.test(w)) && a.classification.label.startsWith('exact CAD'));
  const inch = brep.replace("#20=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );", "#18=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\n#19=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#18);\n#20=( CONVERSION_BASED_UNIT('INCH',#19) LENGTH_UNIT() NAMED_UNIT(#17) );");
  m = await imp('inch.step', inch); ok('conversion-based inch unit', m.units.length === 'in', JSON.stringify(m.units) + m.meta.lengthUnit);
  m = await imp('nounit.stp', brep.replace(/#20=.*\n/, '').replace('(#20,#21)', '(#21)')); ok('no length unit → null', m.units.length === null);
  const tess = head('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }') + `#1=COORDINATES_LIST('',8,(${CV.map((v) => `(${v.map((c) => c.toFixed(1)).join(',')})`).join(',')}));\n#2=TRIANGULATED_FACE('',#1,8,((0.,0.,1.)),$,(),(${CT.slice(0, 6).map((t) => `(${t.map((v) => v + 1).join(',')})`).join(',')}));\n#3=COMPLEX_TRIANGULATED_FACE('',#1,8,((0.,0.,1.)),$,(2,3,7,6,4,8,1,5),((1,2,4,3),(2,5,3,6)),((7,5,6,8,1)));\n#9=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) );\nENDSEC;\nEND-ISO-10303-21;\n`;
  m = await imp('tess.stp', keep('tess.stp', tess)); ok('AP242 detected', m.meta.ap === 'AP242' && m.units.length === 'm');
  ok('AP242 tessellated faces → triangles', m.kind === 'surface' && m.meta.tessellated.faces === 2 && m.triangles.length / 3 === 6 + 4 + 3, `${m.triangles.length / 3} triangles`);
  const at = G.analyse(m); near('tessellated area: 3 cube faces from triangles + 3.5 from strips/fan', at.area, 3 + 3.5, 1e-9);
}

sec('IGES');
{
  const L = (data, s, n) => data.padEnd(72).slice(0, 72) + s + String(n).padStart(7);
  const f8 = (...v) => v.map((x) => String(x).padStart(8)).join('');
  const glob = '1H,,1H;,4Htest,8Htest.igs,7HTestCAD,4Hprep,32,38,6,308,15,4Htest,1.0,1,4HINCH,1,0.01,15H20240101.000000,1E-6,10.0,6Hauthor,3Horg,11,0,15H20240101.000000;';
  const pts = [[0, 0, 0], [2, 0, 0], [2, 3, 0], [0, 3, 4]], out = [L('IGES test file', 'S', 1)];
  for (let i = 0; i * 72 < glob.length; i++) out.push(L(glob.slice(72 * i, 72 * i + 72), 'G', i + 1));
  const D = [], Pl = [];
  const ent = (type, params, status = '00000000') => { const de = 2 * (D.length / 2) + 1; Pl.push(params.padEnd(64).slice(0, 64) + String(de).padStart(8) + 'P' + String(Pl.length + 1).padStart(7)); D.push(f8(type, Pl.length, 0, 0, 0, 0, 0, 0) + status + 'D' + String(D.length + 1).padStart(7)); D.push(f8(type, 0, 0, 1, 0, '', '', '', 0) + 'D' + String(D.length + 1).padStart(7)); };
  for (const p of pts) ent(116, `116,${p.map((c) => c.toFixed(1)).join(',')},0;`);
  ent(110, '110,0.,0.,0.,2.,3.,4.;'); ent(144, '144,1,1,0,3;'); ent(126, '126,1,1,0,0,1,0,0.,0.,1.,1.,1.,1.,50.,60.,0.,70.,80.,0.,0.,1.;', '00010500');
  out.push(...D, ...Pl, L(`S${'1'.padStart(7)}G${String(out.length - 1).padStart(7)}D${String(D.length).padStart(7)}P${String(Pl.length).padStart(7)}`, 'T', 1));
  const src = out.join('\n') + '\n';
  ok('fixture is 80 columns', out.every((l) => l.length === 80));
  const m = await imp('part.igs', keep('part.igs', src)); fmtId('iges', m, 'iges', 'partial');
  ok('global section: units inch + system', m.units.length === 'in' && m.units.source === 'file' && m.meta.originatingSystem === 'TestCAD' && m.meta.author === 'author' && m.meta.fileName === 'test.igs', JSON.stringify({ u: m.units, s: m.meta.originatingSystem, a: m.meta.author, f: m.meta.fileName }));
  ok('directory census', m.meta.census['116 point'] === 4 && m.meta.census['110 line'] === 1 && m.meta.census['144 trimmed parametric surface'] === 1 && m.meta.entityCount === 7, JSON.stringify(m.meta.census));
  const a = G.analyse(m); ok('116 points + 110 line read', m.meta.read.points === 4 && m.meta.read.lines === 1 && a.nVerts === 6 && m.lines.length === 2 && a.bbox.size.join() === '2,3,4', `${a.nVerts} verts, bbox ${a.bbox.size.join()}`);
  ok('parameter-space curve excluded from extents', m.meta.skippedParametricCurves === 1); ok('kind cad-brep + clear warning', m.kind === 'cad-brep' && m.warnings.some((w) => /NOT tessellated/.test(w)));
}

sec('Parasolid text');
{
  const src = '**ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz**************************\n**PARASOLID !"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~0123456789**************************\n**PART1;\nMC=x64;\nMC_MODEL=unknown;\nAPPL=TestCAD 2024;\nFORMAT=text;\nGUISE=transmit;\nDATE=1-may-2024;\n**PART2;\nSCH=SCH_3400201_34000;\nUSFLD_SIZE=0;\n**PART3;\n**END_OF_HEADER*****************************************************************\nT51 : TRANSMIT FILE created by modeller version 3400201 17 SCH_3400201_34000_13006 0\n12 1 12 0 2 0 0 0 0 1e3 1e-8 0 0 0 1 0 3 1 3 4 5 0 6 7 0\n';
  const m = await imp('part.x_t', keep('part.x_t', src)); fmtId('x_t', m, 'xt', 'metadata');
  ok('header fields + schema + version', m.meta.header.APPL === 'TestCAD 2024' && m.meta.header.FORMAT === 'text' && m.meta.schema.startsWith('SCH_3400201') && m.meta.modellerVersion === '3400201', JSON.stringify(m.meta));
  ok('honest: metadata only, no invented census', m.kind === 'metadata-only' && m.meta.entityCensus === null && m.positions.length === 0 && m.warnings.some((w) => /NOT decoded/.test(w)) && m.warnings.some((w) => /pathway/i.test(w)));
}

// ---------- kernel-backed readers (vendored WebAssembly) ----------
// STEP AP214 of a box with real B-rep topology (vertices, line edges, planar faces, closed shell, product, colour)
function stepBox(sx, sy, sz, { name = 'Box', unit = '.MILLI.', rgb = [1, 0, 0] } = {}) {
  const L = []; let id = 100; const add = (s) => { L.push(`#${++id}=${s};`); return `#${id}`; }, f = (v) => (Number.isInteger(v) ? v + '.' : String(v));
  const P = CV.map((v) => [v[0] * sx, v[1] * sy, v[2] * sz]), pt = (p) => add(`CARTESIAN_POINT('',(${p.map(f).join(',')}))`), dir = (d) => add(`DIRECTION('',(${d.map(f).join(',')}))`);
  const vx = P.map((p) => add(`VERTEX_POINT('',${pt(p)})`)), edges = new Map();
  const edge = (a, b) => { const k = Math.min(a, b) + '_' + Math.max(a, b); if (!edges.has(k)) { const i = Math.min(a, b), j = Math.max(a, b), d = P[j].map((c, q) => c - P[i][q]), l = Math.hypot(...d); edges.set(k, add(`EDGE_CURVE('',${vx[i]},${vx[j]},${add(`LINE('',${pt(P[i])},${add(`VECTOR('',${dir(d.map((c) => c / l))},${f(l)})`)})`)},.T.)`)); } return add(`ORIENTED_EDGE('',*,*,${edges.get(k)},${a < b ? '.T.' : '.F.'})`); };
  const faces = CQ.map((q) => {
    const u = P[q[1]].map((c, i) => c - P[q[0]][i]), w = P[q[3]].map((c, i) => c - P[q[0]][i]), n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]], nl = Math.hypot(...n), ul = Math.hypot(...u);
    const plane = add(`PLANE('',${add(`AXIS2_PLACEMENT_3D('',${pt(P[q[0]])},${dir(n.map((c) => c / nl))},${dir(u.map((c) => c / ul))})`)})`);
    return add(`ADVANCED_FACE('',(${add(`FACE_OUTER_BOUND('',${add(`EDGE_LOOP('',(${q.map((a, i) => edge(a, q[(i + 1) % 4])).join(',')}))`)},.T.)`)}),${plane},.T.)`);
  });
  const solid = add(`MANIFOLD_SOLID_BREP('${name}',${add(`CLOSED_SHELL('',(${faces.join(',')}))`)})`);
  const lu = add(`( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(${unit},.METRE.) )`), au = add('( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )'), su = add('( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() )');
  const unc = add(`UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),${lu},'distance_accuracy_value','')`);
  const ctx = add(`( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((${unc})) GLOBAL_UNIT_ASSIGNED_CONTEXT((${lu},${au},${su})) REPRESENTATION_CONTEXT('','') )`);
  const absr = add(`ADVANCED_BREP_SHAPE_REPRESENTATION('${name}',(${solid}),${ctx})`);
  const app = add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')"); add(`APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2010,${app})`);
  const prod = add(`PRODUCT('${name}','${name}','',(${add(`PRODUCT_CONTEXT('',${app},'mechanical')`)}))`); add(`PRODUCT_RELATED_PRODUCT_CATEGORY('part','',(${prod}))`);
  const pd = add(`PRODUCT_DEFINITION('design','',${add(`PRODUCT_DEFINITION_FORMATION('','',${prod})`)},${add(`PRODUCT_DEFINITION_CONTEXT('part definition',${app},'design')`)})`);
  add(`SHAPE_DEFINITION_REPRESENTATION(${add(`PRODUCT_DEFINITION_SHAPE('','',${pd})`)},${absr})`);
  if (rgb) add(`STYLED_ITEM('',(${add(`PRESENTATION_STYLE_ASSIGNMENT((${add(`SURFACE_STYLE_USAGE(.BOTH.,${add(`SURFACE_SIDE_STYLE('',(${add(`SURFACE_STYLE_FILL_AREA(${add(`FILL_AREA_STYLE('',(${add(`FILL_AREA_STYLE_COLOUR('',${add(`COLOUR_RGB('',${rgb.map(f).join(',')})`)})`)}))`)})`)}))`)})`)}))`)}),${solid})`);
  return `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('box'),'2;1');\nFILE_NAME('box.stp','2024-05-01T10:00:00',('A. Engineer'),('Test Org'),'TestPre 1.0','TestCAD 2024','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n${L.join('\n')}\nENDSEC;\nEND-ISO-10303-21;\n`;
}
// IGES file of rational B-spline surfaces (entity 128), one parameter record per surface
function igesSurfaces(records, units = '2,2HMM') {
  const L8 = (data, s, n) => data.padEnd(72).slice(0, 72) + s + String(n).padStart(7), f8 = (...v) => v.map((x) => String(x).padStart(8)).join('');
  const glob = `1H,,1H;,4Htest,8Htest.igs,7HTestCAD,4Hprep,32,38,6,308,15,4Htest,1.0,${units},1,0.01,15H20240101.000000,1E-6,10.0,6Hauthor,3Horg,11,0,15H20240101.000000;`;
  const out = [L8('IGES test file', 'S', 1)], D = [], Pl = [];
  for (let i = 0; i * 72 < glob.length; i++) out.push(L8(glob.slice(72 * i, 72 * i + 72), 'G', i + 1));
  for (const par of records) {
    const de = D.length + 1, first = Pl.length + 1, toks = par.split(','); let line = '', n = 0;
    const flush = () => { Pl.push(line.padEnd(64) + String(de).padStart(8) + 'P' + String(Pl.length + 1).padStart(7)); line = ''; n++; };
    toks.forEach((t, i) => { const piece = t + (i < toks.length - 1 ? ',' : ''); if (line.length + piece.length > 64) flush(); line += piece; }); flush();
    D.push(f8(128, first, 0, 0, 0, 0, 0, 0) + '00000000D' + String(D.length + 1).padStart(7)); D.push(f8(128, 0, 0, n, 0, '', '', '', 0) + 'D' + String(D.length + 1).padStart(7));
  }
  out.push(...D, ...Pl, L8(`S${'1'.padStart(7)}G${String(out.length - 1).padStart(7)}D${String(D.length).padStart(7)}P${String(Pl.length).padStart(7)}`, 'T', 1));
  return out.join('\n') + '\n';
}
const bilinear = (q, P) => `128,1,1,1,1,0,0,1,0,0,0.,0.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,${[q[0], q[1], q[3], q[2]].map((v) => P[v].map((x) => x.toFixed(1)).join(',')).join(',')},0.,1.,0.,1.;`;
// one planar rectangular face a × b in OpenCASCADE's text BREP format (shapes are numbered from the end)
function brepRect(a, b) {
  const V = [[0, 0, 0], [a, 0, 0], [a, b, 0], [0, b, 0]], E = [[0, 1, a, [1, 0, 0]], [1, 2, b, [0, 1, 0]], [2, 3, a, [-1, 0, 0]], [3, 0, b, [0, -1, 0]]], N = 10, ref = (k) => N - k;
  const ve = V.map((p) => `Ve\n1e-07\n${p.join(' ')}\n0 0\n\n0101101\n*`), ed = E.map(([s, e, l], i) => `Ed\n 1e-07 1 1 0\n1  ${i + 1} 0 0 ${l}\n0\n\n0101000\n+${ref(s)} 0 -${ref(e)} 0 *`);
  return `DBRep_DrawableShape\n\nCASCADE Topology V1, (c) Matra-Datavision\nLocations 0\nCurve2ds 0\nCurves 4\n${E.map(([s, , , d]) => `1 ${V[s].join(' ')} ${d.join(' ')} `).join('\n')}\nPolygon3D 0\nPolygonOnTriangulations 0\nSurfaces 1\n1 0 0 0 0 0 1 1 0 0 0 1 0 \nTriangulations 0\n\nTShapes ${N}\n${ve.join('\n')}\n${ed.join('\n')}\nWi\n\n0101100\n+${ref(4)} 0 +${ref(5)} 0 +${ref(6)} 0 +${ref(7)} 0 *\nFa\n0  1e-07 1 0\n\n0111000\n+${ref(8)} 0 *\n\n+1 0 \n`;
}

sec('STEP / IGES / BREP exact geometry (OpenCASCADE)');
{
  const src = stepBox(20, 30, 40, { name: 'Wing rib', rgb: [1, 0, 0] });
  let m = await imp('box.stp', keep('box.stp', src)), a = G.analyse(m); fmtId('box', m, 'step', 'native');
  ok('B-rep faces tessellated: 6 planar faces → 12 triangles', m.kind === 'surface' && a.nTris === 12 && m.meta.tessellation.faces === 6 && m.meta.tessellation.bodies === 1, `${m.kind}, ${a.nTris} triangles, ${JSON.stringify(m.warnings)}`);
  ok('bbox 20 × 30 × 40', a.bbox.size.every((s, i) => Math.abs(s - [20, 30, 40][i]) < 1e-9), JSON.stringify(a.bbox.size)); near('area', a.area, 2 * (20 * 30 + 30 * 40 + 20 * 40), 1e-9); near('volume', a.volume, 24000, 1e-9);
  ok('closed, outward', a.watertight && !a.inwardNormals && a.components === 1); ok('units from the text parser kept', m.units.length === 'mm' && m.units.source === 'file');
  ok('body group named after the product, with colour', m.groups.length === 1 && m.groups[0].name === 'Wing rib' && m.groups[0].kind === 'solid' && m.groups[0].count === 12 && m.groups[0].faces === 6 && m.groups[0].color?.[0] === 1 && m.groups[0].color[2] === 0, JSON.stringify(m.groups));
  ok('per-face triangle ranges', m.meta.brepFaces.length === 6 && m.meta.brepFaces[1].first === 2 && m.meta.brepFaces[1].last === 3 && m.elements[0].face.length === 12 && m.elements[0].face[11] === 5);
  ok('header + census merged from the text parser', m.meta.ap === 'AP214' && m.meta.originatingSystem === 'TestCAD 2024' && m.meta.census.ADVANCED_FACE === 6 && m.meta.census.MANIFOLD_SOLID_BREP === 1 && m.meta.products[0] === 'Wing rib' && /OpenCASCADE/.test(m.meta.kernel));
  ok('no "not tessellated" warning when the kernel succeeded', !m.warnings.some((w) => /NOT tessellated/.test(w)), m.warnings.join(' | '));
  const mp = G.massProperties(m, 1); ok('mass properties of the STEP solid', Math.abs(mp.volume - 24000) < 1e-6 && mp.closed && Math.abs(mp.cg[0] - 10) < 1e-9 && Math.abs(mp.cg[2] - 20) < 1e-9);
  const sl = G.slice(m, { axis: 'z', value: 12 }); ok('slice of the STEP solid', sl.length === 1 && Math.abs(G.sectionMetrics(sl[0]).area - 600) < 1e-9);
  m = await imp('box_m.stp', stepBox(0.5, 0.5, 2, { unit: '$' })); a = G.analyse(m); ok('metre-unit STEP keeps file coordinates', m.units.length === 'm' && Math.abs(a.volume - 0.5) < 1e-12 && Math.abs(a.bbox.size[2] - 2) < 1e-12, `${m.units.length} ${a.volume}`);
  m = await imp('box.stp', src, { wasm: false }); ok('opts.wasm = false → text-level fallback, support downgraded, clear warning', m.kind === 'cad-brep' && m.support === 'partial' && m.triangles.length === 0 && m.positions.length > 0 && m.warnings.some((w) => /switched off/.test(w)) && m.warnings.some((w) => /NOT tessellated/.test(w)), `${m.kind}/${m.support}`);
  m = await imp('broken.stp', src.replace(/#\d+=CLOSED_SHELL[^\n]*\n/, '')); ok('kernel finds nothing usable → graceful fallback with reason', m.kind === 'cad-brep' && m.support === 'partial' && m.warnings.some((w) => /could not be tessellated/.test(w)), `${m.kind} ${m.warnings.join(' | ')}`);

  const P10 = CV.map((v) => v.map((c) => c * 10)), ig = igesSurfaces(CQ.map((q) => bilinear(q, P10)));
  m = await imp('faces.igs', keep('faces.igs', ig)); a = G.analyse(m); fmtId('iges surfaces', m, 'iges', 'native');
  ok('six B-spline faces tessellated', m.kind === 'surface' && a.nTris === 12 && m.groups.length === 6 && m.groups[0].kind === 'surface', `${m.kind} ${a.nTris} ${m.groups.length} ${m.warnings.join('|')}`); near('IGES area', a.area, 600, 1e-9); ok('IGES bbox + units', a.bbox.size.every((s) => Math.abs(s - 10) < 1e-9) && m.units.length === 'mm');
  ok('IGES census kept from the text parser', m.meta.census['128 rational B-spline surface'] === 6 && m.meta.originatingSystem === 'TestCAD');
  const hi = G.heal(m); ok('unsewn IGES faces weld into a closed box', G.analyse(hi.model).watertight && Math.abs(G.analyse(hi.model).volume - 1000) < 1e-6);
  // a curved (parabolic) degree 2 × 1 patch: tessellation density follows the deflection options
  const curved = igesSurfaces(['128,2,1,2,1,0,0,1,0,0,0.,0.,0.,1.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,10.,0.,10.,20.,0.,0.,0.,10.,0.,10.,10.,10.,20.,10.,0.,0.,1.,0.,1.;']);
  const coarse = await imp('curved.igs', curved, { linearDeflection: 0.05, angularDeflection: 1 }), fine = await imp('curved.igs', curved, { linearDeflection: 0.0005, angularDeflection: 0.05 });
  let arc = 0; for (let i = 0; i < 20000; i++) { const t = (i + 0.5) / 20000; arc += Math.hypot(20, 20 - 40 * t) / 20000; }
  ok('finer deflection → more triangles', fine.triangles.length > 2 * coarse.triangles.length && coarse.triangles.length >= 6, `${coarse.triangles.length / 3} → ${fine.triangles.length / 3}`);
  ok('fine tessellation converges on the exact area', Math.abs(G.analyse(fine).area - 10 * arc) < 1e-3 * 10 * arc && Math.abs(G.analyse(fine).area - 10 * arc) < Math.abs(G.analyse(coarse).area - 10 * arc), `${G.analyse(coarse).area} / ${G.analyse(fine).area} vs ${10 * arc}`);
  ok('deflection recorded in the metadata', fine.meta.tessellation.linearDeflection === 0.0005 && fine.meta.tessellation.angularDeflection === 0.05);

  m = await imp('face.brep', keep('face.brep', brepRect(2, 3))); a = G.analyse(m); fmtId('brep', m, 'brep', 'native');
  ok('BREP face tessellated', m.kind === 'surface' && m.groups[0].kind === 'surface' && a.nTris === 2 && Math.abs(a.area - 6) < 1e-12 && a.bbox.size[0] === 2 && a.bbox.size[1] === 3, `${m.kind} ${a.nTris} ${a.area}`); ok('BREP has no units → null + note', m.units.length === null && m.warnings.some((w) => /does not state a length unit/.test(w)));
  ok('BREP detected by content', G.detectFormat('x.txt', enc(brepRect(1, 1))).id === 'brep');
  m = await imp('face.brep', brepRect(2, 3), { wasm: false }); ok('BREP without the kernel → metadata-only + pathway', m.kind === 'metadata-only' && m.support === 'partial' && m.meta.topologyVersion === 1);
}

sec('LAZ (LASzip)');
{
  // 500 points on a 10 × 10 × 5 lattice (x = 100 + 0.1 i, y = 200 + 0.1 j, z = 0.25 k), LAS 1.2 format 1, written with laspy/lazrs
  const b64 = 'TEFTRgAAAAAAAAAAAAAAAAAAAAAAAAAAAQJPVEhFUgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGxhc3B5IDIuNy4wAAAAAAAAAAAAAAAAAAAAAAAAAAAAGgHqB+MARwEAAAEAAACBHAD0AQAAAAAAAAAAAAAAAAAAAAAAAAAAAAD8qfHSTWJQP/yp8dJNYlA//Knx0k1iUD8AAAAAAABZQAAAAAAAAGlAAAAAAAAAAACamZmZmTlZQAAAAAAAAFlAzczMzMwcaUAAAAAAAABpQAAAAAAAAPA/AAAAAAAAAAAAAGxhc3ppcCBlbmNvZGVkAAC8Vi4AaHR0cDovL2xhc3ppcC5vcmcAAAAAAAAAAAAAAAAAAAACAAAAAgIAAAAAAABQwwAA/////////////////////wIABgAUAAIABwAIAAIA/QIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQCUOwGUkaZQ06TbGD+YXysigGw5c7LuWV2LA4VX4IOYMl6ceDAZ48TYhCn4OAYCkpkNbRqM+9G4F80ZAnjHFwHXAGn2R31q07SrvLw6tzpFIW+L03dLGhav7s+0D9kvPke629a5YkpUzOD/57PHW2/5Gpld23NhUdO4kiXDpDNmyC/2F2vPSdmZ7bgmNL6GQ8oL1QYG+ARZpcTVB4nFmO2dnLtjVvOghdHnmRgBCuSxCsUAOTLFjdQ0jv6okDdB40laainEI5JEGxjG3qkrn6tic1ufM3/DvBi0aAwVf+nRzYFCN/AkiKSu9rN782f0O25gwCiVztuF9XrP4chmPpozRhblig9oBowZnVI6vknp4zSoaPq46j4bjnsgGreSmXD7A71Ir7plFWgeSAIbg/iMW35svAYXCY8E7AlS7eJhOEk9HOdUTkwxfUW9iAnGD7UiByEcM7LxNMqjbaWM6NC9RVOJJpQmJ0FstnghXYY+ydagUkmF4UlDoqvh/zmHHkQosGzJGtBpanM3bxEA7AAAAAAAAAAEAAABMUAAAAA==';
  let laz = null; try { laz = Uint8Array.from(Buffer.from(b64, 'base64')); } catch { laz = null; }
  ok('fixture decodes', laz && laz.length > 300 && laz[0] === 0x4c, laz ? String(laz.length) : 'null');
  ok('LAZ detected by content', G.detectFormat('scan.las', laz).id === 'laz', G.detectFormat('scan.las', laz).note);
  let m = await imp('grid.laz', keep('grid.laz', laz)), a = G.analyse(m); fmtId('laz', m, 'laz', 'native');
  ok('all 500 points decompressed', m.kind === 'pointcloud' && a.nVerts === 500 && m.meta.pointCount === 500 && /LASzip/.test(m.meta.decoder), `${m.kind} ${a.nVerts} ${m.warnings.join('|')}`);
  ok('coordinates exact (scale/offset applied)', Math.abs(a.bbox.min[0] - 100) < 1e-9 && Math.abs(a.bbox.max[0] - 100.9) < 1e-9 && Math.abs(a.bbox.min[1] - 200) < 1e-9 && Math.abs(a.bbox.max[1] - 200.9) < 1e-9 && Math.abs(a.bbox.max[2] - 1) < 1e-9, JSON.stringify(a.bbox));
  ok('point 137 = (100.7, 200.3, 0.25)', Math.abs(m.positions[3 * 137] - 100.7) < 1e-9 && Math.abs(m.positions[3 * 137 + 1] - 200.3) < 1e-9 && Math.abs(m.positions[3 * 137 + 2] - 0.25) < 1e-9);
  m = await imp('grid.laz', laz, { maxPoints: 50 }); ok('subsampled to opts.maxPoints by stride', m.positions.length === 150 && m.meta.stride === 10 && m.warnings.some((w) => /subsampled/.test(w)));
  m = await imp('grid.laz', laz, { wasm: false }); ok('decoder switched off → header-only metadata', m.kind === 'metadata-only' && m.meta.pointCount === 500);
  const bad = laz.slice(); for (let i = 400; i < bad.length; i++) bad[i] ^= 0x5a; m = await imp('bad.laz', bad); ok('corrupt compressed data → no throw', Array.isArray(m.warnings) && m.positions.length % 3 === 0);
  m = await imp('grid.laz', laz); ok('decoder still works after a corrupt file', m.positions.length === 1500);
}

// ---------- HDF5 / NetCDF containers ----------
const h5 = await loadH5();
let h5n = 0;
/** Build an HDF5 file in the kernel's in-memory file system and return its bytes. */
const h5file = (fill) => { const p = `/fixture_${++h5n}.h5`, f = new h5.File(p, 'w'); try { fill(f); } finally { f.close(); } const b = h5.FS.readFile(p).slice(); h5.FS.unlink(p); return b; };
const chars = (s, n = s.length) => Int8Array.from({ length: n }, (_, i) => (i < s.length ? s.charCodeAt(i) : 32));
/** One CGNS node: an HDF5 group with name/label/type attributes and a " data" dataset. */
const cg = (parent, name, label, type, data) => { const g = parent.create_group(name); g.create_attribute('name', name); g.create_attribute('label', label); g.create_attribute('type', type); if (data !== undefined) g.create_dataset({ name: ' data', data }); return g; };
const cgCoords = (zone, pts) => { const gc = cg(zone, 'GridCoordinates', 'GridCoordinates_t', 'MT'); ['X', 'Y', 'Z'].forEach((ax, k) => cg(gc, 'Coordinate' + ax, 'DataArray_t', 'R8', Float64Array.from(pts, (p) => p[k]))); };
const cgBase = (f, units) => { cg(f, 'CGNSLibraryVersion', 'CGNSLibraryVersion_t', 'R4', new Float32Array([4.2])); const b = cg(f, 'Base', 'CGNSBase_t', 'I4', new Int32Array([3, 3])); if (units) cg(b, 'DimensionalUnits', 'DimensionalUnits_t', 'C1', chars(['Kilogram', units, 'Second', 'Kelvin', 'Radian'].map((w) => w.padEnd(32)).join(''))); return b; };

sec('CGNS (HDF5)');
{
  const unstructured = h5file((f) => {
    const base = cgBase(f, 'Meter'); cg(base, 'FAR', 'Family_t', 'MT');
    const z = cg(base, 'Fluid', 'Zone_t', 'I4', new Int32Array([8, 6, 0])); cg(z, 'ZoneType', 'ZoneType_t', 'C1', chars('Unstructured')); cgCoords(z, CV);
    const e1 = cg(z, 'Interior', 'Elements_t', 'I4', new Int32Array([10, 0])); cg(e1, 'ElementRange', 'IndexRange_t', 'I4', new Int32Array([1, 6])); cg(e1, 'ElementConnectivity', 'DataArray_t', 'I4', Int32Array.from(TETS.flat(), (v) => v + 1));
    const e2 = cg(z, 'Skin', 'Elements_t', 'I4', new Int32Array([5, 0])); cg(e2, 'ElementRange', 'IndexRange_t', 'I4', new Int32Array([7, 18])); cg(e2, 'ElementConnectivity', 'DataArray_t', 'I4', Int32Array.from(BT.flat(), (v) => v + 1));
    const zbc = cg(z, 'ZoneBC', 'ZoneBC_t', 'MT');
    const b1 = cg(zbc, 'wall', 'BC_t', 'C1', chars('BCWall')); cg(b1, 'PointRange', 'IndexRange_t', 'I4', new Int32Array([7, 16])); cg(b1, 'GridLocation', 'GridLocation_t', 'C1', chars('FaceCenter'));
    const b2 = cg(zbc, 'outer', 'BC_t', 'C1', chars('BCFarfield')); cg(b2, 'PointList', 'IndexArray_t', 'I4', new Int32Array([17, 18])); cg(b2, 'GridLocation', 'GridLocation_t', 'C1', chars('FaceCenter')); cg(b2, 'FamilyName', 'FamilyName_t', 'C1', chars('FAR'));
    const b3 = cg(zbc, 'probe', 'BC_t', 'C1', chars('BCGeneral')); cg(b3, 'PointList', 'IndexArray_t', 'I4', new Int32Array([1, 2]));
  });
  ok('fixture is a real HDF5 file', unstructured[0] === 0x89 && unstructured.length > 2000);
  let m = await imp('cube.cgns', keep('cube.cgns', unstructured)); fmtId('unstructured', m, 'cgns', 'partial'); await cubeCheck('TETRA_4 + TRI_3 sections', m, { nv: 8 });
  ok('elements', m.kind === 'volume-mesh' && m.elements.find((e) => e.type === 'tet4')?.count === 6 && m.elements.find((e) => e.type === 'tri3')?.count === 12 && m.triangles.length === 36, m.elements.map((e) => e.type + e.count).join());
  ok('BC by PointRange → group', m.groups.find((g) => g.name === 'wall')?.count === 10 && m.groups.find((g) => g.name === 'wall').bcType === 'BCWall', JSON.stringify(m.groups));
  ok('BC by PointList with a family → group named after the family', m.groups.find((g) => g.name === 'FAR')?.count === 2 && m.groups.find((g) => g.name === 'FAR').bcName === 'outer');
  ok('volume section → zone group', m.groups.find((g) => g.name === 'Interior')?.kind === 'zone' && m.groups.find((g) => g.name === 'Interior').count === 6);
  ok('vertex-located BC listed, not mapped', m.meta.boundaryConditions.length === 3 && m.warnings.some((w) => /defined on vertices/.test(w)));
  ok('DimensionalUnits → metres; version; families', m.units.length === 'm' && m.units.source === 'file' && Math.abs(m.meta.cgnsVersion - 4.2) < 1e-6 && m.meta.families[0] === 'FAR' && m.meta.zones[0].sections.length === 2, JSON.stringify(m.units) + m.meta.cgnsVersion);
  ok('classified as CFD volume mesh', G.analyse(m).classification.label === 'CFD volume mesh'); ok('quality: no inverted cells', G.meshQuality(m).perType.find((p) => p.type === 'tet4').negJacobian === 0);
  ok('schema taken from the content, not the extension', (await imp('mystery.h5', unstructured)).format === 'cgns' && (await imp('wrong.med', unstructured)).format === 'cgns');
  const structured = h5file((f) => { const base = cgBase(f, null), z = cg(base, 'Block', 'Zone_t', 'I4', new Int32Array([3, 2, 2, 2, 1, 1, 0, 0, 0])); cg(z, 'ZoneType', 'ZoneType_t', 'C1', chars('Structured')); const pts = []; for (let k = 0; k < 2; k++) for (let j = 0; j < 2; j++) for (let i = 0; i < 3; i++) pts.push([i * 0.5, j, k]); cgCoords(z, pts); });
  m = await imp('block.cgns', structured); await cubeCheck('structured zone', m, { nv: 12 }); ok('structured zone → 2 hexahedra', m.elements[0].type === 'hex8' && m.elements[0].count === 2 && m.units.length === null && m.meta.zones[0].dimensions.join() === '3,2,2');
  const mixed = h5file((f) => { const z = cg(cgBase(f, 'Millimeter'), 'Z', 'Zone_t', 'I4', new Int32Array([8, 1, 0])); cg(z, 'ZoneType', 'ZoneType_t', 'C1', chars('Unstructured')); cgCoords(z, CV); const e = cg(z, 'Mixed', 'Elements_t', 'I4', new Int32Array([20, 0])); cg(e, 'ElementRange', 'IndexRange_t', 'I4', new Int32Array([1, 3])); cg(e, 'ElementConnectivity', 'DataArray_t', 'I4', new Int32Array([17, 1, 2, 3, 4, 5, 6, 7, 8, 7, 1, 4, 3, 2, 3, 1, 2])); cg(e, 'ElementStartOffset', 'DataArray_t', 'I4', new Int32Array([0, 9, 14, 17])); });
  m = await imp('mixed.cgns', mixed); await cubeCheck('MIXED section', m, { nv: 8 }); ok('MIXED: hex + quad + bar, millimetres', m.elements.map((e) => e.type + ':' + e.count).sort().join() === 'hex8:1,line2:1,quad4:1' && m.units.length === 'mm', m.elements.map((e) => e.type + ':' + e.count).join());
  const poly = h5file((f) => { const z = cg(cgBase(f, 'Meter'), 'Poly', 'Zone_t', 'I8', new BigInt64Array([8n, 1n, 0n])); cg(z, 'ZoneType', 'ZoneType_t', 'C1', chars('Unstructured')); cgCoords(z, CV); const e = cg(z, 'Faces', 'Elements_t', 'I4', new Int32Array([22, 0])); cg(e, 'ElementRange', 'IndexRange_t', 'I4', new Int32Array([1, 6])); cg(e, 'ElementConnectivity', 'DataArray_t', 'I8', BigInt64Array.from(CQ.flat(), (v) => BigInt(v + 1))); cg(e, 'ElementStartOffset', 'DataArray_t', 'I8', new BigInt64Array([0n, 4n, 8n, 12n, 16n, 20n, 24n])); const c = cg(z, 'Cells', 'Elements_t', 'I4', new Int32Array([23, 0])); cg(c, 'ElementRange', 'IndexRange_t', 'I4', new Int32Array([7, 7])); cg(c, 'ElementConnectivity', 'DataArray_t', 'I4', new Int32Array([1, 2, 3, 4, 5, 6])); cg(c, 'ElementStartOffset', 'DataArray_t', 'I4', new Int32Array([0, 6])); });
  m = await imp('poly.cgns', poly); await cubeCheck('NGON_n / NFACE_n boundary (64-bit integers)', m, { nv: 8 }); ok('polyhedral zone: surface only, said clearly', m.kind === 'surface-mesh' && m.meta.polyhedral && m.meta.zones[0].polyhedralCells === 1 && m.warnings.some((w) => /polyhedral cells are not converted/.test(w)));
  m = await imp('old.cgns', cat(enc('@(#)ADF Database Version B02012>'), new Uint8Array(300))); ok('legacy ADF CGNS stays metadata-only with the conversion route', m.kind === 'metadata-only' && m.meta.container === 'ADF' && m.warnings.some((w) => /cgnsconvert/.test(w)), JSON.stringify(m.meta) + m.warnings.join('|'));
  m = await imp('cube.cgns', unstructured, { wasm: false }); ok('kernel switched off → metadata-only, no throw', m.kind === 'metadata-only' && m.warnings.some((w) => /switched off/.test(w)));
}

sec('MED (Salome)');
{
  const med = h5file((f) => {
    const info = f.create_group('INFOS_GENERALES'); info.create_attribute('MAJ', new Int32Array([4])); info.create_attribute('MIN', new Int32Array([1])); info.create_attribute('REL', new Int32Array([0]));
    const mesh = f.create_group('ENS_MAA').create_group('wingbox'); mesh.create_attribute('DIM', new Int32Array([3])); mesh.create_attribute('ESP', new Int32Array([3])); mesh.create_attribute('UNI', 'mm              mm              mm              '); mesh.create_attribute('DES', 'test mesh');
    const step = mesh.create_group('-0000000000000000001-0000000000000000001'), noe = step.create_group('NOE'), mai = step.create_group('MAI');
    noe.create_dataset({ name: 'COO', data: Float64Array.from([0, 1, 2].flatMap((k) => CV.map((v) => v[k]))) });
    const te = mai.create_group('TE4'), medTets = TETS.map((t) => [t[1], t[0], t[2], t[3]]);              // MED handedness: opposite to the platform convention
    te.create_dataset({ name: 'NOD', data: Int32Array.from([0, 1, 2, 3].flatMap((k) => medTets.map((t) => t[k] + 1))) }); te.create_dataset({ name: 'FAM', data: new Int32Array(6).fill(-1) });
    const tr = mai.create_group('TR3'); tr.create_dataset({ name: 'NOD', data: Int32Array.from([0, 1, 2].flatMap((k) => BT.map((t) => t[k] + 1))) }); tr.create_dataset({ name: 'FAM', data: Int32Array.from(BT, (_, i) => (i < 4 ? -2 : 0)) });
    mai.create_group('PO1').create_dataset({ name: 'NOD', data: new Int32Array([1]) });
    const fam = f.create_group('FAS').create_group('wingbox').create_group('ELEME');
    for (const [key, numId, gname] of [['FAM_-1_body', -1, 'body'], ['FAM_-2_root', -2, 'root rib']]) { const g = fam.create_group(key); g.create_attribute('NUM', new Int32Array([numId])); g.create_group('GRO').create_dataset({ name: 'NOM', data: chars(gname, 80) }); }
  });
  const m = await imp('wingbox.med', keep('wingbox.med', med)); fmtId('med', m, 'med', 'partial'); await cubeCheck('TE4 + TR3 (axis-by-axis coordinates, node-by-node connectivity)', m, { nv: 8 });
  ok('elements + kind', m.elements.find((e) => e.type === 'tet4')?.count === 6 && m.elements.find((e) => e.type === 'tri3')?.count === 12 && m.kind === 'structural-mesh');
  ok('families → group names', m.groups.find((g) => g.name === 'body')?.count === 6 && m.groups.find((g) => g.name === 'root rib')?.count === 4, JSON.stringify(m.groups));
  ok('units from UNI, version, description', m.units.length === 'mm' && m.meta.medVersion === '4.1.0' && m.meta.meshes[0].description === 'test mesh' && m.meta.meshes[0].elements.TE4 === 6, JSON.stringify(m.meta.meshes) + m.meta.medVersion);
  ok('MED handedness converted and logged', G.meshQuality(m).perType.find((p) => p.type === 'tet4').negJacobian === 0 && !G.analyse(m).inwardNormals && m.warnings.some((w) => /renumbered/.test(w)) && m.meta.reorderedBlocks[0] === '6 tet4');
  ok('unmapped type reported', m.warnings.some((w) => /PO1/.test(w)));
}

sec('Exodus II');
{
  // minimal classic NetCDF (CDF-1) writer
  const cdf = ({ dims, gatts = {}, vars }) => {
    const TY = { char: 2, int: 4, double: 6 }, SZ = { char: 1, int: 4, double: 8 }, chunks = []; const u32 = (v) => { const b = new Uint8Array(4); new DataView(b.buffer).setUint32(0, v); chunks.push(b); };
    const name = (s) => { const b = enc(s); u32(b.length); chunks.push(b, new Uint8Array((4 - (b.length % 4)) % 4)); };
    const atts = (o) => { const k = Object.keys(o); if (!k.length) { u32(0); u32(0); return; } u32(0x0c); u32(k.length); for (const n of k) { name(n); const v = o[n]; if (typeof v === 'string') { const b = enc(v); u32(2); u32(b.length); chunks.push(b, new Uint8Array((4 - (b.length % 4)) % 4)); } else { u32(6); u32(1); const b = new Uint8Array(8); new DataView(b.buffer).setFloat64(0, v); chunks.push(b); } } };
    const dimNames = dims.map((d) => d[0]), size = (v) => v.dims.reduce((s, d) => s * dims[dimNames.indexOf(d)][1], 1) * SZ[v.type], pad4 = (n) => n + ((4 - (n % 4)) % 4);
    const header = (begins) => { chunks.length = 0; chunks.push(enc('CDF'), Uint8Array.of(1)); u32(0); u32(0x0a); u32(dims.length); for (const [n, s] of dims) { name(n); u32(s); } atts(gatts); u32(0x0b); u32(vars.length); vars.forEach((v, i) => { name(v.name); u32(v.dims.length); for (const d of v.dims) u32(dimNames.indexOf(d)); atts(v.atts || {}); u32(TY[v.type]); u32(pad4(size(v))); u32(begins[i]); }); return cat(...chunks); };
    const h0 = header(vars.map(() => 0)), begins = []; let o = h0.length; for (const v of vars) { begins.push(o); o += pad4(size(v)); }
    const out = new Uint8Array(o); out.set(header(begins)); const dv = new DataView(out.buffer);
    vars.forEach((v, i) => { let p = begins[i]; if (v.type === 'char') out.set(enc(v.data), p); else for (const x of v.data) { if (v.type === 'int') dv.setInt32(p, x); else dv.setFloat64(p, x); p += SZ[v.type]; } });
    return out;
  };
  const exo = cdf({
    dims: [['len_name', 33], ['time_step', 0], ['num_dim', 3], ['num_nodes', 8], ['num_elem', 1], ['num_el_blk', 1], ['num_el_in_blk1', 1], ['num_nod_per_el1', 8], ['num_side_sets', 2], ['num_side_ss1', 4], ['num_side_ss2', 2], ['num_node_sets', 1], ['num_nod_ns1', 4]],
    gatts: { title: 'cube database', api_version: 8.11, version: 8.11 },
    vars: [
      { name: 'coordx', dims: ['num_nodes'], type: 'double', data: CV.map((v) => v[0]) }, { name: 'coordy', dims: ['num_nodes'], type: 'double', data: CV.map((v) => v[1]) }, { name: 'coordz', dims: ['num_nodes'], type: 'double', data: CV.map((v) => v[2]) },
      { name: 'eb_prop1', dims: ['num_el_blk'], type: 'int', data: [10], atts: { name: 'ID' } }, { name: 'eb_names', dims: ['num_el_blk', 'len_name'], type: 'char', data: 'core' },
      { name: 'connect1', dims: ['num_el_in_blk1', 'num_nod_per_el1'], type: 'int', data: [1, 2, 3, 4, 5, 6, 7, 8], atts: { elem_type: 'HEX8' } },
      { name: 'ss_prop1', dims: ['num_side_sets'], type: 'int', data: [5, 6] }, { name: 'ss_names', dims: ['num_side_sets', 'len_name'], type: 'char', data: 'sides'.padEnd(33, '\0') + 'caps' },
      { name: 'elem_ss1', dims: ['num_side_ss1'], type: 'int', data: [1, 1, 1, 1] }, { name: 'side_ss1', dims: ['num_side_ss1'], type: 'int', data: [1, 2, 3, 4] }, { name: 'elem_ss2', dims: ['num_side_ss2'], type: 'int', data: [1, 1] }, { name: 'side_ss2', dims: ['num_side_ss2'], type: 'int', data: [5, 6] },
      { name: 'ns_prop1', dims: ['num_node_sets'], type: 'int', data: [7] }, { name: 'ns_names', dims: ['num_node_sets', 'len_name'], type: 'char', data: 'clamped' }, { name: 'node_ns1', dims: ['num_nod_ns1'], type: 'int', data: [1, 2, 3, 4] },
    ],
  });
  ok('classic NetCDF detected as Exodus', G.detectFormat('cube.exo', exo).id === 'exodus');
  let m = await imp('cube.exo', keep('cube.exo', exo)); fmtId('classic NetCDF', m, 'exodus', 'partial'); await cubeCheck('HEX8 block', m, { nv: 8 });
  ok('block id + name → zone', m.groups.find((g) => g.kind === 'zone')?.name === 'core' && m.groups.find((g) => g.kind === 'zone').tag === 10 && m.elements.find((e) => e.type === 'hex8')?.count === 1 && m.kind === 'volume-mesh', JSON.stringify(m.groups));
  ok('side sets → boundary faces with names (all six sides distinct)', m.groups.find((g) => g.name === 'sides')?.count === 4 && m.groups.find((g) => g.name === 'caps')?.count === 2 && m.elements.find((e) => e.type === 'quad4')?.count === 6 && m.triangles.length === 36, JSON.stringify(m.groups));
  const caps = m.elements.find((e) => e.type === 'quad4'), zs = []; for (let e = 0; e < 6; e++) if (m.groups[caps.group[e]].name === 'caps') zs.push([0, 1, 2, 3].map((k) => m.positions[3 * caps.conn[4 * e + k] + 2]).join(''));
  ok('side numbering: sides 5 and 6 of a HEX are the bottom and top faces', zs.sort().join() === '0000,1111', zs.join());
  ok('node set census + header', m.groups.find((g) => g.name === 'clamped')?.nodes === 4 && m.meta.title === 'cube database' && m.meta.container === 'NetCDF classic (CDF-1)' && m.meta.blocks[0].elemType === 'HEX8' && m.units.length === null);
  const exo4 = h5file((f) => {
    f.create_attribute('title', 'tet cube'); for (const [n, s] of [['num_dim', 3], ['num_nodes', 8], ['num_elem', 6], ['num_el_blk', 1]]) f.create_dataset({ name: n, data: new Float32Array(s) });
    ['x', 'y', 'z'].forEach((ax, k) => f.create_dataset({ name: 'coord' + ax, data: Float64Array.from(CV, (v) => v[k]) }));
    f.create_dataset({ name: 'connect1', data: Int32Array.from(TETS.flat(), (v) => v + 1), shape: [6, 4] }).create_attribute('elem_type', 'TETRA'); f.create_dataset({ name: 'eb_prop1', data: new Int32Array([3]) }); f.create_dataset({ name: 'eb_names', data: ['fuel tank'], dtype: 'S33' });
  });
  m = await imp('tets.e', keep('tets.e', exo4)); fmtId('NetCDF-4', m, 'exodus', 'partial'); await cubeCheck('NetCDF-4 / HDF5 TETRA block', m, { nv: 8 });
  ok('HDF5 variant: block, name, title', m.elements[0].type === 'tet4' && m.elements[0].count === 6 && m.groups[0].name === 'fuel tank' && m.meta.title === 'tet cube' && m.meta.container === 'NetCDF-4 / HDF5', JSON.stringify(m.groups) + m.meta.title);
  m = await imp('plain.nc', cdf({ dims: [['x', 2]], vars: [{ name: 'temperature', dims: ['x'], type: 'double', data: [1, 2] }] })); ok('NetCDF that is not Exodus → metadata-only with reason', m.kind === 'metadata-only' && m.warnings.some((w) => /not an Exodus II database/.test(w)));
}

sec('Fluent CFF (.msh.h5)');
{
  const cff = h5file((f) => {
    const mesh = f.create_group('meshes').create_group('1'); mesh.create_attribute('dimension', new Int32Array([3])); mesh.create_attribute('nodeCount', new Int32Array([8])); mesh.create_attribute('faceCount', new Int32Array([7])); mesh.create_attribute('cellCount', new Int32Array([1]));
    mesh.create_group('nodes').create_group('coords').create_dataset({ name: '1', data: Float64Array.from(CV.flat()), shape: [8, 3] });
    const faces = mesh.create_group('faces'), zt = faces.create_group('zoneTopology');
    zt.create_dataset({ name: 'id', data: new Int32Array([3, 4, 9]) }); zt.create_dataset({ name: 'minId', data: new Int32Array([1, 5, 7]) }); zt.create_dataset({ name: 'maxId', data: new Int32Array([4, 6, 7]) }); zt.create_dataset({ name: 'zoneType', data: new Int32Array([3, 10, 2]) }); zt.create_dataset({ name: 'name', data: 'wing-wall;inlet;interior-fluid' });
    const sec1 = faces.create_group('nodes').create_group('1'); sec1.create_dataset({ name: 'nnodes', data: new Uint8Array([4, 4, 4, 4, 4, 4, 3]) }); sec1.create_dataset({ name: 'nodes', data: Uint32Array.from([...CQ.flat(), 0, 2, 6], (v) => v + 1) });
    faces.create_group('c0').create_dataset({ name: '1', data: new Uint32Array(7).fill(1) }); faces.create_group('c1').create_dataset({ name: '1', data: new Uint32Array([0, 0, 0, 0, 0, 0, 1]) });
    const cz = mesh.create_group('cells').create_group('zoneTopology'); cz.create_dataset({ name: 'id', data: new Int32Array([2]) }); cz.create_dataset({ name: 'minId', data: new Int32Array([1]) }); cz.create_dataset({ name: 'maxId', data: new Int32Array([1]) }); cz.create_dataset({ name: 'name', data: 'fluid' });
  });
  const m = await imp('case.msh.h5', keep('case.msh.h5', cff)); fmtId('cff', m, 'fluent-h5', 'partial'); await cubeCheck('boundary faces', m, { nv: 8 });
  ok('face zones → boundary groups; internal face left out', m.groups.find((g) => g.name === 'wing-wall')?.count === 4 && m.groups.find((g) => g.name === 'inlet')?.count === 2 && m.meta.internalFaces === 1 && m.groups.find((g) => g.kind === 'zone')?.name === 'fluid', JSON.stringify(m.groups));
  ok('honest: surface only', m.kind === 'surface-mesh' && m.meta.cells === 1 && m.warnings.some((w) => /NOT reconstructed/.test(w)));
  const other = h5file((f) => { f.create_group('results').create_dataset({ name: 'pressure', data: new Float64Array([1, 2, 3]) }); });
  const g = await imp('data.h5', other); ok('unknown HDF5 schema → object tree listed, metadata-only', g.format === 'hdf5' && g.kind === 'metadata-only' && g.meta.tree.some((l) => /\/results\/pressure\s+\[3\]/.test(l)), JSON.stringify(g.meta.tree));
}

// ---------- GeoTIFF ----------
sec('GeoTIFF terrain');
{
  /** TIFF LZW encoder (MSB-first codes, early change) for the fixtures. */
  const lzwEncode = (src) => {
    const out = []; let buf = 0, cnt = 0, bits = 9, next = 258, table = new Map();
    const put = (code) => { buf = (buf << bits) | code; cnt += bits; while (cnt >= 8) { out.push((buf >> (cnt - 8)) & 255); cnt -= 8; } buf &= (1 << cnt) - 1; };
    put(256); let w = -1;
    for (const c of src) {
      if (w < 0) { w = c; continue; }
      const key = w * 256 + c, hit = table.get(key);
      if (hit !== undefined) { w = hit; continue; }
      put(w); table.set(key, next++); if (next > (1 << bits) - 1 && bits < 12) bits++;
      if (next >= 4094) { put(256); table = new Map(); bits = 9; next = 258; }
      w = c;
    }
    if (w >= 0) put(w); put(257); if (cnt) out.push((buf << (8 - cnt)) & 255);
    return Uint8Array.from(out);
  };
  /** Write a single-band (or RGB) TIFF. data: typed array of W·H·S samples. */
  const tiff = ({ W, H, data, fmt = 3, le = true, comp = 1, pred = 1, tile = 0, rps = 8, S = 1, photo = 1, scale = null, tie = null, keys = null, nodata = null }) => {
    const Bs = data.BYTES_PER_ELEMENT, sample = (dv, p, v) => { if (fmt === 3) { if (Bs === 4) dv.setFloat32(p, v, le); else dv.setFloat64(p, v, le); } else if (Bs === 1) dv.setUint8(p, v & 255); else if (Bs === 2) dv.setUint16(p, v & 65535, le); else dv.setUint32(p, v >>> 0, le); };
    const block = (r0, c0, bh, bw) => {
      const raw = new Uint8Array(bh * bw * S * Bs), dv = new DataView(raw.buffer);
      for (let r = 0; r < bh; r++) for (let c = 0; c < bw; c++) for (let s = 0; s < S; s++) { const rr = r0 + r, cc = c0 + c; sample(dv, ((r * bw + c) * S + s) * Bs, rr < H && cc < W ? data[(rr * W + cc) * S + s] : 0); }
      if (pred === 2) for (let r = 0; r < bh; r++) for (let i = bw * S - 1; i >= S; i--) { const p = (r * bw * S + i) * Bs; if (Bs === 1) raw[p] = (raw[p] - raw[p - S]) & 255; else if (Bs === 2) dv.setUint16(p, (dv.getUint16(p, le) - dv.getUint16(p - 2 * S, le)) & 65535, le); else dv.setUint32(p, (dv.getUint32(p, le) - dv.getUint32(p - 4 * S, le)) >>> 0, le); }
      if (pred === 3) for (let r = 0; r < bh; r++) { const n = bw * S, row = new Uint8Array(n * Bs), o = r * n * Bs; for (let i = 0; i < n; i++) for (let q = 0; q < Bs; q++) row[q * n + i] = raw[o + i * Bs + (le ? Bs - 1 - q : q)]; for (let i = n * Bs - 1; i >= S; i--) row[i] = (row[i] - row[i - S]) & 255; raw.set(row, o); }
      return comp === 5 ? lzwEncode(raw) : comp === 8 ? new Uint8Array(zlib.deflateSync(raw)) : raw;
    };
    const blocks = []; if (tile) for (let r = 0; r < H; r += tile) for (let c = 0; c < W; c += tile) blocks.push(block(r, c, tile, tile)); else for (let r = 0; r < H; r += rps) blocks.push(block(r, 0, Math.min(rps, H - r), W));
    const entries = [], extra = []; let extraLen = 0;
    const E = (tag, type, vals) => entries.push({ tag, type, vals });
    E(256, 4, [W]); E(257, 4, [H]); E(258, 3, new Array(S).fill(Bs * 8)); E(259, 3, [comp]); E(262, 3, [photo]); E(277, 3, [S]); E(339, 3, new Array(S).fill(fmt)); if (pred !== 1) E(317, 3, [pred]);
    if (tile) { E(322, 4, [tile]); E(323, 4, [tile]); E(324, 4, blocks.map(() => 0)); E(325, 4, blocks.map((b) => b.length)); } else { E(278, 4, [rps]); E(273, 4, blocks.map(() => 0)); E(279, 4, blocks.map((b) => b.length)); }
    if (scale) E(33550, 12, scale); if (tie) E(33922, 12, tie); if (keys) E(34735, 3, [1, 1, 0, keys.length, ...keys.flatMap(([k, v]) => [k, 0, 1, v])]); if (nodata !== null) E(42113, 2, [...enc(String(nodata) + '\0')]);
    entries.sort((a, b) => a.tag - b.tag);
    const SZ = { 2: 1, 3: 2, 4: 4, 12: 8 }, ifdAt = 8, ifdLen = 2 + 12 * entries.length + 4; let dataAt = ifdAt + ifdLen;
    for (const e of entries) { const n = e.vals.length * SZ[e.type]; if (n > 4) { e.at = dataAt + extraLen; extraLen += n + (n % 2); } }
    let off = dataAt + extraLen; const offs = blocks.map((b) => { const o = off; off += b.length; return o; });
    for (const e of entries) if (e.tag === 273 || e.tag === 324) e.vals = offs;
    const out = new Uint8Array(off), dv = new DataView(out.buffer); out.set(le ? [0x49, 0x49] : [0x4d, 0x4d]); dv.setUint16(2, 42, le); dv.setUint32(4, ifdAt, le); dv.setUint16(ifdAt, entries.length, le);
    const W1 = (p, type, v) => { if (type === 3) dv.setUint16(p, v, le); else if (type === 4) dv.setUint32(p, v, le); else if (type === 12) dv.setFloat64(p, v, le); else dv.setUint8(p, v); };
    entries.forEach((e, i) => { const p = ifdAt + 2 + 12 * i; dv.setUint16(p, e.tag, le); dv.setUint16(p + 2, e.type, le); dv.setUint32(p + 4, e.vals.length, le); if (e.at) { dv.setUint32(p + 8, e.at, le); e.vals.forEach((v, k) => W1(e.at + k * SZ[e.type], e.type, v)); } else e.vals.forEach((v, k) => W1(p + 8 + k * SZ[e.type], e.type, v)); });
    blocks.forEach((b, i) => out.set(b, offs[i])); void extra;
    return out;
  };
  const W = 41, H = 31, z = (c, r) => 100 + 2 * c + 3 * r, plane = Float32Array.from({ length: W * H }, (_, i) => z(i % W, Math.floor(i / W)));
  const geo = { scale: [10, 10, 0], tie: [0, 0, 0, 500000, 4100000, 0], keys: [[1024, 1], [1025, 2], [3072, 32633], [3076, 9001]] };
  const exactArea = 400 * 300 * Math.sqrt(1 + 0.2 * 0.2 + 0.3 * 0.3);
  const check = async (label, bytes, { units = 'm' } = {}) => {
    const m = await imp('dem.tif', bytes), a = G.analyse(m);
    ok(`${label}: height field read`, m.format === 'geotiff' && m.support === 'partial' && m.kind === 'surface' && a.nVerts === W * H && a.nTris === 2 * (W - 1) * (H - 1), `${m.kind} ${a.nVerts}/${a.nTris} ${m.warnings.join(' | ')}`);
    ok(`${label}: georeferenced extents`, Math.abs(a.bbox.min[0] - 500000) < 1e-6 && Math.abs(a.bbox.max[0] - 500400) < 1e-6 && Math.abs(a.bbox.max[1] - 4100000) < 1e-6 && Math.abs(a.bbox.min[1] - 4099700) < 1e-6, JSON.stringify(a.bbox));
    ok(`${label}: elevations exact`, Math.abs(a.bbox.min[2] - 100) < 1e-9 && Math.abs(a.bbox.max[2] - z(W - 1, H - 1)) < 1e-9 && m.meta.zRange[0] === 100); near(`${label}: surface area of the tilted plane`, a.area, exactArea, 1e-9);
    ok(`${label}: units`, m.units.length === units, JSON.stringify(m.units)); return m;
  };
  let m = await check('float32 strips, uncompressed', keep('dem.tif', tiff({ W, H, data: plane, ...geo })));
  ok('georeferencing in meta', m.meta.geoTiff && m.meta.pixelScale[0] === 10 && m.meta.tiePoints[3] === 500000 && m.meta.geoKeys.projectedCRS === 32633 && m.meta.geoKeys.linearUnits === 9001 && m.meta.width === W && m.meta.compression === 'none', JSON.stringify(m.meta.geoKeys));
  ok('upward-facing surface', G.massProperties(m, 1).closed === false && (() => { const P = m.positions, T = m.triangles, u = [P[3 * T[1]] - P[3 * T[0]], P[3 * T[1] + 1] - P[3 * T[0] + 1]], v = [P[3 * T[2]] - P[3 * T[0]], P[3 * T[2] + 1] - P[3 * T[0] + 1]]; return u[0] * v[1] - u[1] * v[0] > 0; })());
  ok('warned that the band is taken as elevation', m.warnings.some((w) => /interpreted as elevation/.test(w)));
  await check('uint16 LZW + horizontal predictor', keep('dem_lzw.tif', tiff({ W, H, data: Uint16Array.from(plane), fmt: 1, comp: 5, pred: 2, ...geo })));
  await check('float32 Deflate + floating-point predictor', keep('dem_def.tif', tiff({ W, H, data: plane, comp: 8, pred: 3, ...geo })));
  await check('float64 LZW, big-endian', tiff({ W, H, data: Float64Array.from(plane), comp: 5, le: false, ...geo }));
  await check('int16 tiles (16 × 16), big-endian', keep('dem_tiled.tif', tiff({ W, H, data: Int16Array.from(plane), fmt: 2, tile: 16, le: false, ...geo })));
  await check('uint32 tiles, Deflate + predictor', tiff({ W, H, data: Uint32Array.from(plane), fmt: 1, tile: 16, comp: 8, pred: 2, ...geo }));
  await check('uint8-range LZW strips of one row', tiff({ W, H, data: plane, comp: 5, rps: 1, ...geo }));
  await check('feet', tiff({ W, H, data: plane, ...geo, keys: [[1024, 1], [1025, 2], [3076, 9002]] }), { units: 'ft' });
  const noisy = Uint16Array.from({ length: 200 * 150 }, (_, i) => (i * 7919) % 3001); m = await imp('noise.tif', tiff({ W: 200, H: 150, data: noisy, fmt: 1, comp: 5, rps: 150 }), { maxGrid: 2000 });
  ok('LZW with code-width growth and table resets decodes exactly', m.positions.length === 3 * 200 * 150 && (() => { for (let i = 0; i < noisy.length; i++) if (m.positions[3 * i + 2] !== noisy[i]) return false; return true; })());
  ok('no georeferencing → pixel coordinates + warning, units null', m.units.length === null && m.meta.georeferenced === false && m.warnings.some((w) => /pixel indices/.test(w)) && G.analyse(m).bbox.size[0] === 199);
  const holed = plane.slice(); holed[5 * W + 5] = -9999; m = await imp('hole.tif', tiff({ W, H, data: holed, nodata: -9999, ...geo })); ok('no-data sample left out with its four cells', m.positions.length / 3 === W * H - 1 && m.triangles.length / 3 === 2 * (W - 1) * (H - 1) - 8 && m.meta.noData === -9999 && m.meta.noDataSamples === 1 && m.meta.zRange[0] === 100);
  m = await imp('big.tif', tiff({ W: 400, H: 300, data: Float32Array.from({ length: 120000 }, (_, i) => (i % 400) * 0.5), comp: 8, rps: 16, ...geo })); let a = G.analyse(m);
  ok('large raster subsampled to opts.maxGrid (default 300)', m.meta.grid.stride === 2 && m.meta.grid.nx === 200 && m.meta.grid.ny === 150 && a.nVerts === 30000 && m.warnings.some((w) => /sampled every 2/.test(w)) && Math.abs(a.bbox.max[2] - 199) < 1e-9, JSON.stringify(m.meta.grid));
  m = await imp('big.tif', tiff({ W: 400, H: 300, data: new Float32Array(120000), ...geo }), { maxGrid: 50 }); ok('opts.maxGrid respected', m.meta.grid.nx <= 50 && m.meta.grid.ny <= 50);
  m = await imp('geo.tif', tiff({ W, H, data: plane, scale: [0.001, 0.001, 0], tie: [0, 0, 0, 12.5, 45.2, 0], keys: [[1024, 2], [2048, 4326]] })); ok('geographic CRS → degrees warning, units not set', m.units.length === null && m.warnings.some((w) => /degrees/.test(w)) && m.kind === 'surface');
  m = await imp('photo.tif', tiff({ W: 8, H: 8, data: new Uint8Array(192), fmt: 1, S: 3, photo: 2 })); ok('RGB imagery is not mistaken for terrain', m.kind === 'metadata-only' && m.meta.width === 8 && m.warnings.some((w) => /imagery/.test(w)) && m.warnings.some((w) => /Conversion pathway/.test(w)));
  const jpeg = tiff({ W, H, data: plane, ...geo }); for (let i = 10; i < 10 + 12 * 14; i += 12) if (new DataView(jpeg.buffer).getUint16(i, true) === 259) new DataView(jpeg.buffer).setUint16(i + 8, 7, true);
  m = await imp('jpeg.tif', jpeg); ok('unsupported compression → metadata-only naming it', m.kind === 'metadata-only' && m.meta.compression === 'JPEG' && m.warnings.some((w) => /JPEG/.test(w)));
}

sec('E57 header + XML');
{
  const xml = enc(`<?xml version="1.0"?><e57Root type="Structure"><data3D type="Vector"><vectorChild type="Structure"><name type="String"><![CDATA[Hangar scan 1]]></name><cartesianBounds type="Structure"><xMinimum type="Float">-2.5</xMinimum><xMaximum type="Float">7.5</xMaximum><yMinimum type="Float">0</yMinimum><yMaximum type="Float">4</yMaximum><zMinimum type="Float">-1</zMinimum><zMaximum type="Float">3</zMaximum></cartesianBounds><points type="CompressedVector" fileOffset="48" recordCount="123456"><prototype type="Structure"><cartesianX type="Float"/><cartesianY type="Float"/><cartesianZ type="Float"/><intensity type="Float"/></prototype></points></vectorChild></data3D>${' '.repeat(700)}</e57Root>`);
  const start = 700, b = new Uint8Array(4096), d = new DataView(b.buffer); b.set(enc('ASTM-E57')); d.setUint32(8, 1, true); d.setBigUint64(16, 4096n, true); d.setBigUint64(24, BigInt(start), true); d.setBigUint64(32, BigInt(xml.length), true); d.setBigUint64(40, 1024n, true);
  for (let p = start, o = 0; o < xml.length;) { const end = (Math.floor(p / 1024) + 1) * 1024 - 4, k = Math.min(end - p, xml.length - o); b.set(xml.subarray(o, o + k), p); o += k; b.fill(0xee, end, end + 4); p = end + 4; }   // page checksums interrupt the XML
  const m = await imp('scan.e57', keep('scan.e57', b)); fmtId('e57', m, 'e57', 'metadata');
  ok('XML section read across page checksums', m.meta.xmlRead && m.meta.scans.length === 1 && m.meta.scans[0].points === 123456 && m.meta.scans[0].name === 'Hangar scan 1' && m.meta.totalPoints === 123456, JSON.stringify(m.meta.scans));
  ok('bounds + fields', m.meta.cartesianBounds.min.join() === '-2.5,0,-1' && m.meta.cartesianBounds.max.join() === '7.5,4,3' && m.meta.fields.includes('cartesianX') && m.meta.fields.includes('intensity')); ok('still honest: no points', m.kind === 'metadata-only' && m.positions.length === 0 && m.warnings.some((w) => /not decoded/.test(w)));
}

sec('simplifyForSolver');
{
  const sph = uvSphere(500, 96, 48); sph.units = { length: 'mm', source: 'file' };
  const s = G.simplifyForSolver(sph, 3000), nT = s.triangles.length / 3;
  ok('plain arrays', Array.isArray(s.positions) && Array.isArray(s.triangles) && s.positions.length % 3 === 0 && s.triangles.length % 3 === 0); ok('triangle budget met and well used', nT <= 3000 && nT > 1200 && sph.triangles.length / 3 === 9024, `${nT} triangles`);
  ok('indices valid, no degenerate triangles', s.triangles.every((v) => Number.isInteger(v) && v >= 0 && v < s.positions.length / 3) && (() => { for (let t = 0; t < s.triangles.length; t += 3) if (s.triangles[t] === s.triangles[t + 1] || s.triangles[t + 1] === s.triangles[t + 2] || s.triangles[t] === s.triangles[t + 2]) return false; return true; })());
  const sm = mk([], []); sm.positions = Float64Array.from(s.positions); sm.triangles = Uint32Array.from(s.triangles); const a = G.analyse(sm);
  ok('converted to metres', s.inMetres === true && s.scaleApplied === 1e-3 && a.bbox.size.every((d) => Math.abs(d - 1) < 0.06), JSON.stringify(a.bbox.size)); near('shape kept: sphere area within 5 %', a.area, Math.PI, 5e-2);
  ok('closed sphere stays closed', a.watertight && Math.abs(a.volume - Math.PI / 6) < 0.05 * Math.PI / 6, `${a.boundaryEdges} boundary edges, V = ${a.volume}`);
  const tight = G.simplifyForSolver(sph, 200); ok('smaller budget respected', tight.triangles.length / 3 <= 200 && tight.triangles.length / 3 > 40, String(tight.triangles.length / 3));
  const c = G.simplifyForSolver(cube(), 3000); ok('small model passes through (welded), units unknown → unscaled', c.triangles.length === 36 && c.positions.length === 24 && c.inMetres === false && c.scaleApplied === 1);
  const stl = await imp('c.stl', G.exportModel(cube(), 'stl')); ok('STL soup is welded on the way', G.simplifyForSolver(stl).positions.length === 24);
  ok('volume mesh → its boundary surface', G.simplifyForSolver(tetCube(), 100).triangles.length === 36); ok('nothing to simplify → empty arrays', G.simplifyForSolver(mk(CV, [], { kind: 'pointcloud' })).triangles.length === 0);
  const def = G.simplifyForSolver(uvSphere(1, 200, 100)); ok('default budget is 3000', def.triangles.length / 3 <= 3000 && def.sourceTriangles === 39600);
}

// ---------- metadata-only signatures ----------
sec('metadata-only formats');
{
  const pad = (b, n = 256) => cat(b, new Uint8Array(Math.max(0, n - b.length)).fill(7));
  const hdf5 = pad(Uint8Array.from([0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a, 2]));
  const e57 = (() => { const b = new Uint8Array(96), d = new DataView(b.buffer); b.set(enc('ASTM-E57')); d.setUint32(8, 1, true); d.setBigUint64(16, 96n, true); d.setBigUint64(24, 48n, true); d.setBigUint64(32, 10n, true); d.setBigUint64(40, 1024n, true); return b; })();
  const tif = (() => { const b = new Uint8Array(200), d = new DataView(b.buffer); b.set([0x49, 0x49, 0x2a, 0]); d.setUint32(4, 8, true); d.setUint16(8, 4, true); const E = (i, tag, type, cnt, val) => { const o = 10 + 12 * i; d.setUint16(o, tag, true); d.setUint16(o + 2, type, true); d.setUint32(o + 4, cnt, true); if (type === 3 && cnt === 1) d.setUint16(o + 8, val, true); else d.setUint32(o + 8, val, true); }; E(0, 256, 4, 1, 640); E(1, 257, 4, 1, 480); E(2, 33550, 12, 3, 100); E(3, 34735, 3, 8, 140); [30, 30, 0].forEach((v, i) => d.setFloat64(100 + 8 * i, v, true)); [1, 1, 0, 1, 3076, 0, 1, 9001].forEach((v, i) => d.setUint16(140 + 2 * i, v, true)); return b; })();
  const sat = '700 0 1 0\n@7 TestCAD @11 ACIS 7.0 NT @24 Wed May 01 10:00:00 2024\n25.4 9.9999999999999995e-07 1e-10\nbody $-1 $1 $-1 $-1 #\nlump $-1 $-1 $2 $0 #\nplane-surface $-1 0 0 0 0 0 1 #\nface $3 $-1 #\nface $4 $-1 #\nEnd-of-ACIS-data\n';
  const cases = [
    ['mesh.cgns', hdf5, 'cgns'], ['model.med', hdf5, 'med'], ['case.msh.h5', hdf5, 'fluent-h5'], ['data.h5', hdf5, 'hdf5'], ['out.exo', hdf5, 'exodus'],
    ['out.e', pad(cat(enc('CDF'), Uint8Array.from([1, 0, 0, 0, 5]))), 'exodus'], ['old.cgns', pad(enc('@(#)ADF Database Version B02012>')), 'cgns'],
    ['part.jt', pad(enc('Version 9.5 JT  DM 4.0 \n')), 'jt'], ['part.x_b', pad(cat(enc('PS'), new Uint8Array(6), enc('?: TRANSMIT FILE'))), 'xb'],
    ['part.sat', enc(sat), 'acis'], ['part.sab', pad(enc('ACIS BinaryFile\x04\0\0\0')), 'acis'], ['scan.e57', e57, 'e57'], ['dem.tif', tif, 'geotiff'],
    ['flow.plt', pad(enc('#!TDV112')), 'tecplot-bin'], ['flow.szplt', pad(enc('#!SZPLT 1.0')), 'tecplot-bin'],
    ['wing.CATPart', pad(enc('V5_CFV2\0\0')), 'nativecad'], ['wing.sldprt', pad(Uint8Array.from([0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1])), 'nativecad'], ['part.prt', pad(Uint8Array.from([1, 2, 3, 0, 0, 5, 200, 201])), 'nativecad'],
  ];
  for (const [name, bytes, id] of cases) {
    truncatable.push([name, bytes]);
    const d = G.detectFormat(name, bytes); ok(`${name} → ${id}`, d.id === id, `got ${d.id} (${d.note})`);
    const m = await imp(name, bytes);
    const decodable = ['cgns', 'med', 'fluent-h5', 'exodus', 'geotiff'].includes(id);    // formats with a real reader: a bare signature still yields no geometry
    ok(`${name}: metadata-only model, never pretends`, m.format === id && m.kind === 'metadata-only' && (decodable || m.support === 'metadata') && m.positions.length === 0 && m.triangles.length === 0 && m.units !== undefined, `${m.format}/${m.kind}/${m.support}`);
    ok(`${name}: warning names the conversion pathway`, m.warnings.some((w) => /Conversion pathway: \S+/.test(w)), m.warnings.join(' | '));
    const a = G.analyse(m); ok(`${name}: analyse is safe on it`, a.nVerts === 0 && a.classification.label.includes('metadata'));
  }
  ok('HDF5 signature wins over a misleading extension', G.detectFormat('mesh.stl', hdf5).id === 'hdf5');
  let m = await imp('part.sat', sat); ok('SAT: header + census + units from the scale', m.meta.product === 'TestCAD' && m.meta.acisVersion === 'ACIS 7.0 NT' && m.meta.census.face === 2 && m.meta.census['plane-surface'] === 1 && m.units.length === 'in', JSON.stringify(m.meta));
  m = await imp('dem.tif', tif); ok('TIFF without strip data: size + pixel scale + GeoKeys still reported', m.kind === 'metadata-only' && m.meta.width === 640 && m.meta.height === 480 && m.meta.geoTiff && m.meta.pixelScale[0] === 30 && m.meta.geoKeys.linearUnits === 9001, JSON.stringify(m.meta));
  m = await imp('scan.e57', e57); ok('E57 header', m.meta.version === '1.0' && m.meta.xmlOffset === 48 && m.meta.pageSize === 1024 && m.meta.scans.length === 0);
  m = await imp('mesh.cgns', hdf5); ok('corrupt HDF5: container + superblock version reported, reason given', m.meta.container === 'HDF5' && m.meta.superblockVersion === 2 && m.warnings.some((w) => /could not open/.test(w)), JSON.stringify(m.meta) + m.warnings.join('|'));
  m = await imp('part.jt', pad(enc('Version 9.5 JT  DM 4.0 \n'))); ok('JT version', m.meta.version === '9.5');
}

// ---------- detection ----------
sec('detectFormat');
{
  const d = (n, s) => G.detectFormat(n, typeof s === 'string' ? enc(s) : s);
  const gm = '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n', fl = '(0 "fluent")\n(2 3)\n(10 (0 1 8 0 3))\n', cdbm = '/COM,ANSYS RELEASE\n/PREP7\nNBLOCK,6,SOLID,8,8\n';
  ok('.msh Gmsh', d('a.msh', gm).id === 'gmsh'); ok('.msh Fluent', d('a.msh', fl).id === 'fluent'); ok('.msh holding ANSYS archive commands', d('a.msh', cdbm).id === 'cdb');
  ok('.dat Nastran', d('a.dat', 'BEGIN BULK\nGRID,1,,0.,0.,0.\n').id === 'nastran'); ok('.dat Nastran with TITLE in case control', d('a.dat', 'SOL 101\nCEND\nTITLE = wing\nBEGIN BULK\nGRID           1              0.      0.      0.\n').id === 'nastran');
  ok('.dat Tecplot', d('a.dat', 'TITLE = "x"\nVARIABLES = "X" "Y"\nZONE N=3, E=1\n').id === 'tecplot'); ok('.dat numeric table', d('a.dat', '0.0 0.0 0.0\n1.0 0.0 0.0\n1.0 1.0 0.0\n').id === 'xyz');
  ok('.xyz point cloud', d('a.xyz', '0.5 0.25 0.1\n1.5 0.25 0.1\n1.5 1.0 0.3\n0.1 0.2 0.3\n').id === 'xyz');
  ok('.xyz Plot3D single block', d('a.xyz', '2 2 2\n0 1 0 1 0 1 0 1\n0 0 1 1 0 0 1 1\n0 0 0 0 1 1 1 1\n').id === 'plot3d'); ok('.xyz Plot3D multi-block', d('a.xyz', '1\n2 2 1\n0 1 0 1\n0 0 1 1\n0 0 0 0\n').id === 'plot3d');
  ok('.xyz integer point rows are not Plot3D', d('a.xyz', '2 2 2\n3 3 3\n4 4 4\n5 5 5\n6 6 6\n').id === 'xyz', d('a.xyz', '2 2 2\n3 3 3\n4 4 4\n5 5 5\n6 6 6\n').id);
  ok('.vtk legacy', d('a.vtk', '# vtk DataFile Version 3.0\n').id === 'vtk'); ok('.vtk holding XML', d('a.vtk', '<?xml version="1.0"?>\n<VTKFile type="PolyData">').id === 'vtkxml');
  ok('.stl ASCII', d('a.stl', 'solid x\nfacet normal 0 0 1\n').id === 'stl'); ok('.stl holding OBJ', d('a.stl', 'v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n').id === 'obj');
  ok('content beats extension, with a note', /content decides/.test(d('a.stl', 'v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n').note));
  ok('.inp Abaqus', d('a.inp', '*HEADING\nx\n*NODE\n').id === 'abaqus'); ok('.unv', d('a.unv', '    -1\n  2411\n').id === 'unv'); ok('.su2', d('a.su2', '% c\nNDIME= 3\n').id === 'su2');
  ok('.stp', d('a.stp', 'ISO-10303-21;\nHEADER;').id === 'step'); ok('.cdb', d('a.cdb', cdbm).id === 'cdb'); ok('glTF JSON', d('a.gltf', '{"asset":{"version":"2.0"}}').id === 'gltf');
  ok('.ply', d('a.ply', 'ply\nformat ascii 1.0\n').id === 'ply'); ok('.off', d('a.off', 'OFF\n8 6 12\n').id === 'off'); ok('.pcd', d('a.pcd', '# .PCD v0.7\nVERSION 0.7\nFIELDS x y z\n').id === 'pcd');
  ok('.q Plot3D solution by extension, low confidence', d('a.q', Uint8Array.from([0, 0, 0, 12, 1, 2, 3, 0, 0, 0])).id === 'plot3dq' && d('a.q', Uint8Array.from([0, 0, 0, 12, 1, 2, 3, 0, 0, 0])).confidence < 0.5);
  const lo = d('a.stp', 'hello world this is not step'); ok('extension-only match has low confidence and says so', lo.id === 'step' && lo.confidence <= 0.3 && /extension only/.test(lo.note), JSON.stringify(lo));
  ok('unknown → unknown', d('a.qqq', 'hello world').id === 'unknown' && d('a.qqq', 'hello world').confidence === 0); ok('empty → unknown', d('a.stl', new Uint8Array(0)).id === 'unknown');
  ok('confidence in 0..1', [d('a.msh', gm), d('a.dat', '0 0 0\n1 1 1\n'), lo].every((r) => r.confidence >= 0 && r.confidence <= 1 && typeof r.note === 'string'));
  let threw = 0; for (const [n, b] of [['a.qqq', enc('hello world')], ['a.stl', new Uint8Array(0)]]) try { await G.importFile(n, b); } catch { threw++; }
  ok('importFile throws only for unknown / empty input', threw === 2);
  const forced = await imp('weird.bin', G.exportModel(cube(), 'stl'), { format: 'stl' }); ok('opts.format forces a reader', forced.format === 'stl' && forced.triangles.length === 36);
}

// ---------- analysis ----------
sec('heal');
{
  // cube as triangle soup with jittered duplicates, one flipped facet, one degenerate facet and one stray vertex
  const V = [], T = []; let s = 7; const r = () => ((s = (s * 16807) % 2147483647) / 2147483647 - 0.5) * 2e-9;
  CT.forEach((t, i) => { const base = V.length; for (const v of t) V.push(CV[v].map((c) => c + r())); T.push(i === 3 ? [base, base + 2, base + 1] : [base, base + 1, base + 2]); });
  V.push([0, 0, 0], [1e-9, 0, 0], [1, 1, 1]); T.push([36, 37, 38]); V.push([9, 9, 9]);
  const dirty = mk(V, T), before = G.analyse(dirty), snapshot = Array.from(dirty.positions).join() + '|' + Array.from(dirty.triangles).join();
  ok('dirty model is not watertight before healing', before.watertight === false && before.boundaryEdges > 0);
  const h = G.heal(dirty, { weldTol: 1e-6 }), a = G.analyse(h.model), step = (k) => h.log.find((l) => l.step === k)?.detail ?? '';
  ok('input untouched', Array.from(dirty.positions).join() + '|' + Array.from(dirty.triangles).join() === snapshot && dirty.log.length === 0); ok('changed flag', h.changed === true);
  ok('weld count: 40 → 8 + stray', /^tolerance 1\.000e-6 .*: 31 vertex/.test(step('weld')), step('weld'));
  const moved = +(/moved = ([\d.e+-]+)/.exec(step('weld'))?.[1] ?? NaN); ok('max vertex movement reported and below tolerance', moved > 0 && moved <= 1e-6, String(moved));
  ok('degenerate count', /^1 degenerate .* 0 duplicate/.test(step('degenerate')), step('degenerate')); ok('orientation count', /^1 triangle\(s\) reversed/.test(step('orientation')), step('orientation'));
  ok('unreferenced count', /^1 unreferenced/.test(step('unreferenced')), step('unreferenced')); ok('dimension change reported', /largest change/.test(step('dimensions')));
  ok('healed: 8 vertices, 12 triangles', a.nVerts === 8 && a.nTris === 12, `${a.nVerts}/${a.nTris}`); ok('healed: watertight, no conflicts', a.watertight && a.orientationConflicts === 0 && a.duplicateVerts === 0 && a.degenerateTris === 0);
  near('healed: volume', a.volume, 1, 1e-6); near('healed: area', a.area, 6, 1e-6); ok('healed: outward normals', a.inwardNormals === false);
  ok('provenance appended to the new model', h.model.log.filter((l) => l.step.startsWith('heal:')).length === h.log.length);
  const inside = cube(); for (let i = 0; i < inside.triangles.length; i += 3) { const b = inside.triangles[i + 1]; inside.triangles[i + 1] = inside.triangles[i + 2]; inside.triangles[i + 2] = b; }
  ok('inside-out cube flagged by analyse', G.analyse(inside).inwardNormals === true && G.analyse(inside).volume === 1);
  const h2 = G.heal(inside); ok('inside-out cube turned outward', G.analyse(h2.model).inwardNormals === false && /12 triangle/.test(h2.log.find((l) => l.step === 'orientation').detail) && /inside-out/.test(h2.log.find((l) => l.step === 'orientation').detail));
  const h3 = G.heal(cube()); ok('clean model: nothing changed', h3.changed === false && /no change/.test(h3.log.at(-1).detail));
  const stl = await imp('c.stl', G.exportModel(cube(), 'stl')), h4 = G.heal(stl); ok('STL soup welds to 8 shared vertices', h4.model.positions.length === 24 && /28 vertex/.test(h4.log[0].detail) && /= 0\.000e\+0/.test(h4.log[0].detail), h4.log[0].detail);
  const h5 = G.heal(tetCube()); ok('volume mesh: elements kept, boundary recomputed', h5.model.elements[0].count === 6 && h5.model.triangles.length === 36 && /skipped/.test(h5.log.find((l) => l.step === 'orientation').detail));
  const pc = mk(CV, [], { kind: 'pointcloud' }), h6 = G.heal(pc); ok('point cloud keeps its points', h6.model.positions.length === 24);
}

sec('transform');
{
  const c = cube(), t = G.transform(c, { scale: 0.001, translate: [1, 0, 0], units: 'm' }), a = G.analyse(t);
  ok('input untouched + log appended', c.positions[3] === 1 && c.log.length === 0 && t.log.at(-1).step === 'transform'); ok('units recorded', t.units.length === 'm' && t.units.source === 'user');
  near('scaled volume', a.volume, 1e-9, 1e-9); near('translated', a.bbox.min[0], 1, 1e-12);
  const y = G.analyse(G.transform(mk([[0, 2, 0], [1, 2, 0], [0, 2, 3]], [[0, 1, 2]]), { swapYZ: true })); ok('swapYZ: Y-up → Z-up', Math.abs(y.bbox.max[2] - 2) < 1e-12 && Math.abs(y.bbox.min[1] + 3) < 1e-12, JSON.stringify(y.bbox));
  const f = G.analyse(G.transform(cube(), { flipX: true })); ok('flipX keeps outward orientation', f.inwardNormals === false && f.volume === 1 && f.bbox.min[0] === -1);
  const ft = G.transform(tetCube(), { flipX: true }); ok('flipX keeps positive tets', G.meshQuality(ft).perType[0].negJacobian === 0);
  const r = G.analyse(G.transform(cube(), { rotateDeg: [0, 0, 90] })); ok('rotate 90° about z', Math.abs(r.bbox.min[0] + 1) < 1e-12 && Math.abs(r.bbox.max[1] - 1) < 1e-12 && Math.abs(r.volume - 1) < 1e-12);
}

sec('massProperties');
{
  const c = G.massProperties(G.transform(cube(), { translate: [5, -3, 2] }), 2700);
  near('cube volume', c.volume, 1, 1e-12); near('cube mass', c.mass, 2700, 1e-12); ok('cube cg', Math.abs(c.cg[0] - 5.5) < 1e-12 && Math.abs(c.cg[1] + 2.5) < 1e-12 && Math.abs(c.cg[2] - 2.5) < 1e-12, c.cg.join()); ok('cube closed', c.closed === true);
  near('cube Ixx = m/6', c.inertia[0][0], 2700 / 6, 1e-10); near('cube Izz = m/6', c.inertia[2][2], 2700 / 6, 1e-10); ok('cube products of inertia vanish', Math.abs(c.inertia[0][1]) < 1e-9 && Math.abs(c.inertia[1][2]) < 1e-9);
  const R = 2, s = G.massProperties(uvSphere(R, 96, 48), 1), Vex = (4 / 3) * Math.PI * R ** 3;
  near('sphere volume → 4/3 π r³ (tessellation error)', s.volume, Vex, 5e-3); ok('sphere volume below exact (inscribed)', s.volume < Vex);
  near('sphere inertia = 2/5 m r²', s.inertia[0][0], 0.4 * s.mass * R * R, 5e-3); near('sphere inertia isotropic', s.inertia[2][2], s.inertia[0][0], 5e-3); ok('sphere cg at origin', Math.hypot(...s.cg) < 1e-9);
  const coarse = G.massProperties(uvSphere(R, 24, 12), 1); ok('refinement reduces the volume error', Math.abs(coarse.volume - Vex) > Math.abs(s.volume - Vex));
  const open = mk(CV, CT.slice(0, 10)); ok('open surface reported as not closed', G.massProperties(open, 1).closed === false && G.analyse(open).volume === null && G.analyse(open).boundaryEdges === 4);
  const as = G.analyse(uvSphere(1, 48, 24)); ok('sphere: one component, watertight, generic part', as.components === 1 && as.watertight && as.classification.label === 'generic part', as.classification.label);
}

sec('slice');
{
  const shoelace = (l) => { let s = 0; for (let i = 0; i + 1 < l.length; i++) s += l[i][0] * l[i + 1][1] - l[i + 1][0] * l[i][1]; return Math.abs(s) / 2; };
  let L = G.slice(cube(), { axis: 'z', value: 0.3 }); ok('cube: one closed loop', L.length === 1 && L[0][0][0] === L[0].at(-1)[0] && L[0][0][1] === L[0].at(-1)[1], JSON.stringify(L.map((l) => l.length))); near('cube: loop area', shoelace(L[0]), 1, 1e-12);
  L = G.slice(cube(), { axis: 'x', value: 0.5 }); near('cube x-plane area', shoelace(L[0]), 1, 1e-12);
  L = G.slice(cube(), { axis: 'y', value: 0 }); ok('plane through a face: still a closed loop', L.length === 1 && Math.abs(shoelace(L[0]) - 1) < 1e-12, JSON.stringify(L));
  ok('plane that misses → no loops', G.slice(cube(), { axis: 'z', value: 2 }).length === 0);
  const stl = await imp('c.stl', G.exportModel(cube(), 'stl')); L = G.slice(stl, { axis: 'z', value: 0.5 }); ok('unwelded STL soup still chains into one loop', L.length === 1 && Math.abs(shoelace(L[0]) - 1) < 1e-12);
  const cyl = cylinder(0.5, 2, 180); L = G.slice(cyl, { axis: 'z', value: 1 }); ok('cylinder: one loop', L.length === 1); near('cylinder: section area → π r²', shoelace(L[0]), Math.PI * 0.25, 1e-3);
  L = G.slice(cyl, { axis: 'x', value: 0.1 }); near('cylinder lengthwise section area', shoelace(L[0]), 2 * 2 * Math.sqrt(0.25 - 0.01), 1e-3);
  const two = mk([...CV, ...CV.map((v) => [v[0] + 3, v[1], v[2]])], [...CT, ...CT.map((t) => t.map((v) => v + 8))]); ok('two bodies → two loops', G.slice(two, { axis: 'z', value: 0.5 }).length === 2 && G.analyse(two).components === 2);
  const openS = G.slice(mk(CV, CT.slice(2)), { axis: 'x', value: 0.5 }); ok('open surface → open polyline', openS.length === 1 && (openS[0][0][0] !== openS[0].at(-1)[0] || openS[0][0][1] !== openS[0].at(-1)[1]));
  // extruded NACA 0012 wing → slice → metrics
  const prof = naca(0.12, 60), n = prof.length, WV = [], WT = [];
  for (const y of [0, 4]) for (const p of prof) WV.push([p[0], y, p[1]]);
  for (let i = 0; i < n; i++) { const j = (i + 1) % n; WT.push([i, j, n + j], [i, n + j, n + i]); }
  for (let i = 1; i + 1 < n; i++) { WT.push([0, i + 1, i]); WT.push([n, n + i, n + i + 1]); }
  const wing = G.heal(mk(WV, WT)).model, aw = G.analyse(wing); ok('extruded wing is watertight', aw.watertight, `${aw.boundaryEdges} boundary edges`);
  ok('wing classified as a lifting surface', aw.classification.label === 'wing / lifting surface' && aw.classification.confidence > 0 && aw.classification.confidence <= 1 && aw.classification.why.length > 10, JSON.stringify(aw.classification));
  ok('dims', Math.abs(aw.dims.length_x - 1) < 1e-9 && Math.abs(aw.dims.span_y - 4) < 1e-9 && Math.abs(aw.dims.height_z - 0.12) < 2e-3);
  const sm = G.sectionMetrics(G.slice(wing, { axis: 'y', value: 2 })[0]); near('sliced wing section t/c', sm.tc, 0.12, 5e-3); near('sliced wing chord', sm.chord, 1, 1e-6); ok('slice loop closed', sm.closed === true);
}

sec('sectionMetrics');
{
  const m = G.sectionMetrics(naca(0.12, 200));
  near('NACA 0012 t/c', m.tc, 0.12, 2e-3); near('chord', m.chord, 1, 1e-9); near('thickness', m.thickness, 0.12, 2e-3); ok('max thickness near 30 % chord', Math.abs(m.xThickness - 0.3) < 0.02, String(m.xThickness));
  ok('symmetric: camber ≈ 0', Math.abs(m.camber) < 1e-6); ok('LE / TE', Math.abs(m.le[0]) < 1e-9 && Math.abs(m.te[0] - 1) < 1e-9); near('area ≈ 0.685 t c²', m.area, 0.0822, 1e-2); ok('perimeter slightly above 2 chords', m.perimeter > 2 && m.perimeter < 2.1);
  const c4 = G.sectionMetrics(naca(0.12, 200, 0.04)); near('cambered: max camber 4 %', c4.camberRatio, 0.04, 2e-2); near('cambered: t/c unchanged', c4.tc, 0.12, 5e-3);
  const cw = G.sectionMetrics(naca(0.12, 200).reverse()); near('clockwise input gives the same t/c', cw.tc, m.tc, 1e-12);
  const closedLoop = [...naca(0.12, 50), naca(0.12, 50)[0]]; ok('explicitly closed loop flagged closed', G.sectionMetrics(closedLoop).closed === true && m.closed === false);
  ok('degenerate input does not throw', G.sectionMetrics([]).chord === 0 && G.sectionMetrics([[0, 0]]).chord === 0 && G.sectionMetrics([[0, 0], [1, 0]]).chord === 1);
}

sec('meshSection + sectionProperties');
{
  const sq = [[0, 0], [1, 0], [1, 1], [0, 1]];
  let m = G.meshSection([sq], { h: 0.04 }), p = G.sectionProperties(m);
  ok('square mesh quality', m.quality.minAngle > 20 && m.quality.meanAspect < 1.3 && m.quality.nTris === m.tris.length && m.quality.areaError < 1e-9, JSON.stringify(m.quality));
  ok('mesh shape', m.nodes.every((q) => q.length === 2) && m.tris.every((t) => t.length === 3 && t.every((v) => v >= 0 && v < m.nodes.length)));
  near('square A', p.A, 1, 1e-9); near('square centroid', p.cx, 0.5, 1e-9); near('square Ixx = 1/12', p.Ixx, 1 / 12, 1e-9); near('square Iyy', p.Iyy, 1 / 12, 1e-9); ok('square Ixy = 0', Math.abs(p.Ixy) < 1e-9);
  near('square J ≈ 0.1406 a⁴', p.J, 0.140577, 6e-3); ok('linear FE converges from below', p.J < 0.140577); ok('dof + note', p.dof > 100 && /Prandtl/.test(p.note));
  const pcw = G.sectionProperties(G.meshSection([sq.slice().reverse()], { h: 0.04 })); near('clockwise input: same J', pcw.J, p.J, 1e-9);
  const circ = Array.from({ length: 360 }, (_, i) => [2 + 0.5 * Math.cos((i * Math.PI) / 180), -1 + 0.5 * Math.sin((i * Math.PI) / 180)]);
  m = G.meshSection([circ], { h: 0.02 }); p = G.sectionProperties(m);
  near('circle J = π r⁴/2', p.J, (Math.PI * 0.5 ** 4) / 2, 5e-3); near('circle A', p.A, Math.PI * 0.25, 1e-3); near('circle Ixx = π r⁴/4', p.Ixx, (Math.PI * 0.5 ** 4) / 4, 2e-3); ok('circle centroid (offset loop)', Math.abs(p.cx - 2) < 1e-6 && Math.abs(p.cy + 1) < 1e-6); near('circle I1 = I2', p.I1, p.I2, 1e-3);
  const ea = 1, eb = 0.5, ell = Array.from({ length: 400 }, (_, i) => [ea * Math.cos((i * Math.PI) / 200), eb * Math.sin((i * Math.PI) / 200)]);
  p = G.sectionProperties(G.meshSection([ell], { h: 0.025 })); near('ellipse J = π a³b³/(a²+b²)', p.J, (Math.PI * ea ** 3 * eb ** 3) / (ea * ea + eb * eb), 5e-3); near('ellipse Ixx = π a b³/4', p.Ixx, (Math.PI * ea * eb ** 3) / 4, 2e-3); near('ellipse Iyy = π a³ b/4', p.Iyy, (Math.PI * ea ** 3 * eb) / 4, 2e-3);
  ok('ellipse principal axes', Math.abs(p.I1 - p.Iyy) < 1e-6 && Math.abs(p.I2 - p.Ixx) < 1e-6 && Math.abs(Math.abs(p.theta_p) - Math.PI / 2) < 1e-3, `theta ${p.theta_p}`);
  const b = 2, t = 0.1; m = G.meshSection([[[0, 0], [b, 0], [b, t], [0, t]]], { h: t / 6 }); p = G.sectionProperties(m);
  near('thin rectangle J → b t³/3 (1 − 0.63 t/b)', p.J, ((b * t ** 3) / 3) * (1 - 0.63 * t / b), 1.5e-2); near('thin rectangle A', p.A, b * t, 1e-9);
  const rot = [[0, 0], [b, 0], [b, t], [0, t]].map(([x, y]) => [x * Math.cos(0.5) - y * Math.sin(0.5), x * Math.sin(0.5) + y * Math.cos(0.5)]), pr = G.sectionProperties(G.meshSection([rot], { h: t / 4 }));
  near('rotated strip: principal angle (axis of I1 is across the strip, modulo π)', ((pr.theta_p % Math.PI) + Math.PI) % Math.PI, 0.5 + Math.PI / 2, 1e-6); near('rotated strip: I2 = b t³/12', pr.I2, (b * t ** 3) / 12, 1e-6); near('rotated strip: I1 = t b³/12', pr.I1, (t * b ** 3) / 12, 1e-6);
  const ro = 1, ri = 0.6, C = (r) => Array.from({ length: 360 }, (_, i) => [r * Math.cos((i * Math.PI) / 180), r * Math.sin((i * Math.PI) / 180)]);
  m = G.meshSection([C(ro), C(ri).reverse()], { h: 0.03 }); p = G.sectionProperties(m);
  near('tube A (hole respected)', p.A, Math.PI * (ro * ro - ri * ri), 1e-3); near('tube J = π/2 (ro⁴ − ri⁴)', p.J, (Math.PI / 2) * (ro ** 4 - ri ** 4), 5e-3); ok('hollow treatment stated in note', p.holes === 1 && /Hollow section/.test(p.note));
  p = G.sectionProperties(G.meshSection([naca(0.12, 100)], { h: 0.008 })); near('aerofoil section area from the mesh', p.A, G.sectionMetrics(naca(0.12, 100)).area, 2e-3); ok('aerofoil J positive and below polar moment', p.J > 0 && p.J < p.Ixx + p.Iyy);
  ok('single loop accepted without wrapping', G.meshSection(sq, { h: 0.2 }).tris.length > 10); ok('default h', G.meshSection([sq]).tris.length > 100);
  ok('degenerate input → empty mesh, no throw', G.meshSection([], {}).tris.length === 0 && G.meshSection([[[0, 0], [1, 1]]], {}).tris.length === 0 && G.sectionProperties({ nodes: [], tris: [] }).A === 0);
  const coarseThin = G.meshSection([[[0, 0], [b, 0], [b, t], [0, t]]], { h: 0.5 }); ok('too-coarse h on a thin strip still returns a valid mesh', coarseThin.tris.length > 0 && Number.isFinite(coarseThin.quality.areaError));
}

sec('sectionConvergence');
{
  const sq = [[0, 0], [1, 0], [1, 1], [0, 1]], c = G.sectionConvergence([sq], [0.1, 0.025, 0.05]);
  ok('rows sorted fine → coarse', c.rows.length === 3 && c.rows[0].h === 0.025 && c.rows[2].h === 0.1 && c.rows[0].nTris > c.rows[1].nTris && c.rows[1].nTris > c.rows[2].nTris);
  ok('J increases monotonically with refinement', c.rows[0].J > c.rows[1].J && c.rows[1].J > c.rows[2].J && c.gci.monotonic);
  ok('observed order ≈ 2', c.gci.p > 1.6 && c.gci.p < 2.4, `p = ${c.gci.p}`); near('Richardson-extrapolated J', c.gci.fExact, 0.140577, 1e-3); ok('GCI small and finite', c.gci.gciFine > 0 && c.gci.gciFine < 0.02);
  ok('areas exact on every grid', c.rows.every((r) => Math.abs(r.A - 1) < 1e-9)); ok('fewer than three sizes → no GCI', G.sectionConvergence([sq], [0.1, 0.2]).gci === null);
}

sec('meshQuality');
{
  const s3 = Math.sqrt(3) / 2, eq = mk([[0, 0, 0], [1, 0, 0], [0.5, s3, 0], [1.5, s3, 0]], [[0, 1, 2], [1, 3, 2]]), q = G.meshQuality(eq);
  ok('equilateral: good', q.grade === 'good' && q.issues.length === 0, JSON.stringify(q.issues)); near('equilateral aspect = 1', q.worst.aspect, 1, 1e-9); ok('equilateral skew = 0, angles 60°', q.worst.skew < 1e-9 && Math.abs(q.worst.minAngle - 60) < 1e-9 && Math.abs(q.perType[0].maxAngle.max - 60) < 1e-9);
  ok('works on bare triangles + stat shape', q.perType.length === 1 && q.perType[0].type === 'tri3' && q.perType[0].count === 2 && q.perType[0].aspect.hist.centers.length === q.perType[0].aspect.hist.counts.length && q.perType[0].skew.mean >= 0);
  const sl = G.meshQuality(mk([[0, 0, 0], [1, 0, 0], [0.5, 0.004, 0], [0.5, s3, 0]], [[0, 1, 2], [0, 1, 3]]));
  ok('sliver: poor or unusable, with issues', ['poor', 'unusable'].includes(sl.grade) && sl.issues.length > 0 && sl.worst.minAngle < 1 && sl.worst.skew > 0.98 && sl.worst.aspect > 50, JSON.stringify(sl.worst));
  const tq = G.meshQuality(tetCube()); ok('Kuhn tets: no negative Jacobians, acceptable', tq.perType[0].type === 'tet4' && tq.perType[0].negJacobian === 0 && tq.perType[0].count === 6 && ['good', 'acceptable'].includes(tq.grade), tq.grade + JSON.stringify(tq.worst));
  const inv = tetCube(); [inv.elements[0].conn[0], inv.elements[0].conn[1]] = [inv.elements[0].conn[1], inv.elements[0].conn[0]];
  const iq = G.meshQuality(inv); ok('one inverted tet → unusable + issue', iq.perType[0].negJacobian === 1 && iq.grade === 'unusable' && iq.issues.some((s) => /inverted/.test(s)));
  const hexM = mk(CV, [], { kind: 'volume-mesh', elements: [{ type: 'hex8', nodesPer: 8, count: 1, conn: Uint32Array.from([0, 1, 2, 3, 4, 5, 6, 7]) }] }), hq = G.meshQuality(hexM);
  ok('unit hex: perfect', hq.grade === 'good' && hq.worst.skew < 1e-12 && hq.worst.aspect === 1 && hq.perType[0].negJacobian === 0 && hq.worst.minAngle === 90);
  const reg = mk([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], [], { elements: [{ type: 'tet4', nodesPer: 4, count: 1, conn: Uint32Array.from([0, 1, 2, 3]) }] }), rq = G.meshQuality(reg);
  near('regular tet aspect = 1', rq.worst.aspect, 1, 1e-9); ok('regular tet skew 0', rq.worst.skew < 1e-9);
  ok('nothing to evaluate → unusable with reason', G.meshQuality(mk(CV, [], { kind: 'pointcloud' })).grade === 'unusable');
  const an = G.analyse(tetCube()); ok('tet cube: analysed through its boundary', an.watertight && an.volume === 1 && an.nTris === 12 && ['structural mesh', 'CFD volume mesh'].includes(an.classification.label));
}

sec('exportModel');
{
  const tc = tetCube();
  for (const [fmt, name] of [['stl', 'e.stl'], ['obj', 'e.obj'], ['vtk', 'e.vtk'], ['msh', 'e.msh'], ['su2', 'e.su2'], ['json', 'e.json']]) {
    const txt = G.exportModel(cube(), fmt); ok(`${fmt}: returns a string`, typeof txt === 'string' && txt.length > 50);
    const m = await imp(name, txt); await cubeCheck(`${fmt} surface round trip`, m, { nv: fmt === 'stl' ? 36 : 8 });
    if (fmt !== 'stl' && fmt !== 'obj') { const v = await imp(name, G.exportModel(tc, fmt)); await cubeCheck(`${fmt} volume round trip`, v, { nv: 8 }); ok(`${fmt}: tets preserved`, v.elements.find((e) => e.type === 'tet4')?.count === 6, v.elements.map((e) => e.type + e.count).join()); }
  }
  const j = await imp('e.json', G.exportModel({ ...tc, units: { length: 'mm', source: 'file' } }, 'json')); ok('json keeps units and groups', j.units.length === 'mm' && j.groups[0].name === 'solid body' && j.elements[0].group[5] === 0 && j.format === 'json');
  let threw = false; try { G.exportModel(cube(), 'dwg'); } catch { threw = true; } ok('unknown export format rejected', threw);
  const precise = mk([[0.1 + 0.2, 1 / 3, 1e-17], [1, 0, 0], [0, 1, 0]], [[0, 1, 2]]), back = await imp('p.vtk', G.exportModel(precise, 'vtk')); ok('coordinates survive at full precision', back.positions[0] === 0.1 + 0.2 && back.positions[1] === 1 / 3 && back.positions[2] === 1e-17);
}

// ---------- robustness ----------
sec('malformed input');
{
  let cases = 0, bad = 0; const errs = [];
  const tryOne = async (name, bytes, opts) => {
    cases++;
    try { const m = await G.importFile(name, bytes, { kernelTimeout: 15000, ...opts }); G.analyse(m); G.meshQuality(m); if (!Array.isArray(m.warnings) || !Array.isArray(m.log) || !(m.positions instanceof Float64Array) || !(m.triangles instanceof Uint32Array) || m.positions.length % 3 || m.triangles.some((v) => v >= m.positions.length / 3)) throw new Error('inconsistent model'); }
    catch (e) { if (!/not recognised|is empty/.test(e.message)) { bad++; if (errs.length < 8) errs.push(`${name} (${bytes.length} B): ${e.constructor.name}: ${e.message}`); } }
  };
  let seed = 99; const rnd = () => ((seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff);
  for (const [name, b] of truncatable) {
    for (const f of [0.999, 0.9, 0.75, 0.5, 0.3, 0.12, 0.03]) await tryOne(name, b.slice(0, Math.max(1, Math.floor(b.length * f))));
    for (let k = 0; k < 4; k++) { const c = b.slice(); for (let i = 0; i < 1 + c.length / 60; i++) c[Math.floor(rnd() * c.length)] = Math.floor(rnd() * 256); await tryOne(name, c); }      // random byte corruption
    const c = b.slice(); for (let i = 0; i < c.length; i++) if (c[i] >= 49 && c[i] <= 56 && rnd() < 0.15) c[i] = 57; await tryOne(name, c);                                                    // digits inflated to 9 (bad counts / indices)
  }
  for (const id of G.FORMATS.map((f) => f.id)) for (const junk of [enc('garbage in, garbage out\n1 2 3\n'), Uint8Array.from({ length: 300 }, (_, i) => (i * 37) & 255), enc('9999999999 9999999999 9999999999\n'), new Uint8Array(3)]) await tryOne('junk.bin', junk, { format: id });
  const huge = new Uint8Array(84); new DataView(huge.buffer).setUint32(80, 0xffffffff, true); await tryOne('huge.stl', huge);
  await tryOne('huge.ply', enc('ply\nformat binary_little_endian 1.0\nelement vertex 4000000000\nproperty float x\nproperty float y\nproperty float z\nend_header\n'));
  await tryOne('huge.msh', enc('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n99999999999\n1 0 0 0\n$EndNodes\n'));
  await tryOne('huge.p3d', enc('100000 100000 100000\n0 0 0\n')); await tryOne('huge.off', enc('OFF\n4000000000 1 0\n0 0 0\n'));
  await tryOne('deep.x3d', enc('<X3D>' + '<a>'.repeat(5000)));
  await tryOne('bad.glb', cat(enc('glTF'), Uint8Array.from([2, 0, 0, 0, 255, 255, 255, 127, 255, 255, 255, 127, 0x4a, 0x53, 0x4f, 0x4e])));
  await tryOne('bad.3mf', cat(Uint8Array.from([0x50, 0x4b, 3, 4]), new Uint8Array(100)));
  ok(`${cases} truncated / corrupted / hostile inputs handled without an unhandled error`, bad === 0, `${bad} failures: ${errs.join(' || ')}`);
  const m = await imp('t.stl', G.exportModel(cube(), 'stl').slice(0, 700)); ok('truncated ASCII STL: partial facets dropped with a warning', m.triangles.length % 3 === 0 && m.triangles.length < 36);
  const big = 'solid big\n' + 'facet normal 0 0 1\n outer loop\n  vertex 0 0 0\n  vertex 1 0 0\n  vertex 0 1 0\n endloop\nendfacet\n'.repeat(200000) + 'endsolid big\n', t0 = Date.now(), mb = await imp('big.stl', big);
  ok(`large ASCII STL (${(big.length / 1e6).toFixed(0)} MB, 200k facets) parsed by scanning`, mb.triangles.length === 600000, `${mb.triangles.length / 3} facets in ${Date.now() - t0} ms`);
}

// ---------- summary ----------
console.log(`geometry: ${pass} passed, ${fail} failed`);
if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
