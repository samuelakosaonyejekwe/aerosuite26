// Readers backed by the vendored WebAssembly kernels: exact CAD geometry (STEP, IGES, OpenCASCADE BREP) is
// tessellated by OpenCASCADE, and LAZ point clouds are decoded by LASzip. The text parsers stay in charge of
// what they read reliably (header, schema, units, census); the kernel adds the faces. When a kernel cannot be
// loaded the readers fall back to the text-only result and say so.

import { Mesh, guard, str, LIMITS } from './parsers-util.js';
import { readSTEP, readIGES } from './parsers-cad.js';
import { lasHeader } from './parsers-surface.js';
import { loadOcct, loadLazPerf, dropKernel, nodeImport } from './parsers-wasm.js';

const OCCT_UNIT = { mm: 'millimeter', cm: 'centimeter', m: 'meter', in: 'inch', ft: 'foot' };

/** Tessellation parameters: opts.linearDeflection (ratio of the bounding box unless linearDeflectionType is 'absolute_value') and opts.angularDeflection (rad). */
function occtParams(opts, unit) {
  const lin = Number.isFinite(+opts.linearDeflection) && +opts.linearDeflection > 0 ? +opts.linearDeflection : 0.001, ang = Number.isFinite(+opts.angularDeflection) && +opts.angularDeflection > 0 ? +opts.angularDeflection : 0.5;
  return { linearUnit: OCCT_UNIT[unit] || 'millimeter', linearDeflectionType: opts.linearDeflectionType === 'absolute_value' ? 'absolute_value' : 'bounding_box_ratio', linearDeflection: lin, angularDeflection: ang };
}

/** Convert the kernel result ({ root, meshes: [{ name, color, brep_faces, position, index }] }) into element blocks: one group per body, face ranges and colours in the metadata. */
function occtToPart(res, bodyKind) {
  const M = new Mesh(), faces = [], faceOf = [], names = new Map();
  const walk = (node, path, depth) => {
    if (!node || depth > 200) return;
    const here = node.name ? (path ? `${path} / ${node.name}` : node.name) : path;
    (node.meshes || []).forEach((mi, k) => { if (!names.has(mi)) names.set(mi, node.meshes.length > 1 ? `${here || 'body'} · ${k + 1}` : here); });
    for (const c of node.children || []) walk(c, here, depth + 1);
  };
  walk(res.root, '', 0);
  let skipped = 0;
  (res.meshes || []).forEach((ms, mi) => {
    const P = ms.position, I = ms.index; if (!P || !I || !P.length || !I.length) { skipped++; return; }
    const base = M.nNodes, nv = Math.floor(P.length / 3), nt = Math.floor(I.length / 3);
    guard(base + nv, 'tessellation vertex', LIMITS.verts);
    for (let i = 0; i < nv; i++) M.node(P[3 * i], P[3 * i + 1], P[3 * i + 2]);
    const g = M.group(bodyKind, M.groups.length, names.get(mi) || ms.name || `body ${mi + 1}`), t0 = faceOf.length;
    M.groups[g].color = ms.color ?? null; M.groups[g].faces = (ms.brep_faces || []).length;
    const tf = new Int32Array(nt).fill(-1);
    for (const f of ms.brep_faces || []) { const id = faces.length; faces.push({ body: g, first: t0 + f.first, last: t0 + f.last, color: f.color ?? null }); for (let t = f.first; t <= f.last && t < nt; t++) tf[t] = id; }
    for (let t = 0; t < nt; t++) { const a = I[3 * t], b = I[3 * t + 1], c = I[3 * t + 2]; if (a < nv && b < nv && c < nv) { M.elem('tri3', [base + a, base + b, base + c], g); faceOf.push(tf[t]); } }
  });
  const out = M.result({ kind: 'surface' });
  if (out.elements[0]) out.elements[0].face = Int32Array.from(faceOf);     // B-rep face index per triangle (see meta.brepFaces)
  return { out, faces, skipped };
}

// ---- kernel worker with a watchdog -------------------------------------------------------------------------
// OpenCASCADE runs synchronously inside WebAssembly and cannot be interrupted, and a damaged file can keep it
// busy indefinitely. It therefore runs in a worker that is terminated when opts.kernelTimeout (ms) expires.
const IS_NODE = typeof process === 'object' && !!process?.versions?.node && typeof window === 'undefined' && typeof importScripts === 'undefined';
let worker = null, queue = Promise.resolve(), jobId = 0;
async function spawnWorker() {
  const u = new URL('./parsers-occt-worker.js', import.meta.url);
  if (IS_NODE) {
    const { Worker: W } = await nodeImport('worker_threads'), w = new W(u), h = { cb: null };
    w.on('message', (m) => h.cb?.(m)); w.on('error', (e) => h.cb?.({ ok: false, fatal: true, error: `kernel worker failed: ${e?.message || e}` })); w.unref();
    return { post: (m) => w.postMessage(m), kill: () => w.terminate(), busy: (on) => (on ? w.ref() : w.unref()), h };
  }
  const w = new Worker(u, { type: 'module' }), h = { cb: null };
  w.onmessage = (e) => h.cb?.(e.data); w.onerror = (e) => { e.preventDefault?.(); h.cb?.({ ok: false, fatal: true, error: `kernel worker failed: ${e?.message || 'it could not be started'}` }); };
  return { post: (m) => w.postMessage(m), kill: () => w.terminate(), busy() {}, h };
}
function inWorker(kind, bytes, params, timeout) {
  const job = queue.then(async () => {
    if (!worker) worker = await spawnWorker();
    const w = worker, id = ++jobId;
    return new Promise((resolve, reject) => {
      const done = (fn, v, kill) => { clearTimeout(timer); w.h.cb = null; if (kill) { w.kill(); if (worker === w) worker = null; } else w.busy(false); fn(v); };
      const timer = setTimeout(() => done(reject, new Error(`the OpenCASCADE kernel did not finish within ${Math.round(timeout / 1000)} s and was stopped (damaged file, or a model too heavy for the chosen deflection; raise opts.kernelTimeout to wait longer)`), true), timeout);
      w.h.cb = (m) => { if (m.id !== undefined && m.id !== id) return; if (m.ok) done(resolve, m, false); else done(reject, new Error(m.error), !!m.fatal); };
      w.busy(true); w.post({ id, kind, bytes, params });
    });
  });
  queue = job.catch(() => {});
  return job;
}
/** The BREP text reader of the kernel does not survive incomplete input: require every section and the closing root reference. */
function brepLooksComplete(b) {
  const s = str(b, 0, b.length), order = ['Locations', 'Curve2ds', 'Curves', 'Polygon3D', 'PolygonOnTriangulations', 'Surfaces', 'Triangulations', 'TShapes']; let at = 0;
  if (!/CASCADE Topology V\d/.test(s.slice(0, 400))) return false;
  const lines = s.split(/\r?\n/).length;
  for (const k of order) { const m = new RegExp(`^${k} (\\d+)\\s*$`, 'm').exec(s.slice(at)); if (!m || +m[1] > lines) return false; at += m.index + m[0].length; }
  const n = +/^TShapes (\d+)\s*$/m.exec(s)[1], body = s.slice(s.search(/^TShapes \d+\s*$/m));
  return n > 0 && (body.match(/^(Ve|Ed|Wi|Fa|Sh|So|CS|Co)\s*$/gm) || []).length === n && (body.match(/\*\s*$/gm) || []).length === n && /\n[+-]\d+ \d+\s*$/.test(s.trimEnd() + ' ');
}

/** Run the kernel on `bytes`; returns the tessellated part or throws with a reason. */
async function tessellate(kind, bytes, ctx, unit, bodyKind) {
  const params = occtParams(ctx.opts, unit), timeout = Math.max(1000, +ctx.opts.kernelTimeout || 180000);
  if (kind === 'brep' && !brepLooksComplete(bytes)) throw new Error('the BREP file is incomplete or damaged (a section or the closing shape reference is missing)');
  let res;
  try { res = await inWorker(kind, bytes, params, timeout); }
  catch (e) {
    if (!/worker failed|Worker is not|not defined|Invalid URL|Failed to construct/i.test(String(e?.message))) throw e;
    // no worker available in this runtime (e.g. a single-file build): translate on the calling thread, without a watchdog
    worker = null;
    const occt = await loadOcct(); let r;
    try { r = kind === 'step' ? occt.ReadStepFile(bytes, params) : kind === 'iges' ? occt.ReadIgesFile(bytes, params) : occt.ReadBrepFile(bytes, params); }
    catch (e2) { dropKernel('OpenCASCADE'); throw new Error(`the OpenCASCADE kernel failed on this file (${typeof e2 === 'number' || /exception/i.test(String(e2?.message)) ? 'internal exception' : e2?.message || e2})`); }
    if (!r || !r.success) throw new Error('the OpenCASCADE kernel could not translate this file');
    res = { root: r.root, meshes: (r.meshes || []).map((m) => ({ name: m.name ?? '', color: m.color ?? null, brep_faces: m.brep_faces || [], position: m.attributes?.position?.array || [], index: m.index?.array || [] })) };
  }
  const { out, faces, skipped } = occtToPart(res, bodyKind);
  if (!out.positions.length || !out.elements.length) throw new Error('the OpenCASCADE kernel found no faces to tessellate (wireframe-only or empty model)');
  if (skipped) ctx.warn(`${skipped} body/bodies produced no triangles and were skipped.`);
  return { out, faces, params };
}

/** Merge the text parser's metadata/units with the kernel's tessellation, or fall back to the text result. */
async function withKernel(kind, label, b, ctx, text) {
  const hold = [], base = text ? text(b, { ...ctx, warn: (m) => hold.push(m) }) : { kind: 'metadata-only', meta: {} };
  const flush = () => { for (const m of hold) ctx.warn(m); };
  if (base.kind === 'surface' && !(base.meta?.census?.ADVANCED_FACE || base.meta?.census?.MANIFOLD_SOLID_BREP)) { flush(); return base; }   // AP242 tessellation only: nothing for the kernel to add
  if (ctx.opts.wasm === false) { flush(); ctx.warn(`${label}: exact-geometry tessellation was switched off (opts.wasm = false); only the text-level content is shown.`); return { ...base, support: 'partial' }; }
  try {
    // bodies are called solids only when the file itself declares solid topology
    const census = base.meta?.census || {}, solid = kind === 'step' ? !!(census.MANIFOLD_SOLID_BREP || census.BREP_WITH_VOIDS || census.FACETED_BREP) : kind === 'iges' ? !!census['186 manifold solid B-rep object'] : !!base.meta?.hasSolids;
    const { out, faces, params } = await tessellate(kind, b, ctx, base.units?.length, solid ? 'solid' : 'surface');
    for (const m of hold) if (!/NOT tessellated|point cloud|indicative only/.test(m)) ctx.warn(m);
    const meta = { ...(base.meta || {}), kernel: 'OpenCASCADE (occt-import-js 0.0.23)', tessellation: { ...params, bodies: out.groups.length, faces: faces.length, triangles: out.elements[0].count }, brepFaces: faces.length <= 20000 ? faces : faces.slice(0, 20000) };
    if (faces.length > 20000) ctx.warn(`The model has ${faces.length} B-rep faces; only the first 20000 face ranges are listed in the metadata.`);
    if (!base.units?.length) ctx.warn(`${label}: the file does not state a length unit; coordinates are as written in the file (the kernel assumes millimetres when none is declared).`);
    return { ...out, units: base.units, meta };
  } catch (e) {
    flush();
    ctx.warn(`${label}: exact faces could not be tessellated — ${e.message}. Falling back to the text-level content (header, census, explicit points and polylines).`);
    return { ...base, support: 'partial' };
  }
}

export const readSTEPKernel = (b, ctx) => withKernel('step', 'STEP', b, ctx, readSTEP);
export const readIGESKernel = (b, ctx) => withKernel('iges', 'IGES', b, ctx, readIGES);
/** OpenCASCADE native BREP text: header version by text, geometry by the kernel. The format carries no units. */
export const readBREP = (b, ctx) => withKernel('brep', 'BREP', b, ctx, (bytes) => {
  const head = str(bytes, 0, Math.min(bytes.length, 400)), v = /CASCADE Topology V(\d+)/.exec(head);
  return { kind: 'metadata-only', meta: { topologyVersion: v ? +v[1] : null, drawable: /DBRep_DrawableShape/.test(head), hasSolids: bytes.length < 64e6 && /^So\s*$/m.test(str(bytes, 0, bytes.length)) } };
});

/** LAZ: header by the LAS reader, point records by LASzip; every stride-th point is kept up to opts.maxPoints. */
export async function readLAZ(b, ctx) {
  const h = lasHeader(b), meta = { ...h, compressed: true };
  if (ctx.opts.wasm === false) throw Object.assign(new Error('LAZ decoding was switched off (opts.wasm = false)'), { meta });
  let L; try { L = await loadLazPerf(); } catch (e) { throw Object.assign(e, { meta }); }
  const fp = L._malloc(b.length); if (!fp) throw Object.assign(new Error('not enough WebAssembly memory to hold the LAZ file'), { meta });
  let lz = null, pp = 0;
  try {
    L.HEAPU8.set(b, fp);
    lz = new L.LASZip(); lz.open(fp, b.length);
    const n = Math.min(guard(lz.getCount() || h.pointCount, 'LAZ point', 4e9), h.pointCount || Infinity), len = lz.getPointLength();
    if (!(len >= 12)) throw new Error('LAZ point record length is not plausible');
    const maxPoints = Math.max(1, Math.floor(ctx.opts.maxPoints ?? 200000)), stride = Math.max(1, Math.ceil(n / maxPoints)), m = Math.ceil(n / stride), pos = new Float64Array(3 * m);
    pp = L._malloc(len);
    for (let i = 0, k = 0; i < n; i++) {
      lz.getPoint(pp);                                    // records are entropy-coded in sequence, so every point is decoded
      if (i % stride) continue;
      const dv = new DataView(L.HEAPU8.buffer, pp, 12);
      pos[3 * k] = dv.getInt32(0, true) * h.scale[0] + h.offset[0]; pos[3 * k + 1] = dv.getInt32(4, true) * h.scale[1] + h.offset[1]; pos[3 * k + 2] = dv.getInt32(8, true) * h.scale[2] + h.offset[2]; k++;
    }
    if (stride > 1) ctx.warn(`Point cloud subsampled: every ${stride}th of ${n} points was kept (${m} points; raise opts.maxPoints to keep more).`);
    return { positions: pos, kind: 'pointcloud', meta: { ...meta, decoder: 'LASzip (laz-perf 0.0.7)', pointsRead: m, stride, unitNote: 'LAS/LAZ units are defined by the coordinate reference system records, which are not interpreted; confirm the units.' } };
  } catch (e) {
    const internal = typeof e === 'number' || /exception|abort|unreachable|out of bounds/i.test(String(e?.message));
    if (internal) { dropKernel('LASzip'); lz = null; pp = 0; L = null; }                 // the instance may be left inconsistent: start a fresh one next time
    throw Object.assign(new Error(internal ? 'the LASzip decoder rejected the compressed point data (corrupt, truncated or not a LASzip-compressed LAS file)' : e?.message || String(e)), { meta });
  } finally { if (L) { try { if (lz) lz.delete(); } catch { /* decoder already gone */ } if (pp) L._free(pp); L._free(fp); } }
}
