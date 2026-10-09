// Readers backed by the vendored WebAssembly kernels: exact CAD geometry (STEP, IGES, OpenCASCADE BREP) is
// tessellated by OpenCASCADE, and LAZ point clouds are decoded by LASzip. The text parsers stay in charge of
// what they read reliably (header, schema, units, census); the kernel adds the faces. When a kernel cannot be
// loaded the readers fall back to the text-only result and say so.

import { Mesh, guard, str, LIMITS, M4 } from './parsers-util.js';
import { readSTEP, readIGES } from './parsers-cad.js';
import { lasHeader } from './parsers-surface.js';
import { loadOcct, loadLazPerf, dropKernel, nodeImport } from './parsers-wasm.js';
import { parseXT, xtToStep } from './parsers-xt.js';
import { parseACIS, acisToStep } from './parsers-sat.js';
import { parseJT, jtPlacements, jtShapes, jtMul } from './parsers-jt.js';
import { readDXFNative } from './parsers-dxf.js';
import { decodeShapeLOD, jtCodecUse, jtOptions } from './parsers-jtshape.js';
import { StepWriter } from './parsers-brep.js';
import { unitFromWord } from './parsers-util.js';

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

const UTF8 = new TextEncoder();
const listCounts = (o) => Object.entries(o).map(([k, n]) => `${n} × ${k}`).join(', ');
/**
 * Common tail of the kernel-format translators: tessellate the generated STEP text and attach the translator's
 * statistics. `tr` = { step, stats, unit, meta, label, solid }.
 */
async function finishTranslation(tr, ctx) {
  const { stats, label } = tr, skipped = Object.values(stats.skippedFaces).reduce((s, n) => s + n, 0), approx = Object.values(stats.approximatedCurves).reduce((s, n) => s + n, 0);
  const meta = { ...tr.meta, translation: stats };
  if (skipped) ctx.warn(`${label}: ${skipped} of ${stats.faces} face(s) could not be translated and are missing from the result (${listCounts(stats.skippedFaces)}).`);
  if (approx) ctx.warn(`${label}: ${approx} edge curve(s) are approximate — ${tr.curveNote} (${listCounts(stats.approximatedCurves)}); the faces they bound are trimmed to that accuracy.`);
  if (!tr.step) throw Object.assign(new Error(stats.faces ? 'none of the faces could be translated' : 'the file holds no faces (wire, point or empty bodies only)'), { meta });
  if (ctx.opts.wasm === false) throw Object.assign(new Error('the entity graph was read, but tessellation was switched off (opts.wasm = false)'), { meta });
  let res;
  try { res = await tessellate('step', UTF8.encode(tr.step), ctx, tr.unit, tr.solid ? 'solid' : 'surface'); }
  catch (e) { throw Object.assign(new Error(`the entity graph was read (${stats.facesTranslated} faces translated to STEP), but the kernel could not tessellate it — ${e.message}`), { meta }); }
  const { out, faces, params } = res;
  return { ...out, units: tr.units, meta: { ...meta, kernel: 'OpenCASCADE (occt-import-js 0.0.23) via in-memory STEP', tessellation: { ...params, bodies: out.groups.length, faces: faces.length, triangles: out.elements[0].count }, brepFaces: faces.slice(0, 20000) }, support: skipped ? 'partial' : undefined };
}

/** Parasolid XT (text or binary): schema-driven node stream → B-rep graph → STEP → kernel tessellation. */
export async function readXT(b, ctx) {
  const parsed = parseXT(b), tr = xtToStep(parsed), solid = [...parsed.nodes.values()].some((n) => n.kind === 'BODY' && n.f.body_type === 1);
  const meta = { header: parsed.meta.header, encoding: parsed.meta.encoding, modellerVersion: parsed.meta.modellerVersion, schema: parsed.meta.schema, embeddedSchema: parsed.meta.embeddedSchema, nodeCount: parsed.meta.nodeCount, bodyNames: tr.stats.bodyNames, unitNote: 'Parasolid models are defined in metres (kernel convention: size box 1000, linear resolution 1e-8); the XT file carries no unit field.' };
  const inf = parsed.meta.inferredSchema;
  if (inf) {
    meta.inferredSchema = inf;
    ctx.warn(`Parasolid XT: the file uses schema ${inf.version} without embedding it. Its node layouts were inferred from the documented base schema and from schemas learned from other versions, and accepted because the whole file then parses to its terminator with every topological pointer on a node of the right type and no other candidate layout does${inf.guessedLayouts.length ? `; ${inf.guessedLayouts.join(', ')} had to be guessed` : ''}${inf.undecidedFields.length ? `; ${inf.undecidedFields.join(', ')} could not be decided and are ignored` : ''}. Check a known dimension.`);
  }
  return finishTranslation({ ...tr, meta, label: 'Parasolid XT', curveNote: 'intersection and surface-parameter curves were replaced by polylines through their defining points', unit: 'm', units: { length: 'm', source: 'file' }, solid }, ctx);
}

const PER_METRE = { m: 1, mm: 1000, cm: 100, in: 1 / 0.0254, ft: 1 / 0.3048 };
/**
 * JT: container and scene graph natively. Per part, exact geometry from an embedded XT B-Rep segment (through the
 * Parasolid translator and the kernel) is preferred; parts without one - or all parts when opts.jtTessellation is
 * true, or when the kernel is unavailable - come from the file's own tessellation (Shape LOD segments).
 * Everything is placed by the scene-graph transforms in the file's model unit.
 */
export async function readJT(b, ctx) {
  const jt = await parseJT(b), meta = { ...jt.meta }, length = unitFromWord(jt.meta.units), k = PER_METRE[length] ?? 1;
  const force = ctx.opts.jtTessellation === true, segById = new Map(jt.segments.map((s) => [s.id, s]));
  if (!jt.lsg) ctx.warn(`JT ${meta.version}: the scene graph of this file could not be decoded, so part names, units and assembly placements are unavailable; shapes are shown in their own coordinate systems.`);
  for (const key of Object.keys(jtCodecUse)) delete jtCodecUse[key];
  jtOptions.allowUnverified = ctx.opts.jtUnverified === true;

  // 1. exact geometry from XT B-Rep segments
  let exact = null, xtParts = 0, xtFailed = '';
  if (jt.xtStreams.length && !force && ctx.opts.wasm !== false) {
    const W = new StepWriter(), stats = { bodies: 0, faces: 0, facesTranslated: 0, edges: 0, skippedFaces: {}, approximatedCurves: {}, instances: 0, wireOnlyBodies: 0, bodyNames: [] };
    const placed = jt.lsg ? jtPlacements(jt.lsg, jt.xtStreams.map((x) => x.segment)) : [];
    for (const x of jt.xtStreams) {
      // XT data is in metres; the scene graph places parts in the file's model unit
      let pls = placed.filter((q) => q.segment === x.segment).map((q) => { const M = q.matrix, r = [[M[0], M[4], M[8]], [M[1], M[5], M[9]], [M[2], M[6], M[10]]], det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1]) - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0]) + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]); return { xf: { r, t: [M[12] / k, M[13] / k, M[14] / k], s: k, mirror: det < 0 }, name: q.name }; });
      if (!pls.length) pls = [{ xf: { r: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0], s: k, mirror: false }, name: null }];
      try { xtToStep(parseXT(x.bytes), { writer: W, stats, placements: pls }); xtParts++; } catch (e) { ctx.warn(`JT: an XT B-Rep segment could not be decoded — ${e.message}.`); }
    }
    if (stats.facesTranslated) {
      try { exact = await finishTranslation({ step: W.finish({ unit: '$', tol: 1e-6 * k, source: 'JT / XT B-Rep' }), stats, meta: {}, label: 'JT (XT B-Rep)', curveNote: 'intersection and surface-parameter curves were replaced by polylines through their defining points', unit: 'm', units: null, solid: true }, ctx); }
      catch (e) { xtFailed = e.message; }
    }
  }

  // 2. tessellation for everything that has no exact geometry
  const cache = new Map(), failures = {}, inst = []; let decoded = 0, skippedUnderXT = 0, nVert = 0, nTri = 0;
  const lodWanted = ctx.opts.lod === 'coarsest' ? 99 : Math.max(0, +ctx.opts.lod | 0);
  const top = jt.lsg ? jtShapes(jt.lsg, lodWanted) : jt.segments.filter((s) => s.type === 6 || s.type === 7).map((s) => ({ segment: s.id, name: null, matrix: null, underXT: false }));
  // assemblies may keep each part in a file of its own: those come from opts.companions
  const shapes = [], extFiles = new Map(), missing = new Set(), comp = ctx.opts.companions || [], baseName = (n) => String(n).split(/[\\/]/).pop().toLowerCase();
  for (const sh of top) {
    if (!sh.external) { shapes.push({ ...sh, src: jt }); continue; }
    const key = baseName(sh.external); let child = extFiles.get(key);
    if (child === undefined) {
      const hit = comp.find((c) => baseName(c.name) === key); child = null;
      if (!hit) missing.add(sh.external);
      else { try { const cj = await parseJT(hit.bytes); child = { jt: cj, segs: new Map(cj.segments.map((s) => [s.id, s])), shapes: cj.lsg ? jtShapes(cj.lsg, lodWanted) : cj.segments.filter((s) => s.type === 6 || s.type === 7).map((s) => ({ segment: s.id, name: null, matrix: null })) }; if (cj.meta.units && jt.meta.units && cj.meta.units !== jt.meta.units) ctx.warn(`JT: the part file "${sh.external}" is in ${cj.meta.units} while the assembly is in ${jt.meta.units}; it is placed without unit conversion.`); } catch (e) { failures[`part file ${sh.external}: ${e.message}`] = 1; } }
      extFiles.set(key, child);
    }
    if (child) for (const c of child.shapes) if (!c.external) shapes.push({ segment: c.segment, name: sh.name || c.name, matrix: c.matrix && sh.matrix ? jtMul(c.matrix, sh.matrix) : sh.matrix || c.matrix, underXT: false, src: child.jt, segs: child.segs, key });
  }
  meta.externalParts = { referenced: extFiles.size, supplied: [...extFiles.values()].filter(Boolean).length };
  if (missing.size) ctx.warn(`JT: this assembly keeps its parts in ${extFiles.size} separate file(s), of which ${missing.size} were not supplied and are missing from the result (${[...missing].slice(0, 6).join(', ')}${missing.size > 6 ? ', …' : ''}); pass them in opts.companions.`);
  for (const sh of shapes) {
    if (sh.underXT && exact) { skippedUnderXT++; continue; }
    const src = sh.src, seg = (sh.segs || segById).get(sh.segment); if (!seg) continue;
    const ck = (sh.key || '') + '/' + sh.segment; let lod = cache.get(ck);
    if (lod === undefined) {
      try { const d = await src.payload(seg); lod = await decodeShapeLOD(d, src.le, src.major, src.guidAt(d, 4)); }
      catch (e) { lod = null; const why = e.message || String(e); failures[why] = (failures[why] || 0) + 1; }
      cache.set(ck, lod);
    }
    if (!lod || !lod.triangles.length) continue;
    inst.push({ lod, m: sh.matrix, name: sh.name }); nVert += lod.positions.length / 3; nTri += lod.triangles.length / 3; decoded++;
    guard(nVert, 'JT vertex', LIMITS.verts); guard(nTri, 'JT triangle', LIMITS.elems);
  }
  // instances are written straight into typed arrays (one pass, no per-triangle allocation)
  let tess = null;
  if (inst.length) {
    const positions = new Float64Array(3 * nVert), conn = new Uint32Array(3 * nTri), group = new Int32Array(nTri), groups = []; let vo = 0, to = 0;
    for (const { lod, m, name } of inst) {
      const P = lod.positions, T = lod.triangles, nv = P.length / 3, nt = T.length / 3, g = groups.length; groups.push({ id: g, name: name || `shape ${g + 1}`, kind: 'component', count: nt, tag: g });
      if (m) { const a = m[0], b = m[1], c = m[2], d = m[4], e = m[5], f = m[6], h = m[8], k2 = m[9], l = m[10], tx = m[12], ty = m[13], tz = m[14]; for (let i = 0, o = 3 * vo; i < 3 * nv; i += 3, o += 3) { const x = P[i], y = P[i + 1], z = P[i + 2]; positions[o] = x * a + y * d + z * h + tx; positions[o + 1] = x * b + y * e + z * k2 + ty; positions[o + 2] = x * c + y * f + z * l + tz; } }
      else positions.set(P, 3 * vo);
      const flip = m ? m[0] * (m[5] * m[10] - m[6] * m[9]) - m[1] * (m[4] * m[10] - m[6] * m[8]) + m[2] * (m[4] * m[9] - m[5] * m[8]) < 0 : false;
      for (let i = 0, o = 3 * to; i < 3 * nt; i += 3, o += 3) { conn[o] = T[i] + vo; conn[o + 1] = T[flip ? i + 2 : i + 1] + vo; conn[o + 2] = T[flip ? i + 1 : i + 2] + vo; }
      group.fill(g, to, to + nt); vo += nv; to += nt;
    }
    tess = { positions, elements: [{ type: 'tri3', nodesPer: 3, count: nTri, conn, group }], groups };
  }
  const nFail = Object.values(failures).reduce((s, n) => s + n, 0);
  meta.shapeDecoding = { shapesUsed: decoded, segmentsDecoded: [...cache.values()].filter(Boolean).length, segmentsFailed: nFail, partsFromXT: exact ? xtParts : 0, codecs: { ...jtCodecUse } };
  if (nFail) ctx.warn(`JT: ${nFail} Shape LOD segment(s) could not be decoded and are missing (${listCounts(failures)}).`);
  const soft = ['CDP1 arithmetic', 'JT8 quantised vertex data'].filter((c) => jtCodecUse[c]);
  if (soft.length && decoded) ctx.warn(`JT: this file uses ${soft.join(' and ')}, which were checked on real files only for self-consistency (coherent, consistently oriented meshes) — no exact reference was available; quantised coordinates are accurate to one quantisation step.`);
  if (jtOptions.allowUnverified && decoded) ctx.warn('JT: decoding paths that could not be verified against real data were enabled (opts.jtUnverified); check the result against a known dimension.');
  if (xtFailed) ctx.warn(`JT: the embedded XT B-Rep could not be tessellated (${xtFailed}); the file's own tessellation is used instead.`);
  if (exact && decoded) ctx.warn(`JT: ${xtParts} part(s) come from their exact XT B-Rep (tessellated here) and ${decoded} shape(s) without B-Rep from the file's tessellation.`);
  else if (exact) ctx.warn(`JT: geometry is the exact XT B-Rep of ${xtParts} part(s), tessellated here; pass opts.jtTessellation = true to use the file's own tessellation instead.`);
  if (jt.meta.jtBrepSegments) ctx.warn(`JT: ${jt.meta.jtBrepSegments} JT B-Rep segment(s) are not decoded.`);
  const units = { length: length || null, source: length ? 'file' : null };
  if (!jt.lsg && exact && !tess) { units.length = 'm'; units.source = 'file'; }      // undecoded scene graph: XT geometry stays in metres
  if (jt.lsg && jt.meta.units && !length) ctx.warn(`JT: the model unit "${jt.meta.units}" has no platform equivalent; assign the units manually.`);
  if (!exact && !tess) throw Object.assign(new Error(nFail ? 'no shape could be decoded' : missing.size ? `the assembly holds no geometry of its own: its ${missing.size} part file(s) must be supplied in opts.companions` : `the file holds no XT B-Rep and no tri-strip tessellation that could be read (${meta.shapeSegments} shape segment(s))`), { meta });
  if (exact && !tess) return { ...exact, units, meta: { ...meta, ...exact.meta } };
  if (!exact) return { ...tess, kind: 'surface', units, meta, support: nFail ? 'partial' : undefined };
  // merge: exact bodies first, then tessellated shapes
  const n0 = exact.positions.length / 3, positions = new Float64Array(exact.positions.length + tess.positions.length); positions.set(exact.positions); positions.set(tess.positions, exact.positions.length);
  const e0 = exact.elements[0], e1 = tess.elements[0], conn = new Uint32Array(e0.conn.length + e1.conn.length), group = new Int32Array(e0.count + e1.count); conn.set(e0.conn); for (let i = 0; i < e1.conn.length; i++) conn[e0.conn.length + i] = e1.conn[i] + n0;
  group.set(e0.group || new Int32Array(e0.count).fill(-1)); for (let i = 0; i < e1.count; i++) group[e0.count + i] = e1.group[i] + exact.groups.length;
  return { positions, elements: [{ type: 'tri3', nodesPer: 3, count: e0.count + e1.count, conn, group }], groups: [...exact.groups, ...tess.groups], kind: 'surface', units, meta: { ...meta, ...exact.meta }, support: exact.support };
}

/** ACIS SAT / SAB: entity records → B-rep graph → STEP → kernel tessellation. */
export async function readACIS(b, ctx) {
  const parsed = parseACIS(b), h = parsed.header, tr = acisToStep(parsed);
  const length = { 1: 'mm', 10: 'cm', 1000: 'm', 25.4: 'in', 304.8: 'ft' }[h.unitsMM] ?? null;
  const census = {}; for (const r of parsed.records.values()) census[r.type] = (census[r.type] || 0) + 1;
  const meta = { encoding: parsed.encoding, saveVersion: h.version, product: h.product, acisVersion: h.acisVersion, date: h.date, millimetresPerUnit: h.unitsMM, resabs: h.resabs, records: parsed.records.size, census, bodyNames: tr.stats.bodyNames };
  const nApprox = Object.values(tr.stats.approximatedSurfaces).reduce((s, n) => s + n, 0);
  if (nApprox) ctx.warn(`ACIS: ${nApprox} procedural spline surface(s) (${listCounts(tr.stats.approximatedSurfaces)}) are represented by the approximating B-spline stored in the file, accurate to the file's fit tolerance rather than exact.`);
  if (h.unitsMM && !length) ctx.warn(`ACIS: the header gives ${h.unitsMM} mm per model unit, which has no platform unit; assign the units manually.`);
  return finishTranslation({ ...tr, meta, label: 'ACIS', curveNote: 'procedural curves are represented by the approximating B-spline stored in the file', unit: 'mm', units: { length, source: length ? 'file' : null }, solid: true }, ctx);
}

/**
 * DXF: entities natively (faces, wires, blocks, layers, units); the ACIS data of 3DSOLID / BODY / REGION / SURFACE
 * entities (obfuscated SAT text up to R2010, binary SAB in ACDSDATA from R2013) goes through the ACIS translator and
 * the kernel, one solid at a time, and is placed by the INSERT transforms it sits under.
 */
export async function readDXF(b, ctx) {
  const d = readDXFNative(b, ctx.opts), meta = { ...d.meta, acis: { solids: d.acis.length, translated: 0, failed: 0, sat: 0, sab: 0, saveVersions: {} } };
  for (const w of [...new Set(d.warnings)]) ctx.warn(`DXF: ${w}`);
  const extra = [], failures = {}, notes = {};
  if (d.acis.length && ctx.opts.wasm === false) ctx.warn(`DXF: ${d.acis.length} ACIS solid(s) / surface(s) were found but not tessellated (opts.wasm = false).`);
  else for (const a of d.acis) {
    if (a.sab) meta.acis.sab++; else meta.acis.sat++;
    try {
      const parsed = parseACIS(a.sab || UTF8.encode(a.sat)), tr = acisToStep(parsed); meta.acis.saveVersions[parsed.header.version] = (meta.acis.saveVersions[parsed.header.version] || 0) + 1;
      const res = await finishTranslation({ ...tr, meta: {}, label: 'ACIS', curveNote: 'procedural curves are represented by the approximating B-spline stored in the file', unit: 'mm', units: null, solid: a.type === '3DSOLID' }, { opts: ctx.opts, warn: (w) => { notes[w] = (notes[w] || 0) + 1; } });
      const P = res.positions.slice(), m = a.matrix, e = res.elements.find((x) => x.type === 'tri3'); if (!e) throw new Error('no triangles');
      if (m) for (let i = 0; i < P.length; i += 3) { const q = M4.apply(m, P[i], P[i + 1], P[i + 2]); P[i] = q[0]; P[i + 1] = q[1]; P[i + 2] = q[2]; }
      const conn = e.conn.slice(); if (m && M4.det3(m) < 0) for (let i = 0; i < conn.length; i += 3) { const t = conn[i + 1]; conn[i + 1] = conn[i + 2]; conn[i + 2] = t; }
      extra.push({ P, conn, name: `${a.type} ${a.handle || extra.length + 1}`, layer: a.layer }); meta.acis.translated++;
    } catch (e) { meta.acis.failed++; const why = String(e.message || e).replace(/^.* — /, ''); failures[why] = (failures[why] || 0) + 1; }
  }
  if (meta.acis.failed) ctx.warn(`DXF: ${meta.acis.failed} of ${d.acis.length} ACIS solid(s) / surface(s) could not be translated and are missing (${listCounts(failures)}).`);
  const nn = Object.entries(notes); if (nn.length) ctx.warn(`DXF ACIS solids: ${nn.slice(0, 4).map(([w, n]) => `${n > 1 ? `[${n} solids] ` : ''}${w}`).join(' ')}`);
  if (!extra.length) { if (!d.part.positions.length) throw Object.assign(new Error(d.acis.length ? 'the drawing holds only ACIS solids, none of which could be translated' : 'the drawing holds no faces, wires or solids that could be read'), { meta }); return { ...d.part, units: d.units, meta, support: meta.acis.failed ? 'partial' : undefined }; }
  // merge the tessellated solids behind the native entities
  let nv = d.part.positions.length / 3, nt = 0; const n0 = nv; for (const x of extra) { nv += x.P.length / 3; nt += x.conn.length / 3; }
  const positions = new Float64Array(3 * nv); positions.set(d.part.positions);
  const elements = d.part.elements.map((e) => ({ ...e, group: e.group || new Int32Array(e.count).fill(-1) })), groups = d.part.groups.slice();
  let tri = elements.find((e) => e.type === 'tri3'); const old = tri ? tri.count : 0, conn = new Uint32Array(3 * (old + nt)), group = new Int32Array(old + nt);
  if (tri) { conn.set(tri.conn); group.set(tri.group); } else { tri = { type: 'tri3', nodesPer: 3 }; elements.push(tri); }
  let vo = n0, to = old;
  for (const x of extra) { positions.set(x.P, 3 * vo); for (let i = 0; i < x.conn.length; i++) conn[3 * to + i] = x.conn[i] + vo; const g = groups.length; groups.push({ id: g, name: x.name, kind: 'body', count: x.conn.length / 3, tag: x.layer }); group.fill(g, to, to + x.conn.length / 3); vo += x.P.length / 3; to += x.conn.length / 3; }
  tri.conn = conn; tri.group = group; tri.count = old + nt;
  return { positions, elements, groups, kind: 'surface', units: d.units, meta, support: meta.acis.failed ? 'partial' : undefined };
}

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
