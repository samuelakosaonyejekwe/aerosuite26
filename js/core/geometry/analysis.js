// Inspection, healing, measurement, slicing and export of imported models.
// All functions are pure: they return new objects and never modify the model they are given.
// Topology (watertightness, manifoldness, orientation) is evaluated on vertices identified by EXACT coincidence,
// so formats that store every facet separately (STL) are judged fairly without moving any vertex; actual
// welding, re-orientation and removal of facets only ever happens in heal(), which logs what it changed.

import * as N from '../numerics.js';
import { ET, FACES, EDGES, shapeOf, displayTriangles, displayLines } from './parsers-util.js';

// ---------- shared helpers ----------
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const len = (a) => Math.hypot(a[0], a[1], a[2]);
const pt = (P, i) => [P[3 * i], P[3 * i + 1], P[3 * i + 2]];
const clone = (v) => (typeof structuredClone === 'function' ? structuredClone(v) : JSON.parse(JSON.stringify(v)));

function bboxOf(P) {
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < P.length; i += 3) for (let k = 0; k < 3; k++) { const v = P[i + k]; if (v < min[k]) min[k] = v; if (v > max[k]) max[k] = v; }
  if (!P.length || !min.every(Number.isFinite) || !max.every(Number.isFinite)) return { min: [0, 0, 0], max: [0, 0, 0], size: [0, 0, 0] };
  return { min, max, size: [max[0] - min[0], max[1] - min[1], max[2] - min[2]] };
}

/**
 * Map every vertex to a representative within `tol` (tol = 0: exactly coincident vertices only).
 * Sweep along a skewed direction so axis-aligned grids do not pile up on one key.
 * Returns { rep, merged, maxMove }.
 */
function weldMap(P, tol) {
  const n = P.length / 3, rep = new Uint32Array(n), key = new Float64Array(n), order = new Uint32Array(n);
  const d = [0.5426, 0.709, 0.4505]; let scale = 0;
  for (let i = 0; i < n; i++) { key[i] = d[0] * P[3 * i] + d[1] * P[3 * i + 1] + d[2] * P[3 * i + 2]; order[i] = i; const a = Math.abs(key[i]); if (a > scale) scale = a; }
  order.sort((a, b) => key[a] - key[b] || a - b);
  const win = tol + 1e-12 * scale, tol2 = tol * tol; let merged = 0, maxMove = 0;
  for (let a = 0; a < n; a++) {
    const i = order[a], x = P[3 * i], y = P[3 * i + 1], z = P[3 * i + 2]; rep[i] = i;
    for (let c = a - 1; c >= 0; c--) {
      const j = order[c]; if (key[i] - key[j] > win) break;
      if (rep[j] !== j) continue;
      const dx = x - P[3 * j], dy = y - P[3 * j + 1], dz = z - P[3 * j + 2], d2 = dx * dx + dy * dy + dz * dz;
      if (d2 <= tol2) { rep[i] = j; merged++; if (d2 > maxMove) maxMove = d2; break; }
    }
  }
  return { rep, merged, maxMove: Math.sqrt(maxMove) };
}

/** Edge census of a triangle list on canonical vertex ids. */
function topology(P, T) {
  const nV = P.length / 3, nT = T.length / 3, { rep, merged } = weldMap(P, 0);
  const edges = new Map(), parent = new Int32Array(nV).fill(-1);
  const find = (i) => { let r = i; while (parent[r] >= 0) r = parent[r]; while (parent[i] >= 0) { const nx = parent[i]; parent[i] = r; i = nx; } return r; };
  const used = new Uint8Array(nV);
  let degenerate = 0, area = 0;
  for (let t = 0; t < nT; t++) {
    const a = rep[T[3 * t]], b = rep[T[3 * t + 1]], c = rep[T[3 * t + 2]];
    if (a === b || b === c || a === c) { degenerate++; continue; }
    const u = sub(pt(P, b), pt(P, a)), v = sub(pt(P, c), pt(P, a)), w = sub(pt(P, c), pt(P, b)), cr = cross(u, v), A2 = len(cr), l2 = Math.max(dot(u, u), dot(v, v), dot(w, w));
    area += 0.5 * A2;
    if (A2 <= 1e-12 * l2) degenerate++;                 // zero-area sliver with distinct vertices: counted, but still part of the topology
    const tri = [a, b, c];
    for (let e = 0; e < 3; e++) {
      const p = tri[e], q = tri[(e + 1) % 3], k = p < q ? p * nV + q : q * nV + p;
      edges.set(k, (edges.get(k) || 0) + (p < q ? 1 : 65536));
    }
    used[a] = used[b] = used[c] = 1;
    const ra = find(a), rb = find(b); if (ra !== rb) parent[rb] = ra;
    const rc = find(c), r2 = find(a); if (rc !== r2) parent[rc] = r2;
  }
  let boundary = 0, nonManifold = 0, conflicts = 0;
  for (const v of edges.values()) { const f = v & 65535, bk = Math.floor(v / 65536), tot = f + bk; if (tot === 1) boundary++; else if (tot > 2) nonManifold++; else if (f !== bk) conflicts++; }
  let components = 0; for (let i = 0; i < nV; i++) if (used[i] && parent[i] < 0) components++;
  return { rep, duplicateVerts: merged, degenerate, boundary, nonManifold, conflicts, components, area, closed: nT - degenerate > 0 && boundary === 0 && nonManifold === 0 && conflicts === 0 };
}

/** Volume integrals of a triangulated surface about `o` (divergence theorem): signed volume, first and second moments. */
function volumeIntegrals(P, T, o) {
  let V = 0; const m1 = [0, 0, 0], C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let t = 0; t < T.length; t += 3) {
    const a = sub(pt(P, T[t]), o), b = sub(pt(P, T[t + 1]), o), c = sub(pt(P, T[t + 2]), o), d = dot(a, cross(b, c));
    V += d / 6;
    for (let i = 0; i < 3; i++) {
      const si = a[i] + b[i] + c[i]; m1[i] += (d / 24) * si;
      for (let j = 0; j < 3; j++) C[i][j] += (d / 120) * (a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + si * (a[j] + b[j] + c[j]));
    }
  }
  return { V, m1, C };
}

/** Deterministic, area-weighted sample of points on the surface (or of the raw points when there are no triangles). */
function samplePoints(P, T, m = 4000) {
  const out = []; let s = 1234567;
  const rnd = () => ((s = (s * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff);
  const nT = T.length / 3;
  if (!nT) { const n = P.length / 3, st = Math.max(1, Math.floor(n / m)); for (let i = 0; i < n; i += st) out.push(pt(P, i)); return out; }
  const cum = new Float64Array(nT); let tot = 0;
  for (let t = 0; t < nT; t++) { const a = pt(P, T[3 * t]); tot += 0.5 * len(cross(sub(pt(P, T[3 * t + 1]), a), sub(pt(P, T[3 * t + 2]), a))); cum[t] = tot; }
  if (!(tot > 0)) return out;
  for (let k = 0, t = 0; k < m; k++) {
    const target = ((k + 0.5) / m) * tot; while (t < nT - 1 && cum[t] < target) t++;
    const a = pt(P, T[3 * t]), b = pt(P, T[3 * t + 1]), c = pt(P, T[3 * t + 2]); let u = rnd(), v = rnd(); if (u + v > 1) { u = 1 - u; v = 1 - v; }
    out.push([a[0] + u * (b[0] - a[0]) + v * (c[0] - a[0]), a[1] + u * (b[1] - a[1]) + v * (c[1] - a[1]), a[2] + u * (b[2] - a[2]) + v * (c[2] - a[2])]);
  }
  return out;
}

const AX = ['x', 'y', 'z'];
/** Heuristic classification from element content, names, aspect ratios, planform fill and mirror symmetry. */
function classify(model, bb, topo) {
  const C = (label, confidence, why) => ({ label, confidence: Math.round(confidence * 100) / 100, why });
  const kind = model.kind, has = (d) => model.elements.some((e) => ET[e.type]?.dim === d);
  if (kind === 'metadata-only') return C('unreadable here (metadata only)', 0.95, 'The file was recognised but its geometry is not decoded; see the conversion pathway.');
  if (kind === 'cad-brep') return C('exact CAD geometry (B-rep, not tessellated)', 0.9, 'Exact CAD entities were found; only their points/polylines are available here.');
  if (kind === 'pointcloud') return C('point cloud / scan', model.lines.length ? 0.6 : 0.9, model.lines.length ? 'Points joined in file order with no faces — a planar coordinate table or section curve.' : 'Vertices without any connectivity.');
  if (kind === 'structural-mesh') return C('structural mesh', 0.85, `Finite-element deck with ${model.elements.map((e) => `${e.count} ${e.type}`).join(', ')}${model.groups.some((g) => g.kind === 'property' || g.kind === 'material') ? ' and property/material assignments' : ''}.`);
  if (has(3)) {
    const cfd = ['su2', 'fluent', 'openfoam', 'plot3d', 'tecplot'].includes(model.format), bnd = model.groups.filter((g) => g.kind === 'boundary').length;
    if (cfd || bnd >= 2 || topo.components >= 2) return C('CFD volume mesh', cfd ? 0.85 : 0.6, cfd ? 'Volume cells in a flow-solver mesh format with boundary markers.' : bnd >= 2 ? `Volume cells with ${bnd} named boundary groups, as a fluid domain would have.` : `The volume boundary has ${topo.components} separate closed parts (e.g. far-field plus an immersed body).`);
    return C('structural mesh', 0.45, 'Volume cells with a single closed boundary and no boundary markers: most likely a solid part, but it could equally be a simple fluid domain.');
  }
  if (model.elements.length && model.elements.every((e) => ET[e.type]?.dim === 1)) return C('structural mesh', 0.6, 'Only one-dimensional (beam/rod) elements.');
  // surface: geometry-based guess
  const names = (model.name + ' ' + model.groups.map((g) => g.name).join(' ')).toLowerCase();
  const hint = [['complete aircraft', /aircraft|airplane|aeroplane|airframe|\bjet\b|uav|drone|glider/], ['rotor or propeller blade', /blade|rotor|prop(eller)?\b/], ['wing / lifting surface', /wing|tailplane|stabili[sz]er|\bfin\b|canard|aileron|flap|aerofoil|airfoil/], ['fuselage / body of revolution', /fuselage|nacelle|body|fairing|nose|pod\b/]].find(([, re]) => re.test(names));
  const order = [0, 1, 2].sort((a, b) => bb.size[b] - bb.size[a]), L = bb.size[order[0]], Mi = bb.size[order[1]], S = bb.size[order[2]];
  if (!(L > 0)) return C('generic part', 0.2, 'Degenerate extents.');
  const pts = samplePoints(model.positions, model.triangles), G = 14, plan = new Uint8Array(G * G), vox = new Set(), q = (v, k, g) => Math.min(g - 1, Math.max(0, Math.floor(((v - bb.min[k]) / (bb.size[k] || 1)) * g)));
  for (const p of pts) { plan[q(p[order[0]], order[0], G) * G + q(p[order[1]], order[1], G)] = 1; vox.add(q(p[0], 0, 8) * 64 + q(p[1], 1, 8) * 8 + q(p[2], 2, 8)); }
  const fill = plan.reduce((s, v) => s + v, 0) / (G * G);
  let sym = 0, symAxis = 0;
  for (let ax = 0; ax < 3; ax++) { let inter = 0; for (const v of vox) { const c = [Math.floor(v / 64), Math.floor(v / 8) % 8, v % 8]; c[ax] = 7 - c[ax]; if (vox.has(c[0] * 64 + c[1] * 8 + c[2])) inter++; } const s = vox.size ? inter / vox.size : 0; if (s > sym) { sym = s; symAxis = ax; } }
  const ar = L / (Mi || 1e-300), thin = S / (Mi || 1e-300), f2 = (v) => v.toFixed(2);
  const facts = `extents ${AX[order[0]]}:${AX[order[1]]}:${AX[order[2]]} = 1 : ${f2(Mi / L)} : ${f2(S / L)}, planform fill ${f2(fill)}, mirror symmetry ${f2(sym)} about ${AX[symAxis]}`;
  let guess;
  if (Mi / L >= 0.45 && S / L <= 0.45 && fill <= 0.6 && sym >= 0.7 && topo.area > 0 && model.triangles.length / 3 >= 200) guess = ['complete aircraft', 0.55, `Two comparable long dimensions, a flat overall envelope, a sparse cross-shaped planform and a plane of symmetry (${facts}).`];
  else if (ar >= 7 && thin <= 0.6) guess = ['rotor or propeller blade', 0.5, `Very slender and thin (${facts}).`];
  else if (ar >= 1.5 && thin <= 0.3 && fill >= 0.45) guess = ['wing / lifting surface', 0.55, `One long dimension with a thin section and a filled planform (${facts}).`];
  else if (ar >= 2.5 && thin >= 0.6) guess = ['fuselage / body of revolution', 0.5, `Slender with a near-equal cross-section in the other two directions (${facts}).`];
  else if (thin <= 0.2 && fill >= 0.6) guess = ['wing / lifting surface', 0.3, `Thin and plate-like with a low aspect ratio (${facts}).`];
  else guess = ['generic part', 0.4, `No aircraft-component signature (${facts}).`];
  if (hint) return hint[0] === guess[0] ? C(guess[0], Math.min(0.9, guess[1] + 0.3), `${guess[2]} The file/group names agree.`) : C(hint[0], 0.5, `Named in the file/group names; geometry alone suggested "${guess[0]}" (${facts}).`);
  return C(...guess);
}

/** Geometric and topological census of a model, with a best-effort classification the user can overrule. */
export function analyse(model) {
  const P = model.positions, T = model.triangles, bb = bboxOf(P), nV = P.length / 3, nT = T.length / 3;
  const tp = topology(P, T), centre = [0, 1, 2].map((k) => 0.5 * (bb.min[k] + bb.max[k]));
  let volume = null, centroid = centre.slice(), inwardNormals = false;
  if (tp.closed) {
    const { V, m1 } = volumeIntegrals(P, T, centre);
    if (V !== 0) { volume = Math.abs(V); inwardNormals = V < 0; centroid = [0, 1, 2].map((k) => centre[k] + m1[k] / V); }
  } else if (nT && tp.area > 0) {
    const c = [0, 0, 0];
    for (let t = 0; t < nT; t++) { const a = pt(P, T[3 * t]), b = pt(P, T[3 * t + 1]), d = pt(P, T[3 * t + 2]), A = 0.5 * len(cross(sub(b, a), sub(d, a))); for (let k = 0; k < 3; k++) c[k] += (A * (a[k] + b[k] + d[k])) / 3; }
    centroid = c.map((v) => v / tp.area);
  } else if (nV) { const c = [0, 0, 0]; for (let i = 0; i < P.length; i += 3) for (let k = 0; k < 3; k++) c[k] += P[i + k]; centroid = c.map((v) => v / nV); }
  return {
    bbox: bb, nVerts: nV, nTris: nT, area: tp.area, volume, centroid, watertight: tp.closed || (nT - tp.degenerate > 0 && tp.boundary === 0 && tp.nonManifold === 0),
    boundaryEdges: tp.boundary, nonManifoldEdges: tp.nonManifold, duplicateVerts: tp.duplicateVerts, degenerateTris: tp.degenerate, orientationConflicts: tp.conflicts, components: tp.components,
    inwardNormals, classification: classify(model, bb, tp), dims: { length_x: bb.size[0], span_y: bb.size[1], height_z: bb.size[2] },
  };
}

// ---------- mesh quality ----------
const angleAt = (p, a, b) => { const u = sub(a, p), v = sub(b, p), l = len(u) * len(v); return l > 0 ? (Math.acos(Math.max(-1, Math.min(1, dot(u, v) / l))) * 180) / Math.PI : 0; };
const det3 = (a, b, c) => dot(a, cross(b, c));
function stat(a) {
  let min = Infinity, max = -Infinity, s = 0; for (const v of a) { if (v < min) min = v; if (v > max) max = v; s += v; }
  if (!a.length) return { min: 0, max: 0, mean: 0, hist: { centers: [], counts: [] } };
  const h = N.histogram(a, 12);
  return { min, max, mean: s / a.length, hist: { centers: h.centers, counts: h.counts } };
}
/** Minimum scaled corner Jacobian sign of a volume element: > 0 valid, ≤ 0 inverted or collapsed. */
function minJacobian(shape, p) {
  if (shape === 'tet') return det3(sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
  let m = Infinity; const J = (c, a, b, d) => { const v = det3(sub(p[a], p[c]), sub(p[b], p[c]), sub(p[d], p[c])); if (v < m) m = v; };
  if (shape === 'hex') for (let c = 0; c < 4; c++) { J(c, (c + 1) % 4, (c + 3) % 4, c + 4); J(c + 4, 4 + ((c + 3) % 4), 4 + ((c + 1) % 4), c); }
  else if (shape === 'wedge') for (let c = 0; c < 3; c++) { J(c, (c + 1) % 3, (c + 2) % 3, c + 3); J(c + 3, 3 + ((c + 2) % 3), 3 + ((c + 1) % 3), c); }
  else if (shape === 'pyr') for (let c = 0; c < 4; c++) J(c, (c + 1) % 4, (c + 3) % 4, 4);
  return m;
}
/**
 * Element quality. Aspect ratio is normalised to 1 for the ideal element (triangle: l_max·perimeter/(4√3·A);
 * tetrahedron: l_max/(2√6·r_in); other shapes: longest/shortest edge). Skewness is the equiangle skew of the
 * face corner angles (ideal 60° on triangular faces, 90° on quadrilateral faces). Quadratic elements are
 * judged on their corner nodes. Works on the display triangles when the model has no explicit elements.
 */
export function meshQuality(model) {
  const P = model.positions, issues = [], perType = [];
  let blocks = model.elements.filter((e) => ET[e.type] && ET[e.type].dim >= 2);
  if (!blocks.length && model.triangles.length) blocks = [{ type: 'tri3', nodesPer: 3, count: model.triangles.length / 3, conn: model.triangles }];
  const worst = { aspect: 0, skew: 0, minAngle: 180 }; let negTotal = 0, collapsed = 0;
  for (const el of blocks) {
    const info = ET[el.type], shape = shapeOf(el.type), nc = info.c, n = el.count;
    const faces = info.dim === 3 ? FACES[shape] : [Array.from({ length: nc }, (_, i) => i)], edges = EDGES[shape];
    const aspect = new Float64Array(n), skew = new Float64Array(n), minA = new Float64Array(n), maxA = new Float64Array(n); let neg = 0;
    for (let e = 0; e < n; e++) {
      const p = []; for (let k = 0; k < nc; k++) p.push(pt(P, el.conn[e * el.nodesPer + k]));
      let lmin = Infinity, lmax = 0, per = 0; for (const [a, b] of edges) { const l = len(sub(p[a], p[b])); if (l < lmin) lmin = l; if (l > lmax) lmax = l; per += l; }
      let amin = 180, amax = 0, sk = 0;
      for (const f of faces) { const m = f.length, ideal = m === 3 ? 60 : 90; for (let k = 0; k < m; k++) { const a = angleAt(p[f[k]], p[f[(k + 1) % m]], p[f[(k + m - 1) % m]]); if (a < amin) amin = a; if (a > amax) amax = a; const s = Math.max((a - ideal) / (180 - ideal), (ideal - a) / ideal); if (s > sk) sk = s; } }
      let ar;
      if (shape === 'tri') { const A = 0.5 * len(cross(sub(p[1], p[0]), sub(p[2], p[0]))); ar = A > 0 ? (lmax * per) / (4 * Math.sqrt(3) * A) : Infinity; }
      else if (shape === 'tet') { let Af = 0; for (const f of faces) Af += 0.5 * len(cross(sub(p[f[1]], p[f[0]]), sub(p[f[2]], p[f[0]]))); const V = Math.abs(det3(sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]))) / 6; ar = V > 0 ? (lmax * Af) / (2 * Math.sqrt(6) * 3 * V) : Infinity; }
      else ar = lmin > 0 ? lmax / lmin : Infinity;
      if (!(ar < 1e9)) { ar = 1e9; collapsed++; }
      if (info.dim === 3) { if (!(minJacobian(shape, p) > 0)) neg++; }
      else if (shape === 'quad') { const nn = [0, 1, 2, 3].map((k) => cross(sub(p[(k + 1) % 4], p[k]), sub(p[(k + 3) % 4], p[k]))), mean = nn.reduce((s, v) => [s[0] + v[0], s[1] + v[1], s[2] + v[2]], [0, 0, 0]); if (nn.some((v) => dot(v, mean) <= 0)) neg++; }
      aspect[e] = ar; skew[e] = Math.min(1, sk); minA[e] = amin; maxA[e] = amax;
    }
    const rec = { type: el.type, count: n, aspect: stat(aspect), skew: stat(skew), minAngle: stat(minA), maxAngle: stat(maxA), negJacobian: neg };
    perType.push(rec);
    worst.aspect = Math.max(worst.aspect, rec.aspect.max); worst.skew = Math.max(worst.skew, rec.skew.max); worst.minAngle = Math.min(worst.minAngle, rec.minAngle.min);
    if (neg && neg === n && info.dim === 3 && n > 1) issues.push(`All ${n} ${el.type} elements have negative orientation: this is a node-ordering convention difference (or a mirrored model), not ${n} inverted cells. Mirror/renumber before handing to a solver.`);
    else if (neg) { negTotal += neg; issues.push(`${neg} of ${n} ${el.type} element(s) are inverted or non-convex (non-positive Jacobian).`); }
    if (info.n > info.c) issues.push(`${el.type}: quality is evaluated on the corner nodes only (mid-side node placement is not checked).`);
  }
  if (!perType.length) return { perType, worst: { aspect: 0, skew: 0, minAngle: 0 }, grade: 'unusable', issues: ['The model has no surface or volume elements to evaluate.'] };
  if (collapsed) issues.push(`${collapsed} element(s) are collapsed (zero area/volume or a zero-length edge).`);
  if (worst.skew > 0.95) issues.push(`Worst equiangle skewness ${worst.skew.toFixed(3)} (> 0.95): effectively degenerate cells.`);
  else if (worst.skew > 0.85) issues.push(`Worst equiangle skewness ${worst.skew.toFixed(3)} (> 0.85): poor cells that can stall or destabilise a solver.`);
  if (worst.minAngle < 10) issues.push(`Smallest corner angle ${worst.minAngle.toFixed(2)}° (< 10°).`);
  if (worst.aspect > 100) issues.push(`Largest aspect ratio ${worst.aspect.toExponential(2)}: acceptable only in deliberately stretched boundary-layer or shell regions.`);
  const grade = negTotal || collapsed || worst.skew > 0.98 || worst.minAngle < 1 ? 'unusable' : worst.skew > 0.85 || worst.minAngle < 10 ? 'poor' : worst.skew > 0.6 || worst.minAngle < 25 || worst.aspect > 20 ? 'acceptable' : 'good';
  return { perType, worst, grade, issues };
}

// ---------- model copying ----------
function copyModel(m) {
  return {
    ...m, positions: m.positions.slice(), triangles: m.triangles.slice(), lines: m.lines.slice(),
    elements: m.elements.map((e) => ({ ...e, conn: e.conn.slice(), ...(e.group ? { group: e.group.slice() } : {}) })),
    groups: m.groups.map((g) => ({ ...g })), units: { ...m.units }, meta: clone(m.meta), log: m.log.map((l) => ({ ...l })), warnings: m.warnings.slice(),
  };
}
/** Group id per display triangle when the triangles map one-to-one onto the surface elements (else null). */
function triangleGroups(m) {
  const g = []; let any = false;
  for (const el of m.elements) { const info = ET[el.type]; if (!info || info.dim !== 2) continue; const per = info.c === 3 ? 1 : 2; for (let e = 0; e < el.count; e++) { const v = el.group ? el.group[e] : -1; if (v >= 0) any = true; for (let k = 0; k < per; k++) g.push(v); } }
  return any && g.length === m.triangles.length / 3 ? Int32Array.from(g) : null;
}

/**
 * Explicit, logged repair. Welds vertices closer than weldTol (the surviving vertex keeps its position, so the
 * largest displacement is bounded by weldTol and reported), removes degenerate and duplicate facets, makes the
 * facet orientation consistent (outward for closed parts) and drops unreferenced vertices.
 * Surface models are rebuilt as triangles; models with volume elements are only welded/compacted.
 */
export function heal(model, { weldTol = null, removeDegenerate = true, fixOrientation = true, removeUnreferenced = true } = {}) {
  const m = copyModel(model), log = [], add = (step, detail) => log.push({ step, detail });
  const bb0 = bboxOf(m.positions), diag = len(bb0.size), tol = weldTol === null || weldTol === undefined ? 1e-6 * diag : Math.max(0, +weldTol || 0);
  const hasVol = m.elements.some((e) => ET[e.type]?.dim === 3), hasConn = m.triangles.length || m.lines.length || m.elements.length;
  let changed = false, P = m.positions, T = m.triangles, tg = hasVol ? null : triangleGroups(m);
  const hadSurfElems = !hasVol && m.elements.some((e) => ET[e.type]?.dim === 2), nonTri = !hasVol && m.elements.some((e) => ET[e.type]?.dim === 2 && e.type !== 'tri3');

  // vertices that nothing references before any repair (reported separately from welded duplicates)
  let stray = 0;
  if (hasConn) { const u0 = new Uint8Array(P.length / 3); for (const v of T) u0[v] = 1; for (const v of m.lines) u0[v] = 1; for (const el of m.elements) for (const v of el.conn) u0[v] = 1; for (const v of u0) if (!v) stray++; }

  // 1. weld
  const w = weldMap(P, tol);
  if (w.merged) {
    changed = true;
    for (let i = 0; i < T.length; i++) T[i] = w.rep[T[i]];
    for (let i = 0; i < m.lines.length; i++) m.lines[i] = w.rep[m.lines[i]];
    for (const el of m.elements) for (let i = 0; i < el.conn.length; i++) el.conn[i] = w.rep[el.conn[i]];
  }
  add('weld', `tolerance ${tol.toExponential(3)} (${weldTol === null || weldTol === undefined ? '1e-6 × bounding-box diagonal' : 'absolute, user-set'}): ${w.merged} vertex/vertices merged onto an existing vertex; largest distance any vertex moved = ${w.maxMove.toExponential(3)}`);

  // 2. degenerate and duplicate facets
  if (removeDegenerate && !hasVol && T.length) {
    const nT = T.length / 3, keep = new Uint8Array(nT), seen = new Set(), nV = P.length / 3; let deg = 0, dup = 0;
    for (let t = 0; t < nT; t++) {
      const a = T[3 * t], b = T[3 * t + 1], c = T[3 * t + 2];
      const u = sub(pt(P, b), pt(P, a)), v = sub(pt(P, c), pt(P, a)), x = sub(pt(P, c), pt(P, b)), l2 = Math.max(dot(u, u), dot(v, v), dot(x, x));
      if (a === b || b === c || a === c || len(cross(u, v)) <= 1e-12 * l2) { deg++; continue; }
      const s = [a, b, c].sort((p, q) => p - q), k = nV < 200000 ? (s[0] * nV + s[1]) * nV + s[2] : s.join('_');
      if (seen.has(k)) { dup++; continue; }
      seen.add(k); keep[t] = 1;
    }
    if (deg || dup) {
      changed = true; const nk = nT - deg - dup, T2 = new Uint32Array(3 * nk), g2 = tg ? new Int32Array(nk) : null; let o = 0;
      for (let t = 0; t < nT; t++) if (keep[t]) { T2[3 * o] = T[3 * t]; T2[3 * o + 1] = T[3 * t + 1]; T2[3 * o + 2] = T[3 * t + 2]; if (g2) g2[o] = tg[t]; o++; }
      T = T2; tg = g2;
    }
    add('degenerate', `${deg} degenerate (zero-area or repeated-vertex) triangle(s) and ${dup} duplicate triangle(s) removed`);
  } else if (removeDegenerate && hasVol) add('degenerate', 'skipped: the model has volume elements, which are never deleted by heal()');

  // 3. orientation
  if (fixOrientation && !hasVol && T.length) {
    const nT = T.length / 3, nV = P.length / 3, cnt = new Map(), first = new Map(), nbr = new Int32Array(3 * nT).fill(-1);
    const ek = (p, q) => (p < q ? p * nV + q : q * nV + p);
    for (let t = 0; t < nT; t++) for (let e = 0; e < 3; e++) { const k = ek(T[3 * t + e], T[3 * t + (e + 1) % 3]); cnt.set(k, (cnt.get(k) || 0) + 1); }
    for (let t = 0; t < nT; t++) for (let e = 0; e < 3; e++) {
      const k = ek(T[3 * t + e], T[3 * t + (e + 1) % 3]); if (cnt.get(k) !== 2) continue;
      const f = first.get(k); if (f === undefined) first.set(k, 3 * t + e); else { nbr[3 * t + e] = f; nbr[f] = 3 * t + e; }
    }
    const state = new Int8Array(nT), flipTri = (t) => { const b = T[3 * t + 1]; T[3 * t + 1] = T[3 * t + 2]; T[3 * t + 2] = b; const a = nbr[3 * t], c = nbr[3 * t + 2]; nbr[3 * t] = c; nbr[3 * t + 2] = a; if (a >= 0) nbr[a] = 3 * t + 2; if (c >= 0) nbr[c] = 3 * t; };
    let flipped = 0, turned = 0, unorientable = 0, comps = 0;
    for (let seed = 0; seed < nT; seed++) {
      if (state[seed]) continue;
      comps++; const comp = [seed], flips = new Set(); state[seed] = 1; let open = false;
      for (let h = 0; h < comp.length; h++) {
        const t = comp[h];
        for (let e = 0; e < 3; e++) {
          const o = nbr[3 * t + e]; if (o < 0) { if (cnt.get(ek(T[3 * t + e], T[3 * t + (e + 1) % 3])) === 1) open = true; continue; }
          const u = Math.floor(o / 3), same = T[o] === T[3 * t + e];   // same start vertex on the shared edge ⇒ same direction ⇒ inconsistent
          if (!state[u]) { state[u] = 1; if (same) { flipTri(u); flips.add(u); } comp.push(u); } else if (same) unorientable++;
        }
      }
      let flipAll = false;
      if (!open) { let V = 0; for (const t of comp) V += det3(pt(P, T[3 * t]), pt(P, T[3 * t + 1]), pt(P, T[3 * t + 2])); flipAll = V < 0; if (flipAll) turned++; }
      else flipAll = flips.size > comp.length / 2;           // open sheet: keep the majority orientation
      if (flipAll) for (const t of comp) flipTri(t);
      flipped += flipAll ? comp.length - flips.size : flips.size;
    }
    if (flipped) changed = true;
    add('orientation', `${flipped} triangle(s) reversed for a consistent orientation across ${comps} connected part(s)${turned ? `; ${turned} closed part(s) were inside-out and now face outward` : ''}${unorientable ? `; ${unorientable} edge conflict(s) could not be resolved (non-orientable or non-manifold region)` : ''}`);
  } else if (fixOrientation && hasVol) add('orientation', 'skipped: boundary faces of volume elements follow the element node ordering');

  // rebuild the surface elements from the healed triangles
  if (!hasVol) {
    if (hadSurfElems) {
      m.elements = m.elements.filter((e) => ET[e.type]?.dim !== 2);
      if (T.length) m.elements.push({ type: 'tri3', nodesPer: 3, count: T.length / 3, conn: T.slice(), ...(tg ? { group: tg } : {}) });
      if (tg) { const c = new Map(); for (const g of tg) c.set(g, (c.get(g) || 0) + 1); for (const g of m.groups) if (c.has(g.id) || g.kind !== 'material') g.count = c.get(g.id) || 0; }
      if (nonTri && changed) add('elements', `quadrilateral/quadratic surface elements were replaced by the healed linear triangles${tg ? ' (group membership kept)' : ''}`);
    }
    m.triangles = T;
  }

  // 4. unreferenced vertices
  if (removeUnreferenced && hasConn) {
    const nV = P.length / 3, used = new Uint8Array(nV);
    for (const v of m.triangles) used[v] = 1; for (const v of m.lines) used[v] = 1; for (const el of m.elements) for (const v of el.conn) used[v] = 1;
    const map = new Int32Array(nV).fill(-1); let k = 0; for (let i = 0; i < nV; i++) if (used[i]) map[i] = k++;
    if (k < nV) {
      changed = true; const P2 = new Float64Array(3 * k);
      for (let i = 0; i < nV; i++) if (map[i] >= 0) { P2[3 * map[i]] = P[3 * i]; P2[3 * map[i] + 1] = P[3 * i + 1]; P2[3 * map[i] + 2] = P[3 * i + 2]; }
      const re = (a) => { for (let i = 0; i < a.length; i++) a[i] = map[a[i]]; };
      re(m.triangles); re(m.lines); for (const el of m.elements) re(el.conn);
      P = P2;
    }
    add('unreferenced', `${Math.min(stray, nV - k)} unreferenced vertex/vertices and ${Math.max(0, nV - k - stray)} welded-away or orphaned duplicate(s) removed (${k} remain)`);
  } else if (removeUnreferenced) add('unreferenced', 'skipped: a point cloud has no connectivity, so every vertex is kept');
  m.positions = P;
  if (hasVol) { m.triangles = displayTriangles(m.elements, P.length / 3); m.lines = displayLines(m.elements); }
  const bb1 = bboxOf(m.positions), dMax = Math.max(...[0, 1, 2].map((k) => Math.abs(bb1.size[k] - bb0.size[k])));
  add('dimensions', `bounding-box size before ${bb0.size.map((v) => +v.toPrecision(9)).join(' × ')}, after ${bb1.size.map((v) => +v.toPrecision(9)).join(' × ')} (largest change ${dMax.toExponential(2)})`);
  add('result', changed ? 'model changed by healing' : 'no change was necessary');
  for (const l of log) m.log.push({ step: `heal:${l.step}`, detail: l.detail });
  return { model: m, log, changed };
}

/**
 * Affine transform, applied in the order scale → swapYZ → flipX → rotate (about x, then y, then z) → translate.
 * swapYZ is the proper rotation that takes a Y-up model to Z-up: (x, y, z) → (x, −z, y).
 * flipX mirrors in the x = 0 plane; element/facet node order is reversed so orientation stays valid.
 * Pass `units` (e.g. 'm') to record the length unit that the transformed coordinates are in.
 */
export function transform(model, { scale = 1, swapYZ = false, flipX = false, translate = [0, 0, 0], rotateDeg = [0, 0, 0], units } = {}) {
  const m = copyModel(model), P = m.positions, s = Number.isFinite(+scale) && +scale !== 0 ? +scale : 1, r = rotateDeg.map((d) => ((+d || 0) * Math.PI) / 180);
  const [cx, sx, cy, sy, cz, sz] = [Math.cos(r[0]), Math.sin(r[0]), Math.cos(r[1]), Math.sin(r[1]), Math.cos(r[2]), Math.sin(r[2])];
  for (let i = 0; i < P.length; i += 3) {
    let x = P[i] * s, y = P[i + 1] * s, z = P[i + 2] * s, t;
    if (swapYZ) { t = y; y = -z; z = t; }
    if (flipX) x = -x;
    t = y * cx - z * sx; z = y * sx + z * cx; y = t;
    t = x * cy + z * sy; z = -x * sy + z * cy; x = t;
    t = x * cz - y * sz; y = x * sz + y * cz; x = t;
    P[i] = x + (+translate[0] || 0); P[i + 1] = y + (+translate[1] || 0); P[i + 2] = z + (+translate[2] || 0);
  }
  if ((s < 0) !== !!flipX) {                              // net reflection: restore a positive orientation
    const T = m.triangles; for (let t = 0; t < T.length; t += 3) { const b = T[t + 1]; T[t + 1] = T[t + 2]; T[t + 2] = b; }
    const SWAP = { tri: [[1, 2]], quad: [[1, 3]], tet: [[0, 1]], hex: [[0, 4], [1, 5], [2, 6], [3, 7]], wedge: [[0, 3], [1, 4], [2, 5]], pyr: [[1, 3]] };
    for (const el of m.elements) { const sw = SWAP[shapeOf(el.type)]; if (!sw) continue; if (ET[el.type].n > ET[el.type].c) m.warnings.push(`Mirroring: ${el.type} corner nodes were reordered but mid-side nodes keep their slots; re-export from the source system if element orientation matters.`); for (let e = 0; e < el.count; e++) for (const [a, b] of sw) { const o = e * el.nodesPer, v = el.conn[o + a]; el.conn[o + a] = el.conn[o + b]; el.conn[o + b] = v; } }
  }
  if (units) m.units = { length: units, source: 'user' };
  m.log.push({ step: 'transform', detail: `scale ${s}${swapYZ ? ', Y-up → Z-up (x, y, z) → (x, −z, y)' : ''}${flipX ? ', mirrored in x' : ''}, rotate [${rotateDeg.join(', ')}]° about x/y/z, translate [${translate.join(', ')}]${units ? `; length unit set to ${units}` : ''}` });
  return m;
}

/** Volume, mass, centre of gravity and inertia tensor about the cg, by divergence-theorem integration over the surface. */
export function massProperties(model, density = 1) {
  const P = model.positions, T = model.triangles, bb = bboxOf(P), o = [0, 1, 2].map((k) => 0.5 * (bb.min[k] + bb.max[k])), tp = topology(P, T);
  let { V, m1, C } = volumeIntegrals(P, T, o);
  if (V < 0) { V = -V; m1 = m1.map((v) => -v); C = C.map((r) => r.map((v) => -v)); }
  if (!(V > 0)) return { volume: 0, mass: 0, cg: o, inertia: [[0, 0, 0], [0, 0, 0], [0, 0, 0]], closed: false };
  const c = m1.map((v) => v / V), cov = C.map((r, i) => r.map((v, j) => v - V * c[i] * c[j])), tr = cov[0][0] + cov[1][1] + cov[2][2];
  const inertia = cov.map((r, i) => r.map((v, j) => density * ((i === j ? tr : 0) - v)));
  return { volume: V, mass: density * V, cg: [o[0] + c[0], o[1] + c[1], o[2] + c[2]], inertia, closed: tp.closed };
}

/**
 * Intersect the surface with the plane axis = value. Returns polylines in the plane's 2-D coordinates:
 * x-plane → (y, z), y-plane → (x, z), z-plane → (x, y). Closed loops repeat their first point at the end.
 * Vertices lying exactly on the plane are treated as just above it (just below it if that yields nothing),
 * so the chaining is purely topological.
 */
export function slice(model, { axis = 'x', value = 0 } = {}) {
  const P = model.positions, T = model.triangles, ax = AX.indexOf(axis); if (ax < 0) throw new Error(`slice axis must be 'x', 'y' or 'z'`);
  const [iu, iv] = ax === 0 ? [1, 2] : ax === 1 ? [0, 2] : [0, 1], nV = P.length / 3, { rep } = weldMap(P, 0);
  let loops = cutPlane(P, T, rep, ax, iu, iv, value, true);
  if (!loops.length) loops = cutPlane(P, T, rep, ax, iu, iv, value, false);   // plane lying exactly on the lowest face: take on-plane vertices as below instead
  const bb = bboxOf(P), eps = 1e-12 * len(bb.size);
  // drop zero-length steps created where the plane passes exactly through vertices
  return loops.map((l) => l.filter((q, i) => i === 0 || Math.hypot(q[0] - l[i - 1][0], q[1] - l[i - 1][1]) > eps)).filter((l) => l.length >= 2);
}
function cutPlane(P, T, rep, ax, iu, iv, value, onPlaneIsAbove) {
  const nV = P.length / 3, pts = new Map(), segs = [], above = (i) => { const d = P[3 * i + ax] - value; return onPlaneIsAbove ? d >= 0 : d > 0; };
  const cut = (a, b) => {
    const k = a < b ? a * nV + b : b * nV + a; let id = pts.get(k);
    if (id === undefined) { const da = P[3 * a + ax] - value, db = P[3 * b + ax] - value, t = da / (da - db); id = { u: P[3 * a + iu] + t * (P[3 * b + iu] - P[3 * a + iu]), v: P[3 * a + iv] + t * (P[3 * b + iv] - P[3 * a + iv]), s: [] }; pts.set(k, id); }
    return id;
  };
  for (let t = 0; t < T.length; t += 3) {
    const v = [rep[T[t]], rep[T[t + 1]], rep[T[t + 2]]]; if (v[0] === v[1] || v[1] === v[2] || v[0] === v[2]) continue;
    const up = v.map(above), nUp = up[0] + up[1] + up[2]; if (nUp === 0 || nUp === 3) continue;
    const ends = []; for (let e = 0; e < 3; e++) if (up[e] !== up[(e + 1) % 3]) ends.push(cut(v[e], v[(e + 1) % 3]));
    const sg = { a: ends[0], b: ends[1], done: false }; ends[0].s.push(sg); ends[1].s.push(sg); segs.push(sg);
  }
  const loops = [], walk = (start, sg) => {
    const line = [[start.u, start.v]]; let p = start;
    for (let guardN = 0; guardN <= segs.length && sg; guardN++) { sg.done = true; p = sg.a === p ? sg.b : sg.a; line.push([p.u, p.v]); sg = p.s.find((x) => !x.done); }
    return line;
  };
  for (const p of pts.values()) if (p.s.length % 2 === 1) for (const sg of p.s) if (!sg.done) loops.push(walk(p, sg));   // open polylines start at free ends
  for (const sg of segs) if (!sg.done) loops.push(walk(sg.a, sg));
  return loops;
}

// ---------- ASCII writers ----------
const GMSH_ID = { line2: 1, tri3: 2, quad4: 3, tet4: 4, hex8: 5, wedge6: 6, pyr5: 7, line3: 8, tri6: 9, quad9: 10, tet10: 11, hex27: 12, wedge18: 13, pyr14: 14, quad8: 16, hex20: 17, wedge15: 18, pyr13: 19 };
const VTK_ID = { line2: 3, tri3: 5, quad4: 9, tet4: 10, hex8: 12, wedge6: 13, pyr5: 14, line3: 21, tri6: 22, quad8: 23, tet10: 24, hex20: 25, wedge15: 26, pyr13: 27, quad9: 28, hex27: 29 };
/** Element blocks to write: the explicit elements, or the display triangles when there are none. */
function blocksOf(m) { return m.elements.length ? m.elements : m.triangles.length ? [{ type: 'tri3', nodesPer: 3, count: m.triangles.length / 3, conn: m.triangles }] : []; }

/** Serialise a model as ASCII: 'stl' | 'obj' | 'vtk' (legacy) | 'msh' (Gmsh 2.2) | 'su2' | 'json' (lossless platform format). */
export function exportModel(model, format) {
  const P = model.positions, T = model.triangles, nV = P.length / 3, out = [], name = String(model.name || 'model').replace(/[\r\n]+/g, ' ');
  const xyz = (i) => `${P[3 * i]} ${P[3 * i + 1]} ${P[3 * i + 2]}`;
  switch (format) {
    case 'stl': {
      const sn = name.replace(/\s+/g, '_'); out.push(`solid ${sn}`);
      for (let t = 0; t < T.length; t += 3) {
        const a = pt(P, T[t]), n = cross(sub(pt(P, T[t + 1]), a), sub(pt(P, T[t + 2]), a)), l = len(n) || 1;
        out.push(`facet normal ${n[0] / l} ${n[1] / l} ${n[2] / l}`, ' outer loop', `  vertex ${xyz(T[t])}`, `  vertex ${xyz(T[t + 1])}`, `  vertex ${xyz(T[t + 2])}`, ' endloop', 'endfacet');
      }
      out.push(`endsolid ${sn}`); break;
    }
    case 'obj': {
      out.push(`# ${name}`, `o ${name.replace(/\s+/g, '_')}`);
      for (let i = 0; i < nV; i++) out.push(`v ${xyz(i)}`);
      const tg = triangleGroups(model); let cur = -2;
      for (let t = 0; t < T.length; t += 3) { if (tg && tg[t / 3] !== cur) { cur = tg[t / 3]; out.push(`g ${(model.groups[cur]?.name ?? 'ungrouped').replace(/\s+/g, '_')}`); } out.push(`f ${T[t] + 1} ${T[t + 1] + 1} ${T[t + 2] + 1}`); }
      for (let i = 0; i < model.lines.length; i += 2) out.push(`l ${model.lines[i] + 1} ${model.lines[i + 1] + 1}`);
      break;
    }
    case 'vtk': {
      out.push('# vtk DataFile Version 3.0', name.slice(0, 250), 'ASCII');
      const blocks = blocksOf(model).filter((b) => VTK_ID[b.type]), poly = blocks.every((b) => b.type === 'tri3');
      out.push(`DATASET ${poly ? 'POLYDATA' : 'UNSTRUCTURED_GRID'}`, `POINTS ${nV} double`);
      for (let i = 0; i < nV; i++) out.push(xyz(i));
      const row = (b, e) => { let s = String(b.nodesPer); for (let k = 0; k < b.nodesPer; k++) s += ' ' + b.conn[e * b.nodesPer + k]; return s; };
      const nc = blocks.reduce((s, b) => s + b.count, 0), sz = blocks.reduce((s, b) => s + b.count * (b.nodesPer + 1), 0);
      if (poly) { if (nc) { out.push(`POLYGONS ${nc} ${sz}`); for (const b of blocks) for (let e = 0; e < b.count; e++) out.push(row(b, e)); } else if (nV) { out.push(`VERTICES ${nV} ${2 * nV}`); for (let i = 0; i < nV; i++) out.push(`1 ${i}`); } }
      else { out.push(`CELLS ${nc} ${sz}`); for (const b of blocks) for (let e = 0; e < b.count; e++) out.push(row(b, e)); out.push(`CELL_TYPES ${nc}`); for (const b of blocks) for (let e = 0; e < b.count; e++) out.push(String(VTK_ID[b.type])); }
      break;
    }
    case 'msh': {
      const blocks = blocksOf(model).filter((b) => GMSH_ID[b.type]), gdim = new Map();
      for (const b of blocks) if (b.group) for (const g of b.group) if (g >= 0 && !gdim.has(g)) gdim.set(g, ET[b.type].dim);
      out.push('$MeshFormat', '2.2 0 8', '$EndMeshFormat');
      if (gdim.size) { out.push('$PhysicalNames', String(gdim.size)); for (const [g, d] of gdim) out.push(`${d} ${g + 1} "${String(model.groups[g]?.name ?? `group ${g}`).replace(/"/g, "'")}"`); out.push('$EndPhysicalNames'); }
      out.push('$Nodes', String(nV)); for (let i = 0; i < nV; i++) out.push(`${i + 1} ${xyz(i)}`); out.push('$EndNodes');
      out.push('$Elements', String(blocks.reduce((s, b) => s + b.count, 0))); let id = 0;
      for (const b of blocks) for (let e = 0; e < b.count; e++) { const g = b.group && b.group[e] >= 0 ? b.group[e] + 1 : 0; let s = `${++id} ${GMSH_ID[b.type]} 2 ${g} ${g}`; for (let k = 0; k < b.nodesPer; k++) s += ' ' + (b.conn[e * b.nodesPer + k] + 1); out.push(s); }
      out.push('$EndElements'); break;
    }
    case 'su2': {
      const all = blocksOf(model).filter((b) => VTK_ID[b.type] && ET[b.type].n === ET[b.type].c), vol = all.filter((b) => ET[b.type].dim === 3), surf = all.filter((b) => ET[b.type].dim === 2), line = all.filter((b) => ET[b.type].dim === 1);
      let planar = !vol.length && surf.length > 0; if (planar) for (let i = 2; i < P.length; i += 3) if (P[i] !== 0) { planar = false; break; }
      const interior = vol.length ? vol : surf, markerBlocks = vol.length ? surf : planar ? line : [];
      const row = (b, e) => { let s = String(VTK_ID[b.type]); for (let k = 0; k < b.nodesPer; k++) s += ' ' + b.conn[e * b.nodesPer + k]; return s; };
      out.push(`% ${name}`, `NDIME= ${planar ? 2 : 3}`, `NELEM= ${interior.reduce((s, b) => s + b.count, 0)}`);
      let k = 0; for (const b of interior) for (let e = 0; e < b.count; e++) out.push(`${row(b, e)} ${k++}`);
      out.push(`NPOIN= ${nV}`); for (let i = 0; i < nV; i++) out.push(planar ? `${P[3 * i]} ${P[3 * i + 1]} ${i}` : `${xyz(i)} ${i}`);
      const markers = new Map();
      for (const b of markerBlocks) for (let e = 0; e < b.count; e++) { const g = b.group ? b.group[e] : -1, nm = g >= 0 ? String(model.groups[g]?.name ?? `marker_${g}`).replace(/\s+/g, '_') : 'unassigned'; if (!markers.has(nm)) markers.set(nm, []); markers.get(nm).push(row(b, e)); }
      if (vol.length && !markers.size && T.length) { const rows = []; for (let t = 0; t < T.length; t += 3) rows.push(`5 ${T[t]} ${T[t + 1]} ${T[t + 2]}`); markers.set('boundary', rows); }
      out.push(`NMARK= ${markers.size}`); for (const [nm, rows] of markers) { out.push(`MARKER_TAG= ${nm}`, `MARKER_ELEMS= ${rows.length}`); for (const r of rows) out.push(r); }
      break;
    }
    case 'json':
      return JSON.stringify({ aerosuiteGeometry: 1, name: model.name, format: model.format, formatName: model.formatName, kind: model.kind, units: model.units, positions: Array.from(P), triangles: Array.from(T), lines: Array.from(model.lines), elements: model.elements.map((e) => ({ type: e.type, nodesPer: e.nodesPer, count: e.count, conn: Array.from(e.conn), ...(e.group ? { group: Array.from(e.group) } : {}) })), groups: model.groups, meta: model.meta, log: model.log, warnings: model.warnings });
    default: throw new Error(`exportModel: unknown format "${format}" (use stl, obj, vtk, msh, su2 or json)`);
  }
  return out.join('\n') + '\n';
}
