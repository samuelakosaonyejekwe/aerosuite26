// Planar sections: aerofoil-style measures of a slice loop, a quality 2-D triangulation of the enclosed region,
// section properties (area, second moments, Saint-Venant torsion constant by finite elements) and a
// mesh-convergence study on the torsion constant.

import * as N from '../numerics.js';

const area2 = (a, b, c) => (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]);
const shoelace = (l) => { let s = 0; for (let i = 0, n = l.length; i < n; i++) { const a = l[i], b = l[(i + 1) % n]; s += a[0] * b[1] - b[0] * a[1]; } return 0.5 * s; };

/** Remove the repeated closing point and consecutive duplicates; returns { pts, closedExplicitly, gap }. */
function cleanLoop(loop) {
  const src = (loop || []).filter((p) => Array.isArray(p) && Number.isFinite(p[0]) && Number.isFinite(p[1])).map((p) => [p[0], p[1]]);
  if (src.length < 2) return { pts: src, closedExplicitly: false, gap: 0 };
  let size = 0; for (const p of src) size = Math.max(size, Math.abs(p[0] - src[0][0]), Math.abs(p[1] - src[0][1]));
  const eps = 1e-9 * (size || 1), gap = dist(src[0], src[src.length - 1]), pts = [src[0]];
  for (let i = 1; i < src.length; i++) if (dist(src[i], pts[pts.length - 1]) > eps) pts.push(src[i]);
  const closedExplicitly = pts.length > 2 && dist(pts[0], pts[pts.length - 1]) <= eps;
  if (closedExplicitly) pts.pop();
  return { pts, closedExplicitly, gap };
}

/**
 * Aerofoil-like measures of one section loop. The chord is the longest distance between two points of the loop;
 * the end with the larger u is taken as the trailing edge (x-aft convention). Thickness and camber are measured
 * perpendicular to the chord line at 400 stations. Lengths are in the loop's units; tc, camberRatio and
 * xThickness (chordwise position of maximum thickness) are fractions of chord.
 */
export function sectionMetrics(loop) {
  const { pts, closedExplicitly, gap } = cleanLoop(loop), n = pts.length;
  const out = { chord: 0, thickness: 0, tc: 0, camber: 0, camberRatio: 0, xThickness: 0, area: 0, perimeter: 0, le: [0, 0], te: [0, 0], closed: closedExplicitly, gap };
  if (n < 2) return out;
  // farthest pair (on a subsample for very long loops)
  const st = Math.max(1, Math.floor(n / 2500)); let best = -1, ia = 0, ib = 0;
  for (let i = 0; i < n; i += st) for (let j = i + st; j < n; j += st) { const d = (pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2; if (d > best) { best = d; ia = i; ib = j; } }
  let le = pts[ia], te = pts[ib]; if (le[0] > te[0] || (le[0] === te[0] && le[1] > te[1])) [le, te] = [te, le];
  const c = dist(le, te); out.chord = c; out.le = le.slice(); out.te = te.slice();
  for (let i = 0; i < n - (closedExplicitly || n > 2 ? 0 : 1); i++) out.perimeter += dist(pts[i], pts[(i + 1) % n]);
  if (n < 3 || !(c > 0)) return out;
  out.area = Math.abs(shoelace(pts));
  const ex = [(te[0] - le[0]) / c, (te[1] - le[1]) / c], q = pts.map((p) => { const dx = p[0] - le[0], dy = p[1] - le[1]; return [dx * ex[0] + dy * ex[1], -dx * ex[1] + dy * ex[0]]; });
  const M = 400; let tmax = 0, xt = 0, cmax = 0;
  for (let k = 0; k < M; k++) {
    const x = (c * (k + 0.5)) / M; let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < n; i++) { const a = q[i], b = q[(i + 1) % n]; if ((a[0] - x) * (b[0] - x) > 0 || a[0] === b[0]) continue; const y = a[1] + ((x - a[0]) / (b[0] - a[0])) * (b[1] - a[1]); if (y < lo) lo = y; if (y > hi) hi = y; }
    if (!(hi >= lo)) continue;
    if (hi - lo > tmax) { tmax = hi - lo; xt = x / c; }
    if (Math.abs(0.5 * (hi + lo)) > Math.abs(cmax)) cmax = 0.5 * (hi + lo);
  }
  out.thickness = tmax; out.tc = tmax / c; out.camber = cmax; out.camberRatio = cmax / c; out.xThickness = xt;
  return out;
}

/** Resample a closed polygon at spacing ≈ h along its own edges, keeping vertices where the outline turns sharply. */
function resample(pts, h) {
  const n = pts.length, corners = [];
  for (let i = 0; i < n; i++) {
    const a = pts[(i + n - 1) % n], b = pts[i], c = pts[(i + 1) % n], u = [b[0] - a[0], b[1] - a[1]], v = [c[0] - b[0], c[1] - b[1]];
    const ang = Math.abs(Math.atan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1]));
    if (ang > 0.5) corners.push(i);                         // > ~29° turn: a real corner (trailing edge, plate corner …)
  }
  const out = [], span = (i0, i1, closedSpan) => {
    const path = [pts[i0 % n]]; for (let i = i0 + 1; i <= i1; i++) path.push(pts[i % n]);
    let L = 0; const s = [0]; for (let i = 1; i < path.length; i++) { L += dist(path[i], path[i - 1]); s.push(L); }
    const m = Math.max(closedSpan ? 3 : 1, Math.round(L / h));
    for (let k = 0, j = 1; k < m; k++) { const t = (k * L) / m; while (j < s.length - 1 && s[j] < t) j++; const f = s[j] > s[j - 1] ? (t - s[j - 1]) / (s[j] - s[j - 1]) : 0; out.push([path[j - 1][0] + f * (path[j][0] - path[j - 1][0]), path[j - 1][1] + f * (path[j][1] - path[j - 1][1])]); }
  };
  if (!corners.length) span(0, n, true);
  else for (let k = 0; k < corners.length; k++) span(corners[k], k + 1 < corners.length ? corners[k + 1] : corners[0] + n, false);
  return out;
}

/** Even-odd point-in-region test over a set of closed polygons. */
function inside(loops, x, y) {
  let c = false;
  for (const l of loops) for (let i = 0, n = l.length, j = n - 1; i < n; j = i++) { const a = l[i], b = l[j]; if ((a[1] > y) !== (b[1] > y) && x < a[0] + ((y - a[1]) / (b[1] - a[1])) * (b[0] - a[0])) c = !c; }
  return c;
}

/**
 * Delaunay triangulation (Bowyer–Watson) of points sorted in x: triangles whose circumcircle lies wholly to the
 * left of the sweep position are retired, so each insertion only tests the active front.
 */
function delaunay(P) {
  const n = P.length; if (n < 3) return [];
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity; for (const p of P) { if (p[0] < x0) x0 = p[0]; if (p[0] > x1) x1 = p[0]; if (p[1] < y0) y0 = p[1]; if (p[1] > y1) y1 = p[1]; }
  const d = Math.max(x1 - x0, y1 - y0) || 1, mx = 0.5 * (x0 + x1), my = 0.5 * (y0 + y1);
  const X = new Float64Array(n + 3), Y = new Float64Array(n + 3);
  for (let i = 0; i < n; i++) { X[i] = (P[i][0] - mx) / d; Y[i] = (P[i][1] - my) / d; }                // normalised to a unit box
  X[n] = -30; Y[n] = -15; X[n + 1] = 0; Y[n + 1] = 30; X[n + 2] = 30; Y[n + 2] = -15;                    // super-triangle
  const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => X[a] - X[b] || Y[a] - Y[b]);
  const circ = (i, j, k) => {
    const ax = X[i], ay = Y[i], bx = X[j] - ax, by = Y[j] - ay, cx = X[k] - ax, cy = Y[k] - ay, D = 2 * (bx * cy - by * cx);
    if (D === 0) return { i, j, k, x: 0, y: 0, r2: Infinity };
    const b2 = bx * bx + by * by, c2 = cx * cx + cy * cy, ux = (cy * b2 - by * c2) / D, uy = (bx * c2 - cx * b2) / D;
    return { i, j, k, x: ax + ux, y: ay + uy, r2: ux * ux + uy * uy };
  };
  let open = [circ(n, n + 1, n + 2)]; const done = [], edges = [];
  for (const i of order) {
    const px = X[i], py = Y[i]; edges.length = 0;
    for (let j = open.length - 1; j >= 0; j--) {
      const c = open[j], dx = px - c.x;
      if (dx > 0 && dx * dx > c.r2) { done.push(c); open[j] = open[open.length - 1]; open.pop(); continue; }
      const dy = py - c.y; if (dx * dx + dy * dy >= c.r2) continue;
      edges.push(c.i, c.j, c.j, c.k, c.k, c.i); open[j] = open[open.length - 1]; open.pop();
    }
    // cavity boundary = edges that occur once
    const seen = new Map();
    for (let e = 0; e < edges.length; e += 2) { const a = edges[e], b = edges[e + 1], k = a < b ? a * (n + 3) + b : b * (n + 3) + a; seen.set(k, seen.has(k) ? null : [a, b]); }
    for (const e of seen.values()) if (e) open.push(circ(e[0], e[1], i));
  }
  const out = [];
  for (const c of done.concat(open)) if (c.i < n && c.j < n && c.k < n) out.push([c.i, c.j, c.k]);
  return out;
}

function triQuality(nodes, tris) {
  let minAngle = 180, asp = 0;
  for (const [i, j, k] of tris) {
    const a = nodes[i], b = nodes[j], c = nodes[k], la = dist(b, c), lb = dist(a, c), lc = dist(a, b), A = 0.5 * Math.abs(area2(a, b, c));
    const ang = (o, p, q) => (Math.acos(Math.max(-1, Math.min(1, (p * p + q * q - o * o) / (2 * p * q || 1e-300)))) * 180) / Math.PI;
    minAngle = Math.min(minAngle, ang(la, lb, lc), ang(lb, la, lc), ang(lc, la, lb));
    asp += A > 0 ? (Math.max(la, lb, lc) * (la + lb + lc)) / (4 * Math.sqrt(3) * A) : 1e9;
  }
  return { minAngle: tris.length ? minAngle : 0, meanAspect: tris.length ? asp / tris.length : 0, nTris: tris.length };
}

/**
 * Triangulate the region enclosed by `loops` (even-odd rule: a loop inside another is a hole; input winding is
 * irrelevant). Boundary points are the loops resampled at spacing h (sharp corners kept), interior points a
 * triangular lattice of pitch h; the Delaunay triangles whose centroid lies in the region are kept and the
 * interior nodes are Laplacian-smoothed. `quality.areaError` is the relative difference between the meshed
 * area and the polygon area - a non-zero value means h is too coarse for a thin or sharp feature.
 */
export function meshSection(loops, { h, smooth = 3 } = {}) {
  const polys = (Array.isArray(loops) && loops.length && typeof loops[0]?.[0] === 'number' ? [loops] : loops || []).map((l) => cleanLoop(l).pts).filter((l) => l.length >= 3 && Math.abs(shoelace(l)) > 0);
  const empty = { nodes: [], tris: [], quality: { minAngle: 0, meanAspect: 0, nTris: 0, areaError: 1, h: h ?? 0, boundaryNodes: 0 } };
  if (!polys.length) return empty;
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity; for (const l of polys) for (const p of l) { if (p[0] < x0) x0 = p[0]; if (p[0] > x1) x1 = p[0]; if (p[1] < y0) y0 = p[1]; if (p[1] > y1) y1 = p[1]; }
  // net area by the even-odd rule: outermost loops add, loops nested at odd depth subtract
  const depth = polys.map((l, i) => polys.reduce((dp, o, j) => dp + (j !== i && inside([o], l[0][0], l[0][1]) ? 1 : 0), 0));
  const polyArea = polys.reduce((s, l, i) => s + (depth[i] % 2 ? -1 : 1) * Math.abs(shoelace(l)), 0);
  if (!(polyArea > 0)) return empty;
  let hh = Number.isFinite(h) && h > 0 ? h : Math.sqrt(polyArea) / 16;
  const cap = 250000; if (polyArea / (0.433 * hh * hh) > cap) hh = Math.sqrt(polyArea / (0.433 * cap));   // keep the triangle count bounded
  const nodes = []; for (const l of polys) for (const p of resample(l, hh)) nodes.push(p);
  const nb = nodes.length;
  // interior lattice by scanline, kept clear of the boundary samples
  const cell = hh, gx = (v) => Math.floor((v - x0) / cell), grid = new Map(), gk = (i, j) => i * 2097152 + j;
  for (let i = 0; i < nb; i++) { const k = gk(gx(nodes[i][0]), Math.floor((nodes[i][1] - y0) / cell)); (grid.get(k) || grid.set(k, []).get(k)).push(i); }
  const clear2 = (0.62 * hh) ** 2, dy = hh * Math.sqrt(3) / 2;
  for (let r = 0, y = y0 + 0.5 * dy; y < y1; r++, y += dy) {
    const xs = []; for (const l of polys) for (let i = 0, n = l.length, j = n - 1; i < n; j = i++) { const a = l[i], b = l[j]; if ((a[1] > y) !== (b[1] > y)) xs.push(a[0] + ((y - a[1]) / (b[1] - a[1])) * (b[0] - a[0])); }
    xs.sort((a, b) => a - b);
    for (let s = 0; s + 1 < xs.length; s += 2) {
      const first = Math.ceil((xs[s] - x0 - (r % 2 ? 0.5 * hh : 0)) / hh);
      for (let k = first; ; k++) {
        const x = x0 + (r % 2 ? 0.5 * hh : 0) + k * hh; if (x >= xs[s + 1]) break;
        const ci = gx(x), cj = Math.floor((y - y0) / cell); let ok = true;
        for (let a = ci - 1; a <= ci + 1 && ok; a++) for (let b = cj - 1; b <= cj + 1 && ok; b++) for (const q of grid.get(gk(a, b)) || []) if ((nodes[q][0] - x) ** 2 + (nodes[q][1] - y) ** 2 < clear2) { ok = false; break; }
        if (ok) nodes.push([x, y]);
      }
    }
  }
  // triangulate a copy perturbed by ~1e-7 h so that cocircular/collinear ties cannot occur; the mesh keeps the true coordinates
  let seed = 987654321; const rnd = () => ((seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff - 0.5);
  const tris = [];
  for (const t of delaunay(nodes.map((p) => [p[0] + 2e-7 * hh * rnd(), p[1] + 2e-7 * hh * rnd()]))) {
    const a = nodes[t[0]], b = nodes[t[1]], c = nodes[t[2]], A2 = area2(a, b, c);
    if (Math.abs(A2) < 1e-9 * hh * hh || !inside(polys, (a[0] + b[0] + c[0]) / 3, (a[1] + b[1] + c[1]) / 3)) continue;
    tris.push(A2 > 0 ? t : [t[0], t[2], t[1]]);
  }
  // Laplacian smoothing of interior nodes (a move is rejected if it would invert a neighbouring triangle)
  const nbrs = nodes.map(() => new Set()), inc = nodes.map(() => []);
  tris.forEach((t, ti) => { for (let k = 0; k < 3; k++) { nbrs[t[k]].add(t[(k + 1) % 3]); nbrs[t[k]].add(t[(k + 2) % 3]); inc[t[k]].push(ti); } });
  for (let it = 0; it < smooth; it++) for (let i = nb; i < nodes.length; i++) {
    if (!nbrs[i].size) continue;
    let sx = 0, sy = 0; for (const j of nbrs[i]) { sx += nodes[j][0]; sy += nodes[j][1]; }
    const old = nodes[i]; nodes[i] = [sx / nbrs[i].size, sy / nbrs[i].size];
    for (const ti of inc[i]) { const t = tris[ti]; if (area2(nodes[t[0]], nodes[t[1]], nodes[t[2]]) <= 1e-9 * hh * hh) { nodes[i] = old; break; } }
  }
  // drop nodes that ended up in no triangle
  const map = new Int32Array(nodes.length).fill(-1), outNodes = []; let usedB = 0;
  for (const t of tris) for (const v of t) if (map[v] < 0) { map[v] = 0; }
  for (let i = 0; i < nodes.length; i++) if (map[i] === 0) { map[i] = outNodes.length; outNodes.push(nodes[i]); if (i < nb) usedB++; }
  const outTris = tris.map((t) => [map[t[0]], map[t[1]], map[t[2]]]);
  let meshArea = 0; for (const t of outTris) meshArea += 0.5 * area2(outNodes[t[0]], outNodes[t[1]], outNodes[t[2]]);
  return { nodes: outNodes, tris: outTris, quality: { ...triQuality(outNodes, outTris), areaError: Math.abs(meshArea - polyArea) / polyArea, h: hh, boundaryNodes: usedB } };
}

/** Jacobi-preconditioned conjugate gradients on a CSR matrix. */
function cgSolve(n, rowPtr, col, val, b) {
  const x = new Float64Array(n), r = Float64Array.from(b), z = new Float64Array(n), p = new Float64Array(n), q = new Float64Array(n), dinv = new Float64Array(n);
  for (let i = 0; i < n; i++) for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) if (col[k] === i) dinv[i] = 1 / val[k];
  let rz = 0, b2 = 0; for (let i = 0; i < n; i++) { z[i] = r[i] * dinv[i]; p[i] = z[i]; rz += r[i] * z[i]; b2 += b[i] * b[i]; }
  let iters = 0, res = Math.sqrt(b2);
  for (; iters < 20 * n + 100 && b2 > 0; iters++) {
    let pq = 0; for (let i = 0; i < n; i++) { let s = 0; for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) s += val[k] * p[col[k]]; q[i] = s; pq += p[i] * s; }
    if (!(pq > 0)) break;
    const al = rz / pq; let r2 = 0, rzn = 0;
    for (let i = 0; i < n; i++) { x[i] += al * p[i]; r[i] -= al * q[i]; r2 += r[i] * r[i]; z[i] = r[i] * dinv[i]; rzn += r[i] * z[i]; }
    res = Math.sqrt(r2); if (r2 <= 1e-24 * b2) { iters++; break; }
    const be = rzn / rz; rz = rzn; for (let i = 0; i < n; i++) p[i] = z[i] + be * p[i];
  }
  return { x, iters, residual: b2 > 0 ? res / Math.sqrt(b2) : 0 };
}

/**
 * Section properties of a triangulated region. Area and second moments are exact for the mesh polygon.
 * J is the Saint-Venant torsion constant from the Prandtl stress function on linear triangles:
 * ∇²φ = −2 in the section, φ = 0 on the outer boundary, J = 2∫φ dA. For hollow (multiply connected) sections
 * φ is held at one unknown constant per hole, with the hole area entering the load and J (J += 2·φ_k·A_k);
 * the returned `note` states which treatment was used. Linear elements make J converge from below at O(h²).
 * Ixx = ∫y² dA and Iyy = ∫x² dA about the centroid; theta_p (rad, from +x) is the direction of the axis about
 * which the second moment is I1, the larger principal value.
 */
export function sectionProperties(mesh) {
  const nodes = mesh?.nodes || [], nn = nodes.length, tris = (mesh?.tris || []).map((t) => (area2(nodes[t[0]], nodes[t[1]], nodes[t[2]]) >= 0 ? t : [t[0], t[2], t[1]]));
  const zero = { A: 0, cx: 0, cy: 0, Ixx: 0, Iyy: 0, Ixy: 0, I1: 0, I2: 0, theta_p: 0, J: 0, dof: 0, note: 'empty mesh' };
  if (!tris.length) return zero;
  let A = 0, Sx = 0, Sy = 0, Jxx = 0, Jyy = 0, Jxy = 0;
  for (const [i, j, k] of tris) {
    const [x1, y1] = nodes[i], [x2, y2] = nodes[j], [x3, y3] = nodes[k], a = 0.5 * area2(nodes[i], nodes[j], nodes[k]), sx = x1 + x2 + x3, sy = y1 + y2 + y3;
    A += a; Sx += (a * sx) / 3; Sy += (a * sy) / 3;
    Jxx += (a / 12) * (y1 * y1 + y2 * y2 + y3 * y3 + sy * sy); Jyy += (a / 12) * (x1 * x1 + x2 * x2 + x3 * x3 + sx * sx); Jxy += (a / 12) * (x1 * y1 + x2 * y2 + x3 * y3 + sx * sy);
  }
  if (!(A > 0)) return zero;
  const cx = Sx / A, cy = Sy / A, Ixx = Jxx - A * cy * cy, Iyy = Jyy - A * cx * cx, Ixy = Jxy - A * cx * cy;
  const avg = 0.5 * (Ixx + Iyy), rad = Math.hypot(0.5 * (Ixx - Iyy), Ixy), theta_p = 0.5 * Math.atan2(-2 * Ixy, Ixx - Iyy);
  // boundary loops of the mesh: edges used by exactly one triangle
  const ecount = new Map(), ek = (a, b) => (a < b ? a * nn + b : b * nn + a);
  for (const t of tris) for (let e = 0; e < 3; e++) { const k = ek(t[e], t[(e + 1) % 3]); ecount.set(k, (ecount.get(k) || 0) + 1); }
  const next = new Map(); for (const t of tris) for (let e = 0; e < 3; e++) { const a = t[e], b = t[(e + 1) % 3]; if (ecount.get(ek(a, b)) === 1) (next.get(a) || next.set(a, []).get(a)).push(b); }
  const bloops = [], visited = new Set();
  for (const start of next.keys()) {
    if (visited.has(start)) continue;
    const loop = []; let v = start;
    for (let g = 0; g <= nn && v !== undefined && !visited.has(v); g++) { visited.add(v); loop.push(v); v = (next.get(v) || []).find((w) => !visited.has(w)); }
    bloops.push(loop);
  }
  const areas = bloops.map((l) => Math.abs(shoelace(l.map((i) => nodes[i])))), outer = areas.indexOf(Math.max(...areas));
  // degrees of freedom: interior nodes, plus one per hole
  const dofOf = new Int32Array(nn).fill(-2); let nd = 0; const holes = [];
  bloops.forEach((l, b) => { if (b === outer) { for (const i of l) dofOf[i] = -1; } else { const d = nd++; holes.push({ dof: d, area: areas[b] }); for (const i of l) if (dofOf[i] === -2) dofOf[i] = d; } });
  for (let i = 0; i < nn; i++) if (dofOf[i] === -2) dofOf[i] = nd++;
  let J = 0, solver = { iters: 0, residual: 0 };
  if (nd > 0) {
    const rows = Array.from({ length: nd }, () => new Map()), f = new Float64Array(nd);
    for (const h of holes) f[h.dof] += 2 * h.area;
    for (const t of tris) {
      const p = t.map((i) => nodes[i]), a = 0.5 * area2(p[0], p[1], p[2]); if (!(a > 0)) continue;
      const bq = [p[1][1] - p[2][1], p[2][1] - p[0][1], p[0][1] - p[1][1]], cq = [p[2][0] - p[1][0], p[0][0] - p[2][0], p[1][0] - p[0][0]];
      for (let i = 0; i < 3; i++) { const di = dofOf[t[i]]; if (di < 0) continue; f[di] += (2 * a) / 3; for (let j = 0; j < 3; j++) { const dj = dofOf[t[j]]; if (dj < 0) continue; rows[di].set(dj, (rows[di].get(dj) || 0) + (bq[i] * bq[j] + cq[i] * cq[j]) / (4 * a)); } }
    }
    const rowPtr = new Int32Array(nd + 1); for (let i = 0; i < nd; i++) rowPtr[i + 1] = rowPtr[i] + rows[i].size;
    const col = new Int32Array(rowPtr[nd]), val = new Float64Array(rowPtr[nd]);
    for (let i = 0, k = 0; i < nd; i++) for (const [c, v] of rows[i]) { col[k] = c; val[k++] = v; }
    solver = cgSolve(nd, rowPtr, col, val, f);
    const phi = (i) => (dofOf[i] < 0 ? 0 : solver.x[dofOf[i]]);
    for (const t of tris) J += (2 * 0.5 * area2(nodes[t[0]], nodes[t[1]], nodes[t[2]]) * (phi(t[0]) + phi(t[1]) + phi(t[2]))) / 3;
    for (const h of holes) J += 2 * solver.x[h.dof] * h.area;
  }
  const note = holes.length
    ? `Hollow section with ${holes.length} hole(s): the stress function is held at an unknown constant on each hole boundary and the hole areas enter the load and J (multiply-connected Prandtl formulation); no thin-wall simplification is made. Linear triangles: J converges from below at O(h²).`
    : 'Solid section: Prandtl stress function with φ = 0 on the boundary, linear triangles (J converges from below at O(h²)).';
  return { A, cx, cy, Ixx, Iyy, Ixy, I1: avg + rad, I2: avg - rad, theta_p, J, dof: nd, note, holes: holes.length, cgIterations: solver.iters, cgResidual: solver.residual };
}

/** Mesh-convergence study of J: one row per element size (sorted fine → coarse) and the GCI of the three finest. */
export function sectionConvergence(loops, hs) {
  const rows = [...hs].filter((h) => Number.isFinite(h) && h > 0).sort((a, b) => a - b).map((h) => { const m = meshSection(loops, { h }), p = sectionProperties(m); return { h, nTris: m.tris.length, J: p.J, A: p.A }; });
  const gci = rows.length >= 3 ? N.gci(rows.slice(0, 3).map((r) => r.h), rows.slice(0, 3).map((r) => r.J)) : null;
  return { rows, gci };
}
