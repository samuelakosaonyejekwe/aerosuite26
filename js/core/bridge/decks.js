// High-fidelity bridge — deck generators. Turns the shared case (and optionally an imported geometry
// or mesh model) into complete, ready-to-run input decks for open-source solvers:
//   SU2 (.cfg + .su2 mesh), OpenFOAM (case directory), CalculiX / Abaqus (.inp), Gmsh (.geo).
// Pure computation: no DOM, no network, no file access. Runs unchanged in the browser and in Node.
// SI units throughout (m, kg, s, N, Pa, K). Axes: x aft, y spanwise (starboard), z up.

import { isa, GAMMA, R_AIR, G0 } from '../atmosphere.js';
import { derived } from '../case.js';
import { METALS } from '../../data/materials.js';
import { exportModel } from '../geometry/analysis.js';

const RAD = Math.PI / 180;
const f = (v, n = 9) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '0');
const slug = (s) => String(s || 'case').toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-+|-+$/g, '').slice(0, 40) || 'case';
const ident = (s) => String(s).replace(/[^A-Za-z0-9_]+/g, '_').replace(/^_+|_+$/g, '') || 'set';

// ---------------------------------------------------------------------------------------------
// SHA-256 (synchronous, for traceability hashes in the manifest)
// ---------------------------------------------------------------------------------------------
const K256 = (() => { const k = [], isP = (n) => { for (let i = 2; i * i <= n; i++) if (n % i === 0) return false; return true; }; for (let n = 2; k.length < 64; n++) if (isP(n)) k.push(n); return k; })();
const KC = K256.map((p) => Math.floor((Math.cbrt(p) % 1) * 2 ** 32) >>> 0), HC = K256.slice(0, 8).map((p) => Math.floor((Math.sqrt(p) % 1) * 2 ** 32) >>> 0);
export function sha256(input) {
  const msg = typeof input === 'string' ? new TextEncoder().encode(input) : input, l = msg.length, n = ((l + 9 + 63) >> 6) << 6, b = new Uint8Array(n);
  b.set(msg); b[l] = 0x80;
  const dv = new DataView(b.buffer); dv.setUint32(n - 8, Math.floor((l * 8) / 2 ** 32)); dv.setUint32(n - 4, (l * 8) >>> 0);
  const H = HC.slice(), w = new Uint32Array(64), rr = (x, k) => (x >>> k) | (x << (32 - k));
  for (let o = 0; o < n; o += 64) {
    for (let i = 0; i < 16; i++) w[i] = dv.getUint32(o + 4 * i);
    for (let i = 16; i < 64; i++) { const a = w[i - 15], c = w[i - 2]; w[i] = (w[i - 16] + (rr(a, 7) ^ rr(a, 18) ^ (a >>> 3)) + w[i - 7] + (rr(c, 17) ^ rr(c, 19) ^ (c >>> 10))) >>> 0; }
    let [a, bb, c, d, e, ff, g, h] = H;
    for (let i = 0; i < 64; i++) {
      const t1 = (h + (rr(e, 6) ^ rr(e, 11) ^ rr(e, 25)) + ((e & ff) ^ (~e & g)) + KC[i] + w[i]) >>> 0, t2 = ((rr(a, 2) ^ rr(a, 13) ^ rr(a, 22)) + ((a & bb) ^ (a & c) ^ (bb & c))) >>> 0;
      h = g; g = ff; ff = e; e = (d + t1) >>> 0; d = c; c = bb; bb = a; a = (t1 + t2) >>> 0;
    }
    H[0] = (H[0] + a) >>> 0; H[1] = (H[1] + bb) >>> 0; H[2] = (H[2] + c) >>> 0; H[3] = (H[3] + d) >>> 0; H[4] = (H[4] + e) >>> 0; H[5] = (H[5] + ff) >>> 0; H[6] = (H[6] + g) >>> 0; H[7] = (H[7] + h) >>> 0;
  }
  return H.map((x) => x.toString(16).padStart(8, '0')).join('');
}
/** Deterministic JSON (sorted keys) so the same case always hashes to the same value. */
export function canonicalJson(v) {
  if (v === null || typeof v !== 'object') return JSON.stringify(v ?? null);
  if (Array.isArray(v)) return `[${v.map(canonicalJson).join(',')}]`;
  return `{${Object.keys(v).sort().filter((k) => v[k] !== undefined).map((k) => `${JSON.stringify(k)}:${canonicalJson(v[k])}`).join(',')}}`;
}
export const caseHash = (c) => sha256(canonicalJson(c));

// ---------------------------------------------------------------------------------------------
// Flow state and wall spacing
// ---------------------------------------------------------------------------------------------
/** Freestream and reference quantities derived from the shared case. */
export function flowState(c) {
  const d = derived(c), atm = isa(c.atm.alt_m, c.atm.dISA_K), V = Math.max(0.1, c.flight.V_ms);
  const hasWing = c.wing.S_m2 > 0 && c.wing.b_m > 0;
  const chord = hasWing ? d.mac : c.rotor.chord_m > 0 ? c.rotor.chord_m : 1;
  const body = !hasWing && c.fuselage.len_m > 0;
  return {
    V, alpha_deg: c.flight.alpha_deg, beta_deg: c.flight.beta_deg, T: atm.T, p: atm.p, rho: atm.rho, mu: atm.mu, nu: atm.nu, a: atm.a,
    mach: V / atm.a, chord, reynolds: (atm.rho * V * chord) / atm.mu, q: 0.5 * atm.rho * V * V,
    S: hasWing ? c.wing.S_m2 : body ? (Math.PI * c.fuselage.dia_m ** 2) / 4 : chord, b: hasWing ? c.wing.b_m : 0, hasWing,
    refLength: hasWing ? d.mac : body ? c.fuselage.len_m : chord, turb: Math.max(1e-4, c.atm.turb_intensity || 0.01), d,
  };
}
/** Height of the first node off the wall for a target y+ (flat-plate skin-friction estimate, Cf = 0.026 Re^-1/7). */
export function wallSpacing(reynolds, chord, yplus = 1) {
  const Re = Math.max(1e3, reynolds), cf = 0.026 / Re ** (1 / 7);
  return (yplus * chord) / (Re * Math.sqrt(cf / 2));
}

// ---------------------------------------------------------------------------------------------
// NACA 4-digit section
// ---------------------------------------------------------------------------------------------
export function parseNaca(code, tcFallback = 0.12) {
  const s = String(code ?? '').replace(/\D/g, '');
  if (s.length === 4 && Number(s.slice(2)) > 0) return { m: Number(s[0]) / 100, p: Number(s[1]) / 10, t: Number(s.slice(2)) / 100, code: s, exact: true };
  const t = Math.min(0.4, Math.max(0.04, tcFallback || 0.12));
  return { m: 0, p: 0, t, code: `00${String(Math.round(t * 100)).padStart(2, '0')}`, exact: false };
}
/** Upper and lower surface points of a closed-trailing-edge NACA 4-digit section at chord fraction x. */
export function nacaAt(a, x) {
  const yt = 5 * a.t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1036 * x ** 4);
  let yc = 0, dy = 0;
  if (a.m > 0 && a.p > 0) { if (x < a.p) { yc = (a.m / a.p ** 2) * (2 * a.p * x - x * x); dy = ((2 * a.m) / a.p ** 2) * (a.p - x); } else { yc = (a.m / (1 - a.p) ** 2) * (1 - 2 * a.p + 2 * a.p * x - x * x); dy = ((2 * a.m) / (1 - a.p) ** 2) * (a.p - x); } }
  const th = Math.atan(dy), s = Math.sin(th), co = Math.cos(th);
  return { upper: [x - yt * s, yc + yt * co], lower: [x + yt * s, yc - yt * co], yc, yt };
}
/** Closed loop of 2n+1 points, clockwise from the trailing edge: lower surface TE→LE, upper surface LE→TE. */
export function nacaLoop(a, n) {
  const pts = [];
  for (let i = 0; i <= n; i++) { const x = 0.5 * (1 + Math.cos((Math.PI * i) / n)); pts.push(i === 0 ? [1, 0] : i === n ? [0, 0] : nacaAt(a, x).lower); }
  for (let i = 1; i <= n; i++) { const x = 0.5 * (1 - Math.cos((Math.PI * i) / n)); pts.push(i === n ? [1, 0] : nacaAt(a, x).upper); }
  return pts;
}

// ---------------------------------------------------------------------------------------------
// Structured C- and O-grids around the section, written as SU2 native meshes
// ---------------------------------------------------------------------------------------------
/** Growth ratio r of an n-cell geometric progression whose first cell is the fraction s1 of the total. */
export function geometricRatio(s1, n) {
  if (!(s1 > 0) || n * s1 >= 1) return 1;
  let lo = 1 + 1e-12, hi = 50;
  for (let i = 0; i < 200; i++) { const r = 0.5 * (lo + hi); if ((r ** n - 1) / (r - 1) > 1 / s1) hi = r; else lo = r; }
  return 0.5 * (lo + hi);
}
const geometricNodes = (s1, n) => { const r = geometricRatio(s1, n), s = [0]; if (r === 1) { for (let j = 1; j <= n; j++) s.push(j / n); return { s, r }; } for (let j = 1; j <= n; j++) s.push((s1 * (r ** j - 1)) / (r - 1)); s[n] = 1; return { s, r }; };
const cellsForGrowth = (s1, g) => Math.max(8, Math.ceil(Math.log(1 + (g - 1) / s1) / Math.log(g)));

export const RESOLUTIONS = { coarse: { nSurf: 64, nWake: 32, growth: 1.22 }, medium: { nSurf: 100, nWake: 48, growth: 1.16 }, fine: { nSurf: 160, nWake: 72, growth: 1.12 } };

/**
 * Structured grid around a NACA 4-digit section.
 * type 'C': wake cut downstream of the trailing edge (best for viscous flow); 'O': closed rings (compact, good for Euler).
 * Lengths are in metres; firstCell is the height of the first node off the wall.
 */
export function airfoilGrid({ naca = '0012', tc = 0.12, chord = 1, type = 'C', nSurf = 100, nWake = 48, nNormal = null, growth = 1.16, firstCell = 1e-5, radius = 40, wake = 30 } = {}) {
  const sec = parseNaca(naca, tc), n = Math.max(12, Math.round(nSurf)), P = nacaLoop(sec, n), isC = type !== 'O';
  const R = Math.max(5, radius), s1 = Math.min(0.05, Math.max(1e-9, firstCell / chord / R));
  const nj = Math.max(8, Math.round(nNormal || cellsForGrowth(s1, growth))), { s: sj, r: ratio } = geometricNodes(s1, nj);
  const nRing = isC ? 2 * n + 1 : 2 * n;                               // surface nodes carrying a grid line
  // arc-length fraction
  const t = [0]; for (let i = 1; i <= 2 * n; i++) t.push(t[i - 1] + Math.hypot(P[i][0] - P[i - 1][0], P[i][1] - P[i - 1][1]));
  const arc = t[2 * n]; for (let i = 0; i <= 2 * n; i++) t[i] /= arc;
  // outward normals (clockwise traversal → outward is the tangent rotated +90°), smoothed so the fan at the trailing edge spreads
  let nrm = [];
  for (let i = 0; i < nRing; i++) {
    if (isC && i === 0) { nrm.push([0, -1]); continue; } if (isC && i === 2 * n) { nrm.push([0, 1]); continue; } if (!isC && i === 0) { nrm.push([1, 0]); continue; }
    const a = P[i - 1], b = P[i + 1], tx = b[0] - a[0], ty = b[1] - a[1], l = Math.hypot(tx, ty) || 1;
    nrm.push([-ty / l, tx / l]);
  }
  const passes = isC ? 4 : 3;
  for (let k = 0; k < passes; k++) {
    const nx = nrm.map((v) => v.slice());
    for (let i = 0; i < nRing; i++) {
      if (isC && (i === 0 || i === nRing - 1)) continue;
      const a = nrm[(i - 1 + nRing) % nRing], b = nrm[(i + 1) % nRing], x = 0.25 * a[0] + 0.5 * nrm[i][0] + 0.25 * b[0], y = 0.25 * a[1] + 0.5 * nrm[i][1] + 0.25 * b[1], l = Math.hypot(x, y) || 1;
      nx[i] = [x / l, y / l];
    }
    nrm = nx;
  }
  const C0 = isC ? [1, 0] : [0.5, 0];
  const thArc = (i) => (isC ? -Math.PI / 2 - Math.PI * t[i] : -2 * Math.PI * t[i]);
  const thRay = [];
  for (let i = 0; i < nRing; i++) {
    const ax = P[i][0] - C0[0], ay = P[i][1] - C0[1], b = ax * nrm[i][0] + ay * nrm[i][1], s = -b + Math.sqrt(Math.max(0, b * b - (ax * ax + ay * ay - R * R)));
    let th = Math.atan2(ay + s * nrm[i][1], ax + s * nrm[i][0]); const ref = thArc(i);
    th += 2 * Math.PI * Math.round((ref - th) / (2 * Math.PI)); thRay.push(th);
  }
  // wake stations (C-grid): geometric growth from the trailing-edge surface spacing
  const nW = isC ? Math.max(4, Math.round(nWake)) : 0, dTE = Math.hypot(P[1][0] - P[0][0], P[1][1] - P[0][1]);
  const xw = isC ? geometricNodes(Math.min(0.5, dTE / wake), nW).s.map((s) => 1 + s * wake) : [];
  const ni = isC ? 2 * nW + 2 * n + 1 : 2 * n;                         // grid lines in the wrap-around direction
  const build = (w) => {
    // outer-boundary angle of every grid line: follow the wall normal where the normals turn steadily, but always advance by
    // at least a share of the arc-length step so lines never bunch or cross (concave lower surface, trailing-edge closure)
    const total = isC ? Math.PI : 2 * Math.PI, last = isC ? nRing - 1 : nRing, dec = [];
    for (let i = 1; i <= last; i++) { const ray = i < nRing ? thRay[i - 1] - thRay[i] : thRay[nRing - 1] - (thRay[0] - 2 * Math.PI), arcStep = total * (t[i] - t[i - 1]); dec.push(Math.max((1 - w) * ray + w * arcStep, 0.25 * arcStep)); }
    const sum = dec.reduce((a, b) => a + b, 0), th = [isC ? -Math.PI / 2 : 0];
    for (let i = 1; i < nRing; i++) th.push(th[i - 1] - (dec[i - 1] * total) / sum);
    const X = new Float64Array(ni * (nj + 1)), Y = new Float64Array(ni * (nj + 1));
    const put = (i, j, x, y) => { X[i * (nj + 1) + j] = x; Y[i * (nj + 1) + j] = y; };
    for (let q = 0; q < nRing; q++) {
      const A = P[q], O = [C0[0] + R * Math.cos(th[q]), C0[1] + R * Math.sin(th[q])], dx = O[0] - A[0], dy = O[1] - A[1];
      const i = isC ? nW + q : q;
      for (let j = 0; j <= nj; j++) put(i, j, A[0] + sj[j] * dx, A[1] + sj[j] * dy);
    }
    for (let m = 1; m <= nW; m++) for (let j = 0; j <= nj; j++) { put(nW - m, j, xw[m], -sj[j] * R); put(nW + 2 * n + m, j, xw[m], sj[j] * R); }
    return { X, Y };
  };
  // node numbering with the wake cut (and the trailing-edge point) merged
  const nid = new Int32Array(ni * (nj + 1)); let nNodes = 0;
  for (let i = 0; i < ni; i++) for (let j = 0; j <= nj; j++) { const twin = isC && j === 0 && i >= nW + 2 * n; nid[i * (nj + 1) + j] = twin ? -1 : nNodes++; }
  if (isC) for (let i = nW + 2 * n; i < ni; i++) nid[i * (nj + 1)] = nid[(ni - 1 - i) * (nj + 1)];
  const id = (i, j) => nid[(i % ni) * (nj + 1) + j], nCellsI = isC ? ni - 1 : ni;
  const minArea = (g) => {
    let mn = Infinity;
    for (let i = 0; i < nCellsI; i++) for (let j = 0; j < nj; j++) {
      const a = i * (nj + 1) + j, b = ((i + 1) % ni) * (nj + 1) + j, c = b + 1, d = a + 1;
      const ar = 0.5 * ((g.X[c] - g.X[a]) * (g.Y[d] - g.Y[b]) - (g.X[d] - g.X[b]) * (g.Y[c] - g.Y[a]));
      if (ar < mn) mn = ar;
    }
    return mn;
  };
  let grid = null, used = null;
  for (const w of [0, 0.25, 0.5, 0.75, 1]) { const g = build(w); if (minArea(g) > 0) { grid = g; used = w; break; } }
  if (!grid) throw new Error('The section grid could not be generated without inverted cells. Use fewer wall-normal cells or a larger first-cell height.');
  const nodes = new Float64Array(2 * nNodes);
  for (let i = 0; i < ni; i++) for (let j = 0; j <= nj; j++) { const q = nid[i * (nj + 1) + j]; if (isC && j === 0 && i >= nW + 2 * n) continue; nodes[2 * q] = grid.X[i * (nj + 1) + j] * chord; nodes[2 * q + 1] = grid.Y[i * (nj + 1) + j] * chord; }
  const quads = new Int32Array(4 * nCellsI * nj); let e = 0;
  for (let i = 0; i < nCellsI; i++) for (let j = 0; j < nj; j++) { quads[e++] = id(i, j); quads[e++] = id(i + 1, j); quads[e++] = id(i + 1, j + 1); quads[e++] = id(i, j + 1); }
  const airfoil = [], farfield = [];
  if (isC) {
    for (let i = nW; i < nW + 2 * n; i++) airfoil.push([id(i, 0), id(i + 1, 0)]);
    for (let i = 0; i < ni - 1; i++) farfield.push([id(i, nj), id(i + 1, nj)]);
    for (let j = 0; j < nj; j++) { farfield.push([id(0, j), id(0, j + 1)]); farfield.push([id(ni - 1, j), id(ni - 1, j + 1)]); }
  } else for (let i = 0; i < ni; i++) { airfoil.push([id(i, 0), id(i + 1, 0)]); farfield.push([id(i, nj), id(i + 1, nj)]); }
  return { dim: 2, nodes, quads, markers: { airfoil, farfield }, info: { type: isC ? 'C' : 'O', naca: sec.code, exactSection: sec.exact, nSurf: n, nWake: nW, nNormal: nj, cells: quads.length / 4, points: nNodes, firstCell: sj[1] * R * chord, growth: ratio, radius_chords: R, wake_chords: isC ? wake : 0, arcBlend: used, chord } };
}

/** Extrude a 2-D section grid along the span (y) into hexahedra with a symmetry plane at each end. */
export function extrudeGrid(g, span, nSpan = 8) {
  const n2 = g.nodes.length / 2, ns = Math.max(1, Math.round(nSpan)), nodes = new Float64Array(3 * n2 * (ns + 1));
  for (let k = 0; k <= ns; k++) for (let i = 0; i < n2; i++) { const o = 3 * (k * n2 + i); nodes[o] = g.nodes[2 * i]; nodes[o + 1] = (span * k) / ns; nodes[o + 2] = g.nodes[2 * i + 1]; }
  const nq = g.quads.length / 4, hexes = new Int32Array(8 * nq * ns); let e = 0;
  for (let k = 0; k < ns; k++) for (let q = 0; q < nq; q++) { const a = g.quads[4 * q], b = g.quads[4 * q + 1], c = g.quads[4 * q + 2], d = g.quads[4 * q + 3]; for (const off of [k * n2, (k + 1) * n2]) { hexes[e++] = a + off; hexes[e++] = d + off; hexes[e++] = c + off; hexes[e++] = b + off; } }
  const sweep = (edges) => { const out = []; for (let k = 0; k < ns; k++) for (const [a, b] of edges) out.push([a + k * n2, b + k * n2, b + (k + 1) * n2, a + (k + 1) * n2]); return out; };
  const cap = (k) => { const out = []; for (let q = 0; q < nq; q++) out.push([0, 1, 2, 3].map((m) => g.quads[4 * q + m] + k * n2)); return out; };
  return { dim: 3, nodes, hexes, markers: { airfoil: sweep(g.markers.airfoil), farfield: sweep(g.markers.farfield), sym_root: cap(0), sym_tip: cap(ns) }, info: { ...g.info, span, nSpan: ns, cells: nq * ns, points: n2 * (ns + 1) } };
}

/** Serialise a section grid (2-D quads) or an extruded grid (hexahedra) in SU2 native ASCII format. */
export function su2MeshText(g) {
  const out = [`% AeroSuite 26 high-fidelity bridge — ${g.info.type}-grid, NACA ${g.info.naca}, chord ${f(g.info.chord, 6)} m`, `NDIME= ${g.dim}`];
  if (g.dim === 2) {
    const nq = g.quads.length / 4; out.push(`NELEM= ${nq}`);
    for (let q = 0; q < nq; q++) out.push(`9 ${g.quads[4 * q]} ${g.quads[4 * q + 1]} ${g.quads[4 * q + 2]} ${g.quads[4 * q + 3]} ${q}`);
    const np = g.nodes.length / 2; out.push(`NPOIN= ${np}`);
    for (let i = 0; i < np; i++) out.push(`${f(g.nodes[2 * i], 15)} ${f(g.nodes[2 * i + 1], 15)} ${i}`);
  } else {
    const nh = g.hexes.length / 8; out.push(`NELEM= ${nh}`);
    for (let q = 0; q < nh; q++) { let s = '12'; for (let m = 0; m < 8; m++) s += ' ' + g.hexes[8 * q + m]; out.push(`${s} ${q}`); }
    const np = g.nodes.length / 3; out.push(`NPOIN= ${np}`);
    for (let i = 0; i < np; i++) out.push(`${f(g.nodes[3 * i], 15)} ${f(g.nodes[3 * i + 1], 15)} ${f(g.nodes[3 * i + 2], 15)} ${i}`);
  }
  const names = Object.keys(g.markers); out.push(`NMARK= ${names.length}`);
  for (const nm of names) { const els = g.markers[nm]; out.push(`MARKER_TAG= ${nm}`, `MARKER_ELEMS= ${els.length}`); for (const el of els) out.push(`${el.length === 2 ? 3 : 9} ${el.join(' ')}`); }
  return out.join('\n') + '\n';
}

// ---------------------------------------------------------------------------------------------
// Settings: auto-filled from the case, reviewed and edited by the user
// ---------------------------------------------------------------------------------------------
export const SOLVERS = {
  su2: { label: 'SU2', what: 'Body-fitted RANS / Euler flow solver (compressible and incompressible)', site: 'https://su2code.github.io', analyses: { 'rans-sa': 'RANS, Spalart–Allmaras', 'rans-sst': 'RANS, k-ω SST', euler: 'Euler (inviscid)' } },
  openfoam: { label: 'OpenFOAM', what: 'General finite-volume CFD on a snappyHexMesh mesh of the 3-D shape', site: 'https://www.openfoam.com', analyses: { simpleFoam: 'Steady incompressible RANS (simpleFoam)', pimpleFoam: 'Large-eddy simulation with the WALE model (pimpleFoam)', rhoSimpleFoam: 'Steady compressible RANS (rhoSimpleFoam)' } },
  calculix: { label: 'CalculiX', what: 'Detailed shell / solid finite elements (decks also read by Abaqus)', site: 'https://www.calculix.de', analyses: { all: 'Static + buckling + modal', static: 'Linear static', buckle: 'Linear buckling', modal: 'Natural frequencies', dynamic: 'Dynamic response to a load pulse (implicit)', explicit: 'Dynamic response to a load pulse (explicit)' } },
};

const F = (key, label, unit, o = {}) => ({ key, label, unit, ...o });
export const SETTING_FIELDS = {
  su2: [
    F('regime', 'Flow model', '', { type: 'select', options: ['compressible', 'incompressible'], group: 'Flow', help: 'Incompressible is recommended below Mach 0.3; compressible above' }),
    F('mesh', 'Mesh', '', { type: 'select', options: ['C-grid', 'O-grid', 'imported'], group: 'Mesh', help: 'C-grid: section mesh with a wake cut (viscous). O-grid: closed rings (inviscid). Imported: the volume mesh selected below' }),
    F('dim', 'Dimensions', '', { type: 'select', options: ['2-D section', '3-D extruded'], group: 'Mesh', help: '3-D extruded sweeps the section grid along the span between two symmetry planes' }),
    F('resolution', 'Resolution', '', { type: 'select', options: Object.keys(RESOLUTIONS), group: 'Mesh' }),
    F('naca', 'NACA 4-digit section', '', { type: 'text', group: 'Mesh' }),
    F('chord_m', 'Chord', 'm', { group: 'Mesh', min: 1e-3 }), F('span_m', 'Extruded span', 'm', { group: 'Mesh', min: 1e-3, help: 'Used by the 3-D extruded mesh only' }),
    F('yplus', 'Target y+', '-', { group: 'Mesh', min: 0.05, max: 300, help: '1 resolves the viscous sublayer (needed without wall functions)' }),
    F('farfield_chords', 'Far-field distance', 'chords', { group: 'Mesh', min: 5, max: 500 }),
    F('mach', 'Mach number', '-', { group: 'Flow', min: 0.001, max: 5 }), F('aoa_deg', 'Angle of attack', 'deg', { group: 'Flow', min: -30, max: 30 }),
    F('reynolds', 'Reynolds number (on chord)', '-', { group: 'Flow', min: 1 }), F('T_K', 'Static temperature', 'K', { group: 'Flow', min: 100 }),
    F('p_Pa', 'Static pressure', 'Pa', { group: 'Flow', min: 1 }), F('V_ms', 'Speed', 'm/s', { group: 'Flow', min: 0.01 }),
    F('ref_area_m2', 'Reference area', 'm²', { group: 'Reference', min: 1e-9, help: 'For a 2-D section this is the chord × 1 m' }), F('ref_length_m', 'Reference length', 'm', { group: 'Reference', min: 1e-9 }),
    F('iter', 'Iterations', '', { group: 'Numerics', min: 1, step: 1 }), F('cfl', 'CFL number', '-', { group: 'Numerics', min: 0.1 }),
    F('mg_levels', 'Multigrid levels', '', { group: 'Numerics', min: 0, max: 4, step: 1, help: '0 is the most robust on stretched viscous grids' }),
    F('conv_log10', 'Stop at residual (log10)', '', { group: 'Numerics', min: -16, max: -2 }),
    F('wall_markers', 'Wall markers (imported mesh)', '', { type: 'text', group: 'Imported mesh', help: 'Comma-separated boundary names' }),
    F('far_markers', 'Far-field markers (imported mesh)', '', { type: 'text', group: 'Imported mesh' }), F('sym_markers', 'Symmetry markers (imported mesh)', '', { type: 'text', group: 'Imported mesh' }),
  ],
  openfoam: [
    F('turbulence', 'Turbulence model', '', { type: 'select', options: ['kOmegaSST', 'SpalartAllmaras'], group: 'Flow', help: 'The LES route always uses the WALE sub-grid model' }),
    F('geometry', 'Shape to mesh', '', { type: 'select', options: ['wing', 'fuselage', 'imported'], group: 'Mesh', help: 'Wing and fuselage are generated from the case; imported uses the surface selected below' }),
    F('V_ms', 'Speed', 'm/s', { group: 'Flow', min: 0.01 }), F('aoa_deg', 'Angle of attack', 'deg', { group: 'Flow', min: -30, max: 30 }),
    F('rho', 'Density', 'kg/m³', { group: 'Flow', min: 1e-4 }), F('nu', 'Kinematic viscosity', 'm²/s', { group: 'Flow', min: 1e-9 }),
    F('T_K', 'Static temperature', 'K', { group: 'Flow', min: 100 }), F('p_Pa', 'Static pressure', 'Pa', { group: 'Flow', min: 1 }),
    F('turb_intensity', 'Turbulence intensity', '-', { group: 'Flow', min: 1e-5, max: 0.5 }),
    F('ref_area_m2', 'Reference area', 'm²', { group: 'Reference', min: 1e-9 }), F('ref_length_m', 'Reference length', 'm', { group: 'Reference', min: 1e-9 }),
    F('cells_per_chord', 'Surface cells per reference length', '', { group: 'Mesh', min: 8, max: 512, step: 1, help: 'Sets the finest snappyHexMesh level; cost grows with its square' }),
    F('layers', 'Prism layers on the wall', '', { group: 'Mesh', min: 0, max: 20, step: 1, help: '0 switches layer addition off (fastest, wall functions carry the boundary layer)' }),
    F('iterations', 'Iterations (steady)', '', { group: 'Numerics', min: 1, step: 1 }), F('end_time_s', 'Simulated time (LES)', 's', { group: 'Numerics', min: 1e-6 }),
  ],
  calculix: [
    F('source', 'Model', '', { type: 'select', options: ['wing box', 'imported'], group: 'Model', help: 'Wing box: shell model of skins, spars and ribs generated from the case. Imported: the structural mesh selected below' }),
    F('element', 'Shell element', '', { type: 'select', options: ['S8R', 'S4'], group: 'Model', help: 'S8R (8-node, reduced integration) is the accurate CalculiX shell; S4 is the 4-node alternative' }),
    F('material', 'Material', '', { type: 'select', options: Object.keys(METALS), group: 'Model' }),
    F('t_skin_mm', 'Skin thickness', 'mm', { group: 'Model', min: 0.05 }), F('t_spar_mm', 'Spar web thickness', 'mm', { group: 'Model', min: 0.05 }), F('t_rib_mm', 'Rib thickness', 'mm', { group: 'Model', min: 0.05 }),
    F('load_factor', 'Load factor', 'g', { group: 'Loads', min: -10, max: 15, help: 'The lift on one wing is this multiple of half the weight' }),
    F('lift_1g_N', 'Lift on this member at 1 g', 'N', { group: 'Loads', min: 0 }),
    F('n_span', 'Elements along the span', '', { group: 'Mesh', min: 4, max: 400, step: 1 }), F('n_chord', 'Elements across the box', '', { group: 'Mesh', min: 2, max: 60, step: 1 }),
    F('n_height', 'Elements through the depth', '', { group: 'Mesh', min: 1, max: 30, step: 1 }), F('rib_every', 'Rib every … span elements', '', { group: 'Mesh', min: 1, max: 100, step: 1 }),
    F('n_modes', 'Natural frequencies to find', '', { group: 'Numerics', min: 1, max: 50, step: 1 }), F('n_buckle', 'Buckling modes to find', '', { group: 'Numerics', min: 1, max: 20, step: 1 }),
    F('dyn_time_s', 'Load pulse duration', 's', { group: 'Numerics', min: 1e-6, help: 'Length of the 1 − cos pulse in the dynamic analyses; the response is followed for twice this time' }),
  ],
};

const markerKind = (name) => (/sym/i.test(name) ? 'sym' : /far|free|inlet|outlet|inflow|outflow|open|infinity|outer/i.test(name) ? 'far' : 'wall');
/** Boundary-marker names of an imported model as the SU2 exporter will write them. */
export function modelMarkers(model) {
  if (!model) return [];
  try { return [...exportModel(model, 'su2').matchAll(/^MARKER_TAG= (.+)$/gm)].map((m) => m[1].trim()); } catch { return []; }
}
export const modelHasVolume = (model) => !!model?.elements?.some((e) => /^(tet|hex|wedge|pyr)/.test(e.type));
export const modelHasSurface = (model) => (model?.triangles?.length || 0) > 0;
export const modelHasStructure = (model) => !!model?.elements?.some((e) => /^(tri|quad|tet|hex|wedge)/.test(e.type));

/** Settings for a solver, filled from the case (and from published native results where they help). */
export function defaultSettings(c, solver, { analysis = null, model = null, up = {} } = {}) {
  const fs = flowState(c), d = fs.d;
  if (solver === 'su2') {
    const a = analysis || 'rans-sa', vol = modelHasVolume(model), names = vol ? modelMarkers(model) : [];
    const pick = (k) => names.filter((nm) => markerKind(nm) === k).join(', ');
    return {
      analysis: a, regime: fs.mach < 0.3 ? 'incompressible' : 'compressible', mesh: vol ? 'imported' : a === 'euler' ? 'O-grid' : 'C-grid', dim: '2-D section', resolution: 'medium',
      naca: parseNaca(c.wing.airfoil, c.wing.tc).code, chord_m: Number(fs.chord.toPrecision(6)), span_m: Number((fs.hasWing ? Math.min(fs.b / 2, fs.chord) : fs.chord).toPrecision(4)), yplus: 1, farfield_chords: 40,
      mach: Number(fs.mach.toPrecision(5)), aoa_deg: fs.alpha_deg, reynolds: Number(fs.reynolds.toPrecision(5)), T_K: Number(fs.T.toPrecision(6)), p_Pa: Number(fs.p.toPrecision(6)), V_ms: fs.V,
      ref_area_m2: vol ? Number(fs.S.toPrecision(6)) : Number(fs.chord.toPrecision(6)), ref_length_m: Number((vol ? fs.refLength : fs.chord).toPrecision(6)),
      iter: a === 'euler' ? 2000 : 6000, cfl: a === 'euler' ? 5 : 10, mg_levels: 0, conv_log10: -8, wall_markers: pick('wall'), far_markers: pick('far'), sym_markers: pick('sym'),
    };
  }
  if (solver === 'openfoam') {
    const a = analysis || 'simpleFoam', surf = modelHasSurface(model);
    return {
      analysis: a, turbulence: 'kOmegaSST', geometry: surf ? 'imported' : fs.hasWing ? 'wing' : 'fuselage', V_ms: fs.V, aoa_deg: fs.alpha_deg, rho: Number(fs.rho.toPrecision(6)), nu: Number(fs.nu.toPrecision(5)),
      T_K: Number(fs.T.toPrecision(6)), p_Pa: Number(fs.p.toPrecision(6)), turb_intensity: fs.turb, ref_area_m2: Number(fs.S.toPrecision(6)), ref_length_m: Number(fs.refLength.toPrecision(6)),
      cells_per_chord: 48, layers: 3, iterations: 1500, end_time_s: Number(((20 * fs.refLength) / fs.V).toPrecision(3)),
    };
  }
  if (solver === 'calculix') {
    const s = c.struct, heli = !fs.hasWing, ch = heli ? c.rotor.chord_m || 0.3 : d.c_root;
    return {
      analysis: analysis || 'all', source: modelHasStructure(model) && !modelHasVolumeOnlyFluid(model) ? 'imported' : 'wing box', element: 'S8R', material: METALS[s.material] ? s.material : 'Al 2024-T3',
      t_skin_mm: Math.min(s.t_skin_mm, Math.max(0.4, 0.6 * ch)), t_spar_mm: Math.min(s.t_spar_mm, Math.max(0.5, ch)), t_rib_mm: Math.min(s.t_skin_mm, Math.max(0.4, 0.6 * ch)),
      load_factor: up.performance?.n_limit ?? c.aero.n_pos, lift_1g_N: Number((d.W / 2).toPrecision(6)), n_span: 40, n_chord: 6, n_height: 2, rib_every: 4, n_modes: 6, n_buckle: 4, dyn_time_s: 0.05,
    };
  }
  throw new Error(`Unknown solver "${solver}".`);
}
const modelHasVolumeOnlyFluid = (model) => /cfd|fluid|flow/i.test(model?.kind || '') && modelHasVolume(model);

// ---------------------------------------------------------------------------------------------
// SU2
// ---------------------------------------------------------------------------------------------
const list = (s) => String(s || '').split(',').map((x) => x.trim()).filter(Boolean);
/** SU2 configuration file text for the given settings and boundary-marker roles. */
export function su2Config(s, { meshFile = 'mesh.su2', walls = ['airfoil'], far = ['farfield'], sym = [], nDim = 2, title = '' } = {}) {
  const euler = s.analysis === 'euler', sst = s.analysis === 'rans-sst', inc = s.regime === 'incompressible';
  const solver = inc ? (euler ? 'INC_EULER' : 'INC_RANS') : euler ? 'EULER' : 'RANS';
  const rho = s.p_Pa / (R_AIR * s.T_K), a = Math.sqrt(GAMMA * R_AIR * s.T_K), V = inc ? s.V_ms : s.mach * a, mu = (rho * V * s.ref_length_m) / Math.max(1, s.reynolds), al = s.aoa_deg * RAD;
  const vel = nDim === 2 ? [V * Math.cos(al), V * Math.sin(al), 0] : [V * Math.cos(al), 0, V * Math.sin(al)];
  const L = [], sec = (t) => L.push('', `% ---- ${t} ----`), kv = (k, v, note) => L.push(`${k}= ${v}${note ? `   % ${note}` : ''}`), tup = (arr) => `( ${arr.join(', ')} )`;
  L.push(`% SU2 configuration generated by the AeroSuite 26 high-fidelity bridge`, `% ${title}`, `% Written for SU2 version 8. Edit freely; every value below came from the shared case or the bridge settings.`);
  sec('Problem definition');
  kv('SOLVER', solver); if (!euler) kv('KIND_TURB_MODEL', sst ? 'SST' : 'SA'); kv('MATH_PROBLEM', 'DIRECT'); kv('RESTART_SOL', 'NO');
  if (inc) {
    sec('Incompressible freestream');
    kv('INC_DENSITY_MODEL', 'CONSTANT'); kv('INC_ENERGY_EQUATION', 'NO'); kv('INC_DENSITY_INIT', f(rho, 7), 'kg/m3'); kv('INC_VELOCITY_INIT', tup(vel.map((v) => f(v, 8))), `m/s — ${f(V, 6)} m/s at ${f(s.aoa_deg, 5)} deg`);
    kv('INC_TEMPERATURE_INIT', f(s.T_K, 7), 'K'); kv('INC_NONDIM', 'INITIAL_VALUES');
    if (!euler) { kv('VISCOSITY_MODEL', 'CONSTANT_VISCOSITY'); kv('MU_CONSTANT', f(mu, 7), `Pa s — gives Reynolds number ${f(s.reynolds, 5)} on ${f(s.ref_length_m, 5)} m`); }
  } else {
    sec('Compressible freestream');
    kv('MACH_NUMBER', f(s.mach, 7)); kv('AOA', f(s.aoa_deg, 7), 'deg'); kv('SIDESLIP_ANGLE', '0.0'); kv('FREESTREAM_TEMPERATURE', f(s.T_K, 7), 'K'); kv('FREESTREAM_PRESSURE', f(s.p_Pa, 8), 'Pa');
    if (!euler) { kv('INIT_OPTION', 'REYNOLDS'); kv('FREESTREAM_OPTION', 'TEMPERATURE_FS'); kv('REYNOLDS_NUMBER', f(s.reynolds, 7)); kv('REYNOLDS_LENGTH', f(s.ref_length_m, 7), 'm'); }
    kv('REF_DIMENSIONALIZATION', euler ? 'DIMENSIONAL' : 'FREESTREAM_VEL_EQ_MACH');
  }
  if (!euler) { if (sst) { kv('FREESTREAM_TURBULENCEINTENSITY', '0.001'); kv('FREESTREAM_TURB2LAMVISCRATIO', '10.0'); } else kv('FREESTREAM_NU_FACTOR', '3.0'); }
  sec('Reference values for the force coefficients');
  kv('REF_ORIGIN_MOMENT_X', f(0.25 * s.ref_length_m, 7), 'quarter chord'); kv('REF_ORIGIN_MOMENT_Y', '0.0'); kv('REF_ORIGIN_MOMENT_Z', '0.0'); kv('REF_LENGTH', f(s.ref_length_m, 7), 'm'); kv('REF_AREA', f(s.ref_area_m2, 7), nDim === 2 ? 'm2 per metre of span' : 'm2');
  sec('Boundary conditions');
  if (euler) kv('MARKER_EULER', tup(walls)); else kv('MARKER_HEATFLUX', tup(walls.flatMap((w) => [w, '0.0'])), 'adiabatic no-slip wall');
  kv('MARKER_FAR', tup(far)); if (sym.length) kv('MARKER_SYM', tup(sym)); kv('MARKER_PLOTTING', tup(walls)); kv('MARKER_MONITORING', tup(walls));
  sec('Numerical method');
  kv('NUM_METHOD_GRAD', inc ? 'WEIGHTED_LEAST_SQUARES' : 'GREEN_GAUSS'); kv('CFL_NUMBER', f(s.cfl, 5));
  kv('CFL_ADAPT', 'YES'); kv('CFL_ADAPT_PARAM', tup(['0.8', '1.1', f(Math.min(s.cfl, 1), 4), f(s.cfl * 10, 5)]), 'factor down, factor up, minimum, maximum');
  kv('ITER', String(Math.round(s.iter)));
  kv('LINEAR_SOLVER', 'FGMRES'); kv('LINEAR_SOLVER_PREC', 'ILU'); kv('LINEAR_SOLVER_ERROR', '1E-6'); kv('LINEAR_SOLVER_ITER', '10');
  kv('MGLEVEL', String(Math.round(s.mg_levels))); if (s.mg_levels > 0) { kv('MGCYCLE', 'W_CYCLE'); kv('MG_PRE_SMOOTH', tup(['1', '2', '3', '3'])); kv('MG_POST_SMOOTH', tup(['0', '0', '0', '0'])); kv('MG_CORRECTION_SMOOTH', tup(['0', '0', '0', '0'])); kv('MG_DAMP_RESTRICTION', '0.75'); kv('MG_DAMP_PROLONGATION', '0.75'); }
  if (inc) { kv('CONV_NUM_METHOD_FLOW', 'FDS'); kv('MUSCL_FLOW', 'YES'); kv('SLOPE_LIMITER_FLOW', 'VENKATAKRISHNAN', 'the limiter keeps the start-up from freestream stable'); kv('VENKAT_LIMITER_COEFF', '0.01'); }
  else if (s.mach >= 0.6 || euler) { kv('CONV_NUM_METHOD_FLOW', 'JST'); kv('JST_SENSOR_COEFF', tup(['0.5', '0.02'])); }
  else { kv('CONV_NUM_METHOD_FLOW', 'ROE'); kv('MUSCL_FLOW', 'YES'); kv('SLOPE_LIMITER_FLOW', 'VENKATAKRISHNAN'); kv('VENKAT_LIMITER_COEFF', '0.05'); }
  kv('TIME_DISCRE_FLOW', 'EULER_IMPLICIT');
  if (!euler) { kv('CONV_NUM_METHOD_TURB', 'SCALAR_UPWIND'); kv('MUSCL_TURB', 'NO'); kv('TIME_DISCRE_TURB', 'EULER_IMPLICIT'); kv('CFL_REDUCTION_TURB', '1.0'); }
  sec('Convergence');
  kv('CONV_FIELD', inc ? 'RMS_PRESSURE' : 'RMS_DENSITY'); kv('CONV_RESIDUAL_MINVAL', f(s.conv_log10, 4), 'log10 of the residual'); kv('CONV_STARTITER', '10');
  sec('Files');
  kv('MESH_FILENAME', meshFile); kv('MESH_FORMAT', 'SU2'); kv('TABULAR_FORMAT', 'CSV'); kv('CONV_FILENAME', 'history'); kv('RESTART_FILENAME', 'restart_flow.dat'); kv('SOLUTION_FILENAME', 'solution_flow.dat');
  kv('VOLUME_FILENAME', 'flow'); kv('SURFACE_FILENAME', 'surface_flow'); kv('OUTPUT_FILES', tup(['RESTART', 'PARAVIEW', 'SURFACE_PARAVIEW', 'SURFACE_CSV'])); kv('OUTPUT_WRT_FREQ', String(Math.max(1, Math.round(s.iter))));
  const turbRes = euler ? [] : sst ? ['RMS_TKE', 'RMS_DISSIPATION'] : ['RMS_NU_TILDE'];
  kv('SCREEN_OUTPUT', tup(['INNER_ITER', 'WALL_TIME', inc ? 'RMS_PRESSURE' : 'RMS_DENSITY', ...turbRes, 'LIFT', 'DRAG'])); kv('SCREEN_WRT_FREQ_INNER', '10');
  kv('HISTORY_OUTPUT', tup(['ITER', 'RMS_RES', 'AERO_COEFF'])); kv('HISTORY_WRT_FREQ_INNER', '1');
  return L.join('\n') + '\n';
}

function su2Case(c, s, model) {
  const files = [], notes = [], warnings = [], title = `${c.meta.name || 'case'} — ${SOLVERS.su2.analyses[s.analysis]}, ${s.regime}`;
  let walls, far, sym = [], nDim, meshInfo;
  if (s.mesh === 'imported') {
    if (!modelHasVolume(model) && !(model?.elements || []).some((e) => /^(tri|quad)/.test(e.type))) throw new Error('Select an imported volume mesh (tetrahedra, hexahedra, prisms or pyramids) or a planar 2-D mesh, or switch the mesh to C-grid or O-grid.');
    const text = exportModel(model, 'su2'), names = [...text.matchAll(/^MARKER_TAG= (.+)$/gm)].map((m) => m[1].trim());
    walls = list(s.wall_markers); far = list(s.far_markers); sym = list(s.sym_markers); nDim = Number(/NDIME= (\d)/.exec(text)?.[1] || 3);
    const unknown = [...walls, ...far, ...sym].filter((nm) => !names.includes(nm)), unassigned = names.filter((nm) => ![...walls, ...far, ...sym].includes(nm));
    if (unknown.length) throw new Error(`These boundary names are not in the imported mesh: ${unknown.join(', ')}. The mesh has: ${names.join(', ') || 'no named boundaries'}.`);
    if (!walls.length || !far.length) throw new Error(`Name at least one wall marker and one far-field marker. The imported mesh has: ${names.join(', ') || 'no named boundaries'}.`);
    if (unassigned.length) throw new Error(`Every boundary needs a role. Not yet assigned: ${unassigned.join(', ')}.`);
    files.push({ path: 'mesh.su2', text }); meshInfo = { source: 'imported', name: model.name, markers: names, nDim };
    notes.push(`Mesh: imported volume mesh "${model.name}" exported in SU2 format; boundary roles as listed in the configuration.`);
  } else {
    const res = RESOLUTIONS[s.resolution] || RESOLUTIONS.medium, euler = s.analysis === 'euler';
    const first = euler ? 2e-3 * s.chord_m : wallSpacing(s.reynolds, s.chord_m, s.yplus);
    const g2 = airfoilGrid({ naca: s.naca, tc: c.wing.tc, chord: s.chord_m, type: s.mesh === 'O-grid' ? 'O' : 'C', nSurf: res.nSurf, nWake: res.nWake, growth: euler ? 1.12 : res.growth, firstCell: first, radius: s.farfield_chords, wake: Math.max(10, 0.75 * s.farfield_chords) });
    const g = s.dim === '3-D extruded' ? extrudeGrid(g2, s.span_m, Math.max(4, Math.round(res.nSurf / 12))) : g2;
    walls = ['airfoil']; far = ['farfield']; sym = g.dim === 3 ? ['sym_root', 'sym_tip'] : []; nDim = g.dim;
    files.push({ path: 'mesh.su2', text: su2MeshText(g) }); meshInfo = { source: 'generated', ...g.info };
    notes.push(`Mesh: structured ${g.info.type}-grid around NACA ${g.info.naca}, ${g.info.cells} cells, first node ${g.info.firstCell.toExponential(2)} m off the wall (target y+ ${euler ? 'not applicable' : s.yplus}), growth ratio ${g.info.growth.toFixed(3)}.`);
    if (!g.info.exactSection) warnings.push(`"${c.wing.airfoil}" is not a NACA 4-digit code, so a symmetric NACA ${g.info.naca} of the case thickness ratio was used.`);
    if (g.info.growth > 1.3) warnings.push(`The wall-normal growth ratio is ${g.info.growth.toFixed(2)}; choose a finer resolution for a production run.`);
    if (g.dim === 2) notes.push('This is a 2-D section analysis: its coefficients are per unit span and are section values, not whole-aircraft values.');
  }
  const cfg = { ...s }; if (nDim === 3 && s.mesh !== 'imported') cfg.ref_area_m2 = s.chord_m * s.span_m; else if (s.mesh !== 'imported') cfg.ref_area_m2 = s.chord_m;
  if (s.mesh !== 'imported') cfg.ref_length_m = s.chord_m;
  if (s.regime === 'compressible' && s.mach < 0.2) warnings.push(`Mach ${s.mach} is low for the compressible solver, which converges slowly there. The incompressible model is recommended.`);
  if (s.regime === 'incompressible' && s.mach > 0.35) warnings.push(`Mach ${s.mach} is too high for the incompressible model; compressibility changes lift and drag noticeably above about Mach 0.3.`);
  files.unshift({ path: 'config.cfg', text: su2Config(cfg, { walls, far, sym, nDim, title }) });
  files.push({ path: 'run.sh', exec: true, text: ['#!/usr/bin/env bash', '# Runs the SU2 case in this folder. SU2_CFD must be on the PATH (https://su2code.github.io/download.html).', 'set -euo pipefail', 'cd "$(dirname "$0")"', 'NP="${NP:-2}"', 'if SU2_CFD --help 2>&1 | grep -q -- "--threads"; then', '  SU2_CFD -t "$NP" config.cfg 2>&1 | tee log.su2', 'elif [ "$NP" -gt 1 ] && command -v mpirun >/dev/null 2>&1; then', '  mpirun -n "$NP" SU2_CFD config.cfg 2>&1 | tee log.su2', 'else', '  SU2_CFD config.cfg 2>&1 | tee log.su2', 'fi', ''].join('\n') });
  return { files, notes, warnings, meshInfo, expected: ['history.csv', 'surface_flow.csv', 'flow.vtu', 'surface_flow.vtu', 'restart_flow.dat', 'log.su2'], commands: ['bash run.sh'] };
}

// ---------------------------------------------------------------------------------------------
// Parametric surfaces (ASCII STL) for snappyHexMesh
// ---------------------------------------------------------------------------------------------
/** Wing stations used by the STL, the Gmsh script and the wing-box model. */
function wingPlan(c) {
  const d = derived(c), w = c.wing, half = w.b_m / 2, cr = d.c_root, ct = d.c_tip, sw = Math.tan((w.sweep_deg || 0) * RAD), di = Math.tan((w.dihedral_deg || 0) * RAD);
  const at = (eta) => { const ch = cr + (ct - cr) * eta, y = half * eta; return { y, chord: ch, xle: y * sw + 0.25 * cr - 0.25 * ch, z: y * di, twist: (w.twist_deg || 0) * eta }; };
  return { half, cr, ct, at, sec: parseNaca(w.airfoil, w.tc), mac: d.mac };
}
const stlText = (name, tris) => { const o = [`solid ${name}`]; for (const [a, b, c] of tris) { const ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2], vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2]; let nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx; const l = Math.hypot(nx, ny, nz) || 1; nx /= l; ny /= l; nz /= l; o.push(` facet normal ${f(nx, 6)} ${f(ny, 6)} ${f(nz, 6)}`, '  outer loop', `   vertex ${f(a[0])} ${f(a[1])} ${f(a[2])}`, `   vertex ${f(b[0])} ${f(b[1])} ${f(b[2])}`, `   vertex ${f(c[0])} ${f(c[1])} ${f(c[2])}`, '  endloop', ' endfacet'); } o.push(`endsolid ${name}`); return o.join('\n') + '\n'; };
const signedVolume = (tris) => { let v = 0; for (const [a, b, c] of tris) v += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6; return v; };

/** Closed, outward-oriented triangles of the full-span wing (both halves), tips capped. */
export function wingTriangles(c, { nAirfoil = 40, nSpan = 24 } = {}) {
  const pl = wingPlan(c), loop = nacaLoop(pl.sec, nAirfoil).slice(0, -1), m = loop.length, rings = [];
  for (let k = -nSpan; k <= nSpan; k++) {
    const eta = Math.abs(k) / nSpan, st = pl.at(0.5 * (1 - Math.cos(Math.PI * eta)) * 0.5 + 0.5 * eta), tw = st.twist * RAD, cs = Math.cos(tw), sn = Math.sin(tw);
    rings.push(loop.map(([x, z]) => { const dx = (x - 0.25) * st.chord, dz = z * st.chord; return [st.xle + 0.25 * st.chord + dx * cs + dz * sn, Math.sign(k) * st.y, st.z - dx * sn + dz * cs]; }));
  }
  const tris = [], sub = (p, q) => [p[0] - q[0], p[1] - q[1], p[2] - q[2]], nrm = (t) => { const u = sub(t[1], t[0]), v = sub(t[2], t[0]); return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]; };
  for (let k = 0; k < rings.length - 1; k++) for (let i = 0; i < m; i++) { const a = rings[k][i], b = rings[k][(i + 1) % m], cc = rings[k + 1][(i + 1) % m], dd = rings[k + 1][i]; tris.push([a, b, cc], [a, cc, dd]); }
  // make the side walls point outwards: test the mid-chord upper-surface facet of the first strip, whose normal must point up
  if (nrm(tris[2 * Math.round(1.5 * nAirfoil)])[2] < 0) for (const t of tris) t.reverse();
  const cap = (ring, sy) => {
    for (let i = 0; i < m / 2; i++) {
      const a = ring[i], b = ring[i + 1], cc = ring[m - i - 1], dd = ring[(m - i) % m], add = (t) => { const ny = nrm(t)[1]; if (ny === 0) return; tris.push(ny * sy > 0 ? t : [t[0], t[2], t[1]]); };
      if (i === 0) add([a, b, cc]); else if (i === m / 2 - 1) add([a, b, dd]); else { add([a, b, cc]); add([a, cc, dd]); }
    }
  };
  cap(rings[0], -1); cap(rings[rings.length - 1], 1);
  return tris;
}
/** Closed triangles of a fuselage-like body of revolution (ellipsoidal nose, cylindrical cabin, tapered tail). */
export function fuselageTriangles(c, { nAx = 48, nCirc = 32 } = {}) {
  const Lf = c.fuselage.len_m, r0 = c.fuselage.dia_m / 2, ln = Math.min(0.3 * Lf, 1.6 * c.fuselage.dia_m), lt = Math.min(0.4 * Lf, 3 * c.fuselage.dia_m), tris = [];
  const rad = (x) => (x < ln ? r0 * Math.sqrt(Math.max(0, 1 - ((ln - x) / ln) ** 2)) : x > Lf - lt ? r0 * Math.sqrt(Math.max(0, 1 - ((x - (Lf - lt)) / lt) ** 2)) ** 1.3 : r0);
  const xs = []; for (let i = 0; i <= nAx; i++) xs.push(0.5 * Lf * (1 - Math.cos((Math.PI * i) / nAx)));
  const ring = (i) => { const x = xs[i], r = rad(x), o = []; for (let j = 0; j < nCirc; j++) { const th = (2 * Math.PI * j) / nCirc; o.push([x - 0.3 * Lf, r * Math.cos(th), r * Math.sin(th)]); } return o; };
  let prev = null;
  for (let i = 1; i < nAx; i++) {
    const cur = ring(i);
    if (!prev) { const tip = [xs[0] - 0.3 * Lf, 0, 0]; for (let j = 0; j < nCirc; j++) tris.push([tip, cur[(j + 1) % nCirc], cur[j]]); }
    else for (let j = 0; j < nCirc; j++) { const a = prev[j], b = prev[(j + 1) % nCirc], cc = cur[(j + 1) % nCirc], dd = cur[j]; tris.push([a, b, cc], [a, cc, dd]); }
    prev = cur;
  }
  const tail = [xs[nAx] - 0.3 * Lf, 0, 0]; for (let j = 0; j < nCirc; j++) tris.push([tail, prev[j], prev[(j + 1) % nCirc]]);
  if (signedVolume(tris) < 0) for (const t of tris) t.reverse();
  return tris;
}
const bboxOf = (pts) => { const mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity]; for (const p of pts) for (let k = 0; k < 3; k++) { if (p[k] < mn[k]) mn[k] = p[k]; if (p[k] > mx[k]) mx[k] = p[k]; } return { min: mn, max: mx, size: mx.map((v, k) => v - mn[k]) }; };

// ---------------------------------------------------------------------------------------------
// OpenFOAM (written for the openfoam.com / ESI releases, v2312 and later)
// ---------------------------------------------------------------------------------------------
const foamHead = (cls, object, location) => ['/*--------------------------------*- C++ -*----------------------------------*\\', '| Generated by the AeroSuite 26 high-fidelity bridge for OpenFOAM (openfoam.com) |', '\\*---------------------------------------------------------------------------*/', 'FoamFile', '{', '    version     2.0;', '    format      ascii;', `    class       ${cls};`, location ? `    location    "${location}";` : null, `    object      ${object};`, '}', '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //', ''].filter((x) => x !== null).join('\n');
const vec = (a) => `(${a.map((v) => f(v, 8)).join(' ')})`;

function openfoamCase(c, s, model) {
  const files = [], notes = [], warnings = [], app = s.analysis, les = app === 'pimpleFoam', comp = app === 'rhoSimpleFoam';
  // --- surface
  let stl, bb, featureHint = 0;
  if (s.geometry === 'imported') {
    if (!modelHasSurface(model)) throw new Error('Select an imported surface (STL, OBJ, STEP tessellation…) or switch the shape to wing or fuselage.');
    stl = exportModel(model, 'stl').replace(/^solid .*$/m, 'solid body').replace(/^endsolid .*$/m, 'endsolid body');
    const P = model.positions, pts = []; for (let i = 0; i < P.length; i += 3) pts.push([P[i], P[i + 1], P[i + 2]]); bb = bboxOf(pts);
    notes.push(`Surface: imported model "${model.name}" written as constant/triSurface/body.stl. It must be closed (watertight) and in metres.`);
    if (model.units?.length && model.units.length !== 'm') warnings.push(`The imported surface is in ${model.units.length}. Convert it to metres on the Geometry page first.`);
  } else {
    const wing = s.geometry === 'wing';
    if (wing && !(c.wing.S_m2 > 0 && c.wing.b_m > 0)) throw new Error('This case has no wing. Choose the fuselage or an imported surface.');
    if (!wing && !(c.fuselage.len_m > 0 && c.fuselage.dia_m > 0)) throw new Error('This case has no fuselage dimensions. Choose the wing or an imported surface.');
    const tris = wing ? wingTriangles(c) : fuselageTriangles(c); stl = stlText('body', tris); bb = bboxOf(tris.flat()); featureHint = wing ? derived(c).c_tip : c.fuselage.dia_m;
    notes.push(wing ? `Surface: full-span wing lofted from the case planform (NACA ${parseNaca(c.wing.airfoil, c.wing.tc).code}, taper, quarter-chord sweep, dihedral and twist), closed at the tips.` : 'Surface: body of revolution with the case fuselage length and diameter (rounded nose, tapered tail).');
  }
  files.push({ path: 'constant/triSurface/body.stl', text: stl });
  // --- domain and refinement levels
  // snappyHexMesh only sees a body that the lines between neighbouring background cell centres cross, so the cells around the
  // body must be smaller than its smallest planform dimension. The block is therefore uniform (cubic cells of size d0) in a
  // core around the body and stretched geometrically towards the far boundaries.
  const Lref = s.ref_length_m, ext = Math.max(...bb.size), ctr = bb.min.map((v, k) => 0.5 * (v + bb.max[k])), sorted = bb.size.slice().sort((a, b) => a - b), feature = Math.min(featureHint || sorted[1], sorted[1] || ext, Lref);
  const dSurf = Lref / Math.max(8, s.cells_per_chord), level = Math.max(1, Math.min(8, Math.floor(Math.log2((0.5 * feature) / dSurf)))), d0 = dSurf * 2 ** level, margin = Math.max(Math.min(Lref, 2 * feature), 3 * d0);
  const far = [[3 * ext, 6 * ext], [2.5 * ext, 2.5 * ext], [2.5 * ext, 2.5 * ext]], dom = { min: [0, 0, 0], max: [0, 0, 0] }, core = { min: [0, 0, 0], max: [0, 0, 0] }, nCells = [0, 0, 0], grading = [];
  for (let k = 0; k < 3; k++) {
    const nc = Math.max(4, Math.ceil((bb.size[k] + 2 * margin) / d0)), mid = 0.5 * (bb.min[k] + bb.max[k]); core.min[k] = mid - 0.5 * nc * d0; core.max[k] = mid + 0.5 * nc * d0;
    const side = far[k].map((dist) => { const len = Math.max(4 * d0, dist - margin), n = cellsForGrowth(d0 / len, 1.3), r = geometricRatio(d0 / len, n); return { len, n, expand: r ** (n - 1) }; });
    dom.min[k] = core.min[k] - side[0].len; dom.max[k] = core.max[k] + side[1].len; nCells[k] = side[0].n + nc + side[1].n;
    const Ltot = dom.max[k] - dom.min[k];
    grading.push(`((${f(side[0].len / Ltot, 9)} ${f(side[0].n / nCells[k], 9)} ${f(1 / side[0].expand, 7)}) (${f((nc * d0) / Ltot, 9)} ${f(nc / nCells[k], 9)} 1) (${f(side[1].len / Ltot, 9)} ${f(side[1].n / nCells[k], 9)} ${f(side[1].expand, 7)}))`);
  }
  const nBackground = nCells[0] * nCells[1] * nCells[2];
  if (nBackground > 600000) warnings.push(`The background mesh alone has ${nBackground} cells because the body has a small feature (${feature.toPrecision(3)} m) relative to its size. Lower "surface cells per reference length" for a first run.`);
  const al = s.aoa_deg * RAD, U = [s.V_ms * Math.cos(al), 0, s.V_ms * Math.sin(al)], I = s.turb_intensity, k0 = 1.5 * (I * s.V_ms) ** 2, ratio = 10, om0 = k0 / (s.nu * ratio), nut0 = les ? 0 : s.nu * ratio, nuT0 = 3 * s.nu;
  const sa = !les && s.turbulence === 'SpalartAllmaras', dt = Number(((0.5 * dSurf) / s.V_ms).toPrecision(3)), steps = Math.max(1, Math.round(s.iterations)), mach = s.V_ms / Math.sqrt(GAMMA * R_AIR * s.T_K);
  if (!comp && mach > 0.35) warnings.push(`The speed is Mach ${mach.toFixed(2)}; use the compressible route (rhoSimpleFoam) above about Mach 0.3.`);
  if (comp && mach >= 0.7) warnings.push(`At Mach ${mach.toFixed(2)} the flow is transonic. The pressure-based rhoSimpleFoam start-up from freestream is fragile there, especially on a coarse mesh: if it stops with a floating-point error, start from a lower speed or first-order schemes, or use the SU2 compressible route, which is built for this regime.`);
  if (les) warnings.push('A wall-resolved or wall-modelled large-eddy simulation needs a far finer mesh and many flow-through times. Treat this deck as a correct starting point and refine it on a cluster.');
  // --- system/controlDict
  const fo = [
    '    forces1', '    {', '        type            forces;', '        libs            (forces);', '        writeControl    timeStep;', `        writeInterval   ${les ? 10 : 1};`, '        patches         ("body.*");',
    comp ? '        rho             rho;' : `        rho             rhoInf;\n        rhoInf          ${f(s.rho, 7)};`, `        CofR            ${vec([bb.min[0] + 0.25 * Lref, ctr[1], ctr[2]])};`, '    }',
    '    forceCoeffs1', '    {', '        type            forceCoeffs;', '        libs            (forces);', '        writeControl    timeStep;', `        writeInterval   ${les ? 10 : 1};`, '        patches         ("body.*");',
    comp ? '        rho             rho;' : '        rho             rhoInf;', `        rhoInf          ${f(s.rho, 7)};`, `        liftDir         ${vec([-Math.sin(al), 0, Math.cos(al)])};`, `        dragDir         ${vec([Math.cos(al), 0, Math.sin(al)])};`,
    `        CofR            ${vec([bb.min[0] + 0.25 * Lref, ctr[1], ctr[2]])};`, '        pitchAxis       (0 1 0);', `        magUInf         ${f(s.V_ms, 7)};`, `        lRef            ${f(Lref, 7)};`, `        Aref            ${f(s.ref_area_m2, 7)};`, '    }',
    '    solverInfo1', '    {', '        type            solverInfo;', '        libs            (utilityFunctionObjects);', `        fields          (U p${comp ? ' e' : ''}${les ? '' : sa ? ' nuTilda' : ' k omega'});`, '        writeResidualFields no;', '    }',
    '    yPlus1', '    {', '        type            yPlus;', '        libs            (fieldFunctionObjects);', '        writeControl    writeTime;', '    }',
  ];
  files.push({ path: 'system/controlDict', text: foamHead('dictionary', 'controlDict', 'system') + [`application     ${app};`, 'startFrom       latestTime;', 'startTime       0;', 'stopAt          endTime;', `endTime         ${les ? f(s.end_time_s, 6) : steps};`, `deltaT          ${les ? dt : 1};`, les ? 'writeControl    adjustable;' : 'writeControl    timeStep;', `writeInterval   ${les ? f(s.end_time_s / 10, 4) : Math.max(1, Math.round(steps / 3))};`, 'purgeWrite      2;', 'writeFormat     binary;', 'writePrecision  8;', 'writeCompression off;', 'timeFormat      general;', 'timePrecision   8;', 'runTimeModifiable true;', les ? 'adjustTimeStep  yes;\nmaxCo           0.8;' : null, '', 'functions', '{', ...fo, '}', ''].filter((x) => x !== null).join('\n') });
  // --- schemes and solution control
  const turbDiv = ['    div(phi,k)      bounded Gauss upwind;', '    div(phi,omega)  bounded Gauss upwind;', '    div(phi,nuTilda) bounded Gauss upwind;'];
  const schemes = les
    ? ['ddtSchemes', '{', '    default         backward;', '}', 'gradSchemes', '{', '    default         Gauss linear;', '}', 'divSchemes', '{', '    default         none;', '    div(phi,U)      Gauss LUST grad(U);', '    div(phi,k)      Gauss limitedLinear 1;', '    div(phi,omega)  Gauss limitedLinear 1;', '    div(phi,nuTilda) Gauss limitedLinear 1;', '    div((nuEff*dev2(T(grad(U))))) Gauss linear;', '}']
    : comp
      ? ['ddtSchemes', '{', '    default         steadyState;', '}', 'gradSchemes', '{', '    default         Gauss linear;', '    grad(U)         cellLimited Gauss linear 1;', '}', 'divSchemes', '{', '    default         bounded Gauss upwind;', '    div(phi,U)      bounded Gauss linearUpwind grad(U);', '    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;', '}']
      : ['ddtSchemes', '{', '    default         steadyState;', '}', 'gradSchemes', '{', '    default         Gauss linear;', '    grad(U)         cellLimited Gauss linear 1;', '}', 'divSchemes', '{', '    default         none;', '    div(phi,U)      bounded Gauss linearUpwind grad(U);', ...turbDiv, '    div((nuEff*dev2(T(grad(U))))) Gauss linear;', '}'];
  files.push({ path: 'system/fvSchemes', text: foamHead('dictionary', 'fvSchemes', 'system') + [...schemes, 'laplacianSchemes', '{', '    default         Gauss linear limited 0.5;', '}', 'interpolationSchemes', '{', '    default         linear;', '}', 'snGradSchemes', '{', '    default         limited 0.5;', '}', 'wallDist', '{', '    method          meshWave;', '}', ''].join('\n') });
  const gamg = (tol, rel) => ['        solver          GAMG;', '        smoother        GaussSeidel;', `        tolerance       ${tol};`, `        relTol          ${rel};`];
  const smooth = (tol, rel) => ['        solver          smoothSolver;', '        smoother        symGaussSeidel;', `        tolerance       ${tol};`, `        relTol          ${rel};`, '        nSweeps         2;'];
  const sol = ['solvers', '{', '    p', '    {', ...gamg('1e-7', les ? '0.05' : '0.05'), '    }'];
  if (les) sol.push('    pFinal', '    {', '        $p;', '        relTol          0;', '    }', '    "(U|k|omega|nuTilda)"', '    {', ...smooth('1e-7', '0.1'), '    }', '    "(U|k|omega|nuTilda)Final"', '    {', '        $U;', '        relTol          0;', '    }', '}', 'PIMPLE', '{', '    nOuterCorrectors 1;', '    nCorrectors     2;', '    nNonOrthogonalCorrectors 1;', '}');
  else {
    sol.push(`    "(U|k|omega|nuTilda${comp ? '|e' : ''})"`, '    {', ...smooth('1e-8', '0.1'), '    }', '}', 'SIMPLE', '{', '    nNonOrthogonalCorrectors 1;', `    consistent      ${comp ? 'no' : 'yes'};`);
    if (comp) sol.push('    transonic       no;', '    pMinFactor      0.2;', '    pMaxFactor      3;');
    sol.push('    residualControl', '    {', '        p               1e-5;', '        U               1e-6;', `        "(k|omega|nuTilda${comp ? '|e' : ''})" 1e-6;`, '    }', '}', 'relaxationFactors', '{', '    fields', '    {', `        p               ${comp ? 0.3 : 1};`, comp ? '        rho             0.05;' : null, '    }', '    equations', '    {', `        U               ${comp ? 0.5 : 0.9};`, comp ? '        e               0.5;' : null, `        "(k|omega|nuTilda)" ${comp ? 0.5 : 0.7};`, '    }', '}');
  }
  if (!les && !comp) sol.push('potentialFlow', '{', '    nNonOrthogonalCorrectors 5;', '}');
  files.push({ path: 'system/fvSolution', text: foamHead('dictionary', 'fvSolution', 'system') + sol.filter((x) => x !== null).join('\n') + '\n' });
  // --- mesh dictionaries
  const [x0, y0, z0] = dom.min, [x1, y1, z1] = dom.max;
  files.push({ path: 'system/blockMeshDict', text: foamHead('dictionary', 'blockMeshDict', 'system') + ['scale 1;', '', 'vertices', '(', ...[[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0], [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]].map((p) => `    ${vec(p)}`), ');', '', 'blocks', '(', `    hex (0 1 2 3 4 5 6 7) (${nCells.join(' ')}) simpleGrading`, '    (', ...grading.map((g) => `        ${g}`), '    )', ');', '', 'edges', '(', ');', '', 'boundary', '(', '    freestream', '    {', '        type patch;', '        faces', '        (', '            (0 4 7 3)', '            (1 2 6 5)', '            (0 1 5 4)', '            (3 7 6 2)', '            (0 3 2 1)', '            (4 5 6 7)', '        );', '    }', ');', ''].join('\n') });
  const wakeBox = { min: [bb.min[0] - 0.5 * Lref, bb.min[1] - 0.5 * Lref, bb.min[2] - 1 * Lref], max: [bb.max[0] + 4 * Lref, bb.max[1] + 0.5 * Lref, bb.max[2] + 1 * Lref] }, nLay = Math.max(0, Math.round(s.layers));
  const inMesh = [core.min[0] + 0.37 * d0, core.min[1] + 0.41 * d0, core.min[2] + 0.43 * d0];               // inside the first core cell: in the fluid, off every cell face
  files.push({ path: 'system/snappyHexMeshDict', text: foamHead('dictionary', 'snappyHexMeshDict', 'system') + [
    'castellatedMesh true;', 'snap            true;', `addLayers       ${nLay > 0 ? 'true' : 'false'};`, '', 'geometry', '{', '    body.stl', '    {', '        type triSurfaceMesh;', '        name body;', '    }', '    wake', '    {', '        type box;', `        min  ${vec(wakeBox.min)};`, `        max  ${vec(wakeBox.max)};`, '    }', '}', '',
    'castellatedMeshControls', '{', '    maxLocalCells   4000000;', '    maxGlobalCells  12000000;', '    minRefinementCells 10;', '    maxLoadUnbalance 0.10;', '    nCellsBetweenLevels 3;', '    features        ();', '    refinementSurfaces', '    {', '        body', '        {', `            level (${Math.max(1, level - 1)} ${level});`, '            patchInfo', '            {', '                type wall;', '            }', '        }', '    }', '    resolveFeatureAngle 30;', '    refinementRegions', '    {', '        wake', '        {', '            mode inside;', `            levels ((1e15 ${Math.max(1, level - 3)}));`, '        }', '    }', `    locationInMesh ${vec(inMesh)};`, '    allowFreeStandingZoneFaces true;', '}', '',
    'snapControls', '{', '    nSmoothPatch    3;', '    tolerance       2.0;', '    nSolveIter      50;', '    nRelaxIter      5;', '    nFeatureSnapIter 10;', '    implicitFeatureSnap true;', '    explicitFeatureSnap false;', '    multiRegionFeatureSnap false;', '}', '',
    'addLayersControls', '{', '    relativeSizes   true;', '    layers', '    {', '        "body.*"', '        {', `            nSurfaceLayers ${Math.max(1, nLay)};`, '        }', '    }', '    expansionRatio  1.2;', '    finalLayerThickness 0.5;', '    minThickness    0.1;', '    nGrow           0;', '    featureAngle    130;', '    slipFeatureAngle 30;', '    nRelaxIter      5;', '    nSmoothSurfaceNormals 1;', '    nSmoothNormals  3;', '    nSmoothThickness 10;', '    maxFaceThicknessRatio 0.5;', '    maxThicknessToMedialRatio 0.3;', '    minMedialAxisAngle 90;', '    nBufferCellsNoExtrude 0;', '    nLayerIter      50;', '    nRelaxedIter    20;', '}', '',
    'meshQualityControls', '{', '    maxNonOrtho     65;', '    maxBoundarySkewness 20;', '    maxInternalSkewness 4;', '    maxConcave      80;', '    minVol          1e-13;', '    minTetQuality   1e-15;', '    minArea         -1;', '    minTwist        0.02;', '    minDeterminant  0.001;', '    minFaceWeight   0.05;', '    minVolRatio     0.01;', '    minTriangleTwist -1;', '    nSmoothScale    4;', '    errorReduction  0.75;', '    relaxed', '    {', '        maxNonOrtho     75;', '    }', '}', '', 'writeFlags      ();', 'mergeTolerance  1e-6;', ''].join('\n') });
  files.push({ path: 'system/decomposeParDict', text: foamHead('dictionary', 'decomposeParDict', 'system') + ['numberOfSubdomains 4;', 'method          scotch;', ''].join('\n') });
  // --- constant
  if (comp) files.push({ path: 'constant/thermophysicalProperties', text: foamHead('dictionary', 'thermophysicalProperties', 'constant') + ['thermoType', '{', '    type            hePsiThermo;', '    mixture         pureMixture;', '    transport       sutherland;', '    thermo          hConst;', '    equationOfState perfectGas;', '    specie          specie;', '    energy          sensibleInternalEnergy;', '}', 'mixture', '{', '    specie', '    {', '        molWeight       28.9647;', '    }', '    thermodynamics', '    {', '        Cp              1004.7;', '        Hf              0;', '    }', '    transport', '    {', '        As              1.458e-06;', '        Ts              110.4;', '    }', '}', ''].join('\n') });
  if (comp) files.push({ path: 'system/fvOptions', text: foamHead('dictionary', 'fvOptions', 'system') + ['// keeps the first iterations from freestream bounded; inactive once the flow has settled', 'limitT', '{', '    type            limitTemperature;', '    active          yes;', '    selectionMode   all;', `    min             ${f(0.5 * s.T_K, 5)};`, `    max             ${f(2.5 * s.T_K, 5)};`, '}', ''].join('\n') });
  files.push({ path: 'constant/transportProperties', text: foamHead('dictionary', 'transportProperties', 'constant') + ['transportModel  Newtonian;', `nu              ${f(s.nu, 7)};`, ''].join('\n') });
  files.push({ path: 'constant/turbulenceProperties', text: foamHead('dictionary', 'turbulenceProperties', 'constant') + (les
    ? ['simulationType LES;', 'LES', '{', '    LESModel        WALE;', '    turbulence      on;', '    printCoeffs     on;', '    delta           cubeRootVol;', '    cubeRootVolCoeffs', '    {', '        deltaCoeff      1;', '    }', '}', '']
    : ['simulationType RAS;', 'RAS', '{', `    RASModel        ${sa ? 'SpalartAllmaras' : 'kOmegaSST'};`, '    turbulence      on;', '    printCoeffs     on;', '}', '']).join('\n') });
  // --- initial and boundary fields
  const field = (name, cls, dims, internal, freestream, wall) => files.push({ path: `0/${name}`, text: foamHead(cls, name, '0') + [`dimensions      ${dims};`, `internalField   uniform ${internal};`, '', 'boundaryField', '{', '    freestream', '    {', ...freestream.map((l) => `        ${l}`), '    }', '    "body.*"', '    {', ...wall.map((l) => `        ${l}`), '    }', '}', ''].join('\n') });
  field('U', 'volVectorField', '[0 1 -1 0 0 0 0]', vec(U), ['type            freestreamVelocity;', 'freestreamValue $internalField;'], ['type            noSlip;']);
  field('p', 'volScalarField', comp ? '[1 -1 -2 0 0 0 0]' : '[0 2 -2 0 0 0 0]', comp ? f(s.p_Pa, 8) : '0', ['type            freestreamPressure;', 'freestreamValue $internalField;'], ['type            zeroGradient;']);
  field('nut', 'volScalarField', '[0 2 -1 0 0 0 0]', f(nut0, 6), ['type            calculated;', 'value           $internalField;'], ['type            nutUSpaldingWallFunction;', 'value           uniform 0;']);
  field('nuTilda', 'volScalarField', '[0 2 -1 0 0 0 0]', f(nuT0, 6), ['type            inletOutlet;', 'inletValue      $internalField;', 'value           $internalField;'], ['type            fixedValue;', 'value           uniform 0;']);
  field('k', 'volScalarField', '[0 2 -2 0 0 0 0]', f(k0, 6), ['type            inletOutlet;', 'inletValue      $internalField;', 'value           $internalField;'], ['type            kqRWallFunction;', 'value           $internalField;']);
  field('omega', 'volScalarField', '[0 0 -1 0 0 0 0]', f(om0, 6), ['type            inletOutlet;', 'inletValue      $internalField;', 'value           $internalField;'], ['type            omegaWallFunction;', 'value           $internalField;']);
  if (comp) {
    field('T', 'volScalarField', '[0 0 0 1 0 0 0]', f(s.T_K, 7), ['type            inletOutlet;', 'inletValue      $internalField;', 'value           $internalField;'], ['type            zeroGradient;']);
    field('alphat', 'volScalarField', '[1 -1 -1 0 0 0 0]', '0', ['type            calculated;', 'value           $internalField;'], ['type            compressible::alphatWallFunction;', 'value           uniform 0;']);
  }
  files.push({ path: 'case.foam', text: '' });
  files.push({ path: 'run.sh', exec: true, text: ['#!/usr/bin/env bash', '# Meshes and runs the OpenFOAM case in this folder. Source the OpenFOAM environment first, for example:', '#   source /usr/lib/openfoam/openfoam2412/etc/bashrc', 'set -eo pipefail', 'cd "$(dirname "$0")"', 'NP="${NP:-1}"', 'blockMesh 2>&1 | tee log.blockMesh', 'snappyHexMesh -overwrite 2>&1 | tee log.snappyHexMesh', 'checkMesh 2>&1 | tee log.checkMesh || true', 'if [ "$NP" -gt 1 ]; then', '  sed -i "s/^numberOfSubdomains.*/numberOfSubdomains $NP;/" system/decomposeParDict', '  decomposePar -force 2>&1 | tee log.decomposePar', `  mpirun --oversubscribe -np "$NP" ${app} -parallel 2>&1 | tee log.${app}`, '  reconstructPar -latestTime 2>&1 | tee log.reconstructPar', 'else', `  ${app} 2>&1 | tee log.${app}`, 'fi', ''].join('\n') });
  const meshInfo = { source: s.geometry, background: nCells, backgroundCells: nBackground, coreCell_m: d0, surfaceLevel: level, surfaceCell_m: d0 / 2 ** level, domain: dom, core, layers: nLay };
  notes.push(`Mesh: background block of ${nCells.join(' × ')} cells, cubic (${d0.toPrecision(3)} m) around the body and stretched towards the far boundaries; snappyHexMesh refines ${level} times to ${(d0 / 2 ** level).toPrecision(3)} m on the surface${nLay ? ` and adds ${nLay} prism layers` : ''}. Wall functions (Spalding) carry the boundary layer.`);
  notes.push(`Force coefficients use reference area ${f(s.ref_area_m2, 5)} m² and length ${f(Lref, 5)} m; lift and drag are resolved along and across the freestream at ${f(s.aoa_deg, 4)}°.`);
  return { files, notes, warnings, meshInfo, expected: ['postProcessing/forceCoeffs1/0/coefficient.dat', 'postProcessing/forces1/0/force.dat', 'postProcessing/solverInfo1/0/solverInfo.dat', `log.${app}`, 'log.snappyHexMesh', 'log.checkMesh'], commands: ['source /usr/lib/openfoam/openfoam2412/etc/bashrc', 'bash run.sh'] };
}

// ---------------------------------------------------------------------------------------------
// CalculiX / Abaqus: shell wing box (skins, spars, ribs) and imported structural meshes
// ---------------------------------------------------------------------------------------------
/** Running lift [N/m] at span fraction eta on a half-wing carrying halfLift in total (Schrenk: mean of elliptic and planform loading). */
const schrenk = (eta, chordRatio, halfLift, half) => (halfLift / half) * 0.5 * (chordRatio + (4 / Math.PI) * Math.sqrt(Math.max(0, 1 - eta * eta)));

/**
 * Shell mesh of the wing box of one half-wing: upper and lower skins, front and rear spar webs, and ribs.
 * The box is the rectangular cell used by the native structures suite (width = box_chord_frac × chord, depth =
 * box_height_frac × thickness ratio × chord) with the front spar at 15 % chord.
 */
export function wingBoxMesh(c, s) {
  const fs = flowState(c), d = fs.d, heli = !fs.hasWing;
  const half = heli ? Math.max(0.5, c.rotor.R_m * (1 - c.rotor.hinge_offset) || 1) : c.wing.b_m / 2, cr = heli ? c.rotor.chord_m || 0.3 : d.c_root, ct = heli ? cr : d.c_tip;
  const tc = heli ? 0.12 : c.wing.tc, bf = heli ? 0.35 : c.struct.box_chord_frac, hf = c.struct.box_height_frac, sw = heli ? 0 : Math.tan((c.wing.sweep_deg || 0) * RAD), xfs = 0.15;
  const ns = Math.max(4, Math.round(s.n_span)), nc = Math.max(2, Math.round(s.n_chord)), nh = Math.max(1, Math.round(s.n_height)), ribEvery = Math.max(1, Math.round(s.rib_every)), quad = s.element === 'S8R';
  const nodes = [], key = new Map(), isRib = (k) => k > 0 && (k % ribEvery === 0 || k === ns);
  const node = (k, a, b) => {
    const id = `${k},${a},${b}`; let q = key.get(id); if (q !== undefined) return q;
    const eta = k / ns, y = half * eta, ch = cr + (ct - cr) * eta, xle = y * sw + 0.25 * cr - 0.25 * ch, w = bf * ch, h = hf * tc * ch;
    q = nodes.length; nodes.push([xle + xfs * ch + (w * a) / nc, y, h * (b / nh - 0.5)]); key.set(id, q); return q;
  };
  const sets = { SKIN_U: [], SKIN_L: [], SPAR_F: [], SPAR_R: [], RIBS: [] }, elems = [];
  const add = (set, q) => { sets[set].push(elems.length); elems.push(q); };
  for (let k = 0; k < ns; k++) {
    for (let a = 0; a < nc; a++) { add('SKIN_L', [node(k, a, 0), node(k, a + 1, 0), node(k + 1, a + 1, 0), node(k + 1, a, 0)]); add('SKIN_U', [node(k, a, nh), node(k + 1, a, nh), node(k + 1, a + 1, nh), node(k, a + 1, nh)]); }
    for (let b = 0; b < nh; b++) { add('SPAR_F', [node(k, 0, b), node(k + 1, 0, b), node(k + 1, 0, b + 1), node(k, 0, b + 1)]); add('SPAR_R', [node(k, nc, b), node(k, nc, b + 1), node(k + 1, nc, b + 1), node(k + 1, nc, b)]); }
  }
  for (let k = 1; k <= ns; k++) if (isRib(k)) for (let a = 0; a < nc; a++) for (let b = 0; b < nh; b++) add('RIBS', [node(k, a, b), node(k, a + 1, b), node(k, a + 1, b + 1), node(k, a, b + 1)]);
  const nCorner = nodes.length, corner = (k, a, b) => key.get(`${k},${a},${b}`);
  if (quad) { const mid = new Map(); for (const q of elems) { for (let m = 0; m < 4; m++) { const a = q[m], b = q[(m + 1) % 4], id = a < b ? `${a}_${b}` : `${b}_${a}`; let mm = mid.get(id); if (mm === undefined) { mm = nodes.length; nodes.push([0, 1, 2].map((x) => 0.5 * (nodes[a][x] + nodes[b][x]))); mid.set(id, mm); } q.push(mm); } } }
  const root = [], tip = []; for (let i = 0; i < nodes.length; i++) { if (Math.abs(nodes[i][1]) < 1e-9 * half) root.push(i); if (Math.abs(nodes[i][1] - half) < 1e-9 * half) tip.push(i); }
  // loads: Schrenk lift × load factor, lumped to the four spar-cap nodes of each station so the resultant acts at 25 % chord
  const Lhalf = s.load_factor * s.lift_1g_N, cbar = 0.5 * (cr + ct), run = (eta) => schrenk(eta, (cr + (ct - cr) * eta) / cbar, Lhalf, half), loads = [];
  const rear = Math.min(1, Math.max(0, (0.25 - xfs) / bf)); let total = 0; const stations = [];
  for (let k = 1; k <= ns; k++) {
    const lo = (k - 0.5) / ns, hi = Math.min(1, (k + 0.5) / ns), m = 8; let F = 0;
    for (let i = 0; i < m; i++) { const e0 = lo + ((hi - lo) * i) / m, e1 = lo + ((hi - lo) * (i + 1)) / m; F += 0.5 * (run(e0) + run(e1)) * (e1 - e0) * half; }
    stations.push(F); total += F;
  }
  const first = (() => { let F = 0; const m = 8, hi = 0.5 / ns; for (let i = 0; i < m; i++) { const e0 = (hi * i) / m, e1 = (hi * (i + 1)) / m; F += 0.5 * (run(e0) + run(e1)) * (e1 - e0) * half; } return F; })();
  const scale = total + first > 0 ? Lhalf / (total + first) : 0;                      // the root half-cell load goes straight into the clamp
  for (let k = 1; k <= ns; k++) { const F = stations[k - 1] * scale; loads.push([corner(k, 0, 0), 0.5 * F * (1 - rear)], [corner(k, 0, nh), 0.5 * F * (1 - rear)], [corner(k, nc, 0), 0.5 * F * rear], [corner(k, nc, nh), 0.5 * F * rear]); }
  const tipRef = corner(ns, 0, nh);
  return { nodes, elems, sets, nCorner, root, tip, tipRef, loads, type: quad ? 'S8R' : 'S4', info: { half, c_root: cr, c_tip: ct, box_w_root: bf * cr, box_h_root: hf * tc * cr, nSpan: ns, nChord: nc, nHeight: nh, ribs: Array.from({ length: ns }, (_, i) => i + 1).filter(isRib).length, nodes: nodes.length, elements: elems.length, appliedLoad_N: loads.reduce((a, l) => a + l[1], 0), halfLift_N: Lhalf, rootReaction_N: first * scale } };
}

const rows = (ids, per = 12) => { const out = []; for (let i = 0; i < ids.length; i += per) out.push(ids.slice(i, i + per).join(', ')); return out; };
function inpMaterial(name) {
  const m = METALS[name] || METALS['Al 2024-T3'];
  return ['*MATERIAL, NAME=MAT', '*ELASTIC', `${f(m.E, 8)}, ${f(m.nu, 6)}`, '*DENSITY', `${f(m.rho, 8)}`];
}
/** The analysis steps as separate decks that include the model file. */
function inpSteps(s, which, { model = 'model.inp', loads, tipSet = 'TIP', rootSet = 'ROOT', refSet = 'TIPREF', stableDt = null }) {
  const cload = ['*CLOAD', ...loads.map(([n, F, dof = 3]) => `${n + 1}, ${dof}, ${f(F, 8)}`)], head = (t) => ['*HEADING', `AeroSuite 26 high-fidelity bridge — ${t}`, `*INCLUDE, INPUT=${model}`];
  const out = {};
  if (which.includes('static')) out['static.inp'] = [...head('linear static, lift distribution × load factor'), '*STEP', '*STATIC', ...cload, '*NODE FILE', 'U', '*EL FILE', 'S', `*NODE PRINT, NSET=${tipSet}`, 'U', `*NODE PRINT, NSET=${rootSet}, TOTALS=ONLY`, 'RF', '*END STEP', ''].join('\n');
  if (which.includes('buckle')) out['buckle.inp'] = [...head('linear buckling under the same load; the factors multiply the applied load'), '*STEP', '*BUCKLE', `${Math.round(s.n_buckle)}`, ...cload, '*NODE FILE', 'U', '*END STEP', ''].join('\n');
  if (which.includes('modal')) out['modal.inp'] = [...head('natural frequencies and mode shapes, root clamped, structural mass only'), '*STEP', '*FREQUENCY', `${Math.round(s.n_modes)}`, '*NODE FILE', 'U', '*END STEP', ''].join('\n');
  if (which.includes('dynamic') || which.includes('explicit')) {
    const T = s.dyn_time_s, amp = []; for (let i = 0; i <= 20; i++) { const t = (T * i) / 20; amp.push(`${f(t, 6)}, ${f(0.5 * (1 - Math.cos((2 * Math.PI * t) / T)), 6)}`); }
    amp.push(`${f(2 * T, 6)}, 0`);
    const gust = ['*AMPLITUDE, NAME=GUST', ...amp], load = [cload[0] + ', AMPLITUDE=GUST', ...cload.slice(1)];
    // implicit (Hilber–Hughes–Taylor) transient: 100 fixed increments per pulse, practical for gust and landing-type loads
    if (which.includes('dynamic')) out['dynamic.inp'] = [...head('implicit dynamic response to a 1 − cos pulse of the same peak load'), ...gust, '*STEP, INC=1000000', '*DYNAMIC, DIRECT', `${f(T / 100, 6)}, ${f(2 * T, 6)}`, ...load, `*NODE PRINT, NSET=${refSet}, FREQUENCY=1`, 'U', '*END STEP', ''].join('\n');
    // explicit central differences: the solver sets the stable increment itself (the value given is only the first guess)
    if (which.includes('explicit')) out['explicit.inp'] = [...head('explicit dynamic response to a 1 − cos pulse of the same peak load'), ...gust, '*STEP, NLGEOM, INC=1000000000', '*DYNAMIC, EXPLICIT', `${f(stableDt || T / 1e5, 4)}, ${f(2 * T, 6)}`, ...load, `*NODE PRINT, NSET=${refSet}, FREQUENCY=2000`, 'U', '*END STEP', ''].join('\n');
  }
  return out;
}
const ABAQUS_TYPE = { tri3: 'S3', quad4: 'S4', tri6: 'S6', quad8: 'S8R', tet4: 'C3D4', tet10: 'C3D10', hex8: 'C3D8', hex20: 'C3D20R', wedge6: 'C3D6', wedge15: 'C3D15' };

function calculixCase(c, s, model) {
  const files = [], notes = [], warnings = [], which = s.analysis === 'all' ? ['static', 'buckle', 'modal'] : [s.analysis], mat = METALS[s.material] || METALS['Al 2024-T3'];
  let stepOpts, meshInfo;
  if (s.source === 'imported') {
    const blocks = (model?.elements || []).filter((b) => ABAQUS_TYPE[b.type]);
    if (!blocks.length) throw new Error('Select an imported structural mesh with shell (triangle, quadrilateral) or solid (tetrahedron, hexahedron, wedge) elements, or switch the model to wing box.');
    const solid = blocks.some((b) => b.type.startsWith('tet') || b.type.startsWith('hex') || b.type.startsWith('wedge')), use = solid ? blocks.filter((b) => !/^(tri|quad)/.test(b.type)) : blocks;
    const P = model.positions, nV = P.length / 3, used = new Uint8Array(nV); for (const b of use) for (let i = 0; i < b.conn.length; i++) used[b.conn[i]] = 1;
    const L = ['** Imported mesh written by the AeroSuite 26 high-fidelity bridge', `** Source: ${String(model.name).replace(/[\r\n]+/g, ' ')}`, '*NODE, NSET=NALL'];
    const mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < nV; i++) if (used[i]) { L.push(`${i + 1}, ${f(P[3 * i])}, ${f(P[3 * i + 1])}, ${f(P[3 * i + 2])}`); for (let k = 0; k < 3; k++) { if (P[3 * i + k] < mn[k]) mn[k] = P[3 * i + k]; if (P[3 * i + k] > mx[k]) mx[k] = P[3 * i + k]; } }
    let eid = 0; const shellSets = [], solidSets = [], groupEls = new Map();
    for (const b of use) {
      const set = `E_${ABAQUS_TYPE[b.type]}`; (b.type.startsWith('tri') || b.type.startsWith('quad') ? shellSets : solidSets).push(set);
      L.push(`*ELEMENT, TYPE=${ABAQUS_TYPE[b.type]}, ELSET=${set}`);
      for (let e = 0; e < b.count; e++) { const ids = []; for (let k = 0; k < b.nodesPer; k++) ids.push(b.conn[e * b.nodesPer + k] + 1); eid++; const r = rows([eid, ...ids], 10); L.push(r.join(',\n')); const g = b.group ? b.group[e] : -1; if (g >= 0) { if (!groupEls.has(g)) groupEls.set(g, []); groupEls.get(g).push(eid); } }
    }
    for (const [g, ids] of groupEls) L.push(`*ELSET, ELSET=G_${ident(model.groups[g]?.name ?? g)}`, ...rows(ids));
    const size = mx.map((v, k) => v - mn[k]), ax = size.indexOf(Math.max(...size)), tol = 0.01 * size[ax], root = [], rest = [];
    for (let i = 0; i < nV; i++) if (used[i]) (P[3 * i + ax] <= mn[ax] + tol ? root : rest).push(i);
    let far = rest[0] ?? root[0]; for (const i of rest) if (P[3 * i + ax] > P[3 * far + ax]) far = i;
    const tip = rest.filter((i) => P[3 * i + ax] >= mx[ax] - tol);
    L.push('*NSET, NSET=ROOT', ...rows(root.map((i) => i + 1)), '*NSET, NSET=TIP', ...rows((tip.length ? tip : [far]).map((i) => i + 1)), '*NSET, NSET=TIPREF', `${far + 1}`, ...inpMaterial(s.material));
    for (const set of shellSets) L.push(`*SHELL SECTION, ELSET=${set}, MATERIAL=MAT`, f(s.t_skin_mm / 1e3, 6));
    for (const set of solidSets) L.push(`*SOLID SECTION, ELSET=${set}, MATERIAL=MAT`);
    L.push('*BOUNDARY', `ROOT, 1, ${solid ? 3 : 6}`, '');
    files.push({ path: 'model.inp', text: L.join('\n') });
    const Ftot = s.load_factor * s.lift_1g_N, loads = rest.map((i) => [i, Ftot / Math.max(1, rest.length)]);
    stepOpts = { loads }; meshInfo = { source: 'imported', name: model.name, nodes: root.length + rest.length, elements: eid, clampAxis: 'xyz'[ax], clamped: root.length };
    notes.push(`Model: imported mesh "${model.name}" (${eid} elements). Nodes within 1 % of the low end of its longest axis (${'xyz'[ax]}) are clamped; ${f(Ftot, 5)} N is spread evenly over the remaining nodes in +z. Replace these with your real supports and loads before relying on the answer.`);
    warnings.push('Loads and supports on an imported mesh are generic placeholders. Check the *BOUNDARY and *CLOAD blocks against your design loads.');
  } else {
    const m = wingBoxMesh(c, s), L = ['** Wing-box shell model written by the AeroSuite 26 high-fidelity bridge', `** Half span ${f(m.info.half, 6)} m, root box ${f(m.info.box_w_root, 5)} m × ${f(m.info.box_h_root, 5)} m, material ${s.material}`, '** Units: m, kg, s, N, Pa', '*NODE, NSET=NALL'];
    m.nodes.forEach((p, i) => L.push(`${i + 1}, ${f(p[0])}, ${f(p[1])}, ${f(p[2])}`));
    for (const [set, ids] of Object.entries(m.sets)) { if (!ids.length) continue; L.push(`*ELEMENT, TYPE=${m.type}, ELSET=${set}`); for (const e of ids) L.push(`${e + 1}, ${m.elems[e].map((n) => n + 1).join(', ')}`); }
    L.push('*NSET, NSET=ROOT', ...rows(m.root.map((i) => i + 1)), '*NSET, NSET=TIP', ...rows(m.tip.map((i) => i + 1)), '*NSET, NSET=TIPREF', `${m.tipRef + 1}`, ...inpMaterial(s.material));
    const sec = (set, t) => (m.sets[set].length ? [`*SHELL SECTION, ELSET=${set}, MATERIAL=MAT`, f(t / 1e3, 6)] : []);
    L.push(...sec('SKIN_U', s.t_skin_mm), ...sec('SKIN_L', s.t_skin_mm), ...sec('SPAR_F', s.t_spar_mm), ...sec('SPAR_R', s.t_spar_mm), ...sec('RIBS', s.t_rib_mm), '*BOUNDARY', 'ROOT, 1, 6', '');
    files.push({ path: 'model.inp', text: L.join('\n') });
    const hmin = Math.min(m.info.half / m.info.nSpan, (c.struct.box_chord_frac * m.info.c_tip) / m.info.nChord, s.t_skin_mm / 1e3) / (m.type === 'S8R' ? 2 : 1), cw = Math.sqrt(mat.E / mat.rho);
    stepOpts = { loads: m.loads, stableDt: Number(((0.5 * hmin) / cw).toPrecision(2)) }; meshInfo = { source: 'wing box', ...m.info, element: m.type };
    notes.push(`Model: wing-box shell mesh, ${m.info.elements} ${m.type} elements and ${m.info.nodes} nodes — two skins, two spar webs and ${m.info.ribs} ribs, clamped at the root.`);
    notes.push(`Load: Schrenk lift distribution for ${f(m.info.halfLift_N, 5)} N on the half-wing (${f(s.load_factor, 4)} g), applied at the spar caps so its resultant acts at 25 % chord. ${f(m.info.appliedLoad_N, 5)} N acts on the mesh; the rest is the root half-cell that goes straight into the clamp.`);
    notes.push('Not included: stringers, fuel and engine inertia relief, dihedral and twist. The native structures suite uses the same rectangular cell, so differences mainly show plate, shear-lag and rib effects.');
  }
  const steps = inpSteps(s, which, stepOpts);
  for (const [name, text] of Object.entries(steps)) files.push({ path: name, text });
  if (which.includes('explicit')) warnings.push('Explicit time stepping is limited by the thinnest shell: CalculiX expands shells to solids and takes increments of a few nanoseconds, so even a millisecond of response needs hours. Use it for short impact events or run the same deck in Abaqus/Explicit with S4R elements; for gust and landing loads choose the implicit dynamic analysis.');
  const jobs = Object.keys(steps).map((n) => n.replace(/\.inp$/, ''));
  files.push({ path: 'run.sh', exec: true, text: ['#!/usr/bin/env bash', '# Runs every CalculiX job in this folder. ccx must be on the PATH (Ubuntu: sudo apt-get install calculix-ccx).', 'set -euo pipefail', 'cd "$(dirname "$0")"', 'export OMP_NUM_THREADS="${NP:-2}"', 'CCX="${CCX:-ccx}"', `for job in ${jobs.join(' ')}; do`, '  "$CCX" -i "$job" 2>&1 | tee "log.$job"', 'done', ''].join('\n') });
  return { files, notes, warnings, meshInfo, expected: jobs.flatMap((j) => [`${j}.dat`, `${j}.frd`, `${j}.sta`, `log.${j}`]), commands: ['bash run.sh'] };
}

// ---------------------------------------------------------------------------------------------
// Gmsh: parametric half-model (wing, optional fuselage) inside a far-field box
// ---------------------------------------------------------------------------------------------
export function gmshGeo(c, { nSection = 30, farfield = 12 } = {}) {
  const fs = flowState(c);
  if (!fs.hasWing) throw new Error('The Gmsh script needs a wing in the case.');
  const pl = wingPlan(c), hasFus = c.fuselage.len_m > 0 && c.fuselage.dia_m > 0, L = [];
  const xs = []; for (let i = 0; i <= nSection; i++) xs.push(0.5 * (1 - Math.cos((Math.PI * i) / nSection)));
  L.push('// Parametric half-model written by the AeroSuite 26 high-fidelity bridge.', '// Open in Gmsh (https://gmsh.info) 4.8 or later. Mesh:  gmsh aircraft.geo -3 -format su2 -o aircraft.su2', '// Physical groups become boundary names: wall, symmetry, farfield. Units: metres; x aft, y starboard, z up.', 'SetFactory("OpenCASCADE");', '',
    `hWall = ${f(pl.mac / 40, 5)};   // target size on the wing surface [m] — lower for a finer mesh`, `hFar = ${f((farfield * pl.half) / 6, 5)};    // size at the far field [m]`, `grow = ${f(2 * pl.half, 5)};    // distance over which the size grows to hFar [m]`, '');
  let p = 0, cv = 0; const loops = [];
  const section = (eta, yOverride) => {
    const st = pl.at(eta), tw = st.twist * RAD, cs = Math.cos(tw), sn = Math.sin(tw), y = yOverride ?? st.y, zoff = yOverride != null ? 0 : st.z;
    const pt = ([x, z]) => { const dx = (x - 0.25) * st.chord, dz = z * st.chord; L.push(`Point(${++p}) = {${f(st.xle + 0.25 * st.chord + dx * cs + dz * sn, 9)}, ${f(y, 9)}, ${f(zoff - dx * sn + dz * cs, 9)}};`); return p; };
    const te = pt([1, 0]), up = [], lo = [];
    for (let i = nSection - 1; i >= 1; i--) up.push(pt(nacaAt(pl.sec, xs[i]).upper));
    const le = pt([0, 0]);
    for (let i = 1; i < nSection; i++) lo.push(pt(nacaAt(pl.sec, xs[i]).lower));
    L.push(`Spline(${++cv}) = {${[te, ...up, le].join(', ')}};`, `Spline(${++cv}) = {${[le, ...lo, te].join(', ')}};`, `Curve Loop(${loops.length + 1}) = {${cv - 1}, ${cv}};`); loops.push(loops.length + 1);
  };
  L.push('// wing sections: inboard extension (so the root passes cleanly through the symmetry plane), root, tip');
  section(0, -0.05 * pl.half); section(0); section(1);
  L.push('', `Ruled ThruSections(1) = {${loops.join(', ')}};`, 'body[] = {1};');
  if (hasFus) {
    const Lf = c.fuselage.len_m, r = c.fuselage.dia_m / 2, xc = 0.25 * pl.cr + 0.08 * Lf;
    L.push('', '// fuselage: ellipsoid of the case length and diameter', 'Sphere(2) = {0, 0, 0, 1};', `Dilate {{0, 0, 0}, {${f(Lf / 2, 7)}, ${f(r, 7)}, ${f(r, 7)}}} { Volume{2}; }`, `Translate {${f(xc, 7)}, 0, 0} { Volume{2}; }`, 'u[] = BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; };', 'body[] = u[];');
  }
  const D = farfield * pl.half, x0 = -D, eps = 1e-3 * pl.half;
  L.push('', '// far-field half box and the fluid volume', `Box(100) = {${f(x0, 7)}, 0, ${f(-D, 7)}, ${f(2.5 * D, 7)}, ${f(D, 7)}, ${f(2 * D, 7)}};`, 'fluid[] = BooleanDifference{ Volume{100}; Delete; }{ Volume{body[]}; Delete; };', '',
    `e = ${f(eps, 4)};`, `sym[] = Surface In BoundingBox{${f(x0, 7)}-e, -e, ${f(-D, 7)}-e, ${f(x0 + 2.5 * D, 7)}+e, e, ${f(D, 7)}+e};`,
    `far[] = Surface In BoundingBox{${f(x0, 7)}-e, -e, ${f(-D, 7)}-e, ${f(x0, 7)}+e, ${f(D, 7)}+e, ${f(D, 7)}+e};`,
    `far[] += Surface In BoundingBox{${f(x0 + 2.5 * D, 7)}-e, -e, ${f(-D, 7)}-e, ${f(x0 + 2.5 * D, 7)}+e, ${f(D, 7)}+e, ${f(D, 7)}+e};`,
    `far[] += Surface In BoundingBox{${f(x0, 7)}-e, ${f(D, 7)}-e, ${f(-D, 7)}-e, ${f(x0 + 2.5 * D, 7)}+e, ${f(D, 7)}+e, ${f(D, 7)}+e};`,
    `far[] += Surface In BoundingBox{${f(x0, 7)}-e, -e, ${f(-D, 7)}-e, ${f(x0 + 2.5 * D, 7)}+e, ${f(D, 7)}+e, ${f(-D, 7)}+e};`,
    `far[] += Surface In BoundingBox{${f(x0, 7)}-e, -e, ${f(D, 7)}-e, ${f(x0 + 2.5 * D, 7)}+e, ${f(D, 7)}+e, ${f(D, 7)}+e};`,
    'wall[] = Surface{:};', 'wall[] -= sym[];', 'wall[] -= far[];', '', 'Physical Surface("wall") = {wall[]};', 'Physical Surface("symmetry") = {sym[]};', 'Physical Surface("farfield") = {far[]};', 'Physical Volume("fluid") = {fluid[]};', '',
    '// mesh size: fine on the wall, growing with distance', 'Field[1] = Distance;', 'Field[1].SurfacesList = {wall[]};', 'Field[1].Sampling = 60;', 'Field[2] = Threshold;', 'Field[2].InField = 1;', 'Field[2].SizeMin = hWall;', 'Field[2].SizeMax = hFar;', 'Field[2].DistMin = 0;', 'Field[2].DistMax = grow;', 'Background Field = 2;',
    'Mesh.MeshSizeExtendFromBoundary = 0;', 'Mesh.MeshSizeFromPoints = 0;', 'Mesh.MeshSizeFromCurvature = 0;', 'Mesh.Algorithm = 6;', 'Mesh.Algorithm3D = 1;', '');
  return L.join('\n');
}

// ---------------------------------------------------------------------------------------------
// Whole case: files + manifest
// ---------------------------------------------------------------------------------------------
const README = (b) => [`${b.title}`, '='.repeat(Math.min(78, b.title.length)), '', 'Prepared by the AeroSuite 26 high-fidelity bridge. Everything needed to run is in this folder.', '', 'Run it', '------', ...b.commands.map((x) => `  ${x}`), '', 'Then bring these files back to the bridge page (drop them, or the whole folder zipped):', ...b.expected.map((x) => `  ${x}`), '', 'What is in the deck', '-------------------', ...b.notes.map((x) => `- ${x}`), ...(b.warnings.length ? ['', 'Check before relying on the result', '----------------------------------', ...b.warnings.map((x) => `- ${x}`)] : []), '', `Case hash (SHA-256): ${b.hash}`, ''].join('\n');

/**
 * Build a complete solver case.
 * @returns {{name, solver, analysis, files: {path, text, exec?}[], manifest, notes, warnings, commands, expected, meshInfo, settings}}
 */
export function buildCase(c, { solver, settings = null, model = null, appVersion = '', up = {}, includeGeo = true, created = null } = {}) {
  if (!SOLVERS[solver]) throw new Error(`Unknown solver "${solver}".`);
  const s = { ...defaultSettings(c, solver, { analysis: settings?.analysis, model, up }), ...(settings || {}) };
  if (!SOLVERS[solver].analyses[s.analysis]) throw new Error(`Unknown ${SOLVERS[solver].label} analysis "${s.analysis}".`);
  for (const fd of SETTING_FIELDS[solver]) { if (fd.type) continue; const v = s[fd.key]; if (typeof v !== 'number' || !Number.isFinite(v)) throw new Error(`"${fd.label}" needs a number.`); if (fd.min != null && v < fd.min) throw new Error(`"${fd.label}" must be at least ${fd.min}.`); if (fd.max != null && v > fd.max) throw new Error(`"${fd.label}" must be at most ${fd.max}.`); }
  const part = solver === 'su2' ? su2Case(c, s, model) : solver === 'openfoam' ? openfoamCase(c, s, model) : calculixCase(c, s, model);
  const hash = caseHash(c), name = `${slug(c.meta.name)}-${solver}-${slug(s.analysis)}`, title = `${c.meta.name || 'Case'} — ${SOLVERS[solver].label}: ${SOLVERS[solver].analyses[s.analysis]}`;
  const files = part.files.slice();
  if (includeGeo && solver !== 'calculix' && c.wing.S_m2 > 0 && c.wing.b_m > 0) { files.push({ path: 'geometry/aircraft.geo', text: gmshGeo(c) }); part.notes.push('Extra: geometry/aircraft.geo is a parametric Gmsh script of the wing (and fuselage) in a far-field box, for meshing the 3-D shape at any resolution.'); }
  files.push({ path: 'case.json', text: JSON.stringify(c, null, 1) });
  files.push({ path: 'README.txt', text: README({ title, commands: part.commands, expected: part.expected, notes: part.notes, warnings: part.warnings, hash }) });
  const manifest = {
    app: 'AeroSuite 26', bridge: 1, appVersion: appVersion || 'unknown', created: created || new Date().toISOString(), name, title, caseName: c.meta.name || '', caseHash: hash,
    solver, solverLabel: SOLVERS[solver].label, analysis: s.analysis, settings: s, mesh: part.meshInfo, model: model ? { name: model.name, format: model.formatName || model.format || '' } : null,
    reference: { ...(({ V, mach, reynolds, rho, mu, T, p, S, b, refLength, chord }) => ({ V_ms: V, mach, reynolds, rho, mu, T_K: T, p_Pa: p, S_m2: S, b_m: b, refLength_m: refLength, chord_m: chord }))(flowState(c)), weight_N: derived(c).W, g0: G0 },
    run: part.commands, decks: files.map((x) => ({ path: x.path, bytes: new TextEncoder().encode(x.text).length, sha256: sha256(x.text) })), expectedOutputs: part.expected, notes: part.notes, warnings: part.warnings,
  };
  files.push({ path: 'manifest.json', text: JSON.stringify(manifest, null, 1) });
  return { name, title, solver, analysis: s.analysis, files, manifest, notes: part.notes, warnings: part.warnings, commands: part.commands, expected: part.expected, meshInfo: part.meshInfo, settings: s };
}
