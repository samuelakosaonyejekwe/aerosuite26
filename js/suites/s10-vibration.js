// Suite 10 — Vibration, Modal and Rotor Dynamics.
// Beam finite-element modal analysis of the wing and of rotating blades (fan diagram), forced response and
// isolation, rotor-bearing whirl with gyroscopics (Campbell, critical speeds, unbalance orbit), drivetrain
// torsional vibration (matrix eigen-solution cross-checked by Holzer), random vibration and ground resonance.

import * as N from '../core/numerics.js';
import { G0 } from '../core/atmosphere.js';
import { METALS } from '../data/materials.js';

const TAU = 2 * Math.PI, rads = (rpm) => (rpm * TAU) / 60, fin = Number.isFinite;
const thin = (a, n = 300) => { const s = Math.max(1, Math.ceil(a.length / n)); return a.filter((_, k) => k % s === 0 || k === a.length - 1); };

// ---- shared numerics ------------------------------------------------------------------------
/** K x = λ M x with K scaled to O(1) so the kernel's absolute Jacobi tolerance acts as a relative one. */
function geig(K, M) {
  let s = 0; for (let k = 0; k < K.length; k++) s = Math.max(s, Math.abs(K[k][k] / M[k][k]));
  s = s || 1;
  const e = N.eigGenSym(N.mscale(K, 1 / s), M);
  return { values: e.values.map((v) => v * s), vectors: e.vectors };
}
/**
 * Euler–Bernoulli beam modes by Hermitian finite elements (w, w' per node), properties sampled at element mid-points.
 * o: root 'clamped'|'pinned', kRoot (rotational spring at a pinned root), omega (spin rate), r0 (root radius from the
 * spin axis), inplane (lead-lag: centrifugal softening −Ω²M), tipMass, points [{x, m}].
 * Returns angular frequencies w[], nodal deflection shapes[][] (M-orthonormal), x[], effective-mass fractions.
 */
export function beamModes(n, L, EIf, mf, o = {}) {
  const nd = 2 * (n + 1), K = N.zeros(nd), M = N.zeros(nd), l = L / n, Om = o.omega || 0, r0 = o.r0 || 0, T = new Array(n).fill(0);
  if (Om) { // element-average centrifugal tension, accumulated from the tip inboard
    let acc = (o.tipMass || 0) * Om * Om * (r0 + L);
    for (let e = n - 1; e >= 0; e--) { const x1 = e * l, x2 = x1 + l, m = mf(x1 + l / 2) * Om * Om; T[e] = acc + m * ((r0 * l) / 2 + (x2 * x2 * l - (x2 ** 3 - x1 ** 3) / 3) / (2 * l)); acc += m * l * (r0 + x1 + l / 2); }
  }
  const ke = [[12, 6 * l, -12, 6 * l], [6 * l, 4 * l * l, -6 * l, 2 * l * l], [-12, -6 * l, 12, -6 * l], [6 * l, 2 * l * l, -6 * l, 4 * l * l]];
  const me = [[156, 22 * l, 54, -13 * l], [22 * l, 4 * l * l, 13 * l, -3 * l * l], [54, 13 * l, 156, -22 * l], [-13 * l, -3 * l * l, -22 * l, 4 * l * l]];
  const kg = [[36, 3 * l, -36, 3 * l], [3 * l, 4 * l * l, -3 * l, -l * l], [-36, -3 * l, 36, -3 * l], [3 * l, -l * l, -3 * l, 4 * l * l]];
  for (let e = 0; e < n; e++) {
    const xm = (e + 0.5) * l, a = EIf(xm) / l ** 3, b = (mf(xm) * l) / 420, g = T[e] / (30 * l);
    for (let p = 0; p < 4; p++) for (let q = 0; q < 4; q++) { K[2 * e + p][2 * e + q] += a * ke[p][q] + g * kg[p][q]; M[2 * e + p][2 * e + q] += b * me[p][q]; }
  }
  if (o.tipMass) M[nd - 2][nd - 2] += o.tipMass;
  for (const p of o.points || []) if (p.m > 0) { const k = 2 * N.clamp(Math.round(p.x / l), 1, n); M[k][k] += p.m; }
  if (o.kRoot) K[1][1] += o.kRoot;
  if (o.inplane && Om) for (let p = 0; p < nd; p++) for (let q = 0; q < nd; q++) K[p][q] -= Om * Om * M[p][q];
  const free = o.root === 'pinned' ? N.range(nd - 1, (k) => k + 1) : N.range(nd - 2, (k) => k + 2);
  const sub = (A) => free.map((p) => free.map((q) => A[p][q])), Mr = sub(M), e = geig(sub(K), Mr);
  const r = free.map((k) => (k % 2 ? 0 : 1)), Mrv = N.matvec(Mr, r), mTot = N.dot(r, Mrv);
  const shapes = e.vectors.map((v) => { const w = new Array(n + 1).fill(0); free.forEach((k, j) => { if (k % 2 === 0) w[k / 2] = v[j]; }); return w; });
  return { w: e.values.map((v) => Math.sqrt(Math.max(v, 0))), lam: e.values, shapes, x: N.range(n + 1, (k) => k * l), meff: e.vectors.map((v) => N.dot(v, Mrv) ** 2 / mTot), vectors: e.vectors, Mr };
}
/** Fixed–free torsion rod modes with linear two-node elements. Returns angular frequencies and twist shapes. */
export function torsionModes(n, L, GJf, If) {
  const K = N.zeros(n), M = N.zeros(n), l = L / n; // node 0 clamped; unknowns are nodes 1..n
  for (let e = 0; e < n; e++) {
    const xm = (e + 0.5) * l, k = GJf(xm) / l, m = (If(xm) * l) / 6, a = e - 1, b = e;
    if (a >= 0) { K[a][a] += k; K[a][b] -= k; K[b][a] -= k; M[a][a] += 2 * m; M[a][b] += m; M[b][a] += m; }
    K[b][b] += k; M[b][b] += 2 * m;
  }
  const e = geig(K, M);
  return { w: e.values.map((v) => Math.sqrt(Math.max(v, 0))), shapes: e.vectors.map((v) => [0, ...v]), x: N.range(n + 1, (k) => k * l) };
}
/**
 * Eigenvalues of a general real or complex matrix (entries numbers or [re, im]) by shifted QR on the Hessenberg form.
 * Same algorithm as the kernel's eig(), on flat arrays without per-operation allocation, for the many small solves here.
 */
export function ceig(A) {
  const n = A.length, hr = new Float64Array(n * n), hi = new Float64Array(n * n), cr = new Float64Array(n), ci = new Float64Array(n), sr = new Float64Array(n), si = new Float64Array(n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { const v = A[i][j]; if (typeof v === 'number') hr[i * n + j] = v; else { hr[i * n + j] = v[0]; hi[i * n + j] = v[1]; } }
  const ab = (x, y) => Math.sqrt(x * x + y * y);
  const swap = (a, b) => { let t = hr[a]; hr[a] = hr[b]; hr[b] = t; t = hi[a]; hi[a] = hi[b]; hi[b] = t; };
  for (let m = 1; m < n - 1; m++) { // Hessenberg reduction by stabilised elementary similarity transformations
    let p = m, mx = 0;
    for (let i = m; i < n; i++) { const v = ab(hr[i * n + m - 1], hi[i * n + m - 1]); if (v > mx) { mx = v; p = i; } }
    if (mx === 0) continue;
    if (p !== m) { for (let j = 0; j < n; j++) swap(p * n + j, m * n + j); for (let j = 0; j < n; j++) swap(j * n + p, j * n + m); }
    const pr = hr[m * n + m - 1], pi = hi[m * n + m - 1], pd = pr * pr + pi * pi;
    for (let i = m + 1; i < n; i++) {
      const ar = hr[i * n + m - 1], ai = hi[i * n + m - 1]; if (ar === 0 && ai === 0) continue;
      const fr = (ar * pr + ai * pi) / pd, fi = (ai * pr - ar * pi) / pd;
      for (let j = 0; j < n; j++) { const br = hr[m * n + j], bi = hi[m * n + j]; hr[i * n + j] -= fr * br - fi * bi; hi[i * n + j] -= fr * bi + fi * br; }
      for (let j = 0; j < n; j++) { const br = hr[j * n + i], bi = hi[j * n + i]; hr[j * n + m] += fr * br - fi * bi; hi[j * n + m] += fr * bi + fi * br; }
    }
  }
  const vals = []; let nn = n, iter = 0;
  while (nn > 0) {
    if (nn === 1) { vals.push([hr[0], hi[0]]); break; }
    let l = nn - 1;
    for (; l > 0; l--) { const d = ab(hr[(l - 1) * n + l - 1], hi[(l - 1) * n + l - 1]) + ab(hr[l * n + l], hi[l * n + l]) || 1; if (ab(hr[l * n + l - 1], hi[l * n + l - 1]) < 1e-14 * d) { hr[l * n + l - 1] = 0; hi[l * n + l - 1] = 0; break; } }
    if (l === nn - 1) { vals.push([hr[l * n + l], hi[l * n + l]]); nn--; iter = 0; continue; }
    if (++iter > 500) throw new Error('ceig: QR iteration did not converge');
    const q = (nn - 2) * n + nn - 2, ar = hr[q], ai = hi[q], br = hr[q + 1], bi = hi[q + 1], er = hr[q + n], ei = hi[q + n], dr = hr[q + n + 1], di = hi[q + n + 1];
    let mr, mi; // Wilkinson shift (exceptional shift every eleventh iteration)
    if (iter % 11 === 10) { mr = ab(er, ei) + (nn > 2 ? ab(hr[q - 1], hi[q - 1]) : 0); mi = 0; }
    else {
      const tr = ar + dr, ti = ai + di, xr = tr * tr - ti * ti - 4 * (ar * dr - ai * di - br * er + bi * ei), xi = 2 * tr * ti - 4 * (ar * di + ai * dr - br * ei - bi * er);
      const r = ab(xr, xi), zr = Math.sqrt((r + xr) / 2), zi = (xi < 0 ? -1 : 1) * Math.sqrt(Math.max(0, (r - xr) / 2));
      const l1r = (tr + zr) / 2, l1i = (ti + zi) / 2, l2r = (tr - zr) / 2, l2i = (ti - zi) / 2;
      if (ab(l1r - dr, l1i - di) < ab(l2r - dr, l2i - di)) { mr = l1r; mi = l1i; } else { mr = l2r; mi = l2i; }
    }
    for (let i = l; i < nn; i++) { hr[i * n + i] -= mr; hi[i * n + i] -= mi; }
    for (let k = l; k < nn - 1; k++) { // Givens rotations from the left
      const xr = hr[k * n + k], xi = hi[k * n + k], yr = hr[(k + 1) * n + k], yi = hi[(k + 1) * n + k], r = Math.sqrt(xr * xr + xi * xi + yr * yr + yi * yi);
      if (r === 0) { cr[k] = NaN; continue; }
      const c1 = xr / r, c2 = xi / r, s1 = yr / r, s2 = yi / r; cr[k] = c1; ci[k] = c2; sr[k] = s1; si[k] = s2;
      for (let j = k; j < nn; j++) {
        const a = k * n + j, b = a + n, ur = hr[a], ui = hi[a], vr = hr[b], vi = hi[b];
        hr[a] = c1 * ur + c2 * ui + s1 * vr + s2 * vi; hi[a] = c1 * ui - c2 * ur + s1 * vi - s2 * vr;
        hr[b] = c1 * vr - c2 * vi - s1 * ur + s2 * ui; hi[b] = c1 * vi + c2 * vr - s1 * ui - s2 * ur;
      }
    }
    for (let k = l; k < nn - 1; k++) { // and from the right
      if (Number.isNaN(cr[k])) continue;
      const c1 = cr[k], c2 = ci[k], s1 = sr[k], s2 = si[k], top = Math.min(k + 2, nn - 1);
      for (let i = l; i <= top; i++) {
        const a = i * n + k, ur = hr[a], ui = hi[a], vr = hr[a + 1], vi = hi[a + 1];
        hr[a] = ur * c1 - ui * c2 + vr * s1 - vi * s2; hi[a] = ur * c2 + ui * c1 + vr * s2 + vi * s1;
        hr[a + 1] = vr * c1 + vi * c2 - ur * s1 - ui * s2; hi[a + 1] = vi * c1 - vr * c2 - ui * s1 + ur * s2;
      }
    }
    for (let i = l; i < nn; i++) { hr[i * n + i] += mr; hi[i * n + i] += mi; }
  }
  return vals.map((v) => [v[0], Math.abs(v[1]) < 1e-12 * (1 + Math.abs(v[0])) ? 0 : v[1]]).sort((p, q) => p[1] - q[1] || p[0] - q[0]);
}
const unit = (s) => { const m = N.amax(s.map(Math.abs)) || 1, sg = s[s.length - 1] < 0 ? -1 : 1; return s.map((v) => (sg * v) / m); };
/** Complex harmonic solve (K − ω²M + iωC) X = F for real matrices. */
const harmonic = (M, C, K, F, w) => N.csolve(K.map((r, p) => r.map((k, q) => [k - w * w * M[p][q], w * C[p][q]])), F);

// ---- preliminary structural estimates from the shared case ----------------------------------
/** Thin-walled two-spar wing-box stiffness at the root and the wing mass model (preliminary, constant gauge). */
function wingDefaults(c, up, d) {
  if (!(c.wing.S_m2 > 0)) return {};
  const mat = METALS[c.struct.material] || METALS['Al 2024-T3'], cr = d.c_root, w = c.struct.box_chord_frac * cr, h = c.struct.box_height_frac * c.wing.tc * cr;
  // same gauge caps as the structures suite (0.06% and 0.1% of the root chord) so that small vehicles are not given airliner gauges
  const ts = Math.min(c.struct.t_skin_mm, Math.max(0.4, 0.6 * cr)) / 1e3, tw = Math.min(c.struct.t_spar_mm, Math.max(0.5, cr)) / 1e3;
  const EI = mat.E * (2 * 1.6 * w * ts * (h / 2) ** 2 + (2 * tw * h ** 3) / 12); // stringers smeared as 60% extra skin area
  const GJ = (mat.G * 4 * (w * h) ** 2) / ((2 * w) / ts + (2 * h) / tw);
  return {
    semi_span: c.wing.b_m / 2 / Math.cos(N.rad(c.wing.sweep_deg)), c_root: cr, taper: c.wing.taper,
    EI_root: up.fea?.EI_root_Nm2 ?? EI, GJ_root: up.fea?.GJ_root_Nm2 ?? GJ,
    m_wing: up.fea?.wing_struct_mass_kg > 0 ? up.fea.wing_struct_mass_kg / 2 + 0.02 * c.mass.mtow_kg : 0.06 * c.mass.mtow_kg, m_fuel: 0.3 * c.mass.fuel_kg, zeta: c.struct.zeta,
  };
}
/** Rotor blade (or propeller blade) properties; stiffnesses are section-shape estimates to be replaced by blade data. */
function bladeDefaults(c) {
  const r = c.rotor, heli = c.meta.type === 'helicopter', rot = r.R_m > 0, R = rot ? r.R_m : c.prop.prop_dia_m / 2;
  if (!(R > 0)) return {};
  const chord = rot && r.chord_m > 0 ? r.chord_m : 0.14 * R, t = 0.12 * chord, E = 70e9;
  return {
    R, chord, n_blades: (rot ? r.n_blades : c.prop.n_blades) || 2, rpm: (rot ? r.rpm : c.prop.rpm) || undefined,
    m_blade: rot && r.blade_mass_kg > 0 ? r.blade_mass_kg : 163 * chord * chord * R, e_root: heli ? r.hinge_offset : 0.12, root: heli ? 'hinged' : 'cantilever',
    EI_flap: E * 0.011 * chord * t ** 3, EI_lag: E * 0.0135 * t * chord ** 3,
  };
}
const hasRotor = (c) => (c.rotor.R_m > 0 || c.prop.prop_dia_m > 0 ? true : 'This analysis needs a rotor or propeller; the case has neither.');
const hasWing = (c) => (c.wing.S_m2 > 0 ? true : 'This analysis needs a wing; use the rotating-blade analysis for rotor and propeller blades.');

// ---- 1. wing modal analysis -----------------------------------------------------------------
const WING_INPUTS = [
  { key: 'semi_span', label: 'Structural semi-span', unit: 'm', default: 17, min: 0.1, group: 'Geometry', help: 'Root to tip along the elastic axis (semi-span / cos sweep)' },
  { key: 'c_root', label: 'Root chord', unit: 'm', default: 5.8, min: 0.02, group: 'Geometry' },
  { key: 'taper', label: 'Taper ratio', unit: '-', default: 0.3, min: 0.05, max: 1, group: 'Geometry' },
  { key: 'EI_root', label: 'Root bending stiffness EI', unit: 'N·m²', default: 1e8, min: 1e-3, group: 'Structure', help: 'From the structures suite when available, else a two-spar box estimate' },
  { key: 'GJ_root', label: 'Root torsional stiffness GJ', unit: 'N·m²', default: 1.2e8, min: 1e-3, group: 'Structure' },
  { key: 'stiff_exp', label: 'Stiffness taper exponent', unit: '-', default: 3, min: 0, max: 5, group: 'Structure', help: 'EI, GJ ∝ chord^p: 3 for a constant-gauge box, 4 when gauge scales with chord' },
  { key: 'm_wing', label: 'Semi-wing mass (structure and systems)', unit: 'kg', default: 4700, min: 1e-4, group: 'Mass', help: 'Distributed ∝ chord²' },
  { key: 'm_fuel', label: 'Fuel in the semi-wing', unit: 'kg', default: 5600, min: 0, group: 'Mass', help: 'Distributed ∝ chord²; lowers every frequency' },
  { key: 'm_point', label: 'Concentrated mass (engine, pod, tip tank)', unit: 'kg', default: 0, min: 0, group: 'Mass' },
  { key: 'eta_point', label: 'Concentrated mass position', unit: '-', default: 0.33, min: 0.05, max: 1, group: 'Mass', help: 'Fraction of the semi-span' },
  { key: 'r_alpha', label: 'Section radius of gyration / chord', unit: '-', default: 0.25, min: 0.1, max: 0.5, group: 'Mass', help: 'About the elastic axis; 0.23–0.30 for wings' },
];
function wingBeam(i, n) {
  const ch = (x) => i.c_root * (1 - (1 - i.taper) * (x / i.semi_span)), s = (x) => (ch(x) / i.c_root) ** i.stiff_exp;
  const m0 = (3 * (i.m_wing + i.m_fuel)) / (i.semi_span * (1 + i.taper + i.taper ** 2)), m = (x) => m0 * (ch(x) / i.c_root) ** 2;
  const b = beamModes(n, i.semi_span, (x) => i.EI_root * s(x), m, { points: [{ x: i.eta_point * i.semi_span, m: i.m_point }] });
  const t = torsionModes(n, i.semi_span, (x) => i.GJ_root * s(x), (x) => m(x) * (i.r_alpha * ch(x)) ** 2);
  return { b, t, m0 };
}
const wingModal = {
  id: 'wing', title: 'Wing bending and torsion modes (beam finite elements)', fidelity: 'numerical',
  summary: 'Natural frequencies and mode shapes of the cantilevered wing from a tapered beam and torsion-rod finite-element model, with modal effective mass.',
  equations: ['Eigenvalue equations', 'Undamped free vibration equations', 'Modal orthogonality relations', 'Euler–Lagrange equations'],
  applicable: hasWing,
  inputs: [...WING_INPUTS, { key: 'nEl', label: 'Beam elements', unit: '', default: 20, min: 3, max: 96, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => wingDefaults(c, up, d),
  run(i) {
    const n = Math.round(i.nEl), { b, t } = wingBeam(i, n), f = b.w.map((w) => w / TAU), ft = t.w.map((w) => w / TAU), eta = b.x.map((x) => x / i.semi_span);
    // modal orthogonality: largest off-diagonal of ΦᵀMΦ for the first six modes
    let ortho = 0; for (let p = 0; p < Math.min(6, f.length); p++) { const Mv = N.matvec(b.Mr, b.vectors[p]); for (let q = 0; q < p; q++) ortho = Math.max(ortho, Math.abs(N.dot(b.vectors[q], Mv))); }
    const warnings = [];
    if (ft[0] / f[0] < 2.5) warnings.push('First torsion is less than 2.5 times first bending: bending–torsion flutter coupling is likely; check Suite 3.');
    if (i.semi_span / i.c_root < 2) warnings.push('Low aspect ratio: beam theory misses chordwise (plate) modes; use a shell model.');
    return {
      kpis: [
        { key: 'f1_Hz', label: 'First bending frequency', value: f[0], unit: 'Hz' }, { key: 'f2_Hz', label: 'Second bending frequency', value: f[1], unit: 'Hz' }, { key: 'f3_Hz', label: 'Third bending frequency', value: f[2], unit: 'Hz' },
        { key: 'f_torsion_Hz', label: 'First torsion frequency', value: ft[0], unit: 'Hz' }, { key: 'f_torsion2_Hz', label: 'Second torsion frequency', value: ft[1], unit: 'Hz' },
        { key: 'freq_ratio_tb', label: 'Torsion / bending frequency ratio', value: ft[0] / f[0], unit: '-', status: ft[0] / f[0] > 2.5 ? 'ok' : 'warn', note: 'Well-separated modes (ratio above about 2.5) delay classical flutter' },
        { key: 'meff_1_pct', label: 'Effective mass of mode 1', value: 100 * b.meff[0], unit: '%' },
        { key: 'ortho_err', label: 'Mass-orthogonality error', value: ortho, unit: '-', status: ortho < 1e-6 ? 'ok' : 'warn', note: 'Largest off-diagonal of ΦᵀMΦ' },
      ],
      plots: [
        { type: 'line', title: 'Bending mode shapes', xlabel: 'Span fraction [-]', ylabel: 'Normalised deflection [-]', series: [0, 1, 2, 3].filter((k) => k < f.length).map((k) => ({ name: `Mode ${k + 1}: ${f[k].toPrecision(3)} Hz`, x: eta, y: unit(b.shapes[k]) })) },
        { type: 'line', title: 'Torsion mode shapes', xlabel: 'Span fraction [-]', ylabel: 'Normalised twist [-]', series: [0, 1, 2].filter((k) => k < ft.length).map((k) => ({ name: `Mode ${k + 1}: ${ft[k].toPrecision(3)} Hz`, x: eta, y: unit(t.shapes[k]) })) },
        { type: 'bar', title: 'Modal effective mass (vertical base motion)', ylabel: 'Share of wing mass [%]', categories: ['1', '2', '3', '4', '5'].slice(0, Math.min(5, f.length)), series: [{ name: 'Effective mass', y: b.meff.slice(0, 5).map((v) => 100 * v) }] },
      ],
      tables: [{ title: 'Natural frequencies', columns: ['Mode', 'Bending [Hz]', 'Torsion [Hz]'], rows: [0, 1, 2, 3].filter((k) => k < ft.length).map((k) => [k + 1, f[k], ft[k]]) }],
      outputs: { wing_modes_Hz: f.slice(0, 6), wing_torsion_Hz: ft.slice(0, 3) },
      warnings, models: ['Finite element modal model (Hermitian beam, linear torsion rod)', 'Consistent mass matrices', 'Jacobi generalised eigen-solver'],
      assumptions: ['Wing clamped at the root; fuselage flexibility and free-free aircraft modes are not included', 'Bending and torsion are uncoupled (no mass or elastic offset coupling)', 'Stiffness ∝ chord^p and mass ∝ chord² along the span', 'No shear deformation or rotary inertia (slender beam)', 'Default stiffness (when the structures suite has not run), mass fractions, stiffness exponent and radius of gyration are typical estimates'],
    };
  },
  convergence: { param: 'nEl', label: 'Beam elements', levels: [4, 8, 16, 32], metric: 'f1_Hz' },
  calibration: { params: [{ key: 'EI_root', min: 1e3, max: 1e11 }, { key: 'm_fuel', min: 0, max: 1e5 }], sweep: 'm_fuel', target: 'f1_Hz', note: 'Ground-vibration-test first bending frequency at several fuel states' },
  verify() {
    const EI = 2e6, m = 12, L = 5, GJ = 8e5, Ia = 1.5, b = beamModes(32, L, () => EI, () => m), t = torsionModes(64, L, () => GJ, () => Ia), k = Math.sqrt(EI / (m * L ** 4));
    return [
      N.check('Uniform cantilever mode 1 (β₁L = 1.8751)', b.w[0], 1.875104 ** 2 * k, 1e-5, 'Euler–Bernoulli exact'),
      N.check('Uniform cantilever mode 2 (β₂L = 4.6941)', b.w[1], 4.694091 ** 2 * k, 1e-4, 'Euler–Bernoulli exact'),
      N.check('Uniform cantilever mode 3 (β₃L = 7.8548)', b.w[2], 7.854757 ** 2 * k, 5e-4, 'Euler–Bernoulli exact'),
      N.check('Fixed–free torsion rod mode 1', t.w[0], (Math.PI / (2 * L)) * Math.sqrt(GJ / Ia), 2e-4, 'Exact (π/2L)·sqrt(GJ/Iα)'),
      N.check('Effective masses sum to the free mass', N.sum(b.meff), 1, 1e-6, 'Completeness of the modal basis'),
      N.check('Default box stiffness of a small wing uses the capped gauges (0.4 mm skin, 0.5 mm web at a 0.5 m chord)', wingDefaults({ wing: { S_m2: 1, b_m: 4, tc: 0.12, taper: 1, sweep_deg: 0 }, struct: { material: 'Al 2024-T3', box_chord_frac: 0.5, box_height_frac: 1, t_skin_mm: 3, t_spar_mm: 5, zeta: 0.02 }, mass: { mtow_kg: 20, fuel_kg: 0 } }, {}, { c_root: 0.5 }).EI_root, METALS['Al 2024-T3'].E * (2 * 1.6 * 0.25 * 0.0004 * 0.03 ** 2 + (2 * 0.0005 * 0.06 ** 3) / 12), 1e-12, 'Thin-walled box second moment of area with the gauge caps of the structures suite'),
    ];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (o.freq_ratio_tb < 2.5) out.push({ severity: 'warn', title: 'Bending and torsion frequencies are close', detail: `Torsion/bending ratio is ${o.freq_ratio_tb.toFixed(2)}.`, action: 'Raise torsional stiffness (closed box, ±45° plies, thicker skins) or move mass forward; confirm the flutter margin in Suite 3.', basis: 'Frequency separation rule for classical bending–torsion flutter' });
    out.push({ severity: 'info', title: 'Confirm with a ground vibration test', detail: `Predicted first bending ${o.f1_Hz.toFixed(2)} Hz and first torsion ${o.f_torsion_Hz.toFixed(2)} Hz are beam-model estimates.`, action: 'Calibrate EI, GJ and the mass distribution against measured modes before using the frequencies for flutter clearance.', basis: 'Model updating practice' });
    return out;
  },
};

// ---- 2. rotating blade: Southwell stiffening and fan diagram ----------------------------------
const bladeFreqs = (i, Om, n) => {
  const r0 = i.e_root * i.R, L = i.R - r0, m = () => i.m_blade / L, root = i.root === 'hinged' ? 'pinned' : 'clamped';
  const flap = beamModes(n, L, () => i.EI_flap, m, { omega: Om, r0, root, tipMass: i.m_tip });
  const lag = beamModes(n, L, () => i.EI_lag, m, { omega: Om, r0, root, tipMass: i.m_tip, inplane: true });
  return { flap, lag };
};
const blade = {
  id: 'blade', title: 'Rotating blade frequencies and fan (Southwell) diagram', fidelity: 'numerical',
  summary: 'Flap and lead-lag bending frequencies of a rotor or propeller blade against rotor speed, with centrifugal stiffening, per-rev excitation lines and resonance crossings.',
  equations: ['Eigenvalue equations', 'Campbell frequency relations', 'Linear multi-degree-of-freedom vibration equations'],
  applicable: hasRotor,
  inputs: [
    { key: 'R', label: 'Blade tip radius', unit: 'm', default: 8, min: 0.02, group: 'Blade' },
    { key: 'e_root', label: 'Root or hinge radius / R', unit: '-', default: 0.05, min: 0, max: 0.5, group: 'Blade', help: 'Flap/lag hinge offset for an articulated rotor, hub cut-out for a hingeless blade' },
    { key: 'root', label: 'Root attachment', type: 'select', options: ['hinged', 'cantilever'], default: 'hinged', group: 'Blade', help: 'Hinged = articulated rotor; cantilever = hingeless rotor or propeller' },
    { key: 'm_blade', label: 'Blade mass', unit: 'kg', default: 110, min: 1e-4, group: 'Blade' },
    { key: 'm_tip', label: 'Tip mass', unit: 'kg', default: 0, min: 0, group: 'Blade' },
    { key: 'EI_flap', label: 'Flapwise stiffness EI', unit: 'N·m²', default: 9e4, min: 1e-6, group: 'Blade', help: 'Default is a section-shape estimate; enter measured blade data' },
    { key: 'EI_lag', label: 'Chordwise (lag) stiffness EI', unit: 'N·m²', default: 7e6, min: 1e-6, group: 'Blade' },
    { key: 'n_blades', label: 'Number of blades', unit: '', default: 4, min: 2, max: 12, step: 1, discrete: true, group: 'Blade' },
    { key: 'rpm', label: 'Operating speed', unit: 'rpm', default: 258, min: 1, group: 'Operation' },
    { key: 'rpm_max_frac', label: 'Diagram top speed / operating', unit: '-', default: 1.25, min: 1, max: 2, group: 'Operation' },
    { key: 'margin_req', label: 'Required frequency separation', unit: '%', default: 5, min: 1, max: 30, group: 'Limits', help: 'Typical design practice keeps blade modes 5–15% clear of each excitation harmonic at operating speed' },
    { key: 'nEl', label: 'Beam elements', unit: '', default: 10, min: 4, max: 48, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c) => bladeDefaults(c),
  run(i, ctx) {
    const n = Math.round(i.nEl), Om0 = rads(i.rpm), nb = Math.round(i.n_blades), F = [[], [], []], Lg = [[], []];
    const sp = [...N.linspace(0, i.rpm, 9), ...N.linspace(i.rpm, i.rpm * i.rpm_max_frac, 3).slice(1)], all = sp.map((s) => bladeFreqs(i, rads(s), n)); // grid contains rest and operating speed
    for (const r of all) { for (let k = 0; k < 3; k++) F[k].push(r.flap.w[k] / TAU); for (let k = 0; k < 2; k++) Lg[k].push(r.lag.w[k] / TAU); }
    const op = all[8], nr = all[0], rev = Om0 / TAU, perRev = (w) => w / Om0;
    const modes = [...[0, 1, 2].map((k) => ({ name: `Flap ${k + 1}`, f: op.flap.w[k] / TAU, f0: nr.flap.w[k] / TAU, hist: F[k] })), ...[0, 1].map((k) => ({ name: `Lag ${k + 1}`, f: op.lag.w[k] / TAU, f0: nr.lag.w[k] / TAU, hist: Lg[k] }))];
    // separation from every integer harmonic at operating speed, and crossings inside the speed range
    const nH = Math.max(8, nb + 2), rows = [], cross = []; let worst = { sep: Infinity };
    for (const m of modes) {
      const h = Math.max(1, Math.round(m.f / rev)), sep = Math.abs(m.f - h * rev) / (h * rev);
      rows.push([m.name, m.f0, m.f, m.f / rev, (m.f ** 2 - m.f0 ** 2) / rev ** 2, `${h}/rev`, 100 * sep]);
      // assessed: harmonics up to nH; the first flap mode of a hinged blade sits near 1/rev by nature and is aerodynamically well damped
      if (m.f / rev > 0.6 && h <= nH && !(m.name === 'Flap 1' && i.root === 'hinged') && sep < worst.sep) worst = { sep, name: m.name, h };
      for (let hh = 1; hh <= nH; hh++) for (let k = 1; k < sp.length; k++) { const g0 = m.hist[k - 1] - (hh * sp[k - 1]) / 60, g1 = m.hist[k] - (hh * sp[k]) / 60; if (g0 * g1 < 0) cross.push([m.name, `${hh}/rev`, sp[k - 1] + ((sp[k] - sp[k - 1]) * g0) / (g0 - g1), (100 * (sp[k - 1] + ((sp[k] - sp[k - 1]) * g0) / (g0 - g1))) / i.rpm]); }
    }
    const rays = { x: [], y: [] }; for (let h = 1; h <= nH; h++) { rays.x.push(0, sp[sp.length - 1], NaN); rays.y.push(0, (h * sp[sp.length - 1]) / 60, NaN); }
    const top = 1.15 * Math.max(F[2][sp.length - 1], Lg[1][sp.length - 1]); rays.y = rays.y.map((v) => (v > top ? top : v)); rays.x = rays.x.map((v, k) => (k % 3 === 1 ? Math.min(v, (top * 60) / (Math.floor(k / 3) + 1)) : v));
    const lagRev = perRev(op.lag.w[0]), warnings = [];
    if (!fin(worst.sep)) worst = { sep: 1, name: 'No mode', h: nH };
    if (100 * worst.sep < i.margin_req) warnings.push(`${worst.name} is within ${(100 * worst.sep).toFixed(1)}% of ${worst.h}/rev at operating speed: resonant amplification of blade loads is expected.`);
    if (lagRev < 1 && lagRev > 0.05) warnings.push('First lag frequency is below 1/rev (soft in-plane): check ground and air resonance.');
    const pure = ctx?.case && !(ctx.case.wing.S_m2 > 0);
    return {
      kpis: [
        { key: 'blade_flap1_rev', label: 'First flap frequency', value: perRev(op.flap.w[0]), unit: '/rev' },
        { key: 'blade_flap2_rev', label: 'Second flap frequency', value: perRev(op.flap.w[1]), unit: '/rev' },
        { key: 'blade_flap3_rev', label: 'Third flap frequency', value: perRev(op.flap.w[2]), unit: '/rev' },
        { key: 'blade_lag1_rev', label: 'First lag frequency', value: lagRev, unit: '/rev', status: Math.abs(lagRev - 1) > 0.1 ? 'ok' : 'warn', note: 'Keep clear of 1/rev' },
        { key: 'blade_f1_Hz', label: 'First flap frequency at operating speed', value: op.flap.w[0] / TAU, unit: 'Hz' },
        { key: 'southwell_K1', label: 'Southwell coefficient, first flap', value: (op.flap.w[0] ** 2 - nr.flap.w[0] ** 2) / Om0 ** 2, unit: '-', note: 'ω² = ω₀² + K·Ω²' },
        { key: 'min_separation_pct', label: 'Smallest separation from a harmonic', value: 100 * worst.sep, unit: '%', status: 100 * worst.sep >= i.margin_req ? 'ok' : 'bad', note: `${worst.name} against ${worst.h}/rev` },
        { key: 'blade_pass_Hz', label: 'Blade-passing frequency', value: nb * rev, unit: 'Hz' },
      ],
      plots: [{ type: 'line', title: 'Fan diagram', xlabel: 'Rotor speed [rpm]', ylabel: 'Frequency [Hz]', series: [...modes.slice(0, 4).map((m) => ({ name: m.name, x: sp, y: m.hist })), { name: `1 to ${nH}/rev`, x: rays.x, y: rays.y, style: 'dash' }], annotations: [{ x: i.rpm, label: 'Operating' }] },
        { type: 'line', title: 'Blade mode shapes at operating speed', xlabel: 'Radius / R [-]', ylabel: 'Normalised deflection [-]', series: [0, 1, 2].map((k) => ({ name: `Flap ${k + 1}`, x: op.flap.x.map((x) => (x + i.e_root * i.R) / i.R), y: unit(op.flap.shapes[k]) })) }],
      tables: [{ title: 'Blade modes at operating speed', columns: ['Mode', 'Non-rotating [Hz]', 'Rotating [Hz]', 'Per rev', 'Southwell K', 'Nearest harmonic', 'Separation [%]'], rows },
        { title: 'Resonance crossings', columns: ['Mode', 'Harmonic', 'Speed [rpm]', '% of operating'], rows: cross.slice(0, 40) }],
      outputs: { ...(pure ? { f1_Hz: op.flap.w[0] / TAU, f2_Hz: op.flap.w[1] / TAU, f3_Hz: op.flap.w[2] / TAU } : {}), lag_freq_rev: lagRev, nearest_harmonic: worst.h },
      warnings, models: ['Finite element beam with centrifugal geometric stiffness', 'Southwell relation', 'Fan (Campbell) diagram'],
      assumptions: ['Uniform blade; flap, lag and torsion uncoupled (no twist, pitch or Coriolis coupling)', 'Lead-lag includes the −Ω² centrifugal softening of in-plane motion', 'No aerodynamic damping or stiffness', 'Default stiffnesses and blade mass are section-shape estimates (typical values), so default frequency placements are indicative only: enter measured blade data before acting on a resonance finding'],
    };
  },
  convergence: { param: 'nEl', label: 'Beam elements', levels: [4, 6, 9, 14, 21], metric: 'blade_flap2_rev' },
  calibration: { params: [{ key: 'EI_flap', min: 1e-3, max: 1e9 }, { key: 'EI_lag', min: 1e-3, max: 1e10 }], sweep: 'rpm', target: 'blade_f1_Hz', note: 'Measured blade frequencies from a whirl-tower or rap test against rotor speed' },
  verify() {
    const b = { R: 6, e_root: 0, root: 'hinged', m_blade: 60, m_tip: 0, EI_flap: 5e4, EI_lag: 2e6 }, Om = 30, r = bladeFreqs(b, Om, 24), c = bladeFreqs({ ...b, root: 'cantilever' }, 0, 24);
    const o = bladeFreqs({ ...b, e_root: 0.06, EI_lag: 1 }, Om, 24);
    return [
      N.check('Centrally hinged blade flaps at exactly 1/rev', r.flap.w[0] / Om, 1, 1e-4, 'Rigid flapping about a central hinge: ω = Ω'),
      N.check('Centrally hinged blade has zero lag frequency', r.lag.w[0] / Om, 0, 2e-2, 'Rigid lead-lag about a central hinge'),
      N.check('Non-rotating cantilever blade first mode', c.flap.w[0], 1.875104 ** 2 * Math.sqrt(5e4 / (10 * 6 ** 4)), 1e-5, 'Euler–Bernoulli exact'),
      N.check('Offset-hinge rigid lag frequency sqrt(3e/2(1−e))', o.lag.w[0] / Om, Math.sqrt((1.5 * 0.06) / 0.94), 5e-3, 'Rigid uniform blade, Johnson, Helicopter Theory'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    // harmonics up to N + 1 per rev carry the significant airloads and feed the N/rev hub loads; higher ones are weak, and a ±5% band around them covers most of the frequency axis
    const strong = o.nearest_harmonic <= Math.round(i.n_blades) + 1;
    if (o.min_separation_pct < i.margin_req) out.push({ severity: strong ? 'warn' : 'advise', title: 'Blade mode too close to an excitation harmonic', detail: res.kpis.find((k) => k.key === 'min_separation_pct').note + `: separation ${o.min_separation_pct.toFixed(1)}% against ${i.margin_req}% required` + (strong ? '.' : '; this harmonic is above N + 1 per rev, where the airload excitation is weak.'), action: 'Retune with tip or spanwise mass, change the stiffness distribution or shift the operating speed; then re-check blade loads and hub vibration in Suite 6.', basis: 'Fan-diagram frequency placement' });
    if (o.blade_lag1_rev < 1) out.push({ severity: 'advise', title: 'Soft in-plane rotor', detail: `First lag mode is ${o.blade_lag1_rev.toFixed(2)}/rev.`, action: 'Run the ground-resonance analysis and size the lag dampers and landing-gear damping.', basis: 'Coleman–Feingold mechanical instability' });
    return out;
  },
};

// ---- 3. forced response, isolation and transient --------------------------------------------
const forced = {
  id: 'forced', title: 'Forced response, transmissibility and isolator design', fidelity: 'analytical',
  summary: 'Harmonic response of an isolated engine or equipment mass on a flexible support: frequency response functions, transmitted force, isolation efficiency and a run-up transient.',
  equations: ['Forced harmonic response equations', 'Linear multi-degree-of-freedom vibration equations', 'Rayleigh damping equations', "Newton's second law"],
  inputs: [
    { key: 'm1', label: 'Isolated mass (engine, gearbox, equipment)', unit: 'kg', default: 2200, min: 1e-4, group: 'System' },
    { key: 'f_iso', label: 'Isolator natural frequency', unit: 'Hz', default: 8, min: 0.2, group: 'System', help: 'Mounted natural frequency on a rigid base; lower gives better isolation but more static deflection' },
    { key: 'zeta_iso', label: 'Isolator damping ratio', unit: '-', default: 0.08, min: 0.001, max: 1, group: 'System', help: '0.05–0.10 elastomer, 0.15–0.3 wire rope or friction-damped' },
    { key: 'm2', label: 'Support modal mass', unit: 'kg', default: 6000, min: 1e-4, group: 'System', help: 'Effective mass of the pylon or airframe mode at the mount' },
    { key: 'f_sup', label: 'Support natural frequency', unit: 'Hz', default: 14, min: 0.2, group: 'System' },
    { key: 'zeta_sup', label: 'Support damping ratio', unit: '-', default: 0.02, min: 0.001, max: 0.5, group: 'System' },
    { key: 'rpm', label: 'Operating speed', unit: 'rpm', default: 4500, min: 1, group: 'Excitation' },
    { key: 'order', label: 'Excitation order', unit: '/rev', default: 1, min: 0.5, max: 60, group: 'Excitation', help: '1 for unbalance, blade count for blade-passing' },
    { key: 'F0', label: 'Force amplitude at operating speed', unit: 'N', default: 2000, min: 0, group: 'Excitation' },
    { key: 'iso_target', label: 'Target isolation', unit: '%', default: 80, min: 10, max: 99, group: 'Limits' },
    { key: 't_run', label: 'Run-up time', unit: 's', default: 8, min: 0.2, max: 120, group: 'Transient' },
    { key: 'nSteps', label: 'Time steps per forcing period', unit: '', default: 16, min: 8, max: 200, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c) => {
    const p = c.prop, m1 = p.type === 'turbofan' || p.type === 'turbojet' ? p.T0_N / 54 : p.type === 'piston' ? p.P0_W / 1300 : p.type === 'electric' ? p.P0_W / 5000 : p.P0_W / 4000, rot = c.rotor.rpm > 0;
    const rpm = rot ? c.rotor.rpm : p.rpm > 0 ? p.rpm : 4500, order = rot ? c.rotor.n_blades : p.rpm > 0 ? p.n_blades || 1 : 1, fex = (rpm * order) / 60;
    return { m1: Math.max(m1, 0.05), m2: 0.08 * c.mass.mtow_kg, rpm, order, f_iso: fex / 3.5, f_sup: fex / 1.8, F0: 0.02 * Math.max(m1, 0.05) * G0, zeta_sup: c.struct.zeta, t_run: N.clamp(600 / fex, 0.5, 20) };
  },
  run(i) {
    const w1 = TAU * i.f_iso, w2 = TAU * i.f_sup, k1 = i.m1 * w1 * w1, c1 = 2 * i.zeta_iso * i.m1 * w1, k2 = i.m2 * w2 * w2, c2 = 2 * i.zeta_sup * i.m2 * w2;
    const M = [[i.m1, 0], [0, i.m2]], K = [[k1, -k1], [-k1, k1 + k2]], Cm = [[c1, -c1], [-c1, c1 + c2]], wf = rads(i.rpm) * i.order;
    const resp = (w) => { const X = harmonic(M, Cm, K, [1, 0], w), d = N.C.sub(X[0], X[1]), ft = N.C.mul([k1, w * c1], d), fg = N.C.mul([k2, w * c2], X[1]); return { x1: N.C.abs(X[0]), x2: N.C.abs(X[1]), tr: N.C.abs(ft), trg: N.C.abs(fg), rel: N.C.abs(d) }; };
    const fs = N.logspace(Math.min(i.f_iso, i.f_sup) / 5, Math.max(wf / TAU, i.f_sup) * 3, 260), R = fs.map((f) => resp(TAU * f)), at = resp(wf), r = wf / w1;
    const trRigid = (rr, z) => Math.sqrt((1 + (2 * z * rr) ** 2) / ((1 - rr * rr) ** 2 + (2 * z * rr) ** 2));
    const eg = N.eigGenSym(K, M), wn = eg.values.map(Math.sqrt), zt = i.zeta_iso, ray = { a: (2 * zt * wn[0] * wn[1]) / (wn[0] + wn[1]), b: (2 * zt) / (wn[0] + wn[1]) };
    // isolator frequency that meets the target on a rigid base: T = target at the forcing frequency
    const tgt = 1 - i.iso_target / 100, fReq = N.findRoot((f) => trRigid(wf / (TAU * f), i.zeta_iso) - tgt, wf / TAU / 60, wf / TAU / 1.45, 80) || NaN;
    // run-up through resonance by Newmark: force ∝ speed², linear speed ramp then dwell
    const dt = TAU / wf / Math.round(i.nSteps), tEnd = i.t_run * 1.25, ns = Math.min(40000, Math.ceil(tEnd / dt)), dtn = tEnd / ns;
    const ph = (t) => (t < i.t_run ? (0.5 * wf * t * t) / i.t_run : wf * (t - i.t_run / 2)), sf = (t) => Math.min(1, t / i.t_run) ** 2;
    const nm = N.newmark(M, Cm, K, (t) => [i.F0 * sf(t) * Math.sin(ph(t)), 0], [0, 0], [0, 0], dtn, ns), rel = nm.u.map((u) => u[0] - u[1]);
    const peakT = N.amax(rel.map(Math.abs)), steady = i.F0 * at.rel, st = (i.m1 * G0) / k1, iso = 1 - at.trg, warnings = [];
    // envelope for plotting
    const blk = Math.max(1, Math.floor(ns / 300)), te = [], env = []; for (let k = 0; k < ns; k += blk) { let m = 0; for (let j = k; j < Math.min(ns, k + blk); j++) m = Math.max(m, Math.abs(rel[j])); te.push(nm.t[k]); env.push(m * 1e3); }
    if (r < Math.SQRT2) warnings.push('The forcing frequency is below √2 times the isolator frequency: the mount amplifies rather than isolates.');
    if (Math.abs(wf / wn[0] - 1) < 0.15 || Math.abs(wf / wn[1] - 1) < 0.15) warnings.push('The forcing frequency is within 15% of a coupled natural frequency.');
    if (st > 0.02) warnings.push(`Static deflection of ${(st * 1e3).toFixed(1)} mm is large: check mount travel, snubbing and engine/driveline alignment.`);
    return {
      kpis: [
        { key: 'transmissibility', label: 'Force transmitted to the airframe / applied', value: at.trg, unit: '-', status: iso * 100 >= i.iso_target ? 'ok' : iso > 0 ? 'warn' : 'bad' },
        { key: 'isolation_pct', label: 'Isolation efficiency', value: 100 * iso, unit: '%' },
        { key: 'tr_rigid_base', label: 'Rigid-base transmissibility', value: trRigid(r, i.zeta_iso), unit: '-' },
        { key: 'freq_ratio', label: 'Forcing / isolator frequency', value: r, unit: '-', status: r > Math.SQRT2 ? 'ok' : 'bad', note: 'Isolation requires a ratio above √2' },
        { key: 'f_coupled1_Hz', label: 'First coupled natural frequency', value: wn[0] / TAU, unit: 'Hz' }, { key: 'f_coupled2_Hz', label: 'Second coupled natural frequency', value: wn[1] / TAU, unit: 'Hz' },
        { key: 'static_defl_mm', label: 'Isolator static deflection (1 g)', value: st * 1e3, unit: 'mm' },
        { key: 'iso_stiffness_Npm', label: 'Isolator stiffness (total)', value: k1, unit: 'N/m' },
        { key: 'f_iso_required_Hz', label: 'Isolator frequency for the target', value: fReq, unit: 'Hz' },
        { key: 'x_steady_mm', label: 'Steady isolator deflection amplitude', value: steady * 1e3, unit: 'mm' },
        { key: 'x_runup_mm', label: 'Peak isolator deflection in run-up', value: peakT * 1e3, unit: 'mm', note: 'Passing through the mount resonance' },
        { key: 'a_sup_g', label: 'Support acceleration amplitude', value: (i.F0 * at.x2 * wf * wf) / G0, unit: 'g' },
        { key: 'rayleigh_alpha', label: 'Rayleigh α for the isolator damping ratio', value: ray.a, unit: '1/s' }, { key: 'rayleigh_beta', label: 'Rayleigh β', value: ray.b, unit: 's' },
      ],
      plots: [
        { type: 'line', title: 'Frequency response functions (receptance)', xlabel: 'Frequency [Hz]', ylabel: 'Displacement per unit force [m/N]', xlog: true, ylog: true, series: [{ name: 'Isolated mass', x: fs, y: R.map((v) => v.x1) }, { name: 'Support', x: fs, y: R.map((v) => v.x2) }], annotations: [{ x: wf / TAU, label: 'Forcing' }] },
        { type: 'line', title: 'Force transmissibility', xlabel: 'Frequency [Hz]', ylabel: 'Transmitted / applied force [-]', xlog: true, ylog: true, series: [{ name: 'Through the isolator', x: fs, y: R.map((v) => v.tr) }, { name: 'Into the airframe', x: fs, y: R.map((v) => v.trg) }, { name: 'Rigid-base theory', x: fs, y: fs.map((f) => trRigid(f / i.f_iso, i.zeta_iso)), style: 'dash' }], annotations: [{ x: wf / TAU, label: 'Forcing' }, { y: 1, label: 'No isolation' }] },
        { type: 'line', title: 'Run-up transient: isolator deflection envelope', xlabel: 'Time [s]', ylabel: 'Deflection amplitude [mm]', series: [{ name: 'Newmark envelope', x: te, y: env }], annotations: [{ y: steady * 1e3, label: 'Steady state' }] },
      ],
      warnings, models: ['Two-degree-of-freedom lumped-mass model', 'Harmonic response by complex solution', 'Newmark average-acceleration transient'],
      assumptions: ['Linear isolators with viscous damping', 'Single translational direction; rocking modes of the mounted mass are not included', 'Excitation force ∝ speed² during the run-up', 'Default masses, mount and support frequencies, damping and force amplitude are class-level estimates (typical values), not data for a specific installation'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps per forcing period', levels: [8, 16, 32, 64], metric: 'x_runup_mm' },
  calibration: { params: [{ key: 'zeta_iso', min: 0.005, max: 0.6 }, { key: 'f_iso', min: 0.5, max: 200 }], sweep: 'rpm', target: 'transmissibility', note: 'Measured transmissibility from a shaker or run-up survey' },
  verify() {
    const b = { m1: 100, f_iso: 10, zeta_iso: 0.1, m2: 1e9, f_sup: 5000, zeta_sup: 0.01, order: 1, F0: 1, iso_target: 80, t_run: 2, nSteps: 20 };
    const a = N.kv(forced.run({ ...b, rpm: 600 * Math.SQRT2 })), k = 100 * (TAU * 10) ** 2, x = harmonic([[100]], [[2 * 0.1 * 100 * TAU * 10]], [[k]], [1], TAU * 10);
    const nm = N.newmark([[1]], [[0]], [[400]], () => [8], [0], [0], 0.002, 400);
    return [
      N.check('Transmissibility is 1 at r = √2 for any damping', a.tr_rigid_base, 1, 1e-9, 'Den Hartog, Mechanical Vibrations'),
      N.check('Two-mass model recovers the rigid-base limit', a.transmissibility, 1, 2e-3, 'Limit of an infinitely heavy, stiff support'),
      N.check('Resonant magnification = 1/(2ζ)', N.C.abs(x[0]) * k, 5, 1e-9, 'Single-degree-of-freedom exact'),
      N.check('Step response overshoots to twice the static value', N.amax(nm.u.map((u) => u[0])), 0.04, 2e-3, 'Undamped SDOF step response'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.isolation_pct < i.iso_target) out.push({ severity: o.isolation_pct < 0 ? 'critical' : 'warn', title: 'Isolation target not met', detail: `Isolation is ${o.isolation_pct.toFixed(0)}% against ${i.iso_target}% required.`, action: fin(o.f_iso_required_Hz) ? `Soften the mounts to about ${o.f_iso_required_Hz.toFixed(1)} Hz, or stiffen the support so that it does not resonate near the forcing frequency.` : 'Soften the mounts or reduce damping; consider a tuned absorber.', basis: 'Transmissibility of a damped isolator' });
    if (o.x_runup_mm > 3 * o.x_steady_mm && o.x_runup_mm > 1) out.push({ severity: 'advise', title: 'Large excursion while passing through the mount resonance', detail: `Run-up peak ${o.x_runup_mm.toFixed(2)} mm against ${o.x_steady_mm.toFixed(2)} mm steady.`, action: 'Accelerate quickly through the resonance, add damping or snubbers. Lower vibration also lengthens fatigue life of brackets and equipment.', basis: 'Transient passage through resonance' });
    return out;
  },
};

// ---- 4. rotor-bearing dynamics: whirl, Campbell diagram, unbalance --------------------------
function rotorSys(i) {
  const L = i.L, a = i.pos * L, b = L - a, EI = (i.E * Math.PI * (i.d_o ** 4 - (i.d_o * i.d_ratio) ** 4)) / 64, m = i.m_disc, Ip = (m * i.r_disc ** 2) / 2, Id = (m * i.r_disc ** 2) / 4 + (m * i.t_disc ** 2) / 12;
  const stiff = (kb) => { // 2×2 flexibility at the disc (deflection, slope) of a simply supported shaft on two equal bearings
    const f11 = (a * a * b * b) / (3 * EI * L) + (a * a + b * b) / (kb * L * L), f12 = (a * b * (b - a)) / (3 * EI * L) + (a - b) / (kb * L * L), f22 = (a * a - a * b + b * b) / (3 * EI * L) + 2 / (kb * L * L), dt = f11 * f22 - f12 * f12;
    return [f22 / dt, -f12 / dt, f11 / dt];
  };
  const kx = stiff(i.kb_x), ky = stiff(i.kb_y);
  // q = [u, v, θ, ψ]: translations along x, y and rotations about x, y; spin about +z
  const K = [[kx[0], 0, 0, kx[1]], [0, ky[0], -ky[1], 0], [0, -ky[1], ky[2], 0], [kx[1], 0, 0, kx[2]]], M = [[m, 0, 0, 0], [0, m, 0, 0], [0, 0, Id, 0], [0, 0, 0, Id]];
  const G = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, Ip], [0, 0, -Ip, 0]];
  const Cd = [2 * i.zeta * Math.sqrt(kx[0] * m), 2 * i.zeta * Math.sqrt(ky[0] * m), 2 * i.zeta * Math.sqrt(ky[2] * Id), 2 * i.zeta * Math.sqrt(kx[2] * Id)];
  return { K, M, G, Cd, kx, ky, m, Ip, Id, a, b, L, EI };
}
/** Undamped whirl modes at spin speed Om: [{w, dir}] ascending in frequency; dir = +1 forward, −1 backward. */
function whirl(s, Om) {
  const A = N.zeros(8);
  for (let p = 0; p < 4; p++) { A[p][p + 4] = 1; for (let q = 0; q < 4; q++) { A[p + 4][q] = -s.K[p][q] / s.M[p][p]; A[p + 4][q + 4] = (-Om * s.G[p][q]) / s.M[p][p]; } }
  const ws = ceig(A).filter((l) => l[1] > 1e-9).map((l) => l[1]).sort((x, y) => x - y);
  return ws.map((w) => {
    const D = s.K.map((r, p) => r.map((k, q) => [k - w * w * 1.0000002 * s.M[p][q], w * Om * s.G[p][q]]));
    let v = [[1, 0], [0.3, 0.2], [0.1, -0.2], [0.2, 0.1]];
    try { for (let it = 0; it < 2; it++) { v = N.csolve(D, v); const nv = Math.sqrt(N.sum(v.map((z) => z[0] * z[0] + z[1] * z[1]))); v = v.map((z) => N.C.scale(z, 1 / nv)); } } catch { /* exactly singular: keep the last iterate */ }
    const sgn = N.C.mul(N.C.conj(v[0]), v[1])[1]; // Im(conj(U)·V) < 0 when y lags x by 90°, i.e. whirl in the spin direction
    return { w, dir: Math.abs(sgn) < 1e-9 ? 0 : sgn < 0 ? 1 : -1 };
  }).map((m, k, all) => (all.some((o, j) => j !== k && Math.abs(o.w - m.w) < 1e-7 * m.w) ? { w: m.w, dir: 0 } : m)); // coincident pair: direction undefined
}
function unbalance(s, Om, e) {
  const D = s.K.map((r, p) => r.map((k, q) => [k - Om * Om * s.M[p][q], Om * Om * s.G[p][q] + (p === q ? Om * s.Cd[p] : 0)]));
  const f = s.m * e * Om * Om, Q = N.csolve(D, [[f, 0], [0, -f], [0, 0], [0, 0]]);
  // semi-axes of the ellipse traced by a rotating vector (x, y) from its forward and backward circular components
  const axes = (x, y) => { const f = N.C.abs(N.C.scale(N.C.add(x, N.C.mul([0, 1], y)), 0.5)), b = N.C.abs(N.C.scale(N.C.sub(x, N.C.mul([0, 1], y)), 0.5)); return { major: f + b, minor: Math.abs(f - b), backward: b > f }; };
  const orb = axes(Q[0], Q[1]);
  const Fx = N.C.add(N.C.scale(Q[0], s.kx[0]), N.C.scale(Q[3], s.kx[1])), My = N.C.add(N.C.scale(Q[0], s.kx[1]), N.C.scale(Q[3], s.kx[2]));
  const Fy = N.C.sub(N.C.scale(Q[1], s.ky[0]), N.C.scale(Q[2], s.ky[1])), Mx = N.C.sub(N.C.scale(Q[2], s.ky[2]), N.C.scale(Q[1], s.ky[1]));
  // bearing reactions rotate with the shaft: the load amplitude is the semi-major axis of the force ellipse at each bearing
  const r1 = axes(N.C.sub(N.C.scale(Fx, s.b / s.L), N.C.scale(My, 1 / s.L)), N.C.add(N.C.scale(Fy, s.b / s.L), N.C.scale(Mx, 1 / s.L))).major;
  const r2 = axes(N.C.add(N.C.scale(Fx, s.a / s.L), N.C.scale(My, 1 / s.L)), N.C.sub(N.C.scale(Fy, s.a / s.L), N.C.scale(Mx, 1 / s.L))).major;
  return { Q, major: orb.major, minor: orb.minor, backward: orb.backward, bearing: Math.max(r1, r2) };
}
const shaft = {
  id: 'shaft', title: 'Rotor-bearing whirl, Campbell diagram and unbalance response', fidelity: 'reduced-order',
  summary: 'A disc on a flexible shaft and flexible bearings with gyroscopic coupling: forward and backward whirl frequencies, critical speeds, unbalance response, orbit and bearing loads.',
  equations: ['Jeffcott rotor equations', 'Rotating shaft equations', 'Gyroscopic equations', 'Campbell frequency relations', 'Forced harmonic response equations'],
  inputs: [
    { key: 'm_disc', label: 'Disc mass (fan, propeller, turbine wheel)', unit: 'kg', default: 240, min: 1e-4, group: 'Rotor' },
    { key: 'r_disc', label: 'Disc radius of gyration ×√2 (equivalent disc radius)', unit: 'm', default: 0.6, min: 1e-3, group: 'Rotor', help: 'Polar inertia = m·r²/2' },
    { key: 't_disc', label: 'Disc axial thickness', unit: 'm', default: 0.15, min: 0, group: 'Rotor' },
    { key: 'L', label: 'Bearing span', unit: 'm', default: 1.5, min: 0.01, group: 'Shaft' },
    { key: 'pos', label: 'Disc position / span', unit: '-', default: 0.35, min: 0.05, max: 0.95, group: 'Shaft', help: '0.5 is the classical Jeffcott rotor with no gyroscopic coupling' },
    { key: 'd_o', label: 'Shaft outer diameter', unit: 'm', default: 0.1, min: 1e-4, group: 'Shaft', help: 'Default is sized for the first rigid-bearing critical speed at 1.4 × operating speed; enter the real shaft' },
    { key: 'd_ratio', label: 'Bore / outer diameter', unit: '-', default: 0.6, min: 0, max: 0.95, group: 'Shaft' },
    { key: 'E', label: 'Shaft Young\'s modulus', unit: 'Pa', default: 200e9, min: 1e9, group: 'Shaft', help: '200 GPa for low-alloy steel (MIL-HDBK-5J), about 71 GPa for aluminium, 110 GPa for titanium' },
    { key: 'kb_x', label: 'Bearing stiffness, horizontal (each)', unit: 'N/m', default: 2e8, min: 1e3, group: 'Bearings' },
    { key: 'kb_y', label: 'Bearing stiffness, vertical (each)', unit: 'N/m', default: 3e8, min: 1e3, group: 'Bearings', help: 'Unequal stiffness splits the critical speeds and makes the orbit elliptical' },
    { key: 'zeta', label: 'Support damping ratio', unit: '-', default: 0.03, min: 1e-4, max: 0.5, group: 'Bearings', help: '0.01–0.03 rolling bearings, 0.05–0.15 with squeeze-film dampers' },
    { key: 'rpm', label: 'Operating speed', unit: 'rpm', default: 4500, min: 1, group: 'Operation' },
    { key: 'rpm_max_frac', label: 'Diagram top speed / operating', unit: '-', default: 2, min: 1.05, max: 6, group: 'Operation' },
    { key: 'G_grade', label: 'Balance quality grade G', unit: 'mm/s', default: 6.3, min: 0.4, max: 630, group: 'Unbalance', help: 'ISO balance quality e·Ω: G2.5 turbines, G6.3 fans and general machinery, G16 propeller shafts' },
    { key: 'sep_req', label: 'Required critical-speed separation', unit: '%', default: 15, min: 1, max: 50, group: 'Limits' },
    { key: 'v_limit', label: 'Vibration velocity alert level', unit: 'mm/s RMS', default: 4.5, min: 0.1, max: 50, group: 'Limits', help: 'Set from the applicable machinery standard or engine manual' },
    { key: 'nSpeeds', label: 'Speed points', unit: '', default: 200, min: 40, max: 2000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c) => {
    const p = c.prop, jet = p.type === 'turbofan' || p.type === 'turbojet', heli = c.meta.type === 'helicopter', D = p.prop_dia_m > 0 ? p.prop_dia_m : c.rotor.R_m > 0 && !heli ? 2 * c.rotor.R_m : 0;
    // representative rotor: fan (jets), tail-rotor drive shaft segment (helicopter), propeller on its shaft (others)
    const m = jet ? 0.002 * p.T0_N : heli ? 6 : Math.max(0.02, 6 * D ** 2.5), r = jet ? 0.07 * Math.sqrt(p.T0_N / 1000) : heli ? 0.06 : Math.max(0.02, 0.25 * D), rpm = jet ? 4500 : heli ? 4100 : p.rpm || c.rotor.rpm || 2500;
    const L = jet ? 2 * r : heli ? 1.4 : Math.max(0.03, 0.35 * D), pos = heli ? 0.5 : 0.35, a = pos * L, b = L - a, ratio = heli ? 0.9 : 0.6, E = heli ? METALS['Al 7075-T6'].E : METALS['Steel 4340 (QT)'].E;
    const k = m * (1.4 * rads(rpm)) ** 2, EI = (k * a * a * b * b) / (3 * L), d = ((64 * EI) / (E * Math.PI * (1 - ratio ** 4))) ** 0.25, kb = 12 * k;
    return { m_disc: m, r_disc: r, t_disc: 0.25 * r, L, pos, d_o: d, d_ratio: ratio, E, kb_x: kb, kb_y: 1.5 * kb, rpm, G_grade: jet || p.type === 'turboshaft' ? 2.5 : 6.3 };
  },
  run(i) {
    const s = rotorSys(i), Om0 = rads(i.rpm), top = Om0 * i.rpm_max_frac, e = (i.G_grade * 1e-3) / Om0, ns = Math.round(i.nSpeeds);
    // Campbell diagram and 1× crossings on each frequency-ordered branch
    const sp = N.linspace(top / 24, top, 24), W = sp.map((o) => whirl(s, o)), fw = { x: [], y: [] }, bw = { x: [], y: [] }, crit = [];
    W.forEach((modes, k) => modes.forEach((m) => { const t = m.dir >= 0 ? fw : bw; t.x.push((sp[k] * 60) / TAU); t.y.push(m.w / TAU); }));
    for (let j = 0; j < 4; j++) for (let k = 1; k < sp.length; k++) {
      const g = (o) => (whirl(s, o)[j]?.w ?? NaN) - o, g0 = W[k - 1][j]?.w - sp[k - 1], g1 = W[k][j]?.w - sp[k];
      if (fin(g0) && fin(g1) && g0 * g1 < 0) { const o = N.brent(g, sp[k - 1], sp[k], 1e-7 * top); crit.push({ rpm: (o * 60) / TAU, dir: whirl(s, o)[j].dir }); }
    }
    const below = W[0].filter((m) => m.w < sp[0]).length; // criticals below the first grid point are not resolved
    crit.sort((p, q) => p.rpm - q.rpm);
    // unbalance response sweep
    const os = N.linspace(top / ns, top, ns), U = os.map((o) => unbalance(s, o, e)), amp = U.map((u) => u.major), peaks = [];
    for (let k = 1; k < ns - 1; k++) if (amp[k] > amp[k - 1] && amp[k] >= amp[k + 1]) { const o = N.goldenSection((x) => -unbalance(s, x, e).major, os[k - 1], os[k + 1], 1e-7); peaks.push({ rpm: (o * 60) / TAU, amp: unbalance(s, o, e).major }); }
    const fwd = crit.filter((c) => c.dir >= 0), nc1 = fwd[0]?.rpm ?? NaN, at = unbalance(s, Om0, e), seps = [...fwd.map((c) => c.rpm), ...peaks.map((p) => p.rpm)].map((r) => Math.abs(r / i.rpm - 1)); // backward crossings count only where they show as a response peak
    const sep = seps.length ? 100 * N.amin(seps) : 100 * (i.rpm_max_frac - 1), vel = (at.major * Om0 * 1e3) / Math.SQRT2, jeff = (Math.sqrt(s.kx[0] / s.m) * 60) / TAU;
    const th = N.linspace(0, TAU, 73), orbit = (Q) => ({ x: th.map((t) => 1e6 * (Q[0][0] * Math.cos(t) - Q[0][1] * Math.sin(t))), y: th.map((t) => 1e6 * (Q[1][0] * Math.cos(t) - Q[1][1] * Math.sin(t))) });
    const pk = peaks[0] ? unbalance(s, rads(peaks[0].rpm), e) : at, warnings = [];
    if (!fin(nc1)) warnings.push('No forward critical speed was found inside the diagram range; raise the top speed to locate it.');
    if (below) warnings.push('A whirl mode lies below the lowest speed analysed; very soft supports are outside the resolved range.');
    if (sep < i.sep_req) warnings.push(`Operating speed is within ${sep.toFixed(1)}% of a critical speed.`);
    if (fin(nc1) && i.rpm > nc1) warnings.push('The rotor runs super-critically: it must pass through a critical speed on every start and stop, so support damping is essential.');
    if (i.L / i.d_o < 4) warnings.push('Short, stiff shaft: shear deformation and bearing rotational stiffness, both neglected, become significant.');
    return {
      kpis: [
        { key: 'crit_speed_rpm', label: 'First forward critical speed', value: nc1, unit: 'rpm' },
        { key: 'crit_speed2_rpm', label: 'Second forward critical speed', value: fwd[1]?.rpm ?? NaN, unit: 'rpm' },
        { key: 'crit_back_rpm', label: 'First backward 1× crossing', value: crit.find((c) => c.dir < 0)?.rpm ?? NaN, unit: 'rpm', note: 'Excited by unbalance only when the supports are anisotropic' },
        { key: 'jeffcott_rpm', label: 'Non-gyroscopic estimate sqrt(k/m)', value: jeff, unit: 'rpm' },
        { key: 'crit_separation_pct', label: 'Separation of operating speed from the nearest critical', value: sep, unit: '%', status: sep >= i.sep_req ? 'ok' : 'bad', note: `Requirement ${i.sep_req}%` },
        { key: 'e_perm_um', label: 'Permissible mass eccentricity', value: e * 1e6, unit: 'µm', note: `Grade G${i.G_grade} at operating speed` },
        { key: 'U_perm_gmm', label: 'Permissible residual unbalance', value: e * s.m * 1e6, unit: 'g·mm' },
        { key: 'orbit_op_um', label: 'Orbit semi-major axis at operating speed', value: at.major * 1e6, unit: 'µm' },
        { key: 'orbit_peak_um', label: 'Largest orbit passing the first critical', value: (peaks[0]?.amp ?? NaN) * 1e6, unit: 'µm' },
        { key: 'vib_velocity_mms', label: 'Vibration velocity at operating speed', value: vel, unit: 'mm/s RMS', status: vel <= i.v_limit ? 'ok' : 'warn' },
        { key: 'bearing_load_N', label: 'Dynamic bearing load at operating speed', value: at.bearing, unit: 'N', note: 'Amplitude of the rotating reaction at the more heavily loaded bearing' },
        { key: 'unbalance_force_N', label: 'Unbalance force at operating speed', value: s.m * e * Om0 * Om0, unit: 'N' },
      ].filter((k) => fin(k.value) || k.key === 'crit_speed_rpm'),
      plots: [
        { type: 'line', title: 'Campbell diagram', xlabel: 'Rotor speed [rpm]', ylabel: 'Whirl frequency [Hz]', series: [{ name: 'Forward whirl', x: fw.x, y: fw.y, style: 'points' }, { name: 'Backward whirl', x: bw.x, y: bw.y, style: 'points' }, { name: '1× (synchronous)', x: [0, (top * 60) / TAU], y: [0, top / TAU], style: 'dash' }], annotations: [{ x: i.rpm, label: 'Operating' }] },
        { type: 'line', title: 'Unbalance response', xlabel: 'Rotor speed [rpm]', ylabel: 'Orbit semi-axis [µm]', ylog: true, series: [{ name: 'Semi-major axis', x: thin(os.map((o) => (o * 60) / TAU)), y: thin(amp.map((v) => v * 1e6)) }, { name: 'Semi-minor axis', x: thin(os.map((o) => (o * 60) / TAU)), y: thin(U.map((u) => Math.max(u.minor * 1e6, 1e-6))) }], annotations: [{ x: i.rpm, label: 'Operating' }] },
        { type: 'line', title: 'Disc orbit', xlabel: 'Horizontal displacement [µm]', ylabel: 'Vertical displacement [µm]', equalAspect: true, series: [{ name: 'At operating speed', ...orbit(at.Q) }, { name: 'At the first response peak', ...orbit(pk.Q) }] },
        { type: 'line', title: 'Dynamic bearing load', xlabel: 'Rotor speed [rpm]', ylabel: 'Bearing load amplitude [N]', ylog: true, series: [{ name: 'More heavily loaded bearing', x: thin(os.map((o) => (o * 60) / TAU)), y: thin(U.map((u) => Math.max(u.bearing, 1e-9))) }] },
      ],
      tables: [{ title: 'Synchronous (1×) crossings and response peaks', columns: ['Kind', 'Speed [rpm]', '% of operating'], rows: [...crit.map((c) => [c.dir >= 0 ? 'Forward critical' : 'Backward crossing', c.rpm, (100 * c.rpm) / i.rpm]), ...peaks.map((p) => ['Unbalance response peak', p.rpm, (100 * p.rpm) / i.rpm])] }],
      warnings, models: ['Jeffcott rotor generalised to an off-centre disc with gyroscopic moments (4 degrees of freedom)', 'Campbell diagram model', 'Rotor-bearing model with anisotropic support stiffness', 'ISO balance quality grade G = e·Ω'],
      assumptions: ['Rigid disc on a massless, uniform, simply supported shaft; two identical bearings', 'Bearing stiffness is constant with speed; no cross-coupling, oil-film or seal forces', 'Viscous damping applied at the disc as an equivalent support damping ratio', 'Default rotor dimensions are representative for the vehicle class, not a specific engine'],
    };
  },
  convergence: { param: 'nSpeeds', label: 'Speed points in the response sweep', levels: [50, 100, 200, 400], metric: 'orbit_peak_um' },
  calibration: { params: [{ key: 'kb_x', min: 1e4, max: 1e11 }, { key: 'kb_y', min: 1e4, max: 1e11 }, { key: 'zeta', min: 1e-3, max: 0.4 }], sweep: 'rpm', target: 'orbit_op_um', note: 'Measured shaft orbit or casing vibration against speed from a run-up or balancing run' },
  verify() {
    const b = { m_disc: 20, r_disc: 0.2, t_disc: 0, L: 1, d_o: 0.04, d_ratio: 0, E: 2e11, kb_x: 1e14, kb_y: 1e14, zeta: 1e-4, rpm: 3000, rpm_max_frac: 2, G_grade: 6.3, sep_req: 15, v_limit: 4.5, nSpeeds: 200 };
    const j = N.kv(shaft.run({ ...b, pos: 0.5 })), kJ = (48 * 2e11 * Math.PI * 0.04 ** 4) / 64;
    // off-centre thin disc: forward whirl must satisfy the complex-coordinate characteristic equation
    const s = rotorSys({ ...b, pos: 0.3 }), Om = 250, m = whirl(s, Om), f = m.filter((x) => x.dir > 0), det = (w) => (s.kx[0] - s.m * w * w) * (s.kx[2] - s.Id * w * w + s.Ip * Om * w) - s.kx[1] ** 2;
    const sc = s.kx[0] * s.kx[2], o = N.kv(shaft.run({ ...b, pos: 0.3 })), cr = Math.sqrt(N.brent((x) => (s.kx[0] - s.m * x) * (s.kx[2] + (s.Ip - s.Id) * x) - s.kx[1] ** 2, 1, s.kx[0] / s.m * 1.5));
    return [
      N.check('Jeffcott critical speed sqrt(48EI/L³m)', j.crit_speed_rpm, (Math.sqrt(kJ / 20) * 60) / TAU, 1e-5, 'Jeffcott (1919)'),
      N.check('Forward whirl satisfies the gyroscopic characteristic equation', det(f[0].w) / sc, 0, 1e-6, 'det[k11−mω², k12; k12, k22−Id·ω²+Ip·Ω·ω] = 0'),
      N.check('Two forward and two backward whirl modes', f.length + 10 * m.filter((x) => x.dir < 0).length, 22, 1e-9, 'Isotropic rotor with one disc'),
      N.check('Forward critical with gyroscopic stiffening', o.crit_speed_rpm, (cr * 60) / TAU, 1e-5, 'Synchronous whirl: effective inertia Id − Ip'),
      N.check('Mid-span disc on rigid bearings: each bearing carries half the transmitted force m·e·Ω²/(1 − r²)', (2 * unbalance(rotorSys({ ...b, pos: 0.5 }), 30, 1e-5).bearing) / (20 * 1e-5 * 900), 1 / (1 - 900 / (kJ / 20)), 1e-6, 'Jeffcott rotor below the critical speed'),
      N.check('Unbalance orbit tends to e at high speed (self-centring)', unbalance(rotorSys({ ...b, pos: 0.5 }), 6000, 1e-5).major / 1e-5, 1 / (1 - (kJ / 20) / 6000 ** 2), 1e-4, 'Jeffcott response r²/(r²−1)'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.crit_separation_pct < i.sep_req) out.push({ severity: 'critical', title: 'Operating speed is too close to a critical speed', detail: `Separation is ${o.crit_separation_pct.toFixed(1)}% against ${i.sep_req}% required; first forward critical at ${fin(o.crit_speed_rpm) ? o.crit_speed_rpm.toFixed(0) : '—'} rpm.`, action: 'Change shaft stiffness or bearing span, move the disc, or retune the bearing support stiffness; add a squeeze-film damper if the rotor must run near or through the critical.', basis: 'Critical-speed separation margin' });
    if (o.vib_velocity_mms > i.v_limit) out.push({ severity: 'warn', title: 'Vibration at operating speed exceeds the alert level', detail: `${o.vib_velocity_mms.toFixed(2)} mm/s RMS for a G${i.G_grade} balance grade.`, action: `Balance to a finer grade: residual unbalance must fall below ${(o.U_perm_gmm * i.v_limit / o.vib_velocity_mms).toFixed(0)} g·mm. Lower vibration also extends bearing life.`, basis: 'Unbalance response at the stated balance quality grade' });
    if (fin(o.crit_speed_rpm) && i.rpm > o.crit_speed_rpm) out.push({ severity: 'advise', title: 'Super-critical rotor', detail: `Peak orbit while passing the critical is ${fin(o.orbit_peak_um) ? o.orbit_peak_um.toFixed(0) : '—'} µm.`, action: 'Check tip and seal clearances against the transit orbit and specify the support damping; avoid dwelling near the critical speed.', basis: 'Unbalance response through resonance' });
    return out;
  },
};

// ---- 5. drivetrain torsional vibration ------------------------------------------------------
const holzer = (J, k, w) => { let th = 1, T = 0; for (let j = 0; j < J.length; j++) { T += J[j] * w * w * th; if (j < k.length) th -= T / k[j]; } return T; };
const torsional = {
  id: 'torsional', title: 'Drivetrain torsional vibration', fidelity: 'reduced-order',
  summary: 'Engine–gearbox–rotor (or propeller) torsional natural frequencies and mode shapes, interference with engine and blade-passing orders, and the dynamic torque in the shafts.',
  equations: ['Torsional vibration equations', 'Eigenvalue equations', 'Forced harmonic response equations', 'Campbell frequency relations'],
  inputs: [
    { key: 'J_eng', label: 'Engine / motor rotor inertia', unit: 'kg·m²', default: 0.12, min: 1e-9, group: 'Inertias', help: 'At engine speed; all engines combined' },
    { key: 'J_gbx', label: 'Gearbox inertia (at output speed)', unit: 'kg·m²', default: 40, min: 1e-9, group: 'Inertias' },
    { key: 'J_rot', label: 'Rotor / propeller / fan inertia', unit: 'kg·m²', default: 10000, min: 1e-9, group: 'Inertias' },
    { key: 'ratio', label: 'Reduction ratio (engine / output speed)', unit: '-', default: 80, min: 0.05, group: 'Drive' },
    { key: 'k_in', label: 'Input shaft stiffness (at engine speed)', unit: 'N·m/rad', default: 2e4, min: 1e-6, group: 'Drive' },
    { key: 'k_out', label: 'Output shaft stiffness', unit: 'N·m/rad', default: 3.5e6, min: 1e-6, group: 'Drive', help: 'Default sizes the shaft for the rated torque at the stated shear stress' },
    { key: 'zeta', label: 'Modal damping ratio', unit: '-', default: 0.02, min: 1e-4, max: 0.5, group: 'Drive', help: '0.01–0.03 for steel shafting; more with elastomeric couplings or lag dampers' },
    { key: 'rpm', label: 'Output operating speed', unit: 'rpm', default: 258, min: 1, group: 'Operation' },
    { key: 'order_out', label: 'Output-side excitation order', unit: '/rev', default: 4, min: 0.5, max: 60, group: 'Excitation', help: 'Blade count for rotors and propellers' },
    { key: 'order_eng', label: 'Engine-side excitation order', unit: '/rev of engine', default: 1, min: 0.5, max: 12, group: 'Excitation', help: '1 for turbines and motors; half the cylinder count for four-stroke pistons' },
    { key: 'T_osc', label: 'Oscillatory torque at the rotor', unit: 'N·m', default: 5000, min: 0, group: 'Excitation' },
    { key: 'sep_req', label: 'Required frequency separation', unit: '%', default: 15, min: 1, max: 50, group: 'Limits' },
  ],
  defaults: (c) => {
    const p = c.prop, r = c.rotor, heli = c.meta.type === 'helicopter', jet = p.type === 'turbofan' || p.type === 'turbojet', P = Math.max(1, p.P0_W * p.n_eng);
    if (jet) { const rf = 0.07 * Math.sqrt(p.T0_N / 1000), Jf = 0.5 * 0.002 * p.T0_N * rf * rf, Pf = (p.T0_N * 150) / 0.8, Q = Pf / rads(4500), d = ((16 * Q) / (Math.PI * 1.2e8)) ** (1 / 3); return { J_eng: 0.35 * Jf, J_gbx: 0.01 * Jf, J_rot: Jf, ratio: 1, rpm: 4500, order_out: 22, order_eng: 1, k_out: (8e10 * Math.PI * d ** 4) / 32 / (4 * rf), k_in: 50 * (8e10 * Math.PI * d ** 4) / 32 / (4 * rf), T_osc: 0.02 * Q }; }
    const b = bladeDefaults(c), nb = b.n_blades || 2, Jb = (nb * (b.m_blade || 1) * (b.R || 1) ** 2) / 3, nRot = heli ? 1 : p.n_eng, rpm = b.rpm || 2000, Q = (heli ? P : P / nRot) / rads(rpm), d = ((16 * Q) / (Math.PI * 1.2e8)) ** (1 / 3), Ls = Math.max(0.05, 0.15 * (b.R || 1));
    const ratio = p.type === 'turboprop' || p.type === 'turboshaft' ? Math.max(1, 20000 / rpm) : 1, Pe = heli ? P : P / nRot, Je = p.type === 'piston' ? 2e-6 * Pe : p.type === 'electric' ? 2.5e-8 * Pe ** 1.2 : 4e-8 * Pe, ko = (8e10 * Math.PI * d ** 4) / 32 / Ls;
    return { J_eng: Je, J_gbx: 0.004 * Jb + 1e-9, J_rot: Jb, ratio, rpm, order_out: nb, order_eng: p.type === 'piston' ? 2 : 1, k_out: ko, k_in: (3 * ko) / ratio ** 2, T_osc: 0.05 * Q };
  },
  run(i) {
    const n2 = i.ratio ** 2, J = [i.J_eng * n2, i.J_gbx, i.J_rot], k = [i.k_in * n2, i.k_out]; // referred to the output shaft
    const K = [[k[0], -k[0], 0], [-k[0], k[0] + k[1], -k[1]], [0, -k[1], k[1]]], M = [[J[0], 0, 0], [0, J[1], 0], [0, 0, J[2]]], e = geig(K, M), wn = e.values.map((v) => Math.sqrt(Math.max(v, 0)));
    const Om = rads(i.rpm), ex = [{ name: `${i.order_out}/rev of the output`, w: i.order_out * Om, ord: i.order_out }, { name: `${i.order_eng}/rev of the engine`, w: i.order_eng * i.ratio * Om, ord: i.order_eng * i.ratio }];
    const rows = []; let sep = Infinity;
    for (const m of [1, 2]) for (const x of ex) { const s = Math.abs(wn[m] / x.w - 1); sep = Math.min(sep, s); rows.push([`Mode ${m}`, wn[m] / TAU, x.name, x.w / TAU, 100 * s, (wn[m] / x.ord) * (60 / TAU)]); }
    // steady response to the rotor oscillatory torque with stiffness-proportional damping
    const beta = (2 * i.zeta) / (wn[1] || 1), Cm = N.mscale(K, beta), torque = (w) => { const X = harmonic(M, Cm, K, [0, 0, i.T_osc], w); return [N.C.abs(N.C.sub(X[0], X[1])) * k[0], N.C.abs(N.C.sub(X[1], X[2])) * k[1]]; };
    const tq = torque(ex[0].w), fs = N.linspace(wn[1] / TAU / 5, Math.max(wn[2] / TAU, ex[0].w / TAU) * 1.4, 300), TQ = fs.map((f) => torque(TAU * f));
    const hw = N.linspace(0.02 * wn[1], 1.25 * wn[2], 240), norm = J[0] * wn[2] ** 2, warnings = [];
    if (100 * sep < i.sep_req) warnings.push(`A torsional mode is within ${(100 * sep).toFixed(1)}% of an excitation order at operating speed.`);
    if (wn[1] / TAU < 8 && i.J_rot > 100) warnings.push('The first torsional mode is low enough to interact with the engine fuel-control (governor) loop: check torsional stability with the engine controller.');
    const sh = (v) => { const m = N.amax(v.map(Math.abs)) || 1; return v.map((x) => x / m); };
    return {
      kpis: [
        { key: 'f_tors1_Hz', label: 'First torsional frequency', value: wn[1] / TAU, unit: 'Hz' }, { key: 'f_tors2_Hz', label: 'Second torsional frequency', value: wn[2] / TAU, unit: 'Hz' },
        { key: 'tors_separation_pct', label: 'Smallest separation from an excitation order', value: 100 * sep, unit: '%', status: 100 * sep >= i.sep_req ? 'ok' : 'bad' },
        { key: 'tors_resonant_rpm', label: 'Output speed for mode 1 at the output order', value: (wn[1] / i.order_out) * (60 / TAU), unit: 'rpm' },
        { key: 'torque_out_Nm', label: 'Dynamic torque in the output shaft', value: tq[1], unit: 'N·m' },
        { key: 'torque_in_Nm', label: 'Dynamic torque in the input shaft (at engine speed)', value: tq[0] / i.ratio, unit: 'N·m' },
        { key: 'torque_magnifier', label: 'Output-shaft dynamic magnifier', value: i.T_osc ? tq[1] / i.T_osc : 0, unit: '-', status: !i.T_osc || tq[1] / i.T_osc < 2 ? 'ok' : 'warn' },
        { key: 'J_referred_ratio', label: 'Referred engine inertia / rotor inertia', value: J[0] / J[2], unit: '-' },
      ],
      plots: [
        { type: 'line', title: 'Torsional mode shapes (referred to the output shaft)', xlabel: 'Station (1 engine, 2 gearbox, 3 rotor)', ylabel: 'Normalised twist [-]', series: [1, 2].map((m) => ({ name: `Mode ${m}: ${(wn[m] / TAU).toPrecision(3)} Hz`, x: [1, 2, 3], y: sh(e.vectors[m]), style: 'line+points' })) },
        { type: 'line', title: 'Shaft dynamic torque per unit rotor torque', xlabel: 'Excitation frequency [Hz]', ylabel: 'Torque magnifier [-]', ylog: true, series: [{ name: 'Output shaft', x: fs, y: TQ.map((t) => Math.max(t[1] / (i.T_osc || 1), 1e-9)) }, { name: 'Input shaft (referred)', x: fs, y: TQ.map((t) => Math.max(t[0] / (i.T_osc || 1), 1e-9)) }], annotations: ex.map((x) => ({ x: x.w / TAU, label: x.name })) },
        { type: 'line', title: 'Holzer residual torque', xlabel: 'Trial frequency [Hz]', ylabel: 'Residual torque / J₁ω₂² [-]', series: [{ name: 'Residual', x: hw.map((w) => w / TAU), y: hw.map((w) => holzer(J, k, w) / norm) }], annotations: [{ y: 0, label: 'Natural frequencies at the zero crossings' }] },
      ],
      tables: [{ title: 'Interference with excitation orders', columns: ['Mode', 'Frequency [Hz]', 'Excitation', 'Excitation [Hz]', 'Separation [%]', 'Resonant output speed [rpm]'], rows }],
      warnings, models: ['Three-inertia lumped torsional model with gear ratio referral', 'Matrix eigen-solution, cross-checked by the Holzer method', 'Harmonic response with stiffness-proportional damping'],
      assumptions: ['Rigid gears and blades; blade lead-lag flexibility (which lowers the first mode of helicopter drivetrains) is not included', 'Linear shafts and couplings, no backlash', 'Default inertias and stiffnesses are sizing-rule estimates'],
    };
  },
  calibration: { params: [{ key: 'k_out', min: 1, max: 1e10 }, { key: 'zeta', min: 1e-3, max: 0.3 }], sweep: 'J_rot', target: 'f_tors1_Hz', note: 'Measured torsional frequencies from a strain-gauged shaft survey' },
  verify() {
    const r = N.kv(torsional.run({ J_eng: 2, J_gbx: 1e-9, J_rot: 5, ratio: 1, k_in: 1e12, k_out: 4000, zeta: 0.02, rpm: 1000, order_out: 2, order_eng: 1, T_osc: 0, sep_req: 15 })), w = Math.sqrt((4000 * 7) / 10);
    const J = [3, 2, 6], k = [5000, 9000], e = geig([[5000, -5000, 0], [-5000, 14000, -9000], [0, -9000, 9000]], [[3, 0, 0], [0, 2, 0], [0, 0, 6]]);
    return [
      N.check('Two-inertia free–free frequency sqrt(k(J₁+J₂)/J₁J₂)', r.f_tors1_Hz * TAU, w, 1e-6, 'Exact'),
      N.check('Holzer residual vanishes at the matrix eigenfrequency', holzer(J, k, Math.sqrt(e.values[1])) / (3 * e.values[1]), 0, 1e-8, 'Holzer (1921)'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.tors_separation_pct < i.sep_req) out.push({ severity: 'warn', title: 'Torsional resonance near operating speed', detail: `Separation ${o.tors_separation_pct.toFixed(1)}% against ${i.sep_req}% required.`, action: 'Retune with a softer or stiffer coupling, change shaft diameter, or add a torsional damper; avoid continuous operation at the resonant speed.', basis: 'Torsional interference diagram' });
    if (o.torque_magnifier > 2) out.push({ severity: 'advise', title: 'High dynamic torque in the output shaft', detail: `Magnifier ${o.torque_magnifier.toFixed(1)} on the applied oscillatory torque.`, action: 'Carry the dynamic torque into the shaft and gear fatigue assessment (Suites 9 and 18).', basis: 'Forced torsional response' });
    return out;
  },
};

// ---- 6. random vibration --------------------------------------------------------------------
const psdProfile = (i) => (f) => (f < i.f_lo || f > i.f_hi ? 0 : f < i.f1 ? i.W0 * (f / i.f1) ** (i.slope_up / 3.0103) : f <= i.f2 ? i.W0 : i.W0 * (f / i.f2) ** (-i.slope_dn / 3.0103));
const random = {
  id: 'random', title: 'Random vibration response of mounted equipment', fidelity: 'analytical',
  summary: 'Response of a mounted item to a broadband acceleration spectrum: response spectrum, RMS and 3σ acceleration, relative displacement and Miles\' estimate.',
  equations: ['Linear multi-degree-of-freedom vibration equations', 'Forced harmonic response equations', "Newton's second law"],
  inputs: [
    { key: 'fn', label: 'Equipment natural frequency', unit: 'Hz', default: 120, min: 1, group: 'Equipment' },
    { key: 'zeta', label: 'Damping ratio', unit: '-', default: 0.05, min: 0.002, max: 0.7, group: 'Equipment', help: 'Q = 1/(2ζ); Q of 10 (ζ = 0.05) is a common assumption' },
    { key: 'mass', label: 'Equipment mass', unit: 'kg', default: 5, min: 1e-4, group: 'Equipment' },
    { key: 'W0', label: 'Plateau acceleration spectral density', unit: 'g²/Hz', default: 0.04, min: 1e-8, group: 'Input spectrum', help: 'Take the level and shape from the applicable environmental test standard for the equipment zone' },
    { key: 'f_lo', label: 'Lowest frequency', unit: 'Hz', default: 10, min: 0.1, group: 'Input spectrum' },
    { key: 'f1', label: 'Plateau start', unit: 'Hz', default: 40, min: 0.1, group: 'Input spectrum' },
    { key: 'f2', label: 'Plateau end', unit: 'Hz', default: 500, min: 0.2, group: 'Input spectrum' },
    { key: 'f_hi', label: 'Highest frequency', unit: 'Hz', default: 2000, min: 1, group: 'Input spectrum' },
    { key: 'slope_up', label: 'Rising slope', unit: 'dB/oct', default: 6, min: 0, max: 24, group: 'Input spectrum' },
    { key: 'slope_dn', label: 'Falling slope', unit: 'dB/oct', default: 6, min: 0, max: 24, group: 'Input spectrum' },
    { key: 'hours', label: 'Exposure time', unit: 'h', default: 1, min: 0.01, group: 'Input spectrum' },
    { key: 'nFreq', label: 'Frequency points', unit: '', default: 2000, min: 100, max: 40000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c) => ({ zeta: Math.max(0.02, 2.5 * c.struct.zeta) }),
  run(i) {
    const W = psdProfile(i), n = Math.round(i.nFreq), lo = Math.min(i.f_lo, i.fn / 4), f = N.logspace(lo, i.f_hi, n), z = i.zeta;
    const H2 = (r) => (1 + (2 * z * r) ** 2) / ((1 - r * r) ** 2 + (2 * z * r) ** 2), Win = f.map(W), Wout = f.map((x, k) => Win[k] * H2(x / i.fn));
    const rel = f.map((x, k) => (Win[k] * G0 * G0) / ((TAU * i.fn) ** 4 * ((1 - (x / i.fn) ** 2) ** 2 + (2 * z * x / i.fn) ** 2))); // relative displacement PSD [m²/Hz]
    const gin = Math.sqrt(N.trapz(f, Win)), gout = Math.sqrt(N.trapz(f, Wout)), xr = Math.sqrt(N.trapz(f, rel)), miles = Math.sqrt((Math.PI / 2) * i.fn * (1 / (2 * z)) * W(i.fn));
    const nu0 = Math.sqrt(N.trapz(f, f.map((x, k) => x * x * Wout[k])) / (gout * gout || 1)), cycles = nu0 * i.hours * 3600, peak = Math.sqrt(2 * Math.log(Math.max(cycles, 2)));
    const warnings = [];
    if (!(W(i.fn) > 0)) warnings.push('The natural frequency lies outside the excitation band: Miles\' equation gives zero and the response is quasi-static.');
    if (i.f1 < i.f_lo || i.f2 > i.f_hi || i.f2 < i.f1) warnings.push('Spectrum break frequencies are out of order; check the input spectrum.');
    if (z > 0.2) warnings.push('Miles\' equation assumes light damping; with ζ above 0.2 use the integrated value.');
    return {
      kpis: [
        { key: 'grms_in', label: 'Input acceleration', value: gin, unit: 'g RMS' },
        { key: 'grms_out', label: 'Response acceleration (integrated)', value: gout, unit: 'g RMS' },
        { key: 'grms_miles', label: 'Response acceleration (Miles)', value: miles, unit: 'g RMS', note: 'sqrt(π/2 · fn · Q · W(fn)); exact only for a flat spectrum and light damping' },
        { key: 'g_3sigma', label: '3σ response acceleration', value: 3 * gout, unit: 'g', note: 'Exceeded 0.27% of the time for a Gaussian response' },
        { key: 'load_3sigma_N', label: '3σ inertial load on the mounts', value: 3 * gout * G0 * i.mass, unit: 'N' },
        { key: 'x_rel_3sigma_mm', label: '3σ relative displacement (sway space)', value: 3 * xr * 1e3, unit: 'mm' },
        { key: 'amplification', label: 'RMS amplification', value: gout / gin, unit: '-' },
        { key: 'nu0_Hz', label: 'Apparent response frequency', value: nu0, unit: 'Hz', note: 'Positive zero-crossing rate' },
        { key: 'cycles', label: 'Response cycles in the exposure', value: cycles, unit: '-' },
        { key: 'peak_factor', label: 'Expected largest peak / RMS', value: peak, unit: '-', note: 'Rayleigh peaks: sqrt(2·ln N)' },
      ],
      plots: [
        { type: 'line', title: 'Acceleration spectral density', xlabel: 'Frequency [Hz]', ylabel: 'ASD [g²/Hz]', xlog: true, ylog: true, series: [{ name: 'Input', x: thin(f), y: thin(Win.map((v) => Math.max(v, 1e-12))) }, { name: 'Response', x: thin(f), y: thin(Wout.map((v) => Math.max(v, 1e-12))) }], annotations: [{ x: i.fn, label: 'fn' }] },
        { type: 'line', title: 'Cumulative response RMS', xlabel: 'Frequency [Hz]', ylabel: 'Cumulative acceleration [g RMS]', xlog: true, series: [{ name: 'Response', x: thin(f), y: thin(N.cumtrapz(f, Wout).map(Math.sqrt)) }, { name: 'Input', x: thin(f), y: thin(N.cumtrapz(f, Win).map(Math.sqrt)) }] },
      ],
      warnings, models: ['Single-degree-of-freedom base-excitation transmissibility', "Miles' equation", 'Gaussian narrow-band response with Rayleigh peaks'],
      assumptions: ['Stationary Gaussian excitation', 'One dominant mode of the mounted item', 'Linear response: no rattling, snubbing or mount non-linearity', 'The default input spectrum is illustrative: take the level and shape from the applicable environmental test standard for the equipment zone'],
    };
  },
  convergence: { param: 'nFreq', label: 'Frequency points', levels: [250, 500, 1000, 2000, 4000], metric: 'grms_out' },
  calibration: { params: [{ key: 'zeta', min: 0.005, max: 0.5 }, { key: 'fn', min: 1, max: 5000 }], sweep: 'W0', target: 'grms_out', note: 'Measured response g RMS from a random vibration test at several input levels' },
  verify() {
    const b = { fn: 100, zeta: 0.02, mass: 1, W0: 0.01, f_lo: 0.01, f1: 0.01, f2: 1e5, f_hi: 1e5, slope_up: 0, slope_dn: 0, hours: 1, nFreq: 40000 }, r = N.kv(random.run(b)), ex = Math.sqrt(((Math.PI * 100) / (4 * 0.02)) * (1 + 4 * 0.02 ** 2) * 0.01);
    return [
      N.check('White-noise response equals the exact integral πfn(1+4ζ²)/(4ζ)·W', r.grms_out, ex, 2e-3, 'Crandall & Mark, Random Vibration in Mechanical Systems'),
      N.check("Miles' equation within 4ζ² of the exact white-noise result", r.grms_miles, ex, 1e-3, 'Miles (1954)'),
      N.check('Positive zero-crossing rate tends to fn for light damping', N.kv(random.run({ ...b, f_hi: 400, f2: 400 })).nu0_Hz, 100, 2e-2, 'Narrow-band process'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.amplification > 3) out.push({ severity: 'advise', title: 'Mounted resonance amplifies the environment', detail: `Response is ${o.amplification.toFixed(1)} times the input RMS; 3σ load ${o.load_3sigma_N.toFixed(0)} N.`, action: 'Move the natural frequency out of the plateau, add damping (Q of 5 halves the response relative to Q of 20), or isolate. Size brackets and fasteners for the 3σ load and check fatigue in Suite 9.', basis: "Miles' equation" });
    out.push({ severity: 'info', title: 'Reserve sway space', detail: `3σ relative displacement is ${o.x_rel_3sigma_mm.toFixed(2)} mm.`, action: 'Keep at least this clearance to neighbouring structure and harnesses to avoid impact and chafing.', basis: 'Relative-displacement response spectrum' });
    return out;
  },
};

// ---- 7. helicopter ground resonance ---------------------------------------------------------
function grModel(i) {
  const eR = i.e_lag * i.R, l = i.R - eR, S = (i.m_blade * l) / 2, I = (i.m_blade * l * l) / 3, nb = Math.round(i.n_blades), h = nb / 2;
  const wx = TAU * i.f_x, wy = TAU * i.f_y, Mx = i.M_x + nb * i.m_blade, My = i.M_y + nb * i.m_blade, cx = 2 * i.zeta_x * wx * Mx, cy = 2 * i.zeta_y * wy * My;
  const wz = (Om) => Math.sqrt((i.k_lag + eR * S * Om * Om) / I);
  /** Eigenvalues of hub (x, y) + cyclic lag (ζ1c, ζ1s) in the non-rotating frame. */
  const eig = (Om, cz, coup = 1) => {
    const s = coup * h * S, wz2 = wz(Om) ** 2;
    const M = [[Mx, 0, 0, -s], [0, My, s, 0], [0, s, h * I, 0], [-s, 0, 0, h * I]], Cm = [[cx, 0, 0, 0], [0, cy, 0, 0], [0, 0, h * cz, 2 * h * I * Om], [0, 0, -2 * h * I * Om, h * cz]];
    const K = [[Mx * wx * wx, 0, 0, 0], [0, My * wy * wy, 0, 0], [0, 0, h * I * (wz2 - Om * Om), h * cz * Om], [0, 0, -h * cz * Om, h * I * (wz2 - Om * Om)]];
    const Mi = N.inv(M), A = N.zeros(8), MK = N.matmul(Mi, K), MC = N.matmul(Mi, Cm);
    for (let p = 0; p < 4; p++) { A[p][p + 4] = 1; for (let q = 0; q < 4; q++) { A[p + 4][q] = -MK[p][q]; A[p + 4][q + 4] = -MC[p][q]; } }
    return ceig(A);
  };
  // Deutsch criterion for each support mode at its resonance centre Ω − ωζ = ω_support
  const deutsch = (w, c) => {
    const g = (Om) => Om - wz(Om) - w, top = 40 * w + 10; if (!(g(top) > 0)) return null;
    const Om = N.brent(g, 0, top, 1e-10), nu = wz(Om) / Om, req = (nb / 4) * S * S * w * w * ((1 - nu) / nu);
    return { Om, nu, req, czReq: req / c };
  };
  return { S, I, nb, wx, wy, cx, cy, wz, eig, deutsch, crit: 2 * I };
}
const groundRes = {
  id: 'groundres', title: 'Helicopter ground resonance (Coleman model)', fidelity: 'reduced-order',
  summary: 'Mechanical instability of a rotor with lag hinges on a flexible undercarriage: regressing lag mode against the airframe-on-gear modes, the damping needed to suppress it and the stability map against rotor speed.',
  equations: ['Euler–Lagrange equations', 'Damped free vibration equations', 'Eigenvalue equations', 'Linear multi-degree-of-freedom vibration equations'],
  applicable: (c) => (c.meta.type === 'helicopter' && c.rotor.R_m > 0 && c.rotor.n_blades >= 3 ? true : 'Ground resonance applies to helicopters with three or more lag-hinged (soft in-plane) blades; two-bladed teetering rotors and stiff in-plane propellers are outside this model.'),
  inputs: [
    { key: 'R', label: 'Rotor radius', unit: 'm', default: 8.18, min: 0.1, group: 'Rotor' },
    { key: 'n_blades', label: 'Number of blades', unit: '', default: 4, min: 3, max: 9, step: 1, discrete: true, group: 'Rotor' },
    { key: 'm_blade', label: 'Blade mass', unit: 'kg', default: 116, min: 0.01, group: 'Rotor' },
    { key: 'e_lag', label: 'Lag hinge offset / R', unit: '-', default: 0.047, min: 0.005, max: 0.3, group: 'Rotor' },
    { key: 'k_lag', label: 'Lag spring stiffness', unit: 'N·m/rad', default: 0, min: 0, group: 'Rotor', help: 'Elastomeric bearing or flexure stiffness about the lag hinge' },
    { key: 'zeta_lag', label: 'Lag damper ratio', unit: '-', default: 0.25, min: 0, max: 2, group: 'Damping', help: 'Fraction of critical lag damping at operating speed; the damper constant is reported' },
    { key: 'rpm', label: 'Operating rotor speed', unit: 'rpm', default: 258, min: 1, group: 'Rotor' },
    { key: 'M_x', label: 'Effective airframe mass at the hub, lateral/roll mode', unit: 'kg', default: 2800, min: 1, group: 'Airframe on gear' },
    { key: 'f_x', label: 'Lateral/roll mode frequency on the gear', unit: 'Hz', default: 1.6, min: 0.1, group: 'Airframe on gear' },
    { key: 'zeta_x', label: 'Lateral/roll mode damping ratio', unit: '-', default: 0.12, min: 0, max: 1, group: 'Damping', help: 'Oleo and tyre damping; 0.05–0.2' },
    { key: 'M_y', label: 'Effective airframe mass at the hub, longitudinal/pitch mode', unit: 'kg', default: 6000, min: 1, group: 'Airframe on gear' },
    { key: 'f_y', label: 'Longitudinal/pitch mode frequency on the gear', unit: 'Hz', default: 3.9, min: 0.1, group: 'Airframe on gear' },
    { key: 'zeta_y', label: 'Longitudinal/pitch mode damping ratio', unit: '-', default: 0.12, min: 0, max: 1, group: 'Damping' },
    { key: 'rpm_max_frac', label: 'Top speed / operating', unit: '-', default: 1.15, min: 1, max: 2, group: 'Rotor' },
    { key: 'nSpeeds', label: 'Rotor speed points', unit: '', default: 60, min: 20, max: 600, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c) => {
    const m = c.mass.mtow_kg, hHub = 0.2 * c.rotor.R_m + c.fuselage.dia_m, hCg = 0.6 * c.fuselage.dia_m, kt = c.gear.tyre_k_Npm * c.gear.n_main * c.gear.tyres_per_strut, wb = c.gear.wheelbase_m;
    const Ir = c.mass.Ixx + m * hCg * hCg, Ip = c.mass.Iyy + m * hCg * hCg, Kr = kt * (c.gear.track_m / 2) ** 2, Kp = kt * (0.2 * wb) ** 2 + 0.5 * kt * (0.8 * wb) ** 2;
    return { R: c.rotor.R_m, n_blades: c.rotor.n_blades, m_blade: c.rotor.blade_mass_kg || undefined, e_lag: c.rotor.hinge_offset || undefined, rpm: c.rotor.rpm, zeta_x: 0.2, zeta_y: 0.2, zeta_lag: 0.6, M_x: Ir / hHub ** 2, f_x: Math.sqrt(Kr / Ir) / TAU, M_y: Ip / hHub ** 2, f_y: Math.sqrt(Kp / Ip) / TAU };
  },
  run(i) {
    const g = grModel(i), Om0 = rads(i.rpm), cz = i.zeta_lag * g.crit * g.wz(Om0), ns = Math.round(i.nSpeeds), top = Om0 * i.rpm_max_frac, sp = N.linspace(top / ns, top, ns);
    const grow = sp.map((o) => N.amax(g.eig(o, cz).map((l) => l[0]))), k = N.argmax(grow), worst = grow[k], pct = sp.map((o) => (100 * o) / Om0);
    const modes = [{ name: 'Lateral/roll', d: g.deutsch(g.wx, g.cx) }, { name: 'Longitudinal/pitch', d: g.deutsch(g.wy, g.cy) }].filter((m) => m.d);
    const inRange = modes.filter((m) => m.d.Om <= top), ref = (inRange.length ? inRange : modes).slice().sort((p, q) => q.d.czReq - p.d.czReq)[0];
    // stability map over rotor speed and lag damping; it also brackets the smallest stabilising lag damper constant
    const cRef = g.crit * g.wz(Om0), zl = N.linspace(0, Math.max(2 * i.zeta_lag, 0.6), 12), spm = N.linspace(top / 18, top, 18), map = zl.map((z) => spm.map((o) => N.amax(g.eig(o, z * cRef).map((l) => l[0]))));
    const rowMax = map.map((r) => N.amax(r)), jU = rowMax.findLastIndex((v) => v > 1e-9);
    let czMin = 0;
    if (jU >= 0) { // refine on the speeds around the least stable column of the last unstable row
      const kc = N.argmax(map[jU]), loS = spm[Math.max(0, kc - 1)], hiS = spm[Math.min(spm.length - 1, kc + 1)], spc = N.linspace(loS, hiS, 7), gmax = (c) => N.amax(spc.map((o) => N.amax(g.eig(o, c).map((l) => l[0]))));
      let lo = zl[jU] * cRef, hi = jU + 1 < zl.length ? zl[jU + 1] * cRef : 2 * lo + cRef; for (let it = 0; it < 6 && gmax(hi) > 0; it++) { lo = hi; hi *= 2; }
      czMin = gmax(hi) > 0 ? Infinity : gmax(lo) <= 0 ? lo : N.brent(gmax, lo, hi, 2e-3 * hi);
    }
    const margin = czMin > 0 ? cz / czMin - 1 : Infinity, nu0 = g.wz(Om0) / Om0, warnings = [];
    if (!(czMin > 0)) warnings.push('The system is stable in the whole speed range even without lag dampers, so the damping margin is unbounded.');
    // unstable band
    let band = [NaN, NaN]; for (let j = 0; j < ns; j++) if (grow[j] > 1e-9) { if (!fin(band[0])) band[0] = pct[j]; band[1] = pct[j]; }
    if (worst > 1e-9) warnings.push(`The rotor–airframe system is unstable between ${band[0].toFixed(0)}% and ${band[1].toFixed(0)}% rotor speed (fastest growth doubles amplitude in ${(Math.LN2 / worst).toFixed(2)} s).`);
    if (nu0 >= 1) warnings.push('The lag frequency is above 1/rev (stiff in-plane): the regressing lag mode cannot coalesce with a support mode and ground resonance does not occur.');
    if (!inRange.length && modes.length) warnings.push('No resonance centre lies inside the rotor speed range; the margin is quoted for the nearest support mode.');
    return {
      kpis: [
        { key: 'gr_margin', label: 'Ground-resonance damping margin', value: margin, unit: '-', status: margin > 0.25 ? 'ok' : margin > 0 ? 'warn' : 'bad', note: 'Lag damper constant / smallest constant that keeps every eigenvalue stable over the speed range − 1' },
        { key: 'gr_growth_1s', label: 'Largest eigenvalue real part in the speed range', value: worst, unit: '1/s', status: worst <= 1e-9 ? 'ok' : 'bad', note: 'Positive means divergent oscillation' },
        { key: 'gr_worst_speed_pct', label: 'Rotor speed of least stability', value: pct[k], unit: '% NR' },
        { key: 'lag_freq_rev', label: 'Lag frequency at operating speed', value: nu0, unit: '/rev', status: nu0 < 1 ? 'warn' : 'ok', note: 'Below 1/rev = soft in-plane, susceptible' },
        { key: 'res_speed_x_pct', label: 'Resonance centre, lateral/roll mode', value: modes[0]?.name === 'Lateral/roll' ? (100 * modes[0].d.Om) / Om0 : NaN, unit: '% NR' },
        { key: 'res_speed_y_pct', label: 'Resonance centre, longitudinal/pitch mode', value: (100 * (modes.find((m) => m.name !== 'Lateral/roll')?.d.Om ?? NaN)) / Om0, unit: '% NR' },
        { key: 'c_lag_Nms', label: 'Lag damper constant (per blade)', value: cz, unit: 'N·m·s/rad' },
        { key: 'c_lag_req_Nms', label: 'Lag damper constant required (eigenvalue analysis)', value: czMin, unit: 'N·m·s/rad' },
        { key: 'c_lag_deutsch_Nms', label: 'Lag damper constant required (Deutsch criterion)', value: ref ? ref.d.czReq : 0, unit: 'N·m·s/rad', note: ref ? `${ref.name} mode at its resonance centre; a light-damping estimate that is optimistic for heavily damped lag motion` : '' },
        { key: 'damping_product', label: 'Damping product c_support · c_lag', value: ref ? (ref.name === 'Lateral/roll' ? g.cx : g.cy) * cz : 0, unit: 'N²·s²' },
        { key: 'damping_product_req', label: 'Deutsch required damping product', value: ref ? ref.d.req : 0, unit: 'N²·s²' },
      ],
      plots: [
        { type: 'line', title: 'Stability against rotor speed', xlabel: 'Rotor speed [% NR]', ylabel: 'Largest real part [1/s]', series: [{ name: 'With the lag damper', x: pct, y: grow }, { name: 'Lag damper removed', x: spm.map((o) => (100 * o) / Om0), y: map[0], style: 'dash' }], annotations: [{ y: 0, label: 'Stability boundary' }, { x: 100, label: 'NR' }] },
        { type: 'line', title: 'Coleman diagram (uncoupled frequencies)', xlabel: 'Rotor speed [% NR]', ylabel: 'Frequency in the fixed frame [Hz]', series: [{ name: 'Regressing lag |Ω − ωζ|', x: pct, y: sp.map((o) => Math.abs(o - g.wz(o)) / TAU) }, { name: 'Progressing lag Ω + ωζ', x: pct, y: sp.map((o) => (o + g.wz(o)) / TAU) }, { name: 'Lateral/roll mode', x: pct, y: sp.map(() => i.f_x), style: 'dash' }, { name: 'Longitudinal/pitch mode', x: pct, y: sp.map(() => i.f_y), style: 'dash' }] },
        { type: 'heat', title: 'Stability map: growth rate', xlabel: 'Rotor speed [% NR]', ylabel: 'Lag damper ratio [-]', zlabel: 'Largest real part [1/s]', x: spm.map((o) => (100 * o) / Om0), y: zl, z: map, contours: 10, diverging: true },
      ],
      warnings, models: ['Ground resonance model (Coleman–Feingold): hub translation with cyclic lag in multiblade coordinates', 'Deutsch damping criterion', 'Eigenvalue stability analysis'],
      assumptions: ['Three or more identical blades, isotropic rotor: constant coefficients in the fixed frame', 'Rigid blades with lag hinge and linear viscous lag damper; uniform blade mass', 'Two uncoupled airframe-on-gear modes represented by effective mass, frequency and damping at the hub', 'No aerodynamics (air resonance is not covered)', 'Default airframe modes are rigid-body-on-tyres estimates without oleo flexibility; default damping ratios are typical values'],
    };
  },
  convergence: { param: 'nSpeeds', label: 'Rotor speed points', levels: [15, 30, 60, 120], metric: 'gr_growth_1s' },
  calibration: { params: [{ key: 'zeta_x', min: 0, max: 0.6 }, { key: 'zeta_lag', min: 0, max: 1.5 }, { key: 'f_x', min: 0.2, max: 20 }], sweep: 'rpm', target: 'gr_growth_1s', note: 'Measured modal damping of the regressing lag mode from ground run-up shake tests' },
  verify() {
    const b = { R: 8, n_blades: 4, m_blade: 10, e_lag: 0.05, k_lag: 0, zeta_lag: 0, rpm: 250, M_x: 3000, f_x: 2, zeta_x: 0.03, M_y: 1e7, f_y: 300, zeta_y: 0.01, rpm_max_frac: 1.2, nSpeeds: 40 };
    const g = grModel(b), Om = 20, l = g.eig(Om, 0, 0).map((v) => Math.abs(v[1])).sort((p, q) => p - q), wz = g.wz(Om), d = g.deutsch(g.wx, g.cx);
    // numerical neutral-stability lag damping at the resonance centre against the Deutsch requirement
    const czN = N.brent((cz) => N.amax(N.linspace(0.9 * d.Om, 1.1 * d.Om, 41).map((o) => N.amax(g.eig(o, cz).map((v) => v[0])))), 0.05 * d.czReq, 5 * d.czReq, 1e-6 * d.czReq);
    return [
      N.check('Uncoupled regressing lag frequency Ω − ωζ', l[2], Om - wz, 1e-6, 'Multiblade coordinate transformation'),
      N.check('Uncoupled progressing lag frequency Ω + ωζ', N.amax(l.filter((v) => v < 100)), Om + wz, 1e-6, 'Multiblade coordinate transformation'),
      N.check('Neutral-stability lag damping agrees with the Deutsch criterion (light damping)', czN, d.czReq, 0.02, 'Deutsch (1946): c_x·c_ζ = (N/4)·S²·ω_x²·(1−ν)/ν, a two-mode light-damping approximation'),
    ];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (o.gr_growth_1s > 1e-9 || o.gr_margin < 0) out.push({ severity: 'critical', title: 'Ground resonance is predicted', detail: `Damping margin ${fin(o.gr_margin) ? (100 * o.gr_margin).toFixed(0) + '%' : 'n/a'}; least stable at ${o.gr_worst_speed_pct.toFixed(0)}% rotor speed.`, action: `Raise the lag damper constant above ${fin(o.c_lag_req_Nms) ? o.c_lag_req_Nms.toFixed(0) : 'the analysed range of'} N·m·s/rad per blade, increase oleo/tyre damping, or move the airframe-on-gear frequency away from the regressing lag frequency. Check with one damper failed and with low tyre and oleo pressure.`, basis: 'Deutsch criterion and eigenvalue analysis of the Coleman model' });
    else if (o.gr_margin < 0.25) out.push({ severity: 'warn', title: 'Thin ground-resonance damping margin', detail: `Margin is ${(100 * o.gr_margin).toFixed(0)}%.`, action: 'Damper wear, cold elastomer, a failed damper or a soft tyre can remove this margin; increase damping or demonstrate stability for those cases.', basis: 'Eigenvalue analysis of the Coleman model' });
    return out;
  },
};

export default {
  id: 'vibration', n: 10,
  tagline: 'Where the natural frequencies, critical speeds and resonances are, and how much margin and damping separates them from what excites them.',
  analyses: [blade, wingModal, forced, shaft, torsional, random, groundRes],
  consumes: [{ from: 'fea', keys: ['EI_root_Nm2', 'GJ_root_Nm2', 'wing_struct_mass_kg'], why: 'Wing stiffness and structural mass for the modal model' }],
  provides: [
    { key: 'f1_Hz', label: 'First bending frequency', unit: 'Hz' }, { key: 'f2_Hz', label: 'Second bending frequency', unit: 'Hz' }, { key: 'f3_Hz', label: 'Third bending frequency', unit: 'Hz' },
    { key: 'f_torsion_Hz', label: 'First torsion frequency', unit: 'Hz' }, { key: 'crit_speed_rpm', label: 'First critical speed', unit: 'rpm' }, { key: 'gr_margin', label: 'Ground-resonance damping margin', unit: '-' },
  ],
  handoff: [
    { model: 'Full-aircraft finite element modal model (free-free, fuselage and empennage modes, engine pylons)', why: 'Needs a detailed 3-D structural model and mass distribution; the native model is a cantilever beam', tool: 'General-purpose FE solver with a normal-modes solution, updated from a ground vibration test' },
    { model: 'Flexible multibody and rotor–gearbox coupled vibration', why: 'Gear mesh stiffness, housing flexibility and contact non-linearity are beyond a lumped model', tool: 'Multibody dynamics / dedicated gearbox dynamics code' },
    { model: 'Multi-disc rotor-bearing models with oil-film, seal and squeeze-film coefficients', why: 'Speed-dependent cross-coupled bearing coefficients and distributed shaft mass require a rotor FE code and bearing data', tool: 'Rotordynamics FE package' },
    { model: 'Air resonance and coupled rotor–fuselage aeromechanical stability', why: 'Requires blade flap–lag–torsion aerodynamics coupled with body modes', tool: 'Comprehensive rotorcraft analysis code' },
    { model: 'Non-linear resonance and contact vibration (freeplay, rub, friction dampers)', why: 'Needs non-linear time-domain or harmonic-balance solvers with test-derived joint models', tool: 'Non-linear dynamics FE / harmonic balance tools' },
    { model: 'Operational modal identification', why: 'Works on measured response data, not on a design model', tool: 'Operational modal analysis software with flight or ground test data' },
  ],
};
