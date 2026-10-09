// Suite 1 — Aerodynamics and Computational Fluid Dynamics.
// Native solvers: 2-D linear-vortex panel method with integral boundary layer, 3-D vortex-lattice and lifting-line
// wing, component drag build-up with Korn wave drag, 1-D finite-volume Euler (shock tube) with an exact Riemann
// solver, shock/expansion relations, a 2-D incompressible Navier–Stokes projection solver (lid-driven cavity) and a
// 1-D Spalart–Allmaras wall-turbulence solve with a first-cell-height calculator, a body-fitted 2-D RANS solver for the wing
// section (core/solvers/rans2d.js) coupled strip by strip to the vortex lattice for viscous wing forces, and — through the native 3-D
// Navier–Stokes kernel in core/solvers/cfd3d.js — DNS (Taylor–Green vortex), LES (turbulent channel), RANS with the
// Spalart–Allmaras model around an immersed wing, fuselage or imported surface, and unsteady bluff-body flow.
// What still needs an external solver (body-fitted 3-D, compressible, moving meshes) is in `handoff`.

import * as N from '../core/numerics.js';
import { isa, G0, GAMMA } from '../core/atmosphere.js';
import { createSolver, wallFriction } from '../core/solvers/cfd3d.js';
import { createRans2d, cGrid, boxGrid, plateGrid, resampleSection, wallSpacing } from '../core/solvers/rans2d.js';

const PI = Math.PI;
const fixedWing = (c) => (c.wing.S_m2 > 0 && c.wing.b_m > 0 ? true : 'This analysis needs a lifting wing; the current case is a pure rotorcraft. Use the airfoil analysis for the blade section and Suite 6 for the rotor.');
const even = (n, lo, hi) => 2 * Math.round(N.clamp(n, lo, hi) / 2);
const thin = (x, y, m = 300) => { const s = Math.max(1, Math.ceil(x.length / m)); return s === 1 ? [x, y] : [x.filter((_, k) => k % s === 0), y.filter((_, k) => k % s === 0)]; };
const ATM = [
  { key: 'alt_m', label: 'Altitude', unit: 'm', default: 0, min: -500, max: 25000, group: 'Flow', help: 'Pressure altitude of the analysis point (ISA)' },
  { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -60, max: 50, group: 'Flow' },
];

// ---- NACA 4- and 5-digit section geometry -------------------------------------------------------
/** Parse a NACA designation. Returns thickness, camber-line function x -> [yc, dyc/dx] and a validity flag. */
export function parseNaca(code, tcOverride = 0) {
  const s = String(code ?? '').replace(/[^0-9]/g, '');
  let m = 0, p = 0, t = 0.12, five = null, ok = true;
  if (s.length === 4) { m = +s[0] / 100; p = +s[1] / 10; t = +s.slice(2) / 100; if (!(m > 0 && p > 0)) { m = 0; p = 0; } }
  else if (s.length === 5 && s[2] === '0' && +s[1] >= 1 && +s[1] <= 5) { const k = +s[1] - 1; five = { r: [0.058, 0.126, 0.2025, 0.29, 0.391][k], k1: ([361.4, 51.64, 15.957, 6.643, 3.23][k] * (+s[0] * 0.15)) / 0.3 }; t = +s.slice(3) / 100; }
  else ok = false;
  if (tcOverride > 0) t = tcOverride;
  if (!(t > 0.01)) t = 0.12;
  t = Math.min(t, 0.4);
  const camber = (x) => {
    if (five) { const { r, k1 } = five; return x < r ? [(k1 / 6) * (x ** 3 - 3 * r * x * x + r * r * (3 - r) * x), (k1 / 6) * (3 * x * x - 6 * r * x + r * r * (3 - r))] : [((k1 * r ** 3) / 6) * (1 - x), (-k1 * r ** 3) / 6]; }
    if (!m) return [0, 0];
    return x < p ? [(m / (p * p)) * (2 * p * x - x * x), ((2 * m) / (p * p)) * (p - x)] : [(m / (1 - p) ** 2) * (1 - 2 * p + 2 * p * x - x * x), ((2 * m) / (1 - p) ** 2) * (p - x)];
  };
  const thick = (x) => 5 * t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1036 * x ** 4); // closed trailing edge
  return { t, ok, camber, thick, cambered: !!(m || five), name: ok ? `NACA ${s}` : 'NACA 0012 (fallback)', rLE: 1.1019 * t * t };
}
/** Surface nodes, clockwise from the trailing edge along the lower surface (cosine spacing), n panels. */
function nacaNodes(af, n) {
  const X = [], Y = [];
  for (let k = 0; k <= n; k++) {
    const x = 0.5 * (1 + Math.cos((2 * PI * k) / n)), [yc, dyc] = af.camber(x), yt = af.thick(x), th = Math.atan(dyc), up = k > n / 2 ? 1 : -1;
    X.push(x - up * yt * Math.sin(th)); Y.push(yc + up * yt * Math.cos(th));
  }
  X[0] = X[n] = 1; Y[0] = Y[n] = 0; X[n / 2] = 0; Y[n / 2] = 0;
  return { X, Y };
}
/** Thin-airfoil theory (Glauert): zero-lift angle [rad] and quarter-chord moment of a camber line. */
export function thinAirfoil(af) {
  const dz = (th) => af.camber(0.5 * (1 - Math.cos(th)))[1];
  const a0L = (-1 / PI) * N.simpson((th) => dz(th) * (Math.cos(th) - 1), 0, PI, 2000);
  const A1 = (2 / PI) * N.simpson((th) => dz(th) * Math.cos(th), 0, PI, 2000), A2 = (2 / PI) * N.simpson((th) => dz(th) * Math.cos(2 * th), 0, PI, 2000);
  return { alphaL0: a0L, cm_c4: (PI / 4) * (A2 - A1) };
}

// ---- linear-strength vortex panel method (Kuethe & Chow) ---------------------------------------
const _k = [0, 0, 0, 0];
/** Influence of panel j (linear vortex, end strengths γj, γj+1) at (x, y) resolved along a direction with sine si, cosine ci. */
function coef(P, x, y, si, ci, j) {
  const dx = x - P.X[j], dy = y - P.Y[j], sj = P.sn[j], cj = P.cs[j], Sj = P.S[j];
  const A = -dx * cj - dy * sj, B = dx * dx + dy * dy, C = si * cj - ci * sj, D = ci * cj + si * sj, E = dx * sj - dy * cj;
  const F = Math.log(1 + (Sj * (Sj + 2 * A)) / B), G = Math.atan2(E * Sj, B + A * Sj);
  const s2 = 2 * sj * cj, c2 = cj * cj - sj * sj, sI = si * c2 - ci * s2, cI = ci * c2 + si * s2, Pp = dx * sI + dy * cI, Q = dx * cI - dy * sI;
  _k[1] = D + (0.5 * Q * F) / Sj - ((A * C + D * E) * G) / Sj; _k[0] = 0.5 * D * F + C * G - _k[1];
  _k[3] = C + (0.5 * Pp * F) / Sj + ((A * D - C * E) * G) / Sj; _k[2] = 0.5 * C * F - D * G - _k[3];
}
/** Factorise the panel system once; the solution for any α is a linear combination of the α = 0° and 90° solutions. */
export function panelSolve(X, Y) {
  const M = X.length - 1, P = { M, X, Y, xc: [], yc: [], S: [], sn: [], cs: [] };
  for (let j = 0; j < M; j++) { const dx = X[j + 1] - X[j], dy = Y[j + 1] - Y[j], S = Math.hypot(dx, dy); P.S.push(S); P.sn.push(dy / S); P.cs.push(dx / S); P.xc.push(0.5 * (X[j] + X[j + 1])); P.yc.push(0.5 * (Y[j] + Y[j + 1])); }
  const AN = N.zeros(M + 1);
  for (let i = 0; i < M; i++) for (let j = 0; j < M; j++) {
    if (i === j) { _k[0] = -1; _k[1] = 1; } else coef(P, P.xc[i], P.yc[i], P.sn[i], P.cs[i], j);
    AN[i][j] += _k[0]; AN[i][j + 1] += _k[1];
  }
  AN[M][0] = 1; AN[M][M] = 1; // Kutta condition: zero net vortex strength at the trailing edge
  const f = N.lu(AN), gc = N.luSolve(f, [...P.sn, 0]), gs = N.luSolve(f, [...P.cs.map((v) => -v), 0]);
  P.gc = gc; P.gs = gs; P.sn0 = [0]; for (let j = 0; j < M; j++) P.sn0.push(P.sn0[j] + P.S[j]); // arc length at the nodes
  return P;
}
/**
 * Inviscid solution at angle of attack al [rad]. The vortex-sheet strength at a node equals the surface speed there
 * (zero velocity inside the body), which is second-order accurate on the true surface; Cp and loads use these nodal values.
 */
export function panelAt(P, al) {
  const ca = Math.cos(al), sa = Math.sin(al), g = P.gc.map((v, j) => v * ca + P.gs[j] * sa), vt = g.map((v) => 2 * PI * v);
  return { g, vt, cl: 4 * PI * circ(P, g), ...loads(P, vt, al) };
}
const circ = (P, g) => { let G = 0; for (let j = 0; j < P.M; j++) G += 0.5 * (g[j] + g[j + 1]) * P.S[j]; return G; };
/** Pressure loads from nodal surface speeds: Cp, pressure lift and drag, normal force and quarter-chord moment. */
function loads(P, vt, al) {
  const cp = vt.map((v) => 1 - v * v), ca = Math.cos(al), sa = Math.sin(al); let fx = 0, fy = 0, m4 = 0;
  for (let j = 0; j < P.M; j++) { const c = 0.5 * (cp[j] + cp[j + 1]) * P.S[j]; fx += c * P.sn[j]; fy -= c * P.cs[j]; m4 += c * ((P.xc[j] - 0.25) * P.cs[j] + P.yc[j] * P.sn[j]); }
  return { cp, clp: fy * ca - fx * sa, cdp: fx * ca + fy * sa, cn: fy, cm4: m4 };
}
/** Velocity (u, v)/V∞ at an off-body point. */
function panelField(P, g, ca, sa, x, y) {
  let u = ca, v = sa;
  for (let j = 0; j < P.M; j++) { coef(P, x, y, 0, 1, j); u += _k[2] * g[j] + _k[3] * g[j + 1]; v += _k[0] * g[j] + _k[1] * g[j + 1]; }
  return [u, v];
}
const inside = (X, Y, x, y) => { let c = false; for (let i = 0, j = X.length - 1; i < X.length; j = i++) if (Y[i] > y !== Y[j] > y && x < ((X[j] - X[i]) * (y - Y[i])) / (Y[j] - Y[i]) + X[i]) c = !c; return c; };

// ---- integral boundary layer: Thwaites (laminar), Michel (transition), Head (turbulent) --------
const headH1 = (H) => (H <= 1.6 ? 3.3 + 0.8234 * (H - 1.1) ** -1.287 : 3.3 + 1.5501 * (H - 0.6778) ** -3.064);
const headH = (H1) => (H1 <= 3.32 ? 4 : H1 >= 5.3 ? 1.1 + ((H1 - 3.3) / 0.8234) ** (-1 / 1.287) : 0.6778 + ((H1 - 3.3) / 1.5501) ** (-1 / 3.064));
/**
 * March one surface from the stagnation point. st = {s, ue, x}: arc length/c, edge speed/V∞ and chordwise station.
 * Returns momentum thickness θ/c, shape factor, cf, transition and separation stations and the Squire–Young drag.
 */
export function boundaryLayer(st, Re, xtr = 1, Hsep = 2.2, xEnd = 0.99) {
  const { s, ue, x } = st, n = s.length, th = new Array(n).fill(0), H = new Array(n).fill(2.59), cf = new Array(n).fill(0);
  const dU = (j) => (j === 0 ? (ue[1] - ue[0]) / (s[1] - s[0]) : j === n - 1 ? (ue[j] - ue[j - 1]) / (s[j] - s[j - 1]) : (ue[j + 1] - ue[j - 1]) / (s[j + 1] - s[j - 1]));
  const stag = ue[0] <= 1e-6;
  let I = 0, jT = -1, lamSep = false, jSep = -1;
  if (stag) { th[0] = Math.sqrt(0.075 / (Re * Math.max(dU(0), 1e-9))); H[0] = 2.24; }
  for (let j = 1; j < n; j++) {
    const ds = s[j] - s[j - 1], u = Math.max(ue[j], 1e-4);
    I += j === 1 && stag ? (u ** 5 * ds) / 6 : 0.5 * (u ** 5 + ue[j - 1] ** 5) * ds;
    const t = Math.sqrt((0.45 * I) / (Re * u ** 6)), lam = N.clamp(t * t * Re * dU(j), -0.1, 0.1), Ret = Re * u * t, Rex = Math.max(Re * u * s[j], 1);
    th[j] = t; H[j] = lam >= 0 ? 2.61 - 3.75 * lam + 5.24 * lam * lam : 2.088 + 0.0731 / (lam + 0.14);
    cf[j] = (2 * (lam >= 0 ? 0.22 + 1.57 * lam - 1.8 * lam * lam : 0.22 + 1.402 * lam + (0.018 * lam) / (lam + 0.107))) / Math.max(Ret, 1e-9);
    if (lam <= -0.09) { lamSep = true; jT = j; break; }
    if (Ret > 1.174 * (1 + 22400 / Rex) * Rex ** 0.46 || x[j] >= xtr) { jT = j; break; }
  }
  if (jT >= 0) {
    let t = th[jT], H1 = headH1(1.4), sep = false;
    for (let j = jT + 1; j < n; j++) {
      const ds = s[j] - s[j - 1], dud = (ue[j] - ue[j - 1]) / ds, ns = 4, h = ds / ns;
      const rhs = (tt, hh, u) => {
        const Hh = sep ? Hsep : headH(hh), c = sep ? 0 : 0.246 * 10 ** (-0.678 * Hh) * Math.max(Re * u * tt, 1) ** -0.268, dt = c / 2 - ((Hh + 2) * tt * dud) / u;
        return [dt, sep ? 0 : (0.0306 * (Math.max(hh, 3.05) - 3) ** -0.6169) / tt - hh * (dt / tt + dud / u), c];
      };
      let c = 0;
      for (let k = 0; k < ns; k++) {
        const u0 = Math.max(ue[j - 1] + dud * h * k, 0.02), um = Math.max(u0 + 0.5 * dud * h, 0.02), k1 = rhs(t, H1, u0), k2 = rhs(Math.max(t + 0.5 * h * k1[0], 1e-9), H1 + 0.5 * h * k1[1], um);
        t = Math.max(t + h * k2[0], 1e-9); H1 = Math.max(H1 + h * k2[1], 3.0); c = k2[2];
        if (!sep && headH(H1) >= Hsep) { sep = true; jSep = j; }
      }
      th[j] = t; H[j] = sep ? Hsep : headH(H1); cf[j] = sep ? 0 : c;
    }
  }
  let jE = n - 1; while (jE > 1 && x[jE] > xEnd) jE--;
  const xSep = jSep >= 0 && x[jSep] < 0.97 ? x[jSep] : 1;
  return { th, H, cf, dstar: th.map((v, j) => v * H[j]), jT, lamSep, xTr: jT >= 0 ? x[jT] : 1, xSep, dsE: th[jE] * H[jE], cd: 2 * th[jE] * Math.max(ue[jE], 0.05) ** ((H[jE] + 5) / 2) };
}
/** Split the panel solution into upper and lower surface stations starting at the stagnation point. */
function surfaces(P, vt) {
  let k = -1;
  for (let i = 1; i < P.M - 1; i++) if (vt[i] < 0 && vt[i + 1] >= 0 && (k < 0 || P.X[i] < P.X[k])) k = i;
  if (k < 0) k = P.M / 2 - 1;
  const f = N.clamp(vt[k + 1] !== vt[k] ? -vt[k] / (vt[k + 1] - vt[k]) : 0.5, 0.02, 0.98), s0 = P.sn0[k] + f * P.S[k], x0 = P.X[k] + f * (P.X[k + 1] - P.X[k]);
  const up = { s: [0], ue: [0], x: [x0] }, lo = { s: [0], ue: [0], x: [x0] };
  for (let i = k + 1; i <= P.M; i++) { up.s.push(P.sn0[i] - s0); up.ue.push(Math.max(vt[i], 1e-4)); up.x.push(P.X[i]); }
  for (let i = k; i >= 0; i--) { lo.s.push(s0 - P.sn0[i]); lo.ue.push(Math.max(-vt[i], 1e-4)); lo.x.push(P.X[i]); }
  return { up, lo, xStag: x0, k };
}
/**
 * One viscous section point: inviscid panel solution, both boundary layers (one-way coupled), Squire–Young drag.
 * Viscous lift = inviscid lift × Kirchhoff separated-flow factor − decambering by the trailing-edge displacement
 * thickness difference (the displacement surface tilts the effective camber line by (δ*u − δ*l)/c).
 */
function sectionPoint(P, al, o) {
  const inv = panelAt(P, al), sf = surfaces(P, inv.vt), bu = boundaryLayer(sf.up, o.Re, o.xtU, o.Hsep), bl = boundaryLayer(sf.lo, o.Re, o.xtL, o.Hsep);
  const f = inv.cl >= 0 ? bu.xSep : bl.xSep, kir = ((1 + Math.sqrt(f)) / 2) ** 2;
  return { inv, sf, bu, bl, f, cl: inv.cl * kir - o.kd * o.cla * (bu.dsE - bl.dsE), cd: bu.cd + bl.cd };
}
const cpKT = (cp, M) => { const b = Math.sqrt(1 - M * M), den = b + ((M * M) / (1 + b)) * (cp / 2); return den > 0.02 ? cp / den : -1e3; };
const cpCrit = (M, g = GAMMA) => (2 / (g * M * M)) * (((2 + (g - 1) * M * M) / (g + 1)) ** (g / (g - 1)) - 1);

const airfoil = {
  id: 'airfoil', title: '2-D airfoil: panel method with integral boundary layer', fidelity: 'numerical',
  summary: 'Pressure distribution, lift, moment, profile drag, transition, separation and stall of a NACA section from a linear-vortex panel solution coupled one-way to laminar and turbulent integral boundary layers.',
  equations: ['Laplace equation', 'Potential-flow equation', 'Kutta–Joukowski lift theorem', 'Circulation theorem', 'Bernoulli equation within its validity limits', 'Thin-airfoil theory', 'Prandtl–Glauert compressibility correction', 'Prandtl boundary-layer equations', 'Panel methods'],
  inputs: [
    { key: 'airfoil', label: 'NACA section', type: 'text', default: '2412', group: 'Geometry', help: '4-digit (2412, 0012, 4415) or 5-digit (23012) designation' },
    { key: 'tc', label: 'Thickness ratio override', unit: '-', default: 0, min: 0, max: 0.4, group: 'Geometry', help: '0 uses the thickness in the designation' },
    { key: 'chord', label: 'Chord', unit: 'm', default: 1.5, min: 0.005, group: 'Geometry' },
    { key: 'sweep_deg', label: 'Sweep (simple sweep theory)', unit: 'deg', default: 0, min: 0, max: 60, group: 'Geometry', help: 'The section sees V·cos Λ and M·cos Λ' },
    { key: 'role', label: 'Section role', type: 'select', options: ['Wing section', 'Rotor or propeller blade section'], default: 'Wing section', group: 'Geometry', help: 'Blade sections of wingless vehicles publish their section coefficients as the vehicle lift data' },
    { key: 'alpha_deg', label: 'Angle of attack', unit: 'deg', default: 3, min: -15, max: 25, group: 'Flow' },
    { key: 'V', label: 'Freestream speed', unit: 'm/s', default: 70, min: 1, group: 'Flow' },
    ...ATM,
    { key: 'xtr_upper', label: 'Forced transition, upper', unit: 'x/c', default: 1, min: 0.01, max: 1, group: 'Boundary layer', help: '1 = free transition (Michel). 0.05 represents a rough or contaminated leading edge' },
    { key: 'xtr_lower', label: 'Forced transition, lower', unit: 'x/c', default: 1, min: 0.01, max: 1, group: 'Boundary layer' },
    { key: 'H_sep', label: 'Turbulent separation shape factor', unit: '-', default: 2.2, min: 1.8, max: 3, group: 'Boundary layer', help: 'Head-method separation criterion; 1.8–2.4 in the literature' },
    { key: 'k_decamber', label: 'Decambering factor', unit: '-', default: 1, min: 0, max: 2, group: 'Boundary layer', help: 'Scales the lift lost to trailing-edge displacement thickness; 0 switches the correction off' },
    { key: 'nPanels', label: 'Panels', unit: '', default: 160, min: 20, max: 400, step: 2, discrete: true, group: 'Numerics' },
    { key: 'field', label: 'Compute the velocity field map', type: 'bool', default: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const base = { alt_m: c.atm.alt_m, dISA: c.atm.dISA_K };
    if (c.wing.S_m2 > 0) return { ...base, airfoil: c.wing.airfoil, tc: c.wing.tc, chord: d.mac || undefined, sweep_deg: Math.abs(c.wing.sweep_deg), alpha_deg: c.flight.alpha_deg, V: c.flight.V_ms, role: 'Wing section' };
    if (!(c.rotor.R_m > 0)) return base;
    const a = isa(c.atm.alt_m, c.atm.dISA_K), nR = c.meta.type === 'helicopter' ? 1 : Math.max(1, c.prop.n_eng), CT = d.W / nR / (a.rho * d.A_disk * d.v_tip ** 2), cl = d.solidity > 0 ? (6 * CT) / d.solidity : 0.5;
    return { ...base, airfoil: '0012', tc: 0.12, chord: c.rotor.chord_m, sweep_deg: 0, V: 0.75 * d.v_tip, alpha_deg: N.clamp(N.deg(cl / 6.2), 1, 10), role: 'Rotor or propeller blade section' };
  },
  run(i) {
    const a = isa(i.alt_m, i.dISA), cosL = Math.cos(N.rad(i.sweep_deg)), Vn = i.V * cosL, M = Vn / a.a, Mc = Math.min(M, 0.95), Re = (Vn * i.chord) / a.nu, beta = Math.sqrt(1 - Mc * Mc);
    const af = parseNaca(i.airfoil, i.tc), g = nacaNodes(af, even(i.nPanels, 20, 400)), P = panelSolve(g.X, g.Y), al = N.rad(i.alpha_deg);
    // lift-curve slope and zero-lift angle from the exactly linear (in sin α, cos α) inviscid solution
    const c0 = panelAt(P, 0), c1 = panelAt(P, N.rad(1)), cla = (c1.cl - c0.cl) / N.rad(1), aL0 = -c0.cl / cla, zl = panelAt(P, aL0), ta = thinAirfoil(af);
    const o = { Re, xtU: i.xtr_upper, xtL: i.xtr_lower, Hsep: i.H_sep, kd: i.k_decamber, cla }, pt = sectionPoint(P, al, o), inv = pt.inv;
    // low-speed viscous polar and stall (Kirchhoff lift on the Head-method separation point)
    const As = N.range(29, (k) => k - 5), pol = As.map((ad) => (ad === i.alpha_deg ? pt : sectionPoint(P, N.rad(ad), o)));
    let km = N.argmax(pol.map((p) => p.cl)), clmax = pol[km].cl, aStall = As[km];
    if (km > 0 && km < As.length - 1) { const y0 = pol[km - 1].cl, y1 = pol[km].cl, y2 = pol[km + 1].cl, dn = y0 - 2 * y1 + y2; if (dn < 0) { const dx = (0.5 * (y0 - y2)) / dn; aStall += dx; clmax = y1 - 0.25 * (y0 - y2) * dx; } }
    // compressibility: Kármán–Tsien on Cp, Prandtl–Glauert on integrated loads
    const cpC = inv.cp.map((v) => cpKT(v, Mc)), cpMinI = N.amin(inv.cp), cpMin = cpKT(cpMinI, Mc), kMin = N.argmin(inv.cp);
    const Mcr = N.findRoot((m) => cpKT(cpMinI, m) - cpCrit(m), 0.05, 0.99, 94, 1e-6);
    const cl = pt.cl / beta, clInv = inv.cl / beta, cm4 = inv.cm4 / beta, xcp = Math.abs(inv.cn) > 1e-4 ? 0.25 - inv.cm4 / inv.cn : NaN;
    const warnings = [];
    if (!af.ok) warnings.push(`"${i.airfoil}" is not a supported NACA 4-digit or non-reflex 5-digit designation; NACA 0012 thickness form was used.`);
    if (M > Mcr) warnings.push(`The section Mach number ${M.toFixed(3)} exceeds the critical Mach number ${Mcr.toFixed(3)}: there is supersonic flow and a shock on the section. Kármán–Tsien is not valid here and wave drag is not included; NACA 4/5-digit sections are also not representative of supercritical transonic sections.`);
    if (M > 0.95) warnings.push('Supersonic section Mach number: use the shock–expansion analysis instead.');
    if (Re < 1e5) warnings.push(`Chord Reynolds number ${Re.toExponential(2)} is below 1e5: laminar separation bubbles dominate and the Michel/Head closures are outside their calibrated range.`);
    if (pt.f < 1) warnings.push(`Turbulent separation is predicted at x/c = ${pt.f.toFixed(2)} on the suction side; lift uses the Kirchhoff separated-flow model and drag is under-predicted because the pressure drag of the separated wake is not resolved.`);
    if (Number.isNaN(xcp)) warnings.push('Centre of pressure is undefined at zero normal force.');
    if (km === As.length - 1) warnings.push('No lift maximum was found below 23°: the stall estimate is a lower bound.');
    const up = pt.bu, lo = pt.bl, half = P.M / 2, xu = g.X.slice(half), xl = g.X.slice(0, half + 1).reverse();
    const plots = [
      { type: 'line', title: `Pressure distribution, ${af.name} at α = ${i.alpha_deg}°`, xlabel: 'x/c [-]', ylabel: '−Cp [-]', series: [{ name: 'Upper surface', x: xu, y: cpC.slice(half).map((v) => -v) }, { name: 'Lower surface', x: xl, y: cpC.slice(0, half + 1).reverse().map((v) => -v) }], annotations: M > 0.3 ? [{ y: -cpCrit(Mc), label: 'Sonic (−Cp*)' }] : [] },
      { type: 'line', title: 'Skin-friction coefficient', xlabel: 'x/c [-]', ylabel: 'cf [-]', series: [{ name: 'Upper surface', x: pt.sf.up.x.slice(1), y: up.cf.slice(1) }, { name: 'Lower surface', x: pt.sf.lo.x.slice(1), y: lo.cf.slice(1) }], annotations: [{ x: up.xTr, label: 'Transition (upper)' }] },
      { type: 'line', title: 'Boundary-layer thicknesses, upper surface', xlabel: 'x/c [-]', ylabel: 'Thickness / chord ×10³ [-]', series: [{ name: 'Displacement δ*', x: pt.sf.up.x, y: up.dstar.map((v) => 1e3 * v) }, { name: 'Momentum θ', x: pt.sf.up.x, y: up.th.map((v) => 1e3 * v) }] },
      { type: 'line', title: 'Low-speed lift curve', xlabel: 'Angle of attack [deg]', ylabel: 'cl [-]', series: [{ name: 'Inviscid panel', x: As, y: pol.map((p) => p.inv.cl), style: 'dash' }, { name: 'With separation (Kirchhoff)', x: As, y: pol.map((p) => p.cl) }], annotations: [{ x: aStall, label: 'Stall' }] },
      { type: 'line', title: 'Low-speed profile-drag polar', xlabel: 'cd [-]', ylabel: 'cl [-]', series: [{ name: 'Squire–Young', x: pol.map((p) => p.cd), y: pol.map((p) => p.cl) }] },
    ];
    if (i.field) {
      const xs = N.linspace(-0.5, 1.5, 61), ys = N.linspace(-0.6, 0.6, 37), ca = Math.cos(al), sa = Math.sin(al);
      const z = ys.map((y) => xs.map((x) => { if (inside(g.X, g.Y, x, y)) return 0; const [u, v] = panelField(P, inv.g, ca, sa, x, y); return Math.min(Math.hypot(u, v), 3); }));
      plots.push({ type: 'heat', title: 'Velocity magnitude around the section (incompressible potential flow)', xlabel: 'x/c [-]', ylabel: 'y/c [-]', zlabel: 'V/V∞ [-]', x: xs, y: ys, z, contours: 14, equalAspect: true, overlay: [{ name: af.name, x: g.X, y: g.Y }] });
    } else plots.push({ type: 'line', title: af.name, xlabel: 'x/c [-]', ylabel: 'y/c [-]', equalAspect: true, series: [{ name: 'Section', x: g.X, y: g.Y }] });
    const rotor = i.role !== 'Wing section', cd0 = pol[N.argmin(pol.map((p) => Math.abs(p.inv.cl)))].cd, ldm = N.amax(pol.map((p) => p.cl / p.cd));
    return {
      kpis: [
        { key: 'cl', label: 'Section lift coefficient', value: cl, unit: '-' },
        { key: 'cd', label: 'Section profile drag coefficient', value: pt.cd, unit: '-', note: 'Squire–Young; friction plus form drag of attached flow' },
        { key: 'cm_c4', label: 'Moment about quarter chord', value: cm4, unit: '-' },
        { key: 'cl_inviscid', label: 'Inviscid lift coefficient', value: clInv, unit: '-' },
        { key: 'cla_2d_per_rad', label: 'Section lift-curve slope', value: cla / beta, unit: '1/rad' },
        { key: 'alpha_L0_deg', label: 'Zero-lift angle', value: N.deg(aL0), unit: 'deg' },
        { key: 'cp_min', label: 'Minimum pressure coefficient', value: cpMin, unit: '-', note: `at x/c = ${g.X[kMin].toFixed(3)}` },
        { key: 'x_cp_frac', label: 'Centre of pressure', value: xcp, unit: 'x/c' },
        { key: 'M_crit', label: 'Critical Mach number', value: Mcr, unit: '-', status: M < Mcr ? 'ok' : 'bad', note: 'Section Mach must stay below this for shock-free flow' },
        { key: 'x_tr_upper', label: 'Transition, upper surface', value: up.xTr, unit: 'x/c', note: up.lamSep ? 'laminar separation bubble' : up.jT >= 0 ? 'natural or forced' : 'laminar to the trailing edge' },
        { key: 'x_tr_lower', label: 'Transition, lower surface', value: lo.xTr, unit: 'x/c' },
        { key: 'x_sep', label: 'Turbulent separation, suction side', value: pt.f, unit: 'x/c', status: pt.f >= 1 ? 'ok' : pt.f > 0.8 ? 'warn' : 'bad', note: '1 = attached to the trailing edge' },
        { key: 'clmax_section', label: 'Section maximum lift (low speed)', value: clmax, unit: '-', note: 'Kirchhoff model on the Head separation point: an estimate, not a measurement' },
        { key: 'alpha_stall_deg', label: 'Section stall angle', value: aStall, unit: 'deg' },
        { key: 'ld_section_max', label: 'Best section lift-to-drag', value: ldm, unit: '-' },
        { key: 'Re_chord', label: 'Chord Reynolds number', value: Re, unit: '-' },
        { key: 'mach_section', label: 'Section Mach number', value: M, unit: '-' },
      ],
      plots,
      tables: [{ title: 'Low-speed section polar', columns: ['α [deg]', 'cl inviscid', 'cl', 'cd', 'x_tr upper', 'x_sep'], rows: pol.filter((_, k) => k % 2 === 1).map((p, k) => [As[2 * k + 1], +p.inv.cl.toFixed(4), +p.cl.toFixed(4), +p.cd.toFixed(5), +p.bu.xTr.toFixed(3), +p.f.toFixed(3)]) }],
      outputs: { Cm_ac: zl.cm4 / beta, cd0_section: cd0, alpha_L0_thin_deg: N.deg(ta.alphaL0), cm_c4_thin: ta.cm_c4, cl_pressure: inv.clp / beta, cd_pressure_inviscid: inv.cdp, ...(rotor ? { CL: cl, CD: pt.cd, CLa_per_rad: cla / beta, CLmax: clmax, LD_max: ldm } : {}) },
      warnings,
      models: ['Linear-strength vortex panel method with Kutta condition', 'Thwaites laminar integral method', 'Michel transition criterion', 'Head entrainment method with Ludwieg–Tillmann skin friction', 'Squire–Young profile drag', 'Kirchhoff trailing-edge separation lift model', 'Kármán–Tsien and Prandtl–Glauert compressibility corrections'],
      assumptions: ['Steady, 2-D, irrotational outer flow; boundary layer coupled one way (no displacement-thickness feedback)', 'A laminar separation bubble is treated as immediate transition', 'Boundary layer marched with the incompressible edge velocity', 'Stall, post-separation drag and maximum lift are engineering estimates', rotor ? 'Blade section evaluated at 75% radius in hover; its coefficients are published as the vehicle lift data because there is no wing' : 'Simple sweep theory for the swept-wing section'],
    };
  },
  convergence: { param: 'nPanels', label: 'Panels on the section', levels: [40, 80, 160, 320], metric: 'cl_inviscid' },
  calibration: { params: [{ key: 'H_sep', min: 1.8, max: 2.8 }, { key: 'k_decamber', min: 0, max: 2 }, { key: 'xtr_upper', min: 0.02, max: 1 }], sweep: 'alpha_deg', target: 'cl', note: 'Supply a measured lift curve (α, cl) at the same Reynolds number; the separation criterion and transition location are the uncertain closures.' },
  verify() {
    // circular cylinder: exact Cp = 1 - 4 sin²θ, zero lift and drag
    const n = 80, X = N.range(n + 1, (k) => 0.5 + 0.5 * Math.cos((2 * PI * k) / n)), Y = N.range(n + 1, (k) => -0.5 * Math.sin((2 * PI * k) / n)), cyl = panelAt(panelSolve(X, Y), 0);
    const P = panelSolve(...Object.values(nacaNodes(parseNaca('0012'), 200))), p4 = panelAt(P, N.rad(4)), s = (panelAt(P, N.rad(1)).cl - panelAt(P, 0).cl) / N.rad(1);
    const P24 = panelSolve(...Object.values(nacaNodes(parseNaca('2412'), 200))), z0 = panelAt(P24, 0), z1 = panelAt(P24, N.rad(1)), ta = thinAirfoil(parseNaca('2412'));
    const Re = 1e6, xs = N.linspace(0, 1, 201), lam = boundaryLayer({ s: xs, ue: xs.map(() => 1), x: xs }, Re, 2, 2.2, 1), tur = boundaryLayer({ s: xs, ue: xs.map(() => 1), x: xs.map((v) => v + 1) }, 1e7, 1, 2.2, 3);
    return [
      N.check('Cylinder minimum Cp = −3', N.amin(cyl.cp), -3, 1e-4, 'Exact potential flow about a circle, Cp = 1 − 4 sin²θ'),
      N.check('Cylinder lift is zero', cyl.cl, 0, 1e-8, 'Symmetry'),
      N.check('Inviscid drag is zero (d’Alembert)', p4.cdp, 0, 1.5e-3, 'd’Alembert paradox; pressure-integration error of the panel discretisation'),
      N.check('Pressure-integrated lift equals Kutta–Joukowski lift', p4.clp, p4.cl, 5e-3, 'Kutta–Joukowski theorem'),
      N.check('NACA 0012 lift slope ≈ 2π(1 + 0.77 t/c)', s, 2 * PI * (1 + 0.77 * 0.12), 0.015, 'Thickness-corrected thin-airfoil slope (approximate closed form)'),
      N.check('Thin-airfoil zero-lift angle of NACA 2412', N.deg(ta.alphaL0), -2.0775, 2e-3, 'Glauert integral, closed form for the 4-digit camber line'),
      N.check('Panel zero-lift angle approaches thin-airfoil theory', N.deg(-z0.cl / ((z1.cl - z0.cl) / N.rad(1))), N.deg(ta.alphaL0), 0.08, 'Thin-airfoil theory; thickness effect is second order'),
      N.check('Laminar flat-plate drag (Thwaites vs Blasius)', lam.cd, 1.328 / Math.sqrt(Re), 0.015, 'Blasius; Thwaites’ method is accurate to about 1%'),
      N.check('Turbulent flat-plate drag (Head vs Prandtl–Schlichting)', tur.cd, 0.455 / Math.log10(1e7) ** 2.58, 0.08, 'Prandtl–Schlichting correlation; Head with Ludwieg–Tillmann friction agrees within 5–8%'),
    ];
  },
  validation: [{ name: 'Thin-airfoil zero-lift angle, NACA 2412', source: 'Glauert thin-airfoil theory (closed form for the NACA 4-digit mean line)', inputs: { airfoil: '2412', tc: 0, V: 50, alt_m: 0, sweep_deg: 0, field: false }, sweep: { key: 'alpha_deg', values: [0] }, target: 'alpha_L0_thin_deg', observed: [-2.0775], tol_pct: 1 }],
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.mach_section > o.M_crit) out.push({ severity: 'warn', title: 'Section is operating above its critical Mach number', detail: `M = ${o.mach_section.toFixed(3)} against Mcrit = ${o.M_crit.toFixed(3)} at cl = ${o.cl.toFixed(2)}. Wave drag rises steeply beyond this point and every drag count is fuel and CO₂.`, action: 'Reduce thickness or section lift, add sweep, or move to a supercritical section and confirm with a transonic RANS solution.', basis: 'Critical pressure coefficient with Kármán–Tsien correction' });
    if (o.x_sep < 0.9) out.push({ severity: o.x_sep < 0.6 ? 'critical' : 'warn', title: 'Trailing-edge separation at the analysis point', detail: `Separation at x/c = ${o.x_sep.toFixed(2)}; stall is predicted near α = ${o.alpha_stall_deg.toFixed(1)}°.`, action: 'Lower the angle of attack, or use the stall margin in Suite 5 to set approach speeds; do not rely on attached-flow loads here.', basis: 'Head shape-factor separation criterion' });
    if (o.x_tr_upper > 0.3 && i.xtr_upper >= 1) out.push({ severity: 'info', title: 'Natural laminar flow is available on this section', detail: `Free transition at ${(100 * o.x_tr_upper).toFixed(0)}% chord (upper) and ${(100 * o.x_tr_lower).toFixed(0)}% (lower) gives cd = ${o.cd.toFixed(4)}.`, action: 'Protect it: keep the leading edge clean and smooth (insects, rain, de-icing fluid, paint steps). Re-run with forced transition at 5% chord to see the fuel penalty of a contaminated surface.', basis: 'Michel transition criterion' });
    if (o.Re_chord < 2e5) out.push({ severity: 'advise', title: 'Low Reynolds number section', detail: `Re = ${o.Re_chord.toExponential(2)}. Laminar separation bubbles control performance and noise.`, action: 'Prefer sections designed for low Reynolds number and consider a boundary-layer trip; validate with wind-tunnel data.', basis: 'Validity range of the Michel and Head correlations' });
    return out;
  },
};

// ---- 3-D wing: vortex-lattice method and Prandtl lifting line ----------------------------------
/** Upwash at (xm, ym) from a planar unit-strength horseshoe with bound segment 1→2 and legs trailing to x = +∞ (Biot–Savart). */
function horseshoe(xm, ym, x1, y1, x2, y2) {
  const a = xm - x1, b = ym - y1, c = xm - x2, d = ym - y2, r1 = Math.hypot(a, b), r2 = Math.hypot(c, d), den = a * d - c * b;
  const bound = Math.abs(den) < 1e-12 ? 0 : (((x2 - x1) * a + (y2 - y1) * b) / r1 - ((x2 - x1) * c + (y2 - y1) * d) / r2) / den;
  return (bound + (1 + a / r1) / (y1 - ym) - (1 + c / r2) / (y2 - ym)) / (4 * PI);
}
/**
 * Planar vortex-lattice solution of a trapezoidal wing. g = {b, cr, taper, tanL (quarter-chord), twist (tip, rad)}.
 * Compressibility by the Prandtl–Glauert (Göthert) stretching x → x/β. Returns unit-α and zero-α loadings.
 */
export function vlm(g, nS, nC, beta = 1, slope = () => 0, cosine = false) {
  const half = g.b / 2, n2 = 2 * nS, nP = n2 * nC, yE = cosine ? N.range(n2 + 1, (k) => -half * Math.cos((PI * k) / n2)) : N.linspace(-half, half, n2 + 1);
  const ym = N.range(n2, (j) => (cosine ? -half * Math.cos((PI * (j + 0.5)) / n2) : 0.5 * (yE[j] + yE[j + 1]))), dy = N.range(n2, (j) => yE[j + 1] - yE[j]);
  const chord = (y) => g.cr * (1 - ((1 - g.taper) * Math.abs(y)) / half), xle = (y) => g.cr / 4 + Math.abs(y) * g.tanL - chord(y) / 4;
  const pn = [];
  for (let js = 0; js < n2; js++) for (let ic = 0; ic < nC; ic++) {
    const ya = yE[js], yb = yE[js + 1], cm = chord(ym[js]);
    pn.push({ js, ic, ya, yb, ym: ym[js], cm, xa: xle(ya) + (chord(ya) * (ic + 0.25)) / nC, xb: xle(yb) + (chord(yb) * (ic + 0.25)) / nC, xc: xle(ym[js]) + (cm * (ic + 0.75)) / nC, xf: (ic + 0.75) / nC });
  }
  const A = N.zeros(nP);
  for (let i = 0; i < nP; i++) for (let j = 0; j < nP; j++) A[i][j] = horseshoe(pn[i].xc / beta, pn[i].ym, pn[j].xa / beta, pn[j].ya, pn[j].xb / beta, pn[j].yb);
  const f = N.lu(A), Ga = N.luSolve(f, pn.map(() => -1)), G0 = N.luSolve(f, pn.map((p) => slope(p.xf) - (g.twist * Math.abs(p.ym)) / half));
  const S = 0.5 * g.cr * (1 + g.taper) * g.b, AR = (g.b * g.b) / S;
  const post = (G) => {
    const strip = new Array(n2).fill(0); let mom = 0; // moment about x = 0 (root leading edge), nose-up positive, per unit q
    pn.forEach((p, k) => { strip[p.js] += G[k]; mom -= 2 * G[k] * dy[p.js] * 0.5 * (p.xa + p.xb); });
    return { strip, CL: (2 * N.dot(strip, dy)) / S, mom };
  };
  const trefftz = (strip) => { let D = 0; for (let i = 0; i < n2; i++) { let w = 0; for (let j = 0; j < n2; j++) w += (strip[j] / (2 * PI)) * (1 / (yE[j] - ym[i]) - 1 / (yE[j + 1] - ym[i])); D -= strip[i] * w * dy[i]; } return D / S; };
  return { pn, yE, ym, dy, S, AR, chord, xle, Ga, G0, a: post(Ga), z: post(G0), trefftz, n2, nC, solve: (rhs) => N.luSolve(f, rhs) };
}
/** Prandtl lifting line by Glauert's Fourier series (symmetric loading). cOverB(θ), alpha(θ) [rad], a0 section slope. */
export function liftingLine(AR, cOverB, alpha, a0 = 2 * PI, nT = 40) {
  const th = N.range(nT, (k) => ((k + 1) * PI) / (2 * nT)), A = th.map((t) => N.range(nT, (j) => { const n = 2 * j + 1, mu = (a0 * cOverB(t)) / 4; return Math.sin(n * t) * (Math.sin(t) + n * mu); }));
  const An = N.solve(A, th.map((t) => ((a0 * cOverB(t)) / 4) * alpha(t) * Math.sin(t))), sum = An.reduce((s, v, j) => s + (2 * j + 1) * v * v, 0);
  return { CL: PI * AR * An[0], CDi: PI * AR * sum, e: sum > 0 ? (An[0] * An[0]) / sum : 1, An };
}
const wingGeom = (i) => { const cr = (2 * i.S) / (i.b * (1 + i.taper)), lam = i.taper; return { b: i.b, cr, taper: lam, tanL: Math.tan(N.rad(i.sweep_deg)), twist: N.rad(i.twist_deg), mac: ((2 / 3) * cr * (1 + lam + lam * lam)) / (1 + lam), ymac: (i.b / 6) * ((1 + 2 * lam) / (1 + lam)) }; };

const wing = {
  id: 'wing', title: '3-D wing: vortex lattice and lifting line', fidelity: 'numerical',
  summary: 'Lift-curve slope, span loading, induced drag, span efficiency, aerodynamic centre and stall onset of a swept, tapered, twisted wing from a vortex-lattice solution, cross-checked with Prandtl lifting-line theory.',
  equations: ['Prandtl lifting-line theory', 'Lifting-surface integral equations', 'Kutta–Joukowski lift theorem', 'Circulation theorem', 'Prandtl–Glauert compressibility correction', 'Vortex lattice models', 'Lifting-line models'],
  applicable: fixedWing,
  inputs: [
    { key: 'S', label: 'Wing area', unit: 'm²', default: 16.2, min: 0.01, group: 'Geometry' },
    { key: 'b', label: 'Span', unit: 'm', default: 11, min: 0.05, group: 'Geometry' },
    { key: 'taper', label: 'Taper ratio', unit: '-', default: 0.6, min: 0.05, max: 1, group: 'Geometry', help: 'Tip chord / root chord' },
    { key: 'sweep_deg', label: 'Quarter-chord sweep', unit: 'deg', default: 0, min: -30, max: 60, group: 'Geometry' },
    { key: 'twist_deg', label: 'Tip twist (washout negative)', unit: 'deg', default: -2, min: -10, max: 5, group: 'Geometry' },
    { key: 'airfoil', label: 'NACA section', type: 'text', default: '2412', group: 'Geometry', help: 'Supplies the camber line of the lattice' },
    { key: 'clmax_section', label: 'Section maximum lift', unit: '-', default: 1.5, min: 0.3, max: 3, group: 'Geometry', help: 'From the airfoil analysis when it has been run' },
    { key: 'alpha_deg', label: 'Wing angle of attack', unit: 'deg', default: 4, min: -10, max: 25, group: 'Flow', help: 'Root-chord incidence to the freestream' },
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 60, min: 1, group: 'Flow' },
    ...ATM,
    { key: 'nSpan', label: 'Spanwise panels per half wing', unit: '', default: 16, min: 2, max: 60, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nChord', label: 'Chordwise panels', unit: '', default: 4, min: 1, max: 10, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up) => ({ S: c.wing.S_m2, b: c.wing.b_m, taper: c.wing.taper, sweep_deg: c.wing.sweep_deg, twist_deg: c.wing.twist_deg, airfoil: c.wing.airfoil, clmax_section: up.cfd?.clmax_section, alpha_deg: c.flight.alpha_deg + c.wing.incidence_deg, V: c.flight.V_ms, alt_m: c.atm.alt_m, dISA: c.atm.dISA_K }),
  run(i) {
    const a = isa(i.alt_m, i.dISA), M = i.V / a.a, beta = Math.sqrt(1 - Math.min(M, 0.9) ** 2), g = wingGeom(i), af = parseNaca(i.airfoil), al = N.rad(i.alpha_deg);
    const nS = Math.round(N.clamp(i.nSpan, 2, 60)), nC = Math.round(N.clamp(i.nChord, 1, 10)), L = vlm(g, nS, nC, beta, (xf) => af.camber(xf)[1], true);
    const CLa = L.a.CL, CL0 = L.z.CL, CL = CLa * al + CL0, strip = L.a.strip.map((v, j) => v * al + L.z.strip[j]), CDi = L.trefftz(strip), e = CDi > 1e-12 ? (CL * CL) / (PI * L.AR * CDi) : 1;
    // aerodynamic centre and moment about it (moments taken about the root leading edge, then transferred)
    const xac = -L.a.mom / (CLa * L.S), Cmac = (L.z.mom + CL0 * L.S * xac) / (L.S * g.mac), xLEmac = g.cr / 4 + g.ymac * g.tanL - g.mac / 4;
    // section lift and first stall: local cl reaches the section limit (reduced by cos Λ for sweep)
    const cl = strip.map((v, j) => (2 * v) / L.chord(L.ym[j])), clA = L.a.strip.map((v, j) => (2 * v) / L.chord(L.ym[j])), clZ = L.z.strip.map((v, j) => (2 * v) / L.chord(L.ym[j]));
    const lim = i.clmax_section * Math.cos(N.rad(i.sweep_deg)), aSt = clA.map((v, j) => (lim - clZ[j]) / v), js = N.argmin(aSt), CLmax = CLa * aSt[js] + CL0, eta = L.ym.map((y) => (2 * y) / i.b);
    const etaCp = N.sum(strip.map((v, j) => v * L.dy[j] * Math.abs(eta[j]))) / (N.dot(strip, L.dy) || 1);
    // lifting line (unswept) for comparison
    const cb = (t) => (g.cr / g.b) * (1 - (1 - g.taper) * Math.abs(Math.cos(t))), a0 = (2 * PI) / beta, aL0 = thinAirfoil(af).alphaL0;
    const ll1 = liftingLine(L.AR, cb, () => 1, a0), llp = liftingLine(L.AR, cb, (t) => al - aL0 + g.twist * Math.abs(Math.cos(t)), a0);
    const ell = strip.map((_, j) => ((4 * CL * L.S) / (PI * i.b)) * Math.sqrt(Math.max(0, 1 - eta[j] ** 2)) / 2);
    // lattice ΔCp map
    const nodes = [], tris = [], cnt = [], val = [], id = (js2, ic) => js2 * (nC + 1) + ic;
    for (let j = 0; j <= L.n2; j++) for (let k = 0; k <= nC; k++) { nodes.push([L.yE[j], -(L.xle(L.yE[j]) + (L.chord(L.yE[j]) * k) / nC)]); cnt.push(0); val.push(0); }
    L.pn.forEach((p, k) => { const dcp = (2 * (L.Ga[k] * al + L.G0[k]) * nC) / p.cm, c4 = [id(p.js, p.ic), id(p.js + 1, p.ic), id(p.js + 1, p.ic + 1), id(p.js, p.ic + 1)]; tris.push([c4[0], c4[1], c4[2]], [c4[0], c4[2], c4[3]]); for (const q of c4) { val[q] += dcp; cnt[q]++; } });
    const warnings = [];
    if (M > 0.7) warnings.push(`Mach ${M.toFixed(2)}: the Prandtl–Glauert rule is a linear subsonic correction and does not capture shocks or transonic loading changes.`);
    if (CL > CLmax) warnings.push(`The wing is beyond the predicted stall onset (CL ${CL.toFixed(2)} > ${CLmax.toFixed(2)}); the linear lattice loads are not valid here.`);
    if (Math.abs(i.sweep_deg) > 10) warnings.push('Lifting-line results are shown for reference only: classical lifting-line theory does not model sweep.');
    return {
      kpis: [
        { key: 'CL_wing', label: 'Wing lift coefficient', value: CL, unit: '-', status: CL <= CLmax ? 'ok' : 'bad' },
        { key: 'CDi', label: 'Induced drag coefficient', value: CDi, unit: '-', note: 'Trefftz-plane integration' },
        { key: 'CLa_per_rad', label: 'Wing lift-curve slope', value: CLa, unit: '1/rad' },
        { key: 'e_span', label: 'Span efficiency (inviscid)', value: e, unit: '-', status: e > 0.9 ? 'ok' : 'warn', note: 'Oswald efficiency including viscous and fuselage effects comes from the drag build-up' },
        { key: 'alpha_L0_wing_deg', label: 'Wing zero-lift angle', value: N.deg(-CL0 / CLa), unit: 'deg' },
        { key: 'CLmax_wing', label: 'Clean wing maximum lift (stall onset)', value: CLmax, unit: '-', note: 'Critical-section method: first section to reach clmax·cos Λ' },
        { key: 'alpha_stall_deg', label: 'Wing angle at stall onset', value: N.deg(aSt[js]), unit: 'deg' },
        { key: 'eta_stall', label: 'Spanwise station of first stall', value: Math.abs(eta[js]), unit: '2y/b', status: Math.abs(eta[js]) < 0.7 ? 'ok' : 'warn', note: 'Stall should start inboard of the ailerons (< 0.7)' },
        { key: 'x_ac_frac_mac', label: 'Aerodynamic centre', value: (xac - xLEmac) / g.mac, unit: 'x/MAC' },
        { key: 'Cm_ac', label: 'Moment about the aerodynamic centre', value: Cmac, unit: '-' },
        { key: 'eta_cp', label: 'Spanwise centre of lift', value: etaCp, unit: '2y/b', note: 'Elliptic loading gives 0.424; sets root bending moment' },
        { key: 'CLa_llt_per_rad', label: 'Lifting-line lift-curve slope', value: ll1.CL, unit: '1/rad' },
        { key: 'e_llt', label: 'Lifting-line span efficiency', value: llp.e, unit: '-' },
        { key: 'aspect_ratio', label: 'Aspect ratio', value: L.AR, unit: '-' },
      ],
      plots: [
        { type: 'line', title: 'Span loading', xlabel: 'Span station 2y/b [-]', ylabel: 'cl·c [m]', series: [{ name: 'Vortex lattice', x: eta, y: strip.map((v) => 2 * v) }, { name: 'Elliptic, same lift', x: eta, y: ell.map((v) => 2 * v), style: 'dash' }] },
        { type: 'line', title: 'Section lift coefficient and stall limit', xlabel: 'Span station 2y/b [-]', ylabel: 'cl [-]', series: [{ name: 'At this α', x: eta, y: cl }, { name: 'At stall onset', x: eta, y: clA.map((v, j) => v * aSt[js] + clZ[j]) }, { name: 'Section limit', x: eta, y: eta.map(() => lim), style: 'dash' }] },
        { type: 'tri', title: 'Lifting-pressure distribution ΔCp on the planform', xlabel: 'y [m]', ylabel: '−x [m]', zlabel: 'ΔCp [-]', nodes, tris, values: val.map((v, k) => v / (cnt[k] || 1)), equalAspect: true, edges: true },
        { type: 'line', title: 'Wing lift curve and induced-drag polar', xlabel: 'Angle of attack [deg]', ylabel: 'CL [-]', series: [{ name: 'Vortex lattice', x: [N.deg(-CL0 / CLa), N.deg(aSt[js])], y: [0, CLmax] }, { name: 'Lifting line', x: [N.deg(-CL0 / CLa), N.deg(aSt[js])], y: [0, ll1.CL * (aSt[js] + CL0 / CLa)], style: 'dash' }], annotations: [{ x: N.deg(aSt[js]), label: 'Stall onset' }] },
      ],
      outputs: { k_inviscid: 1 / (PI * L.AR * e) },
      warnings,
      models: ['Vortex-lattice method (horseshoe vortices, quarter-chord/three-quarter-chord rule, cosine spanwise spacing)', 'Trefftz-plane induced drag', 'Prandtl–Glauert (Göthert) compressibility transformation', 'Prandtl lifting-line theory (Glauert series)', 'Critical-section stall onset'],
      assumptions: ['Thin, planar lifting surface: dihedral, thickness, fuselage and nacelles are not represented', 'Linear attached flow; rigid wake trailing along the body axis', 'Camber line of the NACA section applied at every span station', 'Section maximum lift reduced by cos Λ for sweep (empirical)'],
    };
  },
  convergence: { param: 'nSpan', label: 'Spanwise panels per half wing', levels: [6, 12, 24, 48], metric: 'CLa_per_rad' },
  verify() {
    const L = vlm({ b: 5, cr: 1, taper: 1, tanL: 1, twist: 0 }, 4, 1), cb = (t) => (4 / (PI * 8)) * Math.sin(t), ll = liftingLine(8, cb, () => 1, 2 * PI, 30);
    const hi = vlm({ b: 20, cr: 1, taper: 1, tanL: 0, twist: 0 }, 40, 1, 1, () => 0, true), l20 = liftingLine(20, () => 1 / 20, () => 1, 2 * PI, 60);
    return [
      N.check('VLM lift slope, 45° swept AR 5 wing, 4 panels', L.a.CL, 3.443, 3e-3, 'Bertin & Cummings, Aerodynamics for Engineers, vortex-lattice worked example (CL = 1.096π·α)'),
      N.check('Lifting line, elliptic wing: CLα = a0/(1 + a0/πAR)', ll.CL, (2 * PI) / (1 + 2 / 8), 1e-6, 'Prandtl closed form'),
      N.check('Lifting line, elliptic wing: span efficiency = 1', ll.e, 1, 1e-6, 'Minimum induced drag theorem'),
      N.check('VLM approaches lifting line at high aspect ratio', hi.a.CL, l20.CL, 0.03, 'Consistency of two independent methods at AR 20 (lifting-surface slope is 2–3% below lifting line)'),
      N.check('Trefftz-plane induced drag ≥ elliptic minimum', hi.trefftz(hi.a.strip) / (hi.a.CL ** 2 / (PI * 20)) >= 1 ? 1 : 0, 1, 1e-12, 'Munk minimum induced drag'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.eta_stall > 0.7) out.push({ severity: 'warn', title: 'Stall starts outboard', detail: `First stall at ${(100 * o.eta_stall).toFixed(0)}% semispan, in the aileron region: roll control degrades and a wing drop is likely.`, action: 'Add washout, reduce taper, or use a higher-clmax tip section; re-run until stall starts inboard of about 60% semispan.', basis: 'Critical-section stall analysis' });
    if (o.e_span < 0.93) out.push({ severity: 'advise', title: 'Span loading departs from elliptic', detail: `Span efficiency ${o.e_span.toFixed(3)}; each 1% of efficiency is roughly 1% of induced drag, about 0.4% of cruise fuel and CO₂.`, action: 'Tune taper and twist (taper near 0.4–0.5 unswept) or study a span extension or winglet in Suite 23, trading against root bending moment.', basis: 'Munk minimum induced drag' });
    if (i.sweep_deg > 20 && o.eta_cp > 0.44) out.push({ severity: 'info', title: 'Swept wing carries load outboard', detail: `Centre of lift at ${(100 * o.eta_cp).toFixed(1)}% semispan.`, action: 'Check root bending and aeroelastic wash-out in Suites 2 and 3; flexible twist will unload the tip in flight.', basis: 'Span load centroid' });
    return out;
  },
};

// ---- component drag build-up --------------------------------------------------------------------
const cfLam = (Re) => 1.328 / Math.sqrt(Re);
const cfTurb = (Re, M = 0) => 0.455 / (Math.log10(Re) ** 2.58 * (1 + 0.144 * M * M) ** 0.65);
/** Mean flat-plate skin friction with a laminar run xl (fraction of length) and a roughness cut-off Reynolds number. */
export function cfMixed(Re, M, xl, lOverK) {
  Re = Math.max(Re, 1e3);
  if (xl >= 1) return cfLam(Re);
  const ct = cfTurb(Math.min(Re, lOverK > 0 ? 38.21 * lOverK ** 1.053 : Infinity), M);
  if (!(xl > 0)) return ct;
  const Rx = Math.max(Re * xl, 1e3);
  return Math.max(ct - xl * (cfTurb(Rx, M) - cfLam(Rx)), cfLam(Re));
}
/** Korn equation with Lock's fourth-power wave-drag rise. */
export function korn(M, CL, tc, cosL, kappa) {
  const Mdd = kappa / cosL - tc / cosL ** 2 - CL / (10 * cosL ** 3), Mcr = Mdd - (0.1 / 80) ** (1 / 3);
  return { Mdd, Mcr, CDw: M > Mcr ? 20 * (M - Mcr) ** 4 : 0 };
}
function buildup(i) {
  const a = isa(i.alt_m, i.dISA), M = i.V / a.a, q = 0.5 * a.rho * i.V ** 2, AR = (i.b * i.b) / i.S, lam = i.taper, cr = (2 * i.S) / (i.b * (1 + lam)), mac = ((2 / 3) * cr * (1 + lam + lam * lam)) / (1 + lam), cosL = Math.cos(N.rad(i.sweep_deg));
  const comp = [], Mf = Math.max(M, 0.2);
  const add = (name, Swet, L, FF, Q, xl) => { if (!(Swet > 0 && L > 0)) return; const Re = (i.V * L) / a.nu, Cf = cfMixed(Re, M, xl, L / i.ks_m); comp.push({ name, Swet, Re, Cf, FF, Q, CD: (Cf * FF * Q * Swet) / i.S }); };
  const ffS = (tc, cl) => (1 + (0.6 / i.xc_max) * tc + 100 * tc ** 4) * 1.34 * Mf ** 0.18 * cl ** 0.28;
  const Sexp = Math.max(0.5 * i.S, i.S - i.fus_D * cr * (1 - ((1 - lam) * i.fus_D) / (2 * i.b)));
  add('Wing', Sexp * (1.977 + 0.52 * i.tc), mac, ffS(i.tc, cosL), i.Q_wing, i.lam_wing);
  add('Horizontal tail', i.S_h * (1.977 + 0.52 * i.tc_tail), i.S_h / Math.max(i.b_h, 1e-6), ffS(i.tc_tail, 1), i.Q_tail, i.lam_wing);
  add('Vertical tail', i.S_v * (1.977 + 0.52 * i.tc_tail), i.S_v / Math.max(i.b_v, 1e-6), ffS(i.tc_tail, 1), i.Q_tail, i.lam_wing);
  const lf = Math.max(i.fus_L / Math.max(i.fus_D, 1e-6), 2.5);
  add('Fuselage', PI * i.fus_D * i.fus_L * (1 - 2 / lf) ** (2 / 3) * (1 + 1 / (lf * lf)), i.fus_L, 1 + 60 / lf ** 3 + lf / 400, 1, i.lam_fus);
  if (i.n_nac > 0) add('Nacelles / pods', i.n_nac * PI * i.nac_D * i.nac_L, i.nac_L, 1 + 0.35 / Math.max(i.nac_L / Math.max(i.nac_D, 1e-6), 1), i.Q_nac, 0);
  const clean = N.sum(comp.map((c) => c.CD)), CDexc = (clean * i.exc_pct) / 100, CDextra = i.dq_extra / i.S, CD0 = clean + CDexc + CDextra;
  const s = 1 - 2 * (i.fus_D / i.b) ** 2, e = 1 / (1 / (i.e_span * Math.max(s, 0.5)) + i.K_visc * PI * AR * CD0), k = 1 / (PI * AR * e);
  const CL = (i.n_load * i.mass_kg * G0) / (q * i.S), w = korn(M, CL, i.tc, cosL, i.kappa), CD = CD0 + k * CL * CL + w.CDw;
  const Swet = N.sum(comp.map((c) => c.Swet));
  return { a, M, q, AR, mac, cosL, comp, clean, CDexc, CDextra, CD0, e, k, CL, w, CD, Swet };
}
const dragBuildup = {
  id: 'buildup', title: 'Aircraft drag build-up and polar', fidelity: 'reduced-order',
  summary: 'Zero-lift drag from the wetted area, skin friction, form and interference factors of each component, plus induced and transonic wave drag, giving the aircraft polar, maximum lift-to-drag ratio and maximum lift.',
  equations: ['Prandtl boundary-layer equations', 'Prandtl–Glauert compressibility correction', 'Wing–fuselage interference models', 'Bernoulli equation within its validity limits'],
  applicable: fixedWing,
  inputs: [
    { key: 'mass_kg', label: 'Aircraft mass', unit: 'kg', default: 1111, min: 0.1, group: 'Flight condition' },
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 62, min: 1, group: 'Flight condition' },
    { key: 'n_load', label: 'Load factor', unit: 'g', default: 1, min: 0, max: 9, group: 'Flight condition' },
    ...ATM.map((f) => ({ ...f, group: 'Flight condition' })),
    { key: 'S', label: 'Wing reference area', unit: 'm²', default: 16.2, min: 0.01, group: 'Wing' },
    { key: 'b', label: 'Span', unit: 'm', default: 11, min: 0.05, group: 'Wing' },
    { key: 'taper', label: 'Taper ratio', unit: '-', default: 0.69, min: 0.05, max: 1, group: 'Wing' },
    { key: 'sweep_deg', label: 'Quarter-chord sweep', unit: 'deg', default: 0, min: 0, max: 60, group: 'Wing' },
    { key: 'tc', label: 'Wing thickness ratio', unit: '-', default: 0.12, min: 0.03, max: 0.25, group: 'Wing' },
    { key: 'xc_max', label: 'Chordwise position of maximum thickness', unit: 'x/c', default: 0.3, min: 0.15, max: 0.6, group: 'Wing', help: '0.3 for NACA 4-digit, 0.4–0.5 for laminar and supercritical sections' },
    { key: 'kappa', label: 'Korn technology factor', unit: '-', default: 0.87, min: 0.8, max: 0.97, group: 'Wing', help: '0.87 conventional sections, 0.95 supercritical' },
    { key: 'e_span', label: 'Inviscid span efficiency', unit: '-', default: 0.97, min: 0.5, max: 1.2, group: 'Wing', help: 'From the wing analysis when it has been run' },
    { key: 'clmax_section', label: 'Section maximum lift', unit: '-', default: 1.5, min: 0.3, max: 3, group: 'Wing' },
    { key: 'CLmax_wing', label: 'Wing maximum lift (0 = estimate)', unit: '-', default: 0, min: 0, max: 4, group: 'Wing', help: '0 uses 0.9·clmax·cos Λ; the wing analysis supplies a computed value' },
    { key: 'S_h', label: 'Horizontal tail area', unit: 'm²', default: 3.35, min: 0, group: 'Tail' },
    { key: 'b_h', label: 'Horizontal tail span', unit: 'm', default: 3.45, min: 0, group: 'Tail' },
    { key: 'S_v', label: 'Vertical tail area', unit: 'm²', default: 1.74, min: 0, group: 'Tail' },
    { key: 'b_v', label: 'Vertical tail height', unit: 'm', default: 1.5, min: 0, group: 'Tail' },
    { key: 'tc_tail', label: 'Tail thickness ratio', unit: '-', default: 0.1, min: 0.03, max: 0.2, group: 'Tail' },
    { key: 'fus_L', label: 'Fuselage length', unit: 'm', default: 8.3, min: 0, group: 'Bodies' },
    { key: 'fus_D', label: 'Fuselage diameter', unit: 'm', default: 1.2, min: 0, group: 'Bodies' },
    { key: 'n_nac', label: 'Nacelles or pods', unit: '', default: 0, min: 0, max: 16, step: 1, discrete: true, group: 'Bodies' },
    { key: 'nac_D', label: 'Nacelle diameter', unit: 'm', default: 1, min: 0.01, group: 'Bodies', help: 'Default is a class estimate from engine thrust or power; replace with the real figure' },
    { key: 'nac_L', label: 'Nacelle length', unit: 'm', default: 2, min: 0.01, group: 'Bodies' },
    { key: 'dq_extra', label: 'Extra drag area D/q', unit: 'm²', default: 0, min: 0, group: 'Bodies', help: 'Fixed landing gear, struts, booms, stopped rotors, antennas' },
    { key: 'lam_wing', label: 'Laminar fraction, wing and tails', unit: '-', default: 0.05, min: 0, max: 0.7, group: 'Surface', help: '0–0.1 metal transports, 0.3–0.5 smooth composite natural-laminar-flow wings' },
    { key: 'lam_fus', label: 'Laminar fraction, fuselage', unit: '-', default: 0, min: 0, max: 0.5, group: 'Surface' },
    { key: 'ks_m', label: 'Equivalent sand roughness', unit: 'm', default: 6.3e-6, min: 1e-7, max: 1e-3, group: 'Surface', help: '6e-6 smooth paint, 1e-5 camouflage paint, 4e-6 production sheet metal, 0.5e-6 polished composite' },
    { key: 'exc_pct', label: 'Excrescence and leakage allowance', unit: '%', default: 6, min: 0, max: 30, group: 'Surface', help: 'Gaps, steps, rivets, antennas, leaks: 3–5% modern transports, 5–10% light aircraft' },
    { key: 'Q_wing', label: 'Wing interference factor', unit: '-', default: 1.0, min: 0.9, max: 1.5, group: 'Interference' },
    { key: 'Q_tail', label: 'Tail interference factor', unit: '-', default: 1.04, min: 0.9, max: 1.5, group: 'Interference' },
    { key: 'Q_nac', label: 'Nacelle interference factor', unit: '-', default: 1.3, min: 0.9, max: 2, group: 'Interference', help: '1.5 mounted directly on the wing, 1.3 less than one diameter away, 1.0 beyond that' },
    { key: 'K_visc', label: 'Lift-dependent viscous drag factor', unit: '-', default: 0.38, min: 0, max: 1, group: 'Interference', help: 'Share of CD0 that grows with CL² in the Oswald efficiency' },
  ],
  defaults: (c, up, d) => {
    const p = c.prop, jet = p.type === 'turbofan' || p.type === 'turbojet', mdot = jet ? p.T0_N / (1000 / (1 + (p.type === 'turbojet' ? 0 : p.bpr)) ** 0.6) : 0;
    const nacD = jet ? 1.25 * Math.sqrt((4 * mdot) / (PI * 195)) : 0.15 * Math.max(p.P0_W / 1e3, 1) ** 0.25, nNac = jet || p.type === 'turboprop' || p.n_eng > 1 ? p.n_eng : 0;
    const fixedGear = c.aero.Vmo_ms > 0 && c.aero.Vmo_ms < 75 && c.mass.mtow_kg > 200 && c.gear.type !== 'skid';
    return {
      mass_kg: c.mass.mtow_kg, V: c.flight.V_ms, n_load: c.flight.n_load, alt_m: c.atm.alt_m, dISA: c.atm.dISA_K, S: c.wing.S_m2, b: c.wing.b_m, taper: c.wing.taper, sweep_deg: Math.abs(c.wing.sweep_deg), tc: c.wing.tc,
      kappa: c.aero.Mmo > 0.7 ? 0.95 : 0.87, xc_max: c.aero.Mmo > 0.7 ? 0.4 : 0.3, e_span: up.cfd?.e_span, clmax_section: up.cfd?.clmax_section, CLmax_wing: up.cfd?.CLmax_wing,
      S_h: c.htail.S_m2, b_h: c.htail.b_m, S_v: c.vtail.S_m2, b_v: c.vtail.b_m, tc_tail: c.htail.tc, fus_L: c.fuselage.len_m, fus_D: c.fuselage.dia_m, n_nac: nNac, nac_D: nNac ? nacD : undefined, nac_L: nNac ? (jet ? 1.7 : 3.5) * nacD : undefined,
      dq_extra: fixedGear ? 0.5 * (c.gear.n_main * c.gear.tyres_per_strut + 1) * 0.35 * c.gear.tyre_dia_m ** 2 : undefined, lam_wing: c.mass.mtow_kg < 5700 && c.aero.Vmo_ms < 60 ? 0.15 : undefined,
    };
  },
  run(i) {
    const r = buildup(i), { M, CD0, k, e, CL, CD, w } = r, LDmax = 1 / (2 * Math.sqrt(CD0 * k)), CLmd = Math.sqrt(CD0 / k);
    const CLmax = i.CLmax_wing > 0 ? i.CLmax_wing : 0.9 * i.clmax_section * r.cosL, cnt = (v) => +(1e4 * v).toFixed(1);
    const CLs = N.linspace(0, CLmax, 60), Ms = N.linspace(Math.min(0.3, M), Math.min(0.95, Math.max(0.9, M + 0.1)), 60);
    const warnings = [];
    if (M > w.Mdd) warnings.push(`Mach ${M.toFixed(3)} is beyond the drag-divergence Mach number ${w.Mdd.toFixed(3)} at CL = ${CL.toFixed(2)}: wave drag is ${cnt(w.CDw)} counts and rising steeply; the Korn/Lock model is only a trend there.`);
    if (M > 0.95) warnings.push('Supersonic flight: this subsonic build-up does not apply.');
    if (CL > CLmax) warnings.push(`Required CL ${CL.toFixed(2)} exceeds the clean maximum ${CLmax.toFixed(2)}: this speed, mass and altitude cannot be flown clean.`);
    if (r.Swet / i.S < 2.5 || r.Swet / i.S > 9) warnings.push(`Wetted-area ratio ${(r.Swet / i.S).toFixed(1)} is outside the usual 3–7 range for aeroplanes; check the component dimensions.`);
    return {
      kpis: [
        { key: 'CD0', label: 'Zero-lift drag coefficient', value: CD0, unit: '-', note: `${cnt(CD0)} drag counts` },
        { key: 'k_induced', label: 'Induced drag factor k', unit: '-', value: k },
        { key: 'e_oswald', label: 'Oswald efficiency', value: e, unit: '-', status: e > 0.7 ? 'ok' : 'warn' },
        { key: 'CL', label: 'Lift coefficient at this point', value: CL, unit: '-', status: CL <= CLmax ? 'ok' : 'bad' },
        { key: 'CD', label: 'Drag coefficient at this point', value: CD, unit: '-' },
        { key: 'LD', label: 'Lift-to-drag ratio at this point', value: CL / CD, unit: '-', status: CL / CD > 0.85 * LDmax ? 'ok' : 'warn', note: 'Within 15% of the maximum is good cruise matching' },
        { key: 'LD_max', label: 'Maximum lift-to-drag ratio', value: LDmax, unit: '-' },
        { key: 'CL_md', label: 'Lift coefficient for best L/D', value: CLmd, unit: '-' },
        { key: 'CD_wave', label: 'Wave drag coefficient', value: w.CDw, unit: '-', status: w.CDw < 0.002 ? 'ok' : 'warn', note: 'Above about 20 counts the aircraft is past drag divergence' },
        { key: 'M_dd', label: 'Drag-divergence Mach number', value: w.Mdd, unit: '-' },
        { key: 'CLmax', label: 'Clean maximum lift coefficient', value: CLmax, unit: '-' },
        { key: 'drag_N', label: 'Drag at this point', value: CD * r.q * i.S, unit: 'N' },
        { key: 'f_flat_plate_m2', label: 'Equivalent flat-plate area', value: CD0 * i.S, unit: 'm²' },
        { key: 'Swet_over_S', label: 'Wetted-area ratio', value: r.Swet / i.S, unit: '-' },
        { key: 'Cfe', label: 'Equivalent skin-friction coefficient', value: (CD0 * i.S) / r.Swet, unit: '-', note: 'Typically 0.0026–0.0035 jets, 0.0045–0.0055 light singles' },
        { key: 'mach', label: 'Flight Mach number', value: M, unit: '-' },
      ],
      plots: [
        { type: 'bar', title: 'Zero-lift drag build-up', ylabel: 'Drag counts (CD × 10⁴) [-]', categories: [...r.comp.map((c) => c.name), 'Excrescence', 'Extra D/q'], series: [{ name: 'CD0 contribution', y: [...r.comp.map((c) => 1e4 * c.CD), 1e4 * r.CDexc, 1e4 * r.CDextra] }] },
        { type: 'line', title: 'Drag polar', xlabel: 'CD [-]', ylabel: 'CL [-]', series: [{ name: 'Subcritical parabolic', x: CLs.map((c) => CD0 + k * c * c), y: CLs, style: 'dash' }, { name: `At Mach ${M.toFixed(2)} with wave drag`, x: CLs.map((c) => CD0 + k * c * c + korn(M, c, i.tc, r.cosL, i.kappa).CDw), y: CLs }], annotations: [{ y: CL, label: 'This point' }] },
        { type: 'line', title: 'Lift-to-drag ratio', xlabel: 'CL [-]', ylabel: 'L/D [-]', series: [{ name: `Mach ${M.toFixed(2)}`, x: CLs.slice(1), y: CLs.slice(1).map((c) => c / (CD0 + k * c * c + korn(M, c, i.tc, r.cosL, i.kappa).CDw)) }], annotations: [{ x: CL, label: 'This point' }, { x: CLmd, label: 'Best L/D' }] },
        { type: 'line', title: 'Drag rise with Mach number at constant lift', xlabel: 'Mach number [-]', ylabel: 'CD [-]', series: [{ name: `CL = ${CL.toFixed(2)}`, x: Ms, y: Ms.map((m) => CD0 + k * CL * CL + korn(m, CL, i.tc, r.cosL, i.kappa).CDw) }], annotations: [{ x: w.Mdd, label: 'Drag divergence' }] },
      ],
      tables: [{ title: 'Component build-up', columns: ['Component', 'Wetted area [m²]', 'Re [million]', 'Cf ×10³', 'Form factor', 'Interference', 'CD0 [counts]'], rows: [...r.comp.map((c) => [c.name, +c.Swet.toFixed(2), +(c.Re / 1e6).toFixed(2), +(1e3 * c.Cf).toFixed(3), +c.FF.toFixed(3), c.Q, cnt(c.CD)]), ['Excrescence and leakage', '', '', '', '', '', cnt(r.CDexc)], ['Extra drag area', '', '', '', '', '', cnt(r.CDextra)], ['Total', +r.Swet.toFixed(2), '', '', '', '', cnt(CD0)]] }],
      warnings,
      models: ['Component build-up: flat-plate skin friction × form factor × interference (empirical)', 'Blasius laminar and Prandtl–Schlichting turbulent skin friction with compressibility and roughness cut-off', 'Oswald efficiency from inviscid span efficiency, fuselage correction and lift-dependent viscous drag (empirical)', 'Korn equation with Lock fourth-power wave-drag rise (empirical)'],
      assumptions: ['Clean configuration: no flap, slat, spoiler or landing-gear drag unless entered as extra drag area', 'No trim drag, no propeller slipstream or jet interference', 'CD0 evaluated at the flight Reynolds and Mach number and held constant along the polar', 'Maximum lift is the clean wing value without tail download'],
    };
  },
  calibration: { params: [{ key: 'exc_pct', min: 0, max: 30 }, { key: 'lam_wing', min: 0, max: 0.6 }, { key: 'K_visc', min: 0, max: 0.8 }, { key: 'kappa', min: 0.84, max: 0.97 }], sweep: 'mass_kg', target: 'CD', note: 'Supply flight-test or wind-tunnel drag (CD against mass or CL at fixed speed); excrescence, laminar extent, viscous lift-dependent drag and the Korn factor are the uncertain constants.' },
  verify() {
    const base = { mass_kg: 1000, V: 50, n_load: 1, alt_m: 0, dISA: 0, S: 10, b: 10, taper: 1, sweep_deg: 0, tc: 0.12, xc_max: 0.3, kappa: 0.87, e_span: 1, clmax_section: 1.5, CLmax_wing: 0, S_h: 0, b_h: 0, S_v: 0, b_v: 0, tc_tail: 0.1, fus_L: 0, fus_D: 0, n_nac: 0, nac_D: 1, nac_L: 2, dq_extra: 0, lam_wing: 0, lam_fus: 0, ks_m: 1e-9, exc_pct: 0, Q_wing: 1, Q_tail: 1, Q_nac: 1, K_visc: 0 };
    const o = N.kv(dragBuildup.run(base)), r = buildup(base), kw = korn(0.75, 0, 0.12, 1, 0.87), h = 1e-5;
    return [
      N.check('Laminar plate: Cf = 1.328/√Re', cfMixed(4e5, 0, 1, 0), 1.328 / Math.sqrt(4e5), 1e-12, 'Blasius'),
      N.check('Wing-only CD0 = Cf·FF·Swet/S', o.CD0, (r.comp[0].Cf * r.comp[0].FF * r.comp[0].Swet) / 10, 1e-12, 'Definition of the build-up'),
      N.check('Ideal polar: k = 1/(π·AR)', o.k_induced, 1 / (PI * 10), 1e-12, 'Elliptic loading'),
      N.check('L/D max = 1/(2√(CD0·k))', o.LD_max, 0.5 / Math.sqrt(o.CD0 * o.k_induced), 1e-12, 'Parabolic polar'),
      N.check('Korn drag-divergence Mach (κ 0.87, t/c 0.12, CL 0)', kw.Mdd, 0.75, 1e-12, 'Korn equation'),
      N.check('dCD/dM = 0.1 at drag divergence', (korn(0.75 + h, 0, 0.12, 1, 0.87).CDw - korn(0.75 - h, 0, 0.12, 1, 0.87).CDw) / (2 * h), 0.1, 1e-6, 'Definition of drag divergence in Lock’s rise'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], big = res.tables[0].rows.slice(0, -3).sort((a, b) => b[6] - a[6])[0];
    if (big) out.push({ severity: 'info', title: `${big[0]} is the largest zero-lift drag item`, detail: `${big[6]} of ${(1e4 * o.CD0).toFixed(1)} counts. At cruise one drag count is about ${(100 / (1e4 * o.CD)).toFixed(2)}% of fuel burn and CO₂.`, action: 'Target wetted area and surface quality of this component first; then excrescences (gaps, steps, antennas).', basis: 'Component drag build-up' });
    if (o.LD < 0.85 * o.LD_max) out.push({ severity: 'advise', title: 'Cruise point is away from best lift-to-drag', detail: `L/D = ${o.LD.toFixed(1)} against a maximum of ${o.LD_max.toFixed(1)} (CL ${o.CL.toFixed(2)} versus ${o.CL_md.toFixed(2)} for best L/D).`, action: o.CL < o.CL_md ? 'Fly higher or slower, or reduce wing area in Suite 23; jets cruise best slightly below CL for maximum L/D.' : 'Fly faster or lower, or increase wing area.', basis: 'Breguet range: fuel burn is inversely proportional to L/D' });
    if (o.CD_wave > 0.002) out.push({ severity: 'warn', title: 'Operating beyond drag divergence', detail: `Wave drag is ${(1e4 * o.CD_wave).toFixed(0)} counts at Mach ${o.mach.toFixed(2)} (Mdd = ${o.M_dd.toFixed(3)}).`, action: 'Reduce cruise Mach, add sweep, thin the wing or use supercritical sections; confirm with transonic CFD and tunnel data.', basis: 'Korn equation' });
    if (i.lam_wing < 0.2 && o.mach < 0.7) out.push({ severity: 'info', title: 'Laminar-flow potential', detail: 'The wing is assumed almost fully turbulent.', action: 'Re-run with a laminar fraction of 0.3–0.5 to size the benefit of natural laminar flow against its surface-quality and contamination demands.', basis: 'Skin-friction build-up' });
    return out;
  },
};

// ---- body-fitted 2-D RANS of the section and the viscous-coupled wing ---------------------------
const WALLS = ['Wall functions (first cell at y⁺ ≈ 50)', 'Wall-resolved (first cell at y⁺ ≈ 1)'];
const RANS_GRID = [
  { key: 'wall', label: 'Wall treatment', type: 'select', options: WALLS, default: WALLS[0], group: 'Numerics', help: 'Wall functions need far fewer cells and are the fast default; wall-resolved grids need about 48 or more wall-normal cells and are more reliable near separation' },
  { key: 'nSurf', label: 'Cells around the section', unit: '', default: 56, min: 32, max: 320, step: 2, discrete: true, group: 'Numerics', help: '48–64 fast look, 96–128 for drag within a few percent, 192 or more for reference work' },
  { key: 'nNormal', label: 'Wall-normal cells', unit: '', default: 22, min: 16, max: 128, step: 1, discrete: true, group: 'Numerics', help: '20–24 with wall functions; 48–64 wall-resolved' },
  { key: 'nWake', label: 'Cells along the wake', unit: '', default: 12, min: 6, max: 96, step: 1, discrete: true, group: 'Numerics' },
  { key: 'radius', label: 'Far-field distance', unit: 'chords', default: 30, min: 10, max: 200, group: 'Numerics', help: 'A point-vortex correction is applied at the far field, so 20–30 chords is enough' },
];
const gridOf = (i) => { const wfn = i.wall !== WALLS[1]; return { wfn, yplus: wfn ? 50 : 1, nSurf: even(i.nSurf, 32, 320), nNormal: Math.round(N.clamp(i.nNormal, wfn ? 16 : 32, 128)), nWake: Math.round(N.clamp(i.nWake, 6, 96)), radius: N.clamp(i.radius, 10, 200) }; };
/** C-grid and solver for a NACA section at chord Reynolds number Re and angle of attack al [rad]. */
function sectionSolver(af, Re, al, G, clInit = 0, extra = {}) {
  const nd = nacaNodes(af, 400), sec = resampleSection(nd.X, nd.Y, G.nSurf), ReS = N.clamp(Re, 2e4, 2e9);
  const grid = cGrid({ X: sec.X, Y: sec.Y, nWake: G.nWake, nj: G.nNormal, h1: wallSpacing(ReS, G.yplus), radius: G.radius });
  return { grid, sec, s: createRans2d(grid, { Re: ReS, alpha: al, wallFunction: G.wfn, clInit, ...extra }) };
}
/** Iterate until lift and drag stop changing (relative change over the last ten iterations below tol) or the cap is reached. */
function settle(s, cap, tol = 5e-4) {
  const h = s.history, n0 = h.cl.length; let change = Infinity;
  while (h.cl.length - n0 < cap && !s.diverged) {
    s.iterate(Math.min(5, cap - (h.cl.length - n0))); const n = h.cl.length;
    if (n - n0 >= 15) { change = Math.max(Math.abs(h.cl[n - 1] - h.cl[n - 11]) / Math.max(0.2, Math.abs(h.cl[n - 1])), (0.5 * Math.abs(h.cd[n - 1] - h.cd[n - 11])) / Math.max(1e-4, Math.abs(h.cd[n - 1]))); if (change < tol) break; }
  }
  return { change, its: h.cl.length - n0, converged: change < 2e-3 && !s.diverged };
}
const sepPoint = (sf, upper) => { const n = sf.x.length, h = n >> 1; let x = 1; for (let m = upper ? h : h - 1; upper ? m < n : m >= 0; m += upper ? 1 : -1) if (sf.cf[m] < 0 && sf.x[m] > 0.03) { x = sf.x[m]; break; } return x; };
const polarCache = new Map();
/**
 * Section polar by body-fitted RANS: the angles of attack [deg, ascending] are run in sequence, each warm-started from the one
 * before. Results are cached on section, Reynolds number (three significant figures) and grid, so wings sharing a section reuse them.
 */
export function ransPolar(af, Re, alphas, G, cap = 40, progress = null) {
  const ReK = Number(Re.toPrecision(3)), key = [af.name, af.t.toFixed(4), ReK, alphas.join(','), G.wfn, G.nSurf, G.nNormal, G.nWake, G.radius, cap].join('|');
  if (polarCache.has(key)) return { ...polarCache.get(key), cached: true };
  const cla = 2 * PI * (1 + 0.77 * af.t), a0 = thinAirfoil(af).alphaL0, { grid, s } = sectionSolver(af, ReK, N.rad(alphas[0]), G, 0.92 * cla * (N.rad(alphas[0]) - a0));
  const P = { alpha: alphas.slice(), cl: [], cd: [], cdp: [], cdf: [], cm: [], conv: [], xsep: [], yplus: [], its: 0, cells: grid.ni * grid.nj, Re: ReK, quality: grid.quality };
  alphas.forEach((a, k) => {
    if (k) s.setAlpha(N.rad(a));
    const st = settle(s, cap, 1e-3), h = s.history, n = h.cl.length, F = s.forces(), sf = s.surface(), m = st.converged ? 1 : Math.min(10, st.its);
    // an unsteady (stalled) point is represented by the mean of its last iterations
    P.cl.push(N.mean(h.cl.slice(n - m))); const cd = N.mean(h.cd.slice(n - m)); P.cd.push(cd); P.cdf.push(F.cdf); P.cdp.push(cd - F.cdf); P.cm.push(F.cm); P.conv.push(st.converged); P.xsep.push(sepPoint(sf, P.cl[k] >= 0)); P.yplus.push(N.mean(sf.yplus)); P.its += st.its;
    progress?.((k + 1) / alphas.length, `Section polar: α = ${a}°`);
  });
  const km = N.argmax(P.cl); [P.aStall, P.clmax] = vertex(P.alpha, P.cl, km); P.kMax = km; P.stallInside = km < alphas.length - 1;
  if (polarCache.size > 24) polarCache.delete(polarCache.keys().next().value);
  polarCache.set(key, P);
  return P;
}
/** Polar look-up at α [deg]: linear between points, extrapolated along the first interval below the table, held above it. */
const polarAt = (P, key, a) => (a < P.alpha[0] ? P[key][0] + ((P[key][1] - P[key][0]) * (a - P.alpha[0])) / (P.alpha[1] - P.alpha[0]) : N.interp1(P.alpha, P[key], a));
const yPlusOk = (yp, wfn) => (wfn ? yp >= 15 && yp <= 300 : yp <= 5);
const nearField = (s, grid, names, win = [-0.6, 1.9, -0.75, 0.75]) => {
  // cell-centre lattice inside a window around the section, as a triangle mesh for contour plots
  const { ni, nj } = grid, id = new Int32Array(ni * nj).fill(-1), nodes = [], tris = [], F = names.map((n) => s.field(n)), vals = names.map(() => []);
  for (let q = 0; q < ni * nj; q++) if (s.xc[q] > win[0] && s.xc[q] < win[1] && s.yc[q] > win[2] && s.yc[q] < win[3]) { id[q] = nodes.length; nodes.push([s.xc[q], s.yc[q]]); F.forEach((f, m) => vals[m].push(f[q])); }
  for (let i = 0; i < ni - 1; i++) for (let j = 0; j < nj - 1; j++) { const a = id[i * nj + j], b = id[(i + 1) * nj + j], c = id[(i + 1) * nj + j + 1], d = id[i * nj + j + 1]; if (a >= 0 && b >= 0 && c >= 0 && d >= 0) tris.push([a, b, c], [a, c, d]); }
  return { nodes, tris, vals };
};
const gridLines = (grid, win = [-0.25, 1.3, -0.35, 0.35]) => {
  const { ni, nj, x, y } = grid, N1 = nj + 1, X = [], Y = [], inW = (q) => x[q] > win[0] && x[q] < win[1] && y[q] > win[2] && y[q] < win[3];
  for (let i = 0; i <= ni; i += 1) { let on = false; for (let j = 0; j <= nj; j++) { const q = i * N1 + j; if (inW(q)) { X.push(x[q]); Y.push(y[q]); on = true; } else if (on) break; } if (on) { X.push(NaN); Y.push(NaN); } }
  for (let j = 0; j <= nj; j += 2) { let on = false; for (let i = 0; i <= ni; i++) { const q = i * N1 + j; if (inW(q)) { X.push(x[q]); Y.push(y[q]); on = true; } else if (on) { X.push(NaN); Y.push(NaN); on = false; } } if (on) { X.push(NaN); Y.push(NaN); } }
  return { X, Y };
};
const TMR_SRC = 'NASA Turbulence Modeling Resource, 2D NACA 0012 airfoil validation, SA model: CFL3D on the 897 × 257 grid, M = 0.15, Re = 6 million, fully turbulent (https://tmbwg.github.io/turbmodels/naca0012_val_sa.html, formerly turbmodels.larc.nasa.gov; read 2026-10-09)';

// chord for Re = 6 million at sea level and 51 m/s (M = 0.15); the 'accurate' grid setting
const TMR_IN = { airfoil: '0012', tc: 0, chord: 1.7185, sweep_deg: 0, V: 51, alt_m: 0, dISA: 0, wall: WALLS[1], nSurf: 128, nNormal: 64, nWake: 24, radius: 30, max_iter: 600 };
// the same condition at 10 m/s, so that the Prandtl–Glauert factor is 1.0004 and the incompressible SU2 solution is compared like for like
const SU2_IN = { ...TMR_IN, chord: 8.764, V: 10 };
const SU2_SRC = 'SU2 8.5.0 (official linux64 binary), incompressible RANS with Spalart–Allmaras, on the C-grid and configuration written by this project\'s High-fidelity bridge (medium resolution: 28 416 cells, y⁺ ≈ 1, far field 40 chords), 1 600 iterations, lift steady to 1e-4; run 2026-10-09. SU2 lift on this grid is 2–3.5% below the present solver and the NASA reference trend';
const rans2d = {
  id: 'rans2d', title: '2-D section: body-fitted RANS (Spalart–Allmaras)', fidelity: 'numerical',
  summary: 'Reynolds-averaged Navier–Stokes solution of the flow around the wing section on a body-fitted C-grid with a boundary-layer mesh: lift, pressure and friction drag, moment, surface pressure and skin friction, separation, the wall-layer profile and the flow field, compared with the panel and boundary-layer method.',
  equations: ['Reynolds-averaged Navier–Stokes equations', 'Navier–Stokes equations', 'Continuity equation', 'Conservation of momentum equation', 'Spalart–Allmaras', 'Wall-resolved or wall-modelled turbulence simulations', 'Finite-volume methods', 'Prandtl–Glauert compressibility correction'],
  inputs: [
    { key: 'airfoil', label: 'NACA section', type: 'text', default: '2412', group: 'Geometry', help: '4-digit (2412, 0012, 4415) or 5-digit (23012) designation' },
    { key: 'tc', label: 'Thickness ratio override', unit: '-', default: 0, min: 0, max: 0.4, group: 'Geometry', help: '0 uses the thickness in the designation' },
    { key: 'chord', label: 'Chord', unit: 'm', default: 1.5, min: 0.005, group: 'Geometry' },
    { key: 'sweep_deg', label: 'Sweep (simple sweep theory)', unit: 'deg', default: 0, min: 0, max: 60, group: 'Geometry', help: 'The section sees V·cos Λ and M·cos Λ' },
    { key: 'alpha_deg', label: 'Angle of attack', unit: 'deg', default: 3, min: -10, max: 22, group: 'Flow' },
    { key: 'V', label: 'Freestream speed', unit: 'm/s', default: 70, min: 1, group: 'Flow' },
    ...ATM,
    ...RANS_GRID,
    { key: 'max_iter', label: 'Maximum iterations', unit: '', default: 70, min: 10, max: 5000, step: 1, discrete: true, group: 'Numerics', help: 'The run stops earlier once lift and drag are steady. 70 suits the default grid; use 300–600 on fine grids' },
    { key: 'x_profile', label: 'Station of the wall-layer profile', unit: 'x/c', default: 0.5, min: 0.05, max: 0.98, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const base = { alt_m: c.atm.alt_m, dISA: c.atm.dISA_K };
    if (c.wing.S_m2 > 0) return { ...base, airfoil: c.wing.airfoil, tc: c.wing.tc, chord: d.mac || undefined, sweep_deg: Math.abs(c.wing.sweep_deg), alpha_deg: c.flight.alpha_deg, V: c.flight.V_ms };
    if (!(c.rotor.R_m > 0)) return base;
    return { ...base, airfoil: '0012', tc: 0.12, chord: c.rotor.chord_m, sweep_deg: 0, V: 0.75 * d.v_tip, alpha_deg: 4 };
  },
  run(i, ctx) {
    const a = isa(i.alt_m, i.dISA), cosL = Math.cos(N.rad(i.sweep_deg)), Vn = i.V * cosL, M = Vn / a.a, Mc = Math.min(M, 0.9), Re = (Vn * i.chord) / a.nu, beta = Math.sqrt(1 - Mc * Mc);
    const af = parseNaca(i.airfoil, i.tc), G = gridOf(i), al = N.rad(i.alpha_deg), cap = Math.round(N.clamp(i.max_iter, 10, 5000));
    // panel and integral boundary-layer solution: initial circulation and the comparison
    const nd = nacaNodes(af, 120), Pn = panelSolve(nd.X, nd.Y), c0 = panelAt(Pn, 0), c1 = panelAt(Pn, N.rad(1)), cla = (c1.cl - c0.cl) / N.rad(1), bo = { Re: Math.max(Re, 1e4), xtU: 1, xtL: 1, Hsep: 2.2, kd: 1, cla }, ibl = sectionPoint(Pn, al, bo), iblT = sectionPoint(Pn, al, { ...bo, xtU: 0.03, xtL: 0.03 });
    const { grid, s } = sectionSolver(af, Re, al, G, 0.93 * ibl.inv.cl), h = s.history;
    ctx?.progress?.(0.05, 'Grid generated');
    let st = { change: Infinity, its: 0, converged: false };
    for (let done = 0; done < cap && !st.converged && !s.diverged;) { const n = Math.min(25, cap - done); st = settle(s, n); done += st.its; ctx?.progress?.(0.05 + (0.9 * done) / cap, `${done} iterations, residual down ${s.residualDrop.toFixed(1)} orders`); if (st.its < n) break; }
    const F = s.forces(), sf = s.surface(), nW = sf.x.length, half = nW >> 1, ypM = N.mean(sf.yplus), ypX = N.amax(sf.yplus), up = F.cl >= 0, xsep = sepPoint(sf, up), pr = s.profile(N.clamp(i.x_profile, 0.05, 0.98), up) || { yplus: [1], uplus: [1], utau: 0 };
    const cpC = sf.cp.map((v) => cpKT(v, Mc)), cl = F.cl / beta, cm = F.cm / beta, drop = s.residualDrop, q = grid.quality, ypGood = yPlusOk(ypM, G.wfn), conv = st.converged && !s.diverged;
    const warnings = [];
    if (!af.ok) warnings.push(`"${i.airfoil}" is not a supported NACA 4-digit or non-reflex 5-digit designation; NACA 0012 thickness form was used.`);
    if (s.diverged) warnings.push('The iteration diverged. Use a finer grid or the other wall treatment.');
    else if (!conv) warnings.push(`Not fully converged — raise the iterations: lift and drag still change by ${(100 * st.change).toPrecision(2)}% over the last ten of ${h.cl.length} iterations (residual down ${drop.toFixed(1)} orders). ${xsep < 0.9 ? 'The flow is separated; a steady RANS solution may not exist here and the values are means of an oscillating iteration.' : ''}`);
    if (!ypGood) warnings.push(G.wfn ? `Mean y⁺ of the first cells is ${ypM.toFixed(0)}: the log-law wall function needs roughly 20–200. ${ypM < 15 ? 'At this low Reynolds number choose the wall-resolved treatment.' : 'Use more wall-normal cells.'}` : `Mean y⁺ of the first cells is ${ypM.toFixed(1)}; a wall-resolved solution needs about 1 (at most 5).`);
    if (G.nSurf < 96) warnings.push(`Fast-look grid (${grid.ni} × ${grid.nj} = ${grid.ni * grid.nj} cells). Against the NASA reference solutions for the NACA 0012 at Re = 6 million grids of this size give drag about 7% high at zero lift and 7–12% high at α = 10°, with lift 3–4% low; near stall the drag error exceeds 25%. For numbers to quote use 128 cells around the section with 64 wall-normal cells, wall-resolved, and 300–600 iterations (drag within about 1–3% up to 10°), and run the mesh-convergence study.`);
    if (G.wfn && xsep < 0.95) warnings.push('Wall functions assume an attached equilibrium boundary layer: with separation present the wall-resolved treatment is more trustworthy.');
    if (M > 0.3) warnings.push(`Section Mach number ${M.toFixed(2)}: the flow solution is incompressible. Lift and moment are multiplied by the Prandtl–Glauert factor 1/√(1 − M²) = ${(1 / beta).toFixed(3)} and Cp is corrected with the Kármán–Tsien rule; drag is not corrected. ${M > 0.6 ? 'Above about Mach 0.6 shocks appear on the section and this correction is not valid: use a compressible solver (High-fidelity bridge).' : ''}`);
    if (Re < 2e5) warnings.push(`Chord Reynolds number ${Re.toExponential(2)}: the Spalart–Allmaras model is run fully turbulent and has no transition model, so laminar separation bubbles and laminar drag buckets are not represented.`);
    if (i.sweep_deg > 0) warnings.push('Simple sweep theory: the section is solved in the flow normal to the sweep line; coefficients refer to that normal dynamic pressure.');
    const xs = sf.x, lo = (A) => A.slice(0, half).reverse(), hi = (A) => A.slice(half), [itT, resT] = thin(N.range(h.res.length, (k) => k + 1), h.res.map((v) => Math.max(v, 1e-16)), 300), [, clT] = thin(h.res, h.cl, 300), [, cdT] = thin(h.res, h.cd.map((v) => 100 * v), 300);
    const nf = nearField(s, grid, ['speed', 'p', 'nut']), gl = gridLines(grid), ov = [{ name: af.name, x: nd.X, y: nd.Y }], k0 = Math.min(5, itT.length - 1);
    const logY = pr.yplus.filter((v) => v > 0.2 && v < 3000), logU = logY.map((v) => (v < 11.06 ? v : Math.log(v) / 0.41 + 5));
    return {
      kpis: [
        { key: 'cl_rans', label: 'Section lift coefficient', value: cl, unit: '-', status: conv ? 'ok' : 'warn', note: M > 0.3 ? 'incompressible solution × Prandtl–Glauert factor' : '' },
        { key: 'cd_rans', label: 'Section drag coefficient', value: F.cd, unit: '-', status: conv && ypGood ? 'ok' : 'warn', note: `${(1e4 * F.cd).toFixed(1)} counts, fully turbulent` },
        { key: 'cd_pressure_rans', label: 'Pressure (form) drag', value: F.cdp, unit: '-' },
        { key: 'cd_friction_rans', label: 'Friction drag', value: F.cdf, unit: '-' },
        { key: 'cm_c4_rans', label: 'Moment about quarter chord', value: cm, unit: '-', note: 'nose-up positive' },
        { key: 'ld_rans', label: 'Section lift-to-drag ratio', value: cl / F.cd, unit: '-' },
        { key: 'cp_min_rans', label: 'Minimum pressure coefficient', value: N.amin(cpC), unit: '-' },
        { key: 'x_sep_rans', label: 'Separation, suction side', value: xsep, unit: 'x/c', status: xsep >= 0.98 ? 'ok' : xsep > 0.8 ? 'warn' : 'bad', note: '1 = attached to the trailing edge (first point of negative skin friction)' },
        { key: 'y_plus_mean', label: 'Mean y⁺ of the first cells', value: ypM, unit: '-', status: ypGood ? 'ok' : 'warn', note: G.wfn ? 'wall functions: 20–200 wanted' : 'wall-resolved: about 1 wanted' },
        { key: 'y_plus_max', label: 'Maximum y⁺ of the first cells', value: ypX, unit: '-' },
        { key: 'residual_drop', label: 'Residual reduction', value: drop, unit: 'orders', status: conv ? 'ok' : 'warn' },
        { key: 'force_change_pct', label: 'Change of lift and drag over the last ten iterations', value: Number.isFinite(st.change) ? 100 * st.change : 100, unit: '%', status: conv ? 'ok' : 'warn', note: 'converged below 0.2%' },
        { key: 'iterations', label: 'Iterations', value: h.cl.length, unit: '' },
        { key: 'cl_panel_bl', label: 'Lift, panel + boundary layer', value: ibl.cl / beta, unit: '-', note: 'existing airfoil analysis, free transition' },
        { key: 'cd_panel_bl', label: 'Drag, panel + boundary layer (free transition)', value: ibl.cd, unit: '-', note: 'Squire–Young with Michel transition: lower than a fully turbulent solution when laminar flow is present' },
        { key: 'cd_panel_bl_turb', label: 'Drag, panel + boundary layer (transition at 3% chord)', value: iblT.cd, unit: '-', note: 'the like-for-like comparison with the fully turbulent RANS solution' },
        { key: 'Re_chord', label: 'Chord Reynolds number', value: Re, unit: '-' },
        { key: 'mach_section', label: 'Section Mach number', value: M, unit: '-' },
        { key: 'cells', label: 'Grid cells', value: grid.ni * grid.nj, unit: '', note: `${grid.ni} × ${grid.nj}` },
        { key: 'min_orthogonality_deg', label: 'Smallest grid-line angle', value: q.minOrthogonality_deg, unit: 'deg', status: q.minOrthogonality_deg > 30 ? 'ok' : 'warn', note: `${q.minOrthogonalityWall_deg.toFixed(0)}° in the four layers next to the wall` },
        { key: 'growth_bl', label: 'Wall-normal growth ratio in the boundary layer', value: grid.growth, unit: '-', status: grid.growth <= 1.3 ? 'ok' : 'warn', note: `largest anywhere ${q.maxGrowth.toFixed(2)}; 1.2–1.3 wanted in the boundary layer` },
        { key: 'aspect_ratio_max', label: 'Largest cell aspect ratio', value: q.maxAspectRatio, unit: '-' },
      ],
      plots: [
        { type: 'line', title: `Pressure distribution, ${af.name} at α = ${i.alpha_deg}°`, xlabel: 'x/c [-]', ylabel: '−Cp [-]', series: [{ name: 'RANS, upper surface', x: hi(xs), y: hi(cpC).map((v) => -v) }, { name: 'RANS, lower surface', x: lo(xs), y: lo(cpC).map((v) => -v) }, { name: 'Inviscid panel method', x: nd.X, y: ibl.inv.cp.map((v) => -cpKT(v, Mc)), style: 'dash' }] },
        { type: 'line', title: 'Skin-friction coefficient', xlabel: 'x/c [-]', ylabel: 'cf [-]', series: [{ name: 'RANS, upper surface', x: hi(xs), y: hi(sf.cf) }, { name: 'RANS, lower surface', x: lo(xs), y: lo(sf.cf) }, { name: 'Integral method, upper (free transition)', x: ibl.sf.up.x.slice(1), y: ibl.bu.cf.slice(1), style: 'dash' }], annotations: [{ y: 0, label: 'Separation below this line' }] },
        { type: 'line', title: `First-cell y⁺ (target ${G.yplus})`, xlabel: 'x/c [-]', ylabel: 'y⁺ [-]', series: [{ name: 'Upper surface', x: hi(xs), y: hi(sf.yplus) }, { name: 'Lower surface', x: lo(xs), y: lo(sf.yplus) }] },
        { type: 'line', title: `Wall-layer profile at x/c = ${(pr.x ?? i.x_profile).toFixed(2)}, ${up ? 'upper' : 'lower'} surface`, xlabel: 'y⁺ [-]', ylabel: 'u⁺ [-]', xlog: true, series: [{ name: 'RANS', x: pr.yplus, y: pr.uplus, style: 'line+points' }, { name: 'u⁺ = y⁺ and u⁺ = ln(y⁺)/0.41 + 5.0', x: logY, y: logU, style: 'dash' }] },
        { type: 'line', title: 'Residual history', xlabel: 'Iteration [-]', ylabel: 'Residual / initial residual [-]', ylog: true, series: [{ name: 'Mean-flow residual', x: itT, y: resT }] },
        { type: 'line', title: 'Force history', xlabel: 'Iteration [-]', ylabel: 'cl and 100·cd [-]', series: [{ name: 'cl', x: itT.slice(k0), y: clT.slice(k0) }, { name: '100 · cd', x: itT.slice(k0), y: cdT.slice(k0) }] },
        { type: 'line', title: `C-grid near the section (${grid.ni} × ${grid.nj} cells, every second wall-parallel line)`, xlabel: 'x/c [-]', ylabel: 'y/c [-]', equalAspect: true, series: [{ name: 'Grid lines', x: gl.X, y: gl.Y }, { name: af.name, x: nd.X, y: nd.Y }] },
        { type: 'tri', title: 'Velocity magnitude', xlabel: 'x/c [-]', ylabel: 'y/c [-]', zlabel: '|V|/V∞ [-]', nodes: nf.nodes, tris: nf.tris, values: nf.vals[0], equalAspect: true, overlay: ov },
        { type: 'tri', title: 'Pressure coefficient (incompressible)', xlabel: 'x/c [-]', ylabel: 'y/c [-]', zlabel: 'Cp [-]', nodes: nf.nodes, tris: nf.tris, values: nf.vals[1].map((v) => N.clamp(2 * v, -3, 1)), equalAspect: true, diverging: true, overlay: ov },
        { type: 'tri', title: 'Eddy-viscosity ratio', xlabel: 'x/c [-]', ylabel: 'y/c [-]', zlabel: 'ν_t/ν [-]', nodes: nf.nodes, tris: nf.tris, values: nf.vals[2], equalAspect: true, overlay: ov },
      ],
      tables: [{ title: 'Body-fitted RANS against the panel and boundary-layer method', columns: ['Quantity', 'RANS (fully turbulent)', 'Panel + boundary layer, free transition', 'Panel + boundary layer, transition at 3%'], rows: [['cl', +cl.toFixed(4), +(ibl.cl / beta).toFixed(4), +(iblT.cl / beta).toFixed(4)], ['cd', +F.cd.toFixed(5), +ibl.cd.toFixed(5), +iblT.cd.toFixed(5)], ['cm about c/4', +cm.toFixed(4), +(ibl.inv.cm4 / beta).toFixed(4), +(iblT.inv.cm4 / beta).toFixed(4)], ['Separation x/c', +xsep.toFixed(3), +ibl.f.toFixed(3), +iblT.f.toFixed(3)]] }],
      outputs: { rans2d_converged: conv ? 1 : 0, cd_counts_rans: 1e4 * F.cd },
      warnings,
      models: ['Steady incompressible Reynolds-averaged Navier–Stokes equations by artificial compressibility, cell-centred finite volume on a body-fitted C-grid', 'Spalart–Allmaras one-equation model, standard constants, no trip term (fully turbulent), negative-ν̃ continuation; free-stream ν̃ = 3ν', G.wfn ? 'Log-law wall function in the first cell (κ = 0.41, B = 5.0), ν̃ = κ·u_τ·y there' : 'No-slip wall resolved to the viscous sublayer', 'Roe flux-difference splitting with third-order upwind-biased (QUICK, κ = 1/2) reconstruction weighted for the grid stretching; first-order upwind turbulence convection', 'Implicit line-relaxation time march with local time steps and residual-driven CFL', 'Far field: characteristic boundary with point-vortex correction', 'Prandtl–Glauert (loads) and Kármán–Tsien (Cp) compressibility corrections', 'Comparison: linear-vortex panel method with Thwaites / Michel / Head boundary layer'],
      assumptions: ['Steady, two-dimensional, incompressible, fully turbulent flow; no transition model', 'Closed (sharp) trailing edge as in the NACA thickness form used throughout this suite', 'Algebraic C-grid: wall-normal lines blended into rays, wake cut along the chord line', i.sweep_deg > 0 ? 'Simple sweep theory for the swept-wing section' : 'Unswept section'],
    };
  },
  convergence: { param: 'nSurf', label: 'Cells around the section', levels: [40, 56, 80, 112], metric: 'cd_rans' },
  calibration: { params: [{ key: 'alpha_deg', min: -5, max: 15 }], sweep: 'alpha_deg', target: 'cl_rans', note: 'The turbulence model has no free constants here. Supply a measured lift curve or drag polar at the same Reynolds number to quantify the model-form error (transition, stall); a constant angle offset absorbs tunnel flow-angularity.' },
  verify() {
    const lin = (a, b, n) => N.range(n + 1, (k) => a + ((b - a) * k) / n);
    // 1. free stream on a randomly distorted grid
    const rnd = N.rng(7), gF = boxGrid({ xf: lin(0, 1, 16), yf: lin(0, 1, 12), distort: () => [0.02 * (rnd() - 0.5), 0.025 * (rnd() - 0.5)] }), sF = createRans2d(gF, { Re: 100, alpha: 0.3, vortex: false, noDamp: true }), rF = sF.residualOnly();
    // 2. observed order on the Kovasznay flow (exact Navier–Stokes solution), distorted grids
    const lam = 20 - Math.sqrt(400 + 4 * PI * PI), kov = (x, y) => { const e = Math.exp(lam * x); return [0.5 * (1 - e * e), 1 - e * Math.cos(2 * PI * y), (lam / (2 * PI)) * e * Math.sin(2 * PI * y)]; };
    const err = [16, 32].map((n) => { const g = boxGrid({ xf: lin(-0.5, 1, n), yf: lin(-0.5, 0.5, n), distort: (x, y) => [(0.15 / n) * Math.sin(7 * x + 3 * y), (0.15 / n) * Math.cos(5 * x - 4 * y)] }), s = createRans2d(g, { Re: 40, turbulent: false, vortex: false, farState: kov, init: kov, beta: 4, kappa: 1 / 3 }); s.iterate(250, 1e-9); let e = 0, v = 0; for (let q = 0; q < n * n; q++) { const d = s.state.U[q] - kov(s.xc[q], s.yc[q])[1]; e += d * d * s.vol[q]; v += s.vol[q]; } return Math.sqrt(e / v); });
    // 3. laminar flat plate against Blasius
    const ReL = 2e5, sL = createRans2d(plateGrid({ ni: 48, nj: 32, h1: 2e-4, height: 0.5, rBL: 1.15 }), { Re: ReL, turbulent: false }); sL.iterate(300, 1e-7);
    const fL = sL.surface(), kL = fL.x.findIndex((x) => x > 0.5);
    // 4. turbulent flat plate (wall-resolved and wall-function grids) against White's correlation and the log law
    const ReT = 5e6, white = (x) => 0.455 / Math.log(0.06 * ReT * x) ** 2, at = (f, x) => f.cf[f.x.findIndex((v) => v > x)] / white(f.x[f.x.findIndex((v) => v > x)]);
    const sT = createRans2d(plateGrid({ ni: 40, nj: 36, h1: wallSpacing(ReT, 1), height: 0.5 }), { Re: ReT }); sT.iterate(260, 1e-6);
    const sW = createRans2d(plateGrid({ ni: 40, nj: 20, h1: wallSpacing(ReT, 60), height: 0.5 }), { Re: ReT, wallFunction: true }); sW.iterate(150, 1e-6);
    const pT = sT.profile(0.9), jL = pT.yplus.findIndex((v) => v > 80), uLog = Math.log(pT.yplus[jL]) / 0.41 + 5;
    // 5. symmetric section at zero incidence
    const z = sectionSolver(parseNaca('0012'), 6e6, 0, { wfn: true, yplus: 50, nSurf: 48, nNormal: 20, nWake: 10, radius: 30 }).s; z.iterate(60);
    return [
      N.check('Free stream preserved on a distorted grid', 1 + rF, 1, 1e-11, 'Uniform flow is an exact solution; the finite-volume metrics are closed'),
      N.check('Observed order of accuracy, Kovasznay flow on distorted grids', Math.log2(err[0] / err[1]), 1.75, 0.2, 'Exact Navier–Stokes solution (Kovasznay 1948), Re = 40, grids 16² and 32²; design order 2 in the interior, the weakly imposed first-order boundaries lower the observed order to about 1.5–1.9'),
      N.check('Laminar flat plate: cf at Re_x = 1.0e5', fL.cf[kL], 0.664 / Math.sqrt(ReL * fL.x[kL]), 0.02, 'Blasius'),
      N.check('Laminar flat plate: drag coefficient', sL.forces().cd, 1.328 / Math.sqrt(ReL), 0.02, 'Blasius, Re_L = 2e5'),
      N.check('Turbulent flat plate, wall-resolved: cf at Re_x = 4.5e6', at(sT.surface(), 0.9), 1, 0.08, 'White: cf = 0.455/ln²(0.06 Re_x). The Spalart–Allmaras model sits about 5% below this correlation at this Reynolds number'),
      N.check('Turbulent flat plate, wall functions: cf at Re_x = 4.5e6', at(sW.surface(), 0.9), 1, 0.08, 'White: cf = 0.455/ln²(0.06 Re_x)'),
      N.check('Turbulent flat plate: u⁺ in the log layer', pT.uplus[jL], uLog, 0.04, 'Law of the wall u⁺ = ln(y⁺)/0.41 + 5.0'),
      N.check('Symmetric section at zero incidence carries no lift', 1 + z.forces().cl, 1, 2e-3, 'Symmetry'),
    ];
  },
  validation: [
    { name: 'NACA 0012 drag at Re = 6 million, α = 0°, 10°, 15° (fully turbulent Spalart–Allmaras)', source: TMR_SRC + '. The reference section is the sharp-trailing-edge NACA 0012 rescaled to unit chord, about 0.9% thinner than the closed-trailing-edge form used here; the reference is a compressible solution at M = 0.15 with the far field 500 chords away', inputs: { ...TMR_IN }, sweep: { key: 'alpha_deg', values: [0, 10, 15] }, target: 'cd_rans', observed: [0.00819, 0.01231, 0.02124], tol_pct: 12 },
    { name: 'NACA 0012 lift at Re = 6 million, α = 10°, 15° (fully turbulent Spalart–Allmaras)', source: TMR_SRC, inputs: { ...TMR_IN }, sweep: { key: 'alpha_deg', values: [10, 15] }, target: 'cl_rans', observed: [1.0909, 1.5461], tol_pct: 4 },
    { name: 'NACA 0012 lift at Re = 6 million, α = 4°, 8°: cross-check with SU2', source: SU2_SRC, inputs: { ...SU2_IN }, sweep: { key: 'alpha_deg', values: [4, 8] }, target: 'cl_rans', observed: [0.4319, 0.8386], tol_pct: 5 },
    { name: 'NACA 2412 lift at Re = 6 million, α = 0°, 4°, 8°: cross-check with SU2', source: SU2_SRC, inputs: { ...SU2_IN, airfoil: '2412' }, sweep: { key: 'alpha_deg', values: [0, 4, 8] }, target: 'cl_rans', observed: [0.2205, 0.6497, 1.0487], tol_pct: 5 },
    { name: 'NACA 0012 and 2412 drag at zero incidence: cross-check with SU2', source: SU2_SRC + '. Drag is compared at α = 0° only: on this medium grid SU2 itself is 7.6% above the NASA reference for the NACA 0012 at α = 0° and its drag error grows quickly with lift (its fine-grid value at α = 4° is 12% below its medium-grid value), so the lifting-case SU2 drag is not a reference', inputs: { ...SU2_IN }, sweep: { key: 'airfoil', values: ['0012', '2412'] }, target: 'cd_rans', observed: [0.008816, 0.009005], tol_pct: 10 },
  ],
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!o.rans2d_converged) out.push({ severity: 'warn', title: 'Solution not converged', detail: `Forces still change by ${o.force_change_pct.toPrecision(2)}% per ten iterations after ${o.iterations} iterations.`, action: 'Raise the maximum iterations (300–600 on fine grids). If the flow is separated, treat the value as the mean of an unsteady flow and confirm with the wall-resolved treatment.', basis: 'Iterative convergence of the steady RANS equations' });
    if (o.x_sep_rans < 0.9) out.push({ severity: o.x_sep_rans < 0.6 ? 'critical' : 'warn', title: 'Separated flow on the suction side', detail: `Skin friction changes sign at x/c = ${o.x_sep_rans.toFixed(2)}.`, action: 'Lower the angle of attack or accept the loss; near and beyond maximum lift a steady one-equation RANS solution is a trend, not a prediction. Use the wing analysis below for the stall margin.', basis: 'Sign of the wall shear stress' });
    if (i.nSurf < 96) out.push({ severity: 'advise', title: 'Refine before quoting drag', detail: `${o.cells} cells is the one-second preview grid.`, action: 'Set 128 cells around the section, 64 wall-normal cells, wall-resolved, 400 iterations (about 15–30 s), and check with the mesh-convergence tab.', basis: 'Grid study against the NASA Turbulence Modeling Resource NACA 0012 case' });
    if (o.cd_panel_bl < 0.8 * o.cd_rans) out.push({ severity: 'info', title: 'Laminar flow would save drag', detail: `Fully turbulent: ${o.cd_counts_rans.toFixed(0)} counts. With natural transition the integral method gives ${(1e4 * o.cd_panel_bl).toFixed(0)} counts.`, action: 'The difference is the value of keeping the leading edge clean and smooth; every count of section drag is fuel and CO₂.', basis: 'Comparison of fully turbulent RANS with the free-transition integral boundary layer' });
    return out;
  },
};

/**
 * Viscous-coupled wing: vortex lattice with sectional decambering (α-shift method of van Dam / Mukherjee–Gopalarathnam).
 * Every spanwise strip is rotated by δ so that its lattice lift equals the viscous section lift at its effective angle
 * α_e = cl/a0 + α_L0 + δ. sec(j, α_e [rad]) returns the strip's { cl, cd, cm }. Returns one record per wing angle [rad].
 */
export function viscousWing(g, nS, nC, beta, af, a0, sec, alphas) {
  const L = vlm(g, nS, nC, beta, (xf) => af.camber(xf)[1], true), n2 = L.n2, aL0 = thinAirfoil(af).alphaL0, ch = L.ym.map((y) => L.chord(y));
  // strip circulation produced by a unit decambering of the symmetric strip pair m
  const D = N.range(nS, (m) => { const G = L.solve(L.pn.map((p) => (p.js === nS + m || p.js === nS - 1 - m ? 1 : 0))), st = new Array(n2).fill(0); L.pn.forEach((p, k) => { st[p.js] += G[k]; }); return st; });
  const dl = new Array(nS).fill(0), out = [], clOf = sec.cl || ((j, e) => sec(j, e).cl);
  for (const al of alphas) {
    let strip = [], cl = [], ae = [], res = 1, it = 0;
    for (; it < 150; it++) {
      strip = L.a.strip.map((v, j) => v * al + L.z.strip[j]); for (let m = 0; m < nS; m++) if (dl[m]) for (let j = 0; j < n2; j++) strip[j] += D[m][j] * dl[m];
      cl = strip.map((v, j) => (2 * v) / ch[j]); res = 0; ae = [];
      for (let m = 0; m < nS; m++) { const j = nS + m, e = cl[j] / a0 + aL0 + dl[m], r = cl[j] - clOf(j, e); ae.push(e); res = Math.max(res, Math.abs(r)); dl[m] += (0.35 * r) / a0; }
      if (res < 2e-5) break;
    }
    const half = N.range(nS, (m) => nS + m), sc = half.map((j, m) => sec(j, ae[m])), S2 = L.S / 2, CL = N.sum(half.map((j) => cl[j] * ch[j] * L.dy[j])) / S2, CDp = N.sum(half.map((j, m) => sc[m].cd * ch[j] * L.dy[j])) / S2;
    const xr = g.cr / 4 + g.ymac * g.tanL, CM = N.sum(half.map((j, m) => (sc[m].cm * ch[j] - cl[j] * (g.cr / 4 + L.ym[j] * g.tanL - xr)) * ch[j] * L.dy[j])) / (S2 * g.mac);
    out.push({ alpha: al, CL, CDi: L.trefftz(strip), CDp, CM, cl: half.map((j) => cl[j]), cd: sc.map((v) => v.cd), ae, eta: half.map((j) => (2 * L.ym[j]) / g.b), strip: half.map((j) => strip[j]), residual: res, iterations: it });
  }
  return { L, pts: out };
}

const wingRans = {
  id: 'wingrans', title: '3-D wing with viscous sections: vortex lattice coupled to section RANS', fidelity: 'numerical',
  summary: 'Viscous lift, drag and stall of the wing at flight Reynolds number: a vortex lattice carries the three-dimensional induced flow and every spanwise strip takes its lift, drag and moment from a body-fitted RANS polar of the wing section at its effective angle of attack. Gives the wing polar, span loading, drag breakdown, maximum lift and where the stall starts.',
  equations: ['Reynolds-averaged Navier–Stokes equations', 'Spalart–Allmaras', 'Prandtl lifting-line theory', 'Lifting-surface integral equations', 'Kutta–Joukowski lift theorem', 'Vortex lattice models', 'Prandtl–Glauert compressibility correction'],
  applicable: fixedWing,
  inputs: [
    { key: 'S', label: 'Wing area', unit: 'm²', default: 16.2, min: 0.01, group: 'Geometry' },
    { key: 'b', label: 'Span', unit: 'm', default: 11, min: 0.05, group: 'Geometry' },
    { key: 'taper', label: 'Taper ratio', unit: '-', default: 0.6, min: 0.05, max: 1, group: 'Geometry' },
    { key: 'sweep_deg', label: 'Quarter-chord sweep', unit: 'deg', default: 0, min: -30, max: 45, group: 'Geometry' },
    { key: 'twist_deg', label: 'Tip twist (washout negative)', unit: 'deg', default: -2, min: -10, max: 5, group: 'Geometry' },
    { key: 'airfoil', label: 'NACA section', type: 'text', default: '2412', group: 'Geometry' },
    { key: 'tc', label: 'Thickness ratio override', unit: '-', default: 0, min: 0, max: 0.4, group: 'Geometry', help: '0 uses the thickness in the designation' },
    { key: 'alpha_deg', label: 'Wing angle of attack', unit: 'deg', default: 4, min: -6, max: 22, group: 'Flow', help: 'Root-chord incidence to the freestream' },
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 60, min: 1, group: 'Flow' },
    ...ATM,
    { key: 'CL_flight', label: 'Lift coefficient in flight (0 = use the angle of attack)', unit: '-', default: 0, min: 0, max: 3, group: 'Aircraft', help: 'Weight / (q·S): the published lift and drag are taken at this point of the polar' },
    { key: 'CD_other', label: 'Drag of everything except the wing', unit: '-', default: 0.012, min: 0, max: 0.2, group: 'Aircraft', help: 'Fuselage, tails, nacelles, excrescences, extra drag area and wave drag from the drag build-up; added to the wing drag for the aircraft polar' },
    { key: 'CD_wing_buildup', label: 'Wing profile drag in the drag build-up (reference)', unit: '-', default: 0, min: 0, max: 0.1, group: 'Aircraft', help: 'Flat-plate friction × form factor, shown beside the RANS value' },
    { key: 'alpha_table', label: 'Section angles of attack', type: 'text', default: '0, 6, 11, 14, 17, 20', group: 'Numerics', help: 'Angles [deg] of the section RANS polar, ascending. More points near stall sharpen the maximum lift: e.g. -4, 0, 4, 8, 10, 12, 14, 15, 16, 17, 18, 20' },
    { key: 'nRe', label: 'Reynolds numbers along the span', unit: '', default: 1, min: 1, max: 3, step: 1, discrete: true, group: 'Numerics', help: '1: one polar at the mean-chord Reynolds number, friction scaled to each strip. 2–3: separate polars between root and tip, interpolated in log Re' },
    ...RANS_GRID.map((f) => (f.key === 'nSurf' ? { ...f, default: 48 } : f.key === 'nNormal' ? { ...f, default: 20 } : f.key === 'nWake' ? { ...f, default: 10 } : f)),
    { key: 'iter_alpha', label: 'Iterations per section angle (cap)', unit: '', default: 120, min: 10, max: 2000, step: 1, discrete: true, group: 'Numerics', help: 'Each angle starts from the one before and stops once lift and drag are steady; attached flow needs 40–70 iterations on the default grid, 150–300 on a fine one' },
    { key: 'nSpan', label: 'Spanwise strips per half wing', unit: '', default: 12, min: 4, max: 40, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const bi = {}; for (const f of dragBuildup.inputs) bi[f.key] = f.default; const src = dragBuildup.defaults(c, up, d); for (const k of Object.keys(src)) if (src[k] !== undefined && Number.isFinite(src[k])) bi[k] = src[k];
    let other, wingBu, CLf;
    if (bi.S > 0 && bi.b > 0) { const r = buildup(bi), w = r.comp.find((q) => q.name === 'Wing'), wcd = w ? w.CD * (1 + bi.exc_pct / 100) : 0; other = r.CD0 - wcd + r.w.CDw; wingBu = wcd; CLf = r.CL; }
    return { S: c.wing.S_m2, b: c.wing.b_m, taper: c.wing.taper, sweep_deg: c.wing.sweep_deg, twist_deg: c.wing.twist_deg, airfoil: c.wing.airfoil, tc: c.wing.tc, alpha_deg: c.flight.alpha_deg + c.wing.incidence_deg, V: c.flight.V_ms, alt_m: c.atm.alt_m, dISA: c.atm.dISA_K, CD_other: other, CD_wing_buildup: wingBu, CL_flight: CLf };
  },
  run(i, ctx) {
    const a = isa(i.alt_m, i.dISA), M = i.V / a.a, g = wingGeom(i), af = parseNaca(i.airfoil, i.tc), G = gridOf(i), cosL = Math.cos(N.rad(i.sweep_deg)), Mn = M * cosL, bn = Math.sqrt(1 - Math.min(Mn, 0.9) ** 2), beta = Math.sqrt(1 - Math.min(M, 0.9) ** 2);
    const cbar = i.S / i.b, ReN = (c) => (i.V * cosL * c * cosL) / a.nu, nS = Math.round(N.clamp(i.nSpan, 4, 40));
    let tab = String(i.alpha_table).split(/[,;\s]+/).map(Number).filter(Number.isFinite); tab = [...new Set(tab)].sort((p, q) => p - q).filter((v) => v >= -12 && v <= 24); if (tab.length < 3) tab = [0, 6, 11, 14, 17, 20];
    const nRe = Math.round(N.clamp(i.nRe, 1, 3)), ct = g.cr * g.taper, cRe = nRe === 1 ? [g.mac] : nRe === 2 ? [ct, g.cr] : [ct, g.mac, g.cr], cap = Math.round(N.clamp(i.iter_alpha, 10, 2000));
    const pol = cRe.map((c, k) => ransPolar(af, ReN(c), tab, G, cap, (f, msg) => ctx?.progress?.((0.9 * (k + f)) / cRe.length, msg))), lr = pol.map((p) => Math.log(p.Re)), P0 = pol[nRe === 3 ? 1 : 0];
    // section data of a strip: simple sweep theory, Prandtl–Glauert on lift and moment, friction scaled to the strip Reynolds number when only one polar is used
    const a0 = (2 * PI * cosL) / bn, cs = [], look = (key, an, lRe) => { if (nRe === 1) return polarAt(P0, key, an); const v = pol.map((p) => polarAt(p, key, an)); return N.interp1(lr, v, lRe); };
    const sec = (j, ae) => {
      if (!cs[j]) { const c = g.cr * (1 - (1 - g.taper) * Math.abs(cs.eta[j])), Re = ReN(c); cs[j] = { lRe: Math.log(Re), kf: nRe === 1 ? cfTurb(Math.max(Re, 1e4)) / cfTurb(Math.max(P0.Re, 1e4)) : 1 }; }
      const an = N.deg(ae) / cosL, q = cs[j];
      return { cl: (look('cl', an, q.lRe) * cosL * cosL) / bn, cd: look('cdp', an, q.lRe) * cosL ** 3 + look('cdf', an, q.lRe) * q.kf, cm: (look('cm', an, q.lRe) * cosL * cosL) / bn };
    };
    sec.cl = (j, ae) => { if (!cs[j]) sec(j, ae); return (look('cl', N.deg(ae) / cosL, cs[j].lRe) * cosL * cosL) / bn; };
    cs.eta = N.range(2 * nS, (j) => -Math.cos((PI * (j + 0.5)) / (2 * nS)));
    const aW = N.range(27, (k) => k - 4), W = viscousWing(g, nS, 3, beta, af, a0, sec, [...aW.map(N.rad), N.rad(i.alpha_deg)]), pts = W.pts.slice(0, aW.length), here = W.pts[aW.length], L = W.L;
    ctx?.progress?.(0.97, 'Wing solved');
    const CLs = pts.map((p) => p.CL), km = N.argmax(CLs), [aSt, CLmax] = vertex(aW, CLs, km), pS = pts[km], secMaxN = (P0.clmax * cosL * cosL) / bn;
    // stall onset: first strip whose effective angle reaches the section stall angle as the wing angle increases
    const aSecStall = N.rad(P0.aStall) * cosL; let onset = null; for (const p of pts) { const m = p.ae.findIndex((e) => e >= aSecStall); if (m >= 0) { onset = { alpha: N.deg(p.alpha), eta: p.eta[N.argmax(p.ae)], CL: p.CL }; break; } }
    const etaSt = onset ? onset.eta : pS.eta[N.argmax(pS.cl)];
    const CDw = (p) => p.CDp + p.CDi, lin = pts.filter((p, k) => aW[k] >= -2 && aW[k] <= 4), CLa = (lin[lin.length - 1].CL - lin[0].CL) / (lin[lin.length - 1].alpha - lin[0].alpha), CLaV = L.a.CL;
    const up = pts.slice(0, km + 1), LD = up.map((p) => (p.CL > 0.05 ? p.CL / (CDw(p) + i.CD_other) : 0)), kL = N.argmax(LD), LDmax = LD[kL];
    // flight point on the polar (below maximum lift)
    const CLup = up.map((p) => p.CL), mono = CLup.every((v, k) => k === 0 || v > CLup[k - 1]), useCL = i.CL_flight > 0 && mono && i.CL_flight < CLmax && i.CL_flight > CLup[0], at = (f) => (useCL ? N.interp1(CLup, up.map(f), i.CL_flight) : f(here));
    const CLf = useCL ? i.CL_flight : here.CL, CDpf = at((p) => p.CDp), CDif = at((p) => p.CDi), CDf = CDpf + CDif + i.CD_other, aF = useCL ? N.interp1(CLup, up.map((p) => N.deg(p.alpha)), i.CL_flight) : i.alpha_deg, eV = CDif > 1e-9 ? (CLf * CLf) / (PI * L.AR * CDif) : 1;
    // sanity gate for publishing as the suite's primary viscous values
    const used = P0.alpha.map((v, k) => k <= P0.kMax), convAll = pol.every((p) => p.conv.every((c, k) => c || k > p.kMax)), ypM = N.mean(P0.yplus.filter((_, k) => used[k])), ypGood = yPlusOk(ypM, G.wfn), cplOk = up.every((p) => p.residual < 1e-3), stallIn = P0.stallInside && km < aW.length - 1, machOk = Mn <= 0.6, reOk = P0.Re >= 3e5;
    const gate = convAll && ypGood && cplOk && stallIn && machOk && reOk && CLmax > 0.3 && CLa > 0, why = [];
    if (!convAll) why.push('section polar points below maximum lift did not converge within the iteration cap');
    if (!ypGood) why.push(`mean first-cell y⁺ = ${ypM.toFixed(G.wfn ? 0 : 1)} is outside the range of the wall treatment`);
    if (!cplOk) why.push('the strip coupling did not converge');
    if (!stallIn) why.push('no lift maximum inside the table of section angles (extend the table)');
    if (!machOk) why.push(`normal Mach number ${Mn.toFixed(2)} is above 0.6, beyond an incompressible section solution with Prandtl–Glauert scaling`);
    if (!reOk) why.push(`section Reynolds number ${P0.Re.toExponential(2)} is below 3e5, where a fully turbulent model without transition is not representative`);
    const warnings = [];
    if (!af.ok) warnings.push(`"${i.airfoil}" is not a supported NACA designation; NACA 0012 thickness form was used.`);
    warnings.push(gate ? 'Sanity gate passed (section polar converged below maximum lift, y⁺ in range, coupling converged, subcritical): the lift, drag, maximum lift and best lift-to-drag ratio at the flight point are published as the suite\'s viscous values, in place of the handbook build-up values.' : `Sanity gate NOT passed: ${why.join('; ')}. The values below are shown for inspection only; the published CL, CD, CLmax and L/D of the suite remain those of the drag build-up.`);
    if (G.nSurf < 96) warnings.push(`Fast-look section grid (${P0.cells} cells, ${tab.length} angles). On this grid section drag is typically 5–20% high and maximum lift is resolved only to the spacing of the angle table. For design numbers use 96–128 cells around the section, 48–64 wall-normal cells wall-resolved, 150–300 iterations per angle and a table with 1° steps through the stall (one to a few minutes; the polar is cached, so later runs with the same section and Reynolds number are immediate).`);
    if (Math.abs(i.sweep_deg) > 1) warnings.push('Simple sweep theory: sections are solved in the flow normal to the quarter-chord line. Spanwise boundary-layer flow, which thickens the tip boundary layer of a swept wing and promotes tip stall, is not modelled; the stall of a swept wing is therefore optimistic.');
    if (M > 0.3) warnings.push(`Mach ${M.toFixed(2)}: the lattice uses the Prandtl–Glauert transformation and section lift is scaled by 1/√(1 − M_n²) = ${(1 / bn).toFixed(2)}; section drag is the incompressible value. Wave drag enters only through the build-up term.`);
    if (nRe === 1 && g.taper < 0.7) warnings.push(`One polar at the mean-chord Reynolds number ${P0.Re.toExponential(2)} is used for all strips (tip ${ReN(ct).toExponential(2)}); friction drag is scaled to the strip Reynolds number but maximum lift is not. Set "Reynolds numbers along the span" to 2 or 3 to compute tip and root polars.`);
    if (i.CL_flight > 0 && !useCL) warnings.push(`The flight lift coefficient ${i.CL_flight.toFixed(2)} is outside the attached part of the computed polar (maximum ${CLmax.toFixed(2)}); the values at the entered angle of attack are shown instead.`);
    warnings.push('Wing alone: no fuselage carry-over, nacelles, flaps or tail; planar wake. Post-stall points use steady RANS section data and are a trend only.');
    const cnt = (v) => +(1e4 * v).toFixed(1), eta = pS.eta, pF = useCL ? pts[N.argmin(pts.map((p) => Math.abs(p.CL - CLf)))] : here, CLv = (ad) => L.a.CL * N.rad(ad) + L.z.CL;
    const kpis = [
      { key: 'CL_wingrans', label: useCL ? 'Lift coefficient at the flight point' : 'Wing lift coefficient at this angle', value: CLf, unit: '-', status: CLf < 0.9 * CLmax ? 'ok' : 'warn', note: `α = ${aF.toFixed(2)}°` },
      { key: 'CD_wingrans', label: 'Aircraft drag coefficient (viscous wing + other components)', value: CDf, unit: '-', status: gate ? 'ok' : 'warn', note: `${cnt(CDf)} counts` },
      { key: 'CD_wing_profile', label: 'Wing profile drag (integrated section RANS)', value: CDpf, unit: '-', note: `${cnt(CDpf)} counts; drag build-up has ${cnt(i.CD_wing_buildup)}` },
      { key: 'CD_wing_induced', label: 'Induced drag (Trefftz plane)', value: CDif, unit: '-', note: `span efficiency ${eV.toFixed(3)}` },
      { key: 'LD_wingrans', label: 'Lift-to-drag ratio at this point', value: CLf / CDf, unit: '-' },
      { key: 'LD_max_wingrans', label: 'Maximum lift-to-drag ratio', value: LDmax, unit: '-', note: `at CL = ${up[kL].CL.toFixed(2)}` },
      { key: 'CLmax_wingrans', label: 'Clean wing maximum lift', value: CLmax, unit: '-', status: stallIn ? 'ok' : 'warn', note: `section maximum ${secMaxN.toFixed(2)} (normal to the sweep line: ${P0.clmax.toFixed(2)})` },
      { key: 'alpha_stall_wingrans_deg', label: 'Wing angle at maximum lift', value: aSt, unit: 'deg' },
      { key: 'eta_stall_wingrans', label: 'Spanwise station where stall starts', value: Math.abs(etaSt), unit: '2y/b', status: Math.abs(etaSt) < 0.7 ? 'ok' : 'warn', note: onset ? `first strip to reach the section stall angle, at α = ${onset.alpha.toFixed(0)}°` : 'strip with the highest lift at maximum lift' },
      { key: 'CLa_wingrans_per_rad', label: 'Viscous wing lift-curve slope', value: CLa, unit: '1/rad' },
      { key: 'CLa_vlm_per_rad', label: 'Inviscid vortex-lattice slope (reference)', value: CLaV, unit: '1/rad' },
      { key: 'CLa_ratio', label: 'Viscous / inviscid lift slope', value: CLa / CLaV, unit: '-', status: Math.abs(CLa / CLaV - 1) < 0.08 ? 'ok' : 'warn', note: 'viscous sections lose a few percent of slope to the boundary layer' },
      { key: 'CM_wingrans', label: 'Wing pitching moment about the quarter chord of the mean chord', value: at((p) => p.CM), unit: '-' },
      { key: 'gate_passed', label: 'Sanity gate for publishing', value: gate ? 1 : 0, unit: '', status: gate ? 'ok' : 'warn', note: gate ? 'passed' : why[0] || '' },
      { key: 'y_plus_mean', label: 'Mean first-cell y⁺ of the section solutions', value: ypM, unit: '-', status: ypGood ? 'ok' : 'warn' },
      { key: 'Re_section', label: 'Section Reynolds number (normal to the sweep line)', value: P0.Re, unit: '-' },
      { key: 'rans_iterations', label: 'RANS iterations in this run', value: pol.reduce((s2, p) => s2 + (p.cached ? 0 : p.its), 0), unit: '', note: pol.every((p) => p.cached) ? 'polar taken from the cache' : `${P0.cells} cells per section grid` },
    ];
    return {
      kpis,
      plots: [
        { type: 'line', title: 'Wing lift curve', xlabel: 'Angle of attack [deg]', ylabel: 'CL [-]', series: [{ name: 'Viscous sections (RANS)', x: aW, y: CLs }, { name: 'Inviscid vortex lattice', x: aW.filter((v) => v <= 14), y: aW.filter((v) => v <= 14).map(CLv), style: 'dash' }], annotations: [{ x: aSt, label: 'Maximum lift' }] },
        { type: 'line', title: 'Drag polar', xlabel: 'CD [-]', ylabel: 'CL [-]', series: [{ name: 'Aircraft (wing + other components)', x: up.map((p) => CDw(p) + i.CD_other), y: CLup }, { name: 'Wing alone', x: up.map(CDw), y: CLup, style: 'dash' }, { name: 'Wing profile drag only', x: up.map((p) => p.CDp), y: CLup, style: 'dash' }], annotations: [{ y: CLf, label: 'Flight point' }] },
        { type: 'line', title: 'Section lift along the span', xlabel: 'Span station 2y/b [-]', ylabel: 'cl [-]', series: [{ name: `Flight point, α = ${N.deg(pF.alpha).toFixed(1)}°`, x: eta, y: pF.cl }, { name: `Maximum lift, α = ${aW[km]}°`, x: eta, y: pS.cl }, { name: 'Section maximum', x: eta, y: eta.map(() => secMaxN), style: 'dash' }] },
        { type: 'line', title: 'Section profile drag along the span', xlabel: 'Span station 2y/b [-]', ylabel: 'cd [-]', series: [{ name: `Flight point`, x: eta, y: pF.cd }, { name: 'Maximum lift', x: eta, y: pS.cd }] },
        { type: 'line', title: 'Span loading at the flight point', xlabel: 'Span station 2y/b [-]', ylabel: 'cl·c [m]', series: [{ name: 'Viscous', x: eta, y: pF.strip.map((v) => 2 * v) }, { name: 'Elliptic, same lift', x: eta, y: eta.map((e) => ((4 * pF.CL * L.S) / (PI * i.b)) * Math.sqrt(Math.max(0, 1 - e * e))), style: 'dash' }] },
        { type: 'line', title: `Section polar from RANS, ${af.name}, Re = ${P0.Re.toExponential(2)}`, xlabel: 'Angle of attack [deg]', ylabel: 'cl and 50·cd [-]', series: [{ name: 'cl', x: P0.alpha, y: P0.cl, style: 'line+points' }, { name: '50 · cd', x: P0.alpha, y: P0.cd.map((v) => 50 * v), style: 'line+points' }], annotations: [{ x: P0.aStall, label: 'Section stall' }] },
        { type: 'bar', title: 'Drag breakdown at the flight point', ylabel: 'Drag counts (CD × 10⁴) [-]', categories: ['Wing profile (RANS)', 'Wing profile (build-up)', 'Induced', 'Other components', 'Total'], series: [{ name: 'CD', y: [cnt(CDpf), cnt(i.CD_wing_buildup), cnt(CDif), cnt(i.CD_other), cnt(CDf)] }] },
      ],
      tables: [
        { title: 'Wing polar', columns: ['α [deg]', 'CL', 'CD profile', 'CD induced', 'CD aircraft', 'L/D', 'CM'], rows: pts.filter((_, k) => k % 2 === 0 && k <= km + 2).map((p) => [N.deg(p.alpha).toFixed(0), +p.CL.toFixed(4), +p.CDp.toFixed(5), +p.CDi.toFixed(5), +(CDw(p) + i.CD_other).toFixed(5), +(p.CL / (CDw(p) + i.CD_other)).toFixed(2), +p.CM.toFixed(4)]) },
        { title: 'Section polar (body-fitted RANS, normal to the sweep line)', columns: ['α [deg]', 'cl', 'cd', 'cd pressure', 'cd friction', 'cm c/4', 'Separation x/c', 'y⁺', 'Converged'], rows: P0.alpha.map((v, k) => [v, +P0.cl[k].toFixed(4), +P0.cd[k].toFixed(5), +P0.cdp[k].toFixed(5), +P0.cdf[k].toFixed(5), +P0.cm[k].toFixed(4), +P0.xsep[k].toFixed(2), +P0.yplus[k].toFixed(1), P0.conv[k] ? 'yes' : 'no']) },
      ],
      outputs: gate ? { CL: CLf, CD: CDf, CLmax, LD_max: LDmax } : {},
      warnings,
      models: ['Vortex-lattice method with Prandtl–Glauert transformation and Trefftz-plane induced drag', 'Sectional decambering (α-shift) coupling: each strip is rotated until its lattice lift equals the viscous section lift at its effective angle of attack (van Dam; Mukherjee & Gopalarathnam)', `Section data: body-fitted incompressible RANS, Spalart–Allmaras, fully turbulent, ${G.wfn ? 'wall functions' : 'wall-resolved'}, ${tab.length} angles of attack at ${nRe} Reynolds number${nRe > 1 ? 's' : ''}`, 'Simple sweep theory for the section flow; Prandtl–Glauert scaling of section lift and moment', 'Aircraft polar: viscous wing + non-wing drag from the component build-up'],
      assumptions: ['Each strip behaves as a two-dimensional section at its effective angle of attack (high aspect ratio, no strong spanwise flow)', 'Same section along the span; rigid, planar wing; clean configuration', 'Profile drag integrated strip by strip; induced drag from the viscous span loading', nRe === 1 ? 'Friction part of the section drag scaled to the strip Reynolds number with the flat-plate law; pressure part and maximum lift taken at the mean-chord Reynolds number' : 'Section data interpolated in the logarithm of the strip Reynolds number'],
    };
  },
  convergence: { param: 'nSpan', label: 'Spanwise strips per half wing', levels: [6, 12, 24], metric: 'CLa_wingrans_per_rad' },
  calibration: { params: [{ key: 'CD_other', min: 0, max: 0.05 }], sweep: 'alpha_deg', target: 'CL_wingrans', note: 'Supply a measured wing or aircraft lift curve; the drag of the non-wing components is the uncertain handbook term in the aircraft polar.' },
  verify() {
    // coupling alone, with analytic section data: a section slope of 0.9 × 2π must reproduce lifting-line theory, and a lift cap must cap the wing
    const g = wingGeom({ S: 12, b: 12, taper: 1, sweep_deg: 0, twist_deg: 0 }), af = parseNaca('0012'), k = 0.9;
    const W = viscousWing(g, 16, 2, 1, af, 2 * PI, (j, ae) => ({ cl: k * 2 * PI * ae, cd: 0.01, cm: 0 }), [0, N.rad(2)]), slope = (W.pts[1].CL - W.pts[0].CL) / N.rad(2), ll = liftingLine(12, () => 1 / 12, () => 1, k * 2 * PI, 60), ll1 = liftingLine(12, () => 1 / 12, () => 1, 2 * PI, 60);
    const Wi = viscousWing(g, 16, 2, 1, af, 2 * PI, (j, ae) => ({ cl: 2 * PI * ae, cd: 0.01, cm: 0 }), [N.rad(3)]);
    const Wc = viscousWing(g, 12, 2, 1, af, 2 * PI, (j, ae) => ({ cl: Math.min(2 * PI * ae, 1.2), cd: 0.01, cm: 0 }), N.range(10, (q) => N.rad(2 * q + 2)));
    return [
      N.check('Coupled wing, section slope 0.9 × 2π: loss of wing lift slope against lifting line', slope / W.L.a.CL, ll.CL / ll1.CL, 0.01, 'Prandtl lifting line with section slopes 0.9 × 2π and 2π (rectangular wing, aspect ratio 12): the ratio isolates the coupling from the 2–3% difference between lattice and lifting line'),
      N.check('Inviscid section data return the vortex-lattice lift', Wi.pts[0].CL, Wi.L.a.CL * N.rad(3), 1e-3, 'Zero decambering when the section follows thin-airfoil theory'),
      N.check('Uniform section drag integrates to the same wing profile drag', Wi.pts[0].CDp, 0.01, 1e-9, 'Area-weighted strip integration'),
      N.check('A section lift limit of 1.2 limits the wing below it', N.amax(Wc.pts.map((p) => p.CL)) < 1.2 && N.amax(Wc.pts.map((p) => p.CL)) > 1.0 ? 1 : 0, 1, 1e-12, 'Wing maximum lift is below the section maximum because the span loading is not uniform'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!o.gate_passed) out.push({ severity: 'warn', title: 'Viscous wing result not published', detail: 'The run did not pass its sanity gate, so the rest of the app keeps the handbook drag build-up values.', action: 'Raise "Iterations per section angle", extend the table of section angles through the stall, or switch the wall treatment as the warning explains; then re-run.', basis: 'Internal gate: converged section polar, y⁺ in range, coupling converged, subcritical normal Mach number' });
    if (o.eta_stall_wingrans > 0.7) out.push({ severity: 'warn', title: 'Stall starts outboard', detail: `The first strip to stall is at ${(100 * o.eta_stall_wingrans).toFixed(0)}% semispan, in the aileron region.`, action: 'Add washout, reduce taper or use a higher-lift tip section; re-run until the stall starts inboard of about 60% semispan.', basis: 'Viscous strip analysis' });
    if (o.CL_wingrans > 0.85 * o.CLmax_wingrans) out.push({ severity: 'warn', title: 'Little stall margin at this point', detail: `CL = ${o.CL_wingrans.toFixed(2)} against a clean maximum of ${o.CLmax_wingrans.toFixed(2)}.`, action: 'Fly faster, lower the mass or deploy high-lift devices (outside this clean-wing model).', basis: 'CL / CLmax' });
    if (i.CD_wing_buildup > 0 && Math.abs(o.CD_wing_profile / i.CD_wing_buildup - 1) > 0.3) out.push({ severity: 'advise', title: 'Wing profile drag differs from the handbook build-up', detail: `RANS strips: ${(1e4 * o.CD_wing_profile).toFixed(0)} counts; flat-plate build-up: ${(1e4 * i.CD_wing_buildup).toFixed(0)} counts.`, action: 'Part of the difference is real (lift-dependent profile drag, fully turbulent sections, no laminar run) and part is grid error on the fast-look grid: refine the section grid before changing the build-up factors.', basis: 'Cross-check of two independent methods' });
    out.push({ severity: 'info', title: 'What this method adds', detail: `Maximum L/D ${o.LD_max_wingrans.toFixed(1)}, clean CLmax ${o.CLmax_wingrans.toFixed(2)} with viscous sections at flight Reynolds number.`, action: 'Use it to trade twist, taper and section choice for stall behaviour and cruise drag; one drag count is roughly ' + (100 / (1e4 * o.CD_wingrans)).toFixed(2) + '% of fuel burn and CO₂. Confirm final numbers with 3-D body-fitted RANS from the High-fidelity bridge.', basis: 'Fidelity ladder' });
    return out;
  },
};

// ---- 1-D compressible Euler: exact Riemann solver and finite-volume scheme ---------------------
/** Exact solution of the Riemann problem for an ideal gas (Toro). States {r, u, p}. Returns {ps, us, sample(ξ)}. */
export function riemannExact(L, R, g = GAMMA) {
  const aL = Math.sqrt((g * L.p) / L.r), aR = Math.sqrt((g * R.p) / R.r), g1 = (g - 1) / (2 * g), g2 = (g + 1) / (2 * g), g5 = 2 / (g + 1), g6 = (g - 1) / (g + 1);
  if ((2 * (aL + aR)) / (g - 1) <= R.u - L.u) return null; // vacuum is generated
  const fK = (p, K, a) => (p > K.p ? (p - K.p) * Math.sqrt(g5 / K.r / (p + g6 * K.p)) : ((2 * a) / (g - 1)) * ((p / K.p) ** g1 - 1));
  const F = (p) => fK(p, L, aL) + fK(p, R, aR) + R.u - L.u;
  let hi = Math.max(L.p, R.p); while (F(hi) < 0) hi *= 2;
  const ps = N.brent(F, 1e-12 * Math.min(L.p, R.p), hi, 1e-13 * hi, 300), us = 0.5 * (L.u + R.u) + 0.5 * (fK(ps, R, aR) - fK(ps, L, aL));
  const side = (K, a, sgn, xi) => { // sgn = -1 left, +1 right
    const pr = ps / K.p;
    if (ps > K.p) { const S = K.u + sgn * a * Math.sqrt(g2 * pr + g1); return sgn * (xi - S) > 0 ? K : { r: (K.r * (pr + g6)) / (g6 * pr + 1), u: us, p: ps }; }
    const as = a * pr ** g1, head = K.u + sgn * a, tail = us + sgn * as;
    if (sgn * (xi - head) > 0) return K;
    if (sgn * (xi - tail) < 0) return { r: K.r * pr ** (1 / g), u: us, p: ps };
    const c = g5 * (a - sgn * 0.5 * (g - 1) * (K.u - xi));
    return { r: K.r * (c / a) ** (2 / (g - 1)), u: g5 * (-sgn * a + 0.5 * (g - 1) * K.u + xi), p: K.p * (c / a) ** ((2 * g) / (g - 1)) };
  };
  const shock = (K, a, sgn) => (ps > K.p ? K.u + sgn * a * Math.sqrt(g2 * (ps / K.p) + g1) : NaN);
  return { ps, us, aL, aR, SL: shock(L, aL, -1), SR: shock(R, aR, 1), rsL: side(L, aL, -1, us - 1e-9 * (aL + aR)).r, rsR: side(R, aR, 1, us + 1e-9 * (aL + aR)).r, sample: (xi) => (xi <= us ? side(L, aL, -1, xi) : side(R, aR, 1, xi)) };
}
const _f = [0, 0, 0];
/** Numerical flux at an interface from left/right primitive states. */
function eulerFlux(kind, g, rl, ul, pl, rr, ur, pr) {
  const El = pl / (g - 1) + 0.5 * rl * ul * ul, Er = pr / (g - 1) + 0.5 * rr * ur * ur, al = Math.sqrt((g * pl) / rl), ar = Math.sqrt((g * pr) / rr);
  const Hl = (El + pl) / rl, Hr = (Er + pr) / rr, w = Math.sqrt(rr / rl), ut = (ul + w * ur) / (1 + w), Ht = (Hl + w * Hr) / (1 + w), at = Math.sqrt(Math.max((g - 1) * (Ht - 0.5 * ut * ut), 1e-12));
  const fl0 = rl * ul, fl1 = rl * ul * ul + pl, fl2 = ul * (El + pl), fr0 = rr * ur, fr1 = rr * ur * ur + pr, fr2 = ur * (Er + pr);
  if (kind === 'Rusanov') { const s = Math.max(Math.abs(ul) + al, Math.abs(ur) + ar); _f[0] = 0.5 * (fl0 + fr0) - 0.5 * s * (rr - rl); _f[1] = 0.5 * (fl1 + fr1) - 0.5 * s * (rr * ur - rl * ul); _f[2] = 0.5 * (fl2 + fr2) - 0.5 * s * (Er - El); return; }
  if (kind === 'Roe') {
    const rt = Math.sqrt(rl * rr), dr = rr - rl, du = ur - ul, dp = pr - pl, fix = (l) => { const d = 0.1 * at; return Math.abs(l) < d ? 0.5 * (l * l / d + d) : Math.abs(l); }; // Harten entropy fix
    const a1 = ((dp - rt * at * du) / (2 * at * at)) * fix(ut - at), a2 = (dr - dp / (at * at)) * fix(ut), a3 = ((dp + rt * at * du) / (2 * at * at)) * fix(ut + at);
    _f[0] = 0.5 * (fl0 + fr0) - 0.5 * (a1 + a2 + a3); _f[1] = 0.5 * (fl1 + fr1) - 0.5 * (a1 * (ut - at) + a2 * ut + a3 * (ut + at)); _f[2] = 0.5 * (fl2 + fr2) - 0.5 * (a1 * (Ht - ut * at) + 0.5 * a2 * ut * ut + a3 * (Ht + ut * at)); return;
  }
  const SL = Math.min(ul - al, ut - at), SR = Math.max(ur + ar, ut + at); // HLLC (Toro) with Einfeldt wave-speed estimates
  if (SL >= 0) { _f[0] = fl0; _f[1] = fl1; _f[2] = fl2; return; }
  if (SR <= 0) { _f[0] = fr0; _f[1] = fr1; _f[2] = fr2; return; }
  const Ss = (pr - pl + rl * ul * (SL - ul) - rr * ur * (SR - ur)) / (rl * (SL - ul) - rr * (SR - ur));
  if (Ss >= 0) { const c = (rl * (SL - ul)) / (SL - Ss); _f[0] = fl0 + SL * (c - rl); _f[1] = fl1 + SL * (c * Ss - rl * ul); _f[2] = fl2 + SL * (c * (El / rl + (Ss - ul) * (Ss + pl / (rl * (SL - ul)))) - El); }
  else { const c = (rr * (SR - ur)) / (SR - Ss); _f[0] = fr0 + SR * (c - rr); _f[1] = fr1 + SR * (c * Ss - rr * ur); _f[2] = fr2 + SR * (c * (Er / rr + (Ss - ur) * (Ss + pr / (rr * (SR - ur)))) - Er); }
}
const LIM = { 'MUSCL minmod': (a, b) => (a * b <= 0 ? 0 : Math.abs(a) < Math.abs(b) ? a : b), 'MUSCL van Leer': (a, b) => (a * b <= 0 ? 0 : (2 * a * b) / (a + b)), 'First order': () => 0 };
/** Finite-volume MUSCL–Hancock-type scheme (SSP-RK2 in time) for the 1-D Euler equations on [0, len] with transmissive ends. */
export function eulerFV(o) {
  const { n, g, len } = o, dx = len / n, nt = n + 4, lim = LIM[o.recon] || LIM['MUSCL minmod'];
  let U = [new Float64Array(nt), new Float64Array(nt), new Float64Array(nt)];
  for (let k = 0; k < nt; k++) { const x = (k - 1.5) * dx, K = x < o.x0 ? o.L : o.R; U[0][k] = K.r; U[1][k] = K.r * K.u; U[2][k] = K.p / (g - 1) + 0.5 * K.r * K.u * K.u; }
  const W = [new Float64Array(nt), new Float64Array(nt), new Float64Array(nt)], sl = [new Float64Array(nt), new Float64Array(nt), new Float64Array(nt)], F = [new Float64Array(nt), new Float64Array(nt), new Float64Array(nt)];
  const prim = (Q) => { let smax = 0; for (let k = 0; k < nt; k++) { const r = Math.max(Q[0][k], 1e-12), u = Q[1][k] / r, p = Math.max((g - 1) * (Q[2][k] - 0.5 * r * u * u), 1e-12); W[0][k] = r; W[1][k] = u; W[2][k] = p; smax = Math.max(smax, Math.abs(u) + Math.sqrt((g * p) / r)); } return smax; };
  const rhs = (Q, dt, out) => {
    for (let v = 0; v < 3; v++) { Q[v][0] = Q[v][1] = Q[v][2]; Q[v][nt - 1] = Q[v][nt - 2] = Q[v][nt - 3]; }
    prim(Q);
    for (let v = 0; v < 3; v++) for (let k = 1; k < nt - 1; k++) sl[v][k] = lim(W[v][k] - W[v][k - 1], W[v][k + 1] - W[v][k]);
    for (let k = 1; k < nt - 2; k++) { // interface k+1/2
      eulerFlux(o.flux, g, Math.max(W[0][k] + 0.5 * sl[0][k], 1e-12), W[1][k] + 0.5 * sl[1][k], Math.max(W[2][k] + 0.5 * sl[2][k], 1e-12), Math.max(W[0][k + 1] - 0.5 * sl[0][k + 1], 1e-12), W[1][k + 1] - 0.5 * sl[1][k + 1], Math.max(W[2][k + 1] - 0.5 * sl[2][k + 1], 1e-12));
      F[0][k] = _f[0]; F[1][k] = _f[1]; F[2][k] = _f[2];
    }
    for (let v = 0; v < 3; v++) for (let k = 2; k < nt - 2; k++) out[v][k] = Q[v][k] - (dt / dx) * (F[v][k] - F[v][k - 1]);
  };
  const U1 = [new Float64Array(nt), new Float64Array(nt), new Float64Array(nt)], U2 = [new Float64Array(nt), new Float64Array(nt), new Float64Array(nt)];
  let t = 0, steps = 0;
  while (t < o.t * (1 - 1e-12) && steps < 200000) {
    const dt = Math.min((o.cfl * dx) / prim(U), o.t - t);
    rhs(U, dt, U1); rhs(U1, dt, U2);
    for (let v = 0; v < 3; v++) for (let k = 2; k < nt - 2; k++) U[v][k] = 0.5 * (U[v][k] + U2[v][k]);
    t += dt; steps++;
  }
  prim(U);
  const x = N.range(n, (k) => (k + 0.5) * dx), pick = (A) => Array.from(A.subarray(2, nt - 2));
  return { x, r: pick(W[0]), u: pick(W[1]), p: pick(W[2]), E: pick(U[2]), steps, dx };
}
const shockTube = {
  id: 'shocktube', title: '1-D Euler finite-volume solver (shock tube)', fidelity: 'numerical',
  summary: 'Unsteady compressible flow in a tube after a diaphragm bursts, solved with a shock-capturing finite-volume scheme and compared cell by cell with the exact Riemann solution.',
  equations: ['Compressible Euler equations', 'Continuity equation', 'Conservation of momentum equation', 'Total energy equation', 'Ideal gas equation of state', 'Rankine–Hugoniot jump conditions', 'Reynolds transport theorem', 'Compressible finite-volume solvers', 'High-resolution shock-capturing schemes'],
  inputs: [
    { key: 'pL', label: 'Left (driver) pressure', unit: 'Pa', default: 101325, min: 1, group: 'Initial state' },
    { key: 'rhoL', label: 'Left density', unit: 'kg/m³', default: 1.225, min: 1e-4, group: 'Initial state' },
    { key: 'uL', label: 'Left velocity', unit: 'm/s', default: 0, group: 'Initial state' },
    { key: 'pR', label: 'Right (driven) pressure', unit: 'Pa', default: 10132.5, min: 1, group: 'Initial state' },
    { key: 'rhoR', label: 'Right density', unit: 'kg/m³', default: 0.153125, min: 1e-4, group: 'Initial state' },
    { key: 'uR', label: 'Right velocity', unit: 'm/s', default: 0, group: 'Initial state' },
    { key: 'gamma', label: 'Ratio of specific heats', unit: '-', default: 1.4, min: 1.05, max: 1.67, group: 'Gas' },
    { key: 'len', label: 'Tube length', unit: 'm', default: 1, min: 0.01, group: 'Domain' },
    { key: 'x0_frac', label: 'Diaphragm position', unit: 'x/L', default: 0.5, min: 0.1, max: 0.9, group: 'Domain' },
    { key: 't_end', label: 'End time', unit: 's', default: 6.954e-4, min: 1e-7, group: 'Domain', help: 'Keep the waves inside the tube; the default is the Sod problem at t = 0.2 in reference units' },
    { key: 'flux', label: 'Flux function', type: 'select', options: ['HLLC', 'Roe', 'Rusanov'], default: 'HLLC', group: 'Numerics' },
    { key: 'recon', label: 'Reconstruction', type: 'select', options: ['MUSCL minmod', 'MUSCL van Leer', 'First order'], default: 'MUSCL minmod', group: 'Numerics' },
    { key: 'cfl', label: 'CFL number', unit: '-', default: 0.5, min: 0.05, max: 0.9, group: 'Numerics' },
    { key: 'nCells', label: 'Cells', unit: '', default: 200, min: 20, max: 4000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: () => ({}),
  run(i) {
    const g = i.gamma, L = { r: i.rhoL, u: i.uL, p: i.pL }, R = { r: i.rhoR, u: i.uR, p: i.pR }, x0 = i.x0_frac * i.len, n = Math.round(N.clamp(i.nCells, 20, 4000));
    const ex = riemannExact(L, R, g), s = eulerFV({ n, g, len: i.len, x0, L, R, t: i.t_end, cfl: N.clamp(i.cfl, 0.05, 0.9), flux: i.flux, recon: i.recon }), warnings = [];
    const Rg = 287.05, xe = s.x, e = ex ? xe.map((x) => ex.sample((x - x0) / i.t_end)) : null;
    const l1 = e ? N.sum(s.r.map((v, k) => Math.abs(v - e[k].r))) * s.dx / (N.sum(e.map((v) => Math.abs(v.r))) * s.dx) : NaN;
    const mass0 = L.r * x0 + R.r * (i.len - x0), mass = N.sum(s.r) * s.dx, en0 = (L.p / (g - 1) + 0.5 * L.r * L.u ** 2) * x0 + (R.p / (g - 1) + 0.5 * R.r * R.u ** 2) * (i.len - x0);
    const inFlow = (K) => K.r * K.u * i.t_end, massErr = (mass - (mass0 + inFlow(L) - inFlow(R))) / mass0;
    const ke = N.sum(s.r.map((v, k) => 0.5 * v * s.u[k] ** 2)) * s.dx;
    // numerical plateau between contact and right wave, and captured shock position (half-jump in pressure)
    let pStarN = NaN, uStarN = NaN, xShN = NaN, xShE = NaN, Ms = 1;
    if (ex) {
      const xc = x0 + ex.us * i.t_end, xr = x0 + (Number.isFinite(ex.SR) ? ex.SR : ex.us + ex.aR) * i.t_end, a = xc + 0.3 * (xr - xc), b = xc + 0.7 * (xr - xc), idx = xe.map((x, k) => k).filter((k) => xe[k] > a && xe[k] < b);
      if (idx.length) { pStarN = N.mean(idx.map((k) => s.p[k])); uStarN = N.mean(idx.map((k) => s.u[k])); }
      if (Number.isFinite(ex.SR)) { xShE = x0 + ex.SR * i.t_end; const pm = 0.5 * (ex.ps + R.p); for (let k = n - 2; k > 0; k--) if (s.p[k] >= pm && s.p[k + 1] < pm) { xShN = xe[k] + ((s.p[k] - pm) / (s.p[k] - s.p[k + 1])) * s.dx; break; } Ms = (ex.SR - R.u) / ex.aR; }
      else if (Number.isFinite(ex.SL)) Ms = (L.u - ex.SL) / ex.aL;
      if ((Number.isFinite(ex.SR) ? xShE : x0 + (R.u + ex.aR) * i.t_end) > i.len || x0 + Math.min(L.u - ex.aL, Number.isFinite(ex.SL) ? ex.SL : 0) * i.t_end < 0) warnings.push('A wave has reached the end of the tube: the transmissive boundary lets it leave, so the comparison with the infinite-domain exact solution is only valid inside the tube.');
      if (!Number.isFinite(ex.SR)) warnings.push('The right-running wave is an expansion, not a shock: shock position is undefined for these states.');
    } else warnings.push('The initial states separate fast enough to create a vacuum; the exact reference solution is not evaluated.');
    const T = (r, p) => p / (Rg * r), d = (y) => thin(xe, y, 400);
    const exS = e ? [{ name: 'Exact Riemann solution', style: 'dash' }] : [];
    const pl = (title, ylabel, num, exv) => ({ type: 'line', title, xlabel: 'x [m]', ylabel, series: [{ name: `${i.flux}, ${i.recon}`, x: d(num)[0], y: d(num)[1] }, ...exS.map((q) => ({ ...q, x: d(exv)[0], y: d(exv)[1] }))] });
    return {
      kpis: [
        { key: 'p_star_Pa', label: 'Pressure behind the shock (computed)', value: pStarN, unit: 'Pa' },
        { key: 'p_star_exact_Pa', label: 'Pressure behind the shock (exact)', value: ex ? ex.ps : NaN, unit: 'Pa' },
        { key: 'u_star_ms', label: 'Contact velocity (computed)', value: uStarN, unit: 'm/s' },
        { key: 'u_star_exact_ms', label: 'Contact velocity (exact)', value: ex ? ex.us : NaN, unit: 'm/s' },
        { key: 'shock_mach', label: 'Shock Mach number', value: Ms, unit: '-', note: Ms > 1 ? 'relative to the gas ahead' : 'no shock forms' },
        { key: 'x_shock_m', label: 'Shock position (captured)', value: xShN, unit: 'm' },
        { key: 'x_shock_exact_m', label: 'Shock position (exact)', value: xShE, unit: 'm' },
        { key: 'T_post_shock_K', label: 'Temperature behind the shock (air)', value: ex ? T(ex.rsR, ex.ps) : NaN, unit: 'K' },
        { key: 'l1_error_rho', label: 'Relative L1 density error', value: l1, unit: '-', status: l1 < 0.02 ? 'ok' : 'warn' },
        { key: 'mass_error', label: 'Mass conservation error', value: massErr, unit: '-', status: Math.abs(massErr) < 1e-9 ? 'ok' : 'warn' },
        { key: 'ke_J_m2', label: 'Kinetic energy in the tube', value: ke, unit: 'J/m²' },
        { key: 'time_steps', label: 'Time steps', value: s.steps, unit: '' },
      ],
      plots: [pl('Density', 'Density [kg/m³]', s.r, e?.map((v) => v.r)), pl('Velocity', 'Velocity [m/s]', s.u, e?.map((v) => v.u)), pl('Pressure', 'Pressure [kPa]', s.p.map((v) => v / 1e3), e?.map((v) => v.p / 1e3)), pl('Temperature (air)', 'Temperature [K]', s.r.map((v, k) => T(v, s.p[k])), e?.map((v) => T(v.r, v.p)))],
      outputs: { energy_error: (N.sum(s.E) * s.dx - (en0 + ((L.p / (g - 1) + 0.5 * L.r * L.u ** 2 + L.p) * L.u - (R.p / (g - 1) + 0.5 * R.r * R.u ** 2 + R.p) * R.u) * i.t_end)) / en0 },
      warnings,
      models: [`Finite-volume Godunov-type scheme with ${i.flux} flux`, `${i.recon} reconstruction of primitive variables`, 'Two-stage strong-stability-preserving Runge–Kutta time integration', 'Exact Riemann solver (reference)'],
      assumptions: ['Inviscid, adiabatic, calorically perfect gas', 'One-dimensional flow; transmissive tube ends', 'Temperatures use the gas constant of air'],
    };
  },
  convergence: { param: 'nCells', label: 'Finite-volume cells', levels: [50, 100, 200, 400, 800], metric: 'ke_J_m2' },
  verify() {
    const sod = { pL: 1, rhoL: 1, uL: 0, pR: 0.1, rhoR: 0.125, uR: 0, gamma: 1.4, len: 1, x0_frac: 0.5, t_end: 0.2, flux: 'HLLC', recon: 'MUSCL minmod', cfl: 0.5, nCells: 400 };
    const ex = riemannExact({ r: 1, u: 0, p: 1 }, { r: 0.125, u: 0, p: 0.1 }), o = N.kv(shockTube.run(sod)), roe = N.kv(shockTube.run({ ...sod, flux: 'Roe' })), t123 = riemannExact({ r: 1, u: -2, p: 0.4 }, { r: 1, u: 2, p: 0.4 });
    return [
      N.check('Exact Riemann solver: Sod p*', ex.ps, 0.30313, 2e-5, 'Toro, Riemann Solvers and Numerical Methods, Test 1'),
      N.check('Exact Riemann solver: Sod u*', ex.us, 0.92745, 2e-5, 'Toro, Test 1'),
      N.check('Exact Riemann solver: Sod ρ*L', ex.rsL, 0.42632, 2e-5, 'Toro, Test 1'),
      N.check('Exact Riemann solver: Sod ρ*R', ex.rsR, 0.26557, 2e-5, 'Toro, Test 1'),
      N.check('Exact Riemann solver: two-rarefaction p*', t123.ps, 0.00189, 5e-3, 'Toro, Test 2'),
      N.check('HLLC plateau pressure equals exact p*', o.p_star_Pa, ex.ps, 2e-3, 'Exact Riemann solution'),
      N.check('Captured shock position (HLLC)', o.x_shock_m, o.x_shock_exact_m, 4e-3, 'Rankine–Hugoniot shock speed; one cell is 0.25% of the tube'),
      N.check('Captured shock position (Roe)', roe.x_shock_m, o.x_shock_exact_m, 4e-3, 'Rankine–Hugoniot shock speed'),
      N.check('Mass is conserved to round-off', o.mass_error, 0, 1e-12, 'Conservative finite-volume update'),
      N.check('L1 density error at 400 cells', o.l1_error_rho < 0.006 ? 1 : 0, 1, 1e-12, 'First-order convergence at discontinuities'),
    ];
  },
  validation: [{ name: 'Sod shock tube, star-region pressure', source: 'Exact Riemann solution (Sod 1978; Toro Test 1): p* = 0.30313', inputs: { pL: 1, rhoL: 1, uL: 0, pR: 0.1, rhoR: 0.125, uR: 0, gamma: 1.4, len: 1, x0_frac: 0.5, t_end: 0.2 }, sweep: { key: 'nCells', values: [200, 400] }, target: 'p_star_Pa', observed: [0.30313, 0.30313], tol_pct: 0.5 }],
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.l1_error_rho > 0.02) out.push({ severity: 'advise', title: 'Discontinuities are smeared', detail: `Relative L1 density error is ${(100 * o.l1_error_rho).toFixed(2)}%.`, action: 'Double the number of cells or switch from first-order to MUSCL reconstruction; run the convergence study to quantify the discretisation error.', basis: 'Comparison with the exact Riemann solution' });
    if (i.recon === 'First order') out.push({ severity: 'info', title: 'First-order scheme selected', detail: 'Monotone but very diffusive: contacts spread over many cells.', action: 'Use MUSCL with a limiter for production accuracy; keep first order only as a robustness fallback.', basis: 'Godunov’s theorem' });
    if (o.T_post_shock_K > 1500) out.push({ severity: 'warn', title: 'High post-shock temperature', detail: `${o.T_post_shock_K.toFixed(0)} K: vibrational excitation and dissociation make γ variable.`, action: 'Treat results as indicative and use an equilibrium or finite-rate real-gas solver.', basis: 'Calorically perfect gas limit (about 800–1000 K for air)' });
    return out;
  },
};

// ---- shock and expansion relations ---------------------------------------------------------------
const nShock = (M, g) => { const m2 = M * M, p = 1 + ((2 * g) / (g + 1)) * (m2 - 1), r = ((g + 1) * m2) / ((g - 1) * m2 + 2); return { p, r, T: p / r, M2: Math.sqrt((1 + 0.5 * (g - 1) * m2) / (g * m2 - 0.5 * (g - 1))), p0: r ** (g / (g - 1)) * (1 / p) ** (1 / (g - 1)) }; };
const thetaOf = (M, b, g) => Math.atan((2 * (M * M * Math.sin(b) ** 2 - 1)) / (Math.tan(b) * (M * M * (g + Math.cos(2 * b)) + 2)));
const pm = (M, g) => { const k = Math.sqrt((g + 1) / (g - 1)); return k * Math.atan(Math.sqrt(M * M - 1) / k) - Math.atan(Math.sqrt(M * M - 1)); };
const pmInv = (nu, g) => N.brent((M) => pm(M, g) - nu, 1, 60, 1e-12);
/** Oblique shock for deflection th [rad]: weak and strong wave angles, or detached. */
export function oblique(M, th, g = GAMMA) {
  const mu = Math.asin(1 / M), bMax = N.goldenSection((b) => -thetaOf(M, b, g), mu, PI / 2, 1e-10), thMax = thetaOf(M, bMax, g);
  if (th > thMax) return { detached: true, thMax, bMax };
  const bw = N.brent((b) => thetaOf(M, b, g) - th, mu + 1e-12, bMax, 1e-13), bs = N.brent((b) => thetaOf(M, b, g) - th, bMax, PI / 2, 1e-13), sh = nShock(M * Math.sin(bw), g);
  return { detached: false, thMax, bMax, beta: bw, betaStrong: bs, M2: sh.M2 / Math.sin(bw - th), p: sh.p, r: sh.r, T: sh.T, p0: sh.p0 };
}
/** Turn a supersonic stream through δ [rad] (compression positive): pressure ratio and downstream Mach number. */
function turn(M, d, g) {
  if (Math.abs(d) < 1e-12) return { M2: M, p: 1, kind: 'none' };
  if (d > 0) { const o = oblique(M, d, g); if (o.detached) { const n = nShock(M, g); return { M2: n.M2, p: n.p, kind: 'detached' }; } return { M2: o.M2, p: o.p, kind: 'shock' }; }
  const nuMax = (PI / 2) * (Math.sqrt((g + 1) / (g - 1)) - 1), nu2 = Math.min(pm(M, g) - d, nuMax * 0.999), M2 = pmInv(nu2, g), t = (m) => 1 + 0.5 * (g - 1) * m * m;
  return { M2, p: (t(M) / t(M2)) ** (g / (g - 1)), kind: 'expansion' };
}
/** Shock–expansion lift and drag of a symmetric double-wedge (diamond) section. */
function diamond(M, al, tc, g) {
  const eps = Math.atan(tc), uf = turn(M, eps - al, g), ur = M > 1 && uf.M2 > 1 ? turn(uf.M2, -2 * eps, g) : { p: 1, M2: uf.M2 }, lf = turn(M, eps + al, g), lr = lf.M2 > 1 ? turn(lf.M2, -2 * eps, g) : { p: 1, M2: lf.M2 };
  const p = [uf.p, uf.p * ur.p, lf.p, lf.p * lr.p], q = 0.5 * g * M * M, fy = 0.5 * (p[2] + p[3] - p[0] - p[1]), fx = 0.5 * Math.tan(eps) * (p[0] + p[2] - p[1] - p[3]);
  return { cl: (fy * Math.cos(al) - fx * Math.sin(al)) / q, cd: (fy * Math.sin(al) + fx * Math.cos(al)) / q, p, detached: uf.kind === 'detached' || lf.kind === 'detached' || uf.M2 <= 1 || lf.M2 <= 1 };
}
const shocks = {
  id: 'shocks', title: 'Shock waves, expansions and supersonic sections', fidelity: 'analytical',
  summary: 'Normal and oblique shock jumps, the θ–β–M relation, Prandtl–Meyer expansion and the lift and wave drag of a thin supersonic section by shock–expansion theory, for intakes, wedges, control surfaces and local supersonic pockets.',
  equations: ['Rankine–Hugoniot jump conditions', 'Normal and oblique shock relations', 'Ideal gas equation of state', 'Compressible Euler equations', 'Entropy transport equation', 'Compressible potential-flow equation'],
  inputs: [
    { key: 'M1', label: 'Upstream Mach number', unit: '-', default: 2, min: 1.01, max: 10, group: 'Flow', help: 'For a transonic wing use the local Mach number ahead of the shock (1.1–1.4)' },
    { key: 'theta_deg', label: 'Wedge or ramp deflection', unit: 'deg', default: 10, min: 0, max: 45, group: 'Geometry' },
    { key: 'turn_deg', label: 'Expansion corner angle', unit: 'deg', default: 10, min: 0, max: 60, group: 'Geometry' },
    { key: 'alpha_deg', label: 'Section angle of attack', unit: 'deg', default: 4, min: 0, max: 20, group: 'Geometry' },
    { key: 'tc', label: 'Diamond-section thickness ratio', unit: '-', default: 0.05, min: 0, max: 0.2, group: 'Geometry', help: '0 gives the flat plate' },
    { key: 'gamma', label: 'Ratio of specific heats', unit: '-', default: 1.4, min: 1.05, max: 1.67, group: 'Flow' },
    ...ATM,
  ],
  defaults: (c) => ({ alt_m: c.atm.alt_m, dISA: c.atm.dISA_K }),
  run(i) {
    const g = i.gamma, M = Math.max(i.M1, 1.001), a = isa(i.alt_m, i.dISA), th = N.rad(i.theta_deg), ns = nShock(M, g), ob = oblique(M, th, g), ex = turn(M, -N.rad(i.turn_deg), g), al = N.rad(i.alpha_deg);
    const dm = diamond(M, al, i.tc, g), B = Math.sqrt(M * M - 1), clA = (4 * al) / B, cdA = (4 * (al * al + i.tc * i.tc)) / B, T0 = a.T * (1 + 0.5 * (g - 1) * M * M), warnings = [];
    if (ob.detached) warnings.push(`A ${i.theta_deg}° deflection exceeds the maximum ${N.deg(ob.thMax).toFixed(2)}° for an attached shock at Mach ${M.toFixed(2)}: the shock detaches into a bow wave, which needs a field solution.`);
    if (dm.detached) warnings.push('The leading-edge shock of the section is detached or the flow behind it is subsonic: shock–expansion theory does not apply; the plotted value uses a normal-shock pressure as a bound.');
    if (ns.T * a.T > 1000) warnings.push('Post-shock temperature exceeds about 1000 K: γ is no longer constant and real-gas effects reduce the temperature rise.');
    if (ns.p > 1.5 && M < 1.6) warnings.push(`Normal-shock pressure ratio ${ns.p.toFixed(2)} at M = ${M.toFixed(2)}: a turbulent boundary layer separates at the shock foot when the upstream Mach number exceeds about 1.3 (shock–boundary-layer interaction is not modelled).`);
    const bs = N.linspace(0.5, 89.9, 180).map(N.rad), tb = [1.5, 2, 3, 5, M].filter((m, k, arr) => arr.indexOf(m) === k).slice(0, 5).map((m) => { const pts = bs.filter((b) => b > Math.asin(1 / m)); return { name: `M = ${m.toFixed(2)}`, x: pts.map((b) => N.deg(thetaOf(m, b, g))), y: pts.map(N.deg) }; });
    const Ms = N.linspace(1.02, Math.max(4, M * 1.2), 80), As = N.linspace(0, 12, 25);
    return {
      kpis: [
        { key: 'p2_p1_normal', label: 'Normal shock pressure ratio', value: ns.p, unit: '-' },
        { key: 'M2_normal', label: 'Mach number behind normal shock', value: ns.M2, unit: '-' },
        { key: 'p02_p01_normal', label: 'Normal shock total-pressure recovery', value: ns.p0, unit: '-', status: ns.p0 > 0.9 ? 'ok' : 'warn', note: 'Intake recovery: below 0.9 calls for oblique-shock compression' },
        { key: 'beta_weak_deg', label: 'Oblique shock angle (weak)', value: ob.detached ? NaN : N.deg(ob.beta), unit: 'deg' },
        { key: 'p2_p1_oblique', label: 'Oblique shock pressure ratio', value: ob.detached ? NaN : ob.p, unit: '-' },
        { key: 'M2_oblique', label: 'Mach number behind oblique shock', value: ob.detached ? NaN : ob.M2, unit: '-' },
        { key: 'p02_p01_oblique', label: 'Oblique shock total-pressure recovery', value: ob.detached ? NaN : ob.p0, unit: '-' },
        { key: 'theta_max_deg', label: 'Maximum attached-shock deflection', value: N.deg(ob.thMax), unit: 'deg', status: ob.detached ? 'bad' : 'ok' },
        { key: 'nu_deg', label: 'Prandtl–Meyer angle of the stream', value: N.deg(pm(M, g)), unit: 'deg' },
        { key: 'M2_expansion', label: 'Mach number after expansion', value: ex.M2, unit: '-' },
        { key: 'p2_p1_expansion', label: 'Expansion pressure ratio', value: ex.p, unit: '-' },
        { key: 'cl_supersonic', label: 'Section lift (shock–expansion)', value: dm.cl, unit: '-' },
        { key: 'cd_wave', label: 'Section wave drag (shock–expansion)', value: dm.cd, unit: '-' },
        { key: 'cl_ackeret', label: 'Section lift (Ackeret linear theory)', value: clA, unit: '-' },
        { key: 'ld_supersonic', label: 'Inviscid section lift-to-drag', value: dm.cd > 0 ? dm.cl / dm.cd : NaN, unit: '-' },
        { key: 'T0_K', label: 'Stagnation temperature', value: T0, unit: 'K', note: 'Upper bound of skin recovery temperature' },
        { key: 'p_pitot_Pa', label: 'Pitot pressure behind the bow shock', value: a.p * ns.p0 * (1 + 0.5 * (g - 1) * M * M) ** (g / (g - 1)), unit: 'Pa' },
      ],
      plots: [
        { type: 'line', title: 'θ–β–M: shock angle against flow deflection', xlabel: 'Deflection θ [deg]', ylabel: 'Shock angle β [deg]', series: tb, annotations: [{ x: i.theta_deg, label: 'This wedge' }] },
        { type: 'line', title: 'Normal shock ratios', xlabel: 'Upstream Mach number [-]', ylabel: 'Ratio [-]', series: [{ name: 'Total-pressure recovery p02/p01', x: Ms, y: Ms.map((m) => nShock(m, g).p0) }, { name: 'Downstream Mach number', x: Ms, y: Ms.map((m) => nShock(m, g).M2) }, { name: 'Density ratio ÷ 6', x: Ms, y: Ms.map((m) => nShock(m, g).r / 6) }], annotations: [{ x: M, label: 'This case' }] },
        { type: 'line', title: 'Supersonic section lift', xlabel: 'Angle of attack [deg]', ylabel: 'cl [-]', series: [{ name: 'Shock–expansion', x: As, y: As.map((ad) => diamond(M, N.rad(ad), i.tc, g).cl) }, { name: 'Ackeret linear theory', x: As, y: As.map((ad) => (4 * N.rad(ad)) / B), style: 'dash' }] },
        { type: 'line', title: 'Supersonic section wave-drag polar', xlabel: 'cd [-]', ylabel: 'cl [-]', series: [{ name: 'Shock–expansion', x: As.map((ad) => diamond(M, N.rad(ad), i.tc, g).cd), y: As.map((ad) => diamond(M, N.rad(ad), i.tc, g).cl) }, { name: 'Ackeret linear theory', x: As.map((ad) => (4 * (N.rad(ad) ** 2 + i.tc ** 2)) / B), y: As.map((ad) => (4 * N.rad(ad)) / B), style: 'dash' }] },
      ],
      outputs: { cd_ackeret: cdA, T2_normal_K: a.T * ns.T, p2_normal_Pa: a.p * ns.p },
      warnings, models: ['Rankine–Hugoniot normal and oblique shock relations', 'θ–β–M relation (weak and strong branches)', 'Prandtl–Meyer expansion', 'Shock–expansion theory for a double-wedge section', 'Ackeret linearised supersonic theory'],
      assumptions: ['Steady inviscid flow of a calorically perfect gas', 'Attached, plane waves; no wave interaction or reflection downstream of the section', 'Friction drag is not included in the section drag'],
    };
  },
  verify() {
    const o = N.kv(shocks.run({ M1: 2, theta_deg: 10, turn_deg: 10, alpha_deg: 0.5, tc: 0, gamma: 1.4, alt_m: 0, dISA: 0 }));
    return [
      N.check('Normal shock M = 2: p2/p1', o.p2_p1_normal, 4.5, 1e-12, 'NACA Report 1135'),
      N.check('Normal shock M = 2: M2', o.M2_normal, 0.57735, 1e-5, 'NACA Report 1135'),
      N.check('Normal shock M = 2: p02/p01', o.p02_p01_normal, 0.72087, 1e-5, 'NACA Report 1135'),
      N.check('Oblique shock M = 2, θ = 10°: β', o.beta_weak_deg, 39.314, 1e-4, 'NACA Report 1135 θ–β–M'),
      N.check('Oblique shock M = 2, θ = 10°: p2/p1', o.p2_p1_oblique, 1.7066, 2e-4, 'NACA Report 1135'),
      N.check('Prandtl–Meyer angle at M = 2', o.nu_deg, 26.3798, 1e-5, 'NACA Report 1135'),
      N.check('Flat plate at small α approaches Ackeret', o.cl_supersonic, o.cl_ackeret, 0.01, 'Linearised supersonic theory, cl = 4α/√(M²−1)'),
    ];
  },
  validation: [{ name: 'Oblique shock angle at Mach 2', source: 'NACA Report 1135 (exact θ–β–M relation, γ = 1.4)', inputs: { M1: 2, gamma: 1.4 }, sweep: { key: 'theta_deg', values: [10, 20] }, target: 'beta_weak_deg', observed: [39.31, 53.42], tol_pct: 0.2 }],
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.p02_p01_normal < 0.9) out.push({ severity: 'advise', title: 'A single normal shock wastes total pressure', detail: `Recovery ${(100 * o.p02_p01_normal).toFixed(1)}% at Mach ${i.M1}; each 1% of intake recovery is roughly 1–1.5% of thrust and fuel.`, action: 'Use one or more oblique shocks (ramp or cone) ahead of the terminal normal shock.', basis: 'Rankine–Hugoniot total-pressure loss' });
    if (i.M1 > 1.3 && i.M1 < 1.6) out.push({ severity: 'warn', title: 'Shock strong enough to separate the boundary layer', detail: `Upstream Mach ${i.M1}: shock-induced separation, buffet and a drag rise are likely on a wing.`, action: 'Keep the local Mach number ahead of the wing shock below about 1.2–1.3 through section design; assess buffet onset with unsteady RANS or tunnel tests.', basis: 'Shock–boundary-layer interaction criterion (empirical, M ≈ 1.3)' });
    if (Number.isFinite(o.ld_supersonic) && i.tc > 0) out.push({ severity: 'info', title: 'Thickness wave drag', detail: `Thickness contributes cd ≈ ${((4 * i.tc * i.tc) / Math.sqrt(i.M1 ** 2 - 1)).toFixed(4)}, proportional to (t/c)².`, action: 'Halving thickness quarters the thickness wave drag; trade against structural depth in Suite 2.', basis: 'Ackeret theory' });
    return out;
  },
};

// ---- 2-D incompressible Navier–Stokes: projection method on a staggered grid ---------------------
/**
 * Direct solver for the cell-centred Poisson equation ∇²p = rhs with homogeneous Neumann walls on an n×n grid.
 * The discrete Laplacian is diagonalised exactly by the cosine modes cos(πk(i+½)/n), so the projection is exact to round-off.
 * Returns solve(P, rhs) working on flat row-major Float64Arrays.
 */
export function poissonNeumann(n, h) {
  const Q = new Float64Array(n * n), lam = new Float64Array(n), T = new Float64Array(n * n), W = new Float64Array(n * n);
  for (let k = 0; k < n; k++) { lam[k] = (2 * Math.cos((PI * k) / n) - 2) / (h * h); for (let i = 0; i < n; i++) Q[k * n + i] = Math.sqrt((k ? 2 : 1) / n) * Math.cos((PI * k * (i + 0.5)) / n); }
  const Qt = new Float64Array(n * n); for (let k = 0; k < n; k++) for (let i = 0; i < n; i++) Qt[i * n + k] = Q[k * n + i];
  const mul = (A, B, C) => { // C = A·B, loop order chosen for contiguous memory access
    C.fill(0);
    for (let i = 0; i < n; i++) { const io = i * n; for (let k = 0; k < n; k++) { const a = A[io + k], ko = k * n; for (let j = 0; j < n; j++) C[io + j] += a * B[ko + j]; } }
  };
  return (P, rhs) => {
    mul(Q, rhs, T); mul(T, Qt, W);
    for (let k = 0; k < n; k++) for (let l = 0; l < n; l++) W[k * n + l] = k || l ? W[k * n + l] / (lam[k] + lam[l]) : 0;
    mul(Qt, W, T); mul(T, Q, P);
  };
}
/** Lid-driven cavity (unit square, unit lid speed) by Chorin's projection method on a MAC grid. */
export function cavity(n, Re, cfl, tMax, tol, upw, progress) {
  const h = 1 / n, nu = 1 / Re, w = n + 2, dt = cfl * Math.min(h, 0.25 * h * h * Re, 2 / Re);
  const U = new Float64Array((n + 1) * w), V = new Float64Array(w * (n + 1)), Us = new Float64Array(U.length), Vs = new Float64Array(V.length), P = new Float64Array(n * n), rhs = new Float64Array(n * n);
  const iu = (i, j) => i * w + j + 1, iv = (i, j) => (i + 1) * (n + 1) + j; // u(i,j): i = 0..n, j = -1..n; v(i,j): i = -1..n, j = 0..n
  const hist = { t: [], r: [] }, solveP = poissonNeumann(n, h); let t = 0, steps = 0, res = 1;
  const maxSteps = Math.ceil(tMax / dt), every = Math.max(1, Math.floor(maxSteps / 250));
  while (steps < maxSteps) {
    for (let i = 0; i <= n; i++) { U[iu(i, -1)] = -U[iu(i, 0)]; U[iu(i, n)] = 2 - U[iu(i, n - 1)]; }
    for (let j = 0; j <= n; j++) { V[iv(-1, j)] = -V[iv(0, j)]; V[iv(n, j)] = -V[iv(n - 1, j)]; }
    Us.set(U); Vs.set(V);
    for (let i = 1; i < n; i++) for (let j = 0; j < n; j++) {
      const k = iu(i, j), u = U[k], ue = 0.5 * (u + U[k + w]), uw = 0.5 * (U[k - w] + u), un = 0.5 * (u + U[k + 1]), us = 0.5 * (U[k - 1] + u), vn = 0.5 * (V[iv(i - 1, j + 1)] + V[iv(i, j + 1)]), vs = 0.5 * (V[iv(i - 1, j)] + V[iv(i, j)]);
      const adv = (ue * ue - uw * uw + vn * un - vs * us + upw * 0.5 * (Math.abs(ue) * (u - U[k + w]) - Math.abs(uw) * (U[k - w] - u) + Math.abs(vn) * (u - U[k + 1]) - Math.abs(vs) * (U[k - 1] - u))) / h;
      Us[k] = u + dt * (-adv + (nu * (U[k + w] + U[k - w] + U[k + 1] + U[k - 1] - 4 * u)) / (h * h));
    }
    for (let i = 0; i < n; i++) for (let j = 1; j < n; j++) {
      const k = iv(i, j), v = V[k], m = n + 1, vn = 0.5 * (v + V[k + 1]), vs = 0.5 * (V[k - 1] + v), ve = 0.5 * (v + V[k + m]), vw = 0.5 * (V[k - m] + v), ue = 0.5 * (U[iu(i + 1, j - 1)] + U[iu(i + 1, j)]), uw = 0.5 * (U[iu(i, j - 1)] + U[iu(i, j)]);
      const adv = (vn * vn - vs * vs + ue * ve - uw * vw + upw * 0.5 * (Math.abs(vn) * (v - V[k + 1]) - Math.abs(vs) * (V[k - 1] - v) + Math.abs(ue) * (v - V[k + m]) - Math.abs(uw) * (V[k - m] - v))) / h;
      Vs[k] = v + dt * (-adv + (nu * (V[k + m] + V[k - m] + V[k + 1] + V[k - 1] - 4 * v)) / (h * h));
    }
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) rhs[j * n + i] = (Us[iu(i + 1, j)] - Us[iu(i, j)] + Vs[iv(i, j + 1)] - Vs[iv(i, j)]) / (h * dt);
    solveP(P, rhs);
    res = 0;
    for (let i = 1; i < n; i++) for (let j = 0; j < n; j++) { const k = iu(i, j), un = Us[k] - (dt * (P[j * n + i] - P[j * n + i - 1])) / h; res = Math.max(res, Math.abs(un - U[k])); U[k] = un; }
    for (let i = 0; i < n; i++) for (let j = 1; j < n; j++) { const k = iv(i, j), vn = Vs[k] - (dt * (P[j * n + i] - P[(j - 1) * n + i])) / h; res = Math.max(res, Math.abs(vn - V[k])); V[k] = vn; }
    res /= dt; t += dt; steps++;
    if (steps % every === 0 || res < tol) { hist.t.push(t); hist.r.push(Math.max(res, 1e-16)); progress?.(steps / maxSteps, `t = ${t.toFixed(2)}`); }
    if (res < tol || !Number.isFinite(res)) break;
  }
  let div = 0; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) div = Math.max(div, Math.abs(U[iu(i + 1, j)] - U[iu(i, j)] + V[iv(i, j + 1)] - V[iv(i, j)]) / h);
  const u = (i, j) => U[iu(i, j)], v = (i, j) => V[iv(i, j)];
  // stream function at grid nodes (ψ = 0 on the walls), cell-centre speed and pressure
  const psi = N.range(n + 1, () => new Array(n + 1).fill(0));
  for (let i = 0; i <= n; i++) for (let j = 1; j <= n; j++) psi[j][i] = psi[j - 1][i] + u(i, j - 1) * h;
  const speed = N.range(n, (j) => N.range(n, (i) => Math.hypot(0.5 * (u(i, j) + u(i + 1, j)), 0.5 * (v(i, j) + v(i, j + 1)))));
  return { n, h, dt, t, steps, res, div, hist, psi, speed, uc: N.range(n, (j) => u(n / 2, j)), vc: N.range(n, (i) => v(i, n / 2)), p: N.range(n, (j) => N.range(n, (i) => P[j * n + i])) };
}
/** Parabolic refinement of a discrete extremum: returns [x, y]. */
const vertex = (xs, ys, k) => { if (k <= 0 || k >= ys.length - 1) return [xs[k], ys[k]]; const d = ys[k - 1] - 2 * ys[k] + ys[k + 1], s = d ? (0.5 * (ys[k - 1] - ys[k + 1])) / d : 0; return [xs[k] + s * (xs[k + 1] - xs[k]), ys[k] - 0.25 * (ys[k - 1] - ys[k + 1]) * s]; };
const cavityFlow = {
  id: 'cavity', title: '2-D Navier–Stokes field solver (lid-driven cavity)', fidelity: 'numerical',
  summary: 'Viscous incompressible flow in a square cavity driven by a moving wall, time-marched to steady state with a pressure-projection finite-volume method. It is the classic benchmark for Navier–Stokes solvers and a model of the recirculation in landing-gear, weapon and flap-cove cavities.',
  equations: ['Navier–Stokes equations', 'Continuity equation', 'Conservation of momentum equation', 'Vorticity transport equation', 'Reynolds transport theorem'],
  inputs: [
    { key: 'Re', label: 'Reynolds number U·L/ν', unit: '-', default: 100, min: 1, max: 5000, group: 'Flow', help: 'Based on lid speed and cavity width. Laminar steady solutions exist up to about 7500' },
    { key: 'L_m', label: 'Cavity width', unit: 'm', default: 0.5, min: 1e-4, group: 'Flow', help: 'Only used to give dimensional time and speed scales' },
    { key: 'n', label: 'Cells per side', unit: '', default: 24, min: 8, max: 96, step: 2, discrete: true, group: 'Numerics' },
    { key: 'cfl', label: 'Time-step safety factor', unit: '-', default: 0.8, min: 0.05, max: 0.95, group: 'Numerics', help: 'Fraction of the advective, diffusive and cell-Reynolds stability limits' },
    { key: 't_max', label: 'Maximum simulated time', unit: 'L/U', default: 40, min: 1, max: 400, group: 'Numerics' },
    { key: 'tol', label: 'Steady-state residual', unit: 'U²/L', default: 1e-4, min: 1e-9, max: 1e-2, group: 'Numerics', help: 'Stop when the largest velocity change per unit time falls below this' },
    { key: 'upwind', label: 'Upwind blending', unit: '-', default: 0, min: 0, max: 1, group: 'Numerics', help: '0 = second-order central; raise towards 1 if the cell Reynolds number exceeds 2 and wiggles appear' },
  ],
  defaults: () => ({}),
  run(i, ctx) {
    const n = even(i.n, 8, 96), nu = isa(ctx?.case?.atm?.alt_m ?? 0).nu, r = cavity(n, i.Re, N.clamp(i.cfl, 0.05, 0.95), i.t_max, i.tol, N.clamp(i.upwind, 0, 1), ctx?.progress);
    const yc = N.range(n, (j) => (j + 0.5) * r.h), xn = N.linspace(0, 1, n + 1), [yMin, uMin] = vertex(yc, r.uc, N.argmin(r.uc)), [xvMax, vMax] = vertex(yc, r.vc, N.argmax(r.vc)), [xvMin, vMin] = vertex(yc, r.vc, N.argmin(r.vc));
    let jm = 0, im = 0; for (let j = 0; j <= n; j++) for (let k = 0; k <= n; k++) if (r.psi[j][k] < r.psi[jm][im]) { jm = j; im = k; }
    const [xv, p1] = vertex(xn, r.psi[jm], im), [yv, p2] = vertex(xn, r.psi.map((row) => row[im]), jm), psiMin = Math.min(p1, p2), U = (i.Re * nu) / i.L_m, warnings = [];
    if (!Number.isFinite(r.res)) warnings.push('The time march diverged: reduce the time-step safety factor or add upwind blending.');
    else if (r.res > i.tol) warnings.push(`Steady state was not reached within t = ${i.t_max} L/U (residual ${r.res.toExponential(2)}): increase the maximum time. At high Reynolds number the flow may be genuinely unsteady.`);
    if (i.Re * r.h > 2 && i.upwind < 0.2) warnings.push(`Cell Reynolds number is ${(i.Re * r.h).toFixed(1)} (> 2): central differencing can produce grid-scale oscillations near the lid corners; refine the grid or add upwind blending.`);
    if (i.Re > 3000) warnings.push('Above Re ≈ 3000 secondary eddies need at least 96 cells per side to resolve; treat this grid as qualitative.');
    return {
      kpis: [
        { key: 'u_min', label: 'Minimum u on the vertical centreline', value: uMin, unit: 'U', note: `at y/L = ${yMin.toFixed(4)}` },
        { key: 'y_u_min', label: 'Height of minimum u', value: yMin, unit: 'L' },
        { key: 'v_max', label: 'Maximum v on the horizontal centreline', value: vMax, unit: 'U', note: `at x/L = ${xvMax.toFixed(4)}` },
        { key: 'v_min', label: 'Minimum v on the horizontal centreline', value: vMin, unit: 'U', note: `at x/L = ${xvMin.toFixed(4)}` },
        { key: 'psi_min', label: 'Primary-vortex stream function', value: psiMin, unit: 'U·L' },
        { key: 'x_vortex', label: 'Primary vortex centre x', value: xv, unit: 'L' },
        { key: 'y_vortex', label: 'Primary vortex centre y', value: yv, unit: 'L' },
        { key: 'residual', label: 'Final steady-state residual', value: r.res, unit: 'U²/L', status: r.res <= i.tol ? 'ok' : 'warn' },
        { key: 'div_max', label: 'Largest cell divergence', value: r.div, unit: 'U/L', status: r.div < 1e-8 ? 'ok' : 'warn', note: 'Mass conservation after projection' },
        { key: 'steps', label: 'Time steps', value: r.steps, unit: '' },
        { key: 'U_lid_ms', label: 'Equivalent lid speed in air', value: U, unit: 'm/s', note: `for a ${i.L_m} m cavity at the case altitude` },
      ],
      plots: [
        { type: 'heat', title: `Stream function, Re = ${i.Re}`, xlabel: 'x/L [-]', ylabel: 'y/L [-]', zlabel: 'ψ/(U·L) [-]', x: xn, y: xn, z: r.psi, contours: 16, equalAspect: true },
        { type: 'heat', title: 'Velocity magnitude', xlabel: 'x/L [-]', ylabel: 'y/L [-]', zlabel: '|V|/U [-]', x: yc, y: yc, z: r.speed, contours: 12, equalAspect: true },
        { type: 'line', title: 'Centreline velocity profiles', xlabel: 'Position along the centreline [-]', ylabel: 'Velocity / U [-]', series: [{ name: 'u(y) at x = 0.5', x: [0, ...yc, 1], y: [0, ...r.uc, 1] }, { name: 'v(x) at y = 0.5', x: [0, ...yc, 1], y: [0, ...r.vc, 0] }] },
        { type: 'line', title: 'Convergence to steady state', xlabel: 'Time [L/U]', ylabel: 'Residual max|Δu|/Δt [U²/L]', ylog: true, series: [{ name: 'Residual', x: r.hist.t, y: r.hist.r }], annotations: [{ y: i.tol, label: 'Tolerance' }] },
        { type: 'heat', title: 'Pressure', xlabel: 'x/L [-]', ylabel: 'y/L [-]', zlabel: 'p/(ρU²) [-]', x: yc, y: yc, z: r.p.map((row) => row.map((v) => N.clamp(v - r.p[n >> 1][n >> 1], -0.3, 0.5))), contours: 14, equalAspect: true, diverging: true },
      ],
      warnings,
      models: ['Incompressible Navier–Stokes, primitive variables', 'Chorin projection method, explicit Euler in time', 'Staggered (MAC) finite-volume grid, second-order central differences with optional upwind blending', 'Direct pressure Poisson solver by discrete cosine eigen-decomposition'],
      assumptions: ['Laminar, two-dimensional, constant-property flow', 'Unit square with no-slip walls and a uniformly moving lid', 'Pressure contours are clipped near the singular lid corners'],
    };
  },
  convergence: { param: 'n', label: 'Cells per side', levels: [16, 24, 32, 40], metric: 'u_min' },
  verify() {
    const n = 32, h = 1 / n, P = new Float64Array(n * n), rhs = new Float64Array(n * n); let e = 0;
    for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) rhs[j * n + k] = -2 * PI * PI * Math.cos(PI * (k + 0.5) * h) * Math.cos(PI * (j + 0.5) * h);
    poissonNeumann(n, h)(P, rhs);
    for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) e = Math.max(e, Math.abs(P[j * n + k] - Math.cos(PI * (k + 0.5) * h) * Math.cos(PI * (j + 0.5) * h)));
    const o = N.kv(cavityFlow.run({ Re: 100, L_m: 0.5, n: 32, cfl: 0.8, t_max: 40, tol: 2e-5, upwind: 0 }));
    return [
      N.check('Poisson solver, manufactured solution cos(πx)cos(πy)', 1 + e, 1, 2e-3, 'Second-order truncation error π²h²/12 at h = 1/32'),
      N.check('Cavity Re = 100: minimum centreline u', o.u_min, -0.2109, 0.03, 'Ghia, Ghia & Shin (1982) benchmark; 32² second-order grid'),
      N.check('Cavity Re = 100: primary-vortex stream function', o.psi_min, -0.103423, 0.03, 'Ghia, Ghia & Shin (1982)'),
      N.check('Cavity Re = 100: vortex centre height', o.y_vortex, 0.7344, 0.03, 'Ghia, Ghia & Shin (1982)'),
      N.check('Discrete mass conservation', 1 + o.div_max, 1, 1e-10, 'Projection enforces a divergence-free field'),
    ];
  },
  validation: [{ name: 'Lid-driven cavity, minimum centreline velocity', source: 'Ghia, Ghia & Shin, J. Comput. Phys. 48 (1982): u_min = −0.2109 (Re 100), −0.3829 (Re 1000)', inputs: { n: 64, t_max: 60, upwind: 0 }, sweep: { key: 'Re', values: [100] }, target: 'u_min', observed: [-0.2109], tol_pct: 2 }],
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.residual > i.tol) out.push({ severity: 'warn', title: 'Solution is not converged in time', detail: `Residual ${o.residual.toExponential(2)} against a tolerance of ${i.tol}.`, action: 'Increase the maximum time; do not read forces or profiles from an unconverged field.', basis: 'Iterative convergence' });
    out.push({ severity: 'info', title: 'Quantify the grid error before using the numbers', detail: `Centreline minimum ${o.u_min.toFixed(4)} on a ${Math.round(i.n)}² grid.`, action: 'Run the convergence study: three grids give the observed order and a grid-convergence index. The same discipline applies to any external RANS solution brought into this platform.', basis: 'Richardson extrapolation / GCI (Roache)' });
    if (i.Re > 1000) out.push({ severity: 'advise', title: 'Open cavities are noise and drag sources', detail: 'At flight Reynolds numbers cavity shear layers are turbulent and can lock into acoustic (Rossiter) tones.', action: 'For gear bays and gaps, assess tones in Suite 11 and consider fairings or leading-edge spoilers; they reduce both airframe noise and drag.', basis: 'Cavity flow physics' });
    return out;
  },
};

// ---- wall turbulence: first-cell height and a 1-D turbulence-model solve ------------------------
/**
 * Fully developed turbulent channel flow in wall units, 0 <= y+ <= Reτ, total stress (1 + νt+) du+/dy+ = 1 − y+/Reτ.
 * model: 'Spalart–Allmaras' solves the ν̃ transport equation; otherwise the van Driest / Nikuradse mixing length.
 */
export function wallProfile(ReTau, model, n, kappa = 0.41, Aplus = 26) {
  const y1 = Math.min(0.2, ReTau / (4 * n)), r = N.brent((q) => (y1 * (q ** n - 1)) / (q - 1) - ReTau, 1 + 1e-9, 3, 1e-13), y = N.range(n + 1, (j) => (y1 * (r ** j - 1)) / (r - 1));
  y[n] = ReTau;
  const tau = y.map((v) => 1 - v / ReTau), nut = new Array(n + 1).fill(0); let iters = 0, change = 0;
  const ml = y.map((v, j) => { const e = 1 - v / ReTau, l = ReTau * (0.14 - 0.08 * e * e - 0.06 * e ** 4) * (kappa / 0.4) * (1 - Math.exp(-v / Aplus)); return l * l * (l > 0 ? (-1 + Math.sqrt(1 + 4 * l * l * tau[j])) / (2 * l * l) : 0); });
  if (model !== 'Spalart–Allmaras') for (let j = 0; j <= n; j++) nut[j] = ml[j];
  else {
    const cb1 = 0.1355, sg = 2 / 3, cb2 = 0.622, cv1 = 7.1, cw1 = cb1 / (kappa * kappa) + (1 + cb2) / sg, cw2 = 0.3, cw3 = 2;
    let nt = y.map((v) => kappa * v * (1 - (0.5 * v) / ReTau));
    const a = new Array(n + 1).fill(0), b = new Array(n + 1).fill(1), c = new Array(n + 1).fill(0), d = new Array(n + 1).fill(0);
    for (iters = 0; iters < 4000; iters++) {
      for (let j = 1; j < n; j++) {
        const chi = Math.max(nt[j], 1e-12), fv1 = chi ** 3 / (chi ** 3 + cv1 ** 3), fv2 = 1 - chi / (1 + chi * fv1), S = tau[j] / (1 + chi * fv1), k2d2 = kappa * kappa * y[j] * y[j];
        const St = Math.max(S + (chi * fv2) / k2d2, 0.3 * S), rr = Math.min(chi / (St * k2d2), 10), gg = rr + cw2 * (rr ** 6 - rr), fw = gg * ((1 + cw3 ** 6) / (gg ** 6 + cw3 ** 6)) ** (1 / 6);
        const hm = y[j] - y[j - 1], hp = y[j + 1] - y[j], hc = 0.5 * (hm + hp), // diffusion in the form (1+cb2)∇·((1+ν̃)∇ν̃) − cb2(1+ν̃)∇²ν̃, which keeps the tridiagonal operator diagonally dominant
        dm = ((1 + cb2) * (1 + 0.5 * (nt[j] + nt[j - 1])) - cb2 * (1 + nt[j])) / (sg * hm * hc), dp = ((1 + cb2) * (1 + 0.5 * (nt[j] + nt[j + 1])) - cb2 * (1 + nt[j])) / (sg * hp * hc);
        const dest = (cw1 * fw * chi) / (y[j] * y[j]), dti = cb1 * St + dest; // pseudo-time step scaled on the local source time scale
        a[j] = -dm; c[j] = -dp; b[j] = dm + dp + 2 * dest + dti; d[j] = cb1 * St * chi + dest * chi + dti * chi;
      }
      b[0] = 1; c[0] = 0; d[0] = 0; a[n] = -1; b[n] = 1; d[n] = 0; // wall value and centreline symmetry
      const nn = N.solveTridiag(a, b, c, d); change = 0;
      for (let j = 0; j <= n; j++) { const v = Math.max(nn[j], 0); change = Math.max(change, Math.abs(v - nt[j]) / (1 + nt[j])); nt[j] += 0.5 * (v - nt[j]); } // under-relaxed update
      if (change < 1e-10) break;
    }
    for (let j = 0; j <= n; j++) { const chi = nt[j]; nut[j] = (chi * chi ** 3) / (chi ** 3 + cv1 ** 3 || 1); }
  }
  const S = tau.map((v, j) => v / (1 + nut[j])), u = N.cumtrapz(y, S), Ub = N.trapz(y, u) / ReTau;
  return { y, u, nut, S, Ub, Ucl: u[n], iters, change, ml };
}
const cfWhite = (Rex) => 0.455 / Math.log(0.06 * Rex) ** 2;
const wallTurb = {
  id: 'wallturb', title: 'Turbulent wall layer: first-cell height and log-law solve', fidelity: 'numerical',
  summary: 'Skin friction, boundary-layer thickness and the near-wall mesh spacing needed for a chosen y⁺, together with a one-dimensional Spalart–Allmaras (or mixing-length) solution of the turbulent wall layer that reproduces the viscous sublayer and the logarithmic law.',
  equations: ['Reynolds-averaged Navier–Stokes equations', 'Prandtl boundary-layer equations', 'Boundary-layer equations', 'Turbulent kinetic energy transport equation', 'Spalart–Allmaras'],
  inputs: [
    { key: 'V', label: 'Freestream speed', unit: 'm/s', default: 70, min: 0.5, group: 'Flow' },
    { key: 'x_m', label: 'Distance from the leading edge', unit: 'm', default: 1.5, min: 1e-3, group: 'Flow', help: 'Use the mean chord for a wing or the fuselage length for a body' },
    ...ATM,
    { key: 'yplus', label: 'Target y⁺ of the first cell centre', unit: '-', default: 1, min: 0.1, max: 500, group: 'Mesh', help: '≈1 for wall-resolved SA or k–ω SST; 30–300 with wall functions' },
    { key: 'growth', label: 'Prism-layer growth ratio', unit: '-', default: 1.2, min: 1.05, max: 1.5, group: 'Mesh' },
    { key: 'model', label: 'Turbulence closure', type: 'select', options: ['Spalart–Allmaras', 'Mixing length (van Driest)'], default: 'Spalart–Allmaras', group: 'Model' },
    { key: 'Re_tau', label: 'Friction Reynolds number (0 = from the boundary layer)', unit: '-', default: 0, min: 0, max: 50000, group: 'Model', help: 'δ·uτ/ν of the 1-D solve' },
    { key: 'kappa', label: 'Von Kármán constant', unit: '-', default: 0.41, min: 0.38, max: 0.43, group: 'Model' },
    { key: 'A_plus', label: 'Van Driest damping constant', unit: '-', default: 26, min: 20, max: 30, group: 'Model', help: 'Mixing-length model only' },
    { key: 'nPts', label: 'Grid points across the layer', unit: '', default: 120, min: 30, max: 600, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => ({ V: c.wing.S_m2 > 0 ? c.flight.V_ms : Math.max(c.flight.V_ms, 0.75 * d.v_tip) || undefined, x_m: (c.wing.S_m2 > 0 ? d.mac : c.rotor.chord_m) || undefined, alt_m: c.atm.alt_m, dISA: c.atm.dISA_K }),
  run(i) {
    const a = isa(i.alt_m, i.dISA), Rex = Math.max((i.V * i.x_m) / a.nu, 10), turb = Rex > 5e5, cf = turb ? cfWhite(Rex) : 0.664 / Math.sqrt(Rex), tw = 0.5 * a.rho * i.V ** 2 * cf, ut = Math.sqrt(tw / a.rho);
    const y1 = (2 * i.yplus * a.nu) / ut, d99 = turb ? (0.37 * i.x_m) / Rex ** 0.2 : (5 * i.x_m) / Math.sqrt(Rex), nLay = Math.ceil(Math.log(1 + (d99 * (i.growth - 1)) / y1) / Math.log(i.growth));
    const ReT = N.clamp(i.Re_tau > 0 ? i.Re_tau : (d99 * ut) / a.nu, 100, 50000), n = Math.round(N.clamp(i.nPts, 30, 600)), w = wallProfile(ReT, i.model, n, i.kappa, i.A_plus);
    const lg = w.y.map((v, j) => j).filter((j) => w.y[j] > 40 && w.y[j] < 0.15 * ReT), fit = lg.length > 3 ? N.polyfit(lg.map((j) => Math.log(w.y[j])), lg.map((j) => w.u[j]), 1) : [NaN, NaN];
    const k1 = w.y.findIndex((v) => v >= 100), u100 = k1 > 0 ? N.interp1(w.y, w.u, 100) : NaN, js = w.y.findIndex((v) => v > 0.5), warnings = [];
    if (!turb) warnings.push(`Re_x = ${Rex.toExponential(2)} is below 5×10⁵: the boundary layer is probably laminar here and the mesh estimate uses the Blasius skin friction; the wall-turbulence profile is shown for a tripped layer.`);
    if (i.yplus > 5 && i.yplus < 30) warnings.push('A first-cell y⁺ between 5 and 30 lies in the buffer layer, where neither wall-resolved nor wall-function treatments are accurate.');
    if (i.model === 'Spalart–Allmaras' && w.change > 1e-6) warnings.push(`The Spalart–Allmaras iteration stopped at a relative change of ${w.change.toExponential(1)}; add grid points.`);
    if (lg.length <= 3) warnings.push('The friction Reynolds number is too low for a distinct logarithmic region; slope and intercept are not fitted.');
    const yp = N.logspace(0.3, Math.max(ReT, 100), 60), ypS = yp.filter((v) => v < 12), ypL = yp.filter((v) => v > 8);
    const xs = N.logspace(Math.max(i.x_m / 200, 1e-4), i.x_m, 60), cfx = xs.map((x) => { const R = (i.V * x) / a.nu; return R > 5e5 ? cfWhite(R) : 0.664 / Math.sqrt(R); });
    return {
      kpis: [
        { key: 'Re_x', label: 'Reynolds number at the station', value: Rex, unit: '-' },
        { key: 'cf_local', label: 'Local skin-friction coefficient', value: cf, unit: '-', note: turb ? 'White turbulent flat-plate correlation (empirical)' : 'Blasius laminar' },
        { key: 'tau_w_Pa', label: 'Wall shear stress', value: tw, unit: 'Pa' },
        { key: 'u_tau_ms', label: 'Friction velocity', value: ut, unit: 'm/s' },
        { key: 'first_cell_height_m', label: `First-cell height for y⁺ = ${i.yplus}`, value: y1, unit: 'm', note: 'Cell height; the centre is at half of it' },
        { key: 'delta99_m', label: 'Boundary-layer thickness', value: d99, unit: 'm' },
        { key: 'n_prism_layers', label: 'Prism layers to cover the boundary layer', value: nLay, unit: '', status: nLay <= 60 ? 'ok' : 'warn' },
        { key: 'y_wallfn_m', label: 'First-cell height for wall functions (y⁺ = 50)', value: (100 * a.nu) / ut, unit: 'm' },
        { key: 'Re_tau', label: 'Friction Reynolds number of the solve', value: ReT, unit: '-' },
        { key: 'kappa_fit', label: 'Fitted log-law slope constant κ', value: 1 / fit[1], unit: '-', status: Math.abs(1 / fit[1] - 0.41) < 0.03 ? 'ok' : 'warn' },
        { key: 'B_fit', label: 'Fitted log-law intercept B', value: fit[0], unit: '-', note: 'Accepted range 5.0–5.5' },
        { key: 'u_plus_100', label: 'u⁺ at y⁺ = 100', value: u100, unit: '-' },
        { key: 'sublayer_ratio', label: 'u⁺/y⁺ in the viscous sublayer', value: js > 0 ? w.u[js] / w.y[js] : NaN, unit: '-' },
        { key: 'Ucl_plus', label: 'Centreline velocity U⁺', value: w.Ucl, unit: '-' },
        { key: 'cf_channel', label: 'Channel friction coefficient 2/U_b⁺²', value: 2 / w.Ub ** 2, unit: '-' },
      ],
      plots: [
        { type: 'line', title: 'Law of the wall', xlabel: 'y⁺ [-]', ylabel: 'u⁺ [-]', xlog: true, series: [{ name: i.model, x: w.y.slice(1), y: w.u.slice(1) }, { name: 'Sublayer u⁺ = y⁺', x: ypS, y: ypS, style: 'dash' }, { name: `Log law (κ = ${i.kappa}, B = 5.0)`, x: ypL, y: ypL.map((v) => Math.log(v) / i.kappa + 5.0), style: 'dash' }], annotations: [{ x: i.yplus, label: 'First cell' }] },
        { type: 'line', title: 'Eddy-viscosity ratio across the layer', xlabel: 'y⁺ [-]', ylabel: 'νt/ν [-]', xlog: true, series: [{ name: i.model, x: w.y.slice(1), y: w.nut.slice(1) }, { name: 'Mixing length (reference)', x: w.y.slice(1), y: w.ml.slice(1), style: 'dash' }] },
        { type: 'line', title: 'Flat-plate skin friction along the surface', xlabel: 'Distance from leading edge [m]', ylabel: 'cf [-]', xlog: true, ylog: true, series: [{ name: 'Local cf', x: xs, y: cfx }], annotations: [{ x: i.x_m, label: 'Station' }] },
      ],
      tables: [{ title: 'Near-wall mesh guide at this station', columns: ['Wall treatment', 'Target y⁺', 'First-cell height [µm]'], rows: [['Wall-resolved (SA, k–ω SST, transition models)', '≈ 1', +((2 * a.nu) / ut * 1e6).toFixed(2)], ['Wall-resolved LES', '< 1 (Δx⁺ ≈ 50, Δz⁺ ≈ 15)', +((1.6 * a.nu) / ut * 1e6).toFixed(2)], ['Wall functions / wall-modelled LES', '30–300', `${((60 * a.nu) / ut * 1e6).toFixed(0)}–${((600 * a.nu) / ut * 1e6).toFixed(0)}`]] }],
      outputs: { sa_iterations: w.iters },
      warnings,
      models: [i.model === 'Spalart–Allmaras' ? 'Spalart–Allmaras one-equation model (standard constants, no trip term), 1-D fully developed wall layer' : 'Prandtl mixing length with van Driest damping and Nikuradse outer length', 'Flat-plate skin-friction correlations for the mesh estimate (empirical)', 'Seventh-power-law boundary-layer thickness (empirical)'],
      assumptions: ['Zero pressure gradient, smooth wall, fully turbulent layer', 'The 1-D solve is a fully developed channel at the friction Reynolds number of the boundary layer: its wall region is universal, its outer region is not that of an aircraft boundary layer', 'First-cell height is the full cell height with the centre at the target y⁺'],
    };
  },
  convergence: { param: 'nPts', label: 'Grid points across the wall layer', levels: [40, 80, 160, 320], metric: 'Ucl_plus' },
  calibration: { params: [{ key: 'kappa', min: 0.38, max: 0.43 }, { key: 'A_plus', min: 22, max: 30 }], sweep: 'Re_tau', target: 'Ucl_plus', note: 'Supply measured centreline or edge velocity U⁺ against friction Reynolds number. Only the mixing-length constants are open to calibration; the Spalart–Allmaras constants are universal and should not be tuned to one case.' },
  verify() {
    const sa = wallProfile(5000, 'Spalart–Allmaras', 200), mlp = wallProfile(5000, 'Mixing length (van Driest)', 200), at = (w, yp) => N.interp1(w.y, w.u, yp);
    const slope = (w) => (at(w, 400) - at(w, 100)) / Math.log(4), bal = Math.max(...sa.y.map((v, j) => Math.abs((1 + sa.nut[j]) * sa.S[j] - (1 - v / 5000))));
    return [
      N.check('Viscous sublayer u⁺ = y⁺ (SA)', at(sa, 1), 1, 5e-3, 'Exact near-wall limit'),
      N.check('Log-law slope 1/κ (SA)', slope(sa), 1 / 0.41, 0.03, 'Spalart–Allmaras reproduces the log law by construction'),
      N.check('u⁺ at y⁺ = 100 (SA)', at(sa, 100), Math.log(100) / 0.41 + 5.0, 0.03, 'Log law with B = 5.0 (accepted scatter 5.0–5.5)'),
      N.check('Log-law slope 1/κ (mixing length)', slope(mlp), 1 / 0.41, 0.04, 'Prandtl mixing length'),
      N.check('u⁺ at y⁺ = 100 (van Driest)', at(mlp, 100), Math.log(100) / 0.41 + 5.2, 0.04, 'van Driest profile, B ≈ 5.2'),
      N.check('Total-stress balance residual', 1 + bal, 1, 1e-12, 'Integrated momentum equation'),
      N.check('SA eddy viscosity ν̃ = κ·y in the log layer', N.interp1(sa.y, sa.nut, 150) / (0.41 * 150), 1, 0.06, 'Spalart–Allmaras log-layer solution (νt ≈ ν̃ for χ ≫ 7)'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    out.push({ severity: 'info', title: 'Near-wall mesh for an external RANS run', detail: `First cell ${(o.first_cell_height_m * 1e6).toFixed(1)} µm for y⁺ = ${i.yplus}, ${o.n_prism_layers} layers at growth ${i.growth} to cover δ = ${(o.delta99_m * 1e3).toFixed(1)} mm.`, action: 'Use these values to set the prism layer when preparing the hand-off mesh; check y⁺ after the first solution and adapt.', basis: 'Flat-plate estimate of wall shear' });
    if (o.n_prism_layers > 45) out.push({ severity: 'advise', title: 'Wall-resolved mesh is expensive at this Reynolds number', detail: `${o.n_prism_layers} prism layers are needed.`, action: 'Consider wall functions (y⁺ 30–300) for attached-flow drag studies and keep wall resolution for separation, transition and heat transfer.', basis: 'Mesh cost versus near-wall accuracy' });
    if (i.yplus > 5 && i.yplus < 30) out.push({ severity: 'warn', title: 'First cell in the buffer layer', detail: `y⁺ = ${i.yplus}.`, action: 'Move the first cell to y⁺ ≈ 1 or above 30.', basis: 'Law of the wall' });
    return out;
  },
};

// ---- native 3-D Navier–Stokes: DNS, LES and RANS on Cartesian immersed-boundary grids ----------
const TWO_PI = 2 * PI;
const cellsOf = (n, lo = 8, hi = 128) => Math.round(N.clamp(n, lo, hi));
const tailOf = (a, f = 0.25) => a.slice(Math.max(0, Math.min(a.length - 1, Math.floor(a.length * (1 - f)))));
const heatOf = (title, sl, zlabel, o = {}) => ({ type: 'heat', title, xlabel: o.xlabel || 'x [-]', ylabel: o.ylabel || 'y [-]', zlabel, x: sl.x, y: sl.y, z: o.map ? sl.z.map((r) => r.map(o.map)) : sl.z, contours: o.contours || 14, equalAspect: o.equalAspect !== false, ...(o.diverging ? { diverging: true } : {}), ...(o.overlay && o.overlay.x.length ? { overlay: [{ name: o.overlayName || 'Body outline', x: o.overlay.x, y: o.overlay.y }] } : {}) });
/** Keep the last result of an expensive, input-determined run so that re-opening it (or running it for another aircraft with identical inputs) is free. */
function memoLast(fn) {
  let key = null, val = null;
  return (i, ctx, extra = '') => { const k = JSON.stringify(i) + extra; if (k !== key || !val) { val = fn(i, ctx); key = k; } return typeof structuredClone === 'function' ? structuredClone(val) : JSON.parse(JSON.stringify(val)); };
}
const COARSE_NOTE = (what, n) => `The default grid (${what}) is deliberately small so that the run returns in about a second on any device. Treat it as a first look: raise "${n}" (up to 128) and use the mesh-convergence tab before quoting numbers.`;
const TGV_FLOWS = ['3-D Taylor–Green vortex', '2-D Taylor–Green vortex (exact solution)'];
const tgv3 = (x, y, z) => [Math.sin(x) * Math.cos(y) * Math.cos(z), -Math.cos(x) * Math.sin(y) * Math.cos(z), 0];
const tgv2 = (x, y) => [Math.sin(x) * Math.cos(y), -Math.cos(x) * Math.sin(y), 0];
/** Taylor–Green vortex in a (2π)³ periodic box with unit velocity and length scales; fixed time step that lands exactly on tEnd. */
function runTgv(o) {
  const n = o.n, h = TWO_PI / n, nu = 1 / o.Re, rk = o.time !== 'ab2', cfl = N.clamp(o.cfl ?? (rk ? 1 : 0.3), 0.05, rk ? 1.5 : 0.5);
  const steps = Math.max(2, Math.ceil(o.tEnd / Math.min((cfl * h) / 1.6, ((rk ? 0.5 : 0.2) * h * h) / (3 * nu))));
  const s = createSolver({ n: [n, n, n], L: [TWO_PI, TWO_PI, TWO_PI], nu, init: o.two ? tgv2 : tgv3, time: rk ? 'rk3' : 'ab2', dt: o.tEnd / steps, model: o.model || 'dns', sgs: o.sgs, Cs: o.Cs, Cw: o.Cw, vanDriest: false, historyEvery: 1 });
  s.run(steps, o.progress ? (f) => o.progress(f, `t = ${(f * o.tEnd).toFixed(2)} L/U`) : null, Math.max(1, Math.round(steps / 20)));
  return s;
}
/** Root-mean-square velocity error against the exact decaying 2-D Taylor–Green solution. */
function tgv2Error(s, nu) {
  const g = s.grid, f = s.fields, dec = Math.exp(-2 * nu * s.time); let e = 0;
  for (let k = 0; k < g.nz; k++) for (let j = 0; j < g.ny; j++) for (let i = 0; i < g.nx; i++) { const I = s.index(i, j, k); e += (f.u[I] - tgv2(g.xf[i + 1], g.y[j])[0] * dec) ** 2 + (f.v[I] - tgv2(g.x[i], g.yf[j + 1])[1] * dec) ** 2 + f.w[I] ** 2; }
  return Math.sqrt(e / (g.nx * g.ny * g.nz));
}
/** Kinetic-energy budget from a solver history: −dE/dt between samples against the mean resolved + modelled dissipation. */
function energyBudget(H) {
  const t = [], dE = [], eps = []; let num = 0, den = 0;
  for (let k = 0; k + 1 < H.t.length; k++) { const dt = H.t[k + 1] - H.t[k]; if (!(dt > 0)) continue; t.push(0.5 * (H.t[k] + H.t[k + 1])); dE.push((H.ke[k] - H.ke[k + 1]) / dt); eps.push(0.5 * (H.dissipation[k] + H.dissipation[k + 1] + H.sgs[k] + H.sgs[k + 1])); num += H.ke[k] - H.ke[k + 1]; den += eps[eps.length - 1] * dt; }
  return { t, dE, eps, ratio: den > 0 ? num / den : NaN };
}

const dns3d = {
  id: 'dns3d', title: '3-D DNS: Taylor–Green vortex in a periodic box', fidelity: 'numerical',
  summary: 'Direct numerical simulation of the three-dimensional incompressible Navier–Stokes equations with no turbulence model: a Taylor–Green vortex stretches, breaks down and decays in a periodic box. Energy decay, dissipation rate, vorticity and Q-criterion fields are computed and the solver is verified against the exact viscous solution and the energy budget.',
  equations: ['Navier–Stokes equations', 'Continuity equation', 'Conservation of momentum equation', 'Vorticity transport equation', 'DNS'],
  inputs: [
    { key: 'flow', label: 'Flow', type: 'select', options: TGV_FLOWS, default: TGV_FLOWS[0], group: 'Flow', help: 'The 2-D vortex has an exact solution (pure viscous decay) and measures the solver error directly; the 3-D vortex transitions to turbulence at high Reynolds number' },
    { key: 'Re', label: 'Reynolds number U·L/ν', unit: '-', default: 100, min: 0.01, max: 5000, group: 'Flow', help: 'Box side is 2πL. 100 is laminar and resolvable on 24³–32³; 1600 is the standard turbulence benchmark and needs 256³ or more for a true DNS' },
    { key: 't_end', label: 'Simulated time', unit: 'L/U', default: 5, min: 0.05, max: 40, group: 'Flow', help: 'The 3-D dissipation peak is near t ≈ 5 at Re 100 and t ≈ 9 at Re 1600' },
    { key: 'L_m', label: 'Length scale L', unit: 'm', default: 0.05, min: 1e-4, group: 'Flow', help: 'Only used to give the equivalent speed in air at the case altitude' },
    { key: 'n', label: 'Cells per side', unit: '', default: 16, min: 8, max: 128, step: 1, discrete: true, group: 'Numerics', help: 'Powers of two (16, 32, 64, 128) use the FFT pressure solver and are fastest' },
    { key: 'scheme', label: 'Time integration', type: 'select', options: ['Runge–Kutta 3', 'Adams–Bashforth 2'], default: 'Runge–Kutta 3', group: 'Numerics' },
    { key: 'cfl', label: 'CFL number', unit: '-', default: 1, min: 0.05, max: 1.5, group: 'Numerics', help: 'Runge–Kutta is stable to about 1.5; Adams–Bashforth is limited to 0.5 internally' },
  ],
  defaults: () => ({}),
  run(i, ctx) { return dnsMemo(i, ctx, String(ctx?.case?.atm?.alt_m ?? 0) + '/' + String(ctx?.case?.atm?.dISA_K ?? 0)); },
  runNow(i, ctx) {
    const n = cellsOf(i.n), two = i.flow === TGV_FLOWS[1], Re = Math.max(i.Re, 1e-6), nu = 1 / Re, E0 = two ? 0.25 : 0.125;
    const s = runTgv({ n, Re, tEnd: i.t_end, two, time: i.scheme === 'Adams–Bashforth 2' ? 'ab2' : 'rk3', cfl: i.cfl, progress: ctx?.progress });
    const H = s.history, d = s.diagnostics(), B = energyBudget(H), kMax = N.argmax(H.dissipation), epsMax = H.dissipation[kMax], eta = (nu ** 3 / Math.max(epsMax, 1e-300)) ** 0.25, h = TWO_PI / n, hEta = h / eta;
    const budgetErr = N.amax(B.dE.map((v, k) => Math.abs(v - B.eps[k]))) / Math.max(epsMax, 1e-300), Eex = E0 * Math.exp(-4 * nu * s.time), err2 = two ? tgv2Error(s, nu) : NaN;
    const nuAir = isa(ctx?.case?.atm?.alt_m ?? 0, ctx?.case?.atm?.dISA_K ?? 0).nu, warnings = [];
    if (s.diverged || !Number.isFinite(d.ke)) warnings.push('The time march diverged: lower the CFL number.');
    if (n <= 20) warnings.push(COARSE_NOTE(`${n}³ cells`, 'Cells per side'));
    if (hEta > 2.1) warnings.push(`Not a resolved DNS at this Reynolds number: the cell size is ${hEta.toFixed(1)} Kolmogorov lengths (a DNS needs about 2 or less, i.e. roughly ${Math.ceil((n * hEta) / 2.1)}³ cells here). The energy-conserving central scheme has no numerical dissipation, so unresolved energy piles up at the grid scale and the dissipation peak is under-predicted. Use the LES analysis or refine.`);
    if (!two && Re >= 1000 && n < 128) warnings.push('For orientation only: published pseudo-spectral DNS of the Re = 1600 Taylor–Green vortex (512³) gives a peak dissipation of about 0.0127 U³/L near t ≈ 9 L/U. A second-order scheme on a grid this coarse cannot reproduce it.');
    if (!two && i.t_end < 3) warnings.push('The run stops before vortex stretching has built up the small scales: the dissipation maximum has not been reached.');
    const [tt, ke] = thin(H.t, H.ke, 300), [, ep] = thin(H.t, H.dissipation, 300), [tb, db] = thin(B.t, B.dE, 300), k0 = 0, sv = s.slice('z', k0, 'vorticity', { max: 80 }), sq = s.slice('z', k0, 'q', { max: 80 }), so = s.slice('y', Math.floor(n / 4), 'wy', { max: 80 });
    const kpis = [
      { key: 'E_final', label: 'Kinetic energy at the end', value: d.ke, unit: 'U²', note: `${(100 * d.ke / E0).toFixed(1)}% of the initial energy` },
      ...(two ? [
        { key: 'E_exact', label: 'Exact kinetic energy', value: Eex, unit: 'U²', note: '¼·exp(−4νt)' },
        { key: 'ke_error', label: 'Relative energy error', value: Math.abs(d.ke - Eex) / Eex, unit: '-', status: Math.abs(d.ke - Eex) / Eex < 0.01 ? 'ok' : 'warn' },
        { key: 'u_error_rms', label: 'RMS velocity error against the exact solution', value: err2, unit: 'U' },
      ] : [
        { key: 'eps0_ratio', label: 'Initial dissipation / exact value ¾ν', value: H.dissipation[0] / (0.75 * nu), unit: '-', note: 'Differs from 1 by the second-order truncation error of the grid' },
      ]),
      { key: 'eps_max', label: 'Peak dissipation rate', value: epsMax, unit: 'U³/L', note: `at t = ${H.t[kMax].toFixed(2)} L/U${kMax === H.t.length - 1 ? ' (still rising at the end of the run)' : ''}` },
      { key: 't_eps_max', label: 'Time of peak dissipation', value: H.t[kMax], unit: 'L/U' },
      { key: 'enstrophy_max', label: 'Peak enstrophy', value: N.amax(H.enstrophy), unit: 'U²/L²' },
      { key: 'budget_ratio', label: 'Energy lost / time-integrated dissipation', value: B.ratio, unit: '-', status: Math.abs(B.ratio - 1) < 0.01 ? 'ok' : 'warn', note: 'Exactly 1 for the continuous equations: dE/dt = −ε' },
      { key: 'budget_error', label: 'Largest energy-budget imbalance', value: budgetErr, unit: 'of peak ε', status: budgetErr < 0.02 ? 'ok' : 'warn' },
      { key: 'dx_over_eta', label: 'Cell size / Kolmogorov length', value: hEta, unit: '-', status: hEta <= 2.1 ? 'ok' : hEta < 4 ? 'warn' : 'bad', note: 'A resolved DNS needs about 2 or less' },
      { key: 'div_max', label: 'Largest cell divergence', value: d.divMax, unit: 'U/L', status: d.divMax < 1e-9 ? 'ok' : 'warn', note: 'Mass conservation after the exact projection' },
      { key: 'steps', label: 'Time steps', value: d.step, unit: '' },
      { key: 'dt', label: 'Time step', value: d.dt, unit: 'L/U' },
      { key: 'cells', label: 'Grid cells', value: n ** 3, unit: '' },
      { key: 'U_equiv_ms', label: 'Equivalent velocity scale in air', value: (Re * nuAir) / i.L_m, unit: 'm/s', note: `for L = ${i.L_m} m at the case altitude` },
    ];
    return {
      kpis,
      plots: [
        { type: 'line', title: 'Kinetic-energy decay', xlabel: 'Time [L/U]', ylabel: 'Kinetic energy [U²]', series: [{ name: `DNS ${n}³`, x: tt, y: ke }, two ? { name: 'Exact ¼·exp(−4νt)', x: tt, y: tt.map((t) => 0.25 * Math.exp(-4 * nu * t)), style: 'dash' } : { name: 'Viscous decay of the initial mode alone, ⅛·exp(−6νt)', x: tt, y: tt.map((t) => 0.125 * Math.exp(-6 * nu * t)), style: 'dash' }] },
        { type: 'line', title: 'Dissipation rate and energy budget', xlabel: 'Time [L/U]', ylabel: 'Rate [U³/L]', series: [{ name: 'ε = ν⟨|∇u|²⟩', x: tt, y: ep }, { name: '−dE/dt', x: tb, y: db, style: 'dash' }] },
        heatOf(`Vorticity magnitude, plane z = ${sv.position.toFixed(2)} L at t = ${s.time.toFixed(2)}`, sv, '|ω| [U/L]', { xlabel: 'x [L]', ylabel: 'y [L]' }),
        heatOf(`Q-criterion, plane z = ${sq.position.toFixed(2)} L`, sq, 'Q [U²/L²]', { xlabel: 'x [L]', ylabel: 'y [L]', diverging: true }),
        heatOf(`Vorticity component ω_y, plane y = ${so.position.toFixed(2)} L`, so, 'ω_y [U/L]', { xlabel: 'x [L]', ylabel: 'z [L]', diverging: true }),
      ],
      tables: [{ title: 'Run summary', columns: ['Quantity', 'Value'], rows: [['Grid', `${n} × ${n} × ${n}`], ['Box', '2πL × 2πL × 2πL, periodic'], ['Time integration', i.scheme], ['Cell size h / L', +h.toFixed(4)], ['Kolmogorov length η / L (at peak ε)', +eta.toFixed(4)], ['Pressure solver', (n & (n - 1)) === 0 ? 'FFT (three directions)' : 'dense Fourier transform (use a power of two for speed)']] }],
      warnings,
      models: ['Incompressible Navier–Stokes equations, no turbulence model (DNS)', 'Staggered (MAC) grid, second-order energy-conserving central differences', i.scheme === 'Adams–Bashforth 2' ? 'Adams–Bashforth 2 fractional-step projection' : 'Three-stage low-storage Runge–Kutta (Wray) fractional-step projection', 'Exact discrete pressure projection by Fourier transforms'],
      assumptions: ['Constant-property incompressible flow in a triply periodic box', 'Dissipation is the discrete ν⟨|∇u|²⟩, which equals ν⟨ω²⟩ for periodic flow', 'Second-order accuracy in space: about 2 to 3 times more points per direction are needed than with a spectral method for the same resolved scales'],
    };
  },
  convergence: { param: 'n', label: 'Cells per side', levels: [12, 16, 24, 32], metric: 'E_final', hOf: (n) => TWO_PI / n },
  verify() {
    const nu = 0.05, a = runTgv({ n: 8, Re: 1 / nu, tEnd: 1, two: true }), b = runTgv({ n: 16, Re: 1 / nu, tEnd: 1, two: true }), ea = tgv2Error(a, nu), eb = tgv2Error(b, nu), db = b.diagnostics();
    const st = runTgv({ n: 16, Re: 0.02, tEnd: 0.002 }), ds = st.diagnostics(), mod = (Math.sin(PI / 16) / (PI / 16)) ** 2; // Stokes limit: every initial mode has |k|² = 3
    const t3 = runTgv({ n: 16, Re: 100, tEnd: 1.5 }), B = energyBudget(t3.history), ab = runTgv({ n: 16, Re: 1 / nu, tEnd: 1, two: true, time: 'ab2' });
    return [
      N.check('2-D Taylor–Green: observed spatial order of accuracy', Math.log2(ea / eb), 2, 0.08, 'Exact Navier–Stokes solution u = sin x cos y·exp(−2νt); errors on 8³ and 16³'),
      N.check('2-D Taylor–Green: kinetic energy at t = 1 (16³)', db.ke, 0.25 * Math.exp(-4 * nu), 5e-3, 'Exact decay ¼·exp(−4νt); remaining error is the second-order modified wavenumber'),
      N.check('2-D Taylor–Green with Adams–Bashforth 2', ab.diagnostics().ke, db.ke, 2e-4, 'Two independent time integrators agree'),
      N.check('3-D Taylor–Green, Stokes limit: energy decay exp(−6νt)', ds.ke, 0.125 * Math.exp(-6 * 50 * 0.002 * mod), 2e-3, 'Linear viscous decay with the discrete (modified) wavenumber of the 16³ grid'),
      N.check('3-D Taylor–Green: initial dissipation ε₀ = ¾ν', t3.history.dissipation[0] / (0.75 * 0.01), mod, 1e-6, 'Analytical enstrophy ⟨ω²⟩ = ¾ of the initial field, times the modified-wavenumber factor'),
      N.check('Energy budget dE/dt = −ε (3-D, Re 100)', B.ratio, 1, 2e-3, 'Kinetic-energy equation; the central scheme adds no numerical dissipation'),
      N.check('Discrete mass conservation in 3-D', 1 + t3.diagnostics().divMax, 1, 1e-11, 'Exact projection'),
    ];
  },
  validation: [{ name: '2-D Taylor–Green vortex, kinetic energy at t = 2', source: 'Exact solution of the Navier–Stokes equations (Taylor 1923): E = ¼·exp(−4νt)', inputs: { flow: TGV_FLOWS[1], t_end: 2, n: 32 }, sweep: { key: 'Re', values: [10, 20, 100] }, target: 'E_final', observed: [0.25 * Math.exp(-0.8), 0.25 * Math.exp(-0.4), 0.25 * Math.exp(-0.08)], tol_pct: 1 }],
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.dx_over_eta > 2.1) out.push({ severity: o.dx_over_eta > 4 ? 'warn' : 'advise', title: 'The smallest eddies are not resolved', detail: `Cell size is ${o.dx_over_eta.toFixed(1)} Kolmogorov lengths on ${Math.round(i.n)}³ cells.`, action: `Refine to about ${Math.ceil((i.n * o.dx_over_eta) / 2.1)} cells per side for a DNS, lower the Reynolds number, or switch to the LES analysis, which models the unresolved dissipation.`, basis: 'Kolmogorov scale η = (ν³/ε)^¼ from the computed peak dissipation' });
    if (Math.abs(o.budget_ratio - 1) > 0.01) out.push({ severity: 'advise', title: 'Energy budget does not close', detail: `Energy lost is ${(100 * o.budget_ratio).toFixed(1)}% of the integrated dissipation.`, action: 'Reduce the CFL number: the imbalance is the time-integration error.', basis: 'Kinetic-energy equation' });
    out.push({ severity: 'info', title: 'Use this case to qualify any flow solver', detail: 'The Taylor–Green vortex is the standard test of whether a scheme conserves energy and how much numerical dissipation it adds.', action: 'Run the mesh-convergence study on the 2-D case to confirm second-order accuracy, then repeat the 3-D case on two grids; apply the same discipline to external solvers through the High-fidelity bridge.', basis: 'Code verification practice (AIAA G-077, ASME V&V 20)' });
    return out;
  },
};
const dnsMemo = memoLast((i, ctx) => dns3d.runNow(i, ctx));

// ---- LES: turbulent channel flow and Taylor–Green vortex ----------------------------------------
const LES_FLOWS = ['Turbulent channel flow', 'Taylor–Green vortex'], SGS_MODELS = ['WALE', 'Smagorinsky with van Driest damping', 'None (implicit, under-resolved DNS)'];
const reichardt = (yp) => Math.log(1 + 0.41 * yp) / 0.41 + 7.8 * (1 - Math.exp(-yp / 11) - (yp / 11) * Math.exp(-yp / 3));
const tanhFaces = (n, beta) => N.range(n + 1, (j) => (beta > 0 ? 1 + Math.tanh(beta * ((2 * j) / n - 1)) / Math.tanh(beta) : (2 * j) / n));
const cfDean = (Reb) => 0.073 * Reb ** -0.25;
const sgsCfg = (name) => (name === SGS_MODELS[0] ? { model: 'les', sgs: 'wale' } : name === SGS_MODELS[1] ? { model: 'les', sgs: 'smagorinsky' } : { model: 'dns' });
/**
 * Fully developed channel flow between two walls at y = 0 and 2δ, periodic in x and z, driven at constant mass flux.
 * Units: bulk velocity U_b = 1, half height δ = 1, so ν = 2/Re_b. Returns the solver, wall-shear history and averages.
 */
function runChannel(o) {
  const n = o.n, nu = 2 / o.Reb, ut0 = Math.sqrt(cfDean(o.Reb) / 2), ReT0 = ut0 / nu, Lx = o.Lx, Lz = o.Lz;
  let beta = 0;
  if (!o.wallModel) { const f = (b) => 0.5 * (1 - Math.tanh(b * (1 - 2 / n)) / Math.tanh(b)) * ReT0 - o.y1; beta = f(0.3) < 0 ? 0 : f(3.2) > 0 ? 3.2 : N.brent(f, 0.3, 3.2, 1e-6); }
  const U0 = (y) => ut0 * reichardt((Math.min(y, 2 - y) * ut0) / nu), amp = o.amp ?? 0.15;
  const s = createSolver({
    n: [o.nx || n, n, o.nz || n], L: [Lx, 2, Lz], yFaces: beta > 0 ? tanhFaces(n, beta) : undefined, nu, bc: { y: 'wall' }, forcing: { bulk: 1 }, time: 'rk3', cfl: o.cfl ?? 1, wallModel: !!o.wallModel, ...sgsCfg(o.sgs), Cs: o.Cs, Cw: o.Cw,
    // mean turbulent profile plus streaks, streamwise vortices and a sinuous wave: transition completes within about 20 δ/U_b
    init: (x, y, z) => { const e = 1 - (y - 1) ** 2, a = (TWO_PI * z) / Lz, b = (TWO_PI * x) / Lx; return [U0(y) + amp * e * Math.cos(a) * (1 + 0.5 * Math.sin(b)), 0.6 * amp * e * e * Math.sin(a) * (1 + 0.6 * Math.cos(b)), 0.6 * amp * e * Math.sin(b) * Math.sin(PI * y)]; },
    perturb: { amplitude: 0.1, seed: 7 }, stats: { start: o.tEnd * (1 - o.avg) },
  });
  const hist = { t: [], cf: [] }; let tauSum = 0, nAvg = 0, steps = 0, sgsSum = 0;
  while (s.time < o.tEnd && steps < 200000 && !s.diverged) {
    s.run(1); steps++;
    const m = s.monitor(), tau = 0.5 * (m.wallShear[2] + m.wallShear[3]);
    if (s.time >= o.tEnd * (1 - o.avg)) { tauSum += tau; sgsSum += m.sgsDissipation; nAvg++; }
    if (steps % 2 === 0) { hist.t.push(s.time); hist.cf.push(2 * tau); }
    if (o.progress && steps % 20 === 0) o.progress(Math.min(s.time / o.tEnd, 1), `t = ${s.time.toFixed(1)} δ/U_b`);
  }
  return { s, nu, beta, hist, steps, tauW: nAvg ? tauSum / nAvg : NaN, sgs: nAvg ? sgsSum / nAvg : 0, ReT0, Lx, Lz };
}

const les3d = {
  id: 'les3d', title: '3-D LES: turbulent channel flow and Taylor–Green vortex', fidelity: 'numerical',
  summary: 'Large-eddy simulation: the large turbulent eddies are computed in three dimensions and time, and only the eddies smaller than the grid are modelled (Smagorinsky or WALE). The default case is fully developed turbulent flow between two walls, giving the mean velocity profile against the law of the wall, the Reynolds stresses and the friction coefficient.',
  equations: ['Navier–Stokes equations', 'Continuity equation', 'Conservation of momentum equation', 'Turbulent kinetic energy transport equation', 'LES subgrid-scale models', 'WALE', 'Wall-modelled LES formulations', 'Wall-resolved or wall-modelled turbulence simulations'],
  inputs: [
    { key: 'flow', label: 'Flow', type: 'select', options: LES_FLOWS, default: LES_FLOWS[0], group: 'Flow' },
    { key: 'Re', label: 'Reynolds number', unit: '-', default: 5600, min: 1000, max: 2e5, group: 'Flow', help: 'Channel: bulk velocity × full height / ν (5600 gives Re_τ ≈ 180, the classic DNS case). Taylor–Green: U·L/ν (1600 is the standard benchmark)' },
    { key: 'sgs', label: 'Subgrid-scale model', type: 'select', options: SGS_MODELS, default: SGS_MODELS[0], group: 'Model' },
    { key: 'Cs', label: 'Smagorinsky constant', unit: '-', default: 0.1, min: 0.05, max: 0.25, group: 'Model', help: '0.1 for wall-bounded shear flow, 0.17 for decaying isotropic turbulence' },
    { key: 'Cw', label: 'WALE constant', unit: '-', default: 0.325, min: 0.2, max: 0.6, group: 'Model', help: '0.325 corresponds to Cs = 0.1; 0.5 to Cs = 0.17' },
    { key: 'wall_model', label: 'Equilibrium (log-law) wall model', type: 'bool', default: false, group: 'Model', help: 'Channel only. Needed above Re ≈ 10⁴, where the near-wall eddies cannot be resolved; the grid is then uniform in y' },
    { key: 't_end', label: 'Simulated time', unit: 'δ/U_b or L/U', default: 20, min: 2, max: 2000, group: 'Flow', help: 'Channel: transition from the perturbed initial field takes about 20; statistics need several hundred for 2% accuracy' },
    { key: 'avg_frac', label: 'Fraction of the run used for statistics', unit: '-', default: 0.5, min: 0.1, max: 0.9, group: 'Flow' },
    { key: 'Lx', label: 'Channel length', unit: 'δ', default: 3.1416, min: 1.5, max: 12.6, group: 'Domain', help: 'π (small box) to 2π or 4π. A longer box needs proportionally more cells' },
    { key: 'Lz', label: 'Channel width', unit: 'δ', default: 1.5708, min: 0.8, max: 6.3, group: 'Domain' },
    { key: 'y1_plus', label: 'Target y⁺ of the first cell centre', unit: '-', default: 2, min: 0.5, max: 10, group: 'Numerics', help: 'Sets the wall-normal grid stretching when no wall model is used' },
    { key: 'n', label: 'Cells per direction', unit: '', default: 12, min: 8, max: 128, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: () => ({}),
  run(i, ctx) { return lesMemo(i, ctx); },
  convergence: { param: 'n', label: 'Cells per direction', levels: [12, 16, 20, 24], metric: 'grid_metric' },
  verify() {
    // laminar Poiseuille flow through the LES code path: WALE must return zero eddy viscosity in pure shear
    const G = 1, nu = 0.1, lam = createSolver({ n: [4, 12, 4], L: [2, 2, 2], yFaces: tanhFaces(12, 1.2), nu, bc: { y: 'wall' }, forcing: { gradient: [G, 0, 0] }, model: 'les', sgs: 'wale', time: 'ab2', diffusion: 'implicit', dt: 0.05, init: (x, y) => [(G / (2 * nu)) * y * (2 - y), 0, 0] });
    lam.run(1500); const dl = lam.diagnostics(), m = lam.monitor();
    // uniform shear du/dy = 1 between a fixed and a moving wall: Smagorinsky νt = (Cs·Δ)²·|S| with |S| = 1
    const sh = createSolver({ n: [4, 8, 4], L: [1, 2, 1], nu: 1, bc: { y: ['wall', { type: 'wall', velocity: [2, 0, 0] }] }, model: 'les', sgs: 'smagorinsky', Cs: 0.17, vanDriest: false, time: 'ab2', dt: 1e-4, init: (x, y) => [y, 0, 0] });
    sh.step(); const dlt = Math.cbrt(0.25 * 0.25 * 0.25);
    const mf = createSolver({ n: [4, 16, 4], L: [2, 2, 2], nu: 0.1, bc: { y: 'wall' }, forcing: { bulk: 1 }, init: [1, 0, 0] }); mf.run(400); const dm = mf.diagnostics();
    const ut = wallFriction(17.5, 0.01, 1.5e-5), yp = (0.01 * ut) / 1.5e-5;
    return [
      N.check('Laminar channel: wall shear balances the pressure gradient', 0.5 * (m.wallShear[2] + m.wallShear[3]), G, 1e-6, 'Integral momentum balance τ_w = −δ·dp/dx (stretched grid, implicit diffusion)'),
      N.check('Laminar channel: bulk velocity', dl.bulk, G / (3 * nu), 0.03, 'Poiseuille solution U_b = δ²(−dp/dx)/(3ν); second-order grid error'),
      N.check('WALE eddy viscosity vanishes in pure shear', 1 + m.nutMax / nu, 1, 1e-9, 'Nicoud & Ducros (1999): the WALE operator is zero for laminar shear'),
      N.check('Smagorinsky eddy viscosity in uniform shear', sh.monitor().nutMax, (0.17 * dlt) ** 2, 1e-9, 'νt = (Cs·Δ)²·√(2 S:S) with du/dy = 1'),
      N.check('Constant-mass-flux forcing: pressure gradient 3νU_b/δ²', dm.gradient, 3 * 0.1, 0.03, 'Poiseuille flow at fixed bulk velocity (16 cells across)'),
      N.check('Wall model inverts the log law', ut * (Math.log(yp) / 0.41 + 5.2), 17.5, 1e-8, 'u⁺ = ln(y⁺)/0.41 + 5.2'),
      N.check('Discrete mass conservation with walls', 1 + dl.divMax, 1, 1e-10, 'Exact projection (cosine/Fourier transforms and tridiagonal solve)'),
    ];
  },
  validation: [{ name: 'Channel friction coefficient at Re_b = 5600 (Re_τ ≈ 180)', source: 'Kim, Moin & Moser, J. Fluid Mech. 177 (1987) DNS: Cf = 8.18×10⁻³; Dean (1978) correlation 0.073·Re_b^−0.25 = 8.44×10⁻³', inputs: { flow: LES_FLOWS[0], n: 32, t_end: 300, sgs: SGS_MODELS[0], wall_model: false }, sweep: { key: 'Re', values: [5600] }, target: 'Cf', observed: [8.18e-3], tol_pct: 15 }],
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (i.flow === LES_FLOWS[1]) { out.push({ severity: 'info', title: 'Share of dissipation carried by the model', detail: `The subgrid model supplies ${(100 * o.sgs_fraction).toFixed(0)}% of the dissipation at its peak.`, action: 'A well-resolved LES keeps this below about 20–30%; above that the result depends on the model constant. Refine, or compare WALE and Smagorinsky.', basis: 'Resolved versus modelled dissipation' }); return out; }
    if (o.eddy_turnovers < 5) out.push({ severity: 'warn', title: 'Statistics are not converged in time', detail: `Averaged over ${o.eddy_turnovers.toFixed(1)} eddy turnover times δ/u_τ; friction and stresses still wander by 5–10% at this length.`, action: 'Increase the simulated time until at least 10 turnovers are averaged (about 150 δ/U_b at Re_τ 180), then check that Cf no longer drifts.', basis: 'Sampling error of turbulence statistics' });
    if (!i.wall_model && (o.dx_plus > 60 || o.dz_plus > 25)) out.push({ severity: 'warn', title: 'Near-wall streaks are under-resolved', detail: `Δx⁺ = ${o.dx_plus.toFixed(0)}, Δz⁺ = ${o.dz_plus.toFixed(0)} (wall-resolved LES needs about 50 and 15–20).`, action: 'Add cells, shrink the box, or switch the wall model on and accept modelled wall shear.', basis: 'Wall-resolved LES resolution guidelines' });
    if (Math.abs(o.Cf_err_pct) > 10) out.push({ severity: 'advise', title: 'Friction differs from the reference', detail: `Cf = ${o.Cf.toExponential(3)} against ${o.Cf_ref.toExponential(3)} (${o.Cf_err_pct.toFixed(0)}%).`, action: 'On coarse grids second-order LES typically under-predicts friction by 10–20% (log-layer mismatch). Refine in the wall-parallel directions first and run the mesh-convergence study.', basis: 'Dean correlation / Kim, Moin & Moser DNS' });
    out.push({ severity: 'info', title: 'What this run can and cannot replace', detail: 'A periodic channel is the calibration case for LES. It shows that the solver sustains and resolves wall turbulence; it is not an aircraft.', action: 'For separated flow, buffet or airframe noise on the real geometry export a body-fitted wall-modelled LES case from the High-fidelity bridge; use the y⁺ and Δ⁺ figures here to size that mesh.', basis: 'Resolution requirements of scale-resolving simulation' });
    return out;
  },
};
const lesMemo = memoLast((i, ctx) => (i.flow === LES_FLOWS[1] ? lesTgv(i, ctx) : lesChannel(i, ctx)));
function lesChannel(i, ctx) {
  const n = cellsOf(i.n), Reb = i.Re, wm = !!i.wall_model, avg = N.clamp(i.avg_frac, 0.1, 0.9);
  const r = runChannel({ n, Reb, tEnd: i.t_end, avg, Lx: i.Lx, Lz: i.Lz, y1: i.y1_plus, wallModel: wm, sgs: i.sgs, Cs: i.Cs, Cw: i.Cw, progress: ctx?.progress });
  const s = r.s, nu = r.nu, S = s.statistics(), d = s.diagnostics(), ut = Math.sqrt(Math.max(r.tauW, 1e-300)), ReT = ut / nu, Cf = 2 * r.tauW, CfRef = cfDean(Reb), warnings = [];
  const half = n >> 1, g = s.grid, fold = (A, sgn = 1) => N.range(half, (j) => 0.5 * (A[j] + sgn * A[n - 1 - j]));
  const yP = N.range(half, (j) => g.y[j] * ReT), U = S ? fold(S.U) : N.range(half, () => NaN), uP = U.map((v) => v / ut), rms = (A) => (S ? fold(A).map((v) => Math.sqrt(Math.max(v, 0)) / ut) : U), uu = rms(S?.uu), vv = rms(S?.vv), ww = rms(S?.ww), uv = S ? fold(S.uv, -1).map((v) => -v / (ut * ut)) : U;
  const nut = S ? fold(S.nut).map((v) => v / nu) : U, visc = uP.map((v, j) => (j === 0 ? v / yP[0] : (v - uP[j - 1]) / (yP[j] - yP[j - 1]))), Ucl = S ? 0.5 * (S.U[half - 1] + S.U[n - half]) : NaN;
  const lg = N.range(half).filter((j) => yP[j] > 30 && g.y[j] < 0.5), logDev = lg.length ? N.mean(lg.map((j) => Math.abs(uP[j] - (Math.log(yP[j]) / 0.41 + 5.2)))) : NaN, kU = N.argmax(uu);
  const dxP = (r.Lx / g.nx) * ReT, dzP = (r.Lz / g.nz) * ReT, y1 = g.y[0] * ReT, turnovers = (i.t_end * avg * ut), sgsFrac = d.dissipation + r.sgs > 0 ? r.sgs / (d.dissipation + r.sgs) : 0, cfT = tailOf(r.hist.cf, avg), cfWander = N.std(cfT) / (N.mean(cfT) || 1);
  if (s.diverged) warnings.push('The time march diverged.');
  if (n <= 20) warnings.push(COARSE_NOTE(`${g.nx} × ${g.ny} × ${g.nz} cells in a ${r.Lx.toFixed(2)}δ × 2δ × ${r.Lz.toFixed(2)}δ box`, 'Cells per direction'));
  warnings.push(`Statistics are averaged over ${(i.t_end * avg).toFixed(0)} δ/U_b = ${turnovers.toFixed(1)} eddy turnover times δ/u_τ (${S ? S.samples : 0} samples); the wall shear fluctuates by ±${(100 * cfWander).toFixed(0)}% within that window. Published channel statistics use 10 or more turnovers: expect several percent of sampling error in every profile here.`);
  if (!wm) warnings.push(`Resolution in wall units: Δx⁺ = ${dxP.toFixed(0)}, Δz⁺ = ${dzP.toFixed(0)}, first cell centre at y⁺ = ${y1.toFixed(1)}, ${n} cells across the channel. Wall-resolved LES needs about Δx⁺ ≤ 50, Δz⁺ ≤ 20 and y⁺ ≈ 1${dxP > 60 || dzP > 25 || y1 > 3 ? ': this grid is coarser, so the near-wall streaks are only partly captured and the friction is typically 10–20% low' : ''}.`);
  else warnings.push(`Wall-modelled LES: the wall shear comes from the log law applied at the first cell centre (y⁺ = ${y1.toFixed(0)}), not from resolved near-wall eddies. ${y1 < 30 ? 'The first cell lies below y⁺ = 30, inside the buffer layer, where the log law does not hold: coarsen the wall-normal grid or switch the wall model off. ' : ''}Coarse wall-modelled grids show the known log-layer mismatch (velocity too high above the first cells).`);
  if (r.Lx * ReT < 300 || r.Lz * ReT < 100) warnings.push('The box is smaller than the minimal flow unit (about 300 × 100 wall units): turbulence may not sustain itself.');
  else if (r.Lx < 6 || r.Lz < 3) warnings.push(`The periodic box (${r.Lx.toFixed(2)}δ × ${r.Lz.toFixed(2)}δ) is smaller than the 2πδ × πδ used for reference data: the largest outer-layer structures are constrained, which mainly affects the profiles near the centreline.`);
  if (ReT > 1000 && !wm) warnings.push(`Re_τ ≈ ${ReT.toFixed(0)}: wall-resolved LES at this Reynolds number needs far more cells than this grid has. Switch the wall model on.`);
  if (i.sgs === SGS_MODELS[2]) warnings.push('No subgrid model: the dispersion error of the second-order scheme acts as an uncontrolled implicit model. This often gives good friction on coarse grids by cancellation of errors and is not a converged result.');
  if (Cf < 0.5 * CfRef) warnings.push('The friction is less than half the turbulent value: the flow has laminarised or has not yet transitioned. Run longer, enlarge the box or refine the grid.');
  const yl = N.logspace(Math.max(0.5, 0.5 * y1), Math.max(ReT, 20), 50), sx = s.slice('z', Math.floor(g.nz / 2), 'u', { max: 80 }), jS = N.argmin(yP.map((v) => Math.abs(v - 15))), sw = s.slice('y', jS, 'u', { max: 80 }), [th, ch] = thin(r.hist.t, r.hist.cf, 300);
  return {
    kpis: [
      { key: 'Re_tau', label: 'Friction Reynolds number u_τδ/ν', value: ReT, unit: '-', note: `Dean correlation gives ${r.ReT0.toFixed(0)}` },
      { key: 'Cf', label: 'Friction coefficient 2τ_w/U_b²', value: Cf, unit: '-', status: Math.abs(Cf / CfRef - 1) < 0.1 ? 'ok' : Math.abs(Cf / CfRef - 1) < 0.25 ? 'warn' : 'bad' },
      { key: 'Cf_ref', label: 'Friction coefficient, Dean correlation', value: CfRef, unit: '-', note: '0.073·Re_b^−0.25 (empirical; DNS at Re_b 5600 gives 8.18×10⁻³)' },
      { key: 'Cf_err_pct', label: 'Friction error against the correlation', value: 100 * (Cf / CfRef - 1), unit: '%' },
      { key: 'Ucl_over_Ub', label: 'Centreline / bulk velocity', value: Ucl, unit: '-', note: 'DNS at Re_τ 180: 1.16' },
      { key: 'Ucl_plus', label: 'Centreline velocity U⁺', value: Ucl / ut, unit: '-', note: 'DNS at Re_τ 180: 18.2' },
      { key: 'loglaw_dev', label: 'Mean deviation from the log law (y⁺ > 30, y < 0.5δ)', value: logDev, unit: 'u_τ', status: logDev < 1 ? 'ok' : 'warn' },
      { key: 'urms_peak_plus', label: 'Peak streamwise fluctuation u′⁺', value: uu[kU], unit: '-', note: `at y⁺ = ${yP[kU].toFixed(0)}; DNS at Re_τ 180: 2.66 at y⁺ ≈ 15` },
      { key: 'uv_peak_plus', label: 'Peak Reynolds shear stress −u′v′⁺', value: N.amax(uv), unit: '-', note: 'DNS at Re_τ 180: 0.72' },
      { key: 'y1_plus', label: 'First cell centre y⁺', value: y1, unit: '-', status: wm ? (y1 >= 30 ? 'ok' : 'warn') : y1 <= 3 ? 'ok' : 'warn' },
      { key: 'dx_plus', label: 'Streamwise cell size Δx⁺', value: dxP, unit: '-', status: wm || dxP <= 60 ? 'ok' : 'warn' },
      { key: 'dz_plus', label: 'Spanwise cell size Δz⁺', value: dzP, unit: '-', status: wm || dzP <= 25 ? 'ok' : 'warn' },
      { key: 'eddy_turnovers', label: 'Averaging time', value: turnovers, unit: 'δ/u_τ', status: turnovers >= 10 ? 'ok' : 'warn' },
      { key: 'sgs_fraction', label: 'Modelled share of dissipation', value: sgsFrac, unit: '-', note: 'Subgrid dissipation / (subgrid + resolved)' },
      { key: 'nut_max_ratio', label: 'Peak eddy viscosity ν_t/ν', value: d.nutMax / nu, unit: '-' },
      { key: 'div_max', label: 'Largest cell divergence', value: d.divMax, unit: 'U_b/δ', status: d.divMax < 1e-9 ? 'ok' : 'warn' },
      { key: 'steps', label: 'Time steps', value: r.steps, unit: '' },
      { key: 'cells', label: 'Grid cells', value: g.nx * g.ny * g.nz, unit: '' },
    ],
    plots: [
      { type: 'line', title: 'Mean velocity in wall units', xlabel: 'y⁺ [-]', ylabel: 'U⁺ [-]', xlog: true, series: [{ name: `LES ${g.nx}×${g.ny}×${g.nz}`, x: yP, y: uP }, { name: 'Sublayer U⁺ = y⁺', x: yl.filter((v) => v < 12), y: yl.filter((v) => v < 12), style: 'dash' }, { name: 'Log law ln(y⁺)/0.41 + 5.2', x: yl.filter((v) => v > 8), y: yl.filter((v) => v > 8).map((v) => Math.log(v) / 0.41 + 5.2), style: 'dash' }] },
      { type: 'line', title: 'Turbulence intensities and shear stress', xlabel: 'y⁺ [-]', ylabel: 'Fluctuation / u_τ, stress / u_τ² [-]', series: [{ name: 'u′⁺ streamwise', x: yP, y: uu }, { name: 'v′⁺ wall-normal', x: yP, y: vv }, { name: 'w′⁺ spanwise', x: yP, y: ww }, { name: '−u′v′⁺ resolved', x: yP, y: uv }, { name: 'Total stress 1 − y/δ', x: yP, y: yP.map((v) => 1 - v / ReT), style: 'dash' }] },
      { type: 'line', title: 'Friction coefficient history', xlabel: 'Time [δ/U_b]', ylabel: 'Cf [-]', series: [{ name: 'Instantaneous, both walls', x: th, y: ch }], annotations: [{ y: CfRef, label: 'Dean correlation' }, { x: i.t_end * (1 - avg), label: 'Averaging starts' }] },
      heatOf(`Instantaneous streamwise velocity, x–y plane at t = ${s.time.toFixed(0)}`, sx, 'u/U_b [-]', { xlabel: 'x/δ [-]', ylabel: 'y/δ [-]' }),
      heatOf(`Near-wall streaks: u in the plane y⁺ ≈ ${yP[jS].toFixed(0)}`, sw, 'u/U_b [-]', { xlabel: 'x/δ [-]', ylabel: 'z/δ [-]' }),
      { type: 'line', title: 'Stress balance: viscous, resolved and modelled', xlabel: 'y⁺ [-]', ylabel: 'Stress / u_τ² [-]', series: [{ name: 'Viscous dU⁺/dy⁺', x: yP, y: visc }, { name: 'Resolved −u′v′⁺', x: yP, y: uv }, { name: 'Subgrid ν_t/ν · dU⁺/dy⁺', x: yP, y: nut.map((v, j) => v * visc[j]) }, { name: 'Sum', x: yP, y: uv.map((v, j) => v + (1 + nut[j]) * visc[j]) }, { name: 'Exact total 1 − y/δ', x: yP, y: yP.map((v) => 1 - v / ReT), style: 'dash' }] },
    ],
    tables: [{ title: 'Mean profile (lower and upper halves averaged)', columns: ['y/δ', 'y⁺', 'U⁺', 'Log law', 'u′⁺', 'v′⁺', 'w′⁺', '−u′v′⁺', 'ν_t/ν'], rows: N.range(half, (j) => [+g.y[j].toFixed(4), +yP[j].toFixed(1), +uP[j].toFixed(2), +(Math.log(yP[j]) / 0.41 + 5.2).toFixed(2), +uu[j].toFixed(2), +vv[j].toFixed(2), +ww[j].toFixed(2), +uv[j].toFixed(3), +nut[j].toFixed(2)]) }],
    outputs: { beta_stretch: r.beta, grid_metric: Cf },
    warnings,
    models: ['Filtered incompressible Navier–Stokes equations (LES), constant mass flux', i.sgs === SGS_MODELS[0] ? `WALE subgrid model, Cw = ${i.Cw}` : i.sgs === SGS_MODELS[1] ? `Smagorinsky subgrid model, Cs = ${i.Cs}, van Driest wall damping (A⁺ = 25)` : 'No subgrid model', wm ? 'Equilibrium log-law wall model (κ = 0.41, B = 5.2) imposing the wall shear stress' : 'Wall-resolved: no-slip walls with tanh-stretched wall-normal grid', 'Staggered grid, second-order energy-conserving central differences, three-stage Runge–Kutta projection', 'Pressure by Fourier transforms in x and z and a tridiagonal solve in y'],
    assumptions: ['Fully developed, statistically steady flow; periodic in the streamwise and spanwise directions', 'Initial field is a turbulent mean profile with large streak and vortex perturbations; the first part of the run is discarded', 'Statistics are plane and time averages of cell-centred velocities, which slightly damps the fluctuation levels', 'Smooth walls, constant properties, incompressible'],
  };
}
function lesTgv(i, ctx) {
  const n = cellsOf(i.n), Re = i.Re, nu = 1 / Re, cfgS = sgsCfg(i.sgs), tEnd = Math.min(i.t_end, 40);
  const s = runTgv({ n, Re, tEnd, model: cfgS.model, sgs: cfgS.sgs, Cs: i.Cs > 0.12 ? i.Cs : 0.17, Cw: i.Cw > 0.4 ? i.Cw : 0.5, progress: ctx?.progress });
  const H = s.history, d = s.diagnostics(), B = energyBudget(H), tot = H.dissipation.map((v, k) => v + H.sgs[k]), kM = N.argmax(tot), kD = N.argmax(B.dE), sgsFrac = tot[kM] > 0 ? H.sgs[kM] / tot[kM] : 0, warnings = [];
  if (s.diverged || !Number.isFinite(d.ke)) warnings.push('The time march diverged.');
  if (n <= 20) warnings.push(COARSE_NOTE(`${n}³ cells`, 'Cells per direction'));
  if (tEnd < i.t_end) warnings.push('The Taylor–Green run is limited to 40 L/U.');
  if (tEnd < 10) warnings.push(`The run ends at t = ${tEnd} L/U; at Re 1600 the dissipation peak is near t ≈ 9 and the decay phase follows.`);
  if (Math.abs(Re - 1600) < 1) warnings.push(`Reference: pseudo-spectral DNS at Re = 1600 gives a peak dissipation of about 0.0127 U³/L near t ≈ 9. This LES gives ${N.amax(B.dE).toExponential(2)} at t = ${B.t[kD]?.toFixed(1)}. Coarse second-order LES typically peaks early and low because the small scales that do the dissipating are not on the grid.`);
  if (sgsFrac > 0.5) warnings.push(`The subgrid model carries ${(100 * sgsFrac).toFixed(0)}% of the dissipation at the peak: this is a very coarse LES and the result depends strongly on the model constant.`);
  const [tt, ke] = thin(H.t, H.ke, 300), [, er] = thin(H.t, H.dissipation, 300), [, es] = thin(H.t, H.sgs, 300), [tb, db] = thin(B.t, B.dE, 300), sv = s.slice('z', 0, 'vorticity', { max: 80 }), sn = s.slice('z', 0, 'nut', { max: 80 }), sq = s.slice('z', 0, 'q', { max: 80 });
  return {
    kpis: [
      { key: 'eps_total_peak', label: 'Peak total dissipation rate', value: tot[kM], unit: 'U³/L', note: `resolved + subgrid, at t = ${H.t[kM].toFixed(1)} L/U` },
      { key: 'eps_peak', label: 'Peak energy decay rate −dE/dt', value: N.amax(B.dE), unit: 'U³/L', note: Math.abs(Re - 1600) < 1 ? 'DNS reference ≈ 0.0127 near t ≈ 9' : '' },
      { key: 't_eps_peak', label: 'Time of peak decay rate', value: B.t[kD] ?? NaN, unit: 'L/U' },
      { key: 'sgs_fraction', label: 'Modelled share of dissipation at the peak', value: sgsFrac, unit: '-', status: sgsFrac < 0.3 ? 'ok' : 'warn' },
      { key: 'E_final', label: 'Kinetic energy at the end', value: d.ke, unit: 'U²', note: `${(100 * d.ke / 0.125).toFixed(0)}% of the initial energy` },
      { key: 'budget_ratio', label: 'Energy lost / integrated (resolved + subgrid) dissipation', value: B.ratio, unit: '-', status: Math.abs(B.ratio - 1) < 0.03 ? 'ok' : 'warn' },
      { key: 'nut_max_ratio', label: 'Peak eddy viscosity ν_t/ν', value: d.nutMax / nu, unit: '-' },
      { key: 'div_max', label: 'Largest cell divergence', value: d.divMax, unit: 'U/L', status: d.divMax < 1e-9 ? 'ok' : 'warn' },
      { key: 'steps', label: 'Time steps', value: d.step, unit: '' },
      { key: 'cells', label: 'Grid cells', value: n ** 3, unit: '' },
    ],
    plots: [
      { type: 'line', title: 'Dissipation: resolved, subgrid and total', xlabel: 'Time [L/U]', ylabel: 'Rate [U³/L]', series: [{ name: 'Resolved ν⟨|∇u|²⟩', x: tt, y: er }, { name: 'Subgrid ⟨2ν_t S:S⟩', x: tt, y: es }, { name: '−dE/dt', x: tb, y: db, style: 'dash' }] },
      { type: 'line', title: 'Kinetic-energy decay', xlabel: 'Time [L/U]', ylabel: 'Kinetic energy [U²]', series: [{ name: `LES ${n}³`, x: tt, y: ke }] },
      heatOf(`Vorticity magnitude, plane z = ${sv.position.toFixed(2)} L at t = ${s.time.toFixed(1)}`, sv, '|ω| [U/L]', { xlabel: 'x [L]', ylabel: 'y [L]' }),
      heatOf('Q-criterion in the same plane', sq, 'Q [U²/L²]', { xlabel: 'x [L]', ylabel: 'y [L]', diverging: true }),
      heatOf('Eddy-viscosity ratio in the same plane', sn, 'ν_t/ν [-]', { xlabel: 'x [L]', ylabel: 'y [L]' }),
    ],
    outputs: { grid_metric: tot[kM] },
    warnings,
    models: ['Filtered incompressible Navier–Stokes equations (LES) in a periodic box', cfgS.model === 'dns' ? 'No subgrid model' : cfgS.sgs === 'wale' ? 'WALE subgrid model' : 'Smagorinsky subgrid model without wall damping', 'Second-order energy-conserving central differences, three-stage Runge–Kutta projection, FFT pressure solver'],
    assumptions: ['Triply periodic box of side 2πL, unit velocity scale', 'The central scheme adds no numerical dissipation, so the energy budget closes with resolved plus subgrid dissipation', 'Model constants below the free-turbulence values (Cs 0.17, Cw 0.5) are raised to them for this wall-free flow'],
  };
}

// ---- RANS of a 3-D body by immersed boundary ----------------------------------------------------
const RANS_BODIES = ['Wing from the case', 'Fuselage (ellipsoid)', 'Sphere', 'Imported surface'], RANS_MODELS = ['Spalart–Allmaras', 'Mixing length', 'None (laminar / implicit)'];
const UNIT_M = { m: 1, mm: 1e-3, cm: 1e-2, in: 0.0254, ft: 0.3048 };
/** Signed distance to a closed polygon (negative inside); X, Y are the nodes with the first repeated at the end. */
export function polySdf(X, Y) {
  const n = X.length - 1;
  return (x, y) => {
    let d2 = 1e30, ins = false;
    for (let i = 0, j = n - 1; i < n; j = i++) {
      const ex = X[j] - X[i], ey = Y[j] - Y[i], wx = x - X[i], wy = y - Y[i], t = N.clamp((wx * ex + wy * ey) / (ex * ex + ey * ey || 1e-30), 0, 1), bx = wx - ex * t, by = wy - ey * t, q = bx * bx + by * by;
      if (q < d2) d2 = q;
      if (Y[i] > y !== Y[j] > y && x < (ex * (y - Y[i])) / ey + X[i]) ins = !ins;
    }
    return (ins ? -1 : 1) * Math.sqrt(d2);
  };
}
const sinhFaces = (n, H, beta) => N.range(n + 1, (j) => (H * Math.sinh(beta * ((2 * j) / n - 1))) / Math.sinh(beta));
/** UV-sphere triangle surface (used by the verification of the voxeliser). */
function sphereMesh(c, r, nu = 32, nv = 16) {
  const P = [], T = [];
  for (let j = 0; j <= nv; j++) for (let k = 0; k < nu; k++) { const th = (PI * j) / nv, ph = (TWO_PI * k) / nu; P.push(c[0] + r * Math.sin(th) * Math.cos(ph), c[1] + r * Math.cos(th), c[2] + r * Math.sin(th) * Math.sin(ph)); }
  for (let j = 0; j < nv; j++) for (let k = 0; k < nu; k++) { const a = j * nu + k, b = j * nu + ((k + 1) % nu), c2 = a + nu, d = b + nu; T.push(a, b, c2, b, d, c2); }
  return { positions: P, triangles: T };
}
/** Geometry, reference quantities and domain of the body to be flown. Lengths in metres; the solver runs with unit speed. */
function ransBody(i, ctx, n) {
  const al = N.rad(i.alpha_deg), ca = Math.cos(al), sa = Math.sin(al), notes = [];
  let kind = i.body, shape = ctx?.case?.shape;
  if (kind === RANS_BODIES[3]) {
    const P = shape?.positions, ok = P && P.length >= 9;
    if (!ok) { notes.push('No imported surface is attached to the case (send one from the geometry workbench); the parametric body was used instead.'); kind = i.S > 0 && i.b > 0 ? RANS_BODIES[0] : i.fus_L > 0 && i.fus_D > 0 ? RANS_BODIES[1] : RANS_BODIES[2]; }
  }
  if (kind === RANS_BODIES[0] && !(i.S > 0 && i.b > 0)) { notes.push('The case has no wing; the fuselage ellipsoid was used instead.'); kind = RANS_BODIES[1]; }
  if (kind === RANS_BODIES[1] && !(i.fus_L > 0 && i.fus_D > 0)) { notes.push('No fuselage dimensions; a sphere was used instead.'); kind = RANS_BODIES[2]; }
  if (kind === RANS_BODIES[0]) {
    const g = wingGeom(i), s = i.b / 2, af = parseNaca(i.airfoil, i.tc), nd = nacaNodes(af, 64), sec = polySdf(nd.X, nd.Y), ct = g.cr * g.taper, cbar = i.S / i.b;
    const chord = (z) => g.cr * (1 - ((1 - g.taper) * z) / s), xle = (z) => g.cr / 4 + z * g.tanL - chord(z) / 4;
    const sdf = (x, y, z) => {
      const zz = Math.abs(z), zc = Math.min(zz, s), c = chord(zc), xb = x * ca - y * sa, yb = x * sa + y * ca; let xi = (xb - xle(zc)) / c, et = yb / c;
      if (g.twist) { const tw = (g.twist * zc) / s, cw = Math.cos(tw), sw = Math.sin(tw), dq = xi - 0.25; xi = 0.25 + dq * cw - et * sw; et = dq * sw + et * cw; }
      const d2 = c * (xi < -0.4 || xi > 1.4 || Math.abs(et) > 0.5 ? Math.hypot(Math.max(-xi, xi - 1, 0), Math.max(Math.abs(et) - 0.2, 0)) + 0.05 : sec(xi, et)), dz = zz - s;
      return dz <= 0 ? Math.max(d2, dz) : d2 > 0 ? Math.hypot(d2, dz) : dz;
    };
    const xmin = Math.min(0, xle(s)), xmax = Math.max(g.cr, xle(s) + ct), H = 2.2 * g.cr, Lz = s + Math.max(cbar, 0.2 * s), x0 = xmin - g.cr, Lx = xmax + 2 * g.cr - x0;
    const zs = 0.3 * s, cs = chord(zs);
    return { kind, name: `${af.name} wing, half model with a symmetry plane`, body: { sdf }, half: true, Sref: i.S / 2, Lref: cbar, Lflow: g.cr, notes, n: [2 * n, n, n], L: [Lx, 2 * H, Lz], origin: [x0, -H, 0], yFaces: sinhFaces(n, H, 3.5), zSlice: zs, xTE: (xle(zs) + cs) * ca, xLE: xle(zs) * ca, cSlice: cs * ca, thick: af.t * cbar, len: cbar, g, af, s, cr: g.cr, ct, blockage: (af.t * i.S / 2 + i.S / 2 * Math.abs(sa)) / (2 * H * Lz), bc: { x: [{ type: 'inflow', velocity: [1, 0, 0] }, 'outflow'], y: 'slip', z: 'slip' } };
  }
  let body, bb, Sref, Lref, name, thick;
  if (kind === RANS_BODIES[3]) {
    const P = shape.positions, flat = typeof P[0] === 'number', nV = flat ? P.length / 3 : P.length, sc = UNIT_M[shape.units] || 1, Q = new Float64Array(3 * nV); bb = [1e30, -1e30, 1e30, -1e30, 1e30, -1e30];
    let cx = 0, cy = 0; for (let k = 0; k < nV; k++) { cx += (flat ? P[3 * k] : P[k][0]) * sc; cy += (flat ? P[3 * k + 1] : P[k][1]) * sc; } cx /= nV; cy /= nV;
    for (let k = 0; k < nV; k++) {
      const x = (flat ? P[3 * k] : P[k][0]) * sc - cx, y = (flat ? P[3 * k + 1] : P[k][1]) * sc - cy, z = (flat ? P[3 * k + 2] : P[k][2]) * sc;
      Q[3 * k] = cx + x * ca + y * sa; Q[3 * k + 1] = cy - x * sa + y * ca; Q[3 * k + 2] = z; // nose-up rotation about the spanwise axis through the centroid
      for (let a = 0; a < 3; a++) { bb[2 * a] = Math.min(bb[2 * a], Q[3 * k + a]); bb[2 * a + 1] = Math.max(bb[2 * a + 1], Q[3 * k + a]); }
    }
    const T = shape.triangles; body = { positions: Q, triangles: T && T.length ? (typeof T[0] === 'number' ? T : T.flat()) : null };
    Lref = bb[1] - bb[0]; Sref = i.S_ref > 0 ? i.S_ref : (bb[1] - bb[0]) * (bb[5] - bb[4]); thick = Math.min(bb[3] - bb[2], bb[5] - bb[4]);
    name = `Imported surface (${nV} vertices, ${body.triangles ? Math.floor(body.triangles.length / 3) : Math.floor(nV / 3)} triangles)`;
    if (!(i.S_ref > 0)) notes.push('No reference area was given for the imported surface: coefficients use the planform bounding box (length × span). Enter the true reference area to compare with other data.');
    notes.push('Imported surface axes are taken as x downstream, y up, z spanwise; the surface is rotated nose-up by the angle of attack about its centroid. A surface that is not watertight is voxelised by majority vote of three ray directions and may leak.');
  } else {
    const sph = kind === RANS_BODIES[2], D = sph ? (i.fus_D > 0 ? i.fus_D : 1) : i.fus_D, Lb = sph ? D : i.fus_L, a = Lb / 2, r = D / 2;
    body = { sdf: (x, y, z) => { const xb = x * ca - y * sa, yb = x * sa + y * ca; if (sph) return Math.hypot(xb, yb, z) - r; const k0 = Math.hypot(xb / a, yb / r, z / r), k1 = Math.hypot(xb / (a * a), yb / (r * r), z / (r * r)); return k1 > 0 ? (k0 * (k0 - 1)) / k1 : -r; } };
    const ex = Math.hypot(a * ca, r * sa), ey = Math.hypot(a * sa, r * ca); bb = [-ex, ex, -ey, ey, -r, r];
    Lref = Lb; Sref = (PI * D * D) / 4; thick = D; name = sph ? `Sphere, diameter ${D} m` : `Ellipsoid of revolution ${Lb} m × ${D} m`;
  }
  const lx = bb[1] - bb[0], ly = bb[3] - bb[2], lz = bb[5] - bb[4], lc = Math.max(ly, lz), cl = N.clamp(i.clearance, 0.5, 4) * lc, up = Math.max(0.6 * lx, 1.5 * lc), dn = Math.max(1.2 * lx, 3 * lc);
  const Ly = ly + 2 * cl, Lz = lz + 2 * cl, yc = 0.5 * (bb[2] + bb[3]), flatBody = ly < 0.4 * lz;
  return { kind, name, body, half: false, Sref, Lref, Lflow: lx, notes, n: [2 * n, n, n], L: [lx + up + dn, Ly, Lz], origin: [bb[0] - up, yc - Ly / 2, 0.5 * (bb[4] + bb[5]) - Lz / 2], yFaces: flatBody ? sinhFaces(n, Ly / 2, 2.5).map((v) => v + yc) : undefined, zSlice: 0.5 * (bb[4] + bb[5]), xTE: bb[1], xLE: bb[0], cSlice: lx, thick, len: lx, blockage: (ly * lz * (kind === RANS_BODIES[3] ? 0.5 : PI / 4)) / (Ly * Lz), bc: { x: [{ type: 'inflow', velocity: [1, 0, 0] }, 'outflow'], y: 'slip', z: 'slip' } };
}

const rans3d = {
  id: 'rans3d', title: '3-D flow-field visualisation: immersed-boundary RANS over a wing, fuselage or imported shape', fidelity: 'numerical',
  summary: 'A picture of the three-dimensional flow, not a force prediction. The body (the wing of the case, a fuselage ellipsoid, or a surface imported from the geometry workbench) is immersed in a Cartesian grid and the Reynolds-averaged Navier–Stokes equations with the Spalart–Allmaras model are marched to a steady state: pressure and velocity fields, the wake and the tip vortex. A Cartesian grid that a browser can afford cannot resolve a flight-Reynolds-number boundary layer, so its lift is typically far too low and its drag far too high; viscous wing forces come from the viscous-section wing analysis, which is shown beside it.',
  equations: ['Reynolds-averaged Navier–Stokes equations', 'Navier–Stokes equations', 'Continuity equation', 'Conservation of momentum equation', 'Spalart–Allmaras', 'Immersed-boundary formulations', 'Wall-resolved or wall-modelled turbulence simulations'],
  inputs: [
    { key: 'body', label: 'Body', type: 'select', options: RANS_BODIES, default: RANS_BODIES[0], group: 'Geometry', help: 'An imported surface is used automatically when one has been sent from the geometry workbench' },
    { key: 'S', label: 'Wing area', unit: 'm²', default: 16.2, min: 0, group: 'Geometry' },
    { key: 'b', label: 'Span', unit: 'm', default: 11, min: 0, group: 'Geometry' },
    { key: 'taper', label: 'Taper ratio', unit: '-', default: 0.6, min: 0.05, max: 1, group: 'Geometry' },
    { key: 'sweep_deg', label: 'Quarter-chord sweep', unit: 'deg', default: 0, min: -30, max: 60, group: 'Geometry' },
    { key: 'twist_deg', label: 'Tip twist (washout negative)', unit: 'deg', default: -2, min: -10, max: 5, group: 'Geometry' },
    { key: 'airfoil', label: 'NACA section', type: 'text', default: '2412', group: 'Geometry' },
    { key: 'tc', label: 'Thickness ratio override', unit: '-', default: 0, min: 0, max: 0.4, group: 'Geometry', help: '0 uses the thickness in the designation' },
    { key: 'fus_L', label: 'Fuselage length', unit: 'm', default: 8.3, min: 0, group: 'Geometry' },
    { key: 'fus_D', label: 'Fuselage or sphere diameter', unit: 'm', default: 1.2, min: 0, group: 'Geometry' },
    { key: 'S_ref', label: 'Reference area for an imported surface', unit: 'm²', default: 0, min: 0, group: 'Geometry', help: '0 uses the planform bounding box' },
    { key: 'alpha_deg', label: 'Angle of attack', unit: 'deg', default: 4, min: -10, max: 25, group: 'Flow' },
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 60, min: 0.5, group: 'Flow' },
    ...ATM,
    { key: 'turb', label: 'Turbulence model', type: 'select', options: RANS_MODELS, default: RANS_MODELS[0], group: 'Model' },
    { key: 'wall_model', label: 'Slip wall with log-law wall stress (experimental)', type: 'bool', default: false, group: 'Model', help: 'Off: no-slip immersed wall, robust but with a numerical boundary layer one cell thick. On: only the wall-normal velocity is removed and the log-law shear stress is applied; this needs near-cubic cells at the surface and can give wrong forces on strongly stretched or very coarse grids' },
    { key: 'quick', label: 'Upwind-biased (QUICK) fraction of the convection scheme', unit: '-', default: 1, min: 0, max: 1, group: 'Numerics', help: '1 = QUICK, 0 = central. Central differencing is unstable on coarse grids at high Reynolds number' },
    { key: 'clearance', label: 'Domain clearance around a body', unit: 'body widths', default: 1.25, min: 0.5, max: 4, group: 'Numerics', help: 'Fuselage, sphere and imported shapes; the wing domain is set from the chord and span' },
    { key: 't_flow', label: 'Simulated time', unit: 'body lengths', default: 4, min: 1, max: 100, group: 'Numerics', help: 'Distance flown, in root chords (wing) or body lengths. The starting vortex must leave the domain: 6–10 for a wing' },
    { key: 'max_steps', label: 'Maximum time steps', unit: '', default: 120, min: 20, max: 50000, step: 1, discrete: true, group: 'Numerics' },
    { key: 'CL_wingrans_ref', label: 'Viscous-section wing lift (reference)', unit: '-', default: 0, min: -2, max: 4, group: 'Reference', help: 'Filled from the viscous-section wing analysis when it has been run; 0 = not available' },
    { key: 'CD_wing_ref', label: 'Viscous-section wing drag, profile + induced (reference)', unit: '-', default: 0, min: 0, max: 1, group: 'Reference' },
    { key: 'n', label: 'Cells across the domain (y and z); twice as many along x', unit: '', default: 12, min: 8, max: 64, step: 1, discrete: true, group: 'Numerics', help: '16 gives 32 × 16 × 16; 64 gives 128 × 64 × 64' },
  ],
  defaults: (c, up) => {
    const hasShape = !!(c.shape && c.shape.positions && c.shape.positions.length >= 9), wingOk = c.wing.S_m2 > 0 && c.wing.b_m > 0, fusOk = c.fuselage.len_m > 0 && c.fuselage.dia_m > 0;
    return { body: hasShape ? RANS_BODIES[3] : wingOk ? RANS_BODIES[0] : fusOk ? RANS_BODIES[1] : RANS_BODIES[2], S: c.wing.S_m2, b: c.wing.b_m, taper: c.wing.taper, sweep_deg: c.wing.sweep_deg, twist_deg: c.wing.twist_deg, airfoil: c.wing.airfoil, tc: c.wing.tc, fus_L: c.fuselage.len_m || undefined, fus_D: c.fuselage.dia_m || undefined, alpha_deg: wingOk ? c.flight.alpha_deg + c.wing.incidence_deg : 0, V: c.flight.V_ms, alt_m: c.atm.alt_m, dISA: c.atm.dISA_K, CL_wingrans_ref: wingOk ? up.cfd?.CL_wingrans : undefined, CD_wing_ref: wingOk && up.cfd?.CD_wing_profile > 0 ? up.cfd.CD_wing_profile + up.cfd.CD_wing_induced : undefined };
  },
  run(i, ctx) {
    const a = isa(i.alt_m, i.dISA), n = cellsOf(i.n, 8, 64), B = ransBody(i, ctx, n), nuS = a.nu / i.V, Re = B.Lref / nuS, M = i.V / a.a, lam = i.turb === RANS_MODELS[2];
    const s = createSolver({ n: B.n, L: B.L, origin: B.origin, yFaces: B.yFaces, nu: nuS, bc: B.bc, body: B.body, init: [1, 0, 0], model: lam ? 'dns' : 'rans', rans: i.turb === RANS_MODELS[1] ? 'mixing-length' : 'sa', mixingLengthMax: 0.09 * B.thick, wallModel: !!i.wall_model && !lam, scheme: 'quick', blend: N.clamp(i.quick, 0, 1), time: 'rk3', cfl: 1, historyEvery: 2, residual: true });
    const tEnd = i.t_flow * B.Lflow, maxSteps = Math.round(N.clamp(i.max_steps, 20, 50000));
    while (s.time < tEnd && s.stepCount < maxSteps && !s.diverged) { s.run(1); if (s.stepCount % 10 === 0) ctx?.progress?.(Math.min(s.time / tEnd, s.stepCount / maxSteps), `${(s.time / B.Lflow).toFixed(1)} lengths flown`); }
    const H = s.history, d = s.diagnostics(), g = s.grid, q = 0.5 * B.Sref, tc = H.t.map((t) => t / B.Lflow), clH = H.fy.map((v) => v / q), cdH = H.fx.map((v) => v / q), tl = (A) => N.mean(tailOf(A, 0.25));
    const CL = tl(clH), CD = tl(cdH), drift = (A) => { const k = A.length, m1 = N.mean(A.slice(Math.floor(0.75 * k))), m0 = N.mean(A.slice(Math.floor(0.5 * k), Math.floor(0.75 * k))); return Math.abs(m1 - m0) / Math.max(Math.abs(m1), 1e-12); }, drCL = Math.abs(CL) > 0.02 ? drift(clH) : 0, drCD = drift(cdH);
    const kz = N.argmin(g.z.map((z) => Math.abs(z - B.zSlice))), full = { max: 4096 }, sp = s.slice('z', kz, 'p', full), sdS = s.slice('z', kz, 'sd', full), su = s.slice('z', kz, 'u', full), pInf = N.mean(sp.z.map((r) => r[0])), cpOf = (v) => 2 * (v - pInf);
    const ol = s.outline('z', kz), lim = { max: 80 }, dxc = g.dx, dyc = N.amin(g.dy), cellsLen = B.len / dxc, cellsThk = B.thick / Math.max(dyc, B.half ? 0 : g.dz), warnings = [...B.notes];
    // surface pressure along the body in the plotted plane: first fluid cell above and below
    const xs = [], cpU = [], cpL = [];
    for (let c = 0; c < g.nx; c++) { let jt = -1, jb = -1; for (let j = g.ny - 2; j >= 1; j--) if (sdS.z[j][c] <= 0) { jt = j + 1; break; } for (let j = 1; j < g.ny - 1; j++) if (sdS.z[j][c] <= 0) { jb = j - 1; break; } if (jt > 0 && jb >= 0) { xs.push((g.x[c] - B.xLE) / B.cSlice); cpU.push(cpOf(sp.z[jt][c])); cpL.push(cpOf(sp.z[jb][c])); } }
    // wake survey one reference length behind the trailing edge (or as far as the domain allows)
    const xW = Math.min(B.xTE + B.Lref, g.x[g.nx - 3]), iW = N.argmin(g.x.map((x) => Math.abs(x - xW))), wake = su.z.map((r) => r[iW]), deficit = 1 - N.amin(wake), sw = s.slice('x', iW, 'wx', lim);
    // reference values for orientation
    let CLref = NaN, CDref = NaN;
    if (B.kind === RANS_BODIES[0]) { const L = vlm(B.g, 12, 3, 1, (xf) => B.af.camber(xf)[1], true), al = N.rad(i.alpha_deg), strip = L.a.strip.map((v, j) => v * al + L.z.strip[j]); CLref = L.a.CL * al + L.z.CL; CDref = L.trefftz(strip) + 2 * cfTurb(Math.max(Re, 1e4), M) * (1 + 2 * B.af.t + 60 * B.af.t ** 4); }
    else if (B.kind === RANS_BODIES[2]) CDref = Re > 3.5e5 ? 0.2 : 0.47;
    const hasWR = B.kind === RANS_BODIES[0] && i.CD_wing_ref > 0;
    warnings.unshift(`Flow-field visualisation only. The integrated forces of this immersed-boundary grid are not force-accurate (lift is typically a third to a half of the correct value and drag several times too high) and are not published to other suites. ${B.half ? (hasWR ? `For the wing use the viscous-section analysis: CL = ${i.CL_wingrans_ref.toFixed(3)}, wing CD = ${i.CD_wing_ref.toFixed(4)}.` : 'For wing forces run the viscous-section wing analysis (vortex lattice coupled to body-fitted section RANS) and the vortex-lattice analysis.') : 'For forces on this body use the drag build-up or a body-fitted solver through the High-fidelity bridge.'}`);
    if (s.diverged || !Number.isFinite(CL)) warnings.push('The time march diverged: raise the upwind fraction, switch the wall model on or refine the grid.');
    if (n <= 20) warnings.push(COARSE_NOTE(`${g.nx} × ${g.ny} × ${g.nz} cells`, 'Cells across the domain'));
    warnings.push(`Resolution: ${cellsLen.toFixed(1)} cells along the ${B.half ? 'mean chord' : 'body length'}${B.half ? ` (${(B.cr / dxc).toFixed(1)} at the root, ${(B.ct / dxc).toFixed(1)} at the tip)` : ''} and ${cellsThk.toFixed(1)} across the ${B.half ? 'maximum thickness' : 'smallest body dimension'}. ${cellsLen < 32 || cellsThk < 6 ? 'The immersed boundary places the wall to within about one cell, so at this resolution the leading-edge suction, the trailing edge and the boundary layer are smeared: the numerical boundary layer is about one cell thick, the flow separates early, the lift can be less than half the true value and the drag several times too high. Read this run as a flow visualisation and a trend, not as a drag prediction.' : 'This is enough for pressure lift within roughly 10%; friction drag still relies entirely on the wall model.'}`);
    warnings.push(`Wall layer: the first fluid cells sit at y⁺ ≈ ${d.yPlus.mean.toFixed(0)} on average (maximum ${d.yPlus.max.toFixed(0)}). ${lam ? 'No turbulence model is active.' : d.yPlus.mean < 5 ? 'The viscous sublayer is resolved.' : d.yPlus.mean <= 300 ? 'The viscous and buffer layers are not resolved; the wall shear comes from the log-law wall model, which is valid for attached flow at this y⁺ and not for separation or transition.' : 'This is far outside the log layer (y⁺ ≲ 300): the wall model is extrapolated and the friction drag, separation location and maximum lift are not predictive. A body-fitted mesh with y⁺ ≈ 1 or a wall-function mesh is needed for those: see the wall-layer analysis and the High-fidelity bridge.'}`);
    if (Math.max(drCL, drCD) > 0.02 || s.time < 0.98 * tEnd) warnings.push(`Not converged to a steady state: the mean ${drCL > drCD ? 'lift' : 'drag'} still changes by ${(100 * Math.max(drCL, drCD)).toFixed(1)}% between the last two quarters of the run${s.time < 0.98 * tEnd ? ` and the step limit stopped the run at ${(s.time / B.Lflow).toFixed(1)} of ${i.t_flow} lengths` : ''}. Increase the simulated time and the maximum number of steps.`);
    if (B.blockage > 0.05) warnings.push(`The body blocks ${(100 * B.blockage).toFixed(1)}% of the domain cross-section between slip walls; forces are raised by roughly twice that fraction. Increase the clearance.`);
    else if (B.half) warnings.push('The slip walls above and below the wing are 2.2 root chords away; like wind-tunnel walls they raise the lift slightly (a few percent).');
    if (M > 0.3) warnings.push(`Flight Mach number ${M.toFixed(2)}: this solver is incompressible. Compressibility raises the lift slope by about 1/√(1 − M²) = ${(1 / Math.sqrt(Math.max(1 - M * M, 0.05))).toFixed(2)}${M > 0.7 ? ' and shocks are absent altogether; transonic flow needs a compressible RANS solver (High-fidelity bridge)' : ''}.`);
    if (B.half) warnings.push('Wing only: no fuselage, nacelles, tail or dihedral; the section is the NACA shape at every station.');
    if (lam) warnings.push('No turbulence model: at this Reynolds number the result is an under-resolved, numerically damped solution, not a laminar flow.');
    const cpClip = (v) => N.clamp(cpOf(v), -2.5, 1.2), ovl = { overlay: ol, xlabel: 'x [m]', ylabel: 'y [m]' }, [tcT, clT] = thin(tc, clH, 300), [, cdT] = thin(tc, cdH, 300), [, rsT] = thin(tc, H.residual.map((v) => Math.max(v, 1e-16)), 300);
    const plots = [
      { type: 'line', title: 'Force convergence (iteration monitor; not force-accurate at this grid)', xlabel: `Distance flown [${B.half ? 'root chords' : 'body lengths'}]`, ylabel: 'Coefficient [-]', series: [{ name: 'CL', x: tcT.slice(1), y: clT.slice(1) }, { name: 'CD', x: tcT.slice(1), y: cdT.slice(1) }], annotations: [...(Number.isFinite(CLref) ? [{ y: CLref, label: 'CL, vortex lattice' }] : []), ...(hasWR ? [{ y: i.CL_wingrans_ref, label: 'CL, viscous-section wing' }] : [])] },
      { type: 'line', title: 'Residual history', xlabel: `Distance flown [${B.half ? 'root chords' : 'body lengths'}]`, ylabel: 'RMS |∂u/∂t| · L/V² [-]', ylog: true, series: [{ name: 'Momentum residual', x: tcT.slice(1), y: rsT.slice(1).map((v) => v * B.Lflow) }] },
      heatOf(`Pressure coefficient, plane z = ${g.z[kz].toFixed(2)} m`, s.slice('z', kz, 'p', lim), 'Cp [-]', { ...ovl, map: cpClip, diverging: true }),
      heatOf(`Velocity magnitude, plane z = ${g.z[kz].toFixed(2)} m`, s.slice('z', kz, 'speed', lim), '|V|/V∞ [-]', ovl),
      ...(lam ? [] : [heatOf('Eddy-viscosity ratio in the same plane', s.slice('z', kz, 'nut', lim), 'ν_t/ν [-]', ovl)]),
      { type: 'line', title: `Surface pressure in the plane z = ${g.z[kz].toFixed(2)} m (first fluid cells)`, xlabel: B.half ? 'x/c [-]' : 'x/L [-]', ylabel: '−Cp [-]', series: [{ name: 'Upper side', x: xs, y: cpU.map((v) => -v) }, { name: 'Lower side', x: xs, y: cpL.map((v) => -v) }] },
      { type: 'line', title: `Wake velocity profile at x = ${g.x[iW].toFixed(2)} m`, xlabel: 'u/V∞ [-]', ylabel: 'y [m]', series: [{ name: 'Streamwise velocity', x: wake, y: su.y }] },
      heatOf(`Streamwise vorticity in the cross-plane x = ${g.x[iW].toFixed(2)} m${B.half ? ' (tip vortex)' : ''}`, sw, 'ω_x · m·s⁻¹/V∞ [1/m]', { xlabel: 'z [m]', ylabel: 'y [m]', diverging: true, overlay: s.outline('x', iW) }),
    ];
    return {
      kpis: [
        { key: 'cells_streamwise', label: B.half ? 'Cells along the mean chord' : 'Cells along the body', value: cellsLen, unit: '', status: cellsLen >= 32 ? 'ok' : cellsLen >= 12 ? 'warn' : 'bad' },
        { key: 'cells_thickness', label: 'Cells across the thickness', value: cellsThk, unit: '', status: cellsThk >= 6 ? 'ok' : cellsThk >= 3 ? 'warn' : 'bad' },
        { key: 'wake_deficit', label: 'Peak wake velocity deficit', value: deficit, unit: 'V∞', note: `one reference length behind the body` },
        ...(hasWR ? [{ key: 'CL_wingrans_ref', label: 'Lift coefficient, viscous-section wing analysis (use this)', value: i.CL_wingrans_ref, unit: '-', status: 'ok', note: 'vortex lattice coupled to body-fitted section RANS, at the flight point' }, { key: 'CD_wing_ref', label: 'Wing drag coefficient, viscous-section wing analysis (use this)', value: i.CD_wing_ref, unit: '-', status: 'ok', note: 'profile + induced' }] : []),
        ...(Number.isFinite(CLref) ? [{ key: 'CL_vlm_ref', label: 'Lift coefficient, vortex lattice (reference)', value: CLref, unit: '-', note: 'Inviscid, incompressible, same planform, twist and camber' }] : []),
        { key: 'CL_rans3d', label: 'Integrated lift coefficient of this grid — not force-accurate', value: CL, unit: '-', status: 'warn', note: `not force-accurate at this grid: a ${cellsLen.toFixed(0)}-cell chord cannot carry the leading-edge suction; expect far too little lift. Mean of the last quarter of the run` },
        { key: 'CD_rans3d', label: 'Integrated drag coefficient of this grid — not force-accurate', value: CD, unit: '-', status: 'warn', note: 'not force-accurate at this grid: dominated by numerical diffusion and a one-cell boundary layer; expect several times the true drag' },
        ...(Number.isFinite(CLref) ? [{ key: 'CL_ratio', label: 'Immersed-boundary lift / vortex-lattice lift', value: CL / CLref, unit: '-', status: 'warn', note: 'a measure of how much lift the grid loses, not a viscous effect' }] : []),
        ...(Number.isFinite(CDref) ? [{ key: 'CD_ref_est', label: B.half ? 'Drag coefficient, handbook estimate (reference)' : 'Sphere drag coefficient, textbook (reference)', value: CDref, unit: '-', note: B.half ? 'Induced (Trefftz) plus flat-plate friction × form factor (empirical)' : Re > 3.5e5 ? 'supercritical, ≈ 0.2 (empirical)' : 'subcritical, ≈ 0.47 (empirical)' }] : []),
        { key: 'force_drift_pct', label: 'Change of mean force over the last quarter', value: 100 * Math.max(drCL, drCD), unit: '%', status: Math.max(drCL, drCD) < 0.02 ? 'ok' : 'warn' },
        { key: 'residual', label: 'Final momentum residual', value: d.residual * B.Lflow, unit: 'V²/L' },
        { key: 'y_plus_mean', label: 'Mean y⁺ of the first fluid cells', value: d.yPlus.mean, unit: '-', status: d.yPlus.mean < 5 || (d.yPlus.mean >= 30 && d.yPlus.mean <= 300) ? 'ok' : 'warn', note: d.yPlus.mean < 5 ? 'wall-resolved' : d.yPlus.mean <= 300 ? 'log layer: wall-modelled' : 'beyond the log layer: wall model extrapolated' },
        { key: 'y_plus_max', label: 'Maximum y⁺ of the first fluid cells', value: d.yPlus.max, unit: '-' },
        { key: 'nut_max_ratio', label: 'Peak eddy viscosity ν_t/ν', value: d.nutMax / nuS, unit: '-' },
        { key: 'Re_ref', label: 'Reynolds number on the reference length', value: Re, unit: '-' },
        { key: 'blockage_pct', label: 'Domain blockage', value: 100 * B.blockage, unit: '%', status: B.blockage < 0.05 ? 'ok' : 'warn' },
        { key: 'div_max', label: 'Largest cell divergence', value: d.divMax * B.Lflow, unit: 'V/L', status: d.divMax * B.Lflow < 1e-8 ? 'ok' : 'warn' },
        { key: 'steps', label: 'Time steps', value: d.step, unit: '' },
        { key: 'cells', label: 'Grid cells', value: g.nx * g.ny * g.nz, unit: '', note: `${s.solidCells} inside the body` },
      ],
      plots,
      tables: [{ title: 'Grid and body', columns: ['Quantity', 'Value'], rows: [['Body', B.name], ['Grid', `${g.nx} × ${g.ny} × ${g.nz}`], ['Domain [m]', `${B.L[0].toFixed(2)} × ${B.L[1].toFixed(2)} × ${B.L[2].toFixed(2)}`], ['Cell size Δx, Δy (min), Δz [m]', `${dxc.toPrecision(3)}, ${dyc.toPrecision(3)}, ${g.dz.toPrecision(3)}`], ['Reference area [m²]', +B.Sref.toPrecision(5)], ['Reference length [m]', +B.Lref.toPrecision(5)], ['Boundaries', B.half ? 'inflow, convective outflow, slip top and bottom, symmetry plane at the root, slip beyond the tip' : 'inflow, convective outflow, slip side walls'], ['Distance flown', `${(s.time / B.Lflow).toFixed(2)} lengths in ${d.step} steps`]] }],
      outputs: { CL_history_last: clH[clH.length - 1], CD_history_last: cdH[cdH.length - 1] },
      warnings,
      models: [lam ? 'Incompressible Navier–Stokes without a turbulence model' : i.turb === RANS_MODELS[1] ? 'Reynolds-averaged Navier–Stokes with a van Driest mixing-length eddy viscosity' : 'Reynolds-averaged Navier–Stokes with the Spalart–Allmaras one-equation model (standard constants, no trip term, fully turbulent)', 'Direct-forcing immersed boundary on a Cartesian staggered grid; forces from the momentum removed by the forcing', i.wall_model && !lam ? 'Equilibrium log-law wall model in the first fluid cells (κ = 0.41, B = 5.2)' : 'No wall model', `Convection: ${(100 * N.clamp(i.quick, 0, 1)).toFixed(0)}% QUICK, remainder central; three-stage Runge–Kutta projection marched to steady state`, 'Pressure by cosine transforms in x and z and a tridiagonal solve in y', B.kind === RANS_BODIES[3] ? 'Triangle surface voxelised by ray casting along three axes, signed distance by fast sweeping' : 'Analytical signed-distance description of the body'],
      assumptions: ['Incompressible, constant-property air; steady freestream along x', 'The wall is represented to within about one cell (first-order immersed boundary); there is no body-fitted boundary-layer mesh', 'Uniform inflow, convective outflow, slip (symmetry) side boundaries', 'Wall distance for the turbulence model is the distance to the immersed surface', B.half ? 'Half model with a symmetry plane at the wing root; coefficients use half the wing area' : 'Whole body in the domain'],
    };
  },
  convergence: { param: 'n', label: 'Cells across the domain', levels: [12, 16, 20, 24], metric: 'CL_rans3d' },
  verify() {
    // voxeliser: triangle sphere against the analytical signed distance
    const c = [0.5, 0.5, 0.5], r = 0.3, nv = 20, vs = createSolver({ n: [nv, nv, nv], L: [1, 1, 1], nu: 1, body: sphereMesh(c, r, 48, 24), dt: 1e-4 }), gv = vs.grid, sdv = vs.fields.sd; let eS = 0, cS = 0;
    for (let k = 0; k < nv; k++) for (let j = 0; j < nv; j++) for (let q = 0; q < nv; q++) { const ex = Math.hypot(gv.x[q] - c[0], gv.y[j] - c[1], gv.z[k] - c[2]) - r; if (Math.abs(ex) < 2 / nv) { eS += Math.abs(sdv[vs.index(q, j, k)] - ex); cS++; } }
    // free-stream preservation through inflow, outflow and slip boundaries
    const fs = createSolver({ n: [8, 8, 8], L: [2, 1, 1], nu: 1e-3, bc: { x: [{ type: 'inflow', velocity: [1, 0, 0] }, 'outflow'], y: 'slip', z: 'slip' }, init: [1, 0, 0], model: 'rans', scheme: 'quick' }); fs.run(5); let eF = 0; for (let k = 0; k < 8; k++) for (let j = 0; j < 8; j++) for (let q = 0; q < 8; q++) eF = Math.max(eF, Math.abs(fs.fields.u[fs.index(q, j, k)] - 1));
    // Spalart–Allmaras in a channel against the independent 1-D solver of the wall-layer analysis (wall units: u_τ = δ = 1)
    const ReT = 400, ny = 32, w1 = wallProfile(ReT, 'Spalart–Allmaras', 160), yF = tanhFaces(ny, 2.3), U1 = (y) => N.interp1(w1.y, w1.u, Math.min(y, 2 - y) * ReT);
    const ch = createSolver({ n: [4, ny, 4], L: [16, 2, 4], yFaces: yF, nu: 1 / ReT, bc: { y: 'wall' }, forcing: { gradient: [1, 0, 0] }, model: 'rans', rans: 'sa', time: 'ab2', diffusion: 'implicit', dt: 0.02, init: (x, y) => [0.9 * reichardt(Math.min(y, 2 - y) * ReT), 0, 0], nuTildeInit: (x, y) => { const dd = Math.min(y, 2 - y); return 0.41 * dd * (1 - 0.5 * dd) * 0.8; } });
    ch.run(1500); const dc = ch.diagnostics(), mc = ch.monitor();
    const sec = polySdf(...Object.values(nacaNodes(parseNaca('0012'), 64))); let area = 0; for (let a = 0; a < 400; a++) for (let b = 0; b < 80; b++) if (sec((a + 0.5) / 400, -0.1 + (b + 0.5) * 0.0025) < 0) area += 0.0025 / 400;
    return [
      N.check('Voxelised triangle sphere: solid volume', vs.solidCells / nv ** 3, (4 / 3) * PI * r ** 3, 0.04, 'Ray-cast inside test on a 20³ grid against the sphere volume (cell-count sampling error)'),
      N.check('Voxelised triangle sphere: signed-distance error near the surface', 1 + (eS / cS) * nv, 1, 0.15, 'Mean |d − d_exact| within two cells of the surface, in cell sizes'),
      N.check('Free stream is preserved', 1 + eF, 1, 1e-11, 'Uniform flow is an exact solution with inflow, convective outflow and slip walls'),
      N.check('Spalart–Allmaras channel: bulk velocity U_b⁺ against the 1-D solver', dc.bulk, w1.Ub, 0.03, 'Same model solved independently in one dimension (wall-layer analysis), Re_τ = 400'),
      N.check('Spalart–Allmaras channel: wall shear balances the pressure gradient', 0.5 * (mc.wallShear[2] + mc.wallShear[3]), 1, 0.03, 'Integral momentum balance; the residual is the remaining unsteadiness after 30 δ/u_τ'),
      N.check('NACA 0012 section area from the signed-distance polygon', area, 0.68085 * 0.12, 0.01, 'Integral of the closed-trailing-edge NACA thickness form: 0.6809·t·c²'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.cells_streamwise < 32 || o.cells_thickness < 6) out.push({ severity: 'warn', title: 'Grid is too coarse for quantitative forces', detail: `${o.cells_streamwise.toFixed(0)} cells along the body and ${o.cells_thickness.toFixed(1)} across its thickness.`, action: `Raise the cell count (the drag needs at least 32 cells per chord and 6 across the thickness; ${Math.ceil((i.n * 32) / Math.max(o.cells_streamwise, 1))} or more across the domain here) and run the mesh-convergence study. Until then use the panel, lattice and drag build-up analyses for numbers and this run for the flow picture.`, basis: 'Immersed-boundary resolution requirement' });
    if (o.y_plus_mean > 300) out.push({ severity: 'warn', title: 'Wall layer is not resolved or validly modelled', detail: `First fluid cells at y⁺ ≈ ${o.y_plus_mean.toFixed(0)}.`, action: 'Friction drag, separation and stall need a body-fitted mesh: take the first-cell height from the wall-layer analysis and export an SU2 or OpenFOAM case from the High-fidelity bridge.', basis: 'Law of the wall: the log layer ends near y⁺ ≈ 300 (0.15δ)' });
    if (o.force_drift_pct > 2) out.push({ severity: 'advise', title: 'Run longer', detail: `Mean force still drifting by ${o.force_drift_pct.toFixed(1)}%.`, action: 'Increase the simulated time and the step limit until the force history is flat.', basis: 'Iterative convergence' });
    if (o.CL_ratio !== undefined && Math.abs(o.CL_ratio - 1) > 0.2) out.push({ severity: 'advise', title: 'Lift differs from the vortex-lattice value', detail: `CL = ${o.CL_rans3d.toFixed(3)} against ${o.CL_vlm_ref.toFixed(3)} (ratio ${o.CL_ratio.toFixed(2)}).`, action: 'On a coarse immersed-boundary grid the difference is discretisation error, not physics. Refine before interpreting it as viscous or separation loss.', basis: 'Cross-check of two independent methods' });
    out.push({ severity: 'info', title: 'Where this solver stops', detail: 'Cartesian immersed-boundary RANS gives the three-dimensional pressure field, wake and vortex system of any closed shape with no meshing step.', action: 'For certification-grade drag, maximum lift, transonic flow or high-lift devices, use a body-fitted RANS mesh: the High-fidelity bridge writes ready-to-run SU2 and OpenFOAM cases for the same geometry and flight condition.', basis: 'Fidelity ladder' });
    return out;
  },
};

// ---- unsteady bluff-body flow: sphere and circular cylinder ------------------------------------
const BLUFF_BODIES = ['Sphere', 'Circular cylinder (spanwise periodic)'], BLUFF_MODELS = ['DNS (no model)', 'LES (WALE)'];
/** Schiller–Naumann sphere drag, valid for Re < 800; reduces to Stokes' 24/Re. */
export const cdSchillerNaumann = (Re) => (24 / Re) * (1 + 0.15 * Re ** 0.687);
/** Roshko's Strouhal number of a circular cylinder (laminar shedding 50 < Re < 150, irregular range 300 < Re < 2000). */
export const stRoshko = (Re) => (Re < 47 ? 0 : Re <= 200 ? 0.212 * (1 - 21.2 / Re) : 0.212 * (1 - 12.7 / Re));
/** Dominant frequency of a signal from its mean-crossings; returns NaN with fewer than three full cycles. */
export function crossingFrequency(t, y) {
  const m = N.mean(y), tc = [];
  for (let k = 1; k < y.length; k++) if (y[k - 1] - m < 0 && y[k] - m >= 0) tc.push(t[k - 1] + ((m - y[k - 1]) / (y[k] - y[k - 1])) * (t[k] - t[k - 1]));
  return tc.length >= 4 ? (tc.length - 1) / (tc[tc.length - 1] - tc[0]) : NaN;
}
/** Flow past a sphere or a spanwise-periodic cylinder of unit diameter in a uniform stream of unit speed. */
function runBluff(o) {
  const n = o.n, cyl = o.cyl, Hd = o.height, nu = 1 / o.Re, les = o.les;
  const dims = cyl ? [4 * n, n, 4] : [2 * n, n, n], L = cyl ? [2 * Hd, Hd, (4 * 2 * Hd) / dims[0]] : [2 * Hd, Hd, Hd], xc = cyl ? 0.25 * L[0] : 0.3 * L[0];
  const s = createSolver({
    n: dims, L, origin: [-xc, -Hd / 2, -L[2] / 2], yFaces: cyl ? sinhFaces(n, Hd / 2, 3) : undefined, nu, bc: { x: [{ type: 'inflow', velocity: [1, 0, 0] }, 'outflow'], y: 'slip', z: cyl ? 'periodic' : 'slip' },
    body: { sdf: cyl ? (x, y) => Math.hypot(x, y) - 0.5 : (x, y, z) => Math.hypot(x, y, z) - 0.5 }, init: [1, 0, 0], perturb: cyl ? { amplitude: 0.05, seed: 11, shape: (x, y) => Math.exp(-((x - 1) ** 2 + y * y)) } : undefined,
    model: les ? 'les' : 'dns', sgs: 'wale', Cw: 0.5, time: 'rk3', cfl: o.cfl ?? 1, historyEvery: 1, residual: true,
  });
  let steps = 0;
  while (s.time < o.tEnd && steps < 400000 && !s.diverged) { s.run(1); steps++; if (o.progress && steps % 10 === 0) o.progress(Math.min(s.time / o.tEnd, 1), `t = ${s.time.toFixed(1)} D/U`); }
  return { s, L, dims, area: cyl ? L[2] : PI / 4, blockage: cyl ? 1 / Hd : PI / 4 / (Hd * Hd), steps };
}

const bluff3d = {
  id: 'bluff3d', title: '3-D unsteady flow past a sphere or cylinder (DNS / LES validation)', fidelity: 'numerical',
  summary: 'Time-accurate Navier–Stokes solution of the flow around a sphere or a circular cylinder with the immersed-boundary method. The drag coefficient is compared with the Schiller–Naumann sphere law and the vortex-shedding frequency of the cylinder with the measured Strouhal number: this is the validation case for the immersed-boundary force prediction.',
  equations: ['Navier–Stokes equations', 'Continuity equation', 'Conservation of momentum equation', 'Vorticity transport equation', 'Immersed-boundary formulations', 'DNS', 'LES subgrid-scale models', 'WALE'],
  inputs: [
    { key: 'body', label: 'Body', type: 'select', options: BLUFF_BODIES, default: BLUFF_BODIES[0], group: 'Flow' },
    { key: 'Re', label: 'Reynolds number U·D/ν', unit: '-', default: 40, min: 2, max: 5000, group: 'Flow', help: 'Sphere: steady axisymmetric wake below about 210, unsteady above 270. Cylinder: vortex shedding above 47; 100–200 is the laminar shedding range' },
    { key: 'D_m', label: 'Diameter', unit: 'm', default: 0.05, min: 1e-5, group: 'Flow', help: 'Gives the equivalent speed and drag force in air at the case altitude' },
    ...ATM,
    { key: 'model', label: 'Turbulence treatment', type: 'select', options: BLUFF_MODELS, default: BLUFF_MODELS[0], group: 'Model', help: 'Use LES above Re ≈ 500, where the wake is turbulent and the grid cannot resolve it' },
    { key: 't_end', label: 'Simulated time', unit: 'D/U', default: 8, min: 1, max: 2000, group: 'Numerics', help: 'A steady sphere wake settles in 8–15; cylinder shedding needs 100–200 to develop and be measured' },
    { key: 'height', label: 'Domain height', unit: 'D', default: 3.2, min: 2.5, max: 24, group: 'Numerics', help: 'Sphere: the domain is 2H × H × H. Cylinder: 2H × H, stretched towards the wake axis; use 8 or more' },
    { key: 'n', label: 'Cells across the domain height', unit: '', default: 12, min: 8, max: 64, step: 1, discrete: true, group: 'Numerics', help: 'Sphere: 2n × n × n cells. Cylinder: 4n × n × 4' },
  ],
  defaults: (c) => ({ alt_m: c.atm.alt_m, dISA: c.atm.dISA_K }),
  run(i, ctx) { return bluffMemo(i, ctx); },
  runNow(i, ctx) {
    const a = isa(i.alt_m, i.dISA), n = cellsOf(i.n, 8, 64), cyl = i.body === BLUFF_BODIES[1], Re = i.Re, Hd = N.clamp(i.height, 2.5, 24), les = i.model === BLUFF_MODELS[1];
    const r = runBluff({ n, cyl, Re, height: Hd, tEnd: i.t_end, les, progress: ctx?.progress }), s = r.s, H = s.history, d = s.diagnostics(), g = s.grid, q = 0.5 * r.area, warnings = [];
    const cdH = H.fx.map((v) => v / q), clH = H.fy.map((v) => v / q), k0 = Math.floor((cyl ? 0.5 : 0.75) * H.t.length), tT = H.t.slice(k0), CD = N.mean(cdH.slice(k0)), CLm = N.mean(clH.slice(k0)), CLrms = N.std(clH.slice(k0)), St = crossingFrequency(tT, clH.slice(k0));
    const shed = Number.isFinite(St) && CLrms > 0.01, CDc = CD * (1 - r.blockage) ** 2, CDref = cyl ? NaN : cdSchillerNaumann(Re), StRef = cyl ? stRoshko(Re) : NaN, dPerCell = 1 / Math.max(g.dx, N.amin(g.dy)), V = (Re * a.nu) / i.D_m;
    const kz = Math.floor(g.nz / 2), su = s.slice('z', kz, 'u', { max: 4096 }), uc = N.range(g.nx, (c) => 0.5 * (su.z[(g.ny - 1) >> 1][c] + su.z[g.ny >> 1][c])), ol = s.outline('z', kz);
    let xr = NaN; for (let c = 0; c + 1 < g.nx; c++) if (g.x[c] > 0.5 && uc[c] < 0 && uc[c + 1] >= 0) xr = g.x[c] + ((0 - uc[c]) / (uc[c + 1] - uc[c])) * g.dx - 0.5;
    if (s.diverged || !Number.isFinite(CD)) warnings.push('The time march diverged.');
    if (n <= 20) warnings.push(COARSE_NOTE(`${g.nx} × ${g.ny} × ${g.nz} cells`, 'Cells across the domain height'));
    warnings.push(`Resolution: ${dPerCell.toFixed(1)} cells per diameter; the boundary layer is about D/√Re = ${(1 / Math.sqrt(Re)).toFixed(2)} D thick, i.e. ${(dPerCell / Math.sqrt(Re)).toFixed(1)} cells. ${dPerCell / Math.sqrt(Re) < 1.5 ? 'With fewer than about two cells in the boundary layer the immersed-boundary drag is typically 10–25% high.' : ''} Blockage is ${(100 * r.blockage).toFixed(1)}% between slip walls; the corrected coefficient uses the continuity estimate CD·(1 − blockage)², which is approximate.`);
    if (cyl && Re > 47 && !shed) warnings.push(`No periodic vortex shedding was detected in ${i.t_end} D/U. Shedding behind a cylinder takes 50–150 D/U to grow from the initial disturbance: run for 150–200 D/U to measure the Strouhal number${Hd < 8 ? ', and use a domain height of 8 D or more' : ''}.`);
    if (cyl && Re > 190) warnings.push('Above Re ≈ 190 the real cylinder wake becomes three-dimensional (mode A and B instabilities). With four spanwise cells this is a quasi-two-dimensional solution, which over-predicts the fluctuating lift and the drag.');
    if (!cyl && Re > 800) warnings.push('The Schiller–Naumann drag law is only valid below Re ≈ 800; above that the sphere drag coefficient levels off near 0.4–0.5.');
    if (!cyl && Re > 270) warnings.push('Above Re ≈ 270 the sphere wake sheds hairpin vortices: the forces are unsteady and the averages need a long run (50 D/U or more).');
    if (!les && Re > 500) warnings.push('At this Reynolds number the wake is turbulent and far from resolved on this grid: switch to LES and refine.');
    if (i.t_end < 8 && !cyl) warnings.push('The run is shorter than the 8–10 D/U the wake needs to become steady: the drag is still falling from its impulsive-start value.');
    const lim = { max: 80 }, ovl = { overlay: ol, xlabel: 'x/D [-]', ylabel: 'y/D [-]' }, [tt, cdT] = thin(H.t, cdH, 400), [, clT] = thin(H.t, clH, 400), pInf = N.mean(s.slice('z', kz, 'p', { max: 4096 }).z.map((row) => row[0]));
    return {
      kpis: [
        { key: 'CD', label: 'Drag coefficient', value: CD, unit: '-', note: `mean over the last ${cyl ? 'half' : 'quarter'} of the run; frontal area ${cyl ? 'D × span' : 'πD²/4'}` },
        { key: 'CD_corrected', label: 'Drag coefficient corrected for blockage', value: CDc, unit: '-', status: cyl || Math.abs(CDc / CDref - 1) < 0.15 ? 'ok' : 'warn' },
        ...(cyl ? [
          { key: 'St', label: 'Strouhal number f·D/U', value: shed ? St : NaN, unit: '-', status: shed && Math.abs(St / StRef - 1) < 0.1 ? 'ok' : 'warn', note: shed ? 'from mean-crossings of the lift' : 'no shedding detected' },
          { key: 'St_ref', label: 'Strouhal number, Roshko correlation', value: StRef, unit: '-', note: '0.212·(1 − 21.2/Re) for 50 < Re < 200 (empirical); about 0.165 at Re 100 and 0.19–0.20 at Re 200' },
          { key: 'CL_rms', label: 'Fluctuating lift coefficient (rms)', value: CLrms, unit: '-' },
        ] : [
          { key: 'CD_ref', label: 'Drag coefficient, Schiller–Naumann', value: CDref, unit: '-', note: '24/Re·(1 + 0.15·Re^0.687), Re < 800 (empirical fit to the standard drag curve)' },
          { key: 'CD_err_pct', label: 'Error of the corrected drag', value: 100 * (CDc / CDref - 1), unit: '%', status: Math.abs(CDc / CDref - 1) < 0.15 ? 'ok' : 'warn' },
          { key: 'CD_stokes', label: 'Stokes drag coefficient 24/Re', value: 24 / Re, unit: '-', note: 'Creeping-flow limit, exact for Re ≪ 1' },
          { key: 'CL_mean', label: 'Mean side-force coefficient', value: CLm, unit: '-', note: 'Zero by symmetry while the wake is axisymmetric' },
        ]),
        { key: 'x_recirc', label: 'Recirculation length behind the body', value: Number.isFinite(xr) ? xr : 0, unit: 'D', note: Number.isFinite(xr) ? 'from the rear of the body to the wake stagnation point' : 'no reversed flow on the wake axis' },
        { key: 'cells_per_D', label: 'Cells per diameter', value: dPerCell, unit: '', status: dPerCell >= 16 ? 'ok' : dPerCell >= 6 ? 'warn' : 'bad' },
        { key: 'blockage_pct', label: 'Blockage', value: 100 * r.blockage, unit: '%', status: r.blockage < 0.05 ? 'ok' : 'warn' },
        { key: 'drag_N', label: cyl ? 'Drag per metre of span in air' : 'Drag force in air', value: CD * 0.5 * a.rho * V * V * (cyl ? i.D_m : (PI * i.D_m ** 2) / 4), unit: cyl ? 'N/m' : 'N', note: `at ${V.toPrecision(3)} m/s for D = ${i.D_m} m` },
        { key: 'shedding_Hz', label: 'Shedding frequency in air', value: shed ? (St * V) / i.D_m : 0, unit: 'Hz', note: shed ? '' : 'steady wake' },
        { key: 'residual', label: 'Final unsteadiness RMS |∂u/∂t|', value: d.residual, unit: 'U²/D' },
        { key: 'div_max', label: 'Largest cell divergence', value: d.divMax, unit: 'U/D', status: d.divMax < 1e-8 ? 'ok' : 'warn' },
        { key: 'steps', label: 'Time steps', value: r.steps, unit: '' },
        { key: 'cells', label: 'Grid cells', value: g.nx * g.ny * g.nz, unit: '', note: `${s.solidCells} inside the body` },
      ],
      plots: [
        { type: 'line', title: 'Force history', xlabel: 'Time [D/U]', ylabel: 'Coefficient [-]', series: [{ name: 'CD', x: tt.slice(1), y: cdT.slice(1) }, { name: cyl ? 'CL' : 'Side force', x: tt.slice(1), y: clT.slice(1) }], annotations: cyl ? [] : [{ y: CDref, label: 'Schiller–Naumann' }] },
        heatOf(`Velocity magnitude, centre plane at t = ${s.time.toFixed(1)} D/U`, s.slice('z', kz, 'speed', lim), '|V|/U [-]', ovl),
        heatOf('Spanwise vorticity, centre plane', s.slice('z', kz, 'wz', lim), 'ω_z·D/U [-]', { ...ovl, diverging: true, map: (v) => N.clamp(v, -6, 6) }),
        heatOf('Pressure coefficient, centre plane', s.slice('z', kz, 'p', lim), 'Cp [-]', { ...ovl, diverging: true, map: (v) => N.clamp(2 * (v - pInf), -1.5, 1.2) }),
        { type: 'line', title: 'Streamwise velocity on the wake axis', xlabel: 'x/D [-]', ylabel: 'u/U [-]', series: [{ name: 'Centreline velocity', x: g.x, y: uc }], annotations: Number.isFinite(xr) ? [{ x: xr + 0.5, label: 'Wake stagnation point' }] : [] },
      ],
      tables: [{ title: 'Reference data for this case', columns: ['Quantity', 'Computed', 'Reference', 'Source'], rows: cyl ? [['Strouhal number', shed ? +St.toFixed(4) : 'not detected', +StRef.toFixed(4), 'Roshko (1954) correlation'], ['Drag coefficient', +CDc.toFixed(3), Re >= 80 && Re <= 250 ? '1.3–1.4' : '—', 'Measurements and 2-D simulations, Re 100–200']] : [['Drag coefficient (blockage-corrected)', +CDc.toFixed(3), +CDref.toFixed(3), 'Schiller & Naumann (1933)'], ['Drag coefficient (raw)', +CD.toFixed(3), +CDref.toFixed(3), 'unbounded stream']] }],
      warnings,
      models: [les ? 'Filtered incompressible Navier–Stokes with the WALE subgrid model' : 'Incompressible Navier–Stokes without a turbulence model', 'Direct-forcing immersed boundary on a Cartesian staggered grid; force from the momentum removed by the forcing', 'Second-order central differences, three-stage Runge–Kutta projection', cyl ? 'Pressure by cosine transform in x, Fourier transform in z and a tridiagonal solve in y' : 'Pressure by cosine transforms in x and z and a tridiagonal solve in y'],
      assumptions: ['Uniform inflow, convective outflow, slip walls at the sides' + (cyl ? ', periodic along the span (four cells: quasi-two-dimensional)' : ''), 'Impulsive start from uniform flow' + (cyl ? ' with a small seeded disturbance behind the cylinder to trigger shedding' : ''), 'The surface is located to within about one cell', 'Blockage correction by continuity only'],
    };
  },
  convergence: { param: 'n', label: 'Cells across the domain height', levels: [12, 16, 20, 24], metric: 'CD_corrected' },
  verify() {
    // steady momentum balance in a periodic array of spheres: the body force on the whole box equals the force on the sphere
    const G = 1, pa = createSolver({ n: [12, 12, 12], L: [1, 1, 1], nu: 1, forcing: { gradient: [G, 0, 0] }, body: { sdf: (x, y, z) => Math.hypot(x - 0.5, y - 0.5, z - 0.5) - 0.3 }, time: 'ab2', dt: 4e-4 }); pa.run(2000);
    const f = pa.monitor().force, ts = N.linspace(0, 40, 801), r = runBluff({ n: 16, cyl: false, Re: 40, height: 4, tEnd: 8 }), Hh = r.s.history, cd = N.mean(tailOf(Hh.fx, 0.2)) / (0.5 * r.area) * (1 - r.blockage) ** 2;
    return [
      N.check('Force on a sphere in a periodic array equals the driving body force', f[0], G, 2e-3, 'Global momentum balance at steady state: F = (−dp/dx)·V_box'),
      N.check('Periodic array: no side force', 1 + Math.hypot(f[1], f[2]), 1, 1e-9, 'Symmetry'),
      N.check('Sphere drag at Re = 40 on a coarse grid (4 cells per diameter)', cd, cdSchillerNaumann(40), 0.15, 'Schiller–Naumann 1.74; blockage-corrected immersed-boundary force'),
      N.check('Schiller–Naumann tends to Stokes drag', cdSchillerNaumann(1e-3) * 1e-3 / 24, 1, 2e-3, 'Stokes (1851): CD = 24/Re'),
      N.check('Frequency estimator on a known signal', crossingFrequency(ts, ts.map((t) => Math.sin(TWO_PI * 0.2 * t + 0.3) + 0.5)), 0.2, 1e-3, 'Mean-crossing period of sin(2π·0.2·t)'),
      N.check('Roshko Strouhal number at Re = 100', stRoshko(100), 0.167, 0.01, 'Roshko (1954); measurements give 0.164–0.167'),
    ];
  },
  validation: [
    { name: 'Sphere drag, steady axisymmetric regime', source: 'Schiller & Naumann (1933) fit to the standard drag curve: CD = 24/Re·(1 + 0.15·Re^0.687)', inputs: { body: BLUFF_BODIES[0], n: 32, t_end: 12, height: 5 }, sweep: { key: 'Re', values: [20, 40, 100] }, target: 'CD_corrected', observed: [cdSchillerNaumann(20), cdSchillerNaumann(40), cdSchillerNaumann(100)], tol_pct: 12 },
    { name: 'Cylinder vortex-shedding frequency', source: 'Roshko (1954), Williamson (1989): St ≈ 0.165 at Re = 100, 0.18 at Re = 150', inputs: { body: BLUFF_BODIES[1], n: 32, t_end: 200, height: 10 }, sweep: { key: 'Re', values: [100, 150] }, target: 'St', observed: [0.165, 0.184], tol_pct: 10 },
  ],
  recommend(res, i) {
    const o = res.outputs, out = [], cyl = i.body === BLUFF_BODIES[1];
    if (o.cells_per_D < 8) out.push({ severity: 'advise', title: 'Coarse body resolution', detail: `${o.cells_per_D.toFixed(1)} cells per diameter.`, action: 'Use 16 or more cells per diameter for a drag error below about 5% and run the mesh-convergence study; the immersed boundary converges at first order.', basis: 'Immersed-boundary accuracy' });
    if (!cyl && Math.abs(o.CD_err_pct) > 15) out.push({ severity: 'warn', title: 'Drag differs from the reference', detail: `Corrected CD = ${o.CD_corrected.toFixed(3)} against ${o.CD_ref.toFixed(3)}.`, action: 'Refine the grid, enlarge the domain and run longer; the difference at this resolution is numerical.', basis: 'Schiller–Naumann drag law' });
    if (cyl && o.shedding_Hz > 0) out.push({ severity: 'info', title: 'Vortex shedding excites structures and makes tones', detail: `Shedding at ${o.shedding_Hz.toFixed(1)} Hz for this diameter and speed (St = ${o.St.toFixed(3)}).`, action: 'Keep the frequency away from the natural frequencies of struts, antennas, probes and landing-gear legs (Suite 10); the same frequency is the Aeolian tone in Suite 11. Helical strakes or a fairing suppress it.', basis: 'Strouhal relation f = St·U/D' });
    if (o.blockage_pct > 5) out.push({ severity: 'advise', title: 'Domain is tight', detail: `Blockage ${o.blockage_pct.toFixed(1)}%.`, action: 'Increase the domain height (and the cell count with it) until the corrected and raw coefficients agree within a few percent.', basis: 'Wall interference' });
    return out;
  },
};

const bluffMemo = memoLast((i, ctx) => bluff3d.runNow(i, ctx));

export default {
  id: 'cfd', n: 1,
  tagline: 'How much lift and drag the aircraft makes, where the flow separates or shocks, and how fine a mesh the answer needs.',
  analyses: [airfoil, wing, dragBuildup, rans2d, wingRans, shockTube, shocks, cavityFlow, wallTurb, dns3d, les3d, rans3d, bluff3d],
  consumes: [],
  provides: [
    { key: 'CL', label: 'Lift coefficient', unit: '-' }, { key: 'CD', label: 'Drag coefficient', unit: '-' }, { key: 'CD0', label: 'Zero-lift drag coefficient', unit: '-' },
    { key: 'k_induced', label: 'Induced drag factor', unit: '-' }, { key: 'CLa_per_rad', label: 'Lift-curve slope', unit: '1/rad' }, { key: 'CLmax', label: 'Maximum lift coefficient', unit: '-' },
    { key: 'LD_max', label: 'Maximum L/D', unit: '-' }, { key: 'Cm_ac', label: 'Moment about aerodynamic centre', unit: '-' }, { key: 'e_oswald', label: 'Oswald efficiency', unit: '-' },
    { key: 'cp_min', label: 'Minimum pressure coefficient', unit: '-' }, { key: 'x_cp_frac', label: 'Centre of pressure', unit: 'x/c' },
  ],
  handoff: [
    { model: 'Body-fitted, wall-resolved 3-D RANS of the complete aircraft at flight Reynolds number (realizable k–ε, k–ω SST, Reynolds-stress, transition SST, γ–Reθ) on 10⁷–10⁸ cells', why: 'Natively, viscous wing forces come from body-fitted 2-D section RANS (Spalart–Allmaras, wall-resolved or wall functions) coupled strip by strip to the vortex lattice: that is a quasi-3-D method without spanwise boundary-layer flow, junctions, nacelles or shocks. 3-D RANS (Spalart–Allmaras, mixing length), LES and DNS are solved natively only on Cartesian immersed-boundary grids, for flow visualisation, up to the resolution the device allows (about 128³ cells). Such a grid places the wall to within one cell and cannot carry a y⁺ ≈ 1 boundary-layer mesh on the real surface, so drag to a few counts, maximum lift, separation onset and transition still need a body-fitted mesh and hours of parallel computing. Only the Spalart–Allmaras and mixing-length closures are implemented', tool: 'Finite-volume RANS solver (SU2, OpenFOAM, CFL3D, Fluent, STAR-CCM+). The High-fidelity bridge page of this app exports ready-to-run SU2 and OpenFOAM cases for the same geometry and flight condition; the wall-layer analysis gives the first-cell height' },
    { model: 'Industrial wall-resolved and wall-modelled LES, DES / DDES / IDDES and dynamic Smagorinsky on the real geometry; DNS at flight Reynolds number', why: 'LES (Smagorinsky, WALE, equilibrium wall model) and DNS run natively for canonical flows and immersed bodies at modest Reynolds number. Scale-resolving simulation of an aircraft needs 10⁸–10¹¹ body-fitted cells and 10⁵–10⁷ time steps; hybrid RANS–LES and the dynamic procedure are not implemented', tool: 'HPC scale-resolving solver (OpenFOAM, PyFR, CharLES, Nek5000), prepared through the High-fidelity bridge' },
    { model: 'Compressible transonic 3-D Euler / RANS with shock–boundary-layer interaction and buffet; Favre-averaged equations', why: 'The native 3-D solver is incompressible. Shocks on the wing are represented only by the Korn equation and the 1-D shock relations; the panel and lattice methods are linear subsonic', tool: 'Compressible RANS solver (SU2, OpenFOAM rhoSimpleFoam / HiSA) from the High-fidelity bridge, validated on NASA Common Research Model data' },
    { model: 'High-lift (slats, flaps, spoilers) aerodynamics', why: 'Multi-element gaps and confluent boundary layers are far below the cell size of a Cartesian immersed-boundary grid and beyond single-element integral methods', tool: 'Body-fitted RANS validated against the High-Lift Prediction Workshop cases' },
    { model: 'Moving, overset and ALE meshes; actuator-line and actuator-disc sources; free-vortex rotor wakes', why: 'The immersed boundary is implemented for stationary rigid bodies only; rotor wakes are treated by Suites 6 and 8', tool: 'Overset RANS (OVERFLOW, HELIOS, OpenFOAM overset) or free-wake codes (CAMRAD II, CHARM)' },
    { model: 'Coupled CFD–structure displacement and Eulerian–Lagrangian particle transport', why: 'Static and dynamic aeroelasticity are solved with reduced-order aerodynamics in Suite 3; droplet transport in Suite 13 uses potential flow', tool: 'Coupled CFD–CSD frameworks; icing CFD (LEWICE3D, FENSAP-ICE)' },
    { model: 'Real-gas equations of state and high-enthalpy flow', why: 'All compressible models here assume a calorically perfect gas', tool: 'Equilibrium or finite-rate chemistry solvers (DPLR, US3D)' },
    { model: 'Body-fitted volume meshing with prism layers (CGNS, SU2, OpenFOAM meshes) and watertight CAD repair', why: 'An imported triangle surface is flown directly by voxelising it into the Cartesian grid, which needs no mesher; boundary-layer volume meshes are not generated in the browser', tool: 'Gmsh, Pointwise or snappyHexMesh, driven by the case files from the High-fidelity bridge' },
    { model: 'Data-assisted / physics-informed turbulence closures', why: 'Requires training data and model governance outside the scope of an offline tool', tool: 'Research frameworks coupled to a RANS solver' },
  ],
};
