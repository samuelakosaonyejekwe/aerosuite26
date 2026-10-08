// Suite 1 — Aerodynamics and Computational Fluid Dynamics.
// Native solvers: 2-D linear-vortex panel method with integral boundary layer, 3-D vortex-lattice and lifting-line
// wing, component drag build-up with Korn wave drag, 1-D finite-volume Euler (shock tube) with an exact Riemann
// solver, shock/expansion relations, a 2-D incompressible Navier–Stokes projection solver (lid-driven cavity) and a
// 1-D Spalart–Allmaras wall-turbulence solve with a first-cell-height calculator.
// Volume-mesh RANS/LES/DNS of the complete aircraft are not solved here: see `handoff` at the end of the file.

import * as N from '../core/numerics.js';
import { isa, G0, GAMMA } from '../core/atmosphere.js';

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
  return { pn, yE, ym, dy, S, AR, chord, xle, Ga, G0, a: post(Ga), z: post(G0), trefftz, n2, nC };
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

export default {
  id: 'cfd', n: 1,
  tagline: 'How much lift and drag the aircraft makes, where the flow separates or shocks, and how fine a mesh the answer needs.',
  analyses: [airfoil, wing, dragBuildup, shockTube, shocks, cavityFlow, wallTurb],
  consumes: [],
  provides: [
    { key: 'CL', label: 'Lift coefficient', unit: '-' }, { key: 'CD', label: 'Drag coefficient', unit: '-' }, { key: 'CD0', label: 'Zero-lift drag coefficient', unit: '-' },
    { key: 'k_induced', label: 'Induced drag factor', unit: '-' }, { key: 'CLa_per_rad', label: 'Lift-curve slope', unit: '1/rad' }, { key: 'CLmax', label: 'Maximum lift coefficient', unit: '-' },
    { key: 'LD_max', label: 'Maximum L/D', unit: '-' }, { key: 'Cm_ac', label: 'Moment about aerodynamic centre', unit: '-' }, { key: 'e_oswald', label: 'Oswald efficiency', unit: '-' },
    { key: 'cp_min', label: 'Minimum pressure coefficient', unit: '-' }, { key: 'x_cp_frac', label: 'Centre of pressure', unit: 'x/c' },
  ],
  handoff: [
    { model: '3-D RANS on the full aircraft (k–ε, realizable k–ε, k–ω, k–ω SST, Reynolds-stress, transition SST, γ–Reθ)', why: 'Needs a body-fitted volume mesh of 10⁷–10⁸ cells and hours of parallel computing; only the 1-D Spalart–Allmaras wall layer and 2-D laminar Navier–Stokes are solved here', tool: 'Finite-volume RANS solver (SU2, OpenFOAM, CFL3D, Fluent, STAR-CCM+) using the mesh guidance from the wall-layer analysis' },
    { model: 'LES, wall-modelled LES, DES / DDES / IDDES, dynamic Smagorinsky, WALE, DNS', why: 'Scale-resolving simulation needs 10⁸–10¹¹ cells and 10⁵–10⁷ time steps', tool: 'HPC scale-resolving solver (OpenFOAM, PyFR, CharLES, Nek5000)' },
    { model: 'Transonic full-potential / Euler / RANS with shock–boundary-layer interaction and buffet', why: 'The panel and lattice methods are linear subsonic; shocks on the wing are represented only by the Korn equation and the 1-D shock relations', tool: 'Transonic RANS solver, validated on NASA Common Research Model data' },
    { model: 'High-lift (slats, flaps, spoilers) aerodynamics', why: 'Multi-element confluent boundary layers and separation are beyond single-element integral methods', tool: 'RANS validated against the High-Lift Prediction Workshop cases' },
    { model: 'Moving, overset and ALE meshes; actuator-line and actuator-disc sources; free-vortex rotor wakes; immersed boundaries', why: 'No general volume-mesh infrastructure in the browser; rotor wakes are treated by Suites 6 and 8', tool: 'Overset RANS (OVERFLOW, HELIOS) or free-wake codes (CAMRAD II, CHARM)' },
    { model: 'Coupled CFD–structure displacement and Eulerian–Lagrangian particle transport', why: 'Static and dynamic aeroelasticity are solved with reduced-order aerodynamics in Suite 3; droplet transport in Suite 13 uses potential flow', tool: 'Coupled CFD–CSD frameworks; icing CFD (LEWICE3D, FENSAP-ICE)' },
    { model: 'Real-gas equations of state and high-enthalpy flow', why: 'All compressible models here assume a calorically perfect gas', tool: 'Equilibrium or finite-rate chemistry solvers (DPLR, US3D)' },
    { model: 'CAD and mesh import (STEP, IGES, STL, CGNS, SU2, OpenFOAM) with watertight repair and volume meshing', why: 'Geometry is parametric (planform and NACA sections); a geometry kernel and mesher cannot be shipped dependency-free', tool: 'CAD kernel with Gmsh, Pointwise or snappyHexMesh' },
    { model: 'Data-assisted / physics-informed turbulence closures', why: 'Requires training data and model governance outside the scope of an offline tool', tool: 'Research frameworks coupled to a RANS solver' },
  ],
};
