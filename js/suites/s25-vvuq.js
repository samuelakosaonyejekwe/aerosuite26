// Suite 25 — Numerical Verification, Experimental Validation and Uncertainty Quantification.
// The V&V / UQ laboratory: manufactured-solution code verification with observed order, a grid-convergence (GCI)
// calculator, validation metrics with validation uncertainty, uncertainty propagation by Monte Carlo, Latin
// hypercube, polynomial chaos and stochastic collocation, Sobol and Morris sensitivity, Bayesian calibration by
// MCMC, and a multifidelity discrepancy model with a reproducibility record.
// (Every analysis of every other suite already gets generic convergence, sensitivity, UQ, calibration and
// validation studies from the application core; this suite holds the self-contained experiments and data tools.)

import * as N from '../core/numerics.js';
import { G0 } from '../core/atmosphere.js';

// ---- text parsing ---------------------------------------------------------------------------
/** Parse rows of numbers ("a, b, c" per line or separated by ";"). Returns rows with at least minCols finite numbers. */
function parseTable(text, minCols) {
  const rows = []; let bad = 0;
  for (const line of String(text || '').split(/[;\n]+/)) { if (!line.trim()) continue; const tok = line.split(/[,\s]+/).filter(Boolean), v = tok.map(Number); if (v.length >= minCols && v.slice(0, minCols).every(Number.isFinite)) rows.push(v); else bad++; }
  return { rows, bad };
}
/** FNV-1a 32-bit hash of a string as 8 hex digits. */
function fnv1a(str) { let h = 0x811c9dc5; for (let k = 0; k < str.length; k++) { h ^= str.charCodeAt(k) & 0xff; h = Math.imul(h, 0x01000193); } return (h >>> 0).toString(16).padStart(8, '0'); }
const stableJson = (o) => JSON.stringify(Object.keys(o).sort().map((k) => [k, typeof o[k] === 'number' ? Number(o[k].toPrecision(12)) : o[k]]));
const rms = (a) => Math.sqrt(N.mean(a.map((v) => v * v)));
const VERSION = 'vvuq-suite 1.0 / mulberry32 RNG / Halton-shift QMC';

// ---- manufactured-solution solvers -----------------------------------------------------------
/** 1-D steady advection–diffusion a·u' − ν·u'' = s on [0,1], manufactured u = sin(2πx) + x². */
function mmsAdvDiff(n, a, nu, scheme) {
  const h = 1 / n, um = (x) => Math.sin(2 * Math.PI * x) + x * x, src = (x) => a * (2 * Math.PI * Math.cos(2 * Math.PI * x) + 2 * x) - nu * (-4 * Math.PI * Math.PI * Math.sin(2 * Math.PI * x) + 2);
  const m = n - 1, lo = new Array(m), di = new Array(m), up = new Array(m), rhs = new Array(m), upw = scheme.startsWith('upwind');
  for (let j = 0; j < m; j++) {
    const x = (j + 1) * h, d = nu / (h * h);
    if (upw) { const ap = Math.max(a, 0) / h, am = Math.min(a, 0) / h; lo[j] = -d - ap; di[j] = 2 * d + ap - am; up[j] = -d + am; } else { lo[j] = -d - a / (2 * h); di[j] = 2 * d; up[j] = -d + a / (2 * h); }
    rhs[j] = src(x); if (j === 0) rhs[j] -= lo[j] * um(0); if (j === m - 1) rhs[j] -= up[j] * um(1);
  }
  const u = N.solveTridiag(lo, di, up, rhs), xs = N.range(m, (j) => (j + 1) * h), err = u.map((v, j) => v - um(xs[j]));
  return { h, x: [0, ...xs, 1], u: [um(0), ...u, um(1)], exact: [um(0), ...xs.map(um), um(1)], err, functional: h * (N.sum(u) + 0.5 * (um(0) + um(1))), resid: [] };
}
/** 2-D Poisson −∇²u = s on the unit square, manufactured u = sin(πx)·sin(2πy) + x·y, five-point stencil, conjugate gradients. */
function mmsPoisson2D(n, tol = 1e-11) {
  const h = 1 / n, m = n - 1, M = m * m, um = (x, y) => Math.sin(Math.PI * x) * Math.sin(2 * Math.PI * y) + x * y, src = (x, y) => 5 * Math.PI * Math.PI * Math.sin(Math.PI * x) * Math.sin(2 * Math.PI * y);
  const b = new Array(M), id = (i, j) => j * m + i;
  for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) { const x = (i + 1) * h, y = (j + 1) * h; let v = h * h * src(x, y); if (i === 0) v += um(0, y); if (i === m - 1) v += um(1, y); if (j === 0) v += um(x, 0); if (j === m - 1) v += um(x, 1); b[id(i, j)] = v; }
  const Ax = (p, out) => { for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) { const k = id(i, j); out[k] = 4 * p[k] - (i > 0 ? p[k - 1] : 0) - (i < m - 1 ? p[k + 1] : 0) - (j > 0 ? p[k - m] : 0) - (j < m - 1 ? p[k + m] : 0); } };
  const u = new Array(M).fill(0), r = b.slice(), p = b.slice(), Ap = new Array(M); let rr = N.dot(r, r); const r0 = Math.sqrt(rr) || 1, resid = [1];
  for (let it = 0; it < 4 * M && Math.sqrt(rr) / r0 > tol; it++) { Ax(p, Ap); const al = rr / N.dot(p, Ap); for (let k = 0; k < M; k++) { u[k] += al * p[k]; r[k] -= al * Ap[k]; } const rn = N.dot(r, r), be = rn / rr; for (let k = 0; k < M; k++) p[k] = r[k] + be * p[k]; rr = rn; resid.push(Math.sqrt(rr) / r0); }
  const err = new Array(M); let fsum = 0; for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) { const k = id(i, j); err[k] = u[k] - um((i + 1) * h, (j + 1) * h); fsum += u[k]; }
  // mean of u over the square by the trapezoid rule (boundary values are known exactly)
  let bsum = 0; for (let k = 0; k <= n; k++) { const t = k * h, w = k === 0 || k === n ? 0.5 : 1; bsum += w * (um(t, 0) + um(t, 1)) * 0.5 + (k > 0 && k < n ? 0.5 * (um(0, t) + um(1, t)) : 0); }
  const jm = Math.floor(m / 2);
  return { h, x: N.range(m, (i) => (i + 1) * h), u: N.range(m, (i) => u[id(i, jm)]), exact: N.range(m, (i) => um((i + 1) * h, (jm + 1) * h)), err, functional: h * h * (fsum + bsum), resid, field: N.range(m, (j) => N.range(m, (i) => err[id(i, j)])) };
}
/** Time integration of the semi-discrete heat equation on a fixed grid; the exact semi-discrete solution isolates the temporal error. */
function mmsTime(nt, method, nu = 1, T = 0.1, nx = 12) {
  const h = 1 / nx, m = nx - 1, dt = T / nt, c = nu / (h * h), lam = -c * 4 * Math.sin((Math.PI * h) / 2) ** 2, xs = N.range(m, (j) => (j + 1) * h), F = (u) => u.map((v, j) => c * ((j > 0 ? u[j - 1] : 0) - 2 * v + (j < m - 1 ? u[j + 1] : 0)));
  let u = xs.map((x) => Math.sin(Math.PI * x)); const ax = (a, b, s) => a.map((v, j) => v + s * b[j]);
  const imp = (th) => { const lo = new Array(m).fill(-th * dt * c), di = new Array(m).fill(1 + 2 * th * dt * c), rhs = ax(u, F(u), (1 - th) * dt); return N.solveTridiag(lo, di, lo, rhs); };
  for (let s = 0; s < nt; s++) {
    if (method === 'Euler explicit') u = ax(u, F(u), dt);
    else if (method === 'RK2 (Heun)') { const k1 = F(u), k2 = F(ax(u, k1, dt)); u = u.map((v, j) => v + 0.5 * dt * (k1[j] + k2[j])); }
    else if (method === 'RK4') { const k1 = F(u), k2 = F(ax(u, k1, dt / 2)), k3 = F(ax(u, k2, dt / 2)), k4 = F(ax(u, k3, dt)); u = u.map((v, j) => v + (dt / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])); }
    else u = imp(method === 'Crank–Nicolson' ? 0.5 : 1);
  }
  const ex = xs.map((x) => Math.exp(lam * T) * Math.sin(Math.PI * x));
  return { h: dt, x: xs, u, exact: ex, err: u.map((v, j) => v - ex[j]), functional: u[Math.floor(m / 2)], resid: [], stable: Math.abs(lam) * dt };
}
const FORMAL = { 'central (2nd order)': 2, 'upwind (1st order)': 1, 'Euler explicit': 1, 'Euler implicit': 1, 'RK2 (Heun)': 2, 'Crank–Nicolson': 2, 'RK4': 4 };

// ---- quadrature rules and orthonormal polynomials ---------------------------------------------
/** Gauss rule for a probability measure by Golub–Welsch: 'normal' (probabilists' Hermite) or 'uniform' on [−1, 1] (Legendre). */
function gaussRule(n, kind) {
  if (n === 1) return { x: [0], w: [1] };
  const J = N.zeros(n); for (let k = 1; k < n; k++) J[k][k - 1] = J[k - 1][k] = kind === 'normal' ? Math.sqrt(k) : k / Math.sqrt(4 * k * k - 1);
  const e = N.eigSym(J); return { x: e.values, w: e.vectors.map((v) => v[0] * v[0]) };
}
/** Orthonormal polynomial values ψ₀…ψ_p at ξ (Hermite for normal, Legendre for uniform). */
function orthoPoly(xi, p, kind) {
  const v = [1]; if (p >= 1) v.push(xi);
  if (kind === 'normal') { for (let n = 1; n < p; n++) v.push(xi * v[n] - n * v[n - 1]); let f = 1; return v.map((q, n) => { if (n > 0) f *= n; return q / Math.sqrt(f); }); }
  for (let n = 1; n < p; n++) v.push(((2 * n + 1) * xi * v[n] - n * v[n - 1]) / (n + 1));
  return v.map((q, n) => q * Math.sqrt(2 * n + 1));
}
function multiIndex(d, p) { const out = []; const rec = (j, left, cur) => { if (j === d) return out.push(cur.slice()); for (let k = 0; k <= left; k++) { cur[j] = k; rec(j + 1, left - k, cur); } }; rec(0, p, new Array(d).fill(0)); return out.sort((a, b) => N.sum(a) - N.sum(b)); }
/** Polynomial chaos by least-squares regression. xis: samples in standard variables; y: model values. */
function pceFit(xis, y, p, kind) {
  const d = xis[0].length, idx = multiIndex(d, p), Phi = xis.map((xi) => { const pv = xi.map((x) => orthoPoly(x, p, kind)); return idx.map((a) => a.reduce((s, k, j) => s * pv[j][k], 1)); }), c = N.lstsq(Phi, y, 1e-12);
  const varT = N.sum(c.slice(1).map((v) => v * v)), S1 = N.range(d, (j) => N.sum(idx.map((a, k) => (k > 0 && a.every((q, m) => (m === j ? q > 0 : q === 0)) ? c[k] * c[k] : 0))) / (varT || 1)), ST = N.range(d, (j) => N.sum(idx.map((a, k) => (a[j] > 0 ? c[k] * c[k] : 0))) / (varT || 1));
  return { mean: c[0], sd: Math.sqrt(varT), S1, ST, terms: idx.length, c, predict: (xi) => { const pv = xi.map((x) => orthoPoly(x, p, kind)); return N.sum(idx.map((a, k) => c[k] * a.reduce((s, q, j) => s * pv[j][q], 1))); } };
}
/** Tensor-grid stochastic collocation: mean and standard deviation by Gauss quadrature with q points per dimension. */
function collocate(f, d, q, kind) {
  const r = gaussRule(q, kind); let m1 = 0, m2 = 0, n = 0; const xi = new Array(d);
  const rec = (j, w) => { if (j === d) { const v = f(xi); m1 += w * v; m2 += w * v * v; n++; return; } for (let k = 0; k < q; k++) { xi[j] = r.x[k]; rec(j + 1, w * r.w[k]); } };
  rec(0, 1); return { mean: m1, sd: Math.sqrt(Math.max(0, m2 - m1 * m1)), n };
}
// ---- test models ------------------------------------------------------------------------------
// each returns { d, kind, names, f(ξ) in standard variables, exact: {mean, sd, S1, ST} | null }
function testModel(i) {
  if (i.model === 'Ishigami function') {
    const a = 7, b = 0.1, P = Math.PI, V = a * a / 8 + (b * P ** 4) / 5 + (b * b * P ** 8) / 18 + 0.5, V1 = 0.5 * (1 + (b * P ** 4) / 5) ** 2, V2 = a * a / 8, V13 = b * b * P ** 8 * (1 / 18 - 1 / 50);
    return { d: 3, kind: 'uniform', names: ['x1', 'x2', 'x3'], unit: '-', f: (xi) => { const x = xi.map((v) => P * v); return Math.sin(x[0]) + a * Math.sin(x[1]) ** 2 + b * x[2] ** 4 * Math.sin(x[0]); }, exact: { mean: a / 2, sd: Math.sqrt(V), S1: [V1 / V, V2 / V, 0], ST: [(V1 + V13) / V, V2 / V, V13 / V] } };
  }
  if (i.model === 'Polynomial test') return { d: 2, kind: 'normal', names: ['ξ1', 'ξ2'], unit: '-', f: (xi) => xi[0] * xi[0] + xi[0] * xi[1], exact: { mean: 1, sd: Math.sqrt(3), S1: [2 / 3, 0], ST: [1, 1 / 3] } };
  // Breguet range R = (V / (c·g))·(L/D)·ln(1/(1 − fuel fraction)) with normal L/D, TSFC and fuel fraction
  const mu = [i.LD, i.tsfc, i.fuel_frac], sd = [i.LD * i.cov_LD, i.tsfc * i.cov_tsfc, i.fuel_frac * i.cov_ff];
  return { d: 3, kind: 'normal', names: ['Lift-to-drag ratio', 'Specific fuel consumption', 'Fuel fraction'], unit: 'km', exact: null, f: (xi) => { const LD = mu[0] + sd[0] * xi[0], c = Math.max(1e-9, mu[1] + sd[1] * xi[1]), ff = N.clamp(mu[2] + sd[2] * xi[2], 1e-6, 0.95); return ((i.V / (c * G0)) * LD * Math.log(1 / (1 - ff))) / 1e3; } };
}
const toStd = (kind, u) => (kind === 'normal' ? N.normInv(N.clamp(u, 1e-12, 1 - 1e-12)) : 2 * u - 1);
const MODEL_IN = [
  { key: 'model', label: 'Model', type: 'select', options: ['Breguet range', 'Ishigami function', 'Polynomial test'], default: 'Breguet range', group: 'Model', help: 'Breguet range with uncertain L/D, fuel consumption and fuel fraction; the Ishigami and polynomial functions have exact statistics for checking the methods' },
  { key: 'V', label: 'Cruise speed', unit: 'm/s', default: 231, min: 1, group: 'Breguet inputs' }, { key: 'LD', label: 'Mean lift-to-drag ratio', unit: '-', default: 17, min: 1, max: 60, group: 'Breguet inputs' }, { key: 'cov_LD', label: 'L/D coefficient of variation', unit: '-', default: 0.04, min: 0, max: 0.3, group: 'Breguet inputs' },
  { key: 'tsfc', label: 'Mean thrust-specific fuel consumption', unit: 'kg/N/s', default: 1.62e-5, min: 1e-7, group: 'Breguet inputs', help: 'For propeller aircraft use BSFC·V/η' }, { key: 'cov_tsfc', label: 'Fuel consumption coefficient of variation', unit: '-', default: 0.03, min: 0, max: 0.3, group: 'Breguet inputs' },
  { key: 'fuel_frac', label: 'Mean cruise fuel fraction', unit: '-', default: 0.18, min: 0.005, max: 0.6, group: 'Breguet inputs' }, { key: 'cov_ff', label: 'Fuel fraction coefficient of variation', unit: '-', default: 0.03, min: 0, max: 0.3, group: 'Breguet inputs' },
];
const modelDefaults = (c, up, d) => { const V = c.mission.cruise_V_ms || c.flight.V_ms || 60, jet = c.prop.type === 'turbofan' || c.prop.type === 'turbojet', LDm = c.wing.S_m2 > 0 && d.k_induced > 0 ? 1 / (2 * Math.sqrt(c.aero.CD0 * d.k_induced)) : 4.5;
  return { V, LD: up.performance?.LD_max ? 0.92 * up.performance.LD_max : 0.92 * LDm, tsfc: jet ? up.propulsion?.tsfc_kg_Ns ?? c.prop.tsfc_kg_Ns : (c.prop.bsfc_kg_Ws * V) / (c.prop.eta_prop || 0.8), fuel_frac: c.mass.fuel_kg > 0 ? N.clamp((0.8 * c.mass.fuel_kg) / c.mass.mtow_kg, 0.01, 0.5) : undefined }; };
/** Shifted Halton points: n × d on [0,1). */
function halton(n, d, u) { const pr = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37], sh = N.range(d, () => u()); return N.range(n, (k) => N.range(d, (j) => { let f = 1, r = 0, q = k + 1; const b = pr[j % pr.length]; while (q > 0) { f /= b; r += f * (q % b); q = Math.floor(q / b); } return (r + sh[j]) % 1; })); }

// ---- analyses -------------------------------------------------------------------------------
const mms = {
  id: 'mms', title: 'Code verification by manufactured solutions', fidelity: 'numerical',
  summary: 'Solves model equations whose exact solution is known by construction, refines the grid or time step, and checks that the error falls at the rate the scheme promises. A solver that passes is doing the mathematics it claims to.',
  equations: ['Taylor series truncation error analysis', 'Numerical consistency relations', 'Stability analysis equations', 'Lax equivalence theorem', 'Richardson extrapolation equation', 'Residual convergence equations'],
  inputs: [
    { key: 'problem', label: 'Verification problem', type: 'select', options: ['1-D advection–diffusion (space)', '2-D Poisson (space)', 'Heat equation (time integrator)'], default: '1-D advection–diffusion (space)', group: 'Problem' },
    { key: 'scheme', label: 'Convection scheme (1-D problem)', type: 'select', options: ['central (2nd order)', 'upwind (1st order)'], default: 'central (2nd order)', group: 'Problem' },
    { key: 'integrator', label: 'Time integrator (heat problem)', type: 'select', options: ['Euler explicit', 'Euler implicit', 'RK2 (Heun)', 'Crank–Nicolson', 'RK4'], default: 'Crank–Nicolson', group: 'Problem' },
    { key: 'a', label: 'Advection speed', unit: 'm/s', default: 1, min: -50, max: 50, group: 'Problem', help: '0 gives the pure Poisson problem' }, { key: 'nu', label: 'Diffusivity', unit: 'm²/s', default: 0.1, min: 1e-4, max: 10, group: 'Problem' },
    { key: 'nBase', label: 'Cells (or steps) on the coarsest level', unit: '', default: 8, min: 4, max: 128, step: 1, discrete: true, group: 'Numerics' }, { key: 'nLevels', label: 'Refinement levels (×2 each)', unit: '', default: 5, min: 3, max: 8, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: () => ({}),
  run(i) {
    const warnings = [], kind = i.problem.startsWith('1-D') ? 0 : i.problem.startsWith('2-D') ? 1 : 2, nL = Math.round(i.nLevels), levels = [], formal = kind === 0 ? FORMAL[i.scheme] : kind === 1 ? 2 : FORMAL[i.integrator];
    const nb = Math.round(i.nBase), n0 = kind === 2 ? Math.max(nb, 32) : nb, maxN = kind === 1 ? 96 : 1 << 20;
    for (let l = 0; l < nL; l++) { const n = n0 * 2 ** l; if (n > maxN) { warnings.push(`Refinement stopped at ${n / 2} cells per side to keep the 2-D solve fast.`); break; } const r = kind === 0 ? mmsAdvDiff(n, i.a, i.nu, i.scheme) : kind === 1 ? mmsPoisson2D(n) : mmsTime(n, i.integrator); levels.push({ n, ...r, L1: N.mean(r.err.map(Math.abs)), L2: rms(r.err), Li: N.amax(r.err.map(Math.abs)) }); }
    const ord = (k, key) => (k > 0 && levels[k][key] > 0 && levels[k - 1][key] > 0 ? Math.log(levels[k - 1][key] / levels[k][key]) / Math.log(2) : NaN), L = levels.length, pObs = ord(L - 1, 'L2'), fin = levels[L - 1];
    const g = L >= 3 ? N.gci([levels[L - 1].h, levels[L - 2].h, levels[L - 3].h], [levels[L - 1].functional, levels[L - 2].functional, levels[L - 3].functional]) : null, hs = levels.map((q) => q.h), ref = hs.map((h) => levels[0].L2 * (h / hs[0]) ** formal);
    const pe = kind === 0 ? Math.abs(i.a) * levels[0].h / i.nu : 0;
    if (kind === 0 && !i.scheme.startsWith('upwind') && pe > 2) warnings.push(`Cell Péclet number is ${pe.toFixed(1)} on the coarsest grid: the central scheme oscillates above 2, so the coarse levels are outside the asymptotic range.`);
    if (kind === 2 && ['Euler explicit', 'RK2 (Heun)', 'RK4'].includes(i.integrator) && levels[0].stable > (i.integrator === 'RK4' ? 2.78 : 2)) warnings.push('The coarsest time step is outside the stability limit of this explicit integrator; errors there are meaningless.');
    if (fin.L2 < 1e-12) warnings.push('The finest-level error is at round-off: the last observed order is not reliable.');
    const okOrder = Math.abs(pObs - formal) < 0.15 * formal;
    return {
      kpis: [
        { key: 'mms_observed_order', label: 'Observed order of accuracy (L2, finest pair)', value: pObs, unit: '-', status: okOrder ? 'ok' : 'bad', note: `Formal order ${formal}` },
        { key: 'mms_formal_order', label: 'Formal order of the scheme', value: formal, unit: '-' }, { key: 'mms_order_deviation', label: 'Observed minus formal order', value: pObs - formal, unit: '-', status: okOrder ? 'ok' : 'bad' },
        { key: 'mms_l2_fine', label: 'L2 error on the finest level', value: fin.L2, unit: '-' }, { key: 'mms_linf_fine', label: 'L∞ error on the finest level', value: fin.Li, unit: '-' },
        { key: 'functional', label: 'Solution functional on the finest level', value: fin.functional, unit: '-', note: kind === 2 ? 'Mid-point value at the end time' : 'Domain mean of the solution' },
        { key: 'mms_gci_fine_pct', label: 'GCI of the functional (finest three levels)', value: g ? 100 * g.gciFine : NaN, unit: '%' }, { key: 'mms_order_functional', label: 'Observed order of the functional', value: g ? g.p : NaN, unit: '-' },
        { key: 'code_verified', label: 'Code verification passed', value: okOrder ? 1 : 0, unit: '', status: okOrder ? 'ok' : 'bad', note: 'Observed order within 15% of formal' },
        ...(kind === 1 ? [{ key: 'cg_iterations', label: 'Conjugate-gradient iterations (finest)', value: fin.resid.length - 1, unit: '' }, { key: 'iter_residual_drop_decades', label: 'Iterative residual reduction', value: -Math.log10(Math.max(fin.resid[fin.resid.length - 1], 1e-300)), unit: 'decades' }] : []),
      ],
      plots: [
        { type: 'line', title: 'Discretisation error versus resolution', xlabel: kind === 2 ? 'Time step [s]' : 'Cell size [m]', ylabel: 'Error norm [-]', xlog: true, ylog: true, series: [{ name: 'L1', x: hs, y: levels.map((q) => Math.max(q.L1, 1e-300)), style: 'line+points' }, { name: 'L2', x: hs, y: levels.map((q) => Math.max(q.L2, 1e-300)), style: 'line+points' }, { name: 'L∞', x: hs, y: levels.map((q) => Math.max(q.Li, 1e-300)), style: 'line+points' }, { name: `Formal slope ${formal}`, x: hs, y: ref.map((v) => Math.max(v, 1e-300)), style: 'dash' }] },
        { type: 'line', title: 'Numerical and manufactured solution (coarsest level)', xlabel: 'x [m]', ylabel: 'u [-]', series: [{ name: 'Numerical', x: levels[0].x, y: levels[0].u, style: 'line+points' }, { name: 'Manufactured (exact)', x: levels[0].x, y: levels[0].exact, style: 'dash' }] },
        ...(kind === 1 ? [{ type: 'line', title: 'Iterative convergence of the linear solver (finest level)', xlabel: 'Iteration [-]', ylabel: 'Relative residual [-]', ylog: true, series: [{ name: 'Conjugate gradients', x: fin.resid.map((_, k) => k).filter((k) => k % Math.max(1, Math.floor(fin.resid.length / 300)) === 0), y: fin.resid.filter((_, k) => k % Math.max(1, Math.floor(fin.resid.length / 300)) === 0).map((v) => Math.max(v, 1e-300)) }] }, { type: 'heat', title: 'Error field on the coarsest level', xlabel: 'x [m]', ylabel: 'y [m]', zlabel: 'u − u exact [-]', x: levels[0].x, y: levels[0].x, z: levels[0].field, contours: 10, diverging: true, equalAspect: true }] : []),
      ],
      tables: [{ title: 'Error norms and observed order', columns: ['Cells / steps', kind === 2 ? 'Δt' : 'h', 'L1', 'L2', 'L∞', 'Order L1', 'Order L2', 'Order L∞'], rows: levels.map((q, k) => [q.n, q.h, q.L1, q.L2, q.Li, ord(k, 'L1'), ord(k, 'L2'), ord(k, 'Li')].map((v) => (Number.isNaN(v) ? '—' : v))) }],
      warnings, models: ['Method of manufactured solutions', kind === 0 ? `Finite differences: ${i.scheme} convection, central diffusion` : kind === 1 ? 'Five-point Laplacian with conjugate-gradient solver' : `Method of lines with ${i.integrator}`, 'Observed order from successive error ratios', 'Grid convergence index of a solution functional'],
      assumptions: ['Uniform grids refined by a factor of two', 'Dirichlet boundary values taken from the manufactured solution', kind === 2 ? 'The reference is the exact solution of the semi-discrete system, so only the temporal error is measured' : 'Linear systems solved to a tolerance far below the discretisation error'],
    };
  },
  convergence: { param: 'nBase', label: 'Cells on the coarsest level', levels: [4, 8, 16, 32], metric: 'functional' },
  verify() {
    const b = { problem: '1-D advection–diffusion (space)', scheme: 'central (2nd order)', integrator: 'RK4', a: 1, nu: 0.1, nBase: 16, nLevels: 5 }, o = (x) => N.kv(mms.run({ ...b, ...x })).mms_observed_order;
    return [
      N.check('Central scheme converges at second order', o({}), 2, 0.02, 'Truncation error O(h²)'),
      N.check('Upwind scheme converges at first order', o({ scheme: 'upwind (1st order)', nBase: 64 }), 1, 0.05, 'Truncation error O(h)'),
      N.check('2-D five-point Poisson converges at second order', o({ problem: '2-D Poisson (space)', nBase: 16, nLevels: 3 }), 2, 0.03, 'Truncation error O(h²)'),
      N.check('RK4 converges at fourth order in time', o({ problem: 'Heat equation (time integrator)', nBase: 32, nLevels: 4 }), 4, 0.03, 'Classical Runge–Kutta'),
      N.check('Crank–Nicolson converges at second order in time', o({ problem: 'Heat equation (time integrator)', integrator: 'Crank–Nicolson' }), 2, 0.02, 'Trapezoidal rule'),
      N.check('Implicit Euler converges at first order in time', o({ problem: 'Heat equation (time integrator)', integrator: 'Euler implicit', nBase: 64 }), 1, 0.05, 'Backward Euler'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs;
    if (o.code_verified) return [{ severity: 'info', title: 'The solver reproduces its formal order', detail: `Observed ${o.mms_observed_order.toFixed(2)} against formal ${o.mms_formal_order}.`, action: 'Apply the same test to any new discretisation before using it for engineering results; then estimate the error of real solutions with the grid-convergence calculator. Verified coarse grids save computing energy.', basis: 'Order-of-accuracy test (the most rigorous code-verification criterion)' }];
    return [{ severity: 'warn', title: 'Observed order does not match the formal order', detail: `Observed ${Number(o.mms_observed_order).toFixed(2)} against formal ${o.mms_formal_order}. ${res.warnings.join(' ')}`, action: i.problem.startsWith('1-D') ? 'Refine until the cell Péclet number is below 2, or add levels: the coarse grids are not in the asymptotic range.' : 'Add refinement levels or reduce the time step into the stable, asymptotic range. If the mismatch persists on fine grids, the implementation has an error.', basis: 'Order-of-accuracy test' }];
  },
};

/** Least-squares fit f = f0 + α·h^p over several grids (golden search on p, linear solve for f0 and α). */
function lsqOrder(h, f) {
  const fit = (p) => { const A = h.map((v) => [1, v ** p]), c = N.lstsq(A, f), r = f.map((v, k) => v - c[0] - c[1] * h[k] ** p); return { c, sse: N.dot(r, r) }; };
  const p = N.goldenSection((q) => fit(q).sse, 0.3, 8, 1e-10), r = fit(p);
  return { p, f0: r.c[0], alpha: r.c[1], sd: Math.sqrt(r.sse / Math.max(1, h.length - 3)) };
}
const gciCalc = {
  id: 'gci', title: 'Solution verification: Richardson extrapolation and grid convergence index', fidelity: 'analytical',
  summary: 'Takes a result computed on three or more grids (or time steps) and estimates the observed order, the grid-independent value and the numerical uncertainty band; also grades an iterative residual history.',
  equations: ['Richardson extrapolation equation', 'Grid Convergence Index formulation', 'Residual convergence equations', 'Taylor series truncation error analysis'],
  inputs: [
    { key: 'grids', label: 'Grid size and result, one grid per line', type: 'text', default: '1, 0.97050\n2, 0.96854\n4, 0.96178', group: 'Grid study', help: '"h, value" with h a representative cell size (or 1/N^(1/dimension)) in any consistent unit; three or more grids. The default is a synthetic three-grid illustration.' },
    { key: 'Fs', label: 'Safety factor', unit: '-', default: 1.25, min: 1, max: 3, group: 'Grid study', help: '1.25 for three or more grids with an observed order; 3 for two grids' },
    { key: 'p_formal', label: 'Formal order of the scheme', unit: '-', default: 2, min: 0.5, max: 6, group: 'Grid study' },
    { key: 'residuals', label: 'Iterative residual history (optional)', type: 'text', default: '1, 0.31, 0.12, 0.043, 0.017, 6.1e-3, 2.3e-3, 8.8e-4, 3.3e-4, 1.2e-4, 4.7e-5, 1.8e-5, 6.6e-6', group: 'Iterative convergence', help: 'Residual norms in iteration order. The default list is synthetic.' },
    { key: 'resid_target', label: 'Required residual reduction', unit: 'decades', default: 5, min: 1, max: 14, group: 'Iterative convergence' },
  ],
  defaults: () => ({}),
  run(i) {
    const warnings = [], t = parseTable(i.grids, 2); let rows = t.rows.filter((r) => r[0] > 0).sort((a, b) => a[0] - b[0]);
    if (t.bad) warnings.push(`${t.bad} line(s) of the grid table could not be read and were skipped.`);
    if (rows.length < 3) { warnings.push('Fewer than three valid grids were entered; the built-in synthetic example was used.'); rows = [[1, 0.9705], [2, 0.96854], [4, 0.96178]]; }
    const h = rows.map((r) => r[0]), f = rows.map((r) => r[1]), g = N.gci(h.slice(0, 3), f.slice(0, 3), i.Fs), osc = !g.monotonic;
    let p = g.p, fEx = g.fExact, gf = g.gciFine, gc = g.gciCoarse;
    if (osc || !Number.isFinite(p) || p <= 0) { p = i.p_formal; const r21 = h[1] / h[0]; fEx = f[0] + (f[0] - f[1]) / (r21 ** p - 1); gf = (3 * Math.abs((f[1] - f[0]) / (f[0] || 1e-300))) / (r21 ** p - 1); gc = (3 * Math.abs((f[2] - f[1]) / (f[1] || 1e-300))) / ((h[2] / h[1]) ** p - 1); warnings.push('Convergence is oscillatory or the differences vanish: the observed order is not defined. The formal order with a safety factor of 3 was used instead.'); }
    else if (p > 2 * i.p_formal || p < 0.5 * i.p_formal) warnings.push(`Observed order ${p.toFixed(2)} is far from the formal order ${i.p_formal}: the grids are probably outside the asymptotic range; treat the extrapolation with caution.`);
    const asym = osc ? NaN : g.asymptotic, ls = rows.length > 3 ? lsqOrder(h, f) : null, U = Math.abs(gf * f[0]);
    if (Number.isFinite(asym) && Math.abs(asym - 1) > 0.1) warnings.push('The asymptotic-range check differs from 1 by more than 10%: refine further before relying on the error band.');
    if (Math.max(h[1] / h[0], h[2] / h[1]) < 1.3) warnings.push('Refinement ratios below 1.3 make the observed order sensitive to noise.');
    const kp = [
      { key: 'observed_order', label: 'Observed order of accuracy', value: p, unit: '-', status: osc ? 'warn' : Math.abs(p - i.p_formal) < 0.5 * i.p_formal ? 'ok' : 'warn', note: osc ? 'Formal order used (oscillatory convergence)' : `Formal order ${i.p_formal}` },
      { key: 'f_extrapolated', label: 'Richardson-extrapolated value', value: fEx, unit: '' }, { key: 'gci_fine_pct', label: 'GCI on the fine grid', value: 100 * gf, unit: '%', status: gf < 0.01 ? 'ok' : gf < 0.05 ? 'warn' : 'bad' },
      { key: 'gci_coarse_pct', label: 'GCI on the medium grid', value: 100 * gc, unit: '%' }, { key: 'numerical_uncertainty', label: 'Numerical uncertainty of the fine-grid value (±)', value: U, unit: '' },
      { key: 'asymptotic_ratio', label: 'Asymptotic-range check', value: asym, unit: '-', status: Math.abs(asym - 1) < 0.1 ? 'ok' : 'warn', note: 'GCI₂₃ / (rᵖ·GCI₁₂); 1 means asymptotic' },
      { key: 'monotonic', label: 'Monotonic convergence', value: osc ? 0 : 1, unit: '', status: osc ? 'warn' : 'ok' }, { key: 'fine_grid_error_pct', label: 'Estimated error of the fine-grid value', value: 100 * Math.abs((f[0] - fEx) / (fEx || 1e-300)), unit: '%' },
    ];
    if (osc) kp.splice(5, 1);
    if (ls) kp.push({ key: 'lsq_order', label: 'Least-squares observed order (all grids)', value: ls.p, unit: '-' }, { key: 'lsq_f0', label: 'Least-squares extrapolated value', value: ls.f0, unit: '' }, { key: 'lsq_uncertainty_pct', label: 'Least-squares numerical uncertainty', value: (100 * (1.25 * Math.abs(f[0] - ls.f0) + ls.sd)) / Math.abs(f[0] || 1e-300), unit: '%', note: '1.25·|f₁ − f₀| plus the fit scatter' });
    const hh = N.linspace(0, h[h.length - 1], 40), plots = [{ type: 'line', title: 'Grid convergence with extrapolation and uncertainty band', xlabel: 'Representative cell size h', ylabel: 'Result', series: [{ name: 'Computed', x: h, y: f, style: 'points' }, { name: 'Fit f₀ + α·hᵖ', x: hh, y: hh.map((v) => fEx + ((f[0] - fEx) * v ** p) / h[0] ** p) }, { name: 'Fine grid + GCI', x: [0, h[0]], y: [f[0] + U, f[0] + U], style: 'dash' }, { name: 'Fine grid − GCI', x: [0, h[0]], y: [f[0] - U, f[0] - U], style: 'dash' }], annotations: [{ y: fEx, label: 'Extrapolated' }] }];
    // iterative convergence
    const rs = String(i.residuals || '').split(/[\s,;]+/).filter(Boolean).map(Number), res = rs.filter((v) => Number.isFinite(v) && v > 0);
    if (rs.length !== res.length) warnings.push(`${rs.length - res.length} residual value(s) could not be read as positive numbers.`);
    if (res.length >= 4) {
      const drop = Math.log10(res[0] / res[res.length - 1]), half = res.slice(Math.floor(res.length / 2)), rate = (half[half.length - 1] / half[0]) ** (1 / (half.length - 1)), more = rate < 1 && drop < i.resid_target ? Math.ceil(((i.resid_target - drop) * Math.log(10)) / -Math.log(rate)) : 0, stall = rate > 0.98;
      kp.push({ key: 'resid_drop_decades', label: 'Residual reduction achieved', value: drop, unit: 'decades', status: drop >= i.resid_target ? 'ok' : 'warn', note: `Target ${i.resid_target}` }, { key: 'resid_rate', label: 'Asymptotic convergence factor per iteration', value: rate, unit: '-', status: stall ? 'bad' : 'ok' }, { key: 'iterations_to_target', label: 'Further iterations to reach the target', value: stall && drop < i.resid_target ? Infinity : more, unit: '' },
        { key: 'iter_error_ratio', label: 'Remaining iterative error / last change', value: rate < 1 ? rate / (1 - rate) : Infinity, unit: '-', note: 'Geometric-series estimate: error ≈ ρ/(1 − ρ) × last update' });
      if (stall) warnings.push('The residual has stalled (convergence factor above 0.98): iterative error may contaminate the grid study.');
      plots.push({ type: 'line', title: 'Iterative residual history', xlabel: 'Iteration [-]', ylabel: 'Residual [-]', ylog: true, series: [{ name: 'Residual', x: res.map((_, k) => k), y: res }], annotations: [{ y: res[0] * 10 ** -i.resid_target, label: 'Target' }] });
    }
    return { kpis: kp, plots, tables: [{ title: 'Grid study', columns: ['Grid', 'h', 'Result', 'Change from next finer', 'Refinement ratio'], rows: rows.map((r, k) => [k + 1, r[0], r[1], k ? r[1] - rows[k - 1][1] : 0, k ? r[0] / rows[k - 1][0] : 1]) }],
      warnings, models: ['Generalised Richardson extrapolation with observed order (non-uniform refinement ratios allowed)', 'Grid Convergence Index (Roache)', 'Least-squares order fit for more than three grids', 'Geometric iterative-convergence estimate'],
      assumptions: ['The three finest grids are used for the GCI', 'Grids are geometrically similar and in the asymptotic range', 'Iterative and round-off errors are much smaller than the discretisation error', 'The built-in default numbers are a synthetic illustration'] };
  },
  verify() {
    const o = N.kv(gciCalc.run({ grids: '1, 0.97050; 2, 0.96854; 4, 0.96178', Fs: 1.25, p_formal: 2, residuals: '', resid_target: 5 }));
    const ex = (h) => 3 + 0.2 * h ** 1.7, q = N.kv(gciCalc.run({ grids: [0.5, 1, 1.5, 3, 4].map((h) => `${h}, ${ex(h)}`).join('\n'), Fs: 1.25, p_formal: 2, residuals: '1,0.5,0.25,0.125,0.0625', resid_target: 3 }));
    return [
      N.check('Observed order ln(ε₃₂/ε₂₁)/ln r', o.observed_order, Math.log(0.00676 / 0.00196) / Math.log(2), 1e-9, 'Richardson extrapolation, constant ratio r = 2'),
      N.check('Extrapolated value f₁ + (f₁ − f₂)/(rᵖ − 1)', o.f_extrapolated, 0.9705 + 0.00196 / (0.00676 / 0.00196 - 1), 1e-10, 'Richardson extrapolation'),
      N.check('Fine-grid GCI = Fs·|ε|/(rᵖ − 1)', o.gci_fine_pct, (125 * (0.00196 / 0.9705)) / (0.00676 / 0.00196 - 1), 1e-9, 'Roache (1994)'),
      N.check('Non-uniform ratios: exact order recovered', q.observed_order, 1.7, 1e-6, 'Manufactured f = 3 + 0.2·h^1.7'),
      N.check('Least-squares fit recovers the grid-independent value', q.lsq_f0, 3, 1e-6, 'Manufactured data'),
      N.check('Geometric residual history: factor 0.5', q.resid_rate, 0.5, 1e-12, 'Definition'),
    ];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (!o.monotonic) out.push({ severity: 'warn', title: 'Convergence is not monotonic', detail: 'The result oscillates between grids.', action: 'Add a finer grid, check iterative convergence on each grid and make sure the grids are refined uniformly; the reported band uses the conservative safety factor of 3.', basis: 'Grid Convergence Index procedure' });
    if (o.gci_fine_pct > 5) out.push({ severity: 'warn', title: 'Numerical uncertainty is large', detail: `GCI ${o.gci_fine_pct.toFixed(1)}% on the fine grid.`, action: 'Refine further: the numerical error would dominate any comparison with test data.', basis: 'Fine-grid GCI' });
    else if (o.gci_fine_pct < 0.5 && o.gci_coarse_pct < 2) out.push({ severity: 'advise', title: 'The medium grid is probably sufficient', detail: `GCI is ${o.gci_coarse_pct.toFixed(2)}% on the medium grid and ${o.gci_fine_pct.toFixed(2)}% on the fine grid.`, action: 'Run production cases on the medium grid and carry its GCI as the numerical uncertainty: similar accuracy for a fraction of the computing time and energy.', basis: 'Mesh selection from the GCI' });
    if (o.resid_drop_decades !== undefined && o.resid_drop_decades < 3) out.push({ severity: 'warn', title: 'Iterative convergence is shallow', detail: `Residuals fell ${o.resid_drop_decades.toFixed(1)} decades.`, action: `Continue for about ${Number.isFinite(o.iterations_to_target) ? o.iterations_to_target : 'many'} more iterations, or improve the solver settings, before using the result in a grid study.`, basis: 'Iterative error must be well below discretisation error' });
    out.push({ severity: 'info', title: 'Carry the uncertainty forward', detail: `Fine-grid value ± ${Number(o.numerical_uncertainty).toPrecision(3)} (numerical).`, action: 'Enter this as the numerical uncertainty in the validation analysis so that model error is judged against what the data and the numerics can resolve.', basis: 'ASME V&V 20 approach' });
    return out;
  },
};

// synthetic illustration: lift coefficient against angle of attack, "prediction" versus "measurement" with measurement uncertainty
const VAL_DEFAULT = '0, 0.210, 0.203, 0.012\n2, 0.425, 0.431, 0.012\n4, 0.640, 0.652, 0.012\n6, 0.852, 0.861, 0.014\n8, 1.058, 1.049, 0.016\n10, 1.251, 1.214, 0.018\n12, 1.415, 1.342, 0.022\n14, 1.520, 1.401, 0.030';
/** Area between the empirical CDFs of two samples (Minkowski L1 / area validation metric). */
function areaMetric(a, b) {
  const A = a.slice().sort((x, y) => x - y), B = b.slice().sort((x, y) => x - y), pts = [...new Set([...A, ...B])].sort((x, y) => x - y); let area = 0, ia = 0, ib = 0;
  for (let k = 0; k < pts.length - 1; k++) { while (ia < A.length && A[ia] <= pts[k]) ia++; while (ib < B.length && B[ib] <= pts[k]) ib++; area += Math.abs(ia / A.length - ib / B.length) * (pts[k + 1] - pts[k]); }
  return area;
}
const validation = {
  id: 'validation', title: 'Validation metrics and validation uncertainty', fidelity: 'analytical',
  summary: 'Compares predictions with measurements: error statistics, the area between their distributions, and a point-by-point test of whether the difference is larger than the combined numerical, input and measurement uncertainty. Ends with a pass/fail against your accuracy requirement.',
  equations: ['Error propagation equations', 'Statistical confidence interval equations', 'Least-squares calibration equations'],
  inputs: [
    { key: 'data', label: 'Data: x, prediction, observation[, measurement uncertainty]', type: 'text', default: VAL_DEFAULT, group: 'Data', help: 'One point per line. The fourth column is the standard uncertainty of the measurement. The default dataset is SYNTHETIC and only illustrates the method.' },
    { key: 'u_num_pct', label: 'Numerical standard uncertainty', unit: '% of prediction', default: 0.5, min: 0, max: 50, group: 'Uncertainties', help: 'About GCI/1.15 from the grid-convergence calculator' },
    { key: 'u_input_pct', label: 'Input-parameter standard uncertainty', unit: '% of prediction', default: 1, min: 0, max: 50, group: 'Uncertainties', help: 'From an uncertainty-propagation study of the simulation inputs' },
    { key: 'u_d_pct', label: 'Measurement standard uncertainty when not in the data', unit: '% of observation', default: 2, min: 0, max: 50, group: 'Uncertainties' },
    { key: 'k_cov', label: 'Coverage factor', unit: '-', default: 2, min: 1, max: 3, group: 'Uncertainties', help: '2 for about 95%' },
    { key: 'tol_pct', label: 'Accuracy requirement', unit: '% of observation', default: 5, min: 0.1, max: 100, group: 'Acceptance' },
    { key: 'share_req', label: 'Share of points that must meet it', unit: '-', default: 0.9, min: 0.5, max: 1, group: 'Acceptance' },
  ],
  defaults: (c, up) => ({ u_num_pct: up.vvuq?.gci_fine_pct > 0 && up.vvuq.gci_fine_pct < 50 ? up.vvuq.gci_fine_pct / 1.15 : undefined }),
  run(i) {
    const warnings = [], t = parseTable(i.data, 3); let rows = t.rows, synthetic = String(i.data).trim() === VAL_DEFAULT;
    if (t.bad) warnings.push(`${t.bad} data line(s) could not be read and were skipped.`);
    if (rows.length < 3) { rows = parseTable(VAL_DEFAULT, 3).rows; synthetic = true; warnings.push('Fewer than three valid data points; the built-in synthetic dataset was used.'); }
    if (synthetic) warnings.push('The dataset is synthetic. Results illustrate the method and say nothing about any real model.');
    rows = rows.slice().sort((a, b) => a[0] - b[0]);
    const x = rows.map((r) => r[0]), S = rows.map((r) => r[1]), D = rows.map((r) => r[2]), uD = rows.map((r) => (Number.isFinite(r[3]) && r[3] >= 0 ? r[3] : (i.u_d_pct / 100) * Math.abs(r[2]))), n = rows.length, m = N.errorMetrics(S, D);
    const E = S.map((s, k) => s - D[k]), uv = S.map((s, k) => Math.hypot((i.u_num_pct / 100) * s, (i.u_input_pct / 100) * s, uD[k])), U = uv.map((v) => i.k_cov * v), within = E.filter((e, k) => Math.abs(e) <= U[k]).length, tolOk = E.filter((e, k) => Math.abs(e) <= (i.tol_pct / 100) * Math.abs(D[k])).length;
    const sdE = N.std(E), ci = n > 1 ? (1.96 * sdE) / Math.sqrt(n) : NaN, area = areaMetric(S, D), accept = tolOk / n >= i.share_req, worst = N.argmax(E.map((e, k) => Math.abs(e) / Math.max(Math.abs(D[k]), 1e-300))), resolvable = E.filter((e, k) => U[k] <= (i.tol_pct / 100) * Math.abs(D[k])).length;
    if (resolvable < n) warnings.push(`At ${n - resolvable} point(s) the validation uncertainty is larger than the accuracy requirement: the data cannot confirm that accuracy there, whatever the model does.`);
    if (Math.abs(m.bias) > 2 * ci && n >= 5) warnings.push('The mean error differs from zero by more than its 95% confidence interval: the model has a systematic bias over this range.');
    if (n < 8) warnings.push('Fewer than eight points: statistics such as R² and the confidence interval are indicative only.');
    return {
      kpis: [
        { key: 'validation_rmse', label: 'Root-mean-square error', value: m.rmse, unit: '' }, { key: 'validation_mae', label: 'Mean absolute error', value: m.mae, unit: '' },
        { key: 'validation_bias', label: 'Mean error (bias)', value: m.bias, unit: '', note: `95% confidence ± ${Number(ci).toPrecision(2)}` }, { key: 'validation_mape_pct', label: 'Mean absolute percentage error', value: m.mape, unit: '%' },
        { key: 'validation_r2', label: 'Coefficient of determination R²', value: m.r2, unit: '-', status: m.r2 > 0.9 ? 'ok' : 'warn' }, { key: 'area_metric', label: 'Area validation metric', value: area, unit: '', note: 'Area between the two empirical distributions; same unit as the data' },
        { key: 'mean_E_over_uval', label: 'Mean |E| / u_val', value: N.mean(E.map((e, k) => Math.abs(e) / Math.max(uv[k], 1e-300))), unit: '-', note: 'Above the coverage factor, model error is detectable' },
        { key: 'share_within_uval_pct', label: 'Points where |E| ≤ k·u_val', value: (100 * within) / n, unit: '%', note: 'Agreement to within what can be resolved' },
        { key: 'share_within_tol_pct', label: 'Points meeting the accuracy requirement', value: (100 * tolOk) / n, unit: '%', status: accept ? 'ok' : 'bad' },
        { key: 'validation_accepted', label: 'Accuracy requirement met', value: accept ? 1 : 0, unit: '', status: accept ? 'ok' : 'bad', note: `${(100 * i.share_req).toFixed(0)}% of points within ${i.tol_pct}%` },
        { key: 'max_error_pct', label: 'Largest relative error', value: (100 * Math.abs(E[worst])) / Math.max(Math.abs(D[worst]), 1e-300), unit: '%', note: `At x = ${x[worst]}` },
      ],
      plots: [
        { type: 'line', title: 'Prediction and observation', xlabel: 'x', ylabel: 'Quantity of interest', series: [{ name: 'Prediction', x, y: S, style: 'line+points' }, { name: 'Observation', x, y: D, style: 'points' }, { name: 'Observation + k·u_D', x, y: D.map((d, k) => d + i.k_cov * uD[k]), style: 'dash' }, { name: 'Observation − k·u_D', x, y: D.map((d, k) => d - i.k_cov * uD[k]), style: 'dash' }] },
        { type: 'line', title: 'Comparison error against validation uncertainty', xlabel: 'x', ylabel: 'E = prediction − observation', series: [{ name: 'Comparison error E', x, y: E, style: 'line+points' }, { name: '+k·u_val', x, y: U, style: 'dash' }, { name: '−k·u_val', x, y: U.map((v) => -v), style: 'dash' }], annotations: [{ y: 0, label: 'Perfect agreement' }] },
        { type: 'line', title: 'Parity plot', xlabel: 'Observation', ylabel: 'Prediction', series: [{ name: 'Data', x: D, y: S, style: 'points' }, { name: 'Perfect', x: [N.amin(D), N.amax(D)], y: [N.amin(D), N.amax(D)], style: 'dash' }] },
        { type: 'line', title: 'Empirical distributions (area metric)', xlabel: 'Quantity of interest', ylabel: 'Cumulative probability [-]', series: [{ name: 'Prediction', x: S.slice().sort((a, b) => a - b), y: N.range(n, (k) => (k + 1) / n), style: 'step' }, { name: 'Observation', x: D.slice().sort((a, b) => a - b), y: N.range(n, (k) => (k + 1) / n), style: 'step' }] },
      ],
      tables: [{ title: 'Point-by-point assessment', columns: ['x', 'Prediction S', 'Observation D', 'E = S − D', 'E [% of D]', 'u_val', '|E| ≤ k·u_val', 'Model error interval low', 'Model error interval high'], rows: rows.map((_, k) => [x[k], S[k], D[k], E[k], (100 * E[k]) / (D[k] || 1e-300), uv[k], Math.abs(E[k]) <= U[k] ? 'yes' : 'no', E[k] - U[k], E[k] + U[k]]) }],
      outputs: { n_points: n, synthetic_data: synthetic ? 1 : 0 },
      warnings, models: ['Error statistics (RMSE, MAE, bias, MAPE, R²)', 'Area validation metric between empirical CDFs', 'Validation comparison error E and validation uncertainty u_val = √(u_num² + u_input² + u_D²)', 'Acceptance against a stated accuracy requirement'],
      assumptions: ['Numerical, input and measurement uncertainties are independent', 'Predictions and observations are paired at the same conditions', 'Model error lies in E ± k·u_val; when |E| ≤ k·u_val it cannot be distinguished from zero', 'Validation holds only inside the range of conditions covered by the data'],
    };
  },
  verify() {
    const o = N.kv(validation.run({ data: '0,1,2\n1,3,2\n2,2,4\n3,6,4', u_num_pct: 0, u_input_pct: 0, u_d_pct: 0, k_cov: 2, tol_pct: 5, share_req: 0.9 }));
    return [
      N.check('RMSE of errors (−1, 1, −2, 2)', o.validation_rmse, Math.sqrt(10 / 4), 1e-12, 'Definition'),
      N.check('Bias of errors (−1, 1, −2, 2)', o.validation_bias + 1, 1, 1e-12, 'Definition'),
      N.check('Area metric of a constant shift equals the shift', areaMetric([1, 2, 3, 4], [1.5, 2.5, 3.5, 4.5]), 0.5, 1e-12, 'Area between translated CDFs'),
      N.check('Area metric of identical samples is zero', areaMetric([3, 1, 2], [1, 2, 3]) + 1, 1, 1e-12, 'Identity'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.synthetic_data) out.push({ severity: 'info', title: 'Replace the synthetic data', detail: 'The built-in dataset only demonstrates the calculations.', action: 'Paste paired predictions and measurements (wind tunnel, structural test, flight test) with their measurement uncertainties.', basis: 'Validation requires independent experimental evidence' });
    // the built-in dataset is an illustration with a deliberate model error near its upper end: it cannot raise a warning about the user's model
    if (!o.validation_accepted) out.push({ severity: o.synthetic_data ? 'info' : 'warn', title: o.synthetic_data ? 'The synthetic example misses its accuracy requirement (as designed)' : 'The accuracy requirement is not met', detail: `${o.share_within_tol_pct.toFixed(0)}% of points are within ${i.tol_pct}%; the largest error is ${o.max_error_pct.toFixed(1)}%.`, action: o.share_within_uval_pct < 60 ? 'The errors exceed what the uncertainties can explain: the model form is deficient in part of the range. Restrict the validated domain, improve the model, or calibrate on separate data and re-validate.' : 'Errors are mostly inside the validation uncertainty: tighter measurements or finer grids are needed before the model can be shown to meet this requirement.', basis: 'Comparison error against validation uncertainty and requirement' });
    else out.push({ severity: 'info', title: 'Requirement met over the tested range', detail: `RMSE ${o.validation_rmse.toPrecision(3)}, bias ${o.validation_bias.toPrecision(2)}.`, action: 'Record the range of conditions covered. Do not extrapolate the validation outside it.', basis: 'Validation domain' });
    if (o.mean_E_over_uval > i.k_cov) out.push({ severity: 'advise', title: 'Model error is detectable', detail: `|E| averages ${o.mean_E_over_uval.toFixed(1)} validation uncertainties.`, action: 'Use the Bayesian calibration or the multifidelity discrepancy analysis to quantify and correct it, keeping a separate dataset for the final validation.', basis: '|E| > k·u_val implies model-form error' });
    return out;
  },
};

const propagation = {
  id: 'propagation', title: 'Uncertainty propagation: Monte Carlo, Latin hypercube, polynomial chaos and collocation', fidelity: 'numerical',
  summary: 'Pushes input uncertainty through a model by four methods and compares the mean and spread they give and how many model runs each needs.',
  equations: ['Error propagation equations', 'Polynomial chaos expansion equations', 'Stochastic collocation equations', 'Statistical confidence interval equations'],
  inputs: [...MODEL_IN,
    { key: 'nSamples', label: 'Monte Carlo and Latin hypercube samples', unit: '', default: 2000, min: 64, max: 200000, step: 1, discrete: true, group: 'Numerics' },
    { key: 'pce_order', label: 'Polynomial chaos order', unit: '', default: 3, min: 1, max: 8, step: 1, discrete: true, group: 'Numerics', help: 'Raise to 6–8 for the Ishigami function' },
    { key: 'pce_ratio', label: 'Regression samples per polynomial term', unit: '-', default: 3, min: 1.2, max: 10, group: 'Numerics' },
    { key: 'colloc_pts', label: 'Collocation points per dimension', unit: '', default: 4, min: 2, max: 12, step: 1, discrete: true, group: 'Numerics' },
    { key: 'seed', label: 'Random seed', unit: '', default: 25, min: 1, max: 1e9, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: modelDefaults,
  run(i) {
    const M = testModel(i), n = Math.round(i.nSamples), u = N.rng(i.seed), warnings = [], st = (a) => ({ mean: N.mean(a), sd: N.std(a) });
    const mcX = N.range(n, () => N.range(M.d, () => toStd(M.kind, u()))), mcY = mcX.map(M.f), lhY = N.lhs(n, M.d, u).map((r) => M.f(r.map((p) => toStd(M.kind, p)))), mc = st(mcY), lh = st(lhY);
    const p = Math.round(i.pce_order), terms = multiIndex(M.d, p).length, nP = Math.ceil(i.pce_ratio * terms), pX = N.lhs(nP, M.d, u).map((r) => r.map((q) => toStd(M.kind, q))), pce = pceFit(pX, pX.map(M.f), p, M.kind), sc = collocate(M.f, M.d, Math.round(i.colloc_pts), M.kind);
    const ref = M.exact || (() => { const r = collocate(M.f, M.d, 12, M.kind); return { mean: r.mean, sd: r.sd }; })(), err = (v, r) => 100 * Math.abs(v / r - 1);
    // convergence of the sampled mean with sample count (prefixes of the same streams; LHS re-drawn per size)
    const ns = [50, 100, 200, 400, 800, 1600, 3200, 6400].filter((v) => v <= n), eMc = ns.map((k) => Math.abs(N.mean(mcY.slice(0, k)) - ref.mean) / ref.sd), u2 = N.rng(i.seed + 1), eLh = ns.map((k) => Math.abs(N.mean(N.lhs(k, M.d, u2).map((r) => M.f(r.map((q) => toStd(M.kind, q))))) - ref.mean) / ref.sd);
    const val = mcX.slice(0, 200), pErr = rms(val.map((x, k) => pce.predict(x) - mcY[k])) / (ref.sd || 1), hist = N.histogram(mcY, 30);
    if (pErr > 0.1) warnings.push(`The polynomial chaos surrogate misses ${(100 * pErr).toFixed(0)}% of the output standard deviation on independent samples: raise the order (the model is strongly non-polynomial) or use sampling.`);
    if (!M.exact) warnings.push('No closed-form statistics exist for this model: the reference is a 12-point-per-dimension Gauss quadrature.');
    if (i.model === 'Breguet range') warnings.push('Input distributions are assumed independent normals with the coefficients of variation you set; correlated inputs need a transformation first.');
    const rows = [['Monte Carlo', n, mc.mean, mc.sd], ['Latin hypercube', n, lh.mean, lh.sd], [`Polynomial chaos (order ${p}, regression)`, nP, pce.mean, pce.sd], [`Stochastic collocation (${Math.round(i.colloc_pts)} points/dim)`, sc.n, sc.mean, sc.sd], [M.exact ? 'Exact' : 'Reference quadrature', M.exact ? 0 : 12 ** M.d, ref.mean, ref.sd]];
    return {
      kpis: [
        { key: 'mean_ref', label: 'Reference mean', value: ref.mean, unit: M.unit }, { key: 'sd_ref', label: 'Reference standard deviation', value: ref.sd, unit: M.unit }, { key: 'cov_out_pct', label: 'Output coefficient of variation', value: (100 * ref.sd) / Math.abs(ref.mean || 1), unit: '%' },
        { key: 'mc_mean', label: 'Monte Carlo mean', value: mc.mean, unit: M.unit, note: `± ${(1.96 * mc.sd / Math.sqrt(n)).toPrecision(2)} (95%)` }, { key: 'mc_sd_err_pct', label: 'Monte Carlo error in standard deviation', value: err(mc.sd, ref.sd), unit: '%' },
        { key: 'lhs_mean_err_sd', label: 'Latin hypercube mean error / σ', value: Math.abs(lh.mean - ref.mean) / ref.sd, unit: '-' }, { key: 'mc_mean_err_sd', label: 'Monte Carlo mean error / σ', value: Math.abs(mc.mean - ref.mean) / ref.sd, unit: '-' },
        { key: 'pce_mean', label: 'Polynomial chaos mean', value: pce.mean, unit: M.unit }, { key: 'pce_sd_err_pct', label: 'Polynomial chaos error in standard deviation', value: err(pce.sd, ref.sd), unit: '%', status: err(pce.sd, ref.sd) < 5 ? 'ok' : 'warn' },
        { key: 'pce_runs', label: 'Model runs for polynomial chaos', value: nP, unit: '' }, { key: 'pce_surrogate_err', label: 'Polynomial chaos surrogate error / σ', value: pErr, unit: '-', status: pErr < 0.1 ? 'ok' : 'warn' },
        { key: 'sc_sd_err_pct', label: 'Collocation error in standard deviation', value: err(sc.sd, ref.sd), unit: '%', status: err(sc.sd, ref.sd) < 5 ? 'ok' : 'warn' }, { key: 'sc_runs', label: 'Model runs for collocation', value: sc.n, unit: '' },
        { key: 'p05', label: '5th percentile (Monte Carlo)', value: N.quantile(mcY, 0.05), unit: M.unit }, { key: 'p95', label: '95th percentile (Monte Carlo)', value: N.quantile(mcY, 0.95), unit: M.unit },
      ],
      plots: [
        { type: 'bar', title: 'Output distribution (Monte Carlo)', ylabel: 'Samples [-]', categories: hist.centers.map((v) => v.toPrecision(4)), series: [{ name: 'Count', y: hist.counts }] },
        { type: 'line', title: 'Error of the sampled mean versus sample count', xlabel: 'Samples [-]', ylabel: '|mean error| / σ [-]', xlog: true, ylog: true, series: [{ name: 'Monte Carlo', x: ns, y: eMc.map((v) => Math.max(v, 1e-8)), style: 'line+points' }, { name: 'Latin hypercube', x: ns, y: eLh.map((v) => Math.max(v, 1e-8)), style: 'line+points' }, { name: '1/√N', x: ns, y: ns.map((k) => 1 / Math.sqrt(k)), style: 'dash' }] },
        { type: 'bar', title: 'Sensitivity indices from the polynomial chaos coefficients', ylabel: 'Share of variance [-]', categories: M.names, series: [{ name: 'First-order', y: pce.S1 }, { name: 'Total', y: pce.ST }] },
      ],
      tables: [{ title: 'Method comparison', columns: ['Method', 'Model runs', 'Mean', 'Standard deviation', 'Mean error [% of σ]', 'σ error [%]'], rows: rows.map((r) => [...r, (100 * Math.abs(r[2] - ref.mean)) / ref.sd, err(r[3], ref.sd)]) }],
      warnings, models: ['Monte Carlo sampling (seeded)', 'Latin hypercube sampling', `Polynomial chaos expansion: ${M.kind === 'normal' ? 'Hermite' : 'Legendre'} basis, total degree, least-squares regression`, 'Stochastic collocation on a tensor Gauss grid (Golub–Welsch nodes)'],
      assumptions: ['Inputs are independent', M.kind === 'normal' ? 'Inputs are normal' : 'Inputs are uniform on [−π, π]', 'Tensor grids grow as (points)^dimensions: practical only for a handful of uncertain inputs'],
    };
  },
  convergence: { param: 'nSamples', label: 'Samples', levels: [250, 500, 1000, 2000, 4000], metric: 'mc_mean', hOf: (n) => 1 / Math.sqrt(n) },
  verify() {
    const gh = gaussRule(5, 'normal'), gl = gaussRule(4, 'uniform'), b = { ...Object.fromEntries(propagation.inputs.map((f) => [f.key, f.default])) };
    const pl = N.kv(propagation.run({ ...b, model: 'Polynomial test', pce_order: 2 })), ish = N.kv(propagation.run({ ...b, model: 'Ishigami function', pce_order: 8, colloc_pts: 10, nSamples: 20000 }));
    return [
      N.check('Gauss–Hermite integrates the fourth normal moment E[ξ⁴] = 3', N.sum(gh.x.map((x, k) => gh.w[k] * x ** 4)), 3, 1e-10, 'Moments of the standard normal'),
      N.check('Gauss–Legendre integrates E[ξ⁶] = 1/7 on [−1, 1]', N.sum(gl.x.map((x, k) => gl.w[k] * x ** 6)), 1 / 7, 1e-10, 'Uniform moments'),
      N.check('Polynomial chaos is exact for a polynomial: mean of ξ₁² + ξ₁ξ₂', pl.pce_mean, 1, 1e-8, 'Exact projection'),
      N.check('Polynomial chaos variance of ξ₁² + ξ₁ξ₂ is 3', pl.pce_sd_err_pct + 1, 1, 1e-6, 'Var = 2 + 1'),
      N.check('Ishigami mean a/2', ish.mean_ref, 3.5, 1e-12, 'Analytical (a = 7, b = 0.1)'),
      N.check('Monte Carlo reproduces the Ishigami standard deviation', ish.mc_sd_err_pct + 100, 100, 0.02, 'Analytical variance 13.8446; sampling error ≈ 1%'),
      N.check('Collocation with 10 points reproduces the Ishigami standard deviation', ish.sc_sd_err_pct + 100, 100, 1e-3, 'Gauss–Legendre exactness for smooth integrands'),
      N.check('Order-8 polynomial chaos reproduces the Ishigami standard deviation', ish.pce_sd_err_pct + 100, 100, 0.03, 'Spectral convergence'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    out.push({ severity: 'info', title: `Output uncertainty: ${o.cov_out_pct.toFixed(1)}% coefficient of variation`, detail: `90% interval ${o.p05.toPrecision(4)} … ${o.p95.toPrecision(4)}.`, action: 'Quote results with this interval, and design to the unfavourable percentile rather than the mean where a shortfall matters (range, payload, margins).', basis: 'Propagated input uncertainty' });
    if (o.pce_sd_err_pct < 2 && o.pce_surrogate_err < 0.1) out.push({ severity: 'advise', title: 'Polynomial chaos is the efficient choice here', detail: `${o.pce_runs} runs match the reference spread to ${o.pce_sd_err_pct.toFixed(2)}%, against ${Math.round(i.nSamples)} for sampling.`, action: 'For expensive solvers (CFD, FE) use a low-order polynomial chaos or collocation with a few inputs: the same statistics for a small fraction of the computing time and energy.', basis: 'Method comparison' });
    else out.push({ severity: 'advise', title: 'Use sampling for this model', detail: `Polynomial chaos at order ${Math.round(i.pce_order)} is off by ${o.pce_sd_err_pct.toFixed(1)}% in standard deviation.`, action: 'Raise the order if the model is smooth, otherwise prefer Latin hypercube sampling, which converges faster than plain Monte Carlo for the mean.', basis: 'Surrogate error on independent samples' });
    return out;
  },
};

/** Saltelli first-order and Jansen total-effect estimators on standard-variable samples A, B. */
function sobolIndices(f, A, B) {
  const n = A.length, d = A[0].length, fA = A.map(f), fB = B.map(f), all = fA.concat(fB), V = N.variance(all), S1 = [], ST = [];
  for (let j = 0; j < d; j++) { let s1 = 0, st = 0; for (let k = 0; k < n; k++) { const x = A[k].slice(); x[j] = B[k][j]; const fab = f(x); s1 += fB[k] * (fab - fA[k]); st += (fA[k] - fab) ** 2; } S1.push(s1 / n / V); ST.push(st / (2 * n) / V); }
  return { S1, ST, V, mean: N.mean(all), evals: n * (d + 2) };
}
/** Morris elementary-effects screening on the unit cube with p = 4 levels. */
function morris(f01, d, r, u) {
  const delta = 2 / 3, ee = N.range(d, () => []);
  for (let t = 0; t < r; t++) { let x = N.range(d, () => Math.floor(u() * 2) / 3), y = f01(x); const order = N.range(d).sort(() => u() - 0.5); for (const j of order) { const xn = x.slice(); xn[j] = x[j] + delta; const yn = f01(xn); ee[j].push((yn - y) / delta); x = xn; y = yn; } }
  return { muStar: ee.map((e) => N.mean(e.map(Math.abs))), sigma: ee.map(N.std) };
}
const sobol = {
  id: 'sobol', title: 'Global sensitivity: Sobol indices and Morris screening', fidelity: 'numerical',
  summary: 'Splits the variance of a model output into the share caused by each uncertain input alone and with its interactions, so effort goes to the inputs that matter; a cheap Morris screening is shown alongside.',
  equations: ['Sensitivity derivative equations', 'Error propagation equations', 'Statistical confidence interval equations'],
  inputs: [...MODEL_IN,
    { key: 'nBase', label: 'Base samples N (model runs = N·(d + 2))', unit: '', default: 4096, min: 64, max: 262144, step: 1, discrete: true, group: 'Numerics' },
    { key: 'sampler', label: 'Sampling', type: 'select', options: ['Halton quasi-random', 'pseudo-random'], default: 'Halton quasi-random', group: 'Numerics' },
    { key: 'morris_r', label: 'Morris trajectories', unit: '', default: 40, min: 4, max: 2000, step: 1, discrete: true, group: 'Numerics' },
    { key: 'seed', label: 'Random seed', unit: '', default: 5, min: 1, max: 1e9, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: modelDefaults,
  run(i) {
    const M = testModel(i), n = Math.round(i.nBase), u = N.rng(i.seed), d = M.d, warnings = [];
    const U = i.sampler.startsWith('Halton') ? halton(n, 2 * d, u) : N.range(n, () => N.range(2 * d, () => u())), A = U.map((r) => r.slice(0, d).map((p) => toStd(M.kind, p))), B = U.map((r) => r.slice(d).map((p) => toStd(M.kind, p))), s = sobolIndices(M.f, A, B);
    const f01 = (x) => M.f(x.map((p) => toStd(M.kind, N.clamp(0.02 + 0.96 * p, 1e-6, 1 - 1e-6)))), mo = morris(f01, d, Math.round(i.morris_r), u), inter = s.ST.map((v, j) => v - s.S1[j]), sum1 = N.sum(s.S1);
    const ns = [64, 128, 256, 512, 1024, 2048, 4096, 8192].filter((v) => v <= n), conv = ns.map((k) => sobolIndices(M.f, A.slice(0, k), B.slice(0, k)).S1), ex = M.exact, maxErr = ex ? Math.max(...s.S1.map((v, j) => Math.abs(v - ex.S1[j])), ...s.ST.map((v, j) => Math.abs(v - ex.ST[j]))) : NaN;
    if (s.S1.some((v) => v < -0.03) || sum1 > 1.05) warnings.push('Some indices are outside [0, 1] by more than sampling noise: increase the base sample count.');
    if (1 - sum1 > 0.1) warnings.push(`First-order indices sum to ${sum1.toFixed(2)}: about ${(100 * (1 - sum1)).toFixed(0)}% of the variance comes from interactions, so one-at-a-time studies will mislead.`);
    const rank = N.range(d).sort((a, b) => s.ST[b] - s.ST[a]);
    return {
      kpis: [
        ...M.names.map((nm, j) => ({ key: `S1_${j + 1}`, label: `First-order index: ${nm}`, value: s.S1[j], unit: '-' })), ...M.names.map((nm, j) => ({ key: `ST_${j + 1}`, label: `Total index: ${nm}`, value: s.ST[j], unit: '-' })),
        { key: 'sum_first_order', label: 'Sum of first-order indices', value: sum1, unit: '-', note: '1 for an additive model' }, { key: 'interaction_share', label: 'Variance share from interactions', value: Math.max(0, 1 - sum1), unit: '-' },
        { key: 'output_variance', label: 'Output variance', value: s.V, unit: `${M.unit}²` }, { key: 'model_runs', label: 'Model runs', value: s.evals, unit: '' },
        ...(ex ? [{ key: 'max_index_error', label: 'Largest deviation from the exact indices', value: maxErr, unit: '-', status: maxErr < 0.03 ? 'ok' : 'warn' }] : []),
      ],
      plots: [
        { type: 'bar', title: 'Sobol indices', ylabel: 'Share of output variance [-]', categories: M.names, series: [{ name: 'First-order', y: s.S1 }, { name: 'Total', y: s.ST }, ...(ex ? [{ name: 'Exact first-order', y: ex.S1 }, { name: 'Exact total', y: ex.ST }] : [])] },
        { type: 'line', title: 'Convergence of the first-order indices', xlabel: 'Base samples N [-]', ylabel: 'First-order index [-]', xlog: true, series: M.names.map((nm, j) => ({ name: nm, x: ns, y: conv.map((c) => c[j]), style: 'line+points' })) },
        { type: 'line', title: 'Morris screening: mean absolute effect versus spread', xlabel: 'μ* (overall influence)', ylabel: 'σ (non-linearity or interaction)', series: M.names.map((nm, j) => ({ name: nm, x: [mo.muStar[j]], y: [mo.sigma[j]], style: 'points' })) },
      ],
      tables: [{ title: 'Sensitivity ranking (by total index)', columns: ['Rank', 'Input', 'First-order S', 'Total ST', 'Interaction ST − S', 'Morris μ*', 'Morris σ'], rows: rank.map((j, k) => [k + 1, M.names[j], s.S1[j], s.ST[j], inter[j], mo.muStar[j], mo.sigma[j]]) }],
      warnings, models: ['Sobol variance decomposition', 'Saltelli (2010) first-order estimator and Jansen total-effect estimator', i.sampler.startsWith('Halton') ? 'Randomly shifted Halton quasi-random sampling' : 'Pseudo-random sampling', 'Morris elementary-effects screening (4 levels)'],
      assumptions: ['Inputs are independent', 'Indices describe variance over the stated input distributions, not local derivatives', 'Morris screening covers the central 96% of each input range'],
    };
  },
  convergence: { param: 'nBase', label: 'Base samples', levels: [256, 512, 1024, 2048, 4096], metric: 'S1_1', hOf: (n) => 1 / Math.sqrt(n) },
  verify() {
    const b = Object.fromEntries(sobol.inputs.map((f) => [f.key, f.default])), o = N.kv(sobol.run({ ...b, model: 'Ishigami function', nBase: 16384 }));
    const u = N.rng(3), U = halton(4096, 4, u), lin = sobolIndices((x) => 3 * x[0] + 4 * x[1], U.map((r) => r.slice(0, 2).map((p) => N.normInv(p))), U.map((r) => r.slice(2).map((p) => N.normInv(p))));
    return [
      N.check('Ishigami S₁ = 0.3139', o.S1_1, 0.3139, 0.02, 'Analytical indices for a = 7, b = 0.1'), N.check('Ishigami S₂ = 0.4424', o.S1_2, 0.4424, 0.02, 'Analytical'),
      N.check('Ishigami S₃ = 0', o.S1_3 + 1, 1, 0.01, 'Analytical'), N.check('Ishigami total index ST₁ = 0.5576', o.ST_1, 0.5576, 0.02, 'Analytical'),
      N.check('Ishigami total index ST₃ = 0.2437', o.ST_3, 0.2437, 0.02, 'Analytical'),
      N.check('Linear model 3ξ₁ + 4ξ₂: S₁ = 9/25', lin.S1[0], 0.36, 0.02, 'Additive model: Sᵢ = aᵢ²/Σa²'), N.check('Additive model: total equals first-order', lin.ST[1], 0.64, 0.02, 'No interactions'),
    ];
  },
  recommend(res) {
    const o = res.outputs, rows = res.tables[0].rows, out = [];
    out.push({ severity: 'advise', title: `Reduce uncertainty in "${rows[0][1]}" first`, detail: `It drives ${(100 * rows[0][3]).toFixed(0)}% of the output variance (total index); "${rows[rows.length - 1][1]}" drives ${(100 * rows[rows.length - 1][3]).toFixed(0)}%.`, action: `Spend test and analysis effort on ${rows[0][1].toLowerCase()}. Inputs with a total index below about 0.05 can be fixed at nominal values, which shrinks every later study.`, basis: 'Sobol total-effect indices (factor prioritisation and fixing)' });
    if (o.interaction_share > 0.1) out.push({ severity: 'info', title: 'Interactions matter', detail: `${(100 * o.interaction_share).toFixed(0)}% of the variance is not explained by any input alone.`, action: 'Do not rely on one-at-a-time sensitivity; vary inputs together (as done here) when judging robustness.', basis: 'Sum of first-order indices below 1' });
    return out;
  },
};

// Bayesian calibration ---------------------------------------------------------------------------
const CAL_MODELS = {
  'Drag polar CD = p1 + p2·CL²': { f: (p, x) => p[0] + p[1] * x * x, names: ['CD0', 'k'], linear: (x) => [1, x * x] },
  'Power law y = p1·x^p2': { f: (p, x) => p[0] * Math.max(x, 1e-300) ** p[1], names: ['p1', 'p2'] },
  'Exponential decay y = p1·exp(−p2·x)': { f: (p, x) => p[0] * Math.exp(-p[1] * x), names: ['p1', 'p2'] },
};
// synthetic drag-polar "measurements": CD0 = 0.021, k = 0.045 plus noise of standard deviation 0.0006
const CAL_DEFAULT = '0.10, 0.02168\n0.20, 0.02251\n0.30, 0.02553\n0.40, 0.02795\n0.50, 0.03171\n0.60, 0.03776\n0.70, 0.04262\n0.80, 0.05034\n0.90, 0.05703\n1.00, 0.06629';
const bayesCal = {
  id: 'bayes', title: 'Model calibration: weighted least squares and Bayesian MCMC', fidelity: 'numerical',
  summary: 'Fits up to two parameters of a chosen model to data, first by weighted least squares with parameter covariance, then by Metropolis–Hastings sampling of the posterior to give credible intervals, a predictive band and convergence diagnostics.',
  equations: ['Least-squares calibration equations', 'Maximum likelihood estimation equations', 'Bayesian inference equations', 'Statistical confidence interval equations'],
  inputs: [
    { key: 'model', label: 'Model', type: 'select', options: Object.keys(CAL_MODELS), default: 'Drag polar CD = p1 + p2·CL²', group: 'Model' },
    { key: 'data', label: 'Data: x, y[, standard uncertainty]', type: 'text', default: CAL_DEFAULT, group: 'Data', help: 'One point per line. The default drag-polar dataset is SYNTHETIC (generated from CD0 = 0.021, k = 0.045 with noise).' },
    { key: 'sigma', label: 'Measurement standard uncertainty (0 = estimate from the fit)', unit: 'unit of y', default: 0.0006, min: 0, group: 'Data' },
    { key: 'p1_min', label: 'Prior lower bound, parameter 1', unit: '', default: 0, group: 'Priors (uniform)' }, { key: 'p1_max', label: 'Prior upper bound, parameter 1', unit: '', default: 0.1, group: 'Priors (uniform)' },
    { key: 'p2_min', label: 'Prior lower bound, parameter 2', unit: '', default: 0, group: 'Priors (uniform)' }, { key: 'p2_max', label: 'Prior upper bound, parameter 2', unit: '', default: 0.2, group: 'Priors (uniform)' },
    { key: 'nSamples', label: 'MCMC samples per chain', unit: '', default: 3000, min: 300, max: 200000, step: 1, discrete: true, group: 'Numerics' }, { key: 'nChains', label: 'Chains', unit: '', default: 3, min: 2, max: 8, step: 1, discrete: true, group: 'Numerics' },
    { key: 'seed', label: 'Random seed', unit: '', default: 77, min: 1, max: 1e9, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: () => ({}),
  run(i) {
    const warnings = [], md = CAL_MODELS[i.model] || CAL_MODELS['Drag polar CD = p1 + p2·CL²'], t = parseTable(i.data, 2); let rows = t.rows, synthetic = String(i.data).trim() === CAL_DEFAULT;
    if (t.bad) warnings.push(`${t.bad} data line(s) could not be read and were skipped.`);
    if (rows.length < 3) { rows = parseTable(CAL_DEFAULT, 2).rows; synthetic = true; warnings.push('Fewer than three valid data points; the built-in synthetic dataset was used.'); }
    if (synthetic) warnings.push('The dataset is synthetic; fitted values only illustrate the method.');
    const x = rows.map((r) => r[0]), y = rows.map((r) => r[1]), n = x.length, lo = [Math.min(i.p1_min, i.p1_max), Math.min(i.p2_min, i.p2_max)], hi = [Math.max(i.p1_min, i.p1_max), Math.max(i.p2_min, i.p2_max)], mid = lo.map((v, j) => 0.5 * (v + hi[j]));
    // weighted least squares / maximum likelihood (Gaussian errors)
    const sg0 = rows.map((r) => (Number.isFinite(r[2]) && r[2] > 0 ? r[2] : i.sigma > 0 ? i.sigma : 1)), resid = (p, sg = sg0) => x.map((v, k) => (md.f(p, v) - y[k]) / sg[k]);
    let p0 = mid.map((v, j) => v || (hi[j] - lo[j]) * 0.5 || 1); if (md.linear) p0 = N.lstsq(x.map((v, k) => md.linear(v).map((q) => q / sg0[k])), y.map((v, k) => v / sg0[k]));
    const lm = N.levenbergMarquardt(resid, p0, { maxIter: 200 }), pw = lm.p, dof = Math.max(1, n - 2), s2 = lm.cost / dof, known = i.sigma > 0 || rows.every((r) => Number.isFinite(r[2]) && r[2] > 0), sg = known ? sg0 : sg0.map(() => Math.sqrt(s2));
    // covariance: (JᵀJ)⁻¹ of the weighted residuals for known σ, times the residual variance when σ is estimated from the fit
    let cov = null; try { const J = N.zeros(n, 2); for (let j = 0; j < 2; j++) { const h = 1e-6 * Math.max(1e-3, Math.abs(pw[j])), a = pw.slice(), b = pw.slice(); a[j] += h; b[j] -= h; const ra = resid(a), rb = resid(b); for (let k = 0; k < n; k++) J[k][j] = (ra[k] - rb[k]) / (2 * h); } cov = N.mscale(N.inv(N.matmul(N.transpose(J), J)), known ? 1 : s2); } catch { cov = null; }
    if (!cov || cov.some((r) => r.some((v) => !Number.isFinite(v)))) { cov = [[(1e-2 * (Math.abs(pw[0]) || 1)) ** 2, 0], [0, (1e-2 * (Math.abs(pw[1]) || 1)) ** 2]]; warnings.push('The least-squares covariance is singular: the parameters are not separately identifiable from these data. A fallback proposal was used.'); }
    const se = [Math.sqrt(Math.abs(cov[0][0])), Math.sqrt(Math.abs(cov[1][1]))], rho = cov[0][1] / (se[0] * se[1] || 1), chi2 = known ? lm.cost : dof;
    if (pw.some((v, j) => v < lo[j] || v > hi[j])) warnings.push('The least-squares estimate lies outside the prior bounds: the prior truncates the posterior. Check the bounds for the selected model.');
    // Metropolis–Hastings with a correlated Gaussian proposal scaled by 2.38²/d
    const logPost = (p) => { if (p.some((v, j) => v < lo[j] || v > hi[j])) return -Infinity; let s = 0; for (let k = 0; k < n; k++) s += ((md.f(p, x[k]) - y[k]) / sg[k]) ** 2; return -0.5 * s; };
    let Lc; try { Lc = N.cholesky(cov.map((r, a) => r.map((v, b) => v * (2.38 * 2.38) / 2 + (a === b ? 1e-300 : 0)))); } catch { Lc = [[1.68 * se[0], 0], [0, 1.68 * se[1]]]; }
    const u = N.rng(i.seed), nS = Math.round(i.nSamples), nC = Math.round(i.nChains), burn = Math.floor(0.25 * nS), chains = []; let acc = 0, tot = 0;
    for (let c = 0; c < nC; c++) {
      let p = pw.map((v, j) => N.clamp(v + 2.5 * se[j] * (c % 2 ? 1 : -1) * (1 + 0.3 * c) * (j ? -1 : 1), lo[j], hi[j])), lp = logPost(p); if (!Number.isFinite(lp)) { p = pw.map((v, j) => N.clamp(v, lo[j], hi[j])); lp = logPost(p); } const ch = [];
      for (let s = 0; s < nS; s++) { const z = [N.randn(u), N.randn(u)], q = [p[0] + Lc[0][0] * z[0], p[1] + Lc[1][0] * z[0] + Lc[1][1] * z[1]], lq = logPost(q); tot++; if (Math.log(u() + 1e-300) < lq - lp) { p = q; lp = lq; acc++; } if (s >= burn) ch.push(p); }
      chains.push(ch);
    }
    const all = chains.flat(), post = [0, 1].map((j) => { const a = all.map((p) => p[j]); return { mean: N.mean(a), sd: N.std(a), lo: N.quantile(a, 0.025), hi: N.quantile(a, 0.975), a }; });
    const rhat = [0, 1].map((j) => { const L = chains[0].length, ms = chains.map((c) => N.mean(c.map((p) => p[j]))), W = N.mean(chains.map((c) => N.variance(c.map((p) => p[j])))), Bv = L * N.variance(ms); return W > 0 ? Math.sqrt(((L - 1) / L * W + Bv / L) / W) : 1; });
    const xs = N.linspace(N.amin(x), N.amax(x), 40), thin = all.filter((_, k) => k % Math.max(1, Math.floor(all.length / 400)) === 0), band = xs.map((v) => { const q = thin.map((p) => md.f(p, v)); return [N.quantile(q, 0.025), N.quantile(q, 0.975)]; }), sgm = N.mean(sg);
    const fit = N.errorMetrics(x.map((v) => md.f(post.map((q) => q.mean), v)), y), ar = acc / tot, h1 = N.histogram(post[0].a, 30), h2 = N.histogram(post[1].a, 30);
    if (Math.max(...rhat) > 1.1) warnings.push('Gelman–Rubin statistic above 1.1: the chains have not mixed; run more samples.');
    if (ar < 0.1 || ar > 0.6) warnings.push(`Acceptance rate ${(100 * ar).toFixed(0)}% is outside the efficient 15–50% band; the proposal is poorly scaled for this posterior.`);
    if (Math.abs(rho) > 0.95) warnings.push(`Parameters are ${(100 * Math.abs(rho)).toFixed(0)}% correlated: the data constrain a combination of them much better than each one. Add data over a wider range of x.`);
    if (known && chi2 / dof > 3) warnings.push(`Reduced χ² is ${(chi2 / dof).toFixed(1)}: the misfit is larger than the stated measurement uncertainty explains — model-form error, or the uncertainty is understated.`);
    return {
      kpis: [
        { key: 'p1_mean', label: `Posterior mean ${md.names[0]}`, value: post[0].mean, unit: '' }, { key: 'p1_lo', label: `${md.names[0]} 95% credible lower`, value: post[0].lo, unit: '' }, { key: 'p1_hi', label: `${md.names[0]} 95% credible upper`, value: post[0].hi, unit: '' },
        { key: 'p2_mean', label: `Posterior mean ${md.names[1]}`, value: post[1].mean, unit: '' }, { key: 'p2_lo', label: `${md.names[1]} 95% credible lower`, value: post[1].lo, unit: '' }, { key: 'p2_hi', label: `${md.names[1]} 95% credible upper`, value: post[1].hi, unit: '' },
        { key: 'p1_wls', label: `Least-squares ${md.names[0]}`, value: pw[0], unit: '', note: `± ${se[0].toPrecision(2)} (1σ)` }, { key: 'p2_wls', label: `Least-squares ${md.names[1]}`, value: pw[1], unit: '', note: `± ${se[1].toPrecision(2)} (1σ)` },
        { key: 'p1_se', label: `Standard error ${md.names[0]}`, value: se[0], unit: '' }, { key: 'p2_se', label: `Standard error ${md.names[1]}`, value: se[1], unit: '' },
        { key: 'param_correlation', label: 'Parameter correlation', value: rho, unit: '-', status: Math.abs(rho) < 0.95 ? 'ok' : 'warn', note: 'Identifiability: |ρ| near 1 means the pair is not separately determined' },
        { key: 'reduced_chi2', label: 'Reduced χ²', value: chi2 / dof, unit: '-', status: !known || chi2 / dof < 3 ? 'ok' : 'warn', note: known ? '≈ 1 when the model and the stated uncertainty are consistent' : 'Fixed at 1: σ was estimated from the residuals' },
        { key: 'sigma_used', label: 'Measurement standard uncertainty used', value: sgm, unit: '' }, { key: 'fit_rmse', label: 'RMSE of the calibrated model', value: fit.rmse, unit: '' },
        { key: 'acceptance_rate', label: 'MCMC acceptance rate', value: ar, unit: '-', status: ar > 0.1 && ar < 0.6 ? 'ok' : 'warn' }, { key: 'rhat_max', label: 'Gelman–Rubin R̂ (largest)', value: Math.max(...rhat), unit: '-', status: Math.max(...rhat) < 1.1 ? 'ok' : 'warn' },
      ],
      plots: [
        { type: 'line', title: 'Calibrated model with 95% posterior band', xlabel: 'x', ylabel: 'y', series: [{ name: 'Data', x, y, style: 'points' }, { name: 'Posterior mean model', x: xs, y: xs.map((v) => md.f(post.map((q) => q.mean), v)) }, { name: 'Parameter band low', x: xs, y: band.map((b) => b[0]), style: 'dash' }, { name: 'Parameter band high', x: xs, y: band.map((b) => b[1]), style: 'dash' }, { name: 'Predictive low (with noise)', x: xs, y: band.map((b) => b[0] - 1.96 * sgm), style: 'dash' }, { name: 'Predictive high (with noise)', x: xs, y: band.map((b) => b[1] + 1.96 * sgm), style: 'dash' }] },
        { type: 'bar', title: `Posterior of ${md.names[0]}`, ylabel: 'Samples [-]', categories: h1.centers.map((v) => v.toPrecision(4)), series: [{ name: md.names[0], y: h1.counts }] },
        { type: 'bar', title: `Posterior of ${md.names[1]}`, ylabel: 'Samples [-]', categories: h2.centers.map((v) => v.toPrecision(4)), series: [{ name: md.names[1], y: h2.counts }] },
        { type: 'line', title: 'Chain traces (first parameter)', xlabel: 'Sample after burn-in [-]', ylabel: md.names[0], series: chains.slice(0, 4).map((c, k) => { const st = Math.max(1, Math.floor(c.length / 200)); return { name: `Chain ${k + 1}`, x: c.map((_, j) => j).filter((j) => j % st === 0), y: c.filter((_, j) => j % st === 0).map((p) => p[0]) }; }) },
      ],
      tables: [{ title: 'Parameter estimates', columns: ['Parameter', 'Least squares', 'Standard error', 'Posterior mean', 'Posterior sd', '2.5%', '97.5%', 'R̂'], rows: [0, 1].map((j) => [md.names[j], pw[j], se[j], post[j].mean, post[j].sd, post[j].lo, post[j].hi, rhat[j]]) },
        { title: 'Parameter correlation matrix (least squares)', columns: ['', md.names[0], md.names[1]], rows: [[md.names[0], 1, rho], [md.names[1], rho, 1]] }],
      outputs: { synthetic_data: synthetic ? 1 : 0 },
      warnings, models: ['Weighted least squares by Levenberg–Marquardt (maximum likelihood for Gaussian errors)', 'Parameter covariance (JᵀWJ)⁻¹ and correlation', 'Metropolis–Hastings MCMC with correlated Gaussian proposal, 25% burn-in', 'Gelman–Rubin potential scale reduction'],
      assumptions: ['Independent Gaussian measurement errors', 'Uniform priors between the stated bounds', 'When no uncertainty is given it is estimated from the residuals, which assumes the model form is right', 'Calibration does not validate: keep separate data for validation'],
    };
  },
  convergence: { param: 'nSamples', label: 'MCMC samples per chain', levels: [750, 1500, 3000, 6000], metric: 'p1_mean', hOf: (n) => 1 / Math.sqrt(n) },
  verify() {
    const b = Object.fromEntries(bayesCal.inputs.map((f) => [f.key, f.default])), o = N.kv(bayesCal.run({ ...b, nSamples: 20000 })), r = parseTable(CAL_DEFAULT, 2).rows, A = r.map((q) => [1, q[0] * q[0]]), ex = N.solve(N.matmul(N.transpose(A), A), N.matvec(N.transpose(A), r.map((q) => q[1]))), Ci = N.inv(N.matmul(N.transpose(A), A));
    const pl = N.kv(bayesCal.run({ ...b, model: 'Power law y = p1·x^p2', data: [1, 2, 3, 4, 5].map((x) => `${x}, ${2 * x ** 1.5}`).join('\n'), sigma: 0.01, p1_min: 0, p1_max: 10, p2_min: 0, p2_max: 5, nSamples: 600 }));
    return [
      N.check('Least squares equals the normal-equation solution (CD0)', o.p1_wls, ex[0], 1e-6, 'Linear least squares'), N.check('Least squares equals the normal-equation solution (k)', o.p2_wls, ex[1], 1e-6, 'Linear least squares'),
      N.check('Standard error = σ·√[(AᵀA)⁻¹]₁₁', o.p1_se, 0.0006 * Math.sqrt(Ci[0][0]), 1e-4, 'Gaussian linear model'),
      N.check('Posterior mean equals least squares for a linear Gaussian model with a flat prior', o.p1_mean, o.p1_wls, (0.1 * o.p1_se) / o.p1_wls, 'Conjugate result; MCMC error below 0.1 standard errors'),
      N.check('Posterior standard deviation equals the standard error', (o.p2_hi - o.p2_lo) / 3.92, o.p2_se, 0.06, 'Gaussian posterior'),
      N.check('Nonlinear fit recovers a power-law exponent', pl.p2_wls, 1.5, 1e-6, 'Noise-free manufactured data'),
      N.check('Standard error with known σ does not depend on the misfit (noise-free data)', pl.p2_se, (() => { const xs = [1, 2, 3, 4, 5], Jm = xs.map((x) => [x ** 1.5 / 0.01, (2 * x ** 1.5 * Math.log(x)) / 0.01]); return Math.sqrt(N.inv(N.matmul(N.transpose(Jm), Jm))[1][1]); })(), 1e-5, 'σ·√[(JᵀJ)⁻¹]₂₂ with the analytical Jacobian of p₁·x^p₂'),
    ];
  },
  recommend(res) {
    const o = res.outputs, out = [], t = res.tables[0].rows;
    out.push({ severity: 'info', title: 'Use the interval, not just the best fit', detail: `${t[0][0]} = ${Number(t[0][3]).toPrecision(4)} (95%: ${Number(t[0][5]).toPrecision(4)} … ${Number(t[0][6]).toPrecision(4)}); ${t[1][0]} = ${Number(t[1][3]).toPrecision(4)} (95%: ${Number(t[1][5]).toPrecision(4)} … ${Number(t[1][6]).toPrecision(4)}).`, action: 'Feed the posterior ranges into the uncertainty studies of the suite that uses these parameters (for a drag polar: performance, mission and economics), so fuel and cost predictions carry honest bands.', basis: 'Bayesian posterior' });
    if (Math.abs(o.param_correlation) > 0.95) out.push({ severity: 'advise', title: 'Parameters are poorly identifiable', detail: `Correlation ${o.param_correlation.toFixed(3)}.`, action: 'Collect data over a wider range of the independent variable, or fix one parameter from independent evidence.', basis: 'Parameter correlation matrix' });
    if (o.reduced_chi2 > 3) out.push({ severity: 'warn', title: 'The model does not fit within measurement uncertainty', detail: `Reduced χ² = ${o.reduced_chi2.toFixed(1)}.`, action: 'Look at the residual pattern: a trend means the model form is missing physics; use the discrepancy analysis rather than forcing the parameters to absorb it.', basis: 'Goodness of fit' });
    if (o.rhat_max > 1.1) out.push({ severity: 'warn', title: 'MCMC has not converged', detail: `R̂ = ${o.rhat_max.toFixed(2)}.`, action: 'Increase the number of samples until R̂ is below 1.1 and the traces overlap.', basis: 'Gelman–Rubin diagnostic' });
    return out;
  },
};

// Multifidelity ------------------------------------------------------------------------------------
/** 1-D Gaussian process (zero mean, squared-exponential kernel) with likelihood search over the length scale. */
function gp1(x, y, noise, ell0) {
  const n = x.length, span = N.amax(x) - N.amin(x) || 1, build = (ell, s2) => { const K = N.range(n, (a) => N.range(n, (b) => s2 * Math.exp(-0.5 * ((x[a] - x[b]) / ell) ** 2) + (a === b ? noise * noise + 1e-10 * s2 : 0))); let L; try { L = N.cholesky(K); } catch { return null; } const z = new Array(n); for (let a = 0; a < n; a++) { let s = y[a]; for (let k = 0; k < a; k++) s -= L[a][k] * z[k]; z[a] = s / L[a][a]; } const al = new Array(n); for (let a = n - 1; a >= 0; a--) { let s = z[a]; for (let k = a + 1; k < n; k++) s -= L[k][a] * al[k]; al[a] = s / L[a][a]; } let ld = 0; for (let a = 0; a < n; a++) ld += Math.log(L[a][a]); return { ell, s2, L, al, lnL: -0.5 * N.dot(y, al) - ld }; };
  const s2 = Math.max(N.mean(y.map((v) => v * v)), 1e-30); let best = null; for (const e of ell0 ? [ell0] : N.logspace(0.05 * span, 2 * span, 14)) { const m = build(e, s2); if (m && (!best || m.lnL > best.lnL)) best = m; }
  if (!best) return { predict: () => ({ m: 0, s: Math.sqrt(s2) }), ell: span };
  return { ell: best.ell, predict: (q) => { const k = x.map((v) => best.s2 * Math.exp(-0.5 * ((q - v) / best.ell) ** 2)), w = new Array(n); for (let a = 0; a < n; a++) { let s = k[a]; for (let b = 0; b < a; b++) s -= best.L[a][b] * w[b]; w[a] = s / best.L[a][a]; } return { m: N.dot(k, best.al), s: Math.sqrt(Math.max(0, best.s2 - N.dot(w, w))) }; } };
}
// synthetic illustration: low-fidelity linear lift curve against a few "high-fidelity" points that bend towards stall
const MF_DEFAULT = '0, 0.200, 0.205\n1, 0.310\n2, 0.420\n3, 0.530, 0.522\n4, 0.640\n5, 0.750\n6, 0.860, 0.828\n7, 0.970\n8, 1.080\n9, 1.190, 1.105\n10, 1.300\n11, 1.410\n12, 1.520, 1.318\n13, 1.630\n14, 1.740, 1.385';
const multifid = {
  id: 'multifidelity', title: 'Multifidelity correction, model discrepancy and reproducibility record', fidelity: 'numerical',
  summary: 'Corrects a cheap model with a handful of expensive results or measurements: a scale factor plus a Gaussian-process discrepancy that also reports its own uncertainty. Compares simple additive and multiplicative bridges, and stamps the study with a reproducibility record.',
  equations: ['Multifidelity model calibration equations', 'Gaussian process discrepancy modelling', 'Bayesian model updating', 'Least-squares calibration equations'],
  inputs: [
    { key: 'data', label: 'Data: x, low-fidelity[, high-fidelity]', type: 'text', default: MF_DEFAULT, group: 'Data', help: 'One x per line with the low-fidelity value; add the high-fidelity value or measurement where available. The default is a SYNTHETIC lift-curve illustration.' },
    { key: 'noise', label: 'Standard uncertainty of the high-fidelity values', unit: 'unit of y', default: 0.005, min: 0, group: 'Data' },
    { key: 'fit_rho', label: 'Fit a scale factor ρ on the low-fidelity model', type: 'bool', default: true, group: 'Model', help: 'Off = pure additive discrepancy (ρ = 1)' },
    { key: 'seed', label: 'Random seed recorded with the study', unit: '', default: 2025, min: 1, max: 1e9, step: 1, discrete: true, group: 'Reproducibility' },
    { key: 'note', label: 'Study note for the record', type: 'text', default: '', group: 'Reproducibility' },
  ],
  defaults: () => ({}),
  run(i, ctx) {
    const warnings = [], t = parseTable(i.data, 2); let rows = t.rows, synthetic = String(i.data).trim() === MF_DEFAULT;
    if (t.bad) warnings.push(`${t.bad} data line(s) could not be read and were skipped.`);
    if (rows.filter((r) => Number.isFinite(r[2])).length < 3 || rows.length < 4) { rows = parseTable(MF_DEFAULT, 2).rows; synthetic = true; warnings.push('At least four points, three of them with high-fidelity values, are needed; the built-in synthetic dataset was used.'); }
    if (synthetic) warnings.push('The dataset is synthetic; it only illustrates the method.');
    rows = rows.slice().sort((a, b) => a[0] - b[0]);
    const x = rows.map((r) => r[0]), lf = rows.map((r) => r[1]), H = rows.map((r, k) => k).filter((k) => Number.isFinite(rows[k][2])), xh = H.map((k) => x[k]), lh = H.map((k) => lf[k]), yh = H.map((k) => rows[k][2]), nh = H.length;
    const rhoOf = (l, y) => (i.fit_rho ? N.dot(l, y) / (N.dot(l, l) || 1) : 1);
    const bridges = {
      'Low fidelity only': () => (q, l) => ({ m: l, s: 0 }), 'Additive constant': (xs, l, y) => { const d = N.mean(y.map((v, k) => v - l[k])); return (q, lq) => ({ m: lq + d, s: 0 }); }, 'Multiplicative factor': (xs, l, y) => { const r = N.dot(l, y) / (N.dot(l, l) || 1); return (q, lq) => ({ m: r * lq, s: 0 }); },
      'Scale + Gaussian-process discrepancy': (xs, l, y, ell) => { const r = rhoOf(l, y), g = gp1(xs, y.map((v, k) => v - r * l[k]), i.noise, ell); return Object.assign((q, lq) => { const p = g.predict(q); return { m: r * lq + p.m, s: p.s }; }, { rho: r, ell: g.ell }); },
    };
    const full = bridges['Scale + Gaussian-process discrepancy'](xh, lh, yh), names = Object.keys(bridges), loo = names.map((nm) => rms(N.range(nh, (k) => { const keep = N.range(nh).filter((j) => j !== k), f = bridges[nm](keep.map((j) => xh[j]), keep.map((j) => lh[j]), keep.map((j) => yh[j]), full.ell); return f(xh[k], lh[k]).m - yh[k]; })));
    const pred = x.map((q, k) => full(q, lf[k])), best = N.argmin(loo), gain = 100 * (1 - loo[3] / (loo[0] || 1e-300)), smax = Math.max(...pred.map((p) => p.s)), kmax = N.argmax(pred.map((p) => p.s));
    if (nh < 5) warnings.push('Fewer than five high-fidelity points: the leave-one-out comparison is rough.');
    if (x[0] < xh[0] || x[x.length - 1] > xh[nh - 1]) warnings.push('Some prediction points lie outside the range of the high-fidelity data: the correction is an extrapolation there and its uncertainty band widens to the prior.');
    const record = { inputs_hash: fnv1a(stableJson(i)), case_hash: ctx?.case ? fnv1a(JSON.stringify(ctx.case)) : 'no case', seed: i.seed, software: VERSION, data_points: rows.length, high_fidelity_points: nh };
    return {
      kpis: [
        { key: 'rho_scale', label: 'Scale factor ρ on the low-fidelity model', value: full.rho, unit: '-' }, { key: 'gp_length', label: 'Discrepancy correlation length', value: full.ell, unit: 'unit of x' },
        { key: 'loo_rmse_lf', label: 'Leave-one-out RMSE: low fidelity only', value: loo[0], unit: '' }, { key: 'loo_rmse_add', label: 'Leave-one-out RMSE: additive constant', value: loo[1], unit: '' }, { key: 'loo_rmse_mult', label: 'Leave-one-out RMSE: multiplicative factor', value: loo[2], unit: '' },
        { key: 'loo_rmse_gp', label: 'Leave-one-out RMSE: scale + discrepancy', value: loo[3], unit: '', status: loo[3] <= loo[0] ? 'ok' : 'warn' }, { key: 'error_reduction_pct', label: 'Error reduction against the uncorrected model', value: gain, unit: '%', status: gain > 0 ? 'ok' : 'warn' },
        { key: 'max_pred_sd', label: 'Largest predictive standard deviation', value: smax, unit: '', note: `At x = ${x[kmax]}` }, { key: 'best_bridge', label: 'Best bridge (1 none, 2 additive, 3 multiplicative, 4 discrepancy)', value: best + 1, unit: '' },
        { key: 'inputs_hash_int', label: 'Reproducibility hash of the inputs (integer form)', value: parseInt(record.inputs_hash, 16), unit: '', note: `Hex ${record.inputs_hash}` },
      ],
      plots: [
        { type: 'line', title: 'Multifidelity prediction', xlabel: 'x', ylabel: 'y', series: [{ name: 'Low fidelity', x, y: lf, style: 'dash' }, { name: 'High-fidelity data', x: xh, y: yh, style: 'points' }, { name: 'Corrected prediction', x, y: pred.map((p) => p.m) }, { name: '+2σ', x, y: pred.map((p) => p.m + 2 * p.s), style: 'dash' }, { name: '−2σ', x, y: pred.map((p) => p.m - 2 * p.s), style: 'dash' }] },
        { type: 'line', title: 'Model discrepancy δ(x) = high fidelity − ρ·low fidelity', xlabel: 'x', ylabel: 'Discrepancy', series: [{ name: 'Gaussian-process mean', x, y: pred.map((p, k) => p.m - full.rho * lf[k]) }, { name: 'Observed', x: xh, y: yh.map((v, k) => v - full.rho * lh[k]), style: 'points' }], annotations: [{ y: 0, label: 'No discrepancy' }] },
        { type: 'bar', title: 'Leave-one-out error of each bridge', ylabel: 'RMSE', categories: names, series: [{ name: 'Leave-one-out RMSE', y: loo }] },
      ],
      tables: [{ title: 'Predictions', columns: ['x', 'Low fidelity', 'High fidelity', 'Corrected', 'Standard deviation'], rows: rows.map((r, k) => [x[k], lf[k], Number.isFinite(r[2]) ? r[2] : '—', pred[k].m, pred[k].s]) },
        { title: 'Reproducibility record', columns: ['Item', 'Value'], rows: [['Inputs hash (FNV-1a, 32 bit)', record.inputs_hash], ['Case hash', record.case_hash], ['Seed', record.seed], ['Software', record.software], ['Data points', record.data_points], ['High-fidelity points', record.high_fidelity_points], ['Note', String(i.note || '')]] }],
      outputs: { synthetic_data: synthetic ? 1 : 0, inputs_hash: record.inputs_hash },
      warnings, models: ['Scale factor by least squares (Kennedy–O’Hagan autoregressive form, simplified)', 'Gaussian-process discrepancy with squared-exponential kernel and likelihood length-scale search', 'Additive and multiplicative bridge functions', 'Leave-one-out cross-validation', 'FNV-1a input hash for the reproducibility record'],
      assumptions: ['The scale factor is fitted first and the discrepancy on the remainder (not a joint Bayesian inference)', 'The discrepancy is smooth and stationary', 'The length scale is frozen during cross-validation', 'High-fidelity values carry independent Gaussian uncertainty'],
    };
  },
  verify() {
    const xs = N.linspace(0, 10, 11), lf = xs.map((x) => 0.1 * x + 0.2), txt = xs.map((x, k) => `${x}, ${lf[k]}${k % 2 === 0 ? `, ${2 * lf[k]}` : ''}`).join('\n'), r = multifid.run({ data: txt, noise: 1e-6, fit_rho: true, seed: 1, note: '' }), o = N.kv(r);
    return [
      N.check('FNV-1a hash of the empty string is the offset basis', parseInt(fnv1a(''), 16), 0x811c9dc5, 0, 'FNV-1a definition'), N.check('FNV-1a hash of "a"', parseInt(fnv1a('a'), 16), 0xe40c292c, 0, 'FNV-1a reference vector'),
      N.check('Scale factor recovered when high fidelity = 2 × low fidelity', o.rho_scale, 2, 1e-9, 'Manufactured data'),
      N.check('Corrected prediction at an unseen point', r.tables[0].rows[3][3], 2 * lf[3], 1e-6, 'Manufactured data'),
      N.check('Gaussian process interpolates its data', gp1([0, 1, 2, 3], [0, 0.8, 0.9, 0.1], 0).predict(2).m, 0.9, 1e-5, 'Interpolation property'),
    ];
  },
  recommend(res) {
    const o = res.outputs, out = [], nm = res.plots[2].categories[o.best_bridge - 1];
    out.push({ severity: o.error_reduction_pct > 20 ? 'advise' : 'info', title: `Best correction: ${nm}`, detail: `Scale + discrepancy cuts the leave-one-out error by ${o.error_reduction_pct.toFixed(0)}% against the uncorrected low-fidelity model.`, action: o.best_bridge === 4 ? 'Use the corrected model for trade studies and reserve the expensive solver or test for the points where the predictive band is widest.' : 'A simple bridge is enough here; the discrepancy has no resolvable structure. Keep the simpler model.', basis: 'Leave-one-out cross-validation' });
    out.push({ severity: 'advise', title: 'Place the next expensive run where the band is widest', detail: res.kpis.find((k) => k.key === 'max_pred_sd').note + ` (σ = ${Number(o.max_pred_sd).toPrecision(2)}).`, action: 'Adding one high-fidelity point there reduces uncertainty most per unit of computing or test cost and energy.', basis: 'Gaussian-process predictive variance' });
    out.push({ severity: 'info', title: 'Archive the reproducibility record', detail: `Inputs hash ${o.inputs_hash}.`, action: 'Store the hash, seed and software line with the results; the same inputs and seed reproduce every number in this suite exactly.', basis: 'Reproducibility of numerical evidence' });
    return out;
  },
};

export default {
  id: 'vvuq', n: 25,
  tagline: 'How far the numbers can be trusted: is the code right, is the grid fine enough, does the model match test data, and how uncertain is the answer.',
  analyses: [mms, gciCalc, validation, propagation, sobol, bayesCal, multifid],
  consumes: [
    { from: 'performance', keys: ['LD_max'], why: 'Mean lift-to-drag ratio for the Breguet uncertainty model' },
    { from: 'propulsion', keys: ['tsfc_kg_Ns'], why: 'Mean fuel consumption for the Breguet uncertainty model' },
  ],
  provides: [{ key: 'gci_fine_pct', label: 'Fine-grid GCI', unit: '%' }, { key: 'observed_order', label: 'Observed order of accuracy', unit: '-' }, { key: 'validation_rmse', label: 'Validation RMSE', unit: '' }],
  handoff: [
    { model: 'Verification of external CFD / FE codes on their own grids', why: 'Manufactured-solution tests here exercise built-in model equations; production solvers must be verified in their own environment', tool: 'Solver-specific MMS and regression suites; then paste the grid results into the GCI calculator' },
    { model: 'Full Kennedy–O’Hagan Bayesian calibration with joint inference of parameters, discrepancy and hyper-parameters', why: 'Needs high-dimensional MCMC over model runs; a sequential scale-then-discrepancy fit is provided', tool: 'Dedicated Bayesian calibration frameworks' },
    { model: 'Sparse-grid collocation and adaptive polynomial chaos in many dimensions', why: 'Tensor grids and total-degree regression are practical only for a few inputs', tool: 'UQ toolkits with sparse grids and adaptive bases' },
    { model: 'Experimental–computational data assimilation (Kalman-type filters on field data)', why: 'Requires time-resolved measurement streams and solver state access', tool: 'Data-assimilation frameworks coupled to the solver' },
    { model: 'Formal credibility assessment and validation hierarchy management', why: 'An organisational process with independent experimental evidence at material, component, subsystem and aircraft level', tool: 'ASME V&V 10/20/40-style programme with a managed evidence database' },
  ],
};
