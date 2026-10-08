// Suite 23 — Multidisciplinary Design Analysis and Optimisation (MDAO).
// A self-contained low-order aircraft model (drag build-up with Korn wave drag, wing bending-material sizing with
// inertia relief, rubber engine, Breguet mission, field and approach constraints) calibrated to the shared case,
// wrapped in MDA solvers, a constrained optimiser, NSGA-II, surrogates with expected improvement, coupled
// sensitivities and robust / reliability-based design. Rotorcraft get a rotor-sizing optimisation.
// The model ranks design trends; mass and cost coefficients are illustrative and must be calibrated before use.

import * as N from '../core/numerics.js';
import { isa, G0, RHO0 } from '../core/atmosphere.js';
import { METALS, BATTERIES } from '../data/materials.js';

// ---- disciplines ----------------------------------------------------------------------------
/** AR and sweep trend common to Raymer's Oswald-efficiency correlations; used as a ratio to the baseline value. */
const eShape = (x) => Math.max(0.05, 1 - 0.045 * x.AR ** 0.68) * Math.cos(N.rad(x.sweep)) ** 0.15;
/** Aerodynamics at weight W [N]: wetted-area drag build-up, induced drag and Korn-equation wave drag. */
function aeroAt(P, x, at, W) {
  const M = P.mach, V = M * at.a, q = 0.5 * at.rho * V * V, cosL = Math.cos(N.rad(x.sweep)), b = Math.sqrt(x.AR * x.S), mac = x.S / b;
  const Re = Math.max(1e4, (at.rho * V * mac) / at.mu), cf = 0.455 / (Math.log10(Re) ** 2.58 * (1 + 0.144 * M * M) ** 0.65);
  const FF = (1 + 2 * x.tc + 100 * x.tc ** 4) * 1.34 * Math.max(M, 0.2) ** 0.18 * cosL ** 0.28, fWing = cf * FF * 2.04 * x.S;
  const CD0 = ((P.fRest || 0) + fWing) / x.S, e = N.clamp((P.e_ref * eShape(x)) / (P.eShape0 || eShape(x)), 0.3, 0.98), k = 1 / (Math.PI * x.AR * e), CL = W / (q * x.S);
  const Mdd = P.kappa / cosL - x.tc / cosL ** 2 - CL / (10 * cosL ** 3), Mcr = Mdd - Math.cbrt(0.1 / 80), CDw = M > Mcr ? 20 * (M - Mcr) ** 4 : 0, CD = CD0 + k * CL * CL + CDw;
  return { V, q, b, CL, CD0, e, k, CDw, CD, LD: CL / CD, Mdd, fWing };
}
/**
 * Wing mass [kg] from bending-material sizing: elliptic air load at ultimate load factor, relieved by the distributed
 * wing and wing-fuel inertia, integrated to shear and bending moment; caps sized to the allowable stress on the
 * swept structural span, webs to shear; non-optimum factor and an area-proportional secondary mass on top.
 */
function wingMass(P, x, mtow, mRelief) {
  const b = Math.sqrt(x.AR * x.S), s = b / 2, cosL = Math.cos(N.rad(x.sweep)), cr = (2 * x.S) / (b * (1 + P.taper)), L = P.n_ult * mtow * G0, n = P.nSta, tau = 0.55 * P.sigma_allow;
  let V = 0, M = 0, mass = 0, pP = 0, fP = 0, etaP = 1;
  for (let j = n; j >= 0; j--) { // cosine-spaced stations, tip to root
    const eta = Math.sin((Math.PI * j) / (2 * n)), c = cr * (1 - (1 - P.taper) * eta), dy = (etaP - eta) * s;
    const p = ((4 * L) / (Math.PI * b)) * Math.sqrt(Math.max(0, 1 - eta * eta)) - (P.n_ult * G0 * mRelief * c) / x.S, Vn = V + 0.5 * (p + pP) * dy;
    M += 0.5 * (V + Vn) * dy; V = Vn;
    const h = Math.max(1e-4, P.box_h_frac * x.tc * c), f = (P.rho_mat * ((2 * Math.abs(M)) / (cosL * P.sigma_allow * h) + Math.abs(V) / tau)) / cosL;
    mass += 0.5 * (f + fP) * dy; pP = p; fP = f; etaP = eta;
  }
  return { m: 2 * mass * P.kNo + P.k_sec * x.S, box: 2 * mass, Mroot: M };
}
const engMass = (P, mtow) => (P.jet ? P.eng_kg_N * P.tw * mtow * G0 : P.eng_kg_W * P.pw * mtow);
/** Fuel volume of the wing box [m³] (Torenbeek's estimate) plus any centre / fuselage tank. */
const tankVol = (P, x) => (0.54 * x.S * x.S * x.tc * (1 + P.taper + P.taper ** 2)) / (Math.sqrt(x.AR * x.S) * (1 + P.taper) ** 2) + P.tank_extra_m3;
function thrustCruise(P, at, V, mtow) {
  if (P.jet) { const M = V / at.a, d0 = at.delta * (1 + 0.2 * M * M) ** 3.5, r = Math.sqrt(M); return P.tw * mtow * G0 * Math.max(0.05, P.bpr >= 1.5 ? d0 * (1 - 0.49 * r) : d0 * (1 - 0.3 * r)); }
  const lap = P.ptype === 'piston' ? Math.max(0, (at.sigma - 0.1325) / 0.8675) : P.electric ? 1 : at.sigma ** 0.7;
  return (P.eta_prop * P.pw * mtow * lap) / Math.max(V, 1);
}
function tofl(P, x, mtow, kTO) {
  const ws = (mtow * G0) / x.S, Vlo = 1.1 * Math.sqrt((2 * ws) / (RHO0 * P.CLmax_to)), W = mtow * G0;
  const T = P.jet ? P.tw * W * (1 - (P.bpr >= 1.5 ? 0.49 : 0.3) * Math.sqrt((0.7 * Vlo) / 340.3)) : Math.min(P.tw > 0 ? P.tw * W : Infinity, (P.eta_prop * P.pw * mtow) / (0.7 * Vlo));
  return (kTO * 1.21 * ws) / (G0 * RHO0 * P.CLmax_to * Math.max(1e-3, T / W - 0.03));
}
const hoverEnergy = (P, mtow) => (P.disk_area > 0 && P.hover_min > 0 ? ((mtow * G0) ** 1.5 / Math.sqrt(2 * RHO0 * P.disk_area) / 0.7 / P.eta_elec) * P.hover_min * 60 : 0);

/** Calibrate the model to the baseline design and reference masses. Returns the parameter set P. */
function makeModel(i) {
  const P = { ...i, range_m: i.range_km * 1e3, electric: i.ptype === 'electric', jet: i.ptype === 'turbofan' || i.ptype === 'turbojet', nSta: Math.max(4, Math.round(i.nSta || 12)), espec: i.batt_Wh_kg * 3600, warn: [], kNo: i.k_no, dISA: 0 };
  const x0 = (P.x0 = { AR: i.AR0, S: i.S0, sweep: i.sweep0, tc: i.tc0, h: i.h0 }), at = isa(x0.h, 0);
  P.eShape0 = eShape(x0); P.fRest = 0;
  const fW = aeroAt(P, x0, at, i.mtow_ref * G0).fWing; P.fRest = i.CD0_ref * i.S0 - fW;
  if (P.fRest < 0.25 * i.CD0_ref * i.S0) { P.fRest = 0.25 * i.CD0_ref * i.S0; P.warn.push('The reference CD0 is lower than the wing friction estimate allows; the non-wing drag area was floored at 25% of the total and the model CD0 is higher than the reference.'); }
  const fuelRef = Math.max(0, i.mtow_ref - i.oew_ref - i.payload), relief = (mw) => mw + (P.electric ? 0 : 0.5 * Math.min(fuelRef, P.fuel_rho * tankVol(P, x0)));
  let w0 = 0.1 * i.mtow_ref;
  if (i.wing_mass_ref > 0) { const box = wingMass({ ...P, kNo: 1, k_sec: 0 }, x0, i.mtow_ref, relief(i.wing_mass_ref)).m, k = (i.wing_mass_ref - i.k_sec * i.S0) / box; P.kNo = N.clamp(k, 0.8, 5); if (k !== P.kNo) P.warn.push(`The reference wing mass implies a non-optimum factor of ${k.toFixed(2)}, outside 0.8–5; it was limited, so the model wing mass differs from the reference.`); }
  for (let k = 0; k < 12; k++) w0 = wingMass(P, x0, i.mtow_ref, relief(w0)).m;
  P.wing0 = w0; P.eng0 = engMass(P, i.mtow_ref); P.batt0 = P.electric ? (i.batt_ref_kWh * 3.6e6) / P.espec : 0; P.fuelRef = P.electric ? P.batt0 : fuelRef;
  P.mFixed = i.oew_ref - w0 - P.eng0 - P.batt0;
  if (P.mFixed < 0.15 * i.oew_ref) { P.mFixed = 0.15 * i.oew_ref; P.warn.push('Wing, engine and battery estimates exceed the reference empty mass; the fixed mass was floored at 15% of it. Check the mass coefficients.'); }
  const t1 = tofl(P, x0, i.mtow_ref, 1); P.kTO = i.tofl_ref > 0 ? N.clamp(i.tofl_ref / t1, 0.8, 4) : 1.5;
  return P;
}
/** Coupled disciplines for one design x. Coupling vector y = [take-off mass, wing mass, fuel or battery mass]. */
function couple(P, x) {
  const at = isa(x.h, P.dISA), cap = P.electric ? Infinity : P.fuel_rho * tankVol(P, x);
  const wing = (y) => wingMass(P, x, y[0], y[1] + (P.electric ? 0 : 0.5 * Math.min(y[2], cap)));
  const energy = (mtow, mE) => {
    if (P.electric) { const A = aeroAt(P, x, at, mtow * G0); return (((mtow * G0 * P.range_m) / (A.LD * P.eta_prop * P.eta_elec * P.allow_frac) + hoverEnergy(P, mtow)) * (1 + P.reserve_frac)) / P.espec; }
    const A = aeroAt(P, x, at, (mtow - (0.5 * mE) / (1 + P.reserve_frac)) * G0), c = P.jet ? P.tsfc : (P.bsfc * A.V) / P.eta_prop;
    return mtow * (1 - P.allow_frac * Math.exp((-P.range_m * G0 * c) / (A.V * A.LD))) * (1 + P.reserve_frac);
  };
  const sum = (mtow, mw, mE) => P.mFixed + mw + engMass(P, mtow) + P.payload + mE;
  return {
    y0: [P.mtow_ref, P.wing0, P.fuelRef || 0.1 * P.mtow_ref],
    jacobi: (y) => [sum(y[0], y[1], y[2]), wing(y).m, energy(y[0], y[2])],
    seidel: (y) => { const mw = wing(y).m, mE = energy(y[0], y[2]); return [sum(y[0], mw, mE), mw, mE]; },
    report(y) {
      const [mtow, mw, mE] = y, mis = P.electric ? 0 : mE / (1 + P.reserve_frac), A = aeroAt(P, x, at, mtow * G0), D = (mtow * G0) / A.LD, Tav = thrustCruise(P, at, A.V, mtow), tf = tofl(P, x, mtow, P.kTO);
      const vapp = 1.3 * Math.sqrt((2 * (mtow - 0.8 * mis) * G0) / (RHO0 * x.S * P.CLmax_land)), Amid = aeroAt(P, x, at, (mtow - 0.5 * mis) * G0);
      const g = [1 - tf / P.runway_m, 1 - vapp / P.vapp_max, 1 - A.b / P.span_max, 1 - A.CL / P.CL_buffet, Tav / D - 1, P.electric ? 1 : 1 - mE / cap];
      return { mtow, wing: mw, energy_kg: mE, mission_kg: mis, eng: engMass(P, mtow), oew: mtow - P.payload - (P.electric ? 0 : mE), kWh: P.electric ? (mE * P.espec) / 3.6e6 : 0, LD: Amid.LD, CL: A.CL, CD0: A.CD0, CDw: Amid.CDw, e: A.e, b: A.b, tofl: tf, vapp, cap, thrust_margin: Tav / D - 1, g, Mroot: wing(y).Mroot };
    },
  };
}
const CON_NAMES = ['Take-off field length', 'Approach speed', 'Wing span limit', 'Cruise lift coefficient (buffet)', 'Cruise thrust', 'Fuel volume'];
/** Fixed-point (Gauss–Seidel or Jacobi) or Newton solution of y = G(y). Returns {y, iters, hist, converged}. */
function solveMDA(G, y0, { method = 'seidel', tol = 1e-9, maxIter = 200, Gj } = {}) {
  let y = y0.slice(); const hist = [], scale = Math.abs(y0[0]) || 1, n = y.length;
  if (method !== 'newton') {
    for (let it = 1; it <= maxIter; it++) { const yn = G(y), r = Math.sqrt(N.sum(yn.map((v, j) => (v - y[j]) ** 2))) / scale; hist.push(r); y = yn; if (!(r < 1e6) || !Number.isFinite(r)) return { y, iters: it, hist, converged: false }; if (r < tol) return { y, iters: it, hist, converged: true }; }
    return { y, iters: maxIter, hist, converged: false };
  }
  const R = (v) => { const g = (Gj || G)(v); return v.map((q, j) => q - g[j]); };
  for (let it = 1; it <= Math.min(maxIter, 40); it++) {
    const r = R(y), nr = N.norm(r) / scale; hist.push(nr);
    if (!Number.isFinite(nr) || nr > 1e6) return { y, iters: it, hist, converged: false };
    if (nr < tol) return { y, iters: it - 1, hist, converged: true };
    const J = N.zeros(n); for (let j = 0; j < n; j++) { const h = 1e-6 * Math.max(1, Math.abs(y[j])), yp = y.slice(); yp[j] += h; const rp = R(yp); for (let k = 0; k < n; k++) J[k][j] = (rp[k] - r[k]) / h; }
    let dy; try { dy = N.solve(J, r.map((v) => -v)); } catch { return { y, iters: it, hist, converged: false }; }
    y = y.map((v, j) => v + dy[j]);
  }
  return { y, iters: 40, hist, converged: false };
}
/** Converged design evaluation. */
function evalDesign(P, x, tol = 1e-8) {
  const c = couple(P, x), s = solveMDA(c.seidel, c.y0, { tol, maxIter: 120 });
  if (!s.converged || !(s.y[0] > 0) || s.y[0] > 20 * P.mtow_ref) return { ok: false, iters: s.iters };
  return { ok: true, iters: s.iters, y: s.y, ...c.report(s.y) };
}
const OBJ = {
  fuel: { label: 'Mission energy mass (fuel incl. reserves, or battery)', unit: 'kg', f: (o) => o.energy_kg },
  mtow: { label: 'Take-off mass', unit: 'kg', f: (o) => o.mtow },
  doc: { label: 'Cost proxy per flight', unit: 'USD', f: (o, P) => (P.electric ? P.elec_usd_kWh * o.kWh : P.fuel_usd_kg * o.mission_kg) + P.c_oew * o.oew },
};

// ---- optimisers -----------------------------------------------------------------------------
/**
 * Augmented-Lagrangian minimisation on the unit box. fn(u) -> { f, g: [g_j], ok } with g_j ≥ 0 feasible.
 * Inner unconstrained solves by Nelder–Mead, then a projected-gradient polish with central differences.
 */
function alOptimise(fn, u0, { outer = 5, inner = 80, mu0 = 10, polish = 8 } = {}) {
  const n = u0.length, clip = (u) => u.map((v) => N.clamp(v, 0, 1)); let evals = 0, mu = mu0, u = clip(u0);
  const ev = (v) => { evals++; const r = fn(clip(v)); return r.ok === false ? { f: 1e3, g: r.g || [], ok: false } : r; };
  let r = ev(u), lam = r.g.map(() => 0), viol = Math.max(0, ...r.g.map((g) => -g)); const hist = [{ f: r.f, viol, evals }];
  const L = (v) => { const q = ev(v); let s = q.f, pen = 0; if (q.ok === false) s += 1e3; for (let j = 0; j < q.g.length; j++) { const t = Math.max(0, lam[j] - mu * q.g[j]); s += (t * t - lam[j] * lam[j]) / (2 * mu); } for (const x of v) pen += Math.max(0, -x, x - 1) ** 2; return s + 1e3 * pen; };
  for (let k = 0; k < outer; k++) {
    const z = N.nelderMead((w) => L(w.map((v) => v - 1)), u.map((v) => v + 1), { maxIter: inner, tol: 1e-11, step: 0.06 }); // shifted by 1 so the simplex step is absolute
    u = clip(z.x.map((v) => v - 1)); r = ev(u);
    const vNew = Math.max(0, ...r.g.map((g) => -g)); lam = lam.map((l, j) => Math.max(0, l - mu * r.g[j]));
    if (vNew > 0.25 * viol && vNew > 1e-6) mu *= 4; viol = vNew; hist.push({ f: r.f, viol, evals });
  }
  for (let k = 0; k < polish; k++) { // projected gradient on the augmented Lagrangian
    const h = 1e-4, g = u.map((_, j) => { const a = u.slice(), b = u.slice(); a[j] = Math.min(1, u[j] + h); b[j] = Math.max(0, u[j] - h); return (L(a) - L(b)) / (a[j] - b[j]); }), L0 = L(u);
    let t = 0.05 / (N.norm(g) || 1), ok = false;
    for (let ls = 0; ls < 12; ls++) { const cand = clip(u.map((v, j) => v - t * g[j])); if (L(cand) < L0 - 1e-14) { u = cand; ok = true; break; } t *= 0.4; }
    if (!ok) break;
  }
  r = ev(u); hist.push({ f: r.f, viol: Math.max(0, ...r.g.map((g) => -g)), evals });
  return { u, f: r.f, g: r.g, lam, hist, evals, mu };
}
/** KKT check at u: least-squares multipliers on the active set and the stationarity residual on free variables. */
function kkt(fn, u, gTol = 3e-3) {
  const n = u.length, h = 1e-4, base = fn(u), gradF = new Array(n), gradG = base.g.map(() => new Array(n));
  for (let j = 0; j < n; j++) { const a = u.slice(), b = u.slice(); a[j] = Math.min(1, u[j] + h); b[j] = Math.max(0, u[j] - h); const ra = fn(a), rb = fn(b), d = a[j] - b[j]; gradF[j] = (ra.f - rb.f) / d; base.g.forEach((_, c) => (gradG[c][j] = (ra.g[c] - rb.g[c]) / d)); }
  const free = N.range(n).filter((j) => !((u[j] < 1e-4 && gradF[j] > 0) || (u[j] > 1 - 1e-4 && gradF[j] < 0)));
  let act = N.range(base.g.length).filter((c) => base.g[c] < gTol), lam = [];
  for (let pass = 0; pass < 6 && act.length && free.length; pass++) {
    try { lam = N.lstsq(free.map((j) => act.map((c) => gradG[c][j])), free.map((j) => gradF[j]), 1e-10); } catch { lam = act.map(() => 0); }
    const neg = lam.findIndex((l) => l < -1e-6); if (neg < 0) break; act = act.filter((_, k) => k !== neg); lam = [];
  }
  if (!act.length || !free.length) lam = act.map(() => 0);
  const res = free.map((j) => gradF[j] - N.sum(act.map((c, k) => lam[k] * gradG[c][j]))), nf = N.norm(free.map((j) => gradF[j]));
  return { gradF, gradG, lambda: base.g.map((_, c) => { const k = act.indexOf(c); return k >= 0 ? Math.max(0, lam[k]) : 0; }), residual: free.length ? N.norm(res) / Math.max(nf, 1e-9) : 0, active: act, free };
}
/** NSGA-II on the unit box for two objectives. fn(u) -> { f: [f1, f2], viol ≥ 0 }. Returns the final non-dominated set. */
function nsga2(fn, n, { pop = 40, gens = 30, seed = 1 } = {}) {
  const u = N.rng(seed), mk = (x) => { const r = fn(x); return { x, f: r.f, viol: r.viol, rank: 0, crowd: 0, data: r.data }; };
  const dom = (a, b) => (a.viol > 0 || b.viol > 0 ? a.viol < b.viol : a.f[0] <= b.f[0] && a.f[1] <= b.f[1] && (a.f[0] < b.f[0] || a.f[1] < b.f[1]));
  const rankAll = (Q) => {
    let rest = Q.slice(), rank = 0; const fronts = [];
    while (rest.length) { const F = rest.filter((a) => !rest.some((b) => b !== a && dom(b, a))); if (!F.length) { fronts.push(rest); rest.forEach((a) => (a.rank = rank)); break; } F.forEach((a) => { a.rank = rank; a.crowd = 0; }); for (const m of [0, 1]) { F.sort((a, b) => a.f[m] - b.f[m]); const span = F[F.length - 1].f[m] - F[0].f[m] || 1; F[0].crowd = F[F.length - 1].crowd = Infinity; for (let j = 1; j < F.length - 1; j++) F[j].crowd += (F[j + 1].f[m] - F[j - 1].f[m]) / span; } fronts.push(F); rest = rest.filter((a) => !F.includes(a)); rank++; }
    return fronts;
  };
  const better = (a, b) => a.rank < b.rank || (a.rank === b.rank && a.crowd > b.crowd), pick = (Q) => { const a = Q[Math.floor(u() * Q.length)], b = Q[Math.floor(u() * Q.length)]; return better(a, b) ? a : b; };
  let Q = N.lhs(pop, n, u).map(mk); rankAll(Q);
  for (let g = 0; g < gens; g++) {
    const kids = [];
    while (kids.length < pop) {
      const p1 = pick(Q).x, p2 = pick(Q).x, c1 = p1.slice(), c2 = p2.slice();
      for (let j = 0; j < n; j++) { // simulated binary crossover (η = 15) and polynomial mutation (η = 20)
        if (u() < 0.9) { const r = u(), bq = r < 0.5 ? (2 * r) ** (1 / 16) : (1 / (2 * (1 - r))) ** (1 / 16); c1[j] = 0.5 * ((1 + bq) * p1[j] + (1 - bq) * p2[j]); c2[j] = 0.5 * ((1 - bq) * p1[j] + (1 + bq) * p2[j]); }
        for (const c of [c1, c2]) { if (u() < 1 / n) { const r = u(); c[j] += r < 0.5 ? (2 * r) ** (1 / 21) - 1 : 1 - (2 * (1 - r)) ** (1 / 21); } c[j] = N.clamp(c[j], 0, 1); }
      }
      kids.push(mk(c1), mk(c2));
    }
    const fronts = rankAll(Q.concat(kids)); Q = [];
    for (const F of fronts) { if (Q.length + F.length <= pop) Q.push(...F); else { F.sort((a, b) => b.crowd - a.crowd); Q.push(...F.slice(0, pop - Q.length)); break; } }
  }
  return rankAll(Q)[0].filter((a) => a.viol <= 0).sort((a, b) => a.f[0] - b.f[0]);
}
// Gaussian process (ordinary Kriging, anisotropic squared-exponential kernel) on inputs scaled to the unit box
const cholSolve = (L, b) => { const n = b.length, y = new Array(n), x = new Array(n); for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i][k] * y[k]; y[i] = s / L[i][i]; } for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k]; x[i] = s / L[i][i]; } return x; };
function gpFit(X, yRaw, ellFixed) {
  const n = X.length, d = X[0].length, ym = N.mean(yRaw), ys = N.std(yRaw) || 1, y = yRaw.map((v) => (v - ym) / ys), one = new Array(n).fill(1);
  const kern = (a, b, ell) => { let s = 0; for (let j = 0; j < d; j++) s += ((a[j] - b[j]) / ell[j]) ** 2; return Math.exp(-0.5 * s); };
  const build = (ell) => { const K = N.range(n, (i) => N.range(n, (j) => kern(X[i], X[j], ell) + (i === j ? 1e-8 : 0))); let L; try { L = N.cholesky(K); } catch { return null; } const Ky = cholSolve(L, y), K1 = cholSolve(L, one), mu = N.dot(one, Ky) / N.dot(one, K1), al = Ky.map((v, j) => v - mu * K1[j]), s2 = Math.max(1e-300, N.dot(y.map((v) => v - mu), al) / n); let ld = 0; for (let j = 0; j < n; j++) ld += Math.log(L[j][j]); return { ell, L, mu, al, s2, K1, lnL: -0.5 * n * Math.log(s2) - ld }; };
  let best = ellFixed ? build(ellFixed) : null;
  if (!best) { const grid = N.logspace(0.12, 3, 7), rec = (j, ell) => { if (j === d) { const m = build(ell.slice()); if (m && (!best || m.lnL > best.lnL)) best = m; return; } for (const g of grid) { ell[j] = g; rec(j + 1, ell); } }; rec(0, new Array(d)); }
  const m = best, denom = N.dot(one, m.K1);
  const predict = (x) => { const ks = X.map((p) => kern(x, p, m.ell)), Kk = cholSolve(m.L, ks), mean = m.mu + N.dot(ks, m.al), v = m.s2 * Math.max(0, 1 - N.dot(ks, Kk) + (1 - N.dot(one, Kk)) ** 2 / denom); return { m: ym + ys * mean, s: ys * Math.sqrt(v) }; };
  const diag = N.range(n, (j) => cholSolve(m.L, N.range(n, (k) => (k === j ? 1 : 0)))[j]), loo = m.al.map((a, j) => (ys * a) / (diag[j] - (m.K1[j] * m.K1[j]) / denom)); // leave-one-out residuals in closed form (ordinary Kriging)
  return { predict, loo, ell: m.ell, lnL: m.lnL };
}
const expImp = (fmin, m, s) => { if (s < 1e-12) return Math.max(0, fmin - m); const z = (fmin - m) / s; return (fmin - m) * N.normCdf(z) + (s * Math.exp(-0.5 * z * z)) / Math.sqrt(2 * Math.PI); };
const quadBasis = (x) => { const r = [1, ...x]; for (let a = 0; a < x.length; a++) for (let b = a; b < x.length; b++) r.push(x[a] * x[b]); return r; };
/** First-order reliability method (Hasofer–Lind / Rackwitz–Fiessler) in standard-normal space. g(u) ≤ 0 is failure. */
function form(g, n, maxIter = 40) {
  let u = new Array(n).fill(0), beta = 0, grad = [], it = 0;
  for (; it < maxIter; it++) {
    const g0 = g(u); grad = u.map((_, j) => { const a = u.slice(), b = u.slice(); a[j] += 1e-3; b[j] -= 1e-3; return (g(a) - g(b)) / 2e-3; });
    const n2 = N.dot(grad, grad); if (!(n2 > 1e-30)) break;
    const un = grad.map((d) => ((N.dot(grad, u) - g0) / n2) * d), dv = N.norm(N.vadd(un, u, -1)); u = un; beta = N.norm(u) * (g(new Array(n).fill(0)) >= 0 ? 1 : -1);
    if (dv < 1e-6) break;
  }
  const ng = N.norm(grad) || 1;
  return { beta, pf: N.normCdf(-beta), u, alpha: grad.map((d) => -d / ng), iters: it + 1 };
}

// ---- shared inputs --------------------------------------------------------------------------
const MODEL = [
  { key: 'AR0', label: 'Baseline aspect ratio', unit: '-', default: 9.5, min: 3, max: 30, group: 'Baseline design' },
  { key: 'S0', label: 'Baseline wing area', unit: 'm²', default: 122.6, min: 0.05, group: 'Baseline design' },
  { key: 'sweep0', label: 'Baseline quarter-chord sweep', unit: 'deg', default: 25, min: 0, max: 50, group: 'Baseline design' },
  { key: 'tc0', label: 'Baseline thickness ratio', unit: '-', default: 0.115, min: 0.05, max: 0.25, group: 'Baseline design' },
  { key: 'taper', label: 'Taper ratio', unit: '-', default: 0.24, min: 0.1, max: 1, group: 'Baseline design' },
  { key: 'h0', label: 'Baseline cruise altitude', unit: 'm', default: 10668, min: 0, max: 15000, group: 'Baseline design' },
  { key: 'mach', label: 'Cruise Mach number', unit: '-', default: 0.78, min: 0.03, max: 0.9, group: 'Mission' },
  { key: 'range_km', label: 'Design range', unit: 'km', default: 4500, min: 1, group: 'Mission' },
  { key: 'payload', label: 'Design payload', unit: 'kg', default: 16600, min: 0, group: 'Mission' },
  { key: 'allow_frac', label: 'Non-cruise mass fraction', unit: '-', default: 0.97, min: 0.8, max: 1, group: 'Mission', help: 'Mass ratio for taxi, take-off, climb and descent allowances (energy ratio for electric)' },
  { key: 'reserve_frac', label: 'Reserves as a fraction of mission fuel', unit: '-', default: 0.2, min: 0, max: 1, group: 'Mission' },
  { key: 'mtow_ref', label: 'Reference take-off mass', unit: 'kg', default: 78000, min: 1, group: 'Reference masses', help: 'Used only to calibrate the fixed mass and start the iteration' },
  { key: 'oew_ref', label: 'Reference operating empty mass', unit: 'kg', default: 42600, min: 0.5, group: 'Reference masses' },
  { key: 'wing_mass_ref', label: 'Reference wing mass (0 = model estimate)', unit: 'kg', default: 0, min: 0, group: 'Reference masses', help: 'From Suite 2 when available; calibrates the non-optimum factor' },
  { key: 'CD0_ref', label: 'Baseline zero-lift drag coefficient', unit: '-', default: 0.0205, min: 0.004, max: 0.2, group: 'Aerodynamics', help: 'From Suite 1 when available' },
  { key: 'e_ref', label: 'Baseline Oswald efficiency', unit: '-', default: 0.8, min: 0.4, max: 0.98, group: 'Aerodynamics' },
  { key: 'kappa', label: 'Korn aerofoil technology factor', unit: '-', default: 0.95, min: 0.8, max: 0.97, group: 'Aerodynamics', help: '0.87 conventional sections, 0.95 supercritical' },
  { key: 'CL_buffet', label: 'Maximum cruise lift coefficient', unit: '-', default: 0.7, min: 0.2, max: 1.6, group: 'Aerodynamics', help: 'Buffet or stall-margin limit at start of cruise' },
  { key: 'CLmax_to', label: 'CLmax take-off', unit: '-', default: 2.2, min: 0.5, max: 4, group: 'Aerodynamics' },
  { key: 'CLmax_land', label: 'CLmax landing', unit: '-', default: 2.9, min: 0.5, max: 4.5, group: 'Aerodynamics' },
  { key: 'ptype', label: 'Powerplant type', type: 'select', options: ['turbofan', 'turbojet', 'turboprop', 'piston', 'electric'], default: 'turbofan', group: 'Propulsion' },
  { key: 'tw', label: 'Static thrust-to-weight ratio (held)', unit: '-', default: 0.314, min: 0, max: 2, group: 'Propulsion', help: 'The engine is rubberised so this ratio stays constant as the mass changes' },
  { key: 'pw', label: 'Power-to-mass ratio (held)', unit: 'W/kg', default: 0, min: 0, max: 3000, group: 'Propulsion' },
  { key: 'bpr', label: 'Bypass ratio', unit: '-', default: 5.7, min: 0, max: 15, group: 'Propulsion' },
  { key: 'tsfc', label: 'Cruise TSFC', unit: 'kg/N/s', default: 1.62e-5, min: 0, group: 'Propulsion', help: 'From Suite 7 when available' },
  { key: 'bsfc', label: 'BSFC', unit: 'kg/W/s', default: 8e-8, min: 0, group: 'Propulsion' },
  { key: 'eta_prop', label: 'Propeller efficiency', unit: '-', default: 0.82, min: 0.3, max: 0.95, group: 'Propulsion' },
  { key: 'eng_kg_N', label: 'Installed engine mass per unit thrust', unit: 'kg/N', default: 0.026, min: 0, max: 0.2, group: 'Mass coefficients (illustrative)' },
  { key: 'eng_kg_W', label: 'Installed engine mass per unit power', unit: 'kg/W', default: 3.5e-4, min: 0, max: 0.01, group: 'Mass coefficients (illustrative)' },
  { key: 'sigma_allow', label: 'Allowable bending stress at ultimate load', unit: 'Pa', default: 3.1e8, min: 2e7, max: 2e9, group: 'Structure', help: 'Working stress after buckling and fatigue knock-downs; illustrative default is 65% of the material ultimate strength' },
  { key: 'rho_mat', label: 'Structural material density', unit: 'kg/m³', default: 2780, min: 500, max: 9000, group: 'Structure' },
  { key: 'n_ult', label: 'Ultimate load factor', unit: 'g', default: 3.75, min: 1.5, max: 15, group: 'Structure' },
  { key: 'box_h_frac', label: 'Effective box depth / section thickness', unit: '-', default: 0.85, min: 0.4, max: 1, group: 'Structure' },
  { key: 'k_no', label: 'Non-optimum factor on bending material', unit: '-', default: 1.6, min: 1, max: 5, group: 'Mass coefficients (illustrative)', help: 'Joints, ribs, cut-outs, minimum gauge' },
  { key: 'k_sec', label: 'Secondary wing mass per unit area', unit: 'kg/m²', default: 16, min: 0, max: 80, group: 'Mass coefficients (illustrative)', help: 'Leading and trailing edges, high-lift devices, controls' },
  { key: 'fuel_rho', label: 'Fuel density', unit: 'kg/m³', default: 804, min: 60, max: 900, group: 'Constraints' },
  { key: 'tank_extra_m3', label: 'Tank volume outside the wing box', unit: 'm³', default: 0, min: 0, group: 'Constraints' },
  { key: 'runway_m', label: 'Take-off field length limit', unit: 'm', default: 2500, min: 50, group: 'Constraints' },
  { key: 'tofl_ref', label: 'Reference take-off distance (0 = uncalibrated)', unit: 'm', default: 0, min: 0, group: 'Constraints', help: 'From Suite 5 when available; calibrates the field-length factor' },
  { key: 'vapp_max', label: 'Approach speed limit', unit: 'm/s', default: 72, min: 5, max: 120, group: 'Constraints' },
  { key: 'span_max', label: 'Wing span limit', unit: 'm', default: 36, min: 0.5, max: 100, group: 'Constraints', help: 'Gate / aerodrome code limit: 15, 24, 36, 52, 65 or 80 m' },
  { key: 'batt_Wh_kg', label: 'Battery pack specific energy', unit: 'Wh/kg', default: 190, min: 50, max: 800, group: 'Electric' },
  { key: 'batt_ref_kWh', label: 'Reference battery energy', unit: 'kWh', default: 0, min: 0, group: 'Electric' },
  { key: 'eta_elec', label: 'Battery-to-shaft efficiency', unit: '-', default: 0.9, min: 0.5, max: 0.99, group: 'Electric' },
  { key: 'hover_min', label: 'Hover time per mission', unit: 'min', default: 0, min: 0, max: 30, group: 'Electric' },
  { key: 'disk_area', label: 'Total lifting disk area', unit: 'm²', default: 0, min: 0, group: 'Electric' },
  { key: 'fuel_usd_kg', label: 'Fuel price', unit: 'USD/kg', default: 0.85, min: 0, group: 'Cost hook' },
  { key: 'elec_usd_kWh', label: 'Electricity price', unit: 'USD/kWh', default: 0.14, min: 0, group: 'Cost hook' },
  { key: 'c_oew', label: 'Ownership and maintenance cost per kg of empty mass per flight', unit: 'USD/kg', default: 0.15, min: 0, group: 'Cost hook', help: 'Illustrative: about 10% of price per year divided by flights per year and empty mass' },
];
const BOUNDS = [
  { key: 'objective', label: 'Objective', type: 'select', options: ['fuel', 'mtow', 'doc'], default: 'fuel', group: 'Optimisation', help: 'fuel = mission fuel with reserves (battery mass if electric); mtow = take-off mass; doc = energy cost plus empty-mass cost per flight' },
  { key: 'free_vars', label: 'Design variables', type: 'select', options: ['AR, S', 'AR, S, altitude', 'AR, S, sweep, t/c', 'AR, S, sweep, t/c, altitude'], default: 'AR, S, sweep, t/c, altitude', group: 'Optimisation' },
  { key: 'AR_min', label: 'Aspect ratio lower bound', unit: '-', default: 6, min: 2, max: 30, group: 'Bounds' },
  { key: 'AR_max', label: 'Aspect ratio upper bound', unit: '-', default: 14, min: 3, max: 40, group: 'Bounds' },
  { key: 'S_band', label: 'Wing area band (± fraction of baseline)', unit: '-', default: 0.35, min: 0.02, max: 0.8, group: 'Bounds' },
  { key: 'sweep_max', label: 'Sweep upper bound', unit: 'deg', default: 38, min: 0, max: 55, group: 'Bounds' },
  { key: 'tc_min', label: 't/c lower bound', unit: '-', default: 0.09, min: 0.04, max: 0.2, group: 'Bounds' },
  { key: 'tc_max', label: 't/c upper bound', unit: '-', default: 0.16, min: 0.06, max: 0.25, group: 'Bounds' },
  { key: 'dh_m', label: 'Cruise altitude band (±)', unit: 'm', default: 1500, min: 0, max: 5000, group: 'Bounds' },
];
const NUM = [{ key: 'nSta', label: 'Spanwise stations in wing sizing', unit: '', default: 12, min: 4, max: 400, step: 1, discrete: true, group: 'Numerics' }];
const gateSpan = (b) => [15, 24, 36, 52, 65, 80].find((v) => v >= b * 1.001) || 1.2 * b;
const modelDefaults = (c, up, d) => {
  const m = c.mass, W = m.mtow_kg * G0, S = c.wing.S_m2 || undefined, at = isa(c.mission.cruise_alt_m || c.atm.alt_m, 0), V = c.mission.cruise_V_ms || c.flight.V_ms, elec = c.prop.type === 'electric';
  const mat = METALS[c.struct.material] || METALS['Al 2024-T3'], pack = 0.75 * (BATTERIES[c.systems.batt_chem]?.wh_kg || 250), b = c.wing.b_m;
  const vapp = S ? 1.3 * Math.sqrt((2 * (m.mtow_kg - 0.8 * m.fuel_kg) * G0) / (RHO0 * S * c.aero.CLmax_land)) : undefined, range = c.mission.range_km * 1e3;
  const ptype = ['turbofan', 'turbojet', 'turboprop', 'piston', 'electric'].includes(c.prop.type) ? c.prop.type : 'turboprop', CLc = S ? W / (0.5 * at.rho * V * V * S) : 0.5;
  return {
    AR0: d.AR || undefined, S0: S, sweep0: c.wing.sweep_deg, tc0: c.wing.tc, taper: N.clamp(c.wing.taper, 0.1, 1), h0: c.mission.cruise_alt_m || c.atm.alt_m, mach: N.clamp(c.mission.cruise_mach || V / at.a, 0.03, 0.9),
    range_km: c.mission.range_km || undefined, payload: m.payload_kg, allow_frac: elec ? 0.95 : 0.97, reserve_frac: range > 0 ? N.clamp(0.05 + (c.mission.reserve_min * 60 * V + (elec ? 0 : c.mission.alternate_km * 1e3)) / range, 0.05, 0.6) : undefined,
    mtow_ref: m.mtow_kg, oew_ref: m.oew_kg, wing_mass_ref: up.fea?.wing_struct_mass_kg > 0 && up.fea.wing_struct_mass_kg < 0.3 * m.mtow_kg ? up.fea.wing_struct_mass_kg : undefined,
    CD0_ref: up.cfd?.CD0 ?? c.aero.CD0, e_ref: N.clamp(up.cfd?.e_oswald ?? c.aero.e, 0.4, 0.98), kappa: c.aero.Mmo > 0.7 ? 0.95 : 0.87, CL_buffet: Math.max(0.7, 1.15 * CLc, c.aero.Mmo > 0.7 ? 0 : 0.7 * c.aero.CLmax_clean), CLmax_to: c.aero.CLmax_to, CLmax_land: c.aero.CLmax_land,
    ptype, tw: d.T_total / W, pw: d.P_total / m.mtow_kg, bpr: c.prop.bpr, tsfc: up.propulsion?.tsfc_kg_Ns ?? c.prop.tsfc_kg_Ns, bsfc: c.prop.bsfc_kg_Ws, eta_prop: c.prop.eta_prop,
    eng_kg_W: ptype === 'piston' ? 1.2e-3 : elec ? 3e-4 : 3.5e-4, sigma_allow: 0.65 * mat.Su, rho_mat: mat.rho, n_ult: c.aero.n_pos * c.struct.sf_ultimate, box_h_frac: c.struct.box_height_frac,
    k_sec: S ? Math.max(4, (0.025 * m.mtow_kg) / S) : undefined, fuel_rho: c.prop.fuel === 'Avgas 100LL' ? 720 : 804, tank_extra_m3: S && !elec ? Math.max(0, (1.05 * m.fuel_kg) / (c.prop.fuel === 'Avgas 100LL' ? 720 : 804) - (0.54 * S * S * c.wing.tc * (1 + c.wing.taper + c.wing.taper ** 2)) / (b * (1 + c.wing.taper) ** 2)) : 0,
    runway_m: c.site.runway_len_m, tofl_ref: up.performance?.tofl_m > 0 && Number.isFinite(up.performance.tofl_m) ? up.performance.tofl_m : undefined, vapp_max: vapp ? Math.ceil(1.03 * vapp) : undefined, span_max: b ? gateSpan(b) : undefined,
    batt_Wh_kg: pack, batt_ref_kWh: elec ? c.systems.batt_kWh : 0, hover_min: elec ? c.mission.hover_min : 0, disk_area: elec && c.rotor.R_m > 0 ? c.prop.n_eng * Math.PI * c.rotor.R_m ** 2 : 0,
    fuel_usd_kg: c.econ.fuel_usd_kg, elec_usd_kWh: c.econ.elec_usd_kWh, c_oew: c.econ.price_usd > 0 && c.econ.cycles_yr > 0 && m.oew_kg > 0 ? (0.1 * c.econ.price_usd) / c.econ.cycles_yr / m.oew_kg : undefined,
    AR_min: d.AR ? Math.max(3, Math.floor(0.6 * d.AR)) : undefined, AR_max: d.AR ? Math.ceil(1.5 * d.AR) : undefined, sweep_max: Math.max(5, c.wing.sweep_deg + 12), tc_min: Math.max(0.06, c.wing.tc - 0.035), tc_max: c.wing.tc + 0.04,
    dh_m: Math.min(1500, 0.5 * (c.mission.cruise_alt_m || 1000)), free_vars: ptype === 'turbofan' || ptype === 'turbojet' ? 'AR, S, sweep, t/c, altitude' : 'AR, S, altitude',
  };
};
const fixedWing = (c) => (c.wing.S_m2 > 0 ? true : 'This analysis sizes a wing; the current case is a pure rotorcraft. Use the rotor sizing optimisation instead.');
/** Design-variable table for a chosen set, with bounds. */
function varsOf(i, set) {
  const all = { AR: ['AR', 'Aspect ratio', '-', Math.min(i.AR_min, i.AR0), Math.max(i.AR_max, i.AR0)], S: ['S', 'Wing area', 'm²', i.S0 * (1 - i.S_band), i.S0 * (1 + i.S_band)], sweep: ['sweep', 'Sweep', 'deg', 0, Math.max(i.sweep_max, i.sweep0)], 'tc': ['tc', 'Thickness ratio', '-', Math.min(i.tc_min, i.tc0), Math.max(i.tc_max, i.tc0)], h: ['h', 'Cruise altitude', 'm', Math.max(0, i.h0 - i.dh_m), i.h0 + i.dh_m] };
  const keys = (set || i.free_vars).split(',').map((s) => s.trim()).map((s) => (s === 't/c' ? 'tc' : s === 'altitude' ? 'h' : s));
  return keys.map((k) => all[k]).filter((v) => v && v[4] > v[3]).map(([key, label, unit, lo, hi]) => ({ key, label, unit, lo, hi }));
}
const toX = (P, vars, u) => { const x = { ...P.x0 }; vars.forEach((v, j) => (x[v.key] = v.lo + (v.hi - v.lo) * u[j])); return x; };
const toU = (P, vars) => vars.map((v) => N.clamp((P.x0[v.key] - v.lo) / (v.hi - v.lo), 0, 1));
const ILLUS = 'Mass and cost coefficients are illustrative defaults; calibrate them to a known aircraft of the same class before trusting absolute numbers';
const MODELS = ['Wetted-area drag build-up (flat-plate friction with form factor)', 'Oswald efficiency trend with aspect ratio and sweep', 'Korn equation with Lock fourth-power wave drag', 'Wing bending-material sizing with inertia relief', 'Rubber engine at constant thrust or power loading', 'Breguet range with allowances and reserves'];
const ASSUME = ['Non-wing drag area, fixed equipment mass and tail sizes do not change with the wing', 'Elliptic span loading; ultimate manoeuvre case sizes the wing box', 'Constant cruise Mach and specific fuel consumption', 'Field length from a ground-roll correlation with a calibrated factor', ILLUS];

// ---- analyses -------------------------------------------------------------------------------
const mda = {
  id: 'mda', title: 'Multidisciplinary analysis: mass–fuel–wing coupling', fidelity: 'reduced-order',
  summary: 'Solves the loop in which take-off mass sets wing loads and fuel burn, which set wing and fuel mass, which set take-off mass again. Compares Gauss–Seidel fixed-point iteration with Newton’s method and reports the coupling strengths.',
  equations: ['Multidisciplinary coupling equations', 'Newton optimisation equations'],
  applicable: fixedWing,
  inputs: [...MODEL, ...NUM],
  defaults: modelDefaults,
  run(i) {
    const P = makeModel(i), c = couple(P, P.x0), gs = solveMDA(c.seidel, c.y0, { tol: 1e-10 }), jc = solveMDA(c.jacobi, c.y0, { tol: 1e-10, method: 'jacobi' }), nw = solveMDA(c.jacobi, c.y0, { method: 'newton', tol: 1e-10 });
    const warnings = [...P.warn];
    if (!gs.converged) { warnings.push('The mass loop diverges: with these inputs every extra kilogram needs more than a kilogram of wing, engine and fuel. Reduce range or payload, or improve L/D or fuel consumption.'); return { kpis: [{ key: 'mda_converged', label: 'MDA converged', value: 0, unit: '', status: 'bad' }], plots: [{ type: 'line', title: 'Coupling residual history', xlabel: 'Iteration [-]', ylabel: 'Relative residual [-]', ylog: true, series: [{ name: 'Gauss–Seidel', x: gs.hist.map((_, k) => k + 1), y: gs.hist.map((v) => (Number.isFinite(v) ? Math.max(v, 1e-16) : 1e16)) }] }], warnings, models: MODELS, assumptions: ASSUME }; }
    const o = c.report(gs.y), J = N.zeros(3), yb = gs.y;
    for (let j = 0; j < 3; j++) { const h = 1e-5 * yb[j], a = yb.slice(), b = yb.slice(); a[j] += h; b[j] -= h; const ga = c.jacobi(a), gb = c.jacobi(b); for (let k = 0; k < 3; k++) J[k][j] = (ga[k] - gb[k]) / (2 * h); }
    const rate = gs.hist.length > 4 ? gs.hist[gs.hist.length - 2] / gs.hist[gs.hist.length - 3] : NaN, rho = Math.max(...N.eig(J).map((v) => Math.hypot(v[0], v[1]))), gmin = Math.min(...o.g), fuelGap = P.electric ? (o.kWh / Math.max(i.batt_ref_kWh, 1e-9) - 1) * 100 : (o.energy_kg / Math.max(P.fuelRef, 1e-9) - 1) * 100;
    const it = (h) => h.map((_, k) => k + 1), fl = (h) => h.map((v) => Math.max(v, 1e-16));
    if (Math.abs(o.mtow / i.mtow_ref - 1) > 0.15) warnings.push(`The converged take-off mass differs from the reference by ${(100 * (o.mtow / i.mtow_ref - 1)).toFixed(0)}%: the design range and payload are not consistent with the reference masses under this model, or the drag and fuel-consumption inputs need calibration.`);
    o.g.forEach((g, k) => { if (g < 0) warnings.push(`Baseline violates the constraint "${CON_NAMES[k]}" by ${(-100 * g).toFixed(1)}%.`); });
    return {
      kpis: [
        { key: 'mtow_kg', label: 'Converged take-off mass', value: o.mtow, unit: 'kg' },
        { key: 'mtow_vs_ref_pct', label: 'Difference from reference take-off mass', value: 100 * (o.mtow / i.mtow_ref - 1), unit: '%', status: Math.abs(o.mtow / i.mtow_ref - 1) < 0.15 ? 'ok' : 'warn' },
        { key: 'wing_mass_kg', label: 'Wing mass', value: o.wing, unit: 'kg', note: `${(100 * o.wing / o.mtow).toFixed(1)}% of take-off mass` },
        { key: 'energy_mass_kg', label: P.electric ? 'Battery mass for the mission' : 'Fuel incl. reserves', value: o.energy_kg, unit: 'kg' },
        { key: 'oew_kg', label: 'Operating empty mass', value: o.oew, unit: 'kg' },
        { key: 'LD_cruise', label: 'Mid-cruise lift-to-drag ratio', value: o.LD, unit: '-' },
        { key: 'CL_cruise', label: 'Start-of-cruise lift coefficient', value: o.CL, unit: '-' },
        { key: 'CD_wave', label: 'Wave drag coefficient', value: o.CDw, unit: '-' },
        { key: 'gs_iterations', label: 'Gauss–Seidel iterations', value: gs.iters, unit: '' },
        { key: 'newton_iterations', label: 'Newton iterations', value: nw.iters, unit: '', status: nw.converged ? 'ok' : 'warn' },
        { key: 'coupling_spectral_radius', label: 'Spectral radius of the coupling Jacobian', value: rho, unit: '-', status: rho < 0.7 ? 'ok' : rho < 1 ? 'warn' : 'bad', note: 'Jacobi iteration converges only below 1; the closer to 1, the larger the mass snowball' },
        { key: 'growth_factor', label: 'Mass growth factor', value: 1 / Math.max(1e-9, 1 - J[0][0] - J[0][1] * J[1][0] - J[0][2] * J[2][0]), unit: 'kg/kg', note: 'Take-off mass added per kg of extra fixed mass (first-order)' },
        { key: 'min_constraint_margin_pct', label: 'Smallest constraint margin', value: 100 * gmin, unit: '%', status: gmin >= 0 ? 'ok' : 'bad', note: CON_NAMES[N.argmin(o.g)] },
      ],
      plots: [
        { type: 'line', title: 'Coupling residual history', xlabel: 'Iteration [-]', ylabel: 'Relative residual [-]', ylog: true, series: [{ name: 'Gauss–Seidel', x: it(gs.hist), y: fl(gs.hist), style: 'line+points' }, { name: 'Jacobi', x: it(jc.hist), y: fl(jc.hist), style: 'dash' }, { name: 'Newton', x: it(nw.hist), y: fl(nw.hist), style: 'line+points' }] },
        { type: 'bar', title: 'Take-off mass breakdown', ylabel: 'Mass [kg]', categories: ['Reference', 'Converged'], stacked: true, series: [{ name: 'Fixed', y: [P.mFixed, P.mFixed] }, { name: 'Wing', y: [P.wing0, o.wing] }, { name: 'Engines', y: [P.eng0, o.eng] }, { name: 'Payload', y: [i.payload, i.payload] }, { name: P.electric ? 'Battery' : 'Fuel', y: [P.fuelRef, o.energy_kg] }] },
      ],
      tables: [
        { title: 'Coupling structure (N2 view)', columns: ['Discipline', 'Reads', 'Returns'], rows: [['Aerodynamics', 'Wing geometry, cruise weight (take-off mass, fuel)', 'L/D, CD0, wave drag, CL'], ['Structures', 'Take-off mass, wing mass and wing fuel (inertia relief)', 'Wing mass'], ['Propulsion', 'Take-off mass (constant thrust or power loading)', 'Engine mass, cruise thrust available'], ['Mission', 'Take-off mass, L/D, specific fuel consumption', 'Mission and reserve fuel or battery mass'], ['Mass roll-up', 'Wing, engine, fuel, payload, fixed mass', 'Take-off mass']] },
        { title: 'Coupling Jacobian ∂(new)/∂(old) at the solution', columns: ['Output \\ input', 'Take-off mass', 'Wing mass', 'Energy mass'], rows: ['Take-off mass', 'Wing mass', 'Energy mass'].map((n, k) => [n, ...J[k]]) },
        { title: 'Constraint margins at the baseline', columns: ['Constraint', 'Margin [%]', 'Status'], rows: CON_NAMES.map((n, k) => [n, 100 * o.g[k], o.g[k] >= 0 ? 'met' : 'VIOLATED']) },
      ],
      outputs: { tofl_model_m: o.tofl, vapp_model_ms: o.vapp, span_m: o.b, fuel_vs_ref_pct: fuelGap, gs_rate: rate },
      warnings, models: [...MODELS, 'Gauss–Seidel and Jacobi fixed-point iteration', 'Newton iteration with finite-difference Jacobian'], assumptions: ASSUME,
    };
  },
  convergence: { param: 'nSta', label: 'Spanwise stations', levels: [6, 12, 24, 48, 96], metric: 'wing_mass_kg' },
  verify() {
    // linear contraction y = A y + b has the exact solution (I − A)⁻¹ b
    const A = [[0.2, 0.3], [0.1, 0.4]], b = [1, 2], G = (y) => [A[0][0] * y[0] + A[0][1] * y[1] + b[0], A[1][0] * y[0] + A[1][1] * y[1] + b[1]], ex = N.solve([[0.8, -0.3], [-0.1, 0.6]], b);
    const fp = solveMDA(G, [0, 0], { tol: 1e-13, method: 'jacobi' }), nw = solveMDA(G, [0, 0], { method: 'newton', tol: 1e-13 });
    const P = { n_ult: 1, taper: 0.5, nSta: 200, sigma_allow: 3e8, rho_mat: 2800, box_h_frac: 0.85, kNo: 1, k_sec: 0 }, x = { AR: 8, S: 50, sweep: 0, tc: 0.12 }, w = wingMass(P, x, 1000, 0);
    const i = Object.fromEntries(MODEL.concat(NUM).map((f) => [f.key, f.default])), o = N.kv(mda.run(i)), r = mda.run(i);
    return [
      N.check('Fixed-point iteration reaches (I − A)⁻¹b', fp.y[0], ex[0], 1e-10, 'Linear system'),
      N.check('Newton reaches the same point in one step', nw.y[1], ex[1], 1e-9, 'Linear system'),
      N.check('Root bending moment of an elliptic load = (L/2)·(4/3π)·(b/2)', w.Mroot, ((1000 * G0) / 2) * (4 / (3 * Math.PI)) * 10, 2e-4, 'Centroid of a quarter ellipse'),
      N.check('Mass roll-up closes at the converged point', r.plots[1].series.reduce((s, q) => s + q.y[1], 0), o.mtow_kg, 1e-7, 'Sum of components equals take-off mass'),
      N.check('Korn drag-divergence Mach for κ=0.95, t/c=0.115, sweep 25°, CL=0.6', aeroAt({ mach: 0.78, kappa: 0.95, e_ref: 0.8 }, { AR: 9.5, S: 100, sweep: 25, tc: 0.115 }, isa(10668), 0.6 * 0.5 * isa(10668).rho * (0.78 * isa(10668).a) ** 2 * 100).Mdd, 0.95 / Math.cos(N.rad(25)) - 0.115 / Math.cos(N.rad(25)) ** 2 - 0.06 / Math.cos(N.rad(25)) ** 3, 1e-12, 'Korn equation with simple sweep theory (Mason)'),
    ];
  },
  calibration: { params: [{ key: 'k_no', min: 1, max: 5 }, { key: 'CD0_ref', min: 0.005, max: 0.1 }], sweep: 'range_km', target: 'mtow_kg', note: 'Supply known take-off masses of the same airframe family at several design ranges.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!o.mtow_kg) return [{ severity: 'critical', title: 'The design does not close', detail: 'The mass iteration diverges.', action: 'Shorten the design range, cut payload, or improve L/D and specific fuel consumption until the loop converges.', basis: 'Fixed-point divergence of the sizing loop' }];
    out.push({ severity: o.growth_factor > 3 ? 'advise' : 'info', title: `Every kilogram saved removes about ${o.growth_factor.toFixed(1)} kg from take-off mass`, detail: `Coupling spectral radius ${o.coupling_spectral_radius.toFixed(2)}.`, action: 'Prioritise fixed-mass reductions early: they cascade through wing, engine and fuel. Each kilogram of fuel avoided is 3.16 kg of CO₂ per flight.', basis: 'Mass growth factor from the coupling Jacobian' });
    if (o.min_constraint_margin_pct < 0) out.push({ severity: 'warn', title: 'The baseline violates a design constraint', detail: res.kpis.find((k) => k.key === 'min_constraint_margin_pct').note + ` is exceeded by ${(-o.min_constraint_margin_pct).toFixed(1)}%.`, action: 'Run the constrained optimisation to find the nearest feasible design, or relax the limit if it is not a hard requirement.', basis: 'Constraint margins' });
    if (o.CD_wave > 0.001) out.push({ severity: 'advise', title: 'Wave drag is significant at the cruise point', detail: `CD wave = ${o.CD_wave.toFixed(4)}.`, action: 'More sweep, a thinner or supercritical section, or a slightly lower cruise Mach removes it; let the optimiser trade these against wing mass.', basis: 'Korn equation' });
    return out;
  },
};

/** Shared optimisation driver used by the single-objective, robust and Pareto analyses. */
function designProblem(i, P, vars, objKey, skip = []) {
  const base = evalDesign(P, P.x0), f0 = base.ok ? OBJ[objKey].f(base, P) : 1, cache = new Map();
  const fn = (u) => { const key = u.map((v) => v.toFixed(10)).join(','); if (cache.has(key)) return cache.get(key); const o = evalDesign(P, toX(P, vars, u)), r = o.ok ? { f: OBJ[objKey].f(o, P) / f0, g: o.g.filter((_, k) => !skip.includes(k)), ok: true, o } : { f: 1e3, g: CON_NAMES.filter((_, k) => !skip.includes(k)).map(() => -1), ok: false }; if (cache.size < 20000) cache.set(key, r); return r; };
  return { base, f0, fn };
}
const opt = {
  id: 'opt', title: 'Constrained design optimisation', fidelity: 'reduced-order',
  summary: 'Minimises fuel, take-off mass or a cost proxy over wing aspect ratio, area, sweep, thickness and cruise altitude, subject to field length, approach speed, span, buffet, thrust and fuel-volume limits. Reports which constraints bind and what relaxing each one would be worth.',
  equations: ['Constrained nonlinear optimisation equations', 'Lagrange multiplier equations', 'Karush–Kuhn–Tucker optimality conditions', 'Gradient-based optimisation equations', 'Multidisciplinary coupling equations'],
  applicable: fixedWing,
  inputs: [...BOUNDS, ...MODEL, ...NUM],
  defaults: modelDefaults,
  run(i, ctx) {
    const P = makeModel(i), vars = varsOf(i), warnings = [...P.warn], pr = designProblem(i, P, vars, i.objective), base = pr.base;
    if (!base.ok) return { kpis: [{ key: 'opt_converged', label: 'Optimisation ran', value: 0, unit: '', status: 'bad' }], tables: [{ title: 'Status', columns: ['Message'], rows: [['The baseline mass loop diverges; fix the baseline in the MDA analysis first.']] }], warnings: [...warnings, 'The baseline design does not close, so no optimisation was attempted.'], models: MODELS, assumptions: ASSUME };
    ctx?.progress?.(0.1, 'Augmented-Lagrangian search');
    const r = alOptimise(pr.fn, toU(P, vars)), best = pr.fn(r.u), o = best.o || base, k = kkt(pr.fn, r.u), x = toX(P, vars, r.u), feas = Math.min(...best.g) > -2e-3;
    const fObj = OBJ[i.objective], fB = fObj.f(base, P), fO = best.ok ? fObj.f(o, P) : NaN, baseFeas = Math.min(...base.g) >= -1e-9;
    if (!feas) warnings.push('No design inside the bounds satisfies every constraint; the result is the least-infeasible compromise. Widen the bounds or relax a limit.');
    if (!baseFeas) warnings.push('The baseline violates at least one constraint, so the improvement figure compares an infeasible baseline with a feasible optimum and can be negative.');
    if (k.residual > 0.05 && feas) warnings.push('The KKT stationarity residual is above 5%: the point is close to, but not exactly at, a local optimum (derivative-free search tolerance).');
    vars.forEach((v, j) => { if (r.u[j] < 1e-3 || r.u[j] > 1 - 1e-3) warnings.push(`${v.label} sits on its ${r.u[j] < 0.5 ? 'lower' : 'upper'} bound (${x[v.key].toPrecision(4)} ${v.unit}); the bound, not the physics, limits the result.`); });
    warnings.push(`${ILLUS}.`);
    const row = (n, a, b, u) => [n, a, b, u, a ? (100 * (b - a)) / Math.abs(a) : 0];
    return {
      kpis: [
        { key: 'opt_mtow_kg', label: 'Optimised take-off mass', value: o.mtow, unit: 'kg' },
        { key: 'opt_AR', label: 'Optimised aspect ratio', value: x.AR, unit: '-' },
        { key: 'opt_S_m2', label: 'Optimised wing area', value: x.S, unit: 'm²' },
        { key: 'opt_fuel_kg', label: 'Optimised fuel incl. reserves', value: P.electric ? 0 : o.energy_kg, unit: 'kg', note: P.electric ? 'Electric: see battery mass' : '' },
        { key: 'opt_improvement_pct', label: `Improvement in ${i.objective}`, value: 100 * (1 - fO / fB), unit: '%', status: feas ? 'ok' : 'warn' },
        { key: 'opt_objective', label: `Objective: ${fObj.label}`, value: fO, unit: fObj.unit },
        { key: 'opt_sweep_deg', label: 'Optimised sweep', value: x.sweep, unit: 'deg' },
        { key: 'opt_tc', label: 'Optimised thickness ratio', value: x.tc, unit: '-' },
        { key: 'opt_cruise_alt_m', label: 'Optimised cruise altitude', value: x.h, unit: 'm' },
        { key: 'opt_span_m', label: 'Optimised span', value: o.b, unit: 'm' },
        { key: 'opt_LD', label: 'Optimised cruise L/D', value: o.LD, unit: '-' },
        { key: 'opt_battery_kg', label: 'Optimised battery mass', value: P.electric ? o.energy_kg : 0, unit: 'kg' },
        { key: 'kkt_residual', label: 'KKT stationarity residual', value: k.residual, unit: '-', status: k.residual < 0.05 ? 'ok' : 'warn', note: 'Norm of ∇f − Σλ∇g on free variables, relative to ∇f' },
        { key: 'n_active', label: 'Active constraints', value: best.g.filter((g) => g < 3e-3).length, unit: '' },
        { key: 'opt_feasible', label: 'Optimum is feasible', value: feas ? 1 : 0, unit: '', status: feas ? 'ok' : 'bad' },
        { key: 'model_evaluations', label: 'Design evaluations', value: r.evals, unit: '' },
      ],
      plots: [
        { type: 'line', title: 'Optimisation history', xlabel: 'Design evaluations [-]', ylabel: 'Objective relative to baseline [-]', series: [{ name: 'Objective', x: r.hist.map((h) => h.evals), y: r.hist.map((h) => Math.min(h.f, 3)), style: 'line+points' }] },
        { type: 'line', title: 'Constraint violation history', xlabel: 'Design evaluations [-]', ylabel: 'Largest violation [-]', series: [{ name: 'Max violation', x: r.hist.map((h) => h.evals), y: r.hist.map((h) => h.viol), style: 'line+points' }] },
        { type: 'bar', title: 'Constraint margins before and after', ylabel: 'Margin [%]', categories: CON_NAMES, series: [{ name: 'Baseline', y: base.g.map((g) => 100 * g) }, { name: 'Optimised', y: best.g.map((g) => 100 * g) }] },
      ],
      tables: [
        { title: 'Baseline versus optimised design', columns: ['Quantity', 'Baseline', 'Optimised', 'Unit', 'Change [%]'], rows: [...vars.map((v) => row(v.label, P.x0[v.key], x[v.key], v.unit)), row('Take-off mass', base.mtow, o.mtow, 'kg'), row('Wing mass', base.wing, o.wing, 'kg'), row(P.electric ? 'Battery mass' : 'Fuel incl. reserves', base.energy_kg, o.energy_kg, 'kg'), row('Operating empty mass', base.oew, o.oew, 'kg'), row('Cruise L/D', base.LD, o.LD, '-'), row('Span', base.b, o.b, 'm'), row('Take-off field length', base.tofl, o.tofl, 'm'), row('Approach speed', base.vapp, o.vapp, 'm/s'), row(fObj.label, fB, fO, fObj.unit)] },
        { title: 'Constraints and Lagrange multipliers', columns: ['Constraint', 'Margin [%]', 'Active', 'Multiplier λ', 'Meaning'], rows: CON_NAMES.map((n, c) => [n, 100 * best.g[c], best.g[c] < 3e-3 ? 'yes' : '', k.lambda[c], k.lambda[c] > 1e-4 ? `Relaxing this limit by 1% improves the objective by about ${k.lambda[c].toFixed(2)}%` : '']) },
      ],
      outputs: { opt_x: vars.map((v) => x[v.key]), opt_vars: vars.map((v) => v.key).join(',') },
      warnings, models: [...MODELS, 'Augmented Lagrangian (Powell–Hestenes–Rockafellar) with Nelder–Mead inner solves', 'Projected-gradient polish with central differences', 'Least-squares Lagrange multipliers on the active set'], assumptions: [...ASSUME, 'Local optimum from a single start at the baseline design'],
    };
  },
  verify() {
    // min (2u₁)² + (2u₂)² subject to 2u₁ + 2u₂ ≥ 1 → u = (0.25, 0.25), f = 0.5, λ = 1
    const fn = (u) => ({ f: 4 * u[0] * u[0] + 4 * u[1] * u[1], g: [2 * u[0] + 2 * u[1] - 1] }), r = alOptimise(fn, [0.8, 0.6]), k = kkt(fn, r.u);
    const f2 = (u) => ({ f: (u[0] - 0.3) ** 2 + (u[1] - 2) ** 2, g: [1] }), r2 = alOptimise(f2, [0.5, 0.5]);
    return [
      N.check('Constrained quadratic optimum f* = 0.5', r.f, 0.5, 1e-4, 'KKT solution by hand'),
      N.check('Optimiser location u₁ = 0.25', r.u[0], 0.25, 2e-3, 'KKT solution by hand'),
      N.check('Lagrange multiplier λ = 1', k.lambda[0], 1, 5e-3, '∇f = λ∇g'),
      N.check('Bound-constrained optimum lands on the box face u₂ = 1', r2.u[1], 1, 1e-9, 'Projection onto the box'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.opt_converged === 0) return [{ severity: 'critical', title: 'Baseline does not close', detail: 'No optimisation was possible.', action: 'Fix the baseline with the MDA analysis first.', basis: 'Mass-loop divergence' }];
    const act = res.tables[1].rows.filter((r) => r[3] > 1e-4).sort((a, b) => b[3] - a[3]);
    if (o.opt_feasible && o.opt_improvement_pct > 0.5) out.push({ severity: 'advise', title: `${o.opt_improvement_pct.toFixed(1)}% lower ${i.objective} is available`, detail: `Aspect ratio ${o.opt_AR.toFixed(2)}, wing area ${o.opt_S_m2.toFixed(1)} m², sweep ${o.opt_sweep_deg.toFixed(1)}°, t/c ${o.opt_tc.toFixed(3)}, cruise at ${o.opt_cruise_alt_m.toFixed(0)} m.`, action: 'Carry the optimised geometry into the case and re-run Suites 1, 2, 3 and 5 to confirm drag, wing mass, flutter margin and field performance with the higher-fidelity models. Fuel saved scales directly to CO₂ and energy cost.', basis: 'Constrained optimum of the low-order multidisciplinary model' });
    else if (o.opt_feasible) out.push({ severity: 'info', title: 'The baseline is already close to the optimum of this model', detail: `Improvement ${o.opt_improvement_pct.toFixed(2)}%.`, action: 'Look for gains outside the wing planform: engine efficiency, fixed mass, or the mission requirements themselves.', basis: 'Constrained optimum' });
    if (act.length) out.push({ severity: 'advise', title: `${act[0][0]} is the most valuable constraint to relax`, detail: act.map((r) => `${r[0]}: λ = ${Number(r[3]).toFixed(2)}`).join('; ') + '.', action: `Each 1% of relief on "${act[0][0]}" is worth about ${Number(act[0][3]).toFixed(2)}% of the objective. Question the requirement behind it before spending design effort elsewhere.`, basis: 'Lagrange multipliers as shadow prices' });
    if (!o.opt_feasible) out.push({ severity: 'warn', title: 'No feasible design inside the bounds', detail: 'At least one constraint stays violated.', action: 'Widen the variable bounds, add high-lift capability (CLmax) or thrust, or relax the binding requirement.', basis: 'Constraint violation at the returned point' });
    return out;
  },
};

const pareto = {
  id: 'pareto', title: 'Multi-objective trade-off: Pareto front by NSGA-II', fidelity: 'reduced-order',
  summary: 'Evolves a population of designs towards the set where one objective cannot improve without worsening the other, and marks the knee point that gives the best balance.',
  equations: ['Pareto optimality relations', 'Constrained nonlinear optimisation equations', 'Multidisciplinary coupling equations'],
  applicable: fixedWing,
  inputs: [
    { key: 'pair', label: 'Objective pair', type: 'select', options: ['fuel vs take-off mass', 'fuel vs field length', 'CO2 vs cost'], default: 'fuel vs take-off mass', group: 'Optimisation', help: 'Field length is a proxy for airport access and, through lower take-off thrust, community noise' },
    ...BOUNDS.filter((f) => f.key !== 'objective'), ...MODEL, ...NUM,
    { key: 'pop', label: 'Population size', unit: '', default: 28, min: 8, max: 200, step: 1, discrete: true, group: 'Numerics' },
    { key: 'gens', label: 'Generations', unit: '', default: 20, min: 2, max: 300, step: 1, discrete: true, group: 'Numerics' },
    { key: 'seed', label: 'Random seed', unit: '', default: 23, min: 1, max: 1e9, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: modelDefaults,
  run(i) {
    const P = makeModel(i), vars = varsOf(i), warnings = [...P.warn], base = evalDesign(P, P.x0); let pair = i.pair;
    if (!base.ok) return { kpis: [{ key: 'n_pareto', label: 'Pareto designs', value: 0, unit: '', status: 'bad' }], tables: [{ title: 'Status', columns: ['Message'], rows: [['The baseline mass loop diverges.']] }], warnings: [...warnings, 'The baseline design does not close.'], models: MODELS, assumptions: ASSUME };
    if (P.electric && pair === 'CO2 vs cost') { pair = 'fuel vs take-off mass'; warnings.push('In-flight CO₂ is zero for a battery aircraft, so the battery-mass versus take-off-mass pair was used instead.'); }
    const D = { 'fuel vs take-off mass': { f: [(o) => o.energy_kg, (o) => o.mtow], lab: [P.electric ? 'Battery mass [kg]' : 'Fuel incl. reserves [kg]', 'Take-off mass [kg]'], skip: [] }, 'fuel vs field length': { f: [(o) => o.energy_kg, (o) => o.tofl], lab: [P.electric ? 'Battery mass [kg]' : 'Fuel incl. reserves [kg]', 'Take-off field length [m]'], skip: [0] }, 'CO2 vs cost': { f: [(o) => 3.16 * o.mission_kg, (o) => OBJ.doc.f(o, P)], lab: ['Mission CO₂ [kg]', 'Cost proxy per flight [USD]'], skip: [] } }[pair];
    const fn = (u) => { const o = evalDesign(P, toX(P, vars, u)); if (!o.ok) return { f: [1e30, 1e30], viol: 10 }; const g = o.g.filter((_, k) => !D.skip.includes(k)); return { f: D.f.map((q) => q(o)), viol: N.sum(g.map((v) => Math.max(0, -v))), data: o }; };
    const front = nsga2(fn, vars.length, { pop: Math.round(i.pop), gens: Math.round(i.gens), seed: i.seed }), fb = D.f.map((q) => q(base));
    if (front.length < 3) { warnings.push('Fewer than three feasible non-dominated designs were found; the objectives may not conflict inside these bounds, or the constraints leave almost no feasible space.'); }
    if (!front.length) return { kpis: [{ key: 'n_pareto', label: 'Pareto designs', value: 0, unit: '', status: 'bad' }], tables: [{ title: 'Status', columns: ['Message'], rows: [['No feasible design found inside the bounds.']] }], warnings, models: MODELS, assumptions: ASSUME };
    const f1 = front.map((a) => a.f[0]), f2 = front.map((a) => a.f[1]), s1 = N.amax(f1) - N.amin(f1) || 1, s2 = N.amax(f2) - N.amin(f2) || 1, A = front[0], B = front[front.length - 1];
    // knee: largest normalised distance below the line joining the two extreme designs
    const dist = front.map((a) => { const x = (a.f[0] - A.f[0]) / s1, y = (a.f[1] - A.f[1]) / s2, bx = (B.f[0] - A.f[0]) / s1, by = (B.f[1] - A.f[1]) / s2; return (bx * y - by * x) / (Math.hypot(bx, by) || 1); }), kn = front[N.argmin(dist)], xk = toX(P, vars, kn.x);
    const dominated = front.some((a) => a.f[0] <= fb[0] && a.f[1] <= fb[1] && (a.f[0] < fb[0] || a.f[1] < fb[1]));
    warnings.push(`${ILLUS}.`);
    return {
      kpis: [
        { key: 'n_pareto', label: 'Non-dominated feasible designs', value: front.length, unit: '' },
        { key: 'knee_f1', label: `Knee: ${D.lab[0]}`, value: kn.f[0], unit: '' }, { key: 'knee_f2', label: `Knee: ${D.lab[1]}`, value: kn.f[1], unit: '' },
        { key: 'knee_AR', label: 'Knee aspect ratio', value: xk.AR, unit: '-' }, { key: 'knee_S_m2', label: 'Knee wing area', value: xk.S, unit: 'm²' },
        { key: 'f1_range_pct', label: 'Spread of objective 1 along the front', value: (100 * s1) / Math.abs(N.amin(f1) || 1), unit: '%' }, { key: 'f2_range_pct', label: 'Spread of objective 2 along the front', value: (100 * s2) / Math.abs(N.amin(f2) || 1), unit: '%' },
        { key: 'baseline_dominated', label: 'Baseline is dominated', value: dominated ? 1 : 0, unit: '', status: dominated ? 'warn' : 'ok', note: 'A front design beats the baseline on both objectives' },
        { key: 'knee_f1_vs_base_pct', label: 'Knee versus baseline, objective 1', value: 100 * (kn.f[0] / fb[0] - 1), unit: '%' }, { key: 'knee_f2_vs_base_pct', label: 'Knee versus baseline, objective 2', value: 100 * (kn.f[1] / fb[1] - 1), unit: '%' },
      ],
      plots: [{ type: 'line', title: `Pareto front: ${pair}`, xlabel: D.lab[0], ylabel: D.lab[1], series: [{ name: 'Pareto front', x: f1, y: f2, style: 'line+points' }, { name: 'Knee point', x: [kn.f[0]], y: [kn.f[1]], style: 'points' }, { name: 'Baseline', x: [fb[0]], y: [fb[1]], style: 'points' }] },
        { type: 'line', title: 'Design variables along the front', xlabel: D.lab[0], ylabel: 'Normalised variable (0 = lower bound, 1 = upper bound) [-]', series: vars.map((v, j) => ({ name: v.label, x: f1, y: front.map((a) => a.x[j]), style: 'line+points' })) }],
      tables: [{ title: 'Pareto-optimal designs', columns: [D.lab[0], D.lab[1], ...vars.map((v) => `${v.label} [${v.unit}]`), 'Take-off mass [kg]', 'L/D [-]'], rows: front.map((a) => { const x = toX(P, vars, a.x); return [a.f[0], a.f[1], ...vars.map((v) => x[v.key]), a.data.mtow, a.data.LD]; }) }],
      warnings, models: [...MODELS, 'NSGA-II: non-dominated sorting, crowding distance, SBX crossover, polynomial mutation, constrained domination'], assumptions: [...ASSUME, 'Knee = point farthest below the chord joining the two extreme designs in normalised objectives'],
    };
  },
  verify() {
    // Schaffer problem on x = 4u − 1: f1 = x², f2 = (x − 2)²; Pareto set 0 ≤ x ≤ 2 where √f1 + √f2 = 2
    const fr = nsga2((u) => { const x = 4 * u[0] - 1; return { f: [x * x, (x - 2) ** 2], viol: 0 }; }, 1, { pop: 30, gens: 25, seed: 3 });
    const worst = Math.max(...fr.map((a) => Math.abs(Math.sqrt(a.f[0]) + Math.sqrt(a.f[1]) - 2)));
    const cons = nsga2((u) => ({ f: [u[0], 1 - u[0]], viol: Math.max(0, 0.3 - u[0]) }), 1, { pop: 20, gens: 20, seed: 4 });
    return [
      N.check('All Schaffer front points satisfy √f1 + √f2 = 2', worst + 2, 2, 1e-3, 'Analytical Pareto set (Schaffer 1985); finite-population tolerance'),
      N.check('Front spans the whole set (min f1 → 0)', fr[0].f[0] + 1, 1, 2e-3, 'Extreme point x = 0'),
      N.check('Constrained domination keeps the front feasible (min u ≥ 0.3)', Math.min(...cons.map((a) => a.x[0])) >= 0.3 ? 1 : 0, 1, 0, 'Deb constrained-domination rule'),
    ];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (!o.n_pareto) return [{ severity: 'warn', title: 'No feasible trade-off found', detail: 'The search found no feasible design.', action: 'Relax bounds or constraints and repeat.', basis: 'Empty feasible set' }];
    out.push({ severity: 'advise', title: 'Use the knee design as the balanced candidate', detail: `Knee at aspect ratio ${o.knee_AR.toFixed(2)} and wing area ${o.knee_S_m2.toFixed(1)} m²: ${o.knee_f1_vs_base_pct.toFixed(1)}% on objective 1 and ${o.knee_f2_vs_base_pct.toFixed(1)}% on objective 2 against the baseline.`, action: 'Moving away from the knee in either direction buys little of one objective for a lot of the other. Pick the end of the front only if one objective has a hard target.', basis: 'Knee point of the Pareto front' });
    if (o.baseline_dominated) out.push({ severity: 'advise', title: 'The baseline is not Pareto-optimal', detail: 'At least one design is better on both objectives at once.', action: 'Adopt a front design: this is a free improvement within the model, including lower fuel burn and CO₂.', basis: 'Pareto dominance' });
    if (o.f1_range_pct < 1 || o.f2_range_pct < 1) out.push({ severity: 'info', title: 'The two objectives barely conflict here', detail: `Spreads along the front are ${o.f1_range_pct.toFixed(2)}% and ${o.f2_range_pct.toFixed(2)}%.`, action: 'A single-objective optimisation is sufficient for this pair.', basis: 'Extent of the Pareto front' });
    return out;
  },
};

const VAR_OPTS = ['AR', 'S', 'sweep', 't/c', 'altitude'];
const surrogate = {
  id: 'surrogate', title: 'Design-space map, surrogate models and efficient global optimisation', fidelity: 'numerical',
  summary: 'Maps the objective and the feasible region over two design variables, fits a quadratic response surface and a Gaussian-process (Kriging) model to a Latin-hypercube sample, checks them by leave-one-out error, and uses expected improvement to find the optimum with few model runs.',
  equations: ['Surrogate-assisted optimisation formulations', 'Bayesian optimisation equations', 'Constrained nonlinear optimisation equations'],
  applicable: fixedWing,
  inputs: [
    { key: 'xvar', label: 'Horizontal variable', type: 'select', options: VAR_OPTS, default: 'AR', group: 'Optimisation' }, { key: 'yvar', label: 'Vertical variable', type: 'select', options: VAR_OPTS, default: 'S', group: 'Optimisation' },
    ...BOUNDS.filter((f) => f.key !== 'free_vars'), ...MODEL, ...NUM,
    { key: 'nDoe', label: 'Latin-hypercube samples', unit: '', default: 24, min: 8, max: 120, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nInfill', label: 'Expected-improvement infill points', unit: '', default: 6, min: 0, max: 40, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nGrid', label: 'Map resolution per axis', unit: '', default: 13, min: 5, max: 61, step: 1, discrete: true, group: 'Numerics' },
    { key: 'seed', label: 'Random seed', unit: '', default: 7, min: 1, max: 1e9, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: modelDefaults,
  run(i) {
    const P = makeModel(i), warnings = [...P.warn], yv = i.yvar === i.xvar ? VAR_OPTS.find((v) => v !== i.xvar) : i.yvar, vars = varsOf(i, `${i.xvar}, ${yv}`);
    if (yv !== i.yvar) warnings.push('The two map variables were identical; the vertical one was changed.');
    const pr = designProblem(i, P, vars, i.objective), fail = (m) => ({ kpis: [{ key: 'gp_loo_rmse_pct', label: 'Kriging leave-one-out error', value: 0, unit: '%', status: 'bad' }], tables: [{ title: 'Status', columns: ['Message'], rows: [[m]] }], warnings: [...warnings, m], models: MODELS, assumptions: ASSUME });
    if (!pr.base.ok) return fail('The baseline mass loop diverges.'); if (vars.length < 2) return fail('The chosen variables have no range between their bounds.');
    const ng = Math.round(i.nGrid), gu = N.linspace(0, 1, ng), gx = gu.map((u) => vars[0].lo + u * (vars[0].hi - vars[0].lo)), gy = gu.map((u) => vars[1].lo + u * (vars[1].hi - vars[1].lo));
    const Z = N.zeros(ng), Gm = N.zeros(ng); let bestG = { f: Infinity };
    for (let b = 0; b < ng; b++) for (let a = 0; a < ng; a++) { const r = pr.fn([gu[a], gu[b]]); Z[b][a] = r.ok ? r.f : NaN; Gm[b][a] = r.ok ? Math.min(...r.g) : -1; if (r.ok && Gm[b][a] >= 0 && r.f < bestG.f) bestG = { f: r.f, u: [gu[a], gu[b]] }; }
    const zmax = N.amax(Z.flat().filter(Number.isFinite)); for (const r of Z) for (let a = 0; a < ng; a++) if (!Number.isFinite(r[a])) r[a] = zmax;
    const bx = [], by = []; // feasibility boundary: sign changes of the worst constraint along grid lines
    for (let a = 0; a < ng; a++) for (let b = 1; b < ng; b++) if (Gm[b][a] * Gm[b - 1][a] < 0) { const t = Gm[b - 1][a] / (Gm[b - 1][a] - Gm[b][a]); bx.push(gx[a]); by.push(gy[b - 1] + t * (gy[b] - gy[b - 1])); }
    const ord = N.range(bx.length).sort((p, q) => bx[p] - bx[q]);
    // DOE and surrogates on the objective and on the worst constraint margin
    const u = N.rng(i.seed), X = N.lhs(Math.round(i.nDoe), 2, u), R = X.map((x) => pr.fn(x)), keep = N.range(X.length).filter((k) => R[k].ok), Xs = keep.map((k) => X[k]), F = keep.map((k) => R[k].f), Gw = keep.map((k) => Math.min(...R[k].g));
    if (Xs.length < 8) return fail('Too few sample designs converged to fit a surrogate.');
    const beta = N.lstsq(Xs.map(quadBasis), F, 1e-12), rsm = (x) => N.dot(quadBasis(x), beta), H = Xs.map(quadBasis), Hinv = N.inv(N.matmul(N.transpose(H), H).map((r, a) => r.map((v, b) => v + (a === b ? 1e-12 : 0))));
    const looR = Xs.map((x, k) => (F[k] - rsm(x)) / Math.max(1e-9, 1 - N.dot(H[k], N.matvec(Hinv, H[k])))), gp = gpFit(Xs, F), rm = (a) => Math.sqrt(N.mean(a.map((v) => v * v))), fbar = Math.abs(N.mean(F));
    let gpc = gpFit(Xs, Gw), gpo = gp; const Xa = Xs.slice(), Fa = F.slice(), Ga = Gw.slice(), infill = [], cand = []; for (let b = 0; b < 25; b++) for (let a = 0; a < 25; a++) cand.push([a / 24, b / 24]);
    for (let k = 0; k < Math.round(i.nInfill); k++) { // constrained expected improvement: EI × probability of feasibility
      const feas = Fa.filter((_, j) => Ga[j] >= 0), fmin = feas.length ? Math.min(...feas) : Math.max(...Fa); let bv = -1, bu = null;
      for (const c of cand) { if (Xa.some((p) => Math.hypot(p[0] - c[0], p[1] - c[1]) < 1e-6)) continue; const po = gpo.predict(c), pc = gpc.predict(c), v = (feas.length ? expImp(fmin, po.m, po.s) : 1) * N.normCdf(pc.m / Math.max(pc.s, 1e-9)); if (v > bv) { bv = v; bu = c; } }
      if (!bu || bv < 1e-12) break; const r = pr.fn(bu); if (!r.ok) { cand.splice(cand.indexOf(bu), 1); continue; }
      Xa.push(bu); Fa.push(r.f); Ga.push(Math.min(...r.g)); infill.push(bu); gpo = gpFit(Xa, Fa, gpo.ell); gpc = gpFit(Xa, Ga, gpc.ell);
    }
    const fe = Fa.map((f, j) => (Ga[j] >= 0 ? f : Infinity)), jb = N.argmin(fe), ego = Number.isFinite(fe[jb]) ? { f: fe[jb], u: Xa[jb] } : null, ph = (uu, v) => v.lo + uu * (v.hi - v.lo);
    if (!ego) warnings.push('No feasible design was sampled; the map shows where the constraints are violated.');
    if (rm(looR) / fbar > 0.02) warnings.push('The quadratic response surface misses more than 2% in leave-one-out error: the objective is not quadratic over this range (wave-drag rise or constraint-driven behaviour). Prefer the Kriging model.');
    return {
      kpis: [
        { key: 'gp_loo_rmse_pct', label: 'Kriging leave-one-out RMS error', value: (100 * rm(gp.loo)) / fbar, unit: '%', status: rm(gp.loo) / fbar < 0.01 ? 'ok' : 'warn' },
        { key: 'rsm_loo_rmse_pct', label: 'Quadratic surface leave-one-out RMS error', value: (100 * rm(looR)) / fbar, unit: '%', status: rm(looR) / fbar < 0.02 ? 'ok' : 'warn' },
        { key: 'ego_best_rel', label: 'Best feasible objective found by infill (relative to baseline)', value: ego ? ego.f : NaN, unit: '-' },
        { key: 'grid_best_rel', label: 'Best feasible objective on the map grid', value: Number.isFinite(bestG.f) ? bestG.f : NaN, unit: '-' },
        { key: 'ego_gap_pct', label: 'Infill optimum versus map optimum', value: ego && Number.isFinite(bestG.f) ? 100 * (ego.f / bestG.f - 1) : NaN, unit: '%', note: `${Xa.length} model runs against ${ng * ng} for the map` },
        { key: 'ego_x', label: `Best ${vars[0].label}`, value: ego ? ph(ego.u[0], vars[0]) : NaN, unit: vars[0].unit }, { key: 'ego_y', label: `Best ${vars[1].label}`, value: ego ? ph(ego.u[1], vars[1]) : NaN, unit: vars[1].unit },
        { key: 'feasible_fraction_pct', label: 'Feasible share of the design space', value: (100 * Gm.flat().filter((g) => g >= 0).length) / (ng * ng), unit: '%' },
        { key: 'gp_length_x', label: 'Kriging correlation length, horizontal', value: gp.ell[0], unit: '-', note: 'In units of the variable range; large = weak influence' }, { key: 'gp_length_y', label: 'Kriging correlation length, vertical', value: gp.ell[1], unit: '-' },
      ],
      plots: [
        { type: 'heat', title: `${OBJ[i.objective].label} relative to baseline`, xlabel: `${vars[0].label} [${vars[0].unit}]`, ylabel: `${vars[1].label} [${vars[1].unit}]`, zlabel: 'Objective / baseline [-]', x: gx, y: gy, z: Z, contours: 12, overlay: [{ name: 'Feasibility boundary', x: ord.map((k) => bx[k]), y: ord.map((k) => by[k]) }, { name: 'Latin-hypercube samples', x: Xs.map((p) => ph(p[0], vars[0])), y: Xs.map((p) => ph(p[1], vars[1])), style: 'points' }, { name: 'Infill points', x: infill.map((p) => ph(p[0], vars[0])), y: infill.map((p) => ph(p[1], vars[1])), style: 'points' }] },
        { type: 'heat', title: 'Worst constraint margin (negative = infeasible)', xlabel: `${vars[0].label} [${vars[0].unit}]`, ylabel: `${vars[1].label} [${vars[1].unit}]`, zlabel: 'Margin [-]', x: gx, y: gy, z: Gm, contours: 10, diverging: true },
        { type: 'line', title: 'Leave-one-out prediction against the model', xlabel: 'Model objective [-]', ylabel: 'Surrogate prediction [-]', series: [{ name: 'Kriging', x: F, y: F.map((f, k) => f - gp.loo[k]), style: 'points' }, { name: 'Quadratic surface', x: F, y: F.map((f, k) => f - looR[k]), style: 'points' }, { name: 'Perfect', x: [N.amin(F), N.amax(F)], y: [N.amin(F), N.amax(F)], style: 'dash' }] },
      ],
      tables: [{ title: 'Quadratic response surface (variables scaled 0…1)', columns: ['Term', 'Coefficient'], rows: ['1', 'x', 'y', 'x²', 'x·y', 'y²'].map((t, k) => [t, beta[k]]) }],
      warnings, models: [...MODELS, 'Latin hypercube design of experiments', 'Quadratic response surface by least squares', 'Ordinary Kriging with anisotropic squared-exponential kernel, likelihood grid search', 'Constrained expected improvement'], assumptions: [...ASSUME, 'Other design variables are held at the baseline', 'Kriging hyper-parameters are frozen during infill'],
    };
  },
  verify() {
    const u = N.rng(1), X = N.lhs(14, 2, u), q = (x) => 1 + 2 * x[0] - 3 * x[1] + 0.5 * x[0] * x[0] + 4 * x[0] * x[1] - x[1] * x[1], b = N.lstsq(X.map(quadBasis), X.map(q));
    const gp = gpFit(X, X.map((x) => Math.sin(3 * x[0]) + x[1] * x[1])), t = [0.37, 0.61];
    return [
      N.check('Quadratic surface recovers the x·y coefficient', b[4], 4, 1e-7, 'Exact quadratic data'),
      N.check('Kriging interpolates its training data', gp.predict(X[3]).m, Math.sin(3 * X[3][0]) + X[3][1] ** 2, 5e-4, 'Interpolation property (1e-8 nugget)'),
      N.check('Kriging predicts a smooth function between samples', gp.predict(t).m, Math.sin(3 * t[0]) + t[1] * t[1], 0.02, '14 samples of a smooth 2-D function'),
      N.check('Expected improvement at μ = f_min, σ = 1 equals φ(0)', expImp(0, 0, 1), 1 / Math.sqrt(2 * Math.PI), 1e-12, 'Jones, Schonlau & Welch (1998)'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.gp_loo_rmse_pct === 0 && res.kpis.length === 1) return [{ severity: 'warn', title: 'Surrogate study could not run', detail: res.warnings[res.warnings.length - 1], action: 'Fix the baseline or widen the variable bounds.', basis: 'Model failure' }];
    if (Number.isFinite(o.ego_best_rel)) out.push({ severity: o.ego_best_rel < 0.995 ? 'advise' : 'info', title: `Best design in this plane: ${(100 * (1 - o.ego_best_rel)).toFixed(1)}% better than baseline`, detail: `${res.kpis[5].label} = ${o.ego_x.toPrecision(4)}, ${res.kpis[6].label} = ${o.ego_y.toPrecision(4)}; found with a fraction of the runs the full map needed (gap to map optimum ${Number.isFinite(o.ego_gap_pct) ? o.ego_gap_pct.toFixed(2) : '—'}%).`, action: 'When each evaluation is a CFD or FE run, use this sample-fit-infill loop instead of a grid: the same answer for far less computing energy.', basis: 'Constrained expected improvement on a Kriging model' });
    if (o.feasible_fraction_pct < 25) out.push({ severity: 'advise', title: 'The feasible region is small', detail: `${o.feasible_fraction_pct.toFixed(0)}% of this design plane meets all constraints.`, action: 'Read the margin map to see which limit cuts the space, then review that requirement or the variable bounds.', basis: 'Feasibility map' });
    if (o.gp_length_x > 2.5 || o.gp_length_y > 2.5) out.push({ severity: 'info', title: 'One variable has little influence in this range', detail: `Correlation lengths ${o.gp_length_x.toFixed(2)} and ${o.gp_length_y.toFixed(2)} of the range.`, action: `Drop the weak variable (${o.gp_length_x > o.gp_length_y ? i.xvar : i.yvar}) from the optimisation and spend the design freedom elsewhere.`, basis: 'Kriging length scales as a screening measure' });
    return out;
  },
};

const PARAMS = [
  ['Zero-lift drag', (P, r) => { P.fRest += r * P.CD0_ref * P.S0; }], ['Specific fuel consumption', (P, r) => { P.tsfc *= 1 + r; P.bsfc *= 1 + r; if (P.electric) P.eta_elec /= 1 + r; }], ['Allowable stress', (P, r) => { P.sigma_allow *= 1 + r; }],
  ['Payload', (P, r) => { P.payload *= 1 + r; }], ['Range', (P, r) => { P.range_m *= 1 + r; }], ['Oswald efficiency', (P, r) => { P.e_ref *= 1 + r; }], ['Fixed mass', (P, r) => { P.mFixed *= 1 + r; }],
];
const sens = {
  id: 'sens', title: 'Design sensitivities: finite differences and coupled adjoint', fidelity: 'numerical',
  summary: 'Computes how the objective responds to every design variable and key technology parameter, compares finite-difference derivatives with the coupled direct and adjoint methods, and shows how the finite-difference error depends on step size.',
  equations: ['Adjoint sensitivity equations', 'Gradient-based optimisation equations', 'Multidisciplinary coupling equations'],
  applicable: fixedWing,
  inputs: [...BOUNDS.filter((f) => f.key === 'objective'), ...MODEL, ...NUM,
    { key: 'fd_step', label: 'Relative finite-difference step', unit: '-', default: 1e-3, min: 1e-8, max: 0.2, group: 'Numerics' }],
  defaults: modelDefaults,
  run(i) {
    const P = makeModel(i), warnings = [...P.warn], obj = OBJ[i.objective], keys = ['AR', 'S', 'sweep', 'tc', 'h'], labels = ['Aspect ratio', 'Wing area', 'Sweep', 'Thickness ratio', 'Cruise altitude'];
    const F = (Pm, x, tol = 1e-13) => { const c = couple(Pm, x), s = solveMDA(c.jacobi, c.y0, { method: 'newton', tol }); return s.converged ? obj.f(c.report(s.y), Pm) : NaN; }, f0 = F(P, P.x0);
    if (!Number.isFinite(f0)) return { kpis: [{ key: 'adjoint_fd_max_diff_pct', label: 'Adjoint versus finite difference', value: 0, unit: '%', status: 'bad' }], tables: [{ title: 'Status', columns: ['Message'], rows: [['The baseline mass loop diverges.']] }], warnings: [...warnings, 'The baseline design does not close.'], models: MODELS, assumptions: ASSUME };
    const scale = (k) => (k === 'sweep' ? Math.max(Math.abs(P.x0.sweep), 10) : k === 'h' ? Math.max(P.x0.h, 1000) : P.x0[k]), xp = (k, d) => ({ ...P.x0, [k]: P.x0[k] + d });
    const fd = keys.map((k) => { const d = i.fd_step * scale(k); return (F(P, xp(k, d)) - F(P, xp(k, -d))) / (2 * d); });
    // coupled derivatives from partials of R(y, x) = y − G(y, x) and of f(y, x) at the converged state
    const c0 = couple(P, P.x0), y = solveMDA(c0.jacobi, c0.y0, { method: 'newton', tol: 1e-13 }).y, fy = (Pm, x, yy) => obj.f(couple(Pm, x).report(yy), Pm), Rf = (x, yy) => { const g = couple(P, x).jacobi(yy); return yy.map((v, j) => v - g[j]); };
    const Ry = N.zeros(3), dfdy = new Array(3);
    for (let j = 0; j < 3; j++) { const h = 1e-6 * y[j], a = y.slice(), b = y.slice(); a[j] += h; b[j] -= h; const ra = Rf(P.x0, a), rb = Rf(P.x0, b); for (let k = 0; k < 3; k++) Ry[k][j] = (ra[k] - rb[k]) / (2 * h); dfdy[j] = (fy(P, P.x0, a) - fy(P, P.x0, b)) / (2 * h); }
    const psi = N.solve(N.transpose(Ry), dfdy), direct = [], adjoint = [];
    for (const k of keys) { const d = 1e-6 * scale(k), ra = Rf(xp(k, d), y), rb = Rf(xp(k, -d), y), Rx = ra.map((v, j) => (v - rb[j]) / (2 * d)), fx = (fy(P, xp(k, d), y) - fy(P, xp(k, -d), y)) / (2 * d); adjoint.push(fx - N.dot(psi, Rx)); direct.push(fx + N.dot(dfdy, N.solve(Ry, Rx.map((v) => -v)))); }
    const diff = Math.max(...keys.map((_, j) => Math.abs(fd[j] - adjoint[j]) / Math.max(Math.abs(adjoint[j]), 1e-12 * Math.abs(f0)))), elasX = keys.map((k, j) => (adjoint[j] * (k === 'sweep' || k === 'h' ? scale(k) : P.x0[k])) / f0);
    const elasP = PARAMS.map(([, ap]) => { const a = { ...P }, b = { ...P }; ap(a, 0.01); ap(b, -0.01); return (F(a, P.x0) - F(b, P.x0)) / (0.02 * f0); });
    const steps = N.logspace(1e-8, 1e-1, 15), errFd = steps.map((s) => { const d = s * P.x0.AR; return Math.abs((F(P, xp('AR', d), 1e-9) - F(P, xp('AR', -d), 1e-9)) / (2 * d) - adjoint[0]) / Math.abs(adjoint[0]); }), errFw = steps.map((s) => { const d = s * P.x0.AR; return Math.abs((F(P, xp('AR', d), 1e-9) - f0) / d - adjoint[0]) / Math.abs(adjoint[0]); });
    const all = [...labels.map((n, j) => [n, elasX[j]]), ...PARAMS.map(([n], j) => [n, elasP[j]])].filter((r) => Number.isFinite(r[1])).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1]));
    if (diff > 1e-3) warnings.push('Finite-difference and adjoint derivatives differ by more than 0.1%: the finite-difference step is too large (truncation) or too small (iteration noise). See the step-size plot.');
    return {
      kpis: [
        { key: 'adjoint_fd_max_diff_pct', label: 'Largest adjoint versus finite-difference gap', value: 100 * diff, unit: '%', status: diff < 1e-3 ? 'ok' : 'warn' },
        { key: 'direct_adjoint_diff', label: 'Direct versus adjoint gap', value: Math.max(...keys.map((_, j) => Math.abs(direct[j] - adjoint[j]) / Math.max(Math.abs(adjoint[j]), 1e-12 * Math.abs(f0)))), unit: '-' },
        { key: 'best_fd_step', label: 'Best central-difference step (relative)', value: steps[N.argmin(errFd)], unit: '-' },
        ...keys.map((k, j) => ({ key: `elasticity_${k}`, label: `Elasticity to ${labels[j].toLowerCase()}`, value: elasX[j], unit: '%/%', note: k === 'sweep' || k === 'h' ? 'Per 1% of the reference scale' : '' })),
        { key: 'top_driver_elasticity', label: `Strongest driver: ${all[0][0]}`, value: all[0][1], unit: '%/%' },
      ],
      plots: [
        { type: 'bar', title: `Tornado: elasticity of ${i.objective} (per cent change per 1% change)`, ylabel: 'Elasticity [%/%]', categories: all.map((r) => r[0]), series: [{ name: 'Elasticity', y: all.map((r) => r[1]) }] },
        { type: 'line', title: 'Finite-difference error versus step size (derivative with respect to aspect ratio)', xlabel: 'Relative step [-]', ylabel: 'Relative error against the adjoint [-]', xlog: true, ylog: true, series: [{ name: 'Central difference', x: steps, y: errFd.map((v) => Math.max(v, 1e-16)), style: 'line+points' }, { name: 'Forward difference', x: steps, y: errFw.map((v) => Math.max(v, 1e-16)), style: 'line+points' }], annotations: [{ x: i.fd_step, label: 'Step used' }] },
      ],
      tables: [{ title: 'Total derivatives of the objective', columns: ['Variable', 'Finite difference', 'Coupled direct', 'Coupled adjoint', 'Unit'], rows: keys.map((k, j) => [labels[j], fd[j], direct[j], adjoint[j], `${obj.unit} per ${['-', 'm²', 'deg', '-', 'm'][j]}`]) },
        { title: 'Adjoint variables (sensitivity of the objective to each coupling residual)', columns: ['Residual', 'ψ'], rows: ['Take-off mass', 'Wing mass', 'Energy mass'].map((n, j) => [n, psi[j]]) }],
      warnings, models: [...MODELS, 'Central and forward finite differences on the converged analysis', 'Coupled direct method: (∂R/∂y)·dy/dx = −∂R/∂x', 'Coupled adjoint method: (∂R/∂y)ᵀψ = (∂f/∂y)ᵀ, df/dx = ∂f/∂x − ψᵀ∂R/∂x'],
      assumptions: [...ASSUME, 'Partial derivatives of the explicit discipline equations are taken by central differences at the converged state', 'Technology-parameter elasticities freeze the baseline calibration'],
    };
  },
  verify() {
    const i = Object.fromEntries([...BOUNDS, ...MODEL, ...NUM].map((f) => [f.key, f.default])); i.fd_step = 1e-4; const o = N.kv(sens.run(i));
    const d = (Math.sin(1 + 1e-5) - Math.sin(1 - 1e-5)) / 2e-5;
    return [
      N.check('Coupled adjoint equals finite difference of the converged analysis', o.adjoint_fd_max_diff_pct / 100 + 1, 1, 2e-5, 'Implicit function theorem'),
      N.check('Direct and adjoint methods agree', o.direct_adjoint_diff + 1, 1, 1e-8, 'Transpose identity'),
      N.check('Central difference of sin at 1', d, Math.cos(1), 1e-9, 'Second-order truncation h²/6'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], rows = res.plots[0];
    if (!rows) return [];
    const top = rows.categories.slice(0, 3).map((n, j) => `${n} (${rows.series[0].y[j].toFixed(2)})`);
    out.push({ severity: 'advise', title: `Biggest levers on ${i.objective}: ${rows.categories[0]}`, detail: `Elasticities: ${top.join(', ')}. A value of 1 means a 1% change moves the objective by 1%.`, action: `Direct technology effort at ${rows.categories[0].toLowerCase()} first. Drag and fuel-consumption gains translate one-for-one into CO₂.`, basis: 'Normalised total derivatives' });
    out.push({ severity: 'info', title: 'Use adjoint-style derivatives when variables are many', detail: `Adjoint and finite differences agree to ${o.adjoint_fd_max_diff_pct.toExponential(1)}%; the adjoint needs one linear solve regardless of the number of variables.`, action: `For finite differences keep the relative step near ${o.best_fd_step.toExponential(0)}; smaller steps amplify solver noise, larger ones add truncation error.`, basis: 'Step-size study' });
    return out;
  },
};

const LIMIT = ['Take-off field length', 'Approach speed', 'Cruise thrust', 'Fuel volume'];
const robust = {
  id: 'robust', title: 'Robust and reliability-based design', fidelity: 'numerical',
  summary: 'Re-optimises aspect ratio and wing area when drag, fuel consumption, material strength and payload are uncertain, minimising the mean plus a multiple of the standard deviation, and computes the reliability index of a chosen constraint by the first-order reliability method.',
  equations: ['Reliability-based multidisciplinary optimisation equations', 'Constrained nonlinear optimisation equations', 'Multidisciplinary coupling equations'],
  applicable: fixedWing,
  inputs: [...BOUNDS.filter((f) => ['objective', 'AR_min', 'AR_max', 'S_band'].includes(f.key)), ...MODEL, ...NUM,
    { key: 'cov_cd0', label: 'Uncertainty of zero-lift drag (CoV)', unit: '-', default: 0.05, min: 0, max: 0.3, group: 'Uncertainty' },
    { key: 'cov_sfc', label: 'Uncertainty of fuel consumption (CoV)', unit: '-', default: 0.03, min: 0, max: 0.3, group: 'Uncertainty' },
    { key: 'cov_sigma', label: 'Uncertainty of allowable stress (CoV)', unit: '-', default: 0.05, min: 0, max: 0.3, group: 'Uncertainty' },
    { key: 'cov_payload', label: 'Uncertainty of payload (CoV)', unit: '-', default: 0.05, min: 0, max: 0.3, group: 'Uncertainty' },
    { key: 'cov_clmax', label: 'Uncertainty of CLmax (CoV)', unit: '-', default: 0.04, min: 0, max: 0.3, group: 'Uncertainty', help: 'Used by the reliability analysis' },
    { key: 'k_sigma', label: 'Robustness weight k (mean + k·σ)', unit: '-', default: 2, min: 0, max: 6, group: 'Uncertainty', help: 'Also the number of standard deviations of margin demanded of each constraint' },
    { key: 'limit_state', label: 'Constraint for the reliability index', type: 'select', options: LIMIT, default: 'Take-off field length', group: 'Uncertainty' },
    { key: 'nSamples', label: 'Samples per design', unit: '', default: 8, min: 4, max: 200, step: 1, discrete: true, group: 'Numerics' },
    { key: 'seed', label: 'Random seed', unit: '', default: 11, min: 1, max: 1e9, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: modelDefaults,
  run(i, ctx) {
    const P = makeModel(i), warnings = [...P.warn], vars = varsOf(i, 'AR, S'), obj = OBJ[i.objective], base = evalDesign(P, P.x0);
    if (!base.ok || vars.length < 2) return { kpis: [{ key: 'robust_objective', label: 'Robust objective', value: 0, unit: '', status: 'bad' }], tables: [{ title: 'Status', columns: ['Message'], rows: [['The baseline mass loop diverges or the variable bounds are empty.']] }], warnings: [...warnings, 'The robust study could not run.'], models: MODELS, assumptions: ASSUME };
    const covs = [i.cov_cd0, i.cov_sfc, i.cov_sigma, i.cov_payload, i.cov_clmax], apply = (z) => { const Q = { ...P }; PARAMS[0][1](Q, covs[0] * z[0]); PARAMS[1][1](Q, covs[1] * z[1]); PARAMS[2][1](Q, covs[2] * z[2]); PARAMS[3][1](Q, covs[3] * z[3]); if (z.length > 4) { Q.CLmax_to = P.CLmax_to * (1 + covs[4] * z[4]); Q.CLmax_land = P.CLmax_land * (1 + covs[4] * z[4]); } return Q; };
    const u = N.rng(i.seed), Zs = N.lhs(Math.round(i.nSamples), 5, u).map((r) => r.map((p) => N.clamp(N.normInv(p), -3, 3))), Ps = Zs.map(apply), f0 = obj.f(base, P), tol = 1e-8; // common random numbers across designs
    const stats = (x) => { const fs = [], gs = base.g.map(() => []); for (const Q of Ps) { const o = evalDesign(Q, x, tol); if (!o.ok) return null; fs.push(obj.f(o, Q) / f0); o.g.forEach((g, k) => gs[k].push(g)); } return { mean: N.mean(fs), sd: N.std(fs), gm: gs.map(N.mean), gs: gs.map(N.std) }; };
    const det = designProblem(i, P, vars, i.objective), rd = alOptimise(det.fn, toU(P, vars), { outer: 4, inner: 45, polish: 4 }); ctx?.progress?.(0.3, 'Robust search');
    const rfn = (uu) => { const s = stats(toX(P, vars, uu)); return s ? { f: s.mean + i.k_sigma * s.sd, g: s.gm.map((m, k) => m - i.k_sigma * s.gs[k]), ok: true } : { f: 1e3, g: base.g.map(() => -1), ok: false }; };
    const rr = alOptimise(rfn, rd.u, { outer: 3, inner: 28, polish: 2 }), xd = toX(P, vars, rd.u), xr = toX(P, vars, rr.u), sB = stats(P.x0), sD = stats(xd), sR = stats(xr), feasR = Math.min(...rr.g) > -3e-3;
    // FORM on the chosen constraint at the deterministic optimum
    const ci = [0, 1, 4, 5][LIMIT.indexOf(i.limit_state)], gf = (x) => (z) => { const o = evalDesign(apply(z), x, 1e-10); return o.ok ? o.g[ci] : -1; }, fD = form(gf(xd), 5, 12), fR = form(gf(xr), 5, 12), names = ['Zero-lift drag', 'Fuel consumption', 'Allowable stress', 'Payload', 'CLmax'];
    if (!feasR) warnings.push(`No design keeps every constraint ${i.k_sigma} standard deviations inside its limit; the robust result is the best compromise. Reduce k, the uncertainties or relax a limit.`);
    if (P.electric && ci === 5) warnings.push('Fuel volume does not apply to a battery aircraft; its reliability index is not meaningful.');
    warnings.push(`Uncertainty magnitudes are assumptions to be replaced with evidence; ${Math.round(i.nSamples)} samples give only a rough standard deviation.`);
    const pct = (v) => 100 * (v - 1);
    return {
      kpis: [
        { key: 'robust_objective', label: 'Robust objective mean + k·σ (relative to nominal baseline)', value: sR ? sR.mean + i.k_sigma * sR.sd : NaN, unit: '-', status: feasR ? 'ok' : 'warn' },
        { key: 'robust_AR', label: 'Robust aspect ratio', value: xr.AR, unit: '-' }, { key: 'robust_S_m2', label: 'Robust wing area', value: xr.S, unit: 'm²' },
        { key: 'det_AR', label: 'Deterministic-optimum aspect ratio', value: xd.AR, unit: '-' }, { key: 'det_S_m2', label: 'Deterministic-optimum wing area', value: xd.S, unit: 'm²' },
        { key: 'robust_mean_pct', label: 'Robust design: mean objective versus baseline', value: sR ? pct(sR.mean) : NaN, unit: '%' }, { key: 'robust_sd_pct', label: 'Robust design: standard deviation', value: sR ? 100 * sR.sd : NaN, unit: '%' },
        { key: 'price_of_robustness_pct', label: 'Price of robustness (mean objective)', value: sR && sD ? 100 * (sR.mean - sD.mean) : NaN, unit: '%-points', note: 'Robust minus deterministic optimum' },
        { key: 'beta_det', label: `Reliability index β of "${i.limit_state}" at the deterministic optimum`, value: fD.beta, unit: '-', status: fD.beta >= 2 ? 'ok' : fD.beta >= 1 ? 'warn' : 'bad' },
        { key: 'pf_det', label: 'Probability of violating it', value: fD.pf, unit: '-' },
        { key: 'beta_robust', label: 'Reliability index β at the robust design', value: fR.beta, unit: '-', status: fR.beta >= 2 ? 'ok' : fR.beta >= 1 ? 'warn' : 'bad' },
        { key: 'pf_robust', label: 'Probability of violating it (robust)', value: fR.pf, unit: '-' },
      ],
      plots: [
        { type: 'bar', title: 'Mean and spread of the objective', ylabel: 'Objective relative to nominal baseline [-]', categories: ['Baseline', 'Deterministic optimum', 'Robust optimum'], series: [{ name: 'Mean', y: [sB, sD, sR].map((s) => (s ? s.mean : 0)) }, { name: `Mean + ${i.k_sigma}σ`, y: [sB, sD, sR].map((s) => (s ? s.mean + i.k_sigma * s.sd : 0)) }] },
        { type: 'bar', title: `FORM importance factors α² for "${i.limit_state}"`, ylabel: 'Share of variance at the design point [-]', categories: names, series: [{ name: 'Deterministic optimum', y: fD.alpha.map((a) => a * a) }, { name: 'Robust optimum', y: fR.alpha.map((a) => a * a) }] },
      ],
      tables: [{ title: 'Constraint margins under uncertainty (mean − k·σ, per cent)', columns: ['Constraint', 'Baseline', 'Deterministic optimum', 'Robust optimum'], rows: CON_NAMES.map((n, k) => [n, ...[sB, sD, sR].map((s) => (s ? 100 * (s.gm[k] - i.k_sigma * s.gs[k]) : NaN))]) },
        { title: 'Most probable failure point (standard-normal coordinates)', columns: ['Variable', 'u* deterministic', 'u* robust'], rows: names.map((n, k) => [n, fD.u[k], fR.u[k]]) }],
      warnings, models: [...MODELS, 'Latin-hypercube sampling with common random numbers', 'Mean + k·σ robust objective with moment-based constraint margins', 'First-order reliability method (Hasofer–Lind / Rackwitz–Fiessler)'], assumptions: [...ASSUME, 'Independent normal uncertainties truncated at ±3σ for sampling', 'Only aspect ratio and wing area are re-optimised; other variables stay at the baseline'],
    };
  },
  verify() {
    // linear limit state g = R − S with R ~ N(10, 1.5), S ~ N(6, 2): β = 4 / 2.5 exactly
    const f = form((z) => 10 + 1.5 * z[0] - (6 + 2 * z[1]), 2), q = form((z) => 3 - z[0] - 0.1 * z[0] * z[0] * 0 - z[1], 2);
    return [
      N.check('FORM reliability index of a linear normal limit state', f.beta, 4 / Math.sqrt(1.5 ** 2 + 2 ** 2), 1e-9, 'Cornell index (exact for linear Gaussian)'),
      N.check('Failure probability Φ(−β)', f.pf, N.normCdf(-1.6), 1e-9, 'Definition'),
      N.check('Design point lies on the limit state: β = 3/√2', q.beta, 3 / Math.SQRT2, 1e-9, 'Distance from the origin to a plane'),
      N.check('Importance factors sum to one', N.sum(f.alpha.map((a) => a * a)), 1, 1e-12, 'Unit direction vector'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!o.robust_AR) return [];
    if (o.beta_det < 1.5) out.push({ severity: 'warn', title: 'The deterministic optimum is fragile', detail: `"${i.limit_state}" has reliability index β = ${o.beta_det.toFixed(2)}: about ${(100 * o.pf_det).toFixed(0)}% chance of violation once uncertainty is counted, because the optimiser parks the design exactly on the limit.`, action: `Adopt the robust design (aspect ratio ${o.robust_AR.toFixed(2)}, area ${o.robust_S_m2.toFixed(1)} m²) with β = ${o.beta_robust.toFixed(2)}, at a mean penalty of ${o.price_of_robustness_pct.toFixed(2)} percentage points.`, basis: 'First-order reliability method' });
    else out.push({ severity: 'info', title: 'The chosen constraint has reasonable reliability', detail: `β = ${o.beta_det.toFixed(2)} at the deterministic optimum.`, action: 'Check the other limit states with the selector before concluding.', basis: 'First-order reliability method' });
    const a = res.plots[1], top = N.argmax(a.series[0].y);
    out.push({ severity: 'advise', title: `Reduce uncertainty in ${a.categories[top].toLowerCase()} first`, detail: `It carries ${(100 * a.series[0].y[top]).toFixed(0)}% of the variance of "${i.limit_state}".`, action: 'Targeted testing or higher-fidelity analysis on this one quantity shrinks the margin you must carry, which saves mass and fuel directly.', basis: 'FORM importance factors α²' });
    return out;
  },
};

/** Rotor sizing model: hover and cruise power from momentum theory, component masses, gross-mass loop. */
function rotorDesign(i, x, mFixed) {
  const [R, sig, vt] = x, A = Math.PI * R * R, nR = Math.max(1, Math.round(i.n_rotors)), a = isa(i.alt_m, i.dISA), rho = a.rho, elec = i.ptype === 'electric';
  const pw = (m) => {
    const T = (m * G0) / nR, Ph = (nR * ((i.kappa * T ** 1.5) / Math.sqrt(2 * rho * A) + (rho * A * vt ** 3 * sig * i.cd0) / 8)) / i.eta_mech, mu = i.V / vt;
    const vh = Math.sqrt(T / (2 * rho * A)); let vi = vh; for (let k = 0; k < 12; k++) vi = (vh * vh) / Math.sqrt(i.V * i.V + vi * vi); // Glauert forward-flight inflow (level disk)
    const Pc = (nR * (i.kappa * T * vi + ((rho * A * vt ** 3 * sig * i.cd0) / 8) * (1 + 4.65 * mu * mu)) + 0.5 * rho * i.V ** 3 * i.f_m2) / i.eta_mech;
    return { T, Ph, Pc, mu };
  };
  const parts = (m) => {
    const p = pw(m), blade = i.blade_kg_m2 * sig * A * nR, hub = i.hub_frac * blade * ((vt * vt) / R) / ((i.vt_ref * i.vt_ref) / i.R_ref), Q = (p.Ph * i.eta_mech) / nR / (vt / R), drive = nR * i.drive_k * Q ** i.drive_n, pp = (i.pp_kg_kW * i.power_margin * Math.max(p.Ph, p.Pc)) / 1e3;
    const E = p.Ph * i.hover_min * 60 + (p.Pc * i.range_km * 1e3) / Math.max(i.V, 1), en = ((elec ? E / (i.eta_elec * i.batt_Wh_kg * 3600) : i.bsfc * E) * (1 + i.reserve_frac));
    return { ...p, blade, hub, drive, pp, energy: en, E, sum: blade + hub + drive + pp + en };
  };
  if (mFixed === undefined) return parts(i.mtow_ref);
  const s = solveMDA((y) => [mFixed + i.payload + parts(y[0]).sum], [i.mtow_ref], { tol: 1e-9, maxIter: 80, method: 'jacobi' });
  if (!s.converged || !(s.y[0] > 0) || s.y[0] > 20 * i.mtow_ref) return { ok: false };
  const p = parts(s.y[0]), ctSig = p.T / (rho * A * vt * vt * sig);
  return { ok: true, m: s.y[0], ...p, ctSig, DL: p.T / A, g: [1 - ctSig / i.ct_sigma_max, 1 - (vt + i.V) / a.a / i.mtip_max, 1 - p.mu / i.mu_max, 1 - R / i.R_max], iters: s.iters };
}
const RCON = ['Blade loading CT/σ', 'Advancing-tip Mach number', 'Advance ratio', 'Rotor radius limit'];
const rotor = {
  id: 'rotor', title: 'Rotor sizing optimisation for rotorcraft and multirotors', fidelity: 'reduced-order',
  summary: 'Chooses rotor radius, solidity and tip speed to minimise gross mass, hover power or mission energy, with blade-loading, tip-Mach, advance-ratio and size limits, closing the loop between power, drive-train, energy and gross mass.',
  equations: ['Constrained nonlinear optimisation equations', 'Karush–Kuhn–Tucker optimality conditions', 'Multidisciplinary coupling equations'],
  applicable: (c) => (c.rotor.R_m > 0 && c.meta.type !== 'aeroplane' ? true : 'This analysis needs a lifting rotor.'),
  inputs: [
    { key: 'objective', label: 'Objective', type: 'select', options: ['gross mass', 'hover power', 'mission energy'], default: 'gross mass', group: 'Optimisation' },
    { key: 'has_wing', label: 'Vehicle also has a wing (lift-plus-cruise)', type: 'bool', default: false, group: 'Baseline', help: 'When set, results are published as rotor_opt_… so that the wing optimisation remains the aircraft-level result' },
    { key: 'mtow_ref', label: 'Reference gross mass', unit: 'kg', default: 9980, min: 0.1, group: 'Baseline' }, { key: 'payload', label: 'Payload', unit: 'kg', default: 3600, min: 0, group: 'Baseline' },
    { key: 'n_rotors', label: 'Lifting rotors', unit: '', default: 1, min: 1, max: 16, step: 1, discrete: true, group: 'Baseline' },
    { key: 'R_ref', label: 'Baseline rotor radius', unit: 'm', default: 8.18, min: 0.02, group: 'Baseline' }, { key: 'sig_ref', label: 'Baseline solidity', unit: '-', default: 0.0825, min: 0.01, max: 0.4, group: 'Baseline' }, { key: 'vt_ref', label: 'Baseline tip speed', unit: 'm/s', default: 221, min: 20, max: 280, group: 'Baseline' },
    { key: 'V', label: 'Cruise speed', unit: 'm/s', default: 72, min: 0, max: 120, group: 'Mission' }, { key: 'range_km', label: 'Range', unit: 'km', default: 500, min: 0, group: 'Mission' }, { key: 'hover_min', label: 'Hover time', unit: 'min', default: 10, min: 0, group: 'Mission' },
    { key: 'reserve_frac', label: 'Energy reserve fraction', unit: '-', default: 0.15, min: 0, max: 1, group: 'Mission' }, { key: 'alt_m', label: 'Design altitude', unit: 'm', default: 500, min: -500, max: 6000, group: 'Mission' }, { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -40, max: 40, group: 'Mission' },
    { key: 'cd0', label: 'Blade profile drag coefficient', unit: '-', default: 0.0095, min: 0.004, max: 0.05, group: 'Aerodynamics' }, { key: 'kappa', label: 'Induced power factor', unit: '-', default: 1.15, min: 1, max: 1.6, group: 'Aerodynamics' }, { key: 'f_m2', label: 'Equivalent flat-plate area', unit: 'm²', default: 3.35, min: 0, group: 'Aerodynamics' },
    { key: 'eta_mech', label: 'Transmission efficiency incl. anti-torque', unit: '-', default: 0.85, min: 0.5, max: 1, group: 'Propulsion' }, { key: 'ptype', label: 'Powerplant', type: 'select', options: ['turboshaft', 'piston', 'electric'], default: 'turboshaft', group: 'Propulsion' },
    { key: 'bsfc', label: 'BSFC', unit: 'kg/W/s', default: 7.8e-8, min: 0, group: 'Propulsion' }, { key: 'batt_Wh_kg', label: 'Battery pack specific energy', unit: 'Wh/kg', default: 190, min: 50, max: 800, group: 'Propulsion' }, { key: 'eta_elec', label: 'Battery-to-shaft efficiency', unit: '-', default: 0.88, min: 0.5, max: 0.99, group: 'Propulsion' },
    { key: 'power_margin', label: 'Installed power / required power', unit: '-', default: 1.3, min: 1, max: 3, group: 'Propulsion' },
    { key: 'blade_kg_m2', label: 'Blade mass per unit blade area', unit: 'kg/m²', default: 27, min: 0.2, max: 80, group: 'Mass coefficients (illustrative)' }, { key: 'hub_frac', label: 'Hub and controls mass / blade mass at baseline', unit: '-', default: 0.6, min: 0, max: 3, group: 'Mass coefficients (illustrative)', help: 'Scaled with centrifugal acceleration Vtip²/R' },
    { key: 'drive_k', label: 'Drive mass coefficient', unit: 'kg/(N·m)ⁿ', default: 0.12, min: 0, max: 5, group: 'Mass coefficients (illustrative)', help: 'Gearbox or direct-drive motor mass = k·(rotor torque)ⁿ' }, { key: 'drive_n', label: 'Drive mass torque exponent n', unit: '-', default: 0.8, min: 0.5, max: 1.2, group: 'Mass coefficients (illustrative)' },
    { key: 'pp_kg_kW', label: 'Powerplant mass per installed kW', unit: 'kg/kW', default: 0.22, min: 0.02, max: 3, group: 'Mass coefficients (illustrative)' },
    { key: 'ct_sigma_max', label: 'Blade loading limit CT/σ', unit: '-', default: 0.12, min: 0.04, max: 0.2, group: 'Constraints', help: 'Hover value leaving stall margin for manoeuvre; 0.10–0.14 typical' }, { key: 'mtip_max', label: 'Advancing-tip Mach limit', unit: '-', default: 0.88, min: 0.3, max: 0.98, group: 'Constraints' },
    { key: 'mu_max', label: 'Advance-ratio limit', unit: '-', default: 0.4, min: 0.05, max: 0.8, group: 'Constraints' }, { key: 'R_max', label: 'Largest allowed rotor radius', unit: 'm', default: 9.5, min: 0.03, group: 'Constraints' },
  ],
  defaults: (c, up, d) => {
    const heli = c.meta.type === 'helicopter', elec = c.prop.type === 'electric', nR = heli ? 1 : c.prop.n_eng, R = c.rotor.R_m, bladeA = c.rotor.n_blades * c.rotor.chord_m * R, vt = d.v_tip || 200, wing = c.wing.S_m2 > 0;
    return { has_wing: wing, mtow_ref: c.mass.mtow_kg, payload: c.mass.payload_kg, n_rotors: nR, R_ref: R, sig_ref: d.solidity || undefined, vt_ref: vt, V: wing ? 0 : c.mission.cruise_V_ms || c.flight.V_ms, range_km: wing ? 0 : c.mission.range_km, hover_min: Math.max(c.mission.hover_min, wing ? 1 : 0), alt_m: c.mission.cruise_alt_m || c.atm.alt_m, dISA: c.atm.dISA_K,
      cd0: c.rotor.cd0, f_m2: c.rotor.flat_plate_m2, eta_mech: heli ? 0.85 : 0.92, ptype: ['turboshaft', 'piston', 'electric'].includes(c.prop.type) ? c.prop.type : 'turboshaft', bsfc: c.prop.bsfc_kg_Ws, batt_Wh_kg: 0.75 * (BATTERIES[c.systems.batt_chem]?.wh_kg || 250),
      blade_kg_m2: c.rotor.blade_mass_kg > 0 && bladeA > 0 ? (c.rotor.blade_mass_kg * c.rotor.n_blades) / bladeA : undefined, drive_k: elec ? 0.15 : 0.12, drive_n: elec ? 1 : 0.8, pp_kg_kW: elec ? 0.25 : c.prop.type === 'piston' ? 0.9 : 0.22, R_max: 1.2 * R, reserve_frac: elec ? 0.2 : 0.15,
      ct_sigma_max: Math.max(0.12, d.solidity > 0 ? 1.05 * ((c.mass.mtow_kg * G0) / nR) / (isa(c.mission.cruise_alt_m || 0).rho * Math.PI * R * R * vt * vt * d.solidity) : 0) };
  },
  run(i) {
    const warnings = [], x0 = [i.R_ref, i.sig_ref, i.vt_ref], b0 = rotorDesign(i, x0), mFixed = i.mtow_ref - i.payload - b0.sum, lo = [0.6 * i.R_ref, 0.5 * i.sig_ref, 0.7 * i.vt_ref], hi = [Math.max(i.R_max, i.R_ref), 1.8 * i.sig_ref, Math.min(1.25 * i.vt_ref, 260)];
    const fail = (m) => ({ kpis: [{ key: 'rotor_opt_ran', label: 'Rotor optimisation ran', value: 0, unit: '', status: 'bad' }], tables: [{ title: 'Status', columns: ['Message'], rows: [[m]] }], warnings: [m], models: ['Momentum-theory rotor sizing'], assumptions: [ILLUS] });
    if (mFixed < 0.1 * i.mtow_ref) return fail('The illustrative component masses already exceed the reference gross mass less payload; lower the mass coefficients or check the reference masses.');
    const base = rotorDesign(i, x0, mFixed); if (!base.ok) return fail('The baseline gross-mass loop diverges for this mission.');
    const fo = { 'gross mass': (o) => o.m, 'hover power': (o) => o.Ph, 'mission energy': (o) => o.E }[i.objective], f0 = fo(base), X = (u) => u.map((v, j) => lo[j] + v * (hi[j] - lo[j]));
    const fn = (u) => { const o = rotorDesign(i, X(u), mFixed); return o.ok ? { f: fo(o) / f0, g: o.g, ok: true, o } : { f: 1e3, g: RCON.map(() => -1), ok: false }; };
    const r = alOptimise(fn, x0.map((v, j) => N.clamp((v - lo[j]) / (hi[j] - lo[j]), 0, 1)), { outer: 4, inner: 60, polish: 6 }), best = fn(r.u), o = best.o || base, x = X(r.u), k = kkt(fn, r.u), feas = Math.min(...best.g) > -2e-3, elec = i.ptype === 'electric', pre = i.has_wing ? 'rotor_' : '';
    if (!feas) warnings.push('No rotor inside the bounds meets every limit; the least-infeasible design is shown.');
    base.g.forEach((g, c) => { if (g < 0) warnings.push(`The baseline rotor violates "${RCON[c]}" by ${(-100 * g).toFixed(1)}%.`); });
    warnings.push(`${ILLUS}. Aeroelastic, acoustic and autorotation requirements are not constraints here; check the result in Suites 6, 10 and 11.`);
    const DLs = N.linspace(0.5 * base.DL, 2 * base.DL, 30), swp = DLs.map((dl) => { const Rr = Math.sqrt((i.mtow_ref * G0) / (Math.max(1, Math.round(i.n_rotors)) * Math.PI * dl)), q = rotorDesign(i, [Rr, x[1], x[2]], mFixed); return q.ok ? q : null; }), row = (n, a, b, u) => [n, a, b, u, a ? (100 * (b - a)) / Math.abs(a) : 0];
    return {
      kpis: [
        { key: `${pre}opt_mtow_kg`, label: 'Optimised gross mass', value: o.m, unit: 'kg' },
        { key: `${pre}opt_fuel_kg`, label: 'Optimised mission fuel incl. reserve', value: elec ? 0 : o.energy, unit: 'kg', note: elec ? 'Electric: see battery mass' : '' },
        { key: `${pre}opt_improvement_pct`, label: `Improvement in ${i.objective}`, value: 100 * (1 - fo(o) / f0), unit: '%', status: feas ? 'ok' : 'warn' },
        { key: 'opt_R_m', label: 'Optimised rotor radius', value: x[0], unit: 'm' }, { key: 'opt_solidity', label: 'Optimised solidity', value: x[1], unit: '-' }, { key: 'opt_tip_speed_ms', label: 'Optimised tip speed', value: x[2], unit: 'm/s' },
        { key: 'opt_disk_loading_Pa', label: 'Optimised disk loading', value: o.DL, unit: 'N/m²' }, { key: 'opt_hover_power_W', label: 'Optimised hover power', value: o.Ph, unit: 'W' }, { key: 'opt_cruise_power_W', label: 'Optimised cruise power', value: o.Pc, unit: 'W' },
        { key: 'opt_ct_sigma', label: 'Blade loading CT/σ', value: o.ctSig, unit: '-', status: o.ctSig <= i.ct_sigma_max * 1.002 ? 'ok' : 'bad' }, { key: 'opt_battery_kg', label: 'Optimised battery mass', value: elec ? o.energy : 0, unit: 'kg' },
        { key: 'kkt_residual', label: 'KKT stationarity residual', value: k.residual, unit: '-', status: k.residual < 0.05 ? 'ok' : 'warn' }, { key: 'opt_feasible', label: 'Optimum is feasible', value: feas ? 1 : 0, unit: '', status: feas ? 'ok' : 'bad' },
      ],
      plots: [
        { type: 'line', title: 'Gross mass and hover power versus disk loading (other variables at optimum)', xlabel: 'Disk loading [N/m²]', ylabel: 'Relative to baseline [-]', series: [{ name: 'Gross mass', x: DLs, y: swp.map((q) => (q ? q.m / base.m : NaN)) }, { name: 'Hover power', x: DLs, y: swp.map((q) => (q ? q.Ph / base.Ph : NaN)) }], annotations: [{ x: o.DL, label: 'Optimum' }, { x: base.DL, label: 'Baseline' }] },
        { type: 'bar', title: 'Mission-dependent mass breakdown', ylabel: 'Mass [kg]', categories: ['Baseline', 'Optimised'], stacked: true, series: [{ name: 'Blades', y: [base.blade, o.blade] }, { name: 'Hub and controls', y: [base.hub, o.hub] }, { name: 'Drive', y: [base.drive, o.drive] }, { name: 'Powerplant', y: [base.pp, o.pp] }, { name: elec ? 'Battery' : 'Fuel', y: [base.energy, o.energy] }] },
        { type: 'bar', title: 'Constraint margins', ylabel: 'Margin [%]', categories: RCON, series: [{ name: 'Baseline', y: base.g.map((g) => 100 * g) }, { name: 'Optimised', y: best.g.map((g) => 100 * g) }] },
      ],
      tables: [{ title: 'Baseline versus optimised rotor', columns: ['Quantity', 'Baseline', 'Optimised', 'Unit', 'Change [%]'], rows: [row('Rotor radius', x0[0], x[0], 'm'), row('Solidity', x0[1], x[1], '-'), row('Tip speed', x0[2], x[2], 'm/s'), row('Gross mass', base.m, o.m, 'kg'), row('Hover power', base.Ph, o.Ph, 'W'), row('Cruise power', base.Pc, o.Pc, 'W'), row('Disk loading', base.DL, o.DL, 'N/m²'), row('CT/σ', base.ctSig, o.ctSig, '-'), row(elec ? 'Battery mass' : 'Fuel mass', base.energy, o.energy, 'kg')] },
        { title: 'Constraints and Lagrange multipliers', columns: ['Constraint', 'Margin [%]', 'Multiplier λ'], rows: RCON.map((n, c) => [n, 100 * best.g[c], k.lambda[c]]) }],
      warnings, models: ['Momentum theory with induced-power factor and uniform profile power', 'Glauert forward-flight inflow', 'Torque-based drive mass, power-based powerplant mass, area-based blade mass', 'Augmented Lagrangian with Nelder–Mead and projected-gradient polish'],
      assumptions: ['Fixed (non-rotor) mass is the reference gross mass less payload and the baseline component estimates', 'Hub mass scales with blade mass and centrifugal acceleration', 'No blade stall, compressibility drag rise or rotor–rotor interference', ILLUS],
    };
  },
  verify() {
    const i = Object.fromEntries(rotor.inputs.map((f) => [f.key, f.default])), q = rotorDesign({ ...i, cd0: 0, kappa: 1, eta_mech: 1, alt_m: 0, V: 0 }, [8, 0.08, 200]), T = i.mtow_ref * G0;
    const r = rotor.run(i), o = N.kv(r), sum = r.plots[1].series.reduce((s, p) => s + p.y[1], 0), fx = i.mtow_ref - i.payload - r.plots[1].series.reduce((s, p) => s + p.y[0], 0);
    return [
      N.check('Ideal hover power T^1.5/√(2ρA)', q.Ph, T ** 1.5 / Math.sqrt(2 * isa(0).rho * Math.PI * 64), 1e-9, 'Rankine–Froude momentum theory'),
      N.check('Gross-mass loop closes at the optimum', fx + i.payload + sum, o.opt_mtow_kg, 1e-7, 'Sum of components'),
      N.check('Optimised rotor respects the blade-loading limit', Math.min(o.opt_ct_sigma, i.ct_sigma_max) / o.opt_ct_sigma, 1, 3e-3, 'Constraint enforcement'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.rotor_opt_ran === 0) return [{ severity: 'warn', title: 'Rotor optimisation could not run', detail: res.warnings[0], action: 'Check the reference masses and the illustrative mass coefficients.', basis: 'Mass closure' }];
    const imp = o.opt_improvement_pct ?? o.rotor_opt_improvement_pct;
    out.push({ severity: imp > 1 ? 'advise' : 'info', title: `${imp.toFixed(1)}% lower ${i.objective} from re-sizing the rotor`, detail: `Radius ${o.opt_R_m.toPrecision(3)} m, solidity ${o.opt_solidity.toFixed(3)}, tip speed ${o.opt_tip_speed_ms.toFixed(0)} m/s; disk loading ${o.opt_disk_loading_Pa.toFixed(0)} N/m².`, action: 'A larger, slower rotor cuts induced power, fuel or battery energy and noise, at the price of blade, hub and gearbox mass; confirm with Suite 6 (trim, autorotation) and Suite 11 (noise).', basis: 'Constrained optimum of the momentum-theory sizing model' });
    const act = res.tables[1].rows.filter((r) => r[2] > 1e-4).sort((a, b) => b[2] - a[2]);
    if (act.length) out.push({ severity: 'advise', title: `Binding limit: ${act[0][0]}`, detail: act.map((r) => `${r[0]} (λ = ${Number(r[2]).toFixed(2)})`).join('; ') + '.', action: `A 1% relaxation of "${act[0][0]}" is worth about ${Number(act[0][2]).toFixed(2)}% of the objective; see whether aerofoil or planform technology can supply it.`, basis: 'Lagrange multipliers' });
    return out;
  },
};

export default {
  id: 'mdao', n: 23,
  tagline: 'Which combination of wing, rotor and cruise condition gives the best aircraft once aerodynamics, structure, propulsion and mission are traded together.',
  analyses: [mda, opt, pareto, surrogate, sens, robust, rotor],
  consumes: [
    { from: 'cfd', keys: ['CD0', 'e_oswald'], why: 'Calibrates the drag build-up' },
    { from: 'fea', keys: ['wing_struct_mass_kg'], why: 'Calibrates the wing non-optimum factor' },
    { from: 'propulsion', keys: ['tsfc_kg_Ns'], why: 'Cruise fuel consumption' },
    { from: 'performance', keys: ['tofl_m'], why: 'Calibrates the take-off field-length factor' },
  ],
  provides: [
    { key: 'opt_mtow_kg', label: 'Optimised take-off mass', unit: 'kg' }, { key: 'opt_AR', label: 'Optimised aspect ratio', unit: '-' }, { key: 'opt_S_m2', label: 'Optimised wing area', unit: 'm²' },
    { key: 'opt_fuel_kg', label: 'Optimised fuel', unit: 'kg' }, { key: 'opt_improvement_pct', label: 'Objective improvement', unit: '%' },
  ],
  handoff: [
    { model: 'High-fidelity MDAO with CFD and FEA in the loop and discrete adjoints', why: 'Needs 3-D flow and structural solvers, mesh deformation and adjoint-capable codes; here the disciplines are low-order and partials are finite-differenced', tool: 'Adjoint CFD/FEA MDO frameworks' },
    { model: 'Topology optimisation (SIMP / level-set)', why: 'Not implemented in this suite; needs a dedicated finite-element density solver and filtering', tool: 'Structural topology-optimisation software, then Suite 2 for verification' },
    { model: 'Distributed MDO architectures (IDF, collaborative optimisation, analytical target cascading)', why: 'Only the multidisciplinary-feasible (MDF) formulation is implemented, which suits a model this cheap', tool: 'MDO architecture frameworks' },
    { model: 'Particle swarm and differential evolution', why: 'NSGA-II and augmented-Lagrangian Nelder–Mead cover the global and local roles', tool: 'General optimisation libraries' },
    { model: 'Full SQP with quasi-Newton Hessian', why: 'A projected-gradient polish follows the derivative-free search instead', tool: 'SQP / interior-point NLP solvers' },
  ],
};
