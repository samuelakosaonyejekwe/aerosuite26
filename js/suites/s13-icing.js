// Suite 13 — Aircraft Icing and Ice Protection.
// Lagrangian droplet trajectories and collection efficiency on the leading edge, Messinger mass and energy balance
// with runback and time-stepped ice growth, empirical aerodynamic penalties, anti-icing and de-icing power,
// an icing-condition envelope map, and rotor / propeller blade icing with centrifugal shedding.
// The leading edge is represented by its equivalent cylinder; 3-D icing CFD is listed in `handoff`.

import * as N from '../core/numerics.js';
import { isa, G0, R_AIR, sutherland, kAir, CP_AIR, PR_AIR, pSat } from '../core/atmosphere.js';

const PI = Math.PI, TF = 273.15, RHO_W = 1000, C_W = 4218, C_I = 2050, L_F = 3.34e5, L_V = 2.5e6, L_S = 2.834e6, REC = 0.85, SC = 0.6;
const pSatIce = (T) => 611.2 * Math.exp((22.46 * (T - TF)) / (T - 0.53)); // Magnus form over ice
const LANGMUIR_D = { w: [0.05, 0.1, 0.2, 0.3, 0.2, 0.1, 0.05], r: [0.31, 0.52, 0.71, 1, 1.37, 1.74, 2.22] };
const rLE = (tc, c) => 1.1019 * tc * tc * c; // NACA 4-digit leading-edge radius
/** Air state at the icing condition. */
const air = (alt, T_C) => { const p = isa(alt).p, T = T_C + TF, rho = p / (R_AIR * T), mu = sutherland(T); return { p, T, rho, mu, k: kAir(T) }; };

// ---- droplet inertia and trajectories ------------------------------------------------------------
/** Langmuir–Blodgett inertia parameters for droplet diameter d [m] on a cylinder of radius R [m]. */
export function inertia(V, d, R, a) {
  const K = (RHO_W * d * d * V) / (18 * a.mu * R), Red = (a.rho * V * d) / a.mu, lam = 1 / (0.8388 + 0.001483 * Red + 0.1847 * Math.sqrt(Red)), K0 = 0.125 + lam * (K - 0.125);
  const x = K0 > 0.125 ? 1.4 * (K0 - 0.125) ** 0.84 : 0;
  return { K, Red, lam, K0, beta0: x / (1 + x) };
}
/**
 * Droplet trajectories in potential flow about a unit cylinder (lengths in radii, speeds in V∞).
 * dv/dt = (CD·Re/24)·(u − v)/K with the Schiller–Naumann drag law (Stokes drag when Red = 0).
 * Returns local collection efficiency β(θ) = dy0/ds, total efficiency E and the impingement limit θmax [rad].
 */
export function trajectories(K, Red, nTraj = 24) {
  const shoot = (y0) => {
    const x0 = -50, r2 = x0 * x0 + y0 * y0; let s = [x0, y0, 1 - (x0 * x0 - y0 * y0) / (r2 * r2), (-2 * x0 * y0) / (r2 * r2)];
    const f = (q) => { const r = q[0] * q[0] + q[1] * q[1], r4 = r * r, du = 1 - (q[0] * q[0] - q[1] * q[1]) / r4 - q[2], dv = (-2 * q[0] * q[1]) / r4 - q[3], c = (1 + 0.15 * (Red * Math.hypot(du, dv)) ** 0.687) / K; return [q[2], q[3], c * du, c * dv]; };
    for (let n = 0; n < 40000; n++) {
      const r = Math.hypot(s[0], s[1]), dt = Math.min(0.5 * K, Math.max(0.15 * (r - 1), 0.002)); // explicit stability of the drag term and resolution near the surface
      const k1 = f(s), k2 = f(s.map((v, j) => v + 0.5 * dt * k1[j])), k3 = f(s.map((v, j) => v + 0.5 * dt * k2[j])), k4 = f(s.map((v, j) => v + dt * k3[j])), sn = s.map((v, j) => v + (dt / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])), rn = Math.hypot(sn[0], sn[1]);
      if (rn <= 1) { const fr = (r - 1) / (r - rn), x = s[0] + fr * (sn[0] - s[0]), y = s[1] + fr * (sn[1] - s[1]); return Math.atan2(y, -x); }
      if (sn[0] > 0.05 || (r - 1 < 1e-5 && Math.hypot(sn[2], sn[3]) < 1e-4)) return NaN;
      s = sn;
    }
    return NaN;
  };
  if (Number.isNaN(shoot(1e-4))) return { theta: [0], beta: [0], E: 0, thetaMax: 0, beta0: 0 };
  let lo = 1e-4, hi = 1.05; for (let k = 0; k < 22; k++) { const m = 0.5 * (lo + hi); if (Number.isNaN(shoot(m))) hi = m; else lo = m; }
  const y0 = N.linspace(0, lo, nTraj + 1), th = y0.map((y, k) => (k ? shoot(y) : 0)), theta = [], beta = [];
  for (let k = 1; k <= nTraj; k++) { if (!(th[k] > th[k - 1])) continue; theta.push(0.5 * (th[k] + th[k - 1])); beta.push((y0[k] - y0[k - 1]) / (th[k] - th[k - 1])); }
  const b0 = beta.length > 1 ? beta[0] + ((beta[0] - beta[1]) * theta[0]) / (theta[1] - theta[0]) : beta[0] || 0;
  return { theta: [0, ...theta, th[nTraj]], beta: [Math.min(b0, 1), ...beta, 0], E: lo, thetaMax: th[nTraj], beta0: Math.min(b0, 1) };
}
/** Collection efficiency for a droplet spectrum (monodisperse at the MVD or the seven-bin Langmuir D distribution). */
function collection(V, mvd, R, a, dist, nTraj) {
  const bins = dist === 'Langmuir D' ? LANGMUIR_D.r.map((r, k) => ({ d: mvd * r, w: LANGMUIR_D.w[k] })) : [{ d: mvd, w: 1 }];
  const parts = bins.map((b) => { const q = inertia(V, b.d, R, a); return { ...b, q, t: trajectories(q.K, q.Red, nTraj) }; }), thMax = Math.max(...parts.map((p) => p.t.thetaMax));
  const fn = (th) => N.sum(parts.map((p) => (th <= p.t.thetaMax ? p.w * N.interp1(p.t.theta, p.t.beta, th) : 0)));
  return { fn, E: N.sum(parts.map((p) => p.w * p.t.E)), thetaMax: thMax, beta0: fn(0), mid: inertia(V, mvd, R, a), parts };
}
/** Convective coefficient on the front of a cylinder (laminar Frössling distribution), with a roughness multiplier. */
const hCyl = (th, V, D, a, kh = 1) => (a.k / D) * 1.14 * Math.sqrt((a.rho * V * D) / a.mu) * PR_AIR ** 0.4 * Math.max(1 - (th / (PI / 2)) ** 3, 0.25) * kh;

// ---- Messinger mass and energy balance at one surface station -----------------------------------
/**
 * o = {beta, mIn (runback entering, kg/m²/s), h, V, T (static, K), p, lwc (kg/m³), rh, q (heater flux W/m²)}.
 * Returns freezing fraction n, surface temperature, ice, evaporation and runback mass fluxes and the heat terms at the solution.
 */
export function messinger(o) {
  const { h, V, T, p } = o, mImp = o.beta * o.lwc * V, mTot = mImp + o.mIn, Trec = T + (REC * V * V) / (2 * CP_AIR), hG = (h / CP_AIR) * (PR_AIR / SC) ** (2 / 3), q = o.q || 0;
  const evap = (Ts, ice) => Math.max(0, (hG * 0.622 * ((ice ? pSatIce(Ts) : pSat(Ts)) - o.rh * pSat(T))) / p);
  const out = (n, Ts, mE) => { const mIce = mTot > 0 ? Math.max(0, Math.min(n * mTot, mTot - mE)) : 0; return { n, Ts, Trec, mImp, mTot, mIce, mEvap: mE, mOut: Math.max(0, mTot - mIce - mE), qConv: h * (Ts - Trec), qEvap: mE * (Ts < TF ? L_S : L_V), qSens: mImp * C_W * (Math.min(Ts, TF) - T), qKin: 0.5 * mImp * V * V, qFreeze: mIce * L_F }; };
  if (!(mTot > 0)) return out(0, Trec + q / h, 0); // dry surface: adiabatic-wall recovery temperature plus heater rise
  const base = (Ts) => h * (Ts - Trec) + mImp * C_W * (Math.min(Ts, TF) - T) - 0.5 * mImp * V * V - q;
  const mE0 = Math.min(evap(TF, false), mTot), n0 = (base(TF) + mE0 * L_V) / (mTot * L_F);
  if (n0 >= 0 && n0 <= 1) return out(n0, TF, Math.min(mE0, (1 - n0) * mTot + 1e-30));
  if (n0 > 1) { // rime: everything freezes on impact and the ice cools below 0 °C
    const g = (Ts) => base(Ts) + Math.min(evap(Ts, true), mTot) * L_S - mTot * L_F - mTot * C_I * (TF - Ts), Ts = N.findRoot(g, Math.min(T, Trec) - 80, TF, 80, 1e-9);
    const t = Number.isFinite(Ts) ? Ts : TF; return out(1, t, Math.min(evap(t, true), mTot));
  }
  // no freezing: the surface runs wet above 0 °C (or dries out if the heater evaporates all of the water)
  const g = (Ts) => base(Ts) + (mImp * C_W + o.mIn * C_W) * Math.max(Ts - TF, 0) + Math.min(evap(Ts, false), mTot) * L_V, Ts = N.findRoot(g, TF, TF + 400, 200, 1e-9), t = Number.isFinite(Ts) ? Ts : TF;
  return out(0, t, Math.min(evap(t, false), mTot));
}

// ---- shared inputs ---------------------------------------------------------------------------------
const CLOUD = [
  { key: 'V', label: 'True airspeed', unit: 'm/s', default: 90, min: 5, max: 300, group: 'Flight condition' },
  { key: 'alt_m', label: 'Pressure altitude', unit: 'm', default: 3000, min: 0, max: 12000, group: 'Flight condition' },
  { key: 'T_C', label: 'Static air temperature', unit: '°C', default: -10, min: -40, max: 5, group: 'Cloud', help: 'Supercooled cloud exists mostly between 0 and −20 °C, rarely below −30 °C' },
  { key: 'lwc', label: 'Liquid water content', unit: 'g/m³', default: 0.5, min: 0, max: 5, group: 'Cloud', help: 'Typically 0.1–0.8 in layer cloud and up to about 3 in convective cloud; take the design value from the certification envelope you are working to' },
  { key: 'mvd', label: 'Median volumetric diameter', unit: 'µm', default: 20, min: 5, max: 1000, group: 'Cloud', help: '15–40 µm in ordinary cloud; above about 50 µm is supercooled large droplet (SLD) icing' },
  { key: 'chord', label: 'Chord', unit: 'm', default: 1.5, min: 0.01, group: 'Geometry' },
  { key: 'tc', label: 'Thickness ratio', unit: '-', default: 0.12, min: 0.04, max: 0.3, group: 'Geometry', help: 'Sets the leading-edge radius 1.1019·(t/c)²·c' },
  { key: 'd_le', label: 'Leading-edge diameter override', unit: 'm', default: 0, min: 0, group: 'Geometry', help: '0 uses the NACA leading-edge radius; enter a diameter for a strut, probe or cylinder' },
];
const cloudDefaults = (c, d) => {
  const wing = c.wing.S_m2 > 0, alt = Math.min(c.atm.alt_m, 5000), V0 = wing ? c.flight.V_ms : 0.75 * d.v_tip;
  return { V: (wing ? V0 * Math.sqrt(isa(c.atm.alt_m).rho / isa(alt).rho) : V0) || undefined, alt_m: alt, T_C: N.clamp(isa(alt, c.atm.dISA_K).T - TF, -25, -5), chord: (wing ? d.mac : c.rotor.chord_m) || undefined, tc: wing ? c.wing.tc : 0.12 };
};
const radius = (i) => (i.d_le > 0 ? i.d_le / 2 : rLE(i.tc, i.chord));
const sldNote = (i, w) => { if (i.mvd > 50) w.push(`MVD ${i.mvd} µm is in the supercooled-large-droplet range: droplet splashing, break-up and impingement well aft of the protected zone are not modelled, so catch on the leading edge is over-predicted and aft ice is missed.`); };

const impingement = {
  id: 'impingement', title: 'Droplet trajectories and collection efficiency', fidelity: 'numerical',
  summary: 'Tracks supercooled droplets through the flow around the leading edge to find how much of the cloud water strikes the surface, where, and how far back the ice can form.',
  equations: ['Lagrangian particle motion equations', 'Droplet impingement equations', 'Lagrangian droplet tracking models', 'Supercooled large droplet models'],
  inputs: [...CLOUD,
    { key: 'dist', label: 'Droplet size distribution', type: 'select', options: ['Monodisperse (MVD)', 'Langmuir D'], default: 'Langmuir D', group: 'Cloud' },
    { key: 'nTraj', label: 'Trajectories per droplet size', unit: '', default: 24, min: 6, max: 200, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => cloudDefaults(c, d),
  run(i) {
    const a = air(i.alt_m, i.T_C), R = radius(i), col = collection(i.V, i.mvd * 1e-6, R, a, i.dist, Math.round(N.clamp(i.nTraj, 6, 200))), q = col.mid, warnings = [];
    const th = N.linspace(0, Math.max(col.thetaMax, 0.02), 50), catchRate = col.E * (i.lwc / 1e3) * i.V * 2 * R;
    sldNote(i, warnings);
    if (col.E === 0) warnings.push('No droplets reach the surface: the inertia parameter is below the critical value and the droplets follow the air around the leading edge.');
    if (i.V / Math.sqrt(1.4 * R_AIR * a.T) > 0.5) warnings.push('Mach number above 0.5: the incompressible potential flow used for the trajectories under-estimates the local flow acceleration.');
    const ks = N.logspace(0.13, 100, 40).map((k) => 1.4 * (k - 0.125) ** 0.84), mono = col.parts.find((p) => p.d === i.mvd * 1e-6) || col.parts[0];
    const paths = [0.25, 0.5, 0.75, 1].map((f) => { const y0 = f * mono.t.E; return { name: `y0 = ${(y0 * R * 1e3).toFixed(2)} mm`, y0 }; });
    return {
      kpis: [
        { key: 'K_inertia', label: 'Inertia parameter K', value: q.K, unit: '-', note: 'Droplets impinge only when K > 1/8 (Stokes drag)' },
        { key: 'K0_modified', label: 'Modified inertia parameter K0', value: q.K0, unit: '-' },
        { key: 'Re_droplet', label: 'Droplet Reynolds number', value: q.Red, unit: '-' },
        { key: 'beta0', label: 'Stagnation collection efficiency', value: col.beta0, unit: '-' },
        { key: 'beta0_correlation', label: 'Stagnation efficiency, Langmuir–Blodgett fit', value: q.beta0, unit: '-', note: 'Empirical curve fit at the MVD, for comparison' },
        { key: 'E_total', label: 'Total collection efficiency', value: col.E, unit: '-' },
        { key: 'theta_limit_deg', label: 'Impingement limit angle', value: N.deg(col.thetaMax), unit: 'deg' },
        { key: 's_limit_mm', label: 'Impingement limit, surface distance', value: col.thetaMax * R * 1e3, unit: 'mm', note: 'Measured from the stagnation line on each surface' },
        { key: 'catch_kg_m_s', label: 'Water catch per unit span', value: catchRate, unit: 'kg/(m·s)' },
        { key: 'r_le_mm', label: 'Leading-edge radius', value: R * 1e3, unit: 'mm' },
        { key: 'sld_flag', label: 'Supercooled large droplets', value: i.mvd > 50 ? 1 : 0, unit: '', status: i.mvd > 50 ? 'warn' : 'ok', note: '1 when MVD exceeds 50 µm' },
      ],
      plots: [
        { type: 'line', title: 'Local collection efficiency', xlabel: 'Surface angle from stagnation [deg]', ylabel: 'β [-]', series: [{ name: i.dist, x: th.map(N.deg), y: th.map(col.fn) }, ...(col.parts.length > 1 ? [{ name: 'MVD only', x: th.map(N.deg), y: th.map((t) => (t <= mono.t.thetaMax ? N.interp1(mono.t.theta, mono.t.beta, t) : 0)), style: 'dash' }] : [])], annotations: [{ x: N.deg(col.thetaMax), label: 'Impingement limit' }] },
        { type: 'bar', title: 'Collection efficiency by droplet size', ylabel: 'Total collection efficiency [-]', categories: col.parts.map((p) => `${(p.d * 1e6).toFixed(0)} µm`), series: [{ name: 'E', y: col.parts.map((p) => p.t.E) }] },
        { type: 'line', title: 'Stagnation collection efficiency against inertia', xlabel: 'Modified inertia parameter K0 [-]', ylabel: 'β0 [-]', xlog: true, series: [{ name: 'Langmuir–Blodgett fit', x: N.logspace(0.13, 100, 40), y: ks.map((x) => x / (1 + x)) }, { name: 'Trajectories, each droplet size', x: col.parts.map((p) => Math.max(p.q.K0, 0.13)), y: col.parts.map((p) => p.t.beta0), style: 'points' }] },
      ],
      tables: [{ title: 'Droplet bins', columns: ['Diameter [µm]', 'LWC share', 'K', 'K0', 'β0', 'E', 'Limit [deg]'], rows: col.parts.map((p) => [+(p.d * 1e6).toFixed(1), p.w, +p.q.K.toFixed(3), +p.q.K0.toFixed(3), +p.t.beta0.toFixed(3), +p.t.E.toFixed(3), +N.deg(p.t.thetaMax).toFixed(1)]) }],
      outputs: { catch_start_y: paths[3].y0 },
      warnings,
      models: ['Lagrangian droplet tracking (RK4) in potential flow about the leading-edge cylinder', 'Schiller–Naumann sphere drag', i.dist === 'Langmuir D' ? 'Langmuir D droplet spectrum (7 bins)' : 'Monodisperse droplets at the MVD', 'Langmuir–Blodgett modified inertia parameter (comparison)'],
      assumptions: ['The leading edge is replaced by a cylinder of the leading-edge radius at zero incidence', 'Droplets are rigid spheres released at the local air velocity fifty radii upstream; gravity is neglected', 'No splashing, bouncing or break-up'],
    };
  },
  convergence: { param: 'nTraj', label: 'Trajectories per droplet size', levels: [8, 16, 32, 64], metric: 'beta0' },
  verify() {
    const heavy = trajectories(1e5, 0, 40), sub = trajectories(0.11, 0, 12), sup = trajectories(0.5, 0, 24), m = sup.theta.length >> 1;
    return [
      N.check('Ballistic limit: total efficiency → 1', heavy.E, 1, 2e-3, 'Straight-line droplets, K → ∞'),
      N.check('Ballistic limit: β(θ) = cos θ', N.interp1(heavy.theta, heavy.beta, PI / 3), 0.5, 5e-3, 'Geometric projection'),
      N.check('No impingement below K = 1/8 (Stokes drag)', sub.E, 0, 1e-9, 'Langmuir & Blodgett critical inertia parameter'),
      N.check('Total efficiency equals ∫β dθ', N.trapz(sup.theta, sup.beta), sup.E, 0.01, 'Mass conservation between far-field and surface'),
      N.check('β decreases monotonically from stagnation', sup.beta[m] < sup.beta[0] && sup.beta[m] > 0 ? 1 : 0, 1, 1e-12, 'Physical ordering'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (i.mvd > 50) out.push({ severity: 'warn', title: 'Supercooled large droplet conditions', detail: `MVD ${i.mvd} µm. Ice forms aft of conventional protected zones as ridges that cause large lift and control losses.`, action: 'Treat as outside ordinary cloud-icing protection: check the aircraft is approved for these conditions, plan an immediate exit, and analyse aft impingement with a splashing-capable 3-D code.', basis: 'Large-droplet impingement physics' });
    out.push({ severity: 'info', title: 'Size the protected zone from the impingement limit', detail: `Droplets strike up to ${o.s_limit_mm.toFixed(1)} mm (${o.theta_limit_deg.toFixed(0)}°) from the stagnation line at zero incidence.`, action: 'Extend heaters or boots beyond this limit with allowance for the incidence range and for runback; a smaller heated area cuts electrical or bleed demand and therefore fuel burn.', basis: 'Trajectory impingement limits' });
    if (o.E_total > 0.6) out.push({ severity: 'advise', title: 'Small leading edge collects efficiently', detail: `Total collection efficiency ${(100 * o.E_total).toFixed(0)}%: thin sections, tails, probes and rotor blades ice faster than the wing.`, action: 'Check tailplane and blade icing first; they limit the aircraft before the wing does.', basis: 'Inertia parameter scaling' });
    return out;
  },
};

// ---- time-stepped Messinger accretion ---------------------------------------------------------------
function accrete(i) {
  const a = air(i.alt_m, i.T_C), R0 = radius(i), nS = Math.round(N.clamp(i.nSteps, 1, 40)), nSt = Math.round(N.clamp(i.nStations, 8, 120)), th = N.linspace(0, PI / 2, nSt), dth = th[1] - th[0], dt = (i.time_min * 60) / nS;
  const tk = new Array(nSt).fill(0), lwc = i.lwc / 1e3; let st = [], col = null, bal = { imp: 0, ice: 0, evap: 0, run: 0 };
  for (let s = 0; s < nS; s++) {
    const R = R0 + tk[0]; col = collection(i.V, i.mvd * 1e-6, R, a, i.dist, Math.round(N.clamp(i.nTraj, 6, 100)));
    let mIn = 0; st = [];
    for (let j = 0; j < nSt; j++) {
      const m = messinger({ beta: col.fn(th[j]), mIn, h: hCyl(th[j], i.V, 2 * R, a, i.k_h), V: i.V, T: a.T, p: a.p, lwc, rh: i.rh, q: i.q_heat });
      st.push(m); mIn = j === 0 ? m.mOut / 2 : m.mOut; tk[j] += (m.mIce / i.rho_ice) * dt;
      const w = (j === 0 ? 1 : 2) * R * dth * dt; bal.imp += m.mImp * w; bal.ice += m.mIce * w; bal.evap += m.mEvap * w; if (j === nSt - 1) bal.run += m.mOut * 2 * R * dth * dt;
    }
  }
  const mass = N.sum(tk.map((t, j) => (j === 0 ? 1 : 2) * i.rho_ice * t * (R0 + t / 2) * dth));
  return { a, R0, th, tk, st, col, bal, mass };
}
const accretion = {
  id: 'accretion', title: 'Ice accretion: Messinger mass and energy balance', fidelity: 'numerical',
  summary: 'Balances impinging water, runback, freezing, evaporation, convection and heating at stations around the leading edge, and grows the ice in time steps to give thickness, shape, mass and whether the ice is rime or glaze.',
  equations: ['Messinger icing energy balance equations', 'Mass conservation equations for water films', 'Energy conservation equations for ice and liquid water', 'Stefan phase-change equations', 'Convective heat-transfer equations', 'Evaporation equations', 'Ice growth equations', 'Runback water models', 'Rime and glaze ice models'],
  inputs: [...CLOUD,
    { key: 'rh', label: 'Relative humidity of the airstream', unit: '-', default: 1, min: 0, max: 1, group: 'Cloud', help: '1 inside cloud' },
    { key: 'time_min', label: 'Exposure time', unit: 'min', default: 10, min: 0.1, max: 120, group: 'Flight condition', help: 'Time in cloud without ice protection or since the last de-icing cycle' },
    { key: 'rho_ice', label: 'Ice density', unit: 'kg/m³', default: 900, min: 300, max: 917, group: 'Model', help: '917 clear glaze, 850–900 typical, down to 300–600 for feathery rime' },
    { key: 'k_h', label: 'Roughness heat-transfer multiplier', unit: '-', default: 1.5, min: 1, max: 4, group: 'Model', help: 'Ratio of rough-ice to smooth laminar-cylinder heat transfer: 1 clean surface, 1.5–3 on accreted ice' },
    { key: 'q_heat', label: 'Surface heater flux', unit: 'W/m²', default: 0, min: 0, max: 1e5, group: 'Model', help: '0 for an unprotected surface' },
    { key: 'dist', label: 'Droplet size distribution', type: 'select', options: ['Monodisperse (MVD)', 'Langmuir D'], default: 'Monodisperse (MVD)', group: 'Cloud' },
    { key: 'nSteps', label: 'Time steps (geometry updates)', unit: '', default: 4, min: 1, max: 40, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nStations', label: 'Surface stations per side', unit: '', default: 31, min: 8, max: 120, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nTraj', label: 'Trajectories per droplet size', unit: '', default: 16, min: 6, max: 100, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => cloudDefaults(c, d),
  run(i) {
    const r = accrete(i), s0 = r.st[0], jm = N.argmax(r.tk), deg = r.th.map(N.deg), warnings = [], n0 = s0.n, kind = s0.mIce <= 0 ? 'no ice' : n0 >= 0.999 ? 'rime' : n0 > 0.6 ? 'mixed' : 'glaze';
    sldNote(i, warnings);
    if (s0.mIce <= 0 && s0.mImp > 0) warnings.push('No ice forms at the stagnation line: kinetic heating, the heater or a temperature near 0 °C keeps the surface wet. Check the runback for freezing further aft.');
    if (r.tk[0] > 0.5 * r.R0) warnings.push('The ice is thick compared with the leading-edge radius: the cylinder model no longer represents the flow, and horn growth alters impingement and heat transfer. Use more time steps and treat the shape as indicative.');
    if (jm > 0 && kind === 'glaze') warnings.push(`Maximum thickness is ${deg[jm].toFixed(0)}° off the stagnation line: a horn is forming from runback water, the most damaging shape for lift and drag.`);
    const res = r.bal.imp > 0 ? (r.bal.imp - r.bal.ice - r.bal.evap - r.bal.run) / r.bal.imp : 0, rr = (t) => (r.R0 + t) / r.R0;
    const sx = [...r.th.map((t, j) => -rr(r.tk[j]) * Math.cos(t)).reverse(), ...r.th.map((t, j) => -rr(r.tk[j]) * Math.cos(t))], sy = [...r.th.map((t, j) => -rr(r.tk[j]) * Math.sin(t)).reverse(), ...r.th.map((t, j) => rr(r.tk[j]) * Math.sin(t))];
    const cx = N.linspace(-PI / 2, PI / 2, 60);
    return {
      kpis: [
        { key: 'ice_thickness_mm', label: 'Ice thickness at the stagnation line', value: r.tk[0] * 1e3, unit: 'mm' },
        { key: 'ice_max_mm', label: 'Maximum ice thickness', value: r.tk[jm] * 1e3, unit: 'mm', note: `at ${deg[jm].toFixed(0)}° from stagnation` },
        { key: 'ice_mass_kg_m', label: 'Ice mass per unit span', value: r.mass, unit: 'kg/m' },
        { key: 'freezing_fraction', label: 'Freezing fraction at stagnation', value: n0, unit: '-', note: `${kind}; 1 = rime, below about 0.6 = glaze with runback` },
        { key: 'T_surface_C', label: 'Surface temperature at stagnation', value: s0.Ts - TF, unit: '°C' },
        { key: 'T_recovery_C', label: 'Adiabatic-wall recovery temperature', value: s0.Trec - TF, unit: '°C', status: s0.Trec < TF ? 'warn' : 'ok', note: 'Above 0 °C kinetic heating alone prevents icing' },
        { key: 'accretion_rate_mm_min', label: 'Stagnation growth rate', value: (s0.mIce / i.rho_ice) * 6e4, unit: 'mm/min' },
        { key: 'beta0', label: 'Stagnation collection efficiency', value: r.col.beta0, unit: '-' },
        { key: 'runback_kg_m_s', label: 'Runback leaving the leading edge', value: 2 * r.st[r.st.length - 1].mOut * r.R0 * (r.th[1] - r.th[0]), unit: 'kg/(m·s)', note: 'Can refreeze aft of a heated zone' },
        { key: 'mass_balance_error', label: 'Water mass balance residual', value: res, unit: '-', status: Math.abs(res) < 1e-9 ? 'ok' : 'warn' },
        { key: 'h_stag_Wm2K', label: 'Stagnation heat-transfer coefficient', value: hCyl(0, i.V, 2 * (r.R0 + r.tk[0]), r.a, i.k_h), unit: 'W/(m²·K)' },
      ],
      plots: [
        { type: 'line', title: 'Ice thickness around the leading edge', xlabel: 'Surface angle from stagnation [deg]', ylabel: 'Ice thickness [mm]', series: [{ name: `After ${i.time_min} min`, x: deg, y: r.tk.map((v) => v * 1e3) }] },
        { type: 'line', title: 'Freezing fraction and collection efficiency (final step)', xlabel: 'Surface angle from stagnation [deg]', ylabel: 'Fraction [-]', series: [{ name: 'Freezing fraction', x: deg, y: r.st.map((m) => (m.mTot > 0 ? m.n : 0)) }, { name: 'Collection efficiency β', x: deg, y: r.th.map(r.col.fn) }] },
        { type: 'line', title: 'Ice shape on the leading-edge cylinder', xlabel: 'x / R [-]', ylabel: 'y / R [-]', equalAspect: true, series: [{ name: 'Ice surface', x: sx, y: sy }, { name: 'Clean leading edge', x: cx.map((t) => -Math.cos(t)), y: cx.map(Math.sin), style: 'dash' }] },
        { type: 'bar', title: 'Stagnation heat balance', ylabel: 'Heat flux [kW/m²]', categories: ['Freezing (release)', 'Kinetic (release)', 'Heater', 'Convection', 'Evaporation', 'Droplet warming'], series: [{ name: 'Heat flux', y: [s0.qFreeze, s0.qKin, i.q_heat, s0.qConv, s0.qEvap, s0.qSens].map((v) => v / 1e3) }] },
        { type: 'line', title: 'Surface temperature', xlabel: 'Surface angle from stagnation [deg]', ylabel: 'Temperature [°C]', series: [{ name: 'Surface', x: deg, y: r.st.map((m) => m.Ts - TF) }], annotations: [{ y: 0, label: 'Freezing' }] },
      ],
      outputs: { ice_type_code: kind === 'rime' ? 1 : kind === 'mixed' ? 2 : kind === 'glaze' ? 3 : 0 },
      warnings,
      models: ['Messinger control-volume mass and energy balance with runback', 'Lagrangian droplet impingement on the growing leading-edge cylinder', 'Frössling laminar cylinder heat transfer with a roughness multiplier (empirical)', 'Chilton–Colburn heat–mass transfer analogy for evaporation and sublimation', 'Multi-step quasi-steady ice growth'],
      assumptions: ['Leading edge treated as a cylinder at zero incidence; ice grows normal to the original surface', 'Unfrozen water runs back along the surface at 0 °C; no shedding, splashing or film dynamics', 'Constant ice density; conduction into the structure neglected', 'Recovery factor 0.85 and freestream speed in the kinetic terms'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps (geometry updates)', levels: [1, 2, 4, 8], metric: 'ice_thickness_mm' },
  calibration: { params: [{ key: 'k_h', min: 1, max: 4 }, { key: 'rho_ice', min: 400, max: 917 }], sweep: 'time_min', target: 'ice_thickness_mm', note: 'Supply icing-tunnel thickness against exposure time (or temperature); the rough-surface heat-transfer multiplier and ice density are the uncertain parameters.' },
  verify() {
    const base = { V: 80, alt_m: 0, T_C: -25, lwc: 0.2, mvd: 20, chord: 1, tc: 0.12, d_le: 0.05, rh: 1, time_min: 5, rho_ice: 900, k_h: 1, q_heat: 0, dist: 'Monodisperse (MVD)', nSteps: 1, nStations: 16, nTraj: 16 };
    const rime = accretion.run(base), o = N.kv(rime), warm = N.kv(accretion.run({ ...base, T_C: 4 })), dry = messinger({ beta: 0, mIn: 0, h: 200, V: 100, T: 263.15, p: 9e4, lwc: 0, rh: 1, q: 0 }), hot = messinger({ beta: 0, mIn: 0, h: 200, V: 100, T: 263.15, p: 9e4, lwc: 0, rh: 1, q: 4000 });
    return [
      N.check('Rime limit: thickness = β0·LWC·V·t/ρ_ice', o.ice_thickness_mm, (o.beta0 * 0.2e-3 * 80 * 300) / 900 * 1e3, 0.03, 'All impinging water freezes; sublimation removes about 2%'),
      N.check('Rime limit: freezing fraction = 1', o.freezing_fraction, 1, 1e-12, 'Messinger model, cold limit'),
      N.check('Non-freezing condition: no ice above 0 °C', warm.ice_thickness_mm, 0, 1e-12, 'Limiting behaviour'),
      N.check('Dry surface: wall at recovery temperature', dry.Ts, 263.15 + (0.85 * 1e4) / (2 * CP_AIR), 1e-12, 'Adiabatic wall'),
      N.check('Dry heated surface: ΔT = q/h', hot.Ts - dry.Ts, 20, 1e-10, 'Newton cooling'),
      N.check('Water mass is conserved', 1 + o.mass_balance_error, 1, 1e-10, 'Impinged = frozen + evaporated + runback'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.ice_thickness_mm > 0 && o.freezing_fraction < 0.6) out.push({ severity: 'warn', title: 'Glaze ice with runback', detail: `Freezing fraction ${o.freezing_fraction.toFixed(2)}: water runs back and freezes as horns (${o.ice_max_mm.toFixed(1)} mm maximum).`, action: 'Expect the largest lift and drag penalties here. Change altitude to colder air (rime is more benign) or out of cloud; make sure ice protection covers the runback zone.', basis: 'Messinger freezing fraction' });
    if (o.accretion_rate_mm_min > 1) out.push({ severity: 'warn', title: 'High accretion rate', detail: `${o.accretion_rate_mm_min.toFixed(2)} mm/min at the leading edge, ${o.ice_mass_kg_m.toFixed(2)} kg per metre of span after ${i.time_min} min.`, action: 'Shorten de-icing cycles, and limit exposure: every minute in these conditions adds drag and therefore fuel burn.', basis: 'Accretion rate (severity classes are operational, not fixed by this model)' });
    if (o.T_recovery_C > 0) out.push({ severity: 'info', title: 'Kinetic heating protects this surface', detail: `Recovery temperature is ${o.T_recovery_C.toFixed(1)} °C at ${i.V} m/s.`, action: 'Slower surfaces (inboard rotor blade, tailplane on approach) are not protected: analyse them at their own speed.', basis: 'Adiabatic-wall recovery temperature' });
    return out;
  },
};

// ---- empirical aerodynamic penalties ----------------------------------------------------------------
/** Empirical section penalties from an ice (or roughness) height ratio k/c and freezing fraction. */
const penalty = (kc, n, A, m, ch, cr) => ({ dcl: kc > 0 ? Math.min(0.6, A * kc ** m) : 0, dcd: kc > 0 ? cr + ch * kc * (1 - 0.7 * N.clamp(n, 0, 1)) : 0 });
const degradation = {
  id: 'degradation', title: 'Aerodynamic penalties of the ice (empirical)', fidelity: 'reduced-order',
  summary: 'Converts the ice thickness into estimated losses of maximum lift and increases in drag, stall speed, thrust and fuel, using simple empirical trends that must be calibrated against icing-tunnel or flight data.',
  equations: ['Iced-airfoil aerodynamic models', 'Aerodynamic–icing feedback equations'],
  applicable: (c) => (c.wing.S_m2 > 0 ? true : 'This analysis needs a wing; rotor and propeller penalties are in the blade-icing analysis.'),
  inputs: [
    { key: 'ice_mm', label: 'Ice thickness on the leading edge', unit: 'mm', default: 8, min: 0, max: 150, group: 'Ice', help: 'From the accretion analysis' },
    { key: 'freezing_fraction', label: 'Freezing fraction', unit: '-', default: 0.5, min: 0, max: 1, group: 'Ice', help: '1 = streamlined rime, low values = glaze horns' },
    { key: 'ice_mass_kg_m', label: 'Ice mass per unit span', unit: 'kg/m', default: 0.5, min: 0, group: 'Ice' },
    { key: 'span_frac', label: 'Unprotected fraction of wing span', unit: '-', default: 1, min: 0, max: 1, group: 'Ice', help: '1 = no ice protection or system failed' },
    { key: 'chord', label: 'Mean chord', unit: 'm', default: 1.5, min: 0.01, group: 'Aircraft' },
    { key: 'S', label: 'Wing area', unit: 'm²', default: 16.2, min: 0.01, group: 'Aircraft' },
    { key: 'b', label: 'Span', unit: 'm', default: 11, min: 0.05, group: 'Aircraft' },
    { key: 'tail_frac', label: 'Tail area / wing area', unit: '-', default: 0.3, min: 0, max: 1, group: 'Aircraft' },
    { key: 'mass_kg', label: 'Aircraft mass', unit: 'kg', default: 1111, min: 0.1, group: 'Aircraft' },
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 62, min: 1, group: 'Flight condition' },
    { key: 'alt_m', label: 'Altitude', unit: 'm', default: 2400, min: 0, max: 12000, group: 'Flight condition' },
    { key: 'CLmax', label: 'Clean maximum lift coefficient', unit: '-', default: 1.5, min: 0.3, group: 'Aerodynamics' },
    { key: 'CD0', label: 'Clean zero-lift drag coefficient', unit: '-', default: 0.027, min: 0.003, group: 'Aerodynamics' },
    { key: 'k', label: 'Induced drag factor', unit: '-', default: 0.05, min: 0.005, group: 'Aerodynamics' },
    { key: 'sfc_power', label: 'Fuel per unit of extra thrust power', unit: 'kg/(W·s)', default: 1e-7, min: 0, group: 'Aerodynamics', help: 'BSFC ÷ propeller efficiency, or TSFC ÷ V for jets; 0 for electric aircraft' },
    { key: 'A_cl', label: 'Lift-loss coefficient A', unit: '-', default: 1.8, min: 0.2, max: 5, group: 'Empirical model', help: 'ΔCLmax/CLmax = A·(k/c)^m. Defaults give about 8% at k/c = 10⁻⁴, 18% at 10⁻³ and 39% at 10⁻²' },
    { key: 'm_cl', label: 'Lift-loss exponent m', unit: '-', default: 0.333, min: 0.1, max: 1, group: 'Empirical model' },
    { key: 'cd_horn', label: 'Ice-shape drag coefficient', unit: '-', default: 1.0, min: 0, max: 3, group: 'Empirical model', help: 'Section Δcd = roughness term + this × (k/c) × (1 − 0.7·freezing fraction)' },
    { key: 'dcd_rough', label: 'Roughness drag increment', unit: '-', default: 0.002, min: 0, max: 0.01, group: 'Empirical model', help: 'Section drag added by a rough leading edge even when the ice is thin' },
  ],
  defaults: (c, up, d) => { const p = c.prop, jet = p.type === 'turbofan' || p.type === 'turbojet', alt = Math.min(c.atm.alt_m, 5000), V = c.flight.V_ms * Math.sqrt(isa(c.atm.alt_m).rho / isa(alt).rho);
    return { ice_mm: up.icing?.ice_max_mm ?? up.icing?.ice_thickness_mm, freezing_fraction: up.icing?.freezing_fraction, ice_mass_kg_m: up.icing?.ice_mass_kg_m, chord: d.mac || undefined, S: c.wing.S_m2, b: c.wing.b_m, tail_frac: (c.htail.S_m2 + c.vtail.S_m2) / c.wing.S_m2, mass_kg: c.mass.mtow_kg, V, alt_m: alt, CLmax: up.cfd?.CLmax ?? c.aero.CLmax_clean, CD0: up.cfd?.CD0 ?? c.aero.CD0, k: up.cfd?.k_induced ?? (d.k_induced || undefined), sfc_power: p.type === 'electric' ? 0 : jet ? p.tsfc_kg_Ns / Math.max(V, 1) : p.bsfc_kg_Ws / Math.max(p.eta_prop, 0.3) }; },
  run(i) {
    const a = isa(i.alt_m), q = 0.5 * a.rho * i.V ** 2, W = i.mass_kg * G0, kc = (i.ice_mm / 1e3) / i.chord, pn = penalty(kc, i.freezing_fraction, i.A_cl, i.m_cl, i.cd_horn, i.dcd_rough);
    const dCLmax = pn.dcl * Math.min(1, i.span_frac * 1.5), dCD = pn.dcd * i.span_frac * (1 + i.tail_frac), CL = W / (q * i.S), CDc = i.CD0 + i.k * CL * CL, CDi = CDc + dCD;
    const CLmaxI = i.CLmax * (1 - dCLmax), Vs = Math.sqrt((2 * W) / (a.rho * i.S * i.CLmax)), VsI = Math.sqrt((2 * W) / (a.rho * i.S * CLmaxI)), dP = dCD * q * i.S * i.V, iceMass = i.ice_mass_kg_m * i.b * i.span_frac * (1 + 0.5 * i.tail_frac);
    const hs = N.linspace(0, Math.max(30, i.ice_mm * 2), 50), warnings = ['These penalties come from a simple empirical correlation, not from a flow solution of the iced shape. Use them for trend and margin awareness only, and calibrate the coefficients with icing-tunnel or flight data for the section concerned.'];
    if (i.V < 1.2 * VsI) warnings.push(`Speed is only ${(i.V / VsI).toFixed(2)} × the iced stall speed: stall warning and protection systems calibrated for the clean wing may not activate in time.`);
    if (CL > CLmaxI) warnings.push('The iced wing cannot support the aircraft at this speed.');
    return {
      kpis: [
        { key: 'dCLmax_pct', label: 'Loss of maximum lift', value: 100 * dCLmax, unit: '%', status: dCLmax < 0.1 ? 'ok' : dCLmax < 0.25 ? 'warn' : 'bad' },
        { key: 'dCD_pct', label: 'Drag increase at this flight condition', value: (100 * dCD) / CDc, unit: '%', status: dCD / CDc < 0.2 ? 'ok' : dCD / CDc < 0.6 ? 'warn' : 'bad' },
        { key: 'CLmax_iced', label: 'Iced maximum lift coefficient', value: CLmaxI, unit: '-' },
        { key: 'V_stall_iced_ms', label: 'Iced stall speed', value: VsI, unit: 'm/s', note: `clean ${Vs.toFixed(1)} m/s` },
        { key: 'dV_stall_pct', label: 'Stall speed increase', value: 100 * (VsI / Vs - 1), unit: '%' },
        { key: 'stall_margin', label: 'Speed / iced stall speed', value: i.V / VsI, unit: '-', status: i.V / VsI > 1.3 ? 'ok' : i.V / VsI > 1.15 ? 'warn' : 'bad', note: 'Keep at least 1.3 in icing' },
        { key: 'LD_iced', label: 'Iced lift-to-drag ratio', value: CL / CDi, unit: '-', note: `clean ${(CL / CDc).toFixed(1)}` },
        { key: 'extra_power_W', label: 'Extra thrust power required', value: dP, unit: 'W' },
        { key: 'fuel_penalty_kg_h', label: 'Extra fuel burn', value: dP * i.sfc_power * 3600, unit: 'kg/h' },
        { key: 'co2_penalty_kg_h', label: 'Extra CO₂', value: dP * i.sfc_power * 3600 * 3.16, unit: 'kg/h' },
        { key: 'ice_mass_total_kg', label: 'Ice mass on wing and tail', value: iceMass, unit: 'kg' },
        { key: 'k_over_c', label: 'Ice height ratio k/c', value: kc, unit: '-' },
      ],
      plots: [
        { type: 'line', title: 'Estimated penalties against ice thickness (empirical)', xlabel: 'Ice thickness [mm]', ylabel: 'Change [%]', series: [{ name: 'Maximum-lift loss', x: hs, y: hs.map((h) => 100 * penalty(h / 1e3 / i.chord, i.freezing_fraction, i.A_cl, i.m_cl, i.cd_horn, i.dcd_rough).dcl * Math.min(1, i.span_frac * 1.5)) }, { name: 'Stall-speed increase', x: hs, y: hs.map((h) => 100 * (1 / Math.sqrt(1 - penalty(h / 1e3 / i.chord, i.freezing_fraction, i.A_cl, i.m_cl, i.cd_horn, i.dcd_rough).dcl * Math.min(1, i.span_frac * 1.5)) - 1)) }, { name: 'Drag increase', x: hs, y: hs.map((h) => (100 * penalty(h / 1e3 / i.chord, i.freezing_fraction, i.A_cl, i.m_cl, i.cd_horn, i.dcd_rough).dcd * i.span_frac * (1 + i.tail_frac)) / CDc) }], annotations: [{ x: i.ice_mm, label: 'This case' }] },
        { type: 'bar', title: 'Drag coefficient, clean and iced', ylabel: 'CD [-]', categories: ['Clean', 'Iced'], stacked: true, series: [{ name: 'Zero-lift', y: [i.CD0, i.CD0] }, { name: 'Induced', y: [i.k * CL * CL, i.k * CL * CL] }, { name: 'Ice', y: [0, dCD] }] },
      ],
      warnings,
      models: ['Power-law maximum-lift loss against ice height ratio (empirical, Brumby-type trend; default coefficients indicative)', 'Protuberance drag of the ice shape plus a roughness increment (empirical)', 'Parabolic drag polar for the clean aircraft'],
      assumptions: ['Penalties scale with the unprotected span fraction; the tail carries the same section penalty', 'Lift-curve slope and trim changes, tailplane stall and control hinge-moment effects are not modelled', 'Added ice weight is reported but not included in the stall speed'],
    };
  },
  calibration: { params: [{ key: 'A_cl', min: 0.5, max: 4 }, { key: 'm_cl', min: 0.15, max: 0.6 }, { key: 'cd_horn', min: 0.2, max: 2.5 }, { key: 'dcd_rough', min: 0, max: 0.006 }], sweep: 'ice_mm', target: 'dCLmax_pct', note: 'Supply measured maximum-lift loss (or drag rise, target dCD_pct) against ice thickness from icing-tunnel or flight tests of the same section.' },
  verify() {
    const b = { ice_mm: 10, freezing_fraction: 0, ice_mass_kg_m: 1, span_frac: 1, chord: 1, S: 10, b: 10, tail_frac: 0, mass_kg: 1000, V: 60, alt_m: 0, CLmax: 1.5, CD0: 0.03, k: 0.05, sfc_power: 1e-7, A_cl: 1.8, m_cl: 1 / 3, cd_horn: 1, dcd_rough: 0 };
    const o = N.kv(degradation.run(b)), z = N.kv(degradation.run({ ...b, ice_mm: 0 }));
    return [
      N.check('Stall-speed ratio = 1/√(1 − ΔCLmax)', 1 + o.dV_stall_pct / 100, 1 / Math.sqrt(1 - o.dCLmax_pct / 100), 1e-12, 'Lift equation'),
      N.check('Extra power = ΔCD·q·S·V', o.extra_power_W, 0.01 * 0.5 * isa(0).rho * 3600 * 10 * 60, 1e-12, 'Drag power'),
      N.check('Clean limit: no penalty without ice', z.dCLmax_pct + z.dCD_pct, 0, 1e-12, 'Limiting behaviour'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.stall_margin < 1.3) out.push({ severity: o.stall_margin < 1.15 ? 'critical' : 'warn', title: 'Reduced stall margin with ice', detail: `Iced stall speed ${o.V_stall_iced_ms.toFixed(1)} m/s (+${o.dV_stall_pct.toFixed(0)}%); margin ${o.stall_margin.toFixed(2)}.`, action: 'Raise minimum manoeuvring and approach speeds, limit bank angle and flap extension, and leave the icing layer.', basis: '1.3 × stall speed margin with the estimated iced CLmax' });
    if (o.dCD_pct > 20) out.push({ severity: 'advise', title: 'Significant drag and fuel penalty', detail: `Drag +${o.dCD_pct.toFixed(0)}%, about ${o.fuel_penalty_kg_h.toFixed(1)} kg/h more fuel (${o.co2_penalty_kg_h.toFixed(1)} kg/h CO₂) and ${o.ice_mass_total_kg.toFixed(0)} kg of ice.`, action: 'Compare this with the ice-protection power in the protection analysis: running the system is normally the lower-fuel option; route or level changes to avoid icing save both.', basis: 'Drag power and specific fuel consumption' });
    out.push({ severity: 'info', title: 'Empirical estimate', detail: 'The lift and drag penalties are correlations with large scatter between ice shapes.', action: 'Calibrate the four coefficients against tunnel data, or hand the ice shape to a RANS solver for the iced-section polar.', basis: 'Model limitation' });
    return out;
  },
};

// ---- ice protection power -------------------------------------------------------------------------------
function protect(i) {
  const a = air(i.alt_m, i.T_C), R = radius(i), q = inertia(i.V, i.mvd * 1e-6, R, a), tr = trajectories(q.K, q.Red, 16), sH = (i.heated_pct / 100) * i.chord, A1 = 2 * sH, lwc = i.lwc / 1e3;
  const Res = (a.rho * i.V * sH) / a.mu, hPlate = (a.k / sH) * 0.037 * Res ** 0.8 * PR_AIR ** (1 / 3), hC = hCyl(0.6, i.V, 2 * R, a, 1), h = Math.max(hPlate, hC) * i.k_h;
  const mCatch = tr.E * lwc * i.V * 2 * R, mA = mCatch / A1, Trec = a.T + (REC * i.V ** 2) / (2 * CP_AIR), hG = (h / CP_AIR) * (PR_AIR / SC) ** (2 / 3), evap = (Ts) => Math.max(0, (hG * 0.622 * (pSat(Ts) - pSat(a.T))) / a.p);
  const Ts = i.T_surf_C + TF, mE = Math.min(evap(Ts), mA), qRW = Math.max(0, h * (Ts - Trec) + mA * (C_W * (Ts - a.T) - 0.5 * i.V ** 2) + mE * L_V);
  const TsE = mA > 0 ? N.findRoot((t) => evap(t) - mA, TF, 372, 100, 1e-9) : TF, tE = Number.isFinite(TsE) ? Math.max(TsE, Ts) : 372, qEV = Math.max(0, h * (tE - Trec) + mA * (C_W * (tE - a.T) - 0.5 * i.V ** 2 + L_V));
  return { a, R, q, tr, A1, h, mCatch, mA, Trec, qRW, qEV, TsE: tE, mE, runback: Math.max(0, mA - mE) * A1 };
}
const protection = {
  id: 'protection', title: 'Anti-icing and de-icing power', fidelity: 'reduced-order',
  summary: 'Heat needed to keep the protected leading edge free of ice, either running wet or fully evaporative, the average power of a cyclic electro-thermal de-icing system, the equivalent bleed-air flow, and what it costs in fuel or battery energy.',
  equations: ['Convective heat-transfer equations', 'Evaporation equations', 'Energy conservation equations for ice and liquid water', 'Electrothermal anti-icing models', 'Electrothermal ice protection–thermal flow equations', 'Pneumatic de-icing models'],
  inputs: [...CLOUD,
    { key: 'heated_pct', label: 'Heated extent per surface', unit: '% chord', default: 8, min: 1, max: 30, group: 'Protected zone', help: 'Chordwise heated length on the upper and on the lower surface' },
    { key: 'span_prot', label: 'Protected span (all surfaces)', unit: 'm', default: 8, min: 0.05, group: 'Protected zone', help: 'Total protected leading-edge length: wing, tail, blades' },
    { key: 'T_surf_C', label: 'Running-wet surface temperature', unit: '°C', default: 5, min: 1, max: 60, group: 'Protected zone' },
    { key: 'k_h', label: 'Heat-transfer multiplier', unit: '-', default: 1.2, min: 1, max: 3, group: 'Protected zone', help: 'Allowance for roughness and unsteady flow on the heated zone' },
    { key: 'eta_heater', label: 'Heater to surface efficiency', unit: '-', default: 0.85, min: 0.3, max: 1, group: 'System', help: 'Fraction of electrical or bleed heat that reaches the outer surface' },
    { key: 'q_deice', label: 'De-icing heater flux', unit: 'W/m²', default: 30000, min: 2000, max: 80000, group: 'De-icing cycle', help: 'Cycled zones are heated intensely for a short time' },
    { key: 't_on', label: 'Heater on-time per zone', unit: 's', default: 8, min: 1, max: 120, group: 'De-icing cycle' },
    { key: 't_cycle', label: 'Cycle period', unit: 's', default: 120, min: 10, max: 600, group: 'De-icing cycle' },
    { key: 'ps_frac', label: 'Parting-strip fraction of heated area', unit: '-', default: 0.1, min: 0, max: 0.5, group: 'De-icing cycle', help: 'Continuously heated strip that splits the ice cap' },
    { key: 'rho_ice', label: 'Ice density', unit: 'kg/m³', default: 900, min: 300, max: 917, group: 'De-icing cycle' },
    { key: 'dT_bleed', label: 'Bleed-air temperature drop in the piccolo', unit: 'K', default: 120, min: 20, max: 300, group: 'System' },
    { key: 'sfc_elec', label: 'Fuel per unit of power off-take', unit: 'kg/(W·s)', default: 9e-8, min: 0, group: 'System', help: 'About 8–10 × 10⁻⁸ for shaft power from a turbine engine; 0 for battery aircraft' },
    { key: 'batt_kWh', label: 'Battery energy (electric aircraft)', unit: 'kWh', default: 0, min: 0, group: 'System' }],
  defaults: (c, up, d) => { const rot = !(c.wing.S_m2 > 0);
    return { ...cloudDefaults(c, d), span_prot: rot ? 0.8 * c.rotor.R_m * c.rotor.n_blades * (c.meta.type === 'helicopter' ? 1 : c.prop.n_eng) || undefined : 0.75 * c.wing.b_m + 0.8 * c.htail.b_m, sfc_elec: c.prop.type === 'electric' ? 0 : undefined, batt_kWh: c.prop.type === 'electric' ? c.systems.batt_kWh : undefined, heated_pct: rot ? 12 : undefined }; },
  run(i) {
    const r = protect(i), A = r.A1 * i.span_prot, Prw = (r.qRW * A) / i.eta_heater, Pev = (r.qEV * A) / i.eta_heater, duty = N.clamp(i.t_on / i.t_cycle, 0, 1);
    const Pde = ((i.q_deice * A * (1 - i.ps_frac) * duty + r.qRW * A * i.ps_frac) / i.eta_heater), rate = (r.tr.beta0 * (i.lwc / 1e3) * i.V) / i.rho_ice, tIC = rate * (i.t_cycle - i.t_on);
    const eNeed = i.rho_ice * (0.5 * tIC * C_I * Math.max(TF - r.a.T, 0) + 1e-4 * (L_F + C_I * Math.max(TF - r.a.T, 0))), eHave = i.q_deice * i.t_on * i.eta_heater, bleed = Prw * i.eta_heater / (0.6 * CP_AIR * i.dT_bleed), warnings = [];
    sldNote(i, warnings);
    if (r.runback > 0) warnings.push(`Running wet, ${(r.runback * i.span_prot * 3600).toFixed(1)} kg/h of water leaves the heated zone and can refreeze behind it as a runback ridge.`);
    if (r.qEV > 4e4) warnings.push('The fully evaporative heat flux exceeds 40 kW/m², beyond what electro-thermal mats normally deliver; this condition is a running-wet or de-icing case.');
    if (eHave < eNeed) warnings.push('The de-icing pulse delivers less energy than is needed to warm the ice layer and melt the bond line: lengthen the on-time or raise the heater flux.');
    if (r.Trec > TF) warnings.push('The recovery temperature is above freezing: no ice protection heat is required at this speed and temperature.');
    const Ts = N.linspace(-30, -1, 30), sw = Ts.map((t) => protect({ ...i, T_C: t })), fuel = (P) => P * i.sfc_elec * 3600;
    return {
      kpis: [
        { key: 'antiice_power_W', label: 'Running-wet anti-icing power', value: Prw, unit: 'W' },
        { key: 'evap_power_W', label: 'Fully evaporative anti-icing power', value: Pev, unit: 'W' },
        { key: 'deice_power_W', label: 'Cyclic de-icing average power', value: Pde, unit: 'W', status: Pde < Prw ? 'ok' : 'warn', note: 'Includes the parting strip' },
        { key: 'q_running_wet_Wm2', label: 'Running-wet surface heat flux', value: r.qRW, unit: 'W/m²' },
        { key: 'q_evaporative_Wm2', label: 'Evaporative surface heat flux', value: r.qEV, unit: 'W/m²', note: `surface at ${(r.TsE - TF).toFixed(0)} °C` },
        { key: 'heated_area_m2', label: 'Heated area', value: A, unit: 'm²' },
        { key: 'water_catch_kg_h', label: 'Water catch on the protected span', value: r.mCatch * i.span_prot * 3600, unit: 'kg/h' },
        { key: 'runback_kg_h', label: 'Runback when running wet', value: r.runback * i.span_prot * 3600, unit: 'kg/h', status: r.runback > 0 ? 'warn' : 'ok' },
        { key: 'intercycle_ice_mm', label: 'Inter-cycle ice thickness', value: tIC * 1e3, unit: 'mm', status: tIC < 0.003 ? 'ok' : 'warn', note: 'Rime growth between heater pulses; keep to a few millimetres' },
        { key: 'deice_energy_ratio', label: 'De-icing pulse energy / required', value: eHave / eNeed, unit: '-', status: eHave >= eNeed ? 'ok' : 'bad' },
        { key: 'bleed_kgs', label: 'Equivalent bleed-air flow (running wet)', value: bleed, unit: 'kg/s' },
        { key: 'fuel_penalty_kg_h', label: 'Fuel for running-wet power', value: fuel(Prw), unit: 'kg/h' },
        { key: 'co2_kg_h', label: 'CO₂ for running-wet power', value: fuel(Prw) * 3.16, unit: 'kg/h' },
        { key: 'battery_pct_per_h', label: 'Battery energy used per hour (running wet)', value: i.batt_kWh > 0 ? (100 * Prw) / 1e3 / i.batt_kWh : 0, unit: '%/h', status: i.batt_kWh > 0 && Prw / 1e3 / i.batt_kWh > 0.15 ? 'warn' : 'ok' },
      ],
      plots: [
        { type: 'line', title: 'Anti-icing heat flux against air temperature', xlabel: 'Static air temperature [°C]', ylabel: 'Surface heat flux [kW/m²]', series: [{ name: 'Running wet', x: Ts, y: sw.map((s) => s.qRW / 1e3) }, { name: 'Fully evaporative', x: Ts, y: sw.map((s) => s.qEV / 1e3) }], annotations: [{ x: i.T_C, label: 'This case' }] },
        { type: 'bar', title: 'Ice-protection power by operating mode', ylabel: 'Power [kW]', categories: ['Cyclic de-icing', 'Running wet', 'Fully evaporative'], series: [{ name: 'Average power', y: [Pde / 1e3, Prw / 1e3, Pev / 1e3] }] },
      ],
      outputs: { h_protected_Wm2K: r.h },
      warnings,
      models: ['Steady surface energy balance: convection, droplet sensible heating, kinetic heating and evaporation', 'Chilton–Colburn analogy for evaporation', 'Flat-plate turbulent or cylinder laminar heat transfer, whichever is larger (empirical correlations)', 'Duty-cycle average for electro-thermal de-icing with a parting strip'],
      assumptions: ['Heated zone treated as one control volume at uniform temperature', 'Water catch from the leading-edge cylinder at zero incidence', 'De-icing bond-line criterion: warm half the ice layer and melt 0.1 mm at the interface; skin thermal mass neglected', 'Bleed estimate assumes 60% piccolo effectiveness'],
    };
  },
  calibration: { params: [{ key: 'k_h', min: 1, max: 3 }, { key: 'eta_heater', min: 0.5, max: 1 }], sweep: 'T_C', target: 'q_running_wet_Wm2', note: 'Supply measured heater flux needed to hold the surface temperature against air temperature in an icing tunnel.' },
  verify() {
    const b = { V: 80, alt_m: 0, T_C: -10, lwc: 0, mvd: 20, chord: 1, tc: 0.12, d_le: 0, heated_pct: 10, span_prot: 5, T_surf_C: 5, k_h: 1, eta_heater: 1, q_deice: 20000, t_on: 10, t_cycle: 100, ps_frac: 0, rho_ice: 900, dT_bleed: 120, sfc_elec: 0, batt_kWh: 0 };
    const dry = protection.run(b), o = N.kv(dry), h = dry.outputs.h_protected_Wm2K, Trec = 263.15 + (0.85 * 6400) / (2 * CP_AIR), wet = N.kv(protection.run({ ...b, lwc: 0.5 }));
    return [
      N.check('Dry air: power = h·A·(Ts − Trec)', o.antiice_power_W, h * 1 * (278.15 - Trec), 1e-10, 'Newton cooling of the heated zone'),
      N.check('De-icing average power = q·A·t_on/t_cycle', o.deice_power_W, 20000 * 1 * 0.1, 1e-12, 'Duty cycle'),
      N.check('Evaporative power exceeds running-wet power', wet.evap_power_W > wet.antiice_power_W ? 1 : 0, 1, 1e-12, 'Latent heat of the full catch'),
      N.check('Wet running adds droplet heating and evaporation', wet.antiice_power_W > o.antiice_power_W ? 1 : 0, 1, 1e-12, 'Energy balance ordering'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    out.push({ severity: 'info', title: 'Choose the protection mode by energy', detail: `Cyclic de-icing ${(o.deice_power_W / 1e3).toFixed(1)} kW, running wet ${(o.antiice_power_W / 1e3).toFixed(1)} kW, fully evaporative ${(o.evap_power_W / 1e3).toFixed(1)} kW.`, action: 'Use evaporative anti-icing only where no runback or shed ice can be tolerated (engine intakes, probes); de-ice the wing and tail cyclically to cut power off-take, fuel and CO₂.', basis: 'Surface energy balance' });
    if (o.runback_kg_h > 0) out.push({ severity: 'advise', title: 'Runback water leaves the heated zone', detail: `${o.runback_kg_h.toFixed(1)} kg/h.`, action: 'Extend the heated zone, raise the surface temperature, or accept and periodically shed a runback ridge; check its aerodynamic effect.', basis: 'Water mass balance' });
    if (o.intercycle_ice_mm > 3) out.push({ severity: 'warn', title: 'Thick inter-cycle ice', detail: `${o.intercycle_ice_mm.toFixed(1)} mm builds between pulses.`, action: 'Shorten the cycle period; inter-cycle ice still costs lift and drag.', basis: 'Accretion rate × off-time' });
    if (i.batt_kWh > 0 && o.battery_pct_per_h > 15) out.push({ severity: 'warn', title: 'Ice protection is a major battery load', detail: `${o.battery_pct_per_h.toFixed(0)}% of the battery per hour when running wet.`, action: 'Prefer cyclic de-icing, hydrophobic or low-adhesion coatings and operational avoidance; include this load in the Suite 19 energy budget and reserves.', basis: 'Energy budget' });
    return out;
  },
};

// ---- icing-condition envelope map -------------------------------------------------------------------------
const stagPoint = (i, T_C, lwc_g) => { const a = air(i.alt_m, T_C), R = radius(i), q = inertia(i.V, i.mvd * 1e-6, R, a); return { q, m: messinger({ beta: q.beta0, mIn: 0, h: hCyl(0, i.V, 2 * R, a, i.k_h), V: i.V, T: a.T, p: a.p, lwc: lwc_g / 1e3, rh: 1, q: 0 }) }; };
const envelope = {
  id: 'envelope', title: 'Icing severity map over temperature and water content', fidelity: 'reduced-order',
  summary: 'Sweeps the stagnation-line energy balance over air temperature and cloud liquid water content to map where ice forms, how fast, and where rime gives way to glaze, so an operating or certification envelope can be overlaid.',
  equations: ['Messinger icing energy balance equations', 'Droplet impingement equations', 'Evaporation equations', 'Ice growth equations'],
  inputs: [...CLOUD,
    { key: 'k_h', label: 'Roughness heat-transfer multiplier', unit: '-', default: 1.5, min: 1, max: 4, group: 'Model' },
    { key: 'rho_ice', label: 'Ice density', unit: 'kg/m³', default: 900, min: 300, max: 917, group: 'Model' },
    { key: 'lwc_max', label: 'Largest water content on the map', unit: 'g/m³', default: 3, min: 0.2, max: 5, group: 'Numerics' },
    { key: 'T_min_C', label: 'Coldest temperature on the map', unit: '°C', default: -30, min: -40, max: -5, group: 'Numerics' },
    { key: 'nGrid', label: 'Map points per axis', unit: '', default: 25, min: 8, max: 60, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => cloudDefaults(c, d),
  run(i) {
    const n = Math.round(N.clamp(i.nGrid, 8, 60)), Ts = N.linspace(i.T_min_C, 2, n), Ls = N.linspace(i.lwc_max / n, i.lwc_max, n), rate = (m) => (m.mIce / i.rho_ice) * 6e4;
    const grid = Ls.map((l) => Ts.map((t) => stagPoint(i, t, l).m)), op = stagPoint(i, i.T_C, i.lwc), warnings = [];
    // Ludlam limit: water content above which not all impinging water can freeze, at each temperature
    const lud = Ts.map((t) => { const f = (l) => { const m = stagPoint(i, t, l).m; return (m.mIce > 0 ? m.n : 0) - 0.999; }; return f(0.01) > 0 ? N.findRoot(f, 0.01, 20, 400, 1e-6) : 0; });
    const ludOp = N.interp1(Ts, lud.map((v) => (Number.isFinite(v) ? v : 20)), i.T_C), tRime = N.findRoot((t) => { const m = stagPoint(i, t, i.lwc).m; return (m.mIce > 0 ? m.n : 0) - 0.999; }, -60, 2, 124, 1e-6);
    const tOn = N.findRoot((t) => stagPoint(i, t, i.lwc).m.mIce - 1e-12, i.T_min_C, 5, 70, 1e-7), rmax = Math.max(...grid.flat().map(rate));
    sldNote(i, warnings);
    if (op.m.mIce <= 0) warnings.push('No ice forms at the operating point.');
    const j = Ts.map((t, k) => k).filter((k) => Number.isFinite(lud[k]) && lud[k] > 0 && lud[k] <= i.lwc_max);
    return {
      kpis: [
        { key: 'rate_mm_min', label: 'Growth rate at the operating point', value: rate(op.m), unit: 'mm/min' },
        { key: 'freezing_fraction_op', label: 'Freezing fraction at the operating point', value: op.m.mIce > 0 ? op.m.n : 0, unit: '-' },
        { key: 'ludlam_lwc_gm3', label: 'Rime limit water content at this temperature', value: ludOp, unit: 'g/m³', note: 'Above this (Ludlam limit) the surface runs wet and glaze forms' },
        { key: 'T_rime_C', label: 'Temperature below which the ice is rime', value: tRime, unit: '°C', note: 'at this water content' },
        { key: 'T_onset_C', label: 'Warmest temperature with ice', value: tOn, unit: '°C', note: 'Set by kinetic heating and evaporative cooling' },
        { key: 'rate_max_mm_min', label: 'Largest growth rate on the map', value: rmax, unit: 'mm/min' },
        { key: 'time_to_5mm_min', label: 'Time to 5 mm of ice', value: rate(op.m) > 0 ? 5 / rate(op.m) : Infinity, unit: 'min', status: rate(op.m) > 0 && 5 / rate(op.m) < 5 ? 'warn' : 'ok' },
        { key: 'beta0', label: 'Stagnation collection efficiency (fit)', value: op.q.beta0, unit: '-' },
      ],
      plots: [
        { type: 'heat', title: 'Stagnation ice growth rate', xlabel: 'Static air temperature [°C]', ylabel: 'Liquid water content [g/m³]', zlabel: 'Growth rate [mm/min]', x: Ts, y: Ls, z: grid.map((row) => row.map(rate)), contours: 12, overlay: [{ name: 'Rime / glaze boundary (Ludlam limit)', x: j.map((k) => Ts[k]), y: j.map((k) => lud[k]) }, { name: 'Operating point', x: [i.T_C], y: [i.lwc] }] },
        { type: 'heat', title: 'Freezing fraction', xlabel: 'Static air temperature [°C]', ylabel: 'Liquid water content [g/m³]', zlabel: 'Freezing fraction [-]', x: Ts, y: Ls, z: grid.map((row) => row.map((m) => (m.mIce > 0 ? m.n : 0))), contours: 10, overlay: [{ name: 'Operating point', x: [i.T_C], y: [i.lwc] }] },
        { type: 'line', title: 'Growth rate against temperature at this water content', xlabel: 'Static air temperature [°C]', ylabel: 'Growth rate [mm/min]', series: [{ name: `LWC ${i.lwc} g/m³`, x: Ts, y: Ts.map((t) => rate(stagPoint(i, t, i.lwc).m)) }], annotations: [{ x: i.T_C, label: 'This case' }] },
      ],
      warnings,
      models: ['Stagnation-line Messinger balance', 'Langmuir–Blodgett stagnation collection efficiency (empirical fit, at the MVD)', 'Frössling stagnation heat transfer with a roughness multiplier'],
      assumptions: ['Clean leading-edge geometry (no growth feedback) and saturated air', 'Single droplet size at the MVD', 'The certification envelope itself is not built in: overlay the applicable design LWC–temperature–MVD points yourself'],
    };
  },
  verify() {
    const i = { V: 80, alt_m: 0, T_C: -20, lwc: 0.2, mvd: 20, chord: 1, tc: 0.12, d_le: 0.05, k_h: 1, rho_ice: 900, lwc_max: 2, T_min_C: -30, nGrid: 10 }, o = N.kv(envelope.run(i)), q = inertia(80, 20e-6, 0.025, air(0, -20));
    const atL = stagPoint(i, -20, o.ludlam_lwc_gm3).m;
    return [
      N.check('Rime growth rate = β0·LWC·V/ρ_ice', o.rate_mm_min, ((q.beta0 * 0.2e-3 * 80) / 900) * 6e4, 0.04, 'Cold limit of the Messinger balance; sublimation removes 2–3%'),
      N.check('Freezing fraction is 1 at the Ludlam limit', atL.n, 1, 3e-3, 'Definition of the limit'),
      N.check('No ice when the recovery temperature exceeds 0 °C', stagPoint({ ...i, V: 150 }, -5, 0.3).m.mIce, 0, 1e-12, 'Kinetic heating: Trec = −5 + 9.5 °C'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.freezing_fraction_op > 0 && o.freezing_fraction_op < 0.999) out.push({ severity: 'advise', title: 'Operating point is in the glaze regime', detail: `Water content ${i.lwc} g/m³ exceeds the rime limit ${o.ludlam_lwc_gm3.toFixed(2)} g/m³ at ${i.T_C} °C.`, action: `Below about ${Number.isFinite(o.T_rime_C) ? o.T_rime_C.toFixed(0) : '−40'} °C the same cloud gives rime, which is less harmful; above ${Number.isFinite(o.T_onset_C) ? o.T_onset_C.toFixed(1) : '0'} °C no ice forms. Use the map to choose an exit altitude.`, basis: 'Ludlam limit' });
    if (o.time_to_5mm_min < 5) out.push({ severity: 'warn', title: 'Rapid ice build-up', detail: `5 mm in ${o.time_to_5mm_min.toFixed(1)} min.`, action: 'Activate ice protection before entering cloud and minimise time in the layer.', basis: 'Accretion rate' });
    out.push({ severity: 'info', title: 'Overlay the design envelope', detail: 'The map covers the entered temperature and water-content range for one droplet size and speed.', action: 'Check each design point of the applicable certification icing envelope (water content falls as droplet size and cloud extent rise) and repeat at holding, climb and approach speeds.', basis: 'Certification practice' });
    return out;
  },
};

// ---- rotor and propeller blade icing with centrifugal shedding ---------------------------------------------
const rotorIce = {
  id: 'rotorice', title: 'Rotor and propeller blade icing and shedding', fidelity: 'reduced-order',
  summary: 'Ice growth along a rotating blade, where speed, kinetic heating and collection efficiency all change with radius; the radius where centrifugal force sheds the ice; and the resulting power rise and out-of-balance force.',
  equations: ['Messinger icing energy balance equations', 'Droplet impingement equations', 'Convective heat-transfer equations', 'Ice growth equations', 'Iced-airfoil aerodynamic models'],
  applicable: (c) => (c.rotor.R_m > 0 || c.prop.prop_dia_m > 0 ? true : 'No rotor or propeller is defined; turbofan fan and intake icing is not modelled natively.'),
  inputs: [
    { key: 'R', label: 'Blade radius', unit: 'm', default: 8, min: 0.05, group: 'Rotor' },
    { key: 'chord', label: 'Blade chord', unit: 'm', default: 0.5, min: 0.005, group: 'Rotor' },
    { key: 'tc', label: 'Blade thickness ratio', unit: '-', default: 0.12, min: 0.04, max: 0.25, group: 'Rotor' },
    { key: 'rpm', label: 'Rotational speed', unit: 'rpm', default: 260, min: 10, group: 'Rotor' },
    { key: 'n_blades', label: 'Blades', unit: '', default: 4, min: 2, max: 12, step: 1, discrete: true, group: 'Rotor' },
    { key: 'cd0', label: 'Clean blade profile drag coefficient', unit: '-', default: 0.01, min: 0.004, max: 0.05, group: 'Rotor' },
    ...CLOUD.filter((f) => ['alt_m', 'T_C', 'lwc', 'mvd'].includes(f.key)),
    { key: 'time_min', label: 'Exposure time', unit: 'min', default: 5, min: 0.1, max: 60, group: 'Flight condition' },
    { key: 'tau_adh', label: 'Ice adhesion shear strength', unit: 'Pa', default: 3e5, min: 1e4, max: 2e6, group: 'Model', help: 'Strongly dependent on temperature, surface and ice type: roughly 0.1–0.6 MPa on metal, much lower on icephobic coatings' },
    { key: 'rho_ice', label: 'Ice density', unit: 'kg/m³', default: 900, min: 300, max: 917, group: 'Model' },
    { key: 'k_h', label: 'Roughness heat-transfer multiplier', unit: '-', default: 1.5, min: 1, max: 4, group: 'Model' },
    { key: 'cd_horn', label: 'Ice-shape drag coefficient', unit: '-', default: 1.0, min: 0, max: 3, group: 'Model', help: 'Empirical: Δcd = 0.002 + this × (k/c) × (1 − 0.7·freezing fraction)' },
    { key: 'nR', label: 'Radial stations', unit: '', default: 30, min: 8, max: 200, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c) => { const rot = c.rotor.R_m > 0, alt = Math.min(c.atm.alt_m, 5000), R = rot ? c.rotor.R_m : c.prop.prop_dia_m / 2;
    return { R, chord: rot ? c.rotor.chord_m : 0.15 * R, rpm: (rot ? c.rotor.rpm : c.prop.rpm) || undefined, n_blades: (rot ? c.rotor.n_blades : c.prop.n_blades) || undefined, cd0: rot ? c.rotor.cd0 : undefined, alt_m: alt, T_C: N.clamp(isa(alt, c.atm.dISA_K).T - TF, -25, -5) }; },
  run(i) {
    const Om = (i.rpm * 2 * PI) / 60, n = Math.round(N.clamp(i.nR, 8, 200)), rr = N.linspace(0.15, 1, n).map((x) => x * i.R), Rle = rLE(i.tc, i.chord), a = air(i.alt_m, i.T_C), tSec = i.time_min * 60;
    const st = rr.map((r) => { const V = Math.max(Om * r, 1), q = inertia(V, i.mvd * 1e-6, Rle, a), m = messinger({ beta: q.beta0, mIn: 0, h: hCyl(0, V, 2 * Rle, a, i.k_h), V, T: a.T, p: a.p, lwc: i.lwc / 1e3, rh: 1, q: 0 }), rate = m.mIce / i.rho_ice, tShed = i.tau_adh / (i.rho_ice * Om * Om * r);
      const tk = Math.min(rate * tSec, tShed), pn = penalty(tk / i.chord, m.mIce > 0 ? m.n : 1, 1.8, 1 / 3, i.cd_horn, 0.002); return { r, V, m, rate, tShed, tk, shed: rate * tSec >= tShed, tToShed: rate > 0 ? tShed / rate : Infinity, pn }; });
    const w3 = st.map((s) => s.r ** 3), dP = N.trapz(rr, st.map((s, j) => s.pn.dcd * w3[j])) / (i.cd0 * N.trapz(rr, w3)), iced = st.filter((s) => s.tk > 0), free = st.find((s) => s.rate <= 0 && s.m.Trec > TF), shedSt = st.filter((s) => s.shed);
    const tFirst = Math.min(...st.map((s) => s.tToShed)), jMax = N.argmax(st.map((s) => s.tk)), wIce = 0.08 * i.chord, mBlade = N.trapz(rr, st.map((s) => i.rho_ice * s.tk * wIce));
    // out-of-balance force if one blade sheds its outboard ice (beyond the first shedding radius) and the others do not
    const jS = st.findIndex((s) => s.tToShed === tFirst), Fimb = Number.isFinite(tFirst) && jS >= 0 ? N.trapz(rr.slice(jS), st.slice(jS).map((s) => i.rho_ice * Math.min(s.rate * tFirst, s.tShed) * wIce * Om * Om * s.r)) : 0, warnings = ['Blade drag and lift penalties use the same empirical correlation as the wing analysis; calibrate before relying on the power figure.'];
    if (i.mvd > 50) warnings.push('Supercooled large droplets: impingement aft of the blade leading-edge protection is not modelled.');
    if (shedSt.length) warnings.push(`Ice sheds outboard of r/R = ${(shedSt[0].r / i.R).toFixed(2)} within the exposure time. Shedding is rarely symmetric: expect vibration and ice impact on the fuselage, tail rotor or pusher propeller.`);
    if (Om * i.R / Math.sqrt(1.4 * R_AIR * a.T) > 0.75) warnings.push('Tip Mach number above 0.75: compressibility raises the true tip recovery temperature and heat transfer beyond this incompressible estimate.');
    return {
      kpis: [
        { key: 'ice_max_mm', label: 'Maximum blade ice thickness', value: st[jMax].tk * 1e3, unit: 'mm', note: `at r/R = ${(st[jMax].r / i.R).toFixed(2)}` },
        { key: 'r_ice_free_frac', label: 'Radius beyond which kinetic heating prevents ice', value: free ? free.r / i.R : 1, unit: 'r/R', note: '1 = iced to the tip' },
        { key: 'r_shed_frac', label: 'Innermost shedding radius', value: shedSt.length ? shedSt[0].r / i.R : 1, unit: 'r/R', note: '1 = no shedding within the exposure' },
        { key: 't_first_shed_min', label: 'Time to first shedding', value: tFirst / 60, unit: 'min' },
        { key: 'shed_thickness_tip_mm', label: 'Shedding thickness at the tip', value: st[n - 1].tShed * 1e3, unit: 'mm' },
        { key: 'profile_power_increase_pct', label: 'Profile power increase', value: 100 * dP, unit: '%', status: dP < 0.2 ? 'ok' : dP < 0.6 ? 'warn' : 'bad', note: 'Empirical' },
        { key: 'ice_mass_blade_kg', label: 'Ice mass per blade', value: mBlade, unit: 'kg' },
        { key: 'imbalance_force_N', label: 'Out-of-balance force after one-blade shed', value: Fimb, unit: 'N', note: 'Rotating 1/rev hub force' },
        { key: 'iced_span_frac', label: 'Fraction of stations carrying ice', value: iced.length / n, unit: '-' },
        { key: 'T_recovery_tip_C', label: 'Tip recovery temperature', value: st[n - 1].m.Trec - TF, unit: '°C' },
        { key: 'tip_speed_ms', label: 'Tip speed', value: Om * i.R, unit: 'm/s' },
      ],
      plots: [
        { type: 'line', title: 'Ice thickness along the blade', xlabel: 'Radius r/R [-]', ylabel: 'Thickness [mm]', series: [{ name: `Ice after ${i.time_min} min`, x: rr.map((r) => r / i.R), y: st.map((s) => s.tk * 1e3) }, { name: 'Shedding thickness', x: rr.map((r) => r / i.R), y: st.map((s) => Math.min(s.tShed * 1e3, 4 * Math.max(st[jMax].tk * 1e3, 1))), style: 'dash' }] },
        { type: 'line', title: 'Freezing fraction and collection efficiency', xlabel: 'Radius r/R [-]', ylabel: 'Fraction [-]', series: [{ name: 'Freezing fraction', x: rr.map((r) => r / i.R), y: st.map((s) => (s.m.mIce > 0 ? s.m.n : 0)) }, { name: 'Collection efficiency β0', x: rr.map((r) => r / i.R), y: st.map((s) => s.m.mImp / Math.max((i.lwc / 1e3) * s.V, 1e-30)) }] },
        { type: 'line', title: 'Leading-edge temperatures', xlabel: 'Radius r/R [-]', ylabel: 'Temperature [°C]', series: [{ name: 'Recovery temperature', x: rr.map((r) => r / i.R), y: st.map((s) => s.m.Trec - TF) }, { name: 'Surface temperature', x: rr.map((r) => r / i.R), y: st.map((s) => s.m.Ts - TF) }], annotations: [{ y: 0, label: 'Freezing' }] },
      ],
      outputs: { dCD_pct: 100 * dP, dCLmax_pct: 100 * Math.max(...st.map((s) => s.pn.dcl)) },
      warnings,
      models: ['Stagnation-line Messinger balance at each radius with the local blade speed', 'Langmuir–Blodgett stagnation collection efficiency (empirical fit)', 'Centrifugal shedding when ρ_ice·t·Ω²·r exceeds the adhesion shear strength', 'Empirical iced-section drag penalty weighted by r³ for profile power'],
      assumptions: ['Hover or static operation: section speed is Ω·r; forward speed and inflow are neglected', 'Ice cohesion between neighbouring stations is ignored, so shedding is local and conservative in time', 'Ice strip width is 8% of chord for the mass and imbalance estimates', 'Adhesion strength is a single user value'],
    };
  },
  convergence: { param: 'nR', label: 'Radial stations', levels: [10, 20, 40, 80], metric: 'ice_mass_blade_kg' },
  calibration: { params: [{ key: 'tau_adh', min: 3e4, max: 1e6 }, { key: 'k_h', min: 1, max: 4 }, { key: 'cd_horn', min: 0.2, max: 2.5 }], sweep: 'time_min', target: 'r_shed_frac', note: 'Supply observed shedding radius (or torque rise, target profile_power_increase_pct) against exposure time from a rotor icing test.' },
  verify() {
    const i = { R: 5, chord: 0.4, tc: 0.12, rpm: 300, n_blades: 4, cd0: 0.01, alt_m: 0, T_C: -20, lwc: 0.3, mvd: 20, time_min: 2, tau_adh: 3e5, rho_ice: 900, k_h: 1, cd_horn: 1, nR: 40 }, o = N.kv(rotorIce.run(i)), Om = 10 * PI;
    return [
      N.check('Tip shedding thickness = τ/(ρ·Ω²·R)', o.shed_thickness_tip_mm, (3e5 / (900 * Om * Om * 5)) * 1e3, 1e-12, 'Centrifugal force balance on the ice layer'),
      N.check('Tip recovery temperature = T + r·V²/2cp', o.T_recovery_tip_C, -20 + (0.85 * (Om * 5) ** 2) / (2 * CP_AIR), 1e-10, 'Adiabatic wall'),
      N.check('Tip speed', o.tip_speed_ms, Om * 5, 1e-12, 'Kinematics'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.r_shed_frac < 1) out.push({ severity: 'warn', title: 'Self-shedding of blade ice', detail: `First shed after ${o.t_first_shed_min.toFixed(1)} min, outboard of ${(100 * o.r_shed_frac).toFixed(0)}% radius; a one-blade shed gives about ${(o.imbalance_force_N / 1e3).toFixed(1)} kN of rotating imbalance.`, action: 'Check hub and mount loads in Suite 10, protect structure in the shed-ice path, and use blade heating in symmetric cycles so that blades shed together.', basis: 'Centrifugal force against adhesion strength' });
    if (o.profile_power_increase_pct > 20) out.push({ severity: 'warn', title: 'Large power rise in icing', detail: `Profile power +${o.profile_power_increase_pct.toFixed(0)}% (empirical).`, action: 'Confirm the power margin in Suite 6; torque rise is the pilot’s main cue of rotor ice. More power is more fuel or battery energy, so plan the shortest exposure.', basis: 'Iced-section drag weighted by r³' });
    if (o.r_ice_free_frac < 1) out.push({ severity: 'info', title: 'Kinetic heating keeps the outer blade clear', detail: `No ice beyond ${(100 * o.r_ice_free_frac).toFixed(0)}% radius at ${i.T_C} °C.`, action: 'Heater mats can stop short of this radius at warm icing temperatures, but check the coldest design point, where ice reaches the tip.', basis: 'Recovery temperature above 0 °C' });
    return out;
  },
};

export default {
  id: 'icing', n: 13,
  tagline: 'How fast ice builds, what shape it takes, what it costs in lift, drag and power, and how much energy it takes to keep it off.',
  analyses: [impingement, accretion, degradation, protection, envelope, rotorIce],
  consumes: [{ from: 'cfd', keys: ['CLmax', 'CD0', 'k_induced'], why: 'Clean aerodynamic baseline for the ice penalties' }],
  provides: [
    { key: 'ice_thickness_mm', label: 'Ice thickness', unit: 'mm' }, { key: 'ice_mass_kg_m', label: 'Ice mass per unit span', unit: 'kg/m' }, { key: 'dCLmax_pct', label: 'Maximum-lift loss', unit: '%' },
    { key: 'dCD_pct', label: 'Drag increase', unit: '%' }, { key: 'antiice_power_W', label: 'Anti-icing power', unit: 'W' }, { key: 'freezing_fraction', label: 'Freezing fraction', unit: '-' },
  ],
  handoff: [
    { model: 'Coupled CFD–droplet–accretion on the real 3-D geometry (Navier–Stokes flow field, Eulerian droplet transport)', why: 'Needs a volume mesh and flow solution around the airframe; here the leading edge is an equivalent cylinder in potential flow', tool: 'Icing CFD suites (LEWICE / LEWICE3D, FENSAP-ICE, ONERA IGLOO)' },
    { model: 'Extended Messinger and shallow-water film models with conduction into the ice and structure', why: 'Transient film and conduction solution on a surface mesh; the classical quasi-steady Messinger balance is used instead', tool: 'FENSAP-ICE ICE3D, LEWICE thermal modules' },
    { model: 'Ice-shape feedback on the flow with mesh deformation or remeshing (multi-shot horns)', why: 'Only the leading-edge radius is updated between time steps; horn aerodynamics are not resolved', tool: 'Multi-shot icing CFD with automatic remeshing' },
    { model: 'Supercooled large droplet splashing, bounce and break-up', why: 'Empirical splashing models are calibrated on specific tunnel datasets and need a full flow field; only an SLD warning is given', tool: 'LEWICE 3.x SLD models, FENSAP-ICE DROP3D with SLD' },
    { model: 'Iced-airfoil polars and stall behaviour', why: 'Massively separated flow behind horns needs RANS or scale-resolving simulation; an empirical penalty is used here', tool: 'RANS / DDES on the predicted ice shape, icing-tunnel tests' },
    { model: 'Engine and intake icing, ice-crystal icing, windshield and probe icing', why: 'Internal-flow and mixed-phase physics are outside the leading-edge model', tool: 'Engine icing codes (GlennICE) and rig tests' },
    { model: 'Transient electro-thermal de-icing with conduction through the laminate and bond-line melting; pneumatic boot mechanics', why: 'Needs a multilayer transient thermal–mechanical model; a duty-cycle energy criterion is used', tool: 'Transient thermal FE (Suite 12 for 1-D conduction; LEWICE electro-thermal, ANSYS) and system rig tests' },
    { model: 'Rotor icing with aeroelastic coupling and cohesive ice fracture', why: 'Blade dynamics and ice fracture mechanics are not coupled here; shedding uses a local adhesion criterion', tool: 'Comprehensive rotor codes coupled to icing CFD; icing whirl-rig tests' },
  ],
};
