// Suite 16 — Flight Control Systems and Autopilot.
// Classical pitch-attitude loop with actuator, sampling and delay (Bode, Nyquist, margins, step metrics, gain
// schedule), LQR/LQG design, cascaded autopilot simulation, stability augmentation, control allocation with
// failures, and a multirotor/helicopter hover cascade.

import * as N from '../core/numerics.js';
import { isa, G0 } from '../core/atmosphere.js';

const TAU = 2 * Math.PI, fin = Number.isFinite, rad = N.rad, deg = N.deg, C = N.C;
const thin = (a, n = 300) => { const s = Math.max(1, Math.ceil(a.length / n)); return a.filter((_, k) => k % s === 0); };

// ---- control-theory helpers -------------------------------------------------------------------
/** SISO frequency response c·(jωI − A)⁻¹·b as [re, im]. */
export function freqResp(A, b, c, w) { const x = N.csolve(A.map((r, p) => r.map((v, q) => [-v, p === q ? w : 0])), b); let s = [0, 0]; for (let k = 0; k < c.length; k++) if (c[k]) s = C.add(s, C.scale(x[k], c[k])); return s; }
/** Gain and phase margins of a loop function L(ω) → [re, im], scanned on a log grid with unwrapped phase. */
export function margins(L, w1 = 1e-3, w2 = 1e3, n = 700) {
  const w = N.logspace(w1, w2, n), mag = [], ph = []; let prev = null, off = 0;
  for (const x of w) { const l = L(x); let p = C.arg(l); if (prev !== null) { while (p + off - prev > Math.PI) off -= TAU; while (p + off - prev < -Math.PI) off += TAU; } p += off; prev = p; mag.push(C.abs(l)); ph.push(p); }
  // reference the phase so that the first gain crossover is read against the nearest −180° − 360°·k line
  let pm = Infinity, wgc = NaN, gm = Infinity, wpc = NaN; const lg = mag.map(Math.log);
  for (let k = 1; k < n; k++) if ((lg[k - 1] > 0) !== (lg[k] > 0)) { const t = lg[k - 1] / (lg[k - 1] - lg[k]), p = ph[k - 1] + t * (ph[k] - ph[k - 1]), m = deg(p - (-Math.PI + TAU * Math.round((p + Math.PI) / TAU))) ; const pmk = Math.abs(m) > 180 ? 360 - Math.abs(m) : Math.abs(m), sgn = m >= 0 ? 1 : -1; if (pmk * 1 < Math.abs(pm)) { pm = sgn * pmk; wgc = Math.exp(Math.log(w[k - 1]) + t * Math.log(w[k] / w[k - 1])); } }
  for (let k = 1; k < n; k++) { const a = (ph[k - 1] + Math.PI) / TAU, b = (ph[k] + Math.PI) / TAU; if (Math.floor(a) !== Math.floor(b) || Number.isInteger(b)) { const tgt = Math.max(Math.floor(a), Math.floor(b)), t = (tgt - a) / (b - a || 1e-300), g = -(lg[k - 1] + t * (lg[k] - lg[k - 1])) * (20 / Math.LN10); if (Math.abs(g) < Math.abs(gm)) { gm = g; wpc = Math.exp(Math.log(w[k - 1]) + t * Math.log(w[k] / w[k - 1])); } } }
  return { w, mag, ph, pm_deg: pm, wgc, gm_dB: gm, wpc };
}
/** Rise time (10–90%), settling time into ±band (2% by default), overshoot and steady-state error of a step response towards ref. */
export function stepMetrics(t, y, ref, band = 0.02) {
  const n = t.length, y0 = y[0], d = ref - y0, f = (fr) => { for (let k = 0; k < n; k++) if ((y[k] - y0) / d >= fr) return t[k]; return NaN; };
  let ts = t[0]; for (let k = n - 1; k >= 0; k--) if (Math.abs(y[k] - ref) > band * Math.abs(d)) { ts = k < n - 1 ? t[k + 1] : NaN; break; }
  const pk = d > 0 ? N.amax(y) : N.amin(y);
  return { rise: f(0.9) - f(0.1), settle: ts, overshoot: Math.max(0, (100 * (pk - ref)) / d), sse: (y[n - 1] - ref) / d };
}
const eigMax = (A) => N.amax(N.eig(A).map((l) => l[0]));
const modes = (A) => N.eig(A).map((l) => ({ re: l[0], im: l[1], wn: Math.hypot(l[0], l[1]), zeta: -l[0] / (Math.hypot(l[0], l[1]) || 1) }));
/** Fastest oscillatory mode (the short period of a conventional aeroplane). */
const fastMode = (A) => modes(A).filter((m) => m.im > 1e-9).sort((p, q) => q.wn - p.wn)[0];
const leastDamped = (A, wmin = 0) => modes(A).filter((m) => m.im > 1e-9 && m.wn > wmin).sort((p, q) => p.zeta - q.zeta)[0];

// ---- plant models ---------------------------------------------------------------------------
const PLANT = [
  { key: 'vehicle', label: 'Plant type', type: 'select', options: ['aeroplane', 'rotorcraft'], default: 'aeroplane', group: 'Plant' },
  { key: 'use_upstream', label: 'Use the flight-dynamics suite model when available', type: 'bool', default: true, group: 'Plant' },
  { key: 'mass_kg', label: 'Mass', unit: 'kg', default: 70000, min: 0.01, group: 'Plant' }, { key: 'Iyy', label: 'Pitch inertia', unit: 'kg·m²', default: 3.8e6, min: 1e-6, group: 'Plant' },
  { key: 'S', label: 'Wing area', unit: 'm²', default: 122, min: 0.01, group: 'Aeroplane' }, { key: 'mac', label: 'Mean aerodynamic chord', unit: 'm', default: 4.2, min: 0.01, group: 'Aeroplane' },
  { key: 'V', label: 'True airspeed', unit: 'm/s', default: 231, min: 1, group: 'Aeroplane' }, { key: 'alt_m', label: 'Altitude', unit: 'm', default: 10668, min: -500, max: 25000, group: 'Plant' },
  { key: 'CLa', label: 'Lift-curve slope CLα', unit: '1/rad', default: 6, min: 1, max: 9, group: 'Aeroplane' }, { key: 'Cma', label: 'Pitch stiffness Cmα', unit: '1/rad', default: -1.5, min: -10, max: 3, group: 'Aeroplane' },
  { key: 'Cmq', label: 'Pitch damping Cmq', unit: '1/rad', default: -40, min: -200, max: 0, group: 'Aeroplane' }, { key: 'CLde', label: 'Elevator lift CLδe', unit: '1/rad', default: 0.45, min: 0, max: 3, group: 'Aeroplane' },
  { key: 'Cmde', label: 'Elevator power Cmδe', unit: '1/rad', default: -2.5, min: -10, max: -0.01, group: 'Aeroplane' }, { key: 'CD0', label: 'Zero-lift drag coefficient', unit: '-', default: 0.02, min: 1e-6, group: 'Aeroplane' },
  { key: 'k_ind', label: 'Induced drag factor', unit: '-', default: 0.042, min: 0, group: 'Aeroplane' }, { key: 'T_max', label: 'Maximum thrust at this condition', unit: 'N', default: 90000, min: 0, group: 'Aeroplane' },
  { key: 'Xu', label: 'Hover speed damping Xu', unit: '1/s', default: -0.02, min: -5, max: 0, group: 'Rotorcraft' }, { key: 'Mu', label: 'Hover speed stability Mu', unit: 'rad/(s·m)', default: 0.016, min: -2, max: 5, group: 'Rotorcraft' },
  { key: 'Mq', label: 'Hover pitch damping Mq', unit: '1/s', default: -0.6, min: -50, max: 0, group: 'Rotorcraft' }, { key: 'Mc', label: 'Pitch control power', unit: 'rad/s² per unit control', default: 8.5, min: 0.01, group: 'Rotorcraft' },
  { key: 'Zw', label: 'Heave damping Zw', unit: '1/s', default: -0.25, min: -10, max: 0, group: 'Rotorcraft' }, { key: 'Zc', label: 'Heave control power', unit: 'm/s² per unit control', default: -9.81, min: -500, max: -0.01, group: 'Rotorcraft' },
];
const slope = (A, sw, be2 = 1) => (A > 0 ? (TAU * A) / (2 + Math.sqrt(4 + A * A * be2 * (1 + Math.tan(sw) ** 2 / be2))) : TAU);
const tauFlap = (E) => { const th = Math.acos(N.clamp(2 * E - 1, -1, 1)); return 0.85 * (1 - (th - Math.sin(th)) / Math.PI); };
/** Compact preliminary derivative estimate from the shared case (the detailed build-up lives in Suite 4). */
function estimate(c, d) {
  const at = isa(c.atm.alt_m, c.atm.dISA_K), V = c.flight.V_ms, M = V / at.a, be2 = Math.max(1 - M * M, 0.19), sw = rad(c.wing.sweep_deg), S = c.wing.S_m2, b = c.wing.b_m, mac = d.mac, lam = c.wing.taper;
  const aw = slope(d.AR, sw, be2), ah = slope(c.htail.b_m ** 2 / Math.max(c.htail.S_m2, 1e-9), sw, be2), av = slope((1.55 * c.vtail.b_m ** 2) / Math.max(c.vtail.S_m2, 1e-9), sw, be2), eps = Math.min(0.7, (2 * aw) / (Math.PI * d.AR)), VH = (c.htail.S_m2 * c.htail.arm_m) / (S * mac);
  const vol = 0.7 * (Math.PI / 4) * c.fuselage.dia_m ** 2 * c.fuselage.len_m, CLa = aw + 0.9 * ah * (c.htail.S_m2 / S) * (1 - eps), Cma = aw * (c.mass.cg_pct_mac / 100 - 0.25) + (1.8 * vol) / (S * mac) - 0.9 * ah * VH * (1 - eps), te = tauFlap(c.controls.Se_Sh), tr = tauFlap(c.controls.Sr_Sv);
  const q = 0.5 * at.rho * V * V, CL = d.W / (q * S), zv = 0.5 * c.vtail.b_m + 0.25 * c.fuselage.dia_m, CYv = -av * (c.vtail.S_m2 / S), lv = c.vtail.arm_m, yb = (1 + 2 * lam) / (3 * (1 + lam));
  const isP = c.prop.type !== 'turbofan' && c.prop.type !== 'turbojet', Tm = isP && d.P_total > 0 ? Math.min(d.T_total > 0 ? d.T_total * at.sigma ** 0.7 : Infinity, (c.prop.eta_prop * d.P_total * (c.prop.type === 'electric' ? 1 : at.sigma ** 0.7)) / V) : d.T_total * at.sigma ** 0.7;
  return { at, V, q, CL, CLa, Cma, Cmq: -2.2 * 0.9 * ah * VH * (c.htail.arm_m / mac), CLde: 0.9 * ah * (c.htail.S_m2 / S) * te, Cmde: -0.9 * ah * VH * te, T_max: Tm,
    CYb: CYv - (2 * (Math.PI / 4) * c.fuselage.dia_m ** 2) / S, Cnb: av * (c.vtail.S_m2 * lv) / (S * b) - (1.8 * vol) / (S * b), Clb: 0.75 * (-(aw * rad(c.wing.dihedral_deg) * yb) / 2 - CL * Math.tan(sw) * yb) + CYv * (zv / b),
    Clp: (-aw * (1 + 3 * lam)) / (12 * (1 + lam)), Cnp: -CL / 8, Clr: CL / 4 - 2 * (lv / b) * (zv / b) * CYv, Cnr: 2 * CYv * (lv / b) ** 2 - c.aero.CD0 / 4, CYr: -2 * CYv * (lv / b),
    CYdr: av * tr * (c.vtail.S_m2 / S), Cndr: (-av * tr * c.vtail.S_m2 * lv) / (S * b), Cldr: av * tr * (c.vtail.S_m2 / S) * (zv / b) };
}
/** Compact hover derivative estimate (see Suite 4 for the model). */
function hoverEst(c) {
  const heli = c.meta.type === 'helicopter', r = c.rotor, m = c.mass.mtow_kg, n = heli ? 1 : Math.max(1, c.prop.n_eng), rho = isa(c.atm.alt_m).rho, T = (m * G0) / n, A = Math.PI * r.R_m ** 2, Om = ((r.rpm || 300) * TAU) / 60, vt = Om * r.R_m;
  const CT = T / (rho * A * vt * vt), lam = Math.sqrt(CT / 2), sig = ((r.n_blades || 2) * (r.chord_m || 0.1 * r.R_m)) / (Math.PI * r.R_m), a = r.cla, dTdw = (2 * a * sig * rho * A * vt * lam) / (16 * lam + a * sig), L = Math.max(c.fuselage.len_m, 4 * r.R_m, 0.3);
  const Iyy = c.mass.Iyy > 0 ? c.mass.Iyy : 0.02 * m * L * L, Ixx = c.mass.Ixx > 0 ? c.mass.Ixx : 0.02 * m * L * L, Izz = c.mass.Izz > 0 ? c.mass.Izz : 0.035 * m * L * L, Q = (1.15 * T * lam * vt + (rho * A * vt ** 3 * sig * 0.01) / 8) / Om;
  if (heli) { const h = 0.75 * c.fuselage.dia_m, da = ((16 * CT) / (sig * a) + 2 * lam) / vt, Kh = T * h + (r.n_blades / 2) * r.hinge_offset * r.R_m * Om * Om * ((r.blade_mass_kg * r.R_m) / 2); return { heli, n, m, Iyy, Ixx, Izz, Xu: -G0 * da, Mu: (Kh * da) / Iyy, Mq: (-Kh * 16) / (r.lock * Om) / Iyy, Mc: Kh / Iyy, Zw: -dTdw / m, Zc: -((rho * A * vt * vt * sig * a) / 6 / (1 + (a * sig) / (16 * lam))) / m, arm: r.tr_arm_m || 1.2 * r.R_m, T, Q, tau_m: 0.1 }; }
  const arm = c.wing.b_m > 0 ? 0.35 * c.wing.b_m : 1.6 * r.R_m, h = 0.15 * r.R_m, kd = 0.25;
  return { heli, n, m, Iyy, Ixx, Izz, Xu: -kd, Mu: (m * kd * h) / Iyy, Mq: (-dTdw * n * arm * arm) / 2 / Iyy, Mc: (m * G0 * arm) / (Math.SQRT2 * Iyy), Zw: (-n * dTdw) / m, Zc: -G0, arm, T, Q, tau_m: N.clamp(0.25 * r.R_m, 0.03, 0.3) };
}
const plantDefaults = (c, up, d) => {
  if (c.wing.S_m2 > 0) { const e = estimate(c, d), f = up.flightdyn || {}; return { vehicle: 'aeroplane', mass_kg: c.mass.mtow_kg, Iyy: c.mass.Iyy > 0 ? c.mass.Iyy : c.mass.mtow_kg * (0.18 * c.fuselage.len_m) ** 2, S: c.wing.S_m2, mac: d.mac, V: c.flight.V_ms, alt_m: c.atm.alt_m, CLa: f.CLa_total ?? e.CLa, Cma: f.Cm_alpha ?? e.Cma, Cmq: f.Cmq ?? e.Cmq, CLde: e.CLde, Cmde: f.Cmde ?? e.Cmde, CD0: up.cfd?.CD0 ?? c.aero.CD0, k_ind: d.k_induced || undefined, T_max: e.T_max }; }
  const h = hoverEst(c); return { vehicle: 'rotorcraft', mass_kg: h.m, Iyy: h.Iyy, alt_m: c.atm.alt_m, Xu: h.Xu, Mu: h.Mu, Mq: h.Mq, Mc: h.Mc, Zw: h.Zw, Zc: h.Zc };
};
/** Longitudinal plant x = [u, w, q, θ]; B columns: [pitch control, thrust/collective]. */
function plant(i, ctx, V = i.V) {
  const f = ctx?.up?.flightdyn, ok = (A, B) => Array.isArray(A) && A.length === 4 && Array.isArray(B) && B.length === 4 && B[0].length >= 1;
  if (i.vehicle === 'rotorcraft') {
    if (i.use_upstream && f && ok(f.A_hover_long, f.B_hover_long)) return { A: f.A_hover_long, B: f.B_hover_long, u0: 0, src: 'Suite 4 hover model' };
    return { A: [[i.Xu, 0, 0, -G0], [0, i.Zw, 0, 0], [i.Mu, 0, i.Mq, 0], [0, 0, 1, 0]], B: [[0, 0], [0, i.Zc], [i.Mc, 0], [0, 0]], u0: 0, src: 'built-in hover model' };
  }
  if (i.use_upstream && f && ok(f.A_long, f.B_long) && V === i.V && Math.abs((f.u0_ms ?? V) - V) < 1e-6) return { A: f.A_long, B: f.B_long, u0: V, src: 'Suite 4 linear model' };
  const at = isa(i.alt_m), q = 0.5 * at.rho * V * V, QS = q * i.S, m = i.mass_kg, CL = (m * G0) / QS, CD = i.CD0 + i.k_ind * CL * CL, c = i.mac;
  const Xu = (-2 * CD * QS) / (m * V), Xw = (-(2 * i.k_ind * CL * i.CLa - CL) * QS) / (m * V), Zu = (-2 * CL * QS) / (m * V), Zw = (-(i.CLa + CD) * QS) / (m * V), Mw = (i.Cma * QS * c) / (V * i.Iyy), Mq = (i.Cmq * (c / (2 * V)) * QS * c) / i.Iyy;
  return { A: [[Xu, Xw, 0, -G0], [Zu, Zw, V, 0], [0, Mw, Mq, 0], [0, 0, 1, 0]], B: [[0, i.T_max / m], [(-i.CLde * QS) / m, 0], [(i.Cmde * QS * c) / i.Iyy, 0], [0, 0]], u0: V, src: 'built-in model from the listed derivatives' };
}
/** PID gains that put the gain crossover of C·G at wc with phase margin pm (Ki = Kp·wc/10). G is the loop without the controller. */
function tunePID(G, wc, pmDeg) { const g = G(wc), D = C.div(C.exp([0, -Math.PI + rad(pmDeg)]), g), Kp = D[0], Kd = (D[1] + Kp / 10) / wc; return { Kp, Ki: (Kp * wc) / 10, Kd }; }

// ---- 1. classical pitch loop ------------------------------------------------------------------
const LOOP_IN = [
  { key: 'autotune', label: 'Tune the PID for the target crossover and phase margin', type: 'bool', default: true, group: 'Controller' },
  { key: 'wc', label: 'Target crossover frequency', unit: 'rad/s', default: 2, min: 0.01, max: 200, group: 'Controller', help: 'About 1–2 times the short-period frequency, and well below the actuator bandwidth' },
  { key: 'pm_target', label: 'Target phase margin', unit: 'deg', default: 55, min: 20, max: 80, group: 'Controller' },
  { key: 'Kp', label: 'Proportional gain (manual)', unit: 'rad/rad', default: -1, group: 'Controller', help: 'Negative for a conventional elevator (trailing edge down pitches nose down)' },
  { key: 'Ki', label: 'Integral gain (manual)', unit: '1/s', default: -0.2, group: 'Controller' }, { key: 'Kd', label: 'Pitch-rate gain (manual)', unit: 's', default: -0.5, group: 'Controller' },
  { key: 'w_act', label: 'Actuator natural frequency', unit: 'rad/s', default: 30, min: 1, max: 1000, group: 'Actuator', help: '20–40 rad/s hydraulic servo, 40–100 electromechanical, motor/ESC lag for multirotors' },
  { key: 'z_act', label: 'Actuator damping ratio', unit: '-', default: 0.7, min: 0.2, max: 2, group: 'Actuator' },
  { key: 'rate_lim', label: 'Actuator rate limit', unit: 'deg/s (or %/s)', default: 60, min: 0.1, max: 5000, group: 'Actuator' }, { key: 'pos_lim', label: 'Actuator travel limit', unit: 'deg (or %)', default: 25, min: 0.1, max: 100, group: 'Actuator' },
  { key: 'Ts', label: 'Controller sample time', unit: 's', default: 0.02, min: 1e-4, max: 0.5, group: 'Digital implementation' }, { key: 'delay', label: 'Computation and transport delay', unit: 's', default: 0.02, min: 0, max: 0.5, group: 'Digital implementation' },
  { key: 'step_deg', label: 'Pitch-attitude step command', unit: 'deg', default: 2, min: 0.01, max: 30, group: 'Test input' },
  { key: 'gm_req', label: 'Required gain margin', unit: 'dB', default: 6, min: 1, max: 20, group: 'Limits', help: '6 dB and 45° are the customary minimum flight-control margins' }, { key: 'pm_req', label: 'Required phase margin', unit: 'deg', default: 45, min: 10, max: 80, group: 'Limits' },
];
const loopDefaults = (c, up, d) => {
  const p = plantDefaults(c, up, d), rot = p.vehicle === 'rotorcraft', h = rot ? hoverEst(c) : null, i = { ...Object.fromEntries(PLANT.map((f) => [f.key, f.default])), ...p }, sp = rot ? null : fastMode(plant(i, null).A);
  const wa = rot ? (h.heli ? 25 : 1 / h.tau_m) : 30, wc = rot ? Math.min(0.25 * wa, Math.max(2.5, 3 * Math.cbrt(G0 * Math.abs(p.Mu)))) : N.clamp(1.3 * (sp?.wn ?? 2), 0.5, 0.2 * wa);
  return { ...p, w_act: wa, wc, rate_lim: rot ? (h.heli ? 100 : 2000) : c.controls.rate_max_dps, pos_lim: rot ? (h.heli ? 20 : 100) : c.controls.de_max_deg, Ts: rot && !h.heli ? 0.004 : 0.02, delay: rot && !h.heli ? 0.004 : 0.02 };
};
/** Loop without the controller: θ/δ plant × second-order actuator × pure delay (Ts/2 for the sample-and-hold plus computation delay). */
const openLoop = (P, i) => (w) => { const g = freqResp(P.A, P.B.map((r) => r[0]), [0, 0, 0, 1], w), act = C.div([i.w_act ** 2, 0], [i.w_act ** 2 - w * w, 2 * i.z_act * i.w_act * w]); return C.mul(C.mul(g, act), C.exp([0, -w * (i.Ts / 2 + i.delay)])); };
const pidResp = (g, w) => [g.Kp, g.Kd * w - g.Ki / w];
/** Discrete-time closed-loop simulation with rate- and travel-limited actuator. Controls are in the plant's control unit (rad or fraction). */
function loopSim(P, i, g, ref, tEnd, limScale) {
  const b = P.B.map((r) => r[0]), dt = Math.min(i.Ts / 2, 0.2 / i.w_act), nSub = Math.max(1, Math.round(i.Ts / dt)), h = i.Ts / nSub, n = Math.ceil(tEnd / i.Ts), nd = Math.round(i.delay / i.Ts), buf = [];
  const rl = i.rate_lim * limScale, pl = i.pos_lim * limScale; let x = [0, 0, 0, 0], da = 0, dv = 0, integ = 0, sat = 0, wasSat = false, tSat = 0, umax = 0; const t = [], th = [], de = [];
  for (let k = 0; k <= n; k++) {
    buf.push([x[3], x[2]]); const [thm, qm] = buf[Math.max(0, buf.length - 1 - nd)], e = ref - thm; let u = g.Kp * e + g.Ki * integ - g.Kd * qm; const us = N.clamp(u, -pl, pl);
    if (us === u) integ += e * i.Ts; // conditional integration as anti-windup
    u = us; t.push(k * i.Ts); th.push(x[3]); de.push(da); umax = Math.max(umax, Math.abs(da)); let s = false;
    for (let j = 0; j < nSub; j++) {
      let acc = i.w_act ** 2 * (u - da) - 2 * i.z_act * i.w_act * dv; dv += acc * h; if (Math.abs(dv) > rl) { dv = N.clamp(dv, -rl, rl); s = true; } da += dv * h; if (Math.abs(da) > pl) { da = N.clamp(da, -pl, pl); dv = 0; s = true; }
      const f = (y) => N.vadd(N.matvec(P.A, y), b, da), k1 = f(x), k2 = f(N.vadd(x, k1, h / 2)), k3 = f(N.vadd(x, k2, h / 2)), k4 = f(N.vadd(x, k3, h)); x = x.map((v, q) => v + (h / 6) * (k1[q] + 2 * k2[q] + 2 * k3[q] + k4[q]));
    }
    if (s) { tSat += i.Ts; if (!wasSat) sat++; } wasSat = s; if (!fin(x[3]) || Math.abs(x[3]) > 1e3) break;
  }
  return { t, th, de, sat, tSat, umax };
}
function loopDesign(P, i) {
  const G = openLoop(P, i), g = i.autotune ? tunePID(G, i.wc, i.pm_target) : { Kp: i.Kp, Ki: i.Ki, Kd: i.Kd }, L = (w) => C.mul(pidResp(g, w), G(w)), lo = Math.min(0.01, i.wc / 300), hi = Math.max(10 * i.w_act, 40 * i.wc);
  return { G, g, L, m: margins(L, lo, hi, 700) };
}
const pitchLoop = {
  id: 'pitchloop', title: 'Pitch-attitude loop: PID with actuator, sampling and delay', fidelity: 'reduced-order',
  summary: 'Classical design of the pitch-attitude hold loop: loop shaping of a PID with pitch-rate feedback, Bode and Nyquist plots, gain and phase margins, bandwidth, a sampled-data step response with actuator limits, and a gain schedule over speed.',
  equations: ['Transfer function equations', 'Proportional–integral–derivative control equations', 'Frequency response equations', 'Actuator dynamic equations', 'State-space equations', 'Root locus relations'],
  inputs: [...LOOP_IN, ...PLANT], defaults: loopDefaults,
  run(i, ctx) {
    const P = plant(i, ctx), rot = i.vehicle === 'rotorcraft', ls = rot ? 0.01 : rad(1), { G, g, L, m } = loopDesign(P, i), warnings = [];
    const ref = rad(i.step_deg), wEst = fin(m.wgc) ? m.wgc : i.wc, sim = loopSim(P, i, g, ref, N.clamp(150 / wEst, 15, 120), ls), sm = stepMetrics(sim.t, sim.th, ref, 0.05), diverged = !fin(sim.th[sim.th.length - 1]) || Math.abs(sim.th[sim.th.length - 1]) > 50 * ref;
    // closed-loop bandwidth: |T| = |L/(1+L)| falls through −3 dB
    const Tm = m.w.map((w) => { const l = L(w); return C.abs(C.div(l, C.add([1, 0], l))); }); let bw = NaN; for (let k = 1; k < Tm.length; k++) if (Tm[k - 1] >= Math.SQRT1_2 && Tm[k] < Math.SQRT1_2) bw = m.w[k - 1] + ((Tm[k - 1] - Math.SQRT1_2) / (Tm[k - 1] - Tm[k])) * (m.w[k] - m.w[k - 1]);
    const Ms = N.amax(m.w.map((w) => 1 / C.abs(C.add([1, 0], L(w))))), phLoss = deg(wEst * (i.Ts / 2 + i.delay)), dmaxFree = (i.rate_lim * ls) / wEst, nyq = m.w.map(L).filter((l) => C.abs(l) < 6);
    // root locus on loop gain, with a first-order Padé delay: closed-loop poles for gain factors 0.1…4
    const b = P.B.map((r) => r[0]), td = i.Ts / 2 + i.delay, cl = (kf) => { // states: plant(4), actuator(2), integrator, Padé
      const A = N.zeros(8), K = [g.Kp * kf, g.Ki * kf, g.Kd * kf]; for (let p = 0; p < 4; p++) { for (let q = 0; q < 4; q++) A[p][q] = P.A[p][q]; A[p][4] = b[p]; }
      A[4][5] = 1; A[5][4] = -(i.w_act ** 2); A[5][5] = -2 * i.z_act * i.w_act; A[6][3] = -1; // ė_int = −θ (zero reference)
      // u = −Kp·θ + Ki·z − Kd·q through Padé: v̇ = (2/td)(u − v)·… output y = 2v − u
      const u = [0, 0, -K[2], -K[0], 0, 0, K[1], 0], a = td > 1e-6 ? 2 / td : 1e6; for (let q = 0; q < 8; q++) { A[7][q] += a * u[q]; A[5][q] += i.w_act ** 2 * -u[q]; } A[7][7] -= a; A[5][7] += i.w_act ** 2 * 2; return N.eig(A); };
    const kf = [0.1, 0.2, 0.35, 0.5, 0.7, 1, 1.4, 2, 2.8, 4], loc = { x: [], y: [] }; for (const k of kf) for (const l of cl(k)) if (l[0] > -3 * i.w_act) { loc.x.push(l[0]); loc.y.push(l[1]); } const nom = cl(1), stable = N.amax(nom.map((l) => l[0])) < 0;
    // gain schedule across speed (aeroplane, built-in model)
    const sched = rot ? [] : [0.6, 0.8, 1, 1.2, 1.4].map((f) => { try { const Pv = plant({ ...i, use_upstream: false }, null, f * i.V), d = loopDesign(Pv, i); return [f * i.V, d.g.Kp, d.g.Ki, d.g.Kd, d.m.gm_dB, d.m.pm_deg]; } catch { return null; } }).filter(Boolean);
    if (!stable || diverged) warnings.push('The closed loop is unstable with these gains.');
    if (i.wc > 0.3 * i.w_act) warnings.push('The crossover is above 30% of the actuator bandwidth: actuator phase lag dominates the margins.');
    if (phLoss > 15) warnings.push(`Sampling and delay cost ${phLoss.toFixed(0)}° of phase at crossover; raise the sample rate or reduce latency.`);
    if (stable && !diverged && !fin(sm.settle)) warnings.push('The attitude did not settle inside the 5% band within the simulated time: a slow speed mode keeps drifting under attitude hold, so an outer speed or position loop is needed.');
    if (sim.sat) warnings.push(`The actuator reached a rate or travel limit ${sim.sat} time(s) during the step (${sim.tSat.toFixed(2)} s in total).`);
    if (rot && N.amax(modes(P.A).map((x) => x.re)) > 0) warnings.push('The plant is open-loop unstable: the gain margin has a lower bound as well (reducing the gain destabilises the loop).');
    const gmOk = Math.abs(m.gm_dB) >= i.gm_req, pmOk = m.pm_deg >= i.pm_req;
    return {
      kpis: [
        { key: 'gm_dB', label: 'Gain margin', value: Math.abs(m.gm_dB), unit: 'dB', status: stable && gmOk ? 'ok' : 'bad', note: `Requirement ${i.gm_req} dB; nearest −180° crossing at ${fin(m.wpc) ? m.wpc.toPrecision(3) : '—'} rad/s${m.gm_dB < 0 ? ' (gain-reduction margin)' : ''}` },
        { key: 'pm_deg', label: 'Phase margin', value: m.pm_deg, unit: 'deg', status: stable && pmOk ? 'ok' : 'bad', note: `Requirement ${i.pm_req}°` },
        { key: 'crossover_rads', label: 'Gain crossover frequency', value: m.wgc, unit: 'rad/s' }, { key: 'bandwidth_rads', label: 'Closed-loop bandwidth (−3 dB)', value: bw, unit: 'rad/s' },
        { key: 'delay_margin_s', label: 'Delay margin', value: fin(m.wgc) ? rad(m.pm_deg) / m.wgc : NaN, unit: 's' }, { key: 'Ms', label: 'Peak sensitivity', value: Ms, unit: '-', status: Ms < 2 ? 'ok' : 'warn', note: 'Below about 2 for robust loops' },
        { key: 'rise_s', label: 'Rise time (10–90%)', value: sm.rise, unit: 's' }, { key: 'settling_s', label: 'Settling time (5% band)', value: diverged ? NaN : sm.settle, unit: 's' }, { key: 'overshoot_pct', label: 'Overshoot', value: diverged ? NaN : sm.overshoot, unit: '%', status: sm.overshoot < 25 ? 'ok' : 'warn' },
        { key: 'Kp', label: 'Proportional gain', value: g.Kp, unit: '-' }, { key: 'Ki', label: 'Integral gain', value: g.Ki, unit: '1/s' }, { key: 'Kd', label: 'Pitch-rate gain', value: g.Kd, unit: 's' },
        { key: 'phase_loss_deg', label: 'Phase lost to sampling and delay', value: phLoss, unit: 'deg' }, { key: 'saturation_events', label: 'Actuator saturation events', value: sim.sat, unit: '-', status: sim.sat ? 'warn' : 'ok' },
        { key: 'peak_control', label: 'Peak actuator deflection', value: sim.umax / ls, unit: rot ? '%' : 'deg' },
        { key: 'pio_amp', label: 'Largest sinusoidal command at crossover before rate limiting', value: dmaxFree / ls, unit: rot ? '%' : 'deg', note: 'Heuristic rate-limit onset indicator: rate limit / crossover frequency' },
      ].filter((k) => fin(k.value) || ['gm_dB', 'pm_deg', 'settling_s', 'overshoot_pct', 'bandwidth_rads'].includes(k.key)),
      plots: [
        { type: 'line', title: 'Bode magnitude of the loop', xlabel: 'Frequency [rad/s]', ylabel: 'Magnitude [dB]', xlog: true, series: [{ name: '|L|', x: thin(m.w), y: thin(m.mag.map((v) => 20 * Math.log10(Math.max(v, 1e-12)))) }, { name: 'Closed loop |T|', x: thin(m.w), y: thin(Tm.map((v) => 20 * Math.log10(Math.max(v, 1e-12)))), style: 'dash' }], annotations: [{ y: 0, label: '0 dB' }] },
        { type: 'line', title: 'Bode phase of the loop', xlabel: 'Frequency [rad/s]', ylabel: 'Phase [deg]', xlog: true, series: [{ name: 'arg L', x: thin(m.w), y: thin(m.ph.map(deg)) }], annotations: [{ y: -180, label: '−180°' }] },
        { type: 'line', title: 'Nyquist diagram', xlabel: 'Real', ylabel: 'Imaginary', equalAspect: true, series: [{ name: 'L(jω)', x: thin(nyq.map((l) => l[0])), y: thin(nyq.map((l) => l[1])) }, { name: 'Critical point', x: [-1], y: [0], style: 'points' }] },
        { type: 'line', title: 'Pitch-attitude step response', xlabel: 'Time [s]', ylabel: 'Pitch attitude [deg]', series: [{ name: 'θ', x: thin(sim.t), y: thin(sim.th.map((v) => N.clamp(deg(v), -1e3, 1e3))) }], annotations: [{ y: i.step_deg, label: 'Command' }] },
        { type: 'line', title: 'Actuator deflection', xlabel: 'Time [s]', ylabel: rot ? 'Control [%]' : 'Deflection [deg]', series: [{ name: 'Actuator', x: thin(sim.t), y: thin(sim.de.map((v) => v / ls)) }] },
        { type: 'line', title: 'Root locus on loop gain (0.1× to 4×)', xlabel: 'Real part [1/s]', ylabel: 'Imaginary part [rad/s]', series: [{ name: 'Closed-loop poles', x: loc.x, y: loc.y, style: 'points' }, { name: 'Nominal gain', x: nom.filter((l) => l[0] > -3 * i.w_act).map((l) => l[0]), y: nom.filter((l) => l[0] > -3 * i.w_act).map((l) => l[1]), style: 'points' }], annotations: [{ x: 0, label: 'Stability boundary' }] },
      ],
      tables: sched.length ? [{ title: 'Gain schedule against airspeed (built-in model)', columns: ['TAS [m/s]', 'Kp', 'Ki [1/s]', 'Kd [s]', 'Gain margin [dB]', 'Phase margin [deg]'], rows: sched }] : [],
      outputs: { Kp_theta: g.Kp, Ki_theta: g.Ki, Kd_theta: g.Kd },
      warnings, models: [`Plant: ${P.src}`, 'PID controller with pitch-rate (derivative-on-measurement) feedback', 'Second-order actuator with rate and travel limits', 'Sample-and-hold and delay as a pure time delay in the frequency domain, exact discrete update in the time simulation', 'Gain-scheduled controller table'],
      assumptions: ['Linear plant about the trim condition', 'Ideal attitude and rate sensors (sensor dynamics and noise are in Suite 17)', 'Root locus uses a first-order Padé approximation of the delay', 'Actuator bandwidth and damping, sample time, delay and the 6 dB / 45° margin requirement are customary values, not data for a specific system'],
    };
  },
  calibration: { params: [{ key: 'w_act', min: 2, max: 300 }, { key: 'delay', min: 0, max: 0.3 }, { key: 'Cmde', min: -8, max: -0.05 }], sweep: 'wc', target: 'pm_deg', note: 'Measured loop frequency response (frequency sweeps on the iron bird or in flight) to identify actuator bandwidth, latency and control power' },
  verify() {
    const m1 = margins((w) => C.div([1, 0], C.mul([0, w], C.mul([1, w], [2, w]))), 1e-2, 1e2, 2000), m2 = margins((w) => C.div([1, 0], C.mul([0, w], [1, w])), 1e-2, 1e2, 2000);
    // standard second-order step response: overshoot exp(−πζ/√(1−ζ²)), 2% settling ≈ time when the envelope reaches 2%
    const z = 0.4, wn = 3, t = N.linspace(0, 12, 6001), wd = wn * Math.sqrt(1 - z * z), y = t.map((x) => 1 - Math.exp(-z * wn * x) * (Math.cos(wd * x) + ((z * wn) / wd) * Math.sin(wd * x))), sm = stepMetrics(t, y, 1);
    const b = Object.fromEntries(pitchLoop.inputs.map((f) => [f.key, f.default])), o = N.kv(pitchLoop.run({ ...b, Ts: 1e-3, delay: 0 }));
    const g = freqResp([[0, 1], [-4, -0.4]], [0, 1], [1, 0], 2);
    return [N.check('Gain margin of 1/(s(s+1)(s+2)) is 6', 10 ** (m1.gm_dB / 20), 6, 2e-3, 'Routh: critical gain 6'), N.check('Phase margin of 1/(s(s+1)) is 51.83°', m2.pm_deg, 51.827, 1e-3, 'Exact crossover ω² = (√5 − 1)/2'),
      N.check('Second-order overshoot exp(−πζ/√(1−ζ²))', sm.overshoot, 100 * Math.exp((-Math.PI * z) / Math.sqrt(1 - z * z)), 2e-3, 'Standard second-order system'),
      N.check('State-space frequency response at resonance', C.abs(g), 1 / 0.8, 1e-9, '1/(ωn² − ω² + 2ζωn·jω) at ω = ωn'),
      N.check('Tuned loop achieves the requested phase margin', o.pm_deg, 55, 0.02, 'Loop-shaping condition L(jωc) = exp(j(−180° + PM))'), N.check('…at the requested crossover', o.crossover_rads, 2, 0.01, 'Loop-shaping condition')];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.gm_dB < i.gm_req || o.pm_deg < i.pm_req) out.push({ severity: 'critical', title: 'Stability margins below the requirement', detail: `Gain margin ${o.gm_dB.toFixed(1)} dB, phase margin ${o.pm_deg.toFixed(0)}°.`, action: 'Lower the crossover frequency, add lead (more pitch-rate gain), reduce latency or use a faster actuator.', basis: `${i.gm_req} dB / ${i.pm_req}° margin requirement` });
    if (o.saturation_events > 0) out.push({ severity: 'warn', title: 'Actuator saturation during the manoeuvre', detail: `${o.saturation_events} event(s); peak deflection ${o.peak_control.toFixed(1)}.`, action: 'Rate saturation adds phase lag and can trigger pilot-induced oscillation: reduce loop gain or command rate, or raise the actuator rate capability.', basis: 'Rate-limit onset (heuristic)' });
    if (o.phase_loss_deg > 15) out.push({ severity: 'advise', title: 'Digital implementation is eating the phase margin', detail: `${o.phase_loss_deg.toFixed(0)}° lost at crossover.`, action: 'Sample at least 20–30 times the crossover frequency and minimise sensor-to-actuator latency.', basis: 'Phase lag ω·(Ts/2 + delay)' });
    return out;
  },
};

// ---- 2. LQR / LQG ---------------------------------------------------------------------------
const lqrAn = {
  id: 'lqr', title: 'Optimal control: LQR state feedback and LQG with a Kalman filter', fidelity: 'reduced-order',
  summary: 'Linear-quadratic regulator on the longitudinal model with Bryson weighting, a steady-state Kalman filter from pitch-rate and attitude measurements, closed-loop poles and the disturbance recovery compared with the open loop.',
  equations: ['Riccati equations', 'Linear quadratic regulator equations', 'Kalman filter equations', 'State-space equations', 'Controllability and observability equations', 'Lyapunov stability equations'],
  inputs: [
    { key: 'u_max', label: 'Acceptable speed deviation', unit: 'm/s', default: 5, min: 0.01, group: 'Weights (Bryson)' }, { key: 'w_max', label: 'Acceptable vertical-speed deviation', unit: 'm/s', default: 3, min: 0.01, group: 'Weights (Bryson)' },
    { key: 'q_max_dps', label: 'Acceptable pitch rate', unit: 'deg/s', default: 5, min: 0.01, group: 'Weights (Bryson)' }, { key: 'th_max_deg', label: 'Acceptable pitch attitude deviation', unit: 'deg', default: 2, min: 0.01, group: 'Weights (Bryson)' },
    { key: 'c1_max', label: 'Acceptable pitch-control deflection', unit: 'deg (or %)', default: 5, min: 0.01, group: 'Weights (Bryson)' }, { key: 'c2_max', label: 'Acceptable thrust / collective change', unit: '%', default: 20, min: 0.01, group: 'Weights (Bryson)' },
    { key: 'sigma_gust', label: 'Process noise: vertical gust RMS', unit: 'm/s', default: 1.5, min: 0.01, group: 'Kalman filter' }, { key: 'sigma_q_dps', label: 'Pitch-rate sensor noise', unit: 'deg/s', default: 0.1, min: 1e-4, group: 'Kalman filter' },
    { key: 'sigma_th_deg', label: 'Attitude sensor noise', unit: 'deg', default: 0.2, min: 1e-4, group: 'Kalman filter' }, { key: 'th0_deg', label: 'Initial pitch disturbance', unit: 'deg', default: 3, min: 0.01, max: 20, group: 'Test input' },
    ...PLANT],
  defaults: (c, up, d) => { const p = plantDefaults(c, up, d); return p.vehicle === 'rotorcraft' ? { ...p, u_max: 1, w_max: 0.5, c1_max: 3 } : p; },
  run(i, ctx) {
    const P = plant(i, ctx), rot = i.vehicle === 'rotorcraft', cs = rot ? 0.01 : rad(1), warnings = [];
    const Q = N.zeros(4), R = [[1 / (i.c1_max * cs) ** 2, 0], [0, 1 / (i.c2_max * 0.01) ** 2]]; [i.u_max, i.w_max, rad(i.q_max_dps), rad(i.th_max_deg)].forEach((v, k) => (Q[k][k] = 1 / (v * v)));
    const { K, P: Pr } = N.lqr(P.A, P.B, Q, R), Acl = N.madd(P.A, N.matmul(P.B, K), -1);
    // Kalman filter: measurements q and θ; process noise enters like a vertical gust (second column of A)
    const Cm = [[0, 0, 1, 0], [0, 0, 0, 1]], Gw = P.A.map((r) => [-r[1]]), Wn = N.mscale(N.matmul(Gw, N.transpose(Gw)), i.sigma_gust ** 2), Wr = Wn.map((r, p) => r.map((v, q) => v + (p === q ? 1e-8 : 0))), Vn = [[rad(i.sigma_q_dps) ** 2, 0], [0, rad(i.sigma_th_deg) ** 2]];
    const Pe = N.care(N.transpose(P.A), N.transpose(Cm), Wr, Vn), Lk = N.matmul(N.matmul(Pe, N.transpose(Cm)), N.inv(Vn)), Aest = N.madd(P.A, N.matmul(Lk, Cm), -1);
    // residual of the Riccati equation as a self-check
    const At = N.transpose(P.A), Ric = N.madd(N.madd(N.madd(N.matmul(At, Pr), N.matmul(Pr, P.A)), N.matmul(N.matmul(N.matmul(Pr, P.B), N.inv(R)), N.matmul(N.transpose(P.B), Pr)), -1), Q), ricErr = N.amax(Ric.flat().map(Math.abs)) / N.amax(Q.flat().map(Math.abs));
    // simulations from an initial pitch disturbance: open loop, LQR, LQG (estimator starts at zero)
    const x0 = [0, 0, 0, rad(i.th0_deg)], tEnd = Math.max(10, 8 / Math.max(0.05, -eigMax(Acl))), ns = Math.ceil(N.clamp(2 * tEnd * N.amax([...modes(Acl), ...modes(Aest), ...modes(P.A)].map((m) => m.wn)), 1500, 40000)), ol = N.rk4((t, x) => N.matvec(P.A, x), 0, x0, tEnd, ns), cl = N.rk4((t, x) => N.matvec(Acl, x), 0, x0, tEnd, ns);
    const A8 = N.zeros(8), BK = N.matmul(P.B, K), LC = N.matmul(Lk, Cm); for (let p = 0; p < 4; p++) for (let q = 0; q < 4; q++) { A8[p][q] = P.A[p][q]; A8[p][4 + q] = -BK[p][q]; A8[4 + p][q] = LC[p][q]; A8[4 + p][4 + q] = P.A[p][q] - BK[p][q] - LC[p][q]; }
    const lg = N.rk4((t, x) => N.matvec(A8, x), 0, [...x0, 0, 0, 0, 0], tEnd, ns), ecl = modes(Acl), eol = modes(P.A), ee = modes(Aest), cost = N.dot(x0, N.matvec(Pr, x0));
    const uPk = N.amax(cl.y.map((x) => Math.abs(N.dot(K[0], x)))) / cs, lim = (v) => N.clamp(v, -1e4, 1e4);
    if (eigMax(P.A) > 0) warnings.push('The open-loop plant is unstable; the regulator stabilises it.');
    if (ricErr > 1e-6) warnings.push('The Riccati solution has a large residual; the model may be badly scaled or not stabilisable.');
    if (N.amax(ee.map((m) => m.wn)) > 60) warnings.push('The Kalman filter has very fast poles: the assumed sensor noise is small relative to the process noise, so it will pass sensor noise to the controls.');
    return {
      kpis: [
        { key: 'lqr_min_damping', label: 'Lowest closed-loop damping ratio', value: N.amin(ecl.map((m) => (m.wn > 1e-9 ? m.zeta : 1))), unit: '-', status: 'ok' }, { key: 'ol_min_damping', label: 'Lowest open-loop damping ratio', value: N.amin(eol.map((m) => (m.wn > 1e-9 ? m.zeta : 1))), unit: '-' },
        { key: 'lqr_slowest_s', label: 'Slowest closed-loop time constant', value: -1 / eigMax(Acl), unit: 's' }, { key: 'est_slowest_s', label: 'Slowest estimator time constant', value: -1 / eigMax(Aest), unit: 's' },
        { key: 'K_theta', label: 'Attitude gain on the pitch control', value: K[0][3], unit: rot ? '1/rad' : 'rad/rad' }, { key: 'K_q', label: 'Pitch-rate gain on the pitch control', value: K[0][2], unit: 's' },
        { key: 'peak_control', label: 'Peak pitch control in the recovery', value: uPk, unit: rot ? '%' : 'deg', status: uPk < 3 * i.c1_max ? 'ok' : 'warn' },
        { key: 'cost_J', label: 'Quadratic cost of the recovery', value: cost, unit: '-' }, { key: 'riccati_residual', label: 'Riccati equation residual', value: ricErr, unit: '-', status: ricErr < 1e-6 ? 'ok' : 'warn' },
        { key: 'est_theta_rms_deg', label: 'Steady-state attitude estimate error', value: deg(Math.sqrt(Math.max(Pe[3][3], 0))), unit: 'deg' },
      ],
      plots: [
        { type: 'line', title: 'Recovery from a pitch disturbance', xlabel: 'Time [s]', ylabel: 'Pitch attitude [deg]', series: [{ name: 'Open loop', x: thin(ol.t), y: thin(ol.y.map((x) => lim(deg(x[3])))) }, { name: 'LQR (full state)', x: thin(cl.t), y: thin(cl.y.map((x) => deg(x[3]))) }, { name: 'LQG (estimated state)', x: thin(lg.t), y: thin(lg.y.map((x) => lim(deg(x[3])))) }] },
        { type: 'line', title: 'Pitch control effort', xlabel: 'Time [s]', ylabel: rot ? 'Control [%]' : 'Deflection [deg]', series: [{ name: 'LQR', x: thin(cl.t), y: thin(cl.y.map((x) => -N.dot(K[0], x) / cs)) }, { name: 'LQG', x: thin(lg.t), y: thin(lg.y.map((x) => lim(-N.dot(K[0], x.slice(4)) / cs))) }] },
        { type: 'line', title: 'Pole map', xlabel: 'Real part [1/s]', ylabel: 'Imaginary part [rad/s]', series: [{ name: 'Open loop', x: eol.map((m) => m.re), y: eol.map((m) => m.im), style: 'points' }, { name: 'LQR closed loop', x: ecl.map((m) => m.re), y: ecl.map((m) => m.im), style: 'points' }, { name: 'Kalman filter', x: ee.map((m) => m.re), y: ee.map((m) => m.im), style: 'points' }], annotations: [{ x: 0, label: 'Stability boundary' }] },
      ],
      tables: [{ title: 'State-feedback gain K (u = −K·x, x = [u, w, q, θ])', columns: ['Control', 'u', 'w', 'q', 'θ'], rows: K.map((r, k) => [k ? 'Thrust / collective' : 'Pitch control', ...r]) }, { title: 'Kalman gain L (columns: q, θ measurements)', columns: ['State', 'from q', 'from θ'], rows: Lk.map((r, k) => [['u', 'w', 'q', 'θ'][k], ...r]) }],
      warnings, models: [`Plant: ${P.src}`, 'LQR controller (algebraic Riccati equation by the matrix sign function)', 'LQG controller: steady-state Kalman filter by the dual Riccati equation', 'Separation principle: regulator and estimator poles designed independently'],
      assumptions: ['Linear time-invariant plant, all states weighted by the Bryson rule (1 / acceptable value²)', 'No actuator dynamics or limits in this design model: check the result in the classical loop analysis', 'LQR has guaranteed margins with full-state feedback; LQG has none, so its loop margins must be verified', 'The Bryson tolerances and sensor-noise levels are illustrative defaults: set them from the actual requirements and sensors'],
    };
  },
  verify() {
    const r = N.lqr([[0, 1], [0, 0]], [[0], [1]], [[1, 0], [0, 1]], [[1]]), b = Object.fromEntries(lqrAn.inputs.map((f) => [f.key, f.default])), o = N.kv(lqrAn.run(b));
    return [N.check('Double integrator LQR gain k1 = 1', r.K[0][0], 1, 1e-8, 'Closed-form Riccati solution, Q = I, R = 1'), N.check('Double integrator LQR gain k2 = √3', r.K[0][1], Math.sqrt(3), 1e-8, 'Closed-form Riccati solution'),
      N.check('Riccati matrix P11 = √3', r.P[0][0], Math.sqrt(3), 1e-8, 'Closed-form Riccati solution'), N.check('Riccati residual of the aircraft design', o.riccati_residual, 0, 1e-7, "A'P + PA − PBR⁻¹B'P + Q = 0")];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.peak_control > 3 * i.c1_max) out.push({ severity: 'advise', title: 'Control effort is large for this disturbance', detail: `Peak ${o.peak_control.toFixed(1)} against an acceptable ${i.c1_max}.`, action: 'Increase the control weighting (smaller acceptable deflection) or relax the state tolerances.', basis: 'Bryson rule trade-off' });
    out.push({ severity: 'info', title: 'Verify robustness of the LQG loop', detail: 'Optimal gains do not guarantee stability margins once the Kalman filter and actuator are in the loop.', action: 'Check gain and phase margins with the classical loop analysis, or apply loop-transfer recovery.', basis: 'Doyle (1978): LQG has no guaranteed margins' });
    return out;
  },
};

// ---- 3. autopilot ---------------------------------------------------------------------------
const autopilot = {
  id: 'autopilot', title: 'Autopilot: altitude, airspeed and heading hold', fidelity: 'numerical',
  summary: 'Time simulation of cascaded altitude-hold (through pitch attitude), autothrottle and heading-hold (through bank angle) loops with turbulence, sensor noise and actuator limits: tracking errors and saturation events.',
  equations: ['State-space equations', 'Proportional–integral–derivative control equations', 'Actuator dynamic equations', 'Transfer function equations'],
  applicable: (c) => (c.wing.S_m2 > 0 ? true : 'This autopilot is for fixed-wing flight; use the hover cascade analysis for rotorcraft and multirotors.'),
  inputs: [
    { key: 'dh', label: 'Altitude step command', unit: 'm', default: 100, min: -1000, max: 1000, group: 'Commands' }, { key: 'dV', label: 'Airspeed step command', unit: 'm/s', default: 0, min: -30, max: 30, group: 'Commands' },
    { key: 'dpsi_deg', label: 'Heading step command', unit: 'deg', default: 30, min: -180, max: 180, group: 'Commands' },
    { key: 'w_pitch', label: 'Pitch-loop crossover', unit: 'rad/s', default: 2, min: 0.05, max: 50, group: 'Design bandwidths' }, { key: 'w_alt_frac', label: 'Altitude-loop bandwidth / pitch-loop', unit: '-', default: 0.06, min: 0.02, max: 0.3, group: 'Design bandwidths', help: 'Outer loops should be several times slower than the loop they command; the vertical-speed loop sits between them at 0.3 × the pitch loop' },
    { key: 'vs_lim', label: 'Vertical-speed command limit', unit: 'm/s', default: 8, min: 0.5, max: 40, group: 'Limits' },
    { key: 'w_speed', label: 'Autothrottle bandwidth', unit: 'rad/s', default: 0.15, min: 0.01, max: 2, group: 'Design bandwidths' }, { key: 'w_hdg', label: 'Heading-loop bandwidth', unit: 'rad/s', default: 0.12, min: 0.01, max: 1, group: 'Design bandwidths' },
    { key: 'tau_roll', label: 'Closed-loop bank response time constant', unit: 's', default: 1, min: 0.1, max: 10, group: 'Design bandwidths' },
    { key: 'theta_lim_deg', label: 'Pitch command limit', unit: 'deg', default: 8, min: 1, max: 30, group: 'Limits' }, { key: 'bank_lim_deg', label: 'Bank command limit', unit: 'deg', default: 25, min: 5, max: 60, group: 'Limits' },
    { key: 'de_lim_deg', label: 'Elevator travel', unit: 'deg', default: 25, min: 1, max: 40, group: 'Limits' }, { key: 'de_rate_dps', label: 'Elevator rate limit', unit: 'deg/s', default: 60, min: 1, max: 500, group: 'Limits' },
    { key: 'tau_eng', label: 'Engine thrust time constant', unit: 's', default: 3, min: 0.05, max: 15, group: 'Limits', help: 'About 2–5 s turbofan, 0.5–1 s propeller, 0.1 s electric' },
    { key: 'throttle_trim', label: 'Trim throttle', unit: '-', default: 0.5, min: 0, max: 1, group: 'Limits' },
    { key: 'sigma_w', label: 'Vertical gust RMS', unit: 'm/s', default: 1, min: 0, max: 10, group: 'Disturbances' }, { key: 'L_w', label: 'Turbulence scale', unit: 'm', default: 533, min: 5, max: 3000, group: 'Disturbances' },
    { key: 'sigma_h', label: 'Altitude sensor noise', unit: 'm', default: 1, min: 0, max: 20, group: 'Disturbances' }, { key: 'seed', label: 'Random seed', unit: '', default: 7, min: 1, max: 1e9, step: 1, discrete: true, group: 'Disturbances' },
    { key: 't_end', label: 'Simulated time', unit: 's', default: 120, min: 10, max: 900, group: 'Numerics' }, { key: 'dt', label: 'Control and integration step', unit: 's', default: 0.02, min: 0.001, max: 0.2, group: 'Numerics' },
    ...PLANT.filter((f) => f.group !== 'Rotorcraft' && f.key !== 'vehicle')],
  defaults: (c, up, d) => {
    const p = plantDefaults(c, up, d), i = { ...Object.fromEntries(PLANT.map((f) => [f.key, f.default])), ...p }, sp = fastMode(plant(i, null).A), at = isa(c.atm.alt_m), q = 0.5 * at.rho * c.flight.V_ms ** 2, CL = d.W / (q * c.wing.S_m2);
    const thr = p.T_max > 0 ? N.clamp((q * c.wing.S_m2 * (i.CD0 + i.k_ind * CL * CL)) / p.T_max, 0.02, 0.95) : 0.5, pt = c.prop.type;
    return { ...p, vehicle: undefined, w_pitch: N.clamp(1.3 * (sp?.wn ?? 2), 0.5, 6), de_lim_deg: c.controls.de_max_deg, de_rate_dps: c.controls.rate_max_dps, tau_eng: pt === 'electric' ? 0.15 : pt === 'piston' || pt === 'turboprop' ? 0.8 : 3, throttle_trim: thr, sigma_w: Math.max(0.3, 0.5 * c.atm.turb_intensity * c.flight.V_ms), dh: c.mass.mtow_kg > 5700 ? 100 : 30, vs_lim: N.clamp(0.04 * c.flight.V_ms, 1, 10), t_end: 120 };
  },
  run(i, ctx) {
    const pi = { ...i, vehicle: 'aeroplane' }, P = plant(pi, ctx), u0 = P.u0, b1 = P.B.map((r) => r[0]), b2 = P.B.map((r) => r[1] ?? 0), warnings = [];
    const act = { w_act: 30, z_act: 0.7, Ts: i.dt, delay: 0 }, g = tunePID(openLoop(P, act), i.w_pitch, 60), wh = i.w_alt_frac * i.w_pitch, wvs = 0.3 * i.w_pitch, Tg = 1 / Math.max(-P.A[1][1], 1e-3), Kh = (wvs * Tg) / u0, Khi = wvs / u0, XdT = b2[0] || 1e-9, Kv = i.w_speed / XdT, Kvi = (Kv * i.w_speed) / 4, Kpsi = (i.w_hdg * u0) / G0;
    const n = Math.round(i.t_end / i.dt), dt = i.dt, rng = N.rng(i.seed), tau = i.L_w / u0, thLim = rad(i.theta_lim_deg), deLim = rad(i.de_lim_deg), deRate = rad(i.de_rate_dps), bkLim = rad(i.bank_lim_deg), psiRef = rad(i.dpsi_deg);
    // vertical speed follows pitch attitude through the flight-path lag Tγ = −1/Zw: ḣ/θ = u0/(Tγ·s + 1). The PI zero cancels that lag, leaving the loop u0·Kh/(Tγ·s) with crossover ωvs
    let hf = 0, x = [0, 0, 0, 0], h = 0, de = 0, thr = 0, ih = 0, ith = 0, iv = 0, wg = 0, phi = 0, psi = 0, satE = 0, satT = 0, wasE = false, wasT = false; const T = [], H = [], U = [], TH = [], DE = [], TR = [], PS = [], PH = [], t0 = 5;
    for (let k = 0; k <= n; k++) {
      const t = k * dt, hRef = t >= t0 ? i.dh : 0, vRef = t >= t0 ? i.dV : 0, pRef = t >= t0 ? psiRef : 0; hf += (h + i.sigma_h * N.randn(rng) - hf) * Math.min(1, 4 * wh * dt); const hm = hf; // first-order altitude filter at four times the loop bandwidth
      // outer loops
      const vsC = N.clamp(wh * (hRef - hm), -i.vs_lim, i.vs_lim), eh = vsC - (u0 * x[3] - x[1]), thU = Kh * eh + Khi * ih, thC = N.clamp(thU, -thLim, thLim); if (thU === thC) ih += eh * dt;
      const eth = thC - x[3], deU = g.Kp * eth + g.Ki * ith - g.Kd * x[2], deC = N.clamp(deU, -deLim, deLim); if (deU === deC) ith += eth * dt;
      const ev = vRef - x[0], trU = Kv * ev + Kvi * iv, trC = N.clamp(trU, -i.throttle_trim, 1 - i.throttle_trim); if (trU === trC) iv += ev * dt;
      const phC = N.clamp(Kpsi * (pRef - psi), -bkLim, bkLim);
      // actuators: rate-limited elevator servo, first-order engine
      const dd = N.clamp((deC - de) * 30, -deRate, deRate), sE = Math.abs(deU) > deLim || Math.abs((deC - de) * 30) > 1.0001 * deRate, sT = trU !== trC; de += dd * dt; thr += ((trC - thr) / i.tau_eng) * dt;
      if (sE && !wasE) satE++; wasE = sE; if (sT && !wasT) satT++; wasT = sT;
      T.push(t); H.push(h); U.push(x[0]); TH.push(x[3]); DE.push(de); TR.push(thr + i.throttle_trim); PS.push(psi); PH.push(phi);
      // plant with gust (first-order Dryden-like filter, exact discrete update), RK4
      wg = wg * Math.exp(-dt / tau) + i.sigma_w * Math.sqrt(1 - Math.exp((-2 * dt) / tau)) * N.randn(rng);
      const f = (y) => { const d = N.matvec(P.A, y); for (let p = 0; p < 4; p++) d[p] += b1[p] * de + b2[p] * thr - P.A[p][1] * wg; return d; }, k1 = f(x), k2 = f(N.vadd(x, k1, dt / 2)), k3 = f(N.vadd(x, k2, dt / 2)), k4 = f(N.vadd(x, k3, dt));
      h += (u0 * x[3] - x[1]) * dt; x = x.map((v, q) => v + (dt / 6) * (k1[q] + 2 * k2[q] + 2 * k3[q] + k4[q])); phi += ((phC - phi) / i.tau_roll) * dt; psi += ((G0 * Math.tan(phi)) / u0) * dt;
      if (!fin(h) || Math.abs(h) > 1e6) { warnings.push('The simulation diverged: the chosen bandwidths are not stable for this aircraft.'); break; }
    }
    const k0 = Math.round(t0 / dt), seg = (a) => a.slice(k0), sh = i.dh ? stepMetrics(seg(T), seg(H), i.dh) : null, sp = i.dpsi_deg ? stepMetrics(seg(T), seg(PS), psiRef) : null, tail = Math.round(0.7 * T.length);
    const rms = (a, ref) => Math.sqrt(N.mean(a.slice(tail).map((v) => (v - ref) ** 2))), hR = rms(H, i.dh), vR = rms(U, i.dV), deR = deg(N.std(DE.slice(tail)));
    if (satE) warnings.push(`The elevator hit a travel or rate limit ${satE} time(s).`); if (satT) warnings.push(`The throttle saturated ${satT} time(s): the climb or speed change is thrust-limited.`);
    if (i.w_alt_frac > 0.15) warnings.push('The altitude loop is less than twice as slow as the vertical-speed loop: expect interaction and overshoot.');
    return {
      kpis: [
        { key: 'alt_rms_m', label: 'Altitude-hold RMS error in turbulence', value: hR, unit: 'm', status: hR < Math.max(5, 3 * i.sigma_h) ? 'ok' : 'warn' }, { key: 'alt_rise_s', label: 'Altitude capture rise time', value: sh ? sh.rise : 0, unit: 's' },
        { key: 'alt_overshoot_pct', label: 'Altitude overshoot', value: sh ? sh.overshoot : 0, unit: '%', status: !sh || sh.overshoot < 15 ? 'ok' : 'warn' }, { key: 'speed_rms_ms', label: 'Airspeed-hold RMS error', value: vR, unit: 'm/s' },
        { key: 'speed_dip_ms', label: 'Largest airspeed excursion', value: N.amax(U.map((v, k) => Math.abs(v - (T[k] >= t0 ? i.dV : 0)))), unit: 'm/s' },
        { key: 'hdg_settle_s', label: 'Heading capture settling time', value: sp ? sp.settle - t0 : 0, unit: 's' }, { key: 'hdg_overshoot_pct', label: 'Heading overshoot', value: sp ? sp.overshoot : 0, unit: '%' },
        { key: 'elevator_rms_deg', label: 'Elevator activity (RMS)', value: deR, unit: 'deg' }, { key: 'saturation_events', label: 'Elevator saturation events', value: satE, unit: '-', status: satE ? 'warn' : 'ok' }, { key: 'throttle_sat_events', label: 'Throttle saturation events', value: satT, unit: '-' },
        { key: 'K_vs', label: 'Vertical-speed gain', value: deg(Kh), unit: 'deg per m/s' }, { key: 'w_alt_rads', label: 'Altitude-loop bandwidth', value: wh, unit: 'rad/s' }, { key: 'K_throttle', label: 'Autothrottle gain', value: Kv, unit: '1/(m/s)' }, { key: 'K_heading', label: 'Heading gain', value: Kpsi, unit: 'rad bank/rad' },
      ].filter((k) => fin(k.value)),
      plots: [
        { type: 'line', title: 'Altitude tracking', xlabel: 'Time [s]', ylabel: 'Altitude change [m]', series: [{ name: 'Altitude', x: thin(T), y: thin(H) }, { name: 'Command', x: [0, t0, t0, i.t_end], y: [0, 0, i.dh, i.dh], style: 'dash' }] },
        { type: 'line', title: 'Airspeed', xlabel: 'Time [s]', ylabel: 'Airspeed change [m/s]', series: [{ name: 'Airspeed', x: thin(T), y: thin(U) }, { name: 'Command', x: [0, t0, t0, i.t_end], y: [0, 0, i.dV, i.dV], style: 'dash' }] },
        { type: 'line', title: 'Heading and bank', xlabel: 'Time [s]', ylabel: 'Angle [deg]', series: [{ name: 'Heading', x: thin(T), y: thin(PS.map(deg)) }, { name: 'Bank', x: thin(T), y: thin(PH.map(deg)) }] },
        { type: 'line', title: 'Pitch attitude and elevator', xlabel: 'Time [s]', ylabel: 'Angle [deg]', series: [{ name: 'Pitch attitude', x: thin(T), y: thin(TH.map(deg)) }, { name: 'Elevator', x: thin(T), y: thin(DE.map(deg)) }] },
        { type: 'line', title: 'Throttle', xlabel: 'Time [s]', ylabel: 'Throttle [-]', series: [{ name: 'Throttle', x: thin(T), y: thin(TR) }] },
      ],
      warnings, models: [`Plant: ${P.src}`, 'Autopilot model: cascaded altitude → vertical speed → pitch attitude → elevator, PI autothrottle, heading → bank', 'Gains derived from the chosen loop bandwidths; vertical-speed PI with its zero on the flight-path lag −1/Zw', 'Seeded first-order gust and Gaussian altitude-sensor noise with a first-order measurement filter', 'Kinematic coordinated-turn heading model with a first-order bank response'],
      assumptions: ['Linear longitudinal dynamics; lateral axis reduced to turn kinematics (the lateral modes are in Suite 4 and the stability-augmentation analysis)', 'Elevator servo as a 30 rad/s rate-limited lag; first-order engine response', 'No mode logic, flight-director or envelope protection', 'Loop-bandwidth ratios, command limits, servo bandwidth and engine time constants are typical design values, not data for a specific autopilot'],
    };
  },
  verify() {
    const b = Object.fromEntries(autopilot.inputs.map((f) => [f.key, f.default])), o = N.kv(autopilot.run({ ...b, sigma_w: 0, sigma_h: 0, t_end: 200 })), tau = 1, wh = 0.1, u0 = 231, Kp = (wh * u0) / G0;
    // heading loop: ψ̇ = g·φ/V, φ̇ = (Kψ(ψref − ψ) − φ)/τ → second order with ωn² = ωh/τ, ζ = 1/(2·sqrt(ωh·τ))
    const z = 1 / (2 * Math.sqrt(wh * tau)), r = N.rk4((t, y) => [(G0 * y[1]) / u0, (Kp * (0.01 - y[0]) - y[1]) / tau], 0, [0, 0], 120, 6000), sm = stepMetrics(r.t, r.y.map((y) => y[0]), 0.01);
    return [N.check('Altitude command is captured with zero steady-state error', o.alt_rms_m / 100 + 1, 1, 2e-3, 'Type-1 altitude loop'),
      N.check('Altitude capture is overdamped: ζ = ½·sqrt(ωvs/ωh) = 1.12', o.alt_overshoot_pct / 100 + 1, 1, 0.03, 'Altitude loop around a first-order vertical-speed loop: s² + ωvs·s + ωvs·ωh, no overshoot for ζ ≥ 1 (tolerance for the pitch-loop dynamics)'), N.check('Heading-loop overshoot of the equivalent second-order system', sm.overshoot + 1, 1 + (z < 1 ? 100 * Math.exp((-Math.PI * z) / Math.sqrt(1 - z * z)) : 0), 1e-3, 'ζ = 1/(2·sqrt(ωh·τ)) for small bank angles')];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (o.saturation_events > 0 || o.throttle_sat_events > 0) out.push({ severity: 'advise', title: 'Control saturation during the manoeuvre', detail: `Elevator ${o.saturation_events}, throttle ${o.throttle_sat_events} event(s).`, action: 'Rate-limit the commands (altitude capture profile), reduce outer-loop bandwidth or accept a thrust-limited climb.', basis: 'Actuator limits' });
    if (o.alt_overshoot_pct > 15) out.push({ severity: 'warn', title: 'Altitude capture overshoots', detail: `${o.alt_overshoot_pct.toFixed(0)}% overshoot.`, action: 'Lower the altitude-loop bandwidth ratio or the vertical-speed limit.', basis: 'Loop separation rule' });
    if (o.elevator_rms_deg > 1) out.push({ severity: 'info', title: 'High elevator activity in turbulence', detail: `${o.elevator_rms_deg.toFixed(2)}° RMS.`, action: 'Filter the altitude sensor or reduce bandwidth: actuator duty cycle drives wear and hydraulic/electrical power.', basis: 'Actuator duty cycle' });
    return out;
  },
};

// ---- 4. stability augmentation ----------------------------------------------------------------
const sas = {
  id: 'sas', title: 'Stability augmentation: yaw damper with washout and pitch damper', fidelity: 'reduced-order',
  summary: 'How yaw-rate feedback to the rudder (through a washout filter) raises Dutch-roll damping and how pitch-rate feedback to the elevator raises short-period damping, with the gains that meet a damping target.',
  equations: ['State-space equations', 'Root locus relations', 'Transfer function equations', 'Lyapunov stability equations'],
  applicable: (c) => (c.wing.S_m2 > 0 ? true : 'Yaw and pitch dampers here are for fixed-wing aircraft; rotorcraft rate damping is covered by the hover cascade.'),
  inputs: [
    { key: 'zeta_dr_target', label: 'Dutch-roll damping target', unit: '-', default: 0.3, min: 0.05, max: 0.9, group: 'Targets' }, { key: 'zeta_sp_target', label: 'Short-period damping target', unit: '-', default: 0.6, min: 0.2, max: 1, group: 'Targets' },
    { key: 'tau_wo', label: 'Washout time constant', unit: 's', default: 3, min: 0.2, max: 20, group: 'Yaw damper', help: 'Passes Dutch-roll yaw rate but not the steady yaw rate of a turn; about 1–4 s' },
    { key: 'K_r_max', label: 'Largest yaw-damper gain considered', unit: 's', default: 3, min: 0.1, max: 20, group: 'Yaw damper' }, { key: 'K_q_max', label: 'Largest pitch-damper gain considered', unit: 's', default: 2, min: 0.05, max: 20, group: 'Pitch damper' },
    { key: 'b', label: 'Wing span', unit: 'm', default: 34, min: 0.05, group: 'Lateral plant' }, { key: 'Ixx', label: 'Roll inertia', unit: 'kg·m²', default: 1.28e6, min: 1e-6, group: 'Lateral plant' }, { key: 'Izz', label: 'Yaw inertia', unit: 'kg·m²', default: 4.9e6, min: 1e-6, group: 'Lateral plant' },
    { key: 'CYb', label: 'CYβ', unit: '1/rad', default: -0.8, max: 0, group: 'Lateral plant' }, { key: 'Clb', label: 'Clβ', unit: '1/rad', default: -0.2, group: 'Lateral plant' }, { key: 'Cnb', label: 'Cnβ', unit: '1/rad', default: 0.13, group: 'Lateral plant' },
    { key: 'Clp', label: 'Clp', unit: '1/rad', default: -0.5, max: 0, group: 'Lateral plant' }, { key: 'Cnp', label: 'Cnp', unit: '1/rad', default: -0.07, group: 'Lateral plant' }, { key: 'Clr', label: 'Clr', unit: '1/rad', default: 0.2, group: 'Lateral plant' },
    { key: 'Cnr', label: 'Cnr', unit: '1/rad', default: -0.25, max: 0, group: 'Lateral plant' }, { key: 'CYdr', label: 'CYδr', unit: '1/rad', default: 0.3, group: 'Lateral plant' }, { key: 'Cldr', label: 'Clδr', unit: '1/rad', default: 0.02, group: 'Lateral plant' },
    { key: 'Cndr', label: 'Cnδr', unit: '1/rad', default: -0.12, max: 0, group: 'Lateral plant' },
    ...PLANT.filter((f) => f.group !== 'Rotorcraft' && f.key !== 'vehicle')],
  defaults: (c, up, d) => { const e = estimate(c, d), m = c.mass.mtow_kg; return { ...plantDefaults(c, up, d), vehicle: undefined, b: c.wing.b_m, Ixx: c.mass.Ixx > 0 ? c.mass.Ixx : m * (0.125 * c.wing.b_m) ** 2, Izz: c.mass.Izz > 0 ? c.mass.Izz : m * (0.1 * (c.wing.b_m + c.fuselage.len_m)) ** 2, CYb: e.CYb, Clb: e.Clb, Cnb: up.flightdyn?.Cn_beta ?? e.Cnb, Clp: e.Clp, Cnp: e.Cnp, Clr: e.Clr, Cnr: e.Cnr, CYdr: e.CYdr, Cldr: e.Cldr, Cndr: e.Cndr }; },
  run(i, ctx) {
    const pi = { ...i, vehicle: 'aeroplane' }, Pl = plant(pi, ctx), f = ctx?.up?.flightdyn, warnings = []; let A, Bm, src = 'built-in model from the listed derivatives';
    if (i.use_upstream && f && Array.isArray(f.A_lat) && f.A_lat.length === 4 && Array.isArray(f.B_lat)) { A = f.A_lat; Bm = f.B_lat.map((r) => r[1]); src = 'Suite 4 linear model'; }
    else { const at = isa(i.alt_m), QS = 0.5 * at.rho * i.V * i.V * i.S, k = i.b / (2 * i.V), L = (c) => (QS * i.b * c) / i.Ixx, Nn = (c) => (QS * i.b * c) / i.Izz, Y = (c) => (QS * c) / (i.mass_kg * i.V);
      A = [[Y(i.CYb), 0, -1, G0 / i.V], [L(i.Clb), L(i.Clp) * k, L(i.Clr) * k, 0], [Nn(i.Cnb), Nn(i.Cnp) * k, Nn(i.Cnr) * k, 0], [0, 1, 0, 0]]; Bm = [Y(i.CYdr), L(i.Cldr), Nn(i.Cndr), 0]; }
    // yaw damper: δr = K·(r − x_w), ẋ_w = (r − x_w)/τ  → 5 states
    const yd = (K) => { const M = N.zeros(5); for (let p = 0; p < 4; p++) { for (let q = 0; q < 4; q++) M[p][q] = A[p][q]; M[p][2] += Bm[p] * K; M[p][4] -= Bm[p] * K; } M[4][2] = 1 / i.tau_wo; M[4][4] = -1 / i.tau_wo; return M; };
    const dr0 = leastDamped(A, 0.1), sgn = Bm[2] < 0 ? 1 : -1, Ks = N.linspace(0, i.K_r_max, 31), zd = Ks.map((K) => leastDamped(yd(sgn * K), 0.1)?.zeta ?? 1), loc = { x: [], y: [] }; for (const K of Ks) for (const l of N.eig(yd(sgn * K))) { loc.x.push(l[0]); loc.y.push(l[1]); }
    let Kr = NaN; for (let k = 1; k < Ks.length; k++) if (zd[k - 1] < i.zeta_dr_target && zd[k] >= i.zeta_dr_target) { Kr = Ks[k - 1] + ((i.zeta_dr_target - zd[k - 1]) / (zd[k] - zd[k - 1])) * (Ks[k] - Ks[k - 1]); break; }
    if (zd[0] >= i.zeta_dr_target) Kr = 0; const KrUse = fin(Kr) ? Kr : Ks[N.argmax(zd)], Ayd = yd(sgn * KrUse), dr1 = leastDamped(Ayd, 0.1);
    // pitch damper: δe = Kq·q
    const b1 = Pl.B.map((r) => r[0]), pd = (K) => Pl.A.map((r, p) => r.map((v, q) => v + (q === 2 ? b1[p] * K : 0))), sp0 = fastMode(Pl.A), sq = b1[2] < 0 ? 1 : -1, Kqs = N.linspace(0, i.K_q_max, 31), zs = Kqs.map((K) => fastMode(pd(sq * K))?.zeta ?? 1);
    let Kq = NaN; for (let k = 1; k < Kqs.length; k++) if (zs[k - 1] < i.zeta_sp_target && zs[k] >= i.zeta_sp_target) { Kq = Kqs[k - 1] + ((i.zeta_sp_target - zs[k - 1]) / (zs[k] - zs[k - 1])) * (Kqs[k] - Kqs[k - 1]); break; }
    if (zs[0] >= i.zeta_sp_target) Kq = 0; const KqUse = fin(Kq) ? Kq : Kqs[N.argmax(zs)], sp1 = fastMode(pd(sq * KqUse));
    // responses to an initial sideslip of 2°
    const x0 = [rad(2), 0, 0, 0], tE = dr0 ? Math.min(60, 6 * (TAU / dr0.im)) : 30, o = N.rk4((t, x) => N.matvec(A, x), 0, x0, tE, 1200), c = N.rk4((t, x) => N.matvec(Ayd, x), 0, [...x0, 0], tE, 1200), lim = (v) => N.clamp(v, -90, 90);
    if (!fin(Kr)) warnings.push(`The Dutch-roll damping target cannot be reached within the gain range; the best achievable is ${N.amax(zd).toFixed(2)}.`);
    if (!fin(Kq)) warnings.push(`The short-period damping target cannot be reached within the gain range; the best achievable is ${N.amax(zs).toFixed(2)}.`);
    if (eigMax(Ayd) > 0.02) warnings.push('A slow real root remains unstable with the yaw damper (spiral mode): it is unaffected by washed-out yaw-rate feedback and needs a wing-leveller.');
    return {
      kpis: [
        { key: 'dr_zeta_open', label: 'Dutch-roll damping, bare airframe', value: dr0 ? dr0.zeta : NaN, unit: '-' }, { key: 'dr_zeta_sas', label: 'Dutch-roll damping with the yaw damper', value: dr1 ? dr1.zeta : 1, unit: '-', status: (dr1 ? dr1.zeta : 1) >= i.zeta_dr_target - 1e-3 ? 'ok' : 'warn' },
        { key: 'K_yaw_damper', label: 'Yaw-damper gain', value: sgn * KrUse, unit: 'rad per rad/s' }, { key: 'dr_omega_sas', label: 'Dutch-roll frequency with the yaw damper', value: dr1 ? dr1.wn : NaN, unit: 'rad/s' },
        { key: 'sp_zeta_open', label: 'Short-period damping, bare airframe', value: sp0 ? sp0.zeta : NaN, unit: '-' }, { key: 'sp_zeta_sas', label: 'Short-period damping with the pitch damper', value: sp1 ? sp1.zeta : 1, unit: '-', status: (sp1 ? sp1.zeta : 1) >= i.zeta_sp_target - 1e-3 ? 'ok' : 'warn' },
        { key: 'K_pitch_damper', label: 'Pitch-damper gain', value: sq * KqUse, unit: 'rad per rad/s' },
        { key: 'rudder_per_dps', label: 'Rudder used per deg/s of yaw rate', value: KrUse, unit: 'deg/(deg/s)' },
      ].filter((k) => fin(k.value)),
      plots: [
        { type: 'line', title: 'Damping against feedback gain', xlabel: 'Gain magnitude [s]', ylabel: 'Damping ratio [-]', series: [{ name: 'Dutch roll (yaw damper)', x: Ks, y: zd }, { name: 'Short period (pitch damper)', x: Kqs, y: zs }], annotations: [{ y: i.zeta_dr_target, label: 'Dutch-roll target' }] },
        { type: 'line', title: 'Yaw-damper root locus', xlabel: 'Real part [1/s]', ylabel: 'Imaginary part [rad/s]', series: [{ name: 'Closed-loop poles', x: loc.x, y: loc.y, style: 'points' }], annotations: [{ x: 0, label: 'Stability boundary' }] },
        { type: 'line', title: 'Sideslip after a 2° disturbance', xlabel: 'Time [s]', ylabel: 'Sideslip [deg]', series: [{ name: 'Bare airframe', x: thin(o.t), y: thin(o.y.map((x) => lim(deg(x[0])))) }, { name: 'Yaw damper on', x: thin(c.t), y: thin(c.y.map((x) => lim(deg(x[0])))) }] },
      ],
      warnings, models: [`Lateral plant: ${src}`, `Longitudinal plant: ${Pl.src}`, 'Stability augmentation model: yaw-rate feedback through a first-order washout; pitch-rate feedback', 'Gain selection by root-locus sweep to a damping target'],
      assumptions: ['Ideal rate gyros and actuators (include them with the classical loop analysis)', 'Single-loop designs; no aileron–rudder interconnect or turn coordination', 'Damping targets, washout time constant and gain ranges are typical design values'],
    };
  },
  verify() {
    const b = Object.fromEntries(sas.inputs.map((f) => [f.key, f.default])), o = N.kv(sas.run({ ...b, zeta_dr_target: 0.3 }));
    // rate feedback on a second-order system: s² + (2ζω + bK)s + ω² → ζ' = ζ + bK/(2ω)
    const w = 2, z = 0.1, bb = 1.5, K = 0.8, A = [[0, 1], [-w * w, -2 * z * w - bb * K]], m = modes(A)[0];
    return [N.check('Rate feedback raises second-order damping by b·K/(2ω)', m.zeta, z + (bb * K) / (2 * w), 1e-9, 'Closed-form characteristic polynomial'), N.check('Selected yaw-damper gain meets the damping target', o.dr_zeta_sas, 0.3, 0.02, 'Root-locus interpolation')];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.dr_zeta_open < 0.08) out.push({ severity: o.dr_zeta_sas >= i.zeta_dr_target - 1e-3 ? 'advise' : 'warn', title: 'Bare-airframe Dutch roll needs a yaw damper', detail: `Damping rises from ${o.dr_zeta_open.toFixed(3)} to ${o.dr_zeta_sas.toFixed(2)} with a gain of ${o.K_yaw_damper.toFixed(2)}.`, action: 'Treat the yaw damper as flight-critical if the bare airframe is below the minimum acceptable damping: assess its failure probability in Suite 22.', basis: 'Dutch-roll damping below the MIL-F-8785C Table VI Level 1 Category B minimum of 0.08' });
    if (o.sp_zeta_sas < i.zeta_sp_target - 1e-3) out.push({ severity: 'warn', title: 'Pitch damper cannot reach the target', detail: `Best short-period damping ${o.sp_zeta_sas.toFixed(2)}.`, action: 'Add angle-of-attack or normal-acceleration feedback, or increase tail volume.', basis: 'Short-period damping target' });
    return out;
  },
};

// ---- 5. control allocation ------------------------------------------------------------------
/** Weighted pseudo-inverse with redistribution: saturated effectors are frozen at their limit and the rest re-solved. */
export function allocate(B, v, w, lo, hi) {
  const m = B[0].length, n = B.length, u = new Array(m).fill(0), free = new Array(m).fill(true); let iter = 0;
  for (; iter < m + 1; iter++) {
    const idx = N.range(m).filter((k) => free[k]); if (!idx.length) break;
    const res = v.map((x, p) => x - N.sum(N.range(m).filter((k) => !free[k]).map((k) => B[p][k] * u[k]))), Bf = B.map((r) => idx.map((k) => r[k])), Wi = idx.map((k) => 1 / w[k]);
    const G = N.range(n, (p) => N.range(n, (q) => N.sum(idx.map((_, k) => Bf[p][k] * Wi[k] * Bf[q][k])) + (p === q ? 1e-12 : 0)));
    let lam; try { lam = N.solve(G, res); } catch { lam = N.lstsq(G, res, 1e-9); }
    let sat = false; idx.forEach((k, j) => { const x = Wi[j] * N.sum(Bf.map((r, p) => r[j] * lam[p])); u[k] = x; });
    for (const k of idx) if (u[k] > hi[k] + 1e-12 || u[k] < lo[k] - 1e-12) { u[k] = N.clamp(u[k], lo[k], hi[k]); free[k] = false; sat = true; }
    if (!sat) break;
  }
  const got = B.map((r) => N.dot(r, u));
  return { u, got, err: Math.sqrt(N.sum(got.map((g, p) => (g - v[p]) ** 2))), nSat: free.filter((x) => !x).length, iter };
}
const allocation = {
  id: 'allocation', title: 'Control allocation and effector failure', fidelity: 'analytical',
  summary: 'Distribution of commanded moments (and thrust) over redundant control effectors by a weighted pseudo-inverse with saturation redistribution, the attainable moment envelope, and re-allocation after an effector failure.',
  equations: ['State-space equations', 'Controllability and observability equations'],
  inputs: [
    { key: 'layout', label: 'Effector layout', type: 'select', options: ['aeroplane surfaces', 'multirotor', 'helicopter'], default: 'aeroplane surfaces', group: 'Effectors' },
    { key: 'n_rotors', label: 'Number of rotors (multirotor)', unit: '', default: 4, min: 3, max: 12, step: 1, discrete: true, group: 'Effectors' },
    { key: 'roll_cmd', label: 'Roll moment demand', unit: '% of nominal authority', default: 30, min: -150, max: 150, group: 'Demand' }, { key: 'pitch_cmd', label: 'Pitch moment demand', unit: '% of nominal authority', default: 40, min: -150, max: 150, group: 'Demand' },
    { key: 'yaw_cmd', label: 'Yaw moment demand', unit: '% of nominal authority', default: 20, min: -150, max: 150, group: 'Demand' }, { key: 'thrust_cmd', label: 'Thrust demand (rotorcraft)', unit: '% of hover thrust', default: 100, min: 0, max: 200, group: 'Demand' },
    { key: 'failed', label: 'Failed effector number (0 = none)', unit: '', default: 1, min: 0, max: 12, step: 1, discrete: true, group: 'Failure', help: 'Aeroplane: 1–2 elevators, 3–4 ailerons, 5 rudder, 6–7 spoilers. Multirotor: rotor index, counted from the front-right going clockwise' },
    { key: 'fail_mode', label: 'Failure mode', type: 'select', options: ['floating / lost', 'jammed at current position'], default: 'floating / lost', group: 'Failure' },
    { key: 'jam_pos', label: 'Jam position', unit: 'fraction of travel', default: 0.5, min: -1, max: 1, group: 'Failure' },
    { key: 'thrust_ratio', label: 'Maximum / hover thrust per rotor', unit: '-', default: 2, min: 1.05, max: 5, group: 'Effectors' },
    { key: 'w_spoiler', label: 'Spoiler usage weight', unit: '-', default: 4, min: 0.1, max: 100, group: 'Effectors', help: 'Higher weight makes the allocator avoid that effector (drag penalty)' },
    { key: 'kq_arm', label: 'Rotor torque / (thrust × arm)', unit: '-', default: 0.05, min: 0.005, max: 0.5, group: 'Effectors' },
  ],
  defaults: (c, up, d) => { const heli = c.meta.type === 'helicopter', mr = !heli && !(c.wing.S_m2 > 0) && c.rotor.R_m > 0; if (!heli && !mr) return { layout: 'aeroplane surfaces' }; const h = hoverEst(c); return heli ? { layout: 'helicopter' } : { layout: 'multirotor', n_rotors: Math.max(3, c.prop.n_eng), thrust_ratio: Math.max(1.1, (c.prop.T0_N * c.prop.n_eng) / d.W) || 2, kq_arm: N.clamp(h.Q / (h.T * h.arm), 0.005, 0.5) }; },
  run(i) {
    let B, names, lo, hi, w, axes, v, trim; const warnings = [];
    if (i.layout === 'multirotor') {
      const n = Math.round(i.n_rotors); names = N.range(n, (k) => `Rotor ${k + 1}`); axes = ['Thrust', 'Roll', 'Pitch', 'Yaw'];
      // rotor k at azimuth ψ (from the nose, clockwise seen from above), alternating spin; effector = thrust / hover thrust per rotor
      const az = N.range(n, (k) => ((k + 0.5) * TAU) / n); B = [N.range(n, () => 1 / n), az.map((a) => -Math.sin(a) / n), az.map((a) => Math.cos(a) / n), N.range(n, (k) => ((k % 2 ? -1 : 1) * i.kq_arm) / n)];
      lo = names.map(() => 0); hi = names.map(() => i.thrust_ratio); w = names.map(() => 1); trim = names.map(() => 1);
      const nomM = 0.25, nomN = 0.5 * i.kq_arm; // nominal authority: a quarter of hover thrust × arm in roll and pitch, half the hover reaction torque in yaw
      v = [i.thrust_cmd / 100, (nomM * i.roll_cmd) / 100, (nomM * i.pitch_cmd) / 100, (nomN * i.yaw_cmd) / 100];
    } else if (i.layout === 'helicopter') {
      names = ['Collective', 'Longitudinal cyclic', 'Lateral cyclic', 'Tail-rotor pedal']; axes = ['Thrust', 'Roll', 'Pitch', 'Yaw']; B = [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]; lo = [-1, -1, -1, -1]; hi = [1, 1, 1, 1]; w = [1, 1, 1, 1]; trim = [0, 0, 0, 0];
      v = [(i.thrust_cmd - 100) / 100, i.roll_cmd / 100, i.pitch_cmd / 100, i.yaw_cmd / 100];
    } else {
      names = ['Left elevator', 'Right elevator', 'Left aileron', 'Right aileron', 'Rudder', 'Left spoiler', 'Right spoiler']; axes = ['Roll', 'Pitch', 'Yaw'];
      // rows: roll, pitch, yaw moment per unit normalised deflection (nominal authority of each axis = 1 with primary effectors at full travel)
      B = [[-0.08, 0.08, -0.5, 0.5, 0.05, -0.35, 0.35], [0.5, 0.5, 0, 0, 0, 0.04, 0.04], [0, 0, 0.04, -0.04, 1, -0.12, 0.12]]; lo = [-1, -1, -1, -1, -1, 0, 0]; hi = [1, 1, 1, 1, 1, 1, 1]; w = [1, 1, 1, 1, 1, i.w_spoiler, i.w_spoiler]; trim = names.map(() => 0);
      v = [i.roll_cmd / 100, i.pitch_cmd / 100, i.yaw_cmd / 100];
    }
    const m = names.length, nom = allocate(B, v, w, lo, hi), fk = Math.round(i.failed) - 1, has = fk >= 0 && fk < m, lo2 = lo.slice(), hi2 = hi.slice();
    if (has) { const p = i.fail_mode.startsWith('jam') ? N.clamp(i.jam_pos, lo[fk], hi[fk]) : 0; lo2[fk] = hi2[fk] = p; }
    const fl = has ? allocate(B, v, w, lo2, hi2) : nom, vn = N.norm(v) || 1;
    // attainable roll–pitch moment envelope at the demanded thrust and zero yaw moment: largest scale factor per direction
    const iR = axes.indexOf('Roll'), iP = axes.indexOf('Pitch'), ths = N.linspace(0, TAU, 49), reach = (l, h) => ths.map((th) => { const dir = v.map((x, p) => (p === iR ? Math.cos(th) : p === iP ? Math.sin(th) : axes[p] === 'Thrust' ? x : 0)); let a = 0, b = 4; for (let k = 0; k < 22; k++) { const mid = 0.5 * (a + b), d = dir.map((x, p) => (p === iR || p === iP ? x * mid : x)); if (allocate(B, d, w, l, h).err < 1e-6) a = mid; else b = mid; } return a; });
    const rN = reach(lo, hi), rF = has ? reach(lo2, hi2) : rN, area = (r) => 0.5 * N.sum(r.slice(1).map((x, k) => x * r[k] * Math.sin(ths[k + 1] - ths[k]))), loss = 1 - area(rF) / (area(rN) || 1);
    const rk = (() => { let r = 0; const G = N.range(B.length, (p) => N.range(B.length, (q) => N.dot(B[p].map((x, k) => (has && k === fk ? 0 : x)), B[q].map((x, k) => (has && k === fk ? 0 : x))))), e = N.eigSym(G).values; for (const x of e) if (x > 1e-9 * N.amax(e)) r++; return r; })();
    if (nom.err > 1e-6 * vn) warnings.push('The demand exceeds the attainable set even without a failure: effectors saturate and the achieved moments fall short.');
    if (has && fl.err > 1e-6 * vn) warnings.push(`With ${names[fk]} failed the demand cannot be met: shortfall ${(100 * fl.err / vn).toFixed(0)}% of the demand.`);
    if (rk < B.length) warnings.push(`After the failure the effectors span only ${rk} of ${B.length} control axes: one axis is uncontrollable${i.layout === 'multirotor' && m < 6 ? ' (a quadrotor cannot keep yaw control with a rotor out)' : ''}.`);
    if (i.layout === 'helicopter') warnings.push('A conventional helicopter has no redundant effectors: each control maps to one axis, so any control failure removes that axis.');
    return {
      kpis: [
        { key: 'alloc_error_pct', label: 'Demand shortfall, no failure', value: (100 * nom.err) / vn, unit: '%', status: nom.err < 1e-6 * vn ? 'ok' : 'warn' }, { key: 'alloc_error_fail_pct', label: 'Demand shortfall with the failure', value: (100 * fl.err) / vn, unit: '%', status: fl.err < 1e-6 * vn ? 'ok' : 'bad' },
        { key: 'saturated', label: 'Saturated effectors, no failure', value: nom.nSat, unit: '-' }, { key: 'saturated_fail', label: 'Saturated or failed effectors after re-allocation', value: fl.nSat, unit: '-' },
        { key: 'envelope_loss_pct', label: 'Roll–pitch moment envelope lost to the failure', value: 100 * loss, unit: '%', status: loss < 0.5 ? 'ok' : 'warn' }, { key: 'control_rank', label: 'Controllable axes after the failure', value: rk, unit: `of ${B.length}`, status: rk === B.length ? 'ok' : 'bad' },
        { key: 'effort', label: 'Weighted control effort, no failure', value: Math.sqrt(N.sum(nom.u.map((x, k) => w[k] * (x - trim[k]) ** 2))), unit: '-' }, { key: 'effort_fail', label: 'Weighted control effort with the failure', value: Math.sqrt(N.sum(fl.u.map((x, k) => w[k] * (x - trim[k]) ** 2))), unit: '-' },
        { key: 'redundancy', label: 'Redundant effectors', value: m - B.length, unit: '-' },
      ],
      plots: [
        { type: 'bar', title: 'Effector commands', ylabel: 'Command [fraction of travel or of hover thrust]', categories: names, series: [{ name: 'No failure', y: nom.u }, { name: has ? `${names[fk]} failed` : 'No failure (repeat)', y: fl.u }] },
        { type: 'bar', title: 'Demanded and achieved control', ylabel: 'Normalised moment / thrust [-]', categories: axes, series: [{ name: 'Demand', y: v }, { name: 'Achieved', y: nom.got }, { name: 'Achieved with failure', y: fl.got }] },
        { type: 'line', title: 'Attainable roll–pitch moment envelope', xlabel: 'Roll moment [normalised]', ylabel: 'Pitch moment [normalised]', equalAspect: true, series: [{ name: 'All effectors', x: rN.map((r, k) => r * Math.cos(ths[k])), y: rN.map((r, k) => r * Math.sin(ths[k])) }, { name: 'After the failure', x: rF.map((r, k) => r * Math.cos(ths[k])), y: rF.map((r, k) => r * Math.sin(ths[k])) }, { name: 'Demand', x: [v[iR]], y: [v[iP]], style: 'points' }] },
      ],
      tables: [{ title: 'Control effectiveness matrix B (axes × effectors)', columns: ['Axis', ...names], rows: B.map((r, p) => [axes[p], ...r]) }],
      warnings, models: ['Weighted pseudo-inverse allocation u = W⁻¹Bᵀ(BW⁻¹Bᵀ)⁻¹v', 'Redistributed pseudo-inverse for saturation handling', 'Fault-tolerant re-allocation with the failed effector removed or frozen', 'Attainable moment set by directional search'],
      assumptions: ['Linear, decoupled effectiveness with normalised units; the aeroplane matrix is an illustrative transport layout, not data for a specific aircraft: replace it with identified effectiveness', 'The failure is postulated by the input, not predicted: its probability belongs to Suite 22', 'Static allocation: effector dynamics and rate limits are not considered', 'Redistribution is not guaranteed to find the true optimum on the boundary of the attainable set (a constrained QP would)'],
    };
  },
  verify() {
    const a = allocate([[1, 1]], [1], [1, 1], [-1, -1], [1, 1]), s = allocate([[1, 1]], [1], [1, 1], [-1, -1], [0.2, 1]), wv = allocate([[1, 1]], [1], [1, 3], [-9, -9], [9, 9]);
    const b = Object.fromEntries(allocation.inputs.map((f) => [f.key, f.default])), o = N.kv(allocation.run({ ...b, failed: 0 }));
    return [N.check('Minimum-norm split of a single demand over two equal effectors', a.u[0], 0.5, 1e-9, 'Pseudo-inverse'), N.check('Saturated effector is frozen and the remainder redistributed', s.u[1], 0.8, 1e-9, 'Redistributed pseudo-inverse'),
      N.check('Weights shift effort inversely: u1/u2 = w2/w1', wv.u[0] / wv.u[1], 3, 1e-9, 'Weighted least-norm solution'), N.check('Feasible demand is met exactly', o.alloc_error_pct + 1, 1, 1e-8, 'B·u = v')];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.control_rank < (i.layout === 'aeroplane surfaces' ? 3 : 4)) out.push(i.layout === 'helicopter'
      ? { severity: 'warn', title: 'Each helicopter control channel is a single path', detail: `With the postulated loss of one channel only ${o.control_rank} of 4 axes remain controllable: a conventional helicopter has no redundant effectors.`, action: 'Make the loss of a whole channel extremely improbable by design — duplicated hydraulic supplies and actuators, dual load paths in the control runs and swashplate — and show it in Suite 22.', basis: 'Rank of the control effectiveness matrix for the postulated failure' }
      : { severity: 'critical', title: 'Loss of control in one axis after the failure', detail: `With the postulated failure only ${o.control_rank} axes remain controllable.`, action: i.layout === 'multirotor' ? 'Use six or more rotors for single-rotor-out capability, or accept a yaw-spinning emergency landing mode for a quadrotor.' : 'Add a redundant effector for that axis or split the surface with independent actuators.', basis: 'Rank of the control effectiveness matrix for the postulated failure' });
    else if (o.alloc_error_fail_pct > 1) out.push({ severity: 'warn', title: 'Demand not achievable after the failure', detail: `Shortfall ${o.alloc_error_fail_pct.toFixed(0)}%; ${o.envelope_loss_pct.toFixed(0)}% of the roll–pitch envelope is lost.`, action: 'Restrict the manoeuvre envelope in the degraded mode and prioritise axes (pitch and roll before yaw).', basis: 'Attainable moment set' });
    else if (o.envelope_loss_pct > 30) out.push({ severity: 'advise', title: 'Reduced control envelope after the failure', detail: `${o.envelope_loss_pct.toFixed(0)}% of the roll–pitch envelope is lost.`, action: 'Feed the degraded authority into the failure-case handling assessment (Suite 22).', basis: 'Attainable moment set' });
    return out;
  },
};

// ---- 6. hover cascade for multirotors and helicopters -------------------------------------------
const hoverCascade = {
  id: 'hovercascade', title: 'Hover position, altitude and attitude cascade', fidelity: 'numerical',
  summary: 'Cascaded rate, attitude, velocity and position loops with an altitude loop for a multirotor or helicopter in hover, simulated at the controller sample rate with motor or actuator lag, tilt and thrust limits, gusts and sensor noise.',
  equations: ['Proportional–integral–derivative control equations', 'State-space equations', 'Actuator dynamic equations', 'Transfer function equations'],
  applicable: (c) => (c.rotor.R_m > 0 ? true : 'This analysis needs lifting rotors.'),
  inputs: [
    { key: 'dx', label: 'Position step', unit: 'm', default: 5, min: -200, max: 200, group: 'Commands' }, { key: 'dz', label: 'Altitude step', unit: 'm', default: 2, min: -100, max: 100, group: 'Commands' },
    { key: 'w_rate', label: 'Rate-loop bandwidth', unit: 'rad/s', default: 20, min: 0.2, max: 200, group: 'Design bandwidths', help: 'Keep below about 80% of the actuator or motor bandwidth 1/τ' }, { key: 'sep', label: 'Bandwidth separation between nested loops', unit: '-', default: 3, min: 2, max: 10, group: 'Design bandwidths' },
    { key: 'w_alt', label: 'Altitude-loop bandwidth', unit: 'rad/s', default: 1.5, min: 0.05, max: 10, group: 'Design bandwidths' },
    { key: 'tau_m', label: 'Motor / actuator time constant', unit: 's', default: 0.05, min: 0.002, max: 1, group: 'Actuation' }, { key: 'u_lim', label: 'Pitch-control authority', unit: 'fraction', default: 0.5, min: 0.02, max: 1, group: 'Actuation', help: 'Share of the control range available for pitch moment' },
    { key: 'thrust_ratio', label: 'Maximum / hover thrust', unit: '-', default: 2, min: 1.05, max: 5, group: 'Actuation' }, { key: 'tilt_lim_deg', label: 'Tilt command limit', unit: 'deg', default: 25, min: 2, max: 60, group: 'Actuation' },
    { key: 'Ts', label: 'Controller sample time', unit: 's', default: 0.004, min: 2e-4, max: 0.1, group: 'Digital implementation' },
    { key: 'gust', label: 'Horizontal gust RMS', unit: 'm/s', default: 1, min: 0, max: 15, group: 'Disturbances' }, { key: 'sigma_gyro_dps', label: 'Rate-gyro noise', unit: 'deg/s', default: 0.3, min: 0, max: 10, group: 'Disturbances' },
    { key: 'seed', label: 'Random seed', unit: '', default: 11, min: 1, max: 1e9, step: 1, discrete: true, group: 'Disturbances' }, { key: 't_end', label: 'Simulated time', unit: 's', default: 12, min: 2, max: 120, group: 'Numerics' },
    ...PLANT.filter((f) => ['Xu', 'Mu', 'Mq', 'Mc', 'Zw', 'Zc'].includes(f.key)).map((f) => ({ ...f, group: 'Hover plant', default: { Xu: -0.25, Mu: 0.4, Mq: -1.7, Mc: 125, Zw: -0.5, Zc: -9.81 }[f.key] }))],
  defaults: (c, up, d) => { const h = hoverEst(c), wr = Math.min(0.8 / h.tau_m, h.heli ? 6 : 30); return { Xu: h.Xu, Mu: h.Mu, Mq: h.Mq, Mc: h.Mc, Zw: h.Zw, Zc: h.Zc, tau_m: h.tau_m, w_rate: wr, w_alt: h.heli ? 0.8 : 1.5, Ts: h.heli ? 0.01 : 0.004, thrust_ratio: h.heli ? 1.3 : Math.max(1.1, (c.prop.T0_N * c.prop.n_eng) / d.W) || 2, u_lim: h.heli ? 0.15 : 0.5, tilt_lim_deg: h.heli ? 15 : 25, dx: h.heli ? 20 : 5, dz: h.heli ? 5 : 2, t_end: N.clamp(170 / wr, 12, 110), gust: Math.max(0.5, c.site.gust_ms || 1) }; },
  run(i) {
    const wA = i.w_rate / i.sep, wV = wA / i.sep, wX = wV / i.sep, Kq = i.w_rate / i.Mc, Kw = 3 * i.w_alt, Kwi = (Kw * i.w_alt) / 2, tiltL = rad(i.tilt_lim_deg), cmax = i.thrust_ratio - 1, vLim = (tiltL * G0) / wV, warnings = [];
    const Ts = i.Ts, nSub = Math.max(1, Math.ceil(Ts / Math.min(i.tau_m / 5, 0.002))), h = Ts / nSub, n = Math.round(i.t_end / Ts), rng = N.rng(i.seed), tg = 2;
    // states: u, q, θ, x, w (down +), z (up +), pitch actuator, collective actuator
    // thrust can change by −80% to +(ratio − 1) of the hover thrust, i.e. the collective by that share of g/|Zc| in its own unit (fraction of hover thrust for a multirotor, rad for a helicopter)
    const cUnit = G0 / Math.abs(i.Zc); let aUp = 0;
    let u = 0, q = 0, th = 0, x = 0, w = 0, z = 0, da = 0, dc = 0, iz = 0, gst = 0, sat = 0, was = false; const T = [], X = [], Z = [], TH = [], DA = [], DC = [], t0 = 1;
    for (let k = 0; k <= n; k++) {
      const t = k * Ts, xr = t >= t0 ? i.dx : 0, zr = t >= t0 ? i.dz : 0, qm = q + rad(i.sigma_gyro_dps) * N.randn(rng);
      const vC = N.clamp(wX * (xr - x), -vLim, vLim), aC = wV * (vC - u), thC = N.clamp(-aC / G0, -tiltL, tiltL), qC = wA * (thC - th), uU = Kq * (qC - qm), uC = N.clamp(uU, -i.u_lim, i.u_lim);
      // altitude: climb-rate command from the height error, PI on the climb-rate error; Zc < 0 so positive collective accelerates upward
      const wC = N.clamp(i.w_alt * (zr - z), -5, 5), ew = wC + w, cU = (Kw * ew + Kwi * iz) / Math.abs(i.Zc), cC = N.clamp(cU, -0.8 * cUnit, cmax * cUnit); if (cU === cC) iz += ew * Ts;
      const s = uU !== uC || cU !== cC; if (s && !was) sat++; was = s;
      T.push(t); X.push(x); Z.push(z); TH.push(th); DA.push(da); DC.push(dc);
      gst = gst * Math.exp(-Ts / tg) + i.gust * Math.sqrt(1 - Math.exp((-2 * Ts) / tg)) * N.randn(rng);
      for (let j = 0; j < nSub; j++) { da += ((uC - da) / i.tau_m) * h; dc += ((cC - dc) / i.tau_m) * h; const ur = u - gst, ud = i.Xu * ur - G0 * th, qd = i.Mu * ur + i.Mq * q + i.Mc * da, wd = i.Zw * w + i.Zc * dc; aUp = Math.max(aUp, -wd); u += ud * h; q += qd * h; th += q * h; x += u * h; w += wd * h; z += -w * h; }
      if (!fin(x) || Math.abs(th) > 10) { warnings.push('The simulation diverged: reduce the rate-loop bandwidth relative to the actuator lag or the sample time.'); break; }
    }
    const k0 = Math.round(t0 / Ts), sx = stepMetrics(T.slice(k0), X.slice(k0), i.dx || 1e-9), sz = stepMetrics(T.slice(k0), Z.slice(k0), i.dz || 1e-9), tail = Math.round(0.75 * T.length), rms = (a, r) => Math.sqrt(N.mean(a.slice(tail).map((v) => (v - r) ** 2)));
    const phm = deg(Math.atan2(1, i.w_rate * i.tau_m)) - deg(i.w_rate * Ts * 1.5), tauQ = 1 / (i.Mc * Kq - i.Mq);
    if (i.w_rate * i.tau_m > 1) warnings.push('The rate-loop bandwidth is above the actuator bandwidth: little phase margin remains.');
    if (i.Mu > 0 && Math.cbrt(G0 * i.Mu) > wA / 2) warnings.push('The attitude loop is not much faster than the unstable hover oscillation: actuation is too slow for tight position control (large rotors need collective-pitch rather than speed control).');
    if (sat) warnings.push(`Control saturation occurred ${sat} time(s) (pitch authority or thrust limit).`);
    if (N.amax(TH.map(Math.abs)) > rad(20)) warnings.push('Tilt exceeded 20°: the small-angle hover model underestimates the height loss and coupling.');
    return {
      kpis: [
        { key: 'pos_rise_s', label: 'Position step rise time', value: sx.rise, unit: 's' }, { key: 'pos_overshoot_pct', label: 'Position overshoot', value: sx.overshoot, unit: '%', status: sx.overshoot < 20 ? 'ok' : 'warn' }, { key: 'pos_settle_s', label: 'Position settling time', value: sx.settle - t0, unit: 's' },
        { key: 'alt_rise_s', label: 'Altitude step rise time', value: sz.rise, unit: 's' }, { key: 'alt_overshoot_pct', label: 'Altitude overshoot', value: sz.overshoot, unit: '%', status: sz.overshoot < 20 ? 'ok' : 'warn' },
        { key: 'pos_hold_rms_m', label: 'Position-hold RMS error in gusts', value: rms(X, i.dx), unit: 'm' }, { key: 'alt_hold_rms_m', label: 'Altitude-hold RMS error', value: rms(Z, i.dz), unit: 'm' },
        { key: 'tilt_max_deg', label: 'Largest tilt angle', value: deg(N.amax(TH.map(Math.abs))), unit: 'deg' }, { key: 'climb_accel_max_g', label: 'Largest upward acceleration', value: aUp / G0, unit: 'g', note: `Thrust limit allows ${(cmax).toFixed(2)} g` }, { key: 'saturation_events', label: 'Control saturation events', value: sat, unit: '-', status: sat ? 'warn' : 'ok' },
        { key: 'rate_gain', label: 'Rate-loop gain', value: Kq, unit: 'control per rad/s' }, { key: 'w_att_rads', label: 'Attitude-loop bandwidth', value: wA, unit: 'rad/s' }, { key: 'w_pos_rads', label: 'Position-loop bandwidth', value: wX, unit: 'rad/s' },
        { key: 'rate_loop_pm_deg', label: 'Rate-loop phase margin estimate', value: phm, unit: 'deg', status: phm > 45 ? 'ok' : 'warn', note: 'atan(1/(ω·τ)) minus sample-and-compute delay of 1.5·Ts' },
        { key: 'rate_tau_s', label: 'Ideal rate-loop time constant', value: tauQ, unit: 's' },
      ].filter((k) => fin(k.value)),
      plots: [
        { type: 'line', title: 'Position and altitude tracking', xlabel: 'Time [s]', ylabel: 'Displacement [m]', series: [{ name: 'Position', x: thin(T), y: thin(X) }, { name: 'Altitude', x: thin(T), y: thin(Z) }, { name: 'Position command', x: [0, t0, t0, i.t_end], y: [0, 0, i.dx, i.dx], style: 'dash' }, { name: 'Altitude command', x: [0, t0, t0, i.t_end], y: [0, 0, i.dz, i.dz], style: 'dash' }] },
        { type: 'line', title: 'Pitch attitude', xlabel: 'Time [s]', ylabel: 'Pitch attitude [deg]', series: [{ name: 'θ', x: thin(T), y: thin(TH.map(deg)) }] },
        { type: 'line', title: 'Control activity', xlabel: 'Time [s]', ylabel: 'Control [fraction]', series: [{ name: 'Pitch control', x: thin(T), y: thin(DA) }, { name: 'Collective / thrust change', x: thin(T), y: thin(DC) }] },
      ],
      warnings, models: ['Cascaded P–P–P–P position/velocity/attitude/rate loops with bandwidth separation', 'PI climb-rate loop for altitude', 'First-order motor or swashplate actuator lag', 'Hover speed–attitude plant (Xu, Mu, Mq) with seeded gusts and gyro noise, discrete controller update'],
      assumptions: ['One horizontal axis plus heave; small-angle thrust tilt', 'Thrust responds through a single first-order lag; no battery sag or rotor inflow dynamics', 'Ideal position and attitude estimates apart from gyro noise (estimation is in Suite 17)', 'Bandwidth separation, tilt, authority and thrust limits and the motor time constant are typical values, not data for a specific vehicle'],
    };
  },
  verify() {
    // rate loop alone with an ideal actuator: q → qc with time constant 1/(Mc·Kq − Mq)
    const Mc = 100, Mq = -2, wr = 10, Kq = wr / Mc, r = N.rk4((t, y) => [Mq * y[0] + Mc * Kq * (1 - y[0])], 0, [0], 1, 2000), tau = 1 / (Mc * Kq - Mq), yss = (Mc * Kq) / (Mc * Kq - Mq), k63 = r.y.findIndex((y) => y[0] >= 0.632121 * yss);
    const b = Object.fromEntries(hoverCascade.inputs.map((f) => [f.key, f.default])), o = N.kv(hoverCascade.run({ ...b, gust: 0, sigma_gyro_dps: 0, t_end: 25 }));
    const lim = N.kv(hoverCascade.run({ ...b, gust: 0, sigma_gyro_dps: 0, dx: 0, dz: 30, Zc: -60, Zw: 0, thrust_ratio: 1.2, t_end: 6 }));
    return [N.check('Rate-loop time constant 1/(Mc·Kq − Mq)', r.t[k63], tau, 0.01, 'First-order closed loop'),
      N.check('Climb acceleration is limited to (T/W − 1)·g whatever the collective unit', lim.climb_accel_max_g, 0.2, 1e-6, 'Thrust limit: a = (T_max − W)/m'), N.check('Position step is captured without steady-state error', o.pos_hold_rms_m / 5 + 1, 1, 0.01, 'Type-1 position loop'), N.check('Altitude step is captured', o.alt_hold_rms_m / 2 + 1, 1, 0.01, 'Integral action in the climb-rate loop')];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (o.rate_loop_pm_deg < 45) out.push({ severity: 'warn', title: 'Rate loop is too fast for the actuation', detail: `Estimated phase margin ${o.rate_loop_pm_deg.toFixed(0)}°.`, action: 'Lower the rate-loop bandwidth, use faster motors/ESCs or a shorter sample time.', basis: 'Actuator lag and sampling delay' });
    if (o.saturation_events > 0) out.push({ severity: 'advise', title: 'Thrust or moment saturation', detail: `${o.saturation_events} event(s) in the manoeuvre.`, action: 'Raise the thrust-to-weight ratio or limit the commanded acceleration; saturation in gusts is the usual cause of loss of control.', basis: 'Control authority' });
    return out;
  },
};

export default {
  id: 'control', n: 16,
  tagline: 'Do the control laws keep the aircraft stable with margin, track what is commanded, and stay in control when an actuator saturates or fails?',
  analyses: [pitchLoop, lqrAn, autopilot, sas, allocation, hoverCascade],
  consumes: [
    { from: 'flightdyn', keys: ['A_long', 'B_long', 'A_lat', 'B_lat', 'A_hover_long', 'B_hover_long', 'CLa_total', 'Cm_alpha', 'Cn_beta'], why: 'Linear aircraft models and derivatives for control design' },
    { from: 'cfd', keys: ['CD0'], why: 'Drag for the speed dynamics of the built-in model' },
  ],
  provides: [{ key: 'gm_dB', label: 'Gain margin', unit: 'dB' }, { key: 'pm_deg', label: 'Phase margin', unit: 'deg' }, { key: 'settling_s', label: 'Settling time', unit: 's' }, { key: 'overshoot_pct', label: 'Overshoot', unit: '%' }, { key: 'bandwidth_rads', label: 'Closed-loop bandwidth', unit: 'rad/s' }],
  handoff: [
    { model: 'Robust H-infinity and μ-synthesis controllers', why: 'Need iterative Riccati/LMI synthesis with weighting-function design and structured-uncertainty analysis', tool: 'Robust-control design toolbox' },
    { model: 'Non-linear model predictive control and non-linear dynamic inversion', why: 'Require an on-line constrained optimiser or a validated non-linear aerodynamic model across the envelope', tool: 'Optimal-control / MPC framework with the full flight-dynamics model' },
    { model: 'Adaptive control', why: 'Stability proofs and tuning depend on the full non-linear plant and failure scenarios', tool: 'Dedicated adaptive-control design and simulation environment' },
    { model: 'Fault detection, isolation and redundancy management logic', why: 'Discrete monitoring logic and voting must be modelled as a hybrid system with the real architecture', tool: 'Model-based design tool with state machines; safety assessment in Suite 22' },
    { model: 'Software-, processor- and hardware-in-the-loop testing', why: 'Needs the real flight-control computer, actuators and a real-time simulator', tool: 'Iron-bird / HIL rig' },
    { model: 'Actuator hinge-moment and power sizing', why: 'Hinge moments and hydraulic or electrical power are sized in the systems suites', tool: 'Suites 18 and 19' },
  ],
};
