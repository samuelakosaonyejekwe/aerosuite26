// Suite 4 — Flight Dynamics, Stability and Control.
// Preliminary stability-derivative build-up from geometry, longitudinal trim, linear longitudinal and
// lateral-directional modes with handling-qualities checks, a quaternion six-degree-of-freedom simulation,
// turbulence response (Dryden / von Kármán) and hover dynamics of helicopters and multirotors.

import * as N from '../core/numerics.js';
import { isa, G0 } from '../core/atmosphere.js';

const TAU = 2 * Math.PI, fin = Number.isFinite, rad = N.rad, deg = N.deg;
const isPropType = (t) => t !== 'turbofan' && t !== 'turbojet';
const thin = (a, n = 300) => { const s = Math.max(1, Math.ceil(a.length / n)); return a.filter((_, k) => k % s === 0); };

// ---- fixed-wing aerodynamic model -------------------------------------------------------------
const FW = [
  { key: 'mass_kg', label: 'Aircraft mass', unit: 'kg', default: 70000, min: 0.01, group: 'Mass' },
  { key: 'Ixx', label: 'Roll inertia Ixx', unit: 'kg·m²', default: 1.2e6, min: 1e-6, group: 'Mass' }, { key: 'Iyy', label: 'Pitch inertia Iyy', unit: 'kg·m²', default: 3.8e6, min: 1e-6, group: 'Mass' },
  { key: 'Izz', label: 'Yaw inertia Izz', unit: 'kg·m²', default: 4.9e6, min: 1e-6, group: 'Mass' }, { key: 'Ixz', label: 'Product of inertia Ixz', unit: 'kg·m²', default: 1e5, group: 'Mass' },
  { key: 'cg', label: 'Centre of gravity', unit: 'MAC', default: 0.25, min: -0.5, max: 1, group: 'Mass', help: 'Fraction of the mean aerodynamic chord aft of its leading edge' },
  { key: 'S', label: 'Wing area', unit: 'm²', default: 122, min: 0.01, group: 'Wing' }, { key: 'b', label: 'Wing span', unit: 'm', default: 34, min: 0.05, group: 'Wing' },
  { key: 'mac', label: 'Mean aerodynamic chord', unit: 'm', default: 4.2, min: 0.01, group: 'Wing' }, { key: 'taper', label: 'Taper ratio', unit: '-', default: 0.3, min: 0.05, max: 1, group: 'Wing' },
  { key: 'sweep_deg', label: 'Quarter-chord sweep', unit: 'deg', default: 25, min: -30, max: 60, group: 'Wing' }, { key: 'dihedral_deg', label: 'Dihedral', unit: 'deg', default: 5, min: -10, max: 15, group: 'Wing' },
  { key: 'alpha0_deg', label: 'Body angle of attack at zero lift', unit: 'deg', default: -4, min: -15, max: 5, group: 'Wing', help: 'Aerofoil zero-lift angle minus the wing incidence' },
  { key: 'CLa_w', label: 'Wing lift-curve slope (0 = estimate)', unit: '1/rad', default: 0, min: 0, max: 8, group: 'Wing', help: 'From the aerodynamics suite when available' },
  { key: 'St', label: 'Horizontal tail area', unit: 'm²', default: 31, min: 0, group: 'Tail' }, { key: 'bt', label: 'Horizontal tail span', unit: 'm', default: 12.4, min: 0.01, group: 'Tail' },
  { key: 'lt', label: 'Horizontal tail arm (CG to tail a.c.)', unit: 'm', default: 17.7, min: 0, group: 'Tail' },
  { key: 'Sv', label: 'Vertical tail area', unit: 'm²', default: 21.5, min: 0, group: 'Tail' }, { key: 'bv', label: 'Vertical tail height', unit: 'm', default: 5.9, min: 0.01, group: 'Tail' },
  { key: 'lv', label: 'Vertical tail arm', unit: 'm', default: 16.2, min: 0, group: 'Tail' },
  { key: 'fus_len', label: 'Fuselage length', unit: 'm', default: 37.6, min: 0.01, group: 'Fuselage' }, { key: 'fus_dia', label: 'Fuselage diameter', unit: 'm', default: 3.95, min: 0.001, group: 'Fuselage' },
  { key: 'CD0', label: 'Zero-lift drag coefficient', unit: '-', default: 0.02, min: 1e-7, group: 'Aerodynamics' }, { key: 'e', label: 'Oswald efficiency', unit: '-', default: 0.8, min: 0.3, group: 'Aerodynamics' },
  { key: 'Cm0', label: 'Zero-lift pitching moment', unit: '-', default: -0.05, min: -0.5, max: 0.5, group: 'Aerodynamics', help: 'Whole aircraft, controls neutral' },
  { key: 'CLmax', label: 'Maximum lift coefficient', unit: '-', default: 1.5, min: 0.3, group: 'Aerodynamics' },
  { key: 'Se_Sh', label: 'Elevator / tail chord ratio', unit: '-', default: 0.3, min: 0.05, max: 1, group: 'Controls' }, { key: 'Sa_S', label: 'Aileron area / wing area', unit: '-', default: 0.06, min: 0.005, max: 0.3, group: 'Controls' },
  { key: 'Sr_Sv', label: 'Rudder / fin chord ratio', unit: '-', default: 0.3, min: 0.05, max: 1, group: 'Controls' },
  { key: 'ptype', label: 'Powerplant type', type: 'select', options: ['turbofan', 'turbojet', 'turboprop', 'turboshaft', 'piston', 'electric'], default: 'turbofan', group: 'Propulsion' },
  { key: 'T_static', label: 'Total static thrust', unit: 'N', default: 240000, min: 0, group: 'Propulsion' }, { key: 'P_total', label: 'Total shaft power', unit: 'W', default: 0, min: 0, group: 'Propulsion' },
  { key: 'eta_prop', label: 'Propeller efficiency', unit: '-', default: 0.82, min: 0.1, max: 0.95, group: 'Propulsion' },
  { key: 'V', label: 'True airspeed', unit: 'm/s', default: 231, min: 1, group: 'Flight condition' }, { key: 'alt_m', label: 'Altitude', unit: 'm', default: 10668, min: -500, max: 25000, group: 'Flight condition' },
  { key: 'compress', label: 'Mach effect on lift (Prandtl–Glauert)', type: 'bool', default: true, group: 'Aerodynamics' },
  { key: 'k_Cma', label: 'Calibration factor on Cmα', unit: '-', default: 1, min: 0.1, max: 3, group: 'Calibration' }, { key: 'k_Cmq', label: 'Calibration factor on Cmq', unit: '-', default: 1, min: 0.1, max: 3, group: 'Calibration' },
  { key: 'k_Clp', label: 'Calibration factor on Clp', unit: '-', default: 1, min: 0.1, max: 3, group: 'Calibration' }, { key: 'k_Cnb', label: 'Calibration factor on Cnβ', unit: '-', default: 1, min: 0.1, max: 3, group: 'Calibration' },
];
const fwDefaults = (c, up, d) => {
  if (!(c.wing.S_m2 > 0)) return {};
  const camber = Number(String(c.wing.airfoil)[0]) || 0, m = c.mass.mtow_kg, k = (v, est) => (v > 0 ? v : est);
  return {
    mass_kg: m, Ixx: k(c.mass.Ixx, m * (0.25 * c.wing.b_m / 2) ** 2), Iyy: k(c.mass.Iyy, m * (0.36 * c.fuselage.len_m / 2) ** 2), Izz: k(c.mass.Izz, m * (0.4 * (c.wing.b_m + c.fuselage.len_m) / 4) ** 2), Ixz: c.mass.Ixz,
    cg: c.mass.cg_pct_mac / 100, S: c.wing.S_m2, b: c.wing.b_m, mac: d.mac, taper: c.wing.taper, sweep_deg: c.wing.sweep_deg, dihedral_deg: c.wing.dihedral_deg, alpha0_deg: -1.04 * camber - c.wing.incidence_deg,
    CLa_w: undefined, St: c.htail.S_m2, bt: c.htail.b_m || undefined, lt: c.htail.arm_m, Sv: c.vtail.S_m2, bv: c.vtail.b_m || undefined, lv: c.vtail.arm_m, fus_len: c.fuselage.len_m || undefined, fus_dia: c.fuselage.dia_m || undefined,
    CD0: up.cfd?.CD0 ?? c.aero.CD0, e: up.cfd?.e_oswald ?? c.aero.e, Cm0: c.aero.Cm0, CLmax: up.cfd?.CLmax ?? c.aero.CLmax_clean, Se_Sh: c.controls.Se_Sh, Sa_S: c.controls.Sa_S, Sr_Sv: c.controls.Sr_Sv,
    ptype: c.prop.type, T_static: up.propulsion?.thrust_static_N ? up.propulsion.thrust_static_N * c.prop.n_eng : d.T_total, P_total: d.P_total, eta_prop: up.propeller?.eta_prop ?? c.prop.eta_prop, V: c.flight.V_ms, alt_m: c.atm.alt_m,
  };
};
const hasWing = (c) => (c.wing.S_m2 > 0 ? true : 'This analysis needs a wing and tail; use the hover-dynamics and six-degree-of-freedom analyses for rotorcraft and multirotors.');
const tauFlap = (E) => { const th = Math.acos(N.clamp(2 * E - 1, -1, 1)); return 0.85 * (1 - (th - Math.sin(th)) / Math.PI); }; // thin-aerofoil effectiveness with a 15% viscous loss
const munk = (f) => N.interp1([2, 4, 6, 8, 10, 12, 20], [0.55, 0.77, 0.86, 0.91, 0.94, 0.955, 0.98], f); // apparent-mass factor k2 − k1 of a slender body
/** Thrust available at full throttle [N]. */
const thrustMax = (i, V, sigma) => (isPropType(i.ptype) && i.P_total > 0 ? Math.min(i.T_static > 0 ? i.T_static * sigma ** 0.7 : Infinity, (i.eta_prop * i.P_total * (i.ptype === 'electric' ? 1 : sigma ** 0.7)) / Math.max(V, 1)) : i.T_static * sigma ** 0.7);

/** Non-dimensional stability and control derivatives (per rad) by a DATCOM-style empirical build-up. */
export function derivs(i, V = i.V) {
  const at = isa(i.alt_m), M = V / at.a, be2 = i.compress ? Math.max(1 - M * M, 0.19) : 1, q = 0.5 * at.rho * V * V, AR = (i.b * i.b) / i.S, sw = rad(i.sweep_deg), lam = i.taper, c = i.mac;
  const slope = (A, s) => (TAU * A) / (2 + Math.sqrt(4 + A * A * be2 * (1 + Math.tan(s) ** 2 / be2)));
  const aw = i.CLa_w > 0 ? i.CLa_w : slope(AR, sw), at_ = slope((i.bt * i.bt) / Math.max(i.St, 1e-9), sw), av = slope((1.55 * i.bv * i.bv) / Math.max(i.Sv, 1e-9), sw), eta = 0.9;
  const eps = Math.min(0.7, (2 * aw) / (Math.PI * AR)), VH = (i.St * i.lt) / (i.S * c), k21 = munk(i.fus_len / i.fus_dia), vol = 0.7 * (Math.PI / 4) * i.fus_dia ** 2 * i.fus_len;
  const CLa = aw + eta * at_ * (i.St / i.S) * (1 - eps), Cma_f = (2 * k21 * vol) / (i.S * c), Cma = i.k_Cma * (aw * (i.cg - 0.25) + Cma_f - eta * at_ * VH * (1 - eps)), hn = i.cg - Cma / CLa;
  const CLq = 2 * eta * at_ * VH, Cmq = -2.2 * i.k_Cmq * eta * at_ * VH * (i.lt / c), CLad = CLq * eps, Cmad = -2 * eta * at_ * VH * (i.lt / c) * eps, te = tauFlap(i.Se_Sh), tr = tauFlap(i.Sr_Sv);
  const CL = (i.mass_kg * G0) / (q * i.S), kI = 1 / (Math.PI * AR * i.e), CD = i.CD0 + kI * CL * CL, zv = 0.5 * i.bv + 0.25 * i.fus_dia, gam = rad(i.dihedral_deg), ybar = (1 + 2 * lam) / (3 * (1 + lam));
  const CYb_v = -av * (i.Sv / i.S), CYb = CYb_v - (2 * (Math.PI / 4) * i.fus_dia ** 2) / i.S, Cnb = i.k_Cnb * (av * (i.Sv * i.lv) / (i.S * i.b) - (2 * k21 * vol) / (i.S * i.b));
  const Clb = 0.75 * (-(aw * gam * ybar) / 2 - CL * Math.tan(sw) * ybar) + CYb_v * (zv / i.b), // strip theory with a 0.75 finite-span relief factor
     Clp = (-i.k_Clp * aw * (1 + 3 * lam)) / (12 * (1 + lam)), Cnp = -CL / 8;
  const Clr = CL / 4 - 2 * (i.lv / i.b) * (zv / i.b) * CYb_v, Cnr = 2 * CYb_v * (i.lv / i.b) ** 2 - i.CD0 / 4, CYr = -2 * CYb_v * (i.lv / i.b);
  // aileron over 60–95% semi-span: Clδa = 2·aw·τ/(S·b)·∫c·y·dy; the aileron chord ratio is its area over the wing area of that strip
  const cr = (2 * i.S) / (i.b * (1 + lam)), y1 = 0.3 * i.b, y2 = 0.475 * i.b, icy = (y) => cr * (y * y / 2 - ((1 - lam) * 2 * y ** 3) / (3 * i.b)), ic = (y) => cr * (y - ((1 - lam) * y * y) / i.b);
  const ta = tauFlap(N.clamp((i.Sa_S * i.S) / (2 * (ic(y2) - ic(y1))), 0.05, 0.6)), Clda = ((2 * aw * ta) / (i.S * i.b)) * (icy(y2) - icy(y1));
  return {
    M, q, AR, kI, CL, CD, CLa, Cma, CLq, Cmq, CLad, Cmad, CLde: eta * at_ * (i.St / i.S) * te, Cmde: -eta * at_ * VH * te, CDa: 2 * kI * CL * CLa, CLu: i.compress ? (M * M / (1 - Math.min(M * M, 0.81))) * CL : 0,
    hn, sm: hn - i.cg, hm: hn - (Cmq * at.rho * i.S * c) / (4 * i.mass_kg), aw, at: at_, av, eps, VH, Cma_f, CYb, Cnb, Clb, Clp, Cnp, Clr, Cnr, CYr, Clda, Cnda: -0.2 * CL * Clda, CYdr: av * tr * (i.Sv / i.S), Cndr: -av * tr * (i.Sv * i.lv) / (i.S * i.b), Cldr: av * tr * (i.Sv / i.S) * (zv / i.b), atm: at,
  };
}
/** The same aircraft with the CG at h [MAC]: the tail arms are measured from the CG, so they change with it and the neutral point stays put. */
const withCg = (i, h) => ({ ...i, cg: h, lt: Math.max(i.lt - (h - i.cg) * i.mac, 0), lv: Math.max(i.lv - (h - i.cg) * i.mac, 0) });
/** Level-flight trim from the linear lift and moment equations: α (body), δe, thrust, throttle. */
function trimAt(i, V) {
  const D = derivs(i, V), det = D.CLa * D.Cmde - D.CLde * D.Cma, a0 = rad(i.alpha0_deg);
  const da = (D.CL * D.Cmde + D.CLde * i.Cm0) / det, de = (-D.CLa * i.Cm0 - D.Cma * D.CL) / det, T = D.q * i.S * D.CD, Tmax = thrustMax(i, V, D.atm.sigma);
  return { D, alpha: a0 + da, de, T, throttle: Tmax > 0 ? T / Tmax : NaN, Tmax };
}
/** Longitudinal state matrices, x = [u, w, q, θ], inputs [δe, throttle], stability axes. */
function longSS(i) {
  const t = trimAt(i, i.V), D = t.D, u0 = i.V, m = i.mass_kg, QS = D.q * i.S, c = i.mac;
  const Xu = (-2 * D.CD * QS) / (m * u0) - (isPropType(i.ptype) ? t.T / (m * u0) : 0), Xw = (-(D.CDa - D.CL) * QS) / (m * u0), Zu = (-(D.CLu + 2 * D.CL) * QS) / (m * u0), Zw = (-(D.CLa + D.CD) * QS) / (m * u0);
  const Zq = (-D.CLq * (c / (2 * u0)) * QS) / m, Mw = (D.Cma * QS * c) / (u0 * i.Iyy), Mwd = (D.Cmad * (c / (2 * u0)) * QS * c) / (u0 * i.Iyy), Mq = (D.Cmq * (c / (2 * u0)) * QS * c) / i.Iyy;
  const Zde = (-D.CLde * QS) / m, Mde = (D.Cmde * QS * c) / i.Iyy, XdT = t.Tmax / m;
  const A = [[Xu, Xw, 0, -G0], [Zu, Zw, u0 + Zq, 0], [Mwd * Zu, Mw + Mwd * Zw, Mq + Mwd * (u0 + Zq), 0], [0, 0, 1, 0]], B = [[0, XdT], [Zde, 0], [Mde + Mwd * Zde, 0], [0, 0]];
  return { A, B, t, D, d: { Xu, Xw, Zu, Zw, Zq, Mw, Mwd, Mq, Zde, Mde, XdT } };
}
/** Lateral-directional state matrices, x = [β, p, r, φ], inputs [δa, δr]. */
function latSS(i) {
  const D = derivs(i), u0 = i.V, QS = D.q * i.S, b = i.b, m = i.mass_kg, k = b / (2 * u0), G = 1 - (i.Ixz * i.Ixz) / (i.Ixx * i.Izz);
  const L = (C) => (QS * b * C) / i.Ixx, Nn = (C) => (QS * b * C) / i.Izz, st = (l, n) => [(l + (i.Ixz / i.Ixx) * n) / G, (n + (i.Ixz / i.Izz) * l) / G];
  const [Lb, Nb] = st(L(D.Clb), Nn(D.Cnb)), [Lp, Np] = st(L(D.Clp) * k, Nn(D.Cnp) * k), [Lr, Nr] = st(L(D.Clr) * k, Nn(D.Cnr) * k), [Lda, Nda] = st(L(D.Clda), Nn(D.Cnda)), [Ldr, Ndr] = st(L(D.Cldr), Nn(D.Cndr));
  const Yb = (QS * D.CYb) / m, Yr = (QS * D.CYr * k) / m, Ydr = (QS * D.CYdr) / m;
  const A = [[Yb / u0, 0, -(1 - Yr / u0), G0 / u0], [Lb, Lp, Lr, 0], [Nb, Np, Nr, 0], [0, 1, 0, 0]], B = [[0, Ydr / u0], [Lda, Ldr], [Nda, Ndr], [0, 0]];
  return { A, B, D, d: { Yb, Yr, Lb, Lp, Lr, Nb, Np, Nr, Lda, Nda, Ldr, Ndr } };
}
const modeOf = (l) => ({ re: l[0], im: l[1], wn: Math.hypot(l[0], l[1]), zeta: -l[0] / (Math.hypot(l[0], l[1]) || 1) });
/** Linear response ẋ = A·x + B·u(t) by RK4; returns {t, x}. */
function linSim(A, B, u, tEnd, n, x0) { const r = N.rk4((t, x) => N.vadd(N.matvec(A, x), N.matvec(B, u(t))), 0, x0 || A.map(() => 0), tEnd, n); return { t: r.t, x: r.y }; }
/** RK4 step count that keeps |λ|·h ≤ 0.5 for the fastest root (stability limit 2.78), and at least n. */
const simSteps = (modes, tEnd, n) => Math.ceil(N.clamp(2 * tEnd * N.amax(modes.map((m) => m.wn)), n, 30000));
/** Mode shape [p, r, φ] per unit sideslip for the eigenvalue l = {re, im}: rows 2–4 of (A − λI)·v = 0 with β = 1. */
function modeShape(A, l) {
  const lam = [l.re, l.im], M = [1, 2, 3].map((r) => [1, 2, 3].map((c) => (r === c ? N.C.sub([A[r][c], 0], lam) : [A[r][c], 0])));
  try { return N.csolve(M, [1, 2, 3].map((r) => [-A[r][0], 0])); } catch { return [1, 2, 3].map(() => [NaN, 0]); }
}
const doublet = (amp, t0, w) => (t) => (t >= t0 && t < t0 + w ? amp : t >= t0 + w && t < t0 + 2 * w ? -amp : 0);

// ---- 1. derivatives -------------------------------------------------------------------------
const derivAn = {
  id: 'derivs', title: 'Stability and control derivatives, neutral point and static margin', fidelity: 'analytical',
  summary: 'A preliminary empirical build-up of the longitudinal and lateral-directional derivatives from geometry, with the neutral point, manoeuvre point and static margin.',
  equations: ['Stability derivative equations', 'Aerodynamic force and moment equations', 'Aircraft trim equations'],
  applicable: hasWing, inputs: FW, defaults: fwDefaults,
  run(i) {
    const D = derivs(i), cgs = N.linspace(i.cg - 0.2, i.cg + 0.35, 30), warnings = [];
    if (D.sm < 0.05) warnings.push(D.sm < 0 ? 'The centre of gravity is aft of the neutral point: the aircraft is statically unstable in pitch.' : 'Static margin is below 5% MAC.');
    if (D.Cnb <= 0) warnings.push('Cnβ is not positive: the aircraft is directionally unstable; the fin is too small for this fuselage.');
    if (D.M > 0.7) warnings.push('Above Mach 0.7 the Prandtl–Glauert lift-slope correction and constant downwash are increasingly inaccurate.');
    const rows = [['CLα', D.CLa, 'Lift-curve slope'], ['Cmα', D.Cma, 'Pitch stiffness'], ['CLq', D.CLq, ''], ['Cmq', D.Cmq, 'Pitch damping'], ['CLα̇', D.CLad, ''], ['Cmα̇', D.Cmad, ''], ['CLδe', D.CLde, ''], ['Cmδe', D.Cmde, 'Elevator power'],
      ['CYβ', D.CYb, ''], ['Clβ', D.Clb, 'Dihedral effect'], ['Cnβ', D.Cnb, 'Weathercock stability'], ['Clp', D.Clp, 'Roll damping'], ['Cnp', D.Cnp, ''], ['Clr', D.Clr, ''], ['Cnr', D.Cnr, 'Yaw damping'], ['CYr', D.CYr, ''],
      ['Clδa', D.Clda, 'Aileron power'], ['Cnδa', D.Cnda, 'Adverse yaw'], ['CYδr', D.CYdr, ''], ['Clδr', D.Cldr, ''], ['Cnδr', D.Cndr, 'Rudder power']];
    return {
      kpis: [
        { key: 'static_margin', label: 'Static margin', value: D.sm, unit: 'MAC', status: D.sm > 0.05 && D.sm < 0.4 ? 'ok' : D.sm > 0 ? 'warn' : 'bad', note: 'Neutral point minus CG; 5–25% MAC is usual for manually flown aircraft' },
        { key: 'neutral_point', label: 'Stick-fixed neutral point', value: D.hn, unit: 'MAC' }, { key: 'manoeuvre_point', label: 'Stick-fixed manoeuvre point', value: D.hm, unit: 'MAC' },
        { key: 'CLa', label: 'Lift-curve slope CLα', value: D.CLa, unit: '1/rad' }, { key: 'Cma', label: 'Pitch stiffness Cmα', value: D.Cma, unit: '1/rad', status: D.Cma < 0 ? 'ok' : 'bad' },
        { key: 'Cmq', label: 'Pitch damping Cmq', value: D.Cmq, unit: '1/rad' }, { key: 'Cmde', label: 'Elevator power Cmδe', value: D.Cmde, unit: '1/rad' },
        { key: 'Clb', label: 'Dihedral effect Clβ', value: D.Clb, unit: '1/rad', status: D.Clb < 0 ? 'ok' : 'warn' }, { key: 'Cnb', label: 'Weathercock stability Cnβ', value: D.Cnb, unit: '1/rad', status: D.Cnb > 0 ? 'ok' : 'bad' },
        { key: 'Clp', label: 'Roll damping Clp', value: D.Clp, unit: '1/rad' }, { key: 'Cnr', label: 'Yaw damping Cnr', value: D.Cnr, unit: '1/rad' },
        { key: 'tail_volume', label: 'Horizontal tail volume coefficient', value: D.VH, unit: '-' }, { key: 'downwash_grad', label: 'Downwash gradient dε/dα', value: D.eps, unit: '-' },
      ],
      plots: [
        { type: 'line', title: 'Pitch stiffness against CG position', xlabel: 'CG position [MAC]', ylabel: 'Cmα [1/rad]', series: [{ name: 'Cmα', x: cgs, y: cgs.map((h) => derivs(withCg(i, h)).Cma) }], annotations: [{ y: 0, label: 'Neutral' }, { x: i.cg, label: 'CG' }, { x: D.hn, label: 'Neutral point' }] },
        { type: 'bar', title: 'Neutral-point build-up', ylabel: 'Contribution [MAC]', categories: ['Wing a.c.', 'Fuselage', 'Horizontal tail', 'Neutral point'], series: [{ name: 'Position / shift', y: [0.25, -D.Cma_f / D.CLa, D.hn - 0.25 + D.Cma_f / D.CLa, D.hn] }] },
      ],
      tables: [{ title: 'Derivatives (per radian)', columns: ['Derivative', 'Value', 'Meaning'], rows }],
      outputs: { CLa_total: D.CLa, Cm_alpha: D.Cma, Cl_delta_a: D.Clda, Cn_beta: D.Cnb },
      warnings, models: ['Empirical component build-up (DATCOM-style, preliminary)', 'Helmbold–Prandtl–Glauert lift-curve slope', 'Munk slender-body fuselage moment', 'Thin-aerofoil control effectiveness with a viscous factor', 'Strip-theory roll damping and dihedral effect'],
      assumptions: ['Rigid aircraft, attached flow, linear aerodynamics', 'Downwash gradient 2·CLα/(π·AR) limited to 0.7; tail efficiency 0.9', 'Power, nacelle, flap and ground effects are not included', 'Ailerons span 60–95% of the semi-span; their chord ratio follows from the aileron area', 'Empirical factors (tail efficiency 0.9, 15% viscous loss of control effectiveness, 0.75 dihedral relief, Cnp = −CL/8, adverse yaw −0.2·CL·Clδa, fuselage volume 0.7 × its enclosing cylinder) are typical textbook values, not data for a specific aircraft', 'Expect ±10–20% on primary and ±50% on cross derivatives until calibrated'],
    };
  },
  calibration: { params: [{ key: 'k_Cma', min: 0.3, max: 2 }, { key: 'k_Cmq', min: 0.3, max: 2 }], sweep: 'cg', target: 'Cma', note: 'Wind-tunnel or flight-identified Cmα at several CG positions' },
  verify() {
    const b = Object.fromEntries(FW.map((f) => [f.key, f.default])), D = derivs({ ...b, V: 60, alt_m: 0, compress: false, sweep_deg: 0, b: 1e4, S: 1e4, mac: 1, St: 1e-9, Sv: 1e-9, fus_dia: 1e-6, bt: 1, bv: 1 }), D2 = derivs(withCg(b, 0.4)), D3 = derivs(b), rect = { ...b, taper: 1, CLa_w: 5 };
    return [N.check('Lift-curve slope tends to 2π at very high aspect ratio', D.CLa, TAU, 2e-3, 'Thin-aerofoil theory'),
      N.check('Cmα shifts by CLα·Δh with CG movement', D2.Cma - D3.Cma, D3.CLa * 0.15, 1e-9, 'dCm/dCL = h − hn, tail arm measured from the CG'),
      N.check('Neutral point does not move with the CG', D2.hn, D3.hn, 1e-12, 'The neutral point is a property of the airframe'),
      N.check('Rectangular-wing strip roll damping −a/6', derivs(rect).Clp, -5 / 6, 1e-9, 'Strip theory: −(a/12)(1+3λ)/(1+λ)'),
      N.check('Aileron power of a rectangular wing follows the aileron area', derivs({ ...rect, Sa_S: 0.35 * 0.4 }).Clda, 5 * tauFlap(0.4) * (0.475 ** 2 - 0.3 ** 2), 1e-12, 'Strip theory, ailerons over 60–95% semi-span with a 40% chord ratio: a·τ·(η2² − η1²)/4')];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.static_margin < 0.05) out.push({ severity: o.static_margin < 0 ? 'critical' : 'warn', title: 'Insufficient static margin', detail: `Static margin ${(100 * o.static_margin).toFixed(1)}% MAC with the CG at ${(100 * i.cg).toFixed(0)}%.`, action: 'Move the CG forward, enlarge the tail or lengthen the tail arm; a relaxed-stability design needs the augmentation of Suite 16.', basis: 'Stick-fixed static stability' });
    if (o.static_margin > 0.35) out.push({ severity: 'advise', title: 'Large static margin', detail: `${(100 * o.static_margin).toFixed(0)}% MAC costs trim drag and elevator authority.`, action: 'An aft CG shift (fuel or payload) reduces tail download, trim drag and fuel burn.', basis: 'Trim drag' });
    if (o.Cnb <= 0.03) out.push({ severity: 'warn', title: 'Weak directional stability', detail: `Cnβ = ${o.Cnb.toFixed(3)} per rad.`, action: 'Increase fin area or arm, or add a dorsal/ventral fin.', basis: 'Rule of thumb: a minimum Cnβ of about 0.03–0.06 per rad is typical (not a sourced requirement)' });
    return out;
  },
};

// ---- 2. trim --------------------------------------------------------------------------------
const trimAn = {
  id: 'trim', title: 'Longitudinal trim across the speed range', fidelity: 'analytical',
  summary: 'Angle of attack, elevator angle and thrust for steady level flight from stall to maximum speed, with the elevator angle per g and the CG range it allows.',
  equations: ['Aircraft trim equations', 'Aerodynamic force and moment equations', 'Flight-path kinematic equations'],
  applicable: hasWing,
  inputs: [...FW, { key: 'de_max_deg', label: 'Elevator travel', unit: 'deg', default: 25, min: 5, max: 40, group: 'Controls' }, { key: 'n_limit', label: 'Limit load factor', unit: 'g', default: 2.5, min: 1.5, max: 9, group: 'Controls' }],
  defaults: (c, up, d) => ({ ...fwDefaults(c, up, d), de_max_deg: c.controls.de_max_deg, n_limit: c.aero.n_pos }),
  run(i) {
    const at = isa(i.alt_m), Vs = Math.sqrt((2 * i.mass_kg * G0) / (at.rho * i.S * i.CLmax)), Vtop = Math.max(1.2 * Vs, Math.min(Math.max(1.6 * i.V, 2.5 * Vs), i.compress ? 0.9 * at.a : Infinity)), Vs_ = N.linspace(1.05 * Vs, Vtop, 40), T = Vs_.map((V) => trimAt(i, V)), t0 = trimAt(i, i.V), D = t0.D;
    // steady pull-up: elevator per g
    const dn = 1, qh = (dn * G0 * i.mac) / (2 * i.V * i.V), sol = N.solve([[D.CLa, D.CLde], [D.Cma, D.Cmde]], [dn * D.CL - D.CLq * qh, -D.Cmq * qh]), dePerG = sol[1];
    const deLim = rad(i.de_max_deg), warnings = [];
    // forward CG limit: elevator needed to trim at CLmax, aft limit: 5% static margin
    const fwd = N.findRoot((h) => { const j = withCg(i, h), Dj = derivs(j, 1.05 * Vs), det = Dj.CLa * Dj.Cmde - Dj.CLde * Dj.Cma; return (-Dj.CLa * i.Cm0 - Dj.Cma * Dj.CL) / det + 0.9 * deLim; }, D.hn - 1.5, D.hn, 60);
    if (Math.abs(t0.de) > 0.6 * deLim) warnings.push('Trim at this condition uses more than 60% of elevator travel; a trimmable stabiliser would normally carry this.');
    if (t0.throttle > 1) warnings.push(`Thrust required is ${(100 * t0.throttle).toFixed(0)}% of the estimated thrust available: level flight ${t0.throttle > 1.1 ? 'cannot be' : 'is at best marginally'} sustained at this speed, altitude and mass.`);
    if (i.V < 1.1 * Vs) warnings.push('The flight speed is within 10% of the stall speed; linear aerodynamics underestimates the angle of attack.');
    const nReach = Math.min(i.n_limit, i.CLmax / D.CL); // below the manoeuvre speed the wing stalls before the limit load factor
    if (nReach > 1 && Math.abs(t0.de + dePerG * (nReach - 1)) > deLim) warnings.push(`Elevator travel is insufficient to reach ${nReach < i.n_limit ? 'maximum lift' : 'the limit load factor'} (${nReach.toFixed(2)} g) at this speed.`);
    return {
      kpis: [
        { key: 'trim_alpha_deg', label: 'Trim angle of attack (body)', value: deg(t0.alpha), unit: 'deg' }, { key: 'trim_de_deg', label: 'Trim elevator angle', value: deg(t0.de), unit: 'deg', status: Math.abs(t0.de) < 0.6 * deLim ? 'ok' : 'warn', note: 'Trailing edge down positive' },
        { key: 'trim_CL', label: 'Trim lift coefficient', value: D.CL, unit: '-' }, { key: 'trim_thrust_N', label: 'Thrust required', value: t0.T, unit: 'N' },
        { key: 'trim_throttle', label: 'Throttle setting', value: t0.throttle, unit: '-', status: t0.throttle <= 1 ? 'ok' : t0.throttle <= 1.1 ? 'warn' : 'bad', note: 'Thrust available from a generic σ^0.7 altitude lapse, good to about ±10%' },
        { key: 'de_per_g_deg', label: 'Elevator angle per g', value: deg(dePerG), unit: 'deg/g', status: dePerG < 0 ? 'ok' : 'bad', note: 'Negative (stick back to pull g) for manoeuvre stability' },
        { key: 'manoeuvre_margin', label: 'Manoeuvre margin', value: D.hm - i.cg, unit: 'MAC' },
        { key: 'cg_fwd_limit', label: 'Forward CG limit (trim at CLmax)', value: fwd, unit: 'MAC' }, { key: 'cg_aft_limit', label: 'Aft CG limit (5% static margin)', value: D.hn - 0.05, unit: 'MAC' },
        { key: 'V_stall_trim_ms', label: 'Stall speed at this mass and altitude', value: Vs, unit: 'm/s' },
      ].filter((k) => fin(k.value) || k.key.startsWith('trim_')),
      plots: [
        { type: 'line', title: 'Trim angles', xlabel: 'True airspeed [m/s]', ylabel: 'Angle [deg]', series: [{ name: 'Angle of attack', x: Vs_, y: T.map((t) => deg(t.alpha)) }, { name: 'Elevator', x: Vs_, y: T.map((t) => deg(t.de)) }], annotations: [{ x: i.V, label: 'Flight point' }] },
        { type: 'line', title: 'Thrust required and available', xlabel: 'True airspeed [m/s]', ylabel: 'Thrust [kN]', series: [{ name: 'Required', x: Vs_, y: T.map((t) => t.T / 1e3) }, { name: 'Available', x: Vs_, y: T.map((t) => t.Tmax / 1e3) }] },
        { type: 'line', title: 'Elevator to trim against lift coefficient', xlabel: 'CL [-]', ylabel: 'Elevator angle [deg]', series: [0, -0.1, 0.1].map((dh) => ({ name: `CG ${(100 * (i.cg + dh)).toFixed(0)}% MAC`, x: T.map((t) => t.D.CL), y: Vs_.map((V) => deg(trimAt(withCg(i, i.cg + dh), V).de)) })) },
      ],
      warnings, models: ['Linear lift and pitching-moment trim', 'Parabolic drag polar', 'Steady pull-up manoeuvre equations'], assumptions: ['Level flight, thrust line through the CG', 'Fixed stabiliser: all trim is carried by the elevator', 'Derivatives from the preliminary build-up', 'Thrust available: static thrust (jets) or shaft power × propeller efficiency (propellers) with a generic σ^0.7 altitude lapse and none for electric motors — a typical value, not engine data; it under-predicts flat-rated turboprops'],
    };
  },
  verify() {
    const b = { ...Object.fromEntries(FW.map((f) => [f.key, f.default])), de_max_deg: 25, n_limit: 2.5 }, t = trimAt(b, b.V), D = t.D;
    return [N.check('Trim satisfies CL = W/(qS)', D.CLa * (t.alpha - rad(b.alpha0_deg)) + D.CLde * t.de, D.CL, 1e-10, 'Lift equation'), N.check('Trim satisfies Cm = 0', b.Cm0 + D.Cma * (t.alpha - rad(b.alpha0_deg)) + D.Cmde * t.de + 1, 1, 1e-10, 'Moment equation')];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.trim_throttle > 1) out.push({ severity: o.trim_throttle > 1.1 ? 'critical' : 'warn', title: o.trim_throttle > 1.1 ? 'Not enough thrust to hold this flight condition' : 'Thrust is marginal at this flight condition', detail: `Throttle required is ${(100 * o.trim_throttle).toFixed(0)}% of the estimated thrust available.${o.trim_throttle > 1.1 ? '' : ' That is inside the ±10% accuracy of the generic altitude-lapse model, so the condition is marginal rather than shown to be infeasible.'}`, action: 'Reduce speed, altitude or mass, or replace the lapse estimate with the installed thrust or power from Suite 7.', basis: 'Thrust–drag balance with a generic σ^0.7 thrust/power lapse' });
    if (Math.abs(o.trim_de_deg) > 0.6 * i.de_max_deg) out.push({ severity: 'advise', title: 'High trim elevator', detail: `${o.trim_de_deg.toFixed(1)}° of ${i.de_max_deg}° travel.`, action: 'Set stabiliser incidence or move the CG aft to cut trim drag and keep control authority.', basis: 'Trim equations' });
    if (o.de_per_g_deg >= 0) out.push({ severity: 'critical', title: 'No manoeuvre stability', detail: 'The CG is aft of the manoeuvre point.', action: 'Move the CG forward or add pitch damping (Suite 16).', basis: 'Stick-fixed manoeuvre point' });
    return out;
  },
};

// ---- 3. longitudinal modes --------------------------------------------------------------------
const longAn = {
  id: 'long', title: 'Longitudinal dynamics: short period and phugoid', fidelity: 'reduced-order',
  summary: 'The linearised longitudinal state equations about trim: eigenvalues, short-period and phugoid frequency and damping against handling-qualities limits, root locus with CG and the response to an elevator doublet.',
  equations: ['Longitudinal linearised state equations', 'Characteristic eigenvalue equations', 'Short-period approximation equations', 'Phugoid approximation equations', 'Stability derivative equations'],
  applicable: hasWing,
  inputs: [...FW,
    { key: 'sp_zeta_min', label: 'Short-period damping, lower limit', unit: '-', default: 0.3, min: 0, max: 1, group: 'Handling criteria', help: 'MIL-F-8785C Table IV, Level 1: 0.30 (Category B), 0.35 (Categories A and C); Level 2: 0.20 and 0.25' },
    { key: 'sp_zeta_max', label: 'Short-period damping, upper limit', unit: '-', default: 2, min: 0.5, max: 5, group: 'Handling criteria', help: 'MIL-F-8785C Table IV, Level 1: 2.0 (Category B), 1.3 (Categories A and C)' },
    { key: 'ph_zeta_min', label: 'Phugoid damping, lower limit', unit: '-', default: 0.04, min: -0.5, max: 1, group: 'Handling criteria', help: 'MIL-F-8785C para 3.2.1.2, Level 1: 0.04 (Level 2: 0)' },
    { key: 'cap_min', label: 'Control anticipation parameter, lower limit', unit: '1/(g·s²)', default: 0.085, min: 0, max: 5, group: 'Handling criteria', help: 'MIL-F-8785C para 3.2.2.1.1, Level 1: 0.085 (Category B), 0.28 (Category A)' },
    { key: 'cap_max', label: 'Control anticipation parameter, upper limit', unit: '1/(g·s²)', default: 3.6, min: 0.1, max: 20, group: 'Handling criteria', help: 'MIL-F-8785C Level 1: 3.6, read from the Category A figure and taken to apply to Category B as well' },
    { key: 'de_doublet_deg', label: 'Elevator doublet amplitude', unit: 'deg', default: 2, min: 0.1, max: 15, group: 'Manoeuvre' }],
  defaults: fwDefaults,
  run(i) {
    const s = longSS(i), ev = N.eig(s.A).map(modeOf), osc = ev.filter((m) => m.im > 1e-9).sort((p, q) => q.wn - p.wn), d = s.d, u0 = i.V, warnings = [];
    const sp = osc[0], ph = osc.length > 1 ? osc[1] : null, spA = { wn: Math.sqrt(Math.max(d.Zw * d.Mq - u0 * d.Mw, 0)), z: 0 }; spA.z = -(d.Mq + d.Zw + d.Mwd * u0) / (2 * spA.wn || 1);
    const phA = { wn: (Math.SQRT2 * G0) / u0, z: -d.Xu / (2 * ((Math.SQRT2 * G0) / u0)) }, nza = (u0 * -d.Zw) / G0, cap = sp ? sp.wn ** 2 / nza : NaN;
    if (!sp) warnings.push('No oscillatory short-period mode: the pitch dynamics are overdamped or statically unstable (real roots).');
    if (!ph) warnings.push('The phugoid has degenerated into real roots.');
    if (ev.some((m) => m.re > 1e-9)) warnings.push(`An unstable longitudinal root exists (time to double ${(Math.LN2 / N.amax(ev.map((m) => m.re))).toFixed(1)} s).`);
    const tEnd = ph ? Math.min(3 * (TAU / ph.im), 600) : 60, r = linSim(s.A, s.B, (t) => [doublet(rad(i.de_doublet_deg), 1, 1.5)(t), 0], tEnd, simSteps(ev, tEnd, 3000)), tS = sp ? Math.min(tEnd, 8 * (TAU / sp.im)) : 20;
    const cgs = N.linspace(i.cg - 0.15, i.cg + 0.4, 23), loc = { x: [], y: [] }; for (const h of cgs) for (const l of N.eig(longSS(withCg(i, h)).A)) { loc.x.push(l[0]); loc.y.push(l[1]); }
    const ok = (v, lo, hi) => (fin(v) && v >= lo && v <= hi ? 'ok' : 'warn'), kS = thin(r.t.map((t, k) => k).filter((k) => r.t[k] <= tS), 400);
    return {
      kpis: [
        { key: 'sp_omega_rads', label: 'Short-period natural frequency', value: sp ? sp.wn : NaN, unit: 'rad/s' }, { key: 'sp_zeta', label: 'Short-period damping ratio', value: sp ? sp.zeta : NaN, unit: '-', status: ok(sp?.zeta, i.sp_zeta_min, i.sp_zeta_max), note: `Criterion ${i.sp_zeta_min}–${i.sp_zeta_max}` },
        { key: 'sp_period_s', label: 'Short-period period', value: sp ? TAU / sp.im : NaN, unit: 's' },
        { key: 'ph_period_s', label: 'Phugoid period', value: ph ? TAU / ph.im : NaN, unit: 's' }, { key: 'ph_zeta', label: 'Phugoid damping ratio', value: ph ? ph.zeta : NaN, unit: '-', status: ph && ph.zeta >= i.ph_zeta_min ? 'ok' : 'warn', note: `Criterion ≥ ${i.ph_zeta_min}` },
        { key: 'cap', label: 'Control anticipation parameter', value: cap, unit: '1/(g·s²)', status: ok(cap, i.cap_min, i.cap_max), note: `ωsp²/(n/α); criterion ${i.cap_min}–${i.cap_max}` },
        { key: 'n_alpha', label: 'Load factor per angle of attack n/α', value: nza, unit: 'g/rad' },
        { key: 'sp_omega_approx', label: 'Short-period frequency, two-state approximation', value: spA.wn, unit: 'rad/s' }, { key: 'sp_zeta_approx', label: 'Short-period damping, approximation', value: spA.z, unit: '-' },
        { key: 'ph_period_lanchester_s', label: 'Phugoid period, Lanchester π√2·V/g', value: TAU / phA.wn, unit: 's' }, { key: 'ph_zeta_approx', label: 'Phugoid damping, approximation', value: phA.z, unit: '-', note: '≈ 1/(√2·L/D)' },
      ].filter((k) => fin(k.value) || ['sp_omega_rads', 'sp_zeta', 'ph_period_s', 'ph_zeta'].includes(k.key)),
      plots: [
        { type: 'line', title: 'Longitudinal roots as the CG moves aft', xlabel: 'Real part [1/s]', ylabel: 'Imaginary part [rad/s]', series: [{ name: `CG ${(100 * cgs[0]).toFixed(0)}% → ${(100 * cgs[22]).toFixed(0)}% MAC`, x: loc.x, y: loc.y, style: 'points' }, { name: 'Current CG', x: ev.map((m) => m.re), y: ev.map((m) => m.im), style: 'points' }], annotations: [{ x: 0, label: 'Stability boundary' }] },
        { type: 'line', title: 'Short-term response to an elevator doublet', xlabel: 'Time [s]', ylabel: 'Angle or rate [deg, deg/s]', series: [{ name: 'Angle of attack', x: kS.map((k) => r.t[k]), y: kS.map((k) => deg(r.x[k][1] / u0)) }, { name: 'Pitch rate', x: kS.map((k) => r.t[k]), y: kS.map((k) => deg(r.x[k][2])) }, { name: 'Pitch attitude', x: kS.map((k) => r.t[k]), y: kS.map((k) => deg(r.x[k][3])) }] },
        { type: 'line', title: 'Long-term response: phugoid', xlabel: 'Time [s]', ylabel: 'Speed change [m/s]', series: [{ name: 'Speed perturbation', x: thin(r.t), y: thin(r.x.map((x) => x[0])) }] },
      ],
      tables: [{ title: 'Eigenvalues', columns: ['Real [1/s]', 'Imaginary [rad/s]', 'ωn [rad/s]', 'ζ [-]'], rows: ev.map((m) => [m.re, m.im, m.wn, m.zeta]) }, { title: 'Dimensional derivatives', columns: ['Derivative', 'Value'], rows: Object.entries(d).map(([k, v]) => [k, v]) }],
      outputs: { A_long: s.A, B_long: s.B, u0_ms: u0 },
      warnings, models: ['Linear state-space model, states u, w, q, θ in stability axes', 'Eigenvalue analysis (QR)', 'Short-period and Lanchester phugoid approximations for comparison'],
      assumptions: ['Small perturbations about steady level flight', 'Constant-thrust jets, constant-power propellers for the speed derivative', 'Handling criteria are inputs; the defaults are the MIL-F-8785C Level 1 limits for Category B (cruise) flight phases — a military flying-qualities specification used here as a general yardstick, so select the class and flight-phase category that apply'],
    };
  },
  calibration: { params: [{ key: 'k_Cma', min: 0.3, max: 2 }, { key: 'k_Cmq', min: 0.3, max: 2.5 }], sweep: 'V', target: 'sp_omega_rads', note: 'Flight-test short-period frequency and damping from elevator doublets at several speeds' },
  verify() {
    const b = { ...Object.fromEntries(longAn.inputs.map((f) => [f.key, f.default])), V: 120, alt_m: 0, compress: false, CD0: 1e-7, e: 1e7, ptype: 'turbojet', cg: -0.6, k_Cmq: 0.1 }, o = N.kv(longAn.run(b)), s = longSS(b), l = N.eig(s.A), tr = N.sum(l.map((x) => x[0]));
    return [N.check('Drag-free phugoid period equals Lanchester π√2·V/g', o.ph_period_s, (Math.PI * Math.SQRT2 * 120) / G0, 0.01, 'Lanchester (1908)'),
      N.check('Sum of eigenvalues equals the trace of A', tr, s.A[0][0] + s.A[1][1] + s.A[2][2], 1e-8, 'Linear algebra identity'),
      N.check('Short-period two-state approximation', o.sp_omega_approx, o.sp_omega_rads, 0.02, 'Etkin & Reid, Dynamics of Flight')];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!fin(o.sp_zeta) || o.sp_zeta < i.sp_zeta_min) out.push({ severity: fin(o.sp_zeta) && o.sp_zeta >= 0.2 ? 'advise' : 'warn', title: 'Short-period damping below the criterion', detail: fin(o.sp_zeta) ? `ζ = ${o.sp_zeta.toFixed(2)} against ${i.sp_zeta_min}.` : 'No oscillatory short period was found.', action: 'Increase tail volume or add a pitch damper (Suite 16, stability augmentation).', basis: 'Short-period damping criterion (input; MIL-F-8785C Table IV: Level 1 default, Level 2 minimum 0.20 in Category B)' });
    if (fin(o.ph_zeta) && o.ph_zeta < i.ph_zeta_min) out.push({ severity: 'advise', title: 'Lightly damped phugoid', detail: `ζ = ${o.ph_zeta.toFixed(3)}; period ${o.ph_period_s.toFixed(0)} s.`, action: 'Normal for clean, efficient aircraft; an altitude- or speed-hold autopilot removes the workload.', basis: 'Phugoid damping criterion (input)' });
    if (fin(o.cap) && (o.cap < i.cap_min || o.cap > i.cap_max)) out.push({ severity: 'advise', title: 'Control anticipation parameter outside the band', detail: `CAP = ${o.cap.toFixed(2)}.`, action: 'Adjust static margin or pitch inertia, or shape the response with the flight control law.', basis: 'CAP criterion (input)' });
    return out;
  },
};

// ---- 4. lateral-directional modes -------------------------------------------------------------
const latAn = {
  id: 'lat', title: 'Lateral-directional dynamics: Dutch roll, roll and spiral', fidelity: 'reduced-order',
  summary: 'The lateral-directional state equations: Dutch-roll frequency and damping, roll time constant and spiral stability against handling-qualities limits, with responses to aileron and rudder inputs.',
  equations: ['Lateral-directional linearised state equations', 'Dutch-roll equations', 'Roll subsidence equations', 'Spiral stability equations', 'Characteristic eigenvalue equations'],
  applicable: hasWing,
  inputs: [...FW,
    { key: 'dr_zeta_min', label: 'Dutch-roll damping, lower limit', unit: '-', default: 0.08, min: 0, max: 1, group: 'Handling criteria', help: 'MIL-F-8785C Table VI, Level 1: 0.08 (Category B), 0.19 (Category A); Level 2: 0.02' },
    { key: 'dr_zw_min', label: 'Dutch-roll ζ·ωn, lower limit', unit: 'rad/s', default: 0.15, min: 0, max: 2, group: 'Handling criteria', help: 'MIL-F-8785C Table VI, Level 1: 0.15 (Category B), 0.35 (Category A); Level 2: 0.05' },
    { key: 'dr_w_min', label: 'Dutch-roll frequency, lower limit', unit: 'rad/s', default: 0.4, min: 0, max: 5, group: 'Handling criteria', help: 'MIL-F-8785C Table VI: 0.4 (1.0 for Class I and IV aircraft in Category A)' },
    { key: 'roll_tau_max', label: 'Roll time constant, upper limit', unit: 's', default: 1.4, min: 0.1, max: 10, group: 'Handling criteria', help: 'MIL-F-8785C Table VII, Level 1: 1.4 s (1.0 s for Class I and IV aircraft in Category A and Class I, II-C and IV in Category C)' },
    { key: 'spiral_T2_min', label: 'Spiral time to double, lower limit', unit: 's', default: 20, min: 1, max: 200, group: 'Handling criteria', help: 'MIL-F-8785C Table VIII, Level 1: 20 s (Category B), 12 s (Categories A and C)' },
    { key: 'da_deg', label: 'Aileron step', unit: 'deg', default: 5, min: 0.1, max: 30, group: 'Manoeuvre' }, { key: 'dr_deg', label: 'Rudder doublet', unit: 'deg', default: 5, min: 0.1, max: 30, group: 'Manoeuvre' }],
  defaults: fwDefaults,
  run(i) {
    const s = latSS(i), ev = N.eig(s.A).map(modeOf), d = s.d, u0 = i.V, dr = ev.filter((m) => m.im > 1e-9).sort((p, q) => q.wn - p.wn)[0], real = ev.filter((m) => Math.abs(m.im) <= 1e-9).sort((p, q) => p.re - q.re), warnings = [];
    const roll = real[0], spiral = real[real.length - 1] !== roll ? real[real.length - 1] : null, tauR = roll ? -1 / roll.re : NaN, T2 = spiral ? Math.LN2 / spiral.re : NaN;
    const drA = Math.sqrt(Math.max((d.Yb * d.Nr - d.Nb * d.Yr + u0 * d.Nb) / u0, 0)), spA = (d.Lb * d.Nr - d.Lr * d.Nb) / (d.Lb || 1e-12);
    if (!dr) warnings.push('No oscillatory Dutch-roll mode was found (real roots): check directional stability Cnβ.');
    if (dr && dr.zeta < 0) warnings.push(`The bare-airframe Dutch roll is ${Math.LN2 / dr.re < 120 ? `unstable (time to double ${(Math.LN2 / dr.re).toFixed(0)} s)` : 'practically undamped'}: it needs a yaw damper (Suite 16).`);
    if (real.length < 2) warnings.push('Roll and spiral have coupled into a lateral phugoid (roll–spiral oscillation).');
    const ra = linSim(s.A, s.B, (t) => [t > 0.5 ? rad(i.da_deg) : 0, 0], 12, simSteps(ev, 12, 1200)), tR = dr ? Math.min(60, 6 * (TAU / dr.im)) : 20, rr = linSim(s.A, s.B, (t) => [0, doublet(rad(i.dr_deg), 0.5, 1)(t)], tR, simSteps(ev, tR, 1500));
    const pss = (-d.Lda / (d.Lp || -1e-12)) * rad(i.da_deg), phiBeta = dr ? N.C.abs(modeShape(s.A, dr)[2]) : NaN;
    const st = (c) => (c ? 'ok' : 'warn');
    return {
      kpis: [
        { key: 'dr_omega_rads', label: 'Dutch-roll natural frequency', value: dr ? dr.wn : NaN, unit: 'rad/s', status: st(dr && dr.wn >= i.dr_w_min), note: `Criterion ≥ ${i.dr_w_min}` },
        { key: 'dr_zeta', label: 'Dutch-roll damping ratio', value: dr ? dr.zeta : NaN, unit: '-', status: st(dr && dr.zeta >= i.dr_zeta_min && dr.zeta * dr.wn >= i.dr_zw_min), note: `Criteria ζ ≥ ${i.dr_zeta_min}, ζωn ≥ ${i.dr_zw_min}` },
        { key: 'dr_period_s', label: 'Dutch-roll period', value: dr ? TAU / dr.im : NaN, unit: 's' },
        { key: 'roll_tau_s', label: 'Roll-mode time constant', value: tauR, unit: 's', status: st(tauR > 0 && tauR <= i.roll_tau_max), note: `Criterion ≤ ${i.roll_tau_max} s` },
        { key: 'spiral_T2_s', label: 'Spiral time to double amplitude', value: T2, unit: 's', status: st(!(T2 > 0) || T2 >= i.spiral_T2_min), note: `Negative = stable (time to half). Criterion: not faster than ${i.spiral_T2_min} s to double` },
        { key: 'phi_beta_ratio', label: 'Dutch-roll |φ/β| ratio', value: phiBeta, unit: '-' },
        { key: 'roll_rate_ss_dps', label: 'Steady roll rate for the aileron step', value: deg(pss), unit: 'deg/s' }, { key: 'pb_2V', label: 'Helix angle pb/2V', value: (pss * i.b) / (2 * u0), unit: 'rad', note: 'About 0.07 is a common roll-performance target' },
        { key: 'dr_omega_approx', label: 'Dutch-roll frequency, approximation', value: drA, unit: 'rad/s' }, { key: 'spiral_root_approx', label: 'Spiral root, approximation', value: spA, unit: '1/s' },
      ].filter((k) => fin(k.value) || ['dr_omega_rads', 'dr_zeta', 'roll_tau_s', 'spiral_T2_s'].includes(k.key)),
      plots: [
        { type: 'line', title: 'Lateral-directional roots', xlabel: 'Real part [1/s]', ylabel: 'Imaginary part [rad/s]', series: [{ name: 'Eigenvalues', x: ev.map((m) => m.re), y: ev.map((m) => m.im), style: 'points' }], annotations: [{ x: 0, label: 'Stability boundary' }] },
        { type: 'line', title: 'Response to an aileron step', xlabel: 'Time [s]', ylabel: 'Angle or rate [deg, deg/s]', series: [{ name: 'Roll rate', x: thin(ra.t), y: thin(ra.x.map((x) => deg(x[1]))) }, { name: 'Bank angle', x: thin(ra.t), y: thin(ra.x.map((x) => N.clamp(deg(x[3]), -720, 720))) }, { name: 'Sideslip', x: thin(ra.t), y: thin(ra.x.map((x) => deg(x[0]))) }] },
        { type: 'line', title: 'Response to a rudder doublet', xlabel: 'Time [s]', ylabel: 'Angle or rate [deg, deg/s]', series: [{ name: 'Sideslip', x: thin(rr.t), y: thin(rr.x.map((x) => deg(x[0]))) }, { name: 'Yaw rate', x: thin(rr.t), y: thin(rr.x.map((x) => deg(x[2]))) }, { name: 'Roll rate', x: thin(rr.t), y: thin(rr.x.map((x) => deg(x[1]))) }] },
      ],
      tables: [{ title: 'Eigenvalues', columns: ['Real [1/s]', 'Imaginary [rad/s]', 'ωn [rad/s]', 'ζ [-]'], rows: ev.map((m) => [m.re, m.im, m.wn, m.zeta]) }, { title: 'Dimensional derivatives', columns: ['Derivative', 'Value'], rows: Object.entries(d).map(([k, v]) => [k, v]) }],
      outputs: { A_lat: s.A, B_lat: s.B },
      warnings, models: ['Linear state-space model, states β, p, r, φ', 'Eigenvalue analysis (QR)', 'Dutch-roll and spiral approximations for comparison'],
      assumptions: ['Small perturbations about wings-level flight; product of inertia included', 'Derivatives from the preliminary build-up: cross derivatives carry the largest uncertainty', 'Handling criteria are inputs; the defaults are the MIL-F-8785C Level 1 limits for Category B (cruise) flight phases (Tables VI–VIII), used here as a general yardstick', 'A bare-airframe result: a yaw damper or other augmentation is not included'],
    };
  },
  calibration: { params: [{ key: 'k_Cnb', min: 0.3, max: 2.5 }, { key: 'k_Clp', min: 0.3, max: 2 }], sweep: 'V', target: 'dr_omega_rads', note: 'Flight-identified Dutch-roll frequency and roll time constant from rudder doublets and aileron steps' },
  verify() {
    const b = Object.fromEntries(latAn.inputs.map((f) => [f.key, f.default])), s = latSS(b), l = N.eig(s.A), o = N.kv(latAn.run(b)), C = N.C;
    // Dutch-roll eigenvector: the sideslip row, which is not used to build the mode shape, must be satisfied too
    const dr = l.map(modeOf).filter((m) => m.im > 1e-9).sort((p, q) => q.wn - p.wn)[0], v = modeShape(s.A, dr), row0 = v.reduce((a, x, k) => C.add(a, C.scale(x, s.A[0][k + 1])), [s.A[0][0] - dr.re, -dr.im]);
    // pure rolling motion: with all coupling removed the roll root is Lp exactly
    const A1 = s.A.map((r) => r.slice()); A1[1][0] = A1[1][2] = 0; const lr = N.eig(A1).map((x) => x[0]);
    return [N.check('Sum of eigenvalues equals the trace of A', N.sum(l.map((x) => x[0])), s.A[0][0] + s.A[1][1] + s.A[2][2], 1e-8, 'Linear algebra identity'),
      N.check('Uncoupled roll root equals Lp', lr.reduce((p, q) => (Math.abs(q - A1[1][1]) < Math.abs(p - A1[1][1]) ? q : p)), A1[1][1], 1e-8, 'Single-degree-of-freedom roll subsidence τ = −1/Lp'),
      N.check('Dutch-roll mode shape satisfies the sideslip equation', C.abs(row0) / dr.wn, 0, 1e-9, '(A − λI)·v = 0'),
      N.check('Dutch-roll approximation within its expected accuracy', o.dr_omega_approx, o.dr_omega_rads, 0.25, 'Nelson, Flight Stability and Automatic Control (approximation of modest accuracy)')];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!fin(o.dr_zeta)) out.push({ severity: 'warn', title: 'No Dutch-roll oscillation found', detail: 'The lateral-directional roots are all real.', action: 'Check the directional stability Cnβ and the inertia data.', basis: 'Eigenvalues of the lateral-directional model' });
    else if (o.dr_zeta < i.dr_zeta_min || o.dr_zeta * o.dr_omega_rads < i.dr_zw_min || o.dr_omega_rads < i.dr_w_min) {
      // bare-airframe result: Level 2 minima and the divergence rate decide how much the aircraft depends on its yaw damper
      const w = o.dr_omega_rads, zw = o.dr_zeta * w, lvl2 = o.dr_zeta >= 0.02 && zw >= 0.05 && w >= 0.4, t2 = zw < 0 ? Math.LN2 / -zw : Infinity;
      out.push({ severity: t2 < 4 ? 'critical' : lvl2 ? 'advise' : 'warn', title: 'Bare-airframe Dutch roll needs a yaw damper', detail: `ζ = ${o.dr_zeta.toFixed(3)}, ωn = ${w.toFixed(2)} rad/s against ζ ≥ ${i.dr_zeta_min}, ζωn ≥ ${i.dr_zw_min} rad/s and ωn ≥ ${i.dr_w_min} rad/s${!fin(t2) ? '' : t2 < 120 ? `; the oscillation doubles in ${t2.toFixed(0)} s` : '; the oscillation is practically undamped'}. Low Dutch-roll damping is usual for swept-wing aircraft at altitude and is corrected by yaw-rate feedback rather than by redesign.`, action: 'Design the yaw damper in Suite 16 (stability augmentation) and treat it as required equipment in Suite 22; if it cannot be relied on, a larger fin or less dihedral effect raises the bare-airframe damping. The cross derivatives behind this estimate carry about ±50% uncertainty: calibrate them before drawing conclusions.', basis: 'Dutch-roll criteria (inputs; the defaults are the MIL-F-8785C Table VI Level 1 Category B minima, Level 2 being ζ ≥ 0.02 and ζωn ≥ 0.05 rad/s)' });
    }
    if (o.roll_tau_s > i.roll_tau_max) out.push({ severity: o.roll_tau_s > 3 ? 'warn' : 'advise', title: 'Sluggish roll response', detail: `Roll time constant ${o.roll_tau_s.toFixed(2)} s against ${i.roll_tau_max} s.`, action: 'Reduce roll inertia (fuel and stores inboard) or add roll-rate feedback.', basis: 'Roll-mode criterion (input; MIL-F-8785C Table VII: Level 1 default, Level 2 limit 3.0 s)' });
    if (o.spiral_T2_s > 0 && o.spiral_T2_s < i.spiral_T2_min) out.push({ severity: o.spiral_T2_s < 8 ? 'warn' : 'advise', title: o.spiral_T2_s < 8 ? 'Rapid spiral divergence' : 'Spiral divergence faster than the Level 1 criterion', detail: `Time to double ${o.spiral_T2_s.toFixed(1)} s against ${i.spiral_T2_min} s.`, action: 'Increase dihedral effect or reduce fin area; a wing-leveller autopilot mode suppresses it, and an aircraft that is always flown through an autopilot is not affected.', basis: 'Spiral criterion (input; MIL-F-8785C Table VIII: Level 1 default, Level 2 limit 8 s)' });
    return out;
  },
};

// ---- rigid-body core and hover model ----------------------------------------------------------
/** Body → NED direction cosine matrix from quaternion [q0, q1, q2, q3]. */
const dcm = (q) => { const [a, b, c, d] = q; return [[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)], [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)], [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d]]; };
const euler = (q) => { const R = dcm(q); return [Math.atan2(R[2][1], R[2][2]), -Math.asin(N.clamp(R[2][0], -1, 1)), Math.atan2(R[1][0], R[0][0])]; };
/**
 * Newton–Euler rigid-body derivative. x = [pN, pE, pD, u, v, w, p, q, r, q0..q3, W]; fm(x, t) returns body force and moment
 * excluding gravity. W integrates the power of those forces so that energy can be audited.
 */
function rigidBody(m, I, Ii, fm, kq = 1) {
  return (t, x) => {
    const v = [x[3], x[4], x[5]], w = [x[6], x[7], x[8]], q = [x[9], x[10], x[11], x[12]], R = dcm(q), { F, M } = fm(x, t, R), g = [R[2][0] * G0, R[2][1] * G0, R[2][2] * G0];
    const wxv = N.cross(w, v), Iw = N.matvec(I, w), wd = N.matvec(Ii, N.vadd(M, N.cross(w, Iw), -1)), e = kq * (1 - N.dot(q, q));
    return [...N.matvec(R, v), F[0] / m + g[0] - wxv[0], F[1] / m + g[1] - wxv[1], F[2] / m + g[2] - wxv[2], ...wd,
      0.5 * (-q[1] * w[0] - q[2] * w[1] - q[3] * w[2]) + e * q[0], 0.5 * (q[0] * w[0] + q[2] * w[2] - q[3] * w[1]) + e * q[1], 0.5 * (q[0] * w[1] + q[3] * w[0] - q[1] * w[2]) + e * q[2], 0.5 * (q[0] * w[2] + q[1] * w[1] - q[2] * w[0]) + e * q[3], N.dot(F, v) + N.dot(M, w)];
  };
}
const energy = (m, I, x) => 0.5 * m * (x[3] ** 2 + x[4] ** 2 + x[5] ** 2) + 0.5 * N.dot([x[6], x[7], x[8]], N.matvec(I, [x[6], x[7], x[8]])) + m * G0 * -x[2];
const HOV = [
  { key: 'kind', label: 'Vehicle kind', type: 'select', options: ['helicopter', 'multirotor'], default: 'helicopter', group: 'Vehicle' },
  { key: 'mass_kg', label: 'Mass', unit: 'kg', default: 9980, min: 0.01, group: 'Mass' }, { key: 'Ixx', label: 'Roll inertia', unit: 'kg·m²', default: 6300, min: 1e-6, group: 'Mass' },
  { key: 'Iyy', label: 'Pitch inertia', unit: 'kg·m²', default: 52000, min: 1e-6, group: 'Mass' }, { key: 'Izz', label: 'Yaw inertia', unit: 'kg·m²', default: 50000, min: 1e-6, group: 'Mass' },
  { key: 'R', label: 'Rotor radius', unit: 'm', default: 8.18, min: 0.02, group: 'Rotor' }, { key: 'n_rotors', label: 'Lifting rotors', unit: '', default: 1, min: 1, max: 16, step: 1, discrete: true, group: 'Rotor' },
  { key: 'n_blades', label: 'Blades per rotor', unit: '', default: 4, min: 2, max: 9, step: 1, discrete: true, group: 'Rotor' }, { key: 'chord', label: 'Blade chord', unit: 'm', default: 0.53, min: 0.002, group: 'Rotor' },
  { key: 'rpm', label: 'Rotor speed', unit: 'rpm', default: 258, min: 10, group: 'Rotor' }, { key: 'cla', label: 'Blade lift-curve slope', unit: '1/rad', default: 5.73, min: 2, max: 7, group: 'Rotor' },
  { key: 'lock', label: 'Lock number', unit: '-', default: 8, min: 0.5, max: 20, group: 'Rotor' }, { key: 'e_hinge', label: 'Flap hinge offset / R', unit: '-', default: 0.047, min: 0, max: 0.3, group: 'Rotor' },
  { key: 'm_blade', label: 'Blade mass', unit: 'kg', default: 116, min: 1e-4, group: 'Rotor' },
  { key: 'h_rotor', label: 'Rotor plane above the CG', unit: 'm', default: 1.8, min: -2, max: 6, group: 'Geometry' }, { key: 'arm', label: 'Rotor arm (multirotor) / tail-rotor arm', unit: 'm', default: 9.9, min: 0.02, group: 'Geometry' },
  { key: 'tr_R', label: 'Tail rotor radius (helicopter)', unit: 'm', default: 1.68, min: 0.01, group: 'Geometry' },
  { key: 'k_drag', label: 'Rotor in-plane drag damping (multirotor)', unit: '1/s', default: 0.25, min: 0, max: 3, group: 'Aerodynamics', help: 'Translational damping −Xu from rotor H-force and blade flapping; 0.1–0.5 1/s is typical of small multirotors' },
  { key: 'alt_m', label: 'Altitude', unit: 'm', default: 0, min: -500, max: 8000, group: 'Flight condition' },
];
const hovDefaults = (c, up, d) => {
  const heli = c.meta.type === 'helicopter', m = c.mass.mtow_kg, n = heli ? 1 : Math.max(1, c.prop.n_eng), L = Math.max(c.fuselage.len_m, 4 * c.rotor.R_m, 0.3), k = (v, e) => (v > 0 ? v : e);
  return { kind: heli ? 'helicopter' : 'multirotor', mass_kg: m, Ixx: k(c.mass.Ixx, 0.02 * m * L * L), Iyy: k(c.mass.Iyy, 0.02 * m * L * L), Izz: k(c.mass.Izz, 0.035 * m * L * L), R: c.rotor.R_m || undefined, n_rotors: n, n_blades: c.rotor.n_blades || 2, chord: c.rotor.chord_m || undefined, rpm: c.rotor.rpm || undefined,
    cla: c.rotor.cla, lock: c.rotor.lock, e_hinge: heli ? c.rotor.hinge_offset : 0, m_blade: c.rotor.blade_mass_kg || undefined, h_rotor: heli ? 0.75 * c.fuselage.dia_m : 0.15 * c.rotor.R_m, arm: (heli ? c.rotor.tr_arm_m || 1.2 * c.rotor.R_m : c.wing.b_m > 0 ? 0.35 * c.wing.b_m : 1.6 * c.rotor.R_m) || undefined, tr_R: c.rotor.tr_R_m || 0.2 * c.rotor.R_m || undefined, alt_m: c.atm.alt_m };
};
const hasRotor = (c) => (c.rotor.R_m > 0 ? true : 'This analysis needs lifting rotors; the case has none.');
/** Hover stability and control derivatives (body axes, per unit mass or inertia). Controls: [longitudinal, lateral, collective, yaw]. */
export function hoverDerivs(i) {
  const rho = isa(i.alt_m).rho, W = i.mass_kg * G0, T = W / i.n_rotors, A = Math.PI * i.R ** 2, Om = (i.rpm * TAU) / 60, vt = Om * i.R, CT = T / (rho * A * vt * vt), lam = Math.sqrt(CT / 2), sig = (i.n_blades * i.chord) / (Math.PI * i.R), a = i.cla;
  const dTdw = (2 * a * sig * rho * A * vt * lam) / (16 * lam + a * sig), Zw = (-i.n_rotors * dTdw) / i.mass_kg, P = i.n_rotors * (1.15 * T * lam * vt + (rho * A * vt ** 3 * sig * 0.01) / 8), Q = P / Om / i.n_rotors;
  const dTdth = (rho * A * vt * vt * sig * a) / 6 / (1 + (a * sig) / (16 * lam)); // thrust per rad of collective with inflow lag-free correction
  let d;
  if (i.kind === 'helicopter') {
    const da1du = ((16 * CT) / (sig * a) + 2 * lam) / vt, da1dq = 16 / (i.lock * Om), Kh = T * i.h_rotor + (i.n_blades / 2) * i.e_hinge * i.R * Om * Om * ((i.m_blade * i.R) / 2);
    const At = Math.PI * i.tr_R ** 2, Ttr = Q / i.arm, vtt = 210, lt = Math.sqrt(Ttr / (2 * rho * At)) / vtt, st = 0.15, dTt = (2 * a * st * rho * At * vtt * lt) / (16 * lt + a * st);
    d = { Xu: -G0 * da1du, Mu: (Kh * da1du) / i.Iyy, Mq: (-Kh * da1dq) / i.Iyy, Yv: -G0 * da1du - dTt / i.mass_kg, Lv: (-Kh * da1du) / i.Ixx, Lp: (-Kh * da1dq) / i.Ixx, Nr: (-dTt * i.arm * i.arm) / i.Izz, Nv: (dTt * i.arm) / i.Izz,
      Mc: Kh / i.Iyy, Lc: Kh / i.Ixx, Xc: G0, Zc: -dTdth / i.mass_kg, Nc: (((rho * At * vtt * vtt * st * a) / 6) * i.arm) / i.Izz, Kh, da1du };
  } else {
    const sx = (i.n_rotors * i.arm * i.arm) / 2, kd = i.k_drag; // Σ of squared lever arms about each horizontal axis for a symmetric layout
    d = { Xu: -kd, Mu: (i.mass_kg * kd * i.h_rotor) / i.Iyy, Mq: (-dTdw * sx) / i.Iyy, Yv: -kd, Lv: (-i.mass_kg * kd * i.h_rotor) / i.Ixx, Lp: (-dTdw * sx) / i.Ixx, Nr: (-2 * i.n_rotors * Q) / (Om * i.Izz), Nv: 0,
      Mc: (W * i.arm) / (Math.SQRT2 * i.Iyy), Lc: (W * i.arm) / (Math.SQRT2 * i.Ixx), Xc: 0, Zc: -G0, Nc: (1.5 * i.n_rotors * Q) / i.Izz, Kh: 0, da1du: 0 };
  }
  const Alon = [[d.Xu, 0, 0, -G0], [0, Zw, 0, 0], [d.Mu, 0, d.Mq, 0], [0, 0, 1, 0]], Blon = [[-d.Xc, 0], [0, d.Zc], [d.Mc, 0], [0, 0]]; // x = [u, w, q, θ], inputs [longitudinal, collective]
  const Alat = [[d.Yv, 0, G0, 0], [d.Lv, d.Lp, 0, 0], [0, 1, 0, 0], [d.Nv, 0, 0, d.Nr]], Blat = [[d.Xc, 0], [d.Lc, 0], [0, 0], [0, d.Nc]]; // x = [v, p, φ, r], inputs [lateral, yaw]
  return { ...d, Zw, CT, lam, sig, P, Q, Om, T, Alon, Blon, Alat, Blat };
}

// ---- 5. six-degree-of-freedom simulation ------------------------------------------------------
const MAN = ['pitch doublet', 'roll step', 'yaw doublet', 'power cut', 'none'];
const sixdof = {
  id: 'sixdof', title: 'Six-degree-of-freedom flight simulation', fidelity: 'numerical',
  summary: 'Non-linear rigid-body simulation with quaternion attitude from a trimmed start: response to control doublets, steps or a power cut, with trajectory, attitude and energy-balance checks.',
  equations: ['Newton–Euler rigid-body equations', 'Six-degree-of-freedom aircraft equations of motion', 'Quaternion attitude propagation equations', 'Direction cosine matrix relations', 'Flight-path kinematic equations', 'Aerodynamic force and moment equations'],
  inputs: [
    { key: 'manoeuvre', label: 'Manoeuvre', type: 'select', options: MAN, default: 'pitch doublet', group: 'Manoeuvre' },
    { key: 'amp', label: 'Control amplitude', unit: 'deg (or % for rotor controls)', default: 2, min: 0, max: 30, group: 'Manoeuvre', help: 'Surface deflection for aeroplanes; percent of control moment authority for rotorcraft' },
    { key: 't_end', label: 'Simulated time', unit: 's', default: 40, min: 1, max: 600, group: 'Manoeuvre' },
    { key: 'gamma0_deg', label: 'Initial flight-path angle', unit: 'deg', default: 0, min: -20, max: 20, group: 'Initial conditions' }, { key: 'psi0_deg', label: 'Initial heading', unit: 'deg', default: 0, min: -180, max: 360, group: 'Initial conditions' },
    { key: 'rtol', label: 'Integrator relative tolerance', unit: '-', default: 1e-6, min: 1e-10, max: 1e-3, group: 'Numerics' },
    ...FW.map((f) => ({ ...f, group: 'Aeroplane · ' + f.group })), ...HOV.filter((f) => !FW.some((g) => g.key === f.key)).map((f) => ({ ...f, group: 'Rotorcraft · ' + f.group })),
    { key: 'rotary', label: 'Use the rotorcraft hover model', type: 'bool', default: false, group: 'Manoeuvre' },
  ],
  defaults: (c, up, d) => (c.wing.S_m2 > 0 ? { ...hovDefaults(c, up, d), ...fwDefaults(c, up, d), rotary: false } : { ...hovDefaults(c, up, d), rotary: true, t_end: 6, amp: c.meta.type === 'helicopter' ? 2 : 0.3 }),
  run(i) {
    const m = i.mass_kg, I = [[i.Ixx, 0, i.rotary ? 0 : -i.Ixz], [0, i.Iyy, 0], [i.rotary ? 0 : -i.Ixz, 0, i.Izz]], Ii = N.inv(I), a = rad(i.amp), warnings = [];
    const pulse = (name) => (i.manoeuvre === name ? (name === 'roll step' ? (t) => (t > 1 ? a : 0) : doublet(a, 1, name === 'yaw doublet' ? 1.5 : 1)) : () => 0), uE = pulse('pitch doublet'), uA = pulse('roll step'), uR = pulse('yaw doublet'), cut = (t) => (i.manoeuvre === 'power cut' && t > 1 ? 0 : 1);
    let fm, x0, trimRes = 0, V0;
    if (i.rotary) {
      const h = hoverDerivs(i), W = m * G0, s = 0.01 * i.amp / (a || 1); // rotor controls in fractions of authority
      fm = (x, t) => ({ F: [m * (h.Xu * x[3] - h.Xc * s * uE(t)), m * (h.Yv * x[4] + h.Xc * s * uA(t)), -W * cut(t) + m * h.Zw * x[5]], M: [i.Ixx * (h.Lv * x[4] + h.Lp * x[6] + h.Lc * s * uA(t)), i.Iyy * (h.Mu * x[3] + h.Mq * x[7] + h.Mc * s * uE(t)), i.Izz * (h.Nv * x[4] + h.Nr * x[8] + h.Nc * s * uR(t))] });
      x0 = [0, 0, -Math.max(i.alt_m, 50), 0, 0, 0, 0, 0, 0, Math.cos(rad(i.psi0_deg) / 2), 0, 0, Math.sin(rad(i.psi0_deg) / 2), 0]; V0 = 0;
    } else {
      const D = derivs(i), a0 = rad(i.alpha0_deg), c = i.mac, b = i.b, gam = rad(i.gamma0_deg);
      const aeroFM = (x, de, da, dr, thr) => {
        const V = Math.max(Math.hypot(x[3], x[4], x[5]), 0.1), al = Math.atan2(x[5], x[3]), be = Math.asin(N.clamp(x[4] / V, -1, 1)), at = isa(-x[2]), qd = 0.5 * at.rho * V * V, ph = (x[6] * b) / (2 * V), qh = (x[7] * c) / (2 * V), rh = (x[8] * b) / (2 * V);
        const CLl = D.CLa * (al - a0) + D.CLq * qh + D.CLde * de, CL = i.CLmax * Math.tanh(CLl / i.CLmax), CD = i.CD0 + D.kI * CL * CL + (Math.abs(CLl) > i.CLmax ? 0.5 * Math.sin(al) ** 2 : 0), CY = D.CYb * be + D.CYr * rh + D.CYdr * dr, T = thr * thrustMax(i, V, at.sigma), ca = Math.cos(al), sa = Math.sin(al);
        return { F: [qd * i.S * (-CD * ca + CL * sa) + T, qd * i.S * CY, qd * i.S * (-CD * sa - CL * ca)], M: [qd * i.S * b * (D.Clb * be + D.Clp * ph + D.Clr * rh + D.Clda * da + D.Cldr * dr), qd * i.S * c * (i.Cm0 + D.Cma * (al - a0) + D.Cmq * qh + D.Cmde * de), qd * i.S * b * (D.Cnb * be + D.Cnp * ph + D.Cnr * rh + D.Cnda * da + D.Cndr * dr)], al, be };
      };
      // trim: body-axis force and pitching-moment balance for steady flight along γ
      const res = (p) => { const th = p[0] + gam, f = aeroFM([0, 0, -i.alt_m, i.V * Math.cos(p[0]), 0, i.V * Math.sin(p[0]), 0, 0, 0], p[1], 0, 0, p[2]); return [f.F[0] / (m * G0) - Math.sin(th), f.F[2] / (m * G0) + Math.cos(th), f.M[1] / (D.q * i.S * c)]; };
      const t0 = trimAt(i, i.V), sol = N.fsolve(res, [t0.alpha, t0.de, N.clamp(t0.throttle || 0.5, 0, 1.5)]), [al, de0, thr0] = sol.x, th0 = al + gam, ps = rad(i.psi0_deg); trimRes = sol.residual;
      if (!sol.converged) warnings.push('The non-linear trim did not converge; the simulation starts out of equilibrium.');
      if (thr0 > 1) warnings.push('Trim needs more than full throttle: the aircraft will decelerate or descend.');
      const cy = Math.cos(ps / 2), sy = Math.sin(ps / 2), cp = Math.cos(th0 / 2), sp = Math.sin(th0 / 2);
      x0 = [0, 0, -i.alt_m, i.V * Math.cos(al), 0, i.V * Math.sin(al), 0, 0, 0, cy * cp, -sy * sp, cy * sp, sy * cp, 0]; V0 = i.V;
      fm = (x, t) => aeroFM(x, de0 + uE(t), uA(t), uR(t), thr0 * cut(t));
    }
    const f = rigidBody(m, I, Ii, fm), sim = N.rk45(f, 0, x0, i.t_end, { rtol: i.rtol, atol: 1e-9, hmax: 0.05, maxSteps: 60000 }), Y = sim.y, T = sim.t, n = T.length;
    const E0 = energy(m, I, x0), eErr = N.amax(Y.map((x) => Math.abs(energy(m, I, x) - E0 - x[13]))) / Math.max(E0, m * G0 * 10, 1), qErr = N.amax(Y.map((x) => Math.abs(Math.hypot(x[9], x[10], x[11], x[12]) - 1)));
    const eu = Y.map((x) => euler([x[9], x[10], x[11], x[12]])), Vt = Y.map((x) => Math.hypot(x[3], x[4], x[5])), al = Y.map((x) => (Vt[0] > 1 ? deg(Math.atan2(x[5], x[3])) : 0));
    const nz = Y.map((x, k) => { const d = fm(x, T[k], null); return -d.F[2] / (m * G0); }), xe = Y[n - 1], ix = thin(N.range(n)), g = (arr) => ix.map((k) => arr[k]);
    if (T[n - 1] < i.t_end - 1e-6) warnings.push('The integration stopped early (step limit): the motion is diverging rapidly.');
    if (!i.rotary && N.amax(al.map(Math.abs)) > 15) warnings.push('Angle of attack exceeded 15°: the simple stall model (lift saturation with a drag rise) is only indicative there.');
    if (i.rotary) warnings.push('The rotorcraft model uses linear hover derivatives with non-linear kinematics: it is valid near hover and for small attitude changes only.');
    if (-xe[2] < 0) warnings.push('The vehicle descended below the reference datum (zero altitude) during the simulation.');
    return {
      kpis: [
        { key: 'trim_residual', label: 'Trim residual', value: trimRes, unit: '-', status: trimRes < 1e-6 ? 'ok' : 'warn' },
        { key: 'energy_err', label: 'Energy-balance error', value: eErr, unit: '-', status: eErr < 1e-3 ? 'ok' : 'warn', note: 'max |E − E0 − work of applied forces| / E0' },
        { key: 'quat_err', label: 'Quaternion norm error', value: qErr, unit: '-', status: qErr < 1e-6 ? 'ok' : 'warn' },
        { key: 'steps', label: 'Integration steps', value: sim.steps, unit: '-' },
        { key: 'alt_change_m', label: 'Altitude change', value: -xe[2] + x0[2], unit: 'm' }, { key: 'speed_end_ms', label: 'Final speed', value: Vt[n - 1], unit: 'm/s' },
        { key: 'nz_max', label: 'Peak normal load factor', value: N.amax(nz), unit: 'g' }, { key: 'nz_min', label: 'Lowest normal load factor', value: N.amin(nz), unit: 'g' },
        { key: 'bank_max_deg', label: 'Largest bank angle', value: deg(N.amax(eu.map((e) => Math.abs(e[0])))), unit: 'deg' }, { key: 'pitch_max_deg', label: 'Largest pitch attitude', value: deg(N.amax(eu.map((e) => Math.abs(e[1])))), unit: 'deg' },
        { key: 'range_m', label: 'Ground distance covered', value: Math.hypot(xe[0], xe[1]), unit: 'm' },
      ],
      plots: [
        { type: 'line', title: 'Attitude', xlabel: 'Time [s]', ylabel: 'Euler angle [deg]', series: [{ name: 'Bank φ', x: g(T), y: g(eu.map((e) => deg(e[0]))) }, { name: 'Pitch θ', x: g(T), y: g(eu.map((e) => deg(e[1]))) }, { name: 'Heading ψ', x: g(T), y: g(eu.map((e) => deg(e[2]))) }] },
        { type: 'line', title: 'Angular rates', xlabel: 'Time [s]', ylabel: 'Rate [deg/s]', series: [{ name: 'p', x: g(T), y: g(Y.map((x) => deg(x[6]))) }, { name: 'q', x: g(T), y: g(Y.map((x) => deg(x[7]))) }, { name: 'r', x: g(T), y: g(Y.map((x) => deg(x[8]))) }] },
        { type: 'line', title: 'Speed', xlabel: 'Time [s]', ylabel: 'Speed [m/s]', series: [{ name: 'True airspeed', x: g(T), y: g(Vt) }] },
        { type: 'line', title: 'Altitude', xlabel: 'Time [s]', ylabel: 'Altitude [m]', series: [{ name: 'Altitude', x: g(T), y: g(Y.map((x) => -x[2])) }] },
        { type: 'line', title: 'Ground track', xlabel: 'East [m]', ylabel: 'North [m]', series: [{ name: 'Track', x: g(Y.map((x) => x[1])), y: g(Y.map((x) => x[0])) }] },
        { type: 'line', title: 'Normal load factor', xlabel: 'Time [s]', ylabel: 'Load factor [g]', series: [{ name: 'nz', x: g(T), y: g(nz) }] },
      ],
      outputs: { V0_ms: V0 },
      warnings, models: ['Non-linear six-degree-of-freedom rigid-body model', 'Quaternion attitude with norm-drift correction', 'Dormand–Prince RK45 adaptive integration', i.rotary ? 'Hover-derivative force and moment model' : 'Derivative-based aerodynamic model with lift saturation, parabolic polar and ISA density'],
      assumptions: ['Flat, non-rotating Earth; no wind', 'Constant mass and inertia', 'Thrust along the body x-axis through the CG; no engine dynamics or gyroscopic moments', 'Controls move instantly (actuators are in Suite 16)', 'Aerodynamic and hover derivatives come from the preliminary build-ups, whose empirical factors are typical values rather than data for a specific aircraft'],
    };
  },
  verify() {
    // torque-free tumbling of an asymmetric body: kinetic energy, angular momentum and quaternion norm are conserved
    const I = [[2, 0, 0], [0, 3, 0], [0, 0, 5]], f = rigidBody(1, I, N.inv(I), () => ({ F: [0, 0, 0], M: [0, 0, 0] }), 0), x0 = [0, 0, 0, 0, 0, 0, 1.2, 0.05, 0.8, 1, 0, 0, 0, 0], r = N.rk45(f, 0, x0, 20, { rtol: 1e-10, atol: 1e-12 }), x = r.y[r.y.length - 1];
    const H = (y) => { const w = [y[6], y[7], y[8]], Rm = dcm([y[9], y[10], y[11], y[12]]); return N.matvec(Rm, N.matvec(I, w)); }, T = (y) => 0.5 * (2 * y[6] ** 2 + 3 * y[7] ** 2 + 5 * y[8] ** 2), H0 = H(x0), H1 = H(x);
    // free fall from rest: z = g·t²/2
    const ff = N.rk45(rigidBody(3, I, N.inv(I), () => ({ F: [0, 0, 0], M: [0, 0, 0] })), 0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], 3, { rtol: 1e-10 });
    const b = Object.fromEntries(sixdof.inputs.map((q) => [q.key, q.default])), o = N.kv(sixdof.run({ ...b, manoeuvre: 'none', t_end: 20 }));
    return [N.check('Torque-free rotation conserves kinetic energy', T(x), T(x0), 1e-7, 'Euler equations'), N.check('Inertial angular momentum is conserved', Math.hypot(H1[0] - H0[0], H1[1] - H0[1], H1[2] - H0[2]) / N.norm(H0), 0, 1e-6, 'Euler equations with quaternion kinematics'),
      N.check('Quaternion norm stays at 1 without correction', Math.hypot(x[9], x[10], x[11], x[12]), 1, 1e-7, 'Kinematic invariant'), N.check('Free fall z = g·t²/2', ff.y[ff.y.length - 1][2], 0.5 * G0 * 9, 1e-8, 'Newton'),
      N.check('Trimmed flight holds altitude with controls fixed', o.alt_change_m / 1000 + 1, 1, 2e-3, 'Equilibrium of the non-linear model'), N.check('Energy audit of the full aircraft model', o.energy_err, 0, 1e-4, 'Work–energy theorem')];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (o.energy_err > 1e-3 || o.quat_err > 1e-6) out.push({ severity: 'warn', title: 'Integration accuracy is marginal', detail: `Energy error ${o.energy_err.toExponential(1)}, quaternion error ${o.quat_err.toExponential(1)}.`, action: 'Tighten the integrator tolerance.', basis: 'Conservation checks' });
    if (o.nz_max > 2.5 || o.nz_min < 0) out.push({ severity: 'advise', title: 'Large load-factor excursion', detail: `nz between ${o.nz_min.toFixed(2)} and ${o.nz_max.toFixed(2)} g.`, action: 'Compare with the V–n envelope of Suite 5 and reduce the control amplitude or add envelope protection (Suite 16).', basis: 'Manoeuvre envelope' });
    return out;
  },
};

// ---- 6. turbulence response -------------------------------------------------------------------
const turb = {
  id: 'turbulence', title: 'Response to atmospheric turbulence (Dryden and von Kármán)', fidelity: 'reduced-order',
  summary: 'Vertical-gust response of the longitudinal model: RMS normal acceleration and pitch motion from the Dryden shaping filter (exact covariance and a seeded time history) and from the von Kármán spectrum.',
  equations: ['Longitudinal linearised state equations', 'Stochastic gust-response equations', 'Aerodynamic force and moment equations'],
  applicable: hasWing,
  inputs: [...FW,
    { key: 'sigma_w', label: 'RMS vertical gust velocity', unit: 'm/s', default: 1.5, min: 0, max: 20, group: 'Turbulence', help: 'An operating (ride-quality) level, not a design gust intensity: about 0.5–1.5 light, 1.5–3 moderate, 3–6 severe (typical)' },
    { key: 'L_w', label: 'Turbulence scale length', unit: 'm', default: 533, min: 5, max: 5000, group: 'Turbulence', help: 'Typical values above about 600 m: 533 m (1750 ft) for the Dryden form; 762 m (2500 ft) for the von Kármán form, the scale also set by 14 CFR 25.341(b)' },
    { key: 't_sim', label: 'Time-history length', unit: 's', default: 120, min: 10, max: 1200, group: 'Turbulence' }, { key: 'seed', label: 'Random seed', unit: '', default: 2024, min: 1, max: 1e9, step: 1, discrete: true, group: 'Turbulence' },
    { key: 'nFreq', label: 'Frequency points', unit: '', default: 300, min: 50, max: 4000, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => ({ ...fwDefaults(c, up, d), sigma_w: Math.max(0.5, c.atm.turb_intensity * c.flight.V_ms) }),
  run(i) {
    const s = longSS(i), A = s.A, u0 = i.V, tau = i.L_w / u0, K = i.sigma_w * Math.sqrt(tau), warnings = [];
    // augmented system: aircraft + Dryden filter states (g1, g2), white-noise input of unit intensity
    const Aa = N.zeros(6), cw = [K / (tau * tau), (K * Math.sqrt(3)) / tau];
    for (let p = 0; p < 4; p++) { for (let q = 0; q < 4; q++) Aa[p][q] = A[p][q]; Aa[p][4] = -A[p][1] * cw[0]; Aa[p][5] = -A[p][1] * cw[1]; }
    Aa[4][5] = 1; Aa[5][4] = -1 / (tau * tau); Aa[5][5] = -2 / tau;
    const Cn = [A[1][0], A[1][1], A[1][2] - u0, A[1][3], -A[1][1] * cw[0], -A[1][1] * cw[1]].map((v) => -v / G0), Q = N.zeros(6); Q[5][5] = 1;
    const stable = N.amax(N.eig(A).map((l) => l[0])) < 0; let sn = NaN, sq = NaN, sth = NaN;
    if (stable) { const P = N.lyap(Aa, Q), v = (c) => Math.sqrt(Math.max(N.dot(c, N.matvec(P, c)), 0)); sn = v(Cn); sq = v([0, 0, 1, 0, 0, 0]); sth = v([0, 0, 0, 1, 0, 0]); } else warnings.push('The open-loop aircraft has an unstable root, so the stationary response variance is unbounded; only the time history is meaningful.');
    // von Kármán: frequency-domain integration of |H(iω)|²·Φ(ω)
    const om = N.logspace(1e-3 / tau, 3e2 / tau, Math.round(i.nFreq)), bg = A.map((r) => -r[1]);
    const H = om.map((w) => { let x; try { x = N.csolve(A.map((r, p) => r.map((v, q) => [-v, p === q ? w : 0])), bg); } catch { return 0; } const az = N.C.sub(N.C.add(N.C.add(N.C.scale(x[0], A[1][0]), N.C.scale(x[1], A[1][1])), N.C.add(N.C.scale(x[2], A[1][2] - u0), N.C.scale(x[3], A[1][3]))), [A[1][1], 0]); return N.C.abs(az) / G0; });
    const vk = (w) => { const y = 1.339 * tau * w; return ((i.sigma_w ** 2 * tau) / Math.PI) * (1 + (8 / 3) * y * y) / (1 + y * y) ** (11 / 6); }, dry = (w) => { const y = tau * w; return ((i.sigma_w ** 2 * tau) / Math.PI) * (1 + 3 * y * y) / (1 + y * y) ** 2; };
    const Pv = om.map((w, k) => H[k] ** 2 * vk(w)), Pd = om.map((w, k) => H[k] ** 2 * dry(w)), snV = Math.sqrt(N.trapz(om, Pv)), snD = Math.sqrt(N.trapz(om, Pd)), N0 = Math.sqrt(N.trapz(om, om.map((w, k) => w * w * Pv[k])) / (snV * snV || 1)) / TAU;
    // seeded time history (RK4 with zero-order-hold white noise)
    const dt = Math.min(0.02, tau / 20), ns = Math.min(20000, Math.round(i.t_sim / dt)), u = N.rng(i.seed); let x = new Array(6).fill(0), nk = 0; const f = (xx) => { const d = N.matvec(Aa, xx); d[5] += nk; return d; }, th = [], nzs = [], wg = [];
    for (let k = 0; k < ns; k++) { nk = N.randn(u) / Math.sqrt(dt); const k1 = f(x), k2 = f(N.vadd(x, k1, dt / 2)), k3 = f(N.vadd(x, k2, dt / 2)), k4 = f(N.vadd(x, k3, dt)); x = x.map((v, j) => v + (dt / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])); if (!fin(x[0]) || Math.abs(x[0]) > 1e6) break; th.push(k * dt); nzs.push(N.dot(Cn, x)); wg.push(cw[0] * x[4] + cw[1] * x[5]); }
    const snT = N.std(nzs), ref = fin(sn) ? sn : snT;
    if (ns * dt < i.t_sim - 1e-9) warnings.push(`The time history is limited to ${ns} steps: ${(ns * dt).toFixed(0)} s of the requested ${i.t_sim} s were simulated.`);
    return {
      kpis: [
        { key: 'sigma_nz', label: 'RMS normal acceleration (Dryden, exact covariance)', value: sn, unit: 'g', status: !fin(sn) || sn < 0.2 ? 'ok' : 'warn' },
        { key: 'sigma_nz_vk', label: 'RMS normal acceleration (von Kármán)', value: snV, unit: 'g' }, { key: 'sigma_nz_dryden_freq', label: 'RMS normal acceleration (Dryden, spectrum integration)', value: snD, unit: 'g' },
        { key: 'sigma_nz_sim', label: 'RMS normal acceleration (time history)', value: snT, unit: 'g' }, { key: 'nz_3sigma', label: '3σ normal-acceleration increment', value: 3 * ref, unit: 'g' },
        { key: 'A_bar', label: 'Gust response factor σn/σw', value: i.sigma_w ? ref / i.sigma_w : 0, unit: 'g/(m/s)' }, { key: 'N0_Hz', label: 'Characteristic frequency N0', value: N0, unit: 'Hz' },
        { key: 'sigma_q_dps', label: 'RMS pitch rate', value: deg(sq), unit: 'deg/s' }, { key: 'sigma_theta_deg', label: 'RMS pitch attitude', value: deg(sth), unit: 'deg' },
        { key: 'sigma_w_sim', label: 'RMS of the generated gust', value: N.std(wg), unit: 'm/s', note: `Target ${i.sigma_w}` },
      ].filter((k) => fin(k.value) || k.key === 'sigma_nz'),
      plots: [
        { type: 'line', title: 'Normal-acceleration spectrum', xlabel: 'Frequency [rad/s]', ylabel: 'PSD [g²/(rad/s)]', xlog: true, ylog: true, series: [{ name: 'von Kármán', x: om, y: Pv.map((v) => Math.max(v, 1e-30)) }, { name: 'Dryden', x: om, y: Pd.map((v) => Math.max(v, 1e-30)) }] },
        { type: 'line', title: 'Gust velocity spectra', xlabel: 'Frequency [rad/s]', ylabel: 'PSD [(m/s)²/(rad/s)]', xlog: true, ylog: true, series: [{ name: 'von Kármán', x: om, y: om.map((w) => Math.max(vk(w), 1e-30)) }, { name: 'Dryden', x: om, y: om.map((w) => Math.max(dry(w), 1e-30)) }] },
        { type: 'line', title: 'Simulated turbulence encounter', xlabel: 'Time [s]', ylabel: 'Load-factor increment [g]', series: [{ name: 'Δnz', x: thin(th, 400), y: thin(nzs, 400) }] },
        { type: 'line', title: 'Generated vertical gust', xlabel: 'Time [s]', ylabel: 'Gust velocity [m/s]', series: [{ name: 'w gust', x: thin(th, 400), y: thin(wg, 400) }] },
      ],
      warnings, models: ['Dryden turbulence model: second-order shaping filter driven by seeded white noise', 'von Kármán turbulence model: spectrum integration through the aircraft transfer function', 'Lyapunov covariance equation for the exact stationary response', 'Rice characteristic frequency'],
      assumptions: ['Vertical gust only, uniform over the span and without pitch-gust (penetration) effect', 'Quasi-steady aerodynamics, rigid aircraft: wing bending and unsteady lift build-up (Suite 3) are not included', 'Frozen, isotropic, Gaussian turbulence', 'The default gust RMS and scale are typical operating values for ride quality, not certification design intensities'],
    };
  },
  convergence: { param: 'nFreq', label: 'Frequency points', levels: [75, 150, 300, 600], metric: 'sigma_nz_vk' },
  verify() {
    const tau = 4, sg = 2, K = sg * Math.sqrt(tau), P = N.lyap([[0, 1], [-1 / 16, -0.5]], [[0, 0], [0, 1]]), c = [K / 16, (K * Math.sqrt(3)) / tau], b = Object.fromEntries(turb.inputs.map((f) => [f.key, f.default])), o = N.kv(turb.run({ ...b, nFreq: 1200 }));
    return [N.check('Dryden filter output variance equals σ²', N.dot(c, N.matvec(P, c)), sg * sg, 1e-9, 'Lyapunov equation for the shaping filter'),
      N.check('Dryden: covariance and spectrum integration agree', o.sigma_nz_dryden_freq, o.sigma_nz, 0.01, 'Parseval / Wiener–Khinchin'),
      N.check('Seeded time history reproduces the gust RMS', o.sigma_w_sim, 1.5, 0.12, 'Sampling scatter of a 120 s record')];
  },
  recommend(res) {
    const o = res.outputs, s = fin(o.sigma_nz) ? o.sigma_nz : o.sigma_nz_sim, out = [];
    if (s > 0.2) out.push({ severity: 'advise', title: 'Rough ride in this turbulence', detail: `RMS normal acceleration ${s.toFixed(2)} g (3σ ${(3 * s).toFixed(2)} g).`, action: 'A higher wing loading, lower lift-curve slope (sweep) or active gust-load alleviation (Suite 16) reduces the response; feed the exceedance spectrum to Suite 9.', basis: 'Power-spectral gust response' });
    return out;
  },
};

// ---- 7. hover dynamics ------------------------------------------------------------------------
const hover = {
  id: 'hover', title: 'Hover dynamics of helicopters and multirotors', fidelity: 'reduced-order',
  summary: 'Linearised hover model: speed and attitude derivatives, the unstable pitch and roll oscillations, heave and yaw damping, control power and the open-loop response to a cyclic pulse.',
  equations: ['Longitudinal and lateral-directional linearised state equations', 'Characteristic eigenvalue equations', 'Stability derivative equations', 'Newton–Euler rigid-body equations'],
  applicable: hasRotor, inputs: [...HOV, { key: 'pulse_pct', label: 'Control pulse (share of 1 rad cyclic or of full differential thrust)', unit: '%', default: 1, min: 0.01, max: 20, group: 'Manoeuvre' }],
  defaults: hovDefaults,
  run(i) {
    const h = hoverDerivs(i), el = N.eig(h.Alon).map(modeOf), ea = N.eig(h.Alat).map(modeOf), osc = (e) => e.find((m) => m.im > 1e-9), lo = osc(el), la = osc(ea), warnings = [];
    const t2 = (m) => (m && m.re > 0 ? Math.LN2 / m.re : NaN), tEnd = lo ? Math.min(2.5 * (TAU / lo.im), 40) : 10, r = linSim(h.Alon, h.Blon, (t) => [t > 0.5 && t < 1.5 ? i.pulse_pct / 100 : 0, 0], tEnd, simSteps(el, tEnd, 1500));
    if (lo && lo.re > 0) warnings.push(`The hover pitch oscillation is unstable (period ${(TAU / lo.im).toFixed(1)} s, doubles in ${t2(lo).toFixed(1)} s): normal for unaugmented rotorcraft, it must be stabilised by the pilot or by attitude feedback.`);
    if (h.CT / h.sig > 0.16) warnings.push('Blade loading CT/σ is above 0.16: the rotor is close to stall and the linear derivatives are optimistic.');
    if (i.kind === 'multirotor') warnings.push('Multirotor derivatives use a linear rotor-drag model and fixed-pitch thrust sensitivity; motor and ESC dynamics are in Suite 16.');
    const gMu = G0 * h.Mu, cub = Math.cbrt(Math.max(gMu, 0));
    return {
      kpis: [
        { key: 'hover_pitch_period_s', label: 'Pitch oscillation period', value: lo ? TAU / lo.im : NaN, unit: 's' }, { key: 'hover_pitch_T2_s', label: 'Pitch oscillation time to double', value: t2(lo), unit: 's', status: !(t2(lo) < 3) || i.kind === 'multirotor' ? 'ok' : 'warn', note: i.kind === 'multirotor' ? 'An unstable open loop is normal for multirotors: the rate and attitude loops close it' : 'Faster than about 3 s is hard to fly without augmentation (typical piloting experience, not a sourced limit)' },
        { key: 'hover_roll_period_s', label: 'Roll oscillation period', value: la ? TAU / la.im : NaN, unit: 's' }, { key: 'hover_roll_T2_s', label: 'Roll oscillation time to double', value: t2(la), unit: 's' },
        { key: 'Xu', label: 'Speed damping Xu', value: h.Xu, unit: '1/s' }, { key: 'Mu', label: 'Speed stability Mu', value: h.Mu, unit: 'rad/(s·m)' }, { key: 'Mq', label: 'Pitch damping Mq', value: h.Mq, unit: '1/s' }, { key: 'Lp', label: 'Roll damping Lp', value: h.Lp, unit: '1/s' },
        { key: 'Zw', label: 'Heave damping Zw', value: h.Zw, unit: '1/s' }, { key: 'heave_tau_s', label: 'Heave time constant', value: -1 / h.Zw, unit: 's' }, { key: 'Nr', label: 'Yaw damping Nr', value: h.Nr, unit: '1/s' },
        { key: 'pitch_ctrl_power', label: 'Pitch control power', value: h.Mc, unit: 'rad/s² per unit control' }, { key: 'pitch_bandwidth_rads', label: 'Pitch rate-response bandwidth −Mq', value: -h.Mq, unit: 'rad/s' },
        { key: 'CT', label: 'Thrust coefficient', value: h.CT, unit: '-' }, { key: 'blade_loading', label: 'Blade loading CT/σ', value: h.CT / h.sig, unit: '-', status: h.CT / h.sig < 0.16 ? 'ok' : 'warn' },
        { key: 'period_limit_s', label: 'Period in the limit of no damping, 2π/(√3/2·(g·Mu)^⅓)', value: cub ? TAU / ((Math.sqrt(3) / 2) * cub) : NaN, unit: 's' },
      ].filter((k) => fin(k.value)),
      plots: [
        { type: 'line', title: 'Hover roots', xlabel: 'Real part [1/s]', ylabel: 'Imaginary part [rad/s]', series: [{ name: 'Longitudinal', x: el.map((m) => m.re), y: el.map((m) => m.im), style: 'points' }, { name: 'Lateral-directional', x: ea.map((m) => m.re), y: ea.map((m) => m.im), style: 'points' }], annotations: [{ x: 0, label: 'Stability boundary' }] },
        { type: 'line', title: 'Open-loop response to a 1 s pitch-control pulse', xlabel: 'Time [s]', ylabel: 'Pitch attitude [deg]', series: [{ name: 'θ', x: thin(r.t), y: thin(r.x.map((x) => N.clamp(deg(x[3]), -90, 90))) }] },
        { type: 'line', title: 'Speed response', xlabel: 'Time [s]', ylabel: 'Forward speed [m/s]', series: [{ name: 'u', x: thin(r.t), y: thin(r.x.map((x) => N.clamp(x[0], -100, 100))) }] },
      ],
      tables: [{ title: 'Hover eigenvalues', columns: ['Axis', 'Real [1/s]', 'Imaginary [rad/s]'], rows: [...el.map((m) => ['Longitudinal', m.re, m.im]), ...ea.map((m) => ['Lateral-directional', m.re, m.im])] }],
      outputs: { A_hover_long: h.Alon, B_hover_long: h.Blon, A_hover_lat: h.Alat, B_hover_lat: h.Blat },
      warnings, models: ['Helicopter rotorcraft flight dynamics model: quasi-steady flapping (tip-path-plane) derivatives in hover', 'Momentum and blade-element thrust sensitivity for heave damping', 'Multirotor rigid-body model with differential-thrust moments and rotor-drag speed damping'],
      assumptions: ['Hover, small perturbations; rotor flapping and inflow respond instantly', 'Longitudinal and lateral motions uncoupled; no fuselage or tailplane aerodynamics', 'Helicopter tail rotor represented by a thrust sensitivity at a nominal 210 m/s tip speed and 0.15 solidity', 'Induced-power factor 1.15, blade profile drag coefficient 0.01 and the multirotor in-plane drag damping are typical values, not data for a specific rotor'],
    };
  },
  calibration: { params: [{ key: 'k_drag', min: 0, max: 2 }, { key: 'lock', min: 1, max: 16 }, { key: 'h_rotor', min: 0, max: 5 }], sweep: 'mass_kg', target: 'hover_pitch_period_s', note: 'Flight-identified hover derivatives (frequency sweeps) or the measured period of the hover oscillation' },
  verify() {
    const Mu = 0.02, A = [[0, 0, 0, -G0], [0, -0.3, 0, 0], [Mu, 0, 0, 0], [0, 0, 1, 0]], o = N.eig(A).find((l) => l[1] > 1e-9), c = Math.cbrt(G0 * Mu);
    const h = hoverDerivs(Object.fromEntries(HOV.map((f) => [f.key, f.default])));
    return [N.check('Undamped hover oscillation frequency (√3/2)·(g·Mu)^⅓', o[1], (Math.sqrt(3) / 2) * c, 1e-9, 'Roots of s³ + g·Mu = 0'), N.check('…and its growth rate ½·(g·Mu)^⅓', o[0], 0.5 * c, 1e-9, 'Roots of s³ + g·Mu = 0'),
      N.check('Hover inflow ratio λ = sqrt(CT/2)', h.lam, Math.sqrt(h.CT / 2), 1e-12, 'Momentum theory')];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], multi = i.kind === 'multirotor';
    if (o.hover_pitch_T2_s < 3 || o.hover_roll_T2_s < 3) out.push({ severity: multi ? 'info' : 'warn', title: multi ? 'Hover is unstable without attitude feedback' : 'Fast hover instability', detail: `Time to double: pitch ${fin(o.hover_pitch_T2_s) ? o.hover_pitch_T2_s.toFixed(1) : '—'} s, roll ${fin(o.hover_roll_T2_s) ? o.hover_roll_T2_s.toFixed(1) : '—'} s.`, action: multi ? 'Normal for a multirotor, which is always flown through rate and attitude loops: design them in Suite 16 (hover cascade) with a bandwidth well above the unstable root.' : 'Provide attitude or rate stabilisation (Suite 16); a stabiliser bar, larger hinge offset or lower rotor height above the CG also help.', basis: 'Hover speed–attitude oscillation' });
    if (o.heave_tau_s > 5) out.push({ severity: 'advise', title: 'Slow heave response', detail: `Heave time constant ${o.heave_tau_s.toFixed(1)} s.`, action: 'Expect altitude drift in gusts; add altitude hold.', basis: 'Heave damping Zw' });
    return out;
  },
};

export default {
  id: 'flightdyn', n: 4,
  tagline: 'Is the aircraft stable, how does it respond to the controls and to turbulence, and does it fly the way a pilot or autopilot needs?',
  analyses: [derivAn, trimAn, longAn, latAn, sixdof, turb, hover],
  consumes: [
    { from: 'cfd', keys: ['CD0', 'e_oswald', 'CLmax'], why: 'Drag polar and lift limit' },
    { from: 'propulsion', keys: ['thrust_static_N'], why: 'Installed thrust for trim and simulation' }, { from: 'propeller', keys: ['eta_prop'], why: 'Propulsive efficiency' },
  ],
  provides: [
    { key: 'static_margin', label: 'Static margin', unit: 'MAC' }, { key: 'sp_omega_rads', label: 'Short-period frequency', unit: 'rad/s' }, { key: 'sp_zeta', label: 'Short-period damping', unit: '-' },
    { key: 'ph_period_s', label: 'Phugoid period', unit: 's' }, { key: 'ph_zeta', label: 'Phugoid damping', unit: '-' }, { key: 'dr_omega_rads', label: 'Dutch-roll frequency', unit: 'rad/s' }, { key: 'dr_zeta', label: 'Dutch-roll damping', unit: '-' },
    { key: 'roll_tau_s', label: 'Roll time constant', unit: 's' }, { key: 'spiral_T2_s', label: 'Spiral time to double', unit: 's' }, { key: 'trim_alpha_deg', label: 'Trim angle of attack', unit: 'deg' }, { key: 'trim_de_deg', label: 'Trim elevator', unit: 'deg' },
  ],
  handoff: [
    { model: 'Aerodynamic coefficient lookup models from CFD or wind-tunnel databases', why: 'The native model uses constant derivatives from an empirical build-up; non-linear tables need test or CFD data', tool: 'Aerodynamic database from CFD / wind tunnel, loaded into a flight-simulation framework' },
    { model: 'Non-linear stall, post-stall and spin models', why: 'Depend on measured high-angle-of-attack and rotary-balance data', tool: 'Flight-test-validated simulation with high-alpha aerodynamic data' },
    { model: 'Comprehensive helicopter flight dynamics (blade-element rotor with flapping and inflow dynamics, forward flight and transition)', why: 'Needs individual-blade or tip-path-plane dynamics with dynamic inflow and interference models', tool: 'Rotorcraft flight-dynamics code' },
    { model: 'Pilot-in-the-loop models and handling-qualities ratings', why: 'Require a piloted simulator or validated pilot models with task definitions', tool: 'Piloted simulation facility' },
    { model: 'Flight dynamics coupled with structural flexibility and with ground contact', why: 'Flexible-aircraft and landing dynamics are separate coupled models (Suites 3 and 14)', tool: 'Coupled flight-loads simulation' },
  ],
};
