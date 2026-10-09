// Suite 8 — Propeller and Rotor Performance.
// Blade-element momentum theory in axial flight (static to cruise) with Prandtl tip and hub loss on a generic
// section polar, performance maps and constant-speed schedules, minimum-induced-loss design (Adkins–Liebeck),
// actuator-disk and ducted-fan momentum limits, installation effects and electric motor–propeller matching.

import * as N from '../core/numerics.js';
import { isa, G0 } from '../core/atmosphere.js';

// ---- shared physics -------------------------------------------------------------------------
const PROP_T = 'propeller (generic planform, helical twist)', ROTOR_T = 'rotor (constant chord, linear twist)';
const hasProp = (c) => (c.prop.prop_dia_m > 0 || c.rotor.R_m > 0 ? true : 'No propeller or rotor is defined for this vehicle (turbofan or turbojet). The fan is part of the engine cycle in Suite 7.');
const nProps = (c) => (c.meta.type === 'helicopter' ? 1 : Math.max(1, c.prop.n_eng));
const fmt = (v, n = 2) => (Number.isFinite(v) ? v.toFixed(n) : 'n/a');
/** Drop KPIs flagged opt that do not exist for this condition (e.g. propulsive efficiency at rest). */
const keep = (list) => list.filter((k) => !(k.opt && !Number.isFinite(k.value))).map(({ opt, ...k }) => k);

/**
 * Generic parametric section polar (not a specific aerofoil): linear lift with Prandtl–Glauert slope, smooth
 * saturation at clmax, flat-plate behaviour beyond about 23°, parabolic profile drag with Reynolds scaling,
 * separation drag and a Lock-type fourth-power drag rise above the critical Mach number.
 * Returns a shared scratch array [cl, cd] that must be read immediately.
 */
function polar(q, al, Re, M) {
  if (q.ideal) { POL[0] = q.cl0 + q.cla * al; POL[1] = 0; return POL; }
  const Mc = q.mach ? Math.min(M, 0.85) : 0, lin = q.cl0 + (q.cla * al) / Math.sqrt(1 - Mc * Mc), att = q.clmax * Math.tanh(lin / q.clmax);
  const w = 1 / (1 + Math.exp((Math.abs(al) - 0.4) / 0.05)), s = Math.sin(al), cl = w * att + (1 - w) * 1.1 * Math.sin(2 * al);
  const ex = Math.max(0, Math.abs(lin) - 0.85 * q.clmax), fRe = N.clamp((1e6 / Math.max(Re, 1e3)) ** 0.2, 0.6, 5) * (Re < 1e5 ? (1e5 / Math.max(Re, 1e3)) ** 0.3 : 1);
  const dM = q.mach ? M - (q.Mcrit - 0.1 * Math.abs(cl)) : 0, d = att - q.cl0;
  POL[0] = cl; POL[1] = q.cdmin * fRe + 0.01 * d * d + 0.12 * ex * ex + (1 - w) * 2 * s * s + (dM > 0 ? 20 * dM ** 4 : 0);
  return POL;
}
const POL = [0, 0];
const polarOf = (i, extra = {}) => ({ cla: i.cla ?? 5.9, cl0: i.cl0 ?? 0.3, clmax: i.clmax ?? 1.3, cdmin: i.cdmin ?? 0.008, Mcrit: i.Mcrit ?? 0.75, mach: true, ideal: false, ...extra });

/** Blade geometry at mid-point stations: chord [m] and twist relative to the 75% station [rad]. */
function bladeGeom(i) {
  const R = i.D / 2, n = Math.round(i.nR), xh = N.clamp(i.hub_ratio, 0.02, 0.6), dx = (1 - xh) / n, rotor = i.blade_type === ROTOR_T, x = [], c = [], tw = [];
  const sh = (v) => Math.sqrt(v * (1.1 - v)), lt = Math.max(0.02, i.J_twist || 0.8) / Math.PI; // tan(helix angle) at the tip for the twist design advance ratio
  for (let j = 0; j < n; j++) {
    const v = xh + (j + 0.5) * dx; x.push(v);
    c.push(i.c75_R * R * (rotor ? 1 : sh(v) / sh(0.75)));
    tw.push(rotor ? N.rad(i.twist_deg) * (v - 0.75) : Math.atan(lt / v) - Math.atan(lt / 0.75));
  }
  return { R, Nb: Math.max(1, Math.round(i.Nb)), xh, dx, x, c, tw, tipLoss: i.tip_loss !== false };
}
/**
 * Blade-element momentum solution in axial flight. Unknown per annulus is the inflow angle φ; the residual
 * V·(sinφ·cosφ + σ′ct/4F) − Ωr·(sin²φ − σ′cn/4F) = 0 combines axial and angular momentum with the blade forces
 * and stays regular at V = 0. op = {V, Om, rho, a, nu}; b75 is the pitch at 75% radius [rad].
 */
function bemt(g, q, op, b75) {
  const n = g.x.length, R = g.R, o = { T: 0, Q: 0, Tm: 0, Qm: 0, Pax: 0, Psw: 0, Ppr: 0, al: [], cl: [], dT: [], dQ: [], va: [], vt: [], phi: [], clMax: -9, Mmax: 0, failed: 0 };
  for (let j = 0; j < n; j++) {
    const x = g.x[j], r = x * R, c = g.c[j], b = b75 + g.tw[j], sig = (g.Nb * c) / (2 * Math.PI * r), Ur = op.Om * r, Wg = Math.hypot(op.V, Ur), Re = (Wg * c) / op.nu, M = Wg / op.a;
    let cn = 0, ct = 0, F = 1, cl = 0, cd = 0;
    const res = (ph) => {
      const s = Math.sin(ph), cs = Math.cos(ph), p = polar(q, b - ph, Re, M); cl = p[0]; cd = p[1]; cn = cl * cs - cd * s; ct = cl * s + cd * cs;
      F = 1;
      if (g.tipLoss) { const sa = Math.max(Math.abs(s), 1e-6); F = (2 / Math.PI) * Math.acos(Math.min(1, Math.exp((-g.Nb * (1 - x)) / (2 * x * sa)))) * (2 / Math.PI) * Math.acos(Math.min(1, Math.exp((-g.Nb * (x - g.xh)) / (2 * g.xh * sa)))); if (F < 1e-4) F = 1e-4; }
      return op.V * (s * cs + (sig * ct) / (4 * F)) - Ur * (s * s - (sig * cn) / (4 * F));
    };
    const ph = N.findRoot(res, 1e-4, Math.PI / 2 - 1e-4, 18, 1e-11);
    if (!Number.isFinite(ph)) { o.failed++; o.al.push(NaN); o.cl.push(0); o.dT.push(0); o.dQ.push(0); o.va.push(0); o.vt.push(0); o.phi.push(NaN); continue; }
    res(ph);
    const s = Math.sin(ph), cs = Math.cos(ph), W = Ur / Math.max(1e-6, cs + (sig * ct) / (4 * F * s)), qd = 0.5 * op.rho * W * W * g.Nb * c, dr = g.dx * R;
    const dT = qd * cn * dr, dFt = qd * ct * dr, va = W * s - op.V, vt = Ur - W * cs;
    o.T += dT; o.Q += dFt * r; o.Tm += 4 * Math.PI * r * op.rho * F * (op.V + va) * va * dr; o.Qm += 4 * Math.PI * r * r * op.rho * F * (op.V + va) * vt * dr;
    o.Pax += dT * va; o.Psw += dFt * vt; o.Ppr += qd * cd * dr * W;
    o.al.push(b - ph); o.cl.push(cl); o.dT.push(dT / dr); o.dQ.push((dFt * r) / dr); o.va.push(va); o.vt.push(vt); o.phi.push(ph);
    if (cl > o.clMax) o.clMax = cl; if (W / op.a > o.Mmax) o.Mmax = W / op.a;
  }
  o.P = o.Q * op.Om; o.eta = o.P > 0 && op.V > 0 ? (o.T * op.V) / o.P : 0; o.b75 = b75;
  return o;
}
/** Operating state from analysis inputs. */
function opOf(i, V = i.V, rpm = i.rpm) {
  const at = isa(i.alt_m ?? 0, i.dISA ?? 0), Om = (rpm * 2 * Math.PI) / 60, n = rpm / 60;
  return { V: Math.max(0, V), Om, n, rho: at.rho, a: at.a, nu: at.nu, D: i.D, A: (Math.PI * i.D * i.D) / 4, J: n > 0 ? V / (n * i.D) : 0, Mtip: Math.hypot(V, (Om * i.D) / 2) / at.a };
}
/** Pitch at 75% radius that gives the target thrust (key 'T') or power (key 'P'); NaN when unreachable. */
function solvePitch(g, q, op, key, target) {
  const p0 = Math.atan2(op.V, 0.75 * op.Om * g.R), f = (b) => bemt(g, q, op, b)[key] - target;
  return N.findRoot(f, p0 - N.rad(4), p0 + N.rad(32), 12, 1e-9);
}
const coeffs = (r, op) => ({ CT: r.T / (op.rho * op.n ** 2 * op.D ** 4), CP: r.P / (op.rho * op.n ** 3 * op.D ** 5), CQ: r.Q / (op.rho * op.n ** 2 * op.D ** 5) });
const idealP = (T, op) => (T > 0 ? T * (op.V / 2 + Math.sqrt((op.V * op.V) / 4 + T / (2 * op.rho * op.A))) : 0); // momentum-theory power T·(V + vi)

const GEOM = [
  { key: 'D', label: 'Propeller / rotor diameter', unit: 'm', default: 3.93, min: 0.05, group: 'Geometry' },
  { key: 'Nb', label: 'Blades', unit: '', default: 6, min: 1, max: 16, step: 1, discrete: true, group: 'Geometry' },
  { key: 'hub_ratio', label: 'Hub (spinner) radius ratio', unit: '-', default: 0.18, min: 0.02, max: 0.6, group: 'Geometry' },
  { key: 'blade_type', label: 'Blade planform and twist law', type: 'select', options: [PROP_T, ROTOR_T], default: PROP_T, group: 'Geometry', help: 'Generic distributions used when the real blade drawing is not available' },
  { key: 'c75_R', label: 'Chord at 75% radius / radius', unit: '-', default: 0.14, min: 0.01, max: 0.5, group: 'Geometry', help: '0.10–0.18 for propellers, 0.05–0.08 for helicopter rotors' },
  { key: 'J_twist', label: 'Advance ratio the helical twist is laid out for', unit: '-', default: 1.5, min: 0.05, max: 5, group: 'Geometry', help: 'Propeller law only: constant geometric pitch at this J' },
  { key: 'twist_deg', label: 'Linear twist, root to tip (rotor law)', unit: 'deg', default: -12, min: -45, max: 10, group: 'Geometry' },
  { key: 'rpm', label: 'Rotational speed', unit: 'rpm', default: 1200, min: 10, group: 'Operating point' },
];
const SECTION = [
  { key: 'cla', label: 'Section lift-curve slope', unit: '1/rad', default: 5.9, min: 3, max: 7, group: 'Blade section (generic polar)' },
  { key: 'cl0', label: 'Lift coefficient at zero incidence', unit: '-', default: 0.3, min: 0, max: 0.8, group: 'Blade section (generic polar)', help: 'Camber: 0 symmetric, 0.3–0.5 typical propeller sections' },
  { key: 'clmax', label: 'Maximum lift coefficient', unit: '-', default: 1.3, min: 0.6, max: 2, group: 'Blade section (generic polar)' },
  { key: 'cdmin', label: 'Minimum profile drag at Re = 10⁶', unit: '-', default: 0.008, min: 0, max: 0.05, group: 'Blade section (generic polar)', help: 'Scaled with Reynolds number inside the model' },
  { key: 'Mcrit', label: 'Critical Mach number at zero lift', unit: '-', default: 0.75, min: 0.5, max: 0.95, group: 'Blade section (generic polar)' },
];
const POINT = [
  { key: 'V', label: 'Flight speed (axial)', unit: 'm/s', default: 140, min: 0, max: 300, group: 'Operating point', help: '0 for static thrust or hover' },
  { key: 'alt_m', label: 'Altitude', unit: 'm', default: 6100, min: -500, max: 20000, group: 'Operating point' },
  { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -60, max: 50, group: 'Operating point' },
];
const NR = { key: 'nR', label: 'Radial stations', unit: '', default: 24, min: 6, max: 400, step: 1, discrete: true, group: 'Numerics' };
const T_REQ = { key: 'T_req', label: 'Required thrust per propeller', unit: 'N', default: 7400, min: 0, group: 'Operating point', help: 'Cruise drag (or weight in hover) divided by the number of propellers' };

function geomDefaults(c, d) {
  const useRotor = !(c.prop.prop_dia_m > 0) || d.isRotary, D = useRotor ? 2 * c.rotor.R_m : c.prop.prop_dia_m, rpm = (useRotor ? c.rotor.rpm : c.prop.rpm) || c.rotor.rpm || c.prop.rpm;
  const chordKnown = c.rotor.chord_m > 0 && c.rotor.R_m > 0 && (useRotor || Math.abs(2 * c.rotor.R_m - D) < 0.05 * D);
  return { D: D || undefined, Nb: (useRotor ? c.rotor.n_blades : c.prop.n_blades) || c.rotor.n_blades || undefined, rpm: rpm || undefined, blade_type: useRotor ? ROTOR_T : PROP_T, c75_R: chordKnown ? c.rotor.chord_m / c.rotor.R_m : undefined, twist_deg: c.rotor.twist_deg, hub_ratio: useRotor ? 0.12 : 0.18, cdmin: useRotor ? c.rotor.cd0 : undefined, cl0: useRotor ? 0.15 : undefined };
}
function pointDefaults(c, up, d) {
  const g = geomDefaults(c, d), rot = d.isRotary, alt = rot ? c.atm.alt_m : c.mission.cruise_alt_m || c.atm.alt_m, V = rot ? 0 : c.mission.cruise_V_ms || c.flight.V_ms, at = isa(alt, c.atm.dISA_K), n = nProps(c);
  let T = d.W / n;
  if (!rot) { const qS = 0.5 * at.rho * V * V * c.wing.S_m2, CL = d.W / qS; T = (qS * (c.aero.CD0 + d.k_induced * CL * CL)) / n; }
  // shaft power at the flight point: from the propulsion suite, otherwise the sea-level rating with a simple altitude lapse
  const sg = at.sigma, lapse = c.prop.type === 'electric' ? 1 : c.prop.type === 'piston' ? Math.max(0, (sg - 0.1325) / 0.8675) : sg ** 0.7;
  const Pav = (up.propulsion?.P_shaft_W > 0 ? up.propulsion.P_shaft_W : d.P_total * lapse) / n;
  return { ...g, V, alt_m: alt, dISA: c.atm.dISA_K, T_req: T, P_avail: Pav || undefined, J_twist: g.rpm && g.D && V > 0 ? N.clamp(V / ((g.rpm / 60) * g.D), 0.2, 4) : undefined, n_props: n };
}

// ---- (a) operating point --------------------------------------------------------------------
const point = {
  id: 'bemt', title: 'Blade-element momentum analysis at the operating point', fidelity: 'numerical',
  summary: 'Solves each blade annulus for its inflow angle so that blade forces match the momentum given to the air, with tip and hub losses, and sets blade pitch to deliver the required thrust or absorb the available power.',
  equations: ['Blade element momentum equations', 'Propeller thrust and torque equations', 'Induced velocity equations', 'Prandtl tip-loss correction (tip and hub)', 'Propulsive efficiency equations', 'Rotor figure-of-merit relations'],
  applicable: hasProp,
  inputs: [...GEOM, ...POINT,
    { key: 'mode', label: 'Blade pitch setting', type: 'select', options: ['match required thrust', 'absorb available power', 'fixed pitch'], default: 'match required thrust', group: 'Operating point' },
    T_REQ,
    { key: 'P_avail', label: 'Shaft power available per propeller', unit: 'W', default: 1.3e6, min: 0, group: 'Operating point', help: 'From the propulsion suite when available' },
    { key: 'beta75_deg', label: 'Blade pitch at 75% radius (fixed-pitch mode)', unit: 'deg', default: 35, min: -5, max: 85, group: 'Operating point' },
    { key: 'n_props', label: 'Number of propellers', unit: '', default: 2, min: 1, max: 32, step: 1, discrete: true, group: 'Operating point' },
    { key: 'tip_loss', label: 'Prandtl tip and hub loss', type: 'bool', default: true, group: 'Model' },
    ...SECTION, NR],
  defaults: (c, up, d) => pointDefaults(c, up, d),
  run(i) {
    const g = bladeGeom(i), q = polarOf(i), op = opOf(i), warnings = [], stat = op.V < 0.5;
    let b = i.mode === 'fixed pitch' ? N.rad(i.beta75_deg) : i.mode === 'absorb available power' ? solvePitch(g, q, op, 'P', i.P_avail) : solvePitch(g, q, op, 'T', i.T_req), matched = Number.isFinite(b);
    if (!matched) { b = N.goldenSection((t) => -bemt(g, q, op, t).T, Math.atan2(op.V, 0.75 * op.Om * g.R), Math.atan2(op.V, 0.75 * op.Om * g.R) + N.rad(30), 1e-5); warnings.push('The target cannot be reached at any blade pitch at this speed and rpm: the blades stall first. Results are shown at maximum thrust.'); }
    const r = bemt(g, q, op, b), k = coeffs(r, op), Pi = idealP(r.T, op), FMop = r.P > 0 ? Pi / r.P : 0;
    // static point at the same rpm: constant-speed propeller absorbing the available power, or the same pitch if no power is given
    let rs = r, staticNote = stat ? 'Same as the operating point' : i.P_avail > 0 && i.mode !== 'fixed pitch' ? 'Pitch reset to absorb the available power' : 'At the operating-point pitch';
    if (!stat) { const os = { ...op, V: 0, J: 0 }, bs = i.P_avail > 0 && i.mode !== 'fixed pitch' ? solvePitch(g, q, os, 'P', i.P_avail) : b; if (Number.isFinite(bs)) rs = bemt(g, q, os, bs); else { rs = bemt(g, q, os, N.goldenSection((t) => -bemt(g, q, os, t).T, N.rad(4), N.rad(32), 1e-4)); staticNote = 'Pitch for maximum static thrust (the available power cannot be absorbed at rest at this rpm)'; } }
    const FM = stat ? FMop : rs.P > 0 ? rs.T ** 1.5 / Math.sqrt(2 * op.rho * op.A) / rs.P : 0, etaI = op.V > 0 && r.T > 0 ? (r.T * op.V) / Pi : NaN;
    const AF = (1e5 / 16) * N.sum(g.x.map((x, j) => (g.c[j] / i.D) * x ** 3 * g.dx)), sol = N.sum(g.c.map((c) => (g.Nb * c * g.dx * g.R) / (Math.PI * g.R * g.R)));
    if (r.clMax > 0.92 * q.clmax) warnings.push(`Peak section lift coefficient ${r.clMax.toFixed(2)} is close to the stall value ${q.clmax.toFixed(2)}: part of the blade is stalled or about to stall.`);
    if (r.Mmax > q.Mcrit + 0.1) warnings.push(`Peak section Mach number ${r.Mmax.toFixed(2)} is well above the critical value: the generic drag-rise term dominates and efficiency and noise predictions are uncertain.`);
    if (r.failed) warnings.push(`${r.failed} annuli had no blade-element momentum solution (negative pitch or windmilling) and were left unloaded.`);
    if (i.P_avail > 0 && r.P > 1.02 * i.P_avail) warnings.push(`Power absorbed (${(r.P / 1e3).toFixed(1)} kW) exceeds the shaft power available (${(i.P_avail / 1e3).toFixed(1)} kW).`);
    if (stat) warnings.push('Static or hover condition: propulsive efficiency is undefined and is not published; use the figure of merit.');
    const kp = [
      { key: 'prop_thrust_N', label: 'Thrust per propeller', value: r.T, unit: 'N', status: matched ? 'ok' : 'bad' },
      { key: 'prop_power_W', label: 'Shaft power absorbed per propeller', value: r.P, unit: 'W', status: !(i.P_avail > 0) || r.P <= 1.02 * i.P_avail ? 'ok' : 'bad' },
      { key: 'prop_torque_Nm', label: 'Shaft torque', value: r.Q, unit: 'N m' },
      { key: 'CT_prop', label: 'Thrust coefficient T/(ρn²D⁴)', value: k.CT, unit: '-' },
      { key: 'CP_prop', label: 'Power coefficient P/(ρn³D⁵)', value: k.CP, unit: '-' },
      { key: 'J_adv', label: 'Advance ratio V/(nD)', value: op.J, unit: '-' },
      { key: 'beta75_op_deg', label: 'Blade pitch at 75% radius', value: N.deg(b), unit: 'deg' },
      { key: 'tip_mach', label: 'Helical tip Mach number', value: op.Mtip, unit: '-', status: op.Mtip < 0.85 ? 'ok' : op.Mtip < 0.95 ? 'warn' : 'bad', note: 'Propellers cruise at about 0.75–0.85; efficiency falls and noise rises quickly above' },
      { key: 'FM_prop', label: stat ? 'Figure of merit' : 'Static figure of merit at this rpm', value: FM, unit: '-', status: !stat ? undefined : FM > 0.6 ? 'ok' : 'warn', note: stat ? 'Lifting rotors reach 0.65–0.80; small low-Reynolds rotors 0.5–0.65' : 'For information: a propeller laid out for cruise is not judged on its static figure of merit' },
      { key: 'static_thrust_N', label: 'Static thrust per propeller at this rpm', value: rs.T, unit: 'N', note: staticNote },
      { key: 'eta_ideal', label: 'Ideal (actuator-disk) efficiency at this thrust', value: etaI, unit: '-', opt: true },
      { key: 'power_margin_pct', label: 'Shaft power margin', value: i.P_avail > 0 ? 100 * (1 - r.P / i.P_avail) : NaN, unit: '%', opt: true },
      { key: 'disk_loading_Pa', label: 'Disk loading', value: r.T / op.A, unit: 'N/m²' },
      { key: 'loss_axial_pct', label: 'Axial induced loss', value: (100 * r.Pax) / (r.P || 1), unit: '% of shaft power' },
      { key: 'loss_swirl_pct', label: 'Swirl loss', value: (100 * r.Psw) / (r.P || 1), unit: '% of shaft power' },
      { key: 'loss_profile_pct', label: 'Profile loss', value: (100 * r.Ppr) / (r.P || 1), unit: '% of shaft power' },
      { key: 'cl_peak', label: 'Peak section lift coefficient', value: r.clMax, unit: '-', status: r.clMax < 0.92 * q.clmax ? 'ok' : 'warn' },
      { key: 'activity_factor', label: 'Blade activity factor', value: AF, unit: '-', note: '80–150 per blade is typical for propellers' },
      { key: 'solidity', label: 'Solidity', value: sol, unit: '-' },
    ];
    if (!stat) kp.splice(2, 0, { key: 'eta_prop', label: 'Propulsive efficiency T·V/P', value: r.eta, unit: '-', status: r.eta > 0.75 ? 'ok' : r.eta > 0.6 ? 'warn' : 'bad' });
    return {
      kpis: keep(kp),
      plots: [
        { type: 'line', title: 'Radial thrust loading', xlabel: 'Radius r/R [-]', ylabel: 'Thrust per unit span, all blades [N/m]', series: [{ name: 'dT/dr', x: g.x, y: r.dT }] },
        { type: 'line', title: 'Radial torque loading', xlabel: 'Radius r/R [-]', ylabel: 'Torque per unit span [N m/m]', series: [{ name: 'dQ/dr', x: g.x, y: r.dQ }] },
        { type: 'line', title: 'Section incidence and pitch', xlabel: 'Radius r/R [-]', ylabel: 'Angle [deg]', series: [{ name: 'Incidence α', x: g.x, y: r.al.map(N.deg) }, { name: 'Blade pitch β', x: g.x, y: g.tw.map((t) => N.deg(b + t)) }, { name: 'Inflow angle φ', x: g.x, y: r.phi.map(N.deg), style: 'dash' }] },
        { type: 'line', title: 'Induced velocities at the disk', xlabel: 'Radius r/R [-]', ylabel: 'Velocity [m/s]', series: [{ name: 'Axial', x: g.x, y: r.va }, { name: 'Swirl (tangential)', x: g.x, y: r.vt }] },
        { type: 'bar', title: 'Where the shaft power goes', ylabel: 'Power [kW]', categories: ['Useful thrust power', 'Axial induced', 'Swirl', 'Profile drag'], series: [{ name: 'Power', y: [r.T * op.V, r.Pax, r.Psw, r.Ppr].map((v) => v / 1e3) }] },
        { type: 'line', title: 'Blade chord distribution used', xlabel: 'Radius r/R [-]', ylabel: 'Chord [m]', series: [{ name: 'Chord', x: g.x, y: g.c }] },
      ],
      outputs: { thrust_total_N: r.T * i.n_props, power_total_W: r.P * i.n_props, momentum_residual_T: Math.abs(r.Tm - r.T) / (Math.abs(r.T) || 1), momentum_residual_Q: Math.abs(r.Qm - r.Q) / (Math.abs(r.Q) || 1), energy_residual: Math.abs(r.P - r.T * op.V - r.Pax - r.Psw - r.Ppr) / (r.P || 1) },
      warnings, models: ['Blade-element momentum theory with swirl, solved on the inflow angle', 'Prandtl tip- and hub-loss factors', 'Generic parametric section polar with stall, Reynolds and Mach corrections', i.blade_type === ROTOR_T ? 'Constant-chord, linearly twisted blade' : 'Generic propeller planform c ∝ sqrt(x(1.1 − x)) with constant-pitch helical twist'],
      assumptions: ['Axial, uniform inflow; isolated propeller; rigid blades', 'Blade chord and twist are generic distributions scaled from diameter, blade count and the 75% chord, not the actual blade', 'Section data are a generic polar, not tables for a specific aerofoil', 'Wake contraction and radial flow neglected; annuli are independent', 'Without the propulsion suite the shaft power available is the sea-level rating with a simple density lapse (σ^0.7 turbine, Gagg–Ferrar piston, none electric)'],
    };
  },
  convergence: { param: 'nR', label: 'Radial stations', levels: [8, 16, 32, 64, 128], metric: 'prop_power_W' },
  calibration: { params: [{ key: 'cdmin', min: 0.004, max: 0.03 }, { key: 'cl0', min: 0, max: 0.8 }, { key: 'c75_R', min: 0.03, max: 0.3 }], sweep: 'V', target: 'prop_power_W', note: 'Dynamometer or wind-tunnel shaft power (or thrust) versus speed at known rpm and pitch.' },
  verify() {
    const i = { D: 2, Nb: 3, hub_ratio: 0.15, blade_type: PROP_T, c75_R: 0.14, J_twist: 0.8, twist_deg: 0, rpm: 2400, V: 60, alt_m: 0, dISA: 0, mode: 'match required thrust', T_req: 1500, P_avail: 0, beta75_deg: 25, n_props: 1, tip_loss: true, nR: 30 };
    const r = point.run(i), o = N.kv(r), g = bladeGeom(i), op = opOf(i), f = bemt(g, polarOf(i, { ideal: true }), op, N.rad(22));
    // frictionless annulus efficiency (1 − a′)/(1 + a) summed over the blade
    let Pe = 0; f.dT.forEach((dT, j) => { const a = f.va[j] / op.V, ap = f.vt[j] / (op.Om * g.x[j] * g.R); Pe += (dT * g.dx * g.R * op.V * (1 + a)) / (1 - ap); });
    return [
      N.check('Pitch solve returns the required thrust', o.prop_thrust_N, 1500, 1e-6, 'Trim condition'),
      N.check('Annulus momentum thrust equals blade-element thrust', r.outputs.momentum_residual_T + 1, 1, 1e-8, 'Axial momentum theorem with tip loss'),
      N.check('Annulus angular momentum equals blade-element torque', r.outputs.momentum_residual_Q + 1, 1, 1e-8, 'Angular momentum theorem'),
      N.check('Shaft power = thrust power + axial + swirl + profile losses', r.outputs.energy_residual + 1, 1, 1e-10, 'Energy conservation in the velocity triangle'),
      N.check('Frictionless blade: P = Σ dT·V·(1 + a)/(1 − a′)', f.P, Pe, 1e-8, 'Glauert general momentum theory, η = (1 − a′)/(1 + a)'),
      N.check('η = J·CT/CP', o.eta_prop, (o.J_adv * o.CT_prop) / o.CP_prop, 1e-12, 'Definition of propeller coefficients'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], stat = !(i.V > 0.5);
    if (!stat && o.eta_prop < 0.75) out.push({ severity: 'advise', title: 'Propulsive efficiency is low at this point', detail: `η = ${o.eta_prop.toFixed(2)} against an ideal ${fmt(o.eta_ideal)}; losses: axial ${o.loss_axial_pct.toFixed(0)}%, swirl ${o.loss_swirl_pct.toFixed(0)}%, profile ${o.loss_profile_pct.toFixed(0)}% of shaft power.`, action: o.loss_profile_pct > o.loss_axial_pct ? 'The propeller is lightly loaded for its blade area and tip speed: reduce rpm or blade area in cruise (see the rpm study in the map analysis). Every point of efficiency is a point of cruise fuel or battery energy.' : 'The disk is heavily loaded: a larger diameter lowers induced loss.', basis: 'Power breakdown from blade-element momentum theory' });
    if (o.tip_mach > 0.85) out.push({ severity: 'warn', title: 'High helical tip Mach number', detail: `Tip Mach ${o.tip_mach.toFixed(2)}.`, action: 'Reduce rpm or diameter, or use thin swept tips; this also reduces propeller noise markedly (Suite 11).', basis: 'Section drag divergence near Mach 0.85–0.9' });
    if (stat) out.push({ severity: o.FM_prop < 0.6 ? 'advise' : 'info', title: 'Hover / static efficiency', detail: `Figure of merit ${o.FM_prop.toFixed(2)} at a disk loading of ${o.disk_loading_Pa.toFixed(0)} N/m²; ${(o.prop_thrust_N / o.prop_power_W * 1e3).toFixed(1)} N of thrust per kW.`, action: 'Lower disk loading (larger diameter) and lower tip speed raise thrust per kW and cut hover energy; check blade loading and stall margin as diameter grows.', basis: 'FM = T^1.5 / (P·sqrt(2ρA))' });
    if (Number.isFinite(o.power_margin_pct) && o.power_margin_pct < 0) out.push({ severity: 'critical', title: 'Required thrust needs more power than is available', detail: `Shortfall ${(-o.power_margin_pct).toFixed(0)}%.`, action: 'Reduce speed or mass, or raise installed power; re-check the drag estimate from Suite 1.', basis: 'Shaft power balance' });
    return out;
  },
};

// ---- (b) performance map --------------------------------------------------------------------
const map = {
  id: 'map', title: 'Performance map, constant-speed schedule and rpm study', fidelity: 'numerical',
  summary: 'Efficiency, thrust and power coefficients against advance ratio for several blade pitch settings, the thrust a constant-speed propeller delivers from static to maximum speed at rated power, and the best rpm for the cruise thrust.',
  equations: ['Blade element momentum equations', 'Propulsive efficiency equations', 'Propeller thrust and torque equations'],
  applicable: hasProp,
  inputs: [...GEOM, ...POINT, T_REQ,
    { key: 'P_rated', label: 'Shaft power available per propeller at this altitude', unit: 'W', default: 1.2e6, min: 0, group: 'Operating point' },
    { key: 'beta_min_deg', label: 'Lowest pitch setting', unit: 'deg', default: 20, min: 0, max: 80, group: 'Map' },
    { key: 'beta_max_deg', label: 'Highest pitch setting', unit: 'deg', default: 50, min: 5, max: 85, group: 'Map' },
    { key: 'n_pitch', label: 'Pitch settings', unit: '', default: 4, min: 2, max: 5, step: 1, discrete: true, group: 'Map' },
    { key: 'nJ', label: 'Advance-ratio points', unit: '', default: 22, min: 8, max: 80, step: 1, discrete: true, group: 'Numerics' },
    { key: 'tip_loss', label: 'Prandtl tip and hub loss', type: 'bool', default: true, group: 'Model' },
    ...SECTION, NR],
  defaults: (c, up, d) => {
    const p = pointDefaults(c, up, d), J = p.J_twist || 0, bc = N.deg(Math.atan(J / (0.75 * Math.PI))) + 4, lo = Math.max(8, Math.round(bc - 16));
    return { ...p, V: d.isRotary ? 0 : p.V, P_rated: p.P_avail * (c.meta.type === 'helicopter' ? 0.85 : 1) || undefined, beta_min_deg: lo, beta_max_deg: Math.max(lo + 12, Math.round(bc + 6)) }; // helicopter: 15% of shaft power to tail rotor, transmission and accessories
  },
  run(i, ctx) {
    const g = bladeGeom(i), q = polarOf(i), op0 = opOf(i), nJ = Math.round(i.nJ), betas = N.linspace(i.beta_min_deg, i.beta_max_deg, Math.round(i.n_pitch)), warnings = [];
    const Jmax = Math.max(0.3, 1.15 * Math.PI * 0.75 * Math.tan(N.rad(i.beta_max_deg))), Js = N.linspace(0.02, Jmax, nJ), se = [], sT = [], sP = []; let best = { eta: 0, J: NaN, beta: NaN };
    betas.forEach((bd, kb) => {
      ctx?.progress?.(kb / (betas.length + 2), 'Map');
      const e = [], ct = [], cp = [];
      for (const J of Js) {
        const op = opOf(i, J * op0.n * i.D), r = bemt(g, q, op, N.rad(bd)), k = coeffs(r, op), ok = r.T > 0 && r.P > 0;
        e.push(ok ? r.eta : NaN); ct.push(ok ? k.CT : NaN); cp.push(ok ? k.CP : NaN);
        if (ok && r.eta > best.eta) best = { eta: r.eta, J, beta: bd };
      }
      se.push({ name: `β75 = ${bd.toFixed(0)}°`, x: Js, y: e }); sT.push({ name: `β75 = ${bd.toFixed(0)}°`, x: Js, y: ct }); sP.push({ name: `β75 = ${bd.toFixed(0)}°`, x: Js, y: cp });
    });
    // constant-speed propeller at rated power
    const Vtop = Math.max(10, 1.25 * i.V, 0.35 * op0.n * i.D), Vs = N.linspace(0, Vtop, 14), cs = { T: [], eta: [], beta: [] }, hasP = i.P_rated > 0;
    for (const V of Vs) { const op = opOf(i, V), b = hasP ? solvePitch(g, q, op, 'P', i.P_rated) : NaN, r = Number.isFinite(b) ? bemt(g, q, op, b) : null; cs.T.push(r ? r.T : NaN); cs.eta.push(r ? r.eta : NaN); cs.beta.push(r ? N.deg(b) : NaN); }
    if (hasP && cs.T.some(Number.isNaN)) warnings.push('At some speeds no blade pitch absorbs the rated power at this rpm (the blades stall or the propeller is too small); those points are omitted.');
    if (!hasP) warnings.push('No rated shaft power is given: the constant-speed schedule is skipped.');
    // rpm study at the required thrust and flight speed
    const fr = N.linspace(0.5, 1.1, 9), rp = fr.map((f) => { const op = opOf(i, i.V, f * i.rpm), b = solvePitch(g, q, op, 'T', i.T_req), r = Number.isFinite(b) ? bemt(g, q, op, b) : null; return r ? (i.V > 0.5 ? r.eta : idealP(r.T, op) / r.P) : NaN; });
    const fin = rp.map((v) => (Number.isFinite(v) ? v : -1)), kb = N.argmax(fin), stat = !(i.V > 0.5), mid = bemt(g, q, { ...op0, V: 0, J: 0 }, N.rad(betas[0])), km = coeffs(mid, op0);
    // thrust at rated power at the flight speed itself (not interpolated across points where the governor has no solution)
    const Tst = cs.T[0], bV = hasP ? solvePitch(g, q, op0, 'P', i.P_rated) : NaN, Tcr = Number.isFinite(bV) ? bemt(g, q, op0, bV).T : NaN;
    return {
      kpis: keep([
        { key: 'eta_map_max', label: 'Peak efficiency on the map', value: best.eta, unit: '-' },
        { key: 'J_eta_max', label: 'Advance ratio at peak efficiency', value: best.J, unit: '-' },
        { key: 'beta_eta_max_deg', label: 'Pitch setting at peak efficiency', value: best.beta, unit: 'deg' },
        { key: 'static_thrust_rated_N', label: 'Static thrust at rated power (constant speed)', value: Tst, unit: 'N', opt: true },
        { key: 'thrust_rated_at_V_N', label: 'Thrust at rated power at the flight speed', value: Tcr, unit: 'N', status: Tcr >= i.T_req ? 'ok' : 'bad', note: `Required: ${i.T_req.toFixed(0)} N`, opt: true },
        { key: 'beta_static_deg', label: 'Constant-speed pitch, static', value: cs.beta[0], unit: 'deg', opt: true },
        { key: 'beta_top_deg', label: 'Constant-speed pitch at top of range', value: cs.beta[cs.beta.length - 1], unit: 'deg', opt: true },
        { key: 'rpm_best', label: stat ? 'Best rpm for hover figure of merit' : 'Best rpm for efficiency at the required thrust', value: fin[kb] > 0 ? fr[kb] * i.rpm : NaN, unit: 'rpm', opt: true },
        { key: 'eta_rpm_best', label: stat ? 'Figure of merit at the best rpm' : 'Efficiency at the best rpm', value: fin[kb] > 0 ? fin[kb] : NaN, unit: '-', opt: true },
        { key: 'eta_rpm_nominal', label: stat ? 'Figure of merit at nominal rpm' : 'Efficiency at nominal rpm', value: N.interp1(fr, fin, 1) > 0 ? N.interp1(fr, fin, 1) : NaN, unit: '-', opt: true },
        { key: 'CT_static', label: `Static thrust coefficient at β75 = ${betas[0].toFixed(0)}°`, value: km.CT, unit: '-' },
        { key: 'CP_static', label: `Static power coefficient at β75 = ${betas[0].toFixed(0)}°`, value: km.CP, unit: '-' },
      ]),
      plots: [
        { type: 'line', title: 'Efficiency map', xlabel: 'Advance ratio J = V/(nD) [-]', ylabel: 'Propulsive efficiency [-]', series: se },
        { type: 'line', title: 'Thrust coefficient map', xlabel: 'Advance ratio J [-]', ylabel: 'CT [-]', series: sT },
        { type: 'line', title: 'Power coefficient map', xlabel: 'Advance ratio J [-]', ylabel: 'CP [-]', series: sP },
        { type: 'line', title: 'Constant-speed propeller at rated power: thrust', xlabel: 'Flight speed [m/s]', ylabel: 'Thrust [N]', series: [{ name: 'Thrust available', x: Vs, y: cs.T }, { name: 'Required at flight speed', x: [i.V], y: [i.T_req], style: 'points' }] },
        { type: 'line', title: 'Constant-speed propeller at rated power: pitch schedule', xlabel: 'Flight speed [m/s]', ylabel: 'Blade pitch at 75% radius [deg]', series: [{ name: 'β75', x: Vs, y: cs.beta }] },
        { type: 'line', title: stat ? 'Figure of merit versus rpm at the required thrust' : 'Efficiency versus rpm at the required thrust', xlabel: 'Rotational speed [rpm]', ylabel: stat ? 'Figure of merit [-]' : 'Propulsive efficiency [-]', series: [{ name: 'Pitch re-trimmed at each rpm', x: fr.map((f) => f * i.rpm), y: rp, style: 'line+points' }] },
      ],
      warnings, models: ['Blade-element momentum theory at fixed rpm (Mach and Reynolds effects vary along each curve)', 'Constant-speed governor: pitch solved to absorb rated power', 'Generic section polar and blade geometry'],
      assumptions: ['Map computed at the entered rpm and altitude; it is not a universal J-only map because compressibility is included', 'Rated power is taken as constant with speed; without the propulsion suite it is the sea-level rating with a simple density lapse (σ^0.7 turbine, Gagg–Ferrar piston, none electric)', 'Negative-thrust (windmilling) points are omitted'],
    };
  },
  convergence: { param: 'nR', label: 'Radial stations', levels: [8, 16, 32, 64], metric: 'eta_map_max' },
  verify() {
    // With a Reynolds- and Mach-independent polar the coefficients depend on J and pitch only.
    const i = { D: 1.8, Nb: 2, hub_ratio: 0.15, blade_type: PROP_T, c75_R: 0.13, J_twist: 0.7, twist_deg: 0, alt_m: 0, dISA: 0, nR: 20, tip_loss: true }, g = bladeGeom(i), q = polarOf(i, { ideal: true });
    const a = opOf({ ...i, V: 30, rpm: 2000 }), b = opOf({ ...i, V: 45, rpm: 3000 }), ra = bemt(g, q, a, N.rad(22)), rb = bemt(g, q, b, N.rad(22));
    const big = bemt(bladeGeom({ ...i, D: 3.6 }), q, opOf({ ...i, D: 3.6, V: 30, rpm: 1000 }), N.rad(22));
    // constant-speed thrust at the flight speed against the operating-point analysis absorbing the same power
    const j = { ...Object.fromEntries(map.inputs.map((f) => [f.key, f.default])), D: 2, Nb: 3, c75_R: 0.14, J_twist: 0.8, rpm: 2400, V: 60, alt_m: 0, T_req: 1500, P_rated: 1.5e5, beta_min_deg: 15, beta_max_deg: 35, n_pitch: 2, nJ: 8, nR: 12 };
    const mj = N.kv(map.run(j)), pj = N.kv(point.run({ ...Object.fromEntries(point.inputs.map((f) => [f.key, f.default])), ...j, mode: 'absorb available power', P_avail: 1.5e5 }));
    // default shaft power when the propulsion suite has not run: sea-level rating with the density lapse of the engine type
    const cs = { meta: { type: 'aeroplane' }, atm: { alt_m: 6000, dISA_K: 0 }, flight: { V_ms: 120 }, mission: { cruise_alt_m: 6000, cruise_V_ms: 120 }, wing: { S_m2: 60 }, aero: { CD0: 0.025 }, prop: { type: 'turboprop', n_eng: 2, P0_W: 1.5e6, prop_dia_m: 3.9, n_blades: 6, rpm: 1000 }, rotor: { R_m: 0, chord_m: 0, rpm: 0, n_blades: 0 } };
    const dd = { W: 2e5, k_induced: 0.04, P_total: 3e6, isRotary: false }, sg = isa(6000).sigma, pd = pointDefaults(cs, {}, dd), pp = pointDefaults({ ...cs, prop: { ...cs.prop, type: 'piston' } }, {}, dd), pe = pointDefaults({ ...cs, prop: { ...cs.prop, type: 'electric' } }, {}, dd), pu = pointDefaults(cs, { propulsion: { P_shaft_W: 2.2e6 } }, dd);
    return [
      N.check('Constant-speed thrust at the flight speed equals the operating-point solution at the same power', mj.thrust_rated_at_V_N, pj.prop_thrust_N, 1e-6, 'Same blade-element momentum solution'),
      N.check('Default shaft power at altitude, turbine: P0·σ^0.7', pd.P_avail, 1.5e6 * sg ** 0.7, 1e-12, 'Density lapse of the shaft rating'),
      N.check('Default shaft power at altitude, piston: Gagg–Ferrar', pp.P_avail, (1.5e6 * (sg - 0.1325)) / 0.8675, 1e-12, 'Gagg–Ferrar altitude lapse'),
      N.check('Default shaft power at altitude, electric: no lapse', pe.P_avail, 1.5e6, 1e-12, 'Motor rating independent of air density'),
      N.check('Shaft power from the propulsion suite is used as published', pu.P_avail, 1.1e6, 1e-12, 'All engines at the flight point ÷ number of propellers'),
      N.check('Similarity: CT unchanged at the same J and pitch (rpm ×1.5)', coeffs(rb, b).CT, coeffs(ra, a).CT, 1e-9, 'Dimensional analysis of the propeller'),
      N.check('Similarity: CP unchanged at the same J and pitch (rpm ×1.5)', coeffs(rb, b).CP, coeffs(ra, a).CP, 1e-9, 'Dimensional analysis of the propeller'),
      N.check('Similarity: thrust scales with n²D⁴ (diameter ×2, rpm ÷2)', big.T, ra.T * 4, 1e-9, 'Dimensional analysis of the propeller'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (Number.isFinite(o.rpm_best) && Number.isFinite(o.eta_rpm_nominal) && o.eta_rpm_best > o.eta_rpm_nominal + 0.015) out.push({ severity: 'advise', title: 'A different rpm is more efficient at this thrust', detail: `${o.rpm_best.toFixed(0)} rpm gives ${o.eta_rpm_best.toFixed(3)} against ${o.eta_rpm_nominal.toFixed(3)} at ${i.rpm.toFixed(0)} rpm.`, action: 'Schedule a lower cruise rpm (constant-speed governor or motor speed command). The gain is a direct saving in fuel or battery energy and usually lowers noise as well.', basis: 'Blade-element momentum rpm sweep at constant thrust' });
    if (Number.isFinite(o.thrust_rated_at_V_N) && o.thrust_rated_at_V_N < i.T_req) out.push({ severity: 'critical', title: 'Rated power does not give the required thrust at flight speed', detail: `${o.thrust_rated_at_V_N.toFixed(0)} N available against ${i.T_req.toFixed(0)} N required.`, action: 'Increase rated power or diameter, or reduce drag.', basis: 'Constant-speed thrust at rated power' });
    out.push({ severity: 'info', title: 'Pitch range needed', detail: `A constant-speed propeller needs about ${fmt(o.beta_static_deg, 0)}° (static) to ${fmt(o.beta_top_deg, 0)}° (top speed) at 75% radius; a fixed-pitch propeller must compromise between take-off thrust and cruise efficiency.`, action: 'Use the efficiency map to choose a fixed pitch, or specify the governor range from the schedule.', basis: 'Rated-power pitch schedule' });
    return out;
  },
};

// ---- (c) minimum-induced-loss design --------------------------------------------------------
/** Adkins–Liebeck design iteration. d = {R, Nb, xh, n, V, Om, rho, nu, a, cl, key 'T'|'P', target}; returns geometry and coefficients. */
function milDesign(d, q) {
  const lam = d.V / (d.Om * d.R), dx = (1 - d.xh) / d.n, A = Math.PI * d.R * d.R, Tc = d.key === 'T' ? (2 * d.target) / (d.rho * d.V * d.V * A) : 0, Pc = d.key === 'P' ? (2 * d.target) / (d.rho * d.V ** 3 * A) : 0;
  // incidence that gives the design lift coefficient on the polar at each station's Mach and Reynolds number
  const alOf = (Re, M) => { const r = N.findRoot((a) => polar(q, a, Re, M)[0] - d.cl, -0.3, 0.45, 15, 1e-10); return Number.isFinite(r) ? r : (d.cl - q.cl0) / q.cla; };
  let zeta = 0, out = null;
  for (let it = 0; it < 80; it++) {
    const x = [], c = [], beta = [], G = []; let al75 = 0, I1 = 0, I2 = 0, J1 = 0, J2 = 0;
    const phiT = Math.atan(lam * (1 + zeta / 2));
    for (let j = 0; j < d.n; j++) {
      const xi = d.xh + (j + 0.5) * dx, f = ((d.Nb / 2) * (1 - xi)) / Math.sin(phiT), F = (2 / Math.PI) * Math.acos(Math.min(1, Math.exp(-f))), phi = Math.atan(Math.tan(phiT) / xi), xs = (d.Om * xi * d.R) / d.V;
      const Gj = F * xs * Math.cos(phi) * Math.sin(phi), Wc = (4 * Math.PI * lam * Gj * d.V * d.R * zeta) / (d.cl * d.Nb);
      const Wg = Math.hypot(d.V, d.Om * xi * d.R), Re = Math.max(Wc, 1e-9) / d.nu, al = alOf(Re, Wg / d.a), p = polar(q, al, Re, Wg / d.a), eps = p[1] / Math.max(p[0], 1e-9), tp = Math.tan(phi);
      const a = (zeta / 2) * Math.cos(phi) ** 2 * (1 - eps * tp), W = (d.V * (1 + a)) / Math.sin(phi);
      x.push(xi); c.push(Wc / W); beta.push(al + phi); G.push(Gj); if (xi <= 0.76) al75 = al;
      const i1 = 4 * xi * Gj * (1 - eps * tp), j1 = 4 * xi * Gj * (1 + eps / tp);
      I1 += i1 * dx; I2 += lam * (i1 / (2 * xi)) * (1 + eps / tp) * Math.sin(phi) * Math.cos(phi) * dx; J1 += j1 * dx; J2 += (j1 / 2) * (1 - eps * tp) * Math.cos(phi) ** 2 * dx;
    }
    let zn, tc, pc;
    if (d.key === 'T') { const disc = (I1 / (2 * I2)) ** 2 - Tc / I2; zn = I1 / (2 * I2) - Math.sqrt(Math.max(0, disc)); tc = Tc; pc = J1 * zn + J2 * zn * zn; }
    else { zn = -J1 / (2 * J2) + Math.sqrt((J1 / (2 * J2)) ** 2 + Pc / J2); pc = Pc; tc = I1 * zn - I2 * zn * zn; }
    out = { x, c, beta, G, zeta: zn, Tc: tc, Pc: pc, eta: tc / pc, T: 0.5 * tc * d.rho * d.V * d.V * A, P: 0.5 * pc * d.rho * d.V ** 3 * A, alpha: al75, iterations: it + 1, converged: Math.abs(zn - zeta) < 1e-9 * Math.max(1, Math.abs(zn)) };
    if (out.converged && it > 0) break;
    zeta = zn;
  }
  return out;
}

const design = {
  id: 'design', title: 'Minimum-induced-loss blade design', fidelity: 'numerical',
  summary: 'Designs the chord and twist distribution that gives the required thrust (or absorbs the given power) with the least induced loss, following the Betz condition with Prandtl tip loss and section drag (Adkins–Liebeck), then checks the design with the blade-element momentum solver.',
  equations: ['Betz minimum-energy-loss condition', 'Circulation equations (Goldstein / Prandtl approximation)', 'Prandtl tip-loss correction', 'Blade element momentum equations', 'Propulsive efficiency equations'],
  applicable: hasProp,
  inputs: [...GEOM.filter((f) => ['D', 'Nb', 'hub_ratio', 'rpm'].includes(f.key)), ...POINT,
    { key: 'spec', label: 'Design requirement', type: 'select', options: ['thrust', 'power'], default: 'thrust', group: 'Operating point' },
    T_REQ,
    { key: 'P_design', label: 'Design shaft power (power requirement)', unit: 'W', default: 1.2e6, min: 0, group: 'Operating point' },
    { key: 'cl_design', label: 'Design section lift coefficient', unit: '-', default: 0.6, min: 0.2, max: 1.2, group: 'Blade section (generic polar)', help: 'Near the best lift-to-drag ratio of the section' },
    { key: 'c75_R', label: 'Baseline chord at 75% radius / radius (comparison blade)', unit: '-', default: 0.14, min: 0.01, max: 0.5, group: 'Geometry' },
    ...SECTION, NR],
  defaults: (c, up, d) => {
    const p = pointDefaults(c, up, d), at = isa(p.alt_m, p.dISA), out = { D: p.D, Nb: p.Nb, rpm: p.rpm, hub_ratio: p.hub_ratio, alt_m: p.alt_m, dISA: p.dISA, T_req: p.T_req, c75_R: p.c75_R, cdmin: p.cdmin, V: p.V };
    if (!(p.V > 0) && p.D) out.V = Math.max(2, Math.sqrt(p.T_req / (2 * at.rho * Math.PI * (p.D / 2) ** 2))); // hovering rotors: design at an axial climb equal to the hover induced velocity
    out.P_design = out.V > 0 ? (p.T_req * out.V) / 0.8 : undefined;
    return out;
  },
  run(i) {
    const op = opOf(i), q = polarOf(i), n = Math.round(i.nR), warnings = [], V = Math.max(i.V, 0.5);
    const d = milDesign({ R: i.D / 2, Nb: Math.max(1, Math.round(i.Nb)), xh: N.clamp(i.hub_ratio, 0.02, 0.6), n, V, Om: op.Om, rho: op.rho, nu: op.nu, a: op.a, cl: i.cl_design, key: i.spec === 'power' ? 'P' : 'T', target: i.spec === 'power' ? i.P_design : i.T_req }, q);
    if (!d.converged || !Number.isFinite(d.eta)) warnings.push('The design iteration did not converge: the disk is too heavily loaded for this diameter and speed (no real solution for the displacement velocity).');
    if (i.V < 0.5) warnings.push('The minimum-induced-loss formulation needs forward (axial) speed; 0.5 m/s was used. For hover-optimised rotors design at a small climb speed.');
    // check the designed blade with the analysis solver
    const g = { R: i.D / 2, Nb: Math.max(1, Math.round(i.Nb)), xh: N.clamp(i.hub_ratio, 0.02, 0.6), dx: (1 - N.clamp(i.hub_ratio, 0.02, 0.6)) / n, x: d.x, c: d.c, tw: d.beta, tipLoss: true }, opv = { ...op, V };
    const chk = bemt(g, q, opv, 0), etaI = (d.T * V) / idealP(d.T, opv);
    // baseline generic blade at the same thrust
    const gb = bladeGeom({ ...i, blade_type: PROP_T, J_twist: N.clamp(V / (op.n * i.D), 0.2, 4), twist_deg: 0, tip_loss: true }), bb = solvePitch(gb, q, opv, 'T', d.T), base = Number.isFinite(bb) ? bemt(gb, q, opv, bb) : null;
    const circ = d.G.map((G) => (2 * Math.PI * V * V * d.zeta * G) / (g.Nb * op.Om)), k75 = N.argmin(d.x.map((x) => Math.abs(x - 0.75)));
    const AF = (1e5 / 16) * N.sum(d.x.map((x, j) => (d.c[j] / i.D) * x ** 3 * g.dx)), sol = N.sum(d.c.map((c) => (g.Nb * c * g.dx * g.R) / (Math.PI * g.R * g.R)));
    if (N.amax(d.c) > 0.6 * g.R) warnings.push('The optimum chord is very wide (over 60% of the radius): raise the design lift coefficient, rpm or blade count.');
    return {
      kpis: [
        { key: 'eta_design', label: 'Design efficiency (design theory)', value: d.eta, unit: '-' },
        { key: 'eta_design_bemt', label: 'Efficiency of the designed blade by BEMT analysis', value: chk.eta, unit: '-' },
        { key: 'eta_ideal_disk', label: 'Actuator-disk ideal efficiency', value: etaI, unit: '-' },
        { key: 'eta_baseline', label: 'Generic baseline blade at the same thrust', value: base ? base.eta : NaN, unit: '-' },
        { key: 'eta_gain_pct', label: 'Efficiency gain over the generic blade', value: base ? 100 * (chk.eta - base.eta) : NaN, unit: 'points' },
        { key: 'T_design_N', label: 'Design thrust', value: d.T, unit: 'N' },
        { key: 'P_design_out_W', label: 'Design shaft power', value: d.P, unit: 'W' },
        { key: 'T_check_N', label: 'Thrust of the designed blade by BEMT analysis', value: chk.T, unit: 'N', status: Math.abs(chk.T / d.T - 1) < 0.05 ? 'ok' : 'warn' },
        { key: 'zeta_displacement', label: 'Displacement velocity ratio ζ = v′/V', value: d.zeta, unit: '-' },
        { key: 'beta75_design_deg', label: 'Design pitch at 75% radius', value: N.deg(d.beta[k75]), unit: 'deg' },
        { key: 'chord_max_m', label: 'Maximum chord', value: N.amax(d.c), unit: 'm' },
        { key: 'activity_factor_design', label: 'Activity factor per blade', value: AF, unit: '-' },
        { key: 'solidity_design', label: 'Solidity', value: sol, unit: '-' },
        { key: 'design_iterations', label: 'Design iterations', value: d.iterations, unit: '' },
      ],
      plots: [
        { type: 'line', title: 'Optimum chord distribution', xlabel: 'Radius r/R [-]', ylabel: 'Chord [m]', series: [{ name: 'Minimum induced loss', x: d.x, y: d.c }, { name: 'Generic baseline', x: gb.x, y: gb.c, style: 'dash' }] },
        { type: 'line', title: 'Optimum blade pitch distribution', xlabel: 'Radius r/R [-]', ylabel: 'Blade angle β [deg]', series: [{ name: 'Minimum induced loss', x: d.x, y: d.beta.map(N.deg) }, ...(base ? [{ name: 'Generic baseline (trimmed)', x: gb.x, y: gb.tw.map((t) => N.deg(bb + t)), style: 'dash' }] : [])] },
        { type: 'line', title: 'Bound circulation per blade', xlabel: 'Radius r/R [-]', ylabel: 'Circulation [m²/s]', series: [{ name: 'Γ', x: d.x, y: circ }] },
        { type: 'line', title: 'Thrust loading of the designed blade (BEMT check)', xlabel: 'Radius r/R [-]', ylabel: 'Thrust per unit span [N/m]', series: [{ name: 'Designed blade', x: d.x, y: chk.dT }, ...(base ? [{ name: 'Generic baseline', x: gb.x, y: base.dT, style: 'dash' }] : [])] },
      ],
      tables: [{ title: 'Designed blade geometry', columns: ['r/R', 'Chord [m]', 'Blade angle [deg]'], rows: d.x.map((x, j) => [x, d.c[j], N.deg(d.beta[j])]).filter((_, j) => j % Math.max(1, Math.floor(n / 12)) === 0) }],
      warnings, models: ['Adkins–Liebeck minimum-induced-loss design (Betz condition, Prandtl tip loss, viscous terms)', 'Blade-element momentum check of the designed geometry', 'Generic section polar at the design lift coefficient'],
      assumptions: ['Every section works at the same design lift coefficient', 'Prandtl approximation to the Goldstein circulation function', 'Single design point; off-design behaviour should be checked with the map analysis', 'No structural, noise or manufacturing constraints on chord and twist'],
    };
  },
  convergence: { param: 'nR', label: 'Radial stations', levels: [10, 20, 40, 80], metric: 'eta_design' },
  verify() {
    const i = { D: 2, Nb: 3, hub_ratio: 0.15, rpm: 2400, V: 60, alt_m: 0, dISA: 0, spec: 'thrust', T_req: 1500, P_design: 0, cl_design: 0.7, c75_R: 0.14, nR: 60 }, o = N.kv(design.run(i));
    const p = N.kv(design.run({ ...i, spec: 'power', P_design: o.P_design_out_W }));
    // many frictionless blades at low advance ratio approach the actuator-disk limit
    const op = opOf({ ...i, V: 20, rpm: 3000 }), m = milDesign({ R: 1, Nb: 60, xh: 0.05, n: 200, V: 20, Om: op.Om, rho: op.rho, nu: op.nu, a: op.a, cl: 0.7, key: 'T', target: 300 }, polarOf(i, { ideal: true }));
    return [
      N.check('Blade-element momentum analysis of the designed blade returns the design thrust', o.T_check_N, 1500, 0.02, 'Consistency of design and analysis (Adkins & Liebeck 1994)'),
      N.check('Analysis efficiency of the designed blade equals the design efficiency', o.eta_design_bemt, o.eta_design, 0.01, 'Consistency of design and analysis'),
      N.check('Power-specified design reproduces the thrust-specified design', p.T_design_N, 1500, 1e-6, 'Duality of the thrust and power constraints'),
      N.check('Many frictionless blades at low J: η → 2/(1 + sqrt(1 + Tc))', m.eta, 2 / (1 + Math.sqrt(1 + m.Tc)), 0.01, 'Froude actuator-disk limit (swirl and tip loss vanish)'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (Number.isFinite(o.eta_gain_pct)) out.push({ severity: o.eta_gain_pct > 1.5 ? 'advise' : 'info', title: 'Optimised chord and twist versus the generic blade', detail: `The minimum-induced-loss blade reaches ${o.eta_design_bemt.toFixed(3)} against ${fmt(o.eta_baseline, 3)} for the generic blade (${o.eta_gain_pct.toFixed(1)} points); the actuator-disk limit is ${o.eta_ideal_disk.toFixed(3)}.`, action: 'Adopt the optimum distribution as the starting geometry and check take-off and climb with the map analysis; each efficiency point saves about one percent of cruise energy and CO₂.', basis: 'Betz minimum-energy-loss condition' });
    if (o.eta_ideal_disk - o.eta_design > 0.12) out.push({ severity: 'advise', title: 'Profile and swirl losses dominate', detail: `Design efficiency is ${(100 * (o.eta_ideal_disk - o.eta_design)).toFixed(0)} points below the momentum limit.`, action: 'Lower the rpm or blade count, or raise the design lift coefficient, to cut blade area and profile drag.', basis: 'Comparison with actuator-disk efficiency' });
    return out;
  },
};

// ---- (d) momentum limits, ducted fan, swirl recovery ----------------------------------------
/** Actuator-disk state for thrust T on area A at speed V: induced velocity, far-wake velocity, ideal power. */
function diskState(T, V, rho, A) {
  const vi = -V / 2 + Math.sqrt((V * V) / 4 + T / (2 * rho * A));
  return { vi, Vw: V + 2 * vi, P: T * (V + vi), eta: V > 0 ? V / (V + vi) : 0, Tc: V > 0 ? T / (0.5 * rho * V * V * A) : Infinity, contraction: Math.sqrt((V + vi) / (V + 2 * vi)) };
}
/** Ducted fan with exit area ratio σd = A_exit/A_disk: jet velocity and ideal power for total thrust T. */
function ductState(T, V, rho, A, sd) {
  const Ve = V / 2 + Math.sqrt((V * V) / 4 + T / (rho * A * sd)), md = rho * A * sd * Ve;
  return { Ve, P: 0.5 * md * (Ve * Ve - V * V), eta: V > 0 ? (2 * V) / (V + Ve) : 0, fanShare: V > 0 ? NaN : 1 / (2 * sd) };
}

const disk = {
  id: 'disk', title: 'Momentum-theory limits, ducted fan and swirl recovery', fidelity: 'analytical',
  summary: 'The best any propulsor of this size could do: ideal (Froude) efficiency, slipstream speed and contraction from actuator-disk theory, the benefit of a duct as a function of its diffuser ratio, and an estimate of the swirl energy a contra-rotating pair could recover.',
  equations: ['Actuator-disk momentum equations', 'Froude momentum theory', 'Wake contraction relations', 'Rotor figure-of-merit relations', 'Ducted propulsor momentum model', 'Angular momentum (swirl) relations'],
  applicable: hasProp,
  inputs: [GEOM[0], GEOM[2], GEOM[7], ...POINT, T_REQ,
    { key: 'P_actual', label: 'Actual shaft power per propeller (for comparison)', unit: 'W', default: 1.2e6, min: 0, group: 'Operating point', help: 'From the blade-element analysis when it has run, otherwise a first estimate (case efficiency, or figure of merit 0.7 at rest)' },
    { key: 'sigma_d', label: 'Duct exit area / disk area σd', unit: '-', default: 1.0, min: 0.5, max: 1.6, group: 'Duct', help: '0.5 reproduces an open propeller at rest; 1.0 straight duct; >1 diffusing duct (separation limits it to about 1.2–1.3)' },
    { key: 'swirl_recovery', label: 'Swirl recovery by a contra-rotating rear row', unit: '-', default: 0.7, min: 0, max: 1, group: 'Calibration', help: 'Empirical fraction of the single-rotation swirl loss recovered' }],
  defaults: (c, up, d) => { const p = pointDefaults(c, up, d); return { D: p.D, hub_ratio: p.hub_ratio, rpm: p.rpm, V: p.V, alt_m: p.alt_m, dISA: p.dISA, T_req: p.T_req, P_actual: up.propeller?.prop_power_W ?? (p.V > 0 ? (p.T_req * p.V) / (c.prop.eta_prop || 0.8) : p.D ? p.T_req ** 1.5 / Math.sqrt(2 * isa(p.alt_m, p.dISA).rho * Math.PI * (p.D / 2) ** 2) / 0.7 : undefined) }; },
  run(i) {
    const op = opOf(i), A = op.A * (1 - i.hub_ratio ** 2), T = i.T_req, s = diskState(T, op.V, op.rho, A), du = ductState(T, op.V, op.rho, A, i.sigma_d), stat = op.V < 0.5, warnings = [];
    // swirl of a uniformly loaded (constant-circulation) single rotor: Γ = 2Q/(ρ·Ua·(R² − rh²)); KE flux = ρ·Ua·Γ²·ln(R/rh)/(4π)
    const R = i.D / 2, rh = R * i.hub_ratio, P = i.P_actual > 0 ? i.P_actual : s.P, Q = P / op.Om, Ua = op.V + s.vi, Gam = (2 * Q) / (op.rho * Ua * (R * R - rh * rh)), Psw = (op.rho * Ua * Gam * Gam * Math.log(R / rh)) / (4 * Math.PI), fsw = Psw / P;
    const zs = N.linspace(-2, 4, 61), rw = zs.map((z) => Math.sqrt((op.V + s.vi) / (op.V + s.vi * (1 + z / Math.sqrt(1 + z * z))))), vz = zs.map((z) => op.V + s.vi * (1 + z / Math.sqrt(1 + z * z)));
    const tcs = N.logspace(1e-2, 1e2, 60), sds = N.linspace(0.5, 1.4, 19), Pd = sds.map((sd) => ductState(T, op.V, op.rho, A, sd).P / s.P);
    const FM = stat && i.P_actual > 0 ? s.P / i.P_actual : NaN, etaRel = !stat && i.P_actual > 0 ? (T * op.V) / i.P_actual / s.eta : NaN;
    if (i.sigma_d > 1.25) warnings.push('Diffuser area ratios above about 1.25 separate in practice unless the diffuser is long or boundary-layer control is used; duct drag and weight are not included.');
    if (i.P_actual > 0 && i.P_actual < s.P) warnings.push('The entered actual power is below the momentum-theory minimum for this thrust and diameter: check the inputs.');
    return {
      kpis: keep([
        { key: 'eta_froude', label: 'Ideal (Froude) propulsive efficiency', value: stat ? NaN : s.eta, unit: '-', opt: true },
        { key: 'P_ideal_W', label: 'Ideal power for this thrust', value: s.P, unit: 'W' },
        { key: 'thrust_loading_Tc', label: 'Thrust loading T/(q·A)', value: stat ? NaN : s.Tc, unit: '-', opt: true },
        { key: 'v_induced_disk_ms', label: 'Induced velocity at the disk', value: s.vi, unit: 'm/s' },
        { key: 'V_slipstream_ms', label: 'Fully developed slipstream velocity', value: s.Vw, unit: 'm/s' },
        { key: 'wake_contraction', label: 'Slipstream radius / disk radius', value: s.contraction, unit: '-', note: '0.707 at rest' },
        { key: 'FM_from_power', label: 'Figure of merit from the actual power', value: FM, unit: '-', opt: true },
        { key: 'eta_relative', label: 'Actual efficiency / ideal efficiency', value: etaRel, unit: '-', note: 'Good propellers reach 0.85–0.92 of the momentum limit', opt: true },
        { key: 'duct_power_ratio', label: 'Ducted / open ideal power at the same thrust and disk area', value: du.P / s.P, unit: '-' },
        { key: 'duct_thrust_ratio_static', label: 'Ducted / open static thrust at the same power', value: (2 * i.sigma_d) ** (1 / 3), unit: '-' },
        { key: 'duct_fan_thrust_share', label: 'Share of static thrust carried by the fan', value: 1 / (2 * i.sigma_d), unit: '-', note: 'The remainder acts on the duct lip and diffuser' },
        { key: 'duct_jet_velocity_ms', label: 'Duct exit velocity', value: du.Ve, unit: 'm/s' },
        { key: 'swirl_loss_frac', label: 'Swirl kinetic-energy loss of a single rotor', value: fsw, unit: '-', note: 'Uniform-loading (constant circulation) estimate' },
        { key: 'cr_eta_gain', label: 'Efficiency gain from contra-rotation (estimate)', value: fsw * i.swirl_recovery, unit: '-', note: 'Empirical recovery fraction applied to the swirl loss' },
      ]),
      plots: [
        { type: 'line', title: 'Ideal efficiency versus thrust loading', xlabel: 'Thrust loading Tc = T/(qA) [-]', ylabel: 'Ideal efficiency [-]', xlog: true, series: [{ name: '2/(1 + sqrt(1 + Tc))', x: tcs, y: tcs.map((t) => 2 / (1 + Math.sqrt(1 + t))) }, ...(stat ? [] : [{ name: 'Operating point', x: [s.Tc], y: [s.eta], style: 'points' }])] },
        { type: 'line', title: 'Slipstream contraction', xlabel: 'Axial distance z/R (downstream +) [-]', ylabel: 'Stream-tube radius r/R [-]', series: [{ name: 'Slipstream boundary', x: zs, y: rw }], annotations: [{ x: 0, label: 'Disk' }] },
        { type: 'line', title: 'Axial velocity along the slipstream axis', xlabel: 'Axial distance z/R (downstream +) [-]', ylabel: 'Velocity [m/s]', series: [{ name: 'V + v(z)', x: zs, y: vz }] },
        { type: 'line', title: 'Ducted fan: ideal power relative to the open propeller', xlabel: 'Exit area ratio σd [-]', ylabel: 'Power ratio at equal thrust and disk area [-]', series: [{ name: 'P ducted / P open', x: sds, y: Pd }], annotations: [{ x: i.sigma_d, label: 'Selected' }] },
      ],
      outputs: stat ? {} : { q_ratio: (s.Vw / op.V) ** 2 },
      warnings, models: ['Rankine–Froude actuator disk', 'Axial induced-velocity development v(z) = vi·(1 + z/sqrt(R² + z²))', 'Ducted-fan momentum model with prescribed exit area', 'Constant-circulation swirl estimate with an empirical contra-rotation recovery factor'],
      assumptions: ['Uniform, inviscid, incompressible flow through the disk', 'Duct friction, lip separation and duct drag and weight are not included', 'Contra-rotation benefit is an estimate: real gains depend on row spacing, loading split and acoustics'],
    };
  },
  calibration: { params: [{ key: 'swirl_recovery', min: 0, max: 1 }], sweep: 'V', target: 'cr_eta_gain', note: 'Measured efficiency difference between single- and contra-rotating configurations.' },
  verify() {
    const i = { D: 2, hub_ratio: 0.02, rpm: 2000, V: 50, alt_m: 0, dISA: 0, T_req: 2000, P_actual: 0, sigma_d: 0.5, swirl_recovery: 0.7 }, o = N.kv(disk.run(i)), r = disk.run(i), rho = isa(0).rho, A = Math.PI * (1 - 0.0004);
    const vi = N.brent((v) => 2 * rho * A * (50 + v) * v - 2000, 0, 100, 1e-13), Tc = 2000 / (0.5 * rho * 2500 * A), s = N.kv(disk.run({ ...i, V: 0 }));
    return [
      N.check('Ideal efficiency = 2/(1 + sqrt(1 + Tc))', o.eta_froude, 2 / (1 + Math.sqrt(1 + Tc)), 1e-12, 'Froude momentum theory'),
      N.check('Induced velocity satisfies T = 2ρA(V + vi)·vi', o.v_induced_disk_ms, vi, 1e-9, 'Actuator-disk momentum equation solved independently'),
      N.check('Slipstream dynamic-pressure ratio = 1 + Tc', r.outputs.q_ratio, 1 + Tc, 1e-12, 'Momentum theory identity'),
      N.check('Static slipstream contraction = 1/√2', s.wake_contraction, Math.SQRT1_2, 1e-12, 'Continuity with far-wake velocity 2·vi'),
      N.check('Static ideal power = T^1.5/sqrt(2ρA)', s.P_ideal_W, 2000 ** 1.5 / Math.sqrt(2 * rho * A), 1e-12, 'Rankine–Froude hover power'),
      N.check('Duct with σd = 0.5 at rest equals the open propeller', s.duct_power_ratio, 1, 1e-12, 'Open-rotor far-wake area is half the disk area'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (Number.isFinite(o.eta_froude) && o.eta_froude < 0.9) out.push({ severity: 'advise', title: 'The disk is heavily loaded', detail: `Even an ideal propulsor of this diameter cannot exceed ${(100 * o.eta_froude).toFixed(1)}% efficiency at this thrust (Tc = ${o.thrust_loading_Tc.toFixed(2)}).`, action: 'Increase diameter or the number of propulsors (distributed propulsion) to lower disk loading; this is the main lever on induced loss, energy use and noise.', basis: 'Froude efficiency' });
    if (o.duct_power_ratio < 0.95) out.push({ severity: 'info', title: 'A duct would reduce ideal power', detail: `With σd = ${i.sigma_d.toFixed(2)} the ideal power is ${(100 * o.duct_power_ratio).toFixed(0)}% of the open propeller at equal thrust and diameter.`, action: 'Worth considering for static and low-speed thrust (eVTOL lift fans, UAVs); weigh it against duct drag in cruise, weight and tip-clearance tolerances.', basis: 'Ducted-fan momentum theory' });
    if (o.swirl_loss_frac > 0.04) out.push({ severity: 'info', title: 'Noticeable swirl loss', detail: `About ${(100 * o.swirl_loss_frac).toFixed(1)}% of shaft power leaves as slipstream rotation.`, action: `Contra-rotation or stator vanes could recover roughly ${(100 * o.cr_eta_gain).toFixed(1)} efficiency points (estimate); a higher rpm at lower torque also reduces swirl.`, basis: 'Angular-momentum balance for constant circulation' });
    return out;
  },
};

// ---- (e) installation effects ---------------------------------------------------------------
const install = {
  id: 'install', title: 'Installation effects: slipstream, scrubbing, blockage and compressibility', fidelity: 'reduced-order',
  summary: 'What the airframe sees behind the propeller and what the propeller loses to the airframe: slipstream dynamic pressure on the wing, extra friction drag on wetted surfaces, drag of bodies standing in the wash, and the efficiency lost to compressibility at the blade tips.',
  equations: ['Actuator-disk momentum equations', 'Wake contraction relations', 'Propeller–wing interaction (slipstream dynamic-pressure method)', 'Blade element momentum equations'],
  applicable: hasProp,
  inputs: [...GEOM, ...POINT, T_REQ,
    { key: 'S_wing', label: 'Wing area', unit: 'm²', default: 61, min: 0, group: 'Airframe' },
    { key: 'b_wing', label: 'Wing span', unit: 'm', default: 27, min: 0, group: 'Airframe' },
    { key: 'n_wing_props', label: 'Propellers blowing the wing', unit: '', default: 2, min: 0, max: 32, step: 1, discrete: true, group: 'Airframe' },
    { key: 'CL_wing', label: 'Wing lift coefficient at this point', unit: '-', default: 0.55, min: 0, max: 3, group: 'Airframe' },
    { key: 'S_scrub', label: 'Wetted area washed by one slipstream', unit: 'm²', default: 20, min: 0, group: 'Airframe', help: 'Nacelle plus the wing strip behind the propeller' },
    { key: 'Cf', label: 'Skin-friction coefficient of washed surfaces', unit: '-', default: 0.004, min: 0.001, max: 0.01, group: 'Calibration' },
    { key: 'S_normal', label: 'Bluff area standing across one slipstream', unit: 'm²', default: 0, min: 0, group: 'Airframe', help: 'Fuselage or boom planform under a lifting rotor; 0 for tractor propellers' },
    { key: 'CD_normal', label: 'Drag coefficient of that bluff area', unit: '-', default: 1.0, min: 0.2, max: 2, group: 'Calibration' },
    ...SECTION, NR],
  defaults: (c, up, d) => {
    const p = pointDefaults(c, up, d), at = isa(p.alt_m, p.dISA);
    return { ...p, S_wing: c.wing.S_m2, b_wing: c.wing.b_m, n_wing_props: c.wing.S_m2 > 0 && !d.isRotary ? p.n_props : 0, CL_wing: c.wing.S_m2 > 0 && p.V > 0 ? d.W / (0.5 * at.rho * p.V * p.V * c.wing.S_m2) : 0,
      S_scrub: d.isRotary ? 0 : 2.2 * p.D * (d.mac || p.D) + 1.5 * p.D * p.D, S_normal: d.isRotary ? (c.fuselage.len_m * c.fuselage.dia_m * 0.5) / p.n_props : 0 };
  },
  run(i) {
    const g = bladeGeom(i), q = polarOf(i), op = opOf(i), A = op.A, s = diskState(i.T_req, op.V, op.rho, A), stat = op.V < 0.5, warnings = [];
    const qInf = 0.5 * op.rho * op.V * op.V, qS = 0.5 * op.rho * s.Vw * s.Vw, Dw = i.D * s.contraction; // contracted slipstream diameter
    const scrub = (qS - qInf) * i.Cf * i.S_scrub, bluff = qS * i.CD_normal * i.S_normal, Tinst = i.T_req - scrub - bluff;
    // wing immersed strip: chord × slipstream diameter per propeller, lift scales with local dynamic pressure
    const cbar = i.b_wing > 0 ? i.S_wing / i.b_wing : 0, Sblown = Math.min(i.S_wing, i.n_wing_props * Dw * cbar), dCL = !stat && i.S_wing > 0 ? i.CL_wing * (Sblown / i.S_wing) * (qS / qInf - 1) : 0;
    // compressibility: same thrust with and without Mach effects
    const b = solvePitch(g, q, op, 'T', i.T_req), r = Number.isFinite(b) ? bemt(g, q, op, b) : null, qi = polarOf(i, { mach: false }), bi = solvePitch(g, qi, op, 'T', i.T_req), ri = Number.isFinite(bi) ? bemt(g, qi, op, bi) : null;
    const dComp = r && ri ? (r.P - ri.P) / r.P : NaN, eta = r && !stat ? r.eta : NaN, etaInst = r && !stat ? (Tinst * op.V) / r.P : NaN;
    if (!r) warnings.push('The propeller cannot deliver the required thrust at this point; compressibility and installed-efficiency figures are not available.');
    if (Tinst < 0.8 * i.T_req) warnings.push('Installation losses exceed 20% of thrust: check the washed and bluff areas entered.');
    const Ms = N.linspace(0.4, 1.0, 13), pen = Ms.map((M) => { const rpm = (Math.sqrt(Math.max(1e-6, (M * op.a) ** 2 - op.V ** 2)) / (i.D / 2)) * (60 / (2 * Math.PI)), o2 = opOf(i, op.V, rpm), b2 = solvePitch(g, q, o2, 'T', i.T_req), r2 = Number.isFinite(b2) ? bemt(g, q, o2, b2) : null; return r2 ? (stat ? idealP(r2.T, o2) / r2.P : r2.eta) : NaN; });
    return {
      kpis: keep([
        { key: 'q_ratio_slipstream', label: 'Slipstream / free-stream dynamic pressure', value: stat ? NaN : qS / qInf, unit: '-', opt: true },
        { key: 'q_slipstream_Pa', label: 'Slipstream dynamic pressure', value: qS, unit: 'Pa' },
        { key: 'V_slip_ms', label: 'Slipstream velocity', value: s.Vw, unit: 'm/s' },
        { key: 'slipstream_dia_m', label: 'Contracted slipstream diameter', value: Dw, unit: 'm' },
        { key: 'blown_area_frac', label: 'Wing area inside the slipstreams', value: i.S_wing > 0 ? Sblown / i.S_wing : 0, unit: '-' },
        { key: 'dCL_blown', label: 'Wing lift-coefficient increment from blowing', value: dCL, unit: '-', note: 'Dynamic-pressure scaling on the immersed strip; swirl and upwash effects not included' },
        { key: 'scrub_drag_N', label: 'Scrubbing drag per propeller', value: scrub, unit: 'N' },
        { key: 'bluff_drag_N', label: 'Download / blockage drag per propeller', value: bluff, unit: 'N' },
        { key: 'thrust_installed_N', label: 'Net installed thrust per propeller', value: Tinst, unit: 'N', status: Tinst > 0.9 * i.T_req ? 'ok' : 'warn' },
        { key: 'install_loss_pct', label: 'Installation thrust loss', value: (100 * (i.T_req - Tinst)) / (i.T_req || 1), unit: '%' },
        { key: 'eta_isolated', label: 'Isolated propeller efficiency', value: eta, unit: '-', opt: true },
        { key: 'eta_installed', label: 'Installed propulsive efficiency', value: etaInst, unit: '-', opt: true },
        { key: 'tip_mach_helical', label: 'Helical tip Mach number', value: op.Mtip, unit: '-', status: op.Mtip < 0.85 ? 'ok' : 'warn' },
        { key: 'compress_power_pct', label: 'Extra shaft power due to compressibility', value: 100 * dComp, unit: '%', note: 'Same thrust with and without the Mach terms of the generic polar', opt: true },
      ]),
      plots: [
        { type: 'bar', title: 'Thrust accounting per propeller', ylabel: 'Force [N]', categories: ['Isolated thrust', 'Scrubbing drag', 'Download / blockage', 'Net installed'], series: [{ name: 'Force', y: [i.T_req, -scrub, -bluff, Tinst] }] },
        { type: 'line', title: stat ? 'Figure of merit versus tip Mach number at the required thrust' : 'Efficiency versus helical tip Mach number at the required thrust', xlabel: 'Helical tip Mach number [-]', ylabel: stat ? 'Figure of merit [-]' : 'Propulsive efficiency [-]', series: [{ name: 'rpm varied, pitch re-trimmed', x: Ms, y: pen, style: 'line+points' }], annotations: [{ x: op.Mtip, label: 'Operating point' }] },
      ],
      warnings, models: ['Actuator-disk slipstream (fully developed velocity and contraction)', 'Scrubbing drag from the dynamic-pressure rise on washed wetted area', 'Bluff-body drag of surfaces across the slipstream', 'Blade-element momentum theory with and without compressibility terms'],
      assumptions: ['Slipstream fully developed and uniform at the airframe', 'Wing lift in the slipstream scales with local dynamic pressure only (no swirl, no change of circulation distribution)', 'Upstream blockage of the nacelle on the propeller inflow is not modelled', 'Washed and bluff areas are user estimates; defaults are rough proportions of the case geometry'],
    };
  },
  calibration: { params: [{ key: 'Cf', min: 0.001, max: 0.01 }, { key: 'CD_normal', min: 0.2, max: 2 }], sweep: 'V', target: 'thrust_installed_N', note: 'Installed versus isolated thrust from powered wind-tunnel tests, or hover download measurements.' },
  verify() {
    const i = { D: 2, Nb: 3, hub_ratio: 0.15, blade_type: PROP_T, c75_R: 0.14, J_twist: 0.8, twist_deg: 0, rpm: 2400, V: 60, alt_m: 0, dISA: 0, T_req: 1500, S_wing: 16, b_wing: 10, n_wing_props: 1, CL_wing: 0.4, S_scrub: 4, Cf: 0.004, S_normal: 0.1, CD_normal: 1, nR: 16 };
    const o = N.kv(install.run(i)), rho = isa(0).rho, Tc = 1500 / (0.5 * rho * 3600 * Math.PI);
    return [
      N.check('Slipstream dynamic-pressure ratio = 1 + T/(qA)', o.q_ratio_slipstream, 1 + Tc, 1e-12, 'Momentum theory identity'),
      N.check('Scrubbing drag = Tc·q·Cf·S', o.scrub_drag_N, Tc * 0.5 * rho * 3600 * 0.004 * 4, 1e-12, 'Dynamic-pressure rise on the washed area'),
      N.check('Mass flow conserved through the contraction', o.slipstream_dia_m ** 2 * o.V_slip_ms, 4 * (60 + (o.V_slip_ms - 60) / 2), 1e-12, 'Continuity between disk and far wake'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.install_loss_pct > 5) out.push({ severity: 'advise', title: 'Installation losses are significant', detail: `${o.install_loss_pct.toFixed(1)}% of thrust is lost (${o.scrub_drag_N.toFixed(0)} N scrubbing, ${o.bluff_drag_N.toFixed(0)} N download/blockage per propeller).`, action: i.S_normal > 0 ? 'Slim or fair the structure under the rotor, or move rotors outboard of the fuselage; hover download is paid for directly in power and battery energy.' : 'Reduce the wetted area inside the slipstream (nacelle size, propeller position).', basis: 'Slipstream dynamic pressure on immersed surfaces' });
    if (o.dCL_blown > 0.05) out.push({ severity: 'info', title: 'Useful blown-wing lift', detail: `Slipstream adds about ΔCL = ${o.dCL_blown.toFixed(2)} at this thrust.`, action: 'Account for it in take-off and landing performance (Suite 5) and for the power-off case, where it disappears.', basis: 'Dynamic-pressure scaling of the immersed wing strip' });
    if (Number.isFinite(o.compress_power_pct) && o.compress_power_pct > 3) out.push({ severity: 'warn', title: 'Compressibility is costing power', detail: `${o.compress_power_pct.toFixed(1)}% more shaft power at tip Mach ${o.tip_mach_helical.toFixed(2)}.`, action: 'Lower the rpm in cruise or use thinner, swept tips; see the efficiency-versus-tip-Mach curve.', basis: 'Generic drag-rise model' });
    return out;
  },
};

// ---- (f) electric motor matching ------------------------------------------------------------
/** Equilibrium speed of a DC-equivalent brushless motor (Kv [rad/s/V], R, I0) at terminal voltage Vm against a load torque Qp(ω). */
function matchMotor(m, Vm, Qp) {
  const Qm = (w) => ((Vm - w / m.Kv) / m.R - m.I0) / m.Kv, w0 = Math.max(1e-6, (Vm - m.I0 * m.R) * m.Kv); // no-load speed
  const w = N.findRoot((x) => Qm(x) - Qp(x), 0.02 * w0, w0, 12, 1e-12 * w0 + 1e-13);
  if (!Number.isFinite(w)) return null;
  const I = (Vm - w / m.Kv) / m.R, Q = Qm(w), Pe = Vm * I, Ps = Q * w;
  return { w, I, Q, Pe, Ps, eta: Pe > 0 ? Ps / Pe : 0 };
}

const motor = {
  id: 'motor', title: 'Electric motor–propeller matching', fidelity: 'reduced-order',
  summary: 'Finds the speed at which the motor torque line meets the propeller torque curve for each throttle setting, giving thrust, current, motor and system efficiency, and the thrust margin at full throttle.',
  equations: ['Electric motor–propeller performance coupling', 'Motor torque–speed line (back-EMF, resistance, no-load current)', 'Blade element momentum equations', 'Shaft torque balance'],
  applicable: (c) => (hasProp(c) !== true ? hasProp(c) : c.prop.type === 'electric' ? true : 'Motor matching applies to electrically driven propellers; this vehicle uses a combustion engine (see Suite 7).'),
  inputs: [...GEOM, ...POINT, T_REQ,
    { key: 'V_bus', label: 'DC bus voltage', unit: 'V', default: 44.4, min: 3, group: 'Motor and drive' },
    { key: 'Kv', label: 'Motor speed constant Kv', unit: 'rpm/V', default: 150, min: 0.5, group: 'Motor and drive', help: 'No-load rpm per volt; torque constant Kt = 60/(2π·Kv) N·m/A' },
    { key: 'R_m', label: 'Winding resistance (line equivalent)', unit: 'Ω', default: 0.08, min: 1e-5, group: 'Motor and drive' },
    { key: 'I0', label: 'No-load current', unit: 'A', default: 0.5, min: 0, group: 'Motor and drive' },
    { key: 'I_max', label: 'Continuous current limit', unit: 'A', default: 25, min: 0.1, group: 'Motor and drive' },
    { key: 'eta_esc', label: 'Inverter / speed-controller efficiency', unit: '-', default: 0.97, min: 0.8, max: 1, group: 'Motor and drive' },
    { key: 'tip_loss', label: 'Prandtl tip and hub loss', type: 'bool', default: true, group: 'Model' },
    ...SECTION, NR],
  defaults: (c, up, d) => {
    // generic motor sized from the case rating when no motor data sheet is available
    const p = pointDefaults(c, up, d), Vb = c.systems.bus_V || 48, P0 = c.prop.P0_W || 1000, Ir = P0 / (0.92 * Vb);
    return { ...p, V_bus: Vb, Kv: p.rpm ? (1.3 * p.rpm) / Vb : undefined, R_m: (0.035 * Vb) / Ir, I0: 0.025 * Ir, I_max: 1.15 * Ir };
  },
  run(i) {
    const g = bladeGeom(i), q = polarOf(i), warnings = [], m = { Kv: (i.Kv * 2 * Math.PI) / 60, R: i.R_m, I0: i.I0 }, stat = i.V < 0.5;
    // fixed pitch chosen so that the propeller gives the required thrust at the nominal rpm
    const op0 = opOf(i); let b = solvePitch(g, q, op0, 'T', i.T_req);
    if (!Number.isFinite(b)) { b = Math.atan2(op0.V, 0.75 * op0.Om * g.R) + N.rad(12); warnings.push('The propeller cannot give the required thrust at nominal rpm; a pitch 12° above the helix angle is used.'); }
    // propeller torque and thrust from the blade-element solver at this pitch (memoised on speed)
    const wMax = i.V_bus * m.Kv, memo = new Map(), pr = (w) => { let r = memo.get(w); if (!r) { r = bemt(g, q, opOf(i, i.V, (w * 60) / (2 * Math.PI)), b); memo.set(w, r); } return r; };
    const Qp = (w) => Math.max(0, pr(w).Q), Tp = (w) => pr(w).T;
    const at = (u) => { const r = matchMotor(m, u * i.V_bus, Qp); return r ? { ...r, u, T: Tp(r.w), Pbus: r.Pe / i.eta_esc } : null; };
    const us = N.linspace(0.2, 1, 9), sw = us.map(at), full = sw[sw.length - 1];
    let uo = N.findRoot((u) => (at(u)?.T ?? -1) - i.T_req, 0.1, 1, 9, 1e-7), sat = !Number.isFinite(uo);
    if (sat) { uo = 1; warnings.push('Full throttle does not give the required thrust: the motor–propeller combination is saturated (speed constant too low or pitch too coarse for this voltage).'); }
    const o = at(uo) || full, rpm = (o.w * 60) / (2 * Math.PI), etaP = stat ? NaN : (o.T * i.V) / o.Ps, etaSys = stat ? NaN : (o.T * i.V) / o.Pbus;
    const op = opOf(i, i.V, rpm), FM = idealP(o.T, op) / o.Ps;
    if (o.I > i.I_max) warnings.push(`Operating current ${o.I.toFixed(1)} A exceeds the continuous limit ${i.I_max.toFixed(1)} A.`);
    if (full && full.I > 1.5 * i.I_max) warnings.push('Full-throttle current is more than 1.5 times the continuous limit: the controller must limit current, so the full-throttle point is not sustainable.');
    const ok = sw.filter(Boolean), kq = N.linspace(0.1 * wMax, wMax, 16);
    return {
      kpis: keep([
        { key: 'rpm_op', label: 'Operating speed', value: rpm, unit: 'rpm' },
        { key: 'throttle_op', label: 'Throttle (voltage ratio) for the required thrust', value: uo, unit: '-', status: sat ? 'bad' : uo < 0.8 ? 'ok' : 'warn', note: 'Leave margin for control and battery sag' },
        { key: 'motor_current_A', label: 'Motor current', value: o.I, unit: 'A', status: o.I <= i.I_max ? 'ok' : 'bad' },
        { key: 'P_elec_W', label: 'DC bus power per motor', value: o.Pbus, unit: 'W' },
        { key: 'P_shaft_motor_W', label: 'Shaft power', value: o.Ps, unit: 'W' },
        { key: 'motor_eff_match', label: 'Motor efficiency at the operating point', value: o.eta, unit: '-', status: o.eta > 0.85 ? 'ok' : 'warn' },
        { key: 'eta_prop_match', label: 'Propeller efficiency at the matched speed', value: etaP, unit: '-', opt: true },
        { key: 'FM_match', label: 'Propeller figure of merit (momentum ideal / shaft power)', value: FM, unit: '-' },
        { key: 'eta_system', label: 'Bus-to-thrust-power efficiency', value: etaSys, unit: '-', opt: true },
        { key: 'thrust_per_kW', label: 'Thrust per bus kilowatt', value: (o.T / o.Pbus) * 1e3, unit: 'N/kW' },
        { key: 'thrust_op_N', label: 'Thrust at the operating point', value: o.T, unit: 'N' },
        { key: 'thrust_max_N', label: 'Thrust at full throttle', value: full ? full.T : NaN, unit: 'N' },
        { key: 'thrust_margin', label: 'Full-throttle thrust / required thrust', value: full ? full.T / (i.T_req || 1) : NaN, unit: '-', status: full && full.T / (i.T_req || 1) > (stat ? 1.6 : 1.2) ? 'ok' : 'warn', note: stat ? 'Multirotors need about 1.8–2 for control authority' : 'Climb and gust margin' },
        { key: 'current_max_A', label: 'Current at full throttle', value: full ? full.I : NaN, unit: 'A' },
        { key: 'beta75_match_deg', label: 'Fixed blade pitch at 75% radius', value: N.deg(b), unit: 'deg' },
      ]),
      plots: [
        { type: 'line', title: 'Torque balance: motor lines and propeller curve', xlabel: 'Speed [rpm]', ylabel: 'Torque [N m]', series: [{ name: 'Propeller', x: kq.map((w) => (w * 60) / (2 * Math.PI)), y: kq.map(Qp) }, ...[0.5, 0.75, 1].map((u) => ({ name: `Motor, ${(100 * u).toFixed(0)}% throttle`, x: kq.map((w) => (w * 60) / (2 * Math.PI)), y: kq.map((w) => { const v = ((u * i.V_bus - w / m.Kv) / m.R - m.I0) / m.Kv; return v >= 0 ? v : NaN; }), style: 'dash' })), { name: 'Operating point', x: [rpm], y: [o.Q], style: 'points' }] },
        { type: 'line', title: 'Thrust versus throttle', xlabel: 'Throttle [-]', ylabel: 'Thrust [N]', series: [{ name: 'Thrust', x: ok.map((r) => r.u), y: ok.map((r) => r.T) }], annotations: [{ y: i.T_req, label: 'Required' }] },
        { type: 'line', title: 'Current versus throttle', xlabel: 'Throttle [-]', ylabel: 'Motor current [A]', series: [{ name: 'Current', x: ok.map((r) => r.u), y: ok.map((r) => r.I) }], annotations: [{ y: i.I_max, label: 'Continuous limit' }] },
        { type: 'line', title: 'Motor efficiency versus throttle', xlabel: 'Throttle [-]', ylabel: 'Efficiency [-]', series: [{ name: 'Motor', x: ok.map((r) => r.u), y: ok.map((r) => r.eta) }] },
      ],
      warnings, models: ['Three-constant brushless motor model (Kv, winding resistance, no-load current)', 'Blade-element momentum propeller torque and thrust evaluated at each speed at fixed pitch', 'Constant-efficiency inverter'],
      assumptions: ['Steady state; no field weakening, saturation or temperature rise of the winding (Suite 19 has the dq-axis motor model)', 'No-load current constant with speed', 'Default motor constants are generic values sized from the case rating, not a specific motor', 'Battery sag is represented only through the bus voltage entered'],
    };
  },
  convergence: { param: 'nR', label: 'Radial stations', levels: [8, 16, 32, 64], metric: 'P_elec_W' },
  calibration: { params: [{ key: 'R_m', min: 1e-4, max: 5 }, { key: 'I0', min: 0, max: 50 }, { key: 'Kv', min: 1, max: 5000 }], sweep: 'V_bus', target: 'motor_current_A', note: 'Static thrust-stand measurements of current (and rpm) versus applied voltage.' },
  verify() {
    // Quadratic load Q = k·ω²: the intersection is the root of a quadratic.
    const m = { Kv: 30, R: 0.05, I0: 1 }, Vm = 40, k = 2e-5, r = matchMotor(m, Vm, (w) => k * w * w);
    const a = k * m.Kv, bq = 1 / (m.Kv * m.R), c = -(Vm / m.R - m.I0), w = (-bq + Math.sqrt(bq * bq - 4 * a * c)) / (2 * a);
    return [
      N.check('Motor–load equilibrium speed for Q = k·ω²', r.w, w, 1e-9, 'Closed-form root of Kv·k·ω² + ω/(Kv·R) − (V/R − I0) = 0'),
      N.check('Electrical power = shaft power + copper loss + no-load loss', r.Pe, r.Ps + r.I * r.I * m.R + (r.w / m.Kv) * m.I0, 1e-12, 'Energy conservation in the motor model'),
      N.check('Kirchhoff voltage law: V = I·R + ω/Kv', r.I * m.R + r.w / m.Kv, Vm, 1e-12, 'Motor circuit equation'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.thrust_margin < (i.V < 0.5 ? 1.6 : 1.2)) out.push({ severity: 'warn', title: 'Small thrust margin at full throttle', detail: `Maximum thrust is ${fmt(o.thrust_margin)} × the requirement; throttle at the operating point is ${(100 * o.throttle_op).toFixed(0)}%.`, action: 'Choose a higher Kv or bus voltage, or a finer pitch / smaller propeller, so the torque lines cross at higher speed; check the current limit.', basis: 'Motor–propeller torque balance' });
    if (o.motor_eff_match < 0.85) out.push({ severity: 'advise', title: 'Motor works away from its efficient region', detail: `Motor efficiency ${(100 * o.motor_eff_match).toFixed(0)}% at ${o.motor_current_A.toFixed(1)} A.`, action: 'Peak efficiency occurs near I = sqrt(I0·V/R); select a motor whose efficient current matches cruise or hover, or adjust the propeller to move the operating point. Each point of efficiency extends endurance by about one percent.', basis: 'Copper and no-load loss balance' });
    out.push({ severity: 'info', title: 'System efficiency', detail: i.V < 0.5 ? `${o.thrust_per_kW.toFixed(1)} N of thrust per bus kW in hover (figure of merit ${fmt(o.FM_match)}, motor ${(100 * o.motor_eff_match).toFixed(0)}%).` : `${(100 * o.eta_system).toFixed(0)}% of bus power becomes thrust power (propeller ${(100 * o.eta_prop_match).toFixed(0)}%, motor ${(100 * o.motor_eff_match).toFixed(0)}%, inverter ${(100 * i.eta_esc).toFixed(0)}%).`, action: 'Pass the bus power to Suite 19 for battery sizing and endurance.', basis: 'Chain efficiency' });
    return out;
  },
};

export default {
  id: 'propeller', n: 8,
  tagline: 'How much thrust the propeller or rotor makes for the power it absorbs, how close that is to the physical limit, and how to design and match it.',
  analyses: [point, map, design, disk, install, motor],
  consumes: [{ from: 'propulsion', keys: ['P_shaft_W'], why: 'Shaft power available at the flight point' }],
  provides: [
    { key: 'prop_thrust_N', label: 'Thrust per propeller', unit: 'N' }, { key: 'prop_power_W', label: 'Shaft power per propeller', unit: 'W' }, { key: 'eta_prop', label: 'Propulsive efficiency', unit: '-' },
    { key: 'CT_prop', label: 'Thrust coefficient', unit: '-' }, { key: 'CP_prop', label: 'Power coefficient', unit: '-' }, { key: 'J_adv', label: 'Advance ratio', unit: '-' },
    { key: 'FM_prop', label: 'Figure of merit', unit: '-' }, { key: 'tip_mach', label: 'Helical tip Mach number', unit: '-' },
  ],
  handoff: [
    { model: 'Lifting-line and vortex-lattice propeller models with Goldstein circulation', why: 'The Prandtl tip-loss approximation is used instead of the exact helical-wake solution', tool: 'Lifting-line / vortex-lattice propeller design code' },
    { model: 'Free-vortex wake and coupled blade-element–vortex-wake analysis', why: 'Wake roll-up and non-axial inflow need time-marching wake methods', tool: 'Free-wake or vortex-particle solver' },
    { model: 'Actuator-line and full rotating-blade CFD, propeller CFD–structure coupling', why: 'Three-dimensional viscous compressible flow on moving meshes is beyond in-browser cost', tool: 'RANS/URANS with sliding or overset meshes, coupled FE blade model' },
    { model: 'Section-specific aerofoil tables', why: 'A generic parametric polar is used; real blades need measured or computed tables versus Reynolds and Mach number', tool: 'Aerofoil analysis (panel/boundary-layer or CFD) or wind-tunnel data' },
    { model: 'Propeller at incidence, 1P loads and propeller–wing/airframe interaction in detail', why: 'Only axial flow and slipstream dynamic-pressure effects are modelled', tool: 'Panel or CFD installation analysis; powered wind-tunnel test' },
    { model: 'Contra-rotating and ducted propulsor design', why: 'Only momentum-level estimates with an empirical swirl-recovery factor are given', tool: 'Dedicated ducted-fan / open-rotor design codes and rig tests' },
  ],
};
