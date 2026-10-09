// Suite 15 — Crashworthiness and Impact Mechanics.
// Lumped-mass crash and seat-pulse model with crush zone, stroking seat and spinal (DRI) response; crush-tube
// energy absorbers; head strike and HIC; bird strike by hydrodynamic theory; explicit 1-D elastic–plastic wave
// propagation with Johnson–Cook / Cowper–Symonds rate effects; water impact of a wedge; emergency-landing inertia loads.
// Explicit 3-D finite elements (native kernel): fuselage barrel-section drop test, full-aircraft hybrid crash model
// (masses, nonlinear beams, crush springs) and projectile impact on a panel of solid elements with erosion.

import * as N from '../core/numerics.js';
import { G0 } from '../core/atmosphere.js';
import { METALS } from '../data/materials.js';
import { createFE, flowStress } from '../core/solvers/explicitfe.js';

// ---- shared helpers -------------------------------------------------------------------------
const MATS = Object.keys(METALS);
const mat = (name) => METALS[name] || METALS['Al 2024-T3'];
const MAT = { key: 'material', label: 'Material', type: 'select', options: MATS, default: 'Al 2024-T3', group: 'Material', help: 'Typical handbook values, not design allowables' };
const kpi = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const thin = (a, max = 300) => { if (a.length <= max) return a.slice(); const o = []; for (let k = 0; k < max; k++) o.push(a[Math.round((k * (a.length - 1)) / (max - 1))]); return o; };
const DATA_NOTE = 'Material data are typical handbook values, not design allowables';
const DRI_W = 52.9, DRI_Z = 0.224; // Dynamic Response Index spinal model: natural frequency [rad/s] and damping ratio
const occupants = (c) => (c.mission.pax > 0 ? c.mission.pax + (c.mission.pax > 19 ? 2 : 1) : 0);
/**
 * Vertical seat dynamic-test condition of the aircraft class: velocity change dv [m/s], minimum peak floor deceleration G [g]
 * and latest time of the peak t_rise [s]. Transport aeroplanes: 35 ft/s, 14 g, 0.08 s (14 CFR 25.562(b)(1)); rotorcraft:
 * 30 ft/s, 30 g, 0.031 s (14 CFR 27.562(b)(1), 29.562 alike); small aeroplanes: 31 ft/s, 19 g, 0.05 s for the first row
 * (former 14 CFR 23.562). Powered-lift (eVTOL) vehicles have no rule of their own in this table and take the rotorcraft set.
 */
function seatTest(c) {
  if (c.meta.type === 'helicopter' || c.meta.type === 'evtol') return { cls: 'rotorcraft', dv: 9.14, G: 30, t_rise: 0.031 };
  if (c.mass.mtow_kg >= 5700 || c.mission.pax > 19) return { cls: 'transport aeroplane', dv: 10.67, G: 14, t_rise: 0.08 };
  return { cls: 'small aeroplane', dv: 9.45, G: 19, t_rise: 0.05 };
}
/** Default seat stroke [m]: the class-typical value, raised to 1.3 × the stroke a seat held at nEa g needs when the floor stops within dStop from the speed v. */
const seatStroke = (v, dStop, nEa, base) => N.clamp(1.3 * ((v * v) / (2 * G0 * (nEa - 1)) - dStop), base, 0.4);
const UNSOURCED = 'Seat energy-absorber limit load (12–14.5 g), seat stroke, subfloor crush strength and the Cowper–Symonds, Abramowicz–Jones and bird-Hugoniot constants offered as defaults are typical textbook values that have not been checked against a primary source here: treat them as illustrative and replace them with test data';
/** Cowper–Symonds dynamic strength factor 1 + (ε̇/D)^(1/q). */
const cowperSymonds = (rate, D, q) => 1 + (Math.max(rate, 0) / D) ** (1 / q);
/** Johnson–Cook flow stress [Pa]; reference strain rate 1/s, reference temperature 293 K. */
const johnsonCook = (jc, ep, rate, T = 293) => (jc.A + jc.B * Math.max(ep, 0) ** jc.n) * (1 + jc.C * Math.log(Math.max(rate, 1))) * (1 - N.clamp((T - 293) / (jc.Tm - 293), 0, 1) ** jc.m);
/** Johnson–Cook constants: database values where present, otherwise a simple fit through yield and ultimate (labelled as such by callers). */
const jcOf = (m) => m.JC || { A: m.Sy, B: (m.Su * (1 + m.ef * 0.3) - m.Sy) / (0.3 * m.ef) ** 0.5, n: 0.5, C: 0.01, m: 1, Tm: 900, fitted: true };
/** Elastic–perfectly-plastic link with optional bottoming stop. Returns force for deformation d (compression positive) and updates state s = {dp, wp}. */
function eppLink(s, d, k, Fy, dmax, kb, tension) {
  let F = k * (d - s.dp);
  if (F > Fy) { const ndp = d - Fy / k; s.wp += Fy * (ndp - s.dp); s.dp = ndp; F = Fy; }
  else if (F < -Fy && tension) { const ndp = d + Fy / k; s.wp += Fy * (s.dp - ndp); s.dp = ndp; F = -Fy; }
  if (!tension && F < 0) F = 0;
  s.el = (F * F) / (2 * k);
  if (d > dmax) { F += kb * (d - dmax); s.el += 0.5 * kb * (d - dmax) ** 2; }
  return F;
}
/** Head Injury Criterion over a window [s] from a uniformly sampled acceleration history in g. */
export function hic(t, ag, window) {
  const n = t.length, dt = (t[n - 1] - t[0]) / (n - 1), I = [0]; for (let k = 1; k < n; k++) I.push(I[k - 1] + 0.5 * (ag[k] + ag[k - 1]) * dt);
  const w = Math.max(1, Math.round(window / dt)); let best = 0;
  for (let a = 0; a < n - 1; a++) for (let b = a + 1; b <= Math.min(n - 1, a + w); b++) { const T = (b - a) * dt, av = (I[b] - I[a]) / T; if (av > 0) { const h = T * av ** 2.5; if (h > best) best = h; } }
  return best;
}

// ---- 1. lumped-mass crash / seat pulse ------------------------------------------------------
const PULSE_BASE = { mode: 'vertical impact', v0: 10.67, m_air: 400, m_seat: 55, m_torso: 34, F_crush: 96000, d_crush: 0.3, k_crush: 9.6e6, F_ea: 12700, s_ea: 0.3, k_seat: 2.5e6, G_pulse: 14, t_rise: 0.08, shape: 'triangular', t_end: 0.3, occupied: true, nSteps: 4000 };
const pulse = {
  id: 'pulse', title: 'Crash pulse, seat stroke and spinal load (lumped-mass model)', fidelity: 'numerical',
  summary: 'A chain of airframe, seat and upper-body masses hits the ground through a crushable subfloor, or is driven by a seat-test floor pulse. Explicit time integration gives the crash pulse, structural crush, energy absorbed, seat stroke, the Dynamic Response Index and a lumbar-load estimate.',
  equations: ['Conservation of linear momentum', 'Conservation of energy', 'Transient structural dynamics equations', 'Crushing energy absorption equations', 'Occupant dynamics equations', 'Plasticity constitutive equations'],
  inputs: [
    { key: 'mode', label: 'Scenario', type: 'select', options: ['vertical impact', 'seat test pulse'], default: 'vertical impact', group: 'Scenario', help: 'Vertical impact: the whole chain falls onto the ground. Seat test pulse: the floor follows a prescribed deceleration pulse' },
    { key: 'v0', label: 'Vertical impact velocity', unit: 'm/s', default: 10.67, min: 0.1, max: 40, group: 'Scenario', help: 'Seat-test vertical velocity changes: 10.67 m/s (35 ft/s) for transport aeroplanes, 9.14 m/s (30 ft/s) for rotorcraft, 9.45 m/s (31 ft/s) in the former small-aeroplane rule; military rotorcraft design to 12.8 m/s (42 ft/s). Used by the vertical-impact scenario; the seat-test scenario takes its velocity change from the pulse' },
    { key: 'G_pulse', label: 'Seat-test pulse peak', unit: 'g', default: 14, min: 1, max: 100, group: 'Scenario', help: 'Mainly-vertical seat test: 14 g for transport aeroplanes, 30 g for rotorcraft, 19 g (first row) or 15 g for small aeroplanes; the longitudinal tests use 16 g (transport) and 18.4 g (rotorcraft). Keep peak, rise time and velocity change as one set and confirm against your certification basis' },
    { key: 't_rise', label: 'Seat-test pulse rise time', unit: 's', default: 0.08, min: 0.002, max: 1, group: 'Scenario', help: 'Latest time of the peak: 0.08 s transport, 0.031 s rotorcraft, 0.05 s (first row) or 0.06 s small aeroplanes. A symmetric triangle gives a velocity change of peak × g × rise time' },
    { key: 'shape', label: 'Seat-test pulse shape', type: 'select', options: ['triangular', 'step'], default: 'triangular', group: 'Scenario' },
    { key: 'm_air', label: 'Airframe mass per seat', unit: 'kg', default: 400, min: 0.1, group: 'Masses' },
    { key: 'm_seat', label: 'Seat moving mass + lower body', unit: 'kg', default: 55, min: 0.01, group: 'Masses' },
    { key: 'm_torso', label: 'Upper-body mass on the spine', unit: 'kg', default: 34, min: 0.01, group: 'Masses', help: 'About 34 kg for a 77 kg occupant' },
    { key: 'occupied', label: 'Occupied (evaluate injury metrics)', type: 'bool', default: true, group: 'Masses' },
    { key: 'F_crush', label: 'Subfloor crush plateau force per seat', unit: 'N', default: 96000, min: 1, group: 'Crush zone', help: 'Mean crushing force of the structure below one seat' },
    { key: 'd_crush', label: 'Available crush depth', unit: 'm', default: 0.3, min: 0.001, group: 'Crush zone', help: 'Subfloor structure plus landing-gear stroke' },
    { key: 'k_crush', label: 'Initial crush stiffness', unit: 'N/m', default: 9.6e6, min: 1, group: 'Crush zone' },
    { key: 'F_ea', label: 'Seat energy-absorber limit load', unit: 'N', default: 12700, min: 1, group: 'Seat', help: 'Stroking seats are commonly set at about 12–14.5 g on the effective occupant mass; a rigid seat has a very high value' },
    { key: 's_ea', label: 'Available seat stroke', unit: 'm', default: 0.3, min: 0, group: 'Seat' },
    { key: 'k_seat', label: 'Seat and cushion stiffness', unit: 'N/m', default: 2.5e6, min: 1, group: 'Seat' },
    { key: 't_end', label: 'Simulated time', unit: 's', default: 0.3, min: 0.01, max: 5, group: 'Numerics' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 4000, min: 200, max: 400000, step: 1, discrete: true, group: 'Numerics', help: 'Raised automatically when the explicit stability limit requires it' },
  ],
  defaults(c) {
    const occ = occupants(c), has = occ > 0, mOcc = has ? 77 : Math.max(0.05, c.mass.payload_kg || 0.2 * c.mass.mtow_kg), seats = Math.max(1, occ), mAir = Math.max(0.2 * mOcc, c.mass.mtow_kg / seats - mOcc - (has ? 12 : 0));
    const mSeat = has ? 55 : 0.6 * mOcc, mTor = has ? 34 : 0.4 * mOcc, mTot = mAir + mSeat + mTor, rot = c.meta.type === 'helicopter' || c.meta.type === 'evtol', d = N.clamp(0.12 * (c.fuselage.dia_m || 1), 0.02, 0.5) + (c.gear.stroke_m || 0), st = seatTest(c), v0 = has ? st.dv : Math.max(2, 2 * c.gear.sink_ms); // subfloor depth plus landing-gear stroke
    // plateau force set so the crush zone can absorb the impact energy in about 80% of the available depth, but not below 8 g
    const Fc = mTot * Math.max(8 * G0, (v0 * v0) / (2 * 0.8 * d) + G0), Fea = (has ? 12 : 30) * G0 * (mSeat + mTor);
    return { v0, G_pulse: has ? st.G : undefined, t_rise: has ? st.t_rise : undefined, m_air: mAir, m_seat: mSeat, m_torso: mTor, occupied: has, F_crush: Fc, d_crush: d, k_crush: Fc / Math.max(0.005, 0.05 * d), F_ea: Fea, s_ea: !has ? 0.02 : seatStroke(v0, 0.65 * d, 12, rot ? 0.3 : 0.15), k_seat: Math.max(2.5e6 * (mSeat + mTor) / 89, 100) };
  },
  run(i, ctx) {
    const sled = i.mode === 'seat test pulse', kS = i.m_torso * DRI_W ** 2, cS = 2 * DRI_Z * i.m_torso * DRI_W, kb1 = 20 * i.k_crush, kb2 = 20 * i.k_seat, warnings = [];
    const wmax = Math.sqrt(Math.max(kb1 / i.m_air, kb2 * (1 / i.m_air + 1 / i.m_seat), kS * (1 / i.m_seat + 1 / i.m_torso) * (1 + 2 * DRI_Z) ** 2)), dtc = 1.6 / wmax, n = Math.max(Math.round(i.nSteps), Math.ceil(i.t_end / dtc)), dt = i.t_end / n;
    if (n > Math.round(i.nSteps)) warnings.push(`Time steps raised from ${Math.round(i.nSteps)} to ${n} to respect the explicit stability limit (Δt ≤ 0.8·2/ω_max).`);
    const tp = sled ? (i.shape === 'step' ? i.t_end : 2 * i.t_rise) : 0, aFloor = (t) => (!sled ? 0 : i.shape === 'step' ? -i.G_pulse * G0 : t >= tp ? 0 : -i.G_pulse * G0 * (t < i.t_rise ? t / i.t_rise : 2 - t / i.t_rise));
    const dv = sled ? (i.shape === 'step' ? i.G_pulse * G0 * i.t_end : i.G_pulse * G0 * i.t_rise) : i.v0;
    // positions measured downwards; x1 airframe, x2 seat pan, x3 upper body. Seat test starts in static equilibrium under gravity.
    const st2 = sled ? ((i.m_seat + i.m_torso) * G0) / i.k_seat : 0, st3 = sled ? (i.m_torso * G0) / kS : 0, x = [0, st2, st2 + st3], v = [dv, dv, dv], s1 = { dp: 0, wp: 0, el: 0 }, s2 = { dp: 0, wp: 0, el: 0 };
    const T = [], A1 = [], A2 = [], A3 = [], Fc = [], Dc = [], Fl = []; let stopped = false, wDamp = 0, wExt = 0, pk1 = 0, pk2 = 0, crush = 0, stroke = 0, dSp = 0, lum = 0, eMax = 0;
    const m = [i.m_air, i.m_seat, i.m_torso], E0 = 0.5 * N.sum(m) * dv * dv + 0.5 * i.k_seat * st2 * st2 + 0.5 * kS * st3 * st3, every = Math.max(1, Math.floor(n / 1500));
    for (let s = 0; s <= n; s++) {
      const t = s * dt, F1 = sled ? 0 : eppLink(s1, x[0], i.k_crush, i.F_crush, i.d_crush, kb1, false), F2 = eppLink(s2, x[1] - x[0], i.k_seat, i.F_ea, i.s_ea + i.F_ea / i.k_seat, kb2, true), dsp = x[2] - x[1], F3 = kS * dsp + cS * (v[2] - v[1]);
      const a = [sled ? aFloor(t) : G0 + (F2 - F1) / m[0], G0 + (F3 - F2) / m[1], G0 - F3 / m[2]];
      crush = Math.max(crush, x[0]); stroke = Math.max(stroke, s2.dp); dSp = Math.max(dSp, dsp - st3); lum = Math.max(lum, F3); pk1 = Math.max(pk1, sled ? -a[0] / G0 : (F1 - F2) / (m[0] * G0)); pk2 = Math.max(pk2, (F2 - F3) / (m[1] * G0) );
      if (s % every === 0 || s === n) { T.push(t); A1.push(sled ? -a[0] / G0 : (F1 - F2) / (m[0] * G0)); A2.push(-(a[1] - G0) / G0); A3.push(-(a[2] - G0) / G0); Fc.push(F1); Dc.push(x[0]); Fl.push(F3); }
      if (v[0] <= 0) stopped = true;
      if (s === n) break;
      // symplectic (semi-implicit) Euler: velocities first, then positions
      for (let k = 0; k < 3; k++) { const vn = v[k] + dt * a[k]; if (!(sled && k === 0)) wExt += m[k] * G0 * vn * dt; v[k] = vn; x[k] += dt * vn; }
      wDamp += cS * (v[2] - v[1]) ** 2 * dt;
      const ke = 0.5 * ((sled ? 0 : m[0] * v[0] ** 2) + m[1] * v[1] ** 2 + m[2] * v[2] ** 2), el = s1.el + s2.el + 0.5 * kS * (x[2] - x[1]) ** 2;
      if (!sled) eMax = Math.max(eMax, Math.abs(ke + el + s1.wp + s2.wp + wDamp - E0 - wExt));
      if (s % 2000 === 0) ctx?.progress?.(s / n, 'Integrating');
    }
    const absorbed = s1.wp + s2.wp, Ein = sled ? NaN : E0, dri = (DRI_W ** 2 * dSp) / G0, bottomC = !sled && crush > i.d_crush, bottomS = stroke > i.s_ea && i.s_ea > 0;
    if (bottomC) warnings.push('The crush zone bottoms out: the remaining energy goes into the stiff structure above and the floor deceleration spikes.');
    if (stroke > i.s_ea) warnings.push(i.s_ea > 0 ? 'The seat uses all of its stroke and bottoms out.' : 'The seat limit load is exceeded but no stroke is available.');
    if (!i.occupied) warnings.push('No occupants in this case: the "upper body" is the payload on a 52.9 rad/s mount and injury metrics are not meaningful.');
    if (!sled && !stopped) warnings.push('The airframe is still moving downwards at the end of the simulated time; extend it.');
    const lumLim = 6672, driSt = !i.occupied ? undefined : dri <= 18 ? 'ok' : dri <= 22.8 ? 'warn' : 'bad';
    return {
      kpis: [
        kpi('peak_g', sled ? 'Floor pulse peak' : 'Airframe (floor) peak deceleration', pk1, 'g'),
        kpi('crush_m', 'Structural crush', Math.min(crush, sled ? 0 : Infinity), 'm', sled ? undefined : bottomC ? 'bad' : crush > 0.85 * i.d_crush ? 'warn' : 'ok', `Available ${i.d_crush.toFixed(3)} m`),
        kpi('absorbed_J', 'Energy absorbed by crush zone and seat', absorbed, 'J'),
        kpi('absorbed_frac', 'Share of impact energy absorbed plastically', sled ? NaN : absorbed / (E0 + wExt), '-'),
        kpi('seat_stroke_m', 'Seat stroke used', stroke, 'm', bottomS ? (i.occupied ? 'bad' : 'warn') : 'ok', `Available ${i.s_ea.toFixed(3)} m`),
        kpi('pelvis_peak_g', 'Seat-pan (pelvis) peak deceleration', pk2, 'g'),
        kpi('DRI', 'Dynamic Response Index', dri, '-', driSt, 'About 18 corresponds to roughly 5% and 22.8 to 50% probability of spinal injury in ejection-seat experience'),
        kpi('lumbar_N', 'Lumbar load estimate (spine spring + damper)', lum, 'N', !i.occupied ? undefined : lum <= lumLim ? 'ok' : 'bad', 'Seat-test pass criterion is 6672 N (1500 lbf); this lumped estimate is indicative only'),
        kpi('dv_ms', 'Velocity change', dv, 'm/s'),
        kpi('energy_err', 'Energy balance error', sled ? 0 : eMax / E0, '-', sled || eMax / E0 < 0.01 ? 'ok' : 'warn'),
        kpi('n_steps', 'Time steps used', n, '-'),
      ],
      plots: [
        { type: 'line', title: 'Deceleration histories', xlabel: 'Time [ms]', ylabel: 'Deceleration [g]', series: [{ name: 'Airframe floor', x: T.map((t) => t * 1e3), y: A1 }, { name: 'Seat pan', x: T.map((t) => t * 1e3), y: A2 }, { name: 'Upper body', x: T.map((t) => t * 1e3), y: A3 }] },
        ...(sled ? [] : [{ type: 'line', title: 'Crush force versus crush distance', xlabel: 'Crush [mm]', ylabel: 'Force [kN]', series: [{ name: 'Subfloor', x: Dc.map((d) => d * 1e3), y: Fc.map((f) => f / 1e3) }], annotations: [{ x: i.d_crush * 1e3, label: 'Bottoming' }] }]),
        { type: 'line', title: 'Spinal load', xlabel: 'Time [ms]', ylabel: 'Force [kN]', series: [{ name: 'Lumbar load estimate', x: T.map((t) => t * 1e3), y: Fl.map((f) => f / 1e3) }], annotations: i.occupied ? [{ y: lumLim / 1e3, label: '6672 N' }] : [] },
      ],
      outputs: { energy_in_J: Ein },
      warnings,
      models: ['Three-mass chain with elastic–perfectly-plastic crush zone and seat energy absorber', 'Dynamic Response Index spinal model (52.9 rad/s, ζ = 0.224)', 'Explicit symplectic Euler integration with stability-limited step'],
      assumptions: ['Purely vertical, one-dimensional impact on a rigid surface; no pitch, roll or horizontal velocity', 'Constant crush plateau force; bottoming modelled as a 20× stiffer stop', 'The occupant is a single upper-body mass on the DRI spring–damper; no restraint, flail or head strike', 'Injury thresholds are indicative; the lumbar load of a lumped model is not a substitute for a dummy test', UNSOURCED],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [2000, 4000, 8000, 16000, 32000], metric: 'DRI' },
  verify() {
    // (a) step floor pulse on a stiff seat: DRI = G·(1 + exp(−ζπ/√(1−ζ²)))
    const st = N.kv(pulse.run({ ...PULSE_BASE, mode: 'seat test pulse', shape: 'step', G_pulse: 10, F_ea: 1e9, k_seat: 2e9, m_seat: 1000, t_end: 0.25, nSteps: 20000 }));
    // (b) rigid chain on a plateau: m·v²/2 + m·g·δ = F·(δ − F/2k)
    const mT = 489, F = 96000, k = 9.6e7, v0 = 6, im = N.kv(pulse.run({ ...PULSE_BASE, v0, m_air: 488.99, m_seat: 0.005, m_torso: 0.005, k_crush: k, d_crush: 5, F_ea: 1e9, k_seat: 5e4, t_end: 0.25, nSteps: 20000 })), dExact = (0.5 * mT * v0 * v0 + (F * F) / (2 * k)) / (F - mT * G0);
    const sl = pulse.run({ ...PULSE_BASE, mode: 'seat test pulse', d_crush: 0.01 }), slk = N.kv(sl);
    return [
      N.check('Seat-test pulse: velocity change of the triangular pulse = G·g·t_rise', slk.dv_ms, 14 * G0 * 0.08, 1e-12, 'Area of a symmetric triangle; 14 g peaking at 0.08 s gives 10.98 m/s, above the 10.67 m/s (35 ft/s) required of the transport seat test'),
      N.check('Seat-test pulse: the crush zone is not part of the scenario (no bottoming reported)', sl.warnings.filter((w) => w.includes('crush zone')).length, 0, 0, 'The floor motion is prescribed'),
      N.check('DRI for a step input: G·(1 + overshoot of a damped oscillator)', st.DRI, 10 * (1 + Math.exp((-DRI_Z * Math.PI) / Math.sqrt(1 - DRI_Z ** 2))), 2e-3, 'Step response of a single-degree-of-freedom system (Stech & Payne DRI model)'),
      N.check('Crush distance from the energy balance on a plateau', im.crush_m, dExact, 2e-3, 'Work–energy theorem for an (almost) single mass'),
      N.check('Energy is conserved by the integrator', im.energy_err, 0, 5e-3, 'Kinetic + elastic + plastic + damped = initial + gravity work'),
      N.check('Peak deceleration = plateau force / weight for a single mass', im.peak_g, F / (mT * G0), 2e-3, 'Newton’s second law'),
    ];
  },
  calibration: { params: [{ key: 'F_crush', min: 100, max: 1e7 }, { key: 'k_crush', min: 1e3, max: 1e10 }], sweep: 'v0', target: 'peak_g', note: 'Drop-test floor decelerations at several impact velocities calibrate the subfloor crush force and stiffness.' },
  recommend(res, i) {
    const o = res.outputs, out = [], sled = i.mode === 'seat test pulse';
    if (!sled && o.crush_m > i.d_crush) out.push({ severity: 'critical', title: 'Crush zone bottoms out', detail: `Needs more than the ${i.d_crush.toFixed(2)} m available at ${i.v0} m/s.`, action: `Raise the crush force towards ${(o.energy_in_J / (0.8 * i.d_crush) / 1e3).toFixed(0)} kN per seat, deepen the subfloor, or add an energy-absorbing landing gear stage.`, basis: 'Energy absorbed = plateau force × crush distance' });
    if (i.occupied && o.DRI > 18) out.push({ severity: o.DRI > 22.8 ? 'critical' : 'warn', title: 'High spinal injury risk', detail: `DRI ${o.DRI.toFixed(1)}, lumbar estimate ${(o.lumbar_N / 1e3).toFixed(1)} kN.`, action: o.seat_stroke_m >= i.s_ea ? 'Provide (more) seat stroke: about ' + ((i.v0 * i.v0) / (2 * G0 * 13.5)).toFixed(2) + ' m is needed to hold 14.5 g at this velocity.' : 'Lower the seat energy-absorber limit load or the subfloor crush force.', basis: 'DRI ≈ 18 for about 5% spinal injury probability; 6672 N lumbar limit in seat tests' });
    if (i.occupied && o.DRI <= 18 && !sled) out.push({ severity: 'info', title: 'Occupant loads are within tolerance for this pulse', detail: `DRI ${o.DRI.toFixed(1)} with ${(100 * o.seat_stroke_m).toFixed(0)} cm of seat stroke and ${(100 * o.crush_m).toFixed(0)} cm of crush.`, action: 'Unused crush depth or stroke is carried mass on every flight: trimming it saves fuel or energy, but keep margin for heavier occupants and combined-velocity impacts.', basis: 'Crashworthiness system trade' });
    return out;
  },
};

// ---- 2. crush-tube energy absorber ----------------------------------------------------------
/** Static mean crushing force [N] of a thin-walled tube. D: diameter or side length. */
function meanCrush(shape, method, s0, D, t) {
  if (shape === 'square') return method === 'Abramowicz–Jones' ? 13.06 * s0 * t ** (5 / 3) * D ** (1 / 3) : 9.56 * s0 * t ** (5 / 3) * D ** (1 / 3);
  return method === 'Abramowicz–Jones' ? (s0 * t * t / 4) * (20.79 * Math.sqrt(D / t) + 11.9) : (2 * (Math.PI * t) ** 1.5 * Math.sqrt(D / 2) * s0) / 3 ** 0.25;
}
const absorber = {
  id: 'absorber', title: 'Crush-tube energy absorber sizing', fidelity: 'analytical',
  summary: 'Mean progressive-crushing force of circular or square thin-walled metal tubes, with a strain-rate correction, and the number and length of tubes needed to stop a mass within a deceleration limit.',
  equations: ['Crushing energy absorption equations', 'Conservation of energy', 'Strain-rate-dependent material equations'],
  inputs: [
    { key: 'mass', label: 'Mass to be arrested', unit: 'kg', default: 500, min: 0.01, group: 'Requirement' },
    { key: 'v0', label: 'Impact velocity', unit: 'm/s', default: 9.14, min: 0.1, max: 60, group: 'Requirement', help: '9.14 m/s (30 ft/s) is the rotorcraft seat-test velocity change; transport aeroplanes use 10.67 m/s (35 ft/s)' },
    { key: 'G_lim', label: 'Deceleration limit', unit: 'g', default: 20, min: 1.5, max: 200, group: 'Requirement' },
    { key: 'shape', label: 'Tube section', type: 'select', options: ['circular', 'square'], default: 'circular', group: 'Tube' },
    { key: 'D', label: 'Diameter or side length', unit: 'm', default: 0.06, min: 0.003, group: 'Tube' },
    { key: 't', label: 'Wall thickness', unit: 'm', default: 0.0015, min: 1e-4, group: 'Tube' },
    { key: 'L', label: 'Tube length', unit: 'm', default: 0.3, min: 0.01, group: 'Tube' },
    { key: 'n_tubes', label: 'Number of tubes', unit: '', default: 4, min: 1, max: 1000, step: 1, discrete: true, group: 'Tube' },
    { key: 'eff', label: 'Stroke efficiency', unit: '-', default: 0.72, min: 0.3, max: 0.9, group: 'Tube', help: 'Usable crush / length before the folds lock up; about 0.7–0.75 for metal tubes (empirical)' },
    MAT,
    { key: 'method', label: 'Mean-force formula', type: 'select', options: ['Abramowicz–Jones', 'Alexander / Wierzbicki–Abramowicz'], default: 'Abramowicz–Jones', group: 'Method' },
    { key: 'cs_D', label: 'Cowper–Symonds D', unit: '1/s', default: 6500, min: 1, group: 'Method', help: 'Typical textbook values (not verified against a primary source here): about 6500 1/s for aluminium alloys and 40.4 1/s for mild steel; high-strength steels and titanium are far less rate-sensitive' },
    { key: 'cs_q', label: 'Cowper–Symonds q', unit: '-', default: 4, min: 1, max: 20, group: 'Method', help: '4 for aluminium, 5 for mild steel' },
  ],
  defaults(c) {
    const occ = occupants(c), m = c.mass.mtow_kg, v0 = occ > 0 ? seatTest(c).dv : Math.max(2, 2 * c.gear.sink_ms), mm = mat(c.struct.material), n = 4, G = 20, s0 = 0.5 * (mm.Sy + mm.Su);
    // tube scaled so that the set gives about 90% of the limit deceleration (Pm ≈ 6·σ0·t·√(D·t) with D/t = 40); fewer, larger tubes when they would be slender
    const stroke = (v0 * v0) / (2 * G0 * (0.9 * G - 1)), base = { mass: m, v0, G_lim: G, shape: 'circular', eff: 0.72, material: c.struct.material, method: 'Abramowicz–Jones', cs_D: mm.rho < 3500 ? 6500 : 1e5, cs_q: 4 };
    let o = base;
    for (const nT of [n, 2, 1]) {
      let t = Math.max(2e-4, Math.sqrt((0.9 * G * m * G0) / (nT * 6 * s0 * Math.sqrt(40)) / 1.2));
      let g = 0;
      for (let k = 0; k < 8; k++) { o = { ...base, n_tubes: nT, t, D: 40 * t, L: Math.max(3 * 40 * t, (1.2 * stroke) / 0.72) }; g = N.kv(absorber.run(o)).decel_g; if (!(g > 0)) break; t = Math.max(2e-4, t * ((0.9 * G) / g) ** 0.6); }
      if (o.L <= 6 * o.D && g <= G) break;
    }
    return o;
  },
  run(i) {
    const m = mat(i.material), s0 = 0.5 * (m.Sy + m.Su), n = Math.round(i.n_tubes), warnings = [], rate = i.v0 / (2 * i.D), dyn = cowperSymonds(rate, i.cs_D, i.cs_q);
    const PmS = meanCrush(i.shape, i.method, s0, i.D, i.t), Pm = PmS * dyn, Ftot = n * Pm, W = i.mass * G0, G = Ftot / W, stroke = i.eff * i.L;
    const sNeed = Ftot > W ? (0.5 * i.mass * i.v0 ** 2) / (Ftot - W) : Infinity, sLim = i.v0 ** 2 / (2 * G0 * (i.G_lim - 1)), A = i.shape === 'square' ? 4 * i.D * i.t : Math.PI * i.D * i.t, mTube = m.rho * A * i.L, sea = (Pm * stroke) / mTube;
    const nReq = Math.ceil((0.5 * i.mass * i.v0 ** 2 + W * stroke) / (Pm * stroke)), Pmax = s0 * dyn * A, euler = (Math.PI ** 2 * m.E * (i.shape === 'square' ? (2 / 3) * i.D ** 3 * i.t : (Math.PI / 8) * i.D ** 3 * i.t)) / (4 * i.L * i.L);
    if (i.D / i.t > 100 || i.D / i.t < 10) warnings.push('D/t outside about 10–100: the folding-mode formulas are outside their test range (thick tubes split or crush as a solid, very thin ones fold irregularly).');
    if (euler < 2 * Pmax) warnings.push('The tube is slender (Euler load of the fixed–free tube is below twice the peak crush load): it may bend globally instead of folding progressively. Shorten it or increase the section.');
    if (sNeed > stroke) warnings.push('The usable stroke is too short to absorb the impact energy: the absorber bottoms out.');
    if (G > i.G_lim) warnings.push('The crushing force exceeds the deceleration limit.');
    const Gs = N.linspace(Math.max(2, 0.3 * i.G_lim), 2.5 * i.G_lim, 50), x = [0, 0.02 * stroke, 0.04 * stroke, stroke, stroke * 1.02].map((v) => v * 1e3);
    return {
      kpis: [
        kpi('Pm_N', 'Mean crushing force per tube (dynamic)', Pm, 'N'), kpi('Pm_static_N', 'Mean crushing force per tube (static)', PmS, 'N'),
        kpi('dyn_factor', 'Strain-rate factor', dyn, '-', undefined, `Characteristic strain rate ≈ ${rate.toFixed(0)} 1/s (order-of-magnitude estimate v/2D)`),
        kpi('decel_g', 'Deceleration of the arrested mass', G, 'g', G <= i.G_lim ? 'ok' : 'bad', `Limit ${i.G_lim} g`),
        kpi('stroke_needed_m', 'Stroke needed', sNeed, 'm', sNeed <= stroke ? 'ok' : 'bad', `Usable ${stroke.toFixed(3)} m`),
        kpi('stroke_at_limit_m', 'Minimum stroke at the deceleration limit', sLim, 'm'),
        kpi('energy_capacity_J', 'Energy capacity of the absorber set', Ftot * stroke, 'J'), kpi('energy_demand_J', 'Kinetic energy to absorb', 0.5 * i.mass * i.v0 ** 2, 'J'),
        kpi('n_required', 'Tubes required for the energy', nReq, '-', nReq <= n ? 'ok' : 'warn'),
        kpi('sea_Jkg', 'Specific energy absorption', sea, 'J/kg'), kpi('absorber_mass_kg', 'Mass of the tube set', n * mTube, 'kg'),
        kpi('peak_to_mean', 'Initial peak / mean force estimate', Pmax / Pm, '-', undefined, 'Squash load over mean load; a trigger (dent or chamfer) is needed to cut the initial peak'),
      ],
      plots: [
        { type: 'line', title: 'Idealised crush force–stroke curve (set of tubes)', xlabel: 'Stroke [mm]', ylabel: 'Force [kN]', series: [{ name: 'Mean-force idealisation', x, y: [0, (n * Pmax) / 1e3, Ftot / 1e3, Ftot / 1e3, (3 * Ftot) / 1e3] }], annotations: [{ y: (i.G_lim * W) / 1e3, label: 'Deceleration limit' }] },
        { type: 'line', title: 'Stroke required versus deceleration level', xlabel: 'Deceleration [g]', ylabel: 'Stroke [m]', series: [{ name: 'v²/(2g(G−1))', x: Gs, y: Gs.map((g) => i.v0 ** 2 / (2 * G0 * (g - 1))) }, { name: 'This design', x: [G], y: [Number.isFinite(sNeed) ? sNeed : 0], style: 'points' }], annotations: [{ y: stroke, label: 'Usable stroke' }] },
      ],
      warnings,
      models: [i.shape === 'circular' ? (i.method === 'Abramowicz–Jones' ? 'Abramowicz–Jones axisymmetric crushing (1984, empirical-theoretical)' : 'Alexander concertina-mode theory (1960)') : i.method === 'Abramowicz–Jones' ? 'Abramowicz–Jones square-tube crushing (1984)' : 'Wierzbicki–Abramowicz square-tube theory (1983)', 'Cowper–Symonds strain-rate factor on the flow stress'],
      assumptions: ['Progressive axial folding at a constant mean force; flow stress = mean of yield and ultimate', 'Vertical impact: the weight continues to act during the stroke', 'Stroke efficiency and strain-rate constants are empirical inputs', 'Composite or honeycomb absorbers need their own test-derived crush stress', UNSOURCED, DATA_NOTE],
    };
  },
  verify() {
    const b = { mass: 500, v0: 9, G_lim: 20, shape: 'circular', D: 0.06, t: 0.0015, L: 0.3, n_tubes: 4, eff: 0.72, material: 'Al 2024-T3', method: 'Abramowicz–Jones', cs_D: 1e30, cs_q: 4 }, r = N.kv(absorber.run(b));
    return [
      N.check('Alexander mean force ≈ 6·σ0·t·√(D·t)', meanCrush('circular', 'Alexander', 1, 0.06, 0.0015), 6 * 0.0015 * Math.sqrt(0.06 * 0.0015), 5e-3, 'Alexander (1960)'),
      N.check('Abramowicz–Jones square tube: Pm/M0 = 52.22·(b/t)^(1/3)', meanCrush('square', 'Abramowicz–Jones', 4, 0.05, 0.002) / 0.002 ** 2, 52.22 * 25 ** (1 / 3), 5e-4, 'Abramowicz & Jones (1984)'),
      N.check('Work–energy: (F − W)·stroke = ½·m·v²', (r.decel_g - 1) * 500 * G0 * r.stroke_needed_m, 0.5 * 500 * 81, 1e-10, 'Work–energy theorem'),
      N.check('Cowper–Symonds factor is 2 at ε̇ = D', cowperSymonds(6500, 6500, 4), 2, 1e-12, 'Cowper & Symonds (1957)'),
    ];
  },
  calibration: { params: [{ key: 'eff', min: 0.3, max: 0.9 }, { key: 'cs_D', min: 1, max: 1e6 }], sweep: 'v0', target: 'Pm_N', note: 'Measured mean crush forces from static and drop-tower tests calibrate the strain-rate constant; crush length gives the stroke efficiency.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.decel_g > i.G_lim) out.push({ severity: 'warn', title: 'Absorber is too stiff', detail: `${o.decel_g.toFixed(1)} g against a ${i.G_lim} g limit.`, action: 'Use fewer or thinner tubes: mean force scales with t^1.5–1.67.', basis: 'Mean crushing force / weight' });
    if (o.stroke_needed_m > i.eff * i.L) out.push({ severity: 'critical', title: 'Insufficient stroke', detail: `${o.stroke_needed_m.toFixed(2)} m needed, ${(i.eff * i.L).toFixed(2)} m usable.`, action: `Lengthen the tubes to at least ${(o.stroke_needed_m / i.eff).toFixed(2)} m or raise the force towards the limit (${o.n_required} tubes of this size absorb the energy).`, basis: 'Energy balance' });
    out.push({ severity: 'info', title: 'Mass efficiency', detail: `Specific energy absorption ${(o.sea_Jkg / 1e3).toFixed(1)} kJ/kg; absorber mass ${o.absorber_mass_kg.toFixed(2)} kg.`, action: 'Metal tubes typically reach 15–30 kJ/kg and well-designed composite absorbers more; a higher value directly reduces carried mass and trip energy.', basis: 'Specific energy absorption' });
    return out;
  },
};

// ---- 3. head strike and HIC -----------------------------------------------------------------
function headImpact(i, n) {
  // head mass decelerated by a pad: linear elastic up to the foam plateau force, densification at 70% compression
  const Fp = i.sig_pad * i.A_contact, k = (i.E_pad * i.A_contact) / i.h_pad, xd = 0.7 * i.h_pad, kd = 30 * k, wmax = Math.sqrt(kd / i.m_head), tEnd = Math.max(0.03, (6 * i.h_pad) / i.v), nn = Math.max(n, Math.ceil(tEnd / (0.2 / wmax))), dt = tEnd / nn;
  let x = 0, v = i.v, xp = 0, xmax = 0, end = false; const t = [], a = [];
  for (let s = 0; s <= nn; s++) {
    let F = k * (x - xp); if (F > Fp) { xp = x - Fp / k; F = Fp; } if (F < 0) F = 0; if (x > xd) F += kd * (x - xd);
    if (x < xp && v < 0) end = true;
    t.push(s * dt); a.push(end ? 0 : F / (i.m_head * G0)); xmax = Math.max(xmax, x);
    v -= (dt * (end ? 0 : F)) / i.m_head; x += dt * v;
  }
  return { t, a, xmax, Fp, k, bottom: xmax > xd, rebound: -v };
}
const HEAD_BASE = { v: 6, m_head: 4.5, h_pad: 0.05, sig_pad: 3.5e5, E_pad: 8e6, A_contact: 0.008, nSteps: 3000 };
const head = {
  id: 'head', title: 'Head strike on a padded surface and HIC', fidelity: 'numerical',
  summary: 'A head form strikes an energy-absorbing pad (seat back, bulkhead or panel). The deceleration pulse is integrated and the Head Injury Criterion is evaluated over 15 ms and 36 ms windows.',
  equations: ['Conservation of linear momentum', 'Conservation of energy', 'Impact contact equations', 'Occupant dynamics equations'],
  applicable: (c) => (occupants(c) > 0 ? true : 'No occupants are defined for this vehicle.'),
  inputs: [
    { key: 'v', label: 'Head velocity relative to the surface', unit: 'm/s', default: 6, min: 0.2, max: 30, group: 'Impact', help: 'From a seat-test occupant simulation; several m/s is typical of a head strike on the seat ahead' },
    { key: 'm_head', label: 'Effective head mass', unit: 'kg', default: 4.5, min: 1, max: 10, group: 'Impact', help: '4.5 kg is the customary 50th-percentile head-form mass' },
    { key: 'h_pad', label: 'Padding thickness', unit: 'm', default: 0.05, min: 0.002, group: 'Padding' },
    { key: 'sig_pad', label: 'Padding crush (plateau) stress', unit: 'Pa', default: 3.5e5, min: 1e3, group: 'Padding', help: 'Dynamic plateau stress of the foam or honeycomb; rate- and temperature-dependent, from test' },
    { key: 'E_pad', label: 'Padding initial modulus', unit: 'Pa', default: 8e6, min: 1e4, group: 'Padding' },
    { key: 'A_contact', label: 'Contact area', unit: 'm²', default: 0.008, min: 1e-4, group: 'Padding', help: 'Effective area engaged by the head (grows with indentation in reality)' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 3000, min: 200, max: 200000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c) { const v = c.meta.type === 'helicopter' ? 7 : 6, need = (0.5 * 4.5 * v * v) / (3.5e5 * 0.008); return { v, h_pad: Math.max(0.05, Math.ceil((1.25 * need) / 0.7 * 200) / 200) }; },
  run(i) {
    const r = headImpact(i, Math.round(i.nSteps)), ns = Math.min(r.t.length, 1500), tu = N.linspace(0, r.t[r.t.length - 1], ns), au = tu.map((t) => N.interp1(r.t, r.a, t)), h15 = hic(tu, au, 0.015), h36 = hic(tu, au, 0.036), pk = N.amax(r.a), warnings = [];
    const sNeed = (0.5 * i.m_head * i.v * i.v) / r.Fp, dur = tu[au.reduce((last, a, k) => (a > 0 ? k : last), 0)];
    if (r.bottom) warnings.push('The padding bottoms out (densifies): the deceleration spikes. Use thicker padding or a higher plateau stress.');
    warnings.push('HIC from a one-dimensional head-form model is a screening value; certification uses an instrumented dummy in a dynamic seat test.');
    return {
      kpis: [
        kpi('HIC', 'Head Injury Criterion (15 ms)', h15, '-', h15 <= 700 ? 'ok' : h15 <= 1000 ? 'warn' : 'bad', 'HIC 1000 is the pass limit in aircraft seat tests (36 ms window); 700 is the automotive HIC15 limit'),
        kpi('HIC36', 'Head Injury Criterion (36 ms)', h36, '-', h36 <= 1000 ? 'ok' : 'bad', 'Limit 1000'),
        kpi('head_peak_g', 'Peak head deceleration', pk, 'g'),
        kpi('pad_crush_m', 'Padding crush', r.xmax, 'm', r.bottom ? 'bad' : 'ok', `Densification at ${(0.7 * i.h_pad).toFixed(3)} m`),
        kpi('pad_needed_m', 'Crush needed at the plateau force', sNeed, 'm'),
        kpi('pulse_ms', 'Pulse duration', dur * 1e3, 'ms'),
        kpi('plateau_g', 'Plateau deceleration', r.Fp / (i.m_head * G0), 'g'),
      ],
      plots: [{ type: 'line', title: 'Head deceleration pulse', xlabel: 'Time [ms]', ylabel: 'Deceleration [g]', series: [{ name: 'Head', x: thin(tu.map((t) => t * 1e3)), y: thin(au) }] }],
      warnings,
      models: ['Single mass on an elastic–plastic pad with densification', 'HIC = max (t₂−t₁)·(mean acceleration)^2.5 over 15 and 36 ms windows'],
      assumptions: ['Normal impact, rigid backing structure, constant contact area', 'Rate-independent plateau stress supplied by the user', 'No neck load, rotation or glancing contact'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [500, 1000, 2000, 4000, 8000], metric: 'HIC' },
  verify() {
    const t = N.linspace(0, 0.03, 1501), sq = t.map((x) => (x >= 0.005 && x < 0.015 ? 60 : 0)), el = headImpact({ ...HEAD_BASE, sig_pad: 1e12, h_pad: 10, E_pad: 8e6 * 200, v: 2 }, 20000), w = Math.sqrt(el.k / 4.5);
    const pl = N.kv(head.run({ ...HEAD_BASE, E_pad: 5e8, nSteps: 20000 })), Tp = 6 / (pl.plateau_g * G0);
    return [
      N.check('HIC of a square pulse = T·A^2.5', hic(t, sq, 0.015), 0.01 * 60 ** 2.5, 5e-3, 'Definition (sampling error of the pulse edges only)'),
      N.check('Elastic pad: peak deceleration v·ω/g', N.amax(el.a), (2 * w) / G0, 1e-3, 'Half-sine response of a mass on a linear spring'),
      N.check('Elastic pad: head rebounds at the impact speed', el.rebound, 2, 1e-3, 'Energy conservation'),
      N.check('Stiff pad on a plateau: HIC → (v/a)·a^2.5', pl.HIC, Tp * pl.plateau_g ** 2.5, 0.02, 'Constant-deceleration pulse'),
    ];
  },
  calibration: { params: [{ key: 'sig_pad', min: 1e3, max: 1e7 }, { key: 'A_contact', min: 1e-4, max: 0.05 }], sweep: 'v', target: 'head_peak_g', note: 'Head-form drop-test peak decelerations at several speeds calibrate the pad plateau stress and effective contact area.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.HIC36 > 1000 || o.HIC > 700) out.push({ severity: o.HIC36 > 1000 ? 'critical' : 'warn', title: 'Head injury criterion exceeded or marginal', detail: `HIC15 ${o.HIC.toFixed(0)}, HIC36 ${o.HIC36.toFixed(0)}, peak ${o.head_peak_g.toFixed(0)} g.`, action: o.pad_crush_m > 0.7 * i.h_pad ? `Padding bottoms out: at least ${(o.pad_needed_m / 0.7 * 1e3).toFixed(0)} mm is needed at this plateau stress.` : 'Lower the plateau stress (softer foam, breakable seat-back feature) so the head stops over a longer distance, or increase seat pitch / add an upper-torso restraint to avoid the contact.', basis: 'HIC ≤ 1000 in dynamic seat tests' });
    else out.push({ severity: 'info', title: 'Head strike is within the criterion', detail: `HIC15 ${o.HIC.toFixed(0)} using ${(o.pad_crush_m * 1e3).toFixed(0)} mm of ${(i.h_pad * 1e3).toFixed(0)} mm padding.`, action: 'Padding thickness beyond what is crushed is unused mass and volume.', basis: 'HIC' });
    return out;
  },
};

// ---- 4. bird strike -------------------------------------------------------------------------
const bird = {
  id: 'bird', title: 'Bird strike loads (hydrodynamic theory)', fidelity: 'analytical',
  summary: 'Treats the bird as a fluid cylinder: shock (Hugoniot) and steady stagnation pressures, impact duration, average and peak force from the momentum transferred, and an energy-based screening of skin penetration.',
  equations: ['Conservation of linear momentum', 'Conservation of energy', 'Impact contact equations', 'Fracture criteria'],
  inputs: [
    { key: 'm_bird', label: 'Bird mass', unit: 'kg', default: 1.81, min: 0.01, max: 10, group: 'Bird', help: '1.81 kg (4 lb) is the usual transport-aircraft bird; 3.6 kg (8 lb) is used for some empennage requirements. Confirm against your certification basis' },
    { key: 'rho_bird', label: 'Bird density', unit: 'kg/m³', default: 950, min: 300, max: 1100, group: 'Bird', help: 'About 950 kg/m³ for gelatine substitutes with porosity' },
    { key: 'LD', label: 'Bird length / diameter', unit: '-', default: 2, min: 1, max: 4, group: 'Bird' },
    { key: 'V', label: 'Impact speed', unit: 'm/s', default: 180, min: 5, max: 400, group: 'Impact', help: 'Design cruise speed at sea level for transport aircraft' },
    { key: 'angle_deg', label: 'Angle between flight path and surface', unit: 'deg', default: 90, min: 5, max: 90, group: 'Impact', help: '90° = normal impact; leading edges and windshields are oblique' },
    { key: 'c0', label: 'Bird-material sound speed', unit: 'm/s', default: 1482, min: 300, max: 2000, group: 'Bird', help: 'Water value; porosity lowers the effective shock speed' },
    { key: 'k_h', label: 'Hugoniot slope', unit: '-', default: 2.0, min: 1, max: 3, group: 'Bird', help: 'Shock speed = c₀ + k·particle speed (≈ 2 for water)' },
    { key: 't_skin', label: 'Target skin thickness', unit: 'm', default: 0.003, min: 1e-4, group: 'Target' },
    MAT,
    { key: 'k_area', label: 'Deforming area / bird footprint (empirical)', unit: '-', default: 3, min: 1, max: 20, group: 'Target', help: 'Ratio of the skin area that stretches plastically to the bird cross-section; calibrate against gun tests' },
    { key: 'eps_f', label: 'Skin biaxial failure strain', unit: '-', default: 0.1, min: 0.01, max: 0.5, group: 'Target' },
    { key: 'eta', label: 'Fraction of normal kinetic energy absorbed by the skin (empirical)', unit: '-', default: 0.5, min: 0.05, max: 1, group: 'Target', help: 'The rest remains as bird flow energy and structural motion' },
  ],
  defaults(c, up) { const V = c.aero.Vmo_ms || c.flight.V_ms; return { V: Math.max(5, V), m_bird: c.mass.mtow_kg < 100 ? 0.45 : undefined, t_skin: up.fea?.t_skin_root_m ? Math.max(0.0008, 0.3 * up.fea.t_skin_root_m) : Math.min(c.struct.t_skin_mm, Math.max(0.5, 0.6 * (c.fuselage.dia_m || 1))) / 1e3, material: c.struct.material, angle_deg: 60 }; },
  run(i) {
    const m = mat(i.material), sinA = Math.sin(N.rad(i.angle_deg)), Vn = i.V * sinA, vol = i.m_bird / i.rho_bird, D = ((4 * vol) / (Math.PI * (i.LD - 1 / 3))) ** (1 / 3), L = i.LD * D, warnings = [];
    const Us = i.c0 + i.k_h * Vn, pH = i.rho_bird * Us * Vn, pS = 0.5 * i.rho_bird * Vn * Vn, Ab = (Math.PI * D * D) / 4, Leff = L + D / Math.tan(N.rad(Math.max(i.angle_deg, 5))) * (i.angle_deg < 89.9 ? 1 : 0), dur = Leff / i.V, imp = i.m_bird * Vn, Favg = imp / dur, Fpk = 2 * Favg;
    const KEn = 0.5 * i.m_bird * Vn * Vn, flow = 0.5 * (m.Sy + m.Su), cap = flow * i.eps_f * i.t_skin * ((i.k_area * Ab) / sinA), rf = cap / (i.eta * KEn), tReq = i.t_skin / rf, Vpen = i.V * Math.sqrt(rf);
    if (Vn < 50) warnings.push('Normal velocity below about 50 m/s: the bird no longer behaves as a fluid and the hydrodynamic theory over-simplifies the load.');
    if (rf < 1) warnings.push('The energy screening predicts skin rupture: the structure behind (spar, bulkhead, systems) must withstand the residual bird.');
    warnings.push('Penetration is an empirical energy screening with user factors; substantiation needs gun tests or validated explicit simulation.');
    const Vs = N.linspace(0.3 * i.V, 1.5 * i.V, 40);
    return {
      kpis: [
        kpi('bird_force_N', 'Peak impact force (triangular pulse)', Fpk, 'N'), kpi('bird_force_avg_N', 'Average impact force', Favg, 'N'),
        kpi('duration_s', 'Impact duration (squash-up time)', dur, 's'), kpi('impulse_Ns', 'Normal impulse', imp, 'N·s'),
        kpi('p_hugoniot_Pa', 'Initial shock (Hugoniot) pressure', pH, 'Pa'), kpi('p_stagnation_Pa', 'Steady-flow stagnation pressure', pS, 'Pa'),
        kpi('ke_normal_J', 'Kinetic energy of the normal velocity component', KEn, 'J'), kpi('ke_total_J', 'Total bird kinetic energy', 0.5 * i.m_bird * i.V ** 2, 'J'),
        kpi('bird_dia_m', 'Bird diameter', D, 'm'), kpi('bird_len_m', 'Bird length', L, 'm'),
        kpi('penetration_RF', 'Skin energy reserve factor (empirical screening)', rf, '-', rf >= 1 ? 'ok' : 'warn', 'Uncalibrated default factors: indicative only. ' + 'Membrane plastic work to rupture / absorbed share of normal kinetic energy'),
        kpi('t_required_m', 'Skin thickness for a reserve factor of 1', tReq, 'm'), kpi('V_penetration_ms', 'Estimated penetration speed', Vpen, 'm/s'),
      ],
      plots: [
        { type: 'line', title: 'Impact force versus speed', xlabel: 'Impact speed [m/s]', ylabel: 'Force [kN]', series: [{ name: 'Peak', x: Vs, y: Vs.map((v) => (2 * i.m_bird * v * v * sinA) / Leff / 1e3) }, { name: 'Average', x: Vs, y: Vs.map((v) => (i.m_bird * v * v * sinA) / Leff / 1e3) }], annotations: [{ x: i.V, label: 'Design' }] },
        { type: 'line', title: 'Impact pressures versus speed', xlabel: 'Impact speed [m/s]', ylabel: 'Pressure [MPa]', ylog: true, series: [{ name: 'Hugoniot (shock)', x: Vs, y: Vs.map((v) => (i.rho_bird * (i.c0 + i.k_h * v * sinA) * v * sinA) / 1e6) }, { name: 'Stagnation (steady flow)', x: Vs, y: Vs.map((v) => (0.5 * i.rho_bird * (v * sinA) ** 2) / 1e6) }] },
        { type: 'line', title: 'Idealised force pulse', xlabel: 'Time [ms]', ylabel: 'Force [kN]', series: [{ name: 'Triangular pulse of equal impulse', x: [0, (dur * 1e3) / 4, dur * 1e3], y: [0, Fpk / 1e3, 0] }] },
      ],
      warnings,
      models: ['Hydrodynamic (fluid cylinder) bird model after Wilbeck', 'Linear shock Hugoniot Us = c₀ + k·up', 'Momentum-based force with a triangular pulse', 'Membrane plastic-work penetration screening (empirical)'],
      assumptions: ['Rigid target for the load calculation; target compliance lowers the real peak force', 'Only the normal momentum is transferred; the tangential part slides off', 'Bird geometry: cylinder with hemispherical ends', 'Bird sound speed 1482 m/s, Hugoniot slope 2 and length/diameter 2 are the customary water-like values, typical rather than verified against a primary source here; the penetration factors (deforming area, absorbed energy share, failure strain) are illustrative until calibrated against gun tests', DATA_NOTE],
    };
  },
  verify() {
    const b = { m_bird: 1.81, rho_bird: 950, LD: 2, V: 150, angle_deg: 90, c0: 1482, k_h: 2, t_skin: 0.003, material: 'Al 2024-T3', k_area: 3, eps_f: 0.1, eta: 0.5 }, r = N.kv(bird.run(b)), D = r.bird_dia_m;
    const lo = N.kv(bird.run({ ...b, V: 5.0001 }));
    return [
      N.check('Impulse = average force × duration = m·V', r.bird_force_avg_N * r.duration_s, 1.81 * 150, 1e-12, 'Momentum conservation'),
      N.check('Bird volume from the cylinder with hemispherical ends', (Math.PI * D ** 3) / 6 + (Math.PI * D * D * D) / 4, 1.81 / 950, 1e-10, 'Geometry: sphere + cylinder of length D'),
      N.check('Hugoniot pressure → acoustic ρ·c₀·V at low speed', lo.p_hugoniot_Pa / (950 * 5.0001), 1482 + 2 * 5.0001, 1e-12, 'Rankine–Hugoniot momentum jump'),
      N.check('Stagnation pressure ½ρV²', r.p_stagnation_Pa, 0.5 * 950 * 150 ** 2, 1e-12, 'Bernoulli'),
    ];
  },
  calibration: { params: [{ key: 'k_area', min: 1, max: 20 }, { key: 'eta', min: 0.05, max: 1 }], sweep: 't_skin', target: 'V_penetration_ms', note: 'Gun-test penetration speeds versus skin thickness calibrate the two empirical penetration factors.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.penetration_RF < 1) out.push({ severity: 'advise', title: 'Skin penetration is likely', detail: `Reserve factor ${o.penetration_RF.toFixed(2)}; about ${(o.t_required_m * 1e3).toFixed(1)} mm of this alloy would be needed.`, action: 'Rather than thickening the whole skin, add a splitter or deflector, use a tougher leading-edge material (fibre-metal laminate), or protect the structure and systems behind the skin.', basis: 'Empirical membrane-energy screening with uncalibrated default factors; the requirement is continued safe flight after the strike, not an intact skin' });
    out.push({ severity: 'info', title: 'Load for the supporting structure', detail: `Peak force ${(o.bird_force_N / 1e3).toFixed(0)} kN over ${(o.duration_s * 1e3).toFixed(1)} ms.`, action: 'Apply this pulse to the rib, spar or frame model in Suite 2 or 10; oblique surfaces cut the normal load with sin(angle), so sweep is a mass-free protection lever.', basis: 'Momentum transfer' });
    return out;
  },
};

// ---- 5. 1-D elastic–plastic wave propagation ------------------------------------------------
/** Explicit FE bar (linear elements, lumped mass, central differences) striking a rigid wall with velocity v0. */
function barImpact(i, mm, ctx) {
  const n = Math.max(4, Math.round(i.nElem)), h = i.L / n, c = Math.sqrt(mm.E / mm.rho), dt = (N.clamp(i.cfl, 0.05, 1) * h) / c, tEnd = (i.t_factor * 2 * i.L) / c, ns = Math.ceil(tEnd / dt), jc = jcOf(mm), model = i.model;
  const u = new Array(n + 1).fill(0), v = new Array(n + 1).fill(-i.v0), sig = new Array(n).fill(0), ep = new Array(n).fill(0), eps = new Array(n).fill(0), Tm = new Array(n).fill(293), ms = N.range(n + 1, (j) => mm.rho * h * (j === 0 || j === n ? 0.5 : 1));
  const sy = (e, rate, T) => (model === 'elastic' ? Infinity : model === 'bilinear' ? mm.Sy + i.H * e : johnsonCook(jc, e, rate, T));
  let contact = true, tRel = NaN, wp = 0, pkWall = 0, eErr = 0; const E0 = 0.5 * mm.rho * i.L * i.v0 ** 2, rows = [], tRow = [], wallT = [], wallS = [], freeV = [], everyR = Math.max(1, Math.floor(ns / 60)), everyW = Math.max(1, Math.floor(ns / 300)), cols = Math.min(n, 80);
  let t = 0, dts = dt;
  for (let s = 0; s <= ns; s++) {
    // stress update from the current strain
    for (let e = 0; e < n; e++) {
      const en = (u[e + 1] - u[e]) / h, de = en - eps[e], rate = Math.abs(de) / dts; eps[e] = en; let st = sig[e] + mm.E * de;
      if (model !== 'elastic') {
        let y = sy(ep[e], rate, Tm[e]);
        if (Math.abs(st) > y) { let dp = 0; for (let k = 0; k < 5; k++) { const f = Math.abs(st) - mm.E * dp - sy(ep[e] + dp, rate, Tm[e]); if (k && Math.abs(f) < 1e-7 * y) break; const dy = (sy(ep[e] + dp + 1e-6, rate, Tm[e]) - sy(ep[e] + dp, rate, Tm[e])) / 1e-6; dp = Math.max(0, dp + f / (mm.E + Math.max(dy, 0))); }
          st -= Math.sign(st) * mm.E * dp; ep[e] += dp; wp += Math.abs(st) * dp * h; if (model === 'Johnson–Cook') Tm[e] += (0.9 * Math.abs(st) * dp) / (mm.rho * mm.cp); }
      }
      sig[e] = st;
    }
    if (contact && sig[0] > 0 && s > 2) { contact = false; tRel = t; dts = Math.min(dt, (0.98 * h) / c); } // a free–free lumped bar is only marginally stable at Courant 1
    pkWall = Math.max(pkWall, -sig[0]);
    let ke = 0, el = 0; for (let j = 0; j <= n; j++) ke += 0.5 * ms[j] * v[j] * v[j]; for (let e = 0; e < n; e++) el += (0.5 * sig[e] * sig[e] * h) / mm.E;
    if (s > 0) eErr = Math.max(eErr, Math.abs(ke + el + wp - E0) / E0);
    if (s % everyW === 0) { wallT.push(t); wallS.push(contact ? -sig[0] : 0); freeV.push(v[n]); }
    if (s % everyR === 0 && rows.length < 61) { tRow.push(t); rows.push(N.range(cols, (k) => sig[Math.min(n - 1, Math.floor(((k + 0.5) * n) / cols))] / 1e6)); }
    if (s === ns) break;
    for (let j = 0; j <= n; j++) { const f = (j < n ? sig[j] : 0) - (j > 0 ? sig[j - 1] : 0); v[j] += (dts * f) / ms[j]; }
    if (contact) v[0] = 0;
    for (let j = 0; j <= n; j++) u[j] += dts * v[j];
    t += dts;
    if (s % 500 === 0) ctx?.progress?.(s / ns, 'Wave propagation');
  }
  // permanent shortening from the plastic strains alone (elastic vibration continues after rebound)
  return { c, dt, ns, n, tRel, pkWall, eErr, wp, E0, short: N.sum(ep) * h, epMax: N.amax(ep), dT: N.amax(Tm) - 293, vOut: N.sum(v.map((x, j) => x * ms[j])) / (mm.rho * i.L), rows, tRow, cols, wallT, wallS, freeV, jc, ep, h };
}
const WAVE_BASE = { material: 'Al 2024-T3', L: 0.1, v0: 10, model: 'elastic', H: 1e9, t_factor: 3, cs_D: 6500, cs_q: 4, cfl: 1, nElem: 200 };
const wave = {
  id: 'wave', title: 'Stress-wave impact of a bar with rate-dependent plasticity', fidelity: 'numerical',
  summary: 'An explicit finite-element bar strikes a rigid wall (Taylor-type impact in one dimension): elastic and plastic wave fronts, impact stress, contact time, permanent shortening, adiabatic heating and energy balance, with Johnson–Cook and Cowper–Symonds flow-stress curves.',
  equations: ['Conservation of linear momentum', 'Conservation of energy', 'Transient structural dynamics equations', 'Plasticity constitutive equations', 'Strain-rate-dependent material equations', 'Impact contact equations'],
  inputs: [
    MAT,
    { key: 'L', label: 'Bar length', unit: 'm', default: 0.1, min: 0.001, group: 'Impact' },
    { key: 'v0', label: 'Impact velocity', unit: 'm/s', default: 120, min: 0.01, max: 600, group: 'Impact' },
    { key: 'model', label: 'Material model', type: 'select', options: ['Johnson–Cook', 'bilinear', 'elastic'], default: 'Johnson–Cook', group: 'Material' },
    { key: 'H', label: 'Plastic hardening modulus (bilinear model)', unit: 'Pa', default: 1e9, min: 0, group: 'Material' },
    { key: 'cs_D', label: 'Cowper–Symonds D (comparison curve)', unit: '1/s', default: 6500, min: 1, group: 'Material' },
    { key: 'cs_q', label: 'Cowper–Symonds q (comparison curve)', unit: '-', default: 4, min: 1, max: 20, group: 'Material' },
    { key: 't_factor', label: 'Simulated time / wave transit time 2L/c', unit: '-', default: 3, min: 0.5, max: 20, group: 'Numerics' },
    { key: 'cfl', label: 'Courant number', unit: '-', default: 0.95, min: 0.05, max: 1, group: 'Numerics', help: 'Explicit stability requires ≤ 1; exactly 1 reproduces elastic waves without dispersion on a uniform mesh' },
    { key: 'nElem', label: 'Elements', unit: '', default: 200, min: 10, max: 4000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c) { const m = mat(c.struct.material), vy = m.Sy / (m.rho * Math.sqrt(m.E / m.rho)); return { material: c.struct.material, v0: Math.max(4 * vy, c.aero.Vmo_ms || c.flight.V_ms || 50), cs_D: m.rho < 3500 ? 6500 : 1e5 }; },
  run(i, ctx) {
    const m = mat(i.material), r = barImpact(i, m, ctx), sEl = m.rho * r.c * i.v0, vy = m.Sy / (m.rho * r.c), warnings = [], jc = r.jc;
    if (i.model === 'Johnson–Cook' && jc.fitted) warnings.push('No Johnson–Cook constants are tabulated for this material: a simple fit through yield and ultimate strength with C = 0.01 was used. Replace with split-Hopkinson-bar data.');
    if (r.epMax > 0.5) warnings.push('Plastic strain above 50%: the small-strain, uniaxial-stress bar ignores mushrooming (radial flow) of a real Taylor specimen and overstates axial stress.');
    if (!Number.isFinite(r.tRel)) warnings.push('The bar has not rebounded within the simulated time.');
    const ee = N.linspace(0, 0.3, 40), rates = [1, 1e2, 1e4], xs = N.range(r.cols, (k) => ((k + 0.5) * i.L * 1e3) / r.cols);
    return {
      kpis: [
        kpi('wave_speed_ms', 'Elastic bar wave speed √(E/ρ)', r.c, 'm/s'),
        kpi('sigma_elastic_Pa', 'Elastic impact stress ρ·c·v', sEl, 'Pa', undefined, sEl < m.Sy ? 'Below yield: the impact is elastic' : 'Exceeds yield: a plastic wave forms'),
        kpi('sigma_wall_Pa', 'Peak stress at the wall (simulated)', r.pkWall, 'Pa'),
        kpi('v_yield_ms', 'Impact velocity that first causes yield', vy, 'm/s'),
        kpi('contact_s', 'Contact duration', r.tRel, 's', undefined, `Elastic value 2L/c = ${((2 * i.L) / r.c * 1e6).toFixed(1)} µs`),
        kpi('shortening_m', 'Permanent shortening', r.short, 'm'), kpi('final_length_ratio', 'Final / initial length', 1 - r.short / i.L, '-'),
        kpi('eps_p_max', 'Peak plastic strain', r.epMax, '-'), kpi('dT_max_K', 'Peak adiabatic temperature rise', r.dT, 'K'),
        kpi('plastic_work_frac', 'Share of kinetic energy dissipated plastically', r.wp / r.E0, '-'),
        kpi('rebound_ms', 'Mean rebound velocity', r.vOut, 'm/s'),
        kpi('energy_err', 'Energy balance error', r.eErr, '-', r.eErr < 0.01 ? 'ok' : 'warn'),
        kpi('dt_s', 'Stable time step used', r.dt, 's'),
      ],
      plots: [
        { type: 'heat', title: 'Stress wave diagram (compression negative)', xlabel: 'Position from the wall [mm]', ylabel: 'Time [µs]', zlabel: 'Stress [MPa]', x: xs, y: r.tRow.map((t) => t * 1e6), z: r.rows, diverging: true },
        { type: 'line', title: 'Wall stress history', xlabel: 'Time [µs]', ylabel: 'Compressive stress [MPa]', series: [{ name: 'Contact stress', x: r.wallT.map((t) => t * 1e6), y: r.wallS.map((s) => s / 1e6) }], annotations: [{ y: sEl / 1e6, label: 'ρ·c·v' }] },
        { type: 'line', title: 'Free-end velocity history', xlabel: 'Time [µs]', ylabel: 'Velocity [m/s]', series: [{ name: 'Free end', x: r.wallT.map((t) => t * 1e6), y: r.freeV }] },
        { type: 'line', title: 'Dynamic flow-stress curves', xlabel: 'Plastic strain [-]', ylabel: 'Flow stress [MPa]', series: [...rates.map((rt) => ({ name: `Johnson–Cook, ${rt.toExponential(0)} 1/s`, x: ee, y: ee.map((e) => johnsonCook(jc, e, rt) / 1e6) })), ...[1e2, 1e4].map((rt) => ({ name: `Cowper–Symonds × JC static, ${rt.toExponential(0)} 1/s`, x: ee, y: ee.map((e) => ((jc.A + jc.B * e ** jc.n) * cowperSymonds(rt, i.cs_D, i.cs_q)) / 1e6), style: 'dash' }))] },
      ],
      warnings,
      models: [`Explicit central-difference bar FE (${r.n} elements, ${r.ns} steps, Courant ${i.cfl})`, i.model === 'Johnson–Cook' ? 'Johnson–Cook flow stress with adiabatic heating (Taylor–Quinney 0.9)' : i.model === 'bilinear' ? 'Rate-independent bilinear plasticity' : 'Linear elastic', 'Rigid frictionless wall with release on tension'],
      assumptions: ['Uniaxial stress, small-strain kinematics, no radial inertia', 'Strain rate in the flow stress taken as the total element strain rate', 'No fracture or erosion criterion', 'Cowper–Symonds constants of the comparison curve are typical textbook values (aluminium D = 6500 1/s, q = 4), illustrative only', DATA_NOTE],
    };
  },
  convergence: { param: 'nElem', label: 'Elements', levels: [25, 50, 100, 200], metric: 'shortening_m' },
  verify() {
    const m = mat('Al 2024-T3'), c = Math.sqrt(m.E / m.rho), el = N.kv(wave.run({ ...WAVE_BASE })), H = 2e9, v0 = 60, pl = N.kv(wave.run({ ...WAVE_BASE, v0, model: 'bilinear', H, cfl: 0.9, nElem: 400, t_factor: 0.8 }));
    const Et = (m.E * H) / (m.E + H), vy = m.Sy / (m.rho * c), sPl = m.Sy + m.rho * Math.sqrt(Et / m.rho) * (v0 - vy);
    return [
      N.check('Wave speed c = √(E/ρ)', el.wave_speed_ms, c, 1e-12, 'D’Alembert solution of the 1-D wave equation'),
      N.check('Elastic impact stress σ = ρ·c·v', el.sigma_wall_Pa, m.rho * c * 10, 1e-6, 'Momentum jump across the wave front (exact at Courant number 1)'),
      N.check('Contact duration 2L/c', el.contact_s, 0.2 / c, 0.01, 'Reflection of the release wave from the free end'),
      N.check('Elastic rebound at the impact speed', el.rebound_ms, 10, 5e-3, 'Energy and momentum conservation (the half-element mass at the wall node is arrested inelastically: O(1/n))'),
      N.check('Elastic energy conservation', el.energy_err, 0, 0.02, 'Explicit central differences; errors limited to the arrested wall-node mass, 1/(2n), and the time-step change at release'),
      N.check('Plastic wave stress σy + ρ·c_p·(v − v_y), bilinear hardening', pl.sigma_wall_Pa, sPl, 0.02, 'von Kármán–Taylor rate-independent plastic wave theory'),
      N.check('Energy balance with plasticity', pl.energy_err, 0, 5e-3, 'Kinetic + elastic + plastic work = initial kinetic energy'),
    ];
  },
  calibration: { params: [{ key: 'H', min: 0, max: 2e10 }], sweep: 'v0', target: 'final_length_ratio', note: 'Taylor-impact final lengths at several velocities calibrate the hardening (and, externally, the Johnson–Cook constants in the material table).' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.sigma_elastic_Pa > mat(i.material).Sy) out.push({ severity: 'info', title: 'Impact velocity exceeds the elastic limit', detail: `Yield starts at ${o.v_yield_ms.toFixed(1)} m/s; at ${i.v0} m/s ${(100 * o.plastic_work_frac).toFixed(0)}% of the kinetic energy is dissipated plastically.`, action: 'For hail, debris, tool-drop or blade-fragment assessments use this velocity threshold as a first screen, then hand the case to an explicit FE solver with a fracture criterion.', basis: 'σ = ρ·c·v against yield' });
    if (o.dT_max_K > 100) out.push({ severity: 'advise', title: 'Significant adiabatic heating', detail: `Temperature rise up to ${o.dT_max_K.toFixed(0)} K.`, action: 'Thermal softening and shear banding become important; use temperature-dependent material data.', basis: 'Taylor–Quinney plastic work heating' });
    return out;
  },
};

// ---- 6. water impact (ditching) -------------------------------------------------------------
function wedgeImpact(i, n) {
  const tb = Math.tan(N.rad(i.beta_deg)), gam = i.theory === 'Wagner' ? Math.PI / 2 : 1, mu = (Math.PI * i.rho_w * i.L_hull) / (2 * i.mass), zc = (i.b_max * tb) / gam, tEnd = (4 * zc) / i.v0, g = i.gravity ? G0 : 0;
  // momentum form: d/dt[(M + ma)·V] = M·g, with ma = ½πρ c² L and c = γ z / tanβ (capped at the chine)
  const cOf = (z) => Math.min((gam * Math.max(z, 0)) / tb, i.b_max), f = (t, y) => { const c = cOf(y[0]), V = y[1] / (i.mass * (1 + mu * c * c)); return [V, i.mass * g]; };
  const r = N.rk4(f, 0, [0, i.mass * i.v0], tEnd, n), T = r.t, Z = r.y.map((y) => y[0]), C = Z.map(cOf), V = r.y.map((y, k) => y[1] / (i.mass * (1 + mu * C[k] * C[k])));
  // force on the hull = M·g − M·dV/dt = (dma/dt·V + μc²·M·g)/(1 + μc²); the added mass stops growing once the chine is wet
  const F = V.map((v, k) => (C[k] < i.b_max ? (Math.PI * i.rho_w * i.L_hull * C[k] * (gam / tb) * v * v + mu * C[k] * C[k] * i.mass * g) / (1 + mu * C[k] * C[k]) : 0)), kp = N.argmax(F);
  return { T, Z, V, F, C, kp, tb, gam, mu, zc };
}
const DITCH_BASE = { mass: 5000, v0: 3, beta_deg: 20, L_hull: 5, b_max: 5, rho_w: 1025, theory: 'von Kármán', gravity: false, nSteps: 2000 };
const ditch = {
  id: 'ditch', title: 'Water impact of the hull bottom (ditching)', fidelity: 'reduced-order',
  summary: 'Vertical water entry of a wedge-shaped hull by added-mass momentum theory (von Kármán, with optional Wagner wetting correction): impact force and deceleration histories, peak pressure and the effect of deadrise angle.',
  equations: ['Conservation of linear momentum', 'Conservation of energy', 'Impact contact equations'],
  inputs: [
    { key: 'mass', label: 'Aircraft mass', unit: 'kg', default: 70000, min: 0.01, group: 'Impact' },
    { key: 'v0', label: 'Vertical velocity at water contact', unit: 'm/s', default: 1.5, min: 0.05, max: 30, group: 'Impact', help: 'A controlled ditching touches down at 1–2 m/s; a crash can be far higher' },
    { key: 'beta_deg', label: 'Deadrise angle of the bottom', unit: 'deg', default: 20, min: 3, max: 60, group: 'Hull', help: 'Angle between the bottom and the water; a round fuselage belly is locally shallow (5–15°), which gives high loads' },
    { key: 'L_hull', label: 'Length of bottom entering the water', unit: 'm', default: 12, min: 0.01, group: 'Hull' },
    { key: 'b_max', label: 'Half-breadth at the chine', unit: 'm', default: 1.8, min: 0.005, group: 'Hull', help: 'Wetted half-width at which the bottom stops widening' },
    { key: 'rho_w', label: 'Water density', unit: 'kg/m³', default: 1025, min: 900, max: 1100, group: 'Impact' },
    { key: 'theory', label: 'Wetted-width theory', type: 'select', options: ['von Kármán', 'Wagner'], default: 'Wagner', group: 'Method', help: 'Wagner includes the pile-up of water (wetted width × π/2) and gives higher, more realistic loads for small deadrise' },
    { key: 'gravity', label: 'Include weight during the impact', type: 'bool', default: true, group: 'Method' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 2000, min: 50, max: 200000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c) => ({ mass: c.mass.mtow_kg, v0: Math.max(0.5, 0.5 * c.gear.sink_ms), L_hull: Math.max(0.05, 0.35 * c.fuselage.len_m), b_max: Math.max(0.01, 0.45 * c.fuselage.dia_m) }),
  run(i) {
    const r = wedgeImpact(i, Math.max(50, Math.round(i.nSteps))), kp = r.kp, Fp = r.F[kp], warnings = [], Vp = r.V[kp], pMax = 0.5 * i.rho_w * Vp * Vp * (1 + Math.PI ** 2 / (4 * r.tb * r.tb)), pMean = Fp / (2 * Math.max(r.C[kp], 1e-9) * i.L_hull);
    const chine = r.C[kp] >= i.b_max * 0.999 || r.F.slice(kp + 1).every((f) => f === 0) && r.C[Math.min(kp + 1, r.C.length - 1)] >= i.b_max;
    if (chine) warnings.push('The chine is immersed before the momentum-theory force peak: the peak load is set by the hull breadth, and later loads (not modelled) come from steady planing and buoyancy.');
    if (i.beta_deg < 8) warnings.push('Deadrise below about 8°: air cushioning and hull flexibility limit the real pressures; rigid-wedge theory over-predicts.');
    const bs = N.linspace(5, 45, 30), k = thin(N.range(r.T.length), 300);
    return {
      kpis: [
        kpi('ditch_force_N', 'Peak water impact force', Fp, 'N'),
        kpi('ditch_peak_g', 'Peak deceleration', Fp / (i.mass * G0), 'g', Fp / (i.mass * G0) < 8 ? 'ok' : 'warn'),
        kpi('t_peak_s', 'Time of peak force', r.T[kp], 's'), kpi('depth_peak_m', 'Immersion at peak force', r.Z[kp], 'm'), kpi('wetted_half_width_m', 'Wetted half-width at peak', r.C[kp], 'm'),
        kpi('p_peak_Pa', 'Peak local pressure (Wagner spray-root)', pMax, 'Pa', undefined, '½ρV²(1 + π²/(4·tan²β)) at the velocity when the force peaks'),
        kpi('p_mean_Pa', 'Mean pressure on the wetted bottom at peak', pMean, 'Pa'),
        kpi('v_after_ms', 'Vertical velocity at chine immersion or end', r.V[r.V.length - 1], 'm/s'),
        kpi('added_mass_ratio', 'Added mass / aircraft mass at the chine', r.mu * i.b_max ** 2, '-'),
      ],
      plots: [
        { type: 'line', title: 'Water impact force history', xlabel: 'Time [ms]', ylabel: 'Force [kN]', series: [{ name: i.theory, x: k.map((j) => r.T[j] * 1e3), y: k.map((j) => r.F[j] / 1e3) }] },
        { type: 'line', title: 'Vertical velocity during entry', xlabel: 'Time [ms]', ylabel: 'Velocity [m/s]', series: [{ name: 'Velocity', x: k.map((j) => r.T[j] * 1e3), y: k.map((j) => r.V[j]) }] },
        { type: 'line', title: 'Peak deceleration versus deadrise angle', xlabel: 'Deadrise angle [deg]', ylabel: 'Peak deceleration [g]', series: [{ name: i.theory, x: bs, y: bs.map((b) => { const q = wedgeImpact({ ...i, beta_deg: b }, 400); return q.F[q.kp] / (i.mass * G0); }) }], annotations: [{ x: i.beta_deg, label: 'Design' }] },
      ],
      warnings,
      models: [`${i.theory} added-mass momentum theory for a rigid wedge`, 'RK4 integration of the momentum equation', 'Wagner peak pressure estimate'],
      assumptions: ['Vertical entry only: forward speed, trim, suction on the rear fuselage and porpoising are not modelled', 'Two-dimensional strip added mass with no aspect-ratio (end-flow) correction: conservative for short hulls', 'Rigid hull; no buoyancy, spray or air-cushion effects'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [100, 200, 400, 800, 1600], metric: 'ditch_force_N' },
  verify() {
    const r = N.kv(ditch.run(DITCH_BASE)), tb = Math.tan(N.rad(20)), mu = (Math.PI * 1025 * 5) / (2 * 5000), cs = 1 / Math.sqrt(5 * mu), Fx = (Math.PI * 1025 * 5 * 9 * cs * (5 / 6) ** 3) / tb, w = N.kv(ditch.run({ ...DITCH_BASE, theory: 'Wagner' }));
    return [
      N.check('von Kármán peak force (5/6)³·πρL·V₀²·c*/tanβ', r.ditch_force_N, Fx, 1e-3, 'von Kármán (1929) momentum theory, maximum at μc² = 1/5'),
      N.check('Wetted half-width at peak c* = 1/√(5μ)', r.wetted_half_width_m, cs, 0.02, 'von Kármán (1929)'),
      N.check('Immersion at peak z* = c*·tanβ', r.depth_peak_m, cs * tb, 0.02, 'von Kármán (1929); tolerance reflects the flat maximum sampled on the time grid'),
      N.check('Wagner wetting raises the peak force by π/2', w.ditch_force_N / r.ditch_force_N, Math.PI / 2, 2e-3, 'Wagner (1932) wetted-width factor'),
    ];
  },
  calibration: { params: [{ key: 'beta_deg', min: 3, max: 60 }], sweep: 'v0', target: 'ditch_peak_g', note: 'Model-scale drop or ditching test decelerations versus sink speed calibrate the effective deadrise angle of a non-wedge bottom.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.ditch_peak_g > 8) out.push({ severity: 'warn', title: 'High water-impact deceleration', detail: `${o.ditch_peak_g.toFixed(1)} g at ${i.v0} m/s sink speed.`, action: 'Reduce touchdown sink rate (procedures, flap setting), or increase local deadrise / add a keel strake: peak load falls roughly with 1/tan(deadrise).', basis: 'Added-mass momentum theory' });
    out.push({ severity: 'info', title: 'Bottom skin pressure', detail: `Peak ${(o.p_peak_Pa / 1e3).toFixed(0)} kPa, mean ${(o.p_mean_Pa / 1e3).toFixed(0)} kPa.`, action: 'Check lower-fuselage skin panels and frames for this pressure in Suite 2; bottom rupture governs flotation time and evacuation.', basis: 'Wagner pressure' });
    return out;
  },
};

// ---- 7. emergency-landing inertia loads -----------------------------------------------------
const emergency = {
  id: 'emergency', title: 'Emergency-landing inertia loads on an item of mass', fidelity: 'analytical',
  summary: 'Static ultimate inertia load factors applied to a seat, battery pack, tank or equipment item on a four-bolt attachment: bolt tension and shear for each direction, combined utilisation and the governing case.',
  equations: ['Conservation of linear momentum', 'Conservation of angular momentum'],
  inputs: [
    { key: 'mass', label: 'Item mass (with contents / occupant)', unit: 'kg', default: 100, min: 0.001, group: 'Item' },
    { key: 'h_cg', label: 'Centre of gravity above the attachment plane', unit: 'm', default: 0.4, min: 0, group: 'Item' },
    { key: 'lx', label: 'Bolt spacing fore–aft', unit: 'm', default: 0.5, min: 0.005, group: 'Attachment' },
    { key: 'ly', label: 'Bolt spacing lateral', unit: 'm', default: 0.45, min: 0.005, group: 'Attachment' },
    { key: 'd_bolt', label: 'Bolt diameter', unit: 'm', default: 0.00635, min: 0.001, group: 'Attachment' },
    { key: 'Su_bolt', label: 'Bolt ultimate tensile strength', unit: 'Pa', default: 860e6, min: 1e7, group: 'Attachment', help: 'About 860 MPa (125 ksi) for standard aerospace steel bolts' },
    { key: 'fitting', label: 'Fitting factor', unit: '-', default: 1.33, min: 1, max: 2, group: 'Attachment', help: '1.33 is customary for seat and restraint attachments' },
    { key: 'g_fwd', label: 'Forward load factor', unit: 'g', default: 9, min: 0, max: 60, group: 'Ultimate inertia factors', help: 'Static emergency-landing values: transport and light aeroplanes 9 g forward, 3 g up; rotorcraft 16 g forward, 4 g up, 8 g side, 20 g down. Confirm every factor against the applicable certification basis' },
    { key: 'g_up', label: 'Upward load factor', unit: 'g', default: 3, min: 0, max: 60, group: 'Ultimate inertia factors' },
    { key: 'g_side', label: 'Sideward load factor', unit: 'g', default: 3, min: 0, max: 60, group: 'Ultimate inertia factors', help: 'Transport: 3 g on airframe, 4 g on seats; light aeroplanes 1.5 g' },
    { key: 'g_down', label: 'Downward load factor', unit: 'g', default: 6, min: 0, max: 60, group: 'Ultimate inertia factors' },
    { key: 'g_aft', label: 'Rearward load factor', unit: 'g', default: 1.5, min: 0, max: 60, group: 'Ultimate inertia factors' },
  ],
  defaults(c) {
    const rot = c.meta.type === 'helicopter' || c.meta.type === 'evtol', light = c.mass.mtow_kg < 5700, occ = occupants(c) > 0, elec = c.prop.type === 'electric', m = occ ? (elec ? Math.max(10, (c.systems.batt_kWh * 1000) / 200 / 4) : 77 + 12) : Math.max(0.05, c.mass.payload_kg || 0.1 * c.mass.mtow_kg);
    const s = N.clamp(0.12 * m ** (1 / 3), 0.03, 1.2), o = { mass: m, h_cg: 0.6 * s, lx: s, ly: 0.9 * s, d_bolt: N.clamp(0.0005 * m ** 0.5 * (rot ? 1.5 : 1), 0.002, 0.02) };
    return rot ? { ...o, g_fwd: 16, g_up: 4, g_side: 8, g_down: 20, g_aft: 1.5 } : light ? { ...o, g_side: 1.5 } : o;
  },
  run(i) {
    const W = i.mass * G0, A = (Math.PI * i.d_bolt ** 2) / 4 * 0.75, Pt = i.Su_bolt * A, Ps = 0.6 * i.Su_bolt * (Math.PI * i.d_bolt ** 2) / 4, warnings = [];
    // each case: horizontal force (x or y) at height h and vertical force (positive = pulls the item off the floor)
    const cases = [['Forward', i.g_fwd, 0, 0], ['Rearward', -i.g_aft, 0, 0], ['Sideward', 0, i.g_side, 0], ['Upward', 0, 0, i.g_up], ['Downward', 0, 0, -i.g_down], ['Forward + sideward (resultant)', i.g_fwd * Math.SQRT1_2, i.g_side * Math.SQRT1_2, 0]];
    const rows = cases.map(([name, gx, gy, gz]) => {
      const Fx = gx * W * i.fitting, Fy = gy * W * i.fitting, Fz = gz * W * i.fitting, T = Fz / 4 + Math.abs(Fx) * i.h_cg / (2 * i.lx) + Math.abs(Fy) * i.h_cg / (2 * i.ly), S = Math.hypot(Fx, Fy) / 4, Rt = Math.max(T, 0) / Pt, Rs = S / Ps;
      // interaction Rs³ + Rt² = 1: utilisation as the load multiplier λ solving (λRs)³ + (λRt)² = 1
      let lo = 0, hi = 1e6; for (let k = 0; k < 80; k++) { const mid = 0.5 * (lo + hi); if ((mid * Rs) ** 3 + (mid * Rt) ** 2 < 1) lo = mid; else hi = mid; }
      return { name, T, S, Rt, Rs, util: Rs + Rt > 0 ? 1 / lo : 0, comp: Math.max(0, -Fz / 4 + Math.abs(Fx) * i.h_cg / (2 * i.lx) + Math.abs(Fy) * i.h_cg / (2 * i.ly)) };
    });
    const g = rows[N.argmax(rows.map((r) => r.util))], mos = g.util > 0 ? 1 / g.util - 1 : Infinity, dReq = i.d_bolt * Math.sqrt(Math.max(g.util, 1e-12));
    if (i.h_cg > 1.5 * Math.min(i.lx, i.ly)) warnings.push('The centre of gravity is high relative to the bolt pattern: overturning dominates and the floor structure sees large local pull-out loads.');
    warnings.push('The load factors are static ultimate values entered by the user; dynamic seat tests and the complete certification basis are not replaced by this check.');
    return {
      kpis: [
        kpi('bolt_util', `Governing bolt utilisation (${g.name})`, g.util, '-', g.util <= 0.87 ? 'ok' : g.util <= 1 ? 'warn' : 'bad', 'Interaction Rs³ + Rt² = 1 with the fitting factor applied'),
        kpi('mos_attach', 'Margin of safety of the attachment', mos, '-', mos >= 0.15 ? 'ok' : mos >= 0 ? 'warn' : 'bad'),
        kpi('bolt_tension_N', 'Peak bolt tension', Math.max(...rows.map((r) => r.T)), 'N'), kpi('bolt_shear_N', 'Peak bolt shear', Math.max(...rows.map((r) => r.S)), 'N'),
        kpi('floor_comp_N', 'Peak compressive load into the floor at one fitting', Math.max(...rows.map((r) => r.comp)), 'N'),
        kpi('d_bolt_req_m', 'Bolt diameter for unit utilisation (approx.)', dReq, 'm'),
        kpi('inertia_fwd_N', 'Forward inertia force (ultimate, unfactored)', i.g_fwd * W, 'N'),
        kpi('bolt_Pt_N', 'Bolt tensile strength', Pt, 'N'), kpi('bolt_Ps_N', 'Bolt shear strength', Ps, 'N'),
      ],
      plots: [{ type: 'bar', title: 'Bolt utilisation by load case', ylabel: 'Utilisation [-]', categories: rows.map((r) => r.name), series: [{ name: 'Utilisation', y: rows.map((r) => r.util) }] }],
      tables: [{ title: 'Load cases (per most-loaded bolt, fitting factor included)', columns: ['Case', 'Tension [N]', 'Shear [N]', 'Rt', 'Rs', 'Utilisation'], rows: rows.map((r) => [r.name, r.T, r.S, r.Rt, r.Rs, r.util]) }],
      warnings,
      models: ['Rigid item on a rectangular four-bolt pattern: overturning couple shared by two bolts, shear shared equally', 'Bolt interaction Rs³ + Rt² = 1'],
      assumptions: ['Tensile stress area taken as 75% of the shank area; shear strength 0.6 × tensile strength', 'Rigid floor and item; no prying, preload or load redistribution after yielding', 'Each direction applied separately plus one combined forward–side resultant'],
    };
  },
  verify() {
    const b = { mass: 100, h_cg: 0, lx: 0.5, ly: 0.5, d_bolt: 0.01, Su_bolt: 800e6, fitting: 1, g_fwd: 9, g_up: 3, g_side: 0, g_down: 0, g_aft: 0 }, r = N.kv(emergency.run(b)), r2 = N.kv(emergency.run({ ...b, h_cg: 0.5, g_up: 0 }));
    return [
      N.check('Upward case: bolt tension = n·W/4', r.bolt_tension_N, (3 * 100 * G0) / 4, 1e-12, 'Statics'),
      N.check('Forward case: bolt shear = n·W/4', r.bolt_shear_N, (9 * 100 * G0) / 4, 1e-12, 'Statics'),
      N.check('Overturning: tension = F·h/(2·lx)', r2.bolt_tension_N, (9 * 100 * G0 * 0.5) / (2 * 0.5), 1e-12, 'Moment equilibrium about the rear bolt line'),
      N.check('Pure tension utilisation = T/Pt', N.kv(emergency.run({ ...b, g_fwd: 0 })).bolt_util, ((3 * 100 * G0) / 4) / (800e6 * 0.75 * Math.PI * 0.01 ** 2 / 4), 1e-9, 'Interaction curve end point'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.bolt_util > 1) out.push({ severity: 'critical', title: 'Attachment fails under emergency-landing loads', detail: `Utilisation ${o.bolt_util.toFixed(2)}.`, action: `Use bolts of about ${(o.d_bolt_req_m * 1e3 * 1.1).toFixed(1)} mm, widen the bolt pattern or lower the centre of gravity; check the floor fitting for ${(o.bolt_tension_N / 1e3).toFixed(1)} kN pull-out.`, basis: 'Ultimate inertia load factors with fitting factor' });
    else out.push({ severity: 'info', title: 'Attachment retains the item', detail: `Utilisation ${o.bolt_util.toFixed(2)}, margin ${o.mos_attach.toFixed(2)}.`, action: 'Carry the reactions into the floor structure check. For battery packs and tanks, add crush and penetration protection: retention alone does not prevent fire or leakage.', basis: 'Static emergency-landing conditions' });
    return out;
  },
};

// ---- explicit 3-D finite-element analyses (kernel in ../core/solvers/explicitfe.js) ---------
/** Zero-phase two-pole Butterworth low-pass of a uniformly sampled signal (forward and backward pass). */
function lowpass(y, fs, fc) {
  if (!(fc > 0) || fc >= 0.45 * fs || y.length < 4) return y.slice();
  const c = 1 / Math.tan((Math.PI * fc) / fs), a0 = 1 / (1 + Math.SQRT2 * c + c * c), b1 = 2 * (1 - c * c) * a0, b2 = (1 - Math.SQRT2 * c + c * c) * a0;
  const pass = (x) => { const o = new Array(x.length); let x1 = x[0], x2 = x[0], y1 = x[0], y2 = x[0]; for (let k = 0; k < x.length; k++) { const v = a0 * (x[k] + 2 * x1 + x2) - b1 * y1 - b2 * y2; x2 = x1; x1 = x[k]; y2 = y1; y1 = v; o[k] = v; } return o; };
  return pass(pass(y).reverse()).reverse();
}
/** Uniform-time recorder of scalar channels during an explicit run. */
function recorder(tEnd, n, chans) {
  const dt = tEnd / n, t = [], data = chans.map(() => []); let k = 0;
  return { dt, t, data, take(fe) { while (k <= n && fe.t >= k * dt - 1e-12 * tEnd) { t.push(k * dt); for (let j = 0; j < chans.length; j++) data[j].push(chans[j](fe)); k++; } } };
}
/** Central-difference derivative of a uniformly sampled history. */
const ddt = (v, dt) => v.map((_, k) => { const a = Math.max(k - 1, 0), b = Math.min(k + 1, v.length - 1); return b > a ? (v[b] - v[a]) / ((b - a) * dt) : 0; });
/** Area, second moment, plastic modulus and depth of a stack of rectangles [width, thickness] listed from the skin side inwards. */
function stackSection(rects) {
  let A = 0, Q = 0, y = 0; const ps = rects.map(([b, t]) => { const p = { b, t, y0: y }; y += t; A += b * t; Q += b * t * (p.y0 + t / 2); return p; });
  const yc = Q / A; let I = 0, acc = 0, yp = y / 2, Zp = 0;
  for (const p of ps) I += (p.b * p.t ** 3) / 12 + p.b * p.t * (p.y0 + p.t / 2 - yc) ** 2;
  for (const p of ps) { if (acc + p.b * p.t >= A / 2) { yp = p.y0 + (A / 2 - acc) / p.b; break; } acc += p.b * p.t; }
  for (const p of ps) { const lo = p.y0, hi = p.y0 + p.t; Zp += hi <= yp || lo >= yp ? p.b * p.t * Math.abs((lo + hi) / 2 - yp) : (p.b * ((yp - lo) ** 2 + (hi - yp) ** 2)) / 2; }
  return { A, I, Zp, yc, h: y };
}
/** Seat and occupant on a floor node: seat mass on an elastic–plastic energy absorber, upper body on the DRI spine spring–damper. */
function addOccupant(fe, floor, pos, up, o) {
  const seat = fe.node(pos[0] + 0.5 * o.hgt * up[0], pos[1] + 0.5 * o.hgt * up[1], pos[2] + 0.5 * o.hgt * up[2]), torso = fe.node(pos[0] + o.hgt * up[0], pos[1] + o.hgt * up[1], pos[2] + o.hgt * up[2]);
  fe.mass(seat, o.m_seat); fe.mass(torso, o.m_torso);
  const kS = o.m_torso * DRI_W ** 2;
  return { seat, torso, sSeat: fe.spring(floor, seat, { k: o.k_seat, Fyc: o.F_ea, Fyt: o.F_ea, dmax: o.s_ea + o.F_ea / o.k_seat, kb: 20 * o.k_seat, kLat: o.k_seat }, up), sSpine: fe.spring(seat, torso, { k: kS, c: 2 * DRI_Z * o.m_torso * DRI_W, kLat: 4 * kS }, up) };
}
const seatSplit = (m, occupied, g, s_ea) => ({ m_seat: (occupied ? 55 / 89 : 0.6) * m, m_torso: (occupied ? 34 / 89 : 0.4) * m, F_ea: g * G0 * m, s_ea, k_seat: Math.max((2.5e6 * m) / 89, 100) });
/** Structured block of hexahedra; map transforms the regular grid coordinates. */
function hexBlock(fe, nx, ny, nz, Lx, Ly, Lz, mt, map = (p) => p) {
  const id = (i, j, k) => (i * (ny + 1) + j) * (nz + 1) + k, nd = [], els = [];
  for (let i = 0; i <= nx; i++) for (let j = 0; j <= ny; j++) for (let k = 0; k <= nz; k++) nd.push(fe.node(...map([(i * Lx) / nx, (j * Ly) / ny, (k * Lz) / nz])));
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) for (let k = 0; k < nz; k++) els.push(fe.hex([id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k), id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)].map((q) => nd[q]), mt));
  return { id: (i, j, k) => nd[id(i, j, k)], els, nd };
}
const RES_NOTE = 'Default resolution is deliberately coarse so the analysis runs in about a second in the browser; raise the mesh and step settings in the Numerics group (and run the convergence study) before using the numbers.';

// ---- 8. fuselage barrel-section drop test (explicit FE) -------------------------------------
function barrelSim(i, ctx) {
  const m = mat(i.material), R = i.D / 2, n = Math.max(8, 2 * Math.round(i.nRing / 2)), hf = i.h_frame, tf = i.t_frame, bf = 0.3 * hf, capH = m.Su / m.Sy;
  const sf = stackSection([[Math.min(i.pitch, 30 * i.t_skin), i.t_skin], [bf, tf], [tf, Math.max(hf - 2 * tf, tf)], [bf, tf]]), hb = 1.5 * hf, tb = 1.25 * tf, sb = stackSection([[0.4 * hb, tb], [tb, Math.max(hb - 2 * tb, tb)], [0.4 * hb, tb]]);
  const fe = createFE({ gravity: [0, 0, -G0] }), ro = N.rad(i.roll_deg), cr = Math.cos(ro), sr = Math.sin(ro), P = (x, z) => [x * cr + z * sr, 0, -x * sr + z * cr], pts = [];
  const node = (x, z) => { const p = P(x, z); pts.push(p); return pts.length - 1; };
  // ring meshed in arcs between the keel, the strut feet, the floor attachments and the crown, so the geometry does not depend on the mesh
  const thf = Math.acos(1 - 2 * i.floor_frac), ths = i.strut_frac * thf, xf = R * Math.sin(thf), zf = -R * Math.cos(thf), H = R + zf, xk = i.struts ? R * Math.sin(ths) : 0.5 * xf, Lr = (2 * Math.PI * R) / n;
  const brk = i.struts ? [0, ths, thf, Math.PI] : [0, thf, Math.PI], half = [0]; for (let s = 0; s + 1 < brk.length; s++) { const ns = Math.max(1, Math.round((n / 2) * (brk[s + 1] - brk[s]) / Math.PI)); for (let q = 1; q <= ns; q++) half.push(brk[s] + ((brk[s + 1] - brk[s]) * q) / ns); }
  const angs = [...half, ...half.slice(1, -1).reverse().map((a) => -a)], ring = angs.map((a) => node(R * Math.sin(a), -R * Math.cos(a))), nr = ring.length, at = (a) => ring[N.argmin(angs.map((q) => Math.abs(q - a)))];
  // floor beam through the key points: ring attachment, seat tracks / strut tops, centre
  const keys = [-xf, -xk, 0, xk, xf], floor = [at(-thf)], keyNode = [at(-thf)];
  for (let s = 0; s < 4; s++) { const ns = Math.max(1, Math.round((keys[s + 1] - keys[s]) / Lr)); for (let q = 1; q <= ns; q++) { const last = s === 3 && q === ns, nd = last ? at(thf) : node(keys[s] + ((keys[s + 1] - keys[s]) * q) / ns, zf); floor.push(nd); if (q === ns) keyNode.push(nd); } }
  const zmin = N.amin(pts.map((p) => p[2])), gap = 1e-4 * i.D; for (const p of pts) { p[2] += gap - zmin; fe.node(p[0], p[1], p[2]); }
  const base = { E: m.E, G: m.G, rho: m.rho, hard: (capH - 1) / 0.2, cap: capH, thetaSoft: i.theta_soft, residual: 0.25 }, beams = [];
  const secR = { ...base, A: sf.A, Iy: sf.I, Iz: sf.I, J: sf.I, Np: m.Sy * sf.A, Mpy: i.hinge_eff * m.Sy * sf.Zp, Mpz: i.hinge_eff * m.Sy * sf.Zp }, secB = { ...base, A: sb.A, Iy: sb.I, Iz: sb.I, J: sb.I, Np: m.Sy * sb.A, Mpy: i.hinge_eff * m.Sy * sb.Zp, Mpz: i.hinge_eff * m.Sy * sb.Zp };
  for (let k = 0; k < nr; k++) beams.push({ e: fe.beam(ring[k], ring[(k + 1) % nr], secR, [0, 1, 0]), w: hf, part: 'Frame' });
  for (let k = 0; k + 1 < floor.length; k++) beams.push({ e: fe.beam(floor[k], floor[k + 1], secB, [0, 1, 0]), w: hb, part: 'Floor beam' });
  if (i.struts) for (const [a, b] of [[keyNode[1], at(-ths)], [keyNode[3], at(ths)]]) {
    const Ls = Math.hypot(pts[a][0] - pts[b][0], pts[a][2] - pts[b][2]), As = 0.3 * sf.A, rs = 0.3 * hf, Is = As * rs * rs;
    beams.push({ e: fe.beam(a, b, { ...base, A: As, Iy: Is, Iz: Is, J: Is, Np: i.hinge_eff * Math.min(m.Sy * As, (Math.PI ** 2 * m.E * Is) / (Ls * Ls)), Mpy: 0.6 * m.Sy * As * rs, Mpz: 0.6 * m.Sy * As * rs }, [0, 1, 0]), w: 0.6 * hf, part: 'Strut' });
  }
  for (const k of ring) fe.mass(k, i.m_shell / nr);
  const fl = floor.slice(1, -1); for (const k of fl) fe.mass(k, i.m_floor / fl.length);
  const up = P(0, 1), occ = [1, 3].map((q) => addOccupant(fe, keyNode[q], pts[keyNode[q]], up, { ...seatSplit(i.m_occ / 2, i.occupied, i.seat_g, i.s_ea), hgt: 0.12 * i.D }));
  const pl = fe.plane({ mu: i.mu }); for (const k of [...ring, ...fl]) pl.contact(k);
  for (let k = 0; k < pts.length + 4; k++) fe.fix(k, [0, 1, 0, 1, 0, 1]);
  fe.velocity(-1, [0, 0, -i.v0]); fe.init();
  const c = keyNode[2], z0 = fe.x[3 * c + 2], avg = (a, b, arr) => 0.5 * (arr[3 * a + 2] + arr[3 * b + 2]), wpB = () => { let s = 0; for (const b of beams) s += fe.sec[b.e].wp || 0; return s; }, E = fe.energy;
  const rec = recorder(i.t_end, 600, [(f) => avg(keyNode[1], keyNode[3], f.v), (f) => avg(occ[0].seat, occ[1].seat, f.v), (f) => avg(occ[0].torso, occ[1].torso, f.v), (f) => z0 - f.x[3 * c + 2], () => E.kinetic, () => E.internal, () => wpB() + E.crush, () => E.contact, () => E.external, () => E.total, () => pl.force]);
  const snaps = [], nSnap = 30, take = () => snaps.push({ t: fe.t, x: Array.from(fe.x), k: Array.from(fe.beamKappa) }); let tReb = NaN, crush = 0, errMax = 0; rec.take(fe); take();
  fe.run(i.t_end, () => {
    rec.take(fe); if (fe.t >= (snaps.length * i.t_end) / nSnap) take();
    const d = z0 - fe.x[3 * c + 2]; if (d > crush) crush = d; if (E.error > errMax) errMax = E.error;
    if (!(tReb > 0) && d > 2 * gap && fe.v[3 * c + 2] >= 0) tReb = fe.t;
    if (fe.nStep % 500 === 0) ctx?.progress?.(fe.t / i.t_end, 'Explicit integration');
    return tReb > 0 && fe.t > 1.5 * tReb && occ.every((o) => fe.springs[o.sSpine].d < 0.95 * fe.springs[o.sSpine].dMaxSeen); // floor has rebounded and the spinal response is past its peak
  });
  take();
  return { fe, rec, snaps, beams, occ, pl, tReb, crush: Math.max(0, crush - gap), errMax, H, n, sf, sb, m, keyNode, ring, floor, nNodes: pts.length + 4, wpB: wpB() };
}
/** Triangulated ribbons of the beam mesh for a 'tri' plot: each element is a strip of its section depth, coloured by the hinge strain estimate. */
function ribbons(r, snap, wMin) {
  const nodes = [], tris = [], values = [], fe = r.fe;
  for (const b of r.beams) {
    if (!fe.beamAlive[b.e]) continue;
    const n1 = fe.beamN[2 * b.e], n2 = fe.beamN[2 * b.e + 1], x1 = snap.x[3 * n1], z1 = snap.x[3 * n1 + 2], x2 = snap.x[3 * n2], z2 = snap.x[3 * n2 + 2], L = Math.hypot(x2 - x1, z2 - z1) || 1, w = Math.max(b.w, wMin) / 2, nx = (-(z2 - z1) / L) * w, nz = ((x2 - x1) / L) * w, q = nodes.length;
    nodes.push([x1 + nx, z1 + nz], [x1 - nx, z1 - nz], [x2 + nx, z2 + nz], [x2 - nx, z2 - nz]); tris.push([q, q + 1, q + 2], [q + 1, q + 3, q + 2]);
    const e1 = 0.5 * snap.k[2 * b.e], e2 = 0.5 * snap.k[2 * b.e + 1]; values.push(e1, e1, e2, e2);
  }
  return { nodes, tris, values };
}
const BARREL_BASE = { D: 3.95, floor_frac: 0.35, pitch: 0.5, h_frame: 0.1, t_frame: 0.002, t_skin: 0.0016, struts: true, strut_frac: 0.55, material: 'Al 2024-T3', hinge_eff: 0.7, theta_soft: 0.6, m_occ: 270, m_floor: 60, m_shell: 90, occupied: true, seat_g: 12, s_ea: 0.15, v0: 9.1, roll_deg: 0, mu: 0.3, filter_Hz: 100, t_end: 0.4, nRing: 24 };
const barrel3d = {
  id: 'barrel3d', title: 'Fuselage barrel-section drop test (explicit 3-D finite elements)', fidelity: 'numerical',
  summary: 'One frame bay of the fuselage barrel — ring frame with its effective skin, floor beam, cargo-bay struts and seated occupants — is dropped onto rigid ground in the explicit finite-element kernel. Plastic hinges form and the lower lobe crushes; the run gives deformed shapes, floor and occupant deceleration histories, crush distance, absorbed energy, the energy balance and the Dynamic Response Index.',
  equations: ['Conservation of linear momentum', 'Conservation of angular momentum', 'Conservation of energy', 'Transient structural dynamics equations', 'Plasticity constitutive equations', 'Impact contact equations', 'Crushing energy absorption equations', 'Occupant dynamics equations'],
  applicable: (c) => (c.fuselage.dia_m > 0 ? true : 'No fuselage diameter is defined for this vehicle.'),
  inputs: [
    { key: 'D', label: 'Fuselage diameter', unit: 'm', default: 3.95, min: 0.05, max: 10, group: 'Geometry' },
    { key: 'floor_frac', label: 'Cabin floor height above the belly / diameter', unit: '-', default: 0.35, min: 0.1, max: 0.5, group: 'Geometry', help: 'The crushable lower lobe is the structure below the floor' },
    { key: 'pitch', label: 'Frame pitch (length of the bay modelled)', unit: 'm', default: 0.5, min: 0.005, max: 2, group: 'Geometry' },
    { key: 'h_frame', label: 'Frame depth', unit: 'm', default: 0.1, min: 0.002, max: 0.5, group: 'Geometry' },
    { key: 't_frame', label: 'Frame web and flange thickness', unit: 'm', default: 0.002, min: 1e-4, max: 0.02, group: 'Geometry', help: 'Channel section with flanges of 0.3 × depth' },
    { key: 't_skin', label: 'Skin thickness', unit: 'm', default: 0.0016, min: 1e-4, max: 0.02, group: 'Geometry', help: 'An effective width of 30 thicknesses works with the frame' },
    { key: 'struts', label: 'Cargo-bay struts under the floor beam', type: 'bool', default: true, group: 'Geometry' },
    { key: 'strut_frac', label: 'Strut foot position along the lower arc', unit: '-', default: 0.55, min: 0.1, max: 0.9, group: 'Geometry', help: 'Fraction of the arc from the keel to the floor attachment' },
    MAT,
    { key: 'hinge_eff', label: 'Plastic-moment efficiency of the thin-walled sections', unit: '-', default: 0.7, min: 0.1, max: 1.2, group: 'Material', help: 'Local buckling and crippling stop thin frames reaching the full plastic moment; calibrate against component bending tests' },
    { key: 'theta_soft', label: 'Hinge rotation at frame fracture', unit: 'rad', default: 0.6, min: 0.05, max: 3, group: 'Material', help: 'Beyond this plastic rotation a hinge keeps 25% of its moment (cracked frame held by the skin)' },
    { key: 'm_occ', label: 'Occupants and seats carried by this bay', unit: 'kg', default: 270, min: 0.001, group: 'Masses' },
    { key: 'm_floor', label: 'Floor structure, systems and cargo on the floor beam', unit: 'kg', default: 60, min: 0.001, group: 'Masses' },
    { key: 'm_shell', label: 'Non-modelled mass on the shell (skin, stringers, lining, bins)', unit: 'kg', default: 90, min: 0.001, group: 'Masses' },
    { key: 'occupied', label: 'Occupied (evaluate injury metrics)', type: 'bool', default: true, group: 'Masses' },
    { key: 'seat_g', label: 'Seat energy-absorber limit load factor', unit: 'g', default: 12, min: 1, max: 200, group: 'Seat', help: 'On the occupant-plus-seat mass; a rigid seat has a very high value' },
    { key: 's_ea', label: 'Available seat stroke', unit: 'm', default: 0.15, min: 0, max: 1, group: 'Seat' },
    { key: 'v0', label: 'Vertical impact velocity', unit: 'm/s', default: 9.14, min: 0.1, max: 40, group: 'Impact', help: '9.14 m/s (30 ft/s) is the rotorcraft seat-test velocity change and a customary section drop-test velocity; the transport seat test uses 10.67 m/s (35 ft/s)' },
    { key: 'roll_deg', label: 'Roll angle at impact', unit: 'deg', default: 0, min: -45, max: 45, group: 'Impact' },
    { key: 'mu', label: 'Ground friction coefficient', unit: '-', default: 0.3, min: 0, max: 1.5, group: 'Impact' },
    { key: 'filter_Hz', label: 'Low-pass filter for deceleration histories', unit: 'Hz', default: 100, min: 1, max: 1e5, group: 'Numerics', help: 'Explicit nodal accelerations are noisy; peak g depends on this cut-off, as it does in a test' },
    { key: 't_end', label: 'Maximum simulated time', unit: 's', default: 0.4, min: 1e-4, max: 5, group: 'Numerics', help: 'The run stops earlier once the floor has rebounded and the spinal response has passed its peak' },
    { key: 'nRing', label: 'Beam elements around the frame', unit: '', default: 24, min: 12, max: 240, step: 2, discrete: true, group: 'Numerics' },
  ],
  defaults(c) {
    const D = c.fuselage.dia_m, occ = occupants(c), has = occ > 0, rot = c.meta.type === 'helicopter' || c.meta.type === 'evtol', pitch = N.clamp(0.13 * D, 0.01, 0.53), cabin = Math.max(D, (rot ? 0.3 : occ > 19 ? 0.72 : 0.35) * c.fuselage.len_m);
    const shell = Math.max(1e-3, (0.22 * (c.mass.oew_kg || 0.5 * c.mass.mtow_kg) * pitch) / Math.max(c.fuselage.len_m, D)), v0 = has ? 9.14 : Math.max(2, 2 * c.gear.sink_ms), ff = rot ? 0.2 : 0.35; // 30 ft/s section drop for occupied cabins
    return { D, floor_frac: ff, pitch, h_frame: Math.max(0.004, 0.025 * D), t_frame: Math.max(4e-4, 5e-4 * D), t_skin: N.clamp(4e-4 * D, 4e-4, 3e-3), material: c.struct.material, m_occ: has ? (occ * 89 * pitch) / cabin : Math.max(1e-3, ((c.mass.payload_kg || 0.2 * c.mass.mtow_kg) * pitch) / cabin), m_floor: 0.4 * shell, m_shell: 0.6 * shell, occupied: has, seat_g: has ? 12 : 30, s_ea: !has ? 0.02 * D : seatStroke(v0, 0.4 * ff * D, 12, rot ? 0.3 : 0.15), v0, filter_Hz: (15 * v0) / (ff * D), t_end: (2.2 * ff * D) / v0 + 0.08 };
  },
  run(i, ctx) {
    const r = barrelSim(i, ctx), fe = r.fe, E = fe.energy, rec = r.rec, T = rec.t, fs = 1 / rec.dt, g = (v) => lowpass(ddt(v, rec.dt), fs, i.filter_Hz).map((a) => a / G0 + 1), aF = g(rec.data[0]), aS = g(rec.data[1]), aT = g(rec.data[2]), warnings = [RES_NOTE];
    const pk = N.amax(aF), pkS = N.amax(aS), dri = (DRI_W ** 2 * Math.max(...r.occ.map((o) => fe.springs[o.sSpine].dMaxSeen))) / G0, stroke = Math.max(...r.occ.map((o) => fe.springs[o.sSeat].dpMax)), Ein = E.initial + E.external, absorbed = r.wpB + E.crush + E.friction;
    const kap = Array.from(fe.beamKappa), hinges = r.beams.filter((b) => Math.max(kap[2 * b.e], kap[2 * b.e + 1]) > 0.02), kMax = N.amax(kap), broken = r.beams.filter((b) => Math.max(kap[2 * b.e], kap[2 * b.e + 1]) >= i.theta_soft).length;
    if (!(r.tReb > 0)) warnings.push('The floor is still moving downwards at the end of the simulated time: extend it. Crush and absorbed energy are lower bounds.');
    if (r.crush > 0.9 * r.H) warnings.push('The lower lobe is crushed almost flat: the cabin floor reaches the ground and the loads are then governed by the floor structure, which this model treats as a single beam.');
    if (broken) warnings.push(`${broken} frame or floor element(s) rotated beyond the fracture rotation and carry only the residual moment.`);
    if (!i.occupied) warnings.push('No occupants in this case: the "upper body" is the payload on a 52.9 rad/s mount and injury metrics are not meaningful.');
    if (r.errMax > 0.05) warnings.push(`The energy balance error reached ${(100 * r.errMax).toFixed(1)}%: treat the result with caution and refine the mesh.`);
    if (i.h_frame > (0.5 * Math.PI * i.D) / r.n) warnings.push('The frame elements are short compared with the frame depth: Euler–Bernoulli beam theory neglects the shear deformation that matters there.');
    const pick = (t) => r.snaps[N.argmin(r.snaps.map((s) => Math.abs(s.t - t)))], tc = T[N.argmax(rec.data[3])], shots = [pick(0.4 * tc), pick(tc), r.snaps[r.snaps.length - 1]], wMin = 0.012 * i.D, xs = shots[0].x.filter((_, k) => k % 3 === 0), xr = [N.amin(xs) - 0.1 * i.D, N.amax(xs) + 0.1 * i.D];
    const ms = T.map((t) => t * 1e3), th = (a) => thin(a, 300), driSt = !i.occupied ? undefined : dri <= 18 ? 'ok' : dri <= 22.8 ? 'warn' : 'bad';
    return {
      kpis: [
        kpi('barrel_peak_g', 'Floor (seat-track) peak deceleration', pk, 'g', undefined, `Low-pass filtered at ${i.filter_Hz.toFixed(0)} Hz`),
        kpi('barrel_crush_m', 'Crush of the lower lobe (floor centre)', r.crush, 'm', r.crush > 0.9 * r.H ? 'bad' : r.crush > 0.7 * r.H ? 'warn' : 'ok', `Floor is ${r.H.toFixed(3)} m above the belly`),
        kpi('barrel_absorbed_J', 'Energy absorbed (plastic hinges, seats, friction)', absorbed, 'J'),
        kpi('barrel_absorbed_frac', 'Share of the impact energy absorbed', absorbed / Ein, '-'),
        kpi('barrel_DRI', 'Dynamic Response Index', dri, '-', driSt, 'About 18 corresponds to roughly 5% and 22.8 to 50% probability of spinal injury'),
        kpi('barrel_pelvis_g', 'Seat-pan peak deceleration', pkS, 'g'), kpi('barrel_seat_stroke_m', i.occupied ? 'Seat stroke used' : 'Payload-mount stroke used', stroke, 'm', stroke > i.s_ea && i.s_ea > 0 ? (i.occupied ? 'bad' : 'warn') : 'ok', `Available ${i.s_ea.toFixed(3)} m`),
        kpi('barrel_hinges', 'Elements with plastic hinges', hinges.length, '-'), kpi('barrel_hinge_rot_rad', 'Largest plastic hinge rotation', kMax, 'rad', kMax >= i.theta_soft ? 'warn' : 'ok', `Fracture rotation ${i.theta_soft} rad`),
        kpi('barrel_ground_force_N', 'Peak ground reaction', N.amax(rec.data[10]), 'N'), kpi('barrel_t_stop_s', 'Time to maximum crush', tc, 's'),
        kpi('barrel_energy_err', 'Energy balance error', r.errMax, '-', r.errMax < 0.03 ? 'ok' : 'warn', 'Kinetic + internal + contact − gravity work against the initial energy'),
        kpi('barrel_mass_kg', 'Mass of the modelled bay', fe.totalMass, 'kg'), kpi('barrel_elements', 'Beam elements', fe.nb, '-'), kpi('barrel_steps', 'Time steps', fe.nStep, '-'), kpi('barrel_dt_s', 'Stable time step', fe.dt, 's'),
      ],
      plots: [
        ...shots.map((s, k) => ({ type: 'tri', title: `Deformed section at ${(s.t * 1e3).toFixed(1)} ms${k === 1 ? ' (maximum crush)' : k === 2 ? ' (end of run)' : ''}`, xlabel: 'Lateral position [m]', ylabel: 'Height [m]', zlabel: 'Plastic strain at the hinge (outer-fibre estimate) [-]', ...ribbons(r, s, wMin), equalAspect: true, edges: false, overlay: [{ name: 'Ground', x: xr, y: [0, 0] }, { name: 'Seats and occupants', x: r.occ.flatMap((o) => [s.x[3 * o.seat], s.x[3 * o.torso], NaN]), y: r.occ.flatMap((o) => [s.x[3 * o.seat + 2], s.x[3 * o.torso + 2], NaN]) }] })),
        { type: 'line', title: 'Deceleration histories', xlabel: 'Time [ms]', ylabel: 'Deceleration [g]', series: [{ name: 'Floor at the seat tracks', x: th(ms), y: th(aF) }, { name: 'Seat pan', x: th(ms), y: th(aS) }, { name: 'Upper body', x: th(ms), y: th(aT) }] },
        { type: 'line', title: 'Energy balance', xlabel: 'Time [ms]', ylabel: 'Energy [kJ]', series: [['Kinetic', 4], ['Internal (elastic + plastic)', 5], ['Plastic work', 6], ['Contact and friction', 7], ['Gravity work', 8], ['Total (should stay constant)', 9]].map(([name, k]) => ({ name, x: th(ms), y: th(rec.data[k]).map((e) => e / 1e3), ...(k === 9 ? { style: 'dash' } : {}) })) },
        { type: 'line', title: 'Crush of the lower lobe', xlabel: 'Time [ms]', ylabel: 'Floor-centre displacement [mm]', series: [{ name: 'Crush', x: th(ms), y: th(rec.data[3]).map((d) => d * 1e3) }], annotations: [{ y: r.H * 1e3, label: 'Floor reaches the ground' }] },
      ],
      tables: [{ title: 'Plastic hinges', columns: ['Member', 'Element', 'Plastic rotation [rad]', 'First yield [ms]'], rows: hinges.sort((a, b) => Math.max(kap[2 * b.e], kap[2 * b.e + 1]) - Math.max(kap[2 * a.e], kap[2 * a.e + 1])).slice(0, 12).map((b) => [b.part, b.e, Math.max(kap[2 * b.e], kap[2 * b.e + 1]), fe.beamYieldT[b.e] * 1e3]) }],
      outputs: { barrel_energy_in_J: Ein },
      warnings,
      models: [`Explicit central-difference finite elements: ${fe.nb} co-rotational beam elements, ${r.nNodes} nodes, ${fe.nStep} steps at the Courant limit`, 'Elastic–plastic stress resultants with plastic hinges, axial-force interaction, hardening to the ultimate strength and post-fracture residual moment', 'Penalty contact with Coulomb friction on rigid ground', 'Seat energy absorber and Dynamic Response Index spinal model (52.9 rad/s, ζ = 0.224) carried on the floor beam'],
      assumptions: ['One frame bay with planar deformation (no fore–aft motion): every frame of the barrel behaves alike and the skin acts only as an effective flange of the frame', 'Thin-walled sections are represented by their stress resultants; local buckling enters through the plastic-moment efficiency and the fracture rotation, both test-calibrated inputs', 'No contact between structural members (floor against frames, struts against skin) and no cargo, luggage or fuel below the floor', 'Rigid ground; the deceleration peaks depend on the stated filter cut-off', DATA_NOTE],
    };
  },
  convergence: { param: 'nRing', label: 'Beam elements around the frame', levels: [16, 24, 32, 48], metric: 'barrel_crush_m' },
  calibration: { params: [{ key: 'hinge_eff', min: 0.2, max: 1.2 }, { key: 'theta_soft', min: 0.05, max: 2 }], sweep: 'v0', target: 'barrel_crush_m', note: 'Section drop-test crush distances (or floor decelerations) at one or more impact velocities calibrate the hinge efficiency and fracture rotation.' },
  verify: () => kernelChecks(),
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.barrel_crush_m > 0.9 * i.floor_frac * i.D) out.push({ severity: 'critical', title: 'Lower lobe collapses completely', detail: `Crush ${o.barrel_crush_m.toFixed(2)} m of about ${(i.floor_frac * i.D).toFixed(2)} m available at ${i.v0} m/s.`, action: 'Add energy-absorbing struts or a crushable sub-floor (corrugated webs, foam or composite absorbers sized in the crush-tube analysis), or thicken the lower frames between the keel and the strut feet.', basis: 'Energy absorbed by plastic hinges against the impact energy' });
    if (i.occupied && o.barrel_DRI > 18) out.push({ severity: o.barrel_DRI > 22.8 ? 'critical' : 'warn', title: 'High spinal injury risk', detail: `DRI ${o.barrel_DRI.toFixed(1)}, floor peak ${o.barrel_peak_g.toFixed(0)} g.`, action: o.barrel_seat_stroke_m >= i.s_ea ? 'Provide more seat stroke or soften the structure below the floor so that it crushes at a lower, longer plateau.' : 'Lower the seat energy-absorber limit load, or reduce the stiffness of the strut load path that transmits the ground reaction straight to the seat tracks.', basis: 'DRI ≈ 18 for about 5% spinal injury probability' });
    else if (i.occupied) out.push({ severity: 'info', title: 'Occupant loads are within tolerance for this drop', detail: `DRI ${o.barrel_DRI.toFixed(1)} with ${(100 * o.barrel_crush_m).toFixed(0)} cm of structural crush.`, action: 'Check the off-nominal cases next: rolled impact, heavier occupants and a loaded cargo bay, which stiffens the lower lobe and raises the floor pulse.', basis: 'Crashworthiness system trade' });
    out.push({ severity: 'advise', title: 'Resolution and calibration', detail: `${o.barrel_elements} elements, energy balance error ${(100 * o.barrel_energy_err).toFixed(2)}%.`, action: 'Run the convergence study, then calibrate the hinge efficiency and fracture rotation against a component or section test before relying on the peak deceleration.', basis: 'Verification and validation practice for explicit crash models' });
    return out;
  },
};

// ---- 9. full-aircraft hybrid crash model (lumped masses, nonlinear beams, crush springs) ----
const CONFIGS = ['aeroplane', 'helicopter', 'multirotor'];
function aircraftSim(i, ctx) {
  const m = mat(i.material), L = i.len_m, D = i.dia_m, n = Math.max(3, Math.min(Math.round(i.nStations), 1 + Math.floor(L / (0.75 * D)))), M = i.mass_kg, capH = m.Su / m.Sy, cfg = i.config, heli = cfg === 'helicopter', multi = cfg === 'multirotor', wing = !multi && i.span_m > 0;
  const th = N.rad(i.pitch_deg), ph = N.rad(i.roll_deg), ct = Math.cos(th), st = Math.sin(th), cp = Math.cos(ph), sp = Math.sin(ph);
  // body axes: x forward, y left, z up. Roll (right side down) about x, then pitch (nose up) about y.
  const Rb = (p) => { const y = p[1] * cp - p[2] * sp, z = p[1] * sp + p[2] * cp; return [p[0] * ct - z * st, y, p[0] * st + z * ct]; };
  const P = [], tgt = [], node = (p, mass) => { P.push(p); tgt.push(mass); return P.length - 1; }, beams = [], fe = createFE({ gravity: [0, 0, -G0] });
  // ---- masses
  const mEngW = wing ? i.n_eng_wing * i.m_engine : 0, mArms = multi ? i.n_eng_wing * i.m_engine : 0, mWing = wing ? i.wing_mass_frac * M + i.m_fuel_wing : 0, mTail = multi ? 0 : 0.02 * M, mOh = heli ? i.m_overhead : 0;
  const mFus = Math.max(0.15 * M, M - mEngW - mArms - mWing - mTail - mOh - i.m_nose), wts = N.range(n, (k) => 0.4 + Math.sin((Math.PI * k) / (n - 1))), wS = N.sum(wts);
  const iw = N.clamp(Math.round(0.45 * (n - 1)), 0, n - 1), stn = N.range(n, (k) => node([L / 2 - (L * k) / (n - 1), 0, 0], (mFus * wts[k]) / wS + (k === 0 ? i.m_nose : 0) + (k === n - 1 ? 0.5 * mTail : 0) + (k === iw ? 0.3 * mWing : 0)));
  const sup = stn.map((s) => tgt[s]), X = (k) => P[stn[k]][0], sy = m.Sy;
  // ---- fuselage beams
  const Af = Math.PI * D * i.t_fus, If = Math.PI * (D / 2) ** 3 * i.t_fus, MpF = i.buckle_eff * sy * D * D * i.t_fus, base = { E: m.E, G: m.G, rho: m.rho, hard: (capH - 1) / 0.1, cap: capH };
  const secF = { ...base, A: Af, Iy: If, Iz: If, J: 2 * If, Np: sy * Af, Mpy: MpF, Mpz: MpF, thetaSoft: i.theta_break, residual: 0.2 }, add = (a, b, sec, ref, name, cabin) => beams.push({ a, b, sec, ref, name, cabin });
  for (let k = 0; k + 1 < n; k++) add(stn[k], stn[k + 1], secF, [0, 1, 0], `Fuselage, station ${k + 1}–${k + 2}`, true);
  // ---- wing or rotor arms
  const lat = []; // nodes of the lateral members for plotting, per side
  if (wing) {
    const b2 = i.span_m / 2, tn = Math.tan(N.rad(i.sweep_deg)), zw = (i.wing_pos === 'low' ? -0.3 : 0.4) * D, dih = i.wing_pos === 'low' ? Math.tan(N.rad(5)) : 0, Iw = i.EI_wing / m.E;
    const secW = (f) => ({ ...base, A: f * i.A_wing, Iy: f * f * Iw, Iz: 4 * f * f * Iw, J: f * f * Iw, Np: sy * f * i.A_wing, Mpy: f * Math.sqrt(f) * i.Mp_wing, Mpz: 3 * f * Math.sqrt(f) * i.Mp_wing, thetaSoft: 0.15, residual: 0.1 });
    for (const s of [1, -1]) {
      const eng = i.n_eng_wing >= 2 ? i.m_engine : 0, eng2 = i.n_eng_wing >= 4 ? i.m_engine : 0;
      const mid = node([X(iw) - 0.35 * b2 * tn, s * 0.35 * b2, zw + 0.35 * b2 * dih], 0.2 * mWing + eng), out = node([X(iw) - 0.85 * b2 * tn, s * 0.85 * b2, zw + 0.85 * b2 * dih], 0.15 * mWing + eng2);
      add(stn[iw], mid, secW(1), [1, 0, 0], `${s > 0 ? 'Left' : 'Right'} wing, inboard`, false); add(mid, out, secW(0.4), [1, 0, 0], `${s > 0 ? 'Left' : 'Right'} wing, outboard`, false);
      lat.push({ s, mid, out, eng }); sup[iw] += 0.35 * mWing + eng + eng2;
    }
  }
  if (multi) {
    const a = i.arm_m, da = 0.1 * a, ta = 0.06 * da, Aa = Math.PI * da * ta, Ia = (Math.PI * da ** 3 * ta) / 8, secA = { ...base, A: Aa, Iy: Ia, Iz: Ia, J: 2 * Ia, Np: sy * Aa, Mpy: sy * da * da * ta, Mpz: sy * da * da * ta, thetaSoft: 0.3, residual: 0.1 }, nr = Math.max(2, Math.round(i.n_eng_wing));
    for (let q = 0; q < nr; q++) { const az = (2 * Math.PI * (q + 0.5)) / nr, tip = node([X(iw) + a * Math.cos(az), a * Math.sin(az), 0.1 * D], i.m_engine); add(stn[iw], tip, secA, [0, 0, 1], `Rotor arm ${q + 1}`, false); lat.push({ s: Math.sin(az) >= 0 ? 1 : -1, tip }); sup[iw] += i.m_engine; }
  }
  // ---- tail fin and overhead rotor / transmission mass
  let fin = -1, oh = -1; const ohLegs = [];
  if (!multi) { fin = node([X(n - 1), 0, 0.5 * D + i.fin_h], 0.5 * mTail); add(stn[n - 1], fin, { ...secF, A: 0.1 * Af, Iy: 0.02 * If, Iz: 0.02 * If, J: 0.04 * If, Np: 0.1 * sy * Af, Mpy: 0.03 * MpF, Mpz: 0.03 * MpF }, [1, 0, 0], 'Fin', false); sup[n - 1] += 0.5 * mTail; }
  if (heli && mOh > 0) {
    const j = N.clamp(Math.round(0.4 * (n - 1)), 1, n - 2); oh = node([X(j), 0, 0.85 * D], mOh); sup[j] += mOh;
    const Am = (i.mount_g * mOh * G0) / (1.5 * sy), rm = 0.04 * D, secM = { ...base, A: Am, Iy: Am * rm * rm, Iz: Am * rm * rm, J: 2 * Am * rm * rm, Np: sy * Am, Mpy: sy * Am * rm, Mpz: sy * Am * rm, epsFail: 0.1, thetaSoft: 0.3, residual: 0.2 };
    for (const k of [j - 1, j, j + 1]) { add(stn[k], oh, secM, [0, 1, 0], `Transmission mount to station ${k + 1}`, false); ohLegs.push(k); }
  }
  // ---- attitude, then shift so that the lowest contact point just touches
  const shape = (k) => (k === 0 ? 0.6 : k === n - 1 ? 0.5 : 1), cps = [], hw = heli || multi;
  const dc = (k) => i.d_crush * shape(k), cArea = ((L / (n - 1)) * D) / 2;
  for (let k = 0; k < n; k++) for (const s of [1, -1]) { const F = (i.crush_g * G0 * sup[k]) / 2, kk = F / (0.05 * dc(k)); cps.push({ node: stn[k], off: [0, s * 0.3 * D, -0.4 * D * shape(k)], law: i.surface === 'water' ? null : { k: kk, Fy: F, dmax: dc(k), kb: 20 * kk }, area: cArea, name: `Belly at station ${k + 1}`, kind: 'belly', stn: k }); }
  if (i.gear === 'extended') {
    const wheels = i.gear_type !== 'skid', aux = N.clamp(i.gear_type === 'tailwheel' ? n - 2 : 1, 0, n - 1), auxI = aux === iw ? (i.gear_type === 'tailwheel' ? n - 1 : 0) : aux, gl = { failAtMax: true, dmax: 1.2 * i.gear_stroke };
    for (const s of [1, -1]) { const F = 0.45 * i.gear_g * M * G0; cps.push({ node: stn[iw], off: [0, s * i.track_m / 2, -(0.5 * D + i.gear_len)], law: { ...gl, k: F / (0.35 * i.gear_stroke), Fy: F }, mu: wheels ? 0.05 : undefined, name: `${s > 0 ? 'Left' : 'Right'} main gear`, kind: 'gear' }); }
    const Fa = 0.2 * i.gear_g * M * G0; cps.push({ node: stn[auxI], off: [0, 0, -(0.5 * D * shape(auxI) + i.gear_len)], law: { ...gl, k: Fa / (0.35 * i.gear_stroke), Fy: Fa }, mu: wheels ? 0.05 : undefined, name: i.gear_type === 'tailwheel' ? 'Tail gear' : 'Nose gear', kind: 'gear' });
  }
  for (const w of lat) {
    if (w.tip !== undefined) { cps.push({ node: w.tip, off: [0, 0, -0.15 * D], name: 'Rotor arm tip', kind: 'tip' }); continue; }
    if (w.eng > 0 && i.wing_pos === 'low') { const F = 15 * G0 * w.eng, dn = 0.25 * D; cps.push({ node: w.mid, off: [0, 0, -0.35 * D], law: { k: F / (0.05 * dn), Fy: F, dmax: dn, kb: (20 * F) / (0.05 * dn) }, name: `${w.s > 0 ? 'Left' : 'Right'} engine nacelle`, kind: 'nacelle' }); }
    cps.push({ node: w.mid, name: `${w.s > 0 ? 'Left' : 'Right'} wing`, kind: 'tip' }, { node: w.out, name: `${w.s > 0 ? 'Left' : 'Right'} wing tip`, kind: 'tip' });
  }
  const up = Rb([0, 0, 1]); let zmin = Infinity;
  for (let k = 0; k < P.length; k++) P[k] = Rb(P[k]);
  for (const c of cps) { if (c.off) c.off = Rb(c.off); zmin = Math.min(zmin, P[c.node][2] + (c.off ? c.off[2] : 0)); }
  const gap = 1e-4 * D; for (const p of P) { p[2] += gap - zmin; fe.node(p[0], p[1], p[2]); }
  // ---- beams (their own mass is taken out of the lumped masses so the total is preserved)
  const share = new Array(P.length).fill(0), deg = new Array(P.length).fill(0); for (const b of beams) { deg[b.a]++; deg[b.b]++; }
  for (const b of beams) { const Lb = Math.hypot(P[b.b][0] - P[b.a][0], P[b.b][1] - P[b.a][1], P[b.b][2] - P[b.a][2]), mb = Math.min(m.rho * b.sec.A * Lb, tgt[b.a] / deg[b.a], tgt[b.b] / deg[b.b]); b.e = fe.beam(b.a, b.b, { ...b.sec, rho: mb / (b.sec.A * Lb) }, Rb(b.ref)); share[b.a] += 0.5 * mb; share[b.b] += 0.5 * mb; }
  // ---- occupants (one representative seat per chosen station)
  const pick = heli ? [N.clamp(Math.round(0.4 * (n - 1)) - 1, 0, n - 1), N.clamp(Math.round(0.4 * (n - 1)), 0, n - 1)] : multi || !i.occupied ? [iw] : [Math.round(0.15 * (n - 1)), iw, Math.round(0.75 * (n - 1))];
  const occK = [...new Set(pick)], mo = i.occupied ? 89 : Math.min(0.2 * M, Math.max(1e-3 * M, i.m_payload)), occ = occK.map((k) => { const mk = Math.min(mo, 0.5 * tgt[stn[k]]); tgt[stn[k]] -= mk; return { k, ...addOccupant(fe, stn[k], P[stn[k]], up, { ...seatSplit(mk, i.occupied, i.seat_g, i.s_ea), hgt: 0.25 * D }) }; });
  for (let k = 0; k < P.length; k++) fe.mass(k, Math.max(tgt[k] - share[k], 0.02 * tgt[k]));
  const pl = fe.plane({ mu: i.mu, fluid: i.surface === 'water' ? { rho: 1000, Cd: 1 } : null }); for (const c of cps) c.id = pl.contact(c.node, c);
  fe.velocity(-1, [i.v_fwd, 0, -i.v_sink]); fe.init();
  const E = fe.energy, mid = stn[iw], x0 = fe.x[3 * mid], wpB = () => { let s = 0; for (let e = 0; e < fe.nb; e++) s += fe.sec[e].wp || 0; return s; };
  const chans = [() => E.kinetic, () => E.internal, () => wpB() + E.crush, () => E.friction, () => E.external, () => E.total, () => pl.force, (f) => N.deg(Math.asin(N.clamp(f.R[9 * mid + 6], -1, 1))), (f) => f.x[3 * mid] - x0, (f) => f.v[3 * mid]];
  for (const o of occ) chans.push((f) => f.v[3 * stn[o.k] + 2], (f) => f.v[3 * stn[o.k]], (f) => f.v[3 * stn[o.k] + 1]);
  const rec = recorder(i.t_end, 500, chans), snaps = [], nSnap = 16, take = () => snaps.push({ t: fe.t, x: Array.from(fe.x) }); let errMax = 0, ohMin = Infinity, ohL0 = 0;
  const ohDist = () => (oh < 0 ? 0 : Math.min(...ohLegs.map((k) => Math.hypot(fe.x[3 * oh] - fe.x[3 * stn[k]], fe.x[3 * oh + 1] - fe.x[3 * stn[k] + 1], fe.x[3 * oh + 2] - fe.x[3 * stn[k] + 2]))));
  ohL0 = ohDist(); rec.take(fe); take();
  fe.run(i.t_end, () => { rec.take(fe); if (fe.t >= (snaps.length * i.t_end) / nSnap) take(); if (E.error > errMax) errMax = E.error; if (oh >= 0) ohMin = Math.min(ohMin, ohDist()); if (fe.nStep % 500 === 0) ctx?.progress?.(fe.t / i.t_end, 'Explicit integration'); });
  take();
  return { fe, rec, snaps, beams, occ, pl, cps, stn, lat, fin, oh, iw, n, errMax, ohLoss: oh >= 0 ? Math.max(0, 1 - ohMin / ohL0) : 0, wpB: wpB(), dc, mid, nNodes: P.length + 2 * occ.length };
}
const AC_BASE = { config: 'aeroplane', mass_kg: 78000, len_m: 37.6, dia_m: 3.95, span_m: 34.1, sweep_deg: 25, wing_pos: 'low', n_eng_wing: 2, m_engine: 2400, m_fuel_wing: 9000, wing_mass_frac: 0.11, m_nose: 0, m_overhead: 0, arm_m: 0.5, fin_h: 5.9, mount_g: 20, occupied: true, m_payload: 0, material: 'Al 2024-T3', t_fus: 0.0026, buckle_eff: 0.5, theta_break: 0.12, EI_wing: 1.5e9, A_wing: 0.05, Mp_wing: 2e7, d_crush: 0.47, crush_g: 10, gear: 'extended', gear_type: 'tricycle', gear_len: 1.2, gear_stroke: 0.45, gear_g: 2.5, track_m: 7.6, seat_g: 12, s_ea: 0.15, v_sink: 9.1, v_fwd: 60, pitch_deg: 5, roll_deg: 0, mu: 0.4, surface: 'rigid ground', filter_Hz: 60, t_end: 0.8, nStations: 7 };
const aircraft3d = {
  id: 'aircraft3d', title: 'Full-aircraft crash model (masses, nonlinear beams and crush springs)', fidelity: 'numerical',
  summary: 'A hybrid crash model of the whole vehicle in the KRASH tradition, generated from the case: fuselage stations, wings with engine and fuel masses (or rotor arms, or a transmission above the cabin), tail, landing gear, crushable belly and seated occupants, joined by elastic–plastic beams. It strikes the ground with forward and sink speed at a pitch and roll attitude; the run gives station decelerations, structural failure locations, the survivable-volume check, slide-out distance and deformed shapes.',
  equations: ['Conservation of linear momentum', 'Conservation of angular momentum', 'Conservation of energy', 'Transient structural dynamics equations', 'Plasticity constitutive equations', 'Impact contact equations', 'Crushing energy absorption equations', 'Occupant dynamics equations'],
  applicable: (c) => (c.fuselage.len_m > 0 && c.fuselage.dia_m > 0 ? true : 'No fuselage dimensions are defined for this vehicle.'),
  inputs: [
    { key: 'config', label: 'Configuration', type: 'select', options: CONFIGS, default: 'aeroplane', group: 'Vehicle', help: 'Aeroplane: wing with engines and fuel. Helicopter: rotor and transmission mass above the cabin. Multirotor: motor masses on arms' },
    { key: 'mass_kg', label: 'Aircraft mass at impact', unit: 'kg', default: 78000, min: 0.05, group: 'Vehicle' },
    { key: 'len_m', label: 'Fuselage length', unit: 'm', default: 37.6, min: 0.05, group: 'Vehicle' },
    { key: 'dia_m', label: 'Fuselage diameter', unit: 'm', default: 3.95, min: 0.02, group: 'Vehicle' },
    { key: 'span_m', label: 'Wing span (0 = no wing)', unit: 'm', default: 34.1, min: 0, group: 'Vehicle' },
    { key: 'sweep_deg', label: 'Wing sweep', unit: 'deg', default: 25, min: -30, max: 60, group: 'Vehicle' },
    { key: 'wing_pos', label: 'Wing position', type: 'select', options: ['low', 'high'], default: 'low', group: 'Vehicle', help: 'A low wing with underslung engines meets the ground before the belly' },
    { key: 'n_eng_wing', label: 'Engines on the wing (or rotors on arms)', unit: '', default: 2, min: 0, max: 8, step: 1, discrete: true, group: 'Vehicle' },
    { key: 'm_engine', label: 'Mass of each wing engine or arm motor', unit: 'kg', default: 2400, min: 0, group: 'Vehicle' },
    { key: 'm_fuel_wing', label: 'Fuel or battery mass in the wing', unit: 'kg', default: 9000, min: 0, group: 'Vehicle' },
    { key: 'wing_mass_frac', label: 'Wing structure / aircraft mass', unit: '-', default: 0.11, min: 0.01, max: 0.3, group: 'Vehicle' },
    { key: 'm_nose', label: 'Concentrated mass at the nose (engine, battery, sensor)', unit: 'kg', default: 0, min: 0, group: 'Vehicle' },
    { key: 'm_overhead', label: 'Rotor, hub and transmission mass above the cabin', unit: 'kg', default: 0, min: 0, group: 'Vehicle', help: 'Helicopter configuration only' },
    { key: 'mount_g', label: 'Transmission mount strength', unit: 'g', default: 20, min: 1, max: 100, group: 'Vehicle', help: 'Vertical load factor at which the mounts yield' },
    { key: 'arm_m', label: 'Rotor arm length (multirotor)', unit: 'm', default: 0.5, min: 0.01, group: 'Vehicle' },
    { key: 'fin_h', label: 'Fin height', unit: 'm', default: 5.9, min: 0, group: 'Vehicle' },
    { key: 'occupied', label: 'Occupied (evaluate injury metrics)', type: 'bool', default: true, group: 'Vehicle' },
    { key: 'm_payload', label: 'Payload mass on the monitored mount (unoccupied vehicles)', unit: 'kg', default: 0, min: 0, group: 'Vehicle' },
    MAT,
    { key: 't_fus', label: 'Equivalent fuselage shell thickness', unit: 'm', default: 0.0026, min: 1e-4, max: 0.05, group: 'Structure', help: 'Skin plus smeared stringers' },
    { key: 'buckle_eff', label: 'Fuselage plastic-moment efficiency', unit: '-', default: 0.5, min: 0.05, max: 1, group: 'Structure', help: 'Shell buckling knock-down on the fully plastic bending moment of the tube; calibrate against test or detailed analysis' },
    { key: 'theta_break', label: 'Hinge rotation at fuselage break', unit: 'rad', default: 0.12, min: 0.01, max: 1, group: 'Structure' },
    { key: 'EI_wing', label: 'Wing root bending stiffness', unit: 'N·m²', default: 1.5e9, min: 1e-3, group: 'Structure', help: 'From the structures suite when available' },
    { key: 'A_wing', label: 'Wing box material area at the root', unit: 'm²', default: 0.05, min: 1e-7, group: 'Structure' },
    { key: 'Mp_wing', label: 'Wing root plastic bending moment', unit: 'N·m', default: 2e7, min: 1e-3, group: 'Structure' },
    { key: 'd_crush', label: 'Crushable depth below the floor', unit: 'm', default: 0.47, min: 0.001, group: 'Crush', help: 'Subfloor structure available before the stiff cabin floor is reached' },
    { key: 'crush_g', label: 'Belly crush strength (load factor on the supported mass)', unit: 'g', default: 10, min: 0.5, max: 200, group: 'Crush', help: 'Plateau force of the lower fuselage divided by the weight it carries; test-derived in practice' },
    { key: 'gear', label: 'Landing gear', type: 'select', options: ['extended', 'retracted'], default: 'extended', group: 'Crush' },
    { key: 'gear_type', label: 'Gear arrangement', type: 'select', options: ['tricycle', 'tailwheel', 'skid'], default: 'tricycle', group: 'Crush' },
    { key: 'gear_len', label: 'Gear length below the belly', unit: 'm', default: 1.2, min: 0.001, group: 'Crush' },
    { key: 'gear_stroke', label: 'Gear stroke', unit: 'm', default: 0.45, min: 0.001, group: 'Crush' },
    { key: 'gear_g', label: 'Gear reaction factor at the stroke plateau', unit: 'g', default: 2.5, min: 0.2, max: 30, group: 'Crush', help: 'Total gear load / weight while stroking; the gear collapses when the stroke is used up' },
    { key: 'track_m', label: 'Main gear track', unit: 'm', default: 7.6, min: 0.01, group: 'Crush' },
    { key: 'seat_g', label: 'Seat energy-absorber limit load factor', unit: 'g', default: 12, min: 1, max: 200, group: 'Seat' },
    { key: 's_ea', label: 'Available seat stroke', unit: 'm', default: 0.15, min: 0, max: 1, group: 'Seat' },
    { key: 'v_sink', label: 'Sink speed at impact', unit: 'm/s', default: 9.1, min: 0.1, max: 40, group: 'Impact', help: 'Rotorcraft default: the 9.14 m/s (30 ft/s) seat-test velocity change. Aeroplane default from the case: twice the landing-gear design sink rate (four times its design energy), an illustrative severe-but-survivable emergency landing rather than a regulatory value' },
    { key: 'v_fwd', label: 'Forward speed at impact', unit: 'm/s', default: 60, min: 0, max: 150, group: 'Impact' },
    { key: 'pitch_deg', label: 'Pitch attitude (nose up +)', unit: 'deg', default: 5, min: -30, max: 30, group: 'Impact' },
    { key: 'roll_deg', label: 'Roll attitude (right side down +)', unit: 'deg', default: 0, min: -30, max: 30, group: 'Impact' },
    { key: 'mu', label: 'Friction coefficient of structure on the surface', unit: '-', default: 0.4, min: 0, max: 1.5, group: 'Impact', help: 'About 0.35–0.5 for metal on concrete or soil; rolling wheels use 0.05 until the gear collapses' },
    { key: 'surface', label: 'Impact surface', type: 'select', options: ['rigid ground', 'water'], default: 'rigid ground', group: 'Impact', help: 'Water: a pressure surface (½ρV² plus buoyancy on the belly area) instead of crush springs' },
    { key: 'filter_Hz', label: 'Low-pass filter for deceleration histories', unit: 'Hz', default: 60, min: 1, max: 1e5, group: 'Numerics' },
    { key: 't_end', label: 'Simulated time', unit: 's', default: 0.8, min: 1e-3, max: 10, group: 'Numerics' },
    { key: 'nStations', label: 'Fuselage mass stations', unit: '', default: 7, min: 3, max: 41, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c, up) {
    const M = c.mass.mtow_kg, D = c.fuselage.dia_m, L = c.fuselage.len_m, occ = occupants(c), has = occ > 0, heli = c.meta.type === 'helicopter', multi = !heli && !(c.wing.S_m2 > 0), mm = mat(c.struct.material), p = c.prop, rot = heli || c.meta.type === 'evtol';
    const mEng = p.type === 'turbofan' ? p.T0_N / (5 * G0) : p.type === 'electric' ? (p.P0_W || 1e3 * c.systems.motor_kW) / 2500 + 0.02 * M / Math.max(1, p.n_eng) : (p.P0_W / (p.type === 'piston' ? 900 : 3500)) * 1.4;
    const onWing = !heli && !multi && p.n_eng >= 2, d = N.clamp(0.12 * D, 0.01, 0.5), v = has && heli ? seatTest(c).dv : has ? Math.min(seatTest(c).dv, 2 * c.gear.sink_ms) : Math.max(2, 2 * c.gear.sink_ms), S = c.wing.S_m2, Vs = S > 0 ? Math.sqrt((2 * M * G0) / (1.225 * S * c.aero.CLmax_land)) : 0;
    const cr = S > 0 ? (2 * S) / (c.wing.b_m * (1 + c.wing.taper)) : 0, w = c.struct.box_chord_frac * cr, h = c.struct.box_height_frac * c.wing.tc * cr, ts = c.struct.t_skin_mm / 1e3, tw = c.struct.t_spar_mm / 1e3, EI = mm.E * (2 * 1.6 * w * ts * (h / 2) ** 2 + (2 * tw * h ** 3) / 12);
    const gl = c.gear.type === 'skid' ? Math.max(0.3 * D, 2 * c.gear.stroke_m) : Math.max(0.2 * D, c.gear.stroke_m + 0.6 * (c.gear.tyre_dia_m || 0.2 * D)), te = (2.5 * (gl + d + 0.5 * L * Math.sin(N.rad(5)) + 0.3 * D)) / v + 0.08;
    return {
      config: heli ? 'helicopter' : multi ? 'multirotor' : 'aeroplane', mass_kg: M, len_m: L, dia_m: D, span_m: multi || heli ? 0 : c.wing.b_m, sweep_deg: c.wing.sweep_deg, wing_pos: p.type === 'turbofan' || c.meta.type === 'evtol' ? 'low' : 'high',
      n_eng_wing: multi ? p.n_eng : onWing ? Math.min(4, 2 * Math.floor(p.n_eng / 2)) : 0, m_engine: Math.max(1e-3 * M, mEng), m_fuel_wing: S > 0 ? 0.5 * (c.mass.fuel_kg || 0) : 0, m_nose: !heli && !multi && !onWing ? Math.min(0.25 * M, mEng * p.n_eng) : 0,
      m_overhead: heli ? 0.12 * M + c.rotor.n_blades * c.rotor.blade_mass_kg : 0, arm_m: multi ? Math.max(0.6 * L, 1.4 * (c.rotor.R_m || 0.5 * p.prop_dia_m)) : undefined, fin_h: Math.max(0.2 * D, c.vtail.b_m || 0.5 * D), occupied: has, m_payload: has ? 0 : c.mass.payload_kg || 0.1 * M,
      material: c.struct.material, t_fus: N.clamp(6.5e-4 * D, 5e-4, 6e-3), EI_wing: S > 0 ? up.fea?.EI_root_Nm2 ?? EI : undefined, A_wing: S > 0 ? 2 * 1.6 * w * ts + 2 * h * tw : undefined, Mp_wing: S > 0 ? 0.6 * mm.Sy * 1.6 * w * ts * h : undefined,
      d_crush: d, crush_g: Math.max(8, (v * v) / (1.6 * d * G0) + 1), gear_type: c.gear.type, gear_len: gl, gear_stroke: Math.max(1e-3, c.gear.stroke_m), gear_g: rot ? 4 : 2.5, track_m: Math.max(0.1 * D, c.gear.track_m), seat_g: has ? 12 : 30, s_ea: !has ? 0.02 * D : seatStroke(v, 0.8 * (d + Math.max(1e-3, c.gear.stroke_m)), 12, rot ? 0.3 : 0.15),
      v_sink: v, v_fwd: heli || multi ? Math.min(15, c.flight.V_ms) : N.clamp(1.1 * (up.performance?.V_stall_ms ?? Vs), 5, 90), pitch_deg: heli || multi ? 0 : 5, filter_Hz: N.clamp((6 * v) / (d + 0.3 * gl), 20, 5000), t_end: te,
    };
  },
  run(i, ctx) {
    const r = aircraftSim(i, ctx), fe = r.fe, E = fe.energy, rec = r.rec, T = rec.t, fs = 1 / rec.dt, lp = (v) => lowpass(ddt(v, rec.dt), fs, i.filter_Hz), warnings = [RES_NOTE], th = (a) => thin(a, 250), ms = T.map((t) => t * 1e3);
    const st = r.occ.map((o, q) => { const az = lp(rec.data[10 + 3 * q]).map((a) => a / G0 + 1), ax = lp(rec.data[11 + 3 * q]).map((a) => -a / G0), ay = lp(rec.data[12 + 3 * q]).map((a) => a / G0), res = az.map((a, k) => Math.hypot(a, ax[k], ay[k]));
      return { k: o.k, az, ax, pkV: N.amax(az), pkL: N.amax(ax), pk: N.amax(res), dri: (DRI_W ** 2 * fe.springs[o.sSpine].dMaxSeen) / G0, stroke: fe.springs[o.sSeat].dpMax }; });
    const kap = Array.from(fe.beamKappa), rows = [], pts = r.pl.pts, belly = r.cps.filter((c) => c.kind === 'belly'), crushOf = (c) => Math.max(pts[c.id].dp, 0);
    let cabinBreak = false;
    for (const b of r.beams) { const k = Math.max(kap[2 * b.e], kap[2 * b.e + 1]), dead = !fe.beamAlive[b.e]; if (dead || k > 0.01) { const brk = k >= b.sec.thetaSoft; if (brk && b.cabin && r.occ.some((o) => Math.abs(o.k - (b.a + 0.5)) < 1.6)) cabinBreak = true; rows.push([b.name, dead ? 'Fracture (member lost)' : brk ? 'Hinge beyond break rotation' : 'Plastic hinge', k, (dead ? fe.beamFailT[b.e] : fe.beamYieldT[b.e]) * 1e3]); } }
    for (const c of r.cps) { const p = pts[c.id]; if (p.failed) rows.push([c.name, 'Gear collapsed', p.dMax, p.tFail * 1e3]); else if (c.law && p.dMax > c.law.dmax) rows.push([c.name, c.kind === 'belly' ? 'Crush depth used up (bottomed)' : 'Bottomed', p.dMax, NaN]); }
    const crush = belly.length && i.surface !== 'water' ? Math.max(...belly.map(crushOf)) : 0, intr = i.surface === 'water' ? 0 : Math.max(0, ...belly.filter((c) => r.occ.some((o) => o.k === c.stn)).map((c) => pts[c.id].dMax - r.dc(c.stn))), loss = Math.max(intr / Math.max(i.dia_m - i.d_crush, 1e-9), r.ohLoss), ok = loss < 0.15 && !cabinBreak;
    const vxEnd = rec.data[9][rec.data[9].length - 1], dist = rec.data[8][rec.data[8].length - 1], slide = dist + (i.mu > 0 ? (vxEnd * Math.abs(vxEnd)) / (2 * i.mu * G0) : Infinity), Ein = E.initial + E.external, absorbed = r.wpB + E.crush + E.friction, pk = N.amax(st.map((s) => s.pk)), dri = N.amax(st.map((s) => s.dri)), vz = rec.data[10][rec.data[10].length - 1];
    if (r.n < Math.round(i.nStations)) warnings.push(`The fuselage is stubby: ${r.n} mass stations are used instead of ${Math.round(i.nStations)} so that each beam element stays longer than three quarters of a diameter.`);
    if (i.config === 'aeroplane' && !(i.span_m > 0)) warnings.push('Aeroplane configuration without a wing span: only the fuselage and tail are modelled.');
    if (Math.abs(vz) > 0.25 * i.v_sink) warnings.push('Vertical motion has not settled by the end of the run (the vehicle is still crushing, bouncing or slapping down): extend the simulated time.');
    if (r.errMax > 0.05) warnings.push(`The energy balance error reached ${(100 * r.errMax).toFixed(1)}%: reduce the time step by raising the number of stations or shortening the run.`);
    if (rows.some((q) => /bottomed/i.test(q[1]))) warnings.push('The crushable depth is used up at one or more stations. The deceleration peaks and DRI after bottoming are set by the assumed bottoming stiffness (20 × the crush stiffness) and by how the mass is lumped into stations, so read them as "exceeds tolerance" rather than as values; crush distance, energy absorbed and slide-out are far less sensitive.');
    if (i.surface === 'water') warnings.push('Water is a pressure surface (½ρV² and buoyancy on the belly area): no planing lift, suction, cavity or skin rupture. Use the ditching analysis for the wedge-entry loads.');
    if (!i.occupied) warnings.push('No occupants in this case: the monitored mass is the payload on a 52.9 rad/s mount and injury metrics are not meaningful.');
    if (i.mu > 0 && vxEnd > 1) warnings.push(`The vehicle is still sliding at ${vxEnd.toFixed(1)} m/s when the run ends: the slide-out distance adds a constant-friction extrapolation of ${((vxEnd * vxEnd) / (2 * i.mu * G0)).toFixed(1)} m and ignores obstacles, ploughing and break-up.`);
    // deformed shapes in a frame moving with the centre station (x) and fixed to the ground (z)
    const pickS = [0, 0.15, 0.35, 0.6, 1].map((f) => r.snaps[Math.min(r.snaps.length - 1, Math.round(f * (r.snaps.length - 1)))]), sideOf = (s) => { const x = [], z = [], xm = s.x[3 * r.mid], put = (k) => { x.push(s.x[3 * k] - xm); z.push(s.x[3 * k + 2]); }, brk = () => { x.push(NaN); z.push(NaN); };
      r.stn.forEach(put); if (r.fin >= 0) put(r.fin); if (r.oh >= 0) { brk(); put(r.stn[N.clamp(Math.round(0.4 * (r.n - 1)) - 1, 0, r.n - 1)]); put(r.oh); put(r.stn[N.clamp(Math.round(0.4 * (r.n - 1)) + 1, 0, r.n - 1)]); } return { x, y: z }; };
    const frontOf = (s) => { const y = [], z = [], put = (k) => { y.push(s.x[3 * k + 1]); z.push(s.x[3 * k + 2]); }, L = r.lat.filter((w) => w.s > 0), Rr = r.lat.filter((w) => w.s < 0);
      if (i.config === 'multirotor') { for (const w of r.lat) { put(r.mid); put(w.tip); y.push(NaN); z.push(NaN); } } else { for (const w of L) { put(w.out); put(w.mid); } put(r.mid); for (const w of Rr) { put(w.mid); put(w.out); } } return { x: y, y: z }; };
    const span = Math.max(i.len_m, i.span_m, 2 * i.arm_m * (i.config === 'multirotor' ? 1 : 0)), ground = (a) => ({ name: 'Ground', x: [-0.55 * a, 0.55 * a], y: [0, 0], style: 'dash' });
    const names = st.map((s) => `Station ${s.k + 1}`);
    return {
      kpis: [
        kpi('ac_peak_g', 'Peak resultant deceleration at the monitored stations', pk, 'g', undefined, `Low-pass filtered at ${i.filter_Hz.toFixed(0)} Hz`),
        kpi('ac_peak_vert_g', 'Peak vertical deceleration', N.amax(st.map((s) => s.pkV)), 'g'), kpi('ac_peak_long_g', 'Peak longitudinal deceleration', N.amax(st.map((s) => s.pkL)), 'g'),
        kpi('ac_DRI', 'Highest Dynamic Response Index', dri, '-', !i.occupied ? undefined : dri <= 18 ? 'ok' : dri <= 22.8 ? 'warn' : 'bad', 'About 18 corresponds to roughly 5% and 22.8 to 50% probability of spinal injury'),
        kpi('ac_crush_m', 'Largest belly crush', crush, 'm', crush > i.d_crush ? 'bad' : crush > 0.85 * i.d_crush ? 'warn' : 'ok', `Crushable depth ${i.d_crush.toFixed(3)} m`),
        kpi('ac_volume_loss_pct', 'Loss of cabin height at the occupied stations', 100 * loss, '%', loss < 0.15 ? 'ok' : 'bad', 'Floor intrusion after the crush depth is used up, or descent of the overhead mass; 15% is the customary survivable-volume limit'),
        kpi('ac_survivable', 'Survivable volume maintained (1 = yes)', ok ? 1 : 0, '-', ok ? 'ok' : 'bad', cabinBreak ? 'A fuselage break lies in the occupied cabin' : 'No fuselage break in the occupied cabin'),
        kpi('ac_failures', 'Structural failure locations', rows.length, '-'), kpi('ac_slide_m', 'Slide-out distance', slide, 'm', undefined, 'Distance travelled in the run plus v²/(2μg) for the remaining speed'),
        kpi('ac_absorbed_J', 'Energy absorbed (crush, plastic hinges, friction)', absorbed, 'J'), kpi('ac_absorbed_frac', 'Share of the impact energy absorbed in the run', absorbed / Ein, '-'),
        kpi('ac_seat_stroke_m', i.occupied ? 'Largest seat stroke' : 'Largest payload-mount stroke', N.amax(st.map((s) => s.stroke)), 'm', N.amax(st.map((s) => s.stroke)) > i.s_ea && i.s_ea > 0 ? (i.occupied ? 'bad' : 'warn') : 'ok', `Available ${i.s_ea.toFixed(3)} m`),
        kpi('ac_ground_force_N', 'Peak ground reaction', N.amax(rec.data[6]), 'N'), kpi('ac_v_fwd_end_ms', 'Forward speed at the end of the run', vxEnd, 'm/s'),
        kpi('ac_energy_err', 'Energy balance error', r.errMax, '-', r.errMax < 0.03 ? 'ok' : 'warn', 'Kinetic + internal + contact − gravity work against the initial energy'),
        kpi('ac_mass_kg', 'Mass of the model', fe.totalMass, 'kg'), kpi('ac_elements', 'Beam elements', fe.nb, '-'), kpi('ac_steps', 'Time steps', fe.nStep, '-'),
      ],
      plots: [
        { type: 'line', title: 'Deformed shape, side view (frame moving with the centre station)', xlabel: 'Fore–aft position [m]', ylabel: 'Height [m]', equalAspect: true, series: [...pickS.map((s) => ({ name: `${(s.t * 1e3).toFixed(0)} ms`, ...sideOf(s), style: 'line+points' })), ground(i.len_m)] },
        ...(r.lat.length ? [{ type: 'line', title: 'Deformed shape, front view', xlabel: 'Lateral position [m]', ylabel: 'Height [m]', equalAspect: true, series: [...pickS.map((s) => ({ name: `${(s.t * 1e3).toFixed(0)} ms`, ...frontOf(s), style: 'line+points' })), ground(span)] }] : []),
        { type: 'line', title: 'Vertical deceleration at the monitored stations', xlabel: 'Time [ms]', ylabel: 'Deceleration [g]', series: st.map((s, q) => ({ name: names[q], x: th(ms), y: th(s.az) })) },
        { type: 'line', title: 'Longitudinal deceleration at the monitored stations', xlabel: 'Time [ms]', ylabel: 'Deceleration [g]', series: st.map((s, q) => ({ name: names[q], x: th(ms), y: th(s.ax) })) },
        { type: 'line', title: 'Energy balance', xlabel: 'Time [ms]', ylabel: 'Energy [kJ]', series: [['Kinetic', 0], ['Internal (elastic + plastic)', 1], ['Crush and plastic work', 2], ['Friction', 3], ['Gravity work', 4], ['Total (should stay constant)', 5]].map(([name, k]) => ({ name, x: th(ms), y: th(rec.data[k]).map((e) => e / 1e3), ...(k === 5 ? { style: 'dash' } : {}) })) },
        { type: 'line', title: 'Pitch attitude of the centre fuselage', xlabel: 'Time [ms]', ylabel: 'Pitch change [deg]', series: [{ name: 'Nose up +', x: th(ms), y: th(rec.data[7]) }] },
      ],
      tables: [
        { title: 'Structural failure locations', columns: ['Member', 'Event', 'Plastic rotation [rad] or crush [m]', 'Time [ms]'], rows: rows.length ? rows : [['None', 'No yielding, bottoming or collapse', 0, NaN]] },
        { title: 'Monitored stations', columns: ['Station', 'Peak vertical [g]', 'Peak longitudinal [g]', 'DRI', 'Seat stroke [m]'], rows: st.map((s) => [s.k + 1, s.pkV, s.pkL, s.dri, s.stroke]) },
      ],
      outputs: { ac_energy_in_J: Ein },
      warnings,
      models: [`Hybrid crash model: ${r.nNodes} mass points, ${fe.nb} co-rotational elastic–plastic beams, ${r.cps.length} external crush springs or contact points, ${fe.nStep} explicit steps`, 'Plastic hinges with axial-force interaction; fuselage break and wing failure as residual-moment hinges, transmission mounts with axial fracture', i.surface === 'water' ? 'Water as a pressure surface with drag and buoyancy' : 'Elastic–plastic crush springs with bottoming on rigid ground, Coulomb friction, collapsing landing gear', 'Seat energy absorber and Dynamic Response Index spinal model at the monitored stations'],
      assumptions: ['Mass and stiffness distributions are generated from overall dimensions and class-typical fractions, not from a drawing: use it to rank scenarios and locate weak points, not to certify', 'Crush springs act along body-fixed offsets and carry test-derived plateau strengths supplied as inputs', 'Seats stroke along the initial body-vertical direction; no restraint, flail or head strike', 'Rigid, flat surface without obstacles; no fuel spillage, fire or post-crash break-up dynamics', DATA_NOTE],
    };
  },
  convergence: { param: 'nStations', label: 'Fuselage mass stations', levels: [5, 7, 9, 13], metric: 'ac_crush_m' },
  calibration: { params: [{ key: 'crush_g', min: 1, max: 100 }, { key: 'buckle_eff', min: 0.1, max: 1 }, { key: 'mu', min: 0.05, max: 1 }], sweep: 'v_sink', target: 'ac_peak_vert_g', note: 'Full-scale or section drop-test floor decelerations at several sink speeds calibrate the belly crush strength; slide-out distances calibrate the friction coefficient.' },
  verify: () => aircraftChecks(),
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!o.ac_survivable) out.push({ severity: 'critical', title: 'Survivable volume is lost', detail: `Cabin height loss ${o.ac_volume_loss_pct.toFixed(0)}% with ${o.ac_failures} failure location(s).`, action: i.config === 'helicopter' ? 'Strengthen the transmission and engine mounts for the vertical crash load factor, and deepen the subfloor so that it does not bottom out.' : 'Add crushable depth below the floor or raise its crush strength, and reinforce the fuselage at the break location (usually just aft of the wing box).', basis: 'Occupied volume reduced by no more than about 15%' });
    if (i.occupied && o.ac_DRI > 18) out.push({ severity: o.ac_DRI > 22.8 ? 'critical' : 'warn', title: 'High spinal injury risk', detail: `DRI ${o.ac_DRI.toFixed(1)}, vertical peak ${o.ac_peak_vert_g.toFixed(0)} g.`, action: 'Use stroking seats, lower the belly crush strength towards the level that just avoids bottoming, and keep the landing gear as the first energy-absorbing stage.', basis: 'DRI ≈ 18 for about 5% spinal injury probability' });
    if (o.ac_crush_m > i.d_crush) out.push({ severity: 'warn', title: 'Belly crush depth is used up', detail: `${o.ac_crush_m.toFixed(2)} m of crush against ${i.d_crush.toFixed(2)} m available.`, action: `Raise the crush strength towards ${(o.ac_energy_in_J / (0.8 * i.d_crush * i.mass_kg * G0)).toFixed(0)} g or deepen the subfloor; check the pitched and rolled cases, which concentrate the load on a few stations.`, basis: 'Energy absorbed = plateau force × crush distance' });
    out.push({ severity: 'info', title: 'Scenario summary', detail: `Slide-out about ${Number.isFinite(o.ac_slide_m) ? o.ac_slide_m.toFixed(0) : '∞'} m, ${(100 * o.ac_absorbed_frac).toFixed(0)}% of the energy absorbed in the run, energy balance error ${(100 * o.ac_energy_err).toFixed(2)}%.`, action: 'Repeat for the envelope of sink speed, pitch and roll, on soft ground and water, then pass the governing case to a detailed explicit model through the High-fidelity bridge.', basis: 'Crash scenario envelope' });
    return out;
  },
};

// ---- 10. rigid or soft-equivalent impactor on a plate of solid elements ---------------------
const PROJ = { 'bird (soft-body equivalent)': { rho: 950, soft: 'stagnation' }, 'hail (ice)': { rho: 900, soft: 'crush' }, 'hard fragment (steel)': { rho: 7850, soft: null } };
const RATE = ['Johnson–Cook', 'Cowper–Symonds', 'rate-independent'];
function impactSim(i, ctx, o = {}) {
  const m = mat(i.material), jc = jcOf(m), pr = PROJ[i.projectile] || PROJ['hard fragment (steel)'], Rs = ((3 * i.m_proj) / (4 * Math.PI * pr.rho)) ** (1 / 3), half = i.a_plate / 2, ne = Math.max(2, Math.round(i.nEl)), h = half / ne;
  const K = m.E / (3 * (1 - 2 * m.nu)), G = m.E / (2 * (1 + m.nu)), cs = Math.sqrt(Math.max(K + (4 / 3) * G, 3 * K) / m.rho), tSim = (i.t_factor * 2 * Rs) / i.V;
  // thickness handling: true thickness when affordable, otherwise one membrane-equivalent layer of scaled thickness
  let nz = Math.max(1, Math.round(i.nThk)), sc = 1; const need = tSim / ((0.9 * i.t_plate) / nz / cs);
  if (need > i.max_steps && !o.noScale) { nz = 1; sc = N.clamp(((tSim / i.max_steps) * cs) / 0.9 / i.t_plate, 1, Math.max(1, h / i.t_plate)); }
  const Tp = sc * i.t_plate, fe = createFE({ hourglass: 0.05 / (sc * sc) }), rate = i.rate_model === 'Johnson–Cook' ? { C: jc.C } : i.rate_model === 'Cowper–Symonds' ? { csD: i.cs_D, csQ: i.cs_q } : {};
  const mt = fe.material({ E: m.E / sc, nu: m.nu, rho: m.rho / sc, A: o.elastic ? Infinity : jc.A / sc, B: jc.B / sc, n: jc.n, m: jc.m, Tm: i.rate_model === 'Johnson–Cook' ? jc.Tm : 0, cp: m.cp, efail: o.elastic ? 0 : i.eps_fail, ...rate });
  const Rc = i.R_curv, blk = hexBlock(fe, ne, ne, nz, half, half, Tp, mt, Rc > 0 ? (p) => { const r = Rc + p[2] - Tp, a = p[0] / Rc; return [r * Math.sin(a), p[1], r * Math.cos(a) - Rc + Tp]; } : undefined);
  const top = [], caps = [], pc = pr.soft === 'stagnation' ? 0.5 * pr.rho * i.V * i.V : pr.soft === 'crush' ? i.p_crush : Infinity;
  for (let a = 0; a <= ne; a++) for (let b = 0; b <= ne; b++) {
    for (let k = 0; k <= nz; k++) { const nd = blk.id(a, b, k), edge = a === ne || b === ne; if (edge && i.clamped && !o.free) fe.fix(nd, [1, 1, 1]); else if (a === 0 || b === 0) fe.fix(nd, [a === 0 ? 1 : 0, b === 0 ? 1 : 0, 0]); }
    top.push(blk.id(a, b, nz)); caps.push(pc * h * h * (a === 0 || a === ne ? 0.5 : 1) * (b === 0 || b === ne ? 0.5 : 1));
  }
  const gap = 1e-3 * Rs, imp = fe.impactor({ c: [0, 0, Tp + Rs + gap], v: [0, 0, -i.V], R: Rs, mass: i.m_proj / 4, free: [0, 0, 1], nodes: Number.isFinite(pc) ? top : blk.nd, ...(Number.isFinite(pc) ? { fcaps: caps } : {}) });
  fe.init();
  const cN = blk.id(0, 0, 0), E = fe.energy, rec = recorder(tSim, 300, [() => -imp.v[2], () => 4 * imp.force, (f) => -(f.x[3 * cN + 2] - f.x0[3 * cN + 2]), () => E.kinetic, () => E.internal + E.viscous, () => E.plastic, () => E.hourglass, () => E.contact, () => E.total]);
  let errMax = 0, Fpk = 0, impulse = 0, tSep = NaN, hit = false, defl = 0; rec.take(fe);
  fe.run(tSim, () => {
    rec.take(fe); if (E.error > errMax) errMax = E.error; const F = 4 * imp.force; if (F > Fpk) Fpk = F; impulse += F * fe.dtLast; const d = -(fe.x[3 * cN + 2] - fe.x0[3 * cN + 2]); if (d > defl && fe.hexAlive[blk.els[0]]) defl = d;
    if (F > 0) hit = true; else if (hit && !(tSep > 0) && imp.v[2] > 0) tSep = fe.t;
    if (fe.nStep % 200 === 0) ctx?.progress?.(fe.t / tSim, 'Explicit integration');
    return (tSep > 0 && fe.t > 1.15 * tSep) || imp.c[2] < -3 * Rs - Tp || (hit && F === 0 && d < 0.9 * defl && imp.c[2] < -Rs); // rebounded, passed through, or passed with the dent recovering
  }, o.maxSteps || 40 * i.max_steps);
  // through-thickness erosion of any column of elements = perforation
  let holes = 0, epMax = 0; const col = (a, b, k) => blk.els[(a * ne + b) * nz + k];
  for (let a = 0; a < ne; a++) for (let b = 0; b < ne; b++) { let dead = 0; for (let k = 0; k < nz; k++) { const e = col(a, b, k); if (!fe.hexAlive[e]) dead++; if (fe.ep[e] > epMax) epMax = fe.ep[e]; } if (dead === nz) holes++; }
  return { fe, rec, imp, blk, Rs, sc, nz, ne, h, Tp, tSim, errMax, Fpk, impulse, tSep, defl, holes, epMax, pc, col, jc, need };
}
const IMPACT_BASE = { projectile: 'hard fragment (steel)', m_proj: 0.025, V: 100, material: 'Al 2024-T3', t_plate: 0.003, a_plate: 0.11, R_curv: 0, clamped: true, eps_fail: 0.18, rate_model: 'Johnson–Cook', cs_D: 6500, cs_q: 4, p_crush: 1e7, t_factor: 2.5, max_steps: 1500, nThk: 2, nEl: 6 };
const impact3d = {
  id: 'impact3d', title: 'Projectile impact on a skin panel (explicit 3-D solid elements with erosion)', fidelity: 'numerical',
  summary: 'A sphere — a hard fragment, a hailstone or a soft-body equivalent of a bird — strikes the centre of a flat or curved panel meshed with solid elements. Johnson–Cook plasticity and erosion at a failure strain decide whether the panel is perforated; the run gives the residual velocity, contact force, dent depth and plastic strain, next to the analytical energy screening.',
  equations: ['Conservation of linear momentum', 'Conservation of energy', 'Transient structural dynamics equations', 'Plasticity constitutive equations', 'Strain-rate-dependent material equations', 'Impact contact equations', 'Fracture criteria'],
  inputs: [
    { key: 'projectile', label: 'Projectile', type: 'select', options: Object.keys(PROJ), default: 'bird (soft-body equivalent)', group: 'Projectile', help: 'Hard fragment: rigid sphere. Hail: rigid sphere whose contact pressure is capped at the ice crushing strength. Bird: sphere whose contact pressure is capped at the stagnation pressure ½ρV²' },
    { key: 'm_proj', label: 'Projectile mass', unit: 'kg', default: 1.81, min: 1e-5, max: 20, group: 'Projectile' },
    { key: 'V', label: 'Impact speed (normal to the panel)', unit: 'm/s', default: 180, min: 1, max: 1500, group: 'Projectile' },
    { key: 'p_crush', label: 'Ice crushing strength (hail only)', unit: 'Pa', default: 1e7, min: 1e5, max: 1e9, group: 'Projectile', help: 'Order of 10 MPa at impact rates; scatter is large' },
    MAT,
    { key: 't_plate', label: 'Panel thickness', unit: 'm', default: 0.003, min: 1e-4, max: 0.1, group: 'Panel' },
    { key: 'a_plate', label: 'Panel side length', unit: 'm', default: 0.6, min: 0.005, max: 5, group: 'Panel', help: 'Square bay between supports; at least three projectile diameters' },
    { key: 'R_curv', label: 'Radius of curvature (0 = flat)', unit: 'm', default: 0, min: 0, max: 100, group: 'Panel', help: 'Cylindrical panel, convex towards the projectile (leading edge or fuselage skin)' },
    { key: 'clamped', label: 'Edges clamped', type: 'bool', default: true, group: 'Panel' },
    { key: 'eps_fail', label: 'Equivalent plastic strain at failure (erosion)', unit: '-', default: 0.18, min: 0.01, max: 2, group: 'Material', help: 'Mesh-dependent erosion strain; calibrate against a penetration test at the same element size' },
    { key: 'rate_model', label: 'Strain-rate model', type: 'select', options: RATE, default: 'Johnson–Cook', group: 'Material' },
    { key: 'cs_D', label: 'Cowper–Symonds D', unit: '1/s', default: 6500, min: 1, group: 'Material' },
    { key: 'cs_q', label: 'Cowper–Symonds q', unit: '-', default: 4, min: 1, max: 20, group: 'Material' },
    { key: 't_factor', label: 'Simulated time / projectile transit time 2R/V', unit: '-', default: 2.5, min: 0.5, max: 50, group: 'Numerics', help: 'The run stops earlier once the projectile has left the panel and the dent is rebounding' },
    { key: 'max_steps', label: 'Time-step budget for the true thickness', unit: '', default: 1500, min: 200, max: 2e6, step: 1, discrete: true, group: 'Numerics', help: 'If resolving the true thickness needs more steps, one layer of scaled thickness with scaled modulus, strength and density is used (same membrane stiffness, strength and mass per area)' },
    { key: 'nThk', label: 'Elements through the thickness', unit: '', default: 2, min: 1, max: 8, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nEl', label: 'Elements along the half side', unit: '', default: 6, min: 4, max: 40, step: 1, discrete: true, group: 'Numerics', help: 'A quarter of the panel is modelled using the two symmetry planes' },
  ],
  defaults(c, up) {
    const mb = c.mass.mtow_kg < 100 ? 0.45 : 1.81, R = ((3 * mb) / (4 * Math.PI * 950)) ** (1 / 3), mm = mat(c.struct.material);
    return { m_proj: mb, V: Math.max(5, c.aero.Vmo_ms || c.flight.V_ms), material: c.struct.material, t_plate: up.fea?.t_skin_root_m ? Math.max(0.0008, 0.3 * up.fea.t_skin_root_m) : Math.min(c.struct.t_skin_mm, Math.max(0.5, 0.6 * (c.fuselage.dia_m || 1))) / 1e3, a_plate: 6 * R, cs_D: mm.rho < 3500 ? 6500 : 1e5 };
  },
  run(i, ctx) {
    const r = impactSim(i, ctx), fe = r.fe, E = fe.energy, rec = r.rec, m = mat(i.material), warnings = [RES_NOTE], soft = Number.isFinite(r.pc), vEnd = -r.imp.v[2], KE0 = 0.5 * i.m_proj * i.V ** 2, D = 2 * r.Rs;
    const perf = r.holes > 0, passed = !soft && r.imp.c[2] < -r.Rs, vRes = perf || passed ? Math.max(0, vEnd) : 0, reb = vEnd < 0 ? -vEnd : 0, dp = i.m_proj * (i.V - vEnd);
    // analytical screening by the membrane-energy model of the bird-strike analysis
    const scr = N.kv(bird.run({ m_bird: i.m_proj, rho_bird: (PROJ[i.projectile] || PROJ['hard fragment (steel)']).rho, LD: 1, V: i.V, angle_deg: 90, c0: 1482, k_h: 2, t_skin: i.t_plate, material: i.material, k_area: 3, eps_f: i.eps_fail, eta: 0.5 }));
    if (r.sc > 1) warnings.push(`Resolving the true ${(i.t_plate * 1e3).toFixed(2)} mm thickness would need about ${Math.round(r.need)} time steps; one layer ${r.sc.toFixed(1)}× thicker with modulus, strength and density divided by ${r.sc.toFixed(1)} is used instead. Membrane stiffness, strength and mass per area are preserved, but bending and through-thickness (plugging) response are not: raise the step budget for a thick or bending-dominated panel.`);
    if (r.h > r.Rs / 1.5) warnings.push('The in-plane element size exceeds two thirds of the projectile radius: the node-based contact is too coarse; refine the mesh or shrink the panel.');
    if (i.a_plate < 2.99 * D) warnings.push('The panel is smaller than three projectile diameters: the clamped edges stiffen the response.');
    if (soft) warnings.push('Soft-body equivalent: the projectile is a sphere whose contact pressure is capped (stagnation pressure for a bird, crushing strength for ice). It transfers a realistic impulse over a realistic footprint but has no spreading flow or initial shock (Hugoniot) peak; use an SPH or Eulerian bird for substantiation.');
    if (i.material && jcOf(m).fitted && i.rate_model === 'Johnson–Cook') warnings.push('No Johnson–Cook constants are tabulated for this material: a simple fit through yield and ultimate strength with C = 0.01 was used.');
    if (E.hourglass > 0.5 * Math.max(E.internal, 1e-300)) warnings.push(`The hourglass stabilisation carries ${(100 * E.hourglass / (E.hourglass + E.internal)).toFixed(0)}% of the deformation energy. In thin one-point elements it stands in for plate bending and for the local dent under a contact only a few nodes wide (it is scaled to the plate-bending stiffness and capped at the plastic moment), so the dent depth is indicative: refine the in-plane mesh and use four or more elements through the thickness for a resolved answer.`);
    if (r.errMax > 0.05) warnings.push(`The energy balance error reached ${(100 * r.errMax).toFixed(1)}%.`);
    if (!perf && !passed && !(r.tSep > 0) && !soft) warnings.push('The projectile is still in contact at the end of the run: extend the simulated time.');
    // deformed section in the symmetry plane y = 0, coloured by plastic strain
    const nodes = [], tris = [], values = [], ne = r.ne, nz = r.nz, nid = (a, k) => a * (nz + 1) + k, cnt = new Array((ne + 1) * (nz + 1)).fill(0), val = new Array((ne + 1) * (nz + 1)).fill(0);
    for (let a = 0; a <= ne; a++) for (let k = 0; k <= nz; k++) { const nd = r.blk.id(a, 0, k); nodes.push([fe.x[3 * nd] * 1e3, fe.x[3 * nd + 2] * 1e3]); }
    for (let a = 0; a < ne; a++) for (let k = 0; k < nz; k++) { const e = r.col(a, 0, k); if (!fe.hexAlive[e]) continue; const q = [nid(a, k), nid(a + 1, k), nid(a + 1, k + 1), nid(a, k + 1)]; tris.push([q[0], q[1], q[2]], [q[0], q[2], q[3]]); for (const j of q) { val[j] += fe.ep[e]; cnt[j]++; } }
    for (let j = 0; j < val.length; j++) values.push(cnt[j] ? val[j] / cnt[j] : 0);
    const arc = N.linspace(-Math.PI / 2, 0, 24), ms = rec.t.map((t) => t * 1e3), xs = N.range(ne, (a) => (a + 0.5) * r.h * 1e3), map = N.range(ne, (b) => N.range(ne, (a) => { let v = 0; for (let k = 0; k < nz; k++) v = Math.max(v, fe.ep[r.col(a, b, k)]); return v; }));
    return {
      kpis: [
        kpi('impact_perforated', 'Panel perforated (1 = yes)', perf || passed ? 1 : 0, '-', perf || passed ? 'bad' : 'ok', 'Erosion through the full thickness at the stated failure strain'),
        kpi('impact_v_residual_ms', 'Residual projectile velocity', vRes, 'm/s'), kpi('impact_v_rebound_ms', 'Rebound velocity', reb, 'm/s'),
        kpi('impact_force_N', 'Peak contact force', r.Fpk, 'N'), kpi('impact_impulse_Ns', 'Impulse delivered to the panel', r.impulse, 'N·s', undefined, `Projectile momentum ${(i.m_proj * i.V).toFixed(2)} N·s`),
        kpi('impact_defl_m', 'Peak deflection at the impact point', r.defl, 'm'), kpi('impact_eps_p_max', 'Peak equivalent plastic strain', r.epMax, '-', r.epMax < 0.7 * i.eps_fail ? 'ok' : 'warn', `Failure strain ${i.eps_fail}`),
        kpi('impact_eroded', 'Eroded elements (quarter model)', fe.erodedHex, '-'), kpi('impact_plastic_frac', 'Share of the projectile energy dissipated in the panel', E.plastic / KE0, '-'),
        kpi('impact_dia_m', 'Projectile diameter', D, 'm'), kpi('impact_ke_J', 'Projectile kinetic energy', KE0, 'J'),
        kpi('impact_screen_RF', 'Energy-screening reserve factor (analytical)', scr.penetration_RF, '-', scr.penetration_RF >= 1 ? 'ok' : 'warn', 'Membrane-energy screening of the bird-strike analysis with its default factors; below 1 predicts rupture'),
        kpi('impact_screen_agree', 'Screening agrees with the simulation (1 = yes)', (scr.penetration_RF < 1) === (perf || passed) ? 1 : 0, '-'),
        kpi('impact_thickness_scale', 'Thickness scale of the equivalent layer', r.sc, '-', r.sc > 1 ? 'warn' : 'ok', '1 = true thickness resolved'),
        kpi('impact_hourglass_frac', 'Share of deformation energy carried by the hourglass stabilisation', E.hourglass / Math.max(E.hourglass + E.internal, 1e-300), '-', E.hourglass < 0.5 * E.internal ? 'ok' : 'warn', 'Stands in for bending of thin one-point elements; small values mean the bending is resolved by the mesh'),
        kpi('impact_energy_err', 'Energy balance error', r.errMax, '-', r.errMax < 0.03 ? 'ok' : 'warn'), kpi('impact_elements', 'Solid elements (quarter model)', fe.nh, '-'), kpi('impact_steps', 'Time steps', fe.nStep, '-'),
      ],
      plots: [
        { type: 'tri', title: `Deformed panel section in the symmetry plane at ${(fe.t * 1e3).toFixed(2)} ms`, xlabel: 'Distance from the impact point [mm]', ylabel: 'Height [mm]', zlabel: 'Equivalent plastic strain [-]', nodes, tris, values, equalAspect: r.defl > 0.05 * i.a_plate, edges: true, overlay: [{ name: 'Projectile', x: arc.map((a) => r.Rs * Math.cos(a) * 1e3), y: arc.map((a) => (r.imp.c[2] + r.Rs * Math.sin(a)) * 1e3) }] },
        { type: 'heat', title: 'Plastic strain over the quarter panel', xlabel: 'x from the impact point [mm]', ylabel: 'y from the impact point [mm]', zlabel: 'Equivalent plastic strain [-]', x: xs, y: xs, z: map, equalAspect: true },
        { type: 'line', title: 'Projectile velocity', xlabel: 'Time [ms]', ylabel: 'Velocity towards the panel [m/s]', series: [{ name: 'Projectile', x: ms, y: rec.data[0] }] },
        { type: 'line', title: 'Contact force', xlabel: 'Time [ms]', ylabel: 'Force [kN]', series: [{ name: 'Simulation', x: ms, y: rec.data[1].map((f) => f / 1e3) }], annotations: soft ? [{ y: scr.bird_force_avg_N / 1e3, label: 'Hydrodynamic average' }] : [] },
        { type: 'line', title: 'Energy balance (quarter model)', xlabel: 'Time [ms]', ylabel: 'Energy [J]', series: [['Kinetic', 3], ['Internal', 4], ['Plastic work', 5], ['Hourglass', 6], ['Contact', 7], ['Total (should stay constant)', 8]].map(([name, k]) => ({ name, x: ms, y: rec.data[k], ...(k === 8 ? { style: 'dash' } : {}) })) },
      ],
      outputs: { impact_momentum_change_Ns: dp },
      warnings,
      models: [`Explicit finite elements: ${fe.nh} one-point hexahedra with Flanagan–Belytschko hourglass control (quarter model, ${r.nz} through the thickness), ${fe.nStep} steps`, `J2 plasticity with ${i.rate_model} flow stress${i.rate_model === 'Johnson–Cook' ? ' and adiabatic heating' : ''}, erosion at an equivalent plastic strain of ${i.eps_fail}`, soft ? 'Rigid sphere with capped contact pressure as a soft-body equivalent' : 'Rigid sphere with penalty contact', 'Membrane-energy penetration screening for comparison'],
      assumptions: ['Normal impact at the panel centre; two symmetry planes', 'Node-to-sphere penalty contact without friction', 'Erosion at a single failure strain: no stress-triaxiality or shear-band criterion, and the result depends on the element size', 'With one or two one-point elements through the thickness, plate bending is carried largely by the hourglass stabilisation (scaled to the plate-bending stiffness, capped at the plastic moment)', 'The projectile does not deform; soft bodies are represented only by the pressure cap', DATA_NOTE],
    };
  },
  convergence: { param: 'nEl', label: 'Elements along the half side', levels: [4, 6, 8, 10], metric: 'impact_defl_m' },
  calibration: { params: [{ key: 'eps_fail', min: 0.02, max: 1 }], sweep: 'V', target: 'impact_v_residual_ms', note: 'Residual velocities from gas-gun penetration tests at several impact speeds calibrate the erosion strain for the element size in use.' },
  verify() {
    // free elastic panel struck by a rigid sphere: momentum and energy are conserved
    const b = { ...IMPACT_BASE, clamped: false, V: 20, t_factor: 4 }, r = impactSim(b, null, { free: true, elastic: true, noScale: true, maxSteps: 3000 }), p = r.fe.momentum()[2], E = r.fe.energy;
    const pl = impactSim({ ...IMPACT_BASE, V: 120 }, null, { noScale: true, maxSteps: 3000 });
    return [
      N.check('Momentum of sphere + free panel is conserved', p, (-b.m_proj / 4) * b.V, 1e-9, 'Newton’s third law in the penalty contact'),
      N.check('Impulse on the panel = momentum lost by the sphere', r.impulse / 4, (b.m_proj / 4) * (b.V + r.imp.v[2]), 0.01, 'Impulse–momentum theorem (time-centred sum of the contact force)'),
      N.check('Energy balance, elastic impact', r.errMax, 0, 0.02, 'Kinetic + internal + hourglass + contact = initial kinetic energy'),
      N.check('Energy balance with Johnson–Cook plasticity and erosion', pl.errMax, 0, 0.03, 'Same balance; eroded elements keep their internal energy'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.impact_perforated) out.push({ severity: 'warn', title: 'The panel is perforated', detail: `Residual velocity ${o.impact_v_residual_ms.toFixed(0)} m/s; analytical reserve factor ${o.impact_screen_RF.toFixed(2)}.`, action: 'Protect what lies behind (spar, systems, fuel, occupants) for the residual projectile, or raise the panel thickness, use a tougher alloy or a fibre-metal laminate, or add curvature and sweep to make the impact oblique.', basis: 'Erosion through the thickness at the calibrated failure strain' });
    else out.push({ severity: 'info', title: 'The panel contains the projectile', detail: `Dent ${(o.impact_defl_m * 1e3).toFixed(1)} mm, peak plastic strain ${o.impact_eps_p_max.toFixed(3)} against ${i.eps_fail}.`, action: o.impact_eps_p_max > 0.7 * i.eps_fail ? 'The margin on strain is small: repeat with a finer mesh and a calibrated failure strain before taking credit.' : 'Check the dent against the aerodynamic smoothness and residual-strength limits, and the supporting ribs and frames for the contact force.', basis: 'Plastic strain against the failure strain' });
    if (!o.impact_screen_agree) out.push({ severity: 'advise', title: 'Simulation and analytical screening disagree', detail: `Screening reserve factor ${o.impact_screen_RF.toFixed(2)}.`, action: 'The screening uses uncalibrated default factors and the simulation a mesh-dependent erosion strain: a penetration test near this speed is needed to settle it.', basis: 'Model-form uncertainty' });
    return out;
  },
};

// ---- verification of the explicit kernel against exact solutions ----------------------------
function kernelChecks() {
  const out = [];
  { // elastic wave in a bar of hexahedra striking a rigid wall
    const Em = 70e9, rho = 2700, v0 = 5, n = 80, c = Math.sqrt(Em / rho), fe = createFE(), b = hexBlock(fe, n, 1, 1, 1, 1 / n, 1 / n, fe.material({ E: Em, nu: 0, rho, A: Infinity }));
    for (let a = 0; a <= n; a++) for (let j = 0; j < 2; j++) for (let k = 0; k < 2; k++) { if (a === 0) fe.fix(b.id(a, j, k), [1, 0, 0]); else fe.velocity(b.id(a, j, k), [-v0, 0, 0]); }
    fe.init(); const tA = [NaN, NaN], st = [n / 4, (3 * n) / 4], prev = [0, 0], lvl = 0.5 * rho * c * v0;
    fe.run(0.9 / c, () => st.forEach((e, q) => { const s = -fe.sig[6 * b.els[e]]; if (!(tA[q] >= 0) && s >= lvl) tA[q] = fe.t - (fe.dtLast * (s - lvl)) / (s - prev[q]); prev[q] = s; }));
    out.push(N.check('Elastic wave speed in a bar of hexahedra, √(E/ρ)', 0.5 / (tA[1] - tA[0]), c, 0.01, 'D’Alembert solution; half-amplitude arrival at two stations (Poisson ratio 0)'),
      N.check('Stress behind the wave front, σ = ρ·c·v', -fe.sig[6 * b.els[4]], rho * c * v0, 5e-3, 'Momentum jump across the wave front'),
      N.check('Energy balance during wave propagation', fe.energy.error, 0, 0.01, 'Kinetic + internal + hourglass = initial kinetic energy'));
  }
  { // patch test on a distorted mesh
    const Em = 70e9, nu = 0.3, fe = createFE(); let s = 3; const rnd = () => { s = (s * 16807) % 2147483647; return s / 2147483647 - 0.5; }, b = hexBlock(fe, 3, 3, 3, 1, 1, 1, fe.material({ E: Em, nu, rho: 2700, A: Infinity }), (p) => p.map((q) => (q > 1e-9 && q < 1 - 1e-9 ? q + 0.1 * rnd() : q)));
    const Lg = [[1e-3, 2e-4, -3e-4], [5e-4, -2e-3, 1e-4], [0, 4e-4, 1.5e-3]], lam = (Em * nu) / ((1 + nu) * (1 - 2 * nu)), mu = Em / (2 * (1 + nu)), sxx = lam * (Lg[0][0] + Lg[1][1] + Lg[2][2]) + 2 * mu * Lg[0][0];
    fe.init(); for (let q = 0; q < fe.nn; q++) for (let k = 0; k < 3; k++) fe.v[3 * q + k] = Lg[k][0] * fe.x[3 * q] + Lg[k][1] * fe.x[3 * q + 1] + Lg[k][2] * fe.x[3 * q + 2];
    const f = fe.internalForces(1); let es = 0, fi = 0, fb = 0;
    for (let e = 0; e < fe.nh; e++) es = Math.max(es, Math.abs(fe.sig[6 * e] / sxx - 1));
    for (let q = 0; q < f.length; q++) fb = Math.max(fb, Math.abs(f[q]));
    for (let a = 1; a < 3; a++) for (let j = 1; j < 3; j++) for (let k = 1; k < 3; k++) for (let d = 0; d < 3; d++) fi = Math.max(fi, Math.abs(f[3 * b.id(a, j, k) + d]));
    out.push(N.check('Patch test: uniform stress in every distorted hexahedron', 1 + es, 1, 1e-10, 'Constant-strain patch test'), N.check('Patch test: interior nodal forces vanish (relative to boundary forces)', fi / fb, 0, 1e-10, 'Equilibrium of a uniform stress field'));
  }
  { // cantilevers by dynamic relaxation
    const Em = 70e9, L = 2, bw = 0.05, h = 0.1, I = (bw * h ** 3) / 12, A = bw * h, Pt = 1000, n = 8, fe = createFE(), nd = N.range(n + 1, (k) => fe.node((k * L) / n, 0, 0));
    for (let k = 0; k < n; k++) fe.beam(nd[k], nd[k + 1], { E: Em, rho: 2700, A, Iy: I, Iz: (h * bw ** 3) / 12 }, [0, 1, 0]);
    fe.fix(nd[0]); fe.load(nd[n], [0, 0, -Pt]); fe.init(); const w1 = 3.516 * Math.sqrt((Em * I) / (2700 * A * L ** 4)); fe.damping = 2 * w1; fe.run(12 / w1);
    out.push(N.check('Cantilever of beam elements: tip deflection PL³/3EI', -fe.x[3 * nd[n] + 2], (Pt * L ** 3) / (3 * Em * I), 2e-3, 'Euler–Bernoulli beam theory, reached by dynamic relaxation'));
    const L2 = 1, nx = 20, nz = 4, f2 = createFE(), b = hexBlock(f2, nx, 1, nz, L2, bw, h, f2.material({ E: Em, nu: 0, rho: 2700, A: Infinity }));
    for (let j = 0; j < 2; j++) for (let k = 0; k <= nz; k++) { f2.fix(b.id(0, j, k), [1, 1, 1]); f2.load(b.id(nx, j, k), [0, 0, -Pt / (2 * (nz + 1))]); }
    f2.init(); const w2 = 3.516 * Math.sqrt((Em * I) / (2700 * A * L2 ** 4)); f2.damping = 2 * w2; f2.run(10 / w2);
    out.push(N.check('Cantilever of hexahedra (20 × 1 × 4): tip deflection PL³/3EI + 6PL/5GA', h / 2 - f2.x[3 * b.id(nx, 0, nz / 2) + 2], (Pt * L2 ** 3) / (3 * Em * I) + (6 * Pt * L2) / (5 * (Em / 2) * A), 0.03, 'Timoshenko beam theory; one-point elements with four layers are within a few per cent'),
      N.check('Hourglass energy stays small (bending of one-point elements)', f2.energy.hourglass / f2.energy.internal, 0, 0.1, 'Flanagan–Belytschko control: below 10% of the internal energy'));
  }
  { // thin ring pulled across a diameter
    const Em = 70e9, a = 0.04, A = a * a, I = a ** 4 / 12, n = 32, Pt = 100, fe = createFE(), nd = N.range(n, (k) => fe.node(Math.sin((2 * Math.PI * k) / n), 0, -Math.cos((2 * Math.PI * k) / n)));
    for (let k = 0; k < n; k++) { fe.beam(nd[k], nd[(k + 1) % n], { E: Em, rho: 2700, A, Iy: I, Iz: I }, [0, 1, 0]); fe.fix(nd[k], [0, 1, 0, 1, 0, 1]); }
    fe.load(nd[0], [0, 0, -Pt]); fe.load(nd[n / 2], [0, 0, Pt]); fe.init(); const w = 2.683 * Math.sqrt((Em * I) / (2700 * A)); fe.damping = 2 * w; fe.run(14 / w);
    out.push(N.check('Ring of beam elements loaded across a diameter: (π/4 − 2/π)·PR³/EI', fe.x[3 * nd[n / 2] + 2] - fe.x[3 * nd[0] + 2] - 2, ((Math.PI / 4 - 2 / Math.PI) * Pt) / (Em * I), 0.02, 'Castigliano solution for a thin ring; 32 straight elements'));
  }
  { // free fall and bounce of an elastic block
    const h0 = 0.05, fe = createFE({ gravity: [0, 0, -G0] }), b = hexBlock(fe, 2, 2, 2, 0.2, 0.2, 0.2, fe.material({ E: 2e9, nu: 0.3, rho: 1000, A: Infinity }), (p) => [p[0], p[1], p[2] + h0]), pl = fe.plane({});
    for (const k of b.nd) pl.contact(k); fe.init(); let tHit = NaN, vHit = 0, eMax = 0; const top = b.id(1, 1, 2);
    fe.run(0.16, () => { if (!(tHit >= 0) && pl.force > 0) { tHit = fe.t - fe.dtLast; vHit = -fe.v[3 * top + 2]; } eMax = Math.max(eMax, fe.energy.error); });
    out.push(N.check('Free fall: time to impact √(2h/g)', tHit, Math.sqrt((2 * h0) / G0), 2e-3, 'Uniform acceleration'), N.check('Free fall: impact velocity √(2gh)', vHit, Math.sqrt(2 * G0 * h0), 2e-3, 'Energy conservation'),
      N.check('Energy balance through impact and rebound', eMax, 0, 0.03, 'Kinetic + internal + hourglass + contact − gravity work = constant'));
  }
  { // Taylor-type impact of a plastic bar: conservation checks
    const mm = METALS['Al 2024-T3'], jc = mm.JC, fe = createFE(), b = hexBlock(fe, 4, 4, 10, 0.008, 0.008, 0.024, fe.material({ E: mm.E, nu: mm.nu, rho: mm.rho, A: jc.A, B: jc.B, n: jc.n, C: jc.C, m: jc.m, Tm: jc.Tm, cp: mm.cp })), pl = fe.plane({});
    for (const k of b.nd) { pl.contact(k); fe.velocity(k, [0, 0, -150]); } fe.init(); const p0 = fe.momentum()[2]; let imp = 0, eMax = 0;
    fe.run(6e-5, () => { imp += pl.force * fe.dtLast; eMax = Math.max(eMax, fe.energy.error); });
    out.push(N.check('Plastic bar impact: impulse of the wall = change of momentum', imp, fe.momentum()[2] - p0, 0.01, 'Impulse–momentum theorem'), N.check('Plastic bar impact: energy balance with Johnson–Cook plasticity', eMax, 0, 0.02, 'Kinetic + internal (elastic + plastic) + hourglass + contact = initial kinetic energy'),
      N.check('Plastic bar impact: adiabatic heating stores the Taylor–Quinney share of the plastic work', b.els.reduce((q, e) => q + mm.rho * mm.cp * fe.hexVol[e] * (fe.temp[e] - 293), 0), 0.9 * fe.energy.plastic, 0.02, 'Σ ρ·c_p·V·ΔT = β·W_p with β = 0.9 (element volumes change by the elastic strain only)'),
      N.check('Johnson–Cook thermal softening: flow stress halves midway to melting (m = 1)', flowStress(fe.mats[fe.material({ E: 70e9, rho: 2700, A: 300e6, m: 1, Tm: 893, T0: 293 })], 0, 0, 593), 150e6, 1e-12, 'σ = A·(1 − T*^m) with T* = (T − T0)/(Tm − T0)'),
      N.check('Plastic bar impact: hourglass energy fraction', fe.energy.hourglass / fe.energy.internal, 0, 0.05, 'Hourglass energy below 5% of internal energy'), N.check('Plastic bar impact: plastic work is most of the dissipated energy', (fe.energy.plastic + fe.energy.kinetic + fe.energy.viscous) / fe.energy.initial, 1, 0.1, 'Energy partition: remaining share is elastic, hourglass and contact energy'));
  }
  return out;
}
function aircraftChecks() {
  const out = [];
  { // single mass on an elastic–plastic crush spring: crush distance from the work–energy theorem
    const mT = 489, F = 96000, k = 9.6e7, v0 = 6, fe = createFE({ gravity: [0, 0, -G0], dtMax: 1e-5 }), nd = fe.node(0, 0, 0), pl = fe.plane({}); fe.mass(nd, mT); const c = pl.contact(nd, { law: { k, Fy: F, dmax: 5, kb: 0 } }); fe.velocity(nd, [0, 0, -v0]); fe.init(); fe.run(0.06);
    out.push(N.check('Mass on a crush spring: crush distance from the energy balance', pl.pts[c].dMax, (0.5 * mT * v0 * v0 + (F * F) / (2 * k)) / (F - mT * G0), 2e-3, 'Work–energy theorem'), N.check('Mass on a crush spring: energy balance', fe.energy.error, 0, 2e-3, 'Kinetic + contact work − gravity work = initial kinetic energy'));
  }
  { // Coulomb friction: sliding distance v²/(2μg)
    const mS = 100, k = 1e6, mu = 0.4, v0 = 10, fe = createFE({ gravity: [0, 0, -G0], dtMax: 5e-4 }), nd = fe.node(0, 0, -(mS * G0) / k), pl = fe.plane({ mu }); fe.mass(nd, mS); pl.contact(nd, { law: { k, Fy: 1e12, dmax: 1, kb: 0 } }); fe.velocity(nd, [v0, 0, 0]); fe.init(); fe.run(1.2 * v0 / (mu * G0));
    out.push(N.check('Sliding distance under Coulomb friction, v²/(2μg)', fe.x[0], (v0 * v0) / (2 * mu * G0), 2e-3, 'Constant deceleration μg'), N.check('Friction work equals the kinetic energy lost', fe.energy.friction, 0.5 * mS * v0 * v0, 5e-3, 'Work–energy theorem'));
  }
  { // free spinning and translating beam: large rigid rotation without spurious straining
    const om = 5, fe = createFE({}), a = fe.node(-1, 0, 0), b = fe.node(1, 0, 0); fe.beam(a, b, { E: 70e9, rho: 2700, A: 1e-3, Iy: 1e-6, Iz: 1e-6 }, [0, 0, 1]); fe.mass(a, 10); fe.mass(b, 10); fe.velocity(a, [3, -om, 0]); fe.velocity(b, [3, om, 0]); fe.init();
    for (const q of [a, b]) { fe.w[3 * q + 2] = om; fe.energy.initial += 0.5 * fe.inertia[q] * om * om; }
    const Lz = () => { let s = 0; for (const q of [a, b]) s += fe.mass_[q] * (fe.x[3 * q] * fe.v[3 * q + 1] - fe.x[3 * q + 1] * fe.v[3 * q]) + fe.inertia[q] * fe.w[3 * q + 2]; return s; }, L0 = Lz(); fe.run((6 * Math.PI) / om);
    const len = Math.hypot(fe.x[3] - fe.x[0], fe.x[4] - fe.x[1]);
    out.push(N.check('Spinning beam: angular momentum conserved over three revolutions', Lz(), L0, 1e-3, 'Free rigid-body motion: orbital plus nodal spin angular momentum about the origin'),
      N.check('Spinning beam: kinetic energy conserved', fe.energy.kinetic + fe.energy.internal, fe.energy.initial, 1e-3, 'No spurious strain energy from the co-rotational formulation'), N.check('Spinning beam: length preserved (centrifugal stretch only)', len, 2, 2e-5, 'Rigid rotation; the centrifugal stretch is of order 1e-5 m'));
  }
  const r = N.kv(aircraft3d.run(AC_BASE));
  out.push(N.check('Full-aircraft model: energy balance error of the reference crash', r.ac_energy_err, 0, 0.03, 'Global energy balance of the explicit solution'), N.check('Full-aircraft model: generated masses add up to the aircraft mass', r.ac_mass_kg, AC_BASE.mass_kg, 1e-6, 'Mass bookkeeping of the model generator'));
  return out;
}

export default {
  id: 'crash', n: 15,
  tagline: 'Will the occupants survive the impact, and what does a bird, the ground or the water do to the structure?',
  analyses: [pulse, absorber, head, bird, wave, ditch, emergency, barrel3d, aircraft3d, impact3d],
  consumes: [{ from: 'fea', keys: ['t_skin_root_m', 'EI_root_Nm2'], why: 'Skin gauge for the bird-strike and projectile screening; wing bending stiffness of the full-aircraft crash model' }, { from: 'performance', keys: ['V_stall_ms'], why: 'Forward speed of the full-aircraft crash scenario (1.1 × stall speed)' }],
  provides: [
    { key: 'peak_g', label: 'Peak crash deceleration', unit: 'g' }, { key: 'crush_m', label: 'Structural crush', unit: 'm' }, { key: 'absorbed_J', label: 'Energy absorbed', unit: 'J' },
    { key: 'bird_force_N', label: 'Bird strike peak force', unit: 'N' }, { key: 'HIC', label: 'Head Injury Criterion', unit: '-' },
  ],
  handoff: [
    { model: 'Detailed full-airframe crash models with millions of shell and solid elements', why: 'The native explicit kernel solves a frame-bay barrel model, a hybrid mass–beam–spring model of the whole aircraft and small solid-element panels in about a second; detailed shell meshes with self-contact, rivet and joint failure and local buckling need 10⁵–10⁷ elements and hours of computing. The High-fidelity bridge page exports a ready-to-run case for them', tool: 'Explicit FE (LS-DYNA / Radioss / Abaqus Explicit / PAM-CRASH class)' },
    { model: 'Anthropomorphic dummy, restraint and airbag models', why: 'Multibody or FE dummies with validated joints, belts, airbags and contact are required for certification injury metrics; a DRI spine on a stroking seat and a head form are used here. The High-fidelity bridge page exports the seat pulse and model set-up', tool: 'MADYMO / LS-DYNA dummy models and dynamic seat tests' },
    { model: 'Smoothed particle hydrodynamics and ALE bird, hail and water models', why: 'A fluid-like projectile with large deformation coupled to a deforming target; natively the bird is hydrodynamic theory or a pressure-capped rigid sphere, and water is a pressure surface. The High-fidelity bridge page exports the impact case', tool: 'Explicit FE with SPH/ALE/CEL' },
    { model: 'Composite progressive crush and delamination under impact', why: 'Crush stress of composite absorbers is test-derived and mesh-sensitive in simulation; the native material model is isotropic metal plasticity', tool: 'Explicit FE with composite damage models; component crush tests' },
    { model: 'Calibrated fracture, fragmentation and blade-off containment', why: 'Native erosion uses a single equivalent-plastic-strain limit on a coarse mesh; triaxiality-dependent failure, fragment clouds and rotating-blade containment need calibrated criteria and fine 3-D meshes', tool: 'Explicit FE with Johnson–Cook damage or equivalent; rig tests' },
    { model: 'Fuel-tank slosh, rupture and post-crash fire', why: 'Coupled fluid–structure impact outside the scope of the native models', tool: 'Explicit FSI simulation; tank drop tests' },
    { model: 'Terrain, obstacle and soft-soil interaction, break-up and tumbling', why: 'The full-aircraft model strikes a flat rigid or water-like surface with friction; ploughing, obstacles and post-break-up dynamics are not modelled', tool: 'Full-vehicle explicit FE crash model with soil and obstacle models' },
  ],
};
