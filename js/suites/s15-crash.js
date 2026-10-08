// Suite 15 — Crashworthiness and Impact Mechanics.
// Lumped-mass crash and seat-pulse model with crush zone, stroking seat and spinal (DRI) response; crush-tube
// energy absorbers; head strike and HIC; bird strike by hydrodynamic theory; explicit 1-D elastic–plastic wave
// propagation with Johnson–Cook / Cowper–Symonds rate effects; water impact of a wedge; emergency-landing inertia loads.

import * as N from '../core/numerics.js';
import { G0 } from '../core/atmosphere.js';
import { METALS } from '../data/materials.js';

// ---- shared helpers -------------------------------------------------------------------------
const MATS = Object.keys(METALS);
const mat = (name) => METALS[name] || METALS['Al 2024-T3'];
const MAT = { key: 'material', label: 'Material', type: 'select', options: MATS, default: 'Al 2024-T3', group: 'Material', help: 'Typical handbook values, not design allowables' };
const kpi = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const thin = (a, max = 300) => { if (a.length <= max) return a.slice(); const o = []; for (let k = 0; k < max; k++) o.push(a[Math.round((k * (a.length - 1)) / (max - 1))]); return o; };
const DATA_NOTE = 'Material data are typical handbook values, not design allowables';
const DRI_W = 52.9, DRI_Z = 0.224; // Dynamic Response Index spinal model: natural frequency [rad/s] and damping ratio
const occupants = (c) => (c.mission.pax > 0 ? c.mission.pax + (c.mission.pax > 19 ? 2 : 1) : 0);
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
const PULSE_BASE = { mode: 'vertical impact', v0: 9.1, m_air: 400, m_seat: 55, m_torso: 34, F_crush: 96000, d_crush: 0.3, k_crush: 9.6e6, F_ea: 12700, s_ea: 0.3, k_seat: 2.5e6, G_pulse: 14, t_rise: 0.08, shape: 'triangular', t_end: 0.3, occupied: true, nSteps: 4000 };
const pulse = {
  id: 'pulse', title: 'Crash pulse, seat stroke and spinal load (lumped-mass model)', fidelity: 'numerical',
  summary: 'A chain of airframe, seat and upper-body masses hits the ground through a crushable subfloor, or is driven by a seat-test floor pulse. Explicit time integration gives the crash pulse, structural crush, energy absorbed, seat stroke, the Dynamic Response Index and a lumbar-load estimate.',
  equations: ['Conservation of linear momentum', 'Conservation of energy', 'Transient structural dynamics equations', 'Crushing energy absorption equations', 'Occupant dynamics equations', 'Plasticity constitutive equations'],
  inputs: [
    { key: 'mode', label: 'Scenario', type: 'select', options: ['vertical impact', 'seat test pulse'], default: 'vertical impact', group: 'Scenario', help: 'Vertical impact: the whole chain falls onto the ground. Seat test pulse: the floor follows a prescribed deceleration pulse' },
    { key: 'v0', label: 'Vertical impact velocity', unit: 'm/s', default: 9.1, min: 0.1, max: 40, group: 'Scenario', help: '9.1 m/s (30 ft/s) is the civil seat-test vertical velocity change; military rotorcraft design to 12.8 m/s (42 ft/s)' },
    { key: 'G_pulse', label: 'Seat-test pulse peak', unit: 'g', default: 14, min: 1, max: 100, group: 'Scenario', help: 'Transport seat tests use a 14 g mainly-vertical and a 16 g longitudinal triangular pulse; confirm against your certification basis' },
    { key: 't_rise', label: 'Seat-test pulse rise time', unit: 's', default: 0.08, min: 0.002, max: 1, group: 'Scenario' },
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
    const mSeat = has ? 55 : 0.6 * mOcc, mTor = has ? 34 : 0.4 * mOcc, mTot = mAir + mSeat + mTor, rot = c.meta.type === 'helicopter' || c.meta.type === 'evtol', d = N.clamp(0.12 * (c.fuselage.dia_m || 1), 0.02, 0.5) + (c.gear.stroke_m || 0), v0 = has ? 9.1 : Math.max(2, 2 * c.gear.sink_ms); // subfloor depth plus landing-gear stroke
    // plateau force set so the crush zone can absorb the impact energy in about 80% of the available depth, but not below 8 g
    const Fc = mTot * Math.max(8 * G0, (v0 * v0) / (2 * 0.8 * d) + G0), Fea = (has ? 12 : 30) * G0 * (mSeat + mTor);
    return { v0, m_air: mAir, m_seat: mSeat, m_torso: mTor, occupied: has, F_crush: Fc, d_crush: d, k_crush: Fc / Math.max(0.005, 0.05 * d), F_ea: Fea, s_ea: !has ? 0.02 : rot ? 0.3 : 0.15, k_seat: Math.max(2.5e6 * (mSeat + mTor) / 89, 100) };
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
    const absorbed = s1.wp + s2.wp, Ein = sled ? NaN : E0, dri = (DRI_W ** 2 * dSp) / G0, bottomC = crush > i.d_crush, bottomS = stroke > i.s_ea && i.s_ea > 0;
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
        kpi('seat_stroke_m', 'Seat stroke used', stroke, 'm', bottomS ? 'bad' : 'ok', `Available ${i.s_ea.toFixed(3)} m`),
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
      assumptions: ['Purely vertical, one-dimensional impact on a rigid surface; no pitch, roll or horizontal velocity', 'Constant crush plateau force; bottoming modelled as a 20× stiffer stop', 'The occupant is a single upper-body mass on the DRI spring–damper; no restraint, flail or head strike', 'Injury thresholds are indicative; the lumbar load of a lumped model is not a substitute for a dummy test'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [2000, 4000, 8000, 16000, 32000], metric: 'DRI' },
  verify() {
    // (a) step floor pulse on a stiff seat: DRI = G·(1 + exp(−ζπ/√(1−ζ²)))
    const st = N.kv(pulse.run({ ...PULSE_BASE, mode: 'seat test pulse', shape: 'step', G_pulse: 10, F_ea: 1e9, k_seat: 2e9, m_seat: 1000, t_end: 0.25, nSteps: 20000 }));
    // (b) rigid chain on a plateau: m·v²/2 + m·g·δ = F·(δ − F/2k)
    const mT = 489, F = 96000, k = 9.6e7, v0 = 6, im = N.kv(pulse.run({ ...PULSE_BASE, v0, m_air: 488.99, m_seat: 0.005, m_torso: 0.005, k_crush: k, d_crush: 5, F_ea: 1e9, k_seat: 5e4, t_end: 0.25, nSteps: 20000 })), dExact = (0.5 * mT * v0 * v0 + (F * F) / (2 * k)) / (F - mT * G0);
    return [
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
    { key: 'v0', label: 'Impact velocity', unit: 'm/s', default: 9.1, min: 0.1, max: 60, group: 'Requirement' },
    { key: 'G_lim', label: 'Deceleration limit', unit: 'g', default: 20, min: 1.5, max: 200, group: 'Requirement' },
    { key: 'shape', label: 'Tube section', type: 'select', options: ['circular', 'square'], default: 'circular', group: 'Tube' },
    { key: 'D', label: 'Diameter or side length', unit: 'm', default: 0.06, min: 0.003, group: 'Tube' },
    { key: 't', label: 'Wall thickness', unit: 'm', default: 0.0015, min: 1e-4, group: 'Tube' },
    { key: 'L', label: 'Tube length', unit: 'm', default: 0.3, min: 0.01, group: 'Tube' },
    { key: 'n_tubes', label: 'Number of tubes', unit: '', default: 4, min: 1, max: 1000, step: 1, discrete: true, group: 'Tube' },
    { key: 'eff', label: 'Stroke efficiency', unit: '-', default: 0.72, min: 0.3, max: 0.9, group: 'Tube', help: 'Usable crush / length before the folds lock up; about 0.7–0.75 for metal tubes (empirical)' },
    MAT,
    { key: 'method', label: 'Mean-force formula', type: 'select', options: ['Abramowicz–Jones', 'Alexander / Wierzbicki–Abramowicz'], default: 'Abramowicz–Jones', group: 'Method' },
    { key: 'cs_D', label: 'Cowper–Symonds D', unit: '1/s', default: 6500, min: 1, group: 'Method', help: 'About 6500 1/s for aluminium alloys and 40.4 1/s for mild steel; high-strength steels and titanium are far less rate-sensitive' },
    { key: 'cs_q', label: 'Cowper–Symonds q', unit: '-', default: 4, min: 1, max: 20, group: 'Method', help: '4 for aluminium, 5 for mild steel' },
  ],
  defaults(c) {
    const occ = occupants(c), m = c.mass.mtow_kg, v0 = occ > 0 ? 9.1 : Math.max(2, 2 * c.gear.sink_ms), mm = mat(c.struct.material), n = 4, G = 20, s0 = 0.5 * (mm.Sy + mm.Su);
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
      assumptions: ['Progressive axial folding at a constant mean force; flow stress = mean of yield and ultimate', 'Vertical impact: the weight continues to act during the stroke', 'Stroke efficiency and strain-rate constants are empirical inputs', 'Composite or honeycomb absorbers need their own test-derived crush stress', DATA_NOTE],
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
      assumptions: ['Rigid target for the load calculation; target compliance lowers the real peak force', 'Only the normal momentum is transferred; the tangential part slides off', 'Bird geometry: cylinder with hemispherical ends', DATA_NOTE],
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
    if (o.penetration_RF < 1) out.push({ severity: 'warn', title: 'Skin penetration is likely', detail: `Reserve factor ${o.penetration_RF.toFixed(2)}; about ${(o.t_required_m * 1e3).toFixed(1)} mm of this alloy would be needed.`, action: 'Rather than thickening the whole skin, add a splitter or deflector, use a tougher leading-edge material (fibre-metal laminate), or protect the structure and systems behind the skin.', basis: 'Empirical membrane-energy screening' });
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
        kpi('sigma_elastic_Pa', 'Elastic impact stress ρ·c·v', sEl, 'Pa', sEl < m.Sy ? 'ok' : 'warn', sEl < m.Sy ? 'Below yield' : 'Exceeds yield: plastic wave forms'),
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
      assumptions: ['Uniaxial stress, small-strain kinematics, no radial inertia', 'Strain rate in the flow stress taken as the total element strain rate', 'No fracture or erosion criterion', DATA_NOTE],
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

export default {
  id: 'crash', n: 15,
  tagline: 'Will the occupants survive the impact, and what does a bird, the ground or the water do to the structure?',
  analyses: [pulse, absorber, head, bird, wave, ditch, emergency],
  consumes: [{ from: 'fea', keys: ['t_skin_root_m'], why: 'Skin gauge for the bird-strike penetration screening' }],
  provides: [
    { key: 'peak_g', label: 'Peak crash deceleration', unit: 'g' }, { key: 'crush_m', label: 'Structural crush', unit: 'm' }, { key: 'absorbed_J', label: 'Energy absorbed', unit: 'J' },
    { key: 'bird_force_N', label: 'Bird strike peak force', unit: 'N' }, { key: 'HIC', label: 'Head Injury Criterion', unit: '-' },
  ],
  handoff: [
    { model: 'Explicit nonlinear finite-element crash simulation of the airframe', why: 'Needs detailed shell/solid meshes, contact, material failure and millions of small time steps; here the structure is a lumped crush element', tool: 'Explicit FE (LS-DYNA / Radioss / Abaqus Explicit / PAM-CRASH class)' },
    { model: 'Smoothed particle hydrodynamics and ALE bird, hail and water models', why: 'Fluid-like projectile with large deformation coupled to a deforming target; only hydrodynamic theory and momentum pulses are computed', tool: 'Explicit FE with SPH/ALE/CEL' },
    { model: 'Anthropomorphic dummy and restraint (belt, airbag) models', why: 'Multibody or FE dummies with validated joints and contact are required for certification injury metrics; a DRI spine and a head form are used here', tool: 'MADYMO / LS-DYNA dummy models and dynamic seat tests' },
    { model: 'Composite progressive crush and delamination under impact', why: 'Crush stress of composite absorbers is test-derived and mesh-sensitive in simulation', tool: 'Explicit FE with composite damage models; component crush tests' },
    { model: 'Fracture, element erosion, fragmentation and blade-off containment', why: 'Requires calibrated failure criteria and 3-D explicit analysis; the bar model has no fracture', tool: 'Explicit FE with Johnson–Cook damage or equivalent; rig tests' },
    { model: 'Fuel-tank slosh, rupture and post-crash fire', why: 'Coupled fluid–structure impact outside the scope of lumped models', tool: 'Explicit FSI simulation; tank drop tests' },
    { model: 'Combined-velocity, pitched or rolled impacts and terrain interaction', why: 'The lumped model is purely vertical', tool: 'Full-vehicle multibody or explicit FE crash model' },
  ],
};
