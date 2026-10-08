// Suite 14 — Landing Gear, Ground Dynamics and Impact.
// Oleo-pneumatic drop simulation, braking with wheel-slip dynamics and anti-skid, ground loads and stability,
// taxi over rough runways, nose-wheel shimmy, and skid-gear landings for rotorcraft and multirotors.

import * as N from '../core/numerics.js';
import { isa, G0, RHO0 } from '../core/atmosphere.js';
import { METALS, FLUIDS } from '../data/materials.js';

// ---- shared helpers --------------------------------------------------------------------------
const num = (key, label, unit, def, min, max, group, help, x) => ({ key, label, unit, default: def, min, max, group, ...(help ? { help } : {}), ...x });
const sel = (key, label, options, def, group, help) => ({ key, label, type: 'select', options, default: def, group, ...(help ? { help } : {}) });
const kp = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const thin = (a, n = 300) => { if (a.length <= n) return a; const s = (a.length - 1) / (n - 1); return N.range(n, (j) => a[Math.round(j * s)]); };
const wheeled = (c) => (c.gear.type !== 'skid' ? true : 'This aircraft has skid gear: use the skid-landing and ground-stability analyses.');
const landMass = (c) => Math.max(0.5 * c.mass.mtow_kg, c.mass.mtow_kg - 0.8 * c.mass.fuel_kg);
/** Touchdown ground speed estimate [m/s]: 1.15·VS in landing configuration, or a run-on landing for rotorcraft. */
const vTouch = (c, up) => (c.wing.S_m2 > 0 ? 1.15 * (up?.performance?.V_stall_ms && c.aero.CLmax_clean > 0 ? up.performance.V_stall_ms * Math.sqrt((c.aero.CLmax_clean / c.aero.CLmax_land) * (landMass(c) / c.mass.mtow_kg)) : Math.sqrt((2 * landMass(c) * G0) / (RHO0 * c.wing.S_m2 * c.aero.CLmax_land))) : 18);
const headwind = (c) => { const s = c.site || {}; return (s.wind_ms || 0) * Math.cos(N.rad((s.wind_dir_deg || 0) - (s.runway_heading_deg ?? s.wind_dir_deg ?? 0))); };
const OIL_RHO = FLUIDS['MIL-PRF-5606 hydraulic'].rho;

// ---- oleo-pneumatic strut model shared by the drop and taxi analyses ---------------------------
const STRUT_INPUTS = [
  num('m_strut', 'Mass carried by this strut', 'kg', 33000, 0.01, 5e5, 'Aircraft', 'Landing mass divided by the number of main struts (two-point landing)'),
  num('unsprung_frac', 'Unsprung mass / mass on the strut', '-', 0.015, 0.001, 0.3, 'Aircraft', 'Wheels, tyres, brakes and axle'),
  num('stroke', 'Strut stroke available', 'm', 0.45, 0.005, 2, 'Shock strut'),
  num('p_static', 'Strut air pressure at static load', 'Pa', 10e6, 5e4, 3e7, 'Shock strut', 'About 10 MPa (1500 psi) on transport oleos, lower on light aircraft'),
  num('r_s', 'Compression ratio static / extended', '-', 4, 1.2, 12, 'Shock strut', 'Typical 4:1 (isothermal)'), num('r_c', 'Compression ratio compressed / static', '-', 3, 1.2, 10, 'Shock strut', 'Typical 3:1 (isothermal)'),
  num('n_poly', 'Polytropic exponent during the stroke', '-', 1.3, 1, 1.4, 'Shock strut', '1.35 gas separated from oil, about 1.1 when gas and oil mix'),
  num('A_o', 'Net orifice area (0 = size automatically)', 'm²', 0, 0, 0.1, 'Shock strut'), num('Cd', 'Orifice discharge coefficient', '-', 0.8, 0.4, 1, 'Shock strut'),
  num('k_reb', 'Recoil damping / compression damping', '-', 4, 1, 50, 'Shock strut', 'Recoil (snubber) orifice is smaller to stop the aircraft bouncing'),
  num('k_tyre', 'Tyre vertical stiffness (all tyres on the strut)', 'N/m', 4.4e6, 100, 1e9, 'Tyre'), num('zeta_tyre', 'Tyre damping ratio', '-', 0.03, 0, 0.5, 'Tyre'),
];
function strutGeom(i) {
  const W = i.m_strut * G0, Aa = W / i.p_static, V0 = (Aa * i.stroke) / (1 - 1 / (i.r_c * i.r_s)), p0 = i.p_static / i.r_s, m2 = i.unsprung_frac * i.m_strut, m1 = i.m_strut - m2, kt = i.k_tyre, kS = 30 * kt, n = i.n_poly, sLim = (0.999 * V0) / Aa;
  const Fa = (s) => p0 * Aa * (V0 / (V0 - Aa * N.clamp(s, 0, sLim))) ** n, Ua = (s) => { const x = V0 / (V0 - Aa * N.clamp(s, 0, sLim)); return Math.abs(n - 1) < 1e-9 ? p0 * V0 * Math.log(x) : ((p0 * V0) / (n - 1)) * (x ** (n - 1) - 1); };
  const Ch = (Ao) => (OIL_RHO * Aa ** 3) / (2 * (i.Cd * Ao) ** 2), ct = 2 * i.zeta_tyre * Math.sqrt(kt * i.m_strut), cS = Math.sqrt(kS * m2);
  // total strut force: air spring + orifice damping, with stiff extension and bottoming stops
  const Fs = (s, sd, C) => Fa(s) + C * sd * Math.abs(sd) * (sd < 0 ? i.k_reb : 1) + (s < 0 ? kS * s + cS * sd : s > i.stroke ? kS * (s - i.stroke) : 0);
  return { W, Aa, V0, p0, m1, m2, kt, kS, cS, ct, Fa, Ua, Ch, Fs, sStatic: (V0 * (1 - 1 / i.r_s)) / Aa, AoGuess: (Vs) => Math.sqrt((OIL_RHO * Aa ** 3 * Vs * Vs) / (2 * 1.1 * W)) / i.Cd };
}
/** 2-DOF drop: sprung mass on the strut, unsprung mass on the tyre. z positive down from first tyre contact. */
function dropSim(i, gm, Ao, Vs, nSteps, tEnd = i.t_end) {
  const C = gm.Ch(Ao), L = i.lift_frac * gm.W, { m1, m2, kt, ct } = gm, s0 = -gm.p0 * gm.Aa / gm.kS;
  const tyre = (d, dd) => (d > 0 ? Math.max(0, kt * d + ct * dd) : 0);
  const r = N.rk4((t, y) => { const F = gm.Fs(y[0] - y[2], y[1] - y[3], C); return [y[1], (m1 * G0 - L - F) / m1, y[3], (m2 * G0 + F - tyre(y[2], y[3])) / m2]; }, 0, [s0, Vs, 0, Vs], tEnd, nSteps);
  const s = r.y.map((y) => y[0] - y[2]), sd = r.y.map((y) => y[1] - y[3]), Fst = s.map((v, k) => gm.Fs(v, sd[k], C)), Ft = r.y.map((y) => tyre(y[2], y[3])), d = r.y.map((y) => y[2]);
  const kF = N.argmax(Ft), kS = N.argmax(s), sMax = s[kS], FsMax = N.amax(Fst), a1 = Fst.map((F) => (F + L - m1 * G0) / (m1 * G0));
  // energy audit: initial kinetic + work of (weight − lift) = kinetic + stored + dissipated
  const Ph = s.map((v, k) => C * Math.abs(sd[k]) ** 3 * (sd[k] < 0 ? i.k_reb : 1) + (v < 0 ? gm.cS * sd[k] ** 2 : 0) + (Ft[k] - (d[k] > 0 ? kt * d[k] : 0)) * r.y[k][3]), Ed = N.trapz(r.t, Ph), e = r.y.length - 1, ye = r.y[e];
  const KE = (y) => 0.5 * m1 * y[1] ** 2 + 0.5 * m2 * y[3] ** 2, st = gm.Ua(s[e]) + (d[e] > 0 ? 0.5 * kt * d[e] ** 2 : 0) + (s[e] < 0 ? 0.5 * gm.kS * s[e] ** 2 : s[e] > i.stroke ? 0.5 * gm.kS * (s[e] - i.stroke) ** 2 : 0) - 0.5 * gm.kS * s0 * s0;
  const Wg = (m1 * G0 - L) * (ye[0] - s0) + m2 * G0 * ye[2], eErr = Math.abs(KE(r.y[0]) + Wg - KE(ye) - st - Ed) / KE(r.y[0]);
  const work = N.trapz(s.slice(0, kS + 1), Fst.slice(0, kS + 1));
  return { t: r.t, s, Fst, Ft, d, a1, FtMax: Ft[kF], tPeak: r.t[kF], sMax, FsMax, dMax: N.amax(d), eff: work / Math.max(1e-12, FsMax * sMax), Nr: Ft[kF] / gm.W, eErr, kS, Ehyd: Ed, bottomed: sMax > i.stroke * 1.0005, C };
}
function sizeOrifice(i, gm, Vs, nSteps) {
  if (i.A_o > 0) return i.A_o;
  const f = (lg) => dropSim(i, gm, 10 ** lg, Vs, nSteps).FtMax;
  return 10 ** N.goldenSection(f, Math.log10(gm.Aa * 3e-4), Math.log10(gm.Aa * 0.3), 1e-3);
}
const strutDefaults = (c) => {
  const m = landMass(c) / Math.max(1, c.gear.n_main), small = c.mass.mtow_kg < 150, light = c.mass.mtow_kg <= 5700;
  return { m_strut: m, stroke: c.gear.stroke_m, p_static: small ? 0.8e6 : light ? 4e6 : 10e6, k_tyre: c.gear.tyre_k_Npm * Math.max(1, c.gear.tyres_per_strut), unsprung_frac: light ? 0.03 : 0.015 };
};

// ---- 1. landing impact ------------------------------------------------------------------------
const drop = {
  id: 'drop', title: 'Landing impact: oleo-pneumatic strut drop simulation', fidelity: 'numerical',
  summary: 'Simulates touchdown on one main gear: a polytropic air spring and orifice damper in series with the tyre absorb the sink energy. Gives the ground load factor, stroke used, the load–stroke curve and strut efficiency, and compares them with the classical energy method.',
  equations: ['Newton–Euler equations', 'Oleo-pneumatic shock absorber equations', 'Hydraulic damping equations', 'Landing impact energy equations', 'Tyre force equations'],
  applicable: wheeled,
  inputs: [
    ...STRUT_INPUTS,
    num('sink', 'Sink speed at touchdown', 'm/s', 3.05, 0.05, 15, 'Touchdown', '3.05 m/s (10 ft/s) is the usual transport limit case; 1.8 m/s (6 ft/s) at maximum take-off mass'),
    num('lift_frac', 'Wing or rotor lift / weight during the impact', '-', 1, 0, 1.2, 'Touchdown', '1.0 for transport aeroplanes; 2/3 is customary for rotorcraft and light aircraft drop tests'),
    num('N_design', 'Design ground reaction factor', '-', 1.5, 0.5, 8, 'Limits', 'Typically 1.2–1.5 for large transports, about 2 for commuter aircraft and helicopters, 2–3 for light aircraft and more for small unmanned aircraft (ground load / static load)'),
    num('strut_eff_ref', 'Strut efficiency assumed by the energy method', '-', 0.8, 0.3, 0.95, 'Limits'),
    num('t_end', 'Simulated time', 's', 0.9, 0.02, 10, 'Numerics'), num('nSteps', 'Time steps', '', 4000, 400, 400000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => ({ ...strutDefaults(c), sink: c.gear.sink_ms, lift_frac: c.meta.type === 'helicopter' || c.mass.mtow_kg <= 5700 ? 0.667 : 1, N_design: c.mass.mtow_kg > 30000 ? 1.5 : c.mass.mtow_kg > 5700 ? 2 : c.mass.mtow_kg > 150 ? 3 : 5, strut_eff_ref: c.gear.strut_eff, t_end: N.clamp((7 * c.gear.stroke_m) / Math.max(0.2, c.gear.sink_ms), 0.15, 3) }),
  run(i) {
    const gm = strutGeom(i), n = Math.round(i.nSteps), Ao = sizeOrifice(i, gm, i.sink, Math.max(400, Math.round(n / 2))), r = dropSim(i, gm, Ao, i.sink, n), warnings = [];
    // energy method (Currey): V²/2g + (1 − L)(S + St) = ηs·N·S + ηt·N·St with a tyre efficiency of 0.47
    const St = (r.Nr * gm.W) / gm.kt, kl = 1 - i.lift_frac, Sreq = (i.sink ** 2 / (2 * G0) - (0.47 * r.Nr - kl) * St) / Math.max(1e-9, i.strut_eff_ref * r.Nr - kl), nz = i.lift_frac + N.amax(r.Fst) / gm.W;
    const sinks = [0.5, 0.75, 1, 1.25, 1.5].map((f) => f * i.sink), sw = sinks.map((v) => dropSim(i, gm, Ao, v, n));
    if (r.bottomed) warnings.push('The strut bottoms out: the remaining energy goes into the structure. Increase stroke, static pressure ratio or damping.');
    if (r.eErr > 0.01) warnings.push('Energy balance closure is poorer than 1%: increase the number of time steps.');
    if (r.Nr > i.N_design) warnings.push(`Ground reaction factor ${r.Nr.toFixed(2)} exceeds the design value ${i.N_design.toFixed(2)}: this is a hard-landing exceedance requiring inspection.`);
    if (r.t[N.argmax(r.s)] > 0.95 * i.t_end) warnings.push('Maximum stroke occurs at the end of the simulated time: lengthen the simulation.');
    const k1 = Math.min(r.s.length - 1, r.kS + Math.round(0.02 * r.s.length));
    return {
      kpis: [
        kp('gear_load_N', 'Peak ground load on one main gear', r.FtMax, 'N'), kp('gear_load_factor', 'Ground reaction factor (peak load / static load)', r.Nr, '-', r.Nr <= i.N_design ? 'ok' : 'bad', `Design ${i.N_design.toFixed(2)}`),
        kp('stroke_used_m', 'Strut stroke used', r.sMax, 'm', r.sMax <= 0.95 * i.stroke ? 'ok' : r.bottomed ? 'bad' : 'warn', `${((100 * r.sMax) / i.stroke).toFixed(0)}% of ${i.stroke.toFixed(3)} m`),
        kp('strut_efficiency', 'Strut efficiency', r.eff, '-', r.eff > 0.75 ? 'ok' : 'warn', 'Work absorbed / (peak force × stroke used)'), kp('nz_cg', 'Peak load factor at the mass', nz, 'g'),
        kp('tyre_deflection_m', 'Peak tyre deflection', r.dMax, 'm'), kp('strut_force_max_N', 'Peak strut force', r.FsMax, 'N'), kp('orifice_area_m2', i.A_o > 0 ? 'Orifice area (specified)' : 'Orifice area (sized for minimum peak load)', Ao, 'm²'),
        kp('piston_area_m2', 'Pneumatic piston area', gm.Aa, 'm²'), kp('stroke_static_m', 'Static strut position', gm.sStatic, 'm'), kp('energy_hydraulic_J', 'Energy dissipated', r.Ehyd, 'J'),
        kp('stroke_energy_method_m', 'Stroke by the energy method at the same load factor', Sreq, 'm', undefined, `Assumed strut efficiency ${i.strut_eff_ref}`), kp('t_peak_s', 'Time of peak ground load', r.tPeak, 's'),
        kp('energy_error', 'Energy balance closure error', r.eErr, '-', r.eErr < 0.01 ? 'ok' : 'warn'),
      ],
      plots: [
        { type: 'line', title: 'Loads during the impact', xlabel: 'Time [s]', ylabel: 'Force [kN]', series: [{ name: 'Ground (tyre) load', x: thin(r.t), y: thin(r.Ft).map((v) => v / 1e3) }, { name: 'Strut force', x: thin(r.t), y: thin(r.Fst).map((v) => v / 1e3) }], annotations: [{ y: gm.W / 1e3, label: 'Static load' }] },
        { type: 'line', title: 'Strut load–stroke curve', xlabel: 'Stroke [m]', ylabel: 'Strut force [kN]', series: [{ name: 'Dynamic (compression and start of recoil)', x: thin(r.s.slice(0, k1 + 1)), y: thin(r.Fst.slice(0, k1 + 1)).map((v) => v / 1e3) }, { name: 'Air spring only', x: N.linspace(0, i.stroke, 40), y: N.linspace(0, i.stroke, 40).map((s) => gm.Fa(s) / 1e3), style: 'dash' }], annotations: [{ x: i.stroke, label: 'Stroke available' }] },
        { type: 'line', title: 'Strut stroke and tyre deflection', xlabel: 'Time [s]', ylabel: 'Displacement [m]', series: [{ name: 'Strut stroke', x: thin(r.t), y: thin(r.s) }, { name: 'Tyre deflection', x: thin(r.t), y: thin(r.d) }] },
        { type: 'line', title: 'Ground reaction factor versus sink speed', xlabel: 'Sink speed [m/s]', ylabel: 'Reaction factor [-]', series: [{ name: 'Simulated', x: sinks, y: sw.map((v) => v.Nr), style: 'line+points' }], annotations: [{ y: i.N_design, label: 'Design' }, { x: i.sink, label: 'Design sink speed' }] },
      ],
      warnings, models: ['Polytropic air spring sized from static pressure and compression ratios', 'Velocity-squared orifice damping with stronger recoil damping', 'Linear tyre spring with light damping', 'Two-degree-of-freedom drop, fixed-step RK4', 'Orifice area by golden-section search for minimum peak ground load'],
      assumptions: ['Vertical motion of one strut with a fixed share of aircraft mass; no pitch, spin-up or spring-back drag loads', 'Constant orifice (no metering pin), no seal or bearing friction', 'Lift constant during the stroke', 'Strut and airframe rigid apart from the shock absorber and tyre'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [1000, 2000, 4000, 8000], metric: 'gear_load_N' },
  verify() {
    // light unsprung mass on a very stiff tyre, no damping, lift = weight: ½·m·V² = U_air(s) + F_air(s)²/(2·kt) at maximum stroke
    const b = { m_strut: 1000, unsprung_frac: 0.002, stroke: 0.3, p_static: 5e6, r_s: 4, r_c: 3, n_poly: 1.3, A_o: 1, Cd: 0.8, k_reb: 1, k_tyre: 4e7, zeta_tyre: 0, sink: 2, lift_frac: 1, N_design: 3, strut_eff_ref: 0.8, t_end: 0.25, nSteps: 60000 };
    const gm = strutGeom(b), r = dropSim(b, gm, 1, 2, 60000), sRef = N.brent((s) => gm.Ua(s) + gm.Fa(s) ** 2 / (2 * gm.kt) - 0.5 * 1000 * 4, 0, 0.3, 1e-12);
    const bd = { ...b, k_tyre: 4e5, unsprung_frac: 0.02, zeta_tyre: 0.03 }, d = dropSim(bd, strutGeom(bd), 2e-5, 2, 8000);
    return [
      N.check('Maximum stroke from the energy integral (undamped)', r.sMax, sRef, 2e-3, 'Energy conservation with a polytropic air spring in series with the tyre'),
      N.check('Static position from the compression ratios', gm.Fa(gm.sStatic) / gm.W, 4 ** 0.3, 1e-9, 'Isothermal static pressure ratio 4, polytropic exponent 1.3'),
      N.check('Energy audit of a damped drop', d.eErr, 0, 5e-3, 'Kinetic + potential = stored + dissipated'),
    ];
  },
  calibration: { params: [{ key: 'Cd', min: 0.4, max: 1 }, { key: 'n_poly', min: 1, max: 1.4 }, { key: 'k_tyre', min: 100, max: 1e8 }], sweep: 'sink', target: 'gear_load_factor', note: 'Supply measured reaction factor against sink speed from landing-gear drop tests.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.gear_load_factor > i.N_design) out.push({ severity: 'critical', title: 'Landing load exceeds the design reaction factor', detail: `${o.gear_load_factor.toFixed(2)} against ${i.N_design.toFixed(2)} at ${i.sink.toFixed(2)} m/s sink.`, action: 'Lengthen the stroke, re-tune the orifice (or add a metering pin), or lower the sink-speed requirement; pass the peak load to Suite 2 and the exceedance to Suite 9.', basis: 'Limit landing condition' });
    if (o.stroke_used_m > 0.95 * i.stroke) out.push({ severity: 'warn', title: 'Little stroke reserve', detail: `${((100 * o.stroke_used_m) / i.stroke).toFixed(0)}% of the stroke is used.`, action: 'Reserve-energy landings (about 1.2 times the limit sink speed) need remaining stroke: raise the air-spring compression ratio or the stroke.', basis: 'Reserve energy absorption' });
    if (o.strut_efficiency < 0.7) out.push({ severity: 'advise', title: 'Strut efficiency is low', detail: `η = ${o.strut_efficiency.toFixed(2)}; good oleos reach 0.8–0.9.`, action: 'A metering pin or better orifice sizing flattens the load–stroke curve, lowering peak load and therefore gear and attachment mass.', basis: 'Load–stroke efficiency' });
    else out.push({ severity: 'info', title: 'Peak load sets structural mass', detail: `Peak ground load ${(o.gear_load_N / 1e3).toFixed(1)} kN per main gear at η = ${o.strut_efficiency.toFixed(2)}.`, action: 'Every reduction in reaction factor lowers gear, wing-attachment and fuselage loads; use this load in Suite 2 and Suite 23 mass trades.', basis: 'Landing load factor' });
    return out;
  },
};

// ---- 2. braking and stopping -------------------------------------------------------------------
// Burckhardt tyre–road friction coefficients (generic automotive data, used here only for the curve shape)
const SURF = { 'Dry': [1.2801, 23.99, 0.52, 0.75], 'Wet': [0.857, 33.822, 0.347, 0.45], 'Snow': [0.1946, 94.129, 0.0646, 0.2], 'Ice': [0.05, 306.39, 0, 0.05] }, SURFS = Object.keys(SURF);
const BMODE = ['Anti-skid (slip regulation)', 'No anti-skid (full brake pressure)', 'Ideal (peak friction throughout)'], HS = { 'Carbon': 1300, 'Steel': 520 }; // mean specific heat [J/kg/K], typical
function muCurve(i) {
  const [c1, c2, c3] = SURF[i.surface] || SURF.Dry, raw = (s) => c1 * (1 - Math.exp(-c2 * s)) - c3 * s, sp = Math.min(1, Math.log((c1 * c2) / Math.max(c3, 1e-9)) / c2), k = i.mu_peak / raw(sp);
  return { mu: (s) => k * raw(N.clamp(s, 0, 1)), sp, slope0: k * c1 * c2, lock: k * raw(1) };
}
function brakeSim(i, mode = i.mode) {
  const mc = muCurve(i), W = i.mass * G0, rho = i.rho, r = i.r_tyre, nb = Math.max(1, Math.round(i.n_brakes)), V0 = Math.max(0.1, i.V_td - i.headwind), st = 0.8 * mc.sp, band = 0.4 * mc.sp, tauB = 0.03;
  const f = (y) => {
    const [V, w, , Tb] = y, Va = V + i.headwind, q = 0.5 * rho * Va * Math.abs(Va), Nn = Math.max(0, W - q * i.CLS), Fz = (i.W_frac * Nn) / nb, slip = N.clamp((V - w * r) / Math.max(V, 0.05), 0, 1);
    let Fx, wd, Tc;
    if (mode === BMODE[2]) { Fx = Math.min(i.T_max / r, mc.mu(mc.sp) * Fz); wd = (V - 0) * 0; Tc = Fx * r; return [(-nb * Fx - q * i.CDS - i.mu_roll * Nn * (1 - i.W_frac)) / i.mass, 0, V, 0, nb * Fx * V / nb]; }
    Fx = mc.mu(slip) * Fz; Tc = mode === BMODE[0] ? i.T_max * N.clamp(1 - (slip - st) / band, 0, 1) : i.T_max;
    wd = (Fx * r - Tb) / i.I_wheel; if (w <= 0 && wd < 0) wd = 0;
    return [(-nb * Fx - q * i.CDS - i.mu_roll * Nn * (1 - i.W_frac)) / i.mass, wd, V, (Tc - Tb) / tauB, Tb * Math.max(w, 0)];
  };
  // state [V, wheel speed, distance, brake torque, energy per brake]; step bound by the wheel-slip time constant
  let y = [V0, V0 / r, 0, 0, 0], t = 0, steps = 0, muSum = 0, tOut = 0; const T = [0], Vh = [V0], Sh = [0], Xh = [0], Th = [0], Vc = Math.max(1.5, 0.03 * V0), nRes = Math.max(0.5, i.nRes), dtOut = V0 / (0.15 * G0) / 300;
  while (y[0] > Vc && steps++ < 400000) {
    const Fz = (i.W_frac * W) / nb, lam = Math.max((mc.slope0 * Fz * r * r) / (i.I_wheel * y[0]), 1 / tauB, 0.2 * G0 / y[0]), h = mode === BMODE[2] ? Math.min(0.02, 1 / (nRes * 20)) : 1 / (nRes * lam);
    const k1 = f(y), k2 = f(N.vadd(y, k1, h / 2)), k3 = f(N.vadd(y, k2, h / 2)), k4 = f(N.vadd(y, k3, h));
    y = y.map((v, j) => v + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])); if (y[1] < 0) y[1] = 0; t += h;
    const slip = N.clamp((y[0] - y[1] * r) / y[0], 0, 1); muSum += h * (mode === BMODE[2] ? mc.mu(mc.sp) : mc.mu(slip));
    if (t - tOut >= dtOut) { tOut = t; T.push(t); Vh.push(y[0]); Sh.push(slip); Xh.push(y[2]); Th.push(y[3]); }
  }
  // finish the last metres at the current deceleration
  const dec = Math.max(0.05, -f(y)[0]), tRem = y[0] / dec, x = y[2] + y[0] ** 2 / (2 * dec), Eb = y[4] + (mode === BMODE[2] ? 0 : 0.5 * (y[3] * y[1]) * tRem);
  T.push(t + tRem); Vh.push(0); Sh.push(Sh[Sh.length - 1]); Xh.push(x); Th.push(y[3]);
  return { T, Vh, Sh, Xh, Th, dist: x, time: t + tRem, Eb, muMean: muSum / Math.max(t, 1e-9), mc, V0, steps, decMean: V0 / (t + tRem) };
}
const brake = {
  id: 'brake', title: 'Braking, anti-skid and stopping distance', fidelity: 'numerical',
  summary: 'Ground roll from touchdown to stop with wheel-slip dynamics on dry, wet, snow or icy runways, with and without anti-skid, plus brake energy, heat-sink temperature and the hydroplaning speed.',
  equations: ['Braking dynamics equations', 'Wheel rotational dynamics equations', 'Tyre force equations', 'Coulomb friction law', 'Newton–Euler equations'],
  applicable: wheeled,
  inputs: [
    num('mass', 'Landing mass', 'kg', 66000, 0.1, 6e5, 'Aircraft'), num('V_td', 'Touchdown airspeed', 'm/s', 68, 2, 120, 'Aircraft'), num('headwind', 'Headwind component', 'm/s', 0, -15, 30, 'Runway', 'Filled from the site weather; negative is a tailwind'),
    num('rho', 'Air density', 'kg/m³', 1.225, 0.3, 1.5, 'Runway'), num('CDS', 'Drag area on the ground (CD·S, spoilers out)', 'm²', 12, 0, 500, 'Aircraft'), num('CLS', 'Residual lift area on the ground (CL·S)', 'm²', 12, -100, 2000, 'Aircraft'),
    num('W_frac', 'Share of weight on braked wheels', '-', 0.9, 0.1, 1, 'Aircraft'), num('mu_roll', 'Rolling friction of unbraked wheels', '-', 0.02, 0, 0.3, 'Runway'),
    sel('surface', 'Runway surface condition', SURFS, 'Dry', 'Runway', 'Sets the shape of the friction–slip curve (generic Burckhardt coefficients)'),
    num('mu_peak', 'Peak tyre friction coefficient', '-', 0.7, 0.02, 1.2, 'Runway', 'Aircraft tyres: dry ≈ 0.6–0.8, wet ≈ 0.3–0.5, compacted snow ≈ 0.2, ice ≈ 0.05 (typical)'),
    sel('mode', 'Brake control', BMODE, BMODE[0], 'Brakes'), num('n_brakes', 'Braked wheels', '', 4, 1, 40, 'Brakes', '', { step: 1, discrete: true }),
    num('T_max', 'Brake torque applied per wheel', 'N·m', 28000, 0.001, 5e6, 'Brakes', 'Default corresponds to a medium autobrake setting (about 0.3 g); maximum-effort torque is roughly twice this'), num('r_tyre', 'Tyre rolling radius', 'm', 0.56, 0.01, 2, 'Brakes'), num('I_wheel', 'Wheel, tyre and brake rotor inertia', 'kg·m²', 20, 1e-7, 5000, 'Brakes'),
    sel('heat_sink', 'Brake heat-sink material', Object.keys(HS), 'Carbon', 'Brakes', 'Carbon stores about 2.5 times more heat per kg than steel and wears more slowly'),
    num('m_hs', 'Heat-sink mass per brake', 'kg', 60, 1e-4, 2000, 'Brakes'), num('T0', 'Initial brake temperature', 'K', 320, 230, 800, 'Brakes'), num('T_limit', 'Heat-sink temperature limit', 'K', 1150, 400, 2500, 'Brakes', 'Fuse plugs release near 450–470 K at the wheel rim; heat-sink limits depend on the material (typical values)'),
    num('p_tyre', 'Tyre inflation pressure', 'Pa', 1.4e6, 3e4, 2.5e6, 'Brakes', 'For the hydroplaning speed'), num('runway_m', 'Runway length available for the ground roll', 'm', 2000, 20, 6000, 'Runway'),
    num('nRes', 'Time steps per wheel-slip time constant', '', 2, 1, 32, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c, up, d) => {
    const m = landMass(c), V = vTouch(c, up), s = c.site || {}, mb = s.runway_mu_brake ?? 0.4, surface = mb >= 0.3 ? 'Dry' : mb >= 0.16 ? 'Wet' : mb >= 0.08 ? 'Snow' : 'Ice', r = 0.48 * (c.gear.tyre_dia_m || 0.3), nb = Math.max(1, c.gear.n_main * c.gear.tyres_per_strut), S = c.wing.S_m2 || 0;
    const rho = isa(s.elev_m ?? 0, (s.T_C ?? 15) - (15 - 0.0065 * (s.elev_m || 0))).rho, Fz = (0.9 * m * G0) / nb, mw = 0.004 * c.mass.mtow_kg, big = c.mass.mtow_kg > 5700;
    return { mass: m, V_td: V, headwind: headwind(c), rho, CDS: S ? S * (c.aero.CD0 * 3 + 0.05) : c.rotor.flat_plate_m2 || 0.5, CLS: S ? 0.1 * S : 0, surface, mu_peak: N.clamp(SURF[surface][3] * (surface === 'Dry' ? mb / 0.4 : 1), 0.03, 1), n_brakes: nb, T_max: 0.35 * Fz * r, r_tyre: r, I_wheel: Math.max(1e-6, 0.6 * mw * r * r), heat_sink: big ? 'Carbon' : 'Steel', m_hs: Math.max(1e-3, (big ? 8.5e-4 : 1.6e-3) * c.mass.mtow_kg * (4 / Math.max(nb, 2))), T0: (s.T_C ?? 15) + 273.15 + 30, p_tyre: big ? 1.4e6 : c.mass.mtow_kg > 600 ? 3e5 : 1.5e5, runway_m: 0.6 * (s.runway_len_m || 2500), W_frac: c.gear.type === 'tailwheel' ? 0.85 : 0.9 };
  },
  run(i) {
    const r = brakeSim(i), alt = brakeSim(i, i.mode === BMODE[1] ? BMODE[0] : BMODE[1]), ideal = brakeSim(i, BMODE[2]), cp = HS[i.heat_sink] || 1300, dT = r.Eb / (i.m_hs * cp), Tb = i.T0 + dT, warnings = [];
    const Vh = 9 * 0.514444 * Math.sqrt(i.p_tyre / 6894.757), Vh2 = (7.7 / 9) * Vh, KE = 0.5 * i.mass * r.V0 ** 2, nb = Math.round(i.n_brakes), eff = r.muMean / i.mu_peak, lockI = i.T_max > r.mc.mu(r.mc.sp) * ((i.W_frac * i.mass * G0) / nb) * i.r_tyre;
    if (i.mode === BMODE[1] && lockI) warnings.push('Without anti-skid the commanded torque exceeds the peak tyre friction: the wheels lock, braking falls to the sliding value and the tyres will flat-spot or burst.');
    if (i.surface !== 'Dry' && r.V0 > Vh2) warnings.push(`Touchdown speed is above the hydroplaning speed of a non-rotating tyre (${Vh2.toFixed(0)} m/s): with standing water the wheels may not spin up and braking is lost until the aircraft slows.`);
    if (r.dist > i.runway_m) warnings.push('The braked ground roll exceeds the runway length available.');
    if (Tb > i.T_limit) warnings.push('Brake heat-sink temperature exceeds the stated limit: expect fuse-plug release, brake fade or fire risk; a cooling period is needed before the next take-off.');
    if (i.headwind < -5) warnings.push('Tailwind above 5 m/s is beyond normal landing limits.');
    const alts = i.mode === BMODE[1] ? 'With anti-skid' : 'Without anti-skid';
    return {
      kpis: [
        kp('stop_dist_m', 'Braked ground roll', r.dist, 'm', r.dist <= 0.85 * i.runway_m ? 'ok' : r.dist <= i.runway_m ? 'warn' : 'bad', `Available ${i.runway_m.toFixed(0)} m`), kp('stop_time_s', 'Time to stop', r.time, 's'),
        kp('decel_mean_g', 'Mean deceleration', r.decMean / G0, 'g'), kp('antiskid_eff', 'Mean friction used / peak friction', eff, '-', !lockI ? undefined : eff > 0.85 ? 'ok' : 'warn', lockI ? 'Friction-limited stop' : 'Torque-limited stop: the tyres are not at their friction limit'),
        kp('stop_dist_alt_m', `Ground roll ${alts.toLowerCase()}`, alt.dist, 'm'), kp('stop_dist_ideal_m', 'Ground roll at peak friction (ideal)', ideal.dist, 'm'),
        kp('brake_energy_J', 'Energy absorbed per brake', r.Eb, 'J'), kp('brake_energy_frac', 'Share of kinetic energy taken by the brakes', (nb * r.Eb) / KE, '-'),
        kp('brake_temp_K', 'Heat-sink temperature after the stop', Tb, 'K', Tb < 0.85 * i.T_limit ? 'ok' : Tb < i.T_limit ? 'warn' : 'bad', `Limit ${i.T_limit.toFixed(0)} K`), kp('brake_dT_K', 'Heat-sink temperature rise', dT, 'K'),
        kp('energy_loading_J_kg', 'Heat-sink energy loading', r.Eb / i.m_hs, 'J/kg'), kp('hydroplane_ms', 'Hydroplaning speed (rotating tyre)', Vh, 'm/s', undefined, `${(Vh / 0.514444).toFixed(0)} kt; ${(Vh2 / 0.514444).toFixed(0)} kt for a non-rotating tyre`),
        kp('slip_opt', 'Slip ratio at peak friction', r.mc.sp, '-'), kp('mu_locked', 'Locked-wheel friction coefficient', r.mc.lock, '-'),
      ],
      plots: [
        { type: 'line', title: 'Speed versus distance', xlabel: 'Distance from brake application [m]', ylabel: 'Ground speed [m/s]', series: [{ name: i.mode, x: r.Xh, y: r.Vh }, { name: alts, x: alt.Xh, y: alt.Vh }, { name: 'Ideal peak friction', x: ideal.Xh, y: ideal.Vh, style: 'dash' }], annotations: [{ x: i.runway_m, label: 'Runway available' }] },
        { type: 'line', title: 'Wheel slip during the stop', xlabel: 'Time [s]', ylabel: 'Slip ratio [-]', series: [{ name: i.mode, x: r.T, y: r.Sh }], annotations: [{ y: r.mc.sp, label: 'Peak friction' }] },
        { type: 'line', title: 'Tyre friction versus slip', xlabel: 'Slip ratio [-]', ylabel: 'Friction coefficient [-]', series: SURFS.map((s) => { const m = muCurve({ surface: s, mu_peak: s === i.surface ? i.mu_peak : SURF[s][3] }), x = N.linspace(0, 1, 60); return { name: s + (s === i.surface ? ' (selected)' : ''), x, y: x.map(m.mu) }; }) },
      ],
      warnings, models: ['Point-mass ground roll with aerodynamic drag and residual lift', 'Wheel spin dynamics with Burckhardt friction–slip curve scaled to the stated peak (generic curve shape)', 'Proportional slip-limiting anti-skid with first-order brake lag', 'Lumped heat-sink temperature rise', 'Horne hydroplaning speed 9·√p (kt, psi)'],
      assumptions: ['Brakes applied at touchdown speed on all braked wheels equally; no reverse thrust', 'Constant vertical load split between braked and unbraked wheels', 'Friction curve independent of speed; no hydroplaning dynamics in the roll itself', 'All brake energy stays in the heat sink during the stop (no cooling)'],
    };
  },
  convergence: { param: 'nRes', label: 'Steps per wheel-slip time constant', levels: [1, 2, 4, 8], metric: 'stop_dist_m' },
  verify() {
    const b = { mass: 10000, V_td: 50, headwind: 0, rho: 1.225, CDS: 0, CLS: 0, W_frac: 1, mu_roll: 0, surface: 'Dry', mu_peak: 0.6, mode: BMODE[2], n_brakes: 2, T_max: 1e9, r_tyre: 0.4, I_wheel: 5, heat_sink: 'Steel', m_hs: 50, T0: 300, T_limit: 1000, p_tyre: 200 * 6894.757, runway_m: 1000, nRes: 2 };
    const id = brakeSim(b), lk = brakeSim({ ...b, mode: BMODE[1] }), as = brakeSim({ ...b, mode: BMODE[0], T_max: 0.9 * 0.6 * 5000 * G0 * 0.4 }), o = N.kv(brake.run(b));
    return [
      N.check('Ideal braking distance V²/(2·μ·g)', id.dist, 2500 / (2 * 0.6 * G0), 2e-3, 'Constant-deceleration kinematics'),
      N.check('Locked-wheel distance V²/(2·μ_lock·g)', lk.dist, 2500 / (2 * lk.mc.lock * G0), 0.01, 'Sliding friction after rapid lock-up'),
      N.check('Torque-limited stop: brakes absorb the kinetic energy less tyre slip work', (2 * as.Eb) / (0.5 * 10000 * 2500), 1 - N.mean(as.Sh.slice(2, -2)), 0.02, 'Energy split between brake and tyre contact patch'),
      N.check('Horne hydroplaning speed at 200 psi', o.hydroplane_ms / 0.514444, 9 * Math.sqrt(200), 1e-9, 'Horne & Dreher, NASA TN D-2056'),
    ];
  },
  calibration: { params: [{ key: 'mu_peak', min: 0.02, max: 1.2 }, { key: 'CDS', min: 0, max: 200 }], sweep: 'V_td', target: 'stop_dist_m', note: 'Supply measured braked ground roll against brake-application speed from landing or rejected take-off tests.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.stop_dist_m > i.runway_m) out.push({ severity: 'critical', title: 'Cannot stop on the runway available', detail: `Ground roll ${o.stop_dist_m.toFixed(0)} m against ${i.runway_m.toFixed(0)} m.`, action: 'Reduce landing mass or approach speed, use reverse thrust or a longer/drier runway; update Suite 5 landing distances.', basis: 'Stopping distance' });
    if (i.mode === BMODE[1] && o.stop_dist_alt_m < 0.9 * o.stop_dist_m) out.push({ severity: 'warn', title: 'Anti-skid would shorten the stop and save tyres', detail: `${o.stop_dist_alt_m.toFixed(0)} m with anti-skid against ${o.stop_dist_m.toFixed(0)} m without.`, action: 'Fit slip-regulating anti-skid: it holds the tyre near peak friction and prevents flat spots and bursts.', basis: 'Friction–slip curve' });
    if (o.brake_temp_K > 0.85 * i.T_limit) out.push({ severity: o.brake_temp_K > i.T_limit ? 'critical' : 'warn', title: 'Brakes run hot', detail: `Heat sink reaches ${o.brake_temp_K.toFixed(0)} K (energy loading ${(o.energy_loading_J_kg / 1e6).toFixed(2)} MJ/kg).`, action: 'Increase heat-sink mass, use reverse thrust or allow a longer roll with lighter braking; plan brake cooling time before the next departure and pass the temperature to Suite 12.', basis: 'Heat-sink temperature limit' });
    out.push({ severity: 'info', title: 'Brake wear and energy', detail: `${((o.brake_energy_J * Math.round(i.n_brakes)) / 1e6).toFixed(1)} MJ is turned into heat on each landing.`, action: i.heat_sink === 'Carbon' ? 'Carbon brakes wear mainly per application: fewer, firmer applications and idle reverse extend life. Electric taxi or regenerative wheel motors can recover part of this energy on electric aircraft.' : 'Steel brakes wear with energy absorbed: use aerodynamic braking and reverse thrust where available; carbon heat sinks save mass and last longer.', basis: 'Brake energy and wear behaviour' });
    return out;
  },
};

// ---- 3. ground loads and stability -------------------------------------------------------------
const LAY = ['tricycle', 'tailwheel', 'skid'], STEER = ['Engaged / locked', 'Free castoring'];
function groundCalc(i) {
  const W = i.mass * G0, B = i.wheelbase, hT = i.track / 2, h = i.h_cg, tri = i.layout === LAY[0], tw = i.layout === LAY[1];
  // aux = nose wheel (tricycle), tail wheel (tailwheel) or the rear contact points (skid)
  const aMain = i.aux_frac * B, aAux = B - aMain, Faux = W * i.aux_frac, Fmain = W - Faux, del = Math.atan(hT / B), dTip = aAux * Math.sin(del);
  const overturn = i.layout === LAY[2] ? Math.atan(h / hT) : Math.atan(h / Math.max(dTip, 1e-9)), tip = Math.atan(aMain / h);
  const FauxBrk = tri ? (W * (aMain + i.decel_g * h)) / B : tw ? Math.max(0, (W * (aMain - i.decel_g * h)) / B) : Faux, nTy = Math.max(1, i.n_main * i.tyres), Fty = Fmain / nTy;
  // linear bicycle model: front axle at +af, rear at −ar from the CG
  const Cm = i.c_alpha * Fmain, Ca = i.steer === STEER[0] ? i.aux_ca_ratio * i.c_alpha * Faux : 0, af = tri ? aAux : aMain, ar = tri ? aMain : aAux, Cf = tri ? Ca : Cm, Cr = tri ? Cm : Ca;
  const A = (U) => [[-(Cf + Cr) / (i.mass * U), -U - (af * Cf - ar * Cr) / (i.mass * U)], [-(af * Cf - ar * Cr) / (i.Izz * U), -(af * af * Cf + ar * ar * Cr) / (i.Izz * U)]];
  const eig2 = (U) => { const a = A(U), tr = a[0][0] + a[1][1], dt = a[0][0] * a[1][1] - a[0][1] * a[1][0], disc = tr * tr - 4 * dt; return { re: disc >= 0 ? (tr + Math.sqrt(disc)) / 2 : tr / 2, im: disc >= 0 ? 0 : Math.sqrt(-disc) / 2, det: dt, tr }; };
  const ov = af * Cf - ar * Cr, over = Math.abs(ov) < 1e-9 * (af * Cf + ar * Cr) ? 0 : ov, Ucrit = over > 0 && Cf > 0 && Cr > 0 ? Math.sqrt((Cf * Cr * B * B) / (i.mass * over)) : over > 0 ? 0 : Infinity;
  return { W, aMain, aAux, Faux, Fmain, overturn, tip, FauxBrk, Fty, nTy, eig2, Ucrit, Cf, Cr, af, ar, tri, tw, hT, over };
}
const ground = {
  id: 'ground', title: 'Ground loads, tip-over and directional stability', fidelity: 'analytical',
  summary: 'Static and braking wheel loads from the centre-of-gravity position, tip-back and turnover angles, tyre and pavement loading, crosswind side load, and the tendency to ground-loop from a linear directional-stability model.',
  equations: ['Newton–Euler equations', 'Tyre force equations', 'Coulomb friction law', 'Multibody equations of motion'],
  inputs: [
    sel('layout', 'Gear layout', LAY, LAY[0], 'Geometry'), num('mass', 'Aircraft mass', 'kg', 78000, 0.1, 6e5, 'Aircraft'), num('wheelbase', 'Wheelbase (or skid contact length)', 'm', 12.6, 0.05, 40, 'Geometry'), num('track', 'Track', 'm', 7.6, 0.05, 20, 'Geometry'),
    num('aux_frac', 'Static load share on the nose or tail wheel', '-', 0.1, 0.02, 0.6, 'Geometry', 'Nose wheel 0.08–0.15 for good steering without overload; tail wheel 0.05–0.1; skids about 0.5'),
    num('h_cg', 'CG height above the ground', 'm', 3.2, 0.02, 10, 'Geometry'), num('n_main', 'Main gear struts', '', 2, 1, 6, 'Tyres', '', { step: 1, discrete: true }), num('tyres', 'Tyres per main strut', '', 2, 1, 8, 'Tyres', '', { step: 1, discrete: true }),
    num('tyre_rated_N', 'Tyre rated static load', 'N', 2.1e5, 1, 1e6, 'Tyres', 'From the tyre catalogue'), num('p_tyre', 'Tyre pressure', 'Pa', 1.4e6, 3e4, 2.5e6, 'Tyres'),
    num('decel_g', 'Braking deceleration', 'g', 0.31, 0, 0.8, 'Loads', '0.31 g (10 ft/s²) is a customary dynamic braking case'), num('tail_down_deg', 'Tail-down (rotation) angle available', 'deg', 12, 0, 30, 'Geometry'),
    num('crosswind', 'Crosswind component', 'm/s', 5, 0, 30, 'Loads'), num('V_td', 'Touchdown speed', 'm/s', 68, 2, 120, 'Loads'), num('mu_side', 'Tyre side-friction coefficient', '-', 0.5, 0.05, 1, 'Loads'),
    num('Izz', 'Yaw inertia', 'kg·m²', 4.9e6, 1e-4, 1e9, 'Directional stability'), num('c_alpha', 'Tyre cornering stiffness per unit vertical load', '1/rad', 4, 0.5, 20, 'Directional stability', 'Aircraft tyres roughly 3–6 per radian (typical, falls with load and speed)'),
    sel('steer', 'Nose or tail wheel', STEER, STEER[0], 'Directional stability', 'A free-castoring wheel gives no side force'),
    num('aux_ca_ratio', 'Nose/tail tyre cornering coefficient relative to the main tyres', '-', 1, 0, 5, 'Directional stability', 'With equal coefficients the layout is neutral-steer; a relatively stiffer front axle gives oversteer and a critical speed'),
  ],
  defaults: (c, up) => {
    const lay = LAY.includes(c.gear.type) ? c.gear.type : LAY[0], skid = lay === LAY[2], dia = c.fuselage.dia_m || 0.3, h = skid ? Math.max(0.05, 0.5 * dia) : 0.5 * (c.gear.tyre_dia_m || 0.2) + c.gear.stroke_m + 0.55 * dia, len = c.fuselage.len_m || 1, m = c.mass.mtow_kg, af = skid ? 0.5 : lay === LAY[1] ? 0.07 : m > 5700 ? 0.1 : 0.18;
    const nT = Math.max(1, c.gear.n_main * c.gear.tyres_per_strut);
    return { layout: lay, mass: m, wheelbase: c.gear.wheelbase_m, track: c.gear.track_m, aux_frac: af, h_cg: h, n_main: c.gear.n_main, tyres: c.gear.tyres_per_strut, tyre_rated_N: (1.25 * m * G0 * (1 - af)) / nT, p_tyre: m > 5700 ? 1.4e6 : m > 600 ? 3e5 : 1.5e5, crosswind: Math.abs((c.site?.wind_ms || 0) * Math.sin(N.rad((c.site?.wind_dir_deg || 0) - (c.site?.runway_heading_deg ?? c.site?.wind_dir_deg ?? 0)))) || 0.15 * vTouch(c, up), V_td: vTouch(c, up), Izz: c.mass.Izz || m * (0.25 * len) ** 2, steer: lay === LAY[1] ? STEER[1] : STEER[0], tail_down_deg: lay === LAY[0] ? 12 : 0 };
  },
  run(i) {
    const r = groundCalc(i), skid = i.layout === LAY[2], warnings = [], ovd = N.deg(r.overturn), tipd = N.deg(r.tip), crab = N.deg(Math.asin(N.clamp(i.crosswind / i.V_td, -1, 1))), side = i.mu_side * r.Fmain / Math.max(1, i.n_main);
    const tyM = i.tyre_rated_N / r.Fty - 1, area = r.Fty / i.p_tyre, rollM = r.hT / i.h_cg - i.mu_side, Us = N.linspace(Math.max(1, 0.05 * i.V_td), 1.3 * i.V_td, 50), ev = Us.map(r.eig2), e0 = r.eig2(i.V_td), t2 = e0.re > 1e-9 ? Math.log(2) / e0.re : Infinity;
    if (!skid && ovd > 63) warnings.push(`Turnover angle ${ovd.toFixed(0)}° exceeds the customary 63° limit: the aircraft may roll over in a sharp turn. Widen the track or lower the CG.`);
    if (r.tri && tipd < i.tail_down_deg) warnings.push('The CG lies behind the main-wheel contact point at the tail-down attitude: the aircraft can tip onto its tail.');
    if (r.tri && (i.aux_frac < 0.06 || i.aux_frac > 0.2)) warnings.push('Nose-wheel static load is outside about 6–20% of weight: too little gives poor steering, too much makes rotation difficult and overloads the nose gear.');
    if (r.tw && i.decel_g > r.aMain / i.h_cg) warnings.push('Braking deceleration exceeds the nose-over limit a/h: the tail will lift under hard braking.');
    if (tyM < 0) warnings.push('Static main tyre load exceeds the rated tyre load.');
    if (rollM < 0) warnings.push('Side friction can overturn the aircraft before the tyres slide (track too narrow for the CG height).');
    if (e0.re > 1e-9) warnings.push(`The ground roll is directionally unstable at touchdown speed (time to double ${t2.toFixed(2)} s): a ground loop develops unless corrected by rudder, steering or differential braking.`);
    return {
      kpis: [
        kp('main_load_N', 'Static load per main strut', r.Fmain / Math.max(1, i.n_main), 'N'), kp('aux_load_N', r.tri ? 'Static nose-gear load' : r.tw ? 'Static tail-wheel load' : 'Static rear-contact load', r.Faux, 'N'),
        kp('aux_load_braking_N', r.tri ? 'Nose-gear load under braking' : 'Tail-wheel load under braking', r.FauxBrk, 'N'), kp('tyre_load_N', 'Static load per main tyre', r.Fty, 'N', tyM >= 0.07 ? 'ok' : tyM >= 0 ? 'warn' : 'bad'),
        kp('tyre_margin', 'Main tyre load margin', tyM, '-', tyM >= 0.07 ? 'ok' : tyM >= 0 ? 'warn' : 'bad', 'About 7% growth margin is customary'), kp('contact_area_m2', 'Tyre contact area', area, 'm²'), kp('contact_pressure_Pa', 'Pavement contact pressure', i.p_tyre, 'Pa'),
        kp('turnover_deg', skid ? 'Lateral rollover angle of the CG' : 'Turnover angle', ovd, 'deg', skid || ovd <= 63 ? 'ok' : 'bad', skid ? 'atan(h / half-track)' : 'Customary limit 63° for land-based aircraft'),
        kp('tipback_deg', r.tw ? 'Nose-over angle (CG behind main wheels)' : 'Tip-back angle (main wheels behind CG)', tipd, 'deg', r.tri ? (tipd >= i.tail_down_deg ? 'ok' : 'bad') : undefined),
        kp('static_roll_margin', 'Rollover margin under side friction (half-track/h − μ)', rollM, '-', rollM > 0 ? 'ok' : 'bad'), kp('crab_deg', 'Crab angle in the crosswind', crab, 'deg'), kp('side_load_N', 'Friction-limited side load per main gear', side, 'N'),
        kp('dir_eig_re', 'Directional mode growth rate at touchdown speed', e0.re, '1/s', e0.re <= 1e-9 ? 'ok' : 'warn', e0.re > 1e-9 ? 'Unstable (ground-loop tendency)' : 'Stable'), kp('t_double_s', 'Time to double heading error', Number.isFinite(t2) ? t2 : 0, 's', undefined, Number.isFinite(t2) ? '' : 'Stable: no divergence'),
        kp('U_crit_ms', 'Speed above which the roll is unstable', Number.isFinite(r.Ucrit) ? r.Ucrit : 0, 'm/s', undefined, Number.isFinite(r.Ucrit) ? '' : 'Stable at all speeds (understeer)'),
      ],
      plots: [
        { type: 'line', title: 'Directional stability of the ground roll', xlabel: 'Ground speed [m/s]', ylabel: 'Largest real part of the eigenvalues [1/s]', series: [{ name: 'Growth rate (positive = unstable)', x: Us, y: ev.map((e) => e.re) }], annotations: [{ y: 0, label: 'Neutral' }, { x: i.V_td, label: 'Touchdown' }] },
        { type: 'bar', title: 'Wheel loads', ylabel: 'Load [kN]', categories: ['Each main strut (static)', r.tri ? 'Nose gear (static)' : 'Tail / rear (static)', r.tri ? 'Nose gear (braking)' : 'Tail / rear (braking)'], series: [{ name: 'Vertical load', y: [r.Fmain / Math.max(1, i.n_main) / 1e3, r.Faux / 1e3, r.FauxBrk / 1e3] }] },
        { type: 'line', title: 'Gear footprint and centre of gravity', xlabel: 'Lateral position [m]', ylabel: 'Longitudinal position (forward +) [m]', equalAspect: true, series: [{ name: 'Contact points', x: skid ? [-r.hT, -r.hT, r.hT, r.hT, -r.hT] : r.tri ? [-r.hT, 0, r.hT, -r.hT] : [-r.hT, 0, r.hT, -r.hT], y: skid ? [r.aAux, -r.aMain, -r.aMain, r.aAux, r.aAux] : r.tri ? [-r.aMain, r.aAux, -r.aMain, -r.aMain] : [r.aMain, -r.aAux, r.aMain, r.aMain], style: 'line+points' }, { name: 'Centre of gravity', x: [0], y: [0], style: 'points' }] },
      ],
      warnings, models: ['Rigid-body statics for wheel loads with braking load transfer', 'Turnover and tip-back geometry', 'Linear two-degree-of-freedom (bicycle) directional model with tyre cornering stiffness', 'Friction-limited side load'],
      assumptions: ['Level, rigid ground and rigid gear', 'Cornering stiffness proportional to vertical load; aerodynamic yaw stiffness and damping neglected (conservative at low speed)', 'Pavement loading reported as contact pressure and single-wheel load only', 'Skid gear treated as four contact points'],
    };
  },
  verify() {
    const b = { layout: LAY[1], mass: 1000, wheelbase: 5, track: 2, aux_frac: 0.1, h_cg: 1, n_main: 2, tyres: 1, tyre_rated_N: 6000, p_tyre: 2e5, decel_g: 0.3, tail_down_deg: 0, crosswind: 5, V_td: 25, mu_side: 0.5, Izz: 2000, c_alpha: 4, steer: STEER[1], aux_ca_ratio: 1 };
    const r = groundCalc(b), W = 1000 * G0, t = groundCalc({ ...b, layout: LAY[0], steer: STEER[0] }), o = groundCalc({ ...b, layout: LAY[0], steer: STEER[0], aux_ca_ratio: 2 }), nt = groundCalc({ ...b, layout: LAY[0], steer: STEER[0] }), e = N.eig([[r.eig2(20).tr, 0], [0, 0]]);
    return [
      N.check('Moment balance about the main wheels', r.Faux * 5, W * 0.5, 1e-12, 'Statics'),
      N.check('Tailwheel with free castor: det(A) = −a·Cα/Izz (always a saddle: ground loop)', r.eig2(20).det, (-0.5 * 4 * 0.9 * W) / 2000, 1e-10, 'Linear bicycle model'),
      N.check('Tricycle nose load under braking W(a + (ax/g)·h)/B', t.FauxBrk, (W * (0.5 + 0.3)) / 5, 1e-12, 'Statics with inertia load at the CG'),
      N.check('Oversteer critical speed: det(A) = 0', o.eig2(o.Ucrit).det / Math.abs(o.eig2(0.5 * o.Ucrit).det), 0, 1e-9, 'U² = Cf·Cr·L²/(m(a·Cf − b·Cr))'),
      N.check('Load-proportional cornering stiffness gives neutral steer', nt.over / (nt.af * nt.Cf), 0, 1e-12, 'a·Cf = b·Cr when Cα ∝ static load'),
      N.check('Trace equals the sum of eigenvalues', e[0][0] + e[1][0], r.eig2(20).tr, 1e-9, 'Matrix identity'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (i.layout !== LAY[2] && o.turnover_deg > 63) out.push({ severity: 'critical', title: 'Turnover angle too large', detail: `${o.turnover_deg.toFixed(0)}° against the customary 63°.`, action: 'Increase the track, move the main gear outboard or lower the CG.', basis: 'Ground turnover criterion' });
    if (o.dir_eig_re > 1e-9) out.push({ severity: i.layout === LAY[1] ? 'advise' : 'warn', title: 'Ground-loop tendency', detail: `Directional divergence doubles in ${o.t_double_s.toFixed(2)} s at touchdown speed.`, action: i.layout === LAY[1] ? 'Inherent in tailwheel layouts: lock the tail wheel for landing, provide effective rudder and differential braking, and keep the CG close behind the main wheels.' : 'Increase nose-wheel steering stiffness or move the main gear aft.', basis: 'Linear directional stability of the ground roll' });
    if (o.tyre_margin < 0.07) out.push({ severity: o.tyre_margin < 0 ? 'critical' : 'advise', title: 'Tyre load margin is small', detail: `Static tyre load is ${(100 * (1 - o.tyre_margin / (1 + o.tyre_margin))).toFixed(0)}% of rating.`, action: 'Select the next tyre size or ply rating; allow for mass growth.', basis: 'Tyre rated load' });
    out.push({ severity: 'info', title: 'Pavement loading', detail: `Single-wheel load ${(o.tyre_load_N / 1e3).toFixed(1)} kN at ${(i.p_tyre / 1e6).toFixed(2)} MPa contact pressure.`, action: 'Lower tyre pressure and more wheels reduce pavement damage and widen the range of usable airfields; a formal ACN/PCN (ACR/PCR) assessment needs the standard pavement method.', basis: 'Contact pressure' });
    return out;
  },
};

// ---- 4. taxi over runway roughness -------------------------------------------------------------
const PROF = ['(1 − cos) bump', 'Random roughness (ISO 8608 spectrum)'], RCLASS = { 'A (very good)': 16e-6, 'B (good)': 64e-6, 'C (average)': 256e-6, 'D (poor)': 1024e-6 }, RNAMES = Object.keys(RCLASS);
/** Seeded random profile with displacement PSD Gd(n) = Gd0·(n/n0)^−2 between n1 and n2 [cycles/m]; exact band energies. */
function roughProfile(cls, seed, n1 = 0.02, n2 = 4, M = 80) {
  const G0d = RCLASS[cls] ?? 64e-6, n0 = 0.1, u = N.rng(seed), ed = N.logspace(n1, n2, M + 1), comp = N.range(M, (k) => ({ n: Math.sqrt(ed[k] * ed[k + 1]), a: Math.sqrt(2 * G0d * n0 * n0 * (1 / ed[k] - 1 / ed[k + 1])), ph: 2 * Math.PI * u() }));
  return { y: (x) => { let s = 0; for (const c of comp) s += c.a * Math.sin(2 * Math.PI * c.n * x + c.ph); return s; }, rms: Math.sqrt(N.sum(comp.map((c) => 0.5 * c.a * c.a))), rmsExact: Math.sqrt(G0d * n0 * n0 * (1 / n1 - 1 / n2)) };
}
function taxiSim(i) {
  const gm = strutGeom(i), C = gm.Ch(i.A_o > 0 ? i.A_o : gm.AoGuess(3)), { m1, m2, kt, ct } = gm, sEq = N.brent((s) => gm.Fa(s) - m1 * G0, 0, i.stroke, 1e-12), dEq = (i.m_strut * G0) / kt;
  const rp = roughProfile(i.rclass, Math.round(i.seed)), x0 = 2, prof = i.profile === PROF[0] ? (x) => (x > x0 && x < x0 + i.bump_L ? 0.5 * i.bump_H * (1 - Math.cos((2 * Math.PI * (x - x0)) / i.bump_L)) : 0) : (x) => rp.y(x) * Math.min(1, x / 3);
  const dy = (x) => (prof(x + 1e-4) - prof(x - 1e-4)) / 2e-4, nt = Math.round(i.nSteps), T = i.dist / i.V;
  // displacements up from static equilibrium: z1 sprung, z2 unsprung
  const forces = (t, y) => { const x = i.V * t, s = sEq - (y[0] - y[2]), sd = -(y[1] - y[3]), F = gm.Fs(s, sd, C), d = dEq - (y[2] - prof(x)), dd = -(y[3] - dy(x) * i.V), Ft = d > 0 ? Math.max(0, kt * d + ct * dd) : 0; return { F, Ft }; };
  const r = N.rk4((t, y) => { const { F, Ft } = forces(t, y); return [y[1], (F - m1 * G0) / m1, y[3], (Ft - F - m2 * G0) / m2]; }, 0, [0, 0, 0, 0], T, nt);
  const acc = [], Ft = [], strk = [];
  r.t.forEach((t, k) => { const f = forces(t, r.y[k]); acc.push((f.F - m1 * G0) / (m1 * G0)); Ft.push(f.Ft); strk.push(sEq - (r.y[k][0] - r.y[k][2])); });
  const ka = (i.n_poly * gm.Fa(sEq) * gm.Aa) / (gm.V0 - gm.Aa * sEq), fB = Math.sqrt((ka * kt) / (ka + kt) / m1) / (2 * Math.PI), fW = Math.sqrt((ka + kt) / m2) / (2 * Math.PI);
  return { t: r.t, x: r.t.map((t) => i.V * t), z1: r.y.map((y) => y[0]), acc, Ft, strk, prof: r.t.map((t) => prof(i.V * t)), gm, sEq, fB, fW, rp, dt: T / nt, ka };
}
const taxi = {
  id: 'taxi', title: 'Taxi and ground roll over runway roughness', fidelity: 'numerical',
  summary: 'Response of one gear leg and its share of the aircraft to a runway bump or to random roughness: vertical acceleration felt in the airframe, dynamic wheel loads and strut travel.',
  equations: ['Newton–Euler equations', 'Oleo-pneumatic shock absorber equations', 'Hydraulic damping equations', 'Tyre force equations', 'Multibody equations of motion'],
  applicable: wheeled,
  inputs: [
    ...STRUT_INPUTS.map((f) => (f.key === 'm_strut' ? { ...f, help: 'Taxi mass share on this strut' } : f)),
    num('V', 'Ground speed', 'm/s', 15, 0.5, 100, 'Runway'), sel('profile', 'Runway profile', PROF, PROF[1], 'Runway'),
    sel('rclass', 'Roughness class', RNAMES, RNAMES[1], 'Runway', 'ISO 8608 road classes used as a generic scale; good paved runways are class A–B'), num('seed', 'Random seed', '', 7, 1, 1e6, 'Runway', '', { step: 1, discrete: true }),
    num('bump_H', 'Bump height', 'm', 0.03, 0, 0.5, 'Runway'), num('bump_L', 'Bump length', 'm', 6, 0.1, 200, 'Runway', 'Wavelengths near speed/frequency of the bounce mode are the most severe'),
    num('dist', 'Distance simulated', 'm', 150, 5, 3000, 'Numerics'), num('nSteps', 'Time steps', '', 8000, 500, 400000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c, up) => ({ ...strutDefaults(c), m_strut: (0.92 * c.mass.mtow_kg) / Math.max(1, c.gear.n_main), V: N.clamp(0.25 * vTouch(c, up), 3, 20), dist: N.clamp(10 * 0.25 * vTouch(c, up), 40, 200), bump_H: N.clamp(0.25 * c.gear.stroke_m, 0.005, 0.05) }),
  run(i) {
    const r = taxiSim(i), k0 = Math.round(0.1 * r.t.length), a = r.acc.slice(k0), rms = Math.sqrt(N.mean(a.map((v) => v * v))), pk = N.amax(r.acc.map(Math.abs)), Fst = r.gm.W, dlf = N.amax(r.Ft) / Fst, warnings = [];
    const sp = N.spectrum(a, r.dt), kp1 = N.argmax(sp.amp.slice(1)) + 1, sMin = N.amin(r.strk), sMax = N.amax(r.strk), off = N.amin(r.Ft.slice(k0)) <= 0;
    if (sMax > i.stroke) warnings.push('The strut bottoms out on this profile.');
    if (sMin < 0) warnings.push('The strut tops out (fully extends) on this profile.');
    if (off) warnings.push('The tyre leaves the ground: loads on recontact are high and steering and braking are interrupted.');
    if (r.dt * 2 * Math.PI * r.fW > 0.5) warnings.push('The time step is coarse for the wheel-hop mode; increase the number of steps.');
    const nf = Math.min(sp.f.length, 200);
    return {
      kpis: [
        kp('acc_rms_g', 'RMS vertical acceleration of the airframe', rms, 'g', rms < 0.15 ? 'ok' : 'warn'), kp('acc_peak_g', 'Peak vertical acceleration', pk, 'g', pk < 0.6 ? 'ok' : 'warn'),
        kp('dyn_load_factor', 'Peak ground load / static load', dlf, '-', dlf < 1.7 ? 'ok' : 'warn'), kp('ground_load_max_N', 'Peak ground load', N.amax(r.Ft), 'N'), kp('ground_load_min_N', 'Lowest ground load', N.amin(r.Ft.slice(k0)), 'N', off ? 'warn' : 'ok'),
        kp('strut_travel_m', 'Strut travel range', sMax - sMin, 'm'), kp('stroke_static_m', 'Static strut position', r.sEq, 'm'), kp('f_bounce_Hz', 'Body bounce frequency (linearised)', r.fB, 'Hz'), kp('f_wheel_hop_Hz', 'Wheel-hop frequency (linearised)', r.fW, 'Hz'),
        kp('f_response_Hz', 'Dominant response frequency', sp.f[kp1], 'Hz'), kp('profile_rms_m', 'Profile RMS height', i.profile === PROF[0] ? i.bump_H / Math.sqrt(8) * Math.sqrt(i.bump_L / i.dist) * Math.sqrt(3) : r.rp.rms, 'm'), kp('critical_wavelength_m', 'Wavelength exciting the bounce mode at this speed', i.V / r.fB, 'm'),
      ],
      plots: [
        { type: 'line', title: 'Airframe vertical acceleration', xlabel: 'Distance [m]', ylabel: 'Acceleration [g]', series: [{ name: 'Sprung mass', x: thin(r.x, 400), y: thin(r.acc, 400) }] },
        { type: 'line', title: 'Runway profile and airframe displacement', xlabel: 'Distance [m]', ylabel: 'Height [mm]', series: [{ name: 'Runway profile', x: thin(r.x, 400), y: thin(r.prof, 400).map((v) => v * 1e3) }, { name: 'Airframe', x: thin(r.x, 400), y: thin(r.z1, 400).map((v) => v * 1e3) }] },
        { type: 'line', title: 'Ground load', xlabel: 'Distance [m]', ylabel: 'Load [kN]', series: [{ name: 'Tyre load', x: thin(r.x, 400), y: thin(r.Ft, 400).map((v) => v / 1e3) }], annotations: [{ y: Fst / 1e3, label: 'Static' }] },
        { type: 'line', title: 'Acceleration spectrum', xlabel: 'Frequency [Hz]', ylabel: 'Amplitude [g]', series: [{ name: 'Sprung-mass acceleration', x: sp.f.slice(0, nf), y: sp.amp.slice(0, nf) }], annotations: [{ x: r.fB, label: 'Bounce' }] },
      ],
      warnings, models: ['Two-degree-of-freedom leg model with the nonlinear oleo and linear tyre', 'Point-contact tyre following the profile', 'Sum-of-sines random profile from a −2 slope displacement spectrum (ISO 8608 form), seeded', 'Fixed-step RK4'],
      assumptions: ['One leg with a fixed share of aircraft mass: no pitch or roll coupling, no fuselage flexibility (cockpit acceleration is usually higher)', 'No strut seal friction: real oleos can stay locked over small roughness, leaving only the tyre as a spring', 'Constant speed, no lift', 'Road-class roughness scale used as a generic stand-in for measured runway profiles'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [2000, 4000, 8000, 16000], metric: 'acc_rms_g' },
  verify() {
    const b = { m_strut: 1000, unsprung_frac: 0.03, stroke: 0.2, p_static: 4e6, r_s: 4, r_c: 3, n_poly: 1.3, A_o: 0, Cd: 0.8, k_reb: 4, k_tyre: 4e5, zeta_tyre: 0.03, V: 2, profile: PROF[0], rclass: RNAMES[1], seed: 3, bump_H: 0.05, bump_L: 60, dist: 70, nSteps: 20000 };
    const q = taxiSim(b), f = taxiSim({ ...b, bump_H: 0 }), rp = roughProfile(RNAMES[1], 5), xs = N.linspace(0, 5000, 400001), yr = xs.map(rp.y);
    return [
      N.check('Long bump taken quasi-statically: airframe rises by the bump height', N.amax(q.z1), 0.05, 0.02, 'Static limit of the transmissibility'),
      N.check('Flat runway: aircraft stays in equilibrium', N.amax(f.acc.map(Math.abs)), 0, 1e-9, 'Static equilibrium'),
      N.check('Random profile RMS equals the integral of its spectrum', rp.rms, rp.rmsExact, 1e-9, 'σ² = ∫Gd(n)dn = Gd0·n0²(1/n1 − 1/n2)'),
      N.check('Sampled profile RMS over 5 km', Math.sqrt(N.mean(yr.map((v) => v * v))), rp.rmsExact, 0.08, 'Finite-length sample of the same spectrum'),
    ];
  },
  calibration: { params: [{ key: 'k_tyre', min: 100, max: 1e8 }, { key: 'zeta_tyre', min: 0, max: 0.3 }, { key: 'Cd', min: 0.4, max: 1 }], sweep: 'V', target: 'acc_rms_g', note: 'Supply measured RMS vertical acceleration against taxi speed from instrumented taxi tests on a surveyed runway.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.acc_rms_g > 0.15) out.push({ severity: 'advise', title: 'Rough ride during taxi', detail: `${o.acc_rms_g.toFixed(2)} g RMS, peak ${o.acc_peak_g.toFixed(2)} g.`, action: `Avoid dwelling near the speed where the dominant runway wavelength equals ${o.critical_wavelength_m.toFixed(1)} m; a softer strut air curve or two-stage oleo reduces taxi loads.`, basis: 'Bounce-mode resonance' });
    if (o.dyn_load_factor > 1.5) out.push({ severity: 'advise', title: 'Taxi loads contribute to fatigue', detail: `Ground load cycles up to ${o.dyn_load_factor.toFixed(2)} times static.`, action: 'Include the taxi load spectrum in Suite 9; ground loads are a large part of landing-gear and wing-root fatigue damage on short-haul aircraft.', basis: 'Ground–air–ground fatigue spectrum' });
    if (!out.length) out.push({ severity: 'info', title: 'Benign taxi response', detail: `${o.acc_rms_g.toFixed(3)} g RMS on a class ${i.rclass[0]} surface.`, action: 'Re-check on the roughest runway in the intended network and at high-speed turn-off.', basis: 'Roughness response' });
    return out;
  },
};

// ---- 5. shimmy --------------------------------------------------------------------------------
/** Third-order castor shimmy model (straight-tangent tyre with relaxation length): states [ψ, ψ', α]. */
function shimmySys(i, V = i.V, e = i.trail, kd = i.k_damp) {
  const CF = i.c_Fa * i.Fz, Ma = e * CF + i.t_p * CF, D = kd + i.kappa / Math.max(V, 0.1), I = i.I_z, sg = i.sigma;
  const A = [[0, 1, 0], [-i.c_tors / I, -D / I, -Ma / I], [V / sg, (e - i.a_c) / sg, -V / sg]];
  // characteristic polynomial a3 s³ + a2 s² + a1 s + a0 and the Routh–Hurwitz margin
  const a3 = (I * sg) / V, a2 = I + (D * sg) / V, a1 = D + (i.c_tors * sg) / V + (Ma * (e - i.a_c)) / V, a0 = i.c_tors + Ma;
  return { A, a3, a2, a1, a0, hurwitz: a2 * a1 - a3 * a0, Ma, CF, stableRH: a1 > 0 && a0 > 0 && a2 * a1 - a3 * a0 > 0 };
}
function shimmyEig(i, V, e, kd) {
  const s = shimmySys(i, V, e, kd), ev = N.eig(s.A); let re = -Infinity, im = 0, osc = null;
  for (const l of ev) { if (l[0] > re) re = l[0]; if (l[1] > 1e-9 && (!osc || l[0] > osc[0])) osc = l; }
  if (osc) im = osc[1];
  return { re, f: im / (2 * Math.PI), zeta: osc ? -osc[0] / Math.hypot(osc[0], osc[1]) : 1, s };
}
const shimmy = {
  id: 'shimmy', title: 'Nose- or tail-wheel shimmy stability', fidelity: 'reduced-order',
  summary: 'Stability of the castoring wheel against self-excited yaw oscillation from a linear model with tyre relaxation: growth rate and frequency against speed, the influence of mechanical trail, and the shimmy-damper rate needed.',
  equations: ['Newton–Euler equations', 'Tyre force equations', 'Wheel rotational dynamics equations', 'Multibody equations of motion'],
  applicable: wheeled,
  inputs: [
    num('V', 'Ground speed', 'm/s', 60, 0.5, 120, 'Operation'), num('Fz', 'Vertical load on the wheel', 'N', 70000, 1, 1e6, 'Operation'),
    num('I_z', 'Yaw inertia of the swivelling parts about the steering axis', 'kg·m²', 6, 1e-7, 1e4, 'Gear'), num('trail', 'Mechanical trail', 'm', 0.06, -0.5, 1, 'Gear', 'Wheel axle behind the steering axis'),
    num('c_tors', 'Torsional stiffness about the steering axis', 'N·m/rad', 1.5e5, 0, 1e9, 'Gear', 'Torque links and steering actuator'), num('k_damp', 'Shimmy damper rate', 'N·m·s/rad', 300, 0, 1e7, 'Gear'),
    num('a_c', 'Tyre contact half-length', 'm', 0.1, 0.001, 1, 'Tyre'), num('sigma', 'Tyre relaxation length', 'm', 0.3, 0.005, 3, 'Tyre', 'Typically about 3 contact half-lengths'),
    num('c_Fa', 'Cornering stiffness per unit load', '1/rad', 6, 0.5, 30, 'Tyre', 'Nose tyres at low load are typically 5–10 per radian (generic)'), num('t_p', 'Pneumatic trail', 'm', 0.05, 0, 0.5, 'Tyre', 'About half the contact half-length'),
    num('kappa', 'Tread-width damping constant', 'N·m²/rad', 100, 0, 1e6, 'Tyre', 'Yaw damping from the finite tread width, moment = κ/V·ψ′'),
    num('V_max', 'Top of the speed range', 'm/s', 90, 2, 150, 'Numerics'), num('nGrid', 'Grid points per axis of the stability map', '', 40, 10, 120, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c, up) => {
    const m = c.mass.mtow_kg, Rn = 0.3 * (c.gear.tyre_dia_m || 0.2), a = 0.33 * Rn, Fz = (c.gear.type === 'tailwheel' ? 0.07 : 0.1) * m * G0, Iz = Math.max(1e-7, 0.002 * m * (0.6 * Rn) ** 2), V = vTouch(c, up);
    const b = { V, Fz, I_z: Iz, trail: 0.3 * a, c_tors: Iz * (2 * Math.PI * 25) ** 2, a_c: a, sigma: 3 * a, c_Fa: 6, t_p: 0.5 * a, kappa: 0.15 * a * a * 6 * Fz, k_damp: 0 };
    // damper sized at 1.5 times the largest rate needed anywhere in the speed range (see run)
    let need = 0; for (const v of N.linspace(2, 1.3 * V, 30)) { const f = (k) => shimmySys(b, v, b.trail, k).hurwitz; if (f(0) < 0) need = Math.max(need, N.brent(f, 0, 1e9 * Iz + 1e6, 1e-6)); }
    return { ...b, k_damp: 1.5 * need, V_max: 1.3 * V };
  },
  run(i) {
    const e0 = shimmyEig(i, i.V, i.trail, i.k_damp), Vs = N.linspace(Math.max(1, 0.02 * i.V_max), i.V_max, 60), sw = Vs.map((v) => shimmyEig(i, v, i.trail, i.k_damp)), sw0 = Vs.map((v) => shimmyEig(i, v, i.trail, 0)), warnings = [];
    const need = (v) => { const f = (k) => shimmySys(i, v, i.trail, k).hurwitz; return f(0) >= 0 ? 0 : N.brent(f, 0, 1e9 * i.I_z + 1e6, 1e-6); }, kNeed = N.amax(Vs.map(need)), worst = N.argmax(sw.map((s) => s.re)), unst = sw.filter((s) => s.re > 1e-9).length;
    const ng = Math.round(i.nGrid), es = N.linspace(0, 1.5 * (i.a_c + i.sigma), ng), vg = N.linspace(Math.max(1, 0.03 * i.V_max), i.V_max, ng), Z = es.map((e) => vg.map((v) => N.clamp(shimmyEig(i, v, e, i.k_damp).re, -30, 30)));
    if (e0.re > 1e-9) warnings.push(`Shimmy is unstable at ${i.V.toFixed(0)} m/s: oscillation at ${e0.f.toFixed(1)} Hz grows with a time constant of ${(1 / e0.re).toFixed(2)} s.`);
    else if (unst) warnings.push(`Stable at the stated speed but unstable in part of the speed range (worst near ${Vs[worst].toFixed(0)} m/s).`);
    if (i.trail > 0 && i.trail < i.a_c + i.sigma && i.k_damp < kNeed) warnings.push('Mechanical trail lies in the unstable band 0 < e < a + σ of the undamped model and the damper is below the rate required.');
    return {
      kpis: [
        kp('shimmy_growth', 'Growth rate of the least stable mode', e0.re, '1/s', e0.re < -1e-9 ? 'ok' : 'bad', e0.re > 1e-9 ? 'Unstable' : 'Stable'), kp('shimmy_zeta', 'Damping ratio of the shimmy mode', e0.zeta, '-', e0.zeta > 0.05 ? 'ok' : e0.zeta > 0 ? 'warn' : 'bad'),
        kp('shimmy_freq_Hz', 'Shimmy frequency', e0.f, 'Hz'), kp('shimmy_wavelength_m', 'Shimmy wavelength on the runway', e0.f > 0 ? i.V / e0.f : 0, 'm'),
        kp('damper_required', 'Damper rate for stability over the speed range', kNeed, 'N·m·s/rad', i.k_damp >= 1.2 * kNeed ? 'ok' : i.k_damp >= kNeed ? 'warn' : 'bad', `Fitted ${i.k_damp.toFixed(1)}`),
        kp('damper_margin', 'Damper rate fitted / required', kNeed > 0 ? i.k_damp / kNeed : 99, '-'), kp('trail_critical_m', 'Trail above which the undamped castor is stable (a + σ)', i.a_c + i.sigma, 'm'),
        kp('worst_speed_ms', 'Least stable speed in the range', Vs[worst], 'm/s'), kp('worst_growth', 'Growth rate at the least stable speed', sw[worst].re, '1/s', sw[worst].re < -1e-9 ? 'ok' : 'bad'), kp('unstable_fraction', 'Share of the speed range that is unstable', unst / Vs.length, '-'),
      ],
      plots: [
        { type: 'line', title: 'Shimmy mode growth rate versus speed', xlabel: 'Ground speed [m/s]', ylabel: 'Real part of the eigenvalue [1/s]', series: [{ name: 'With the damper fitted', x: Vs, y: sw.map((s) => s.re) }, { name: 'No damper', x: Vs, y: sw0.map((s) => s.re), style: 'dash' }], annotations: [{ y: 0, label: 'Stability boundary' }, { x: i.V, label: 'Speed' }] },
        { type: 'line', title: 'Shimmy frequency versus speed', xlabel: 'Ground speed [m/s]', ylabel: 'Frequency [Hz]', series: [{ name: 'Oscillatory mode', x: Vs, y: sw.map((s) => s.f) }] },
        { type: 'heat', title: 'Stability map: growth rate against speed and trail', xlabel: 'Ground speed [m/s]', ylabel: 'Mechanical trail [m]', zlabel: 'Growth rate [1/s] (positive = unstable)', x: vg, y: es, z: Z, contours: 12, diverging: true, overlay: [{ name: 'Fitted trail', x: [vg[0], vg[ng - 1]], y: [i.trail, i.trail] }] },
      ],
      warnings, models: ['Castor yaw dynamics about the steering axis with torsional stiffness and viscous damper', 'Straight-tangent stretched-string tyre: first-order relaxation of the slip angle with contact-length lead', 'Linear cornering force and aligning moment; tread-width damping κ/V', 'Eigenvalues of the 3×3 state matrix; damper requirement from the Routh–Hurwitz condition'],
      assumptions: ['Small angles, constant speed and vertical load', 'Rigid gear in lateral bending (no lateral–torsional coupling) and no free play', 'Tyre data are generic; measured cornering stiffness and relaxation length should replace them', 'Linear damper (real shimmy dampers are velocity-squared with friction)'],
    };
  },
  verify() {
    const b = { V: 30, Fz: 10000, I_z: 1, trail: 0.4, c_tors: 0, k_damp: 0, a_c: 0.1, sigma: 0.3, c_Fa: 6, t_p: 0, kappa: 0, V_max: 60, nGrid: 10 };
    const e = shimmyEig(b, 30, 0.4, 0), wRef = Math.sqrt((0.4 * 6 * 10000) / 1), st = shimmyEig(b, 30, 0.5, 0), un = shimmyEig(b, 30, 0.3, 0), d = { ...b, c_tors: 5e3, k_damp: 20, kappa: 50, t_p: 0.05, trail: 0.15 }, sd = shimmySys(d), ed = N.eig(sd.A);
    // the eigenvalues must be roots of the characteristic cubic
    const l = ed[ed.length - 1], p = N.C.add(N.C.add(N.C.scale(N.C.mul(l, N.C.mul(l, l)), sd.a3), N.C.scale(N.C.mul(l, l), sd.a2)), N.C.add(N.C.scale(l, sd.a1), [sd.a0, 0]));
    return [
      N.check('Undamped stability boundary at trail = a + σ', e.re / (2 * Math.PI * e.f), 0, 1e-8, 'Routh–Hurwitz: a2·a1 − a3·a0 ∝ (e − a − σ)'),
      N.check('Frequency on the boundary ω² = e·C_F/I', 2 * Math.PI * e.f, wRef, 1e-8, 'Roots of the characteristic cubic on the boundary'),
      N.check('Stable for trail above a + σ', st.re < 0 ? 1 : 0, 1, 1e-12, 'von Schlippe–Dietrich / Moreland-type result'),
      N.check('Unstable for trail between 0 and a + σ', un.re > 0 ? 1 : 0, 1, 1e-12, 'Same criterion'),
      N.check('Eigenvalue satisfies the characteristic cubic', N.C.abs(p) / (sd.a0 + Math.abs(sd.a1) * N.C.abs(l)), 0, 1e-9, 'det(sI − A) = 0'),
    ];
  },
  calibration: { params: [{ key: 'sigma', min: 0.01, max: 2 }, { key: 'c_Fa', min: 1, max: 20 }, { key: 'kappa', min: 0, max: 1e5 }], sweep: 'V', target: 'shimmy_freq_Hz', note: 'Supply measured shimmy frequency against speed from dynamometer or taxi tests to fit the tyre relaxation length and cornering stiffness.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.worst_growth > 1e-9) out.push({ severity: o.shimmy_growth > 1e-9 ? 'critical' : 'warn', title: 'Shimmy instability', detail: `Unstable over ${(100 * o.unstable_fraction).toFixed(0)}% of the speed range; damper needed ≥ ${o.damper_required.toFixed(1)} N·m·s/rad, fitted ${i.k_damp.toFixed(1)}.`, action: 'Fit or enlarge the shimmy damper, use twin co-rotating wheels, increase torsional stiffness and remove free play, or move the trail out of the 0 … a + σ band.', basis: 'Routh–Hurwitz stability of the castor–tyre system' });
    else if (o.damper_margin < 1.5) out.push({ severity: 'advise', title: 'Shimmy stability relies on a thin damper margin', detail: `Damper is ${o.damper_margin.toFixed(2)} times the minimum.`, action: 'Damper wear, air in the damper and torque-link free play erode this margin in service: specify inspection intervals and a wear limit (Suite 22).', basis: 'Stability margin' });
    else out.push({ severity: 'info', title: 'Castor is stable across the speed range', detail: `Shimmy-mode damping ratio ${o.shimmy_zeta.toFixed(2)} at ${i.V.toFixed(0)} m/s.`, action: 'Confirm with measured tyre data and include gear lateral bending in a multibody model before freezing the design.', basis: 'Eigenvalue analysis' });
    return out;
  },
};

// ---- 6. skid-gear landing ----------------------------------------------------------------------
function skidGeom(i) {
  const m = METALS[i.mat] || METALS['Al 7075-T6'], D = i.D_tube, d = D * (1 - 2 * i.t_ratio), I = (Math.PI * (D ** 4 - d ** 4)) / 64, n = Math.max(1, Math.round(i.n_legs));
  const k = (n * 3 * m.E * I) / i.ell ** 3, Mp = (m.Sy * (D ** 3 - d ** 3)) / 6, My = (m.Sy * I) / (D / 2), Pp = (n * Mp) / i.ell, Py = (n * My) / i.ell;
  return { m, k, Pp, Py, dy: Pp / k, I, Mp, mass: n * m.rho * (Math.PI / 4) * (D * D - d * d) * i.ell };
}
/** Drop on elastic–perfectly-plastic cross tubes: state [deflection, velocity], plastic set tracked explicitly. */
function skidSim(i, gm, V, nSteps = i.nSteps) {
  const W = i.mass * G0, L = i.lift_frac * W, nt = Math.round(nSteps), tEnd = i.t_end, h = tEnd / nt; let z = 0, v = V, dp = 0, zMax = 0, Fmax = 0, Eabs = 0, Epk = 0;
  const T = [0], Z = [0], F = [0], force = (zz) => N.clamp(gm.k * (zz - dp), 0, gm.Pp);
  for (let s = 1; s <= nt; s++) {
    // velocity-Verlet with a return-mapping update of the plastic deflection
    const f0 = force(z), a0 = (W - L - f0) / i.mass, zn = z + v * h + 0.5 * a0 * h * h; if (gm.k * (zn - dp) > gm.Pp) dp = zn - gm.Pp / gm.k;
    const f1 = force(zn), a1 = (W - L - f1) / i.mass; Eabs += 0.5 * (f0 + f1) * (zn - z); v += 0.5 * (a0 + a1) * h; z = zn;
    if (z > zMax) { zMax = z; Epk = Eabs; } if (f1 > Fmax) Fmax = f1; T.push(s * h); Z.push(z); F.push(f1);
    if (v < 0 && f1 <= 0 && s > 10) break;
  }
  return { T, Z, F, zMax, Fmax, set: dp, Eabs: Epk, Eplastic: Eabs, plastic: dp > 0 };
}
const skidPeak = (i, gm, V) => { const W = i.mass * G0, net = W * (1 - i.lift_frac), KE = 0.5 * i.mass * V * V, de = (net + Math.sqrt(net * net + 2 * gm.k * KE)) / gm.k; return de <= gm.dy ? de : gm.Pp > net ? (KE + gm.Pp ** 2 / (2 * gm.k)) / (gm.Pp - net) : Infinity; };
const skid = {
  id: 'skid', title: 'Skid-gear landing: energy absorption by cross-tube bending', fidelity: 'reduced-order',
  summary: 'Vertical landing on skid gear whose cross tubes bend elastically and then plastically: peak load factor, deflection, permanent set and the reserve-energy (hard landing) case.',
  equations: ['Landing impact energy equations', 'Newton–Euler equations', 'Structural impact equations'],
  applicable: (c) => (c.gear.type === 'skid' ? true : 'This aircraft has wheeled gear: use the oleo-pneumatic drop analysis.'),
  inputs: [
    num('mass', 'Landing mass', 'kg', 6.3, 0.01, 2e4, 'Aircraft'), num('sink', 'Sink speed at touchdown', 'm/s', 1.5, 0.05, 10, 'Touchdown', 'About 2 m/s (6.5 ft/s) limit drop for small rotorcraft'),
    num('lift_frac', 'Rotor lift / weight during the impact', '-', 0.667, 0, 1, 'Touchdown', 'Two-thirds is customary for rotorcraft limit drop tests'),
    num('n_legs', 'Cross-tube ends carrying load', '', 4, 2, 8, 'Gear', 'Two cross tubes, each bending at both ends', { step: 1, discrete: true }), num('ell', 'Bending arm from fuselage attachment to skid', 'm', 0.16, 0.01, 3, 'Gear'),
    num('D_tube', 'Cross-tube outer diameter', 'm', 0.012, 0.001, 0.3, 'Gear'), num('t_ratio', 'Wall thickness / diameter', '-', 0.1, 0.02, 0.5, 'Gear'), sel('mat', 'Tube material', Object.keys(METALS), 'Al 7075-T6', 'Gear'),
    num('clearance', 'Deflection available before the fuselage touches', 'm', 0.08, 0.005, 2, 'Limits'), num('reserve', 'Reserve-energy factor on drop energy', '-', 1.5, 1, 3, 'Limits', '1.5 × limit energy is the customary reserve-energy drop'),
    num('t_end', 'Simulated time', 's', 0.5, 0.01, 5, 'Numerics'), num('nSteps', 'Time steps', '', 4000, 200, 400000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => {
    const W = c.mass.mtow_kg * G0, ell = Math.max(0.03, 0.8 * (c.gear.track_m || 0.3)), m = METALS['Al 7075-T6'], tr = 0.1, D = Math.cbrt((6 * ((3 * W) / 4) * ell) / (m.Sy * (1 - (1 - 2 * tr) ** 3)));
    return { mass: c.mass.mtow_kg, sink: c.gear.sink_ms, ell, D_tube: D, t_ratio: tr, clearance: Math.max(0.03, 0.4 * ell), t_end: N.clamp((40 * (c.gear.stroke_m || 0.03)) / Math.max(0.2, c.gear.sink_ms), 0.1, 3) };
  },
  run(i) {
    const gm = skidGeom(i), W = i.mass * G0, r = skidSim(i, gm, i.sink), Vr = i.sink * Math.sqrt(i.reserve), rr = skidSim(i, gm, Vr), nz = r.Fmax / W, warnings = [], an = skidPeak(i, gm, i.sink);
    if (r.plastic) warnings.push(`The cross tubes yield in the limit landing and take a permanent set of ${(r.set * 1e3).toFixed(1)} mm: inspection and replacement are needed.`);
    if (rr.zMax > i.clearance) warnings.push('In the reserve-energy landing the gear deflects beyond the available clearance: the fuselage or payload strikes the ground.');
    if (gm.Pp <= W * (1 - i.lift_frac)) warnings.push('The plastic collapse load is below the net weight: the gear cannot arrest the descent.');
    const ds = N.linspace(0, Math.max(1.5 * rr.zMax, 2 * gm.dy), 50), sinks = N.linspace(0.3 * i.sink, 2 * i.sink, 25);
    return {
      kpis: [
        kp('gear_load_N', 'Peak ground load per skid', r.Fmax / 2, 'N'), kp('gear_load_factor', 'Ground reaction factor (peak load / weight)', nz, '-'), kp('nz_cg', 'Peak load factor at the mass', nz + i.lift_frac, 'g'),
        kp('stroke_used_m', 'Peak gear deflection', r.zMax, 'm', r.zMax <= i.clearance ? 'ok' : 'bad', `Clearance ${i.clearance.toFixed(3)} m`), kp('permanent_set_m', 'Permanent set after the landing', r.set, 'm', r.plastic ? 'warn' : 'ok'),
        kp('absorbed_J', 'Energy absorbed by the gear at peak deflection', r.Eabs, 'J'), kp('plastic_work_J', 'Energy dissipated plastically', r.Eplastic, 'J'), kp('yield_load_factor', 'Reaction factor at first yield', gm.Py / W, '-'), kp('collapse_load_factor', 'Reaction factor at full plastic hinge', gm.Pp / W, '-'),
        kp('stiffness_Npm', 'Vertical gear stiffness', gm.k, 'N/m'), kp('deflection_reserve_m', 'Deflection in the reserve-energy landing', rr.zMax, 'm', rr.zMax <= i.clearance ? 'ok' : 'bad'),
        kp('sink_first_yield_ms', 'Sink speed at which the tubes start to yield plastically', Math.sqrt(Math.max(0, (gm.k * gm.dy ** 2 - 2 * W * (1 - i.lift_frac) * gm.dy) / i.mass)), 'm/s'), kp('tube_mass_kg', 'Mass of the bending lengths', gm.mass, 'kg'), kp('deflection_energy_m', 'Peak deflection from the energy balance', an, 'm'),
      ],
      plots: [
        { type: 'line', title: 'Gear load and deflection', xlabel: 'Time [s]', ylabel: 'Ground load / weight [-]', series: [{ name: 'Limit landing', x: thin(r.T), y: thin(r.F).map((v) => v / W) }, { name: 'Reserve-energy landing', x: thin(rr.T), y: thin(rr.F).map((v) => v / W) }] },
        { type: 'line', title: 'Load–deflection characteristic', xlabel: 'Deflection [mm]', ylabel: 'Total ground load [N]', series: [{ name: 'Elastic–perfectly-plastic cross tubes', x: ds.map((v) => v * 1e3), y: ds.map((v) => Math.min(gm.k * v, gm.Pp)) }], annotations: [{ x: r.zMax * 1e3, label: 'Limit landing' }, { x: i.clearance * 1e3, label: 'Clearance' }] },
        { type: 'line', title: 'Peak deflection versus sink speed', xlabel: 'Sink speed [m/s]', ylabel: 'Deflection [mm]', series: [{ name: 'Energy balance', x: sinks, y: sinks.map((v) => Math.min(skidPeak(i, gm, v), 10 * i.clearance) * 1e3) }], annotations: [{ y: i.clearance * 1e3, label: 'Clearance' }, { x: i.sink, label: 'Limit' }, { x: Vr, label: 'Reserve' }] },
      ],
      warnings, models: ['Cantilever cross-tube bending, elastic to the full plastic moment then perfectly plastic', 'Single-degree-of-freedom drop with constant rotor lift, velocity-Verlet integration with plastic return mapping', 'Closed-form energy balance for the peak deflection'],
      assumptions: ['Vertical, level landing shared equally by all cross-tube ends', 'No skid spreading friction, strain hardening or rate effects', 'Rigid ground and fuselage', 'Handbook yield strength, not a design allowable'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [500, 1000, 2000, 4000], metric: 'stroke_used_m' },
  verify() {
    const b = { mass: 1000, sink: 2, lift_frac: 0.667, n_legs: 4, ell: 0.5, D_tube: 0.05, t_ratio: 0.1, mat: 'Al 7075-T6', clearance: 0.3, reserve: 1.5, t_end: 0.6, nSteps: 20000 }, gm = skidGeom(b), r = skidSim(b, gm, 2), e = skidSim(b, gm, 0.3);
    const net = 1000 * G0 * (1 - 0.667), de = (net + Math.sqrt(net * net + 2 * gm.k * 0.5 * 1000 * 0.09)) / gm.k;
    return [
      N.check('Plastic landing: peak deflection from the energy balance', r.zMax, (0.5 * 1000 * 4 + gm.Pp ** 2 / (2 * gm.k)) / (gm.Pp - net), 1e-3, '½mV² + (W − L)δ = P_p·δ − P_p²/2k'),
      N.check('Elastic landing: peak deflection', e.zMax, de, 1e-3, '½mV² + (W − L)δ = ½kδ²'),
      N.check('Plastic moment of a thin tube σy(D³ − d³)/6', gm.Mp, (METALS['Al 7075-T6'].Sy * (0.05 ** 3 - 0.04 ** 3)) / 6, 1e-12, 'Plastic section modulus'),
      N.check('Energy absorbed equals drop energy at maximum deflection', r.Eabs, 0.5 * 1000 * 4 + net * r.zMax, 1e-3, 'Work–energy theorem'),
    ];
  },
  calibration: { params: [{ key: 'ell', min: 0.01, max: 3 }, { key: 't_ratio', min: 0.02, max: 0.5 }], sweep: 'sink', target: 'stroke_used_m', note: 'Supply measured peak gear deflection against drop speed from skid-gear drop tests.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.deflection_reserve_m > i.clearance) out.push({ severity: 'critical', title: 'Insufficient reserve energy absorption', detail: `Reserve-energy deflection ${(o.deflection_reserve_m * 1e3).toFixed(0)} mm exceeds the ${(i.clearance * 1e3).toFixed(0)} mm clearance.`, action: 'Use larger or thicker cross tubes, a more ductile alloy, or more ground clearance; for payload protection pass the residual velocity to Suite 15.', basis: 'Reserve-energy drop (1.5 × limit energy)' });
    if (o.permanent_set_m > 0) out.push({ severity: 'warn', title: 'Limit landing causes permanent set', detail: `${(o.permanent_set_m * 1e3).toFixed(1)} mm of permanent deflection.`, action: `Keep the limit landing elastic: yielding starts at about ${o.sink_first_yield_ms.toFixed(2)} m/s sink. Define a hard-landing inspection threshold from this value.`, basis: 'Elastic limit of the cross tubes' });
    else out.push({ severity: 'info', title: 'Limit landing is elastic', detail: `Peak reaction ${o.gear_load_factor.toFixed(2)} × weight with no permanent set; yielding begins near ${o.sink_first_yield_ms.toFixed(2)} m/s.`, action: 'A stiffer gear than necessary raises landing loads on the airframe and payload; trim tube size to balance load factor against deflection and mass.', basis: 'Load factor versus deflection trade' });
    return out;
  },
};

export default {
  id: 'gear', n: 14,
  tagline: 'How hard the aircraft lands, how far it takes to stop, and whether it stays stable and within limits on the ground.',
  analyses: [drop, brake, ground, taxi, shimmy, skid],
  consumes: [{ from: 'performance', keys: ['V_stall_ms'], why: 'Touchdown speed for braking, crosswind and shimmy conditions' }],
  provides: [
    { key: 'gear_load_N', label: 'Peak gear load', unit: 'N' }, { key: 'gear_load_factor', label: 'Ground reaction factor', unit: '-' }, { key: 'stroke_used_m', label: 'Stroke used', unit: 'm' },
    { key: 'stop_dist_m', label: 'Braked ground roll', unit: 'm' }, { key: 'brake_temp_K', label: 'Brake temperature', unit: 'K' },
  ],
  handoff: [
    { model: 'Full multibody landing-gear model coupled to six-degree-of-freedom aircraft dynamics', why: 'Needs articulated gear kinematics, spin-up/spring-back, pitch and roll coupling and asymmetric landings; here one leg is simulated vertically', tool: 'Multibody dynamics software' },
    { model: 'Coupled multibody–finite-element landing impact with flexible airframe', why: 'Dynamic landing loads on a flexible wing and fuselage need modal or FE structural models', tool: 'FE / multibody co-simulation; Suite 2 for static stressing with the loads from this suite' },
    { model: 'Detailed tyre contact models (Hertzian contact patch, deformable tyre FE, hydroplaning fluid film)', why: 'The tyre is a linear spring with a generic friction–slip curve', tool: 'Tyre test data fitted to semi-empirical tyre models, explicit FE for tyre impact' },
    { model: 'Brake thermal–mechanical model (disc stack temperature field, fade, wear, cooling)', why: 'A single lumped heat-sink temperature is computed', tool: 'Brake dynamometer data and thermal FE; Suite 12 for transient conduction' },
    { model: 'Nonlinear shimmy with free play, friction and gear lateral–torsional flexibility', why: 'The native model is linear with a rigid leg', tool: 'Multibody shimmy analysis with measured tyre data' },
    { model: 'Hard-landing and crash impact beyond the gear stroke (structural collapse)', why: 'Explicit nonlinear structural impact is outside the gear model', tool: 'Suite 15 (crashworthiness) and explicit FE codes' },
    { model: 'Pavement strength classification (ACN/PCN, ACR/PCR)', why: 'Requires the standardised layered-pavement procedure', tool: 'ICAO pavement design software' },
  ],
};
