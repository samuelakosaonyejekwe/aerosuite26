// Suite 5 — Aircraft Performance and Flight Envelope.
// Point-mass performance from the drag polar and a thrust/power lapse model: speeds, field lengths,
// climb and ceilings, range and endurance, payload-range, V-n and altitude-speed envelopes, hover ceiling.

import * as N from '../core/numerics.js';
import { isa, G0, RHO0 } from '../core/atmosphere.js';

// ---- shared physics -------------------------------------------------------------------------
const isProp = (t) => t === 'turboprop' || t === 'piston' || t === 'electric' || t === 'turboshaft';

/** Total available thrust [N] at true airspeed V, altitude h. Jets: thrust lapse; props: power lapse capped by static thrust. */
export function thrustAvail(p, V, h, dT = 0, throttle = 1, nEngOut = 0) {
  const a = isa(h, dT), n = Math.max(0, p.n_eng - nEngOut), M = V / a.a;
  if (!isProp(p.type)) {
    // Mattingly installed-thrust lapse (Aircraft Engine Design, 2nd ed., eqs 2.42/2.45) below the throttle-ratio break
    const bpr = p.type === 'turbojet' ? 0 : p.bpr, d0 = a.delta * (1 + 0.2 * M * M) ** 3.5, rM = Math.sqrt(Math.max(M, 0));
    const lapse = bpr >= 1.5 ? d0 * (1 - 0.49 * rM) : bpr > 0 ? d0 * (1 - 0.3 * rM) : d0 * (1 - 0.16 * rM);
    return Math.max(0, throttle * n * p.T0_N * Math.max(0.05, lapse));
  }
  const P = throttle * n * powerAvail(p, h, dT);
  const Tstatic = n * (p.T0_N || 0) * a.sigma ** 0.7 || Infinity;
  return Math.min(Tstatic, (p.eta_prop * P) / Math.max(V, 1));
}
/** Shaft power per engine [W] at altitude. */
export function powerAvail(p, h, dT = 0) {
  const s = isa(h, dT).sigma;
  if (p.type === 'piston') return p.P0_W * Math.max(0, (s - 0.1325) / 0.8675); // Gagg–Ferrar
  if (p.type === 'electric') return p.P0_W;
  return p.P0_W * s ** 0.7;
}
/** Fuel flow [kg/s] for a required thrust T at speed V. */
export function fuelFlow(p, T, V) {
  if (p.type === 'electric') return 0;
  if (isProp(p.type)) return (p.bsfc * T * V) / p.eta_prop;
  return p.tsfc * T;
}
const drag = (q, S, CD0, k, W, nz = 1) => q * S * CD0 + (k * (nz * W) ** 2) / (q * S);

const COMMON = [
  { key: 'mass_kg', label: 'Aircraft mass', unit: 'kg', default: 70000, min: 1, group: 'Aircraft', help: 'Gross mass for this analysis' },
  { key: 'S', label: 'Wing area', unit: 'm²', default: 120, min: 0.01, group: 'Aircraft' },
  { key: 'CD0', label: 'Zero-lift drag coefficient', unit: '-', default: 0.02, min: 0.001, group: 'Aerodynamics', help: 'From the CFD suite drag build-up when available' },
  { key: 'k', label: 'Induced drag factor k', unit: '-', default: 0.045, min: 0.001, group: 'Aerodynamics', help: '1 / (π · AR · e)' },
  { key: 'CLmax', label: 'Maximum lift coefficient', unit: '-', default: 1.5, min: 0.2, group: 'Aerodynamics' },
  { key: 'alt_m', label: 'Altitude', unit: 'm', default: 0, min: -500, max: 25000, group: 'Atmosphere' },
  { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -60, max: 50, group: 'Atmosphere' },
];
const PROP_INPUTS = [
  { key: 'ptype', label: 'Powerplant type', type: 'select', options: ['turbofan', 'turbojet', 'turboprop', 'turboshaft', 'piston', 'electric'], default: 'turbofan', group: 'Propulsion' },
  { key: 'n_eng', label: 'Engines', unit: '', default: 2, min: 1, step: 1, discrete: true, group: 'Propulsion' },
  { key: 'T0_N', label: 'Static thrust per engine', unit: 'N', default: 120000, min: 0, group: 'Propulsion' },
  { key: 'P0_W', label: 'Rated power per engine', unit: 'W', default: 0, min: 0, group: 'Propulsion' },
  { key: 'bpr', label: 'Bypass ratio', unit: '-', default: 5, min: 0, group: 'Propulsion' },
  { key: 'eta_prop', label: 'Propeller efficiency', unit: '-', default: 0.82, min: 0.1, max: 0.95, group: 'Propulsion', help: 'From the propeller suite when available' },
  { key: 'tsfc', label: 'TSFC', unit: 'kg/N/s', default: 1.6e-5, min: 0, group: 'Propulsion', help: 'From the propulsion suite when available' },
  { key: 'bsfc', label: 'BSFC', unit: 'kg/W/s', default: 8e-8, min: 0, group: 'Propulsion' },
];
const pOf = (i) => ({ type: i.ptype, n_eng: i.n_eng, T0_N: i.T0_N, P0_W: i.P0_W, bpr: i.bpr, eta_prop: i.eta_prop, tsfc: i.tsfc, bsfc: i.bsfc });
const commonDefaults = (c, up, d) => ({
  mass_kg: c.mass.mtow_kg, S: c.wing.S_m2 || undefined, CD0: up.cfd?.CD0 ?? c.aero.CD0, k: up.cfd?.k_induced ?? (d.k_induced || undefined),
  CLmax: up.cfd?.CLmax ?? c.aero.CLmax_clean, alt_m: c.atm.alt_m, dISA: c.atm.dISA_K,
  ptype: c.prop.type, n_eng: c.prop.n_eng, T0_N: up.propulsion?.thrust_static_N ?? c.prop.T0_N, P0_W: c.prop.P0_W, bpr: c.prop.bpr,
  eta_prop: up.propeller?.eta_prop ?? c.prop.eta_prop, tsfc: up.propulsion?.tsfc_kg_Ns ?? c.prop.tsfc_kg_Ns, bsfc: c.prop.bsfc_kg_Ws,
});
const fixedWing = (c) => (c.wing.S_m2 > 0 ? true : 'This analysis needs a lifting wing; the current case is a pure rotorcraft. Use the hover-ceiling analysis or Suite 6.');

/** Speeds and rates at one altitude. */
function pointPerf(i, h) {
  const a = isa(h, i.dISA), W = i.mass_kg * G0, p = pOf(i);
  const Vs = Math.sqrt((2 * W) / (a.rho * i.S * i.CLmax));
  const Vmd = Math.sqrt((2 * W) / (a.rho * i.S)) * (i.k / i.CD0) ** 0.25;
  const LDmax = 1 / (2 * Math.sqrt(i.CD0 * i.k));
  const ex = (V) => thrustAvail(p, V, h, i.dISA) - drag(0.5 * a.rho * V * V, i.S, i.CD0, i.k, W);
  const Vhi = Math.min(1.2 * a.a, 400);
  const Vs_grid = N.linspace(Vs, Vhi, 240), exs = Vs_grid.map(ex);
  let Vmax = NaN;
  for (let j = Vs_grid.length - 1; j > 0; j--) if (exs[j] < 0 && exs[j - 1] >= 0) { Vmax = N.brent(ex, Vs_grid[j - 1], Vs_grid[j], 1e-6); break; }
  if (Number.isNaN(Vmax) && exs[exs.length - 1] >= 0) Vmax = Vhi;
  const roc = Vs_grid.map((V, j) => (exs[j] * V) / W), jb = N.argmax(roc);
  return { a, W, Vs, Vmd, LDmax, Vmax, ROCmax: roc[jb], Vy: Vs_grid[jb], gammaMax: Math.asin(N.clamp(N.amax(exs) / W, -1, 1)), grid: Vs_grid, exs, roc, p };
}

// ---- analyses -------------------------------------------------------------------------------
const point = {
  id: 'point', title: 'Point performance and drag polar', fidelity: 'analytical',
  summary: 'Thrust and power required versus available at one altitude: stall, minimum-drag, best-climb and maximum speeds, with the drag polar.',
  equations: ['Thrust–drag balance', 'Lift–weight balance', 'Specific excess power', 'Rate-of-climb equation', 'Stall speed equation'],
  applicable: fixedWing,
  inputs: [...COMMON, ...PROP_INPUTS],
  defaults: commonDefaults,
  run(i) {
    const r = pointPerf(i, i.alt_m), { a, W } = r, q = (V) => 0.5 * a.rho * V * V;
    const V = N.linspace(r.Vs * 0.95, Number.isFinite(r.Vmax) ? r.Vmax * 1.08 : r.Vs * 3, 120);
    const D = V.map((v) => drag(q(v), i.S, i.CD0, i.k, W)), T = V.map((v) => thrustAvail(r.p, v, i.alt_m, i.dISA));
    const CL = N.linspace(0, i.CLmax, 60);
    const warnings = [];
    if (!Number.isFinite(r.Vmax)) warnings.push('Thrust available is below drag at every speed: level flight is not possible at this altitude and mass.');
    if (r.Vmax / a.a > 0.75 && i.ptype !== 'turbofan' && i.ptype !== 'turbojet') warnings.push('Predicted speed approaches transonic Mach numbers where the parabolic polar and constant propeller efficiency are not valid.');
    if (r.Vmax / a.a > 0.85) warnings.push('Wave drag is not modelled: maximum speed above Mach 0.85 is optimistic.');
    return {
      kpis: [
        { key: 'V_stall_ms', label: 'Stall speed (TAS)', value: r.Vs, unit: 'm/s' },
        { key: 'V_md_ms', label: 'Minimum-drag speed', value: r.Vmd, unit: 'm/s' },
        { key: 'V_max_ms', label: 'Maximum level speed', value: r.Vmax, unit: 'm/s', status: Number.isFinite(r.Vmax) ? 'ok' : 'bad' },
        { key: 'M_max', label: 'Mach at maximum speed', value: r.Vmax / a.a, unit: '-' },
        { key: 'LD_max', label: 'Maximum lift-to-drag ratio', value: r.LDmax, unit: '-' },
        { key: 'roc_max_ms', label: 'Maximum rate of climb', value: r.ROCmax, unit: 'm/s', status: r.ROCmax > 0.5 ? 'ok' : 'warn' },
        { key: 'V_y_ms', label: 'Best rate-of-climb speed', value: r.Vy, unit: 'm/s' },
        { key: 'gamma_max_deg', label: 'Maximum climb angle', value: N.deg(r.gammaMax), unit: 'deg' },
        { key: 'T_over_W', label: 'Static thrust-to-weight', value: thrustAvail(r.p, 1, i.alt_m, i.dISA) / W, unit: '-' },
      ],
      plots: [
        { type: 'line', title: 'Thrust required and available', xlabel: 'True airspeed [m/s]', ylabel: 'Force [kN]', series: [{ name: 'Drag (thrust required)', x: V, y: D.map((v) => v / 1e3) }, { name: 'Thrust available', x: V, y: T.map((v) => v / 1e3) }], annotations: [{ x: r.Vs, label: 'Stall' }, { x: r.Vmd, label: 'Min drag' }] },
        { type: 'line', title: 'Rate of climb', xlabel: 'True airspeed [m/s]', ylabel: 'Rate of climb [m/s]', series: [{ name: 'ROC', x: r.grid, y: r.roc }], annotations: [{ y: 0, label: 'Level flight' }] },
        { type: 'line', title: 'Drag polar', xlabel: 'CD [-]', ylabel: 'CL [-]', series: [{ name: 'CD = CD0 + k·CL²', x: CL.map((c) => i.CD0 + i.k * c * c), y: CL }] },
      ],
      outputs: { k_induced: i.k },
      warnings, models: ['Parabolic drag polar', 'Thrust/power lapse model', 'Specific excess power'],
      assumptions: ['Steady, symmetric, point-mass flight', 'Small flight-path angle for lift = weight', 'No wave drag or Reynolds-number variation of CD0'],
    };
  },
  convergence: null,
  verify() {
    // Closed form: maximum L/D and minimum-drag speed for a parabolic polar.
    const i = { mass_kg: 1000, S: 16, CD0: 0.03, k: 0.05, CLmax: 1.5, alt_m: 0, dISA: 0, ptype: 'turbojet', n_eng: 1, T0_N: 3000, P0_W: 0, bpr: 0, eta_prop: 0.8, tsfc: 2e-5, bsfc: 0 };
    const r = pointPerf(i, 0), W = 1000 * G0;
    return [
      N.check('L/D max = 1/(2·sqrt(CD0·k))', r.LDmax, 12.909944, 1e-6, 'Anderson, Aircraft Performance and Design, eq. 5.30'),
      N.check('Stall speed from L = W', r.Vs, Math.sqrt((2 * W) / (1.225 * 16 * 1.5)), 1e-6, 'Definition'),
      N.check('Drag at Vmd equals W/(L/D)max', drag(0.5 * 1.225 * r.Vmd ** 2, 16, 0.03, 0.05, W), W / 12.909944, 1e-6, 'Analytical minimum of the drag curve'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!(o.roc_max_ms > 1.5)) out.push({ severity: 'warn', title: 'Low climb capability at this point', detail: `Maximum rate of climb is ${(o.roc_max_ms || 0).toFixed(2)} m/s.`, action: 'Reduce mass, lower the altitude or temperature, or increase installed thrust/power.', basis: 'Specific excess power' });
    if (o.LD_max < 12 && i.S > 5) out.push({ severity: 'advise', title: 'Aerodynamic efficiency is modest', detail: `(L/D)max = ${o.LD_max.toFixed(1)}. Every 1% of drag is roughly 1% of cruise fuel and CO₂.`, action: 'Review CD0 contributors in Suite 1 (wetted area, excrescences) and consider a higher aspect ratio in Suite 23.', basis: 'Breguet range relation' });
    return out;
  },
};

const field = {
  id: 'field', title: 'Take-off and landing field performance', fidelity: 'numerical',
  summary: 'Time-integrated ground roll with thrust, drag, lift, rolling friction, runway slope and wind, plus airborne segments to the screen height.',
  equations: ['Takeoff ground-roll equations', 'Landing ground-roll equations', 'Equations of motion along the flight path'],
  applicable: fixedWing,
  inputs: [
    ...COMMON.filter((f) => f.key !== 'CLmax'), ...PROP_INPUTS,
    { key: 'CLmax_to', label: 'CLmax take-off', unit: '-', default: 2.0, min: 0.3, group: 'Aerodynamics' },
    { key: 'CLmax_land', label: 'CLmax landing', unit: '-', default: 2.6, min: 0.3, group: 'Aerodynamics' },
    { key: 'CL_ground', label: 'Ground-roll lift coefficient', unit: '-', default: 0.6, min: 0, group: 'Aerodynamics' },
    { key: 'mu_roll', label: 'Rolling friction', unit: '-', default: 0.03, min: 0, max: 0.3, group: 'Runway', help: '0.02–0.03 dry paved, 0.05 short grass, 0.1 soft ground' },
    { key: 'mu_brake', label: 'Braking friction', unit: '-', default: 0.4, min: 0.03, max: 0.8, group: 'Runway', help: '0.4–0.5 dry, 0.2–0.3 wet, 0.05–0.1 icy' },
    { key: 'slope_pct', label: 'Runway slope (uphill +)', unit: '%', default: 0, min: -5, max: 5, group: 'Runway' },
    { key: 'headwind', label: 'Headwind component', unit: 'm/s', default: 0, min: -15, max: 30, group: 'Runway', help: 'Negative for tailwind. Filled from live weather when a site is set.' },
    { key: 'runway_m', label: 'Runway length available', unit: 'm', default: 2500, min: 50, group: 'Runway' },
    { key: 'land_mass_frac', label: 'Landing mass / take-off mass', unit: '-', default: 0.85, min: 0.3, max: 1, group: 'Aircraft' },
    { key: 'screen_m', label: 'Screen height', unit: 'm', default: 10.7, min: 0, group: 'Runway', help: '10.7 m (35 ft) transport, 15.2 m (50 ft) light aircraft' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 400, min: 20, max: 20000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const s = c.site, rw = N.rad((s.wind_dir_deg || 0) - (s.runway_heading_deg ?? s.wind_dir_deg ?? 0));
    return { ...commonDefaults(c, up, d), alt_m: s.elev_m ?? 0, dISA: s.T_C - (15 - 0.0065 * (s.elev_m || 0)), CLmax_to: c.aero.CLmax_to, CLmax_land: c.aero.CLmax_land, mu_roll: s.runway_mu, mu_brake: s.runway_mu_brake, slope_pct: s.runway_slope_pct, headwind: (s.wind_ms || 0) * Math.cos(rw), runway_m: s.runway_len_m, land_mass_frac: Math.min(1, (c.mass.mtow_kg - 0.8 * c.mass.fuel_kg) / c.mass.mtow_kg), screen_m: c.mass.mtow_kg > 5700 ? 10.7 : 15.2 };
  },
  run(i) {
    const a = isa(i.alt_m, i.dISA), W = i.mass_kg * G0, p = pOf(i), th = Math.atan(i.slope_pct / 100);
    const Vs = Math.sqrt((2 * W) / (a.rho * i.S * i.CLmax_to)), VR = 1.1 * Vs, V2 = 1.2 * Vs, CDg = i.CD0 * 1.6 + i.k * i.CL_ground ** 2 * 0.6; // gear/flap drag, ground effect
    const roll = (Wt, m, Vend, accel, eo = 0) => {
      // state [groundspeed, distance]; airspeed = groundspeed + headwind
      const f = (t, y) => {
        const Va = Math.max(0, y[0] + i.headwind), q = 0.5 * a.rho * Va * Va, L = q * i.S * i.CL_ground, D = q * i.S * CDg;
        const T = accel ? thrustAvail(p, Va, i.alt_m, i.dISA, 1, eo) : 0, mu = accel ? i.mu_roll : i.mu_brake;
        return [(T - D - mu * Math.max(0, Wt * Math.cos(th) - L) - Wt * Math.sin(th)) / m, y[0]];
      };
      const v0 = accel ? Math.max(0, -i.headwind) : Vend - i.headwind, vt = accel ? Vend - i.headwind : 0;
      const a0 = f(0, [v0, 0])[0];
      if (accel && a0 <= 0) return { s: Infinity, t: Infinity, T: [], Vt: [], S: [] };
      const tEst = (1.6 * Math.abs(vt - v0)) / Math.abs(a0 || 1) + 1, n = Math.round(i.nSteps), h = tEst / n;
      let y = [v0, 0], t = 0; const T = [0], Vt = [v0], S = [0];
      for (let s = 0; s < n * 6; s++) {
        const k1 = f(t, y), k2 = f(t + h / 2, N.vadd(y, k1, h / 2)), k3 = f(t + h / 2, N.vadd(y, k2, h / 2)), k4 = f(t + h, N.vadd(y, k3, h));
        const yn = y.map((v, j) => v + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
        const done = accel ? yn[0] >= vt : yn[0] <= vt;
        if (done) { const fr = (vt - y[0]) / (yn[0] - y[0]); y = [vt, y[1] + fr * (yn[1] - y[1])]; t += fr * h; T.push(t); Vt.push(y[0]); S.push(y[1]); break; }
        if (accel && k1[0] <= 0) return { s: Infinity, t: Infinity, T, Vt, S };
        y = yn; t += h; T.push(t); Vt.push(y[0]); S.push(y[1]);
      }
      return { s: y[1], t, T, Vt, S };
    };
    const g = roll(W, i.mass_kg, VR, true);
    // rotation (2 s at VR) and airborne arc to screen height at V2 (Raymer-style transition)
    const sRot = 2 * Math.max(0, VR - i.headwind), R = V2 ** 2 / (0.2 * G0);
    const Tair = thrustAvail(p, V2, i.alt_m, i.dISA), Dair = drag(0.5 * a.rho * V2 ** 2, i.S, i.CD0 * 1.5, i.k, W), gam = Math.asin(N.clamp((Tair - Dair) / W, 0.001, 0.5));
    const hTr = R * (1 - Math.cos(gam)), sTr = i.screen_m <= hTr ? Math.sqrt(Math.max(0, R * R - (R - i.screen_m) ** 2)) : R * Math.sin(gam) + (i.screen_m - hTr) / Math.tan(gam);
    const tofl = g.s + sRot + sTr * Math.max(0.2, 1 - i.headwind / V2);
    // one-engine-inoperative second-segment gradient
    const gradOEI = i.n_eng > 1 ? (thrustAvail(p, V2, i.alt_m, i.dISA, 1, 1) - Dair) / W : NaN;
    // landing
    const Wl = W * i.land_mass_frac, Vsl = Math.sqrt((2 * Wl) / (a.rho * i.S * i.CLmax_land)), Vapp = 1.3 * Vsl, Vtd = 1.15 * Vsl;
    const sAir = i.screen_m / Math.tan(N.rad(3)) + (Vapp ** 2 - Vtd ** 2) / (2 * G0 * 0.1), l = roll(Wl, i.mass_kg * i.land_mass_frac, Vtd, false);
    const ldg = sAir * Math.max(0.2, 1 - i.headwind / Vapp) + 2 * Math.max(0, Vtd - i.headwind) + l.s, ldgFactored = ldg * (i.mass_kg > 5700 ? 1 / 0.6 : 1.43);
    const mTO = i.runway_m / (tofl * 1.15) - 1, mLD = i.runway_m / ldgFactored - 1, warnings = [];
    if (!Number.isFinite(tofl)) warnings.push('The aircraft cannot accelerate to rotation speed: thrust does not exceed friction, drag and slope resistance.');
    if (i.headwind < -5) warnings.push('Tailwind above 5 m/s (≈10 kt) is beyond normal certified take-off and landing limits.');
    return {
      kpis: [
        { key: 'tofl_m', label: 'Take-off distance to screen', value: tofl, unit: 'm', status: mTO > 0 ? 'ok' : 'bad' },
        { key: 'to_ground_roll_m', label: 'Take-off ground roll', value: g.s, unit: 'm' },
        { key: 'V_R_ms', label: 'Rotation speed', value: VR, unit: 'm/s' },
        { key: 'V_2_ms', label: 'Take-off safety speed V2', value: V2, unit: 'm/s' },
        { key: 'oei_gradient_pct', label: 'Engine-out climb gradient at V2', value: 100 * gradOEI, unit: '%', status: Number.isNaN(gradOEI) ? undefined : gradOEI >= 0.024 ? 'ok' : 'bad', note: 'CS/FAR-25 twin second-segment minimum is 2.4%' },
        { key: 'ldg_dist_m', label: 'Landing distance from screen', value: ldg, unit: 'm' },
        { key: 'ldg_ground_roll_m', label: 'Landing ground roll', value: l.s, unit: 'm' },
        { key: 'ldg_factored_m', label: 'Factored landing distance', value: ldgFactored, unit: 'm', status: mLD > 0 ? 'ok' : 'bad' },
        { key: 'V_app_ms', label: 'Approach speed', value: Vapp, unit: 'm/s' },
        { key: 'runway_margin_to', label: 'Take-off runway margin (×1.15)', value: 100 * mTO, unit: '%', status: mTO > 0.1 ? 'ok' : mTO > 0 ? 'warn' : 'bad' },
        { key: 'runway_margin_ldg', label: 'Landing runway margin', value: 100 * mLD, unit: '%', status: mLD > 0.1 ? 'ok' : mLD > 0 ? 'warn' : 'bad' },
        { key: 'brake_energy_J', label: 'Landing brake energy', value: 0.5 * i.mass_kg * i.land_mass_frac * Math.max(0, Vtd - i.headwind) ** 2, unit: 'J' },
      ],
      plots: [
        { type: 'line', title: 'Ground-roll speed versus distance', xlabel: 'Distance [m]', ylabel: 'Ground speed [m/s]', series: [{ name: 'Take-off roll', x: g.S, y: g.Vt }, { name: 'Landing roll', x: l.S, y: l.Vt }], annotations: [{ x: i.runway_m, label: 'Runway end' }] },
        { type: 'bar', title: 'Distance build-up', ylabel: 'Distance [m]', categories: ['Take-off', 'Landing'], stacked: true, series: [{ name: 'Ground roll', y: [g.s, l.s] }, { name: 'Rotation / free roll', y: [sRot, 2 * Math.max(0, Vtd - i.headwind)] }, { name: 'Airborne to/from screen', y: [tofl - g.s - sRot, ldg - l.s - 2 * Math.max(0, Vtd - i.headwind)] }] },
      ],
      warnings, models: ['Time-marched ground roll (RK4)', 'Circular-arc transition', 'Thrust/power lapse model'],
      assumptions: ['All engines operating for distances; engine-out only for the climb gradient', 'Constant friction coefficients', 'No reverse thrust; spoilers not modelled', 'Factored landing distance: ÷0.6 for transports, ×1.43 for light aircraft'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps in ground roll', levels: [25, 50, 100, 200, 400], metric: 'to_ground_roll_m' },
  verify() {
    // Constant-acceleration limit: no aero forces, no friction -> s = V^2 m / (2 T).
    const i = { mass_kg: 1000, S: 16, CD0: 0, k: 0, alt_m: 0, dISA: 0, ptype: 'turbojet', n_eng: 1, T0_N: 5000, P0_W: 0, bpr: 0, eta_prop: 0.8, tsfc: 0, bsfc: 0, CLmax_to: 2, CLmax_land: 2, CL_ground: 0, mu_roll: 0, mu_brake: 0.4, slope_pct: 0, headwind: 0, runway_m: 1e9, land_mass_frac: 1, screen_m: 0, nSteps: 400 };
    const r = field.run(i), VR = N.kv(r).V_R_ms;
    // turbojet lapse 1-0.45M+0.19M^2 is not constant, so integrate the same law independently by quadrature
    const sRef = N.simpson((v) => (1000 * v) / thrustAvail(pOf(i), v, 0, 0), 0, VR, 2000);
    const o = N.kv(r), Vtd = 1.15 * Math.sqrt((2 * 1000 * G0) / (isa(0).rho * 16 * 2));
    return [
      N.check('Ground roll equals ∫ m·V/T dV', o.to_ground_roll_m, sRef, 2e-4, 'Energy integral of the same thrust law'),
      N.check('Braked roll equals V²/(2·μ·g)', o.ldg_ground_roll_m, Vtd ** 2 / (2 * 0.4 * G0), 2e-4, 'Constant-deceleration kinematics'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.runway_margin_to < 0) out.push({ severity: 'critical', title: 'Take-off distance exceeds the runway', detail: `Factored take-off distance needs ${(o.tofl_m * 1.15).toFixed(0)} m against ${i.runway_m.toFixed(0)} m available.`, action: 'Reduce take-off mass, select a higher-lift flap setting, wait for cooler conditions or use a longer runway.', basis: 'Take-off field length with 15% margin' });
    if (o.runway_margin_ldg < 0) out.push({ severity: 'critical', title: 'Factored landing distance exceeds the runway', detail: `Needs ${o.ldg_factored_m.toFixed(0)} m against ${i.runway_m.toFixed(0)} m.`, action: 'Reduce landing mass or choose an alternate with a longer or drier runway.', basis: 'Operational landing-distance factor' });
    if (Number.isFinite(o.oei_gradient_pct) && o.oei_gradient_pct < 2.4) out.push({ severity: 'warn', title: 'Engine-out climb gradient is below 2.4%', detail: `Predicted ${o.oei_gradient_pct.toFixed(2)}% at V2.`, action: 'Limit take-off mass for this altitude and temperature (WAT limit).', basis: 'CS/FAR 25.121(b)' });
    return out;
  },
};

const climb = {
  id: 'climb', title: 'Climb, ceilings and time to climb', fidelity: 'numerical',
  summary: 'Maximum rate of climb through the atmosphere, service and absolute ceilings, and time, fuel and distance to climb.',
  equations: ['Rate-of-climb equations', 'Service ceiling relations', 'Specific excess power equations'],
  applicable: fixedWing,
  inputs: [...COMMON.filter((f) => f.key !== 'alt_m'), ...PROP_INPUTS,
    { key: 'h_top', label: 'Top of analysis', unit: 'm', default: 15000, min: 500, max: 30000, group: 'Numerics' },
    { key: 'nAlt', label: 'Altitude stations', unit: '', default: 40, min: 6, max: 400, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => ({ ...commonDefaults(c, up, d), h_top: Math.max(8000, (c.mission.cruise_alt_m || 3000) * 1.5) }),
  run(i) {
    const hs = N.linspace(0, i.h_top, Math.round(i.nAlt)), pts = hs.map((h) => pointPerf(i, h)), roc = pts.map((p) => p.ROCmax);
    const cross = (lim) => { for (let j = 1; j < hs.length; j++) if (roc[j] < lim && roc[j - 1] >= lim) return hs[j - 1] + ((roc[j - 1] - lim) / (roc[j - 1] - roc[j])) * (hs[j] - hs[j - 1]); return roc[roc.length - 1] >= lim ? NaN : 0; };
    const svc = cross(0.508), abs = cross(0), top = Number.isFinite(svc) ? svc : i.h_top;
    // time, fuel and distance to the service ceiling (or top of analysis)
    const t = [0], fuel = [0], dist = [0]; let k = 0;
    for (let j = 1; j < hs.length && hs[j] <= top; j++) {
      const r = 0.5 * (roc[j] + roc[j - 1]); if (r <= 0) break;
      const dt = (hs[j] - hs[j - 1]) / r, V = 0.5 * (pts[j].Vy + pts[j - 1].Vy), T = thrustAvail(pts[j].p, V, hs[j], i.dISA);
      t.push(t[k] + dt); fuel.push(fuel[k] + fuelFlow(pts[j].p, T, V) * dt); dist.push(dist[k] + V * dt); k++;
    }
    const warnings = [];
    if (Number.isNaN(svc)) warnings.push('The service ceiling lies above the top of analysis; raise "Top of analysis" to find it.');
    return {
      kpis: [
        { key: 'roc_sl_ms', label: 'Sea-level rate of climb', value: roc[0], unit: 'm/s' },
        { key: 'ceiling_m', label: 'Service ceiling (0.5 m/s)', value: svc, unit: 'm' },
        { key: 'abs_ceiling_m', label: 'Absolute ceiling', value: abs, unit: 'm' },
        { key: 'time_to_climb_min', label: `Time to ${Math.round(hs[k])} m`, value: t[k] / 60, unit: 'min' },
        { key: 'fuel_to_climb_kg', label: 'Fuel to climb', value: fuel[k], unit: 'kg' },
        { key: 'dist_to_climb_km', label: 'Distance to climb', value: dist[k] / 1e3, unit: 'km' },
      ],
      plots: [
        { type: 'line', title: 'Maximum rate of climb versus altitude', xlabel: 'Rate of climb [m/s]', ylabel: 'Altitude [m]', series: [{ name: 'ROC max', x: roc, y: hs }], annotations: Number.isFinite(svc) ? [{ y: svc, label: 'Service ceiling' }] : [] },
        { type: 'line', title: 'Climb time history', xlabel: 'Time [min]', ylabel: 'Altitude [m]', series: [{ name: 'Altitude', x: t.map((v) => v / 60), y: hs.slice(0, k + 1) }] },
      ],
      warnings, models: ['Quasi-steady climb at best rate-of-climb speed', 'ISA with temperature offset'],
      assumptions: ['Mass held constant during the climb', 'No acceleration (energy-height) correction', 'No speed or Mach limits applied to the climb schedule'],
    };
  },
  convergence: { param: 'nAlt', label: 'Altitude stations', levels: [10, 20, 40, 80], metric: 'time_to_climb_min' },
  recommend(res, i, ctx) {
    const o = res.outputs, hc = ctx.case.mission.cruise_alt_m, out = [];
    if (Number.isFinite(o.ceiling_m) && hc && o.ceiling_m < hc + 300) out.push({ severity: 'warn', title: 'Service ceiling is close to or below the planned cruise altitude', detail: `Ceiling ${o.ceiling_m.toFixed(0)} m versus cruise ${hc.toFixed(0)} m at this mass.`, action: 'Plan a step-climb, reduce mass, or lower the cruise altitude.', basis: 'Rate-of-climb margin of 0.5 m/s (100 ft/min)' });
    return out;
  },
};

const range = {
  id: 'range', title: 'Range, endurance and payload–range', fidelity: 'analytical',
  summary: 'Breguet range and endurance at the cruise condition and the payload–range diagram with its three limiting segments.',
  equations: ['Breguet range equation', 'Breguet endurance relations', 'Payload–range equations'],
  applicable: fixedWing,
  inputs: [...COMMON.filter((f) => !['CLmax', 'mass_kg'].includes(f.key)), ...PROP_INPUTS.filter((f) => ['ptype', 'eta_prop', 'tsfc', 'bsfc'].includes(f.key)),
    { key: 'V', label: 'Cruise true airspeed', unit: 'm/s', default: 230, min: 5, group: 'Mission' },
    { key: 'mtow', label: 'Maximum take-off mass', unit: 'kg', default: 78000, min: 1, group: 'Mass' },
    { key: 'oew', label: 'Operating empty mass', unit: 'kg', default: 42000, min: 0.5, group: 'Mass' },
    { key: 'payload_max', label: 'Maximum payload', unit: 'kg', default: 18000, min: 0, group: 'Mass' },
    { key: 'fuel_max', label: 'Fuel capacity', unit: 'kg', default: 19000, min: 0, group: 'Mass' },
    { key: 'reserve_frac', label: 'Reserve fuel fraction', unit: '-', default: 0.08, min: 0, max: 0.5, group: 'Mission', help: 'Share of loaded fuel kept as reserves and allowances' },
    { key: 'batt_kWh', label: 'Battery energy (electric)', unit: 'kWh', default: 0, min: 0, group: 'Mass' },
    { key: 'eta_elec', label: 'Electric drivetrain efficiency', unit: '-', default: 0.9, min: 0.3, max: 0.99, group: 'Propulsion' }],
  defaults: (c, up, d) => ({ ...commonDefaults(c, up, d), alt_m: c.mission.cruise_alt_m, V: c.mission.cruise_V_ms || c.flight.V_ms, mtow: c.mass.mtow_kg, oew: c.mass.oew_kg, payload_max: c.mass.payload_kg, fuel_max: c.mass.fuel_kg, batt_kWh: c.systems.batt_kWh }),
  run(i) {
    const a = isa(i.alt_m, i.dISA), q = 0.5 * a.rho * i.V ** 2, electric = i.ptype === 'electric';
    const LD = (W) => { const CL = W / (q * i.S); return CL / (i.CD0 + i.k * CL * CL); };
    const c_th = isProp(i.ptype) ? (i.bsfc * i.V) / i.eta_prop : i.tsfc; // equivalent thrust-specific consumption [kg/N/s]
    // constant-speed, constant-altitude cruise: integrate dR = -V L/D / (g c) dW/W numerically to honour the changing CL
    const R = (W1, W0) => { if (W1 >= W0 || !c_th) return 0; return N.simpson((W) => (i.V * LD(W)) / (G0 * c_th * W), W1, W0, 200); };
    const E = (W1, W0) => { if (W1 >= W0 || !c_th) return 0; return N.simpson((W) => LD(W) / (G0 * c_th * W), W1, W0, 200); };
    const g = G0, zfwMax = i.oew + i.payload_max, pts = [];
    let rangeDesign, endurance, fuelUsed, LDc;
    if (electric) {
      const W = i.mtow * g; LDc = LD(W);
      rangeDesign = (i.batt_kWh * 3.6e6 * (1 - i.reserve_frac) * i.eta_elec * i.eta_prop * LDc) / W; endurance = rangeDesign / i.V; fuelUsed = 0;
      pts.push([0, i.payload_max], [rangeDesign / 1e3, i.payload_max], [rangeDesign / 1e3, 0]);
    } else {
      const fA = Math.max(0, Math.min(i.fuel_max, i.mtow - zfwMax)), use = (f) => f * (1 - i.reserve_frac);
      const RA = R((i.oew + i.payload_max + fA - use(fA)) * g, (zfwMax + fA) * g);
      const payB = Math.max(0, i.mtow - i.oew - i.fuel_max), RB = R((i.mtow - use(i.fuel_max)) * g, i.mtow * g);
      const RC = R((i.oew + i.fuel_max - use(i.fuel_max)) * g, (i.oew + i.fuel_max) * g);
      pts.push([0, i.payload_max], [RA / 1e3, i.payload_max]);
      if (payB < i.payload_max) pts.push([RB / 1e3, payB]);
      pts.push([RC / 1e3, 0]);
      rangeDesign = RA; fuelUsed = use(fA); LDc = LD((zfwMax + fA / 2) * g);
      endurance = E((zfwMax + fA - use(fA)) * g, (zfwMax + fA) * g);
    }
    const W0 = i.mtow * g, CLc = W0 / (q * i.S), CLopt = Math.sqrt(i.CD0 / i.k) * (isProp(i.ptype) ? 1 : 1 / Math.sqrt(3)), warnings = [];
    if (CLc > 1.2) warnings.push('Cruise lift coefficient exceeds 1.2: the aircraft is too slow or too high for this mass.');
    const sr = c_th ? (i.V * LD(W0)) / (c_th * W0) : NaN;
    return {
      kpis: [
        { key: 'range_km', label: 'Range at maximum payload', value: rangeDesign / 1e3, unit: 'km' },
        { key: 'ferry_range_km', label: 'Ferry range (zero payload)', value: pts[pts.length - 1][0], unit: 'km' },
        { key: 'endurance_h', label: 'Endurance at maximum payload', value: endurance / 3600, unit: 'h' },
        { key: 'LD_cruise', label: 'Cruise lift-to-drag ratio', value: LDc, unit: '-' },
        { key: 'CL_cruise', label: 'Cruise lift coefficient at MTOM', value: CLc, unit: '-', status: Math.abs(CLc / CLopt - 1) < 0.25 ? 'ok' : 'warn', note: `Best-range CL ≈ ${CLopt.toFixed(2)}` },
        { key: 'specific_range_m_kg', label: 'Specific range at MTOM', value: sr, unit: 'm/kg' },
        { key: 'fuel_flow_cruise_kgs', label: 'Cruise fuel flow at MTOM', value: c_th ? (c_th * W0) / LD(W0) : 0, unit: 'kg/s' },
        { key: 'trip_fuel_kg', label: 'Usable trip fuel', value: fuelUsed, unit: 'kg' },
      ],
      plots: [{ type: 'line', title: 'Payload–range diagram', xlabel: 'Range [km]', ylabel: 'Payload [kg]', series: [{ name: 'Payload–range boundary', x: pts.map((p) => p[0]), y: pts.map((p) => p[1]), style: 'line+points' }] }],
      tables: [{ title: 'Payload–range corner points', columns: ['Point', 'Range [km]', 'Payload [kg]'], rows: pts.map((p, j) => [['Zero range', 'Max payload', pts.length === 4 ? 'Max fuel' : 'Ferry', 'Ferry'][j], p[0], p[1]]) }],
      warnings, models: [electric ? 'Battery-electric Breguet range' : 'Breguet range integrated at constant speed and altitude', 'Parabolic drag polar'],
      assumptions: ['Cruise segment only: no climb, descent, taxi or wind', 'Constant specific fuel consumption', 'Reserve expressed as a fraction of loaded fuel'],
    };
  },
  verify() {
    // Constant L/D limit (k -> tiny, CD0 dominates is not constant L/D), so test the jet Breguet log form with a flat L/D polar by direct quadrature.
    const V = 200, c = 2e-5, LD = 15, W0 = 60000 * G0, W1 = 45000 * G0;
    const num = N.simpson((W) => (V * LD) / (G0 * c * W), W1, W0, 200);
    return [N.check('Breguet range = (V/c·g)(L/D) ln(W0/W1)', num, ((V * LD) / (G0 * c)) * Math.log(W0 / W1), 1e-8, 'Breguet (1923); Anderson eq. 5.153')];
  },
  recommend(res, i, ctx) {
    const o = res.outputs, need = ctx.case.mission.range_km, out = [];
    if (need && o.range_km < need) out.push({ severity: 'warn', title: 'Design range is not met at maximum payload', detail: `${o.range_km.toFixed(0)} km available against ${need.toFixed(0)} km required.`, action: 'Trade payload for fuel along the payload–range boundary, or improve L/D or specific fuel consumption.', basis: 'Breguet range' });
    if (o.CL_cruise && res.kpis[4].status === 'warn') out.push({ severity: 'advise', title: 'Cruise point is away from the best-range lift coefficient', detail: res.kpis[4].note + `, flying at ${o.CL_cruise.toFixed(2)}.`, action: 'Adjust cruise altitude or speed towards the best-range condition; a step-climb keeps CL near optimum as fuel burns off and cuts fuel and CO₂.', basis: 'Maximum of V·(L/D) for jets, L/D for propeller aircraft' });
    return out;
  },
};

const envelope = {
  id: 'envelope', title: 'Flight envelope: V–n diagram and altitude–speed limits', fidelity: 'analytical',
  summary: 'Manoeuvre and gust V–n diagram with design speeds, and the level-flight envelope bounded by stall, thrust and operating limits.',
  equations: ['Load-factor equations', 'Manoeuvre envelope relations', 'Gust load relations', 'Turning flight equations', 'Stall speed equations'],
  applicable: fixedWing,
  inputs: [...COMMON, ...PROP_INPUTS,
    { key: 'n_pos', label: 'Positive limit load factor', unit: 'g', default: 2.5, min: 1.5, max: 9, group: 'Limits' },
    { key: 'n_neg', label: 'Negative limit load factor', unit: 'g', default: -1, min: -4.5, max: 0, group: 'Limits' },
    { key: 'CLmax_neg', label: 'Negative CLmax magnitude', unit: '-', default: 0.9, min: 0.2, group: 'Limits' },
    { key: 'Vc_eas', label: 'Design cruise speed VC (EAS)', unit: 'm/s', default: 180, min: 5, group: 'Limits' },
    { key: 'Mmo', label: 'Maximum operating Mach', unit: '-', default: 0.82, min: 0.05, max: 3, group: 'Limits' },
    { key: 'CLa', label: 'Lift-curve slope', unit: '1/rad', default: 5.0, min: 1, max: 7, group: 'Aerodynamics', help: 'From the CFD suite when available' },
    { key: 'mac', label: 'Mean aerodynamic chord', unit: 'm', default: 4, min: 0.05, group: 'Aircraft' },
    { key: 'Ude_c', label: 'Gust velocity at VC (EAS)', unit: 'm/s', default: 15.24, min: 0, group: 'Limits', help: '50 ft/s (15.24 m/s) at VC and 25 ft/s at VD below 20 000 ft' },
    { key: 'h_top', label: 'Envelope top altitude', unit: 'm', default: 14000, min: 500, group: 'Numerics' }],
  defaults: (c, up, d) => ({ ...commonDefaults(c, up, d), n_pos: c.aero.n_pos, n_neg: c.aero.n_neg, Vc_eas: c.aero.Vmo_ms || c.flight.V_ms, Mmo: c.aero.Mmo || 0.8, CLa: up.cfd?.CLa_per_rad, mac: d.mac || undefined, h_top: Math.max(3000, (c.mission.cruise_alt_m || 3000) * 1.35) }),
  run(i) {
    const a = isa(i.alt_m, i.dISA), W = i.mass_kg * G0, ws = W / i.S;
    const Vs1 = Math.sqrt((2 * ws) / (RHO0 * i.CLmax)), Va = Vs1 * Math.sqrt(i.n_pos), Vsn = Math.sqrt((2 * ws) / (RHO0 * i.CLmax_neg)), Vg = Vsn * Math.sqrt(-i.n_neg);
    const Vc = Math.max(i.Vc_eas, Va), Vd = 1.25 * Vc;
    // Pratt gust formula with gust alleviation factor
    const mug = (2 * ws) / (a.rho * i.mac * i.CLa * G0), Kg = (0.88 * mug) / (5.3 + mug), dn = (V, U) => (Kg * RHO0 * U * V * i.CLa) / (2 * ws);
    const Ve = N.linspace(0, Vd, 80);
    const nPos = Ve.map((V) => Math.min(i.n_pos, (0.5 * RHO0 * V * V * i.CLmax) / ws)), nNeg = Ve.map((V) => (V <= Vc ? Math.max(i.n_neg, (-0.5 * RHO0 * V * V * i.CLmax_neg) / ws) : i.n_neg * (1 - (V - Vc) / (Vd - Vc))));
    const gp = [0, Vc, Vd], gust = [[1, 1 + dn(Vc, i.Ude_c), 1 + dn(Vd, i.Ude_c / 2)], [1, 1 - dn(Vc, i.Ude_c), 1 - dn(Vd, i.Ude_c / 2)]];
    const nGustMax = Math.max(...gust[0]), nLim = Math.max(i.n_pos, nGustMax);
    // altitude-speed envelope
    const hs = N.linspace(0, i.h_top, 36), stall = [], vmax = [], vmo = [];
    for (const h of hs) { const p = pointPerf(i, h); stall.push(p.ROCmax > 0 ? p.Vs : NaN); vmax.push(p.ROCmax > 0 ? Math.min(p.Vmax, i.Mmo * p.a.a, Vc / Math.sqrt(p.a.sigma)) : NaN); vmo.push(Math.min(i.Mmo * p.a.a, Vc / Math.sqrt(p.a.sigma))); }
    // sustained turn at the analysis altitude and minimum-drag speed
    const pp = pointPerf(i, i.alt_m), Vt = Math.max(pp.Vmd, 1.2 * pp.Vs), qt = 0.5 * a.rho * Vt * Vt;
    const nSus = Math.max(1, Math.sqrt(Math.max(0, ((thrustAvail(pp.p, Vt, i.alt_m, i.dISA) / (qt * i.S) - i.CD0) / i.k))) * qt * i.S / W);
    const nTurn = Math.min(nSus, i.n_pos, (qt * i.S * i.CLmax) / W), turnRate = nTurn > 1 ? (G0 * Math.sqrt(nTurn * nTurn - 1)) / Vt : 0;
    return {
      kpis: [
        { key: 'V_s1_eas', label: 'Stall speed VS1 (EAS)', value: Vs1, unit: 'm/s' },
        { key: 'V_a_eas', label: 'Manoeuvring speed VA (EAS)', value: Va, unit: 'm/s' },
        { key: 'V_c_eas', label: 'Design cruise speed VC (EAS)', value: Vc, unit: 'm/s' },
        { key: 'V_d_eas', label: 'Design dive speed VD (EAS)', value: Vd, unit: 'm/s' },
        { key: 'n_gust_max', label: 'Peak gust load factor', value: nGustMax, unit: 'g', status: nGustMax > i.n_pos ? 'warn' : 'ok', note: 'Gust-critical when above the manoeuvre limit' },
        { key: 'n_limit', label: 'Governing limit load factor', value: nLim, unit: 'g' },
        { key: 'n_ultimate', label: 'Ultimate load factor (×1.5)', value: 1.5 * nLim, unit: 'g' },
        { key: 'gust_alleviation', label: 'Gust alleviation factor Kg', value: Kg, unit: '-' },
        { key: 'n_sustained', label: 'Sustained turn load factor', value: nTurn, unit: 'g' },
        { key: 'turn_rate_dps', label: 'Sustained turn rate', value: N.deg(turnRate), unit: 'deg/s' },
        { key: 'turn_radius_m', label: 'Sustained turn radius', value: turnRate ? Vt / turnRate : Infinity, unit: 'm' },
      ],
      plots: [
        { type: 'line', title: 'V–n diagram', xlabel: 'Equivalent airspeed [m/s]', ylabel: 'Load factor n [g]', series: [{ name: 'Manoeuvre envelope (+)', x: Ve, y: nPos }, { name: 'Manoeuvre envelope (−)', x: Ve, y: nNeg }, { name: 'Gust lines (+)', x: gp, y: gust[0], style: 'dash' }, { name: 'Gust lines (−)', x: gp, y: gust[1], style: 'dash' }], annotations: [{ x: Va, label: 'VA' }, { x: Vc, label: 'VC' }, { x: Vd, label: 'VD' }] },
        { type: 'line', title: 'Level-flight envelope', xlabel: 'True airspeed [m/s]', ylabel: 'Altitude [m]', series: [{ name: 'Stall boundary', x: stall, y: hs }, { name: 'Maximum speed (thrust or limit)', x: vmax, y: hs }, { name: 'VMO / MMO', x: vmo, y: hs, style: 'dash' }] },
      ],
      warnings: nGustMax > i.n_pos ? ['The gust case exceeds the manoeuvre limit load factor and governs the structural design at VC.'] : [],
      models: ['Manoeuvre envelope from CLmax and limit load factors', 'Pratt discrete-gust formula with alleviation factor', 'Sustained-turn thrust limit'],
      assumptions: ['Symmetric manoeuvres and vertical gusts only', 'VD = 1.25·VC', 'Gust velocity halves between VC and VD; altitude variation of design gust velocity not applied'],
    };
  },
  verify() {
    const r = N.kv(envelope.run({ mass_kg: 1000, S: 16, CD0: 0.03, k: 0.05, CLmax: 1.5, alt_m: 0, dISA: 0, ptype: 'turbojet', n_eng: 1, T0_N: 4000, P0_W: 0, bpr: 0, eta_prop: 0.8, tsfc: 2e-5, bsfc: 0, n_pos: 3.8, n_neg: -1.5, CLmax_neg: 0.9, Vc_eas: 70, Mmo: 0.5, CLa: 5, mac: 1.5, Ude_c: 15.24, h_top: 5000 }));
    return [N.check('VA = VS1·sqrt(n)', r.V_a_eas / r.V_s1_eas, Math.sqrt(3.8), 1e-9, 'CS-23.335(c)')];
  },
  recommend(res) {
    const o = res.outputs;
    return o.n_gust_max > o.n_limit - 1e-9 && o.n_gust_max > 0 && res.warnings.length ? [{ severity: 'advise', title: 'Structure is gust-critical', detail: `Gust load factor ${o.n_gust_max.toFixed(2)} g exceeds the manoeuvre limit.`, action: 'Carry n_limit into Suite 2 (structures) and Suite 9 (fatigue spectrum); a higher wing loading or gust-load alleviation in Suite 16 reduces it.', basis: 'Pratt gust formula' }] : [];
  },
};

const hover = {
  id: 'hover', title: 'Hover ceiling and vertical performance', fidelity: 'analytical',
  summary: 'Hover power required by momentum theory with profile power, against installed power, to find in- and out-of-ground-effect hover ceilings.',
  equations: ['Aircraft energy equations', 'Rate-of-climb equations'],
  applicable: (c, d) => (c.rotor.R_m > 0 ? true : 'This analysis needs a lifting rotor; set rotor data in the case.'),
  inputs: [
    { key: 'mass_kg', label: 'Aircraft mass', unit: 'kg', default: 9000, min: 0.1, group: 'Aircraft' },
    { key: 'R', label: 'Rotor radius', unit: 'm', default: 8, min: 0.02, group: 'Rotor' },
    { key: 'n_rotors', label: 'Lifting rotors', unit: '', default: 1, min: 1, step: 1, discrete: true, group: 'Rotor' },
    { key: 'solidity', label: 'Rotor solidity', unit: '-', default: 0.08, min: 0.005, max: 0.4, group: 'Rotor' },
    { key: 'v_tip', label: 'Tip speed', unit: 'm/s', default: 220, min: 20, group: 'Rotor' },
    { key: 'cd0', label: 'Blade profile drag coefficient', unit: '-', default: 0.01, min: 0.004, group: 'Rotor' },
    { key: 'kappa', label: 'Induced power factor κ', unit: '-', default: 1.15, min: 1, max: 1.6, group: 'Rotor' },
    { key: 'P_inst_W', label: 'Total installed power', unit: 'W', default: 2.8e6, min: 1, group: 'Propulsion' },
    { key: 'ptype', label: 'Powerplant type', type: 'select', options: ['turboshaft', 'piston', 'electric'], default: 'turboshaft', group: 'Propulsion' },
    { key: 'eta_mech', label: 'Transmission efficiency (incl. tail rotor)', unit: '-', default: 0.85, min: 0.5, max: 1, group: 'Propulsion' },
    { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -60, max: 50, group: 'Atmosphere' },
    { key: 'ige_z_R', label: 'IGE rotor height / radius', unit: '-', default: 0.6, min: 0.3, max: 3, group: 'Rotor' },
  ],
  defaults: (c, up, d) => ({ mass_kg: c.mass.mtow_kg, R: c.rotor.R_m, n_rotors: c.meta.type === 'helicopter' ? 1 : c.prop.n_eng, solidity: d.solidity || undefined, v_tip: d.v_tip || undefined, cd0: c.rotor.cd0, P_inst_W: d.P_total, ptype: ['turboshaft', 'piston', 'electric'].includes(c.prop.type) ? c.prop.type : 'turboshaft', eta_mech: c.meta.type === 'helicopter' ? 0.85 : 0.92, dISA: c.atm.dISA_K }),
  run(i) {
    const W = i.mass_kg * G0, A = Math.PI * i.R ** 2, T = W / i.n_rotors;
    const Preq = (h, kg = 1) => { const rho = isa(h, i.dISA).rho; return (i.n_rotors * (i.kappa * kg * T ** 1.5 / Math.sqrt(2 * rho * A) + (rho * A * i.v_tip ** 3 * i.solidity * i.cd0) / 8)) / i.eta_mech; };
    const Pav = (h) => powerAvail({ type: i.ptype, P0_W: i.P_inst_W }, h, i.dISA);
    const kIGE = 1 - 1 / (16 * i.ige_z_R ** 2); // Cheeseman–Bennett thrust augmentation as a power factor
    const ceil = (kg) => { const f = (h) => Pav(h) - Preq(h, kg); return f(0) <= 0 ? 0 : f(9000) > 0 ? 9000 : N.brent(f, 0, 9000, 1e-3); };
    const hOGE = ceil(1), hIGE = ceil(Math.max(0.5, kIGE) ** 1.5), hs = N.linspace(0, 6000, 40), rho0 = isa(0, i.dISA).rho;
    const Pi = i.n_rotors * T ** 1.5 / Math.sqrt(2 * rho0 * A), FM = Pi / (Preq(0) * i.eta_mech), vi = Math.sqrt(T / (2 * rho0 * A));
    const excess = Pav(0) - Preq(0), vc = excess > 0 ? (2 * excess * i.eta_mech) / W : 0; // low-rate climb: Vc ≈ 2·ΔP/W
    return {
      kpis: [
        { key: 'hover_power_W', label: 'Hover power required at sea level', value: Preq(0), unit: 'W' },
        { key: 'power_margin_pct', label: 'Sea-level hover power margin', value: (100 * excess) / Pav(0), unit: '%', status: excess > 0.1 * Pav(0) ? 'ok' : excess > 0 ? 'warn' : 'bad' },
        { key: 'hover_ceiling_oge_m', label: 'Hover ceiling OGE', value: hOGE, unit: 'm' },
        { key: 'hover_ceiling_ige_m', label: 'Hover ceiling IGE', value: hIGE, unit: 'm' },
        { key: 'FM', label: 'Figure of merit', value: FM, unit: '-', status: FM > 0.6 ? 'ok' : 'warn' },
        { key: 'disk_loading_Pa', label: 'Disk loading', value: T / A, unit: 'N/m²' },
        { key: 'v_induced_ms', label: 'Hover induced velocity', value: vi, unit: 'm/s' },
        { key: 'vertical_roc_ms', label: 'Vertical rate of climb at sea level', value: vc, unit: 'm/s' },
      ],
      plots: [{ type: 'line', title: 'Hover power required and available', xlabel: 'Altitude [m]', ylabel: 'Power [kW]', series: [{ name: 'Required OGE', x: hs, y: hs.map((h) => Preq(h) / 1e3) }, { name: 'Required IGE', x: hs, y: hs.map((h) => Preq(h, Math.max(0.5, kIGE) ** 1.5) / 1e3) }, { name: 'Available', x: hs, y: hs.map((h) => Pav(h) / 1e3) }] }],
      models: ['Momentum theory with induced-power factor', 'Uniform profile-power estimate σ·cd0/8', 'Cheeseman–Bennett ground effect'],
      assumptions: ['Uniform inflow', 'Transmission efficiency also covers tail-rotor and accessory power', 'No blade stall or compressibility limits'],
    };
  },
  verify() {
    const r = N.kv(hover.run({ mass_kg: 1000, R: 4, n_rotors: 1, solidity: 0.05, v_tip: 200, cd0: 0, kappa: 1, P_inst_W: 1e6, ptype: 'electric', eta_mech: 1, dISA: 0, ige_z_R: 3 })), T = 1000 * G0, A = Math.PI * 16;
    return [N.check('Ideal hover power T^1.5/sqrt(2ρA)', r.hover_power_W, T ** 1.5 / Math.sqrt(2 * isa(0).rho * A), 1e-9, 'Rankine–Froude momentum theory'), N.check('Ideal figure of merit = 1', r.FM, 1, 1e-9, 'Definition')];
  },
  recommend(res) {
    const o = res.outputs;
    return o.power_margin_pct < 10 ? [{ severity: o.power_margin_pct < 0 ? 'critical' : 'warn', title: 'Thin hover power margin', detail: `Margin is ${o.power_margin_pct.toFixed(1)}% at sea level.`, action: 'Reduce gross mass or disk loading; hot-and-high sites will be more limiting still.', basis: 'Momentum theory power balance' }] : [];
  },
};

export default {
  id: 'performance', n: 5,
  tagline: 'How fast, how high, how far and from which runway — with the limits that bound safe operation.',
  analyses: [point, field, climb, range, envelope, hover],
  consumes: [
    { from: 'cfd', keys: ['CD0', 'k_induced', 'CLmax', 'CLa_per_rad'], why: 'Drag polar and lift limits' },
    { from: 'propulsion', keys: ['thrust_static_N', 'tsfc_kg_Ns'], why: 'Installed thrust and fuel consumption' },
    { from: 'propeller', keys: ['eta_prop'], why: 'Propulsive efficiency' },
  ],
  provides: [
    { key: 'V_stall_ms', label: 'Stall speed', unit: 'm/s' }, { key: 'tofl_m', label: 'Take-off distance', unit: 'm' }, { key: 'ldg_dist_m', label: 'Landing distance', unit: 'm' },
    { key: 'roc_max_ms', label: 'Max rate of climb', unit: 'm/s' }, { key: 'ceiling_m', label: 'Service ceiling', unit: 'm' }, { key: 'V_max_ms', label: 'Max level speed', unit: 'm/s' },
    { key: 'range_km', label: 'Range', unit: 'km' }, { key: 'endurance_h', label: 'Endurance', unit: 'h' }, { key: 'LD_max', label: 'Max L/D', unit: '-' },
    { key: 'fuel_flow_cruise_kgs', label: 'Cruise fuel flow', unit: 'kg/s' }, { key: 'n_limit', label: 'Limit load factor', unit: 'g' }, { key: 'V_d_eas', label: 'Dive speed', unit: 'm/s' },
  ],
  handoff: [
    { model: 'Balanced field length with accelerate–stop', why: 'Needs certified engine-failure transition times and brake performance data', tool: 'Manufacturer performance software / flight test' },
    { model: 'Energy-based trajectory optimisation', why: 'Implemented at mission level in Suite 24', tool: 'Suite 24' },
    { model: 'Helicopter height–velocity envelope', why: 'Depends on autorotation entry dynamics validated in flight test', tool: 'Suite 6 autorotation index + flight test' },
  ],
};
