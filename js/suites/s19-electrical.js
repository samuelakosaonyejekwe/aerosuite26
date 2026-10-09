// Suite 19 — Electrical and Hybrid-Electric Power Systems.
// Battery equivalent-circuit discharge with thermal coupling, ageing and thermal-runaway propagation, PMSM dq-axis
// envelope and efficiency map with inverter losses, DC-equivalent nodal load flow with fault levels, load shedding and
// electrical load analysis, cable sizing versus bus voltage with the Paschen limit, hybrid-electric energy management
// (rule-based and dynamic programming) and averaged DC/DC converter design with the bus load-step transient.

import * as N from '../core/numerics.js';
import { isa, G0 } from '../core/atmosphere.js';
import { FLUIDS, BATTERIES } from '../data/materials.js';

// ---- shared helpers -------------------------------------------------------------------------
const CHEMS = Object.keys(BATTERIES), FUELS = ['Jet A-1', 'Avgas 100LL', 'SAF (HEFA-SPK)', 'Liquid hydrogen'];
const chemOf = (name) => BATTERIES[name] || BATTERIES[CHEMS[0]];
const isElec = (c) => c.prop.type === 'electric';
const RHO_CU = 1.724e-8, ALPHA_CU = 0.00393; // copper resistivity at 20 °C [Ω m] (International Annealed Copper Standard, 58 MS/m) and its usual temperature coefficient [1/K]
/** Down-sample a history to at most m points, keeping the last one. */
const ds = (a, m = 300) => { if (a.length <= m) return a; const k = Math.ceil(a.length / m), o = []; for (let j = 0; j < a.length; j += k) o.push(a[j]); if ((a.length - 1) % k) o.push(a[a.length - 1]); return o; };

/** Class-level power estimates from the shared case (used only to fill defaults; upstream suites override them). */
function est(c, up, d) {
  const a = isa(c.atm.alt_m, c.atm.dISA_K), W = c.mass.mtow_kg * G0, r = c.rotor, n = Math.max(1, c.prop.n_eng), S = c.wing.S_m2;
  const V = c.mission.cruise_V_ms || c.flight.V_ms || 30, nr = c.meta.type === 'helicopter' ? 1 : n, etaP = N.clamp(up.propeller?.eta_prop ?? c.prop.eta_prop ?? 0.8, 0.3, 0.95);
  let Ph = 0, Pc = 0;
  if (r.R_m > 0 && d.v_tip > 0) { // momentum theory hover (κ = 1.15) + σ·cd0/8 profile power; Glauert high-speed induced power in cruise
    const T = W / nr, vh = Math.sqrt(T / (2 * a.rho * d.A_disk)), Pi = 1.15 * T * vh, P0 = (a.rho * d.A_disk * d.v_tip ** 3 * d.solidity * r.cd0) / 8, mu = V / d.v_tip;
    Ph = nr * (Pi + P0);
    if (!(S > 0)) Pc = nr * ((Pi * vh) / Math.hypot(V, vh) + P0 * (1 + 4.65 * mu * mu)) + 0.5 * a.rho * V ** 3 * (r.flat_plate_m2 || 0);
  }
  if (S > 0) { const q = 0.5 * a.rho * V * V, CL = W / (q * S); Pc = (q * S * (c.aero.CD0 + (d.k_induced || 0.05) * CL * CL) * V) / etaP; }
  const rc = up.rotorcraft || {};
  if (r.R_m > 0 && rc.hover_power_W > 0) Ph = rc.hover_power_W;
  if (!(S > 0) && rc.P_cruise_W > 0) Pc = rc.P_cruise_W;
  const genTot = Math.max(50, (c.systems.gen_kVA || 0) * 900 * Math.min(2, n)); // W at 0.9 power factor, at most two channels
  const Ppeak = isElec(c) || d.isRotary ? Math.max(1.05 * Ph, Pc) : d.P_total > 0 ? d.P_total : 2.2 * Pc;
  const tCruise = c.mission.range_km > 0 ? (c.mission.range_km * 1e3) / V / 60 : 60;
  return { electric: isElec(c), Ph, Pc, V, n, genTot, Pess: 0.5 * genTot ** 0.75, Ppeak, tCruise, len: Math.max(c.fuselage.len_m || 0, 2 * r.R_m, 0.5) };
}

// ---- battery equivalent-circuit model -------------------------------------------------------
// Generic open-circuit-voltage shape g(SOC) in [0,1] (layered-oxide-like); scaled between the chemistry voltage limits.
const OCV_S = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], OCV_G = [0, 0.25, 0.375, 0.458, 0.508, 0.542, 0.583, 0.65, 0.725, 0.808, 0.9, 1];
/** Pack parameters from energy, nominal voltage and chemistry. Resistance split 60/25/15 % over R0, R1 (15 s), R2 (300 s). */
function makePack(i) {
  const ch = chemOf(i.chem), Ns = Math.max(1, Math.round(i.V_nom / ch.v_nom)), vlo = ch.v_min + 0.2;
  const ocv1 = (s) => vlo + (ch.v_max - vlo) * N.interp1(OCV_S, OCV_G, N.clamp(s, 0, 1));
  const vmean = N.trapz(OCV_S, OCV_S.map(ocv1)), Q_Ah = (i.E_kWh * 1e3) / (Ns * vmean), R = (i.r_scale * ch.r_mohm_ah * 1e-3 * Ns) / Q_Ah, mass = (i.E_kWh * 1e3) / ch.wh_kg / i.pack_frac;
  return { ch, Ns, Q_Ah, Q: Q_Ah * 3600, R, R0: 0.6 * R, R1: 0.25 * R, tau1: 15, R2: 0.15 * R, tau2: 300, ocv: (s) => Ns * ocv1(s), Vcut: Ns * ch.v_min, mass, Cth: mass * ch.cp, hA: i.hA_W_K, Tamb: i.T_amb_K, kT: i.R_temp_K };
}
/**
 * Integrate the 2-RC equivalent circuit with a lumped thermal node through a list of segments {P | I, t, n} by RK4.
 * State: [SOC, v1, v2, T, E_terminal, E_heat, E_chemical]. Stops at the cut-off voltage, at SOC = socStop or on power collapse.
 */
function ecmSim(pk, segs, y0, socStop = 0) {
  const f = (y, sg) => {
    const rT = pk.kT ? Math.exp(pk.kT * (1 / y[3] - 1 / 298.15)) : 1, R0 = pk.R0 * rT, R1 = pk.R1 * rT, R2 = pk.R2 * rT, e = pk.ocv(y[0]), voc = e - y[1] - y[2];
    let I, bad = false;
    if (sg.I !== undefined) I = sg.I; else { const disc = voc * voc - 4 * R0 * sg.P; bad = disc < 0; I = (voc - Math.sqrt(Math.max(0, disc))) / (2 * R0); }
    const V = voc - I * R0, heat = I * I * R0 + (y[1] * y[1]) / R1 + (y[2] * y[2]) / R2;
    return { dy: [-I / pk.Q, (I * R1 - y[1]) / pk.tau1, (I * R2 - y[2]) / pk.tau2, (heat - pk.hA * (y[3] - pk.Tamb)) / pk.Cth, V * I, heat, e * I], V, I, heat, bad };
  };
  const H = { t: [], V: [], I: [], soc: [], T: [], heat: [] }, rec = (t, y, o) => { H.t.push(t); H.V.push(o.V); H.I.push(o.I); H.soc.push(y[0]); H.T.push(y[3]); H.heat.push(o.heat); };
  let y = y0.slice(), t = 0, stop = '';
  outer: for (const sg of segs) {
    if (!(sg.t > 0)) continue;
    const h = sg.t / sg.n; let o = f(y, sg); rec(t, y, o);
    for (let s = 0; s < sg.n; s++) {
      const k1 = o.dy, k2 = f(N.vadd(y, k1, h / 2), sg).dy, k3 = f(N.vadd(y, k2, h / 2), sg).dy, k4 = f(N.vadd(y, k3, h), sg).dy;
      const yn = y.map((v, j) => v + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])), on = f(yn, sg);
      if (yn[0] <= socStop) { const fr = N.clamp((y[0] - socStop) / (y[0] - yn[0] || 1), 0, 1); y = y.map((v, j) => v + fr * (yn[j] - v)); t += fr * h; rec(t, y, f(y, sg)); stop = socStop > 0 ? 'reserve' : 'empty'; break outer; }
      y = yn; t += h; o = on; rec(t, y, o);
      if (o.bad) { stop = 'collapse'; break outer; }
      if (o.V < pk.Vcut) { stop = 'cutoff'; break outer; }
    }
  }
  return { ...H, y, tEnd: t, stop };
}

const battery = {
  id: 'battery', title: 'Battery pack discharge, voltage sag and temperature', fidelity: 'numerical',
  summary: 'Sizes the pack in series and parallel, then discharges a two-RC equivalent circuit through the mission power profile with a lumped thermal model: state of charge, voltage sag, C-rate, heat and endurance.',
  equations: ['Battery equivalent circuit equations', 'Ohm’s law', 'Kirchhoff’s voltage law', 'Power balance equations', 'Thermal energy equations'],
  inputs: [
    { key: 'chem', label: 'Cell chemistry', type: 'select', options: CHEMS, default: CHEMS[0], group: 'Battery', help: 'Sets specific energy, voltage window, resistance and heat capacity from the material database' },
    { key: 'E_kWh', label: 'Pack energy (rated, new)', unit: 'kWh', default: 160, min: 0.001, group: 'Battery' },
    { key: 'V_nom', label: 'Nominal pack voltage', unit: 'V', default: 800, min: 3, max: 3000, group: 'Battery', help: 'Sets the series cell count' },
    { key: 'cell_Ah', label: 'Cell capacity', unit: 'Ah', default: 5, min: 0.1, max: 500, group: 'Battery', help: 'Only used for the parallel count; 3–5 Ah cylindrical, 50–150 Ah pouch/prismatic' },
    { key: 'pack_frac', label: 'Cell-to-pack mass fraction', unit: '-', default: 0.75, min: 0.4, max: 1, group: 'Battery' },
    { key: 'r_scale', label: 'Resistance scale factor', unit: '-', default: 1, min: 0.05, max: 10, group: 'Battery', help: 'Multiplies the database cell resistance; >1 for aged or cold cells, <1 for power cells. Calibrate against a discharge test' },
    { key: 'R_temp_K', label: 'Resistance temperature coefficient', unit: 'K', default: 1500, min: 0, max: 6000, group: 'Battery', help: 'R(T) = R·exp[k·(1/T − 1/298 K)]; 0 switches the dependence off. Generic value' },
    { key: 'soc0', label: 'Initial state of charge', unit: '-', default: 1, min: 0.05, max: 1, group: 'Initial conditions' },
    { key: 'T0_K', label: 'Initial pack temperature', unit: 'K', default: 298.15, min: 233, max: 333, group: 'Initial conditions' },
    { key: 'P1_W', label: 'Phase 1 shaft power (hover / take-off)', unit: 'W', default: 5e5, min: 0, group: 'Mission', help: 'Applied in two halves: at the start and at the end of the mission' },
    { key: 't1_min', label: 'Phase 1 total time', unit: 'min', default: 2, min: 0, group: 'Mission' },
    { key: 'P2_W', label: 'Phase 2 power (cruise / essential loads)', unit: 'W', default: 1.5e5, min: 0, group: 'Mission' },
    { key: 't2_min', label: 'Phase 2 time', unit: 'min', default: 30, min: 0, group: 'Mission' },
    { key: 'P_aux_W', label: 'Auxiliary electrical load', unit: 'W', default: 3000, min: 0, group: 'Mission', help: 'Avionics, thermal management, payload: drawn at the battery terminals' },
    { key: 'eta_drive', label: 'Battery-to-shaft efficiency', unit: '-', default: 0.9, min: 0.3, max: 1, group: 'Mission', help: 'Cables × inverter × motor; 1 when the powers above are already electrical' },
    { key: 'soc_res', label: 'Reserve state of charge', unit: '-', default: 0.15, min: 0, max: 0.6, group: 'Limits' },
    { key: 'T_amb_K', label: 'Coolant / ambient temperature', unit: 'K', default: 298.15, min: 220, max: 330, group: 'Cooling' },
    { key: 'hA_W_K', label: 'Pack cooling conductance h·A', unit: 'W/K', default: 320, min: 0, group: 'Cooling', help: 'Heat removed per kelvin above coolant; ~2 W/K per kWh for liquid-cooled packs' },
    { key: 'nSteps', label: 'Time steps over the mission', unit: '', default: 400, min: 20, max: 20000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const e = est(c, up, d), E = c.systems.batt_kWh || undefined, ch = chemOf(c.systems.batt_chem), mass = ((E || 1) * 1e3) / ch.wh_kg / 0.75;
    const T = isa(c.atm.alt_m, c.atm.dISA_K).T, area = 6 * (mass / 2000) ** (2 / 3);
    if (!e.electric) return { chem: c.systems.batt_chem, E_kWh: E, V_nom: 28, P1_W: 0, t1_min: 0, P2_W: Math.min(e.Pess, 1400 * (E || 1)), t2_min: 30, P_aux_W: 0, eta_drive: 1, soc_res: 0.1, hA_W_K: 8 * area, T_amb_K: Math.max(273, Math.min(T, 313)) };
    return { chem: c.systems.batt_chem, E_kWh: E, V_nom: c.systems.bus_V || undefined, P1_W: e.Ph, t1_min: e.Ph > 0 ? c.mission.hover_min : 0, P2_W: e.Pc || undefined, t2_min: e.tCruise, P_aux_W: 0.02 * Math.max(e.Pc, e.Ph), hA_W_K: Math.max(25 * area, 2 * (E || 0)), T_amb_K: N.clamp(T, 263, 313) };
  },
  run(i) {
    const pk = makePack(i), n = Math.max(8, Math.round(i.nSteps)), Pe = (P) => P / i.eta_drive + i.P_aux_W, t1 = i.t1_min * 30, t2 = i.t2_min * 60, tt = 2 * t1 + t2 || 1;
    const nOf = (t) => Math.max(4, Math.round((n * t) / tt)), segs = [{ P: Pe(i.P1_W), t: t1, n: nOf(t1) }, { P: Pe(i.P2_W), t: t2, n: nOf(t2) }, { P: Pe(i.P1_W), t: t1, n: nOf(t1) }];
    const y0 = [i.soc0, 0, 0, i.T0_K, 0, 0, 0], r = ecmSim(pk, segs, y0), socEnd = r.y[0], warnings = [];
    // endurance: constant phase-2 power (phase 1 when there is no phase 2) from the initial state down to the reserve
    const Pend = Pe(i.P2_W > 0 ? i.P2_W : i.P1_W), tMax = (2 * i.E_kWh * 3.6e6) / Math.max(Pend, 1e-9);
    const en = Pend > 0 ? ecmSim(pk, [{ P: Pend, t: tMax, n: Math.max(100, n) }], y0, i.soc_res) : { tEnd: Infinity };
    const Imax = N.amax(r.I), Vmin = N.amin(r.V), Tmax = N.amax(r.T), cRate = Imax / pk.Q_Ah, Eout = r.y[4] / 3.6e6;
    const Ecap = (0.5 * pk.tau1 * r.y[1] ** 2) / pk.R1 + (0.5 * pk.tau2 * r.y[2] ** 2) / pk.R2, bal = r.y[6] ? Math.abs(r.y[4] + r.y[5] + Ecap - r.y[6]) / r.y[6] : 0;
    const done = !r.stop;
    if (!done) warnings.push(r.stop === 'collapse' ? 'The pack cannot deliver the demanded power: the terminal voltage collapses (demand exceeds the maximum-power point).' : `The mission is not completed: the pack reaches ${r.stop === 'cutoff' ? 'its cut-off voltage' : 'zero charge'} after ${(r.tEnd / 60).toFixed(1)} of ${(tt / 60).toFixed(1)} min.`);
    if (cRate > 5) warnings.push(`Peak discharge rate is ${cRate.toFixed(1)}C: beyond the continuous rating of typical energy cells; the database resistance and specific energy do not describe power cells.`);
    if (Tmax > 333.15) warnings.push('Pack temperature exceeds 60 °C, a typical operating limit for lithium-ion cells.');
    if (i.T0_K < 273.15) warnings.push('Below 0 °C the generic resistance law is extrapolated and usable capacity loss is not modelled.');
    const tm = ds(r.t).map((v) => v / 60);
    return {
      kpis: [
        { key: 'batt_soc_end', label: 'State of charge at end of mission', value: socEnd, unit: '-', status: !done ? 'bad' : socEnd >= i.soc_res ? 'ok' : 'warn', note: `Reserve ${(100 * i.soc_res).toFixed(0)}%` },
        { key: 'pack_V_min_V', label: 'Minimum pack terminal voltage', value: Vmin, unit: 'V', status: Vmin > pk.Vcut ? 'ok' : 'bad', note: `Cut-off ${pk.Vcut.toFixed(1)} V` },
        { key: 'pack_sag_pct', label: 'Peak voltage sag below open-circuit', value: 100 * N.amax(r.V.map((v, j) => 1 - v / pk.ocv(r.soc[j]))), unit: '%' },
        { key: 'pack_I_peak_A', label: 'Peak pack current', value: Imax, unit: 'A' },
        { key: 'c_rate_peak', label: 'Peak discharge rate', value: cRate, unit: 'C', status: cRate <= 5 ? 'ok' : 'warn', note: 'Energy cells: up to about 3–5C' },
        { key: 'batt_temp_K', label: 'Peak pack temperature', value: Tmax, unit: 'K', status: Tmax <= 333.15 ? 'ok' : 'bad', note: 'Typical limit 333 K (60 °C)' },
        { key: 'energy_used_kWh', label: 'Energy delivered at the terminals', value: Eout, unit: 'kWh' },
        { key: 'batt_heat_kWh', label: 'Heat generated in the pack', value: r.y[5] / 3.6e6, unit: 'kWh' },
        { key: 'batt_heat_peak_W', label: 'Peak pack heat generation', value: N.amax(r.heat), unit: 'W' },
        { key: 'batt_discharge_eff', label: 'Discharge energy efficiency', value: r.y[6] ? r.y[4] / r.y[6] : 1, unit: '-' },
        { key: 'endurance_elec_min', label: 'Endurance to reserve at phase-2 power', value: en.tEnd / 60, unit: 'min' },
        { key: 'n_series', label: 'Cells in series', value: pk.Ns, unit: '' },
        { key: 'n_parallel', label: 'Cells in parallel', value: Math.max(1, Math.ceil(pk.Q_Ah / i.cell_Ah - 1e-9)), unit: '' },
        { key: 'pack_capacity_Ah', label: 'Pack capacity', value: pk.Q_Ah, unit: 'Ah' },
        { key: 'pack_R_ohm', label: 'Pack DC resistance at 25 °C', value: pk.R, unit: 'Ω' },
        { key: 'pack_mass_kg', label: 'Pack mass', value: pk.mass, unit: 'kg' },
        { key: 'energy_balance_err', label: 'Energy accounting residual', value: bal, unit: '-', note: 'Chemical = terminal + heat + RC-stored' },
      ],
      plots: [
        { type: 'line', title: 'Pack voltage', xlabel: 'Time [min]', ylabel: 'Voltage [V]', series: [{ name: 'Terminal voltage', x: tm, y: ds(r.V) }, { name: 'Open-circuit voltage', x: tm, y: ds(r.soc).map(pk.ocv), style: 'dash' }], annotations: [{ y: pk.Vcut, label: 'Cut-off' }] },
        { type: 'line', title: 'State of charge', xlabel: 'Time [min]', ylabel: 'SOC [-]', series: [{ name: 'SOC', x: tm, y: ds(r.soc) }], annotations: [{ y: i.soc_res, label: 'Reserve' }] },
        { type: 'line', title: 'Pack current', xlabel: 'Time [min]', ylabel: 'Current [A]', series: [{ name: 'Discharge current', x: tm, y: ds(r.I) }] },
        { type: 'line', title: 'Pack temperature', xlabel: 'Time [min]', ylabel: 'Temperature [K]', series: [{ name: 'Lumped pack temperature', x: tm, y: ds(r.T) }], annotations: [{ y: 333.15, label: '60 °C limit' }] },
      ],
      tables: [{ title: 'Pack build', columns: ['Quantity', 'Value', 'Unit'], rows: [['Series cells', pk.Ns, '-'], ['Capacity', pk.Q_Ah, 'Ah'], ['R0', pk.R0, 'Ω'], ['R1 (τ = 15 s)', pk.R1, 'Ω'], ['R2 (τ = 300 s)', pk.R2, 'Ω'], ['Thermal mass', pk.Cth, 'J/K'], ['Open-circuit voltage, full', pk.ocv(1), 'V'], ['Cut-off voltage', pk.Vcut, 'V']] }],
      warnings, models: ['Two-RC equivalent-circuit battery model with constant-power terminal load', 'Generic OCV(SOC) curve scaled to the chemistry voltage window', 'Lumped pack thermal model with Joule and polarisation heat', 'RK4 time integration'],
      assumptions: ['Rated energy equals ∫OCV dq, so capacity loss with rate appears only through resistive heat (no Peukert exponent)', 'Resistance split 60/25/15 % between R0 and two RC pairs with 15 s and 300 s time constants (generic)', 'Reversible (entropic) heat and cell-to-cell imbalance neglected', 'Uniform pack temperature', 'Cell data (specific energy, voltage window, resistance, heat capacity) are typical class values from the material database, optimistic against the datasheets checked: a high-power 21700 cell gives about 231 Wh/kg with a 2.5 V cut-off, an 18650 NCA cell 243 Wh/kg', 'Pack mass fraction, cooling conductance and the resistance temperature coefficient are illustrative defaults'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps over the mission', levels: [50, 100, 200, 400], metric: 'batt_temp_K' },
  calibration: { params: [{ key: 'r_scale', min: 0.05, max: 10 }, { key: 'hA_W_K', min: 0, max: 1e5 }], sweep: 'P2_W', target: 'pack_V_min_V', note: 'Supply measured minimum terminal voltage (or peak temperature) from constant-power discharge tests at several power levels.' },
  verify() {
    // Constant-current pulse on a flat-OCV circuit: analytic two-RC response and coulomb counting.
    const pk = { Q: 36000, R0: 0.01, R1: 0.005, tau1: 15, R2: 0.004, tau2: 300, ocv: () => 100, Vcut: 0, Cth: 1e6, hA: 0, Tamb: 300, kT: 0 }, I = 50, t = 60;
    const r = ecmSim(pk, [{ I, t, n: 240 }], [1, 0, 0, 300, 0, 0, 0]), Vex = 100 - I * (0.01 + 0.005 * (1 - Math.exp(-t / 15)) + 0.004 * (1 - Math.exp(-t / 300)));
    const o = N.kv(battery.run({ ...Object.fromEntries(battery.inputs.map((f) => [f.key, f.default])), R_temp_K: 0 }));
    return [
      N.check('2-RC terminal voltage after a 60 s current pulse', r.V[r.V.length - 1], Vex, 1e-8, 'V = OCV − I·[R0 + ΣRk(1 − e^(−t/τk))]'),
      N.check('Coulomb counting', r.y[0], 1 - (I * t) / 36000, 1e-10, 'dSOC/dt = −I/Q'),
      N.check('Energy conservation over the default mission', o.energy_balance_err, 0, 1e-6, 'Chemical energy = terminal energy + heat + energy stored in the RC pairs'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.batt_soc_end < i.soc_res) out.push({ severity: o.batt_soc_end <= 0.02 ? 'critical' : 'warn', title: 'Mission ends below the reserve state of charge', detail: `End SOC ${(100 * o.batt_soc_end).toFixed(1)}% against a ${(100 * i.soc_res).toFixed(0)}% reserve.`, action: 'Shorten the hover or cruise segment, add pack energy, or reduce power demand (mass, disk loading, drag).', basis: 'Energy reserve policy set in the inputs' });
    if (o.batt_temp_K > 318) out.push({ severity: o.batt_temp_K > 333.15 ? 'critical' : 'advise', title: 'Pack runs hot', detail: `Peak ${(o.batt_temp_K - 273.15).toFixed(0)} °C. Sustained operation above about 45 °C accelerates ageing; above 60 °C is outside typical cell limits.`, action: 'Increase cooling conductance, pre-cool before flight, or lower the peak C-rate with more parallel capacity. Pass the heat load to Suite 12.', basis: 'Cell operating window; Arrhenius ageing' });
    if (o.c_rate_peak > 3) out.push({ severity: 'advise', title: 'High discharge rate', detail: `Peak ${o.c_rate_peak.toFixed(1)}C with ${o.pack_sag_pct.toFixed(1)}% voltage sag.`, action: 'Select power-optimised cells (lower resistance, lower specific energy) or raise pack energy; check end-of-life resistance growth in the ageing analysis.', basis: 'I²R loss grows with the square of C-rate' });
    out.push({ severity: 'info', title: 'Energy use and battery life', detail: `${o.energy_used_kWh.toPrecision(3)} kWh per mission at ${(100 * o.batt_discharge_eff).toFixed(1)}% discharge efficiency.`, action: 'Charging from low-carbon electricity sets the operational CO₂; shallow cycling and moderate temperatures extend pack life and cut embodied material use.', basis: 'Energy accounting of this run' });
    return out;
  },
};

// ---- ageing and thermal-runaway propagation -------------------------------------------------
/** Chain of lumped cells with conduction to neighbours, convection and a one-step Arrhenius heat release (rate capped at 5/τ). */
function trSim(p) {
  const n = p.n, f = (t, y) => {
    const dy = new Array(2 * n);
    for (let j = 0; j < n; j++) {
      const T = y[j], x = y[n + j], rx = p.E > 0 ? (Math.min(5, Math.exp(p.EaR * (1 / p.Tonset - 1 / T))) / p.tau) * (1 - x) : 0;
      let q = p.E * rx - p.hA * (T - p.Tamb);
      if (j > 0) q += p.G * (y[j - 1] - T);
      if (j < n - 1) q += p.G * (y[j + 1] - T);
      dy[j] = q / p.C; dy[n + j] = rx;
    }
    return dy;
  };
  return N.rk4(f, 0, [p.T0, ...new Array(n - 1).fill(p.Tamb), ...new Array(n).fill(0)], p.tEnd, p.nSteps);
}
/** Semi-empirical capacity loss after t years: linear cycle fade scaled by depth of discharge, rate and temperature, plus √t calendar fade. */
function fade(i, t) {
  const ch = chemOf(i.chem), arr = (T) => Math.exp((i.Ea_J_mol / 8.314) * (1 / 298.15 - 1 / T));
  const N80 = (ch.cycles_80 * (0.8 / Math.max(i.dod, 0.02)) ** i.k_dod) / (Math.max(1, arr(i.T_cell_K)) * Math.exp(i.k_rate * Math.max(0, i.c_rate - 1)));
  return { cyc: (0.2 * i.cycles_yr * t) / N80, cal: i.k_cal * arr(i.T_store_K) * (0.5 + i.soc_store) * Math.sqrt(t), N80 };
}

const ageing = {
  id: 'ageing', title: 'Battery ageing and thermal-runaway propagation', fidelity: 'reduced-order',
  summary: 'Projects capacity fade from cycling and storage with a labelled semi-empirical law, and simulates whether a single-cell thermal runaway spreads along a row of cells through the inter-cell thermal path.',
  equations: ['Thermal energy equations', 'Battery ageing models (semi-empirical)', 'Battery thermal runaway models (lumped, one-step Arrhenius)'],
  inputs: [
    { key: 'chem', label: 'Cell chemistry', type: 'select', options: CHEMS, default: CHEMS[0], group: 'Battery' },
    { key: 'cycles_yr', label: 'Charge–discharge cycles per year', unit: '1/yr', default: 1500, min: 0, group: 'Duty' },
    { key: 'dod', label: 'Depth of discharge per cycle', unit: '-', default: 0.8, min: 0.02, max: 1, group: 'Duty' },
    { key: 'c_rate', label: 'Mean discharge rate', unit: 'C', default: 1, min: 0.05, max: 20, group: 'Duty' },
    { key: 'T_cell_K', label: 'Cell temperature while cycling', unit: 'K', default: 303, min: 263, max: 343, group: 'Duty' },
    { key: 'T_store_K', label: 'Storage temperature', unit: 'K', default: 298, min: 253, max: 333, group: 'Duty' },
    { key: 'soc_store', label: 'Storage state of charge', unit: '-', default: 0.5, min: 0, max: 1, group: 'Duty' },
    { key: 'years', label: 'Projection horizon', unit: 'yr', default: 10, min: 0.5, max: 40, group: 'Duty' },
    { key: 'k_dod', label: 'Depth-of-discharge exponent', unit: '-', default: 1, min: 0, max: 3, group: 'Ageing model', help: 'Cycle life ∝ (0.8/DoD)^k. Generic; calibrate' },
    { key: 'k_rate', label: 'Rate stress coefficient', unit: '1/C', default: 0.15, min: 0, max: 2, group: 'Ageing model', help: 'Cycle life ÷ exp[k·(C − 1)] above 1C. Generic; calibrate' },
    { key: 'k_cal', label: 'Calendar fade in the first year', unit: '-', default: 0.02, min: 0, max: 0.2, group: 'Ageing model', help: 'Capacity fraction lost in year one at 25 °C and 50 % SOC; grows with √time' },
    { key: 'Ea_J_mol', label: 'Ageing activation energy', unit: 'J/mol', default: 4e4, min: 0, max: 1.2e5, group: 'Ageing model' },
    { key: 'n_cells', label: 'Cells in the row', unit: '', default: 8, min: 1, max: 40, step: 1, discrete: true, group: 'Thermal runaway' },
    { key: 'cell_Ah', label: 'Cell capacity', unit: 'Ah', default: 5, min: 0.1, max: 500, group: 'Thermal runaway' },
    { key: 'soc', label: 'State of charge at the event', unit: '-', default: 1, min: 0, max: 1, group: 'Thermal runaway' },
    { key: 'heat_ratio', label: 'Heat released / stored electrical energy', unit: '-', default: 1.5, min: 0, max: 5, group: 'Thermal runaway', help: 'Total heat release: fractional thermal-runaway calorimetry of 18650 cells at full charge gave 1.3–1.6 (NASA, Walker et al. 2018). Depends on chemistry and state of charge; all of it is kept in the cell body here' },
    { key: 'T_onset_K', label: 'Runaway onset temperature', unit: 'K', default: 423, min: 350, max: 600, group: 'Thermal runaway', help: 'Generic 150 °C; lower for high-nickel cells, higher for LiFePO4' },
    { key: 'tau_s', label: 'Reaction time constant at onset', unit: 's', default: 5, min: 0.5, max: 600, group: 'Thermal runaway' },
    { key: 'G_W_K', label: 'Cell-to-cell thermal conductance', unit: 'W/K', default: 0.3, min: 0, max: 50, group: 'Thermal runaway', help: 'Direct contact ≈ 0.2–1 W/K for small cells; insulating barriers ≈ 0.01–0.05 W/K' },
    { key: 'hA_W_K', label: 'Cell cooling conductance', unit: 'W/K', default: 0.05, min: 0, max: 50, group: 'Thermal runaway' },
    { key: 'T_amb_K', label: 'Initial / coolant temperature', unit: 'K', default: 300, min: 233, max: 340, group: 'Thermal runaway' },
    { key: 't_end_s', label: 'Simulated time', unit: 's', default: 600, min: 10, max: 36000, group: 'Numerics' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 2400, min: 100, max: 100000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const e = est(c, up, d), E = (c.systems.batt_kWh || 0) * 1e3;
    return e.electric ? { chem: c.systems.batt_chem, cycles_yr: c.econ.cycles_yr || undefined, c_rate: E > 0 ? N.clamp(Math.max(e.Pc, 0.3 * e.Ph) / 0.9 / E, 0.1, 10) : undefined, T_cell_K: 303 } : { chem: c.systems.batt_chem, cycles_yr: 50, dod: 0.2, c_rate: 0.5, soc_store: 0.9, T_cell_K: 298 };
  },
  run(i) {
    const ch = chemOf(i.chem), ts = N.linspace(0, i.years, 81), fd = ts.map((t) => fade(i, t)), soh = fd.map((f) => Math.max(0, 1 - f.cyc - f.cal)), warnings = [];
    const g = (t) => { const f = fade(i, t); return 0.8 - (1 - f.cyc - f.cal); }, life = g(60) < 0 ? 60 : N.brent(g, 0, 60, 1e-8), fEnd = fd[80];
    if (life >= 60) warnings.push('The 80 % end-of-life criterion is not reached within 60 years: the result is capped.');
    if (i.T_cell_K < 283) warnings.push('Low-temperature degradation (lithium plating) is not represented: life is optimistic below about 10 °C.');
    // thermal-runaway propagation along a row of cells
    const mCell = (i.cell_Ah * ch.v_nom) / ch.wh_kg, n = Math.round(i.n_cells), E = i.heat_ratio * i.soc * i.cell_Ah * ch.v_nom * 3600;
    const r = trSim({ n, C: mCell * ch.cp, G: i.G_W_K, hA: i.hA_W_K, Tamb: i.T_amb_K, Tonset: i.T_onset_K, E, tau: i.tau_s, EaR: 15600, T0: i.T_onset_K + 30, tEnd: i.t_end_s, nSteps: Math.round(i.nSteps) });
    const trig = N.range(n, (j) => { const k = r.y.findIndex((y) => y[n + j] > 0.5); return k < 0 ? NaN : r.t[k]; }), nTrig = trig.filter(Number.isFinite).length;
    const tLast = N.amax(trig.filter(Number.isFinite).concat(0)), Tpk = N.amax(r.y.map((y) => N.amax(y.slice(0, n)))), idx = ds(N.range(r.t.length), 240);
    const show = n <= 6 ? N.range(n) : [0, 1, 2, Math.floor(n / 2), n - 2, n - 1];
    if (nTrig > 1 && nTrig < n && tLast > 0.8 * i.t_end_s) warnings.push('Propagation is still in progress at the end of the simulated time: extend it.');
    return {
      kpis: [
        { key: 'batt_life_yr', label: 'Life to 80 % capacity', value: life, unit: 'yr' },
        { key: 'cycles_to_eol', label: 'Cycles to 80 % capacity', value: life * i.cycles_yr, unit: '' },
        { key: 'soh_end', label: `Capacity retention after ${i.years} yr`, value: soh[80], unit: '-', status: soh[80] >= 0.8 ? 'ok' : 'warn', note: 'End of life at 0.80' },
        { key: 'fade_cycle', label: 'Cycle fade at the horizon', value: Math.min(1, fEnd.cyc), unit: '-' },
        { key: 'fade_calendar', label: 'Calendar fade at the horizon', value: Math.min(1, fEnd.cal), unit: '-' },
        { key: 'cycle_life_eff', label: 'Effective cycle life at this duty', value: fEnd.N80, unit: 'cycles', note: `Database rating ${ch.cycles_80} cycles at reference conditions` },
        { key: 'tr_cells_triggered', label: 'Cells driven into runaway', value: nTrig, unit: '', status: nTrig <= 1 ? 'ok' : 'warn', note: 'Target: no propagation beyond the initiating cell (conservative screening model)' },
        { key: 'tr_propagation_s', label: 'Time from first to last triggered cell', value: tLast, unit: 's' },
        { key: 'tr_front_cells_min', label: 'Propagation rate', value: tLast > 0 ? ((nTrig - 1) * 60) / tLast : 0, unit: 'cells/min' },
        { key: 'tr_T_peak_K', label: 'Peak cell temperature', value: Tpk, unit: 'K' },
        { key: 'tr_heat_cell_kJ', label: 'Heat released per cell', value: E / 1e3, unit: 'kJ' },
      ],
      plots: [
        { type: 'line', title: 'Capacity retention projection', xlabel: 'Time [yr]', ylabel: 'Capacity / new capacity [-]', series: [{ name: 'Total', x: ts, y: soh }, { name: 'Cycling only', x: ts, y: fd.map((f) => Math.max(0, 1 - f.cyc)), style: 'dash' }, { name: 'Calendar only', x: ts, y: fd.map((f) => Math.max(0, 1 - f.cal)), style: 'dash' }], annotations: [{ y: 0.8, label: 'End of life' }] },
        { type: 'line', title: 'Thermal-runaway propagation along the cell row', xlabel: 'Time [s]', ylabel: 'Cell temperature [K]', series: show.map((j) => ({ name: `Cell ${j + 1}`, x: idx.map((k) => r.t[k]), y: idx.map((k) => r.y[k][j]) })), annotations: [{ y: i.T_onset_K, label: 'Onset' }] },
      ],
      tables: [{ title: 'Runaway trigger times', columns: ['Cell', 'Trigger time [s]'], rows: trig.map((t, j) => [j + 1, Number.isFinite(t) ? t : 'not triggered']) }],
      warnings, models: ['Semi-empirical capacity fade: linear cycle term with depth-of-discharge, rate and Arrhenius factors plus √t calendar term (generic coefficients)', 'Lumped cell chain with one-step Arrhenius heat release for thermal-runaway propagation'],
      assumptions: ['Cycle life from the database is taken at 80 % depth of discharge, 1C and 25 °C; the database values are unverified and optimistic for power cells (a high-power 21700 cell datasheet guarantees ≥ 80 % after only 500 cycles at 1C)', 'Ageing coefficients (activation energy, calendar fade, rate and depth-of-discharge stress) and the 150 °C runaway onset are generic, unsourced values: calibrate them before relying on the result', 'Cycle and calendar fade add linearly; resistance growth and knee-point behaviour are not modelled', 'Runaway: uniform cell temperature, no vent-gas combustion, ejecta or electrical short paths; apparent activation temperature 15 600 K (generic)', 'The first cell is assumed already 30 K above onset', 'All of the runaway heat stays in the cell bodies, whereas in tests much of it leaves with the vented gas and ejecta: cell temperatures and propagation are over-predicted, so the model can show that a design is robust but not that it will fail', 'The default cell-to-cell conductance represents cells in direct contact with no barrier (illustrative)'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps in the runaway simulation', levels: [600, 1200, 2400, 4800], metric: 'tr_T_peak_K' },
  calibration: { params: [{ key: 'k_cal', min: 0, max: 0.2 }, { key: 'k_dod', min: 0, max: 3 }, { key: 'k_rate', min: 0, max: 2 }], sweep: 'years', target: 'soh_end', note: 'Supply measured capacity retention versus time from cell ageing tests at the intended duty.' },
  verify() {
    const base = { n: 1, C: 70, G: 0, hA: 0.1, Tamb: 300, Tonset: 423, E: 0, tau: 5, EaR: 15600, T0: 400, tEnd: 700, nSteps: 700 };
    const cool = trSim(base), adi = trSim({ ...base, hA: 0, E: 70000, T0: 460, tEnd: 200, nSteps: 4000 });
    const i0 = { ...Object.fromEntries(ageing.inputs.map((f) => [f.key, f.default])), cycles_yr: 0 }, ch = chemOf(i0.chem);
    const i1 = { ...i0, cycles_yr: ch.cycles_80, k_cal: 0, T_cell_K: 298.15, dod: 0.8, c_rate: 1 };
    return [
      N.check('Inert cell cooling follows exp(−hA·t/C)', cool.y[700][0], 300 + 100 * Math.exp(-1), 1e-8, 'Lumped-capacitance solution'),
      N.check('Adiabatic runaway temperature rise = E/(m·cp)', adi.y[4000][0], 460 + 1000, 1e-4, 'Energy conservation'),
      N.check('Calendar fade grows with √t', fade(i0, 4).cal / fade(i0, 1).cal, 2, 1e-12, 'Parabolic growth law'),
      N.check('Reference duty reaches 80 % at the rated cycle life', 1 - fade(i1, 1).cyc, 0.8, 1e-12, 'Definition of rated cycle life'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.tr_cells_triggered > 1) out.push({ severity: 'advise', title: 'Containment of a single-cell thermal runaway is not shown', detail: `${o.tr_cells_triggered} of ${Math.round(i.n_cells)} cells are triggered within ${o.tr_propagation_s.toFixed(0)} s with ${i.G_W_K} W/K between cells and ${i.hA_W_K} W/K of cooling per cell. The lumped model keeps all the released heat in the cells, so it errs on the side of propagation: it can show that a design is robust, not that it will fail.`, action: 'Enter the real thermal paths of the pack; add inter-cell thermal barriers (lower the cell-to-cell conductance), spacing or heat sinking and vent paths until the screen passes, then demonstrate non-propagation by test.', basis: 'Design objective: containment of a single-cell failure, to be shown by test; conservative lumped energy balance' });
    else out.push({ severity: 'info', title: 'No propagation predicted with these thermal paths', detail: 'Only the initiating cell goes into runaway in this lumped model.', action: 'Confirm by abuse testing: the model omits vent-gas combustion and electrical short paths.', basis: 'Lumped energy balance' });
    if (o.batt_life_yr < 5) out.push({ severity: 'advise', title: 'Short pack life', detail: `${o.batt_life_yr.toFixed(1)} years (${o.cycles_to_eol.toFixed(0)} cycles) to 80 % capacity; cycling contributes ${(100 * o.fade_cycle).toFixed(0)} points and storage ${(100 * o.fade_calendar).toFixed(0)} points at the horizon.`, action: 'Reduce depth of discharge, charge rate and cell temperature, and store at moderate state of charge; each replacement pack carries embodied emissions and cost (Suite 26).', basis: 'Semi-empirical ageing law (generic coefficients — calibrate)' });
    return out;
  },
};

// ---- permanent-magnet synchronous machine ---------------------------------------------------
/** Machine parameters from ratings: flux linkage and inductances chosen so the voltage limit is met at base speed and rated torque. */
function pmsmDesign(i) {
  const wb = (i.rpm_base * Math.PI) / 30, p = Math.round(i.pole_pairs), Vmax = (i.m_max * i.Vdc) / Math.sqrt(3), Tr = i.P_rated_W / wb;
  const psi = (Vmax * (1 - i.f_cu)) / (p * wb * Math.hypot(1, i.saliency * i.l_pu)), Imax = Tr / (1.5 * p * psi), Ld = (i.l_pu * psi) / Imax;
  return { wb, p, Vmax, Tr, psi, Imax, Ld, Lq: i.saliency * Ld, Rs: (i.f_cu * i.P_rated_W) / (1.5 * Imax * Imax), Pr: i.P_rated_W, f_fe: i.f_fe, f_mech: i.f_mech, flux0: psi * Math.hypot(1, i.saliency * i.l_pu) };
}
/** Steady-state dq operating point for shaft torque T at speed w: MTPA when the voltage allows, otherwise field weakening on the voltage limit. */
function pmsmOp(m, T, w) {
  const we = m.p * w, wr = w / m.wb, Pmech = m.f_mech * m.Pr * (0.4 * wr + 0.6 * wr ** 3), Tem = T + (w > 1e-9 ? Pmech / w : 0);
  const iqOf = (id) => Tem / (1.5 * m.p * (m.psi + (m.Ld - m.Lq) * id));
  const vdq = (id) => { const iq = iqOf(id); return [m.Rs * id - we * m.Lq * iq, m.Rs * iq + we * (m.Ld * id + m.psi)]; };
  const dv = (id) => { const v = vdq(id); return v[0] * v[0] + v[1] * v[1] - m.Vmax * m.Vmax; };
  let id = m.Lq > m.Ld * (1 + 1e-9) ? N.goldenSection((x) => x * x + iqOf(x) ** 2, -m.Imax, 0, 1e-9) : 0, ok = true;
  if (dv(id) > 0) { // move along the constant-torque locus to more negative id until the voltage limit is met
    let x0 = id, found = false;
    for (let k = 1; k <= 48; k++) { const x1 = id + ((-m.Imax - id) * k) / 48; if (dv(x1) <= 0) { id = N.brent(dv, x1, x0, 1e-12 * m.Imax); found = true; break; } x0 = x1; }
    if (!found) ok = false;
  }
  const iq = iqOf(id), I = Math.hypot(id, iq), v = vdq(id), V = Math.hypot(v[0], v[1]);
  if (I > m.Imax * (1 + 1e-9)) ok = false;
  const Pcu = 1.5 * m.Rs * I * I, flux = Math.hypot(m.psi + m.Ld * id, m.Lq * iq) / m.flux0, Pfe = m.f_fe * m.Pr * (0.5 * wr + 0.5 * wr * wr) * flux * flux;
  const Pout = T * w, Pin = Pout + Pcu + Pfe + Pmech;
  return { ok, id, iq, I, V, Pcu, Pfe, Pmech, Pout, Pin, eff: Pin > 0 ? Pout / Pin : 0, pf: N.clamp((Pin - Pfe) / Math.max(1.5 * V * I, 1e-12), 0, 1) };
}
/** Largest feasible shaft torque at speed w (bisection on feasibility). */
function pmsmTmax(m, w) { let lo = 0, hi = 1.6 * m.Tr; for (let k = 0; k < 44; k++) { const md = 0.5 * (lo + hi); if (pmsmOp(m, md, w).ok) lo = md; else hi = md; } return lo; }
/** Two-level inverter losses: conduction (threshold + resistive, equal transistor/diode characteristics) and switching. */
function inverterLoss(i, m, Ipk) { const cond = (6 * i.V_th * Ipk) / Math.PI + 1.5 * (i.V_on_rated / m.Imax) * Ipk * Ipk, sw = (6 * i.f_sw_kHz * 1e3 * i.k_sw * i.Vdc * Ipk) / Math.PI; return { cond, sw, total: cond + sw }; }

const motor = {
  id: 'motor', title: 'PMSM drive: torque–speed envelope and efficiency map', fidelity: 'reduced-order',
  summary: 'Steady-state dq-axis model of a permanent-magnet synchronous machine with current and voltage limits, field weakening, copper, iron and mechanical losses, and an inverter conduction and switching loss model.',
  equations: ['Park transformation equations (dq steady state)', 'Electrical machine equations', 'Faraday’s law of electromagnetic induction (back-EMF)', 'Power balance equations', 'Converter switching equations (averaged loss form)'],
  inputs: [
    { key: 'P_rated_W', label: 'Rated shaft power', unit: 'W', default: 150000, min: 10, group: 'Machine' },
    { key: 'rpm_base', label: 'Base (corner) speed', unit: 'rpm', default: 1100, min: 10, max: 200000, group: 'Machine' },
    { key: 'speed_ratio', label: 'Maximum / base speed', unit: '-', default: 1.8, min: 1.05, max: 6, group: 'Machine' },
    { key: 'pole_pairs', label: 'Pole pairs', unit: '', default: 8, min: 1, max: 60, step: 1, discrete: true, group: 'Machine' },
    { key: 'l_pu', label: 'Per-unit d-axis inductance Ld·Imax/ψ', unit: '-', default: 0.5, min: 0.02, max: 1.5, group: 'Machine', help: '1 gives an ideal infinite constant-power range; 0.2–0.6 is typical' },
    { key: 'saliency', label: 'Saliency ratio Lq/Ld', unit: '-', default: 1, min: 1, max: 4, group: 'Machine', help: '1 for surface magnets; 1.5–3 for interior magnets (adds reluctance torque)' },
    { key: 'f_cu', label: 'Copper loss fraction at rated point', unit: '-', default: 0.025, min: 1e-9, max: 0.2, group: 'Losses' },
    { key: 'f_fe', label: 'Iron loss fraction at rated point', unit: '-', default: 0.015, min: 0, max: 0.2, group: 'Losses', help: 'Split equally between hysteresis (∝ f) and eddy (∝ f²), scaled with flux²' },
    { key: 'f_mech', label: 'Friction and windage fraction at base speed', unit: '-', default: 0.005, min: 0, max: 0.1, group: 'Losses' },
    { key: 'Vdc', label: 'DC-link voltage', unit: 'V', default: 800, min: 5, max: 3000, group: 'Inverter' },
    { key: 'm_max', label: 'Maximum modulation utilisation', unit: '-', default: 0.95, min: 0.5, max: 1, group: 'Inverter', help: 'Fraction of Vdc/√3 (space-vector limit) available as peak phase voltage' },
    { key: 'f_sw_kHz', label: 'Switching frequency', unit: 'kHz', default: 10, min: 1, max: 200, group: 'Inverter' },
    { key: 'V_th', label: 'Device threshold voltage', unit: 'V', default: 0.8, min: 0, max: 3, group: 'Inverter', help: '0 for MOSFET bridges, 0.7–1.2 V for IGBT/diode' },
    { key: 'V_on_rated', label: 'Resistive on-state drop at rated current', unit: 'V', default: 1, min: 0, max: 5, group: 'Inverter' },
    { key: 'k_sw', label: 'Switching energy coefficient', unit: 'J/(V·A)', default: 1e-7, min: 0, max: 2e-6, group: 'Inverter', help: '(Eon+Eoff+Err)/(V·I): ≈ 0.5e-7 SiC, ≈ 3e-7 Si IGBT. Take from the device data sheet' },
    { key: 'P_op_W', label: 'Operating-point shaft power', unit: 'W', default: 83000, min: 0, group: 'Operating point' },
    { key: 'rpm_op', label: 'Operating-point speed', unit: 'rpm', default: 900, min: 1, group: 'Operating point' },
    { key: 'nMap', label: 'Map resolution (speed points)', unit: '', default: 24, min: 6, max: 96, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const e = est(c, up, d);
    if (!e.electric) { const P = Math.max(100, (c.systems.gen_kVA || 0) * 900), hv = c.systems.bus_V > 60; return { P_rated_W: P, rpm_base: 12000, pole_pairs: 3, Vdc: hv ? 270 : 28, V_th: hv ? 0.8 : 0, V_on_rated: hv ? 1 : 0.3, f_sw_kHz: 20, P_op_W: 0.6 * P, rpm_op: 12000 }; }
    const P = c.prop.P0_W || (c.systems.motor_kW || 0) * 1e3 || undefined, rpm = c.prop.rpm || c.rotor.rpm || undefined, Pop = (e.Ph > 0 ? e.Ph : e.Pc) / e.n || undefined, lv = c.systems.bus_V < 100;
    return { P_rated_W: P, rpm_base: rpm, Vdc: c.systems.bus_V || undefined, P_op_W: Pop, rpm_op: P && rpm && Pop ? rpm * N.clamp((Pop / P) ** (1 / 3), 0.3, 1.2) : undefined, pole_pairs: lv ? 7 : 8, V_th: lv ? 0 : 0.8, V_on_rated: lv ? 0.3 : 1, f_sw_kHz: lv ? 24 : 10 };
  },
  run(i) {
    const m = pmsmDesign(i), wMax = i.speed_ratio * m.wb, nx = Math.round(i.nMap), ny = Math.max(6, Math.round(0.8 * nx)), rpm = (w) => (w * 30) / Math.PI, warnings = [];
    const ws = N.linspace(0.08 * m.wb, wMax, nx), Tenv = ws.map((w) => pmsmTmax(m, w)), fr = N.linspace(0.08, 1, ny);
    const effMap = fr.map((f) => ws.map((w, k) => { const o = pmsmOp(m, f * Tenv[k], w), inv = inverterLoss(i, m, o.I).total; return Tenv[k] > 0 ? (100 * o.Pout) / (o.Pin + inv) : 0; }));
    const mEff = fr.map((f) => ws.map((w, k) => (Tenv[k] > 0 ? pmsmOp(m, f * Tenv[k], w).eff : 0))), effPk = N.amax(mEff.map(N.amax));
    // operating point (clamped to the envelope when outside)
    const wOp = Math.min((i.rpm_op * Math.PI) / 30, wMax), TmaxOp = pmsmTmax(m, wOp); let Top = i.P_op_W / wOp;
    if ((i.rpm_op * Math.PI) / 30 > wMax * 1.0001) warnings.push('The operating speed exceeds the maximum speed of the machine: it was limited to the maximum.');
    if (Top > TmaxOp) { warnings.push(`The operating point (${Top.toFixed(1)} N·m) lies outside the torque–speed envelope (${TmaxOp.toFixed(1)} N·m at this speed): results are shown at the envelope.`); Top = TmaxOp; }
    const o = pmsmOp(m, Top, wOp), inv = inverterLoss(i, m, o.I), invEff = o.Pin > 0 ? o.Pin / (o.Pin + inv.total) : 1, PmaxEnv = N.amax(ws.map((w, k) => w * Tenv[k]));
    const fe = (m.p * wOp) / (2 * Math.PI);
    if (i.f_sw_kHz * 1e3 < 10 * ((m.p * wMax) / (2 * Math.PI))) warnings.push('Switching frequency is below ten times the maximum electrical frequency: current ripple and harmonic losses (not modelled) become significant.');
    return {
      kpis: [
        { key: 'motor_eff', label: 'Machine efficiency at the operating point', value: o.eff, unit: '-', status: o.eff > 0.9 ? 'ok' : 'warn' },
        { key: 'inverter_eff', label: 'Inverter efficiency at the operating point', value: invEff, unit: '-' },
        { key: 'drive_eff', label: 'Drive efficiency (inverter × machine)', value: o.eff * invEff, unit: '-' },
        { key: 'drive_losses_W', label: 'Drive losses at the operating point', value: o.Pin - o.Pout + inv.total, unit: 'W' },
        { key: 'T_op_Nm', label: 'Operating torque', value: Top, unit: 'N·m' },
        { key: 'torque_margin_pct', label: 'Torque margin to the envelope', value: TmaxOp > 0 ? 100 * (1 - Top / TmaxOp) : 0, unit: '%', status: Top < 0.999 * TmaxOp ? 'ok' : 'bad' },
        { key: 'I_phase_rms_A', label: 'Phase current (RMS)', value: o.I / Math.SQRT2, unit: 'A' },
        { key: 'id_A', label: 'd-axis current', value: o.id, unit: 'A', note: 'Negative in field weakening' },
        { key: 'power_factor', label: 'Fundamental power factor', value: o.pf, unit: '-' },
        { key: 'f_elec_Hz', label: 'Electrical frequency', value: fe, unit: 'Hz' },
        { key: 'T_rated_Nm', label: 'Rated torque', value: m.Tr, unit: 'N·m' },
        { key: 'P_env_max_W', label: 'Maximum shaft power on the envelope', value: PmaxEnv, unit: 'W' },
        { key: 'P_at_max_speed_W', label: 'Shaft power at maximum speed', value: wMax * Tenv[nx - 1], unit: 'W' },
        { key: 'eff_peak', label: 'Peak machine efficiency on the map', value: effPk, unit: '-' },
      ],
      plots: [
        { type: 'heat', title: 'Drive efficiency map (machine × inverter)', xlabel: 'Speed [rpm]', ylabel: 'Torque [% of envelope torque at that speed]', zlabel: 'Efficiency [%]', x: ws.map(rpm), y: fr.map((f) => 100 * f), z: effMap, contours: 12, overlay: [{ name: 'Operating point', x: [rpm(wOp)], y: [TmaxOp > 0 ? (100 * Top) / TmaxOp : 0] }] },
        { type: 'line', title: 'Torque–speed envelope', xlabel: 'Speed [rpm]', ylabel: 'Torque [N·m]', series: [{ name: 'Maximum torque', x: ws.map(rpm), y: Tenv }, { name: 'Operating point', x: [rpm(wOp)], y: [Top], style: 'points' }], annotations: [{ x: i.rpm_base, label: 'Base speed' }] },
        { type: 'line', title: 'Power–speed envelope', xlabel: 'Speed [rpm]', ylabel: 'Shaft power [kW]', series: [{ name: 'Maximum power', x: ws.map(rpm), y: ws.map((w, k) => (w * Tenv[k]) / 1e3) }] },
        { type: 'bar', title: 'Loss breakdown at the operating point', ylabel: 'Loss [W]', categories: ['Copper', 'Iron', 'Friction & windage', 'Inverter conduction', 'Inverter switching'], series: [{ name: 'Loss', y: [o.Pcu, o.Pfe, o.Pmech, inv.cond, inv.sw] }] },
      ],
      tables: [{ title: 'Equivalent machine parameters (derived from ratings)', columns: ['Parameter', 'Value', 'Unit'], rows: [['Magnet flux linkage ψ', m.psi, 'V·s'], ['Ld', m.Ld, 'H'], ['Lq', m.Lq, 'H'], ['Phase resistance', m.Rs, 'Ω'], ['Peak phase current limit', m.Imax, 'A'], ['Peak phase voltage limit', m.Vmax, 'V']] }],
      outputs: { drive_heat_W: o.Pin - o.Pout + inv.total },
      warnings, models: ['dq-axis PMSM steady-state model with MTPA and voltage-limited field weakening', 'Loss separation: copper I²R, iron (hysteresis ∝ f, eddy ∝ f², flux² scaling), friction and windage', 'Averaged two-level inverter conduction and switching losses'],
      assumptions: ['Machine parameters are back-calculated from ratings, not from an electromagnetic design', 'Constant inductances (no magnetic saturation or cross-coupling) and constant magnet flux (no temperature effect)', 'Iron loss treated as an added input power, not as a braking torque', 'Sinusoidal currents: PWM harmonic and AC winding losses neglected', 'Generator operation of conventional aircraft machines is represented by the equivalent motoring point', 'Loss fractions, per-unit inductance and inverter device parameters are typical values: take them from machine and device data'],
    };
  },
  convergence: { param: 'nMap', label: 'Map speed points', levels: [8, 16, 32, 64], metric: 'eff_peak' },
  calibration: { params: [{ key: 'f_cu', min: 1e-4, max: 0.2 }, { key: 'f_fe', min: 0, max: 0.2 }, { key: 'f_mech', min: 0, max: 0.1 }], sweep: 'P_op_W', target: 'motor_eff', note: 'Supply dynamometer efficiency versus shaft power at a fixed speed (repeat at several speeds).' },
  verify() {
    const i = { ...Object.fromEntries(motor.inputs.map((f) => [f.key, f.default])), f_cu: 1e-9, f_fe: 0, f_mech: 0, l_pu: 0.6 }, m = pmsmDesign(i), w = 2 * m.wb, we = m.p * w;
    const idx = ((m.Vmax / we) ** 2 - m.psi ** 2 - (m.Ld * m.Imax) ** 2) / (2 * m.psi * m.Ld), Tex = 1.5 * m.p * m.psi * Math.sqrt(m.Imax ** 2 - idx ** 2), o = pmsmOp(m, 0.5 * m.Tr, 0.5 * m.wb);
    return [
      N.check('Maximum torque below base speed equals rated torque', pmsmTmax(m, 0.5 * m.wb), m.Tr, 1e-8, 'T = 1.5·p·ψ·Imax (surface-magnet machine, id = 0)'),
      N.check('Field-weakening torque at twice base speed', pmsmTmax(m, w), Tex, 1e-6, 'Intersection of the current circle and the voltage circle (lossless)'),
      N.check('Lossless machine: electrical input = shaft power', o.Pin, 1.5 * m.p * 0.5 * m.wb * m.psi * o.iq, 1e-8, 'P = 1.5·ωe·ψ·iq'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.torque_margin_pct < 10) out.push({ severity: o.torque_margin_pct <= 0.1 ? 'critical' : 'warn', title: 'Operating point is at the edge of the envelope', detail: `Torque margin ${o.torque_margin_pct.toFixed(1)}%.`, action: 'Raise the DC-link voltage or rated power, lower the operating speed, or accept a deeper field-weakening design (higher per-unit inductance).', basis: 'Current and voltage limits of the dq model' });
    if (o.drive_eff < 0.9) out.push({ severity: 'advise', title: 'Drive efficiency is below 90 %', detail: `Machine ${(100 * o.motor_eff).toFixed(1)}% × inverter ${(100 * o.inverter_eff).toFixed(1)}% at ${(i.P_op_W / 1e3).toPrecision(3)} kW; ${o.drive_losses_W.toFixed(0)} W must be rejected as heat.`, action: 'Move the main operating point towards the high-efficiency island (match propeller speed and gearing), and consider SiC devices for part-load switching loss. Each point of efficiency is a point of battery energy.', basis: 'Loss breakdown of this run' });
    else out.push({ severity: 'info', title: 'Drive efficiency', detail: `${(100 * o.drive_eff).toFixed(1)}% at the operating point; peak machine efficiency ${(100 * o.eff_peak).toFixed(1)}%.`, action: 'Pass the drive heat load to the thermal suite and use the drive efficiency in mission energy calculations.', basis: 'Loss breakdown of this run' });
    return out;
  },
};

// ---- distribution network -------------------------------------------------------------------
/** Newton solution of the nodal equations G·V = I_src − P/V with constant-power loads. */
function loadFlow(G, Is, P, V0) {
  const n = Is.length; let V = new Array(n).fill(V0), ok = true;
  const F = (v) => v.map((vi, r) => N.dot(G[r], v) - Is[r] + P[r] / vi);
  for (let it = 0; it < 60; it++) {
    const dV = N.solve(G.map((row, r) => row.map((g, cc) => (r === cc ? g - P[r] / (V[r] * V[r]) : g))), F(V).map((x) => -x));
    V = V.map((v, r) => Math.max(0.02 * V0, v + dV[r]));
    if (N.amax(dV.map(Math.abs)) < 1e-11 * V0) break;
    if (it === 59) ok = false;
  }
  if (N.amin(V) < 0.5 * V0) ok = false;
  return { V, res: N.amax(F(V).map(Math.abs)), ok };
}
const NODES = ['Main bus 1', 'Main bus 2', 'Essential bus', 'Power feeder 1', 'Power feeder 2', 'Avionics bus'];
const PHASES = ['Ground', 'Take-off / hover', 'Climb', 'Cruise', 'Descent', 'Landing'];
// load groups: [name, node, share of connected system load, priority (1 essential … 3 sheddable), demand factor per phase] — generic factors
const GROUPS = [
  ['Avionics and flight controls', 5, 0.1, 1, [0.8, 1, 1, 1, 1, 1]], ['Essential lighting, fuel and engine systems', 2, 0.1, 1, [0.6, 1, 1, 0.9, 0.9, 1]],
  ['Actuation and hydraulic pumps', 3, 0.2, 2, [0.3, 1, 0.6, 0.4, 0.5, 1]], ['Ice and rain protection', 3, 0.15, 2, [0, 0.5, 1, 0.6, 1, 0.5]],
  ['Environmental control and galley', 4, 0.3, 3, [0.8, 0.3, 0.7, 1, 0.8, 0.3]], ['Cabin, lighting and miscellaneous', 4, 0.15, 3, [1, 0.6, 0.8, 1, 0.8, 0.6]],
];
/** Loads [W] for flight phase k: propulsion split over the two power feeders (priority 1) plus the system groups. */
function netLoads(i, k) {
  const pp = [0.03 * i.P_prop_peak_W, i.P_prop_peak_W, Math.min(i.P_prop_peak_W, 1.5 * i.P_prop_cruise_W), i.P_prop_cruise_W, 0.4 * i.P_prop_cruise_W, i.vtol ? i.P_prop_peak_W : 0.5 * i.P_prop_cruise_W][k];
  const L = GROUPS.map((g) => ({ name: g[0], node: g[1], P: i.P_sys_W * g[2] * g[4][k], prio: g[3] }));
  if (i.P_prop_peak_W > 0) L.unshift({ name: 'Propulsion drive A', node: 3, P: pp / 2, prio: 1, prop: true }, { name: 'Propulsion drive B', node: 4, P: pp / 2, prio: 1, prop: true });
  return L;
}
/** Branch resistances from a design current density (so every feeder has the same design voltage drop per metre). */
function netBranches(i) {
  const rho = RHO_CU * (1 + ALPHA_CU * 60) * i.loop_factor * i.J_Amm2 * 1e6;
  const pk = N.range(6, (nd) => N.amax(N.range(6, (k) => N.sum(netLoads(i, k).filter((l) => l.node === nd).map((l) => l.P)))));
  const R = (L, P) => (rho * L * i.V_bus) / Math.max(P, 1e-6 * i.P_src_W + 1e-9), Lf = i.L_feeder_m;
  return [[0, 1, R(0.3 * Lf, i.P_src_W), 'Bus tie'], [0, 2, R(0.3 * Lf, pk[2] + pk[5] || i.P_batt_W), 'Essential feeder'], [0, 3, R(Lf, pk[3]), 'Power feeder 1'], [1, 4, R(Lf, pk[4]), 'Power feeder 2'], [2, 5, R(0.4 * Lf, pk[5]), 'Avionics feeder']];
}
/** Solve one network state. st: {s1, s2, batt, tie, loads, fault (node index or −1)}. */
function netSolve(i, st) {
  const V0 = i.V_bus, G = N.zeros(6), Is = new Array(6).fill(0), P = new Array(6).fill(0), br = netBranches(i), flt = st.fault >= 0;
  for (let k = 0; k < 6; k++) G[k][k] = 1e-12;
  for (const [a, b, R] of br) { if (a === 0 && b === 1 && !st.tie) continue; G[a][a] += 1 / R; G[b][b] += 1 / R; G[a][b] -= 1 / R; G[b][a] -= 1 / R; }
  const Rs = (flt ? V0 / i.fault_mult : (i.src_droop_pct / 100) * V0) * (V0 / i.P_src_W), Rb = (flt ? 0.02 : 0.05) * V0 * (V0 / i.P_batt_W);
  const src = [[0, st.s1, V0, Rs], [1, st.s2, V0, Rs], [2, st.batt, 0.95 * V0, Rb]].filter((s) => s[1]);
  for (const [nd, , Vs, R] of src) { G[nd][nd] += 1 / R; Is[nd] += Vs / R; }
  if (flt) G[st.fault][st.fault] += 1e9; else for (const l of st.loads) P[l.node] += l.P;
  const lf = flt ? { V: N.solve(G, Is), res: 0, ok: true } : loadFlow(G, Is, P, V0), V = lf.V;
  const cable = N.sum(br.map(([a, b, R]) => (a === 0 && b === 1 && !st.tie ? 0 : (V[a] - V[b]) ** 2 / R)));
  const sI = src.map(([nd, , Vs, R]) => (Vs - V[nd]) / R), Pdel = N.sum(src.map(([nd], k) => V[nd] * sI[k])), srcLoss = N.sum(src.map(([, , , R], k) => sI[k] * sI[k] * R));
  return { ...lf, cable, srcLoss, Pdel, Pload: N.sum(P), Ifault: flt ? V[st.fault] * 1e9 : 0 };
}
/** Shed whole loads by priority (3 first, then 2) until demand fits the capacity; then power-limit priority-1 loads if still above. */
function shed(loads, cap) {
  const kept = loads.map((l) => ({ ...l })), dropped = []; let tot = N.sum(kept.map((l) => l.P));
  for (const pr of [3, 2]) for (let k = kept.length - 1; k >= 0 && tot > cap; k--) if (kept[k].prio === pr && kept[k].P > 0) { tot -= kept[k].P; dropped.push(kept[k].name); kept[k].P = 0; }
  const scale = tot > cap ? cap / tot : 1; if (scale < 1) kept.forEach((l) => { l.P *= scale; });
  return { kept, dropped, shedW: N.sum(loads.map((l) => l.P)) - N.sum(kept.map((l) => l.P)), scale };
}
/** Largest drop between a source bus and the loads it feeds [V]; with the tie closed everything hangs on the higher bus. */
const feederDrop = (V, tie) => (tie ? Math.max(V[0], V[1]) - N.amin(V) : Math.max(V[0] - Math.min(V[2], V[3], V[5]), V[1] - V[4]));
const BREAKERS = [1, 2, 3, 5, 7.5, 10, 15, 20, 25, 35, 50, 75, 100, 150, 200, 300, 400, 600, 800, 1000, 1500, 2000];

const network = {
  id: 'network', title: 'Power distribution: load flow, faults, load shedding and load analysis', fidelity: 'numerical',
  summary: 'Nodal (Kirchhoff) analysis of a two-channel bus network with feeders and constant-power loads: bus voltages and losses in every flight phase, fault currents, protection ratings, priority load shedding after loss of a source, and the electrical load analysis with capacity margins.',
  equations: ['Kirchhoff’s current law', 'Kirchhoff’s voltage law', 'Ohm’s law', 'Power balance equations'],
  inputs: [
    { key: 'V_bus', label: 'Bus voltage', unit: 'V', default: 115, min: 5, max: 3000, group: 'Architecture', help: 'AC systems are analysed as a DC-equivalent real-power network' },
    { key: 'n_src', label: 'Main sources (generators or battery strings)', unit: '', default: 2, min: 1, max: 2, step: 1, discrete: true, group: 'Architecture' },
    { key: 'P_src_W', label: 'Rating of each main source', unit: 'W', default: 81000, min: 1, group: 'Architecture' },
    { key: 'P_batt_W', label: 'Emergency battery power', unit: 'W', default: 5000, min: 1, group: 'Architecture' },
    { key: 'src_droop_pct', label: 'Source voltage droop at rated power', unit: '%', default: 0.5, min: 0.01, max: 25, group: 'Architecture', help: 'Regulated generator ≈ 0.5 %; battery string = I·R/V at rated power' },
    { key: 'fault_mult', label: 'Source fault current / rated current', unit: '-', default: 3, min: 1, max: 200, group: 'Protection', help: '≈ 3 for regulated aircraft generators; 10–50 for batteries' },
    { key: 'eta_src', label: 'Source conversion efficiency', unit: '-', default: 0.9, min: 0.5, max: 1, group: 'Architecture', help: 'Generator and rectifier; 1 for battery strings (their loss is the droop)' },
    { key: 'P_sys_W', label: 'Connected system (non-propulsive) load', unit: 'W', default: 97000, min: 0, group: 'Loads' },
    { key: 'P_prop_peak_W', label: 'Propulsion electrical demand, peak', unit: 'W', default: 0, min: 0, group: 'Loads', help: '0 for conventionally powered aircraft' },
    { key: 'P_prop_cruise_W', label: 'Propulsion electrical demand, cruise', unit: 'W', default: 0, min: 0, group: 'Loads' },
    { key: 'vtol', label: 'Peak propulsion power again at landing', type: 'bool', default: false, group: 'Loads' },
    { key: 'eta_drive', label: 'Propulsion drive efficiency', unit: '-', default: 0.93, min: 0.5, max: 1, group: 'Loads', help: 'Inverter × motor, for the drive heat included in total losses' },
    { key: 'L_feeder_m', label: 'Power feeder length', unit: 'm', default: 15, min: 0.1, max: 200, group: 'Cables' },
    { key: 'J_Amm2', label: 'Design current density', unit: 'A/mm²', default: 4, min: 0.5, max: 15, group: 'Cables', help: 'Conductors are sized at this density for their peak connected load. Bundled copper aircraft wire carries about 4–8 A/mm² at mid gauge but only 2.4–4.6 A/mm² in large feeders (AC 43.13-1B Table 11-9)' },
    { key: 'loop_factor', label: 'Conductor length factor', unit: '-', default: 2, min: 1, max: 2, group: 'Cables', help: '2 for DC two-wire (composite airframe); 1 for metallic structure return or the per-phase equivalent of a three-phase AC feeder' },
  ],
  defaults: (c, up, d) => {
    const e = est(c, up, d), V = c.systems.bus_V || undefined, L = 0.5 * e.len;
    if (!e.electric) return { V_bus: V, n_src: Math.min(2, e.n), P_src_W: e.genTot / Math.min(2, e.n), P_batt_W: 1.5 * e.Pess, P_sys_W: 0.6 * e.genTot, L_feeder_m: L, loop_factor: V > 60 ? 1 : 2 };
    const ch = chemOf(c.systems.batt_chem), Pk = Math.max(e.Ph, e.Pc) / 0.93, Psys = Math.max(10, 0.02 * Pk), Pstr = 0.65 * (Pk + Psys), E = (c.systems.batt_kWh || 1) * 1e3;
    const droop = N.clamp((100 * 2 * ch.r_mohm_ah * 1e-3 * Pstr) / (E * ch.v_nom), 0.2, 20); // string of half the pack: R·P/V² = 2·r·P/(E·v_cell)
    return { V_bus: V, n_src: 2, P_src_W: Pstr, P_batt_W: Math.max(20, 0.6 * Psys), src_droop_pct: droop, fault_mult: N.clamp(100 / droop, 3, 200), eta_src: 1, P_sys_W: Psys, P_prop_peak_W: Pk, P_prop_cruise_W: e.Pc / 0.93, vtol: e.Ph > 0, L_feeder_m: L };
  },
  run(i) {
    const two = Math.round(i.n_src) >= 2, cap = (two ? 2 : 1) * i.P_src_W, normal = { s1: true, s2: two, batt: false, tie: !two, fault: -1 }, warnings = [];
    const ph = PHASES.map((_, k) => { const loads = netLoads(i, k), r = netSolve(i, { ...normal, loads }), tot = N.sum(loads.map((l) => l.P)); return { loads, r, tot, prop: N.sum(loads.filter((l) => l.prop).map((l) => l.P)) }; });
    const kPk = N.argmax(ph.map((p) => p.tot)), Vmin = N.amin(ph.map((p) => N.amin(p.r.V))), cr = ph[3], margin = 100 * (1 - ph[kPk].tot / cap), dFeed = (100 * N.amax(ph.map((p) => feederDrop(p.r.V, normal.tie)))) / i.V_bus;
    const losses = cr.r.cable + cr.r.srcLoss + cr.r.Pdel * (1 / i.eta_src - 1) + cr.prop * (1 - i.eta_drive);
    // degraded: one main source lost at the peak-demand phase (or all generation lost when there is only one), then battery only
    const capB = two ? i.P_src_W : i.P_batt_W, sB = shed(ph[kPk].loads, capB), rB = netSolve(i, { s1: false, s2: two, batt: !two, tie: true, fault: -1, loads: sB.kept });
    const ess = ph[3].loads.filter((l) => l.prio === 1 && !l.prop), sC = shed(ess, i.P_batt_W), rC = netSolve(i, { s1: false, s2: false, batt: true, tie: true, fault: -1, loads: sC.kept });
    const fBus = netSolve(i, { ...normal, fault: 0, loads: [] }).Ifault, fFd = netSolve(i, { ...normal, fault: 3, loads: [] }).Ifault;
    const resid = N.amax(ph.map((p) => p.r.res)), balErr = Math.abs(cr.r.Pdel - cr.r.Pload - cr.r.cable) / Math.max(cr.r.Pdel, 1e-9);
    if (ph.some((p) => !p.r.ok)) warnings.push('The load flow does not converge to an acceptable voltage in at least one phase: the feeders or sources cannot carry the constant-power demand (voltage collapse).');
    if (margin < 0) warnings.push(`Peak demand exceeds installed source capacity by ${(-margin).toFixed(0)}% in "${PHASES[kPk]}".`);
    if (sB.scale < 1) warnings.push(`After losing one source, shedding all non-essential loads is not enough: essential and propulsion loads are power-limited to ${(100 * sB.scale).toFixed(0)}%.`);
    if (i.P_prop_peak_W > 0) warnings.push('With all main sources lost the emergency battery supports essential systems only; propulsion is not sustained in that state.');
    if (sC.scale < 1) warnings.push(`Priority-1 loads exceed the emergency battery rating: only ${(100 * sC.scale).toFixed(0)}% can be supplied on battery alone, so a further emergency source (ram-air turbine, APU) or deeper shedding is needed.`);
    const rows = netLoads(i, kPk).filter((l) => l.P > 0).map((l) => { const I = l.P / i.V_bus, rt = BREAKERS.find((b) => b >= 1.25 * I) ?? Math.ceil((1.25 * I) / 500) * 500, If = l.node === 0 || l.node === 1 ? fBus : fFd; return [l.name, NODES[l.node], l.prio, l.P / 1e3, I, rt, If / rt]; });
    const pct = (r) => r.V.map((v) => (100 * v) / i.V_bus);
    return {
      kpis: [
        { key: 'bus_V_min', label: 'Minimum bus voltage, normal operation', value: Vmin, unit: 'V', status: dFeed <= 5 ? 'ok' : dFeed <= 10 ? 'warn' : 'bad', note: 'Lowest node in any flight phase, including source droop; criterion: feeder drop ≤ 5 % of nominal' },
        { key: 'v_drop_max_pct', label: 'Largest voltage drop below nominal, normal operation', value: 100 * (1 - Vmin / i.V_bus), unit: '%', note: 'Source droop plus feeder drop' },
        { key: 'v_drop_feeder_pct', label: 'Largest feeder drop from source bus to load', value: dFeed, unit: '%', status: dFeed <= 5 ? 'ok' : dFeed <= 10 ? 'warn' : 'bad' },
        { key: 'elec_losses_W', label: 'Electrical losses in cruise', value: losses, unit: 'W', note: 'Cables + source droop + source conversion + propulsion drives' },
        { key: 'cable_loss_cruise_W', label: 'Cable loss in cruise', value: cr.r.cable, unit: 'W' },
        { key: 'dist_eff', label: 'Distribution efficiency in cruise', value: cr.r.Pdel > 0 ? cr.r.Pload / (cr.r.Pdel + cr.r.srcLoss) : 1, unit: '-' },
        { key: 'P_demand_peak_W', label: `Peak demand (${PHASES[kPk]})`, value: ph[kPk].tot, unit: 'W' },
        { key: 'gen_margin_pct', label: 'Source capacity margin at peak demand', value: margin, unit: '%', status: margin > 10 ? 'ok' : margin > 0 ? 'warn' : 'bad' },
        { key: 'shed_W', label: 'Load shed after loss of one source', value: sB.shedW, unit: 'W' },
        { key: 'bus_V_min_degraded', label: 'Minimum bus voltage, one source lost', value: N.amin(rB.V), unit: 'V' },
        { key: 'bus_V_min_emergency', label: 'Minimum bus voltage, battery only', value: N.amin(rC.V), unit: 'V' },
        { key: 'ess_supported_pct', label: 'Priority-1 system load supported on battery alone', value: 100 * sC.scale, unit: '%', status: sC.scale >= 1 ? 'ok' : 'warn' },
        { key: 'fault_bus_A', label: 'Bolted fault current at main bus 1', value: fBus, unit: 'A' },
        { key: 'fault_feeder_A', label: 'Bolted fault current at the end of power feeder 1', value: fFd, unit: 'A' },
        { key: 'kcl_residual_A', label: 'Largest nodal current imbalance', value: resid, unit: 'A' },
        { key: 'power_balance_err', label: 'Power balance residual in cruise', value: balErr, unit: '-' },
      ],
      plots: [
        { type: 'bar', title: 'Node voltages by operating state', ylabel: 'Voltage [% of nominal]', categories: NODES, series: [{ name: `Normal, ${PHASES[kPk]}`, y: pct(ph[kPk].r) }, { name: 'One source lost (after shedding)', y: pct(rB) }, { name: 'Battery only (essential loads)', y: pct(rC) }] },
        { type: 'bar', title: 'Electrical load analysis by flight phase', ylabel: 'Demand [kW]', categories: PHASES, stacked: true, series: [{ name: 'Propulsion', y: ph.map((p) => p.prop / 1e3) }, { name: 'Essential (priority 1)', y: ph.map((p) => N.sum(p.loads.filter((l) => l.prio === 1 && !l.prop).map((l) => l.P)) / 1e3) }, { name: 'Important (priority 2)', y: ph.map((p) => N.sum(p.loads.filter((l) => l.prio === 2).map((l) => l.P)) / 1e3) }, { name: 'Sheddable (priority 3)', y: ph.map((p) => N.sum(p.loads.filter((l) => l.prio === 3).map((l) => l.P)) / 1e3) }] },
      ],
      tables: [
        { title: 'Electrical load analysis', columns: ['Phase', 'Propulsion [kW]', 'Systems [kW]', 'Total [kW]', 'Capacity [kW]', 'Margin [%]', 'Min bus voltage [V]', 'Cable loss [W]'], rows: ph.map((p, k) => [PHASES[k], p.prop / 1e3, (p.tot - p.prop) / 1e3, p.tot / 1e3, cap / 1e3, 100 * (1 - p.tot / cap), N.amin(p.r.V), p.r.cable]) },
        { title: `Protection sizing at peak demand (${PHASES[kPk]})`, columns: ['Load', 'Bus', 'Priority', 'Power [kW]', 'Current [A]', 'Protection rating [A]', 'Fault / rating [-]'], rows },
        { title: 'Loads shed after loss of one source', columns: ['Order', 'Load'], rows: sB.dropped.length ? sB.dropped.map((nm, k) => [k + 1, nm]) : [[0, 'none required']] },
      ],
      warnings, models: ['Nodal analysis (Kirchhoff current law) with Newton iteration for constant-power loads', 'Thevenin sources with droop for load flow and a fault-current multiple for short circuits', 'Priority-based load shedding', 'Electrical load analysis with generic demand factors per flight phase'],
      assumptions: ['AC systems are treated as a DC-equivalent real-power network: reactive power, harmonics and unbalance are not modelled', 'Split-bus operation with two sources; the bus tie closes when a source is lost', 'Conductors sized at one design current density; contact and protective-device resistances neglected', 'The emergency battery is shown at bus voltage (any DC/DC conversion is ideal)', 'Protection ratings are the next standard size above 125 % of load current; time–current coordination is not analysed', 'Load groups, demand factors, source droop and the default source, system-load and emergency-battery ratings are generic class-level estimates, not an aircraft load analysis'],
    };
  },
  verify() {
    const Vs = 100, R = 0.2, P = 2000, lf = loadFlow([[1 / 1e-9 + 1 / R, -1 / R], [-1 / R, 1 / R]], [Vs / 1e-9, 0], [0, P], Vs);
    const o = N.kv(network.run({ ...Object.fromEntries(network.inputs.map((f) => [f.key, f.default])), P_prop_peak_W: 60000, P_prop_cruise_W: 30000 }));
    return [
      N.check('Constant-power load on one feeder', lf.V[1], (Vs + Math.sqrt(Vs * Vs - 4 * R * P)) / 2, 1e-9, 'V² − Vs·V + R·P = 0'),
      N.check('Kirchhoff current balance at every node', o.kcl_residual_A / (o.P_demand_peak_W / 115), 0, 1e-8, 'ΣI = 0'),
      N.check('Source power = load power + cable loss', o.power_balance_err, 0, 1e-8, 'Tellegen / energy conservation'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.gen_margin_pct < 10) out.push({ severity: o.gen_margin_pct < 0 ? 'critical' : 'warn', title: 'Thin source capacity margin', detail: `${o.gen_margin_pct.toFixed(1)}% at peak demand of ${(o.P_demand_peak_W / 1e3).toPrecision(3)} kW.`, action: 'Increase source rating, stagger intermittent loads, or lower demand factors with load management.', basis: 'Electrical load analysis: capacity ≥ demand with growth margin' });
    if (o.v_drop_feeder_pct > 5) out.push({ severity: 'warn', title: 'Feeder voltage drop exceeds 5 %', detail: `${o.v_drop_feeder_pct.toFixed(1)}% is lost between the source bus and the furthest load (lowest node ${o.bus_V_min.toFixed(1)} V on a ${i.V_bus} V system).`, action: 'Use larger conductors (lower current density), shorter feeders or a higher distribution voltage; see the cable-sizing analysis.', basis: 'Feeder drop limit assumed at 5 % of nominal (typical design practice)' });
    else if (o.v_drop_max_pct > 5) out.push({ severity: 'advise', title: 'Bus voltage sags with the source', detail: `Lowest node ${o.bus_V_min.toFixed(1)} V on a ${i.V_bus} V system: ${(o.v_drop_max_pct - o.v_drop_feeder_pct).toFixed(1)}% is source droop (battery internal resistance or generator regulation) and ${o.v_drop_feeder_pct.toFixed(1)}% feeder drop.`, action: 'Equipment on an unregulated battery bus must accept the whole voltage window of the pack (see the battery analysis); lower-resistance cells or more parallel capacity reduce the sag at peak power.', basis: 'Source droop at rated power (input)' });
    if (o.shed_W > 0) out.push({ severity: 'advise', title: 'Load shedding is required after loss of one source', detail: `${(o.shed_W / 1e3).toPrecision(3)} kW is shed; see the shedding table for the order.`, action: 'Confirm that every shed load is non-essential for continued safe flight and landing, and feed the logic to Suite 22.', basis: 'Priority-based load management' });
    out.push({ severity: 'info', title: 'Electrical losses become heat and fuel or battery energy', detail: `${(o.elec_losses_W / 1e3).toPrecision(3)} kW in cruise.`, action: 'Pass the loss to the thermal suite; higher voltage and efficient conversion cut both the heat load and the energy drawn.', basis: 'Loss summation of this run' });
    return out;
  },
};

// ---- cable sizing and Paschen limit ---------------------------------------------------------
const WIRE = { Copper: { rho: RHO_CU, alpha: ALPHA_CU, dens: 8960 }, Aluminium: { rho: 2.82e-8, alpha: 0.0039, dens: 2700 } };
// Paschen law in Townsend form for air. A = 15 /(Torr·cm) and B = 365 V/(Torr·cm) are the tabulated values (112.5 /(kPa·cm), 2737.5 V/(kPa·cm));
// the secondary-emission coefficient γ = 0.01 is an illustrative value, which puts the minimum at 305 V against the commonly quoted 327 V
const PA = 15, PB = 365, PG = Math.log(1 + 1 / 0.01), PD_MIN = (Math.E * PG) / PA;
const paschen = (pd) => (pd > PD_MIN ? (PB * pd) / Math.log((PA * pd) / PG) : PB * PD_MIN); // conservative: the minimum is used left of it
function cableSize(i, V) {
  const ac = i.system === 'AC three-phase', w = WIRE[i.material] || WIRE.Copper, a = isa(i.alt_m), I = ac ? i.P_W / (3 * V * 0.9) : i.P_W / V, L = ac ? i.L_m : i.L_m * i.loop_factor, nC = ac ? 3 : 1, rho = w.rho * (1 + w.alpha * (i.T_cond_C - 20));
  const derate = Math.max(0.5, a.sigma ** 0.2), A_amp = I / (i.J_max_Amm2 * 1e6 * derate), A_drop = (rho * L * I) / ((i.dv_pct / 100) * V), A = Math.max(A_amp, A_drop);
  const rc = Math.sqrt(A / Math.PI), tIns = (i.ins_base_mm + (i.ins_mm_kV * V * i.v_peak_factor) / 1e3) * 1e-3, R = (rho * L) / A;
  const mCond = nC * w.dens * A * L, mIns = nC * 1700 * Math.PI * ((rc + tIns) ** 2 - rc * rc) * L;
  return { I, A, R, drop: I * R, loss: nC * I * I * R, mass: mCond + mIns, mCond, mIns, gov: A_drop > A_amp ? 1 : 0, derate };
}

const cable = {
  id: 'cable', title: 'Cable sizing versus bus voltage and the Paschen limit', fidelity: 'analytical',
  summary: 'Sizes a feeder for ampacity and voltage drop, shows how cable mass and loss fall as bus voltage rises, and checks the peak voltage against the Paschen breakdown curve of air at altitude.',
  equations: ['Ohm’s law', 'Kirchhoff’s voltage law', 'Power balance equations', 'Paschen’s law (Townsend form)'],
  inputs: [
    { key: 'P_W', label: 'Transmitted power', unit: 'W', default: 150000, min: 1, group: 'Feeder' },
    { key: 'V_bus', label: 'Bus voltage', unit: 'V', default: 800, min: 5, max: 5000, group: 'Feeder' },
    { key: 'system', label: 'System', type: 'select', options: ['DC', 'AC three-phase'], default: 'DC', group: 'Feeder', help: 'For AC the bus voltage is the phase voltage; three phase conductors, power factor 0.9, no neutral' },
    { key: 'L_m', label: 'Route length', unit: 'm', default: 6, min: 0.05, max: 300, group: 'Feeder' },
    { key: 'loop_factor', label: 'Conductor length factor', unit: '-', default: 2, min: 1, max: 2, group: 'Feeder', help: '2 for two-wire, 1 for structure return' },
    { key: 'material', label: 'Conductor material', type: 'select', options: Object.keys(WIRE), default: 'Copper', group: 'Feeder' },
    { key: 'J_max_Amm2', label: 'Allowable current density at sea level', unit: 'A/mm²', default: 6, min: 0.5, max: 20, group: 'Limits', help: 'A bundled-wire rating: about 4–8 A/mm² for mid-gauge copper wire and 2.4–4.6 A/mm² for large feeders (AC 43.13-1B Table 11-9); a single wire in free air carries several times more. Take it from the applicable wiring standard' },
    { key: 'dv_pct', label: 'Allowable voltage drop', unit: '%', default: 2, min: 0.1, max: 15, group: 'Limits' },
    { key: 'T_cond_C', label: 'Conductor temperature', unit: '°C', default: 90, min: -55, max: 260, group: 'Limits' },
    { key: 'alt_m', label: 'Altitude', unit: 'm', default: 3000, min: 0, max: 25000, group: 'Environment' },
    { key: 'ins_base_mm', label: 'Insulation thickness at low voltage', unit: 'mm', default: 0.25, min: 0, max: 5, group: 'Insulation' },
    { key: 'ins_mm_kV', label: 'Added insulation per kV of peak voltage', unit: 'mm/kV', default: 0.5, min: 0, max: 5, group: 'Insulation', help: 'Generic allowance; real values come from partial-discharge-free insulation design' },
    { key: 'v_peak_factor', label: 'Peak / nominal voltage', unit: '-', default: 1, min: 1, max: 3, group: 'Insulation', help: '1 for DC, √2 for AC phase voltage, higher with switching overshoot' },
    { key: 'gap_mm', label: 'Smallest air gap between live parts', unit: 'mm', default: 1, min: 0.01, max: 100, group: 'Insulation' },
  ],
  defaults: (c, up, d) => {
    const e = est(c, up, d), V = c.systems.bus_V || undefined, alt = Math.max(c.mission.cruise_alt_m || 0, c.atm.alt_m || 0);
    return e.electric ? { P_W: (Math.max(e.Ph, e.Pc) / e.n / 0.93) || undefined, V_bus: V, L_m: 0.5 * e.len, alt_m: alt } : { P_W: Math.max(50, (c.systems.gen_kVA || 0) * 900), V_bus: V, L_m: 0.5 * e.len, alt_m: alt, system: V > 60 ? 'AC three-phase' : 'DC', v_peak_factor: V > 60 ? Math.SQRT2 : 1 };
  },
  run(i) {
    const s = cableSize(i, i.V_bus), s2 = cableSize(i, 2 * i.V_bus), Vs = N.logspace(Math.min(28, i.V_bus), Math.max(3000, i.V_bus), 40), sw = Vs.map((v) => cableSize(i, v));
    const a = isa(i.alt_m), torr = (h) => { const q = isa(h); return (q.p / 133.322) * (293.15 / q.T); }, pd = torr(i.alt_m) * i.gap_mm * 0.1, Vbd = paschen(pd), Vpk = i.V_bus * i.v_peak_factor, Vmin = PB * PD_MIN;
    const hs = N.linspace(0, 20000, 41), pds = N.logspace(0.1, 10 ** 2.5, 60), warnings = [];
    if (Vpk > Vmin) warnings.push(`Peak voltage ${Vpk.toFixed(0)} V is above the Paschen minimum of air (${Vmin.toFixed(0)} V in this model; about 327 V is the commonly quoted measured value): partial discharge is possible wherever an unfavourable gap–pressure product exists. Insulation must be designed and tested partial-discharge-free at altitude.`);
    if (Vpk > Vbd) warnings.push(`Peak voltage exceeds the breakdown voltage of the ${i.gap_mm} mm air gap at ${i.alt_m.toFixed(0)} m (${Vbd.toFixed(0)} V).`);
    return {
      kpis: [
        { key: 'cable_area_mm2', label: 'Conductor cross-section', value: s.A * 1e6, unit: 'mm²', note: s.gov ? 'Governed by voltage drop' : 'Governed by ampacity' },
        { key: 'cable_I_A', label: 'Feeder current', value: s.I, unit: 'A' },
        { key: 'cable_mass_kg', label: 'Cable mass (conductor + insulation)', value: s.mass, unit: 'kg' },
        { key: 'cable_loss_W', label: 'Cable loss', value: s.loss, unit: 'W' },
        { key: 'cable_drop_pct', label: 'Voltage drop', value: (100 * s.drop) / i.V_bus, unit: '%', status: (100 * s.drop) / i.V_bus <= i.dv_pct * 1.0001 ? 'ok' : 'warn' },
        { key: 'cable_mass_2V_kg', label: 'Cable mass at twice the voltage', value: s2.mass, unit: 'kg' },
        { key: 'cable_mass_saving_pct', label: 'Mass saved by doubling the voltage', value: 100 * (1 - s2.mass / s.mass), unit: '%' },
        { key: 'ampacity_derate', label: 'Altitude ampacity factor', value: s.derate, unit: '-', note: 'Generic density-based factor' },
        { key: 'paschen_min_V', label: 'Paschen minimum of air (model)', value: Vmin, unit: 'V' },
        { key: 'V_breakdown_gap_V', label: 'Breakdown voltage of the stated gap at altitude', value: Vbd, unit: 'V' },
        { key: 'pd_margin', label: 'Gap breakdown voltage / peak voltage', value: Vbd / Vpk, unit: '-', status: Vpk < Vmin ? 'ok' : Vpk < Vbd ? 'warn' : 'bad', note: 'Below the Paschen minimum no gap can break down' },
      ],
      plots: [
        { type: 'line', title: 'Cable mass versus bus voltage', xlabel: 'Bus voltage [V]', ylabel: 'Mass [kg]', xlog: true, ylog: true, series: [{ name: 'Total', x: Vs, y: sw.map((q) => q.mass) }, { name: 'Conductor', x: Vs, y: sw.map((q) => q.mCond), style: 'dash' }, { name: 'Insulation', x: Vs, y: sw.map((q) => q.mIns), style: 'dash' }], annotations: [{ x: i.V_bus, label: 'Selected' }] },
        { type: 'line', title: 'Cable loss versus bus voltage', xlabel: 'Bus voltage [V]', ylabel: 'Loss [W]', xlog: true, ylog: true, series: [{ name: 'I²R loss', x: Vs, y: sw.map((q) => q.loss) }], annotations: [{ x: i.V_bus, label: 'Selected' }] },
        { type: 'line', title: 'Paschen curve for air', xlabel: 'Pressure × gap [Torr·cm]', ylabel: 'Breakdown voltage [V]', xlog: true, ylog: true, series: [{ name: 'Townsend-form Paschen law', x: pds, y: pds.map((x) => (x > PD_MIN * 0.42 ? (PB * x) / Math.log((PA * x) / PG) : NaN)).map((v) => (v > 0 && v < 1e5 ? v : NaN)) }, { name: 'This gap at altitude', x: [pd], y: [Vbd], style: 'points' }], annotations: [{ y: Vpk, label: 'Peak voltage' }] },
        { type: 'line', title: 'Breakdown voltage of the stated gap versus altitude', xlabel: 'Altitude [m]', ylabel: 'Breakdown voltage [V]', series: [{ name: `${i.gap_mm} mm gap`, x: hs, y: hs.map((h) => paschen(torr(h) * i.gap_mm * 0.1)) }], annotations: [{ y: Vpk, label: 'Peak voltage' }] },
      ],
      warnings, models: ['Conductor sizing for ampacity and voltage drop', 'Paschen law in Townsend form with the tabulated air constants A = 15 /(Torr·cm), B = 365 V/(Torr·cm) and an illustrative γ = 0.01'],
      assumptions: ['DC, or balanced three-phase AC at 0.9 power factor with resistive drop only; skin effect and inductive drop neglected', 'Ampacity as a current-density limit scaled by (density ratio)^0.2 — a generic stand-in for the altitude and bundle derating curves of the applicable wiring standard (0.88 at 6100 m where the AC 43.13-1B worked example uses 0.91)', 'Current-density, voltage-drop and insulation-thickness defaults are typical values', 'Insulation thickness grows linearly with peak voltage (generic allowance); connectors, shielding and supports are not included', 'Paschen: uniform field between clean electrodes; left of the minimum the minimum is used because longer discharge paths usually exist', `Ambient pressure ${(a.p / 1e3).toFixed(1)} kPa from the standard atmosphere`],
    };
  },
  verify() {
    const i = { ...Object.fromEntries(cable.inputs.map((f) => [f.key, f.default])), J_max_Amm2: 1e3, ins_base_mm: 0, ins_mm_kV: 0 }, s = cableSize(i, 400), s2 = cableSize(i, 800);
    const pm = N.goldenSection((x) => (PB * x) / Math.log((PA * x) / PG), PD_MIN * 0.5, PD_MIN * 4, 1e-12);
    return [
      N.check('Drop-limited conductor meets the allowed drop exactly', s.drop / 400, 0.02, 1e-12, 'ΔV = I·ρ·L/A'),
      N.check('Drop-limited cable mass scales with 1/V²', s.mass / s2.mass, 4, 1e-12, 'A ∝ I/V ∝ 1/V²'),
      N.check('Paschen minimum position', pm, (Math.E * PG) / PA, 1e-6, '(p·d)min = e·ln(1 + 1/γ)/A'),
      N.check('Paschen minimum voltage', (PB * pm) / Math.log((PA * pm) / PG), (Math.E * PG * PB) / PA, 1e-9, 'Vmin = e·B·ln(1 + 1/γ)/A'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.cable_mass_saving_pct > 30 && i.V_bus * i.v_peak_factor < 250) out.push({ severity: 'advise', title: 'A higher bus voltage would cut cable mass substantially', detail: `Doubling the voltage saves ${o.cable_mass_saving_pct.toFixed(0)}% of the ${o.cable_mass_kg.toPrecision(3)} kg feeder and reduces I²R loss.`, action: 'Trade a higher distribution voltage against insulation, connector and protection requirements; lighter cables save energy on every flight.', basis: 'I = P/V; drop-limited area ∝ 1/V²' });
    if (o.pd_margin < 1 || i.V_bus * i.v_peak_factor > o.paschen_min_V) out.push({ severity: o.pd_margin < 1 ? 'critical' : 'warn', title: 'Partial-discharge risk at altitude', detail: `Peak voltage ${(i.V_bus * i.v_peak_factor).toFixed(0)} V against a Paschen minimum near ${o.paschen_min_V.toFixed(0)} V.`, action: 'Specify partial-discharge-free insulation and connector systems, control clearances and creepage, and qualify by partial-discharge inception testing at minimum operating pressure.', basis: 'Paschen law for air; qualitative screening only' });
    if (!out.length) out.push({ severity: 'info', title: 'Feeder sized', detail: `${o.cable_area_mm2.toPrecision(3)} mm² ${i.material.toLowerCase()}, ${o.cable_mass_kg.toPrecision(3)} kg, ${o.cable_loss_W.toPrecision(3)} W loss.`, action: 'Aluminium conductors save mass on long, large feeders at the cost of larger terminations.', basis: 'Ampacity and voltage-drop limits set in the inputs' });
    return out;
  },
};

// ---- hybrid-electric energy management ------------------------------------------------------
const BIG = 1e30;
/** Mission as uniform-ish time steps of junction demand. Returns {dt[], P[] (shaft demand), label[]}. */
function hybMission(i) {
  const pk = i.P_peak_W, segs = [[i.land_peak ? i.t_peak_min / 2 : i.t_peak_min, pk], [i.t_climb_min, i.climb_frac * pk], [i.t_cruise_min, i.P_cruise_W], [i.t_desc_min, i.desc_frac * i.P_cruise_W]];
  if (i.land_peak) segs.push([i.t_peak_min / 2, pk]);
  const tt = N.sum(segs.map((s) => s[0])) || 1, nT = Math.round(i.nT), dt = [], P = [];
  for (const [t, p] of segs) { if (!(t > 0)) continue; const n = Math.max(1, Math.round((nT * t) / tt)); for (let k = 0; k < n; k++) { dt.push((t * 60) / n); P.push(p); } }
  return { dt, P };
}
/** Fuel power [W] of the engine at shaft power P (Willans line: idle fuel fraction f of rated fuel flow). */
const fuelPower = (P, Pr, eta, f, off) => (P <= 1e-9 * Pr ? (off ? 0 : (f * Pr) / eta) : (f * Pr + (1 - f) * P) / eta);
function hybSolve(i) {
  const ms = hybMission(i), nT = ms.dt.length, series = i.arch === 'series', fl = FLUIDS[i.fuel] || FLUIDS['Jet A-1'], E = i.batt_kWh * 3.6e6;
  const eb = series ? 1 : i.eta_drive, ee = series ? i.eta_gen : i.eta_gbx, D = ms.P.map((p) => (series ? p / i.eta_drive : p)); // junction: DC bus (series) or shaft (parallel)
  const Pbase = i.P_peak_W / i.eta_gbx, Pr = i.eng_frac * Pbase, Pb = i.c_rate_max * i.batt_kWh * 1e3, sHi = Math.max(i.soc0, 0.9), fp = (P, R = Pr) => fuelPower(P, R, i.eta_eng, i.f_idle, i.allow_off);
  // one step: terminal battery power u (+ discharge) -> {soc', fuel energy} or null when infeasible
  const step = (k, soc, u) => {
    if (Math.abs(u) > Pb * (1 + 1e-9)) return null;
    const s2 = soc - ((u > 0 ? u / i.eta_batt : u * i.eta_batt) * ms.dt[k]) / E, e = D[k] - (u > 0 ? u * eb : u / eb);
    if (s2 < i.soc_min - 1e-9 || s2 > sHi + 1e-9 || e < -1e-9 * (D[k] + 1)) return null;
    const Pe = Math.max(0, e) / ee; if (Pe > Pr * (1 + 1e-9)) return null;
    return { s2: N.clamp(s2, i.soc_min, sHi), fuel: fp(Pe) * ms.dt[k], Pe };
  };
  const uDef = (k) => Math.max(0, (D[k] - Pr * ee) / eb); // boost the engine cannot supply
  // all-engine baseline
  let fConv = 0; for (let k = 0; k < nT; k++) fConv += fp(ms.P[k] / i.eta_gbx, Pbase) * ms.dt[k];
  // rule-based: mandatory boost, plus a fixed blend while the charge exceeds what later boosts will need
  const need = new Array(nT + 1).fill(0); for (let k = nT - 1; k >= 0; k--) need[k] = need[k + 1] + (uDef(k) * ms.dt[k]) / i.eta_batt / E;
  const rule = { soc: [i.soc0], u: [], Pe: [], fuel: 0, ok: true }; let s = i.soc0;
  for (let k = 0; k < nT; k++) {
    const uAvail = Math.min(Pb, ((s - i.soc_min) * E * i.eta_batt) / ms.dt[k]), blend = s - i.soc_min > need[k + 1] + (i.hyb_frac * D[k] / eb) * ms.dt[k] / i.eta_batt / E ? (i.hyb_frac * D[k]) / eb : 0;
    let u = Math.min(Math.max(uDef(k), blend), Math.max(0, uAvail)), st = step(k, s, u);
    if (!st) { rule.ok = false; u = Math.min(uDef(k), Math.max(0, uAvail)); st = { s2: Math.max(i.soc_min, s - ((u / i.eta_batt) * ms.dt[k]) / E), fuel: fp(Pr) * ms.dt[k], Pe: Pr }; }
    s = st.s2; rule.fuel += st.fuel; rule.soc.push(s); rule.u.push(u); rule.Pe.push(st.Pe);
  }
  // dynamic programming over SOC with interpolated cost-to-go. The grid of step k spans [lo[k], sHi], where lo[k] is the
  // exact lowest SOC from which the rest of the mission is still feasible (mandatory boosts less the charging the engine can do).
  const nS = Math.max(2, Math.round(i.nSoc)), lo = new Array(nT + 1).fill(i.soc_min), uGrid = N.linspace(-Pb, Pb, 41); let feas = true;
  for (let k = nT - 1; k >= 0; k--) {
    const ud = uDef(k), chg = ud > 0 ? 0 : Math.min(Pb, Math.max(0, Pr * ee - D[k]) * eb) * i.eta_batt;
    if (ud > Pb * (1 + 1e-9)) feas = false;
    lo[k] = Math.max(i.soc_min, lo[k + 1] + ((ud / i.eta_batt - chg) * ms.dt[k]) / E);
  }
  if (lo[0] > i.soc0 + 1e-12 || N.amax(lo) > sHi) feas = false;
  const grids = lo.map((l) => N.linspace(Math.min(l, sHi), sHi, nS));
  const interp = (k, J, sv) => { const l = lo[k], h = (sHi - l) / (nS - 1); if (sv < l - 1e-9) return BIG; if (!(h > 0)) return J[0]; const x = N.clamp((sv - l) / h, 0, nS - 1), j = Math.min(nS - 2, Math.floor(x)), w = x - j; return (1 - w) * J[j] + w * J[j + 1]; };
  const toSoc = (k, soc, target) => { const dE = ((soc - target) * E) / ms.dt[k]; return dE >= 0 ? dE * i.eta_batt : dE / i.eta_batt; }; // terminal power that lands on a target SOC
  const best = (k, soc, Jn) => {
    let bu = NaN, bc = BIG, bs = null;
    for (const u of uGrid.concat([0, uDef(k), D[k] / eb, toSoc(k, soc, lo[k + 1]), toSoc(k, soc, sHi)])) { const st = step(k, soc, u); if (!st) continue; const cst = st.fuel + interp(k + 1, Jn, st.s2); if (cst < bc) { bc = cst; bu = u; bs = st; } }
    return { u: bu, cost: bc, st: bs };
  };
  const dp = { soc: [i.soc0], u: [], Pe: [], fuel: 0, ok: feas };
  if (feas) {
    const Js = new Array(nT + 1); Js[nT] = new Array(nS).fill(0);
    for (let k = nT - 1; k >= 0; k--) Js[k] = grids[k].map((sv) => best(k, sv, Js[k + 1]).cost);
    s = i.soc0;
    for (let k = 0; k < nT; k++) { const bq = best(k, s, Js[k + 1]); if (!bq.st || bq.cost >= BIG) { dp.ok = false; break; } s = bq.st.s2; dp.fuel += bq.st.fuel; dp.soc.push(s); dp.u.push(bq.u); dp.Pe.push(bq.st.Pe); }
  }
  const t = [0]; for (const dtk of ms.dt) t.push(t[t.length - 1] + dtk / 60);
  return { ms, D, t, fConv: fConv / fl.LHV, rule, dp, fl, Pr, Pbase, E, Eshaft: N.sum(ms.P.map((p, k) => p * ms.dt[k])) };
}

const hybrid = {
  id: 'hybrid', title: 'Hybrid-electric energy management: rule-based versus optimal', fidelity: 'numerical',
  summary: 'Splits mission power between a (possibly downsized) engine and a battery in a series or parallel architecture, comparing a rule-based strategy with the fuel-optimal split found by dynamic programming over the state of charge.',
  equations: ['Hybrid-electric propulsion power-split equations', 'Power balance equations', 'Gas turbine–electric motor power balance', 'Bellman dynamic-programming recursion'],
  inputs: [
    { key: 'arch', label: 'Architecture', type: 'select', options: ['parallel', 'series'], default: 'parallel', group: 'Architecture', help: 'Parallel: engine and motor both drive the shaft. Series: engine drives a generator only' },
    { key: 'P_peak_W', label: 'Peak shaft power demand', unit: 'W', default: 2e6, min: 10, group: 'Mission' },
    { key: 'P_cruise_W', label: 'Cruise shaft power demand', unit: 'W', default: 1e6, min: 1, group: 'Mission' },
    { key: 't_peak_min', label: 'Time at peak power (take-off / hover)', unit: 'min', default: 2, min: 0, group: 'Mission' },
    { key: 'land_peak', label: 'Split peak time between take-off and landing', type: 'bool', default: false, group: 'Mission', help: 'On for VTOL: half the peak time before and half after the cruise' },
    { key: 't_climb_min', label: 'Climb time', unit: 'min', default: 15, min: 0, group: 'Mission' },
    { key: 'climb_frac', label: 'Climb power / peak power', unit: '-', default: 0.8, min: 0.1, max: 1, group: 'Mission' },
    { key: 't_cruise_min', label: 'Cruise time', unit: 'min', default: 90, min: 0, group: 'Mission' },
    { key: 't_desc_min', label: 'Descent time', unit: 'min', default: 15, min: 0, group: 'Mission' },
    { key: 'desc_frac', label: 'Descent power / cruise power', unit: '-', default: 0.3, min: 0, max: 1, group: 'Mission' },
    { key: 'eng_frac', label: 'Engine rating / all-engine baseline rating', unit: '-', default: 0.8, min: 0.05, max: 1, group: 'Engine', help: 'Below 1 the engine is downsized and the battery must supply the peak' },
    { key: 'eta_eng', label: 'Engine thermal efficiency at rated power', unit: '-', default: 0.3, min: 0.1, max: 0.6, group: 'Engine', help: 'From the propulsion suite when available' },
    { key: 'f_idle', label: 'Idle fuel flow / rated fuel flow', unit: '-', default: 0.2, min: 0, max: 0.6, group: 'Engine', help: 'Willans-line part-load model; ≈ 0.2–0.3 for gas turbines, ≈ 0.1 for piston engines' },
    { key: 'allow_off', label: 'Engine may be shut down in flight', type: 'bool', default: false, group: 'Engine' },
    { key: 'fuel', label: 'Fuel', type: 'select', options: FUELS, default: 'Jet A-1', group: 'Engine' },
    { key: 'batt_kWh', label: 'Battery energy', unit: 'kWh', default: 150, min: 0.001, group: 'Battery' },
    { key: 'soc0', label: 'Initial state of charge', unit: '-', default: 0.95, min: 0.1, max: 1, group: 'Battery' },
    { key: 'soc_min', label: 'Minimum state of charge', unit: '-', default: 0.2, min: 0, max: 0.9, group: 'Battery' },
    { key: 'c_rate_max', label: 'Maximum charge / discharge rate', unit: 'C', default: 4, min: 0.1, max: 30, group: 'Battery' },
    { key: 'eta_batt', label: 'Battery one-way efficiency', unit: '-', default: 0.96, min: 0.5, max: 1, group: 'Battery' },
    { key: 'wh_kg_pack', label: 'Pack specific energy', unit: 'Wh/kg', default: 180, min: 30, max: 600, group: 'Battery' },
    { key: 'eta_drive', label: 'Inverter × motor efficiency', unit: '-', default: 0.93, min: 0.5, max: 1, group: 'Efficiencies' },
    { key: 'eta_gen', label: 'Generator × rectifier efficiency', unit: '-', default: 0.94, min: 0.5, max: 1, group: 'Efficiencies' },
    { key: 'eta_gbx', label: 'Gearbox efficiency', unit: '-', default: 0.98, min: 0.5, max: 1, group: 'Efficiencies' },
    { key: 'hyb_frac', label: 'Rule-based battery share of demand', unit: '-', default: 0.1, min: 0, max: 1, group: 'Strategy', help: 'Blend used while charge exceeds what later peak boosts will need' },
    { key: 'nT', label: 'Time steps', unit: '', default: 80, min: 10, max: 600, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nSoc', label: 'State-of-charge grid points', unit: '', default: 61, min: 5, max: 401, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const e = est(c, up, d), vt = e.Ph > 0, ef = e.electric ? 0.4 : 0.9, pk = e.Ppeak || undefined, alt = c.mission.cruise_alt_m || 1000, t = c.prop.type;
    const etaUp = up.propulsion?.eta_thermal, eta = etaUp > 0.1 && etaUp < 0.6 ? etaUp : t === 'turbofan' || t === 'turbojet' ? 0.42 : 0.3;
    const boost = (1 - ef) * (pk || 0) / 0.98 / 0.93 / 1e3; // kW the battery must add at peak power
    return { P_peak_W: pk, P_cruise_W: Math.min(e.Pc, pk || e.Pc) || undefined, t_peak_min: vt ? Math.max(1, c.mission.hover_min || 2) : 2, land_peak: vt, t_climb_min: Math.max(1, alt / 480), t_cruise_min: e.tCruise, t_desc_min: Math.max(1, alt / 300),
      eng_frac: ef, eta_eng: eta, f_idle: t === 'piston' ? 0.1 : 0.2, fuel: FUELS.includes(c.prop.fuel) ? c.prop.fuel : undefined, batt_kWh: e.electric ? c.systems.batt_kWh || undefined : Math.max(1.4 * boost / 4, 0.002 * c.mass.mtow_kg) || undefined };
  },
  run(i) {
    const r = hybSolve(i), fl = r.fl, warnings = [], fR = r.rule.fuel / fl.LHV, fD = r.dp.ok ? r.dp.fuel / fl.LHV : NaN, mB = (i.batt_kWh * 1e3) / i.wh_kg_pack;
    if (!r.rule.ok) warnings.push('Rule-based strategy: the battery cannot supply the boost the downsized engine needs in at least one step (energy or C-rate limit). Its fuel figure assumes the engine at rated power there and the demand is not met.');
    if (!r.dp.ok) warnings.push('No feasible power split exists for this mission: enlarge the battery, raise its C-rate limit or the engine rating.');
    if (r.dp.ok && r.rule.ok && fD > fR * 1.002) warnings.push('The dynamic-programming result is slightly worse than the rule-based one: refine the SOC grid and time steps.');
    if (i.allow_off) warnings.push('In-flight engine shutdown is permitted in this run: restart reliability and one-engine-inoperative performance must be justified separately.');
    const tm = r.t.slice(0, -1).map((v, k) => 0.5 * (v + r.t[k + 1])), kW = (a) => a.map((v) => v / 1e3), save = 100 * (1 - fD / r.fConv), usedD = r.dp.ok ? (i.soc0 - r.dp.soc[r.dp.soc.length - 1]) * i.batt_kWh : NaN;
    return {
      kpis: [
        { key: 'fuel_conv_kg', label: 'Fuel, all-engine baseline', value: r.fConv, unit: 'kg' },
        { key: 'fuel_rule_kg', label: 'Fuel, hybrid with rule-based split', value: fR, unit: 'kg', status: r.rule.ok ? 'ok' : 'bad' },
        { key: 'fuel_dp_kg', label: 'Fuel, hybrid with optimal split', value: fD, unit: 'kg', status: r.dp.ok ? 'ok' : 'bad' },
        { key: 'fuel_saved_pct', label: 'Fuel saved versus baseline (optimal)', value: save, unit: '%', status: save > 0 ? 'ok' : 'warn', note: 'Battery mass penalty on power demand is not included' },
        { key: 'dp_gain_pct', label: 'Optimal versus rule-based fuel', value: 100 * (1 - fD / fR), unit: '%' },
        { key: 'co2_saved_kg', label: 'Tailpipe CO₂ avoided per mission', value: (r.fConv - fD) * (fl.co2_per_kg || 0), unit: 'kg' },
        { key: 'hyb_soc_end', label: 'Final state of charge (optimal)', value: r.dp.ok ? r.dp.soc[r.dp.soc.length - 1] : NaN, unit: '-' },
        { key: 'hyb_batt_used_kWh', label: 'Battery energy used (optimal)', value: usedD, unit: 'kWh' },
        { key: 'hyb_energy_frac', label: 'Battery share of mission shaft energy', value: (usedD * 3.6e6 * i.eta_batt * i.eta_drive) / r.Eshaft, unit: '-' },
        { key: 'hyb_batt_mass_kg', label: 'Battery pack mass', value: mB, unit: 'kg' },
        { key: 'eng_rated_W', label: 'Hybrid engine rating', value: r.Pr, unit: 'W' },
        { key: 'mission_shaft_kWh', label: 'Mission shaft energy', value: r.Eshaft / 3.6e6, unit: 'kWh' },
      ],
      plots: [
        { type: 'line', title: 'State-of-charge trajectories', xlabel: 'Time [min]', ylabel: 'SOC [-]', series: [{ name: 'Rule-based', x: r.t.slice(0, r.rule.soc.length), y: r.rule.soc }, { name: 'Optimal (dynamic programming)', x: r.t.slice(0, r.dp.soc.length), y: r.dp.soc }], annotations: [{ y: i.soc_min, label: 'Minimum' }] },
        { type: 'line', title: 'Optimal power split', xlabel: 'Time [min]', ylabel: 'Power [kW]', series: [{ name: 'Shaft demand', x: tm, y: kW(r.ms.P), style: 'step' }, { name: 'Engine shaft power', x: tm.slice(0, r.dp.Pe.length), y: kW(r.dp.Pe), style: 'step' }, { name: 'Battery terminal power', x: tm.slice(0, r.dp.u.length), y: kW(r.dp.u), style: 'step' }], annotations: [{ y: r.Pr / 1e3, label: 'Engine rating' }] },
        { type: 'line', title: 'Rule-based power split', xlabel: 'Time [min]', ylabel: 'Power [kW]', series: [{ name: 'Shaft demand', x: tm, y: kW(r.ms.P), style: 'step' }, { name: 'Engine shaft power', x: tm, y: kW(r.rule.Pe), style: 'step' }, { name: 'Battery terminal power', x: tm, y: kW(r.rule.u), style: 'step' }] },
        { type: 'bar', title: 'Mission fuel by strategy', ylabel: 'Fuel [kg]', categories: ['All-engine baseline', 'Hybrid, rule-based', 'Hybrid, optimal'], series: [{ name: 'Fuel', y: [r.fConv, fR, r.dp.ok ? fD : 0] }] },
      ],
      warnings, models: ['Quasi-static power-split model (series or parallel) with constant conversion efficiencies', 'Willans-line engine part-load fuel model', 'Rule-based blend with boost reserve', 'Deterministic dynamic programming over state of charge with interpolated cost-to-go'],
      assumptions: ['Power demand is prescribed: the battery and electrical-machine mass is not fed back into drag or power (do that in Suite 23/24)', 'Rated thermal efficiency is the same for the baseline and the downsized engine (no scale effect), with no altitude lapse of engine power', 'Battery efficiency is constant; voltage sag and thermal limits are in the battery analysis', 'The battery is charged on the ground: fuel saving is tailpipe only and excludes electricity generation', 'Conversion efficiencies, the idle fuel fraction, C-rate limit and pack specific energy are typical values, not data for specific equipment'],
    };
  },
  convergence: { param: 'nSoc', label: 'State-of-charge grid points', levels: [21, 41, 81, 161], metric: 'fuel_dp_kg' },
  verify() {
    // Constant engine efficiency and loss-free conversion: fuel = (demand energy − usable battery energy)/(η·LHV), whatever the split.
    const i = { ...Object.fromEntries(hybrid.inputs.map((f) => [f.key, f.default])), f_idle: 0, eta_batt: 1, eta_drive: 1, eta_gen: 1, eta_gbx: 1, eng_frac: 1, c_rate_max: 20 }, r = hybSolve(i), o = N.kv(hybrid.run(i));
    const ex = (r.Eshaft - (i.soc0 - i.soc_min) * r.E) / (i.eta_eng * r.fl.LHV), d = Object.fromEntries(hybrid.inputs.map((f) => [f.key, f.default])), od = N.kv(hybrid.run(d));
    return [
      N.check('Optimal fuel with constant efficiency equals the energy balance', o.fuel_dp_kg, ex, 1e-9, 'First law: fuel energy·η = demand − battery energy'),
      N.check('All-engine baseline fuel', o.fuel_conv_kg, r.Eshaft / (i.eta_eng * r.fl.LHV), 1e-12, 'Fuel = ∫P dt/(η·LHV)'),
      N.check('Optimal split never burns more than the rule-based split', Math.max(0, od.fuel_dp_kg / od.fuel_rule_kg - 1), 0, 2e-3, 'Optimality of dynamic programming (grid tolerance)'),
    ];
  },
  recommend(res, i, ctx) {
    const o = res.outputs, out = [], mtow = ctx?.case?.mass?.mtow_kg;
    if (!Number.isFinite(o.fuel_dp_kg)) return [{ severity: 'critical', title: 'The hybrid system cannot fly this mission', detail: 'No feasible power split was found.', action: 'Increase battery energy or its C-rate limit, or raise the engine rating fraction.', basis: 'Power and energy feasibility' }];
    out.push({ severity: o.fuel_saved_pct > 3 ? 'advise' : 'info', title: `Hybridisation changes mission fuel by ${(-o.fuel_saved_pct).toFixed(1)}%`, detail: `${o.fuel_conv_kg.toPrecision(4)} kg baseline, ${o.fuel_dp_kg.toPrecision(4)} kg with the optimal split (${o.co2_saved_kg.toPrecision(3)} kg tailpipe CO₂ avoided). The battery supplies ${(100 * o.hyb_energy_frac).toFixed(1)}% of shaft energy.`, action: 'The gain comes from running a smaller engine nearer its best efficiency and from stored grid energy; it shrinks on long missions. Combine with SAF or hydrogen (fuel selector) for larger CO₂ cuts and check life-cycle electricity emissions.', basis: 'Mission fuel integration of this run' });
    const mAdd = o.hyb_batt_mass_kg - ((ctx?.case?.systems?.batt_kWh || 0) * 1e3) / i.wh_kg_pack; // battery beyond what the aircraft already carries
    if (mtow && mAdd > 0.05 * mtow) out.push({ severity: 'warn', title: 'Added battery mass is a large share of take-off mass', detail: `${mAdd.toFixed(0)} kg more battery than the aircraft carries now is ${(100 * mAdd / mtow).toFixed(1)}% of MTOM and is not yet reflected in the power demand.`, action: 'Re-size the aircraft with the battery mass in Suite 23 and re-fly the mission in Suite 24 before trusting the fuel saving.', basis: 'Mass–energy coupling (Breguet)' });
    if (o.dp_gain_pct > 1) out.push({ severity: 'advise', title: 'The energy-management strategy matters', detail: `The optimal split uses ${o.dp_gain_pct.toFixed(1)}% less fuel than the simple rule.`, action: 'Implement a predictive or equivalent-consumption strategy that follows the optimal state-of-charge trajectory shown.', basis: 'Dynamic-programming benchmark' });
    return out;
  },
};

// ---- DC/DC converter and bus transient ------------------------------------------------------
/** Averaged steady state of a buck (Vout < Vin) or boost converter in continuous conduction with conduction and switching losses. */
function dcdc(i) {
  const buck = i.Vout <= i.Vin, f = i.f_sw_kHz * 1e3, L = i.L_uH * 1e-6, Cc = i.C_uF * 1e-6, RL = i.R_L_mohm * 1e-3, Ron = i.R_on_mohm * 1e-3, Io = i.P_out_W / i.Vout;
  let D, IL = Io;
  if (buck) D = (i.Vout + IL * (RL + Ron) + i.V_d) / (i.Vin + i.V_d);
  else { D = 1 - i.Vin / i.Vout; for (let k = 0; k < 80; k++) { IL = Io / Math.max(1 - D, 1e-3); D = N.clamp((i.Vout + i.V_d + IL * (Ron + RL) - i.Vin) / (i.Vout + i.V_d), 0, 0.999); } IL = Io / Math.max(1 - D, 1e-3); }
  const dI = buck ? (i.Vout * (1 - D)) / (L * f) : (i.Vin * D) / (L * f), dV = buck ? dI / (8 * f * Cc) : (Io * D) / (f * Cc);
  const Pcond = (IL * IL + (dI * dI) / 12) * (RL + Ron) + (1 - D) * i.V_d * IL, Psw = 0.5 * (buck ? i.Vin : i.Vout) * IL * i.t_sw_ns * 1e-9 * f; // ½·V·I per linear transition, turn-on plus turn-off each cycle
  return { buck, D, IL, dI, dV, Pcond, Psw, eff: i.P_out_W / (i.P_out_W + Pcond + Psw), ccm: dI / 2 < IL, valid: D > 0 && D < 0.98 };
}
/** Source–line–capacitor bus: L di/dt = Vs − R·i − v, C dv/dt = i − i_load(v). Load step P1 → P2 at t = 0 (resistive or constant power). */
function busStep(i, n) {
  const Vs = i.Vin, R = i.R_line_mohm * 1e-3, L = i.L_line_uH * 1e-6, Cb = i.C_bus_uF * 1e-6, cpl = i.load_type === 'constant power', G1 = i.P_step1_W / Vs ** 2, G2 = i.P_step2_W / Vs ** 2;
  const v0 = cpl ? (Vs + Math.sqrt(Math.max(0, Vs * Vs - 4 * R * i.P_step1_W))) / 2 : Vs / (1 + R * G1), i0 = cpl ? i.P_step1_W / v0 : G1 * v0;
  const load = (v) => (cpl ? i.P_step2_W / Math.max(v, 0.2 * Vs) : G2 * v), r = N.rk4((t, y) => [(Vs - R * y[0] - y[1]) / L, (y[0] - load(y[1])) / Cb], 0, [i0, v0], i.t_end_ms * 1e-3, n);
  const vss = cpl ? (Vs + Math.sqrt(Math.max(0, Vs * Vs - 4 * R * i.P_step2_W))) / 2 : Vs / (1 + R * G2), Geff = cpl ? -i.P_step2_W / (vss * vss) : G2;
  const wn = Math.sqrt(Math.max(1e-30, (1 + R * Geff) / (L * Cb))), zeta = (R / L + Geff / Cb) / (2 * wn);
  return { t: r.t, v: r.y.map((y) => y[1]), i: r.y.map((y) => y[0]), v0, i0, vss, wn, zeta, G2, R, L, Cb, Vs };
}

const converter = {
  id: 'converter', title: 'DC/DC converter design point and bus load-step transient', fidelity: 'numerical',
  summary: 'Average-value buck or boost converter: duty cycle, inductor and capacitor ripple and efficiency; then the response of the supply bus (line resistance and inductance with a bus capacitor) to a load step, including the constant-power-load stability margin.',
  equations: ['Kirchhoff’s voltage law', 'Kirchhoff’s current law', 'Converter switching equations (state-space averaged)', 'Power network differential–algebraic equations (RLC bus)'],
  inputs: [
    { key: 'Vin', label: 'Input (bus) voltage', unit: 'V', default: 270, min: 3, max: 3000, group: 'Converter' },
    { key: 'Vout', label: 'Output voltage', unit: 'V', default: 28, min: 1, max: 3000, group: 'Converter', help: 'Below the input → buck; above → boost' },
    { key: 'P_out_W', label: 'Output power', unit: 'W', default: 1000, min: 0.1, group: 'Converter' },
    { key: 'f_sw_kHz', label: 'Switching frequency', unit: 'kHz', default: 100, min: 1, max: 2000, group: 'Converter' },
    { key: 'L_uH', label: 'Inductance', unit: 'µH', default: 22, min: 0.01, group: 'Converter' },
    { key: 'C_uF', label: 'Output capacitance', unit: 'µF', default: 47, min: 0.01, group: 'Converter' },
    { key: 'R_L_mohm', label: 'Inductor resistance', unit: 'mΩ', default: 5, min: 0, group: 'Losses' },
    { key: 'R_on_mohm', label: 'Switch on-resistance', unit: 'mΩ', default: 10, min: 0, group: 'Losses' },
    { key: 'V_d', label: 'Rectifier forward voltage', unit: 'V', default: 0, min: 0, max: 2, group: 'Losses', help: '0 for a synchronous rectifier (its resistance is then taken equal to the switch)' },
    { key: 't_sw_ns', label: 'Switching transition time (rise + fall)', unit: 'ns', default: 40, min: 0, max: 2000, group: 'Losses' },
    { key: 'R_line_mohm', label: 'Source and line resistance', unit: 'mΩ', default: 700, min: 0.001, group: 'Bus transient' },
    { key: 'L_line_uH', label: 'Line inductance', unit: 'µH', default: 20, min: 0.001, group: 'Bus transient', help: 'About 1 µH per metre of two-wire run' },
    { key: 'C_bus_uF', label: 'Bus capacitance', unit: 'µF', default: 30, min: 0.001, group: 'Bus transient' },
    { key: 'P_step1_W', label: 'Load before the step', unit: 'W', default: 200, min: 0, group: 'Bus transient' },
    { key: 'P_step2_W', label: 'Load after the step', unit: 'W', default: 1000, min: 0.1, group: 'Bus transient' },
    { key: 'load_type', label: 'Load behaviour', type: 'select', options: ['resistive', 'constant power'], default: 'resistive', group: 'Bus transient', help: 'Regulated converters and drives behave as constant-power loads (negative incremental resistance)' },
    { key: 't_end_ms', label: 'Simulated time', unit: 'ms', default: 1, min: 0.001, max: 5000, group: 'Numerics' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 2000, min: 100, max: 200000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const e = est(c, up, d), bus = c.systems.bus_V || 28, Vin = e.electric ? bus : bus > 60 ? 270 : 28, Vout = Vin > 60 ? 28 : 12;
    const P = N.clamp(e.electric ? 0.02 * Math.max(e.Ph, e.Pc) : 0.05 * e.genTot, 10, 5000), f = 1e5, D = Vout / Vin, Io = P / Vout, dI = 0.3 * Io;
    const Rl = (0.01 * Vin * Vin) / P, Ll = Math.max(1, 1.0 * e.len), Cb = Math.max(1, (2e-3 * P) / (Vin * Vin) * 1e6), wn = 1 / Math.sqrt(Ll * 1e-6 * Cb * 1e-6), zeta = N.clamp(Rl / (Ll * 1e-6) / (2 * wn), 0.02, 50);
    return { Vin, Vout, P_out_W: P, L_uH: ((Vout * (1 - D)) / (dI * f)) * 1e6, C_uF: (dI / (8 * f * 0.01 * Vout)) * 1e6, R_L_mohm: (0.01 * Vout / Io) * 1e3, R_on_mohm: (0.01 * Vout / Io) * 1e3, R_line_mohm: Rl * 1e3, L_line_uH: Ll, C_bus_uF: Cb, P_step1_W: 0.2 * P, P_step2_W: P, t_end_ms: (1e3 * Math.max(8 / (Math.min(zeta, 1 / zeta) * wn), 30 * Rl * Cb * 1e-6)) };
  },
  run(i) {
    const s = dcdc(i), n = Math.round(i.nSteps), b = busStep(i, n), warnings = [], vMin = N.amin(b.v), under = (100 * (b.vss - vMin)) / b.vss;
    const band = 0.02 * Math.abs(b.v0 - b.vss) + 1e-12; let kS = 0; for (let k = b.v.length - 1; k >= 0; k--) if (Math.abs(b.v[k] - b.vss) > band) { kS = Math.min(k + 1, b.v.length - 1); break; }
    const cplMargin = (b.Vs * b.Vs * b.R * b.Cb) / (b.L * i.P_step2_W), settled = Math.abs(b.v[b.v.length - 1] - b.vss) <= band, idx = ds(N.range(b.t.length), 400);
    if (!s.valid) warnings.push('The required duty cycle is outside 0–0.98: this conversion ratio is not practical for a single non-isolated stage.');
    if (s.valid && s.D < 0.08) warnings.push('Very small duty cycle: a transformer-isolated or two-stage topology is normally used for this conversion ratio, and switching loss dominates.');
    if (!s.ccm) warnings.push('Inductor ripple exceeds twice the average current: the converter enters discontinuous conduction and the averaged continuous-conduction relations no longer hold.');
    if (i.load_type === 'constant power' && (cplMargin < 1 || b.zeta < 0)) warnings.push('The bus is unstable with a constant-power load: the negative incremental resistance outweighs the line damping. Add bus capacitance or damping.');
    if (!settled) warnings.push('The transient has not settled within the simulated time: extend it.');
    return {
      kpis: [
        { key: 'duty', label: `Duty cycle (${s.buck ? 'buck' : 'boost'})`, value: s.D, unit: '-', status: s.valid ? 'ok' : 'bad' },
        { key: 'conv_eff', label: 'Converter efficiency', value: s.eff, unit: '-', status: s.eff > 0.9 ? 'ok' : 'warn' },
        { key: 'conv_loss_W', label: 'Converter losses', value: s.Pcond + s.Psw, unit: 'W' },
        { key: 'I_L_A', label: 'Average inductor current', value: s.IL, unit: 'A' },
        { key: 'ripple_I_A', label: 'Inductor current ripple (peak-to-peak)', value: s.dI, unit: 'A', status: s.ccm ? 'ok' : 'warn' },
        { key: 'ripple_V_pct', label: 'Output voltage ripple (peak-to-peak)', value: (100 * s.dV) / i.Vout, unit: '%', status: s.dV / i.Vout < 0.02 ? 'ok' : 'warn', note: 'Typical target ≤ 1–2 %' },
        { key: 'step_V_min_V', label: 'Minimum bus voltage after the load step', value: vMin, unit: 'V' },
        { key: 'step_undershoot_pct', label: 'Undershoot below the new steady state', value: Math.max(0, under), unit: '%' },
        { key: 'step_V_final_V', label: 'New steady-state bus voltage', value: b.vss, unit: 'V' },
        { key: 'bus_zeta', label: 'Bus damping ratio (linearised)', value: b.zeta, unit: '-', status: b.zeta > 0.2 ? 'ok' : b.zeta > 0 ? 'warn' : 'bad' },
        { key: 'bus_fn_Hz', label: 'Bus natural frequency', value: b.wn / (2 * Math.PI), unit: 'Hz' },
        { key: 'step_settle_ms', label: 'Settling time (2 % band)', value: b.t[kS] * 1e3, unit: 'ms' },
        { key: 'cpl_margin', label: 'Constant-power-load stability margin V²RC/(L·P)', value: cplMargin, unit: '-', status: cplMargin > 2 ? 'ok' : cplMargin > 1 ? 'warn' : 'bad', note: 'Must exceed 1 for a constant-power load' },
      ],
      plots: [
        { type: 'line', title: 'Bus voltage after the load step', xlabel: 'Time [ms]', ylabel: 'Voltage [V]', series: [{ name: 'Bus voltage', x: idx.map((k) => b.t[k] * 1e3), y: idx.map((k) => b.v[k]) }], annotations: [{ y: b.vss, label: 'New steady state' }] },
        { type: 'line', title: 'Line current after the load step', xlabel: 'Time [ms]', ylabel: 'Current [A]', series: [{ name: 'Line current', x: idx.map((k) => b.t[k] * 1e3), y: idx.map((k) => b.i[k]) }] },
        { type: 'bar', title: 'Converter loss breakdown', ylabel: 'Loss [W]', categories: ['Conduction', 'Switching'], series: [{ name: 'Loss', y: [s.Pcond, s.Psw] }] },
      ],
      warnings, models: ['State-space-averaged buck/boost converter in continuous conduction', 'Conduction (I²R, forward drop) and linear-transition switching loss ½·V·I·(t_rise + t_fall)·f', 'Second-order RLC bus model with resistive or constant-power load (RK4)'],
      assumptions: ['Ideal regulation: the converter holds its output voltage, so only steady-state duty and ripple are computed', 'Capacitor ESR, magnetic core loss and gate-drive loss are neglected', 'Line modelled as a lumped series R–L; the source is an ideal voltage behind it', 'Switching-level waveforms and control-loop dynamics are not simulated', 'Default component and line values are illustrative, sized for about 30% inductor ripple, 1% output ripple and 1% line drop'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps in the bus transient', levels: [250, 500, 1000, 2000], metric: 'step_V_min_V' },
  verify() {
    const i = { ...Object.fromEntries(converter.inputs.map((f) => [f.key, f.default])), R_line_mohm: 200, t_end_ms: 0.2 }, b = busStep(i, 4000), id = dcdc({ ...i, R_L_mohm: 0, R_on_mohm: 0, V_d: 0, t_sw_ns: 0 });
    // analytic underdamped response of v'' + 2ζωn v' + ωn² v = ωn² vss
    const al = b.zeta * b.wn, wd = b.wn * Math.sqrt(1 - b.zeta ** 2), A = b.v0 - b.vss, dv0 = (b.i0 - b.G2 * b.v0) / b.Cb, B = (dv0 + al * A) / wd, t = 0.2e-3;
    const bo = dcdc({ ...i, Vin: 28, Vout: 56, R_L_mohm: 0, R_on_mohm: 0, V_d: 0, t_sw_ns: 0 }), sw = dcdc({ ...i, R_L_mohm: 0, R_on_mohm: 0, V_d: 0, t_sw_ns: 40 });
    return [
      N.check('RLC bus step response at t = 0.2 ms', b.v[4000], b.vss + Math.exp(-al * t) * (A * Math.cos(wd * t) + B * Math.sin(wd * t)), 1e-9, 'Closed-form second-order response'),
      N.check('Ideal buck duty cycle D = Vout/Vin', id.D, 28 / 270, 1e-12, 'Volt-second balance'),
      N.check('Ideal boost duty cycle D = 1 − Vin/Vout', bo.D, 0.5, 1e-12, 'Volt-second balance'),
      N.check('Ideal converter is lossless', id.eff, 1, 1e-12, 'Power balance'),
      N.check('Switching loss ½·Vin·I·(t_rise + t_fall)·f', sw.Psw, 0.5 * 270 * (1000 / 28) * 40e-9 * 1e5, 1e-12, 'Triangular voltage–current overlap of a hard-switched transition'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.cpl_margin < 2) out.push({ severity: i.load_type !== 'constant power' ? 'advise' : o.cpl_margin < 1 ? 'critical' : 'warn', title: 'Low stability margin with constant-power loads', detail: `V²RC/(L·P) = ${o.cpl_margin.toFixed(2)}; ${i.load_type !== 'constant power' ? 'the load in this run is resistive, but ' : ''}regulated converters and motor drives present negative incremental resistance.`, action: 'Increase bus capacitance, add an RC damper, shorten the line, or limit converter input bandwidth; confirm with an impedance-ratio (Middlebrook) assessment.', basis: 'Linearised constant-power-load criterion P < V²·R·C/L' });
    if (o.step_undershoot_pct > 10) out.push({ severity: 'advise', title: 'Large bus undershoot on load application', detail: `${o.step_undershoot_pct.toFixed(1)}% below the new steady state (damping ratio ${o.bus_zeta.toFixed(2)}).`, action: 'Add bus capacitance or soft-start the load, and compare the excursion with the power-quality transient limits of the applicable standard.', basis: 'Second-order RLC response' });
    if (o.conv_eff < 0.95) out.push({ severity: 'advise', title: 'Converter efficiency below 95 %', detail: `${(100 * o.conv_eff).toFixed(1)}% with ${o.conv_loss_W.toPrecision(3)} W of heat.`, action: 'Lower on-resistance or switching frequency, use wide-band-gap devices, or split a large conversion ratio into two stages.', basis: 'Loss breakdown of this run' });
    if (!out.length) out.push({ severity: 'info', title: 'Converter and bus behave well at this point', detail: `Efficiency ${(100 * o.conv_eff).toFixed(1)}%, ripple ${o.ripple_V_pct.toFixed(2)}%, undershoot ${o.step_undershoot_pct.toFixed(1)}%.`, action: 'Re-check at minimum input voltage and maximum load, where current stress and ripple peak.', basis: 'Averaged model' });
    return out;
  },
};

export default {
  id: 'electrical', n: 19,
  tagline: 'Can the electrical system deliver the power — batteries, drives, buses and cables — safely, efficiently and for how long?',
  analyses: [battery, ageing, motor, network, cable, hybrid, converter],
  consumes: [
    { from: 'rotorcraft', keys: ['hover_power_W', 'P_cruise_W'], why: 'Hover and cruise shaft power of rotor-borne electric aircraft' },
    { from: 'propeller', keys: ['eta_prop'], why: 'Cruise shaft power of electric fixed-wing aircraft' },
    { from: 'propulsion', keys: ['eta_thermal'], why: 'Engine efficiency for the hybrid-electric study' },
  ],
  provides: [
    { key: 'batt_soc_end', label: 'End-of-mission state of charge', unit: '-' }, { key: 'batt_temp_K', label: 'Peak battery temperature', unit: 'K' }, { key: 'bus_V_min', label: 'Minimum bus voltage', unit: 'V' },
    { key: 'elec_losses_W', label: 'Electrical losses in cruise', unit: 'W' }, { key: 'endurance_elec_min', label: 'Battery endurance', unit: 'min' }, { key: 'motor_eff', label: 'Machine efficiency', unit: '-' }, { key: 'energy_used_kWh', label: 'Battery energy used', unit: 'kWh' },
  ],
  handoff: [
    { model: 'Doyle–Fuller–Newman and single-particle battery models', why: 'Need electrode-level parameters (diffusivities, kinetics, porosities) that are not available at this level and a coupled PDE solve; the equivalent circuit is used instead', tool: 'Electrochemical battery simulator (pseudo-2D), calibrated to cell tests' },
    { model: 'Inverter and converter switching-level simulation', why: 'Device waveforms, PWM harmonics, dead-time and EMI need sub-microsecond circuit simulation; averaged loss models are used here', tool: 'Circuit simulator (SPICE-class or power-electronics simulator)' },
    { model: 'Maxwell-equation (finite-element) electromagnetic machine analysis', why: 'Saturation, cogging, magnet eddy and AC winding losses depend on the detailed 2-D/3-D geometry', tool: 'Electromagnetic FE package with coupled thermal analysis' },
    { model: 'Induction-machine and wound-field generator models', why: 'Only the permanent-magnet synchronous machine is implemented; generators are represented by an equivalent motoring point', tool: 'Machine design software or equivalent-circuit tools with test data' },
    { model: 'AC three-phase load flow, power quality and harmonics', why: 'The network is solved as a DC-equivalent real-power system; reactive power, unbalance, frequency transients and harmonic distortion are not modelled', tool: 'Power-system analysis tool with aircraft power-quality standard checks' },
    { model: 'Protection coordination and arc-fault behaviour', why: 'Requires device time–current curves and arc models; only bolted-fault levels and ratings are estimated', tool: 'Protection coordination software and laboratory fault testing' },
    { model: 'Detailed thermal-runaway chemistry, venting and 3-D propagation', why: 'Multi-reaction kinetics, gas venting and combustion are test-driven; a lumped one-step model is used', tool: 'Accelerating-rate calorimetry plus 3-D thermal-abuse CFD/FE' },
    { model: 'Electrical–thermal–structural coupling', why: 'Only lumped electro-thermal coupling is included', tool: 'Multiphysics FE; Suite 12 for system-level heat rejection' },
  ],
};
