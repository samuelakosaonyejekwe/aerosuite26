// Suite 18 — Hydraulic, Pneumatic and Mechanical Systems.
// Hydraulic supply network, servo-actuator sizing and dynamics (hydraulic or electromechanical), pressure surge and
// accumulators, compressible duct and orifice flow, gear transmissions, rolling bearings and shafts, lubrication films.

import * as N from '../core/numerics.js';
import { isa, G0, R_AIR, GAMMA } from '../core/atmosphere.js';
import { METALS, FLUIDS } from '../data/materials.js';

// ---- shared helpers --------------------------------------------------------------------------
const num = (key, label, unit, def, min, max, group, help, x) => ({ key, label, unit, default: def, min, max, group, ...(help ? { help } : {}), ...x });
const sel = (key, label, options, def, group, help) => ({ key, label, type: 'select', options, default: def, group, ...(help ? { help } : {}) });
const kp = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const thin = (a, n = 300) => { if (a.length <= n) return a; const s = (a.length - 1) / (n - 1); return N.range(n, (j) => a[Math.round(j * s)]); };
const HYD = ['Skydrol LD-4', 'MIL-PRF-5606 hydraulic'], METAL_NAMES = Object.keys(METALS);
const hasHyd = (c) => (c.systems.hyd_p_Pa > 0 ? true : 'This aircraft has no central hydraulic system (system pressure is 0). Use the actuator analysis in electromechanical mode, and the gearbox, bearing and lubrication analyses.');
const OIL = FLUIDS['MIL-PRF-23699 oil'];
/**
 * Shaft power [W] through one transmission. Suite 7 publishes P_shaft_W for all engines together: a helicopter combines
 * them in one main gearbox, every other layout has one gearbox and shaft line per engine.
 */
const shaftPower = (c, up) => { const n = Math.max(1, c.prop.n_eng), total = up?.propulsion?.P_shaft_W > 0 ? up.propulsion.P_shaft_W : c.prop.P0_W * n; return c.meta.type === 'helicopter' ? total : total / n; };

/** Darcy friction factor: laminar 64/Re, Swamee–Jain turbulent, linear blend through transition (2000–4000). */
export function fDarcy(Re, rr = 0) {
  if (!(Re > 1e-12)) return 0;
  const tur = 0.25 / Math.log10(rr / 3.7 + 5.74 / Math.max(Re, 2000) ** 0.9) ** 2;
  if (Re <= 2000) return 64 / Re;
  if (Re >= 4000) return tur;
  const w = (Re - 2000) / 2000; return (1 - w) * (64 / Re) + w * tur;
}
/** Pressure loss [Pa] in a round line for volume flow Q (Darcy–Weisbach plus minor losses K). */
export function pipeDp(Q, D, L, rho, mu, rough = 0, K = 0) {
  const A = (Math.PI * D * D) / 4, v = Math.abs(Q) / A, Re = (rho * v * D) / mu;
  return N.sign(Q) * (fDarcy(Re, rough / D) * (L / D) + K) * 0.5 * rho * v * v;
}

// ---- 1. hydraulic supply network ---------------------------------------------------------------
const CONS = ['Primary flight controls', 'Secondary controls (flaps, gear)', 'Utility and brakes'];
function hydNet(i) {
  const fl = FLUIDS[i.fluid] || FLUIDS[HYD[0]], rho = fl.rho, mu = fl.mu * i.visc_factor, Qr = i.Q_rated;
  // pressure-compensated variable-displacement pump: slight droop to rated flow, then collapse at full displacement
  const pPump = (Q) => (Q <= Qr ? i.p_set * (1 - (i.droop * Q) / Qr) : i.p_set * (1 - i.droop) * Math.max(0, 1 - (Q - Qr) / (0.02 * Qr)));
  const br = [[i.Qv1, i.u1, i.lf1], [i.Qv2, i.u2, i.lf2], [i.Qv3, i.u3, i.lf3]].map(([Qv, u, lf]) => ({ Kv: i.p_set / 3 / Math.max(Qv * Math.max(u, 1e-9), 1e-15) ** 2, u, pl: lf * i.p_set }));
  const lineB = (Q) => 2 * pipeDp(Q, i.D_b, i.L_b, rho, mu, i.rough, i.K_fit);
  const qBranch = (b, dp) => { const av = dp - b.pl; if (b.u <= 0 || av <= 0) return 0; const g = (Q) => b.Kv * Q * Q + lineB(Q) - av, ub = Math.sqrt(av / b.Kv); return g(ub) <= 0 ? ub : N.brent(g, 0, ub, 1e-13 * ub + 1e-18); };
  const leak = (pm) => (i.Q_leak * Math.max(pm, 0)) / i.p_set;
  const state = (Q) => { const po = pPump(Q), pm = po - pipeDp(Q, i.D_p, i.L_p, rho, mu, i.rough, i.K_fit), ql = Math.min(Q, leak(pm)), pr = i.p_res + pipeDp(Q - ql, i.D_r, i.L_r, rho, mu, i.rough, i.K_fit), qs = br.map((b) => qBranch(b, pm - pr)); return { po, pm, pr, ql, qs, F: N.sum(qs) + leak(pm) - Q }; };
  const Qmax = 1.02 * Qr, Q = state(0).F <= 0 ? 0 : state(Qmax).F >= 0 ? Qmax : N.brent((q) => state(q).F, 0, Qmax, 1e-12 * Qr), s = state(Q);
  // system curve: pump-outlet pressure the consumers would need to pass a given total flow
  const sysCurve = (Qt) => { const g = (dp) => N.sum(br.map((b) => qBranch(b, dp))) - Qt, hi = 3 * i.p_set; if (g(hi) < 0) return NaN; const dp = Qt <= 0 ? N.amin(br.filter((b) => b.u > 0).map((b) => b.pl).concat([hi])) : N.brent(g, 0, hi, 1); return dp + i.p_res + pipeDp(Qt, i.D_r, i.L_r, rho, mu, i.rough, i.K_fit) + pipeDp(Qt, i.D_p, i.L_p, rho, mu, i.rough, i.K_fit); };
  return { ...s, Q, rho, mu, fl, br, pPump, sysCurve, lineB, dpP: s.po - s.pm, dpR: s.pr - i.p_res };
}
const network = {
  id: 'network', title: 'Hydraulic power generation and distribution network', fidelity: 'numerical',
  summary: 'Operating point of a pressure-compensated pump feeding three consumer groups through pressure and return lines: flows, pressures, the pressure-loss budget, pump power and heat rejected to the fluid.',
  equations: ['Continuity equation', 'Bernoulli equation', 'Darcy–Weisbach equation', 'Hagen–Poiseuille equation', 'Orifice flow equations', 'Mechanical power balance equations'],
  applicable: hasHyd,
  inputs: [
    sel('fluid', 'Hydraulic fluid', HYD, HYD[0], 'Fluid', 'Phosphate-ester (Skydrol) is fire-resistant; mineral MIL-PRF-5606 is lighter but flammable'),
    num('visc_factor', 'Viscosity multiplier', '-', 1, 0.2, 500, 'Fluid', '1 at normal operating temperature; roughly 10–50 during a −40 °C cold start (typical)'),
    num('p_set', 'Pump compensator (system) pressure', 'Pa', 20.7e6, 1e5, 60e6, 'Pump', '20.7 MPa (3000 psi) classic, 34.5 MPa (5000 psi) on newer large aircraft'),
    num('Q_rated', 'Pump rated flow', 'm³/s', 2.3e-3, 1e-7, 0.1, 'Pump', '1 L/min = 1.667e-5 m³/s'),
    num('droop', 'Pressure droop at rated flow', '-', 0.03, 0, 0.3, 'Pump'), num('eta_pump', 'Pump overall efficiency', '-', 0.87, 0.3, 0.98, 'Pump'),
    num('p_res', 'Reservoir pressure', 'Pa', 0.35e6, 0, 5e6, 'Pump', 'Bootstrap or air-pressurised reservoirs ≈ 0.3–0.5 MPa'),
    num('Q_leak', 'Internal leakage at system pressure', 'm³/s', 7e-5, 0, 0.01, 'Pump', 'Pump case drain plus servo-valve quiescent leakage'),
    num('D_p', 'Pressure line bore', 'm', 0.016, 0.001, 0.2, 'Lines'), num('L_p', 'Pressure line length', 'm', 15, 0.01, 500, 'Lines'),
    num('D_r', 'Return line bore', 'm', 0.022, 0.001, 0.2, 'Lines'), num('L_r', 'Return line length', 'm', 15, 0.01, 500, 'Lines'),
    num('D_b', 'Consumer branch line bore', 'm', 0.01, 0.0005, 0.2, 'Lines'), num('L_b', 'Consumer branch length (each way)', 'm', 6, 0.01, 500, 'Lines'),
    num('rough', 'Tube wall roughness', 'm', 1.5e-6, 0, 1e-3, 'Lines', 'Drawn tubing ≈ 1.5 µm'), num('K_fit', 'Fitting loss coefficient per line', '-', 4, 0, 200, 'Lines', 'Sum of bends, tees, filters and quick-disconnects'),
    ...CONS.flatMap((nm, j) => [
      num(`Qv${j + 1}`, `${nm}: valve rated flow`, 'm³/s', [1.2e-3, 1.5e-3, 5e-4][j], 1e-8, 0.1, 'Consumers', 'Flow with the valve fully open at a valve drop of one third of system pressure'),
      num(`u${j + 1}`, `${nm}: valve opening`, '-', [0.5, 0.3, 0][j], 0, 1, 'Consumers'), num(`lf${j + 1}`, `${nm}: load pressure / system pressure`, '-', [0.5, 0.6, 0.3][j], 0, 1, 'Consumers')]),
  ],
  defaults: (c) => {
    const p = c.systems.hyd_p_Pa, Q = Math.max(2e-6, 3e-8 * c.mass.mtow_kg), len = Math.max(1, c.fuselage.len_m || 5), dOf = (q, v) => Math.max(0.003, Math.sqrt((4 * q) / (Math.PI * v)));
    return { p_set: p, Q_rated: Q, Q_leak: 0.03 * Q, fluid: c.mass.mtow_kg > 5700 && c.meta.type !== 'helicopter' ? HYD[0] : HYD[1], D_p: dOf(Q, 7), D_r: dOf(Q, 4), D_b: dOf(0.5 * Q, 8), L_p: 0.4 * len, L_r: 0.4 * len, L_b: 0.15 * len, Qv1: 0.55 * Q, Qv2: 0.65 * Q, Qv3: 0.2 * Q };
  },
  run(i) {
    const r = hydNet(i), Ph = (r.po - i.p_res) * r.Q, Ps = Ph / i.eta_pump, useful = N.sum(r.br.map((b, j) => b.pl * r.qs[j])), heat = Ps - useful, dT = r.Q > 0 ? heat / (r.rho * r.Q * r.fl.cp) : 0, warnings = [];
    const j0 = N.argmax(r.qs), dpB = r.lineB(r.qs[j0]), dpV = r.br[j0].Kv * r.qs[j0] ** 2, lineDp = r.dpP + r.dpR + dpB, vP = r.Q / ((Math.PI * i.D_p ** 2) / 4), Re = (r.rho * vP * i.D_p) / r.mu;
    if (r.Q >= i.Q_rated * 0.999) warnings.push('Demand reaches the pump rated flow: the pump is at full displacement and system pressure sags. Add pump capacity, an accumulator for peaks, or priority valves.');
    if (vP > 9) warnings.push(`Pressure-line velocity ${vP.toFixed(1)} m/s exceeds the usual 6–9 m/s guideline: expect high losses and surge pressures.`);
    if (i.visc_factor > 5 && Re < 2000) warnings.push('Cold, viscous fluid: line flow is laminar and losses scale with viscosity; consumer rates are reduced until the fluid warms.');
    const Qs = N.linspace(0, 1.02 * i.Q_rated, 60), lpm = (q) => q * 6e4;
    return {
      kpis: [
        kp('pump_power_W', 'Pump shaft power', Ps, 'W'), kp('hyd_power_W', 'Hydraulic power delivered by the pump', Ph, 'W'),
        kp('pump_flow_m3s', 'Pump flow', r.Q, 'm³/s', r.Q < 0.95 * i.Q_rated ? 'ok' : 'warn', `${lpm(r.Q).toFixed(1)} L/min of ${lpm(i.Q_rated).toFixed(1)} rated`),
        kp('p_pump_Pa', 'Pump outlet pressure', r.po, 'Pa'), kp('p_manifold_Pa', 'Supply manifold pressure', r.pm, 'Pa', r.pm > 0.9 * i.p_set ? 'ok' : 'warn'), kp('p_return_Pa', 'Return manifold pressure', r.pr, 'Pa'),
        kp('line_dp_Pa', 'Line pressure loss on the busiest consumer path', lineDp, 'Pa', lineDp < 0.1 * i.p_set ? 'ok' : 'warn', 'Pressure + branch + return lines; guideline < 10% of system pressure'),
        kp('useful_power_W', 'Power delivered to the loads', useful, 'W'), kp('heat_rejection_W', 'Heat rejected to the fluid', heat, 'W'),
        kp('system_eff', 'Shaft-to-load efficiency', Ps > 0 ? useful / Ps : 0, '-'), kp('fluid_dT_K', 'Fluid temperature rise per pass', dT, 'K'),
        kp('Re_pressure_line', 'Pressure-line Reynolds number', Re, '-'), kp('continuity_residual', 'Flow continuity residual', Math.abs(r.F) / Math.max(i.Q_rated, 1e-15), '-'),
      ],
      plots: [
        { type: 'line', title: 'Pump characteristic and system demand', xlabel: 'Flow [L/min]', ylabel: 'Pressure at pump outlet [MPa]', series: [{ name: 'Pump (pressure-compensated)', x: Qs.map(lpm), y: Qs.map((q) => r.pPump(q) / 1e6) }, { name: 'System demand', x: Qs.map(lpm), y: Qs.map((q) => r.sysCurve(Math.max(0, q - r.ql)) / 1e6) }, { name: 'Operating point', x: [lpm(r.Q)], y: [r.po / 1e6], style: 'points' }] },
        { type: 'bar', title: 'Pressure budget along the busiest consumer path', ylabel: 'Pressure [MPa]', categories: ['Pressure line', 'Branch lines', 'Control valve', 'Load', 'Return line', 'Reservoir'], series: [{ name: 'Pressure used', y: [r.dpP, dpB, dpV, r.br[j0].pl, r.dpR, i.p_res].map((v) => v / 1e6) }] },
      ],
      tables: [{ title: 'Consumer flows', columns: ['Consumer', 'Flow [L/min]', 'Valve drop [MPa]', 'Load pressure [MPa]', 'Load power [kW]'], rows: CONS.map((nm, j) => [nm, lpm(r.qs[j]), (r.br[j].Kv * r.qs[j] ** 2) / 1e6, r.br[j].pl / 1e6, (r.br[j].pl * r.qs[j]) / 1e3]) }],
      warnings, models: ['Pressure-compensated pump characteristic', 'Darcy–Weisbach lines (laminar 64/Re, Swamee–Jain turbulent)', 'Square-law valve orifices', 'Nodal continuity solved by bracketed root finding'],
      assumptions: ['Steady, incompressible, isothermal flow', 'One pump and one lumped consumer per group; identical supply and return branch runs', 'Leakage proportional to pressure', 'All throttling and leakage losses become heat in the fluid', 'Default pump flow, line bores, valve sizes, leakage and loss coefficients are scaled typical values, not data for a specific aircraft'],
    };
  },
  verify() {
    const fl = FLUIDS[HYD[0]], Q = 1e-5, D = 0.01, L = 5, Re = 1e5, rr = 1e-4;
    let f = 0.02; for (let k = 0; k < 60; k++) f = 1 / (-2 * Math.log10(rr / 3.7 + 2.51 / (Re * Math.sqrt(f)))) ** 2;
    // one consumer, no line losses, no leakage, no droop: Q = sqrt((p_set − p_res − p_load)/Kv)
    const i = { fluid: HYD[0], visc_factor: 1, p_set: 20e6, Q_rated: 1e-3, droop: 0, eta_pump: 1, p_res: 0, Q_leak: 0, D_p: 0.5, L_p: 0.01, D_r: 0.5, L_r: 0.01, D_b: 0.5, L_b: 0.01, rough: 0, K_fit: 0, Qv1: 5e-4, u1: 1, lf1: 0.25, Qv2: 1e-4, u2: 0, lf2: 0, Qv3: 1e-4, u3: 0, lf3: 0 }, r = hydNet(i), o = N.kv(network.run(i));
    return [
      N.check('Laminar line loss equals Hagen–Poiseuille 128μLQ/(πD⁴)', pipeDp(Q, D, L, fl.rho, fl.mu), (128 * fl.mu * L * Q) / (Math.PI * D ** 4), 1e-10, 'Hagen–Poiseuille'),
      N.check('Swamee–Jain friction factor against Colebrook', fDarcy(Re, rr), f, 0.015, 'Colebrook–White solved iteratively'),
      N.check('Single-orifice operating point', r.Q, 5e-4 * Math.sqrt(15e6 / (20e6 / 3)), 1e-5, 'Orifice law Q ∝ √Δp'),
      N.check('Power balance: shaft power = load power + heat', o.pump_power_W, o.useful_power_W + o.heat_rejection_W, 1e-12, 'First law'),
      N.check('Flow continuity at the manifold', o.continuity_residual, 0, 1e-8, 'ΣQ = 0'),
    ];
  },
  calibration: { params: [{ key: 'K_fit', min: 0, max: 100 }, { key: 'Q_leak', min: 0, max: 1e-3 }, { key: 'droop', min: 0, max: 0.2 }], sweep: 'u1', target: 'p_manifold_Pa', note: 'Supply measured manifold pressure against valve opening (or demanded flow) from an iron-bird or pump rig test.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.p_manifold_Pa < 0.85 * i.p_set) out.push({ severity: 'warn', title: 'System pressure sags under this demand', detail: `Manifold pressure is ${(o.p_manifold_Pa / 1e6).toFixed(1)} MPa against a ${(i.p_set / 1e6).toFixed(1)} MPa setting.`, action: 'Increase pump displacement, add an accumulator for transient peaks, or enlarge the pressure line.', basis: 'Pump–system operating point' });
    if (o.system_eff < 0.4 && o.pump_power_W > 0) out.push({ severity: 'advise', title: 'Most pump power is throttled into heat', detail: `Only ${(100 * o.system_eff).toFixed(0)}% of shaft power reaches the loads; ${(o.heat_rejection_W / 1e3).toFixed(1)} kW heats the fluid and must be removed by the fuel- or air-cooled heat exchanger.`, action: 'Valve-controlled hydraulics waste power at part load. Electro-hydrostatic or electromechanical actuators draw power on demand and remove the continuous engine off-take, saving fuel and CO₂; pass the heat load to Suite 12.', basis: 'Throttling loss (p_s − p_load)·Q' });
    if (i.fluid === HYD[1]) out.push({ severity: 'info', title: 'Mineral hydraulic fluid selected', detail: 'MIL-PRF-5606 is flammable; fire-resistant synthetic hydrocarbon or phosphate-ester fluids are preferred for transport aircraft.', action: 'Confirm fluid choice against fire-safety, seal compatibility and disposal requirements; phosphate esters need careful handling and waste treatment.', basis: 'Fluid selection' });
    return out;
  },
};

// ---- 2. actuator sizing and dynamics -----------------------------------------------------------
const KIND = ['Hydraulic servo-actuator', 'Electromechanical (ball-screw) actuator'];
/** Valve-controlled double-acting actuator: state [x, v, pA, pB]; u(t,y) is the valve opening in [−1, 1]. */
function hydActRhs(p, u) {
  const dq = (dp) => N.sign(dp) * Math.sqrt((2 * Math.abs(dp)) / p.rho);
  return (t, y) => {
    const [x, v, pA, pB] = y, uu = typeof u === 'function' ? u(t, y) : u, a = p.CdW * Math.abs(uu);
    const QA = uu >= 0 ? a * dq(p.ps - pA) : -a * dq(pA - p.pr), QB = uu >= 0 ? a * dq(pB - p.pr) : -a * dq(p.ps - pB);
    const VA = Math.max(0.2 * p.Vd, p.Vd + p.A * (p.S / 2 + x)), VB = Math.max(0.2 * p.Vd, p.Vd + p.A * (p.S / 2 - x)), ql = p.Cl * (pA - pB);
    return [v, (p.A * (pA - pB) - p.c * v - p.k * x - p.Fext) / p.m, (p.beta / VA) * (QA - p.A * v - ql), (p.beta / VB) * (-QB + p.A * v + ql)];
  };
}
function actSize(i) {
  const hyd = i.kind === KIND[0], Freq = i.HM / i.arm / Math.max(1, i.n_act), Fstall = i.F_margin * Freq, S = 2 * i.arm * Math.sin(N.rad(i.defl_deg)), vReq = N.rad(i.rate_dps) * i.arm;
  const o = { hyd, Freq, Fstall, S, vReq };
  if (hyd) {
    const fl = FLUIDS[i.fluid] || FLUIDS[HYD[0]], dps = i.p_sys - i.p_ret, A = Fstall / dps, vNL = vReq / Math.sqrt(1 - 1 / i.F_margin), CdW = (A * vNL) / Math.sqrt(dps / fl.rho), Vd = i.dead_frac * A * S, Vt = 2 * Vd + A * S;
    const kh = (4 * fl.bulk * A * A) / Vt, m = i.m_load, wh = Math.sqrt(kh / m);
    Object.assign(o, { fl, A, vNL, CdW, Vd, Vt, kh, wh, bore: Math.sqrt((4 * A) / Math.PI / (1 - i.rod_frac ** 2)), Qpk: A * vNL, par: { A, S, Vd, beta: fl.bulk, rho: fl.rho, ps: i.p_sys, pr: i.p_ret, CdW, m, c: i.c_load, k: i.k_load, Fext: i.load_frac * Freq, Cl: (i.leak_frac * A * vNL) / dps } });
  } else {
    const vNL = vReq / (1 - 1 / i.F_margin), wm = (i.motor_rpm * 2 * Math.PI) / 60, G = (wm * i.lead) / (2 * Math.PI * vNL), Tpk = (Fstall * i.lead) / (2 * Math.PI * i.eta_screw * G), Jm = i.Jm_per_Nm * Tpk;
    const mRefl = (Jm * ((2 * Math.PI * G) / i.lead) ** 2) / i.eta_screw;
    Object.assign(o, { vNL, wm, G, Tpk, Jm, mRefl, meq: i.m_load + mRefl, Ppk: (Fstall * vNL) / 4 / (i.eta_screw * i.eta_motor), Fext: i.load_frac * Freq });
  }
  return o;
}
function actSim(i, z) {
  const xc = i.step_frac * z.S / 2, sat = (e) => N.clamp((i.Kv * e) / z.vNL, -1, 1); let r, F;
  if (z.hyd) {
    const p = z.par, pm = 0.5 * (p.ps + p.pr), y0 = [0, 0, pm + p.Fext / (2 * p.A), pm - p.Fext / (2 * p.A)];
    r = N.rk45(hydActRhs(p, (t, y) => sat(xc - y[0])), 0, y0, i.t_end, { rtol: 1e-6, atol: 1e-9, hmax: i.t_end / 200, maxSteps: 60000 }); F = r.y.map((y) => p.A * (y[2] - y[3]));
  } else {
    // PI position loop (integral time 4/Kv, clamped when the drive saturates) commanding motor voltage
    const uOf = (y) => (i.Kv * (xc - y[0] + (i.Kv / 4) * y[2])) / z.vNL, fm = (y) => z.Fstall * N.clamp(N.clamp(uOf(y), -1, 1) - y[1] / z.vNL, -1, 1);
    r = N.rk45((t, y) => [y[1], (fm(y) - i.c_load * y[1] - i.k_load * y[0] - z.Fext) / z.meq, Math.abs(uOf(y)) < 1 ? xc - y[0] : 0], 0, [0, 0, 0], i.t_end, { rtol: 1e-7, atol: 1e-10, hmax: i.t_end / 200, maxSteps: 60000 }); F = r.y.map(fm);
  }
  const x = r.y.map((y) => y[0]), v = r.y.map((y) => y[1]), xe = x[x.length - 1], j90 = x.findIndex((q) => q >= 0.9 * xc);
  return { t: r.t, x, v, F, xc, xe, t90: j90 >= 0 ? r.t[j90] : NaN, os: xc > 0 ? Math.max(0, (N.amax(x) - xc) / xc) : 0, vmax: N.amax(v), pA: z.hyd ? r.y.map((y) => y[2]) : null, pB: z.hyd ? r.y.map((y) => y[3]) : null, steps: r.steps };
}
const actuator = {
  id: 'actuator', title: 'Actuator sizing and step response (hydraulic or electromechanical)', fidelity: 'numerical',
  summary: 'Sizes a flight-control actuator from hinge moment and rate, then simulates a position step with valve flow and fluid compressibility (hydraulic) or motor torque–speed limits and reflected inertia (electromechanical).',
  equations: ['Hydraulic actuator force equations', 'Orifice flow equations', 'Fluid compressibility equations', 'Hydraulic capacitance equations', 'Continuity equation', 'Mechanical power balance equations'],
  inputs: [
    sel('kind', 'Actuator technology', KIND, KIND[0], 'Requirement'),
    num('HM', 'Design hinge moment (total)', 'N·m', 16000, 0.001, 1e7, 'Requirement', 'Estimated from dynamic pressure, surface area and chord when a case is loaded; replace with Suite 1/4 data'),
    num('n_act', 'Actuators sharing the load', '', 2, 1, 8, 'Requirement', '', { step: 1, discrete: true }),
    num('arm', 'Horn (moment) arm', 'm', 0.12, 0.002, 2, 'Requirement'), num('defl_deg', 'Surface travel (± from neutral)', 'deg', 25, 1, 80, 'Requirement'),
    num('rate_dps', 'Required surface rate under design load', 'deg/s', 60, 0.5, 1000, 'Requirement'),
    num('F_margin', 'Stall force / design load', '-', 1.5, 1.05, 5, 'Requirement', '1.5 places the design point at maximum power transfer for a valve-controlled actuator'),
    num('m_load', 'Effective load mass at the actuator', 'kg', 1200, 1e-4, 1e6, 'Load', 'Surface inertia divided by the horn arm squared'), num('k_load', 'Load spring rate (aerodynamic)', 'N/m', 0, 0, 1e9, 'Load'), num('c_load', 'Load damping', 'N·s/m', 5000, 0, 1e8, 'Load'),
    num('load_frac', 'Steady opposing load / design load during the step', '-', 0.5, 0, 0.99, 'Load'),
    num('p_sys', 'Supply pressure', 'Pa', 20.7e6, 1e5, 60e6, 'Hydraulic'), num('p_ret', 'Return pressure', 'Pa', 0.35e6, 0, 5e6, 'Hydraulic'), sel('fluid', 'Hydraulic fluid', HYD, HYD[0], 'Hydraulic'),
    num('rod_frac', 'Rod diameter / bore', '-', 0.5, 0, 0.9, 'Hydraulic', 'Balanced (through-rod) cylinder'), num('dead_frac', 'Dead volume per side / swept volume', '-', 0.1, 0.01, 2, 'Hydraulic'), num('leak_frac', 'Cross-piston leakage at full pressure / peak flow', '-', 0.01, 0, 0.3, 'Hydraulic'),
    num('lead', 'Ball-screw lead', 'm/rev', 0.005, 0.0005, 0.05, 'Electromechanical'), num('eta_screw', 'Screw and gear efficiency', '-', 0.85, 0.3, 0.98, 'Electromechanical'), num('motor_rpm', 'Motor no-load speed', 'rpm', 8000, 100, 60000, 'Electromechanical'),
    num('Jm_per_Nm', 'Motor inertia per unit peak torque', 'kg·m²/(N·m)', 5e-5, 1e-7, 1e-2, 'Electromechanical', 'Typical brushless servo motors ≈ 2–10 ×10⁻⁵'), num('eta_motor', 'Motor and drive efficiency', '-', 0.88, 0.3, 0.99, 'Electromechanical'),
    num('Kv', 'Position-loop gain', '1/s', 20, 0.5, 500, 'Control', 'Velocity constant: commanded rate per unit position error'), num('step_frac', 'Step command / half stroke', '-', 0.2, 0.01, 0.9, 'Control'),
    num('t_end', 'Simulated time', 's', 0.8, 0.02, 30, 'Numerics'),
  ],
  defaults: (c, up, d) => {
    const heli = c.meta.type === 'helicopter', W = d.W, q = 0.5 * 1.225 * (c.aero.Vmo_ms || c.flight.V_ms || 30) ** 2, Sh = c.htail.S_m2, ce = Sh > 0 && c.htail.b_m > 0 ? (0.3 * Sh) / c.htail.b_m : 0;
    // hinge moment: Ch ≈ 0.12 at design deflection (typical, partly balanced surface); rotor servos carry a fraction of weight
    const HM = heli ? 0.3 * W * 0.1 : Sh > 0 ? q * c.controls.Se_Sh * Sh * ce * 0.12 : 0.05 * W * 0.1, arm = heli ? 0.1 : N.clamp(0.25 * ce || 0.03, 0.01, 0.15), n = heli ? 3 : c.mass.mtow_kg > 5700 ? 2 : 1, hyd = c.systems.hyd_p_Pa > 0;
    return { kind: hyd ? KIND[0] : KIND[1], HM, arm, n_act: n, defl_deg: heli ? 15 : c.controls.de_max_deg, rate_dps: c.controls.rate_max_dps, p_sys: hyd ? c.systems.hyd_p_Pa : undefined, fluid: c.mass.mtow_kg > 5700 && !heli ? HYD[0] : HYD[1], m_load: Math.max(0.02, (0.02 * HM) / arm / n), c_load: Math.max(0.05, (4 * HM) / arm / n / 1000) };
  },
  run(i, ctx) {
    const z = actSize(i), s = actSim(i, z), warnings = [], vLoad = z.hyd ? z.vNL * Math.sqrt(1 - z.Freq / z.Fstall) : z.vNL * (1 - z.Freq / z.Fstall), ess = s.xc - s.xe;
    const Fs = N.linspace(0, z.Fstall, 40), env = Fs.map((F) => (z.hyd ? z.vNL * Math.sqrt(1 - F / z.Fstall) : z.vNL * (1 - F / z.Fstall)));
    const kpis = [
      kp('actuator_force_N', 'Stall (maximum) force per actuator', z.Fstall, 'N'), kp('design_load_N', 'Design load per actuator', z.Freq, 'N'),
      kp('actuator_rate_ms', 'No-load piston rate', z.vNL, 'm/s'), kp('rate_at_load_ms', 'Rate under design load', vLoad, 'm/s', vLoad >= 0.999 * z.vReq ? 'ok' : 'bad', `Required ${z.vReq.toFixed(4)} m/s`),
      kp('stroke_m', 'Stroke', z.S, 'm'),
    ];
    if (z.hyd) {
      const zh = 0.5 * (i.c_load / (z.par.m * z.wh)) + 0.1, KvMax = 2 * zh * z.wh;
      kpis.push(kp('piston_area_m2', 'Piston area', z.A, 'm²'), kp('bore_m', 'Cylinder bore', z.bore, 'm'), kp('peak_flow_m3s', 'Peak flow demand', z.Qpk, 'm³/s', undefined, `${(z.Qpk * 6e4).toFixed(1)} L/min`),
        kp('hyd_stiffness_Npm', 'Hydraulic stiffness (centred)', z.kh, 'N/m'), kp('hyd_freq_Hz', 'Hydraulic natural frequency', z.wh / (2 * Math.PI), 'Hz'), kp('peak_power_W', 'Peak hydraulic power demand', i.p_sys * z.Qpk, 'W'),
        kp('Kv_limit', 'Stability limit on loop gain (2ζω)', KvMax, '1/s', i.Kv < 0.7 * KvMax ? 'ok' : 'warn', 'With an assumed 0.1 valve/structural damping ratio'));
      if (i.Kv > KvMax) warnings.push('Loop gain exceeds the 2ζ·ω_h stability guideline: expect a lightly damped or unstable response. Lower the gain or add pressure feedback.');
    } else {
      kpis.push(kp('motor_torque_Nm', 'Motor peak torque', z.Tpk, 'N·m'), kp('gear_ratio', 'Motor-to-screw gear ratio', z.G, '-', z.G >= 0.5 && z.G <= 20 ? 'ok' : 'warn'), kp('reflected_mass_kg', 'Motor inertia reflected to the output', z.mRefl, 'kg', undefined, `${(z.mRefl / Math.max(i.m_load, 1e-9)).toFixed(1)} × load mass`),
        kp('peak_power_W', 'Peak electrical power', z.Ppk, 'W'), kp('hold_torque_Nm', 'Motor torque to hold the steady load', (z.Fext * i.lead) / (2 * Math.PI * i.eta_screw * z.G), 'N·m', undefined, 'Held by motor current (heat) unless a brake or irreversible screw is fitted'));
      if (z.G < 0.5 || z.G > 20) warnings.push('The gear ratio needed to match motor speed and screw lead is outside 0.5–20: change the screw lead or motor speed.');
    }
    kpis.push(kp('t_rise_s', 'Time to 90% of the step', s.t90, 's'), kp('overshoot_pct', 'Overshoot', 100 * s.os, '%', s.os < 0.1 ? 'ok' : 'warn'), kp('ss_error_m', 'Position error at end of run', ess, 'm'), kp('rate_peak_ms', 'Peak rate during the step', s.vmax, 'm/s'));
    if (Number.isNaN(s.t90)) warnings.push('The actuator did not reach 90% of the commanded step within the simulated time.');
    const plots = [
      { type: 'line', title: 'Position step response', xlabel: 'Time [s]', ylabel: 'Position [mm]', series: [{ name: 'Actuator position', x: thin(s.t), y: thin(s.x).map((v) => v * 1e3) }], annotations: [{ y: s.xc * 1e3, label: 'Command' }] },
      { type: 'line', title: 'Actuator output force', xlabel: 'Time [s]', ylabel: 'Force [kN]', series: [{ name: z.hyd ? 'Pressure force A·(pA − pB)' : 'Motor force at the screw', x: thin(s.t), y: thin(s.F).map((v) => v / 1e3) }] },
      { type: 'line', title: 'Load–speed envelope', xlabel: 'Opposing load [kN]', ylabel: 'Rate [m/s]', series: [{ name: 'Maximum rate', x: Fs.map((v) => v / 1e3), y: env }, { name: 'Design point', x: [z.Freq / 1e3], y: [z.vReq], style: 'points' }] },
    ];
    if (z.hyd) plots.push({ type: 'line', title: 'Chamber pressures', xlabel: 'Time [s]', ylabel: 'Pressure [MPa]', series: [{ name: 'Chamber A', x: thin(s.t), y: thin(s.pA).map((v) => v / 1e6) }, { name: 'Chamber B', x: thin(s.t), y: thin(s.pB).map((v) => v / 1e6) }] });
    const noHyd = ctx?.case && !(ctx.case.systems.hyd_p_Pa > 0);
    return {
      kpis, plots, warnings, outputs: noHyd ? { pump_power_W: 0, line_dp_Pa: 0 } : {},
      models: z.hyd ? ['Four-way critical-centre servo-valve orifice flow', 'Compressible chamber volumes (bulk modulus)', 'Mass–spring–damper load, adaptive RK45'] : ['Linear motor torque–speed characteristic with torque saturation', 'Ball-screw kinematics with reflected rotor inertia', 'Mass–spring–damper load, adaptive RK45'],
      assumptions: z.hyd ? ['Balanced double-acting cylinder, constant supply and return pressure', 'Ideal valve (no lag, no overlap), proportional position loop', 'Constant bulk modulus (no entrained air)', 'No seal friction or end-stop contact', 'Default hinge moment, load mass, damping and loop gain are class-level estimates (typical values)'] : ['Rigid screw and nut (no backlash or compliance)', 'Proportional–integral position loop commanding motor voltage', 'Constant efficiency, no thermal limit on motor torque', 'Default hinge moment, load mass, motor and screw data are class-level estimates (typical values)'],
    };
  },
  verify() {
    const b = { kind: KIND[0], HM: 1000, n_act: 1, arm: 0.1, defl_deg: 30, rate_dps: 60, F_margin: 1.5, m_load: 50, k_load: 0, c_load: 0, load_frac: 0, p_sys: 20e6, p_ret: 0, fluid: HYD[0], rod_frac: 0.5, dead_frac: 0.1, leak_frac: 0, lead: 0.005, eta_screw: 0.9, motor_rpm: 6000, Jm_per_Nm: 5e-5, eta_motor: 0.9, Kv: 20, step_frac: 0.2, t_end: 0.5 };
    const z = actSize(b), p = z.par, pm = 10e6;
    // (1) fully open valve, no load: steady rate = Cd·w·sqrt(Δp_s/ρ)/A
    const nl = N.rk45(hydActRhs({ ...p, c: 1 }, 1), 0, [-0.04, 0, pm, pm], 0.05, { rtol: 1e-8, atol: 1e-10 }), vss = nl.y[nl.y.length - 1][1];
    // (2) blocked valve: mass on the trapped oil columns rings at ω_h = sqrt(4βA²/(V_t·m))
    const rg = N.rk45(hydActRhs(p, 0), 0, [0, 1e-4, pm, pm], 0.2, { rtol: 1e-9, atol: 1e-12, hmax: 1e-4 }), zc = [];
    for (let k = 1; k < rg.t.length; k++) if (rg.y[k - 1][1] > 0 && rg.y[k][1] <= 0) zc.push(rg.t[k - 1] + ((rg.t[k] - rg.t[k - 1]) * rg.y[k - 1][1]) / (rg.y[k - 1][1] - rg.y[k][1]));
    const wNum = (2 * Math.PI * (zc.length - 1)) / (zc[zc.length - 1] - zc[0]);
    // (3) electromechanical: full command, no load -> no-load speed; stall force from torque and lead
    const e = { ...b, kind: KIND[1] }, ze = actSize(e), em = N.rk45((t, y) => [y[1], (ze.Fstall * N.clamp(1 - y[1] / ze.vNL, -1, 1)) / ze.meq], 0, [0, 0], 2, { rtol: 1e-8 });
    return [
      N.check('Stall force = Δp·A', z.Fstall, 20e6 * z.A, 1e-12, 'Actuator force equation'),
      N.check('No-load rate of a valve-controlled actuator', vss, z.vNL, 1e-3, 'Merritt, Hydraulic Control Systems: Q = Cd·w·sqrt(p_s/ρ)'),
      N.check('Hydraulic natural frequency from the transient', wNum, z.wh, 2e-3, 'ω_h = sqrt(4βA²/(V_t·m))'),
      N.check('Ball-screw force F = 2π·η·G·T/lead', (2 * Math.PI * e.eta_screw * ze.G * ze.Tpk) / e.lead, ze.Fstall, 1e-12, 'Screw kinematics'),
      N.check('Electromechanical no-load speed from the transient', em.y[em.y.length - 1][1], ze.vNL, 1e-4, 'Torque–speed line'),
    ];
  },
  calibration: { params: [{ key: 'c_load', min: 0, max: 1e7 }, { key: 'leak_frac', min: 0, max: 0.2 }, { key: 'dead_frac', min: 0.01, max: 2 }], sweep: 'Kv', target: 't_rise_s', note: 'Supply measured step-response rise time against loop gain from an actuator test bench (frequency response data constrain damping and stiffness).' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.rate_at_load_ms < 0.999 * N.rad(i.rate_dps) * i.arm) out.push({ severity: 'warn', title: 'Rate requirement not met under load', detail: `${o.rate_at_load_ms.toFixed(3)} m/s available at the design load.`, action: 'Increase the stall-force margin or valve/motor size.', basis: 'Load–speed envelope' });
    if (o.overshoot_pct > 15) out.push({ severity: 'advise', title: 'Lightly damped position response', detail: `${o.overshoot_pct.toFixed(0)}% overshoot.`, action: 'Reduce loop gain, add dynamic pressure feedback (hydraulic) or velocity feedback (electromechanical); check interaction with surface flutter modes in Suite 3.', basis: 'Servo-loop stability' });
    if (i.kind === KIND[0]) out.push({ severity: 'info', title: 'Power-on-demand alternative', detail: `Peak demand ${(o.peak_power_W / 1e3).toFixed(1)} kW is drawn from a continuously pressurised system.`, action: 'An electro-hydrostatic or electromechanical actuator draws power only when moving and removes central hydraulic leakage losses and fluid mass: re-run in electromechanical mode to compare.', basis: 'More-electric actuation trade' });
    else if (o.reflected_mass_kg > 5 * i.m_load) out.push({ severity: 'advise', title: 'Motor inertia dominates the load', detail: `Reflected mass is ${(o.reflected_mass_kg / i.m_load).toFixed(0)} times the load mass.`, action: 'High reflected inertia raises end-stop and jam loads; use a coarser lead or lower gear ratio, and provide jam-tolerance (dual load path or clutch).', basis: 'Reflected inertia J·(2πG/lead)²' });
    return out;
  },
};

// ---- 3. pressure surge and accumulator ---------------------------------------------------------
/** Method of characteristics for a line fed from a constant-pressure source with a closing valve at the far end. */
function moc(i) {
  const fl = FLUIDS[i.fluid] || FLUIDS[HYD[0]], rho = fl.rho, E = (METALS[i.tube] || METALS['Ti-6Al-4V']).E, e = i.wall_mm / 1e3;
  const a = Math.sqrt(fl.bulk / rho / (1 + (fl.bulk * i.D) / (E * e))), A = (Math.PI * i.D ** 2) / 4, Q0 = i.v0 * A, n = Math.max(2, Math.round(i.nReach)), dx = i.L / n, dt = dx / a, B = (rho * a) / A;
  const Re = (rho * i.v0 * i.D) / fl.mu, f = i.fric * fDarcy(Re, 1.5e-6 / i.D), R = (rho * f * dx) / (2 * i.D * A * A), dpf = R * n * Q0 * Q0, dpv0 = i.p_sys - dpf - i.p_back;
  const Tsim = Math.min(i.n_periods * ((4 * i.L) / a) + i.t_close, 30000 * dt), nt = Math.ceil(Tsim / dt);
  let p = N.range(n + 1, (j) => i.p_sys - R * j * Q0 * Q0), Q = new Array(n + 1).fill(Q0);
  const t = [0], pv = [p[n]], pmid = [p[n >> 1]], pmax = p.slice(), pmin = p.slice();
  for (let s = 1; s <= nt; s++) {
    const pn = new Array(n + 1), Qn = new Array(n + 1), tau = i.t_close > 0 ? Math.max(0, 1 - (s * dt) / i.t_close) : 0;
    for (let j = 1; j < n; j++) { const CP = p[j - 1] + B * Q[j - 1] - R * Q[j - 1] * Math.abs(Q[j - 1]), CM = p[j + 1] - B * Q[j + 1] + R * Q[j + 1] * Math.abs(Q[j + 1]); pn[j] = 0.5 * (CP + CM); Qn[j] = (CP - CM) / (2 * B); }
    pn[0] = i.p_sys; Qn[0] = (i.p_sys - (p[1] - B * Q[1] + R * Q[1] * Math.abs(Q[1]))) / B;
    const CP = p[n - 1] + B * Q[n - 1] - R * Q[n - 1] * Math.abs(Q[n - 1]) - i.p_back, Cv = dpv0 > 0 ? (tau * Q0) ** 2 / (2 * dpv0) : 0, disc = (B * Cv) ** 2 + 2 * Cv * CP;
    Qn[n] = Cv > 0 && disc > 0 ? -B * Cv + Math.sqrt(disc) : 0; pn[n] = CP + i.p_back - B * Qn[n];
    p = pn; Q = Qn; t.push(s * dt); pv.push(p[n]); pmid.push(p[n >> 1]);
    for (let j = 0; j <= n; j++) { if (p[j] > pmax[j]) pmax[j] = p[j]; if (p[j] < pmin[j]) pmin[j] = p[j]; }
  }
  return { a, A, Q0, rho, f, Re, dpf, dpv0, t, pv, pmid, pmax, pmin, peak: N.amax(pmax), low: N.amin(pmin), x: N.linspace(0, i.L, n + 1), tc: (2 * i.L) / a, E, e, p0v: i.p_sys - dpf };
}
function accu(i) {
  // charge slowly (isothermal) from pre-charge to p_max, discharge polytropically to p_min
  const p0 = i.pre_frac * i.p_min, VgMax = p0 / i.p_max, dV = VgMax * ((i.p_max / i.p_min) ** (1 / i.n_gas) - 1), V0 = i.V_demand / dV;
  const Vg1 = V0 * VgMax, E = Math.abs(i.n_gas - 1) < 1e-9 ? i.p_max * Vg1 * Math.log(i.p_max / i.p_min) : ((i.p_max * Vg1) / (i.n_gas - 1)) * (1 - (i.p_min / i.p_max) ** ((i.n_gas - 1) / i.n_gas));
  return { p0, V0, usable: dV, E };
}
const surge = {
  id: 'surge', title: 'Pressure surge (water hammer) and accumulator sizing', fidelity: 'numerical',
  summary: 'Pressure waves in a hydraulic line when a valve closes, by the method of characteristics with pipe-wall elasticity, compared with the Joukowsky estimate; plus gas-charged accumulator volume for a required fluid delivery.',
  equations: ['Continuity equation', 'Navier–Stokes equations', 'Fluid compressibility equations', 'Hydraulic capacitance equations', 'Darcy–Weisbach equation'],
  applicable: hasHyd,
  inputs: [
    sel('fluid', 'Hydraulic fluid', HYD, HYD[0], 'Line'), num('p_sys', 'Supply pressure', 'Pa', 20.7e6, 1e5, 60e6, 'Line'), num('p_back', 'Pressure downstream of the valve', 'Pa', 0.35e6, 0, 5e6, 'Line'),
    num('D', 'Line bore', 'm', 0.012, 0.001, 0.2, 'Line'), num('wall_mm', 'Tube wall thickness', 'mm', 0.9, 0.1, 20, 'Line'), sel('tube', 'Tube material', METAL_NAMES, 'Ti-6Al-4V', 'Line'),
    num('L', 'Line length', 'm', 12, 0.1, 500, 'Line'), num('v0', 'Initial flow velocity', 'm/s', 6, 0.01, 40, 'Line'),
    num('t_close', 'Valve closure time', 's', 0.01, 0, 10, 'Valve', 'Fast solenoid and servo-valves close in 5–30 ms'),
    num('fric', 'Line friction multiplier', '-', 1, 0, 20, 'Line', 'Scales the steady Darcy friction; unsteady (frequency-dependent) friction is not modelled'),
    num('V_demand', 'Fluid volume the accumulator must deliver', 'm³', 5e-4, 1e-7, 1, 'Accumulator', '1 litre = 1e-3 m³'), num('p_max', 'Accumulator maximum (system) pressure', 'Pa', 20.7e6, 1e5, 60e6, 'Accumulator'),
    num('p_min', 'Minimum useful pressure', 'Pa', 13.8e6, 5e4, 60e6, 'Accumulator'), num('pre_frac', 'Gas pre-charge / minimum pressure', '-', 0.9, 0.3, 1, 'Accumulator', '0.9 keeps the piston or bladder off its stop'),
    num('n_gas', 'Polytropic exponent on discharge', '-', 1.4, 1, 1.67, 'Accumulator', '1.0 slow (isothermal), 1.4 rapid (adiabatic nitrogen)'),
    num('n_periods', 'Wave periods (4L/a) simulated', '', 4, 1, 50, 'Numerics'), num('nReach', 'Line reaches', '', 40, 4, 800, 'Numerics', '', { step: 1, discrete: true }),
  ],
  // branch line carrying half the pump flow at 6 m/s; where the 3 mm minimum bore governs, the initial velocity is the lower value that this flow really gives
  defaults: (c) => { const p = c.systems.hyd_p_Pa, Q = Math.max(2e-6, 3e-8 * c.mass.mtow_kg), D = Math.max(0.003, Math.sqrt((4 * 0.5 * Q) / (Math.PI * 6))); return { p_sys: p, p_max: p, p_min: 0.67 * p, D, v0: (0.5 * Q) / ((Math.PI * D * D) / 4), wall_mm: Math.max(0.4, (1e3 * 4 * p * D) / (2 * METALS['Ti-6Al-4V'].Sy)), L: Math.max(1, 0.5 * (c.fuselage.len_m || 5)), V_demand: Math.max(1e-6, 0.12 * Q), fluid: c.mass.mtow_kg > 5700 && c.meta.type !== 'helicopter' ? HYD[0] : HYD[1] }; },
  run(i) {
    const r = moc(i), ac = accu(i), jk = r.rho * r.a * i.v0, slow = i.t_close > r.tc, est = slow ? (jk * r.tc) / i.t_close : jk, hoop = (r.peak * i.D) / (2 * r.e), Sy = (METALS[i.tube] || METALS['Ti-6Al-4V']).Sy, warnings = [];
    if (r.dpv0 <= 0) warnings.push('Line friction alone consumes the whole supply pressure at this velocity: reduce velocity or enlarge the bore. The surge result is not meaningful.');
    if (r.low < 0) warnings.push('Pressure falls below zero absolute in the rarefaction: the fluid would cavitate (column separation), which this single-phase model does not represent.');
    if (i.p_min >= i.p_max) warnings.push('Accumulator minimum pressure must be below the maximum pressure.');
    return {
      kpis: [
        kp('wave_speed_ms', 'Pressure-wave speed', r.a, 'm/s'), kp('t_critical_s', 'Critical closure time 2L/a', r.tc, 's', slow ? 'ok' : 'warn', slow ? 'Closure is slower than the wave return' : 'Closure is faster than the wave return: full Joukowsky surge'),
        kp('joukowsky_dp_Pa', 'Joukowsky surge ρ·a·Δv', jk, 'Pa'), kp('surge_estimate_Pa', 'Hand estimate for this closure time', est, 'Pa'),
        kp('surge_dp_Pa', 'Simulated surge above steady pressure at the valve', N.amax(r.pv) - r.p0v, 'Pa'),
        kp('peak_p_Pa', 'Peak line pressure', r.peak, 'Pa', r.peak < 1.35 * i.p_sys ? 'ok' : r.peak < 1.5 * i.p_sys ? 'warn' : 'bad', 'Transients are commonly limited to about 135% of system pressure; proof pressure is typically 150%'),
        kp('peak_ratio', 'Peak pressure / system pressure', r.peak / i.p_sys, '-'), kp('min_p_Pa', 'Minimum line pressure', r.low, 'Pa', r.low > 0 ? 'ok' : 'bad'),
        kp('hoop_stress_Pa', 'Tube hoop stress at peak pressure', hoop, 'Pa', hoop < Sy / 1.5 ? 'ok' : 'warn'), kp('friction_dp_Pa', 'Steady line friction loss', r.dpf, 'Pa'),
        kp('accu_volume_m3', 'Accumulator gas volume required', ac.V0, 'm³', undefined, `${(ac.V0 * 1e3).toFixed(2)} L`), kp('accu_precharge_Pa', 'Gas pre-charge pressure', ac.p0, 'Pa'),
        kp('accu_usable_frac', 'Usable fluid / gas volume', ac.usable, '-'), kp('accu_energy_J', 'Energy delivered by the accumulator', ac.E, 'J'),
      ],
      plots: [
        { type: 'line', title: 'Pressure history', xlabel: 'Time [ms]', ylabel: 'Pressure [MPa]', series: [{ name: 'At the valve', x: thin(r.t, 400).map((v) => v * 1e3), y: thin(r.pv, 400).map((v) => v / 1e6) }, { name: 'Mid-line', x: thin(r.t, 400).map((v) => v * 1e3), y: thin(r.pmid, 400).map((v) => v / 1e6) }], annotations: [{ y: (r.p0v + jk) / 1e6, label: 'Joukowsky' }] },
        { type: 'line', title: 'Pressure envelope along the line', xlabel: 'Distance from the source [m]', ylabel: 'Pressure [MPa]', series: [{ name: 'Maximum', x: thin(r.x), y: thin(r.pmax).map((v) => v / 1e6) }, { name: 'Minimum', x: thin(r.x), y: thin(r.pmin).map((v) => v / 1e6) }] },
      ],
      warnings, models: ['Method of characteristics (fixed grid, Courant number 1)', 'Korteweg wave speed with thin-wall pipe elasticity', 'Joukowsky surge relation', 'Polytropic gas accumulator'],
      assumptions: ['Single-phase liquid, constant bulk modulus (no entrained air)', 'Quasi-steady Darcy friction', 'Constant-pressure source upstream; valve area closes linearly with time', 'Accumulator charged isothermally, discharged polytropically, ideal gas', 'Default line velocity, closure time and the 135% / 150% transient limits are typical hydraulic-system practice, not sourced values'],
    };
  },
  convergence: { param: 'nReach', label: 'Line reaches', levels: [10, 20, 40, 80, 160], metric: 'surge_dp_Pa' },
  verify() {
    const b = { fluid: HYD[0], p_sys: 20e6, p_back: 0, D: 0.01, wall_mm: 1, tube: 'Steel 4340 (QT)', L: 10, v0: 5, t_close: 0, fric: 0, V_demand: 1e-3, p_max: 20e6, p_min: 10e6, pre_frac: 1, n_gas: 1, n_periods: 3, nReach: 20 }, r = moc(b), fl = FLUIDS[HYD[0]], ac = accu(b);
    return [
      N.check('Instant closure, no friction: surge = ρ·a·v0', N.amax(r.pv) - 20e6, fl.rho * r.a * 5, 1e-9, 'Joukowsky (1898)'),
      N.check('Wave speed with elastic wall', r.a, Math.sqrt(fl.bulk / fl.rho / (1 + (fl.bulk * 0.01) / (METALS['Steel 4340 (QT)'].E * 0.001))), 1e-12, 'Korteweg formula'),
      N.check('Wave returns as an equal rarefaction', 20e6 - N.amin(r.pv), fl.rho * r.a * 5, 1e-9, 'Reflection at a constant-pressure source'),
      N.check('Default initial velocity is consistent with the default bore and flow (small system, minimum bore governs)', ((d) => (d.v0 * Math.PI * d.D ** 2) / 4)(surge.defaults({ systems: { hyd_p_Pa: 6.9e6 }, mass: { mtow_kg: 1000 }, fuselage: { len_m: 8 }, meta: { type: 'aeroplane' } })), 0.5 * 3e-5, 1e-12, 'Continuity: v0·A = half the pump flow'),
      N.check('Default initial velocity is 6 m/s when the bore is not clamped', surge.defaults({ systems: { hyd_p_Pa: 20.7e6 }, mass: { mtow_kg: 78000 }, fuselage: { len_m: 37 }, meta: { type: 'aeroplane' } }).v0, 6, 1e-12, 'Line sizing rule'),
      N.check('Isothermal accumulator: ΔV = p0·V0·(1/p_min − 1/p_max)', b.V_demand, 10e6 * ac.V0 * (1 / 10e6 - 1 / 20e6), 1e-12, 'Boyle’s law'),
    ];
  },
  calibration: { params: [{ key: 'fric', min: 0, max: 20 }, { key: 'wall_mm', min: 0.1, max: 10 }], sweep: 't_close', target: 'surge_dp_Pa', note: 'Supply measured surge pressure against valve closure time from a line rig with a fast pressure transducer.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.peak_ratio > 1.35) out.push({ severity: o.peak_ratio > 1.5 ? 'critical' : 'warn', title: 'Surge pressure exceeds the usual transient allowance', detail: `Peak is ${(100 * o.peak_ratio).toFixed(0)}% of system pressure.`, action: `Slow the valve to well beyond 2L/a = ${(o.t_critical_s * 1e3).toFixed(1)} ms, reduce line velocity, or fit a surge damper/accumulator near the valve.`, basis: 'Typical hydraulic transient limit of 135% of system pressure' });
    if (o.accu_usable_frac < 0.15) out.push({ severity: 'advise', title: 'Accumulator is used inefficiently', detail: `Only ${(100 * o.accu_usable_frac).toFixed(0)}% of the gas volume is delivered as fluid.`, action: 'Widen the working pressure band or raise the pre-charge towards the minimum pressure to cut accumulator volume and mass.', basis: 'Polytropic gas law' });
    return out;
  },
};

// ---- 4. pneumatic duct, orifice and reservoir --------------------------------------------------
const g = GAMMA, PROC = ['Isothermal', 'Adiabatic (isentropic)'];
/** Fanno parameter f·Lmax/D (Darcy friction factor, length to the sonic point) and static pressure ratio p/p_sonic. */
export const fanno = (M) => (1 - M * M) / (g * M * M) + ((g + 1) / (2 * g)) * Math.log(((g + 1) * M * M) / (2 + (g - 1) * M * M));
export const fannoP = (M) => (1 / M) * Math.sqrt((g + 1) / (2 + (g - 1) * M * M));
const fannoInv = (F) => (F <= 0 ? 1 : F >= fanno(1e-7) ? 1e-7 : N.brent((M) => fanno(M) - F, 1e-7, 1, 1e-15));
/** Mass flow [kg/s] through an orifice of effective area CdA from stagnation (p0, T0) to back pressure pb. */
export function orificeFlow(CdA, p0, T0, pb) {
  if (pb >= p0) return 0;
  const rc = (2 / (g + 1)) ** (g / (g - 1)), r = Math.max(pb / p0, rc);
  return CdA * p0 * Math.sqrt(((2 * g) / ((g - 1) * R_AIR * T0)) * (r ** (2 / g) - r ** ((g + 1) / g)));
}
function ductFlow(i) {
  const A = (Math.PI * i.D ** 2) / 4, fLD = (i.f * i.L) / i.D + i.K_fit;
  if (i.p_back >= i.p0) return { A, fLD, M1: 0, M2: 0, mdot: 0, p1: i.p0, p2: i.p0, choked: false };
  const M1max = fannoInv(fLD), st = (M1) => { const M2 = fannoInv(fanno(M1) - fLD), p1 = i.p0 / (1 + 0.2 * M1 * M1) ** 3.5; return { M2, p1, p2: (p1 * fannoP(M2)) / fannoP(M1) }; };
  const cm = st(M1max * (1 - 1e-12)), choked = i.p_back <= cm.p2, lo = Math.min(1e-5, 0.5 * M1max), M1 = choked ? M1max : st(lo).p2 <= i.p_back ? lo : N.brent((M) => st(M).p2 - i.p_back, lo, M1max * (1 - 1e-12), 1e-14), s = st(M1), T1 = i.T0 / (1 + 0.2 * M1 * M1);
  return { A, fLD, M1, M2: choked ? 1 : s.M2, p1: s.p1, p2: s.p2, choked, T1, mdot: (s.p1 / (R_AIR * T1)) * M1 * Math.sqrt(g * R_AIR * T1) * A };
}
function blowdown(i) {
  const CdA = (i.Cd * Math.PI * i.d_or ** 2) / 4, m0 = (i.p0 * i.V_tank) / (R_AIR * i.T0), iso = i.process === PROC[0], nt = Math.max(4, Math.round(i.nSteps));
  const Tof = (m) => (iso ? i.T0 : i.T0 * (m / m0) ** (g - 1)), pof = (m) => (m * R_AIR * Tof(m)) / i.V_tank;
  const r = N.rk4((t, y) => [-orificeFlow(CdA, pof(y[0]), Tof(y[0]), i.p_amb)], 0, [m0], i.t_end, nt), m = r.y.map((y) => y[0]);
  return { t: r.t, m, p: m.map(pof), T: m.map(Tof), CdA, m0 };
}
const pneumatic = {
  id: 'pneumatic', title: 'Pneumatic duct, orifice and reservoir flow', fidelity: 'numerical',
  summary: 'Compressible air flow through a duct with friction (Fanno line) including choking, mass flow through an orifice or valve, and the pressure decay of a reservoir discharging through it.',
  equations: ['Compressible gas flow equations', 'Continuity equation', 'Orifice flow equations', 'Darcy–Weisbach equation'],
  inputs: [
    num('p0', 'Supply total pressure', 'Pa', 250e3, 2e4, 4e7, 'Supply', 'Engine bleed ≈ 200–350 kPa after the pressure regulator; storage bottles 10–20 MPa'), num('T0', 'Supply total temperature', 'K', 470, 150, 1200, 'Supply'),
    num('p_back', 'Duct delivery pressure at the consumer', 'Pa', 230e3, 500, 4e7, 'Supply', 'Static pressure held at the duct exit by the downstream valve or consumer'),
    num('p_amb', 'Ambient pressure for orifice and reservoir discharge', 'Pa', 80e3, 500, 4e7, 'Supply'),
    num('D', 'Duct bore', 'm', 0.08, 0.001, 1, 'Duct'), num('L', 'Duct length', 'm', 15, 0, 300, 'Duct'), num('f', 'Darcy friction factor', '-', 0.015, 0, 0.1, 'Duct', 'Smooth ducts at high Reynolds number ≈ 0.012–0.02'),
    num('K_fit', 'Bend and fitting losses (equivalent fL/D)', '-', 2, 0, 200, 'Duct'),
    num('d_or', 'Orifice or valve throat diameter', 'm', 0.02, 0.0002, 1, 'Orifice and reservoir'), num('Cd', 'Discharge coefficient', '-', 0.8, 0.3, 1, 'Orifice and reservoir', 'Sharp-edged ≈ 0.6–0.65 (rising towards 0.85 when choked); well-rounded nozzle ≈ 0.97'),
    num('V_tank', 'Reservoir volume', 'm³', 0.5, 1e-6, 5000, 'Orifice and reservoir'), sel('process', 'Reservoir gas process', PROC, PROC[1], 'Orifice and reservoir', 'Fast discharge is close to adiabatic; slow discharge with wall heat transfer approaches isothermal'),
    num('t_end', 'Discharge time simulated', 's', 20, 0.01, 1e5, 'Numerics'), num('nSteps', 'Time steps', '', 400, 10, 100000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => {
    const amb = isa(c.mission.cruise_alt_m || c.atm.alt_m || 0, c.atm.dISA_K || 0), pSite = (c.site?.p_hPa ?? 1013.25) * 100, Ts = (c.site?.T_C ?? 15) + 273.15, sc = Math.sqrt(Math.max(c.mass.mtow_kg, 1) / 78000);
    if (c.systems.bleed) return { p0: 250e3, T0: 470, p_back: 230e3, p_amb: Math.max(amb.p, isa(c.systems.cabin_alt_m || 0).p), D: N.clamp(0.09 * sc, 0.02, 0.2), L: Math.max(2, 0.4 * (c.fuselage.len_m || 5)), d_or: N.clamp(0.025 * sc, 0.004, 0.08), V_tank: Math.max(0.01, 0.4 * sc * sc), t_end: 30 };
    return { p0: 800e3, T0: Ts, p_back: 740e3, p_amb: pSite, D: N.clamp(0.008 * Math.sqrt(sc) * 3, 0.003, 0.02), L: Math.max(0.3, 0.3 * (c.fuselage.len_m || 1)), d_or: N.clamp(0.002 * Math.sqrt(sc) * 3, 0.0005, 0.006), V_tank: N.clamp(0.02 * sc * sc * 30, 2e-4, 0.05), t_end: 20 };
  },
  run(i) {
    const d = ductFlow(i), b = blowdown(i), mOr = orificeFlow(b.CdA, i.p0, i.T0, i.p_amb), rc = (2 / (g + 1)) ** (g / (g - 1)), last = b.t.length - 1, pEnd = b.p[last], warnings = [];
    const half = 0.5 * (i.p0 + i.p_amb), jh = b.p.findIndex((p) => p <= half), tHalf = jh > 0 ? b.t[jh - 1] + ((b.t[jh] - b.t[jh - 1]) * (b.p[jh - 1] - half)) / (b.p[jh - 1] - b.p[jh]) : NaN;
    if (i.p_back >= i.p0) warnings.push('Duct delivery pressure is not below the supply pressure: there is no duct flow.');
    if (d.choked) warnings.push('The duct is choked at its exit: mass flow no longer responds to lower back pressure, and the exit plane is sonic (noise, shock losses downstream).');
    if (d.M1 > 0.3 && !d.choked) warnings.push(`Duct inlet Mach number is ${d.M1.toFixed(2)}: bleed ducts are normally sized below about Mach 0.2–0.25 to limit losses and noise.`);
    if (Number.isNaN(tHalf)) warnings.push('The reservoir did not fall to the mean of supply and ambient pressure within the simulated time.');
    const pbs = N.linspace(0.05 * i.p0, i.p0, 50);
    return {
      kpis: [
        kp('duct_mdot_kgs', 'Duct mass flow', d.mdot, 'kg/s'), kp('duct_M_in', 'Duct inlet Mach number', d.M1, '-', d.M1 < 0.3 ? 'ok' : 'warn'), kp('duct_M_out', 'Duct exit Mach number', d.M2, '-'),
        kp('duct_choked', 'Duct choked (1 = yes)', d.choked ? 1 : 0, '-', d.choked ? 'warn' : 'ok'), kp('duct_dp_Pa', 'Total-pressure loss along the duct', d.mdot > 0 ? i.p0 - d.p2 * (1 + 0.2 * d.M2 ** 2) ** 3.5 : 0, 'Pa'),
        kp('duct_p_exit_Pa', 'Duct exit static pressure', d.p2, 'Pa'), kp('duct_fLD', 'Friction parameter fL/D incl. fittings', d.fLD, '-'),
        kp('orifice_mdot_kgs', 'Orifice mass flow at supply conditions', mOr, 'kg/s'), kp('orifice_choked', 'Orifice choked (1 = yes)', i.p_amb / i.p0 <= rc ? 1 : 0, '-', undefined, `Critical pressure ratio ${rc.toFixed(4)}`),
        kp('tank_t_half_s', 'Time for the reservoir to lose half its excess pressure', tHalf, 's'), kp('tank_p_end_Pa', 'Reservoir pressure at end', pEnd, 'Pa'), kp('tank_T_end_K', 'Reservoir gas temperature at end', b.T[last], 'K', b.T[last] > 233 ? 'ok' : 'warn'),
        kp('tank_mass_out_kg', 'Air mass discharged', b.m0 - b.m[last], 'kg'),
      ],
      plots: [
        { type: 'line', title: 'Reservoir pressure during discharge', xlabel: 'Time [s]', ylabel: 'Pressure [kPa]', series: [{ name: 'Reservoir pressure', x: thin(b.t), y: thin(b.p).map((v) => v / 1e3) }], annotations: [{ y: i.p_amb / 1e3, label: 'Ambient' }] },
        { type: 'line', title: 'Mass flow versus back pressure', xlabel: 'Back pressure / supply pressure [-]', ylabel: 'Mass flow [kg/s]', series: [{ name: 'Duct with friction', x: pbs.map((p) => p / i.p0), y: pbs.map((p) => ductFlow({ ...i, p_back: p }).mdot) }, { name: 'Orifice', x: pbs.map((p) => p / i.p0), y: pbs.map((p) => orificeFlow(b.CdA, i.p0, i.T0, p)) }], annotations: [{ x: i.p_back / i.p0, label: 'Duct delivery' }, { x: i.p_amb / i.p0, label: 'Orifice ambient' }] },
        { type: 'line', title: 'Reservoir gas temperature during discharge', xlabel: 'Time [s]', ylabel: 'Temperature [K]', series: [{ name: 'Gas temperature', x: thin(b.t), y: thin(b.T) }] },
      ],
      warnings, models: ['Fanno flow (adiabatic constant-area duct with friction) with subsonic inlet and choking', 'Isentropic nozzle from the supply plenum to the duct inlet', 'Compressible orifice equation with critical-ratio choking', 'Lumped reservoir mass balance, RK4'],
      assumptions: ['Perfect gas, γ = 1.4, constant friction factor', 'Adiabatic duct walls; fittings represented as equivalent fL/D', 'Reservoir gas spatially uniform; isothermal or isentropic limit', 'No condensation or icing at low discharge temperatures', 'Default supply conditions, duct and reservoir sizes are illustrative class-level values'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [25, 50, 100, 200, 400], metric: 'tank_p_end_Pa' },
  verify() {
    const b = { p0: 1e6, T0: 300, p_back: 1e5, p_amb: 1e5, D: 0.05, L: 10, f: 0.02, K_fit: 0, d_or: 0.01, Cd: 1, V_tank: 0.1, process: PROC[0], t_end: 8, nSteps: 400 };
    const bd = blowdown(b), CdA = (Math.PI * 1e-4) / 4, Kc = Math.sqrt(g) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1))), tau = 0.1 / (CdA * Kc * Math.sqrt(R_AIR * 300)), d = ductFlow(b);
    // round trip: with the solved inlet Mach the duct length to the exit state must equal fL/D
    return [
      N.check('Fanno parameter fL*/D at M = 0.5', fanno(0.5), 1.06906, 1e-5, 'Standard Fanno tables, γ = 1.4'),
      N.check('Fanno pressure ratio p/p* at M = 0.5', fannoP(0.5), 2.13809, 1e-5, 'Standard Fanno tables, γ = 1.4'),
      N.check('Choked orifice flow 0.04042·p0·A/√T0', orificeFlow(CdA, 1e6, 300, 1e5), (0.040416 * 1e6 * CdA) / Math.sqrt(300), 1e-4, 'Choked-flow constant for air'),
      N.check('Isothermal choked blow-down p = p0·exp(−t/τ)', bd.p[400], 1e6 * Math.exp(-8 / tau), 1e-6, 'Analytical solution while choked'),
      N.check('Duct solution satisfies the Fanno length', fanno(d.M1) - fanno(d.M2), d.fLD, 1e-8, 'Fanno relation'),
    ];
  },
  calibration: { params: [{ key: 'f', min: 0.005, max: 0.06 }, { key: 'Cd', min: 0.5, max: 1 }, { key: 'K_fit', min: 0, max: 50 }], sweep: 'p_back', target: 'duct_mdot_kgs', note: 'Supply measured duct mass flow against back pressure from a flow-bench test.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.duct_dp_Pa > 0.1 * i.p0) out.push({ severity: 'advise', title: 'High pressure loss in the duct', detail: `${(o.duct_dp_Pa / 1e3).toFixed(0)} kPa (${((100 * o.duct_dp_Pa) / i.p0).toFixed(0)}% of supply) is lost to friction.`, action: 'Increase the duct bore or shorten the run: lost bleed pressure has to be made up by extracting from a higher compressor stage, which costs fuel.', basis: 'Fanno-line total-pressure loss' });
    // the isentropic limit is reached by any full blow-down beyond a pressure ratio of about 2: it bounds the cooling, it is not a prediction
    if (o.tank_T_end_K < 233) out.push({ severity: 'advise', title: 'Cold gas during discharge (adiabatic bound)', detail: `Reservoir gas could fall to ${(o.tank_T_end_K - 273.15).toFixed(0)} °C if no heat reached it from the walls; real discharges lie between this and the isothermal case.`, action: 'Dry the gas to a low dew point and check seal and valve materials for low-temperature operation; re-run with the isothermal process to bracket the result.', basis: 'Isentropic expansion T/T0 = (p/p0)^((γ−1)/γ)' });
    return out;
  },
};

// ---- 5. gear transmission ---------------------------------------------------------------------
const FINAL = ['Planetary (sun in, carrier out, ring fixed)', 'Parallel-axis spur pair'];
/** Lewis form factor for 20° full-depth involute teeth (curve fit, generic). */
const lewisY = (z) => Math.max(0.2, 0.484 - 2.87 / z);
/** Hertz peak pressure for two parallel cylinders: load per length Fp, radii R1, R2 (R2 < 0 for an internal contact). */
export function hertzLine(Fp, R1, R2, E1, nu1, E2 = E1, nu2 = nu1) { const Es = 1 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2), Rr = 1 / (1 / R1 + 1 / R2); return Math.sqrt((Fp * Es) / (Math.PI * Rr)); }
function gearChain(i) {
  const m = METALS[i.mat] || METALS['Steel 4340 (QT)'], Rt = Math.max(i.rpm_in, i.rpm_out) / Math.min(i.rpm_in, i.rpm_out), reducer = i.rpm_in >= i.rpm_out, phi = N.rad(i.phi_deg), K = i.Ko * i.Kv * i.Km;
  const nSt = Math.max(1, Math.ceil(Math.log(Rt) / Math.log(i.r_max) - 1e-9)), plan = i.final === FINAL[0] && Rt >= 3, ratios = new Array(nSt).fill(Rt ** (1 / nSt));
  if (plan && ratios[0] < 3) { ratios.fill(nSt > 1 ? (Rt / 3) ** (1 / (nSt - 1)) : Rt); ratios[nSt - 1] = nSt > 1 ? 3 : Rt; }
  if (plan && ratios[nSt - 1] > 7) ratios[nSt - 1] = 7;
  let rpm = reducer ? i.rpm_in : i.rpm_out, P = i.P_W, eta = 1; const stages = [];
  // stages are listed from the high-speed end; torque is evaluated on the high-speed (pinion or sun) member
  ratios.forEach((r, s) => {
    const isPl = plan && s === nSt - 1, w = (rpm * 2 * Math.PI) / 60, T = P / Math.max(w, 1e-9), z1 = Math.round(i.z_pinion), np = isPl ? Math.round(i.n_planets) : 1; let z2, ratio, u, zr = 0;
    if (isPl) { const z0 = Math.round(z1 * (r - 1)), okz = (q) => (q - z1) % 2 === 0 && (q + z1) % np === 0 && q > z1 + 24; zr = z0 + ((z0 - z1) % 2 ? 1 : 0); for (let k = 0; k <= 2 * np + 2; k++) { if (okz(z0 + k)) { zr = z0 + k; break; } if (okz(z0 - k)) { zr = z0 - k; break; } } z2 = (zr - z1) / 2; ratio = 1 + zr / z1; u = z2 / z1; } else { z2 = Math.max(z1, Math.round(z1 * r)); ratio = z2 / z1; u = ratio; }
    const Y = lewisY(z1), mb = Math.cbrt((2 * T * K) / (i.psi * z1 * Y * np * i.sb_allow)), Es = m.E / (2 * (1 - m.nu ** 2));
    const mh = Math.cbrt((4 * T * K * Es * (1 + 1 / u)) / (Math.PI * np * i.psi * z1 * z1 * Math.cos(phi) * Math.sin(phi) * i.sh_allow ** 2));
    const mod = i.module_mm > 0 ? i.module_mm / 1e3 : Math.ceil(Math.max(mb, mh, 3e-4) * 4000) / 4000, b = i.psi * mod, d1 = mod * z1, Ft = (2 * T) / (d1 * np);
    const sb = (Ft * K) / (b * mod * Y), sh = hertzLine((Ft * K) / (b * Math.cos(phi)), (d1 / 2) * Math.sin(phi), ((mod * z2) / 2) * Math.sin(phi), m.E, m.nu);
    const rpmRel = isPl ? rpm * (1 - 1 / ratio) : rpm, v = ((rpmRel * 2 * Math.PI) / 60) * (d1 / 2), loss = isPl ? 2 * i.loss_mesh * (1 - 1 / ratio) + i.loss_brg : i.loss_mesh + i.loss_brg;
    stages.push({ type: isPl ? 'Planetary' : 'Spur pair', ratio, z1, z2, zr, np, mod, b, d1, T, Ft, sb, sh, v, fm: (z1 * rpmRel) / 60, rpm, eff: 1 - loss, assembles: isPl ? Number.isInteger((z1 + zr) / np) : true });
    eta *= 1 - loss; P *= 1 - loss; rpm /= ratio;
  });
  const ratio = stages.reduce((a, s) => a * s.ratio, 1);
  return { stages, eta, ratio, Rt, reducer, m, rpmLow: rpm, Pout: i.P_W * eta, Tout: (i.P_W * eta) / (((reducer ? rpm : stages[0].rpm) * 2 * Math.PI) / 60) };
}
const GB_CASE = (type, vehicle) => ({ meta: { type: vehicle }, prop: { type, n_eng: 2, P0_W: 1e6, rpm: 1200 }, rotor: { rpm: vehicle === 'helicopter' ? 300 : 0 }, systems: { gen_kVA: 30 } });
const gearbox = {
  id: 'gearbox', title: 'Gear transmission: ratios, tooth stresses and efficiency', fidelity: 'analytical',
  summary: 'Lays out a reduction (or accessory) gear train between engine and rotor, propeller or accessory speeds, sizes the teeth for bending and surface contact, and estimates efficiency, mesh frequencies and the heat carried away by the oil.',
  equations: ['Gear kinematic equations', 'Gear contact equations', 'Shaft torsion equations', 'Mechanical power balance equations'],
  inputs: [
    num('P_W', 'Transmitted power', 'W', 1.4e6, 1, 1e8, 'Duty'), num('rpm_in', 'Input (engine or motor) speed', 'rpm', 21000, 10, 200000, 'Duty'), num('rpm_out', 'Output speed', 'rpm', 258, 1, 200000, 'Duty'),
    sel('final', 'Final (low-speed) stage', FINAL, FINAL[0], 'Layout', 'A planetary stage is used only when the overall ratio is at least 3'), num('r_max', 'Maximum ratio per stage', '-', 4.5, 1.5, 8, 'Layout'),
    num('z_pinion', 'Pinion / sun teeth', '', 23, 12, 80, 'Layout', '', { step: 1, discrete: true }), num('n_planets', 'Planets', '', 4, 3, 8, 'Layout', '', { step: 1, discrete: true }),
    num('psi', 'Face width / module', '-', 12, 4, 30, 'Layout'), num('phi_deg', 'Pressure angle', 'deg', 20, 14.5, 30, 'Layout'), num('module_mm', 'Module (0 = size automatically)', 'mm', 0, 0, 50, 'Layout'),
    sel('mat', 'Gear steel (elastic properties)', METAL_NAMES, 'Steel 4340 (QT)', 'Material'),
    num('sb_allow', 'Allowable tooth-root bending stress', 'Pa', 380e6, 5e7, 1.5e9, 'Material', 'Case-carburised aerospace gear steels ≈ 350–450 MPa for long life (typical)'),
    num('sh_allow', 'Allowable contact stress', 'Pa', 1300e6, 2e8, 3e9, 'Material', 'Case-carburised steels ≈ 1200–1500 MPa (typical)'),
    num('Ko', 'Overload factor', '-', 1.25, 1, 3, 'Rating factors', 'AGMA-style application factor (user-supplied)'), num('Kv', 'Dynamic factor', '-', 1.2, 1, 2.5, 'Rating factors', 'Precision ground aerospace gears ≈ 1.1–1.3'), num('Km', 'Load-distribution factor', '-', 1.3, 1, 2.5, 'Rating factors', 'Includes planet load sharing'),
    num('loss_mesh', 'Power loss per mesh', '-', 0.006, 0.001, 0.05, 'Losses', 'Ground, oil-jet lubricated spur/helical meshes ≈ 0.4–1% (typical)'), num('loss_brg', 'Bearing and windage loss per stage', '-', 0.003, 0, 0.05, 'Losses'),
    num('dT_oil', 'Oil temperature rise across the gearbox', 'K', 30, 5, 80, 'Losses'),
  ],
  defaults: (c, up, d) => {
    const t = c.prop.type, heli = c.meta.type === 'helicopter'; let P = shaftPower(c, up), rin, rout = c.rotor.rpm || c.prop.rpm || 0;
    if (t === 'turboshaft') rin = 21000; else if (t === 'turboprop') rin = 20000; else if (t === 'turbofan' || t === 'turbojet') { rin = 15000; rout = 8000; P = Math.max(2e4, 1.5 * (c.systems.gen_kVA || 30) * 1000); } else rin = Math.max(c.prop.rpm || 2700, 1) * (t === 'electric' ? 3 : 2.2);
    if (!rout) rout = rin / 3;
    return { P_W: Math.max(P, 10), rpm_in: rin, rpm_out: rout, final: heli || t === 'turboprop' ? FINAL[0] : FINAL[1], n_planets: heli ? 5 : 4 };
  },
  run(i) {
    const r = gearChain(i), heat = i.P_W * (1 - r.eta), oil = heat / (OIL.rho * OIL.cp * i.dT_oil), ub = r.stages.map((s) => s.sb / i.sb_allow), uh = r.stages.map((s) => s.sh / i.sh_allow), vmax = N.amax(r.stages.map((s) => s.v)), warnings = [];
    if (!r.reducer) warnings.push('Output speed is above input speed: the train is treated as a speed increaser with the same gear geometry (pinion on the high-speed side).');
    if (Math.abs(r.ratio / r.Rt - 1) > 0.03) warnings.push(`Integer tooth counts give an overall ratio of ${r.ratio.toFixed(3)} against ${r.Rt.toFixed(3)} requested; adjust tooth numbers to trim.`);
    if (vmax > 120) warnings.push(`Pitch-line velocity reaches ${vmax.toFixed(0)} m/s: high-speed gearing above roughly 100–150 m/s needs special lubrication, profile accuracy and scuffing checks.`);
    if (r.stages.some((s) => !s.assembles)) warnings.push('The planetary stage does not satisfy the equal-spacing assembly condition (sun + ring teeth divisible by the number of planets); change the sun tooth count or planet number.');
    if (i.module_mm > 0 && (N.amax(ub) > 1 || N.amax(uh) > 1)) warnings.push('The specified module is too small for the allowable stresses.');
    return {
      kpis: [
        kp('gearbox_eff', 'Transmission efficiency', r.eta, '-'), kp('gear_ratio', 'Overall ratio achieved', r.ratio, '-'), kp('n_stages', 'Stages', r.stages.length, ''),
        kp('output_torque_Nm', r.reducer ? 'Output (low-speed) torque' : 'Output (high-speed) torque', r.Tout, 'N·m'), kp('heat_to_oil_W', 'Heat to the oil', heat, 'W'), kp('oil_flow_m3s', 'Oil flow for the temperature rise', oil, 'm³/s', undefined, `${(oil * 6e4).toFixed(2)} L/min`),
        kp('bending_stress_max_Pa', 'Highest tooth-root bending stress', N.amax(r.stages.map((s) => s.sb)), 'Pa', N.amax(ub) <= 1 ? 'ok' : 'bad'), kp('contact_stress_max_Pa', 'Highest contact (Hertz) stress', N.amax(r.stages.map((s) => s.sh)), 'Pa', N.amax(uh) <= 1 ? 'ok' : 'bad'),
        kp('pitch_velocity_ms', 'Highest pitch-line velocity', vmax, 'm/s', vmax < 120 ? 'ok' : 'warn'), kp('mesh_freq_Hz', 'First-stage mesh frequency', r.stages[0].fm, 'Hz'), kp('mesh_freq_last_Hz', 'Final-stage mesh frequency', r.stages[r.stages.length - 1].fm, 'Hz'),
        kp('module_max_mm', 'Largest module', 1e3 * N.amax(r.stages.map((s) => s.mod)), 'mm'),
      ],
      plots: [
        { type: 'bar', title: 'Tooth stress utilisation by stage', ylabel: 'Stress / allowable [-]', categories: r.stages.map((s, j) => `Stage ${j + 1} (${s.type})`), series: [{ name: 'Bending', y: ub }, { name: 'Contact', y: uh }] },
        { type: 'bar', title: 'Gear-mesh excitation frequencies', ylabel: 'Frequency [Hz]', categories: r.stages.map((s, j) => `Stage ${j + 1}`), series: [{ name: 'Mesh frequency', y: r.stages.map((s) => s.fm) }] },
      ],
      tables: [{ title: 'Stage summary (high-speed end first)', columns: ['Stage', 'Type', 'Ratio', 'Pinion/sun teeth', 'Gear/planet teeth', 'Ring teeth', 'Module [mm]', 'Face width [mm]', 'Input speed [rpm]', 'Input torque [N·m]', 'Tooth load [kN]', 'Bending [MPa]', 'Contact [MPa]', 'Pitch velocity [m/s]', 'Mesh freq. [Hz]'],
        rows: r.stages.map((s, j) => [j + 1, s.type, s.ratio, s.z1, s.z2, s.zr, s.mod * 1e3, s.b * 1e3, s.rpm, s.T, s.Ft / 1e3, s.sb / 1e6, s.sh / 1e6, s.v, s.fm]) }],
      outputs: { mesh_freqs_Hz: r.stages.map((s) => s.fm) },
      warnings, models: ['Kinematic ratio split with integer tooth counts', 'Lewis tooth-root bending with user rating factors', 'Hertz line contact at the pitch point', 'Constant fractional loss per mesh and per stage'],
      assumptions: ['Spur geometry, 20°-type full-depth teeth; helical overlap not credited', 'AGMA-style factors Ko, Kv, Km are user inputs, not computed from an AGMA/ISO rating', 'Allowable stresses, mesh and bearing loss fractions and default speeds are typical figures, not certified or sourced data', 'Equal load sharing between planets beyond the load-distribution factor; no scuffing or micropitting check'],
    };
  },
  verify() {
    const b = { P_W: 100e3, rpm_in: 6000, rpm_out: 1500, final: FINAL[0], r_max: 4.5, z_pinion: 24, n_planets: 4, psi: 10, phi_deg: 20, module_mm: 3, mat: 'Steel 4340 (QT)', sb_allow: 4e8, sh_allow: 1.3e9, Ko: 1, Kv: 1, Km: 1, loss_mesh: 0.01, loss_brg: 0, dT_oil: 30 };
    const r = gearChain(b), s = r.stages[0], T = 100e3 / ((6000 * 2 * Math.PI) / 60), Ft = (2 * T) / (0.072 * 4), sp = gearChain({ ...b, final: FINAL[1] }).stages[0];
    return [
      N.check('Planetary ratio 1 + Zr/Zs', s.ratio, 1 + 72 / 24, 1e-12, 'Willis equation, ring fixed'),
      N.check('Lewis bending stress Ft/(b·m·Y)', s.sb, Ft / (0.03 * 0.003 * (0.484 - 2.87 / 24)), 1e-10, 'Lewis (1892)'),
      N.check('Hertz line contact: p = 0.418·sqrt(F′E/R′) for ν = 0.3', hertzLine(1e5, 0.02, 0.03, 200e9, 0.3), 0.418 * Math.sqrt((1e5 * 200e9) / 0.012), 1e-3, 'Roark, cylinders in contact'),
      N.check('Twin turboprop: each gearbox carries half of the published all-engine shaft power', gearbox.defaults(GB_CASE('turboprop', 'aeroplane'), { propulsion: { P_shaft_W: 3e6 } }, {}).P_W, 1.5e6, 1e-12, 'Suite 7 publishes P_shaft_W for all engines together'),
      N.check('Helicopter: the main gearbox combines all engines (no double counting)', gearbox.defaults(GB_CASE('turboshaft', 'helicopter'), { propulsion: { P_shaft_W: 3e6 } }, {}).P_W, 3e6, 1e-12, 'Combining gearbox'),
      N.check('Without upstream data the helicopter gearbox takes rated power × engines', gearbox.defaults(GB_CASE('turboshaft', 'helicopter'), {}, {}).P_W, 2e6, 1e-12, 'Case fallback'),
      N.check('Speed increaser: output torque is on the high-speed shaft, T = η·P/ω_out', (gearChain({ ...b, final: FINAL[1], rpm_in: 1500, rpm_out: 6000 }).Tout * 6000 * 2 * Math.PI) / 60, 100e3 * 0.99, 1e-10, 'Energy conservation'),
      N.check('Spur pair power balance T_out·ω_out = η·P', (gearChain({ ...b, final: FINAL[1] }).Tout * (6000 / sp.ratio) * 2 * Math.PI) / 60, 100e3 * 0.99, 1e-10, 'Energy conservation'),
    ];
  },
  calibration: { params: [{ key: 'loss_mesh', min: 0.001, max: 0.03 }, { key: 'loss_brg', min: 0, max: 0.03 }], sweep: 'P_W', target: 'heat_to_oil_W', note: 'Supply measured oil heat rejection against transmitted power from a gearbox rig test.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    out.push({ severity: 'info', title: 'Transmission losses', detail: `${(100 * (1 - o.gearbox_eff)).toFixed(1)}% of power (${(o.heat_to_oil_W / 1e3).toFixed(1)} kW) becomes heat in the oil.`, action: 'Pass the heat load to the oil cooler sizing in Suite 12. Each 0.1% of mesh loss removed is fuel or battery energy saved for the whole flight; superfinished teeth and low-viscosity oils reduce it.', basis: 'Power balance' });
    if (o.contact_stress_max_Pa > i.sh_allow || o.bending_stress_max_Pa > i.sb_allow) out.push({ severity: 'warn', title: 'Tooth stress exceeds the allowable', detail: `Bending ${(o.bending_stress_max_Pa / 1e6).toFixed(0)} MPa, contact ${(o.contact_stress_max_Pa / 1e6).toFixed(0)} MPa.`, action: 'Increase the module or face width, add planets, or split the torque path.', basis: 'Lewis bending and Hertz contact' });
    out.push({ severity: 'advise', title: 'Avoid mesh-frequency resonances', detail: `Mesh excitation at ${o.mesh_freq_Hz.toFixed(0)} Hz (first stage) and ${o.mesh_freq_last_Hz.toFixed(0)} Hz (final stage).`, action: 'Check these tones and their harmonics against casing and shaft modes in Suite 10 and cabin noise in Suite 11.', basis: 'Gear-mesh excitation' });
    return out;
  },
};

// ---- 6. rolling bearings and shafts -------------------------------------------------------------
const BRG = ['Ball bearing (p = 3)', 'Roller bearing (p = 10/3)'], REL = ['90', '95', '96', '97', '98', '99'], A1 = { 90: 1, 95: 0.64, 96: 0.55, 97: 0.47, 98: 0.37, 99: 0.25 }; // ISO 281:2007 life modification factor for reliability
function brgShaft(i) {
  const pe = i.btype === BRG[0] ? 3 : 10 / 3, P = i.Fa / Math.max(i.Fr, 1e-12) > i.e_lim ? i.X * i.Fr + i.Y * i.Fa : i.Fr, L10 = (i.C / Math.max(P, 1e-12)) ** pe * 1e6 / (60 * i.rpm), a1 = A1[i.rel] ?? 1;
  const Creq = P * ((i.L_target_h * 60 * i.rpm) / (1e6 * a1 * i.a_iso)) ** (1 / pe), m = METALS[i.mat] || METALS['Steel 4340 (QT)'], w = (i.rpm * 2 * Math.PI) / 60, T = i.P_W / w, k = i.bore_ratio;
  const J = (Math.PI * i.D_o ** 4 * (1 - k ** 4)) / 32, tau = (T * i.D_o) / (2 * J), twist = (T * i.L_s) / (m.G * J), wc = (Math.PI / i.L_s) ** 2 * Math.sqrt((m.E * (J / 2)) / (m.rho * ((Math.PI * i.D_o ** 2 * (1 - k * k)) / 4)));
  return { pe, P, L10, a1, Lna: a1 * i.a_iso * L10, Creq, T, tau, twist, wc, Nc: (wc * 60) / (2 * Math.PI), m, dn: i.bore_mm * i.rpm, mass: m.rho * ((Math.PI * i.D_o ** 2 * (1 - k * k)) / 4) * i.L_s };
}
const bearing = {
  id: 'bearing', title: 'Rolling-bearing life and drive-shaft checks', fidelity: 'analytical',
  summary: 'Basic and reliability-adjusted rating life of a rolling bearing from its equivalent dynamic load, the speed (dn) check, and the torsional stress, twist and first bending critical speed of a drive shaft.',
  equations: ['Bearing load equations', 'Shaft torsion equations', 'Mechanical power balance equations'],
  inputs: [
    sel('btype', 'Bearing type', BRG, BRG[0], 'Bearing'), num('C', 'Basic dynamic load rating C', 'N', 60000, 10, 5e7, 'Bearing', 'From the bearing catalogue'),
    num('Fr', 'Radial load', 'N', 8000, 0, 5e7, 'Bearing'), num('Fa', 'Axial load', 'N', 1500, 0, 5e7, 'Bearing'),
    num('X', 'Radial factor X', '-', 0.56, 0, 2, 'Bearing', 'Deep-groove ball bearing typical: X = 0.56, Y ≈ 1.2–2.3 when Fa/Fr > e'), num('Y', 'Axial factor Y', '-', 1.6, 0, 10, 'Bearing'), num('e_lim', 'Threshold e on Fa/Fr', '-', 0.27, 0, 2, 'Bearing'),
    num('rpm', 'Shaft speed', 'rpm', 6000, 1, 300000, 'Duty'), num('bore_mm', 'Bearing bore', 'mm', 50, 2, 1000, 'Bearing'),
    sel('rel', 'Required reliability [%]', REL, '90', 'Life'), num('a_iso', 'Lubrication/cleanliness life factor', '-', 1, 0.05, 50, 'Life', 'a_ISO: below 1 for thin films or contamination, above 1 for clean, well-lubricated bearings'),
    num('L_target_h', 'Target life', 'h', 10000, 10, 5e5, 'Life', 'e.g. time between overhauls'), num('dn_limit', 'Speed limit dn', 'mm·rpm', 1.5e6, 1e5, 4e6, 'Life', 'Typical: grease ≈ 0.5·10⁶, oil-jet ≈ 1–1.5·10⁶, aero-engine main-shaft bearings up to ≈ 2.5–3·10⁶'),
    num('P_W', 'Shaft power', 'W', 500e3, 1, 1e8, 'Shaft'), num('D_o', 'Shaft outer diameter', 'm', 0.06, 0.002, 1, 'Shaft'), num('bore_ratio', 'Shaft inner / outer diameter', '-', 0.8, 0, 0.98, 'Shaft'),
    num('L_s', 'Shaft length between supports', 'm', 1.2, 0.02, 20, 'Shaft'), sel('mat', 'Shaft material', METAL_NAMES, 'Steel 4340 (QT)', 'Shaft'),
  ],
  defaults: (c, up, d) => {
    const heli = c.meta.type === 'helicopter', jet = c.prop.type === 'turbofan' || c.prop.type === 'turbojet', rpm = heli ? c.rotor.tr_rpm * 3.5 || 4000 : jet ? 8000 : c.prop.rpm || c.rotor.rpm || 3000;
    const P = Math.max(5, heli ? 0.1 * d.P_total : jet ? 1.5 * (c.systems.gen_kVA || 30) * 1000 : shaftPower(c, up)), T = P / ((rpm * 2 * Math.PI) / 60);
    // shaft sized for a torsional shear of about 25% of shear yield; bearing loads from a tooth load at an assumed gear radius (estimates)
    const mat = heli ? 'Al 7075-T6' : 'Steel 4340 (QT)', ms = METALS[mat], k = heli ? 0.94 : 0.8, Ds = Math.max(0.004, Math.cbrt((16 * T) / (Math.PI * (1 - k ** 4) * 0.25 * 0.577 * ms.Sy))), L = heli ? N.clamp((c.rotor.tr_arm_m || 6) / 6.5, 0.4, 1.6) : Math.max(0.05, 12 * Ds);
    // diameter is the larger of the torsion size and the size that keeps the first critical speed 35% above running speed
    const Do = Math.max(Ds, (1.35 * ((rpm * 2 * Math.PI) / 60) * L * L) / (Math.PI ** 2 * Math.sqrt((ms.E * (1 + k * k)) / (16 * ms.rho)))), rg = Math.max(0.01, 1.5 * Ds), Fr = Math.max(1, (1.1 * T) / rg), life = c.econ.util_fh_yr > 2000 ? 10000 : 3000, Cq = Fr * ((life * 60 * rpm) / 1e6) ** (1 / 3);
    return { rpm, P_W: P, D_o: Do, bore_ratio: k, L_s: L, mat, Fr, Fa: 0.15 * Fr, C: 1.25 * Cq, bore_mm: Math.max(4, 1e3 * Do), L_target_h: life };
  },
  run(i) {
    const r = brgShaft(i), lifeM = r.Lna / i.L_target_h, ty = 0.577 * r.m.Sy, ms = ty / Math.max(r.tau, 1e-9) - 1, sr = i.rpm / r.Nc, warnings = [];
    if (r.dn > i.dn_limit) warnings.push(`Speed parameter dn = ${r.dn.toExponential(2)} mm·rpm exceeds the stated limit: cage, lubrication and centrifugal ball loads need a high-speed bearing design.`);
    if (sr > 0.75 && sr < 1.3) warnings.push('The shaft runs close to its first bending critical speed (within 25–30%): move the critical speed or add a support/damper.');
    if (sr >= 1.3) warnings.push('The shaft is supercritical: it must pass through the critical speed on every start and needs damping and good balance.');
    if (r.P > 0.5 * i.C) warnings.push('Equivalent load exceeds half the dynamic rating: the life equation is outside its normal range and static capacity may govern.');
    const Ps = N.logspace(Math.max(r.P / 5, 1e-3), r.P * 3, 40);
    return {
      kpis: [
        kp('bearing_L10_h', 'Basic rating life L10', r.L10, 'h'), kp('bearing_Lna_h', `Adjusted life at ${i.rel}% reliability`, r.Lna, 'h', lifeM >= 1 ? 'ok' : 'bad', `Target ${i.L_target_h.toFixed(0)} h`),
        kp('life_margin', 'Adjusted life / target life', lifeM, '-', lifeM >= 1.2 ? 'ok' : lifeM >= 1 ? 'warn' : 'bad'), kp('P_equiv_N', 'Equivalent dynamic load', r.P, 'N'), kp('C_required_N', 'Dynamic rating needed for the target life', r.Creq, 'N'),
        kp('a1', 'Reliability factor a1', r.a1, '-'), kp('dn_mm_rpm', 'Speed parameter dn', r.dn, 'mm·rpm', r.dn <= i.dn_limit ? 'ok' : 'warn'),
        kp('shaft_torque_Nm', 'Shaft torque', r.T, 'N·m'), kp('shaft_shear_Pa', 'Torsional shear stress', r.tau, 'Pa', ms > 1 ? 'ok' : ms > 0 ? 'warn' : 'bad'), kp('shaft_ms_yield', 'Margin on shear yield (0.577·Sy)', ms, '-', ms > 1 ? 'ok' : ms > 0 ? 'warn' : 'bad'),
        kp('shaft_twist_deg', 'Shaft twist', N.deg(r.twist), 'deg'), kp('shaft_crit_rpm', 'First bending critical speed (simply supported)', r.Nc, 'rpm'), kp('speed_ratio', 'Running speed / critical speed', sr, '-', sr < 0.75 ? 'ok' : 'warn'),
        kp('shaft_mass_kg', 'Shaft mass', r.mass, 'kg'),
      ],
      plots: [
        { type: 'line', title: 'Bearing life versus equivalent load', xlabel: 'Equivalent load [kN]', ylabel: 'Life [h]', xlog: true, ylog: true, series: [{ name: 'L10 (90%)', x: Ps.map((p) => p / 1e3), y: Ps.map((p) => ((i.C / p) ** r.pe * 1e6) / (60 * i.rpm)) }, { name: `${i.rel}% reliability`, x: Ps.map((p) => p / 1e3), y: Ps.map((p) => (r.a1 * i.a_iso * (i.C / p) ** r.pe * 1e6) / (60 * i.rpm)) }], annotations: [{ y: i.L_target_h, label: 'Target' }, { x: r.P / 1e3, label: 'Load' }] },
        { type: 'bar', title: 'Life factor against required reliability', ylabel: 'a1 [-]', categories: REL.map((v) => v + '%'), series: [{ name: 'a1 (ISO 281)', y: REL.map((v) => A1[v]) }] },
      ],
      warnings, models: ['ISO 281 basic rating life L10 = (C/P)^p', 'Reliability factor a1 (ISO 281:2007 table)', 'Elementary torsion of a hollow circular shaft', 'Euler–Bernoulli first bending critical speed, simply supported'],
      assumptions: ['Constant load and speed (use a cubic-mean load for a duty spectrum)', 'X, Y and e are catalogue inputs', 'Rigid bearings and uniform shaft; gyroscopic effects, coupling masses and support flexibility ignored (Suite 10 covers rotor dynamics)', 'Typical handbook material properties; default loads, rating, dn limit and shaft size are sizing-rule estimates'],
    };
  },
  verify() {
    const b = { btype: BRG[0], C: 20000, Fr: 10000, Fa: 0, X: 0.56, Y: 1.6, e_lim: 0.27, rpm: 1000, bore_mm: 50, rel: '90', a_iso: 1, L_target_h: 1000, dn_limit: 1e6, P_W: 1e5, D_o: 0.05, bore_ratio: 0, L_s: 1, mat: 'Steel 4340 (QT)' };
    const r = brgShaft(b), rr = brgShaft({ ...b, btype: BRG[1] }), T = 1e5 / ((1000 * 2 * Math.PI) / 60);
    return [
      N.check('Ball bearing C/P = 2: L10 = 8 million revolutions', r.L10, 8e6 / 60000, 1e-12, 'ISO 281'),
      N.check('Roller bearing exponent 10/3', rr.L10, (2 ** (10 / 3) * 1e6) / 60000, 1e-12, 'ISO 281'),
      N.check('Solid shaft shear 16T/(πd³)', r.tau, (16 * T) / (Math.PI * 0.05 ** 3), 1e-12, 'Elementary torsion'),
      N.check('Critical speed (π/L)²·(d/4)·sqrt(E/ρ) for a solid shaft', r.wc, Math.PI ** 2 * (0.05 / 4) * Math.sqrt(METALS['Steel 4340 (QT)'].E / METALS['Steel 4340 (QT)'].rho), 1e-12, 'Pinned–pinned Euler–Bernoulli beam'),
      N.check('Required rating reproduces the target life', brgShaft({ ...b, C: r.Creq }).Lna, 1000, 1e-10, 'Inverse of the life equation'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.life_margin < 1) out.push({ severity: 'warn', title: 'Bearing life is below target', detail: `${o.bearing_Lna_h.toFixed(0)} h at ${i.rel}% reliability against ${i.L_target_h.toFixed(0)} h.`, action: `Select a bearing with C ≥ ${(o.C_required_N / 1e3).toFixed(1)} kN, reduce the load, or improve lubrication cleanliness (a_ISO).`, basis: 'ISO 281 rating life' });
    else if (o.life_margin > 5) out.push({ severity: 'info', title: 'Bearing is generously sized', detail: `Adjusted life is ${o.life_margin.toFixed(1)} times the target.`, action: 'A smaller bearing would save mass; alternatively extend the overhaul interval and feed the life to Suite 22 (reliability) and Suite 26 (maintenance cost).', basis: 'Life margin' });
    if (o.speed_ratio > 0.75) out.push({ severity: 'warn', title: 'Shaft critical-speed margin is insufficient', detail: `Running at ${(100 * o.speed_ratio).toFixed(0)}% of the first critical speed.`, action: 'Shorten the span with an intermediate bearing, increase diameter, use a higher specific-stiffness material (composite tube), or design as supercritical with dampers; confirm in Suite 10.', basis: 'Critical-speed separation of about 25%' });
    return out;
  },
};

// ---- 7. lubrication film ----------------------------------------------------------------------
/** 1-D Reynolds equation for a fixed-incline pad: d/dx(h³ dp/dx) = 6μU dh/dx, p = 0 at both ends. */
function slider(i) {
  const n = Math.max(4, Math.round(i.nNodes)), dx = i.B / n, h1 = i.h2 * i.h_ratio, h = (x) => h1 - ((h1 - i.h2) * x) / i.B, x = N.linspace(0, i.B, n + 1);
  const a = [], b = [], c = [], d = [];
  for (let j = 1; j < n; j++) { const hw = h(x[j] - dx / 2), he = h(x[j] + dx / 2); a.push(hw ** 3); c.push(he ** 3); b.push(-(hw ** 3 + he ** 3)); d.push(6 * i.mu * i.U * (he - hw) * dx); }
  const p = [0, ...N.solveTridiag(a, b, c, d), 0], Wl = N.trapz(x, p), xcp = N.trapz(x, p.map((v, j) => v * x[j])) / Math.max(Wl, 1e-300);
  // flow per unit width at cell faces (must be constant) and shear on the runner
  const q = [], tauI = []; let Fl = 0;
  for (let j = 0; j < n; j++) { const hm = h(x[j] + dx / 2), dp = (p[j + 1] - p[j]) / dx; q.push((i.U * hm) / 2 - (hm ** 3 * dp) / (12 * i.mu)); Fl += ((i.mu * i.U) / hm + (hm / 2) * dp) * dx; tauI.push(dp); }
  const K = i.h_ratio - 1, WlA = K > 1e-9 ? ((6 * i.mu * i.U * i.B ** 2) / (i.h2 ** 2 * K * K)) * (Math.log(1 + K) - (2 * K) / (K + 2)) : 0;
  return { x, p, Wl, WlA, xcp, q, Fl, h: x.map(h), pmax: N.amax(p), qm: N.mean(q), qdev: (N.amax(q) - N.amin(q)) / Math.abs(N.mean(q)) };
}
const lubrication = {
  id: 'lubrication', title: 'Hydrodynamic lubrication film (tilting-pad slider)', fidelity: 'numerical',
  summary: 'Pressure in the converging oil film of a thrust or slider pad from the Reynolds equation: load capacity, film thickness for a given load, friction, power loss and the lubrication regime.',
  equations: ['Navier–Stokes equations', 'Continuity equation', 'Bearing load equations', 'Mechanical power balance equations'],
  inputs: [
    num('B', 'Pad length in the sliding direction', 'm', 0.03, 0.001, 1, 'Geometry'), num('Wd', 'Pad width', 'm', 0.09, 0.001, 2, 'Geometry'),
    num('h2', 'Outlet (minimum) film thickness', 'm', 20e-6, 1e-7, 1e-3, 'Geometry'), num('h_ratio', 'Inlet / outlet film thickness', '-', 2.2, 1.01, 10, 'Geometry', 'Load capacity peaks near 2.2'),
    num('U', 'Sliding speed', 'm/s', 20, 0.01, 300, 'Operation'), num('mu', 'Oil dynamic viscosity at film temperature', 'Pa·s', 0.012, 1e-4, 5, 'Operation', 'MIL-PRF-23699 ≈ 0.025 at 40 °C, ≈ 0.005 at 100 °C (typical)'),
    num('load', 'Applied load on the pad', 'N', 5000, 0, 1e8, 'Operation'), num('Rq_um', 'Combined surface roughness Rq', 'µm', 0.4, 0.01, 20, 'Operation', 'Root-sum-square of both surfaces'),
    num('nNodes', 'Grid intervals', '', 100, 8, 5000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => ({ U: N.clamp(((c.rotor.rpm || c.prop.rpm || 3000) * 2 * Math.PI) / 60 * 0.05, 2, 120) }),
  run(i) {
    const r = slider(i), ar = i.Wd / i.B, Wcap = r.Wl * i.Wd, F = r.Fl * i.Wd, Ploss = F * i.U;
    const hReq = i.load > 0 ? i.h2 * Math.sqrt(Wcap / i.load) : Infinity, lam = Math.min(hReq, 1) / (i.Rq_um * 1e-6), Q = r.qm * i.Wd, dT = Q > 0 ? Ploss / (OIL.rho * OIL.cp * Q) : 0, warnings = [];
    if (ar < 4) warnings.push(`Width/length = ${ar.toFixed(1)}: side leakage, ignored by the 1-D film equation, reduces the real load capacity below the value shown. Use a 2-D Reynolds solution for design.`);
    if (lam < 3) warnings.push(`Film parameter λ = ${lam.toFixed(1)} is below 3: asperity contact (mixed lubrication) and wear are likely at the applied load.`);
    if (dT > 40) warnings.push('Adiabatic oil temperature rise exceeds 40 K: viscosity at the film temperature will be lower than assumed, reducing the film thickness.');
    return {
      kpis: [
        kp('load_capacity_N', 'Load capacity at the stated film thickness', Wcap, 'N', Wcap >= i.load ? 'ok' : 'warn'), kp('p_max_Pa', 'Peak film pressure', r.pmax * (i.load > 0 ? i.load / Wcap : 1), 'Pa', undefined, 'At the applied load'),
        kp('h_min_m', 'Minimum film thickness at the applied load', Number.isFinite(hReq) ? hReq : i.h2, 'm', lam >= 3 ? 'ok' : lam >= 1 ? 'warn' : 'bad'), kp('lambda_ratio', 'Film parameter λ = h_min/Rq', Number.isFinite(lam) ? lam : 1e6, '-', lam >= 3 ? 'ok' : lam >= 1 ? 'warn' : 'bad', 'λ > 3 full film; 1–3 mixed; < 1 boundary'),
        kp('friction_coeff', 'Friction coefficient at the stated film', F / Math.max(Wcap, 1e-12), '-'), kp('power_loss_W', 'Viscous power loss at the stated film', Ploss, 'W'), kp('oil_flow_m3s', 'Oil flow through the film', Q, 'm³/s'),
        kp('oil_dT_K', 'Adiabatic oil temperature rise', dT, 'K', dT < 40 ? 'ok' : 'warn'), kp('x_cp', 'Centre of pressure from the inlet / pad length', r.xcp / i.B, '-', undefined, 'Pivot location for a tilting pad'),
        kp('load_per_width_Nm', 'Load capacity per unit width', r.Wl, 'N/m'), kp('flow_uniformity', 'Flow non-uniformity along the film', r.qdev, '-', r.qdev < 1e-3 ? 'ok' : 'warn'),
      ],
      plots: [
        { type: 'line', title: 'Film pressure', xlabel: 'Distance from the inlet [mm]', ylabel: 'Pressure [MPa]', series: [{ name: 'Pressure', x: thin(r.x).map((v) => v * 1e3), y: thin(r.p).map((v) => v / 1e6) }], annotations: [{ x: r.xcp * 1e3, label: 'Centre of pressure' }] },
        { type: 'line', title: 'Film thickness', xlabel: 'Distance from the inlet [mm]', ylabel: 'Film thickness [µm]', series: [{ name: 'Film', x: thin(r.x).map((v) => v * 1e3), y: thin(r.h).map((v) => v * 1e6) }] },
      ],
      warnings, models: ['1-D Reynolds equation, conservative finite differences with a tridiagonal solve', 'Couette + Poiseuille shear on the runner', 'Film-thickness scaling W ∝ 1/h² at fixed inclination ratio'],
      assumptions: ['Infinitely wide pad (no side leakage), rigid surfaces', 'Isoviscous, incompressible Newtonian oil', 'Fixed inlet/outlet film ratio (a tilting pad adjusts it automatically)', 'No cavitation (converging film only)', 'Default pad size, load, viscosity and roughness are illustrative typical values'],
    };
  },
  convergence: { param: 'nNodes', label: 'Grid intervals', levels: [10, 20, 40, 80, 160], metric: 'load_per_width_Nm' },
  verify() {
    const r = slider({ B: 0.05, Wd: 1, h2: 2e-5, h_ratio: 2.2, U: 10, mu: 0.02, load: 0, Rq_um: 0.4, nNodes: 400 });
    return [
      N.check('Plane-slider load capacity (analytical)', r.Wl, r.WlA, 1e-4, 'Reynolds (1886): W = 6μUB²/(h2²K²)·[ln(1+K) − 2K/(K+2)]'),
      N.check('Maximum load coefficient 0.1602 near K = 1.2', (r.WlA * (2e-5) ** 2) / (0.02 * 10 * 0.05 ** 2), 0.1602, 1e-3, 'Cameron, Basic Lubrication Theory'),
      N.check('Mass conservation along the film', r.qdev, 0, 1e-9, 'Volume flow per unit width is constant'),
    ];
  },
  calibration: { params: [{ key: 'mu', min: 1e-3, max: 0.5 }, { key: 'h_ratio', min: 1.1, max: 5 }], sweep: 'U', target: 'power_loss_W', note: 'Supply measured bearing friction power against speed from a bearing rig to fit the effective film viscosity.' },
  recommend(res) {
    const o = res.outputs, out = [];
    if (o.lambda_ratio < 3) out.push({ severity: o.lambda_ratio < 1 ? 'critical' : 'warn', title: 'Oil film too thin for the load', detail: `λ = ${o.lambda_ratio.toFixed(1)}; minimum film ${(o.h_min_m * 1e6).toFixed(1)} µm.`, action: 'Increase pad area or speed, use a higher-viscosity oil, lower the oil inlet temperature or improve the surface finish.', basis: 'Film parameter λ' });
    else out.push({ severity: 'info', title: 'Full-film lubrication', detail: `λ = ${Math.min(o.lambda_ratio, 999).toFixed(0)}: the surfaces are fully separated and wear is negligible; the price is ${o.power_loss_W.toFixed(0)} W of viscous loss per pad.`, action: 'If the margin is large, a lower-viscosity oil cuts friction loss and oil-cooler load.', basis: 'Reynolds film' });
    return out;
  },
};

export default {
  id: 'hydmech', n: 18,
  tagline: 'Whether pumps, lines, actuators, ducts, gears, bearings and shafts deliver the force, speed and life required — and what power and heat that costs.',
  analyses: [network, actuator, surge, pneumatic, gearbox, bearing, lubrication],
  consumes: [{ from: 'propulsion', keys: ['P_shaft_W'], why: 'Shaft power of all engines, shared between the gearboxes and drive shafts' }],
  provides: [
    { key: 'pump_power_W', label: 'Hydraulic pump power', unit: 'W' }, { key: 'actuator_force_N', label: 'Actuator stall force', unit: 'N' }, { key: 'actuator_rate_ms', label: 'Actuator no-load rate', unit: 'm/s' },
    { key: 'line_dp_Pa', label: 'Hydraulic line pressure loss', unit: 'Pa' }, { key: 'gearbox_eff', label: 'Gearbox efficiency', unit: '-' }, { key: 'bearing_L10_h', label: 'Bearing L10 life', unit: 'h' },
  ],
  handoff: [
    { model: '3-D Navier–Stokes flow in pumps, valves and manifolds', why: 'Component internal flow needs a resolved 3-D mesh; valves and pumps are represented here by characteristic curves and discharge coefficients', tool: 'CFD (finite volume) with moving-mesh or cavitation models' },
    { model: 'Full hydraulic/pneumatic system differential–algebraic network with all consumers, logic and failure cases', why: 'The native network has one pump and three lumped consumers; a complete aircraft circuit needs a component-library system simulator', tool: '1-D multi-domain system simulation tools' },
    { model: 'Flexible-shaft whirl, gear-mesh stiffness dynamics and gearbox–rotor coupling', why: 'Needs rotor-dynamic models with gyroscopics and support flexibility', tool: 'Suite 10 (vibration and rotor dynamics); dedicated drivetrain multibody tools' },
    { model: 'AGMA/ISO 6336 gear rating, tooth contact analysis, scuffing and micropitting', why: 'Standard ratings need detailed tooth geometry, accuracy grades and material certificates; only Lewis bending and pitch-point Hertz contact are computed', tool: 'Gear design and rating software' },
    { model: 'Elastohydrodynamic and thermal lubrication of gear and rolling contacts; 2-D journal and finite-width pad bearings', why: 'Requires coupled elastic deformation, piezo-viscous and thermal film solution', tool: 'EHL / thermo-hydrodynamic bearing codes' },
    { model: 'Fluid–structure coupling of hydraulic lines (pump ripple, clamp loads) and column separation', why: 'The surge model is single-phase with rigid supports', tool: 'Line-dynamics (transfer-matrix / MOC with FSI) codes' },
    { model: 'Bearing and gearbox temperature fields', why: 'Only bulk oil heat load is estimated', tool: 'Thermal network or CHT analysis (Suite 12 for network level)' },
  ],
};
