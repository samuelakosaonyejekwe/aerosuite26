// Suite 7 — Propulsion and Engine Performance.
// Gas-turbine parametric cycle with variable gas properties, simplified matched off-design maps, generic
// compressor map with surge margin, combustor energy balance and emissions, piston fuel-air cycle,
// spool acceleration transient, and electric / hybrid / fuel-cell propulsion chains.
//
// Published coupling keys (one coherent meaning for every vehicle class):
//   thrust_N, P_shaft_W, fuel_flow_kgs, co2_kg_s, heat_rejection_W  -> all engines together, at the flight point
//   thrust_static_N                                               -> per engine, sea-level static
//   tsfc_kg_Ns                                                    -> flight point (equivalent V·BSFC/η_prop for shaft engines, 0 for electric)
//   T4_K                                                          -> hot-section temperature (turbine inlet / peak cylinder gas / motor winding)

import * as N from '../core/numerics.js';
import { isa, T0 as T_SL, P0 as P_SL, R_AIR } from '../core/atmosphere.js';
import { FLUIDS } from '../data/materials.js';

// ---- gas properties ---------------------------------------------------------------------------
// cp(T) [kJ/kg/K] as polynomials in z = T/1000. Dry air: Walsh & Fletcher, Gas Turbine Performance, ch. 3 (200–2000 K).
// CO2, H2O, O2: quartics through tabulated ideal-gas values at 300/600/1000/1500/2000 K. Combustion products are an exact
// ideal-gas mixture: air + fuel·(CO2 + H2O formed − O2 consumed), so any fuel composition (kerosene, SAF, hydrogen) is covered.
const AIR = [0.992313, 0.236688, -1.852148, 6.083152, -8.893933, 7.097112, -3.234725, 0.794571, -0.081873];
const TK = [0.3, 0.6, 1.0, 1.5, 2.0];
const fitCp = (v) => [...N.solve(TK.map((t) => [1, t, t * t, t ** 3, t ** 4]), v), 0, 0, 0, 0];
const CP_CO2 = fitCp([0.846, 1.075, 1.234, 1.326, 1.371]), CP_H2O = fitCp([1.864, 2.016, 2.288, 2.609, 2.836]), CP_O2 = fitCp([0.918, 1.003, 1.09, 1.143, 1.181]);
const R_CO2 = 188.92, R_H2O = 461.52, R_O2 = 259.84, Y_O2 = 0.2314, TREF = 298.15, TMAX = 2200, TMIN = 150, FARADAY = 96485.332;
const FUELS = ['Jet A-1', 'SAF (HEFA-SPK)', 'Liquid hydrogen', 'Avgas 100LL'];

/** Fuel stoichiometry from its CO2 factor: carbon mass fraction -> CO2, the rest is hydrogen -> H2O. All per kg fuel. */
function fuelOf(name) {
  const F = FLUIDS[name] || FLUIDS['Jet A-1'], mC = (F.co2_per_kg * 12.011) / 44.009, mH = 1 - mC;
  const co2 = F.co2_per_kg, h2o = mH * 8.9365, o2 = (co2 * 31.998) / 44.009 + mH * 7.9365;
  const gD = { c: AIR.map((_, k) => co2 * CP_CO2[k] + h2o * CP_H2O[k] - o2 * CP_O2[k]), R: co2 * R_CO2 + h2o * R_H2O - o2 * R_O2 };
  return { name, LHV: F.LHV, co2, h2o, o2, fst: Y_O2 / o2, lc: F.lifecycle_factor ?? 1, gD };
}
/** Gas mixture for fuel-air ratio f (kg fuel burned per kg air). ideal = calorically perfect air (cp 1004.7, γ 1.4). */
function gas(f = 0, fu = null, ideal = false) {
  if (ideal) return { c: [1.0047], R: R_AIR };
  if (!f || !fu) return { c: AIR, R: R_AIR };
  return { c: AIR.map((a, k) => (a + f * fu.gD.c[k]) / (1 + f)), R: (R_AIR + f * fu.gD.R) / (1 + f) };
}
const poly = (c, z, w) => { let s = 0; for (let k = c.length - 1; k >= 0; k--) s = s * z + c[k] / w(k); return s; };
const cp = (g, T) => 1000 * poly(g.c, N.clamp(T, TMIN, TMAX) / 1000, () => 1);
const hRaw = (g, T) => 1e6 * (T / 1000) * poly(g.c, T / 1000, (k) => k + 1);
/** Sensible enthalpy h = ∫cp dT [J/kg] (cp frozen above TMAX where the fits end). */
const h = (g, T) => (T <= TMAX ? hRaw(g, T) : hRaw(g, TMAX) + cp(g, TMAX) * (T - TMAX));
const phiRaw = (g, T) => { const z = T / 1000; let s = 0; for (let k = g.c.length - 1; k >= 1; k--) s = s * z + g.c[k] / k; return 1000 * (g.c[0] * Math.log(z) + s * z); };
/** Entropy function φ = ∫cp dT/T [J/kg/K]; s2 − s1 = φ2 − φ1 − R ln(p2/p1). */
const phi = (g, T) => (T <= TMAX ? phiRaw(g, T) : phiRaw(g, TMAX) + cp(g, TMAX) * Math.log(T / TMAX));
const gam = (g, T) => { const c = cp(g, T); return c / (c - g.R); };
const Tofh = (g, hh, T) => { for (let k = 0; k < 40; k++) { const d = (h(g, T) - hh) / cp(g, T); T = Math.max(100, T - d); if (Math.abs(d) < 1e-7) break; } return T; };
const Tofphi = (g, ph, T) => { for (let k = 0; k < 40; k++) { const d = (ph - phi(g, T)) / cp(g, T); T *= Math.exp(N.clamp(d, -1, 1)); if (Math.abs(d) < 1e-11) break; } return T; };
/** Burner fuel-air ratio from the energy balance h_air(T3) + f·ηb·LHV = (1+f)·h_prod(T4), LHV referenced to 298 K. */
const burnerF = (fu, T3, T4, eta_b) => { const a = gas(); return (h(a, T4) - h(a, T3)) / (eta_b * fu.LHV - (h(fu.gD, T4) - h(fu.gD, TREF))); };
/** Adiabatic flame temperature without dissociation for lean-to-stoichiometric f. */
function flameT(fu, T3, f, eta_b = 1) {
  const a = gas(), rhs = h(a, T3) + f * (eta_b * fu.LHV + h(fu.gD, TREF));
  return N.brent((T) => h(a, T) + f * h(fu.gD, T) - rhs, T3, 6000, 1e-6);
}
/** Convergent or fully-expanding nozzle. Returns effective (gross-thrust) velocity per unit mass flow. */
function nozzle(g, Tt, pt, pa, full, cv = 1) {
  if (!(pt > pa * 1.0001)) return { V: 0, Fg: 0, Ts: Tt, pe: pa, Me: 0, choked: false };
  const ph = phi(g, Tt), ht = h(g, Tt);
  let Ts = Tofphi(g, ph - g.R * Math.log(pt / pa), Tt * (pa / pt) ** 0.26), pe = pa, V = Math.sqrt(Math.max(0, 2 * (ht - h(g, Ts)))), choked = false;
  if (!full && V > Math.sqrt(gam(g, Ts) * g.R * Ts)) { // sonic throat: V = a(T*)
    Ts = N.brent((T) => 2 * (ht - h(g, T)) - gam(g, T) * g.R * T, Ts, Tt, 1e-7); V = Math.sqrt(2 * (ht - h(g, Ts)));
    pe = pt * Math.exp(-(ph - phi(g, Ts)) / g.R); choked = true;
  }
  return { V: cv * V, Fg: cv * V + ((pe - pa) * g.R * Ts) / (pe * V), Ts, pe, Me: V / Math.sqrt(gam(g, Ts) * g.R * Ts), choked };
}
const isShaft = (a) => a === 'turboprop' || a === 'turboshaft';
function inlet(i, alt, M, dISA, air = gas()) {
  const a0 = isa(alt, dISA), V0 = M * a0.a, Tt2 = V0 > 0 ? Tofh(air, h(air, a0.T) + 0.5 * V0 * V0, a0.T * (1 + 0.2 * M * M)) : a0.T;
  const ram = M > 1 ? 1 - 0.075 * (M - 1) ** 1.35 : 1; // MIL-E-5008B shock recovery, 1 < M < 5
  return { a0, V0, Tt2, pt0: a0.p * Math.exp((phi(air, Tt2) - phi(air, a0.T)) / air.R), pt2: a0.p * Math.exp((phi(air, Tt2) - phi(air, a0.T)) / air.R) * i.pi_d * ram, ram };
}

/**
 * Parametric (design-point) cycle per kg/s of core air. o: arch, alt, M, dISA, opr, fpr, bpr, T4, ep_c, ep_f, ep_t, eta_b,
 * pi_b, pi_d, eta_m, eta_gb, bleed, npr, cv, fuel, eta_prop, full (fully-expanding nozzle), ideal (perfect gas, no fuel mass).
 */
function cycle(o) {
  const fu = fuelOf(o.fuel), id = !!o.ideal, air = gas(0, null, id), R = air.R, I = inlet(o, o.alt, o.M, o.dISA, air), { a0, V0, Tt2, pt2 } = I;
  const fan = o.arch === 'turbofan', shaft = isShaft(o.arch), bpr = fan ? o.bpr : 0, fpr = fan ? Math.min(o.fpr, o.opr) : 1, st = [], warn = [];
  const sOf = (g, T, p) => phi(g, T) - phi(g, a0.T) - g.R * Math.log(p / a0.p);
  const add = (id_, name, g, T, p, m) => st.push({ id: id_, name, T, p, m, s: sOf(g, T, p) });
  const comp = (Tt, pt, PR, ep) => { const T2 = Tofphi(air, phi(air, Tt) + (R * Math.log(PR)) / ep, Tt * PR ** (0.2857 / ep)); return { Tt: T2, pt: pt * PR, w: h(air, T2) - h(air, Tt) }; };
  add('0', 'Ambient (static)', air, a0.T, a0.p, 1 + bpr); add('2', 'Compressor / fan face', air, Tt2, pt2, 1 + bpr);
  const s13 = fan ? comp(Tt2, pt2, fpr, o.ep_f) : { Tt: Tt2, pt: pt2, w: 0 };
  if (fan) add('13', 'Fan exit', air, s13.Tt, s13.pt, 1 + bpr);
  const s3 = comp(s13.Tt, s13.pt, o.opr / fpr, o.ep_c); add('3', 'Compressor exit', air, s3.Tt, s3.pt, 1);
  const eps = id ? 0 : o.bleed, T4 = Math.max(o.T4, s3.Tt + 20);
  if (o.T4 < s3.Tt + 20) warn.push('Turbine inlet temperature is at or below the compressor delivery temperature: no heat can be added at this pressure ratio.');
  const f = id ? (cp(air, 300) * (T4 - s3.Tt)) / (o.eta_b * fu.LHV) : burnerF(fu, s3.Tt, T4, o.eta_b);
  const g4 = gas(f, fu, id), pt4 = s3.pt * o.pi_b, m4 = id ? 1 : (1 - eps) * (1 + f); add('4', 'Turbine inlet', g4, T4, pt4, m4);
  // HP turbine drives the core compressor; cooling air (taken at station 3) rejoins after the rotor without doing work
  const wHPT = s3.w / o.eta_m, T45i = Tofh(g4, h(g4, T4) - wHPT / m4, 0.8 * T4), pt45 = pt4 * Math.exp(-(phi(g4, T4) - phi(g4, T45i)) / (o.ep_t * g4.R));
  const fm = id ? 0 : f * (1 - eps), gm = gas(fm, fu, id), m45 = 1 + fm, T45 = eps > 0 ? Tofh(gm, (m4 * h(g4, T45i) + eps * h(air, s3.Tt)) / m45, T45i) : T45i;
  add('45', 'HP turbine exit (mixed)', gm, T45, pt45, m45);
  let T5 = T45, pt5 = pt45, wLPT = 0, wPT = 0;
  if (fan) { wLPT = ((1 + bpr) * s13.w) / o.eta_m; T5 = Tofh(gm, h(gm, T45) - wLPT / m45, 0.8 * T45); pt5 = pt45 * Math.exp(-(phi(gm, T45) - phi(gm, T5)) / (o.ep_t * gm.R)); }
  else if (shaft) {
    pt5 = a0.p * o.npr;
    if (pt45 > pt5) { T5 = Tofphi(gm, phi(gm, T45) - o.ep_t * gm.R * Math.log(pt45 / pt5), 0.8 * T45); wPT = m45 * (h(gm, T45) - h(gm, T5)); } else { pt5 = pt45; warn.push('The gas generator leaves no pressure for the power turbine: shaft power is zero.'); }
  }
  if (fan || shaft) add('5', shaft ? 'Power turbine exit' : 'LP turbine exit', gm, T5, pt5, m45);
  if (pt5 < a0.p) warn.push('Turbine exit pressure is below ambient: the cycle cannot drive its compressors at this temperature and pressure ratio.');
  const n9 = nozzle(gm, T5, pt5, a0.p, !!o.full, o.cv), n19 = fan ? nozzle(air, s13.Tt, s13.pt, a0.p, !!o.full, o.cv) : { Fg: 0, V: 0, Ts: s13.Tt, pe: a0.p, choked: false };
  add('9', 'Core nozzle exit (static)', gm, n9.Ts, n9.pe, m45); if (fan) add('19', 'Bypass nozzle exit (static)', air, n19.Ts, n19.pe, bpr);
  const wshaft = wPT * o.eta_m * o.eta_gb, Fj = m45 * n9.Fg - V0, ff = id ? f : fm, Q = ff * fu.LHV;
  const Fs = shaft ? (o.arch === 'turboprop' ? Fj : 0) : m45 * n9.Fg + bpr * n19.Fg - (1 + bpr) * V0;
  const dKE = 0.5 * (m45 * n9.Fg ** 2 + bpr * n19.Fg ** 2 - (1 + bpr) * V0 * V0);
  const eta_th = shaft ? (wshaft + (o.arch === 'turboprop' ? Math.max(0, dKE) : 0)) / Q : dKE / Q;
  const eta_p = shaft ? o.eta_prop : dKE > 0 ? (Fs * V0) / dKE : 0, eta_o = shaft ? (wshaft * o.eta_prop + Math.max(0, Fs) * V0) / Q : (Fs * V0) / Q;
  const lossMech = (wHPT + wLPT) * (1 - o.eta_m) + wPT * (1 - o.eta_m * o.eta_gb);
  return { st, warn, fu, a0, V0, Tt2, pt2, T3: s3.Tt, p3: s3.pt, T4, f, ff, far4: f, T45, T5, pt5, pt45, pt4, n9, n19, Fs, Fj, wshaft, eta_th, eta_p, eta_o, lossMech,
    w: { fan: s13.w, hpc: s3.w, hpt: wHPT, lpt: wLPT, pt: wPT }, m45, bpr, fpr, opr: o.opr, valid: !warn.length && (shaft ? wshaft > 0 : Fs > 0) };
}
/** Design point = sea-level static ISA at rated turbine temperature; core mass flow sized to the rated thrust or shaft power. */
function design(i) {
  const T4 = i.tit_K / Math.max(1, i.theta_break), cy = cycle({ ...i, alt: 0, M: 0, dISA: 0, T4 });
  return { cy, T4, mdot: isShaft(i.arch) ? i.P0_W / Math.max(cy.wshaft, 1) : i.T0_N / Math.max(cy.Fs, 1) };
}
/**
 * Simplified matched off-design point. Assumptions (stated in every result): choked turbines with constant turbine temperature
 * ratio, so compressor and fan specific work scale with T4; constant polytropic efficiencies and bypass ratio; core flow from the
 * choked HP nozzle (m ∝ pt4/√T4); control holds T4/Tt2 up to the theta break and T4 = limit beyond it. thr scales T4.
 */
function offDesign(i, des, alt, M, dISA, thr = 1) {
  const air = gas(), I = inlet(i, alt, M, dISA), T4 = thr * i.tit_K * Math.min(1, I.Tt2 / T_SL / Math.max(1, i.theta_break)), r = T4 / des.T4, fan = i.arch === 'turbofan';
  const T13 = fan ? Tofh(air, h(air, I.Tt2) + des.cy.w.fan * r, I.Tt2 * 1.1) : I.Tt2, fpr = fan ? Math.exp((i.ep_f * (phi(air, T13) - phi(air, I.Tt2))) / air.R) : 1;
  const T3 = Tofh(air, h(air, T13) + des.cy.w.hpc * r, T13 * 2), opr = Math.max(1.05, fpr * Math.exp((i.ep_c * (phi(air, T3) - phi(air, T13))) / air.R));
  const mdot = des.mdot * ((I.pt2 * opr) / (des.cy.pt2 * des.cy.opr)) * Math.sqrt(des.T4 / T4), cy = cycle({ ...i, alt, M, dISA, opr, fpr, T4 });
  const P = Math.max(0, mdot * cy.wshaft), wf = mdot * cy.ff, sh = isShaft(i.arch);
  const Tprop = sh ? Math.min(staticThrust(i.T0_N, P, i.D_prop, I.a0.rho) * (i.T0_N > 0 ? I.a0.sigma ** 0.7 : 1), (i.eta_prop * P) / Math.max(I.V0, 1)) : 0;
  const T = sh ? Tprop + Math.max(0, mdot * cy.Fs) : Math.max(0, mdot * cy.Fs);
  const Wc = (mdot * Math.sqrt(I.Tt2 / T_SL)) / (I.pt2 / des.cy.pt2) / des.mdot, Nc = Math.sqrt((r * T_SL) / I.Tt2);
  return { cy, mdot, T4, opr, fpr, T, P, wf, tsfc: T > 0 ? wf / T : NaN, bsfc: P > 0 ? wf / P : NaN, Wc, Nc, I };
}
/** Static thrust per engine: the rating when given, otherwise actuator-disk momentum theory with a 0.75 figure of merit. */
const staticThrust = (T0, P, D, rho) => (T0 > 0 ? T0 : (0.75 * P * Math.sqrt(2 * rho * 0.25 * Math.PI * Math.max(D, 0.05) ** 2)) ** (2 / 3));
/** Simple T3–P3 NOx correlation EINOx = a·(P3/1 atm)^n·exp(T3/Ts) [g NO2 per kg fuel]; coefficients are calibratable inputs. */
const einox = (i, T3, p3) => i.nox_a * (p3 / P_SL) ** i.nox_n * Math.exp(T3 / i.nox_Ts);

// ---- shared inputs ----------------------------------------------------------------------------
const FLIGHT = [
  { key: 'alt_m', label: 'Altitude', unit: 'm', default: 10668, min: -500, max: 20000, group: 'Flight condition', help: 'Engine inlet boundary condition: ISA static pressure and temperature' },
  { key: 'mach', label: 'Flight Mach number', unit: '-', default: 0.78, min: 0, max: 3, group: 'Flight condition' },
  { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -60, max: 50, group: 'Flight condition' },
  { key: 'throttle', label: 'Throttle (T4 / rated T4)', unit: '-', default: 1, min: 0.5, max: 1.05, group: 'Flight condition', help: '1 = maximum rating at the flight condition' },
];
const ENGINE = [
  { key: 'arch', label: 'Architecture', type: 'select', options: ['turbofan', 'turbojet', 'turboprop', 'turboshaft'], default: 'turbofan', group: 'Engine', help: 'Two-spool separate-flow turbofan, single-spool turbojet, or gas generator with free power turbine' },
  { key: 'n_eng', label: 'Engines', unit: '', default: 2, min: 1, step: 1, discrete: true, group: 'Engine' },
  { key: 'T0_N', label: 'Rated static thrust per engine', unit: 'N', default: 120000, min: 0, group: 'Engine', help: 'Sizes the core mass flow for jets; propeller static thrust for turboprops (0 = estimate by momentum theory)' },
  { key: 'P0_W', label: 'Rated shaft power per engine', unit: 'W', default: 1.8e6, min: 0, group: 'Engine', help: 'Sizes the core mass flow for turboprop and turboshaft engines' },
  { key: 'opr', label: 'Overall pressure ratio', unit: '-', default: 32, min: 2, max: 70, group: 'Cycle' },
  { key: 'fpr', label: 'Fan pressure ratio', unit: '-', default: 1.6, min: 1.1, max: 4, group: 'Cycle', help: 'Turbofan only; about 1.4–1.7 for bypass ratios of 5–10' },
  { key: 'bpr', label: 'Bypass ratio', unit: '-', default: 5.7, min: 0, max: 20, group: 'Cycle', help: 'Turbofan only' },
  { key: 'tit_K', label: 'Rated turbine inlet temperature T4', unit: 'K', default: 1650, min: 900, max: 2100, group: 'Cycle' },
  { key: 'theta_break', label: 'Throttle ratio (theta break)', unit: '-', default: 1.0, min: 1, max: 1.2, group: 'Cycle', help: 'Inlet total temperature ratio Tt2/288.15 K at which the T4 limit is reached; 1 = flat-rated to ISA sea level' },
  { key: 'ep_c', label: 'Compressor polytropic efficiency', unit: '-', default: 0.9, min: 0.7, max: 1, group: 'Components' },
  { key: 'ep_f', label: 'Fan polytropic efficiency', unit: '-', default: 0.9, min: 0.7, max: 1, group: 'Components' },
  { key: 'ep_t', label: 'Turbine polytropic efficiency', unit: '-', default: 0.89, min: 0.7, max: 1, group: 'Components' },
  { key: 'eta_b', label: 'Combustion efficiency', unit: '-', default: 0.995, min: 0.8, max: 1, group: 'Components' },
  { key: 'pi_b', label: 'Combustor pressure ratio', unit: '-', default: 0.95, min: 0.8, max: 1, group: 'Components' },
  { key: 'pi_d', label: 'Inlet pressure recovery (subsonic)', unit: '-', default: 0.99, min: 0.8, max: 1, group: 'Components', help: 'Above Mach 1 the MIL-E-5008B shock recovery 1 − 0.075·(M − 1)^1.35 is applied in addition' },
  { key: 'eta_m', label: 'Spool mechanical efficiency', unit: '-', default: 0.99, min: 0.9, max: 1, group: 'Components' },
  { key: 'eta_gb', label: 'Reduction gearbox efficiency', unit: '-', default: 0.985, min: 0.9, max: 1, group: 'Components', help: 'Turboprop / turboshaft output gearbox' },
  { key: 'bleed', label: 'Turbine cooling air fraction', unit: '-', default: 0.1, min: 0, max: 0.3, group: 'Components', help: 'Share of core air bypassing the combustor; rejoins after the HP turbine' },
  { key: 'npr', label: 'Exhaust total pressure / ambient (shaft engines)', unit: '-', default: 1.08, min: 1.0, max: 2, group: 'Components', help: 'Pressure left after the power turbine: 1.02–1.05 turboshaft, 1.05–1.2 turboprop' },
  { key: 'cv', label: 'Nozzle velocity coefficient', unit: '-', default: 0.985, min: 0.9, max: 1, group: 'Components' },
  { key: 'fuel', label: 'Fuel', type: 'select', options: FUELS, default: 'Jet A-1', group: 'Fuel' },
  { key: 'eta_prop', label: 'Propeller / rotor propulsive efficiency', unit: '-', default: 0.82, min: 0.1, max: 0.95, group: 'Engine', help: 'Converts shaft power to thrust for shaft engines' },
  { key: 'D_prop', label: 'Propeller / rotor diameter per engine', unit: 'm', default: 3.9, min: 0.05, group: 'Engine', help: 'Only used to estimate static thrust when no rating is given' },
  { key: 'nox_a', label: 'NOx correlation coefficient a', unit: 'g/kg', default: 0.0986, min: 0, max: 1, group: 'Emissions', help: 'EINOx = a·(P3/1 atm)^n·exp(T3/Ts). Calibrate a against engine certification data.' },
  { key: 'nox_n', label: 'NOx pressure exponent n', unit: '-', default: 0.4, min: 0, max: 1, group: 'Emissions' },
  { key: 'nox_Ts', label: 'NOx temperature scale Ts', unit: 'K', default: 194.4, min: 100, max: 400, group: 'Emissions' },
];
const GT = ['turbofan', 'turbojet', 'turboprop', 'turboshaft'];
const isGT = (c) => (GT.includes(c.prop.type) ? true : `Gas-turbine analysis; this case has a ${c.prop.type} powerplant. Use the ${c.prop.type === 'piston' ? 'piston-engine' : 'electric / hybrid chain'} analysis.`);
const fuelName = (c) => (FUELS.includes(c.prop.fuel) ? c.prop.fuel : undefined);
const machOf = (c) => c.mission.cruise_mach || (c.mission.cruise_V_ms || c.flight.V_ms) / isa(c.mission.cruise_alt_m || c.atm.alt_m, c.atm.dISA_K).a;
const diaOf = (c) => c.prop.prop_dia_m || (c.rotor.R_m > 0 ? (2 * c.rotor.R_m) / Math.sqrt(c.meta.type === 'helicopter' ? Math.max(1, c.prop.n_eng) : 1) : undefined);
function gtDefaults(c) {
  const p = c.prop, arch = GT.includes(p.type) ? p.type : 'turbofan', sh = isShaft(arch), small = sh || p.T0_N < 30000;
  return { arch, alt_m: c.mission.cruise_alt_m || c.atm.alt_m, mach: machOf(c), dISA: c.atm.dISA_K, n_eng: p.n_eng, T0_N: sh ? p.T0_N : p.T0_N || undefined, P0_W: p.P0_W || undefined,
    opr: p.opr, bpr: arch === 'turbofan' ? Math.max(0.2, p.bpr) : p.bpr, fpr: N.clamp(1.3 + 1.8 / Math.max(p.bpr, 0.6), 1.25, Math.min(3.5, p.opr / 2)), tit_K: p.tit_K,
    ep_c: small ? 0.84 : 0.9, ep_t: small ? 0.85 : 0.89, pi_b: small ? 0.94 : 0.95, bleed: N.clamp((p.tit_K - 1350) / 3000, 0, 0.2), npr: arch === 'turboshaft' ? 1.04 : 1.1, fuel: fuelName(c), eta_prop: p.eta_prop, D_prop: diaOf(c) };
}
const gtAssume = ['Thermally perfect gas with variable cp(T); no dissociation', 'Polytropic component efficiencies; cooling air rejoins after the HP turbine without doing work', 'Convergent nozzles (choked or unchoked); no installation, nacelle or power off-take losses', 'Core mass flow sized so the sea-level static ISA rating equals the rated thrust or shaft power'];
const odAssume = 'Off-design: choked turbines with constant temperature ratio (compressor work ∝ T4), constant bypass ratio and component efficiencies, core flow ∝ pt4/√T4; control holds T4/Tt2 up to the theta break';
const kp = (key, label, value, unit, x) => ({ key, label, value, unit, ...x });

// ---- (a) gas-turbine cycle ----------------------------------------------------------------------
const cyc = {
  id: 'cycle', title: 'Gas-turbine thermodynamic cycle', fidelity: 'reduced-order',
  summary: 'Station-by-station cycle of a turbojet, two-spool turbofan, turboprop or turboshaft with temperature-dependent gas properties: thrust or shaft power, fuel consumption, efficiencies, emissions and the T–s diagram.',
  equations: ['Brayton cycle equations', 'Conservation of mass, momentum and energy', 'Compressor and turbine energy equations', 'Compressor and turbine efficiency relations', 'Nozzle flow equations', 'Thrust equation', 'Combustion energy equation', 'Shaft power balance'],
  applicable: isGT, inputs: [...FLIGHT, ...ENGINE], defaults: gtDefaults,
  run(i) {
    const des = design(i), od = offDesign(i, des, i.alt_m, i.mach, i.dISA, i.throttle), cy = od.cy, sh = isShaft(i.arch), n = i.n_eng, fu = cy.fu;
    const wf = n * od.wf, T = n * od.T, P = n * od.P, ei = einox(i, cy.T3, cy.p3), Qf = wf * fu.LHV, useful = sh ? P : n * od.mdot * Math.max(0, cy.eta_th) * cy.ff * fu.LHV;
    const tsS = sh ? NaN : des.cy.ff / des.cy.Fs, Tstat = sh ? staticThrust(i.T0_N, i.P0_W, i.D_prop, 1.225) : i.T0_N, warnings = [...new Set([...des.cy.warn, ...cy.warn])];
    const tsfc = sh ? (Number.isFinite(od.bsfc) ? (od.bsfc * Math.max(od.I.V0, 1)) / i.eta_prop : 0) : od.tsfc;
    if (cy.T4 > 1900) warnings.push('Turbine inlet temperature above 1900 K is beyond current cooled-turbine practice and the no-dissociation gas model.');
    if (fu.co2 === 0) warnings.push('Hydrogen fuel: the NOx correlation was derived for kerosene combustors and is indicative only.');
    if (i.mach > 1) warnings.push(`Supersonic inlet: MIL-E-5008B shock pressure recovery ${od.I.ram.toFixed(3)} applied; convergent nozzles under-expand and lose thrust.`);
    if (!(T > 0)) warnings.push('No positive net thrust at this flight condition and throttle.');
    const core = cy.st.filter((s) => s.id !== '19'), byp = cy.st.filter((s) => ['2', '13', '19'].includes(s.id));
    const kpis = [
      kp('thrust_N', 'Net thrust at the flight point (all engines)', T, 'N', { status: T > 0 ? 'ok' : 'bad' }),
      kp('thrust_static_N', 'Sea-level static thrust per engine', Tstat, 'N', { note: sh ? 'Propeller/rotor static thrust: rating or momentum-theory estimate' : 'Design point: equals the rating' }),
      kp('tsfc_kg_Ns', sh ? 'Equivalent thrust-specific fuel consumption' : 'TSFC at the flight point', tsfc, 'kg/N/s', { note: sh ? 'BSFC·V/η_prop' : '' }),
      kp('P_shaft_W', 'Shaft power at the flight point (all engines)', P, 'W'),
      kp('eta_thermal', 'Thermal efficiency', cy.eta_th, '-', { status: cy.eta_th > 0.3 ? 'ok' : 'warn' }), kp('eta_propulsive', 'Propulsive efficiency', cy.eta_p, '-'), kp('eta_overall', 'Overall efficiency', cy.eta_o, '-'),
      kp('fuel_flow_kgs', 'Fuel flow (all engines)', wf, 'kg/s'), kp('T4_K', 'Turbine inlet temperature', cy.T4, 'K', { status: cy.T4 <= 1900 ? 'ok' : 'warn' }),
      kp('T3_K', 'Compressor delivery temperature', cy.T3, 'K', { status: cy.T3 < 1000 ? 'ok' : 'warn', note: 'Disc material limit is about 950–1000 K' }), kp('P3_Pa', 'Compressor delivery pressure', cy.p3, 'Pa'),
      kp('opr_flight', 'Overall pressure ratio at the flight point', od.opr, '-'), kp('EGT_K', 'Exhaust gas total temperature', cy.T5, 'K'),
      kp('EINOx_g_kg', 'NOx emission index', ei, 'g/kg'), kp('co2_kg_s', 'CO₂ emission (all engines)', wf * fu.co2, 'kg/s'), kp('h2o_kg_s', 'Water vapour emission (all engines)', wf * fu.h2o, 'kg/s'),
      kp('heat_rejection_W', 'Heat to oil system (bearing and gearbox losses)', n * od.mdot * cy.lossMech, 'W'), kp('exhaust_heat_W', 'Exhaust and jet residual power', Math.max(0, Qf - useful), 'W'),
      kp('mdot_core_kgs', 'Design core mass flow per engine', des.mdot, 'kg/s'), kp('far', 'Combustor fuel–air ratio', cy.f, '-'),
    ];
    if (sh) kpis.push(kp('bsfc_kg_Ws', 'BSFC at the flight point', od.bsfc, 'kg/W/s'), kp('bsfc_static_kg_Ws', 'BSFC at the sea-level static rating', des.cy.ff / des.cy.wshaft, 'kg/W/s', { note: `${(des.cy.ff / des.cy.wshaft * 3.6e9).toFixed(0)} g/kWh` }), kp('jet_thrust_N', 'Residual jet thrust (all engines)', n * od.mdot * Math.max(0, cy.Fj), 'N'));
    else kpis.push(kp('tsfc_static_kg_Ns', 'TSFC at the sea-level static rating', tsS, 'kg/N/s'), kp('specific_thrust_ms', 'Specific thrust (per total airflow)', cy.Fs / (1 + cy.bpr), 'N·s/kg'), kp('V9_ms', 'Core jet velocity', cy.n9.V, 'm/s'), kp('mdot_total_kgs', 'Design total airflow per engine', des.mdot * (1 + cy.bpr), 'kg/s'));
    return {
      kpis,
      plots: [
        { type: 'line', title: 'Temperature–entropy diagram at the flight point', xlabel: 'Specific entropy above ambient [J/kg/K]', ylabel: 'Temperature [K]', series: [{ name: 'Core stream', x: core.map((s) => s.s), y: core.map((s) => s.T), style: 'line+points' }, ...(byp.length === 3 ? [{ name: 'Bypass stream', x: byp.map((s) => s.s), y: byp.map((s) => s.T), style: 'line+points' }] : [])] },
        { type: 'bar', title: 'Where the fuel energy goes', ylabel: 'Power [MW]', categories: ['All engines'], stacked: true, series: [{ name: sh ? 'Shaft power' : 'Thrust power', y: [(sh ? P : T * od.I.V0) / 1e6] }, { name: sh ? 'Residual jet power' : 'Jet kinetic energy left in the wake', y: [Math.max(0, useful - (sh ? P : T * od.I.V0)) / 1e6] }, { name: 'Exhaust heat and losses', y: [Math.max(0, Qf - useful) / 1e6] }] },
      ],
      tables: [{ title: 'Station table at the flight point (totals unless marked static)', columns: ['Station', 'Location', 'Temperature [K]', 'Pressure [kPa]', 'Flow / core flow [-]', 'Entropy [J/kg/K]'], rows: cy.st.map((s) => [s.id, s.name, +s.T.toFixed(1), +(s.p / 1e3).toFixed(2), +s.m.toFixed(4), +s.s.toFixed(1)]) }],
      outputs: { bsfc_kg_Ws: sh ? od.bsfc : 0, EICO2_g_kg: fu.co2 * 1e3, EIH2O_g_kg: fu.h2o * 1e3 },
      warnings, models: [`${i.arch} parametric cycle with variable cp(T) for air and ${fu.name} products`, 'Simplified matched-engine scaling from the sea-level static design point', 'T3–P3 NOx correlation (empirical, calibratable)', 'CO₂ and H₂O from fuel stoichiometry (exact)'],
      assumptions: [...gtAssume, odAssume, 'Published totals are for all engines at the flight point; static thrust is per engine'],
    };
  },
  verify() {
    const base = { arch: 'turbojet', alt: 0, M: 0, dISA: 0, opr: 12, fpr: 1, bpr: 0, T4: 1400, ep_c: 1, ep_f: 1, ep_t: 1, eta_b: 1, pi_b: 1, pi_d: 1, eta_m: 1, eta_gb: 1, bleed: 0, npr: 1, cv: 1, fuel: 'Jet A-1', eta_prop: 1, full: true, ideal: true };
    const c = cycle(base), g = 1004.7 / (1004.7 - R_AIR), e = (g - 1) / g, T3 = 288.15 * 12 ** e, T5 = 1400 - (T3 - 288.15), V9 = Math.sqrt(2 * 1004.7 * T5 * (1 - (1 / (12 * (T5 / 1400) ** (1 / e))) ** e));
    const nz = nozzle(gas(0, null, true), 800, 3e5, 1e5, false), air = gas();
    const fu = fuelOf('Jet A-1'), f = burnerF(fu, 800, 1600, 1), cpD = (T) => cp(fu.gD, T);
    return [
      N.check('Ideal Brayton thermal efficiency 1 − PR^(−(γ−1)/γ)', c.eta_th, 1 - 12 ** -e, 1e-8, 'Ideal air-standard Brayton cycle'),
      N.check('Ideal turbojet static specific thrust', c.Fs, V9, 1e-8, 'Closed-form ideal cycle analysis'),
      N.check('Choked nozzle exit pressure pt/((γ+1)/2)^(γ/(γ−1))', nz.pe, 3e5 / ((g + 1) / 2) ** (g / (g - 1)), 1e-6, 'Isentropic sonic relation'),
      N.check('Analytic enthalpy integral of cp(T), 288–1600 K', h(air, 1600) - h(air, 288.15), N.simpson((T) => cp(air, T), 288.15, 1600, 400), 1e-8, 'Simpson quadrature of the same polynomial'),
      N.check('Analytic entropy function ∫cp dT/T, 288–1600 K', phi(air, 1600) - phi(air, 288.15), N.simpson((T) => cp(air, T) / T, 288.15, 1600, 400), 1e-8, 'Simpson quadrature'),
      N.check('Burner energy balance closes', f * fu.LHV, N.simpson((T) => cp(air, T), 800, 1600, 400) + f * N.simpson(cpD, TREF, 1600, 400), 1e-8, 'First law, LHV at 298 K'),
    ];
  },
  calibration: { params: [{ key: 'ep_c', min: 0.8, max: 0.94 }, { key: 'ep_t', min: 0.8, max: 0.93 }, { key: 'bleed', min: 0, max: 0.25 }], sweep: 'throttle', target: 'tsfc_kg_Ns', note: 'Supply test-cell or flight TSFC (or BSFC) against throttle setting at a known altitude and Mach number.' },
  recommend(res, i) {
    const o = res.outputs, out = [], fu = fuelOf(i.fuel);
    if (fu.co2 > 0) out.push({ severity: 'advise', title: 'Fuel choice is the largest CO₂ lever', detail: `This point emits ${(o.co2_kg_s * 3600).toFixed(0)} kg CO₂ per hour. ${fu.lc < 1 ? `With ${fu.name} the life-cycle figure is about ${(100 * fu.lc).toFixed(0)}% of that.` : 'A HEFA-type sustainable aviation fuel burns the same in this cycle but cuts life-cycle CO₂ by roughly 75%; hydrogen removes CO₂ at the exhaust but emits about 2.6 times more water vapour per unit of energy.'}`, action: 'Re-run with "SAF (HEFA-SPK)" or "Liquid hydrogen" as the fuel and compare fuel flow, water vapour and NOx in the combustor analysis.', basis: 'Fuel stoichiometry and life-cycle factor in the fuel database' });
    if (o.EINOx_g_kg > 30) out.push({ severity: 'warn', title: 'High NOx emission index', detail: `EINOx ≈ ${o.EINOx_g_kg.toFixed(1)} g/kg at T3 = ${o.T3_K.toFixed(0)} K.`, action: 'A lower pressure ratio or a lean-burn / staged combustor reduces NOx; weigh this against the fuel-burn penalty of the lower pressure ratio.', basis: 'T3–P3 NOx correlation (uncalibrated default)' });
    if (o.T3_K > 950) out.push({ severity: 'warn', title: 'Compressor delivery temperature is at the disc material limit', detail: `T3 = ${o.T3_K.toFixed(0)} K.`, action: 'Reduce the overall pressure ratio or check the last-stage disc and cooling-air temperature in Suite 12.', basis: 'Typical nickel-alloy HP compressor disc limit of 950–1000 K' });
    if (!isShaft(i.arch) && o.eta_propulsive < 0.7 && i.mach > 0.3) out.push({ severity: 'advise', title: 'Propulsive efficiency is the weak link', detail: `η_propulsive = ${o.eta_propulsive.toFixed(2)} against η_thermal = ${o.eta_thermal.toFixed(2)}.`, action: 'A higher bypass ratio with a lower fan pressure ratio lowers jet velocity and fuel burn; check nacelle drag and weight in Suite 23.', basis: 'Froude efficiency 2/(1 + Vj/V0)' });
    return out;
  },
};

// ---- (b) off-design map ---------------------------------------------------------------------------
const map = {
  id: 'map', title: 'Altitude–Mach performance map and throttle characteristic', fidelity: 'reduced-order',
  summary: 'Thrust (or shaft power) and specific fuel consumption across the flight envelope, and the part-throttle fuel-consumption loop at the flight point, from a simplified matched-engine model.',
  equations: ['Brayton cycle equations', 'Thrust equation', 'Shaft power balance', 'Compressor–turbine work matching'],
  applicable: isGT,
  inputs: [...FLIGHT, ...ENGINE,
    { key: 'h_top', label: 'Top altitude of the map', unit: 'm', default: 12500, min: 1000, max: 20000, group: 'Numerics' },
    { key: 'M_top', label: 'Highest Mach number of the map', unit: '-', default: 0.9, min: 0.1, max: 2.5, group: 'Numerics' },
    { key: 'nGrid', label: 'Grid points per axis', unit: '', default: 10, min: 4, max: 40, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c) => { const g = gtDefaults(c); return { ...g, h_top: Math.max(3000, 1.15 * g.alt_m), M_top: Math.max(0.3, Math.min(2.5, 1.15 * g.mach)) }; },
  run(i) {
    const des = design(i), sh = isShaft(i.arch), n = Math.round(i.nGrid), hs = N.linspace(0, i.h_top, n), Ms = N.linspace(0, i.M_top, n), ref = sh ? i.P0_W : i.T0_N;
    const pts = hs.map((hh) => Ms.map((M) => offDesign(i, des, hh, M, i.dISA, 1))), val = (p) => (sh ? p.P : p.T), sfc = (p) => (sh ? p.bsfc * 3.6e9 : p.tsfc * 1e6);
    const here = offDesign(i, des, i.alt_m, i.mach, i.dISA, 1), thr = N.linspace(0.62, 1, 16).map((t) => offDesign(i, des, i.alt_m, i.mach, i.dISA, t)).filter((p) => val(p) > 0 && Number.isFinite(sfc(p)));
    const best = thr.length ? thr[N.argmin(thr.map(sfc))] : here, hot = offDesign(i, des, 0, 0, 15, 1), rows = [0, Math.floor((n - 1) / 2), n - 1];
    const warnings = ['The map comes from a simplified matched-engine model, not from measured component maps: use it for trends and replace it with engine-deck data for performance guarantees.'];
    if (!(val(here) > 0)) warnings.push('The engine produces no useful output at the flight point.');
    return {
      kpis: [
        kp(sh ? 'map_power_W' : 'map_thrust_N', sh ? 'Maximum shaft power at the flight point (per engine)' : 'Maximum net thrust at the flight point (per engine)', val(here), sh ? 'W' : 'N'),
        kp('lapse_ratio', sh ? 'Power lapse: flight point / sea-level static' : 'Thrust lapse: flight point / sea-level static', val(here) / ref, '-'),
        kp(sh ? 'map_bsfc_kg_Ws' : 'map_tsfc_kg_Ns', 'Specific fuel consumption at the flight point, full throttle', sh ? here.bsfc : here.tsfc, sh ? 'kg/W/s' : 'kg/N/s'),
        kp('sfc_best', 'Best part-throttle specific fuel consumption', sfc(best), sh ? 'g/kWh' : 'mg/N/s'), kp('throttle_best', 'Output at the best-consumption throttle / maximum', val(best) / Math.max(val(here), 1e-9), '-'),
        kp('hot_day_ratio', 'ISA+15 K sea-level static output / rating', val(hot) / ref, '-', { note: 'Loss on a hot day with the T4 limit held' }),
        kp('opr_cruise', 'Overall pressure ratio at the flight point', here.opr, '-'), kp('T4_cruise_K', 'Turbine inlet temperature at the flight point', here.T4, 'K'), kp('mdot_cruise_kgs', 'Core mass flow at the flight point', here.mdot, 'kg/s'),
      ],
      plots: [
        { type: 'heat', title: sh ? 'Maximum shaft power / rated power' : 'Maximum net thrust / rated static thrust', xlabel: 'Mach number [-]', ylabel: 'Altitude [m]', zlabel: sh ? 'P / P0 [-]' : 'T / T0 [-]', x: Ms, y: hs, z: pts.map((r) => r.map((p) => val(p) / ref)), contours: 10, overlay: [{ name: 'Flight point', x: [i.mach], y: [i.alt_m] }] },
        { type: 'line', title: 'Specific fuel consumption at full throttle', xlabel: 'Mach number [-]', ylabel: sh ? 'BSFC [g/kWh]' : 'TSFC [mg/N/s]', series: rows.map((r) => ({ name: `${Math.round(hs[r])} m`, x: Ms, y: pts[r].map((p) => (val(p) > 0 ? sfc(p) : NaN)) })) },
        { type: 'line', title: 'Throttle characteristic at the flight point', xlabel: sh ? 'Shaft power per engine [kW]' : 'Net thrust per engine [kN]', ylabel: sh ? 'BSFC [g/kWh]' : 'TSFC [mg/N/s]', series: [{ name: 'Part-throttle loop', x: thr.map((p) => val(p) / 1e3), y: thr.map(sfc), style: 'line+points' }] },
      ],
      tables: [{ title: 'Throttle sweep at the flight point', columns: ['T4 [K]', 'OPR [-]', sh ? 'Power [kW]' : 'Thrust [kN]', sh ? 'BSFC [g/kWh]' : 'TSFC [mg/N/s]', 'Fuel flow [kg/s]', 'Corrected speed [-]'], rows: thr.map((p) => [+p.T4.toFixed(0), +p.opr.toFixed(2), +(val(p) / 1e3).toFixed(2), +sfc(p).toFixed(2), +p.wf.toFixed(4), +p.Nc.toFixed(3)]) }],
      warnings, models: ['Simplified matched-engine off-design model', 'Variable-cp parametric cycle at each grid point'], assumptions: [odAssume, ...gtAssume.slice(0, 3)],
    };
  },
  verify() {
    const i = { ...Object.fromEntries(ENGINE.map((f) => [f.key, f.default])), theta_break: 1 }, des = design(i), d = offDesign(i, des, 0, 0, 0, 1), up = offDesign(i, des, 8000, 0, 0, 1), a = isa(8000);
    return [
      N.check('Off-design model reproduces the design thrust', d.T, i.T0_N, 1e-9, 'Consistency at the design point'),
      N.check('Corrected flow follows m·√θ/δ ∝ pressure ratio at constant T4/Tt2', ((up.mdot * Math.sqrt(a.theta)) / a.delta / des.mdot) * (i.opr / up.opr), 1, 1e-7, 'Choked-turbine flow function m·√T4/pt4 = constant'),
      N.check('Pressure ratio is invariant at constant T4/Tt2', up.opr, i.opr, 2e-2, 'Similarity; exact only for constant cp, so 2% is allowed for the cp change with inlet temperature'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.throttle_best < 0.97) out.push({ severity: 'info', title: 'Best fuel consumption is below full throttle', detail: `Minimum specific fuel consumption occurs at ${(100 * o.throttle_best).toFixed(0)}% of maximum output at this flight point.`, action: 'Size the engine so cruise falls near this throttle setting; an oversized engine cruises too far down the loop and burns more fuel.', basis: 'Part-throttle loop of the matched-engine model' });
    if (o.hot_day_ratio < 0.9) out.push({ severity: 'advise', title: 'Significant hot-day loss', detail: `Sea-level static output falls to ${(100 * o.hot_day_ratio).toFixed(0)}% of the rating at ISA+15 K.`, action: 'Use this ratio for hot-day take-off and hover checks in Suite 5, or flat-rate the engine by raising the theta break.', basis: 'T4-limited operation above the theta break' });
    return out;
  },
};

// ---- generic compressor map (labelled generic; shared by the turbomachinery and transient analyses) -------------------
/** Parametric β-line map normalised to the design point (W = 1, N = 1). Not a measured map. */
function cmap(PRd, SMd, ep, etaD) {
  const k = 3.5 * ep, d0 = 0.08, dc = 0.04, bd = 1 - Math.sqrt(1 - d0 / (d0 + dc)); let PK = PRd * (1 + SMd), p = 2;
  const pk = (n) => (1 + (PK ** (1 / k) - 1) * n * n) ** k; // surge-line pressure ratio: work ∝ N²
  for (let it = 0; it < 60; it++) { p = (0.75 * Math.log(pk(0.8) / PK)) / Math.log(0.8); const nS = (1 / (1 - d0)) ** (1 / p), t = (PRd * (1 + SMd)) ** (1 / k); PK = (1 + (t - 1) / (nS * nS)) ** k; }
  const cb = (1 - (PRd - 1) / (PK - 1)) / (bd * bd), ns = N.linspace(0.3, 1.25, 60), ws = ns.map((n) => (1 - d0) * n ** p), ps = ns.map(pk);
  const at = (n, b) => ({ W: n ** p * (1 - d0 + (d0 + dc) * (1 - (1 - b) ** 2)), PR: 1 + (pk(n) - 1) * (1 - cb * b * b), eta: etaD * Math.max(0.5, 1 - 0.5 * (n - 1) ** 2 - 0.25 * (b - bd) ** 2) });
  return { at, bd, surgePR: (W) => N.interp1(ws, ps, W), sm: (W, PR) => N.interp1(ws, ps, W) / PR - 1, surge: { W: ws, PR: ps }, p };
}
const etaIs = (PR, ep) => (PR ** 0.2857 - 1) / (PR ** (0.2857 / ep) - 1);

// ---- (c) turbomachinery -----------------------------------------------------------------------------
const turbo = {
  id: 'turbo', title: 'Compressor map, surge margin and turbine work matching', fidelity: 'reduced-order',
  summary: 'Stage-count estimates from stage loading, a generic compressor map with the engine operating line and surge margin, and the compressor–turbine power balance on each spool.',
  equations: ['Compressor energy equation', 'Turbine energy equation', 'Compressor and turbine efficiency relations', 'Shaft power balance', 'Euler turbomachinery work (stage loading)'],
  applicable: isGT,
  inputs: [...ENGINE.filter((f) => !f.key.startsWith('nox')),
    { key: 'sm_design', label: 'Design surge margin of the generic map', unit: '-', default: 0.2, min: 0.05, max: 0.4, group: 'Map', help: 'Surge pressure ratio / operating pressure ratio − 1 at constant corrected flow. 0.15–0.25 is typical.' },
    { key: 'psi_c', label: 'Compressor stage loading Δh/U²', unit: '-', default: 0.33, min: 0.15, max: 0.6, group: 'Stages', help: 'Mean-radius work coefficient; 0.3–0.4 for axial stages' },
    { key: 'U_c', label: 'Compressor mean blade speed', unit: 'm/s', default: 380, min: 150, max: 550, group: 'Stages' },
    { key: 'psi_t', label: 'Turbine stage loading Δh/U²', unit: '-', default: 1.6, min: 0.8, max: 3, group: 'Stages' },
    { key: 'U_t', label: 'HP turbine mean blade speed', unit: 'm/s', default: 420, min: 150, max: 600, group: 'Stages' },
    { key: 'U_lp', label: 'LP / power turbine mean blade speed', unit: 'm/s', default: 280, min: 100, max: 500, group: 'Stages' }],
  defaults: gtDefaults,
  run(i) {
    const des = design(i), c = des.cy, m = des.mdot, mp = cmap(i.opr, i.sm_design, i.ep_c, etaIs(i.opr, i.ep_c));
    const line = N.linspace(0.6, 1, 17).map((t) => offDesign(i, des, 0, 0, 0, t)).filter((p) => p.opr > 1.2), sm = line.map((p) => mp.sm(p.Wc, p.opr)), jm = N.argmin(sm);
    const sx = [], sy = []; for (const n of [0.6, 0.7, 0.8, 0.9, 1.0, 1.05]) { for (const b of N.linspace(0, 1, 14)) { const q = mp.at(n, b); sx.push(q.W); sy.push(q.PR); } sx.push(NaN); sy.push(NaN); }
    const st = (dh, psi, U) => Math.max(1, Math.ceil(dh / (psi * U * U) - 1e-9)), nC = st(c.w.hpc, i.psi_c, i.U_c), nH = st(c.w.hpt / c.m45, i.psi_t, i.U_t), wl = c.w.lpt || c.w.pt, nL = wl > 0 ? st(wl / c.m45, i.psi_t, i.U_lp) : 0;
    const Pc = m * c.w.hpc, Ph = m * c.w.hpt, Pf = m * (1 + c.bpr) * c.w.fan, Pl = m * wl, sI = line.filter((p, j) => sm[j] < Infinity).map((p) => p.Wc), warnings = ['The compressor map is a generic parametric map scaled to the design point; it is not the map of any real compressor. Surge margins are indicative.'];
    if (sm[jm] < 0.1) warnings.push('The steady operating line comes within 10% of the generic surge line at part power: handling bleed or variable stators would be needed.');
    return {
      kpis: [
        kp('sm_design_pct', 'Surge margin at the design point', 100 * mp.sm(1, i.opr), '%'), kp('sm_min_pct', 'Minimum steady-state surge margin on the operating line', 100 * sm[jm], '%', { status: sm[jm] > 0.15 ? 'ok' : sm[jm] > 0.08 ? 'warn' : 'bad', note: 'At least 15% is usually kept for transients, distortion and deterioration' }),
        kp('n_stages_comp', i.arch === 'turbofan' ? 'HP compression stages (booster + HPC)' : 'Compressor stages', nC, '-'), kp('n_stages_hpt', 'HP turbine stages', nH, '-'), kp('n_stages_lpt', isShaft(i.arch) ? 'Power turbine stages' : 'LP turbine stages', nL, '-'),
        kp('P_comp_W', 'Core compressor power', Pc, 'W'), kp('P_hpt_W', 'HP turbine power', Ph, 'W'), kp('work_balance_residual', 'HP spool power balance residual', (Ph * i.eta_m - Pc) / Pc, '-'),
        kp('hpt_expansion_ratio', 'HP turbine expansion ratio', c.pt4 / c.pt45, '-'), kp('lpt_expansion_ratio', 'LP / power turbine expansion ratio', c.pt45 / c.pt5, '-'),
        kp('eta_c_isentropic', 'Compressor isentropic efficiency', etaIs(i.opr, i.ep_c), '-', { note: 'Equivalent to the polytropic input at this pressure ratio' }), kp('T45_K', 'Inter-turbine temperature', c.T45, 'K'),
      ],
      plots: [
        { type: 'line', title: 'Generic compressor map with the steady operating line', xlabel: 'Corrected flow / design [-]', ylabel: 'Pressure ratio [-]', series: [{ name: 'Speed lines 60–105%', x: sx, y: sy }, { name: 'Surge line', x: mp.surge.W.filter((w) => w < 1.25), y: mp.surge.PR.filter((_, j) => mp.surge.W[j] < 1.25), style: 'dash' }, { name: 'Operating line', x: line.map((p) => p.Wc), y: line.map((p) => p.opr), style: 'line+points' }, { name: 'Design point', x: [1], y: [i.opr], style: 'points' }] },
        { type: 'line', title: 'Surge margin along the operating line', xlabel: 'Corrected flow / design [-]', ylabel: 'Surge margin [%]', series: [{ name: 'Steady state', x: sI, y: sm.map((v) => 100 * v) }], annotations: [{ y: 15, label: 'Typical minimum' }] },
        { type: 'bar', title: 'Spool power balance at the design point', ylabel: 'Power [MW]', categories: ['HP spool', i.arch === 'turbofan' ? 'LP spool' : 'Output'], series: [{ name: 'Absorbed (compressor, fan or shaft)', y: [Pc / 1e6, (i.arch === 'turbofan' ? Pf : Pl * i.eta_m * i.eta_gb) / 1e6] }, { name: 'Delivered by turbine', y: [Ph / 1e6, Pl / 1e6] }] },
      ],
      tables: [{ title: 'Turbomachinery work matching (design point)', columns: ['Component', 'Specific work [kJ/kg]', 'Power [kW]', 'Stages [-]'], rows: [['Fan (all flow)', +(c.w.fan / 1e3).toFixed(1), +(Pf / 1e3).toFixed(0), i.arch === 'turbofan' ? 1 : 0], ['Core compressor', +(c.w.hpc / 1e3).toFixed(1), +(Pc / 1e3).toFixed(0), nC], ['HP turbine', +(c.w.hpt / c.m45 / 1e3).toFixed(1), +(Ph / 1e3).toFixed(0), nH], ['LP / power turbine', +(wl / c.m45 / 1e3).toFixed(1), +(Pl / 1e3).toFixed(0), nL]] }],
      warnings, models: ['Generic parametric β-line compressor map (labelled generic)', 'Stage count from mean-radius stage loading', 'Matched operating line from the simplified off-design model'],
      assumptions: [odAssume, 'One map represents the whole core compression system', 'Stage loading and blade speed are the same for all stages of a component'],
    };
  },
  verify() {
    const mp = cmap(20, 0.2, 0.9, 0.85), q = mp.at(1, mp.bd), r = N.kv(turbo.run({ ...Object.fromEntries(turbo.inputs.map((f) => [f.key, f.default])) }));
    return [
      N.check('Generic map passes through the design flow', q.W, 1, 1e-12, 'Map construction'), N.check('Generic map passes through the design pressure ratio', q.PR, 20, 1e-10, 'Map construction'),
      N.check('Design surge margin equals the requested value', mp.sm(1, 20), 0.2, 2e-3, 'Map construction (surge line interpolated on 60 points)'),
      N.check('HP turbine power × mechanical efficiency = compressor power', r.work_balance_residual + 1, 1, 1e-9, 'Shaft power balance'),
    ];
  },
  recommend(res) {
    const o = res.outputs;
    return o.sm_min_pct < 15 ? [{ severity: o.sm_min_pct < 8 ? 'warn' : 'advise', title: 'Limited surge margin at part power', detail: `Minimum steady-state margin on the generic map is ${o.sm_min_pct.toFixed(1)}%.`, action: 'Check the acceleration transient analysis, and plan handling bleed or variable stator vanes; confirm with the real compressor map from rig tests.', basis: 'Surge margin at constant corrected flow, 15% guideline' }] : [];
  },
};

// ---- (d) combustor ----------------------------------------------------------------------------------
const combustor = {
  id: 'combustor', title: 'Combustor energy balance, emissions and contrail water', fidelity: 'analytical',
  summary: 'Fuel–air ratio and flame temperature from an energy balance for kerosene, SAF or hydrogen, exact CO₂ and water emission indices, an empirical NOx estimate and the Schmidt–Appleman contrail criterion.',
  equations: ['Combustion energy equation', 'Chemical species conservation (complete-combustion stoichiometry)', 'First law of thermodynamics', 'Ideal and real gas equations of state'],
  applicable: isGT,
  inputs: [
    { key: 'fuel', label: 'Fuel', type: 'select', options: FUELS, default: 'Jet A-1', group: 'Fuel' },
    { key: 'T3', label: 'Combustor inlet temperature T3', unit: 'K', default: 800, min: 300, max: 1100, group: 'Boundary conditions' },
    { key: 'P3', label: 'Combustor inlet pressure P3', unit: 'Pa', default: 3.0e6, min: 5e4, max: 7e6, group: 'Boundary conditions' },
    { key: 'T4', label: 'Combustor exit temperature T4', unit: 'K', default: 1600, min: 700, max: 2200, group: 'Boundary conditions' },
    { key: 'mdot_air', label: 'Combustor air mass flow', unit: 'kg/s', default: 50, min: 0.001, group: 'Boundary conditions', help: 'All engines together' },
    { key: 'eta_b', label: 'Combustion efficiency', unit: '-', default: 0.995, min: 0.8, max: 1, group: 'Combustor' },
    { key: 'phi_pz', label: 'Primary-zone equivalence ratio', unit: '-', default: 1.0, min: 0.3, max: 1, group: 'Combustor', help: 'Lean to stoichiometric only: 1.0 for a conventional rich-burn front end is the upper bound of this model, 0.5–0.7 for lean burn' },
    { key: 'nox_a', label: 'NOx correlation coefficient a', unit: 'g/kg', default: 0.0986, min: 0, max: 1, group: 'Emissions' },
    { key: 'nox_n', label: 'NOx pressure exponent n', unit: '-', default: 0.4, min: 0, max: 1, group: 'Emissions' },
    { key: 'nox_Ts', label: 'NOx temperature scale Ts', unit: 'K', default: 194.4, min: 100, max: 400, group: 'Emissions' },
    { key: 'alt_m', label: 'Altitude (for the contrail criterion)', unit: 'm', default: 10668, min: 0, max: 20000, group: 'Contrail' },
    { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -60, max: 50, group: 'Contrail' },
    { key: 'eta_o', label: 'Overall propulsion efficiency', unit: '-', default: 0.33, min: 0.05, max: 0.6, group: 'Contrail', help: 'Thrust power / fuel power at the flight point' },
  ],
  defaults(c) {
    const i = { ...Object.fromEntries([...FLIGHT, ...ENGINE].map((f) => [f.key, f.default])), ...Object.fromEntries(Object.entries(gtDefaults(c)).filter(([, v]) => v !== undefined)) };
    const od = offDesign(i, design(i), i.alt_m, i.mach, i.dISA, 1), ok = od.cy.eta_o > 0.05 && od.cy.eta_o < 0.6;
    return { fuel: i.fuel, T3: od.cy.T3, P3: od.cy.p3, T4: od.cy.T4, mdot_air: od.mdot * (1 - i.bleed) * i.n_eng, alt_m: i.alt_m, dISA: i.dISA, eta_o: ok ? od.cy.eta_o : undefined };
  },
  run(i) {
    const fu = fuelOf(i.fuel), T4 = Math.max(i.T4, i.T3 + 10), f = burnerF(fu, i.T3, T4, i.eta_b), eq = f / fu.fst, wf = f * i.mdot_air, a = isa(i.alt_m, i.dISA);
    const Tpz = flameT(fu, i.T3, i.phi_pz * fu.fst, i.eta_b), Tst = flameT(fu, i.T3, fu.fst, 1), ei = einox(i, i.T3, i.P3);
    // Schmidt–Appleman: slope of the exhaust mixing line G = EI_H2O·cp·p / (ε·Q·(1−η)); threshold temperature fit of Schumann (1996) at water saturation
    const G = (fu.h2o * 1004 * a.p) / (0.622 * fu.LHV * (1 - i.eta_o)), lg = Math.log(Math.max(G - 0.053, 1e-6)), Tlm = 273.15 - 46.46 + 9.43 * lg + 0.72 * lg * lg;
    const ps = N.linspace(0.2, 1, 25), names = ['Jet A-1', 'SAF (HEFA-SPK)', 'Liquid hydrogen'], fs = names.map(fuelOf), T3s = N.linspace(500, 950, 19), warnings = [];
    if (eq > 1) warnings.push('The requested exit temperature needs more than stoichiometric fuel: it cannot be reached at this inlet temperature.');
    if (Tpz > 2300) warnings.push('Primary-zone temperature above about 2300 K is over-predicted because dissociation is not modelled (typically by 100–250 K); thermal NOx forms rapidly above about 1900 K.');
    if (fu.co2 === 0) warnings.push('Hydrogen: the NOx correlation is a kerosene-combustor correlation; hydrogen burns leaner and hotter locally, so calibrate coefficient a before use.');
    if (i.T4 < i.T3 + 10) warnings.push('Exit temperature was raised to the inlet temperature plus 10 K.');
    return {
      kpis: [
        kp('far_comb', 'Fuel–air ratio', f, '-'), kp('phi_overall', 'Overall equivalence ratio', eq, '-', { status: eq < 1 ? 'ok' : 'bad' }), kp('far_stoich', 'Stoichiometric fuel–air ratio', fu.fst, '-'),
        kp('T_flame_pz_K', 'Primary-zone adiabatic flame temperature', Tpz, 'K', { status: Tpz < 2300 ? 'ok' : 'warn', note: 'No dissociation' }), kp('T_flame_stoich_K', 'Stoichiometric adiabatic flame temperature', Tst, 'K'),
        kp('fuel_flow_comb_kgs', 'Fuel flow', wf, 'kg/s'), kp('heat_release_W', 'Heat release rate', wf * i.eta_b * fu.LHV, 'W'),
        kp('EICO2_g_kg', 'CO₂ emission index', 1e3 * fu.co2, 'g/kg'), kp('EIH2O_g_kg', 'Water emission index', 1e3 * fu.h2o, 'g/kg'), kp('EINOx_comb_g_kg', 'NOx emission index (empirical)', ei, 'g/kg'),
        kp('co2_comb_kg_s', 'CO₂ emission rate', wf * fu.co2, 'kg/s'), kp('co2_lifecycle_kg_s', 'Life-cycle CO₂ rate', wf * fu.co2 * fu.lc, 'kg/s', { note: 'Combustion CO₂ × life-cycle factor of the fuel' }),
        kp('h2o_comb_kg_s', 'Water vapour emission rate', wf * fu.h2o, 'kg/s'), kp('nox_g_s', 'NOx emission rate', wf * ei, 'g/s'),
        kp('contrail_G_Pa_K', 'Contrail mixing-line slope G', G, 'Pa/K'), kp('contrail_T_threshold_K', 'Contrail threshold temperature (saturated air)', Tlm, 'K', { status: a.T < Tlm ? 'warn' : 'ok', note: `Ambient ${a.T.toFixed(1)} K; contrails can form when ambient is colder` }),
        kp('contrail_possible', 'Contrail formation possible', a.T < Tlm ? 1 : 0, '-'),
      ],
      plots: [
        { type: 'line', title: 'Adiabatic flame temperature at this inlet temperature', xlabel: 'Equivalence ratio [-]', ylabel: 'Flame temperature [K]', series: fs.map((F) => ({ name: F.name, x: ps, y: ps.map((p) => flameT(F, i.T3, p * F.fst, 1)) })), annotations: [{ y: 1900, label: 'Rapid thermal NOx' }] },
        { type: 'bar', title: 'Emissions per unit of fuel energy', ylabel: 'Emission [g/MJ]', categories: names, series: [{ name: 'CO₂ (combustion)', y: fs.map((F) => (1e9 * F.co2) / F.LHV) }, { name: 'CO₂ (life-cycle)', y: fs.map((F) => (1e9 * F.co2 * F.lc) / F.LHV) }, { name: 'H₂O', y: fs.map((F) => (1e9 * F.h2o) / F.LHV) }] },
        { type: 'line', title: 'NOx emission index versus combustor inlet temperature', xlabel: 'T3 [K]', ylabel: 'EINOx [g/kg]', series: [1, 0.5].map((s) => ({ name: `P3 = ${(s * i.P3 / 1e5).toFixed(1)} bar`, x: T3s, y: T3s.map((t) => einox(i, t, s * i.P3)) })) },
      ],
      tables: [{ title: 'Complete-combustion stoichiometry per kg of fuel', columns: ['Fuel', 'O₂ consumed [kg]', 'CO₂ [kg]', 'H₂O [kg]', 'Stoichiometric air/fuel [-]', 'LHV [MJ/kg]'], rows: fs.map((F) => [F.name, +F.o2.toFixed(3), +F.co2.toFixed(3), +F.h2o.toFixed(3), +(1 / F.fst).toFixed(2), +(F.LHV / 1e6).toFixed(2)]) }],
      warnings, models: ['First-law combustor balance with variable-cp products', 'Complete-combustion stoichiometry (exact CO₂ and H₂O)', 'T3–P3 NOx correlation (empirical, calibratable coefficient)', 'Schmidt–Appleman contrail criterion with Schumann threshold fit'],
      assumptions: ['No dissociation or finite-rate chemistry: flame temperatures are upper bounds above about 2000 K', 'Fuel enters at 298 K; LHV referenced to 298 K', 'Contrail threshold given for water-saturated ambient air; drier air needs colder temperatures'],
    };
  },
  verify() {
    const k = fuelOf('Jet A-1'), hy = fuelOf('Liquid hydrogen'), r = N.kv(combustor.run(Object.fromEntries(combustor.inputs.map((f) => [f.key, f.default])))), a = gas();
    const Tf = flameT(k, 700, 0.03, 1), g = gas(0.03, k);
    return [
      N.check('Kerosene (C12H23) CO₂ index 44.009·12/(12·12.011 + 23·1.008)', k.co2, (44.009 * 12) / (12 * 12.011 + 23 * 1.008), 3e-3, 'Stoichiometry; fuel database value 3.16 kg/kg'),
      N.check('Mass balance: fuel + O₂ = CO₂ + H₂O', k.co2 + k.h2o - k.o2, 1, 1e-4, 'Conservation of mass'),
      N.check('Hydrogen stoichiometric air/fuel ratio', 1 / hy.fst, (0.5 * 31.998) / 2.016 / Y_O2, 1e-3, '2H₂ + O₂ → 2H₂O with 23.14% O₂ by mass in air'),
      N.check('Hydrogen water index = 18.015/2.016', hy.h2o, 18.015 / 2.016, 1e-3, 'Stoichiometry'),
      N.check('Flame-temperature energy balance closes', 1.03 * (h(g, Tf) - h(g, TREF)), h(a, 700) - h(a, TREF) + 0.03 * k.LHV, 1e-8, 'First law'),
      N.check('Fuel–air ratio reproduces the exit temperature', flameT(k, 800, r.far_comb, 0.995), 1600, 1e-7, 'Inverse of the same balance'),
    ];
  },
  calibration: { params: [{ key: 'nox_a', min: 0.01, max: 0.5 }, { key: 'nox_n', min: 0.2, max: 0.7 }], sweep: 'T3', target: 'EINOx_comb_g_kg', note: 'Supply measured EINOx against combustor inlet temperature (e.g. the four certification thrust settings) for the engine family.' },
  recommend(res, i) {
    const o = res.outputs, fu = fuelOf(i.fuel), out = [];
    if (fu.co2 > 0 && fu.lc === 1) out.push({ severity: 'advise', title: 'Drop-in SAF cuts life-cycle CO₂ without changing the cycle', detail: `Combustion CO₂ is ${(o.co2_comb_kg_s * 3600).toFixed(0)} kg/h. HEFA-SPK has the same exhaust CO₂ but about a quarter of the life-cycle CO₂, with about 2% lower fuel mass flow from its higher heating value.`, action: 'Evaluate SAF blend ratios against fuel cost in Suite 26.', basis: 'Life-cycle factor in the fuel database' });
    if (o.contrail_possible) out.push({ severity: 'info', title: 'Persistent-contrail conditions are possible at this altitude', detail: `Ambient temperature is below the ${o.contrail_T_threshold_K.toFixed(0)} K threshold.${fu.co2 === 0 ? ' Hydrogen raises the threshold because it emits about 2.6 times more water per unit of energy.' : ''}`, action: 'Where humidity forecasts show ice-supersaturated air, a cruise-level change of 600 m is often enough to avoid persistent contrails; weigh against the fuel penalty in Suite 24.', basis: 'Schmidt–Appleman criterion' });
    if (o.T_flame_pz_K > 2200) out.push({ severity: 'advise', title: 'Hot primary zone drives thermal NOx', detail: `Primary-zone flame temperature ≈ ${o.T_flame_pz_K.toFixed(0)} K.`, action: 'A leaner primary zone (equivalence ratio 0.5–0.7) or staged combustion lowers flame temperature and NOx; check lean blow-out margin with a kinetics tool.', basis: 'Zeldovich thermal NOx grows rapidly above about 1900 K' });
    return out;
  },
};

// ---- (e) piston engine --------------------------------------------------------------------------------
/** Otto or Diesel cycle per kg of trapped air. ideal: perfect gas with heat q [J/kg]. Otherwise variable-cp fuel-air cycle. */
function pistonCycle(o) {
  const id = !!o.ideal, fu = o.fu, f = id ? 0 : o.phi * fu.fst, fb = Math.min(f, fu.fst), ga = gas(0, null, id), gp = id ? ga : gas(fb, fu), mm = 1 + f;
  const q = id ? o.q : fb * o.eta_c * fu.LHV * (1 - o.hl), u = (g, T) => h(g, T) - g.R * T, r = o.r;
  // isentrope at varying volume: φ(T) − R ln T + R ln v = const
  const isoV = (g, T1, vr, T) => { const tg = phi(g, T1) - g.R * Math.log(T1) - g.R * Math.log(vr); for (let k = 0; k < 40; k++) { const d = (tg - phi(g, T) + g.R * Math.log(T)) / (cp(g, T) - g.R); T *= Math.exp(N.clamp(d, -1, 1)); if (Math.abs(d) < 1e-12) break; } return T; };
  const T1 = o.T1, p1 = o.p1, v1 = (ga.R * T1) / p1, T2 = isoV(ga, T1, 1 / r, T1 * r ** 0.35), p2 = (ga.R * T2) / (v1 / r);
  let T3 = T2 + q / (mm * 900), rc = 1;
  if (o.diesel) { const tg = (h(ga, T2) - h(ga, TREF) + q) / mm + h(gp, TREF); T3 = Tofh(gp, tg, T3); rc = (mm * gp.R * T3) / (ga.R * T2); }
  else { const tg = (u(ga, T2) - u(ga, TREF) + q) / mm + u(gp, TREF); for (let k = 0; k < 60; k++) { const d = (u(gp, T3) - tg) / (cp(gp, T3) - gp.R); T3 -= d; if (Math.abs(d) < 1e-8) break; } }
  const v3 = (v1 / r) * rc, p3 = (mm * gp.R * T3) / v3, T4 = isoV(gp, T3, v1 / v3, T3 * (v3 / v1) ** 0.28), p4 = (mm * gp.R * T4) / v1;
  const w = (o.diesel ? p3 * (v3 - v1 / r) : 0) + mm * (u(gp, T3) - u(gp, T4)) - (u(ga, T2) - u(ga, T1));
  // p–V loop by quadrature (independent of the state-energy work above)
  const n = Math.max(4, Math.round(o.nPts || 0)), vc = N.linspace(v1, v1 / r, n + 1), pc = vc.map((v, k) => (k === 0 ? p1 : (ga.R * isoV(ga, T1, v / v1, T2)) / v));
  const ve = N.linspace(v3, v1, n + 1), pe = ve.map((v, k) => (k === 0 ? p3 : (mm * gp.R * isoV(gp, T3, v / v3, T4)) / v));
  const wpv = o.nPts ? N.trapz(vc, pc) + (o.diesel ? p3 * (v3 - v1 / r) : 0) + N.trapz(ve, pe) : NaN;
  return { T1, T2, T3, T4, p1, p2, p3, p4, v1, rc, w, wpv, q, f, eta: w / (id ? q : Math.max(f, 1e-12) * fu.LHV), imep: w / (v1 * (1 - 1 / r)), V: [...vc, ...ve, v1].map((v) => v / v1), P: [...pc, ...pe, p1] };
}
const piston = {
  id: 'piston', title: 'Piston engine cycle, altitude lapse and fuel consumption', fidelity: 'reduced-order',
  summary: 'Air-standard and variable-property fuel–air Otto or Diesel cycle with volumetric and mechanical efficiency: power, BSFC, the p–V diagram and power lapse with altitude.',
  equations: ['Otto cycle equations', 'Diesel cycle equations', 'First and second laws of thermodynamics', 'Ideal gas equation of state', 'Combustion energy equation'],
  applicable: (c) => (c.prop.type === 'piston' ? true : `Piston-engine analysis; this case has a ${c.prop.type} powerplant.`),
  inputs: [
    { key: 'cyc', label: 'Cycle', type: 'select', options: ['Otto (spark ignition)', 'Diesel (compression ignition)'], default: 'Otto (spark ignition)', group: 'Engine' },
    { key: 'n_eng', label: 'Engines', unit: '', default: 1, min: 1, step: 1, discrete: true, group: 'Engine' },
    { key: 'P0_W', label: 'Rated sea-level power per engine', unit: 'W', default: 119000, min: 10, group: 'Engine' },
    { key: 'disp_L', label: 'Displacement per engine (0 = size to rated power)', unit: 'L', default: 0, min: 0, max: 100, group: 'Engine' },
    { key: 'rpm', label: 'Crankshaft speed', unit: 'rpm', default: 2700, min: 300, max: 12000, group: 'Engine', help: 'Initial condition for the cycle; four-stroke operation assumed' },
    { key: 'r', label: 'Compression ratio', unit: '-', default: 8.5, min: 5, max: 24, group: 'Engine', help: '7–10 spark ignition, 15–20 Diesel' },
    { key: 'map_ratio', label: 'Manifold pressure / ambient', unit: '-', default: 0.96, min: 0.3, max: 3, group: 'Engine', help: 'Below 1 naturally aspirated (throttle and intake loss); above 1 turbo- or supercharged' },
    { key: 'phi', label: 'Equivalence ratio', unit: '-', default: 1.05, min: 0.4, max: 1.4, group: 'Engine', help: '1.1–1.25 full-rich take-off, 0.9–1.0 lean cruise; Diesel 0.4–0.7' },
    { key: 'eta_v', label: 'Volumetric efficiency', unit: '-', default: 0.85, min: 0.5, max: 1.1, group: 'Efficiencies' },
    { key: 'eta_mech', label: 'Mechanical efficiency at sea level', unit: '-', default: 0.88, min: 0.6, max: 0.98, group: 'Efficiencies', help: 'Brake / indicated power; sets a constant friction mean effective pressure' },
    { key: 'eta_c', label: 'Combustion efficiency', unit: '-', default: 0.97, min: 0.7, max: 1, group: 'Efficiencies' },
    { key: 'hl', label: 'Heat-loss fraction to walls', unit: '-', default: 0.15, min: 0, max: 0.4, group: 'Efficiencies', help: 'Share of the heat release lost to coolant and oil during combustion and expansion' },
    { key: 'diagram', label: 'Indicator-diagram factor', unit: '-', default: 0.9, min: 0.6, max: 1, group: 'Efficiencies', help: 'Real / ideal loop work: finite burn time, valve timing and pumping' },
    { key: 'dT_charge', label: 'Charge heating in the intake', unit: 'K', default: 25, min: 0, max: 120, group: 'Efficiencies' },
    { key: 'fuel', label: 'Fuel', type: 'select', options: ['Avgas 100LL', 'Jet A-1', 'SAF (HEFA-SPK)', 'Liquid hydrogen'], default: 'Avgas 100LL', group: 'Fuel' },
    { key: 'EINOx', label: 'NOx emission index (user-supplied)', unit: 'g/kg', default: 4, min: 0, max: 60, group: 'Fuel', help: 'Not predicted: piston NOx depends strongly on mixture strength. Enter test data.' },
    { key: 'alt_m', label: 'Altitude', unit: 'm', default: 2400, min: -500, max: 12000, group: 'Flight condition' },
    { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -60, max: 50, group: 'Flight condition' },
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 62, min: 0, max: 250, group: 'Flight condition' },
    { key: 'eta_prop', label: 'Propeller efficiency', unit: '-', default: 0.8, min: 0.1, max: 0.95, group: 'Propeller' },
    { key: 'T0_N', label: 'Static thrust per engine (0 = estimate)', unit: 'N', default: 0, min: 0, group: 'Propeller' },
    { key: 'D_prop', label: 'Propeller diameter', unit: 'm', default: 1.9, min: 0.05, group: 'Propeller' },
    { key: 'nPts', label: 'Points per stroke in the p–V loop', unit: '', default: 60, min: 4, max: 2000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c) => ({ n_eng: c.prop.n_eng, P0_W: c.prop.P0_W || undefined, rpm: c.prop.rpm || undefined, fuel: ['Avgas 100LL', 'Jet A-1', 'SAF (HEFA-SPK)', 'Liquid hydrogen'].includes(c.prop.fuel) ? c.prop.fuel : undefined, alt_m: c.mission.cruise_alt_m || c.atm.alt_m, dISA: c.atm.dISA_K, V: c.mission.cruise_V_ms || c.flight.V_ms, eta_prop: c.prop.eta_prop, T0_N: c.prop.T0_N, D_prop: diaOf(c) }),
  run(i) {
    const fu = fuelOf(i.fuel), diesel = i.cyc.startsWith('Diesel'), g = 1.4;
    const at = (alt, nPts = 0) => { const a = isa(alt, i.dISA); return { a, c: pistonCycle({ fu, r: i.r, T1: a.T + i.dT_charge, p1: a.p * i.map_ratio, phi: i.phi, eta_c: i.eta_c, hl: i.hl, diesel, nPts }), rhoM: (a.p * i.map_ratio) / (R_AIR * a.T) }; };
    const sl = at(0), cps = (i.rpm / 120), ind0 = i.eta_v * sl.rhoM * cps * sl.c.w * i.diagram; // indicated power per m³ of displacement at sea level
    const Vd = i.disp_L > 0 ? i.disp_L / 1e3 : i.P0_W / (ind0 * i.eta_mech), Pf = ind0 * Vd * (1 - i.eta_mech); // constant friction power
    const pb = (s) => { const ma = i.eta_v * s.rhoM * cps * Vd; return { ma, Pi: ma * s.c.w * i.diagram, P: Math.max(0, ma * s.c.w * i.diagram - Pf), wf: ma * s.c.f }; };
    const fl = at(i.alt_m, i.nPts), b = pb(fl), b0 = pb(sl), c = fl.c, n = i.n_eng, P = n * b.P, wf = n * b.wf, eta_b = b.P / (b.wf * fu.LHV), bsfc = b.P > 0 ? b.wf / b.P : NaN;
    const etaAir = diesel ? 1 - (i.r ** (1 - g) * (c.rc ** g - 1)) / (g * (c.rc - 1)) : 1 - i.r ** (1 - g);
    const hs = N.linspace(0, 9000, 31), lapse = hs.map((hh) => pb(at(hh)).P / b0.P), gf = hs.map((hh) => Math.max(0, (isa(hh, i.dISA).sigma - 0.1325) / 0.8675) / ((isa(0, i.dISA).sigma - 0.1325) / 0.8675));
    const Ts = staticThrust(i.T0_N, b0.P, i.D_prop, 1.225), T = n * Math.min(Ts * fl.a.sigma ** 0.7, (i.eta_prop * b.P) / Math.max(i.V, 1)), Qf = wf * fu.LHV;
    const Qcool = n * (b.ma * c.f * Math.min(1, 1 / i.phi) * i.eta_c * fu.LHV * i.hl + (b.Pi - b.P)), warnings = ['NOx is a user-supplied emission index: this suite has no predictive NOx model for piston engines.'];
    if (c.T3 > 2400) warnings.push('Peak cycle temperature is over-predicted (no dissociation, instantaneous combustion); treat it as an upper bound.');
    if (!diesel && i.r > 10.5 && i.map_ratio > 0.9) warnings.push('Compression ratio above about 10.5 at high manifold pressure risks detonation with aviation gasoline.');
    if (i.phi > 1) warnings.push('Rich mixture: the excess fuel is carried as inert mass and leaves unburned (CO and H₂ formation is not modelled).');
    if (!(b.P > 0)) warnings.push('Friction exceeds indicated power at this altitude: the engine delivers no brake power.');
    return {
      kpis: [
        kp('P_shaft_W', 'Brake power at altitude (all engines)', P, 'W', { status: P > 0 ? 'ok' : 'bad' }), kp('power_sl_W', 'Sea-level brake power per engine', b0.P, 'W'),
        kp('bsfc_kg_Ws', 'Brake specific fuel consumption', bsfc, 'kg/W/s', { note: `${(bsfc * 3.6e9).toFixed(0)} g/kWh` }), kp('eta_thermal', 'Brake thermal efficiency', eta_b, '-', { status: eta_b > 0.25 ? 'ok' : 'warn' }),
        kp('eta_air_standard', 'Air-standard cycle efficiency', etaAir, '-'), kp('eta_fuel_air', 'Fuel–air cycle efficiency', c.eta, '-'), kp('eta_overall', 'Overall efficiency (× propeller)', eta_b * i.eta_prop, '-'),
        kp('fuel_flow_kgs', 'Fuel flow (all engines)', wf, 'kg/s'), kp('thrust_N', 'Thrust at the flight point (all engines)', T, 'N'), kp('thrust_static_N', 'Sea-level static thrust per engine', Ts, 'N', { note: i.T0_N > 0 ? 'Rating' : 'Momentum-theory estimate' }),
        kp('tsfc_kg_Ns', 'Equivalent thrust-specific fuel consumption', T > 0 ? wf / T : 0, 'kg/N/s'), kp('T4_K', 'Peak cylinder gas temperature', c.T3, 'K', { note: 'Hot-section temperature for this engine type; upper bound' }),
        kp('p_peak_Pa', 'Peak cylinder pressure', c.p3, 'Pa'), kp('egt_K', 'Gas temperature at exhaust-valve opening', c.T4, 'K'), kp('map_Pa', 'Manifold pressure', c.p1, 'Pa'),
        kp('imep_Pa', 'Indicated mean effective pressure', b.Pi / (Vd * cps), 'Pa'), kp('bmep_Pa', 'Brake mean effective pressure', b.P / (Vd * cps), 'Pa'),
        kp('displacement_L', 'Displacement per engine', Vd * 1e3, 'L'), kp('cutoff_ratio', 'Cut-off ratio', c.rc, '-'), kp('work_pv_J_kg', 'Loop work by ∮p dV', c.wpv, 'J/kg', { note: 'Per kg of trapped air; compare with the state-energy work' }), kp('work_cycle_J_kg', 'Loop work from state energies', c.w, 'J/kg'),
        kp('EINOx_g_kg', 'NOx emission index (user-supplied)', i.EINOx, 'g/kg'), kp('co2_kg_s', 'CO₂ emission (all engines)', wf * fu.co2, 'kg/s'), kp('heat_rejection_W', 'Heat to coolant and oil (all engines)', Qcool, 'W'), kp('exhaust_heat_W', 'Exhaust heat and unburned fuel', Math.max(0, Qf - P - Qcool), 'W'),
        kp('lapse_ratio', 'Power at altitude / sea level', b.P / b0.P, '-'),
      ],
      plots: [
        { type: 'line', title: 'Pressure–volume diagram (ideal fuel–air loop)', xlabel: 'Volume / maximum volume [-]', ylabel: 'Pressure [bar]', series: [{ name: `${diesel ? 'Diesel' : 'Otto'} cycle at ${Math.round(i.alt_m)} m`, x: c.V, y: c.P.map((p) => p / 1e5) }] },
        { type: 'line', title: 'Power lapse with altitude', xlabel: 'Altitude [m]', ylabel: 'Power / sea-level power [-]', series: [{ name: 'Cycle with constant friction', x: hs, y: lapse }, { name: 'Gagg–Ferrar (σ − 0.1325)/0.8675', x: hs, y: gf, style: 'dash' }], annotations: [{ x: i.alt_m, label: 'Flight point' }] },
        { type: 'bar', title: 'Fuel energy balance at the flight point', ylabel: 'Power [kW]', categories: ['All engines'], stacked: true, series: [{ name: 'Brake power', y: [P / 1e3] }, { name: 'Coolant, oil and friction', y: [Qcool / 1e3] }, { name: 'Exhaust and unburned fuel', y: [Math.max(0, Qf - P - Qcool) / 1e3] }] },
      ],
      tables: [{ title: 'Cycle states (per kg of trapped air)', columns: ['State', 'Temperature [K]', 'Pressure [bar]'], rows: [['1 Start of compression', c.T1, c.p1], ['2 End of compression', c.T2, c.p2], ['3 End of combustion', c.T3, c.p3], ['4 End of expansion', c.T4, c.p4]].map((r) => [r[0], +r[1].toFixed(1), +(r[2] / 1e5).toFixed(3)]) }],
      outputs: { EICO2_g_kg: fu.co2 * 1e3 },
      warnings, models: [`${diesel ? 'Diesel' : 'Otto'} fuel–air cycle with variable cp(T) and ${fu.name} products`, 'Volumetric efficiency and constant friction mean effective pressure', 'Gagg–Ferrar altitude lapse (shown for comparison)', 'CO₂ from fuel stoichiometry (exact)'],
      assumptions: ['Four-stroke engine; instantaneous combustion; no dissociation; residual gas neglected', 'Heat loss removed from the heat release; the indicator-diagram factor covers finite burn time and pumping', 'Friction power is constant with altitude at fixed speed', 'Manifold-to-ambient pressure ratio is constant with altitude (no critical altitude modelled for turbocharged engines)'],
    };
  },
  convergence: { param: 'nPts', label: 'Points per stroke in the p–V quadrature', levels: [10, 20, 40, 80, 160], metric: 'work_pv_J_kg' },
  verify() {
    const fu = fuelOf('Avgas 100LL'), g = 1004.7 / (1004.7 - R_AIR), o = pistonCycle({ ideal: true, fu, r: 9, T1: 300, p1: 1e5, q: 1.5e6, nPts: 4000 });
    const d = pistonCycle({ ideal: true, diesel: true, fu, r: 18, T1: 300, p1: 1e5, q: 1.2e6 }), fa = pistonCycle({ fu, r: 8.5, T1: 320, p1: 9e4, phi: 0.9, eta_c: 1, hl: 0, nPts: 4000 });
    return [
      N.check('Ideal Otto efficiency 1 − r^(1−γ)', o.eta, 1 - 9 ** (1 - g), 1e-8, 'Air-standard Otto cycle'),
      N.check('Ideal Diesel efficiency with cut-off ratio', d.eta, 1 - (18 ** (1 - g) * (d.rc ** g - 1)) / (g * (d.rc - 1)), 1e-8, 'Air-standard Diesel cycle'),
      N.check('Ideal Diesel cut-off ratio 1 + q/(cp·T2)', d.rc, 1 + 1.2e6 / (1004.7 * 300 * 18 ** (g - 1)), 1e-8, 'Constant-pressure heat addition'),
      N.check('Ideal loop: ∮p dV equals the state-energy work', o.wpv, o.w, 1e-5, 'First law (trapezoidal quadrature, 4000 points)'),
      N.check('Fuel–air loop: ∮p dV equals the state-energy work', fa.wpv, fa.w, 1e-5, 'First law with variable cp'),
    ];
  },
  calibration: { params: [{ key: 'eta_mech', min: 0.75, max: 0.95 }, { key: 'diagram', min: 0.75, max: 1 }, { key: 'eta_v', min: 0.7, max: 1 }], sweep: 'alt_m', target: 'P_shaft_W', note: 'Supply dynamometer or flight-manual brake power against altitude at fixed speed and mixture; add BSFC data to separate friction from breathing.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (i.phi > 1.05) out.push({ severity: 'advise', title: 'Rich mixture wastes fuel in cruise', detail: `At equivalence ratio ${i.phi.toFixed(2)} about ${(100 * (1 - 1 / i.phi)).toFixed(0)}% of the fuel leaves unburned; BSFC is ${(o.bsfc_kg_Ws * 3.6e9).toFixed(0)} g/kWh.`, action: 'Lean to an equivalence ratio of 0.9–1.0 in cruise where cylinder temperatures allow: fuel burn and CO₂ fall almost in proportion.', basis: 'Fuel–air cycle energy balance' });
    if (o.lapse_ratio < 0.7) out.push({ severity: 'info', title: 'Large altitude power loss', detail: `Only ${(100 * o.lapse_ratio).toFixed(0)}% of sea-level power remains at ${Math.round(i.alt_m)} m.`, action: 'Turbocharging (manifold pressure ratio above 1) restores power at altitude; check the propeller match in Suite 8 and climb performance in Suite 5.', basis: 'Charge-density lapse with constant friction' });
    if (i.fuel === 'Avgas 100LL') out.push({ severity: 'advise', title: 'Leaded fuel', detail: 'Avgas 100LL contains tetra-ethyl lead; CO₂ here is ' + (o.co2_kg_s * 3600).toFixed(1) + ' kg/h.', action: 'Consider an unleaded-capable or jet-fuel (Diesel-cycle) engine: the Diesel option also cuts fuel burn through its higher compression ratio.', basis: 'Fuel database and cycle efficiency' });
    return out;
  },
};

// ---- (f) spool transient ------------------------------------------------------------------------------
/** Lumped gas-generator spool: generic compressor map + choked turbine + rotor inertia (constant cp 1005 / 1150 J/kg/K). */
function spoolModel(i) {
  const cpc = 1005, cpt = 1150, mp = cmap(i.opr, i.sm_design, i.ep_c, etaIs(i.opr, i.ep_c)), T2 = T_SL, d = mp.at(1, mp.bd);
  const T3d = T2 * (1 + (i.opr ** 0.2857 - 1) / d.eta), Wd = i.mdot, wfd = (Wd * cpt * (i.tit_K - T3d)) / (i.eta_b * i.LHV), K = (Wd * Math.sqrt(i.tit_K)) / i.opr, kt = (cpc * (T3d - T2)) / i.tit_K;
  /** Quasi-steady gas path at speed n and fuel flow wf: intersect the speed line with the choked-turbine flow function. */
  const point = (n, wf) => {
    const res = (b) => { const q = mp.at(n, b), W = Wd * q.W, T3 = T2 * (1 + (q.PR ** 0.2857 - 1) / q.eta), T4 = T3 + (wf * i.eta_b * i.LHV) / (W * cpt); return { r: (W * Math.sqrt(T4)) / q.PR - K, q, W, T3, T4 }; };
    const r0 = res(0), r1 = res(1); let s = r0, surged = false;
    if (r0.r > 0) surged = true; else if (r1.r < 0) s = r1; else s = res(N.brent((b) => res(b).r, 0, 1, 1e-10));
    const Pc = s.W * cpc * (s.T3 - T2), Pt = kt * s.W * s.T4;
    return { ...s, Pc, Pt, surged, sm: mp.sm(s.q.W, s.q.PR), net: Pt - Pc };
  };
  const steadyN = (wf) => { const v = N.findRoot((n) => point(n, wf).net, 0.35, 1.12, 22, 1e-9); return Number.isFinite(v) ? v : NaN; };
  return { mp, point, steadyN, wfd, Wd };
}
const transient = {
  id: 'transient', title: 'Spool acceleration transient and surge-margin excursion', fidelity: 'numerical',
  summary: 'Time integration of the gas-generator rotor from idle to maximum for a fuel ramp: speed, turbine temperature overshoot and how far the compressor moves towards surge.',
  equations: ['Engine spool dynamics equation J·ω·dω/dt = P_turbine − P_compressor', 'Compressor and turbine energy equations', 'Shaft power balance', 'Conservation of mass (choked turbine flow function)'],
  applicable: isGT,
  inputs: [
    { key: 'opr', label: 'Design pressure ratio', unit: '-', default: 32, min: 2, max: 70, group: 'Engine' },
    { key: 'tit_K', label: 'Design turbine inlet temperature', unit: 'K', default: 1650, min: 900, max: 2100, group: 'Engine' },
    { key: 'mdot', label: 'Design core mass flow', unit: 'kg/s', default: 55, min: 0.01, group: 'Engine' },
    { key: 'ep_c', label: 'Compressor polytropic efficiency', unit: '-', default: 0.9, min: 0.7, max: 1, group: 'Engine' },
    { key: 'eta_b', label: 'Combustion efficiency', unit: '-', default: 0.995, min: 0.8, max: 1, group: 'Engine' },
    { key: 'LHV', label: 'Fuel lower heating value', unit: 'J/kg', default: 43.15e6, min: 1e7, max: 1.3e8, group: 'Engine' },
    { key: 'sm_design', label: 'Design surge margin of the generic map', unit: '-', default: 0.2, min: 0.05, max: 0.4, group: 'Engine' },
    { key: 'J', label: 'Spool polar moment of inertia', unit: 'kg·m²', default: 8, min: 1e-5, group: 'Rotor', help: 'Default is a generic size scaling; replace with the real rotor inertia' },
    { key: 'rpm_des', label: 'Design spool speed', unit: 'rpm', default: 11000, min: 500, max: 200000, group: 'Rotor', help: 'Initial conditions are set from the idle fuel flow' },
    { key: 'idle_frac', label: 'Idle fuel flow / maximum', unit: '-', default: 0.15, min: 0.05, max: 0.6, group: 'Fuel schedule' },
    { key: 't_ramp', label: 'Fuel ramp time idle → maximum', unit: 's', default: 3, min: 0.05, max: 30, group: 'Fuel schedule', help: 'Shorter ramps accelerate faster but push the compressor towards surge and overshoot T4' },
    { key: 't_end', label: 'Simulated time', unit: 's', default: 10, min: 1, max: 60, group: 'Numerics' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 200, min: 20, max: 5000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c) {
    const i = { ...Object.fromEntries(ENGINE.map((f) => [f.key, f.default])), ...Object.fromEntries(Object.entries(gtDefaults(c)).filter(([, v]) => v !== undefined)) }, m = design(i).mdot, r = Math.sqrt(m / (1.225 * 150 * Math.PI * 0.75));
    return { opr: i.opr, tit_K: i.tit_K, mdot: m, ep_c: i.ep_c, LHV: fuelOf(i.fuel).LHV, J: 0.01 * m ** (5 / 3), rpm_des: (420 / r) * (30 / Math.PI) };
  },
  run(i) {
    const M = spoolModel(i), wd = (i.rpm_des * Math.PI) / 30, wf0 = i.idle_frac * M.wfd, nI = M.steadyN(wf0), n0 = Number.isFinite(nI) ? nI : 0.6, nS = Math.round(i.nSteps), t0 = 0.5;
    const fuel = (t) => wf0 + (M.wfd - wf0) * N.clamp((t - t0) / i.t_ramp, 0, 1), sol = N.rk4((t, y) => [M.point(y[0], fuel(t)).net / (i.J * wd * wd * y[0])], 0, [n0], i.t_end, nS);
    const pts = sol.t.map((t, k) => M.point(sol.y[k][0], fuel(t))), n = sol.y.map((y) => y[0]), sm = pts.map((p) => p.sm), jm = N.argmin(sm), T4 = pts.map((p) => p.T4), surged = pts.some((p) => p.surged);
    const tgt = n0 + 0.95 * (1 - n0); let tA = NaN; for (let k = 1; k < n.length; k++) if (n[k] >= tgt && n[k - 1] < tgt) { tA = sol.t[k - 1] + ((tgt - n[k - 1]) / (n[k] - n[k - 1])) * (sol.t[1] - sol.t[0]) - t0; break; }
    const st = N.linspace(i.idle_frac, 1, 8).map((x) => { const ns = M.steadyN(x * M.wfd); return Number.isFinite(ns) ? M.point(ns, x * M.wfd) : null; }).filter(Boolean), ds = Math.max(1, Math.ceil(n.length / 300)), pick = (a) => a.filter((_, k) => k % ds === 0 || k === a.length - 1);
    const sx = [], sy = []; for (const s of [0.7, 0.8, 0.9, 1.0]) { for (const b of N.linspace(0, 1, 12)) { const q = M.mp.at(s, b); sx.push(q.W); sy.push(q.PR); } sx.push(NaN); sy.push(NaN); }
    const warnings = ['Simplified transient: quasi-steady gas path on a generic compressor map, no volume dynamics, heat soak, tip-clearance change or bleed scheduling. Use it to compare fuel schedules, not to certify handling.'];
    if (!Number.isFinite(nI)) warnings.push('No steady idle speed was found for this idle fuel flow; the run starts from 60% speed.');
    if (surged) warnings.push('The compressor reached the generic surge line during the acceleration: the fuel ramp is too fast for this spool inertia.');
    if (!Number.isFinite(tA)) warnings.push('The spool did not reach 95% of the speed change within the simulated time; increase the simulated time.');
    return {
      kpis: [
        kp('t_accel_s', 'Acceleration time to 95% of the speed change', tA, 's', { status: tA < 5 ? 'ok' : 'warn', note: 'Counted from the start of the fuel ramp; about 5 s idle-to-95% thrust is a common handling target' }),
        kp('sm_min_trans_pct', 'Minimum surge margin during the transient', 100 * sm[jm], '%', { status: surged ? 'bad' : sm[jm] > 0.05 ? 'ok' : 'warn' }), kp('sm_excursion_pct', 'Surge-margin loss relative to steady state', 100 * (M.mp.sm(1, i.opr) - sm[jm]), '%'),
        kp('T4_peak_K', 'Peak turbine inlet temperature', N.amax(T4), 'K', { status: N.amax(T4) < i.tit_K + 50 ? 'ok' : 'warn' }), kp('T4_overshoot_K', 'T4 overshoot above the design value', Math.max(0, N.amax(T4) - i.tit_K), 'K'),
        kp('n_idle_pct', 'Idle speed', 100 * n0, '%'), kp('n_final_pct', 'Final speed', 100 * n[n.length - 1], '%'), kp('rotor_energy_J', 'Rotor kinetic energy at design speed', 0.5 * i.J * wd * wd, 'J'), kp('wf_max_kgs', 'Design fuel flow', M.wfd, 'kg/s'),
      ],
      plots: [
        { type: 'line', title: 'Spool speed and fuel flow', xlabel: 'Time [s]', ylabel: 'Fraction of design [%]', series: [{ name: 'Spool speed', x: pick(sol.t), y: pick(n).map((v) => 100 * v) }, { name: 'Fuel flow', x: pick(sol.t), y: pick(sol.t.map((t) => (100 * fuel(t)) / M.wfd)), style: 'dash' }] },
        { type: 'line', title: 'Transient trajectory on the generic compressor map', xlabel: 'Corrected flow / design [-]', ylabel: 'Pressure ratio [-]', series: [{ name: 'Speed lines 70–100%', x: sx, y: sy }, { name: 'Surge line', x: M.mp.surge.W.filter((w) => w < 1.1), y: M.mp.surge.PR.filter((_, k) => M.mp.surge.W[k] < 1.1), style: 'dash' }, { name: 'Steady operating line', x: st.map((p) => p.q.W), y: st.map((p) => p.q.PR) }, { name: 'Acceleration', x: pick(pts.map((p) => p.q.W)), y: pick(pts.map((p) => p.q.PR)) }] },
        { type: 'line', title: 'Surge margin during the acceleration', xlabel: 'Time [s]', ylabel: 'Surge margin [%]', series: [{ name: 'Surge margin', x: pick(sol.t), y: pick(sm).map((v) => 100 * v) }], annotations: [{ y: 0, label: 'Surge' }] },
        { type: 'line', title: 'Turbine inlet temperature', xlabel: 'Time [s]', ylabel: 'T4 [K]', series: [{ name: 'T4', x: pick(sol.t), y: pick(T4) }], annotations: [{ y: i.tit_K, label: 'Design T4' }] },
      ],
      warnings, models: ['Rotor inertia ODE integrated with RK4', 'Generic parametric compressor map (labelled generic)', 'Choked-turbine flow function with constant turbine temperature ratio', 'Linear fuel ramp schedule'],
      assumptions: ['Single gas-generator spool at sea-level static inlet conditions', 'Constant cp of 1005 J/kg/K (compressor) and 1150 J/kg/K (turbine)', 'No power off-take, bleed or variable geometry', 'Default inertia and speed come from a generic size scaling'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [50, 100, 200, 400], metric: 't_accel_s' },
  verify() {
    const i = Object.fromEntries(transient.inputs.map((f) => [f.key, f.default])), M = spoolModel(i), d = M.point(1, M.wfd), r = N.kv(transient.run({ ...i, t_end: 20, nSteps: 600 }));
    // free decay of a pure-inertia rotor under a constant braking power P: ω² = ω0² − 2Pt/J (exact) with the same RK4 driver
    const J = 2, w0 = 1000, Pb = 2e5, s = N.rk4((t, y) => [-Pb / (J * y[0])], 0, [w0], 2, 200);
    return [
      N.check('Design point is an equilibrium: turbine power = compressor power', d.Pt / d.Pc, 1, 1e-9, 'Shaft power balance'),
      N.check('Design fuel flow gives the design pressure ratio', d.q.PR, i.opr, 1e-7, 'Matching of map and turbine flow function'),
      N.check('Acceleration settles at 100% speed', r.n_final_pct, 100, 1e-3, 'Steady state of the spool ODE'),
      N.check('Rotor energy equation ω² = ω0² − 2Pt/J', s.y[200][0], Math.sqrt(w0 * w0 - (2 * Pb * 2) / J), 1e-9, 'Exact solution of J·ω·dω/dt = −P'),
    ];
  },
  calibration: { params: [{ key: 'J', min: 1e-4, max: 500 }, { key: 'idle_frac', min: 0.05, max: 0.4 }], sweep: 't_ramp', target: 't_accel_s', note: 'Supply measured acceleration times against fuel-ramp time from engine test-cell slam accelerations.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.sm_min_trans_pct < 5) out.push({ severity: o.sm_min_trans_pct <= 0 ? 'critical' : 'warn', title: 'Acceleration fuel schedule erodes the surge margin', detail: `Minimum transient margin ${o.sm_min_trans_pct.toFixed(1)}% with a ${i.t_ramp.toFixed(1)} s ramp.`, action: 'Lengthen the fuel ramp or add an acceleration limiter (fuel flow / compressor delivery pressure), transient bleed or variable stators.', basis: 'Surge margin on the generic map' });
    if (o.t_accel_s > 5) out.push({ severity: 'advise', title: 'Slow thrust response', detail: `95% of the speed change takes ${o.t_accel_s.toFixed(1)} s.`, action: 'Raise the idle speed (higher idle fuel flow) or shorten the ramp if surge margin allows; go-around handling normally needs most of the thrust within about 5–8 s.', basis: 'Spool energy balance' });
    if (o.T4_overshoot_K > 50) out.push({ severity: 'advise', title: 'Turbine temperature overshoot shortens hot-section life', detail: `T4 peaks ${o.T4_overshoot_K.toFixed(0)} K above the design value.`, action: 'Rate-limit the fuel ramp near maximum; pass the peak to Suite 9 for creep and thermal-fatigue life.', basis: 'Over-fuelling during acceleration' });
    return out;
  },
};

// ---- (g) electric / hybrid / fuel-cell chain ------------------------------------------------------------
const ARCH_E = ['battery-electric', 'series hybrid', 'parallel hybrid', 'fuel-cell electric'];
/** PEM cell voltage [V] at current density j [mA/cm²]: V = Eoc − A·ln(j) − r·j − m·exp(n·j). */
const cellV = (i, j) => i.fc_Eoc - i.fc_A * Math.log(Math.max(j, 1)) - i.fc_r * j - i.fc_m * Math.exp(i.fc_n * j);
/** Steady power flow from energy source to shaft for one architecture. All engines / motors together. */
function chain(i, arch) {
  const fu = fuelOf(arch === 'fuel-cell electric' ? 'Liquid hydrogen' : i.fuel), Ps = i.n_mot * i.P_mot_W * i.power_frac, H = arch === 'battery-electric' ? 1 : i.Hp;
  const par = arch === 'parallel hybrid', Pmot = par ? H * Ps : Ps, Pm_in = Pmot / (i.eta_mot * i.eta_gb), Pbus = Pm_in / i.eta_inv, Pb_bus = (par ? 1 : H) * Pbus, Pchem = Pb_bus / i.eta_batt;
  let wf = 0, Pfuel = 0, Qfc = 0, Vc = NaN, Peng = 0, lossGen = 0;
  if (arch === 'series hybrid') { const Pg = (1 - H) * Pbus; Peng = Pg / (i.eta_gen * i.eta_rect); lossGen = Peng - Pg; wf = i.bsfc_gen * Peng; }
  if (par) { Peng = ((1 - H) * Ps) / i.eta_gb; lossGen = Peng - (1 - H) * Ps; wf = i.bsfc_gen * Peng; }
  if (arch === 'fuel-cell electric') { const net = ((1 - H) * Pbus) / i.eta_dcdc, gross = net / (1 - i.fc_bop); Vc = cellV(i, i.fc_j); wf = ((gross / Vc) * 2.01588e-3) / (2 * FARADAY) / i.fc_util; Qfc = wf * fu.LHV - gross; lossGen = gross - (1 - H) * Pbus; }
  Pfuel = wf * fu.LHV;
  const lossE = (Pchem - Pb_bus) + (Pbus - Pm_in) + (Pm_in - Pmot), Psrc = Pchem + Pfuel;
  return { arch, fu, Ps, Pmot, Pbus, Pchem, Pb_bus, wf, Pfuel, Psrc, eta: Psrc > 0 ? Ps / Psrc : 0, lossE, lossBatt: Pchem - Pb_bus, lossInv: Pbus - Pm_in, lossMot: Pm_in - Pmot, lossGen, Qfc, Vc, Peng,
    heat: lossE + lossGen + (arch === 'fuel-cell electric' ? Qfc : 0), co2: wf * fu.co2, h2o: wf * fu.h2o };
}
const electric = {
  id: 'electric', title: 'Electric, hybrid-electric and fuel-cell propulsion chain', fidelity: 'analytical',
  summary: 'Power flow from battery, generator or hydrogen fuel cell through inverter, motor and propeller: chain efficiency, losses to be cooled, endurance, fuel and CO₂ for four architectures side by side.',
  equations: ['Power balance equations', 'Gas turbine–electric motor power balance', 'Battery–generator hybrid power-split equations', 'Faraday law of electrolysis (fuel-cell hydrogen consumption)', 'Fuel-cell polarisation equation'],
  inputs: [
    { key: 'arch', label: 'Architecture', type: 'select', options: ARCH_E, default: 'battery-electric', group: 'Architecture' },
    { key: 'publish', label: 'This chain is the aircraft primary propulsion', type: 'bool', default: false, group: 'Architecture', help: 'When set, its results are published as the propulsion outputs used by other suites' },
    { key: 'n_mot', label: 'Propulsion motors', unit: '', default: 2, min: 1, step: 1, discrete: true, group: 'Motors' },
    { key: 'P_mot_W', label: 'Rated shaft power per motor', unit: 'W', default: 150000, min: 1, group: 'Motors' },
    { key: 'power_frac', label: 'Power setting / rated', unit: '-', default: 1, min: 0.05, max: 1.2, group: 'Motors' },
    { key: 'Hp', label: 'Battery share of propulsive power (hybrids)', unit: '-', default: 0.3, min: 0, max: 1, group: 'Architecture', help: 'Degree of hybridisation in power; ignored for battery-electric' },
    { key: 'eta_batt', label: 'Battery discharge efficiency', unit: '-', default: 0.96, min: 0.7, max: 1, group: 'Efficiencies', help: 'Resistive loss inside the pack; Suite 19 resolves it over the mission' },
    { key: 'eta_inv', label: 'Inverter efficiency', unit: '-', default: 0.98, min: 0.8, max: 1, group: 'Efficiencies' },
    { key: 'eta_mot', label: 'Motor efficiency', unit: '-', default: 0.95, min: 0.6, max: 1, group: 'Efficiencies' },
    { key: 'eta_gb', label: 'Gearbox efficiency (1 = direct drive)', unit: '-', default: 1, min: 0.9, max: 1, group: 'Efficiencies' },
    { key: 'eta_gen', label: 'Generator efficiency', unit: '-', default: 0.95, min: 0.7, max: 1, group: 'Efficiencies' },
    { key: 'eta_rect', label: 'Rectifier efficiency', unit: '-', default: 0.98, min: 0.8, max: 1, group: 'Efficiencies' },
    { key: 'eta_dcdc', label: 'Fuel-cell DC/DC converter efficiency', unit: '-', default: 0.97, min: 0.8, max: 1, group: 'Efficiencies' },
    { key: 'eta_prop', label: 'Propeller / rotor propulsive efficiency', unit: '-', default: 0.82, min: 0.1, max: 0.95, group: 'Efficiencies' },
    { key: 'bsfc_gen', label: 'Engine BSFC (hybrids)', unit: 'kg/W/s', default: 8e-8, min: 3e-8, max: 3e-7, group: 'Hybrid engine' },
    { key: 'fuel', label: 'Engine fuel (hybrids)', type: 'select', options: FUELS, default: 'Jet A-1', group: 'Hybrid engine' },
    { key: 'batt_kWh', label: 'Battery energy', unit: 'kWh', default: 160, min: 0, group: 'Battery' },
    { key: 'usable', label: 'Usable fraction of battery energy', unit: '-', default: 0.8, min: 0.1, max: 1, group: 'Battery' },
    { key: 'fc_j', label: 'Fuel-cell current density', unit: 'mA/cm²', default: 600, min: 20, max: 1500, group: 'Fuel cell' },
    { key: 'fc_Eoc', label: 'Open-circuit voltage Eoc', unit: 'V', default: 1.031, min: 0.8, max: 1.25, group: 'Fuel cell', help: 'Representative PEM constants in the Larminie–Dicks polarisation form; calibrate against stack data' },
    { key: 'fc_A', label: 'Activation (Tafel) coefficient A', unit: 'V', default: 0.03, min: 0.005, max: 0.1, group: 'Fuel cell' },
    { key: 'fc_r', label: 'Area-specific resistance r', unit: 'kΩ·cm²', default: 2.45e-4, min: 2e-5, max: 2e-3, group: 'Fuel cell' },
    { key: 'fc_m', label: 'Mass-transport coefficient m', unit: 'V', default: 2.11e-5, min: 0, max: 1e-3, group: 'Fuel cell' },
    { key: 'fc_n', label: 'Mass-transport exponent n', unit: 'cm²/mA', default: 8e-3, min: 1e-3, max: 2e-2, group: 'Fuel cell' },
    { key: 'fc_util', label: 'Hydrogen utilisation', unit: '-', default: 0.97, min: 0.7, max: 1, group: 'Fuel cell' },
    { key: 'fc_bop', label: 'Balance-of-plant parasitic fraction', unit: '-', default: 0.12, min: 0, max: 0.4, group: 'Fuel cell', help: 'Compressor, pumps and cooling fans as a share of gross stack power' },
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 67, min: 0, max: 300, group: 'Flight condition' },
    { key: 'alt_m', label: 'Altitude', unit: 'm', default: 600, min: -500, max: 15000, group: 'Flight condition' },
    { key: 'T0_N', label: 'Static thrust per motor (0 = estimate)', unit: 'N', default: 0, min: 0, group: 'Motors' },
    { key: 'D_prop', label: 'Propeller / rotor diameter', unit: 'm', default: 2.8, min: 0.05, group: 'Motors' },
    { key: 'dT_wind', label: 'Winding temperature rise at rated power', unit: 'K', default: 90, min: 10, max: 180, group: 'Motors', help: 'From the motor data sheet; scaled with the square of the power setting (copper loss)' },
  ],
  defaults(c, up, d) {
    const p = c.prop, el = p.type === 'electric', P = el ? p.P0_W : c.systems.motor_kW > 0 ? c.systems.motor_kW * 1e3 : p.P0_W || (p.T0_N * (c.mission.cruise_V_ms || c.flight.V_ms)) / 3 || undefined;
    return { arch: el ? 'battery-electric' : 'parallel hybrid', publish: el, n_mot: p.n_eng, P_mot_W: P, eta_prop: p.eta_prop, bsfc_gen: p.bsfc_kg_Ws, fuel: fuelName(c), batt_kWh: el ? c.systems.batt_kWh : Math.max(c.systems.batt_kWh, (0.3 * (P || 0) * p.n_eng) / 2e3) || undefined,
      V: c.mission.cruise_V_ms || c.flight.V_ms, alt_m: c.mission.cruise_alt_m || c.atm.alt_m, T0_N: p.T0_N, D_prop: diaOf(c) };
  },
  run(i) {
    const r = chain(i, i.arch), all = ARCH_E.map((a) => chain(i, a)), a = isa(i.alt_m), fc = i.arch === 'fuel-cell electric', E = i.batt_kWh * 3.6e6 * i.usable;
    const Ts = staticThrust(i.T0_N, i.P_mot_W, i.D_prop, 1.225), T = i.n_mot * Math.min(Ts * (i.T0_N > 0 ? a.sigma ** 0.7 : 1), (i.eta_prop * i.P_mot_W * i.power_frac) / Math.max(i.V, 1));
    const end = r.Pchem > 0 ? E / r.Pchem : Infinity, crate = E > 0 ? r.Pchem / (i.batt_kWh * 1e3) : 0, Tw = a.T + i.dT_wind * i.power_frac ** 2;
    const js = N.linspace(20, 1400, 70), vs = js.map((j) => cellV(i, j)), Vfc = cellV(i, i.fc_j), warnings = [];
    if (r.Pchem > 0 && i.batt_kWh <= 0) warnings.push('The battery supplies power but no battery energy is defined: endurance is zero.');
    if (crate > 5) warnings.push(`Battery discharge rate ${crate.toFixed(1)} C exceeds about 5 C: expect strong voltage sag, heating and accelerated ageing (see Suite 19).`);
    if (fc && Vfc < 0.5) warnings.push('Fuel-cell operating point is in the mass-transport region (cell voltage below 0.5 V): reduce current density.');
    if (i.arch !== 'battery-electric') warnings.push('Hybrid and fuel-cell figures use a constant engine BSFC or a representative polarisation curve; they are architecture comparisons, not certified performance.');
    const out = i.publish ? { thrust_N: T, thrust_static_N: Ts, tsfc_kg_Ns: T > 0 ? r.wf / T : 0, P_shaft_W: r.Ps, eta_thermal: r.eta, eta_overall: r.eta * i.eta_prop, fuel_flow_kgs: r.wf, T4_K: Tw, EINOx_g_kg: 0, co2_kg_s: r.co2, heat_rejection_W: r.heat } : {};
    return {
      kpis: [
        kp('el_P_shaft_W', 'Shaft power (all motors and engines)', r.Ps, 'W'), kp('el_P_source_W', 'Power drawn from the energy sources', r.Psrc, 'W'), kp('el_P_batt_W', 'Battery chemical power', r.Pchem, 'W'), kp('el_P_bus_W', 'DC bus power', r.Pbus, 'W'),
        kp('el_eta_chain', 'Source-to-shaft efficiency', r.eta, '-', { status: r.eta > 0.8 || i.arch !== 'battery-electric' ? 'ok' : 'warn' }), kp('el_eta_overall', 'Source-to-thrust efficiency', r.eta * i.eta_prop, '-'),
        kp('el_thrust_N', 'Thrust at the flight point (all motors)', T, 'N'), kp('el_thrust_static_N', 'Static thrust per motor', Ts, 'N', { note: i.T0_N > 0 ? 'Rating' : 'Momentum-theory estimate' }),
        kp('el_fuel_flow_kgs', fc ? 'Hydrogen flow' : 'Fuel flow', r.wf, 'kg/s'), kp('el_co2_kg_s', 'CO₂ emission at the exhaust', r.co2, 'kg/s'), kp('el_h2o_kg_s', 'Water emission', r.h2o, 'kg/s'),
        kp('el_heat_W', 'Heat to be rejected by the cooling system', r.heat, 'W', { note: 'Battery, inverter, motor, generator and fuel-cell losses' }),
        kp('el_endurance_min', 'Battery endurance at this power', end / 60, 'min', { status: end > 600 ? 'ok' : 'warn' }), kp('el_c_rate', 'Battery discharge rate', crate, '1/h', { status: crate <= 5 ? 'ok' : 'warn' }),
        kp('el_T_winding_K', 'Motor winding temperature estimate', Tw, 'K', { status: Tw < 453 ? 'ok' : 'warn', note: 'Class H insulation limit is 180 °C (453 K)' }),
        kp('fc_V_cell', 'Fuel-cell voltage at the operating current density', Vfc, 'V'), kp('fc_eta_lhv', 'Fuel-cell stack efficiency (LHV)', (Vfc / 1.254) * i.fc_util, '-', { note: 'Cell voltage / 1.254 V × utilisation' }), kp('fc_power_density_W_cm2', 'Fuel-cell power density', (Vfc * i.fc_j) / 1e3, 'W/cm²'),
      ],
      plots: [
        { type: 'bar', title: 'Power flow from source to shaft', ylabel: 'Power [kW]', categories: [i.arch], stacked: true, series: [{ name: 'Shaft power', y: [r.Ps / 1e3] }, { name: 'Motor and gearbox loss', y: [r.lossMot / 1e3] }, { name: 'Inverter loss', y: [r.lossInv / 1e3] }, { name: 'Battery loss', y: [r.lossBatt / 1e3] }, { name: 'Generator / converter loss', y: [r.lossGen / 1e3] }, { name: 'Engine or fuel-cell heat', y: [Math.max(0, r.Psrc - r.Ps - r.lossE - r.lossGen) / 1e3] }] },
        { type: 'bar', title: 'Source power for the same shaft power', ylabel: 'Power [kW]', categories: ARCH_E, stacked: true, series: [{ name: 'Battery', y: all.map((x) => x.Pchem / 1e3) }, { name: 'Fuel (LHV)', y: all.map((x) => x.Pfuel / 1e3) }] },
        { type: 'bar', title: 'Exhaust CO₂ by architecture', ylabel: 'CO₂ [kg/h]', categories: ARCH_E, series: [{ name: 'CO₂ at the aircraft', y: all.map((x) => x.co2 * 3600) }] },
        { type: 'line', title: 'Fuel-cell polarisation curve', xlabel: 'Current density [mA/cm²]', ylabel: 'Cell voltage [V]', series: [{ name: 'Cell voltage', x: js, y: vs }], annotations: [{ x: i.fc_j, label: 'Operating point' }] },
      ],
      tables: [{ title: 'Architecture comparison at the same shaft power', columns: ['Architecture', 'Source power [kW]', 'Efficiency [-]', 'Fuel [kg/h]', 'CO₂ [kg/h]', 'Heat to reject [kW]'], rows: all.map((x) => [x.arch, +(x.Psrc / 1e3).toFixed(1), +x.eta.toFixed(3), +(x.wf * 3600).toFixed(2), +(x.co2 * 3600).toFixed(1), +(x.heat / 1e3).toFixed(1)]) }],
      outputs: out,
      warnings, models: ['Steady efficiency-chain power balance', 'Series / parallel power split by degree of hybridisation', 'PEM polarisation curve (Larminie–Dicks form, representative constants)', 'Faraday-law hydrogen consumption (exact)'],
      assumptions: ['Constant component efficiencies at the operating point; no voltage sag or thermal derating (resolved in Suite 19)', 'Hybrid engine at constant BSFC', 'Grid or hydrogen production emissions are outside the aircraft boundary and not counted', 'Winding temperature scales with the square of the power setting from the rated rise'],
    };
  },
  verify() {
    const i = Object.fromEntries(electric.inputs.map((f) => [f.key, f.default])), b = chain(i, 'battery-electric'), fcx = chain({ ...i, Hp: 0, fc_bop: 0, eta_dcdc: 1, fc_util: 1 }, 'fuel-cell electric'), s = chain(i, 'series hybrid');
    return [
      N.check('Battery chain efficiency = product of component efficiencies', b.eta, i.eta_batt * i.eta_inv * i.eta_mot * i.eta_gb, 1e-12, 'Power balance'),
      N.check('Energy conservation: source = shaft + losses', b.Ps + b.lossE, b.Psrc, 1e-12, 'First law'),
      N.check('Faraday law: 1 A through one cell consumes M/(2F) kg/s of hydrogen', (fcx.wf * fcx.Vc) / fcx.Pbus, 2.01588e-3 / (2 * 96485.332), 1e-9, 'Faraday constant 96485.332 C/mol, two electrons per H₂'),
      N.check('Series hybrid bus balance: battery + generator = bus', s.Pb_bus + s.Peng * i.eta_gen * i.eta_rect, s.Pbus, 1e-12, 'Kirchhoff power balance at the DC bus'),
    ];
  },
  calibration: { params: [{ key: 'fc_r', min: 5e-5, max: 1e-3 }, { key: 'fc_A', min: 0.01, max: 0.08 }, { key: 'fc_Eoc', min: 0.9, max: 1.2 }], sweep: 'fc_j', target: 'fc_V_cell', note: 'Supply a measured stack polarisation curve (cell voltage against current density).' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (i.arch === 'battery-electric' && o.el_endurance_min < 30) out.push({ severity: 'warn', title: 'Short battery endurance at this power', detail: `${o.el_endurance_min.toFixed(1)} min of usable energy at ${(o.el_P_batt_W / 1e3).toFixed(0)} kW.`, action: 'Lower the power demand (disk loading, drag), raise battery energy, or compare a series hybrid or fuel-cell range extender in the architecture table.', basis: 'Usable battery energy / chemical power' });
    if (o.el_c_rate > 3) out.push({ severity: 'advise', title: 'High discharge rate shortens battery life', detail: `${o.el_c_rate.toFixed(1)} C continuous.`, action: 'Size the pack for no more than about 3 C in sustained flight phases; check ageing and temperature in Suite 19.', basis: 'Cell cycle-life sensitivity to C-rate' });
    out.push({ severity: 'info', title: 'Zero exhaust CO₂ is not zero life-cycle CO₂', detail: `Battery and fuel-cell chains emit ${i.arch === 'battery-electric' ? 'nothing' : (o.el_co2_kg_s * 3600).toFixed(1) + ' kg/h CO₂'} at the aircraft, but electricity and hydrogen production carry their own footprint.`, action: 'Use the energy figures here with the grid or hydrogen carbon intensity in Suite 26 to compare options on a life-cycle basis.', basis: 'System boundary of this analysis' });
    if (o.el_heat_W > 0.15 * o.el_P_shaft_W) out.push({ severity: 'advise', title: 'Large low-grade heat load', detail: `${(o.el_heat_W / 1e3).toFixed(0)} kW must be rejected at low temperature.`, action: 'Size the thermal-management system in Suite 12; fuel-cell heat at 60–80 °C needs large radiators and adds cooling drag.', basis: 'Chain losses' });
    return out;
  },
};

export default {
  id: 'propulsion', n: 7,
  tagline: 'How much thrust or power the powerplant gives, what it burns and emits, and how it behaves away from its design point.',
  analyses: [cyc, map, turbo, combustor, piston, transient, electric],
  consumes: [],
  provides: [
    { key: 'thrust_N', label: 'Thrust at the flight point (all engines)', unit: 'N' }, { key: 'thrust_static_N', label: 'Static thrust per engine', unit: 'N' }, { key: 'tsfc_kg_Ns', label: 'Thrust-specific fuel consumption', unit: 'kg/N/s' },
    { key: 'P_shaft_W', label: 'Shaft power (all engines)', unit: 'W' }, { key: 'eta_thermal', label: 'Thermal (source-to-shaft) efficiency', unit: '-' }, { key: 'eta_overall', label: 'Overall efficiency', unit: '-' },
    { key: 'fuel_flow_kgs', label: 'Fuel flow (all engines)', unit: 'kg/s' }, { key: 'T4_K', label: 'Hot-section temperature (turbine inlet / peak cylinder gas / motor winding)', unit: 'K' }, { key: 'EINOx_g_kg', label: 'NOx emission index', unit: 'g/kg' },
    { key: 'co2_kg_s', label: 'CO₂ emission rate', unit: 'kg/s' }, { key: 'heat_rejection_W', label: 'Heat to the cooling / oil system', unit: 'W' },
  ],
  handoff: [
    { model: 'Component-map engine deck (full compressor and turbine maps, variable geometry, bleed schedules)', why: 'Needs proprietary rig-measured maps; the built-in map is a generic parametric shape', tool: 'Engine performance deck / gas-turbine simulation code (NPSS-class)' },
    { model: 'Three-dimensional turbomachinery CFD and coupled zero-dimensional cycle–CFD', why: 'Blade-row resolved RANS/LES is far beyond in-browser cost', tool: 'Turbomachinery CFD solver with mixing-plane or sliding interfaces' },
    { model: 'Detailed chemical kinetics, reacting-flow combustor CFD, soot/CO/UHC and equilibrium dissociation', why: 'Stiff multi-species chemistry and spray modelling are not solved; NOx is an empirical T3–P3 correlation', tool: 'Chemical-kinetics solver / reacting-flow CFD; chemical-equilibrium code' },
    { model: 'Turbine blade cooling and conjugate heat transfer', why: 'Only a cooling-air fraction is modelled; metal temperatures need internal-passage geometry', tool: 'Conjugate heat-transfer CFD / Suite 12 for first estimates' },
    { model: 'Inlet distortion and engine–airframe aerodynamic interaction', why: 'Requires the installed flow field and parallel-compressor or body-force models', tool: 'Coupled airframe–inlet CFD with distortion descriptors' },
    { model: 'Full transient engine model (volume dynamics, heat soak, tip clearance, starting)', why: 'The transient here is a single-spool quasi-steady gas path on a generic map', tool: 'Transient engine simulation code validated on test-cell data' },
    { model: 'Crank-angle-resolved piston engine simulation (Wiebe burn, knock, gas exchange, turbocharger matching)', why: 'The piston model is an ideal fuel–air loop with empirical factors', tool: 'One-dimensional engine gas-dynamics code' },
    { model: 'Engine structural models and engine health state estimation', why: 'Structural life is handled in Suites 2, 9 and 10; health estimation needs fleet sensor data', tool: 'FE rotor-dynamics codes; gas-path-analysis / Kalman-filter health monitoring tools' },
    { model: 'Electrochemical fuel-cell and battery models', why: 'The fuel cell is a static polarisation curve; batteries are resolved in Suite 19', tool: 'Suite 19; electrochemical stack simulation' },
  ],
};
