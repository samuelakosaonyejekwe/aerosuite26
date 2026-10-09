// Suite 20 — Fuel Systems and Environmental Control.
// Fuel feed hydraulics and vapour margins, fuel quantity and centre-of-gravity travel, tank sloshing, fuel cold-soak and
// tank inerting, cabin heat/ventilation/humidity balance, pressurisation and decompression, and the air-conditioning pack.

import * as N from '../core/numerics.js';
import { isa, G0, R_AIR, GAMMA, CP_AIR, pSat, pressureAltitude } from '../core/atmosphere.js';
import { FLUIDS } from '../data/materials.js';

// ---- shared helpers --------------------------------------------------------------------------
const num = (key, label, unit, def, min, max, group, help, x) => ({ key, label, unit, default: def, min, max, group, ...(help ? { help } : {}), ...x });
const sel = (key, label, options, def, group, help) => ({ key, label, type: 'select', options, default: def, group, ...(help ? { help } : {}) });
const kp = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const thin = (a, n = 300) => { if (a.length <= n) return a; const s = (a.length - 1) / (n - 1); return N.range(n, (j) => a[Math.round(j * s)]); };
const g = GAMMA, KG = (g - 1) / g;
// Approximate fuel volatility and freezing data: true vapour pressure p_v = pv38·exp(−B(1/T − 1/311 K)).
// pv38 and B are representative figures (they vary with batch and weathering); freezing points are specification maxima.
const FUELX = { 'Jet A-1': { pv38: 1.0e3, B: 4500, Tfr: 226.15 }, 'Avgas 100LL': { pv38: 40e3, B: 3400, Tfr: 215.15 }, 'SAF (HEFA-SPK)': { pv38: 1.0e3, B: 4500, Tfr: 233.15 } };
const FUELS = Object.keys(FUELX), fuelOf = (name) => ({ ...(FLUIDS[name] || FLUIDS['Jet A-1']), ...(FUELX[name] || FUELX['Jet A-1']) });
const pVap = (f, T) => f.pv38 * Math.exp(-f.B * (1 / T - 1 / 311.15));
const hasFuel = (c) => (!(c.mass.fuel_kg > 0) ? 'This aircraft carries no liquid fuel (battery-electric). The cabin and equipment-bay air analysis applies instead.' : c.prop.fuel === 'Liquid hydrogen' ? 'Cryogenic hydrogen storage needs dedicated boil-off and two-phase models that are not solved here.' : true);
const fuelName = (c) => (FUELX[c.prop.fuel] ? c.prop.fuel : undefined);
const siteT = (c) => (c.site?.T_C ?? 15) + 273.15;
const cruise = (c) => { const alt = c.mission.cruise_alt_m || c.atm.alt_m || 0, a = isa(alt, c.atm.dISA_K || 0), V = c.mission.cruise_V_ms || c.flight.V_ms || 30; return { alt, a, V, M: V / a.a }; };
const fuelFlowEst = (c, up, d) => up.performance?.fuel_flow_cruise_kgs ?? up.propulsion?.fuel_flow_kgs ?? (d.T_total > 0 ? 0.22 * c.prop.tsfc_kg_Ns * d.T_total : 0.65 * c.prop.bsfc_kg_Ws * d.P_total);
const fDarcy = (Re, rr = 0) => { if (!(Re > 1e-12)) return 0; const tur = 0.25 / Math.log10(rr / 3.7 + 5.74 / Math.max(Re, 2000) ** 0.9) ** 2; if (Re <= 2000) return 64 / Re; if (Re >= 4000) return tur; const w = (Re - 2000) / 2000; return (1 - w) * (64 / Re) + w * tur; };
const pipeDp = (Q, D, L, rho, mu, rough = 0, K = 0) => { const v = Math.abs(Q) / ((Math.PI * D * D) / 4); return N.sign(Q) * (fDarcy((rho * v * D) / mu, rough / D) * (L / D) + K) * 0.5 * rho * v * v; };
/** Compressible orifice mass flow [kg/s] from (p0, T0) to back pressure pb through effective area CdA. */
const orifice = (CdA, p0, T0, pb) => { if (pb >= p0) return 0; const rc = (2 / (g + 1)) ** (g / (g - 1)), r = Math.max(pb / p0, rc); return CdA * p0 * Math.sqrt(((2 * g) / ((g - 1) * R_AIR * T0)) * (r ** (2 / g) - r ** ((g + 1) / g))); };
const humRatio = (pv, p) => (0.622 * pv) / Math.max(p - pv, 1);
const dewPoint = (pv) => { const L = Math.log(Math.max(pv, 1e-3) / 610.94); return (273.15 * 17.625 - 30.11 * L) / (17.625 - L); };

// ---- 1. fuel feed system ----------------------------------------------------------------------
function feedState(i, alt = i.alt_m, Q = i.mdot / fuelOf(i.fuel).rho, pumpOn = true) {
  const f = fuelOf(i.fuel), rho = f.rho, pv = pVap(f, i.T_fuel), pt = isa(alt).p + i.p_ullage, A = (Math.PI * i.D ** 2) / 4, v = Q / A;
  const pin = pt + rho * G0 * i.h_fuel - 0.5 * rho * v * v * i.K_in, dpp = pumpOn ? i.dp0 * (1 - (Q / i.Q_max) ** 2) : -0.5 * rho * v * v * i.K_pump_off;
  const dpl = pipeDp(Q, i.D, i.L, rho, f.mu, i.rough, i.K_fit), pe = pin + dpp - dpl - rho * G0 * i.dz;
  return { f, rho, pv, pt, pin, dpp, dpl, pe, v, Re: (rho * v * i.D) / f.mu, npsh: (pin - pv) / (rho * G0), margin: pe - pv - i.p_req };
}
/** Largest flow that still meets the engine inlet pressure requirement. */
const feedMaxQ = (i, alt, pumpOn) => { const hi = pumpOn ? i.Q_max : 50 * i.Q_max, fn = (Q) => feedState(i, alt, Q, pumpOn).margin; return fn(1e-12) <= 0 ? 0 : fn(hi) >= 0 ? hi : N.brent(fn, 1e-12, hi, 1e-12 * i.Q_max); };
const feed = {
  id: 'feed', title: 'Engine fuel feed: pump operating point and vapour margins', fidelity: 'numerical',
  summary: 'Pressure delivered to the engine from tank, boost pump, line losses and elevation, the pump suction head against fuel vapour pressure at altitude, and whether gravity alone can feed the engine if the pump fails.',
  equations: ['Darcy–Weisbach equation', 'Bernoulli equation', 'Hydrostatic pressure equations', 'Vapour–liquid equilibrium equations', 'Conservation of mass, momentum and energy equations', 'Incompressible and compressible fluid flow equations'],
  applicable: hasFuel,
  inputs: [
    sel('fuel', 'Fuel', FUELS, FUELS[0], 'Fuel'), num('T_fuel', 'Fuel temperature', 'K', 303, 200, 360, 'Fuel', 'Hot fuel after a ground soak is the critical case for vapour margins'),
    num('mdot', 'Fuel flow demanded per feed line', 'kg/s', 1.2, 1e-6, 50, 'Demand', 'Take-off flow of one engine'), num('alt_m', 'Altitude', 'm', 10668, -500, 20000, 'Demand'),
    num('p_ullage', 'Tank ullage pressure above ambient', 'Pa', 0, -5000, 2e5, 'Tank', '0 for an open-vented tank'), num('h_fuel', 'Fuel head above the pump inlet', 'm', 0.3, 0, 10, 'Tank'),
    num('dz', 'Engine inlet height above the pump', 'm', -0.4, -20, 20, 'Line', 'Negative when the engine is below the tank'),
    num('D', 'Feed line bore', 'm', 0.032, 0.002, 0.3, 'Line'), num('L', 'Feed line length', 'm', 10, 0.05, 200, 'Line'), num('rough', 'Wall roughness', 'm', 5e-6, 0, 1e-3, 'Line'),
    num('K_fit', 'Fittings, valves and filter loss coefficient', '-', 8, 0, 500, 'Line'), num('K_in', 'Pump inlet loss coefficient', '-', 1, 0, 50, 'Line'),
    num('dp0', 'Boost pump shut-off pressure rise', 'Pa', 150e3, 0, 5e6, 'Pump', 'Centrifugal tank boost pumps ≈ 100–250 kPa'), num('Q_max', 'Boost pump flow at zero pressure rise', 'm³/s', 5e-3, 1e-8, 1, 'Pump'),
    num('eta_pump', 'Pump overall efficiency', '-', 0.5, 0.05, 0.9, 'Pump', 'Small centrifugal fuel pumps ≈ 0.3–0.6'), num('K_pump_off', 'Loss coefficient of the stopped pump', '-', 6, 0, 200, 'Pump', 'For the gravity-feed case'),
    num('npsh_req', 'Pump NPSH required', 'm', 1, 0, 50, 'Requirement'), num('p_req', 'Required engine inlet pressure above vapour pressure', 'Pa', 35e3, 0, 5e5, 'Requirement', 'Engine manufacturers typically ask for about 35 kPa (5 psi) above true vapour pressure'),
  ],
  defaults: (c, up, d) => {
    const jet = d.T_total > 0 && !(c.prop.P0_W > 0), md = Math.max(1e-6, jet ? 0.6 * c.prop.tsfc_kg_Ns * c.prop.T0_N : c.prop.bsfc_kg_Ws * c.prop.P0_W), f = fuelOf(c.prop.fuel), Q = md / f.rho, turbine = c.prop.type !== 'piston';
    const heli = c.meta.type === 'helicopter', D = Math.max(0.004, Math.sqrt((4 * Q) / (Math.PI * 1.8)));
    return { fuel: fuelName(c), T_fuel: Math.max(siteT(c), 288.15) + 10, mdot: md, alt_m: cruise(c).alt, D, L: Math.max(0.3, 0.3 * (c.fuselage.len_m || 3)), Q_max: 3 * Q, dp0: turbine ? 150e3 : 35e3, p_req: turbine ? 35e3 : 5e3, h_fuel: N.clamp(0.4 * c.wing.tc * (d.c_root || 0.5), 0.03, 0.6), dz: heli ? 1.5 : c.prop.type === 'piston' ? -0.5 : c.prop.type === 'turboprop' ? 0.3 : -0.4, npsh_req: turbine ? 1 : 0.3 };
  },
  run(i) {
    const s = feedState(i), Q = i.mdot / s.rho, Ppump = (Q * Math.max(s.dpp, 0)) / i.eta_pump, Qmax = feedMaxQ(i, i.alt_m, true), Qg = feedMaxQ(i, i.alt_m, false), nm = s.npsh - i.npsh_req, warnings = [];
    const alts = N.linspace(0, Math.max(3000, 1.3 * i.alt_m), 40), sw = alts.map((h) => feedState(i, h)), sg = alts.map((h) => feedState(i, h, Q, false));
    const hLim = N.findRoot((h) => feedState(i, h, Q, false).margin, 0, 20000, 80), Qs = N.linspace(0, i.Q_max, 50);
    if (s.margin < 0) warnings.push('Engine inlet pressure is below the required margin above vapour pressure: risk of vapour lock or engine-pump cavitation.');
    if (nm < 0) warnings.push('Boost-pump suction head is below NPSH required: the pump will cavitate at this altitude and fuel temperature.');
    if (s.v > 3) warnings.push(`Line velocity ${s.v.toFixed(1)} m/s exceeds the usual 2–3 m/s guideline for fuel lines (pressure loss and static-charge generation).`);
    if (Q > i.Q_max) warnings.push('Demanded flow exceeds the pump zero-head flow: the pump adds no pressure.');
    return {
      kpis: [
        kp('fuel_pump_power_W', 'Boost pump electrical power', Ppump, 'W'), kp('p_engine_inlet_Pa', 'Engine inlet pressure (absolute)', s.pe, 'Pa'),
        kp('inlet_margin_Pa', 'Margin above the engine inlet requirement', s.margin, 'Pa', s.margin > 10e3 ? 'ok' : s.margin > 0 ? 'warn' : 'bad'),
        kp('npsh_available_m', 'Pump NPSH available', s.npsh, 'm', nm > 0.5 ? 'ok' : nm > 0 ? 'warn' : 'bad', `Required ${i.npsh_req.toFixed(1)} m`),
        kp('p_vapour_Pa', 'Fuel true vapour pressure', s.pv, 'Pa'), kp('pump_dp_Pa', 'Pump pressure rise at the demand', s.dpp, 'Pa'), kp('line_dp_Pa', 'Feed line pressure loss', s.dpl, 'Pa'),
        kp('line_velocity_ms', 'Line velocity', s.v, 'm/s', s.v <= 3 ? 'ok' : 'warn'), kp('Re_line', 'Line Reynolds number', s.Re, '-'),
        kp('max_flow_ratio', 'Deliverable flow / demand (pump on)', Qmax / Q, '-', Qmax >= 1.25 * Q ? 'ok' : Qmax >= Q ? 'warn' : 'bad'),
        kp('gravity_flow_ratio', 'Deliverable flow / demand (pump failed)', Qg / Q, '-', Qg >= Q ? 'ok' : 'warn', 'Gravity or suction feed capability'),
        kp('suction_feed_ceiling_m', 'Altitude limit with the boost pump failed', Number.isFinite(hLim) ? hLim : sg[0].margin > 0 ? 20000 : 0, 'm'),
      ],
      plots: [
        { type: 'line', title: 'Engine inlet margin versus altitude', xlabel: 'Altitude [m]', ylabel: 'Pressure above the inlet requirement [kPa]', series: [{ name: 'Boost pump on', x: alts, y: sw.map((v) => v.margin / 1e3) }, { name: 'Boost pump failed', x: alts, y: sg.map((v) => v.margin / 1e3) }], annotations: [{ y: 0, label: 'Limit' }, { x: i.alt_m, label: 'Design altitude' }] },
        { type: 'line', title: 'Pump characteristic and system demand', xlabel: 'Flow [L/min]', ylabel: 'Pressure rise [kPa]', series: [{ name: 'Boost pump', x: Qs.map((q) => q * 6e4), y: Qs.map((q) => (i.dp0 * (1 - (q / i.Q_max) ** 2)) / 1e3) }, { name: 'Needed from the pump', x: Qs.map((q) => q * 6e4), y: Qs.map((q) => { const z = feedState(i, i.alt_m, q); return (z.dpp - z.margin) / 1e3; }) }], annotations: [{ x: Q * 6e4, label: 'Demand' }] },
        { type: 'bar', title: 'Pressure build-up from tank to engine', ylabel: 'Pressure [kPa]', categories: ['Tank surface', 'Fuel head', 'Pump rise', 'Line loss', 'Elevation', 'Engine inlet', 'Vapour pressure + requirement'], series: [{ name: 'Pressure', y: [s.pt, s.rho * G0 * i.h_fuel, s.dpp, -s.dpl, -s.rho * G0 * i.dz, s.pe, s.pv + i.p_req].map((v) => v / 1e3) }] },
      ],
      warnings, models: ['Quadratic centrifugal pump characteristic', 'Darcy–Weisbach line with minor losses', 'Hydrostatic head and vented-tank pressure from the standard atmosphere', 'Exponential (Clausius–Clapeyron type) fuel vapour-pressure model — approximate'],
      assumptions: ['Steady single-phase flow; no dissolved-air evolution', 'Fuel vapour pressure from representative constants; real fuels vary with batch and weathering', 'Tank vented to ambient static pressure plus the stated ullage pressure', '1 g level flight (no manoeuvre or attitude head changes)', 'Pump characteristic, loss coefficients, NPSH and the engine inlet pressure requirement are typical values; the default case combines take-off fuel flow with cruise altitude and warm fuel, which is conservative'],
    };
  },
  verify() {
    const b = { fuel: FUELS[0], T_fuel: 288.15, mdot: 0.5, alt_m: 0, p_ullage: 0, h_fuel: 1, dz: 0, D: 0.02, L: 1e-9, rough: 0, K_fit: 0, K_in: 0, dp0: 1e5, Q_max: 2e-3, eta_pump: 0.5, K_pump_off: 1, npsh_req: 1, p_req: 0 }, f = fuelOf(FUELS[0]);
    // gravity discharge to ambient through an exit loss K = 1: v = sqrt(2·g·h)
    const tor = { ...b, p_req: 101325 - pVap(f, 288.15) }, Qg = feedMaxQ(tor, 0, false), s = feedState({ ...b, L: 5, D: 0.01 }, 0, 1e-5);
    const lin = { ...b, p_req: 101325 + f.rho * G0 * 1 + 0.5e5 - pVap(f, 288.15) }; // requirement placed where the pump must supply exactly 50 kPa
    return [
      N.check('Gravity discharge velocity (Torricelli)', Qg / ((Math.PI * 0.02 ** 2) / 4), Math.sqrt(2 * G0 * 1), 1e-8, 'Torricelli’s law'),
      N.check('Laminar feed-line loss (Hagen–Poiseuille)', s.dpl, (128 * f.mu * 5 * 1e-5) / (Math.PI * 0.01 ** 4), 1e-9, 'Hagen–Poiseuille'),
      N.check('Pump/system intersection for a quadratic pump curve', feedMaxQ(lin, 0, true), 2e-3 * Math.sqrt(0.5), 1e-7, 'Δp = Δp0(1 − (Q/Qmax)²)'),
      N.check('NPSH available = (p_tank + ρgh − p_v)/(ρg)', feedState(b, 0, 0).npsh, (101325 - pVap(f, 288.15)) / (f.rho * G0) + 1, 1e-12, 'Definition'),
    ];
  },
  calibration: { params: [{ key: 'K_fit', min: 0, max: 200 }, { key: 'dp0', min: 1e3, max: 1e6 }, { key: 'Q_max', min: 1e-7, max: 0.1 }], sweep: 'mdot', target: 'p_engine_inlet_Pa', note: 'Supply measured engine inlet pressure against fuel flow from a fuel-system rig or ground run.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.inlet_margin_Pa < 0 || o.npsh_available_m < i.npsh_req) out.push({ severity: 'critical', title: 'Insufficient fuel pressure margin', detail: `Engine inlet margin ${(o.inlet_margin_Pa / 1e3).toFixed(1)} kPa; NPSH available ${o.npsh_available_m.toFixed(2)} m.`, action: 'Raise boost-pump pressure, enlarge the feed line, pressurise the tank, or restrict altitude with hot fuel.', basis: 'Vapour-pressure margin at the engine and pump inlets' });
    if (o.gravity_flow_ratio < 1) out.push({ severity: 'advise', title: 'No gravity or suction feed at this condition', detail: `With the pump failed only ${(100 * o.gravity_flow_ratio).toFixed(0)}% of the demanded flow can be delivered; the limit is about ${o.suction_feed_ceiling_m.toFixed(0)} m.`, action: 'Provide a second pump on an independent supply, or define a suction-feed altitude limitation; carry the pump failure into Suite 22.', basis: 'Feed availability after a single failure' });
    if (i.fuel === FUELS[1]) out.push({ severity: 'info', title: 'Volatile fuel', detail: `Vapour pressure is ${(o.p_vapour_Pa / 1e3).toFixed(0)} kPa at this temperature, far above kerosene.`, action: 'Keep fuel lines away from heat sources and check hot-day climb; unleaded or kerosene-type fuels reduce vapour and lead emissions.', basis: 'Fuel volatility' });
    return out;
  },
};

// ---- 2. fuel quantity and CG travel -----------------------------------------------------------
const LAYOUT = ['Wing tanks + centre tank', 'Fuselage tanks (forward + aft)'], SEQ = ['Centre (or forward) tank first', 'Wing (or aft) tanks first', 'All tanks together'];
/** Fuel centroid of one wing tank against fill fraction: fuel collects inboard (dihedral), cross-section ∝ chord². */
function wingTank(i) {
  const s = i.b / 2, y1 = i.eta1 * s, y2 = i.eta2 * s, n = 80, ys = N.linspace(y1, y2, n), c = (y) => i.c_root * (1 - (1 - i.taper) * (y / s)), tanL = Math.tan(N.rad(i.sweep_deg));
  const x = (y) => i.x_wing + y * tanL + 0.15 * c(y), a = ys.map((y) => c(y) ** 2), V = N.cumtrapz(ys, a), My = N.cumtrapz(ys, a.map((v, j) => v * ys[j])), Mx = N.cumtrapz(ys, a.map((v, j) => v * x(ys[j])));
  return (fill) => { const v = N.clamp(fill, 0, 1) * V[n - 1]; if (v <= 1e-12 * V[n - 1]) return { x: x(y1), y: y1 }; return { x: N.interp1(V, Mx, v) / v, y: N.interp1(V, My, v) / v }; };
}
function cgBurn(i) {
  const wing = i.layout === LAYOUT[0], tank = wing ? wingTank(i) : null, mA0 = i.fuel_kg * i.frac_A, mB0 = i.fuel_kg - mA0, n = Math.max(10, Math.round(i.nSteps)), burnTot = i.fuel_kg * (1 - i.reserve_frac);
  const t = [], mA = [], mL = [], mR = [], xcg = [], ycg = [], roll = [], xf = [];
  for (let k = 0; k <= n; k++) {
    const used = (burnTot * k) / n; let a, b;
    if (i.seq === SEQ[0]) { a = Math.max(0, mA0 - used); b = mB0 - Math.max(0, used - mA0); } else if (i.seq === SEQ[1]) { b = Math.max(0, mB0 - used); a = mA0 - Math.max(0, used - mB0); } else { a = mA0 * (1 - used / i.fuel_kg); b = mB0 * (1 - used / i.fuel_kg); }
    // group B is split left/right (wing layout) with a flow mismatch between the sides
    const ub = mB0 - b, l = wing ? Math.max(0, mB0 / 2 - (ub / 2) * (1 + i.imbalance)) : b / 2, r = wing ? b - l : b / 2;
    let xA = i.x_A, xL = i.x_B, xR = i.x_B, yL = 0, yR = 0;
    if (wing) { const cl = tank(l / Math.max(mB0 / 2, 1e-12)), cr = tank(r / Math.max(mB0 / 2, 1e-12)); xL = cl.x; xR = cr.x; yL = -cl.y; yR = cr.y; }
    const mf = a + l + r, M = i.zfw_kg + mf, mx = a * xA + l * xL + r * xR, my = l * yL + r * yR;
    t.push(used / i.burn_rate); mA.push(a); mL.push(l); mR.push(r); xcg.push(mx / M); ycg.push(my / M); roll.push(my * G0); xf.push(mf > 1e-9 ? mx / mf : 0);
  }
  return { t, mA, mL, mR, xcg, ycg, roll, xf, wing };
}
const cg = {
  id: 'cg', title: 'Fuel quantity, centre-of-gravity travel and lateral balance', fidelity: 'reduced-order',
  summary: 'Tracks tank quantities through the flight for a chosen burn sequence and shows how far the aircraft centre of gravity moves, fore–aft and laterally, as the fuel is used.',
  equations: ['Conservation of mass, momentum and energy equations', 'Hydrostatic pressure equations'],
  applicable: hasFuel,
  inputs: [
    sel('layout', 'Tank layout', LAYOUT, LAYOUT[0], 'Tanks'), num('fuel_kg', 'Fuel on board at start', 'kg', 18800, 0.01, 5e5, 'Tanks'),
    num('frac_A', 'Share of fuel in the centre (or forward) tank', '-', 0.3, 0, 1, 'Tanks'), sel('seq', 'Burn sequence', SEQ, SEQ[0], 'Tanks', 'Centre tank first keeps wing fuel for bending relief'),
    num('x_A', 'Centre/forward tank position aft of the zero-fuel CG', 'm', -0.3, -50, 50, 'Tanks'), num('x_B', 'Aft tank position aft of the zero-fuel CG', 'm', 1.5, -50, 50, 'Tanks', 'Fuselage layout only'),
    num('x_wing', 'Wing root quarter-chord aft of the zero-fuel CG', 'm', -0.5, -50, 50, 'Wing tank geometry'), num('b', 'Wing span', 'm', 34.1, 0.1, 100, 'Wing tank geometry'), num('c_root', 'Root chord', 'm', 5.8, 0.02, 30, 'Wing tank geometry'),
    num('taper', 'Taper ratio', '-', 0.24, 0.05, 1, 'Wing tank geometry'), num('sweep_deg', 'Quarter-chord sweep', 'deg', 25, -10, 60, 'Wing tank geometry'),
    num('eta1', 'Tank inboard end / semi-span', '-', 0.1, 0, 0.9, 'Wing tank geometry'), num('eta2', 'Tank outboard end / semi-span', '-', 0.7, 0.05, 1, 'Wing tank geometry'),
    num('zfw_kg', 'Zero-fuel mass', 'kg', 59200, 0.1, 6e5, 'Aircraft'), num('mac', 'Mean aerodynamic chord', 'm', 4.2, 0.02, 30, 'Aircraft'),
    num('burn_rate', 'Average fuel burn rate', 'kg/s', 0.75, 1e-7, 50, 'Mission'), num('reserve_frac', 'Fuel remaining at landing / fuel at start', '-', 0.12, 0, 0.9, 'Mission'),
    num('imbalance', 'Left/right burn mismatch', '-', 0.03, -0.5, 0.5, 'Mission', 'Fraction by which the left wing tank is drawn down faster than the right'),
    num('cg_range_pct', 'Certified CG range available', '% MAC', 20, 1, 60, 'Limits'), num('nSteps', 'Burn steps', '', 120, 10, 5000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c, up, d) => {
    const wing = c.wing.S_m2 > 0 && c.meta.type !== 'helicopter', len = c.fuselage.len_m || 2, big = c.mass.mtow_kg > 30000;
    const geo = { b: c.wing.b_m || 1, c_root: d.c_root || 0.3, taper: c.wing.taper, sweep_deg: c.wing.sweep_deg, eta1: 0.1, eta2: 0.7, x_wing: 0 }, mac = d.mac || 0.3 * len;
    // wing assumed positioned so that full wing tanks sit slightly aft of the zero-fuel CG; centre tank at 40% root chord
    const xw = wing ? 0.05 * mac - wingTank(geo)(1).x : 0;
    return { layout: wing ? LAYOUT[0] : LAYOUT[1], fuel_kg: c.mass.fuel_kg, frac_A: wing ? (big ? 0.3 : 0) : 0.5, x_A: wing ? xw + 0.15 * geo.c_root : -0.06 * len, x_B: 0.08 * len, x_wing: xw, b: c.wing.b_m || 1, c_root: d.c_root || 0.3, taper: c.wing.taper, sweep_deg: c.wing.sweep_deg, zfw_kg: d.zfw_kg || c.mass.mtow_kg - c.mass.fuel_kg, mac: d.mac || 0.3 * len, burn_rate: Math.max(1e-7, fuelFlowEst(c, up, d)), seq: SEQ[0] };
  },
  run(i) {
    const r = cgBurn(i), shift = N.amax(r.xcg) - N.amin(r.xcg), pct = (100 * shift) / i.mac, ymax = N.amax(r.ycg.map(Math.abs)), rmax = N.amax(r.roll.map(Math.abs)), n = r.t.length - 1, warnings = [];
    if (pct > i.cg_range_pct) warnings.push('Fuel burn alone moves the CG by more than the stated certified range: the loading envelope cannot be respected throughout the flight without transfer.');
    if (i.eta2 <= i.eta1 && r.wing) warnings.push('Wing tank outboard end must be outboard of the inboard end.');
    const th = r.t.map((v) => v / 3600);
    return {
      kpis: [
        kp('fuel_cg_shift_m', 'Aircraft CG travel due to fuel burn', shift, 'm'), kp('cg_travel_pct_mac', 'CG travel', pct, '% MAC', pct < 0.5 * i.cg_range_pct ? 'ok' : pct < i.cg_range_pct ? 'warn' : 'bad', `Range available ${i.cg_range_pct}% MAC`),
        kp('cg_start_m', 'CG at start (aft of zero-fuel CG)', r.xcg[0], 'm'), kp('cg_end_m', 'CG at landing (aft of zero-fuel CG)', r.xcg[n], 'm'),
        kp('fuel_cg_start_m', 'Fuel centroid at start', r.xf[0], 'm'), kp('fuel_cg_end_m', 'Fuel centroid at landing', r.xf[n], 'm'),
        kp('lateral_cg_max_m', 'Largest lateral CG offset', ymax, 'm'), kp('roll_moment_max_Nm', 'Largest rolling moment from fuel imbalance', rmax, 'N·m'),
        kp('imbalance_max_kg', 'Largest left/right fuel difference', N.amax(r.mL.map((v, j) => Math.abs(v - r.mR[j]))), 'kg'),
        kp('flight_time_h', 'Time to burn to the landing fuel', r.t[n] / 3600, 'h'), kp('fuel_end_kg', 'Fuel remaining at landing', r.mA[n] + r.mL[n] + r.mR[n], 'kg'),
      ],
      plots: [
        { type: 'line', title: 'Aircraft CG position through the flight', xlabel: 'Time [h]', ylabel: 'CG aft of the zero-fuel CG [% MAC]', series: [{ name: 'Longitudinal CG', x: th, y: r.xcg.map((v) => (100 * v) / i.mac) }] },
        { type: 'line', title: 'Tank quantities', xlabel: 'Time [h]', ylabel: 'Fuel [kg]', series: [{ name: r.wing ? 'Centre tank' : 'Forward tank', x: th, y: r.mA }, { name: r.wing ? 'Left wing tank' : 'Aft tank (half)', x: th, y: r.mL }, { name: r.wing ? 'Right wing tank' : 'Aft tank (half) ', x: th, y: r.mR }] },
        { type: 'line', title: 'Rolling moment from lateral fuel imbalance', xlabel: 'Time [h]', ylabel: 'Rolling moment [kN·m]', series: [{ name: 'Moment (right wing heavy +)', x: th, y: r.roll.map((v) => v / 1e3) }] },
      ],
      warnings, models: ['Quasi-static tank quantity bookkeeping with sequenced burn', 'Wing tank fuel centroid from chord² area distribution with inboard collection', 'Mass-weighted centre of gravity'],
      assumptions: ['Level 1 g attitude: fuel collects at the inboard (lowest) end of each wing tank', 'Tank cross-section proportional to chord squared between the stated span stations, centroid at 40% chord', 'Positions are measured aft of the zero-fuel centre of gravity', 'No fuel transfer other than the burn sequence; unusable fuel ignored', 'Tank span stations, tank positions, burn mismatch and the certified CG range are illustrative defaults, not data for a specific aircraft'],
    };
  },
  verify() {
    const b = { layout: LAYOUT[1], fuel_kg: 1000, frac_A: 0, seq: SEQ[2], x_A: 0, x_B: 2, x_wing: 0, b: 20, c_root: 2, taper: 1, sweep_deg: 0, eta1: 0.1, eta2: 0.7, zfw_kg: 9000, mac: 2, burn_rate: 1, reserve_frac: 0, imbalance: 0, cg_range_pct: 20, nSteps: 50 };
    const r = cgBurn(b), t = wingTank(b)(0.5), w = cgBurn({ ...b, layout: LAYOUT[0], frac_A: 0.3, seq: SEQ[0], reserve_frac: 0.2, imbalance: 0.1 }), k = 30;
    return [
      N.check('Two-mass CG travel m·x/(M+m)', N.amax(r.xcg) - N.amin(r.xcg), (1000 * 2) / 10000, 1e-12, 'Centre-of-mass definition'),
      N.check('Half-full constant-chord tank centroid', t.y, 1 + 0.5 * 6 * 0.5, 5e-4, 'Prismatic tank filled from the inboard end (80-station cumulative table, linear interpolation)'),
      N.check('Fuel mass conservation during the burn', w.mA[k] + w.mL[k] + w.mR[k], 1000 - (800 * k) / 50, 1e-12, 'Σ tank masses = initial − burned'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.cg_travel_pct_mac > 0.5 * i.cg_range_pct) out.push({ severity: o.cg_travel_pct_mac > i.cg_range_pct ? 'critical' : 'warn', title: 'Fuel burn uses a large share of the CG envelope', detail: `CG moves ${o.cg_travel_pct_mac.toFixed(1)}% MAC of ${i.cg_range_pct}% available.`, action: 'Change the burn sequence, reposition tanks nearer the CG, or add automatic fuel transfer; re-check trim drag and static margin in Suite 4 at both ends of the travel.', basis: 'CG envelope' });
    if (o.roll_moment_max_Nm > 0 && o.imbalance_max_kg > 0.03 * i.fuel_kg) out.push({ severity: 'advise', title: 'Lateral fuel imbalance builds up', detail: `Up to ${o.imbalance_max_kg.toFixed(0)} kg difference (${(o.roll_moment_max_Nm / 1e3).toFixed(1)} kN·m).`, action: 'Provide cross-feed and an imbalance alert; the aileron trim needed adds drag and fuel burn.', basis: 'Lateral balance' });
    out.push({ severity: 'info', title: 'CG position affects cruise efficiency', detail: `CG travels from ${((100 * o.cg_start_m) / i.mac).toFixed(1)}% to ${((100 * o.cg_end_m) / i.mac).toFixed(1)}% MAC relative to the zero-fuel CG.`, action: 'Holding the CG towards the aft limit in cruise (by sequencing or trim-tank transfer) lowers tail download and trim drag, typically saving of the order of 1% fuel on long sectors.', basis: 'Trim drag' });
    return out;
  },
};

// ---- 3. sloshing ------------------------------------------------------------------------------
const SHAPE = ['Rectangular tank', 'Upright cylindrical tank'], EXC = ['1 − cos pulse', 'Step', 'Sinusoid'], XI = [1.8412, 5.3314, 8.5363]; // zeros of J1'
/** Linear potential-flow slosh modes: natural frequency and equivalent sloshing mass fraction (exact for these shapes). */
export function sloshModes(shape, L, h, nModes = 3) {
  return N.range(nModes, (n) => {
    if (shape === SHAPE[0]) { const k = ((2 * n + 1) * Math.PI) / L; return { w: Math.sqrt(G0 * k * Math.tanh(k * h)), mf: (8 * Math.tanh(k * h)) / ((2 * n + 1) ** 3 * Math.PI ** 3 * (h / L)) }; }
    const R = L / 2, x = XI[n]; return { w: Math.sqrt(((G0 * x) / R) * Math.tanh((x * h) / R)), mf: (2 * R * Math.tanh((x * h) / R)) / (x * (x * x - 1) * h) };
  });
}
const excite = (i) => (t) => (i.exc === EXC[1] ? i.a_lat : i.exc === EXC[0] ? (t < i.T_exc ? 0.5 * i.a_lat * (1 - Math.cos((2 * Math.PI * t) / i.T_exc)) : 0) : i.a_lat * Math.sin((2 * Math.PI * t) / i.T_exc));
/** 1-D nonlinear shallow-water equations in the tank frame: MUSCL + HLL, Heun time stepping, reflective walls. */
function shallow(i, acc, init, tEnd = i.t_end) {
  const n = Math.max(8, Math.round(i.nCells)), dx = i.L / n, eps = 1e-6 * i.h, mm = (a, b) => (a * b <= 0 ? 0 : Math.abs(a) < Math.abs(b) ? a : b);
  let h = N.range(n, (j) => i.h * (1 + (init ? init((j + 0.5) * dx) : 0))), q = new Array(n).fill(0), t = 0, steps = 0, dry = false;
  const m0 = N.sum(h), flux = (hl, ul, hr, ur) => {
    const cl = Math.sqrt(G0 * hl), cr = Math.sqrt(G0 * hr), sl = Math.min(ul - cl, ur - cr), sr = Math.max(ul + cl, ur + cr), f0 = hl * ul, f1 = hl * ul * ul + 0.5 * G0 * hl * hl, g0 = hr * ur, g1 = hr * ur * ur + 0.5 * G0 * hr * hr;
    if (sl >= 0) return [f0, f1]; if (sr <= 0) return [g0, g1];
    return [(sr * f0 - sl * g0 + sl * sr * (hr - hl)) / (sr - sl), (sr * f1 - sl * g1 + sl * sr * (hr * ur - hl * ul)) / (sr - sl)];
  };
  const rhs = (hh, qq, a) => {
    const u = hh.map((v, j) => (v > eps ? qq[j] / v : 0)), dh = new Array(n), du = new Array(n), F = new Array(n + 1), dH = new Array(n), dQ = new Array(n);
    for (let j = 0; j < n; j++) { dh[j] = mm(hh[j] - (j > 0 ? hh[j - 1] : hh[0]), (j < n - 1 ? hh[j + 1] : hh[n - 1]) - hh[j]); du[j] = mm(u[j] - (j > 0 ? u[j - 1] : -u[0]), (j < n - 1 ? u[j + 1] : -u[n - 1]) - u[j]); }
    for (let f = 0; f <= n; f++) {
      const hr = f < n ? Math.max(eps, hh[f] - dh[f] / 2) : 0, ur = f < n ? u[f] - du[f] / 2 : 0, hl = f > 0 ? Math.max(eps, hh[f - 1] + dh[f - 1] / 2) : 0, ul = f > 0 ? u[f - 1] + du[f - 1] / 2 : 0;
      F[f] = f === 0 ? flux(hr, -ur, hr, ur) : f === n ? flux(hl, ul, hl, -ul) : flux(hl, ul, hr, ur);
    }
    for (let j = 0; j < n; j++) { dH[j] = -(F[j + 1][0] - F[j][0]) / dx; dQ[j] = -(F[j + 1][1] - F[j][1]) / dx - hh[j] * a; }
    return [dH, dQ];
  };
  const T = [0], hL = [h[0]], hR = [h[n - 1]], Fw = [0];
  while (t < tEnd && steps++ < 40000) {
    let smax = 1e-9; for (let j = 0; j < n; j++) smax = Math.max(smax, Math.abs(h[j] > eps ? q[j] / h[j] : 0) + Math.sqrt(G0 * Math.max(h[j], eps)));
    const dt = Math.min((0.4 * dx) / smax, tEnd - t), [a1, b1] = rhs(h, q, acc(t)), h1 = h.map((v, j) => Math.max(eps, v + dt * a1[j])), q1 = q.map((v, j) => v + dt * b1[j]), [a2, b2] = rhs(h1, q1, acc(t + dt));
    h = h.map((v, j) => { const z = 0.5 * (v + h1[j] + dt * a2[j]); if (z < 10 * eps) dry = true; return Math.max(eps, z); }); q = q.map((v, j) => 0.5 * (v + q1[j] + dt * b2[j])); t += dt;
    T.push(t); hL.push(h[0]); hR.push(h[n - 1]); Fw.push(0.5 * G0 * (h[n - 1] ** 2 - h[0] ** 2));
  }
  return { T, hL, hR, Fw, massErr: Math.abs(N.sum(h) - m0) / m0, dry, steps, h, x: N.range(n, (j) => (j + 0.5) * dx) };
}
const slosh = {
  id: 'slosh', title: 'Fuel sloshing: natural frequencies and manoeuvre response', fidelity: 'numerical',
  summary: 'Sloshing frequencies and equivalent sloshing masses of a partly filled tank from potential-flow theory, the lateral force during a manoeuvre, and a nonlinear shallow-water simulation of the free surface that captures steep waves at low fill.',
  equations: ['Sloshing free-surface equations', 'Conservation of mass, momentum and energy equations', 'Hydrostatic pressure equations'],
  applicable: hasFuel,
  inputs: [
    sel('shape', 'Tank shape', SHAPE, SHAPE[0], 'Tank'), num('L', 'Tank length in the direction of motion (or diameter)', 'm', 2, 0.02, 50, 'Tank'), num('Wd', 'Tank width (rectangular)', 'm', 1.5, 0.02, 50, 'Tank'),
    num('h', 'Fuel depth', 'm', 0.3, 0.002, 20, 'Tank'), sel('fuel', 'Fuel', FUELS, FUELS[0], 'Tank'),
    num('zeta', 'Slosh damping ratio', '-', 0.02, 0, 0.5, 'Tank', 'Smooth bare tank ≈ 0.005–0.01; ring or perforated baffles ≈ 0.05–0.15 (typical test values)'),
    sel('exc', 'Lateral acceleration history', EXC, EXC[0], 'Manoeuvre'), num('a_lat', 'Peak lateral acceleration', 'm/s²', 1, 0, 30, 'Manoeuvre'), num('T_exc', 'Pulse duration or oscillation period', 's', 2, 0.02, 120, 'Manoeuvre'),
    num('f_aircraft', 'Aircraft mode frequency to keep clear of', 'Hz', 0.25, 0, 50, 'Manoeuvre', 'Dutch roll or a structural mode; filled from the flight-dynamics suite when available'),
    num('t_end', 'Simulated time', 's', 10, 0.1, 600, 'Numerics'), num('nCells', 'Shallow-water cells', '', 60, 16, 1000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c, up, d) => {
    const f = fuelOf(c.prop.fuel), L = Math.max(0.08, 0.8 * (c.fuselage.dia_m || 0.3)), Wd = Math.max(0.08, 0.4 * (d.c_root || c.fuselage.dia_m || 0.3)), h = N.clamp((0.25 * c.mass.fuel_kg) / f.rho / (L * Wd), 0.04 * L, 0.5 * L), T1 = (2 * Math.PI) / sloshModes(SHAPE[0], L, h, 1)[0].w;
    return { fuel: fuelName(c), L, Wd, h, T_exc: 2 * T1, t_end: 6 * T1, f_aircraft: up.flightdyn?.dr_omega_rads ? up.flightdyn.dr_omega_rads / (2 * Math.PI) : undefined };
  },
  run(i) {
    const f = fuelOf(i.fuel), rect = i.shape === SHAPE[0], vol = rect ? i.L * i.Wd * i.h : (Math.PI * i.L ** 2 * i.h) / 4, m = f.rho * vol, md = sloshModes(i.shape, i.L, i.h, 3), acc = excite(i), warnings = [];
    // modal response: x_n'' + 2ζω x_n' + ω² x_n = −a(t); tank force F = −m0·a − Σ m_n (a + x_n'') = −m0·a + Σ m_n (2ζω x_n' + ω² x_n)
    const ns = 2400, r = N.rk4((t, y) => md.flatMap((q, n) => [y[2 * n + 1], -acc(t) - 2 * i.zeta * q.w * y[2 * n + 1] - q.w ** 2 * y[2 * n]]), 0, new Array(6).fill(0), i.t_end, ns);
    const m0 = m * (1 - N.sum(md.map((q) => q.mf))), F = r.t.map((t, k) => { let s = -m0 * acc(t); md.forEach((q, n) => { s += q.mf * m * (2 * i.zeta * q.w * r.y[k][2 * n + 1] + q.w ** 2 * r.y[k][2 * n]); }); return s; });
    const Fpk = N.amax(F.map(Math.abs)), daf = i.a_lat > 0 ? Fpk / (m * i.a_lat) : 1, f1 = md[0].w / (2 * Math.PI), sw = shallow(i, acc), wall = N.amax(sw.hL.concat(sw.hR)) - i.h, FswPk = f.rho * (rect ? i.Wd : (Math.PI * i.L) / 4) * N.amax(sw.Fw.map(Math.abs));
    const fsw = Math.sqrt(G0 * i.h) / (2 * i.L), hl = i.h / i.L, sep = i.f_aircraft > 0 ? Math.abs(f1 / i.f_aircraft - 1) : 1;
    if (hl > 0.15) warnings.push(`Depth/length = ${hl.toFixed(2)}: the shallow-water simulation overestimates the wave speed (frequency ${fsw.toFixed(2)} Hz against ${f1.toFixed(2)} Hz exact); rely on the modal results at this fill level.`);
    if (!rect) warnings.push('The shallow-water simulation is planar (one horizontal direction); for the cylindrical tank it is indicative only.');
    if (sw.dry) warnings.push('The tank floor is uncovered during the response: the wave becomes a travelling bore and the linear modal model is no longer valid.');
    if (wall > 0.5 * i.h) warnings.push('Wave amplitude exceeds half the fuel depth: strongly nonlinear sloshing (and roof impact, if the tank is shallow) should be expected.');
    if (sep < 0.2) warnings.push(`First slosh frequency ${f1.toFixed(2)} Hz is within 20% of the stated aircraft mode at ${i.f_aircraft.toFixed(2)} Hz: coupled slosh–airframe response is possible.`);
    const hs = N.linspace(0.02 * i.L, i.L, 40);
    return {
      kpis: [
        kp('slosh_freq_Hz', 'First sloshing frequency', f1, 'Hz', sep >= 0.2 ? 'ok' : 'warn'), kp('slosh_freq2_Hz', 'Second sloshing frequency', md[1].w / (2 * Math.PI), 'Hz'),
        kp('slosh_mass_kg', 'First-mode sloshing mass', md[0].mf * m, 'kg'), kp('slosh_mass_frac', 'First-mode sloshing mass / fuel mass', md[0].mf, '-'), kp('fuel_mass_kg', 'Fuel mass in the tank', m, 'kg'),
        kp('slosh_force_peak_N', 'Peak lateral force on the tank (modal)', Fpk, 'N'), kp('dynamic_factor', 'Peak force / (fuel mass × peak acceleration)', daf, '-', daf < 1.5 ? 'ok' : 'warn'),
        kp('sw_wall_rise_m', 'Peak wave height at the wall (shallow-water)', wall, 'm', wall < 0.5 * i.h ? 'ok' : 'warn'), kp('sw_force_peak_N', 'Peak lateral force (shallow-water)', FswPk, 'N'),
        kp('cg_shift_static_m', 'Fuel CG shift under steady peak acceleration', (i.L ** 2 * i.a_lat) / (12 * i.h * G0) * (rect ? 1 : 0.75), 'm', undefined, 'Free-surface effect, I/V·tanθ'),
        kp('sw_freq_Hz', 'Shallow-water first frequency √(gh)/2L', fsw, 'Hz'), kp('sw_mass_error', 'Shallow-water mass conservation error', sw.massErr, '-', sw.massErr < 1e-9 ? 'ok' : 'warn'),
      ],
      plots: [
        { type: 'line', title: 'Lateral force on the tank', xlabel: 'Time [s]', ylabel: 'Force [N]', series: [{ name: 'Linear modal model', x: thin(r.t), y: thin(F) }, { name: 'Nonlinear shallow water', x: thin(sw.T), y: thin(sw.Fw).map((v) => f.rho * (rect ? i.Wd : (Math.PI * i.L) / 4) * v) }, { name: 'Frozen fuel (−m·a)', x: thin(r.t), y: thin(r.t).map((t) => -m * acc(t)), style: 'dash' }] },
        { type: 'line', title: 'Fuel depth at the tank walls (shallow-water)', xlabel: 'Time [s]', ylabel: 'Depth [m]', series: [{ name: 'Wall at x = 0', x: thin(sw.T), y: thin(sw.hL) }, { name: 'Wall at x = L', x: thin(sw.T), y: thin(sw.hR) }], annotations: [{ y: i.h, label: 'Still level' }] },
        { type: 'line', title: 'First sloshing frequency versus fuel depth', xlabel: 'Fuel depth [m]', ylabel: 'Frequency [Hz]', series: [{ name: 'Exact linear theory', x: hs, y: hs.map((h) => sloshModes(i.shape, i.L, h, 1)[0].w / (2 * Math.PI)) }, { name: 'Shallow-water limit', x: hs, y: hs.map((h) => Math.sqrt(G0 * h) / (2 * i.L) * (rect ? 1 : (2 * XI[0]) / Math.PI)), style: 'dash' }], annotations: [{ x: i.h, label: 'Fill' }, ...(i.f_aircraft > 0 ? [{ y: i.f_aircraft, label: 'Aircraft mode' }] : [])] },
      ],
      warnings, models: ['Linear potential-flow slosh modes: ω² = g·k·tanh(k·h)', 'Equivalent mechanical model (sloshing masses on springs, Graham–Rodriguez / Abramson)', 'Nonlinear shallow-water equations: MUSCL–HLL finite volume, Heun time stepping'],
      assumptions: ['Rigid tank, inviscid fuel with equivalent viscous damping ratio', 'Excitation along one horizontal axis; no roof impact', 'Shallow-water model: hydrostatic pressure, depth-uniform velocity (valid for depth/length below about 0.15)', 'No coupling back to the aircraft motion', 'Tank dimensions are scaled from the fuselage and root chord and the damping ratio is a typical test value'],
    };
  },
  convergence: { param: 'nCells', label: 'Shallow-water cells', levels: [15, 30, 60, 120], metric: 'sw_wall_rise_m' },
  verify() {
    const md = sloshModes(SHAPE[0], 2, 0.5, 3), cyl = sloshModes(SHAPE[1], 2, 1, 3), b = { L: 2, h: 0.1, nCells: 100, t_end: 12 };
    const sw = shallow(b, () => 0, (x) => 0.01 * Math.cos((Math.PI * x) / 2), 12), zc = [];
    for (let k = 1; k < sw.T.length; k++) { const a = sw.hL[k - 1] - 0.1, c = sw.hL[k] - 0.1; if (a > 0 && c <= 0) zc.push(sw.T[k - 1] + ((sw.T[k] - sw.T[k - 1]) * a) / (a - c)); }
    const st = shallow({ ...b, nCells: 60 }, () => 0.5, null, 60), k0 = st.T.findIndex((v) => v > 40), tilt = N.mean(st.hL.slice(k0).map((v, k) => v - st.hR[k0 + k])) / (2 - 2 / 60);
    return [
      N.check('First mode ω² = g·k·tanh(k·h), rectangular', md[0].w, Math.sqrt(((G0 * Math.PI) / 2) * Math.tanh((Math.PI / 2) * 0.5)), 1e-12, 'Lamb, Hydrodynamics §257'),
      N.check('Shallow limit: first-mode sloshing mass fraction → 8/π²', sloshModes(SHAPE[0], 1, 1e-4, 1)[0].mf, 8 / Math.PI ** 2, 1e-6, 'Graham & Rodriguez (1952), h/L → 0'),
      N.check('Shallow limit: all modes together carry the whole liquid mass', N.sum(sloshModes(SHAPE[0], 1, 1e-5, 2000).map((q) => q.mf)), 1, 1e-3, 'Σ 8/((2n−1)²π²) = 1'),
      N.check('Cylindrical first-mode mass fraction (R/2.2h)·tanh(1.84h/R)', cyl[0].mf, (1 / 2.2) * Math.tanh(1.8412), 1e-3, 'Abramson, NASA SP-106'),
      N.check('Shallow-water free oscillation period 2L/√(gh)', (zc[zc.length - 1] - zc[0]) / (zc.length - 1), 4 / Math.sqrt(G0 * 0.1), 0.01, 'Linear long-wave theory'),
      N.check('Steady free-surface slope −a/g under constant acceleration', tilt, 0.5 / G0, 0.03, 'Hydrostatic equilibrium, time-averaged over the residual oscillation'),
      N.check('Shallow-water mass conservation', sw.massErr, 0, 1e-11, 'Finite-volume conservation'),
    ];
  },
  calibration: { params: [{ key: 'zeta', min: 0, max: 0.4 }], sweep: 'T_exc', target: 'slosh_force_peak_N', note: 'Supply measured peak slosh force against excitation period from a shaker-table tank test to fit the damping ratio.' },
  recommend(res, i, ctx) {
    const o = res.outputs, out = [], share = o.slosh_mass_kg / (ctx?.case?.mass?.mtow_kg || Infinity); // coupling strength grows with the sloshing mass relative to the aircraft
    if (i.f_aircraft > 0 && Math.abs(o.slosh_freq_Hz / i.f_aircraft - 1) < 0.2) out.push({ severity: share > 0.05 ? 'warn' : 'advise', title: 'Slosh frequency close to an aircraft mode', detail: `${o.slosh_freq_Hz.toFixed(2)} Hz against ${i.f_aircraft.toFixed(2)} Hz${share > 0 ? `; the sloshing mass is ${(100 * share).toFixed(1)}% of the take-off mass` : ''}. The tank size is a generic default: enter the real free-surface length between ribs or baffles.`, action: 'Add baffles or ribs to shorten the free-surface length (frequency rises as the cell length falls), or check the coupled response in Suite 4.', basis: 'Frequency separation of 20% (rule of thumb); the coupling matters once the sloshing mass exceeds a few percent of the aircraft mass' });
    if (o.dynamic_factor > 1.3) out.push({ severity: 'advise', title: 'Sloshing amplifies the manoeuvre load', detail: `Peak tank force is ${o.dynamic_factor.toFixed(2)} times the frozen-fuel value.`, action: 'Increase damping with perforated baffles, and use the peak force for tank-wall and attachment loads in Suite 2.', basis: 'Equivalent mechanical slosh model' });
    if (o.sw_wall_rise_m > 0.5 * i.h) out.push({ severity: 'advise', title: 'Steep waves at this fill level', detail: `Wall wave height ${o.sw_wall_rise_m.toFixed(2)} m on ${i.h.toFixed(2)} m depth.`, action: 'Check pump-inlet uncovering (use a collector cell with flapper valves) and consider a volume-of-fluid CFD study for impact pressures.', basis: 'Nonlinear shallow-water response' });
    return out;
  },
};

// ---- 4. fuel thermal state and tank inerting ---------------------------------------------------
function fuelThermal(i) {
  const f = fuelOf(i.fuel), nt = Math.max(4, Math.round(i.nSteps)), mMin = 0.02 * i.m_fuel, mOf = (t) => Math.max(mMin, i.m_fuel - i.burn_rate * t), rhoG = i.p_tank / (R_AIR * i.T_ull), Qn = i.mdot_nea / rhoG;
  const Vu = (t) => Math.max(0.01 * i.V_tank, i.V_tank - mOf(t) / f.rho), dVu = (t) => (mOf(t) > mMin ? i.burn_rate / f.rho : 0);
  const r = N.rk4((t, y) => { const qa = Math.max(0, dVu(t) - Qn); return [(i.U * i.A_wet * (i.T_skin - y[0]) + i.Q_return) / (mOf(t) * f.cp), (Qn * (i.x_nea - y[1]) + qa * (0.21 - y[1])) / Vu(t)]; }, 0, [i.T_fuel0, i.x_o2_0], i.t_end, nt);
  return { t: r.t, T: r.y.map((y) => y[0]), x: r.y.map((y) => y[1]), f, m: r.t.map(mOf), Qn, Teq: i.T_skin + i.Q_return / Math.max(i.U * i.A_wet, 1e-12), Vu0: Vu(0) };
}
const cross = (t, y, lim) => { for (let k = 1; k < t.length; k++) if (y[k] <= lim && y[k - 1] > lim) return t[k - 1] + ((t[k] - t[k - 1]) * (y[k - 1] - lim)) / (y[k - 1] - y[k]); return y[0] <= lim ? 0 : NaN; };
const fueltherm = {
  id: 'fueltherm', title: 'Fuel temperature in cruise and tank inerting', fidelity: 'reduced-order',
  summary: 'Cold-soak of tank fuel towards the skin recovery temperature with heat returned from the engine oil coolers, the margin to fuel freezing, and the wash-out of oxygen from the ullage by nitrogen-enriched air.',
  equations: ['Conservation of mass, momentum and energy equations', 'Heat exchanger equations', 'Ideal gas law', 'Ventilation transport equations'],
  applicable: hasFuel,
  inputs: [
    sel('fuel', 'Fuel', FUELS, FUELS[0], 'Fuel'), num('m_fuel', 'Fuel mass in the tank at start', 'kg', 18800, 0.01, 5e5, 'Fuel'), num('burn_rate', 'Fuel burn rate', 'kg/s', 0.75, 0, 50, 'Fuel'),
    num('T_fuel0', 'Fuel temperature at start', 'K', 293, 200, 340, 'Fuel', 'Uplift temperature, filled from the site temperature'),
    num('T_skin', 'Tank skin (recovery) temperature in cruise', 'K', 243, 180, 340, 'Heat transfer', 'From the thermal suite recovery temperature when available'),
    num('U', 'Fuel-to-skin heat-transfer coefficient', 'W/m²/K', 15, 0.5, 200, 'Heat transfer', 'Natural convection in the fuel in series with the external boundary layer; typical 5–30'),
    num('A_wet', 'Tank skin area in contact with fuel', 'm²', 110, 0.01, 5000, 'Heat transfer'), num('Q_return', 'Heat returned to the tank (oil coolers, recirculation)', 'W', 0, 0, 5e6, 'Heat transfer'),
    num('T_margin', 'Required margin above the freezing point', 'K', 3, 0, 20, 'Limits', 'Operators usually hold fuel about 3 K above the specification freezing point'),
    num('V_tank', 'Tank total volume', 'm³', 25, 1e-4, 1000, 'Inerting'), num('p_tank', 'Ullage pressure', 'Pa', 23800, 3000, 2e5, 'Inerting'), num('T_ull', 'Ullage gas temperature', 'K', 260, 180, 340, 'Inerting'),
    num('mdot_nea', 'Nitrogen-enriched air flow', 'kg/s', 0.012, 0, 5, 'Inerting'), num('x_nea', 'Oxygen fraction of the enriched air', '-', 0.05, 0, 0.21, 'Inerting', 'Air-separation modules deliver about 2–12% oxygen depending on flow'),
    num('x_o2_0', 'Initial ullage oxygen fraction', '-', 0.21, 0, 0.21, 'Inerting'), num('x_limit', 'Inert oxygen limit', '-', 0.12, 0.05, 0.2, 'Inerting', 'About 12% by volume at sea level, rising with altitude (commonly used criterion)'),
    num('t_end', 'Cruise time', 's', 18000, 10, 2e5, 'Numerics'), num('nSteps', 'Time steps', '', 400, 20, 100000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c, up, d) => {
    const { a, V, M } = cruise(c), f = fuelOf(c.prop.fuel), rate = Math.max(1e-7, fuelFlowEst(c, up, d)), Vt = (1.08 * c.mass.fuel_kg) / f.rho, Trec = up.thermal?.T_recovery_K ?? a.T * (1 + 0.178 * M * M);
    return { fuel: fuelName(c), m_fuel: c.mass.fuel_kg, burn_rate: rate, T_fuel0: siteT(c), T_skin: Trec, A_wet: Math.max(0.05, c.wing.S_m2 > 0 ? 0.9 * c.wing.S_m2 : 6 * Vt ** (2 / 3) * 0.5), V_tank: Vt, p_tank: a.p, T_ull: 0.5 * (Trec + siteT(c)), mdot_nea: ((Vt / 3000) * a.p) / (R_AIR * 0.5 * (Trec + siteT(c))), t_end: Math.max(600, Math.min((0.85 * c.mass.fuel_kg) / rate, ((c.mission.range_km || 300) * 1e3) / Math.max(V, 1))) };
  },
  run(i) {
    const r = fuelThermal(i), n = r.t.length - 1, Tmin = N.amin(r.T), lim = r.f.Tfr + i.T_margin, tFr = cross(r.t, r.T, lim), tIn = cross(r.t, r.x, i.x_limit), warnings = [];
    const qNeed = i.x_limit > i.x_nea ? ((r.Vu0 + i.V_tank) / 2 / 1200) * Math.log((0.21 - i.x_nea) / (i.x_limit - i.x_nea)) * (i.p_tank / (R_AIR * i.T_ull)) : NaN;
    if (Tmin < lim) warnings.push(`Fuel reaches ${Tmin.toFixed(1)} K, inside the ${i.T_margin} K margin above the ${r.f.Tfr.toFixed(1)} K freezing point: wax crystals can block filters and reduce usable fuel.`);
    if (Number.isNaN(tIn) && r.x[n] > i.x_limit) warnings.push('The ullage does not become inert within the cruise time at this enriched-air flow.');
    if (N.amax(r.T) > 328) warnings.push('Fuel temperature exceeds about 55 °C: vapour pressure and flammability exposure rise, and tank sealant limits should be checked.');
    const th = r.t.map((v) => v / 3600);
    return {
      kpis: [
        kp('T_fuel_end_K', 'Fuel temperature at end of cruise', r.T[n], 'K'), kp('T_fuel_min_K', 'Lowest fuel temperature', Tmin, 'K', Tmin >= lim + 5 ? 'ok' : Tmin >= lim ? 'warn' : 'bad'),
        kp('freeze_margin_K', 'Margin above the freezing point', Tmin - r.f.Tfr, 'K', Tmin >= lim + 5 ? 'ok' : Tmin >= lim ? 'warn' : 'bad', `Freezing point ${r.f.Tfr.toFixed(1)} K (specification maximum)`),
        kp('t_to_freeze_limit_h', 'Time to reach the freezing margin', Number.isNaN(tFr) ? i.t_end / 3600 : tFr / 3600, 'h', Number.isNaN(tFr) ? 'ok' : 'bad', Number.isNaN(tFr) ? 'Not reached within the cruise time' : ''),
        kp('T_fuel_eq_K', 'Equilibrium fuel temperature', r.Teq, 'K'), kp('tau_h', 'Cold-soak time constant at start', (i.m_fuel * r.f.cp) / (i.U * i.A_wet) / 3600, 'h'),
        kp('heat_loss_W', 'Heat lost through the skin at start', i.U * i.A_wet * (i.T_fuel0 - i.T_skin), 'W'),
        kp('o2_end_pct', 'Ullage oxygen at end', 100 * r.x[n], '%', r.x[n] <= i.x_limit ? 'ok' : 'warn', `Limit ${(100 * i.x_limit).toFixed(1)}%`),
        kp('t_inert_min', 'Time to reach the inert limit', Number.isNaN(tIn) ? i.t_end / 60 : tIn / 60, 'min', Number.isNaN(tIn) ? 'warn' : 'ok', Number.isNaN(tIn) ? 'Not reached within the cruise time' : ''),
        kp('nea_flow_20min_kgs', 'Enriched-air flow to inert the mean ullage in 20 min', qNeed, 'kg/s'), kp('nea_volume_flow_m3s', 'Enriched-air volume flow at tank conditions', r.Qn, 'm³/s'),
      ],
      plots: [
        { type: 'line', title: 'Bulk fuel temperature', xlabel: 'Time [h]', ylabel: 'Temperature [K]', series: [{ name: 'Fuel', x: thin(th), y: thin(r.T) }], annotations: [{ y: lim, label: 'Freezing point + margin' }, { y: i.T_skin, label: 'Skin' }] },
        { type: 'line', title: 'Ullage oxygen concentration', xlabel: 'Time [h]', ylabel: 'Oxygen [% by volume]', series: [{ name: 'Ullage oxygen', x: thin(th), y: thin(r.x).map((v) => 100 * v) }], annotations: [{ y: 100 * i.x_limit, label: 'Inert limit' }] },
        { type: 'line', title: 'Fuel remaining', xlabel: 'Time [h]', ylabel: 'Fuel mass [kg]', series: [{ name: 'Fuel in tank', x: thin(th), y: thin(r.m) }] },
      ],
      warnings, models: ['Lumped (well-mixed) fuel energy balance with decreasing mass', 'Well-mixed ullage oxygen balance with enriched-air wash and vent make-up', 'Specification-maximum freezing points'],
      assumptions: ['Uniform bulk fuel temperature (no stratification; fuel next to the skin is colder than the bulk)', 'Constant heat-transfer coefficient and wetted area', 'Ullage at constant pressure and temperature; fuel vapour and oxygen evolution from the fuel neglected', 'Vent make-up is ambient air at 21% oxygen when enriched-air flow is below the ullage growth rate', 'Heat-transfer coefficient, enriched-air flow and purity, the 12% oxygen limit and the 3 K freezing margin are typical or commonly used values; volatility constants are representative'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [25, 50, 100, 200, 400], metric: 'T_fuel_end_K' },
  verify() {
    const b = { fuel: FUELS[0], m_fuel: 1000, burn_rate: 0, T_fuel0: 300, T_skin: 240, U: 20, A_wet: 10, Q_return: 0, T_margin: 3, V_tank: 2.5, p_tank: 50000, T_ull: 250, mdot_nea: 0.002, x_nea: 0.05, x_o2_0: 0.21, x_limit: 0.12, t_end: 4000, nSteps: 200 };
    const r = fuelThermal(b), tau = (1000 * 2010) / 200, Vu = 2.5 - 1000 / 804, Q = 0.002 / (50000 / (R_AIR * 250));
    return [
      N.check('Lumped cold-soak T = Ts + (T0 − Ts)·exp(−t/τ)', r.T[200], 240 + 60 * Math.exp(-4000 / tau), 1e-9, 'Lumped-capacitance solution'),
      N.check('Ullage wash-out x = xn + (x0 − xn)·exp(−Q·t/V)', r.x[200], 0.05 + 0.16 * Math.exp((-Q * 4000) / Vu), 1e-7, 'Well-mixed dilution'),
      N.check('Equilibrium with returned heat: Ts + Q/(UA)', fuelThermal({ ...b, Q_return: 4000 }).Teq, 260, 1e-12, 'Steady energy balance'),
    ];
  },
  calibration: { params: [{ key: 'U', min: 1, max: 100 }, { key: 'A_wet', min: 0.01, max: 2000 }], sweep: 't_end', target: 'T_fuel_end_K', note: 'Supply measured bulk fuel temperature against cruise time from flight test or a cold-soak chamber test.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.freeze_margin_K < i.T_margin) out.push({ severity: 'critical', title: 'Fuel freezing margin is lost', detail: `Margin ${o.freeze_margin_K.toFixed(1)} K after ${o.t_to_freeze_limit_h.toFixed(1)} h.`, action: 'Fly faster or lower (warmer recovery temperature), return more oil-cooler heat to the tank, or load fuel with a lower freezing point.', basis: 'Fuel specification freezing point plus operating margin' });
    else if (o.freeze_margin_K < i.T_margin + 5) out.push({ severity: 'warn', title: 'Fuel is close to its freezing margin', detail: `Lowest bulk temperature ${o.T_fuel_min_K.toFixed(1)} K.`, action: 'Monitor fuel temperature on long cold sectors; plan a descent or Mach increase option.', basis: 'Fuel temperature monitoring practice' });
    if (o.o2_end_pct > 100 * i.x_limit) out.push({ severity: 'advise', title: 'Ullage is not inert', detail: `${o.o2_end_pct.toFixed(1)}% oxygen at the end of cruise.`, action: `Increase the enriched-air flow (about ${Number.isFinite(o.nea_flow_20min_kgs) ? o.nea_flow_20min_kgs.toExponential(2) : '—'} kg/s inerts the mean ullage in 20 minutes) or its purity; note that descent draws in ambient air.`, basis: 'Ullage oxygen balance' });
    if (i.Q_return > 0) out.push({ severity: 'info', title: 'Waste heat is used to protect the fuel', detail: `${(i.Q_return / 1e3).toFixed(1)} kW of oil-cooler heat raises the equilibrium fuel temperature to ${o.T_fuel_eq_K.toFixed(0)} K.`, action: 'Recovering oil heat into the fuel avoids ram-air cooler drag and improves freeze margin at no fuel cost.', basis: 'Energy recovery' });
    return out;
  },
};

// ---- 5. cabin / bay air: heat load, ventilation, CO2 and humidity -----------------------------
function cabinCalc(i) {
  const Qocc = i.n_occ * i.q_occ, Qsol = i.A_glaze * i.tau_glass * i.G_solar, Qwall = i.U_wall * i.A_wall * (i.T_skin - i.T_cab), Q = Qocc + i.Q_equip + Qsol + Qwall;
  const mFresh = Math.max(i.n_occ * i.m_fresh_pp, i.m_vent), mLoad = Q > 0 ? Q / (CP_AIR * Math.max(1, i.T_cab - i.T_supply)) : 0, mSup = Math.max(mLoad, mFresh / Math.max(1e-6, 1 - i.recirc)), rho = i.p_cab / (R_AIR * i.T_cab);
  // CO2 by mole fraction; water by humidity ratio
  const nF = mFresh / 0.028965, nC = (i.n_occ * i.co2_pp) / 0.04401, Ncab = (i.p_cab * i.V_cab) / (8.31446 * i.T_cab), Css = i.co2_out * 1e-6 + (nF > 0 ? nC / nF : 0), tauC = nF > 0 ? Ncab / nF : Infinity;
  const wOut = humRatio(i.rh_out * pSat(i.T_out), i.p_out), wCab = wOut + (mFresh > 0 ? (i.n_occ * i.h2o_pp) / mFresh : 0), pv = (wCab * i.p_cab) / (0.622 + wCab), rh = pv / pSat(i.T_cab), Td = dewPoint(pv);
  const Twall = i.T_cab - (i.U_wall * (i.T_cab - i.T_skin)) / i.h_in;
  return { Qocc, Qsol, Qwall, Q, mFresh, mLoad, mSup, rho, Css, tauC, nF, nC, Ncab, wOut, wCab, rh, Td, Twall, ach: (mSup / rho / i.V_cab) * 3600 };
}
const cabin = {
  id: 'cabin', title: 'Cabin and equipment-bay air: heat load, ventilation, CO₂ and humidity', fidelity: 'reduced-order',
  summary: 'Adds up occupant, equipment, solar and wall heat loads, sizes the supply and fresh-air flows, and predicts carbon-dioxide level, humidity, dew point and wall condensation for a cabin or an equipment bay.',
  equations: ['Ventilation transport equations', 'Psychrometric relations', 'Ideal gas law', 'Conservation of mass, momentum and energy equations', 'Heat exchanger equations'],
  inputs: [
    num('n_occ', 'Occupants', '', 171, 0, 1000, 'Heat loads', 'Passengers plus crew; 0 for an equipment or payload bay', { step: 1, discrete: true }), num('q_occ', 'Sensible heat per occupant', 'W', 100, 30, 300, 'Heat loads', 'Seated adult ≈ 70–100 W sensible'),
    num('Q_equip', 'Equipment and lighting heat', 'W', 12000, 0, 5e6, 'Heat loads'), num('A_glaze', 'Window area exposed to the sun', 'm²', 6, 0, 500, 'Heat loads'), num('tau_glass', 'Window solar transmittance', '-', 0.6, 0, 1, 'Heat loads'),
    num('G_solar', 'Solar irradiance', 'W/m²', 900, 0, 1400, 'Heat loads'), num('U_wall', 'Wall overall heat-transfer coefficient', 'W/m²/K', 1.0, 0.05, 50, 'Heat loads', 'Insulated cabin wall ≈ 0.7–1.5; uninsulated skin ≈ 4–8'),
    num('A_wall', 'Wall area', 'm²', 420, 0.01, 10000, 'Heat loads'), num('T_skin', 'Outside skin temperature', 'K', 318, 180, 400, 'Heat loads', 'Ground hot soak, or recovery temperature in flight'),
    num('T_cab', 'Target cabin or bay temperature', 'K', 297, 250, 350, 'Air supply'), num('T_supply', 'Supply air temperature', 'K', 283, 220, 350, 'Air supply', 'Above about 275 K to avoid icing in the ducts'),
    num('m_fresh_pp', 'Fresh air per occupant', 'kg/s', 0.0042, 0, 0.05, 'Air supply', '0.0042 kg/s (0.55 lb/min) is the usual transport-category design minimum'), num('m_vent', 'Minimum ventilation flow (unoccupied bay)', 'kg/s', 0, 0, 100, 'Air supply'),
    num('recirc', 'Recirculated fraction of the supply', '-', 0.4, 0, 0.9, 'Air supply'), num('V_cab', 'Cabin or bay volume', 'm³', 240, 0.001, 5000, 'Air supply'), num('p_cab', 'Cabin or bay pressure', 'Pa', 101325, 5000, 110000, 'Air supply'),
    num('co2_pp', 'CO₂ generated per occupant', 'kg/s', 1.0e-5, 0, 1e-4, 'Air quality', 'About 0.3 L/min at rest (typical)'), num('co2_out', 'Outside-air CO₂', 'ppm', 420, 250, 2000, 'Air quality'),
    num('h2o_pp', 'Moisture released per occupant', 'kg/s', 1.4e-5, 0, 1e-4, 'Air quality', 'About 50 g/h at rest (typical)'),
    num('T_out', 'Outside air temperature', 'K', 303, 180, 330, 'Outside air', 'Filled from the site weather'), num('p_out', 'Outside air pressure', 'Pa', 101325, 5000, 110000, 'Outside air'), num('rh_out', 'Outside relative humidity', '-', 0.5, 0, 1, 'Outside air'),
    num('h_in', 'Inside wall film coefficient', 'W/m²/K', 8, 1, 50, 'Outside air'), num('t_end', 'Time simulated for the CO₂ build-up', 's', 3600, 10, 1e5, 'Numerics'), num('nSteps', 'Time steps', '', 200, 10, 20000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c, up, d) => {
    const pax = c.mission.pax || 0, occ = pax > 0 ? pax + (pax > 19 ? Math.ceil(pax / 50) + 2 : 1) : 0, len = c.fuselage.len_m || 1, dia = c.fuselage.dia_m || 0.3, Ts = siteT(c), pS = (c.site?.p_hPa ?? 1013.25) * 100;
    const Vc = occ > 0 ? 0.55 * (Math.PI / 4) * dia * dia * len * 0.7 : 0.15 * (Math.PI / 4) * dia * dia * len, Aw = (occ > 0 ? 0.7 : 0.3) * Math.PI * dia * len, { a, M } = cruise(c);
    if (occ > 0) return { n_occ: occ, q_occ: c.systems.pax_heat_W, Q_equip: 40 * occ + 0.02 * (c.systems.gen_kVA || 0) * 1000 * c.prop.n_eng + 300, A_glaze: 0.035 * occ + 1, A_wall: Aw, V_cab: Vc, U_wall: c.mass.mtow_kg > 5700 ? 1 : 3, T_skin: up.thermal?.T_soak_K ?? Ts + 15, T_out: Ts, p_out: pS, rh_out: c.site?.rh ?? 0.5, p_cab: pS, recirc: c.mass.mtow_kg > 5700 ? 0.4 : 0 };
    const Q = Math.max(5, up.electrical?.elec_losses_W ?? (c.prop.type === 'electric' ? 0.06 * d.P_total : 60 + 0.03 * (c.systems.gen_kVA || 0) * 1000)), Trec = a.T * (1 + 0.178 * M * M);
    return { n_occ: 0, Q_equip: Q, A_glaze: 0, A_wall: Aw, V_cab: Math.max(1e-3, Vc), U_wall: 5, T_skin: Trec, T_out: a.T, p_out: a.p, rh_out: c.site?.rh ?? 0.5, p_cab: a.p, T_cab: Trec + 20, T_supply: Trec, m_vent: Q / (CP_AIR * 20) * 0.2, recirc: 0 };
  },
  run(i, ctx) {
    const r = cabinCalc(i), warnings = [], nt = Math.max(10, Math.round(i.nSteps)), C0 = i.co2_out * 1e-6;
    const co = r.nF > 0 || r.nC > 0 ? N.rk4((t, y) => [(r.nF * (C0 - y[0]) + r.nC) / r.Ncab], 0, [C0], i.t_end, nt) : { t: [0, i.t_end], y: [[C0], [C0]] }, ppm = co.y.map((y) => 1e6 * y[0]), cond = r.Twall < r.Td;
    if (r.Css > 5000e-6) warnings.push('Steady CO₂ exceeds 5000 ppm (0.5% by volume), the usual airworthiness limit for occupied compartments.');
    if (cond) warnings.push(`Inner wall surface (${r.Twall.toFixed(1)} K) is below the cabin dew point (${r.Td.toFixed(1)} K): condensation or frost will form on the structure.`);
    if (r.rh < 0.1 && i.n_occ > 0) warnings.push('Cabin relative humidity is below 10%: typical of high-altitude cruise with dry outside air; a comfort rather than a safety issue.');
    if (i.T_supply >= i.T_cab && r.Q > 0) warnings.push('Supply air is not colder than the target temperature: the cooling load cannot be removed.');
    if (r.mLoad > 3 * Math.max(r.mFresh, 1e-9) / Math.max(1e-6, 1 - i.recirc) && i.n_occ > 0) warnings.push('The flow needed for cooling is far above the ventilation minimum: cooling, not air quality, sizes the air supply at this condition.');
    const Tsk = N.linspace(220, 340, 40), noFuel = ctx?.case && !(ctx.case.mass.fuel_kg > 0), unp = ctx?.case && !(ctx.case.fuselage.cabin_dp_Pa > 0), outputs = {};
    if (noFuel) Object.assign(outputs, { fuel_cg_shift_m: 0, fuel_pump_power_W: 0 });
    if (unp) outputs.cabin_alt_m = ctx.case.mission.cruise_alt_m || ctx.case.atm.alt_m || 0;
    if (ctx?.case) outputs.bleed_kgs = ctx.case.systems.bleed ? r.mFresh : 0;
    return {
      kpis: [
        kp('ecs_cooling_W', 'Cooling load', Math.max(0, r.Q), 'W'), kp('ecs_heating_W', 'Heating load', Math.max(0, -r.Q), 'W'), kp('Q_net_W', 'Net heat gain', r.Q, 'W'),
        kp('supply_flow_kgs', 'Supply air flow required', r.mSup, 'kg/s'), kp('fresh_flow_kgs', 'Fresh air flow', r.mFresh, 'kg/s'), kp('air_changes_per_h', 'Air changes per hour', r.ach, '1/h'),
        kp('co2_ss_ppm', 'Steady CO₂ concentration', 1e6 * r.Css, 'ppm', r.Css < 2500e-6 ? 'ok' : r.Css < 5000e-6 ? 'warn' : 'bad', 'Airliner cabins typically run at about 1000–2000 ppm; limit 5000 ppm'),
        kp('co2_end_ppm', 'CO₂ at end of the simulated time', ppm[ppm.length - 1], 'ppm'), kp('co2_tau_min', 'CO₂ build-up time constant', Number.isFinite(r.tauC) ? r.tauC / 60 : 0, 'min'),
        kp('rh_cabin', 'Relative humidity', r.rh, '-', r.rh <= 0.7 ? 'ok' : 'warn'), kp('dew_point_K', 'Dew point', r.Td, 'K'), kp('T_wall_inner_K', 'Inner wall surface temperature', r.Twall, 'K', cond ? 'warn' : 'ok'),
        kp('cabin_pressure_alt_m', 'Pressure altitude of the compartment', pressureAltitude(i.p_cab), 'm'),
      ],
      plots: [
        { type: 'bar', title: 'Heat load budget', ylabel: 'Heat gain [W]', categories: ['Occupants', 'Equipment', 'Solar', 'Walls', 'Net'], series: [{ name: 'Heat gain', y: [r.Qocc, i.Q_equip, r.Qsol, r.Qwall, r.Q] }] },
        { type: 'line', title: 'CO₂ build-up from outside-air level', xlabel: 'Time [min]', ylabel: 'CO₂ [ppm]', series: [{ name: 'Compartment CO₂', x: thin(co.t).map((v) => v / 60), y: thin(ppm) }], annotations: [{ y: 5000, label: 'Limit' }] },
        { type: 'line', title: 'Net heat gain versus outside skin temperature', xlabel: 'Skin temperature [K]', ylabel: 'Net heat gain [kW]', series: [{ name: 'Net load (cooling +, heating −)', x: Tsk, y: Tsk.map((T) => cabinCalc({ ...i, T_skin: T }).Q / 1e3) }], annotations: [{ x: i.T_skin, label: 'Condition' }, { y: 0, label: 'Balance' }] },
      ],
      outputs, warnings, models: ['Steady sensible heat balance', 'Well-mixed CO₂ dilution (mole balance)', 'Psychrometrics: humidity ratio 0.622·p_v/(p − p_v), Magnus saturation pressure', 'Wall surface temperature from series film resistance'],
      assumptions: ['Perfectly mixed compartment air', 'Sensible loads only for the supply sizing; latent load reported through humidity', 'Recirculated air returns at cabin temperature; fan heat neglected', 'Uniform wall coefficient and skin temperature', 'Occupant heat, CO₂ and moisture rates, equipment heat, glazing and wall coefficients are typical values'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [10, 20, 40, 80], metric: 'co2_end_ppm' },
  verify() {
    const b = { n_occ: 100, q_occ: 100, Q_equip: 5000, A_glaze: 10, tau_glass: 0.5, G_solar: 1000, U_wall: 1, A_wall: 300, T_skin: 320, T_cab: 297, T_supply: 283, m_fresh_pp: 0.005, m_vent: 0, recirc: 0.5, V_cab: 200, p_cab: 80000, co2_pp: 1e-5, co2_out: 400, h2o_pp: 0, T_out: 293.15, p_out: 101325, rh_out: 0.5, h_in: 8, t_end: 600, nSteps: 100 };
    const r = cabinCalc(b), o = N.kv(cabin.run(b)), tau = r.Ncab / r.nF, pv = 0.5 * pSat(293.15);
    return [
      N.check('Heat load sum', r.Q, 10000 + 5000 + 5000 + 300 * 23, 1e-12, 'Energy balance'),
      N.check('Steady CO₂ = outside + generation/fresh-air (mole basis)', r.Css, 400e-6 + (100 * 1e-5) / 0.04401 / (0.5 / 0.028965), 1e-12, 'Dilution equation'),
      N.check('CO₂ transient C = Css − (Css − C0)·exp(−t/τ)', o.co2_end_ppm, 1e6 * (r.Css - (r.Css - 400e-6) * Math.exp(-600 / tau)), 1e-8, 'Well-mixed first-order response'),
      N.check('Dew point inverts the saturation curve', pSat(dewPoint(pv)), pv, 1e-10, 'Magnus formula'),
      N.check('Supply flow removes the load: m·cp·ΔT = Q', o.supply_flow_kgs * CP_AIR * 14, r.Q, 1e-12, 'Sensible heat equation'),
    ];
  },
  calibration: { params: [{ key: 'U_wall', min: 0.2, max: 20 }, { key: 'tau_glass', min: 0.1, max: 0.9 }], sweep: 'T_skin', target: 'Q_net_W', note: 'Supply measured cabin heat load (pack duty) against skin temperature from environmental-chamber or ground hot-soak tests.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.co2_ss_ppm > 2500 && i.n_occ > 0) out.push({ severity: o.co2_ss_ppm > 5000 ? 'critical' : 'warn', title: 'Fresh-air rate is low for air quality', detail: `Steady CO₂ ${o.co2_ss_ppm.toFixed(0)} ppm.`, action: 'Increase fresh air per occupant or reduce the recirculated share.', basis: 'CO₂ dilution; 5000 ppm airworthiness limit, typical cabins 1000–2000 ppm' });
    if (o.ecs_cooling_W > 0 && res.kpis.find((k) => k.key === 'Q_net_W') && i.A_glaze * i.tau_glass * i.G_solar > 0.25 * o.ecs_cooling_W) out.push({ severity: 'advise', title: 'Solar gain is a large share of the cooling load', detail: `${((100 * i.A_glaze * i.tau_glass * i.G_solar) / o.ecs_cooling_W).toFixed(0)}% of the load enters through the glazing.`, action: 'Window shades or low-transmittance glazing on the ground reduce pack duty, ground power or APU fuel and pull-down time.', basis: 'Heat load budget' });
    if (i.n_occ > 0 && i.recirc < 0.3 && o.fresh_flow_kgs > 0.5) out.push({ severity: 'advise', title: 'Recirculation would cut outside-air demand', detail: `Only ${(100 * i.recirc).toFixed(0)}% of the supply is recirculated.`, action: 'Filtered recirculation of 40–50% reduces the fresh (bleed or compressed) air flow and its fuel penalty while keeping the fresh-air minimum per occupant.', basis: 'Fresh-air energy cost' });
    if (i.n_occ === 0) out.push({ severity: 'info', title: 'Equipment bay ventilation', detail: `${(o.supply_flow_kgs * 1e3).toFixed(1)} g/s of air removes ${o.ecs_cooling_W.toFixed(0)} W with the stated temperature rise.`, action: 'Size the ram-air inlet for this flow at the slowest flight speed and check hover or ground operation, when there is no ram air, in Suite 12.', basis: 'Sensible heat balance' });
    return out;
  },
};

// ---- 6. pressurisation and decompression -------------------------------------------------------
const PROC = ['Isothermal', 'Adiabatic (isentropic)'];
/** Cabin pressure for an aircraft altitude: proportional cabin-altitude schedule capped by the maximum differential. */
function schedule(i, h) {
  const frac = N.clamp((h - i.field_m) / Math.max(1, i.ceiling_m - i.field_m), 0, 1), hc = Math.min(h, i.field_m + frac * (i.cab_alt_max - i.field_m)), pa = isa(h).p, ps = isa(hc).p, limited = ps - pa > i.dp_max, pc = limited ? pa + i.dp_max : ps;
  return { pa, pc, hc: limited ? pressureAltitude(pc) : hc, dp: pc - pa, limited };
}
function decomp(i, s) {
  const A = (i.Cd * Math.PI * i.d_hole ** 2) / 4 + i.A_leak, m0 = (s.pc * i.V_cab) / (R_AIR * i.T_cab), iso = i.process === PROC[0], nt = Math.max(10, Math.round(i.nSteps));
  const Tof = (m) => (iso ? i.T_cab : i.T_cab * (m / m0) ** (g - 1)), pof = (m) => (m * R_AIR * Tof(m)) / i.V_cab;
  const r = N.rk4((t, y) => [i.m_in - orifice(A, pof(y[0]), Tof(y[0]), s.pa)], 0, [m0], i.t_end, nt), p = r.y.map((y) => pof(y[0]));
  return { t: r.t, p, alt: p.map(pressureAltitude), T: r.y.map((y) => Tof(y[0])), A, m0 };
}
const pressurisation = {
  id: 'pressurisation', title: 'Cabin pressurisation schedule and decompression', fidelity: 'numerical',
  summary: 'Cabin altitude and pressure differential against flight altitude, the outflow-valve flow that holds it, cabin rate of change in a climb, and how quickly the cabin loses pressure through a failed window or door seal.',
  equations: ['Cabin pressure balance equations', 'Ideal gas law', 'Incompressible and compressible fluid flow equations', 'Conservation of mass, momentum and energy equations'],
  applicable: (c) => (c.fuselage.cabin_dp_Pa > 0 ? true : 'This aircraft is unpressurised (cabin pressure differential is 0): cabin altitude equals flight altitude. See the cabin air analysis.'),
  inputs: [
    num('alt_m', 'Flight altitude', 'm', 10668, 0, 20000, 'Schedule'), num('dp_max', 'Maximum cabin pressure differential', 'Pa', 57000, 1000, 100000, 'Schedule'),
    num('cab_alt_max', 'Cabin altitude at the ceiling', 'm', 2400, 0, 4500, 'Schedule', '2438 m (8000 ft) is the maximum in normal operation for transport aircraft (14 CFR 25.841(a)); newer composite fuselages use about 1800 m'), num('ceiling_m', 'Maximum operating altitude', 'm', 12500, 500, 20000, 'Schedule'),
    num('field_m', 'Departure field elevation', 'm', 0, -400, 4500, 'Schedule'), num('roc', 'Aircraft rate of climb', 'm/s', 10, 0.1, 60, 'Schedule'),
    num('V_cab', 'Pressurised volume', 'm³', 330, 0.5, 5000, 'Cabin'), num('T_cab', 'Cabin temperature', 'K', 295, 250, 320, 'Cabin'), num('m_in', 'Air supply (pack) flow', 'kg/s', 0.9, 0, 50, 'Cabin'),
    num('A_leak', 'Effective fuselage leakage area', 'm²', 8e-4, 0, 0.5, 'Cabin', 'Structural leakage through seals and drains'),
    num('r_fus', 'Fuselage radius', 'm', 1.975, 0.1, 5, 'Structure'), num('t_skin_mm', 'Skin thickness', 'mm', 1.6, 0.2, 20, 'Structure'), num('sigma_allow', 'Allowable hoop stress', 'Pa', 100e6, 1e7, 1e9, 'Structure', 'Fatigue-driven: about 80–110 MPa for aluminium pressure cabins (typical)'),
    num('d_hole', 'Decompression opening equivalent diameter', 'm', 0.3, 0.005, 3, 'Decompression', 'A cabin window ≈ 0.3 m; a door seal failure is much smaller'), num('Cd', 'Opening discharge coefficient', '-', 0.75, 0.3, 1, 'Decompression'),
    sel('process', 'Cabin air process', PROC, PROC[1], 'Decompression', 'Rapid decompression is close to adiabatic (fogging, cooling)'),
    num('t_end', 'Decompression time simulated', 's', 120, 0.5, 3600, 'Numerics'), num('nSteps', 'Time steps', '', 600, 20, 100000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => {
    const pax = c.mission.pax || 0, len = c.fuselage.len_m || 5, dia = c.fuselage.dia_m || 1.5, V = 0.75 * (Math.PI / 4) * dia * dia * len, alt = cruise(c).alt;
    return { alt_m: alt, dp_max: c.fuselage.cabin_dp_Pa, cab_alt_max: c.systems.cabin_alt_m || 2400, ceiling_m: Math.max(alt * 1.12, alt + 500), field_m: c.site?.elev_m ?? 0, V_cab: V, m_in: Math.max(0.01, 0.0055 * (pax + 5)), A_leak: 2.5e-6 * V, r_fus: dia / 2, t_skin_mm: Math.max(0.8, Math.min(c.struct.t_skin_mm, 1e3 * (c.fuselage.cabin_dp_Pa * dia) / 2 / 85e6)), d_hole: N.clamp(0.076 * dia, 0.1, 0.4), roc: c.mass.mtow_kg > 30000 ? 10 : 7 };
  },
  run(i) {
    const s = schedule(i, i.alt_m), mLeak = orifice(i.A_leak, s.pc, i.T_cab, s.pa), mOfv = i.m_in - mLeak, Aofv = mOfv > 0 ? mOfv / orifice(1, s.pc, i.T_cab, s.pa) : 0, d = decomp(i, s), warnings = [];
    const hs = N.linspace(i.field_m, Math.max(i.ceiling_m, i.alt_m), 60), sc = hs.map((h) => schedule(i, h)), dh = 50, rate = ((schedule(i, i.alt_m + dh / 2).hc - schedule(i, i.alt_m - dh / 2).hc) / dh) * i.roc, rMax = N.amax(sc.slice(1).map((v, j) => ((v.hc - sc[j].hc) / (hs[j + 1] - hs[j])) * i.roc));
    const hoop = (s.dp * i.r_fus) / (i.t_skin_mm / 1e3), up = (lim) => { for (let k = 1; k < d.t.length; k++) if (d.alt[k] >= lim && d.alt[k - 1] < lim) return d.t[k - 1] + ((d.t[k] - d.t[k - 1]) * (lim - d.alt[k - 1])) / (d.alt[k] - d.alt[k - 1]); return d.alt[0] >= lim ? 0 : NaN; };
    const t10 = up(3048), t15 = up(4572), t25 = up(7620), n = d.t.length - 1, tau = d.m0 / Math.max(1e-12, orifice(d.A, s.pc, i.T_cab, s.pa));
    if (s.limited) warnings.push(`The scheduled cabin altitude cannot be held at this flight altitude: the differential limit governs and cabin altitude rises to ${s.hc.toFixed(0)} m.`);
    if (mOfv < 0) warnings.push('Fuselage leakage exceeds the air supply: the cabin cannot hold pressure at this altitude with the stated inflow.');
    if (s.hc > 2440) warnings.push('Cabin altitude exceeds 2438 m (8000 ft), the maximum for normal operation of transport aircraft (14 CFR 25.841(a)).');
    if (rMax > 2.6) warnings.push(`Cabin climb rate reaches ${rMax.toFixed(1)} m/s during the climb: above the roughly 2.5 m/s (500 ft/min) comfort guideline.`);
    const nan = (v) => (Number.isNaN(v) ? i.t_end : v), nt = (v) => (Number.isNaN(v) ? 'Not reached within the simulated time' : '');
    return {
      kpis: [
        kp('cabin_alt_m', 'Cabin altitude', s.hc, 'm', s.hc <= 2440 ? 'ok' : 'warn', 'Transport-category maximum in normal operation: 2438 m (8000 ft)'), kp('cabin_p_Pa', 'Cabin pressure', s.pc, 'Pa'), kp('cabin_dp_Pa', 'Cabin pressure differential', s.dp, 'Pa', s.limited ? 'warn' : 'ok', `Limit ${(i.dp_max / 1e3).toFixed(1)} kPa`),
        kp('leak_flow_kgs', 'Fuselage leakage flow', mLeak, 'kg/s'), kp('ofv_flow_kgs', 'Outflow-valve flow', mOfv, 'kg/s', mOfv > 0 ? 'ok' : 'bad'), kp('ofv_area_m2', 'Outflow-valve effective area', Aofv, 'm²'),
        kp('cabin_rate_ms', 'Cabin climb rate at this point of the climb', rate, 'm/s'), kp('cabin_rate_max_ms', 'Largest cabin climb rate over the climb', rMax, 'm/s', rMax <= 2.6 ? 'ok' : 'warn', 'Comfort guideline ≈ 2.5 m/s climb, 1.5 m/s descent'),
        kp('hoop_stress_Pa', 'Fuselage skin hoop stress', hoop, 'Pa', hoop <= i.sigma_allow ? 'ok' : 'warn'),
        kp('t_to_3048m_s', 'Decompression: time to 3048 m (10 000 ft) cabin altitude', nan(t10), 's', undefined, nt(t10)), kp('t_to_4572m_s', 'Decompression: time to 4572 m (15 000 ft)', nan(t15), 's', undefined, nt(t15)), kp('t_to_7620m_s', 'Decompression: time to 7620 m (25 000 ft)', nan(t25), 's', undefined, nt(t25)),
        kp('decomp_tau_s', 'Initial decompression time constant', tau, 's'), kp('cabin_alt_end_m', 'Cabin altitude at end of the simulated time', d.alt[n], 'm'), kp('p_cabin_end_Pa', 'Cabin pressure at end of the simulated time', d.p[n], 'Pa'), kp('T_cabin_end_K', 'Cabin air temperature at end', d.T[n], 'K'),
      ],
      plots: [
        { type: 'line', title: 'Cabin altitude schedule', xlabel: 'Flight altitude [m]', ylabel: 'Cabin altitude [m]', series: [{ name: 'Cabin altitude', x: hs, y: sc.map((v) => v.hc) }, { name: 'Unpressurised', x: hs, y: hs, style: 'dash' }], annotations: [{ x: i.alt_m, label: 'Flight altitude' }] },
        { type: 'line', title: 'Cabin pressure differential', xlabel: 'Flight altitude [m]', ylabel: 'Differential [kPa]', series: [{ name: 'Cabin − ambient', x: hs, y: sc.map((v) => v.dp / 1e3) }], annotations: [{ y: i.dp_max / 1e3, label: 'Maximum' }] },
        { type: 'line', title: 'Cabin altitude after the opening appears', xlabel: 'Time [s]', ylabel: 'Cabin altitude [m]', series: [{ name: 'Cabin altitude', x: thin(d.t), y: thin(d.alt) }], annotations: [{ y: 3048, label: '3048 m' }, { y: 4572, label: '4572 m' }, { y: i.alt_m, label: 'Flight altitude' }] },
      ],
      warnings, models: ['Proportional cabin-altitude schedule with differential-pressure limit', 'Compressible (choked/unchoked) orifice flow for leakage, outflow valve and decompression', 'Lumped cabin mass balance, RK4', 'Thin-wall hoop stress Δp·r/t'],
      assumptions: ['Standard-atmosphere ambient pressure', 'Uniform cabin air; packs continue at constant flow during the decompression, outflow valve closed', 'Aircraft holds altitude during the decompression (no emergency descent)', 'Hoop stress ignores frames, stringers and cut-outs', 'The opening is postulated by the input; leakage area, allowable hoop stress, pack flow and the comfort rate limits are typical values'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [50, 100, 200, 400, 800], metric: 'p_cabin_end_Pa' },
  verify() {
    const b = { alt_m: 12000, dp_max: 57000, cab_alt_max: 2400, ceiling_m: 12500, field_m: 0, roc: 10, V_cab: 100, T_cab: 295, m_in: 0, A_leak: 0, r_fus: 2, t_skin_mm: 2, sigma_allow: 1e8, d_hole: 0.1, Cd: 1, process: PROC[0], t_end: 8, nSteps: 400 };
    const s = schedule(b, 12000), d = decomp(b, s), A = (Math.PI * 0.01) / 4, Kc = Math.sqrt(g) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1))), tau = 100 / (A * Kc * Math.sqrt(R_AIR * 295)), lo = schedule({ ...b, dp_max: 30000 }, 12000);
    return [
      N.check('Isothermal choked decompression p = p0·exp(−t/τ)', d.p[400], s.pc * Math.exp(-8 / tau), 1e-7, 'Analytical solution while choked'),
      N.check('Differential limit: cabin pressure = ambient + Δp_max', lo.pc, isa(12000).p + 30000, 1e-12, 'Schedule cap'),
      N.check('Cabin altitude consistent with cabin pressure', isa(lo.hc).p, lo.pc, 1e-6, 'Standard atmosphere inversion'),
      N.check('Hoop stress Δp·r/t', N.kv(pressurisation.run(b)).hoop_stress_Pa, (s.dp * 2) / 0.002, 1e-12, 'Thin-wall cylinder'),
    ];
  },
  calibration: { params: [{ key: 'A_leak', min: 0, max: 0.05 }, { key: 'Cd', min: 0.4, max: 1 }], sweep: 't_end', target: 'p_cabin_end_Pa', note: 'Supply measured cabin pressure against time from a ground leak-down test (packs off) to fit the effective leakage area.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    // the decompression is a postulated failure: it sizes the oxygen system and the emergency descent, and only an excursion that no descent can prevent is a design finding
    if (o.cabin_alt_end_m > 12192) out.push({ severity: 'warn', title: 'Cabin altitude would exceed 12 190 m (40 000 ft) after the decompression', detail: `With a ${i.d_hole.toFixed(2)} m opening the cabin reaches ${o.cabin_alt_end_m.toFixed(0)} m, passing 7620 m after ${o.t_to_7620m_s.toFixed(0)} s.`, action: 'Limit the largest credible opening (smaller windows, structural containment) or the cruise altitude: an emergency descent cannot prevent an excursion this fast. Carry the event into Suite 22.', basis: 'Decompression mass balance; transport-category practice keeps the cabin below 12 190 m (40 000 ft) at all times after a credible failure' });
    else if (o.t_to_4572m_s < i.t_end) out.push({ severity: 'advise', title: 'Decompression sets the oxygen and emergency-descent requirement', detail: `With a ${i.d_hole.toFixed(2)} m opening and the aircraft holding altitude, cabin altitude passes 3048 m after ${o.t_to_3048m_s.toFixed(0)} s and 4572 m after ${o.t_to_4572m_s.toFixed(0)} s${o.t_to_7620m_s < i.t_end ? `, and 7620 m after ${o.t_to_7620m_s.toFixed(0)} s` : ''}.`, action: 'Size oxygen deployment and the emergency-descent profile to these times, so that the cabin is back below 7620 m (25 000 ft) within about two minutes; a smaller opening or a lower cruise altitude lengthens them. Carry the event into Suite 22.', basis: 'Decompression mass balance for a postulated opening' });
    if (o.hoop_stress_Pa > i.sigma_allow) out.push({ severity: 'warn', title: 'Hoop stress exceeds the fatigue allowable', detail: `${(o.hoop_stress_Pa / 1e6).toFixed(0)} MPa per pressurisation cycle.`, action: 'Increase skin gauge or reduce the differential; pass the pressure cycle to Suite 9 for fatigue and crack-growth life.', basis: 'Pressure-cabin fatigue practice' });
    out.push({ severity: 'info', title: 'Cabin altitude versus structure trade', detail: `Cabin altitude ${o.cabin_alt_m.toFixed(0)} m at ${(o.cabin_dp_Pa / 1e3).toFixed(1)} kPa differential.`, action: 'A lower cabin altitude improves passenger comfort but raises the differential, skin stress and structural mass (fuel); the leakage flow of ' + o.leak_flow_kgs.toFixed(2) + ' kg/s is conditioned air thrown away, so tighter sealing saves bleed or compressor power.', basis: 'Pressurisation design trade' });
    return out;
  },
};

// ---- 7. air-conditioning pack -----------------------------------------------------------------
const PACK = ['Bleed air-cycle (bootstrap)', 'Electric air-cycle (bootstrap, ram-air compressor)', 'Vapour-cycle (electric)'];
/** Bootstrap air-cycle machine: primary HX → compressor → secondary HX → turbine, compressor driven by the turbine. */
function bootstrap(i, p1, T1, Tsink) {
  const T2 = T1 - i.eps1 * (T1 - Tsink), p2 = p1 * (1 - i.dp_hx), p5 = i.p_cab * 1.04;
  const st = (PR) => { const T3 = T2 * (1 + (PR ** KG - 1) / i.eta_c), p3 = p2 * PR, T4 = T3 - i.eps2 * (T3 - Tsink), p4 = p3 * (1 - i.dp_hx), er = Math.min(1, p5 / p4), T5 = T4 * (1 - i.eta_t * (1 - er ** KG)); return { T3, p3, T4, p4, T5, res: i.eta_m * (T4 - T5) - (T3 - T2) }; };
  const PRc = st(1).res <= 0 ? 1 : st(8).res > 0 ? 8 : N.brent((x) => st(x).res, 1, 8, 1e-12), s = st(PRc);
  return { T1, p1, T2, p2, ...s, p5, PRc };
}
function packCalc(i) {
  const a = isa(i.alt_m, i.dISA), V = i.M * a.a, Tr = a.T * (1 + 0.2 * i.M ** 2), pr = a.p * (1 + 0.2 * i.M ** 2) ** 3.5 * i.eta_ram, k = KG;
  const elecFuel = (P) => (P / i.eta_gen) * i.sfc_shaft, drag = (m) => m * V * i.tsfc;
  // bleed architecture: air taken at port pressure, throttled and pre-cooled to the pack inlet state
  const pPort = Math.max(i.p_bleed * i.port_ratio, pr * 1.01), Tport = Tr * (1 + ((pPort / pr) ** k - 1) / i.eta_eng_c), bl = bootstrap(i, i.p_bleed, Math.min(Tport, i.T_precool), Tr), Pbleed = i.m_air * CP_AIR * (Tport - Tr);
  // electric architecture: dedicated compressor raises ram air only to the pressure the pack needs
  const pE = Math.max(i.p_cab * i.pr_elec, pr * 1.01), TE = Tr * (1 + ((pE / pr) ** k - 1) / i.eta_cac), el = bootstrap(i, pE, TE, Tr), Pel = (i.m_air * CP_AIR * (TE - Tr)) / i.eta_motor;
  const sup = (b) => { // water: saturated outlet with latent heat release, then limited to 275 K by the anti-ice bypass
    const w = i.w_in, f = (T) => CP_AIR * (T - b.T5) - 2.5e6 * Math.max(0, w - humRatio(pSat(T), b.p5)), Tw = f(b.T5) >= 0 ? b.T5 : N.brent(f, b.T5, b.T5 + 150, 1e-9), Ts = Math.max(Tw, 275);
    return { Tw, Ts, cond: Math.max(0, w - humRatio(pSat(Tw), b.p5)), Q: i.m_air * CP_AIR * (i.T_cab - Ts) };
  };
  const sb = sup(bl), se = sup(el), fuelB = Pbleed * i.sfc_shaft + drag((i.ram_ratio + 1) * i.m_air), fuelE = elecFuel(Pel) + drag((i.ram_ratio + 1) * i.m_air);
  // vapour cycle: fraction of the Carnot coefficient of performance
  const cop = (i.cop_frac * i.T_evap) / Math.max(1, i.T_cond - i.T_evap), Pvc = i.Q_cool / cop + 0.03 * i.Q_cool, fuelV = elecFuel(Pvc) + drag(i.m_air);
  return { a, V, Tr, pr, pPort, Tport, bl, el, sb, se, Pbleed, Pel, fuelB, fuelE, cop, Pvc, fuelV };
}
const aircycle = {
  id: 'aircycle', title: 'Air-conditioning pack: air-cycle machine and bleed versus electric supply', fidelity: 'reduced-order',
  summary: 'Thermodynamic model of a bootstrap air-cycle pack (heat exchangers, compressor, turbine, water condensation), its cooling capacity, and the fuel cost of supplying it with engine bleed, an electric compressor or a vapour-cycle unit.',
  equations: ['Conservation of mass, momentum and energy equations', 'Heat exchanger equations', 'Ideal gas law', 'Psychrometric relations', 'Incompressible and compressible fluid flow equations'],
  applicable: (c) => ((c.mission.pax || 0) > 0 ? true : 'No occupied cabin: equipment cooling by ventilation air is covered by the cabin and equipment-bay analysis.'),
  inputs: [
    sel('mode', 'Pack architecture', PACK, PACK[0], 'Architecture'), num('Q_cool', 'Cooling duty required', 'W', 45000, 0, 5e6, 'Requirement', 'From the cabin heat-load analysis'),
    num('m_air', 'Pack air flow', 'kg/s', 0.9, 1e-4, 50, 'Requirement'), num('T_cab', 'Cabin temperature', 'K', 297, 270, 320, 'Requirement'), num('p_cab', 'Cabin pressure', 'Pa', 75300, 20000, 110000, 'Requirement'),
    num('alt_m', 'Flight altitude', 'm', 10668, 0, 20000, 'Flight condition'), num('M', 'Mach number', '-', 0.78, 0, 0.95, 'Flight condition'), num('dISA', 'ISA deviation', 'K', 0, -60, 50, 'Flight condition'),
    num('w_in', 'Humidity ratio of the intake air', 'kg/kg', 0.0005, 0, 0.04, 'Flight condition', 'Near zero at cruise; 0.01–0.02 on a humid ground day'), num('eta_ram', 'Ram pressure recovery', '-', 0.9, 0.5, 1, 'Flight condition'),
    num('p_bleed', 'Regulated bleed pressure at the pack', 'Pa', 250e3, 5e4, 1e6, 'Bleed supply'), num('port_ratio', 'Engine port pressure / regulated pressure', '-', 1.5, 1, 5, 'Bleed supply', 'Bleed is taken at a higher pressure than needed and throttled'),
    num('T_precool', 'Pre-cooler outlet temperature', 'K', 470, 300, 600, 'Bleed supply'), num('eta_eng_c', 'Engine compressor efficiency to the port', '-', 0.86, 0.6, 0.95, 'Bleed supply'),
    num('pr_elec', 'Electric compressor delivery / cabin pressure', '-', 2.4, 1.1, 6, 'Electric supply'), num('eta_cac', 'Cabin air compressor efficiency', '-', 0.78, 0.5, 0.9, 'Electric supply'), num('eta_motor', 'Motor and drive efficiency', '-', 0.92, 0.6, 0.99, 'Electric supply'),
    num('eta_gen', 'Generation and distribution efficiency', '-', 0.9, 0.5, 0.99, 'Electric supply'),
    num('eps1', 'Primary heat-exchanger effectiveness', '-', 0.8, 0.3, 0.98, 'Air-cycle machine', 'From the thermal suite when available'), num('eps2', 'Secondary heat-exchanger effectiveness', '-', 0.85, 0.3, 0.98, 'Air-cycle machine'),
    num('eta_c', 'ACM compressor efficiency', '-', 0.75, 0.4, 0.9, 'Air-cycle machine'), num('eta_t', 'ACM turbine efficiency', '-', 0.8, 0.4, 0.92, 'Air-cycle machine'), num('eta_m', 'ACM mechanical efficiency', '-', 0.95, 0.6, 1, 'Air-cycle machine'),
    num('dp_hx', 'Pressure loss fraction per heat exchanger', '-', 0.03, 0, 0.2, 'Air-cycle machine'), num('ram_ratio', 'Ram cooling air / pack air', '-', 2, 0, 10, 'Air-cycle machine'),
    num('cop_frac', 'Vapour-cycle COP / Carnot COP', '-', 0.45, 0.1, 0.8, 'Vapour cycle'), num('T_evap', 'Evaporator temperature', 'K', 278, 240, 300, 'Vapour cycle'), num('T_cond', 'Condenser temperature', 'K', 328, 280, 380, 'Vapour cycle'),
    num('sfc_shaft', 'Fuel per unit shaft energy', 'kg/J', 5.2e-8, 0, 3e-7, 'Fuel penalty', '1/(η_thermal·LHV): ≈ 5·10⁻⁸ for a modern turbofan core, ≈ 8·10⁻⁸ for small engines; 0 for battery aircraft'),
    num('tsfc', 'Thrust-specific fuel consumption', 'kg/N/s', 1.6e-5, 0, 1e-4, 'Fuel penalty', 'For the ram-air momentum drag'),
  ],
  defaults: (c, up, d) => {
    const { alt, a, M } = cruise(c), pax = c.mission.pax || 0, turbine = ['turbofan', 'turbojet', 'turboprop', 'turboshaft'].includes(c.prop.type), elec = c.prop.type === 'electric', press = c.fuselage.cabin_dp_Pa > 0;
    const pc = press ? Math.max(isa(c.systems.cabin_alt_m || 2400).p, a.p) : a.p, Q = up.fuelecs?.ecs_cooling_W > 0 ? up.fuelecs.ecs_cooling_W : 130 * (pax + 2) + 800, m = Math.max(0.02, 0.0055 * (pax + 3), (1.1 * Q) / (CP_AIR * 22));
    const tsfc = elec ? 0 : up.propulsion?.tsfc_kg_Ns ?? (d.T_total > 0 ? c.prop.tsfc_kg_Ns : (c.prop.bsfc_kg_Ws * Math.max(c.flight.V_ms, 20)) / c.prop.eta_prop);
    return { mode: turbine && c.systems.bleed ? PACK[0] : turbine ? PACK[1] : PACK[2], Q_cool: Q, m_air: m, p_cab: pc, alt_m: alt, M, dISA: c.atm.dISA_K, eps1: up.thermal?.hx_effectiveness ? N.clamp(up.thermal.hx_effectiveness, 0.5, 0.9) : undefined, sfc_shaft: elec ? 0 : turbine ? (up.propulsion?.eta_thermal ? 1 / (up.propulsion.eta_thermal * 43.15e6) : 5.2e-8) : c.prop.bsfc_kg_Ws, tsfc, T_cond: Math.max(siteT(c), 288) + 25 };
  },
  run(i, ctx) {
    const r = packCalc(i), vc = i.mode === PACK[2], bleed = i.mode === PACK[0], b = bleed ? r.bl : r.el, s = bleed ? r.sb : r.se, warnings = [], elecAc = ctx?.case?.prop?.type === 'electric';
    const Qcap = vc ? i.Q_cool : s.Q, Pin = vc ? r.Pvc : bleed ? r.Pbleed : r.Pel, fuel = vc ? r.fuelV : bleed ? r.fuelB : r.fuelE, mReq = vc ? i.Q_cool / (CP_AIR * Math.max(1, i.T_cab - (i.T_evap + 5))) : i.Q_cool / Math.max(1e-9, CP_AIR * (i.T_cab - s.Ts));
    if (!vc && s.Tw < 275) warnings.push(`Turbine outlet would be ${s.Tw.toFixed(0)} K: bypass (trim) air is mixed in to hold 275 K and prevent icing of the water separator, so not all of the expansion cooling is used.`);
    if (!vc && Qcap < i.Q_cool) warnings.push('Pack cooling capacity at this flow is below the required duty: increase the pack flow or heat-exchanger effectiveness.');
    if (!vc && b.PRc >= 8) warnings.push('Bootstrap compressor pressure ratio reached the solver limit of 8: check the component efficiencies.');
    if (!vc && s.cond > 0) warnings.push(`${(1e3 * s.cond).toFixed(1)} g of water per kg of air condenses at the turbine outlet and must be removed by the water separator.`);
    if (elecAc) warnings.push('Battery-electric aircraft: the fuel penalty is zero; use the electrical input power as the load on the battery in Suite 19.');
    const kpis = [
      kp('ecs_cooling_W', 'Cooling delivered to the cabin', Math.min(Qcap, Math.max(i.Q_cool, 0)) || Qcap, 'W', Qcap >= i.Q_cool ? 'ok' : 'warn', `Capacity ${Qcap.toFixed(0)} W, required ${i.Q_cool.toFixed(0)} W`),
      kp('pack_capacity_W', 'Pack cooling capacity at the stated flow', Qcap, 'W'), kp('bleed_kgs', 'Engine bleed flow', bleed ? i.m_air : 0, 'kg/s'), kp('flow_required_kgs', 'Air flow needed for the duty', mReq, 'kg/s'),
      kp('ecs_power_W', bleed ? 'Engine compression power embodied in the bleed' : 'Electrical input power', Pin, 'W'), kp('cop', 'Coefficient of performance', Pin > 0 ? Qcap / Pin : 0, '-'),
      kp('fuel_penalty_kg_h', 'Fuel penalty of the selected architecture', fuel * 3600, 'kg/h'),
    ];
    if (!vc) kpis.push(kp('T_supply_K', 'Pack outlet (supply) temperature', s.Ts, 'K'), kp('T_turbine_out_K', 'Turbine outlet temperature before trim', s.Tw, 'K'), kp('acm_PR', 'Bootstrap compressor pressure ratio', b.PRc, '-'), kp('water_g_per_kg', 'Water condensed', 1e3 * s.cond, 'g/kg'),
      kp('fuel_bleed_kg_h', 'Fuel penalty with engine bleed', r.fuelB * 3600, 'kg/h'), kp('fuel_electric_kg_h', 'Fuel penalty with an electric compressor', r.fuelE * 3600, 'kg/h'),
      kp('electric_saving_pct', 'Fuel saving of electric over bleed supply', r.fuelB > 0 ? 100 * (1 - r.fuelE / r.fuelB) : 0, '%'), kp('bootstrap_residual', 'Turbine–compressor work balance residual', Math.abs(b.res) / Math.max(1, b.T3 - b.T2), '-'));
    else kpis.push(kp('cop_carnot', 'Carnot COP between evaporator and condenser', i.T_evap / Math.max(1, i.T_cond - i.T_evap), '-'), kp('heat_rejected_W', 'Heat rejected at the condenser', i.Q_cool + r.Pvc, 'W'));
    const st = ['Pack inlet', 'After primary HX', 'After compressor', 'After secondary HX', 'Turbine outlet'], plots = vc
      ? [{ type: 'bar', title: 'Vapour-cycle energy flows', ylabel: 'Power [kW]', categories: ['Cooling duty', 'Electrical input', 'Heat rejected'], series: [{ name: 'Power', y: [i.Q_cool / 1e3, r.Pvc / 1e3, (i.Q_cool + r.Pvc) / 1e3] }] }]
      : [{ type: 'line', title: 'Air temperature through the pack', xlabel: 'Station [-]', ylabel: 'Temperature [K]', series: [{ name: 'Bleed-supplied pack', x: [1, 2, 3, 4, 5], y: [r.bl.T1, r.bl.T2, r.bl.T3, r.bl.T4, r.sb.Tw], style: 'line+points' }, { name: 'Electric-compressor pack', x: [1, 2, 3, 4, 5], y: [r.el.T1, r.el.T2, r.el.T3, r.el.T4, r.se.Tw], style: 'line+points' }], annotations: [{ y: r.Tr, label: 'Ram air (heat sink)' }, { y: i.T_cab, label: 'Cabin' }] },
        { type: 'bar', title: 'Fuel penalty by air-supply architecture', ylabel: 'Fuel [kg/h]', categories: ['Engine bleed', 'Electric compressor'], stacked: true, series: [{ name: 'Compression power', y: [r.Pbleed * i.sfc_shaft * 3600, (r.Pel / i.eta_gen) * i.sfc_shaft * 3600] }, { name: 'Intake momentum drag (cabin + cooling air)', y: [1, 1].map(() => (i.ram_ratio + 1) * i.m_air * r.V * i.tsfc * 3600) }] }];
    return {
      kpis, plots, warnings,
      tables: vc ? [] : [{ title: 'Pack stations (selected architecture)', columns: ['Station', 'Temperature [K]', 'Pressure [kPa]'], rows: st.map((nm, j) => [nm, [b.T1, b.T2, b.T3, b.T4, b.T5][j], [b.p1, b.p2, b.p3, b.p4, b.p5][j] / 1e3]) }],
      models: vc ? ['Vapour-compression cycle as a fraction of the Carnot COP', 'Fuel penalty from electrical power and intake momentum drag'] : ['Bootstrap air-cycle: isentropic-efficiency compressor and turbine, effectiveness heat exchangers', 'Turbine–compressor power balance solved for the compressor pressure ratio', 'Saturated-outlet water condensation with latent heat', 'Fuel penalty: compression work × shaft SFC plus ram-air momentum drag × TSFC'],
      assumptions: ['Calorically perfect air (cp = 1004.7 J/kg/K, γ = 1.4)', 'Heat exchangers reject to ram air at total temperature; fan power neglected', 'Bleed fuel penalty approximated as the engine compression work embodied in the bleed air at a constant shaft SFC', 'All intake momentum of the cabin air and heat-exchanger cooling air is lost (no thrust recovery), for both supply architectures', 'Steady state at one flight condition', 'Component efficiencies, heat-exchanger effectiveness, bleed pressures, ram-air ratio and specific fuel consumptions are typical values, not data for a specific pack or engine'],
    };
  },
  verify() {
    const b = { mode: PACK[0], Q_cool: 30000, m_air: 0.5, T_cab: 297, p_cab: 80000, alt_m: 10000, M: 0.8, dISA: 0, w_in: 0, eta_ram: 1, p_bleed: 250e3, port_ratio: 1.5, T_precool: 470, eta_eng_c: 0.86, pr_elec: 2.4, eta_cac: 0.78, eta_motor: 0.92, eta_gen: 0.9, eps1: 0.8, eps2: 0.85, eta_c: 0.75, eta_t: 0.8, eta_m: 0.95, dp_hx: 0.03, ram_ratio: 2, cop_frac: 1, T_evap: 280, T_cond: 320, sfc_shaft: 5e-8, tsfc: 1.6e-5 };
    const r = packCalc(b), s = r.bl, a = isa(10000), Tr = a.T * 1.128, pr = a.p * 1.128 ** 3.5, id = bootstrap({ ...b, eps1: 1, eps2: 1, eta_c: 1, eta_t: 1, eta_m: 1, dp_hx: 0 }, 250e3, 470, 250);
    return [
      N.check('Bootstrap work balance η_m·(T4 − T5) = T3 − T2', b.eta_m * (s.T4 - s.T5), s.T3 - s.T2, 1e-9, 'Turbine drives the compressor'),
      N.check('Ideal pack (ε = 1, η = 1): turbine work equals compressor work, PR satisfies PR^k − 1 = 1 − (p5/(p1·PR))^k', id.PRc ** KG - 1, 1 - (id.p5 / (250e3 * id.PRc)) ** KG, 1e-9, 'Isentropic relations with both streams cooled to the sink temperature'),
      N.check('Electric compressor power m·cp·T0·(PR^k − 1)/(η_c·η_m)', r.Pel, (0.5 * CP_AIR * Tr * (((2.4 * 80000) / pr) ** KG - 1)) / 0.78 / 0.92, 1e-10, 'Adiabatic compression'),
      N.check('Vapour-cycle COP at Carnot fraction 1', r.cop, 280 / 40, 1e-12, 'Carnot refrigerator'),
    ];
  },
  calibration: { params: [{ key: 'eps2', min: 0.4, max: 0.98 }, { key: 'eta_t', min: 0.5, max: 0.92 }, { key: 'eta_c', min: 0.5, max: 0.9 }], sweep: 'm_air', target: 'pack_capacity_W', note: 'Supply measured pack cooling capacity against pack flow from an environmental-control rig or altitude-chamber test.' },
  recommend(res, i, ctx) {
    const o = res.outputs, out = [];
    if (i.mode !== PACK[2] && o.electric_saving_pct > 5 && i.sfc_shaft > 0) out.push({ severity: i.mode === PACK[0] ? 'advise' : 'info', title: 'Bleed-less (electric) air supply saves fuel', detail: `Electric compression costs ${o.fuel_electric_kg_h.toFixed(0)} kg/h against ${o.fuel_bleed_kg_h.toFixed(0)} kg/h for engine bleed (${o.electric_saving_pct.toFixed(0)}% less), because air is compressed only to the pressure the pack needs instead of being throttled and pre-cooled.`, action: 'Consider electrically driven cabin air compressors; include the larger generators and their mass in Suite 19 and the mission fuel in Suite 24 before deciding.', basis: 'Compression work comparison at equal pack flow' });
    if (i.mode === PACK[2] && o.cop > 0) out.push({ severity: 'info', title: 'Vapour-cycle cooling', detail: `COP ${o.cop.toFixed(1)}: ${(o.ecs_power_W / 1e3).toFixed(1)} kW of electricity for ${(i.Q_cool / 1e3).toFixed(1)} kW of cooling.`, action: 'On battery aircraft this power comes directly out of range: pre-condition the cabin on ground power and use recirculation to cut the duty.', basis: 'Refrigeration COP' });
    if (o.pack_capacity_W < i.Q_cool) out.push({ severity: 'warn', title: 'Pack capacity is below the cabin cooling load', detail: `${(o.pack_capacity_W / 1e3).toFixed(1)} kW available for ${(i.Q_cool / 1e3).toFixed(1)} kW required.`, action: `Raise the pack flow to about ${o.flow_required_kgs.toFixed(2)} kg/s or improve the secondary heat exchanger.`, basis: 'Sensible cooling m·cp·(T_cabin − T_supply)' });
    return out;
  },
};

export default {
  id: 'fuelecs', n: 20,
  tagline: 'Whether fuel reaches the engines with margin, where the fuel mass moves and how cold it gets, and what it takes to keep the cabin pressurised, ventilated and comfortable.',
  analyses: [feed, cg, slosh, fueltherm, cabin, pressurisation, aircycle],
  consumes: [
    { from: 'performance', keys: ['fuel_flow_cruise_kgs'], why: 'Cruise fuel burn for CG travel and tank cold-soak' },
    { from: 'propulsion', keys: ['fuel_flow_kgs', 'tsfc_kg_Ns', 'eta_thermal'], why: 'Fuel demand and the fuel cost of bleed, power and drag' },
    { from: 'flightdyn', keys: ['dr_omega_rads'], why: 'Dutch-roll frequency to keep clear of the slosh frequency' },
    { from: 'thermal', keys: ['T_recovery_K', 'T_soak_K', 'hx_effectiveness'], why: 'Skin temperatures and heat-exchanger effectiveness' },
    { from: 'electrical', keys: ['elec_losses_W'], why: 'Equipment-bay heat load' },
  ],
  provides: [
    { key: 'fuel_cg_shift_m', label: 'CG travel due to fuel', unit: 'm' }, { key: 'slosh_freq_Hz', label: 'First slosh frequency', unit: 'Hz' }, { key: 'fuel_pump_power_W', label: 'Fuel boost pump power', unit: 'W' },
    { key: 'cabin_alt_m', label: 'Cabin altitude', unit: 'm' }, { key: 'ecs_cooling_W', label: 'ECS cooling load', unit: 'W' }, { key: 'bleed_kgs', label: 'Bleed air flow', unit: 'kg/s' },
  ],
  handoff: [
    { model: 'Volume-of-fluid free-surface sloshing (3-D, with baffles, roof impact and breaking waves)', why: 'Needs a 3-D two-phase Navier–Stokes solver on the real tank geometry; here sloshing is linear potential theory plus 1-D shallow water', tool: 'VOF or SPH CFD solver' },
    { model: 'Cabin ventilation CFD (airflow patterns, thermal comfort, contaminant transport)', why: 'Requires a 3-D turbulent buoyant flow solution of the cabin interior; the native model is well-mixed', tool: 'RANS/LES CFD with thermal comfort post-processing' },
    { model: 'Coupled fuel sloshing–aircraft dynamics and tank structural–fluid coupling', why: 'The slosh force is computed for a prescribed tank motion and a rigid tank', tool: 'Flight-dynamics model with slosh pendulum states (Suite 4 hand-off) or coupled CFD–FE' },
    { model: 'Full fuel-system network with transfer, cross-feed, jet pumps, vent and refuel/defuel logic', why: 'One feed line is solved natively; a complete system needs a component-library network simulator', tool: '1-D thermal-fluid system simulation tools' },
    { model: 'Fuel vapour–liquid equilibrium and flammability exposure (multi-component distillation, oxygen evolution)', why: 'A single exponential vapour-pressure fit is used and vapour in the ullage is neglected', tool: 'Fleet-average flammability (Monte Carlo) and multi-component fuel property models' },
    { model: 'Detailed ECS network with turbomachinery maps, valves and control laws; cryogenic hydrogen fuel systems', why: 'The pack uses constant efficiencies and effectiveness at one operating point; hydrogen needs two-phase boil-off models', tool: 'ECS system simulation with component maps; cryogenic tank models' },
  ],
};
