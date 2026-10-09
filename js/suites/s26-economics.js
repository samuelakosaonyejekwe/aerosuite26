// Suite 26 — Aircraft Economics and Techno-Economic Analysis.
// Operating-cost build-up (direct and indirect, cash and ownership kept apart), investment appraisal on yearly
// cash flows (NPV, IRR, MIRR, payback), life-cycle cost with learning curves, lease-versus-buy and replacement
// timing, Monte Carlo financial risk with correlated drivers and a deferral option, route / mission / survey
// economics, and the economics of carbon and alternative energy carriers.
// A few defaults carry a sourced anchor, named in their help text (FAA Economic Values crew and maintenance cost per
// block hour, EIA / BTS jet-fuel price, EASA SAF price ratio, EU ETS allowance price, EIA electricity price and grid
// intensity, US corporate tax rate, IATA route price elasticity). Every other price, charge, rate and cost coefficient
// is an ILLUSTRATIVE placeholder. Replace them with operator accounts, lessor quotes, charge tables and supplier data.

import * as N from '../core/numerics.js';
import { isa, G0 } from '../core/atmosphere.js';
import { FLUIDS, BATTERIES } from '../data/materials.js';

// ---- money helpers --------------------------------------------------------------------------
const pvFactor = (r, n) => (Math.abs(r) < 1e-12 ? n : (1 - (1 + r) ** -n) / r);
/** Level payment that repays `principal` over n periods at rate r (annuity, payments in arrears). */
const annuity = (principal, r, n) => (n > 0 ? principal / pvFactor(r, n) : 0);
const npvOf = (cf, r) => cf.reduce((s, v, t) => s + v / (1 + r) ** t, 0);
/** Internal rate of return with a multiple-root guard. Returns { irr, roots, signChanges }. */
function irrOf(cf) {
  let sc = 0, last = 0; for (const v of cf) { if (v !== 0) { if (last && Math.sign(v) !== last) sc++; last = Math.sign(v); } }
  if (!sc) return { irr: NaN, roots: [], signChanges: 0 };
  const grid = [...N.linspace(-0.95, 1, 196), ...N.linspace(1.05, 10, 60)], roots = []; let f0 = npvOf(cf, grid[0]);
  for (let k = 1; k < grid.length; k++) { const f1 = npvOf(cf, grid[k]); if (Number.isFinite(f0) && Number.isFinite(f1) && f0 * f1 < 0) roots.push(N.brent((r) => npvOf(cf, r), grid[k - 1], grid[k], 1e-12)); f0 = f1; }
  const pos = roots.filter((r) => r > -0.5); return { irr: pos.length ? pos[0] : NaN, roots, signChanges: sc }; // roots below −50% are not reported as a rate of return
}
function mirrOf(cf, rFin, rRe) { const n = cf.length - 1; let pvNeg = 0, fvPos = 0; cf.forEach((v, t) => { if (v < 0) pvNeg += v / (1 + rFin) ** t; else fvPos += v * (1 + rRe) ** (n - t); }); return pvNeg < 0 && fvPos > 0 && n > 0 ? (fvPos / -pvNeg) ** (1 / n) - 1 : NaN; }
/** First time the cumulative series turns non-negative, with linear interpolation inside the year. */
function paybackOf(cf, r = 0) { let cum = 0; for (let t = 0; t < cf.length; t++) { const v = cf[t] / (1 + r) ** t, prev = cum; cum += v; if (t > 0 && cum >= 0 && prev < 0) return t - 1 + -prev / v; if (t === 0 && cum >= 0) return 0; } return Infinity; }
/** Depreciation fractions of the depreciable base per year. */
function depSchedule(method, life, table) {
  const n = Math.max(1, Math.round(life));
  if (method === 'declining balance') { const rate = 2 / n, out = []; let book = 1; for (let t = 0; t < n; t++) { const sl = book / (n - t), d = Math.max(rate * book, sl); out.push(d); book -= d; } return out; } // double-declining with switch to straight line
  if (method === 'custom table') { const v = String(table || '').split(/[\s,;]+/).filter(Boolean).map(Number).filter((x) => Number.isFinite(x) && x >= 0), s = N.sum(v); if (v.length && s > 0) return v.map((x) => x / s); }
  return new Array(n).fill(1 / n);
}
const learnB = (lr) => Math.log(lr) / Math.log(2);
/** Cumulative production cost of Q units. Wright: cumulative average = T1·Q^b. Crawford: unit n costs T1·n^b. */
const cumCost = (T1, Q, lr, model) => { const b = learnB(lr); if (model.startsWith('Wright')) return T1 * Q ** (1 + b); let s = 0; const n = Math.floor(Q); if (n <= 2000) { for (let k = 1; k <= n; k++) s += k ** b; return T1 * s; } return T1 * ((n + 0.5) ** (1 + b) - 0.5 ** (1 + b)) / (1 + b); };
/** American call on an asset paying a continuous yield, Cox–Ross–Rubinstein lattice. */
function crrCall(S, K, r, q, sigma, T, n) {
  if (!(sigma > 0) || !(T > 0)) return Math.max(S - K, 0);
  const dt = T / n, u = Math.exp(sigma * Math.sqrt(dt)), d = 1 / u, p = N.clamp((Math.exp((r - q) * dt) - d) / (u - d), 0, 1), disc = Math.exp(-r * dt), v = new Array(n + 1);
  for (let j = 0; j <= n; j++) v[j] = Math.max(S * u ** j * d ** (n - j) - K, 0);
  for (let s = n - 1; s >= 0; s--) for (let j = 0; j <= s; j++) v[j] = Math.max(disc * (p * v[j + 1] + (1 - p) * v[j]), S * u ** j * d ** (s - j) - K);
  return v[0];
}
const bsCall = (S, K, r, sigma, T) => { const d1 = (Math.log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * Math.sqrt(T)), d2 = d1 - sigma * Math.sqrt(T); return S * N.normCdf(d1) - K * Math.exp(-r * T) * N.normCdf(d2); };

// ---- shared inputs --------------------------------------------------------------------------
const BASE = [
  { key: 'price_usd', label: 'Aircraft price', unit: 'USD', default: 105e6, min: 0, group: 'Capital', help: 'Purchase price including buyer-furnished equipment' },
  { key: 'life_yr', label: 'Service life / appraisal horizon', unit: 'yr', default: 25, min: 1, max: 50, step: 1, discrete: true, group: 'Capital' },
  { key: 'residual_frac', label: 'Residual value at end of life', unit: '-', default: 0.15, min: 0, max: 0.9, group: 'Capital', help: 'Fraction of the purchase price in constant money' },
  { key: 'insurance_pct', label: 'Hull and liability insurance', unit: '%/yr of price', default: 0.35, min: 0, max: 15, group: 'Capital' },
  { key: 'debt_frac', label: 'Share of price financed by debt', unit: '-', default: 0.7, min: 0, max: 1, group: 'Capital' }, { key: 'loan_rate', label: 'Loan interest rate (nominal)', unit: '-', default: 0.055, min: 0, max: 0.4, group: 'Capital' }, { key: 'loan_yr', label: 'Loan term', unit: 'yr', default: 12, min: 1, max: 40, step: 1, discrete: true, group: 'Capital' },
  { key: 'util_fh_yr', label: 'Utilisation', unit: 'FH/yr', default: 3400, min: 1, max: 6500, group: 'Operation' }, { key: 'cycles_yr', label: 'Flights per year', unit: '1/yr', default: 1650, min: 1, group: 'Operation' },
  { key: 'block_kmh', label: 'Block speed', unit: 'km/h', default: 760, min: 1, group: 'Operation', help: 'Stage distance divided by flight time; from the mission simulation when available' },
  { key: 'mtow_t', label: 'Maximum take-off mass', unit: 't', default: 78, min: 0.0001, group: 'Operation', help: 'Basis of landing and navigation charges' },
  { key: 'seats', label: 'Seats', unit: '', default: 165, min: 0, step: 1, discrete: true, group: 'Revenue' }, { key: 'load_factor', label: 'Passenger load factor', unit: '-', default: 0.82, min: 0, max: 1, group: 'Revenue' },
  { key: 'yield_usd_pkm', label: 'Passenger yield', unit: 'USD/pax-km', default: 0.09, min: 0, group: 'Revenue' }, { key: 'cargo_kg', label: 'Cargo carried per flight', unit: 'kg', default: 0, min: 0, group: 'Revenue' }, { key: 'cargo_usd_tkm', label: 'Cargo yield', unit: 'USD/t-km', default: 0.25, min: 0, group: 'Revenue' },
  { key: 'charter_usd_fh', label: 'Charter / contract revenue per flight hour', unit: 'USD/FH', default: 0, min: 0, group: 'Revenue', help: 'For helicopters, UAV services and charter. For a commercial operation without a passenger yield a placeholder is filled: the rate at which the purchase just earns the discount rate, plus 5%. Replace it with the contract rate.' },
  { key: 'fuel_kg_fh', label: 'Fuel burn', unit: 'kg/FH', default: 2450, min: 0, group: 'Energy', help: 'Block fuel ÷ flight time from Suite 24 when available, so that taxi fuel is carried by the flight hours. For scale: US narrow-bodies averaged 822 US gal (about 2490 kg) per block hour in the FAA Economic Values tables.' }, { key: 'fuel_usd_kg', label: 'Fuel price', unit: 'USD/kg', default: 1.12, min: 0, group: 'Energy', help: 'Price of the fuel selected below, refreshed from live data when available. The default is jet fuel at 3.38 USD per US gallon, the EIA US Gulf Coast spot average for January–September 2026 (2025 averaged 0.70 USD/kg): a volatile ex-refinery price, not an into-plane price.' },
  { key: 'fuel', label: 'Fuel', type: 'select', options: ['Jet A-1', 'Avgas 100LL', 'SAF (HEFA-SPK)', 'Liquid hydrogen'], default: 'Jet A-1', group: 'Energy' },
  { key: 'energy_kWh_fh', label: 'Battery energy use', unit: 'kWh/FH', default: 0, min: 0, group: 'Energy' }, { key: 'elec_usd_kWh', label: 'Electricity price', unit: 'USD/kWh', default: 0.14, min: 0, group: 'Energy', help: 'US commercial-sector average was 0.145 USD/kWh in mid-2026 (EIA)' },
  { key: 'saf_blend', label: 'SAF blend share (by mass)', unit: '-', default: 0, min: 0, max: 1, group: 'Energy' }, { key: 'saf_price_ratio', label: 'SAF price ÷ fossil fuel price', unit: '-', default: 2.84, min: 0.5, max: 12, group: 'Energy', help: 'EASA 2024 reference prices: aviation biofuel 2085 EUR/t against 734 EUR/t for conventional fuel (2.84). Synthetic e-fuels were about 10.5 times the conventional price.' },
  { key: 'carbon_usd_t', label: 'Carbon price', unit: 'USD/t CO₂', default: 80, min: 0, group: 'Energy', help: 'Refreshed from live data when available. EU ETS allowances averaged about 83 USD/t in 2025.' }, { key: 'carbon_cover', label: 'Share of emissions that is priced', unit: '-', default: 0.5, min: 0, max: 1, group: 'Energy', help: 'An unsourced modelling choice: it depends on which flights fall inside a pricing scheme and on free allowances. Set 1 when every tonne pays the full price (flights wholly inside a scheme), 0 where no scheme applies.' },
  { key: 'crew_usd_fh', label: 'Crew cost', unit: 'USD/FH', default: 1378, min: 0, group: 'Operating costs', help: 'Flight and cabin crew. Anchor: 1378 USD per block hour for US passenger narrow-bodies of 165,000 lb or more (FAA Economic Values, Table 4-7, year to June 2023). Block hours include taxi, so a rate per flight hour is a few per cent higher. Other classes take the figure entered in the case.' },
  { key: 'maint_usd_fh', label: 'Maintenance cost (airframe, engines, components)', unit: 'USD/FH', default: 954, min: 0, group: 'Operating costs', help: 'Total, including overhaul reserves. The shares below split it; they do not add to it. Anchor: 954 USD per block hour for US passenger narrow-bodies of 165,000 lb or more (FAA Economic Values, Table 4-7, year to June 2023). When the case gives no figure the anchor is scaled with (MTOM / 78 t)^0.6, an illustrative fixed-wing size law; rotorcraft cost several times more.' },
  { key: 'maint_labour_frac', label: 'Labour share of maintenance', unit: '-', default: 0.45, min: 0, max: 1, group: 'Operating costs' }, { key: 'eng_share', label: 'Engine / motor share of maintenance (overhaul reserve)', unit: '-', default: 0.4, min: 0, max: 1, group: 'Operating costs' },
  { key: 'rotor_share', label: 'Rotor-blade and gearbox share of maintenance', unit: '-', default: 0, min: 0, max: 1, group: 'Operating costs', help: 'Helicopter dynamic components overhaul and retirement reserve' },
  { key: 'batt_kWh', label: 'Battery pack energy', unit: 'kWh', default: 0, min: 0, group: 'Operating costs' }, { key: 'batt_usd_kWh', label: 'Battery replacement cost', unit: 'USD/kWh', default: 350, min: 0, group: 'Operating costs', help: 'An assumption for a certified aviation pack, not a market price: no aviation source was found. Automotive packs cost about 139 USD/kWh (US DOE, 2023) and 108 USD/kWh across sectors (BloombergNEF, 2025).' }, { key: 'batt_cycles', label: 'Battery cycle life', unit: 'cycles', default: 1200, min: 1, group: 'Operating costs', help: 'Equivalent full cycles to end of life. A flight that uses part of the pack counts as that part of a cycle.' },
  { key: 'landing_usd_t', label: 'Landing charge', unit: 'USD/t MTOM', default: 10, min: 0, group: 'Charges (illustrative)', help: 'Illustrative: no verified average exists. A high-cost example is New York JFK in 2026 at 8.47 USD per 1000 lb, about 18.7 USD/t.' }, { key: 'nav_usd_100km', label: 'En-route navigation unit rate', unit: 'USD per 100 km·√(MTOM/50 t)', default: 30, min: 0, group: 'Charges (illustrative)', help: 'The charge formula (distance / 100 km × √(MTOM / 50 t) × unit rate) is the EUROCONTROL route-charge formula; the unit rate itself is illustrative and varies by state.' },
  { key: 'handling_usd_cycle', label: 'Ground handling and turnaround per flight', unit: 'USD', default: 990, min: 0, group: 'Charges (illustrative)' }, { key: 'env_usd_cycle', label: 'Noise and emission charges per flight', unit: 'USD', default: 0, min: 0, group: 'Charges (illustrative)' },
  { key: 'ioc_usd_pax', label: 'Passenger service and station cost per passenger', unit: 'USD', default: 8, min: 0, group: 'Indirect costs' }, { key: 'sales_pct_rev', label: 'Sales and distribution', unit: '% of revenue', default: 4, min: 0, max: 40, group: 'Indirect costs' }, { key: 'admin_pct_doc', label: 'General and administrative', unit: '% of direct cost', default: 4, min: 0, max: 40, group: 'Indirect costs' },
  { key: 'discount', label: 'Discount rate (nominal)', unit: '-', default: 0.08, min: 0, max: 0.5, group: 'Finance' }, { key: 'inflation', label: 'General inflation', unit: '-', default: 0.025, min: -0.05, max: 0.5, group: 'Finance', help: 'Refreshed from live data when available' },
  { key: 'tax', label: 'Corporate tax rate', unit: '-', default: 0.21, min: 0, max: 0.6, group: 'Finance', help: 'US federal rate (26 U.S.C. § 11); state and other national rates differ' },
  { key: 'fx_rate', label: 'Exchange rate (local currency per USD)', unit: '-', default: 1, min: 1e-6, group: 'Finance', help: 'Headline results are repeated in the local currency. Filled from the live exchange-rate feed for the case currency when available.' }, { key: 'currency', label: 'Local currency label', type: 'text', default: 'USD', group: 'Finance' },
];
const baseOf = (keys) => BASE.filter((f) => keys.includes(f.key)), baseBut = (drop) => BASE.filter((f) => !drop.includes(f.key));
/** Cruise consumption estimate from the case when no mission result is available: { ff [kg/h], kWh per h }. */
function estBurn(c, d) {
  const V = c.mission.cruise_V_ms || c.flight.V_ms || 50, rho = isa(c.mission.cruise_alt_m || 0, c.atm.dISA_K).rho, m = c.mass.mtow_kg - 0.4 * c.mass.fuel_kg, W = m * G0, elec = c.prop.type === 'electric', jet = c.prop.type === 'turbofan' || c.prop.type === 'turbojet'; let P;
  if (c.wing.S_m2 > 0 && d.k_induced > 0) { const q = 0.5 * rho * V * V, D = q * c.wing.S_m2 * c.aero.CD0 + (d.k_induced * W * W) / (q * c.wing.S_m2); if (jet) return { ff: 1.08 * c.prop.tsfc_kg_Ns * D * 3600, kWh: 0 }; P = (D * V) / (c.prop.eta_prop || 0.8); }
  else if (c.rotor.R_m > 0) { const nR = c.meta.type === 'helicopter' ? 1 : c.prop.n_eng, T = W / nR, A = Math.PI * c.rotor.R_m ** 2, vh = Math.sqrt(T / (2 * rho * A)); let vi = vh; for (let k = 0; k < 40; k++) vi = 0.5 * (vi + (vh * vh) / Math.hypot(V, vi)); P = (nR * (1.15 * T * vi + ((rho * A * d.v_tip ** 3 * d.solidity * c.rotor.cd0) / 8) * (1 + 4.65 * (V / d.v_tip) ** 2)) + 0.5 * rho * V ** 3 * c.rotor.flat_plate_m2) / (c.meta.type === 'helicopter' ? 0.85 : 0.92); }
  else P = 0.5 * d.P_total;
  return elec ? { ff: 0, kWh: (1.1 * P) / 0.9 / 1e3 } : { ff: 1.1 * c.prop.bsfc_kg_Ws * P * 3600, kWh: 0 };
}
// FAA Economic Values, Table 4-7 (Form 41, year to June 2023): US passenger narrow-bodies of 165,000 lb (74.8 t) MTOW or more, USD per block hour
const FAA_NB = { crew: 1378, maint: 954, mtow_t: 74.8 };
/** Maintenance cost for a class without its own figure: the FAA narrow-body anchor scaled with (MTOM / 78 t)^0.6 (illustrative size law). */
const maintForClass = (mtow_t) => FAA_NB.maint * (mtow_t / 78) ** 0.6;
const CHARGE_EFF = 0.92; // grid-to-pack charging efficiency
const declared = (fields) => Object.fromEntries(fields.map((f) => [f.key, f.default]));
/** Keep only the values an analysis declares as inputs. */
const onlyInputs = (vals, an) => Object.fromEntries(an.inputs.map((f) => f.key).filter((k) => k in vals).map((k) => [k, vals[k]]));
const baseDefaults = (c, up, d) => {
  const e = c.econ, elec = c.prop.type === 'electric', heli = c.meta.type === 'helicopter', est = estBurn(c, d), mi = up.mission || {}, V = c.mission.cruise_V_ms || c.flight.V_ms || 50, small = c.mass.mtow_kg < 5700;
  // hours that carry the mission fuel: flight time when the mission publishes it (utilisation is in flight hours), else block time
  const bt = mi.flight_time_h > 0 ? mi.flight_time_h : mi.block_time_h > 0 ? mi.block_time_h : 0, cur = e.currency || 'USD', fx = e.fx_rates?.[cur];
  const b = {
    price_usd: e.price_usd || undefined, life_yr: up.fatigue?.life_fh > 0 && e.util_fh_yr > 0 ? N.clamp(Math.round(Math.min(e.life_yr, up.fatigue.life_fh / e.util_fh_yr)), 3, 50) : e.life_yr, residual_frac: e.residual_frac, insurance_pct: heli ? 1.5 : c.meta.type === 'uav' ? 3 : c.meta.type === 'evtol' ? 2 : small ? 1.2 : 0.35,
    util_fh_yr: e.util_fh_yr, cycles_yr: e.cycles_yr, block_kmh: bt && mi.mission_dist_km > 0 ? mi.mission_dist_km / bt : (c.wing.S_m2 > 0 ? 0.9 : 0.85) * V * 3.6, mtow_t: c.mass.mtow_kg / 1e3, seats: e.seats, load_factor: e.load_factor, yield_usd_pkm: e.yield_usd_pkm,
    fuel_kg_fh: elec ? 0 : bt && mi.block_fuel_kg > 0 ? mi.block_fuel_kg / bt : est.ff, fuel_usd_kg: e.fuel_usd_kg, fuel: FLUIDS[c.prop.fuel]?.LHV ? c.prop.fuel : 'Jet A-1', energy_kWh_fh: elec ? (bt && mi.mission_energy_kWh > 0 ? mi.mission_energy_kWh / bt : est.kWh) / CHARGE_EFF : 0, elec_usd_kWh: e.elec_usd_kWh, carbon_usd_t: e.carbon_usd_t,
    crew_usd_fh: e.crew_usd_fh, maint_usd_fh: e.maint_usd_fh > 0 ? e.maint_usd_fh : maintForClass(c.mass.mtow_kg / 1e3), eng_share: heli ? 0.3 : elec ? 0.15 : 0.4, rotor_share: heli ? 0.25 : 0, batt_kWh: elec ? c.systems.batt_kWh : 0, batt_cycles: BATTERIES[c.systems.batt_chem]?.cycles_80 || 1200,
    landing_usd_t: small ? 12 : 10, nav_usd_100km: small ? 0 : 30, handling_usd_cycle: heli ? 150 : c.meta.type === 'uav' ? 5 : small ? 25 : 6 * e.seats, ioc_usd_pax: e.yield_usd_pkm > 0 ? (small ? 3 : 8) : 0, sales_pct_rev: e.yield_usd_pkm > 0 ? 4 : 2,
    discount: e.discount, inflation: e.inflation, tax: e.tax, currency: cur, fx_rate: cur !== 'USD' && fx > 0 ? fx : undefined, debt_frac: small && e.crew_usd_fh === 0 ? 0 : 0.7,
  };
  // commercial operation without a passenger yield: placeholder charter rate at which the purchase just earns the discount rate, plus 5%
  // (owner-flown aircraft keep zero revenue). A cost-plus rate is the fallback when no such rate exists.
  if (!(e.yield_usd_pkm > 0) && e.crew_usd_fh > 0) {
    const i0 = { ...declared(BASE), ...Object.fromEntries(Object.entries(b).filter(([, v]) => v !== undefined)), ...finExtraDefaults(c), dep_method: 'straight-line', charter_usd_fh: 0, cargo_kg: 0 };
    const costPlus = (1.12 * opCost(i0).tocFh) / (1 - b.sales_pct_rev / 100), npv = (r) => npvOf(cashFlows({ ...i0, charter_usd_fh: r }).proj, i0.discount), hi = 20 * costPlus;
    b.charter_usd_fh = Math.ceil(npv(hi) > 0 && npv(0) < 0 ? 1.05 * N.brent(npv, 0, hi, 1e-6 * hi) : costPlus);
  }
  return b;
};
const PLACEHOLDER = 'Unless their help text names a source, prices, charges, rates and cost shares are illustrative placeholders';

/** Operating-cost model. All money in nominal USD of year 0; per flight hour unless stated. `m` holds scenario multipliers. */
function opCost(i, m = {}) {
  const safOnly = i.fuel === 'SAF (HEFA-SPK)', fl = FLUIDS[i.fuel] || FLUIDS['Jet A-1'], saf = FLUIDS['SAF (HEFA-SPK)'], x = safOnly ? 1 : i.saf_blend, util = i.util_fh_yr * (m.util ?? 1), cyc = i.cycles_yr * (m.util ?? 1), fhc = util / Math.max(cyc, 1e-9), stage = i.block_kmh * fhc;
  const fuelP = i.fuel_usd_kg * (m.fuel ?? 1), hedge = m.hedge ?? 0, pEff = hedge * i.fuel_usd_kg + (1 - hedge) * fuelP, base = safOnly ? FLUIDS['Jet A-1'] : fl;
  // The fuel price is that of the selected fuel. A blend swaps a mass share x of it for SAF at saf_price_ratio × that price, at constant
  // energy (fuel mass scales with the blend heating value). With SAF itself selected the price entered already is the SAF price.
  const lhv = (1 - x) * base.LHV + x * saf.LHV, mfuel = (i.fuel_kg_fh * (safOnly ? saf.LHV : base.LHV)) / lhv, fuel = mfuel * pEff * (safOnly ? 1 : 1 - x + x * i.saf_price_ratio), elecC = i.energy_kWh_fh * i.elec_usd_kWh * (m.fuel ?? 1);
  const co2Tail = (mfuel * ((1 - x) * (base.co2_per_kg ?? 0) + x * saf.co2_per_kg)) / 1e3, co2Life = (mfuel * ((1 - x) * (base.co2_per_kg ?? 0) + x * saf.co2_per_kg * saf.lifecycle_factor)) / 1e3, carbon = co2Life * i.carbon_usd_t * (m.carbon ?? 1) * i.carbon_cover;
  // battery wear: each flight uses the share of a full cycle that its energy is of the pack (a full cycle when the energy use is not given)
  const cycFlight = i.batt_kWh > 0 && i.energy_kWh_fh > 0 ? Math.min(1, (CHARGE_EFF * i.energy_kWh_fh * fhc) / i.batt_kWh) : 1;
  const maint = i.maint_usd_fh * (m.maint ?? 1), batt = i.batt_kWh > 0 ? (i.batt_kWh * i.batt_usd_kWh * cycFlight) / i.batt_cycles / fhc : 0;
  const airport = (i.landing_usd_t * i.mtow_t + i.handling_usd_cycle + i.env_usd_cycle) / fhc, nav = (i.nav_usd_100km * (stage / 100) * Math.sqrt(i.mtow_t / 50)) / fhc, ins = ((i.insurance_pct / 100) * i.price_usd) / util;
  const dep = (i.price_usd * (1 - i.residual_frac)) / i.life_yr / util, interest = (i.loan_rate * i.debt_frac * i.price_usd * (1 + i.residual_frac)) / 2 / util;
  const cats = { 'Fuel / energy': fuel + elecC, 'Carbon': carbon, 'Crew': i.crew_usd_fh, 'Maintenance: airframe labour': maint * (1 - i.eng_share - i.rotor_share) * i.maint_labour_frac, 'Maintenance: airframe material': maint * (1 - i.eng_share - i.rotor_share) * (1 - i.maint_labour_frac), 'Maintenance: engine reserve': maint * i.eng_share, 'Maintenance: rotor and gearbox reserve': maint * i.rotor_share, 'Battery replacement': batt, 'Airport and handling': airport, 'Navigation': nav, 'Insurance': ins, 'Depreciation': dep, 'Interest': interest };
  const doc = N.sum(Object.values(cats)), cash = doc - dep - interest, lf = N.clamp(i.load_factor * (m.demand ?? 1), 0, 1), ask = i.seats * stage * cyc, rpk = ask * lf, pax = i.seats * lf * cyc;
  const revPax = i.yield_usd_pkm * (m.yield ?? 1) * rpk, revCargo = (i.cargo_kg / 1e3) * stage * cyc * i.cargo_usd_tkm, revCh = i.charter_usd_fh * util * (m.demand ?? 1) * (m.yield ?? 1), rev = revPax + revCargo + revCh;
  const iocPax = i.ioc_usd_pax * pax, sales = (i.sales_pct_rev / 100) * rev, admin = (i.admin_pct_doc / 100) * doc * util, ioc = iocPax + sales + admin, toc = doc * util + ioc;
  return { cats, doc, cash, dep, interest, util, cyc, fhc, stage, ask, rpk, pax, lf, rev, revPax, revCargo, revCh, iocPax, sales, admin, ioc, toc, tocFh: toc / util, cashOpex: cash * util + ioc, fuelFh: fuel + elecC, carbonFh: carbon, co2Fh: co2Tail, co2LifeFh: co2Life, maintFh: maint + batt, battFh: batt };
}
/** Yearly after-tax cash-flow model in nominal money. Returns project (unlevered) and equity (levered) flows. */
function cashFlows(i, m = {}) {
  const n = Math.max(1, Math.round(i.life_yr)), o = opCost(i, m), g = 1 + i.inflation, gf = 1 + (i.fuel_escal ?? i.inflation), gc = 1 + (i.carbon_escal ?? i.inflation), dep = depSchedule(i.dep_method || 'straight-line', i.dep_life || n, i.dep_table), basis = i.price_usd * (1 - (i.dep_to_residual === false ? 0 : i.residual_frac));
  // the loan is repaid as an annuity over its own term; whatever is still owed at the end of the horizon is settled in the final year
  const term = Math.max(1, Math.round(i.loan_yr)), loan = i.debt_frac * i.price_usd, pay = annuity(loan, i.loan_rate, term), proj = [-i.price_usd], eq = [-(i.price_usd - loan)], rows = []; let bal = loan, book = i.price_usd, lossP = 0, lossE = 0;
  const fuelEnergy = o.fuelFh * o.util, carbon = o.carbonFh * o.util, other = o.cashOpex - fuelEnergy - carbon - o.sales;
  for (let t = 1; t <= n; t++) {
    const rev = o.rev * g ** t, opex = other * g ** t + fuelEnergy * gf ** t + carbon * gc ** t + o.sales * g ** t, d = (dep[t - 1] || 0) * basis, ebitda = rev - opex, ebit = ebitda - d; book -= d;
    const intr = t <= term ? bal * i.loan_rate : 0; let princ = t <= term ? Math.min(bal, pay - intr) : 0; if (t === n) princ = bal; bal -= princ;
    // tax with loss carry-forward, separately for the unlevered and levered views
    const tx = (profit, loss) => { const taxable = profit - loss; return taxable > 0 ? [i.tax * taxable, 0] : [0, -taxable]; }, [taxP, lp] = tx(ebit, lossP), [taxE, le] = tx(ebit - intr, lossE); lossP = lp; lossE = le;
    let sale = 0; if (t === n) { const sv = i.price_usd * (m.residual ?? i.residual_frac) * g ** n; sale = sv - i.tax * Math.max(0, sv - book); }
    proj.push(ebitda - taxP + sale); eq.push(ebitda - taxE - intr - princ + sale); rows.push({ t, rev, opex, ebitda, dep: d, ebit, tax: taxP, interest: intr, principal: princ, cf: proj[t], cfe: eq[t] });
  }
  return { proj, eq, rows, o, pay, loan };
}const fmt = (v) => (Math.abs(v) >= 1e9 ? `${(v / 1e9).toFixed(2)} bn` : Math.abs(v) >= 1e6 ? `${(v / 1e6).toFixed(2)} M` : Math.abs(v) >= 1e3 ? `${(v / 1e3).toFixed(1)} k` : v.toFixed(0));

// ---- analyses -------------------------------------------------------------------------------
const opcost = {
  id: 'opcost', title: 'Direct and indirect operating cost build-up', fidelity: 'analytical',
  summary: 'Builds cost per flight hour from fuel or energy, carbon, crew, maintenance and reserves, charges, insurance, depreciation and interest, adds indirect costs, and turns it into cost per flight, per seat-kilometre and per passenger, with break-even load factor and prices.',
  equations: ['Direct operating cost equations', 'Indirect operating cost equations', 'Cost per flight hour equations', 'Cost per available seat kilometre equations', 'Revenue per available seat kilometre equations', 'Yield equations', 'Load factor equations', 'Break-even analysis equations', 'Contribution margin equations', 'Depreciation equations', 'Aircraft utilisation equations'],
  inputs: BASE,
  defaults: baseDefaults,
  run(i) {
    const o = opCost(i), warnings = [], names = Object.keys(o.cats), vals = Object.values(o.cats), s = i.sales_pct_rev / 100, hasPax = i.seats > 0 && o.ask > 0;
    const fixed = o.doc * o.util + o.admin - (o.revCargo + o.revCh) * (1 - s), marginPerLf = (i.yield_usd_pkm * o.stage * (1 - s) - i.ioc_usd_pax) * i.seats * o.cyc, belf = hasPax && marginPerLf > 0 ? fixed / marginPerLf : NaN;
    const cask = hasPax ? o.toc / o.ask : NaN, rask = hasPax ? o.rev / o.ask : NaN, ticketBe = o.pax > 0 ? Math.max(0, ((o.doc * o.util + o.admin + o.iocPax) / (1 - s) - o.revCargo - o.revCh) / o.pax) : NaN, charterBe = (o.doc * o.util + o.admin + o.iocPax) / (1 - s) / o.util;
    const tkm = ((i.seats * o.lf * 0.095 + i.cargo_kg / 1e3) * o.stage * o.cyc), profit = o.rev - o.toc, payKm = o.stage * o.cyc;
    if (!hasPax) warnings.push('No seats are defined: seat-kilometre measures and the break-even load factor do not apply and are shown as not-a-number. Use the cost per flight hour and the break-even charter rate.');
    else if (!(i.yield_usd_pkm > 0)) warnings.push('Passenger yield is zero: the break-even load factor is undefined. Revenue, if any, comes from the charter rate.');
    else if (!(marginPerLf > 0)) warnings.push('Each extra passenger costs more to serve than the fare brings in: no load factor breaks even.');
    if (!(o.rev > 0)) warnings.push('No revenue is defined (private or non-commercial operation): the results are costs of ownership and use.');
    if (i.eng_share + i.rotor_share > 1) warnings.push('Engine and rotor maintenance shares add to more than 1; the airframe share is negative. Correct the shares.');
    if (i.cycles_yr > i.util_fh_yr * 12) warnings.push('Average flight time is under five minutes: check utilisation and flights per year.');
    if (i.fuel === 'Liquid hydrogen' && i.fuel_kg_fh > 0 && i.fuel_usd_kg < 2) warnings.push(`The fuel price of ${i.fuel_usd_kg} USD/kg looks like a kerosene price, but the selected fuel is liquid hydrogen: enter the delivered hydrogen price per kilogram.`);
    warnings.push(`${PLACEHOLDER}. Replace them before relying on the totals.`);
    const loc = i.fx_rate !== 1 || (i.currency && i.currency !== 'USD');
    return {
      kpis: [
        { key: 'doc_usd_fh', label: 'Direct operating cost', value: o.doc, unit: 'USD/FH' }, { key: 'cash_doc_usd_fh', label: 'Cash direct operating cost', value: o.cash, unit: 'USD/FH', note: 'Excludes depreciation and interest' },
        { key: 'toc_usd_fh', label: 'Total operating cost (direct + indirect)', value: o.tocFh, unit: 'USD/FH' }, { key: 'ioc_usd_fh', label: 'Indirect operating cost', value: o.ioc / o.util, unit: 'USD/FH' },
        { key: 'cost_per_flight_usd', label: 'Total cost per flight', value: o.toc / o.cyc, unit: 'USD' }, { key: 'cost_per_km_usd', label: 'Total cost per kilometre', value: o.toc / payKm, unit: 'USD/km' },
        { key: 'cask_usd', label: 'Cost per available seat-kilometre (CASK)', value: cask, unit: 'USD/ASK' }, { key: 'rask_usd', label: 'Revenue per available seat-kilometre (RASK)', value: rask, unit: 'USD/ASK', status: !hasPax || !(o.rev > 0) ? undefined : rask >= cask ? 'ok' : 'warn' },
        { key: 'cost_per_pax_km_usd', label: 'Cost per passenger-kilometre', value: o.rpk > 0 ? o.toc / o.rpk : NaN, unit: 'USD/pkm' }, { key: 'cost_per_tkm_usd', label: 'Cost per tonne-kilometre', value: tkm > 0 ? o.toc / tkm : NaN, unit: 'USD/tkm', note: '95 kg per passenger with baggage' },
        { key: 'breakeven_load_factor', label: 'Break-even load factor', value: belf, unit: '-', status: !Number.isFinite(belf) ? undefined : belf <= i.load_factor ? 'ok' : belf <= 1 ? 'warn' : 'bad' },
        { key: 'breakeven_ticket_usd', label: 'Break-even average fare per passenger', value: ticketBe, unit: 'USD', note: 'Fare that covers total cost after cargo and charter revenue; 0 when those already cover it' }, { key: 'breakeven_charter_usd_fh', label: 'Break-even charter rate', value: charterBe, unit: 'USD/FH', note: 'Revenue per flight hour that covers total cost' },
        { key: 'co2_cost_usd_fh', label: 'Carbon cost', value: o.carbonFh, unit: 'USD/FH' }, { key: 'co2_kg_fh', label: 'CO₂ emitted', value: 1e3 * o.co2Fh, unit: 'kg/FH' }, { key: 'co2_g_pkm', label: 'CO₂ per passenger-kilometre', value: o.rpk > 0 ? (1e6 * o.co2Fh * o.util) / o.rpk : NaN, unit: 'g/pkm' },
        { key: 'annual_opex_usd', label: 'Annual cash operating expenditure', value: o.cashOpex, unit: 'USD/yr' }, { key: 'annual_revenue_usd', label: 'Annual revenue', value: o.rev, unit: 'USD/yr' },
        { key: 'operating_profit_usd', label: 'Annual operating result', value: profit, unit: 'USD/yr', status: o.rev > 0 ? (profit >= 0 ? 'ok' : 'bad') : undefined }, { key: 'operating_margin_pct', label: 'Operating margin', value: o.rev > 0 ? (100 * profit) / o.rev : NaN, unit: '%' },
        { key: 'stage_km', label: 'Average stage length', value: o.stage, unit: 'km' }, { key: 'maint_usd_fh_total', label: 'Maintenance incl. battery replacement', value: o.maintFh, unit: 'USD/FH' },
        ...(loc ? [{ key: 'doc_local_fh', label: `Direct operating cost in ${i.currency}`, value: o.doc * i.fx_rate, unit: `${i.currency}/FH` }, { key: 'toc_local_fh', label: `Total operating cost in ${i.currency}`, value: o.tocFh * i.fx_rate, unit: `${i.currency}/FH` }] : []),
      ],
      plots: [
        { type: 'bar', title: 'Cost per flight hour by category', ylabel: 'Cost [USD/FH]', categories: ['Cost build-up'], stacked: true, series: [...[['Fuel / energy', [0]], ['Carbon', [1]], ['Crew', [2]], ['Maintenance and reserves', [3, 4, 5, 6, 7]], ['Charges', [8, 9]]].map(([n, ix]) => ({ name: n, y: [N.sum(ix.map((k) => vals[k]))] })), { name: 'Ownership and indirect', y: [vals[10] + vals[11] + vals[12] + o.ioc / o.util] }] },
        { type: 'bar', title: 'Detailed cost items', ylabel: 'Cost [USD/FH]', categories: [...names, 'Indirect: passenger service', 'Indirect: sales', 'Indirect: administration'], series: [{ name: 'USD per flight hour', y: [...vals, o.iocPax / o.util, o.sales / o.util, o.admin / o.util] }] },
        ...(hasPax && i.yield_usd_pkm > 0 ? [{ type: 'line', title: 'Annual result versus load factor', xlabel: 'Load factor [-]', ylabel: 'Operating result [USD M/yr]', series: [{ name: 'Operating result', x: N.linspace(0.3, 1, 15), y: N.linspace(0.3, 1, 15).map((lf) => { const q = opCost({ ...i, load_factor: lf }); return (q.rev - q.toc) / 1e6; }) }], annotations: [{ y: 0, label: 'Break-even' }, { x: i.load_factor, label: 'Planned' }] }] : []),
      ],
      tables: [{ title: 'Operating cost build-up', columns: ['Item', 'Class', 'USD per FH', 'USD per flight', 'USD per year', `${i.currency || 'local'} per FH`, 'Share of total [%]'], rows: [...names.map((nm, k) => [nm, k >= 11 ? 'Direct (ownership)' : 'Direct (cash)', vals[k], vals[k] * o.fhc, vals[k] * o.util, vals[k] * i.fx_rate, (100 * vals[k] * o.util) / o.toc]), ...[['Passenger service and station', o.iocPax], ['Sales and distribution', o.sales], ['General and administrative', o.admin]].map(([nm, v]) => [nm, 'Indirect', v / o.util, v / o.cyc, v, (v / o.util) * i.fx_rate, (100 * v) / o.toc]), ['Total', '', o.tocFh, o.toc / o.cyc, o.toc, o.tocFh * i.fx_rate, 100]] }],
      outputs: { stage_km: o.stage, ask_per_year: o.ask },
      warnings, models: ['Direct operating cost build-up with separate cash and ownership items', 'Indirect cost model: per passenger, share of revenue, share of direct cost', 'Weight-based landing charge and distance × √mass navigation charge', 'Straight-line depreciation to residual value and average-balance interest', 'Battery replacement from equivalent full cycles per flight', 'Cost–volume–profit break-even'],
      assumptions: ['Maintenance is one total split into airframe labour, airframe material, engine reserve and rotor/gearbox reserve: the shares divide it and never add to it', 'Battery replacement is separate from the maintenance total; a flight wears the pack by the share of a full cycle it uses', 'Depreciation and interest are ownership costs: they are in the direct operating cost but not in the cash cost or the cash-flow model', 'Stage length = block speed × utilisation ÷ flights; all money in today’s USD', 'Hours are flight hours. Rates quoted per block hour (such as the FAA crew and maintenance anchors) are used unchanged, which understates them by the taxi share of a block hour', 'The fuel price is that of the selected fuel; a blend buys its SAF share at the price ratio', PLACEHOLDER],
    };
  },
  verify() {
    const z = { ...Object.fromEntries(BASE.map((f) => [f.key, f.default])), crew_usd_fh: 0, maint_usd_fh: 0, landing_usd_t: 0, nav_usd_100km: 0, handling_usd_cycle: 0, insurance_pct: 0, carbon_usd_t: 0, price_usd: 0, ioc_usd_pax: 0, sales_pct_rev: 0, admin_pct_doc: 0 }, a = N.kv(opcost.run(z));
    const d = Object.fromEntries(BASE.map((f) => [f.key, f.default])), r = opcost.run(d), o = N.kv(r), q = opCost({ ...d, load_factor: o.breakeven_load_factor }), tot = r.tables[0].rows;
    return [
      N.check('Fuel-only cost = burn × price', a.doc_usd_fh, z.fuel_kg_fh * z.fuel_usd_kg, 1e-12, 'Definition'),
      N.check('SAF selected as the fuel: cost = burn × the price entered (no second premium)', opCost({ ...z, fuel: 'SAF (HEFA-SPK)' }).doc, z.fuel_kg_fh * z.fuel_usd_kg, 1e-12, 'The price entered is that of the selected fuel'),
      N.check('30% SAF blend: fuel cost in closed form', opCost({ ...z, saf_blend: 0.3 }).doc, (() => { const J = FLUIDS['Jet A-1'], S = FLUIDS['SAF (HEFA-SPK)']; return ((z.fuel_kg_fh * J.LHV) / (0.7 * J.LHV + 0.3 * S.LHV)) * z.fuel_usd_kg * (0.7 + 0.3 * z.saf_price_ratio); })(), 1e-12, 'Constant energy, SAF share at the price ratio'),
      N.check('Battery wear: a flight using a quarter of the pack costs a quarter of a cycle', opCost({ ...z, fuel_kg_fh: 0, batt_kWh: 100, batt_usd_kWh: 300, batt_cycles: 1000, energy_kWh_fh: 25 / CHARGE_EFF / (z.util_fh_yr / z.cycles_yr) }).battFh, (100 * 300 * 0.25) / 1000 / (z.util_fh_yr / z.cycles_yr), 1e-12, 'Equivalent full cycles'),
      N.check('Battery wear without an energy figure: one full cycle per flight', opCost({ ...z, fuel_kg_fh: 0, batt_kWh: 100, batt_usd_kWh: 300, batt_cycles: 1000 }).battFh, (100 * 300) / 1000 / (z.util_fh_yr / z.cycles_yr), 1e-12, 'Conservative fallback'),
      N.check('Carbon cost = fuel × 3.16 × price', N.kv(opcost.run({ ...z, carbon_usd_t: 100, carbon_cover: 1 })).co2_cost_usd_fh, 2.45 * 3.16 * 100, 1e-12, 'CO₂ emission factor of kerosene'),
      N.check('Straight-line depreciation per flight hour', opCost(d).dep, (d.price_usd * (1 - d.residual_frac)) / d.life_yr / d.util_fh_yr, 1e-12, 'Definition'),
      N.check('Result is zero at the break-even load factor', (q.rev - q.toc) / q.toc + 1, 1, 1e-10, 'Cost–volume–profit'),
      N.check('Cost items sum to the total (no double counting)', N.sum(tot.slice(0, -1).map((x) => x[4])), tot[tot.length - 1][4], 1e-12, 'Additivity'),
      N.check('CASK × ASK = total cost', o.cask_usd * r.outputs.ask_per_year, o.toc_usd_fh * d.util_fh_yr, 1e-12, 'Definition'),
    ];
  },
  calibration: { params: [{ key: 'maint_usd_fh', min: 0, max: 20000 }, { key: 'crew_usd_fh', min: 0, max: 10000 }, { key: 'handling_usd_cycle', min: 0, max: 50000 }], sweep: 'util_fh_yr', target: 'toc_usd_fh', note: 'Supply recorded total cost per flight hour at several utilisation levels (operator accounts) to separate fixed and variable cost.' },
  recommend(res, i) {
    const o = res.outputs, out = [], rows = res.tables[0].rows.slice(0, -1).sort((a, b) => b[4] - a[4]);
    out.push({ severity: 'info', title: `Largest cost: ${rows[0][0]} (${rows[0][6].toFixed(0)}% of total)`, detail: `Then ${rows[1][0]} (${rows[1][6].toFixed(0)}%) and ${rows[2][0]} (${rows[2][6].toFixed(0)}%).`, action: /Fuel/.test(rows[0][0]) ? 'A 1% cut in fuel burn saves about ' + fmt(0.01 * rows[0][4]) + ' USD a year and the same share of CO₂: check cruise speed and altitude (Suite 24), mass (Suite 23) and drag (Suite 1).' : `Focus cost reduction on ${rows[0][0].toLowerCase()} first.`, basis: 'Cost breakdown ranking' });
    if (Number.isFinite(o.breakeven_load_factor) && o.breakeven_load_factor > i.load_factor) out.push({ severity: o.breakeven_load_factor > 1 ? 'critical' : 'warn', title: 'The operation loses money at the planned load factor', detail: `Break-even load factor ${(100 * o.breakeven_load_factor).toFixed(0)}% against ${(100 * i.load_factor).toFixed(0)}% planned; break-even fare ${o.breakeven_ticket_usd.toFixed(0)} USD.`, action: 'Raise yield or load factor, increase utilisation to dilute ownership cost, or cut the dominant cost item.', basis: 'Break-even analysis' });
    if (o.co2_cost_usd_fh > 0.05 * o.doc_usd_fh) out.push({ severity: 'advise', title: 'Carbon is a material cost', detail: `${o.co2_cost_usd_fh.toFixed(0)} USD/FH (${(100 * o.co2_cost_usd_fh / o.doc_usd_fh).toFixed(0)}% of direct cost).`, action: 'Every fuel-saving measure now pays twice. Compare SAF and other options in the sustainability analysis.', basis: 'Emissions–carbon pricing coupling' });
    const fixedShare = (res.tables[0].rows.filter((r) => /Insurance|Depreciation|Interest/.test(r[0])).reduce((s, r) => s + r[6], 0));
    // aircraft inside the class of the sourced anchor: show how far the entered crew and maintenance rates sit from it
    if (i.mtow_t >= FAA_NB.mtow_t && i.seats >= 100 && i.fuel_kg_fh > 0) { const dc = i.crew_usd_fh / FAA_NB.crew - 1, dm = i.maint_usd_fh / FAA_NB.maint - 1; if (Math.abs(dc) > 0.2 || Math.abs(dm) > 0.2) out.push({ severity: 'advise', title: 'Crew or maintenance rate is far from the published narrow-body average', detail: `Crew ${i.crew_usd_fh.toFixed(0)} USD/FH (${dc >= 0 ? '+' : ''}${(100 * dc).toFixed(0)}%) and maintenance ${i.maint_usd_fh.toFixed(0)} USD/FH (${dm >= 0 ? '+' : ''}${(100 * dm).toFixed(0)}%) against ${FAA_NB.crew} and ${FAA_NB.maint} USD per block hour.`, action: 'Check the rates against your own accounts. The published figures are US fleet averages for the year to June 2023 and are per block hour, so per flight hour they are a few per cent higher still.', basis: 'FAA Economic Values, Table 4-7: US passenger narrow-bodies of 165,000 lb or more' }); }
    if (fixedShare > 30) out.push({ severity: 'advise', title: 'Ownership cost dominates: fly more hours', detail: `${fixedShare.toFixed(0)}% of cost is depreciation, interest and insurance.`, action: 'Each 10% more utilisation cuts cost per flight hour by about ' + (0.1 * fixedShare / 1.1).toFixed(1) + '%. Check the turnaround and availability levers in Suites 24 and 22.', basis: 'Fixed-cost dilution' });
    return out;
  },
};

const FIN_EXTRA = [
  { key: 'dep_method', label: 'Tax depreciation method', type: 'select', options: ['straight-line', 'declining balance', 'custom table'], default: 'straight-line', group: 'Tax depreciation' },
  { key: 'dep_life', label: 'Tax depreciation life', unit: 'yr', default: 12, min: 1, max: 50, step: 1, discrete: true, group: 'Tax depreciation' },
  { key: 'dep_table', label: 'Custom depreciation table (% per year)', type: 'text', default: '14.29, 24.49, 17.49, 12.49, 8.93, 8.92, 8.93, 4.46', group: 'Tax depreciation', help: 'Used with "custom table"; the default is the seven-year accelerated (MACRS-style) schedule. Values are normalised to 100%.' },
  { key: 'fuel_escal', label: 'Fuel and energy price escalation', unit: '-', default: 0.03, min: -0.1, max: 0.5, group: 'Escalation' }, { key: 'carbon_escal', label: 'Carbon price escalation', unit: '-', default: 0.05, min: -0.1, max: 0.5, group: 'Escalation' },
];
const finExtraDefaults = (c) => ({ dep_life: Math.min(12, c.econ.life_yr), fuel_escal: c.econ.inflation + 0.005, carbon_escal: c.econ.inflation + 0.025 });
const finDefaults = (c, up, d) => ({ ...baseDefaults(c, up, d), ...finExtraDefaults(c) });
const hasRevenue = (i) => opCost(i).rev > 0;
const invest = {
  id: 'invest', title: 'Investment appraisal: cash flow, NPV, IRR and payback', fidelity: 'analytical',
  summary: 'Projects revenue, operating cost, tax and financing year by year over the service life, and reduces the cash flows to net present value, internal rate of return, payback period and return on investment, for the project and for the equity holder.',
  equations: ['Net present value equation', 'Internal rate of return equation', 'Discounted cash flow equations', 'Payback period equation', 'Discounted payback period equation', 'Return on investment equation', 'Depreciation equations', 'Amortisation equations', 'Annuity equations', 'Loan repayment equations', 'Interest compounding equations', 'Capital expenditure equations', 'Operating expenditure equations', 'Profitability equations'],
  inputs: [...BASE, ...FIN_EXTRA],
  defaults: finDefaults,
  run(i) {
    const c = cashFlows(i), warnings = [], r = i.discount, npv = npvOf(c.proj, r), ir = irrOf(c.proj), npvE = npvOf(c.eq, r), irE = irrOf(c.eq), pb = paybackOf(c.proj), dpb = paybackOf(c.proj, r), n = c.proj.length - 1;
    const realR = (1 + r) / (1 + i.inflation) - 1, roi = i.price_usd > 0 ? (N.sum(c.proj.slice(1)) - i.price_usd) / i.price_usd : NaN, pi = i.price_usd > 0 ? (npv + i.price_usd) / i.price_usd : NaN, cum = N.range(n + 1, (t) => N.sum(c.proj.slice(0, t + 1))), cumD = N.range(n + 1, (t) => npvOf(c.proj.slice(0, t + 1), r));
    if (!(c.o.rev > 0)) warnings.push('No revenue is defined: every yearly cash flow is negative, so IRR and payback do not exist and NPV is the present cost of owning and operating the aircraft.');
    else if (!Number.isFinite(ir.irr)) warnings.push('Operating cash flows are too weak for a meaningful rate of return (no root above −50%): the operation does not cover its cash costs and the purchase.');
    if (ir.signChanges > 1) warnings.push(`The cash-flow series changes sign ${ir.signChanges} times: more than one IRR can exist (${ir.roots.length} found). Rely on NPV and MIRR.`);
    if (!Number.isFinite(pb)) warnings.push('The investment is not paid back within the service life.');
    if (i.loan_yr > i.life_yr) warnings.push('The loan term exceeds the appraisal horizon; the outstanding balance is repaid in the final year.');
    if (r <= i.inflation) warnings.push('The nominal discount rate is at or below inflation: the real discount rate is zero or negative.');
    warnings.push(`${PLACEHOLDER}.`);
    const yrs = N.range(n + 1), loc = i.fx_rate !== 1 || (i.currency && i.currency !== 'USD');
    return {
      kpis: [
        { key: 'npv_usd', label: 'Net present value (project, after tax)', value: npv, unit: 'USD', status: c.o.rev > 0 ? (npv >= 0 ? 'ok' : npv >= -0.02 * i.price_usd ? 'warn' : 'bad') : undefined, note: `Nominal discount rate ${(100 * r).toFixed(1)}%, real ${(100 * realR).toFixed(1)}%` },
        { key: 'irr', label: 'Internal rate of return (project)', value: ir.irr, unit: '-', status: Number.isFinite(ir.irr) ? (ir.irr >= r ? 'ok' : 'warn') : undefined }, { key: 'mirr', label: 'Modified internal rate of return', value: mirrOf(c.proj, r, r), unit: '-' },
        { key: 'payback_yr', label: 'Simple payback period', value: pb, unit: 'yr', status: Number.isFinite(pb) ? 'ok' : c.o.rev > 0 ? 'bad' : undefined }, { key: 'disc_payback_yr', label: 'Discounted payback period', value: dpb, unit: 'yr' },
        { key: 'roi', label: 'Return on investment over the life', value: roi, unit: '-', note: '(Σ net cash − investment) ÷ investment, undiscounted' }, { key: 'profitability_index', label: 'Profitability index', value: pi, unit: '-', note: 'PV of inflows ÷ investment' },
        { key: 'npv_equity_usd', label: 'Net present value to equity (with loan)', value: npvE, unit: 'USD' }, { key: 'irr_equity', label: 'Internal rate of return on equity', value: irE.irr, unit: '-' },
        { key: 'loan_payment_usd_yr', label: 'Annual loan payment', value: c.pay, unit: 'USD/yr' }, { key: 'equity_usd', label: 'Equity invested', value: i.price_usd - c.loan, unit: 'USD' },
        { key: 'ebitda_y1_usd', label: 'Year-1 EBITDA', value: c.rows[0].ebitda, unit: 'USD' }, { key: 'eaa_usd_yr', label: 'Equivalent annual value', value: npv / pvFactor(r, n), unit: 'USD/yr', note: 'NPV spread as a level annual amount' },
        { key: 'real_discount_rate', label: 'Real discount rate (Fisher)', value: realR, unit: '-' },
        ...(loc ? [{ key: 'npv_local', label: `Net present value in ${i.currency}`, value: npv * i.fx_rate, unit: i.currency }] : []),
      ],
      plots: [
        { type: 'line', title: 'Cumulative cash flow (project)', xlabel: 'Year', ylabel: 'Cumulative cash flow [USD M]', series: [{ name: 'Undiscounted', x: yrs, y: cum.map((v) => v / 1e6), style: 'line+points' }, { name: 'Discounted', x: yrs, y: cumD.map((v) => v / 1e6), style: 'line+points' }], annotations: [{ y: 0, label: 'Payback' }] },
        { type: 'bar', title: 'Yearly cash flow', ylabel: 'Cash flow [USD M]', categories: yrs.map(String), series: [{ name: 'Project (unlevered)', y: c.proj.map((v) => v / 1e6) }, { name: 'Equity (levered)', y: c.eq.map((v) => v / 1e6) }] },
        { type: 'line', title: 'Net present value versus discount rate', xlabel: 'Discount rate [-]', ylabel: 'NPV [USD M]', series: [{ name: 'NPV', x: N.linspace(0, 0.25, 26), y: N.linspace(0, 0.25, 26).map((q) => npvOf(c.proj, q) / 1e6) }], annotations: [{ y: 0, label: 'IRR' }, { x: r, label: 'Discount rate' }] },
      ],
      tables: [{ title: 'Cash-flow statement (nominal USD)', columns: ['Year', 'Revenue', 'Cash operating cost', 'EBITDA', 'Tax depreciation', 'Tax (unlevered)', 'Project cash flow', 'Interest', 'Principal', 'Equity cash flow'], rows: [[0, 0, 0, 0, 0, 0, c.proj[0], 0, 0, c.eq[0]], ...c.rows.map((q) => [q.t, q.rev, q.opex, q.ebitda, q.dep, q.tax, q.cf, q.interest, q.principal, q.cfe])] }],
      warnings, models: ['Yearly discounted cash-flow model in nominal money', `Tax depreciation: ${i.dep_method}`, 'Annuity loan amortisation', 'Tax on operating profit with loss carry-forward', 'IRR by bracketing and Brent’s method with multiple-root detection', 'MIRR, simple and discounted payback, profitability index'],
      assumptions: ['Project NPV uses unlevered after-tax cash flows: interest is not deducted, because the discount rate already carries the cost of capital', 'Equity view deducts interest and principal and uses the interest tax shield', 'Revenue and non-fuel costs escalate with general inflation; fuel and carbon with their own rates', 'Residual value is the stated fraction of price, escalated with inflation and taxed on the gain over book value', 'Cash flows fall at year end; the purchase is at year 0', PLACEHOLDER],
    };
  },
  verify() {
    const cf = [-1000, ...new Array(10).fill(200)], ir = irrOf([-100, 110]), two = irrOf([-100, 230, -132]);
    const d = { ...Object.fromEntries([...BASE, ...FIN_EXTRA].map((f) => [f.key, f.default])) }, sch = depSchedule('declining balance', 8), c = cashFlows(d);
    return [
      N.check('NPV of a level annuity: −I + A·(1 − (1+r)⁻ⁿ)/r', npvOf(cf, 0.08), -1000 + (200 * (1 - 1.08 ** -10)) / 0.08, 1e-12, 'Annuity present-value formula'),
      N.check('IRR of (−100, +110) is 10%', ir.irr, 0.1, 1e-9, 'Definition'),
      N.check('Two IRRs of (−100, 230, −132): 10% and 20%', two.roots.length === 2 ? two.roots[0] + two.roots[1] : NaN, 0.3, 1e-8, 'Quadratic with roots 1.1 and 1.2'),
      N.check('Annuity payment P·r/(1 − (1+r)⁻ⁿ)', annuity(1e6, 0.06, 10), (1e6 * 0.06) / (1 - 1.06 ** -10), 1e-12, 'Loan formula'),
      N.check('Loan is fully repaid: principal sums to the amount borrowed', N.sum(c.rows.map((q) => q.principal)), c.loan, 1e-9, 'Amortisation identity'),
      N.check('Loan longer than the horizon: annual payment is the annuity over the loan term', cashFlows({ ...d, life_yr: 5, loan_yr: 10 }).pay, annuity(d.debt_frac * d.price_usd, d.loan_rate, 10), 1e-12, 'Loan formula'),
      N.check('Loan longer than the horizon: the balance is settled in the final year', N.sum(cashFlows({ ...d, life_yr: 5, loan_yr: 10 }).rows.map((q) => q.principal)), d.debt_frac * d.price_usd, 1e-9, 'Balloon = balance of a 10-year annuity after 5 payments'),
      N.check('Balloon equals the annuity balance P·(1+r)⁴ − A·((1+r)⁴ − 1)/r', cashFlows({ ...d, life_yr: 5, loan_yr: 10 }).rows[4].principal, (() => { const P = d.debt_frac * d.price_usd, r = d.loan_rate, A = annuity(P, r, 10), k = (1 + r) ** 4; return P * k - (A * (k - 1)) / r; })(), 1e-9, 'Balance after four payments, all repaid in year 5'),
      N.check('Declining-balance schedule sums to 100%', N.sum(sch), 1, 1e-12, 'Depreciation identity'),
      N.check('Payback of (−1000, 400, 400, 400) is 2.5 years', paybackOf([-1000, 400, 400, 400]), 2.5, 1e-12, 'Linear interpolation within the year'),
      N.check('MIRR of (−100, 0, 121) at 10% is 10%', mirrOf([-100, 0, 121], 0.1, 0.1), 0.1, 1e-12, 'Definition'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!(o.annual_revenue_usd > 0) && !(o.ebitda_y1_usd > 0) && !Number.isFinite(o.irr)) { out.push({ severity: 'info', title: 'This is a cost-of-ownership case', detail: `Present cost of owning and operating over ${Math.round(i.life_yr)} years: ${fmt(-o.npv_usd)} USD (${fmt(-o.eaa_usd_yr)} USD per year equivalent).`, action: 'Enter a charter rate or yield to appraise a commercial operation, or compare this annual figure with renting or fractional ownership.', basis: 'Net present cost' }); return out; }
    // an NPV within 2% of the price either way is break-even at the accuracy of the inputs
    const thin = Math.abs(o.npv_usd) <= 0.02 * i.price_usd;
    if (o.npv_usd < 0) out.push({ severity: thin ? 'advise' : 'warn', title: thin ? 'The investment about earns its cost of capital' : 'The investment does not earn its cost of capital', detail: `NPV ${fmt(o.npv_usd)} USD; IRR ${Number.isFinite(o.irr) ? (100 * o.irr).toFixed(1) + '%' : 'undefined'} against a ${(100 * i.discount).toFixed(1)}% discount rate.`, action: 'Test the levers in the risk analysis tornado: yield and load factor, utilisation, purchase price, and fuel burn. A lower-priced or leased aircraft may clear the hurdle (see lease versus buy).', basis: 'NPV decision rule' });
    else out.push({ severity: 'info', title: 'The investment clears the hurdle rate', detail: `NPV ${fmt(o.npv_usd)} USD, IRR ${(100 * o.irr).toFixed(1)}%, payback ${Number.isFinite(o.payback_yr) ? o.payback_yr.toFixed(1) : '—'} years.`, action: 'Stress it with the risk analysis before committing; check that the result survives a fuel or carbon price shock.', basis: 'NPV decision rule' });
    if (Number.isFinite(o.irr_equity) && Number.isFinite(o.irr) && o.irr < i.loan_rate) out.push({ severity: 'advise', title: 'Debt works against you here', detail: `Project IRR ${(100 * o.irr).toFixed(1)}% is below the ${(100 * i.loan_rate).toFixed(1)}% loan rate, so leverage lowers the equity return.`, action: 'Reduce the debt share or renegotiate the rate.', basis: 'Financial leverage' });
    return out;
  },
};

const lcc = {
  id: 'lcc', title: 'Life-cycle cost, development and production cost with learning', fidelity: 'analytical',
  summary: 'Adds development, production, operation and disposal into a life-cycle cost per aircraft, with a learning curve for unit cost against quantity and the programme break-even quantity.',
  equations: ['Life-cycle cost equations', 'Total cost of ownership equations', 'Learning curve equations', 'Capital expenditure equations', 'Operating expenditure equations', 'Marginal cost equations', 'Average cost equations', 'Discounted cash flow equations'],
  inputs: [...baseBut(['loan_yr', 'tax']),
    { key: 'oew_kg', label: 'Operating empty mass', unit: 'kg', default: 42600, min: 0.1, group: 'Programme (illustrative)' },
    { key: 'dev_a', label: 'Development cost coefficient a', unit: 'USD/kg^b', default: 1e5, min: 0, group: 'Programme (illustrative)', help: 'RDT&E = a · (empty mass)^b. Illustrative power law; calibrate to programmes of the same class.' }, { key: 'dev_b', label: 'Development cost exponent b', unit: '-', default: 1.1, min: 0.3, max: 2, group: 'Programme (illustrative)' },
    { key: 'qty', label: 'Production quantity', unit: 'aircraft', default: 1500, min: 1, max: 100000, step: 1, discrete: true, group: 'Programme (illustrative)' },
    { key: 'mfg_share', label: 'Average production cost ÷ price at that quantity', unit: '-', default: 0.7, min: 0.05, max: 2, group: 'Programme (illustrative)', help: 'Sets the first-unit cost through the learning curve' },
    { key: 'learning', label: 'Learning rate', unit: '-', default: 0.85, min: 0.6, max: 1, group: 'Programme (illustrative)', help: 'Cost ratio for each doubling of quantity. 0.85 is the aerospace value of the NASA Cost Estimating Handbook; RAND airframe data put manufacturing labour near 0.77 and materials near 0.86.' },
    { key: 'learn_model', label: 'Learning-curve model', type: 'select', options: ['Wright (cumulative average)', 'Crawford (unit)'], default: 'Wright (cumulative average)', group: 'Programme (illustrative)' },
    { key: 'disposal_frac', label: 'Disposal cost ÷ price', unit: '-', default: 0.01, min: 0, max: 0.3, group: 'Capital' }, { key: 'maint_growth', label: 'Real growth of maintenance cost with age', unit: '1/yr', default: 0.015, min: 0, max: 0.15, group: 'Operating costs' },
  ],
  defaults: (c, up, d) => onlyInputs({ ...baseDefaults(c, up, d), oew_kg: c.mass.oew_kg, qty: c.mass.mtow_kg > 40000 ? 1500 : c.mass.mtow_kg > 5700 ? 600 : 2000, learning: c.meta.type === 'uav' ? 0.9 : 0.85 }, lcc),
  run(i) {
    const warnings = [], o = opCost(i), n = Math.round(i.life_yr), b = learnB(i.learning), Q = Math.round(i.qty), dev = i.dev_a * i.oew_kg ** i.dev_b, rr = (1 + i.discount) / (1 + i.inflation) - 1;
    // first-unit cost from the average production cost at the planned quantity
    const avgTarget = i.mfg_share * i.price_usd, T1 = avgTarget / (cumCost(1, Q, i.learning, i.learn_model) / Q), total = (q) => cumCost(T1, q, i.learning, i.learn_model), unitAt = (q) => (i.learn_model.startsWith('Wright') ? T1 * (q ** (1 + b) - (q - 1) ** (1 + b)) : T1 * q ** b);
    const qs = [...new Set(N.logspace(1, Math.max(4 * Q, 10), 40).map(Math.round))], be = N.findRoot((q) => i.price_usd * q - dev - total(q), 1, 100 * Q, 400), unitFull = total(Q) / Q + dev / Q;
    // operator life-cycle cost per aircraft in constant (real) money, with real maintenance growth
    const cashNoMaint = o.cashOpex - o.maintFh * o.util, maintY = (t) => o.maintFh * o.util * (1 + i.maint_growth) ** (t - 1), yrs = N.range(n, (t) => t + 1), ops = N.sum(yrs.map((t) => cashNoMaint + maintY(t))), opsPv = N.sum(yrs.map((t) => (cashNoMaint + maintY(t)) / (1 + rr) ** t));
    const resid = i.residual_frac * i.price_usd, disp = i.disposal_frac * i.price_usd, lccU = i.price_usd + ops + disp - resid, lccPv = i.price_usd + opsPv + (disp - resid) / (1 + rr) ** n, fhLife = o.util * n;
    const parts = { Acquisition: i.price_usd, 'Fuel / energy': o.fuelFh * o.util * n, Carbon: o.carbonFh * o.util * n, Crew: i.crew_usd_fh * o.util * n, Maintenance: N.sum(yrs.map(maintY)), 'Charges and insurance': (o.cats['Airport and handling'] + o.cats.Navigation + o.cats.Insurance) * o.util * n, Indirect: o.ioc * n, 'Disposal less residual': disp - resid };
    if (!(i.price_usd > 0)) warnings.push('Aircraft price is zero: acquisition and learning-curve results are meaningless.');
    if (unitFull > i.price_usd) warnings.push(`At ${Q} aircraft the full unit cost (production plus amortised development) is above the price: the programme does not recover its development cost at this quantity.`);
    if (!Number.isFinite(be)) warnings.push('No break-even quantity exists within 100× the planned production: price does not cover the marginal unit cost.');
    warnings.push('The development-cost power law and the production cost share are illustrative; sourced cost-estimating relationships and labour rates are needed for a real estimate.');
    return {
      kpis: [
        { key: 'lcc_usd', label: 'Life-cycle cost per aircraft (constant money, undiscounted)', value: lccU, unit: 'USD' }, { key: 'lcc_pv_usd', label: 'Life-cycle cost per aircraft (present value)', value: lccPv, unit: 'USD', note: `Real discount rate ${(100 * rr).toFixed(1)}%` },
        { key: 'lcc_per_fh_usd', label: 'Life-cycle cost per flight hour', value: lccU / fhLife, unit: 'USD/FH' }, { key: 'acq_share_pct', label: 'Acquisition share of life-cycle cost', value: (100 * (i.price_usd - resid)) / lccU, unit: '%' },
        { key: 'ops_cost_life_usd', label: 'Operating expenditure over the life', value: ops, unit: 'USD' }, { key: 'capex_usd', label: 'Capital expenditure per aircraft', value: i.price_usd, unit: 'USD' },
        { key: 'dev_cost_usd', label: 'Development cost (RDT&E)', value: dev, unit: 'USD' }, { key: 'first_unit_cost_usd', label: 'First-unit production cost', value: T1, unit: 'USD' },
        { key: 'avg_unit_cost_usd', label: `Average production cost over ${Q} aircraft`, value: total(Q) / Q, unit: 'USD' }, { key: 'last_unit_cost_usd', label: `Cost of aircraft number ${Q}`, value: unitAt(Q), unit: 'USD', note: 'Marginal cost' },
        { key: 'unit_cost_full_usd', label: 'Unit cost incl. amortised development', value: unitFull, unit: 'USD', status: unitFull <= i.price_usd ? 'ok' : 'warn' },
        { key: 'breakeven_qty', label: 'Programme break-even quantity', value: be, unit: 'aircraft', status: Number.isFinite(be) && be <= Q ? 'ok' : 'warn' }, { key: 'programme_margin_pct', label: 'Programme margin at planned quantity', value: i.price_usd > 0 ? 100 * (1 - unitFull / i.price_usd) : NaN, unit: '%' },
        ...(i.fx_rate !== 1 || (i.currency && i.currency !== 'USD') ? [{ key: 'lcc_local', label: `Life-cycle cost per aircraft in ${i.currency}`, value: lccU * i.fx_rate, unit: i.currency }] : []),
      ],
      plots: [
        { type: 'bar', title: 'Life-cycle cost breakdown per aircraft', ylabel: 'Cost [USD M]', categories: ['Life-cycle cost'], stacked: true, series: Object.entries(parts).filter(([, v]) => v > 0).map(([nm, v]) => ({ name: nm, y: [v / 1e6] })).slice(0, 6) },
        { type: 'line', title: 'Unit cost versus production quantity', xlabel: 'Cumulative quantity [aircraft]', ylabel: 'Cost per aircraft [USD M]', xlog: true, ylog: true, series: [{ name: 'Average production cost', x: qs, y: qs.map((q) => total(q) / q / 1e6) }, { name: 'Marginal (unit) cost', x: qs, y: qs.map((q) => unitAt(q) / 1e6) }, { name: 'Average incl. development', x: qs, y: qs.map((q) => (total(q) + dev) / q / 1e6) }, { name: 'Price', x: [qs[0], qs[qs.length - 1]], y: [i.price_usd / 1e6, i.price_usd / 1e6], style: 'dash' }], annotations: [{ x: Q, label: 'Planned' }] },
        { type: 'line', title: 'Annual operating cost with ageing', xlabel: 'Year', ylabel: 'Cost [USD M/yr, constant money]', series: [{ name: 'Maintenance', x: yrs, y: yrs.map((t) => maintY(t) / 1e6) }, { name: 'Other cash operating cost', x: yrs, y: yrs.map(() => cashNoMaint / 1e6) }] },
      ],
      tables: [{ title: 'Life-cycle cost per aircraft (constant USD)', columns: ['Element', 'USD', 'Share [%]', 'USD per FH'], rows: [...Object.entries(parts).map(([nm, v]) => [nm, v, (100 * v) / lccU, v / fhLife]), ['Total', lccU, 100, lccU / fhLife]] }],
      warnings, models: ['Life-cycle cost = acquisition + operation + disposal − residual', `Learning curve: ${i.learn_model}`, 'Power-law development cost-estimating relationship (user-parameterised)', 'Programme break-even quantity'],
      assumptions: ['Life-cycle cost is in constant money; the present value uses the real discount rate', 'Ownership financing cost is not added to life-cycle cost (it is in the discount rate)', 'Maintenance grows at a constant real rate with age; other operating costs are constant in real terms', 'First-unit cost is back-calculated from the stated average production cost at the planned quantity', PLACEHOLDER],
    };
  },
  verify() {
    const d = Object.fromEntries(lcc.inputs.map((f) => [f.key, f.default])), o = N.kv(lcc.run(d)), r = lcc.run(d).tables[0].rows;
    return [
      N.check('Wright: doubling quantity multiplies the cumulative average by the learning rate', cumCost(1, 200, 0.85, 'Wright') / 200 / (cumCost(1, 100, 0.85, 'Wright') / 100), 0.85, 1e-12, 'Wright (1936)'),
      N.check('Crawford: unit 8 costs T1·LR³', cumCost(1, 8, 0.8, 'Crawford') - cumCost(1, 7, 0.8, 'Crawford'), 0.8 ** 3, 1e-12, 'Unit learning curve'),
      N.check('Average production cost matches the stated share of price', o.avg_unit_cost_usd, 0.7 * 105e6, 1e-10, 'Calibration identity'),
      N.check('Life-cycle elements sum to the total', N.sum(r.slice(0, -1).map((x) => x[1])), o.lcc_usd, 1e-10, 'Additivity'),
      N.check('Crawford large-quantity integral matches the direct sum', (2000.5 ** (1 + learnB(0.85)) - 0.5 ** (1 + learnB(0.85))) / (1 + learnB(0.85)), cumCost(1, 2000, 0.85, 'Crawford'), 2e-3, 'Midpoint-rule approximation'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], rows = res.tables[0].rows.slice(0, -1).sort((a, b) => b[1] - a[1]);
    out.push({ severity: 'info', title: `Life-cycle cost is ${(o.lcc_usd / Math.max(i.price_usd, 1)).toFixed(1)}× the purchase price`, detail: `Largest elements: ${rows[0][0]} (${rows[0][2].toFixed(0)}%), ${rows[1][0]} (${rows[1][2].toFixed(0)}%).`, action: `Judge design and procurement choices on life-cycle cost, not price: a dearer but more efficient configuration pays back when ${rows[0][0].toLowerCase()} dominates. Lower fuel burn also cuts lifetime CO₂.`, basis: 'Life-cycle cost breakdown' });
    if (o.unit_cost_full_usd > i.price_usd) out.push({ severity: 'warn', title: 'Programme does not break even at the planned quantity', detail: `Full unit cost ${fmt(o.unit_cost_full_usd)} USD against a price of ${fmt(i.price_usd)} USD; break-even at ${Number.isFinite(o.breakeven_qty) ? Math.round(o.breakeven_qty) : 'no'} aircraft.`, action: 'Raise the planned quantity (family variants, shared components), cut development scope, or improve the learning rate through design for manufacture.', basis: 'Cost–volume–profit with learning' });
    return out;
  },
};

const finance = {
  id: 'finance', title: 'Lease versus buy, financing and replacement timing', fidelity: 'analytical',
  summary: 'Compares the after-tax present cost of buying with cash, buying with a loan, a finance lease and an operating lease over the same period, and finds the age at which replacing the aircraft minimises the equivalent annual cost.',
  equations: ['Aircraft financing equations', 'Annuity equations', 'Loan repayment equations', 'Net present value equation', 'Depreciation equations', 'Amortisation equations', 'Total cost of ownership equations'],
  inputs: [...baseOf(['price_usd', 'life_yr', 'residual_frac', 'debt_frac', 'loan_rate', 'loan_yr', 'util_fh_yr', 'maint_usd_fh', 'discount', 'inflation', 'tax', 'fx_rate', 'currency']),
    { key: 'horizon_yr', label: 'Comparison period', unit: 'yr', default: 10, min: 1, max: 40, step: 1, discrete: true, group: 'Lease' },
    { key: 'lease_factor_pct', label: 'Operating lease rate factor', unit: '%/month of price', default: 0.8, min: 0.1, max: 3, group: 'Lease', help: 'Monthly rent as a percentage of aircraft value; illustrative' },
    { key: 'lease_reserve_usd_fh', label: 'Lessor maintenance reserve', unit: 'USD/FH', default: 0, min: 0, group: 'Lease', help: 'Paid on top of your own maintenance; enter only the part not refunded' },
    { key: 'finlease_rate', label: 'Finance lease implicit rate', unit: '-', default: 0.065, min: 0, max: 0.4, group: 'Lease' },
    { key: 'maint_growth', label: 'Real growth of maintenance cost with age', unit: '1/yr', default: 0.04, min: 0, max: 0.3, group: 'Replacement' },
    { key: 'overhaul_usd', label: 'Heavy check / overhaul cost', unit: 'USD', default: 0, min: 0, group: 'Replacement' }, { key: 'overhaul_yr', label: 'Heavy check interval', unit: 'yr', default: 6, min: 1, max: 20, step: 1, discrete: true, group: 'Replacement' },
  ],
  defaults: (c, up, d) => { const b = baseDefaults(c, up, d); return { ...Object.fromEntries(['price_usd', 'life_yr', 'residual_frac', 'debt_frac', 'util_fh_yr', 'maint_usd_fh', 'discount', 'inflation', 'tax', 'fx_rate', 'currency'].map((k) => [k, b[k]])), debt_frac: Math.max(b.debt_frac, 0.5), horizon_yr: Math.min(10, c.econ.life_yr), overhaul_usd: 0.04 * (c.econ.price_usd || 0), lease_reserve_usd_fh: 0 }; },
  run(i) {
    const warnings = [], H = Math.min(Math.round(i.horizon_yr), Math.round(i.life_yr)), r = i.discount, P = i.price_usd, g = 1 + i.inflation, mv = (t) => P * i.residual_frac ** (t / i.life_yr) * g ** t, dep = (P * (1 - i.residual_frac)) / i.life_yr, book = (t) => P - dep * t;
    const saleAfterTax = mv(H) - i.tax * (mv(H) - book(H)), shield = N.sum(N.range(H, (t) => (i.tax * dep) / (1 + r) ** (t + 1)));
    // present cost (positive = cost) of each way of having the aircraft for H years, after tax
    const cash = P - shield - saleAfterTax / (1 + r) ** H;
    const L = i.debt_frac * P, pay = annuity(L, i.loan_rate, Math.round(i.loan_yr)); let bal = L, pvLoan = P - L;
    for (let t = 1; t <= H; t++) { const intr = t <= Math.round(i.loan_yr) ? bal * i.loan_rate : 0, pr = t <= Math.round(i.loan_yr) ? pay - intr : 0; bal -= pr; pvLoan += (intr * (1 - i.tax) + pr) / (1 + r) ** t; }
    const loan = pvLoan + bal / (1 + r) ** H - shield - saleAfterTax / (1 + r) ** H;
    const rent = 12 * (i.lease_factor_pct / 100) * P, opl = N.sum(N.range(H, (t) => ((rent * g ** t + i.lease_reserve_usd_fh * i.util_fh_yr * g ** t) * (1 - i.tax)) / (1 + r) ** t)); // rent in advance, indexed
    const fpay = annuity(P, i.finlease_rate, H); let fb = P, fin = 0; for (let t = 1; t <= H; t++) { const intr = fb * i.finlease_rate; fb -= fpay - intr; fin += (fpay - i.tax * intr) / (1 + r) ** t; } fin += -shield - saleAfterTax / (1 + r) ** H;
    const opts = [['Buy with cash', cash], [`Buy with ${(100 * i.debt_frac).toFixed(0)}% loan`, loan], ['Finance lease', fin], ['Operating lease', opl]], best = N.argmin(opts.map((q) => q[1])), rentPv = N.sum(N.range(H, (t) => (rent * g ** t * (1 - i.tax)) / (1 + r) ** t)), beFactor = rentPv > 0 ? (i.lease_factor_pct * (Math.min(cash, loan, fin) - (opl - rentPv))) / rentPv : NaN; // rent scales with the factor, the reserve does not
    // economic life: minimum equivalent annual cost in real terms
    const rr = (1 + r) / g - 1, m0 = i.maint_usd_fh * i.util_fh_yr, ages = N.range(Math.min(40, Math.max(Math.round(1.6 * i.life_yr), 8)), (k) => k + 1), eac = ages.map((n) => { let pv = P - (P * i.residual_frac ** (n / i.life_yr)) / (1 + rr) ** n; for (let t = 1; t <= n; t++) pv += (m0 * (1 + i.maint_growth) ** (t - 1) + (i.overhaul_usd > 0 && t % Math.round(i.overhaul_yr) === 0 && t < n ? i.overhaul_usd : 0)) / (1 + rr) ** t; return pv / pvFactor(rr, n); }), kb = N.argmin(eac);
    if (i.loan_rate * (1 - i.tax) > r) warnings.push('The after-tax loan rate is above the discount rate, so borrowing adds present cost.');
    if (kb === ages.length - 1) warnings.push('Equivalent annual cost is still falling at the oldest age examined: with this maintenance growth there is no economic reason to replace within the range.');
    warnings.push('Lease rate factor, implicit rate and market-value curve are illustrative; use lessor quotes and appraiser values.');
    return {
      kpis: [
        { key: 'pv_cost_cash_usd', label: 'Present cost: buy with cash', value: cash, unit: 'USD' }, { key: 'pv_cost_loan_usd', label: 'Present cost: buy with loan', value: loan, unit: 'USD' }, { key: 'pv_cost_finlease_usd', label: 'Present cost: finance lease', value: fin, unit: 'USD' }, { key: 'pv_cost_oplease_usd', label: 'Present cost: operating lease', value: opl, unit: 'USD' },
        { key: 'best_option', label: 'Cheapest option (1 cash, 2 loan, 3 finance lease, 4 operating lease)', value: best + 1, unit: '' }, { key: 'lease_vs_buy_usd', label: 'Operating lease minus best purchase route', value: opl - Math.min(cash, loan, fin), unit: 'USD', note: 'Positive = leasing costs more' },
        { key: 'breakeven_lease_factor_pct', label: 'Break-even lease rate factor', value: beFactor, unit: '%/month', note: 'Below this, the operating lease beats the cheapest purchase route' }, { key: 'annual_rent_usd', label: 'First-year operating lease rent', value: rent, unit: 'USD/yr' },
        { key: 'loan_payment_usd_yr', label: 'Annual loan payment', value: pay, unit: 'USD/yr' }, { key: 'market_value_end_usd', label: `Market value after ${H} years`, value: mv(H), unit: 'USD' },
        { key: 'economic_life_yr', label: 'Economic life (minimum equivalent annual cost)', value: ages[kb], unit: 'yr' }, { key: 'eac_min_usd_yr', label: 'Minimum equivalent annual cost', value: eac[kb], unit: 'USD/yr', note: 'Capital plus maintenance, constant money' },
        { key: 'eac_at_life_usd_yr', label: 'Equivalent annual cost at the planned life', value: eac[Math.min(ages.length, Math.round(i.life_yr)) - 1], unit: 'USD/yr' },
      ],
      plots: [
        { type: 'bar', title: `After-tax present cost over ${H} years`, ylabel: 'Present cost [USD M]', categories: opts.map((q) => q[0]), series: [{ name: 'Present cost', y: opts.map((q) => q[1] / 1e6) }] },
        { type: 'line', title: 'Equivalent annual cost versus replacement age', xlabel: 'Age at replacement [yr]', ylabel: 'Equivalent annual cost [USD M/yr]', series: [{ name: 'Capital + maintenance', x: ages, y: eac.map((v) => v / 1e6) }], annotations: [{ x: ages[kb], label: 'Economic life' }, { x: i.life_yr, label: 'Planned life' }] },
        { type: 'line', title: 'Market value and book value', xlabel: 'Age [yr]', ylabel: 'Value [USD M]', series: [{ name: 'Market value (nominal)', x: N.range(Math.round(i.life_yr) + 1), y: N.range(Math.round(i.life_yr) + 1, (t) => mv(t) / 1e6) }, { name: 'Book value', x: N.range(Math.round(i.life_yr) + 1), y: N.range(Math.round(i.life_yr) + 1, (t) => book(t) / 1e6), style: 'dash' }] },
      ],
      outputs: { loan_terms_in_horizon: Math.min(Math.round(i.loan_yr), H) },
      tables: [{ title: 'Financing comparison', columns: ['Option', 'Present cost [USD]', `Present cost [${i.currency || 'local'}]`, 'Versus cheapest [USD]'], rows: opts.map((q) => [q[0], q[1], q[1] * i.fx_rate, q[1] - opts[best][1]]) }],
      warnings, models: ['After-tax present cost of cash purchase, annuity loan, finance lease and operating lease', 'Exponential market-value curve to the residual value', 'Economic-life model: equivalent annual cost with rising maintenance and falling resale value'],
      assumptions: ['All options give the same aircraft and the same operating cost, so only financing, tax shields and end value differ', 'Operating-lease rent is paid in advance and indexed to inflation; purchase options sell at market value at the end of the period', 'Straight-line tax depreciation; tax shields are usable when they arise', 'The economic-life calculation is in constant money at the real discount rate'],
    };
  },
  verify() {
    const b = { ...Object.fromEntries(finance.inputs.map((f) => [f.key, f.default])), tax: 0, inflation: 0, lease_reserve_usd_fh: 0, residual_frac: 1e-9 }, o = N.kv(finance.run({ ...b, horizon_yr: 10, maint_usd_fh: 0, maint_growth: 0, overhaul_usd: 0 }));
    const rent = 12 * 0.008 * b.price_usd, due = rent * pvFactor(b.discount, 10) * (1 + b.discount), f0 = N.kv(finance.run({ ...b, finlease_rate: b.discount, life_yr: 10, horizon_yr: 10 }));
    return [
      N.check('Operating lease present cost = annuity-due value of the rent', o.pv_cost_oplease_usd, due, 1e-12, 'Annuity due: PV = A·(1 − (1+r)⁻ⁿ)/r·(1+r)'),
      N.check('Finance lease at the discount rate costs the same as paying cash', f0.pv_cost_finlease_usd, f0.pv_cost_cash_usd, 1e-9, 'Present value of an annuity at its own rate equals the principal'),
      N.check('At the break-even lease factor the lease costs the same as the cheapest purchase route', (() => { const q = { ...b, lease_reserve_usd_fh: 150, tax: 0.21, inflation: 0.02, residual_frac: 0.15 }, k = N.kv(finance.run(q)), e = N.kv(finance.run({ ...q, lease_factor_pct: k.breakeven_lease_factor_pct })); return e.lease_vs_buy_usd / e.pv_cost_oplease_usd + 1; })(), 1, 1e-9, 'Rent is linear in the factor; the maintenance reserve is not scaled'),
      N.check('Equivalent annual cost with no maintenance = P × capital recovery factor', o.eac_at_life_usd_yr, b.price_usd / pvFactor(b.discount, 25), 1e-6, 'Capital recovery factor'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], t = res.tables[0].rows, best = t[o.best_option - 1];
    out.push({ severity: 'advise', title: `Cheapest route over ${Math.round(i.horizon_yr)} years: ${best[0]}`, detail: t.map((r) => `${r[0]} ${fmt(r[1])}`).join('; ') + ' USD present cost.', action: o.best_option === 4 ? 'Lease, and keep the flexibility to return the aircraft; check return conditions and reserves in the contract.' : `Buy. The operating lease only wins below a rate factor of ${Number(o.breakeven_lease_factor_pct).toFixed(2)}% per month, unless flexibility or balance-sheet limits matter more than cost.`, basis: 'After-tax present cost comparison' });
    if (o.economic_life_yr < i.life_yr - 2) out.push({ severity: 'advise', title: `Replace earlier: economic life is about ${o.economic_life_yr} years`, detail: `Equivalent annual cost is ${fmt(o.eac_min_usd_yr)} USD at ${o.economic_life_yr} years against ${fmt(o.eac_at_life_usd_yr)} USD at the planned ${i.life_yr} years.`, action: 'Plan the fleet roll-over around the economic life; a newer aircraft also brings lower fuel burn and emissions, which this capital-and-maintenance view does not yet credit.', basis: 'Minimum equivalent annual cost' });
    else if (o.economic_life_yr > i.life_yr + 2) out.push({ severity: 'info', title: 'Keeping the aircraft longer lowers annual cost', detail: `Economic life ${o.economic_life_yr} years exceeds the planned ${i.life_yr}.`, action: 'Check the structural life limit (Suite 9) before extending: life extension avoids the material and energy of a new airframe.', basis: 'Minimum equivalent annual cost' });
    return out;
  },
};

const DRIVERS = [['Fuel / energy price', 'fuel'], ['Demand (load factor, charter hours)', 'demand'], ['Yield / rate', 'yield'], ['Maintenance cost', 'maint'], ['Utilisation', 'util'], ['Carbon price', 'carbon']];
const risk = {
  id: 'risk', title: 'Financial risk: Monte Carlo, scenarios, hedging and the option to defer', fidelity: 'numerical',
  summary: 'Samples fuel price, demand, maintenance cost, utilisation, carbon price and residual value (with fuel and demand correlated) to get the distribution of NPV and the chance of loss, ranks the drivers, runs named scenarios, shows what a fuel hedge does, and values the right to wait.',
  equations: ['Stochastic financial risk equations', 'Net present value equation', 'Discounted cash flow equations', 'Fuel consumption–fuel price cost equations', 'Emissions–carbon pricing equations', 'Aircraft reliability–availability–revenue equations'],
  inputs: [...BASE, ...FIN_EXTRA,
    { key: 'vol_fuel', label: 'Fuel price uncertainty (σ of log)', unit: '-', default: 0.25, min: 0, max: 1, group: 'Uncertainty' }, { key: 'sd_demand', label: 'Demand uncertainty (σ)', unit: '-', default: 0.06, min: 0, max: 0.5, group: 'Uncertainty' },
    { key: 'rho_fuel_demand', label: 'Correlation of fuel price and demand', unit: '-', default: -0.3, min: -0.95, max: 0.95, group: 'Uncertainty', help: 'Negative: high fuel prices tend to coincide with weak demand' },
    { key: 'sd_maint', label: 'Maintenance cost uncertainty (σ)', unit: '-', default: 0.12, min: 0, max: 0.6, group: 'Uncertainty' }, { key: 'sd_util', label: 'Utilisation uncertainty (σ)', unit: '-', default: 0.05, min: 0, max: 0.4, group: 'Uncertainty' },
    { key: 'vol_carbon', label: 'Carbon price uncertainty (σ of log)', unit: '-', default: 0.35, min: 0, max: 1.5, group: 'Uncertainty' }, { key: 'sd_residual', label: 'Residual value uncertainty (σ, absolute fraction)', unit: '-', default: 0.05, min: 0, max: 0.4, group: 'Uncertainty' },
    { key: 'availability', label: 'Technical availability', unit: '-', default: 1, min: 0.5, max: 1, group: 'Uncertainty', help: 'From Suite 22: scales achievable utilisation and revenue' },
    { key: 'hedge_frac', label: 'Share of fuel hedged at today’s price', unit: '-', default: 0, min: 0, max: 1, group: 'Hedging and options' },
    { key: 'defer_yr', label: 'Period the investment could be deferred', unit: 'yr', default: 2, min: 0, max: 10, group: 'Hedging and options' }, { key: 'riskfree', label: 'Risk-free rate', unit: '-', default: 0.035, min: 0, max: 0.3, group: 'Hedging and options' },
    { key: 'nSamples', label: 'Monte Carlo samples', unit: '', default: 2000, min: 200, max: 100000, step: 1, discrete: true, group: 'Numerics' }, { key: 'seed', label: 'Random seed', unit: '', default: 26, min: 1, max: 1e9, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => ({ ...finDefaults(c, up, d), availability: up.safety?.dispatch_reliability > 0.5 ? Math.min(1, up.safety.dispatch_reliability) : up.safety?.availability > 0.5 ? up.safety.availability : undefined }),
  run(i) {
    const warnings = [], n = Math.round(i.nSamples), u = N.rng(i.seed), av = i.availability, npvM = (m, h = i.hedge_frac) => npvOf(cashFlows(i, { ...m, util: (m.util ?? 1) * av, hedge: h }).proj, i.discount), base = npvM({});
    const Lc = N.cholesky([[1, i.rho_fuel_demand], [i.rho_fuel_demand, 1]]), S = new Array(n), SH = new Array(n), Fz = new Array(n), Dz = new Array(n), other = 0.5;
    for (let s = 0; s < n; s++) {
      const z = [N.randn(u), N.randn(u)], zf = Lc[0][0] * z[0], zd = Lc[1][0] * z[0] + Lc[1][1] * z[1]; Fz[s] = zf; Dz[s] = zd;
      const m = { fuel: Math.exp(i.vol_fuel * zf - 0.5 * i.vol_fuel ** 2), demand: Math.max(0.3, 1 + i.sd_demand * zd), maint: Math.max(0.3, 1 + i.sd_maint * N.randn(u)), util: Math.max(0.3, 1 + i.sd_util * N.randn(u)), carbon: Math.exp(i.vol_carbon * N.randn(u) - 0.5 * i.vol_carbon ** 2), residual: N.clamp(i.residual_frac + i.sd_residual * N.randn(u), 0, 0.95) };
      S[s] = npvM(m); SH[s] = npvM(m, i.hedge_frac > 0 ? 0 : other);
    }
    const mean = N.mean(S), sd = N.std(S), pLoss = S.filter((v) => v < 0).length / n, var5 = N.quantile(S, 0.05), tail = S.filter((v) => v <= var5), cvar = tail.length ? N.mean(tail) : var5, hist = N.histogram(S, 30), sdAlt = N.std(SH);
    const sds = { fuel: i.vol_fuel, demand: i.sd_demand, yield: i.sd_demand, maint: i.sd_maint, util: i.sd_util, carbon: i.vol_carbon }, tor = DRIVERS.map(([nm, k]) => ({ nm, lo: npvM({ [k]: 1 - sds[k] }), hi: npvM({ [k]: 1 + sds[k] }) })).sort((a, b) => Math.abs(b.hi - b.lo) - Math.abs(a.hi - a.lo));
    const scen = [['Base', {}], ['High fuel (+50%)', { fuel: 1.5 }], ['Low demand (−10%)', { demand: 0.9 }], ['Carbon price ×2.5', { carbon: 2.5 }], ['High fuel and low demand', { fuel: 1.5, demand: 0.9 }], ['Favourable (fuel −20%, demand +5%)', { fuel: 0.8, demand: 1.05 }]].map(([nm, m]) => { const c = cashFlows(i, { ...m, util: av, hedge: i.hedge_frac }), v = npvOf(c.proj, i.discount); return [nm, v, irrOf(c.proj).irr, v - base]; });
    // deferral option: project value V (PV of operating flows) against the investment I; volatility from the simulated NPV spread
    const I = i.price_usd, V = base + I, sig = V > 0 ? N.clamp(sd / V / Math.sqrt(Math.max(1, i.life_yr / 4)), 0.02, 1) : 0, yld = V > 0 ? 1 / pvFactor(i.discount, Math.round(i.life_yr)) : 0, opt = V > 0 && I > 0 ? crrCall(V, I, i.riskfree, yld, sig, i.defer_yr, 200) : 0, wait = Math.max(0, opt - Math.max(base, 0));
    const ns = [200, 400, 800, 1600, 3200, 6400].filter((v) => v <= n), run = ns.map((k) => N.mean(S.slice(0, k)));
    const rev = hasRevenue(i);
    if (!rev) warnings.push('No revenue is defined: NPV is a present cost, "probability of loss" is 100% by construction, and the deferral option has no meaning.');
    if (av < 1) warnings.push(`Technical availability of ${(100 * av).toFixed(1)}% scales utilisation and revenue in every sample.`);
    warnings.push('Volatilities and the correlation are judgement inputs; the project volatility used for the option is derived from the simulated NPV spread and is an approximation.');
    return {
      kpis: [
        { key: 'npv_mean_usd', label: 'Mean NPV', value: mean, unit: 'USD', note: `± ${fmt(1.96 * sd / Math.sqrt(n))} (95% sampling error)` }, { key: 'npv_base_usd', label: 'Deterministic NPV', value: base, unit: 'USD' }, { key: 'npv_sd_usd', label: 'Standard deviation of NPV', value: sd, unit: 'USD' },
        { key: 'p_loss', label: 'Probability that NPV is negative', value: pLoss, unit: '-', status: !rev ? undefined : pLoss < 0.2 ? 'ok' : pLoss < 0.9 ? 'warn' : 'bad', note: 'Return below the discount rate; 50% for a project that just earns it' }, { key: 'var5_usd', label: 'NPV at the 5th percentile (value at risk)', value: var5, unit: 'USD' }, { key: 'cvar5_usd', label: 'Mean of the worst 5% (conditional value at risk)', value: cvar, unit: 'USD' },
        { key: 'npv_p95_usd', label: 'NPV at the 95th percentile', value: N.quantile(S, 0.95), unit: 'USD' }, { key: 'corr_fuel_demand_sample', label: 'Sampled fuel–demand correlation', value: N.corr(Fz, Dz), unit: '-' },
        { key: 'hedge_sd_change_pct', label: i.hedge_frac > 0 ? 'NPV spread reduction from the hedge' : 'NPV spread reduction if 50% of fuel were hedged', value: i.hedge_frac > 0 ? 100 * (1 - sd / sdAlt) : 100 * (1 - sdAlt / sd), unit: '%' },
        { key: 'top_driver_swing_usd', label: `Largest ±1σ swing: ${tor[0].nm}`, value: Math.abs(tor[0].hi - tor[0].lo), unit: 'USD' },
        { key: 'option_defer_usd', label: `Value of the option to invest within ${i.defer_yr} years`, value: opt, unit: 'USD' }, { key: 'value_of_waiting_usd', label: 'Value of waiting over investing now', value: wait, unit: 'USD', note: 'Option value less max(NPV, 0)' }, { key: 'project_vol', label: 'Project value volatility used', value: sig, unit: '1/√yr' },
        ...(i.fx_rate !== 1 || (i.currency && i.currency !== 'USD') ? [{ key: 'npv_mean_local', label: `Mean NPV in ${i.currency}`, value: mean * i.fx_rate, unit: i.currency }] : []),
      ],
      plots: [
        { type: 'bar', title: 'Distribution of NPV', ylabel: 'Samples [-]', categories: hist.centers.map((v) => (v / 1e6).toPrecision(3)), series: [{ name: 'NPV [USD M]', y: hist.counts }] },
        { type: 'bar', title: 'Tornado: NPV change for ±1σ of each driver', ylabel: 'NPV change from base [USD M]', categories: tor.map((t) => t.nm), series: [{ name: '−1σ', y: tor.map((t) => (t.lo - base) / 1e6) }, { name: '+1σ', y: tor.map((t) => (t.hi - base) / 1e6) }] },
        { type: 'bar', title: 'Scenario NPV', ylabel: 'NPV [USD M]', categories: scen.map((q) => q[0]), series: [{ name: 'NPV', y: scen.map((q) => q[1] / 1e6) }] },
        { type: 'line', title: 'Convergence of the mean NPV', xlabel: 'Samples [-]', ylabel: 'Mean NPV [USD M]', xlog: true, series: [{ name: 'Running mean', x: ns, y: run.map((v) => v / 1e6), style: 'line+points' }] },
      ],
      tables: [{ title: 'Scenarios', columns: ['Scenario', 'NPV [USD]', 'IRR [-]', 'Change from base [USD]'], rows: scen.map((q) => [q[0], q[1], Number.isFinite(q[2]) ? q[2] : '—', q[3]]) },
        { title: 'Driver sensitivity (±1σ, one at a time)', columns: ['Driver', 'σ', 'NPV at −1σ', 'NPV at +1σ', 'Swing'], rows: tor.map((t) => [t.nm, sds[DRIVERS.find((d) => d[0] === t.nm)[1]], t.lo, t.hi, Math.abs(t.hi - t.lo)]) }],
      warnings, models: ['Monte Carlo on the yearly cash-flow model (seeded)', 'Cholesky-correlated fuel price and demand; lognormal prices, normal quantities', 'Value at risk and conditional value at risk at 5%', 'One-at-a-time tornado and named scenarios', 'Fuel hedge as a fixed-price share', 'Cox–Ross–Rubinstein binomial lattice for the American option to defer'],
      assumptions: ['Each sample draws one level per driver that holds for the whole life (persistent uncertainty), which overstates year-to-year diversification risk but captures structural risk', 'Hedged fuel is bought at today’s price escalated like the rest', 'The deferral option treats the present value of operating cash flows as the underlying asset, the price as the strike and the annuity yield as the cost of waiting', PLACEHOLDER],
    };
  },
  convergence: { param: 'nSamples', label: 'Monte Carlo samples', levels: [250, 500, 1000, 2000, 4000], metric: 'npv_mean_usd', hOf: (n) => 1 / Math.sqrt(n) },
  verify() {
    const d = Object.fromEntries(risk.inputs.map((f) => [f.key, f.default])), z = N.kv(risk.run({ ...d, vol_fuel: 0, sd_demand: 0, sd_maint: 0, sd_util: 0, vol_carbon: 0, sd_residual: 0, nSamples: 200 })), o = N.kv(risk.run({ ...d, nSamples: 6000 }));
    return [
      N.check('Zero volatility reproduces the deterministic NPV', z.npv_mean_usd, z.npv_base_usd, 1e-9, 'Degenerate distribution'),
      N.check('Sampled correlation matches the input', o.corr_fuel_demand_sample, -0.3, 0.1, 'Cholesky factorisation; sampling error ≈ 0.012'),
      N.check('Binomial American call without yield equals Black–Scholes', crrCall(100, 100, 0.05, 0, 0.2, 1, 400), bsCall(100, 100, 0.05, 0.2, 1), 2e-3, 'No early exercise without dividends; lattice error O(1/n)'),
      N.check('Black–Scholes reference value', bsCall(100, 100, 0.05, 0.2, 1), 10.4506, 1e-4, 'Textbook value for S = K = 100, r = 5%, σ = 20%, T = 1'),
      N.check('Option with zero volatility is its intrinsic value', crrCall(120, 100, 0.05, 0, 0, 1, 50), 20, 1e-12, 'Limit'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], tor = res.tables[1].rows;
    // without revenue every sample is a cost: "probability of loss" carries no information
    if (!hasRevenue(i)) return [{ severity: 'info', title: 'Cost-of-ownership risk', detail: `Present cost of owning and operating ${fmt(-o.npv_mean_usd)} USD on average; 5%–95% range ${fmt(-o.npv_p95_usd)} … ${fmt(-o.var5_usd)} USD.`, action: `The largest driver is ${tor[0][0].toLowerCase()}. Enter a charter rate or yield to appraise a commercial operation.`, basis: 'Monte Carlo distribution of the present cost' }];
    // a project that just earns its cost of capital has a 50% chance of a negative NPV: only a near-certain loss is critical
    if (o.p_loss > 0.3) out.push({ severity: o.p_loss > 0.9 ? 'critical' : o.p_loss > 0.5 ? 'warn' : 'advise', title: o.p_loss > 0.5 ? `${(100 * o.p_loss).toFixed(0)}% chance that the investment earns less than its cost of capital` : `${(100 * o.p_loss).toFixed(0)}% chance of earning less than the cost of capital`, detail: `Mean NPV ${fmt(o.npv_mean_usd)} USD; worst-5% average ${fmt(o.cvar5_usd)} USD.`, action: `The largest driver is ${tor[0][0].toLowerCase()}: secure it contractually where possible (fuel hedge, capacity or power-by-the-hour agreements, guaranteed residual value).`, basis: 'Monte Carlo NPV distribution' });
    else out.push({ severity: 'info', title: `Probability of loss ${(100 * o.p_loss).toFixed(0)}%`, detail: `NPV 5%–95% range ${fmt(o.var5_usd)} … ${fmt(o.npv_p95_usd)} USD.`, action: `Monitor ${tor[0][0].toLowerCase()}, the driver with the largest swing.`, basis: 'Monte Carlo NPV distribution' });
    if (/Fuel/.test(tor[0][0]) || /Fuel/.test(tor[1][0])) out.push({ severity: 'advise', title: 'Fuel price is a leading risk', detail: `A 50% hedge would change the NPV spread by about ${o.hedge_sd_change_pct.toFixed(0)}%.`, action: 'Hedge part of the fuel, and reduce exposure at source: every per cent of fuel burn saved lowers both the cost and its volatility, and the CO₂.', basis: 'Tornado and hedge comparison' });
    if (o.value_of_waiting_usd > 0.02 * i.price_usd && hasRevenue(i)) out.push({ severity: 'advise', title: 'Waiting has value', detail: `The option to invest within ${i.defer_yr} years is worth ${fmt(o.option_defer_usd)} USD, ${fmt(o.value_of_waiting_usd)} USD more than investing today.`, action: 'Secure the option cheaply (delivery slots, purchase rights) rather than committing now, unless competitive or regulatory timing forces the decision.', basis: 'Real-option value by binomial lattice' });
    return out;
  },
};

const route = {
  id: 'route', title: 'Route, mission and fleet economics', fidelity: 'analytical',
  summary: 'For aeroplanes: fare, frequency and profit on a route with price-sensitive demand. For helicopters: cost per mission with standby and availability. For unmanned aircraft: cost per flight hour and per square kilometre surveyed.',
  equations: ['Cost–volume–profit equations', 'Contribution margin equations', 'Marginal cost equations', 'Aircraft utilisation equations', 'Load factor equations', 'Yield equations', 'Fleet scheduling–revenue optimisation equations', 'Aircraft reliability–availability–revenue equations', 'Profitability equations'],
  inputs: [
    { key: 'mode', label: 'Operation type', type: 'select', options: ['Airline route', 'Helicopter mission', 'UAV survey'], default: 'Airline route', group: 'Operation' },
    ...baseBut(['loan_yr', 'discount', 'inflation', 'tax', 'fx_rate', 'currency']),
    { key: 'route_km', label: 'Route length', unit: 'km', default: 1500, min: 1, group: 'Airline route' }, { key: 'demand_ref', label: 'Daily one-way demand at the reference fare', unit: 'pax/day', default: 600, min: 0, group: 'Airline route' },
    { key: 'fare_ref', label: 'Reference fare', unit: 'USD', default: 135, min: 0.01, group: 'Airline route' }, { key: 'elasticity', label: 'Price elasticity of demand', unit: '-', default: -1.4, min: -5, max: -0.1, group: 'Airline route', help: 'Route-level estimate of the IATA / InterVISTAS air-travel demand study (−1.4; national level −0.8). Leisure markets are more elastic than business markets.' },
    { key: 'fare_cap', label: 'Highest sustainable fare ÷ reference fare', unit: '-', default: 1.3, min: 0.3, max: 20, group: 'Airline route', help: 'Competition and substitutes cap the fare; set high for a monopoly route' },
    { key: 'lf_max', label: 'Highest achievable average load factor', unit: '-', default: 0.92, min: 0.3, max: 1, group: 'Airline route' }, { key: 'turn_min', label: 'Turnaround time', unit: 'min', default: 45, min: 5, max: 600, group: 'Airline route' }, { key: 'day_h', label: 'Operating day', unit: 'h', default: 16, min: 1, max: 24, group: 'Airline route' },
    { key: 'mission_fh', label: 'Flight hours per mission', unit: 'FH', default: 1.5, min: 0.05, group: 'Helicopter mission' }, { key: 'missions_yr', label: 'Missions per year', unit: '1/yr', default: 400, min: 1, group: 'Helicopter mission' },
    { key: 'standby_crew_usd_yr', label: 'Standby and specialist crew cost', unit: 'USD/yr', default: 900000, min: 0, group: 'Helicopter mission', help: 'Duty crews, hoist operators, medical crew, base costs not charged by the flight hour' },
    { key: 'availability', label: 'Aircraft availability', unit: '-', default: 0.95, min: 0.3, max: 1, group: 'Helicopter mission', help: 'From Suite 22' }, { key: 'coverage_req', label: 'Required mission coverage', unit: '-', default: 0.98, min: 0.5, max: 0.9999, group: 'Helicopter mission', help: 'Share of call-outs that must be met' },
    { key: 'survey_V_kmh', label: 'Survey ground speed', unit: 'km/h', default: 90, min: 1, group: 'UAV survey' }, { key: 'swath_m', label: 'Effective swath width', unit: 'm', default: 300, min: 1, group: 'UAV survey' }, { key: 'survey_eff', label: 'Survey efficiency (turns, transit, overlap)', unit: '-', default: 0.6, min: 0.05, max: 1, group: 'UAV survey' },
  ],
  defaults: (c, up, d) => { const b = baseDefaults(c, up, d), fhc = c.econ.cycles_yr > 0 ? c.econ.util_fh_yr / c.econ.cycles_yr : 1, stage = b.block_kmh * fhc, heli = c.meta.type === 'helicopter', uav = c.meta.type === 'uav';
    return onlyInputs({ ...b, mode: heli ? 'Helicopter mission' : uav ? 'UAV survey' : 'Airline route', route_km: stage, fare_ref: c.econ.yield_usd_pkm > 0 ? c.econ.yield_usd_pkm * stage : Math.max(1, (1.1 * opCost({ ...Object.fromEntries(BASE.map((f) => [f.key, f.default])), ...Object.fromEntries(Object.entries(b).filter(([, v]) => v !== undefined)) }).tocFh * fhc) / Math.max(1, c.econ.seats * c.econ.load_factor)), demand_ref: Math.max(1, c.econ.seats * c.econ.load_factor * (c.econ.yield_usd_pkm > 0 ? 4 : 1)), turn_min: c.mass.mtow_kg > 40000 ? 45 : c.mass.mtow_kg > 5700 ? 30 : 20,
      mission_fh: heli ? fhc : 1.5, missions_yr: heli ? c.econ.cycles_yr : 400, standby_crew_usd_yr: heli ? 0.04 * (c.econ.price_usd || 0) : 0, availability: up.safety?.availability > 0.3 ? Math.min(up.safety.availability, 0.999) * 0.96 : 0.95, survey_V_kmh: (c.mission.cruise_V_ms || c.flight.V_ms || 20) * 3.6 * 0.9, swath_m: uav ? Math.max(50, 0.3 * (c.mission.cruise_alt_m || 100) * 2) : 300 }, route); },
  run(i) {
    const o = opCost(i), warnings = [], s = i.sales_pct_rev / 100; let kpis, plots, tables, models;
    if (i.mode === 'Airline route') {
      if (!(i.seats > 0)) { warnings.push('No seats are defined: the airline route model does not apply. Select the helicopter or UAV mode, or enter seats.'); }
      const seats = Math.max(i.seats, 1), tBlock = i.route_km / i.block_kmh, fixedFh = o.doc + o.admin / o.util, navFix = (i.nav_usd_100km * (i.route_km / 100) * Math.sqrt(i.mtow_t / 50)), cFlight = (fixedFh - o.cats['Airport and handling'] - o.cats.Navigation) * tBlock + (i.landing_usd_t * i.mtow_t + i.handling_usd_cycle + i.env_usd_cycle) + navFix;
      const Qd = (p) => i.demand_ref * (p / i.fare_ref) ** i.elasticity, day = (f, p) => { const pax = Math.min(Qd(p), f * seats * i.lf_max); return { pax, profit: pax * (p * (1 - s) - i.ioc_usd_pax) - f * cFlight, lf: pax / (f * seats), spill: Qd(p) - pax }; };
      const fmax = Math.max(1, Math.min(40, Math.ceil((3 * i.demand_ref) / (seats * i.lf_max)))), res = N.range(fmax, (k) => { const f = k + 1, p = N.goldenSection((q) => -day(f, q).profit, Math.min(0.2, 0.5 * i.fare_cap) * i.fare_ref, Math.max(i.fare_cap, 0.3) * i.fare_ref, 1e-8), d = day(f, p); return { f, p, ...d }; }), best = res[N.argmax(res.map((r) => r.profit))];
      const acNeeded = (2 * best.f * (tBlock * 60 + i.turn_min)) / (i.day_h * 60), cm = best.p * (1 - s) - i.ioc_usd_pax, fares = N.linspace(0.4 * i.fare_ref, Math.min(2.5, Math.max(1.2, i.fare_cap)) * i.fare_ref, 40), pStar = i.elasticity < -1 ? ((i.ioc_usd_pax / (1 - s)) * i.elasticity) / (1 + i.elasticity) : Infinity;
      if (i.elasticity >= -1) warnings.push('Demand is inelastic (elasticity between −1 and 0): revenue rises without limit as the fare rises in this constant-elasticity model; the optimum sits on the fare cap.');
      if (best.p >= i.fare_cap * i.fare_ref * (1 - 1e-6)) warnings.push('The best fare sits on the fare cap: the result depends on how much pricing power the route really has.');
      if (best.profit < 0) warnings.push('No fare and frequency combination makes this route profitable with this aircraft and demand.');
      kpis = [
        { key: 'route_profit_day_usd', label: 'Best daily route profit (one direction)', value: best.profit, unit: 'USD/day', status: best.profit >= 0 ? 'ok' : 'bad' }, { key: 'opt_fare_usd', label: 'Profit-maximising fare', value: best.p, unit: 'USD' }, { key: 'opt_frequency', label: 'Best daily frequency', value: best.f, unit: 'flights/day' },
        { key: 'route_load_factor', label: 'Load factor at the optimum', value: best.lf, unit: '-' }, { key: 'pax_day', label: 'Passengers carried per day', value: best.pax, unit: 'pax/day' }, { key: 'spill_pax_day', label: 'Demand turned away', value: best.spill, unit: 'pax/day' },
        { key: 'cost_per_flight_usd', label: 'Cost per flight (excluding passenger-variable cost)', value: cFlight, unit: 'USD' }, { key: 'contribution_per_pax_usd', label: 'Contribution margin per passenger', value: cm, unit: 'USD' }, { key: 'breakeven_pax_per_flight', label: 'Break-even passengers per flight', value: cm > 0 ? cFlight / cm : NaN, unit: 'pax' },
        { key: 'aircraft_needed', label: 'Aircraft needed for the rotation', value: acNeeded, unit: 'aircraft', note: 'Both directions, block plus turnaround, in the operating day' }, { key: 'route_margin_pct', label: 'Route margin', value: best.pax > 0 ? (100 * best.profit) / (best.pax * best.p) : NaN, unit: '%' },
        { key: 'annual_route_profit_usd', label: 'Annual route profit (both directions)', value: 2 * 365 * best.profit, unit: 'USD/yr' }, { key: 'monopoly_fare_usd', label: 'Unconstrained profit-maximising fare', value: pStar, unit: 'USD', note: 'c·ε/(1 + ε) for constant elasticity' },
      ];
      plots = [{ type: 'line', title: 'Daily profit versus fare at the best frequency', xlabel: 'Fare [USD]', ylabel: 'Profit [USD k/day]', series: [{ name: `${best.f} flights/day`, x: fares, y: fares.map((p) => day(best.f, p).profit / 1e3) }, ...(best.f > 1 ? [{ name: `${best.f - 1} flights/day`, x: fares, y: fares.map((p) => day(best.f - 1, p).profit / 1e3), style: 'dash' }] : []), { name: `${best.f + 1} flights/day`, x: fares, y: fares.map((p) => day(best.f + 1, p).profit / 1e3), style: 'dash' }], annotations: [{ x: best.p, label: 'Optimum' }, { y: 0, label: 'Break-even' }] },
        { type: 'line', title: 'Profit and load factor versus frequency (best fare at each)', xlabel: 'Flights per day [-]', ylabel: 'Profit [USD k/day]', series: [{ name: 'Profit', x: res.map((r) => r.f), y: res.map((r) => r.profit / 1e3), style: 'line+points' }] },
        { type: 'line', title: 'Demand curve and capacity', xlabel: 'Fare [USD]', ylabel: 'Passengers per day [-]', series: [{ name: 'Demand', x: fares, y: fares.map(Qd) }, { name: 'Capacity at best frequency', x: [fares[0], fares[fares.length - 1]], y: [best.f * seats * i.lf_max, best.f * seats * i.lf_max], style: 'dash' }] }];
      tables = [{ title: 'Frequency options', columns: ['Flights/day', 'Best fare [USD]', 'Passengers/day', 'Load factor', 'Spill [pax/day]', 'Profit [USD/day]'], rows: res.slice(0, 12).map((r) => [r.f, r.p, r.pax, r.lf, r.spill, r.profit]) }];
      models = ['Constant-elasticity demand curve', 'Capacity-constrained fare optimisation by golden-section search', 'Frequency enumeration', 'Rotation count from block and turnaround time'];
    } else if (i.mode === 'Helicopter mission') {
      const varFh = o.cash - o.cats.Insurance, fixedYr = (o.cats.Insurance + o.dep + o.interest) * o.util + i.standby_crew_usd_yr + o.admin, fhYr = i.missions_yr * i.mission_fh, A = N.clamp(i.availability, 0.01, 0.9999);
      const nAc = Math.max(1, Math.ceil(Math.log(1 - i.coverage_req) / Math.log(1 - A) - 1e-9)), cover = 1 - (1 - A) ** nAc, fleetFixed = fixedYr * nAc, total = fleetFixed + varFh * fhYr, perMission = total / (i.missions_yr * cover), ms = N.linspace(0.25 * i.missions_yr, 3 * i.missions_yr, 30);
      if (fhYr > o.util * nAc * 1.001) warnings.push('Mission hours exceed the stated annual utilisation of the fleet: raise utilisation or the number of aircraft.');
      kpis = [
        { key: 'cost_per_mission_usd', label: 'Cost per completed mission', value: perMission, unit: 'USD' }, { key: 'cost_per_op_fh_usd', label: 'Cost per operational flight hour', value: total / (fhYr * cover), unit: 'USD/FH' }, { key: 'variable_cost_usd_fh', label: 'Variable cost per flight hour', value: varFh, unit: 'USD/FH' },
        { key: 'fixed_cost_usd_yr', label: 'Fixed cost of readiness per year (fleet)', value: fleetFixed, unit: 'USD/yr' }, { key: 'standby_share_pct', label: 'Share of cost that is readiness (fixed)', value: (100 * fleetFixed) / total, unit: '%' },
        { key: 'aircraft_needed', label: 'Aircraft needed for the required coverage', value: nAc, unit: 'aircraft', note: `Coverage 1 − (1 − A)ⁿ = ${(100 * cover).toFixed(2)}%` }, { key: 'coverage_achieved', label: 'Mission coverage achieved', value: cover, unit: '-', status: cover >= i.coverage_req ? 'ok' : 'warn' },
        { key: 'missed_missions_yr', label: 'Missions missed per year', value: i.missions_yr * (1 - cover), unit: '1/yr' }, { key: 'annual_contract_cost_usd', label: 'Annual cost of the service', value: total, unit: 'USD/yr' }, { key: 'breakeven_contract_usd_yr', label: 'Break-even annual contract value', value: total / (1 - s), unit: 'USD/yr' },
      ];
      plots = [{ type: 'line', title: 'Cost per mission versus missions per year', xlabel: 'Missions per year [-]', ylabel: 'Cost per mission [USD]', series: [{ name: 'Cost per completed mission', x: ms, y: ms.map((m) => (fleetFixed + varFh * m * i.mission_fh) / (m * cover)) }], annotations: [{ x: i.missions_yr, label: 'Planned' }] },
        { type: 'bar', title: 'Annual cost structure', ylabel: 'Cost [USD M/yr]', categories: ['Service'], stacked: true, series: [{ name: 'Variable flight cost', y: [(varFh * fhYr) / 1e6] }, { name: 'Ownership and insurance', y: [((o.cats.Insurance + o.dep + o.interest) * o.util * nAc) / 1e6] }, { name: 'Standby and specialist crew', y: [(i.standby_crew_usd_yr * nAc) / 1e6] }, { name: 'Administration', y: [(o.admin * nAc) / 1e6] }] },
        { type: 'line', title: 'Coverage versus fleet size', xlabel: 'Aircraft [-]', ylabel: 'Coverage [-]', series: [{ name: 'Coverage', x: [1, 2, 3, 4, 5], y: [1, 2, 3, 4, 5].map((k) => 1 - (1 - A) ** k), style: 'line+points' }], annotations: [{ y: i.coverage_req, label: 'Required' }] }];
      tables = [{ title: 'Mission cost build-up', columns: ['Item', 'USD per year', 'USD per mission'], rows: [['Variable flight cost', varFh * fhYr, (varFh * fhYr) / (i.missions_yr * cover)], ['Readiness: ownership, insurance, administration', fleetFixed - i.standby_crew_usd_yr * nAc, (fleetFixed - i.standby_crew_usd_yr * nAc) / (i.missions_yr * cover)], ['Readiness: standby and specialist crew', i.standby_crew_usd_yr * nAc, (i.standby_crew_usd_yr * nAc) / (i.missions_yr * cover)], ['Total', total, perMission]] }];
      models = ['Fixed (readiness) and variable (flight-hour) cost split', 'Independent-aircraft coverage 1 − (1 − A)ⁿ', 'Cost per completed mission'];
    } else {
      const rate = (i.survey_V_kmh * i.swath_m) / 1e3 * i.survey_eff, perKm2 = rate > 0 ? o.tocFh / rate : NaN, Vs = N.linspace(0.5 * i.survey_V_kmh, 1.5 * i.survey_V_kmh, 20);
      kpis = [
        { key: 'cost_per_km2_usd', label: 'Cost per square kilometre surveyed', value: perKm2, unit: 'USD/km²' }, { key: 'area_rate_km2_h', label: 'Area coverage rate', value: rate, unit: 'km²/FH' }, { key: 'cost_per_fh_usd', label: 'Total cost per flight hour', value: o.tocFh, unit: 'USD/FH' },
        { key: 'cost_per_line_km_usd', label: 'Cost per line-kilometre', value: o.tocFh / (i.survey_V_kmh * i.survey_eff), unit: 'USD/km' }, { key: 'annual_area_km2', label: 'Area surveyed per year', value: rate * o.util, unit: 'km²/yr' },
        { key: 'energy_cost_share_pct', label: 'Energy share of cost', value: (100 * o.fuelFh) / o.tocFh, unit: '%' }, { key: 'crew_cost_share_pct', label: 'Crew share of cost', value: (100 * i.crew_usd_fh) / o.tocFh, unit: '%' },
      ];
      plots = [{ type: 'line', title: 'Cost per square kilometre versus survey speed', xlabel: 'Ground speed [km/h]', ylabel: 'Cost [USD/km²]', series: [{ name: 'At constant cost per hour', x: Vs, y: Vs.map((v) => o.tocFh / ((v * i.swath_m) / 1e3 * i.survey_eff)) }], annotations: [{ x: i.survey_V_kmh, label: 'Planned' }] },
        { type: 'bar', title: 'Cost per flight hour', ylabel: 'Cost [USD/FH]', categories: ['Energy', 'Crew', 'Maintenance', 'Ownership and insurance', 'Other'], series: [{ name: 'USD/FH', y: [o.fuelFh + o.carbonFh, i.crew_usd_fh, o.maintFh, o.dep + o.interest + o.cats.Insurance, o.tocFh - o.fuelFh - o.carbonFh - i.crew_usd_fh - o.maintFh - o.dep - o.interest - o.cats.Insurance] }] }];
      tables = [{ title: 'Survey productivity', columns: ['Quantity', 'Value', 'Unit'], rows: [['Ground speed', i.survey_V_kmh, 'km/h'], ['Swath', i.swath_m, 'm'], ['Efficiency', i.survey_eff, '-'], ['Area rate', rate, 'km²/FH'], ['Cost', perKm2, 'USD/km²']] }];
      models = ['Area rate = speed × swath × efficiency', 'Cost per flight hour from the operating-cost model'];
    }
    warnings.push(`${PLACEHOLDER}. Demand, elasticity and coverage inputs are market and contract data you must supply.`);
    return { kpis, plots, tables, warnings, models, assumptions: ['Costs per flight hour come from the same operating-cost model as the cost build-up', i.mode === 'Airline route' ? 'One fare class and one demand curve per direction; competition enters only through the fare cap; no network or connecting traffic' : i.mode === 'Helicopter mission' ? 'Aircraft are available independently; call-outs arrive one at a time' : 'Survey lines are flown at constant speed; efficiency covers turns, transit and overlap', PLACEHOLDER] };
  },
  verify() {
    const d = { ...Object.fromEntries(route.inputs.map((f) => [f.key, f.default])), mode: 'Airline route', demand_ref: 5, lf_max: 1, fare_cap: 10 }, r = route.run(d), o = N.kv(r), one = r.tables[0].rows[0], s = 0.04;
    const h = N.kv(route.run({ ...d, mode: 'Helicopter mission', availability: 0.9, coverage_req: 0.98 })), u = N.kv(route.run({ ...d, mode: 'UAV survey', survey_V_kmh: 100, swath_m: 200, survey_eff: 0.5 }));
    return [
      N.check('Unconstrained optimum fare = c·ε/(1 + ε)', one[1], ((8 / (1 - s)) * -1.4) / (1 - 1.4), 1e-5, 'Monopoly pricing with constant elasticity (capacity not binding)'),
      N.check('Coverage needs 2 aircraft at 90% availability for 98%', h.aircraft_needed, 2, 0, '1 − 0.1² = 0.99 ≥ 0.98'),
      N.check('Survey area rate = speed × swath × efficiency', u.area_rate_km2_h, 100 * 0.2 * 0.5, 1e-12, 'Definition'),
      N.check('Monopoly fare KPI', o.monopoly_fare_usd, 8 / 0.96 * 3.5, 1e-12, 'Closed form'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (i.mode === 'Airline route') {
      out.push({ severity: o.route_profit_day_usd >= 0 ? 'advise' : 'warn', title: o.route_profit_day_usd >= 0 ? `Fly ${o.opt_frequency} a day at about ${o.opt_fare_usd.toFixed(0)} USD` : 'The route does not pay with this aircraft', detail: `Load factor ${(100 * o.route_load_factor).toFixed(0)}%, ${o.pax_day.toFixed(0)} passengers and ${fmt(o.route_profit_day_usd)} USD a day each way; ${o.aircraft_needed.toFixed(2)} aircraft tied up.`, action: o.route_profit_day_usd >= 0 ? (o.spill_pax_day > 0.1 * o.pax_day ? 'Demand is being turned away: test a larger aircraft or one more frequency.' : 'Hold frequency; use the spare seats for connecting or promotional traffic.') : 'Use a smaller or cheaper-to-operate aircraft, combine with another route, or drop it.', basis: 'Fare and frequency optimisation on a constant-elasticity demand curve' });
      if (o.route_load_factor < 0.6) out.push({ severity: 'advise', title: 'The aircraft is too large for this market', detail: `Best load factor is only ${(100 * o.route_load_factor).toFixed(0)}%.`, action: 'A right-sized aircraft cuts cost and fuel burn per passenger; empty seats carry CO₂ for nothing.', basis: 'Capacity matching' });
    } else if (i.mode === 'Helicopter mission') {
      out.push({ severity: 'info', title: `${fmt(o.cost_per_mission_usd)} USD per mission, ${o.standby_share_pct.toFixed(0)}% of it readiness`, detail: `${o.aircraft_needed} aircraft give ${(100 * o.coverage_achieved).toFixed(1)}% coverage; ${o.missed_missions_yr.toFixed(1)} missions a year would be missed.`, action: o.standby_share_pct > 60 ? 'Readiness dominates: share standby aircraft and crews across bases or contracts, and raise availability (Suite 22) before adding airframes.' : 'Flight-hour cost dominates: work on fuel burn, maintenance cost per hour and mission routing.', basis: 'Fixed and variable cost split with availability' });
      if (o.aircraft_needed > 1) out.push({ severity: 'advise', title: 'Availability drives fleet size', detail: `Availability ${(100 * i.availability).toFixed(1)}% forces ${o.aircraft_needed} aircraft for ${(100 * i.coverage_req).toFixed(1)}% coverage.`, action: 'Each point of availability gained may remove a whole standby aircraft; price the maintenance and spares investment against that.', basis: 'Reliability–availability–revenue coupling' });
    } else out.push({ severity: 'info', title: `${Number(o.cost_per_km2_usd).toFixed(2)} USD per km² surveyed`, detail: `${o.area_rate_km2_h.toFixed(1)} km² per flight hour; crew is ${o.crew_cost_share_pct.toFixed(0)}% and energy ${o.energy_cost_share_pct.toFixed(0)}% of cost.`, action: o.crew_cost_share_pct > 40 ? 'Crew dominates: one operator supervising several aircraft, or longer-endurance sorties, cut cost most.' : 'Raise swath (sensor, altitude) or survey efficiency before raising speed, which costs energy.', basis: 'Productivity-based unit cost' });
    return out;
  },
};

const sustain = {
  id: 'sustain', title: 'Sustainability economics: carbon cost, energy carriers and abatement cost', fidelity: 'analytical',
  summary: 'Compares kerosene, sustainable fuel blends, hydrogen and battery electricity on cost and CO₂ for the same useful energy, gives the cost of each tonne of CO₂ avoided, and projects the carbon bill over the aircraft life.',
  equations: ['Emissions–carbon pricing equations', 'Fuel consumption–fuel price cost equations', 'Hybrid-electric propulsion–energy cost equations', 'Discounted cash flow equations', 'Cost per passenger-kilometre equations'],
  inputs: [...baseOf(['fuel_kg_fh', 'fuel_usd_kg', 'fuel', 'energy_kWh_fh', 'elec_usd_kWh', 'saf_blend', 'saf_price_ratio', 'carbon_usd_t', 'carbon_cover', 'util_fh_yr', 'cycles_yr', 'block_kmh', 'seats', 'load_factor', 'life_yr', 'discount']),
    { key: 'carbon_escal', label: 'Carbon price escalation (nominal)', unit: '-', default: 0.05, min: -0.1, max: 0.5, group: 'Energy' },
    { key: 'h2_usd_kg', label: 'Hydrogen price (delivered, liquid)', unit: 'USD/kg', default: 6, min: 0, group: 'Alternative carriers (illustrative)', help: 'Unsourced placeholder for delivered liquid hydrogen. The only verified figure is the US DOE production-cost target of 1 USD/kg.' }, { key: 'h2_co2_kg_kg', label: 'Hydrogen production CO₂ intensity', unit: 'kg CO₂/kg H₂', default: 1, min: 0, max: 30, group: 'Alternative carriers (illustrative)', help: 'Near 0–1 for electrolysis on renewable power; about 10 for unabated natural-gas reforming' },
    { key: 'h2_energy_penalty', label: 'Hydrogen aircraft energy penalty', unit: '-', default: 1.1, min: 0.8, max: 2, group: 'Alternative carriers (illustrative)', help: 'Extra energy for tank mass and volume' },
    { key: 'eta_fuel', label: 'Fuel-to-thrust-power efficiency', unit: '-', default: 0.35, min: 0.1, max: 0.6, group: 'Alternative carriers (illustrative)', help: 'Overall efficiency of the combustion powerplant; from Suite 7 when available' }, { key: 'eta_batt', label: 'Grid-to-thrust-power efficiency of a battery aircraft', unit: '-', default: 0.72, min: 0.3, max: 0.95, group: 'Alternative carriers (illustrative)' },
    { key: 'grid_g_kWh', label: 'Grid carbon intensity', unit: 'g CO₂/kWh', default: 367, min: 0, max: 1200, group: 'Alternative carriers (illustrative)', help: 'Regional. The default is the 2023 US average of 0.81 lb CO₂ per kWh (EIA); enter the value of the grid that charges the aircraft.' }, { key: 'batt_energy_penalty', label: 'Battery aircraft energy penalty', unit: '-', default: 1.3, min: 0.8, max: 4, group: 'Alternative carriers (illustrative)', help: 'Extra energy for battery mass; mission feasibility is checked in Suite 24, not here' },
  ],
  defaults: (c, up, d) => { const b = baseDefaults(c, up, d); return { ...Object.fromEntries(['fuel_kg_fh', 'fuel_usd_kg', 'fuel', 'energy_kWh_fh', 'elec_usd_kWh', 'carbon_usd_t', 'util_fh_yr', 'cycles_yr', 'block_kmh', 'seats', 'load_factor', 'life_yr', 'discount'].map((k) => [k, b[k]])), carbon_escal: c.econ.inflation + 0.025,
      // a battery aircraft has no combustion powerplant of its own: the comparison then uses a turboprop-class efficiency, not the electric drive's
      eta_fuel: c.prop.type === 'electric' ? 0.28 : N.clamp(up.propulsion?.eta_overall ?? (c.prop.type === 'piston' ? 0.24 : c.prop.type === 'turbofan' ? 0.35 : 0.28), 0.1, 0.6) }; },
  run(i) {
    const warnings = [], jet = i.fuel === 'SAF (HEFA-SPK)' || i.fuel === 'Liquid hydrogen' ? FLUIDS['Jet A-1'] : FLUIDS[i.fuel] || FLUIDS['Jet A-1'], saf = FLUIDS['SAF (HEFA-SPK)'], h2 = FLUIDS['Liquid hydrogen'], elecBase = !(i.fuel_kg_fh > 0) && i.energy_kWh_fh > 0;
    // useful (thrust) energy per flight hour that every carrier must deliver [kWh]
    // (the penalty of a carrier the aircraft already uses is inside its measured consumption, so it is taken out here and put back below)
    const isH2 = i.fuel === 'Liquid hydrogen', isSaf = i.fuel === 'SAF (HEFA-SPK)', Eu = elecBase ? (i.energy_kWh_fh * i.eta_batt) / i.batt_energy_penalty : (i.fuel_kg_fh * (FLUIDS[i.fuel]?.LHV || jet.LHV) * i.eta_fuel) / 3.6e6 / (isH2 ? i.h2_energy_penalty : 1), cp = i.carbon_usd_t * i.carbon_cover / 1e3;
    // price of fossil fuel: the entered price is that of the selected fuel, so it is divided by the SAF ratio when SAF is selected
    const pFossil = isSaf ? i.fuel_usd_kg / i.saf_price_ratio : i.fuel_usd_kg, pH2 = isH2 ? i.fuel_usd_kg : i.h2_usd_kg;
    if (isH2) warnings.push('Liquid hydrogen is the selected fuel: its price is the fuel price entered, and the kerosene options are priced at that same figure per kilogram. Enter a kerosene price in a separate run for a fair fossil comparison.');
    if (!(Eu > 0)) warnings.push('Neither fuel burn nor electrical energy use is defined; all results are zero.');
    const mJet = (Eu * 3.6e6) / (jet.LHV * i.eta_fuel), fuelOpt = (x, nm) => { const lhv = (1 - x) * jet.LHV + x * saf.LHV, m = (Eu * 3.6e6) / (lhv * i.eta_fuel); return { nm, energy: m * pFossil * (1 - x + x * i.saf_price_ratio), co2: m * ((1 - x) * jet.co2_per_kg + x * saf.co2_per_kg * saf.lifecycle_factor), mass: m }; };
    const mH = (Eu * 3.6e6 * i.h2_energy_penalty) / (h2.LHV * i.eta_fuel), eB = (Eu * i.batt_energy_penalty) / i.eta_batt;
    const opts = [fuelOpt(0, jet === FLUIDS['Avgas 100LL'] ? 'Avgas' : 'Fossil kerosene'), ...(i.saf_blend > 0 && i.saf_blend < 1 ? [fuelOpt(i.saf_blend, `SAF blend ${(100 * i.saf_blend).toFixed(0)}%`)] : [fuelOpt(0.3, 'SAF blend 30%')]), fuelOpt(1, 'SAF 100%'), { nm: 'Liquid hydrogen', energy: mH * pH2, co2: mH * i.h2_co2_kg_kg, mass: mH }, { nm: 'Battery electric', energy: eB * i.elec_usd_kWh, co2: (eB * i.grid_g_kWh) / 1e3, mass: 0 }];
    const ref = opts[0]; opts.forEach((q) => { q.carbon = q.co2 * cp; q.total = q.energy + q.carbon; q.abate = ref.co2 - q.co2 > 1e-9 ? (1e3 * (q.energy - ref.energy)) / (ref.co2 - q.co2) : NaN; });
    const cur = elecBase ? opts[4] : i.fuel === 'Liquid hydrogen' ? opts[3] : i.fuel === 'SAF (HEFA-SPK)' ? opts[2] : i.saf_blend > 0 ? opts[1] : opts[0], n = Math.round(i.life_yr), yrs = N.range(n, (t) => t + 1), cy = yrs.map((t) => cur.co2 * i.util_fh_yr * cp * (1 + i.carbon_escal) ** t), pvC = N.sum(cy.map((v, k) => v / (1 + i.discount) ** (k + 1)));
    const fhc = i.util_fh_yr / Math.max(i.cycles_yr, 1e-9), rpk = i.seats * i.load_factor * i.block_kmh, gpk = rpk > 0 ? (1e3 * cur.co2) / rpk : NaN, cheapest = opts.slice(1).filter((q) => Number.isFinite(q.abate)).sort((a, b) => a.abate - b.abate)[0], cps = N.linspace(0, Math.max(400, 3 * i.carbon_usd_t), 30);
    if (!(rpk > 0)) warnings.push('No seats or load factor are defined: CO₂ per passenger-kilometre does not apply.');
    if (elecBase) warnings.push('The aircraft is battery-electric: the fuel options show what the same useful energy would cost and emit with a combustion powerplant of the stated efficiency.');
    warnings.push('Alternative-carrier prices, efficiencies and penalties are illustrative. Hydrogen and battery options change the aircraft itself; whether the mission is feasible is a Suite 23 / 24 question.');
    return {
      kpis: [
        { key: 'co2_kg_fh', label: 'Life-cycle CO₂ of the current energy carrier', value: cur.co2, unit: 'kg/FH' }, { key: 'co2_g_pkm', label: 'CO₂ per passenger-kilometre', value: gpk, unit: 'g/pkm' }, { key: 'co2_t_yr', label: 'CO₂ per year', value: (cur.co2 * i.util_fh_yr) / 1e3, unit: 't/yr' },
        { key: 'energy_cost_usd_fh', label: 'Energy cost of the current carrier', value: cur.energy, unit: 'USD/FH' }, { key: 'carbon_cost_fh_usd', label: 'Carbon cost of the current carrier', value: cur.carbon, unit: 'USD/FH' }, { key: 'energy_cost_per_flight_usd', label: 'Energy and carbon cost per flight', value: cur.total * fhc, unit: 'USD' },
        { key: 'carbon_cost_life_pv_usd', label: 'Present value of the carbon bill over the life', value: pvC, unit: 'USD' }, { key: 'carbon_cost_last_year_usd', label: 'Carbon bill in the final year', value: cy[n - 1] || 0, unit: 'USD/yr' },
        { key: 'abatement_saf_usd_t', label: 'Abatement cost: 100% SAF', value: opts[2].abate, unit: 'USD/t CO₂' }, { key: 'abatement_h2_usd_t', label: 'Abatement cost: hydrogen', value: opts[3].abate, unit: 'USD/t CO₂' }, { key: 'abatement_elec_usd_t', label: 'Abatement cost: battery electric', value: opts[4].abate, unit: 'USD/t CO₂' },
        { key: 'cheapest_abatement_usd_t', label: `Cheapest abatement: ${cheapest ? cheapest.nm : 'none'}`, value: cheapest ? cheapest.abate : NaN, unit: 'USD/t CO₂', status: cheapest && cheapest.abate <= i.carbon_usd_t * i.carbon_cover ? 'ok' : undefined, note: 'Pays for itself when below the effective carbon price' },
        { key: 'useful_energy_kWh_fh', label: 'Useful (thrust) energy required', value: Eu, unit: 'kWh/FH' }, { key: 'fossil_fuel_equiv_kg_fh', label: 'Fossil-fuel equivalent burn', value: mJet, unit: 'kg/FH' },
      ],
      plots: [
        { type: 'bar', title: 'Energy and carbon cost per flight hour by carrier', ylabel: 'Cost [USD/FH]', categories: opts.map((q) => q.nm), stacked: true, series: [{ name: 'Energy', y: opts.map((q) => q.energy) }, { name: 'Carbon', y: opts.map((q) => q.carbon) }] },
        { type: 'bar', title: 'Life-cycle CO₂ per flight hour by carrier', ylabel: 'CO₂ [kg/FH]', categories: opts.map((q) => q.nm), series: [{ name: 'CO₂', y: opts.map((q) => q.co2) }] },
        { type: 'line', title: 'Total energy cost versus carbon price', xlabel: 'Carbon price [USD/t CO₂]', ylabel: 'Energy + carbon cost [USD/FH]', series: opts.map((q) => ({ name: q.nm, x: cps, y: cps.map((p) => q.energy + (q.co2 * p * i.carbon_cover) / 1e3) })), annotations: [{ x: i.carbon_usd_t, label: 'Today' }] },
        { type: 'line', title: 'Carbon bill over the service life', xlabel: 'Year', ylabel: 'Carbon cost [USD k/yr, nominal]', series: [{ name: 'Current carrier', x: yrs, y: cy.map((v) => v / 1e3) }] },
      ],
      tables: [{ title: 'Energy carrier comparison (same useful energy)', columns: ['Carrier', 'Fuel mass [kg/FH]', 'Energy cost [USD/FH]', 'CO₂ [kg/FH]', 'Carbon cost [USD/FH]', 'Total [USD/FH]', 'CO₂ saving [%]', 'Abatement cost [USD/t]', 'Break-even carbon price [USD/t]'], rows: opts.map((q) => [q.nm, q.mass, q.energy, q.co2, q.carbon, q.total, ref.co2 > 0 ? 100 * (1 - q.co2 / ref.co2) : 0, Number.isFinite(q.abate) ? q.abate : '—', Number.isFinite(q.abate) ? q.abate / Math.max(i.carbon_cover, 1e-9) : '—']) }],
      outputs: { current_carrier: cur.nm },
      warnings, models: ['Equal-useful-energy comparison of energy carriers', 'Life-cycle CO₂ factors: combustion factor for fossil fuel, life-cycle factor for SAF, production intensity for hydrogen, grid intensity for electricity', 'Abatement cost = extra energy cost ÷ CO₂ avoided', 'Escalating carbon-price trajectory, discounted'],
      assumptions: ['Fossil fuel is counted at its combustion CO₂ only (upstream emissions excluded), SAF at its life-cycle factor', 'Hydrogen burns in a powerplant of the same efficiency with an energy penalty for tanks', 'Non-CO₂ climate effects (contrails, NOx) are not priced', 'Aircraft capital and maintenance differences between carriers are not included: this is the energy and carbon bill only', PLACEHOLDER],
    };
  },
  verify() {
    const d = Object.fromEntries(sustain.inputs.map((f) => [f.key, f.default])), r = sustain.run(d), o = N.kv(r), J = FLUIDS['Jet A-1'], S = FLUIDS['SAF (HEFA-SPK)'], k = J.LHV / S.LHV;
    return [
      N.check('SAF abatement cost in closed form', o.abatement_saf_usd_t, (1e3 * d.fuel_usd_kg * (d.saf_price_ratio * k - 1)) / (J.co2_per_kg - S.co2_per_kg * S.lifecycle_factor * k), 1e-10, 'Extra cost per kg of kerosene replaced ÷ CO₂ avoided'),
      N.check('Fossil reference reproduces the entered fuel burn', o.fossil_fuel_equiv_kg_fh, d.fuel_kg_fh, 1e-12, 'Energy equivalence'),
      N.check('A hydrogen aircraft is compared at its own entered burn', sustain.run({ ...d, fuel: 'Liquid hydrogen', fuel_kg_fh: 900 }).tables[0].rows[3][1], 900, 1e-12, 'The tank penalty is already inside the measured consumption'),
      N.check('SAF selected: its energy cost is burn × the price entered', N.kv(sustain.run({ ...d, fuel: 'SAF (HEFA-SPK)' })).energy_cost_usd_fh, d.fuel_kg_fh * d.fuel_usd_kg, 1e-12, 'The price entered is that of the selected fuel'),
      N.check('Kerosene CO₂ = 3.16 kg per kg of fuel', o.co2_kg_fh, 3.16 * d.fuel_kg_fh, 1e-12, 'Combustion emission factor'),
      N.check('Carbon bill present value: growing annuity', o.carbon_cost_life_pv_usd, (() => { const a = (3.16 * d.fuel_kg_fh * d.util_fh_yr * d.carbon_usd_t * d.carbon_cover) / 1e3, q = (1 + d.carbon_escal) / (1 + d.discount); return (a * q * (1 - q ** d.life_yr)) / (1 - q); })(), 1e-10, 'Geometric series'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], t = res.tables[0].rows, eff = i.carbon_usd_t * i.carbon_cover;
    const best = res.kpis.find((k) => k.key === 'cheapest_abatement_usd_t').label, dropIn = /SAF/.test(best), same = best.endsWith(o.current_carrier), pays = o.cheapest_abatement_usd_t <= eff;
    if (Number.isFinite(o.cheapest_abatement_usd_t) && same) out.push({ severity: 'info', title: 'The current energy carrier is already the cheapest low-carbon option', detail: `${o.current_carrier}: ${o.cheapest_abatement_usd_t.toFixed(0)} USD per tonne of CO₂ avoided against fossil fuel on the energy bill.`, action: 'Keep the energy source low-carbon (renewable electricity or hydrogen supply contracts): that is where the remaining emissions are.', basis: 'Marginal abatement cost against the carbon price' });
    else if (Number.isFinite(o.cheapest_abatement_usd_t) && pays && !dropIn) out.push({ severity: 'info', title: 'A different energy carrier would cut the energy and carbon bill', detail: best + `; effective carbon price today ${eff.toFixed(0)} USD/t.`, action: 'This is the energy bill only. Hydrogen and battery options need a different aircraft, whose mission feasibility (Suites 23 and 24), capital and maintenance cost are not assessed here: treat it as a reason to study that aircraft, not as a switch you can make.', basis: 'Marginal abatement cost against the carbon price' });
    else if (Number.isFinite(o.cheapest_abatement_usd_t)) out.push({ severity: pays ? 'advise' : 'info', title: pays ? 'A lower-carbon fuel already pays for itself' : `Cheapest abatement costs ${o.cheapest_abatement_usd_t.toFixed(0)} USD per tonne`, detail: best + `; effective carbon price today ${eff.toFixed(0)} USD/t.`, action: pays ? 'Switch as supply allows; lock in offtake contracts before the carbon price rises further.' : `Efficiency first: every kilogram of fuel not burned avoids 3.16 kg of CO₂ at negative cost. This ranking is the energy bill only: hydrogen and battery options need a different aircraft whose feasibility and capital cost are not assessed here. Alternative carriers break even when the carbon price passes ${o.cheapest_abatement_usd_t.toFixed(0)} USD/t or their price premium falls.`, basis: 'Marginal abatement cost against the carbon price' });
    if (o.carbon_cost_life_pv_usd > 0) out.push({ severity: 'info', title: `Carbon bill over the life: ${fmt(o.carbon_cost_life_pv_usd)} USD in present value`, detail: `Rising to ${fmt(o.carbon_cost_last_year_usd)} USD in the final year at ${(100 * i.carbon_escal).toFixed(1)}% escalation.`, action: 'Include this in fleet decisions: it is the budget available for efficiency retrofits, SAF contracts or earlier replacement.', basis: 'Discounted carbon-cost trajectory' });
    const el = t[4]; if (el && el[3] > t[0][3]) out.push({ severity: 'advise', title: 'Battery-electric would emit more on this grid', detail: `${el[3].toFixed(0)} against ${t[0][3].toFixed(0)} kg CO₂ per flight hour at ${i.grid_g_kWh} g/kWh.`, action: 'The climate case for electric flight depends on low-carbon electricity; contract renewable supply for charging.', basis: 'Grid carbon intensity' });
    return out;
  },
};

export default {
  id: 'economics', n: 26,
  tagline: 'What the aircraft costs to buy, fly and keep, whether the investment pays, how risky it is, and which technical or commercial lever is worth the most.',
  analyses: [opcost, invest, lcc, finance, risk, route, sustain],
  consumes: [
    { from: 'mission', keys: ['block_fuel_kg', 'block_time_h', 'mission_energy_kWh'], why: 'Fuel or energy per flight hour and block speed' },
    { from: 'fatigue', keys: ['life_fh'], why: 'Structural life limits the service life' },
    { from: 'safety', keys: ['availability', 'dispatch_reliability'], why: 'Availability scales utilisation, revenue and fleet size' },
    { from: 'propulsion', keys: ['eta_overall'], why: 'Powerplant efficiency for the energy-carrier comparison' },
  ],
  provides: [
    { key: 'doc_usd_fh', label: 'Direct operating cost', unit: 'USD/FH' }, { key: 'cask_usd', label: 'Cost per available seat-kilometre', unit: 'USD' }, { key: 'npv_usd', label: 'Net present value', unit: 'USD' }, { key: 'irr', label: 'Internal rate of return', unit: '-' },
    { key: 'payback_yr', label: 'Payback period', unit: 'yr' }, { key: 'lcc_usd', label: 'Life-cycle cost', unit: 'USD' }, { key: 'co2_cost_usd_fh', label: 'Carbon cost', unit: 'USD/FH' }, { key: 'breakeven_load_factor', label: 'Break-even load factor', unit: '-' },
  ],
  handoff: [
    { model: 'Published parametric cost-estimating relationships (DAPCA-IV, PRICE, SEER and similar)', why: 'Their coefficients, labour rates and dollar-year factors must come from the licensed or published source; a transparent user-set power law is provided instead', tool: 'Parametric cost-estimating tools with sourced databases' },
    { model: 'Airline revenue management, network and fleet-assignment optimisation', why: 'Needs booking data, fare classes, connecting flows and integer programming over the whole network', tool: 'Revenue-management and network-planning systems' },
    { model: 'Detailed maintenance cost modelling (task-level MSG-3 programme, shop-visit forecasting)', why: 'Requires the maintenance planning document, part prices and removal statistics', tool: 'Maintenance cost and engine shop-visit forecasting tools, fed by Suite 22' },
    { model: 'Fuel and carbon price forecasting', why: 'Market forecasting is outside an engineering model; prices are user inputs refreshed from live data, with volatility set by the user', tool: 'Commodity market analysis' },
    { model: 'Tax, accounting and lease-contract detail', why: 'Jurisdiction-specific rules (allowances, lease accounting standards, return conditions) are simplified', tool: 'Corporate finance and tax advice' },
    { model: 'Multi-stage and compound real options', why: 'Only a single deferral option on a binomial lattice is valued', tool: 'Real-options analysis software' },
  ],
};
