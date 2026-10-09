// Suite 22 — Aircraft Safety, Reliability and Failure Analysis.
// Fault trees (minimal cut sets, exact top-event probability, importance, common cause), FMECA and risk matrix,
// reliability block diagrams with exponential/Weibull life, Markov availability and dispatch, event trees,
// Bayesian failure-rate updating, Weibull fitting and the optimal preventive-replacement interval.
// Default probability objectives are the AC 25.1309-1B / AC 23.1309-1E values, the turbine in-flight-shutdown rate is the
// AC 120-42B diversion-time threshold and the common-cause β is the NUREG/CR-5485 generic screening value. Every other
// failure rate, conditional probability and cost in this file is an ILLUSTRATIVE placeholder of generic order of
// magnitude. Replace them with sourced fleet, supplier or handbook data before drawing conclusions.
import * as N from '../core/numerics.js';

// ---- probability helpers --------------------------------------------------------------------
const pop = (m) => { let c = 0; for (; m; m &= m - 1) c++; return c; };
const binom = (n, k) => { let c = 1; for (let j = 1; j <= k; j++) c = (c * (n - k + j)) / j; return c; };
function combos(n, k) { const out = []; const rec = (s, c) => { if (c.length === k) return out.push(c.slice()); for (let j = s; j < n; j++) { c.push(j); rec(j + 1, c); c.pop(); } }; rec(0, []); return out; }
/** Reliability of a k-out-of-n block of identical units each with reliability R. */
const kofn = (k, n, R) => { let s = 0; for (let j = k; j <= n; j++) s += binom(n, j) * R ** j * (1 - R) ** (n - j); return Math.min(1, s); };
/** Weibull survival conditional on survival to `age`: R(t | age) = R(t + age) / R(age). shape = 1 gives the exponential law. */
const weibR = (t, shape, scale, age = 0) => Math.exp(-(((t + age) / scale) ** shape) + (age / scale) ** shape);
const weibScale = (mttf, shape) => mttf / N.gamma(1 + 1 / shape);
/** Regularised lower incomplete gamma function P(a, x) (series / continued fraction). */
function gammP(a, x) {
  if (x <= 0) return 0;
  const gl = N.gammaln(a), pre = Math.exp(-x + a * Math.log(x) - gl);
  if (x < a + 1) { let ap = a, d = 1 / a, s = d; for (let n = 0; n < 800; n++) { ap++; d *= x / ap; s += d; if (Math.abs(d) < Math.abs(s) * 1e-15) break; } return s * pre; }
  let b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
  for (let n = 1; n < 800; n++) { const an = -n * (n - a); b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300; c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300; d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < 1e-15) break; }
  return 1 - pre * h;
}
/** Quantile of a Gamma(shape a, rate b) variable by bisection. */
function gammaQuantile(p, a, b) { let lo = 0, hi = a + 12 * Math.sqrt(a) + 60; for (let j = 0; j < 200; j++) { const m = 0.5 * (lo + hi); if (gammP(a, m) < p) lo = m; else hi = m; } return (0.5 * (lo + hi)) / b; }
const gammaPdf = (x, a, b) => (x <= 0 ? 0 : Math.exp(a * Math.log(b) + (a - 1) * Math.log(x) - b * x - N.gammaln(a)));

// ---- fault-tree engine ----------------------------------------------------------------------
// tree = { top, nodes: { name: { gate: 'AND'|'OR'|'KOFN', k, kids: [names] } | { lam [1/h], t [h] } | { p } } }
const minimise = (S) => { const u = [...new Set(S)].sort((a, b) => pop(a) - pop(b)), out = []; for (const s of u) if (!out.some((m) => (m & s) === m)) out.push(s); return out; };
/** Minimal cut sets by top-down gate expansion (MOCUS) with Boolean absorption; sets are bit masks over the basic events. */
function cutSets(tree) {
  const events = Object.keys(tree.nodes).filter((k) => !tree.nodes[k].gate), idx = Object.fromEntries(events.map((k, j) => [k, j])), memo = {};
  if (events.length > 30) throw new Error('more than 30 basic events');
  const and = (lists) => lists.reduce((acc, L) => minimise(acc.flatMap((a) => L.map((b) => a | b))), [0]);
  const expand = (name, depth) => {
    if (memo[name]) return memo[name];
    const nd = tree.nodes[name];
    if (!nd) throw new Error(`undefined node "${name}"`);
    if (depth > 30) throw new Error('circular gate definition');
    if (!nd.gate) return (memo[name] = [1 << idx[name]]);
    const ks = nd.kids.map((k) => expand(k, depth + 1));
    if (nd.gate === 'OR') return (memo[name] = minimise(ks.flat()));
    if (nd.gate === 'AND') return (memo[name] = and(ks));
    if (!(nd.k >= 1 && nd.k <= ks.length)) throw new Error(`k-of-n gate "${name}" needs 1 ≤ k ≤ n`);
    return (memo[name] = minimise(combos(ks.length, nd.k).flatMap((c) => and(c.map((j) => ks[j])))));
  };
  return { events, sets: expand(tree.top, 0) };
}
/** Exact probability of the union of cut sets for independent basic events (Shannon decomposition, i.e. an implicit BDD). */
function probUnion(sets, p) {
  const rec = (S) => {
    if (!S.length) return 0;
    if (S[0] === 0) return 1;
    const bit = S[0] & -S[0], j = 31 - Math.clz32(bit);
    return p[j] * rec(minimise(S.map((s) => s & ~bit))) + (1 - p[j]) * rec(S.filter((s) => !(s & bit)));
  };
  return rec(minimise(sets));
}
const eventP = (nd) => (nd.p !== undefined ? N.clamp(nd.p, 0, 1) : -Math.expm1(-Math.max(0, nd.lam) * Math.max(0, nd.t)));
function solveTree(tree, full = true) {
  const { events, sets } = cutSets(tree), p = events.map((e) => eventP(tree.nodes[e]));
  const pc = sets.map((s) => events.reduce((q, _, j) => (s & (1 << j) ? q * p[j] : q), 1));
  const pTop = probUnion(sets, p), rare = N.sum(pc), mcub = 1 - pc.reduce((q, v) => q * (1 - v), 1);
  if (!full) return { pTop };
  const imp = events.map((name, j) => {
    const b = 1 << j, p1 = p.slice(), p0 = p.slice(); p1[j] = 1; p0[j] = 0;
    return { name, p: p[j], fv: pTop > 0 ? probUnion(sets.filter((s) => s & b), p) / pTop : 0, birnbaum: probUnion(sets, p1) - probUnion(sets, p0), ccf: !!tree.nodes[name].ccf };
  });
  const ccfMask = events.reduce((m, e, j) => (tree.nodes[e].ccf ? m | (1 << j) : m), 0);
  const cuts = sets.map((s, c) => ({ names: events.filter((_, j) => s & (1 << j)), order: pop(s), p: pc[c] })).sort((a, b) => b.p - a.p);
  const single = sets.filter((s) => pop(s) === 1 && !(s & ccfMask)).map((s) => ({ name: events[31 - Math.clz32(s)], p: p[31 - Math.clz32(s)] }));
  // top event with every common-cause event removed: what the independent failures alone give
  const pInd = probUnion(sets.filter((s) => !(s & ccfMask)), p);
  return { events, p, sets, pTop, pInd, rare, mcub, imp, cuts, single, ccfShare: pTop > 0 ? probUnion(sets.filter((s) => s & ccfMask), p) / pTop : 0, minOrder: sets.length ? Math.min(...sets.map(pop)) : 0 };
}
/** Parse "TOP = OR(A, B); A = AND(e1, e2); B = 2oo3(x, y, z); e1 = 1e-5; e2 = 2e-4 @ 500; x = p:0.01". First name is the top event. */function parseTree(text, tDefault) {
  const nodes = {}; let top = null;
  for (const raw of String(text).split(/[;\n]+/)) {
    const line = raw.trim(); if (!line) continue;
    const m = line.match(/^([\w.-]+)\s*=\s*(.+)$/); if (!m) throw new Error(`cannot read "${line}"`);
    const name = m[1], rhs = m[2].trim(), g = rhs.match(/^(AND|OR|(\d+)oo(\d+))\s*\((.*)\)$/i);
    if (g) {
      const kids = g[4].split(',').map((s) => s.trim()).filter(Boolean); if (!kids.length) throw new Error(`gate "${name}" has no inputs`);
      nodes[name] = g[2] ? { gate: 'KOFN', k: +g[2], kids } : { gate: g[1].toUpperCase(), kids };
      if (g[2] && +g[3] !== kids.length) throw new Error(`gate "${name}" declares ${g[3]} inputs but lists ${kids.length}`);
    } else if (/^p\s*:/i.test(rhs)) { const v = Number(rhs.replace(/^p\s*:/i, '')); if (!(v >= 0 && v <= 1)) throw new Error(`probability of "${name}" must be 0…1`); nodes[name] = { p: v }; }
    else { const [l, t] = rhs.split('@').map((s) => Number(s)); if (!(l >= 0)) throw new Error(`failure rate of "${name}" is not a number`); nodes[name] = { lam: l, t: Number.isFinite(t) ? t : tDefault }; }
    top = top || name;
  }
  if (!top || !nodes[top].gate) throw new Error('the first line must define the top gate');
  return { top, nodes };
}
/** Default tree from the vehicle architecture: loss of propulsive power, primary flight control or electrical power. */
function defaultTree(i) {
  const n = {}, T = i.t_flight, tops = [];
  // A function is lost when fewer than `need` of `cnt` channels remain, i.e. after m = cnt − need + 1 failures. One
  // common-cause event takes those m channels out together at rate λ·β^(m−1) (multiple-Greek-letter model, equal factors).
  const group = (tag, label, cnt, need, lam) => {
    cnt = Math.round(cnt); need = N.clamp(Math.round(need), 1, cnt);
    const kids = N.range(cnt, (j) => `${tag}${j + 1}`), m = cnt - need + 1;
    kids.forEach((k) => (n[k] = { lam: lam * (m > 1 ? 1 - i.beta : 1), t: T }));
    if (m > 1) { n[`${tag}_CCF`] = { lam: lam * i.beta ** (m - 1), t: T, ccf: true }; n[`${tag}_INDEP`] = { gate: 'KOFN', k: m, kids }; n[label] = { gate: 'OR', kids: [`${tag}_INDEP`, `${tag}_CCF`] }; }
    else n[label] = { gate: 'OR', kids };
  };
  group('ENG', 'POWER_LOSS', Math.max(1, i.n_eng), i.eng_need, i.lam_engine);
  n.NO_SAFE_LANDING = { p: i.p_forced }; n.THRUST_CAT = { gate: 'AND', kids: ['POWER_LOSS', 'NO_SAFE_LANDING'] }; tops.push('THRUST_CAT');
  const ctl = [];
  if (i.n_hyd >= 1 && !i.manual_rev) { group('HYD', 'HYD_LOSS', i.n_hyd, 1, i.lam_hyd); ctl.push('HYD_LOSS'); }
  if (i.n_fcc >= 1) { group('FCC', 'FCS_LOSS', i.n_fcc, 1, i.lam_fcc); ctl.push('FCS_LOSS'); }
  n.CTRL_JAM = { lam: i.lam_jam, t: T }; ctl.push('CTRL_JAM');
  n.CONTROL_LOSS = { gate: 'OR', kids: ctl }; tops.push('CONTROL_LOSS');
  if (i.n_gen >= 1) { group('GEN', 'GEN_LOSS', i.n_gen, 1, i.lam_gen); n.BATT_LATENT = { lam: i.lam_batt, t: i.t_latent }; n.ELEC_LOSS = { gate: 'AND', kids: ['GEN_LOSS', 'BATT_LATENT'] }; }
  else { const nb = Math.max(1, Math.round(i.n_batt)); group('BAT', 'ELEC_LOSS', nb, nb > 1 ? nb - 1 : 1, i.lam_batt); } // battery-only aircraft: flight continues with one pack lost when more than one is fitted
  tops.push('ELEC_LOSS');
  // nobody on board: losing the aircraft is catastrophic only if it also strikes a person
  if (i.unmanned) { n.AIRCRAFT_LOST = { gate: 'OR', kids: tops }; n.PERSON_STRUCK = { p: i.p_third }; n.CATASTROPHIC = { gate: 'AND', kids: ['AIRCRAFT_LOST', 'PERSON_STRUCK'] }; }
  else n.CATASTROPHIC = { gate: 'OR', kids: tops };
  return { top: 'CATASTROPHIC', nodes: n };
}

// ---- shared inputs --------------------------------------------------------------------------
const ARCH = [
  { key: 'n_eng', label: 'Engines / motors', unit: '', default: 2, min: 1, max: 12, step: 1, discrete: true, group: 'Architecture' },
  { key: 'eng_need', label: 'Engines needed to continue safe flight', unit: '', default: 1, min: 1, max: 12, step: 1, discrete: true, group: 'Architecture', help: '1 for a transport twin; all of them for a quadrotor without motor-out capability' },
  { key: 'n_hyd', label: 'Independent hydraulic systems', unit: '', default: 3, min: 0, max: 4, step: 1, discrete: true, group: 'Architecture', help: '0 for manual or all-electric controls' },
  { key: 'n_gen', label: 'Generators', unit: '', default: 2, min: 0, max: 6, step: 1, discrete: true, group: 'Architecture' },
  { key: 'n_fcc', label: 'Flight-control computer channels', unit: '', default: 3, min: 0, max: 5, step: 1, discrete: true, group: 'Architecture', help: 'Fly-by-wire channels whose total loss means loss of control; 0 for mechanically signalled controls' },
];
const OCCUPANCY = [
  { key: 'n_batt', label: 'Independent battery packs (battery-only aircraft)', unit: '', default: 1, min: 1, max: 12, step: 1, discrete: true, group: 'Architecture', help: 'Used when there is no generator. With more than one pack the flight continues after losing one.' },
  { key: 'unmanned', label: 'Unmanned aircraft (nobody on board)', type: 'bool', default: false, group: 'Architecture', help: 'Loss of the aircraft is then not catastrophic by itself; effects are taken one severity class lower and the fault tree adds the chance of striking a person.' },
];
const RATES = [
  { key: 'lam_engine', label: 'Engine in-flight shutdown rate', unit: '1/FH', default: 3e-5, min: 0, max: 1, group: 'Failure rates', help: 'Per engine. The turbine default is the AC 120-42B in-flight-shutdown threshold for diversion times above 120 and up to 180 minutes (5e-5 up to 120 minutes, 2e-5 beyond 180 minutes): a regulatory ceiling, not a measured rate. The piston (1e-4) and electric-motor (2e-5) values are unsourced placeholders. Replace with fleet or supplier data.' },
  { key: 'lam_hyd', label: 'Hydraulic system loss rate', unit: '1/FH', default: 1e-4, min: 0, max: 1, group: 'Failure rates (illustrative)' },
  { key: 'lam_gen', label: 'Generator channel loss rate', unit: '1/FH', default: 2e-4, min: 0, max: 1, group: 'Failure rates (illustrative)' },
  { key: 'lam_fcc', label: 'Flight-control channel loss rate', unit: '1/FH', default: 1e-4, min: 0, max: 1, group: 'Failure rates (illustrative)' },
];
const OBJ_HELP = 'AC 25.1309-1B for large aeroplanes: 1e-9 catastrophic, 1e-7 hazardous, 1e-5 major, 1e-3 minor per flight hour ("on the order of"). AC 23.1309-1E Figure 2 relaxes them for small aeroplanes, down to 1e-6 / 1e-5 / 1e-4 / 1e-3 for a single piston engine up to 6000 lb (Class I). Filled from the vehicle class; set the values of your own certification basis.';
const OBJECTIVES = [
  { key: 'obj_cat', label: 'Objective: catastrophic', unit: '1/FH', default: 1e-9, min: 1e-12, max: 1, group: 'Safety objectives', help: OBJ_HELP },
  { key: 'obj_haz', label: 'Objective: hazardous', unit: '1/FH', default: 1e-7, min: 1e-12, max: 1, group: 'Safety objectives' },
  { key: 'obj_maj', label: 'Objective: major', unit: '1/FH', default: 1e-5, min: 1e-12, max: 1, group: 'Safety objectives' },
  { key: 'obj_min', label: 'Objective: minor', unit: '1/FH', default: 1e-3, min: 1e-12, max: 1, group: 'Safety objectives' },
];
// [catastrophic, hazardous, major, minor] per flight hour: AC 25.1309-1B (= AC 23.1309-1E Class IV) and AC 23.1309-1E Figure 2 Classes I–III
const OBJ_25 = [1e-9, 1e-7, 1e-5, 1e-3], OBJ_23 = { I: [1e-6, 1e-5, 1e-4, 1e-3], II: [1e-7, 1e-6, 1e-5, 1e-3], III: [1e-8, 1e-7, 1e-5, 1e-3] };
/** Objectives for the vehicle class. Rotorcraft and eVTOL keep the large-aeroplane values; unmanned aircraft get the Class I values as a placeholder (neither is sourced for those classes). */
const objDefaults = (c) => {
  const lb = c.mass.mtow_kg * 2.20462, single = c.prop.n_eng === 1, recip = c.prop.type === 'piston' || c.prop.type === 'electric';
  const o = c.meta.type === 'uav' ? OBJ_23.I : c.meta.type !== 'aeroplane' || lb > 12500 ? OBJ_25 : lb > 6000 ? OBJ_23.III : single && recip ? OBJ_23.I : OBJ_23.II;
  return { obj_cat: o[0], obj_haz: o[1], obj_maj: o[2], obj_min: o[3] };
};
const schemeOf = (cat) => { const m = (v) => Math.abs(Math.log10(cat / v)) < 1e-6; return m(1e-9) ? 'AC 25.1309-1B large-aeroplane value (AC 23.1309-1E Class IV is the same)' : m(1e-8) ? 'AC 23.1309-1E Class III value' : m(1e-7) ? 'AC 23.1309-1E Class II value' : m(1e-6) ? 'AC 23.1309-1E Class I value' : 'user-defined value'; };
/** The advisory material words every probability band as an order of magnitude: half a decade above an objective is still "of that order". */
const ORDER = Math.sqrt(10);
const flightHours = (c) => { const V = c.mission.cruise_V_ms || c.flight.V_ms || 0, fromRange = V > 0 && c.mission.range_km > 0 ? (c.mission.range_km * 1000) / V / 3600 : 0, perCycle = c.econ.cycles_yr > 0 ? c.econ.util_fh_yr / c.econ.cycles_yr : 0; return Math.max(0.1, perCycle || fromRange || 1); };
const archDefaults = (c) => {
  const elec = c.prop.type === 'electric', rotary = c.wing.S_m2 === 0 && c.rotor.R_m > 0, small = c.mass.mtow_kg < 5700, n = c.prop.n_eng, heli = c.meta.type === 'helicopter';
  // fly-by-wire is assumed for unmanned aircraft, eVTOL and large transports; other classes are taken as mechanically signalled
  const fbw = c.meta.type === 'uav' || c.meta.type === 'evtol' || (c.meta.type === 'aeroplane' && c.mass.mtow_kg > 40000);
  return {
    n_eng: n, eng_need: elec && n >= 6 ? n - 1 : elec && rotary ? n : 1,
    n_hyd: c.systems.hyd_p_Pa > 0 && !small ? (c.mass.mtow_kg > (heli ? 8000 : 40000) ? 3 : 2) : 0, n_gen: c.systems.gen_kVA > 0 ? Math.max(1, Math.min(n, 2)) : 0,
    n_fcc: c.meta.type === 'uav' ? 1 : fbw ? 3 : 0,
    lam_engine: c.prop.type === 'piston' ? 1e-4 : elec ? 2e-5 : 3e-5,
  };
};
const occDefaults = (c) => ({ n_batt: c.prop.type === 'electric' && !(c.systems.gen_kVA > 0) ? (c.meta.type === 'evtol' ? 4 : c.mass.mtow_kg > 600 ? 2 : 1) : 1, unmanned: c.meta.type === 'uav' });
const SEV = ['Catastrophic', 'Hazardous', 'Major', 'Minor'];
const objOf = (i) => [i.obj_cat, i.obj_haz, i.obj_maj, i.obj_min];
const ILLUSTRATIVE = 'Failure rates and conditional probabilities are illustrative generic placeholders, not sourced reliability data (the turbine in-flight-shutdown default is a regulatory threshold and the common-cause factor a generic screening value)';
// ---- analyses -------------------------------------------------------------------------------
const fta = {
  id: 'fta', title: 'Fault tree analysis of the catastrophic top event', fidelity: 'analytical',
  summary: 'Builds a fault tree for loss of power, flight control and electrical supply from the vehicle architecture (or reads your own tree), finds the minimal cut sets and the exact top-event probability, separates what independent failures give from what common cause adds, and ranks every basic event by importance.',
  equations: ['Fault tree probability equations', 'Probability theory', 'Conditional probability equations', 'Total probability theorem', 'Exponential reliability distribution'],
  inputs: [...ARCH, ...OCCUPANCY,
    { key: 'manual_rev', label: 'Flight controls work without hydraulic power (manual reversion)', type: 'bool', default: false, group: 'Architecture', help: 'When set, loss of every hydraulic system is not counted as loss of control' },
    { key: 't_flight', label: 'Average flight duration (exposure time)', unit: 'h', default: 2, min: 0.05, max: 24, group: 'Exposure' },
    { key: 't_latent', label: 'Latent-failure check interval', unit: 'h', default: 50, min: 1, max: 20000, group: 'Exposure', help: 'Exposure time of the standby battery, whose failure stays hidden until the next check. 14 CFR 25.1309(b)(5)(iii) limits the probability of latent failures combined with an active failure to 1/1000: with the default standby failure rate of 2e-5 per hour that allows about 50 h.' },
    ...RATES,
    { key: 'lam_batt', label: 'Battery / standby source failure rate', unit: '1/FH', default: 2e-5, min: 0, max: 1, group: 'Failure rates (illustrative)', help: 'Per pack for a battery-only aircraft' },
    { key: 'lam_jam', label: 'Single-point control jam or disconnect rate', unit: '1/FH', default: 1e-10, min: 0, max: 1, group: 'Failure rates (illustrative)', help: 'A single failure with a catastrophic effect has to be excluded by design; the placeholder is a tenth of the catastrophic objective' },
    { key: 'beta', label: 'Common-cause β-factor', unit: '-', default: 0.05, min: 0, max: 0.5, group: 'Failure rates', help: 'Share of a channel\'s failures that also take out the next redundant channel; m channels are lost together at rate λ·β^(m−1). The default 0.05 is the NUREG/CR-5485 generic screening value for a redundant pair tested at staggered times (0.10 when tested together). It is a nuclear-industry screening number, deliberately conservative: replace it only with a value backed by a common-cause analysis of your own channels.' },
    { key: 'p_forced', label: 'P(aircraft lost | total power loss)', unit: '-', default: 0.2, min: 0, max: 1, group: 'Failure rates (illustrative)', help: 'Chance that the forced landing, ditching or autorotation after total power loss is not survivable' },
    { key: 'p_third', label: 'P(person fatally struck | unmanned aircraft lost)', unit: '-', default: 0.01, min: 0, max: 1, group: 'Failure rates (illustrative)', help: 'Used only for unmanned aircraft. Depends on the population overflown and the size of the aircraft; take it from your operational risk assessment.' },
    ...OBJECTIVES.slice(0, 1),
    { key: 'tree_text', label: 'Custom fault tree (optional)', type: 'text', default: '', group: 'Custom tree', help: 'e.g. TOP = OR(A, C); A = AND(e1, e2); C = 2oo3(x, y, z); e1 = 1e-5; e2 = 2e-4 @ 500; x = p:0.01 — rates per hour, "@" sets an exposure time in hours, "p:" a fixed probability. Leave empty to use the tree built from the architecture.' },
  ],
  defaults: (c) => {
    const a = archDefaults(c), o = objDefaults(c), rotaryElec = c.wing.S_m2 === 0 && c.prop.type === 'electric';
    return { ...a, ...occDefaults(c), manual_rev: c.meta.type === 'aeroplane' && a.n_hyd >= 1 && c.mass.mtow_kg <= 40000, t_flight: flightHours(c), obj_cat: o.obj_cat, lam_jam: 0.1 * o.obj_cat,
      // light aeroplanes are stall-speed limited so that most forced landings are survivable
      p_forced: rotaryElec ? 0.9 : c.meta.type === 'uav' ? 0.5 : c.meta.type === 'helicopter' ? 0.15 : c.meta.type === 'evtol' ? 0.1 : c.mass.mtow_kg < 5700 ? 0.02 : 0.2 };
  },
  run(i) {
    const warnings = []; let tree = null, custom = false;
    if (String(i.tree_text || '').trim()) { try { tree = parseTree(i.tree_text, i.t_flight); solveTree(tree); custom = true; } catch (e) { tree = null; warnings.push(`The custom tree could not be used (${e.message}); the architecture-based default tree was solved instead.`); } }
    tree = tree || defaultTree(i);
    const r = solveTree(tree), perFh = r.pTop / i.t_flight, indFh = r.pInd / i.t_flight, imp = r.imp.slice().sort((a, b) => b.fv - a.fv), top = imp.slice(0, 10), obj = i.obj_cat, e0 = (v) => v.toExponential(0), e2 = (v) => v.toExponential(2);
    const Ts = N.logspace(0.25, 16, 13), pT = Ts.map((T) => { const t = custom ? parseTree(i.tree_text, T) : defaultTree({ ...i, t_flight: T }); return solveTree(t, false).pTop / T; });
    const spf = r.single.filter((e) => e.p / i.t_flight > obj), lat = !custom && tree.nodes.BATT_LATENT ? eventP(tree.nodes.BATT_LATENT) : NaN;
    // independent failures decide the verdict; a total that fails only through the screening β is a flag for common-cause analysis
    const indStatus = indFh <= obj ? 'ok' : indFh <= 100 * obj ? 'warn' : 'bad', totStatus = perFh <= obj ? 'ok' : indStatus === 'bad' ? 'bad' : 'warn';
    if (spf.length) warnings.push(`Single failure(s) lead directly to the top event with a probability above the objective: ${spf.map((e) => e.name).join(', ')}. Fail-safe design practice does not accept a single failure with a catastrophic effect.`);
    if (perFh > obj && indFh <= obj) warnings.push(`Independent failures meet the objective (${e2(indFh)} per FH) but the common-cause events lift the total to ${e2(perFh)} per FH. With a generic screening β this flags the redundant groups for a common-cause analysis; it does not show that the design is unsafe.`);
    else if (r.ccfShare > 0.5) warnings.push('The result is driven by the common-cause assumption (β-factor), not by independent failures: it is only as good as the β value and the independence of the channels.');
    if (r.rare > 0.1) warnings.push('Cut-set probabilities are not small: the rare-event approximation is poor here; the exact value is reported.');
    if (lat > 1e-3 * (1 + 1e-9)) warnings.push(`The standby source has failed unnoticed with probability ${e2(lat)} by the end of the ${i.t_latent} h check interval. For large aeroplanes 14 CFR 25.1309(b)(5)(iii) limits latent failures combined with an active failure to 1/1000: shorten the check interval or monitor the source.`);
    if (i.unmanned && !custom) warnings.push(`Unmanned aircraft: the top event is a fatal injury on the ground, i.e. the aircraft is lost and a person is struck (illustrative probability ${i.p_third}). The objective is a placeholder, because the data registry holds no sourced quantitative objective for unmanned aircraft.`);
    warnings.push(`${ILLUSTRATIVE}. The objective of ${e0(obj)} per flight hour is the ${schemeOf(obj)}; rotorcraft, powered-lift and unmanned aircraft have their own certification bases.`);
    return {
      kpis: [
        { key: 'p_catastrophic_per_fh', label: 'Top-event probability per flight hour, including common-cause screening', value: perFh, unit: '1/FH', status: totStatus, note: `Objective ${e0(obj)} per FH (average over the flight); ${e2(indFh)} without common cause` },
        { key: 'p_independent_per_fh', label: 'Top-event probability from independent failures only', value: indFh, unit: '1/FH', status: indStatus, note: 'Common-cause (β-factor) events removed' },
        { key: 'p_top_per_flight', label: 'Top-event probability per flight (exact)', value: r.pTop, unit: '-' },
        { key: 'p_top_rare_event', label: 'Rare-event approximation Σ P(cut set)', value: r.rare, unit: '-', note: 'Upper bound' },
        { key: 'p_top_mcub', label: 'Minimal-cut-set upper bound', value: r.mcub, unit: '-' },
        { key: 'n_cut_sets', label: 'Minimal cut sets', value: r.sets.length, unit: '' },
        { key: 'min_cut_order', label: 'Smallest cut-set order', value: r.minOrder, unit: '', note: 'Order 1 = one event (single failure or common cause) reaches the top' },
        { key: 'n_single_point_above_obj', label: 'Single failures above the objective', value: spf.length, unit: '', status: spf.length ? 'bad' : 'ok' },
        { key: 'ccf_share_pct', label: 'Common-cause share of the top event', value: 100 * r.ccfShare, unit: '%', status: r.ccfShare > 0.5 ? 'warn' : 'ok' },
        { key: 'objective_margin_decades', label: 'Margin to the objective', value: Math.log10(obj / Math.max(perFh, 1e-300)), unit: 'decades', status: totStatus, note: `${Math.log10(obj / Math.max(indFh, 1e-300)).toFixed(2)} decades without common cause` },
        ...(Number.isFinite(lat) ? [{ key: 'p_latent_standby', label: 'Probability that the standby source has failed unnoticed', value: lat, unit: '-', status: lat <= 1e-3 * (1 + 1e-9) ? 'ok' : 'warn', note: 'At the end of the check interval; limit 1/1000 for large aeroplanes' }] : []),
      ],
      plots: [
        { type: 'bar', title: 'Fussell–Vesely importance of basic events', ylabel: 'Share of top-event probability [-]', categories: top.map((e) => e.name), series: [{ name: 'Fussell–Vesely', y: top.map((e) => e.fv) }] },
        { type: 'line', title: 'Top-event probability per flight hour versus flight duration', xlabel: 'Flight duration [h]', ylabel: 'Probability per flight hour [1/FH]', xlog: true, ylog: true, series: [{ name: 'Top event', x: Ts, y: pT.map((v) => Math.max(v, 1e-300)) }], annotations: [{ y: obj, label: 'Objective' }, { x: i.t_flight, label: 'This flight' }] },
      ],
      tables: [
        { title: 'Minimal cut sets (most probable first)', columns: ['Events', 'Order', 'Probability per flight', 'Share [%]'], rows: r.cuts.slice(0, 25).map((c) => [c.names.join(' · '), c.order, c.p, r.rare > 0 ? (100 * c.p) / r.rare : 0]) },
        { title: 'Basic-event importance', columns: ['Event', 'Probability per flight', 'Fussell–Vesely', 'Birnbaum', 'Common cause'], rows: imp.map((e) => [e.name, e.p, e.fv, e.birnbaum, e.ccf ? 'yes' : '']) },
      ],
      outputs: { top_event: tree.top, custom_tree: custom ? 1 : 0 },
      warnings, models: ['Fault tree with AND / OR / k-of-n gates', 'MOCUS minimal cut sets with Boolean absorption', 'Exact top-event probability by Shannon decomposition', 'β-factor common-cause model (multiple-Greek-letter form with equal factors)', 'Fussell–Vesely and Birnbaum importance'],
      assumptions: ['Basic events are statistically independent apart from the explicit β-factor events', 'Constant failure rates: P = 1 − exp(−λ·t)', 'Probability per flight hour is the per-flight probability divided by the flight duration', 'No repair in flight; latent failures exposed for the full check interval', 'Total loss of electrical power and of powered or computed flight control is taken as catastrophic; with manual reversion, hydraulic loss is not', ILLUSTRATIVE],
    };
  },
  verify() {
    const t1 = solveTree(parseTree('TOP = OR(G, C); G = AND(A, B); A = p:0.1; B = p:0.2; C = p:0.05', 1));
    const t2 = solveTree(parseTree('TOP = 2oo3(A, B, C); A = p:0.1; B = p:0.1; C = p:0.1', 1));
    const t3 = solveTree(parseTree('TOP = AND(G1, G2); G1 = OR(A, B); G2 = OR(A, C); A = p:0.1; B = p:0.2; C = p:0.3', 1)); // shared event: A + B·C
    const b = Object.fromEntries(fta.inputs.map((f) => [f.key, f.default])), q = defaultTree({ ...b, n_eng: 4, eng_need: 3, lam_engine: 1e-4, beta: 0.1 }).nodes;
    const lone = { ...b, n_eng: 1, eng_need: 1, n_hyd: 0, n_gen: 0, n_fcc: 0, n_batt: 1, lam_engine: 1e-3, lam_batt: 0, lam_jam: 0, p_forced: 1, t_flight: 2 };
    const u = N.kv(fta.run({ ...lone, unmanned: true, p_third: 0.25 })), d = N.kv(fta.run(b)), g = solveTree(defaultTree({ ...lone, n_eng: 2, lam_engine: 0.05, beta: 0.2, t_flight: 1 }));
    return [
      N.check('OR(AND(A,B),C) exact probability', t1.pTop, 1 - (1 - 0.02) * (1 - 0.05), 1e-12, 'Hand solution 1 − (1 − pA·pB)(1 − pC)'),
      N.check('2-out-of-3 failure probability 3p² − 2p³', t2.pTop, 3 * 0.01 - 2 * 0.001, 1e-12, 'Binomial'),
      N.check('2-out-of-3 gives three minimal cut sets', t2.sets.length, 3, 0, 'Combinatorics'),
      N.check('Repeated event absorbed: (A+B)(A+C) = A + BC', t3.pTop, 0.1 + 0.9 * 0.06, 1e-12, 'Boolean absorption law'),
      N.check('Birnbaum importance of C in OR(AND(A,B),C)', t1.imp.find((e) => e.name === 'C').birnbaum, 1 - 0.02, 1e-12, '∂P(top)/∂pC'),
      N.check('Rare-event bound is conservative', t1.rare >= t1.pTop ? 1 : 0, 1, 0, 'Boole inequality'),
      N.check('3-of-4 group: two failures defeat it, so the common-cause rate is λ·β', q.ENG_CCF.lam, 1e-4 * 0.1, 1e-12, 'Multiple-Greek-letter model: m failures in common at λ·β^(m−1)'),
      N.check('3-of-4 group is a 2-out-of-4 failure gate', q.ENG_INDEP.k, 2, 0, 'm = n − k + 1'),
      N.check('Redundant pair: independent part is (1 − exp(−λ(1−β)t))²', g.pInd, (1 - Math.exp(-0.04)) ** 2, 1e-12, 'Common-cause event removed'),
      N.check('Redundant pair: total adds the common-cause event', g.pTop, 1 - (1 - (1 - Math.exp(-0.04)) ** 2) * Math.exp(-0.01), 1e-12, 'Union of independent pair failure and λβ event'),
      N.check('Unmanned: P(top) = P(aircraft lost) × P(person struck)', u.p_top_per_flight, 0.25 * (1 - Math.exp(-2e-3)), 1e-12, 'AND gate with an independent conditional event'),
      N.check('Default latent exposure sits at the 1/1000 limit', d.p_latent_standby, 1 - Math.exp(-2e-5 * 50), 1e-12, 'P = 1 − exp(−λ·t) with λ = 2e-5 /h, t = 50 h'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], imp = res.tables[1].rows[0], ind = o.p_independent_per_fh, tot = o.p_catastrophic_per_fh, obj = i.obj_cat, e2 = (v) => v.toExponential(2);
    const lead = `The largest contributor is ${imp[0]} (${(100 * imp[2]).toFixed(0)}% of the total).`, placeholder = i.unmanned ? ' The objective for an unmanned aircraft is a placeholder, so this is a prompt to set your own, not a finding of non-compliance.' : '';
    // without a sourced objective for unmanned aircraft an exceedance cannot be more than advice
    if (ind > ORDER * obj) out.push({ severity: i.unmanned ? 'advise' : ind > 100 * obj ? 'critical' : 'warn', title: 'Independent failures alone exceed the objective', detail: `${e2(ind)} per FH from independent failures (${e2(tot)} with common cause) against ${obj.toExponential(0)}. ${lead}${placeholder}`, action: `Work on ${imp[0]} first: add an independent channel, lower its failure rate, or shorten its exposure time. Re-run to confirm the gain before touching smaller contributors.`, basis: 'Fussell–Vesely importance ranking against the quantitative safety objective' });
    else if (ind > obj) out.push({ severity: 'advise', title: 'Top-event probability is of the same order as the objective', detail: `${e2(ind)} per FH from independent failures against ${obj.toExponential(0)}. ${lead}`, action: 'The objectives are stated as orders of magnitude and these inputs are placeholders: firm up the failure rate and exposure time of the leading contributor before concluding either way.', basis: 'Probability bands defined "on the order of" a value (AC 25.1309-1B; AC 23.1309-1E Figure 2, Note 1)' });
    else if (tot > obj) out.push({ severity: 'advise', title: 'Objective met by independent failures; common cause still to be shown', detail: `${e2(ind)} per FH independent, ${e2(tot)} with the β-factor events (β = ${i.beta}) against ${obj.toExponential(0)}; ${o.ccf_share_pct.toFixed(0)}% of the total is common cause.`, action: 'Treat the β result as a screening flag. Analyse each redundant group for common causes (segregation and zonal separation, dissimilar channels, independent power and software, maintenance that touches every channel in one visit), then replace β with a value you can justify.', basis: 'β-factor screening with generic values (NUREG/CR-5485); credit for redundancy requires demonstrated independence' });
    else if (o.ccf_share_pct > 50) out.push({ severity: 'advise', title: 'Redundancy is limited by common cause', detail: `${o.ccf_share_pct.toFixed(0)}% of the top event comes from β-factor events, so adding more identical channels buys little.`, action: 'Use dissimilar or segregated channels (different power source, routing, software), and review maintenance practices that could affect all channels in one visit.', basis: 'β-factor common-cause model' });
    if (o.n_single_point_above_obj > 0) out.push({ severity: 'warn', title: 'Single-point failure above the objective', detail: `${o.n_single_point_above_obj} single event(s) reach the top event on their own; see the warning list.`, action: 'Add redundancy or a design feature that removes the single failure, or show by design (not probability) that it cannot occur.', basis: 'Fail-safe design concept: no single failure may be catastrophic' });
    if (o.p_latent_standby > 1e-3 * (1 + 1e-9)) out.push({ severity: obj <= 1e-9 * (1 + 1e-6) ? 'warn' : 'advise', title: 'Latent standby failure above 1 in 1000', detail: `The standby source has failed unnoticed with probability ${e2(o.p_latent_standby)} at the end of the ${i.t_latent} h check interval.`, action: 'Shorten the check interval or add monitoring; halving the interval roughly halves the exposure and avoids carrying a dead battery (and its replacement mass) unknowingly.', basis: '14 CFR 25.1309(b)(5)(iii): latent failures combined with an active failure limited to 1/1000 (large aeroplanes)' });
    return out;
  },
};
// Failure modes that are not judged against a probability objective, with the reason shown in the worksheet
const BASIS = { oei: 'design case: flight continues with an engine out', forced: 'design case: forced landing after engine failure', dt: 'shown by fatigue and damage-tolerance evaluation', pr: 'particular risk: shown by design and test' };
const FMEA_ROWS = (i) => {
  const multi = i.n_eng > i.eng_need, r = [], nb = Math.max(1, Math.round(i.n_batt));
  // [item, failure mode, effect, severity index 0..3, item failure rate, mode ratio α, effect probability β, detection 1..10, basis when no probability objective applies]
  r.push(['Engine / motor', 'In-flight shutdown or loss of power', multi ? 'Continued flight on remaining engines, diversion' : 'Forced landing / autorotation', multi ? 2 : 1, i.lam_engine, 0.9, 1, 2, multi ? BASIS.oei : BASIS.forced]);
  r.push(['Engine / motor', 'Uncontained rotor burst', 'Debris damage to structure and systems', 0, i.lam_engine, 0.001, 0.05, 8, BASIS.pr]);
  r.push(['Fuel / energy supply', 'Loss of feed to one engine', multi ? 'Engine runs down; cross-feed, continued flight' : 'Loss of power, forced landing', multi ? 2 : 1, i.lam_engine, 0.05, 1, 5, multi ? BASIS.oei : BASIS.forced]);
  if (i.n_hyd >= 1) { r.push(['Hydraulic system', 'Loss of one system (leak, pump)', i.n_hyd > 1 ? 'Loss of redundancy, gear by alternate means' : 'Loss of powered controls', i.n_hyd > 1 ? 3 : 1, i.lam_hyd, 1, 1, 2, '']); r.push(['Flight-control actuator', 'Jam or runaway', 'Loss of control of one surface', i.n_hyd > 1 ? 2 : 1, i.lam_hyd * 0.02, 0.5, 0.3, 6, '']); }
  if (i.n_gen >= 1) { r.push(['Electrical generation', 'Loss of one generator channel', i.n_gen > 1 ? 'Load shedding' : 'Battery-only flight, land as soon as possible', i.n_gen > 1 ? 3 : 2, i.lam_gen, 1, 1, 2, '']); r.push(['Battery / standby power', 'Capacity loss (latent)', 'No effect until generation is lost', 3, i.lam_batt, 1, 1, 9, '']); }
  else r.push(['Battery pack', 'Loss of one pack or string', nb > 1 ? 'Reduced power and range, land as soon as practicable' : 'Loss of propulsion and control power', nb > 1 ? 2 : 0, i.lam_batt, 1, 1, 4, ''], ['Battery pack', 'Thermal runaway not contained', 'Fire, loss of the aircraft', 0, i.lam_batt, 0.01, 1, 6, BASIS.pr]);
  if (i.n_fcc >= 1) r.push(['Flight-control computer', 'Loss of one channel', i.n_fcc > 1 ? 'Reversion to remaining channels' : 'Loss of computed control', i.n_fcc > 1 ? 3 : 0, i.lam_fcc, 1, 1, 2, ''], ['Flight-control computer', 'Undetected erroneous output', 'Hard-over or oscillatory command', 0, i.lam_fcc, 0.01, i.n_fcc > 1 ? 1e-4 : 0.01, 9, '']);
  r.push(['Air data / navigation sensors', 'Erroneous output', 'Misleading information to crew or autopilot', 1, 5e-5, 0.1, 0.01, 7, '']);
  r.push(['Primary structure', 'Fatigue crack undetected to critical length', 'Structural failure', 0, i.lam_struct, 1, 1, 6, BASIS.dt]);
  if (i.has_rotor) r.push(['Rotor blade / hub', 'Fatigue crack growth to failure', 'Loss of the aircraft', 0, i.lam_struct, 1, 1, 6, BASIS.dt]);
  if (i.has_gearbox) r.push(['Main gearbox', 'Loss of lubrication or gear failure', 'Loss of drive to the rotor', 0, 2e-5, 0.05, 0.3, 5, BASIS.dt], ['Tail rotor / anti-torque drive', 'Loss of drive', 'Loss of yaw control, autorotative landing', 1, 1e-5, 0.3, 0.5, 5, BASIS.dt]);
  if (i.has_gear) r.push(['Landing gear', 'Failure to extend', 'Gear-up landing', 2, 2e-5, 0.2, 0.5, 2, ''], ['Wheel brakes / steering', 'Loss of braking on one side', 'Runway excursion risk', 2, 5e-5, 0.3, 0.1, 3, '']);
  // nobody on board: every effect is one class less severe
  return i.unmanned ? r.map((q) => { const s = q.slice(); s[3] = Math.min(3, s[3] + 1); return s; }) : r;
};
const sevIndex = (s) => { const v = String(s).trim().toLowerCase(), k = SEV.findIndex((n) => n.toLowerCase().startsWith(v.slice(0, 3))); return k >= 0 ? k : N.clamp(Math.round(Number(v)) - 1, 0, 3); };
const fmeca = {
  id: 'fmeca', title: 'FMECA table, hazard classification and risk matrix', fidelity: 'analytical',
  summary: 'Generates a failure modes, effects and criticality table for the systems present on the vehicle, classifies each effect by severity, compares its probability per flight hour with the objective for that class and places it on a risk matrix.',
  equations: ['Failure rate equations', 'Probability theory', 'Conditional probability equations'],
  inputs: [...ARCH, ...OCCUPANCY, ...RATES,
    { key: 'lam_batt', label: 'Battery / standby source failure rate', unit: '1/FH', default: 2e-5, min: 0, max: 1, group: 'Failure rates (illustrative)' },
    { key: 'lam_struct', label: 'Structural failure rate between inspections', unit: '1/FH', default: 1e-9, min: 0, max: 1, group: 'Failure rates (illustrative)', help: 'Shown for completeness only: structural integrity is demonstrated by the fatigue and damage-tolerance evaluation of Suite 9, not by this number' },
    { key: 't_flight', label: 'Average flight duration', unit: 'h', default: 2, min: 0.05, max: 24, group: 'Exposure' },
    { key: 'has_rotor', label: 'Lifting rotors', type: 'bool', default: false, group: 'Architecture' },
    { key: 'has_gearbox', label: 'Main gearbox and tail-rotor drive (conventional helicopter)', type: 'bool', default: false, group: 'Architecture' },
    { key: 'has_gear', label: 'Retractable gear / wheel brakes', type: 'bool', default: true, group: 'Architecture' },
    ...OBJECTIVES,
    { key: 'extra_rows', label: 'Additional failure modes', type: 'text', default: '', group: 'Custom rows', help: 'One per line: item; failure mode; severity (1 catastrophic … 4 minor, or the word); failure rate per FH; mode ratio α; effect probability β; detection 1–10' },
  ],
  defaults: (c) => ({ ...archDefaults(c), ...occDefaults(c), ...objDefaults(c), t_flight: flightHours(c), has_rotor: c.rotor.R_m > 0 && c.meta.type !== 'aeroplane', has_gearbox: c.meta.type === 'helicopter', has_gear: c.gear.type !== 'skid' }),
  run(i) {
    const warnings = [], rows = FMEA_ROWS(i), obj = objOf(i); let badLines = 0;
    for (const raw of String(i.extra_rows || '').split(/\n+/)) {
      const f = raw.split(';').map((s) => s.trim()); if (f.length < 2 && !f[0]) continue;
      const lam = Number(f[3]);
      if (f.length < 4 || !(lam >= 0)) { badLines++; continue; }
      rows.push([f[0], f[1], 'User-defined', sevIndex(f[2]), lam, Number.isFinite(+f[4]) && f[4] !== undefined && f[4] !== '' ? +f[4] : 1, Number.isFinite(+f[5]) && f[5] !== undefined && f[5] !== '' ? +f[5] : 1, N.clamp(Math.round(+f[6]) || 5, 1, 10), '']);
    }
    if (badLines) warnings.push(`${badLines} additional row(s) could not be read (need at least item; mode; severity; failure rate) and were skipped.`);
    const band = (p) => (p > 1e-3 ? 0 : p > 1e-5 ? 1 : p > 1e-7 ? 2 : p > 1e-9 ? 3 : 4), BANDS = ['Frequent > 1e-3', 'Probable', 'Remote', 'Extremely remote', 'Extremely improbable < 1e-9'], half = Math.log10(ORDER);
    const mat = N.zeros(4, 5), full = rows.map((r) => {
      const pFh = r[4] * r[5] * r[6], Cm = pFh * i.t_flight, occ = N.clamp(Math.round(10 + Math.log10(Math.max(pFh, 1e-10))), 1, 10), rpn = (4 - r[3]) * 2.5 * occ * r[7], margin = Math.log10(obj[r[3]] / Math.max(pFh, 1e-300));
      mat[3 - r[3]][4 - band(pFh)]++;
      // a product of rates that lands exactly on the objective must not fail by round-off; within half a decade is "of the same order"
      return { r, pFh, Cm, rpn, margin, verdict: r[8] ? `n/a — ${r[8]}` : margin >= -1e-9 ? 'yes' : margin >= -half ? 'same order' : 'NO' };
    }).sort((a, b) => (a.r[8] ? 1 : 0) - (b.r[8] ? 1 : 0) || a.margin - b.margin);
    const gov = full.filter((f) => !f.r[8]), sumBy = SEV.map((_, s) => N.sum(gov.filter((f) => f.r[3] === s).map((f) => f.pFh))), nBad = gov.filter((f) => f.verdict === 'NO').length, nNear = gov.filter((f) => f.verdict === 'same order').length, worst = gov[0];
    if (full.length > gov.length) warnings.push(`${full.length - gov.length} mode(s) are listed without a compliance verdict: engine failure is a design case (continued flight with an engine out, or a forced landing on a single), its consequence is quantified in the fault tree; structure and rotating dynamic components are shown by fatigue and damage-tolerance evaluation, and rotor burst or battery fire by design and test, not by a probability objective.`);
    if (i.unmanned) warnings.push('Unmanned aircraft: every effect is taken one severity class lower than for an occupied aircraft, and the objectives are placeholders (no sourced scheme for unmanned aircraft in the data registry).');
    warnings.push(`${ILLUSTRATIVE}. Objectives: ${schemeOf(i.obj_cat)}. Severity classes are generic and must come from your own functional hazard assessment.`);
    return {
      kpis: [
        { key: 'n_failure_modes', label: 'Failure modes assessed', value: rows.length, unit: '' },
        { key: 'n_noncompliant', label: 'Modes above their probability objective', value: nBad, unit: '', status: nBad ? (i.unmanned ? 'warn' : 'bad') : nNear ? 'warn' : 'ok', note: nNear ? `${nNear} more within half a decade of it` : 'More than half a decade above' },
        { key: 'worst_margin_decades', label: 'Smallest margin to objective', value: worst.margin, unit: 'decades', status: worst.margin >= -1e-9 ? 'ok' : worst.margin >= -half || i.unmanned ? 'warn' : 'bad', note: `${worst.r[0]}: ${worst.r[1]}` },
        { key: 'p_cat_single_modes_per_fh', label: 'Sum of catastrophic single-mode probabilities', value: sumBy[0], unit: '1/FH', status: sumBy[0] <= i.obj_cat ? 'ok' : 'warn', note: 'Modes judged against a probability objective' },
        { key: 'p_haz_modes_per_fh', label: 'Sum of hazardous mode probabilities', value: sumBy[1], unit: '1/FH' },
        { key: 'p_major_modes_per_fh', label: 'Sum of major mode probabilities', value: sumBy[2], unit: '1/FH' },
        { key: 'max_rpn', label: 'Highest risk priority number', value: Math.max(...full.map((f) => f.rpn)), unit: '-', note: 'Severity × occurrence × detection, each 1–10' },
      ],
      plots: [
        { type: 'heat', title: 'Risk matrix: number of failure modes per cell', xlabel: 'Probability band (1 extremely improbable … 5 frequent)', ylabel: 'Severity (1 minor … 4 catastrophic)', zlabel: 'Modes [-]', x: [1, 2, 3, 4, 5], y: [1, 2, 3, 4], z: mat, contours: 0, overlay: [{ name: 'Acceptability boundary (large-aeroplane bands)', x: [1.5, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5], y: [4.5, 3.5, 3.5, 2.5, 2.5, 1.5, 1.5, 0.5] }] },
        { type: 'bar', title: 'Probability per flight hour by severity class against objectives', ylabel: 'log10 of probability per FH [-]', categories: SEV, series: [{ name: 'Sum of modes', y: sumBy.map((v) => Math.log10(Math.max(v, 1e-15))) }, { name: 'Objective per condition', y: obj.map((v) => Math.log10(v)) }] },
      ],
      tables: [{ title: 'FMECA worksheet (least margin first)', columns: ['Item', 'Failure mode', 'Effect', 'Severity', 'λ item [1/FH]', 'α', 'β', 'P effect [1/FH]', 'Criticality Cm per flight', 'Detection', 'RPN', 'Margin [decades]', 'Compliant'], rows: full.map((f) => [f.r[0], f.r[1], f.r[2], SEV[f.r[3]], f.r[4], f.r[5], f.r[6], f.pFh, f.Cm, f.r[7], f.rpn, f.margin, f.verdict]) },
        { title: 'Probability bands', columns: ['Band', 'Label'], rows: BANDS.map((b, j) => [5 - j, b]) }],
      warnings, models: ['Failure Modes, Effects and Criticality Analysis (mode criticality Cm = β·α·λ·t)', 'Risk priority number', 'Severity–probability risk matrix', 'Hazard classification against quantitative objectives'],
      assumptions: ['Each row is an independent single failure mode of one item; combinations are treated in the fault tree analysis', 'Occurrence rank is 10 + log10(probability per FH), limited to 1–10', 'A mode within half a decade above its objective is reported as of the same order, because the objectives are orders of magnitude', 'The risk-matrix bands and boundary are the large-aeroplane ones whatever objectives are entered', ILLUSTRATIVE],
    };
  },
  verify() {
    const b = { n_eng: 2, eng_need: 1, n_hyd: 0, n_gen: 0, n_fcc: 0, n_batt: 1, unmanned: false, lam_engine: 1e-5, lam_hyd: 0, lam_gen: 0, lam_fcc: 0, lam_batt: 0, lam_struct: 1e-3, t_flight: 3, has_rotor: false, has_gearbox: false, has_gear: false, obj_cat: 1e-9, obj_haz: 1e-7, obj_maj: 1e-5, obj_min: 1e-3 };
    const r = fmeca.run({ ...b, extra_rows: 'Test item; test mode; 2; 4e-6; 0.5; 0.25; 3' }), row = r.tables[0].rows.find((x) => x[0] === 'Test item'), cells = N.sum(r.plots[0].z.map((z) => N.sum(z)));
    // 5e-5 × 0.1 × 0.02 is 1.0000000000000001e-7 in floating point: exactly on the hazardous objective
    const e = fmeca.run({ ...b, extra_rows: 'Edge; on the objective; 2; 5e-5; 0.1; 0.02; 3\nNear; twice the objective; 3; 2e-5; 1; 1; 3\nFar; ten times the objective; 3; 1e-4; 1; 1; 3' }), v = (nm) => e.tables[0].rows.find((x) => x[0] === nm)[12];
    const u = fmeca.run({ ...b, unmanned: true }), st = r.tables[0].rows.find((x) => x[0] === 'Primary structure');
    return [
      N.check('Mode criticality Cm = β·α·λ·t', row[8], 0.25 * 0.5 * 4e-6 * 3, 1e-12, 'MIL-STD-1629A mode criticality number'),
      N.check('Margin = log10(objective / probability)', row[11], Math.log10(1e-7 / 5e-7), 1e-10, 'Definition'),
      N.check('Every mode appears once in the risk matrix', cells, N.kv(r).n_failure_modes, 0, 'Count conservation'),
      N.check('A mode exactly on its objective is compliant', v('Edge') === 'yes' ? 1 : 0, 1, 0, 'No failure by floating-point round-off'),
      N.check('Twice the objective is of the same order, ten times is not', (v('Near') === 'same order' ? 1 : 0) + (v('Far') === 'NO' ? 1 : 0), 2, 0, 'Half-decade band: √10 ≈ 3.16'),
      N.check('Only the mode a decade above its objective is counted as non-compliant', N.kv(e).n_noncompliant, 1, 0, 'Count'),
      N.check('Structure carries no probability verdict even at a high rate', /^n\/a/.test(st[12]) ? 1 : 0, 1, 0, 'Shown by fatigue and damage-tolerance evaluation'),
      N.check('Unmanned aircraft: engine loss on a twin is one class lower (minor)', u.tables[0].rows.find((x) => x[1] === 'In-flight shutdown or loss of power')[3] === 'Minor' ? 1 : 0, 1, 0, 'Severity shift by one class'),
    ];
  },
  recommend(res, i) {
    const out = [], rows = res.tables[0].rows, bad = rows.filter((r) => r[12] === 'NO'), near = rows.filter((r) => r[12] === 'same order'), list = (a) => a.slice(0, 4).map((r) => `${r[0]} — ${r[1]} (${r[3]}, ${Number(r[7]).toExponential(1)}/FH)`).join('; ') + '.';
    // without a sourced objective for unmanned aircraft an exceedance cannot be more than advice
    if (bad.length) out.push({ severity: i.unmanned ? 'advise' : bad.some((r) => r[3] === 'Catastrophic') ? 'critical' : 'warn', title: `${bad.length} failure mode(s) exceed their probability objective`, detail: list(bad) + (i.unmanned ? ' The objectives for an unmanned aircraft are placeholders: set those of your operational approval.' : ''), action: 'For each: reduce the failure rate, add redundancy or monitoring so the effect becomes less severe, or limit exposure. Then confirm the combination logic in the fault tree analysis.', basis: 'Inverse relation between severity and allowable probability' });
    if (near.length) out.push({ severity: 'advise', title: `${near.length} failure mode(s) are of the same order as their objective`, detail: list(near), action: 'These sit within half a decade above the objective, inside the accuracy of placeholder rates. Replace the rate with supplier or fleet data, or confirm the severity class in the functional hazard assessment, before spending design effort.', basis: 'Probability bands defined as orders of magnitude (AC 25.1309-1B; AC 23.1309-1E Figure 2, Note 1)' });
    const latent = rows.filter((r) => r[9] >= 8);
    if (latent.length) out.push({ severity: 'advise', title: 'Poorly detectable failure modes', detail: latent.map((r) => `${r[0]} — ${r[1]}`).join('; ') + '.', action: 'Add built-in test, a scheduled functional check or a condition monitor. Earlier detection also avoids secondary damage and unplanned part replacement, which cuts material use and cost.', basis: 'Detection ranking ≥ 8' });
    out.push({ severity: 'info', title: 'Replace the placeholder data', detail: 'The table is generated with generic illustrative rates and severities.', action: 'Enter rates from supplier reliability reports or fleet records and severities from your functional hazard assessment; use the Bayesian updating analysis to combine prior estimates with service data.', basis: 'Data traceability' });
    return out;
  },
};
/** Series system of k-of-n blocks; each block {name, n, k, mttf, shape}. */
const blocksOf = (i) => [
  { name: 'Propulsion', n: Math.round(i.n_eng), k: N.clamp(Math.round(i.eng_need), 1, Math.round(i.n_eng)), mttf: 1 / Math.max(i.lam_engine, 1e-30), shape: i.shape_mech },
  { name: 'Hydraulics', n: Math.round(i.n_hyd), k: 1, mttf: 1 / Math.max(i.lam_hyd, 1e-30), shape: i.shape_mech },
  { name: 'Electrical generation', n: Math.round(i.n_gen), k: 1, mttf: 1 / Math.max(i.lam_gen, 1e-30), shape: 1 },
  { name: 'Flight-control computing', n: Math.round(i.n_fcc), k: 1, mttf: 1 / Math.max(i.lam_fcc, 1e-30), shape: 1 },
].filter((b) => b.n >= 1).map((b) => ({ ...b, scale: weibScale(b.mttf, b.shape) }));
const sysR = (blocks, t, age) => blocks.reduce((R, b) => R * kofn(b.k, b.n, weibR(t, b.shape, b.scale, age)), 1);
function sysMttf(blocks, age) {
  let tmax = Math.min(...blocks.map((b) => b.mttf)); for (let j = 0; j < 60 && sysR(blocks, tmax, age) > 1e-12; j++) tmax *= 1.6;
  return N.simpson((t) => sysR(blocks, t, age), 0, tmax, 3000);
}
const rbd = {
  id: 'rbd', title: 'Reliability block diagram, redundancy and mission reliability', fidelity: 'numerical',
  summary: 'Treats the vehicle as a series chain of redundant k-out-of-n blocks with exponential or Weibull life, giving mission reliability, mean time to failure and hazard-rate curves, and cross-checks the analytic result with a seeded Monte Carlo simulation.',
  equations: ['Reliability function equations', 'Hazard function equations', 'Exponential reliability distribution', 'Weibull distribution', 'Failure rate equations'],
  inputs: [...ARCH, ...RATES,
    { key: 't_mission', label: 'Mission duration', unit: 'h', default: 2, min: 0.01, max: 1e5, group: 'Exposure' },
    { key: 'age_h', label: 'Accumulated operating hours', unit: 'h', default: 0, min: 0, max: 2e5, group: 'Exposure', help: 'Component age at the start of the mission; matters only when the Weibull shape differs from 1' },
    { key: 'shape_mech', label: 'Weibull shape of mechanical items', unit: '-', default: 1, min: 0.3, max: 6, group: 'Life model', help: '1 = constant hazard (exponential); > 1 wear-out; < 1 infant mortality. Applied to propulsion and hydraulics.' },
    { key: 'bath_infant_frac', label: 'Bathtub: early-failure hazard at 100 h ÷ random rate', unit: '-', default: 1, min: 0, max: 50, group: 'Life model', help: 'Used only for the single-component bathtub illustration' },
    { key: 'bath_wear_life_h', label: 'Bathtub: wear-out characteristic life', unit: 'h', default: 20000, min: 10, max: 1e6, group: 'Life model' },
    { key: 'nSamples', label: 'Monte Carlo samples', unit: '', default: 4000, min: 200, max: 200000, step: 1, discrete: true, group: 'Numerics' },
    { key: 'seed', label: 'Random seed', unit: '', default: 22, min: 1, max: 1e9, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up) => ({ ...archDefaults(c), t_mission: flightHours(c), bath_wear_life_h: up.fatigue?.life_fh > 0 ? N.clamp(up.fatigue.life_fh, 100, 1e6) : undefined, lam_hyd: up.hydmech?.bearing_L10_h > 0 ? N.clamp(0.105 / up.hydmech.bearing_L10_h + 5e-5, 1e-6, 1e-2) : undefined }),
  run(i) {
    const blocks = blocksOf(i), warnings = [], Rm = sysR(blocks, i.t_mission, i.age_h), mttf = sysMttf(blocks, i.age_h);
    // Monte Carlo: sample remaining component lives, block fails at its (n−k+1)-th failure, system at the first block failure
    const u = N.rng(i.seed), ns = Math.round(i.nSamples); let surv = 0, sumT = 0, sumT2 = 0;
    for (let s = 0; s < ns; s++) {
      let tSys = Infinity;
      for (const b of blocks) {
        const lives = N.range(b.n, () => b.scale * (-Math.log(1 - u() * (1 - 1e-12)) + (i.age_h / b.scale) ** b.shape) ** (1 / b.shape) - i.age_h).sort((x, y) => x - y);
        tSys = Math.min(tSys, lives[b.n - b.k]);
      }
      if (tSys > i.t_mission) surv++; sumT += tSys; sumT2 += tSys * tSys;
    }
    const Rmc = surv / ns, mMc = sumT / ns, seM = Math.sqrt(Math.max(0, sumT2 / ns - mMc * mMc) / ns);
    const z = 1.96, wc = (Rmc + (z * z) / (2 * ns)) / (1 + (z * z) / ns), wh = (z * Math.sqrt((Rmc * (1 - Rmc)) / ns + (z * z) / (4 * ns * ns))) / (1 + (z * z) / ns); // Wilson interval
    const ts = N.linspace(0, Math.min(3 * mttf, 20 * Math.max(i.t_mission, mttf / 20)), 80), dt = 1e-3 * i.t_mission;
    const hazSys = -(Math.log(sysR(blocks, i.t_mission + dt, i.age_h)) - Math.log(sysR(blocks, Math.max(0, i.t_mission - dt), i.age_h))) / (i.t_mission + dt - Math.max(0, i.t_mission - dt));
    // bathtub illustration for one mechanical item: early (shape 0.5) + random + wear-out (shape 3.5) competing risks
    const lam0 = i.lam_engine, tb = N.logspace(1, 3 * i.bath_wear_life_h, 90), hz = (t, sh, sc) => (sh / sc) * (t / sc) ** (sh - 1), scI = (0.5 / Math.max(lam0 * i.bath_infant_frac, 1e-30)) ** 2 / 100; // shape-0.5 scale giving h(100 h) = frac·λ
    const hInf = tb.map((t) => (i.bath_infant_frac > 0 ? hz(t, 0.5, scI) : 0)), hWear = tb.map((t) => hz(t, 3.5, i.bath_wear_life_h));
    if (Rm > 0.9999 && ns < 1e5) warnings.push('Mission unreliability is far smaller than 1/samples, so the Monte Carlo mission estimate cannot resolve it; compare the mean time to failure instead, or raise the sample count.');
    if (i.shape_mech !== 1 && i.age_h === 0) warnings.push('A Weibull shape other than 1 is set but component age is zero: enter accumulated hours to see the ageing effect on mission reliability.');
    warnings.push(`${ILLUSTRATIVE}.`);
    return {
      kpis: [
        { key: 'system_reliability', label: 'Mission reliability (analytic)', value: Rm, unit: '-', status: Rm > 0.999 ? 'ok' : Rm > 0.99 ? 'warn' : 'bad', note: 'Probability that no block loses its function during the mission' },
        { key: 'mission_unreliability', label: 'Mission failure probability', value: 1 - Rm, unit: '-' },
        { key: 'mtbf_h', label: 'System mean time to failure', value: mttf, unit: 'h', note: 'Integral of the system reliability function' },
        { key: 'hazard_rate_per_h', label: 'System hazard rate during the mission', value: hazSys, unit: '1/h' },
        { key: 'mc_reliability', label: 'Mission reliability (Monte Carlo)', value: Rmc, unit: '-', note: `95% Wilson interval ${(wc - wh).toFixed(5)} … ${Math.min(1, wc + wh).toFixed(5)}` },
        { key: 'mc_mttf_h', label: 'Mean time to failure (Monte Carlo)', value: mMc, unit: 'h', status: Math.abs(mMc - mttf) <= 3.5 * seM ? 'ok' : 'warn', note: `± ${(1.96 * seM).toPrecision(3)} h (95%)` },
        { key: 'mc_mttf_err_pct', label: 'Monte Carlo versus analytic difference', value: (100 * (mMc - mttf)) / mttf, unit: '%' },
      ],
      plots: [
        { type: 'line', title: 'Reliability versus operating time', xlabel: 'Time [h]', ylabel: 'Reliability [-]', series: [...blocks.map((b) => ({ name: `${b.name} (${b.k}-of-${b.n})`, x: ts, y: ts.map((t) => kofn(b.k, b.n, weibR(t, b.shape, b.scale, i.age_h))) })), { name: 'System (series)', x: ts, y: ts.map((t) => sysR(blocks, t, i.age_h)), style: 'dash' }].slice(-6) },
        { type: 'line', title: 'Bathtub hazard-rate curve of a single mechanical item', xlabel: 'Age [h]', ylabel: 'Hazard rate [1/h]', xlog: true, ylog: true, series: [{ name: 'Early failures (shape 0.5)', x: tb, y: hInf.map((v) => Math.max(v, 1e-300)), style: 'dash' }, { name: 'Random failures', x: tb, y: tb.map(() => Math.max(lam0, 1e-300)), style: 'dash' }, { name: 'Wear-out (shape 3.5)', x: tb, y: hWear.map((v) => Math.max(v, 1e-300)), style: 'dash' }, { name: 'Total hazard', x: tb, y: tb.map((_, j) => Math.max(hInf[j] + lam0 + hWear[j], 1e-300)) }] },
      ],
      tables: [{ title: 'Block summary', columns: ['Block', 'Units n', 'Required k', 'Unit MTTF [h]', 'Weibull shape', 'Unit mission reliability', 'Block mission reliability', 'Block failure probability'], rows: blocks.map((b) => { const R = weibR(i.t_mission, b.shape, b.scale, i.age_h), Rb = kofn(b.k, b.n, R); return [b.name, b.n, b.k, b.mttf, b.shape, R, Rb, 1 - Rb]; }) }],
      warnings, models: ['Series reliability block diagram of k-out-of-n blocks', 'Exponential and two-parameter Weibull life models with conditional (aged) reliability', 'Seeded Monte Carlo life simulation with Wilson confidence interval', 'Competing-risk bathtub hazard illustration'],
      assumptions: ['Units within a block are identical, active and independent (no common cause here; see the fault tree)', 'No repair during the mission', 'Unit failure rates are converted to Weibull scale so that the unit MTTF is preserved', ILLUSTRATIVE],
    };
  },
  convergence: { param: 'nSamples', label: 'Monte Carlo samples', levels: [500, 1000, 2000, 4000, 8000], metric: 'mc_mttf_h', hOf: (n) => 1 / Math.sqrt(n) },
  verify() {
    const lam = 1e-3, b = [{ name: 'x', n: 4, k: 2, mttf: 1 / lam, shape: 1, scale: 1 / lam }], w = [{ name: 'w', n: 1, k: 1, mttf: weibScale(1, 1) * 500 * N.gamma(1 + 1 / 2.5), shape: 2.5, scale: 500 }];
    const r = N.kv(rbd.run({ n_eng: 3, eng_need: 2, n_hyd: 0, n_gen: 0, n_fcc: 0, lam_engine: 1e-3, lam_hyd: 0, lam_gen: 0, lam_fcc: 0, t_mission: 100, age_h: 0, shape_mech: 1, bath_infant_frac: 0, bath_wear_life_h: 1e4, nSamples: 40000, seed: 5 }));
    return [
      N.check('k-of-n exponential MTTF = (1/λ)·Σ 1/j, j = k…n', sysMttf(b, 0), (1 / lam) * (1 / 2 + 1 / 3 + 1 / 4), 1e-6, 'Order statistics of exponential lives'),
      N.check('Weibull MTTF = η·Γ(1 + 1/β)', sysMttf(w, 0), 500 * N.gamma(1.4), 1e-6, 'Weibull mean'),
      N.check('2-of-3 mission reliability 3R² − 2R³', r.system_reliability, 3 * Math.exp(-0.2) - 2 * Math.exp(-0.3), 1e-12, 'Binomial'),
      N.check('Monte Carlo MTTF agrees with the analytic value', r.mc_mttf_h, (1 / 1e-3) * (1 / 2 + 1 / 3), 0.015, 'Sampling error ≈ 0.4% at 40 000 samples'),
      N.check('Monte Carlo mission reliability agrees', r.mc_reliability, r.system_reliability, 0.006, 'Binomial sampling error'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], rows = res.tables[0].rows.slice().sort((a, b) => b[7] - a[7]), w = rows[0];
    if (w && w[7] > 0) out.push({ severity: o.system_reliability < 0.99 ? 'warn' : 'info', title: `${w[0]} is the weakest block`, detail: `It contributes ${(100 * w[7] / Math.max(N.sum(rows.map((r) => r[7])), 1e-300)).toFixed(0)}% of mission unreliability (${w[2]}-of-${w[1]} redundancy, unit MTTF ${Number(w[3]).toPrecision(3)} h).`, action: w[1] === w[2] ? 'There is no spare unit in this block: adding one redundant unit usually cuts its failure probability by orders of magnitude.' : 'Improve the unit failure rate or add a further redundant unit, and check for common-cause exposure in the fault tree.', basis: 'Series system: unreliabilities add' });
    if (i.shape_mech > 1.5) out.push({ severity: 'advise', title: 'Wear-out behaviour: schedule replacement', detail: `Weibull shape ${i.shape_mech} means the hazard rate rises with age.`, action: 'Use the preventive-replacement analysis to set a hard-time or on-condition limit; replacing just before wear-out keeps reliability high without discarding useful life.', basis: 'Increasing hazard rate (Weibull shape > 1)' });
    return out;
  },
};

/** Generator matrix of a k-of-n repairable system; state j = number of failed units. */
function markovQ(n, lam, mu, crews, melAllowed, muMel) {
  const Q = N.zeros(n + 1);
  for (let j = 0; j <= n; j++) {
    if (j < n) { Q[j][j + 1] = (n - j) * lam; Q[j][j] -= Q[j][j + 1]; }
    if (j > 0) { Q[j][j - 1] = j <= melAllowed ? muMel : Math.min(j, crews) * mu; Q[j][j] -= Q[j][j - 1]; }
  }
  return Q;
}
function steadyState(Q) { const n = Q.length, A = N.transpose(Q); A[n - 1] = new Array(n).fill(1); const b = new Array(n).fill(0); b[n - 1] = 1; return N.solve(A, b); }
const markov = {
  id: 'markov', title: 'Markov availability and dispatch reliability of a repairable redundant system', fidelity: 'numerical',
  summary: 'Models a redundant system whose units fail and are repaired as a continuous-time Markov chain: availability over time and in the long run, time between system failures, and the effect of deferring repairs under a minimum-equipment-list allowance on dispatch.',
  equations: ['Markov transition equations', 'Availability equations', 'Repairable system reliability equations', 'Failure rate equations'],
  inputs: [
    { key: 'n_units', label: 'Installed units', unit: '', default: 2, min: 1, max: 8, step: 1, discrete: true, group: 'Architecture' },
    { key: 'k_need', label: 'Units required for the function', unit: '', default: 1, min: 1, max: 8, step: 1, discrete: true, group: 'Architecture' },
    { key: 'lam', label: 'Unit failure rate', unit: '1/h', default: 2e-4, min: 1e-9, max: 10, group: 'Failure and repair (illustrative)' },
    { key: 'mttr_h', label: 'Mean time to repair', unit: 'h', default: 4, min: 0.01, max: 1e4, group: 'Failure and repair (illustrative)', help: 'Active repair including access, fault-finding and test' },
    { key: 'crews', label: 'Simultaneous repair crews', unit: '', default: 1, min: 1, max: 8, step: 1, discrete: true, group: 'Failure and repair (illustrative)' },
    { key: 'mel_allowed', label: 'Failed units allowed at dispatch (MEL)', unit: '', default: 1, min: 0, max: 7, step: 1, discrete: true, group: 'Dispatch', help: '0 = no dispatch with any unit failed' },
    { key: 'mel_interval_h', label: 'Mean time to rectify a deferred defect', unit: 'h', default: 72, min: 0.1, max: 3000, group: 'Dispatch', help: 'Deferred items are repaired at the next convenient opportunity inside the rectification interval' },
    { key: 'nogo_per_100dep', label: 'Other technical no-go events per 100 departures', unit: '-', default: 0.8, min: 0, max: 50, group: 'Dispatch', help: 'Rest-of-aircraft technical delays and cancellations; illustrative, replace with operator data' },
    { key: 't_end', label: 'Simulated time', unit: 'h', default: 200, min: 1, max: 1e5, group: 'Numerics' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 400, min: 20, max: 20000, step: 1, discrete: true, group: 'Numerics' },
  ],
  // the modelled units are the generators when the aircraft has them, otherwise the electric motors (with their own failure rate)
  defaults: (c) => { const motors = !(c.systems.gen_kVA > 0) && c.prop.type === 'electric', n = c.systems.gen_kVA > 0 ? Math.max(1, Math.min(c.prop.n_eng, 2)) : motors ? Math.max(2, Math.min(c.prop.n_eng, 8)) : 1, need = motors ? (n >= 6 ? n - 1 : n) : 1; return { n_units: n, k_need: need, mel_allowed: Math.max(0, Math.min(1, n - need)), lam: motors ? archDefaults(c).lam_engine : undefined }; },
    run(i) {
    const n = Math.round(i.n_units), k = N.clamp(Math.round(i.k_need), 1, n), mu = 1 / i.mttr_h, mel = N.clamp(Math.round(i.mel_allowed), 0, n - k), warnings = [];
    const Q = markovQ(n, i.lam, mu, Math.round(i.crews), mel, 1 / i.mel_interval_h), up = (j) => j <= n - k, pi = steadyState(Q);
    // stiff-safe step count: explicit RK4 needs h·|q|max < 2.78
    const qmax = Math.max(...Q.map((r, j) => -r[j])), nSt = Math.max(Math.round(i.nSteps), Math.ceil((i.t_end * qmax) / 2)), y0 = N.range(n + 1, (j) => (j === 0 ? 1 : 0));
    const sol = N.rk4((t, y) => y.map((_, c) => { let s = 0; for (let r = 0; r <= n; r++) s += y[r] * Q[r][c]; return s; }), 0, y0, i.t_end, nSt), st = Math.max(1, Math.floor(nSt / 300));
    const tt = [], At = [], P0 = [], Pdown = [];
    for (let s = 0; s <= nSt; s += st) { const y = sol.y[s]; tt.push(sol.t[s]); At.push(N.sum(y.filter((_, j) => up(j)))); P0.push(y[0]); Pdown.push(1 - At[At.length - 1]); }
    const yEnd = sol.y[nSt], Aend = N.sum(yEnd.filter((_, j) => up(j))), Ass = N.sum(pi.filter((_, j) => up(j)));
    let freq = 0; for (let j = 0; j <= n; j++) for (let m = 0; m <= n; m++) if (up(j) && !up(m)) freq += pi[j] * Q[j][m];
    // mean time to first system failure: absorbing down states, Q_TT·τ = −1
    const nt = n - k + 1, tau = N.solve(N.range(nt, (a) => N.range(nt, (b) => Q[a][b])), new Array(nt).fill(-1));
    const pDisp = N.sum(pi.filter((_, j) => j <= mel)), dispatch = pDisp * (1 - i.nogo_per_100dep / 100);
    const piNoMel = steadyState(markovQ(n, i.lam, mu, Math.round(i.crews), 0, mu)), availNoMel = N.sum(piNoMel.filter((_, j) => up(j)));
    if (nSt > Math.round(i.nSteps)) warnings.push(`Time steps were raised to ${nSt} to keep the explicit integration stable for this repair rate.`);
    if (mel > 0 && Ass < availNoMel - 1e-12) warnings.push('Deferring repairs under the MEL improves dispatch but lowers the availability of the function itself, because the system flies longer with reduced redundancy.');
    warnings.push(`${ILLUSTRATIVE}; the rest-of-aircraft no-go rate is an operator statistic you must supply.`);
    return {
      kpis: [
        { key: 'availability', label: 'Steady-state availability of the function', value: Ass, unit: '-', status: Ass > 0.999 ? 'ok' : Ass > 0.99 ? 'warn' : 'bad' },
        { key: 'availability_t', label: `Availability at ${i.t_end} h`, value: Aend, unit: '-' },
        { key: 'unavailability', label: 'Steady-state unavailability', value: 1 - Ass, unit: '-' },
        { key: 'mtbf_system_h', label: 'Mean time between system failures', value: freq > 0 ? 1 / freq : Infinity, unit: 'h' },
        { key: 'mttff_h', label: 'Mean time to first system failure', value: tau[0], unit: 'h' },
        { key: 'mdt_h', label: 'Mean system down time per failure', value: freq > 0 ? (1 - Ass) / freq : 0, unit: 'h' },
        { key: 'p_full_redundancy', label: 'Time share with every unit serviceable', value: pi[0], unit: '-' },
        { key: 'dispatch_reliability', label: 'Technical dispatch reliability', value: dispatch, unit: '-', status: dispatch > 0.99 ? 'ok' : dispatch > 0.97 ? 'warn' : 'bad', note: 'P(system dispatchable) × (1 − other no-go events per departure)' },
        { key: 'dispatch_gain_mel_pct', label: 'Dispatch gain from MEL deferral', value: 100 * (pDisp - piNoMel[0]), unit: '%-points' },
      ],
      plots: [
        { type: 'line', title: 'Availability build-down from an all-serviceable start', xlabel: 'Time [h]', ylabel: 'Probability [-]', series: [{ name: 'Function available', x: tt, y: At }, { name: 'All units serviceable', x: tt, y: P0 }], annotations: [{ y: Ass, label: 'Steady state' }] },
        { type: 'bar', title: 'Long-run state probabilities', ylabel: 'log10 of probability [-]', categories: pi.map((_, j) => `${j} failed`), series: [{ name: 'With MEL deferral', y: pi.map((v) => Math.log10(Math.max(v, 1e-30))) }, { name: 'Repair at once', y: piNoMel.map((v) => Math.log10(Math.max(v, 1e-30))) }] },
      ],
      tables: [{ title: 'Markov states', columns: ['Failed units', 'Function', 'Dispatch', 'Steady-state probability', 'Failure rate out [1/h]', 'Repair rate out [1/h]'], rows: pi.map((v, j) => [j, up(j) ? 'available' : 'lost', j <= mel ? 'go' : 'no-go', v, j < n ? Q[j][j + 1] : 0, j > 0 ? Q[j][j - 1] : 0]) }],
      warnings, models: ['Continuous-time Markov chain (birth–death) for a k-out-of-n repairable system', 'Kolmogorov forward equations integrated by RK4', 'Steady state from the balance equations', 'MEL deferred-repair states'],
      assumptions: ['Exponential failure and repair times', 'Surviving units keep operating while others are under repair', 'Deferred defects are rectified at a constant rate 1/interval', 'Dispatch reliability multiplies this system by an independent rest-of-aircraft no-go rate', ILLUSTRATIVE],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [100, 200, 400, 800], metric: 'availability_t' },
  verify() {
    const lam = 0.01, mttr = 20, mu = 1 / mttr, base = { crews: 1, mel_allowed: 0, mel_interval_h: 1, nogo_per_100dep: 0, nSteps: 400 };
    const a = N.kv(markov.run({ ...base, n_units: 1, k_need: 1, lam, mttr_h: mttr, t_end: 30 })), b = N.kv(markov.run({ ...base, n_units: 2, k_need: 1, lam, mttr_h: mttr, crews: 2, t_end: 2000 }));
    return [
      N.check('Two-state steady availability μ/(λ+μ)', a.availability, mu / (lam + mu), 1e-12, 'Closed form'),
      N.check('Two-state transient A(t) = μ/(λ+μ) + λ/(λ+μ)·exp(−(λ+μ)t)', a.availability_t, mu / (lam + mu) + (lam / (lam + mu)) * Math.exp(-(lam + mu) * 30), 1e-8, 'Closed form'),
      N.check('Parallel pair with two crews: 1 − (λ/(λ+μ))²', b.availability, 1 - (lam / (lam + mu)) ** 2, 1e-12, 'Independent repairable units'),
      N.check('Single unit MTBF = 1/λ + 1/μ', a.mtbf_system_h, 1 / lam + mttr, 1e-10, 'Renewal cycle'),
      N.check('Parallel pair mean time to first failure (3λ+μ)/(2λ²)', b.mttff_h, (3 * lam + mu) / (2 * lam * lam), 1e-10, 'Standard result for a 1-of-2 repairable system'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.dispatch_reliability < 0.99) out.push({ severity: 'warn', title: 'Dispatch reliability below 99%', detail: `${(100 * o.dispatch_reliability).toFixed(2)}% of departures leave without a technical hold.`, action: o.dispatch_gain_mel_pct < 0.01 && i.n_units > i.k_need ? 'Allow dispatch with one unit inoperative (subject to an approved MEL) and stock the spare at the main base.' : 'Attack the largest no-go drivers: lower the failure rate, cut repair time with better access and fault isolation, or position spares.', basis: 'Steady-state share of dispatchable states' });
    if (o.availability < 0.999) out.push({ severity: 'advise', title: 'Availability is repair-limited', detail: `Unavailability ${o.unavailability.toExponential(2)} with MTTR ${i.mttr_h} h.`, action: 'Unavailability scales roughly with (λ·MTTR)^(spares + 1): halving MTTR on a duplex system cuts it about four-fold. Feed the availability to Suite 26 to value it.', basis: 'Markov steady state' });
    out.push({ severity: 'info', title: 'Use availability in the economics', detail: `Availability ${o.availability.toFixed(5)} and dispatch reliability ${o.dispatch_reliability.toFixed(4)} are published for Suite 26.`, action: 'Cancelled sorties and standby aircraft are priced there; higher dispatch reliability also avoids ferry and positioning flights and their fuel and CO₂.', basis: 'Reliability–availability–revenue coupling' });
    return out;
  },
};

const eventTree = {
  id: 'eventtree', title: 'Event tree for an engine-fire initiating event', fidelity: 'analytical',
  summary: 'Follows an initiating event through detection, shutdown, two extinguisher shots and the diversion, multiplying branch probabilities to get the frequency and severity of every outcome.',
  equations: ['Event tree probability equations', 'Conditional probability equations', 'Total probability theorem'],
  inputs: [
    { key: 'f_init', label: 'Initiating event frequency (engine fire)', unit: '1/FH', default: 1e-6, min: 0, max: 1, group: 'Initiating event (illustrative)', help: 'Per aircraft flight hour. Placeholder.' },
    { key: 'p_detect', label: 'P(fire detected)', unit: '-', default: 0.999, min: 0, max: 1, group: 'Branch probabilities (illustrative)' },
    { key: 'p_shutdown', label: 'P(engine shut down and isolated)', unit: '-', default: 0.995, min: 0, max: 1, group: 'Branch probabilities (illustrative)' },
    { key: 'p_ext1', label: 'P(first extinguisher shot succeeds)', unit: '-', default: 0.9, min: 0, max: 1, group: 'Branch probabilities (illustrative)' },
    { key: 'p_ext2', label: 'P(second shot succeeds)', unit: '-', default: 0.8, min: 0, max: 1, group: 'Branch probabilities (illustrative)', help: 'Set 0 when only one bottle is fitted' },
    { key: 'p_land', label: 'P(safe landing with a persisting fire)', unit: '-', default: 0.95, min: 0, max: 1, group: 'Branch probabilities (illustrative)', help: 'Depends on diversion time and fire containment; 0.99 is filled for flights under half an hour' },
    { key: 'p_land_undetected', label: 'P(safe landing after an undetected fire)', unit: '-', default: 0.5, min: 0, max: 1, group: 'Branch probabilities (illustrative)' },
    ...OBJECTIVES,
  ],
  defaults: (c) => ({ ...objDefaults(c), f_init: c.prop.type === 'electric' ? 2e-7 : c.prop.type === 'piston' ? 3e-6 : 1e-6, p_ext2: c.mass.mtow_kg > 5700 ? 0.8 : 0, p_ext1: c.mass.mtow_kg > 5700 || c.meta.type === 'helicopter' ? 0.9 : 0, p_land: flightHours(c) < 0.5 ? 0.99 : 0.95 }),
  run(i) {
    const d = i.p_detect, s = i.p_shutdown, e1 = i.p_ext1, e2 = i.p_ext2, l = i.p_land, lu = i.p_land_undetected;
    const seq = [ // [description, conditional probability, severity index]
      ['Detected, shut down, out on first shot', d * s * e1, 2], ['Detected, shut down, out on second shot', d * s * (1 - e1) * e2, 2],
      ['Detected, shut down, fire persists, safe landing', d * s * (1 - e1) * (1 - e2) * l, 1], ['Detected, shut down, fire persists, loss of aircraft', d * s * (1 - e1) * (1 - e2) * (1 - l), 0],
      ['Detected, shutdown fails, safe landing', d * (1 - s) * l, 1], ['Detected, shutdown fails, loss of aircraft', d * (1 - s) * (1 - l), 0],
      ['Not detected, safe landing', (1 - d) * lu, 1], ['Not detected, loss of aircraft', (1 - d) * (1 - lu), 0],
    ];
    const f = seq.map((q) => q[1] * i.f_init), by = SEV.map((_, k) => N.sum(f.filter((_, j) => seq[j][2] === k))), obj = objOf(i), worst = N.argmax(seq.map((q) => (q[2] === 0 ? q[1] : -1)));
    const st = (k) => (by[k] <= obj[k] ? 'ok' : by[k] <= ORDER * obj[k] ? 'warn' : k === 0 ? 'bad' : 'warn');
    return {
      kpis: [
        { key: 'fire_catastrophic_per_fh', label: 'Catastrophic outcome frequency', value: by[0], unit: '1/FH', status: st(0), note: `Objective ${i.obj_cat.toExponential(0)}; within half a decade counts as the same order` },
        { key: 'fire_hazardous_per_fh', label: 'Hazardous outcome frequency', value: by[1], unit: '1/FH', status: st(1) },
        { key: 'fire_major_per_fh', label: 'Major outcome frequency', value: by[2], unit: '1/FH', status: st(2) },
        { key: 'p_cat_given_fire', label: 'P(catastrophic | fire)', value: i.f_init > 0 ? by[0] / i.f_init : N.sum(seq.filter((q) => q[2] === 0).map((q) => q[1])), unit: '-' },
        { key: 'p_sum_check', label: 'Sum of sequence probabilities', value: N.sum(seq.map((q) => q[1])), unit: '-', note: 'Must equal 1' },
      ],
      plots: [{ type: 'bar', title: 'Outcome frequency of each event-tree sequence', ylabel: 'log10 of frequency per FH [-]', categories: seq.map((_, j) => `S${j + 1}`), series: [{ name: 'Sequence frequency', y: f.map((v) => Math.log10(Math.max(v, 1e-30))) }] },
        { type: 'bar', title: 'Outcome frequency by severity against objectives', ylabel: 'log10 of frequency per FH [-]', categories: SEV.slice(0, 3), series: [{ name: 'Event tree', y: by.slice(0, 3).map((v) => Math.log10(Math.max(v, 1e-30))) }, { name: 'Objective', y: obj.slice(0, 3).map((v) => Math.log10(v)) }] }],
      tables: [{ title: 'Event-tree sequences', columns: ['Sequence', 'Path', 'Conditional probability', 'Frequency [1/FH]', 'Severity'], rows: seq.map((q, j) => [`S${j + 1}`, q[0], q[1], f[j], SEV[q[2]]]) }],
      outputs: { dominant_cat_sequence: worst + 1 },
      warnings: [`${ILLUSTRATIVE}. Initiating frequency and branch probabilities depend on the installation (detection loops, bottle count, diversion time) and must be substantiated. Objectives: ${schemeOf(i.obj_cat)}.`],
      models: ['Event tree with five pivotal events', 'Outcome classification against quantitative objectives'],
      assumptions: ['Branch probabilities are conditional on the path taken and otherwise independent', 'One initiating event at a time', 'Outcome severities are fixed per sequence'],
    };
  },
    verify() {
    const r = eventTree.run({ f_init: 1e-6, p_detect: 0.9, p_shutdown: 0.8, p_ext1: 0.5, p_ext2: 0.5, p_land: 0.9, p_land_undetected: 0.5, obj_cat: 1e-9, obj_haz: 1e-7, obj_maj: 1e-5, obj_min: 1e-3 }), o = N.kv(r);
    return [
      N.check('Sequence probabilities sum to one', o.p_sum_check, 1, 1e-14, 'Total probability theorem'),
      N.check('Catastrophic frequency by hand', o.fire_catastrophic_per_fh, 1e-6 * (0.9 * 0.8 * 0.25 * 0.1 + 0.9 * 0.2 * 0.1 + 0.1 * 0.5), 1e-12, 'Product of branch probabilities along each path'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], s = res.tables[0].rows[o.dominant_cat_sequence - 1];
    const near = o.fire_catastrophic_per_fh <= ORDER * i.obj_cat;
    if (o.fire_catastrophic_per_fh > i.obj_cat) out.push({ severity: near ? 'advise' : 'warn', title: near ? 'Fire event tree is of the same order as the catastrophic objective' : 'Fire event tree exceeds the catastrophic objective', detail: `${o.fire_catastrophic_per_fh.toExponential(2)} per FH against ${i.obj_cat.toExponential(0)}${near ? ' (within half a decade, inside the accuracy of placeholder branch probabilities)' : ''}; the dominant path is "${s[1]}".`, action: /Not detected/.test(s[1]) ? 'Improve detection coverage (dual loops, fault monitoring): undetected fire dominates.' : /shutdown fails/.test(s[1]) ? 'Improve isolation reliability (independent fuel and hydraulic shut-off valves).' : 'Add extinguishing capacity or shorten the maximum diversion time.', basis: 'Dominant event-tree sequence' });
    else out.push({ severity: 'info', title: 'Fire sequences meet the objective with these inputs', detail: `Catastrophic frequency ${o.fire_catastrophic_per_fh.toExponential(2)} per FH.`, action: 'Substantiate each branch probability with test or service evidence before relying on it.', basis: 'Event tree against the quantitative objective' });
    return out;
  },
};

const parseList = (text) => { const tok = String(text || '').split(/[\s,;]+/).filter(Boolean), v = tok.map(Number); return { values: v.filter((x) => Number.isFinite(x) && x > 0), bad: v.filter((x) => !(Number.isFinite(x) && x > 0)).length }; };
/** Weibull fit by median-rank regression (Bernard's approximation) on complete failure data. */
function weibullMRR(times) {
  const t = times.slice().sort((a, b) => a - b), n = t.length, X = t.map(Math.log), Y = t.map((_, j) => Math.log(-Math.log(1 - (j + 1 - 0.3) / (n + 0.4))));
  const mx = N.mean(X), my = N.mean(Y); let sxy = 0, sxx = 0; for (let j = 0; j < n; j++) { sxy += (X[j] - mx) * (Y[j] - my); sxx += (X[j] - mx) ** 2; }
  const shape = sxy / sxx, scale = Math.exp(mx - my / shape);
  return { shape, scale, r2: N.corr(X, Y) ** 2, X, Y, t };
}
const bayes = {
  id: 'bayes', title: 'Failure-rate updating from service data and Weibull life fitting', fidelity: 'analytical',
  summary: 'Combines a prior estimate of a failure rate with observed fleet hours and failures (Gamma–Poisson Bayes update) to give a posterior rate with credible interval, and fits a Weibull life distribution to a list of failure times.',
  equations: ['Bayes’ theorem', 'Failure rate equations', 'Weibull distribution', 'Exponential reliability distribution', 'Reliability function equations'],
  inputs: [
    { key: 'prior_lambda', label: 'Prior failure rate', unit: '1/h', default: 1e-4, min: 1e-12, max: 10, group: 'Prior', help: 'Handbook, supplier prediction or similar-equipment value' },
    { key: 'prior_strength', label: 'Prior strength (equivalent failures)', unit: '-', default: 1, min: 0.01, max: 1000, group: 'Prior', help: 'Gamma shape: 0.5–1 for a weak prior, larger when the prior rests on substantial evidence' },
    { key: 'fleet_hours', label: 'Observed fleet operating hours', unit: 'h', default: 50000, min: 0, max: 1e10, group: 'Service data' },
    { key: 'failures', label: 'Observed failures', unit: '', default: 3, min: 0, max: 1e6, step: 1, discrete: true, group: 'Service data' },
    { key: 'cred', label: 'Credible / confidence level', unit: '-', default: 0.9, min: 0.5, max: 0.999, group: 'Service data' },
    { key: 'times_text', label: 'Failure times for the Weibull fit', type: 'text', default: '410, 620, 780, 905, 1010, 1150, 1290, 1420, 1610, 1890', group: 'Weibull data', help: 'Times to failure in hours, separated by commas or spaces. The default list is synthetic and for illustration only. Complete (uncensored) data are assumed.' },
  ],
  defaults: () => ({}),
  run(i) {
    const warnings = [], a0 = i.prior_strength, b0 = a0 / i.prior_lambda, n = Math.round(i.failures), a1 = a0 + n, b1 = b0 + i.fleet_hours, lo = (1 - i.cred) / 2;
    const post = a1 / b1, qLo = gammaQuantile(lo, a1, b1), qHi = gammaQuantile(1 - lo, a1, b1), mle = i.fleet_hours > 0 ? n / i.fleet_hours : NaN, upper = i.fleet_hours > 0 ? gammaQuantile(i.cred, n + 1, i.fleet_hours) : NaN;
    const xs = N.linspace(0, Math.max(gammaQuantile(0.999, a0, b0), gammaQuantile(0.999, a1, b1)) * 1.02, 160).slice(1);
    const { values: times, bad } = parseList(i.times_text); let w = null;
    if (bad) warnings.push(`${bad} entr${bad > 1 ? 'ies' : 'y'} in the failure-time list could not be read as positive numbers and were ignored.`);
    if (times.length >= 3 && N.amax(times) > N.amin(times)) w = weibullMRR(times); else warnings.push('At least three distinct failure times are needed for the Weibull fit; it was skipped.');
    if (w && times.length < 8) warnings.push('Fewer than eight failure times: the fitted shape has wide uncertainty.');
    if (!(i.fleet_hours > 0)) warnings.push('No service hours entered: the posterior equals the prior.');
    const prPredZero = (b1 / (b1 + 1000)) ** a1, kp = [
      { key: 'lambda_post_per_h', label: 'Posterior mean failure rate', value: post, unit: '1/h' },
      { key: 'lambda_lo_per_h', label: `Lower ${(100 * i.cred).toFixed(0)}% credible bound`, value: qLo, unit: '1/h' },
      { key: 'lambda_hi_per_h', label: `Upper ${(100 * i.cred).toFixed(0)}% credible bound`, value: qHi, unit: '1/h' },
      { key: 'mtbf_post_h', label: 'Posterior MTBF', value: 1 / post, unit: 'h' },
      { key: 'lambda_mle_per_h', label: 'Service-data-only estimate (failures / hours)', value: mle, unit: '1/h' },
      { key: 'lambda_upper_classical_per_h', label: `One-sided ${(100 * i.cred).toFixed(0)}% upper confidence bound`, value: upper, unit: '1/h', note: 'χ² bound, valid also for zero failures' },
      { key: 'prior_weight_pct', label: 'Weight of the prior in the posterior', value: (100 * b0) / b1, unit: '%' },
      { key: 'p_no_failure_next_1000h', label: 'P(no failure in the next 1000 h)', value: prPredZero, unit: '-', note: 'Posterior predictive (negative binomial)' },
    ];
    const plots = [{ type: 'line', title: 'Prior and posterior density of the failure rate', xlabel: 'Failure rate [1/h]', ylabel: 'Probability density [h]', series: [{ name: 'Prior', x: xs, y: xs.map((x) => gammaPdf(x, a0, b0)), style: 'dash' }, { name: 'Posterior', x: xs, y: xs.map((x) => gammaPdf(x, a1, b1)) }], annotations: Number.isFinite(mle) ? [{ x: mle, label: 'Data only' }] : [] }];
    const tables = [{ title: 'Gamma–Poisson update', columns: ['Quantity', 'Prior', 'Posterior'], rows: [['Shape a (failures)', a0, a1], ['Rate b (hours)', b0, b1], ['Mean [1/h]', a0 / b0, post], ['Standard deviation [1/h]', Math.sqrt(a0) / b0, Math.sqrt(a1) / b1]] }];
    if (w) {
      const b10 = w.scale * (-Math.log(0.9)) ** (1 / w.shape), fx = [Math.min(...w.X), Math.max(...w.X)];
      kp.push({ key: 'weibull_shape', label: 'Weibull shape β', value: w.shape, unit: '-', note: w.shape > 1.2 ? 'Wear-out' : w.shape < 0.8 ? 'Infant mortality' : 'Near-random failures' }, { key: 'weibull_scale_h', label: 'Weibull characteristic life η', value: w.scale, unit: 'h' },
        { key: 'weibull_mttf_h', label: 'Weibull mean life', value: w.scale * N.gamma(1 + 1 / w.shape), unit: 'h' }, { key: 'weibull_b10_h', label: 'B10 life', value: b10, unit: 'h', note: 'Age by which 10% have failed' }, { key: 'weibull_r2', label: 'Probability-plot fit R²', value: w.r2, unit: '-', status: w.r2 > 0.9 ? 'ok' : 'warn' });
      plots.push({ type: 'line', title: 'Weibull probability plot (median ranks)', xlabel: 'ln(time to failure [h])', ylabel: 'ln(−ln(1 − F)) [-]', series: [{ name: 'Failure data', x: w.X, y: w.Y, style: 'points' }, { name: 'Fitted line', x: fx, y: fx.map((x) => w.shape * (x - Math.log(w.scale))) }] });
      tables.push({ title: 'Failure data and median ranks', columns: ['Rank', 'Time [h]', 'Median rank F'], rows: w.t.map((t, j) => [j + 1, t, (j + 1 - 0.3) / (w.t.length + 0.4)]) });
    }
    return { kpis: kp, plots, tables, warnings, models: ['Gamma–Poisson conjugate Bayesian update', 'Posterior credible interval from the incomplete gamma function', 'Two-parameter Weibull by median-rank regression (Bernard)'],
      assumptions: ['Failures follow a homogeneous Poisson process with a constant but uncertain rate', 'The prior is a Gamma distribution with mean equal to the prior rate', 'Weibull data are complete: no suspensions or censored units', 'The default failure-time list is synthetic'] };
  },
  verify() {
    const o = N.kv(bayes.run({ prior_lambda: 1e-3, prior_strength: 2, fleet_hours: 8000, failures: 6, cred: 0.9, times_text: '' }));
    const beta = 2.5, eta = 1000, n = 12, syn = N.range(n, (j) => eta * (-Math.log(1 - (j + 0.7) / (n + 0.4))) ** (1 / beta)), w = weibullMRR(syn);
    return [
      N.check('Posterior mean (a + n)/(b + T)', o.lambda_post_per_h, 8 / 10000, 1e-12, 'Conjugate Gamma–Poisson'),
      N.check('Gamma quantile reduces to exponential for a = 1', gammaQuantile(0.9, 1, 2), -Math.log(0.1) / 2, 1e-9, 'Exponential quantile'),
      N.check('Regularised incomplete gamma P(3, 2.5)', gammP(3, 2.5), 1 - Math.exp(-2.5) * (1 + 2.5 + 2.5 * 2.5 / 2), 1e-10, 'Closed form for integer shape'),
      N.check('Median-rank regression recovers the Weibull shape', w.shape, beta, 1e-9, 'Data placed exactly on the median ranks'),
      N.check('Median-rank regression recovers the Weibull scale', w.scale, eta, 1e-9, 'Data placed exactly on the median ranks'),
    ];
  },
  calibration: { params: [{ key: 'prior_lambda', min: 1e-9, max: 1 }], sweep: 'fleet_hours', target: 'lambda_post_per_h', note: 'Supply the observed cumulative failure rate at successive fleet-hour milestones to back out a consistent prior.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.prior_weight_pct > 50) out.push({ severity: 'advise', title: 'The result is dominated by the prior', detail: `The prior carries ${o.prior_weight_pct.toFixed(0)}% of the weight.`, action: 'Collect more service hours before relying on the posterior, or justify the prior strength with documented evidence.', basis: 'Ratio of prior to posterior Gamma rate parameter' });
    if (o.lambda_post_per_h > 1.5 * i.prior_lambda) out.push({ severity: 'warn', title: 'Service experience is worse than predicted', detail: `Posterior ${o.lambda_post_per_h.toExponential(2)}/h against a prior of ${i.prior_lambda.toExponential(2)}/h.`, action: 'Update the fault tree and FMECA with the posterior upper bound and review whether safety objectives are still met; investigate root cause.', basis: 'Bayesian update of the failure rate' });
    if (o.weibull_shape > 1.5) out.push({ severity: 'advise', title: 'Failure times show wear-out', detail: `Shape β = ${o.weibull_shape.toFixed(2)}, B10 life ${o.weibull_b10_h.toFixed(0)} h.`, action: 'Set a preventive replacement or overhaul interval with the next analysis (enter this shape and scale). It cuts unscheduled removals and secondary damage.', basis: 'Weibull shape > 1 implies an increasing hazard rate' });
    if (o.weibull_shape < 0.8) out.push({ severity: 'advise', title: 'Failure times show infant mortality', detail: `Shape β = ${o.weibull_shape.toFixed(2)}.`, action: 'Preventive replacement would make reliability worse. Look at manufacturing quality, burn-in and installation errors instead.', basis: 'Weibull shape < 1 implies a decreasing hazard rate' });
    return out;
  },
};

/** Age-replacement cost rate: [Cp·R(T) + Cf·F(T)] / ∫₀ᵀ R dt. */
function ageCost(T, shape, scale, Cp, Cf) { const R = weibR(T, shape, scale), up = N.simpson((t) => weibR(t, shape, scale), 0, T, 400); return { rate: (Cp * R + Cf * (1 - R)) / up, R, up }; }
const maint = {
  id: 'maint', title: 'Optimal preventive-replacement interval', fidelity: 'numerical',
  summary: 'Finds the replacement age that minimises cost per operating hour for a part with Weibull life, trading the cost of planned replacements against the higher cost of in-service failures, and compares it with running to failure.',
  equations: ['Weibull distribution', 'Reliability function equations', 'Hazard function equations', 'Repairable system reliability equations'],
  inputs: [
    { key: 'shape', label: 'Weibull shape β', unit: '-', default: 3, min: 0.5, max: 8, group: 'Life model', help: 'From the Weibull fit; must exceed 1 for preventive replacement to pay' },
    { key: 'scale_h', label: 'Weibull characteristic life η', unit: 'h', default: 20000, min: 1, max: 1e7, group: 'Life model', help: 'Linked to the fatigue life from Suite 9 when available' },
    { key: 'cost_planned', label: 'Cost of a planned replacement', unit: 'USD', default: 20000, min: 0, max: 1e9, group: 'Costs (illustrative)', help: 'Part, labour and planned downtime' },
    { key: 'cost_failure', label: 'Cost of an in-service failure', unit: 'USD', default: 150000, min: 0, max: 1e10, group: 'Costs (illustrative)', help: 'Unscheduled removal, secondary damage, delay or cancellation, recovery' },
    { key: 'r_min', label: 'Minimum reliability at replacement', unit: '-', default: 0, min: 0, max: 0.999999, group: 'Limits', help: 'Optional safety floor: the interval is shortened so that R(T) stays above it. 0 = cost only.' },
    { key: 'insp_interval_h', label: 'Current inspection / replacement interval', unit: 'h', default: 0, min: 0, max: 1e7, group: 'Limits', help: 'For comparison; 0 to ignore. Linked to Suite 9.' },
  ],
  defaults: (c, up) => ({ scale_h: up.fatigue?.life_fh > 0 ? N.clamp(up.fatigue.life_fh, 10, 1e7) : undefined, insp_interval_h: up.fatigue?.inspection_interval_fh > 0 ? up.fatigue.inspection_interval_fh : undefined, cost_planned: Math.max(500, 2e-4 * (c.econ.price_usd || 1e7)), cost_failure: Math.max(3000, 1.5e-3 * (c.econ.price_usd || 1e7)) }),
  run(i) {
    const mttf = i.scale_h * N.gamma(1 + 1 / i.shape), rtf = i.cost_failure / mttf, warnings = [], f = (T) => ageCost(T, i.shape, i.scale_h, i.cost_planned, i.cost_failure).rate;
    const Ts = N.logspace(0.02 * i.scale_h, 4 * i.scale_h, 120), cs = Ts.map(f), jb = N.argmin(cs);
    let Topt = Infinity, cOpt = rtf;
    if (i.shape > 1 && i.cost_failure > i.cost_planned && jb < Ts.length - 1 && cs[jb] < rtf * (1 - 1e-6)) { Topt = N.goldenSection(f, Ts[Math.max(0, jb - 1)], Ts[jb + 1], 1e-9); cOpt = f(Topt); }
    else warnings.push(i.shape <= 1 ? 'With a Weibull shape of 1 or less the hazard rate does not rise with age, so preventive replacement never lowers cost: run to failure (or on condition).' : 'A planned replacement costs as much as a failure, or the optimum lies beyond four characteristic lives: run to failure is the cost optimum.');
    const Tsafe = i.r_min > 0 ? i.scale_h * (-Math.log(i.r_min)) ** (1 / i.shape) : Infinity, Tuse = Math.min(Topt, Tsafe), sel = Number.isFinite(Tuse) ? ageCost(Tuse, i.shape, i.scale_h, i.cost_planned, i.cost_failure) : { rate: rtf, R: 0, up: mttf };
    if (Tsafe < Topt) warnings.push('The reliability floor governs: the interval is shorter than the pure cost optimum.');
    const cur = i.insp_interval_h > 0 ? ageCost(i.insp_interval_h, i.shape, i.scale_h, i.cost_planned, i.cost_failure) : null;
    warnings.push('Replacement and failure costs are illustrative placeholders; enter your own part, labour and disruption costs.');
    return {
      kpis: [
        { key: 'replace_interval_h', label: 'Recommended replacement age', value: Tuse, unit: 'h', note: Number.isFinite(Tuse) ? (Tsafe < Topt ? 'Set by the reliability floor' : 'Minimum cost rate') : 'No finite optimum: run to failure' },
        { key: 'cost_rate_usd_h', label: 'Cost rate at the recommended age', value: sel.rate, unit: 'USD/h' },
        { key: 'cost_rate_rtf_usd_h', label: 'Cost rate when run to failure', value: rtf, unit: 'USD/h' },
        { key: 'saving_pct', label: 'Saving against run to failure', value: 100 * (1 - sel.rate / rtf), unit: '%', status: sel.rate <= rtf ? 'ok' : 'warn' },
        { key: 'reliability_at_replacement', label: 'Reliability at the replacement age', value: Number.isFinite(Tuse) ? sel.R : 0, unit: '-' },
        { key: 'mean_life_used_pct', label: 'Share of mean life used', value: Number.isFinite(Tuse) ? (100 * sel.up) / mttf : 100, unit: '%', note: 'Mean operating time per part ÷ mean life' },
        { key: 'mttf_part_h', label: 'Mean life of the part', value: mttf, unit: 'h' },
        ...(cur ? [{ key: 'cost_rate_current_usd_h', label: 'Cost rate at the current interval', value: cur.rate, unit: 'USD/h', status: cur.rate <= 1.1 * sel.rate ? 'ok' : 'warn' }] : []),
      ],
      plots: [
        { type: 'line', title: 'Cost per operating hour versus replacement age', xlabel: 'Replacement age [h]', ylabel: 'Cost rate [USD/h]', xlog: true, series: [{ name: 'Age-replacement policy', x: Ts, y: cs }, { name: 'Run to failure', x: [Ts[0], Ts[Ts.length - 1]], y: [rtf, rtf], style: 'dash' }], annotations: [...(Number.isFinite(Tuse) ? [{ x: Tuse, label: 'Recommended' }] : []), ...(cur ? [{ x: i.insp_interval_h, label: 'Current' }] : [])] },
        { type: 'line', title: 'Reliability and hazard build-up with age', xlabel: 'Age [h]', ylabel: 'Reliability [-]', series: [{ name: 'R(t)', x: Ts, y: Ts.map((t) => weibR(t, i.shape, i.scale_h)) }], annotations: Number.isFinite(Tuse) ? [{ x: Tuse, label: 'Recommended' }] : [] },
      ],
      warnings, models: ['Age-replacement policy (renewal-reward cost rate)', 'Two-parameter Weibull life', 'Golden-section minimisation'],
      assumptions: ['Replacement restores the part to as-new', 'Costs are constant in real terms; no discounting within one part life', 'Failures are revealed immediately'],
    };
  },
  verify() {
    const sh = 2.5, sc = 1000, Cp = 1, Cf = 8, o = N.kv(maint.run({ shape: sh, scale_h: sc, cost_planned: Cp, cost_failure: Cf, r_min: 0, insp_interval_h: 0 })), T = o.replace_interval_h, a = ageCost(T, sh, sc, Cp, Cf), h = (sh / sc) * (T / sc) ** (sh - 1);
    const e = N.kv(maint.run({ shape: 1, scale_h: sc, cost_planned: Cp, cost_failure: Cf, r_min: 0, insp_interval_h: 0 }));
    return [
      N.check('Optimality condition h(T)·∫R dt − F(T) = Cp/(Cf − Cp)', h * a.up - (1 - a.R), Cp / (Cf - Cp), 1e-5, 'Barlow & Hunter age-replacement first-order condition'),
      N.check('Exponential life: run-to-failure cost rate Cf/MTTF', e.cost_rate_usd_h, Cf / sc, 1e-12, 'Memoryless property: no finite optimum'),
      N.check('Long-interval limit equals Cf/MTTF', ageCost(40 * sc, sh, sc, Cp, Cf).rate, Cf / (sc * N.gamma(1 + 1 / sh)), 1e-5, 'Renewal-reward theorem'),
    ];
  },
  calibration: { params: [{ key: 'shape', min: 0.5, max: 8 }, { key: 'scale_h', min: 1, max: 1e7 }], sweep: 'insp_interval_h', target: 'reliability_at_replacement', note: 'Supply the observed fraction of parts surviving to several removal ages to fit the Weibull life.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (Number.isFinite(o.replace_interval_h)) out.push({ severity: o.saving_pct > 10 ? 'advise' : 'info', title: `Replace at about ${o.replace_interval_h.toPrecision(3)} h`, detail: `Cost rate ${o.cost_rate_usd_h.toPrecision(3)} USD/h, ${o.saving_pct.toFixed(0)}% below run-to-failure, with ${(100 * o.reliability_at_replacement).toFixed(1)}% of parts reaching the limit and ${o.mean_life_used_pct.toFixed(0)}% of mean life used.`, action: 'Adopt this as the hard-time limit, or better, as the threshold for an on-condition inspection so that healthy parts stay in service (less material scrapped). Pass the cost rate to Suite 26.', basis: 'Minimum of the age-replacement cost rate' });
    else out.push({ severity: 'info', title: 'Run to failure or monitor on condition', detail: 'No finite replacement age lowers the cost rate with these inputs.', action: 'Do not replace on a fixed schedule; invest in condition monitoring where the failure effect is severe.', basis: 'Non-increasing hazard or cost ratio ≤ 1' });
    if (Number.isFinite(o.cost_rate_current_usd_h) && o.cost_rate_current_usd_h > 1.1 * o.cost_rate_usd_h) out.push({ severity: 'advise', title: 'The current interval is off the optimum', detail: `Current ${i.insp_interval_h} h costs ${o.cost_rate_current_usd_h.toPrecision(3)} USD/h against ${o.cost_rate_usd_h.toPrecision(3)} USD/h at the optimum.`, action: i.insp_interval_h < o.replace_interval_h ? 'The interval can be extended on cost grounds, provided the damage-tolerance inspection threshold of Suite 9 is still respected.' : 'Shorten the interval: failures are costing more than the replacements saved.', basis: 'Cost-rate comparison' });
    return out;
  },
};

export default {
  id: 'safety', n: 22,
  tagline: 'How likely is a failure, what causes it, how available is the aircraft, and which design or maintenance change buys the most safety.',
  analyses: [fta, fmeca, rbd, markov, eventTree, bayes, maint],
  consumes: [
    { from: 'fatigue', keys: ['life_fh', 'inspection_interval_fh'], why: 'Structural life and inspection interval for the wear-out model and replacement policy' },
    { from: 'hydmech', keys: ['bearing_L10_h'], why: 'Pump bearing life as a hydraulic failure-rate indicator' },
  ],
  provides: [
    { key: 'p_catastrophic_per_fh', label: 'Catastrophic top-event probability', unit: '1/FH' }, { key: 'system_reliability', label: 'Mission reliability', unit: '-' }, { key: 'mtbf_h', label: 'System mean time to failure', unit: 'h' },
    { key: 'availability', label: 'Steady-state availability', unit: '-' }, { key: 'dispatch_reliability', label: 'Dispatch reliability', unit: '-' },
  ],
  handoff: [
    { model: 'Certification-grade system safety assessment (FHA / PSSA / SSA, zonal and particular-risk analysis)', why: 'Requires the real system architecture, validated failure data, independence arguments and authority agreement; the default trees here are generic', tool: 'ARP4761-style safety process with dedicated FTA/FMEA tooling and qualified reliability data' },
    { model: 'Dynamic fault trees and Bayesian networks', why: 'Sequence-dependent gates (spares, priority-AND) and general conditional dependence need state-space or inference engines beyond static cut sets', tool: 'Dynamic FTA / Bayesian-network software' },
    { model: 'Physics-of-failure and stochastic degradation models', why: 'Need component-level stress, material and usage data; only Weibull life and the Suite 9 fatigue link are provided', tool: 'Prognostics and health-management tools fed by Suites 9 and 10' },
    { model: 'Censored-data life analysis (suspensions, interval data)', why: 'Median-rank regression here assumes complete failure data', tool: 'Maximum-likelihood Weibull analysis software' },
  ],
};
