// Suite 9 — Fatigue, Fracture and Damage Tolerance.
// Stress-life and strain-life initiation with mean-stress and notch corrections, rainflow counting with
// Palmgren–Miner damage on a flight spectrum, fatigue crack growth (Paris, Walker, Forman, simplified NASGRO)
// with residual strength and inspection intervals, probabilistic damage tolerance, and elastic–plastic
// fracture measures with a failure assessment diagram.

import * as N from '../core/numerics.js';
import { METALS, designAllowables } from '../data/materials.js';

// ---- shared helpers -------------------------------------------------------------------------
const MATS = Object.keys(METALS);
const mat = (name) => METALS[name] || METALS['Al 2024-T3'];
const MAT = { key: 'material', label: 'Material', type: 'select', options: MATS, default: 'Al 2024-T3', group: 'Material', help: 'Typical handbook fatigue and fracture constants, not design allowables' };
const kpi = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const DATA_NOTE = 'Fatigue and fracture constants are typical mean handbook values, not design allowables; scatter is covered only by the stated factors';
/** Safe-life scatter factor input. AC 23-13A para 2-26: 8 for a life from analysis with S-N data; 4 is the full-scale fatigue-test factor with two specimens (Table 2). */
const SCATTER = { key: 'scatter', label: 'Life scatter factor', unit: '-', default: 8, min: 1, max: 10, group: 'Method', help: 'Safe life = mean life / scatter factor. 8 for a life obtained by analysis from S–N data (AC 23-13A para 2-26); 4 only when the mean life is demonstrated by full-scale fatigue tests of two specimens (AC 23-13A Table 2: 4.96, 4.0, 3.7, 3.54 for one to four specimens). Transport-category practice (AC 25.571-1D) starts from a base factor of 3 on test results' };
const MEAN = ['Morrow', 'Smith–Watson–Topper', 'Goodman', 'Gerber', 'Soderberg', 'Walker', 'None'];
/**
 * 1 g gross stress [Pa] at the critical location, limit load factor and flight hours per flight from the case and upstream suites.
 * Without Suite 2 the limit-load stress is taken as 0.4 × the ultimate design strength, which is what its default cover sizing
 * produces (design stress 0.6·Ftu at ultimate load, ultimate factor 1.5).
 */
function usage(c, up) {
  const n = up.performance?.n_limit ?? c.aero.n_pos ?? 2.5, m = mat(c.struct.material), sLim = up.fea?.sigma_bend_root_Pa ?? up.fea?.sigma_max_Pa ?? 0.4 * designAllowables(m).Su;
  return { n, s1g: sLim / Math.max(n, 1), fh: c.econ.cycles_yr > 0 ? c.econ.util_fh_yr / c.econ.cycles_yr : 1, life_fl: c.econ.cycles_yr * c.econ.life_yr };
}
/** Equivalent fully reversed stress amplitude for amplitude sa and mean sm. */
export function eqAmplitude(m, sa, sm, method, gamma = 0.5) {
  const smax = sm + sa, pos = Math.max(sm, 0), lim = (d) => (d > 1e-9 ? sa / d : Infinity);
  switch (method) {
    case 'Goodman': return lim(1 - pos / m.Su);
    case 'Gerber': return lim(1 - (pos / m.Su) ** 2);
    case 'Soderberg': return lim(1 - pos / m.Sy);
    case 'Morrow': return lim(1 - sm / m.sf);
    case 'Smith–Watson–Topper': return smax > 0 ? Math.sqrt(smax * sa) : 0;
    case 'Walker': return smax > 0 ? smax ** (1 - gamma) * sa ** gamma : 0;
    default: return sa;
  }
}
/** Basquin life in cycles for a fully reversed amplitude: σa = σ'f (2N)^b. */
export const basquin = (m, sar) => (sar <= 0 ? Infinity : sar >= m.sf ? 0.5 : 0.5 * (sar / m.sf) ** (1 / m.b));
/** Peterson notch sensitivity q = 1/(1 + a/r) and fatigue notch factor. */
const peterson = (Kt, r_mm, a_mm) => { const q = r_mm > 0 ? 1 / (1 + a_mm / r_mm) : 1; return { q, Kf: 1 + q * (Kt - 1) }; };
/** Peterson characteristic length [mm]: steels from ultimate strength; 0.51 mm for aluminium alloys, 0.25 mm otherwise (empirical). */
const petersonA = (m) => (m.E > 150e9 ? 0.0254 * (2070e6 / m.Su) ** 1.8 : m.rho < 3000 ? 0.51 : 0.25);

// ---- 1. stress-life -------------------------------------------------------------------------
const SN_BASE = { material: 'Al 2024-T3', s_max: 120e6, R: 0, Kt: 1, r_mm: 3, a_pet_mm: 0.51, k_surf: 1, k_size: 1, k_env: 1, mean: 'Morrow', gamma: 0.5, scatter: 4, fh_per_flight: 1, cyc_per_flight: 1, target_flights: 60000 };
const sn = {
  id: 'sn', title: 'Stress-life (S–N) fatigue with mean-stress and notch corrections', fidelity: 'analytical',
  summary: 'Cycles to crack initiation under constant-amplitude loading from the Basquin curve, with Goodman, Gerber, Soderberg, Morrow, Smith–Watson–Topper or Walker mean-stress correction, Peterson notch sensitivity and surface, size and environment factors.',
  equations: ['Basquin fatigue life equation', 'Morrow mean-stress correction', 'Smith–Watson–Topper relation', 'Goodman relation', 'Gerber relation', 'Soderberg relation'],
  inputs: [
    MAT,
    { key: 's_max', label: 'Maximum nominal stress in the cycle', unit: 'Pa', default: 120e6, group: 'Loading', help: 'Gross-section stress at the once-per-flight peak; defaults to 1.3 × the 1 g stress from Suite 2' },
    { key: 'R', label: 'Stress ratio R = σmin/σmax', unit: '-', default: 0, min: -3, max: 0.95, group: 'Loading', help: '−1 fully reversed, 0 zero-to-tension; wing lower skin ground–air–ground is roughly −0.3 to 0' },
    { key: 'Kt', label: 'Elastic stress-concentration factor', unit: '-', default: 3, min: 1, max: 10, group: 'Notch', help: '3 for an open hole; use the plane-stress FE analysis in Suite 2' },
    { key: 'r_mm', label: 'Notch root radius', unit: 'mm', default: 3, min: 0.01, group: 'Notch' },
    { key: 'a_pet_mm', label: 'Peterson characteristic length', unit: 'mm', default: 0.51, min: 0.001, group: 'Notch', help: 'Empirical: about 0.5 mm for aluminium alloys, 0.06–0.25 mm for high-strength steels' },
    { key: 'k_surf', label: 'Surface-finish factor', unit: '-', default: 0.9, min: 0.2, max: 1.2, group: 'Modifying factors', help: 'Empirical: 1.0 polished, about 0.9 machined, 0.7 or less as-forged' },
    { key: 'k_size', label: 'Size factor', unit: '-', default: 1, min: 0.5, max: 1, group: 'Modifying factors' },
    { key: 'k_env', label: 'Environment (corrosion-fatigue) factor', unit: '-', default: 1, min: 0.2, max: 1, group: 'Modifying factors', help: 'Empirical knock-down on fatigue strength in a corrosive environment; 1 = laboratory air. Needs test data' },
    { key: 'mean', label: 'Mean-stress correction', type: 'select', options: MEAN, default: 'Morrow', group: 'Method' },
    { key: 'gamma', label: 'Walker exponent γ', unit: '-', default: 0.5, min: 0.1, max: 1, group: 'Method', help: '0.5 reproduces Smith–Watson–Topper; fit to S–N data at several R' },
    SCATTER,
    { key: 'cyc_per_flight', label: 'Cycles of this kind per flight', unit: '', default: 1, min: 0.01, group: 'Usage' },
    { key: 'fh_per_flight', label: 'Flight hours per flight', unit: 'h', default: 1.5, min: 0.01, group: 'Usage' },
    { key: 'target_flights', label: 'Design service goal', unit: 'flights', default: 60000, min: 1, group: 'Usage' },
  ],
  defaults(c, up) { const u = usage(c, up), m = mat(c.struct.material); return { material: c.struct.material, s_max: 1.3 * u.s1g, a_pet_mm: petersonA(m), fh_per_flight: u.fh, target_flights: u.life_fl || undefined }; },
  run(i) {
    const m = mat(i.material), { q, Kf } = peterson(i.Kt, i.r_mm, i.a_pet_mm), kk = i.k_surf * i.k_size * i.k_env, warnings = [];
    const sa = (0.5 * i.s_max * (1 - i.R) * Kf) / kk, sm = 0.5 * i.s_max * (1 + i.R), sar = eqAmplitude(m, sa, sm, i.mean, i.gamma), Nf = basquin(m, sar), flights = Nf / i.cyc_per_flight, safe = flights / i.scatter;
    const sTarget = m.sf * (2 * i.target_flights * i.cyc_per_flight * i.scatter) ** m.b, fos = sar > 0 ? sTarget / sar : Infinity;
    if (i.Kt * Math.abs(i.s_max) > m.Sy) warnings.push('The elastic notch stress exceeds yield: the stress-life method misses local plasticity and mean-stress relaxation. Use the strain-life analysis.');
    if (Nf < 1e4) warnings.push('Predicted life is in the low-cycle regime (< 10⁴ cycles), where strain-life is the appropriate method.');
    if (sm + sa >= m.Su) warnings.push('Maximum local stress reaches the ultimate strength: static failure governs.');
    if (i.k_env < 1) warnings.push('The environment factor is an empirical knock-down; corrosion-fatigue interaction is frequency- and time-dependent and must be confirmed by test.');
    const Ns = N.logspace(1e2, 1e9, 60), sms = N.linspace(0, 0.95 * m.Su, 40), sRef = sTarget;
    const haigh = (meth) => sms.map((x) => { let lo = 0, hi = m.Su; for (let k = 0; k < 50; k++) { const mid = 0.5 * (lo + hi); if (eqAmplitude(m, mid, x, meth, i.gamma) < sRef) lo = mid; else hi = mid; } return lo / 1e6; });
    return {
      kpis: [
        kpi('sn_life_cycles', 'Mean cycles to crack initiation', Nf, 'cycles'),
        kpi('sn_life_flights', 'Mean life', flights, 'flights'),
        kpi('sn_safe_life_flights', 'Safe life (mean / scatter factor)', safe, 'flights', safe >= i.target_flights ? 'ok' : safe >= 0.5 * i.target_flights ? 'warn' : 'bad', `Design service goal ${i.target_flights.toFixed(0)} flights`),
        kpi('sn_safe_life_fh', 'Safe life in flight hours', safe * i.fh_per_flight, 'h'),
        kpi('sigma_ar_Pa', 'Equivalent fully reversed amplitude', sar, 'Pa'),
        kpi('sigma_a_local_Pa', 'Local stress amplitude (Kf and factors applied)', sa, 'Pa'), kpi('sigma_mean_Pa', 'Mean stress', sm, 'Pa'),
        kpi('Kf', 'Fatigue notch factor', Kf, '-'), kpi('q_notch', 'Notch sensitivity', q, '-'),
        kpi('stress_reserve', 'Stress reserve factor at the design goal', fos, '-', fos >= 1 ? 'ok' : 'bad', 'Allowable equivalent amplitude at goal × scatter / applied'),
        kpi('S_1e7_Pa', 'Fatigue strength at 10⁷ cycles (R = −1)', m.sf * (2e7) ** m.b, 'Pa'),
      ],
      plots: [
        { type: 'line', title: 'S–N curve (fully reversed, smooth)', xlabel: 'Cycles to failure [-]', ylabel: 'Stress amplitude [MPa]', xlog: true, series: [{ name: 'Basquin mean curve', x: Ns, y: Ns.map((n) => (m.sf * (2 * n) ** m.b) / 1e6) }, { name: 'Operating point (equivalent amplitude)', x: [Math.max(Nf, 1)], y: [sar / 1e6], style: 'points' }], annotations: [{ x: i.target_flights * i.cyc_per_flight * i.scatter, label: 'Goal × scatter' }] },
        { type: 'line', title: 'Constant-life (Haigh) diagram at the factored design goal', xlabel: 'Mean stress [MPa]', ylabel: 'Allowable amplitude [MPa]', series: [...['Goodman', 'Gerber', 'Soderberg', 'Morrow', 'Smith–Watson–Topper'].map((meth) => ({ name: meth, x: sms.map((v) => v / 1e6), y: haigh(meth) })), { name: 'Operating point', x: [sm / 1e6], y: [sa / 1e6], style: 'points' }] },
      ],
      warnings,
      models: ['Basquin stress-life curve', `${i.mean} mean-stress correction`, 'Peterson notch sensitivity (empirical)', 'Surface, size and environment factors (empirical)'],
      assumptions: ['Constant-amplitude loading; life to initiation of an engineering-size crack', 'Notch factor and modifying factors applied to the amplitude only; mean stress taken as nominal', 'No endurance limit is assumed (appropriate for aluminium; conservative for steel and titanium below their limit)', 'Peterson characteristic length and the surface, size and environment factors are typical empirical values, not sourced data', DATA_NOTE],
    };
  },
  verify() {
    const m = mat('Al 2024-T3'), r1 = N.kv(sn.run({ ...SN_BASE, s_max: 200e6, R: -1, mean: 'None' })), r2 = N.kv(sn.run({ ...SN_BASE, s_max: 200e6, R: 0, mean: 'Smith–Watson–Topper' })), r3 = N.kv(sn.run({ ...SN_BASE, s_max: 200e6, R: 0, mean: 'Goodman' }));
    return [
      N.check('Basquin: σ′f·(2N)^b returns the applied amplitude', m.sf * (2 * r1.sn_life_cycles) ** m.b, 200e6, 1e-10, 'Basquin (1910)'),
      N.check('SWT at R = 0: σar = σmax/√2', r2.sigma_ar_Pa, 200e6 / Math.SQRT2, 1e-12, 'Smith, Watson & Topper (1970)'),
      N.check('Goodman: σar = σa/(1 − σm/Su)', r3.sigma_ar_Pa, 100e6 / (1 - 100e6 / m.Su), 1e-12, 'Goodman relation'),
      N.check('Walker with γ = 0.5 equals SWT', eqAmplitude(m, 80e6, 40e6, 'Walker', 0.5), eqAmplitude(m, 80e6, 40e6, 'Smith–Watson–Topper'), 1e-12, 'Walker (1970)'),
      N.check('Peterson: blunt notch is fully sensitive (Kf → Kt)', peterson(3, 1e6, 0.5).Kf, 3, 1e-6, 'Peterson (1959)'),
    ];
  },
  calibration: { params: [{ key: 'gamma', min: 0.1, max: 1 }, { key: 'k_surf', min: 0.2, max: 1.2 }], sweep: 's_max', target: 'sn_life_cycles', note: 'Coupon S–N data (life versus maximum stress at one or more stress ratios) calibrate the Walker exponent and surface factor.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.sn_safe_life_flights < i.target_flights) out.push({ severity: o.stress_reserve < 0.8 ? 'critical' : 'warn', title: 'Safe life is below the design service goal', detail: `${o.sn_safe_life_flights.toExponential(2)} flights against ${i.target_flights.toFixed(0)}; the stress must fall by ${(100 * (1 - o.stress_reserve)).toFixed(0)}%.`, action: `Lower Kt (now ${i.Kt}) with a larger radius or bushing, cold-work or interference-fit the hole, improve surface finish, or reduce the 1 g stress. Life scales with stress to the power ${(-1 / mat(i.material).b).toFixed(0)}.`, basis: 'Basquin curve with scatter factor' });
    else out.push({ severity: 'info', title: 'Life extension potential', detail: `Safe life exceeds the goal by a factor ${(o.sn_safe_life_flights / i.target_flights).toFixed(1)}.`, action: 'A documented margin supports a later life-extension programme: keeping the airframe in service longer amortises its embodied energy and material over more flights.', basis: 'Safe-life margin' });
    if (i.k_env < 1 || i.k_surf < 0.85) out.push({ severity: 'advise', title: 'Surface and environment drive this result', detail: `Combined factor ${(i.k_surf * i.k_size * i.k_env).toFixed(2)}.`, action: 'Protective finishes, sealing and shot-peening are low-mass ways to recover fatigue strength.', basis: 'Empirical modifying factors' });
    return out;
  },
};

// ---- 2. strain-life with Neuber notch plasticity --------------------------------------------
/** Cyclic Ramberg–Osgood constants from strain-life compatibility: n' = b/c, K' = σ'f/ε'f^n'. */
const cyclic = (m) => { const n = m.b / m.c; return { n, K: m.sf / m.ef ** n }; };
/** Neuber: find local stress on σ·ε(σ) = (Kt·S)²/E, with ε from the Ramberg–Osgood curve scaled by f (1 for loading, 2 for ranges). */
function neuber(m, KtS, f) {
  const { n, K } = cyclic(m), tgt = (KtS * KtS) / m.E, g = (s) => s * (s / m.E + f * (s / (f * K)) ** (1 / n)) - tgt;
  if (!(KtS > 0)) return { s: 0, e: 0 };
  const s = N.brent(g, 0, KtS, 1e-9 * KtS + 1e-6);
  return { s, e: s / m.E + f * (s / (f * K)) ** (1 / n) };
}
/** Reversals to failure 2N from strain amplitude, mean and maximum stress. */
function strainLife(m, ea, sm, smax, method) {
  const f = (x) => { const r = 10 ** x; return method === 'Smith–Watson–Topper' ? (m.sf ** 2 / m.E) * r ** (2 * m.b) + m.sf * m.ef * r ** (m.b + m.c) - Math.max(smax, 0) * ea : (Math.max(m.sf - sm, 1) / m.E) * r ** m.b + m.ef * r ** m.c - ea; };
  if (method === 'Smith–Watson–Topper' && smax <= 0) return Infinity;
  if (f(0) <= 0) return 1; if (f(14) >= 0) return Infinity;
  return 10 ** N.brent(f, 0, 14, 1e-10);
}
const en = {
  id: 'en', title: 'Strain-life (ε–N) fatigue with Neuber notch plasticity', fidelity: 'numerical',
  summary: 'Local stress–strain response at a notch from Neuber’s rule on the cyclic stress–strain curve, the stabilised hysteresis loop, and cycles to initiation from the Coffin–Manson–Basquin curve with Morrow or Smith–Watson–Topper mean-stress correction.',
  equations: ['Coffin–Manson equation', 'Basquin fatigue life equation', 'Morrow mean-stress correction', 'Smith–Watson–Topper relation'],
  inputs: [
    MAT,
    { key: 's_max', label: 'Maximum nominal stress', unit: 'Pa', default: 120e6, group: 'Loading' },
    { key: 'R', label: 'Stress ratio', unit: '-', default: 0, min: -3, max: 0.95, group: 'Loading' },
    { key: 'Kt', label: 'Stress-concentration factor', unit: '-', default: 3, min: 1, max: 10, group: 'Notch' },
    { key: 'mean', label: 'Mean-stress correction', type: 'select', options: ['Morrow', 'Smith–Watson–Topper'], default: 'Smith–Watson–Topper', group: 'Method' },
    SCATTER,
    { key: 'cyc_per_flight', label: 'Cycles of this kind per flight', unit: '', default: 1, min: 0.01, group: 'Usage' },
    { key: 'target_flights', label: 'Design service goal', unit: 'flights', default: 60000, min: 1, group: 'Usage' },
  ],
  defaults(c, up) { const u = usage(c, up); return { material: c.struct.material, s_max: 1.3 * u.s1g, target_flights: u.life_fl || undefined }; },
  run(i) {
    const m = mat(i.material), cy = cyclic(m), warnings = [], Smax = Math.abs(i.s_max), dS = Smax * (1 - i.R);
    const pk = neuber(m, i.Kt * Smax, 1), rg = neuber(m, i.Kt * dS, 2), sig = Math.sign(i.s_max) || 1, smaxL = sig * pk.s, sminL = smaxL - sig * rg.s, sm = 0.5 * (smaxL + sminL), ea = rg.e / 2;
    const rev = strainLife(m, ea, sm, Math.max(smaxL, sminL), i.mean), Nf = rev / 2, flights = Nf / i.cyc_per_flight, safe = flights / i.scatter, eap = ea - rg.s / (2 * m.E), rev_t = ((m.ef * m.E) / m.sf) ** (1 / (m.b - m.c));
    if (Smax > m.Sy) warnings.push('Nominal stress exceeds yield: Neuber’s rule assumes net-section elasticity and is no longer valid.');
    if (eap / ea > 0.5) warnings.push('Plastic strain dominates the cycle (low-cycle fatigue): expect a life of only hundreds to thousands of cycles.');
    const revs = N.logspace(10, 1e9, 70), up = N.linspace(0, rg.s, 30), loopS = [], loopE = [];
    for (const s of up) { loopS.push((sminL + sig * s) / 1e6); loopE.push(100 * (sig * (pk.e - rg.e) + sig * (s / m.E + 2 * (s / (2 * cy.K)) ** (1 / cy.n)))); }
    for (const s of up) { loopS.push((smaxL - sig * s) / 1e6); loopE.push(100 * (sig * pk.e - sig * (s / m.E + 2 * (s / (2 * cy.K)) ** (1 / cy.n)))); }
    return {
      kpis: [
        kpi('en_life_cycles', 'Mean cycles to crack initiation', Nf, 'cycles'),
        kpi('en_safe_life_flights', 'Safe life (mean / scatter factor)', safe, 'flights', safe >= i.target_flights ? 'ok' : safe >= 0.5 * i.target_flights ? 'warn' : 'bad', `Design service goal ${i.target_flights.toFixed(0)} flights`),
        kpi('eps_a', 'Local strain amplitude', ea, '-'), kpi('eps_a_plastic', 'Plastic strain amplitude', eap, '-'),
        kpi('sigma_max_local_Pa', 'Local maximum stress', smaxL, 'Pa', undefined, `Elastic estimate Kt·S = ${(i.Kt * i.s_max / 1e6).toFixed(0)} MPa`),
        kpi('sigma_mean_local_Pa', 'Local mean stress', sm, 'Pa', undefined, 'Notch yielding relaxes the mean stress below the elastic value'),
        kpi('transition_reversals', 'Transition life 2Nt', rev_t, 'reversals', undefined, 'Elastic and plastic strain amplitudes are equal'),
        kpi('K_cyclic_Pa', 'Cyclic strength coefficient K′', cy.K, 'Pa'), kpi('n_cyclic', 'Cyclic hardening exponent n′', cy.n, '-'),
      ],
      plots: [
        { type: 'line', title: 'Strain-life curve', xlabel: 'Reversals to failure 2N [-]', ylabel: 'Strain amplitude [-]', xlog: true, ylog: true, series: [{ name: 'Total', x: revs, y: revs.map((r) => (m.sf / m.E) * r ** m.b + m.ef * r ** m.c) }, { name: 'Elastic (Basquin)', x: revs, y: revs.map((r) => (m.sf / m.E) * r ** m.b), style: 'dash' }, { name: 'Plastic (Coffin–Manson)', x: revs, y: revs.map((r) => m.ef * r ** m.c), style: 'dash' }, { name: 'Operating point', x: [Number.isFinite(rev) ? rev : 1e9], y: [ea], style: 'points' }] },
        { type: 'line', title: 'Stabilised local hysteresis loop at the notch', xlabel: 'Local strain [%]', ylabel: 'Local stress [MPa]', series: [{ name: 'Hysteresis loop (Masing)', x: loopE, y: loopS }] },
      ],
      warnings,
      models: ['Coffin–Manson–Basquin strain-life curve', 'Ramberg–Osgood cyclic curve with constants from compatibility (n′ = b/c)', 'Neuber’s rule with Masing behaviour', `${i.mean} mean-stress correction`],
      assumptions: ['Stabilised cyclic response from the first cycle; no cyclic creep or sequence memory beyond one loop', 'Kt used directly in Neuber’s rule (conservative versus Kf)', 'Uniaxial local stress state', DATA_NOTE],
    };
  },
  verify() {
    const m = mat('Al 7075-T6'), rt = ((m.ef * m.E) / m.sf) ** (1 / (m.b - m.c)), nb = neuber(m, 600e6, 1), b = { material: 'Al 7075-T6', s_max: 60e6, R: -1, Kt: 2, mean: 'Morrow', scatter: 1, cyc_per_flight: 1, target_flights: 1 }, lo = N.kv(en.run(b));
    return [
      N.check('Transition life: elastic strain equals plastic strain', (m.sf / m.E) * rt ** m.b, m.ef * rt ** m.c, 1e-10, 'Coffin–Manson–Basquin'),
      N.check('Neuber: local σ·ε = (Kt·S)²/E', nb.s * nb.e, 600e6 ** 2 / m.E, 1e-7, 'Neuber (1961)'),
      N.check('Elastic limit: strain-life life tends to the Basquin life at Kt·S', lo.en_life_cycles, basquin(m, 120e6), 0.02, 'High-cycle limit (plastic strain is negligible but not zero)'),
      N.check('Strain-life curve inverted consistently', (m.sf / m.E) * 1e5 ** m.b + m.ef * 1e5 ** m.c, ((x) => (m.sf / m.E) * x ** m.b + m.ef * x ** m.c)(strainLife(m, (m.sf / m.E) * 1e5 ** m.b + m.ef * 1e5 ** m.c, 0, 0, 'Morrow')), 1e-8, 'Root-finding check'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.en_safe_life_flights < i.target_flights) out.push({ severity: 'warn', title: 'Notch initiation life is below the goal', detail: `Safe life ${o.en_safe_life_flights.toExponential(2)} flights; local strain amplitude ${(100 * o.eps_a).toFixed(2)}%.`, action: 'Reduce Kt or nominal stress; compressive residual stress from cold-working lowers the local mean stress and is very effective at this strain level.', basis: 'Strain-life with Neuber correction' });
    if (o.eps_a_plastic > 0.3 * o.eps_a) out.push({ severity: 'advise', title: 'Significant cyclic plasticity at the notch', detail: `${(100 * o.eps_a_plastic / o.eps_a).toFixed(0)}% of the strain amplitude is plastic.`, action: 'Treat this location as fatigue-critical: plan a damage-tolerance (crack-growth) assessment and an inspection.', basis: 'Low-cycle fatigue regime' });
    return out;
  },
};

// ---- 3. spectrum, rainflow and Miner --------------------------------------------------------
/** ASTM E1049 rainflow counting. Returns [[range, mean, count], ...] with count 1 (full cycle) or 0.5 (half cycle). */
export function rainflow(series) {
  const tp = [];
  for (const v of series) { const n = tp.length; if (n && v === tp[n - 1]) continue; if (n >= 2 && (tp[n - 1] - tp[n - 2]) * (v - tp[n - 1]) > 0) tp[n - 1] = v; else tp.push(v); }
  const st = [], out = [];
  for (const p of tp) {
    st.push(p);
    while (st.length >= 3) {
      const n = st.length, X = Math.abs(st[n - 1] - st[n - 2]), Y = Math.abs(st[n - 2] - st[n - 3]);
      if (X < Y) break;
      if (n === 3) { out.push([Y, 0.5 * (st[0] + st[1]), 0.5]); st.shift(); } // range Y contains the starting point: half cycle
      else { out.push([Y, 0.5 * (st[n - 2] + st[n - 3]), 1]); const last = st.pop(); st.pop(); st.pop(); st.push(last); }
    }
  }
  for (let k = 0; k < st.length - 1; k++) out.push([Math.abs(st[k + 1] - st[k]), 0.5 * (st[k] + st[k + 1]), 0.5]);
  return out;
}
/** Illustrative flight-by-flight load-factor sequence: ground, then gust/manoeuvre pairs with exponentially distributed increments, then ground. */
function flightSequence(i, nFl) {
  const u = N.rng(Math.round(i.seed)), seq = [], per = Math.max(0, i.N0_per_fh * i.fh_per_flight);
  for (let f = 0; f < nFl; f++) {
    seq.push(i.n_ground);
    let k = Math.floor(per); if (u() < per - k) k++;
    for (let j = 0; j < k; j++) { const dn = Math.min(i.dn_max, i.dn_min - i.dn_scale * Math.log(1 - u() * (1 - 1e-12))); seq.push(1 + dn, 1 - dn); }
    seq.push(1);
  }
  seq.push(i.n_ground);
  return seq;
}
const spectrum = {
  id: 'spectrum', title: 'Flight spectrum, rainflow counting and Miner damage', fidelity: 'numerical',
  summary: 'Builds a flight-by-flight load sequence (ground–air–ground plus gust and manoeuvre cycles), counts cycles by the rainflow method and sums Palmgren–Miner damage to give life in flights and flight hours.',
  equations: ['Palmgren–Miner cumulative damage rule', 'Basquin fatigue life equation', 'Morrow mean-stress correction', 'Smith–Watson–Topper relation', 'Goodman relation'],
  inputs: [
    MAT,
    { key: 's_1g', label: 'Gross stress at 1 g', unit: 'Pa', default: 80e6, min: 1e5, group: 'Loading', help: 'Limit-load stress from Suite 2 divided by the limit load factor' },
    { key: 'Kt', label: 'Stress-concentration factor', unit: '-', default: 3, min: 1, max: 10, group: 'Loading' },
    { key: 'n_ground', label: 'Ground load factor equivalent', unit: 'g', default: -0.3, min: -1.5, max: 1, group: 'Spectrum (illustrative)', help: 'Stress on the ground as a fraction of the 1 g flight stress; negative for a wing lower surface' },
    { key: 'N0_per_fh', label: 'Gust/manoeuvre cycles per flight hour', unit: '1/h', default: 12, min: 0, max: 2000, group: 'Spectrum (illustrative)', help: 'Cycles exceeding the smallest counted increment. Replace with a measured or regulatory exceedance spectrum' },
    { key: 'dn_min', label: 'Smallest counted load-factor increment', unit: 'g', default: 0.05, min: 0, max: 1, group: 'Spectrum (illustrative)' },
    { key: 'dn_scale', label: 'Exceedance decay scale', unit: 'g', default: 0.12, min: 0.005, max: 2, group: 'Spectrum (illustrative)', help: 'Exceedances fall by e for each such increment: N(Δn) = N₀·exp(−(Δn − Δn_min)/scale)' },
    { key: 'dn_max', label: 'Truncation increment', unit: 'g', default: 1.5, min: 0.05, max: 10, group: 'Spectrum (illustrative)', help: 'Limit load factor minus 1 is a natural cap' },
    { key: 'fh_per_flight', label: 'Flight hours per flight', unit: 'h', default: 1.5, min: 0.01, group: 'Usage' },
    { key: 'target_flights', label: 'Design service goal', unit: 'flights', default: 60000, min: 1, group: 'Usage' },
    { key: 'mean', label: 'Mean-stress correction', type: 'select', options: MEAN, default: 'Smith–Watson–Topper', group: 'Method' },
    SCATTER,
    { key: 'seed', label: 'Random seed', unit: '', default: 2024, min: 1, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nFlights', label: 'Flights in the simulated block', unit: '', default: 400, min: 10, max: 20000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c, up) { const u = usage(c, up); return { material: c.struct.material, s_1g: u.s1g, dn_max: Math.max(0.2, u.n - 1), fh_per_flight: u.fh, target_flights: u.life_fl || undefined, N0_per_fh: c.meta.type === 'helicopter' ? 60 : undefined }; },
  run(i) {
    const m = mat(i.material), nFl = Math.max(10, Math.round(i.nFlights)), seq = flightSequence(i, nFl), cyc = rainflow(seq), k = i.Kt * i.s_1g, warnings = [];
    let D = 0, Dgag = 0, nCyc = 0; const edges = N.linspace(0, N.amax(cyc.map((c) => c[0])) * 1.0001, 13), dBin = new Array(12).fill(0), cBin = new Array(12).fill(0);
    const mb = N.linspace(N.amin(cyc.map((c) => c[1])) - 1e-9, N.amax(cyc.map((c) => c[1])) + 1e-9, 9), z = N.range(12, () => new Array(8).fill(0));
    for (const [rg, mn, cnt] of cyc) {
      const Nf = basquin(m, eqAmplitude(m, 0.5 * rg * k, mn * k, i.mean)), d = cnt / Nf, b = Math.min(11, Math.floor((rg / edges[12]) * 12)), bm = Math.min(7, Math.max(0, Math.floor(((mn - mb[0]) / (mb[8] - mb[0])) * 8)));
      D += d; nCyc += cnt; dBin[b] += d; cBin[b] += cnt; z[b][bm] += cnt; if (rg > 1 - i.n_ground - 1e-9) Dgag += d;
    }
    const dpf = D / nFl, life = dpf > 0 ? 1 / dpf : Infinity, safe = life / i.scatter, smaxLoc = k * N.amax(seq);
    if (smaxLoc > m.Sy) warnings.push('The largest local elastic stress in the spectrum exceeds yield: linear damage summation ignores the beneficial or harmful residual stresses left by such overloads.');
    warnings.push('The default spectrum is illustrative (exponential exceedance law with a random sequence). Use a measured, regulatory or mission-analysis spectrum for any substantiation.');
    const ex = N.linspace(i.dn_min, i.dn_max, 40), nShow = Math.min(seq.length, 300), ctr = (e) => N.range(e.length - 1, (j) => 0.5 * (e[j] + e[j + 1]));
    return {
      kpis: [
        kpi('damage_per_flight', 'Miner damage per flight', dpf, '1/flight'),
        kpi('life_cycles', 'Mean fatigue life', life, 'flights'),
        kpi('life_fh', 'Mean fatigue life in flight hours', life * i.fh_per_flight, 'h'),
        kpi('safe_life_flights', 'Safe life (mean / scatter factor)', safe, 'flights', safe >= i.target_flights ? 'ok' : safe >= 0.5 * i.target_flights ? 'warn' : 'bad', `Design service goal ${i.target_flights.toFixed(0)} flights`),
        kpi('safe_life_fh', 'Safe life in flight hours', safe * i.fh_per_flight, 'h'),
        kpi('damage_at_goal', 'Miner sum at the design goal (with scatter)', dpf * i.target_flights * i.scatter, '-', dpf * i.target_flights * i.scatter <= 1 ? 'ok' : 'bad', 'Must not exceed 1'),
        kpi('gag_damage_frac', 'Share of damage from ground–air–ground cycles', D > 0 ? Dgag / D : 0, '-'),
        kpi('cycles_per_flight', 'Counted cycles per flight', nCyc / nFl, '-'),
        kpi('sigma_local_max_Pa', 'Largest local elastic stress in the block', smaxLoc, 'Pa'),
      ],
      plots: [
        { type: 'line', title: 'Load-factor sequence (start of the block)', xlabel: 'Turning point [-]', ylabel: 'Load factor [g]', series: [{ name: 'Sequence', x: N.range(nShow), y: seq.slice(0, nShow) }] },
        { type: 'heat', title: 'Rainflow matrix (cycle counts in the block)', xlabel: 'Cycle mean [g]', ylabel: 'Cycle range [g]', zlabel: 'log10(1 + cycles)', x: ctr(mb), y: ctr(edges), z: z.map((r) => r.map((v) => Math.log10(1 + v))) },
        { type: 'bar', title: 'Damage by cycle range', ylabel: 'Share of total damage [%]', categories: ctr(edges).map((v) => v.toFixed(2) + ' g'), series: [{ name: 'Damage', y: dBin.map((v) => (D > 0 ? (100 * v) / D : 0)) }] },
        { type: 'line', title: 'Assumed exceedance spectrum', xlabel: 'Load-factor increment Δn [g]', ylabel: 'Exceedances per flight hour [1/h]', ylog: true, series: [{ name: 'N(Δn)', x: ex, y: ex.map((d) => Math.max(1e-12, i.N0_per_fh * Math.exp(-(d - i.dn_min) / i.dn_scale))) }] },
      ],
      warnings,
      models: ['Seeded random flight-by-flight sequence from an exponential exceedance law (illustrative)', 'ASTM E1049 rainflow counting', `Basquin S–N curve with ${i.mean} mean-stress correction`, 'Palmgren–Miner linear damage'],
      assumptions: ['Stress proportional to load factor; Kt applied to amplitude and mean (local elastic stress)', 'Linear damage accumulation with no load-sequence or overload-retardation effect', 'Failure at Miner sum 1 on the mean curve; scatter handled by the life factor', 'Spectrum parameters (cycles per hour, decay scale, ground load) are illustrative, not sourced: AC 23-13A Appendix 1 or measured exceedance data should replace them', DATA_NOTE],
    };
  },
  convergence: { param: 'nFlights', label: 'Flights in the block', levels: [50, 100, 200, 400, 800], metric: 'damage_per_flight', hOf: (n) => 1 / Math.sqrt(n) }, // sampling study: error falls with 1/√flights
  verify() {
    const rf = rainflow([-2, 1, -3, 5, -1, 3, -4, 4, -2]), cnt = (r) => N.sum(rf.filter((c) => c[0] === r).map((c) => c[2])), m = mat('Al 2024-T3');
    const ca = N.kv(spectrum.run({ material: 'Al 2024-T3', s_1g: 100e6, Kt: 1, n_ground: 0, N0_per_fh: 0, dn_min: 0, dn_scale: 0.1, dn_max: 1, fh_per_flight: 1, target_flights: 1, mean: 'Smith–Watson–Topper', scatter: 1, seed: 1, nFlights: 50 }));
    return [
      N.check('ASTM E1049 example: 0.5 cycle of range 3', cnt(3), 0.5, 0, 'ASTM E1049-85 §5.4.4 rainflow example'),
      N.check('ASTM E1049 example: 1.5 cycles of range 4', cnt(4), 1.5, 0, 'ASTM E1049-85'),
      N.check('ASTM E1049 example: 0.5 cycle of range 6', cnt(6), 0.5, 0, 'ASTM E1049-85'),
      N.check('ASTM E1049 example: 1 cycle of range 8', cnt(8), 1, 0, 'ASTM E1049-85'),
      N.check('ASTM E1049 example: 0.5 cycle of range 9', cnt(9), 0.5, 0, 'ASTM E1049-85'),
      N.check('Constant-amplitude 0–1 g block: damage per flight = 1/N', ca.damage_per_flight, 1 / basquin(m, 100e6 / Math.SQRT2), 1e-9, 'Palmgren–Miner with one cycle per flight'),
    ];
  },
  calibration: { params: [{ key: 'dn_scale', min: 0.005, max: 2 }, { key: 'N0_per_fh', min: 0, max: 2000 }], sweep: 'dn_max', target: 'cycles_per_flight', note: 'Recorded exceedance counts from flight-loads monitoring calibrate the spectrum parameters.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    // below the goal on the mean curve the detail cracks in service; above it, only the scatter-factored safe life is unsubstantiated
    const meanShort = o.life_cycles < i.target_flights;
    if (o.damage_at_goal > 1) out.push({ severity: meanShort ? 'critical' : 'warn', title: meanShort ? 'Mean fatigue life is below the design service goal' : 'Safe life with the scatter factor is below the design service goal', detail: `Miner sum ${o.damage_at_goal.toFixed(2)} at the goal with scatter factor ${i.scatter}; mean life ${o.life_cycles.toFixed(0)} flights, safe life ${o.safe_life_flights.toFixed(0)} flights.`, action: o.gag_damage_frac > 0.5 ? 'Ground–air–ground cycles dominate: reduce the 1 g stress level or Kt at this detail.' : 'Gust and manoeuvre cycles dominate: gust-load alleviation (Suite 16) or a lower wing loading sensitivity reduces damage every flight hour.' + (meanShort ? '' : ' Alternatively substantiate the life by full-scale fatigue test (lower scatter factor) or by damage tolerance with inspections (crack-growth analysis).'), basis: 'Palmgren–Miner sum ≤ 1 with the scatter factor (AC 23-13A: 8 for analysis by S–N data)' });
    else out.push({ severity: 'info', title: 'Life margin and extension', detail: `Safe life ${o.safe_life_flights.toFixed(0)} flights (${o.safe_life_fh.toFixed(0)} h) against a goal of ${i.target_flights.toFixed(0)}.`, action: 'Usage monitoring of actual exceedances lets the fleet consume real rather than assumed damage, typically extending service life and deferring replacement.', basis: 'Individual aircraft tracking' });
    out.push({ severity: 'advise', title: 'Replace the illustrative spectrum', detail: `${o.cycles_per_flight.toFixed(1)} counted cycles per flight from an assumed exponential exceedance law.`, action: 'Upload measured or regulatory exceedance data and calibrate; fatigue life is very sensitive to the tail of the spectrum.', basis: 'Spectrum sensitivity' });
    return out;
  },
};

// ---- 4. crack growth, residual strength and inspection --------------------------------------
const GEOMS = ['centre crack', 'edge crack', 'crack from hole (one side)', 'cracks from hole (both sides)'];
const sec = (x) => 1 / Math.cos(Math.min(x, 1.55));
/** Geometry factor β in K = β·σ·√(π·a). a: half-length (centre), depth (edge) or length from the hole edge; W: full width; r: hole radius. */
export function beta(geom, a, W, r) {
  if (geom === 'edge crack') { const x = Math.min(a / W, 0.7); return 1.122 - 0.231 * x + 10.55 * x * x - 21.71 * x ** 3 + 30.382 * x ** 4; }
  if (geom === 'crack from hole (one side)') return (0.6762 + 0.8733 / (0.3245 + a / r)) * Math.sqrt(sec((Math.PI * (2 * r + a)) / (2 * W)));
  if (geom === 'cracks from hole (both sides)') return (0.9439 + 0.6865 / (0.2772 + a / r)) * Math.sqrt(sec((Math.PI * (r + a)) / W));
  return Math.sqrt(sec((Math.PI * a) / W));
}
/** Largest admissible crack size for the geometry (95% of the ligament, with a further 10% taken off at a hole where the curve fit degrades). */
const aLimit = (geom, W, r) => 0.95 * (geom === 'edge crack' ? 0.7 * W : geom.includes('hole') ? 0.9 * (W / 2 - r) : W / 2); // the hole is central: either crack runs out of plate at W/2 − r
/** Crack growth rate [m/cycle]; ΔK, Kc, ΔKth in MPa√m. */
export function growthRate(law, dK, R, p) {
  if (!(dK > 0) || dK <= p.dKth * (law === 'Paris' ? 0 : 1)) return 0;
  const Kmax = dK / (1 - R);
  if (law === 'Walker') return p.C * (dK * (1 - R) ** (p.gamma - 1)) ** p.m;
  if (law === 'Forman') { const d = (1 - R) * p.Kc - dK; return d > 0 ? (p.C * p.Kc * dK ** p.m) / d : Infinity; }
  if (law === 'NASGRO (simplified)') return Kmax >= p.Kc ? Infinity : (p.C * dK ** p.m * (1 - p.dKth / dK) ** p.pn) / (1 - Kmax / p.Kc) ** p.qn;
  return p.C * dK ** p.m;
}
/** Integrate cycles from a0 to a1 on a geometric grid of n intervals (Simpson in ln a). Returns the a(N) history. */
function growCrack(i, p, a0, a1, n) {
  const R = Math.max(i.R, 0), ds = i.s_max * (1 - R) / 1e6, f = (a) => { const r = growthRate(i.law, beta(i.geom, a, i.W, i.r_hole) * ds * Math.sqrt(Math.PI * a), R, p); return r > 0 ? (Number.isFinite(r) ? a / r : 0) : Infinity; };
  const as = [a0], Ns = [0], h = Math.log(a1 / a0) / n; let tot = 0;
  for (let k = 0; k < n; k++) { const x0 = Math.log(a0) + k * h, s = (h / 6) * (f(Math.exp(x0)) + 4 * f(Math.exp(x0 + h / 2)) + f(Math.exp(x0 + h))); tot += s; as.push(Math.exp(x0 + h)); Ns.push(tot); }
  return { as, Ns, N: tot };
}
const CG_BASE = { material: 'Al 2024-T3', geom: 'centre crack', law: 'Paris', s_max: 100e6, R: 0, W: 100, r_hole: 0.003, a0: 0.001, a_det: 0.003, Kc_MPam: 0, gamma: 0.5, pn: 0.25, qn: 0.25, cyc_per_flight: 1, fh_per_flight: 1, scatter: 2, s_limit: 0, nSteps: 400 };
const crack = {
  id: 'crack', title: 'Crack growth, residual strength and inspection interval', fidelity: 'numerical',
  summary: 'Integrates a fatigue crack-growth law from an initial flaw to the critical size set by fracture toughness, gives the residual-strength curve, and sets the repeat inspection interval from the detectable-to-critical growth period.',
  equations: ['Paris–Erdogan crack growth law', 'Forman crack growth equation', 'Walker crack growth relation', 'Stress intensity factor equations', 'Griffith fracture criterion', 'Irwin fracture mechanics relations'],
  inputs: [
    MAT,
    { key: 'geom', label: 'Crack configuration', type: 'select', options: GEOMS, default: 'crack from hole (one side)', group: 'Geometry' },
    { key: 'W', label: 'Panel width', unit: 'm', default: 0.15, min: 0.005, group: 'Geometry', help: 'Stringer or crack-stopper pitch for a skin panel' },
    { key: 'r_hole', label: 'Hole radius', unit: 'm', default: 0.003, min: 0.0005, group: 'Geometry' },
    { key: 'a0', label: 'Initial flaw size', unit: 'm', default: 0.00127, min: 1e-5, group: 'Flaws', help: 'Crack length from the hole edge / half-length / depth. 1.27 mm (0.05 in) is the initial primary flaw assumed at holes for slow-crack-growth structure (JSSG-2006, USAF damage-tolerant design handbook §1.3.4.1)' },
    { key: 'a_det', label: 'Detectable crack size', unit: 'm', default: 0.005, min: 1e-5, group: 'Flaws', help: 'Depends on the inspection method and access: a few mm for eddy current, 25 mm or more for general visual' },
    { key: 's_max', label: 'Maximum gross stress in the cycle', unit: 'Pa', default: 100e6, min: 1e5, group: 'Loading', help: 'Equivalent once-per-flight stress' },
    { key: 'R', label: 'Stress ratio', unit: '-', default: 0, min: -1, max: 0.9, group: 'Loading', help: 'Negative R is treated as 0 (compressive part assumed not to open the crack)' },
    { key: 's_limit', label: 'Limit-load stress for residual strength (0 = cycle maximum)', unit: 'Pa', default: 0, min: 0, group: 'Loading', help: 'The critical crack is the size at which the structure can just carry this stress' },
    { key: 'cyc_per_flight', label: 'Equivalent cycles per flight', unit: '', default: 1, min: 0.01, group: 'Loading' },
    { key: 'fh_per_flight', label: 'Flight hours per flight', unit: 'h', default: 1.5, min: 0.01, group: 'Loading' },
    { key: 'law', label: 'Growth law', type: 'select', options: ['Paris', 'Walker', 'Forman', 'NASGRO (simplified)'], default: 'Forman', group: 'Growth law' },
    { key: 'Kc_MPam', label: 'Fracture toughness (0 = plane-strain KIc from database)', unit: 'MPa√m', default: 0, min: 0, group: 'Growth law', help: 'Thin sheet has a much higher plane-stress Kc than KIc; using KIc is conservative' },
    { key: 'gamma', label: 'Walker exponent γ', unit: '-', default: 0.5, min: 0.1, max: 1, group: 'Growth law' },
    { key: 'pn', label: 'NASGRO threshold exponent p', unit: '-', default: 0.25, min: 0, max: 2, group: 'Growth law' },
    { key: 'qn', label: 'NASGRO instability exponent q', unit: '-', default: 0.25, min: 0, max: 2, group: 'Growth law' },
    { key: 'scatter', label: 'Scatter factor on the growth period', unit: '-', default: 2, min: 1, max: 5, group: 'Inspection', help: 'Inspection interval = detectable-to-critical period / factor. 2 corresponds to the two-lifetime slow-crack-growth requirement (JSSG-2006) and to the factor of 2 of AC 25.571-1D where inspections are effective; 3 where they are not' },
    { key: 'nSteps', label: 'Integration intervals', unit: '', default: 400, min: 10, max: 20000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c, up) { const u = usage(c, up); return { material: c.struct.material, s_max: 1.3 * u.s1g, s_limit: u.n * u.s1g, fh_per_flight: u.fh }; },
  run(i) {
    const m = mat(i.material), Kc = i.Kc_MPam > 0 ? i.Kc_MPam : m.KIc / 1e6, p = { C: m.parisC, m: m.parisM, dKth: m.dKth, Kc, gamma: i.gamma, pn: i.pn, qn: i.qn }, warnings = [];
    const sRes = Math.max(i.s_limit, i.s_max), Kof = (a, s) => beta(i.geom, a, i.W, i.r_hole) * (s / 1e6) * Math.sqrt(Math.PI * a), aLim = aLimit(i.geom, i.W, i.r_hole), n = Math.max(10, Math.round(i.nSteps));
    const aLo = Math.min(i.a0, aLim) * 1e-3; let ac = aLim, netGov = true;
    if (Kof(aLim, sRes) > Kc) { ac = Kof(aLo, sRes) >= Kc ? aLo : N.brent((a) => Kof(a, sRes) - Kc, aLo, aLim, 1e-12); netGov = false; }
    if (netGov) warnings.push('Fracture toughness is not reached within the panel: the critical size was set at 95% of the ligament, where net-section yield governs.');
    const a0 = Math.min(i.a0, 0.999 * ac), ad = N.clamp(i.a_det, a0, 0.999 * ac), R = Math.max(i.R, 0), dK0 = Kof(a0, i.s_max) * (1 - R);
    if (i.a0 >= ac) warnings.push('The initial flaw is already critical at the residual-strength stress.');
    if (i.a_det >= ac) warnings.push('The detectable crack is not smaller than the critical crack: this detail is not inspectable and needs a safe-life or slow-crack-growth-without-inspection justification.');
    const g = growCrack(i, p, a0, ac, n), gd = growCrack(i, p, ad, ac, n), noGrow = !Number.isFinite(g.N), life = g.N / i.cyc_per_flight, per = gd.N / i.cyc_per_flight, interval = per / i.scatter;
    if (noGrow) warnings.push('ΔK at the initial flaw is below the threshold: no growth is predicted. Thresholds are sensitive to environment and load history, so confirm before relying on this.');
    if (i.R < 0) warnings.push('Negative stress ratio treated as R = 0.');
    const as = N.logspace(a0, ac, 80), resid = as.map((a) => Math.min((Kc / (beta(i.geom, a, i.W, i.r_hole) * Math.sqrt(Math.PI * a))), m.Sy / 1e6));
    const plots = [
      { type: 'line', title: 'Crack length versus flights', xlabel: 'Flights [-]', ylabel: 'Crack size [mm]', series: [{ name: 'Crack growth', x: thin(g.Ns.map((v) => (Number.isFinite(v) ? v / i.cyc_per_flight : 0))), y: thin(g.as.map((a) => a * 1e3)) }], annotations: [{ y: ad * 1e3, label: 'Detectable' }, { y: ac * 1e3, label: 'Critical' }] },
      { type: 'line', title: 'Residual strength', xlabel: 'Crack size [mm]', ylabel: 'Gross failure stress [MPa]', series: [{ name: 'Fracture (K = Kc), capped at yield', x: as.map((a) => a * 1e3), y: resid }, ...(i.geom === 'centre crack' ? [feddersen(Kc, m.Sy / 1e6, i.W)] : [])], annotations: [{ y: sRes / 1e6, label: 'Required' }] },
      { type: 'line', title: 'Growth rate curve', xlabel: 'ΔK [MPa√m]', ylabel: 'da/dN [m/cycle]', xlog: true, ylog: true, series: [{ name: i.law, x: N.logspace(Math.max(1, 0.8 * p.dKth), 0.98 * Kc * (1 - R), 60), y: N.logspace(Math.max(1, 0.8 * p.dKth), 0.98 * Kc * (1 - R), 60).map((d) => Math.max(1e-14, Math.min(1e-2, growthRate(i.law, d, R, p)))) }] },
    ];
    return {
      kpis: [
        kpi('crit_crack_m', 'Critical crack size', ac, 'm', undefined, `At ${(sRes / 1e6).toFixed(0)} MPa with Kc = ${Kc.toFixed(0)} MPa√m`),
        kpi('growth_cycles', 'Flights from initial flaw to critical', life, 'flights'),
        kpi('growth_fh', 'Crack-growth life in flight hours', life * i.fh_per_flight, 'h'),
        kpi('insp_period_flights', 'Flights from detectable to critical', per, 'flights'),
        kpi('inspection_interval_flights', 'Repeat inspection interval', interval, 'flights', interval > 500 ? 'ok' : 'warn', `Detectable-to-critical period / ${i.scatter}`),
        kpi('inspection_interval_fh', 'Repeat inspection interval in flight hours', interval * i.fh_per_flight, 'h'),
        kpi('threshold_insp_flights', 'Inspection threshold (first inspection)', life / i.scatter, 'flights', undefined, 'Initial-flaw-to-critical life / scatter factor'),
        kpi('dK_initial', 'ΔK at the initial flaw', dK0, 'MPa√m', dK0 > p.dKth ? undefined : 'ok', `Threshold ${p.dKth} MPa√m`),
        kpi('K_max_initial', 'Kmax at the initial flaw (cycle)', Kof(a0, i.s_max), 'MPa√m'),
        kpi('beta_initial', 'Geometry factor at the initial flaw', beta(i.geom, a0, i.W, i.r_hole), '-'),
        kpi('resid_strength_det_Pa', 'Residual strength at the detectable size', Math.min((Kc / (beta(i.geom, ad, i.W, i.r_hole) * Math.sqrt(Math.PI * ad))) * 1e6, m.Sy), 'Pa'),
      ],
      plots, warnings,
      models: [`${i.law} growth law` + (i.law === 'Forman' ? ' (coefficient normalised to the Paris constant at R = 0)' : i.law === 'NASGRO (simplified)' ? ' without the crack-closure function' : ''), `${i.geom}: ` + (i.geom === 'edge crack' ? 'Tada polynomial' : i.geom.includes('hole') ? 'Bowie solution curve-fit with secant finite-width correction' : 'secant (Feddersen) finite-width correction'), 'Simpson integration on a geometric crack-length grid'],
      assumptions: ['Through-thickness crack under constant-amplitude equivalent loading; no retardation, no load interaction', 'Linear elastic fracture mechanics; failure at Kmax = Kc under the residual-strength stress', 'Paris constants are quoted at R ≈ 0 in laboratory air; environment and thickness effects need test data', 'Walker and NASGRO exponents are generic placeholders, not sourced values (NASGRO/AFGROW material files should replace them); the default toughness is the plane-strain KIc', DATA_NOTE],
    };
  },
  convergence: { param: 'nSteps', label: 'Integration intervals', levels: [25, 50, 100, 200, 400], metric: 'growth_cycles' },
  verify() {
    const m = mat('Al 2024-T3'), s = 100, a0 = 0.001, r = N.kv(crack.run({ ...CG_BASE })), ac = (m.KIc / 1e6 / s) ** 2 / Math.PI, e = 1 - m.parisM / 2;
    const exact = (ac ** e - a0 ** e) / (m.parisC * (s * Math.sqrt(Math.PI)) ** m.parisM * e);
    return [
      N.check('Critical half-crack a = (Kc/σ)²/π in a wide panel', r.crit_crack_m, ac, 1e-4, 'Irwin: K = σ√(πa)'),
      N.check('Paris life against the closed-form integral', r.growth_cycles, exact, 1e-3, 'Closed-form integration of da/dN = C(Δσ√(πa))^m (finite-width effect < 0.1% at W = 100 m)'),
      N.check('Edge crack β → 1.122 for a short crack', beta('edge crack', 1e-9, 1, 1), 1.122, 1e-6, 'Tada, Paris & Irwin, Stress Analysis of Cracks Handbook'),
      N.check('One-sided hole crack cannot outgrow the ligament W/2 − r', N.kv(crack.run({ ...CG_BASE, geom: 'crack from hole (one side)', W: 0.1, s_max: 10e6 })).crit_crack_m, 0.95 * 0.9 * (0.05 - 0.003), 1e-12, 'Geometry: central hole of radius r in a panel of width W'),
      N.check('Both hole-crack configurations share the same ligament limit', aLimit('crack from hole (one side)', 0.1, 0.003), aLimit('cracks from hole (both sides)', 0.1, 0.003), 1e-12, 'Geometry'),
      N.check('Short crack at a hole: β → 1.12·Kt = 3.36', beta('crack from hole (one side)', 1e-9, 1e3, 0.003), 3.36, 3e-3, 'Bowie (1956); free-edge factor 1.12 × Kirsch Kt = 3'),
      N.check('Forman reduces to Paris when ΔK ≪ Kc at R = 0', growthRate('Forman', 0.01, 0, { C: 1e-11, m: 3, Kc: 1e6, dKth: 0 }), 1e-11 * 0.01 ** 3, 1e-6, 'Forman, Kearney & Engle (1967), normalised form'),
      N.check('Walker with γ = 1 is independent of R', growthRate('Walker', 10, 0.5, { C: 1e-11, m: 3, gamma: 1, dKth: 0 }), 1e-8, 1e-12, 'Walker (1970)'),
    ];
  },
  calibration: { params: [{ key: 'gamma', min: 0.1, max: 1 }, { key: 'Kc_MPam', min: 10, max: 250 }], sweep: 's_max', target: 'growth_cycles', note: 'Measured crack-growth lives of panels at several stress levels or ratios calibrate the Walker exponent and the apparent toughness; for C and m, fit da/dN–ΔK coupon data externally and edit the material constants.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    // the database toughness is the plane-strain KIc, a lower bound for sheet and thin plate: a verdict that rests on it is provisional
    const lowerBound = !(i.Kc_MPam > 0);
    if (i.a_det >= o.crit_crack_m) out.push({ severity: lowerBound ? 'warn' : 'critical', title: lowerBound ? 'Detail is not inspectable with the plane-strain toughness' : 'Detail is not inspectable', detail: `Critical crack ${(o.crit_crack_m * 1e3).toFixed(1)} mm is not larger than the detectable ${(i.a_det * 1e3).toFixed(1)} mm${lowerBound ? ' when the database KIc is used; thin gauges have a higher apparent toughness Kc' : ''}.`, action: (lowerBound ? 'Enter the measured Kc for the actual thickness first. If the result stands: l' : 'L') + 'ower the stress, use a tougher alloy or thinner (plane-stress) gauge, add a crack stopper, or apply a more sensitive inspection method.', basis: 'Damage tolerance: detectable before critical' });
    else if (o.inspection_interval_flights < 500) out.push({ severity: 'warn', title: 'Very short inspection interval', detail: `${o.inspection_interval_flights.toFixed(0)} flights (${o.inspection_interval_fh.toFixed(0)} h).`, action: 'Frequent inspections cost downtime and access damage: reduce the stress level, improve the detectable size (eddy current instead of visual) or redesign the detail.', basis: 'Detectable-to-critical period / scatter factor' });
    else out.push({ severity: 'info', title: 'Inspection programme', detail: `First inspection by ${o.threshold_insp_flights.toFixed(0)} flights, then every ${o.inspection_interval_flights.toFixed(0)} flights (${o.inspection_interval_fh.toFixed(0)} h).`, action: 'Align the interval with a scheduled check; a damage-tolerant detail can stay in service on condition instead of being retired at a fixed life, saving material and cost.', basis: 'Slow-crack-growth damage tolerance' });
    return out;
  },
};
const thin = (a, max = 200) => { if (a.length <= max) return a.slice(); const o = []; for (let k = 0; k < max; k++) o.push(a[Math.round((k * (a.length - 1)) / (max - 1))]); return o; };
/** Feddersen residual-strength construction for a centre-cracked panel: tangents from yield at zero crack and from zero stress at 2a = W. */
function feddersen(Kc, Sy, W) {
  const a1 = Math.min((9 / (4 * Math.PI)) * (Kc / Sy) ** 2, W / 6), a2 = W / 6, sc = (a) => Kc / Math.sqrt(Math.PI * a), x = [0, a1], y = [Sy, a1 < W / 6 ? (2 * Sy) / 3 : sc(a2)];
  for (const a of N.linspace(a1, a2, 12).slice(1)) { x.push(a); y.push(sc(a)); }
  x.push(W / 2); y.push(0);
  return { name: 'Feddersen construction', x: x.map((v) => v * 1e3), y, style: 'dash' };
}

// ---- 5. probabilistic damage tolerance ------------------------------------------------------
const pod = (a, a50, sg) => (a > 0 ? N.normCdf(Math.log(a / a50) / sg) : 0);
const probdt = {
  id: 'probdt', title: 'Probabilistic damage tolerance (Monte Carlo)', fidelity: 'numerical',
  summary: 'Samples the initial flaw size, crack-growth coefficient and toughness to give the distribution of crack-growth life, the probability of failure within the service life, and how a repeat inspection with a given probability of detection reduces it.',
  equations: ['Paris–Erdogan crack growth law', 'Walker crack growth relation', 'Stress intensity factor equations', 'Irwin fracture mechanics relations'],
  inputs: [
    MAT,
    { key: 'geom', label: 'Crack configuration', type: 'select', options: GEOMS, default: 'crack from hole (one side)', group: 'Geometry' },
    { key: 'W', label: 'Panel width', unit: 'm', default: 0.15, min: 0.005, group: 'Geometry' },
    { key: 'r_hole', label: 'Hole radius', unit: 'm', default: 0.003, min: 0.0005, group: 'Geometry' },
    { key: 's_max', label: 'Maximum gross stress per flight', unit: 'Pa', default: 100e6, min: 1e5, group: 'Loading' },
    { key: 'R', label: 'Stress ratio', unit: '-', default: 0, min: 0, max: 0.9, group: 'Loading' },
    { key: 's_limit', label: 'Residual-strength stress (0 = cycle maximum)', unit: 'Pa', default: 0, min: 0, group: 'Loading' },
    { key: 'a0_med', label: 'Median initial flaw size', unit: 'm', default: 0.0003, min: 1e-6, group: 'Distributions', help: 'Equivalent initial flaw size of manufacturing quality; fit to teardown or fractography data' },
    { key: 'a0_sd', label: 'Log-standard deviation of flaw size', unit: '-', default: 0.6, min: 0, max: 2, group: 'Distributions' },
    { key: 'C_sd', label: 'Log-standard deviation of growth coefficient C', unit: '-', default: 0.25, min: 0, max: 1.5, group: 'Distributions', help: 'Material scatter in da/dN; about 0.2–0.3 for aluminium plate' },
    { key: 'Kc_cov', label: 'Coefficient of variation of toughness', unit: '-', default: 0.08, min: 0, max: 0.4, group: 'Distributions' },
    { key: 'gamma', label: 'Walker exponent γ', unit: '-', default: 0.5, min: 0.1, max: 1, group: 'Loading' },
    { key: 'Kc_MPam', label: 'Mean fracture toughness (0 = plane-strain KIc from database)', unit: 'MPa√m', default: 0, min: 0, group: 'Distributions', help: 'Use the same value as in the crack-growth analysis; KIc is conservative for thin gauges' },
    { key: 'life_flights', label: 'Service life', unit: 'flights', default: 60000, min: 1, group: 'Service' },
    { key: 'interval', label: 'Repeat inspection interval', unit: 'flights', default: 6000, min: 1, group: 'Service', help: 'Defaults to the deterministic interval from the crack-growth analysis' },
    { key: 'a50', label: 'Crack size with 50% probability of detection', unit: 'm', default: 0.003, min: 1e-5, group: 'Service' },
    { key: 'pod_sd', label: 'Log-standard deviation of the POD curve', unit: '-', default: 0.5, min: 0.05, max: 2, group: 'Service', help: 'Log-normal POD(a) = Φ(ln(a/a50)/σ); from inspection reliability trials' },
    { key: 'seed', label: 'Random seed', unit: '', default: 7, min: 1, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nSamples', label: 'Monte Carlo samples', unit: '', default: 4000, min: 100, max: 200000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c, up) { const u = usage(c, up); return { material: c.struct.material, s_max: 1.3 * u.s1g, s_limit: u.n * u.s1g, life_flights: u.life_fl || undefined, interval: up.fatigue?.inspection_interval_flights > 1 ? up.fatigue.inspection_interval_flights : undefined }; },
  run(i) {
    const m = mat(i.material), aLim = aLimit(i.geom, i.W, i.r_hole), amin = Math.min(i.a0_med * Math.exp(-5 * i.a0_sd), aLim * 1e-3), ng = 600, as = N.logspace(amin, aLim, ng), sRes = Math.max(i.s_limit, i.s_max), warnings = [];
    // G(a) = ∫ da / ΔK_eff^m is independent of C, so each sample costs two table look-ups
    const dKe = as.map((a) => beta(i.geom, a, i.W, i.r_hole) * (i.s_max / 1e6) * (1 - i.R) * Math.sqrt(Math.PI * a) * (1 - i.R) ** (i.gamma - 1)), inv = dKe.map((d) => d ** -m.parisM), G = N.cumtrapz(as, inv);
    const Kres = as.map((a) => beta(i.geom, a, i.W, i.r_hole) * (sRes / 1e6) * Math.sqrt(Math.PI * a)); for (let k = 1; k < ng; k++) Kres[k] = Math.max(Kres[k], Kres[k - 1] * (1 + 1e-12));
    const u = N.rng(Math.round(i.seed)), n = Math.max(100, Math.round(i.nSamples)), lives = new Array(n), nI = Math.floor(i.life_flights / i.interval), Kc0 = i.Kc_MPam > 0 ? i.Kc_MPam : m.KIc / 1e6;
    const ints = [0.25, 0.5, 1, 2, 4].map((f) => f * i.interval), pfI = ints.map(() => 0); let pf0 = 0, pfIn = 0;
    const missTo = (g0, C, Nf, itvIn) => { const itv = Math.max(itvIn, i.life_flights / 400); let miss = 1; for (let k = 1; k * itv < Math.min(Nf, i.life_flights) && miss > 1e-12; k++) miss *= 1 - pod(N.interp1(G, as, g0 + C * k * itv), i.a50, i.pod_sd); return miss; };
    for (let s = 0; s < n; s++) {
      const a0 = N.clamp(i.a0_med * Math.exp(i.a0_sd * N.randn(u)), amin, aLim), C = m.parisC * Math.exp(i.C_sd * N.randn(u)), Kc = Math.max(0.3 * Kc0, Kc0 * (1 + i.Kc_cov * N.randn(u)));
      const ac = N.interp1(Kres, as, Kc), g0 = N.interp1(as, G, a0), Nf = Math.max(0, (N.interp1(as, G, ac) - g0) / C); lives[s] = Nf;
      if (Nf < i.life_flights) { pf0++; pfIn += missTo(g0, C, Nf, i.interval); ints.forEach((itv, j) => { pfI[j] += missTo(g0, C, Nf, itv); }); }
    }
    pf0 /= n; pfIn /= n;
    const lg = lives.map((v) => Math.log10(Math.max(v, 1))), hist = N.histogram(lg, 30), med = N.quantile(lives, 0.5), p01 = N.quantile(lives, 0.01), ap = N.logspace(i.a50 / 20, i.a50 * 20, 60);
    if (pf0 > 0 && pf0 * n < 20) warnings.push('Fewer than 20 failures were sampled: the failure probability has a large sampling error. Increase the number of samples.');
    if (pf0 === 0) warnings.push(`No failure in ${n} samples: the probability of failure is below about ${(3 / n).toExponential(1)} (95% confidence) and cannot be resolved by plain Monte Carlo.`);
    if (nI > 400) warnings.push('More than 400 inspections in the service life: the simulation inspects at most 400 times (every ' + (i.life_flights / 400).toFixed(0) + ' flights), which is conservative.');
    if (nI < 1) warnings.push('The inspection interval is longer than the service life: no inspection takes place.');
    return {
      kpis: [
        kpi('pof_no_inspection', 'Probability of failure in service life, no inspection', pf0, '-', pf0 < 1e-3 ? 'ok' : pf0 < 0.05 ? 'warn' : 'bad'),
        kpi('pof_with_inspection', 'Probability of failure with repeat inspection', pfIn, '-', pfIn < 1e-3 ? 'ok' : pfIn < 0.01 ? 'warn' : 'bad', 'Per detail over the service life'),
        kpi('risk_reduction', 'Risk reduction factor from inspection', pfIn > 0 ? pf0 / pfIn : 1, '-', undefined, pf0 > 0 ? '' : 'No failures sampled: not resolvable'),
        kpi('life_median_flights', 'Median crack-growth life', med, 'flights'),
        kpi('life_p01_flights', '1st-percentile crack-growth life', p01, 'flights', p01 >= i.life_flights ? 'ok' : 'warn'),
        kpi('life_cov', 'Scatter of life (log10 standard deviation)', N.std(lg), '-'),
        kpi('n_inspections', 'Inspections in the service life', nI, '-'),
        kpi('pof_std_error', 'Sampling standard error of the no-inspection probability', Math.sqrt((pf0 * (1 - pf0)) / n), '-'),
      ],
      plots: [
        { type: 'bar', title: 'Distribution of crack-growth life', ylabel: 'Samples [-]', categories: hist.centers.map((c) => '10^' + c.toFixed(1)), series: [{ name: 'Samples', y: hist.counts }] },
        { type: 'line', title: 'Probability of detection', xlabel: 'Crack size [mm]', ylabel: 'POD [-]', xlog: true, series: [{ name: 'POD(a)', x: ap.map((a) => a * 1e3), y: ap.map((a) => pod(a, i.a50, i.pod_sd)) }] },
        { type: 'line', title: 'Probability of failure versus inspection interval', xlabel: 'Inspection interval [flights]', ylabel: 'Probability of failure in service life [-]', xlog: true, ylog: true, series: [{ name: 'With inspection', x: ints, y: pfI.map((v) => Math.max(v / n, 1e-12)), style: 'line+points' }], annotations: [{ y: Math.max(pf0, 1e-12), label: 'No inspection' }] },
      ],
      warnings,
      models: ['Monte Carlo with seeded generator', 'Walker-corrected Paris growth integrated once as G(a) and scaled per sample', 'Log-normal initial flaw and growth coefficient, normal toughness', 'Log-normal probability-of-detection curve; detected cracks are assumed repaired'],
      assumptions: ['One crack per detail, independent inspections, perfect repair on detection', 'Constant-amplitude equivalent loading; load scatter is not sampled', 'Distribution parameters and the probability-of-detection curve (a50, log-sd) are user inputs: defaults are illustrative, not fleet or inspection-trial data (MIL-HDBK-1823A gives the model form only)', DATA_NOTE],
    };
  },
  convergence: { param: 'nSamples', label: 'Monte Carlo samples', levels: [500, 1000, 2000, 4000, 8000], metric: 'life_median_flights', hOf: (n) => 1 / Math.sqrt(n) },
  verify() {
    const b = { material: 'Al 2024-T3', geom: 'centre crack', W: 100, r_hole: 0.003, s_max: 100e6, R: 0, s_limit: 0, a0_med: 0.001, a0_sd: 0, C_sd: 0, Kc_cov: 0, gamma: 0.5, Kc_MPam: 0, life_flights: 1e9, interval: 1e9, a50: 0.003, pod_sd: 0.5, seed: 1, nSamples: 100 };
    const r = N.kv(probdt.run(b)), det = N.kv(crack.run({ ...CG_BASE, nSteps: 800 }));
    return [
      N.check('Zero scatter reproduces the deterministic crack-growth life', r.life_median_flights, det.growth_cycles, 2e-3, 'Consistency with the crack-growth analysis'),
      N.check('Zero scatter: life has no spread', r.life_cov, 0, 1e-9, 'Limit case'),
      N.check('POD(a50) = 0.5', pod(0.003, 0.003, 0.5), 0.5, 1e-7, 'Definition of the log-normal POD curve'),
      N.check('No scatter and life beyond the crack-growth life: failure is certain', r.pof_no_inspection, 1, 0, 'Limit case'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.pof_with_inspection > 1e-3) out.push({ severity: 'warn', title: 'Residual risk is high even with inspection', detail: `Probability of failure ${o.pof_with_inspection.toExponential(1)} per detail over ${i.life_flights.toFixed(0)} flights.`, action: 'Shorten the interval, use a method with a smaller a50, or reduce stress; the interval plot shows the trade.', basis: 'Risk-based inspection planning' });
    else if (o.pof_no_inspection < 1e-4) out.push({ severity: 'advise', title: 'Inspection adds little here', detail: `Risk without inspection is already ${o.pof_no_inspection.toExponential(1)}.`, action: 'A longer interval or opportunity-based inspection may be justified, reducing maintenance burden and disassembly damage.', basis: 'Risk-based inspection planning' });
    out.push({ severity: 'info', title: 'Calibrate the flaw distribution', detail: `Median flaw ${(i.a0_med * 1e3).toFixed(2)} mm, log-sd ${i.a0_sd}.`, action: 'Teardown or fractography data from the fleet tightens the initial-flaw distribution, which usually dominates the risk estimate.', basis: 'Equivalent initial flaw size method' });
    return out;
  },
};

// ---- 6. elastic–plastic fracture measures and FAD -------------------------------------------
/** R6 Option 1 failure assessment curve. */
const fadCurve = (Lr) => (1 - 0.14 * Lr * Lr) * (0.3 + 0.7 * Math.exp(-0.65 * Lr ** 6));
const epfm = {
  id: 'epfm', title: 'Fracture parameters, plastic zone and failure assessment diagram', fidelity: 'analytical',
  summary: 'Stress intensity, energy release rate (J in small-scale yielding), crack-tip opening displacement by the Irwin and Dugdale models, plastic-zone size, plane-strain validity, and a failure assessment diagram that covers the transition from brittle fracture to plastic collapse.',
  equations: ['Stress intensity factor equations', 'Griffith fracture criterion', 'Irwin fracture mechanics relations', 'J-integral formulation', 'Crack-tip opening displacement relations'],
  inputs: [
    MAT,
    { key: 'geom', label: 'Crack configuration', type: 'select', options: ['centre crack', 'edge crack'], default: 'centre crack', group: 'Geometry' },
    { key: 'a', label: 'Crack size (half-length or depth)', unit: 'm', default: 0.01, min: 1e-5, group: 'Geometry', help: 'Defaults to half the critical half-length of a centre crack at the limit stress with the database KIc' },
    { key: 'W', label: 'Panel width', unit: 'm', default: 0.15, min: 0.005, group: 'Geometry' },
    { key: 'B', label: 'Thickness', unit: 'm', default: 0.003, min: 1e-4, group: 'Geometry' },
    { key: 'sigma', label: 'Gross stress', unit: 'Pa', default: 150e6, min: 0, group: 'Loading', help: 'Limit-load stress from Suite 2' },
    { key: 'Kmat_MPam', label: 'Material toughness (0 = KIc from database)', unit: 'MPa√m', default: 0, min: 0, group: 'Material' },
  ],
  defaults(c, up) {
    // crack assessed by default: half the critical half-length of a centre crack at the limit stress, kept inside the default panel
    const u = usage(c, up), m = mat(c.struct.material), sig = u.n * u.s1g, ac = sig > 0 ? (m.KIc / sig) ** 2 / Math.PI : NaN;
    return { material: c.struct.material, sigma: sig, B: up.fea?.t_skin_root_m ?? Math.min(c.struct.t_skin_mm, 3) / 1e3, a: Number.isFinite(ac) ? N.clamp(0.5 * ac, 2e-4, 0.03) : undefined };
  },
  run(i) {
    const m = mat(i.material), Kmat = i.Kmat_MPam > 0 ? i.Kmat_MPam * 1e6 : m.KIc, warnings = [], a = Math.min(i.a, (i.geom === 'edge crack' ? 0.65 : 0.47) * i.W);
    const Bmin = 2.5 * (Kmat / m.Sy) ** 2, pe = i.B >= Bmin, Ep = pe ? m.E / (1 - m.nu ** 2) : m.E, K = beta(i.geom, a, i.W, 1) * i.sigma * Math.sqrt(Math.PI * a);
    const rp = (1 / ((pe ? 6 : 2) * Math.PI)) * (K / m.Sy) ** 2, Keff = beta(i.geom, a + rp, i.W, 1) * i.sigma * Math.sqrt(Math.PI * (a + rp)), J = (K * K) / Ep, ctodI = (K * K) / (m.Sy * Ep) * (pe ? 0.5 : 1);
    const ratio = Math.min(i.sigma / m.Sy, 0.999), ctodD = ((8 * m.Sy * a) / (Math.PI * m.E)) * Math.log(1 / Math.cos((Math.PI * ratio) / 2)), Jd = m.Sy * ctodD;
    const net = i.geom === 'edge crack' ? 1 - a / i.W : 1 - (2 * a) / i.W, Lr = i.sigma / (net * m.Sy), Kr = K / Kmat, LrMax = (m.Sy + m.Su) / (2 * m.Sy);
    // reserve factor: scale the load until the assessment point reaches the curve or the collapse cut-off
    let rf = Infinity; if (Kr > 0 || Lr > 0) { const f = (s) => (s * Lr > LrMax ? -1 : fadCurve(s * Lr) - s * Kr); let lo = 0, hi = 1; while (f(hi) > 0 && hi < 1e6) { lo = hi; hi *= 2; } for (let k = 0; k < 60; k++) { const mid = 0.5 * (lo + hi); if (f(mid) > 0) lo = mid; else hi = mid; } rf = lo; }
    const sGriff = Kmat / (beta(i.geom, a, i.W, 1) * Math.sqrt(Math.PI * a));
    if (a < i.a) warnings.push('The crack was clipped to the validity range of the geometry factor.');
    if (!pe) warnings.push(`Thickness ${(i.B * 1e3).toFixed(1)} mm is below the plane-strain requirement ${(Bmin * 1e3).toFixed(1)} mm: the real toughness is higher than KIc (conservative), and tearing resistance (R-curve) behaviour should be used for an accurate residual strength.`);
    if (rp > 0.1 * a) warnings.push('The plastic zone exceeds 10% of the crack length: small-scale-yielding LEFM is inaccurate; rely on the failure assessment diagram.');
    const Ls = N.linspace(0, LrMax, 60);
    return {
      kpis: [
        kpi('K_MPam', 'Stress intensity factor K', K / 1e6, 'MPa√m', K < Kmat ? 'ok' : 'bad', `Toughness ${(Kmat / 1e6).toFixed(0)} MPa√m`),
        kpi('K_eff_MPam', 'K with Irwin plastic-zone correction', Keff / 1e6, 'MPa√m'),
        kpi('J_Jm2', 'Energy release rate G = J (small-scale yielding)', J, 'J/m²'), kpi('J_dugdale_Jm2', 'J from the Dugdale strip-yield model', Jd, 'J/m²', undefined, 'Plane stress, infinite plate'),
        kpi('Jc_Jm2', 'Critical energy release rate Kc²/E′', (Kmat * Kmat) / Ep, 'J/m²'),
        kpi('ctod_m', 'Crack-tip opening displacement K²/(m·σy·E′)', ctodI, 'm', undefined, 'Small-scale-yielding estimate with m = 1 in plane stress and 2 in plane strain'), kpi('ctod_dugdale_m', 'Crack-tip opening displacement (Dugdale)', ctodD, 'm'),
        kpi('plastic_zone_m', 'Plastic-zone radius', rp, 'm', rp < 0.1 * a ? 'ok' : 'warn', pe ? 'Plane strain' : 'Plane stress'),
        kpi('B_plane_strain_m', 'Thickness needed for plane strain', Bmin, 'm'),
        kpi('Kr', 'Fracture ratio Kr', Kr, '-'), kpi('Lr', 'Load ratio Lr', Lr, '-'),
        kpi('fad_reserve', 'Reserve factor on load (failure assessment diagram)', rf, '-', rf >= 1.1 ? 'ok' : rf >= 1 ? 'warn' : 'bad', 'Load multiplier to reach the assessment curve; a cracked structure must still carry limit load, i.e. at least 1 when the stress entered is the limit-load stress'),
        kpi('sigma_fracture_Pa', 'LEFM fracture stress for this crack', sGriff, 'Pa'),
      ],
      outputs: { plane_strain: pe ? 1 : 0 },
      plots: [{ type: 'line', title: 'Failure assessment diagram', xlabel: 'Load ratio Lr = σ_net/σ_y [-]', ylabel: 'Fracture ratio Kr = K/K_mat [-]', series: [{ name: 'Assessment curve (R6 Option 1)', x: [...Ls, LrMax], y: [...Ls.map(fadCurve), 0] }, { name: 'Loading line', x: [0, Lr * Math.min(rf, 50)], y: [0, Kr * Math.min(rf, 50)], style: 'dash' }, { name: 'Assessment point', x: [Lr], y: [Kr], style: 'points' }] }],
      warnings,
      models: ['Linear elastic fracture mechanics with secant / Tada geometry factors', 'Irwin plastic-zone correction', 'Dugdale strip-yield model', 'R6 Option 1 failure assessment curve with collapse cut-off at the flow stress'],
      assumptions: ['Through crack, mode I, monotonic load', 'J equals G only in small-scale yielding; no fully plastic (EPRI h-function) J solution is evaluated', 'Reference stress taken as the net-section stress', DATA_NOTE],
    };
  },
  verify() {
    const m = mat('Al 7075-T6'), b = { material: 'Al 7075-T6', geom: 'centre crack', a: 0.005, W: 100, B: 0.001, sigma: 10e6, Kmat_MPam: 0 }, r = N.kv(epfm.run(b)), br = N.kv(epfm.run({ ...b, Kmat_MPam: 4 }));
    return [
      N.check('K = σ√(πa) for a small crack in a wide plate', r.K_MPam, 10 * Math.sqrt(Math.PI * 0.005), 1e-6, 'Irwin (1957)'),
      N.check('Dugdale CTOD → K²/(σy·E) at low stress', r.ctod_dugdale_m, (r.K_MPam * 1e6) ** 2 / (m.Sy * m.E), 2e-4, 'Dugdale (1960) / Burdekin & Stone (1966), small-load expansion'),
      N.check('Assessment curve f(0) = 1', fadCurve(0), 1, 1e-12, 'R6 Option 1'),
      N.check('Brittle limit: reserve factor = K_mat/K', br.fad_reserve, 4 / br.K_MPam, 2e-3, 'Lr → 0 limit of the failure assessment diagram'),
      N.check('G = K²/E in plane stress', r.J_Jm2, (r.K_MPam * 1e6) ** 2 / m.E, 1e-10, 'Irwin relation'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    // with the database KIc in a gauge below the plane-strain thickness the fracture ratio is a lower bound, so a fracture-dominated failure is provisional
    const lowerBound = !(i.Kmat_MPam > 0) && !o.plane_strain && o.Kr > o.Lr;
    if (o.fad_reserve < 1) out.push({ severity: lowerBound ? 'warn' : 'critical', title: lowerBound ? 'Residual strength is short with the plane-strain toughness' : 'The cracked section cannot carry the applied stress', detail: `Reserve factor ${o.fad_reserve.toFixed(2)} (Kr = ${o.Kr.toFixed(2)}, Lr = ${o.Lr.toFixed(2)})${lowerBound ? '; the gauge is below the plane-strain thickness, where the apparent toughness exceeds KIc' : ''}.`, action: lowerBound ? 'Enter the measured Kc (or R-curve result) for this thickness; if the reserve stays below 1, use a tougher material or reduce the crack size limit through inspection.' : o.Kr > o.Lr ? 'Fracture-dominated: use a tougher material or reduce the crack size limit through inspection.' : 'Collapse-dominated: add section or a redundant load path.', basis: 'Failure assessment diagram; residual strength ≥ limit load' });
    else if (o.fad_reserve < 1.1) out.push({ severity: 'advise', title: 'Thin residual-strength margin', detail: `Reserve factor ${o.fad_reserve.toFixed(2)} on the stress entered.`, action: 'Tighten the inspection interval so that cracks are found well before this size.', basis: 'Failure assessment diagram; residual strength ≥ limit load' });
    return out;
  },
};

export default {
  id: 'fatigue', n: 9,
  tagline: 'How long before a crack starts, how fast does it grow, and how often must we look?',
  analyses: [sn, en, spectrum, crack, probdt, epfm],
  consumes: [
    { from: 'fea', keys: ['sigma_bend_root_Pa', 'sigma_max_Pa', 't_skin_root_m'], why: 'Stress level at the fatigue-critical location' },
    { from: 'performance', keys: ['n_limit'], why: 'Limit load factor to convert limit stress to the 1 g stress and to cap the spectrum' },
  ],
  provides: [
    { key: 'life_cycles', label: 'Fatigue life', unit: 'flights' }, { key: 'life_fh', label: 'Fatigue life', unit: 'h' }, { key: 'crit_crack_m', label: 'Critical crack size', unit: 'm' },
    { key: 'inspection_interval_fh', label: 'Inspection interval', unit: 'h' }, { key: 'damage_per_flight', label: 'Damage per flight', unit: '1/flight' },
  ],
  handoff: [
    { model: 'Coupled finite-element fracture mechanics, XFEM and crack-growth remeshing', why: 'Needs 3-D meshes with enriched or adaptively remeshed crack fronts; handbook geometry factors are used instead', tool: 'FE-based fracture codes (FRANC3D / Zencrack / Abaqus XFEM class), NASGRO / AFGROW stress-intensity libraries' },
    { model: 'Cohesive-zone and continuum-damage-mechanics crack simulation', why: 'Requires calibrated traction–separation or damage-evolution laws and nonlinear FE', tool: 'Nonlinear FE with cohesive elements' },
    { model: 'Load-interaction (retardation) models and full NASGRO equation with crack closure', why: 'Willenborg/Wheeler/strip-yield closure models need cycle-by-cycle integration with calibrated parameters; the simplified NASGRO form here omits closure', tool: 'NASGRO / AFGROW' },
    { model: 'Fully plastic J-integral (EPRI h-functions), J–R tearing analysis', why: 'Tabulated h-functions and material R-curves are geometry- and hardening-specific; only small-scale-yielding J, Dugdale and the R6 curve are evaluated', tool: 'EPFM handbooks / elastic–plastic FE' },
    { model: 'Thermomechanical, multiaxial and fretting fatigue', why: 'Critical-plane and TMF models need multiaxial stress histories from FE and dedicated material data', tool: 'FE-based fatigue post-processors (fe-safe / nCode class)' },
    { model: 'Composite progressive fatigue damage', why: 'No general validated closed-form model; laminate static failure is in Suite 21', tool: 'Test-based life factors; progressive damage FE' },
    { model: 'Corrosion-fatigue interaction kinetics', why: 'Environment-, frequency- and alloy-specific; represented only by an empirical knock-down factor and the pit-as-initial-flaw route from Suite 21', tool: 'Environmental fatigue testing' },
  ],
};
