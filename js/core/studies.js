// Generic study engine. Every analysis in every suite gets, for free and in a uniform way:
// input resolution from the shared case and upstream suites, mesh/resolution convergence with GCI,
// sensitivity ranking, uncertainty propagation, parameter calibration and validation against data.

import * as N from './numerics.js';
import { isa } from './atmosphere.js';
import { derived } from './case.js';

export function makeCtx(caseData, up = {}, progress = () => {}) {
  return { case: caseData, d: derived(caseData), up, atm: isa(caseData.atm.alt_m, caseData.atm.dISA_K), progress };
}

/** Whether an analysis applies to this vehicle; returns true or a human-readable reason. */
export function applicable(an, ctx) {
  if (!an.applicable) return true;
  try { return an.applicable(ctx.case, ctx.d); } catch { return true; }
}

/** Build the input set: declared defaults < case/upstream-derived defaults < user overrides. */
export function resolveInputs(an, ctx, overrides = {}) {
  const inp = {};
  for (const f of an.inputs || []) inp[f.key] = f.default;
  let src = {};
  if (an.defaults) { try { src = an.defaults(ctx.case, ctx.up, ctx.d) || {}; } catch { src = {}; } }
  for (const k of Object.keys(src)) if (src[k] !== undefined && src[k] !== null && !(typeof src[k] === 'number' && !Number.isFinite(src[k]))) inp[k] = src[k];
  const linked = Object.keys(src);
  for (const k of Object.keys(overrides)) if (overrides[k] !== undefined && overrides[k] !== '') inp[k] = overrides[k];
  // coerce numeric fields
  for (const f of an.inputs || []) if ((f.type || 'number') === 'number') { const v = Number(inp[f.key]); inp[f.key] = Number.isFinite(v) ? v : f.default; }
  return { inp, linked };
}

/** Run one analysis and normalise its result. */
export async function runAnalysis(an, inp, ctx) {
  const t0 = Date.now();
  const res = (await an.run(inp, ctx)) || {};
  res.kpis = res.kpis || []; res.plots = res.plots || []; res.tables = res.tables || [];
  res.warnings = res.warnings || []; res.models = res.models || []; res.assumptions = res.assumptions || [];
  res.outputs = res.outputs || {};
  for (const k of res.kpis) if (k.key && !(k.key in res.outputs) && typeof k.value === 'number') res.outputs[k.key] = k.value;
  res.elapsed_ms = Date.now() - t0;
  return res;
}
const metricOf = (res, key) => { const v = res.outputs?.[key]; return typeof v === 'number' ? v : NaN; };

/**
 * Resolution/mesh convergence study with Richardson extrapolation and the Grid Convergence Index.
 * Also performs mesh selection: the coarsest level whose estimated discretisation error is within `target`.
 */
export async function convergenceStudy(an, inp, ctx, { levels, metric, target = 0.01 } = {}) {
  const spec = an.convergence;
  if (!spec) throw new Error('This analysis has no discretisation parameter.');
  levels = (levels || spec.levels).slice().sort((a, b) => a - b); metric = metric || spec.metric;
  const hOf = spec.hOf || ((n) => 1 / n), rows = [];
  for (let i = 0; i < levels.length; i++) {
    ctx.progress(i / levels.length, `Resolution ${levels[i]}`);
    const r = await runAnalysis(an, { ...inp, [spec.param]: levels[i] }, ctx);
    rows.push({ level: levels[i], h: hOf(levels[i]), value: metricOf(r, metric), ms: r.elapsed_ms });
  }
  const out = { param: spec.param, label: spec.label || spec.param, metric, rows, target };
  if (rows.length >= 3) {
    const k = rows.length - 1, g = N.gci([rows[k].h, rows[k - 1].h, rows[k - 2].h], [rows[k].value, rows[k - 1].value, rows[k - 2].value]);
    out.gci = g;
    const ref = Number.isFinite(g.fExact) && g.monotonic && g.p > 0.3 && g.p < 12 ? g.fExact : rows[k].value;
    out.reference = ref; out.referenceKind = ref === rows[k].value ? 'finest level' : 'Richardson extrapolation';
    for (const r of rows) r.err = Math.abs((r.value - ref) / (ref || 1e-300));
    const pick = rows.find((r) => r.err <= target) || rows[k];
    out.selected = pick.level; out.selectedReason = pick.err <= target
      ? `Coarsest resolution whose estimated discretisation error (${(100 * pick.err).toPrecision(2)}%) is within the ${(100 * target).toPrecision(2)}% target.`
      : `No tested level met the ${(100 * target).toPrecision(2)}% target; the finest was chosen. Add finer levels.`;
    out.verdict = !g.monotonic ? 'Oscillatory or non-monotonic convergence: treat the error estimate with caution.'
      : Math.abs(g.asymptotic - 1) < 0.1 ? 'Solutions are in the asymptotic range of convergence.' : 'Not yet clearly in the asymptotic range; refine further for a firmer estimate.';
  }
  ctx.progress(1, 'done');
  return out;
}

/** One-at-a-time sweep of a single input. Returns per-value metric table for all numeric outputs. */
export async function sweep(an, inp, ctx, key, values) {
  const rows = [];
  for (let i = 0; i < values.length; i++) {
    ctx.progress(i / values.length, `${key} = ${values[i]}`);
    try { const r = await runAnalysis(an, { ...inp, [key]: values[i] }, ctx); rows.push({ x: values[i], out: numericOutputs(r) }); } catch (e) { rows.push({ x: values[i], out: {}, error: e.message }); }
  }
  ctx.progress(1, 'done');
  return rows;
}
export const numericOutputs = (r) => Object.fromEntries(Object.entries(r.outputs || {}).filter(([, v]) => typeof v === 'number' && Number.isFinite(v)));

/** Local sensitivity: normalised derivative (d metric / metric) / (d input / input) for each numeric input. */
export async function sensitivity(an, inp, ctx, metric, keys, pct = 0.05) {
  const base = metricOf(await runAnalysis(an, inp, ctx), metric), rows = [];
  keys = keys || (an.inputs || []).filter((f) => (f.type || 'number') === 'number' && !f.discrete).map((f) => f.key);
  for (let i = 0; i < keys.length; i++) {
    const k = keys[i], x = inp[k];
    ctx.progress(i / keys.length, k);
    if (typeof x !== 'number' || x === 0) continue;
    try {
      const hi = metricOf(await runAnalysis(an, { ...inp, [k]: x * (1 + pct) }, ctx), metric), lo = metricOf(await runAnalysis(an, { ...inp, [k]: x * (1 - pct) }, ctx), metric);
      if (Number.isFinite(hi) && Number.isFinite(lo)) rows.push({ key: k, lo, hi, elasticity: base ? (hi - lo) / (2 * pct * base) : 0 });
    } catch { /* a perturbed point outside the model's validity is skipped */ }
  }
  rows.sort((a, b) => Math.abs(b.elasticity) - Math.abs(a.elasticity));
  ctx.progress(1, 'done');
  return { metric, base, pct, rows };
}

/**
 * Uncertainty propagation by Latin hypercube or Monte Carlo sampling.
 * vars: [{key, type:'normal'|'uniform'|'lognormal'|'triangular', cov (fraction) | sd | min | max}]
 */
export async function uqStudy(an, inp, ctx, { vars, n = 200, method = 'lhs', seed = 2024, metrics } = {}) {
  const u = N.rng(seed), d = vars.length, U = method === 'lhs' ? N.lhs(n, d, u) : N.range(n, () => N.range(d, () => u()));
  const dist = vars.map((v) => {
    const m = inp[v.key];
    if (v.type === 'uniform') return { type: 'uniform', min: v.min ?? m * (1 - (v.cov || 0.1) * Math.sqrt(3)), max: v.max ?? m * (1 + (v.cov || 0.1) * Math.sqrt(3)) };
    if (v.type === 'triangular') return { type: 'triangular', min: v.min ?? m * 0.8, mode: m, max: v.max ?? m * 1.2 };
    return { type: v.type || 'normal', mean: m, sd: v.sd ?? Math.abs(m) * (v.cov || 0.05) };
  });
  const X = [], Y = {}; let failed = 0;
  for (let i = 0; i < n; i++) {
    if (i % 10 === 0) ctx.progress(i / n, `Sample ${i}/${n}`);
    const x = U[i].map((p, j) => N.sampleDist(dist[j], N.clamp(p, 1e-9, 1 - 1e-9))), s = { ...inp };
    vars.forEach((v, j) => (s[v.key] = x[j]));
    try {
      const o = numericOutputs(await runAnalysis(an, s, ctx));
      X.push(x);
      for (const k of metrics || Object.keys(o)) (Y[k] = Y[k] || []).push(o[k]);
    } catch { failed++; }
  }
  const stats = {};
  for (const [k, arr] of Object.entries(Y)) {
    if (arr.length !== X.length || arr.some((v) => !Number.isFinite(v))) continue;
    const m = N.mean(arr), sd = N.std(arr);
    if (sd === 0) continue;
    stats[k] = {
      mean: m, sd, cov: m ? sd / Math.abs(m) : NaN, p05: N.quantile(arr, 0.05), p50: N.quantile(arr, 0.5), p95: N.quantile(arr, 0.95),
      ci95: 1.96 * sd / Math.sqrt(arr.length), hist: N.histogram(arr, 24),
      // squared rank-free correlation as a first-order importance indicator
      importance: vars.map((v, j) => ({ key: v.key, r: N.corr(X.map((x) => x[j]), arr) })),
    };
  }
  ctx.progress(1, 'done');
  return { n: X.length, failed, method, vars, dist, stats };
}

/**
 * Calibrate chosen inputs so that a swept model output matches measured data.
 * data rows: [sweepValue, observed]. Returns fitted parameters, standard errors and fit quality.
 */
export async function calibrate(an, inp, ctx, { params, sweepKey, target, data }) {
  const xs = data.map((r) => r[0]), obs = data.map((r) => r[1]);
  const p0 = params.map((p) => inp[p.key]), scale = p0.map((v) => Math.abs(v) || 1);
  const predict = async (p) => {
    const s = { ...inp }; params.forEach((q, j) => (s[q.key] = N.clamp(p[j] * scale[j], q.min ?? -Infinity, q.max ?? Infinity)));
    const out = [];
    for (const x of xs) out.push(metricOf(await runAnalysis(an, { ...s, [sweepKey]: x }, ctx), target));
    return out;
  };
  // synchronous optimiser over an async model: evaluate through a small cache-driven loop
  const cache = new Map(), key = (p) => p.map((v) => v.toPrecision(12)).join(',');
  const before = await predict(p0.map((v, j) => v / scale[j]));
  let pending = null;
  const resid = (p) => { const k = key(p); if (cache.has(k)) return cache.get(k); pending = p; throw pending; };
  let fit;
  for (let guard = 0; guard < 5000; guard++) {
    try { fit = N.levenbergMarquardt(resid, p0.map((v, j) => v / scale[j]), { maxIter: 40 }); break; }
    catch (e) { if (e !== pending) throw e; const pr = await predict(pending); cache.set(key(pending), pr.map((v, i) => (Number.isFinite(v) ? v - obs[i] : 1e6))); ctx.progress(Math.min(0.95, guard / 400), 'Fitting'); }
  }
  if (!fit) throw new Error('Calibration did not converge.');
  const after = await predict(fit.p);
  ctx.progress(1, 'done');
  return {
    params: params.map((q, j) => ({ key: q.key, initial: p0[j], fitted: N.clamp(fit.p[j] * scale[j], q.min ?? -Infinity, q.max ?? Infinity), se: fit.cov ? Math.sqrt(Math.abs(fit.cov[j][j])) * scale[j] : NaN })),
    x: xs, observed: obs, before, after, metricsBefore: N.errorMetrics(before, obs), metricsAfter: N.errorMetrics(after, obs), evaluations: cache.size,
  };
}

/** Compare predictions with a reference dataset. vset: {inputs, sweep:{key, values}, target, observed, tol_pct}. */
export async function validate(an, inp, ctx, vset) {
  const base = { ...inp, ...(vset.inputs || {}) }, pred = [];
  for (const x of vset.sweep.values) pred.push(metricOf(await runAnalysis(an, { ...base, [vset.sweep.key]: x }, ctx), vset.target));
  const m = N.errorMetrics(pred, vset.observed), tol = vset.tol_pct ?? 10;
  const within = pred.filter((p, i) => Math.abs(p - vset.observed[i]) <= (tol / 100) * Math.max(Math.abs(vset.observed[i]), 1e-12)).length;
  return { ...vset, predicted: pred, metrics: m, within, pass: within === pred.length };
}

/** Run every self-test declared by a suite's analyses. */
export function verifySuite(def) {
  const out = [];
  for (const an of def.analyses) if (an.verify) {
    let checks;
    try { checks = an.verify(); } catch (e) { checks = [{ name: 'verification threw', actual: NaN, expected: 0, tol: 0, error: NaN, pass: false, ref: e.message }]; }
    for (const c of checks) out.push({ analysis: an.id, ...c });
  }
  return out;
}
