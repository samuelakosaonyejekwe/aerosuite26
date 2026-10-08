// Suite workbench: pick an analysis, set it up, run it, study it, verify it, and read its specification.

import { h, setKids, icon, clear, add, num, btn, badge, card, empty, toast, debounce, ago } from '../dom.js';
import { renderPlot } from '../plots.js';
import { resultPanel, dataTable, recList } from '../results.js';
import { SUITES, suiteMeta } from '../../core/registry.js';
import { state, getOverrides, setOverride, clearOverrides, saveRun, loadRun, on, setSetting } from '../../core/store.js';
import { runJob, cancelAll } from '../../core/runner.js';
import { loadSpec, terms } from '../../core/spec.js';
import { live, SUITE_QUERIES } from '../../core/live.js';
import { replaceHash } from '../../app.js';

const TABS = [['run', 'Set up & run', 'play'], ['mesh', 'Mesh & convergence', 'layers'], ['studies', 'Studies', 'sliders'], ['vv', 'Verification & validation', 'check'], ['spec', 'Specification', 'doc'], ['live', 'Live resources', 'globe']];
const FID = { analytical: 'Closed-form', 'reduced-order': 'Reduced-order model', numerical: 'Numerical solver' };
const parseList = (s) => String(s).split(/[\s,;]+/).map(Number).filter(Number.isFinite);
const parsePairs = (s) => String(s).trim().split(/\n+/).map((l) => l.split(/[\s,;]+/).map(Number)).filter((r) => r.length >= 2 && r.slice(0, 2).every(Number.isFinite)).map((r) => [r[0], r[1]]);

export async function render(root, [suiteId, analysisId, tabId], { setCrumb }) {
  const meta = suiteMeta(suiteId);
  if (!meta) { root.append(empty('warn', 'Unknown suite', 'Pick a suite from the menu.')); return; }
  setCrumb(`${meta.n}. ${meta.short}`);
  const disposers = [], body = h('div');
  let desc, an, tab = TABS.some((t) => t[0] === tabId) ? tabId : 'run', busy = false;
  const describe = async () => { desc = await runJob({ kind: 'describe', suite: suiteId, case: state.case, up: state.up, overridesBy: Object.fromEntries(Object.entries(state.overrides).filter(([k]) => k.startsWith(suiteId + '.')).map(([k, v]) => [k.slice(suiteId.length + 1), v])) }); };
  await describe();
  an = desc.analyses.find((a) => a.id === analysisId) || desc.analyses.find((a) => a.applicable === true) || desc.analyses[0];
  const i = SUITES.findIndex((s) => s.id === suiteId), prev = SUITES[i - 1], next = SUITES[i + 1];
  const chips = h('div', { class: 'chips', role: 'tablist', 'aria-label': 'Analyses' }), tabs = h('div', { class: 'tabs', role: 'tablist' });

  root.append(
    h('div', { class: 'page-h' },
      h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, `Suite ${meta.n} of 26 · ${meta.group}`), h('h1', null, meta.title), h('p', null, desc.tagline)),
      h('div', { class: 'row' },
        h('button', { class: 'arrow', disabled: !prev, title: prev ? `Previous suite: ${prev.short}` : '', 'aria-label': 'Previous suite', onclick: () => (location.hash = `#/suite/${prev.id}`) }, icon('back', 20)),
        h('button', { class: 'arrow', disabled: !next, title: next ? `Next suite: ${next.short}` : '', 'aria-label': 'Next suite', onclick: () => (location.hash = `#/suite/${next.id}`) }, icon('fwd', 20)),
        btn('Run whole suite', () => runAll(), { ic: 'play', title: 'Run every applicable analysis in this suite in order' }))),
    chips, tabs, body);

  const select = (a) => { an = a; replaceHash(`#/suite/${suiteId}/${an.id}/${tab}`); paintChips(); paintBody(); };
  function paintChips() {
    clear(chips);
    add(chips, desc.analyses.map((a) => { const r = state.runs[`${suiteId}.${a.id}`]; return h('button', { class: `chip ${a.id === an.id ? 'on' : ''} ${a.applicable === true ? '' : 'na'}`, role: 'tab', 'aria-selected': a.id === an.id, title: a.applicable === true ? a.summary : a.applicable, onclick: () => select(a) }, r ? h('i', { class: `dot ${r.status}`, style: { display: 'inline-block', marginRight: '7px' } }) : null, a.title); }));
    clear(tabs);
    add(tabs, TABS.map(([id, label, ic]) => h('button', { class: `tab ${id === tab ? 'on' : ''}`, role: 'tab', 'aria-selected': id === tab, onclick: () => { tab = id; replaceHash(`#/suite/${suiteId}/${an.id}/${tab}`); paintChips(); paintBody(); } }, icon(ic, 16), label)));
  }
  const dispose = () => { while (disposers.length) { try { disposers.pop()(); } catch { /* ignore */ } } };
  function paintBody() { dispose(); clear(body); ({ run: tabRun, mesh: tabMesh, studies: tabStudies, vv: tabVV, spec: tabSpec, live: tabLive })[tab](); }
  const job = (kind, opts, onProgress) => runJob({ kind, suite: suiteId, analysis: an.id, case: state.case, up: state.up, overrides: getOverrides(suiteId, an.id), opts }, onProgress);
  const notApplicable = () => (an.applicable === true ? null : h('div', { class: 'note warn' }, icon('info'), h('div', null, h('b', null, 'Not applicable to the current aircraft. '), an.applicable, ' ', h('a', { href: '#/case' }, 'Change the case'), ' or choose another analysis.')));

  async function doRun(showToast = true) {
    const payload = await job('run');
    await saveRun(suiteId, an.id, payload);
    if (showToast && payload.recs.some((r) => r.severity === 'critical')) toast('Run complete — there are critical findings to review.', 'bad');
    return payload;
  }
  async function runAll() {
    if (busy) return; busy = true;
    try { let n = 0; for (const a of desc.analyses) { if (a.applicable !== true) continue; const payload = await runJob({ kind: 'run', suite: suiteId, analysis: a.id, case: state.case, up: state.up, overrides: getOverrides(suiteId, a.id) }); await saveRun(suiteId, a.id, payload); n++; } toast(`${n} analyses finished in ${meta.short}.`, 'ok'); await describe(); paintChips(); paintBody(); }
    catch (e) { toast(e.message, 'bad'); } finally { busy = false; }
  }

  // ---------------- Set up & run ----------------
  function tabRun() {
    const na = notApplicable(), resHost = h('div', { class: 'stack' }), prog = h('progress', { max: 1, value: 0, hidden: true }), msg = h('span', { class: 'muted small' });
    const runBtn = btn('Run analysis', () => go(), { ic: 'play', kind: 'primary big', disabled: !!na });
    const cancelBtn = btn('Stop', () => { cancelAll(); }, { ic: 'stop', kind: 'ghost' }); cancelBtn.hidden = true;
    let panel = null;
    const show = (payload) => { panel?.destroy(); clear(resHost); if (!payload) { resHost.append(empty('play', 'No results yet', na ? 'This analysis does not apply to the current aircraft.' : 'Check the inputs on the left, then press Run. Values tinted blue come from your case or from upstream suites.')); return; } panel = resultPanel(payload, { suiteLabel: meta.short }); const r = state.runs[`${suiteId}.${an.id}`]; resHost.append(h('div', { class: 'row small muted' }, icon('clock', 15), `Last run ${ago(r?.ts)} · ${payload.res.elapsed_ms} ms`), panel.el); };
    disposers.push(() => panel?.destroy());
    async function go() {
      if (busy) return; busy = true; runBtn.disabled = true; prog.hidden = false; cancelBtn.hidden = false; prog.removeAttribute('value'); msg.textContent = 'Solving…';
      try { const p = await job('run', null, (f, m) => { prog.value = f; msg.textContent = m || 'Solving…'; }); await saveRun(suiteId, an.id, p); show(p); msg.textContent = ''; paintChipsOnly(); if (p.recs.some((r) => r.severity === 'critical')) toast('Run complete — critical findings need attention.', 'bad'); }
      catch (e) { msg.textContent = ''; if (e.message !== 'Cancelled') { toast(`The analysis failed: ${e.message}`, 'bad'); resHost.prepend(h('div', { class: 'note bad' }, icon('warn'), h('div', null, h('b', null, 'The solver stopped. '), e.message, ' Check the inputs for values outside their valid range.'))); } }
      finally { busy = false; runBtn.disabled = !!na; prog.hidden = true; cancelBtn.hidden = true; }
    }
    const paintChipsOnly = () => paintChips();
    const autoRun = debounce(() => { if (state.settings.autoRun && !na) go(); }, 500);
    const form = inputForm(an, suiteId, { onChange: autoRun, onReset: async () => { clearOverrides(suiteId, an.id); await describe(); an = desc.analyses.find((a) => a.id === an.id); paintBody(); } });
    body.append(h('div', { class: 'stack' },
      h('div', { class: 'note' }, icon('info'), h('div', null, h('b', null, an.title + '. '), an.summary, ' ', badge(FID[an.fidelity] || an.fidelity, 'accent'))),
      na,
      h('div', { class: 'split' },
        card('Inputs, initial and boundary conditions', h('div', { class: 'stack' }, form,
          h('div', { class: 'row' }, runBtn, cancelBtn, h('label', { class: 'switch', title: 'Re-run automatically whenever an input changes' }, h('input', { type: 'checkbox', checked: !!state.settings.autoRun, onchange: (e) => setSetting('autoRun', e.target.checked) }), 'Auto-run')),
          prog, msg)),
        resHost)));
    loadRun(suiteId, an.id).then(show);
  }

  // ---------------- Mesh & convergence ----------------
  function tabMesh() {
    const c = an.convergence, numerics = an.inputs.filter((f) => f.group === 'Numerics');
    if (!c) { body.append(h('div', { class: 'stack' }, empty('layers', 'This analysis has no mesh or step size', 'It is evaluated in closed form, so there is no discretisation error to estimate. Analyses in this suite marked “Numerical solver” have a resolution parameter you can study here.'), numerics.length ? card('Numerical settings', inputForm({ ...an, inputs: numerics }, suiteId, {})) : null)); return; }
    const out = h('div', { class: 'stack' }), levels = h('input', { class: 'inp', value: c.levels.join(', '), 'aria-label': 'Resolution levels' }), target = h('input', { class: 'inp', type: 'number', value: 1, min: 0.001, step: 0.1, 'aria-label': 'Target error percent' });
    const metric = metricSelect(c.metric), prog = h('progress', { max: 1, value: 0, hidden: true }), plots = [];
    disposers.push(() => plots.forEach((p) => p.destroy()));
    const go = async () => {
      const lv = parseList(levels.value).filter((v) => v > 0); if (lv.length < 3) { toast('Enter at least three resolution levels.', 'bad'); return; }
      prog.hidden = false; clear(out); plots.splice(0).forEach((p) => p.destroy());
      try {
        const s = await job('convergence', { levels: lv, metric: metric.value, target: Number(target.value) / 100 }, (f) => (prog.value = f)), g = s.gci;
        const sel = s.rows.find((r) => r.level === s.selected);
        out.append(
          h('div', { class: `note ${g?.monotonic ? 'ok' : 'warn'}` }, icon(g?.monotonic ? 'check' : 'warn'), h('div', null, h('b', null, `Recommended resolution: ${num(s.selected)}. `), s.selectedReason, ' ', s.verdict || '', ' ', btn('Use this resolution', () => { setOverride(suiteId, an.id, c.param, s.selected); toast(`${c.label} set to ${s.selected} for this analysis.`, 'ok'); }, { kind: 'sm primary' }))),
          h('div', { class: 'kpis' },
            kp('Observed order of accuracy', g?.p), kp('Grid convergence index (finest)', g ? 100 * g.gciFine : NaN, '%'), kp('Extrapolated value', s.reference), kp('Asymptotic-range indicator', g?.asymptotic, '', 'Close to 1 means the solutions are in the asymptotic range'), kp('Cost at recommended level', sel?.ms, 'ms')),
          h('div', { class: 'plots' }));
        const grid = out.lastChild;
        plots.push(renderPlot(grid, { type: 'line', title: `${metric.selectedOptions[0].textContent} versus resolution`, xlabel: c.label, ylabel: metric.value, xlog: true, series: [{ name: 'Computed', x: s.rows.map((r) => r.level), y: s.rows.map((r) => r.value), style: 'line+points' }], annotations: Number.isFinite(s.reference) ? [{ y: s.reference, label: s.referenceKind }] : [] }));
        if (s.rows.some((r) => r.err > 0)) plots.push(renderPlot(grid, { type: 'line', title: 'Estimated discretisation error', xlabel: 'Relative cell size h', ylabel: 'Relative error', xlog: true, ylog: true, series: [{ name: 'Error estimate', x: s.rows.filter((r) => r.err > 0).map((r) => r.h), y: s.rows.filter((r) => r.err > 0).map((r) => r.err), style: 'line+points' }], annotations: [{ y: s.target, label: 'Target' }] }));
        out.append(card(null, dataTable({ title: 'Refinement levels', columns: [c.label, 'Relative size h', metric.value, 'Estimated error [%]', 'Solve time [ms]'], rows: s.rows.map((r) => [r.level, r.h, r.value, r.err != null ? 100 * r.err : null, r.ms]) })));
      } catch (e) { out.append(h('div', { class: 'note bad' }, icon('warn'), h('div', null, e.message))); } finally { prog.hidden = true; }
    };
    body.append(h('div', { class: 'stack' }, notApplicable(),
      card('Mesh sensitivity, testing and selection', h('div', { class: 'stack' },
        h('p', { class: 'muted' }, `The solution is recomputed at several values of “${c.label}”. Richardson extrapolation gives the observed order of accuracy and the Grid Convergence Index (Roache, safety factor 1.25); the coarsest resolution that meets your error target is then selected, so you do not pay for accuracy you do not need.`),
        h('div', { class: 'grid g3' }, lab(`Levels of “${c.label}”`, levels), lab('Quantity to converge', metric), lab('Target discretisation error [%]', target)),
        h('div', { class: 'row' }, btn('Run convergence study', go, { ic: 'play', kind: 'primary', disabled: an.applicable !== true }), prog))),
      numerics.length ? card('Numerical settings for this analysis', inputForm({ ...an, inputs: numerics }, suiteId, {}), { collapsible: true, open: false }) : null, out));
  }
  const kp = (label, value, unit = '', note = '') => h('div', { class: 'kpi', title: note }, h('span', { class: 'l' }, label), h('div', { class: 'v' }, num(value), unit ? h('small', null, unit) : null));
  const lab = (text, control) => h('label', { class: 'stack small', style: { gap: '4px' } }, h('span', { class: 'muted' }, text), control);
  /** <select> of numeric outputs of the last run (falls back to a given key). */
  function metricSelect(preferred) {
    const sel = h('select', { class: 'inp', 'aria-label': 'Output quantity' }), last = state.results[`${suiteId}.${an.id}`];
    const opts = last ? Object.entries(last.res.outputs).filter(([, v]) => typeof v === 'number' && Number.isFinite(v)).map(([k]) => [k, last.res.kpis.find((x) => x.key === k)?.label || k]) : preferred ? [[preferred, preferred]] : [];
    if (preferred && !opts.some((o) => o[0] === preferred)) opts.unshift([preferred, preferred]);
    add(sel, opts.map(([k, l]) => h('option', { value: k, selected: k === preferred }, l)));
    return sel;
  }
  const numericInputs = () => an.inputs.filter((f) => (f.type || 'number') === 'number' && !f.discrete && an.values[f.key] !== 0);
  async function ensureRun() { if (!state.results[`${suiteId}.${an.id}`]) { const p = await loadRun(suiteId, an.id); if (!p) await doRun(false); } }

  // ---------------- Studies: sensitivity, sweep, uncertainty ----------------
  function tabStudies() {
    const na = notApplicable(); if (na) { body.append(na); return; }
    const host = h('div', { class: 'stack' }); body.append(host); host.append(h('div', { class: 'row muted' }, h('i', { class: 'spin' }), 'Preparing…'));
    const plots = []; disposers.push(() => plots.forEach((p) => p.destroy()));
    ensureRun().then(() => {
      clear(host);
      const inputs = numericInputs(), first = state.results[`${suiteId}.${an.id}`]?.res.kpis[0]?.key;
      // sensitivity
      const sM = metricSelect(first), sOut = h('div', { class: 'stack' }), sP = h('progress', { max: 1, value: 0, hidden: true });
      const sens = async () => { sP.hidden = false; clear(sOut); try { const r = await job('sensitivity', { metric: sM.value, pct: 0.05 }, (f) => (sP.value = f)); const top = r.rows.slice(0, 12); if (!top.length) { sOut.append(h('p', { class: 'muted' }, 'No input changes this output.')); return; } const g = h('div', { class: 'plots' }); sOut.append(g); plots.push(renderPlot(g, { type: 'bar', title: 'Sensitivity ranking (elasticity: % change in output per % change in input)', ylabel: 'Elasticity [-]', categories: top.map((x) => an.inputs.find((f) => f.key === x.key)?.label || x.key), series: [{ name: 'Elasticity', y: top.map((x) => x.elasticity) }] })); sOut.append(h('p', { class: 'muted small' }, `The strongest driver is “${an.inputs.find((f) => f.key === top[0].key)?.label}”: a 1% increase changes the output by ${num(top[0].elasticity, 3)}%. Inputs near the top deserve the most care in measurement and the tightest tolerances.`)); } catch (e) { sOut.append(h('div', { class: 'note bad' }, e.message)); } finally { sP.hidden = true; } };
      // sweep
      const wK = h('select', { class: 'inp' }, inputs.map((f) => h('option', { value: f.key }, f.label))), wM = metricSelect(first), wFrom = h('input', { class: 'inp', type: 'number' }), wTo = h('input', { class: 'inp', type: 'number' }), wN = h('input', { class: 'inp', type: 'number', value: 15, min: 3, max: 80 }), wOut = h('div', { class: 'stack' }), wP = h('progress', { max: 1, value: 0, hidden: true });
      const setRange = () => { const f = an.inputs.find((x) => x.key === wK.value), v = an.values[wK.value]; wFrom.value = Number((f.min != null ? Math.max(f.min, v * 0.5) : v * 0.5).toPrecision(4)); wTo.value = Number((f.max != null ? Math.min(f.max, v * 1.5) : v * 1.5).toPrecision(4)); };
      if (inputs.length) setRange(); wK.addEventListener('change', setRange);
      const sweep = async () => { const a = Number(wFrom.value), b = Number(wTo.value), n = Math.max(3, Math.min(80, Number(wN.value) | 0)); wP.hidden = false; clear(wOut); try { const rows = await job('sweep', { key: wK.value, values: Array.from({ length: n }, (_, k) => a + ((b - a) * k) / (n - 1)) }, (f) => (wP.value = f)); const f = an.inputs.find((x) => x.key === wK.value), g = h('div', { class: 'plots' }); wOut.append(g); plots.push(renderPlot(g, { type: 'line', title: `${wM.selectedOptions[0].textContent} versus ${f.label}`, xlabel: `${f.label}${f.unit ? ` [${f.unit}]` : ''}`, ylabel: wM.value, series: [{ name: wM.selectedOptions[0].textContent, x: rows.map((r) => r.x), y: rows.map((r) => r.out[wM.value] ?? NaN), style: 'line+points' }], annotations: [{ x: an.values[wK.value], label: 'Current' }] })); const failed = rows.filter((r) => r.error).length; if (failed) wOut.append(h('p', { class: 'muted small' }, `${failed} point(s) lay outside the model's valid range and are left blank.`)); } catch (e) { wOut.append(h('div', { class: 'note bad' }, e.message)); } finally { wP.hidden = true; } };
      // uncertainty
      const uChecks = inputs.map((f, k) => ({ f, cb: h('input', { type: 'checkbox', checked: k < 4 }), cov: h('input', { class: 'inp', type: 'number', value: 5, min: 0.1, max: 60, step: 0.5, style: { width: '74px' }, 'aria-label': `Variation of ${f.label} in percent` }) }));
      const uN = h('input', { class: 'inp', type: 'number', value: 200, min: 20, max: 5000 }), uMeth = h('select', { class: 'inp' }, h('option', { value: 'lhs' }, 'Latin hypercube'), h('option', { value: 'mc' }, 'Monte Carlo')), uM = metricSelect(first), uOut = h('div', { class: 'stack' }), uP = h('progress', { max: 1, value: 0, hidden: true });
      const uq = async () => { const vars = uChecks.filter((x) => x.cb.checked).map((x) => ({ key: x.f.key, type: 'normal', cov: Number(x.cov.value) / 100 })); if (!vars.length) { toast('Tick at least one uncertain input.', 'bad'); return; } uP.hidden = false; clear(uOut); try { const r = await job('uq', { vars, n: Number(uN.value) | 0, method: uMeth.value, seed: 2024 }, (f) => (uP.value = f)), s = r.stats[uM.value]; if (!s) { uOut.append(h('p', { class: 'muted' }, 'This output does not vary with the chosen inputs.')); return; } uOut.append(h('div', { class: 'kpis' }, kp('Mean', s.mean), kp('Standard deviation', s.sd), kp('Coefficient of variation', 100 * s.cov, '%'), kp('5th percentile', s.p05), kp('Median', s.p50), kp('95th percentile', s.p95), kp('95% confidence half-width of the mean', s.ci95))); const g = h('div', { class: 'plots' }); uOut.append(g); plots.push(renderPlot(g, { type: 'bar', title: `Distribution of ${uM.selectedOptions[0].textContent}`, ylabel: 'Samples', categories: s.hist.centers.map((c) => num(c, 3)), series: [{ name: 'Count', y: s.hist.counts }] })); plots.push(renderPlot(g, { type: 'bar', title: 'Which uncertain input matters most (correlation with the output)', ylabel: 'Correlation coefficient [-]', categories: s.importance.map((x) => an.inputs.find((f) => f.key === x.key)?.label || x.key), series: [{ name: 'Correlation', y: s.importance.map((x) => x.r) }] })); uOut.append(h('p', { class: 'muted small' }, `${r.n} samples (${r.failed} outside the valid range). With these input uncertainties there is a 90% chance the result lies between ${num(s.p05)} and ${num(s.p95)}.`)); } catch (e) { uOut.append(h('div', { class: 'note bad' }, e.message)); } finally { uP.hidden = true; } };
      host.append(
        card('Sensitivity ranking — which inputs drive the answer', h('div', { class: 'stack' }, h('div', { class: 'row' }, lab('Output', sM), btn('Rank inputs', sens, { ic: 'play', kind: 'primary' })), sP, sOut), { collapsible: true }),
        card('Parameter sweep — how the answer changes with one input', h('div', { class: 'stack' }, h('div', { class: 'grid g4' }, lab('Input to vary', wK), lab('From', wFrom), lab('To', wTo), lab('Points', wN), lab('Output', wM)), h('div', { class: 'row' }, btn('Run sweep', sweep, { ic: 'play', kind: 'primary' })), wP, wOut), { collapsible: true }),
        card('Uncertainty propagation — how sure is the answer', h('div', { class: 'stack' },
          h('p', { class: 'muted small' }, 'Tick the inputs that are uncertain and give each a variation (standard deviation as a percentage of its value). The model is then sampled to give the spread of the result.'),
          h('div', { class: 'grid g3' }, uChecks.map((x) => h('label', { class: 'row small' }, x.cb, h('span', { class: 'grow' }, x.f.label), x.cov, '%'))),
          h('div', { class: 'grid g4' }, lab('Samples', uN), lab('Sampling', uMeth), lab('Output', uM)), h('div', { class: 'row' }, btn('Propagate uncertainty', uq, { ic: 'play', kind: 'primary' })), uP, uOut), { collapsible: true }));
    }).catch((e) => { clear(host); host.append(h('div', { class: 'note bad' }, e.message)); });
  }

  // ---------------- Verification, validation, calibration ----------------
  function tabVV() {
    const vOut = h('div', { class: 'stack' }), plots = []; disposers.push(() => plots.forEach((p) => p.destroy()));
    const verify = async () => { clear(vOut); vOut.append(h('div', { class: 'row muted' }, h('i', { class: 'spin' }), 'Running benchmarks…')); try { const checks = await runJob({ kind: 'verify', suite: suiteId }), pass = checks.filter((c) => c.pass).length; clear(vOut); vOut.append(h('div', { class: `note ${pass === checks.length ? 'ok' : 'bad'}` }, icon(pass === checks.length ? 'check' : 'warn'), h('div', null, h('b', null, `${pass} of ${checks.length} code-verification benchmarks pass. `), 'Each compares this implementation with an exact or analytical solution. Passing shows the equations are solved correctly; it does not by itself show that the model represents your aircraft.')), dataTable({ title: 'Verification benchmarks', columns: ['Analysis', 'Benchmark', 'Computed', 'Exact', 'Relative error', 'Tolerance', 'Result', 'Reference'], rows: checks.map((c) => [desc.analyses.find((a) => a.id === c.analysis)?.title || c.analysis, c.name, c.actual, c.expected, c.error, c.tol, c.pass ? 'Pass' : 'FAIL', c.ref || '']) })); } catch (e) { clear(vOut); vOut.append(h('div', { class: 'note bad' }, e.message)); } };
    const inputs = numericInputs(), cal = an.calibration, first = state.results[`${suiteId}.${an.id}`]?.res.kpis[0]?.key;
    // validation against user data
    const dK = h('select', { class: 'inp' }, inputs.map((f) => h('option', { value: f.key, selected: f.key === cal?.sweep }, f.label))), dM = metricSelect(cal?.target || first), dTol = h('input', { class: 'inp', type: 'number', value: 10, min: 0.1 }), dTxt = h('textarea', { class: 'inp', placeholder: 'One measurement per line:  input value, measured output\n0, 0.21\n4, 0.65\n8, 1.08', 'aria-label': 'Measured data' }), dOut = h('div', { class: 'stack' });
    const valRun = async (vset) => { clear(dOut); try { const r = await job('validate', vset); const m = r.metrics, g = h('div', { class: 'plots' }); dOut.append(h('div', { class: `note ${r.pass ? 'ok' : 'warn'}` }, icon(r.pass ? 'check' : 'warn'), h('div', null, h('b', null, `${r.within} of ${r.predicted.length} points agree within ${r.tol_pct ?? 10}%. `), `RMSE ${num(m.rmse)}, bias ${num(m.bias)}, mean absolute percentage error ${num(m.mape, 3)}%, R² ${num(m.r2, 3)}.`, r.source ? ` Source: ${r.source}.` : '')), g); plots.push(renderPlot(g, { type: 'line', title: 'Prediction against reference data', xlabel: an.inputs.find((f) => f.key === r.sweep.key)?.label || r.sweep.key, ylabel: r.target, series: [{ name: 'Model prediction', x: r.sweep.values, y: r.predicted }, { name: 'Reference data', x: r.sweep.values, y: r.observed, style: 'points' }] })); } catch (e) { dOut.append(h('div', { class: 'note bad' }, e.message)); } };
    const valUser = () => { const d = parsePairs(dTxt.value); if (d.length < 2) { toast('Enter at least two lines of “input, measured output”.', 'bad'); return; } valRun({ name: 'User data', sweep: { key: dK.value, values: d.map((r) => r[0]) }, target: dM.value, observed: d.map((r) => r[1]), tol_pct: Number(dTol.value) }); };
    // calibration
    const cParams = (cal?.params?.length ? cal.params.map((p) => ({ ...p, f: an.inputs.find((f) => f.key === p.key) })).filter((p) => p.f) : inputs.slice(0, 6).map((f) => ({ key: f.key, min: f.min, max: f.max, f }))).map((p, k) => ({ ...p, cb: h('input', { type: 'checkbox', checked: cal ? true : k < 1 }) }));
    const cOut = h('div', { class: 'stack' }), cP = h('progress', { max: 1, value: 0, hidden: true });
    const calRun = async () => { const d = parsePairs(dTxt.value), params = cParams.filter((p) => p.cb.checked).map((p) => ({ key: p.key, min: p.min, max: p.max })); if (d.length < params.length + 1 || !params.length) { toast('Calibration needs at least one parameter and more data points than parameters.', 'bad'); return; } cP.hidden = false; clear(cOut); try { const r = await job('calibrate', { params, sweepKey: dK.value, target: dM.value, data: d }, (f) => (cP.value = f)), g = h('div', { class: 'plots' }); cOut.append(h('div', { class: 'note ok' }, icon('check'), h('div', null, h('b', null, `RMSE reduced from ${num(r.metricsBefore.rmse)} to ${num(r.metricsAfter.rmse)}. `), 'Keep a separate dataset for validation: agreement with the data used for fitting does not demonstrate predictive accuracy. ', btn('Apply fitted values', () => { r.params.forEach((p) => setOverride(suiteId, an.id, p.key, Number(p.fitted.toPrecision(6)))); toast('Calibrated values applied to this analysis.', 'ok'); }, { kind: 'sm primary' }))), dataTable({ title: 'Calibrated parameters', columns: ['Parameter', 'Initial', 'Fitted', 'Standard error'], rows: r.params.map((p) => [an.inputs.find((f) => f.key === p.key)?.label || p.key, p.initial, p.fitted, p.se]) }), g); plots.push(renderPlot(g, { type: 'line', title: 'Model before and after calibration', xlabel: an.inputs.find((f) => f.key === dK.value)?.label, ylabel: dM.value, series: [{ name: 'Before', x: r.x, y: r.before, style: 'dash' }, { name: 'After', x: r.x, y: r.after }, { name: 'Measured', x: r.x, y: r.observed, style: 'points' }] })); } catch (e) { cOut.append(h('div', { class: 'note bad' }, e.message)); } finally { cP.hidden = true; } };
    const specHost = h('div', { class: 'spec stack' });
    loadSpec().then((sp) => { const s = sp.suites[meta.n]; add(specHost, [['Calibration', s.calibration], ['Verification', s.verification], ['Validation', s.validation]].map(([t, x]) => h('div', null, h('h4', null, t + ' requirement'), h('p', null, x)))); }).catch(() => {});
    body.append(h('div', { class: 'stack' },
      card('Code verification — is the mathematics solved correctly?', h('div', { class: 'stack' }, h('div', { class: 'row' }, btn('Run verification benchmarks', verify, { ic: 'check', kind: 'primary' }), h('span', { class: 'muted small' }, 'Compares every solver in this suite with exact solutions.')), vOut)),
      an.validation.length ? card('Built-in reference cases', h('div', { class: 'stack' }, an.validation.map((v) => h('div', { class: 'row' }, h('span', { class: 'grow' }, h('b', null, v.name), h('span', { class: 'muted small' }, ` — ${v.source || ''}`)), btn('Compare', () => valRun(v), { kind: 'sm', disabled: an.applicable !== true }))))) : null,
      an.applicable === true ? card('Validate and calibrate against your own measurements', h('div', { class: 'stack' },
        h('p', { class: 'muted small' }, cal?.note || 'Paste test, rig or flight data. Validation compares the model with the data as it stands; calibration adjusts the ticked parameters to fit it. Use different datasets for the two.'),
        h('div', { class: 'grid g3' }, lab('Varied input (first column)', dK), lab('Measured output (second column)', dM), lab('Acceptance tolerance [%]', dTol)), dTxt,
        h('div', { class: 'row' }, btn('Validate', valUser, { ic: 'check', kind: 'primary' }), btn('Calibrate', calRun, { ic: 'sliders' }), btn('Load a CSV file', async () => { const { pickFiles } = await import('../dom.js'); const [f] = await pickFiles({ accept: '.csv,.txt,.dat' }); if (f) dTxt.value = await f.text(); }, { ic: 'upload', kind: 'ghost' })),
        h('div', { class: 'row small' }, h('span', { class: 'muted' }, 'Parameters to calibrate:'), cParams.map((p) => h('label', { class: 'row' }, p.cb, p.f.label))), cP, dOut, cOut)) : notApplicable(),
      card('What the specification requires for this suite', specHost, { collapsible: true, open: false })));
  }

  // ---------------- Specification ----------------
  function tabSpec() {
    const host = h('div', { class: 'stack' }); body.append(host); host.append(h('div', { class: 'row muted' }, h('i', { class: 'spin' }), 'Loading specification…'));
    loadSpec().then((sp) => {
      const s = sp.suites[meta.n], impl = desc.analyses.flatMap((a) => a.equations).map((e) => e.toLowerCase());
      const hit = (t) => { const x = t.toLowerCase().replace(/ equations?| relations?| formulations?/g, ''); return impl.some((e) => e.includes(x.slice(0, 22)) || x.includes(e.replace(/ equations?| relations?| formulations?/g, '').slice(0, 22))); };
      const termBlock = (title, text) => { const ts = terms(text); return h('div', null, h('h4', null, title), ts.length > 2 ? h('div', { class: 'terms', style: { marginTop: '6px' } }, ts.map((t) => h('span', { class: `term ${hit(t) ? 'on' : ''}`, title: hit(t) ? 'Used by an analysis in this suite' : '' }, t))) : h('p', null, text)); };
      clear(host);
      host.append(
        h('div', { class: 'note' }, icon('info'), h('div', null, h('b', null, 'Scope: '), s.scope, '. Items highlighted in green are named by an analysis in this suite as part of what it solves.')),
        card('Governing equations and models', h('div', { class: 'spec stack' }, termBlock('Classical governing equations and formulations', s.classical), termBlock('Hybrid equations and computational formulations', s.hybrid), termBlock('Required computational models', s.models), ...(s.notes || []).map((n) => h('p', null, n)))),
        card('How each analysis here is solved', dataTable({ title: '', columns: ['Analysis', 'Method', 'Equations and models it implements'], rows: desc.analyses.map((a) => [a.title, FID[a.fidelity] || a.fidelity, a.equations.join('; ')]) })),
        desc.handoff.length ? card('Models handed off to external high-fidelity solvers', h('div', { class: 'stack' }, h('p', { class: 'muted small' }, 'These models from the specification are not solved inside the browser. The suite prepares consistent inputs and reduced-order estimates; use the class of tool named here when you need that fidelity, then bring the results back through the calibration and validation tab.'), dataTable({ title: '', columns: ['Model', 'Why it is handed off', 'Recommended tool class'], rows: desc.handoff.map((x) => [x.model, x.why, x.tool]) }))) : null,
        card('Initial and boundary conditions', h('div', { class: 'spec stack' }, h('div', null, h('h4', null, 'Initial conditions'), h('p', null, s.ic)), h('div', null, h('h4', null, 'Boundary conditions'), h('p', null, s.bc))), { collapsible: true, open: false }),
        card('Input and output data', h('div', { class: 'spec stack' }, termBlock('Accepted input data', s.inputs), termBlock('Generated output data', s.outputs)), { collapsible: true, open: false }),
        card('Geometry and mesh formats', h('div', { class: 'spec stack' }, h('p', null, s.geometry), h('div', null, btn('Open the geometry & mesh workbench', () => (location.hash = '#/geometry'), { ic: 'cube' }))), { collapsible: true, open: false }),
        card('Data exchanged with other suites', h('div', { class: 'grid g2' },
          h('div', null, h('h4', null, 'Uses'), desc.consumes.length ? h('ul', { class: 'small' }, desc.consumes.map((c) => h('li', null, h('a', { href: `#/suite/${c.from}` }, suiteMeta(c.from)?.short || c.from), `: ${c.keys.join(', ')} — ${c.why || ''}`))) : h('p', { class: 'muted small' }, 'Only the shared case.')),
          h('div', null, h('h4', null, 'Publishes'), h('ul', { class: 'small' }, desc.provides.map((p) => h('li', null, `${p.label}${p.unit && p.unit !== '-' ? ` [${p.unit}]` : ''}`, state.up[suiteId]?.[p.key] != null ? h('b', null, ` = ${num(state.up[suiteId][p.key])}`) : null))))), { collapsible: true, open: false }));
    }).catch((e) => { clear(host); host.append(h('div', { class: 'note bad' }, e.message)); });
  }

  // ---------------- Live resources ----------------
  function tabLive() {
    const [qLit, qCode] = SUITE_QUERIES[suiteId] || [meta.title, meta.short], lit = h('div', { class: 'stack' }), code = h('div', { class: 'stack' });
    let sort = 'recent';
    const stamp = (r) => h('p', { class: 'muted small' }, r.error ? `Showing saved results (${r.error.toLowerCase()}). ` : '', r.ts ? `Updated ${ago(r.ts)}.` : '');
    const loadLit = async (force) => { clear(lit); lit.append(h('div', { class: 'row muted' }, h('i', { class: 'spin' }), 'Searching…')); const r = await live('literature', { q: qLit, sort }, { force }); clear(lit); if (!r.data) { lit.append(h('div', { class: 'note warn' }, icon('info'), h('div', null, `No connection to the literature index (${r.error}). Results appear here when you are online; they are then kept for offline use.`))); return; } lit.append(stamp(r), ...r.data.items.map((w) => h('div', null, h('a', { href: w.url, target: '_blank', rel: 'noopener noreferrer' }, w.title), ' ', w.oa ? badge('Open access', 'ok') : null, h('div', { class: 'muted small' }, [w.authors, w.venue, w.date, `${w.cited} citations`].filter(Boolean).join(' · '))))); };
    const loadCode = async (force) => { clear(code); code.append(h('div', { class: 'row muted' }, h('i', { class: 'spin' }), 'Searching…')); const r = await live('opensource', { q: qCode }, { force }); clear(code); if (!r.data) { code.append(h('div', { class: 'note warn' }, icon('info'), h('div', null, `No connection to the repository index (${r.error}).`))); return; } code.append(stamp(r), ...r.data.items.map((x) => h('div', null, h('a', { href: x.url, target: '_blank', rel: 'noopener noreferrer' }, x.name), ' ', badge(`★ ${num(x.stars)}`), x.lang ? badge(x.lang) : null, x.license && x.license !== 'NOASSERTION' ? badge(x.license) : null, h('div', { class: 'muted small' }, x.desc, x.pushed ? ` · updated ${x.pushed.slice(0, 10)}` : '')))); };
    body.append(h('div', { class: 'grid g2' },
      card('Research literature', lit, { actions: h('span', { class: 'row' }, h('select', { class: 'inp', style: { width: 'auto' }, 'aria-label': 'Sort', onchange: (e) => { sort = e.target.value; loadLit(); } }, h('option', { value: 'recent' }, 'Newest'), h('option', { value: 'cited' }, 'Most cited')), btn('', () => loadLit(true), { ic: 'refresh', kind: 'sm ghost', title: 'Refresh now' })) }),
      card('Open-source solvers and tools', code, { actions: btn('', () => loadCode(true), { ic: 'refresh', kind: 'sm ghost', title: 'Refresh now' }) })),
      h('p', { class: 'muted small', style: { marginTop: '10px' } }, 'Fetched directly by your browser from OpenAlex and GitHub, refreshed daily and whenever you press refresh, and kept on this device for offline reading. Listings are search results, not endorsements: check each source before relying on it.'));
    loadLit(); loadCode();
  }

  paintChips(); paintBody();
  const off = on('case', debounce(async () => { await describe(); an = desc.analyses.find((a) => a.id === an.id) || desc.analyses[0]; paintChips(); if (tab === 'run') paintBody(); }, 400));
  return () => { off(); dispose(); };
}

/** Input form for one analysis, grouped, with provenance tags and inline range checks. */
function inputForm(an, suiteId, { onChange, onReset }) {
  const ov = getOverrides(suiteId, an.id), groups = new Map();
  for (const f of an.inputs) { const g = f.group || 'Inputs'; if (!groups.has(g)) groups.set(g, []); groups.get(g).push(f); }
  const form = h('div', { class: 'form' });
  let gi = 0;
  for (const [g, fields] of groups) {
    const det = h('details', { class: 'fgroup', open: gi++ < 3 || fields.some((f) => f.key in ov) }, h('summary', null, g));
    for (const f of fields) {
      const type = f.type || 'number', v = an.values[f.key], edited = f.key in ov, linked = an.linked.includes(f.key) && !edited, err = h('div', { class: 'err', hidden: true });
      let ctl;
      const commit = (val) => { setOverride(suiteId, an.id, f.key, val); an.values[f.key] = val; ctl.classList.remove('linked'); ctl.classList.add('edited'); setKids(tagHost, resetTag()); onChange?.(); };
      const check = (val) => { let m = ''; if (type === 'number') { if (!Number.isFinite(val)) m = 'Enter a number.'; else if (f.min != null && val < f.min) m = `Below the minimum of ${num(f.min)}${f.unit && f.unit !== '-' ? ' ' + f.unit : ''}.`; else if (f.max != null && val > f.max) m = `Above the maximum of ${num(f.max)}${f.unit && f.unit !== '-' ? ' ' + f.unit : ''}.`; } err.textContent = m; err.hidden = !m; ctl.classList.toggle('invalid', !!m); return !m; };
      if (type === 'select') ctl = h('select', { class: 'inp', id: `f-${f.key}`, onchange: (e) => commit(e.target.value) }, (f.options || []).map((o) => h('option', { value: o, selected: String(o) === String(v) }, String(o))));
      else if (type === 'bool') ctl = h('input', { type: 'checkbox', id: `f-${f.key}`, checked: !!v, onchange: (e) => commit(e.target.checked) });
      else if (type === 'text') ctl = h('textarea', { class: 'inp', id: `f-${f.key}`, value: v ?? '', spellcheck: false, onchange: (e) => commit(e.target.value) });
      else ctl = h('input', { class: 'inp', id: `f-${f.key}`, type: 'number', inputMode: 'decimal', step: f.step ?? 'any', value: typeof v === 'number' ? Number(v.toPrecision(7)) : v, onchange: (e) => { const val = Number(e.target.value); if (check(val)) commit(f.discrete ? Math.round(val) : val); } , oninput: (e) => check(Number(e.target.value)) });
      if (linked) ctl.classList.add('linked'); if (edited) ctl.classList.add('edited');
      const resetTag = () => h('button', { class: 'tag ed', type: 'button', title: 'You changed this value. Click to return to the value from the case.', onclick: () => { setOverride(suiteId, an.id, f.key, undefined); onReset?.(); } }, 'edited ×');
      const tagHost = h('span', null, edited ? resetTag() : linked ? h('span', { class: 'tag', title: 'Filled automatically from your case, live data or an upstream suite' }, 'linked') : null);
      det.append(h('div', { class: `field ${type === 'text' ? 'wide' : ''}` }, h('label', { for: `f-${f.key}` }, f.label, f.unit && f.unit !== '-' ? h('span', { class: 'u' }, `[${f.unit}]`) : null, tagHost), ctl, f.help ? h('div', { class: 'help' }, f.help) : null, err));
    }
    form.append(det);
  }
  if (onReset && Object.keys(ov).length) form.append(h('div', null, btn('Reset all inputs to case values', onReset, { ic: 'refresh', kind: 'sm ghost' })));
  return form;
}
