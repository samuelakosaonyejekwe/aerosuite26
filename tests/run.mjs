// Headless test harness: `node tests/run.mjs [suiteId ...] [--quick]`
// For every suite it checks the module contract, runs each analysis on every aircraft preset,
// validates result shapes, executes the built-in verification benchmarks and the convergence studies,
// and finally runs the full coupled chain per preset.

import { existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { SUITES, loadSuite, executionOrder } from '../js/core/registry.js';
import { PRESETS } from '../js/core/case.js';
import { makeCtx, resolveInputs, runAnalysis, applicable, verifySuite, convergenceStudy } from '../js/core/studies.js';

const args = process.argv.slice(2), quick = args.includes('--quick'), only = args.filter((a) => !a.startsWith('--'));
const root = fileURLToPath(new URL('../js/suites/', import.meta.url));
let fail = 0, soft = 0, ran = 0;
const bad = (m) => { fail++; console.log('  FAIL ' + m); };
const PLOT_TYPES = ['line', 'bar', 'heat', 'polar', 'tri'];
const numArr = (a) => Array.isArray(a) && a.every((v) => typeof v === 'number');

function checkPlot(p, where) {
  if (!PLOT_TYPES.includes(p.type)) return bad(`${where}: unknown plot type ${p.type}`);
  if (!p.title) bad(`${where}: plot without title`);
  if (p.type === 'line') for (const s of p.series || []) { if (!numArr(s.x) || !numArr(s.y) || s.x.length !== s.y.length) bad(`${where}: line series "${s.name}" x/y mismatch in "${p.title}"`); if (!s.name) bad(`${where}: unnamed series in "${p.title}"`); if (s.x.length && !s.y.some(Number.isFinite)) bad(`${where}: series "${s.name}" in "${p.title}" has no finite values`); if (s.x.some((v) => Number.isNaN(v)) && !s.x.some(Number.isFinite)) bad(`${where}: series "${s.name}" x is all NaN`); }
  if (p.type === 'bar') for (const s of p.series || []) if (!numArr(s.y) || s.y.length !== p.categories.length) bad(`${where}: bar series length in "${p.title}"`);
  if (p.type === 'heat' && !p.z.flat().some(Number.isFinite)) bad(`${where}: heat map "${p.title}" has no finite values`);
  if (p.type === 'heat') { if (!numArr(p.x) || !numArr(p.y) || p.z.length !== p.y.length || p.z[0].length !== p.x.length) bad(`${where}: heat grid shape in "${p.title}"`); }
  if (p.type === 'polar') for (const s of p.series || []) if (!numArr(s.theta_deg) || s.theta_deg.length !== s.r.length) bad(`${where}: polar series in "${p.title}"`);
  if (p.type === 'tri') { if (p.values.length !== p.nodes.length) bad(`${where}: tri values/nodes in "${p.title}"`); }
}

const defs = [];
for (const meta of SUITES) {
  if (!existsSync(root + meta.file + '.js')) { if (!only.length || only.includes(meta.id)) console.log(`-- ${meta.n} ${meta.id}: (module not present)`); continue; }
  let def;
  try { def = await loadSuite(meta.id); } catch (e) { bad(`${meta.id}: import failed: ${e.stack}`); continue; }
  defs.push(def);
  if (only.length && !only.includes(meta.id)) continue;
  console.log(`== ${meta.n} ${meta.id}: ${def.analyses.length} analyses`);
  if (def.id !== meta.id || def.n !== meta.n) bad(`${meta.id}: id/n mismatch`);
  if (!def.analyses?.length) bad(`${meta.id}: no analyses`);
  const ids = new Set();
  for (const an of def.analyses) {
    if (ids.has(an.id)) bad(`${meta.id}: duplicate analysis id ${an.id}`); ids.add(an.id);
    for (const f of ['id', 'title', 'summary', 'fidelity', 'inputs', 'run']) if (!an[f]) bad(`${meta.id}.${an.id}: missing ${f}`);
    const keys = new Set();
    for (const f of an.inputs || []) { if (keys.has(f.key)) bad(`${meta.id}.${an.id}: duplicate input ${f.key}`); keys.add(f.key); if (!f.label) bad(`${meta.id}.${an.id}: input ${f.key} has no label`); if (f.default === undefined) bad(`${meta.id}.${an.id}: input ${f.key} has no default`); }
  }
  // verification benchmarks
  const checks = verifySuite(def);
  for (const c of checks) if (!c.pass) bad(`${meta.id}.${c.analysis} verification "${c.name}": got ${c.actual}, expected ${c.expected} (tol ${c.tol})`);
  console.log(`  verification: ${checks.filter((c) => c.pass).length}/${checks.length} pass`);
  if (!checks.length) bad(`${meta.id}: no verification benchmarks`);
  // every analysis on every preset
  for (const [pk, preset] of Object.entries(PRESETS)) {
    const ctx = makeCtx(preset.data());
    for (const an of def.analyses) {
      const ok = applicable(an, ctx);
      if (ok !== true) continue;
      const where = `${meta.id}.${an.id}[${pk}]`;
      try {
        const { inp } = resolveInputs(an, ctx), t = Date.now(), res = await runAnalysis(an, inp, ctx); ran++;
        if (!res.kpis.length) bad(`${where}: no KPIs`);
        if (!res.plots.length && !res.tables.length) bad(`${where}: no plots or tables`);
        for (const k of res.kpis) { if (typeof k.value !== 'number') bad(`${where}: KPI ${k.key} is not a number`); else if (!Number.isFinite(k.value)) { soft++; if (!quick) console.log(`  note ${where}: KPI ${k.key} = ${k.value}`); } if (!k.key || !k.label) bad(`${where}: KPI missing key/label`); }
        res.plots.forEach((p) => checkPlot(p, where));
        for (const tb of res.tables) if (!tb.columns || tb.rows.some((r) => r.length !== tb.columns.length)) bad(`${where}: table "${tb.title}" shape`);
        try { structuredClone(res); } catch { bad(`${where}: result is not structured-cloneable (functions/classes in result)`); }
        if (Date.now() - t > 4000) bad(`${where}: too slow (${Date.now() - t} ms)`);
        if (an.recommend) { const rec = an.recommend(res, inp, ctx) || []; for (const r of rec) if (!r.title || !r.action || !r.severity) bad(`${where}: malformed recommendation`); }
      } catch (e) { bad(`${where}: threw ${e.stack?.split('\n').slice(0, 3).join(' | ')}`); }
    }
  }
  // convergence studies on the first applicable preset
  if (!quick) for (const an of def.analyses) if (an.convergence) {
    const pk = Object.keys(PRESETS).find((k) => applicable(an, makeCtx(PRESETS[k].data())) === true); if (!pk) continue;
    const ctx = makeCtx(PRESETS[pk].data());
    try { const s = await convergenceStudy(an, resolveInputs(an, ctx).inp, ctx); console.log(`  convergence ${an.id}: order ${s.gci?.p?.toFixed(2)}, GCI ${(100 * s.gci?.gciFine).toPrecision(2)}%, selected ${s.selected}`); if (s.rows.some((r) => !Number.isFinite(r.value))) bad(`${meta.id}.${an.id}: convergence metric not finite`); }
    catch (e) { bad(`${meta.id}.${an.id}: convergence threw ${e.message}`); }
  }
}

// integrated chain
if (!only.length || args.includes('--chain')) {
  const order = executionOrder(defs);
  console.log('== integrated order: ' + order.join(' > '));
  for (const [pk, preset] of Object.entries(PRESETS)) {
    const up = {}, c = preset.data(); let n = 0, t = Date.now();
    for (const id of order) {
      const def = defs.find((d) => d.id === id); up[id] = up[id] || {};
      for (const an of def.analyses) {
        const ctx = makeCtx(c, up); if (applicable(an, ctx) !== true) continue;
        try { const res = await runAnalysis(an, resolveInputs(an, ctx).inp, ctx); Object.assign(up[id], res.outputs); n++; } catch (e) { bad(`chain[${pk}] ${id}.${an.id}: ${e.message}`); }
      }
    }
    for (const def of defs) for (const p of def.provides || []) if (!(p.key in (up[def.id] || {})) && def.analyses.some((a) => applicable(a, makeCtx(c)) === true)) { soft++; if (!quick) console.log(`  note chain[${pk}]: ${def.id} did not publish "${p.key}"`); }
    console.log(`  chain[${pk}]: ${n} analyses in ${Date.now() - t} ms`);
  }
}
console.log(`\n${ran} analysis runs, ${soft} soft notes, ${fail} failures`);
process.exit(fail ? 1 : 0);
