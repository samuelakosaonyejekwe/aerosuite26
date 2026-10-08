// Reports & assurance: traceable engineering report, platform-wide verification evidence, data export.

import { h, icon, clear, add, num, btn, badge, card, toast, downloadText, toCsv, pickFiles } from '../dom.js';
import { SUITES, suiteMeta } from '../../core/registry.js';
import { state, loadRun, exportProject, importProject, clearRuns, on } from '../../core/store.js';
import { runJob } from '../../core/runner.js';
import { dataTable, kpiGrid, recList } from '../results.js';
import { CASE_FIELDS, CASE_SECTIONS } from '../../core/case.js';
import { loadSpec } from '../../core/spec.js';

export async function render(root, _p, { setCrumb }) {
  setCrumb('Reports & assurance');
  const host = h('div'); root.append(host);
  const vOut = h('div', { class: 'stack' }), vProg = h('progress', { max: 1, value: 0, hidden: true }), reportHost = h('div', { class: 'stack' });

  const verifyAll = async () => {
    vProg.hidden = false; clear(vOut);
    try {
      const res = await runJob({ kind: 'verifyAll' }, (f) => (vProg.value = f)), flat = res.flatMap((s) => s.checks.map((c) => ({ ...c, suite: s.suite }))), pass = flat.filter((c) => c.pass).length;
      vOut.append(h('div', { class: `note ${pass === flat.length ? 'ok' : 'bad'}` }, icon(pass === flat.length ? 'check' : 'warn'), h('div', null, h('b', null, `${pass} of ${flat.length} verification benchmarks pass across ${res.length} suites. `), 'Run on this device, just now, against exact and analytical solutions.')),
        dataTable({ title: 'Platform verification evidence', columns: ['Suite', 'Benchmark', 'Computed', 'Exact', 'Relative error', 'Tolerance', 'Result', 'Reference'], rows: flat.map((c) => [`${suiteMeta(c.suite).d}. ${suiteMeta(c.suite).short}`, c.name, c.actual, c.expected, c.error, c.tol, c.pass ? 'Pass' : 'FAIL', c.ref || '']) }, { max: 600 }));
    } catch (e) { vOut.append(h('div', { class: 'note bad' }, e.message)); } finally { vProg.hidden = true; }
  };

  const allKpis = async () => { const rows = []; for (const k of Object.keys(state.runs)) { const [s, a] = k.split('.'), p = await loadRun(s, a); if (!p) continue; for (const x of p.res.kpis) rows.push([suiteMeta(s).d, suiteMeta(s).short, a, x.key, x.label, x.value, x.unit || '', x.status || '', new Date(state.runs[k].ts).toISOString()]); } return rows.sort((x, y) => x[0] - y[0]); };

  const buildReport = async () => {
    clear(reportHost); reportHost.append(h('div', { class: 'row muted' }, h('i', { class: 'spin' }), 'Assembling report…'));
    const c = state.case, keys = Object.keys(state.runs).sort((a, b) => suiteMeta(a.split('.')[0]).d - suiteMeta(b.split('.')[0]).d);
    let spec = null; try { spec = await loadSpec(); } catch { /* optional */ }
    const parts = [h('div', null, h('div', { class: 'eyebrow' }, 'AeroSuite 26 engineering report'), h('h1', null, c.meta.name || 'Unnamed case'), h('p', { class: 'muted' }, `${c.meta.type} · generated ${new Date().toLocaleString()} · site: ${c.site.name} (${c.site.source})`))];
    parts.push(card('Case definition', h('div', { class: 'grid g3' }, Object.entries(CASE_SECTIONS).map(([sec, title]) => h('div', null, h('h4', null, title), h('dl', { class: 'kv small' }, CASE_FIELDS.filter((f) => f[0] === sec).filter((f) => c[sec][f[1]] !== 0 && c[sec][f[1]] !== '').map((f) => [h('dt', null, f[2]), h('dd', null, typeof c[sec][f[1]] === 'number' ? `${num(c[sec][f[1]])} ${f[3] || ''}` : String(c[sec][f[1]]))])))))));
    let lastSuite = null;
    for (const k of keys) {
      const [s, a] = k.split('.'), p = await loadRun(s, a); if (!p) continue; const m = suiteMeta(s);
      if (s !== lastSuite) { parts.push(h('h2', { style: { marginTop: '18px' } }, `${m.d}. ${m.title}`)); if (spec) parts.push(h('p', { class: 'muted small' }, spec.suites[m.n].scope)); lastSuite = s; }
      parts.push(card(a, h('div', { class: 'stack' }, kpiGrid(p.res.kpis),
        p.res.warnings.length ? h('div', { class: 'note warn' }, icon('warn'), h('div', null, p.res.warnings.join(' '))) : null,
        p.recs.length ? recList(p.recs) : null,
        h('p', { class: 'small muted' }, h('b', null, 'Models: '), p.res.models.join('; ') || 'not stated', h('br'), h('b', null, 'Assumptions: '), p.res.assumptions.join('; ') || 'none stated', h('br'), h('b', null, 'Inputs: '), Object.entries(p.inp).filter(([, v]) => typeof v !== 'string' || v.length < 40).map(([key, v]) => `${key} = ${typeof v === 'number' ? num(v) : v}`).join(', ')))));
    }
    if (!keys.length) parts.push(h('p', { class: 'muted' }, 'No analyses have been run for this case yet.'));
    parts.push(h('div', { class: 'note' }, icon('shield'), h('div', null, 'Results are computed predictions from the stated models and inputs. They distinguish calculated values from empirical estimates and assumed inputs, and they are not experimentally validated for this aircraft unless validation data were supplied in the relevant suite. This report is not, on its own, evidence of airworthiness.')));
    clear(reportHost); add(reportHost, parts);
    reportHost.prepend(h('div', { class: 'row no-print' }, btn('Print or save as PDF', () => window.print(), { ic: 'doc', kind: 'primary' })));
  };

  host.append(
    h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Traceable evidence'), h('h1', null, 'Reports & assurance'), h('p', null, 'Export what you have done in forms others can check: a complete project file that reproduces every run, tables of all results, a printable engineering report, and verification evidence for the solvers themselves.'))),
    h('div', { class: 'grid g2 no-print' },
      card('Project and results', h('div', { class: 'stack' },
        h('div', { class: 'row' }, btn('Save project', () => downloadText(`${(state.case.meta.name || 'case').toLowerCase().replace(/[^a-z0-9]+/g, '-')}.aerosuite.json`, JSON.stringify(exportProject(), null, 1), 'application/json'), { ic: 'download', kind: 'primary' }),
          btn('Open project', async () => { const [f] = await pickFiles({ accept: '.json' }); if (!f) return; try { importProject(JSON.parse(await f.text())); toast('Project opened.', 'ok'); } catch (e) { toast(`Could not open: ${e.message}`, 'bad'); } }, { ic: 'upload' }),
          btn('All results (CSV)', async () => { const rows = await allKpis(); if (!rows.length) { toast('There are no results to export yet.', 'info'); return; } downloadText('aerosuite-results.csv', toCsv(['Suite no.', 'Suite', 'Analysis', 'Key', 'Quantity', 'Value', 'Unit', 'Status', 'Run at'], rows), 'text/csv'); }, { ic: 'table' })),
        h('p', { class: 'muted small' }, `${Object.keys(state.runs).length} analyses stored for this case. The project file holds the case, every input override, published outputs and findings, so a colleague can reproduce the study exactly.`),
        h('div', { class: 'row' }, btn('Build printable report', buildReport, { ic: 'doc' }), btn('Clear stored results', async () => { if (confirm('Remove all stored results and findings for this case? Your case data and inputs are kept.')) { await clearRuns(); toast('Stored results cleared.', 'ok'); } }, { ic: 'close', kind: 'ghost' })))),
      card('Solver verification', h('div', { class: 'stack' }, h('p', { class: 'muted small' }, 'Runs every built-in benchmark in all 26 suites on this device and tabulates computed against exact values. Attach this to a report as code-verification evidence. Solution verification (mesh convergence) and validation against test data are done per analysis inside each suite.'), h('div', { class: 'row' }, btn('Verify all suites now', verifyAll, { ic: 'check', kind: 'primary' })), vProg))),
    h('div', { class: 'gap' }), vOut, reportHost);
  const off = on('case', () => {});
  return off;
}
