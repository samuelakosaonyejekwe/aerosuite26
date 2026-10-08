// Shared renderers for analysis results, recommendations and data tables.

import { h, icon, num, badge, card, toCsv, downloadText, btn } from './dom.js';
import { renderPlot } from './plots.js';

const STATUS_WORD = { ok: 'Within limits', warn: 'Check', bad: 'Not acceptable' };

export function kpiGrid(kpis) {
  return h('div', { class: 'kpis' }, kpis.map((k) => h('div', { class: `kpi ${k.status || ''}`, title: k.note || '' },
    h('span', { class: 'l' }, k.label),
    h('div', { class: 'v' }, num(k.value), k.unit && k.unit !== '-' ? h('small', null, k.unit) : null),
    k.status || k.note ? h('div', { class: 's' }, k.status ? h('b', null, STATUS_WORD[k.status]) : null, k.note ? h('span', null, (k.status ? '· ' : '') + k.note) : null) : null)));
}

const SEV = { critical: ['warn', 'Act now'], warn: ['warn', 'Attention'], advise: ['bulb', 'Opportunity'], info: ['info', 'Note'] };
export function recCard(r, origin) {
  const [ic, word] = SEV[r.severity] || SEV.info;
  return h('article', { class: `rec ${r.severity}` }, h('span', { class: 'ic' }, icon(ic, 22)),
    h('h4', null, r.title, ' ', badge(word, r.severity === 'critical' ? 'bad' : r.severity === 'warn' ? 'warn' : 'accent'), origin ? h('a', { class: 'small', href: origin.hash, style: { marginLeft: '8px', fontWeight: 500 } }, origin.label) : null),
    h('div', null, r.detail ? h('p', null, r.detail) : null, r.action ? h('p', { class: 'act' }, h('b', null, 'Recommended: '), r.action) : null, r.basis ? h('p', { class: 'basis' }, 'Basis: ', r.basis) : null));
}
export const recList = (recs, origin) => h('div', { class: 'stack' }, recs.map((r) => recCard(r, origin)));

export function dataTable(t, { max = 200 } = {}) {
  const isNum = (j) => t.rows.every((r) => typeof r[j] === 'number' || r[j] == null || r[j] === '');
  const numeric = t.columns.map((_, j) => isNum(j));
  return h('div', { class: 'stack' },
    h('div', { class: 'row' }, h('h4', { class: 'grow' }, t.title || ''), btn('CSV', () => downloadText(`${(t.title || 'table').toLowerCase().replace(/[^a-z0-9]+/g, '-')}.csv`, toCsv(t.columns, t.rows), 'text/csv'), { ic: 'download', kind: 'sm ghost', title: 'Download this table' })),
    h('div', { class: 'table-wrap' }, h('table', { class: 'data' },
      h('thead', null, h('tr', null, t.columns.map((c, j) => h('th', { class: numeric[j] ? 'num' : '' }, c)))),
      h('tbody', null, t.rows.slice(0, max).map((r) => h('tr', null, r.map((v, j) => h('td', { class: numeric[j] ? 'num' : '' }, typeof v === 'number' ? num(v, 5) : v == null ? '' : String(v)))))))),
    t.rows.length > max ? h('p', { class: 'muted small' }, `Showing ${max} of ${t.rows.length} rows; the CSV has all of them.`) : null);
}

/** Full result panel. Returns { el, destroy }. */
export function resultPanel(payload, { suiteLabel } = {}) {
  const { res, recs = [] } = payload, plots = [], host = h('div', { class: 'stack' });
  if (res.warnings.length) host.append(h('div', { class: 'note warn' }, icon('warn'), h('div', null, h('b', null, res.warnings.length === 1 ? 'Validity note' : 'Validity notes'), h('ul', { style: { margin: '4px 0 0', paddingLeft: '18px' } }, res.warnings.map((w) => h('li', null, w))))));
  host.append(kpiGrid(res.kpis));
  if (recs.length) host.append(card('What this means and what to do', recList(recs), { collapsible: true }));
  if (res.plots.length) { const grid = h('div', { class: 'plots' }); host.append(grid); for (const p of res.plots) { try { plots.push(renderPlot(grid, p)); } catch (e) { grid.append(h('div', { class: 'note bad' }, `Chart "${p.title}" could not be drawn: ${e.message}`)); } } }
  for (const t of res.tables) host.append(card(null, dataTable(t)));
  host.append(card('How this was calculated', h('div', { class: 'grid g3' },
    h('div', null, h('h4', null, 'Models used'), listOr(res.models, 'Not stated')),
    h('div', null, h('h4', null, 'Assumptions'), listOr(res.assumptions, 'None stated')),
    h('div', null, h('h4', null, 'Run record'), h('dl', { class: 'kv' }, h('dt', null, 'Solve time'), h('dd', null, `${res.elapsed_ms} ms`), h('dt', null, 'Inputs from case'), h('dd', null, `${(payload.linked || []).length} linked`), h('dt', null, 'Suite'), h('dd', null, suiteLabel || '')))),
  { collapsible: true, open: false }));
  return { el: host, destroy: () => plots.forEach((p) => p.destroy()) };
}
const listOr = (arr, none) => (arr?.length ? h('ul', { class: 'small', style: { margin: '6px 0 0', paddingLeft: '18px' } }, arr.map((x) => h('li', null, x))) : h('p', { class: 'muted small' }, none));
