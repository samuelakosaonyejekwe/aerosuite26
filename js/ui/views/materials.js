// Materials library: browse, search and compare every material the suites use; see where the numbers
// come from; add your own from test or supplier data; import and export.

import { h, icon, clear, add, num, btn, badge, card, empty, toast, pickFiles, downloadText, toCsv } from '../dom.js';
import { renderPlot } from '../plots.js';
import { state, patchCase, on } from '../../core/store.js';
import { TABLES, CLASS_LABEL, FIELDS, INDICES, isBuiltin, customMaterials, registerCustom, validate, saveCustom, deleteCustom, toCsvRows, parseMaterials } from '../../core/matlib.js';
import { designAllowables } from '../../data/materials.js';

const CASE_SLOT = { metals: ['struct', 'material', 'primary structural metal'], plies: ['struct', 'ply', 'composite ply system'], batteries: ['systems', 'batt_chem', 'battery chemistry'], fluids: ['prop', 'fuel', 'fuel'] };
const SUMMARY = { metals: ['E', 'rho', 'Sy', 'Su', 'KIc'], plies: ['E1', 'E2', 'Xt', 'Xc', 'rho'], fluids: ['rho', 'mu', 'cp', 'LHV'], batteries: ['wh_kg', 'v_nom', 'r_mohm_ah', 'cycles_80'] };
const show = (f, v) => (v == null || v === '' ? '–' : num(v / f.scale));
let extraP = null;
/** Optional sourced reference entries (js/data/materials-extra.json); absent in some builds. */
const loadExtra = () => (extraP ||= (async () => { try { const r = await fetch(new URL('../../data/materials-extra.json', import.meta.url)); return r.ok ? await r.json() : null; } catch { return null; } })());

export async function render(root, [clsArg, nameArg], { setCrumb }) {
  setCrumb('Materials library');
  registerCustom(customMaterials());
  const host = h('div'); root.append(host);
  let cls = TABLES[clsArg] ? clsArg : 'metals', selected = nameArg && TABLES[cls][nameArg] ? nameArg : null, query = '', sortKey = null, sortDir = 1, editing = null;
  const plots = [], killPlots = () => plots.splice(0).forEach((p) => p.destroy());
  const tabs = h('div', { class: 'tabs', role: 'tablist' }), body = h('div', { class: 'stack' });
  const extra = await loadExtra();

  host.append(
    h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Shared by all 26 suites'), h('h1', null, 'Materials library'), h('p', null, 'The metals, composite ply systems, fuels, fluids and battery cells the analyses draw on — with every property, the comparison indices engineers actually use, and room for your own materials from coupon tests or supplier data.')),
      h('div', { class: 'row' }, btn('Add a material', () => { editing = { name: '', base: selected || Object.keys(TABLES[cls])[0], fresh: true }; paint(); }, { ic: 'upload', kind: 'primary' }))),
    tabs, body);

  const names = () => Object.keys(TABLES[cls]).filter((n) => !n.startsWith('__') && n.toLowerCase().includes(query.toLowerCase()));
  function paintTabs() { clear(tabs); add(tabs, Object.keys(TABLES).map((k) => h('button', { class: `tab ${k === cls ? 'on' : ''}`, role: 'tab', 'aria-selected': k === cls, onclick: () => { cls = k; selected = null; editing = null; sortKey = null; query = ''; paintTabs(); paint(); } }, CLASS_LABEL[k], ' ', badge(String(Object.keys(TABLES[k]).length))))); }

  function paint() {
    killPlots(); clear(body);
    const T = TABLES[cls], F = FIELDS[cls], cols = SUMMARY[cls].map((k) => F.find((f) => f.key === k)), slot = CASE_SLOT[cls], inCase = state.case[slot[0]][slot[1]];
    let list = names();
    if (sortKey) list.sort((a, b) => (sortKey === 'name' ? a.localeCompare(b) : ((T[a][sortKey] ?? -Infinity) - (T[b][sortKey] ?? -Infinity))) * sortDir);
    const th = (label, key, numeric) => h('th', { class: numeric ? 'num' : '', style: { cursor: 'pointer' }, title: 'Sort', onclick: () => { sortDir = sortKey === key ? -sortDir : 1; sortKey = key; paint(); } }, label, sortKey === key ? (sortDir > 0 ? ' ▲' : ' ▼') : '');
    const search = h('input', { class: 'inp', type: 'search', placeholder: `Search ${CLASS_LABEL[cls].toLowerCase()}…`, value: query, 'aria-label': 'Search materials', oninput: (e) => { query = e.target.value; const pos = e.target.selectionStart; paint(); const s2 = body.querySelector('input[type=search]'); s2?.focus(); s2?.setSelectionRange(pos, pos); } });

    body.append(card(null, h('div', { class: 'stack' },
      h('div', { class: 'row' }, h('div', { class: 'grow', style: { minWidth: '200px' } }, search),
        btn('Export', () => { const r = toCsvRows(cls); downloadText(`aerosuite-${cls}.csv`, toCsv(r.columns, r.rows), 'text/csv'); }, { ic: 'download', kind: 'sm', title: 'Download this family as a table (SI units)' }),
        btn('Import', importFile, { ic: 'upload', kind: 'sm', title: 'Add materials from a CSV or JSON file' }),
        btn('Template', () => { const r = toCsvRows(cls, [Object.keys(T)[0]]); downloadText(`aerosuite-${cls}-template.csv`, toCsv(r.columns, r.rows), 'text/csv'); }, { ic: 'doc', kind: 'sm ghost', title: 'A one-row example showing the column names and SI units' })),
      list.length ? h('div', { class: 'table-wrap' }, h('table', { class: 'data' },
        h('thead', null, h('tr', null, th('Material', 'name', false), ...cols.map((f) => th(`${f.label.replace(/ \(typical\)/, '')} [${f.unit}]`, f.key, true)), h('th', null, 'Status'))),
        h('tbody', null, list.map((n) => h('tr', { style: { cursor: 'pointer', background: n === selected ? 'var(--accent-soft)' : '' }, tabIndex: 0, onclick: () => { selected = n; editing = null; paint(); }, onkeydown: (e) => { if (e.key === 'Enter') { selected = n; editing = null; paint(); } } },
          h('td', null, h('b', null, n)), ...cols.map((f) => h('td', { class: 'num' }, show(f, T[n][f.key]))),
          h('td', null, isBuiltin(cls, n) ? badge('Built-in') : badge('Yours', 'accent'), ' ', cls === 'metals' && designAllowables(T[n]).design ? badge('Design allowable', 'ok') : null, ' ', n === inCase ? badge('In your aircraft', 'ok') : null))))))
        : empty('search', 'Nothing matches', 'Try a shorter search, or add the material you need.'))));

    if (editing) body.append(editor());
    else if (selected && T[selected]) body.append(detail(selected));
    else body.append(h('div', { class: 'note' }, icon('info'), h('div', null, 'Tap a material to see all its properties. ', h('b', null, `${slot[2][0].toUpperCase()}${slot[2].slice(1)} currently in your aircraft: `), inCase, '.')));

    // comparison chart
    if (list.length > 1) {
      const opts = [...F.filter((f) => list.some((n) => Number.isFinite(T[n][f.key]))).map((f) => ({ id: f.key, label: `${f.label} [${f.unit}]`, fn: (m) => (m[f.key] == null ? NaN : m[f.key] / f.scale) })), ...INDICES[cls].map(([label, unit, fn], i) => ({ id: 'idx' + i, label: `${label} [${unit}]`, fn }))];
      const sel = h('select', { class: 'inp', style: { width: 'auto', maxWidth: '100%' }, 'aria-label': 'Property to compare' }, opts.map((o) => h('option', { value: o.id, selected: o.id === 'idx1' }, o.label))), chart = h('div');
      const draw = () => { killPlots(); clear(chart); const o = opts.find((x) => x.id === sel.value) || opts[0], rows = list.map((n) => [n, o.fn(T[n])]).filter((r) => Number.isFinite(r[1])).sort((a, b) => b[1] - a[1]); if (rows.length) plots.push(renderPlot(chart, { type: 'bar', title: o.label, ylabel: o.label, categories: rows.map((r) => r[0]), series: [{ name: o.label, y: rows.map((r) => r[1]) }] })); };
      sel.addEventListener('change', draw);
      body.append(card('Compare', h('div', { class: 'stack' }, h('div', { class: 'row' }, h('span', { class: 'muted small' }, 'Rank the listed materials by'), sel), chart), { collapsible: true }));
      draw();
    }
    if (extra?.items?.length) body.append(reference());
    body.append(h('div', { class: 'note' }, icon('shield'), h('div', null, h('b', null, 'About the numbers. '), 'Built-in values are typical room-temperature figures for preliminary work. Where a statistically based design allowable has been verified against a handbook it is held separately and used for margins of safety; the material’s page says which basis applies. Fatigue, fracture and thermal constants are typical values unless a source is shown. For substantiation, enter your programme’s own allowables as a material of your own.')));
  }

  // ---- detail ---------------------------------------------------------------------------------
  function detail(n) {
    const T = TABLES[cls], m = T[n], F = FIELDS[cls], groups = [...new Set(F.map((f) => f.group))], slot = CASE_SLOT[cls], mine = !isBuiltin(cls, n);
    const da = cls === 'metals' ? designAllowables(m) : null;
    const canUse = cls !== 'fluids' || m.LHV > 0;
    return card(n, h('div', { class: 'stack' },
      h('div', { class: 'row' }, mine ? badge('Your material', 'accent') : badge('Built-in'), da ? badge(da.design ? 'Margins use the design allowable' : 'Margins use typical strengths', da.design ? 'ok' : 'warn') : null,
        h('span', { class: 'grow' }),
        canUse ? btn('Use in my aircraft', () => { patchCase(slot[0], slot[1], n); toast(`${n} is now the ${slot[2]} of your aircraft. Re-run the suites to see its effect.`, 'ok'); paint(); }, { ic: 'check', kind: 'primary sm', disabled: state.case[slot[0]][slot[1]] === n }) : null,
        btn('Copy and edit', () => { editing = { name: `${n} (my data)`, base: n, fresh: true }; paint(); }, { ic: 'sliders', kind: 'sm', title: 'Start a material of your own from these values' }),
        mine ? btn('Edit', () => { editing = { name: n, base: n, fresh: false }; paint(); }, { kind: 'sm' }) : null,
        mine ? btn('Delete', () => { if (confirm(`Delete “${n}” from your materials?`)) { if (state.case[slot[0]][slot[1]] === n) patchCase(slot[0], slot[1], Object.keys(T).find((k) => isBuiltin(cls, k))); deleteCustom(cls, n); selected = null; toast('Material deleted.', 'ok'); paintTabs(); paint(); } }, { kind: 'sm ghost' }) : null),
      da ? h('p', { class: 'small muted' }, h('b', null, 'Strength basis: '), da.basis) : null, m.note ? h('p', { class: 'small muted' }, h('b', null, 'Note: '), m.note) : null,
      h('div', { class: 'grid g3' }, groups.map((g) => h('div', null, h('h4', null, g), h('dl', { class: 'kv small', style: { marginTop: '6px' } }, F.filter((f) => f.group === g).map((f) => [h('dt', { title: f.help || '' }, f.label), h('dd', null, m[f.key] == null ? h('span', { class: 'muted' }, 'not given') : `${show(f, m[f.key])} ${f.unit === '-' ? '' : f.unit}`)])))),
        h('div', null, h('h4', null, 'Comparison indices'), h('dl', { class: 'kv small', style: { marginTop: '6px' } }, INDICES[cls].map(([label, unit, fn]) => { const v = fn(m); return [h('dt', null, label), h('dd', null, Number.isFinite(v) ? `${num(v)} ${unit === '-' ? '' : unit}` : h('span', { class: 'muted' }, 'n/a'))]; }))),
        cls === 'metals' && m.JC ? h('div', null, h('h4', null, 'High-rate plasticity (Johnson–Cook)'), h('dl', { class: 'kv small', style: { marginTop: '6px' } }, [['A [MPa]', m.JC.A / 1e6], ['B [MPa]', m.JC.B / 1e6], ['n', m.JC.n], ['C', m.JC.C], ['m', m.JC.m], ['Melting point [K]', m.JC.Tm]].map(([k, v]) => [h('dt', null, k), h('dd', null, num(v))]))) : null)));
  }

  // ---- editor ---------------------------------------------------------------------------------
  function editor() {
    const F = FIELDS[cls], base = TABLES[cls][editing.base] || {}, groups = [...new Set(F.map((f) => f.group))], inputs = {}, msg = h('div', { class: 'stack' });
    const nameIn = h('input', { class: 'inp', type: 'text', value: editing.name, placeholder: 'e.g. Al 2024-T3 sheet, my A-basis', 'aria-label': 'Material name', maxLength: 80 });
    const noteIn = h('input', { class: 'inp', type: 'text', value: (editing.fresh ? '' : base.note) || '', placeholder: 'Where the values come from: test report, supplier datasheet, specification…', 'aria-label': 'Source note', maxLength: 300 });
    const read = () => { const m = {}; for (const f of F) { const raw = inputs[f.key].value.trim(); if (raw !== '') m[f.key] = Number(raw) * f.scale; } if (noteIn.value.trim()) m.note = noteIn.value.trim(); if (cls === 'metals' && base.JC) m.JC = base.JC; return m; };
    const save = () => {
      const name = nameIn.value.trim(), m = read(), v = validate(cls, name, m, { editing: editing.fresh ? null : editing.base });
      clear(msg); add(msg, [...v.errors.map((e) => h('div', { class: 'note bad' }, icon('warn'), h('div', null, e))), ...v.warnings.map((w) => h('div', { class: 'note warn' }, icon('info'), h('div', null, w)))]);
      if (v.errors.length) return;
      saveCustom(cls, name, m, editing.fresh ? null : editing.base); selected = name; editing = null; toast(`“${name}” saved. It now appears in every suite’s material list.`, 'ok'); paintTabs(); paint();
    };
    return card(editing.fresh ? 'New material' : `Edit ${editing.base}`, h('div', { class: 'stack' },
      h('p', { class: 'muted small' }, editing.fresh ? `Values start from “${editing.base}”. Replace what you have data for; anything you leave is carried over from that material, so say so in the source note.` : 'Change any value and save.'),
      h('div', { class: 'grid g2' }, h('label', { class: 'stack small', style: { gap: '4px' } }, h('span', { class: 'muted' }, 'Name'), nameIn), h('label', { class: 'stack small', style: { gap: '4px' } }, h('span', { class: 'muted' }, 'Source note'), noteIn)),
      h('div', { class: 'grid g3' }, groups.map((g) => h('div', { class: 'form' }, h('h4', null, g), F.filter((f) => f.group === g).map((f) => { const id = `m-${f.key}`; inputs[f.key] = h('input', { class: 'inp', id, type: 'number', step: 'any', inputMode: 'decimal', value: base[f.key] == null ? '' : Number((base[f.key] / f.scale).toPrecision(7)), placeholder: f.req ? '' : 'optional' }); return h('div', { class: 'field' }, h('label', { for: id }, f.label, f.unit !== '-' ? h('span', { class: 'u' }, `[${f.unit}]`) : null), inputs[f.key], f.help ? h('div', { class: 'help' }, f.help) : null); })))),
      msg, h('div', { class: 'row' }, btn('Save material', save, { ic: 'check', kind: 'primary' }), btn('Cancel', () => { editing = null; paint(); }, { kind: 'ghost' }))));
  }

  // ---- import ---------------------------------------------------------------------------------
  async function importFile() {
    const [f] = await pickFiles({ accept: '.csv,.json,.txt,.tsv' }); if (!f) return;
    try {
      if (f.size > 5e6) throw new Error('The file is larger than 5 MB; a materials table is normally a few kilobytes.');
      const items = parseMaterials(await f.text()); let ok = 0; const problems = [];
      for (const it of items) {
        const complete = { ...(TABLES[it.cls][it.name] && !isBuiltin(it.cls, it.name) ? TABLES[it.cls][it.name] : {}), ...it.data }, v = validate(it.cls, it.name, complete, { editing: it.name });
        if (v.errors.length) problems.push(`${it.name || '(no name)'}: ${v.errors[0]}`); else { saveCustom(it.cls, it.name, complete); ok++; }
      }
      toast(`${ok} material${ok === 1 ? '' : 's'} imported${problems.length ? `; ${problems.length} skipped` : ''}.`, problems.length && !ok ? 'bad' : 'ok', 6000);
      paintTabs(); paint();
      if (problems.length) body.prepend(h('div', { class: 'note warn' }, icon('warn'), h('div', null, h('b', null, 'Not imported: '), h('ul', { style: { margin: '4px 0 0', paddingLeft: '18px' } }, problems.slice(0, 12).map((p) => h('li', null, p))), 'Each material needs every required property of its family, in SI units. Download the template to see the columns.')));
    } catch (e) { toast(`${f.name} could not be read: ${e.message}`, 'bad', 7000); }
  }

  // ---- sourced reference data (optional file) ----------------------------------------------------
  function reference() {
    const items = extra.items.filter((x) => (x.class || 'metals') === cls || (cls === 'metals' && !['composite ply', 'ply'].includes(x.category) && !x.class) || (cls === 'plies' && ['composite ply', 'ply'].includes(x.category)));
    if (!items.length) return null;
    const F = FIELDS[cls], keys = SUMMARY[cls].filter((k) => items.some((x) => Number.isFinite(x[k])));
    return card(`Sourced reference data (${items.length})`, h('div', { class: 'stack' },
      h('p', { class: 'muted small' }, 'Additional materials with values read from the cited public documents. They are reference entries: to use one in the analyses, press “Add to my materials” and complete any property the source does not give.'),
      h('div', { class: 'table-wrap' }, h('table', { class: 'data' }, h('thead', null, h('tr', null, h('th', null, 'Material'), h('th', null, 'Form / condition'), ...keys.map((k) => { const f = F.find((x) => x.key === k); return h('th', { class: 'num' }, `${f.label.replace(/ \(typical\)/, '')} [${f.unit}]`); }), h('th', null, 'Source'), h('th', null, ''))),
        h('tbody', null, items.map((x) => h('tr', null, h('td', null, h('b', null, x.name)), h('td', null, x.form || x.condition || ''), ...keys.map((k) => { const f = F.find((y) => y.key === k); return h('td', { class: 'num' }, show(f, x[k])); }),
          h('td', null, x.source?.url ? h('a', { href: x.source.url, target: '_blank', rel: 'noopener noreferrer' }, x.source.org || x.source.title || 'source') : (x.source?.title || '')),
          h('td', null, btn('Add to my materials', () => { const nearest = Object.keys(TABLES[cls]).filter((n) => isBuiltin(cls, n)).sort((a, b) => Math.abs(Math.log((TABLES[cls][a].rho || 1) / (x.rho || 1))) - Math.abs(Math.log((TABLES[cls][b].rho || 1) / (x.rho || 1))))[0]; const merged = { ...TABLES[cls][nearest] }; for (const f of F) if (Number.isFinite(x[f.key])) merged[f.key] = x[f.key]; delete merged.Sy_A; delete merged.Su_A; delete merged.allow_ref; merged.note = `Sourced values: ${x.source?.title || ''}. Other properties carried over from ${nearest} — replace with data for this material.`.slice(0, 300); TABLES[cls]['__draft__'] = merged; editing = { name: x.name, base: '__draft__', fresh: true }; paint(); delete TABLES[cls]['__draft__']; body.querySelector('.card input[type=text]')?.scrollIntoView({ block: 'center' }); }, { kind: 'sm' })))))))), { collapsible: true, open: false });
  }

  paintTabs(); paint();
  const off = on('materials', (m) => { registerCustom(m); paintTabs(); paint(); });
  return () => { off(); killPlots(); };
}
