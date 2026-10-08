// Integrated run: all 26 suites as one coupled system. Shows who feeds whom, runs the chain in
// dependency order with each suite's outputs handed to the next, and summarises the outcome.

import { h, setKids, icon, clear, num, btn, badge, card, toast } from '../dom.js';
import { SUITES, suiteMeta } from '../../core/registry.js';
import { state, saveRun, on } from '../../core/store.js';
import { runJob, cancelAll } from '../../core/runner.js';
import { dataTable } from '../results.js';

let graphCache = null;
async function couplingGraph() {
  if (graphCache) return graphCache;
  const nodes = [];
  for (const s of SUITES) { try { const d = await runJob({ kind: 'describe', suite: s.id, case: state.case, up: state.up }); nodes.push({ id: s.id, consumes: d.consumes, provides: d.provides, n: d.analyses.length, applicable: d.analyses.filter((a) => a.applicable === true).length }); } catch { nodes.push({ id: s.id, consumes: [], provides: [], n: 0, applicable: 0, missing: true }); } }
  return (graphCache = nodes);
}

export async function render(root, _p, { setCrumb }) {
  setCrumb('Integrated run');
  const host = h('div'); root.append(host);
  const mapHost = h('div'), detail = h('div', { class: 'muted small', style: { minHeight: '3em' } }, 'Hover or tap a suite to see what it receives and what it passes on.');
  const prog = h('progress', { max: 1, value: 0, hidden: true }), msg = h('span', { class: 'muted small' }), out = h('div', { class: 'stack' });
  const runBtn = btn('Run all suites', () => go(), { ic: 'play', kind: 'primary big' }), stopBtn = btn('Stop', () => cancelAll(), { ic: 'stop', kind: 'ghost' }); stopBtn.hidden = true;
  const passes = h('select', { class: 'inp', style: { width: 'auto' }, 'aria-label': 'Coupling passes' }, h('option', { value: 1 }, 'One pass'), h('option', { value: 2, selected: true }, 'Two passes (feeds results back upstream)'), h('option', { value: 3 }, 'Three passes'));
  host.append(
    h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Coupled multidisciplinary analysis'), h('h1', null, 'Integrated run'), h('p', null, 'The suites are connected, not merged: each keeps its own models, but they exchange loads, masses, performance, limits and costs through a shared data bus. Running the chain propagates one change — a heavier wing, a hotter day, a dearer fuel — through every discipline to the bottom line.'))),
    h('div', { class: 'grid g2' },
      card('How the 26 suites feed each other', h('div', { class: 'stack' }, mapHost, detail)),
      card('Run the chain', h('div', { class: 'stack' },
        h('p', { class: 'muted small' }, 'Upstream suites run first (aerodynamics and propulsion before performance, performance before mission, mission before economics). A second pass lets downstream results — optimised mass, fatigue life, availability — flow back into the suites that use them.'),
        h('div', { class: 'row' }, runBtn, stopBtn, passes), prog, msg,
        h('div', { class: 'note' }, icon('info'), h('div', null, 'Each suite uses your saved input overrides. Analyses that do not apply to the current aircraft class are skipped and listed.'))))),
    h('div', { class: 'gap' }), out);

  const nodes = await couplingGraph();
  drawMap();
  function drawMap() {
    const W = 560, H = 560, cx = W / 2, cy = H / 2, R = 215, pos = {};
    SUITES.forEach((s, i) => { const a = (i / SUITES.length) * 2 * Math.PI - Math.PI / 2; pos[s.id] = [cx + R * Math.cos(a), cy + R * Math.sin(a), a]; });
    const svg = h('svg', { class: 'flowmap', viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'Data flow between suites' }), links = [];
    for (const n of nodes) for (const c of n.consumes) if (pos[c.from] && c.from !== n.id) { const [x1, y1] = pos[c.from], [x2, y2] = pos[n.id]; const p = h('path', { class: 'lnk', d: `M${x1},${y1} Q${cx},${cy} ${x2},${y2}` }); links.push({ p, from: c.from, to: n.id, c }); svg.append(p); }
    const focus = (id) => {
      for (const l of links) l.p.classList.toggle('hot', id != null && (l.from === id || l.to === id));
      svg.querySelectorAll('.nd').forEach((g) => g.classList.toggle('hot', g.dataset.id === id));
      if (!id) return;
      const ins = links.filter((l) => l.to === id), outs = links.filter((l) => l.from === id), m = suiteMeta(id);
      setKids(detail, h('b', { style: { color: 'var(--ink)' } }, `${m.d}. ${m.short} `), h('a', { href: `#/suite/${id}` }, 'open'), h('br'),
        'Receives: ', ins.length ? ins.map((l) => `${suiteMeta(l.from).short} (${l.c.keys.slice(0, 3).join(', ')}${l.c.keys.length > 3 ? '…' : ''})`).join('; ') : 'the shared case only', h('br'), 'Feeds: ', outs.length ? outs.map((l) => suiteMeta(l.to).short).join(', ') : 'reports and decision support');
    };
    for (const s of SUITES) {
      const [x, y, a] = pos[s.id], runs = Object.entries(state.runs).filter(([k]) => k.startsWith(s.id + '.')).map(([, v]) => v.status), st = runs.includes('bad') ? 'bad' : runs.includes('warn') ? 'warn' : runs.length ? 'ok' : '';
      const g = h('g', { class: `nd ${st}`, dataset: { id: s.id }, tabindex: 0, onmouseenter: () => focus(s.id), onfocus: () => focus(s.id), onclick: () => focus(s.id), ondblclick: () => (location.hash = `#/suite/${s.id}`) }, h('circle', { cx: x, cy: y, r: 13 }));
      const t = h('text', { x, y: y + 4, 'text-anchor': 'middle', style: { fontWeight: 650, pointerEvents: 'none' } }); t.textContent = s.d; g.append(t);
      const lx = cx + (R + 22) * Math.cos(a), ly = cy + (R + 22) * Math.sin(a), right = Math.cos(a) > 0.15, left = Math.cos(a) < -0.15;
      const lab = h('text', { x: lx, y: ly + 4, 'text-anchor': right ? 'start' : left ? 'end' : 'middle', style: { fontSize: '9.5px', fill: 'var(--ink-2)' } }); lab.textContent = s.short.length > 14 ? s.short.slice(0, 13) + '…' : s.short; g.append(lab);
      svg.append(g);
    }
    setKids(mapHost, svg);
  }

  async function go() {
    runBtn.disabled = true; stopBtn.hidden = false; prog.hidden = false; prog.value = 0; clear(out);
    const t0 = Date.now(), nPass = Number(passes.value); let up = state.up, last = null;
    try {
      for (let p = 0; p < nPass; p++) {
        last = await runJob({ kind: 'integrated', case: state.case, up, overridesAll: state.overrides }, (f, m) => { prog.value = (p + f) / nPass; msg.textContent = `Pass ${p + 1} of ${nPass} · ${m || ''}`; });
        for (const e of last.log) if (e.payload) await saveRun(e.suite, e.analysis, e.payload);
        up = state.up;
      }
      const ok = last.log.filter((e) => e.payload), skipped = last.log.filter((e) => e.skipped), failed = last.log.filter((e) => e.error), recs = ok.flatMap((e) => e.payload.recs), crit = recs.filter((r) => r.severity === 'critical').length;
      msg.textContent = '';
      toast(`Integrated run finished: ${ok.length} analyses in ${((Date.now() - t0) / 1000).toFixed(1)} s.`, crit ? 'bad' : 'ok');
      out.append(
        h('div', { class: 'kpis' }, kp('Analyses solved', ok.length), kp('Skipped (not applicable)', skipped.length), kp('Failed', failed.length, failed.length ? 'bad' : 'ok'), kp('Critical findings', crit, crit ? 'bad' : 'ok'), kp('Findings needing attention', recs.filter((r) => r.severity === 'warn').length), kp('Wall time', (Date.now() - t0) / 1000, '', 's')),
        h('div', { class: 'row' }, btn('See the recommendations', () => (location.hash = '#/decisions'), { ic: 'bulb', kind: 'primary' }), btn('Build the report', () => (location.hash = '#/reports'), { ic: 'doc' })),
        card(null, dataTable({ title: 'Execution order and outcome', columns: ['#', 'Suite', 'Analysis', 'Outcome', 'Time [ms]', 'Warnings', 'Findings'], rows: last.log.map((e) => [suiteMeta(e.suite).d, suiteMeta(e.suite).short, e.title, e.error ? `Failed: ${e.error}` : e.skipped ? 'Skipped — not applicable' : e.payload.res.kpis.some((k) => k.status === 'bad') ? 'Limit exceeded' : e.payload.res.kpis.some((k) => k.status === 'warn') ? 'Check' : 'OK', e.payload?.res.elapsed_ms ?? null, e.payload?.res.warnings.length ?? null, e.payload?.recs.length ?? null]) }, { max: 400 })),
        card('Headline outputs on the data bus', busTable(), { collapsible: true }));
      drawMap();
    } catch (e) { msg.textContent = ''; if (e.message !== 'Cancelled') out.append(h('div', { class: 'note bad' }, icon('warn'), h('div', null, e.message))); else toast('Run stopped.', 'info'); }
    finally { runBtn.disabled = false; stopBtn.hidden = true; prog.hidden = true; }
  }
  const kp = (label, value, st = '', unit = '') => h('div', { class: `kpi ${st}` }, h('span', { class: 'l' }, label), h('div', { class: 'v' }, num(value), unit ? h('small', null, unit) : null));
  function busTable() {
    const rows = [];
    for (const n of nodes) for (const p of n.provides) { const v = state.up[n.id]?.[p.key]; rows.push([`${suiteMeta(n.id).d}. ${suiteMeta(n.id).short}`, p.label, typeof v === 'number' ? v : null, p.unit || '', typeof v === 'number' ? 'published' : 'not available for this aircraft']); }
    return dataTable({ title: '', columns: ['Suite', 'Quantity', 'Value', 'Unit', 'Status'], rows }, { max: 400 });
  }
  if (Object.keys(state.runs).length) out.append(card('Headline outputs on the data bus', busTable(), { collapsible: true, open: false }));
  const off = on('case', () => { graphCache = null; });
  return off;
}
