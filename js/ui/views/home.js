// Overview: where you are, what to do next, and the 26 suites at a glance.

import { h, setKids, icon, num, btn, badge, card, ago } from '../dom.js';
import { SUITES, GROUPS } from '../../core/registry.js';
import { state, on } from '../../core/store.js';
import { derived } from '../../core/case.js';
import { loadSpec } from '../../core/spec.js';

export async function render(root, _p, { setCrumb }) {
  setCrumb('Overview');
  const host = h('div'); root.append(host);
  const paint = async () => {
    const c = state.case, d = derived(c), nRuns = Object.keys(state.runs).length, suitesRun = new Set(Object.keys(state.runs).map((k) => k.split('.')[0])).size;
    const recs = Object.values(state.recs).flat(), crit = recs.filter((r) => r.severity === 'critical').length, warn = recs.filter((r) => r.severity === 'warn').length;
    let spec = null; try { spec = await loadSpec(); } catch { /* offline first load */ }
    setKids(host, 
      h('section', { class: 'hero' },
        h('svg', { class: 'deco', viewBox: '0 0 24 24', fill: 'none', stroke: '#fff', 'stroke-width': 0.6 }, h('path', { d: 'M3 13l7-1 5-7h2l-2 7 5-.5 1.5-2H23l-1 3 1 3h-1.5L20 13.5l-5-.5 2 7h-2l-5-7-7-1z' })),
        h('h1', null, 'One connected workspace for every aircraft engineering analysis.'),
        h('p', null, 'Describe an aeroplane, helicopter, rotorcraft or UAV once. Twenty-six linked suites — from aerodynamics and structures to mission, safety and economics — share that single case, pass results to each other, and turn the numbers into recommendations. It runs entirely on this device, online or offline.'),
        h('div', { class: 'row' }, btn('Set up your aircraft', () => (location.hash = '#/case'), { ic: 'sliders', kind: 'primary big' }), btn('Run all 26 suites', () => (location.hash = '#/integrated'), { ic: 'graph', kind: 'ghost big' }))),
      h('div', { class: 'gap' }),
      h('div', { class: 'steps' },
        h('a', { class: `step ${state.preset ? 'done' : ''}`, href: '#/case' }, h('b', null, 'Define the case'), h('p', null, 'Pick a starting aircraft, edit any value, import a data file, and choose where on Earth it operates.')),
        h('a', { class: `step ${state.geometry.length ? 'done' : ''}`, href: '#/geometry' }, h('b', null, 'Bring geometry or a mesh'), h('p', null, 'Import CAD, tessellated or mesh files from other tools. Inspect, heal, measure and section them.')),
        h('a', { class: `step ${nRuns ? 'done' : ''}`, href: nRuns ? '#/integrated' : '#/suite/cfd' }, h('b', null, 'Run analyses'), h('p', null, 'Work suite by suite, or run the whole coupled chain. Check mesh convergence and uncertainty as you go.')),
        h('a', { class: `step ${recs.length ? 'done' : ''}`, href: '#/decisions' }, h('b', null, 'Decide'), h('p', null, 'Read prioritised recommendations with their engineering basis, sustainability and cost effects.'))),
      h('div', { class: 'gap' }),
      h('div', { class: 'grid g3' },
        card('Current case', h('div', { class: 'stack' }, h('div', null, h('b', null, c.meta.name || 'Unnamed case'), ' ', badge(c.meta.type, 'accent')),
          h('dl', { class: 'kv' }, h('dt', null, 'Take-off mass'), h('dd', null, `${num(c.mass.mtow_kg)} kg`),
            c.wing.S_m2 > 0 ? [h('dt', null, 'Wing'), h('dd', null, `${num(c.wing.S_m2)} m², span ${num(c.wing.b_m)} m, aspect ratio ${num(d.AR, 3)}`)] : null,
            c.rotor.R_m > 0 ? [h('dt', null, 'Rotor'), h('dd', null, `${c.rotor.n_blades} blades, radius ${num(c.rotor.R_m)} m, ${num(c.rotor.rpm)} rpm`)] : null,
            h('dt', null, 'Powerplant'), h('dd', null, `${c.prop.n_eng} × ${c.prop.type}`), h('dt', null, 'Mission'), h('dd', null, `${num(c.mission.range_km)} km at ${num(c.mission.cruise_alt_m)} m`)),
          h('div', null, btn('Edit case', () => (location.hash = '#/case'), { ic: 'sliders', kind: 'sm' })))),
        card('Operating site and live conditions', h('div', { class: 'stack' }, h('div', null, icon('pin', 16), ' ', h('b', null, c.site.name)),
          h('dl', { class: 'kv' }, h('dt', null, 'Elevation'), h('dd', null, `${num(c.site.elev_m)} m`), h('dt', null, 'Temperature'), h('dd', null, `${num(c.site.T_C, 3)} °C`), h('dt', null, 'Pressure'), h('dd', null, `${num(c.site.p_hPa, 5)} hPa`), h('dt', null, 'Wind'), h('dd', null, `${num(c.site.wind_ms, 3)} m/s from ${num(c.site.wind_dir_deg, 3)}°`), h('dt', null, 'Runway'), h('dd', null, `${num(c.site.runway_len_m)} m`), h('dt', null, 'Source'), h('dd', null, c.site.source, c.site.updated ? ` · ${ago(c.site.updated)}` : '')),
          h('div', null, btn(c.site.lat == null ? 'Choose a location' : 'Change location', () => (location.hash = '#/case/site'), { ic: 'globe', kind: 'sm' })))),
        card('Progress', h('div', { class: 'stack' },
          h('dl', { class: 'kv' }, h('dt', null, 'Suites with results'), h('dd', null, `${suitesRun} of 26`), h('dt', null, 'Analyses run'), h('dd', null, String(nRuns)), h('dt', null, 'Findings'), h('dd', null, crit ? badge(`${crit} critical`, 'bad') : null, ' ', warn ? badge(`${warn} attention`, 'warn') : null, !crit && !warn ? (recs.length ? badge(`${recs.length} notes`, 'ok') : 'none yet') : null)),
          h('progress', { max: 26, value: suitesRun }),
          h('div', { class: 'row' }, btn('Decision support', () => (location.hash = '#/decisions'), { ic: 'bulb', kind: 'sm' }), btn('Report', () => (location.hash = '#/reports'), { ic: 'doc', kind: 'sm ghost' }))))),
      ...GROUPS.map((g) => h('section', null, h('h2', { style: { margin: '26px 0 12px' } }, g), h('div', { class: 'grid g3' }, SUITES.filter((s) => s.group === g).map((s) => {
        const runs = Object.entries(state.runs).filter(([k]) => k.startsWith(s.id + '.')), st = runs.some(([, v]) => v.status === 'bad') ? 'bad' : runs.some(([, v]) => v.status === 'warn') ? 'warn' : runs.length ? 'ok' : '';
        return h('a', { class: 'suite-card', href: `#/suite/${s.id}` }, h('div', { class: 'top-r' }, h('span', { class: 'num' }, icon(s.icon, 20)), h('h3', null, s.short), h('span', { class: 'nn' }, `${s.n}/26`)), h('p', null, spec?.suites[s.n]?.scope || s.title), h('div', { class: 'row' }, st ? badge(st === 'ok' ? `${runs.length} run${runs.length > 1 ? 's' : ''} · within limits` : st === 'warn' ? 'Results need a look' : 'Limit exceeded', st) : badge('Not run yet')));
      })))),
      h('div', { class: 'gap' }),
      h('div', { class: 'note' }, icon('shield'), h('div', null, h('b', null, 'Engineering credibility. '), 'Every result states the model, its assumptions and its validity limits, and each suite carries verification benchmarks you can run yourself. Results distinguish calculated values, empirical estimates and assumed inputs. A passing calculation is not airworthiness evidence: validate against test data for your aircraft before any safety-critical or certification use. ', h('a', { href: '#/about' }, 'Capabilities and limits'))));
  };
  await paint();
  const off = on('run', paint), off2 = on('case', paint);
  return () => { off(); off2(); };
}
