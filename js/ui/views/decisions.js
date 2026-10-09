// Decision support: every recommendation from every suite, prioritised, with a cross-discipline
// scorecard and a sustainability view built from the coupled results.

import { h, setKids, icon, clear, add, num, btn, badge, card, empty, downloadText } from '../dom.js';
import { SUITES, suiteMeta } from '../../core/registry.js';
import { state, on } from '../../core/store.js';
import { recCard } from '../results.js';
import { renderPlot } from '../plots.js';
import { FLUIDS } from '../../data/materials.js';
import { status as liveStatus } from '../../core/live.js';
import { idb } from '../../core/store.js';

const ORDER = { critical: 0, warn: 1, advise: 2, info: 3 };
const SCORE = [ // [suite, key, label, unit, good direction]
  ['cfd', 'LD_max', 'Maximum lift-to-drag ratio', '-'], ['performance', 'tofl_m', 'Take-off distance', 'm'], ['performance', 'range_km', 'Range at max payload', 'km'], ['performance', 'ceiling_m', 'Service ceiling', 'm'],
  ['rotorcraft', 'FM', 'Rotor figure of merit', '-'], ['rotorcraft', 'hover_power_W', 'Hover power', 'W'], ['propulsion', 'tsfc_kg_Ns', 'Thrust-specific fuel consumption', 'kg/N/s'], ['propulsion', 'eta_overall', 'Overall propulsive efficiency', '-'],
  ['fea', 'margin_of_safety', 'Structural margin of safety', '-'], ['fea', 'wing_struct_mass_kg', 'Wing structural mass', 'kg'], ['aeroelastic', 'flutter_margin', 'Flutter margin', '-'], ['fatigue', 'life_fh', 'Fatigue life', 'FH'],
  ['flightdyn', 'static_margin', 'Static margin', '-'], ['control', 'pm_deg', 'Phase margin', 'deg'], ['gear', 'gear_load_factor', 'Landing load factor', 'g'], ['crash', 'peak_g', 'Crash peak deceleration', 'g'],
  ['acoustics', 'OASPL_dB', 'Overall sound pressure level', 'dB'], ['thermal', 'thermal_margin_K', 'Thermal margin', 'K'], ['icing', 'antiice_power_W', 'Anti-icing power', 'W'], ['electrical', 'batt_soc_end', 'Battery state of charge at end', '-'],
  ['safety', 'p_independent_per_fh', 'Catastrophic failure probability (independent failures)', '/FH'], ['safety', 'p_catastrophic_per_fh', 'Catastrophic failure probability incl. common-cause screening', '/FH'], ['safety', 'dispatch_reliability', 'Dispatch reliability', '-'], ['mission', 'block_fuel_kg', 'Mission block fuel', 'kg'], ['mission', 'co2_kg', 'Mission CO₂', 'kg'],
  ['economics', 'doc_usd_fh', 'Direct operating cost', 'USD/FH'], ['economics', 'npv_usd', 'Net present value', 'USD'], ['economics', 'irr', 'Internal rate of return', '-'], ['mdao', 'opt_improvement_pct', 'Improvement found by optimisation', '%'],
];

export async function render(root, _p, { setCrumb }) {
  setCrumb('Decision support');
  const host = h('div'); root.append(host);
  let sev = 'all', suite = 'all'; const plots = [];
  const paint = async () => {
    plots.splice(0).forEach((p) => p.destroy());
    const all = [];
    for (const [k, list] of Object.entries(state.recs)) { const [s, a] = k.split('.'); for (const r of list || []) all.push({ ...r, suite: s, analysis: a }); }
    all.sort((x, y) => (ORDER[x.severity] ?? 9) - (ORDER[y.severity] ?? 9) || suiteMeta(x.suite).d - suiteMeta(y.suite).d);
    const shown = all.filter((r) => (sev === 'all' || r.severity === sev) && (suite === 'all' || r.suite === suite)), count = (s) => all.filter((r) => r.severity === s).length;
    const notRun = SUITES.filter((s) => !Object.keys(state.runs).some((k) => k.startsWith(s.id + '.')));
    const up = state.up, c = state.case, fuel = FLUIDS[c.prop.fuel] || FLUIDS['Jet A-1'];
    // sustainability from coupled outputs
    const co2 = up.mission?.co2_kg, blockFuel = up.mission?.block_fuel_kg, energy = up.mission?.mission_energy_kWh, pax = c.mission.pax || c.econ.seats || 0, dist = c.mission.range_km || 0;
    const perPkm = co2 != null && pax && dist ? (1000 * co2) / (pax * (c.econ.load_factor || 1) * dist) : null;
    const grid = (await idb.get(`live.grid.cc:${c.site.country || 'WLD'}`))?.data?.gCO2_kWh;
    const saf = FLUIDS['SAF (HEFA-SPK)'], scen = blockFuel ? [['Current fuel', co2 ?? blockFuel * fuel.co2_per_kg], ['30% sustainable aviation fuel blend', blockFuel * (0.7 * fuel.co2_per_kg + 0.3 * saf.co2_per_kg * saf.lifecycle_factor)], ['100% sustainable aviation fuel', blockFuel * saf.co2_per_kg * saf.lifecycle_factor], ['1% drag reduction', (co2 ?? blockFuel * fuel.co2_per_kg) * 0.99], ['3% lighter empty mass', (co2 ?? blockFuel * fuel.co2_per_kg) * (1 - 0.03 * 0.6 * (c.mass.oew_kg / c.mass.mtow_kg))]] : null;
    const scoreRows = SCORE.filter(([s, k]) => typeof up[s]?.[k] === 'number').map(([s, k, label, unit]) => [`${suiteMeta(s).d}. ${suiteMeta(s).short}`, label, up[s][k], unit]);

    setKids(host, 
      h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'From numbers to action'), h('h1', null, 'Decision support'), h('p', null, 'Each analysis interprets its own results against engineering criteria and says what to do about them. This page gathers those findings across all suites, most urgent first, alongside the sustainability and cost picture.')),
        h('div', { class: 'row' }, btn('Export findings', () => downloadText('aerosuite-findings.json', JSON.stringify({ case: c.meta, generated: new Date().toISOString(), findings: all }, null, 1), 'application/json'), { ic: 'download', disabled: !all.length }))),
      !all.length && !scoreRows.length ? empty('bulb', 'No findings yet', 'Run some analyses, or the whole coupled chain, and the recommendations will be collected here.', btn('Run all suites', () => (location.hash = '#/integrated'), { ic: 'graph', kind: 'primary' })) : null,
      all.length || scoreRows.length ? h('div', { class: 'kpis' },
        kp('Act now', count('critical'), count('critical') ? 'bad' : 'ok'), kp('Needs attention', count('warn'), count('warn') ? 'warn' : 'ok'), kp('Opportunities', count('advise')), kp('Notes', count('info')), kp('Suites with results', 26 - notRun.length, '', ' / 26')) : null,
      h('div', { class: 'gap' }),
      h('div', { class: 'grid g2' },
        scoreRows.length ? card('Cross-discipline scorecard', h('div', { class: 'table-wrap' }, h('table', { class: 'data' }, h('thead', null, h('tr', null, ['Suite', 'Quantity', 'Value', 'Unit'].map((x, i) => h('th', { class: i === 2 ? 'num' : '' }, x)))), h('tbody', null, scoreRows.map((r) => h('tr', null, h('td', null, r[0]), h('td', null, r[1]), h('td', { class: 'num' }, num(r[2])), h('td', null, r[3]))))))) : null,
        card('Sustainability', h('div', { class: 'stack' },
          blockFuel || energy ? h('div', { class: 'kpis' }, blockFuel ? kp('Mission fuel', blockFuel, '', ' kg') : null, co2 != null ? kp('Mission CO₂', co2, '', ' kg') : null, perPkm != null ? kp('CO₂ per passenger-kilometre', perPkm, '', ' g') : null, energy ? kp('Mission energy', energy, '', ' kWh') : null, energy && grid ? kp('Charging CO₂ on today’s GB grid', (energy * grid) / 1000, '', ' kg') : null, up.acoustics?.OASPL_dB != null ? kp('Overall sound pressure level', up.acoustics.OASPL_dB, '', ' dB') : null, up.propulsion?.EINOx_g_kg != null ? kp('NOx emission index', up.propulsion.EINOx_g_kg, '', ' g/kg') : null)
            : h('p', { class: 'muted' }, 'Run the mission suite to quantify fuel, energy and CO₂ for this aircraft.'),
          h('div', { id: 'sus-plot' }),
          h('p', { class: 'muted small' }, 'Scenario bars scale the mission result: sustainable-fuel cases apply the life-cycle factor in the fuel database to combustion CO₂; the drag and mass cases use first-order Breguet sensitivities. Treat them as screening estimates and confirm with Suite 24 and Suite 26.')))),
      h('div', { class: 'gap' }),
      all.length ? card('Findings and recommendations', h('div', { class: 'stack' },
        h('div', { class: 'row' }, h('div', { class: 'chips', style: { paddingBottom: 0 } }, [['all', `All (${all.length})`], ['critical', `Act now (${count('critical')})`], ['warn', `Attention (${count('warn')})`], ['advise', `Opportunities (${count('advise')})`], ['info', `Notes (${count('info')})`]].map(([k, l]) => h('button', { class: `chip ${sev === k ? 'on' : ''}`, onclick: () => { sev = k; paint(); } }, l))),
          h('select', { class: 'inp', style: { width: 'auto' }, 'aria-label': 'Filter by suite', onchange: (e) => { suite = e.target.value; paint(); } }, h('option', { value: 'all' }, 'All suites'), [...new Set(all.map((r) => r.suite))].map((s) => h('option', { value: s, selected: suite === s }, `${suiteMeta(s).d}. ${suiteMeta(s).short}`)))),
        shown.length ? shown.map((r) => recCard(r, { hash: `#/suite/${r.suite}/${r.analysis}`, label: `${suiteMeta(r.suite).d}. ${suiteMeta(r.suite).short} →` })) : h('p', { class: 'muted' }, 'Nothing matches this filter.'))) : null,
      notRun.length && notRun.length < 26 ? card('Suggested next steps', h('div', { class: 'stack' }, h('p', { class: 'muted small' }, 'These suites have no results yet for this case, so their findings are missing from the picture above.'), h('div', { class: 'row' }, notRun.map((s) => h('a', { class: 'chip', href: `#/suite/${s.id}` }, `${s.d}. ${s.short}`)))), { collapsible: true }) : null,
      h('div', { class: 'gap' }),
      h('div', { class: 'note' }, icon('info'), h('div', null, 'Recommendations follow from the models and inputs stated in each analysis and name the criterion they rest on. Prices, rates and reference conditions come from the live feeds where available', Object.values(liveStatus).some((v) => v.ok) ? '' : ' (none reachable yet on this device)', '. They support engineering judgement; they do not replace it.')));
    if (scen) { const el = host.querySelector('#sus-plot'); if (el) plots.push(renderPlot(el, { type: 'bar', title: 'Mission CO₂ under improvement scenarios', ylabel: 'CO₂ per mission [kg]', categories: scen.map((s) => s[0]), series: [{ name: 'CO₂', y: scen.map((s) => s[1]) }] })); }
  };
  const kp = (label, value, st = '', unit = '') => h('div', { class: `kpi ${st}` }, h('span', { class: 'l' }, label), h('div', { class: 'v' }, num(value), unit ? h('small', null, unit) : null));
  await paint();
  const off = on('run', paint);
  return () => { off(); plots.forEach((p) => p.destroy()); };
}
