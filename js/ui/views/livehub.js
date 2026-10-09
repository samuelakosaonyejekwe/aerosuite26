// Live data hub: every external feed, what it is used for, when it was last refreshed, and its values.

import { h, icon, clear, add, num, btn, badge, card, toast, ago } from '../dom.js';
import { CONNECTORS, live, status, refreshSite, refreshGlobal, jetFromBrent, loadSnapshot, seedFromSnapshot, carbonPrice, carbonChoice, config, loadConfig, loadLicences, allowed, attributions } from '../../core/live.js';
import { state, setSetting, on, idb } from '../../core/store.js';
import { renderPlot } from '../plots.js';
import { gridIndex } from '../../core/gridwx.js';
import { dataTable } from '../results.js';

const CLASS_WORD = { 'commercial-ok': 'Commercial use allowed', 'commercial-with-key': 'Commercial use with a key, own server or licence', unclear: 'Unclear — no explicit grant', 'not-allowed': 'Not allowed without consent' };
/** Attribution lines required by the providers of the given registry sources (shown next to the data they belong to). */
export async function attributionLine(ids, lead = 'Data: ') {
  const a = await attributions(ids); if (!a.length) return null;
  return h('p', { class: 'muted small attribution' }, lead, a.map((x, i) => [i ? ' · ' : '', x.url ? h('a', { href: x.url, target: '_blank', rel: 'noopener noreferrer' }, x.text) : x.text]));
}
/**
 * "Data sources and licences": the audited registry (js/data/licences.json) as tables — every external source and
 * dataset with its terms, the quoted basis, its class and whether this deployment uses it; then the third-party
 * software components. Shown on the Live data page and the About page.
 */
/** Printable compliance pack (the content of COMPLIANCE.md): per mode, what is contacted or shipped, on what terms, and the open questions. */
async function compliancePack(reg) {
  const { complianceRows, complianceQuestions } = await import('../../core/compliance.js'), qs = complianceQuestions(reg), host = h('div', { class: 'stack compliance-pack' });
  const tbl = (mode, title) => dataTable({ title, columns: ['Source', 'How it is used', 'Licence or terms', 'Obligations', 'Attribution shown at', 'Sentence relied on', 'Address', 'Read on'], rows: complianceRows(reg, mode).map((r) => [r.name, r.use, r.licence, r.obligations, r.shownAt, r.quote, r.url, r.retrieved]) });
  const print = () => window.print(); // the page's print styles hide the navigation; the open pack prints with the page
  host.append(h('h4', null, `Compliance pack — terms read on ${reg.retrieved}`),
    h('p', { class: 'small' }, 'For each deployment mode: every source that is contacted or shipped, its licence, the obligation, where its attribution is shown, and the sentence relied on with its address. A record for review, not legal advice. The same content is in ', h('a', { href: 'COMPLIANCE.md', target: '_blank', rel: 'noopener' }, 'COMPLIANCE.md'), '.'),
    h('div', null, h('b', null, `Questions that need a decision (${qs.length})`), qs.length ? h('ol', { class: 'small' }, qs.map((q) => h('li', null, h('b', null, q.name), ' — ', q.question))) : h('p', { class: 'small' }, 'None.')),
    tbl('commercial', 'Commercial deployment (no keys)'), tbl('nonCommercial', 'Non-commercial deployment (default configuration)'));
  return h('details', { class: 'compliance' }, h('summary', { style: { cursor: 'pointer' } }, h('b', null, 'Compliance pack'), ' — printable record for legal review'), h('div', { class: 'row', style: { margin: '8px 0' } }, btn('Print', print, { ic: 'doc', kind: 'sm' })), host);
}
export async function licenceSection({ pack = false } = {}) {
  await loadConfig();
  let reg; try { reg = await loadLicences(); } catch (e) { return card('Data sources and licences', h('p', { class: 'muted' }, `The licence registry could not be read (${e.message}).`)); }
  const inUse = async (x) => (x.used === false ? 'No — removed' : (await allowed(x.id)) ? (!config.commercial && x.usedInCommercial !== true ? 'Yes (non-commercial deployment)' : 'Yes') : 'No');
  const rows = []; for (const x of reg.sources) rows.push([x.name, x.usedFor, x.licence, CLASS_WORD[x.class] || x.class, await inUse(x), (x.obligations || []).join('; '), x.attribution || '', x.replacement || '', x.quote ? `“${x.quote}”` : '', x.url || '', x.retrieved]);
  const comp = reg.components.map((c) => [c.name, c.version, c.licence, c.shipped ? 'Shipped with the app' : 'Development tool, not shipped', c.licenceFile, c.modifications, (c.obligations || []).join('; '), c.source]);
  const required = [...new Set(reg.sources.filter((x) => x.attribution && x.used !== false).map((x) => x.attribution))];
  return card('Data sources and licences', h('div', { class: 'stack' },
    h('div', { class: `note ${config.commercial ? 'ok' : ''}` }, icon(config.commercial ? 'check' : 'info'), h('div', null, h('b', null, config.commercial ? 'Commercial deployment. ' : 'Non-commercial deployment. '),
      config.commercial ? 'This copy contacts only sources whose terms allow commercial use, sources for which its operator configured a key or address, and sources its operator has licensed. The others are switched off and say so on their cards.'
        : 'This copy may also use sources whose free tiers are limited to non-commercial use (marked below). An operator who sells or monetises the app sets "commercial": true in config.json, which switches every feed to the compliant set.',
      ` Terms were read on ${reg.retrieved}; each line gives the page and the sentence relied on. This table is a compliance record, not legal advice.`)),
    dataTable({ title: 'External data sources and datasets', columns: ['Source', 'Used for', 'Licence or terms', 'Class', 'Used in this deployment', 'Obligations', 'Attribution shown', 'Replacement when not usable', 'Sentence relied on', 'Terms page', 'Read on'], rows }),
    dataTable({ title: 'Third-party software components', columns: ['Component', 'Version', 'Licence', 'Distribution', 'Licence text', 'Modifications', 'Obligations', 'Source code'], rows: comp }),
    h('p', { class: 'small' }, 'The application itself is proprietary (see ', h('a', { href: 'LICENSE', target: '_blank', rel: 'noopener' }, 'LICENSE'), '). Third-party components keep their own licences: see ', h('a', { href: 'THIRD_PARTY_NOTICES.md', target: '_blank', rel: 'noopener' }, 'THIRD_PARTY_NOTICES.md'), ', which also explains how to replace the LGPL-licensed OpenCASCADE kernel (a separate file loaded at run time).'),
    h('details', null, h('summary', { class: 'small', style: { cursor: 'pointer' } }, 'Attribution statements'), h('ul', { class: 'small' }, required.map((t) => h('li', null, t)))),
    pack ? await compliancePack(reg) : null));
}

const paramsFor = (id) => { const s = state.case.site; if (CONNECTORS[id].needs === 'site') return s.lat == null ? null : { lat: s.lat, lon: s.lon }; if (id === 'grid') return { country: s.country || '' }; if (id === 'macro') return { country: s.country || 'WLD' }; if (id === 'literature' || id === 'opensource') return null; return {}; };
const fmtTtl = (ms) => (ms < 3600e3 ? `${Math.round(ms / 60e3)} min` : ms < 86400e3 ? `${Math.round(ms / 3600e3)} h` : `${Math.round(ms / 86400e3)} d`);

export async function render(root, _p, { setCrumb }) {
  setCrumb('Live data');
  const host = h('div'); root.append(host);
  const plots = [], list = h('div', { class: 'grid g2' }), oilHost = h('div'), snapHost = h('div'), bgHost = h('span'), licHost = h('div', { id: 'licences' });
  let fxNow = null;
  const inCase = (label, value, unit, source) => h('p', { class: 'small in-case' }, h('b', null, `${label} in the case: `), `${num(value)} ${unit}`, source ? h('span', { class: 'muted' }, ` — ${source}`) : h('span', { class: 'muted' }, ' — entered by hand or preset'));
  const preview = (id, d) => {
    if (!d) return h('span', { class: 'muted' }, 'No data yet');
    const kv = (pairs) => h('dl', { class: 'kv small' }, pairs.filter((p) => p[1] != null).map(([k, v]) => [h('dt', null, k), h('dd', null, typeof v === 'number' ? num(v) : String(v))]));
    switch (id) {
      case 'weather': return h('div', { class: 'stack' }, kv([['Provider', d.provider_label], ['Valid', d.time + 'Z'], ['Temperature [°C]', d.T_C], ['Pressure [hPa]', d.p_hPa], ['Wind [m/s]', d.wind_ms], ['Gust [m/s]', d.gust_ms ?? 'not supplied'], [`Freezing level${d.freezing_level_estimated ? ' (ISA estimate)' : ''} [m]`, d.freezing_level_m], ['Levels aloft', (d.aloft || []).length || 'not supplied']]), d.note ? h('p', { class: 'small wx-note' }, d.note) : null);
      case 'climate': return h('div', { class: 'stack' }, kv([['Period', d.from === d.to ? d.from : `${d.from} to ${d.to}`], ['Hot day, 99th pct [°C]', d.hot99_C], ['Cold day, 1st pct [°C]', d.cold01_C], ['Strong wind, 99th pct [m/s]', d.wind99_ms]]), d.basis ? h('p', { class: 'muted small' }, d.basis) : null);
      case 'airquality': return h('div', { class: 'stack' }, kv([['PM10 [µg/m³]', d.pm10], ['PM2.5 [µg/m³]', d.pm2_5], ['Dust [µg/m³]', d.dust], ['Aerosol optical depth', d.aod], ['NO₂ [µg/m³]', d.no2], ['Ozone [µg/m³]', d.o3]]), d.note ? h('p', { class: 'muted small' }, d.note) : null);
      case 'marine': return kv([['Source', d.provider === 'dwd-gwam' ? 'DWD wave model (2° cells)' : d.provider === 'open-meteo' ? 'Open-Meteo' : null], ['Wave height [m]', d.wave_height_m ?? 'inland / none'], ['Wave period [s]', d.wave_period_s], ['Wave direction, from [°]', d.wave_dir_deg]]);
      case 'aerodromes': return kv([['Aerodromes within 40 km', d.fields.length], ['Runways', d.nRunways], ['Nearest', d.fields[0]?.name], ['Longest runway nearby [m]', Math.max(0, ...d.fields.map((f) => f.longest_m || 0)) || null], ['Source', d.source], ['Database date', d.dataset_date]]);
      case 'spaceweather': return kv([['Planetary Kp', d.kp], ['At', d.time], ['GNSS outlook', d.kp >= 7 ? 'severe storm — expect degraded accuracy' : d.kp >= 5 ? 'storm — monitor integrity' : 'quiet']]);
      case 'fx': return kv([['Source', d.source ? `${d.source} euro reference rates; per-USD rates are own calculations` : null], ['Date', d.date], ['EUR per USD', d.rates.EUR], ['GBP per USD', d.rates.GBP], ['JPY per USD', d.rates.JPY], ['Currencies', Object.keys(d.rates).length]]);
      case 'oil': return kv([['Source', 'U.S. Energy Information Administration'], ['Date', d.date], ['Brent [USD/bbl]', d.brent_usd_bbl], ['Jet fuel estimated from Brent [USD/kg]', jetFromBrent(d.brent_usd_bbl, state.settings.fuelCrack)], ['Estimate is used', 'only when the quoted jet-fuel series is missing or old'], ['Annualised volatility', d.vol_annual]]);
      case 'jetfuel': return h('div', { class: 'stack' }, kv([['Price date', d.date], ['US Gulf Coast jet fuel, spot [USD/gal]', d.usd_gal], ['Converted [USD/kg]', d.usd_kg], ['Density used [kg/L]', d.density_kg_l], ['Annualised volatility', d.vol_annual], ['Series', d.series], ['Last scheduled fetch', d._error ? `failed (${d._error}); showing the previous value` : null]]), inCase('Fuel price', state.case.econ.fuel_usd_kg, 'USD/kg', state.case.econ.fuel_source));
      case 'carbon': {
        if (!d.markets) return h('span', { class: 'muted' }, 'No data yet');
        const pref = state.settings.carbonMarket || 'auto', used = carbonChoice(d, state.case.site.country, pref), ids = (d.order || Object.keys(d.markets)).filter((id) => d.markets[id]);
        const usd = (m) => carbonPrice(m, fxNow)?.usd_t ?? null;
        const pick = h('select', { class: 'inp', style: { width: 'auto' }, 'aria-label': 'Carbon market used in the case', onchange: (e) => { setSetting('carbonMarket', e.target.value); refreshGlobal().then(paintList); } },
          h('option', { value: 'auto', selected: pref === 'auto' }, `Automatic: by the site’s country${used ? ` (now ${used.market})` : ''}`), ids.map((id) => h('option', { value: id, selected: pref === id }, d.markets[id].market)));
        return h('div', { class: 'stack' },
          h('div', { class: 'table-wrap' }, h('table', { class: 'data carbon-markets' }, h('thead', null, h('tr', null, ['Market', 'Price', 'USD/t', 'Date', 'Publisher'].map((x, k) => h('th', { class: k === 1 || k === 2 ? 'num' : '' }, x)))),
            h('tbody', null, ids.map((id) => { const m = d.markets[id], u = usd(m); return h('tr', { class: used && used.id === m.id ? 'on' : '' }, h('td', null, used && used.id === m.id ? h('b', null, m.market) : m.market), h('td', { class: 'num' }, `${num(m.price)} ${m.currency}`), h('td', { class: 'num' }, u == null ? '–' : num(u)), h('td', null, m.period || (m.scheme_year ? `year ${m.scheme_year}` : m.date)), h('td', null, m.url ? h('a', { href: m.url, target: '_blank', rel: 'noopener noreferrer' }, m.source) : m.source)); })))),
          h('label', { class: 'row small' }, 'Market used in the case', pick),
          used ? h('p', { class: 'muted small' }, used.instrument, used.lag_note ? ` ${used.lag_note}` : '') : null,
          Object.keys(d.errors || {}).length ? h('p', { class: 'small' }, `Not updated in the last run: ${Object.entries(d.errors).map(([k, v]) => `${k} (${v})`).join('; ')}.`) : null,
          inCase('Carbon price', state.case.econ.carbon_usd_t, 'USD/t CO₂', state.case.econ.carbon_source)); }
      case 'rates': return kv([['US effective federal funds rate [%]', d.us_fed_funds_pct], ['as of', d.us_date], ['Source', d.us_via], ['ECB deposit facility rate [%]', d.ecb_deposit_pct], ['as of ', d.ecb_date]]);
      case 'macro': return kv([['Economy', d.country], [`Inflation ${d.inflation_year || ''} [%]`, d.inflation_pct], [`Lending rate ${d.lending_year || ''} [%]`, d.lending_pct], [`GDP growth ${d.gdp_year || ''} [%]`, d.gdp_growth_pct]]);
      case 'grid': return kv([['Economy', d.country], ['Carbon intensity [gCO₂e/kWh]', d.gCO2_kWh], ['Basis', d.basis], ['Level', d.index], ['National average [gCO₂e/kWh]', d.from ? d.national_g : null]]);
      default: return h('span', { class: 'muted small' }, `${d.items?.length || 0} results`);
    }
  };
  const paintList = async () => {
    await loadConfig();
    const cards = []; fxNow = (await idb.get('live.fx.usd'))?.data || null;
    for (const [id, c] of Object.entries(CONNECTORS)) {
      const st = status[id] || {}, p = paramsFor(id), perSuite = id === 'literature' || id === 'opensource';
      const rec = p ? await idb.get(`live.${id}.${c.key(p)}`) : null, age = rec ? Date.now() - rec.ts : null, fresh = age != null && age < c.ttl;
      cards.push(card(c.title, h('div', { class: 'stack' },
        h('div', { class: 'row small' }, perSuite ? badge('Per suite', 'accent') : !p ? badge('Needs a location', 'warn') : st.policy && !rec ? badge('Not part of this deployment', 'warn') : c.bundled && rec ? badge(rec.data?.osm ? 'Stored in the app + OpenStreetMap' : 'Stored in the app — works offline', 'ok') : fresh ? badge('Fresh', 'ok') : rec ? badge(c.cloud && rec.via !== 'direct' ? 'Older than expected — snapshot not updated' : 'Stale — refreshing when online', 'warn') : st.ok === false ? badge(st.error || 'Unavailable', 'bad') : badge('Not fetched yet'), c.cloud && rec?.via !== 'direct' ? badge('Cloud snapshot only', 'accent') : rec?.via === 'snapshot' ? badge('From cloud snapshot') : null,
          h('span', { class: 'muted' }, rec ? `${c.bundled ? 'Looked up' : rec.via !== 'direct' && (c.cloud || rec.via === 'snapshot') ? 'Fetched by the snapshot job' : 'Updated'} ${ago(rec.ts)}` : '', c.bundled ? '' : c.cloud && rec?.via !== 'direct' ? ' · snapshot rebuilt about every 3 h' : ` · refreshes every ${fmtTtl(c.ttl)}`)),
        h('p', { class: 'muted small' }, c.use), perSuite ? h('p', { class: 'small' }, 'Open the “Live resources” tab inside any suite.') : st.policy && !rec ? h('p', { class: 'small' }, st.error) : preview(id, rec?.data),
        config.attribution ? await attributionLine(c.reg ? c.reg() : []) : null,
        h('div', { class: 'row small' }, h('a', { href: c.home, target: '_blank', rel: 'noopener noreferrer' }, c.provider, ' ', icon('external', 13)),
          p && !perSuite ? btn('Refresh', async () => { const r = await live(id, p, { force: true }); toast(r.error ? `${c.title}: ${r.error}` : `${c.title} updated.`, r.error ? 'bad' : 'ok'); paintList(); }, { ic: 'refresh', kind: 'sm ghost' }) : !p && !perSuite ? btn('Choose location', () => (location.hash = '#/case/site'), { kind: 'sm ghost', ic: 'pin' }) : null))));
    }
    clear(list); add(list, cards);
    const snap = await loadSnapshot().catch(() => null); clear(snapHost);
    if (snap) { const fs = Object.entries(snap.feeds || {}), bad = fs.filter(([, f]) => !f.ok), age = Date.now() - Date.parse(snap.generated); snapHost.append(h('div', { class: 'stack' }, h('div', { class: 'row small' }, badge(age < 8 * 3600e3 ? 'Current' : 'Not updated recently', age < 8 * 3600e3 ? 'ok' : 'warn'), h('span', { class: 'muted' }, `Built ${ago(Date.parse(snap.generated))} · ${fs.length - bad.length} of ${fs.length} feeds fetched${snap.from ? ` · read from ${new URL(snap.from).host}` : ' · saved copy on this device'}`)),
      bad.length ? h('p', { class: 'small' }, `Not fetched in the last run: ${bad.map(([id, f]) => `${CONNECTORS[id]?.title || id} (${f.error})`).join('; ')}. The previous values are kept.`) : null,
      age >= 8 * 3600e3 ? h('p', { class: 'small muted' }, 'The snapshot is rebuilt only by the host that runs the schedule. A copy of the app on another host, or a host whose schedule is switched off, carries an older file; feeds your browser can reach directly are still refreshed.') : null)); }
    else snapHost.append(h('p', { class: 'muted small' }, 'No cloud snapshot could be read (offline at first start, or this copy of the app was opened from a file). Feeds your browser can reach directly still work; the jet-fuel spot price and policy rates need the snapshot, so the case keeps its entered values for those.'));
    if (globalThis.__AEROSUITE_STANDALONE__) { const c = globalThis.__AEROSUITE_ESSENTIAL__?.carries; snapHost.append(h('div', { class: 'note standalone-data' }, icon('info'), h('div', null, h('b', null, 'Single-file copy. '), c ? `This file carries its own data, as of its build: ${num(c.airports)} large and medium airports with their runways, ${num(c.cities)} cities of 100 000 inhabitants or more, the design-temperature grid, national grid emission factors, the tool lists${c.snapshot ? `, and the cloud snapshot of ${c.snapshot.slice(0, 10)} (carbon prices, exchange rates, fuel prices, policy rates)` : ''}. ` : 'This file carries no bundled data. ', 'Smaller airports and places, terrain, and the forecast grids (winds aloft, sea state, air quality) need the app served from a web address; with a connection, the live feeds your browser can reach still update.'))); }
    gridIndex().then((g) => { const P = g.products || {}, label = { wx: 'weather (NOAA GFS)', sea: 'sea state (DWD wave model)', air: 'aerosols (NOAA GEFS-Aerosols)' }; snapHost.append(h('p', { class: 'small grid-status' }, h('b', null, 'Forecast grids: '), Object.entries(label).map(([id, t], k) => `${k ? ' · ' : ''}${t} ${P[id]?.times?.length ? `run ${P[id].cycle.slice(0, 13)}Z, valid to ${P[id].times[P[id].times.length - 1].valid.slice(0, 16)}Z${P[id].ok === false ? ' (not updated in the last run)' : ''}` : 'not available'}`), '. Stored on this device for offline use.')); }).catch(() => snapHost.append(h('p', { class: 'small muted grid-status' }, 'Forecast grids could not be read (offline before they were stored, or the app was opened from a file).')));
    clear(licHost); licHost.append(await licenceSection());
    const oil = (await idb.get('live.oil.brent'))?.data;
    plots.splice(0).forEach((x) => x.destroy()); clear(oilHost);
    if (oil?.history) plots.push(renderPlot(oilHost, { type: 'line', title: 'Brent crude, last 12 months', xlabel: 'Trading weeks ago', ylabel: 'Price [USD/bbl]', series: [{ name: 'Brent', x: oil.history.map((_, i) => i - oil.history.length + 1), y: oil.history.map((r) => r.v) }] }));
  };
  host.append(
    h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Always-fresh inputs'), h('h1', null, 'Live data'), h('p', null, 'Three sources keep the inputs current without depending on any one computer. Airports and runways worldwide are stored inside the app. Weather, exchange rates, oil, space weather, economic indicators, research and code indexes are fetched by your browser straight from their public providers. A snapshot of the location-independent feeds — including a carbon price, a quoted jet-fuel price and policy rates, whose providers do not accept browser requests — is rebuilt on the hosting service about every three hours. Everything is saved on your device for offline use. Each card names its provider and carries the attribution its licence asks for; the full list of sources, terms and licences is at the foot of this page.')),
      h('div', { class: 'row' }, btn('Refresh everything', async () => { toast('Refreshing live data…', 'info'); await seedFromSnapshot({ force: true }).catch(() => 0); await Promise.all([refreshSite({ force: true }), refreshGlobal({ force: true })]); toast('Live data refreshed.', 'ok'); paintList(); }, { ic: 'refresh', kind: 'primary' }))),
    card('Settings', h('div', { class: 'grid g3' },
      h('label', { class: 'switch' }, h('input', { type: 'checkbox', checked: state.settings.autoRefreshLive, onchange: (e) => setSetting('autoRefreshLive', e.target.checked) }), 'Refresh automatically while the app is open'),
      h('label', { class: 'switch' }, h('input', { type: 'checkbox', checked: state.settings.autoApplyLive, onchange: (e) => setSetting('autoApplyLive', e.target.checked) }), 'Write live values into the case'),
      h('label', { class: 'row small' }, 'Jet-fuel refining margin over Brent', h('input', { class: 'inp', type: 'number', style: { width: '80px' }, min: 0, max: 200, step: 1, value: Math.round(state.settings.fuelCrack * 100), onchange: (e) => { setSetting('fuelCrack', Math.max(0, Number(e.target.value) / 100)); refreshGlobal().then(paintList); } }), '%'))),
    h('div', { class: 'gap' }), card('Cloud snapshot', snapHost), h('div', { class: 'gap' }), list, h('div', { class: 'gap' }), oilHost,
    h('div', { class: 'gap' }), licHost,
    h('div', { class: 'gap' }),
    h('div', { class: 'note' }, icon('info'), h('div', null, h('b', null, 'What “live” means here. '), 'A feed is as current as its provider publishes: surface weather every 15 minutes, space weather every 3 hours, exchange rates and central-bank rates each working day, oil and jet-fuel prices each trading day (published a few days in arrears), EU allowance auctions on most working days (the UK ETS figure used by commercial deployments once a year), World Bank indicators yearly. Carbon and jet-fuel values are public reference prices, not quotations for a particular contract or airport.')),
    h('div', { class: 'gap' }),
    h('div', { class: 'note' }, icon('info'), h('div', null, h('b', null, 'While the app is closed. '), 'A web app cannot run on your device while it is fully closed. Two things cover the gap. First, the cloud snapshot is rebuilt on the hosting service about every 3 hours whether or not anyone has the app open (scheduled jobs can start late, and stop if the host disables the schedule), so the first screen after opening shows location-independent data that is normally no more than a few hours old. Second, browsers with periodic background sync — currently the installed app in Chrome and Edge on Android and desktop — also refresh the saved feeds a few times a day at their own discretion; Safari on iPhone, iPad and Mac, and Firefox, do not offer this. Weather for your site depends on the location, so it is fetched when the app is opened, brought back to the foreground, or reconnects, and every few minutes while it stays open. Offline, the last saved values are used and marked as such. ', bgHost)));
  (async () => { let txt = 'On this device: background refresh while closed is not available; data refreshes when the app is opened.'; try { const reg = await Promise.race([navigator.serviceWorker.ready, new Promise((r) => setTimeout(r, 1500))]); if (reg?.periodicSync && (await reg.periodicSync.getTags()).includes('refresh-live')) txt = 'On this device: background refresh while closed is switched on.'; else if (reg?.periodicSync) txt = 'On this device: the browser supports background refresh but has not granted it (it is offered to installed apps that are used regularly).'; } catch { /* keep the default */ } bgHost.textContent = txt; })();
  await paintList();
  const off = on('live', () => paintList());
  return () => { off(); plots.forEach((p) => p.destroy()); };
}
