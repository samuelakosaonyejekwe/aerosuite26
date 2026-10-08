// Live data hub: every external feed, what it is used for, when it was last refreshed, and its values.

import { h, icon, clear, add, num, btn, badge, card, toast, ago } from '../dom.js';
import { CONNECTORS, live, status, refreshSite, refreshGlobal, jetFromBrent } from '../../core/live.js';
import { state, setSetting, on, idb } from '../../core/store.js';
import { renderPlot } from '../plots.js';

const paramsFor = (id) => { const s = state.case.site; if (CONNECTORS[id].needs === 'site') return s.lat == null ? null : { lat: s.lat, lon: s.lon }; if (id === 'macro') return { country: s.country || 'WLD' }; if (id === 'literature' || id === 'opensource') return null; return {}; };
const fmtTtl = (ms) => (ms < 3600e3 ? `${Math.round(ms / 60e3)} min` : ms < 86400e3 ? `${Math.round(ms / 3600e3)} h` : `${Math.round(ms / 86400e3)} d`);

export async function render(root, _p, { setCrumb }) {
  setCrumb('Live data');
  const host = h('div'); root.append(host);
  const plots = [], list = h('div', { class: 'grid g2' }), oilHost = h('div');
  const preview = (id, d) => {
    if (!d) return h('span', { class: 'muted' }, 'No data yet');
    const kv = (pairs) => h('dl', { class: 'kv small' }, pairs.filter((p) => p[1] != null).map(([k, v]) => [h('dt', null, k), h('dd', null, typeof v === 'number' ? num(v) : String(v))]));
    switch (id) {
      case 'weather': return kv([['Observed', d.time + 'Z'], ['Temperature [°C]', d.T_C], ['Pressure [hPa]', d.p_hPa], ['Wind [m/s]', d.wind_ms], ['Gust [m/s]', d.gust_ms], ['Freezing level [m]', d.freezing_level_m], ['Levels aloft', d.aloft.length]]);
      case 'climate': return kv([['Period', `${d.from} to ${d.to}`], ['Hot day, 99th pct [°C]', d.hot99_C], ['Cold day, 1st pct [°C]', d.cold01_C], ['Strong wind, 99th pct [m/s]', d.wind99_ms]]);
      case 'airquality': return kv([['PM10 [µg/m³]', d.pm10], ['PM2.5 [µg/m³]', d.pm2_5], ['Dust [µg/m³]', d.dust], ['Aerosol optical depth', d.aod]]);
      case 'marine': return kv([['Wave height [m]', d.wave_height_m ?? 'inland / none'], ['Wave period [s]', d.wave_period_s]]);
      case 'aerodromes': return kv([['Aerodromes within 40 km', d.fields.length], ['Runways mapped', d.nRunways], ['Nearest', d.fields[0]?.name]]);
      case 'spaceweather': return kv([['Planetary Kp', d.kp], ['At', d.time], ['GNSS outlook', d.kp >= 7 ? 'severe storm — expect degraded accuracy' : d.kp >= 5 ? 'storm — monitor integrity' : 'quiet']]);
      case 'fx': return kv([['Date', d.date], ['EUR per USD', d.rates.EUR], ['GBP per USD', d.rates.GBP], ['JPY per USD', d.rates.JPY], ['Currencies', Object.keys(d.rates).length]]);
      case 'oil': return kv([['Date', d.date], ['Brent [USD/bbl]', d.brent_usd_bbl], ['Derived jet fuel [USD/kg]', jetFromBrent(d.brent_usd_bbl, state.settings.fuelCrack)], ['Annualised volatility', d.vol_annual]]);
      case 'macro': return kv([['Economy', d.country], [`Inflation ${d.inflation_year || ''} [%]`, d.inflation_pct], [`Lending rate ${d.lending_year || ''} [%]`, d.lending_pct], [`GDP growth ${d.gdp_year || ''} [%]`, d.gdp_growth_pct]]);
      case 'grid': return kv([['Carbon intensity [gCO₂/kWh]', d.gCO2_kWh], ['Level', d.index], ['Period from', d.from]]);
      default: return h('span', { class: 'muted small' }, `${d.items?.length || 0} results`);
    }
  };
  const paintList = async () => {
    const cards = [];
    for (const [id, c] of Object.entries(CONNECTORS)) {
      const st = status[id] || {}, p = paramsFor(id), perSuite = id === 'literature' || id === 'opensource';
      const rec = p ? await idb.get(`live.${id}.${c.key(p)}`) : null, age = rec ? Date.now() - rec.ts : null, fresh = age != null && age < c.ttl;
      cards.push(card(c.title, h('div', { class: 'stack' },
        h('div', { class: 'row small' }, perSuite ? badge('Per suite', 'accent') : !p ? badge('Needs a location', 'warn') : fresh ? badge('Fresh', 'ok') : rec ? badge('Stale — refreshing when online', 'warn') : st.ok === false ? badge(st.error || 'Unavailable', 'bad') : badge('Not fetched yet'), h('span', { class: 'muted' }, rec ? `Updated ${ago(rec.ts)}` : '', ` · refreshes every ${fmtTtl(c.ttl)}`)),
        h('p', { class: 'muted small' }, c.use), perSuite ? h('p', { class: 'small' }, 'Open the “Live resources” tab inside any suite.') : preview(id, rec?.data),
        h('div', { class: 'row small' }, h('a', { href: c.home, target: '_blank', rel: 'noopener noreferrer' }, c.provider, ' ', icon('external', 13)),
          p && !perSuite ? btn('Refresh', async () => { const r = await live(id, p, { force: true }); toast(r.error ? `${c.title}: ${r.error}` : `${c.title} updated.`, r.error ? 'bad' : 'ok'); paintList(); }, { ic: 'refresh', kind: 'sm ghost' }) : !p && !perSuite ? btn('Choose location', () => (location.hash = '#/case/site'), { kind: 'sm ghost', ic: 'pin' }) : null))));
    }
    clear(list); add(list, cards);
    const oil = (await idb.get('live.oil.brent'))?.data;
    plots.splice(0).forEach((x) => x.destroy()); clear(oilHost);
    if (oil?.history) plots.push(renderPlot(oilHost, { type: 'line', title: 'Brent crude, last 12 months', xlabel: 'Trading weeks ago', ylabel: 'Price [USD/bbl]', series: [{ name: 'Brent', x: oil.history.map((_, i) => i - oil.history.length + 1), y: oil.history.map((r) => r.v) }] }));
  };
  host.append(
    h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Always-fresh inputs'), h('h1', null, 'Live data'), h('p', null, 'Your browser fetches these feeds directly from their public providers — weather services, mapping, space-weather, central-bank and energy data, research and code indexes. Nothing passes through a server belonging to this app, so freshness never depends on any one machine being switched on. Each feed is saved on your device for offline use and refreshed on its own schedule whenever you are connected.')),
      h('div', { class: 'row' }, btn('Refresh everything', async () => { toast('Refreshing live data…', 'info'); await Promise.all([refreshSite({ force: true }), refreshGlobal({ force: true })]); toast('Live data refreshed.', 'ok'); paintList(); }, { ic: 'refresh', kind: 'primary' }))),
    card('Settings', h('div', { class: 'grid g3' },
      h('label', { class: 'switch' }, h('input', { type: 'checkbox', checked: state.settings.autoRefreshLive, onchange: (e) => setSetting('autoRefreshLive', e.target.checked) }), 'Refresh automatically while the app is open'),
      h('label', { class: 'switch' }, h('input', { type: 'checkbox', checked: state.settings.autoApplyLive, onchange: (e) => setSetting('autoApplyLive', e.target.checked) }), 'Write live values into the case'),
      h('label', { class: 'row small' }, 'Jet-fuel refining margin over Brent', h('input', { class: 'inp', type: 'number', style: { width: '80px' }, min: 0, max: 200, step: 1, value: Math.round(state.settings.fuelCrack * 100), onchange: (e) => { setSetting('fuelCrack', Math.max(0, Number(e.target.value) / 100)); refreshGlobal().then(paintList); } }), '%'))),
    h('div', { class: 'gap' }), list, h('div', { class: 'gap' }), oilHost,
    h('div', { class: 'gap' }),
    h('div', { class: 'note' }, icon('info'), h('div', null, h('b', null, 'What “live” means here. '), 'A feed is as current as its provider publishes: surface weather every 15 minutes, space weather every 3 hours, exchange rates each working day, oil prices each trading day, World Bank indicators yearly. When the installed app is closed, browsers that support background sync refresh the feeds a few times a day; otherwise they refresh the moment you open it. Offline, the last saved values are used and marked as such.')));
  await paintList();
  const off = on('live', () => paintList());
  return () => { off(); plots.forEach((p) => p.destroy()); };
}
