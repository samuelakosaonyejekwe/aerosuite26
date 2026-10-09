// Live data hub: every external feed, what it is used for, when it was last refreshed, and its values.

import { h, icon, clear, add, num, btn, badge, card, toast, ago } from '../dom.js';
import { CONNECTORS, live, status, refreshSite, refreshGlobal, jetFromBrent, loadSnapshot, seedFromSnapshot, fuelPrice, carbonPrice } from '../../core/live.js';
import { state, setSetting, on, idb } from '../../core/store.js';
import { renderPlot } from '../plots.js';

const paramsFor = (id) => { const s = state.case.site; if (CONNECTORS[id].needs === 'site') return s.lat == null ? null : { lat: s.lat, lon: s.lon }; if (id === 'macro') return { country: s.country || 'WLD' }; if (id === 'literature' || id === 'opensource') return null; return {}; };
const fmtTtl = (ms) => (ms < 3600e3 ? `${Math.round(ms / 60e3)} min` : ms < 86400e3 ? `${Math.round(ms / 3600e3)} h` : `${Math.round(ms / 86400e3)} d`);

export async function render(root, _p, { setCrumb }) {
  setCrumb('Live data');
  const host = h('div'); root.append(host);
  const plots = [], list = h('div', { class: 'grid g2' }), oilHost = h('div'), snapHost = h('div'), bgHost = h('span');
  let fxNow = null;
  const inCase = (label, value, unit, source) => h('p', { class: 'small in-case' }, h('b', null, `${label} in the case: `), `${num(value)} ${unit}`, source ? h('span', { class: 'muted' }, ` — ${source}`) : h('span', { class: 'muted' }, ' — entered by hand or preset'));
  const preview = (id, d) => {
    if (!d) return h('span', { class: 'muted' }, 'No data yet');
    const kv = (pairs) => h('dl', { class: 'kv small' }, pairs.filter((p) => p[1] != null).map(([k, v]) => [h('dt', null, k), h('dd', null, typeof v === 'number' ? num(v) : String(v))]));
    switch (id) {
      case 'weather': return kv([['Observed', d.time + 'Z'], ['Temperature [°C]', d.T_C], ['Pressure [hPa]', d.p_hPa], ['Wind [m/s]', d.wind_ms], ['Gust [m/s]', d.gust_ms], ['Freezing level [m]', d.freezing_level_m], ['Levels aloft', d.aloft.length]]);
      case 'climate': return kv([['Period', `${d.from} to ${d.to}`], ['Hot day, 99th pct [°C]', d.hot99_C], ['Cold day, 1st pct [°C]', d.cold01_C], ['Strong wind, 99th pct [m/s]', d.wind99_ms]]);
      case 'airquality': return kv([['PM10 [µg/m³]', d.pm10], ['PM2.5 [µg/m³]', d.pm2_5], ['Dust [µg/m³]', d.dust], ['Aerosol optical depth', d.aod]]);
      case 'marine': return kv([['Wave height [m]', d.wave_height_m ?? 'inland / none'], ['Wave period [s]', d.wave_period_s]]);
      case 'aerodromes': return kv([['Aerodromes within 40 km', d.fields.length], ['Runways', d.nRunways], ['Nearest', d.fields[0]?.name], ['Longest runway nearby [m]', Math.max(0, ...d.fields.map((f) => f.longest_m || 0)) || null], ['Source', d.source], ['Database date', d.dataset_date]]);
      case 'spaceweather': return kv([['Planetary Kp', d.kp], ['At', d.time], ['GNSS outlook', d.kp >= 7 ? 'severe storm — expect degraded accuracy' : d.kp >= 5 ? 'storm — monitor integrity' : 'quiet']]);
      case 'fx': return kv([['Date', d.date], ['EUR per USD', d.rates.EUR], ['GBP per USD', d.rates.GBP], ['JPY per USD', d.rates.JPY], ['Currencies', Object.keys(d.rates).length]]);
      case 'oil': return kv([['Date', d.date], ['Brent [USD/bbl]', d.brent_usd_bbl], ['Jet fuel estimated from Brent [USD/kg]', jetFromBrent(d.brent_usd_bbl, state.settings.fuelCrack)], ['Estimate is used', 'only when the quoted jet-fuel series is missing or old'], ['Annualised volatility', d.vol_annual]]);
      case 'jetfuel': return h('div', { class: 'stack' }, kv([['Price date', d.date], ['US Gulf Coast jet fuel, spot [USD/gal]', d.usd_gal], ['Converted [USD/kg]', d.usd_kg], ['Density used [kg/L]', d.density_kg_l], ['Annualised volatility', d.vol_annual], ['Series', d.series], ['Last scheduled fetch', d._error ? `failed (${d._error}); showing the previous value` : null]]), inCase('Fuel price', state.case.econ.fuel_usd_kg, 'USD/kg', state.case.econ.fuel_source));
      case 'carbon': { const cp = carbonPrice(d, fxNow); return h('div', { class: 'stack' }, kv([['Auction date', d.date], ['EU allowance auction price [EUR/t CO₂]', d.eur_t], ['Converted [USD/t CO₂]', cp?.usd_t ?? d.usd_t], ['Exchange rate [EUR per USD]', fxNow?.rates?.EUR ?? d.eur_per_usd], ['Auction', d.auction], ['Cover ratio', d.cover_ratio], ['Last scheduled fetch', d._error ? `failed (${d._error}); showing the previous value` : null]]), inCase('Carbon price', state.case.econ.carbon_usd_t, 'USD/t CO₂', state.case.econ.carbon_source)); }
      case 'rates': return kv([['US effective federal funds rate [%]', d.us_fed_funds_pct], ['as of', d.us_date], ['ECB deposit facility rate [%]', d.ecb_deposit_pct], ['as of ', d.ecb_date]]);
      case 'macro': return kv([['Economy', d.country], [`Inflation ${d.inflation_year || ''} [%]`, d.inflation_pct], [`Lending rate ${d.lending_year || ''} [%]`, d.lending_pct], [`GDP growth ${d.gdp_year || ''} [%]`, d.gdp_growth_pct]]);
      case 'grid': return kv([['Carbon intensity [gCO₂/kWh]', d.gCO2_kWh], ['Level', d.index], ['Period from', d.from]]);
      default: return h('span', { class: 'muted small' }, `${d.items?.length || 0} results`);
    }
  };
  const paintList = async () => {
    const cards = []; fxNow = (await idb.get('live.fx.usd'))?.data || null;
    for (const [id, c] of Object.entries(CONNECTORS)) {
      const st = status[id] || {}, p = paramsFor(id), perSuite = id === 'literature' || id === 'opensource';
      const rec = p ? await idb.get(`live.${id}.${c.key(p)}`) : null, age = rec ? Date.now() - rec.ts : null, fresh = age != null && age < c.ttl;
      cards.push(card(c.title, h('div', { class: 'stack' },
        h('div', { class: 'row small' }, perSuite ? badge('Per suite', 'accent') : !p ? badge('Needs a location', 'warn') : c.bundled && rec ? badge(rec.data?.osm ? 'Stored in the app + OpenStreetMap' : 'Stored in the app — works offline', 'ok') : fresh ? badge('Fresh', 'ok') : rec ? badge(c.cloud ? 'Older than expected — snapshot not updated' : 'Stale — refreshing when online', 'warn') : st.ok === false ? badge(st.error || 'Unavailable', 'bad') : badge('Not fetched yet'), c.cloud ? badge('Cloud snapshot only', 'accent') : rec?.via === 'snapshot' ? badge('From cloud snapshot') : null,
          h('span', { class: 'muted' }, rec ? `${c.bundled ? 'Looked up' : c.cloud || rec.via === 'snapshot' ? 'Fetched by the snapshot job' : 'Updated'} ${ago(rec.ts)}` : '', c.bundled ? '' : c.cloud ? ' · snapshot rebuilt about every 3 h' : ` · refreshes every ${fmtTtl(c.ttl)}`)),
        h('p', { class: 'muted small' }, c.use), perSuite ? h('p', { class: 'small' }, 'Open the “Live resources” tab inside any suite.') : preview(id, rec?.data),
        h('div', { class: 'row small' }, h('a', { href: c.home, target: '_blank', rel: 'noopener noreferrer' }, c.provider, ' ', icon('external', 13)),
          p && !perSuite ? btn('Refresh', async () => { const r = await live(id, p, { force: true }); toast(r.error ? `${c.title}: ${r.error}` : `${c.title} updated.`, r.error ? 'bad' : 'ok'); paintList(); }, { ic: 'refresh', kind: 'sm ghost' }) : !p && !perSuite ? btn('Choose location', () => (location.hash = '#/case/site'), { kind: 'sm ghost', ic: 'pin' }) : null))));
    }
    clear(list); add(list, cards);
    const snap = await loadSnapshot().catch(() => null); clear(snapHost);
    if (snap) { const fs = Object.entries(snap.feeds || {}), bad = fs.filter(([, f]) => !f.ok), age = Date.now() - Date.parse(snap.generated); snapHost.append(h('div', { class: 'stack' }, h('div', { class: 'row small' }, badge(age < 8 * 3600e3 ? 'Current' : 'Not updated recently', age < 8 * 3600e3 ? 'ok' : 'warn'), h('span', { class: 'muted' }, `Built ${ago(Date.parse(snap.generated))} · ${fs.length - bad.length} of ${fs.length} feeds fetched${snap.from ? ` · read from ${new URL(snap.from).host}` : ' · saved copy on this device'}`)),
      bad.length ? h('p', { class: 'small' }, `Not fetched in the last run: ${bad.map(([id, f]) => `${CONNECTORS[id]?.title || id} (${f.error})`).join('; ')}. The previous values are kept.`) : null,
      age >= 8 * 3600e3 ? h('p', { class: 'small muted' }, 'The snapshot is rebuilt only by the host that runs the schedule. A copy of the app on another host, or a host whose schedule is switched off, carries an older file; feeds your browser can reach directly are still refreshed.') : null)); }
    else snapHost.append(h('p', { class: 'muted small' }, 'No cloud snapshot could be read (offline at first start, or this copy of the app was opened from a file). Feeds your browser can reach directly still work; carbon price, jet-fuel spot price and policy rates need the snapshot, so the case keeps its entered values for those.'));
    const oil = (await idb.get('live.oil.brent'))?.data;
    plots.splice(0).forEach((x) => x.destroy()); clear(oilHost);
    if (oil?.history) plots.push(renderPlot(oilHost, { type: 'line', title: 'Brent crude, last 12 months', xlabel: 'Trading weeks ago', ylabel: 'Price [USD/bbl]', series: [{ name: 'Brent', x: oil.history.map((_, i) => i - oil.history.length + 1), y: oil.history.map((r) => r.v) }] }));
  };
  host.append(
    h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Always-fresh inputs'), h('h1', null, 'Live data'), h('p', null, 'Three sources keep the inputs current without depending on any one computer. Airports and runways worldwide are stored inside the app. Weather, exchange rates, oil, space weather, economic indicators, research and code indexes are fetched by your browser straight from their public providers. A snapshot of the location-independent feeds — including a carbon price, a quoted jet-fuel price and policy rates, whose providers do not accept browser requests — is rebuilt on the hosting service about every three hours. Everything is saved on your device for offline use.')),
      h('div', { class: 'row' }, btn('Refresh everything', async () => { toast('Refreshing live data…', 'info'); await seedFromSnapshot({ force: true }).catch(() => 0); await Promise.all([refreshSite({ force: true }), refreshGlobal({ force: true })]); toast('Live data refreshed.', 'ok'); paintList(); }, { ic: 'refresh', kind: 'primary' }))),
    card('Settings', h('div', { class: 'grid g3' },
      h('label', { class: 'switch' }, h('input', { type: 'checkbox', checked: state.settings.autoRefreshLive, onchange: (e) => setSetting('autoRefreshLive', e.target.checked) }), 'Refresh automatically while the app is open'),
      h('label', { class: 'switch' }, h('input', { type: 'checkbox', checked: state.settings.autoApplyLive, onchange: (e) => setSetting('autoApplyLive', e.target.checked) }), 'Write live values into the case'),
      h('label', { class: 'row small' }, 'Jet-fuel refining margin over Brent', h('input', { class: 'inp', type: 'number', style: { width: '80px' }, min: 0, max: 200, step: 1, value: Math.round(state.settings.fuelCrack * 100), onchange: (e) => { setSetting('fuelCrack', Math.max(0, Number(e.target.value) / 100)); refreshGlobal().then(paintList); } }), '%'))),
    h('div', { class: 'gap' }), card('Cloud snapshot', snapHost), h('div', { class: 'gap' }), list, h('div', { class: 'gap' }), oilHost,
    h('div', { class: 'gap' }),
    h('div', { class: 'note' }, icon('info'), h('div', null, h('b', null, 'What “live” means here. '), 'A feed is as current as its provider publishes: surface weather every 15 minutes, space weather every 3 hours, exchange rates and central-bank rates each working day, oil and jet-fuel prices each trading day (published a few days in arrears), EU allowance auctions on most working days, World Bank indicators yearly. Carbon and jet-fuel values are public reference prices, not quotations for a particular contract or airport.')),
    h('div', { class: 'gap' }),
    h('div', { class: 'note' }, icon('info'), h('div', null, h('b', null, 'While the app is closed. '), 'A web app cannot run on your device while it is fully closed. Two things cover the gap. First, the cloud snapshot is rebuilt on the hosting service about every 3 hours whether or not anyone has the app open (scheduled jobs can start late, and stop if the host disables the schedule), so the first screen after opening shows location-independent data that is normally no more than a few hours old. Second, browsers with periodic background sync — currently the installed app in Chrome and Edge on Android and desktop — also refresh the saved feeds a few times a day at their own discretion; Safari on iPhone, iPad and Mac, and Firefox, do not offer this. Weather for your site depends on the location, so it is fetched when the app is opened, brought back to the foreground, or reconnects, and every few minutes while it stays open. Offline, the last saved values are used and marked as such. ', bgHost)));
  (async () => { let txt = 'On this device: background refresh while closed is not available; data refreshes when the app is opened.'; try { const reg = await Promise.race([navigator.serviceWorker.ready, new Promise((r) => setTimeout(r, 1500))]); if (reg?.periodicSync && (await reg.periodicSync.getTags()).includes('refresh-live')) txt = 'On this device: background refresh while closed is switched on.'; else if (reg?.periodicSync) txt = 'On this device: the browser supports background refresh but has not granted it (it is offered to installed apps that are used regularly).'; } catch { /* keep the default */ } bgHost.textContent = txt; })();
  await paintList();
  const off = on('live', () => paintList());
  return () => { off(); plots.forEach((p) => p.destroy()); };
}
