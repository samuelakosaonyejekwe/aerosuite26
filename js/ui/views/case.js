// Case & input portal: choose or import an aircraft, edit every shared quantity once, and tie the
// case to a real place on Earth so live weather, terrain and runway data flow into the analyses.

import { h, icon, clear, add, num, btn, badge, card, toast, pickFiles, downloadText, debounce, ago, toCsv } from '../dom.js';
import { state, setCase, patchCase, patchCaseMany, on, importProject, exportProject, setSetting } from '../../core/store.js';
import { PRESETS, CASE_FIELDS, CASE_SECTIONS, derived, blankCase } from '../../core/case.js';
import { METALS, PLIES, BATTERIES } from '../../data/materials.js';
import { geocode, searchAirports, live, refreshSite, refreshGlobal } from '../../core/live.js';
import { isa } from '../../core/atmosphere.js';

const OPTS = { METALS: Object.keys(METALS), PLIES: Object.keys(PLIES), BATTERIES: Object.keys(BATTERIES) };
// Rolling friction by surface word. The bundled airport database (OurAirports codes such as ASP, CON, GRS, GVL mapped to
// plain words by tools/fetch-data.mjs) and OpenStreetMap both use these words; water, metal, wood and unknown fall to 0.04.
const SURFACE_MU = { asphalt: 0.03, concrete: 0.03, paved: 0.03, grass: 0.07, gravel: 0.05, dirt: 0.06, unpaved: 0.06, sand: 0.1, compacted: 0.05, ice: 0.02, snow: 0.05 };

/** Parse a case file: JSON (project or case), CSV/TSV "section.key,value", or INI-style "key = value". */
export function parseCaseFile(name, text) {
  const t = text.trim(), unknown = [], patch = {};
  if (t.startsWith('{')) { const o = JSON.parse(t); return { project: o.case ? o : { case: o }, unknown }; }
  const byKey = new Map(); for (const [s, k] of CASE_FIELDS) { byKey.set(`${s}.${k}`, [s, k]); if (!byKey.has(k)) byKey.set(k, [s, k]); else if (byKey.get(k)[0] !== s) byKey.set(k, null); }
  let section = '';
  for (const raw of t.split(/\r?\n/)) {
    const line = raw.trim(); if (!line || line.startsWith('#') || line.startsWith('//')) continue;
    const sec = line.match(/^\[(.+)\]$/); if (sec) { section = sec[1].trim(); continue; }
    const m = line.match(/^"?([A-Za-z_][\w.]*)"?\s*[=,;\t:]\s*"?([^",;\t]*)"?/); if (!m) continue;
    const key = m[1], hit = byKey.get(key) || (section && byKey.get(`${section}.${key}`)) || null;
    if (!hit) { if (!/^(parameter|key|name|field)$/i.test(key)) unknown.push(key); continue; }
    const v = m[2].trim(), n = Number(v);
    (patch[hit[0]] ||= {})[hit[1]] = v !== '' && Number.isFinite(n) ? n : v;
  }
  if (!Object.keys(patch).length) throw new Error('No recognised case fields were found in this file.');
  return { patch, unknown };
}
const caseToCsv = (c) => toCsv(['parameter', 'value', 'unit', 'description'], CASE_FIELDS.map(([s, k, label, unit]) => [`${s}.${k}`, c[s][k], unit, label]));

export async function render(root, [focus], { setCrumb }) {
  setCrumb('Case & input portal');
  const host = h('div'); root.append(host);
  let tab = focus === 'site' ? 'site' : focus === 'data' ? 'data' : 'aircraft';
  const body = h('div'), tabs = h('div', { class: 'tabs' });
  host.append(h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Step 1'), h('h1', null, 'Case & input portal'), h('p', null, 'Everything the 26 suites share is defined here once: the aircraft, its masses, aerodynamics, propulsion, structure, systems, mission, economics — and the place it operates. Each suite then takes what it needs, and lets you override values locally.'))), tabs, body);
  const paintTabs = () => { clear(tabs); add(tabs, [['aircraft', 'Aircraft data', 'sliders'], ['site', 'Location & live conditions', 'globe'], ['data', 'Import & export', 'upload']].map(([id, l, ic]) => h('button', { class: `tab ${tab === id ? 'on' : ''}`, onclick: () => { tab = id; paintTabs(); paint(); } }, icon(ic, 16), l))); };
  const paint = () => { clear(body); ({ aircraft: tabAircraft, site: tabSite, data: tabData })[tab](); };

  // ---------- aircraft ----------
  function tabAircraft() {
    const c = state.case, d = derived(c), checks = [];
    const sum = c.mass.oew_kg + c.mass.payload_kg + c.mass.fuel_kg;
    if (c.mass.mtow_kg && Math.abs(sum - c.mass.mtow_kg) / c.mass.mtow_kg > 0.12) checks.push(['warn', `Empty mass + payload + fuel = ${num(sum)} kg, which differs from the take-off mass of ${num(c.mass.mtow_kg)} kg by ${num((100 * (sum - c.mass.mtow_kg)) / c.mass.mtow_kg, 3)}%. That is normal only if payload and fuel are traded against each other.`]);
    if (c.wing.S_m2 > 0 && (d.AR < 3 || d.AR > 25)) checks.push(['warn', `Wing aspect ratio is ${num(d.AR, 3)}; most aircraft lie between 5 and 14. Check span and area.`]);
    if (c.wing.S_m2 > 0 && d.wing_loading > 9000) checks.push(['warn', `Wing loading is ${num(d.wing_loading)} N/m², above any transport aircraft in service.`]);
    if (c.mission.cruise_V_ms && c.mission.cruise_alt_m != null) { const M = c.mission.cruise_V_ms / isa(c.mission.cruise_alt_m).a; if (c.mission.cruise_mach && Math.abs(M - c.mission.cruise_mach) > 0.04) checks.push(['info', `Cruise speed of ${num(c.mission.cruise_V_ms)} m/s at ${num(c.mission.cruise_alt_m)} m is Mach ${num(M, 3)}, but cruise Mach is entered as ${c.mission.cruise_mach}. The suites use the speed.`]); }
    if (c.prop.type === 'electric' && c.mass.fuel_kg > 0) checks.push(['info', 'The powerplant is electric but fuel mass is not zero.']);
    if (!checks.length) checks.push(['ok', 'The case passes the built-in consistency checks.']);
    body.append(h('div', { class: 'stack' },
      card('Start from a representative aircraft', h('div', { class: 'grid g4' }, Object.entries(PRESETS).map(([k, p]) => h('button', { class: `suite-card ${state.preset === k ? 'on' : ''}`, style: { textAlign: 'left', cursor: 'pointer', borderColor: state.preset === k ? 'var(--accent)' : '' }, onclick: () => { setCase(p.data(), k); toast(`Loaded “${p.label}”. Results from the previous case are kept until you re-run.`, 'ok'); paint(); } }, h('div', { class: 'top-r' }, h('span', { class: 'num' }, icon(p.icon, 20)), h('h3', null, p.label)), state.preset === k ? badge('Current', 'accent') : null))),
        { collapsible: true }),
      h('div', { class: 'grid g2' },
        card('Derived quantities', h('dl', { class: 'kv' },
          c.wing.S_m2 > 0 ? [h('dt', null, 'Aspect ratio'), h('dd', null, num(d.AR, 4)), h('dt', null, 'Mean aerodynamic chord'), h('dd', null, `${num(d.mac)} m`), h('dt', null, 'Root / tip chord'), h('dd', null, `${num(d.c_root)} / ${num(d.c_tip)} m`), h('dt', null, 'Wing loading'), h('dd', null, `${num(d.wing_loading)} N/m²`)] : null,
          d.T_total ? [h('dt', null, 'Static thrust-to-weight'), h('dd', null, num(d.T_total / d.W, 3))] : null, d.P_total ? [h('dt', null, 'Power loading'), h('dd', null, `${num(d.P_total / c.mass.mtow_kg)} W/kg`)] : null,
          c.rotor.R_m > 0 ? [h('dt', null, 'Rotor tip speed'), h('dd', null, `${num(d.v_tip)} m/s`), h('dt', null, 'Solidity'), h('dd', null, num(d.solidity, 3)), h('dt', null, 'Disk loading'), h('dd', null, `${num(d.disk_loading)} N/m²`)] : null,
          h('dt', null, 'Zero-fuel mass'), h('dd', null, `${num(d.zfw_kg)} kg`))),
        card('Consistency checks', h('div', { class: 'stack' }, checks.map(([k, t]) => h('div', { class: `note ${k === 'info' ? '' : k}` }, icon(k === 'ok' ? 'check' : k === 'warn' ? 'warn' : 'info'), h('div', null, t)))))),
      ...Object.entries(CASE_SECTIONS).map(([sec, title], i) => card(title, h('div', { class: 'form' }, CASE_FIELDS.filter((f) => f[0] === sec).map((f) => field(f))), { collapsible: true, open: i < 4 }))));
  }
  const repaintDerived = debounce(() => { if (tab === 'aircraft') { const y = scrollY; paint(); scrollTo(0, y); } }, 900);
  function field([sec, key, label, unit, help, options]) {
    const v = state.case[sec][key], opts = typeof options === 'string' ? OPTS[options] : options, id = `c-${sec}-${key}`;
    const set = (val) => { patchCase(sec, key, val); if (state.preset !== 'custom') { state.preset = 'custom'; } repaintDerived(); };
    const ctl = opts ? h('select', { class: 'inp', id, onchange: (e) => set(e.target.value) }, opts.map((o) => h('option', { value: o, selected: o === v }, o)))
      : typeof v === 'number' || v == null ? h('input', { class: 'inp', id, type: 'number', step: 'any', inputMode: 'decimal', value: v ?? '', onchange: (e) => { const n = Number(e.target.value); if (Number.isFinite(n)) { e.target.classList.remove('invalid'); set(n); } else e.target.classList.add('invalid'); } })
      : h('input', { class: 'inp', id, type: 'text', value: v, onchange: (e) => set(e.target.value) });
    return h('div', { class: 'field' }, h('label', { for: id }, label, unit ? h('span', { class: 'u' }, `[${unit}]`) : null), ctl, help ? h('div', { class: 'help' }, help) : null);
  }

  // ---------- site ----------
  let sr = null, qText = '', drawSearch = () => {}; // place/airport search on screen ({ text, airs, places, geoErr }) and the text being typed: both survive repaints
  function tabSite() {
    const s = state.case.site, results = h('div', { class: 'stack' }), fieldsHost = h('div', { class: 'stack' }), wxHost = h('div', { class: 'stack' });
    const q = h('input', { class: 'inp', type: 'search', placeholder: 'City or place, airport name, ICAO or IATA code…', value: qText, oninput: (e) => { qText = e.target.value; }, 'aria-label': 'Search for a place', onkeydown: (e) => { if (e.key === 'Enter') search(); } });
    const lat = h('input', { class: 'inp', type: 'number', step: 'any', value: s.lat ?? '', placeholder: 'Latitude', 'aria-label': 'Latitude' }), lon = h('input', { class: 'inp', type: 'number', step: 'any', value: s.lon ?? '', placeholder: 'Longitude', 'aria-label': 'Longitude' });
    const choose = async (p) => {
      patchCaseMany('site', { name: p.name + (p.admin ? `, ${p.admin}` : ''), lat: p.lat, lon: p.lon, elev_m: p.elev_m ?? s.elev_m, country: p.country || '' });
      sr = null; qText = ''; q.value = ''; clear(results); toast(`Location set to ${p.name}. Fetching live conditions…`, 'info');
      const r = await refreshSite({ force: true }); refreshGlobal({ force: true }).catch(() => {});
      if (r?.weather?.data?.elev_m != null && p.elev_m == null) patchCaseMany('site', { elev_m: r.weather.data.elev_m });
      if (r?.weather?.error) toast(`Live weather is unavailable (${r.weather.error}); using saved or standard values.`, 'bad');
      paint();
    };
    // The search lives outside this function's lifetime (see `sr` above): the tab is repainted whenever live data
    // changes the case, and an open result list must survive that.
    const draw = () => {
      clear(results); if (!sr) return;
      const a = sr.airs || [], pl = sr.places || [], list = /^[A-Za-z0-9]{3,4}$/.test(sr.text) && a.length ? [...a, ...pl] : [...pl, ...a];
      if (sr.geoErr) results.append(h('div', { class: 'note warn' }, icon('info'), h('div', null, navigator.onLine === false ? `You are offline, so place names cannot be looked up. ${a.length ? 'Airports are found from the database stored in the app.' : 'Search by airport name or ICAO/IATA code, or enter latitude and longitude.'}` : `Place search failed: ${sr.geoErr.message || 'no reply'}${a.length ? '. Airports below come from the database stored in the app.' : ''}`)));
      else if (sr.airs && sr.places && !list.length) results.append(h('p', { class: 'muted' }, 'No place or airport found. Try another spelling, an ICAO/IATA code, or enter coordinates.'));
      add(results, list.map((p) => h('button', { class: `btn ${p.airport ? 'hit-airport' : 'hit-place'}`, style: { justifyContent: 'flex-start' }, onclick: () => choose(p) }, icon(p.airport ? 'jet' : 'pin', 16), h('span', null, h('b', null, p.name), ` — ${p.admin}${p.airport ? ` · ${p.kind}${p.longest_m ? `, longest runway ${num(p.longest_m)} m` : ''}` : ''} · ${num(p.lat, 5)}°, ${num(p.lon, 5)}°${p.elev_m != null ? ` · ${num(p.elev_m)} m` : ''}`))));
      if (!sr.airs || !sr.places) results.append(h('div', { class: 'row muted' }, h('i', { class: 'spin' }), sr.places ? 'Searching airports…' : 'Searching places…'));
    };
    drawSearch = draw;
    const search = () => {
      const text = (qText = q.value).trim(); if (!text) return;
      // airports come from the database stored in the app (instant, works offline); places from the online geocoder
      const cur = (sr = { text, airs: null, places: null, geoErr: null }), upd = () => { if (sr === cur) drawSearch(); };
      upd();
      searchAirports(text).catch(() => []).then((a) => { cur.airs = a; upd(); });
      geocode(text).then((p) => { cur.places = p; }, (e) => { cur.geoErr = e; cur.places = []; }).then(upd);
    };
    const mine = () => { if (!navigator.geolocation) { toast('This device does not offer location.', 'bad'); return; } navigator.geolocation.getCurrentPosition((pos) => choose({ name: 'My location', admin: '', lat: pos.coords.latitude, lon: pos.coords.longitude, elev_m: pos.coords.altitude ?? null }), () => toast('Location permission was not granted.', 'bad'), { timeout: 15000 }); };
    const manual = () => { const a = Number(lat.value), b = Number(lon.value); if (!(Math.abs(a) <= 90) || !(Math.abs(b) <= 180) || lat.value === '' || lon.value === '') { toast('Enter a latitude between −90 and 90 and a longitude between −180 and 180.', 'bad'); return; } choose({ name: `${a.toFixed(4)}°, ${b.toFixed(4)}°`, admin: '', lat: a, lon: b, elev_m: null }); };

    // weather summary + wind/runway diagram
    const std = isa(s.elev_m || 0), dISA = s.T_C + 273.15 - std.T, rho = (s.p_hPa * 100) / (287.05 * (s.T_C + 273.15));
    let da = 0; for (let lo = -1000, hi = 12000, k = 0; k < 50; k++) { da = (lo + hi) / 2; if (isa(da).rho > rho) lo = da; else hi = da; }
    const rel = s.runway_heading_deg == null ? null : ((s.wind_dir_deg - s.runway_heading_deg) * Math.PI) / 180, hw = rel == null ? null : s.wind_ms * Math.cos(rel), xw = rel == null ? null : s.wind_ms * Math.sin(rel);
    const best = hw == null ? null : Math.abs(hw); // use the runway end that faces the wind
    wxHost.append(h('div', { class: 'grid g2' },
      h('div', { class: 'stack' },
        h('dl', { class: 'kv' }, h('dt', null, 'Place'), h('dd', null, s.name), h('dt', null, 'Coordinates'), h('dd', null, s.lat == null ? 'not set' : `${num(s.lat, 6)}°, ${num(s.lon, 6)}°`), h('dt', null, 'Elevation'), h('dd', null, `${num(s.elev_m)} m`),
          h('dt', null, 'Temperature'), h('dd', null, `${num(s.T_C, 3)} °C (ISA ${dISA >= 0 ? '+' : ''}${num(dISA, 3)} K)`), h('dt', null, 'Station pressure'), h('dd', null, `${num(s.p_hPa, 5)} hPa`), h('dt', null, 'Relative humidity'), h('dd', null, `${num(100 * s.rh, 3)} %`),
          h('dt', null, 'Density altitude'), h('dd', null, h('b', null, `${num(da, 3)} m`)), h('dt', null, 'Wind'), h('dd', null, `${num(s.wind_ms, 3)} m/s from ${num(s.wind_dir_deg, 3)}°, gusting ${num(s.gust_ms, 3)} m/s`),
          h('dt', null, 'Precipitation'), h('dd', null, `${num(s.precip_mm_h, 2)} mm/h · runway ${s.runway_state || 'dry'}`), h('dt', null, 'Freezing level'), h('dd', null, `${num(s.freezing_level_m, 3)} m`), h('dt', null, 'Hot-day design temperature'), h('dd', null, `${num(s.design_hot_C, 3)} °C (99th percentile, last 12 months)`),
          h('dt', null, 'Sea state'), h('dd', null, s.wave_height_m ? `${num(s.wave_height_m, 2)} m significant wave height` : 'not applicable / calm'), h('dt', null, 'Geomagnetic Kp'), h('dd', null, num(s.kp_index, 2)), h('dt', null, 'Source'), h('dd', null, s.source, s.updated ? ` · ${ago(s.updated)}` : '')),
        h('div', { class: 'row' }, btn('Refresh now', async () => { const r = await refreshSite({ force: true }); toast(r?.weather?.error ? `Could not refresh: ${r.weather.error}` : r ? 'Live conditions updated.' : 'Choose a location first.', r?.weather?.error ? 'bad' : 'ok'); paint(); }, { ic: 'refresh' }),
          h('label', { class: 'switch' }, h('input', { type: 'checkbox', checked: state.settings.autoApplyLive, onchange: (e) => setSetting('autoApplyLive', e.target.checked) }), 'Write live values into the case automatically'))),
      h('div', { class: 'stack' }, windRose(s),
        rel == null ? h('p', { class: 'muted small' }, 'Pick a runway below to see head- and crosswind components.') : h('div', { class: 'kpis' }, h('div', { class: 'kpi' }, h('span', { class: 'l' }, 'Headwind on the into-wind runway end'), h('div', { class: 'v' }, num(best, 3), h('small', null, 'm/s'))), h('div', { class: `kpi ${Math.abs(xw) > 10 ? 'warn' : 'ok'}` }, h('span', { class: 'l' }, 'Crosswind component'), h('div', { class: 'v' }, num(Math.abs(xw), 3), h('small', null, 'm/s')))))));
    if (s.winds_aloft?.length) wxHost.append(h('details', null, h('summary', { class: 'small', style: { cursor: 'pointer' } }, 'Winds and temperatures aloft'), h('div', { class: 'table-wrap', style: { marginTop: '8px' } }, h('table', { class: 'data' }, h('thead', null, h('tr', null, ['Pressure [hPa]', 'Height [m]', 'Temperature [°C]', 'ISA deviation [K]', 'Wind [m/s]', 'From [°]'].map((x) => h('th', { class: 'num' }, x)))), h('tbody', null, s.winds_aloft.map((a) => h('tr', null, [a.hPa, a.alt_m, a.T_C, a.T_C + 273.15 - isa(a.alt_m).T, a.speed_ms, a.dir_deg].map((v) => h('td', { class: 'num' }, num(v, 4))))))))));

    const applyRunway = (f, rw) => {
      const mu = SURFACE_MU[String(rw.surface).split(/[:;_]/)[0]] ?? 0.04, a = rw.le?.elev_m, b = rw.he?.elev_m, slope = a != null && b != null && rw.len_m > 0 ? Math.round((10000 * (b - a)) / rw.len_m) / 100 : null; // uphill + in the direction of the stated heading
      patchCaseMany('site', { runway_len_m: rw.len_m, runway_heading_deg: rw.heading_deg, runway_surface: rw.surface, runway_mu: mu, ...(slope != null && Math.abs(slope) < 15 ? { runway_slope_pct: slope } : {}), ...(f.ele_m != null ? { elev_m: f.ele_m } : {}), name: `${f.name}${f.icao ? ` (${f.icao})` : ''}` });
      toast(`Runway ${rw.ref || ''} at ${f.name} applied: ${rw.len_m} m, ${rw.surface}${rw.heading_source === 'designator' ? '. Heading taken from the runway number (magnetic, nearest 10°)' : ''}.`, 'ok'); paint();
    };
    const runwayTip = (rw) => [rw.heading_deg != null ? `Heading ${rw.heading_deg}°${rw.heading_source === 'designator' ? ' (from the runway number)' : rw.heading_source === 'surveyed' ? ' true' : ''}` : null, rw.le?.elev_m != null && rw.he?.elev_m != null ? `Threshold elevations ${rw.le.elev_m} m (${rw.le.ident}) and ${rw.he.elev_m} m (${rw.he.ident})` : null, rw.le?.displaced_m ? `${rw.le.ident} threshold displaced ${rw.le.displaced_m} m` : null, rw.he?.displaced_m ? `${rw.he.ident} threshold displaced ${rw.he.displaced_m} m` : null, rw.lighted ? 'Lighted' : null].filter(Boolean).join(' · ');
    const loadFields = async (force, osm = false) => {
      if (s.lat == null) return; clear(fieldsHost); fieldsHost.append(h('div', { class: 'row muted' }, h('i', { class: 'spin' }), osm ? 'Asking OpenStreetMap for additional detail…' : 'Looking up aerodromes within 40 km…'));
      const r = await live('aerodromes', { lat: s.lat, lon: s.lon, ...(osm ? { osm: true } : {}) }, { force: force || osm }); clear(fieldsHost);
      const osmBtn = btn('Add detail from OpenStreetMap', () => loadFields(true, true), { ic: 'globe', kind: 'sm ghost', title: 'Optional: asks the public OpenStreetMap servers for heliports, small strips and missing runway details near this place' });
      if (!r.data) { fieldsHost.append(h('div', { class: 'note warn' }, icon('info'), h('div', null, `Aerodrome data is unavailable (${r.error}). Enter runway length, heading and surface in the fields below.`))); return; }
      if (osm && r.data.osm_error) fieldsHost.append(h('div', { class: 'note warn' }, icon('info'), h('div', null, `OpenStreetMap did not answer (${r.data.osm_error}). Showing the airport database stored in the app.`)));
      if (!r.data.fields.length) { fieldsHost.append(h('p', { class: 'muted' }, 'No aerodrome within 40 km in the airport database stored in the app. Small private strips and most heliports are not in it: try OpenStreetMap, or enter the runway or landing-site data by hand below.'), h('div', { class: 'row' }, osmBtn)); return; }
      add(fieldsHost, r.data.fields.slice(0, 8).map((f) => h('div', { class: 'stack aerodrome', style: { gap: '4px' } }, h('div', null, h('b', null, f.name), ' ', f.icao ? badge(f.icao, 'accent') : null, f.iata ? badge(f.iata) : null, h('span', { class: 'muted small' }, ` · ${f.kind ? `${f.kind} · ` : ''}${num(f.dist_km, 3)} km away${f.ele_m != null ? ` · ${num(f.ele_m)} m` : ''}${f.closed_runways ? ` · ${f.closed_runways} closed runway${f.closed_runways > 1 ? 's' : ''} not shown` : ''}`)),
        f.runways.length ? h('div', { class: 'row' }, f.runways.slice(0, 6).map((rw) => btn(`${rw.ref || 'Runway'} · ${num(rw.len_m)} m${rw.width_m ? ` × ${num(rw.width_m)} m` : ''} · ${rw.surface}`, () => applyRunway(f, rw), { kind: 'sm runway', ic: 'check', title: runwayTip(rw) || null }))) : h('span', { class: 'muted small' }, 'No runway data for this field.'))));
      fieldsHost.append(h('div', { class: 'row' }, osmBtn),
        h('p', { class: 'muted small' }, `Airports and runways: OurAirports (public domain)${r.data.dataset_date ? `, data of ${r.data.dataset_date}` : ''}, stored in the app and available offline.${r.data.osm ? ' Additional detail © OpenStreetMap contributors (ODbL).' : ''} Community-maintained data: confirm declared distances in the official aeronautical publication before operational use.`));
    };
    const siteField = (key, label, unit, help) => { const id = `s-${key}`; return h('div', { class: 'field' }, h('label', { for: id }, label, unit ? h('span', { class: 'u' }, `[${unit}]`) : null), h('input', { class: 'inp', id, type: 'number', step: 'any', value: state.case.site[key] ?? '', onchange: (e) => { const n = Number(e.target.value); if (Number.isFinite(n)) patchCase('site', key, n); } }), help ? h('div', { class: 'help' }, help) : null); };
    body.append(h('div', { class: 'stack' },
      card('Where does the aircraft operate?', h('div', { class: 'stack' },
        h('div', { class: 'row' }, h('div', { class: 'grow', style: { minWidth: '200px' } }, q), btn('Search', search, { ic: 'search', kind: 'primary' }), btn('Use my location', mine, { ic: 'pin' })),
        h('div', { class: 'row' }, h('div', { style: { width: '150px' } }, lat), h('div', { style: { width: '150px' } }, lon), btn('Set coordinates', manual, { kind: 'ghost' })), results)),
      card('Conditions at the site', wxHost),
      card('Aerodromes and runways nearby', fieldsHost, { actions: s.lat != null ? btn('', () => loadFields(true), { ic: 'refresh', kind: 'sm ghost', title: 'Refresh' }) : null }),
      card('Site and runway values used by the suites', h('div', { class: 'form' }, siteField('elev_m', 'Field elevation', 'm'), siteField('T_C', 'Outside air temperature', '°C'), siteField('p_hPa', 'Station pressure', 'hPa'), siteField('wind_ms', 'Wind speed', 'm/s'), siteField('wind_dir_deg', 'Wind direction (from)', 'deg'), siteField('runway_len_m', 'Runway length available', 'm'), siteField('runway_heading_deg', 'Runway heading', 'deg'), siteField('runway_slope_pct', 'Runway slope (uphill +)', '%'), siteField('runway_mu', 'Rolling friction', '-', '0.02–0.03 paved, 0.05 gravel, 0.07 short grass'), siteField('runway_mu_brake', 'Braking friction', '-', '0.4 dry, 0.25 wet, 0.1 ice or freezing precipitation')), { collapsible: true, open: false })));
    draw();
    if (s.lat != null) loadFields(false); else fieldsHost.append(h('p', { class: 'muted' }, 'Choose a location to list nearby aerodromes with runway length, heading and surface.'));
  }

  // ---------- import & export ----------
  function tabData() {
    const report = h('div', { class: 'stack' });
    const take = async (file) => {
      clear(report);
      try {
        if (file.size > 20e6) throw new Error('Case files are small text files; this one is larger than 20 MB. Geometry and mesh files belong in the Geometry & mesh page.');
        const r = parseCaseFile(file.name, await file.text());
        if (r.project) { importProject(r.project); report.append(h('div', { class: 'note ok' }, icon('check'), h('div', null, `Loaded “${state.case.meta.name || file.name}” from ${file.name}.`))); }
        else { const c = state.case; let n = 0; for (const [sec, obj] of Object.entries(r.patch)) { Object.assign(c[sec], obj); n += Object.keys(obj).length; } setCase(c, 'custom'); report.append(h('div', { class: 'note ok' }, icon('check'), h('div', null, `${n} values from ${file.name} were merged into the current case.`))); }
        if (r.unknown.length) report.append(h('div', { class: 'note warn' }, icon('info'), h('div', null, `${r.unknown.length} unrecognised entr${r.unknown.length > 1 ? 'ies were' : 'y was'} ignored: ${r.unknown.slice(0, 12).join(', ')}${r.unknown.length > 12 ? '…' : ''}. Download the template to see the accepted names.`)));
      } catch (e) { report.append(h('div', { class: 'note bad' }, icon('warn'), h('div', null, `${file.name} could not be read: ${e.message}`))); }
    };
    const drop = h('div', { class: 'drop', tabIndex: 0, role: 'button', onclick: async () => { for (const f of await pickFiles({ accept: '.json,.csv,.tsv,.txt,.ini,.dat,.cfg' })) take(f); }, ondragover: (e) => { e.preventDefault(); drop.classList.add('over'); }, ondragleave: () => drop.classList.remove('over'), ondrop: (e) => { e.preventDefault(); drop.classList.remove('over'); [...e.dataTransfer.files].forEach(take); } },
      icon('upload', 30), h('b', null, 'Drop a case or project file here, or click to choose'), h('span', { class: 'muted small' }, 'JSON project or case · CSV or TSV (parameter, value) · key = value text'));
    body.append(h('div', { class: 'stack' },
      card('Import case data', h('div', { class: 'stack' }, drop, report, h('p', { class: 'muted small' }, 'Tabular files use one row per quantity, for example “mass.mtow_kg, 78000”. Values merge into the current case, so a file can hold only what differs. Suite-specific tables (load spectra, laminate stacks, measured data) are entered inside each suite.'))),
      card('Export', h('div', { class: 'row' }, btn('Project (JSON)', () => downloadText(`${slug(state.case.meta.name)}.aerosuite.json`, JSON.stringify(exportProject(), null, 1), 'application/json'), { ic: 'download', kind: 'primary' }), btn('Case table (CSV)', () => downloadText(`${slug(state.case.meta.name)}-case.csv`, caseToCsv(state.case), 'text/csv'), { ic: 'download' }), btn('Blank template (CSV)', () => downloadText('aerosuite-case-template.csv', caseToCsv(blankCase()), 'text/csv'), { ic: 'doc', kind: 'ghost' })))));
  }
  const slug = (s) => String(s || 'case').toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '') || 'case';

  paintTabs(); paint();
  // Live data repaints the location tab; never in the middle of a click (the button under the pointer would vanish)
  let pressedAt = 0; host.addEventListener('pointerdown', () => { pressedAt = Date.now(); }, true);
  const repaintSite = debounce(() => { if (Date.now() - pressedAt < 700) { repaintSite(); return; } if (tab === 'site' && document.activeElement?.tagName !== 'INPUT') paint(); }, 600);
  const off = on('case', repaintSite);
  return off;
}

/** Compass with wind arrow and runway bar. */
function windRose(s) {
  const R = 78, cx = 100, cy = 100, pt = (deg, r) => [cx + r * Math.sin((deg * Math.PI) / 180), cy - r * Math.cos((deg * Math.PI) / 180)];
  const svg = h('svg', { viewBox: '0 0 200 200', class: 'flowmap', style: { maxWidth: '230px', margin: '0 auto' }, role: 'img', 'aria-label': 'Wind and runway directions' }, h('circle', { cx, cy, r: R, fill: 'none', stroke: 'var(--line)', 'stroke-width': 1.5 }));
  for (const [d, l] of [[0, 'N'], [90, 'E'], [180, 'S'], [270, 'W']]) { const [x, y] = pt(d, R + 12); const t = h('text', { x, y: y + 4, 'text-anchor': 'middle' }); t.textContent = l; svg.append(t); }
  if (s.runway_heading_deg != null) { const [x1, y1] = pt(s.runway_heading_deg, R - 8), [x2, y2] = pt(s.runway_heading_deg + 180, R - 8); svg.append(h('line', { x1, y1, x2, y2, stroke: 'var(--ink-2)', 'stroke-width': 9, 'stroke-linecap': 'butt' })); }
  if (s.wind_ms > 0.2) { const [x1, y1] = pt(s.wind_dir_deg, R - 4), [x2, y2] = pt(s.wind_dir_deg, 14), [a1, b1] = pt(s.wind_dir_deg + 16, 30), [a2, b2] = pt(s.wind_dir_deg - 16, 30); svg.append(h('line', { x1, y1, x2, y2, stroke: 'var(--series-1)', 'stroke-width': 3, 'stroke-linecap': 'round' }), h('polyline', { points: `${a1},${b1} ${x2},${y2} ${a2},${b2}`, fill: 'none', stroke: 'var(--series-1)', 'stroke-width': 3, 'stroke-linecap': 'round', 'stroke-linejoin': 'round' })); }
  return svg;
}
