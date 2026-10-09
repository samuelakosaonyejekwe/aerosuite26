// Data checks: `node tests/data.mjs [--browser]`
//  - the bundled airport and runway database (js/data/airports, read through js/core/airports.js);
//  - the cloud snapshot file (data/snapshot.json) and the readers in tools/snapshot.mjs;
//  - unit conversions shared by the app and the snapshot tool.
// With --browser it also opens the app in headless Chromium and checks that runways appear on the location
// page with the OpenStreetMap servers blocked, and again with the network switched off after the first load.

import { readFileSync, readdirSync, statSync } from 'node:fs';
import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import * as airports from '../js/core/airports.js';
import { CONNECTORS, usdKgFromUsdGal, usdKgFromUsdBbl, usdFromEur, jetFromBrent, fuelPrice, carbonPrice, searchAirports, JET_DENSITY_KG_L, L_PER_US_GAL } from '../js/core/live.js';
import { parseCsv, surfaceWord } from '../tools/fetch-data.mjs';
import { FEEDS, parseDateValueCsv, parseEiaDailyHtml, parseEexAuctions, unzip } from '../tools/snapshot.mjs';

const root = fileURLToPath(new URL('..', import.meta.url));
let fail = 0, n = 0;
const ok = (cond, msg) => { n++; if (!cond) { fail++; console.log('  FAIL ' + msg); } };
const near = (a, b, tol) => Math.abs(a - b) <= tol;
const section = (t) => console.log(t);

// ---- airports ------------------------------------------------------------------------------------
section('Airport database');
const TOTAL_BUDGET = 3.0e6, TARGET = 2.5e6, TILE_BUDGET = 200e3;
const dir = root + 'js/data/airports/', files = readdirSync(dir), bytes = Object.fromEntries(files.map((f) => [f, statSync(dir + f).size])), total = Object.values(bytes).reduce((a, b) => a + b, 0);
const meta = await airports.meta();
ok(/^\d{4}-\d{2}-\d{2}$/.test(meta.dataset_date), 'index.json records the dataset date');
ok(meta.source === 'OurAirports' && /public domain/i.test(meta.license), 'index.json records source and licence');
ok(total <= TOTAL_BUDGET, `database is ${(total / 1e6).toFixed(2)} MB, over the ${TOTAL_BUDGET / 1e6} MB limit`);
const tileFiles = files.filter((f) => /^t\d+_\d+\.txt$/.test(f)), largest = Math.max(...tileFiles.map((f) => bytes[f]));
ok(largest <= TILE_BUDGET, `largest tile is ${(largest / 1e3).toFixed(0)} kB, over ${TILE_BUDGET / 1e3} kB`);
ok(tileFiles.length === Object.keys(meta.tiles).length, 'every tile in the index has a file and the reverse');
let lines = 0; for (const f of tileFiles) { const k = readFileSync(dir + f, 'utf8').split('\n').filter((l) => l && l[0] !== '#').length; lines += k; if (k !== meta.tiles[f.slice(1, -4)]) ok(false, `tile ${f} holds ${k} aerodromes, index says ${meta.tiles[f.slice(1, -4)]}`); }
ok(lines === meta.counts.airports && lines > 25000, `tiles hold ${lines} aerodromes; index says ${meta.counts.airports}`);
console.log(`  ${lines} aerodromes, ${meta.counts.runways} runways, ${tileFiles.length} tiles, ${(total / 1e6).toFixed(2)} MB (target ${TARGET / 1e6} MB${total > TARGET ? ' — above target, within limit' : ''}), largest tile ${(largest / 1e3).toFixed(0)} kB, data of ${meta.dataset_date}`);

const lagos = await airports.nearest(6.5774, 3.3212, 40, 25), dnmm = lagos.find((a) => a.ident === 'DNMM');
ok(!!dnmm, 'nearest(6.5774, 3.3212) includes DNMM');
if (dnmm) {
  const len = dnmm.runways.map((r) => r.len_m).sort((a, b) => b - a);
  ok(dnmm.runways.length === 2 && near(len[0], 3900, 60) && near(len[1], 2745, 60), `DNMM runways should be about 3.9 km and 2.7 km, got ${len.join(', ')} m`);
  ok(dnmm.icao === 'DNMM' && dnmm.iata === 'LOS' && dnmm.country === 'NG' && dnmm.type === 'large airport', 'DNMM codes, country and type');
  ok(dnmm.dist_km < 2 && near(dnmm.elev_m, 41, 15), 'DNMM distance and elevation');
  const r = dnmm.runways.find((x) => x.ref === '18R/36L');
  ok(!!r && r.surface === 'asphalt' && r.lighted && !r.closed && near(r.heading_deg, 180, 3) && near(r.he.heading_degT, 0, 3) && r.width_m >= 45, 'DNMM 18R/36L surface, lighting, headings, width');
  ok(!!r && r.le.lat != null && near(airports.distanceKm(r.le, r.he) * 1000, r.len_m, 150), 'DNMM 18R/36L thresholds are a runway length apart');
}
ok(lagos.every((a, i) => i === 0 || a.dist_km >= lagos[i - 1].dist_km) && lagos.every((a) => a.dist_km <= 40), 'nearest() is sorted and within the radius');
const egll = await airports.findByCode('EGLL');
ok(egll?.name === 'London Heathrow Airport' && egll.iata === 'LHR' && egll.country === 'GB' && egll.runways.filter((r) => r.len_m > 3500).length === 2, "findByCode('EGLL')");
ok(egll?.runways.some((r) => r.le.displaced_m > 0 || r.he.displaced_m > 0), 'EGLL carries displaced thresholds');
ok((await airports.findByCode('lhr'))?.ident === 'EGLL' && (await airports.findByCode('JFK'))?.ident === 'KJFK', 'findByCode accepts IATA codes in any case');
ok((await airports.findByCode('ZZZZ9')) === null && (await airports.findByCode('')) === null, 'findByCode returns null for unknown codes');
const nbo = await airports.search('Nairobi');
ok(nbo.some((a) => a.ident === 'HKJK') && nbo.some((a) => a.ident === 'HKNW'), "search('Nairobi') finds Jomo Kenyatta and Wilson");
ok((await airports.search('DNMM'))[0]?.ident === 'DNMM' && (await airports.search('los'))[0]?.ident === 'DNMM', 'search() puts an exact code first');
ok((await airports.search('heathrow'))[0]?.ident === 'EGLL' && (await airports.search('sao paulo')).some((a) => a.ident === 'SBGR'), 'search() by name, ignoring accents');
ok((await airports.search('q')).length === 0, 'search() ignores one-letter queries');
ok((await airports.nearest(-17.7, 179.95, 400, 50)).length > 3 && Array.isArray(await airports.nearest(89.9, 10, 300)) && (await airports.nearest(0, -30, 40)).length === 0, 'nearest() across the antimeridian, at the pole and over open ocean');
ok(airports.reciprocal('09L') === '27R' && airports.reciprocal('36') === '18' && airports.reciprocal('18C') === '36C' && airports.reciprocal('H1') === '', 'runway designator reciprocals');
ok(airports.tileOf(6.5, 3.3).id === '9_18' && airports.tileOf(90, 180).id === '17_35' && airports.tileOf(-90, -180).id === '0_0', 'tile indexing at the edges');
for (const [raw, word] of [['ASP', 'asphalt'], ['ASPH-G', 'asphalt'], ['CON', 'concrete'], ['CONC', 'concrete'], ['ASPH/ CONC', 'paved'], ['GRS', 'grass'], ['TURF', 'grass'], ['GVL', 'gravel'], ['GRVL', 'gravel'], ['DIRT', 'dirt'], ['GRE', 'dirt'], ['WATER', 'water'], ['SAND', 'sand'], ['MATS', 'metal'], ['PEM', 'asphalt'], ['', 'unknown'], ['X', 'unknown']]) ok(surfaceWord(raw) === word, `surface "${raw}" should map to ${word}, got ${surfaceWord(raw)}`);
ok(airports.SURFACES.every((w) => /^[a-z]+$/.test(w)), 'surface words are plain words');
ok(JSON.stringify(parseCsv('a,b\n1,"x, ""y"""\n"line\nbreak",2\n')) === JSON.stringify([{ a: '1', b: 'x, "y"' }, { a: 'line\nbreak', b: '2' }]), 'CSV reader handles quotes and line breaks');

// the connector the location page uses: bundled data, no network at all
const noNet = async () => { throw new Error('network must not be used'); };
const con = await CONNECTORS.aerodromes.load({ lat: 6.5774, lon: 3.3212 }, noNet), f0 = con.fields[0];
ok(con.source === 'OurAirports' && con.osm === false && con.nRunways >= 2 && f0.icao === 'DNMM' && f0.kind === 'large airport' && f0.longest_m === f0.runways[0].len_m, 'aerodromes connector answers from the bundled database without the network');
ok(['name', 'icao', 'iata', 'kind', 'ele_m', 'dist_km', 'lat', 'lon', 'runways', 'longest_m'].every((k) => k in f0) && ['ref', 'len_m', 'surface', 'width_m', 'heading_deg', 'lat', 'lon'].every((k) => k in f0.runways[0]), 'aerodromes connector keeps the shape the case page expects');
const hits = await searchAirports('EGLL');
ok(hits[0]?.airport && hits[0].name === 'London Heathrow Airport' && hits[0].country === 'GB' && Number.isFinite(hits[0].lat), 'searchAirports() returns place-search rows');

// ---- conversions ---------------------------------------------------------------------------------
section('Unit conversions');
ok(near(usdKgFromUsdGal(L_PER_US_GAL * JET_DENSITY_KG_L), 1, 1e-12) && near(usdKgFromUsdGal(4.342), 1.42667, 1e-4), 'USD/gal → USD/kg at 0.804 kg/L (4.342 USD/gal = 1.4267 USD/kg)');
ok(near(usdKgFromUsdBbl(158.987 * 0.804), 1, 1e-12) && near(jetFromBrent(100, 0.2), 120 / 127.8255, 1e-6), 'USD/bbl → USD/kg and the Brent-derived estimate');
ok(near(usdFromEur(85.07, 0.89397), 95.16, 0.005) && near(usdFromEur(1, 1), 1, 1e-12), 'EUR → USD with a rate quoted as EUR per USD');
const today = new Date().toISOString().slice(0, 10), old = new Date(Date.now() - 90 * 86400e3).toISOString().slice(0, 10);
ok(fuelPrice({ usd_kg: 1.4, usd_gal: 4.26, date: today }, { brent_usd_bbl: 100, date: today }, 0.2).quoted === true, 'quoted jet fuel is preferred over the Brent estimate');
ok(fuelPrice({ usd_kg: 1.4, usd_gal: 4.26, date: old }, { brent_usd_bbl: 100, date: today }, 0.2).quoted === false && fuelPrice(null, { brent_usd_bbl: 100, date: today }, 0.2).quoted === false && fuelPrice(null, null, 0.2) === null, 'Brent estimate is the fall-back when the quoted series is missing or old');
ok(near(carbonPrice({ eur_t: 80, date: today, eur_per_usd: 0.9 }, { rates: { EUR: 0.8 }, date: today }).usd_t, 100, 1e-9) && near(carbonPrice({ eur_t: 81, date: today, eur_per_usd: 0.9 }, null).usd_t, 90, 1e-9) && carbonPrice({ eur_t: 80, date: old, eur_per_usd: 0.9 }, null) === null, 'carbon price uses the newest exchange rate and ignores old auctions');

// ---- snapshot readers ----------------------------------------------------------------------------
section('Snapshot readers');
ok(JSON.stringify(parseDateValueCsv('observation_date,X\n2026-01-02,1.5\n2026-01-05,.\n2026-01-06,\n2026-01-07,2\n')) === JSON.stringify([{ t: '2026-01-02', v: 1.5 }, { t: '2026-01-07', v: 2 }]), 'date,value CSV skips missing observations');
const eia = parseEiaDailyHtml("<tr><td class='B6'>&nbsp;&nbsp;2026 Sep-28 to Oct- 2</td> <td class='B3'>4.343</td> <td class='B3'>4.399</td> <td class='B3'></td> <td class='B3'>4.355</td> <td class='B3'>4.307</td></tr>");
ok(eia.length === 4 && eia[0].t === '2026-09-28' && eia[3].t === '2026-10-02' && eia[3].v === 4.307, 'EIA history table reader (week rows, month roll-over, blank days)');
const auctions = parseEexAuctions([{ D: 'Public' }, { B: 'Date', D: 'Auction Name', E: 'Contract', F: 'Status', G: 'Auction Price €/tCO2', L: 'Auction Volume tCO2', V: 'Cover Ratio', Z: 'Zone' }, { B: 46301, D: 'Auction 4. Period CAP3 EU', E: 'T3PA', F: 'successful', G: 84.31, L: 2791500, V: 1.49, Z: 'EU' }, { B: 46303, D: 'Auction 4. Period CAP3 EU', E: 'T3PA', F: 'successful', G: 85.07, L: 2791500, V: 1.54, Z: 'EU' }, { B: 46302, D: 'x', E: 'T3PA', F: 'cancelled', G: 0, Z: 'EU' }]);
ok(auctions.length === 2 && auctions[0].date === '2026-10-08' && auctions[0].eur_t === 85.07 && auctions[1].date === '2026-10-06', 'EEX auction report reader (newest first, spreadsheet dates, unsuccessful auctions skipped)');
let threw = false; try { parseEexAuctions([{ B: 'Date', G: 'Auction Price $/t' }]); } catch { threw = true; } ok(threw, 'a changed price unit in the auction report is refused, not misread');
threw = false; try { unzip(Buffer.from('not a zip file at all, just text.')); } catch { threw = true; } ok(threw, 'unzip rejects non-archives');

// ---- snapshot file -------------------------------------------------------------------------------
section('Cloud snapshot file');
let snap = null; try { snap = JSON.parse(readFileSync(root + 'data/snapshot.json', 'utf8')); } catch (e) { ok(false, 'data/snapshot.json is missing or not JSON: ' + e.message); }
if (snap) {
  ok(snap.schema === 1 && Number.isFinite(Date.parse(snap.generated)) && snap.interval_h === 3, 'snapshot header (schema, generated, interval)');
  ok(JSON.stringify(Object.keys(snap.feeds)) === JSON.stringify(Object.keys(FEEDS)), 'snapshot holds exactly the feeds the tool defines');
  for (const [id, f] of Object.entries(snap.feeds)) {
    ok(typeof f.ok === 'boolean' && typeof f.source === 'string' && f.source.length > 5 && /^https:\/\//.test(f.url) && /^https:\/\//.test(f.home) && typeof f.terms === 'string', `${id}: ok flag, source name, address and terms`);
    ok(f.ok ? f.data && Number.isFinite(f.ts) && Number.isFinite(Date.parse(f.fetched)) && !f.error : typeof f.error === 'string' && f.error.length > 0, `${id}: ${f.ok ? 'timestamp and data' : 'error text'}`);
    ok(!!CONNECTORS[id] && typeof CONNECTORS[id].key(id === 'macro' ? { country: 'WLD' } : {}) === 'string', `${id}: the app has a connector for it`);
  }
  const d = (id) => snap.feeds[id]?.data;
  if (d('fx')) ok(d('fx').base === 'USD' && d('fx').rates.USD === 1 && d('fx').rates.EUR > 0.3 && d('fx').rates.EUR < 3, 'fx: USD base with a plausible EUR rate');
  if (d('oil')) ok(d('oil').brent_usd_bbl > 5 && d('oil').brent_usd_bbl < 500 && /^\d{4}-\d{2}-\d{2}$/.test(d('oil').date), 'oil: Brent price and date');
  if (d('spaceweather')) ok(d('spaceweather').kp >= 0 && d('spaceweather').kp <= 9, 'space weather: Kp in range');
  if (d('macro')) ok(d('macro').country === 'World', 'macro: world aggregate');
  if (d('jetfuel')) ok(near(d('jetfuel').usd_kg, usdKgFromUsdGal(d('jetfuel').usd_gal), 1e-4) && d('jetfuel').density_kg_l === JET_DENSITY_KG_L && /Gulf Coast/.test(d('jetfuel').instrument) && d('jetfuel').history.length > 10, 'jet fuel: USD/kg is consistent with USD/gal and the stated density');
  if (d('carbon')) ok(d('carbon').eur_t > 1 && d('carbon').eur_t < 1000 && near(d('carbon').usd_t, usdFromEur(d('carbon').eur_t, d('carbon').eur_per_usd), 0.006) && /EUA/.test(d('carbon').instrument) && /^\d{4}-\d{2}-\d{2}$/.test(d('carbon').date), 'carbon: USD/t is consistent with EUR/t and the stated exchange rate');
  if (d('rates')) ok([d('rates').us_fed_funds_pct, d('rates').ecb_deposit_pct].some((v) => Number.isFinite(v) && v > -2 && v < 30), 'rates: at least one policy rate');
  console.log(`  generated ${snap.generated}, ${Object.values(snap.feeds).filter((f) => f.ok).length}/${Object.keys(snap.feeds).length} feeds ok`);
}

// ---- browser -------------------------------------------------------------------------------------
if (process.argv.includes('--browser')) {
  section('Browser: runways on the location page without OpenStreetMap, then fully offline');
  const { chromium } = await import('playwright-core');
  const port = 20000 + Math.floor(Math.random() * 20000), base = `http://localhost:${port}/`;
  const server = spawn(process.execPath, [root + 'tools/serve.mjs', String(port)], { stdio: 'ignore' });
  for (let k = 0; k < 40; k++) { try { if ((await fetch(base + 'version.json')).ok) break; } catch { /* not listening yet */ } await new Promise((r) => setTimeout(r, 150)); }
  const browser = await chromium.launch();
  const until = async (fn, ms = 40000) => { const t = Date.now(); for (;;) { const v = await fn(); if (v) return v; if (Date.now() - t > ms) return v; await new Promise((r) => setTimeout(r, 250)); } };
  let step = 'start';
  try {
    const ctx = await browser.newContext({ viewport: { width: 1280, height: 900 } }), page = await ctx.newPage(), errors = [];
    page.on('pageerror', (e) => errors.push(e.message));
    let overpass = 0, snapshotReads = 0;
    await page.route(/overpass|maps\.mail\.ru/, (r) => { overpass++; r.abort(); });
    page.on('request', (r) => { if (/data\/snapshot\.json/.test(r.url())) snapshotReads++; });
    await page.goto(base + '#/case/site'); await page.waitForSelector('#page input[type=search]', { timeout: 30000 });
    step = 'search'; await page.locator('#page input[type=search]').fill('Lagos'); await page.locator('.btn.primary', { hasText: 'Search' }).click();
    const nHits = await until(() => page.locator('.hit-place, .hit-airport').count());
    ok(nHits > 0, 'searching "Lagos" lists results');
    const t0 = Date.now(), nAir = await until(() => page.locator('.hit-airport', { hasText: 'Murtala Muhammed' }).count(), 30000);
    ok(nAir === 1, `the airport serving Lagos is merged into the place results (found ${nAir} after ${Date.now() - t0} ms; list: ${(await page.locator('.card').first().innerText()).replace(/\s+/g, ' ').slice(0, 300)})`);
    await until(() => page.locator('.hit-place').count(), 25000); // the online geocoder may take a moment; airports do not wait for it
    const place = page.locator('.hit-place', { hasText: 'Nigeria' }).first(), pick = (await place.count()) ? place : page.locator('.hit-airport', { hasText: 'Murtala Muhammed' }).first();
    console.log(`  selecting: ${(await pick.innerText()).replace(/\s+/g, ' ').slice(0, 90)}`);
    await pick.click();
    const nRw = await until(() => page.locator('.btn.runway').count(), 90000);
    const labels = nRw ? await page.locator('.btn.runway').allInnerTexts() : [];
    ok(nRw >= 2 && labels.some((t) => /18R\/36L/.test(t)) && labels.some((t) => /18L\/36R/.test(t)), `runway buttons for Lagos appear (got: ${labels.join(' | ') || 'none'})`);
    ok(overpass === 0, `no request was made to an Overpass server (${overpass} attempted)`);
    ok(snapshotReads > 0, 'the app read data/snapshot.json at start-up');
    console.log(`  runways shown with Overpass blocked: ${labels.map((t) => t.replace(/\s+/g, ' ').trim()).join(' | ')}`);
    await page.locator('.btn.runway', { hasText: '18R/36L' }).click();
    ok(await until(async () => (await page.evaluate(() => JSON.parse(localStorage.getItem('aerosuite26.case')).site.runway_len_m)) === 3900, 8000), 'choosing a runway writes its length into the case');
    const site = await page.evaluate(() => JSON.parse(localStorage.getItem('aerosuite26.case')).site), econ = await page.evaluate(() => JSON.parse(localStorage.getItem('aerosuite26.case')).econ);
    ok(site.runway_surface === 'asphalt' && site.runway_mu === 0.03 && near(site.runway_heading_deg, 180, 3), 'surface, rolling friction and heading follow the runway');
    ok(await until(async () => { const e = await page.evaluate(() => JSON.parse(localStorage.getItem('aerosuite26.case')).econ); return /EUA/.test(e.carbon_source || '') && /jet fuel spot|Brent/.test(e.fuel_source || ''); }, 20000), `carbon and fuel prices reach the case from the feeds (carbon: ${econ.carbon_source || 'none yet'})`);
    await page.locator('#page input[type=search]').fill('EGLL'); await page.locator('.btn.primary', { hasText: 'Search' }).click();
    { const found = await until(() => page.locator('.hit-airport', { hasText: 'London Heathrow' }).count(), 60000); ok(found, `an ICAO code typed into the place search finds the airport (input "${await page.locator('#page input[type=search]').inputValue()}"; list: ${(await page.locator('.card').first().innerText()).replace(/\s+/g, ' ').slice(0, 300)})`); }

    // live-data page rows for the snapshot-only feeds
    step = 'live page'; await page.goto(base + '#/live'); await page.waitForSelector('#page .card', { timeout: 20000 });
    ok(await until(() => page.locator('.card', { hasText: 'Carbon price (EU ETS allowance)' }).locator('dd').count(), 20000), 'Live data page shows the carbon price row');
    ok(await until(() => page.locator('.card', { hasText: 'Jet-fuel spot price' }).locator('dd').count(), 20000), 'Live data page shows the jet-fuel row');
    ok((await page.locator('.card', { hasText: 'Cloud snapshot' }).first().innerText()).includes('feeds fetched'), 'Live data page shows the snapshot status');

    // fully offline after the first load: the service worker must have stored the airport tiles
    step = 'offline';
    const want = files.length, counts = () => page.evaluate(async () => { let air = 0, all = 0; for (const c of await caches.keys()) if (c.startsWith('aerosuite-app-')) for (const r of await (await caches.open(c)).keys()) { all++; if (/js\/data\/airports\//.test(r.url)) air++; } return { air, all }; });
    let c = await until(async () => { const v = await counts(); return v.air >= want ? v : null; }, 240000) || (await counts());
    for (let prev = -1; prev !== c.all;) { prev = c.all; await new Promise((r) => setTimeout(r, 2500)); c = await counts(); } // let the last few files land
    ok(c.air >= want, `service worker stored the airport database (${c.air} of ${want} files; ${c.all} files in all)`);
    await ctx.setOffline(true); server.kill();
    await page.goto(base + '#/case/site'); await page.reload(); await page.waitForSelector('#page input[type=search]', { timeout: 60000 });
    await page.locator('input[aria-label=Latitude]').fill('-1.3192'); await page.locator('input[aria-label=Longitude]').fill('36.9278'); await page.locator('.btn', { hasText: 'Set coordinates' }).click();
    const off = await until(async () => ((await page.locator('.aerodrome', { hasText: 'Jomo Kenyatta' }).locator('.btn.runway').count()) ? page.locator('.aerodrome', { hasText: 'Jomo Kenyatta' }).locator('.btn.runway').allInnerTexts() : null), 90000);
    ok(!!off && off.some((t) => /06\/24/.test(t)), `offline: runways for a new location (Nairobi) appear from the stored database (${(off || []).join(' | ') || 'none'})`);
    console.log(`  offline, server stopped — Nairobi runways: ${(off || []).map((t) => t.replace(/\s+/g, ' ').trim()).join(' | ')}`);
    await page.locator('#page input[type=search]').fill('NBO'); await page.locator('.btn.primary', { hasText: 'Search' }).click();
    ok(await until(() => page.locator('.hit-airport', { hasText: 'Jomo Kenyatta' }).count(), 60000), 'offline: airport search by IATA code still works');
    ok(errors.length === 0, `no uncaught exceptions (${errors.join('; ')})`);
    await ctx.close();
  } catch (e) { ok(false, `browser check aborted at step "${step}": ` + e.message.split('\n')[0]); }
  await browser.close(); server.kill();
}

console.log(`${n} checks, ${fail} failed`);
process.exit(fail ? 1 : 0);
