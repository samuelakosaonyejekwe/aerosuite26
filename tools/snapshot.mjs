// Cloud snapshot of every location-independent live feed: `node tools/snapshot.mjs [--out <file>] [--only a,b]`
// Runs on the hosting side on a schedule (.github/workflows/snapshot.yml) and locally. It fetches, without any
// key or account, the feeds the app also reads from the browser (exchange rates, Brent, space weather, World
// Bank world aggregates) plus the ones a browser cannot reach because the provider sends no CORS headers
// (jet-fuel spot price, EU carbon allowance price, policy interest rates), and writes data/snapshot.json.
// A feed that fails keeps its last good value (from the previous snapshot) and is marked ok: false.
//
// Environment: SNAPSHOT_PREV_URL — address of the previously published snapshot.json, used as the last-good
// source when the working copy has none newer.

import { readFile, writeFile, mkdir } from 'node:fs/promises';
import { inflateRawSync } from 'node:zlib';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { CONNECTORS, usdKgFromUsdGal, usdFromEur, JET_DENSITY_KG_L, L_PER_US_GAL } from '../js/core/live.js';

const root = resolve(fileURLToPath(new URL('..', import.meta.url)));
const UA = 'AeroSuite26-snapshot/1.0 (scheduled public-data snapshot; static web app)';
export const INTERVAL_H = 3;

/** fetch with a timeout and one retry; records every address requested. */
function makeGet(log) {
  return async (url, opt = {}) => {
    log.push(url); let last;
    for (let k = 0; k < 2; k++) {
      try { const r = await fetch(url, { ...opt, headers: { 'user-agent': UA, ...(opt.headers || {}) }, signal: AbortSignal.timeout(45000) }); if (r.status >= 500 && k === 0) { last = new Error(`HTTP ${r.status}`); continue; } return r; }
      catch (e) { last = e; await new Promise((res) => setTimeout(res, 1500)); }
    }
    throw new Error(last?.cause?.code || last?.message || 'request failed');
  };
}
const okText = async (r) => { if (!r.ok) throw new Error(`HTTP ${r.status}`); return r.text(); };

// ---- small readers ------------------------------------------------------------------------------
/** "date,value" CSV (FRED graph export) → [{ t, v }], skipping missing observations. */
export function parseDateValueCsv(text) {
  return text.trim().split(/\r?\n/).slice(1).map((l) => l.split(',')).filter((r) => /^\d{4}-\d{2}-\d{2}$/.test(r[0]) && r[1] !== '' && Number.isFinite(Number(r[1]))).map((r) => ({ t: r[0], v: Number(r[1]) }));
}
/** EIA "history" page for a daily series (weeks as rows, Monday–Friday as columns) → [{ t, v }]. */
export function parseEiaDailyHtml(html) {
  const out = [], MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  for (const m of html.matchAll(/<td class='B6'>(?:&nbsp;)*(\d{4}) (\w{3})-\s*(\d+) to[^<]*<\/td>((?:\s*<td class='B3'>[^<]*<\/td>){5})/g)) {
    const start = Date.UTC(Number(m[1]), MON.indexOf(m[2]), Number(m[3])); if (!Number.isFinite(start)) continue;
    [...m[4].matchAll(/<td class='B3'>([^<]*)<\/td>/g)].forEach((c, i) => { const v = Number(c[1]); if (c[1].trim() !== '' && Number.isFinite(v)) out.push({ t: new Date(start + i * 86400e3).toISOString().slice(0, 10), v }); });
  }
  return out;
}
/** Minimal ZIP reader (stored and deflated entries), enough for .xlsx workbooks. */
export function unzip(buf) {
  let e = buf.length - 22; while (e >= 0 && buf.readUInt32LE(e) !== 0x06054b50) e--;
  if (e < 0) throw new Error('not a ZIP archive');
  const out = new Map(); let p = buf.readUInt32LE(e + 16);
  for (let i = 0, n = buf.readUInt16LE(e + 10); i < n; i++) {
    if (buf.readUInt32LE(p) !== 0x02014b50) throw new Error('damaged ZIP directory');
    const method = buf.readUInt16LE(p + 10), size = buf.readUInt32LE(p + 20), nameLen = buf.readUInt16LE(p + 28), extra = buf.readUInt16LE(p + 30), comment = buf.readUInt16LE(p + 32), at = buf.readUInt32LE(p + 42), name = buf.toString('utf8', p + 46, p + 46 + nameLen);
    const start = at + 30 + buf.readUInt16LE(at + 26) + buf.readUInt16LE(at + 28), raw = buf.subarray(start, start + size);
    out.set(name, method === 0 ? raw : inflateRawSync(raw));
    p += 46 + nameLen + extra + comment;
  }
  return out;
}
const xmlText = (s) => s.replace(/<[^>]+>/g, '').replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"').replace(/&apos;/g, "'").replace(/&#(\d+);/g, (_, d) => String.fromCodePoint(Number(d))).replace(/&amp;/g, '&');
/** First worksheet of an .xlsx file → rows of { column letter: value } (numbers as numbers, strings as text). */
export function readXlsxSheet(buf) {
  const z = unzip(buf), shared = z.has('xl/sharedStrings.xml') ? [...z.get('xl/sharedStrings.xml').toString('utf8').matchAll(/<si>([\s\S]*?)<\/si>/g)].map((m) => xmlText(m[1].replace(/<rPh>[\s\S]*?<\/rPh>/g, ''))) : [];
  const name = [...z.keys()].filter((k) => /^xl\/worksheets\/sheet\d+\.xml$/.test(k)).sort()[0]; if (!name) throw new Error('workbook has no worksheet');
  const rows = [];
  for (const r of z.get(name).toString('utf8').matchAll(/<row [^>]*>([\s\S]*?)<\/row>/g)) {
    const row = {};
    for (const c of r[1].matchAll(/<c r="([A-Z]+)\d+"([^>]*?)(?:\/>|>([\s\S]*?)<\/c>)/g)) {
      const v = /<v>([\s\S]*?)<\/v>/.exec(c[3] || ''), inline = /<is>([\s\S]*?)<\/is>/.exec(c[3] || '');
      if (inline) row[c[1]] = xmlText(inline[1]); else if (v) row[c[1]] = /t="s"/.test(c[2]) ? shared[Number(v[1])] : /t="(str|b|e)"/.test(c[2]) ? xmlText(v[1]) : Number(v[1]);
    }
    rows.push(row);
  }
  return rows;
}
const excelDate = (serial) => new Date(Math.round((Math.floor(serial) - 25569) * 86400e3)).toISOString().slice(0, 10);
/** EEX "Emission Spot Primary Market Auction Report" rows → successful EUA auctions, newest first. */
export function parseEexAuctions(rows) {
  const hi = rows.findIndex((r) => Object.values(r).some((v) => typeof v === 'string' && /^Auction Price/i.test(v))); if (hi < 0) throw new Error('auction table not found in the report');
  const col = (re) => Object.keys(rows[hi]).find((k) => typeof rows[hi][k] === 'string' && re.test(rows[hi][k].replace(/\s+/g, ' ')));
  const C = { date: col(/^Date$/i), name: col(/^Auction Name$/i), contract: col(/^Contract$/i), status: col(/^Status$/i), price: col(/^Auction Price/i), volume: col(/^Auction Volume/i), cover: col(/^Cover Ratio$/i), zone: col(/^Zone$/i) };
  if (!C.date || !C.price || !/€\s*\/\s*tCO2/i.test(rows[hi][C.price])) throw new Error('auction report layout changed (price column is not €/tCO2)');
  return rows.slice(hi + 1).filter((r) => typeof r[C.date] === 'number' && typeof r[C.price] === 'number' && r[C.price] > 0 && (!C.status || /^successful$/i.test(String(r[C.status] || ''))))
    .map((r) => ({ date: excelDate(r[C.date]), eur_t: r[C.price], name: String(r[C.name] || ''), contract: String(r[C.contract] || ''), zone: String(r[C.zone] || ''), volume_t: typeof r[C.volume] === 'number' ? r[C.volume] : null, cover_ratio: typeof r[C.cover] === 'number' ? r[C.cover] : null }))
    .sort((a, b) => (a.date < b.date ? 1 : a.date > b.date ? -1 : 0));
}
const annualVol = (px) => { const r = px.slice(1).map((v, i) => Math.log(v / px[i])).filter(Number.isFinite); if (r.length < 20) return null; const m = r.reduce((a, b) => a + b, 0) / r.length; return Math.sqrt(r.reduce((a, b) => a + (b - m) ** 2, 0) / (r.length - 1)) * Math.sqrt(252); };
const weekly = (rows) => rows.filter((_, i) => i % 5 === 0 || i === rows.length - 1);
const iso = (d) => d.toISOString().slice(0, 10);

// ---- feeds --------------------------------------------------------------------------------------
// Each: { source, home, terms, run(get, done) → data }. `done` holds the feeds finished earlier (for FX).
const viaConnector = (id, params, terms) => ({ source: CONNECTORS[id].provider, home: CONNECTORS[id].home, terms, run: (get) => CONNECTORS[id].load(params, get) });
export const FEEDS = {
  fx: viaConnector('fx', {}, 'ECB euro foreign exchange reference rates, free reuse with attribution; served by the open Frankfurter API'),
  oil: viaConnector('oil', {}, 'EIA data are US Government works in the public domain; dataset repackaged under ODC-PDDL 1.0'),
  spaceweather: viaConnector('spaceweather', {}, 'NOAA SWPC products are US Government works in the public domain'),
  macro: viaConnector('macro', { country: 'WLD' }, 'World Bank Open Data, CC BY 4.0'),
  jetfuel: {
    source: 'US Energy Information Administration — U.S. Gulf Coast Kerosene-Type Jet Fuel Spot Price FOB (series EER_EPJK_PF4_RGC_DPG; FRED series DJFUELUSGULF)', home: 'https://www.eia.gov/dnav/pet/hist/EER_EPJK_PF4_RGC_DPGD.htm',
    terms: 'EIA data are US Government works in the public domain (attribution requested). Read through the FRED graph CSV export, with the EIA history page as fall-back.',
    async run(get) {
      let rows = [], via = 'FRED';
      try { rows = parseDateValueCsv(await get(`https://fred.stlouisfed.org/graph/fredgraph.csv?id=DJFUELUSGULF&cosd=${iso(new Date(Date.now() - 400 * 86400e3))}`).then(okText)); if (!rows.length) throw new Error('empty series'); }
      catch { via = 'EIA'; rows = parseEiaDailyHtml(await get('https://www.eia.gov/dnav/pet/hist/EER_EPJK_PF4_RGC_DPGD.htm').then(okText)); }
      const tail = rows.slice(-260), last = tail[tail.length - 1]; if (!last || !(last.v > 0.2 && last.v < 30)) throw new Error('no plausible jet-fuel observation'); // USD/gal sanity window
      return { date: last.t, usd_gal: last.v, usd_kg: Math.round(usdKgFromUsdGal(last.v) * 1e4) / 1e4, density_kg_l: JET_DENSITY_KG_L, litres_per_us_gal: L_PER_US_GAL, instrument: 'U.S. Gulf Coast kerosene-type jet fuel, spot price FOB, daily', series: via === 'FRED' ? 'FRED DJFUELUSGULF' : 'EIA EER_EPJK_PF4_RGC_DPG', via, vol_annual: annualVol(tail.map((r) => r.v)), history: weekly(tail) };
    },
  },
  carbon: {
    source: 'European Energy Exchange (EEX) — Emission Spot Primary Market Auction Report (EU ETS allowance auctions held for the EU Member States)', home: 'https://www.eex.com/en/markets/environmentals/eu-ets1-eu-ets2-auctions/eu-ets1-auctions',
    terms: 'Auction results are published by the auction platform as required by the EU Auctioning Regulation; the report file is marked "Public". Only the most recent clearing price is quoted here, with attribution and a link to the source; the auction dataset itself is not redistributed. Not a licence for EEX exchange market data.',
    async run(get, done) {
      const y = new Date().getUTCFullYear(), file = (yr) => `https://public.eex-group.com/eex/eua-auction-report/emission-spot-primary-market-auction-report-${yr}-data.xlsx`;
      let all = [], used = null;
      for (const yr of [y, y - 1]) { try { const r = await get(file(yr)); if (!r.ok) throw new Error(`HTTP ${r.status}`); const a = parseEexAuctions(readXlsxSheet(Buffer.from(await r.arrayBuffer()))); if (a.length) { all = all.concat(a); used ||= file(yr); if (all.length >= 60) break; } } catch (e) { if (yr === y - 1 && !all.length) throw e; } }
      const eua = all.filter((a) => /EUA|T3PA/i.test(a.contract) || /CAP3/i.test(a.name) || !a.contract), pool = eua.filter((a) => a.zone === 'EU').length ? eua.filter((a) => a.zone === 'EU') : eua, last = pool[0];
      if (!last || !(last.eur_t > 1 && last.eur_t < 1000)) throw new Error('no plausible auction price');
      const fx = done.fx?.data, rate = fx?.rates?.EUR;
      return { date: last.date, eur_t: last.eur_t, usd_t: Number.isFinite(rate) ? Math.round(usdFromEur(last.eur_t, rate) * 100) / 100 : null, eur_per_usd: rate ?? null, fx_date: fx?.date ?? null,
        instrument: 'EU ETS allowance (EUA), phase 4 — clearing price of the primary auction on EEX, EUR per tonne CO₂', auction: last.name, contract: last.contract, zone: last.zone, file: used }; // only the latest clearing price is republished, with attribution — not the auction dataset
    },
  },
  rates: {
    source: 'Board of Governors of the Federal Reserve System (effective federal funds rate, FRED series DFF) and European Central Bank (deposit facility rate, ECB Data Portal series FM.D.U2.EUR.4F.KR.DFR.LEV)', home: 'https://data.ecb.europa.eu/data/datasets/FM/FM.D.U2.EUR.4F.KR.DFR.LEV',
    terms: 'Federal Reserve H.15 data are in the public domain; ECB statistics may be reused free of charge with attribution.',
    async run(get) {
      const out = {};
      try { const r = parseDateValueCsv(await get(`https://fred.stlouisfed.org/graph/fredgraph.csv?id=DFF&cosd=${iso(new Date(Date.now() - 30 * 86400e3))}`).then(okText)).pop(); if (r) Object.assign(out, { us_fed_funds_pct: r.v, us_date: r.t }); } catch (e) { out.us_error = e.message; }
      try { const t = (await get('https://data-api.ecb.europa.eu/service/data/FM/D.U2.EUR.4F.KR.DFR.LEV?lastNObservations=1&format=csvdata').then(okText)).trim().split(/\r?\n/), h = t[0].split(','), row = t[t.length - 1].split(','), v = Number(row[h.indexOf('OBS_VALUE')]), d = row[h.indexOf('TIME_PERIOD')]; if (Number.isFinite(v) && /^\d{4}-\d{2}-\d{2}$/.test(d)) Object.assign(out, { ecb_deposit_pct: v, ecb_date: d }); else throw new Error('unexpected reply'); } catch (e) { out.ecb_error = e.message; }
      if (out.us_fed_funds_pct == null && out.ecb_deposit_pct == null) throw new Error(out.us_error || out.ecb_error || 'no data');
      return out;
    },
  },
};

export async function buildSnapshot({ only = null, previous = null } = {}) {
  const feeds = {};
  for (const [id, f] of Object.entries(FEEDS)) { // in order: later feeds may use earlier ones (carbon uses fx)
    if (only && !only.includes(id)) { if (previous?.feeds?.[id]) feeds[id] = previous.feeds[id]; continue; }
    const urls = [], now = Date.now(), base = { source: f.source, home: f.home, terms: f.terms };
    try { const data = await f.run(makeGet(urls), feeds); feeds[id] = { ok: true, ts: now, fetched: new Date(now).toISOString(), ...base, url: urls[urls.length - 1] || f.home, urls, data }; }
    catch (e) { const old = previous?.feeds?.[id]; feeds[id] = { ok: false, error: String(e.message || e).slice(0, 200), tried: new Date(now).toISOString(), ts: old?.data ? old.ts : null, fetched: old?.data ? old.fetched : null, ...base, url: urls[urls.length - 1] || f.home, urls, data: old?.data ?? null }; }
  }
  return { schema: 1, generated: new Date().toISOString(), generator: 'tools/snapshot.mjs', interval_h: INTERVAL_H, feeds };
}

const isMain = process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url);
if (isMain) {
  const args = process.argv.slice(2), arg = (k) => (args.includes(k) ? args[args.indexOf(k) + 1] : null);
  const out = resolve(arg('--out') || join(root, 'data/snapshot.json')), only = arg('--only')?.split(',') || null;
  const when = (s) => Date.parse(s?.generated) || 0;
  let previous = null;
  try { previous = JSON.parse(await readFile(out, 'utf8')); } catch { /* first run */ }
  if (process.env.SNAPSHOT_PREV_URL) try { const r = await fetch(process.env.SNAPSHOT_PREV_URL, { signal: AbortSignal.timeout(20000), headers: { 'user-agent': UA } }); if (r.ok) { const p = await r.json(); if (p?.feeds && when(p) > when(previous)) previous = p; } } catch { /* not published yet */ }
  const snap = await buildSnapshot({ only, previous });
  await mkdir(dirname(out), { recursive: true });
  await writeFile(out, JSON.stringify(snap, null, 1) + '\n');
  for (const [id, f] of Object.entries(snap.feeds)) console.log(`${f.ok ? 'ok  ' : 'FAIL'} ${id.padEnd(13)} ${f.ok ? '' : f.error + (f.data ? ' — keeping the value fetched ' + f.fetched : ' — no earlier value')}`);
  console.log(`${out}: ${Object.values(snap.feeds).filter((f) => f.ok).length}/${Object.keys(snap.feeds).length} feeds, generated ${snap.generated}`);
}
