// Carbon (emissions-allowance) prices from official publishers whose licence explicitly allows commercial reuse.
// Each market: descriptive fields + fetch(get) → { date, price, currency, instrument, url, lag_note? }, where `get(url, opt)` is a
// fetch-like function supplied by the caller. A fetch throws on any unexpected layout rather than return a doubtful number.
// Run directly (`node tools/carbon-markets.mjs`) to fetch every market and print the result.
import { pathToFileURL } from 'node:url';

const EU_ETS = ['AT', 'BE', 'BG', 'HR', 'CY', 'CZ', 'DK', 'EE', 'FI', 'FR', 'DE', 'GR', 'HU', 'IE', 'IT', 'LV', 'LT', 'LU', 'MT', 'NL', 'PL', 'PT', 'RO', 'SK', 'SI', 'ES', 'SE', 'IS', 'LI', 'NO'];
const MONTHS = ['january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december'];
const iso = (y, m, d) => new Date(Date.UTC(y, m, d)).toISOString().slice(0, 10); // m is 0-based; d = 0 gives the last day of the previous month
const today = () => new Date().toISOString().slice(0, 10);
const body = async (r, url) => { if (!r || !r.ok) throw new Error(`HTTP ${r?.status} for ${url}`); return r.text(); };
const cellText = (s) => s.replace(/<[^>]+>/g, ' ').replace(/&nbsp;|&#160;/g, ' ').replace(/&euro;|&#8364;/g, '€').replace(/&amp;/g, '&').replace(/\s+/g, ' ').trim();
/** Every "12 March 2027" in a text → ISO dates, in order of appearance. */
const longDates = (s) => [...s.matchAll(/(\d{1,2})\s+([A-Za-z]+)\s+(20\d\d)/g)].map((m) => { const k = MONTHS.indexOf(m[2].toLowerCase()); return k < 0 ? null : iso(Number(m[3]), k, Number(m[1])); }).filter(Boolean);
/** RFC 4180 text → rows of cells (quoted cells may hold commas and line breaks). */
function parseCsv(text) {
  const rows = []; let row = [], cell = '', quoted = false;
  const t = text.replace(/^﻿/, '');
  for (let i = 0; i < t.length; i++) {
    const c = t[i];
    if (quoted) { if (c !== '"') cell += c; else if (t[i + 1] === '"') { cell += '"'; i++; } else quoted = false; }
    else if (c === '"') quoted = true;
    else if (c === ',') { row.push(cell); cell = ''; }
    else if (c === '\n') { row.push(cell); rows.push(row); row = []; cell = ''; }
    else if (c !== '\r') cell += c;
  }
  if (cell !== '' || row.length) { row.push(cell); rows.push(row); }
  return rows;
}
const num = (s) => { const v = String(s ?? '').replace(/[\s  $]/g, '').replace(',', '.'); return /^\d+(\.\d+)?$/.test(v) ? Number(v) : NaN; };

// ---- EU ETS: price of CBAM certificates (European Commission) -----------------------------------------------------------------
const CBAM_URL = 'https://taxation-customs.ec.europa.eu/carbon-border-adjustment-mechanism/price-cbam-certificates_en';
/** Last day (Friday) of ISO week `w` of `y`, as yyyy-mm-dd. */
const isoWeekFriday = (y, w) => { const jan4 = new Date(Date.UTC(y, 0, 4)), mon = new Date(jan4.getTime() - ((jan4.getUTCDay() + 6) % 7) * 86400e3 + (w - 1) * 7 * 86400e3); return new Date(mon.getTime() + 4 * 86400e3).toISOString().slice(0, 10); };
/** Every date written in a cell, as yyyy-mm-dd, in order: "5 October 2026", "05/10/2026", "5.10.2026", "2026-10-05". */
function anyDates(text) {
  const found = [];
  for (const m of String(text).matchAll(/(\d{4})-(\d{2})-(\d{2})|(\d{1,2})[./](\d{1,2})[./](20\d\d)/g)) { const [y, mo, d] = m[1] ? [m[1], m[2], m[3]] : [m[6], m[5], m[4]]; if (mo >= 1 && mo <= 12 && d >= 1 && d <= 31) found.push({ at: m.index, d: iso(Number(y), Number(mo) - 1, Number(d)) }); }
  const t = String(text); let from = 0; for (const d of longDates(t)) { const day = String(Number(d.slice(8))), at = t.indexOf(day, from); found.push({ at: at < 0 ? from : at, d }); from = at < 0 ? from : at + 1; }
  return found.sort((a, b) => a.at - b.at).map((x) => x.d);
}
/**
 * The Commission's CBAM certificate price table(s) → [{ period, quarterly, date, published, price }], newest first.
 * Reads the quarterly layout of 2026 (period "Q3 2026") and tolerates the weekly layout announced for 2027, whatever
 * form the week takes: ISO week ("Week 2 2027", "W02 2027", "2027-W02"), a date range in words or figures, or a
 * single date. The price column is the one whose heading says price (in € / EUR); cells may carry a € sign. A row
 * whose price is still empty is skipped; a price that cannot be read is an error, never a guess.
 */
export function parseCbamPrices(html) {
  const out = [];
  for (const t of html.matchAll(/<table\b[\s\S]*?<\/table>/g)) {
    const rows = [...t[0].matchAll(/<tr\b[^>]*>([\s\S]*?)<\/tr>/g)].map((r) => [...r[1].matchAll(/<(t[hd])\b[^>]*>([\s\S]*?)<\/t[hd]>/g)].map((m) => ({ th: m[1] === 'th', text: cellText(m[2]) })));
    const hi = rows.findIndex((r) => r.length >= 2 && r.some((c) => /price/i.test(c.text)) && !r.some((c) => /^[€\s]*\d{1,3}([.,]\d{1,2})?[€\s]*$/.test(c.text))); if (hi < 0) continue;
    const head = rows[hi].map((c) => c.text), euro = (h) => /price/i.test(h) && /(€|EUR|euro)/i.test(h);
    const iPrice = head.some(euro) ? head.findIndex(euro) : head.findIndex((h) => /price/i.test(h) && !/publication|date/i.test(h)), iPub = head.findIndex((h) => /publi/i.test(h));
    if (iPrice < 0) continue;
    const iPeriod = head.findIndex((_, i) => i !== iPrice && i !== iPub); if (iPeriod < 0) continue;
    for (const r of rows.slice(hi + 1)) {
      const c = r.map((x) => x.text); if (c.length !== head.length) continue;
      const raw = c[iPrice].replace(/€|EUR|euros?/gi, '').replace(/[\s  ]/g, ''); if (raw === '' || /^(-|–|—|n\/?a|tbc|tbd)$/i.test(raw)) continue; // a period whose price is not out yet
      if (!/^\d{1,3}([.,]\d{1,4})?$/.test(raw)) throw new Error(`CBAM price table: unreadable price "${c[iPrice]}"`);
      const price = Number(raw.replace(',', '.')), label = c[iPeriod], published = iPub >= 0 ? anyDates(c[iPub])[0] || null : null;
      const q = /^Q([1-4])[\s/-]*(20\d\d)$|^(20\d\d)[\s/-]*Q([1-4])$/i.exec(label), wk = /\b(?:week|wk|w)\s*\.?\s*(\d{1,2})\D{1,6}(20\d\d)|\b(20\d\d)\s*[-/ ]?\s*(?:week|wk|w)\s*\.?\s*(\d{1,2})\b/i.exec(label), inLabel = anyDates(label);
      const date = q ? iso(Number(q[2] || q[3]), Number(q[1] || q[4]) * 3, 0) : wk ? isoWeekFriday(Number(wk[2] || wk[3]), Number(wk[1] || wk[4])) : inLabel.length ? inLabel.sort()[inLabel.length - 1] : published; // quarter → its last day; week → its Friday or last date given
      if (!date || (wk && !(Number(wk[1] || wk[4]) >= 1 && Number(wk[1] || wk[4]) <= 53))) throw new Error(`CBAM price table: cannot date the period "${label}"`);
      out.push({ period: label, quarterly: !!q, date, published, price });
    }
  }
  return out.sort((a, b) => (a.date < b.date ? 1 : a.date > b.date ? -1 : 0));
}

// ---- WCI (California–Québec): joint auction results, Données Québec -----------------------------------------------------------
const QC_DATASET = 'resultats-et-statistiques-des-ventes-aux-encheres-d-unites-d-emission';
const QC_API = `https://www.donneesquebec.ca/recherche/api/3/action/package_show?id=${QC_DATASET}`;
const QC_CSV = 'https://www.donneesquebec.ca/recherche/dataset/4719fd2b-458b-4058-ac3a-12bed98f50b6/resource/10273e0d-0aae-42fd-9366-5fd2181e0a5a/download/donnees-ouvertes-encheres.csv';
/** MELCCFP auction file → [{ date, no, usd, cad, floor_usd, fx }], newest first (current-vintage settlement price of each auction). */
export function parseQuebecAuctions(csv) {
  const rows = parseCsv(csv); if (rows.length < 2) throw new Error('Québec auction file is empty');
  const H = rows[0].map((h) => h.trim()), col = (name) => { const i = H.indexOf(name); if (i < 0) throw new Error(`Québec auction file: column ${name} is missing`); return i; };
  const C = { date: col('Date_vente'), no: col('No_vente'), usd: col('Prix_final_USD_present'), cad: col('Prix_final_CAD_present'), floor: col('Prix_min_USD_present'), fx: col('Taux_change') };
  return rows.slice(1).filter((r) => /^\d{4}-\d{2}-\d{2}$/.test((r[C.date] || '').trim()))
    .map((r) => ({ date: r[C.date].trim(), no: (r[C.no] || '').trim(), usd: num(r[C.usd]), cad: num(r[C.cad]), floor_usd: num(r[C.floor]), fx: num(r[C.fx]) }))
    .sort((a, b) => (a.date < b.date ? 1 : a.date > b.date ? -1 : 0));
}

export const CARBON_MARKETS = {
  'eu-ets-cbam': {
    id: 'eu-ets-cbam',
    market: 'EU ETS (EUA) — average auction clearing price (price of CBAM certificates)',
    scheme: 'EU ETS', countries: EU_ETS, currency: 'EUR', unit: 't CO2e',
    publisher: 'European Commission, Directorate-General for Taxation and Customs Union',
    home: CBAM_URL,
    licence: 'Creative Commons Attribution 4.0 International (CC BY 4.0) — European Commission reuse policy (Decision 2011/833/EU)',
    licenceUrl: 'https://commission.europa.eu/legal-notice_en',
    quote: 'This means that reuse is allowed, provided appropriate credit is given and changes are indicated.',
    attribution: 'Source: European Commission, "Price of CBAM certificates" (© European Union, CC BY 4.0).',
    async fetch(get) {
      const all = parseCbamPrices(await body(await get(CBAM_URL), CBAM_URL));
      if (!all.length) throw new Error('CBAM price table not found or without a published price (page layout changed?)');
      const p = all[0], now = today();
      if (!(p.price >= 5 && p.price <= 500)) throw new Error(`implausible CBAM certificate price ${p.price}`);
      if (p.published && (p.published > now || p.published < p.date.slice(0, 8) + '01')) throw new Error(`CBAM price table: inconsistent dates (${p.period}, published ${p.published})`);
      if (p.date > iso(new Date().getUTCFullYear(), new Date().getUTCMonth() + 1, 0) || p.date < iso(new Date().getUTCFullYear() - 1, new Date().getUTCMonth(), 1)) throw new Error(`CBAM price table: newest period ${p.period} (${p.date}) is not recent`);
      return {
        date: p.date, price: p.price, currency: 'EUR', period: p.period, quarterly: p.quarterly, published: p.published, url: CBAM_URL,
        instrument: `Price of CBAM certificates for ${p.period}: volume-weighted average of the clearing prices of EU ETS allowance (EUA) auctions on the common auction platform over that period, EUR per tonne CO2e (Regulation (EU) 2023/956, Implementing Regulation (EU) 2025/2548)`,
        lag_note: p.quarterly ? 'A quarterly average, published in the first week after the quarter ends; not a daily or spot price. The date is the last day of the quarter. Weekly averages are announced from 2027.' : 'An average over the stated period, published after it ends; not a daily or spot price.',
      };
    },
  },
  'wci-joint-auction': {
    id: 'wci-joint-auction',
    market: 'WCI linked market (California Cap-and-Trade / Québec SPEDE) — joint auction settlement price',
    scheme: 'WCI (California–Québec)', countries: ['US', 'CA'], subdivisions: ['US-CA', 'CA-QC'], currency: 'USD', unit: 't CO2e',
    publisher: 'Gouvernement du Québec, Ministère de l’Environnement, de la Lutte contre les changements climatiques, de la Faune et des Parcs (MELCCFP)',
    home: `https://www.donneesquebec.ca/recherche/dataset/${QC_DATASET}`,
    licence: 'Creative Commons Attribution 4.0 International (CC BY 4.0)',
    licenceUrl: 'https://www.donneesquebec.ca/licence/',
    quote: 'Cette licence permet à d’autres personnes de distribuer, remixer, arranger et adapter votre œuvre, même à des fins commerciales, tant qu’on vous attribue le crédit de la création originale en citant votre nom.',
    attribution: 'Source: Gouvernement du Québec (MELCCFP), « Résultats et statistiques des ventes aux enchères d’unités d’émission de GES », Données Québec, CC BY 4.0.',
    async fetch(get) {
      let url = QC_CSV;
      try { // the portal's catalogue names the current CSV file; the address above is the fall-back
        const res = JSON.parse(await body(await get(QC_API), QC_API)).result.resources.filter((r) => /^csv$/i.test(r.format || '') && /^https:\/\/www\.donneesquebec\.ca\/.+\.csv$/i.test(r.url || ''));
        if (res.length === 1) url = res[0].url;
      } catch { /* keep the fall-back */ }
      const all = parseQuebecAuctions(await body(await get(url), url)), a = all[0], now = today();
      if (!a || all.length < 20) throw new Error('Québec auction file: too few auctions (layout changed?)');
      if (!Number.isFinite(a.usd) || !Number.isFinite(a.floor_usd)) throw new Error(`Québec auction file: no settlement price for the ${a.date} auction`);
      if (!(a.usd >= 5 && a.usd <= 300) || a.usd < a.floor_usd - 0.005) throw new Error(`implausible settlement price ${a.usd} USD (floor ${a.floor_usd})`);
      if (Number.isFinite(a.cad) && Number.isFinite(a.fx) && Math.abs(a.usd * a.fx - a.cad) > 0.05 * a.cad) throw new Error('Québec auction file: USD and CAD settlement prices disagree');
      if (a.date > now || a.date < iso(new Date().getUTCFullYear() - 1, new Date().getUTCMonth(), 1)) throw new Error(`Québec auction file: newest auction ${a.date} is not recent`);
      return {
        date: a.date, price: a.usd, currency: 'USD', price_cad: Number.isFinite(a.cad) ? a.cad : null, auction: a.no, url,
        instrument: `Settlement price of current-vintage allowances at California–Québec joint auction no. ${a.no}, USD per allowance (1 t CO2e)`,
        lag_note: 'Joint auctions are held four times a year (February, May, August, November); the open-data file is updated a few weeks after each auction. Not a secondary-market price.',
      };
    },
  },
};

// ---- run directly: fetch every market and print it ----------------------------------------------------------------------------
if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  const get = (url, opt = {}) => fetch(url, { ...opt, headers: { 'user-agent': 'AeroSuite26-snapshot/1.0', ...(opt.headers || {}) }, signal: AbortSignal.timeout(45000) });
  let failed = 0;
  for (const m of Object.values(CARBON_MARKETS)) {
    try { const d = await m.fetch(get); console.log(`${m.id}\t${d.date}\t${d.price}\t${d.currency}\t${d.url}`); }
    catch (e) { failed++; console.log(`${m.id}\tERROR\t${e?.cause?.code || e?.message || e}`); }
  }
  process.exitCode = failed ? 1 : 0;
}
