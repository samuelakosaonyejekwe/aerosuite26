// Bundled world airport and runway database: `node tools/fetch-data.mjs [--from <dir with the CSV files>]`
// Downloads the public-domain OurAirports datasets, keeps open aerodromes and heliports that carry usable
// data, and writes compact, spatially indexed text tiles (10° × 10°) plus two small look-up files to
// js/data/airports/. The app reads them through js/core/airports.js and needs no network for them.
//
// Source: OurAirports (https://ourairports.com/data/), released to the public domain by its contributors.

import { readFile, writeFile, mkdir, readdir, rm } from 'node:fs/promises';
import { join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { TILE_DEG, TYPES, SURFACES, SUFFIXES, tileOf, fold, reciprocal } from '../js/core/airports.js';

const root = resolve(fileURLToPath(new URL('..', import.meta.url))), outDir = join(root, 'js/data/airports');
const BASE = 'https://davidmegginson.github.io/ourairports-data/';
const args = process.argv.slice(2), from = args.includes('--from') ? args[args.indexOf('--from') + 1] : null;
const FT = 0.3048, MIN_STRIP_M = 400; // uncoded small strips shorter than this are left out to keep the bundle small

/** RFC 4180 CSV → array of objects (quoted fields may hold commas, quotes and line breaks). */
export function parseCsv(text) {
  const rows = []; let row = [], cell = '', q = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (q) { if (ch === '"') { if (text[i + 1] === '"') { cell += '"'; i++; } else q = false; } else cell += ch; }
    else if (ch === '"') q = true;
    else if (ch === ',') { row.push(cell); cell = ''; }
    else if (ch === '\n') { row.push(cell); rows.push(row); row = []; cell = ''; }
    else if (ch !== '\r') cell += ch;
  }
  if (cell !== '' || row.length) { row.push(cell); rows.push(row); }
  const head = rows.shift();
  return rows.filter((r) => r.length === head.length).map((r) => Object.fromEntries(head.map((h, i) => [h, r[i]])));
}

/** OurAirports free-text surface → one of the plain words in SURFACES (never guessed: unclear codes stay "unknown"). */
export function surfaceWord(raw) {
  const s = String(raw || '').toUpperCase().normalize('NFD').replace(/[̀-ͯ]/g, '').trim();
  if (!s) return 'unknown';
  if (/WATER|^WAT\b/.test(s)) return 'water';
  if (/\bICE\b/.test(s)) return 'ice';
  if (/SNOW/.test(s)) return 'snow';
  const asp = /ASP|BIT|TARMAC|MACADAM|\bPEM\b/.test(s), con = /CON|\bPCC\b|CEMENT/.test(s);
  if (asp && con) return 'paved';
  if (asp) return 'asphalt';
  if (con) return 'concrete';
  if (/UNPAVED|UNSEALED/.test(s)) return 'unpaved';
  if (/PAVED|SEALED|\bPER\b|\bPAD\b|\bCOP\b|COMPOSITE|BRICK|ROOF/.test(s)) return 'paved';
  if (/TURF|GRASS|\bGRS\b|\bSOD\b/.test(s)) return 'grass';
  if (/GRAV|GRVL|\bGVL\b|\bGRV\b|PICARRA|STONE|ROCK/.test(s)) return 'gravel';
  if (/LATERITE|\bLAT\b|CLAY|CORAL|CALICHE|TREATED|COMPACT|OILED|MURRAM|HARD/.test(s)) return 'compacted';
  if (/SAND|\bSAN\b/.test(s)) return 'sand';
  if (/DIRT|EARTH|\bGRE\b|SOIL|\bTER\b|GROUND|MUD|SILT|NATURAL/.test(s)) return 'dirt';
  if (/MATS?\b|METAL|\bMET\b|\bPSP\b|STEEL|ALUM/.test(s)) return 'metal';
  if (/WOOD|DECK/.test(s)) return 'wood';
  return 'unknown';
}

async function fetchCsv(name) {
  if (from) return { text: await readFile(join(from, name), 'utf8'), modified: null };
  const r = await fetch(BASE + name, { signal: AbortSignal.timeout(120000) }); if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
  return { text: await r.text(), modified: r.headers.get('last-modified') };
}

const isMain = process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url);
if (isMain) {
  const [a, r] = await Promise.all([fetchCsv('airports.csv'), fetchCsv('runways.csv')]);
  const datasetDate = (args.includes('--date') ? args[args.indexOf('--date') + 1] : a.modified ? new Date(a.modified).toISOString().slice(0, 10) : new Date().toISOString().slice(0, 10));
  const airports = parseCsv(a.text), runways = parseCsv(r.text);
  const num = (v) => (v === '' || v == null || !Number.isFinite(Number(v)) ? null : Number(v));
  const byAirport = new Map(); for (const x of runways) { if (!byAirport.has(x.airport_ident)) byAirport.set(x.airport_ident, []); byAirport.get(x.airport_ident).push(x); }
  const clean = (s) => String(s || '').replace(/[|;~\n\r\t]/g, ' ').replace(/\s+/g, ' ').trim();
  const code = (s) => clean(s).toUpperCase().replace(/[^A-Z0-9-]/g, '');

  const tiles = new Map(), codeLines = new Map(), nameLines = [], count = { airports: 0, runways: 0, byType: {}, bySurface: {} };
  let seenIn = 0;
  for (const ap of airports) {
    const type = { large_airport: 'L', medium_airport: 'M', small_airport: 'S', heliport: 'H', seaplane_base: 'W', balloonport: 'B' }[ap.type]; if (!type) continue; // "closed" and unknown types are dropped
    const lat = num(ap.latitude_deg), lon = num(ap.longitude_deg); if (lat == null || lon == null || Math.abs(lat) > 90 || Math.abs(lon) > 180 || (lat === 0 && lon === 0)) continue;
    seenIn++;
    const ident = code(ap.ident); if (!ident) continue;
    const icao = code(ap.icao_code), iata = code(ap.iata_code).slice(0, 3);
    const rws = (byAirport.get(ap.ident) || []).filter((x) => num(x.length_ft) > 0);
    const open = rws.filter((x) => x.closed !== '1');
    // "usable data": every large and medium airport; anything with an ICAO/IATA code or scheduled service; other
    // aerodromes when they have an open runway of known length (small strips: at least MIN_STRIP_M).
    const coded = !!(icao || iata), notable = type === 'L' || type === 'M' || !!iata || ap.scheduled_service === 'yes';
    const longest = open.reduce((m, x) => Math.max(m, num(x.length_ft) * FT), 0);
    const keep = type === 'L' || type === 'M' || coded || ap.scheduled_service === 'yes' || (type === 'H' ? false : type === 'S' ? longest >= MIN_STRIP_M : open.length > 0);
    if (!keep) continue;

    const t = tileOf(lat, lon), lat0 = t.lat0, lon0 = t.lon0;
    const latI = Math.round((lat - lat0) * 1e4), lonI = Math.round((lon - lon0) * 1e4), latR = lat0 + latI / 1e4, lonR = lon0 + lonI / 1e4;
    const elev = num(ap.elevation_ft) == null ? '' : Math.round(num(ap.elevation_ft) * FT);
    let name = clean(ap.name) || ident; const sfx = SUFFIXES.reduce((best, s, i) => (i > 0 && name.endsWith(s) && name.length > s.length && s.length > SUFFIXES[best].length ? i : best), 0); if (sfx > 0) name = name.slice(0, -SUFFIXES[sfx].length) + '~' + sfx.toString(36);
    const end = (x, p) => { const la = num(x[`${p}_latitude_deg`]), lo = num(x[`${p}_longitude_deg`]), el = num(x[`${p}_elevation_ft`]), dt = num(x[`${p}_displaced_threshold_ft`]); const okPos = la != null && lo != null && Math.abs(la - latR) < 1 && Math.abs(lo - lonR) < 1 && !(la === 0 && lo === 0); return [okPos ? Math.round((la - latR) * 1e5) : '', okPos ? Math.round((lo - lonR) * 1e5) : '', el == null ? '' : Math.round(el * FT), dt ? Math.round(dt * FT) : '']; };
    const hdg = (v) => { const h = num(v); return h == null || h < 0 || h > 360 ? null : Math.round(h * 10) / 10; };
    const rwTxt = rws.sort((x, y) => num(y.length_ft) - num(x.length_ft)).slice(0, 12).map((x) => {
      const sw = surfaceWord(x.surface), le = hdg(x.le_heading_degT), he = hdg(x.he_heading_degT), wid = num(x.width_ft) > 0 ? Math.round(num(x.width_ft) * FT) : '';
      count.bySurface[sw] = (count.bySurface[sw] || 0) + 1; count.runways++;
      const heOut = he == null || (le != null && Math.abs(((he - le + 360) % 360) - 180) < 0.6) ? '' : he; // omitted when it is simply the reciprocal
      const leId = code(x.le_ident), heId = code(x.he_ident), f = [heId === reciprocal(leId) ? leId : `${leId}/${heId}`, Math.round(num(x.length_ft) * FT), wid, SURFACES.indexOf(sw), (x.lighted === '1' ? 1 : 0) + (x.closed === '1' ? 2 : 0), le ?? '', heOut, ...end(x, 'le'), ...end(x, 'he')];
      while (f.length && f[f.length - 1] === '') f.pop();
      return f.join(',');
    }).join(';');
    const cc = clean(ap.iso_country).toUpperCase().slice(0, 2), line = [ident, icao === ident ? '1' : icao, iata, name, latI, lonI, elev, rwTxt].join('|').replace(/\|+$/, '');
    if (!tiles.has(t.id)) tiles.set(t.id, []); tiles.get(t.id).push({ line, group: `${cc}|${type}`, ident });
    if (!codeLines.has(t.id)) codeLines.set(t.id, new Set()); for (const c of [/^[A-Z]{2}-\d+$/.test(ident) ? '' : ident, icao, iata]) if (c) codeLines.get(t.id).add(c); // OurAirports' own serial identifiers are not indexed
    if (notable) nameLines.push({ rank: 'LMSWHB'.indexOf(type), type, line: [ident, t.id, name, fold(clean(ap.municipality)).split(' ').every((w) => ` ${fold(clean(ap.name))} `.includes(` ${w} `)) ? '' : clean(ap.municipality)].join('|').replace(/\|+$/, '') });
    count.airports++; count.byType[TYPES[type]] = (count.byType[TYPES[type]] || 0) + 1;
  }

  await mkdir(outDir, { recursive: true });
  for (const f of await readdir(outDir)) await rm(join(outDir, f));
  let bytes = 0, maxTile = 0; const tileIndex = {};
  const put = async (name, text) => { const b = Buffer.byteLength(text); bytes += b; await writeFile(join(outDir, name), text); return b; };
  for (const [id, list] of [...tiles].sort()) { list.sort((x, y) => (x.group < y.group ? -1 : x.group > y.group ? 1 : x.ident < y.ident ? -1 : 1)); let g = null; const b = await put(`t${id}.txt`, list.map((x) => (x.group !== g ? `#${(g = x.group)}\n` : '') + x.line).join('\n') + '\n'); maxTile = Math.max(maxTile, b); tileIndex[id] = list.length; }
  const codesB = await put('codes.txt', [...codeLines].sort().map(([id, set]) => `${id}:${[...set].sort().join(' ')}`).join('\n') + '\n');
  const namesB = await put('names.txt', (() => { let g = null; return nameLines.sort((x, y) => x.rank - y.rank || (x.line < y.line ? -1 : 1)).map((x) => (x.type !== g ? `#${(g = x.type)}\n` : '') + x.line).join('\n') + '\n'; })());
  const index = {
    schema: 1, source: 'OurAirports', home: 'https://ourairports.com/data/', license: 'Public domain (dedicated by the OurAirports contributors)',
    files: [BASE + 'airports.csv', BASE + 'runways.csv'], dataset_date: datasetDate, generated: new Date().toISOString(),
    selection: `Open large and medium airports; every other open aerodrome or heliport with an ICAO/IATA code or scheduled service; uncoded seaplane bases and balloonports with a runway of known length; uncoded small aerodromes with an open runway of at least ${MIN_STRIP_M} m. Uncoded heliports are omitted. Closed aerodromes are omitted; closed runways are kept and flagged.`,
    tile_deg: TILE_DEG, counts: { input_open: seenIn, airports: count.airports, runways: count.runways, named_in_search: nameLines.length, by_type: count.byType, by_surface: count.bySurface },
    tiles: tileIndex,
  };
  await put('index.json', JSON.stringify(index));
  console.log(`OurAirports ${datasetDate}: ${count.airports} aerodromes, ${count.runways} runways in ${tiles.size} tiles`);
  console.log(`tiles+index ${(bytes / 1e6).toFixed(2)} MB total  ·  largest tile ${(maxTile / 1e3).toFixed(0)} kB  ·  codes ${(codesB / 1e3).toFixed(0)} kB  ·  names ${(namesB / 1e3).toFixed(0)} kB`);
}
