// Bundled world airport and runway database: `node tools/fetch-data.mjs [--from <dir with the CSV files>]`
// Also builds the bundled place-name database (`--places`, GeoNames) and the coarse terrain grid (`--terrain`, NOAA ETOPO 2022): see below.
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

// ---- places: GeoNames cities (CC BY 4.0) → js/data/places/cities.txt -----------------------------
// `node tools/fetch-data.mjs --places [--from <dir with cities15000.zip and admin1CodesASCII.txt>]`
// Every populated place with at least 15 000 inhabitants (and all capitals in the file), plus a second tier of places of 5 000 to
// 15 000 inhabitants under small/ (one file per first letter, read only when a search needs it), grouped by country and
// first-order division so the names of those are written once. Line: name|ascii name (if different)|lat|lon|elevation|population in thousands,
// with latitude and longitude in thousandths of a degree. Read by js/core/places.js; no network at run time.
async function buildPlaces() {
  const { unzip } = await import('./snapshot.mjs'), GEO = 'https://download.geonames.org/export/dump/', dir = join(root, 'js/data/places');
  const grab = async (name) => { if (from) return { buf: await readFile(join(from, name)), modified: null }; const r = await fetch(GEO + name, { signal: AbortSignal.timeout(180000) }); if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`); return { buf: Buffer.from(await r.arrayBuffer()), modified: r.headers.get('last-modified') }; };
  const z = await grab('cities15000.zip'), txt = unzip(z.buf).get('cities15000.txt').toString('utf8'), adm = (await grab('admin1CodesASCII.txt')).buf.toString('utf8');
  const admin = new Map(adm.split('\n').map((l) => l.split('\t')).filter((r) => r.length >= 2).map((r) => [r[0], r[1]]));
  const clean = (v) => String(v || '').replace(/[|#\n\r\t]/g, ' ').replace(/\s+/g, ' ').trim(), groups = new Map(), big = new Set(); let n = 0;
  // second tier, loaded on demand: places of 5 000 to 15 000 inhabitants (cities5000 minus the first tier), split by the first letter of the name
  const small = new Map(), txt5 = unzip((await grab('cities5000.zip')).buf).get('cities5000.txt').toString('utf8'); let nSmall = 0;
  const parse = (l) => { const f = l.split('\t'); if (f.length < 18) return null; const name = clean(f[1]), ascii = clean(f[2]), lat = Number(f[4]), lon = Number(f[5]), cc = f[8], a1 = admin.get(`${cc}.${f[10]}`) || '', pop = Number(f[14]) || 0, elev = f[15] !== '' && Number.isFinite(Number(f[15])) ? Number(f[15]) : Number(f[16]);
    if (!name || !/^[A-Z]{2}$/.test(cc) || !Number.isFinite(lat) || !Number.isFinite(lon) || Math.abs(lat) > 90 || Math.abs(lon) > 180) return null;
    return { id: f[0], ascii, cc, a1: clean(a1), pop, line: `${name}|${ascii && ascii !== name ? ascii : ''}|${Math.round(lat * 1000)}|${Math.round(lon * 1000)}|${Number.isFinite(elev) && elev > -9000 ? Math.round(elev) : ''}|${Math.round(pop / 1000)}` }; };
  for (const l of txt.split('\n')) { const f = l.split('\t'); if (f.length >= 18) big.add(f[0]); }
  for (const l of txt5.split('\n')) { const r = parse(l); if (!r || big.has(r.id)) continue; const b = /^[a-z]/.test(fold(r.ascii)) ? fold(r.ascii)[0] : '0'; if (!small.has(b)) small.set(b, new Map()); const g = small.get(b), key = `${r.cc}|${r.a1}`; if (!g.has(key)) g.set(key, []); g.get(key).push(r); nSmall++; }
  for (const l of txt.split('\n')) {
    const f = l.split('\t'); if (f.length < 18) continue;
    const name = clean(f[1]), ascii = clean(f[2]), lat = Number(f[4]), lon = Number(f[5]), cc = f[8], a1 = admin.get(`${cc}.${f[10]}`) || '', pop = Number(f[14]) || 0, elev = f[15] !== '' && Number.isFinite(Number(f[15])) ? Number(f[15]) : Number(f[16]);
    if (!name || !/^[A-Z]{2}$/.test(cc) || !Number.isFinite(lat) || !Number.isFinite(lon) || Math.abs(lat) > 90 || Math.abs(lon) > 180) continue;
    const key = `${cc}|${clean(a1)}`; if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push({ pop, line: `${name}|${ascii && ascii !== name ? ascii : ''}|${Math.round(lat * 1000)}|${Math.round(lon * 1000)}|${Number.isFinite(elev) && elev > -9000 ? Math.round(elev) : ''}|${Math.round(pop / 1000)}` }); n++;
  }
  if (n < 20000) throw new Error(`only ${n} places read: the GeoNames file layout may have changed`);
  const out = [...groups.keys()].sort().map((k) => `#${k}\n${groups.get(k).sort((a, b) => b.pop - a.pop).map((x) => x.line).join('\n')}`).join('\n') + '\n';
  await rm(dir, { recursive: true, force: true }); await mkdir(dir, { recursive: true });
  await writeFile(join(dir, 'cities.txt'), out);
  await mkdir(join(dir, 'small'), { recursive: true }); let smallBytes = 0; const buckets = {};
  for (const [b, g] of small) { const t = [...g.keys()].sort().map((k) => `#${k}\n${g.get(k).sort((x, y) => y.pop - x.pop).map((x) => x.line).join('\n')}`).join('\n') + '\n'; await writeFile(join(dir, 'small', `${b}.txt`), t); smallBytes += Buffer.byteLength(t); buckets[b] = [...g.values()].reduce((a, v) => a + v.length, 0); }
  if (nSmall < 15000) throw new Error(`only ${nSmall} smaller places read: the GeoNames cities5000 layout may have changed`);
  const date = z.modified ? new Date(z.modified).toISOString().slice(0, 10) : new Date().toISOString().slice(0, 10);
  await writeFile(join(dir, 'index.json'), JSON.stringify({ schema: 1, source: 'GeoNames', home: 'https://www.geonames.org/', files: [GEO + 'cities15000.zip', GEO + 'cities5000.zip', GEO + 'admin1CodesASCII.txt'], license: 'Creative Commons Attribution 4.0 (CC BY 4.0)', attribution: 'Place names: GeoNames (geonames.org), CC BY 4.0', registry: 'geonames', dataset_date: date, generated: new Date().toISOString(),
    selection: 'cities15000: populated places with more than 15 000 inhabitants, and capitals', changes: 'Reduced to name, ASCII name, country, first-order division, position (0.001°), elevation (GeoNames elevation, else its digital-elevation-model value) and population (thousands)', counts: { places: n, groups: groups.size, small: nSmall }, file: 'cities.txt', bytes: Buffer.byteLength(out),
    small: { selection: 'cities5000 without the places above: 5 000 to 15 000 inhabitants', dir: 'small/', by: 'first letter of the ASCII name (0 = other)', buckets, bytes: smallBytes } }, null, 1) + '\n');
  console.log(`GeoNames ${date}: ${n} places in ${groups.size} divisions, ${(Buffer.byteLength(out) / 1e6).toFixed(2)} MB; ${nSmall} smaller places in ${small.size} files, ${(smallBytes / 1e6).toFixed(2)} MB (loaded on demand)`);
}

// ---- terrain: NOAA ETOPO 2022 (public domain) → js/data/terrain/etopo-30m.bin --------------------
// `node tools/fetch-data.mjs --terrain`
// Reads the 60 arc-second surface-elevation grid through the NOAA NCEI THREDDS OPeNDAP service at every sixth
// point (0.1°), and stores for each 0.5° × 0.5° cell the mean of its 25 samples (land and ice surface; sea = 0),
// as int16 metres, north to south from 90°N, west to east from 180°W. Read by js/core/places.js.
async function buildTerrain() {
  const URL0 = 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/60s_surface_elev_netcdf/ETOPO_2022_v1_60s_N90W180_surface.nc', dir = join(root, 'js/data/terrain');
  const NJ = 360, NI = 720, K = 5, mean = new Float64Array(NJ * NI), top = new Float32Array(NJ * NI).fill(-1e9);
  for (let band = 0; band < 6; band++) { // six latitude bands of 30° keep each reply near 4 MB
    const r0 = band * 1800, q = `z[${r0 + 2}:6:${r0 + 1799}][2:6:21599]`, r = await fetch(`${URL0}.dods?${q}`, { signal: AbortSignal.timeout(600000) }); if (!r.ok) throw new Error(`ETOPO band ${band}: HTTP ${r.status}`);
    const buf = Buffer.from(await r.arrayBuffer()), at = buf.indexOf('\nData:\n'), p = at + 7, n = buf.readUInt32BE(p); if (at < 0 || n !== 300 * 3600) throw new Error(`ETOPO band ${band}: unexpected reply (${n} values)`);
    for (let j = 0; j < 300; j++) for (let i = 0; i < 3600; i++) { const v = Math.max(0, buf.readFloatBE(p + 8 + 4 * (j * 3600 + i))), row = NJ - 1 - Math.floor((band * 300 + j) / K), c = row * NI + Math.floor(i / K); mean[c] += v / (K * K); if (v > top[c]) top[c] = v; } // the file runs south to north
    process.stdout.write(` band ${band + 1}/6`);
  }
  console.log();
  const bin = Buffer.alloc(NJ * NI * 2); let mx = 0;
  for (let c = 0; c < NJ * NI; c++) { bin.writeInt16LE(Math.round(mean[c]), 2 * c); mx = Math.max(mx, top[c]); }
  if (!(mx > 7000 && mx < 9000)) throw new Error(`implausible highest elevation ${mx} m`);
  await mkdir(dir, { recursive: true });
  await writeFile(join(dir, 'etopo-30m.bin'), bin);
  await writeFile(join(dir, 'index.json'), JSON.stringify({ schema: 1, source: 'ETOPO 2022 Global Relief Model, 60 arc-second surface elevation', publisher: 'NOAA National Centers for Environmental Information', home: 'https://www.ncei.noaa.gov/products/etopo-global-relief-model', service: URL0, license: 'US Government work, public domain', attribution: 'Terrain: NOAA NCEI ETOPO 2022 (doi:10.25921/fd45-gt74)', registry: 'noaa-etopo', generated: new Date().toISOString(),
    method: 'Every sixth point of the 60 arc-second grid (0.1° spacing); per 0.5° cell the mean of the 25 samples; sea set to 0. Peaks inside a cell are higher than its mean. A coarse terrain model for route profiles and a last-resort site elevation — not for obstacle clearance.', grid: { step_deg: 0.5, nj: NJ, ni: NI, lat_north: 90, lon_west: -180 }, type: 'int16le', fields: ['mean_m'], file: 'etopo-30m.bin', bytes: bin.length, highest_m: Math.round(mx) }, null, 1) + '\n');
  console.log(`ETOPO 2022: ${NJ}×${NI} cells, ${(bin.length / 1e3).toFixed(0)} kB, highest sample ${Math.round(mx)} m`);
}

if (isMain && args.includes('--places')) await buildPlaces();
else if (isMain && args.includes('--terrain')) await buildTerrain();
else if (isMain) {
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
