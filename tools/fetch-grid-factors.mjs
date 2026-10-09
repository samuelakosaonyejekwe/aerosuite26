#!/usr/bin/env node
// Rebuilds js/data/grid-factors.json: national electricity-grid emission factors
// (g CO2e per kWh of electricity generated, latest available year per country).
//
//   node tools/fetch-grid-factors.mjs
//
// Source: Ember, "Yearly Electricity Data" (long-format CSV, no key required).
// Licence: CC BY 4.0 — https://ember-energy.org/creative-commons/
// Dependency-free (Node 20+). The existing file is only replaced when the download
// parses and passes the plausibility checks below.

import { writeFile, mkdir } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const DATA_URL = 'https://files.ember-energy.org/public-downloads/generation/outputs/release_generation_yearly_global.csv';
const OUT = resolve(dirname(fileURLToPath(import.meta.url)), '..', 'js', 'data', 'grid-factors.json');

const G_MIN = 0, G_MAX = 1300;   // plausible range, g/kWh
const MIN_COUNTRIES = 150;
const MIN_YEAR = 2015;           // ignore series that stopped long ago
const PER_LINE = 10;             // countries per output line

// ISO 3166-1 alpha-3 -> alpha-2 (all officially assigned codes, plus XKX = Kosovo, user-assigned).
const ISO3_TO_2 = Object.fromEntries((
  'ABW:AW AFG:AF AGO:AO AIA:AI ALA:AX ALB:AL AND:AD ARE:AE ARG:AR ARM:AM ASM:AS ATA:AQ ATF:TF ATG:AG AUS:AU AUT:AT AZE:AZ ' +
  'BDI:BI BEL:BE BEN:BJ BES:BQ BFA:BF BGD:BD BGR:BG BHR:BH BHS:BS BIH:BA BLM:BL BLR:BY BLZ:BZ BMU:BM BOL:BO BRA:BR BRB:BB BRN:BN BTN:BT BVT:BV BWA:BW ' +
  'CAF:CF CAN:CA CCK:CC CHE:CH CHL:CL CHN:CN CIV:CI CMR:CM COD:CD COG:CG COK:CK COL:CO COM:KM CPV:CV CRI:CR CUB:CU CUW:CW CXR:CX CYM:KY CYP:CY CZE:CZ ' +
  'DEU:DE DJI:DJ DMA:DM DNK:DK DOM:DO DZA:DZ ECU:EC EGY:EG ERI:ER ESH:EH ESP:ES EST:EE ETH:ET FIN:FI FJI:FJ FLK:FK FRA:FR FRO:FO FSM:FM ' +
  'GAB:GA GBR:GB GEO:GE GGY:GG GHA:GH GIB:GI GIN:GN GLP:GP GMB:GM GNB:GW GNQ:GQ GRC:GR GRD:GD GRL:GL GTM:GT GUF:GF GUM:GU GUY:GY ' +
  'HKG:HK HMD:HM HND:HN HRV:HR HTI:HT HUN:HU IDN:ID IMN:IM IND:IN IOT:IO IRL:IE IRN:IR IRQ:IQ ISL:IS ISR:IL ITA:IT JAM:JM JEY:JE JOR:JO JPN:JP ' +
  'KAZ:KZ KEN:KE KGZ:KG KHM:KH KIR:KI KNA:KN KOR:KR KWT:KW LAO:LA LBN:LB LBR:LR LBY:LY LCA:LC LIE:LI LKA:LK LSO:LS LTU:LT LUX:LU LVA:LV ' +
  'MAC:MO MAF:MF MAR:MA MCO:MC MDA:MD MDG:MG MDV:MV MEX:MX MHL:MH MKD:MK MLI:ML MLT:MT MMR:MM MNE:ME MNG:MN MNP:MP MOZ:MZ MRT:MR MSR:MS MTQ:MQ MUS:MU MWI:MW MYS:MY MYT:YT ' +
  'NAM:NA NCL:NC NER:NE NFK:NF NGA:NG NIC:NI NIU:NU NLD:NL NOR:NO NPL:NP NRU:NR NZL:NZ OMN:OM PAK:PK PAN:PA PCN:PN PER:PE PHL:PH PLW:PW PNG:PG POL:PL PRI:PR PRK:KP PRT:PT PRY:PY PSE:PS PYF:PF ' +
  'QAT:QA REU:RE ROU:RO RUS:RU RWA:RW SAU:SA SDN:SD SEN:SN SGP:SG SGS:GS SHN:SH SJM:SJ SLB:SB SLE:SL SLV:SV SMR:SM SOM:SO SPM:PM SRB:RS SSD:SS STP:ST SUR:SR SVK:SK SVN:SI SWE:SE SWZ:SZ SXM:SX SYC:SC SYR:SY ' +
  'TCA:TC TCD:TD TGO:TG THA:TH TJK:TJ TKL:TK TKM:TM TLS:TL TON:TO TTO:TT TUN:TN TUR:TR TUV:TV TWN:TW TZA:TZ UGA:UG UKR:UA UMI:UM URY:UY USA:US UZB:UZ ' +
  'VAT:VA VCT:VC VEN:VE VGB:VG VIR:VI VNM:VN VUT:VU WLF:WF WSM:WS XKX:XK YEM:YE ZAF:ZA ZMB:ZM ZWE:ZW'
).split(' ').map((p) => p.split(':')));

// Minimal RFC 4180 line splitter (quoted fields, doubled quotes). Records never span lines in this file.
function splitCsvLine(line) {
  const out = [];
  let cur = '', q = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (q) {
      if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; }
      else cur += c;
    } else if (c === '"') q = true;
    else if (c === ',') { out.push(cur); cur = ''; }
    else cur += c;
  }
  out.push(cur);
  return out;
}

function fail(msg) {
  console.error(`fetch-grid-factors: ${msg}`);
  console.error('js/data/grid-factors.json was NOT modified.');
  process.exit(1);
}

async function main() {
  console.log(`Downloading ${DATA_URL}`);
  let text;
  try {
    const res = await fetch(DATA_URL, { headers: { 'User-Agent': 'Mozilla/5.0 (compatible; grid-factors-fetch)' }, redirect: 'follow' });
    if (!res.ok) fail(`HTTP ${res.status} ${res.statusText}`);
    text = await res.text();
  } catch (e) {
    fail(`download failed: ${e.message}`);
  }

  const lines = text.split(/\r?\n/);
  const header = splitCsvLine(lines[0].replace(/^﻿/, ''));
  const col = (name) => { const i = header.indexOf(name); if (i < 0) fail(`column "${name}" missing — the source format changed`); return i; };
  const cArea = col('Area'), cIso = col('ISO 3 code'), cYear = col('Year'), cType = col('Area type');
  const cSrc = col('Electricity source'), cG = col('Emissions intensity (gCO2e/kWh)');

  const best = new Map();      // key -> { year, g } : latest year with an in-range value
  const rejected = new Map();  // key -> out-of-range values that were ignored
  let world = null;
  const unmapped = new Set();

  for (let n = 1; n < lines.length; n++) {
    const line = lines[n];
    if (!line || !line.includes('Total generation')) continue;
    const f = splitCsvLine(line);
    if (f[cSrc] !== 'Total generation' || f[cG] === '' || f[cG] == null) continue;
    const year = Number(f[cYear]), g = Number(f[cG]);
    if (!Number.isInteger(year) || !Number.isFinite(g)) continue;

    let key;
    if (f[cType] === 'Country or economy') {
      key = ISO3_TO_2[f[cIso]];
      if (!key) { unmapped.add(`${f[cArea]} (${f[cIso] || 'no code'})`); continue; }
    } else if (f[cArea] === 'World') key = 'WORLD';
    else continue;             // other aggregates (EU, G20, OECD, continents …) are skipped

    if (g < G_MIN || g > G_MAX) {
      rejected.set(key, [...(rejected.get(key) || []), `${year}: ${g}`]);
      continue;
    }
    const prev = best.get(key);
    if (!prev || year > prev.year) best.set(key, { year, g: Math.round(g) });
  }

  world = best.get('WORLD') || null;
  best.delete('WORLD');

  const countries = {};
  const stale = [];
  for (const k of [...best.keys()].sort()) {
    const v = best.get(k);
    if (v.year < MIN_YEAR) { stale.push(`${k} (${v.year})`); continue; }
    countries[k] = v;
  }

  const n = Object.keys(countries).length;
  if (n < MIN_COUNTRIES) fail(`only ${n} countries parsed (expected at least ${MIN_COUNTRIES})`);
  if (!world) fail('no plausible World value found');
  for (const [k, v] of Object.entries({ ...countries, WORLD: world })) {
    if (!(v.g >= G_MIN && v.g <= G_MAX) || !(v.year >= MIN_YEAR && v.year <= new Date().getUTCFullYear())) fail(`implausible entry ${k}: ${JSON.stringify(v)}`);
  }

  const meta = {
    schema: 1,
    source: 'Ember — Yearly Electricity Data (emissions intensity of total electricity generation)',
    publisher: 'Ember (Ember Energy Research CIC)',
    home: 'https://ember-energy.org/data/yearly-electricity-data/',
    dataUrl: DATA_URL,
    licence: 'CC BY 4.0',
    licenceUrl: 'https://ember-energy.org/creative-commons/',
    quote: 'Ember content is released under a Creative Commons Attribution Licence (CC-BY-4.0)',
    attribution: 'Grid emission factors: Ember, Yearly Electricity Data (ember-energy.org), CC BY 4.0. Values rounded to whole g/kWh; latest year per country.',
    retrieved: new Date().toISOString().slice(0, 10),
    unit: 'gCO2e/kWh',
    method: 'Generation-based, lifecycle: Ember\'s "Emissions intensity (gCO2e/kWh)" of total generation, i.e. CO2-equivalent emissions including upstream and supply-chain stages (so hydro, solar, wind and nuclear are non-zero) divided by electricity generated in the country; no import/export or grid-loss adjustment.',
    world,
  };

  // Compact but diff-friendly: metadata one key per line, ~10 countries per line.
  const keys = Object.keys(countries);
  const rows = [];
  for (let i = 0; i < keys.length; i += PER_LINE) {
    rows.push(keys.slice(i, i + PER_LINE).map((k) => `${JSON.stringify(k)}:${JSON.stringify(countries[k])}`).join(','));
  }
  const head = Object.entries(meta).map(([k, v]) => `${JSON.stringify(k)}:${JSON.stringify(v)}`).join(',\n');
  const json = `{\n${head},\n"countries":{\n${rows.join(',\n')}\n}}\n`;
  JSON.parse(json);            // must round-trip

  await mkdir(dirname(OUT), { recursive: true });
  await writeFile(OUT, json, 'utf8');

  const years = {};
  for (const v of Object.values(countries)) years[v.year] = (years[v.year] || 0) + 1;
  console.log(`Wrote ${OUT}`);
  console.log(`  ${n} countries, ${Buffer.byteLength(json)} bytes; World ${world.year}: ${world.g} g/kWh`);
  console.log(`  latest-year distribution: ${JSON.stringify(years)}`);
  if (rejected.size) console.log(`  out-of-range values ignored: ${[...rejected].map(([k, v]) => `${k} [${v.slice(-3).join('; ')}]${countries[k] ? ` -> used ${countries[k].year}` : ' -> country omitted'}`).join(' | ')}`);
  if (stale.length) console.log(`  omitted, no data since ${MIN_YEAR}: ${stale.join(', ')}`);
  if (unmapped.size) console.log(`  no ISO alpha-2 mapping (skipped): ${[...unmapped].join(', ')}`);
}

main().catch((e) => fail(e.stack || e.message));
