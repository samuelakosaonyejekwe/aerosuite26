// Materials catalogue checks: `node tests/materials.mjs`
//  - unit conversions used by the build;
//  - js/data/materials-catalogue.json: schema, per-family ranges, unique names, a source id, locator and grade for
//    every number, coverage lists, and that it is exactly what tools/build-materials.mjs produces from tools/materials-data/;
//  - spot checks against the source documents themselves: when MATERIALS_DOCS names a directory holding the downloaded
//    PDFs (paths as in tools/materials-data/sources.json, "file") and python3 with PyMuPDF is available, the cited pages are read again as plain
//    text — a different route from the table extraction used to transcribe them — and the printed numbers are looked
//    up there. Without the documents that part is skipped and says so.
// No network is used.

import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { spawnSync } from 'node:child_process';
import { build, serialise, validate, toSI, parisCtoSI, sig, parseCsv, UNIT, KEYS, RANGES } from '../tools/build-materials.mjs';

const root = fileURLToPath(new URL('..', import.meta.url)), text = (f) => readFileSync(root + f, 'utf8');
let fail = 0, n = 0;
const ok = (cond, msg) => { n++; if (!cond) { fail++; console.log('  FAIL ' + msg); } };
const near = (a, b, rel = 1e-6) => Math.abs(a - b) <= rel * Math.abs(b);
const section = (t) => console.log(t);

// ---- conversions ---------------------------------------------------------------------------------
section('Unit conversions');
ok(near(toSI(1, 'ksi'), 6.894757e6) && near(toSI(64, 'ksi'), 441.264448e6) && near(toSI(10.5, 'Msi'), 72.3949485e9), 'ksi and Msi to Pa');
ok(near(toSI(0.1, 'lb/in3'), 2767.9905) && near(toSI(0.283, 'lb/in3'), 7833.413), 'lb/in³ to kg/m³');
ok(near(toSI(0.23, 'Btu/(lb F)'), 962.964) && near(toSI(91, 'Btu ft/(h ft2 F)'), 157.496885) && near(toSI(12.8, '1e-6/F'), 23.04e-6), 'specific heat, conductivity and expansion');
ok(near(toSI(30, 'ksi*sqrt(in)'), 32.96529e6) && near(UNIT['ksi*sqrt(in)'], 6.894757e6 * Math.sqrt(0.0254), 1e-6), 'ksi√in to Pa√m, consistent with ksi and inch');
ok(near(toSI(4.5, 'pcf'), 72.0830835) && near(UNIT.pcf, 0.45359237 / 0.3048 ** 3, 1e-6) && near(toSI(405, 'psi'), 2.792376585e6) && near(toSI(0.0053, 'in'), 1.3462e-4), 'pcf, psi and inch');
ok(near(toSI(4.22, 'in-lb/in2'), 739.035) && near(UNIT['in-lb/in2'], 4.4482216 / 0.0254, 1e-6) && near(toSI(0.221, 'kJ/m2'), 221) && near(toSI(0.816, 'specific gravity 60/60F'), 815.197), 'fracture energy and specific gravity');
ok(near(parisCtoSI(1.3e-8, 2.791), (1.3e-8 * 0.0254) / 1.098843 ** 2.791) && near(parisCtoSI(1e-9, 3) * 1.098843 ** 3, 2.54e-11), 'crack-growth coefficient: in/cycle, ksi√in to m/cycle, MPa√m');
{ let threw = false; try { toSI(1, 'furlong'); } catch { threw = true; } ok(threw && sig(72394948500) === 72390000000 && sig(0.000134620) === 0.0001346, 'unknown units are refused; values keep four significant figures'); }
ok(parseCsv('a,b\n"x, y",2\n"say ""hi""",3\n').length === 2 && parseCsv('a,b\n"x, y",2\n')[0].a === 'x, y' && parseCsv('a,b\n"say ""hi""",3\n')[0].a === 'say "hi"', 'CSV reader handles quoted commas and quotes');

// ---- the stored catalogue ------------------------------------------------------------------------
section('Catalogue file');
const cat = JSON.parse(text('js/data/materials-catalogue.json')), items = cat.items, byName = new Map(items.map((i) => [i.name, i])), ids = new Set(cat.sources.map((s) => s.id));
let rebuilt = null; try { rebuilt = build(); } catch (e) { ok(false, 'build failed: ' + e.message); }
ok(rebuilt && serialise(rebuilt) === text('js/data/materials-catalogue.json'), 'the stored file is what tools/build-materials.mjs produces (run `node tools/build-materials.mjs`)');
{ let v = false; try { v = validate(cat); } catch (e) { console.log(e.message.split('\n').slice(0, 6).join('\n')); } ok(v, 'the build validation accepts the stored file'); }
ok(/^\d{4}-\d\d-\d\d$/.test(cat.generated) && cat.sources.length >= 8 && cat.sources.every((s) => s.id && s.title && s.org && /^https:\/\//.test(s.url) && s.licence_note.length > 30) && ids.size === cat.sources.length, 'header and source list');
ok(byName.size === items.length, 'no duplicate names');
const count = (c) => items.filter((i) => i.class === c).length;
ok(count('metals') >= 150 && count('plies') >= 20 && count('cores') >= 10 && count('fluids') >= 5 && count('batteries') >= 1, `counts: metals ${count('metals')}, plies ${count('plies')}, cores ${count('cores')}, fluids ${count('fluids')}, batteries ${count('batteries')}`);
const META = ['name', 'class', 'family', 'form', 'basis', 'allow_ref', 'src', 'locator', 'grade', 'coverage', 'note', 'B_basis'];
for (const it of items) {
  const props = Object.keys(it).filter((k) => !META.includes(k));
  ok(typeof it.name === 'string' && KEYS[it.class] && it.family && it.form && typeof it.note === 'string' && it.note.length > 10 && ['A', 'B', 'S', 'typical'].includes(it.basis), `${it.name}: name, class, family, form, basis, note`);
  ok(props.length > 0 && props.every((k) => ids.has(it.src[k]) && typeof it.locator[k] === 'string' && it.locator[k].length > 5 && ['handbook', 'primary-report', 'datasheet', 'secondary'].includes(it.grade[k])), `${it.name}: every property has a source id that exists, a locator and a grade`);
  ok(props.every((k) => (k === 'JC' ? Object.values(it.JC).every(Number.isFinite) : Number.isFinite(it[k]))), `${it.name}: properties are numbers`);
  ok(Object.keys(it.src).every((k) => props.includes(k)) && it.coverage.sourced.every((k) => props.includes(k)) && it.coverage.missing.every((k) => !props.includes(k)) && KEYS[it.class].every((k) => props.includes(k) || it.coverage.missing.includes(k)), `${it.name}: coverage lists what is sourced and what the app must estimate`);
  if (it.class === 'metals') {
    const R = RANGES[it.family];
    ok(R && it.E >= R.E[0] * 1e9 && it.E <= R.E[1] * 1e9 && it.rho >= R.rho[0] && it.rho <= R.rho[1], `${it.name}: E ${(it.E / 1e9).toFixed(1)} GPa and density ${it.rho} within the ${it.family} range`);
    ok(it.Sy > 0 && it.Sy <= it.Su && it.Su < 2.5e9 && (!it.Sy_A || (it.Sy_A <= it.Sy && it.Su_A <= it.Su && /MIL-HDBK-5J [ABS]-basis.*Table/.test(it.allow_ref))) && (it.nu == null || (it.nu > 0.05 && it.nu < 0.4)), `${it.name}: strengths ordered, allowable named, Poisson ratio sane`);
  }
  if (it.class === 'plies') ok(it.E1 > 2e10 && it.E2 > 0 && it.E2 <= 1.15 * it.E1 && ['Xt', 'Xc', 'Yt', 'Yc', 'S'].some((k) => it[k] > 0) && (it.t == null || (it.t > 5e-5 && it.t < 6e-4)), `${it.name}: ply moduli, a strength and thickness`);
  if (it.class === 'batteries') ok(it.v_min < it.v_nom && it.v_nom < it.v_max && it.wh_kg > 50 && it.wh_kg < 450, `${it.name}: voltage window and specific energy`);
  if (it.class === 'fluids') ok(it.rho > 600 && it.rho < 900 && ['mu', 'k', 'cp', 'LHV'].some((k) => it[k] > 0), `${it.name}: density and one transport or energy property`);
  if (it.class === 'cores') ok(it.rho > 10 && it.rho < 250 && it.tau_W <= it.tau_L && it.Sc_stab > 0, `${it.name}: core density and strengths`);
}
// a handful of values worked by hand from the printed figures
const P = (name) => byName.get(name) || {};
ok(P('Al 2024-T3 sheet').Su_A === sig(64 * 6.894757e6) && P('Al 2024-T3 sheet').Sy_A === sig(47 * 6.894757e6) && P('Al 2024-T3 sheet').rho === 2768 && P('Al 2024-T3 sheet').nu === 0.33, 'Al 2024-T3 sheet: 64 / 47 ksi A-basis, 0.100 lb/in³');
ok(P('Al 7075-T651 plate').KIc === sig(26 * 1.098843e6) && P('Al 7075-T651 plate').parisM === 2.791 && near(P('Al 7075-T651 plate').parisC, 2.538e-10, 1e-3), 'Al 7075-T651 plate: KIc 26 ksi√in, Walker fit at R = 0');
ok(P('Al 2024-T3 sheet').JC?.A === 369e6 && P('Ti-6Al-4V Annealed sheet').JC?.B === 1092e6 && P('Ti-6Al-4V Annealed sheet').JC?.n === 0.93, 'Johnson–Cook constants of 2024-T3 and Ti-6Al-4V');
ok(P('IM7/8552 unidirectional tape (NCAMP, 190 gsm, 35% RC)').GIc === 221 && near(P('IM7/8552 unidirectional tape (NCAMP, 190 gsm, 35% RC)').E1, 22.99 * 6.894757e9, 1e-3), 'IM7/8552: E1 22.99 Msi, GIc 0.221 kJ/m²');
ok(P('Al 5052 honeycomb 1/8-0.0010 (4.5 pcf)').Sc_stab === sig(405 * 6894.757) && P('Al 5056 honeycomb 1/8-0.0020 (8.1 pcf)').tau_L === sig(740 * 6894.757), 'honeycomb core rows');

// ---- intermediate tables -------------------------------------------------------------------------
section('Transcribed tables');
const csv = (f) => parseCsv(text('tools/materials-data/' + f)), h5 = csv('mil-hdbk-5j-design-properties.csv'), recs = csv('records.csv');
const srcFiles = JSON.parse(text('tools/materials-data/sources.json'));
ok(h5.length >= 150 && h5.every((r) => r.name && r.table && Number(r.pdf_page) > 0 && /^\d+ of \d+$/.test(r.column) && ['A', 'S'].includes(r.basis) && Number(r.Ftu_ksi) >= Number(r.Fty_ksi) && r.verified_against_plain_text), 'handbook rows carry table, page, column, basis and the cross-check that passed');
ok(recs.every((r) => ids.has(r.source) && r.locator && r.key && r.printed !== '' && (r.key === 'parisC' || r.unit in UNIT)) && srcFiles.every((s) => s.file), 'records name a known source, a locator and a known unit');
for (const r of h5.filter((_, i) => i % 13 === 0)) { const it = P(r.name); ok(near(it.Su_A, r.Ftu_ksi * 6.894757e6, 1e-3) && near(it.Sy_A, r.Fty_ksi * 6.894757e6, 1e-3) && near(it.E, r.E_1e3ksi * 6.894757e9, 1e-3) && near(it.rho, r.density_lb_in3 * 27679.905, 1e-3), `${r.name}: catalogue values are the conversions of the printed ones`); }

// ---- spot checks against the documents -----------------------------------------------------------
section('Spot checks against the source documents');
const DOCS = process.env.MATERIALS_DOCS || '';
const fileOf = Object.fromEntries(srcFiles.map((s) => [s.id, `${DOCS}/${s.file}`]));
/** Plain text of the given pages, read again from the PDF. Returns null when the document or the reader is missing. */
function pages(id, list) {
  if (!DOCS || !existsSync(fileOf[id])) return null;
  const py = 'import sys, json, pymupdf\nd = pymupdf.open(sys.argv[1])\nprint(json.dumps({p: d[int(p) - 1].get_text() for p in sys.argv[2:]}))';
  const r = spawnSync('python3', ['-P', '-c', py, fileOf[id], ...list.map(String)], { encoding: 'utf8', maxBuffer: 64e6 });
  if (r.status !== 0) return null;
  try { return JSON.parse(r.stdout); } catch { return null; }
}
const tokens = (t) => t.replace(/\.\s\.\s\./g, '...').split(/\s+/).filter(Boolean).map((w) => (/^\d+(\.\d+)?[a-z]{1,2}$/.test(w) ? w.replace(/[a-z]+$/, '') : w));
const hasSeq = (hay, seq) => { for (let i = 0; i + seq.length <= hay.length; i++) { let j = 0; while (j < seq.length && hay[i + j] === seq[j]) j++; if (j === seq.length) return true; } return false; };
let spot = 0, skipped = [];
{ // MIL-HDBK-5J: the whole strength column of an entry, in table order, must appear in the page text
  const pick = []; const seen = new Set();
  for (const r of h5) { if (r.verified_against_plain_text !== 'column' || seen.has(r.table)) continue; seen.add(r.table); pick.push(r); }
  const sample = pick.filter((_, i) => i % 3 === 0).slice(0, 30), txt = pages('mil-hdbk-5j', [...new Set(sample.map((r) => r.pdf_page))]);
  if (!txt) skipped.push('MIL-HDBK-5J');
  else for (const r of sample) {
    const hay = tokens(txt[r.pdf_page]), hay2 = hay.filter((w) => !/^[a-z]$/.test(w)), seq = r.column_values.split(' ').map((v) => v.split('=')[1]).map((v) => (/^[. ]+$/.test(v) ? '...' : v));
    const it = P(r.name), ftu = r.Ftu_ksi, fty = r.Fty_ksi, agrees = near(it.Su_A, ftu * 6.894757e6, 1e-3) && near(it.Sy_A, fty * 6.894757e6, 1e-3);
    spot++; ok((hasSeq(hay, seq) || hasSeq(hay2, seq)) && agrees && seq.includes(ftu) && seq.includes(fty) && hay.includes(r.E_1e3ksi) && hay.includes(r.density_lb_in3), `${r.name}: Table ${r.table} p. ${r.pdf_page}: column ${seq.slice(0, 6).join(' ')} …, E ${r.E_1e3ksi}, ω ${r.density_lb_in3} found in the page; catalogue ${ftu} / ${fty} ksi`);
  }
}
{ // other documents: the printed figures, in their printed neighbourhood
  const flat = (t) => t.replace(/\s+/g, ' ');
  const checks = [
    ['faa-ar-05-15', 83, /C = 0\.130E-07 n = 2\.791/, () => P('Al 7075-T651 plate').parisM === 2.791],
    ['faa-ar-05-15', 39, /C = 0\.167E-08 n = 3\.273/, () => P('Al 2024-T3 sheet').parisM === 3.273],
    ['faa-ar-05-15', 77, /C = 0\.600E-09 n = 3\.613/, () => P('Al 7050-T7451 plate').parisM === 3.613],
    ['faa-ar-05-15', 37, /K1c = 30\.0/, () => P('Al 2024-T3 sheet').KIc === sig(30 * 1.098843e6)],
    ['llnl-ucrl-id-134691', 25, /1098 1092 \.93 \.014 1\.1/, () => P('Ti-6Al-4V Annealed sheet').JC.A === 1098e6],
    ['llnl-ucrl-id-134691', 25, /369 684 \.73 \.0083 1\.7/, () => P('Al 2024-T3 sheet').JC.C === 0.0083],
    ['nasa-tm-2011-217059', 17, /Heat of Combustion \(MJ\/kg\) 43\.3 44\.4 43\.8 44\.1 43\.8/, () => P('JP-8 (AAFEX test fuel)').LHV === 43.3e6],
    ['nasa-tm-2011-217059', 17, /Specific Gravity 0\.816 0\.738 0\.777 0\.763 0\.789/, () => near(P('Fischer-Tropsch SPK from natural gas, Shell (AAFEX FT1)').rho, 0.738 * 999.016, 1e-3)],
    ['nasa-tm-2011-217059', 17, /Flash Point, deg C 46 41 43 42 46/, () => P('JP-8 / FT1 50:50 blend (AAFEX)').flash_C === 43],
    ['nasa-tm-2010-216838', 19, /Average 4\.22 6\.50/, () => near(P('IM7/8552 unidirectional tape (NCAMP, 190 gsm, 35% RC)').GIIc, 4.22 * 175.1268, 1e-3)],
    ['nasa-ntrs-20120016494', 16, /MEAN 0\.237 0\.221 0\.206/, () => P('IM7/8552 unidirectional tape (NCAMP, 190 gsm, 35% RC)').GIc === 221],
    ['ncamp-ncp-rp-2009-028', 37, /\(22\.57\) \(22\.99\) \(24\.00\)/, () => near(P('IM7/8552 unidirectional tape (NCAMP, 190 gsm, 35% RC)').E1, 22.99 * 6.894757e9, 1e-3)],
    ['ncamp-ncp-rp-2009-028', 35, /Mean 324\.62 362\.69 248\.94 9\.29 41\.44/, () => near(P('IM7/8552 unidirectional tape (NCAMP, 190 gsm, 35% RC)').Xc, 248.94 * 6.894757e6, 1e-3)],
    ['molicel-inr21700-p42a', 1, /Gravimetric 230 Wh\/kg/, () => P('Molicel INR-21700-P42A (NMC/NCA 21700 cell)').wh_kg === 230],
    ['molicel-inr21700-p42a', 1, /Nominal 3\.6 V Charge 4\.2 V Discharge 2\.5 V/, () => P('Molicel INR-21700-P42A (NMC/NCA 21700 cell)').v_min === 2.5],
    ['mil-c-7438g', 16, /370 385 405 510/, () => P('Al 5052 honeycomb 1/8-0.0010 (4.5 pcf)').Sc_stab === sig(405 * 6894.757)],
    ['mil-c-7438g', 17, /687 720 740 740/, () => P('Al 5056 honeycomb 1/8-0.0020 (8.1 pcf)').tau_L === sig(740 * 6894.757)],
    ['mil-hdbk-5j', 315, /7075-T651 Plate L-T .{0,40}? 26 /, () => P('Al 7075-T651 plate').KIc === sig(26 * 1.098843e6)],
    ['mil-hdbk-5j', 313, /2024-T351 Plate L-T .{0,40}? 31 /, () => P('Al 2024-T351 plate').KIc === sig(31 * 1.098843e6)],
    ['mil-hdbk-17-2f', 226, /Mean 211 207 199 197 236 232/, () => near(P('T300 15k/976 unidirectional tape (MIL-HDBK-17-2F)').Xt, 211 * 6.894757e6, 1e-3)],
    ['mil-hdbk-17-2f', 226, /Mean 19\.6 19\.3 20\.8 20\.4 22\.6 22\.4/, () => near(P('T300 15k/976 unidirectional tape (MIL-HDBK-17-2F)').E1, 19.6 * 6.894757e9, 1e-3)],
    ['mil-hdbk-17-2f', 228, /Mean 5\.66 4\.73 3\.81 3\.47/, () => near(P('T300 15k/976 unidirectional tape (MIL-HDBK-17-2F)').Yt, 5.66 * 6.894757e6, 1e-3)],
  ];
  const cache = {};
  for (const [id, page, rx, agrees] of checks) {
    cache[id] ??= pages(id, [...new Set(checks.filter((c) => c[0] === id).map((c) => c[1]))]);
    if (!cache[id]) { if (!skipped.includes(id)) skipped.push(id); continue; }
    spot++; ok(rx.test(flat(cache[id][page])) && agrees(), `${id} p. ${page}: ${rx.source.replace(/\\/g, '')} is printed there and the catalogue holds its conversion`);
  }
}
if (skipped.length) console.log(`  (skipped: source documents or PyMuPDF not available for ${skipped.join(', ')}; set MATERIALS_DOCS to the directory holding ${srcFiles.map((s) => s.file).slice(0, 2).join(', ')} …)`);
else ok(spot >= 25, `${spot} values re-read in the documents`);
console.log(`  ${spot} values re-read in the source documents`);

console.log(`${n} checks, ${fail} failed`);
process.exit(fail ? 1 : 0);
