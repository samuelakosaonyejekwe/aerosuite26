// Materials catalogue: `node tools/build-materials.mjs [--check]`
// Rebuilds js/data/materials-catalogue.json from the transcribed tables in tools/materials-data/ and validates every
// item. Each number in those tables was read in the cited document (table, PDF page and column are kept beside it), in
// the units the document prints; this script only converts units, assembles the items and refuses to write anything
// that fails a check. `--check` builds and validates without writing and reports whether the stored file is current.
//
//   mil-hdbk-5j-design-properties.csv   one column of a MIL-HDBK-5J design-property table per row
//   mil-hdbk-5j-fracture-toughness.csv  the handbook's plane-strain toughness tables; "assigned_item" names the entry
//   mil-hdbk-17-2f-lamina.csv           room-temperature lamina means per material system, with the table of each value
//   mil-c-7438g-cores.csv               aluminium honeycomb core rows
//   records.csv                         one property per row for everything else (item, key, printed value, unit, source)
//   sources.json                        the documents

import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('..', import.meta.url)), DATA = root + 'tools/materials-data/', OUT = root + 'js/data/materials-catalogue.json';
export const RETRIEVED = '2026-10-09';

// ---- units ---------------------------------------------------------------------------------------
/** SI value per printed unit. */
export const UNIT = {
  ksi: 6.894757e6, Msi: 6.894757e9, psi: 6894.757, MPa: 1e6, 'ksi*sqrt(in)': 1.098843e6, 'lb/in3': 27679.905, pcf: 16.018463, 'g/cm3': 1000,
  'specific gravity 60/60F': 999.016, 'Btu/(lb F)': 4186.8, 'Btu ft/(h ft2 F)': 1.730735, '1e-6/F': 1.8e-6, in: 0.0254, 'MJ/kg': 1e6, 'kJ/m2': 1000,
  'in-lb/in2': 175.1268, '-': 1, degC: 1, V: 1, 'Wh/kg': 1, 'mOhm*Ah': 1,
};
export function toSI(value, unit) { if (!(unit in UNIT)) throw new Error(`unknown unit "${unit}"`); return Number(value) * UNIT[unit]; }
/** Paris/Walker coefficient from in/cycle with ΔK in ksi√in to m/cycle with ΔK in MPa√m. */
export const parisCtoSI = (C, n) => (Number(C) * 0.0254) / 1.098843 ** Number(n);
/** Round to a number of significant figures (conversions must not add false precision). */
export const sig = (x, n = 4) => (x === 0 ? 0 : Number(Number(x).toPrecision(n)));

// ---- reading -------------------------------------------------------------------------------------
export function parseCsv(text) {
  const rows = []; let row = [], cur = '', q = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (q) { if (ch === '"') { if (text[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += ch; }
    else if (ch === '"') q = true;
    else if (ch === ',') { row.push(cur); cur = ''; }
    else if (ch === '\n') { row.push(cur); rows.push(row); row = []; cur = ''; }
    else if (ch !== '\r') cur += ch;
  }
  if (cur !== '' || row.length) { row.push(cur); rows.push(row); }
  const head = rows.shift();
  return rows.filter((r) => r.length > 1).map((r) => Object.fromEntries(head.map((h, i) => [h, r[i] ?? ''])));
}
const csv = (f) => parseCsv(readFileSync(DATA + f, 'utf8'));
const num = (s) => (s === '' || s == null ? null : Number(s));

// ---- what the app expects per class (js/core/matlib.js FIELDS; the keys the app must estimate when absent) -----
export const KEYS = {
  metals: ['E', 'G', 'nu', 'rho', 'Sy', 'Su', 'KIc', 'parisC', 'parisM', 'dKth', 'sf', 'b', 'ef', 'c', 'k', 'cp', 'alpha'],
  plies: ['E1', 'E2', 'G12', 'nu12', 'rho', 't', 'Xt', 'Xc', 'Yt', 'Yc', 'S', 'GIc', 'GIIc', 'a1', 'a2'],
  fluids: ['rho', 'mu', 'cp', 'k'],
  batteries: ['wh_kg', 'v_nom', 'v_max', 'v_min', 'r_mohm_ah', 'cp', 'cycles_80'],
  cores: ['rho', 'cell_size', 'foil_t', 'Sc_stab', 'tau_L', 'tau_W'],
};
const OPTIONAL = { metals: ['Sy_A', 'Su_A', 'JC'], plies: [], fluids: ['LHV', 'co2_per_kg', 'lifecycle_factor', 'flash_C', 'bulk'], batteries: [], cores: [] };
/** The least an item must carry to be listed at all. */
const CORE = {
  metals: (m) => m.rho > 0 && m.E > 0 && m.Sy > 0 && m.Su > 0,
  plies: (m) => m.E1 > 0 && m.E2 > 0 && ['Xt', 'Xc', 'Yt', 'Yc', 'S'].some((k) => m[k] > 0),
  fluids: (m) => m.rho > 0 && ['mu', 'k', 'cp', 'LHV'].some((k) => m[k] > 0),
  batteries: (m) => m.v_min > 0 && m.v_max > 0 && m.wh_kg > 0,
  cores: (m) => m.rho > 0 && m.Sc_stab > 0,
};

// ---- assembly ------------------------------------------------------------------------------------
export function build() {
  const sources = JSON.parse(readFileSync(DATA + 'sources.json', 'utf8'));
  const items = new Map();
  const put = (it) => { if (items.has(it.name)) throw new Error(`duplicate item "${it.name}"`); items.set(it.name, { src: {}, locator: {}, grade: {}, ...it }); return items.get(it.name); };
  const set = (it, key, value, src, locator, grade) => {
    if (!Number.isFinite(value)) throw new Error(`${it.name}: ${key} is not a number`);
    if (key.startsWith('JC.')) { it.JC = { ...(it.JC || {}), [key.slice(3)]: value }; key = 'JC'; } else it[key] = value;
    it.src[key] = src; it.locator[key] = locator; it.grade[key] = grade;
  };
  const notes = new Map(), note = (it, s) => { if (s) notes.set(it.name, [...(notes.get(it.name) || []), s]); };

  // metals: MIL-HDBK-5J design properties
  for (const r of csv('mil-hdbk-5j-design-properties.csv')) {
    const S = 'mil-hdbk-5j', loc = `Table ${r.table} (PDF p. ${r.pdf_page}), column ${r.column}`;
    const form = [r.form, `condition ${r.condition}`, r.dimensions, r.specification, `${r.basis}-basis`].filter(Boolean).join('; ');
    const it = put({ name: r.name, class: 'metals', family: r.family, form, basis: r.basis });
    const g = 'handbook';
    set(it, 'E', sig(toSI(r.E_1e3ksi, 'Msi')), S, `${loc}, row ${r.E_row === 'E' ? 'E' : 'E, ' + r.E_row}`, g);
    if (r.G_1e3ksi) set(it, 'G', sig(toSI(r.G_1e3ksi, 'Msi')), S, `${loc}, row G`, g);
    if (r.mu) set(it, 'nu', Number(r.mu), S, `${loc}, row µ`, g);
    set(it, 'rho', sig(toSI(r.density_lb_in3, 'lb/in3')), S, `${loc}, row ω`, g);
    const dir = (d) => (d ? `, ${d}` : '');
    set(it, 'Sy_A', sig(toSI(r.Fty_ksi, 'ksi')), S, `${loc}, row Fty${dir(r.Fty_dir)}`, g);
    set(it, 'Su_A', sig(toSI(r.Ftu_ksi, 'ksi')), S, `${loc}, row Ftu${dir(r.Ftu_dir)}`, g);
    set(it, 'Sy', it.Sy_A, S, it.locator.Sy_A, g); set(it, 'Su', it.Su_A, S, it.locator.Su_A, g);
    it.allow_ref = `MIL-HDBK-5J ${r.basis}-basis, ${r.form.toLowerCase()}${r.dimensions ? ', ' + r.dimensions : ''}, Table ${r.table}`;
    if (r.C_btu_lb_F) set(it, 'cp', sig(toSI(r.C_btu_lb_F, 'Btu/(lb F)')), S, `${loc}, row C${r.C_temp ? ' (' + r.C_temp + ')' : ''}`, g);
    if (r.K_btu_hr_ft_F) set(it, 'k', sig(toSI(r.K_btu_hr_ft_F, 'Btu ft/(h ft2 F)')), S, `${loc}, row K${r.K_temp ? ' (' + r.K_temp + ')' : ''}`, g);
    if (r['alpha_1e-6_F']) set(it, 'alpha', sig(toSI(r['alpha_1e-6_F'], '1e-6/F')), S, `${loc}, row α${r.alpha_temp ? ' (' + r.alpha_temp + ')' : ''}`, g);
    note(it, `Sy and Su repeat the handbook ${r.basis}-basis design values (Fty, Ftu${r.Ftu_dir ? ', ' + r.Ftu_dir + ' direction' : ''}): the handbook tabulates no typical strengths.`);
    if (r.E_row === 'Primary') note(it, 'E is the primary modulus of the clad product.');
    if (r.E_row === 'L') note(it, 'E is the longitudinal value; the table lists L and LT separately.');
    if (r.C_temp || r.K_temp || r.alpha_temp) note(it, `Thermal values as printed: ${[r.C_temp && 'C ' + r.C_temp, r.K_temp && 'K ' + r.K_temp, r.alpha_temp && 'α ' + r.alpha_temp].filter(Boolean).join('; ')}.`);
    if (r.footnotes) note(it, `Table notes: ${r.footnotes}`);
  }
  // metals: MIL-HDBK-5J plane-strain fracture toughness (averages "for information only")
  for (const r of csv('mil-hdbk-5j-fracture-toughness.csv')) {
    if (!r.assigned_item) continue;
    const it = items.get(r.assigned_item); if (!it) throw new Error(`fracture toughness assigned to unknown item "${r.assigned_item}"`);
    set(it, 'KIc', sig(toSI(r.KIc_avg_ksi_sqrt_in, 'ksi*sqrt(in)')), 'mil-hdbk-5j', `Table ${r.table} (PDF p. ${r.pdf_page}), ${r.alloy.replace(/^Al /, '')} ${r.condition}, ${r.product_form}, ${r.orientation}, average`, 'handbook');
    note(it, `KIc is the handbook average for ${r.product_form.toLowerCase()} (${r.orientation}, product thickness ${r.product_thickness_in} in, ${r.sample_size} specimens, range ${r.KIc_min_ksi_sqrt_in}–${r.KIc_max_ksi_sqrt_in} ksi√in), given there for information only.`);
  }
  // plies: MIL-HDBK-17-2F
  for (const r of csv('mil-hdbk-17-2f-lamina.csv')) {
    const S = 'mil-hdbk-17-2f', g = 'handbook';
    const it = { name: `${r.material} (MIL-HDBK-17-2F)`, class: 'plies', family: r.family, src: {}, locator: {}, grade: {}, basis: 'typical' };
    it.form = [r.form || 'lamina', r.fiber_volume_pct && `fibre volume ${r.fiber_volume_pct} %`, r.normalized_by && `normalisation: ${r.normalized_by}`, `${r.test_temp_F} °F, ${r.moisture || 'ambient'}`].filter(Boolean).join('; ');
    const B = {}, cls = [];
    for (const [k, col, unit] of [['E1', 'E1_Msi', 'Msi'], ['E2', 'E2_Msi', 'Msi'], ['G12', 'G12_Msi', 'Msi'], ['nu12', 'nu12', '-'], ['Xt', 'Xt_ksi', 'ksi'], ['Xc', 'Xc_ksi', 'ksi'], ['Yt', 'Yt_ksi', 'ksi'], ['Yc', 'Yc_ksi', 'ksi'], ['S', 'S_ksi', 'ksi']]) {
      if (r[col] === '') continue;
      set(it, k, sig(toSI(r[col], unit)), S, `Table ${r[k + '_table']} (PDF p. ${r[k + '_pdf_page']}), first column (${r.test_temp_F} °F), Mean${/^normalized/.test(r[k + '_basis']) ? ', normalised' : ''}`, g);
      cls.push(`${k} ${r[k + '_basis']}`);
      if (r[k + '_Bvalue_ksi']) B[k] = sig(toSI(r[k + '_Bvalue_ksi'], 'ksi'));
    }
    if (r.density_g_cm3) set(it, 'rho', sig(toSI(r.density_g_cm3, 'g/cm3')), S, `Table ${r.header_table}, heading: composite density`, g);
    if (r.ply_thickness_in) set(it, 't', sig(toSI(r.ply_thickness_in, 'in')), S, `Table ${r.header_table}, heading: ${r.cpt_source}`, g);
    if (!CORE.plies(it)) continue; // E1, E2 and one strength are the least a ply entry needs
    put(it); Object.assign(items.get(it.name), it);
    if (Object.keys(B).length) items.get(it.name).B_basis = B;
    note(it, `Room-temperature means; data class and sample per value: ${cls.join('; ')}. Screening and interim data are not design allowables.`);
    if (Object.keys(B).length) note(it, 'B_basis holds the handbook B-values printed for the same columns.');
  }
  // cores: MIL-C-7438G
  for (const r of csv('mil-c-7438g-cores.csv')) {
    const S = 'mil-c-7438g', loc = `Table ${r.table} (PDF p. ${r.pdf_page}), row ${r.nominal_density_pcf} PCF, ${r.cell_size_in} in cell, ${r.foil_thickness_in} in foil`, g = 'handbook';
    const [a, b] = r.cell_size_in.split('/').map(Number);
    const it = put({ name: r.name, class: 'cores', family: 'honeycomb core', form: `${r.foil_alloy} aluminium foil honeycomb, hexagonal cells, room temperature`, basis: 'S' });
    set(it, 'rho', sig(toSI(r.nominal_density_pcf, 'pcf')), S, loc, g); set(it, 'cell_size', sig(toSI(a / b, 'in')), S, loc, g); set(it, 'foil_t', sig(toSI(r.foil_thickness_in, 'in')), S, loc, g);
    set(it, 'Sc_stab', sig(toSI(r.stabilized_compressive_strength_psi, 'psi')), S, loc, g); set(it, 'tau_L', sig(toSI(r.plate_shear_strength_L_psi, 'psi')), S, loc, g); set(it, 'tau_W', sig(toSI(r.plate_shear_strength_W_psi, 'psi')), S, loc, g);
    note(it, 'Specification values for room temperature: stabilised flatwise compressive strength and plate shear strength in the ribbon (L) and transverse (W) directions. The specification states that the values are for test purposes only and are not design allowables. No moduli are tabulated.');
  }
  // everything else, one property per row
  const recs = csv('records.csv'), nOf = (item, src) => recs.find((r) => r.item === item && r.source === src && r.key === 'parisM');
  for (const r of recs) {
    let it = items.get(r.item);
    if (!it) { if (!r.family) throw new Error(`records.csv: "${r.item}" is not an item and has no family`); it = put({ name: r.item, class: r.class, family: r.family, form: r.form, basis: 'typical' }); }
    if (it.class !== r.class) throw new Error(`records.csv: ${r.item} is ${it.class}, not ${r.class}`);
    const v = r.key === 'parisC' ? sig(parisCtoSI(r.printed, nOf(r.item, r.source).printed)) : sig(toSI(r.printed, r.unit));
    if (r.key in it && r.key !== 'KIc') throw new Error(`records.csv: ${r.item} already has ${r.key}`);
    if (r.key === 'KIc' && 'KIc' in it) continue; // the handbook test average is kept
    set(it, r.key, v, r.source, r.locator, r.grade);
    if (r.note) note(it, `${r.key.replace(/\..*/, '')}: ${r.note}`);
  }

  const out = [];
  for (const it of items.values()) {
    const all = [...KEYS[it.class], ...OPTIONAL[it.class]], sourced = all.filter((k) => k in it), missing = KEYS[it.class].filter((k) => !(k in it));
    it.coverage = { sourced, missing };
    it.note = [...new Set(notes.get(it.name) || [])].join(' ');
    const { src, locator, grade, coverage, note: n, ...rest } = it;
    out.push({ ...rest, src, locator, grade, coverage, note: n });
  }
  const cat = { generated: RETRIEVED, about: 'Sourced materials catalogue. Every number was read in the cited document and converted to SI; properties a source does not give are absent (listed under coverage.missing) and are never filled in here. Strengths of metals are handbook design values on the stated basis, repeated under Sy/Su because the handbook gives no typical values. Regenerate with `node tools/build-materials.mjs`.', sources: sources.map(({ file, ...s }) => s), items: out };
  validate(cat);
  return cat;
}

// ---- validation ----------------------------------------------------------------------------------
const GPa = 1e9, MPa = 1e6;
/** Plausible ranges per family: [E min, E max] in GPa, [rho min, rho max] in kg/m³. */
export const RANGES = {
  aluminium: { E: [66, 80], rho: [2600, 2900] }, 'aluminium (cast)': { E: [68, 76], rho: [2600, 2850] },
  magnesium: { E: [42, 47], rho: [1740, 1850] }, 'magnesium (cast)': { E: [42, 47], rho: [1740, 1850] },
  titanium: { E: [95, 125], rho: [4300, 4850] },
  'carbon steel': { E: [190, 210], rho: [7700, 8000] }, 'low-alloy steel': { E: [190, 210], rho: [7700, 8000] }, 'intermediate-alloy steel': { E: [190, 210], rho: [7700, 8000] }, 'high-alloy steel': { E: [180, 210], rho: [7700, 8000] },
  'stainless / PH steel': { E: [175, 205], rho: [7600, 8000] }, 'iron-base superalloy': { E: [195, 210], rho: [7900, 8350] },
  'nickel / cobalt': { E: [195, 245], rho: [8000, 9200] }, beryllium: { E: [280, 300], rho: [1800, 1900] }, copper: { E: [95, 135], rho: [7800, 8400] },
};
export function validate(cat) {
  const bad = [], ids = new Set(cat.sources.map((s) => s.id)), names = new Set();
  const say = (it, m) => bad.push(`${it.name}: ${m}`), within = (it, k, v, lo, hi) => { if (v != null && !Number.isNaN(v) && !(v >= lo && v <= hi)) say(it, `${k} = ${v} outside ${lo}…${hi}`); };
  for (const s of cat.sources) if (!s.id || !s.title || !s.org || !/^https:\/\//.test(s.url) || !s.licence_note) bad.push(`source ${s.id}: incomplete`);
  for (const it of cat.items) {
    if (names.has(it.name)) say(it, 'duplicate name'); names.add(it.name);
    if (!KEYS[it.class]) { say(it, `unknown class ${it.class}`); continue; }
    if (!it.family || !it.form || !it.note) say(it, 'family, form and note are required');
    if (!CORE[it.class](it)) say(it, 'core properties are missing');
    const known = new Set([...KEYS[it.class], ...OPTIONAL[it.class]]);
    for (const [k, v] of Object.entries(it)) {
      if (['name', 'class', 'family', 'form', 'basis', 'allow_ref', 'src', 'locator', 'grade', 'coverage', 'note', 'B_basis'].includes(k)) continue;
      if (!known.has(k)) say(it, `unexpected key ${k}`);
      if (k === 'JC') { for (const x of Object.values(v)) if (!Number.isFinite(x) || x <= 0) say(it, 'JC constants must be positive numbers'); }
      else if (!Number.isFinite(v)) say(it, `${k} is not a number`);
      else if (v <= 0 && k !== 'flash_C') say(it, `${k} must be positive`);
      if (!ids.has(it.src[k])) say(it, `${k} has no known source id`);
      if (!it.locator[k]) say(it, `${k} has no locator`);
      if (!['handbook', 'primary-report', 'datasheet', 'secondary'].includes(it.grade[k])) say(it, `${k} has no grade`);
    }
    for (const k of Object.keys(it.src)) if (!(k in it)) say(it, `source given for absent key ${k}`);
    if (it.coverage.sourced.some((k) => !(k in it)) || it.coverage.missing.some((k) => k in it) || KEYS[it.class].some((k) => !(k in it) !== it.coverage.missing.includes(k))) say(it, 'coverage does not match the item');
    if (it.class === 'metals') {
      const R = RANGES[it.family]; if (!R) { say(it, `no range for family ${it.family}`); continue; }
      within(it, 'E', it.E / GPa, ...R.E); within(it, 'rho', it.rho, ...R.rho); within(it, 'nu', it.nu, 0.05, 0.4); within(it, 'G', it.G / GPa, 0.3 * it.E / GPa, 0.5 * it.E / GPa);
      within(it, 'Su', it.Su / MPa, 50, 2500); within(it, 'Sy', it.Sy / MPa, 30, 2300);
      if (it.Sy > it.Su) say(it, 'Sy exceeds Su'); if (it.Sy_A > it.Sy || it.Su_A > it.Su) say(it, 'design value above Sy/Su');
      if (!['A', 'B', 'S', 'typical'].includes(it.basis)) say(it, 'basis'); if (it.Sy_A && !it.allow_ref) say(it, 'allow_ref missing');
      within(it, 'KIc', it.KIc / 1e6, 10, 250); within(it, 'cp', it.cp, 350, 2000); within(it, 'k', it.k, 5, 250); within(it, 'alpha', it.alpha * 1e6, 4, 30);
      within(it, 'parisM', it.parisM, 2, 5); within(it, 'parisC', it.parisC, 1e-13, 1e-9); if (('parisC' in it) !== ('parisM' in it)) say(it, 'parisC and parisM go together');
    }
    if (it.class === 'plies') {
      within(it, 'E1', it.E1 / GPa, 20, 350); within(it, 'E2', it.E2 / GPa, 5, 80); if (it.E2 > 1.15 * it.E1) say(it, 'E2 above E1'); within(it, 'G12', it.G12 / GPa, 2, 10); within(it, 'nu12', it.nu12, 0.01, 0.5);
      for (const k of ['Xt', 'Xc', 'Yt', 'Yc', 'S']) within(it, k, it[k] / MPa, 20, 3000);
      within(it, 't', it.t * 1e3, 0.05, 0.6); within(it, 'rho', it.rho, 1300, 2100); within(it, 'GIc', it.GIc, 50, 2000); within(it, 'GIIc', it.GIIc, 100, 5000);
      for (const [k, v] of Object.entries(it.B_basis || {})) if (!(v > 0 && v <= it[k])) say(it, `B-basis ${k} above the mean`);
    }
    if (it.class === 'fluids') { within(it, 'rho', it.rho, 600, 900); within(it, 'LHV', it.LHV / 1e6, 40, 46); within(it, 'flash_C', it.flash_C, 30, 110); }
    if (it.class === 'batteries') { if (!(it.v_min < it.v_nom && it.v_nom < it.v_max)) say(it, 'voltages out of order'); within(it, 'wh_kg', it.wh_kg, 50, 450); within(it, 'r_mohm_ah', it.r_mohm_ah, 5, 300); }
    if (it.class === 'cores') { within(it, 'rho', it.rho, 10, 250); within(it, 'cell_size', it.cell_size * 1e3, 1, 15); within(it, 'foil_t', it.foil_t * 1e6, 10, 200); within(it, 'Sc_stab', it.Sc_stab / MPa, 0.05, 20); if (it.tau_W > it.tau_L) say(it, 'W shear above L shear'); }
  }
  if (bad.length) throw new Error(`materials catalogue: ${bad.length} violation(s)\n  ${bad.join('\n  ')}`);
  return true;
}
export const serialise = (cat) => JSON.stringify(cat, null, 1) + '\n';

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) {
  let cat;
  try { cat = build(); } catch (e) { console.error(String(e.message || e)); process.exit(1); }
  const text = serialise(cat), count = {}; for (const it of cat.items) count[it.class] = (count[it.class] || 0) + 1;
  if (process.argv.includes('--check')) {
    let current = false; try { current = readFileSync(OUT, 'utf8') === text; } catch { /* not built yet */ }
    console.log(`${cat.items.length} items valid (${Object.entries(count).map(([k, v]) => `${k} ${v}`).join(', ')}); stored file is ${current ? 'current' : 'OUT OF DATE'}`);
    process.exit(current ? 0 : 1);
  }
  writeFileSync(OUT, text);
  console.log(`wrote js/data/materials-catalogue.json: ${cat.items.length} items (${Object.entries(count).map(([k, v]) => `${k} ${v}`).join(', ')}), ${cat.sources.length} sources`);
}
