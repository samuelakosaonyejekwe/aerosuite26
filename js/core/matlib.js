// Materials library: one description of every property the suites use, plus the user's own materials.
// Built-in data live in js/data/materials.js. Materials a user adds are stored on the device and registered
// into the same tables, so every suite can select them exactly like a built-in one.

import { METALS, PLIES, FLUIDS, BATTERIES } from '../data/materials.js';

export const TABLES = { metals: METALS, plies: PLIES, fluids: FLUIDS, batteries: BATTERIES };
export const CLASS_LABEL = { metals: 'Metals', plies: 'Composite plies', fluids: 'Fuels & fluids', batteries: 'Battery cells' };
const BUILTIN = Object.fromEntries(Object.entries(TABLES).map(([k, t]) => [k, new Set(Object.keys(t))]));
export const isBuiltin = (cls, name) => BUILTIN[cls].has(name);
export const builtinNames = (cls) => [...BUILTIN[cls]];

// [key, label, display unit, SI-per-display-unit, group, required, help]
const F = (key, label, unit, scale, group, req = true, help = '') => ({ key, label, unit, scale, group, req, help });
export const FIELDS = {
  metals: [
    F('E', 'Young’s modulus', 'GPa', 1e9, 'Elastic'), F('G', 'Shear modulus', 'GPa', 1e9, 'Elastic'), F('nu', 'Poisson’s ratio', '-', 1, 'Elastic'), F('rho', 'Density', 'kg/m³', 1, 'Elastic'),
    F('Sy', 'Yield strength (typical)', 'MPa', 1e6, 'Strength'), F('Su', 'Ultimate strength (typical)', 'MPa', 1e6, 'Strength'),
    F('Sy_A', 'Design yield allowable', 'MPa', 1e6, 'Strength', false, 'Statistically based value (A- or B-basis). Used for margins of safety when given.'), F('Su_A', 'Design ultimate allowable', 'MPa', 1e6, 'Strength', false),
    F('KIc', 'Fracture toughness KIc', 'MPa√m', 1e6, 'Fracture'), F('parisC', 'Paris coefficient C', 'm/cycle·(MPa√m)⁻ᵐ', 1, 'Fracture', true, 'da/dN = C·ΔKᵐ with ΔK in MPa√m'), F('parisM', 'Paris exponent m', '-', 1, 'Fracture'), F('dKth', 'Threshold ΔKth', 'MPa√m', 1, 'Fracture'),
    F('sf', 'Fatigue strength coefficient σ′f', 'MPa', 1e6, 'Fatigue', true, 'Basquin: stress amplitude = σ′f·(2N)ᵇ'), F('b', 'Fatigue strength exponent b', '-', 1, 'Fatigue'), F('ef', 'Fatigue ductility coefficient ε′f', '-', 1, 'Fatigue'), F('c', 'Fatigue ductility exponent c', '-', 1, 'Fatigue'),
    F('k', 'Thermal conductivity', 'W/m/K', 1, 'Thermal'), F('cp', 'Specific heat', 'J/kg/K', 1, 'Thermal'), F('alpha', 'Thermal expansion', 'µm/m/K', 1e-6, 'Thermal'),
  ],
  plies: [
    F('E1', 'Fibre-direction modulus E1', 'GPa', 1e9, 'Elastic'), F('E2', 'Transverse modulus E2', 'GPa', 1e9, 'Elastic'), F('G12', 'In-plane shear modulus G12', 'GPa', 1e9, 'Elastic'), F('nu12', 'Major Poisson’s ratio', '-', 1, 'Elastic'), F('rho', 'Density', 'kg/m³', 1, 'Elastic'), F('t', 'Cured ply thickness', 'mm', 1e-3, 'Elastic'),
    F('Xt', 'Fibre tension strength Xt', 'MPa', 1e6, 'Strength'), F('Xc', 'Fibre compression strength Xc', 'MPa', 1e6, 'Strength'), F('Yt', 'Transverse tension strength Yt', 'MPa', 1e6, 'Strength'), F('Yc', 'Transverse compression strength Yc', 'MPa', 1e6, 'Strength'), F('S', 'In-plane shear strength S', 'MPa', 1e6, 'Strength'),
    F('GIc', 'Mode I toughness GIc', 'J/m²', 1, 'Fracture'), F('GIIc', 'Mode II toughness GIIc', 'J/m²', 1, 'Fracture'),
    F('a1', 'Expansion along fibres', 'µm/m/K', 1e-6, 'Thermal'), F('a2', 'Expansion across fibres', 'µm/m/K', 1e-6, 'Thermal'),
  ],
  fluids: [
    F('rho', 'Density', 'kg/m³', 1, 'Fluid'), F('mu', 'Dynamic viscosity', 'mPa·s', 1e-3, 'Fluid'), F('cp', 'Specific heat', 'J/kg/K', 1, 'Fluid'), F('k', 'Thermal conductivity', 'W/m/K', 1, 'Fluid'),
    F('LHV', 'Lower heating value', 'MJ/kg', 1e6, 'Fuel', false), F('co2_per_kg', 'CO₂ per kg burnt', 'kg/kg', 1, 'Fuel', false), F('lifecycle_factor', 'Life-cycle CO₂ factor', '-', 1, 'Fuel', false, 'Share of combustion CO₂ counted on a life-cycle basis'), F('flash_C', 'Flash point', '°C', 1, 'Fuel', false),
    F('bulk', 'Bulk modulus', 'GPa', 1e9, 'Hydraulic', false),
  ],
  batteries: [
    F('wh_kg', 'Specific energy', 'Wh/kg', 1, 'Cell'), F('v_nom', 'Nominal voltage', 'V', 1, 'Cell'), F('v_max', 'Maximum voltage', 'V', 1, 'Cell'), F('v_min', 'Minimum voltage', 'V', 1, 'Cell'),
    F('r_mohm_ah', 'Internal resistance × capacity', 'mΩ·Ah', 1, 'Cell'), F('cp', 'Specific heat', 'J/kg/K', 1, 'Cell'), F('cycles_80', 'Cycles to 80 % capacity', 'cycles', 1, 'Cell'),
  ],
};
/** Derived comparison indices, per class: [label, unit, fn]. */
export const INDICES = {
  metals: [['Specific stiffness E/ρ', 'MN·m/kg', (m) => m.E / m.rho / 1e6], ['Specific strength Su/ρ', 'kN·m/kg', (m) => m.Su / m.rho / 1e3], ['Yield-to-ultimate ratio', '-', (m) => m.Sy / m.Su], ['Thermal diffusivity', 'mm²/s', (m) => (m.k / (m.rho * m.cp)) * 1e6], ['Critical crack size at Sy/2 (through crack)', 'mm', (m) => ((m.KIc / (m.Sy / 2)) ** 2 / Math.PI) * 1e3]],
  plies: [['Specific stiffness E1/ρ', 'MN·m/kg', (m) => m.E1 / m.rho / 1e6], ['Specific strength Xt/ρ', 'kN·m/kg', (m) => m.Xt / m.rho / 1e3], ['Anisotropy E1/E2', '-', (m) => m.E1 / m.E2], ['Areal mass per ply', 'g/m²', (m) => m.rho * m.t * 1e3]],
  fluids: [['Kinematic viscosity', 'mm²/s', (m) => (m.mu / m.rho) * 1e6], ['Energy per litre', 'MJ/L', (m) => (m.LHV ? (m.LHV * m.rho) / 1e9 : NaN)]],
  batteries: [['Energy per cell-volt', 'Ah/kg', (m) => m.wh_kg / m.v_nom], ['Usable voltage window', 'V', (m) => m.v_max - m.v_min]],
};

// ---- user materials (kept on the device) --------------------------------------------------------
const KEY = 'aerosuite26.materials';
const blank = () => ({ metals: {}, plies: {}, fluids: {}, batteries: {} });
export function customMaterials() {
  try { const v = JSON.parse(globalThis.localStorage?.getItem(KEY) || 'null'); return v && typeof v === 'object' ? { ...blank(), ...v } : blank(); } catch { return blank(); }
}
function persist(c) { try { globalThis.localStorage?.setItem(KEY, JSON.stringify(c)); } catch { /* storage blocked */ } }

/** Put user materials into the shared tables (called at start-up, after every change, and inside each solver job). */
export function registerCustom(custom) {
  if (!custom || typeof custom !== 'object') return;
  for (const cls of Object.keys(TABLES)) {
    for (const name of Object.keys(TABLES[cls])) if (!isBuiltin(cls, name) && !(custom[cls] && name in custom[cls])) delete TABLES[cls][name];
    for (const [name, m] of Object.entries(custom[cls] || {})) if (!isBuiltin(cls, name) && m && typeof m === 'object') TABLES[cls][name] = m;
  }
}
/** Check a material before it is saved. Returns { errors: [], warnings: [] }. */
export function validate(cls, name, m, { editing = null } = {}) {
  const errors = [], warnings = [];
  if (!name || !String(name).trim()) errors.push('Give the material a name.');
  else if (isBuiltin(cls, name)) errors.push('That name belongs to a built-in material. Choose a different name.');
  else if (name !== editing && name in TABLES[cls]) errors.push('You already have a material with that name.');
  for (const f of FIELDS[cls]) {
    const v = m[f.key];
    if (v == null || v === '') { if (f.req) errors.push(`${f.label} is required.`); continue; }
    if (!Number.isFinite(v)) errors.push(`${f.label} must be a number.`);
    else if (v <= 0 && !['b', 'c', 'a1', 'a2', 'flash_C', 'co2_per_kg', 'lifecycle_factor'].includes(f.key)) errors.push(`${f.label} must be greater than zero.`);
  }
  if (errors.length) return { errors, warnings };
  if (cls === 'metals') {
    if (m.nu >= 0.5) errors.push('Poisson’s ratio must be below 0.5.');
    if (m.Sy > m.Su) errors.push('Yield strength cannot exceed ultimate strength.');
    const Giso = m.E / (2 * (1 + m.nu)); if (Math.abs(m.G / Giso - 1) > 0.15) warnings.push(`Shear modulus differs by ${Math.round(100 * (m.G / Giso - 1))}% from E/(2(1+ν)) = ${(Giso / 1e9).toFixed(1)} GPa. Fine for anisotropic products; check otherwise.`);
    if (m.b >= 0 || m.c >= 0) errors.push('Fatigue exponents b and c must be negative.');
    if (m.Sy_A && m.Sy_A > m.Sy) warnings.push('The design yield allowable is above the typical yield strength; allowables are normally lower.');
    if (m.Su_A && m.Su_A > m.Su) warnings.push('The design ultimate allowable is above the typical ultimate strength; allowables are normally lower.');
    if ((m.Sy_A > 0) !== (m.Su_A > 0)) warnings.push('Give both design allowables or neither; margins use them only as a pair.');
  }
  if (cls === 'plies') { if (m.E2 > m.E1) warnings.push('E2 is larger than E1: direction 1 should be the fibre direction.'); if (m.nu12 >= 1) errors.push('ν12 must be below 1.'); }
  if (cls === 'batteries' && !(m.v_min < m.v_nom && m.v_nom < m.v_max)) errors.push('Voltages must satisfy minimum < nominal < maximum.');
  return { errors, warnings };
}
export function saveCustom(cls, name, m, previousName = null) {
  const c = customMaterials();
  if (previousName && previousName !== name) delete c[cls][previousName];
  c[cls][name] = { ...m, custom: true };
  persist(c); registerCustom(c);
}
export function deleteCustom(cls, name) { const c = customMaterials(); delete c[cls][name]; persist(c); registerCustom(c); }

// ---- exchange ------------------------------------------------------------------------------------
/** CSV with one material per row, SI units, columns = class, name, then every property key of that class. */
export function toCsvRows(cls, names = Object.keys(TABLES[cls])) {
  const keys = FIELDS[cls].map((f) => f.key);
  return { columns: ['class', 'name', ...keys, 'note'], rows: names.map((n) => ['' + cls, n, ...keys.map((k) => TABLES[cls][n][k] ?? ''), TABLES[cls][n].allow_ref || TABLES[cls][n].note || '']) };
}
/** Parse materials from JSON ({metals:{name:{…}}}) or CSV as written by toCsvRows. Returns [{cls, name, data}]. */
export function parseMaterials(text) {
  const t = text.trim(), out = [];
  if (t.startsWith('{')) {
    const o = JSON.parse(t);
    for (const cls of Object.keys(TABLES)) for (const [name, data] of Object.entries(o[cls] || {})) if (data && typeof data === 'object') out.push({ cls, name, data: pick(cls, data) });
    return out;
  }
  const lines = t.split(/\r?\n/).filter((l) => l.trim() && !l.startsWith('#')), split = (l) => { const r = []; let cur = '', q = false; for (const ch of l) { if (ch === '"') q = !q; else if ((ch === ',' || ch === ';' || ch === '\t') && !q) { r.push(cur); cur = ''; } else cur += ch; } r.push(cur); return r.map((s) => s.trim()); };
  const head = split(lines[0]).map((s) => s.toLowerCase()), iC = head.indexOf('class'), iN = head.indexOf('name');
  if (iN < 0) throw new Error('The first row must name the columns and include “name”.');
  for (const l of lines.slice(1)) {
    const cells = split(l), cls = (iC >= 0 ? cells[iC] : 'metals').toLowerCase(); if (!TABLES[cls]) continue;
    const data = {}; for (const f of FIELDS[cls]) { const j = head.indexOf(f.key.toLowerCase()); if (j >= 0 && cells[j] !== '' && cells[j] != null) data[f.key] = Number(cells[j]); }
    const jn = head.indexOf('note'); if (jn >= 0 && cells[jn]) data.note = cells[jn].slice(0, 300);
    out.push({ cls, name: cells[iN], data });
  }
  return out;
}
function pick(cls, d) { const o = {}; for (const f of FIELDS[cls]) if (d[f.key] != null && d[f.key] !== '') o[f.key] = Number(d[f.key]); if (typeof d.allow_ref === 'string') o.allow_ref = d.allow_ref.slice(0, 300); if (typeof d.note === 'string') o.note = d.note.slice(0, 300); if (cls === 'metals' && d.JC && typeof d.JC === 'object') o.JC = Object.fromEntries(['A', 'B', 'n', 'C', 'm', 'Tm'].filter((k) => Number.isFinite(Number(d.JC[k]))).map((k) => [k, Number(d.JC[k])])); return o; }

/** Extend a suite's material drop-down with the user's materials (used when describing analyses). */
export function extendOptions(options) {
  if (!Array.isArray(options)) return options;
  for (const cls of Object.keys(TABLES)) { const b = builtinNames(cls); if (b.length && b.every((n) => options.includes(n))) return [...new Set([...options, ...Object.keys(TABLES[cls])])]; }
  return options;
}
