// Suite 21 — Materials and Composite Mechanics.
// Classical laminate theory with hygrothermal loads, ply failure criteria (maximum stress, Tsai–Hill, Tsai–Wu,
// Hashin, Puck inter-fibre), first-ply and ply-discount last-ply failure, micromechanics, interlaminar fracture
// and shear, viscoelastic creep / Norton relaxation, moisture diffusion and pitting, and metal yield criteria.

import * as N from '../core/numerics.js';
import { METALS, PLIES } from '../data/materials.js';

// ---- shared helpers -------------------------------------------------------------------------
const PLY_NAMES = Object.keys(PLIES), MATS = Object.keys(METALS);
const plyOf = (name) => PLIES[name] || PLIES['IM7/8552 carbon-epoxy'];
const kpi = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const DATA_NOTE = 'Ply and metal data are typical handbook values, not statistically based design allowables';
const rep = (a, n) => { const o = []; for (let k = 0; k < n; k++) o.push(...a); return o; };

/** Parse a stacking sequence such as "[0/45/-45/90]s", "[0_2/±45/90]s", "[0/90]2s" or "[(0/90)2/45]s" into ply angles [deg]. Returns null when it cannot be read. */
export function parseLayup(str) {
  const tokOf = (t) => {
    let m = /^\((.*)\)(\d*)$/.exec(t); if (m) return rep(seqOf(m[1]), +(m[2] || 1));
    m = /^(±|\+|-)?(\d+(?:\.\d+)?)(?:_(\d+))?$/.exec(t); if (!m) throw new Error('bad token');
    const a = +m[2]; return rep(m[1] === '±' ? [a, -a] : [m[1] === '-' ? -a : a], +(m[3] || 1));
  };
  const seqOf = (s) => { const out = []; let depth = 0, cur = ''; for (const ch of s + '/') { if (ch === '(') depth++; if (ch === ')') depth--; if (ch === '/' && depth === 0) { if (cur) out.push(...tokOf(cur)); cur = ''; } else cur += ch; } return out; };
  try {
    let s = String(str || '').replace(/\s+/g, '').replace(/[−–]/g, '-').replace(/\+\/-|\+-/g, '±'), suf = '';
    const m = /^\[(.*)\]([0-9a-zA-Z]*)$/.exec(s); if (m) { s = m[1]; suf = m[2].toLowerCase(); }
    let seq = seqOf(s);
    for (const tok of suf.match(/\d+|[a-z]/g) || []) { if (/\d/.test(tok)) seq = rep(seq, +tok); else if (tok === 's') seq = [...seq, ...seq.slice().reverse()]; else if (tok !== 't') return null; }
    return seq.length && seq.length <= 400 && seq.every(Number.isFinite) ? seq : null;
  } catch { return null; }
}

/** Reduced stiffnesses [Q11, Q12, Q22, Q66] of an orthotropic ply. */
const Qof = (p) => { const d = 1 - (p.nu12 ** 2 * p.E2) / p.E1; return [p.E1 / d, (p.nu12 * p.E2) / d, p.E2 / d, p.G12]; };
function qbar([Q11, Q12, Q22, Q66], thDeg) {
  const c = Math.cos(N.rad(thDeg)), s = Math.sin(N.rad(thDeg)), c2 = c * c, s2 = s * s;
  const b11 = Q11 * c2 * c2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * s2 * s2, b22 = Q11 * s2 * s2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * c2 * c2;
  const b12 = (Q11 + Q22 - 4 * Q66) * s2 * c2 + Q12 * (s2 * s2 + c2 * c2), b66 = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * s2 * c2 + Q66 * (s2 * s2 + c2 * c2);
  const b16 = (Q11 - Q12 - 2 * Q66) * s * c * c2 + (Q12 - Q22 + 2 * Q66) * s * s2 * c, b26 = (Q11 - Q12 - 2 * Q66) * s * s2 * c + (Q12 - Q22 + 2 * Q66) * s * c * c2;
  return [[b11, b12, b16], [b12, b22, b26], [b16, b26, b66]];
}
/** Free expansion strains in laminate axes (engineering shear) for principal coefficients e1, e2. */
const expBar = (e1, e2, thDeg) => { const c = Math.cos(N.rad(thDeg)), s = Math.sin(N.rad(thDeg)); return [e1 * c * c + e2 * s * s, e1 * s * s + e2 * c * c, 2 * (e1 - e2) * s * c]; };

/** Build the laminate: plies[k] holds the (possibly degraded) properties of ply k, bottom to top. */
function laminate(angles, plies, t) {
  const n = angles.length, h = n * t, A = N.zeros(3), B = N.zeros(3), D = N.zeros(3), NT = [0, 0, 0], MT = [0, 0, 0], NH = [0, 0, 0], MH = [0, 0, 0], lay = [];
  for (let k = 0; k < n; k++) {
    const p = plies[k], z0 = -h / 2 + k * t, z1 = z0 + t, Q = Qof(p), Qb = qbar(Q, angles[k]), al = expBar(p.a1, p.a2, angles[k]), be = expBar(p.b1 || 0, p.b2 || 0, angles[k]);
    const qa = N.matvec(Qb, al), qh = N.matvec(Qb, be);
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) { A[r][c] += Qb[r][c] * t; B[r][c] += (Qb[r][c] * (z1 * z1 - z0 * z0)) / 2; D[r][c] += (Qb[r][c] * (z1 ** 3 - z0 ** 3)) / 3; }
      NT[r] += qa[r] * t; MT[r] += (qa[r] * (z1 * z1 - z0 * z0)) / 2; NH[r] += qh[r] * t; MH[r] += (qh[r] * (z1 * z1 - z0 * z0)) / 2;
    }
    lay.push({ th: angles[k], z0, z1, Q, Qb, p });
  }
  const ABD = N.range(6, (r) => N.range(6, (c) => (r < 3 ? (c < 3 ? A[r][c] : B[r][c - 3]) : c < 3 ? B[r - 3][c] : D[r - 3][c - 3])));
  return { n, h, t, A, B, D, ABD, abd: N.inv(ABD), NT, MT, NH, MH, lay };
}
/** Mid-plane strains and curvatures for loads [Nx,Ny,Nxy,Mx,My,Mxy] plus uniform ΔT and moisture ΔM. */
const lamStrain = (L, load, dT = 0, dM = 0) => N.matvec(L.abd, load.map((v, k) => v + dT * (k < 3 ? L.NT[k] : L.MT[k - 3]) + dM * (k < 3 ? L.NH[k] : L.MH[k - 3])));
/** Stresses in material axes [σ1, σ2, τ12] and laminate axes at height z of ply k. */
function plyStress(L, k, z, e, dT = 0, dM = 0) {
  const { th, Q, p } = L.lay[k], c = Math.cos(N.rad(th)), s = Math.sin(N.rad(th)), ex = e[0] + z * e[3], ey = e[1] + z * e[4], gxy = e[2] + z * e[5];
  const e1 = ex * c * c + ey * s * s + gxy * s * c - p.a1 * dT - (p.b1 || 0) * dM, e2 = ex * s * s + ey * c * c - gxy * s * c - p.a2 * dT - (p.b2 || 0) * dM, g12 = 2 * (ey - ex) * s * c + gxy * (c * c - s * s);
  const s1 = Q[0] * e1 + Q[1] * e2, s2 = Q[1] * e1 + Q[2] * e2, t12 = Q[3] * g12;
  return { m: [s1, s2, t12], sx: s1 * c * c + s2 * s * s - 2 * t12 * s * c, sy: s1 * s * s + s2 * c * c + 2 * t12 * s * c, e1: e1 + p.a1 * dT + (p.b1 || 0) * dM, e2, g12 };
}

// ---- ply failure criteria: each returns { fi, mode } with failure at fi = 1 -------------------
const sgn = (s, t, c) => (s >= 0 ? t : c);
export const CRIT = {
  'Maximum stress': (s1, s2, t, p) => { const a = s1 >= 0 ? s1 / p.Xt : -s1 / p.Xc, b = s2 >= 0 ? s2 / p.Yt : -s2 / p.Yc, c = Math.abs(t) / p.S, fi = Math.max(a, b, c); return { fi, mode: fi === a ? (s1 >= 0 ? 'fibre tension' : 'fibre compression') : fi === b ? (s2 >= 0 ? 'matrix tension' : 'matrix compression') : 'in-plane shear' }; },
  'Tsai–Hill': (s1, s2, t, p) => { const X = sgn(s1, p.Xt, p.Xc), Y = sgn(s2, p.Yt, p.Yc); return { fi: Math.sqrt((s1 / X) ** 2 - (s1 * s2) / (X * X) + (s2 / Y) ** 2 + (t / p.S) ** 2), mode: 'interactive' }; },
  'Tsai–Wu': (s1, s2, t, p) => { const F11 = 1 / (p.Xt * p.Xc), F22 = 1 / (p.Yt * p.Yc); return { fi: (1 / p.Xt - 1 / p.Xc) * s1 + (1 / p.Yt - 1 / p.Yc) * s2 + F11 * s1 * s1 + F22 * s2 * s2 + (t / p.S) ** 2 - Math.sqrt(F11 * F22) * s1 * s2, mode: 'interactive' }; },
  Hashin: (s1, s2, t, p) => { const f = hashinF(s1, s2, t, p), m = hashinM(s1, s2, t, p); return f >= m ? { fi: f, mode: s1 >= 0 ? 'fibre tension' : 'fibre compression' } : { fi: m, mode: s2 >= 0 ? 'matrix tension' : 'matrix compression' }; },
  Puck: (s1, s2, t, p) => {
    const ff = s1 >= 0 ? s1 / p.Xt : -s1 / p.Xc, pp = p.pPlus ?? 0.35, pm = p.pMinus ?? 0.3; let iff, mode;
    if (s2 >= 0) { iff = Math.sqrt((t / p.S) ** 2 + ((1 - (pp * p.Yt) / p.S) * s2 / p.Yt) ** 2) + (pp * s2) / p.S; mode = 'inter-fibre mode A'; }
    else {
      const RA = (p.S / (2 * pm)) * (Math.sqrt(1 + (2 * pm * p.Yc) / p.S) - 1), ptt = (pm * RA) / p.S, tc = p.S * Math.sqrt(1 + 2 * ptt);
      if (Math.abs(s2) * tc <= RA * Math.abs(t)) { iff = (Math.sqrt(t * t + (pm * s2) ** 2) + pm * s2) / p.S; mode = 'inter-fibre mode B'; }
      else { iff = (((t / (2 * (1 + ptt) * p.S)) ** 2 + (s2 / p.Yc) ** 2) * p.Yc) / -s2; mode = 'inter-fibre mode C (wedge)'; }
    }
    return ff >= iff ? { fi: ff, mode: s1 >= 0 ? 'fibre tension' : 'fibre compression' } : { fi: iff, mode };
  },
};
const hashinF = (s1, s2, t, p) => (s1 >= 0 ? (s1 / p.Xt) ** 2 + (t / p.S) ** 2 : (s1 / p.Xc) ** 2);
function hashinM(s1, s2, t, p) { if (s2 >= 0) return (s2 / p.Yt) ** 2 + (t / p.S) ** 2; const ST = (p.stRatio ?? 0.378) * p.Yc; return (s2 / (2 * ST)) ** 2 + ((p.Yc / (2 * ST)) ** 2 - 1) * (s2 / p.Yc) + (t / p.S) ** 2; }
const CRIT_NAMES = Object.keys(CRIT);
/** Reserve factor on the mechanical stress sm with the residual stress sr held fixed: smallest λ with fi(sr + λ·sm) = 1. */
function reserve(fn, sm, sr, p) {
  const f = (l) => fn(sr[0] + l * sm[0], sr[1] + l * sm[1], sr[2] + l * sm[2], p) - 1;
  if (f(0) >= 0) return 0;
  let hi = 1, lo = 0; for (let k = 0; k < 80 && f(hi) < 0; k++) { lo = hi; hi *= 2; }
  if (f(hi) < 0) return Infinity;
  for (let k = 0; k < 48; k++) { const mid = 0.5 * (lo + hi); if (f(mid) < 0) lo = mid; else hi = mid; }
  return 0.5 * (lo + hi);
}
/** Governing reserve factor of a laminate for a criterion: scans the top and bottom of every ply. */
function lamReserve(L, load, dT, dM, name, skip) {
  const em = lamStrain(L, load), er = lamStrain(L, [0, 0, 0, 0, 0, 0], dT, dM), fn = CRIT[name]; let best = { rf: Infinity, ply: -1, mode: '-' }; const per = [];
  for (let k = 0; k < L.n; k++) {
    let rk = Infinity, mk = '-';
    if (!(skip && skip[k])) for (const z of [L.lay[k].z0, L.lay[k].z1]) {
      const sm = plyStress(L, k, z, em).m, sr = plyStress(L, k, z, er, dT, dM).m, p = L.lay[k].p, rf = reserve((a, b, c, q) => fn(a, b, c, q).fi, sm, sr, p);
      if (rf < rk) { rk = rf; mk = Number.isFinite(rf) ? fn(sr[0] + rf * sm[0], sr[1] + rf * sm[1], sr[2] + rf * sm[2], p).mode : '-'; }
    }
    per.push(rk); if (rk < best.rf) best = { rf: rk, ply: k, mode: mk };
  }
  return { ...best, per };
}
/** Laminate from analysis inputs; falls back to a quasi-isotropic stack when the text cannot be parsed. */
function lamFromInputs(i, warnings) {
  let ang = parseLayup(i.layup);
  if (!ang) { ang = [0, 45, -45, 90, 90, -45, 45, 0]; warnings.push(`Could not read the stacking sequence "${i.layup}"; a [0/45/-45/90]s laminate was used instead.`); }
  const base = plyOf(i.ply), p = { ...base, b1: 0, b2: i.beta2 ?? 0, pPlus: i.pPlus, pMinus: i.pMinus, stRatio: i.stRatio }, t = (i.t_ply_mm > 0 ? i.t_ply_mm / 1e3 : base.t) * Math.max(1, Math.round(i.n_block || 1));
  return { ang, p, t, L: laminate(ang, ang.map(() => p), t) };
}
const engConst = (L) => { const a = L.abd, h = L.h; return { Ex: 1 / (h * a[0][0]), Ey: 1 / (h * a[1][1]), Gxy: 1 / (h * a[2][2]), nuxy: -a[0][1] / a[0][0], Exf: 12 / (h ** 3 * a[3][3]), Eyf: 12 / (h ** 3 * a[4][4]) }; };
const LAM_INPUTS = [
  { key: 'layup', label: 'Stacking sequence', type: 'text', default: '[0/45/-45/90]s', group: 'Laminate', help: 'Angles in degrees from the load axis, bottom to top: [0/45/-45/90]s, [0_2/±45/90]s, [0/90]2s, [(0/90)2/45]s' },
  { key: 'ply', label: 'Ply material', type: 'select', options: PLY_NAMES, default: 'IM7/8552 carbon-epoxy', group: 'Laminate', help: 'Typical unidirectional ply data, not design allowables' },
  { key: 't_ply_mm', label: 'Cured ply thickness (0 = database)', unit: 'mm', default: 0, min: 0, max: 2, group: 'Laminate' },
  { key: 'n_block', label: 'Plies per listed angle', unit: '', default: 1, min: 1, max: 200, step: 1, discrete: true, group: 'Laminate', help: 'Scales the whole laminate thickness. In-plane response is identical whether the extra plies are blocked or dispersed; disperse them in the real laminate' },
];
/** Running load [N/m] in a wing cover at limit load, estimated from the case when Suite 2 has not run. */
function coverLoad(c, up, d) {
  if (up.fea?.cover_load_Npm > 0 && c.wing.S_m2 > 0) return up.fea.cover_load_Npm;
  if (!(c.wing.S_m2 > 0)) return undefined;
  const n = up.performance?.n_limit ?? c.aero.n_pos, M = 0.8 * n * (d.W / 2) * 0.42 * (c.wing.b_m / 2), w = c.struct.box_chord_frac * d.c_root, h = c.struct.box_height_frac * c.wing.tc * d.c_root;
  return M / (w * h);
}
/** Plies per angle needed for the first-ply-failure reserve to reach the ultimate factor under Nx (defaults only). */
function blockFor(layup, ply, Nx, sf) {
  const w = [], { L } = lamFromInputs({ layup, ply, t_ply_mm: 0, n_block: 1 }, w), rf = lamReserve(L, [Nx, 0, 0, 0, 0, 0], 0, 0, 'Tsai–Wu').rf;
  return Number.isFinite(rf) && rf > 0 ? N.clamp(Math.ceil((1.05 * sf) / rf), 1, 200) : 1;
}

// ---- 1. classical laminate theory -----------------------------------------------------------
const clt = {
  id: 'clt', title: 'Laminate stiffness and ply stresses (classical laminate theory)', fidelity: 'analytical',
  summary: 'ABD stiffness matrices, engineering constants and directional stiffness of the laminate, then ply-by-ply stresses and reserve factors under in-plane loads, moments, a temperature change and absorbed moisture.',
  equations: ['Generalized Hooke’s law', 'Orthotropic constitutive equations', 'Classical laminate theory equations', 'Tsai–Hill failure criterion', 'Tsai–Wu failure criterion', 'Hashin failure criteria', 'Puck failure criteria'],
  inputs: [...LAM_INPUTS,
    { key: 'Nx', label: 'Running load Nx', unit: 'N/m', default: 2e5, group: 'Loads', help: 'Force per unit width along the 0° direction (limit load)' },
    { key: 'Ny', label: 'Running load Ny', unit: 'N/m', default: 0, group: 'Loads' },
    { key: 'Nxy', label: 'Shear flow Nxy', unit: 'N/m', default: 0, group: 'Loads' },
    { key: 'Mx', label: 'Bending moment Mx', unit: 'N·m/m', default: 0, group: 'Loads' },
    { key: 'My', label: 'Bending moment My', unit: 'N·m/m', default: 0, group: 'Loads' },
    { key: 'Mxy', label: 'Twisting moment Mxy', unit: 'N·m/m', default: 0, group: 'Loads' },
    { key: 'dT', label: 'Temperature change from stress-free state', unit: 'K', default: 0, min: -300, max: 300, group: 'Environment', help: 'About −100 to −160 K from cure to room temperature for 180 °C-cure epoxy' },
    { key: 'dM', label: 'Absorbed moisture (mass fraction)', unit: '-', default: 0, min: 0, max: 0.05, group: 'Environment', help: 'e.g. 0.01 = 1% by weight; see the moisture analysis' },
    { key: 'beta2', label: 'Transverse swelling coefficient β₂', unit: '1/(mass fraction)', default: 0.4, min: 0, max: 1.5, group: 'Environment', help: 'Strain per unit moisture mass fraction; roughly 0.3–0.6 for carbon/epoxy. Measure for your system' },
    { key: 'criterion', label: 'Governing failure criterion', type: 'select', options: CRIT_NAMES, default: 'Tsai–Wu', group: 'Strength' },
    { key: 'sf', label: 'Required reserve factor at limit load', unit: '-', default: 1.5, min: 1, max: 4, group: 'Strength' },
    { key: 'pPlus', label: 'Puck slope p⊥∥(+)', unit: '-', default: 0.35, min: 0.1, max: 0.5, group: 'Strength', help: '0.35 carbon, 0.30 glass (Puck recommendations)' },
    { key: 'pMinus', label: 'Puck slope p⊥∥(−)', unit: '-', default: 0.3, min: 0.1, max: 0.5, group: 'Strength', help: '0.30 carbon, 0.25 glass' },
    { key: 'stRatio', label: 'Hashin transverse shear strength / Yc', unit: '-', default: 0.378, min: 0.2, max: 0.8, group: 'Strength', help: '0.378 corresponds to a 53° fracture plane; 0.5 removes the linear compression term' },
  ],
  defaults(c, up, d) { const Nx = coverLoad(c, up, d); return { layup: c.struct.layup, ply: c.struct.ply, sf: c.struct.sf_ultimate, Nx, n_block: blockFor(c.struct.layup, c.struct.ply, Nx ?? 2e5, c.struct.sf_ultimate) }; },
  run(i) {
    const warnings = [], { ang, p, L } = lamFromInputs(i, warnings), E = engConst(L), load = [i.Nx, i.Ny, i.Nxy, i.Mx, i.My, i.Mxy];
    const e = lamStrain(L, load, i.dT, i.dM), res = Object.fromEntries(CRIT_NAMES.map((nm) => [nm, lamReserve(L, load, i.dT, i.dM, nm)])), gov = res[i.criterion] || res['Tsai–Wu'];
    const scaleB = Math.max(...L.B.flat().map(Math.abs)) / (Math.max(...L.A.flat().map(Math.abs)) * L.h), sym = scaleB < 1e-9, bal = Math.abs(L.A[0][2]) + Math.abs(L.A[1][2]) < 1e-9 * L.A[0][0];
    const frac = (a) => ang.filter((v) => Math.abs(((Math.abs(v) % 180) + 180) % 180 - a) < 1 || (a === 0 && Math.abs(Math.abs(v) % 180 - 180) < 1)).length / ang.length, f0 = frac(0), f90 = frac(90), f45 = ang.filter((v) => Math.abs(Math.abs(v) - 45) < 1).length / ang.length;
    if (!sym) warnings.push('The laminate is not symmetric: extension–bending coupling (B ≠ 0) will warp the part on cool-down from cure.');
    if (!bal) warnings.push('The laminate is not balanced: extension–shear coupling (A16, A26 ≠ 0).');
    if (Math.min(f0, f45, f90) < 0.0999 && ang.length >= 8) warnings.push('Fewer than 10% of plies in one of the 0°, ±45° or 90° directions: common design practice keeps at least 10% in each to cover unexpected load paths.');
    if (gov.rf < 1) warnings.push(`Ply ${gov.ply + 1} (${ang[gov.ply]}°) is predicted to fail (${gov.mode}) below the applied load.`);
    // through-thickness profiles
    const zs = [], sxs = [], s1s = [], s2s = [];
    for (let k = 0; k < L.n; k++) for (const z of [L.lay[k].z0, L.lay[k].z1]) { const s = plyStress(L, k, z, e, i.dT, i.dM); zs.push(z * 1e3); sxs.push(s.sx / 1e6); s1s.push(s.m[0] / 1e6); s2s.push(s.m[1] / 1e6); }
    const a = N.inv(L.A).map((r) => r.map((v) => v * L.h)), th = N.linspace(0, 360, 73), Eth = th.map((d) => { const c = Math.cos(N.rad(d)), s = Math.sin(N.rad(d)); return 1 / (a[0][0] * c ** 4 + (2 * a[0][1] + a[2][2]) * s * s * c * c + a[1][1] * s ** 4 + 2 * a[0][2] * c ** 3 * s + 2 * a[1][2] * c * s ** 3) / 1e9; });
    const Gth = th.map((d) => { const c = Math.cos(N.rad(d)), s = Math.sin(N.rad(d)); return 1 / (4 * (a[0][0] + a[1][1] - 2 * a[0][1]) * s * s * c * c + a[2][2] * (c * c - s * s) ** 2 + 4 * (a[0][2] - a[1][2]) * s * c * (c * c - s * s)) / 1e9; });
    const cte = N.matvec(L.abd, [...L.NT, ...L.MT]), nShow = Math.min(L.n, 40), rfMin = Math.min(...CRIT_NAMES.map((nm) => res[nm].rf)), fmt = (v) => (Number.isFinite(v) ? v : 1e9);
    return {
      kpis: [
        kpi('Ex_Pa', 'In-plane modulus Ex', E.Ex, 'Pa'), kpi('Ey_Pa', 'In-plane modulus Ey', E.Ey, 'Pa'), kpi('Gxy_Pa', 'In-plane shear modulus Gxy', E.Gxy, 'Pa'), kpi('nu_xy', 'Poisson ratio νxy', E.nuxy, '-'),
        kpi('Ex_flex_Pa', 'Flexural modulus (x)', E.Exf, 'Pa'), kpi('laminate_t_m', 'Laminate thickness', L.h, 'm'), kpi('n_plies', 'Number of plies', L.n * Math.max(1, Math.round(i.n_block)), '-'),
        kpi('areal_mass_kgm2', 'Areal mass', p.rho * L.h, 'kg/m²'),
        kpi('min_RF', `Reserve factor (${i.criterion})`, gov.rf, '-', gov.rf >= i.sf ? 'ok' : gov.rf >= 1 ? 'warn' : 'bad', `First-ply failure load / applied load; ${i.sf} required. Critical ply ${gov.ply + 1} (${ang[gov.ply] ?? '-'}°), ${gov.mode}`),
        kpi('RF_lowest_any', 'Lowest reserve factor over all criteria', rfMin, '-'),
        kpi('eps_x', 'Mid-plane strain εx', e[0], '-', Math.abs(e[0]) < 0.006 ? 'ok' : 'warn', 'Design strains for damage-tolerant carbon laminates are typically limited to about 0.4–0.6%'),
        kpi('eps_y', 'Mid-plane strain εy', e[1], '-'), kpi('gamma_xy', 'Mid-plane shear strain γxy', e[2], '-'), kpi('kappa_x', 'Curvature κx', e[3], '1/m'),
        kpi('cte_x', 'Laminate thermal expansion αx', cte[0], '1/K'), kpi('cte_y', 'Laminate thermal expansion αy', cte[1], '1/K'),
        kpi('frac_0', 'Share of 0° plies', f0, '-'), kpi('frac_45', 'Share of ±45° plies', f45, '-'), kpi('frac_90', 'Share of 90° plies', f90, '-'),
      ],
      plots: [
        { type: 'polar', title: 'Directional in-plane stiffness', rlabel: 'Modulus [GPa]', series: [{ name: 'E(θ)', theta_deg: th, r: Eth }, { name: 'G(θ)', theta_deg: th, r: Gth }] },
        { type: 'line', title: 'Through-thickness stress', xlabel: 'Stress [MPa]', ylabel: 'z [mm]', series: [{ name: 'σx (laminate axis)', x: sxs, y: zs }, { name: 'σ1 (fibre direction)', x: s1s, y: zs, style: 'dash' }, { name: 'σ2 (transverse)', x: s2s, y: zs, style: 'dash' }] },
        { type: 'bar', title: 'Reserve factor by ply' + (L.n > nShow ? ` (first ${nShow} plies)` : ''), ylabel: 'Reserve factor [-] (capped at 10)', categories: N.range(nShow, (k) => `${k + 1}: ${ang[k]}°`), series: ['Maximum stress', 'Tsai–Wu', 'Hashin', 'Puck'].map((nm) => ({ name: nm, y: res[nm].per.slice(0, nShow).map((v) => Math.min(10, v)) })) },
      ],
      tables: [
        { title: 'ABD matrix (A [N/m], B [N], D [N·m])', columns: ['', '1', '2', '6', '1 ', '2 ', '6 '], rows: L.ABD.map((r, k) => [['A/B 1', 'A/B 2', 'A/B 6', 'B/D 1', 'B/D 2', 'B/D 6'][k], ...r]) },
        { title: 'Reserve factors by criterion', columns: ['Criterion', 'Reserve factor', 'Critical ply', 'Angle [deg]', 'Mode'], rows: CRIT_NAMES.map((nm) => [nm, fmt(res[nm].rf), res[nm].ply + 1, ang[res[nm].ply] ?? 0, res[nm].mode]) },
      ],
      outputs: { nu_xy: E.nuxy, symmetric: sym ? 1 : 0, balanced: bal ? 1 : 0 },
      warnings,
      models: ['Classical laminate theory (Kirchhoff)', 'Hygrothermal equivalent loads', `Ply failure: ${CRIT_NAMES.join(', ')}`],
      assumptions: ['Thin laminate, perfectly bonded plies, plane stress in each ply, linear elastic to first failure', 'Reserve factors scale the mechanical load with thermal and moisture residual stresses held constant', 'Free-edge interlaminar stresses, holes, impact damage and in-situ strength effects are not included', DATA_NOTE],
    };
  },
  verify() {
    const p = PLIES['T300/5208 carbon-epoxy'], Q = Qof(p), qi = N.kv(clt.run(base({ layup: '[0/45/-45/90]s', ply: 'T300/5208 carbon-epoxy' })));
    const U1 = (3 * Q[0] + 3 * Q[2] + 2 * Q[1] + 4 * Q[3]) / 8, U4 = (Q[0] + Q[2] + 6 * Q[1] - 4 * Q[3]) / 8, U5 = (Q[0] + Q[2] - 2 * Q[1] + 4 * Q[3]) / 8;
    const ud = N.kv(clt.run(base({ layup: '[0]4', ply: 'T300/5208 carbon-epoxy', Nx: 1e5 }))), th = N.kv(clt.run(base({ layup: '[0]8', ply: 'T300/5208 carbon-epoxy', Nx: 0, dT: -150 })));
    const cp = lamFromInputs(base({ layup: '[0/90]s', ply: 'T300/5208 carbon-epoxy' }), []), er = lamStrain(cp.L, [0, 0, 0, 0, 0, 0], -150), fsum = N.sum(N.range(4, (k) => plyStress(cp.L, k, 0, er, -150).sx * cp.L.t));
    const pk = { ...p, pPlus: 0.35, pMinus: 0.3 }, RA = (p.S / 0.6) * (Math.sqrt(1 + (0.6 * p.Yc) / p.S) - 1), tc = p.S * Math.sqrt(1 + (0.6 * RA) / p.S);
    return [
      N.check('Quasi-isotropic Ex = (U1² − U4²)/U1', qi.Ex_Pa, (U1 * U1 - U4 * U4) / U1, 1e-10, 'Tsai & Pagano laminate invariants'),
      N.check('Quasi-isotropic Ex = Ey', qi.Ey_Pa, qi.Ex_Pa, 1e-10, 'In-plane isotropy of π/4 laminates'),
      N.check('Quasi-isotropic Gxy = U5 = E/(2(1+ν))', qi.Gxy_Pa, U5, 1e-10, 'Laminate invariants'),
      N.check('T300/5208 quasi-isotropic modulus ≈ 69.7 GPa', qi.Ex_Pa, 69.7e9, 5e-3, 'Tsai & Hahn, Introduction to Composite Materials (for the tabulated ply data)'),
      N.check('Unidirectional laminate returns E1', ud.Ex_Pa, p.E1, 1e-10, 'Limit case'),
      N.check('Unidirectional 0° tension: RF = Xt·h/Nx (Tsai–Wu)', ud.min_RF, (p.Xt * 4 * p.t) / 1e5, 1e-6, 'Uniaxial strength'),
      N.check('Free thermal expansion of a unidirectional laminate is stress-free (αx = α1)', th.cte_x, p.a1, 1e-9, 'Limit case'),
      N.check('Cross-ply thermal residual stresses self-equilibrate', fsum / (cp.L.h * 1e6), 0, 1e-9, 'Force equilibrium with no applied load'),
      N.check('Parser: [0_2/±45/90]s has 10 plies', parseLayup('[0_2/±45/90]s').length, 10, 0, 'Stacking-sequence notation'),
      N.check('Parser: [0/90]2s mirrors after repeating', parseLayup('[0/90]2s').join(','), '0,90,0,90,90,0,90,0', 0, 'Stacking-sequence notation'),
      N.check('Puck: pure transverse tension fails at Yt', CRIT.Puck(0, p.Yt, 0, pk).fi, 1, 1e-12, 'Puck & Schürmann (1998), mode A'),
      N.check('Puck: modes B and C meet at the transition point', CRIT.Puck(0, -RA * 1.0000001, tc, pk).fi, CRIT.Puck(0, -RA * 0.9999999, tc, pk).fi, 1e-6, 'Puck & Schürmann (1998)'),
      N.check('Puck: pure transverse compression fails at Yc', CRIT.Puck(0, -p.Yc, 0, pk).fi, 1, 1e-12, 'Puck mode C'),
      N.check('Hashin: pure transverse compression fails at Yc', hashinM(0, -p.Yc, 0, { ...p, stRatio: 0.378 }), 1, 1e-12, 'Hashin (1980)'),
      N.check('Tsai–Wu: uniaxial fibre compression fails at Xc', CRIT['Tsai–Wu'](-p.Xc, 0, 0, p).fi, 1, 1e-12, 'Tsai & Wu (1971)'),
    ].map((c) => (typeof c.actual === 'string' ? { ...c, actual: c.actual === c.expected ? 1 : 0, expected: 1, error: c.actual === c.expected ? 0 : 1, pass: c.actual === c.expected } : c));
  },
  validation: [{ name: 'Quasi-isotropic T300/5208 laminate modulus', source: 'Laminate invariants for the tabulated T300/5208 ply (Tsai & Hahn): E = 69.7 GPa independent of ply block count', inputs: { layup: '[0/45/-45/90]s', ply: 'T300/5208 carbon-epoxy', t_ply_mm: 0, Nx: 1e5, Ny: 0, Nxy: 0, Mx: 0, My: 0, Mxy: 0, dT: 0, dM: 0 }, sweep: { key: 'n_block', values: [1, 2, 4] }, target: 'Ex_Pa', observed: [69.7e9, 69.7e9, 69.7e9], tol_pct: 1 }],
  calibration: { params: [{ key: 'beta2', min: 0, max: 1.5 }], sweep: 'dM', target: 'eps_y', note: 'Measured transverse swelling strain of a unidirectional coupon versus moisture content calibrates β₂.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.min_RF < 1) out.push({ severity: 'critical', title: 'First-ply failure below the applied load', detail: `Reserve factor ${o.min_RF.toFixed(2)} (${i.criterion}).`, action: `Increase thickness (about ${Math.ceil((i.sf / Math.max(o.min_RF, 1e-6)) * Math.max(1, i.n_block))} plies per listed angle would reach ${i.sf}), or move plies towards the principal load direction.`, basis: 'First-ply failure at limit load' });
    else if (o.min_RF < i.sf) out.push({ severity: 'warn', title: 'Reserve factor below the required value', detail: `${o.min_RF.toFixed(2)} against ${i.sf}.`, action: 'Add plies in the governing direction or reduce the running load; confirm with the progressive-failure analysis whether the first failure is benign matrix cracking.', basis: 'Ultimate factor on first-ply failure' });
    else if (o.min_RF > 2.5 * i.sf && Number.isFinite(o.min_RF)) out.push({ severity: 'advise', title: 'Laminate is far from first-ply failure', detail: `Reserve factor ${o.min_RF.toFixed(1)}.`, action: 'Strength does not size this laminate. Check stiffness, buckling, bearing and damage-tolerance strain limits before removing plies; each ply removed saves ' + `${(res.kpis.find((k) => k.key === 'areal_mass_kgm2').value / Math.max(1, res.kpis.find((k) => k.key === 'n_plies').value)).toFixed(3)} kg/m².`, basis: 'Reserve factor' });
    if (!o.symmetric) out.push({ severity: 'warn', title: 'Unsymmetric laminate', detail: 'B matrix is non-zero.', action: 'Mirror the stack about the mid-plane to avoid cure warpage and bending under in-plane load.', basis: 'Laminate design rule' });
    if (Math.abs(o.eps_x) > 0.005) out.push({ severity: 'advise', title: 'High design strain', detail: `εx = ${(100 * o.eps_x).toFixed(2)}%.`, action: 'Notched and impact-damaged carbon laminates usually limit design strain to roughly 0.4–0.5%; treat the unnotched reserve factor as optimistic.', basis: 'Damage-tolerance strain cut-off (programme-specific)' });
    return out;
  },
};
const base = (o) => ({ layup: '[0/45/-45/90]s', ply: 'IM7/8552 carbon-epoxy', t_ply_mm: 0, n_block: 1, Nx: 2e5, Ny: 0, Nxy: 0, Mx: 0, My: 0, Mxy: 0, dT: 0, dM: 0, beta2: 0.4, criterion: 'Tsai–Wu', sf: 1.5, pPlus: 0.35, pMinus: 0.3, stRatio: 0.378, ...o });

// ---- 2. first-ply / last-ply failure and envelopes ------------------------------------------
/** Ply-discount progressive failure under proportional loading. Returns the event list and the maximum load multiplier. */
function progressive(ang, p0, t, load, dT, kd) {
  const n = ang.length, st = ang.map(() => ({ m: false, f: false })), events = [];
  let peak = 0, lamPrev = 0;
  for (let guard = 0; guard < 3 * n + 2; guard++) {
    const plies = st.map((s) => (s.f ? { ...p0, E1: p0.E1 * 1e-3, E2: p0.E2 * 1e-3, G12: p0.G12 * 1e-3, nu12: p0.nu12 * 1e-3 } : s.m ? { ...p0, E2: p0.E2 * kd, G12: p0.G12 * kd, nu12: p0.nu12 * kd } : p0));
    let L; try { L = laminate(ang, plies, t); } catch { break; }
    const em = lamStrain(L, load), er = lamStrain(L, [0, 0, 0, 0, 0, 0], dT, 0); let best = { rf: Infinity };
    for (let k = 0; k < n; k++) if (!st[k].f) for (const z of [L.lay[k].z0, L.lay[k].z1]) {
      const sm = plyStress(L, k, z, em).m, sr = plyStress(L, k, z, er, dT).m, rf = reserve(hashinF, sm, sr, p0), rm = st[k].m ? Infinity : reserve(hashinM, sm, sr, p0);
      if (rf < best.rf) best = { rf, k, fibre: true }; if (rm < best.rf) best = { rf: rm, k, fibre: false };
    }
    if (!Number.isFinite(best.rf)) break;
    const lam = Math.max(best.rf, lamPrev), eBefore = lamStrain(L, load.map((v) => v * lam), dT, 0);
    events.push({ lam, raw: best.rf, ply: best.k, angle: ang[best.k], mode: best.fibre ? 'fibre' : 'matrix', ex: eBefore[0], ey: eBefore[1], gxy: eBefore[2], cascade: best.rf < lamPrev });
    if (best.fibre) st[best.k].f = true; else st[best.k].m = true;
    peak = Math.max(peak, best.rf); lamPrev = lam;
    if (best.fibre && best.rf < 0.5 * peak) break; // load has collapsed: the laminate cannot recover
    if (st.every((s) => s.f)) break;
  }
  return { events, lpf: peak, fpf: events.length ? events[0].raw : Infinity };
}
const strength = {
  id: 'strength', title: 'First-ply failure, progressive failure and failure envelope', fidelity: 'numerical',
  summary: 'Load at which the first ply fails, the ply-discount sequence of matrix and fibre failures up to laminate collapse, and the biaxial Nx–Ny failure envelope for several criteria.',
  equations: ['Classical laminate theory equations', 'Tsai–Wu failure criterion', 'Hashin failure criteria', 'Puck failure criteria', 'Continuum damage evolution equations'],
  inputs: [...LAM_INPUTS,
    { key: 'Nx', label: 'Reference load Nx', unit: 'N/m', default: 2e5, group: 'Loads', help: 'Limit running load; the analysis scales the whole reference load proportionally' },
    { key: 'Ny', label: 'Reference load Ny', unit: 'N/m', default: 0, group: 'Loads' },
    { key: 'Nxy', label: 'Reference shear flow Nxy', unit: 'N/m', default: 0, group: 'Loads' },
    { key: 'dT', label: 'Cure cool-down ΔT', unit: 'K', default: 0, min: -300, max: 100, group: 'Environment', help: 'Adds thermal residual stresses, which usually bring matrix cracking forward' },
    { key: 'criterion', label: 'First-ply failure criterion', type: 'select', options: CRIT_NAMES, default: 'Tsai–Wu', group: 'Strength' },
    { key: 'kd', label: 'Matrix-failed stiffness retention', unit: '-', default: 0.1, min: 0.001, max: 1, group: 'Strength', help: 'Factor on E2, G12 and ν12 of a ply after matrix failure (ply-discount model; calibrate against notched/unnotched coupon curves)' },
    { key: 'sf', label: 'Required reserve factor at limit load', unit: '-', default: 1.5, min: 1, max: 4, group: 'Strength' },
    { key: 'stRatio', label: 'Hashin transverse shear strength / Yc', unit: '-', default: 0.378, min: 0.2, max: 0.8, group: 'Strength' },
    { key: 'nDir', label: 'Envelope directions', unit: '', default: 72, min: 16, max: 360, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c, up, d) { const Nx = coverLoad(c, up, d); return { layup: c.struct.layup, ply: c.struct.ply, sf: c.struct.sf_ultimate, Nx, n_block: blockFor(c.struct.layup, c.struct.ply, Nx ?? 2e5, c.struct.sf_ultimate) }; },
  run(i) {
    const warnings = [], { ang, p, t, L } = lamFromInputs(i, warnings), load = [i.Nx, i.Ny, i.Nxy, 0, 0, 0], mag = Math.hypot(i.Nx, i.Ny, i.Nxy);
    const fp = lamReserve(L, load, i.dT, 0, i.criterion), pr = mag > 0 ? progressive(ang, p, t, load, i.dT, i.kd) : { events: [], lpf: Infinity, fpf: Infinity };
    const fx = lamReserve(L, [1, 0, 0, 0, 0, 0], i.dT, 0, i.criterion), cx = lamReserve(L, [-1, 0, 0, 0, 0, 0], i.dT, 0, i.criterion), ux = progressive(ang, p, t, [1, 0, 0, 0, 0, 0], i.dT, i.kd);
    if (!(mag > 0)) warnings.push('No reference load given: reserve factors are infinite; uniaxial strengths are still reported.');
    if (fp.rf === 0) warnings.push('Thermal residual stresses alone exceed the ply strength in this criterion: matrix cracking is predicted before any load is applied.');
    // envelope in the Nx–Ny plane
    const nD = Math.max(16, Math.round(i.nDir)), env = (nm, lpf) => { const x = [], y = [], nn = lpf ? Math.min(nD, 36) : nD; for (let k = 0; k <= nn; k++) { const a = (2 * Math.PI * k) / nn, dir = [Math.cos(a), Math.sin(a), 0, 0, 0, 0], r = lpf ? progressive(ang, p, t, dir, i.dT, i.kd).lpf : lamReserve(L, dir, i.dT, 0, nm).rf, v = Number.isFinite(r) ? r : NaN; x.push((v * dir[0]) / 1e6); y.push((v * dir[1]) / 1e6); } return { x, y }; };
    const series = [{ name: 'Tsai–Wu (first ply)', ...env('Tsai–Wu') }, { name: 'Maximum stress (first ply)', ...env('Maximum stress') }, { name: 'Puck (first ply)', ...env('Puck') }, { name: 'Hashin ply-discount (last ply)', ...env('Hashin', true), style: 'dash' }];
    if (mag > 0) series.push({ name: 'Applied (limit)', x: [i.Nx / 1e6], y: [i.Ny / 1e6], style: 'points' });
    const ev = pr.events, firstF = ev.find((e) => e.mode === 'fibre'), benign = fp.mode.includes('matrix') || fp.mode.includes('inter-fibre') || fp.mode.includes('shear');
    return {
      kpis: [
        kpi('min_RF', `First-ply failure reserve factor (${i.criterion})`, fp.rf, '-', fp.rf >= i.sf ? 'ok' : fp.rf >= 1 ? 'warn' : 'bad', `Ply ${fp.ply + 1} (${ang[fp.ply] ?? '-'}°), ${fp.mode}; ${i.sf} required`),
        kpi('RF_last_ply', 'Last-ply failure reserve factor (ply discount)', pr.lpf, '-', pr.lpf >= i.sf ? 'ok' : pr.lpf >= 1 ? 'warn' : 'bad'),
        kpi('lpf_over_fpf', 'Last-ply / first-ply load', pr.lpf / (pr.fpf || NaN), '-', undefined, 'Above 1: the laminate tolerates matrix cracking before collapse'),
        kpi('fpf_load_Npm', 'First-ply failure load, uniaxial tension Nx', fx.rf, 'N/m'),
        kpi('fpf_comp_Npm', 'First-ply failure load, uniaxial compression Nx', cx.rf, 'N/m', undefined, 'Material strength only: panel buckling is checked in Suite 2'),
        kpi('ult_load_Npm', 'Ultimate tensile load Nx (ply discount)', ux.lpf, 'N/m'),
        kpi('ult_strength_Pa', 'Unnotched tensile strength', ux.lpf / L.h, 'Pa'),
        kpi('fpf_strain', 'Strain at first-ply failure (uniaxial tension)', fx.rf * L.abd[0][0], '-'),
        kpi('laminate_t_m', 'Laminate thickness', L.h, 'm'),
        kpi('n_events', 'Ply failure events to collapse', ev.length, '-'),
      ],
      plots: [
        { type: 'line', title: 'Failure envelope in the Nx–Ny plane', xlabel: 'Nx [MN/m]', ylabel: 'Ny [MN/m]', equalAspect: true, series },
        { type: 'line', title: 'Progressive failure under the reference load path', xlabel: 'Mid-plane strain εx [%]', ylabel: 'Load multiplier [-]', series: [{ name: 'Load path (ply discount)', x: [0, ...ev.map((e) => 100 * e.ex)], y: [0, ...ev.map((e) => e.lam)], style: 'line+points' }], annotations: [{ y: 1, label: 'Limit' }, { y: i.sf, label: 'Ultimate' }] },
      ],
      tables: [{ title: 'Failure sequence under the reference load', columns: ['Event', 'Load multiplier', 'Ply', 'Angle [deg]', 'Mode', 'εx [%]', 'At constant load'], rows: ev.map((e, k) => [k + 1, e.lam, e.ply + 1, e.angle, e.mode, 100 * e.ex, e.cascade ? 'yes' : 'no']) }],
      outputs: { first_fibre_RF: firstF ? firstF.lam : NaN, fpf_benign: benign ? 1 : 0 },
      warnings,
      models: [`First-ply failure by ${i.criterion}`, 'Ply-discount progressive failure with Hashin fibre/matrix modes', 'Proportional loading with constant thermal residual stress'],
      assumptions: ['Sudden stiffness discount at ply failure (no fracture-energy regularisation, no delamination)', 'Load-controlled: plies that fail below the current load fail immediately', 'Unnotched laminate; compression strength excludes buckling and kink-band in-situ effects', DATA_NOTE],
    };
  },
  convergence: { param: 'nDir', label: 'Envelope directions', levels: [18, 36, 72, 144], metric: 'min_RF' },
  verify() {
    const p = PLIES['AS4/3501-6 carbon-epoxy'], b = { layup: '[0]8', ply: 'AS4/3501-6 carbon-epoxy', t_ply_mm: 0, n_block: 1, Nx: 1e5, Ny: 0, Nxy: 0, dT: 0, criterion: 'Maximum stress', kd: 0.1, sf: 1.5, stRatio: 0.378, nDir: 16 };
    const ud = N.kv(strength.run(b)), tr = N.kv(strength.run({ ...b, layup: '[90]8' })), cpI = { ...b, layup: '[0/90]s', criterion: 'Maximum stress' }, cp = N.kv(strength.run(cpI));
    // cross-ply ultimate by netting analysis after full matrix discount: only the 0° plies carry load, with degraded 90° stiffness
    const Lc = lamFromInputs(cpI, []).L, e0 = p.Xt / p.E1, kd = 0.1, Qd = Qof({ ...p, E2: p.E2 * kd, G12: p.G12 * kd, nu12: p.nu12 * kd });
    return [
      N.check('Unidirectional tension strength Xt·h', ud.fpf_load_Npm, p.Xt * 8 * p.t, 1e-6, 'Uniaxial strength'),
      N.check('Unidirectional compression strength Xc·h', ud.fpf_comp_Npm, p.Xc * 8 * p.t, 1e-6, 'Uniaxial strength'),
      N.check('Transverse laminate fails at Yt·h', tr.fpf_load_Npm, p.Yt * 8 * p.t, 1e-6, 'Uniaxial strength'),
      N.check('Cross-ply first-ply failure at the 90° transverse strain', cp.fpf_strain, lamReserve(Lc, [1, 0, 0, 0, 0, 0], 0, 0, 'Maximum stress').rf * Lc.abd[0][0], 1e-12, 'Consistency of strain and load'),
      N.check('Cross-ply ultimate: 0° plies at fibre failure strain, 90° plies discounted', cp.ult_load_Npm, crossPlyUlt(p, Qd, e0), 2e-3, 'Ply-discount (netting) analysis with Poisson coupling solved independently'),
    ];
  },
  calibration: { params: [{ key: 'kd', min: 0.001, max: 1 }], sweep: 'n_block', target: 'ult_strength_Pa', note: 'Unnotched tension coupon strengths calibrate the matrix-failure stiffness retention.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.min_RF < i.sf && o.RF_last_ply >= i.sf && o.fpf_benign) out.push({ severity: 'advise', title: 'First failure is matrix cracking with reserve to collapse', detail: `First-ply RF ${o.min_RF.toFixed(2)}, last-ply RF ${o.RF_last_ply.toFixed(2)}.`, action: 'Acceptable only if matrix cracking at ultimate load is permitted (no fuel-tightness, fatigue or moisture-ingress concern); otherwise add plies or reorient.', basis: 'Progressive ply-discount analysis' });
    if (o.RF_last_ply < i.sf) out.push({ severity: o.RF_last_ply < 1 ? 'critical' : 'warn', title: 'Laminate collapse below the required load', detail: `Last-ply reserve factor ${o.RF_last_ply.toFixed(2)} against ${i.sf}.`, action: 'Add fibre in the load direction or increase thickness.', basis: 'Ultimate load ≥ limit × ultimate factor' });
    if (o.lpf_over_fpf < 1.05) out.push({ severity: 'info', title: 'Fibre-dominated, brittle failure', detail: 'Collapse follows first-ply failure almost immediately.', action: 'Stress concentrations will not redistribute: apply notched (open-hole) allowables and consider ±45° plies for damage tolerance.', basis: 'Last-ply / first-ply load ratio' });
    out.push({ severity: 'info', title: 'Sustainability lever', detail: `Unnotched strength ${(o.ult_strength_Pa / 1e6).toFixed(0)} MPa at ${(plyOf(i.ply).rho / 1000).toFixed(2)} g/cm³.`, action: 'Tailoring ply percentages to the load path typically saves 10–20% of cover mass against a quasi-isotropic lay-up; weigh this against the higher embodied energy and end-of-life cost of carbon fibre in Suite 26.', basis: 'Specific strength; life-cycle trade' });
    return out;
  },
};
/** Independent check: ultimate Nx of [0/90]s with matrix-discounted 90° plies when the 0° plies reach their fibre failure stress (max-stress fibre mode of Hashin with zero shear). */
function crossPlyUlt(p, Qd, e0unused) {
  // unknowns: εx, εy with Ny = 0. 0° ply: σx = Q11εx + Q12εy; 90° ply (discounted): σx = Qd22 εx + Qd12 εy, σy = Qd12 εx + Qd11 εy
  const Q = Qof(p), t = p.t, A11 = 2 * t * (Q[0] + Qd[2]), A12 = 2 * t * (Q[1] + Qd[1]), A22 = 2 * t * (Q[2] + Qd[0]);
  // per unit Nx: εy = −A12/A22 εx, εx = 1/(A11 − A12²/A22); fibre stress in the 0° ply per unit Nx
  const ex = 1 / (A11 - (A12 * A12) / A22), ey = (-A12 / A22) * ex, s1 = Q[0] * ex + Q[1] * ey;
  return p.Xt / s1;
}

// ---- 3. micromechanics ----------------------------------------------------------------------
/** Halpin–Tsai: P/Pm = (1 + ξηVf)/(1 − ηVf), η = (Pf/Pm − 1)/(Pf/Pm + ξ). */
const halpinTsai = (Pf, Pm, Vf, xi) => { const eta = (Pf / Pm - 1) / (Pf / Pm + xi); return (Pm * (1 + xi * eta * Vf)) / (1 - eta * Vf); };
function micro(i, Vf) {
  const Vv = i.Vv, Vm = 1 - Vf - Vv, E1 = i.Ef1 * Vf + i.Em * Vm, nu12 = i.nuf * Vf + i.num * Vm, Gm = i.Em / (2 * (1 + i.num));
  return { E1, nu12, E2_rom: 1 / (Vf / i.Ef2 + (Vm + Vv) / i.Em), E2: halpinTsai(i.Ef2, i.Em, Vf, i.xiE) * (1 - Vv) ** 2, G12_rom: 1 / (Vf / i.Gf + (Vm + Vv) / Gm), G12: halpinTsai(i.Gf, Gm, Vf, i.xiG) * (1 - Vv) ** 2, rho: i.rho_f * Vf + i.rho_m * Vm, Xt: i.Sf * (Vf + (Vm * i.Em) / i.Ef1), Wf: (i.rho_f * Vf) / (i.rho_f * Vf + i.rho_m * Vm), Gm };
}
const micromech = {
  id: 'micromech', title: 'Micromechanics: ply properties from fibre and matrix', fidelity: 'analytical',
  summary: 'Rule-of-mixtures and Halpin–Tsai estimates of ply stiffness, density and fibre-direction strength from constituent properties and fibre volume fraction, with a void knock-down.',
  equations: ['Generalized Hooke’s law', 'Orthotropic constitutive equations'],
  inputs: [
    { key: 'Ef1', label: 'Fibre axial modulus', unit: 'Pa', default: 276e9, min: 1e9, group: 'Fibre', help: '230 GPa standard-modulus carbon, 276–290 GPa intermediate-modulus, 72 GPa E-glass' },
    { key: 'Ef2', label: 'Fibre transverse modulus', unit: 'Pa', default: 15e9, min: 1e9, group: 'Fibre', help: 'Carbon fibres are strongly anisotropic (10–20 GPa); glass is isotropic' },
    { key: 'Gf', label: 'Fibre shear modulus G12', unit: 'Pa', default: 20e9, min: 1e8, group: 'Fibre' },
    { key: 'nuf', label: 'Fibre Poisson ratio', unit: '-', default: 0.2, min: 0, max: 0.45, group: 'Fibre' },
    { key: 'rho_f', label: 'Fibre density', unit: 'kg/m³', default: 1780, min: 500, group: 'Fibre' },
    { key: 'Sf', label: 'Fibre tensile strength', unit: 'Pa', default: 5.0e9, min: 1e8, group: 'Fibre' },
    { key: 'Em', label: 'Matrix modulus', unit: 'Pa', default: 4.0e9, min: 1e8, group: 'Matrix' },
    { key: 'num', label: 'Matrix Poisson ratio', unit: '-', default: 0.35, min: 0.1, max: 0.49, group: 'Matrix' },
    { key: 'rho_m', label: 'Matrix density', unit: 'kg/m³', default: 1300, min: 500, group: 'Matrix' },
    { key: 'Vf', label: 'Fibre volume fraction', unit: '-', default: 0.6, min: 0.05, max: 0.8, group: 'Composite' },
    { key: 'Vv', label: 'Void volume fraction', unit: '-', default: 0, min: 0, max: 0.1, group: 'Composite', help: 'Autoclave prepreg is typically below 1%; out-of-autoclave and wet lay-up can be 2–5%' },
    { key: 'xiE', label: 'Halpin–Tsai ξ for E2', unit: '-', default: 2, min: 0.01, max: 100, group: 'Composite', help: 'Empirical reinforcement factor: 2 for round fibres in a square array' },
    { key: 'xiG', label: 'Halpin–Tsai ξ for G12', unit: '-', default: 1, min: 0.01, max: 100, group: 'Composite' },
    { key: 'ply', label: 'Compare with database ply', type: 'select', options: PLY_NAMES, default: 'IM7/8552 carbon-epoxy', group: 'Composite' },
  ],
  defaults(c) { const glass = /glass/i.test(c.struct.ply), kev = /Kevlar/i.test(c.struct.ply); return { ply: c.struct.ply, ...(glass ? { Ef1: 72e9, Ef2: 72e9, Gf: 30e9, nuf: 0.22, rho_f: 2550, Sf: 3.4e9, Vf: 0.5 } : kev ? { Ef1: 124e9, Ef2: 7e9, Gf: 3e9, nuf: 0.36, rho_f: 1440, Sf: 3.6e9 } : {}) }; },
  run(i) {
    const r = micro(i, i.Vf), db = plyOf(i.ply), vs = N.linspace(0.2, 0.75, 45), cur = vs.map((v) => micro(i, Math.min(v, 1 - i.Vv - 1e-6))), warnings = [];
    if (i.Vf + i.Vv > 0.8) warnings.push('Fibre volume fractions above about 0.7–0.75 cannot be consolidated without dry fibre contact; transverse and shear properties will fall short of these estimates.');
    if (i.Vv > 0.02) warnings.push('Void content above 2%: interlaminar shear strength typically falls several percent per 1% of voids (empirical); the stiffness knock-down used here does not capture that.');
    const pct = (a, b) => (100 * (a - b)) / b;
    return {
      kpis: [
        kpi('E1_Pa', 'Longitudinal modulus E1 (rule of mixtures)', r.E1, 'Pa'), kpi('E2_Pa', 'Transverse modulus E2 (Halpin–Tsai)', r.E2, 'Pa'), kpi('G12_Pa', 'Shear modulus G12 (Halpin–Tsai)', r.G12, 'Pa'), kpi('nu12', 'Major Poisson ratio ν12', r.nu12, '-'),
        kpi('E2_rom_Pa', 'E2 lower bound (inverse rule of mixtures)', r.E2_rom, 'Pa'), kpi('G12_rom_Pa', 'G12 lower bound (inverse rule of mixtures)', r.G12_rom, 'Pa'),
        kpi('rho_kgm3', 'Composite density', r.rho, 'kg/m³'), kpi('Wf', 'Fibre mass fraction', r.Wf, '-'), kpi('Xt_Pa', 'Fibre-direction tensile strength estimate', r.Xt, 'Pa', undefined, 'Upper bound: assumes all fibres reach their mean strength together'),
        kpi('specific_E1', 'Specific stiffness E1/ρ', r.E1 / r.rho, 'm²/s²'),
        kpi('dE1_vs_db_pct', `E1 versus database ply`, pct(r.E1, db.E1), '%'), kpi('dE2_vs_db_pct', 'E2 versus database ply', pct(r.E2, db.E2), '%'), kpi('dG12_vs_db_pct', 'G12 versus database ply', pct(r.G12, db.G12), '%'),
      ],
      plots: [
        { type: 'line', title: 'Fibre-direction modulus versus fibre volume fraction', xlabel: 'Fibre volume fraction [-]', ylabel: 'E1 [GPa]', series: [{ name: 'Rule of mixtures', x: vs, y: cur.map((c) => c.E1 / 1e9) }, { name: 'Database ply', x: [i.Vf], y: [db.E1 / 1e9], style: 'points' }] },
        { type: 'line', title: 'Matrix-dominated moduli versus fibre volume fraction', xlabel: 'Fibre volume fraction [-]', ylabel: 'Modulus [GPa]', series: [{ name: 'E2 Halpin–Tsai', x: vs, y: cur.map((c) => c.E2 / 1e9) }, { name: 'E2 inverse rule of mixtures', x: vs, y: cur.map((c) => c.E2_rom / 1e9), style: 'dash' }, { name: 'G12 Halpin–Tsai', x: vs, y: cur.map((c) => c.G12 / 1e9) }, { name: 'G12 inverse rule of mixtures', x: vs, y: cur.map((c) => c.G12_rom / 1e9), style: 'dash' }, { name: 'Database E2, G12', x: [i.Vf, i.Vf], y: [db.E2 / 1e9, db.G12 / 1e9], style: 'points' }] },
      ],
      warnings,
      models: ['Voigt rule of mixtures (E1, ν12, density, Xt)', 'Reuss inverse rule of mixtures (lower bounds)', 'Halpin–Tsai semi-empirical equations (E2, G12)', 'Void knock-down (1 − Vv)² on matrix-dominated moduli (empirical)'],
      assumptions: ['Perfect fibre–matrix bond, straight aligned fibres, uniform packing', 'Default constituent values are generic; the comparison ply is typical handbook data', 'Strength from micromechanics is an upper-bound estimate, not an allowable'],
    };
  },
  verify() {
    const i = { Ef1: 230e9, Ef2: 15e9, Gf: 20e9, nuf: 0.2, rho_f: 1780, Sf: 4e9, Em: 3.5e9, num: 0.35, rho_m: 1250, Vf: 0.6, Vv: 0, xiE: 2, xiG: 1 }, Gm = 3.5e9 / 2.7;
    return [
      N.check('Rule of mixtures E1 = Ef·Vf + Em·Vm', micro(i, 0.6).E1, 230e9 * 0.6 + 3.5e9 * 0.4, 1e-12, 'Voigt bound'),
      N.check('Halpin–Tsai → matrix modulus at Vf = 0', micro(i, 0).E2, 3.5e9, 1e-12, 'Limit case'),
      N.check('Halpin–Tsai → fibre modulus at Vf = 1', halpinTsai(15e9, 3.5e9, 1, 2), 15e9, 1e-12, 'Limit case'),
      N.check('Halpin–Tsai with ξ → 0 equals the inverse rule of mixtures', halpinTsai(20e9, Gm, 0.6, 1e-12), 1 / (0.6 / 20e9 + 0.4 / Gm), 1e-9, 'Halpin & Kardos (1976)'),
      N.check('Halpin–Tsai with ξ → ∞ equals the rule of mixtures', halpinTsai(20e9, Gm, 0.6, 1e12), 0.6 * 20e9 + 0.4 * Gm, 1e-9, 'Halpin & Kardos (1976)'),
    ];
  },
  calibration: { params: [{ key: 'xiE', min: 0.01, max: 100 }, { key: 'Ef2', min: 1e9, max: 1e11 }], sweep: 'Vf', target: 'E2_Pa', note: 'Measured transverse modulus at several fibre volume fractions calibrates ξ and the fibre transverse modulus.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (Math.abs(o.dE1_vs_db_pct) > 10) out.push({ severity: 'advise', title: 'Constituent data do not reproduce the database ply', detail: `E1 differs by ${o.dE1_vs_db_pct.toFixed(0)}% from ${i.ply}.`, action: 'Set the fibre modulus and volume fraction for the actual prepreg before using these estimates in laminate analysis.', basis: 'Rule of mixtures is accurate to a few percent for E1' });
    if (i.Vv > 0.01) out.push({ severity: 'warn', title: 'Void content needs process control', detail: `${(100 * i.Vv).toFixed(1)}% voids.`, action: 'Improve debulking/cure pressure or accept a matrix-dominated strength knock-down verified by short-beam-shear coupons.', basis: 'Typical acceptance limit of 1–2% voids for primary structure' });
    out.push({ severity: 'info', title: 'Stiffness-to-mass lever', detail: `Specific stiffness ${(o.specific_E1 / 1e6).toFixed(0)} MN·m/kg at Vf = ${i.Vf}.`, action: 'Raising fibre volume fraction from 0.55 to 0.60 adds about 9% to E1 for about 2% more density; mass saved in stiffness-critical parts feeds straight into fuel or battery energy.', basis: 'Rule of mixtures' });
    return out;
  },
};

// ---- 4. interlaminar fracture and shear -----------------------------------------------------
/** Energy release rates [J/m²] of the standard delamination specimens by simple beam theory; h is the arm (half-laminate) thickness. */
function gRates(test, P, a, b, h, E, Lh, c) {
  const k = b * b * h ** 3 * E;
  if (test === 'DCB') return { GI: (12 * P * P * a * a) / k, GII: 0 };
  if (test === 'ENF') return { GI: 0, GII: (9 * P * P * a * a) / (16 * k) };
  const PI = (P * (3 * c - Lh)) / (4 * Lh), PII = (P * (c + Lh)) / Lh;
  return { GI: PI > 0 ? (12 * PI * PI * a * a) / k : 0, GII: (9 * PII * PII * a * a) / (16 * k) };
}
const bk = (GIc, GIIc, mix, eta) => GIc + (GIIc - GIc) * mix ** eta;
/** Interlaminar shear through a laminate strip in cylindrical bending: τxz(z) from equilibrium and the FSDT shear-correction factor. */
function shearProfile(L, G13, G23, Q) {
  const lay = L.lay, E = lay.map((l) => l.Qb[0][0]), EA = N.sum(lay.map((l, k) => E[k] * (l.z1 - l.z0))), zn = N.sum(lay.map((l, k) => (E[k] * (l.z1 * l.z1 - l.z0 * l.z0)) / 2)) / EA;
  const Db = N.sum(lay.map((l, k) => (E[k] * ((l.z1 - zn) ** 3 - (l.z0 - zn) ** 3)) / 3)), zs = [], tau = []; let acc = 0, U = 0, GA = 0, tmax = 0;
  for (let k = lay.length - 1; k >= 0; k--) {
    const l = lay[k], c = Math.cos(N.rad(l.th)), s = Math.sin(N.rad(l.th)), G = 1 / ((c * c) / G13 + (s * s) / G23), f = (z) => acc + (E[k] * ((l.z1 - zn) ** 2 - (z - zn) ** 2)) / 2;
    for (let j = 0; j <= 6; j++) { const z = l.z1 - ((l.z1 - l.z0) * j) / 6, t = (Q * f(z)) / Db; zs.push(z * 1e3); tau.push(t); tmax = Math.max(tmax, Math.abs(t)); }
    for (const [g, w] of [[0.5 - Math.sqrt(0.15), 5 / 18], [0.5, 4 / 9], [0.5 + Math.sqrt(0.15), 5 / 18]]) U += (w * (l.z1 - l.z0) * ((Q * f(l.z0 + g * (l.z1 - l.z0))) / Db) ** 2) / G; // 3-point Gauss: exact for the quartic integrand
    GA += G * (l.z1 - l.z0); acc = f(l.z0);
  }
  return { zs, tau, tmax, k: Q ? (Q * Q) / (U * GA) : NaN, GA, closure: acc / (EA * L.h) };
}
const interlaminar = {
  id: 'interlaminar', title: 'Delamination and interlaminar shear', fidelity: 'analytical',
  summary: 'Mode I, mode II and mixed-mode energy release rates of standard delamination specimens by beam theory with the Benzeggagh–Kenane criterion, the delamination growth load–displacement curve, and the through-thickness shear stress and shear-correction factor of the laminate.',
  equations: ['First-order shear deformation theory', 'Classical laminate theory equations', 'Continuum damage evolution equations'],
  inputs: [...LAM_INPUTS,
    { key: 'test', label: 'Delamination configuration', type: 'select', options: ['DCB', 'ENF', 'MMB'], default: 'DCB', group: 'Delamination', help: 'DCB: mode I opening; ENF: mode II sliding (three-point bend); MMB: mixed mode' },
    { key: 'P', label: 'Applied load', unit: 'N', default: 40, min: 0, group: 'Delamination' },
    { key: 'a0', label: 'Delamination length', unit: 'm', default: 0.05, min: 0.002, group: 'Delamination' },
    { key: 'b', label: 'Specimen width', unit: 'm', default: 0.025, min: 0.002, group: 'Delamination' },
    { key: 'h_arm', label: 'Arm thickness (half laminate)', unit: 'm', default: 0.0015, min: 1e-4, group: 'Delamination' },
    { key: 'Lh', label: 'Half-span (ENF, MMB)', unit: 'm', default: 0.05, min: 0.005, group: 'Delamination' },
    { key: 'c_lever', label: 'MMB lever length', unit: 'm', default: 0.04, min: 0.005, group: 'Delamination' },
    { key: 'GIc', label: 'Mode I toughness GIc (0 = database)', unit: 'J/m²', default: 0, min: 0, group: 'Toughness' },
    { key: 'GIIc', label: 'Mode II toughness GIIc (0 = database)', unit: 'J/m²', default: 0, min: 0, group: 'Toughness' },
    { key: 'eta', label: 'Benzeggagh–Kenane exponent η', unit: '-', default: 2, min: 0.5, max: 5, group: 'Toughness', help: 'Fitted to mixed-mode bending tests; roughly 1.5–2.5 for carbon/epoxy' },
    { key: 'Q', label: 'Transverse shear force per unit width', unit: 'N/m', default: 20000, group: 'Interlaminar shear' },
    { key: 'G23', label: 'Ply transverse shear modulus G23', unit: 'Pa', default: 3.3e9, min: 1e8, group: 'Interlaminar shear', help: 'About E2/(2(1+ν23)) with ν23 ≈ 0.4–0.5' },
  ],
  defaults(c) { const p = plyOf(c.struct.ply); return { layup: c.struct.layup, ply: c.struct.ply, G23: p.E2 / 2.9, h_arm: 12 * p.t, P: 60 }; },
  run(i) {
    const warnings = [], { p, L } = lamFromInputs(i, warnings), GIc = i.GIc || p.GIc, GIIc = i.GIIc || p.GIIc, E = p.E1;
    const g = gRates(i.test, i.P, i.a0, i.b, i.h_arm, E, i.Lh, i.c_lever), GT = g.GI + g.GII, g1 = gRates(i.test, 1, i.a0, i.b, i.h_arm, E, i.Lh, i.c_lever), mix = g1.GI + g1.GII > 0 ? g1.GII / (g1.GI + g1.GII) : 0;
    const Gc = bk(GIc, GIIc, mix, i.eta), fi = GT / Gc, Pc = Math.sqrt(Gc / (g1.GI + g1.GII));
    if (i.test !== 'DCB' && i.a0 > i.Lh) warnings.push('The delamination extends beyond the half-span: the ENF/MMB beam formulas no longer apply.');
    if (i.test === 'MMB' && 3 * i.c_lever < i.Lh) warnings.push('Lever shorter than L/3: the mode I component is closing, so the specimen is in pure mode II with crack-face contact.');
    if (i.a0 / i.h_arm < 10) warnings.push('Short, thick arms (a/h < 10): simple beam theory underestimates compliance; apply the corrected beam theory root-rotation term.');
    // DCB-type propagation curve: load–opening with the delamination growing at G = Gc
    const I = (i.b * i.h_arm ** 3) / 12, as = N.linspace(i.a0, 3 * i.a0, 40), comp = (a) => (i.test === 'DCB' ? (2 * a ** 3) / (3 * E * I) : i.test === 'ENF' ? (2 * i.Lh ** 3 + 3 * Math.min(a, i.Lh) ** 3) / (8 * E * i.b * i.h_arm ** 3) : NaN);
    const pcOf = (a) => { const q = gRates(i.test, 1, a, i.b, i.h_arm, E, i.Lh, i.c_lever); return Math.sqrt(Gc / (q.GI + q.GII)); }, plots = [];
    if (i.test !== 'MMB') { const aa = i.test === 'ENF' ? as.filter((a) => a <= i.Lh) : as; plots.push({ type: 'line', title: 'Load–displacement with delamination growth', xlabel: 'Load-point displacement [mm]', ylabel: 'Load [N]', series: [{ name: 'Elastic loading then growth at G = Gc', x: [0, ...aa.map((a) => pcOf(a) * comp(a) * 1e3)], y: [0, ...aa.map(pcOf)] }], annotations: [{ y: i.P, label: 'Applied' }] }); }
    plots.push({ type: 'line', title: 'Critical load versus delamination length', xlabel: 'Delamination length [mm]', ylabel: 'Critical load [N]', series: [{ name: 'Onset of growth', x: as.map((a) => a * 1e3), y: as.map(pcOf) }], annotations: [{ y: i.P, label: 'Applied' }] });
    const mx = N.linspace(0, 1, 41);
    plots.push({ type: 'line', title: 'Mixed-mode fracture criterion', xlabel: 'Mode mixity GII/(GI + GII) [-]', ylabel: 'Critical energy release rate [J/m²]', series: [{ name: `Benzeggagh–Kenane, η = ${i.eta}`, x: mx, y: mx.map((v) => bk(GIc, GIIc, v, i.eta)) }, { name: 'Applied state', x: [mix], y: [GT], style: 'points' }] });
    // interlaminar shear of the laminate
    const sp = shearProfile(L, p.G12, i.G23, i.Q), rfS = sp.tmax > 0 ? p.S / sp.tmax : Infinity;
    plots.push({ type: 'line', title: 'Interlaminar shear stress through the thickness', xlabel: 'τxz [MPa]', ylabel: 'z [mm]', series: [{ name: 'τxz', x: sp.tau.map((v) => v / 1e6), y: sp.zs }] });
    return {
      kpis: [
        kpi('GI_Jm2', 'Mode I energy release rate', g.GI, 'J/m²'), kpi('GII_Jm2', 'Mode II energy release rate', g.GII, 'J/m²'), kpi('mode_mix', 'Mode mixity GII/GT', mix, '-'),
        kpi('Gc_Jm2', 'Mixed-mode toughness (B-K)', Gc, 'J/m²'),
        kpi('delam_index', 'Delamination growth index GT/Gc', fi, '-', fi < 0.67 ? 'ok' : fi < 1 ? 'warn' : 'bad', 'Growth is predicted at 1'),
        kpi('P_crit_N', 'Critical load for growth', Pc, 'N'),
        kpi('tau_max_Pa', 'Peak interlaminar shear stress', sp.tmax, 'Pa'),
        kpi('RF_ils', 'Interlaminar shear reserve factor', rfS, '-', rfS >= 1.5 ? 'ok' : rfS >= 1 ? 'warn' : 'bad', 'Against the in-plane shear strength S used as a proxy for interlaminar shear strength'),
        kpi('k_shear', 'FSDT shear-correction factor', sp.k, '-', undefined, '5/6 for a homogeneous section'),
        kpi('GA_shear_Npm', 'Transverse shear stiffness k·ΣG·t', sp.k * sp.GA, 'N/m'),
      ],
      plots, warnings,
      models: ['Simple beam theory for DCB, ENF and MMB specimens (Reeder–Crews load partition)', 'Benzeggagh–Kenane mixed-mode criterion', 'Equilibrium-derived interlaminar shear in cylindrical bending with energy-based shear correction'],
      assumptions: ['Unidirectional 0° arms with modulus E1 for the fracture specimens; no root rotation, large-displacement or R-curve (fibre-bridging) effects', 'Self-similar growth between the mid-plane plies', 'Interlaminar shear strength approximated by the in-plane shear strength S', DATA_NOTE],
    };
  },
  verify() {
    const E = 150e9, b = 0.025, h = 0.0015, a = 0.05, P = 50, Lh = 0.05, dC = (C) => (C(a * 1.0001) - C(a * 0.9999)) / (2e-4 * a), I = (b * h ** 3) / 12;
    const gd = gRates('DCB', P, a, b, h, E, Lh, 0.04), ge = gRates('ENF', P, a, b, h, E, Lh, 0.04), gm = gRates('MMB', P, a, b, h, E, Lh, Lh / 3);
    const iso = { E1: 70e9, E2: 70e9, G12: 26.9e9, nu12: 0.3, a1: 0, a2: 0 }, Li = laminate([0, 0, 0, 0], [iso, iso, iso, iso], 0.001), sp = shearProfile(Li, 26.9e9, 26.9e9, 1000);
    return [
      N.check('DCB: G = (P²/2b)·dC/da', gd.GI, ((P * P) / (2 * b)) * dC((x) => (2 * x ** 3) / (3 * E * I)), 1e-6, 'Irwin–Kies compliance relation'),
      N.check('ENF: G = (P²/2b)·dC/da', ge.GII, ((P * P) / (2 * b)) * dC((x) => (2 * Lh ** 3 + 3 * x ** 3) / (8 * E * b * h ** 3)), 1e-6, 'Irwin–Kies compliance relation'),
      N.check('MMB with c = L/3 is pure mode II', gm.GI, 0, 1e-12, 'Reeder & Crews (1990)'),
      N.check('B-K returns GIc in pure mode I and GIIc in pure mode II', bk(200, 600, 0, 2) + bk(200, 600, 1, 2), 800, 1e-12, 'Benzeggagh & Kenane (1996)'),
      N.check('Homogeneous section: peak shear 1.5·Q/h', sp.tmax, (1.5 * 1000) / 0.004, 1e-9, 'Jourawski shear formula'),
      N.check('Homogeneous section: shear-correction factor 5/6', sp.k, 5 / 6, 1e-6, 'Reissner (1945)'),
    ];
  },
  calibration: { params: [{ key: 'eta', min: 0.5, max: 5 }], sweep: 'c_lever', target: 'Gc_Jm2', note: 'Mixed-mode bending toughness measured at several lever lengths (mode mixities) calibrates the B-K exponent.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.delam_index >= 1) out.push({ severity: 'critical', title: 'Delamination growth predicted', detail: `GT/Gc = ${o.delam_index.toFixed(2)} at ${i.P} N.`, action: 'Reduce the peel load path (taper ply drops, add fasteners or z-pins at stiffener run-outs) or select a toughened resin system.', basis: 'Benzeggagh–Kenane criterion' });
    if (o.RF_ils < 1.5) out.push({ severity: 'warn', title: 'Low interlaminar shear reserve', detail: `Reserve factor ${o.RF_ils.toFixed(2)}.`, action: 'Thicken the laminate, spread the transverse load, or move stiff 0° plies away from the mid-plane where shear peaks.', basis: 'τxz against shear strength' });
    out.push({ severity: 'info', title: 'Inspection and life extension', detail: `Critical load for growth ${o.P_crit_N.toFixed(0)} N at the present delamination length.`, action: 'Use this no-growth threshold to set allowable damage limits: structures shown not to grow detectable delaminations can stay in service without repair, avoiding scrap of high-embodied-energy parts.', basis: 'No-growth damage-tolerance approach' });
    return out;
  },
};

// ---- 5. viscoelastic creep and Norton relaxation --------------------------------------------
/** Creep strain history of a Prony-series solid under a stress history sig(t), by the recursive (exponential) algorithm on the creep form. */
function pronyCreep(E0, terms, sig, tEnd, n) {
  // generalised Kelvin form: J(t) = 1/E0 + Σ Jk (1 − exp(−t/τk)); internal strains qk follow dqk/dt = (Jk σ − qk)/τk
  const dt = tEnd / n, q = terms.map(() => 0), t = [0], eps = [sig(0) / E0];
  for (let s = 1; s <= n; s++) {
    const t1 = s * dt, s0 = sig(t1 - dt), s1 = sig(t1);
    terms.forEach((k, j) => { const x = dt / k.tau, e = Math.exp(-x), g = x > 1e-8 ? (1 - e) / x : 1 - x / 2; q[j] = e * q[j] + k.J * (s0 * (g - e) + s1 * (1 - g)); }); // exact for piecewise-linear stress
    t.push(t1); eps.push(s1 / E0 + N.sum(q));
  }
  return { t, eps };
}
const visco = {
  id: 'visco', title: 'Creep and stress relaxation', fidelity: 'numerical',
  summary: 'Matrix-dominated viscoelastic creep of a composite under sustained load (standard linear solid / Prony series with temperature shift) and Norton power-law relaxation of a preloaded metal part such as a bolt.',
  equations: ['Creep constitutive equations', 'Viscoelastic constitutive equations'],
  inputs: [
    { key: 'sigma', label: 'Sustained stress', unit: 'Pa', default: 30e6, group: 'Viscoelastic', help: 'Matrix-dominated stress, e.g. in-plane shear or transverse stress in a ply' },
    { key: 'E0', label: 'Instantaneous modulus', unit: 'Pa', default: 5.29e9, min: 1e6, group: 'Viscoelastic', help: 'G12 or E2 of the ply' },
    { key: 'g1', label: 'Relaxing fraction, term 1', unit: '-', default: 0.1, min: 0, max: 0.9, group: 'Viscoelastic', help: 'Share of the modulus that relaxes with time constant τ₁ (illustrative default: fit to creep tests)' },
    { key: 'tau1', label: 'Relaxation time τ₁ at reference temperature', unit: 's', default: 3.6e5, min: 1e-3, group: 'Viscoelastic' },
    { key: 'g2', label: 'Relaxing fraction, term 2', unit: '-', default: 0.1, min: 0, max: 0.9, group: 'Viscoelastic' },
    { key: 'tau2', label: 'Relaxation time τ₂ at reference temperature', unit: 's', default: 3.6e7, min: 1e-3, group: 'Viscoelastic' },
    { key: 'T_K', label: 'Service temperature', unit: 'K', default: 343, min: 150, max: 700, group: 'Viscoelastic' },
    { key: 'Tref_K', label: 'Reference temperature', unit: 'K', default: 296, min: 150, max: 700, group: 'Viscoelastic' },
    { key: 'Ea', label: 'Activation energy (Arrhenius shift)', unit: 'J/mol', default: 1.0e5, min: 0, group: 'Viscoelastic', help: 'Time–temperature shift below Tg; illustrative default' },
    { key: 't_end_h', label: 'Duration', unit: 'h', default: 10000, min: 0.01, group: 'Viscoelastic' },
    { key: 's0', label: 'Initial preload stress (metal)', unit: 'Pa', default: 400e6, min: 1e5, group: 'Norton creep' },
    { key: 'E_metal', label: 'Metal modulus', unit: 'Pa', default: 113.8e9, min: 1e9, group: 'Norton creep' },
    { key: 'A_n', label: 'Norton coefficient A', unit: '1/(s·MPaⁿ)', default: 1e-24, min: 0, group: 'Norton creep', help: 'ε̇ = A·σⁿ·exp(−Q/RT) with σ in MPa. Strongly material- and temperature-specific: supply test data' },
    { key: 'n_n', label: 'Norton exponent n', unit: '-', default: 5, min: 1.01, max: 20, group: 'Norton creep' },
    { key: 'Q_n', label: 'Creep activation energy Q', unit: 'J/mol', default: 0, min: 0, group: 'Norton creep', help: '0 when A is already given at the service temperature' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 400, min: 20, max: 20000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c) { const p = plyOf(c.struct.ply), m = METALS[c.struct.material] || METALS['Al 2024-T3']; return { E0: p.G12, sigma: 0.3 * p.S, E_metal: m.E, s0: 0.6 * m.Sy }; },
  run(i) {
    const R = 8.314462, aT = Math.exp((i.Ea / R) * (1 / i.T_K - 1 / i.Tref_K)), tEnd = i.t_end_h * 3600, n = Math.max(20, Math.round(i.nSteps)), warnings = [];
    const gs = [[i.g1, i.tau1 * aT], [i.g2, i.tau2 * aT]].filter(([g]) => g > 0), gsum = N.sum(gs.map((g) => g[0]));
    if (gsum >= 0.98) warnings.push('The relaxing fractions sum to almost 1: the long-term modulus is near zero (fluid-like response).');
    // relaxation modulus E(t) = E0 (1 − Σ gk (1 − exp(−t/τk)))
    const Erel = (t) => i.E0 * (1 - N.sum(gs.map(([g, tau]) => g * (1 - Math.exp(-t / tau)))));
    const cr = relaxToCreep(i.E0, gs, i.sigma, tEnd, n), eps0 = i.sigma / i.E0, epsEnd = cr.eps[n], Einf = i.E0 * (1 - gsum);
    // Norton relaxation at fixed total strain: dσ/dt = −E·A·σⁿ·exp(−Q/RT), σ in MPa inside the power law
    const k = i.A_n * Math.exp(-i.Q_n / (R * i.T_K)), rate = (s) => -i.E_metal * k * (Math.max(s, 0) / 1e6) ** i.n_n, tg = N.range(n + 1, (j) => (tEnd * j) / n), sg = [i.s0];
    for (let j = 0; j < n; j++) { const h = tEnd / n, s = sg[j], k1 = rate(s), k2 = rate(s + 0.5 * h * k1), k3 = rate(s + 0.5 * h * k2), k4 = rate(s + h * k3); sg.push(Math.max(0, s + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4))); }
    const sEnd = sg[n], loss = 1 - sEnd / i.s0, s90 = 0.9 * i.s0, c0 = (i.s0 / 1e6) ** (1 - i.n_n), t90 = k > 0 ? (((s90 / 1e6) ** (1 - i.n_n) - c0) * 1e6) / ((i.n_n - 1) * i.E_metal * k) : Infinity;
    if (epsEnd / eps0 > 1.5) warnings.push('Creep strain exceeds 150% of the elastic strain: linear viscoelasticity is doubtful at this stress and temperature.');
    const th = cr.t.map((t) => t / 3600), idx = N.range(Math.min(n + 1, 201), (j) => Math.round((j * n) / Math.min(n, 200)));
    return {
      kpis: [
        kpi('creep_strain', 'Total strain at end of duration', epsEnd, '-'),
        kpi('creep_ratio', 'Strain / instantaneous elastic strain', epsEnd / eps0, '-', epsEnd / eps0 < 1.15 ? 'ok' : 'warn'),
        kpi('E_relaxed_Pa', 'Relaxation modulus at end of duration', Erel(tEnd), 'Pa'),
        kpi('E_inf_Pa', 'Long-term (fully relaxed) modulus', Einf, 'Pa'),
        kpi('shift_factor', 'Time–temperature shift factor a_T', aT, '-', undefined, 'Below 1: creep is faster than at the reference temperature'),
        kpi('preload_end_Pa', 'Metal preload stress at end of duration', sEnd, 'Pa'),
        kpi('preload_loss_pct', 'Preload loss', 100 * loss, '%', loss < 0.1 ? 'ok' : loss < 0.25 ? 'warn' : 'bad'),
        kpi('t_10pct_loss_h', 'Time to 10% preload loss', t90 / 3600, 'h'),
      ],
      plots: [
        { type: 'line', title: 'Viscoelastic creep under sustained stress', xlabel: 'Time [h]', ylabel: 'Strain [%]', series: [{ name: 'Total strain', x: idx.map((j) => th[j]), y: idx.map((j) => 100 * cr.eps[j]) }], annotations: [{ y: 100 * eps0, label: 'Elastic' }, { y: (100 * i.sigma) / Math.max(Einf, 1), label: 'Long-term limit' }] },
        { type: 'line', title: 'Relaxation modulus', xlabel: 'Time [h]', ylabel: 'Modulus [GPa]', xlog: true, series: [{ name: 'E(t)', x: N.logspace(0.01, Math.max(1, i.t_end_h), 60), y: N.logspace(0.01, Math.max(1, i.t_end_h), 60).map((t) => Erel(t * 3600) / 1e9) }] },
        { type: 'line', title: 'Norton-law preload relaxation', xlabel: 'Time [h]', ylabel: 'Stress [MPa]', series: [{ name: 'Preload stress', x: idx.map((j) => tg[j] / 3600), y: idx.map((j) => sg[j] / 1e6) }] },
      ],
      warnings,
      models: ['Prony-series (generalised Maxwell) linear viscoelastic solid', 'Arrhenius time–temperature shift', 'Crank–Nicolson integration of the internal-variable equations on a geometric time grid', 'Norton power-law creep, RK4 relaxation at fixed strain'],
      assumptions: ['Linear, thermorheologically simple viscoelasticity below the glass transition', 'Default Prony and Norton parameters are illustrative placeholders, not material data: fit them to creep or relaxation tests', 'Fibre-direction response is treated as elastic; primary and tertiary creep of metals are not modelled'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [50, 100, 200, 400, 800], metric: 'creep_strain' },
  verify() {
    // standard linear solid: creep ε(t) = σ/E∞ − σ(1/E∞ − 1/E0) exp(−t/τc), τc = τ·E0/E∞
    const E0 = 5e9, g = 0.3, tau = 1000, sig = 1e7, t = 2000, Einf = E0 * (1 - g), exact = sig / Einf - sig * (1 / Einf - 1 / E0) * Math.exp((-t * Einf) / (tau * E0));
    const num = relaxToCreep(E0, [[g, tau]], sig, t, 800).eps.pop(), kel = pronyCreep(E0, [{ J: 1 / Einf - 1 / E0, tau: (tau * E0) / Einf }], () => sig, t, 10).eps.pop();
    const b = { sigma: 1e7, E0: 5e9, g1: 0, tau1: 1, g2: 0, tau2: 1, T_K: 300, Tref_K: 300, Ea: 0, t_end_h: 1000, s0: 300e6, E_metal: 100e9, A_n: 1e-20, n_n: 4, Q_n: 0, nSteps: 400 }, r = N.kv(visco.run(b));
    const sExact = 1e6 * ((300 ** -3) + 3 * (100e9 / 1e6) * 1e-20 * 3.6e6) ** (-1 / 3);
    return [
      N.check('Standard linear solid creep (internal-variable Crank–Nicolson)', num, exact, 1e-4, 'Closed-form Zener creep compliance'),
      N.check('Standard linear solid creep (recursive Kelvin algorithm)', kel, exact, 1e-10, 'Closed-form Zener creep compliance; the recursion is exact for constant stress'),
      N.check('Norton relaxation σ^(1−n) = σ0^(1−n) + (n−1)·E·A·t', r.preload_end_Pa, sExact, 1e-6, 'Closed-form integration of dσ/dt = −E·A·σⁿ'),
      N.check('No relaxing terms: strain stays elastic', r.creep_ratio, 1, 1e-12, 'Limit case'),
    ];
  },
  calibration: { params: [{ key: 'g1', min: 0, max: 0.9 }, { key: 'tau1', min: 1, max: 1e9 }, { key: 'g2', min: 0, max: 0.9 }], sweep: 't_end_h', target: 'creep_strain', note: 'Creep strain versus time from a constant-load test at the service temperature calibrates the Prony terms.' },
  recommend(res) {
    const o = res.outputs, out = [];
    if (o.creep_ratio > 1.15) out.push({ severity: 'advise', title: 'Noticeable matrix creep', detail: `Strain grows to ${(100 * (o.creep_ratio - 1)).toFixed(0)}% above the elastic value.`, action: 'Carry sustained loads in fibre-dominated directions, lower the service temperature margin to Tg, or use the relaxed modulus for long-term deflection and bolted-joint clamp-up.', basis: 'Linear viscoelastic creep' });
    if (o.preload_loss_pct > 10) out.push({ severity: 'warn', title: 'Preload relaxation', detail: `${o.preload_loss_pct.toFixed(0)}% of preload lost; 10% is reached after ${Number.isFinite(o.t_10pct_loss_h) ? o.t_10pct_loss_h.toFixed(0) : '∞'} h.`, action: 'Re-torque at that interval, reduce the preload stress, or choose a more creep-resistant alloy for hot joints.', basis: 'Norton power-law relaxation' });
    out.push({ severity: 'info', title: 'Parameters need test data', detail: 'Default creep parameters are placeholders.', action: 'Upload creep or relaxation curves and use calibration; time–temperature superposition lets short hot tests stand in for long service lives, saving test time and energy.', basis: 'Model calibration' });
    return out;
  },
};
/**
 * Creep strain under constant stress of a generalised Maxwell (Prony) solid, gs = [[g, τ], ...]. The internal stresses q_k obey the
 * linear system q̇_k = (E_k/E0)·Σ q_j/τ_j − q_k/τ_k and ε̇ = Σ(q_j/τ_j)/E0; integrated by Crank–Nicolson on a geometric time grid
 * so that relaxation times far shorter than the duration are still resolved.
 */
function relaxToCreep(E0, gs, sig, tEnd, n) {
  const m = gs.length, t = [0], eps = [sig / E0]; if (!m) { for (let j = 1; j <= n; j++) { t.push((tEnd * j) / n); eps.push(sig / E0); } return { t, eps }; }
  const A = N.range(m, (r) => N.range(m, (c) => gs[r][0] / gs[c][1] - (r === c ? 1 / gs[r][1] : 0))), t0 = Math.min(tEnd / n, Math.min(...gs.map((g) => g[1])) / 20), ratio = (tEnd / t0) ** (1 / Math.max(1, n - 1));
  let q = gs.map(() => 0), e = sig / E0, tc = 0; // at t = 0+ the springs carry σ·g_k each relative to E0·ε0
  q = gs.map(([g]) => g * sig);
  for (let j = 1; j <= n; j++) {
    const tn = j === n ? tEnd : t0 * ratio ** (j - 1), h = tn - tc, L = N.range(m, (r) => N.range(m, (c) => (r === c ? 1 : 0) - 0.5 * h * A[r][c])), rhs = q.map((v, r) => v + 0.5 * h * N.dot(A[r], q)), qn = N.solve(L, rhs);
    e += (h * N.sum(gs.map(([, tau], k) => (0.5 * (q[k] + qn[k])) / tau))) / E0; q = qn; tc = tn; t.push(tn); eps.push(e);
  }
  return { t, eps };
}

// ---- 6. moisture diffusion and corrosion pitting --------------------------------------------
/** 1-D Fickian diffusion through thickness h by Crank–Nicolson; both faces (or one) held at the saturation level. Returns average content history. */
function fick(h, D, c0, cs, tEnd, nz, nt, twoSided) {
  const dz = h / nz, dt = tEnd / nt, r = (D * dt) / (2 * dz * dz), n = nz + 1; let c = new Array(n).fill(c0); c[0] = cs; if (twoSided) c[nz] = cs;
  const avg = (v) => (N.sum(v) - 0.5 * (v[0] + v[nz])) / nz, t = [0], M = [c0], prof = [];
  const a = new Array(n).fill(-r), b = new Array(n).fill(1 + 2 * r), cc = new Array(n).fill(-r);
  b[0] = 1; cc[0] = 0; if (twoSided) { a[nz] = 0; b[nz] = 1; } else { a[nz] = -2 * r; } // insulated far face by mirror node
  for (let s = 1; s <= nt; s++) {
    const d = c.map((v, k) => (k === 0 ? cs : k === nz ? (twoSided ? cs : v + 2 * r * (c[nz - 1] - v)) : v + r * (c[k - 1] - 2 * v + c[k + 1])));
    c = N.solveTridiag(a, b, cc, d); t.push(s * dt); M.push(s === 0 ? c0 : avg(c));
    if (s === Math.round(nt / 20) || s === Math.round(nt / 4) || s === nt) prof.push({ t: s * dt, c: c.slice() });
  }
  return { t, M, prof, z: N.range(n, (k) => k * dz) };
}
/** Series solution for the average moisture content of a plate exposed on both faces. */
const fickExact = (h, D, c0, cs, t) => { let s = 0; for (let j = 0; j < 200; j++) s += Math.exp((-((2 * j + 1) ** 2) * Math.PI ** 2 * D * t) / (h * h)) / (2 * j + 1) ** 2; return cs - (cs - c0) * (8 / Math.PI ** 2) * s; };
const environment = {
  id: 'environment', title: 'Moisture uptake, hot-wet knock-down and corrosion pitting', fidelity: 'numerical',
  summary: 'Fickian moisture diffusion into a laminate over years of service, the resulting glass-transition and strength knock-downs, and power-law growth of corrosion pits in a metal part with the remaining section and equivalent flaw size.',
  equations: ['Continuum damage evolution equations'],
  inputs: [
    { key: 'h_m', label: 'Laminate thickness', unit: 'm', default: 0.004, min: 1e-4, group: 'Moisture', help: 'Defaults to the laminate from the stiffness analysis' },
    { key: 'D_m2s', label: 'Through-thickness diffusivity', unit: 'm²/s', default: 2e-13, min: 1e-16, max: 1e-9, group: 'Moisture', help: 'Order 10⁻¹³ m²/s for carbon/epoxy at room temperature, rising steeply with temperature; measure for your system' },
    { key: 'rh', label: 'Relative humidity of the environment', unit: '-', default: 0.85, min: 0, max: 1, group: 'Moisture' },
    { key: 'M_sat100', label: 'Saturation moisture content at 100% RH', unit: '%', default: 1.5, min: 0, max: 8, group: 'Moisture', help: 'Mass % of the laminate; about 1–2% for aerospace epoxies' },
    { key: 'M0', label: 'Initial moisture content', unit: '%', default: 0, min: 0, max: 8, group: 'Moisture' },
    { key: 'years', label: 'Exposure time', unit: 'yr', default: 20, min: 0.01, max: 100, group: 'Moisture' },
    { key: 'sides', label: 'Exposed faces', type: 'select', options: ['both', 'one'], default: 'both', group: 'Moisture' },
    { key: 'Tg_dry_K', label: 'Dry glass-transition temperature', unit: 'K', default: 473, min: 300, max: 700, group: 'Knock-down (empirical)' },
    { key: 'dTg_per_pct', label: 'Tg loss per 1% moisture', unit: 'K/%', default: 20, min: 0, max: 60, group: 'Knock-down (empirical)', help: 'Empirical; of order 15–30 K per 1% for epoxies. Replace with wet DMA data' },
    { key: 'k_str_per_pct', label: 'Matrix-dominated strength loss per 1% moisture', unit: '%/%', default: 8, min: 0, max: 40, group: 'Knock-down (empirical)', help: 'Illustrative placeholder for hot-wet compression/shear knock-down; replace with coupon data' },
    { key: 'T_service_K', label: 'Maximum service temperature', unit: 'K', default: 355, min: 200, max: 600, group: 'Knock-down (empirical)' },
    { key: 't_metal_m', label: 'Metal part thickness', unit: 'm', default: 0.003, min: 1e-4, group: 'Corrosion pitting (empirical)' },
    { key: 'k_pit', label: 'Pit growth coefficient', unit: 'm/yrⁿ', default: 1.2e-4, min: 0, group: 'Corrosion pitting (empirical)', help: 'Depth = k·tⁿ. Entirely environment- and alloy-specific: fit to inspection findings' },
    { key: 'n_pit', label: 'Pit growth exponent', unit: '-', default: 0.333, min: 0.1, max: 1, group: 'Corrosion pitting (empirical)', help: '1/3 corresponds to constant volumetric dissolution of a hemispherical pit' },
    { key: 'protect_yr', label: 'Protection life before pitting starts', unit: 'yr', default: 8, min: 0, group: 'Corrosion pitting (empirical)', help: 'Life of anodising, primer and sealant' },
    { key: 'nz', label: 'Grid intervals through thickness', unit: '', default: 60, min: 8, max: 600, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c, up) { return { h_m: up.composites?.laminate_t_m, t_metal_m: up.fea?.t_skin_root_m ?? Math.min(c.struct.t_skin_mm, 3) / 1e3, years: c.econ.life_yr, rh: N.clamp(c.site.rh + 0.25, 0.3, 0.95) }; },
  run(i) {
    const tEnd = i.years * 3.15576e7, cs = i.M_sat100 * i.rh, two = i.sides === 'both', nz = Math.max(8, Math.round(i.nz)), nt = Math.max(60, 4 * nz), warnings = [];
    const f = fick(i.h_m, i.D_m2s, i.M0, cs, tEnd, nz, nt, two), Mend = f.M[nt], heff = two ? i.h_m : 2 * i.h_m, t90 = (0.0 + 0.305 * heff * heff) / i.D_m2s; // 90% saturation: D·t/h² ≈ 0.305 for two-sided exposure (from the series solution)
    const Tg = i.Tg_dry_K - i.dTg_per_pct * Mend, margin = Tg - i.T_service_K, ret = Math.max(0, 1 - (i.k_str_per_pct * Mend) / 100);
    if (margin < 28) warnings.push(`Wet glass transition is only ${margin.toFixed(0)} K above the service temperature; a margin of about 28 K (50 °F) is common design practice.`);
    // pitting
    const tp = Math.max(0, i.years - i.protect_yr), dpit = i.k_pit * tp ** i.n_pit, rem = Math.max(0, 1 - dpit / i.t_metal_m), tThru = i.k_pit > 0 ? i.protect_yr + (i.t_metal_m / i.k_pit) ** (1 / i.n_pit) : Infinity, t10 = i.k_pit > 0 ? i.protect_yr + ((0.1 * i.t_metal_m) / i.k_pit) ** (1 / i.n_pit) : Infinity;
    if (rem < 0.9) warnings.push('Pit depth exceeds 10% of the thickness: typical blend-out limits are exceeded and the pit acts as a fatigue crack starter.');
    const yr = f.t.map((t) => t / 3.15576e7), ys = N.linspace(0, i.years, 60);
    return {
      kpis: [
        kpi('moisture_pct', 'Average moisture content at end of exposure', Mend, '%'),
        kpi('saturation_frac', 'Fraction of equilibrium uptake reached', cs > i.M0 ? (Mend - i.M0) / (cs - i.M0) : 1, '-'),
        kpi('M_equilibrium_pct', 'Equilibrium moisture content', cs, '%'),
        kpi('t_90_yr', 'Time to 90% of equilibrium', t90 / 3.15576e7, 'yr'),
        kpi('Tg_wet_K', 'Wet glass-transition temperature', Tg, 'K', margin >= 28 ? 'ok' : margin > 0 ? 'warn' : 'bad', 'Empirical linear knock-down'),
        kpi('Tg_margin_K', 'Margin of wet Tg over service temperature', margin, 'K'),
        kpi('strength_retention', 'Matrix-dominated strength retention', ret, '-', ret > 0.85 ? 'ok' : 'warn', 'Empirical placeholder: replace with hot-wet coupon data'),
        kpi('pit_depth_m', 'Corrosion pit depth at end of exposure', dpit, 'm'),
        kpi('section_remaining', 'Remaining thickness fraction at the pit', rem, '-', rem >= 0.9 ? 'ok' : rem > 0.6 ? 'warn' : 'bad'),
        kpi('t_pit_10pct_yr', 'Time for a pit to reach 10% of thickness', t10, 'yr'),
        kpi('t_perforation_yr', 'Time to perforation', tThru, 'yr'),
      ],
      plots: [
        { type: 'line', title: 'Moisture uptake', xlabel: 'Time [yr]', ylabel: 'Average moisture content [%]', series: [{ name: 'Finite-difference solution', x: thinArr(yr), y: thinArr(f.M) }, ...(two ? [{ name: 'Series solution', x: ys, y: ys.map((y) => fickExact(i.h_m, i.D_m2s, i.M0, cs, y * 3.15576e7)), style: 'dash' }] : [])], annotations: [{ y: cs, label: 'Equilibrium' }] },
        { type: 'line', title: 'Moisture profile through the thickness', xlabel: 'z [mm]', ylabel: 'Moisture content [%]', series: f.prof.map((p) => ({ name: `${(p.t / 3.15576e7).toPrecision(2)} yr`, x: f.z.map((z) => z * 1e3), y: p.c })) },
        { type: 'line', title: 'Corrosion pit growth (empirical power law)', xlabel: 'Time [yr]', ylabel: 'Pit depth [mm]', series: [{ name: 'Pit depth', x: ys, y: ys.map((y) => i.k_pit * Math.max(0, y - i.protect_yr) ** i.n_pit * 1e3) }], annotations: [{ y: i.t_metal_m * 100, label: '10% of thickness' }] },
      ],
      outputs: { equivalent_flaw_m: dpit },
      warnings,
      models: ['1-D Fickian diffusion, Crank–Nicolson finite differences', 'Linear moisture knock-down of Tg and matrix-dominated strength (empirical)', 'Power-law pit growth after a protection life (empirical)'],
      assumptions: ['Constant humidity and temperature (use a time-averaged environment); no moisture cycling, no edge diffusion', 'Constant diffusivity; non-Fickian (two-stage) uptake is not modelled', 'Knock-down and pitting parameters are user-supplied empirical values: the defaults are illustrative only', 'Pit treated as an equivalent initial flaw for Suite 9; galvanic, exfoliation and stress-corrosion mechanisms are not modelled'],
    };
  },
  convergence: { param: 'nz', label: 'Grid intervals through thickness', levels: [10, 20, 40, 80, 160], metric: 'moisture_pct' },
  verify() {
    const h = 0.002, D = 1e-12, t = 2e5, f = fick(h, D, 0, 1, t, 100, 800, true), f1 = fick(h / 2, D, 0, 1, t, 50, 800, false), sm = (0.01 * h * h) / D, fs = fick(h, D, 0, 1, sm, 400, 3200, true);
    return [
      N.check('Two-sided uptake versus the Fourier series solution', f.M[800], fickExact(h, D, 0, 1, t), 2e-4, 'Crank, The Mathematics of Diffusion, eq. 4.18'),
      N.check('One-sided plate of half thickness equals the two-sided plate', f1.M[800], f.M[800], 1e-6, 'Symmetry'),
      N.check('Short-time uptake M/M∞ = 4·√(D·t/(π·h²))', fs.M[3200], 4 * Math.sqrt(0.01 / Math.PI), 5e-3, 'Semi-infinite solution (Shen & Springer, 1976)'),
    ];
  },
  calibration: { params: [{ key: 'D_m2s', min: 1e-16, max: 1e-9 }, { key: 'M_sat100', min: 0, max: 8 }], sweep: 'years', target: 'moisture_pct', note: 'Coupon mass gain versus conditioning time calibrates diffusivity and saturation content (ASTM D5229-type test).' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.Tg_margin_K < 28) out.push({ severity: o.Tg_margin_K < 0 ? 'critical' : 'warn', title: 'Insufficient hot-wet margin', detail: `Wet Tg ${o.Tg_wet_K.toFixed(0)} K against ${i.T_service_K} K service temperature.`, action: 'Select a higher-Tg resin, limit the service temperature (paint colour, heat-source shielding), or seal the laminate against moisture.', basis: 'Material operational limit: wet Tg minus a margin (programme-specific; about 28 K is common)' });
    if (o.section_remaining < 0.9) out.push({ severity: 'warn', title: 'Corrosion pitting consumes the section', detail: `Pit depth ${(o.pit_depth_m * 1e3).toFixed(2)} mm after ${i.years} years.`, action: `Inspect before year ${Number.isFinite(o.t_pit_10pct_yr) ? o.t_pit_10pct_yr.toFixed(0) : '-'}; renew surface protection, improve drainage and sealing, and pass the pit depth to Suite 9 as the initial flaw.`, basis: 'Empirical pit-growth law; blend-out limit of about 10% thickness' });
    out.push({ severity: 'info', title: 'Durability is a sustainability lever', detail: `Laminate reaches ${(100 * o.saturation_frac).toFixed(0)}% of equilibrium moisture in ${i.years} years.`, action: 'Design to end-of-life (saturated) properties from the start: it avoids mid-life restrictions and extends the useful life of high-embodied-energy structure.', basis: 'Environmental design condition' });
    return out;
  },
};
const thinArr = (a, max = 200) => { if (a.length <= max) return a.slice(); const o = []; for (let k = 0; k < max; k++) o.push(a[Math.round((k * (a.length - 1)) / (max - 1))]); return o; };

// ---- 7. metal yield criteria ----------------------------------------------------------------
/** Hill (1948) plane-stress equivalent stress referred to the yield stress in the x (rolling) direction, from Lankford r-values. */
const hill48 = (sx, sy, t, r0, r45, r90) => Math.sqrt(Math.max(0, sx * sx - ((2 * r0) / (1 + r0)) * sx * sy + ((r0 * (1 + r90)) / (r90 * (1 + r0))) * sy * sy + (((r0 + r90) * (2 * r45 + 1)) / (r90 * (1 + r0))) * t * t));
const yieldA = {
  id: 'yield', title: 'Metal yield criteria: von Mises, Tresca and Hill', fidelity: 'analytical',
  summary: 'Equivalent stress and reserve factor of a plane-stress state in a metal part by the von Mises and Tresca criteria and by Hill’s anisotropic criterion for rolled sheet, with the yield loci.',
  equations: ['von Mises yield criterion', 'Hill yield criterion', 'Generalized Hooke’s law'],
  inputs: [
    { key: 'sx', label: 'Stress σx (rolling direction)', unit: 'Pa', default: 200e6, group: 'Stress state', help: 'Defaults to the peak structural stress from Suite 2' },
    { key: 'sy', label: 'Stress σy', unit: 'Pa', default: 0, group: 'Stress state' },
    { key: 'txy', label: 'Shear stress τxy', unit: 'Pa', default: 0, group: 'Stress state' },
    { key: 'material', label: 'Material', type: 'select', options: MATS, default: 'Al 2024-T3', group: 'Material', help: 'Typical values, not design allowables' },
    { key: 'r0', label: 'Lankford coefficient r₀', unit: '-', default: 1, min: 0.2, max: 5, group: 'Anisotropy (Hill)', help: 'Width/thickness plastic strain ratio in the rolling direction; 1 = isotropic. Aluminium sheet is typically 0.5–0.8' },
    { key: 'r45', label: 'Lankford coefficient r₄₅', unit: '-', default: 1, min: 0.2, max: 5, group: 'Anisotropy (Hill)' },
    { key: 'r90', label: 'Lankford coefficient r₉₀', unit: '-', default: 1, min: 0.2, max: 5, group: 'Anisotropy (Hill)' },
  ],
  defaults: (c, up) => ({ material: c.struct.material, sx: up.fea?.sigma_bend_root_Pa ?? up.fea?.sigma_max_Pa, txy: up.fea?.tau_skin_max_Pa }),
  run(i) {
    const m = METALS[i.material] || METALS['Al 2024-T3'], c = 0.5 * (i.sx + i.sy), R = Math.hypot(0.5 * (i.sx - i.sy), i.txy), s1 = c + R, s2 = c - R, warnings = [];
    const vm = Math.sqrt(i.sx ** 2 - i.sx * i.sy + i.sy ** 2 + 3 * i.txy ** 2), tr = Math.max(Math.abs(s1 - s2), Math.abs(s1), Math.abs(s2)), hl = hill48(i.sx, i.sy, i.txy, i.r0, i.r45, i.r90);
    const rf = (s) => (s > 0 ? m.Sy / s : Infinity), rfMin = Math.min(rf(vm), rf(tr), rf(hl));
    if (rfMin < 1) warnings.push('The stress state is outside at least one yield surface: permanent deformation is predicted.');
    const th = N.linspace(0, 2 * Math.PI, 145), locus = (f) => { const x = [], y = []; for (const a of th) { const k = m.Sy / f(Math.cos(a), Math.sin(a)); x.push((k * Math.cos(a)) / 1e6); y.push((k * Math.sin(a)) / 1e6); } return { x, y }; };
    return {
      kpis: [
        kpi('sigma_vm_Pa', 'von Mises equivalent stress', vm, 'Pa'), kpi('sigma_tresca_Pa', 'Tresca equivalent stress (2·τmax)', tr, 'Pa'), kpi('sigma_hill_Pa', 'Hill equivalent stress', hl, 'Pa'),
        kpi('sigma_1_Pa', 'Major principal stress', s1, 'Pa'), kpi('sigma_2_Pa', 'Minor principal stress', s2, 'Pa'), kpi('tau_max_Pa', 'Maximum in-plane shear stress', R, 'Pa'),
        kpi('RF_vm', 'Reserve factor on yield (von Mises)', rf(vm), '-', rf(vm) >= 1.1 ? 'ok' : rf(vm) >= 1 ? 'warn' : 'bad'),
        kpi('RF_tresca', 'Reserve factor on yield (Tresca)', rf(tr), '-'), kpi('RF_hill', 'Reserve factor on yield (Hill)', rf(hl), '-'),
        kpi('RF_yield_min', 'Lowest yield reserve factor', rfMin, '-', rfMin >= 1.1 ? 'ok' : rfMin >= 1 ? 'warn' : 'bad'),
      ],
      plots: [{ type: 'line', title: 'Yield loci in the σx–σy plane (τxy = 0)', xlabel: 'σx [MPa]', ylabel: 'σy [MPa]', equalAspect: true, series: [{ name: 'von Mises', ...locus((a, b) => Math.sqrt(a * a - a * b + b * b)) }, { name: 'Tresca', ...locus((a, b) => Math.max(Math.abs(a - b), Math.abs(a), Math.abs(b))) }, { name: 'Hill 1948', ...locus((a, b) => hill48(a, b, 0, i.r0, i.r45, i.r90)), style: 'dash' }, { name: 'Stress state (direct components)', x: [i.sx / 1e6], y: [i.sy / 1e6], style: 'points' }] }],
      warnings,
      models: ['von Mises (distortion energy)', 'Tresca (maximum shear stress)', 'Hill 1948 quadratic anisotropic criterion in plane stress'],
      assumptions: ['Plane stress, initial yield only (no hardening)', 'Hill criterion normalised to the rolling-direction yield stress; r-values default to isotropy', DATA_NOTE],
    };
  },
  verify() {
    const b = { material: 'Al 2024-T3', r0: 1, r45: 1, r90: 1 }, Sy = METALS['Al 2024-T3'].Sy, sh = N.kv(yieldA.run({ ...b, sx: 0, sy: 0, txy: Sy / Math.sqrt(3) })), bi = N.kv(yieldA.run({ ...b, sx: 100e6, sy: 60e6, txy: 30e6 }));
    return [
      N.check('Pure shear yields at Sy/√3 (von Mises)', sh.RF_vm, 1, 1e-12, 'von Mises (1913)'),
      N.check('Pure shear: Tresca predicts yield at Sy/2', sh.RF_tresca, Math.sqrt(3) / 2, 1e-12, 'Tresca (1864)'),
      N.check('Hill with r = 1 reduces to von Mises', bi.sigma_hill_Pa, bi.sigma_vm_Pa, 1e-12, 'Hill (1948)'),
      N.check('Hill equibiaxial yield = σ0·√((1+r)/2) for planar isotropy r = 2', 1 / hill48(1, 1, 0, 2, 2, 2), Math.sqrt(1.5), 1e-12, 'Hill (1948), normal anisotropy'),
    ];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (o.RF_yield_min < 1) out.push({ severity: 'critical', title: 'Yield at the applied stress', detail: `Lowest reserve factor ${o.RF_yield_min.toFixed(2)}.`, action: 'Reduce the stress or select a higher-strength temper; limit load must not cause detrimental permanent deformation.', basis: 'No yield at limit load' });
    else if (o.RF_tresca < 1.05 && o.RF_vm >= 1.05) out.push({ severity: 'info', title: 'Criteria disagree near yield', detail: 'Tresca is up to 15% more conservative than von Mises in shear-dominated states.', action: 'Use von Mises for ductile aerospace alloys unless your design manual prescribes Tresca.', basis: 'Taylor–Quinney experiments' });
    return out;
  },
};

export default {
  id: 'composites', n: 21,
  tagline: 'How stiff and strong is the laminate, how does it fail, and how do time and environment degrade it?',
  analyses: [clt, strength, micromech, interlaminar, visco, environment, yieldA],
  consumes: [
    { from: 'fea', keys: ['cover_load_Npm', 'sigma_bend_root_Pa', 'tau_skin_max_Pa', 't_skin_root_m'], why: 'Running loads and stresses to check the laminate and metal against' },
    { from: 'performance', keys: ['n_limit'], why: 'Limit load factor when Suite 2 has not run' },
  ],
  provides: [
    { key: 'Ex_Pa', label: 'Laminate modulus Ex', unit: 'Pa' }, { key: 'Ey_Pa', label: 'Laminate modulus Ey', unit: 'Pa' }, { key: 'Gxy_Pa', label: 'Laminate shear modulus', unit: 'Pa' },
    { key: 'fpf_load_Npm', label: 'First-ply failure load', unit: 'N/m' }, { key: 'min_RF', label: 'Minimum reserve factor', unit: '-' }, { key: 'laminate_t_m', label: 'Laminate thickness', unit: 'm' },
  ],
  handoff: [
    { model: 'Cohesive-zone delamination and progressive damage finite elements', why: 'Needs 3-D/shell FE meshes with interface elements, fracture-energy regularisation and explicit or arc-length solvers; here only beam-theory energy release rates and ply discount are solved', tool: 'Nonlinear FE with cohesive elements and continuum damage (Abaqus / LS-DYNA / B2000++ class)' },
    { model: 'Multiscale micromechanics–FE coupling and RVE homogenisation', why: 'Requires microstructure meshes or generalised method of cells; closed-form rule-of-mixtures and Halpin–Tsai are used instead', tool: 'Micromechanics codes (GMC/HFGMC, FE-RVE)' },
    { model: 'Higher-order shear deformation and layer-wise theories, free-edge interlaminar stresses', why: 'Edge boundary-layer stress fields are three-dimensional; only equilibrium-derived τxz and the FSDT shear factor are computed', tool: '3-D or layer-wise FE' },
    { model: 'Full Puck action-plane (3-D) and fibre-kinking criteria (LaRC), in-situ strengths', why: 'Need fracture-plane searches and in-situ data not in the ply database; the 2-D Puck inter-fibre modes are implemented', tool: 'Composite post-processors / user-material FE' },
    { model: 'Impact damage, compression-after-impact and notched (open-hole, bearing) strength', why: 'Empirical or high-fidelity damage models calibrated per laminate family', tool: 'Explicit FE and coupon test programmes' },
    { model: 'Coupled moisture–thermal–mechanical and fatigue-damage evolution', why: 'Only uncoupled Fickian uptake with empirical knock-downs is solved', tool: 'Coupled-field FE; fatigue test programmes' },
    { model: 'Electrochemical corrosion and stress-corrosion cracking', why: 'Mechanism- and environment-specific kinetics; only an empirical pit-growth law is offered', tool: 'Corrosion testing and specialist electrochemical models' },
  ],
};
