// Suite 2 — Structural Mechanics and Finite Element Analysis.
// Beam finite elements for wings, rotor blades and booms (Euler–Bernoulli / Timoshenko, centrifugal geometric
// stiffness), 2-D plane-stress Q4/Q8 elements with a mapped mesh generator, Rayleigh–Ritz plate buckling with
// stringer-column checks, an axisymmetric shell strip for the pressurised fuselage, a 3-D direct-stiffness space
// truss, a corotational large-deflection beam and multi-cell Bredt–Batho torsion.

import * as N from '../core/numerics.js';
import { G0, RHO0 } from '../core/atmosphere.js';
import { METALS, designAllowables } from '../data/materials.js';

// ---- shared helpers -------------------------------------------------------------------------
const MATS = Object.keys(METALS);
const mat = (name) => METALS[name] || METALS['Al 2024-T3'];
const MAT = { key: 'material', label: 'Material', type: 'select', options: MATS, default: 'Al 2024-T3', group: 'Material', help: 'Elastic constants are typical handbook values. Strength margins use the MIL-HDBK-5J design allowable where the database holds a verified one, otherwise the typical strength; each result states which' };
const kpi = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const st3 = (m, warnAt = 0.15) => (m >= warnAt ? 'ok' : m >= 0 ? 'warn' : 'bad');
const thin = (a, max = 200) => { if (a.length <= max) return a.slice(); const o = []; for (let k = 0; k < max; k++) o.push(a[Math.round((k * (a.length - 1)) / (max - 1))]); return o; };
/** Assumption line naming the strength basis used for the margins of material m. */
const matNote = (m) => { const a = designAllowables(m); return a.design ? `Strength margins use ${a.basis}; elastic constants are typical handbook values. Confirm product form, grain direction and thickness against your own allowables` : `Strength margins use ${a.basis}; they are not statistically based design allowables`; };
/** Rough installed engine mass per engine [kg] from rating — an order-of-magnitude statistical estimate used only as a default. */
const engineMass = (p) => (p.type === 'turbofan' || p.type === 'turbojet' ? p.T0_N / (5 * G0) : p.type === 'electric' ? p.P0_W / 5000 : p.P0_W * (p.type === 'piston' ? 1.1e-3 : 4e-4));

/** Symmetric positive-definite banded solve. A[i][k] = K[i][i+k] (k = 0..bw); A and b are overwritten. */
function bandSolve(A, b, bw) {
  const n = b.length;
  for (let i = 0; i < n; i++) {
    const Ai = A[i], d = Ai[0], km = Math.min(bw, n - 1 - i);
    if (!(d > 0)) throw new Error('Stiffness matrix is singular or not positive definite: check supports and stability.');
    for (let k = 1; k <= km; k++) {
      const f = Ai[k] / d; if (f === 0) continue;
      const Aj = A[i + k]; for (let m = k; m <= km; m++) Aj[m - k] -= f * Ai[m];
      b[i + k] -= f * b[i];
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) { const Ai = A[i], km = Math.min(bw, n - 1 - i); let s = b[i]; for (let k = 1; k <= km; k++) s -= Ai[k] * x[i + k]; x[i] = s / Ai[0]; }
  return x;
}
/** Impose a zero displacement on DOF d of a banded system. */
function fixDof(A, b, d, bw) { for (let k = 1; k <= bw; k++) { A[d][k] = 0; if (d - k >= 0) A[d - k][k] = 0; } A[d][0] = 1; b[d] = 0; }

// ---- beam finite elements -------------------------------------------------------------------
const G3 = [0.5 - Math.sqrt(0.15), 0.5, 0.5 + Math.sqrt(0.15)], W3 = [5 / 18, 4 / 9, 5 / 18];
/**
 * Hermite beam FE with optional shear flexibility (Φ-corrected Timoshenko element), axial-tension geometric
 * stiffness Tg(e, ξ), Winkler foundation kf[e] and linearly varying distributed load q (nodal values).
 * DOF per node: [w, θ = dw/dx]. Returns displacements, nodal transverse force V and bending moment M = EI·w''.
 */
function beamSolve(x, EI, kGA, q, P, { fixed = [0, 1], springs = [], Tg = null, kf = null } = {}) {
  const n = x.length - 1, nd = 2 * (n + 1), A = N.range(nd, () => new Float64Array(4)), f = new Float64Array(nd), el = [];
  for (let e = 0; e < n; e++) {
    const L = x[e + 1] - x[e], phi = Number.isFinite(kGA[e]) ? (12 * EI[e]) / (kGA[e] * L * L) : 0, c = EI[e] / ((1 + phi) * L ** 3), L2 = L * L;
    const k = [[12 * c, 6 * L * c, -12 * c, 6 * L * c], [6 * L * c, (4 + phi) * L2 * c, -6 * L * c, (2 - phi) * L2 * c], [-12 * c, -6 * L * c, 12 * c, -6 * L * c], [6 * L * c, (2 - phi) * L2 * c, -6 * L * c, (4 + phi) * L2 * c]];
    if (Tg) for (let g = 0; g < 3; g++) {
      const z = G3[g], T = Tg(e, z), d = [(-6 * z + 6 * z * z) / L, 1 - 4 * z + 3 * z * z, (6 * z - 6 * z * z) / L, -2 * z + 3 * z * z];
      for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) k[a][b] += W3[g] * L * T * d[a] * d[b];
    }
    if (kf) { const m = (kf[e] * L) / 420, M = [[156, 22 * L, 54, -13 * L], [22 * L, 4 * L2, 13 * L, -3 * L2], [54, 13 * L, 156, -22 * L], [-13 * L, -3 * L2, -22 * L, 4 * L2]]; for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) k[a][b] += m * M[a][b]; }
    const q1 = q[e], q2 = q[e + 1], fe = [(L * (7 * q1 + 3 * q2)) / 20, (L2 * (3 * q1 + 2 * q2)) / 60, (L * (3 * q1 + 7 * q2)) / 20, (-L2 * (2 * q1 + 3 * q2)) / 60];
    el.push({ k, fe });
    for (let a = 0; a < 4; a++) { f[2 * e + a] += fe[a]; for (let b = a; b < 4; b++) A[2 * e + a][b - a] += k[a][b]; }
  }
  for (let j = 0; j <= n; j++) f[2 * j] += P[j] || 0;
  for (const [d, ks] of springs) A[d][0] += ks;
  const F = Array.from(f);
  for (const d of fixed) fixDof(A, f, d, 3);
  const u = bandSolve(A, f, 3), V = new Array(n + 1), M = new Array(n + 1);
  for (let e = 0; e < n; e++) {
    const ue = [u[2 * e], u[2 * e + 1], u[2 * e + 2], u[2 * e + 3]], r = el[e].k.map((row, a) => N.dot(row, ue) - el[e].fe[a]);
    V[e] = -r[0]; M[e] = -r[1]; if (e === n - 1) { V[n] = r[2]; M[n] = r[3]; }
  }
  let energy = 0; for (let d = 0; d < nd; d++) energy += 0.5 * F[d] * u[d];
  return { u, w: N.range(n + 1, (j) => u[2 * j]), th: N.range(n + 1, (j) => u[2 * j + 1]), V, M, energy };
}
/** Thin-walled rectangular wing-box section: skin (with smeared stringers) carries bending, webs carry shear, the closed cell carries torque. */
const boxSection = (w, h, ts, tw, ks) => { const te = ts * (1 + ks); return { I: 0.5 * w * te * h * h + (tw * h ** 3) / 6, Aenc: w * h, J: (4 * (w * h) ** 2) / ((2 * w) / ts + (2 * h) / tw), A: 2 * w * te + 2 * h * tw, As: 2 * h * tw }; };

const beam = {
  id: 'beam', title: 'Wing, rotor-blade or boom as a finite-element beam', fidelity: 'numerical',
  summary: 'Spanwise shear, bending and torque, deflection, stress, margins and structural mass of the primary lifting member, with section stiffness built from the wing-box geometry and optional fully-stressed sizing of the covers.',
  equations: ['Euler–Bernoulli beam equations', 'Timoshenko beam equations', 'Static equilibrium equations', 'Hooke’s law', 'Geometric stiffness equations', 'Finite element weak-form equilibrium equations', 'Saint-Venant torsion equations', 'Principle of virtual work'],
  inputs: [
    { key: 'kind', label: 'Member represented', type: 'select', options: ['wing', 'rotor blade', 'boom'], default: 'wing', group: 'Geometry', help: 'Only changes labels; the physics is set by the inputs below' },
    { key: 'L', label: 'Semi-span / member length', unit: 'm', default: 17, min: 0.05, group: 'Geometry', help: 'Root to tip, measured perpendicular to the aircraft centreline' },
    { key: 'c_root', label: 'Root chord', unit: 'm', default: 5.8, min: 0.005, group: 'Geometry' },
    { key: 'c_tip', label: 'Tip chord', unit: 'm', default: 1.4, min: 0.005, group: 'Geometry' },
    { key: 'tc', label: 'Thickness-to-chord ratio', unit: '-', default: 0.115, min: 0.02, max: 1, group: 'Geometry' },
    { key: 'sweep_deg', label: 'Structural-axis sweep', unit: 'deg', default: 25, min: -45, max: 60, group: 'Geometry' },
    { key: 'box_chord_frac', label: 'Box width / chord', unit: '-', default: 0.45, min: 0.05, max: 1, group: 'Section' },
    { key: 'box_height_frac', label: 'Box height / section thickness', unit: '-', default: 0.85, min: 0.3, max: 1, group: 'Section' },
    { key: 't_skin_mm', label: 'Cover skin thickness (minimum gauge when sizing)', unit: 'mm', default: 3, min: 0.1, max: 80, group: 'Section', help: 'Defaults to the case skin gauge, capped at 0.06 % of the root chord so small vehicles are not given airliner gauges' },
    { key: 't_taper', label: 'Tip / root skin thickness ratio', unit: '-', default: 1, min: 0.1, max: 1, group: 'Section', help: 'Only used when thicknesses are as specified' },
    { key: 't_web_mm', label: 'Spar web thickness', unit: 'mm', default: 5, min: 0.1, max: 80, group: 'Section' },
    { key: 'k_str', label: 'Stringer area / skin area', unit: '-', default: 0.5, min: 0, max: 2, group: 'Section', help: 'Smeared stringers add bending material but carry no shear; 0.4–0.6 for transport wings, 0 for a plain tube' },
    { key: 'sizing', label: 'Cover thickness', type: 'select', options: ['sized', 'as specified'], default: 'sized', group: 'Section', help: '"sized" thickens each station until the ultimate bending stress equals the design allowable (fully-stressed design)' },
    { key: 'f_design', label: 'Design stress / ultimate strength', unit: '-', default: 0.6, min: 0.2, max: 1, group: 'Section', help: 'Knock-down on the ultimate strength (the design allowable when the database holds one) kept for compression stability, fatigue and joints when sizing; 0.55–0.7 is typical for metallic covers' },
    MAT,
    { key: 'n_load', label: 'Limit load factor', unit: 'g', default: 2.5, min: -6, max: 12, group: 'Loads', help: 'From the V–n diagram of Suite 5 when available' },
    { key: 'sf', label: 'Ultimate safety factor', unit: '-', default: 1.5, min: 1, max: 3, group: 'Loads' },
    { key: 'L_1g_N', label: 'Distributed lift on this member at 1 g', unit: 'N', default: 382000, min: 0, group: 'Loads', help: 'Half the weight for a wing, thrust per blade for a rotor' },
    { key: 'dist', label: 'Lift distribution', type: 'select', options: ['schrenk', 'elliptic', 'uniform', 'triangular'], default: 'schrenk', group: 'Loads', help: 'Schrenk: mean of chord-proportional and elliptic; triangular (∝ radius) suits rotor blades' },
    { key: 'P_tip_1g_N', label: 'Tip point load at 1 g', unit: 'N', default: 0, group: 'Loads', help: 'Rotor thrust at the end of a boom, winglet load' },
    { key: 'm_fuel_kg', label: 'Fuel carried in this member', unit: 'kg', default: 7500, min: 0, group: 'Loads', help: 'Bending relief; distributed with chord squared' },
    { key: 'm_dist_kg', label: 'Other distributed mass', unit: 'kg', default: 0, min: 0, group: 'Loads', help: 'Non-structural mass spread along the member (blade fairing, systems)' },
    { key: 'm_eng_kg', label: 'Point mass on the member', unit: 'kg', default: 2400, min: 0, group: 'Loads', help: 'Wing-mounted engine and pylon: relieves bending in flight' },
    { key: 'eta_eng', label: 'Point-mass station (fraction of length)', unit: '-', default: 0.33, min: 0, max: 1, group: 'Loads' },
    { key: 'm_tip_kg', label: 'Tip mass', unit: 'kg', default: 0, min: 0, group: 'Loads' },
    { key: 'k_nonopt', label: 'Non-optimum mass factor', unit: '-', default: 1.8, min: 1, max: 4, group: 'Loads', help: 'Empirical: ribs, joints, leading/trailing edges and controls relative to the ideal box mass' },
    { key: 'q_Pa', label: 'Dynamic pressure for the torque case', unit: 'Pa', default: 31000, min: 0, group: 'Loads' },
    { key: 'Cm0', label: 'Section pitching-moment coefficient', unit: '-', default: -0.08, min: -0.5, max: 0.5, group: 'Loads' },
    { key: 'ea_offset', label: 'Elastic axis aft of aerodynamic centre / chord', unit: '-', default: 0.15, min: -0.3, max: 0.5, group: 'Loads' },
    { key: 'omega', label: 'Rotational speed (rotor blades)', unit: 'rad/s', default: 0, min: 0, group: 'Rotation', help: 'Non-zero adds centrifugal tension and its geometric stiffness' },
    { key: 'r0', label: 'Root radius from the rotation axis', unit: 'm', default: 0, min: 0, group: 'Rotation' },
    { key: 'root', label: 'Root attachment', type: 'select', options: ['clamped', 'pinned'], default: 'clamped', group: 'Supports', help: 'Pinned = flap hinge (needs rotation or a root spring to be stable)' },
    { key: 'k_root', label: 'Root rotational spring (pinned root)', unit: 'N·m/rad', default: 0, min: 0, group: 'Supports', help: 'Represents joint or carry-through flexibility instead of an ideal clamp' },
    { key: 'theory', label: 'Beam theory', type: 'select', options: ['Euler–Bernoulli', 'Timoshenko'], default: 'Euler–Bernoulli', group: 'Numerics' },
    { key: 'n_members', label: 'Identical members on the aircraft', unit: '', default: 2, min: 1, step: 1, discrete: true, group: 'Geometry' },
    { key: 'nElem', label: 'Beam elements', unit: '', default: 40, min: 2, max: 400, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c, up, d) {
    const s = c.struct, n = up.performance?.n_limit ?? c.aero.n_pos, base = { material: s.material, n_load: n, sf: s.sf_ultimate };
    if (c.wing.S_m2 > 0) {
      const onWing = c.prop.n_eng >= 2 && c.meta.type !== 'evtol', Vd = up.performance?.V_d_eas ?? 1.25 * (c.aero.Vmo_ms || c.flight.V_ms);
      return { ...base, kind: 'wing', L: c.wing.b_m / 2, c_root: d.c_root, c_tip: d.c_tip, tc: c.wing.tc, sweep_deg: c.wing.sweep_deg, box_chord_frac: s.box_chord_frac, box_height_frac: s.box_height_frac,
        t_skin_mm: Math.min(s.t_skin_mm, Math.max(0.4, 0.6 * d.c_root)), t_web_mm: Math.min(s.t_spar_mm, Math.max(0.5, d.c_root)), L_1g_N: d.W / 2, m_fuel_kg: 0.4 * c.mass.fuel_kg,
        m_eng_kg: onWing ? (engineMass(c.prop) * c.prop.n_eng) / 2 : 0, q_Pa: 0.5 * RHO0 * Vd * Vd, Cm0: c.aero.Cm0, n_members: 2 };
    }
    if (c.meta.type === 'helicopter' && c.rotor.R_m > 0) {
      const r = c.rotor, r0 = r.hinge_offset * r.R_m, L = r.R_m - r0, sec = boxSection(0.35 * r.chord_m, 0.85 * 0.12 * r.chord_m, s.t_skin_mm / 1e3, s.t_spar_mm / 1e3, 0);
      return { ...base, kind: 'rotor blade', L, c_root: r.chord_m, c_tip: r.chord_m, tc: 0.12, sweep_deg: 0, box_chord_frac: 0.35, box_height_frac: 0.85, t_skin_mm: s.t_skin_mm, t_web_mm: s.t_spar_mm, k_str: 0, sizing: 'as specified',
        L_1g_N: d.W / Math.max(1, r.n_blades), dist: 'triangular', m_fuel_kg: 0, m_eng_kg: 0, m_dist_kg: Math.max(0, r.blade_mass_kg - sec.A * mat(s.material).rho * L), k_nonopt: 1, q_Pa: 0, Cm0: 0, ea_offset: 0,
        omega: d.omega, r0, root: r.hinge_offset > 0 ? 'pinned' : 'clamped', n_members: r.n_blades };
    }
    const L = Math.max(0.1, 1.6 * (c.rotor.R_m || 0.5 * c.prop.prop_dia_m || 0.2)), tube = Math.max(0.012, 0.07 * L), nm = Math.max(1, c.prop.n_eng);
    return { ...base, kind: 'boom', L, c_root: tube, c_tip: tube, tc: 1, sweep_deg: 0, box_chord_frac: 1, box_height_frac: 1, t_skin_mm: Math.max(0.5, 40 * tube), t_web_mm: Math.max(0.5, 40 * tube), k_str: 0, sizing: 'as specified',
      L_1g_N: 0, dist: 'uniform', P_tip_1g_N: d.W / nm, m_fuel_kg: 0, m_eng_kg: 0, m_tip_kg: engineMass(c.prop), k_nonopt: 1.2, q_Pa: 0, Cm0: 0, ea_offset: 0, n_members: nm };
  },
  run(i) {
    const m = mat(i.material), al = designAllowables(m), n = Math.max(2, Math.round(i.nElem)), cosL = Math.cos(N.rad(i.sweep_deg)), Ls = i.L / cosL, warnings = [];
    const s = N.linspace(0, Ls, n + 1), eta = s.map((v) => v / Ls), chord = eta.map((e) => i.c_root + (i.c_tip - i.c_root) * e), cm = 0.5 * (i.c_root + i.c_tip);
    const bw = chord.map((c) => i.box_chord_frac * c), bh = chord.map((c) => i.box_height_frac * i.tc * c), tw = i.t_web_mm / 1e3, tmin = i.t_skin_mm / 1e3, sized = i.sizing === 'sized';
    const shape = eta.map((e, j) => (i.dist === 'elliptic' ? Math.sqrt(Math.max(0, 1 - e * e)) : i.dist === 'uniform' ? 1 : i.dist === 'triangular' ? (i.r0 + s[j]) / (i.r0 + Ls) : 0.5 * (chord[j] / cm + (4 / Math.PI) * Math.sqrt(Math.max(0, 1 - e * e)))));
    const lift = shape.map((v) => (i.n_load * i.L_1g_N * v) / N.trapz(s, shape)), c2 = N.trapz(s, chord.map((c) => c * c)), fuel = chord.map((c) => (i.m_fuel_kg * c * c) / c2);
    const jEng = Math.round(i.eta_eng * n), P = new Array(n + 1).fill(0);
    P[jEng] -= i.n_load * G0 * i.m_eng_kg; P[n] += i.n_load * (i.P_tip_1g_N - G0 * i.m_tip_kg);
    let pinned = i.root === 'pinned';
    if (pinned && !(i.omega > 0) && !(i.k_root > 0)) { pinned = false; warnings.push('A pinned root without rotation or a root spring is a mechanism; the root was treated as clamped.'); }
    let ts = eta.map((e) => tmin * (sized ? 1 : 1 + (i.t_taper - 1) * e)), sec, mdist, q, fe, Tn, sol;
    for (let it = 0, nIt = sized ? 8 : 1; it < nIt; it++) {
      sec = s.map((_, j) => boxSection(bw[j], bh[j], ts[j], tw, i.k_str));
      mdist = sec.map((x) => m.rho * x.A * i.k_nonopt + i.m_dist_kg / Ls);
      q = lift.map((l, j) => l - i.n_load * G0 * (mdist[j] + fuel[j]));
      // centrifugal tension T(s) = Ω² ∫ m'(r) r dr + tip mass term, integrated element by element (2-point Gauss on the remainder of each element)
      Tn = new Array(n + 1).fill(0); Tn[n] = i.omega ** 2 * i.m_tip_kg * (i.r0 + Ls);
      const dT = (e, z0) => { const L = s[e + 1] - s[e], h = 1 - z0; let a = 0; for (const g of [0.5 - 0.5 / Math.sqrt(3), 0.5 + 0.5 / Math.sqrt(3)]) { const z = z0 + g * h; a += 0.5 * h * L * (mdist[e] + (mdist[e + 1] - mdist[e]) * z) * (i.r0 + s[e] + L * z); } return i.omega ** 2 * a; };
      for (let e = n - 1; e >= 0; e--) Tn[e] = Tn[e + 1] + dT(e, 0);
      const mid = (a, e) => 0.5 * (a[e] + a[e + 1]);
      fe = N.range(n, (e) => boxSection(mid(bw, e), mid(bh, e), mid(ts, e), tw, i.k_str));
      sol = beamSolve(s, fe.map((x) => m.E * x.I), fe.map((x) => (i.theory === 'Timoshenko' ? m.G * x.As : Infinity)), q, P,
        { fixed: pinned ? [0] : [0, 1], springs: pinned && i.k_root > 0 ? [[1, i.k_root]] : [], Tg: i.omega > 0 ? (e, z) => Tn[e + 1] + dT(e, z) : null });
      if (!sized || it === nIt - 1) break;
      ts = ts.map((_, j) => Math.max(tmin, ((Math.abs(sol.M[j]) * i.sf * bh[j]) / (2 * i.f_design * al.Su) - (tw * bh[j] ** 3) / 6) / (0.5 * bw[j] * bh[j] ** 2 * (1 + i.k_str))));
    }
    // torque about the elastic axis and twist (Bredt–Batho single cell)
    const tq = lift.map((l, j) => l * i.ea_offset * chord[j] + i.Cm0 * i.q_Pa * chord[j] ** 2 * cosL), Tq = N.cumtrapz(s, tq).map((v, _, a) => a[n] - v);
    const twist = N.cumtrapz(s, Tq.map((t, j) => t / (m.G * sec[j].J)));
    const sb = sol.M.map((M, j) => (Math.abs(M) * bh[j]) / (2 * sec[j].I)), sax = Tn.map((T, j) => T / sec[j].A);
    // sol.V is the total transverse force on the cut; under centrifugal tension part of it, T·w′, is carried by the inclined axial force, not by the webs
    const Qs = sol.V.map((V, j) => V - Tn[j] * sol.th[j]);
    const tauS = Tq.map((t, j) => Math.abs(t) / (2 * sec[j].Aenc * ts[j])), tauW = Qs.map((Q, j) => Math.abs(Q) / sec[j].As + Math.abs(Tq[j]) / (2 * sec[j].Aenc * tw));
    const vm = sb.map((b, j) => Math.max(Math.hypot(b + sax[j], Math.sqrt(3) * tauS[j]), Math.sqrt(3) * tauW[j])), jm = N.argmax(vm), sMax = vm[jm];
    const mosU = sMax > 0 ? al.Su / (i.sf * sMax) - 1 : Infinity, mosY = sMax > 0 ? al.Sy / sMax - 1 : Infinity, mos = Math.min(mosU, mosY);
    const mBox = N.trapz(s, sec.map((x) => m.rho * x.A)), mOne = mBox * i.k_nonopt, mAll = mOne * Math.round(i.n_members), tip = sol.w[n], y = s.map((v) => v * cosL);
    const gauge = sized ? ts.filter((t) => t <= tmin * 1.0001).length / (n + 1) : NaN;
    if (Math.abs(tip) / Ls > 0.15 && !pinned) warnings.push(`Tip deflection is ${(100 * Math.abs(tip) / Ls).toFixed(0)}% of the length: small-deflection theory overstates it. Use the large-deflection analysis.`);
    if (bh[0] / Ls > 0.1 && i.theory !== 'Timoshenko') warnings.push('The member is deep relative to its length; switch to Timoshenko theory to include shear deflection.');
    if (Math.abs(i.sweep_deg) > 35) warnings.push('High sweep: bending–torsion coupling of a swept box and the root triangle are not represented by a straight beam.');
    if (mos < 0) warnings.push(`Negative margin of safety at ${(100 * eta[jm]).toFixed(0)}% of the length: the section yields at limit load or fails before ultimate load.`);
    if (i.omega > 0 && !pinned) warnings.push('Clamped rotating blade: a narrow bending boundary layer forms at the root; refine the mesh there or check the hub fitting separately.');
    const k = thin(N.range(n + 1), 161), pick = (a, f = 1) => k.map((j) => a[j] * f), yk = pick(y), wing = i.kind === 'wing';
    const rows = thin(N.range(n + 1), 11).map((j) => [y[j], chord[j], ts[j] * 1e3, (m.E * sec[j].I) / 1e6, (m.G * sec[j].J) / 1e6, sol.V[j] / 1e3, sol.M[j] / 1e3, Tq[j] / 1e3, vm[j] / 1e6]);
    return {
      kpis: [
        kpi('sigma_max_Pa', 'Peak von Mises stress at limit load', sMax, 'Pa', st3(mos), `Yield ${(al.Sy / 1e6).toFixed(0)} MPa, ultimate ${(al.Su / 1e6).toFixed(0)} MPa (${al.design ? 'design allowables' : 'typical strengths'})`),
        kpi('margin_of_safety', 'Governing margin of safety', mos, '-', st3(mos), `${mosU < mosY ? 'Ultimate strength governs' : 'Yield at limit load governs'}; ${al.design ? 'design allowables' : 'typical strengths, not design allowables'}`),
        kpi('tip_deflection_m', 'Tip deflection at limit load', tip, 'm', Math.abs(tip) / Ls > 0.15 && !pinned ? 'warn' : 'ok', pinned ? `Rigid-body coning of ${N.deg(Math.atan2(tip, Ls)).toFixed(1)}° about the hinge plus elastic bending` : ''),
        kpi('tip_twist_deg', 'Tip twist', N.deg(twist[n]), 'deg'),
        kpi('root_moment_Nm', 'Root bending moment', sol.M[0], 'N·m'),
        kpi('root_shear_N', 'Root shear force', sol.V[0], 'N', undefined, i.omega > 0 ? 'Total transverse force (root reaction); the web shear excludes the part carried by the inclined centrifugal tension' : ''),
        kpi('root_torque_Nm', 'Root torque', Tq[0], 'N·m'),
        kpi('sigma_bend_root_Pa', 'Root cover bending stress', sb[0] + sax[0], 'Pa'),
        kpi('tau_skin_max_Pa', 'Peak skin shear stress', N.amax(tauS), 'Pa'),
        kpi('EI_root_Nm2', 'Root bending stiffness EI', m.E * sec[0].I, 'N·m²'),
        kpi('GJ_root_Nm2', 'Root torsional stiffness GJ', m.G * sec[0].J, 'N·m²'),
        kpi('t_skin_root_m', 'Root cover skin thickness', ts[0], 'm'),
        kpi('cover_load_Npm', 'Root cover running load at limit load', (sb[0] + sax[0]) * ts[0] * (1 + i.k_str), 'N/m', undefined, 'Direct stress × smeared cover thickness: the sizing load for a composite cover in Suite 21'),
        kpi('member_mass_kg', `Mass of one ${i.kind}`, mOne, 'kg', undefined, 'Ideal box mass × non-optimum factor (empirical)'),
        kpi('wing_struct_mass_kg', 'Wing structural mass (all wings)', wing ? mAll : 0, 'kg', undefined, wing ? '' : 'No wing on this member'),
        kpi('centrifugal_root_N', 'Root centrifugal tension', Tn[0], 'N'),
        kpi('strain_energy_J', 'Strain energy (½ F·u)', sol.energy, 'J'),
      ],
      plots: [
        { type: 'line', title: 'Applied loading at limit load', xlabel: 'Spanwise position [m]', ylabel: 'Running load [kN/m]', series: [{ name: 'Lift', x: yk, y: pick(lift, 1e-3) }, { name: 'Inertia relief', x: yk, y: k.map((j) => (-i.n_load * G0 * (mdist[j] + fuel[j])) / 1e3) }, { name: 'Net', x: yk, y: pick(q, 1e-3), style: 'dash' }] },
        { type: 'line', title: 'Shear force diagram', xlabel: 'Spanwise position [m]', ylabel: 'Shear force [kN]', series: [{ name: 'Shear', x: yk, y: pick(sol.V, 1e-3) }] },
        { type: 'line', title: 'Bending moment and torque diagrams', xlabel: 'Spanwise position [m]', ylabel: 'Moment [kN·m]', series: [{ name: 'Bending moment', x: yk, y: pick(sol.M, 1e-3) }, { name: 'Torque', x: yk, y: pick(Tq, 1e-3) }] },
        { type: 'line', title: 'Deflected shape at limit load', xlabel: 'Spanwise position [m]', ylabel: 'Deflection [m]', series: [{ name: 'Deflection', x: yk, y: pick(sol.w) }] },
        { type: 'line', title: 'Stress along the member at limit load', xlabel: 'Spanwise position [m]', ylabel: 'Stress [MPa]', series: [{ name: 'Cover direct stress', x: yk, y: k.map((j) => (sb[j] + sax[j]) / 1e6) }, { name: 'Web shear', x: yk, y: pick(tauW, 1e-6) }, { name: 'Skin shear (torque)', x: yk, y: pick(tauS, 1e-6) }, { name: 'von Mises', x: yk, y: pick(vm, 1e-6), style: 'dash' }], annotations: [{ y: al.Sy / 1e6, label: al.design ? 'Yield allowable' : 'Yield (typical)' }] },
        { type: 'line', title: 'Cover skin thickness', xlabel: 'Spanwise position [m]', ylabel: 'Thickness [mm]', series: [{ name: 'Skin', x: yk, y: pick(ts, 1e3), style: sized ? 'line' : 'dash' }] },
      ],
      tables: [{ title: 'Station summary', columns: ['y [m]', 'Chord [m]', 'Skin [mm]', 'EI [MN·m²]', 'GJ [MN·m²]', 'Shear [kN]', 'Moment [kN·m]', 'Torque [kN·m]', 'von Mises [MPa]'], rows }],
      outputs: { span_y_m: thin(y, 41), EI_Nm2: thin(sec.map((x) => m.E * x.I), 41), GJ_Nm2: thin(sec.map((x) => m.G * x.J), 41), mass_per_m_kgm: thin(mdist.map((v, j) => (v + fuel[j]) / cosL), 41), min_gauge_fraction: gauge, web_shear_root_N: Qs[0], allowable_basis: al.design ? 1 : 0 },
      warnings,
      models: [`${i.theory} Hermite beam elements (${n})`, sized ? 'Fully-stressed cover sizing with minimum gauge' : 'Thicknesses as specified', i.dist === 'schrenk' ? 'Schrenk spanwise lift approximation' : `${i.dist} lift distribution`, 'Bredt–Batho single-cell torsion', ...(i.omega > 0 ? ['Centrifugal tension with consistent geometric stiffness'] : [])],
      assumptions: ['Straight beam along the swept structural axis; no bending–torsion coupling', 'Thin-walled rectangular box: covers and smeared stringers carry bending, webs carry shear', 'Structural mass = ideal box mass × empirical non-optimum factor; relief from structure, fuel and point masses scales with load factor', 'Default non-optimum factor, design-stress ratio, stringer ratio, fuel share and engine mass are typical or empirical values, not sourced data', 'Stresses are nominal: no joints, cut-outs, buckling or fatigue knock-downs', matNote(m)],
    };
  },
  convergence: { param: 'nElem', label: 'Beam elements', levels: [5, 10, 20, 40, 80], metric: 'tip_deflection_m' },
  verify() {
    const b = { kind: 'boom', L: 2, c_root: 0.1, c_tip: 0.1, tc: 1, sweep_deg: 0, box_chord_frac: 1, box_height_frac: 0.5, t_skin_mm: 2, t_taper: 1, t_web_mm: 2, k_str: 0, sizing: 'as specified', f_design: 0.75, material: 'Al 2024-T3', n_load: 1, sf: 1.5, L_1g_N: 0, dist: 'uniform', P_tip_1g_N: 100, m_fuel_kg: 0, m_dist_kg: 0, m_eng_kg: 0, eta_eng: 0.3, m_tip_kg: 0, k_nonopt: 1, q_Pa: 0, Cm0: 0, ea_offset: 0, omega: 0, r0: 0, root: 'clamped', k_root: 0, theory: 'Euler–Bernoulli', n_members: 1, nElem: 8 };
    const mm = mat(b.material), sec = boxSection(0.1, 0.05, 0.002, 0.002, 0), EI = mm.E * sec.I, mp = mm.rho * sec.A, P = 100, L = 2, wq = -G0 * mp;
    const r1 = N.kv(beam.run(b)), r2 = N.kv(beam.run({ ...b, theory: 'Timoshenko' })), r3 = N.kv(beam.run({ ...b, L_1g_N: 500, P_tip_1g_N: 0, material: 'Al 2024-T3' }));
    // rotating pinned blade under a load proportional to radius: exact solution is a straight line w = k·r/(m'Ω²) with zero bending
    // (the difference of two lift levels removes the self-weight contribution, which has no straight-line solution)
    const om = 60, kq = (2 * 500) / (L * L), rotRun = (Lf) => beam.run({ ...b, L_1g_N: Lf, P_tip_1g_N: 0, dist: 'triangular', omega: om, root: 'pinned' }), rot = (Lf) => N.kv(rotRun(Lf)).tip_deflection_m;
    const q1 = rotRun(1000), q5 = rotRun(500), dV = N.kv(q1).root_shear_N - N.kv(q5).root_shear_N, dQ = q1.outputs.web_shear_root_N - q5.outputs.web_shear_root_N;
    return [
      N.check('Cantilever tip load: PL³/3EI + self-weight qL⁴/8EI', r1.tip_deflection_m, (P * L ** 3) / (3 * EI) + (wq * L ** 4) / (8 * EI), 1e-9, 'Euler–Bernoulli beam theory (Hermite elements are nodally exact)'),
      N.check('Timoshenko shear deflection adds PL/kGA + qL²/2kGA', r2.tip_deflection_m - r1.tip_deflection_m, (P * L) / (mm.G * sec.As) + (wq * L * L) / (2 * mm.G * sec.As), 1e-6, 'Timoshenko beam theory'),
      N.check('Uniform load root moment qL²/2', r3.root_moment_Nm, ((500 / L + wq) * L * L) / 2, 1e-9, 'Statics'),
      N.check('Strain energy equals ∫M²/2EI', r1.strain_energy_J, N.simpson((x) => (P * (L - x) + (wq * (L - x) ** 2) / 2) ** 2 / (2 * EI), 0, L, 200), 1e-6, 'Clapeyron theorem (cubic interpolation of the quartic self-weight deflection gives O(h⁴) work error)'),
      N.check('Rotating hinged blade cones as a straight line: w_tip = k·R/(m′Ω²)', rot(1000) - rot(500), (kq * L) / (mp * om * om), 1e-8, 'String equation (T w′)′ + f = 0 with T = ½m′Ω²(R² − r²) and f = k·r'),
      N.check('Rotating hinged blade: the root reaction equals the added lift', dV, 500, 1e-8, 'Vertical equilibrium of the straight-line coning solution'),
      N.check('Rotating hinged blade: no web shear in pure coning (reaction carried by the inclined tension)', dQ / dV, 0, 1e-8, 'Q = −(EI·w″)′ = 0 for a straight line; V = Q + T·w′'),
      N.check('Margins use the A-basis allowable when the database holds one', (N.kv(beam.run({ ...b, material: 'Al 2024-T3' })).margin_of_safety + 1) * N.kv(beam.run({ ...b, material: 'Al 2024-T3' })).sigma_max_Pa, Math.min(441e6 / 1.5, 324e6), 1e-12, 'MIL-HDBK-5J 2024-T3 sheet A-basis Ftu 64 ksi, Fty 47 ksi'),
      N.check('Bredt–Batho GJ of the rectangular cell', r1.GJ_root_Nm2, (mm.G * 4 * (0.1 * 0.05) ** 2 * 0.002) / (2 * 0.1 + 2 * 0.05), 1e-12, 'Bredt–Batho formula'),
    ];
  },
  recommend(res, i, ctx) {
    const o = res.outputs, out = [], fuel = ctx.case.mass.fuel_kg, mtow = ctx.case.mass.mtow_kg;
    if (o.margin_of_safety < 0) out.push({ severity: 'critical', title: 'Primary structure fails the static strength check', detail: `Margin of safety ${o.margin_of_safety.toFixed(2)} with ${(o.sigma_max_Pa / 1e6).toFixed(0)} MPa at limit load.`, action: 'Increase cover or web thickness (or select "sized"), deepen the box, reduce the limit load factor through manoeuvre/gust load alleviation, or choose a stronger alloy.', basis: 'No yield at limit load; no failure at limit × ultimate factor' });
    else if (o.margin_of_safety > 1 && i.sizing !== 'sized') out.push({ severity: 'advise', title: 'Structure is lightly stressed', detail: `Margin of safety ${o.margin_of_safety.toFixed(2)}: static strength does not size this member.`, action: 'Check stiffness, buckling and fatigue before removing material; if none governs, thinner gauges save mass.', basis: 'Fully-stressed design principle' });
    if (i.kind === 'wing' && o.min_gauge_fraction > 0.5) out.push({ severity: 'info', title: 'Minimum gauge governs most of the span', detail: `${(100 * o.min_gauge_fraction).toFixed(0)}% of stations sit at the minimum skin gauge.`, action: 'Outboard wing mass is set by handling, buckling and manufacturing limits rather than strength: consider a lower gauge with closer stiffening, or composite covers (Suite 21).', basis: 'Sizing result' });
    if (i.kind === 'wing' && mtow > 0 && o.wing_struct_mass_kg > 0) {
      const dm = 0.05 * o.wing_struct_mass_kg, df = fuel > 0 ? (fuel * dm) / mtow : 0;
      out.push({ severity: 'info', title: 'Mass lever', detail: `Wing structure ≈ ${o.wing_struct_mass_kg.toFixed(0)} kg (${(100 * o.wing_struct_mass_kg / mtow).toFixed(1)}% of MTOM). A 5% saving (${dm.toFixed(0)} kg) is worth about ${df.toFixed(1)} kg of fuel and ${(3.16 * df).toFixed(1)} kg CO₂ per design-range flight${fuel > 0 ? '' : ' (electric: energy saving scales the same way)'}.`, action: 'Trade load alleviation, higher design stress with verified fatigue life, or a deeper box against drag in Suite 23.', basis: 'First-order Breguet sensitivity ΔW_fuel/W_fuel ≈ ΔW/W; 3.16 kg CO₂ per kg kerosene' });
    }
    if (Math.abs(o.tip_twist_deg) > 3) out.push({ severity: 'advise', title: 'Significant elastic twist', detail: `Tip twist ${o.tip_twist_deg.toFixed(1)}° changes the spanwise lift distribution.`, action: 'Pass GJ to Suite 3 for divergence, aileron reversal and flutter checks.', basis: 'Static aeroelastic coupling' });
    return out;
  },
};

// ---- 2-D plane elasticity -------------------------------------------------------------------
const QX = [-1, 1, 1, -1, 0, 1, 0, -1], QY = [-1, -1, 1, 1, -1, 0, 1, 0], R3 = 1 / Math.sqrt(3);
function dShape(nn, xi, et) {
  const dx = new Array(nn), dy = new Array(nn);
  for (let a = 0; a < nn; a++) {
    const xa = QX[a], ya = QY[a];
    if (nn === 4) { dx[a] = 0.25 * xa * (1 + et * ya); dy[a] = 0.25 * ya * (1 + xi * xa); }
    else if (a < 4) { dx[a] = 0.25 * xa * (1 + et * ya) * (2 * xi * xa + et * ya); dy[a] = 0.25 * ya * (1 + xi * xa) * (xi * xa + 2 * et * ya); }
    else if (xa === 0) { dx[a] = -xi * (1 + et * ya); dy[a] = 0.5 * ya * (1 - xi * xi); }
    else { dx[a] = 0.5 * xa * (1 - et * et); dy[a] = -et * (1 + xi * xa); }
  }
  return [dx, dy];
}
/** Cartesian shape-function gradients and Jacobian determinant at (ξ, η). */
function gradQ(xy, xi, et) {
  const nn = xy.length, [dx, dy] = dShape(nn, xi, et); let xx = 0, yx = 0, xe = 0, ye = 0;
  for (let a = 0; a < nn; a++) { xx += dx[a] * xy[a][0]; yx += dx[a] * xy[a][1]; xe += dy[a] * xy[a][0]; ye += dy[a] * xy[a][1]; }
  const det = xx * ye - yx * xe;
  return { det, bx: dx.map((v, a) => (ye * v - yx * dy[a]) / det), by: dx.map((v, a) => (-xe * v + xx * dy[a]) / det) };
}
/**
 * Linear plane-stress / plane-strain FE with 4-node (2×2 Gauss) or 8-node serendipity (3×3, or 2×2 reduced) quadrilaterals.
 * Stresses are sampled at the 2×2 Gauss points, extrapolated to the corners and averaged at the nodes.
 */
function planeFE({ nodes, elems, E, nu, t, planeStrain = false, reduced = false, fixed, f }) {
  const nd = 2 * nodes.length, nn = elems[0].length, c = planeStrain ? E / ((1 + nu) * (1 - 2 * nu)) : E / (1 - nu * nu);
  const D11 = planeStrain ? c * (1 - nu) : c, D12 = c * nu, D33 = planeStrain ? c * (1 - 2 * nu) / 2 : c * (1 - nu) / 2;
  const gp = nn === 8 && !reduced ? [[-Math.sqrt(0.6), 5 / 9], [0, 8 / 9], [Math.sqrt(0.6), 5 / 9]] : [[-R3, 1], [R3, 1]];
  let bw = 0; for (const el of elems) { const lo = N.amin(el), hi = N.amax(el); bw = Math.max(bw, 2 * (hi - lo) + 1); }
  const A = N.range(nd, () => new Float64Array(bw + 1)), b = Float64Array.from(f), kes = [];
  for (const el of elems) {
    const xy = el.map((k) => nodes[k]), ke = N.zeros(2 * nn);
    for (const [xi, wx] of gp) for (const [et, wy] of gp) {
      const g = gradQ(xy, xi, et), w = wx * wy * g.det * t;
      if (!(g.det > 0)) throw new Error('Inverted or degenerate element in the mesh.');
      for (let a = 0; a < nn; a++) for (let q = 0; q < nn; q++) {
        const ax = g.bx[a], ay = g.by[a], qx = g.bx[q], qy = g.by[q];
        ke[2 * a][2 * q] += w * (ax * D11 * qx + ay * D33 * qy); ke[2 * a][2 * q + 1] += w * (ax * D12 * qy + ay * D33 * qx);
        ke[2 * a + 1][2 * q] += w * (ay * D12 * qx + ax * D33 * qy); ke[2 * a + 1][2 * q + 1] += w * (ay * D11 * qy + ax * D33 * qx);
      }
    }
    kes.push(ke);
    for (let a = 0; a < 2 * nn; a++) { const ga = 2 * el[a >> 1] + (a & 1); for (let q = 0; q < 2 * nn; q++) { const gq = 2 * el[q >> 1] + (q & 1); if (gq >= ga) A[ga][gq - ga] += ke[a][q]; } }
  }
  for (const d of fixed) fixDof(A, b, d, bw);
  const u = bandSolve(A, b, bw), resid = f.map((v) => -v), acc = N.range(nodes.length, () => [0, 0, 0, 0]);
  elems.forEach((el, e) => {
    const xy = el.map((k) => nodes[k]), ue = []; for (const k of el) ue.push(u[2 * k], u[2 * k + 1]);
    kes[e].forEach((row, a) => { resid[2 * el[a >> 1] + (a & 1)] += N.dot(row, ue); });
    const sg = [];
    for (const et of [-R3, R3]) for (const xi of [-R3, R3]) {
      const g = gradQ(xy, xi, et); let ex = 0, ey = 0, gxy = 0;
      for (let a = 0; a < nn; a++) { ex += g.bx[a] * ue[2 * a]; ey += g.by[a] * ue[2 * a + 1]; gxy += g.by[a] * ue[2 * a] + g.bx[a] * ue[2 * a + 1]; }
      sg.push([D11 * ex + D12 * ey, D12 * ex + D11 * ey, D33 * gxy]);
    }
    const corner = [0, 1, 2, 3].map((a) => { const X = QX[a] / R3, Y = QY[a] / R3, w = [(1 - X) * (1 - Y) / 4, (1 + X) * (1 - Y) / 4, (1 - X) * (1 + Y) / 4, (1 + X) * (1 + Y) / 4]; return [0, 1, 2].map((s) => w[0] * sg[0][s] + w[1] * sg[1][s] + w[2] * sg[2][s] + w[3] * sg[3][s]); });
    for (let a = 0; a < nn; a++) {
      const v = a < 4 ? corner[a] : [0, 1, 2].map((s) => 0.5 * (corner[a - 4][s] + corner[(a - 3) % 4][s])), r = acc[el[a]];
      r[0] += v[0]; r[1] += v[1]; r[2] += v[2]; r[3]++;
    }
  });
  const sx = acc.map((r) => r[0] / r[3]), sy = acc.map((r) => r[1] / r[3]), txy = acc.map((r) => r[2] / r[3]);
  const sz = planeStrain ? sx.map((v, k) => nu * (v + sy[k])) : sx.map(() => 0);
  const vm = sx.map((v, k) => Math.sqrt(0.5 * ((v - sy[k]) ** 2 + (sy[k] - sz[k]) ** 2 + (sz[k] - v) ** 2) + 3 * txy[k] ** 2));
  return { u, sx, sy, txy, vm, resid, energy: 0.5 * N.dot(f, u), ndof: nd, bw };
}
/** Structured mapped mesh: P(ρ, τ) → [x, y] for ρ, τ in [0, 1]; ρ maps to the element ξ direction. */
function gridMesh(nI, nJ, P, q8) {
  const m = q8 ? 2 : 1, NI = m * nI, NJ = m * nJ, id = N.range(NJ + 1, () => new Array(NI + 1).fill(-1)), nodes = [], elems = [];
  for (let j = 0; j <= NJ; j++) for (let i = 0; i <= NI; i++) { if (q8 && i % 2 && j % 2) continue; id[j][i] = nodes.length; nodes.push(P(i / NI, j / NJ)); }
  for (let j = 0; j < nJ; j++) for (let i = 0; i < nI; i++) {
    const a = m * i, b = m * j;
    elems.push(q8 ? [id[b][a], id[b][a + 2], id[b + 2][a + 2], id[b + 2][a], id[b][a + 1], id[b + 1][a + 2], id[b + 2][a + 1], id[b + 1][a]] : [id[b][a], id[b][a + 1], id[b + 1][a + 1], id[b + 1][a]]);
  }
  return { nodes, elems, id, NI, NJ };
}
/** Consistent nodal forces for a uniform traction (tx, ty) on a straight boundary run of node ids. */
function edgeLoad(f, nodes, ids, q8, tx, ty, t) {
  for (let st = q8 ? 2 : 1, k = 0; k + st < ids.length; k += st) {
    const a = ids[k], c = ids[k + st], len = Math.hypot(nodes[c][0] - nodes[a][0], nodes[c][1] - nodes[a][1]) * t;
    const parts = q8 ? [[a, 1 / 6], [ids[k + 1], 2 / 3], [c, 1 / 6]] : [[a, 0.5], [c, 0.5]];
    for (const [nd, w] of parts) { f[2 * nd] += w * len * tx; f[2 * nd + 1] += w * len * ty; }
  }
}
const quadTris = (elems) => elems.flatMap((e) => (e.length === 4 ? [[e[0], e[1], e[2]], [e[0], e[2], e[3]]] : [[e[0], e[4], e[7]], [e[4], e[1], e[5]], [e[5], e[2], e[6]], [e[6], e[3], e[7]], [e[4], e[5], e[7]], [e[5], e[6], e[7]]]));
/** Quarter model of a rectangular plate (half-length Lx along the load, half-width Ly) with a central elliptical hole (semi-axes a along x, b along y). */
function holeModel(i) {
  const q8 = i.elem !== 'Q4', nI = Math.max(2, Math.round(i.nMesh)), nJ = 2 * Math.ceil(nI / 2), Lx = i.L_m / 2, Ly = i.W_m / 2, b = i.d_m / 2, a = b * i.ellip, g = Math.max(1.0001, i.grade);
  const P = (r, t) => {
    const ph = (t * Math.PI) / 2, xi = a * Math.cos(ph), yi = b * Math.sin(ph), xo = t <= 0.5 ? Lx : Lx * (2 - 2 * t), yo = t <= 0.5 ? Ly * 2 * t : Ly, rg = (g ** r - 1) / (g - 1);
    return [xi + rg * (xo - xi), yi + rg * (yo - yi)];
  };
  const mesh = gridMesh(nI, nJ, P, q8), { nodes, id, NI, NJ } = mesh, f = new Array(2 * nodes.length).fill(0), fixed = [];
  edgeLoad(f, nodes, N.range(NJ / 2 + 1, (j) => id[j][NI]), q8, i.sx_Pa, 0, i.t_m);
  edgeLoad(f, nodes, N.range(NJ / 2 + 1, (j) => id[NJ / 2 + j][NI]), q8, 0, i.sy_Pa, i.t_m);
  for (let k = 0; k <= NI; k++) { fixed.push(2 * id[0][k] + 1, 2 * id[NJ][k]); }
  return { mesh, f, fixed, a, b, Lx, Ly, q8 };
}
const PLANE_BASE = { elem: 'Q8', nMesh: 12, grade: 8, L_m: 0.3, W_m: 0.15, d_m: 0.02, ellip: 1, t_m: 0.003, sx_Pa: 100e6, sy_Pa: 0, material: 'Al 2024-T3', state: 'plane stress' };

const plane = {
  id: 'plane', title: 'Plane-stress finite elements: plate with a hole', fidelity: 'numerical',
  summary: 'Quadrilateral finite elements on a mapped mesh around a circular or elliptical cut-out under uniaxial or biaxial stress: stress contours on the deformed mesh, stress-concentration factor and comparison with the Kirsch and Inglis solutions.',
  equations: ['Linear elasticity equations', 'Strain–displacement compatibility equations', 'Generalized Hooke’s law', 'Finite element weak-form equilibrium equations', 'Minimum potential energy principle'],
  inputs: [
    { key: 'W_m', label: 'Plate width (across the load)', unit: 'm', default: 0.15, min: 0.002, group: 'Geometry' },
    { key: 'L_m', label: 'Plate length (along the load)', unit: 'm', default: 0.3, min: 0.002, group: 'Geometry' },
    { key: 'd_m', label: 'Hole width across the load', unit: 'm', default: 0.02, min: 0.0002, group: 'Geometry', help: 'Diameter for a circular hole' },
    { key: 'ellip', label: 'Hole length along load / width', unit: '-', default: 1, min: 0.1, max: 5, group: 'Geometry', help: '1 = circle; below 1 is a slot across the load (sharper, higher Kt)' },
    { key: 't_m', label: 'Thickness', unit: 'm', default: 0.003, min: 1e-5, group: 'Geometry' },
    { key: 'sx_Pa', label: 'Remote stress along the plate', unit: 'Pa', default: 100e6, group: 'Loads', help: 'Gross-section stress; defaults to the wing cover stress from the beam analysis' },
    { key: 'sy_Pa', label: 'Remote transverse stress', unit: 'Pa', default: 0, group: 'Loads', help: 'Biaxial case, e.g. hoop and axial stress around a fuselage cut-out' },
    MAT,
    { key: 'state', label: 'Stress state', type: 'select', options: ['plane stress', 'plane strain'], default: 'plane stress', group: 'Material' },
    { key: 'elem', label: 'Element', type: 'select', options: ['Q4', 'Q8', 'Q8 reduced'], default: 'Q8', group: 'Numerics', help: 'Q4: bilinear, 2×2 Gauss. Q8: quadratic serendipity, 3×3 Gauss (2×2 when reduced)' },
    { key: 'nMesh', label: 'Elements from hole to edge', unit: '', default: 12, min: 2, max: 32, step: 1, discrete: true, group: 'Numerics', help: 'The same count is used around the quarter hole; cost grows with the fourth power' },
    { key: 'grade', label: 'Radial grading (outer / inner element size)', unit: '-', default: 8, min: 1, max: 60, group: 'Numerics' },
  ],
  defaults: (c, up) => ({ material: c.struct.material, sx_Pa: up.fea?.sigma_bend_root_Pa, t_m: up.fea?.t_skin_root_m }),
  run(i) {
    const m = mat(i.material), mod = holeModel(i), { mesh, a, b, Lx, Ly } = mod, { nodes, id, NI, NJ } = mesh, warnings = [];
    const r = planeFE({ nodes, elems: mesh.elems, E: m.E, nu: m.nu, t: i.t_m, planeStrain: i.state === 'plane strain', reduced: i.elem === 'Q8 reduced', fixed: mod.fixed, f: mod.f });
    const pk = id[NJ][0], pk2 = id[0][0], sRef = Math.max(Math.abs(i.sx_Pa), Math.abs(i.sy_Pa)) || 1, circ = Math.abs(i.ellip - 1) < 1e-9;
    const top = r.sx[pk], side = r.sy[pk2], peak = Math.abs(top) >= Math.abs(side) ? top : side, Ktg = peak / sRef, net = Math.abs(top) >= Math.abs(side) ? 1 - b / Ly : 1 - a / Lx;
    const inf = (1 + 2 / i.ellip) * i.sx_Pa - i.sy_Pa; // Inglis/Kirsch stress at the end of the transverse axis, infinite plate
    const heywood = circ && i.sy_Pa === 0 ? ((2 + (1 - b / Ly) ** 3) / (1 - b / Ly)) * i.sx_Pa : NaN;
    const al = designAllowables(m), vmMax = N.amax(r.vm), mos = al.Sy / vmMax - 1, uEnd = r.u[2 * id[0][NI]];
    let Rx = 0; for (let k = 0; k <= NI; k++) Rx += r.resid[2 * id[NJ][k]];
    const eqErr = i.sx_Pa ? Math.abs(Rx + i.sx_Pa * Ly * i.t_m) / Math.abs(i.sx_Pa * Ly * i.t_m) : 0;
    if (b / Ly > 0.6 || a / Lx > 0.6) warnings.push('The hole removes more than 60% of the section: the mapped mesh is strongly distorted and net-section yielding governs.');
    if (vmMax > al.Sy) warnings.push('Peak elastic stress exceeds yield: the real material yields locally and redistributes; the elastic peak is still the correct input for fatigue (Suite 9).');
    if (i.elem === 'Q4' && i.nMesh < 16) warnings.push('Bilinear Q4 elements converge slowly at a stress raiser; run the convergence study or use Q8.');
    // sections for the line plots
    const lig = N.range(NI + 1, (k) => id[NJ][k]), yl = lig.map((k) => nodes[k][1]), hole = N.range(NJ + 1, (j) => id[j][0]);
    const th = hole.map((k) => N.deg(Math.atan2(nodes[k][1], nodes[k][0]))), hoop = hole.map((k, j) => {
      const ph = (j / NJ) * Math.PI / 2, tx = -a * Math.sin(ph), ty = b * Math.cos(ph), tn = Math.hypot(tx, ty), cx = tx / tn, cy = ty / tn;
      return (r.sx[k] * cx * cx + r.sy[k] * cy * cy + 2 * r.txy[k] * cx * cy) / 1e6;
    });
    const uMax = Math.max(...r.u.map(Math.abs)) || 1, sc = (0.05 * Math.max(Lx, Ly)) / uMax;
    const series1 = [{ name: 'FE σx on the net section', x: yl.map((v) => v * 1e3), y: lig.map((k) => r.sx[k] / 1e6), style: 'line+points' }], series2 = [{ name: 'FE hoop stress', x: th, y: hoop, style: 'line+points' }];
    if (circ) {
      series1.push({ name: 'Kirsch (infinite plate)', x: yl.map((v) => v * 1e3), y: yl.map((y) => (0.5 * i.sx_Pa * (2 + (b / y) ** 2 + 3 * (b / y) ** 4) + 0.5 * i.sy_Pa * ((b / y) ** 2 - 3 * (b / y) ** 4)) / 1e6), style: 'dash' });
      series2.push({ name: 'Kirsch (infinite plate)', x: th, y: th.map((t) => (i.sx_Pa * (1 - 2 * Math.cos(2 * N.rad(t))) + i.sy_Pa * (1 + 2 * Math.cos(2 * N.rad(t)))) / 1e6), style: 'dash' });
    }
    return {
      kpis: [
        kpi('Kt_gross', 'Stress-concentration factor (gross)', Ktg, '-', undefined, 'Peak stress / largest remote stress'),
        kpi('Kt_net', 'Stress-concentration factor (net section)', Ktg * net, '-'),
        kpi('sigma_peak_Pa', 'Peak stress at the hole edge', peak, 'Pa'),
        kpi('sigma_inf_plate_Pa', 'Infinite-plate reference (Kirsch / Inglis)', inf, 'Pa'),
        kpi('sigma_heywood_Pa', 'Finite-width reference (Heywood, empirical)', heywood, 'Pa', undefined, 'Circular hole under uniaxial load only'),
        kpi('vm_max_Pa', 'Peak von Mises stress', vmMax, 'Pa', mos >= 0 ? 'ok' : 'warn', `Yield ${(al.Sy / 1e6).toFixed(0)} MPa (${al.design ? 'design allowable' : 'typical'}); local yielding at a notch is a fatigue, not a static, concern`),
        kpi('mos_local_yield', 'Margin against first local yield', mos, '-'),
        kpi('elongation_m', 'End displacement of the quarter model', uEnd, 'm'),
        kpi('equilibrium_err', 'Reaction vs applied load error', eqErr, '-', eqErr < 1e-8 ? 'ok' : 'warn'),
        kpi('n_dof', 'Degrees of freedom', r.ndof, '-'),
        kpi('strain_energy_J', 'Strain energy of the quarter model', r.energy, 'J'),
      ],
      plots: [
        { type: 'tri', title: `von Mises stress on the deformed quarter model (displacement ×${sc.toPrecision(2)})`, xlabel: 'x [mm]', ylabel: 'y [mm]', zlabel: 'von Mises [MPa]', nodes: nodes.map((p, k) => [(p[0] + sc * r.u[2 * k]) * 1e3, (p[1] + sc * r.u[2 * k + 1]) * 1e3]), tris: quadTris(mesh.elems), values: r.vm.map((v) => v / 1e6), equalAspect: true, edges: mesh.elems.length <= 600 },
        { type: 'line', title: 'Stress across the net section', xlabel: 'Distance from the hole centre [mm]', ylabel: 'Stress σx [MPa]', series: series1 },
        { type: 'line', title: 'Hoop stress around the hole edge', xlabel: 'Angle from the load axis [deg]', ylabel: 'Tangential stress [MPa]', series: series2 },
      ],
      warnings,
      models: [`${i.elem} isoparametric quadrilaterals, ${mesh.elems.length} elements, ${r.ndof} DOF`, 'Banded Cholesky direct solver', 'Gauss-point stress recovery with nodal averaging'],
      assumptions: ['Quarter model with two symmetry planes', 'Linear elastic, small strain', `${i.state}; uniform remote traction on the plate ends`, 'Open (unloaded) hole: pin bearing and fastener load transfer are not modelled', matNote(m)],
    };
  },
  convergence: { param: 'nMesh', label: 'Elements from hole to edge', levels: [4, 8, 16, 24], metric: 'Kt_gross' },
  verify() {
    const out = [], m = mat('Al 2024-T3');
    // patch test: distorted rectangle under uniform tension must return the exact uniform stress and elongation
    for (const type of ['Q4', 'Q8', 'Q8 reduced']) {
      const q8 = type !== 'Q4', P = (r, t) => [2 * (r + 0.12 * Math.sin(Math.PI * r) * Math.sin(Math.PI * t)), t - 0.1 * Math.sin(Math.PI * r) * Math.sin(Math.PI * t)];
      const g = gridMesh(3, 3, P, q8), f = new Array(2 * g.nodes.length).fill(0), fixed = [];
      edgeLoad(f, g.nodes, N.range(g.NJ + 1, (j) => g.id[j][g.NI]), q8, 1e6, 0, 0.01);
      for (let j = 0; j <= g.NJ; j++) fixed.push(2 * g.id[j][0]); fixed.push(2 * g.id[0][0] + 1);
      const r = planeFE({ nodes: g.nodes, elems: g.elems, E: m.E, nu: m.nu, t: 0.01, reduced: type === 'Q8 reduced', fixed, f });
      out.push(N.check(`Patch test ${type}: uniform stress on a distorted mesh`, Math.max(...r.sx.map((v) => Math.abs(v - 1e6))) / 1e6 + Math.max(...r.sy.map(Math.abs)) / 1e6, 0, 1e-8, 'Irons patch test'));
      if (type === 'Q8') out.push(N.check('Patch test Q8: elongation σL/E', r.u[2 * g.id[0][g.NI]], (1e6 * 2) / m.E, 1e-9, 'Hooke’s law'));
    }
    const k1 = N.kv(plane.run({ ...PLANE_BASE, W_m: 0.4, L_m: 0.4, nMesh: 16, grade: 30 }));
    const k2 = N.kv(plane.run({ ...PLANE_BASE, W_m: 0.2, L_m: 0.5, nMesh: 16, grade: 20 }));
    const k3 = N.kv(plane.run({ ...PLANE_BASE, W_m: 0.4, L_m: 0.4, ellip: 0.5, nMesh: 18, grade: 100 }));
    const k4 = N.kv(plane.run({ ...PLANE_BASE, W_m: 0.4, L_m: 0.4, sy_Pa: 100e6, nMesh: 16, grade: 30 }));
    out.push(N.check('Kirsch Kt = 3 (hole 1/20 of the plate width)', k1.Kt_gross, 3, 0.015, 'Kirsch (1898); finite-size effect below 1%'));
    out.push(N.check('Finite width d/W = 0.1: Heywood Ktg = 3.03', k2.Kt_gross, (2 + 0.9 ** 3) / 0.9, 0.015, 'Heywood fit to Howland (empirical, ±1%)'));
    out.push(N.check('Inglis ellipse Kt = 1 + 2b/a = 5', k3.Kt_gross, 5, 0.015, 'Inglis (1913)'));
    out.push(N.check('Equal biaxial stress: Kt = 2', k4.Kt_gross, 2, 0.015, 'Kirsch superposition 3σ₁ − σ₂'));
    out.push(N.check('Global equilibrium of reactions', k2.equilibrium_err, 0, 1e-9, 'Statics'));
    return out;
  },
  validation: [{ name: 'Kirsch stress concentration of a small circular hole', source: 'Kirsch (1898), exact elasticity solution for an infinite plate: Kt = 3', inputs: { ...PLANE_BASE, W_m: 0.4, L_m: 0.4, grade: 30 }, sweep: { key: 'nMesh', values: [12, 16, 20] }, target: 'Kt_gross', observed: [3, 3, 3], tol_pct: 2 }],
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.Kt_gross > 3.5) out.push({ severity: 'warn', title: 'Sharp stress raiser', detail: `Kt = ${o.Kt_gross.toFixed(2)}: fatigue life falls roughly with the 3rd–4th power of local stress.`, action: 'Round the cut-out (larger corner radius, ellipse aligned with the load), add a doubler or thicken the land locally; pass Kt to Suite 9.', basis: 'Inglis: Kt = 1 + 2·(half-width across load)/(half-length along load)' });
    if (o.mos_local_yield < 0) out.push({ severity: 'advise', title: 'Local yielding at the hole edge', detail: `Elastic peak ${(o.sigma_peak_Pa / 1e6).toFixed(0)} MPa exceeds yield.`, action: 'Static strength is governed by net-section stress, but run the strain-life (Neuber) analysis in Suite 9 for durability; cold-working or interference-fit bushes extend life without added mass.', basis: 'Notch plasticity' });
    out.push({ severity: 'info', title: 'Discretisation check', detail: `${o.n_dof} DOF with ${i.elem}.`, action: 'Run the convergence study: the peak stress should change by less than 1–2% between the two finest meshes before it is used for fatigue.', basis: 'Grid convergence index' });
    return out;
  },
};

// ---- panel buckling -------------------------------------------------------------------------
/** Rayleigh–Ritz buckling of a simply supported plate a×b under Nx, Ny (compression positive) and Nxy with an M×M double-sine series. Returns the load multiplier and mode. */
function ritzBuckle(a, b, D, Nx, Ny, Nxy, M) {
  const idx = []; for (let p = 1; p <= M; p++) for (let q = 1; q <= M; q++) idx.push([p, q]);
  const n = idx.length, K = idx.map(([p, q]) => (D * a * b * Math.PI ** 4 * ((p / a) ** 2 + (q / b) ** 2) ** 2) / 4), C = N.zeros(n);
  for (let r = 0; r < n; r++) {
    const [m, nn] = idx[r];
    C[r][r] = ((a * b) / 4) * (Nx * ((m * Math.PI) / a) ** 2 + Ny * ((nn * Math.PI) / b) ** 2) / K[r];
    for (let s = r + 1; s < n; s++) { const [p, q] = idx[s]; if ((m + p) % 2 && (nn + q) % 2) C[r][s] = C[s][r] = (8 * Nxy * m * nn * p * q) / ((p * p - m * m) * (nn * nn - q * q)) / Math.sqrt(K[r] * K[s]); }
  }
  const e = N.eigSym(C), mu = e.values[n - 1];
  return { lambda: mu > 0 ? 1 / mu : Infinity, mode: e.vectors[n - 1].map((v, r) => v / Math.sqrt(K[r])), idx };
}
/** Euler–Johnson column strength [Pa] for slenderness L'/ρ and crippling (cut-off) stress. */
const eulerJohnson = (E, scc, sl) => (sl < Math.PI * Math.sqrt((2 * E) / scc) ? scc * (1 - (scc * sl * sl) / (4 * Math.PI ** 2 * E)) : (Math.PI ** 2 * E) / (sl * sl));

const buckling = {
  id: 'buckling', title: 'Skin–stringer panel buckling', fidelity: 'numerical',
  summary: 'Buckling of a stiffened cover panel under compression and shear: skin buckling between stringers by a Rayleigh–Ritz eigen-solution and classical coefficients, stringer local buckling, and Euler–Johnson column failure of the stringer with its effective skin.',
  equations: ['Kirchhoff–Love plate equations', 'Euler buckling equation', 'Geometric stiffness equations', 'Minimum potential energy principle'],
  inputs: [
    { key: 'a_m', label: 'Rib (or frame) pitch', unit: 'm', default: 0.6, min: 0.02, group: 'Geometry', help: 'Panel length along the load' },
    { key: 'b_m', label: 'Stringer pitch', unit: 'm', default: 0.15, min: 0.01, group: 'Geometry' },
    { key: 't_m', label: 'Skin thickness', unit: 'm', default: 0.003, min: 1e-4, group: 'Geometry' },
    { key: 'h_s', label: 'Stringer web height', unit: 'm', default: 0.04, min: 0.002, group: 'Stringer (Z-section)' },
    { key: 'b_f', label: 'Stringer flange width', unit: 'm', default: 0.02, min: 0.001, group: 'Stringer (Z-section)' },
    { key: 't_s', label: 'Stringer thickness', unit: 'm', default: 0.003, min: 1e-4, group: 'Stringer (Z-section)' },
    { key: 'fixity', label: 'Column end-fixity coefficient', unit: '-', default: 1, min: 1, max: 4, group: 'Geometry', help: '1 = pinned at ribs; 1.5–2 for typical rib restraint (judgement, calibrate by test)' },
    { key: 'edges', label: 'Skin edge support', type: 'select', options: ['simply supported', 'clamped'], default: 'simply supported', group: 'Geometry', help: 'Simply supported is the usual conservative choice between open-section stringers' },
    { key: 'sc_Pa', label: 'Applied compressive stress', unit: 'Pa', default: 150e6, min: 0, group: 'Loads', help: 'Limit-load cover stress; defaults to the beam analysis' },
    { key: 'sy_Pa', label: 'Applied transverse compressive stress', unit: 'Pa', default: 0, group: 'Loads' },
    { key: 'tau_Pa', label: 'Applied shear stress', unit: 'Pa', default: 20e6, min: 0, group: 'Loads' },
    { key: 'sf', label: 'Factor on applied loads', unit: '-', default: 1.5, min: 1, max: 3, group: 'Loads', help: '1.5 checks buckling at ultimate load; 1.0 at limit load' },
    { key: 'postbuckle', label: 'Allow skin to buckle (post-buckled design)', type: 'bool', default: false, group: 'Loads', help: 'Thin fuselage skins are commonly allowed to buckle below ultimate load; wing compression covers usually are not' },
    MAT,
    { key: 'nTerms', label: 'Ritz terms per direction', unit: '', default: 8, min: 2, max: 12, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c, up, d) {
    // a consistent starting panel: stringer pitch chosen so the skin just does not buckle at the factored stress, stringer and rib pitch scaled from it
    const wing = c.wing.S_m2 > 0, f = wing ? up.fea || {} : {}, m = mat(c.struct.material), sf = c.struct.sf_ultimate, sc = f.sigma_bend_root_Pa ?? 150e6;
    const t = f.t_skin_root_m ?? Math.min(c.struct.t_skin_mm, Math.max(0.4, 0.6 * (d.c_root || 5))) / 1e3, Ke = (Math.PI ** 2 * m.E) / (12 * (1 - m.nu ** 2));
    const b = N.clamp(0.85 * t * Math.sqrt((4 * Ke) / (sf * sc)), 0.012, 0.25), hs = Math.max(0.45 * b, 12 * t);
    return { material: c.struct.material, sf, sc_Pa: sc, tau_Pa: f.tau_skin_max_Pa ?? 0.15 * sc, t_m: t, b_m: b, h_s: hs, b_f: 0.4 * hs, t_s: Math.max(t, hs / 14), a_m: N.clamp(wing ? 0.12 * d.c_root : 10 * hs, 2 * b, Math.min(0.75, 9 * hs)) };
  },
  run(i) {
    const m = mat(i.material), Sy = designAllowables(m).Sy, { a_m: a, b_m: b, t_m: t } = i, D = (m.E * t ** 3) / (12 * (1 - m.nu ** 2)), Ke = (Math.PI ** 2 * m.E) / (12 * (1 - m.nu ** 2)), M = Math.max(2, Math.round(i.nTerms)), warnings = [];
    const sc = i.sc_Pa * i.sf, sy = i.sy_Pa * i.sf, tau = i.tau_Pa * i.sf, ar = a / b;
    // classical coefficients
    let kc = Infinity; for (let mm = 1; mm <= 40; mm++) kc = Math.min(kc, (mm / ar + ar / mm) ** 2);
    const clamped = i.edges === 'clamped', kcU = clamped ? (kc * 6.97) / 4 : kc, r1 = Math.max(ar, 1 / ar), bs = Math.min(a, b), ks = clamped ? 8.98 + 5.6 / r1 ** 2 : 5.34 + 4 / r1 ** 2;
    const scr = kcU * Ke * (t / b) ** 2, tcr = ks * Ke * (t / bs) ** 2, Rc = sc / scr, Rs = tau / tcr;
    const lamInt = Rs > 0 ? (-Rc + Math.sqrt(Rc * Rc + 4 * Rs * Rs)) / (2 * Rs * Rs) : Rc > 0 ? 1 / Rc : Infinity;
    const rz = ritzBuckle(a, b, D, sc * t, sy * t, tau * t, M), lamSkinE = clamped ? lamInt : rz.lambda;
    // plasticity cut-off: elastic buckling stress cannot exceed yield
    const vmApp = Math.sqrt(sc * sc - sc * sy + sy * sy + 3 * tau * tau), lamY = vmApp > 0 ? Sy / vmApp : Infinity, lamSkin = Math.min(lamSkinE, lamY);
    if (lamSkinE > lamY) warnings.push('Elastic skin buckling stress exceeds yield: the skin is stocky and the result is capped at yield (no plasticity-corrected buckling curve is applied).');
    // stringer: local buckling of flange (one edge free, k = 0.43) and web (k = 4) taken as the crippling cut-off
    const sFl = 0.43 * Ke * (i.t_s / i.b_f) ** 2, sWeb = 4 * Ke * (i.t_s / i.h_s) ** 2, scc = Math.min(Sy, sFl, sWeb);
    // column of stringer + effective skin (von Kármán–Sechler effective width, iterated on the column stress)
    let scol = scc, we = b, sec;
    for (let it = 0; it < 30; it++) {
      we = lamSkinE * sc >= scol || !(sc > 0) ? b : Math.min(b, 1.7 * t * Math.sqrt(m.E / scol));
      const parts = [[we * t, 0], [i.b_f * i.t_s, t / 2 + i.t_s / 2], [i.h_s * i.t_s, t / 2 + i.h_s / 2], [i.b_f * i.t_s, t / 2 + i.h_s - i.t_s / 2]], own = [(we * t ** 3) / 12, (i.b_f * i.t_s ** 3) / 12, (i.t_s * i.h_s ** 3) / 12, (i.b_f * i.t_s ** 3) / 12];
      const A = N.sum(parts.map((p) => p[0])), zb = N.sum(parts.map((p) => p[0] * p[1])) / A, I = N.sum(parts.map((p, k) => own[k] + p[0] * (p[1] - zb) ** 2));
      sec = { A, I, rho: Math.sqrt(I / A), zb };
      const nw = eulerJohnson(m.E, scc, a / Math.sqrt(i.fixity) / sec.rho);
      if (Math.abs(nw - scol) < 1e-6 * scol) { scol = nw; break; }
      scol = 0.5 * (scol + nw);
    }
    const sl = a / Math.sqrt(i.fixity) / sec.rho, lamCol = sc > 0 ? scol / sc : Infinity, lamLoc = sc > 0 ? scc / sc : Infinity;
    const gov = i.postbuckle ? Math.min(lamCol, lamLoc) : Math.min(lamSkin, lamCol, lamLoc), govName = gov === lamCol ? 'stringer column' : gov === lamLoc ? 'stringer local buckling' : 'skin buckling';
    if (!Number.isFinite(gov)) warnings.push('No compressive or shear load is applied: buckling factors are infinite.');
    if (b / t < 15) warnings.push('Very stocky skin (b/t < 15): thin-plate buckling theory is not the governing failure mode.');
    // mode shape
    const xs = N.linspace(0, a, 41), ys = N.linspace(0, b, 21), z = ys.map((y) => xs.map((x) => { let w = 0; rz.idx.forEach(([p, q], r) => { w += rz.mode[r] * Math.sin((p * Math.PI * x) / a) * Math.sin((q * Math.PI * y) / b); }); return w; }));
    const zm = Math.max(...z.flat().map(Math.abs)) || 1, sls = N.linspace(5, 150, 60);
    const tAr = (b * t + (2 * i.b_f + i.h_s) * i.t_s) / b;
    return {
      kpis: [
        kpi('buckling_factor', `Governing buckling factor (${govName})`, gov, '-', st3(gov - 1, 0.1), `Critical load / (applied × ${i.sf}); must be ≥ 1`),
        kpi('lambda_skin', 'Skin buckling factor', lamSkin, '-', i.postbuckle ? undefined : st3(lamSkin - 1, 0.1), clamped ? 'Classical clamped coefficients (k = 6.97 scaled for aspect ratio; 8.98 + 5.6/(a/b)² in shear) with interaction Rc + Rs² = 1; transverse stress ignored' : 'Rayleigh–Ritz eigenvalue, capped at yield'),
        kpi('lambda_ritz', 'Rayleigh–Ritz eigenvalue (simply supported, elastic)', rz.lambda, '-'),
        kpi('lambda_interaction', 'Classical interaction estimate', lamInt, '-'),
        kpi('lambda_column', 'Stringer column factor', lamCol, '-', st3(lamCol - 1, 0.1)),
        kpi('lambda_local', 'Stringer local-buckling factor', lamLoc, '-', st3(lamLoc - 1, 0.1)),
        kpi('sigma_cr_skin_Pa', 'Skin compression buckling stress', Math.min(scr, Sy), 'Pa'),
        kpi('tau_cr_skin_Pa', 'Skin shear buckling stress', Math.min(tcr, Sy / Math.sqrt(3)), 'Pa'),
        kpi('k_compression', 'Compression buckling coefficient', kcU, '-'),
        kpi('k_shear', 'Shear buckling coefficient', ks, '-'),
        kpi('sigma_column_Pa', 'Column failure stress', scol, 'Pa'),
        kpi('sigma_crippling_Pa', 'Crippling cut-off stress', scc, 'Pa', undefined, 'Lowest of yield and flange/web local buckling (conservative)'),
        kpi('slenderness', 'Column slenderness L′/ρ', sl, '-'),
        kpi('eff_width_m', 'Effective skin width', we, 'm'),
        kpi('t_smeared_m', 'Smeared panel thickness', tAr, 'm'),
      ],
      plots: [
        { type: 'heat', title: 'Critical skin buckling mode (normalised deflection)', xlabel: 'Along the load x [m]', ylabel: 'Across the bay y [m]', zlabel: 'w / w_max', x: xs, y: ys, z: z.map((r) => r.map((v) => v / zm)), contours: 10, diverging: true },
        { type: 'line', title: 'Euler–Johnson column curve', xlabel: 'Slenderness L′/ρ [-]', ylabel: 'Failure stress [MPa]', series: [{ name: 'Euler–Johnson', x: sls, y: sls.map((s) => eulerJohnson(m.E, scc, s) / 1e6) }, { name: 'Euler', x: sls, y: sls.map((s) => Math.min(2 * Sy, (Math.PI ** 2 * m.E) / (s * s)) / 1e6), style: 'dash' }, { name: 'This panel', x: [sl], y: [scol / 1e6], style: 'points' }], annotations: [{ y: sc / 1e6, label: 'Applied' }] },
        { type: 'bar', title: 'Buckling factors by mode', ylabel: 'Critical / applied [-]', categories: ['Skin', 'Stringer local', 'Stringer column'], series: [{ name: 'Factor', y: [lamSkin, lamLoc, lamCol].map((v) => (Number.isFinite(v) ? Math.min(v, 20) : 20)) }] },
      ],
      warnings,
      models: [`Rayleigh–Ritz double-sine series (${M}×${M} terms) with symmetric Jacobi eigen-solution`, 'Classical plate buckling coefficients', 'von Kármán–Sechler effective width (semi-empirical)', 'Euler–Johnson column curve'],
      assumptions: ['Flat, perfect, isotropic panel; no curvature benefit, imperfection knock-down or plasticity correction below yield', 'Stringer torsional restraint of the skin is ignored when simply supported', 'Default stringer and rib pitches and stringer proportions come from a sizing rule (starting point only)', 'Crippling taken as the lower of yield and elastic local buckling instead of an empirical crippling curve', 'Tension yield is used as the compression cut-off (compression yield Fcy of sheet is usually a little lower)', matNote(m)],
    };
  },
  convergence: { param: 'nTerms', label: 'Ritz terms per direction', levels: [2, 4, 6, 8, 10], metric: 'lambda_ritz' },
  verify() {
    const E = 70e9, nu = 0.3, t = 0.002, b = 0.2, D = (E * t ** 3) / (12 * (1 - nu * nu)), Ncr = (k) => (k * Math.PI ** 2 * D) / (b * b);
    const sq = ritzBuckle(b, b, D, Ncr(1), 0, 0, 6), lg = ritzBuckle(3 * b, b, D, Ncr(1), 0, 0, 6), sh = ritzBuckle(b, b, D, 0, 0, Ncr(1), 10), bi = ritzBuckle(b, b, D, Ncr(1), Ncr(1), 0, 4);
    return [
      N.check('Square plate in compression: k = 4', sq.lambda, 4, 1e-9, 'Bryan (1891); Timoshenko & Gere, Theory of Elastic Stability'),
      N.check('Plate a/b = 3 in compression: k = 4 (m = 3)', lg.lambda, 4, 1e-9, 'Timoshenko & Gere'),
      N.check('Square plate in shear: k = 9.34', sh.lambda, 9.34, 5e-3, 'Stein & Neff (1947); Timoshenko & Gere'),
      N.check('Square plate, equal biaxial compression: k = 2', bi.lambda, 2, 1e-9, 'Timoshenko & Gere'),
      N.check('Euler column σ = π²E/(L/ρ)²', eulerJohnson(E, 300e6, 120), (Math.PI ** 2 * E) / 120 ** 2, 1e-12, 'Euler (1757)'),
      N.check('Johnson parabola meets Euler at σcc/2', eulerJohnson(E, 300e6, Math.PI * Math.sqrt((2 * E) / 300e6) * 0.999999), 150e6, 1e-5, 'Johnson parabola tangency'),
    ];
  },
  validation: [{ name: 'Simply supported long plate in compression, k = 4', source: 'Timoshenko & Gere, Theory of Elastic Stability: σcr = 4π²E/(12(1−ν²))·(t/b)²', inputs: { a_m: 0.6, b_m: 0.2, edges: 'simply supported', sc_Pa: 1e6, sy_Pa: 0, tau_Pa: 0, sf: 1, material: 'Al 7075-T6', nTerms: 6, postbuckle: false }, sweep: { key: 't_m', values: [0.001, 0.002, 0.003] }, target: 'lambda_ritz', observed: [0.001, 0.002, 0.003].map((t) => (4 * Math.PI ** 2 * 71.7e9 * (t / 0.2) ** 2) / (12 * (1 - 0.33 ** 2)) / 1e6), tol_pct: 0.1 }],
  calibration: { params: [{ key: 'fixity', min: 1, max: 4 }], sweep: 'a_m', target: 'sigma_column_Pa', note: 'Panel compression test failure stresses versus rib pitch calibrate the end-fixity coefficient.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.buckling_factor < 1) out.push({ severity: 'critical', title: 'Panel buckles before the required load', detail: `Buckling factor ${o.buckling_factor.toFixed(2)} at ${i.sf}× applied load.`, action: o.lambda_column <= o.buckling_factor * 1.001 ? 'Reduce the rib pitch or use a deeper stringer (higher radius of gyration).' : o.lambda_local <= o.buckling_factor * 1.001 ? 'Thicken the stringer or shorten its free flange.' : 'Reduce stringer pitch or thicken the skin: buckling stress scales with (t/b)².', basis: 'Buckling factor ≥ 1 at ultimate load for compression-critical covers' });
    else if (o.buckling_factor > 2 && Number.isFinite(o.buckling_factor)) out.push({ severity: 'advise', title: 'Stability margin is generous', detail: `Buckling factor ${o.buckling_factor.toFixed(2)}.`, action: 'A wider stringer or rib pitch removes parts and mass; re-check fatigue and damage tolerance before thinning.', basis: 'Minimum-mass stiffened-panel design: skin, local and column modes should be near-coincident' });
    if (!i.postbuckle && o.lambda_skin < 1 && o.lambda_column > 1) out.push({ severity: 'advise', title: 'Skin buckles but stringers hold', detail: 'A post-buckled (tension-field) design may be acceptable for shear-dominated fuselage panels.', action: 'Enable the post-buckled option only where skin buckling is permitted, and check permanent set at limit load.', basis: 'Effective-width post-buckling design' });
    return out;
  },
};

// ---- pressurised fuselage -------------------------------------------------------------------
const fuselage = {
  id: 'fuselage', title: 'Pressurised fuselage shell with frames', fidelity: 'numerical',
  summary: 'Hoop and longitudinal stress of the pressure cabin, the bending disturbance where frames restrain the skin (axisymmetric shell strip solved by finite elements), frame stress, thick-wall check and minimum skin gauge.',
  equations: ['Classical shell equations', 'Static equilibrium equations', 'Hooke’s law', 'Finite element weak-form equilibrium equations'],
  applicable: (c) => (c.fuselage.cabin_dp_Pa > 0 && c.fuselage.dia_m > 0 ? true : 'The current case has no cabin pressure differential.'),
  inputs: [
    { key: 'dp_Pa', label: 'Cabin pressure differential', unit: 'Pa', default: 57000, min: 100, group: 'Loads', help: 'Maximum relief-valve setting' },
    { key: 'kp', label: 'Pressure factor (ultimate)', unit: '-', default: 1.33, min: 1, max: 2, group: 'Loads', help: 'Transport practice applies 1.33 on pressure alone before the 1.5 ultimate factor; confirm against your certification basis' },
    { key: 'sf', label: 'Ultimate safety factor', unit: '-', default: 1.5, min: 1, max: 3, group: 'Loads' },
    { key: 'M_Nm', label: 'Fuselage bending moment at this station', unit: 'N·m', default: 0, group: 'Loads', help: 'Limit moment from tail and inertia loads; adds an axial stress M·R/I' },
    { key: 'R_m', label: 'Fuselage radius', unit: 'm', default: 1.975, min: 0.05, group: 'Geometry' },
    { key: 't_m', label: 'Skin thickness', unit: 'm', default: 0.0012, min: 1e-4, group: 'Geometry' },
    { key: 'Lf_m', label: 'Frame pitch', unit: 'm', default: 0.5, min: 0.05, group: 'Geometry' },
    { key: 'Af_m2', label: 'Frame cross-section area', unit: 'm²', default: 2e-4, min: 0, group: 'Geometry' },
    { key: 'Astr_m2', label: 'Total stringer area around the section', unit: 'm²', default: 0.008, min: 0, group: 'Geometry', help: 'Carries part of the longitudinal load' },
    { key: 's_hoop_allow', label: 'Design hoop stress (fatigue-driven)', unit: 'Pa', default: 100e6, min: 1e6, group: 'Material', help: 'Typically 80–110 MPa for aluminium pressure cabins; a durability choice, not a strength limit' },
    MAT,
    { key: 'nElem', label: 'Elements in a half bay', unit: '', default: 40, min: 4, max: 400, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c) {
    const R = c.fuselage.dia_m / 2, p = c.fuselage.cabin_dp_Pa, t = Math.max(0.0008, (p * R) / 100e6);
    return { dp_Pa: p, R_m: R, t_m: t, material: c.struct.material, sf: c.struct.sf_ultimate, Lf_m: N.clamp(0.25 * R, 0.15, 0.55), Af_m2: 0.3 * t * N.clamp(0.25 * R, 0.15, 0.55), Astr_m2: 0.5 * 2 * Math.PI * R * t };
  },
  run(i) {
    const m = mat(i.material), { dp_Pa: p, R_m: R, t_m: t } = i, n = Math.max(4, Math.round(i.nElem)), nu = m.nu, warnings = [];
    const D = (m.E * t ** 3) / (12 * (1 - nu * nu)), x = N.linspace(0, i.Lf_m / 2, n + 1), Ask = 2 * Math.PI * R * t;
    const sxM = (p * Math.PI * R * R) / (Ask + i.Astr_m2), pe = p - (nu * sxM * t) / R; // Poisson contraction from the axial stress reduces the radial load
    const sol = beamSolve(x, new Array(n).fill(D), new Array(n).fill(Infinity), new Array(n + 1).fill(pe), [], { fixed: [1, 2 * n + 1], springs: [[0, (m.E * i.Af_m2) / (2 * R * R)]], kf: new Array(n).fill((m.E * t) / (R * R)) });
    const w = sol.w, Mb = sol.M, hoop = w.map((v) => (m.E * v) / R + nu * sxM), sb = Mb.map((M) => (6 * M) / (t * t));
    const sxB = i.M_Nm / (Math.PI * R * R * (t + i.Astr_m2 / (2 * Math.PI * R))), sxO = sb.map((v) => sxM + sxB + Math.abs(v)), hO = hoop.map((h, j) => h + nu * Math.abs(sb[j]));
    const vm = sxO.map((s, j) => Math.sqrt(s * s - s * hO[j] + hO[j] ** 2)), vmMax = N.amax(vm), hoopU = (p * R) / t, sFrame = (m.E * w[0]) / R;
    const al = designAllowables(m), mos = al.Su / (i.kp * i.sf * vmMax) - 1, mosY = al.Sy / (i.kp * vmMax) - 1, ri = R - t / 2, ro = R + t / 2, lame = (p * (ro * ro + ri * ri)) / (ro * ro - ri * ri);
    const tMin = (p * R) / i.s_hoop_allow, beta = (3 * (1 - nu * nu) / (R * R * t * t)) ** 0.25, mPerM = m.rho * (Ask + i.Astr_m2 + (2 * Math.PI * R * i.Af_m2) / i.Lf_m);
    if (R / t < 10) warnings.push('R/t < 10: thin-shell theory is inaccurate; use the Lamé thick-wall value.');
    if (hoop[n] > i.s_hoop_allow) warnings.push(`Mid-bay hoop stress ${(hoop[n] / 1e6).toFixed(0)} MPa exceeds the design hoop stress: fatigue life and crack-arrest capability of the cabin are at risk.`);
    if (i.Lf_m * beta < 3) warnings.push('Frames are closer than the shell bending length: the disturbances from adjacent frames overlap (this is captured by the model).');
    const xm = x.map((v) => v * 1e3);
    return {
      kpis: [
        kpi('hoop_unstiffened_Pa', 'Hoop stress pR/t (no frames)', hoopU, 'Pa'),
        kpi('hoop_midbay_Pa', 'Hoop stress at mid-bay', hoop[n], 'Pa', hoop[n] <= i.s_hoop_allow ? 'ok' : 'warn', `Design hoop stress ${(i.s_hoop_allow / 1e6).toFixed(0)} MPa`),
        kpi('hoop_frame_Pa', 'Skin hoop stress at the frame', hoop[0], 'Pa'),
        kpi('sigma_long_Pa', 'Longitudinal membrane stress', sxM + sxB, 'Pa'),
        kpi('sigma_bend_frame_Pa', 'Skin bending stress at the frame', Math.abs(sb[0]), 'Pa', undefined, 'Local discontinuity stress: a fatigue driver at the frame shear ties'),
        kpi('sigma_frame_Pa', 'Frame hoop stress', sFrame, 'Pa'),
        kpi('sigma_vm_max_Pa', 'Peak skin von Mises stress', vmMax, 'Pa'),
        kpi('mos_pressure', 'Margin of safety under factored pressure', Math.min(mos, mosY), '-', st3(Math.min(mos, mosY)), al.design ? 'Design allowables' : 'Typical strengths, not design allowables'),
        kpi('radial_growth_m', 'Radial growth at mid-bay', w[n], 'm'),
        kpi('lame_hoop_Pa', 'Thick-wall (Lamé) hoop stress at the inner surface', lame, 'Pa'),
        kpi('t_min_m', 'Minimum skin gauge for the design hoop stress', tMin, 'm', t >= tMin ? 'ok' : 'warn'),
        kpi('shell_mass_kgm', 'Shell mass per metre of cabin', mPerM, 'kg/m'),
      ],
      plots: [
        { type: 'line', title: 'Radial displacement between frame and mid-bay', xlabel: 'Distance from the frame [mm]', ylabel: 'Radial displacement [mm]', series: [{ name: 'w', x: xm, y: w.map((v) => v * 1e3) }], annotations: [{ y: ((hoopU - nu * sxM) * R / m.E) * 1e3, label: 'Unstiffened' }] },
        { type: 'line', title: 'Skin stresses along the bay', xlabel: 'Distance from the frame [mm]', ylabel: 'Stress [MPa]', series: [{ name: 'Hoop (membrane)', x: xm, y: hoop.map((v) => v / 1e6) }, { name: 'Axial, surface (membrane + bending)', x: xm, y: sxO.map((v) => v / 1e6) }, { name: 'Axial bending only', x: xm, y: sb.map((v) => v / 1e6), style: 'dash' }, { name: 'von Mises', x: xm, y: vm.map((v) => v / 1e6), style: 'dash' }] },
      ],
      warnings,
      models: ['Axisymmetric cylindrical shell strip (beam on elastic foundation) with Hermite elements', 'Frame as a ring spring E·A/R²', 'Lamé thick-wall solution'],
      assumptions: ['Circular section, uniform skin, frames attached continuously to the skin', 'Closed-end axial load shared by skin and stringers', 'The 1.33 pressure factor and default frame and stringer proportions are typical transport practice: confirm against the certification basis', 'No cut-outs, doors, lap joints, bulkheads or floor-beam effects', 'The design hoop stress is a typical durability value, not a sourced limit', matNote(m)],
    };
  },
  convergence: { param: 'nElem', label: 'Elements in a half bay', levels: [5, 10, 20, 40, 80], metric: 'sigma_bend_frame_Pa' },
  verify() {
    const b = { dp_Pa: 60000, kp: 1.33, sf: 1.5, M_Nm: 0, R_m: 2, t_m: 0.002, Lf_m: 2, Af_m2: 1e3, Astr_m2: 0, s_hoop_allow: 100e6, material: 'Al 2024-T3', nElem: 200 }, m = mat(b.material), r = N.kv(fuselage.run(b)), pr = (60000 * 2) / 0.002;
    return [
      N.check('Mid-bay hoop stress far from a frame: pR/t', r.hoop_midbay_Pa, pr, 1e-6, 'Membrane theory'),
      N.check('Longitudinal stress pR/2t', r.sigma_long_Pa, pr / 2, 1e-12, 'Membrane theory'),
      N.check('Bending stress at a rigid ring: √(3/(1−ν²))·(1−ν/2)·pR/t', r.sigma_bend_frame_Pa, Math.sqrt(3 / (1 - m.nu ** 2)) * (1 - m.nu / 2) * pr, 2e-4, 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells, §115'),
      N.check('Lamé → pR/t for a thin wall', r.lame_hoop_Pa, pr, 1e-3, 'Lamé (1852)'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.t_min_m > i.t_m) out.push({ severity: 'warn', title: 'Skin is thinner than the durability gauge', detail: `Needs ${(o.t_min_m * 1e3).toFixed(2)} mm for a ${(i.s_hoop_allow / 1e6).toFixed(0)} MPa hoop stress; ${(i.t_m * 1e3).toFixed(2)} mm specified.`, action: 'Thicken the skin, lower the cabin differential (higher cabin altitude), or justify a higher hoop stress with fatigue and crack-arrest evidence in Suite 9.', basis: 'Hoop stress pR/t against the fatigue-driven design value' });
    if (o.mos_pressure < 0) out.push({ severity: 'critical', title: 'Negative static margin under factored pressure', detail: `Margin ${o.mos_pressure.toFixed(2)}.`, action: 'Increase skin thickness or reduce the differential.', basis: 'Ultimate pressure case' });
    if (o.sigma_bend_frame_Pa > 0.5 * o.hoop_midbay_Pa) out.push({ severity: 'advise', title: 'Strong frame "pillowing" stress', detail: `Skin bending at the frame is ${(o.sigma_bend_frame_Pa / 1e6).toFixed(0)} MPa.`, action: 'A lighter (more compliant) frame or shear-tied frame with tear straps reduces the discontinuity; this location needs a fatigue check.', basis: 'Shell discontinuity stress' });
    out.push({ severity: 'info', title: 'Pressure-cabin mass lever', detail: `Shell mass ≈ ${o.shell_mass_kgm.toFixed(1)} kg per metre. Skin gauge scales directly with differential pressure and radius.`, action: 'Each 10% reduction in differential pressure saves about 10% of pressure-sized skin mass, with fuel and CO₂ savings every flight, at the cost of a higher cabin altitude.', basis: 't = pR/σ_hoop' });
    return out;
  },
};

// ---- 3-D space truss (direct stiffness) -----------------------------------------------------
/** Pin-jointed space truss. nodes [[x,y,z]], members [[i,j,A]], fixed node indices, loads [[node, Fx, Fy, Fz]]. */
function truss3d(nodes, members, fixed, loads, E) {
  const nd = 3 * nodes.length, K = N.zeros(nd), f = new Array(nd).fill(0), geo = [];
  for (const [a, b, A] of members) {
    const d = N.vadd(nodes[b], nodes[a], -1), L = N.norm(d), c = d.map((v) => v / L), k = (E * A) / L; geo.push({ L, c, k });
    for (let p = 0; p < 3; p++) for (let q = 0; q < 3; q++) { const v = k * c[p] * c[q]; K[3 * a + p][3 * a + q] += v; K[3 * b + p][3 * b + q] += v; K[3 * a + p][3 * b + q] -= v; K[3 * b + p][3 * a + q] -= v; }
  }
  for (const [nn, fx, fy, fz] of loads) { f[3 * nn] += fx; f[3 * nn + 1] += fy; f[3 * nn + 2] += fz; }
  const free = []; for (let k = 0; k < nd; k++) if (!fixed.includes(Math.floor(k / 3))) free.push(k);
  const uf = N.solve(free.map((r) => free.map((c) => K[r][c])), free.map((r) => f[r])), u = new Array(nd).fill(0);
  free.forEach((d, k) => { u[d] = uf[k]; });
  const force = members.map(([a, b], e) => geo[e].k * N.dot(geo[e].c, [u[3 * b] - u[3 * a], u[3 * b + 1] - u[3 * a + 1], u[3 * b + 2] - u[3 * a + 2]]));
  const R = N.matvec(K, u).map((v, k) => v - f[k]);
  return { u, force, R, geo };
}
/** Engine-mount geometry: 4 firewall nodes (fixed, x = 0) and 4 engine-ring nodes (x = L). */
function mountModel(i) {
  const a = i.a_m / 2, c = i.c_m / 2, sg = [[1, 1], [-1, 1], [-1, -1], [1, -1]], nodes = [...sg.map(([y, z]) => [0, a * y, a * z]), ...sg.map(([y, z]) => [i.L_m, c * y, c * z])];
  const At = Math.PI * (i.D_m - i.tw_m) * i.tw_m, members = [], names = [];
  for (let k = 0; k < 4; k++) { members.push([k, 4 + k, At]); names.push(`Longeron ${k + 1}`); }
  for (let k = 0; k < 4; k++) { members.push([(k + 1) % 4, 4 + k, At]); names.push(`Diagonal ${k + 1}`); }
  for (let k = 0; k < 4; k++) { members.push([4 + k, 4 + ((k + 1) % 4), At]); names.push(`Ring ${k + 1}`); }
  members.push([4, 6, At]); names.push('Ring brace');
  return { nodes, members, names, At };
}
const frame = {
  id: 'frame', title: 'Engine or motor mount space truss', fidelity: 'numerical',
  summary: 'Direct-stiffness solution of a welded-tube mount truss between the firewall and the engine ring under inertia, thrust and torque: member forces, tube stress and Euler buckling utilisation, attachment reactions.',
  equations: ['Static equilibrium equations', 'Hooke’s law', 'Finite element weak-form equilibrium equations', 'Euler buckling equation'],
  inputs: [
    { key: 'm_eng_kg', label: 'Engine / motor mass (with accessories)', unit: 'kg', default: 150, min: 0.01, group: 'Loads' },
    { key: 'n_z', label: 'Vertical limit load factor', unit: 'g', default: 3.8, min: 0, max: 20, group: 'Loads' },
    { key: 'n_y', label: 'Side load factor', unit: 'g', default: 1.33, min: 0, max: 10, group: 'Loads', help: 'Light-aeroplane practice uses about one third of the vertical factor, not less than 1.33; confirm against your basis' },
    { key: 'T_N', label: 'Thrust', unit: 'N', default: 2600, min: 0, group: 'Loads' },
    { key: 'Q_Nm', label: 'Engine torque reaction', unit: 'N·m', default: 420, min: 0, group: 'Loads', help: 'Shaft power / shaft speed × torque factor (≥ 1.25–2 for piston engines)' },
    { key: 'x_cg', label: 'Engine CG ahead of the ring', unit: 'm', default: 0.25, min: 0, group: 'Geometry' },
    { key: 'L_m', label: 'Mount length (firewall to ring)', unit: 'm', default: 0.45, min: 0.02, group: 'Geometry' },
    { key: 'a_m', label: 'Firewall attachment spacing', unit: 'm', default: 0.6, min: 0.02, group: 'Geometry' },
    { key: 'c_m', label: 'Engine-ring attachment spacing', unit: 'm', default: 0.4, min: 0.02, group: 'Geometry' },
    { key: 'D_m', label: 'Tube outer diameter', unit: 'm', default: 0.019, min: 0.002, group: 'Tubes' },
    { key: 'tw_m', label: 'Tube wall thickness', unit: 'm', default: 0.0012, min: 0.0002, group: 'Tubes' },
    { ...MAT, default: 'Steel 4340 (QT)' },
    { key: 'sf', label: 'Ultimate safety factor', unit: '-', default: 1.5, min: 1, max: 3, group: 'Loads' },
  ],
  defaults(c) {
    const p = c.prop, m = Math.max(0.05, engineMass(p)), nz = c.aero.n_pos, s = N.clamp(0.11 * m ** (1 / 3), 0.03, 1.6), om = p.rpm > 0 ? (p.rpm * 2 * Math.PI) / 60 : 0, Q = om && p.P0_W ? p.P0_W / om : 0;
    // tube sized so that the direct stress in a longeron is of order 150 MPa under the pitching couple (starting point only)
    const Fm = (nz * m * G0 * (1 + (2 * 0.6 * s) / (0.8 * s)) + p.T0_N) / 2, D = Math.max(0.004, Math.sqrt((12 * Fm) / (Math.PI * 150e6)));
    const o = { m_eng_kg: m, n_z: nz, n_y: Math.max(1.33, nz / 3), T_N: p.T0_N, Q_Nm: 1.5 * Q, x_cg: 0.6 * s, L_m: s, a_m: 1.2 * s, c_m: 0.8 * s, D_m: D, tw_m: Math.max(0.0005, D / 12), sf: c.struct.sf_ultimate, material: 'Steel 4340 (QT)' };
    for (let k = 0; k < 12; k++) { // resize the tube until the peak utilisation is about 0.7
      const u = N.kv(frame.run(o)).util_max; if (!(u > 0) || (u > 0.6 && u < 0.8)) break;
      o.D_m = Math.max(0.004, o.D_m * N.clamp((u / 0.7) ** 0.4, 0.6, 1.8)); o.tw_m = Math.max(0.0005, o.D_m / 12);
    }
    return o;
  },
  run(i) {
    const m = mat(i.material), al = designAllowables(m), { nodes, members, names, At } = mountModel(i), c = i.c_m / 2, W = i.m_eng_kg * G0, warnings = [];
    // rigid engine: resolve CG forces and moments into the four ring nodes
    const Fz = -i.n_z * W, Fy = i.n_y * W, My = -Fz * i.x_cg, Mz = Fy * i.x_cg, loads = [];
    for (let k = 4; k < 8; k++) { const [, y, z] = nodes[k]; loads.push([k, i.T_N / 4 + (My * z) / (4 * c * c) - (Mz * y) / (4 * c * c), Fy / 4 - (i.Q_Nm * z) / (8 * c * c), Fz / 4 + (i.Q_Nm * y) / (8 * c * c)]); }
    let sol;
    try { sol = truss3d(nodes, members, [0, 1, 2, 3], loads, m.E); } catch { throw new Error('The mount truss is a mechanism for this geometry; check the attachment spacings.'); }
    const I = (Math.PI / 64) * (i.D_m ** 4 - (i.D_m - 2 * i.tw_m) ** 4), stress = sol.force.map((F) => F / At), Pcr = sol.geo.map((g) => (Math.PI ** 2 * m.E * I) / (g.L * g.L));
    const rho = Math.sqrt(I / At), util = sol.force.map((F, e) => { const sAllow = F >= 0 ? al.Su : Math.min(eulerJohnson(m.E, al.Sy, sol.geo[e].L / rho), al.Su); return (i.sf * Math.abs(F)) / (sAllow * At); });
    const je = N.argmax(util), uMax = util[je], mos = 1 / uMax - 1, disp = Math.max(...N.range(4, (k) => Math.hypot(sol.u[12 + 3 * k], sol.u[13 + 3 * k], sol.u[14 + 3 * k])));
    const Rsum = [0, 1, 2].map((p) => N.sum([0, 1, 2, 3].map((k) => sol.R[3 * k + p]))), eq = Math.hypot(Rsum[0] + i.T_N, Rsum[1] + Fy, Rsum[2] + Fz) / Math.max(1, Math.hypot(i.T_N, Fy, Fz));
    const Rmax = Math.max(...[0, 1, 2, 3].map((k) => Math.hypot(sol.R[3 * k], sol.R[3 * k + 1], sol.R[3 * k + 2]))), mass = m.rho * At * N.sum(sol.geo.map((g) => g.L));
    if (i.D_m / i.tw_m > 60) warnings.push('Tube D/t above 60: local wall buckling can precede column buckling.');
    if (uMax > 1) warnings.push(`${names[je]} is over-utilised (${uMax.toFixed(2)}) at ultimate load.`);
    const sc = disp > 0 ? (0.05 * i.L_m) / disp : 0, side = (def) => { const x = [], z = []; for (const [a, b] of members) { for (const k of [a, b]) { x.push(nodes[k][0] + (def ? sc * sol.u[3 * k] : 0)); z.push(nodes[k][2] + (def ? sc * sol.u[3 * k + 2] : 0)); } x.push(NaN); z.push(NaN); } return { x, y: z }; };
    return {
      kpis: [
        kpi('util_max', `Peak member utilisation (${names[je]})`, uMax, '-', uMax <= 0.87 ? 'ok' : uMax <= 1 ? 'warn' : 'bad', 'Ultimate load / (tension strength or Euler–Johnson column strength)'),
        kpi('mos_mount', 'Margin of safety of the mount', mos, '-', st3(mos)),
        kpi('force_max_N', 'Largest member force (limit)', sol.force[N.argmax(sol.force.map(Math.abs))], 'N'),
        kpi('stress_max_Pa', 'Largest member stress (limit)', Math.max(...stress.map(Math.abs)), 'Pa'),
        kpi('reaction_max_N', 'Largest firewall attachment reaction', Rmax, 'N'),
        kpi('ring_disp_m', 'Engine-ring displacement', disp, 'm'),
        kpi('mount_mass_kg', 'Tube mass of the mount', mass, 'kg'),
        kpi('equilibrium_err', 'Reaction vs applied load error', eq, '-', eq < 1e-8 ? 'ok' : 'warn'),
      ],
      plots: [
        { type: 'bar', title: 'Member axial force at limit load (tension positive)', ylabel: 'Force [kN]', categories: names, series: [{ name: 'Axial force', y: sol.force.map((F) => F / 1e3) }] },
        { type: 'bar', title: 'Member utilisation at ultimate load', ylabel: 'Utilisation [-]', categories: names, series: [{ name: 'Utilisation', y: util }] },
        { type: 'line', title: `Side view of the mount (deflection ×${sc.toPrecision(2)})`, xlabel: 'x forward of firewall [m]', ylabel: 'z up [m]', equalAspect: true, series: [{ name: 'Undeformed', ...side(false), style: 'dash' }, { name: 'Deformed', ...side(true) }] },
      ],
      tables: [{ title: 'Members', columns: ['Member', 'Length [m]', 'Force [kN]', 'Stress [MPa]', 'Euler load [kN]', 'Utilisation'], rows: members.map((_, e) => [names[e], sol.geo[e].L, sol.force[e] / 1e3, stress[e] / 1e6, Pcr[e] / 1e3, util[e]]) }],
      warnings,
      models: ['3-D pin-jointed direct-stiffness truss (13 members, 12 free DOF)', 'Euler–Johnson column strength for compression members'],
      assumptions: ['Engine treated as a rigid body loading the four ring nodes', 'Pin joints: welded-joint bending and fatigue at clusters are not modelled', 'Firewall attachments rigid', 'Vertical, side, thrust and torque loads applied together (conservative combination)', 'No fitting factor is applied: fittings not proven by test need at least 1.15 (14 CFR 25.625(a))', 'Default engine mass, mount geometry and side load factor are class-level estimates', matNote(m)],
    };
  },
  verify() {
    // symmetric tripod: three bars from a ring of radius r at z = 0 to an apex at height h under a vertical load P
    const r = 1, h = 2, P = 1000, E = 200e9, A = 1e-4, nodes = [...[0, 1, 2].map((k) => [r * Math.cos((2 * Math.PI * k) / 3), r * Math.sin((2 * Math.PI * k) / 3), 0]), [0, 0, h]];
    const s = truss3d(nodes, [[0, 3, A], [1, 3, A], [2, 3, A]], [0, 1, 2], [[3, 0, 0, -P]], E), Lb = Math.hypot(r, h);
    const km = N.kv(frame.run({ m_eng_kg: 150, n_z: 3.8, n_y: 1.33, T_N: 2600, Q_Nm: 420, x_cg: 0.25, L_m: 0.45, a_m: 0.6, c_m: 0.4, D_m: 0.019, tw_m: 0.0012, material: 'Steel 4340 (QT)', sf: 1.5 }));
    return [
      N.check('Tripod member force −P·L/(3h)', s.force[0], (-P * Lb) / (3 * h), 1e-10, 'Method of joints'),
      N.check('Tripod apex deflection P·L³/(3·E·A·h²)', -s.u[11], (P * Lb ** 3) / (3 * E * A * h * h), 1e-10, 'Unit-load method'),
      N.check('Mount truss: reactions balance applied loads', km.equilibrium_err, 0, 1e-9, 'Statics'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.util_max > 1) out.push({ severity: 'critical', title: 'Mount member overloaded', detail: `Utilisation ${o.util_max.toFixed(2)} at ultimate load.`, action: 'Increase tube diameter first (column strength rises with D²·area), then wall thickness; shortening the mount or widening the firewall spacing lowers member forces.', basis: 'Euler–Johnson column and tension strength at ultimate load' });
    else if (o.util_max < 0.4) out.push({ severity: 'advise', title: 'Mount is lightly loaded', detail: `Peak utilisation ${o.util_max.toFixed(2)}; tube mass ${o.mount_mass_kg.toFixed(2)} kg.`, action: 'Smaller tubes save mass, but check vibration isolation, fatigue at welds and minimum weldable gauge first.', basis: 'Utilisation' });
    out.push({ severity: 'info', title: 'Check dynamic and failure cases', detail: 'Static inertia, thrust and torque only.', action: 'Add gyroscopic, propeller-loss / blade-off and crash-inertia cases, and run engine-mount modes in Suite 10.', basis: 'Engine-mount load cases beyond the static envelope' });
    return out;
  },
};

// ---- geometrically nonlinear beam (corotational) --------------------------------------------
/**
 * 2-D corotational Euler–Bernoulli beam (Crisfield). Cantilever along x clamped at node 0, nodal dead loads Fref
 * ([Fx, Fy, M] per node) scaled in `steps` increments up to lamMax with full Newton iteration.
 */
function corotCantilever(L, EI, EA, Fref, steps, lamMax) {
  const n = EI.length, l0 = L / n, nd = 3 * (n + 1), u = new Array(nd).fill(0), alpha = new Array(n).fill(0), hist = [{ lam: 0, ux: 0, uy: 0, th: 0, it: 0 }];
  let maxIt = 0, ok = true;
  const assemble = (withK) => {
    const A = withK ? N.range(nd, () => new Float64Array(6)) : null, fi = new Array(nd).fill(0);
    for (let e = 0; e < n; e++) {
      const d = 3 * e, du = u[d + 3] - u[d], dx = l0 + du, dy = u[d + 4] - u[d + 1], ln = Math.hypot(dx, dy), ext = (2 * l0 * du + du * du + dy * dy) / (ln + l0); // round-off-safe ln − l0
      let al = Math.atan2(dy, dx); al += 2 * Math.PI * Math.round((alpha[e] - al) / (2 * Math.PI)); alpha[e] = al;
      const c = dx / ln, s = dy / ln, t1 = u[d + 2] - al, t2 = u[d + 5] - al, Nf = (EA * ext) / l0, k = EI[e] / l0, M1 = k * (4 * t1 + 2 * t2), M2 = k * (2 * t1 + 4 * t2);
      const r = [-c, -s, 0, c, s, 0], z = [s, -c, 0, -s, c, 0], b2 = z.map((v, a) => (a === 2 ? 1 : 0) - v / ln), b3 = z.map((v, a) => (a === 5 ? 1 : 0) - v / ln);
      for (let a = 0; a < 6; a++) fi[d + a] += r[a] * Nf + b2[a] * M1 + b3[a] * M2;
      if (withK) for (let a = 0; a < 6; a++) for (let b = a; b < 6; b++) A[d + a][b - a] += (EA / l0) * r[a] * r[b] + k * (4 * b2[a] * b2[b] + 2 * b2[a] * b3[b] + 2 * b3[a] * b2[b] + 4 * b3[a] * b3[b]) + (Nf / ln) * z[a] * z[b] + ((M1 + M2) / (ln * ln)) * (r[a] * z[b] + z[a] * r[b]);
    }
    return { A, fi };
  };
  const fn = N.norm(Fref) || 1;
  for (let st = 1; st <= steps && ok; st++) {
    const lam = (lamMax * st) / steps; let it = 0, conv = false;
    for (; it < 40; it++) {
      const { A, fi } = assemble(true), R = Fref.map((f, k) => lam * f - fi[k]);
      for (const d of [0, 1, 2]) fixDof(A, R, d, 5);
      if (N.norm(R) < 1e-8 * lam * fn + 1e-14) { conv = true; break; }
      let du; try { du = bandSolve(A, R, 5); } catch { break; }
      for (let k = 0; k < nd; k++) u[k] += du[k];
      if (N.norm(du) < 1e-12 * L) { conv = true; break; }
    }
    maxIt = Math.max(maxIt, it); ok = conv;
    if (conv) hist.push({ lam, ux: u[nd - 3], uy: u[nd - 2], th: u[nd - 1], it });
  }
  return { hist, u, converged: ok, maxIt, x: N.range(n + 1, (j) => j * l0 + u[3 * j]), y: N.range(n + 1, (j) => u[3 * j + 1]) };
}
const nonlinear = {
  id: 'nonlinear', title: 'Large-deflection (geometrically nonlinear) cantilever', fidelity: 'numerical',
  summary: 'Corotational beam elements with incremental Newton iteration trace the load–deflection curve of a flexible wing, blade or boom beyond the range where linear theory holds, and show how much linear analysis overstates deflection.',
  equations: ['Nonlinear elasticity equations', 'Green–Lagrange strain relations', 'Geometric stiffness equations', 'Principle of virtual work', 'Finite element weak-form equilibrium equations'],
  applicable: (c) => (c.meta.type === 'helicopter' && !(c.wing.S_m2 > 0) ? 'Rotor blades are held by centrifugal stiffening, which the beam analysis includes; this non-rotating large-deflection model does not apply.' : true),
  inputs: [
    { key: 'L', label: 'Length', unit: 'm', default: 17, min: 0.05, group: 'Geometry' },
    { key: 'EI_root', label: 'Root bending stiffness EI', unit: 'N·m²', default: 4e8, min: 1e-3, group: 'Stiffness', help: 'From the beam analysis when available' },
    { key: 'EI_ratio', label: 'Tip / root EI', unit: '-', default: 0.02, min: 1e-4, max: 1, group: 'Stiffness', help: 'Stiffness varies exponentially between root and tip' },
    { key: 'EA_root', label: 'Root axial stiffness EA', unit: 'N', default: 5e9, min: 1, group: 'Stiffness' },
    { key: 'F_dist_N', label: 'Total distributed transverse load (reference)', unit: 'N', default: 600000, group: 'Loads', help: 'Net limit load on the member, spread like an elliptic lift distribution' },
    { key: 'P_tip_N', label: 'Tip force (reference)', unit: 'N', default: 0, group: 'Loads' },
    { key: 'M_tip_Nm', label: 'Tip moment (reference)', unit: 'N·m', default: 0, group: 'Loads' },
    { key: 'lam_max', label: 'Maximum load multiplier', unit: '-', default: 1.5, min: 0.05, max: 20, group: 'Loads', help: '1.0 = reference (limit) load, 1.5 = ultimate' },
    { key: 'steps', label: 'Load increments', unit: '', default: 15, min: 2, max: 200, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nElem', label: 'Elements', unit: '', default: 24, min: 2, max: 160, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults(c, up, d) {
    const bd = beam.defaults(c, up, d), full = { ...Object.fromEntries(beam.inputs.map((f) => [f.key, f.default])), ...bd }, m = mat(full.material);
    const sec = (ch) => boxSection(full.box_chord_frac * ch, full.box_height_frac * full.tc * ch, full.t_skin_mm / 1e3, full.t_web_mm / 1e3, full.k_str), s0 = sec(full.c_root), s1 = sec(full.c_tip);
    const EI0 = up.fea?.EI_root_Nm2 ?? m.E * s0.I, F = full.n_load * full.L_1g_N * 0.75;
    return { L: full.L / Math.cos(N.rad(full.sweep_deg)), EI_root: EI0, EI_ratio: N.clamp(s1.I / s0.I, 1e-4, 1), EA_root: m.E * s0.A * (EI0 / (m.E * s0.I)), F_dist_N: F, P_tip_N: full.n_load * (full.P_tip_1g_N - G0 * full.m_tip_kg) };
  },
  run(i) {
    const n = Math.max(2, Math.round(i.nElem)), steps = Math.max(2, Math.round(i.steps)), EI = N.range(n, (e) => i.EI_root * i.EI_ratio ** ((e + 0.5) / n)), warnings = [];
    const xs = N.linspace(0, i.L, n + 1), sh = xs.map((x) => Math.sqrt(Math.max(0, 1 - (x / i.L) ** 2))), wts = sh.map((v, j) => v * (j === 0 || j === n ? 0.5 : 1)), ws = N.sum(wts) || 1, F = new Array(3 * (n + 1)).fill(0);
    for (let j = 0; j <= n; j++) F[3 * j + 1] = (i.F_dist_N * wts[j]) / ws;
    F[3 * n + 1] += i.P_tip_N; F[3 * n + 2] += i.M_tip_Nm;
    const nl = corotCantilever(i.L, EI, i.EA_root, F, steps, i.lam_max), last = nl.hist[nl.hist.length - 1];
    // linear reference with the same mesh and lumped loads (Hermite elements)
    const lin = beamSolve(xs, EI, new Array(n).fill(Infinity), new Array(n + 1).fill(0), N.range(n + 1, (j) => F[3 * j + 1]));
    const linTip = lin.w[n] + (i.M_tip_Nm ? linMomentTip(xs, EI, i.M_tip_Nm) : 0);
    const ratio = last.lam && linTip ? last.uy / (linTip * last.lam) : NaN, lams = nl.hist.map((h) => h.lam);
    if (!nl.converged) warnings.push(`Newton iteration stopped converging at load multiplier ${last.lam.toFixed(2)}: use more load increments.`);
    if (Math.abs(last.th) > 1.2) warnings.push('Tip rotation exceeds about 70°: loads are treated as fixed in direction (dead loads), whereas aerodynamic lift follows the surface.');
    return {
      kpis: [
        kpi('tip_deflection_nl_m', 'Tip deflection (nonlinear)', last.uy, 'm'),
        kpi('tip_deflection_lin_m', 'Tip deflection (linear theory)', linTip * last.lam, 'm'),
        kpi('nl_to_lin_ratio', 'Nonlinear / linear deflection', ratio, '-', !Number.isFinite(ratio) ? undefined : Math.abs(1 - ratio) < 0.05 ? 'ok' : 'warn', Number.isFinite(ratio) ? 'Within 5% means linear analysis is adequate' : 'Undefined: no transverse load'),
        kpi('tip_shortening_m', 'Tip axial shortening', -last.ux, 'm'),
        kpi('tip_rotation_deg', 'Tip rotation', N.deg(last.th), 'deg'),
        kpi('deflection_over_L', 'Tip deflection / length', last.uy / i.L, '-'),
        kpi('lam_reached', 'Load multiplier reached', last.lam, '-', nl.converged ? 'ok' : 'bad'),
        kpi('newton_iters_max', 'Most Newton iterations in an increment', nl.maxIt, '-'),
      ],
      plots: [
        { type: 'line', title: 'Load–deflection curve', xlabel: 'Tip deflection [m]', ylabel: 'Load multiplier [-]', series: [{ name: 'Nonlinear (corotational)', x: nl.hist.map((h) => h.uy), y: lams, style: 'line+points' }, { name: 'Linear theory', x: lams.map((l) => l * linTip), y: lams, style: 'dash' }] },
        { type: 'line', title: 'Deformed shape at the final load', xlabel: 'x [m]', ylabel: 'y [m]', equalAspect: true, series: [{ name: 'Nonlinear', x: nl.x, y: nl.y }, { name: 'Linear', x: xs, y: lin.w.map((v) => v * last.lam), style: 'dash' }, { name: 'Unloaded', x: [0, i.L], y: [0, 0], style: 'dash' }] },
        { type: 'line', title: 'Tip shortening', xlabel: 'Load multiplier [-]', ylabel: 'Axial shortening [m]', series: [{ name: 'Shortening', x: lams, y: nl.hist.map((h) => -h.ux) }] },
      ],
      warnings,
      models: [`Corotational 2-D Euler–Bernoulli beam elements (${n})`, `Incremental load control with full Newton–Raphson (${steps} increments)`],
      assumptions: ['Linear elastic material, large displacements and rotations, small strains', 'Dead (fixed-direction) loads; follower aerodynamic loading is not modelled', 'Planar bending without torsion', 'Stiffness varies exponentially from root to tip'],
    };
  },
  convergence: { param: 'nElem', label: 'Elements', levels: [4, 8, 16, 32, 64], metric: 'tip_deflection_nl_m' },
  verify() {
    const L = 1, EI = 1, n = 40, EA = 1e8, tipLoad = (al) => { const F = new Array(3 * (n + 1)).fill(0); F[3 * n + 1] = (al * EI) / (L * L); return corotCantilever(L, new Array(n).fill(EI), EA, F, 10, 1).hist.pop(); };
    // independent elastica reference by shooting on θ'' = −α cos θ, θ(0) = 0, θ'(1) = 0
    const elastica = (al) => {
      const shoot = (k0) => N.rk4((s, y) => [y[1], -al * Math.cos(y[0]), Math.cos(y[0]), Math.sin(y[0])], 0, [0, k0, 0, 0], 1, 400).y.pop();
      const k0 = N.brent((k) => shoot(k)[1], 0, al, 1e-13), y = shoot(k0); return { w: y[3], u: 1 - y[2] };
    };
    const a1 = tipLoad(1), a5 = tipLoad(5), e5 = elastica(5), F = new Array(3 * (n + 1)).fill(0); F[3 * n + 2] = (Math.PI / 2) * EI / L;
    const mq = corotCantilever(L, new Array(n).fill(EI), EA, F, 6, 1).hist.pop();
    const small = N.kv(nonlinear.run({ L: 2, EI_root: 1000, EI_ratio: 1, EA_root: 1e9, F_dist_N: 0, P_tip_N: 1, M_tip_Nm: 0, lam_max: 1, steps: 2, nElem: 10 }));
    return [
      N.check('Elastica PL²/EI = 1: tip deflection 0.3017 L', a1.uy, 0.30172, 1e-3, 'Bisshopp & Drucker (1945); Mattiasson (1981)'),
      N.check('Elastica PL²/EI = 1: tip shortening 0.0564 L', -a1.ux, 0.05643, 3e-3, 'Mattiasson (1981)'),
      N.check('Elastica PL²/EI = 5 vs independent shooting solution (deflection)', a5.uy, e5.w, 1e-3, 'Elastica ODE integrated by RK4 with Brent shooting'),
      N.check('Elastica PL²/EI = 5 vs independent shooting solution (shortening)', -a5.ux, e5.u, 2e-3, 'Elastica ODE'),
      N.check('End moment ML/EI = π/2 bends the beam into a quarter circle', mq.uy, 2 / Math.PI, 1e-3, 'Pure bending: constant curvature M/EI'),
      N.check('Small-load limit recovers PL³/3EI', small.tip_deflection_nl_m, 8 / 3000, 1e-5, 'Linear beam theory'),
    ];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (!Number.isFinite(o.nl_to_lin_ratio)) return [{ severity: 'info', title: 'No transverse load', detail: 'The linear deflection is zero, so the nonlinear / linear ratio is undefined.', action: 'Enter a distributed load, tip force or tip moment.', basis: 'Input check' }];
    if (Math.abs(1 - o.nl_to_lin_ratio) > 0.05) out.push({ severity: 'advise', title: 'Geometric nonlinearity matters', detail: `Linear theory misses the true deflection by ${(100 * Math.abs(1 - o.nl_to_lin_ratio)).toFixed(0)}% (deflection is ${(100 * o.deflection_over_L).toFixed(0)}% of the length).`, action: 'Use nonlinear static aeroelastic loads for this configuration: lift tilts inboard as the wing bends, changing root bending moment, trim and flutter (Suite 3).', basis: 'Elastica theory; linear beam theory is accurate below about 10% deflection' });
    else out.push({ severity: 'info', title: 'Linear analysis is adequate', detail: `Nonlinear and linear deflections agree within ${(100 * Math.abs(1 - o.nl_to_lin_ratio)).toFixed(1)}%.`, action: 'Linear loads and aeroelastic analysis can be used with confidence at this load level.', basis: 'Nonlinear / linear deflection ratio' });
    return out;
  },
};
/** Linear tip deflection of a cantilever with piecewise-constant EI under a tip moment (second moment–area theorem). */
function linMomentTip(xs, EI, M) { let w = 0; const L = xs[xs.length - 1]; for (let e = 0; e < EI.length; e++) w += (M / EI[e]) * (xs[e + 1] - xs[e]) * (L - 0.5 * (xs[e] + xs[e + 1])); return w; }

// ---- torsion of thin-walled closed sections -------------------------------------------------
/** Multi-cell Bredt–Batho torsion. cells: [{A, d (∮ds/t of the whole cell)}], shared[k] = ds/t of the wall between cell k and k+1. */
function multicell(cells, shared, T, G) {
  const n = cells.length, Am = N.zeros(n + 1), b = new Array(n + 1).fill(0);
  for (let k = 0; k < n; k++) { Am[k][k] = cells[k].d; if (k > 0) Am[k][k - 1] = -shared[k - 1]; if (k < n - 1) Am[k][k + 1] = -shared[k]; Am[k][n] = -2 * cells[k].A * G; Am[n][k] = 2 * cells[k].A; }
  b[n] = 1; const x1 = N.solve(Am, b); // unit torque: q and twist rate scale linearly with T
  return { q: x1.slice(0, n).map((v) => v * T), rate: x1[n] * T, GJ: 1 / x1[n] };
}
const torsion = {
  id: 'torsion', title: 'Torsion of the closed wing section (Bredt–Batho)', fidelity: 'analytical',
  summary: 'Shear flows, shear stress, twist rate and torsional stiffness of a one- or two-cell thin-walled section (nose cell plus main box), compared with the same section cut open.',
  equations: ['Saint-Venant torsion equations', 'Static equilibrium equations', 'Hooke’s law'],
  applicable: (c) => (c.wing.S_m2 > 0 || (c.meta.type === 'helicopter' && c.rotor.chord_m > 0) ? true : 'No wing or blade section is defined for this vehicle.'),
  inputs: [
    { key: 'chord', label: 'Chord at the station', unit: 'm', default: 5.8, min: 0.01, group: 'Geometry' },
    { key: 'tc', label: 'Thickness-to-chord ratio', unit: '-', default: 0.115, min: 0.02, max: 0.5, group: 'Geometry' },
    { key: 'x_fs', label: 'Front spar position / chord', unit: '-', default: 0.15, min: 0.02, max: 0.5, group: 'Geometry' },
    { key: 'x_rs', label: 'Rear spar position / chord', unit: '-', default: 0.6, min: 0.2, max: 0.95, group: 'Geometry' },
    { key: 'h_frac', label: 'Box height / section thickness', unit: '-', default: 0.85, min: 0.3, max: 1, group: 'Geometry' },
    { key: 'cells', label: 'Cells carrying torque', type: 'select', options: ['box only', 'nose cell + box'], default: 'nose cell + box', group: 'Geometry' },
    { key: 't_skin', label: 'Box skin thickness', unit: 'm', default: 0.006, min: 1e-4, group: 'Thickness' },
    { key: 't_nose', label: 'Nose skin thickness', unit: 'm', default: 0.002, min: 1e-4, group: 'Thickness' },
    { key: 't_fs', label: 'Front spar web thickness', unit: 'm', default: 0.005, min: 1e-4, group: 'Thickness' },
    { key: 't_rs', label: 'Rear spar web thickness', unit: 'm', default: 0.005, min: 1e-4, group: 'Thickness' },
    { key: 'T_Nm', label: 'Applied torque (limit)', unit: 'N·m', default: 300000, group: 'Loads', help: 'Root torque from the beam analysis when available' },
    { key: 'sf', label: 'Ultimate safety factor', unit: '-', default: 1.5, min: 1, max: 3, group: 'Loads' },
    MAT,
  ],
  defaults(c, up, d) {
    const s = c.struct, heli = !(c.wing.S_m2 > 0), ch = heli ? c.rotor.chord_m : d.c_root, bf = heli ? 0.35 : s.box_chord_frac, tw = heli ? s.t_spar_mm : Math.min(s.t_spar_mm, Math.max(0.5, ch));
    const ts = up.fea?.t_skin_root_m ?? (heli ? s.t_skin_mm : Math.min(s.t_skin_mm, Math.max(0.4, 0.6 * ch))) / 1e3, qd = 0.5 * RHO0 * (1.25 * (c.aero.Vmo_ms || c.flight.V_ms)) ** 2;
    return { chord: ch, tc: heli ? 0.12 : c.wing.tc, x_fs: 0.15, x_rs: 0.15 + bf, h_frac: s.box_height_frac, t_skin: ts, t_nose: Math.max(0.0004, 0.5 * ts), t_fs: tw / 1e3, t_rs: tw / 1e3,
      T_Nm: Math.abs(up.fea?.root_torque_Nm || 0) || (heli ? 0.02 * d.W * ch : Math.abs(c.aero.Cm0) * qd * c.wing.S_m2 * d.mac / 2 + 0.15 * d.mac * c.aero.n_pos * d.W / 2), sf: s.sf_ultimate, material: s.material };
  },
  run(i) {
    const m = mat(i.material), c = i.chord, h = i.h_frac * i.tc * c, w = Math.max(1e-6, (i.x_rs - i.x_fs) * c), an = i.x_fs * c, bn = h / 2, warnings = [];
    // nose cell idealised as a half ellipse (semi-axes x_fs·c and h/2); perimeter by Ramanujan's approximation
    const pn = 0.5 * Math.PI * (3 * (an + bn) - Math.sqrt((3 * an + bn) * (an + 3 * bn))), box = { A: w * h, d: (2 * w) / i.t_skin + h / i.t_fs + h / i.t_rs }, nose = { A: 0.5 * Math.PI * an * bn, d: pn / i.t_nose + h / i.t_fs };
    const two = i.cells !== 'box only', r = two ? multicell([nose, box], [h / i.t_fs], i.T_Nm, m.G) : multicell([box], [], i.T_Nm, m.G), qB = r.q[two ? 1 : 0], qN = two ? r.q[0] : 0;
    const walls = [['Box skins', qB, i.t_skin], ['Rear spar web', qB, i.t_rs], ['Front spar web', qB - qN, i.t_fs], ...(two ? [['Nose skin', qN, i.t_nose]] : [])], tau = walls.map(([, q, t]) => Math.abs(q) / t), jm = N.argmax(tau);
    const al = designAllowables(m), tauAllow = al.Su / Math.sqrt(3), mos = tau[jm] > 0 ? tauAllow / (i.sf * tau[jm]) - 1 : Infinity;
    const Jopen = ((2 * w * i.t_skin ** 3 + h * i.t_fs ** 3 + (two ? pn * i.t_nose ** 3 : 0)) / 3), GJopen = m.G * Jopen, single = (m.G * 4 * box.A ** 2) / box.d;
    // skin shear buckling screen between spars is left to the buckling analysis; flag very thin walls
    if (Math.min(i.t_skin, i.t_nose) / c < 2e-4) warnings.push('Very thin skins relative to chord: shear buckling will limit the usable shear flow (see the buckling analysis).');
    if (i.x_rs <= i.x_fs) warnings.push('Rear spar is not aft of the front spar; the box width was clipped.');
    // outline for the section sketch
    const ph = N.linspace(Math.PI / 2, 1.5 * Math.PI, 31), ox = [...ph.map((p) => an + an * Math.cos(p)), i.x_rs * c, i.x_rs * c, an], oy = [...ph.map((p) => bn * Math.sin(p)), -bn, bn, bn];
    return {
      kpis: [
        kpi('GJ_Nm2', 'Torsional stiffness GJ', r.GJ, 'N·m²'),
        kpi('twist_rate_deg_m', 'Twist rate at limit torque', N.deg(r.rate), 'deg/m'),
        kpi('q_box_Npm', 'Shear flow in the main box', qB, 'N/m'),
        kpi('q_nose_Npm', 'Shear flow in the nose cell', qN, 'N/m'),
        kpi('tau_max_Pa', `Peak shear stress (${walls[jm][0]})`, tau[jm], 'Pa'),
        kpi('mos_shear', 'Margin of safety in shear (ultimate)', mos, '-', st3(mos), `Allowable Su/√3 with Su = ${(al.Su / 1e6).toFixed(0)} MPa (${al.design ? 'design allowable' : 'typical'}); a measured shear allowable Fsu should replace it`),
        kpi('GJ_box_only_Nm2', 'GJ of the main box alone', single, 'N·m²'),
        kpi('nose_torque_share', 'Share of torque carried by the nose cell', two ? (2 * nose.A * qN) / (i.T_Nm || 1) : 0, '-'),
        kpi('GJ_open_Nm2', 'GJ if the section were cut open', GJopen, 'N·m²', undefined, 'Σ b·t³/3: shows why closing the cell matters'),
        kpi('closed_to_open', 'Closed / open stiffness ratio', r.GJ / GJopen, '-'),
      ],
      plots: [
        { type: 'line', title: 'Idealised section', xlabel: 'Chordwise position [m]', ylabel: 'Height [m]', equalAspect: true, series: [{ name: two ? 'Nose cell and box' : 'Box (nose not structural)', x: ox, y: oy }, { name: 'Front spar', x: [an, an], y: [-bn, bn] }] },
        { type: 'bar', title: 'Shear stress by wall at limit torque', ylabel: 'Shear stress [MPa]', categories: walls.map((v) => v[0]), series: [{ name: 'Shear stress', y: tau.map((v) => v / 1e6) }] },
      ],
      tables: [{ title: 'Shear flows', columns: ['Wall', 'Shear flow [kN/m]', 'Thickness [mm]', 'Shear stress [MPa]'], rows: walls.map(([nm, q, t], k) => [nm, q / 1e3, t * 1e3, tau[k] / 1e6]) }],
      warnings,
      models: [two ? 'Two-cell Bredt–Batho torsion with compatibility of twist' : 'Single-cell Bredt–Batho torsion', 'Open-section Saint-Venant torsion constant Σbt³/3'],
      assumptions: ['Thin walls with constant shear flow per wall; free warping (no root restraint stiffening)', 'Nose cell idealised as a half ellipse; trailing-edge structure ignored', 'Default spar positions and nose-skin gauge are typical values', 'Skins assumed not to buckle in shear', matNote(m)],
    };
  },
  verify() {
    const G = 27e9, a = 0.2, t = 0.002, T = 1000, sqr = { A: a * a, d: (4 * a) / t }, two = multicell([sqr, sqr], [a / t], T, G), one = multicell([{ A: 2 * a * a, d: (6 * a) / t }], [], T, G);
    const R = 0.1, tube = multicell([{ A: Math.PI * R * R, d: (2 * Math.PI * R) / t }], [], T, G);
    return [
      N.check('Two equal cells: no net flow in the shared wall, GJ equals the outer cell', two.GJ, one.GJ, 1e-12, 'Symmetry of the Bredt–Batho compatibility equations'),
      N.check('Two equal cells: equal shear flows T/(4A)', two.q[0], T / (4 * a * a), 1e-12, 'Bredt–Batho'),
      N.check('Thin circular tube GJ = 2πR³tG', tube.GJ, 2 * Math.PI * R ** 3 * t * G, 1e-12, 'Saint-Venant torsion, thin-wall limit'),
      N.check('Thin circular tube shear flow T/(2A)', tube.q[0], T / (2 * Math.PI * R * R), 1e-12, 'Bredt–Batho'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.mos_shear < 0) out.push({ severity: 'critical', title: 'Shear stress exceeds the allowable under ultimate torque', detail: `Margin ${o.mos_shear.toFixed(2)}.`, action: 'Thicken the governing wall or enlarge the enclosed area: shear flow falls with 1/(2A).', basis: 'τ = T/(2·A·t) against Su/√3' });
    if (i.cells === 'box only') out.push({ severity: 'info', title: 'Nose cell not counted', detail: 'Only the main box carries torque in this run.', action: 'If the leading edge is structural and continuous, the two-cell option adds stiffness at no mass cost — useful for flutter and aileron reversal margins.', basis: 'Bredt–Batho: GJ ∝ A²' });
    out.push({ severity: 'info', title: 'Keep the cell closed', detail: `Cutting the section open would divide GJ by ${o.closed_to_open.toExponential(1)}.`, action: 'Access panels and cut-outs in the box skins must be load-carrying (stressed doors) or reinforced with a torsion frame.', basis: 'Open- vs closed-section torsion' });
    return out;
  },
};

export default {
  id: 'fea', n: 2,
  tagline: 'Is the structure strong, stiff and stable enough — and how much does it weigh?',
  analyses: [beam, plane, buckling, fuselage, frame, nonlinear, torsion],
  consumes: [{ from: 'performance', keys: ['n_limit', 'V_d_eas'], why: 'Limit load factor and dive speed for the design load case' }],
  provides: [
    { key: 'sigma_max_Pa', label: 'Peak stress at limit load', unit: 'Pa' }, { key: 'tip_deflection_m', label: 'Tip deflection', unit: 'm' }, { key: 'margin_of_safety', label: 'Margin of safety', unit: '-' },
    { key: 'buckling_factor', label: 'Buckling factor', unit: '-' }, { key: 'EI_root_Nm2', label: 'Root bending stiffness', unit: 'N·m²' }, { key: 'GJ_root_Nm2', label: 'Root torsional stiffness', unit: 'N·m²' },
    { key: 'wing_struct_mass_kg', label: 'Wing structural mass', unit: 'kg' },
  ],
  handoff: [
    { model: 'General 3-D solid and shell finite elements on imported CAD/meshes', why: 'Needs a geometry kernel, unstructured meshing and sparse solvers for 10⁵–10⁷ DOF, beyond an in-browser dense/banded kernel', tool: 'General-purpose FE solver (Nastran / Abaqus / ANSYS / CalculiX / Code_Aster class)' },
    { model: 'Contact mechanics, fastener and joint load transfer', why: 'Requires contact search and nonlinear constraint enforcement on real joint geometry', tool: 'Nonlinear implicit FE with contact; fastener-flexibility joint tools' },
    { model: 'Elastic–plastic, hyperelastic and damage constitutive FE', why: 'Only elastic stress with yield-based margins is solved here; path-dependent material integration needs calibrated hardening data', tool: 'Nonlinear implicit/explicit FE solver' },
    { model: 'Post-buckling, von Kármán large-deflection plates and imperfection-sensitive shells', why: 'Requires arc-length continuation on shell meshes with measured imperfections; only bifurcation loads and effective width are computed', tool: 'Nonlinear shell FE with Riks/arc-length control' },
    { model: 'Mindlin–Reissner plate/shell and composite shell models', why: 'Laminate response is covered at ply level in Suite 21; shell FE of built-up structure is not solved natively', tool: 'Shell FE solver with composite property cards' },
    { model: 'Geometrically exact 3-D beams with bending–torsion coupling, multibody–FE coupling', why: 'The corotational beam is planar; coupled 3-D rotor and swept-wing dynamics are external', tool: 'Multibody / comprehensive rotorcraft code, nonlinear aeroelastic beam solver' },
    { model: 'Structural–thermal coupling and multiscale homogenisation', why: 'Thermal fields come from Suite 12 but thermo-elastic FE and RVE homogenisation are not solved natively', tool: 'Coupled-field FE; micromechanics codes' },
    { model: 'Topology-optimised and lattice structures, global aircraft load-path FE', why: 'Needs full-airframe models and optimisation on large meshes', tool: 'Global FE model (GFEM) with structural optimisation' },
  ],
};
