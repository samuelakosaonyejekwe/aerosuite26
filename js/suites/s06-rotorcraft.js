// Suite 6 — Rotorcraft Aerodynamics and Aeromechanics.
// Hover blade-element momentum theory, forward-flight trim with rigid-blade flapping, flap/lag time response and
// hub-load harmonics, Pitt–Peters dynamic inflow, prescribed tip-vortex wake with Biot–Savart induction and
// blade–vortex miss distance, a Leishman–Beddoes-type dynamic-stall demonstrator, and autorotation.
// Conventions: azimuth ψ = 0 over the tail, 90° on the advancing side; flapping β = a0 − a1·cosψ − b1·sinψ.

import * as N from '../core/numerics.js';
import { isa, G0 } from '../core/atmosphere.js';

// ---- shared physics -------------------------------------------------------------------------
const rotorOnly = (c) => (c.rotor.R_m > 0 && c.rotor.rpm > 0 ? true : 'This suite needs a lifting rotor (radius and speed in the case); the current vehicle is a fixed-wing aeroplane. Use Suite 8 for its propeller.');
const nRot = (c) => (c.meta.type === 'helicopter' ? 1 : Math.max(1, c.prop.n_eng));
const TWO_PI = 2 * Math.PI;

/**
 * Sectional loads on a blade element; uT, uP are the tangential and perpendicular velocities over tip speed.
 * Returns a shared scratch array [fz, fx, alpha]: fz = u²(cl·cosφ − cd·sinφ) (thrust direction), fx = u²(cl·sinφ + cd·cosφ) (against rotation).
 * p.linear selects the small-angle, stall-free textbook form used by the closed-form verification cases.
 */
function section(p, th, uT, uP) {
  if (p.linear) { SEC[0] = p.cla * (th * uT * uT - uP * uT); SEC[1] = p.cla * (th * uT * uP - uP * uP) + p.cd0 * uT * uT; SEC[2] = th - uP / (uT || 1e-9); return SEC; }
  const rev = uT < 0, u2 = uT * uT + uP * uP, phi = Math.atan2(uP, rev ? -uT : uT), al = rev ? -th - phi : th - phi; // reverse-flow region: flow meets the trailing edge
  const M = Math.sqrt(u2) * p.Mtip, Mc = M < 0.8 ? M : 0.8, a = p.cla / Math.sqrt(1 - Mc * Mc), as = p.clmax / a, aa = al < 0 ? -al : al, ab = aa < as ? aa : as, dM = M - p.Mdd;
  const cl = aa < as ? a * al : al < 0 ? -p.clmax : p.clmax;
  const cd = p.cd0 + p.cd2 * ab * ab + (aa > as ? Math.min(1.8, 1.5 * (aa - as)) : 0) + (dM > 0 ? 20 * dM * dM * dM * dM : 0);
  const s = Math.sin(phi), c = Math.cos(phi);
  SEC[0] = u2 * (cl * c - cd * s); SEC[1] = (rev ? -u2 : u2) * (cl * s + cd * c); SEC[2] = al;
  return SEC;
}
const SEC = [0, 0, 0]; // scratch result of section(): read immediately, never stored
/** Rotor/section parameter set from analysis inputs. */
function rotorOf(i, extra = {}) {
  const at = isa(i.alt_m ?? 0, i.dISA ?? 0), Om = (i.rpm * TWO_PI) / 60, Vt = Om * i.R, A = Math.PI * i.R * i.R;
  return { R: i.R, Nb: Math.max(1, Math.round(i.Nb)), sigma: (Math.max(1, Math.round(i.Nb)) * i.chord) / (Math.PI * i.R), cla: i.cla, cd0: i.cd0, cd2: i.cd2 ?? 0.3, clmax: i.clmax ?? 1.3, tw: N.rad(i.twist_deg ?? 0), rc: i.root_cut ?? 0.15, Om, Vt, A, rho: at.rho, a: at.a, Mtip: Vt / at.a, Mdd: i.M_dd ?? 0.8, qA: at.rho * A * Vt * Vt, linear: false, ...extra };
}
/** Glauert induced inflow ratio for thrust coefficient CT, advance ratio μ and normal free-stream component λc (positive down through the disk). */
function lamI(CT, mu, lc) {
  if (CT <= 0) return 0;
  const lh = Math.sqrt(CT / 2);
  if (mu < 1e-9 && Math.abs(lc) < 1e-9) return lh;
  const g = (l) => l - CT / (2 * Math.sqrt(mu * mu + (lc + l) ** 2));
  const r = N.findRoot(g, 0, 2 * lh + Math.abs(lc) + 0.05, 40, 1e-13);
  return Number.isFinite(r) ? r : lh;
}
/** Azimuth-integrated rotor loads in the tip-path plane for controls x = [θ75, θ1c, θ1s, a0]. */
function diskLoads(p, s, x, grid) {
  const [th0, c1, s1, a0] = x, nR = s.nR, nP = s.nPsi, dr = (1 - p.rc) / nR; let CT = 0, CQ = 0, CH = 0, M0 = 0, Mc = 0, Ms = 0;
  for (let k = 0; k < nP; k++) {
    const psi = (TWO_PI * k) / nP, sp = Math.sin(psi), cp = Math.cos(psi), th = th0 + c1 * cp + s1 * sp; let Fz = 0, Fx = 0, Q = 0, Mz = 0;
    for (let j = 0; j < nR; j++) {
      const r = p.rc + (j + 0.5) * dr, uT = r + s.mu * sp, uP = s.lam + s.lami * r * (s.kx * cp + s.ky * sp) + s.mu * a0 * cp;
      const q = section(p, th + p.tw * (r - 0.75), uT, uP), al = q[2]; let fz = q[0], fx = q[1];
      if (r > s.B) { fz = 0; fx = p.cd0 * uT * Math.abs(uT); } // Prandtl tip-loss as an effective radius: no lift outboard of B
      Fz += fz * dr; Fx += fx * dr; Q += fx * r * dr; Mz += fz * r * dr;
      if (grid) grid[j][k] = al;
    }
    CT += Fz; CQ += Q; CH += Fx * sp - Fz * a0 * cp; M0 += Mz; Mc += Mz * cp; Ms += Mz * sp;
  }
  const f = p.sigma / (2 * nP);
  return { CT: f * CT, CQ: f * CQ, CH: f * CH, M0: M0 / nP, M1c: (2 * Mc) / nP, M1s: (2 * Ms) / nP };
}
/**
 * Forward-flight trim of one rotor in the tip-path plane: finds θ75, cyclic and coning so that thrust balances weight and
 * drag and the first-harmonic flap moments vanish (centrally hinged rigid blade). o = {V, Wr, Dr, gam (climb angle), kappa,
 * lock, drees, tipLoss, nR, nPsi, passes}.
 */
function trimFF(p, o) {
  let H = 0, out = null, x = null;
  for (let pass = 0; pass < (o.passes ?? 2); pass++) {
    const aT = Math.atan2(o.Dr + H, o.Wr), T = Math.hypot(o.Wr, o.Dr + H), CT = T / p.qA, aS = aT + (o.gam || 0); // disk incidence to the flight path
    const mu = (o.V * Math.cos(aS)) / p.Vt, lc = (o.V * Math.sin(aS)) / p.Vt, li0 = lamI(CT, mu, lc), li = (o.kappa ?? 1) * li0, lam = lc + li;
    const chi = Math.atan2(mu, Math.abs(lam) + 1e-12), dre = o.drees && mu > 1e-3;
    const s = { mu, lam, lami: li, kx: dre ? ((4 / 3) * (1 - Math.cos(chi) - 1.8 * mu * mu)) / Math.sin(chi) : 0, ky: dre ? -2 * mu : 0, nR: o.nR, nPsi: o.nPsi, B: o.tipLoss ? 1 - Math.sqrt(2 * CT) / p.Nb : 1.01 };
    const res = (y) => { const L = diskLoads(p, s, y); return [(L.CT - CT) / CT, L.M1c, L.M1s, y[3] - (o.lock / (2 * p.cla)) * L.M0]; };
    const th = (6 * CT) / (p.sigma * p.cla) + 1.5 * lam, sol = N.fsolve(res, x || [th, 0, (-8 / 3) * mu * th + 2 * mu * lam, (o.lock / 8) * (th - (4 * lam) / 3)], { tol: 1e-8, maxIter: 25 });
    x = sol.x; const L = diskLoads(p, s, x); H = L.CH * p.qA;
    out = { x, s, L, CT, T, aT, aS, mu, lam, li, li0, converged: sol.residual < 1e-6, th75: x[0], th1c: x[1], th1s: x[2], a0: x[3], a1: -x[2], b1: x[1], P: L.CQ * p.qA * p.Vt, Q: L.CQ * p.qA * p.R, H };
  }
  out.aRet = out.th75 + p.tw * 0.25 - out.th1s - Math.atan2(out.lam - out.s.lami * out.s.ky, 1 - out.mu); // Gessow–Myers retreating-tip incidence α(1.0, 270°)
  return out;
}
/** Tail-rotor thrust and power to react main-rotor torque Q at speed V (edgewise Glauert inflow). */
function tailRotor(i, rho, Q, V) {
  if (!(i.tr_R > 0 && i.tr_arm > 0 && i.tr_rpm > 0) || i.n_rotors > 1) return { T: 0, P: 0 };
  const T = Q / i.tr_arm, A = Math.PI * i.tr_R ** 2, Vt = ((i.tr_rpm * TWO_PI) / 60) * i.tr_R, sig = (i.tr_Nb * i.tr_chord) / (Math.PI * i.tr_R), vh2 = T / (2 * rho * A);
  const vi = Math.sqrt((-V * V + Math.sqrt(V ** 4 + 4 * vh2 * vh2)) / 2);
  return { T, P: (i.kappa ?? 1.15) * T * vi + ((rho * A * Vt ** 3 * sig * i.cd0) / 8) * (1 + 4.65 * (V / Vt) ** 2) };
}

const MASS = [
  { key: 'mass_kg', label: 'Gross mass', unit: 'kg', default: 9980, min: 0.05, group: 'Aircraft' },
  { key: 'n_rotors', label: 'Lifting rotors', unit: '', default: 1, min: 1, max: 32, step: 1, discrete: true, group: 'Aircraft', help: '1 for a conventional helicopter, 4 for a quadrotor' },
];
const ROTOR = [
  { key: 'R', label: 'Rotor radius', unit: 'm', default: 8.18, min: 0.02, group: 'Rotor' },
  { key: 'Nb', label: 'Blades per rotor', unit: '', default: 4, min: 1, max: 12, step: 1, discrete: true, group: 'Rotor' },
  { key: 'chord', label: 'Mean blade chord', unit: 'm', default: 0.53, min: 0.002, group: 'Rotor' },
  { key: 'twist_deg', label: 'Linear twist, root to tip', unit: 'deg', default: -16, min: -40, max: 10, group: 'Rotor', help: 'Negative = washout. Typical −8 to −18' },
  { key: 'rpm', label: 'Rotor speed', unit: 'rpm', default: 258, min: 10, group: 'Rotor' },
  { key: 'cla', label: 'Section lift-curve slope', unit: '1/rad', default: 5.73, min: 3, max: 7, group: 'Blade section' },
  { key: 'cd0', label: 'Section profile drag cd0', unit: '-', default: 0.0095, min: 0, max: 0.08, group: 'Blade section', help: '0.008–0.011 full scale; 0.015–0.03 for small low-Reynolds rotors' },
  { key: 'clmax', label: 'Section maximum lift coefficient', unit: '-', default: 1.3, min: 0.5, max: 2, group: 'Blade section' },
];
const ATM = [
  { key: 'alt_m', label: 'Altitude', unit: 'm', default: 500, min: -500, max: 9000, group: 'Atmosphere' },
  { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -60, max: 50, group: 'Atmosphere' },
];
const TAIL = [
  { key: 'tr_R', label: 'Tail-rotor radius', unit: 'm', default: 1.68, min: 0, group: 'Tail rotor', help: '0 for no tail rotor (multirotor, coaxial, tandem)' },
  { key: 'tr_chord', label: 'Tail-rotor chord', unit: 'm', default: 0.247, min: 0, group: 'Tail rotor' },
  { key: 'tr_Nb', label: 'Tail-rotor blades', unit: '', default: 4, min: 0, step: 1, discrete: true, group: 'Tail rotor' },
  { key: 'tr_rpm', label: 'Tail-rotor speed', unit: 'rpm', default: 1190, min: 0, group: 'Tail rotor' },
  { key: 'tr_arm', label: 'Tail-rotor arm', unit: 'm', default: 9.93, min: 0, group: 'Tail rotor' },
];
const KAPPA = { key: 'kappa', label: 'Induced power factor κ', unit: '-', default: 1.15, min: 1, max: 1.6, group: 'Calibration', help: 'Non-uniform inflow, tip loss and wake effects on momentum theory; 1.10–1.20 typical' };
const baseDefaults = (c, up, d) => ({
  mass_kg: c.mass.mtow_kg, n_rotors: nRot(c), R: c.rotor.R_m || undefined, Nb: c.rotor.n_blades || undefined, chord: c.rotor.chord_m || undefined, twist_deg: c.rotor.twist_deg,
  rpm: c.rotor.rpm || undefined, cla: c.rotor.cla, cd0: c.rotor.cd0, alt_m: c.atm.alt_m, dISA: c.atm.dISA_K,
});
const tailDefaults = (c) => (c.meta.type === 'helicopter' && c.rotor.tr_R_m > 0 ? { tr_R: c.rotor.tr_R_m, tr_chord: c.rotor.tr_chord_m, tr_Nb: c.rotor.tr_n_blades, tr_rpm: c.rotor.tr_rpm, tr_arm: c.rotor.tr_arm_m } : { tr_R: 0, tr_chord: 0, tr_Nb: 0, tr_rpm: 0, tr_arm: 0 });
const ffSpeed = (c, d) => (d.isRotary ? c.mission.cruise_V_ms || c.flight.V_ms : Math.min(c.flight.V_ms, 0.15 * d.v_tip));
const wingShare = (c, d, V) => (c.wing.S_m2 > 0 ? Math.min(0.9, (0.5 * isa(c.atm.alt_m, c.atm.dISA_K).rho * V * V * c.wing.S_m2 * 0.7 * c.aero.CLmax_clean) / d.W) : 0);
const fmt = (v, n = 1) => (Number.isFinite(v) ? v.toFixed(n) : 'n/a');

// ---- (a) hover: blade-element momentum theory -----------------------------------------------
/** Hover BEMT with Prandtl tip loss. g = {Nb, sigma, tap, tw, ideal, cla, clmax, cd0, cd2, rc, tipLoss, nR}. */
function hoverBEMT(g, th75) {
  const n = g.nR, dr = (1 - g.rc) / n, r = [], lam = [], dT = [], al = [], Fs = []; let CT = 0, CTm = 0, CPi = 0, CP0 = 0, aMax = -9;
  for (let j = 0; j < n; j++) {
    const x = g.rc + (j + 0.5) * dr, th = g.ideal ? (th75 * 0.75) / x : th75 + g.tw * (x - 0.75), s = (g.sigma * 2 * (1 - (1 - g.tap) * x)) / (1 + g.tap);
    const cl = (l) => N.clamp(g.cla * (th - l / x), -g.clmax, g.clmax);
    const F = (l) => (g.tipLoss ? Math.max(1e-6, (2 / Math.PI) * Math.acos(Math.min(1, Math.exp((-g.Nb * (1 - x)) / (2 * Math.max(l, 1e-9)))))) : 1);
    let lo = 0, hi = 0.6, l = 0;
    if (cl(0) > 0) for (let k = 0; k < 60; k++) { l = 0.5 * (lo + hi); if (8 * F(l) * l * l - s * cl(l) * x > 0) hi = l; else lo = l; } // momentum = blade element at this annulus
    const a = th - l / x, as = g.clmax / g.cla, cd = g.cd0 + g.cd2 * Math.min(Math.abs(a), as) ** 2 + (Math.abs(a) > as ? Math.min(1.8, 1.5 * (Math.abs(a) - as)) : 0), dCT = 0.5 * s * cl(l) * x * x * dr;
    CT += dCT; CTm += 4 * F(l) * l * l * x * dr; CPi += l * dCT; CP0 += 0.5 * s * cd * x ** 3 * dr; aMax = Math.max(aMax, a);
    r.push(x); lam.push(l); dT.push(dCT / dr); al.push(a); Fs.push(F(l));
  }
  return { CT, CTm, CPi, CP0, r, lam, dT, al, F: Fs, aMax };
}
const hoverGeom = (i, p) => ({ Nb: p.Nb, sigma: p.sigma, tap: i.taper, tw: p.tw, ideal: i.twist_law === 'ideal (θtip/r)', cla: i.cla, clmax: i.clmax, cd0: i.cd0, cd2: i.cd2, rc: i.root_cut, tipLoss: i.tip_loss, nR: Math.round(i.nR) });

const hover = {
  id: 'hover', title: 'Hover by blade-element momentum theory', fidelity: 'numerical',
  summary: 'Radial inflow and loading of a hovering rotor from annulus-by-annulus momentum balance with Prandtl tip loss, giving thrust and power coefficients, figure of merit and the collective pitch needed to carry the weight.',
  equations: ['Blade element momentum theory', 'Rotor momentum theory', 'Rotor thrust and torque equations', 'Induced velocity equations', 'Prandtl tip-loss correction', 'Coaxial / tandem overlap interference (momentum theory)'],
  applicable: rotorOnly,
  inputs: [...MASS, ...ROTOR,
    { key: 'taper', label: 'Blade taper ratio (tip/root chord)', unit: '-', default: 1, min: 0.2, max: 1.5, group: 'Rotor' },
    { key: 'twist_law', label: 'Twist law', type: 'select', options: ['linear', 'ideal (θtip/r)'], default: 'linear', group: 'Rotor' },
    { key: 'root_cut', label: 'Root cut-out r/R', unit: '-', default: 0.15, min: 0, max: 0.5, group: 'Rotor' },
    { key: 'cd2', label: 'Profile drag rise d²cd/dα²', unit: '1/rad²', default: 0.3, min: 0, max: 2, group: 'Blade section', help: 'cd = cd0 + cd2·α² below stall' },
    { key: 'tip_loss', label: 'Prandtl tip loss', type: 'bool', default: true, group: 'Model' },
    { key: 'download', label: 'Fuselage download fraction', unit: '-', default: 0.03, min: 0, max: 0.3, group: 'Aircraft', help: 'Vertical drag of the airframe in the rotor wash; 0.02–0.05 for helicopters' },
    { key: 'overlap', label: 'Rotor disk overlap fraction m′', unit: '-', default: 0, min: 0, max: 1, group: 'Model', help: '0 isolated rotors, 1 coaxial; induced power factor κint = sqrt(2/(2 − m′))' },
    { key: 'ki_corr', label: 'Residual induced-power correction', unit: '-', default: 1.05, min: 0.9, max: 1.4, group: 'Calibration', help: 'Wake contraction, swirl and unsteadiness not captured by BEMT; calibrate on hover-stand power' },
    ...TAIL, ...ATM,
    { key: 'nR', label: 'Radial stations', unit: '', default: 60, min: 6, max: 2000, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => ({ ...baseDefaults(c, up, d), ...tailDefaults(c), download: c.meta.type === 'helicopter' ? 0.03 : 0.02 }),
  run(i) {
    const p = rotorOf(i), g = hoverGeom(i, p), W = i.mass_kg * G0, T = (W * (1 + i.download)) / i.n_rotors, CTreq = T / p.qA, warnings = [];
    let th = N.findRoot((t) => hoverBEMT(g, t).CT - CTreq, -0.05, 0.7, 30, 1e-12), trimmed = Number.isFinite(th);
    if (!trimmed) { th = N.goldenSection((t) => -hoverBEMT(g, t).CT, 0, 0.7, 1e-6); warnings.push('The rotor cannot produce the required thrust at any collective: the blades stall first. Results are shown at maximum thrust. Increase solidity, tip speed or radius.'); }
    const b = hoverBEMT(g, th), kint = Math.sqrt(2 / (2 - N.clamp(i.overlap, 0, 1))), CPi = b.CPi * kint * i.ki_corr, CP = CPi + b.CP0, Pu = p.qA * p.Vt;
    const Tr = b.CT * p.qA, Pmr = i.n_rotors * CP * Pu, Pideal = (i.n_rotors * Tr ** 1.5) / Math.sqrt(2 * p.rho * p.A), FM = Pideal / Pmr, vh = Math.sqrt(Tr / (2 * p.rho * p.A));
    const tr = tailRotor(i, p.rho, (CP * Pu) / p.Om, 0), kap = b.CPi / (b.CT ** 1.5 / Math.SQRT2), aS = i.clmax / i.cla;
    if (b.aMax > aS) warnings.push(`Blade sections reach ${fmt(N.deg(b.aMax))}° incidence, beyond the stall angle of ${fmt(N.deg(aS))}°: profile power is rising steeply.`);
    if (p.Mtip > 0.72) warnings.push(`Hover tip Mach number is ${p.Mtip.toFixed(2)}; compressibility drag rise is not included in the hover polar.`);
    if (b.CT / p.sigma > 0.14) warnings.push('Blade loading CT/σ exceeds 0.14: little or no stall margin is left for manoeuvre or gusts.');
    // hover polar: sweep collective
    const ths = N.linspace(0.02, Math.max(th * 1.25, 0.12), 14), pol = ths.map((t) => { const q = hoverBEMT(g, t), cp = q.CPi * kint * i.ki_corr + q.CP0; return [q.CT / p.sigma, q.CT > 0 ? q.CT ** 1.5 / Math.SQRT2 / cp : 0, cp / p.sigma]; });
    return {
      kpis: [
        { key: 'hover_power_W', label: 'Rotor shaft power in hover (all lifting rotors)', value: Pmr, unit: 'W' },
        { key: 'FM', label: 'Figure of merit', value: FM, unit: '-', status: FM > 0.65 ? 'ok' : FM > 0.5 ? 'warn' : 'bad', note: 'Good full-scale rotors reach 0.70–0.80; small rotors 0.55–0.70' },
        { key: 'CT', label: 'Thrust coefficient', value: b.CT, unit: '-' },
        { key: 'CP', label: 'Power coefficient', value: CP, unit: '-' },
        { key: 'CT_sigma', label: 'Blade loading CT/σ', value: b.CT / p.sigma, unit: '-', status: b.CT / p.sigma < 0.12 ? 'ok' : 'warn', note: 'Design values 0.06–0.10; stall-limited near 0.12–0.14' },
        { key: 'collective_deg', label: 'Collective pitch at 75% radius', value: N.deg(th), unit: 'deg', status: trimmed ? 'ok' : 'bad' },
        { key: 'v_induced_ms', label: 'Momentum-theory induced velocity', value: vh, unit: 'm/s' },
        { key: 'kappa_bemt', label: 'Induced power factor from BEMT', value: kap * i.ki_corr, unit: '-', note: 'Ratio of induced power to the uniform-inflow ideal' },
        { key: 'kappa_int', label: 'Rotor–rotor interference factor', value: kint, unit: '-' },
        { key: 'P_induced_W', label: 'Induced power', value: i.n_rotors * CPi * Pu, unit: 'W' },
        { key: 'P_profile_W', label: 'Profile power', value: i.n_rotors * b.CP0 * Pu, unit: 'W' },
        { key: 'P_tail_W', label: 'Tail-rotor power', value: tr.P, unit: 'W' },
        { key: 'hover_power_installed_W', label: 'Hover power including tail rotor', value: Pmr + tr.P, unit: 'W', note: 'Excludes transmission and accessory losses' },
        { key: 'thrust_per_rotor_N', label: 'Thrust per rotor', value: Tr, unit: 'N' },
        { key: 'disk_loading_Pa', label: 'Disk loading', value: Tr / p.A, unit: 'N/m²' },
        { key: 'power_loading_N_W', label: 'Power loading', value: (i.n_rotors * Tr) / Pmr, unit: 'N/W' },
        { key: 'tip_mach_hover', label: 'Tip Mach number', value: p.Mtip, unit: '-' },
        { key: 'alpha_max_deg', label: 'Peak section incidence', value: N.deg(b.aMax), unit: 'deg', status: b.aMax < aS ? 'ok' : 'warn' },
      ],
      plots: [
        { type: 'line', title: 'Radial induced velocity', xlabel: 'Radius r/R [-]', ylabel: 'Induced velocity [m/s]', series: [{ name: 'BEMT', x: b.r, y: b.lam.map((l) => l * p.Vt) }, { name: 'Uniform (momentum)', x: [0, 1], y: [vh, vh], style: 'dash' }] },
        { type: 'line', title: 'Radial thrust loading', xlabel: 'Radius r/R [-]', ylabel: 'Thrust per unit span, all blades [N/m]', series: [{ name: 'dT/dr', x: b.r, y: b.dT.map((v) => (v * p.qA) / p.R) }] },
        { type: 'line', title: 'Section incidence and tip-loss factor', xlabel: 'Radius r/R [-]', ylabel: 'Incidence [deg] / F×10 [-]', series: [{ name: 'Incidence α [deg]', x: b.r, y: b.al.map(N.deg) }, { name: 'Prandtl F × 10', x: b.r, y: b.F.map((v) => 10 * v), style: 'dash' }], annotations: [{ y: N.deg(aS), label: 'Stall' }] },
        { type: 'line', title: 'Hover efficiency versus blade loading', xlabel: 'Blade loading CT/σ [-]', ylabel: 'Figure of merit [-]', series: [{ name: 'FM', x: pol.map((v) => v[0]), y: pol.map((v) => v[1]) }], annotations: [{ x: b.CT / p.sigma, label: 'Operating point' }] },
      ],
      outputs: { solidity: p.sigma, v_tip_ms: p.Vt, bemt_momentum_residual: Math.abs(b.CTm - b.CT) / (Math.abs(b.CT) || 1) },
      warnings, models: ['Blade-element momentum theory (annular)', 'Prandtl tip-loss function', 'Parabolic section polar with hard stall limit (generic)', 'Overlap-area momentum interference factor', 'Momentum tail rotor'],
      assumptions: ['Hover out of ground effect, rigid blades, no swirl or wake contraction', 'Section data are a generic parametric polar, not tables for a specific aerofoil', 'Download applied as a simple thrust increment', 'Figure of merit is for the lifting-rotor system including interference', 'Download fraction, residual induced-power correction and section polar constants are typical values to be calibrated on hover-stand data'],
    };
  },
  convergence: { param: 'nR', label: 'Radial stations', levels: [10, 20, 40, 80, 160], metric: 'hover_power_W' },
  calibration: { params: [{ key: 'ki_corr', min: 0.9, max: 1.4 }, { key: 'cd0', min: 0.005, max: 0.04 }], sweep: 'mass_kg', target: 'hover_power_W', note: 'Hover-stand or flight-test shaft power versus thrust (gross mass) at known density.' },
  verify() {
    const i = { mass_kg: 1000, n_rotors: 1, R: 4, Nb: 3, chord: 0.3, twist_deg: 0, rpm: 400, cla: 5.73, cd0: 0, clmax: 1e9, taper: 1, twist_law: 'ideal (θtip/r)', root_cut: 0, cd2: 0, tip_loss: false, download: 0, overlap: 0, ki_corr: 1, tr_R: 0, tr_chord: 0, tr_Nb: 0, tr_rpm: 0, tr_arm: 0, alt_m: 0, dISA: 0, nR: 40 };
    const r = hover.run(i), o = N.kv(r), T = 1000 * G0, A = Math.PI * 16, rho = isa(0).rho;
    const c = N.kv(hover.run({ ...i, overlap: 1, n_rotors: 2 }));
    const t = hover.run({ ...i, twist_law: 'linear', twist_deg: -10, tip_loss: true, cd0: 0.01, nR: 80 });
    return [
      N.check('Ideal-twist, loss-free rotor: P = T^1.5/sqrt(2ρA)', o.hover_power_W, T ** 1.5 / Math.sqrt(2 * rho * A), 1e-8, 'Rankine–Froude momentum theory'),
      N.check('Ideal figure of merit = 1', o.FM, 1, 1e-8, 'Definition of figure of merit'),
      N.check('Uniform inflow v = sqrt(T/2ρA)', r.plots[0].series[0].y[20], Math.sqrt(T / (2 * rho * A)), 1e-8, 'Momentum theory'),
      N.check('Coaxial momentum interference κint = √2', (c.hover_power_W / 2) / ((T / 2) ** 1.5 / Math.sqrt(2 * rho * A)), Math.SQRT2, 1e-8, 'Momentum theory for two rotors sharing one disk area'),
      N.check('Annulus momentum thrust equals blade-element thrust', t.outputs.bemt_momentum_residual + 1, 1, 1e-9, 'Conservation of momentum in each annulus'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.FM < 0.6) out.push({ severity: 'advise', title: 'Hover efficiency is low', detail: `Figure of merit ${o.FM.toFixed(2)}: ${(100 * o.P_profile_W / o.hover_power_W).toFixed(0)}% of power is profile drag.`, action: o.P_profile_W / o.hover_power_W > 0.3 ? 'Reduce tip speed or solidity, or use a lower-drag section: profile power scales with σ·cd0·Vtip³.' : 'Increase blade twist or taper to even out the inflow, and reduce disk overlap.', basis: 'FM = ideal induced power / actual power' });
    if (o.CT_sigma > 0.12) out.push({ severity: 'warn', title: 'Blade loading is close to the stall limit', detail: `CT/σ = ${o.CT_sigma.toFixed(3)} with peak incidence ${o.alpha_max_deg.toFixed(1)}°.`, action: 'Add blade area (chord or blade count) or tip speed; check hot-and-high conditions, which raise CT at constant weight.', basis: 'Rotor stall boundary CT/σ ≈ 0.12–0.14' });
    out.push({ severity: 'info', title: 'Disk loading sets the hover energy bill', detail: `Induced power is ${(100 * o.P_induced_W / o.hover_power_W).toFixed(0)}% of hover power at ${o.disk_loading_Pa.toFixed(0)} N/m². Ideal power per unit thrust scales with sqrt(disk loading).`, action: 'A 10% larger radius cuts induced power by about 9% and the matching fuel or battery energy; trade it against blade mass and footprint in Suite 23.', basis: 'P/T = sqrt(DL / 2ρ)' });
    if (o.kappa_int > 1.01) out.push({ severity: 'info', title: 'Rotor overlap costs induced power', detail: `Interference factor κint = ${o.kappa_int.toFixed(3)}.`, action: 'Increase rotor spacing or vertical separation; measured coaxial penalties are lower than this same-plane momentum value because the lower rotor ingests a contracted wake.', basis: 'Momentum theory on the overlapped area' });
    return out;
  },
};

// ---- (b) forward flight ---------------------------------------------------------------------
const FF_INPUTS = [...MASS, ...ROTOR,
  { key: 'V', label: 'True airspeed', unit: 'm/s', default: 72, min: 0, max: 150, group: 'Flight' },
  { key: 'gamma_deg', label: 'Flight-path angle (climb +)', unit: 'deg', default: 0, min: -30, max: 30, group: 'Flight' },
  { key: 'f_plate', label: 'Equivalent flat-plate drag area', unit: 'm²', default: 3.35, min: 0, group: 'Aircraft', help: 'Parasite drag D = ½ρV²·f' },
  { key: 'wing_lift_frac', label: 'Share of weight carried by a wing', unit: '-', default: 0, min: 0, max: 0.95, group: 'Aircraft', help: '0 for pure rotorcraft; compound and transitioning eVTOL aircraft offload the rotors' },
  { key: 'lock', label: 'Lock number γ', unit: '-', default: 8.2, min: 0.5, max: 20, group: 'Blade dynamics', help: 'ρ·a·c·R⁴ / blade flap inertia; 5–10 for articulated rotors' },
  { key: 'inflow', label: 'Induced inflow model', type: 'select', options: ['Drees linear', 'uniform'], default: 'Drees linear', group: 'Model' },
  { key: 'tip_loss', label: 'Tip-loss effective radius', type: 'bool', default: true, group: 'Model' },
  { key: 'M_dd', label: 'Section drag-divergence Mach', unit: '-', default: 0.8, min: 0.6, max: 0.95, group: 'Blade section' },
  KAPPA, ...TAIL, ...ATM];
const ffDefaults = (c, up, d) => { const V = ffSpeed(c, d); return { ...baseDefaults(c, up, d), ...tailDefaults(c), V, f_plate: c.rotor.flat_plate_m2 || undefined, wing_lift_frac: wingShare(c, d, V), lock: c.rotor.lock || undefined }; };
const ffOpts = (i, p, V, nR, nPsi, passes = 2) => {
  const Wr = (i.mass_kg * G0 * (1 - (i.wing_lift_frac || 0))) / i.n_rotors, g = N.rad(i.gamma_deg || 0);
  return { V, Wr: Wr * Math.cos(g), Dr: (0.5 * p.rho * V * V * i.f_plate) / i.n_rotors + Wr * Math.sin(g), gam: 0, kappa: i.kappa, lock: i.lock, drees: i.inflow === 'Drees linear', tipLoss: i.tip_loss, nR, nPsi, passes };
};

const forward = {
  id: 'forward', title: 'Forward flight: trim, flapping, power and stall limits', fidelity: 'numerical',
  summary: 'Trims the rotor at speed with Glauert inflow and rigid-blade flapping, integrates blade loads around the azimuth for power, and maps blade incidence to locate retreating-blade stall and advancing-tip compressibility limits.',
  equations: ['Blade element theory', 'Rotor inflow equations (Glauert, Drees)', 'Blade flapping equations', 'Blade feathering equations', 'Rotor thrust and torque equations', 'Actuator-disk momentum equations'],
  applicable: rotorOnly,
  inputs: [...FF_INPUTS,
    { key: 'nR', label: 'Radial stations', unit: '', default: 20, min: 6, max: 200, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nPsi', label: 'Azimuth stations', unit: '', default: 36, min: 12, max: 360, step: 1, discrete: true, group: 'Numerics' }],
  defaults: ffDefaults,
  run(i, ctx) {
    const p = rotorOf(i), nR = Math.round(i.nR), nPsi = Math.round(i.nPsi), warnings = [], W = i.mass_kg * G0;
    const t = trimFF(p, ffOpts(i, p, i.V, nR, nPsi, 3)), tr = tailRotor(i, p.rho, t.Q, i.V), aS = i.clmax / i.cla;
    if (!t.converged) warnings.push('Trim did not converge: the rotor cannot carry this thrust at this speed (stall or reverse flow dominate). Values are indicative only.');
    // incidence map
    const grid = N.range(nR, () => new Array(nPsi).fill(0)); diskLoads(p, t.s, t.x, grid);
    const rs = N.range(nR, (j) => p.rc + ((j + 0.5) * (1 - p.rc)) / nR), ps = N.range(nPsi + 1, (k) => (360 * k) / nPsi);
    let stalled = 0, area = 0; for (let j = 0; j < nR; j++) for (let k = 0; k < nPsi; k++) { area += rs[j]; if (Math.abs(grid[j][k]) > aS && rs[j] + t.mu * Math.sin(N.rad(ps[k])) > 0) stalled += rs[j]; }
    // power curve
    const Vmax = Math.max(8, Math.min(Math.max(1.35 * i.V, 0.22 * p.Vt), 0.45 * p.Vt, 0.95 * p.a - p.Vt)), Vs = N.linspace(0, Vmax, 16), cur = { Pi: [], P0: [], Pp: [], Pt: [], Pe: [], Pb: [] };
    const energy = (V, q, tq) => { const mu = V / p.Vt; return { Pi: i.n_rotors * q.T * q.li * p.Vt, P0: ((i.n_rotors * p.qA * p.Vt * p.sigma * i.cd0) / 8) * (1 + 4.65 * mu * mu), Pp: 0.5 * p.rho * V ** 3 * i.f_plate + W * (1 - i.wing_lift_frac) * V * Math.sin(N.rad(i.gamma_deg)), Pt: tq.P }; };
    Vs.forEach((V, k) => {
      ctx?.progress?.(k / 24, 'Power curve');
      const q = trimFF(p, ffOpts(i, p, V, Math.min(nR, 12), Math.min(nPsi, 24), 2)), tq = tailRotor(i, p.rho, q.Q, V), e = energy(V, q, tq);
      cur.Pi.push(e.Pi); cur.P0.push(e.P0); cur.Pp.push(e.Pp); cur.Pt.push(e.Pt); cur.Pe.push(e.Pi + e.P0 + e.Pp + e.Pt); cur.Pb.push(q.converged ? i.n_rotors * q.P + tq.P : NaN);
    });
    const eOp = energy(i.V, t, tr), PeOp = eOp.Pi + eOp.P0 + eOp.Pp + eOp.Pt;
    const Pc = cur.Pb.map((v, k) => (Number.isFinite(v) ? v : cur.Pe[k])), kbe = N.argmin(Pc), kbr = N.argmin(Pc.map((v, k) => (k ? v / Vs[k] : Infinity)));
    // retreating-blade stall boundary: CT/σ at which α(1.0, 270°) reaches the static stall angle
    const mus = N.linspace(0.08, 0.44, 6), lim = mus.map((mu) => {
      let lo = 0.02, hi = 0.3;
      for (let k = 0; k < 7; k++) { const cts = 0.5 * (lo + hi), V = mu * p.Vt, q = trimFF(p, { ...ffOpts(i, p, V, 8, 16, 1), Wr: cts * p.sigma * p.qA, Dr: (0.5 * p.rho * V * V * i.f_plate) / i.n_rotors }); if (q.aRet > aS || !q.converged) hi = cts; else lo = cts; }
      return 0.5 * (lo + hi);
    });
    const cts = t.CT / p.sigma, ctLim = N.interp1(mus, lim, t.mu), Madv = (p.Vt + i.V) / p.a, Ptot = i.n_rotors * t.P + tr.P;
    if (t.aRet > aS) warnings.push(`Retreating-tip incidence ${fmt(N.deg(t.aRet))}° exceeds the static stall angle ${fmt(N.deg(aS))}°: expect stall flutter, high control loads and vibration.`);
    if (Madv > 0.9) warnings.push(`Advancing-tip Mach number ${Madv.toFixed(2)} is beyond 0.9: strong compressibility drag and noise; the generic drag-rise term is not validated there.`);
    if (t.mu > 0.4) warnings.push('Advance ratio above 0.4: the reverse-flow region is large and the rigid-blade first-harmonic trim loses accuracy.');
    if (i.wing_lift_frac > 0) warnings.push(`A wing carries ${(100 * i.wing_lift_frac).toFixed(0)}% of the weight at this speed; the rotors are analysed in edgewise, partially unloaded flight (transition condition).`);
    return {
      kpis: [
        { key: 'P_cruise_W', label: 'Rotor shaft power at this speed (all rotors + tail rotor)', value: Ptot, unit: 'W', status: t.converged ? 'ok' : 'warn' },
        { key: 'P_energy_W', label: 'Power by the energy method', value: PeOp, unit: 'W', note: 'Induced + profile + parasite + tail rotor; cross-check of the blade-element torque' },
        { key: 'mu_adv', label: 'Advance ratio μ', value: t.mu, unit: '-' },
        { key: 'lambda_i', label: 'Induced inflow ratio', value: t.li, unit: '-' },
        { key: 'v_induced_ff_ms', label: 'Mean induced velocity', value: t.li * p.Vt, unit: 'm/s' },
        { key: 'alpha_tpp_deg', label: 'Tip-path-plane forward tilt', value: N.deg(t.aT), unit: 'deg' },
        { key: 'theta75_ff_deg', label: 'Collective at 75% radius', value: N.deg(t.th75), unit: 'deg' },
        { key: 'cyc_lon_deg', label: 'Longitudinal cyclic B1 (to the tip-path plane)', value: N.deg(-t.th1s), unit: 'deg' },
        { key: 'cyc_lat_deg', label: 'Lateral cyclic A1 (to the tip-path plane)', value: N.deg(-t.th1c), unit: 'deg' },
        { key: 'coning_deg', label: 'Coning angle a0', value: N.deg(t.a0), unit: 'deg', status: t.a0 < N.rad(8) ? 'ok' : 'warn' },
        { key: 'a1_deg', label: 'Longitudinal flapping a1 (rearward, to the no-feathering plane)', value: N.deg(t.a1), unit: 'deg' },
        { key: 'b1_deg', label: 'Lateral flapping b1 (to the no-feathering plane)', value: N.deg(t.b1), unit: 'deg' },
        { key: 'CT_sigma_ff', label: 'Blade loading CT/σ', value: cts, unit: '-', status: cts < ctLim ? 'ok' : 'bad' },
        { key: 'CT_sigma_limit', label: 'Retreating-stall limit on CT/σ at this μ', value: ctLim, unit: '-', note: 'Gessow–Myers criterion with the static stall angle' },
        { key: 'alpha_ret_deg', label: 'Retreating-tip incidence α(1.0, 270°)', value: N.deg(t.aRet), unit: 'deg', status: t.aRet < aS ? 'ok' : 'bad' },
        { key: 'stalled_disk_pct', label: 'Disk area beyond static stall', value: (100 * stalled) / area, unit: '%' },
        { key: 'M_adv_tip', label: 'Advancing-tip Mach number', value: Madv, unit: '-', status: Madv < 0.88 ? 'ok' : Madv < 0.94 ? 'warn' : 'bad', note: 'Cruise values of 0.80–0.88 are usual; drag divergence and impulsive noise grow above about 0.9' },
        { key: 'V_be_ms', label: 'Best-endurance speed (minimum power)', value: Vs[kbe], unit: 'm/s' },
        { key: 'V_br_ms', label: 'Best-range speed (minimum power/speed)', value: Vs[kbr], unit: 'm/s' },
        { key: 'P_min_W', label: 'Minimum power required', value: Pc[kbe], unit: 'W' },
        { key: 'P_tail_ff_W', label: 'Tail-rotor power', value: tr.P, unit: 'W' },
        { key: 'H_force_N', label: 'Rotor in-plane H-force (per rotor)', value: t.H, unit: 'N' },
      ],
      plots: [
        { type: 'line', title: 'Power required versus speed', xlabel: 'True airspeed [m/s]', ylabel: 'Power [kW]', series: [{ name: 'Induced', x: Vs, y: cur.Pi.map((v) => v / 1e3) }, { name: 'Profile', x: Vs, y: cur.P0.map((v) => v / 1e3) }, { name: 'Parasite + climb', x: Vs, y: cur.Pp.map((v) => v / 1e3) }, { name: 'Tail rotor', x: Vs, y: cur.Pt.map((v) => v / 1e3) }, { name: 'Total (energy method)', x: Vs, y: cur.Pe.map((v) => v / 1e3), style: 'dash' }, { name: 'Total (blade-element trim)', x: Vs, y: cur.Pb.map((v) => v / 1e3), style: 'line+points' }], annotations: [{ x: Vs[kbe], label: 'Best endurance' }, { x: Vs[kbr], label: 'Best range' }] },
        { type: 'heat', title: 'Blade section incidence over the disk', xlabel: 'Azimuth ψ [deg] (90° = advancing side)', ylabel: 'Radius r/R [-]', zlabel: 'Incidence [deg]', x: ps, y: rs, z: grid.map((row) => [...row, row[0]].map((a) => N.clamp(N.deg(a), -10, 30))), contours: 12 },
        { type: 'line', title: 'Retreating-blade stall boundary', xlabel: 'Advance ratio μ [-]', ylabel: 'Blade loading CT/σ [-]', series: [{ name: 'Stall limit (tip incidence = static stall)', x: mus, y: lim }, { name: 'Operating point', x: [t.mu], y: [cts], style: 'points' }] },
        { type: 'polar', title: 'Incidence around the azimuth at 85% radius', rlabel: 'Incidence [deg]', series: [{ name: 'α at 0.85R', theta_deg: ps, r: [...grid[Math.min(nR - 1, Math.round(((0.85 - p.rc) / (1 - p.rc)) * nR - 0.5))], grid[Math.min(nR - 1, Math.round(((0.85 - p.rc) / (1 - p.rc)) * nR - 0.5))][0]].map((a) => Math.max(0, N.deg(a))) }] },
      ],
      tables: [{ title: 'Trim solution', columns: ['Quantity', 'Value', 'Unit'], rows: [['Thrust per rotor', t.T, 'N'], ['CT', t.CT, '-'], ['Total inflow ratio λ (TPP)', t.lam, '-'], ['Drees kx', t.s.kx, '-'], ['Shaft torque per rotor', t.Q, 'N m'], ['Tail-rotor thrust', tr.T, 'N']] }],
      outputs: { trim_converged: t.converged ? 1 : 0 },
      warnings, models: ['Blade-element integration over radius and azimuth', 'Glauert momentum inflow with Drees linear gradient', 'Rigid, centrally hinged blade: first-harmonic flap equilibrium', 'Generic section polar with stall limit, Prandtl–Glauert slope and Lock-type drag rise', 'Energy-method power build-up (profile factor 1 + 4.65μ²)'],
      assumptions: ['Steady flight, isolated rotor, no fuselage or rotor–rotor interference in forward flight', 'Hinge offset neglected in trim (see the blade-dynamics analysis)', 'Flapping and cyclic are referred to the tip-path and no-feathering planes; fuselage attitude is not solved', 'Static stall angle used for the stall boundary; dynamic stall delays it (see the dynamic-stall analysis)', 'Induced-power factor κ, section clmax and drag-divergence Mach number are typical values, not data for a specific blade'],
    };
  },
  convergence: { param: 'nR', label: 'Radial stations', levels: [8, 16, 32, 64], metric: 'P_cruise_W' },
  calibration: { params: [{ key: 'kappa', min: 1, max: 1.6 }, { key: 'f_plate', min: 0, max: 20 }, { key: 'cd0', min: 0.005, max: 0.04 }], sweep: 'V', target: 'P_cruise_W', note: 'Flight-test or wind-tunnel shaft power versus airspeed at known mass and density.' },
  verify() {
    // Linear aerodynamics, uniform inflow, untwisted blade, full radius: classical closed-form flapping results.
    const i = { R: 6, Nb: 4, chord: 0.4, twist_deg: 0, rpm: 300, cla: 5.7, cd0: 0.008, alt_m: 0, dISA: 0, root_cut: 0 }, p = rotorOf(i, { linear: true });
    const t = trimFF(p, { V: 0.2 * p.Vt, Wr: 0.007 * p.qA, Dr: 0.0004 * p.qA, kappa: 1, lock: 8, drees: false, tipLoss: false, nR: 200, nPsi: 48, passes: 1 }), mu = t.mu, lam = t.lam, th = t.th75;
    const a1 = ((8 / 3) * th * mu - 2 * mu * lam) / (1 + 1.5 * mu * mu), a0 = (8 / 2) * (th * (0.25 + mu * mu / 4) - (a1 * mu) / 3 - lam / 3), b1 = ((4 / 3) * mu * a0) / (1 + mu * mu / 2);
    const CT = (p.sigma * p.cla / 2) * ((th / 3) * (1 + 1.5 * mu * mu) - (mu * a1) / 2 - lam / 2);
    return [
      N.check('Longitudinal flapping a1 = (8/3·θ·μ − 2μλ)/(1 + 3μ²/2)', t.a1, a1, 2e-4, 'Classical rigid-blade flapping in the tip-path plane (Johnson, Helicopter Theory)'),
      N.check('Lateral flapping b1 = (4/3)·μ·a0/(1 + μ²/2)', t.b1, b1, 2e-4, 'Classical rigid-blade flapping, uniform inflow'),
      N.check('Coning a0 = γ/2·[θ(1+μ²)/4 − μ·a1/3 − λ/3]', t.a0, a0, 2e-4, 'Mean flap-moment balance'),
      N.check('Thrust CT = σa/2·[θ/3·(1+3μ²/2) − μ·a1/2 − λ/2]', t.CT, CT, 2e-4, 'Blade-element thrust integral'),
      N.check('Glauert inflow: λi·sqrt(μ² + λ²) = CT/2', t.li0 * Math.hypot(t.mu, t.lam), t.CT / 2, 1e-9, 'Glauert high-speed momentum theory'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.alpha_ret_deg > N.deg(i.clmax / i.cla) || o.CT_sigma_ff > o.CT_sigma_limit) out.push({ severity: 'critical', title: 'Retreating-blade stall at this flight condition', detail: `CT/σ = ${o.CT_sigma_ff.toFixed(3)} against a limit of ${o.CT_sigma_limit.toFixed(3)} at μ = ${o.mu_adv.toFixed(2)}.`, action: 'Reduce speed, mass or altitude, or add blade area / tip speed. Offloading the rotor with a wing raises the speed limit.', basis: 'Gessow–Myers retreating-tip incidence criterion' });
    else if (o.CT_sigma_ff > 0.85 * o.CT_sigma_limit) out.push({ severity: 'warn', title: 'Thin retreating-stall margin', detail: `Blade loading is ${(100 * o.CT_sigma_ff / o.CT_sigma_limit).toFixed(0)}% of the stall limit.`, action: 'Check manoeuvre load factors and hot-and-high conditions; the limit falls with altitude at constant weight.', basis: 'Retreating-blade stall boundary' });
    if (o.M_adv_tip > 0.85) out.push({ severity: o.M_adv_tip > 0.9 ? 'warn' : 'advise', title: 'Advancing-tip compressibility', detail: `Advancing-tip Mach number ${o.M_adv_tip.toFixed(2)}${o.M_adv_tip > 0.9 ? '' : ': normal for a helicopter at cruise speed, but power and noise rise quickly from here'}.`, action: 'Reduce rotor speed in cruise, use thinner or swept tips; this also lowers high-speed impulsive noise (Suite 11).', basis: 'Section drag divergence: thin tip sections tolerate about Mach 0.85–0.9 before the power rise becomes steep' });
    if (i.V > 1.08 * o.V_br_ms || i.V < 0.85 * o.V_br_ms) out.push({ severity: 'advise', title: 'Cruise speed is away from best range', detail: `Best-range speed is about ${o.V_br_ms.toFixed(0)} m/s and best endurance ${o.V_be_ms.toFixed(0)} m/s; the analysis speed is ${i.V.toFixed(0)} m/s.`, action: 'Fly near the best-range speed to minimise fuel or battery energy and CO₂ per kilometre; loiter and search at the best-endurance speed.', basis: 'Minimum of P/V and of P on the power curve' });
    out.push({ severity: 'info', title: 'Parasite drag drives high-speed power', detail: `Flat-plate area ${i.f_plate.toFixed(2)} m²; parasite power grows with V³.`, action: 'Hub and landing-gear fairings typically remove 10–20% of flat-plate area and the same share of high-speed parasite power.', basis: 'P = ½ρV³f' });
    return out;
  },
};

// ---- (c) blade flap / lag dynamics ----------------------------------------------------------
/** RK4 time response of a rigid blade in flap and lag. Returns per-step arrays for the last revolution and convergence history. */
function flapSim(p, s, c) {
  const e = c.e, nu2 = 1 + (1.5 * e) / (1 - e) + c.kbeta, nz2 = Math.max(1e-4, (1.5 * c.e_lag) / (1 - c.e_lag) + c.kzeta), nz = Math.sqrt(nz2), g2 = c.lock / (2 * p.cla), r0 = Math.max(p.rc, e), nR = c.nR, dr = (1 - r0) / nR;
  const ev = (psi, y) => {
    const sp = Math.sin(psi), cp = Math.cos(psi), th = c.th75 + c.th1c * cp + c.th1s * sp; let Mf = 0, Ml = 0, Sz = 0;
    for (let j = 0; j < nR; j++) {
      const x = r0 + (j + 0.5) * dr, uT = x + s.mu * sp, uP = s.lam + s.lami * x * (s.kx * cp + s.ky * sp) + (x - e) * y[1] + s.mu * y[0] * cp;
      const q = section(p, th + p.tw * (x - 0.75), uT, uP); let fz = q[0], fx = q[1];
      if (x > s.B) { fz = 0; fx = p.cd0 * uT * Math.abs(uT); }
      Mf += (x - e) * fz * dr; Ml += (x - e) * fx * dr; Sz += fz * dr;
    }
    return { d: [y[1], g2 * Mf - nu2 * y[0], y[3], g2 * Ml - 2 * c.zeta_lag * nz * y[3] - nz2 * y[2] - 2 * y[0] * y[1]], Sz }; // Coriolis: flapping up drives the blade forward (lead)
  };
  const nS = c.nStep, sol = N.rk4((t, y) => ev(t, y).d, 0, [0, 0, 0, 0], TWO_PI * c.nRev, nS * c.nRev), last = sol.y.slice(-nS - 1, -1), psi = sol.t.slice(-nS - 1, -1);
  const prev = sol.y.slice(-2 * nS - 1, -nS - 1), perr = c.nRev > 1 ? Math.max(...last.map((y, k) => Math.abs(y[0] - prev[k][0]))) : NaN;
  const hist = N.range(c.nRev, (k) => N.mean(sol.y.slice(k * nS, (k + 1) * nS).map((y) => y[0])));
  const sz = [], bdd = []; last.forEach((y, k) => { const q = ev(psi[k], y); sz.push(q.Sz); bdd.push(q.d[1]); });
  const harm = (arr, n) => { let a = 0, b = 0; arr.forEach((v, k) => { a += v * Math.cos(n * psi[k]); b += v * Math.sin(n * psi[k]); }); return n ? [(2 * a) / nS, (2 * b) / nS] : [a / nS, 0]; };
  const beta = last.map((y) => y[0]), zeta = last.map((y) => y[2]);
  return { psi: psi.map((v) => v - psi[0]), beta, zeta, sz, bdd, harm, a0: harm(beta, 0)[0], a1: -harm(beta, 1)[0], b1: -harm(beta, 1)[1], zeta0: harm(zeta, 0)[0], zeta1: Math.hypot(...harm(zeta, 1)), nu: Math.sqrt(nu2), nz, perr, hist };
}

const bladedyn = {
  id: 'bladedyn', title: 'Blade flap and lead–lag dynamics, hub load harmonics', fidelity: 'numerical',
  summary: 'Integrates the rigid-blade flapping and lead–lag equations over many revolutions until the response repeats, then extracts flapping harmonics, natural frequencies and the vibratory hub force that the blades pass to the airframe.',
  equations: ['Blade flapping equations', 'Blade lead–lag equations', 'Centrifugal stiffening equations', 'Coriolis acceleration relations', 'Euler–Lagrange rotor dynamics equations (rigid blade)', 'Blade element theory'],
  applicable: rotorOnly,
  inputs: [...FF_INPUTS.filter((f) => f.group !== 'Tail rotor'), // the tail rotor plays no part in the blade response
    { key: 'control', label: 'Control setting', type: 'select', options: ['collective only (rotor free to flap)', 'trimmed cyclic (tip-path plane held)'], default: 'collective only (rotor free to flap)', group: 'Flight' },
    { key: 'hinge', label: 'Flap hinge offset e/R', unit: '-', default: 0.047, min: 0, max: 0.3, group: 'Blade dynamics', help: '0 teetering/central hinge, 0.03–0.06 articulated, 0.10–0.15 equivalent for hingeless' },
    { key: 'k_beta', label: 'Flap spring Kβ/(Ib·Ω²)', unit: '-', default: 0, min: 0, max: 2, group: 'Blade dynamics' },
    { key: 'hinge_lag', label: 'Lag hinge offset e/R', unit: '-', default: 0.047, min: 0.005, max: 0.3, group: 'Blade dynamics' },
    { key: 'k_zeta', label: 'Lag spring Kζ/(Ib·Ω²)', unit: '-', default: 0, min: 0, max: 4, group: 'Blade dynamics' },
    { key: 'zeta_lag', label: 'Lag damper ratio', unit: '-', default: 0.3, min: 0, max: 2, group: 'Blade dynamics', help: 'Fraction of critical damping of the lag mode supplied by the lag damper' },
    { key: 'blade_mass', label: 'Blade mass', unit: 'kg', default: 116, min: 0.001, group: 'Blade dynamics' },
    { key: 'nRev', label: 'Rotor revolutions', unit: '', default: 16, min: 3, max: 200, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nStep', label: 'Time steps per revolution', unit: '', default: 72, min: 16, max: 1440, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nR', label: 'Radial stations', unit: '', default: 16, min: 6, max: 200, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => { const { tr_R, tr_chord, tr_Nb, tr_rpm, tr_arm, ...f } = ffDefaults(c, up, d); return { ...f, hinge: c.rotor.hinge_offset, hinge_lag: Math.max(0.01, c.rotor.hinge_offset), blade_mass: c.rotor.blade_mass_kg || undefined }; },
  run(i) {
    const p = rotorOf(i), nR = Math.round(i.nR), nS = Math.round(i.nStep), warnings = [], free = i.control.startsWith('collective');
    const t = trimFF(p, ffOpts(i, p, i.V, nR, 36, 2));
    const s = free ? { ...t.s, lam: t.lam + t.mu * t.a1 } : t.s; // inflow through the control (no-feathering) plane when no cyclic is applied
    const r = flapSim(p, s, { th75: t.th75, th1c: free ? 0 : t.th1c, th1s: free ? 0 : t.th1s, e: i.hinge, kbeta: i.k_beta, e_lag: i.hinge_lag, kzeta: i.k_zeta, zeta_lag: i.zeta_lag, lock: i.lock, nR, nStep: nS, nRev: Math.round(i.nRev) });
    // dimensional blade root vertical shear: airload minus flap inertia
    const qb = 0.5 * p.rho * i.chord * p.R * p.Vt ** 2, Sb = (i.blade_mass * p.R * (1 - i.hinge)) / 2, shear = r.sz.map((v, k) => qb * v - Sb * p.Om ** 2 * r.bdd[k]);
    const nH = Math.min(2 * p.Nb, Math.floor(nS / 2) - 1), amp = N.range(nH, (n) => Math.hypot(...r.harm(shear, n + 1))), hub = p.Nb * Math.hypot(...r.harm(shear, p.Nb)), Tmean = p.Nb * r.harm(shear, 0)[0];
    if (nH < p.Nb) warnings.push(`Only ${nS} steps per revolution: the ${p.Nb}/rev blade-passage harmonic is close to the sampling limit; use more time steps.`);
    const Ib = (i.blade_mass * p.R ** 2 * (1 - i.hinge) ** 3) / 3, hubM = (p.Nb / 2) * (r.nu ** 2 - 1) * Ib * p.Om ** 2 * Math.hypot(r.a1, r.b1);
    const damp = (i.lock / 8) * (1 - i.hinge) ** 3 * (1 + i.hinge / 3), phase = N.deg(Math.atan2(damp, r.nu ** 2 - 1));
    if (!(r.perr < 2e-4)) warnings.push(`The response has not become periodic (revolution-to-revolution flap change ${fmt(N.deg(r.perr), 3)}°): add revolutions or lag damping.`);
    if (Math.abs(r.nz - 1) < 0.1) warnings.push('Lag frequency is within 10% of 1/rev: resonance with the rotor speed.');
    if (r.nz < 1 && i.zeta_lag < 0.05) warnings.push('Soft in-plane rotor with little lag damping: check ground resonance in Suite 10.');
    if (!t.converged) warnings.push('The underlying trim did not converge; controls are indicative.');
    return {
      kpis: [
        { key: 'nu_beta', label: 'Flap frequency ratio νβ', value: r.nu, unit: '/rev', note: '1.00 central hinge, 1.02–1.04 articulated, 1.08–1.15 hingeless' },
        { key: 'nu_zeta', label: 'Lead–lag frequency ratio νζ', value: r.nz, unit: '/rev', status: Math.abs(r.nz - 1) > 0.1 ? 'ok' : 'bad', note: '0.2–0.3 articulated, 0.6–0.8 soft in-plane, >1 stiff in-plane' },
        { key: 'coning_dyn_deg', label: 'Coning from the time response', value: N.deg(r.a0), unit: 'deg' },
        { key: 'a1_dyn_deg', label: 'Longitudinal flapping (rearward)', value: N.deg(r.a1), unit: 'deg' },
        { key: 'b1_dyn_deg', label: 'Lateral flapping', value: N.deg(r.b1), unit: 'deg' },
        { key: 'flap_phase_deg', label: 'Flap response phase lag at 1/rev', value: phase, unit: 'deg', note: '90° for a central hinge; less with hinge offset or flap spring' },
        { key: 'lag_mean_deg', label: 'Mean lag angle', value: N.deg(r.zeta0), unit: 'deg' },
        { key: 'lag_1p_deg', label: '1/rev lag amplitude', value: N.deg(r.zeta1), unit: 'deg' },
        { key: 'hub_vib_N', label: `${p.Nb}/rev vertical hub force amplitude`, value: hub, unit: 'N', note: 'Rigid-blade estimate: elastic blade modes, which usually dominate, are not included' },
        { key: 'hub_vib_g', label: 'Equivalent airframe vibration', value: hub / (i.mass_kg * G0), unit: 'g', status: hub / (i.mass_kg * G0) < 0.05 ? 'ok' : 'warn' },
        { key: 'hub_moment_Nm', label: 'Steady hub moment from tip-path tilt', value: hubM, unit: 'N m' },
        { key: 'thrust_dyn_N', label: 'Mean thrust from root shear', value: Tmean, unit: 'N' },
        { key: 'periodicity_deg', label: 'Revolution-to-revolution flap change', value: N.deg(r.perr), unit: 'deg', status: r.perr < 2e-4 ? 'ok' : 'warn' },
      ],
      plots: [
        { type: 'line', title: 'Flap and lag angle around the azimuth (converged revolution)', xlabel: 'Azimuth ψ [deg]', ylabel: 'Angle [deg]', series: [{ name: 'Flap β', x: r.psi.map(N.deg), y: r.beta.map(N.deg) }, { name: 'Lag ζ', x: r.psi.map(N.deg), y: r.zeta.map(N.deg) }] },
        { type: 'line', title: 'Blade root vertical shear', xlabel: 'Azimuth ψ [deg]', ylabel: 'Shear force [N]', series: [{ name: 'One blade', x: r.psi.map(N.deg), y: shear }] },
        { type: 'bar', title: 'Harmonics of blade root vertical shear', ylabel: 'Amplitude [N]', categories: N.range(nH, (n) => `${n + 1}/rev`), series: [{ name: 'One blade', y: amp }] },
        { type: 'line', title: 'Approach to the periodic solution', xlabel: 'Revolution [-]', ylabel: 'Mean flap angle over the revolution [deg]', series: [{ name: 'Coning', x: N.range(r.hist.length, (k) => k + 1), y: r.hist.map(N.deg), style: 'line+points' }] },
      ],
      warnings, models: ['Rigid blade with offset flap and lag hinges and root springs', 'Blade-element airloads with the generic section polar', 'Classical RK4 time marching to a periodic state', 'Hub filtering: only harmonics that are multiples of the blade count reach the non-rotating frame'],
      assumptions: ['Uniform blade mass; Lock number as entered', 'Prescribed inflow from the trim solution (no dynamic wake feedback)', 'Flap–lag coupling through the Coriolis term only; no torsion or pitch–flap coupling', 'In-plane hub forces and elastic bending are not modelled, so the vibratory hub force is a lower-bound estimate'],
    };
  },
  convergence: { param: 'nStep', label: 'Time steps per revolution', levels: [24, 48, 96, 192], metric: 'a1_dyn_deg' },
  calibration: { params: [{ key: 'lock', min: 1, max: 20 }, { key: 'hinge', min: 0, max: 0.3 }, { key: 'zeta_lag', min: 0, max: 2 }], sweep: 'V', target: 'a1_dyn_deg', note: 'Measured blade flapping harmonics versus airspeed from an instrumented rotor.' },
  verify() {
    const p = rotorOf({ R: 6, Nb: 4, chord: 0.4, twist_deg: 0, rpm: 300, cla: 5.7, cd0: 0.008, alt_m: 0, dISA: 0, root_cut: 0 }, { linear: true });
    const c = { th75: 0.14, th1c: 0, th1s: 0, e: 0, kbeta: 0, e_lag: 0.05, kzeta: 0, zeta_lag: 0.5, lock: 8, nR: 100, nStep: 96, nRev: 14 }, s0 = { mu: 0, lam: 0.05, lami: 0, kx: 0, ky: 0 };
    const h = flapSim(p, s0, c), hc = flapSim(p, s0, { ...c, th1c: 0.03 });
    const t = trimFF(p, { V: 0.15 * p.Vt, Wr: 0.006 * p.qA, Dr: 0.0002 * p.qA, kappa: 1, lock: 8, drees: false, tipLoss: false, nR: 100, nPsi: 48, passes: 1 });
    const f = flapSim(p, { ...t.s, lam: t.lam + t.mu * t.a1 }, { ...c, th75: t.th75 });
    return [
      N.check('Hover coning a0 = γ/8·(θ − 4λ/3)', h.a0, (8 / 8) * (0.14 - (4 * 0.05) / 3), 2e-4, 'Steady solution of the flapping equation'),
      N.check('Flap–feather equivalence: lateral cyclic gives equal flapping 90° later', -hc.b1, 0.03, 2e-4, 'Exact resonance response of a centrally hinged blade in hover'),
      N.check('Time-marched longitudinal flapping equals the harmonic-balance trim', f.a1, t.a1, 0.01, 'Flap–feather equivalence between the no-feathering and tip-path planes'),
      N.check('Flap frequency with hinge offset νβ² = 1 + 3e/(2(1 − e))', flapSim(p, s0, { ...c, e: 0.05, nRev: 3 }).nu ** 2, 1 + 0.075 / 0.95, 1e-12, 'Rigid blade with offset hinge'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.hub_vib_g > 0.05) out.push({ severity: 'warn', title: 'Vibratory hub force is significant', detail: `${o.hub_vib_N.toFixed(0)} N at blade-passage frequency (${o.hub_vib_g.toFixed(3)} g on the gross mass), before any elastic amplification.`, action: 'Consider more blades, hub absorbers or isolation, and keep airframe modes away from the blade-passage frequency (Suite 10). Lower vibration extends component fatigue life (Suite 9).', basis: 'Blade-count filtering of rotating-frame loads' });
    if (Math.abs(o.nu_zeta - 1) < 0.1) out.push({ severity: 'critical', title: 'Lead–lag frequency is near 1/rev', detail: `νζ = ${o.nu_zeta.toFixed(2)}.`, action: 'Change the lag hinge offset or lag stiffness to move the frequency away from the rotor speed.', basis: 'Resonance avoidance' });
    out.push({ severity: 'info', title: 'Hinge offset sets control power and phase', detail: `νβ = ${o.nu_beta.toFixed(3)}, phase lag ${o.flap_phase_deg.toFixed(0)}°, hub moment ${(o.hub_moment_Nm / 1e3).toFixed(1)} kN·m for the present tip-path tilt.`, action: 'Rig the swashplate phase to the computed lag; higher hinge offset gives more hub moment (agility, CG range) at the cost of vibration and gust response.', basis: 'Second-order flap response' });
    return out;
  },
};

// ---- (d) Pitt–Peters dynamic inflow ---------------------------------------------------------
const inflow = {
  id: 'inflow', title: 'Dynamic inflow response (Pitt–Peters three-state)', fidelity: 'reduced-order',
  summary: 'Shows how the rotor wake lags a collective-pitch step: thrust overshoots, then settles as the induced flow builds up through the Pitt–Peters uniform, lateral and longitudinal inflow states.',
  equations: ['Pitt–Peters dynamic inflow model', 'Rotor inflow equations', 'Coupled blade element–dynamic inflow formulation', 'Actuator-disk momentum equations'],
  applicable: rotorOnly,
  inputs: [...MASS, ...ROTOR.filter((f) => !['cd0', 'clmax'].includes(f.key)),
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 0, min: 0, max: 120, group: 'Flight' },
    { key: 'alpha_d_deg', label: 'Disk forward tilt', unit: 'deg', default: 0, min: -10, max: 20, group: 'Flight' },
    { key: 'd_theta_deg', label: 'Collective step', unit: 'deg', default: 1, min: -5, max: 5, group: 'Input' },
    { key: 'dCM', label: 'Aerodynamic pitching-moment step ΔCM', unit: '-', default: 0, min: -0.002, max: 0.002, group: 'Input', help: 'Excites the longitudinal inflow gradient state' },
    { key: 't_end_rev', label: 'Simulated time', unit: 'rev', default: 12, min: 2, max: 100, group: 'Numerics' },
    { key: 'nSteps', label: 'Time steps', unit: '', default: 600, min: 50, max: 20000, step: 1, discrete: true, group: 'Numerics' }, ...ATM],
  defaults: (c, up, d) => { const b = baseDefaults(c, up, d); delete b.cd0; return b; },
  run(i) {
    const p = rotorOf({ ...i, cd0: 0 }), CT0 = (i.mass_kg * G0) / i.n_rotors / p.qA, aD = N.rad(i.alpha_d_deg), mu = (i.V * Math.cos(aD)) / p.Vt, lc = (i.V * Math.sin(aD)) / p.Vt, sa = p.sigma * p.cla;
    const M = [8 / (3 * Math.PI), 16 / (45 * Math.PI), 16 / (45 * Math.PI)]; // apparent-mass terms (uncorrected thrust term)
    const ctOf = (th, l0) => (sa / 2) * ((th / 3) * (1 + 1.5 * mu * mu) - (p.tw * mu * mu) / 8 - (lc + l0) / 2); // quasi-steady blade-element thrust
    const l00 = lamI(CT0, mu, lc), th0 = (3 / (1 + 1.5 * mu * mu)) * ((2 * CT0) / sa + (p.tw * mu * mu) / 8 + (lc + l00) / 2), th1 = th0 + N.rad(i.d_theta_deg);
    const rhs = (th, cm) => (t, y) => {
      const lam = lc + y[0], VT = Math.max(1e-6, Math.hypot(mu, lam)), Vm = (mu * mu + lam * (lam + y[0])) / VT, X = Math.tan(Math.atan2(mu, Math.abs(lam) + 1e-12) / 2), k = (15 * Math.PI * X) / 64;
      // M·λ' + V·L̃⁻¹·λ = C with static gains L̃ = [[1/2, 0, −k], [0, 2(1+X²), 0], [k, 0, 2(1−X²)]]; the (0, c) block is inverted analytically
      const det = (1 - X * X) + k * k, Li = [[(2 * (1 - X * X)) / det, 0, k / det], [0, 1 / (2 * (1 + X * X)), 0], [-k / det, 0, 0.5 / det]];
      const f = [ctOf(th, y[0]), 0, cm], V = [VT, Vm, Vm];
      return [0, 1, 2].map((r) => (f[r] - V[r] * (Li[r][0] * y[0] + Li[r][1] * y[1] + Li[r][2] * y[2])) / M[r]);
    };
    // initial steady state at the old collective (thrust only): march to equilibrium, then apply the step
    const n = Math.round(i.nSteps), tEnd = TWO_PI * i.t_end_rev, pre = N.rk4(rhs(th0, 0), 0, [l00, 0, 0], 60, 600), y0 = pre.y[pre.y.length - 1];
    const sol = N.rk4(rhs(th1, i.dCM), 0, y0, tEnd, n), tt = sol.t.map((t) => t / p.Om), ct = sol.y.map((y) => ctOf(th1, y[0])), ctEnd = ct[n], ctQS = ctOf(th1, y0[0]);
    const peak = i.d_theta_deg >= 0 ? N.amax(ct) : N.amin(ct), over = (peak - ctEnd) / ((ctEnd - ctOf(th0, y0[0])) || 1e-12);
    const l63 = y0[0] + 0.632 * (sol.y[n][0] - y0[0]); let tau = NaN;
    for (let k = 1; k <= n; k++) if ((sol.y[k][0] - l63) * (sol.y[k - 1][0] - l63) <= 0) { tau = tt[k - 1] + ((l63 - sol.y[k - 1][0]) / ((sol.y[k][0] - sol.y[k - 1][0]) || 1e-30)) * (tt[k] - tt[k - 1]); break; }
    const tauHover = M[0] / (4 * Math.sqrt(CT0 / 2) + sa / 4) / p.Om, ds = Math.max(1, Math.floor(n / 300)), pick = (a) => a.filter((_, k) => k % ds === 0);
    return {
      kpis: [
        { key: 'tau_inflow_s', label: 'Inflow build-up time constant (63%)', value: tau, unit: 's' },
        { key: 'tau_inflow_rev', label: 'Time constant in rotor revolutions', value: (tau * p.Om) / TWO_PI, unit: 'rev' },
        { key: 'tau_hover_lin_s', label: 'Linearised hover time constant', value: tauHover, unit: 's', note: 'm11 / (4λh + σa/4) / Ω' },
        { key: 'thrust_overshoot', label: 'Thrust overshoot ratio', value: over, unit: '-', note: 'Peak thrust change divided by the final thrust change, minus one' },
        { key: 'lambda0_final', label: 'Final uniform inflow λ0', value: sol.y[n][0], unit: '-' },
        { key: 'lambda_c_final', label: 'Final longitudinal inflow gradient λc', value: sol.y[n][2], unit: '-' },
        { key: 'inflow_gradient_ratio', label: 'Fore–aft gradient λc/λ0', value: sol.y[n][2] / (sol.y[n][0] || 1e-12), unit: '-' },
        { key: 'dT_final_pct', label: 'Final thrust change', value: 100 * (ctEnd / ctOf(th0, y0[0]) - 1), unit: '%' },
        { key: 'theta_trim_deg', label: 'Initial collective at 75% radius', value: N.deg(th0), unit: 'deg' },
      ],
      plots: [
        { type: 'line', title: 'Thrust response to the control step', xlabel: 'Time [s]', ylabel: 'Thrust per rotor [N]', series: [{ name: 'With dynamic inflow', x: pick(tt), y: pick(ct).map((v) => v * p.qA) }, { name: 'Quasi-static inflow (final value)', x: [tt[0], tt[n]], y: [ctEnd * p.qA, ctEnd * p.qA], style: 'dash' }, { name: 'Frozen inflow (initial jump)', x: [tt[0], tt[n]], y: [ctQS * p.qA, ctQS * p.qA], style: 'dash' }] },
        { type: 'line', title: 'Inflow states', xlabel: 'Time [s]', ylabel: 'Inflow ratio [-]', series: [{ name: 'Uniform λ0', x: pick(tt), y: pick(sol.y.map((y) => y[0])) }, { name: 'Lateral λs', x: pick(tt), y: pick(sol.y.map((y) => y[1])) }, { name: 'Longitudinal λc', x: pick(tt), y: pick(sol.y.map((y) => y[2])) }] },
      ],
      outputs: { lambda0_initial: y0[0], lambda_h: Math.sqrt(CT0 / 2), CT_initial: ctOf(th0, y0[0]), CT_req: CT0 },
      warnings: mu > 0.02 && Math.abs(lc + l00) < 0.3 * Math.sqrt(CT0 / 2) && lc < 0 ? ['Descent with a small net flow through the disk: the momentum-based mass-flow parameter is unreliable near the vortex-ring state.'] : [],
      models: ['Pitt–Peters three-state dynamic inflow with nonlinear mass-flow parameters', 'Quasi-steady linear blade-element thrust', 'RK4 time integration'],
      assumptions: ['Rigid, non-flapping blades: thrust responds instantly to pitch and inflow', 'Apparent-mass term 8/(3π) for the uniform state (uncorrected form)', 'Only thrust and an optional pitching-moment step force the wake; roll moment is zero', 'Higher-harmonic and radial inflow states (Peters–He) are not included'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [75, 150, 300, 600], metric: 'tau_inflow_s' },
  verify() {
    const i = { mass_kg: 5000, n_rotors: 1, R: 6, Nb: 4, chord: 0.4, twist_deg: -8, rpm: 300, cla: 5.7, V: 0, alpha_d_deg: 0, d_theta_deg: 0.02, dCM: 0, t_end_rev: 10, nSteps: 2000, alt_m: 0, dISA: 0 };
    const r = inflow.run(i), o = N.kv(r), f = inflow.run({ ...i, V: 40, d_theta_deg: 0.5 }), fo = N.kv(f), mu = 40 / (rotorOf({ ...i, cd0: 0 }).Vt);
    const ctF = f.outputs.CT_initial * (1 + fo.dT_final_pct / 100);
    return [
      N.check('Steady hover inflow λ0 = sqrt(CT/2)', r.outputs.lambda0_initial, r.outputs.lambda_h, 1e-6, 'Momentum theory (static limit of the Pitt–Peters model)'),
      N.check('Trim collective reproduces the required thrust', r.outputs.CT_initial, r.outputs.CT_req, 1e-6, 'Blade-element thrust equation'),
      N.check('Small-step time constant equals m11/(4λh + σa/4)/Ω', o.tau_inflow_s, o.tau_hover_lin_s, 5e-3, 'Linearisation of the coupled blade-element / inflow equation'),
      N.check('Forward-flight steady state satisfies Glauert: λ0·sqrt(μ² + λ0²) = CT/2', fo.lambda0_final * Math.hypot(mu, fo.lambda0_final), ctF / 2, 1e-4, 'Glauert momentum theory'),
      N.check('Steady fore–aft gradient λc/λ0 = (15π/32)·tan(χ/2)', fo.inflow_gradient_ratio, ((15 * Math.PI) / 32) * Math.tan(Math.atan2(mu, fo.lambda0_final) / 2), 1e-3, 'Static gain of the Pitt–Peters L-matrix'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs;
    return [{ severity: 'info', title: 'Inflow lag shapes the heave response', detail: `The wake needs about ${fmt(o.tau_inflow_s, 2)} s (${fmt(o.tau_inflow_rev, 1)} rev) to adjust; thrust first overshoots by ${fmt(100 * o.thrust_overshoot, 0)}% of the final change.`, action: 'Include this inflow state in the flight-control model (Suite 16) when the collective or thrust loop bandwidth approaches 1/τ; ignoring it over-predicts heave damping and gain margin.', basis: 'Pitt–Peters dynamic inflow' }];
  },
};

// ---- (e) prescribed wake, Biot–Savart, blade–vortex interaction -----------------------------
/** Induced velocity at P from a straight vortex segment A→B of unit circulation with a Scully-type core rc. */
function bsSeg(P, A, B, rc2) {
  const r1 = [P[0] - A[0], P[1] - A[1], P[2] - A[2]], r2 = [P[0] - B[0], P[1] - B[1], P[2] - B[2]], cr = N.cross(r1, r2), n1 = N.norm(r1), n2 = N.norm(r2), c2 = N.dot(cr, cr);
  const L2 = (B[0] - A[0]) ** 2 + (B[1] - A[1]) ** 2 + (B[2] - A[2]) ** 2, den = n1 * n2 * (n1 * n2 + N.dot(r1, r2));
  if (den < 1e-14 || L2 < 1e-20) return [0, 0, 0];
  const h2 = c2 / L2, k = ((n1 + n2) / (4 * Math.PI * den)) * (h2 / (h2 + rc2));
  return [k * cr[0], k * cr[1], k * cr[2]];
}
/** Tip-vortex point (lengths in rotor radii) released from a blade now at azimuth psi0, at wake age w. */
function tipPoint(wk, psi0, w) {
  const az = psi0 - w;
  if (wk.hover) { const r = wk.A + (1 - wk.A) * Math.exp(-wk.Lam * w), w1 = TWO_PI / wk.Nb; return [r * Math.cos(az), r * Math.sin(az), w <= w1 ? wk.k1 * w : wk.k1 * w1 + wk.k2 * (w - w1)]; }
  return [Math.cos(az) + wk.mu * w, Math.sin(az), -wk.lam * w];
}
function wakeModel(p, i, gamDesc) {
  const W = (i.mass_kg * G0) / i.n_rotors, V = i.V, D = (0.5 * p.rho * V * V * i.f_plate) / i.n_rotors, aT = Math.atan2(D - W * Math.sin(gamDesc), W * Math.cos(gamDesc)), T = Math.hypot(W * Math.cos(gamDesc), D - W * Math.sin(gamDesc));
  const CT = T / p.qA, aS = aT, mu = (V * Math.cos(aS)) / p.Vt, lc = (V * Math.sin(aS)) / p.Vt, li = lamI(CT, mu, lc), tw = i.twist_deg;
  return { hover: mu < 0.03, Nb: p.Nb, mu, lam: lc + li, li, CT, aS, A: 0.78, Lam: 0.145 + 27 * CT, k1: -0.25 * (CT / p.sigma + 0.001 * tw), k2: -(1.41 + 0.0141 * tw) * Math.sqrt(CT / 2), G: (TWO_PI * CT) / p.Nb };
}
/** Plan-view crossings of the tip vortices with a blade over one revolution; returns [{psi, r, miss, age}] (miss in radii, + above the blade). */
function bviScan(wk, a0, nPsi = 72, revs = 3) {
  const hits = [], dA = N.rad(5), nA = Math.round((revs * TWO_PI) / dA);
  for (let kb = 0; kb < nPsi; kb++) {
    const pb = (TWO_PI * kb) / nPsi, bx = Math.cos(pb), by = Math.sin(pb);
    for (let j = 0; j < wk.Nb; j++) {
      const psi0 = pb + (TWO_PI * j) / wk.Nb; let P = tipPoint(wk, psi0, 6 * dA);
      for (let k = 7; k <= nA; k++) {
        const Q = tipPoint(wk, psi0, k * dA), dx = Q[0] - P[0], dy = Q[1] - P[1], det = dx * -by + dy * bx;
        if (Math.abs(det) > 1e-12) { const t = (-P[0] * -by + -P[1] * bx) / det, u = (dx * -P[1] - dy * -P[0]) / det; if (t >= 0 && t < 1 && u > 0.3 && u <= 1) hits.push({ psi: pb, r: u, miss: P[2] + t * (Q[2] - P[2]) - u * a0, age: (k - 1 + t) * dA }); }
        P = Q;
      }
    }
  }
  return hits;
}

const wake = {
  id: 'wake', title: 'Tip-vortex wake geometry, induced velocity and blade–vortex interaction', fidelity: 'reduced-order',
  summary: 'Lays out the tip vortices with a prescribed wake (Landgrebe in hover, a skewed helix in forward flight), computes the velocity they induce at the disk by the Biot–Savart law, and locates where vortices pass close to the blades in descent.',
  equations: ['Biot–Savart law', 'Vortex filament equations', 'Kelvin circulation theorem (constant filament strength)', 'Prescribed-wake model (Landgrebe)', 'Rotor inflow equations'],
  applicable: rotorOnly,
  inputs: [...MASS, ...ROTOR.filter((f) => !['cd0', 'clmax', 'cla'].includes(f.key)),
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 35, min: 0, max: 120, group: 'Flight' },
    { key: 'descent_deg', label: 'Descent angle', unit: 'deg', default: 6, min: -10, max: 20, group: 'Flight', help: 'Approach descents of 3–9° bring the wake back into the disk' },
    { key: 'f_plate', label: 'Equivalent flat-plate drag area', unit: 'm²', default: 3.35, min: 0, group: 'Aircraft' },
    { key: 'coning_deg', label: 'Coning angle', unit: 'deg', default: 4, min: 0, max: 12, group: 'Rotor' },
    { key: 'core_c', label: 'Vortex core radius / chord', unit: '-', default: 0.2, min: 0.02, max: 1, group: 'Calibration' },
    { key: 'turns', label: 'Wake revolutions retained', unit: '', default: 6, min: 2, max: 20, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nSeg', label: 'Vortex segments per revolution', unit: '', default: 36, min: 8, max: 360, step: 1, discrete: true, group: 'Numerics' }, ...ATM],
  defaults: (c, up, d) => { const b = baseDefaults(c, up, d); delete b.cd0; delete b.cla; return { ...b, V: d.isRotary ? 0.5 * (c.mission.cruise_V_ms || c.flight.V_ms) : Math.min(c.flight.V_ms, 0.12 * d.v_tip), f_plate: c.rotor.flat_plate_m2 || undefined, coning_deg: up.rotorcraft?.coning_deg }; },
  run(i) {
    const p = rotorOf({ ...i, cla: 5.73, cd0: 0 }), wk = wakeModel(p, i, N.rad(i.descent_deg)), a0 = N.rad(i.coning_deg), nSeg = Math.round(i.nSeg), nT = Math.round(i.turns), dW = TWO_PI / nSeg, rc2 = ((i.core_c * i.chord) / p.R) ** 2;
    // Biot–Savart downwash along the fore–aft diameter, averaged over blade phase
    const xs = N.linspace(-0.92, 0.92, 24), nPh = 4, w = xs.map(() => 0);
    for (let ph = 0; ph < nPh; ph++) for (let j = 0; j < wk.Nb; j++) {
      const psi0 = (TWO_PI * (j + ph / nPh)) / wk.Nb; let A = tipPoint(wk, psi0, 0);
      for (let k = 1; k <= nSeg * nT; k++) { const B = tipPoint(wk, psi0, k * dW); xs.forEach((x, m) => { w[m] -= (wk.G * bsSeg([x, 0, 0], A, B, rc2)[2]) / nPh; }); A = B; }
    }
    const inner = xs.map((x, m) => (Math.abs(x) < 0.7 ? w[m] : NaN)).filter(Number.isFinite), lamBS = N.mean(inner), front = N.mean(xs.map((x, m) => (x < -0.3 && x > -0.8 ? w[m] : NaN)).filter(Number.isFinite)), rear = N.mean(xs.map((x, m) => (x > 0.3 && x < 0.8 ? w[m] : NaN)).filter(Number.isFinite));
    // blade–vortex interactions
    const hits = wk.hover ? [] : bviScan(wk, a0), cR = i.chord / p.R, near = hits.filter((h) => Math.abs(h.miss) < cR), best = hits.length ? hits[N.argmin(hits.map((h) => Math.abs(h.miss)))] : null;
    const degs = N.linspace(-2, 14, 17), sweep = degs.map((dg) => { const q = wakeModel(p, i, N.rad(dg)); if (q.hover) return NaN; const h = bviScan(q, a0, 36, 2).filter((v) => v.psi < Math.PI); return h.length ? N.amin(h.map((v) => Math.abs(v.miss))) / cR : NaN; });
    const fin = sweep.map((v) => (Number.isFinite(v) ? v : Infinity)), kc = N.argmin(fin), warnings = [];
    if (wk.hover) warnings.push('Hover or near-hover: the Landgrebe prescribed wake is used and blade–vortex interaction in descent is not evaluated. Set an airspeed for the descent study.');
    if (wk.lam < 0) warnings.push('Net flow through the disk is upward: the wake is convected above the rotor (steep descent). A rigid wake is a poor model here.');
    // plan and side views at blade azimuth 0
    const top = [], side = []; for (let j = 0; j < wk.Nb; j++) { for (let k = 0; k <= Math.min(nT, 3) * 36; k++) { const P = tipPoint(wk, (TWO_PI * j) / wk.Nb, (k * TWO_PI) / 36); top.push(P); if (j === 0) side.push(P); } top.push([NaN, NaN, NaN]); }
    const circ = N.linspace(0, TWO_PI, 73);
    return {
      kpis: [
        { key: 'bvi_miss_chords', label: 'Closest blade–vortex miss distance', value: best ? Math.abs(best.miss) / cR : NaN, unit: 'chords', status: !best ? undefined : Math.abs(best.miss) / cR > 1 ? 'ok' : 'warn', note: 'Below about one chord the interaction is strong and impulsive' },
        { key: 'bvi_miss_m', label: 'Closest miss distance', value: best ? best.miss * p.R : NaN, unit: 'm', note: 'Positive: vortex above the blade' },
        { key: 'bvi_psi_deg', label: 'Azimuth of the closest interaction', value: best ? N.deg(best.psi) : NaN, unit: 'deg' },
        { key: 'bvi_r_R', label: 'Radius of the closest interaction', value: best ? best.r : NaN, unit: 'r/R' },
        { key: 'bvi_count', label: 'Interactions within one chord per revolution', value: near.length, unit: '', note: 'Counted on a 5° azimuth grid' },
        { key: 'bvi_critical_descent_deg', label: 'Descent angle of closest advancing-side interaction', value: Number.isFinite(fin[kc]) ? degs[kc] : NaN, unit: 'deg' },
        { key: 'lambda_wake', label: 'Biot–Savart mean induced inflow (|x| < 0.7R)', value: lamBS, unit: '-' },
        { key: 'lambda_momentum', label: 'Momentum-theory induced inflow', value: wk.li, unit: '-' },
        { key: 'wake_inflow_ratio', label: 'Wake / momentum inflow ratio', value: lamBS / (wk.li || 1e-12), unit: '-', note: 'Tip vortices only: the inboard vortex sheet and the truncated far wake are missing' },
        { key: 'inflow_rear_front', label: 'Rear / front downwash ratio', value: rear / (front || 1e-12), unit: '-' },
        { key: 'tip_vortex_circ_m2s', label: 'Tip-vortex circulation', value: wk.G * p.Vt * p.R, unit: 'm²/s' },
        { key: 'wake_skew_deg', label: 'Wake skew angle', value: N.deg(Math.atan2(wk.mu, Math.abs(wk.lam) + 1e-12)), unit: 'deg' },
        { key: 'vortex_spacing_R', label: 'Vertical spacing between successive vortices', value: wk.hover ? Math.abs(wk.k2) * (TWO_PI / wk.Nb) : (wk.lam * TWO_PI) / wk.Nb, unit: 'R' },
      ],
      plots: [
        { type: 'line', title: 'Tip-vortex trajectories, plan view', xlabel: 'x/R (aft +) [-]', ylabel: 'y/R (advancing side +) [-]', equalAspect: true, series: [{ name: 'Tip vortices', x: top.map((v) => v[0]), y: top.map((v) => v[1]) }, { name: 'Rotor disk', x: circ.map(Math.cos), y: circ.map(Math.sin), style: 'dash' }, { name: 'Interactions within one chord', x: near.map((h) => h.r * Math.cos(h.psi)), y: near.map((h) => h.r * Math.sin(h.psi)), style: 'points' }] },
        { type: 'line', title: 'Tip-vortex trajectory, side view (one blade)', xlabel: 'x/R (aft +) [-]', ylabel: 'z/R (up +) [-]', series: [{ name: 'Tip vortex', x: side.map((v) => v[0]), y: side.map((v) => v[2]) }, { name: 'Blade tip path (coned)', x: [-1, 0, 1], y: [a0, 0, a0], style: 'dash' }] },
        { type: 'line', title: 'Induced downwash along the fore–aft diameter', xlabel: 'x/R (aft +) [-]', ylabel: 'Induced velocity [m/s]', series: [{ name: 'Biot–Savart (tip vortices)', x: xs, y: w.map((v) => v * p.Vt) }, { name: 'Momentum theory (uniform)', x: [-1, 1], y: [wk.li * p.Vt, wk.li * p.Vt], style: 'dash' }] },
        { type: 'line', title: 'Closest advancing-side miss distance versus descent angle', xlabel: 'Descent angle [deg]', ylabel: 'Miss distance [chords]', series: [{ name: 'Minimum |miss|', x: degs, y: sweep, style: 'line+points' }], annotations: [{ y: 1, label: 'One chord' }] },
      ],
      warnings, models: [wk.hover ? 'Landgrebe prescribed hover wake (tip vortex)' : 'Rigid skewed helical wake convected by the mean inflow', 'Biot–Savart law on straight segments with a Scully-type core', 'Geometric blade–vortex miss-distance search'],
      assumptions: ['Constant-strength tip vortices from uniformly loaded blades, Γ = 2π·CT·ΩR²/Nb', 'No wake distortion, roll-up dynamics or vortex decay; inboard sheet neglected', 'Blades are straight and coned; flapping harmonics and elastic deflection are not applied to the miss distance', 'Wake truncated after the retained revolutions'],
    };
  },
  convergence: { param: 'nSeg', label: 'Vortex segments per revolution', levels: [12, 24, 48, 96], metric: 'lambda_wake' },
  calibration: { params: [{ key: 'core_c', min: 0.02, max: 1 }, { key: 'coning_deg', min: 0, max: 12 }], sweep: 'descent_deg', target: 'bvi_miss_chords', note: 'Measured vortex positions or miss distances from wake imaging (e.g. PIV or laser light sheet) versus descent angle.' },
  verify() {
    // Polygonal vortex ring against the exact axial velocity of a circular ring.
    const n = 360, z = 0.6; let w = 0;
    for (let k = 0; k < n; k++) w += bsSeg([0, 0, z], [Math.cos((TWO_PI * k) / n), Math.sin((TWO_PI * k) / n), 0], [Math.cos((TWO_PI * (k + 1)) / n), Math.sin((TWO_PI * (k + 1)) / n), 0], 0)[2];
    const inf = bsSeg([0, 1, 0], [-1e6, 0, 0], [1e6, 0, 0], 0)[2];
    const wk = { hover: false, Nb: 2, mu: 0.2, lam: 0.03 }, P = tipPoint(wk, 1, 2);
    return [
      N.check('Vortex ring on its axis: w = ΓR²/(2(R² + z²)^1.5)', w, 1 / (2 * (1 + z * z) ** 1.5), 1e-4, 'Exact Biot–Savart integral for a circular ring'),
      N.check('Infinite straight vortex: v = Γ/(2πh)', inf, 1 / TWO_PI, 1e-9, 'Biot–Savart law'),
      N.check('Skewed helix convects at (μ, −λ) per radian of wake age', (P[0] - Math.cos(-1)) / 2 + (-P[2] / 2), 0.2 + 0.03, 1e-12, 'Rigid-wake kinematics'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (Number.isFinite(o.bvi_miss_chords) && o.bvi_miss_chords < 1) out.push({ severity: 'warn', title: 'Close blade–vortex interaction at this descent condition', detail: `A tip vortex passes ${o.bvi_miss_chords.toFixed(2)} chords from a blade near ψ = ${o.bvi_psi_deg.toFixed(0)}°, r/R = ${o.bvi_r_R.toFixed(2)}; ${o.bvi_count} close encounters per revolution.`, action: 'For approach noise abatement fly steeper or shallower than the critical descent angle, or change speed, so the wake clears the disk; quantify the benefit in Suite 11.', basis: 'Geometric miss distance on a rigid wake' });
    if (Number.isFinite(o.bvi_critical_descent_deg)) out.push({ severity: 'info', title: 'BVI-critical descent angle', detail: `At ${i.V.toFixed(0)} m/s the advancing-side miss distance is smallest near a ${o.bvi_critical_descent_deg.toFixed(0)}° descent.`, action: 'Use this to shape low-noise approach profiles; confirm with a free-wake analysis before publishing procedures.', basis: 'Wake convected by μ and λ' });
    return out;
  },
};

// ---- (f) dynamic stall ----------------------------------------------------------------------
/** Leishman–Beddoes-type indicial model in discrete (exponential-recurrence) form. Returns the last cycle. */
function lbLoop(o) {
  const b2 = 1 - o.M * o.M, cna = o.cna / Math.sqrt(b2), nS = o.nStep, ds = TWO_PI / (o.k * nS), a1r = N.rad(o.alpha1), s1 = N.rad(o.S1), s2 = N.rad(o.S2);
  const fsep = (a) => (!o.sep ? 1 : Math.abs(a) <= a1r ? 1 - 0.3 * Math.exp((Math.abs(a) - a1r) / s1) : 0.04 + 0.66 * Math.exp((a1r - Math.abs(a)) / s2));
  const e1 = Math.exp(-o.b1 * b2 * ds), e2 = Math.exp(-o.b2 * b2 * ds), ep = Math.exp(-ds / o.Tp), ev = Math.exp(-ds / o.Tv);
  let X = 0, Y = 0, Dp = 0, Df = 0, CNv = 0, tv = 0, aqP = null, cnpP = 0, fpP = 1, cvP = 0, alP = o.a0, onset = NaN, reatt = NaN;
  const out = { al: [], cn: [], cnc: [], f: [], cnv: [] };
  for (let n = 0; n <= o.nCyc * nS; n++) {
    const s = n * ds, al = o.a0 + o.a1 * Math.sin(o.k * s), d1 = o.a1 * o.k * Math.cos(o.k * s), d2 = -o.a1 * o.k * o.k * Math.sin(o.k * s), aq = al + (o.rate ? d1 : 0); // incidence at 3/4 chord for pitch about 1/4 chord
    if (aqP === null) aqP = aq;
    const dA = aq - aqP; X = X * e1 + o.A1 * dA * Math.sqrt(e1); Y = Y * e2 + o.A2 * dA * Math.sqrt(e2);
    const aE = aq - X - Y, cnc = cna * aE, cni = o.nc ? Math.PI * d1 + (Math.PI / 2) * d2 : 0, cnp = cnc + cni;
    Dp = Dp * ep + (cnp - cnpP) * Math.sqrt(ep); const cnPrime = cnp - Dp, fp = fsep(cnPrime / cna);
    const Tf = tv > 0 && tv <= o.Tvl ? o.Tf / 2 : o.Tf, ef = Math.exp(-ds / Tf); Df = Df * ef + (fp - fpP) * Math.sqrt(ef);
    const fpp = N.clamp(fp - Df, 0, 1), Kf = ((1 + Math.sqrt(fpp)) / 2) ** 2, cv = cnc * (1 - Kf), up = Math.abs(al) > Math.abs(alP);
    if (o.sep && Math.abs(cnPrime) > o.CN1) { if (tv === 0 && n > (o.nCyc - 1) * nS && Number.isNaN(onset)) onset = al; tv += ds; } else if (tv > 0 && !up) { if (fpp > 0.7) { tv = 0; if (n > (o.nCyc - 1) * nS && Number.isNaN(reatt)) reatt = al; } } else if (up && Math.abs(cnPrime) < o.CN1) tv = 0;
    CNv = tv > 0 && tv <= o.Tvl ? CNv * ev + (cv - cvP) * Math.sqrt(ev) : CNv * ev * ev; // vortex lift accumulates while the vortex is over the chord, then decays twice as fast
    if (n >= (o.nCyc - 1) * nS) { out.al.push(al); out.cn.push(cna * Kf * aE + cni + CNv); out.cnc.push(cnc); out.f.push(fpp); out.cnv.push(CNv); }
    aqP = aq; cnpP = cnp; fpP = fp; cvP = cv; alP = al;
  }
  return { ...out, cna, onset, reatt, stat: (a) => cna * ((1 + Math.sqrt(fsep(a))) / 2) ** 2 * a, ds };
}

const dynstall = {
  id: 'dynstall', title: 'Dynamic stall of a pitching blade section', fidelity: 'reduced-order',
  summary: 'A Leishman–Beddoes-type state model of an aerofoil oscillating in pitch: lift lags in attached flow, stall is delayed beyond the static angle, a leading-edge vortex adds lift, and reattachment is late — giving the characteristic lift hysteresis loop of the retreating blade.',
  equations: ['Dynamic stall state equations (Leishman–Beddoes type)', 'Indicial (Wagner-type) unsteady lift', 'Kirchhoff–Helmholtz trailing-edge separation', 'Blade feathering equations (prescribed pitch)'],
  applicable: rotorOnly,
  inputs: [
    { key: 'alpha_mean_deg', label: 'Mean incidence', unit: 'deg', default: 10, min: -5, max: 25, group: 'Motion' },
    { key: 'alpha_amp_deg', label: 'Incidence amplitude', unit: 'deg', default: 8, min: 0.1, max: 15, group: 'Motion' },
    { key: 'k', label: 'Reduced frequency k = ωc/2V', unit: '-', default: 0.08, min: 0.005, max: 0.5, group: 'Motion', help: '1/rev pitching at the retreating blade is typically 0.05–0.15' },
    { key: 'M', label: 'Section Mach number', unit: '-', default: 0.3, min: 0, max: 0.7, group: 'Flow' },
    { key: 'cna', label: 'Incompressible normal-force slope', unit: '1/rad', default: 6.28, min: 4, max: 7, group: 'Section' },
    { key: 'alpha1_deg', label: 'Static break angle α1', unit: 'deg', default: 15.25, min: 6, max: 25, group: 'Section', help: 'Incidence where the separation point reaches 70% chord' },
    { key: 'S1_deg', label: 'Separation fit S1', unit: 'deg', default: 3.0, min: 0.3, max: 10, group: 'Section' },
    { key: 'S2_deg', label: 'Separation fit S2', unit: 'deg', default: 2.3, min: 0.3, max: 10, group: 'Section' },
    { key: 'CN1', label: 'Critical normal force for leading-edge separation', unit: '-', default: 1.45, min: 0.6, max: 2.5, group: 'Section' },
    { key: 'Tp', label: 'Pressure lag Tp', unit: 'semi-chords', default: 1.7, min: 0.1, max: 6, group: 'Time constants' },
    { key: 'Tf', label: 'Separation lag Tf', unit: 'semi-chords', default: 3.0, min: 0.1, max: 12, group: 'Time constants' },
    { key: 'Tv', label: 'Vortex decay Tv', unit: 'semi-chords', default: 6.0, min: 0.5, max: 15, group: 'Time constants' },
    { key: 'Tvl', label: 'Vortex passage time Tvl', unit: 'semi-chords', default: 7.0, min: 2, max: 15, group: 'Time constants' },
    { key: 'nStep', label: 'Steps per cycle', unit: '', default: 240, min: 40, max: 5000, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nCyc', label: 'Cycles', unit: '', default: 4, min: 2, max: 20, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => {
    // 1/rev pitching seen at 75% radius on the retreating side
    const V = Math.max(0.2 * d.v_tip, 0.75 * d.v_tip - ffSpeed(c, d)), a = isa(c.atm.alt_m, c.atm.dISA_K).a;
    return { k: N.clamp((d.omega * c.rotor.chord_m) / (2 * V), 0.005, 0.5), M: N.clamp(V / a, 0, 0.7) };
  },
  run(i) {
    const o = { a0: N.rad(i.alpha_mean_deg), a1: N.rad(i.alpha_amp_deg), k: i.k, M: i.M, cna: i.cna, alpha1: i.alpha1_deg, S1: i.S1_deg, S2: i.S2_deg, CN1: i.CN1, Tp: i.Tp, Tf: i.Tf, Tv: i.Tv, Tvl: i.Tvl, A1: 0.3, b1: 0.14, A2: 0.7, b2: 0.53, nStep: Math.round(i.nStep), nCyc: Math.round(i.nCyc), sep: true, rate: true, nc: true };
    const r = lbLoop(o), att = lbLoop({ ...o, sep: false }), as = N.linspace(Math.min(0, o.a0 - o.a1), o.a0 + o.a1 + N.rad(4), 80), st = as.map(r.stat), cnMaxS = N.amax(N.linspace(0, N.rad(30), 300).map(r.stat)), kmax = N.argmax(r.cn), cnMax = r.cn[kmax];
    let area = 0; for (let k = 1; k < r.al.length; k++) area += 0.5 * (r.cn[k] + r.cn[k - 1]) * (r.al[k] - r.al[k - 1]);
    const stalls = Number.isFinite(r.onset), warnings = [];
    if (!stalls) warnings.push('The section stays below the leading-edge separation criterion over the whole cycle: the loop shows unsteady attached-flow and trailing-edge separation effects only.');
    if (i.M > 0.5) warnings.push('Above Mach 0.5 the time constants and separation parameters change markedly and shock-induced separation appears; recalibrate before use.');
    const ds = Math.max(1, Math.floor(r.al.length / 300)), pick = (a) => a.filter((_, k) => k % ds === 0);
    return {
      kpis: [
        { key: 'CN_max_dyn', label: 'Peak dynamic normal-force coefficient', value: cnMax, unit: '-' },
        { key: 'CN_max_static', label: 'Static maximum (same model)', value: cnMaxS, unit: '-' },
        { key: 'lift_overshoot', label: 'Dynamic / static maximum lift', value: cnMax / cnMaxS, unit: '-', note: 'Typically 1.1–1.8, growing with reduced frequency and stall penetration' },
        { key: 'alpha_CNmax_deg', label: 'Incidence at peak lift', value: N.deg(r.al[kmax]), unit: 'deg' },
        { key: 'stall_onset_deg', label: 'Dynamic stall onset incidence (upstroke)', value: N.deg(r.onset), unit: 'deg', note: `Static leading-edge criterion is met at ${N.deg(i.CN1 / r.cna).toFixed(1)}°` },
        { key: 'stall_delay_deg', label: 'Stall delay beyond the static critical incidence CN1/CNα', value: N.deg(r.onset - i.CN1 / r.cna), unit: 'deg' },
        { key: 'reattach_deg', label: 'Reattachment incidence (downstroke)', value: N.deg(r.reatt), unit: 'deg' },
        { key: 'hysteresis_area', label: 'Lift hysteresis loop area ∮CN dα', value: area, unit: 'rad', note: 'Negative (clockwise) loops are typical of stall' },
        { key: 'f_min', label: 'Minimum attached-flow fraction', value: N.amin(r.f), unit: '-' },
        { key: 'CN_vortex_max', label: 'Peak vortex-lift increment', value: N.amax(r.cnv), unit: '-' },
      ],
      plots: [
        { type: 'line', title: 'Normal-force hysteresis loop', xlabel: 'Incidence α [deg]', ylabel: 'Normal-force coefficient CN [-]', series: [{ name: 'Dynamic (with stall)', x: pick(r.al).map(N.deg), y: pick(r.cn) }, { name: 'Unsteady attached flow', x: pick(att.al).map(N.deg), y: pick(att.cn), style: 'dash' }, { name: 'Static', x: as.map(N.deg), y: st, style: 'dash' }], annotations: Number.isFinite(r.onset) ? [{ x: N.deg(r.onset), label: 'Onset' }] : [] },
        { type: 'line', title: 'Separation point and vortex lift over the cycle', xlabel: 'Cycle phase [deg]', ylabel: 'Attached fraction f / vortex CN [-]', series: [{ name: 'Attached-flow fraction f', x: pick(N.linspace(0, 360, r.al.length)), y: pick(r.f) }, { name: 'Vortex-lift increment', x: pick(N.linspace(0, 360, r.al.length)), y: pick(r.cnv) }] },
      ],
      warnings, models: ['Leishman–Beddoes-type dynamic stall: indicial attached flow, pressure lag, Kirchhoff separation with lagged separation point, vortex-lift accumulation and decay', 'Two-term exponential indicial lift function (A1 = 0.3, b1 = 0.14, A2 = 0.7, b2 = 0.53) with β² compressibility scaling', 'Incompressible apparent-mass (non-circulatory) terms'],
      assumptions: ['Simplified against the published model: non-circulatory loads use the incompressible form, pitching moment and chord force are not computed, and the time-constant switching logic is reduced to a halved Tf during vortex passage', 'Default constants are representative of a NACA 0012-class section near Mach 0.3 and must be recalibrated for other sections and Mach numbers', 'Harmonic pitch about the quarter chord at constant free-stream speed'],
    };
  },
  convergence: { param: 'nStep', label: 'Steps per cycle', levels: [60, 120, 240, 480], metric: 'CN_max_dyn' },
  calibration: { params: [{ key: 'Tp', min: 0.1, max: 6 }, { key: 'Tf', min: 0.1, max: 12 }, { key: 'Tv', min: 0.5, max: 15 }, { key: 'CN1', min: 0.6, max: 2.5 }, { key: 'alpha1_deg', min: 6, max: 25 }], sweep: 'k', target: 'CN_max_dyn', note: 'Oscillating-aerofoil wind-tunnel loops (peak CN versus reduced frequency) for the actual section and Mach number.' },
  verify() {
    // Attached flow, circulatory part only: frequency response of the indicial kernel, and Theodorsen's function with R. T. Jones' coefficients.
    const resp = (A1, b1, A2, b2, k) => {
      const r = lbLoop({ a0: 0, a1: 0.02, k, M: 0, cna: TWO_PI, alpha1: 90, S1: 3, S2: 3, CN1: 99, Tp: 1.7, Tf: 3, Tv: 6, Tvl: 7, A1, b1, A2, b2, nStep: 2000, nCyc: 12, sep: false, rate: false, nc: false });
      let re = 0, im = 0; const n = r.al.length - 1; for (let j = 0; j < n; j++) { const ph = (TWO_PI * j) / n; re += r.cnc[j] * Math.sin(ph); im += r.cnc[j] * Math.cos(ph); }
      return [(2 * re) / n / (TWO_PI * 0.02), (2 * im) / n / (TWO_PI * 0.02)];
    };
    const k = 0.2, C = (A1, b1, A2, b2) => { const t = (A, b) => [(A * k * k) / (k * k + b * b), (A * k * b) / (k * k + b * b)], u = t(A1, b1), v = t(A2, b2); return [1 - u[0] - v[0], -u[1] - v[1]]; };
    const lb = resp(0.3, 0.14, 0.7, 0.53, k), ex = C(0.3, 0.14, 0.7, 0.53), jo = resp(0.165, 0.0455, 0.335, 0.3, k);
    const st = lbLoop({ a0: 0, a1: 0.1, k: 0.1, M: 0, cna: TWO_PI, alpha1: 15, S1: 3, S2: 2.3, CN1: 1.45, Tp: 1.7, Tf: 3, Tv: 6, Tvl: 7, A1: 0.3, b1: 0.14, A2: 0.7, b2: 0.53, nStep: 40, nCyc: 2, sep: true, rate: true, nc: true });
    return [
      N.check('Circulatory lift, in-phase part of the lift-deficiency function', lb[0], ex[0], 1e-4, 'Exact frequency response of the exponential indicial function'),
      N.check('Circulatory lift, quadrature part', lb[1], ex[1], 1e-3, 'Exact frequency response of the exponential indicial function'),
      N.check('Theodorsen F(k = 0.2) with R. T. Jones coefficients', jo[0], 0.7276, 0.025, 'Theodorsen (NACA Report 496): C(0.2) = 0.7276 − 0.1886i; Jones two-term fit is accurate to about 2%'),
      N.check('Theodorsen G(k = 0.2) with R. T. Jones coefficients', jo[1], -0.1886, 0.025, 'Theodorsen (NACA Report 496)'),
      N.check('Static Kirchhoff lift at the break angle: CN = CNα·((1+√0.7)/2)²·α1', st.stat(N.rad(15)), TWO_PI * ((1 + Math.sqrt(0.7)) / 2) ** 2 * N.rad(15), 1e-12, 'Kirchhoff–Helmholtz flat-plate separation'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (Number.isFinite(o.stall_onset_deg)) out.push({ severity: 'advise', title: 'The section enters dynamic stall each cycle', detail: `Onset at ${o.stall_onset_deg.toFixed(1)}° (${o.stall_delay_deg.toFixed(1)}° later than in steady flow), peak CN ${o.CN_max_dyn.toFixed(2)} = ${o.lift_overshoot.toFixed(2)} × static, reattachment near ${fmt(o.reattach_deg)}°.`, action: 'The extra lift is transient and comes with large nose-down pitching moments (not modelled here) that drive pitch-link loads and stall flutter. Limit blade loading or speed, or select a section with a higher dynamic-stall onset; carry the oscillatory loads to Suite 9 for fatigue.', basis: 'Leishman–Beddoes leading-edge separation criterion' });
    else out.push({ severity: 'info', title: 'No leading-edge dynamic stall at this condition', detail: `Peak CN ${o.CN_max_dyn.toFixed(2)} stays below the critical value ${i.CN1.toFixed(2)}.`, action: 'Increase the mean incidence or amplitude to explore the manoeuvre and high-speed corners of the envelope.', basis: 'Critical normal-force criterion CN1' });
    return out;
  },
};

// ---- (g) autorotation, vortex-ring state ----------------------------------------------------
/** Induced velocity ratio vi/vh in vertical flight, x = Vc/vh (climb +): momentum branches outside −2 < x < 0, empirical quartic inside. */
const viAxial = (x) => (x >= 0 ? -x / 2 + Math.sqrt(x * x / 4 + 1) : x <= -2 ? -x / 2 - Math.sqrt(x * x / 4 - 1) : 1 - 1.125 * x - 1.372 * x * x - 1.718 * x ** 3 - 0.655 * x ** 4);

const autorot = {
  id: 'autorot', title: 'Autorotation, rotor energy and vortex-ring state', fidelity: 'analytical',
  summary: 'Power-off descent rate versus forward speed from the energy balance, the stored rotor energy available for the entry and the landing flare, and the descent band in which the rotor works in its own recirculating wake.',
  equations: ['Autorotation energy balance', 'Rotor momentum theory in axial and forward flight', 'Induced velocity equations (empirical descent fit)', 'Rotor angular-momentum equation'],
  applicable: rotorOnly,
  inputs: [...MASS, ...ROTOR.filter((f) => !['twist_deg', 'cla', 'clmax'].includes(f.key)),
    { key: 'blade_mass', label: 'Blade mass', unit: 'kg', default: 116, min: 0.001, group: 'Rotor' },
    { key: 'f_plate', label: 'Equivalent flat-plate drag area', unit: 'm²', default: 3.35, min: 0, group: 'Aircraft' },
    { key: 'omega_min_frac', label: 'Minimum usable rotor speed / normal', unit: '-', default: 0.75, min: 0.4, max: 0.95, group: 'Limits', help: 'Below this the blades stall or cone excessively' },
    { key: 'vrs_lo', label: 'Vortex-ring band, lower descent rate / vh', unit: '-', default: 0.5, min: 0.2, max: 1, group: 'Limits', help: 'Textbook guideline: roughly 0.5–1.5 vh at low forward speed' },
    { key: 'vrs_hi', label: 'Vortex-ring band, upper descent rate / vh', unit: '-', default: 1.5, min: 1, max: 2, group: 'Limits' },
    KAPPA, ...ATM],
  defaults: (c, up, d) => { const b = baseDefaults(c, up, d); delete b.twist_deg; delete b.cla; return { ...b, blade_mass: c.rotor.blade_mass_kg || undefined, f_plate: c.rotor.flat_plate_m2 || undefined }; },
  run(i, ctx) {
    const p = rotorOf({ ...i, cla: 5.73 }), W = i.mass_kg * G0, T = W / i.n_rotors, vh = Math.sqrt(T / (2 * p.rho * p.A)), n = i.n_rotors, warnings = [];
    const P0 = (mu) => ((n * p.qA * p.Vt * p.sigma * i.cd0) / 8) * (1 + 4.65 * mu * mu);
    // forward autorotation: W·Vd = κ·W·vi + P0 + ½ρV³f with Glauert vi for flow up through the disk
    const sink = (V) => {
      let Vd = (P0(V / p.Vt) + 0.5 * p.rho * V ** 3 * i.f_plate) / W + vh, vi = vh;
      for (let k = 0; k < 60; k++) { const g = (v) => v - (vh * vh) / Math.hypot(V, Vd - v), r = N.findRoot(g, 0, (vh * vh) / Math.max(V, 1e-6) + 1e-9, 30, 1e-12); vi = Number.isFinite(r) ? r : vi; const nv = (i.kappa * W * vi + P0(V / p.Vt) + 0.5 * p.rho * V ** 3 * i.f_plate) / W; if (Math.abs(nv - Vd) < 1e-10) { Vd = nv; break; } Vd = 0.5 * (Vd + nv); }
      return { Vd, vi };
    };
    // vertical autorotation: T·(Vc + κ·vi) + P0 = 0 with the empirical descent-state induced velocity
    const xv = N.findRoot((x) => W * vh * (x + i.kappa * viAxial(x)) + P0(0), -4, -0.5, 80, 1e-10), VdVert = -xv * vh;
    const Vs = N.linspace(Math.max(2 * vh, 0.08 * p.Vt), Math.min(0.42 * p.Vt, 0.92 * p.a - p.Vt), 40), sk = Vs.map((V) => sink(V).Vd), kmin = N.argmin(sk), kgl = N.argmax(Vs.map((V, k) => V / sk[k]));
    // stored energy
    const Ib = (i.blade_mass * p.R ** 2) / 3, IR = n * p.Nb * Ib, KE = 0.5 * IR * p.Om ** 2, DL = T / p.A, AI = KE / (W * DL), usable = KE * (1 - i.omega_min_frac ** 2);
    const Ph = n * (i.kappa * T * vh) + P0(0), tEq = usable / Ph, Q0 = Ph / p.Om, tau = (IR * p.Om) / Q0, tDecay = tau * (1 / i.omega_min_frac - 1);
    const flare = usable / (0.5 * i.mass_kg * sk[kmin] ** 2), tt = N.linspace(0, Math.max(2, 2.5 * tDecay), 60);
    if (i.Nb * i.blade_mass < 0.002 * i.mass_kg || n > 2) warnings.push('Small fixed-pitch multirotors and distributed lift rotors cannot usually enter or sustain autorotation; the figures show stored energy and descent physics only. Safe recovery relies on redundancy or a ballistic parachute.');
    if (!Number.isFinite(xv)) warnings.push('No vertical autorotation solution was found in the range of the empirical induced-velocity fit.');
    const AIft = AI / 0.006366, ell = N.linspace(0, Math.PI, 40), heli = n <= 2; // autorotation guide values are helicopter practice; distributed-lift vehicles recover by redundancy
    return {
      kpis: [
        { key: 'autorotation_index', label: 'Autorotation index I·Ω²/(2·W·DL)', value: AI, unit: 'm³/N', status: !heli ? undefined : AIft > 20 ? 'ok' : AIft > 10 ? 'warn' : 'bad', note: `${AIft.toFixed(1)} ft³/lb; helicopter guide values often quoted are about 20 (single engine) and 10 (multi-engine)${heli ? '' : '. Not a criterion for multirotor or distributed-lift vehicles'}` },
        { key: 'rotor_KE_J', label: 'Rotor kinetic energy', value: KE, unit: 'J' },
        { key: 'hover_time_equiv_s', label: 'Equivalent hover time from usable rotor energy', value: tEq, unit: 's', status: !heli ? undefined : tEq > 1.5 ? 'ok' : tEq > 0.8 ? 'warn' : 'bad' },
        { key: 'rotor_decay_time_s', label: 'Time for rotor speed to fall to the minimum (no pilot action, hover)', value: tDecay, unit: 's', status: !heli ? undefined : tDecay > 1 ? 'ok' : 'warn', note: 'Available intervention time after power loss' },
        { key: 'sink_min_ms', label: 'Minimum autorotative descent rate', value: sk[kmin], unit: 'm/s' },
        { key: 'V_sink_min_ms', label: 'Speed for minimum descent rate', value: Vs[kmin], unit: 'm/s' },
        { key: 'glide_ratio', label: 'Best autorotative glide ratio', value: Vs[kgl] / sk[kgl], unit: '-' },
        { key: 'V_glide_ms', label: 'Speed for best glide', value: Vs[kgl], unit: 'm/s' },
        { key: 'sink_vertical_ms', label: 'Vertical autorotation descent rate', value: VdVert, unit: 'm/s', note: `${fmt(-xv, 2)} × hover induced velocity` },
        { key: 'flare_energy_ratio', label: 'Usable rotor energy / descent kinetic energy at minimum sink', value: flare, unit: '-', status: !heli ? undefined : flare > 1 ? 'ok' : 'warn' },
        { key: 'v_hover_induced_ms', label: 'Hover induced velocity vh', value: vh, unit: 'm/s' },
        { key: 'vrs_sink_lo_ms', label: 'Vortex-ring band begins near', value: i.vrs_lo * vh, unit: 'm/s' },
        { key: 'vrs_sink_hi_ms', label: 'Vortex-ring band ends near', value: i.vrs_hi * vh, unit: 'm/s' },
        { key: 'rotor_inertia_kgm2', label: 'Rotor polar inertia (all rotors)', value: IR, unit: 'kg m²' },
      ],
      plots: [
        { type: 'line', title: 'Autorotative descent rate versus forward speed', xlabel: 'True airspeed [m/s]', ylabel: 'Descent rate [m/s]', series: [{ name: 'Steady autorotation (energy balance)', x: Vs, y: sk }, { name: 'Vertical autorotation (empirical inflow)', x: [0], y: [VdVert], style: 'points' }], annotations: [{ x: Vs[kmin], label: 'Min sink' }, { x: Vs[kgl], label: 'Best glide' }] },
        { type: 'line', title: 'Rotor speed decay after power loss in hover (fixed collective)', xlabel: 'Time [s]', ylabel: 'Rotor speed [% of normal]', series: [{ name: 'Ω/Ω0 = 1/(1 + t/τ)', x: tt, y: tt.map((t) => 100 / (1 + t / tau)) }], annotations: [{ y: 100 * i.omega_min_frac, label: 'Minimum usable' }] },
        { type: 'line', title: 'Induced velocity in axial flight and the vortex-ring band', xlabel: 'Climb velocity Vc/vh [-]', ylabel: 'Induced velocity vi/vh [-]', series: [{ name: 'Momentum theory / empirical fit', x: N.linspace(-3, 2, 101), y: N.linspace(-3, 2, 101).map(viAxial) }, { name: 'Indicative vortex-ring region (guideline)', x: [-i.vrs_lo, -i.vrs_lo, -i.vrs_hi, -i.vrs_hi], y: [0, 2.2, 2.2, 0], style: 'dash' }], annotations: [{ x: xv, label: 'Autorotation' }] },
        { type: 'line', title: 'Indicative vortex-ring-state region', xlabel: 'Forward speed Vx/vh [-]', ylabel: 'Descent rate Vd/vh [-]', series: [{ name: 'Avoid: textbook guideline, not a prediction', x: [...ell.map((a) => Math.sin(a)), 0], y: [...ell.map((a) => 0.5 * (i.vrs_lo + i.vrs_hi) - 0.5 * (i.vrs_hi - i.vrs_lo) * Math.cos(a)), i.vrs_lo] }] },
      ],
      outputs: { helicopter_criteria: heli ? 1 : 0, energy_residual: Math.abs(W * sk[kmin] - (i.kappa * W * sink(Vs[kmin]).vi + P0(Vs[kmin] / p.Vt) + 0.5 * p.rho * Vs[kmin] ** 3 * i.f_plate)) / (W * sk[kmin]) },
      warnings, models: ['Energy-balance autorotation with Glauert induced velocity', 'Empirical quartic induced-velocity fit for −2 ≤ Vc/vh ≤ 0 (descent states)', 'Constant-collective rotor speed decay (torque ∝ Ω²)', 'Autorotation index (stored energy per unit weight and disk loading)'],
      assumptions: ['Uniform blade mass: flap inertia m·R²/3 per blade', 'Steady autorotation at normal rotor speed; entry and flare manoeuvres are not simulated', 'The vortex-ring band is a guideline region (momentum theory is invalid there), not a computed stability boundary', 'Profile power from the mean cd0 with the 1 + 4.65μ² factor', 'Minimum usable rotor speed and the vortex-ring band limits are typical textbook values'],
    };
  },
  verify() {
    const i = { mass_kg: 3000, n_rotors: 1, R: 5.5, Nb: 4, chord: 0.35, rpm: 380, cd0: 0, blade_mass: 60, f_plate: 0, omega_min_frac: 0.8, vrs_lo: 0.5, vrs_hi: 1.5, kappa: 1, alt_m: 0, dISA: 0 };
    const r = autorot.run(i), o = N.kv(r), tau = o.rotor_decay_time_s / (1 / 0.8 - 1);
    const sol = N.rk4((t, y) => [-(y[0] * y[0]) / tau], 0, [1], o.rotor_decay_time_s, 400);
    return [
      N.check('Empirical descent fit meets momentum theory at Vc/vh = −2 (vi/vh = 1)', 1 - 1.125 * -2 - 1.372 * 4 - 1.718 * -8 - 0.655 * 16, 1, 0.03, 'Windmill-brake momentum solution'),
      N.check('Ideal vertical autorotation at Vc ≈ −1.77·vh (loss-free rotor)', -o.sink_vertical_ms / o.v_hover_induced_ms, -1.775, 0.01, 'Zero of Vc + vi with the descent fit; textbook ideal autorotation 1.75–1.8 vh'),
      N.check('Rotor speed decay Ω/Ω0 = 1/(1 + t/τ) reaches the minimum at the stated time', sol.y[400][0], 0.8, 1e-8, 'Closed-form solution of I·dΩ/dt = −Q0·(Ω/Ω0)²'),
      N.check('Energy balance closes at the minimum-sink point', r.outputs.energy_residual + 1, 1, 1e-8, 'W·Vd = induced + profile + parasite power'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.autorotation_index / 0.006366 < 10) out.push(o.helicopter_criteria
      ? { severity: 'warn', title: 'Low stored rotor energy', detail: `Autorotation index ${(o.autorotation_index / 0.006366).toFixed(1)} ft³/lb; the rotor decays to the minimum usable speed in ${o.rotor_decay_time_s.toFixed(1)} s if collective is not lowered.`, action: 'Add blade tip mass or rotor speed, or lower disk loading.', basis: 'Autorotation index (helicopter guide values of about 20 ft³/lb single-engine, 10 ft³/lb multi-engine) and rotor-speed decay' }
      : { severity: 'advise', title: 'This configuration cannot rely on autorotation', detail: `Stored rotor energy is small (index ${(o.autorotation_index / 0.006366).toFixed(1)} ft³/lb, ${o.rotor_decay_time_s.toFixed(1)} s to the minimum rotor speed): normal for small or distributed lift rotors, which are not designed to autorotate.`, action: 'Show safe recovery by motor and rotor redundancy with one-rotor-out controllability (Suites 16 and 22), or a ballistic parachute.', basis: 'Stored rotor energy against hover power; helicopter autorotation guide values do not apply to this class' });
    out.push({ severity: 'info', title: 'Power-off glide performance', detail: `Minimum sink ${o.sink_min_ms.toFixed(1)} m/s at ${o.V_sink_min_ms.toFixed(0)} m/s; best glide ${o.glide_ratio.toFixed(1)}:1 at ${o.V_glide_ms.toFixed(0)} m/s; vertical autorotation about ${fmt(o.sink_vertical_ms)} m/s.`, action: 'Use these speeds for the height–velocity diagram and emergency procedures; confirm in flight test.', basis: 'Energy balance in steady autorotation' });
    out.push({ severity: 'advise', title: 'Stay clear of the vortex-ring band in powered descent', detail: `With vh = ${o.v_hover_induced_ms.toFixed(1)} m/s, avoid descent rates between about ${o.vrs_sink_lo_ms.toFixed(1)} and ${o.vrs_sink_hi_ms.toFixed(1)} m/s at forward speeds below about ${o.v_hover_induced_ms.toFixed(0)} m/s.`, action: 'Fly approaches with forward speed above vh or descent rate below the band; high disk loading (eVTOL, multirotor) moves the band to higher, more easily reached descent rates in absolute terms only if vh is large — check the approach profile in Suite 24.', basis: 'Textbook vortex-ring guideline 0.5–1.5 vh' });
    return out;
  },
};

export default {
  id: 'rotorcraft', n: 6,
  tagline: 'How much power the rotor needs, how the blades flap and stall, how the wake behaves, and what happens when the engine stops.',
  analyses: [hover, forward, bladedyn, inflow, wake, dynstall, autorot],
  consumes: [],
  provides: [
    { key: 'hover_power_W', label: 'Hover rotor power', unit: 'W' }, { key: 'FM', label: 'Figure of merit', unit: '-' }, { key: 'CT', label: 'Thrust coefficient', unit: '-' },
    { key: 'v_induced_ms', label: 'Hover induced velocity', unit: 'm/s' }, { key: 'collective_deg', label: 'Hover collective', unit: 'deg' }, { key: 'P_cruise_W', label: 'Rotor power at speed', unit: 'W' },
    { key: 'coning_deg', label: 'Coning angle', unit: 'deg' }, { key: 'a1_deg', label: 'Longitudinal flapping', unit: 'deg' }, { key: 'b1_deg', label: 'Lateral flapping', unit: 'deg' },
    { key: 'hub_vib_N', label: 'N/rev hub force', unit: 'N' }, { key: 'autorotation_index', label: 'Autorotation index', unit: 'm³/N' },
  ],
  handoff: [
    { model: 'Rotor CFD–CSD (coupled Navier–Stokes and structural dynamics)', why: 'Needs moving overset meshes, turbulence modelling and elastic blade coupling far beyond in-browser cost', tool: 'Overset RANS/DES solver coupled to a comprehensive code' },
    { model: 'Free-vortex wake and comprehensive rotorcraft analysis', why: 'Wake distortion, roll-up and full aeroelastic trim are replaced here by a rigid prescribed wake and rigid blades', tool: 'Comprehensive rotorcraft aeromechanics code (free-wake, multibody)' },
    { model: 'Elastic blade (rotating beam) flap–lag–torsion modes', why: 'Only rigid hinged-blade motion is integrated; the vibratory hub loads therefore omit elastic amplification', tool: 'Finite-element rotating-beam / multibody rotor code; Suite 10 for rotating-beam frequencies' },
    { model: 'Peters–He finite-state inflow (higher harmonics and radial shape functions)', why: 'Only the three Pitt–Peters states are solved', tool: 'Comprehensive rotorcraft code or flight-dynamics simulation' },
    { model: 'Full Leishman–Beddoes and ONERA dynamic-stall models with moment and drag', why: 'The demonstrator computes normal force only, with simplified switching logic and incompressible impulsive loads', tool: 'Comprehensive rotorcraft code with section-specific calibrated coefficients' },
    { model: 'Rotor–fuselage, rotor–tail-rotor and rotor–rotor forward-flight interference', why: 'Requires wake–body interaction modelling; only hover download and overlap interference factors are applied', tool: 'Panel/free-wake or CFD interaction analysis' },
    { model: 'Vortex-ring-state boundary prediction', why: 'Shown as a guideline band; the unsteady recirculating flow needs dedicated inflow models validated in flight test', tool: 'Flight test, free-wake or CFD' },
  ],
};
