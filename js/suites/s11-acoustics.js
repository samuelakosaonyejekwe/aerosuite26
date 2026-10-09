// Suite 11 — Aeroacoustics and Noise Prediction.
// Rotor and propeller harmonic noise (Gutin / Garrick–Watkins loading and compact thickness sources, checked against a
// retarded-time Ffowcs Williams–Hawkings point-source integration), source breakdown with jet (Lighthill) and empirical
// broadband models in one-third-octave bands, atmospheric absorption (ISO 9613-1) and A-weighting, a moving-source
// flyover with Doppler shift and ground footprint, cabin noise by mass-law transmission, and a finite-difference
// Helmholtz solver for duct silencers. Computational aeroacoustics on resolved flow fields is listed in `handoff`.

import * as N from '../core/numerics.js';
import { isa, G0 } from '../core/atmosphere.js';
import { METALS } from '../data/materials.js';
import { CERT_NOISE, CERT_NOISE_ISSUE } from '../data/ref/cert-noise.js';

const PI = Math.PI, PREF = 2e-5, dB = (p) => 20 * Math.log10(Math.max(p, 1e-12) / PREF), sumdB = (a) => 10 * Math.log10(a.reduce((s, v) => s + 10 ** (v / 10), 0) || 1e-30);
/** Bessel function of the first kind, integer order, from its integral representation (Simpson; exact to round-off). */
export const besselJ = (n, x) => N.simpson((t) => Math.cos(n * t - x * Math.sin(t)), 0, PI, 2 * Math.ceil(60 + Math.abs(x) + n)) / PI;
/** A-weighting [dB] (IEC 61672). */
export const aWeight = (f) => { const f2 = f * f; return 20 * Math.log10((12194 ** 2 * f2 * f2) / ((f2 + 20.6 ** 2) * Math.sqrt((f2 + 107.7 ** 2) * (f2 + 737.9 ** 2)) * (f2 + 12194 ** 2))) + 2.0; };
/** Pure-tone atmospheric absorption [dB/m], ISO 9613-1. T [K], rh 0–1, p [Pa]. */
export function absorption(f, T = 293.15, rh = 0.7, p = 101325) {
  const pr = 101325, T0 = 293.15, pa = p / pr, tr = T / T0, h = (rh * 100 * 10 ** (-6.8346 * (273.16 / T) ** 1.261 + 4.6151)) / pa;
  const frO = pa * (24 + (4.04e4 * h * (0.02 + h)) / (0.391 + h)), frN = pa * tr ** -0.5 * (9 + 280 * h * Math.exp(-4.17 * (tr ** (-1 / 3) - 1)));
  return 8.686 * f * f * ((1.84e-11 * Math.sqrt(tr)) / pa + tr ** -2.5 * ((0.01275 * Math.exp(-2239.1 / T)) / (frO + (f * f) / frO) + (0.1068 * Math.exp(-3352 / T)) / (frN + (f * f) / frN)));
}
const BANDS = N.range(28, (k) => 1000 * 10 ** ((k - 17) / 10)); // one-third-octave centres, 20 Hz – 10 kHz
const bandOf = (f) => N.clamp(Math.round(10 * Math.log10(f / 1000)) + 17, 0, 27);
/** Speech interference level: arithmetic mean of the 500 Hz, 1, 2 and 4 kHz octave-band levels, each octave summed from its three one-third-octave bands. */
export const speechInterference = (L) => N.mean([14, 17, 20, 23].map((k) => sumdB([L[k - 1], L[k], L[k + 1]])));
/** Generic broadband hump: share of the overall mean-square pressure in each band for a peak frequency fp. */
const hump = (fp) => { const w = BANDS.map((f) => { const x = f / fp; return (x * x) / (1 + x * x) ** 2; }), s = N.sum(w); return w.map((v) => v / s); };

// ---- rotor / propeller harmonic noise -------------------------------------------------------------
/** Radial stations with thrust, torque and blade-volume shares (loading ∝ x²√(1−x), a typical propeller shape). */
function stations(o) {
  const n = Math.max(1, Math.round(o.nRad)), xs = n === 1 ? [0.8] : N.range(n, (j) => 0.2 + (0.8 * (j + 0.5)) / n), w = xs.map((x) => (n === 1 ? 1 : x * x * Math.sqrt(1 - x))), sw = N.sum(w), swx = N.sum(w.map((v, j) => v * xs[j]));
  return xs.map((x, j) => ({ r: x * o.R, dT: (o.T * w[j]) / sw, dQ: (o.Q * w[j] * x) / swx, dVol: (0.685 * o.tc * o.chord ** 2 * 0.8 * o.R) / n }));
}
/**
 * RMS pressure of harmonic m at distance r and angle th (from the forward thrust axis), axial flight Mach M.
 * Loading: Gutin's formula in the Garrick–Watkins forward-flight form, summed over radial stations.
 * Thickness: compact rotating blade volume (monopole) in the same far-field form. The two are in phase quadrature.
 */
export function harmonic(o, st, m, r, th) {
  const b2 = 1 - o.M * o.M, x = r * Math.cos(th), y = r * Math.sin(th), S0 = Math.sqrt(x * x + b2 * y * y), n = m * o.B, k = (n * o.Om) / o.c; let L = 0, Tk = 0;
  for (const s of st) { const J = besselJ(n, (k * y * s.r) / S0); L += (k / (2 * Math.SQRT2 * PI * S0)) * ((s.dT * (o.M + x / S0)) / b2 - (s.dQ * o.c) / (o.Om * s.r * s.r)) * J; Tk += ((Math.SQRT2 * o.B * o.rho * s.dVol * (n * o.Om) ** 2) / (4 * PI * S0)) * J; }
  return { L: Math.abs(L), T: Math.abs(Tk), p: Math.hypot(L, Tk) };
}
/**
 * Retarded-time integration of the Ffowcs Williams–Hawkings equation for compact rotating sources (static hub and
 * observer): Lowson's far-field formula for the blade forces and ρ0 ∂²/∂t²[Vol/(4πr|1−Mr|)] for the blade volume.
 * Returns complex Fourier coefficients (already × B) of harmonics n = mB for loading and thickness.
 */
export function fwhRotating(o, st, obs, nH, nTau = 720) {
  const L = N.range(nH, () => [0, 0]), Tk = N.range(nH, () => [0, 0]);
  for (const s of st) for (let q = 0; q < nTau; q++) {
    const tau = (2 * PI * q) / (nTau * o.Om), ps = o.Om * tau, cp = Math.cos(ps), sp = Math.sin(ps), rx = obs[0], ry = obs[1] - s.r * cp, rz = -s.r * sp, r = Math.hypot(rx, ry, rz), hx = rx / r, hy = ry / r, hz = rz / r;
    const Mt = (o.Om * s.r) / o.c, Mr = Mt * (-hy * sp + hz * cp), hRad = hy * cp + hz * sp, Fb = -s.dT / o.B, Ft = s.dQ / (o.B * s.r), Fr = hx * Fb + Ft * (-hy * sp + hz * cp);
    const pL = ((-Ft * o.Om * hRad) / (1 - Mr) ** 2 + (Fr * -(o.Om * Mt) * hRad) / (1 - Mr) ** 3) / (4 * PI * o.c * r), t = tau + r / o.c;
    for (let m = 1; m <= nH; m++) { const n = m * o.B, ph = -n * o.Om * t, c = Math.cos(ph), sn = Math.sin(ph), wL = (o.B * pL * (1 - Mr)) / nTau, wT = (-o.B * (n * o.Om) ** 2 * o.rho * s.dVol) / (4 * PI * r * nTau); L[m - 1][0] += wL * c; L[m - 1][1] += wL * sn; Tk[m - 1][0] += wT * c; Tk[m - 1][1] += wT * sn; }
  }
  return { L, Tk };
}
const ROTOR = [
  { key: 'T_N', label: 'Thrust per rotor', unit: 'N', default: 3000, min: 0, group: 'Rotor / propeller' },
  { key: 'P_W', label: 'Shaft power per rotor', unit: 'W', default: 120000, min: 0, group: 'Rotor / propeller' },
  { key: 'R', label: 'Tip radius', unit: 'm', default: 0.95, min: 0.02, group: 'Rotor / propeller' },
  { key: 'B', label: 'Blades', unit: '', default: 2, min: 2, max: 12, step: 1, discrete: true, group: 'Rotor / propeller' },
  { key: 'rpm', label: 'Rotational speed', unit: 'rpm', default: 2700, min: 30, group: 'Rotor / propeller' },
  { key: 'chord', label: 'Mean blade chord', unit: 'm', default: 0.15, min: 0.002, group: 'Rotor / propeller' },
  { key: 'tc', label: 'Blade thickness ratio', unit: '-', default: 0.1, min: 0.02, max: 0.25, group: 'Rotor / propeller' },
  { key: 'n_rotors', label: 'Rotors or propellers', unit: '', default: 1, min: 0, max: 16, step: 1, discrete: true, group: 'Rotor / propeller', help: 'Summed as uncorrelated sources (+10·log N)' },
  { key: 'V_axial', label: 'Axial flight speed', unit: 'm/s', default: 0, min: 0, max: 250, group: 'Rotor / propeller', help: '0 for static thrust or hover. Edgewise (helicopter forward) flight is not represented' },
];
const SITE = [
  { key: 'alt_m', label: 'Altitude of the source', unit: 'm', default: 0, min: -500, max: 15000, group: 'Atmosphere' },
  { key: 'T_C', label: 'Air temperature on the path', unit: '°C', default: 15, min: -40, max: 50, group: 'Atmosphere' },
  { key: 'rh', label: 'Relative humidity', unit: '-', default: 0.7, min: 0.05, max: 1, group: 'Atmosphere' },
];
const hasRotor = (c) => c.rotor.R_m > 0 || (c.prop.prop_dia_m > 0 && c.prop.n_blades > 0);
const rotorDefaults = (c, up, d) => {
  if (!hasRotor(c)) return { n_rotors: 0 };
  const a = isa(c.atm.alt_m, c.atm.dISA_K);
  if (c.rotor.R_m > 0) { const nR = c.meta.type === 'helicopter' ? 1 : Math.max(1, c.prop.n_eng), T = d.W / nR; return { T_N: T, P_W: (c.meta.type === 'helicopter' ? up.rotorcraft?.hover_power_W : undefined) ?? T ** 1.5 / Math.sqrt(2 * a.rho * d.A_disk) / 0.7, R: c.rotor.R_m, B: c.rotor.n_blades, rpm: c.rotor.rpm, chord: c.rotor.chord_m, tc: 0.12, n_rotors: nR, V_axial: 0 }; }
  const R = c.prop.prop_dia_m / 2; return { T_N: c.prop.T0_N || undefined, P_W: c.prop.P0_W || undefined, R, B: c.prop.n_blades, rpm: c.prop.rpm, chord: 0.16 * R, tc: 0.08, n_rotors: c.prop.n_eng, V_axial: 0 };
};
const siteDefaults = (c) => ({ alt_m: c.site.elev_m ?? 0, T_C: c.site.T_C, rh: c.site.rh });
const rotorState = (i) => { const a = isa(i.alt_m, i.T_C + 273.15 - isa(i.alt_m).T), Om = (i.rpm * 2 * PI) / 60; return { a, T: i.T_N, Q: i.P_W / Math.max(Om, 1e-6), R: i.R, B: Math.round(i.B), Om, c: a.a, rho: a.rho, chord: i.chord, tc: i.tc, M: Math.min(i.V_axial / a.a, 0.9), nRad: i.nRad ?? 8 }; };

const tonal = {
  id: 'tonal', title: 'Rotor and propeller harmonic noise', fidelity: 'reduced-order',
  summary: 'Sound pressure level of the blade-passing tone and its harmonics from the steady blade loading and blade thickness, with the directivity pattern and the pressure waveform at the observer.',
  equations: ['Ffowcs Williams–Hawkings equation', 'Acoustic wave equation', 'Rotor loading noise models', 'Rotor thickness noise models', 'FW-H acoustic models', 'Rotor aerodynamic loading–acoustic radiation coupling'],
  applicable: (c) => (hasRotor(c) ? true : 'No open rotor or propeller is defined. Ducted-fan tones need a duct-mode model (see hand-off); jet and airframe noise are in the source-breakdown analysis.'),
  inputs: [...ROTOR,
    { key: 'r_obs', label: 'Observer distance from the hub', unit: 'm', default: 150, min: 1, group: 'Observer' },
    { key: 'theta_deg', label: 'Observer angle from the thrust axis', unit: 'deg', default: 105, min: 0, max: 180, group: 'Observer', help: '0 = ahead on the axis, 90 = in the disc plane, > 90 = behind (for a hovering rotor: below the disc)' },
    ...SITE.filter((f) => f.key !== 'rh'), // absorption is not applied to the tone levels here, so humidity is not an input
    { key: 'nHarm', label: 'Harmonics', unit: '', default: 8, min: 1, max: 20, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nRad', label: 'Radial source stations', unit: '', default: 8, min: 1, max: 40, step: 1, discrete: true, group: 'Numerics', help: '1 reproduces the classical effective-radius (0.8R) Gutin model' }],
  defaults: (c, up, d) => ({ ...rotorDefaults(c, up, d), alt_m: c.site.elev_m ?? 0, T_C: c.site.T_C, r_obs: c.rotor.R_m > 0 ? 150 : 300 }),
  run(i) {
    const o = rotorState(i), st = stations(o), nH = Math.round(N.clamp(i.nHarm, 1, 20)), th = N.rad(i.theta_deg), nR = Math.max(1, Math.round(i.n_rotors)), add = 10 * Math.log10(nR), bpf = (o.B * o.Om) / (2 * PI);
    const hs = N.range(nH, (k) => harmonic(o, st, k + 1, i.r_obs, th)), lev = hs.map((h) => dB(h.p) + add), levL = hs.map((h) => dB(h.L) + add), levT = hs.map((h) => dB(h.T) + add), fr = hs.map((_, k) => (k + 1) * bpf);
    const oaspl = sumdB(lev), dBA = sumdB(lev.map((v, k) => v + aWeight(fr[k]))), Mt = (o.Om * o.R) / o.c, Mh = Math.hypot(Mt, o.M), ths = N.linspace(2, 178, 89), dir = ths.map((t) => Math.max(dB(harmonic(o, st, 1, i.r_obs, N.rad(t)).p) + add, 0)), kp = N.argmax(dir);
    const one = stations({ ...o, nRad: 1 }), gut = dB(harmonic(o, one, 1, i.r_obs, th).L) + add, warnings = [];
    if (Mh > 0.9) warnings.push(`Helical tip Mach number ${Mh.toFixed(2)}: the sources are no longer acoustically compact, shocks form at the tip and quadrupole noise appears. Levels are under-predicted.`);
    if (i.r_obs < 4 * i.R) warnings.push('The observer is within four radii of the hub: the far-field formulas used here omit near-field terms.');
    if (o.M > 0.6) warnings.push('High axial Mach number: the Garrick–Watkins form assumes compact chordwise loading.');
    // waveform over one blade passage from the retarded-time integration (static hub)
    const w = fwhRotating({ ...o, M: 0 }, st, [i.r_obs * Math.cos(th), i.r_obs * Math.sin(th), 0], nH), tt = N.linspace(0, 1 / bpf, 120);
    const wave = (C) => tt.map((t) => 2 * N.sum(C.map((c, k) => c[0] * Math.cos(2 * PI * fr[k] * t) - c[1] * Math.sin(2 * PI * fr[k] * t))));
    return {
      kpis: [
        { key: 'bpf_Hz', label: 'Blade-passing frequency', value: bpf, unit: 'Hz' },
        { key: 'SPL_bpf_dB', label: 'Level of the blade-passing tone', value: Math.max(lev[0], 0), unit: 'dB' },
        { key: 'OASPL_tonal_dB', label: 'Overall tonal level', value: Math.max(oaspl, 0), unit: 'dB' },
        { key: 'tonal_dBA', label: 'A-weighted tonal level', value: Math.max(dBA, 0), unit: 'dBA' },
        { key: 'SPL_loading_dB', label: 'Loading noise at BPF', value: Math.max(levL[0], 0), unit: 'dB' },
        { key: 'SPL_thickness_dB', label: 'Thickness noise at BPF', value: Math.max(levT[0], 0), unit: 'dB' },
        { key: 'SPL_gutin_dB', label: 'Gutin effective-radius estimate at BPF', value: Math.max(gut, 0), unit: 'dB', note: 'Single source ring at 0.8R, loading only' },
        { key: 'tip_mach', label: 'Rotational tip Mach number', value: Mt, unit: '-', note: 'Tone level rises steeply with tip Mach: about 0.6–0.7 for helicopter rotors, 0.7–0.8 for propellers at take-off, below 0.6 for quiet designs' },
        { key: 'tip_mach_helical', label: 'Helical tip Mach number', value: Mh, unit: '-', status: Mh < 0.85 ? 'ok' : Mh < 0.95 ? 'warn' : 'bad', note: 'Above about 0.85–0.9 the tip flow is transonic and impulsive noise appears' },
        { key: 'theta_peak_deg', label: 'Direction of the loudest BPF radiation', value: ths[kp], unit: 'deg' },
        { key: 'SPL_peak_dir_dB', label: 'BPF level in the loudest direction', value: dir[kp], unit: 'dB' },
      ],
      plots: [
        { type: 'bar', title: 'Harmonic levels at the observer', ylabel: 'Sound pressure level [dB re 20 µPa]', categories: fr.map((f) => `${f.toFixed(0)} Hz`), series: [{ name: 'Loading', y: levL.map((v) => Math.max(v, 0)) }, { name: 'Thickness', y: levT.map((v) => Math.max(v, 0)) }, { name: 'Total', y: lev.map((v) => Math.max(v, 0)) }] },
        { type: 'polar', title: 'Directivity of the blade-passing tone', rlabel: 'SPL [dB]', series: [{ name: `BPF at ${i.r_obs} m`, theta_deg: ths, r: dir }] },
        { type: 'line', title: 'Acoustic pressure over one blade passage (static, one rotor)', xlabel: 'Time [ms]', ylabel: 'Pressure [Pa]', series: [{ name: 'Loading', x: tt.map((t) => t * 1e3), y: wave(w.L) }, { name: 'Thickness', x: tt.map((t) => t * 1e3), y: wave(w.Tk) }, { name: 'Total', x: tt.map((t) => t * 1e3), y: wave(w.L.map((c, k) => [c[0] + w.Tk[k][0], c[1] + w.Tk[k][1]])) }] },
      ],
      outputs: { harmonic_freq_Hz: fr, harmonic_SPL_dB: lev.map((v) => Math.max(v, 0)) },
      warnings,
      models: ['Gutin steady-loading noise in the Garrick–Watkins axial-flight form, integrated over radial stations', 'Compact-chord thickness (monopole) noise of the rotating blade volume', 'Retarded-time Ffowcs Williams–Hawkings point-source integration (Lowson) for the waveform', 'Uncorrelated summation over rotors'],
      assumptions: ['Steady, azimuthally uniform blade loads: no blade–vortex interaction, inflow distortion or unsteady loading harmonics', 'Acoustically compact chord; far field; free field (no ground reflection or installation effects)', 'Axial flight only; assumed radial load shape x²√(1−x)', 'Atmospheric absorption is applied in the source-breakdown and flyover analyses, not here'],
    };
  },
  convergence: { param: 'nRad', label: 'Radial source stations', levels: [2, 4, 8, 16, 32], metric: 'SPL_bpf_dB' },
  calibration: { params: [{ key: 'T_N', min: 100, max: 1e5 }, { key: 'tc', min: 0.03, max: 0.2 }], sweep: 'rpm', target: 'SPL_bpf_dB', note: 'Supply measured blade-passing-tone level against rotational speed at a known microphone position; effective thrust (loading) and blade thickness are the adjustable source strengths.' },
  verify() {
    const o = { T: 2000, Q: 300, R: 1, B: 2, Om: 200, c: 340, rho: 1.225, chord: 0.15, tc: 0.1, M: 0, nRad: 1 }, st = stations(o), r = 200, out = [];
    for (const td of [60, 120]) { const th = N.rad(td), w = fwhRotating(o, st, [r * Math.cos(th), r * Math.sin(th), 0], 2, 1440), h = harmonic(o, st, 1, r, th), h2 = harmonic(o, st, 2, r, th);
      out.push(N.check(`Gutin loading tone equals FW-H integration, θ = ${td}°`, h.L, Math.SQRT2 * Math.hypot(...w.L[0]), 0.01, 'Retarded-time Lowson point-force solution; far-field terms O(R/r)'));
      out.push(N.check(`Second harmonic loading, θ = ${td}°`, h2.L, Math.SQRT2 * Math.hypot(...w.L[1]), 0.015, 'Retarded-time Lowson point-force solution'));
      out.push(N.check(`Thickness tone equals FW-H monopole integration, θ = ${td}°`, h.T, Math.SQRT2 * Math.hypot(...w.Tk[0]), 0.01, 'ρ0 ∂²/∂t² [Vol / 4πr(1−Mr)] evaluated at retarded time')); }
    out.push(N.check('Bessel J1(1)', besselJ(1, 1), 0.4400505857, 1e-9, 'Abramowitz & Stegun'), N.check('Bessel J4(2.5)', besselJ(4, 2.5), 0.0737819, 1e-5, 'Abramowitz & Stegun'));
    out.push(N.check('Doubling the distance lowers the tone by 6.02 dB', dB(harmonic(o, st, 1, 200, 1).p) - dB(harmonic(o, st, 1, 400, 1).p), 20 * Math.log10(2), 1e-9, 'Spherical spreading'));
    return out;
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.tip_mach_helical > 0.85) out.push({ severity: o.tip_mach_helical > 0.95 ? 'critical' : 'warn', title: 'Transonic blade tips', detail: `Helical tip Mach ${o.tip_mach_helical.toFixed(2)}: shocks and quadrupole noise appear, which this compact-source model does not include; BPF tone ${o.SPL_bpf_dB.toFixed(0)} dB at ${i.r_obs} m is a lower bound.`, action: 'Reduce rotational speed or diameter, or use thin swept tips; recover thrust with more blades or chord.', basis: 'Helical tip Mach number above the range of acoustically compact sources (about 0.85–0.9)' });
    else if (o.tip_mach > 0.65) out.push({ severity: 'advise', title: 'Tip speed is the main lever on tone noise', detail: `Rotational tip Mach ${o.tip_mach.toFixed(2)} (normal for a propeller at take-off power); BPF tone ${o.SPL_bpf_dB.toFixed(0)} dB at ${i.r_obs} m.`, action: 'Reduce rotational speed and recover thrust with more blades, larger diameter or more chord; a 10% tip-speed cut is typically worth several decibels and usually improves static efficiency as well.', basis: 'Bessel-function radiation efficiency J_mB(mB·Mt·sin θ)' });
    if (o.SPL_thickness_dB > o.SPL_loading_dB) out.push({ severity: 'advise', title: 'Thickness noise exceeds loading noise', detail: `${o.SPL_thickness_dB.toFixed(0)} dB against ${o.SPL_loading_dB.toFixed(0)} dB at BPF.`, action: 'Thin the outer blade, sweep the tip, or lower the tip Mach number.', basis: 'Monopole source strength ∝ blade volume × (mBΩ)²' });
    out.push({ severity: 'info', title: 'More blades move energy to higher, weaker harmonics', detail: `BPF is ${o.bpf_Hz.toFixed(0)} Hz with ${Math.round(i.B)} blades.`, action: 'Re-run with one more blade at the same thrust and tip speed: the Bessel order rises and the tone falls, but the frequency moves towards the ear’s most sensitive range, so compare A-weighted levels.', basis: 'Harmonic order mB' });
    return out;
  },
};

// ---- source breakdown in one-third-octave bands -----------------------------------------------------
/** Band levels of each source at distance r, angle th from the flight/thrust axis. i = inputs of the `spectrum` analysis. */
function sources(i, r, th) {
  const o = rotorState(i), a = o.a, out = [], nR = Math.round(i.n_rotors);
  if (nR > 0 && i.T_N > 0 && i.R > 0) {
    const st = stations(o), tone = new Array(28).fill(-99), bpf = (o.B * o.Om) / (2 * PI);
    for (let m = 1; m <= 10 && m * bpf < 11000; m++) { const b = bandOf(m * bpf); tone[b] = sumdB([tone[b], dB(harmonic(o, st, m, r, th).p) + 10 * Math.log10(nR)]); }
    out.push({ name: 'Rotor / propeller tones', L: tone });
    // Schlegel–King–Mull vortex (broadband) noise, empirical: SPL = 10·log(6.1e-27·A_b·V_0.7⁶/1e-16) + 20·log(CL/0.4) at 300 ft, A_b in ft², V in ft/s
    const Ab = o.B * i.chord * i.R * 0.8, V07 = 0.7 * o.Om * i.R, CT = i.T_N / (o.rho * PI * i.R ** 2 * (o.Om * i.R) ** 2), sig = (o.B * i.chord) / (PI * i.R), clb = Math.max((6 * CT) / sig, 0.05);
    const L300 = 10 * Math.log10((6.1e-27 * (Ab / 0.092903) * (V07 / 0.3048) ** 6) / 1e-16) + 20 * Math.log10(clb / 0.4) + i.K_bb + 10 * Math.log10(nR), fp = (0.28 * V07) / (i.tc * i.chord + (i.chord * clb) / 5.7);
    out.push({ name: 'Rotor broadband (empirical)', L: hump(fp).map((w) => L300 + 20 * Math.log10(91.44 / r) + 10 * Math.log10(w)) });
  }
  if (i.n_jets > 0 && i.U_jet > 0) {
    const A = (PI * i.D_jet ** 2) / 4, Ue = Math.max(i.U_jet - i.V_flight, 1), W = i.K_jet * (i.U_jet / a.a) ** 5 * 0.5 * a.rho * A * i.U_jet ** 3 * (Ue / i.U_jet) ** i.n_rel, Mc = (0.62 * i.U_jet) / a.a;
    const D = (t) => Math.max(1 - Mc * Math.cos(t), 0.4) ** -5, Dm = 0.5 * N.simpson((t) => D(t) * Math.sin(t), 0, PI, 200), L = 10 * Math.log10(W / 1e-12) - 10 * Math.log10(4 * PI * r * r) + 10 * Math.log10(D(PI - th) / Dm) + 10 * Math.log10(i.n_jets);
    out.push({ name: 'Jet mixing (Lighthill)', L: hump((0.25 * i.U_jet) / i.D_jet).map((w) => L + 10 * Math.log10(w)), W });
  }
  if (i.S_wing > 0 && i.V_flight > 10) {
    // Fink-type clean-airframe trailing-edge noise (empirical; formula in knots and feet, overhead)
    const cbar = i.S_wing / i.b_wing, dl = 0.37 * cbar * ((i.V_flight * cbar) / a.nu) ** -0.2, L = 50 * Math.log10(i.V_flight / 51.444) + 10 * Math.log10((dl * i.b_wing) / (r * r)) + 101.3 + i.K_af;
    out.push({ name: 'Airframe trailing edge (empirical)', L: hump((0.1 * i.V_flight) / dl).map((w) => L + 10 * Math.log10(w)) });
  }
  return { out, a };
}
/** Total and per-source levels at the receiver with absorption, ground effect and A-weighting. */
function received(i, r, th) {
  const { out } = sources(i, r, th), T = i.T_C + 273.15, p = isa(i.alt_m).p, att = BANDS.map((f) => (i.absorb ? absorption(f, T, i.rh, p) * r : 0) - (i.ground ? 3 : 0));
  const src = out.map((s) => { const L = s.L.map((v, k) => v - att[k]); return { name: s.name, L, oa: sumdB(L), dBA: sumdB(L.map((v, k) => v + aWeight(BANDS[k]))), W: s.W }; });
  const tot = BANDS.map((_, k) => sumdB(src.map((s) => s.L[k])));
  return { src, tot, oa: src.length ? sumdB(tot) : 0, dBA: src.length ? sumdB(tot.map((v, k) => v + aWeight(BANDS[k]))) : 0 };
}
const SRC_INPUTS = [...ROTOR,
  { key: 'n_jets', label: 'Jet engines', unit: '', default: 0, min: 0, max: 8, step: 1, discrete: true, group: 'Jet' },
  { key: 'U_jet', label: 'Jet exhaust velocity', unit: 'm/s', default: 300, min: 0, max: 1200, group: 'Jet', help: 'Fully mixed equivalent velocity. Class estimate from thrust and bypass ratio; replace with cycle data from Suite 7' },
  { key: 'D_jet', label: 'Nozzle equivalent diameter', unit: 'm', default: 1, min: 0.02, group: 'Jet' },
  { key: 'K_jet', label: 'Lighthill acoustic-efficiency constant', unit: '-', default: 1e-4, min: 1e-5, max: 1e-3, group: 'Jet', help: 'Acoustic power = K·M_j⁵ × jet kinetic power; about 10⁻⁴ for subsonic jets' },
  { key: 'n_rel', label: 'Flight-effect exponent', unit: '-', default: 5, min: 0, max: 8, group: 'Jet', help: 'Jet noise scales with ((U_j − V)/U_j)^n in flight (empirical)' },
  { key: 'S_wing', label: 'Wing area', unit: 'm²', default: 0, min: 0, group: 'Airframe' },
  { key: 'b_wing', label: 'Wing span', unit: 'm', default: 10, min: 0.1, group: 'Airframe' },
  { key: 'V_flight', label: 'Flight speed', unit: 'm/s', default: 60, min: 0, max: 300, group: 'Airframe' },
  { key: 'K_bb', label: 'Rotor broadband calibration offset', unit: 'dB', default: 0, min: -15, max: 15, group: 'Calibration' },
  { key: 'K_af', label: 'Airframe noise calibration offset', unit: 'dB', default: 0, min: -15, max: 15, group: 'Calibration', help: '0 = aerodynamically clean; about +8 dB for conventional transports with gear and flaps retracted' },
  ...SITE,
  { key: 'absorb', label: 'Apply atmospheric absorption', type: 'bool', default: true, group: 'Atmosphere' },
  { key: 'ground', label: 'Add 3 dB for ground reflection', type: 'bool', default: true, group: 'Atmosphere', help: 'Energy doubling at a microphone near hard ground' },
];
const srcDefaults = (c, up, d) => {
  const p = c.prop, jet = p.type === 'turbofan' || p.type === 'turbojet', Uj = 1000 / (1 + (p.type === 'turbojet' ? 0 : p.bpr)) ** 0.6, wing = c.wing.S_m2 > 0;
  const Vto = wing ? up.performance?.V_2_ms ?? 1.2 * Math.sqrt((2 * d.W) / (1.225 * c.wing.S_m2 * c.aero.CLmax_to)) : 0; // take-off safety speed V2 = 1.2·VS in the take-off configuration
  return { ...rotorDefaults(c, up, d), ...siteDefaults(c), n_jets: jet ? p.n_eng : 0, U_jet: jet ? Uj : undefined, D_jet: jet ? Math.sqrt((4 * p.T0_N) / (PI * 1.225 * Uj * Uj)) : undefined, S_wing: c.wing.S_m2, b_wing: c.wing.b_m || undefined, V_flight: Vto };
};
const spectrum = {
  id: 'spectrum', title: 'Noise source breakdown and one-third-octave spectrum', fidelity: 'reduced-order',
  summary: 'Adds rotor or propeller tones, rotor broadband noise, jet mixing noise and airframe noise at an observer, in one-third-octave bands with atmospheric absorption and A-weighting, to show which source sets the level people hear.',
  equations: ['Lighthill acoustic analogy', 'Acoustic propagation equations', 'Acoustic intensity equations', 'Acoustic energy equations', 'Broadband noise models', 'Engine noise models', 'Stochastic broadband noise prediction formulations'],
  inputs: [...SRC_INPUTS,
    { key: 'r_obs', label: 'Observer distance', unit: 'm', default: 300, min: 1, group: 'Observer', help: '450 m is the lateral reference distance of the transport noise certification (Annex 16 / FAR 36); 150 m is the helicopter overflight height' },
    { key: 'theta_deg', label: 'Observer angle from the forward axis', unit: 'deg', default: 110, min: 0, max: 180, group: 'Observer', help: 'Jet noise peaks towards the rear (130–150°)' }],
  defaults: (c, up, d) => ({ ...srcDefaults(c, up, d), r_obs: c.mass.mtow_kg > 5700 ? 450 : 150 }),
  run(i) {
    const r = received(i, i.r_obs, N.rad(i.theta_deg)), warnings = [], top = r.src.slice().sort((p, q) => q.dBA - p.dBA)[0], jetS = r.src.find((s) => s.W);
    if (!r.src.length) warnings.push('No noise source is active: enter a rotor, a jet or a wing with a flight speed.');
    if (r.src.some((s) => s.name.includes('empirical'))) warnings.push('Broadband rotor and airframe levels come from empirical scaling laws with a generic spectral shape; expect ±5 dB or more until calibrated against measurements of a similar vehicle.');
    if (i.U_jet / isa(i.alt_m).a > 1.2 && i.n_jets > 0) warnings.push('Jet Mach number above 1.2: shock-associated noise and the U³ supersonic scaling are not modelled.');
    const kMax = N.argmax(r.tot.map((v, k) => v + aWeight(BANDS[k])));
    return {
      kpis: [
        { key: 'OASPL_dB', label: 'Overall sound pressure level', value: r.oa, unit: 'dB' },
        { key: 'SPL_dBA', label: 'A-weighted level', value: r.dBA, unit: 'dBA', status: r.dBA < 85 ? 'ok' : 'warn', note: 'For orientation: about 65 dBA is conversational level; 85 dBA is the usual hearing-protection threshold for sustained exposure. Community limits depend on the airport and the metric' },
        { key: 'f_peak_dBA_Hz', label: 'Band with the highest A-weighted level', value: BANDS[kMax], unit: 'Hz' },
        { key: 'dominant_dBA', label: 'A-weighted level of the dominant source', value: top ? top.dBA : 0, unit: 'dBA', note: top ? top.name : 'none' },
        { key: 'jet_PWL_dB', label: 'Jet sound power level (per engine)', value: jetS ? 10 * Math.log10(jetS.W / 1e-12) : 0, unit: 'dB re 1 pW' },
        { key: 'absorption_1k_dB_km', label: 'Atmospheric absorption at 1 kHz', value: absorption(1000, i.T_C + 273.15, i.rh, isa(i.alt_m).p) * 1000, unit: 'dB/km' },
        { key: 'n_sources', label: 'Active source models', value: r.src.length, unit: '' },
      ],
      plots: [
        { type: 'line', title: `One-third-octave spectrum at ${i.r_obs} m`, xlabel: 'Band centre frequency [Hz]', ylabel: 'Band level [dB re 20 µPa]', xlog: true, series: [...r.src.slice(0, 4).map((s) => ({ name: s.name, x: BANDS, y: s.L.map((v) => (v > 0 ? v : NaN)), style: 'line+points' })), { name: 'Total', x: BANDS, y: r.tot.map((v) => Math.max(v, 0)) }] },
        { type: 'bar', title: 'Source ranking', ylabel: 'Level [dB / dBA]', categories: r.src.length ? r.src.map((s) => s.name) : ['none'], series: [{ name: 'Unweighted', y: r.src.length ? r.src.map((s) => Math.max(s.oa, 0)) : [0] }, { name: 'A-weighted', y: r.src.length ? r.src.map((s) => Math.max(s.dBA, 0)) : [0] }] },
        { type: 'line', title: 'A-weighted total spectrum', xlabel: 'Band centre frequency [Hz]', ylabel: 'A-weighted band level [dBA]', xlog: true, series: [{ name: 'Total', x: BANDS, y: r.tot.map((v, k) => Math.max(v + aWeight(BANDS[k]), 0)) }] },
      ],
      tables: [{ title: 'Sources at the observer', columns: ['Source', 'OASPL [dB]', 'Level [dBA]'], rows: r.src.map((s) => [s.name, +s.oa.toFixed(1), +s.dBA.toFixed(1)]) }],
      warnings,
      models: ['Harmonic rotor noise (see the tonal analysis)', 'Schlegel–King–Mull rotor vortex-noise correlation (empirical)', 'Lighthill U⁸ jet mixing noise with convective-amplification directivity and a relative-velocity flight effect', 'Fink-type trailing-edge airframe noise scaling, V⁵ (empirical)', 'ISO 9613-1 atmospheric absorption; A-weighting; one-third-octave energy summation'],
      assumptions: ['Sources are uncorrelated and add on an energy basis', 'Generic single-hump broadband spectra; tones assigned to the band containing them', 'Homogeneous still atmosphere: no refraction by wind or temperature gradients, no shielding', 'Fan, core, combustion and turbine noise of turbofans, tail-rotor noise and blade–vortex interaction are not modelled, so turbofan and helicopter totals are lower bounds',
        'Unverified constants (illustrative until calibrated): jet acoustic efficiency 10⁻⁴·M⁵, flight-effect exponent 5, jet peak Strouhal number 0.25, and the clean-airframe constant 101.3 dB (one opened reference prints a form equivalent to 107.5 dB for jet aircraft, 6.2 dB higher). The class estimate of jet velocity 1000/(1 + BPR)^0.6 m/s is a rough default',
        'Sourced constants: rotor vortex-noise law and Strouhal number 0.28 (Schlegel, King & Mull as reproduced in JPL TR 32-1462), jet convection factor 0.62 and airframe peak Strouhal number 0.1 (NASA RP-1258), +8 dB from clean to conventional airframes (NASA TM-83199)'],
    };
  },
  calibration: { params: [{ key: 'K_jet', min: 3e-5, max: 4e-4 }, { key: 'K_bb', min: -12, max: 12 }, { key: 'K_af', min: -5, max: 15 }], sweep: 'r_obs', target: 'SPL_dBA', note: 'Supply measured A-weighted (or overall, target OASPL_dB) level against distance or against jet velocity / rotor speed; the three constants set the absolute level of each source.' },
  verify() {
    const b = Object.fromEntries(spectrum.inputs.map((f) => [f.key, f.default])), jet = { ...b, n_rotors: 0, n_jets: 1, U_jet: 250, D_jet: 0.8, S_wing: 0, V_flight: 0, absorb: false, ground: false, r_obs: 200, theta_deg: 90 };
    const k = (x) => N.kv(spectrum.run(x)), j1 = k(jet), j2 = k({ ...jet, r_obs: 400 }), j3 = k({ ...jet, U_jet: 500 });
    return [
      N.check('A-weighting at 1 kHz is 0 dB', 1 + aWeight(1000), 1, 1e-3, 'IEC 61672 definition'),
      N.check('A-weighting at 100 Hz', aWeight(100), -19.1, 3e-3, 'IEC 61672 table'),
      N.check('Atmospheric absorption at 1 kHz, 20 °C, 70% RH', absorption(1000) * 1000, 5.0, 0.02, 'ISO 9613-2 Table 2 (dB/km)'),
      N.check('Atmospheric absorption at 4 kHz, 20 °C, 70% RH', absorption(4000) * 1000, 22.9, 0.02, 'ISO 9613-2 Table 2 (dB/km)'),
      N.check('Atmospheric absorption at 500 Hz, 10 °C, 70% RH', absorption(500, 283.15) * 1000, 1.9, 0.03, 'ISO 9613-2 Table 2 (dB/km)'),
      N.check('Spherical spreading: 6.02 dB per doubling', j1.OASPL_dB - j2.OASPL_dB, 20 * Math.log10(2), 1e-9, 'Inverse-square law'),
      N.check('Lighthill eighth-power law: 24.08 dB per velocity doubling', j3.jet_PWL_dB - j1.jet_PWL_dB, 80 * Math.log10(2), 1e-9, 'Lighthill (1952)'),
      N.check('Band levels sum to the overall level', sumdB(hump(300).map((w) => 90 + 10 * Math.log10(w))), 90, 1e-12, 'Energy summation'),
    ];
  },
  validation: [{ name: 'Atmospheric absorption, 20 °C and 70% humidity', source: 'ISO 9613-2:1996 Table 2: 5.0 dB/km at 1 kHz', inputs: { T_C: 20, rh: 0.7, alt_m: 0 }, sweep: { key: 'r_obs', values: [300] }, target: 'absorption_1k_dB_km', observed: [5.0], tol_pct: 3 }],
  recommend(res, i) {
    const o = res.outputs, t = res.tables[0].rows.slice().sort((a, b) => b[2] - a[2]), out = [];
    if (t.length) out.push({ severity: 'info', title: `${t[0][0]} sets the perceived level`, detail: `${t[0][2]} dBA of ${o.SPL_dBA.toFixed(1)} dBA total at ${i.r_obs} m${t[1] ? `; next is ${t[1][0]} at ${t[1][2]} dBA` : ''}.`, action: 'Work on the top source first: reducing any other source by 10 dB changes the total by less than 0.5 dB once it is 10 dB below the leader.', basis: 'Energy summation of uncorrelated sources' });
    if (i.n_jets > 0 && t.length && t[0][0].startsWith('Jet')) out.push({ severity: 'advise', title: 'Jet velocity is the lever', detail: `Exhaust velocity ${i.U_jet.toFixed(0)} m/s.`, action: 'A higher bypass ratio that lowers jet velocity by 10% cuts jet noise by about 3.7 dB and improves propulsive efficiency, fuel burn and CO₂ at the same time.', basis: 'Lighthill U⁸ law' });
    if (o.SPL_dBA > 85) out.push({ severity: 'warn', title: 'High level at the observer', detail: `${o.SPL_dBA.toFixed(0)} dBA.`, action: 'Increase stand-off distance or height, adjust the operating procedure (reduced thrust or rotor speed), and assess community exposure with the flyover footprint.', basis: 'Occupational and community noise practice (limits vary by jurisdiction)' });
    return out;
  },
};

// ---- moving source: flyover time history, Doppler shift and ground footprint -----------------------------
const flyLevel = (i, r, cosTh, cosPhi, f, att) => i.L_ref + 20 * Math.log10(i.r_ref / r) - 40 * Math.log10(1 - (i.V / i.c0) * cosTh) + (i.directivity === 'Dipole (vertical axis)' ? 20 * Math.log10(Math.max(Math.abs(cosPhi), 0.1)) : 0) - att * Math.max(r - i.r_ref, 0) + aWeight(f);
function footprint(i, att, n) {
  const g = N.rad(i.gamma_deg), fA = aWeight(i.f0), rTh = i.r_ref * 10 ** ((i.L_ref + fA - i.L_th) / 20), hover = i.V < 0.5;
  if (!(rTh > i.h0)) return { area: 0, x: [0, 1], y: [0, 1], z: [[0, 0], [0, 0]], track: 0, rTh };
  const S = hover ? 0 : Math.min(Math.sin(g) > 1e-3 ? (rTh - i.h0) / Math.sin(g) : 30000, 30000), nP = hover ? 1 : 60, path = N.range(nP, (k) => { const s = nP > 1 ? (S * k) / (nP - 1) : 0; return [s * Math.cos(g), i.h0 + s * Math.sin(g)]; });
  const x0 = -rTh, x1 = S * Math.cos(g) + rTh, nx = n, ny = Math.max(4, Math.round(n / 2)), dx = (x1 - x0) / nx, dy = rTh / ny, xs = N.range(nx, (k) => x0 + (k + 0.5) * dx), ys = N.range(ny, (k) => (k + 0.5) * dy); let cells = 0;
  const z = ys.map((y) => xs.map((x) => { let best = -1e9; for (const [px, pz] of path) { const r = Math.hypot(x - px, y, pz), cT = ((x - px) * Math.cos(g) - pz * Math.sin(g)) / r, L = flyLevel(i, r, cT, pz / r, i.f0, att); if (L > best) best = L; } if (best >= i.L_th) cells++; return best; }));
  return { area: 2 * cells * dx * dy, x: xs, y: ys, z, track: S, rTh };
}
const CERT_CLASSES = ['Jet aeroplane', 'Propeller aeroplane', 'Helicopter', 'None'];
/**
 * Certified noise levels of the aircraft closest in maximum take-off mass (EASA type-certificate data sheets for noise).
 * Returns the data set used, the index of the level most comparable to a flyover and the n nearest rows.
 */
export function certNeighbours(cls, mtom, n = 5) {
  const key = cls === 'Jet aeroplane' ? 'jets' : cls === 'Propeller aeroplane' ? (mtom > 8618 ? 'heavyProp' : 'lightProp') : cls === 'Helicopter' ? (mtom <= 3175 ? 'rotor11' : 'rotor8') : null;
  if (!key || !(mtom > 0)) return null;
  const set = CERT_NOISE[key], rows = set.rows.slice().sort((a, b) => Math.abs(Math.log(a[2] / mtom)) - Math.abs(Math.log(b[2] / mtom))).slice(0, n), ref = set.points.length > 1 ? 1 : 0;
  const lev = rows.map((r) => r[set.first + 2 * ref]).sort((a, b) => a - b);
  return { key, set, rows, ref, median: lev[lev.length >> 1], lo: lev[0], hi: lev[lev.length - 1], direct: key === 'rotor11' || key === 'lightProp' };
}
const flyover = {
  id: 'flyover', title: 'Flyover: moving source, Doppler shift and ground footprint', fidelity: 'numerical',
  summary: 'Time history of the level heard on the ground as the aircraft passes, solved at the retarded (emission) time with convective amplification and Doppler shift, the single-event exposure level, and the ground area enclosed by a chosen noise contour during climb-out or hover.',
  equations: ['Ffowcs Williams–Hawkings equation', 'Doppler frequency relations', 'Convected wave equation', 'Acoustic propagation equations', 'Acoustic intensity equations', 'FW-H acoustic models'],
  inputs: [
    { key: 'L_ref', label: 'Source level at the reference distance', unit: 'dB', default: 85, min: 20, max: 160, group: 'Source', help: 'Level of the stationary source at its characteristic frequency. From the source-breakdown analysis it is set so that the A-weighted level matches the full spectrum' },
    { key: 'r_ref', label: 'Reference distance', unit: 'm', default: 150, min: 1, group: 'Source' },
    { key: 'f0', label: 'Characteristic source frequency', unit: 'Hz', default: 200, min: 10, max: 10000, group: 'Source', help: 'Blade-passing frequency or the peak band; used for Doppler shift, absorption and A-weighting' },
    { key: 'directivity', label: 'Source directivity', type: 'select', options: ['Monopole (uniform)', 'Dipole (vertical axis)'], default: 'Monopole (uniform)', group: 'Source', help: 'A lift or thrust dipole radiates mostly below the vehicle' },
    { key: 'V', label: 'Flight speed', unit: 'm/s', default: 70, min: 0, max: 300, group: 'Flight path', help: '0 for hover' },
    { key: 'h', label: 'Height over the microphone', unit: 'm', default: 300, min: 5, group: 'Flight path' },
    { key: 'y_side', label: 'Sideline offset of the microphone', unit: 'm', default: 0, min: 0, group: 'Flight path' },
    { key: 'gamma_deg', label: 'Climb angle for the footprint', unit: 'deg', default: 6, min: 0, max: 90, group: 'Flight path' },
    { key: 'h0', label: 'Height at the start of the footprint track', unit: 'm', default: 15, min: 1, group: 'Flight path', help: 'Hover height when the speed is zero' },
    { key: 'L_th', label: 'Footprint contour level', unit: 'dBA', default: 65, min: 30, max: 110, group: 'Footprint' },
    ...SITE.filter((f) => f.key !== 'alt_m'),
    { key: 'absorb', label: 'Apply atmospheric absorption', type: 'bool', default: true, group: 'Atmosphere' },
    { key: 'mtom_kg', label: 'Maximum take-off mass, to look up certified aircraft', unit: 'kg', default: 0, min: 0, group: 'Reality check', help: '0 = no comparison. Certified levels of the nearest aircraft by mass are listed beside the prediction' },
    { key: 'cert_class', label: 'Certification class for the comparison', type: 'select', options: CERT_CLASSES, default: 'None', group: 'Reality check' },
    { key: 'nGrid', label: 'Footprint grid points along track', unit: '', default: 80, min: 20, max: 300, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => { const wing = c.wing.S_m2 > 0, V = wing ? (1.25 / 1.2) * (up.performance?.V_2_ms ?? 1.2 * Math.sqrt((2 * d.W) / (1.225 * c.wing.S_m2 * c.aero.CLmax_to))) : 0, roc = up.performance?.roc_sl_ms; // initial climb at about V2 + 10 kt; climb rate at the airfield, not at the cruise point
    // source level and characteristic frequency: from the source-breakdown analysis when it has run, otherwise from the same source models on the case
    let dBA = up.acoustics?.SPL_dBA, fpk = up.acoustics?.f_peak_dBA_Hz ?? up.acoustics?.bpf_Hz; const rRef = c.mass.mtow_kg > 5700 ? 450 : 150;
    if (dBA === undefined || !fpk) {
      const si = Object.fromEntries(spectrum.inputs.map((f) => [f.key, f.default])), sd = spectrum.defaults(c, up, d);
      for (const k of Object.keys(sd)) if (sd[k] !== undefined && sd[k] !== null && !(typeof sd[k] === 'number' && !Number.isFinite(sd[k]))) si[k] = sd[k];
      const r = received(si, rRef, N.rad(si.theta_deg));
      if (r.src.length) { dBA = r.dBA; fpk = BANDS[N.argmax(r.tot.map((v, k) => v + aWeight(BANDS[k])))]; } else { dBA = undefined; fpk = undefined; }
    }
    return { L_ref: dBA !== undefined && fpk ? dBA - aWeight(fpk) : undefined, r_ref: dBA !== undefined && fpk ? rRef : undefined, f0: fpk, V: wing ? V : c.flight.V_ms, h: wing ? (c.mass.mtow_kg > 5700 ? 450 : 150) : 150, gamma_deg: wing ? (roc > 0 ? N.clamp(N.deg(Math.asin(Math.min(roc / Math.max(V, 1), 0.5))), 2, 30) : undefined) : 0, h0: wing ? undefined : 150, T_C: c.site.T_C, rh: c.site.rh, mtom_kg: c.mass.mtow_kg,
      cert_class: c.meta.type === 'helicopter' ? 'Helicopter' : c.meta.type !== 'aeroplane' ? 'None' : c.prop.type === 'turbofan' || c.prop.type === 'turbojet' ? 'Jet aeroplane' : 'Propeller aeroplane' }; },
  run(i0) {
    const T = i0.T_C + 273.15, c0 = Math.sqrt(1.4 * 287.05287 * T), i = { ...i0, c0, V: Math.min(i0.V, 0.85 * c0) }, M = i.V / c0, d2 = i.h ** 2 + i.y_side ** 2, d = Math.sqrt(d2), al = (f) => (i.absorb ? absorption(f, T, i.rh) : 0), warnings = [];
    // time history at the microphone: emission time τ from c(t − τ) = |x_obs − x_src(τ)| (uniform motion → quadratic)
    const span = M > 0.01 ? Math.min(60, (12 * d) / i.V) : 10, ts = N.linspace(-span + d / c0, span + d / c0, 241), tau = ts.map((t) => (M > 1e-6 ? (c0 * c0 * t - Math.sqrt(c0 ** 4 * t * t - (c0 * c0 - i.V ** 2) * (c0 * c0 * t * t - d2))) / (c0 * c0 - i.V ** 2) : t - d / c0));
    const hist = tau.map((ta) => { const x = i.V * ta, r = Math.sqrt(x * x + d2), cT = -x / r, f = i.f0 / (1 - M * cT); return { r, cT, f, L: flyLevel(i, r, cT, i.h / r, f, al(f)) }; });
    const LA = hist.map((h) => h.L), km = N.argmax(LA), LAmax = LA[km], dt = ts[1] - ts[0], sel = 10 * Math.log10(N.sum(LA.map((v) => 10 ** (v / 10))) * dt), above = LA.filter((v) => v >= LAmax - 10).length * dt;
    const fp = footprint(i, al(i.f0), Math.round(N.clamp(i.nGrid, 20, 300))), retErr = Math.max(...tau.map((ta, k) => Math.abs(c0 * (ts[k] - ta) - hist[k].r)));
    if (i0.V > i.V) warnings.push('Flight speed was limited to Mach 0.85: the moving-source solution here is subsonic.');
    if (M < 0.01) warnings.push('Hover: the level is constant in time and the exposure level depends only on the chosen duration window.');
    if (fp.area === 0) warnings.push(`The ${i.L_th} dBA contour does not reach the ground from the starting height.`);
    if (fp.track >= 30000) warnings.push('In level or shallow flight the contour is an open strip: the area is for the first 30 km of track only.');
    if (i.gamma_deg > 0 && M >= 0.01 && fp.track > 0) warnings.push('The footprint treats the whole source as one tone at the characteristic frequency; use the band spectrum for certification-type metrics.');
    // reality check: certified levels of real aircraft of similar mass
    const cert = certNeighbours(i.cert_class, i.mtom_kg), certKpi = [], certTab = [];
    if (cert) {
      const u = cert.set.unit, pt = cert.set.points, own = cert.key === 'rotor11' ? sel : LAmax;
      certKpi.push({ key: 'cert_level_ref', label: `Certified ${pt[cert.ref].toLowerCase()} level of comparable aircraft (median of ${cert.rows.length})`, value: cert.median, unit: u, note: `${cert.lo}–${cert.hi} ${u} for ${cert.rows[0][0]} and neighbours in mass; measured at the certification reference point, not at this microphone` });
      if (cert.direct) certKpi.push({ key: 'cert_level_delta_dB', label: `This prediction minus the certified ${cert.key === 'rotor11' ? 'exposure level' : 'maximum level'}`, value: own - cert.median, unit: 'dB', note: cert.key === 'rotor11' ? 'Comparable when the flyover is at 150 m in level flight at high speed' : 'Indicative only: the certification microphone is 2 500 m from brake release, under the climbing aircraft' });
      certTab.push({ title: `Certified noise levels of comparable aircraft: ${cert.set.metric}`, columns: ['Type', 'Engine', 'MTOM [kg]', 'Annex 16 chapter', ...pt.flatMap((q) => [`${q} [${u}]`, `${q} limit [${u}]`])], rows: cert.rows.map((r) => [r[0], r[1], r[2], r[cert.set.first - 1], ...r.slice(cert.set.first)]) });
    }
    return {
      kpis: [...certKpi,
        { key: 'SPL_peak_dBA', label: 'Maximum A-weighted level at the microphone', value: LAmax, unit: 'dBA' },
        { key: 'SEL_dBA', label: 'Sound exposure level', value: sel, unit: 'dBA', note: 'Event energy normalised to 1 s' },
        { key: 'duration_10dB_s', label: 'Time within 10 dB of the maximum', value: above, unit: 's' },
        { key: 'footprint_km2', label: `Ground area inside the ${i.L_th} dBA contour`, value: fp.area / 1e6, unit: 'km²' },
        { key: 'f_approach_Hz', label: 'Received frequency on approach', value: hist[0].f, unit: 'Hz' },
        { key: 'f_recede_Hz', label: 'Received frequency receding', value: hist[hist.length - 1].f, unit: 'Hz' },
        { key: 'doppler_max', label: 'Largest Doppler factor in the record', value: hist[0].f / i.f0, unit: '-', note: `limit 1/(1 − M) = ${(1 / (1 - M)).toFixed(4)}` },
        { key: 't_peak_s', label: 'Time of maximum relative to overhead passage', value: ts[km] - d / c0, unit: 's', note: 'Sound from overhead arrives d/c later; convective amplification moves the peak earlier' },
        { key: 'r_contour_m', label: 'Slant range of the contour level', value: fp.rTh, unit: 'm', note: 'Without absorption' },
        { key: 'mach', label: 'Flight Mach number', value: M, unit: '-' },
        { key: 'retarded_time_residual_m', label: 'Retarded-time equation residual', value: retErr, unit: 'm' },
      ],
      plots: [
        { type: 'line', title: 'Level at the microphone during the flyover', xlabel: 'Time after overhead passage [s]', ylabel: 'A-weighted level [dBA]', series: [{ name: 'Moving source', x: ts.map((t) => t - d / c0), y: LA }], annotations: [{ y: LAmax - 10, label: '10 dB down' }] },
        { type: 'line', title: 'Received frequency (Doppler shift)', xlabel: 'Time after overhead passage [s]', ylabel: 'Frequency [Hz]', series: [{ name: 'Received', x: ts.map((t) => t - d / c0), y: hist.map((h) => h.f) }], annotations: [{ y: i.f0, label: 'Emitted' }] },
        { type: 'heat', title: `Maximum level on the ground (half footprint, contour at ${i.L_th} dBA)`, xlabel: 'Distance along track [m]', ylabel: 'Sideline distance [m]', zlabel: 'L_Amax [dBA]', x: fp.x, y: fp.y, z: fp.z, contours: 12 },
      ],
      tables: certTab,
      warnings,
      models: ['Uniformly moving point source solved at the retarded time (far-field Ffowcs Williams–Hawkings monopole/dipole)', 'Convective amplification (1 − M cos θ)⁻² in pressure and Doppler shift f/(1 − M cos θ)', 'Spherical spreading and ISO 9613-1 absorption at the received frequency', 'Footprint from the maximum level over a straight climb track', ...(cert ? [`Measured reference: ${CERT_NOISE_ISSUE}`] : [])],
      assumptions: ['The source is represented by one overall level and one characteristic frequency', 'Straight flight path, homogeneous still air, flat ground, no ground-reflection correction', 'Footprint uses the emission geometry; propagation delay does not change the maximum level',
        ...(cert ? [cert.direct ? 'Certified levels are for real aircraft of similar mass at the Annex 16 reference conditions; they are a plausibility check, not a prediction for this design' : 'Certified levels are effective perceived noise levels (EPNdB) at the Annex 16 reference points: a tone- and duration-corrected metric that typically runs 10–15 dB above the maximum A-weighted level of the same event, so compare trends and orders of magnitude only'] : [])],
    };
  },
  convergence: { param: 'nGrid', label: 'Footprint grid points along track', levels: [40, 80, 160, 320], metric: 'footprint_km2' },
  verify() {
    const b = { L_ref: 90, r_ref: 100, f0: 500, directivity: 'Monopole (uniform)', V: 100, h: 200, y_side: 0, gamma_deg: 0, h0: 100, L_th: 70, T_C: 15, rh: 0.7, absorb: false, nGrid: 200 }, c0 = Math.sqrt(1.4 * 287.05287 * 288.15);
    const far = N.kv(flyover.run({ ...b, h: 5 })), hov = N.kv(flyover.run({ ...b, V: 0 })), rTh = 100 * 10 ** ((90 + aWeight(500) - 70) / 20), o = N.kv(flyover.run(b));
    return [
      N.check('Doppler factor approaches 1/(1 − M)', far.doppler_max, 1 / (1 - 100 / c0), 2e-3, 'Doppler relation for a source approaching head-on'),
      N.check('Approach × recede frequency product', far.f_approach_Hz * far.f_recede_Hz, 500 ** 2 / (1 - (100 / c0) ** 2), 4e-3, 'f²/(1 − M²) in the head-on limits'),
      N.check('Retarded-time equation satisfied', 1 + o.retarded_time_residual_m, 1, 1e-8, 'c(t − τ) = r(τ)'),
      N.check('Hover footprint area = π(r_c² − h²)', hov.footprint_km2, (PI * (rTh ** 2 - 100 ** 2)) / 1e6, 0.03, 'Spherical spreading; grid quadrature of the contour'),
      N.check('Hover level = L_ref − 20·log(h/r_ref) + A(f)', hov.SPL_peak_dBA, 90 - 20 * Math.log10(2) + aWeight(500), 1e-9, 'Inverse-square law'),
      N.check('Certified-noise lookup: lightest jet in the extract, lateral level', certNeighbours('Jet aeroplane', 2722, 1).rows[0][5], 81.3, 1e-12, 'EASA TCDSN jets issue 53: SF50 with FJ33-5A, 81.3 EPNdB lateral'),
      N.check('Certified-noise lookup: light helicopters use the Chapter 11 overflight exposure level', certNeighbours('Helicopter', 621, 1).median, 77.4, 1e-12, 'EASA TCDSN rotorcraft issue 52: R22 Beta, 77.4 dB(A) SEL'),
      N.check('Neighbours are ordered by mass ratio', (() => { const r = certNeighbours('Jet aeroplane', 78000, 5).rows.map((q) => Math.abs(Math.log(q[2] / 78000))); return r.every((v, k) => !k || v >= r[k - 1]) && r[4] < 0.2 ? 1 : 0; })(), 1, 1e-12, 'Sorting on |ln(MTOM ratio)|; five types within 20% of 78 t'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    out.push({ severity: 'info', title: 'Height is the cheapest noise reduction', detail: `L_Amax ${o.SPL_peak_dBA.toFixed(1)} dBA at ${i.h} m; the ${i.L_th} dBA footprint is ${o.footprint_km2.toFixed(3)} km².`, action: 'Doubling the height over the community lowers the maximum level by about 6 dB (more with absorption); steeper climb-out or approach shrinks the footprint roughly in proportion.', basis: 'Spherical spreading' });
    if (o.footprint_km2 > 1) out.push({ severity: 'advise', title: 'Large noise footprint', detail: `${o.footprint_km2.toFixed(2)} km² inside ${i.L_th} dBA.`, action: 'Evaluate a steeper initial climb, reduced-thrust or reduced-rotor-speed procedures and track changes away from housing; lower source noise usually also means lower tip or jet speeds and lower energy use.', basis: 'Contour area' });
    if (Number.isFinite(o.cert_level_ref)) { const kc = res.kpis.find((k) => k.key === 'cert_level_ref'); out.push({ severity: Number.isFinite(o.cert_level_delta_dB) && Math.abs(o.cert_level_delta_dB) > 10 ? 'advise' : 'info', title: 'Compare with certified aircraft of this size', detail: `Aircraft near ${i.mtom_kg.toFixed(0)} kg: ${kc.label.toLowerCase()} ${o.cert_level_ref} ${kc.unit} (see the table). This prediction gives L_Amax ${o.SPL_peak_dBA.toFixed(1)} dBA and SEL ${o.SEL_dBA.toFixed(1)} dBA at ${i.h} m${Number.isFinite(o.cert_level_delta_dB) ? `, ${o.cert_level_delta_dB >= 0 ? '+' : ''}${o.cert_level_delta_dB.toFixed(1)} dB against the certified figure` : ''}.`, action: 'If the prediction is far from real aircraft of the same class, calibrate the source constants in the source-breakdown analysis before using the footprint; sources not modelled (fan, core, tail rotor, blade–vortex interaction) make the prediction a lower bound.', basis: 'EASA certification noise levels (type-certificate data sheets for noise)' }); }
    if (o.duration_10dB_s > 20) out.push({ severity: 'info', title: 'Long event duration', detail: `${o.duration_10dB_s.toFixed(0)} s within 10 dB of the peak raises the exposure level to ${o.SEL_dBA.toFixed(1)} dBA.`, action: 'Slow, low overflights are judged by exposure, not peak level: trade speed against height.', basis: 'Sound exposure level' });
    return out;
  },
};

// ---- cabin noise by mass-law transmission -------------------------------------------------------------------
/** Field-incidence transmission loss of a single or double wall in each one-third-octave band. */
export function wallTL(f, m1, m2, d, fc, eta, rc, c) {
  const tl = (m) => 10 * Math.log10(1 + ((PI * f * m) / rc) ** 2) - 5, coin = (m) => (f > fc ? 20 * Math.log10((PI * f * m) / rc) + 10 * Math.log10((2 * eta * f) / (PI * fc)) : tl(m));
  if (!(m2 > 0)) return Math.max(Math.min(tl(m1), coin(m1)), 0);
  const f0 = (1 / (2 * PI)) * Math.sqrt(((rc * c) / d) * (1 / m1 + 1 / m2)), fl = c / (2 * PI * d), t1 = Math.max(Math.min(tl(m1), coin(m1)), 0), t2 = Math.max(tl(m2), 0);
  return f < f0 ? Math.max(tl(m1 + m2), 0) : f < fl ? Math.max(t1 + t2 + 20 * Math.log10((4 * PI * f * d) / c), tl(m1 + m2)) : t1 + t2 + 6;
}
const cabin = {
  id: 'cabin', title: 'Cabin noise: exterior excitation and wall transmission', fidelity: 'reduced-order',
  summary: 'Estimates the interior level from boundary-layer pressure fluctuations and an exterior tone on the fuselage, the transmission loss of the skin and trim panel, and the absorption in the cabin.',
  equations: ['Acoustic intensity equations', 'Acoustic energy equations', 'Cabin acoustic transmission models', 'Structural vibration–acoustic coupling'],
  inputs: [
    { key: 'V', label: 'Cruise true airspeed', unit: 'm/s', default: 140, min: 1, group: 'Exterior' },
    { key: 'alt_m', label: 'Cruise altitude', unit: 'm', default: 6000, min: 0, max: 15000, group: 'Exterior' },
    { key: 'x_m', label: 'Distance from the nose to the cabin station', unit: 'm', default: 8, min: 0.1, group: 'Exterior', help: 'Sets the boundary-layer thickness and its peak frequency' },
    { key: 'L_tone', label: 'Exterior tone level on the fuselage', unit: 'dB', default: 0, min: 0, max: 160, group: 'Exterior', help: 'Propeller or rotor blade-passing tone at the skin in cruise; 0 for none. The default extrapolates the far-field tone to 1.5 radii from the hub and is a rough estimate: enter a measured or near-field value when available' },
    { key: 'f_tone', label: 'Tone frequency', unit: 'Hz', default: 100, min: 10, max: 5000, group: 'Exterior' },
    { key: 'material', label: 'Skin material', type: 'select', options: Object.keys(METALS), default: 'Al 2024-T3', group: 'Wall' },
    { key: 't_skin_mm', label: 'Skin thickness', unit: 'mm', default: 1.6, min: 0.3, max: 10, group: 'Wall' },
    { key: 'm_trim', label: 'Trim panel surface density', unit: 'kg/m²', default: 1.5, min: 0, max: 20, group: 'Wall', help: '0 for a bare single wall' },
    { key: 'gap_m', label: 'Insulation gap', unit: 'm', default: 0.08, min: 0.005, max: 0.3, group: 'Wall' },
    { key: 'eta', label: 'Skin damping loss factor', unit: '-', default: 0.02, min: 0.001, max: 0.3, group: 'Wall', help: '0.01 bare metal, 0.05–0.15 with damping tiles' },
    { key: 'alpha_cabin', label: 'Mean cabin absorption coefficient', unit: '-', default: 0.3, min: 0.02, max: 0.9, group: 'Cabin', help: '0.1 bare, 0.3 furnished with seats and carpet, 0.5 well treated' },
    { key: 'cabin_alt_m', label: 'Cabin altitude', unit: 'm', default: 2400, min: 0, max: 5000, group: 'Cabin' },
  ],
  defaults: (c, up, d) => {
    // exterior tone on the skin: the blade-passing tone of the tonal analysis brought in to 1.5 radii from the hub; for a propeller the tone was
    // computed at static take-off thrust, so the loading noise is scaled to the cruise thrust (level ∝ thrust)
    const V = c.mission.cruise_V_ms || c.flight.V_ms, alt = c.mission.cruise_alt_m || c.atm.alt_m, prop = !(c.rotor.R_m > 0) && c.wing.S_m2 > 0 && c.prop.T0_N > 0, qS = 0.5 * isa(alt, c.atm.dISA_K).rho * V * V * c.wing.S_m2;
    const Tcr = prop && qS > 0 ? (qS * (c.aero.CD0 + d.k_induced * (d.W / qS) ** 2)) / Math.max(1, c.prop.n_eng) : 0, dT = Tcr > 0 ? Math.min(0, 20 * Math.log10(Tcr / c.prop.T0_N)) : 0;
    return { V, alt_m: alt, x_m: Math.max(0.3 * c.fuselage.len_m, 0.2) || undefined, L_tone: hasRotor(c) && up.acoustics?.SPL_bpf_dB ? up.acoustics.SPL_bpf_dB + dT + 20 * Math.log10((c.rotor.R_m > 0 ? 150 : 300) / Math.max(1.5 * (c.rotor.R_m || c.prop.prop_dia_m / 2), 0.5)) : undefined, f_tone: up.acoustics?.bpf_Hz, material: c.struct.material in METALS ? c.struct.material : undefined, t_skin_mm: Math.min(c.struct.t_skin_mm, 2.5) || undefined, cabin_alt_m: c.fuselage.cabin_dp_Pa > 0 ? c.systems.cabin_alt_m : Math.min(c.mission.cruise_alt_m || c.atm.alt_m, 5000), m_trim: c.mission.pax > 0 ? undefined : 0 }; },
  run(i) {
    const a = isa(i.alt_m), ci = isa(i.cabin_alt_m), mat = METALS[i.material] || METALS['Al 2024-T3'], t = i.t_skin_mm / 1e3, m1 = mat.rho * t, M = i.V / a.a, q = 0.5 * a.rho * i.V ** 2;
    const prms = (0.006 * q) / (1 + 0.14 * M * M), Rex = Math.max((i.V * i.x_m) / a.nu, 1e4), dstar = (0.046 * i.x_m) / Rex ** 0.2, fp = (0.1 * i.V) / dstar, ext = hump(fp).map((w) => dB(prms) + 10 * Math.log10(w));
    if (i.L_tone > 0) { const b = bandOf(i.f_tone); ext[b] = sumdB([ext[b], i.L_tone]); }
    const rc = 0.5 * (a.rho * a.a + ci.rho * ci.a), cL = Math.sqrt(mat.E / (mat.rho * (1 - mat.nu ** 2))), fc = (ci.a ** 2 * Math.sqrt(3)) / (PI * cL * t), tl = BANDS.map((f) => wallTL(f, m1, i.m_trim, i.gap_m, fc, i.eta, rc, ci.a));
    const room = 10 * Math.log10(0.8 / i.alpha_cabin), inn = ext.map((v, k) => v - tl[k] + room), dBAi = sumdB(inn.map((v, k) => v + aWeight(BANDS[k]))), sil = speechInterference(inn), f0 = i.m_trim > 0 ? (1 / (2 * PI)) * Math.sqrt(((rc * ci.a) / i.gap_m) * (1 / m1 + 1 / i.m_trim)) : NaN;
    const warnings = ['Idealised infinite-panel transmission: frames, stringers, windows, structure-borne paths and leaks are not represented, so real cabins are typically several decibels louder than this estimate.'];
    if (i.L_tone > 0 && i.m_trim > 0 && Math.abs(Math.log2(i.f_tone / f0)) < 0.35) warnings.push(`The exterior tone at ${i.f_tone.toFixed(0)} Hz is close to the double-wall resonance at ${f0.toFixed(0)} Hz, where the trim panel gives no benefit.`);
    if (!(i.m_trim > 0)) warnings.push('Single wall (no trim panel): the double-wall resonance is undefined.');
    if (i.L_tone > 0 && i.f_tone < 250) warnings.push('Low-frequency tones are only weakly attenuated by a light wall: mass law gives about 6 dB per doubling of frequency or mass.');
    const kb = bandOf(i.f_tone);
    return {
      kpis: [
        { key: 'cabin_SPL_dBA', label: 'Cabin level', value: dBAi, unit: 'dBA', status: dBAi < 85 ? 'ok' : dBAi < 95 ? 'warn' : 'bad', note: 'Guide: 75–85 dBA in airliner cabins in cruise; 80 and 85 dBA are the occupational action values for an 8-hour day, above which crews need hearing protection' },
        { key: 'cabin_OASPL_dB', label: 'Cabin overall level', value: sumdB(inn), unit: 'dB' },
        { key: 'exterior_OASPL_dB', label: 'Exterior level on the skin', value: sumdB(ext), unit: 'dB' },
        { key: 'SIL_dB', label: 'Speech interference level', value: sil, unit: 'dB', status: sil < 60 ? 'ok' : sil < 70 ? 'warn' : 'bad', note: 'Mean of the 500 Hz, 1, 2 and 4 kHz octave-band levels; below about 60 dB for normal conversation' },
        { key: 'TL_500_dB', label: 'Wall transmission loss at 500 Hz', value: tl[14], unit: 'dB' },
        { key: 'TL_tone_dB', label: 'Wall transmission loss at the tone', value: tl[kb], unit: 'dB' },
        { key: 'tone_inside_dB', label: 'Tone level inside', value: i.L_tone > 0 ? i.L_tone - tl[kb] + room : 0, unit: 'dB' },
        { key: 'f_mass_air_mass_Hz', label: 'Double-wall resonance', value: f0, unit: 'Hz' },
        { key: 'f_coincidence_Hz', label: 'Skin coincidence frequency', value: fc, unit: 'Hz' },
        { key: 'skin_surface_density', label: 'Skin surface density', value: m1, unit: 'kg/m²' },
        { key: 'p_rms_tbl_Pa', label: 'Boundary-layer wall-pressure fluctuation', value: prms, unit: 'Pa' },
      ],
      plots: [
        { type: 'line', title: 'Exterior and cabin spectra', xlabel: 'Band centre frequency [Hz]', ylabel: 'Band level [dB re 20 µPa]', xlog: true, series: [{ name: 'On the fuselage skin', x: BANDS, y: ext }, { name: 'In the cabin', x: BANDS, y: inn }] },
        { type: 'line', title: 'Wall transmission loss', xlabel: 'Band centre frequency [Hz]', ylabel: 'Transmission loss [dB]', xlog: true, series: [{ name: i.m_trim > 0 ? 'Skin + gap + trim' : 'Skin only', x: BANDS, y: tl }, { name: 'Skin alone (mass law)', x: BANDS, y: BANDS.map((f) => wallTL(f, m1, 0, i.gap_m, fc, i.eta, rc, ci.a)), style: 'dash' }], annotations: Number.isFinite(f0) ? [{ x: f0, label: 'Mass–air–mass' }] : [] },
      ],
      warnings,
      models: ['Turbulent-boundary-layer wall-pressure level p_rms = 0.006·q/(1 + 0.14·M²) (Lowson, empirical) with a generic spectrum', 'Mass-law transmission loss, field incidence (normal incidence − 5 dB)', 'Double-wall mass–air–mass resonance and coincidence dip (infinite-panel theory)', 'Diffuse-field receiving-room relation L_in = L_ext − TL + 10·log(S/A)'],
      assumptions: ['Wall area is 80% of the cabin surface; uniform absorption', 'Characteristic impedance is the mean of the outside and cabin air', 'Airborne path only; engine, ECS and structure-borne noise are not included', 'Skin material properties from the shared materials table', 'Trim mass, insulation gap, damping and absorption defaults are typical values, not data for a specific cabin; a cabin altitude of 2 400 m is the 8 000 ft maximum of CS/FAR 25.841(a)'],
    };
  },
  calibration: { params: [{ key: 'alpha_cabin', min: 0.05, max: 0.8 }, { key: 'eta', min: 0.005, max: 0.2 }, { key: 'm_trim', min: 0, max: 10 }], sweep: 'V', target: 'cabin_SPL_dBA', note: 'Supply measured cabin level against airspeed (or a measured transmission-loss curve); cabin absorption, damping and effective trim mass absorb the idealisations of the wall model.' },
  verify() {
    const rc = 415, f = 2000, a = wallTL(f, 5, 0, 0.1, 1e9, 0.02, rc, 340), b = wallTL(f, 10, 0, 0.1, 1e9, 0.02, rc, 340), o = N.kv(cabin.run({ V: 100, alt_m: 0, x_m: 5, L_tone: 0, f_tone: 100, material: 'Al 2024-T3', t_skin_mm: 2, m_trim: 2, gap_m: 0.1, eta: 0.02, alpha_cabin: 0.3, cabin_alt_m: 0 })), r0 = isa(0);
    return [
      N.check('Normal-incidence mass law', a + 5, 10 * Math.log10(1 + ((PI * f * 5) / rc) ** 2), 1e-12, 'Exact limp-wall transmission, τ = 1/(1 + (πfm/ρc)²)'),
      N.check('Doubling the mass adds 6.02 dB', b - a, 20 * Math.log10(2), 1e-4, 'Mass law (high-frequency limit)'),
      N.check('Mass–air–mass resonance', o.f_mass_air_mass_Hz, Math.sqrt(((r0.rho * r0.a ** 2) / 0.1) * (1 / (2780 * 0.002) + 1 / 2)) / (2 * PI), 1e-10, 'f0 = (1/2π)·√(ρc²/d·(1/m1 + 1/m2))'),
      N.check('Coincidence frequency of a 2 mm aluminium skin', o.f_coincidence_Hz, (r0.a ** 2 * Math.sqrt(3)) / (PI * Math.sqrt(73.1e9 / (2780 * (1 - 0.33 ** 2))) * 0.002), 1e-10, 'Thin-plate bending wave speed equal to the speed of sound'),
      N.check('Speech interference level of a flat one-third-octave spectrum: band level + 10·log 3', speechInterference(BANDS.map(() => 60)), 60 + 10 * Math.log10(3), 1e-12, 'Octave band = energy sum of its three one-third-octave bands'),
      N.check('Boundary-layer pressure level', o.p_rms_tbl_Pa, (0.006 * 0.5 * r0.rho * 1e4) / (1 + 0.14 * (100 / r0.a) ** 2), 1e-12, 'Lowson correlation as implemented'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.cabin_SPL_dBA > 80) out.push({ severity: o.cabin_SPL_dBA > 85 ? 'warn' : 'advise', title: o.cabin_SPL_dBA > 85 ? 'Cabin is loud' : 'Cabin level is at the upper end of normal', detail: `${o.cabin_SPL_dBA.toFixed(0)} dBA, speech interference level ${o.SIL_dB.toFixed(0)} dB${o.cabin_SPL_dBA > 85 ? ': above the 85 dBA action value for a working day, so crews need hearing protection or headsets' : ' (the estimate carries several decibels of uncertainty)'}.`, action: 'Add skin damping, increase the insulation gap, or raise trim mass in the loudest zone only; every kilogram of treatment costs fuel for the life of the aircraft, so target the dominant band first.', basis: 'Mass law and double-wall theory; occupational noise action values of 80 and 85 dB(A) for an 8-hour day (EU Directive 2003/10/EC)' });
    if (i.L_tone > 0 && o.tone_inside_dB > o.cabin_OASPL_dB - 3) out.push({ severity: 'advise', title: 'A propeller or rotor tone dominates the cabin', detail: `${o.tone_inside_dB.toFixed(0)} dB at ${i.f_tone.toFixed(0)} Hz with only ${o.TL_tone_dB.toFixed(0)} dB of wall attenuation.`, action: 'Low-frequency tones are better treated at source (tip clearance, blade count, synchrophasing) or with tuned vibration absorbers and active noise control than with added mass.', basis: 'Mass law: 6 dB per octave' });
    if (Number.isFinite(o.f_mass_air_mass_Hz) && o.f_mass_air_mass_Hz > 150) out.push({ severity: 'info', title: 'Double-wall resonance is in the audible working range', detail: `${o.f_mass_air_mass_Hz.toFixed(0)} Hz.`, action: 'A deeper gap or heavier trim lowers it; keep it at least an octave below the blade-passing frequency.', basis: 'Mass–air–mass resonance' });
    return out;
  },
};

// ---- duct acoustics: finite-difference Helmholtz solver -----------------------------------------------------
/**
 * 1-D Helmholtz (Webster horn) equation d/dx(S dp/dx) + k² S p = 0 by second-order finite differences, unit incident
 * wave at the inlet and an anechoic outlet. S(x) is the duct area. Returns transmitted and reflected amplitudes.
 */
export function helmholtzDuct(S, len, k, n) {
  const h = len / n, ar = new Float64Array(n + 1), br = new Float64Array(n + 1), bi = new Float64Array(n + 1), cr = new Float64Array(n + 1), dr = new Float64Array(n + 1), di = new Float64Array(n + 1);
  for (let j = 0; j <= n; j++) { const sm = S((j - 0.5) * h), sp = S((j + 0.5) * h); ar[j] = sm; cr[j] = sp; br[j] = -(sm + sp) + k * k * h * h * 0.5 * (sm + sp); }
  // radiation conditions by ghost-node elimination: inlet dp/dx − ik p = −2ik, outlet dp/dx = −ik p
  const s0 = S(0), sN = S(len); ar[0] = 0; cr[0] = 2 * s0; br[0] = s0 * (-2 + k * k * h * h); bi[0] = -2 * h * k * s0; di[0] = -4 * h * k * s0; cr[n] = 0; ar[n] = 2 * sN; br[n] = sN * (-2 + k * k * h * h); bi[n] = -2 * h * k * sN;
  // complex Thomas algorithm
  const cpr = new Float64Array(n + 1), cpi = new Float64Array(n + 1), dpr = new Float64Array(n + 1), dpi = new Float64Array(n + 1), pr = new Float64Array(n + 1), pi = new Float64Array(n + 1);
  let mr = br[0], mi = bi[0], den = mr * mr + mi * mi; cpr[0] = (cr[0] * mr) / den; cpi[0] = (-cr[0] * mi) / den; dpr[0] = (dr[0] * mr + di[0] * mi) / den; dpi[0] = (di[0] * mr - dr[0] * mi) / den;
  for (let j = 1; j <= n; j++) { mr = br[j] - ar[j] * cpr[j - 1]; mi = bi[j] - ar[j] * cpi[j - 1]; den = mr * mr + mi * mi; const nr = dr[j] - ar[j] * dpr[j - 1], ni = di[j] - ar[j] * dpi[j - 1]; cpr[j] = (cr[j] * mr) / den; cpi[j] = (-cr[j] * mi) / den; dpr[j] = (nr * mr + ni * mi) / den; dpi[j] = (ni * mr - nr * mi) / den; }
  pr[n] = dpr[n]; pi[n] = dpi[n];
  for (let j = n - 1; j >= 0; j--) { pr[j] = dpr[j] - (cpr[j] * pr[j + 1] - cpi[j] * pi[j + 1]); pi[j] = dpi[j] - (cpr[j] * pi[j + 1] + cpi[j] * pr[j + 1]); }
  return { T: Math.hypot(pr[n], pi[n]), R: Math.hypot(pr[0] - 1, pi[0]), amp: Array.from(pr, (v, j) => Math.hypot(v, pi[j])), h };
}
const duct = {
  id: 'duct', title: 'Duct silencer: finite-difference Helmholtz solver', fidelity: 'numerical',
  summary: 'Solves the acoustic wave equation in the frequency domain along an exhaust or ventilation duct with an expansion chamber, giving the transmission loss against frequency, compared with the exact plane-wave solution.',
  equations: ['Helmholtz equation', 'Acoustic wave equation', 'Acoustic energy equations', 'Acoustic finite element models'],
  inputs: [
    { key: 'D_pipe', label: 'Pipe diameter', unit: 'm', default: 0.06, min: 0.005, group: 'Geometry' },
    { key: 'D_chamber', label: 'Chamber diameter', unit: 'm', default: 0.18, min: 0.005, group: 'Geometry' },
    { key: 'L_chamber', label: 'Chamber length', unit: 'm', default: 0.5, min: 0.02, group: 'Geometry' },
    { key: 'T_gas_K', label: 'Gas temperature', unit: 'K', default: 600, min: 200, max: 1400, group: 'Gas', help: 'About 600–900 K in a piston or APU exhaust, ambient in a ventilation duct' },
    { key: 'f_target', label: 'Tone to be attenuated', unit: 'Hz', default: 180, min: 5, max: 5000, group: 'Target', help: 'Engine firing frequency or blade-passing frequency' },
    { key: 'f_max', label: 'Highest frequency of the sweep', unit: 'Hz', default: 1500, min: 50, max: 8000, group: 'Numerics' },
    { key: 'nNodes', label: 'Grid intervals along the duct', unit: '', default: 240, min: 30, max: 4000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up) => ({ f_target: [up.acoustics?.bpf_Hz, c.prop.n_blades > 0 && c.prop.rpm > 0 ? (c.prop.n_blades * c.prop.rpm) / 60 : 0].find((f) => f > 40), T_gas_K: c.prop.type === 'electric' ? 310 : c.prop.type === 'piston' ? 800 : 700, D_pipe: c.prop.type === 'piston' ? N.clamp(0.02 * (c.prop.P0_W / 1e3) ** 0.25 * 1.0, 0.015, 0.15) : undefined }),
  run(i) {
    const c = Math.sqrt(1.4 * 287.05287 * i.T_gas_K), S1 = (PI * i.D_pipe ** 2) / 4, S2 = (PI * i.D_chamber ** 2) / 4, m = S2 / S1, Lp = i.L_chamber, len = 3 * Lp, n = 3 * Math.max(10, Math.round(i.nNodes / 3));
    const S = (x) => (x > Lp && x < 2 * Lp ? S2 : S1), tlN = (f) => -20 * Math.log10(helmholtzDuct(S, len, (2 * PI * f) / c, n).T), tlA = (f) => 10 * Math.log10(1 + 0.25 * (m - 1 / m) ** 2 * Math.sin((2 * PI * f * Lp) / c) ** 2);
    const fs = N.linspace(i.f_max / 150, i.f_max, 150), num = fs.map(tlN), ana = fs.map(tlA), err = Math.max(...num.map((v, k) => Math.abs(v - ana[k]))), tT = tlN(i.f_target), fcut = (1.8412 * c) / (PI * Math.max(i.D_chamber, i.D_pipe)), sol = helmholtzDuct(S, len, (2 * PI * i.f_target) / c, n), warnings = [];
    if (i.f_max > fcut) warnings.push(`Above ${fcut.toFixed(0)} Hz the first transverse mode of the chamber propagates and the plane-wave (1-D) model over-predicts attenuation.`);
    if (((2 * PI * i.f_max) / c) * sol.h > 0.5) warnings.push('Fewer than about 12 grid points per wavelength at the top of the sweep: increase the grid intervals.');
    if (tT < 5) warnings.push(`Only ${tT.toFixed(1)} dB at the target tone: it lies near a pass band (chamber length close to a multiple of half a wavelength).`);
    return {
      kpis: [
        { key: 'TL_target_dB', label: 'Transmission loss at the target tone', value: tT, unit: 'dB', status: tT > 10 ? 'ok' : tT > 5 ? 'warn' : 'bad' },
        { key: 'TL_target_exact_dB', label: 'Exact plane-wave value', value: tlA(i.f_target), unit: 'dB' },
        { key: 'TL_peak_dB', label: 'Peak transmission loss of the chamber', value: 10 * Math.log10(1 + 0.25 * (m - 1 / m) ** 2), unit: 'dB' },
        { key: 'f_first_peak_Hz', label: 'First attenuation peak', value: c / (4 * Lp), unit: 'Hz', note: 'Chamber a quarter wavelength long' },
        { key: 'f_first_pass_Hz', label: 'First pass band', value: c / (2 * Lp), unit: 'Hz', note: 'No attenuation' },
        { key: 'L_optimum_m', label: 'Chamber length for peak loss at the target', value: c / (4 * i.f_target), unit: 'm' },
        { key: 'area_ratio', label: 'Expansion area ratio', value: m, unit: '-' },
        { key: 'f_cuton_Hz', label: 'Plane-wave validity limit', value: fcut, unit: 'Hz' },
        { key: 'max_error_dB', label: 'Largest difference from the exact solution', value: err, unit: 'dB', status: err < 0.5 ? 'ok' : 'warn' },
        { key: 'energy_balance', label: 'Reflected + transmitted power', value: sol.R ** 2 + sol.T ** 2, unit: '-', note: 'Should equal 1 for a loss-free duct' },
        { key: 'c_gas_ms', label: 'Speed of sound in the gas', value: c, unit: 'm/s' },
      ],
      plots: [
        { type: 'line', title: 'Transmission loss of the expansion chamber', xlabel: 'Frequency [Hz]', ylabel: 'Transmission loss [dB]', series: [{ name: 'Finite-difference Helmholtz', x: fs, y: num }, { name: 'Exact plane-wave theory', x: fs, y: ana, style: 'dash' }], annotations: [{ x: i.f_target, label: 'Target tone' }] },
        { type: 'line', title: `Pressure amplitude along the duct at ${i.f_target.toFixed(0)} Hz`, xlabel: 'Distance along the duct [m]', ylabel: '|p| / incident amplitude [-]', series: [{ name: 'Standing-wave pattern', ...(([x, y]) => ({ x, y }))(((x, y) => { const s = Math.max(1, Math.ceil(x.length / 300)); return [x.filter((_, k) => k % s === 0), y.filter((_, k) => k % s === 0)]; })(N.range(n + 1, (j) => j * sol.h), sol.amp)) }], annotations: [{ x: Lp, label: 'Chamber inlet' }, { x: 2 * Lp, label: 'Chamber outlet' }] },
      ],
      warnings,
      models: ['Frequency-domain Helmholtz (Webster) equation, second-order conservative finite differences', 'Non-reflecting inlet and outlet conditions with a unit incident wave', 'Complex tridiagonal (Thomas) solver', 'Exact transfer-matrix solution of the simple expansion chamber (reference)'],
      assumptions: ['Plane waves, rigid walls, no mean flow, no absorption or viscous losses', 'Uniform gas temperature', 'Inlet and outlet pipes each one chamber length long, both anechoically terminated'],
    };
  },
  convergence: { param: 'nNodes', label: 'Grid intervals along the duct', levels: [60, 120, 240, 480, 960], metric: 'TL_target_dB' },
  verify() {
    const b = { D_pipe: 0.05, D_chamber: 0.2, L_chamber: 0.4, T_gas_K: 300, f_target: 217, f_max: 800, nNodes: 600 }, o = N.kv(duct.run(b)), st = helmholtzDuct(() => 1, 1, 5, 200), c = Math.sqrt(1.4 * 287.05287 * 300);
    return [
      N.check('Transmission loss at the target equals the exact solution', o.TL_target_dB, o.TL_target_exact_dB, 5e-3, 'TL = 10·log[1 + ¼(m − 1/m)² sin²(kL)]; error O((kh)²)'),
      N.check('Largest error across the sweep below 0.2 dB', o.max_error_dB < 0.2 ? 1 : 0, 1, 1e-12, 'Second-order dispersion error'),
      N.check('Acoustic power is conserved (R² + T² = 1)', o.energy_balance, 1, 1e-3, 'Loss-free duct with equal inlet and outlet areas'),
      N.check('Uniform pipe transmits completely', st.T, 1, 1e-4, 'Travelling plane wave'),
      N.check('Quarter-wave peak frequency', o.f_first_peak_Hz, c / 1.6, 1e-12, 'c/(4L)'),
    ];
  },
  validation: [{ name: 'Expansion chamber, area ratio 16, peak transmission loss', source: 'Plane-wave transfer-matrix theory (Davis et al., NACA Report 1192): TL_max = 10·log[1 + ¼(m − 1/m)²] = 18.1 dB', inputs: { D_pipe: 0.05, D_chamber: 0.2, L_chamber: 0.4, T_gas_K: 300, f_max: 800, nNodes: 600 }, sweep: { key: 'f_target', values: [217.0] }, target: 'TL_target_dB', observed: [18.1], tol_pct: 1 }],
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.TL_target_dB < 10) out.push({ severity: 'advise', title: 'Chamber is poorly tuned to the target tone', detail: `${o.TL_target_dB.toFixed(1)} dB at ${i.f_target.toFixed(0)} Hz against a possible ${o.TL_peak_dB.toFixed(1)} dB.`, action: `Set the chamber length near ${o.L_optimum_m.toFixed(2)} m (a quarter wavelength in the hot gas), or an odd multiple, and increase the area ratio if more loss is needed.`, basis: 'Quarter-wave expansion chamber' });
    out.push({ severity: 'info', title: 'Reactive silencers cost little pressure loss', detail: `Area ratio ${o.area_ratio.toFixed(1)} gives up to ${o.TL_peak_dB.toFixed(1)} dB without flow restriction.`, action: 'Prefer a tuned chamber to a restrictive baffle: back-pressure costs engine power and fuel. Add an absorptive lining for the broadband range above the plane-wave limit.', basis: 'Plane-wave duct acoustics' });
    if (o.max_error_dB > 0.5) out.push({ severity: 'warn', title: 'Numerical dispersion', detail: `Finite-difference result departs ${o.max_error_dB.toFixed(2)} dB from the exact solution.`, action: 'Increase the grid intervals (at least 20 points per wavelength at the highest frequency); run the convergence study.', basis: 'Second-order dispersion error' });
    return out;
  },
};

export default {
  id: 'acoustics', n: 11,
  tagline: 'How loud the aircraft is, which source is responsible, what people on the ground and in the cabin hear, and what would make it quieter.',
  analyses: [tonal, spectrum, flyover, cabin, duct],
  consumes: [
    { from: 'rotorcraft', keys: ['hover_power_W'], why: 'Rotor torque for loading noise' },
    { from: 'performance', keys: ['V_2_ms', 'roc_sl_ms'], why: 'Take-off safety speed and sea-level climb rate for the source speed and the footprint climb angle' },
  ],
  provides: [
    { key: 'OASPL_dB', label: 'Overall sound pressure level', unit: 'dB' }, { key: 'SPL_peak_dBA', label: 'Maximum flyover level', unit: 'dBA' }, { key: 'bpf_Hz', label: 'Blade-passing frequency', unit: 'Hz' },
    { key: 'cabin_SPL_dBA', label: 'Cabin level', unit: 'dBA' }, { key: 'footprint_km2', label: 'Noise footprint area', unit: 'km²' },
  ],
  handoff: [
    { model: 'Computational aeroacoustics: linearised Euler / Navier–Stokes and LES coupled to permeable-surface FW-H or Kirchhoff integrals', why: 'Needs time-resolved 3-D flow data (terabytes) from scale-resolving CFD; only compact analytical sources are integrated here', tool: 'LES/DES solver with an FW-H post-processor (PSU-WOPWOP, OpenFOAM libAcoustics, Actran, PowerFLOW)' },
    { model: 'Blade–vortex interaction, unsteady loading harmonics and edgewise-flight rotor noise', why: 'Requires the free wake and blade airloads at every azimuth; this suite uses steady axial loading', tool: 'Comprehensive rotor code + acoustic solver (CAMRAD II or CHARM with PSU-WOPWOP)' },
    { model: 'Turbofan fan, core, turbine and combustion noise with duct modes and liners', why: 'Rotor–stator interaction, cut-off duct modes and liner impedance need engine geometry not in the case', tool: 'NASA ANOPP / ANOPP2, engine-manufacturer methods' },
    { model: 'Certification metrics: tone-corrected perceived noise level and EPNL', why: 'Requires the full noy tables, tone corrections and a band time history at the certification points', tool: 'ANOPP2, SAE ARP 866/1845 procedures, FAA AEDT for contours' },
    { model: 'Refraction by wind and temperature gradients, terrain and ground impedance', why: 'A homogeneous still atmosphere and flat ground are assumed', tool: 'Ray-tracing or parabolic-equation propagation codes; AEDT / IMPACT for airport contours' },
    { model: 'Finite-element / boundary-element cabin acoustics and structural–acoustic coupling', why: 'Needs the fuselage structural model and cabin cavity mesh; infinite-panel transmission is used here', tool: 'Vibro-acoustic FE/BE and statistical energy analysis (Actran, VA One, COMSOL)' },
    { model: 'Broadband self-noise by turbulence-resolved source models (full BPM, Amiet with measured spectra)', why: 'Depends on section boundary-layer data at every radius and calibrated spectral functions; an empirical overall-level law is used', tool: 'NAFNoise, ANOPP2 self-noise modules, acoustic wind-tunnel tests' },
  ],
};
