// International Standard Atmosphere (ISO 2533 / US 1976) to 47 km with temperature offset,
// plus air transport properties and airspeed conversions shared by every suite.

export const G0 = 9.80665, R_AIR = 287.05287, GAMMA = 1.4, P0 = 101325, T0 = 288.15, RHO0 = 1.225, A0 = 340.294;
const LAYERS = [
  [0, -0.0065, 288.15, 101325],
  [11000, 0, 216.65, 22632.06],
  [20000, 0.001, 216.65, 5474.889],
  [32000, 0.0028, 228.65, 868.0187],
  [47000, 0, 270.65, 110.9063],
];

/**
 * Atmospheric state at geopotential altitude h [m] with ISA temperature deviation dT [K].
 * Returns {T, p, rho, a, mu, nu, sigma, delta, theta}.
 */
export function isa(h, dT = 0) {
  h = Math.max(-2000, Math.min(47000, h));
  let L = LAYERS[0];
  for (const l of LAYERS) if (h >= l[0]) L = l;
  const [hb, lapse, Tb, pb] = L;
  const Tstd = Tb + lapse * (h - hb);
  const p = lapse === 0 ? pb * Math.exp((-G0 * (h - hb)) / (R_AIR * Tb)) : pb * (Tstd / Tb) ** (-G0 / (lapse * R_AIR));
  const T = Tstd + dT, rho = p / (R_AIR * T), a = Math.sqrt(GAMMA * R_AIR * T), mu = sutherland(T);
  return { T, p, rho, a, mu, nu: mu / rho, sigma: rho / RHO0, delta: p / P0, theta: T / T0 };
}
/** Dynamic viscosity of air by Sutherland's law [Pa s]. */
export const sutherland = (T) => (1.458e-6 * T ** 1.5) / (T + 110.4);
/** Thermal conductivity of air [W/m/K] (Sutherland-type fit). */
export const kAir = (T) => (2.64638e-3 * T ** 1.5) / (T + 245.4 * 10 ** (-12 / T));
export const CP_AIR = 1004.7, PR_AIR = 0.71;

/** Pressure altitude [m] from static pressure [Pa] (troposphere/stratosphere). */
export function pressureAltitude(p) {
  if (p > 22632.06) return (T0 / 0.0065) * (1 - (p / P0) ** ((0.0065 * R_AIR) / G0));
  return 11000 - ((R_AIR * 216.65) / G0) * Math.log(p / 22632.06);
}
/** Density altitude [m] from actual density. */
export function densityAltitude(rho) {
  let lo = -2000, hi = 30000;
  for (let i = 0; i < 60; i++) { const m = 0.5 * (lo + hi); if (isa(m).rho > rho) lo = m; else hi = m; }
  return 0.5 * (lo + hi);
}
/** Saturation vapour pressure over water [Pa] (Magnus/Tetens), T in K. */
export const pSat = (T) => 610.94 * Math.exp((17.625 * (T - 273.15)) / (T - 30.11));
/** Moist-air density [kg/m3] from p [Pa], T [K] and relative humidity (0-1). */
export function moistDensity(p, T, rh = 0) { const pv = rh * pSat(T); return (p - pv) / (R_AIR * T) + pv / (461.495 * T); }

export const tasFromMach = (M, h, dT = 0) => M * isa(h, dT).a;
export const machFromTas = (V, h, dT = 0) => V / isa(h, dT).a;
export const easFromTas = (V, h, dT = 0) => V * Math.sqrt(isa(h, dT).sigma);
/** Calibrated airspeed from Mach and altitude (subsonic compressible pitot relation). */
export function casFromMach(M, h, dT = 0) {
  const { p } = isa(h, dT), qc = p * ((1 + 0.2 * M * M) ** 3.5 - 1);
  return A0 * Math.sqrt(5 * ((qc / P0 + 1) ** (2 / 7) - 1));
}
export const dynPressure = (V, h, dT = 0) => 0.5 * isa(h, dT).rho * V * V;
export const reynolds = (V, L, h, dT = 0) => (V * L) / isa(h, dT).nu;

/** Isentropic ratios at Mach M: {T0_T, p0_p, rho0_rho}. */
export function isentropic(M, g = GAMMA) { const t = 1 + 0.5 * (g - 1) * M * M; return { T0_T: t, p0_p: t ** (g / (g - 1)), rho0_rho: t ** (1 / (g - 1)) }; }
/** Normal-shock jump relations for upstream Mach M1 > 1. */
export function normalShock(M1, g = GAMMA) {
  const m2 = M1 * M1, p2_p1 = 1 + ((2 * g) / (g + 1)) * (m2 - 1), rho2_rho1 = ((g + 1) * m2) / ((g - 1) * m2 + 2);
  const M2 = Math.sqrt((1 + 0.5 * (g - 1) * m2) / (g * m2 - 0.5 * (g - 1)));
  const T2_T1 = p2_p1 / rho2_rho1;
  const p02_p01 = (((g + 1) * m2) / ((g - 1) * m2 + 2)) ** (g / (g - 1)) * ((g + 1) / (2 * g * m2 - (g - 1))) ** (1 / (g - 1));
  return { M2, p2_p1, rho2_rho1, T2_T1, p02_p01 };
}
