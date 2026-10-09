// Shared material database. Values are typical room-temperature handbook-class figures for
// preliminary analysis; they are NOT certified design allowables. Replace with programme-specific,
// traceable allowables (e.g. MMPDS / CMH-17 basis values) before any substantiation work.
// Units: E,G,Sy,Su [Pa]; rho [kg/m3]; k [W/m/K]; cp [J/kg/K]; alpha [1/K]; KIc [Pa sqrt(m)].
// Fatigue: Basquin S = sf*(2N)^b (stress amplitude); Coffin-Manson ef, c; Paris da/dN = C (dK)^m with dK in MPa sqrt(m), da/dN in m/cycle.
//
// Sy, Su are TYPICAL strengths. Where a statistically based design allowable was verified against MIL-HDBK-5J
// (see js/data/sources.json, ids mmpds-*), it is held separately as Sy_A (Fty) and Su_A (Ftu) [Pa], longitudinal grain,
// with allow_ref naming the basis, product form and table. Strength margins should use designAllowables(); fatigue,
// fracture and energy-absorption models keep the typical values. Elastic constants of the low-alloy steels and of
// Ti-6Al-4V are the MIL-HDBK-5J values; the 300M strengths equal its S-basis values for the 280 ksi condition. Fatigue,
// fracture, thermal and Johnson-Cook constants are unsourced typical values.

export const METALS = {
  'Al 2024-T3': { E: 73.1e9, G: 28e9, nu: 0.33, rho: 2780, Sy: 345e6, Su: 483e6, Sy_A: 324e6, Su_A: 441e6, allow_ref: 'MIL-HDBK-5J A-basis, bare sheet 0.010–0.128 in, Table 3.2.3.0(b1)', k: 121, cp: 875, alpha: 23.2e-6, KIc: 37e6, sf: 850e6, b: -0.086, ef: 0.22, c: -0.59, parisC: 5e-11, parisM: 3.0, dKth: 3.0, JC: { A: 369e6, B: 684e6, n: 0.73, C: 0.0083, m: 1.7, Tm: 775 } },
  'Al 7075-T6': { E: 71.7e9, G: 26.9e9, nu: 0.33, rho: 2810, Sy: 503e6, Su: 572e6, Sy_A: 483e6, Su_A: 538e6, allow_ref: 'MIL-HDBK-5J A-basis, bare sheet 0.040–0.125 in, Table 3.7.6.0(b1)', k: 130, cp: 960, alpha: 23.6e-6, KIc: 29e6, sf: 1050e6, b: -0.09, ef: 0.19, c: -0.52, parisC: 2.7e-11, parisM: 3.3, dKth: 2.5, JC: { A: 520e6, B: 477e6, n: 0.52, C: 0.001, m: 1.61, Tm: 893 } },
  'Al 7050-T7451': { E: 71.7e9, G: 26.9e9, nu: 0.33, rho: 2830, Sy: 469e6, Su: 524e6, Sy_A: 441e6, Su_A: 510e6, allow_ref: 'MIL-HDBK-5J S/A-basis design values, plate 0.250–1.500 in, Table 3.7.4.0(b1)', k: 157, cp: 860, alpha: 23.5e-6, KIc: 35e6, sf: 950e6, b: -0.09, ef: 0.2, c: -0.55, parisC: 3e-11, parisM: 3.1, dKth: 2.5 },
  'Al-Li 2195': { E: 76e9, G: 28.5e9, nu: 0.33, rho: 2700, Sy: 560e6, Su: 600e6, k: 90, cp: 900, alpha: 22e-6, KIc: 33e6, sf: 1000e6, b: -0.09, ef: 0.15, c: -0.55, parisC: 3e-11, parisM: 3.2, dKth: 2.5 },
  'Ti-6Al-4V': { E: 110.3e9, G: 42.7e9, nu: 0.31, rho: 4430, Sy: 880e6, Su: 950e6, Sy_A: 827e6, Su_A: 896e6, allow_ref: 'MIL-HDBK-5J A-basis, annealed plate 0.1875–2.000 in, Table 5.4.1.0(b)', k: 6.7, cp: 526, alpha: 8.6e-6, KIc: 75e6, sf: 1700e6, b: -0.095, ef: 0.8, c: -0.7, parisC: 1e-11, parisM: 3.2, dKth: 4.0, JC: { A: 1098e6, B: 1092e6, n: 0.93, C: 0.014, m: 1.1, Tm: 1878 } },
  'Steel 4340 (QT)': { E: 200e9, G: 75.8e9, nu: 0.32, rho: 7850, Sy: 1470e6, Su: 1720e6, k: 44.5, cp: 475, alpha: 12.3e-6, KIc: 60e6, sf: 2000e6, b: -0.091, ef: 0.48, c: -0.6, parisC: 5e-12, parisM: 3.0, dKth: 5.0, JC: { A: 792e6, B: 510e6, n: 0.26, C: 0.014, m: 1.03, Tm: 1793 } },
  'Steel 300M': { E: 200e9, G: 75.8e9, nu: 0.32, rho: 7870, Sy: 1586e6, Su: 1931e6, k: 37, cp: 448, alpha: 11.3e-6, KIc: 57e6, sf: 2100e6, b: -0.09, ef: 0.4, c: -0.6, parisC: 5e-12, parisM: 3.0, dKth: 5.0 },
  'Inconel 718': { E: 200e9, G: 77e9, nu: 0.29, rho: 8190, Sy: 1100e6, Su: 1375e6, Sy_A: 1000e6, Su_A: 1241e6, allow_ref: 'MIL-HDBK-5J A-basis, solution-treated and aged sheet 0.010–0.187 in, Table 6.3.5.0(b)', k: 11.4, cp: 435, alpha: 13e-6, KIc: 100e6, sf: 2200e6, b: -0.1, ef: 0.5, c: -0.65, parisC: 4e-12, parisM: 3.1, dKth: 6.0 },
  'Mg AZ31B': { E: 45e9, G: 17e9, nu: 0.35, rho: 1770, Sy: 200e6, Su: 260e6, k: 96, cp: 1000, alpha: 26e-6, KIc: 28e6, sf: 450e6, b: -0.1, ef: 0.15, c: -0.55, parisC: 2e-10, parisM: 3.2, dKth: 1.5 },
};

// Unidirectional ply properties: E1,E2,G12 [Pa], nu12, strengths Xt,Xc,Yt,Yc,S [Pa], ply thickness t [m],
// fracture energies GIc/GIIc [J/m2], CTE a1,a2 [1/K]. Typical published characterisation sets, not basis values; the
// IM7/8552 strengths are lower than the manufacturer's typical datasheet values and the AS4/3501-6 set is unverified.
export const PLIES = {
  'T300/5208 carbon-epoxy': { E1: 181e9, E2: 10.3e9, G12: 7.17e9, nu12: 0.28, Xt: 1500e6, Xc: 1500e6, Yt: 40e6, Yc: 246e6, S: 68e6, rho: 1600, t: 0.125e-3, GIc: 200, GIIc: 600, a1: 0.02e-6, a2: 22.5e-6 },
  'AS4/3501-6 carbon-epoxy': { E1: 142e9, E2: 10.3e9, G12: 7.2e9, nu12: 0.27, Xt: 2280e6, Xc: 1440e6, Yt: 57e6, Yc: 228e6, S: 71e6, rho: 1580, t: 0.125e-3, GIc: 190, GIIc: 570, a1: -0.9e-6, a2: 27e-6 },
  'IM7/8552 carbon-epoxy': { E1: 171e9, E2: 9.08e9, G12: 5.29e9, nu12: 0.32, Xt: 2326e6, Xc: 1200e6, Yt: 62e6, Yc: 200e6, S: 92e6, rho: 1570, t: 0.131e-3, GIc: 277, GIIc: 788, a1: -0.1e-6, a2: 31e-6 },
  'E-glass/epoxy': { E1: 38.6e9, E2: 8.27e9, G12: 4.14e9, nu12: 0.26, Xt: 1062e6, Xc: 610e6, Yt: 31e6, Yc: 118e6, S: 72e6, rho: 1800, t: 0.15e-3, GIc: 300, GIIc: 1200, a1: 8.6e-6, a2: 22.1e-6 },
  'Kevlar 49/epoxy': { E1: 76e9, E2: 5.5e9, G12: 2.3e9, nu12: 0.34, Xt: 1400e6, Xc: 235e6, Yt: 12e6, Yc: 53e6, S: 34e6, rho: 1460, t: 0.125e-3, GIc: 250, GIIc: 900, a1: -4e-6, a2: 79e-6 },
};

export const FLUIDS = {
  'Jet A-1': { rho: 804, mu: 1.6e-3, cp: 2010, k: 0.115, LHV: 43.15e6, co2_per_kg: 3.16, flash_C: 38 },
  'Avgas 100LL': { rho: 720, mu: 0.5e-3, cp: 2200, k: 0.11, LHV: 43.5e6, co2_per_kg: 3.1 },
  'SAF (HEFA-SPK)': { rho: 760, mu: 1.5e-3, cp: 2050, k: 0.115, LHV: 44.1e6, co2_per_kg: 3.16, lifecycle_factor: 0.25 },
  'Liquid hydrogen': { rho: 70.8, mu: 1.3e-5, cp: 9700, k: 0.1, LHV: 119.96e6, co2_per_kg: 0 },
  'MIL-PRF-5606 hydraulic': { rho: 860, mu: 0.0132, cp: 1900, k: 0.13, bulk: 1.4e9 },
  'Skydrol LD-4': { rho: 1000, mu: 0.011, cp: 1750, k: 0.13, bulk: 1.5e9 },
  'MIL-PRF-23699 oil': { rho: 1000, mu: 0.025, cp: 2000, k: 0.15 },
};

export const BATTERIES = {
  'Li-ion NMC (cell)': { wh_kg: 250, v_nom: 3.65, v_max: 4.2, v_min: 2.8, r_mohm_ah: 60, cp: 1000, cycles_80: 1200 },
  'Li-ion NCA (cell)': { wh_kg: 265, v_nom: 3.6, v_max: 4.2, v_min: 2.7, r_mohm_ah: 55, cp: 1000, cycles_80: 1000 },
  'LiFePO4 (cell)': { wh_kg: 160, v_nom: 3.2, v_max: 3.65, v_min: 2.5, r_mohm_ah: 80, cp: 1100, cycles_80: 3000 },
  'Li-metal / solid-state (projected)': { wh_kg: 400, v_nom: 3.8, v_max: 4.3, v_min: 3.0, r_mohm_ah: 70, cp: 1000, cycles_80: 800 },
};

export const metalNames = () => Object.keys(METALS);
/**
 * Static strengths to use in a margin of safety: the verified design allowable when the database holds one, otherwise the
 * typical value. Returns { Sy, Su, design (bool), basis (text naming which was used) }.
 */
export function designAllowables(m) {
  const design = m.Sy_A > 0 && m.Su_A > 0, Sy = design ? m.Sy_A : m.Sy, Su = design ? m.Su_A : m.Su, v = `Fty ${(Sy / 1e6).toFixed(0)} MPa, Ftu ${(Su / 1e6).toFixed(0)} MPa`;
  return { Sy, Su, design, basis: design ? `${m.allow_ref}: ${v}` : `typical handbook strengths (${v}); no verified design allowable is held for this material` };
}
export const plyNames = () => Object.keys(PLIES);
