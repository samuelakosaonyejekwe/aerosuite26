// The shared aircraft case: one data model that all 26 suites read from, so a quantity is entered
// once and propagates everywhere. Presets are representative class-level figures (not any specific
// certified type) meant as editable starting points.

const base = () => ({
  meta: { name: '', type: 'aeroplane', notes: '' },
  // site: where the aircraft operates; weather is filled by the live-data connectors
  site: { name: 'Sea-level standard day', lat: null, lon: null, elev_m: 0, runway_len_m: 2500, runway_slope_pct: 0, runway_mu: 0.04, runway_mu_brake: 0.4, T_C: 15, p_hPa: 1013.25, rh: 0.5, wind_ms: 0, wind_dir_deg: 0, gust_ms: 0, runway_heading_deg: null, runway_surface: 'asphalt', precip_mm_h: 0, visibility_m: 10000, cloud_pct: 0, freezing_level_m: 2300, kp_index: 2, wave_height_m: 0, design_hot_C: 30, winds_aloft: [], country: '', source: 'ISA default', updated: null },
  atm: { alt_m: 0, dISA_K: 0, turb_intensity: 0.01 },
  flight: { V_ms: 70, alpha_deg: 3, beta_deg: 0, n_load: 1 },
  mass: { mtow_kg: 0, oew_kg: 0, payload_kg: 0, fuel_kg: 0, cg_pct_mac: 25, Ixx: 0, Iyy: 0, Izz: 0, Ixz: 0 },
  wing: { S_m2: 0, b_m: 0, taper: 0.4, sweep_deg: 0, twist_deg: -2, dihedral_deg: 3, airfoil: '2412', tc: 0.12, incidence_deg: 2 },
  htail: { S_m2: 0, b_m: 0, arm_m: 0, tc: 0.1 },
  vtail: { S_m2: 0, b_m: 0, arm_m: 0 },
  fuselage: { len_m: 0, dia_m: 0, cabin_dp_Pa: 0 },
  aero: { CD0: 0.025, e: 0.8, CLmax_clean: 1.5, CLmax_to: 1.9, CLmax_land: 2.3, Cm0: -0.05, n_pos: 2.5, n_neg: -1, Vmo_ms: 0, Mmo: 0 },
  controls: { Se_Sh: 0.3, Sa_S: 0.06, Sr_Sv: 0.3, de_max_deg: 25, da_max_deg: 20, dr_max_deg: 25, rate_max_dps: 60 },
  prop: { type: 'turbofan', n_eng: 2, T0_N: 0, P0_W: 0, bpr: 5, opr: 30, tit_K: 1600, tsfc_kg_Ns: 1.6e-5, bsfc_kg_Ws: 8e-8, prop_dia_m: 0, n_blades: 0, rpm: 0, eta_prop: 0.82, fuel: 'Jet A-1' },
  rotor: { R_m: 0, n_blades: 0, chord_m: 0, twist_deg: -10, rpm: 0, cla: 5.73, cd0: 0.01, hinge_offset: 0.04, lock: 8, blade_mass_kg: 0, tr_R_m: 0, tr_chord_m: 0, tr_n_blades: 0, tr_rpm: 0, tr_arm_m: 0, shaft_tilt_deg: 3, flat_plate_m2: 0 },
  struct: { material: 'Al 2024-T3', box_chord_frac: 0.45, box_height_frac: 0.85, t_skin_mm: 3, t_spar_mm: 5, sf_ultimate: 1.5, zeta: 0.02, layup: '[0/45/-45/90]s', ply: 'IM7/8552 carbon-epoxy' },
  gear: { type: 'tricycle', n_main: 2, tyres_per_strut: 2, stroke_m: 0.4, sink_ms: 3.05, strut_eff: 0.8, tyre_dia_m: 1.1, tyre_k_Npm: 2e6, wheelbase_m: 12, track_m: 7.6 },
  systems: { hyd_p_Pa: 20.7e6, bus_V: 115, gen_kVA: 90, batt_kWh: 2, batt_chem: 'Li-ion NMC (cell)', motor_kW: 0, cabin_alt_m: 2400, pax_heat_W: 100, bleed: true },
  mission: { range_km: 0, cruise_alt_m: 0, cruise_mach: 0, cruise_V_ms: 0, reserve_min: 45, alternate_km: 185, pax: 0, hover_min: 0, profile: 'A-to-B' },
  econ: { price_usd: 0, fuel_usd_kg: 0.85, elec_usd_kWh: 0.14, util_fh_yr: 3000, cycles_yr: 1500, seats: 0, load_factor: 0.82, yield_usd_pkm: 0.09, crew_usd_fh: 600, maint_usd_fh: 0, life_yr: 25, residual_frac: 0.15, discount: 0.08, inflation: 0.025, tax: 0.21, carbon_usd_t: 80, currency: 'USD' },
});

const merge = (t, s) => { for (const k of Object.keys(s)) { if (s[k] && typeof s[k] === 'object' && !Array.isArray(s[k])) t[k] = merge(t[k] || {}, s[k]); else t[k] = s[k]; } return t; };
const make = (o) => merge(base(), o);

export const PRESETS = {
  narrowbody: {
    label: 'Narrow-body jet airliner (150-180 seat class)', icon: 'jet',
    data: () => make({
      meta: { name: 'Narrow-body airliner', type: 'aeroplane' },
      atm: { alt_m: 10668 }, flight: { V_ms: 231, alpha_deg: 2.5 },
      mass: { mtow_kg: 78000, oew_kg: 42600, payload_kg: 16600, fuel_kg: 18800, Ixx: 1.28e6, Iyy: 3.8e6, Izz: 4.9e6, Ixz: 1.0e5 },
      wing: { S_m2: 122.6, b_m: 34.1, taper: 0.24, sweep_deg: 25, twist_deg: -3, dihedral_deg: 5, airfoil: '2412', tc: 0.115 },
      htail: { S_m2: 31, b_m: 12.45, arm_m: 17.7 }, vtail: { S_m2: 21.5, b_m: 5.87, arm_m: 16.2 },
      fuselage: { len_m: 37.6, dia_m: 3.95, cabin_dp_Pa: 57000 },
      aero: { CD0: 0.0205, e: 0.8, CLmax_clean: 1.5, CLmax_to: 2.2, CLmax_land: 2.9, Cm0: -0.08, n_pos: 2.5, n_neg: -1, Vmo_ms: 180, Mmo: 0.82 },
      prop: { type: 'turbofan', n_eng: 2, T0_N: 120000, bpr: 5.7, opr: 32, tit_K: 1650, tsfc_kg_Ns: 1.62e-5 },
      gear: { stroke_m: 0.45, tyre_dia_m: 1.17, wheelbase_m: 12.6, track_m: 7.6, tyre_k_Npm: 2.2e6 },
      systems: { hyd_p_Pa: 20.7e6, bus_V: 115, gen_kVA: 90 },
      mission: { range_km: 4500, cruise_alt_m: 10668, cruise_mach: 0.78, cruise_V_ms: 231, pax: 165 },
      econ: { price_usd: 105e6, seats: 165, util_fh_yr: 3400, cycles_yr: 1650, crew_usd_fh: 950, maint_usd_fh: 1150 },
    }),
  },
  turboprop: {
    label: 'Regional twin turboprop (70 seat class)', icon: 'prop',
    data: () => make({
      meta: { name: 'Regional turboprop', type: 'aeroplane' },
      atm: { alt_m: 6100 }, flight: { V_ms: 140, alpha_deg: 2 },
      mass: { mtow_kg: 23000, oew_kg: 13300, payload_kg: 7400, fuel_kg: 5000, Ixx: 2.5e5, Iyy: 4.6e5, Izz: 6.6e5, Ixz: 2e4 },
      wing: { S_m2: 61, b_m: 27.05, taper: 0.59, sweep_deg: 3, twist_deg: -2, dihedral_deg: 2.5, airfoil: '4318', tc: 0.18 },
      htail: { S_m2: 11.7, b_m: 7.3, arm_m: 13.5 }, vtail: { S_m2: 12.5, b_m: 4.5, arm_m: 12.2 },
      fuselage: { len_m: 27.2, dia_m: 2.87, cabin_dp_Pa: 41000 },
      aero: { CD0: 0.027, e: 0.82, CLmax_clean: 1.6, CLmax_to: 2.1, CLmax_land: 2.7, n_pos: 2.5, n_neg: -1, Vmo_ms: 128, Mmo: 0.55 },
      prop: { type: 'turboprop', n_eng: 2, P0_W: 1.846e6, T0_N: 33000, opr: 15, tit_K: 1450, bsfc_kg_Ws: 7.9e-8, prop_dia_m: 3.93, n_blades: 6, rpm: 1200, eta_prop: 0.86 },
      gear: { stroke_m: 0.35, tyre_dia_m: 0.87, wheelbase_m: 10.8, track_m: 4.1, tyre_k_Npm: 1.4e6 },
      mission: { range_km: 1300, cruise_alt_m: 6100, cruise_mach: 0.44, cruise_V_ms: 140, pax: 70 },
      econ: { price_usd: 27e6, seats: 70, util_fh_yr: 2600, cycles_yr: 2800, crew_usd_fh: 520, maint_usd_fh: 620, yield_usd_pkm: 0.17 },
    }),
  },
  ga: {
    label: 'Light piston single (4 seat general aviation)', icon: 'ga',
    data: () => make({
      meta: { name: 'Light piston single', type: 'aeroplane' },
      atm: { alt_m: 2400 }, flight: { V_ms: 62, alpha_deg: 2 },
      mass: { mtow_kg: 1111, oew_kg: 767, payload_kg: 200, fuel_kg: 144, Ixx: 1285, Iyy: 1825, Izz: 2667, Ixz: 0 },
      wing: { S_m2: 16.2, b_m: 11, taper: 0.69, sweep_deg: 0, twist_deg: -3, dihedral_deg: 1.7, airfoil: '2412', tc: 0.12 },
      htail: { S_m2: 3.35, b_m: 3.45, arm_m: 4.6 }, vtail: { S_m2: 1.74, b_m: 1.5, arm_m: 4.7 },
      fuselage: { len_m: 8.28, dia_m: 1.2, cabin_dp_Pa: 0 },
      aero: { CD0: 0.031, e: 0.75, CLmax_clean: 1.52, CLmax_to: 1.7, CLmax_land: 2.1, n_pos: 3.8, n_neg: -1.52, Vmo_ms: 63, Mmo: 0.25 },
      prop: { type: 'piston', n_eng: 1, P0_W: 119000, T0_N: 2600, bsfc_kg_Ws: 7.6e-8, prop_dia_m: 1.9, n_blades: 2, rpm: 2700, eta_prop: 0.8, fuel: 'Avgas 100LL' },
      gear: { n_main: 2, tyres_per_strut: 1, stroke_m: 0.18, tyre_dia_m: 0.44, wheelbase_m: 1.65, track_m: 2.5, tyre_k_Npm: 2.5e5 },
      systems: { hyd_p_Pa: 6.9e6, bus_V: 28, gen_kVA: 1.7, batt_kWh: 0.6, bleed: false, cabin_alt_m: 0 },
      mission: { range_km: 1100, cruise_alt_m: 2400, cruise_mach: 0.19, cruise_V_ms: 62, pax: 3 },
      econ: { price_usd: 450000, seats: 4, util_fh_yr: 350, cycles_yr: 400, crew_usd_fh: 0, maint_usd_fh: 45, fuel_usd_kg: 2.2, yield_usd_pkm: 0, load_factor: 0.6 },
    }),
  },
  helicopter: {
    label: 'Medium twin-turbine utility helicopter', icon: 'heli',
    data: () => make({
      meta: { name: 'Medium utility helicopter', type: 'helicopter' },
      atm: { alt_m: 500 }, flight: { V_ms: 70, alpha_deg: -3 },
      mass: { mtow_kg: 9980, oew_kg: 5220, payload_kg: 3600, fuel_kg: 1100, Ixx: 6300, Iyy: 52000, Izz: 50000, Ixz: 2200 },
      wing: { S_m2: 0, b_m: 0 }, htail: { S_m2: 4.2, b_m: 4.4, arm_m: 9 }, vtail: { S_m2: 3, b_m: 2.5, arm_m: 9.9 },
      fuselage: { len_m: 15.3, dia_m: 2.4, cabin_dp_Pa: 0 },
      aero: { CD0: 0.02, e: 0.8, n_pos: 3.5, n_neg: -1, Vmo_ms: 98, Mmo: 0.3 },
      prop: { type: 'turboshaft', n_eng: 2, P0_W: 1.41e6, T0_N: 0, opr: 17, tit_K: 1500, bsfc_kg_Ws: 7.8e-8 },
      rotor: { R_m: 8.18, n_blades: 4, chord_m: 0.53, twist_deg: -16, rpm: 258, cla: 5.73, cd0: 0.0095, hinge_offset: 0.047, lock: 8.2, blade_mass_kg: 116, tr_R_m: 1.68, tr_chord_m: 0.247, tr_n_blades: 4, tr_rpm: 1190, tr_arm_m: 9.93, shaft_tilt_deg: 3, flat_plate_m2: 3.35 },
      gear: { type: 'tailwheel', stroke_m: 0.3, sink_ms: 2.44, tyre_dia_m: 0.66, wheelbase_m: 8.8, track_m: 2.7, tyre_k_Npm: 8e5 },
      systems: { hyd_p_Pa: 20.7e6, bus_V: 115, gen_kVA: 45, bleed: true, cabin_alt_m: 0 },
      mission: { range_km: 500, cruise_alt_m: 500, cruise_V_ms: 72, cruise_mach: 0.21, pax: 12, hover_min: 10, profile: 'Offshore transport' },
      econ: { price_usd: 21e6, seats: 12, util_fh_yr: 600, cycles_yr: 700, crew_usd_fh: 480, maint_usd_fh: 1900, yield_usd_pkm: 0, load_factor: 0.7 },
    }),
  },
  fixedUav: {
    label: 'Fixed-wing tactical UAV (25 kg class)', icon: 'uav',
    data: () => make({
      meta: { name: 'Fixed-wing UAV', type: 'uav' },
      atm: { alt_m: 1000 }, flight: { V_ms: 25, alpha_deg: 3 },
      mass: { mtow_kg: 25, oew_kg: 15, payload_kg: 5, fuel_kg: 5, Ixx: 3.5, Iyy: 4.2, Izz: 7.2, Ixz: 0.1 },
      wing: { S_m2: 1.1, b_m: 3.3, taper: 0.6, sweep_deg: 2, twist_deg: -2, dihedral_deg: 3, airfoil: '4412', tc: 0.12 },
      htail: { S_m2: 0.18, b_m: 0.85, arm_m: 1.25 }, vtail: { S_m2: 0.11, b_m: 0.4, arm_m: 1.25 },
      fuselage: { len_m: 2.2, dia_m: 0.25 },
      aero: { CD0: 0.03, e: 0.8, CLmax_clean: 1.4, CLmax_to: 1.4, CLmax_land: 1.6, n_pos: 3.8, n_neg: -1.5, Vmo_ms: 40, Mmo: 0.15 },
      prop: { type: 'piston', n_eng: 1, P0_W: 2800, T0_N: 90, bsfc_kg_Ws: 1.4e-7, prop_dia_m: 0.56, n_blades: 2, rpm: 6500, eta_prop: 0.7, fuel: 'Avgas 100LL' },
      gear: { n_main: 2, tyres_per_strut: 1, stroke_m: 0.05, sink_ms: 2, tyre_dia_m: 0.1, wheelbase_m: 0.7, track_m: 0.6, tyre_k_Npm: 3e4 },
      systems: { hyd_p_Pa: 0, bus_V: 28, gen_kVA: 0.3, batt_kWh: 0.3, bleed: false, cabin_alt_m: 0 },
      mission: { range_km: 800, cruise_alt_m: 1000, cruise_V_ms: 25, cruise_mach: 0.075, pax: 0, profile: 'Survey / loiter' },
      econ: { price_usd: 160000, seats: 0, util_fh_yr: 800, cycles_yr: 300, crew_usd_fh: 90, maint_usd_fh: 30, fuel_usd_kg: 2.2, yield_usd_pkm: 0 },
    }),
  },
  multirotor: {
    label: 'Electric quadrotor UAV (multirotor)', icon: 'quad',
    data: () => make({
      meta: { name: 'Electric quadrotor', type: 'uav' },
      atm: { alt_m: 100 }, flight: { V_ms: 12, alpha_deg: -8 },
      mass: { mtow_kg: 6.3, oew_kg: 3.4, payload_kg: 1, fuel_kg: 0, Ixx: 0.12, Iyy: 0.12, Izz: 0.21, Ixz: 0 },
      wing: { S_m2: 0, b_m: 0 }, fuselage: { len_m: 0.5, dia_m: 0.25 },
      aero: { CD0: 0.05, n_pos: 3, n_neg: 0, Vmo_ms: 23, Mmo: 0.07 },
      prop: { type: 'electric', n_eng: 4, P0_W: 500, T0_N: 32, prop_dia_m: 0.43, n_blades: 2, rpm: 5200, eta_prop: 0.7 },
      rotor: { R_m: 0.215, n_blades: 2, chord_m: 0.03, twist_deg: -12, rpm: 5200, cd0: 0.02, flat_plate_m2: 0.03, blade_mass_kg: 0.02, lock: 3 },
      gear: { type: 'skid', stroke_m: 0.03, sink_ms: 1.5, wheelbase_m: 0.4, track_m: 0.4, tyre_k_Npm: 2e4 },
      systems: { hyd_p_Pa: 0, bus_V: 44.4, gen_kVA: 0, batt_kWh: 0.26, motor_kW: 0.5, bleed: false, cabin_alt_m: 0 },
      mission: { range_km: 6, cruise_alt_m: 100, cruise_V_ms: 12, cruise_mach: 0.035, hover_min: 3, pax: 0, profile: 'Inspection' },
      econ: { price_usd: 14000, seats: 0, util_fh_yr: 400, cycles_yr: 1200, crew_usd_fh: 60, maint_usd_fh: 6, yield_usd_pkm: 0 },
    }),
  },
  evtol: {
    label: 'Battery-electric eVTOL air taxi (lift + cruise)', icon: 'evtol',
    data: () => make({
      meta: { name: 'eVTOL air taxi', type: 'evtol' },
      atm: { alt_m: 600 }, flight: { V_ms: 67, alpha_deg: 3 },
      mass: { mtow_kg: 2400, oew_kg: 1950, payload_kg: 450, fuel_kg: 0, Ixx: 4200, Iyy: 3900, Izz: 7600, Ixz: 100 },
      wing: { S_m2: 13.5, b_m: 11.5, taper: 0.6, sweep_deg: 0, twist_deg: -2, dihedral_deg: 2, airfoil: '4415', tc: 0.15 },
      htail: { S_m2: 2.6, b_m: 3.4, arm_m: 4.4 }, vtail: { S_m2: 1.9, b_m: 1.5, arm_m: 4.4 },
      fuselage: { len_m: 7.3, dia_m: 1.5 },
      aero: { CD0: 0.035, e: 0.8, CLmax_clean: 1.5, CLmax_to: 1.7, CLmax_land: 1.9, n_pos: 3, n_neg: -1, Vmo_ms: 85, Mmo: 0.26 },
      prop: { type: 'electric', n_eng: 6, P0_W: 150000, T0_N: 5200, prop_dia_m: 2.8, n_blades: 5, rpm: 1100, eta_prop: 0.83 },
      rotor: { R_m: 1.4, n_blades: 5, chord_m: 0.16, twist_deg: -18, rpm: 1100, cd0: 0.012, flat_plate_m2: 0.9, blade_mass_kg: 3.5, lock: 5 },
      struct: { material: 'Al 7075-T6' },
      gear: { n_main: 2, tyres_per_strut: 1, stroke_m: 0.15, sink_ms: 2.4, tyre_dia_m: 0.4, wheelbase_m: 3, track_m: 2.4, tyre_k_Npm: 4e5 },
      systems: { hyd_p_Pa: 0, bus_V: 800, gen_kVA: 0, batt_kWh: 160, motor_kW: 150, bleed: false, cabin_alt_m: 0 },
      mission: { range_km: 120, cruise_alt_m: 600, cruise_V_ms: 67, cruise_mach: 0.2, hover_min: 2, pax: 4, reserve_min: 20, profile: 'Urban air taxi' },
      econ: { price_usd: 3.5e6, seats: 4, util_fh_yr: 1800, cycles_yr: 6000, crew_usd_fh: 140, maint_usd_fh: 180, yield_usd_pkm: 1.9, load_factor: 0.65 },
    }),
  },
};

export const defaultCase = () => PRESETS.narrowbody.data();
export const blankCase = base;

/** Fill in any fields missing from an imported or older case with baseline values. */
export function normaliseCase(c) { return merge(base(), c || {}); }

/** Derived geometric and mass quantities used by many suites. */
export function derived(c) {
  const w = c.wing, S = w.S_m2, b = w.b_m, AR = S > 0 ? (b * b) / S : 0, lam = w.taper;
  const cr = S > 0 ? (2 * S) / (b * (1 + lam)) : 0, ct = cr * lam;
  const mac = cr ? ((2 / 3) * cr * (1 + lam + lam * lam)) / (1 + lam) : 0;
  const W = c.mass.mtow_kg * 9.80665;
  const r = c.rotor, A_disk = Math.PI * r.R_m * r.R_m, omega = (r.rpm * 2 * Math.PI) / 60;
  return {
    AR, c_root: cr, c_tip: ct, mac, W, wing_loading: S > 0 ? W / S : 0,
    k_induced: AR > 0 ? 1 / (Math.PI * AR * c.aero.e) : 0,
    y_mac: (b / 6) * ((1 + 2 * lam) / (1 + lam)),
    A_disk, omega, v_tip: omega * r.R_m, solidity: A_disk ? (r.n_blades * r.chord_m) / (Math.PI * r.R_m) : 0,
    disk_loading: A_disk ? W / (A_disk * (c.meta.type === 'helicopter' ? 1 : Math.max(1, c.prop.n_eng))) : 0,
    isRotary: c.meta.type === 'helicopter' || (c.wing.S_m2 === 0 && c.rotor.R_m > 0),
    T_total: c.prop.T0_N * c.prop.n_eng, P_total: c.prop.P0_W * c.prop.n_eng,
    zfw_kg: c.mass.oew_kg + c.mass.payload_kg,
  };
}

/** Field metadata for the case editor: [section, key, label, unit, help]. */
export const CASE_FIELDS = [
  ['meta', 'name', 'Case name', '', 'A short name for this study'],
  ['meta', 'type', 'Vehicle class', '', 'aeroplane, helicopter, uav or evtol', ['aeroplane', 'helicopter', 'uav', 'evtol']],
  ['atm', 'alt_m', 'Analysis altitude', 'm', 'Pressure altitude of the main analysis point'],
  ['atm', 'dISA_K', 'ISA temperature deviation', 'K', 'Outside air temperature minus the standard-day value'],
  ['atm', 'turb_intensity', 'Turbulence intensity', '-', 'RMS gust velocity divided by airspeed'],
  ['flight', 'V_ms', 'True airspeed', 'm/s', 'Speed of the main analysis point'],
  ['flight', 'alpha_deg', 'Angle of attack', 'deg', ''],
  ['flight', 'beta_deg', 'Sideslip angle', 'deg', ''],
  ['flight', 'n_load', 'Load factor', 'g', 'Lift divided by weight'],
  ['mass', 'mtow_kg', 'Maximum take-off mass', 'kg', ''],
  ['mass', 'oew_kg', 'Operating empty mass', 'kg', ''],
  ['mass', 'payload_kg', 'Payload', 'kg', ''],
  ['mass', 'fuel_kg', 'Fuel (or 0 for electric)', 'kg', ''],
  ['mass', 'cg_pct_mac', 'Centre of gravity', '% MAC', 'Longitudinal CG position as a percentage of mean aerodynamic chord'],
  ['mass', 'Ixx', 'Roll inertia Ixx', 'kg m²', ''], ['mass', 'Iyy', 'Pitch inertia Iyy', 'kg m²', ''], ['mass', 'Izz', 'Yaw inertia Izz', 'kg m²', ''], ['mass', 'Ixz', 'Product of inertia Ixz', 'kg m²', ''],
  ['wing', 'S_m2', 'Wing reference area', 'm²', 'Set 0 for pure rotorcraft'], ['wing', 'b_m', 'Wing span', 'm', ''], ['wing', 'taper', 'Taper ratio', '-', 'Tip chord / root chord'],
  ['wing', 'sweep_deg', 'Quarter-chord sweep', 'deg', ''], ['wing', 'twist_deg', 'Tip twist (washout negative)', 'deg', ''], ['wing', 'dihedral_deg', 'Dihedral', 'deg', ''],
  ['wing', 'airfoil', 'NACA 4-digit section', '', 'e.g. 2412, 0012, 4415'], ['wing', 'tc', 'Thickness ratio', '-', ''], ['wing', 'incidence_deg', 'Wing incidence', 'deg', ''],
  ['htail', 'S_m2', 'Horizontal tail area', 'm²', ''], ['htail', 'b_m', 'Horizontal tail span', 'm', ''], ['htail', 'arm_m', 'Horizontal tail arm', 'm', 'CG to tail aerodynamic centre'],
  ['vtail', 'S_m2', 'Vertical tail area', 'm²', ''], ['vtail', 'b_m', 'Vertical tail height', 'm', ''], ['vtail', 'arm_m', 'Vertical tail arm', 'm', ''],
  ['fuselage', 'len_m', 'Fuselage length', 'm', ''], ['fuselage', 'dia_m', 'Fuselage diameter', 'm', ''], ['fuselage', 'cabin_dp_Pa', 'Cabin pressure differential', 'Pa', ''],
  ['aero', 'CD0', 'Zero-lift drag coefficient', '-', ''], ['aero', 'e', 'Oswald efficiency', '-', ''], ['aero', 'CLmax_clean', 'CLmax clean', '-', ''], ['aero', 'CLmax_to', 'CLmax take-off', '-', ''], ['aero', 'CLmax_land', 'CLmax landing', '-', ''],
  ['aero', 'Cm0', 'Zero-lift pitching moment', '-', ''], ['aero', 'n_pos', 'Positive limit load factor', 'g', ''], ['aero', 'n_neg', 'Negative limit load factor', 'g', ''], ['aero', 'Vmo_ms', 'Max operating speed (EAS)', 'm/s', ''], ['aero', 'Mmo', 'Max operating Mach', '-', ''],
  ['prop', 'type', 'Powerplant type', '', '', ['turbofan', 'turbojet', 'turboprop', 'turboshaft', 'piston', 'electric']],
  ['prop', 'n_eng', 'Number of engines / motors', '', ''], ['prop', 'T0_N', 'Static thrust per engine', 'N', ''], ['prop', 'P0_W', 'Rated shaft power per engine', 'W', ''],
  ['prop', 'bpr', 'Bypass ratio', '-', ''], ['prop', 'opr', 'Overall pressure ratio', '-', ''], ['prop', 'tit_K', 'Turbine inlet temperature', 'K', ''],
  ['prop', 'tsfc_kg_Ns', 'Thrust-specific fuel consumption', 'kg/N/s', ''], ['prop', 'bsfc_kg_Ws', 'Power-specific fuel consumption', 'kg/W/s', ''],
  ['prop', 'prop_dia_m', 'Propeller diameter', 'm', ''], ['prop', 'n_blades', 'Propeller blades', '', ''], ['prop', 'rpm', 'Propeller speed', 'rpm', ''], ['prop', 'eta_prop', 'Propeller efficiency', '-', ''],
  ['prop', 'fuel', 'Fuel', '', '', ['Jet A-1', 'Avgas 100LL', 'SAF (HEFA-SPK)', 'Liquid hydrogen']],
  ['rotor', 'R_m', 'Main rotor radius', 'm', ''], ['rotor', 'n_blades', 'Main rotor blades', '', ''], ['rotor', 'chord_m', 'Blade chord', 'm', ''], ['rotor', 'twist_deg', 'Blade linear twist', 'deg', ''], ['rotor', 'rpm', 'Rotor speed', 'rpm', ''],
  ['rotor', 'cla', 'Blade lift-curve slope', '1/rad', ''], ['rotor', 'cd0', 'Blade profile drag coefficient', '-', ''], ['rotor', 'hinge_offset', 'Flap hinge offset e/R', '-', ''], ['rotor', 'lock', 'Lock number', '-', ''], ['rotor', 'blade_mass_kg', 'Blade mass', 'kg', ''],
  ['rotor', 'tr_R_m', 'Tail rotor radius', 'm', ''], ['rotor', 'tr_chord_m', 'Tail rotor chord', 'm', ''], ['rotor', 'tr_n_blades', 'Tail rotor blades', '', ''], ['rotor', 'tr_rpm', 'Tail rotor speed', 'rpm', ''], ['rotor', 'tr_arm_m', 'Tail rotor arm', 'm', ''],
  ['rotor', 'shaft_tilt_deg', 'Shaft forward tilt', 'deg', ''], ['rotor', 'flat_plate_m2', 'Equivalent flat-plate area', 'm²', 'Fuselage parasite drag area'],
  ['struct', 'material', 'Primary structural metal', '', '', 'METALS'], ['struct', 'ply', 'Composite ply system', '', '', 'PLIES'], ['struct', 'layup', 'Laminate stacking sequence', '', 'e.g. [0/45/-45/90]s'],
  ['struct', 'box_chord_frac', 'Wing-box chord fraction', '-', ''], ['struct', 'box_height_frac', 'Wing-box height / max thickness', '-', ''], ['struct', 't_skin_mm', 'Skin thickness', 'mm', ''], ['struct', 't_spar_mm', 'Spar web thickness', 'mm', ''],
  ['struct', 'sf_ultimate', 'Ultimate safety factor', '-', ''], ['struct', 'zeta', 'Structural damping ratio', '-', ''],
  ['gear', 'type', 'Gear layout', '', '', ['tricycle', 'tailwheel', 'skid']], ['gear', 'n_main', 'Main gear struts', '', ''], ['gear', 'tyres_per_strut', 'Tyres per strut', '', ''], ['gear', 'stroke_m', 'Shock strut stroke', 'm', ''],
  ['gear', 'sink_ms', 'Design sink rate', 'm/s', ''], ['gear', 'strut_eff', 'Strut efficiency', '-', ''], ['gear', 'tyre_dia_m', 'Tyre diameter', 'm', ''], ['gear', 'tyre_k_Npm', 'Tyre vertical stiffness', 'N/m', ''], ['gear', 'wheelbase_m', 'Wheelbase', 'm', ''], ['gear', 'track_m', 'Track', 'm', ''],
  ['systems', 'hyd_p_Pa', 'Hydraulic system pressure', 'Pa', ''], ['systems', 'bus_V', 'Main bus voltage', 'V', ''], ['systems', 'gen_kVA', 'Generator rating (each)', 'kVA', ''], ['systems', 'batt_kWh', 'Battery energy', 'kWh', ''],
  ['systems', 'batt_chem', 'Battery chemistry', '', '', 'BATTERIES'], ['systems', 'motor_kW', 'Propulsion motor power (each)', 'kW', ''], ['systems', 'cabin_alt_m', 'Max cabin altitude', 'm', ''], ['systems', 'pax_heat_W', 'Heat per occupant', 'W', ''],
  ['mission', 'range_km', 'Design range', 'km', ''], ['mission', 'cruise_alt_m', 'Cruise altitude', 'm', ''], ['mission', 'cruise_mach', 'Cruise Mach', '-', ''], ['mission', 'cruise_V_ms', 'Cruise true airspeed', 'm/s', ''],
  ['mission', 'reserve_min', 'Final reserve', 'min', ''], ['mission', 'alternate_km', 'Alternate distance', 'km', ''], ['mission', 'pax', 'Passengers', '', ''], ['mission', 'hover_min', 'Hover time', 'min', ''],
  ['econ', 'price_usd', 'Aircraft price', 'USD', ''], ['econ', 'fuel_usd_kg', 'Fuel price', 'USD/kg', 'Updated by the live fuel-price connector when available'], ['econ', 'elec_usd_kWh', 'Electricity price', 'USD/kWh', ''],
  ['econ', 'util_fh_yr', 'Utilisation', 'FH/yr', ''], ['econ', 'cycles_yr', 'Flight cycles', '1/yr', ''], ['econ', 'seats', 'Seats', '', ''], ['econ', 'load_factor', 'Load factor', '-', ''], ['econ', 'yield_usd_pkm', 'Passenger yield', 'USD/pax-km', ''],
  ['econ', 'crew_usd_fh', 'Crew cost', 'USD/FH', ''], ['econ', 'maint_usd_fh', 'Maintenance cost', 'USD/FH', ''], ['econ', 'life_yr', 'Service life', 'yr', ''], ['econ', 'residual_frac', 'Residual value fraction', '-', ''],
  ['econ', 'discount', 'Discount rate', '-', ''], ['econ', 'inflation', 'Inflation rate', '-', 'Updated by the live World Bank connector when available'], ['econ', 'tax', 'Tax rate', '-', ''], ['econ', 'carbon_usd_t', 'Carbon price', 'USD/t CO₂', ''],
];
export const CASE_SECTIONS = { meta: 'Identity', atm: 'Atmosphere', flight: 'Flight condition', mass: 'Mass & inertia', wing: 'Wing', htail: 'Horizontal tail', vtail: 'Vertical tail', fuselage: 'Fuselage', aero: 'Aerodynamic data', prop: 'Propulsion', rotor: 'Rotor system', struct: 'Structure & materials', gear: 'Landing gear', systems: 'Systems', mission: 'Mission', econ: 'Economics' };
