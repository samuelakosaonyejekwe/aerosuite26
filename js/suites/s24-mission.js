// Suite 24 — Aircraft Mission and Operational Simulation.
// Point-mass (3-DOF, quasi-steady) mission integration by phase with fuel accounting and reserves, great-circle
// routing with winds and terrain, cruise and payload–range trades, energy-state climb optimisation, helicopter
// out-and-back and search missions, electric / eVTOL mission energy and stochastic fuel planning.

import * as N from '../core/numerics.js';
import { isa, G0, RHO0 } from '../core/atmosphere.js';
import { FLUIDS, BATTERIES } from '../data/materials.js';

// ---- shared physics -------------------------------------------------------------------------
const isProp = (t) => t === 'turboprop' || t === 'piston' || t === 'electric' || t === 'turboshaft';
const powerAvail = (p, h, dT = 0) => { const s = isa(h, dT).sigma; return p.ptype === 'piston' ? p.P0_W * Math.max(0, (s - 0.1325) / 0.8675) : p.ptype === 'electric' ? p.P0_W : p.P0_W * s ** 0.7; };
/** Total thrust available [N] (jet thrust lapse after Mattingly; propeller power lapse capped by static thrust). */
function thrustAvail(p, V, h, dT = 0) {
  const a = isa(h, dT), M = V / a.a;
  if (!isProp(p.ptype)) { const bpr = p.ptype === 'turbojet' ? 0 : p.bpr, d0 = a.delta * (1 + 0.2 * M * M) ** 3.5, r = Math.sqrt(Math.max(M, 0)); return p.n_eng * p.T0_N * Math.max(0.05, bpr >= 1.5 ? d0 * (1 - 0.49 * r) : bpr > 0 ? d0 * (1 - 0.3 * r) : d0 * (1 - 0.16 * r)); }
  return Math.min(p.n_eng * (p.T0_N || 0) * a.sigma ** 0.7 || Infinity, (p.eta_prop * p.n_eng * powerAvail(p, h, dT)) / Math.max(V, 1));
}
/** Fuel flow [kg/s] for thrust T at speed V. */
const fuelFlow = (p, T, V) => (p.ptype === 'electric' ? 0 : isProp(p.ptype) ? (p.bsfc * T * V) / p.eta_prop : p.tsfc * T);
const maxFuelFlow = (p) => (isProp(p.ptype) ? p.bsfc * p.P0_W * p.n_eng : p.tsfc * p.T0_N * p.n_eng);
const dragOf = (p, m, V, rho) => { const q = 0.5 * rho * V * V; return q * p.S * p.CD0 + (p.k * (m * G0) ** 2) / (q * p.S); };
/** Level-flight fuel flow [kg/s] at mass m. */
const ffLevel = (p, m, h, V, dT) => fuelFlow(p, dragOf(p, m, V, isa(h, dT).rho), V);
const vMinDrag = (p, m, rho) => Math.sqrt((2 * m * G0) / (rho * p.S)) * (p.k / p.CD0) ** 0.25;
/** Cruise at constant altitude and true airspeed over a ground distance [m]; midpoint integration in distance. Returns fuel and time. */
function cruiseSeg(p, m0, dist, h, V, headwind, dT, n) {
  const gs = Math.max(1, V - headwind), dx = dist / n; let m = m0;
  for (let j = 0; j < n; j++) { const k1 = ffLevel(p, m, h, V, dT) / gs, k2 = ffLevel(p, m - 0.5 * k1 * dx, h, V, dT) / gs; m -= k2 * dx; }
  return { fuel: m0 - m, time: dist / gs };
}
/** Rotor power required in level flight [W] for all rotors: Glauert induced, profile with advance-ratio growth, parasite. */
function rotorPower(r, m, V, rho) {
  const T = (m * G0) / r.nR, A = Math.PI * r.R * r.R, vh = Math.sqrt(T / (2 * rho * A)); let vi = vh;
  for (let k = 0; k < 40; k++) { const vn = (vh * vh) / Math.sqrt(V * V + vi * vi); if (Math.abs(vn - vi) < 1e-10) { vi = vn; break; } vi = 0.5 * (vi + vn); }
  const mu = V / r.vt, Pi = r.nR * r.kappa * T * vi, Po = ((r.nR * rho * A * r.vt ** 3 * r.sig * r.cd0) / 8) * (1 + 4.65 * mu * mu), Pp = 0.5 * rho * V ** 3 * r.f;
  return { P: (Pi + Po + Pp) / r.eta, Pi, Po, Pp, vi, mu };
}
/** Speeds for minimum power (best endurance) and maximum V/P (best range, with headwind). */
function rotorSpeeds(r, m, rho, wind = 0, Vmax = 110) {
  const Vs = N.linspace(1, Vmax, 110), P = Vs.map((V) => rotorPower(r, m, V, rho).P), jb = N.argmin(P), jr = N.argmax(Vs.map((V, j) => Math.max(0, V - wind) / P[j]));
  return { Vs, P, Vbe: Vs[jb], Pbe: P[jb], Vbr: Vs[jr], Pbr: P[jr] };
}
// ---- geodesy --------------------------------------------------------------------------------
const R_E = 6371008.8, WGS = { a: 6378137, f: 1 / 298.257223563 };
function haversine(p, q) { const a = N.rad(p[0]), b = N.rad(q[0]), dl = N.rad(q[1] - p[1]), s = Math.sin((b - a) / 2) ** 2 + Math.cos(a) * Math.cos(b) * Math.sin(dl / 2) ** 2; return 2 * R_E * Math.asin(Math.min(1, Math.sqrt(s))); }
const bearing = (p, q) => { const a = N.rad(p[0]), b = N.rad(q[0]), dl = N.rad(q[1] - p[1]); return (N.deg(Math.atan2(Math.sin(dl) * Math.cos(b), Math.cos(a) * Math.sin(b) - Math.sin(a) * Math.cos(b) * Math.cos(dl))) + 360) % 360; };
/** Destination on the sphere from a start point, initial bearing [deg] and distance [m]. */
function direct(p, brg, d) { const a = N.rad(p[0]), t = N.rad(brg), dr = d / R_E, lat = Math.asin(Math.sin(a) * Math.cos(dr) + Math.cos(a) * Math.sin(dr) * Math.cos(t)); return [N.deg(lat), ((p[1] + N.deg(Math.atan2(Math.sin(t) * Math.sin(dr) * Math.cos(a), Math.cos(dr) - Math.sin(a) * Math.sin(lat))) + 540) % 360) - 180]; }
/** Vincenty inverse geodesic on WGS-84 [m]; falls back to the haversine value for near-antipodal points. */
function vincenty(p, q) {
  const { a, f } = WGS, b = a * (1 - f), U1 = Math.atan((1 - f) * Math.tan(N.rad(p[0]))), U2 = Math.atan((1 - f) * Math.tan(N.rad(q[0]))), L = N.rad(q[1] - p[1]), sU1 = Math.sin(U1), cU1 = Math.cos(U1), sU2 = Math.sin(U2), cU2 = Math.cos(U2);
  let lam = L, sS = 0, cS = 1, sig = 0, cSqA = 1, c2m = 0;
  for (let it = 0; it < 200; it++) {
    const sl = Math.sin(lam), cl = Math.cos(lam); sS = Math.hypot(cU2 * sl, cU1 * sU2 - sU1 * cU2 * cl); if (sS === 0) return 0;
    cS = sU1 * sU2 + cU1 * cU2 * cl; sig = Math.atan2(sS, cS); const sA = (cU1 * cU2 * sl) / sS; cSqA = 1 - sA * sA; c2m = cSqA ? cS - (2 * sU1 * sU2) / cSqA : 0;
    const C = (f / 16) * cSqA * (4 + f * (4 - 3 * cSqA)), ln = L + (1 - C) * f * sA * (sig + C * sS * (c2m + C * cS * (-1 + 2 * c2m * c2m)));
    if (Math.abs(ln - lam) < 1e-12) { lam = ln; break; } lam = ln; if (it === 199) return haversine(p, q);
  }
  const u2 = (cSqA * (a * a - b * b)) / (b * b), A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2))), B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
  return b * A * (sig - B * sS * (c2m + (B / 4) * (cS * (-1 + 2 * c2m * c2m) - (B / 6) * c2m * (-3 + 4 * sS * sS) * (-3 + 4 * c2m * c2m))));
}
/** Wind triangle: ground speed and wind-correction angle for a course [deg], TAS and wind FROM direction [deg] at speed ws. */
function windTriangle(course, V, wdir, ws) { const rel = N.rad(wdir - course), xw = ws * Math.sin(rel), hw = ws * Math.cos(rel), s = N.clamp(xw / Math.max(V, 1e-9), -1, 1); return { gs: Math.max(0, V * Math.sqrt(1 - s * s) - hw), wca: N.deg(Math.asin(s)), headwind: hw, crosswind: xw }; }
const parsePts = (text) => { const out = []; let bad = 0; for (const t of String(text || '').split(/[;\n]+/)) { if (!t.trim()) continue; const v = t.split(/[,\s]+/).filter(Boolean).map(Number); if (v.length >= 2 && Math.abs(v[0]) <= 90 && Math.abs(v[1]) <= 180) out.push([v[0], v[1]]); else bad++; } return { pts: out, bad }; };
const parseNums = (text) => { const v = String(text || '').split(/[\s,;]+/).filter(Boolean).map(Number); return { vals: v.filter(Number.isFinite), bad: v.filter((x) => !Number.isFinite(x)).length }; };

// ---- shared inputs --------------------------------------------------------------------------
const AC = [
  { key: 'mtow', label: 'Maximum take-off mass', unit: 'kg', default: 78000, min: 0.1, group: 'Mass' },
  { key: 'oew', label: 'Operating empty mass', unit: 'kg', default: 42600, min: 0.05, group: 'Mass' },
  { key: 'payload', label: 'Payload', unit: 'kg', default: 16600, min: 0, group: 'Mass' },
  { key: 'fuel_load', label: 'Fuel loaded', unit: 'kg', default: 18800, min: 0, group: 'Mass' },
  { key: 'S', label: 'Wing area', unit: 'm²', default: 122.6, min: 0.01, group: 'Aerodynamics' },
  { key: 'CD0', label: 'Zero-lift drag coefficient', unit: '-', default: 0.0205, min: 0.001, group: 'Aerodynamics', help: 'From Suite 1 when available' },
  { key: 'k', label: 'Induced drag factor k', unit: '-', default: 0.042, min: 0.001, group: 'Aerodynamics' },
  { key: 'CLmax', label: 'Maximum lift coefficient (clean)', unit: '-', default: 1.5, min: 0.2, group: 'Aerodynamics' },
  { key: 'ptype', label: 'Powerplant type', type: 'select', options: ['turbofan', 'turbojet', 'turboprop', 'piston'], default: 'turbofan', group: 'Propulsion' },
  { key: 'n_eng', label: 'Engines', unit: '', default: 2, min: 1, step: 1, discrete: true, group: 'Propulsion' },
  { key: 'T0_N', label: 'Static thrust per engine', unit: 'N', default: 120000, min: 0, group: 'Propulsion' },
  { key: 'P0_W', label: 'Rated power per engine', unit: 'W', default: 0, min: 0, group: 'Propulsion' },
  { key: 'bpr', label: 'Bypass ratio', unit: '-', default: 5.7, min: 0, group: 'Propulsion' },
  { key: 'eta_prop', label: 'Propeller efficiency', unit: '-', default: 0.82, min: 0.1, max: 0.95, group: 'Propulsion' },
  { key: 'tsfc', label: 'TSFC', unit: 'kg/N/s', default: 1.62e-5, min: 0, group: 'Propulsion', help: 'From Suite 7 when available' },
  { key: 'bsfc', label: 'BSFC', unit: 'kg/W/s', default: 8e-8, min: 0, group: 'Propulsion' },
  { key: 'fuel', label: 'Fuel', type: 'select', options: ['Jet A-1', 'Avgas 100LL', 'SAF (HEFA-SPK)', 'Liquid hydrogen'], default: 'Jet A-1', group: 'Propulsion' },
];
const FLT = [
  { key: 'range_km', label: 'Stage length', unit: 'km', default: 4500, min: 1, group: 'Mission' },
  { key: 'cruise_alt', label: 'Cruise altitude', unit: 'm', default: 10668, min: 0, max: 16000, group: 'Mission' },
  { key: 'cruise_V', label: 'Cruise true airspeed', unit: 'm/s', default: 231, min: 5, group: 'Mission' },
  { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -40, max: 40, group: 'Atmosphere' },
];
const jetLike = (t) => t === 'turbofan' || t === 'turbojet';
const acDefaults = (c, up, d) => ({
  mtow: c.mass.mtow_kg, oew: c.mass.oew_kg, payload: c.mass.payload_kg, fuel_load: Math.max(0, Math.min(c.mass.fuel_kg, c.mass.mtow_kg - c.mass.oew_kg - c.mass.payload_kg)), S: c.wing.S_m2 || undefined, CD0: up.cfd?.CD0 ?? c.aero.CD0, k: up.cfd?.k_induced ?? (d.k_induced || undefined), CLmax: up.cfd?.CLmax ?? c.aero.CLmax_clean,
  ptype: ['turbofan', 'turbojet', 'turboprop', 'piston'].includes(c.prop.type) ? c.prop.type : 'turboprop', n_eng: c.prop.n_eng, T0_N: up.propulsion?.thrust_static_N ?? c.prop.T0_N, P0_W: c.prop.P0_W, bpr: c.prop.bpr, eta_prop: up.propeller?.eta_prop ?? c.prop.eta_prop, tsfc: up.propulsion?.tsfc_kg_Ns ?? c.prop.tsfc_kg_Ns, bsfc: c.prop.bsfc_kg_Ws, fuel: FLUIDS[c.prop.fuel]?.LHV ? c.prop.fuel : 'Jet A-1',
  range_km: c.mission.range_km || undefined, cruise_alt: c.mission.cruise_alt_m, cruise_V: c.mission.cruise_V_ms || c.flight.V_ms, dISA: c.atm.dISA_K,
});
const fuelWing = (c) => (c.wing.S_m2 > 0 && c.prop.type !== 'electric' && c.mass.fuel_kg > 0 ? true : 'This analysis is for fuel-burning fixed-wing aircraft. Use the electric mission or helicopter mission analysis for this vehicle.');
const ROTOR_IN = [
  { key: 'R', label: 'Rotor radius', unit: 'm', default: 8.18, min: 0.02, group: 'Rotor' }, { key: 'nR', label: 'Lifting rotors', unit: '', default: 1, min: 1, max: 16, step: 1, discrete: true, group: 'Rotor' },
  { key: 'sig', label: 'Rotor solidity', unit: '-', default: 0.0825, min: 0.005, max: 0.4, group: 'Rotor' }, { key: 'vt', label: 'Tip speed', unit: 'm/s', default: 221, min: 20, group: 'Rotor' },
  { key: 'cd0', label: 'Blade profile drag coefficient', unit: '-', default: 0.0095, min: 0.004, group: 'Rotor' }, { key: 'kappa', label: 'Induced power factor', unit: '-', default: 1.15, min: 1, max: 1.6, group: 'Rotor' },
  { key: 'f', label: 'Equivalent flat-plate area', unit: 'm²', default: 3.35, min: 0, group: 'Rotor' }, { key: 'eta', label: 'Transmission efficiency incl. anti-torque and accessories', unit: '-', default: 0.85, min: 0.5, max: 1, group: 'Rotor' },
];
const rotorDefaults = (c, d) => ({ R: c.rotor.R_m || undefined, nR: c.meta.type === 'helicopter' ? 1 : c.prop.n_eng, sig: d.solidity || undefined, vt: d.v_tip || undefined, cd0: c.rotor.cd0, f: c.rotor.flat_plate_m2 || undefined, eta: c.meta.type === 'helicopter' ? 0.85 : 0.92 });
/** Quick cruise consumption estimate straight from the case: fuel flow [kg/s], electrical power [W], speed [m/s]. */
function estCruise(c, d) {
  const V = c.mission.cruise_V_ms || c.flight.V_ms || 50, h = c.mission.cruise_alt_m || 0, rho = isa(h, c.atm.dISA_K).rho, m = c.mass.mtow_kg - 0.4 * c.mass.fuel_kg, elec = c.prop.type === 'electric'; let P;
  if (c.wing.S_m2 > 0 && d.k_induced > 0) { const D = dragOf({ S: c.wing.S_m2, CD0: c.aero.CD0, k: d.k_induced }, m, V, rho); if (!isProp(c.prop.type)) return { ff: c.prop.tsfc_kg_Ns * D, P: 0, V, elec: false }; P = (D * V) / c.prop.eta_prop; }
  else if (c.rotor.R_m > 0) P = rotorPower({ ...rotorDefaults(c, d), R: c.rotor.R_m, sig: d.solidity, vt: d.v_tip, f: c.rotor.flat_plate_m2, kappa: 1.15 }, m, V, rho).P; else P = 0.5 * d.P_total;
  return elec ? { ff: 0, P: P / 0.9, V, elec } : { ff: c.prop.bsfc_kg_Ws * P, P: 0, V, elec };
}

// ---- analyses -------------------------------------------------------------------------------
/** Full mission for a fuel-burning fixed-wing aircraft. Returns phase table and histories. */
function flyMission(i) {
  const p = i, n = Math.max(8, Math.round(i.nSteps)), fl = FLUIDS[i.fuel] || FLUIDS['Jet A-1'], warn = [], ffMax = maxFuelFlow(p), ffIdle = i.idle_frac * ffMax;
  const tom = i.oew + i.payload + i.fuel_load, H = { t: [0], h: [i.elev_dep], V: [0], m: [tom], ff: [ffIdle], x: [0] }, ph = [];
  let m = tom, t = 0, x = 0; const log = (name, fuel, time, dist) => { ph.push({ name, fuel, time, dist }); };
  const rec = (h, V, ff) => { H.t.push(t); H.h.push(h); H.V.push(V); H.m.push(m); H.ff.push(ff); H.x.push(x); };
  const sched = (h) => Math.min(i.climb_cas / Math.sqrt(isa(h, i.dISA).sigma), i.cruise_V);
  // taxi-out and take-off (full power to 450 m above the field at the initial climb speed)
  let f = ffIdle * i.taxi_min * 60; m -= f; t += i.taxi_min * 60; log('Taxi-out', f, i.taxi_min * 60, 0); rec(i.elev_dep, 0, ffIdle);
  const hStart = i.elev_dep + 450; f = ffMax * i.to_min * 60; m -= f; t += i.to_min * 60; x += 0.5 * sched(hStart) * i.to_min * 60; log('Take-off and initial climb', f, i.to_min * 60, 0.5 * sched(hStart) * i.to_min * 60); rec(hStart, sched(hStart), ffMax);
  // climb at the speed schedule with climb thrust; RK2 in time, energy correction for the accelerating TAS
  const rate = (h, mm) => { const V = sched(h), a = isa(h, i.dISA), T = i.climb_rating * thrustAvail(p, V, h, i.dISA), D = dragOf(p, mm, V, a.rho), dVdh = (sched(h + 50) - sched(h - 50)) / 100; return { roc: ((T - D) * V) / (mm * G0) / (1 + (V / G0) * dVdh), V, T, ff: fuelFlow(p, T, V) }; };
  let h = hStart, hc = Math.max(i.cruise_alt, hStart), r0 = rate(h, m); const dtc = Math.max(2, ((hc - h) / Math.max(0.5 * r0.roc, 0.5)) / n), c0 = { m, t, x };
  for (let s = 0; s < 6 * n && h < hc - 1e-6; s++) {
    const k1 = rate(h, m); if (k1.roc < i.roc_min) { warn.push(`Climb capability falls below ${i.roc_min} m/s at ${h.toFixed(0)} m, short of the planned ${hc.toFixed(0)} m; the cruise was flown at the altitude reached.`); hc = h; break; }
    const k2 = rate(h + 0.5 * k1.roc * dtc, m - 0.5 * k1.ff * dtc), dt = Math.min(dtc, (hc - h) / Math.max(k2.roc, 1e-6));
    h += k2.roc * dt; m -= k2.ff * dt; t += dt; x += Math.max(0, k2.V - i.wind_climb) * dt; rec(h, k2.V, k2.ff);
  }
  log('Climb', c0.m - m, t - c0.t, x - c0.x);
  // descent geometry first (fixed flight-path angle, idle thrust), so the cruise distance is known
  const hEnd = i.elev_dest + 450, nd = 16, gam = N.rad(i.desc_angle); let dDesc = 0, tDesc = 0; const dseg = [];
  for (let j = 0; j < nd; j++) { const hm = hc - ((j + 0.5) * (hc - hEnd)) / nd, V = sched(hm), dt = Math.max(0, hc - hEnd) / nd / (V * Math.sin(gam)); dseg.push({ hm, V, dt }); tDesc += dt; dDesc += Math.max(0, V * Math.cos(gam) - i.wind_desc) * dt; }
  let dCr = i.range_km * 1e3 - x - dDesc; if (dCr < 0) { warn.push('The stage is shorter than the climb and descent distances: the aircraft never reaches the planned cruise altitude on a real flight. Lower the cruise altitude for this stage.'); dCr = 0; }
  // cruise, with optional step climbs of 600 m when the heavier-altitude specific range is beaten and climb margin exists
  const c1 = { m, t, x }, nc = n, dx = dCr / nc, gs = Math.max(1, i.cruise_V - i.wind_cruise); let steps = 0, hcr = hc;
  for (let j = 0; j < nc && dCr > 0; j++) {
    if (i.step_climb && hcr + 600 <= i.alt_max && ffLevel(p, m, hcr + 600, i.cruise_V, i.dISA) < 0.995 * ffLevel(p, m, hcr, i.cruise_V, i.dISA)) { const a = isa(hcr + 600, i.dISA), V = i.cruise_V; if (((i.climb_rating * thrustAvail(p, V, hcr + 600, i.dISA) - dragOf(p, m, V, a.rho)) * V) / (m * G0) > 1.5 && (m * G0) / (0.5 * a.rho * V * V * p.S) < 0.9 * i.CLmax * 0.6) { m -= isProp(p.ptype) ? (p.bsfc * m * G0 * 600) / p.eta_prop : (p.tsfc * m * G0 * 600) / V; hcr += 600; steps++; } }
    const k1 = ffLevel(p, m, hcr, i.cruise_V, i.dISA) / gs, k2 = ffLevel(p, m - 0.5 * k1 * dx, hcr, i.cruise_V, i.dISA) / gs; m -= k2 * dx; t += dx / gs; x += dx; rec(hcr, i.cruise_V, k2 * gs);
  }
  log('Cruise', c1.m - m, t - c1.t, x - c1.x);
  const c2 = { m, t, x }, sc = (hcr - hEnd) / Math.max(1e-9, hc - hEnd);
  for (const s of dseg) { const dt = s.dt * (hc > hEnd ? sc : 0), hm = hEnd + (s.hm - hEnd) * sc; m -= ffIdle * 1.3 * dt; t += dt; x += Math.max(0, s.V * Math.cos(gam) - i.wind_desc) * dt; rec(hm, s.V, ffIdle * 1.3); }
  log('Descent', c2.m - m, t - c2.t, x - c2.x);
  f = 2.5 * ffIdle * i.appr_min * 60; m -= f; t += i.appr_min * 60; log('Approach and landing', f, i.appr_min * 60, 0); rec(i.elev_dest, 0, 2.5 * ffIdle);
  const trip = tom - m - ph[0].fuel, tAir = t - ph[0].time; f = 0.5 * ffIdle * i.taxi_min * 60; m -= f; t += 0.5 * i.taxi_min * 60; log('Taxi-in', f, 0.5 * i.taxi_min * 60, 0); rec(i.elev_dest, 0, ffIdle);
  // reserves evaluated at landing mass: contingency, alternate cruise, final hold at minimum-drag speed 450 m above the alternate
  const mL = m, hAlt = Math.min(hc, i.alt_alt), aA = isa(hAlt, i.dISA), Va = Math.min(i.cruise_V, Math.max(1.3 * vMinDrag(p, mL, aA.rho), 0.75 * i.cruise_V)), alt = i.alternate_km > 0 ? cruiseSeg(p, mL + 1, i.alternate_km * 1e3, hAlt, Va, 0, i.dISA, 12).fuel * 1.1 : 0;
  const aH = isa(i.elev_dest + 450, i.dISA), Vh = Math.max(vMinDrag(p, mL, aH.rho), 1.25 * Math.sqrt((2 * mL * G0) / (aH.rho * p.S * i.CLmax))), hold = ffLevel(p, mL, i.elev_dest + 450, Vh, i.dISA) * i.hold_min * 60, cont = (i.contingency_pct / 100) * trip;
  const block = tom - m;
  return { ph, H, tom, trip, tAir, block, tBlock: t, dist: x, reserves: { cont, alt, hold, total: cont + alt + hold }, margin: i.fuel_load - block - (cont + alt + hold), hc, hcr, steps, warn, fl, mL, CLc: (c1.m * G0) / (0.5 * isa(hc, i.dISA).rho * i.cruise_V ** 2 * p.S) };
}
const profile = {
  id: 'profile', title: 'Full mission simulation with fuel, time and reserves', fidelity: 'numerical',
  summary: 'Flies the stage phase by phase — taxi, take-off, climb on a speed schedule, cruise with optional step climbs, descent, approach and taxi-in — integrating mass and fuel, then checks the fuel left against contingency, alternate and holding reserves.',
  equations: ['Three-degree-of-freedom trajectory equations', 'Aircraft energy equations', 'Fuel consumption equations', 'Breguet range equation', 'Flight path kinematic equations', 'Atmospheric equations', 'Mission constraint equations'],
  applicable: fuelWing,
  inputs: [...AC, ...FLT,
    { key: 'climb_cas', label: 'Climb / descent calibrated airspeed', unit: 'm/s', default: 150, min: 5, group: 'Mission', help: 'Held until the cruise true airspeed is reached' },
    { key: 'climb_rating', label: 'Climb thrust / thrust available at altitude', unit: '-', default: 1, min: 0.5, max: 1, group: 'Mission' },
    { key: 'roc_min', label: 'Minimum acceptable rate of climb', unit: 'm/s', default: 0.5, min: 0.2, max: 5, group: 'Mission' },
    { key: 'step_climb', label: 'Allow step climbs in cruise', type: 'bool', default: false, group: 'Mission' }, { key: 'alt_max', label: 'Maximum operating altitude', unit: 'm', default: 12500, min: 100, max: 18000, group: 'Mission' },
    { key: 'desc_angle', label: 'Descent path angle', unit: 'deg', default: 3, min: 1, max: 8, group: 'Mission' },
    { key: 'taxi_min', label: 'Taxi-out time (taxi-in is half)', unit: 'min', default: 12, min: 0, max: 60, group: 'Ground and terminal' }, { key: 'to_min', label: 'Take-off and initial climb time at full power', unit: 'min', default: 1.5, min: 0.2, max: 5, group: 'Ground and terminal' },
    { key: 'appr_min', label: 'Approach and landing time', unit: 'min', default: 5, min: 0, max: 20, group: 'Ground and terminal' }, { key: 'idle_frac', label: 'Idle fuel flow / maximum fuel flow', unit: '-', default: 0.07, min: 0.01, max: 0.3, group: 'Ground and terminal' },
    { key: 'elev_dep', label: 'Departure elevation', unit: 'm', default: 0, min: -400, max: 4500, group: 'Atmosphere' }, { key: 'elev_dest', label: 'Destination elevation', unit: 'm', default: 0, min: -400, max: 4500, group: 'Atmosphere' },
    { key: 'wind_climb', label: 'Headwind in climb', unit: 'm/s', default: 0, min: -80, max: 80, group: 'Winds', help: 'Negative for tailwind; filled from live winds aloft when available' }, { key: 'wind_cruise', label: 'Headwind in cruise', unit: 'm/s', default: 0, min: -100, max: 100, group: 'Winds' }, { key: 'wind_desc', label: 'Headwind in descent', unit: 'm/s', default: 0, min: -80, max: 80, group: 'Winds' },
    { key: 'contingency_pct', label: 'Contingency fuel', unit: '% of trip', default: 5, min: 0, max: 20, group: 'Reserves' }, { key: 'alternate_km', label: 'Distance to alternate', unit: 'km', default: 185, min: 0, group: 'Reserves' },
    { key: 'alt_alt', label: 'Alternate cruise altitude', unit: 'm', default: 6000, min: 100, max: 12000, group: 'Reserves' }, { key: 'hold_min', label: 'Final reserve holding time', unit: 'min', default: 30, min: 0, max: 120, group: 'Reserves', help: '30 min for turbine aeroplanes and 45 min for piston aeroplanes are widely used rules' },
    { key: 'pax', label: 'Passengers', unit: '', default: 165, min: 0, step: 1, discrete: true, group: 'Mass' },
    { key: 'nSteps', label: 'Integration steps per phase', unit: '', default: 60, min: 8, max: 4000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => { const s = isa(c.mission.cruise_alt_m, 0).sigma, V = c.mission.cruise_V_ms || c.flight.V_ms, eas = V * Math.sqrt(s), vs = c.wing.S_m2 > 0 ? Math.sqrt((2 * c.mass.mtow_kg * G0) / (RHO0 * c.wing.S_m2 * c.aero.CLmax_clean)) : 0, small = c.mass.mtow_kg < 5700;
    const vmd = c.wing.S_m2 > 0 && d.k_induced > 0 ? Math.sqrt((2 * c.mass.mtow_kg * G0) / (RHO0 * c.wing.S_m2)) * (d.k_induced / c.aero.CD0) ** 0.25 : eas;
    return { ...acDefaults(c, up, d), climb_cas: Math.max(1.25 * vs, isProp(c.prop.type) ? 1.1 * vmd : Math.min(c.aero.Vmo_ms ? 0.85 * c.aero.Vmo_ms : eas, 1.05 * eas)), alt_max: Math.max(c.mission.cruise_alt_m * 1.17, c.mission.cruise_alt_m + 600), taxi_min: small ? 5 : 12, to_min: small ? 1 : 1.5, appr_min: small ? 3 : 5, idle_frac: c.prop.type === 'piston' ? 0.12 : 0.07,
      elev_dep: c.site.elev_m ?? 0, elev_dest: 0, alternate_km: c.mission.alternate_km, alt_alt: Math.min(6000, 0.6 * c.mission.cruise_alt_m + 300), hold_min: c.mission.reserve_min, pax: c.mission.pax, contingency_pct: 5 }; },
  run(i) {
    const r = flyMission(i), warnings = [...r.warn], co2 = r.block * (r.fl.co2_per_kg ?? 3.16), kWh = (r.block * r.fl.LHV) / 3.6e6, over = r.tom - i.mtow, feas = r.margin >= 0 && over <= 1e-6 && r.hc >= i.cruise_alt - 1 && r.ph[3].dist > 0;
    if (over > 1e-6) warnings.push(`Take-off mass exceeds the maximum by ${over.toFixed(0)} kg: offload payload or fuel.`);
    if (r.margin < 0) warnings.push(`Fuel on board is ${(-r.margin).toFixed(0)} kg short of trip fuel plus reserves.`);
    if (r.CLc > 0.6 * i.CLmax && jetLike(i.ptype)) warnings.push('Start-of-cruise lift coefficient is high for the cruise Mach number: buffet margin may be thin at this altitude and mass.');
    if (r.ph[3].dist > 0 && thrustAvail(i, i.cruise_V, r.hc, i.dISA) < dragOf(i, r.tom - r.ph[0].fuel - r.ph[1].fuel - r.ph[2].fuel, i.cruise_V, isa(r.hc, i.dISA).rho)) warnings.push('Drag at the start of cruise exceeds the thrust available at this altitude and speed: the aircraft would have to cruise slower or lower until fuel burns off. Fuel is still computed at the planned speed.');
    if (i.cruise_V / isa(i.cruise_alt, i.dISA).a > 0.82) warnings.push('Cruise Mach exceeds 0.82: wave drag is not in the parabolic polar, so fuel burn is optimistic.');
    const st = Math.max(1, Math.floor(r.H.t.length / 300)), ds = (a) => a.filter((_, j) => j % st === 0 || j === a.length - 1), tm = ds(r.H.t).map((v) => v / 60), pk = i.pax > 0 ? i.pax * r.dist / 1e3 : 0;
    return {
      kpis: [
        { key: 'block_fuel_kg', label: 'Block fuel', value: r.block, unit: 'kg' }, { key: 'trip_fuel_kg', label: 'Trip fuel (take-off to landing)', value: r.trip, unit: 'kg' },
        { key: 'block_time_h', label: 'Block time', value: r.tBlock / 3600, unit: 'h' }, { key: 'flight_time_h', label: 'Flight time', value: r.tAir / 3600, unit: 'h' },
        { key: 'reserve_required_kg', label: 'Reserves required', value: r.reserves.total, unit: 'kg', note: 'Contingency + alternate + final hold' },
        { key: 'reserve_margin_kg', label: 'Fuel margin above reserves', value: r.margin, unit: 'kg', status: r.margin >= 0 ? 'ok' : 'bad' },
        { key: 'mission_feasible', label: 'Mission feasible', value: feas ? 1 : 0, unit: '', status: feas ? 'ok' : 'bad', note: 'Mass limit, fuel with reserves, cruise altitude reached' },
        { key: 'co2_kg', label: 'CO₂ emitted (block)', value: co2, unit: 'kg' }, { key: 'mission_energy_kWh', label: 'Fuel energy used (block)', value: kWh, unit: 'kWh' },
        { key: 'tom_kg', label: 'Take-off mass', value: r.tom, unit: 'kg', status: over <= 1e-6 ? 'ok' : 'bad' }, { key: 'landing_mass_kg', label: 'Landing mass', value: r.mL, unit: 'kg' },
        { key: 'cruise_alt_final_m', label: 'Final cruise altitude', value: r.hcr, unit: 'm', note: r.steps ? `${r.steps} step climb(s)` : '' },
        { key: 'fuel_per_km_kg', label: 'Block fuel per kilometre', value: r.block / (r.dist / 1e3), unit: 'kg/km' },
        { key: 'co2_per_pax_km_g', label: 'CO₂ per passenger-kilometre', value: pk ? (1e3 * co2) / pk : 0, unit: 'g/pkm', note: pk ? 'All seats in the passenger count occupied' : 'No passengers defined' },
        { key: 'block_speed_ms', label: 'Block speed', value: r.dist / r.tBlock, unit: 'm/s' },
      ],
      plots: [
        { type: 'line', title: 'Altitude profile', xlabel: 'Time [min]', ylabel: 'Altitude [m]', series: [{ name: 'Altitude', x: tm, y: ds(r.H.h) }] },
        { type: 'line', title: 'True airspeed', xlabel: 'Time [min]', ylabel: 'True airspeed [m/s]', series: [{ name: 'TAS', x: tm, y: ds(r.H.V) }] },
        { type: 'line', title: 'Aircraft mass', xlabel: 'Time [min]', ylabel: 'Mass [kg]', series: [{ name: 'Mass', x: tm, y: ds(r.H.m) }], annotations: [{ y: i.oew + i.payload + r.reserves.total, label: 'Zero-fuel mass + reserves' }] },
        { type: 'line', title: 'Fuel flow', xlabel: 'Time [min]', ylabel: 'Fuel flow [kg/s]', series: [{ name: 'Fuel flow', x: tm, y: ds(r.H.ff), style: 'step' }] },
        { type: 'bar', title: 'Fuel by phase and reserves', ylabel: 'Fuel [kg]', categories: [...r.ph.map((q) => q.name), 'Contingency', 'Alternate', 'Final hold'], series: [{ name: 'Fuel', y: [...r.ph.map((q) => q.fuel), r.reserves.cont, r.reserves.alt, r.reserves.hold] }] },
      ],
      tables: [{ title: 'Mission phases', columns: ['Phase', 'Fuel [kg]', 'Time [min]', 'Distance [km]'], rows: [...r.ph.map((q) => [q.name, q.fuel, q.time / 60, q.dist / 1e3]), ['Block total', r.block, r.tBlock / 60, r.dist / 1e3]] }],
      outputs: { cruise_fuel_flow_kgs: r.ph[3].time > 0 ? r.ph[3].fuel / r.ph[3].time : 0, mission_dist_km: r.dist / 1e3 },
      warnings, models: ['Quasi-steady point-mass flight in the vertical plane', 'Parabolic drag polar', 'Thrust and power lapse with altitude and Mach number', 'Constant-CAS / constant-TAS speed schedule with acceleration correction', 'ISA with temperature offset'],
      assumptions: ['Constant specific fuel consumption in every phase', 'Idle-thrust descent on a fixed path angle; fuel flow 1.3× ground idle', 'Step climbs are instantaneous with the potential-energy fuel added', 'Reserves computed at landing mass; alternate includes 10% for the missed approach and climb', 'Winds are constant head- or tailwind components per phase'],
    };
  },
  convergence: { param: 'nSteps', label: 'Integration steps per phase', levels: [15, 30, 60, 120, 240], metric: 'block_fuel_kg' },
  verify() {
    const p = { ptype: 'turbojet', n_eng: 1, T0_N: 1e5, P0_W: 0, bpr: 0, eta_prop: 0.8, tsfc: 2e-5, bsfc: 0, S: 100, CD0: 0.02, k: 0.045 }, h = 9000, V = 220, m0 = 60000, dist = 3e6, q = 0.5 * isa(h).rho * V * V;
    const num = cruiseSeg(p, m0, dist, h, V, 0, 0, 400), m1 = m0 - num.fuel, c = Math.sqrt(p.k / p.CD0) / (q * p.S), exact = (V / (p.tsfc * G0 * Math.sqrt(p.CD0 * p.k))) * (Math.atan(m0 * G0 * c) - Math.atan(m1 * G0 * c));
    const i = Object.fromEntries(profile.inputs.map((f) => [f.key, f.default])), r = flyMission(i);
    return [
      N.check('Constant-altitude, constant-speed cruise range (arctangent Breguet form)', exact, dist, 1e-6, 'Closed-form integral of dR = −V·dW/(c·g·D) for a parabolic polar'),
      N.check('Fuel mass conservation: phase fuels sum to the mass change', N.sum(r.ph.map((x) => x.fuel)), r.tom - r.H.m[r.H.m.length - 1], 1e-10, 'Conservation of mass'),
      N.check('Phase distances sum to the stage length', N.sum(r.ph.map((x) => x.dist)), i.range_km * 1e3, 1e-9, 'Kinematic closure'),
      N.check('Headwind lengthens cruise time by V/(V − w)', cruiseSeg(p, m0, 1e6, h, V, 20, 0, 50).time, 1e6 / 200, 1e-12, 'Ground speed = TAS − headwind'),
    ];
  },
  calibration: { params: [{ key: 'CD0', min: 0.008, max: 0.08 }, { key: 'tsfc', min: 8e-6, max: 4e-5 }, { key: 'bsfc', min: 4e-8, max: 2e-7 }], sweep: 'range_km', target: 'block_fuel_kg', note: 'Supply recorded block fuel for several stage lengths at a known payload (operator or flight-recorder data).' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!o.mission_feasible) out.push({ severity: 'critical', title: 'The mission is not feasible as planned', detail: res.warnings.slice(0, 2).join(' '), action: o.reserve_margin_kg < 0 ? `Load ${(-o.reserve_margin_kg).toFixed(0)} kg more fuel if mass allows, otherwise offload payload, plan a technical stop or choose a closer alternate.` : 'Reduce take-off mass or lower the cruise altitude.', basis: 'Fuel policy: trip + contingency + alternate + final reserve; maximum take-off mass' });
    else if (o.reserve_margin_kg > 0.08 * i.fuel_load) out.push({ severity: 'advise', title: 'More fuel is loaded than this stage needs', detail: `${o.reserve_margin_kg.toFixed(0)} kg above the required reserves.`, action: `Carrying surplus fuel burns fuel: unloading it would save roughly ${(0.03 * o.reserve_margin_kg * o.flight_time_h).toFixed(0)} kg per flight (about 3% of the extra mass per flight hour) and the matching CO₂, unless tankering is cheaper.`, basis: 'Cost of weight' });
    if (!i.step_climb && o.flight_time_h > 3 && jetLike(i.ptype)) out.push({ severity: 'advise', title: 'Try step climbs', detail: 'On a long stage the optimum altitude rises as fuel burns off.', action: 'Enable step climbs and compare block fuel; 1–2% is typical on long sectors, a direct CO₂ saving.', basis: 'Specific range versus altitude and mass' });
    if (i.wind_cruise > 15) out.push({ severity: 'info', title: 'Strong headwind in cruise', detail: `${i.wind_cruise} m/s headwind.`, action: 'Check other flight levels and lateral tracks in the route analysis; flying slightly faster than still-air best-range speed pays in a headwind.', basis: 'Wind effect on ground specific range' });
    return out;
  },
};

const route = {
  id: 'route', title: 'Route, winds, terrain clearance and diversion distance', fidelity: 'analytical',
  summary: 'Turns a list of waypoints into great-circle legs with distances and bearings, applies a wind to get ground speeds and times, checks cruise altitude against a terrain profile, and finds the farthest point from any suitable airport.',
  equations: ['Flight path kinematic equations', 'Mission constraint equations', 'Atmospheric equations'],
  inputs: [
    { key: 'waypoints', label: 'Waypoints (lat, lon; lat, lon; …)', type: 'text', default: '48.0, 2.0; 52.2, 21.0; 55.6, 37.3', group: 'Route', help: 'Decimal degrees, north and east positive. The default is an illustrative route; it is replaced by the case site and design range when a site is set.' },
    { key: 'alternates', label: 'En-route alternates (lat, lon; …)', type: 'text', default: '', group: 'Route', help: 'Departure and destination are always counted as suitable airports' },
    { key: 'tas', label: 'Cruise true airspeed', unit: 'm/s', default: 231, min: 1, group: 'Flight' }, { key: 'cruise_alt', label: 'Cruise altitude', unit: 'm', default: 10668, min: 0, group: 'Flight' },
    { key: 'wind_dir', label: 'Wind direction (from)', unit: 'deg', default: 270, min: 0, max: 360, group: 'Wind', help: 'Filled from live winds aloft when available' }, { key: 'wind_ms', label: 'Wind speed', unit: 'm/s', default: 0, min: 0, max: 120, group: 'Wind' },
    { key: 'ff_kgs', label: 'Cruise fuel flow', unit: 'kg/s', default: 0.72, min: 0, group: 'Consumption', help: 'From Suite 5 or the mission simulation' }, { key: 'power_kW', label: 'Cruise electrical power (electric aircraft)', unit: 'kW', default: 0, min: 0, group: 'Consumption' },
    { key: 'terrain', label: 'Terrain elevations along the route', type: 'text', default: '', group: 'Terrain', help: 'Metres, evenly spaced from departure to destination; from live elevation data or typed in. Leave empty to skip the check.' },
    { key: 'clearance_min', label: 'Required terrain clearance', unit: 'm', default: 600, min: 0, max: 3000, group: 'Terrain', help: '300 m is common over flat terrain, 600 m in mountainous areas' },
    { key: 'oei_alt', label: 'Engine-out drift-down altitude (0 = skip)', unit: 'm', default: 0, min: 0, group: 'Terrain' },
    { key: 'v_div', label: 'Diversion speed (one engine inoperative)', unit: 'm/s', default: 200, min: 1, group: 'Diversion' }, { key: 'div_min', label: 'Diversion time limit', unit: 'min', default: 60, min: 5, max: 420, group: 'Diversion', help: '60 min for twins without extended-diversion approval' },
  ],
  defaults: (c, up, d) => { const e = estCruise(c, d), V = c.mission.cruise_V_ms || c.flight.V_ms, s = c.site, has = typeof s.lat === 'number' && typeof s.lon === 'number', R = Math.min(c.mission.range_km || 100, 15000) * 1e3, o = has ? [s.lat, s.lon] : [48, 2], mid = direct(o, 70, R / 2), dst = direct(o, 75, R), f = (p) => `${p[0].toFixed(4)}, ${p[1].toFixed(4)}`;
    return { waypoints: [o, mid, dst].map(f).join('; '), tas: V, cruise_alt: c.mission.cruise_alt_m, wind_dir: s.wind_dir_deg || 270, wind_ms: 0, ff_kgs: up.performance?.fuel_flow_cruise_kgs ?? e.ff, power_kW: e.P / 1e3, v_div: 0.85 * V, clearance_min: c.mission.cruise_alt_m < 1500 ? 150 : 600 }; },
  run(i) {
    const warnings = [], wp = parsePts(i.waypoints), al = parsePts(i.alternates); let pts = wp.pts;
    if (wp.bad) warnings.push(`${wp.bad} waypoint entr${wp.bad > 1 ? 'ies' : 'y'} could not be read and were skipped.`); if (al.bad) warnings.push(`${al.bad} alternate entr${al.bad > 1 ? 'ies' : 'y'} could not be read.`);
    if (pts.length < 2) { pts = [[48, 2], [52.2, 21]]; warnings.push('Fewer than two valid waypoints: an illustrative two-point route was used.'); }
    const legs = []; let dist = 0, time = 0;
    for (let j = 1; j < pts.length; j++) { const dv = vincenty(pts[j - 1], pts[j]), brg = bearing(pts[j - 1], pts[j]), w = windTriangle(brg, i.tas, i.wind_dir, i.wind_ms); legs.push({ d: dv, dh: haversine(pts[j - 1], pts[j]), brg, ...w, t: w.gs > 0.5 ? dv / w.gs : Infinity }); dist += dv; time += legs[j - 1].t; }
    if (!Number.isFinite(time)) warnings.push('The wind stops all progress on at least one leg (ground speed near zero).');
    // densified track for the map and the diversion scan
    const track = []; legs.forEach((l, j) => { const n = Math.max(2, Math.min(60, Math.ceil(l.d / 50e3))); const a = N.rad(pts[j][0]), b = N.rad(pts[j + 1][0]), la = N.rad(pts[j][1]), lb = N.rad(pts[j + 1][1]), dg = l.dh / R_E; for (let s = j ? 1 : 0; s <= n; s++) { const fr = s / n; if (dg < 1e-9) { track.push({ p: pts[j] }); continue; } const A = Math.sin((1 - fr) * dg) / Math.sin(dg), B = Math.sin(fr * dg) / Math.sin(dg), X = A * Math.cos(a) * Math.cos(la) + B * Math.cos(b) * Math.cos(lb), Y = A * Math.cos(a) * Math.sin(la) + B * Math.cos(b) * Math.sin(lb), Z = A * Math.sin(a) + B * Math.sin(b); track.push({ p: [N.deg(Math.atan2(Z, Math.hypot(X, Y))), N.deg(Math.atan2(Y, X))] }); } });
    const airports = [pts[0], pts[pts.length - 1], ...al.pts], far = track.map((q) => Math.min(...airports.map((ap) => haversine(q.p, ap)))), jf = N.argmax(far), dMax = far[jf], dLim = i.v_div * i.div_min * 60;
    const ter = parseNums(i.terrain), hasT = ter.vals.length >= 2, tMax = hasT ? N.amax(ter.vals) : NaN, clr = hasT ? i.cruise_alt - tMax : NaN, clrO = hasT && i.oei_alt > 0 ? i.oei_alt - tMax : NaN;
    if (ter.bad) warnings.push(`${ter.bad} terrain value(s) could not be read.`);
    if (!hasT) warnings.push('No terrain profile supplied: terrain clearance was not checked.');
    if (dist > 0.6 * Math.PI * R_E) warnings.push('Very long legs: a single wind vector for the whole route is a crude assumption.');
    const fuel = i.ff_kgs * time, kWh = (i.power_kW * time) / 3600, gcd = vincenty(pts[0], pts[pts.length - 1]), cum = N.cumtrapz(N.range(legs.length + 1), [0, ...legs.map((l) => l.d)]).map((_, j) => N.sum(legs.slice(0, j).map((l) => l.d)) / 1e3);
    return {
      kpis: [
        { key: 'route_dist_km', label: 'Route distance (WGS-84 geodesic)', value: dist / 1e3, unit: 'km' }, { key: 'direct_dist_km', label: 'Direct departure–destination distance', value: gcd / 1e3, unit: 'km' },
        { key: 'route_extension_pct', label: 'Route extension over the direct track', value: gcd > 0 ? 100 * (dist / gcd - 1) : 0, unit: '%', status: dist / Math.max(gcd, 1) < 1.08 ? 'ok' : 'warn' },
        { key: 'route_time_h', label: 'En-route time with wind', value: time / 3600, unit: 'h' }, { key: 'mean_gs_ms', label: 'Mean ground speed', value: time > 0 ? dist / time : 0, unit: 'm/s' },
        { key: 'wind_time_penalty_pct', label: 'Time change due to wind', value: 100 * (time / (dist / i.tas) - 1), unit: '%' },
        { key: 'route_fuel_kg', label: 'Cruise fuel for the route', value: fuel, unit: 'kg' }, { key: 'route_energy_kWh', label: 'Cruise electrical energy for the route', value: kWh, unit: 'kWh' },
        { key: 'max_wca_deg', label: 'Largest wind-correction angle', value: Math.max(...legs.map((l) => Math.abs(l.wca))), unit: 'deg' },
        ...(hasT ? [{ key: 'terrain_clearance_m', label: 'Smallest terrain clearance in cruise', value: clr, unit: 'm', status: clr >= i.clearance_min ? 'ok' : 'bad', note: `Highest terrain ${tMax.toFixed(0)} m` }] : []),
        ...(Number.isFinite(clrO) ? [{ key: 'oei_clearance_m', label: 'Terrain clearance at drift-down altitude', value: clrO, unit: 'm', status: clrO >= i.clearance_min ? 'ok' : 'bad' }] : []),
        { key: 'max_diversion_km', label: 'Farthest point from a suitable airport', value: dMax / 1e3, unit: 'km', status: dMax <= dLim ? 'ok' : 'warn' },
        { key: 'diversion_time_min', label: 'Diversion time from that point', value: dMax / i.v_div / 60, unit: 'min', status: dMax <= dLim ? 'ok' : 'warn', note: `Limit ${i.div_min} min` },
      ],
      plots: [
        { type: 'line', title: 'Route map (equirectangular)', xlabel: 'Longitude [deg]', ylabel: 'Latitude [deg]', series: [{ name: 'Great-circle track', x: track.map((q) => q.p[1]), y: track.map((q) => q.p[0]) }, { name: 'Waypoints', x: pts.map((q) => q[1]), y: pts.map((q) => q[0]), style: 'points' }, ...(al.pts.length ? [{ name: 'Alternates', x: al.pts.map((q) => q[1]), y: al.pts.map((q) => q[0]), style: 'points' }] : []), { name: 'Farthest from an airport', x: [track[jf].p[1]], y: [track[jf].p[0]], style: 'points' }] },
        { type: 'bar', title: 'Ground speed by leg', ylabel: 'Speed [m/s]', categories: legs.map((_, j) => `Leg ${j + 1}`), series: [{ name: 'Ground speed', y: legs.map((l) => l.gs) }, { name: 'True airspeed', y: legs.map(() => i.tas) }] },
        ...(hasT ? [{ type: 'line', title: 'Terrain and cruise altitude', xlabel: 'Distance [km]', ylabel: 'Altitude [m]', series: [{ name: 'Terrain', x: N.linspace(0, dist / 1e3, ter.vals.length), y: ter.vals }, { name: 'Cruise altitude', x: [0, dist / 1e3], y: [i.cruise_alt, i.cruise_alt] }, { name: 'Terrain + required clearance', x: N.linspace(0, dist / 1e3, ter.vals.length), y: ter.vals.map((v) => v + i.clearance_min), style: 'dash' }] }] : []),
      ],
      tables: [{ title: 'Legs', columns: ['Leg', 'Distance [km]', 'Sphere distance [km]', 'Initial bearing [deg]', 'Headwind [m/s]', 'Crosswind [m/s]', 'Wind correction [deg]', 'Ground speed [m/s]', 'Time [min]', 'Cumulative [km]'], rows: legs.map((l, j) => [j + 1, l.d / 1e3, l.dh / 1e3, l.brg, l.headwind, l.crosswind, l.wca, l.gs, l.t / 60, cum[j + 1]]) }],
      warnings, models: ['Vincenty inverse geodesic on the WGS-84 ellipsoid', 'Haversine great-circle distance and initial bearing', 'Wind triangle', 'Great-circle interpolation for the track', 'Nearest-airport diversion scan'],
      assumptions: ['One wind vector applies to the whole route at cruise level', 'Cruise fuel only: add climb, descent and reserves with the mission simulation', 'Terrain values are evenly spaced along the route', 'Diversion distance is still-air great-circle distance to the nearest listed airport'],
    };
  },
  verify() {
    const w = windTriangle(90, 100, 90, 20), x = windTriangle(0, 100, 90, 30);
    return [
      N.check('Haversine: quarter of the equator', haversine([0, 0], [0, 90]), (Math.PI * R_E) / 2, 1e-12, 'Great-circle arc'),
      N.check('Vincenty: one degree of longitude on the equator = a·π/180', vincenty([0, 0], [0, 1]), (WGS.a * Math.PI) / 180, 1e-10, 'Equatorial geodesic'),
      N.check('Vincenty: meridian quadrant of WGS-84', vincenty([0, 0], [90, 0]), 10001965.729, 1e-9, 'WGS-84 meridian quarter length'),
      N.check('Initial bearing due east along the equator', bearing([0, 10], [0, 20]), 90, 1e-12, 'Definition'),
      N.check('Pure headwind: GS = TAS − wind', w.gs, 80, 1e-12, 'Wind triangle'),
      N.check('Pure crosswind: GS = √(TAS² − wind²)', x.gs, Math.sqrt(100 * 100 - 900), 1e-12, 'Wind triangle'),
      N.check('Direct and inverse problems agree', haversine([10, 20], direct([10, 20], 37, 5e5)), 5e5, 1e-9, 'Spherical consistency'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.terrain_clearance_m !== undefined && o.terrain_clearance_m < i.clearance_min) out.push({ severity: 'critical', title: 'Cruise altitude does not clear terrain by the required margin', detail: `${o.terrain_clearance_m.toFixed(0)} m against ${i.clearance_min} m.`, action: 'Raise the cruise level or re-route around the high ground.', basis: 'Minimum obstacle clearance' });
    if (o.oei_clearance_m !== undefined && o.oei_clearance_m < i.clearance_min) out.push({ severity: 'warn', title: 'Engine-out drift-down does not clear terrain', detail: `${o.oei_clearance_m.toFixed(0)} m at the drift-down altitude.`, action: 'Define escape routes and decision points, or limit take-off mass to raise the engine-out ceiling.', basis: 'En-route one-engine-inoperative obstacle clearance' });
    if (o.diversion_time_min > i.div_min) out.push({ severity: 'warn', title: 'Part of the route is beyond the diversion time limit', detail: `${o.diversion_time_min.toFixed(0)} min from the nearest listed airport at ${i.v_div} m/s.`, action: 'Add suitable en-route alternates, re-route closer to airports, or operate under an extended-diversion-time approval.', basis: 'Diversion time rule' });
    if (o.route_extension_pct > 5) out.push({ severity: 'advise', title: 'The route is noticeably longer than the direct track', detail: `${o.route_extension_pct.toFixed(1)}% extension.`, action: `A more direct routing would save about ${(o.route_fuel_kg * o.route_extension_pct / (100 + o.route_extension_pct)).toFixed(0)} kg of fuel and its CO₂ on each flight.`, basis: 'Great-circle distance' });
    if (o.wind_time_penalty_pct > 5) out.push({ severity: 'info', title: 'Wind adds significant time', detail: `${o.wind_time_penalty_pct.toFixed(1)}% longer than still air.`, action: 'Compare other levels and tracks with live winds; carry the extra trip fuel.', basis: 'Wind triangle' });
    return out;
  },
};

/** Specific range [m/kg] of a fuel aircraft. */
const specRange = (p, m, h, V, dT) => V / ffLevel(p, m, h, V, dT);
const trade = {
  id: 'trade', title: 'Cruise optimisation, cost index and payload–range', fidelity: 'analytical',
  summary: 'Maps specific range over altitude and speed to find the maximum-range and long-range cruise points, shows how the economic speed moves with the cost index, and builds the payload–range diagram with reserves.',
  equations: ['Breguet range equation', 'Fuel consumption equations', 'Mission constraint equations', 'Optimal control equations'],
  applicable: fuelWing,
  inputs: [...AC, ...FLT,
    { key: 'mach_max', label: 'Maximum operating Mach', unit: '-', default: 0.82, min: 0.05, max: 0.95, group: 'Limits' }, { key: 'alt_max', label: 'Maximum operating altitude', unit: 'm', default: 12500, min: 500, max: 18000, group: 'Limits' },
    { key: 'cost_index', label: 'Cost index (time cost ÷ fuel cost)', unit: 'kg/min', default: 30, min: 0, max: 500, group: 'Economics', help: '0 = minimum fuel; high values favour speed' },
    { key: 'reserve_kg', label: 'Reserve and allowance fuel', unit: 'kg', default: 3200, min: 0, group: 'Mission', help: 'Contingency, alternate, hold, taxi' }, { key: 'climb_allow', label: 'Climb and descent fuel allowance', unit: '-', default: 0.03, min: 0, max: 0.2, group: 'Mission', help: 'Fraction of take-off mass' },
    { key: 'payload_max', label: 'Maximum payload', unit: 'kg', default: 16600, min: 0, group: 'Mass' }, { key: 'fuel_cap', label: 'Fuel capacity', unit: 'kg', default: 18800, min: 0, group: 'Mass' },
    { key: 'nGrid', label: 'Grid points per axis', unit: '', default: 31, min: 9, max: 121, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => { const e = estCruise(c, d), small = c.mass.mtow_kg < 5700; return { ...acDefaults(c, up, d), mach_max: c.aero.Mmo || 0.6, alt_max: Math.max(c.mission.cruise_alt_m * 1.17, c.mission.cruise_alt_m + 600, 1500), cost_index: small ? 2e-4 * c.mass.mtow_kg : 3.5e-4 * c.mass.mtow_kg, reserve_kg: e.ff * (c.mission.reserve_min * 60 + (c.mission.alternate_km * 1e3) / e.V) + 0.03 * c.mass.fuel_kg, payload_max: c.mass.payload_kg, fuel_cap: c.mass.fuel_kg }; },
  run(i) {
    const p = i, ng = Math.round(i.nGrid), mMid = Math.min(i.mtow, i.oew + i.payload + i.fuel_load) - 0.4 * i.fuel_load, hs = N.linspace(0, i.alt_max, ng), warnings = [];
    const lim = (h) => { const a = isa(h, i.dISA), Vs = 1.2 * Math.sqrt((2 * mMid * G0) / (a.rho * p.S * i.CLmax)), Vt = N.findRoot((V) => thrustAvail(p, V, h, i.dISA) - dragOf(p, mMid, V, a.rho), vMinDrag(p, mMid, a.rho), 1.2 * a.a, 60); return { a, Vs, Vhi: Math.min(i.mach_max * a.a, Number.isFinite(Vt) ? Vt : thrustAvail(p, vMinDrag(p, mMid, a.rho), h, i.dISA) > dragOf(p, mMid, vMinDrag(p, mMid, a.rho), a.rho) ? i.mach_max * a.a : 0) }; };
    const Vall = hs.map(lim), Vlo = Math.min(...Vall.map((l) => l.Vs)), Vhi = Math.max(...Vall.map((l) => l.Vhi), Vlo * 1.5), Vs = N.linspace(Vlo, Vhi, ng), SR = hs.map((h, a) => Vs.map((V) => (V >= Vall[a].Vs && V <= Vall[a].Vhi ? specRange(p, mMid, h, V, i.dISA) / 1e3 : 0)));
    // per-altitude optimum by golden section inside the admissible speed band
    const per = hs.map((h, a) => { const l = Vall[a]; if (l.Vhi <= l.Vs) return null; const V = N.goldenSection((v) => -specRange(p, mMid, h, v, i.dISA), l.Vs, l.Vhi, 1e-7), sr = specRange(p, mMid, h, V, i.dISA); let Vl = V; if (V < l.Vhi - 1e-6) { const f = (v) => specRange(p, mMid, h, v, i.dISA) - 0.99 * sr; Vl = f(l.Vhi) < 0 ? N.brent(f, V, l.Vhi, 1e-7) : l.Vhi; } return { h, V, sr, Vl, limited: V >= l.Vhi - 1e-4 }; });
    const ok = per.filter(Boolean); if (!ok.length) { return { kpis: [{ key: 'sr_max_km_kg', label: 'Best specific range', value: 0, unit: 'km/kg', status: 'bad' }], tables: [{ title: 'Status', columns: ['Message'], rows: [['No altitude offers level flight between stall margin and the thrust or Mach limit at this mass.']] }], warnings: ['Level cruise is not possible at this mass with the given thrust and limits.'], models: ['Parabolic drag polar'], assumptions: [] }; }
    const best = ok[N.argmax(ok.map((q) => q.sr))], srNow = specRange(p, mMid, i.cruise_alt, i.cruise_V, i.dISA), hb = best.h;
    // cost index: minimise (fuel flow + CI) / V at the best altitude
    const lb = lim(hb), econ = (ci) => N.goldenSection((v) => (ffLevel(p, mMid, hb, v, i.dISA) + ci / 60) / v, lb.Vs, lb.Vhi, 1e-7), CIs = N.linspace(0, Math.max(4 * i.cost_index, 1e-3), 25), Vci = CIs.map(econ), Ve = econ(i.cost_index);
    // payload–range with reserves: range from integrating specific range over the cruise burn at the planned cruise condition
    const rng = (tom, fuel) => { const use = fuel - i.reserve_kg - i.climb_allow * tom; if (use <= 0) return 0; const m0 = tom * (1 - 0.6 * i.climb_allow); return N.simpson((m) => specRange(p, m, i.cruise_alt, i.cruise_V, i.dISA), m0 - use, m0, 60) + 0.6 * i.climb_allow * tom * srNow * 0.5; };
    const fA = Math.max(0, Math.min(i.fuel_cap, i.mtow - i.oew - i.payload_max)), payB = Math.max(0, Math.min(i.payload_max, i.mtow - i.oew - i.fuel_cap)), pr = [[0, i.payload_max], [rng(i.oew + i.payload_max + fA, fA) / 1e3, i.payload_max]];
    if (payB < i.payload_max) pr.push([rng(i.oew + payB + i.fuel_cap, i.fuel_cap) / 1e3, payB]); pr.push([rng(i.oew + i.fuel_cap, i.fuel_cap) / 1e3, 0]);
    if (best.limited) warnings.push('The best-range speed is capped by the Mach or thrust limit at the optimum altitude.');
    if (jetLike(i.ptype)) warnings.push('Wave drag is not in the polar: the optimum for jets tends to sit on the Mach limit; treat speeds close to it as an upper bound.');
    const a0 = isa(hb, i.dISA).a;
    return {
      kpis: [
        { key: 'sr_max_km_kg', label: 'Best specific range', value: best.sr / 1e3, unit: 'km/kg' }, { key: 'opt_alt_m', label: 'Best-range altitude', value: hb, unit: 'm', note: hb >= i.alt_max - 1 ? 'At the altitude limit' : '' },
        { key: 'mrc_V_ms', label: 'Maximum-range cruise speed', value: best.V, unit: 'm/s' }, { key: 'mrc_mach', label: 'Maximum-range cruise Mach', value: best.V / a0, unit: '-' },
        { key: 'lrc_V_ms', label: 'Long-range cruise speed (99% of best range)', value: best.Vl, unit: 'm/s' }, { key: 'lrc_mach', label: 'Long-range cruise Mach', value: best.Vl / a0, unit: '-' },
        { key: 'econ_V_ms', label: 'Economic speed at the cost index', value: Ve, unit: 'm/s' }, { key: 'econ_fuel_penalty_pct', label: 'Fuel penalty of the economic speed', value: 100 * (best.sr / specRange(p, mMid, hb, Ve, i.dISA) - 1), unit: '%' },
        { key: 'sr_planned_km_kg', label: 'Specific range at the planned cruise', value: srNow / 1e3, unit: 'km/kg' }, { key: 'planned_vs_best_pct', label: 'Planned cruise versus best', value: 100 * (srNow / best.sr - 1), unit: '%', status: srNow > 0.96 * best.sr ? 'ok' : 'warn' },
        { key: 'range_max_payload_km', label: 'Range with maximum payload', value: pr[1][0], unit: 'km' }, { key: 'ferry_range_km', label: 'Ferry range', value: pr[pr.length - 1][0], unit: 'km' },
      ],
      plots: [
        { type: 'heat', title: 'Specific range over speed and altitude', xlabel: 'True airspeed [m/s]', ylabel: 'Altitude [m]', zlabel: 'Specific range [km/kg]', x: Vs, y: hs, z: SR, contours: 12, overlay: [{ name: 'Maximum-range speed', x: ok.map((q) => q.V), y: ok.map((q) => q.h) }, { name: 'Long-range cruise', x: ok.map((q) => q.Vl), y: ok.map((q) => q.h) }, { name: 'Planned cruise', x: [i.cruise_V], y: [i.cruise_alt], style: 'points' }] },
        { type: 'line', title: 'Economic speed versus cost index', xlabel: 'Cost index [kg/min]', ylabel: 'True airspeed [m/s]', series: [{ name: 'ECON speed', x: CIs, y: Vci }], annotations: [{ x: i.cost_index, label: 'Selected' }, { y: best.V, label: 'MRC' }] },
        { type: 'line', title: 'Payload–range diagram with reserves', xlabel: 'Range [km]', ylabel: 'Payload [kg]', series: [{ name: 'Payload–range boundary', x: pr.map((q) => q[0]), y: pr.map((q) => q[1]), style: 'line+points' }] },
      ],
      tables: [{ title: 'Optimum cruise by altitude', columns: ['Altitude [m]', 'MRC speed [m/s]', 'LRC speed [m/s]', 'Specific range [km/kg]', 'Speed limited'], rows: ok.filter((_, j) => j % Math.max(1, Math.floor(ok.length / 12)) === 0).map((q) => [q.h, q.V, q.Vl, q.sr / 1e3, q.limited ? 'yes' : '']) }],
      warnings, models: ['Specific range from the parabolic polar and constant specific fuel consumption', 'Golden-section search for maximum-range and economic speeds', 'Payload–range by integrating specific range over the burn'],
      assumptions: ['Evaluated at a mid-mission mass', 'Lower speed bound is 1.2× stall; upper bound is the Mach limit or thrust-limited speed', 'Payload–range holds the planned cruise altitude and speed, with fixed reserve fuel and a climb allowance'],
    };
  },
  verify() {
    const p = { ptype: 'turbojet', n_eng: 2, T0_N: 4e5, P0_W: 0, bpr: 0, eta_prop: 0.8, tsfc: 2e-5, bsfc: 0, S: 100, CD0: 0.02, k: 0.05 }, m = 50000, h = 8000, rho = isa(h).rho;
    const V = N.goldenSection((v) => -specRange(p, m, h, v, 0), 80, 330, 1e-10), q = { ...p, ptype: 'piston', P0_W: 1e6, bsfc: 8e-8 }, Vp = N.goldenSection((v) => -specRange(q, m, h, v, 0), 60, 300, 1e-10);
    return [
      N.check('Jet maximum-range speed = Vmd·3^¼', V, vMinDrag(p, m, rho) * 3 ** 0.25, 1e-6, 'CL = √(CD0/3k) for maximum V/D'),
      N.check('Propeller maximum-range speed = minimum-drag speed', Vp, vMinDrag(p, m, rho), 1e-6, 'Maximum L/D'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!o.sr_max_km_kg) return [{ severity: 'critical', title: 'No level cruise possible', detail: res.warnings[0], action: 'Reduce mass or check thrust and drag inputs.', basis: 'Thrust–drag balance' }];
    if (o.planned_vs_best_pct < -4) out.push({ severity: 'advise', title: 'The planned cruise point wastes fuel', detail: `Specific range is ${(-o.planned_vs_best_pct).toFixed(1)}% below the best available (${o.opt_alt_m.toFixed(0)} m, ${o.mrc_V_ms.toFixed(0)} m/s).`, action: `Fly closer to ${o.opt_alt_m.toFixed(0)} m at the long-range cruise speed of ${o.lrc_V_ms.toFixed(0)} m/s: the same percentage comes off fuel cost and CO₂.`, basis: 'Specific-range map' });
    out.push({ severity: 'info', title: 'Speed is a cost decision', detail: `At cost index ${i.cost_index.toFixed(0)} kg/min the economic speed is ${o.econ_V_ms.toFixed(0)} m/s, costing ${o.econ_fuel_penalty_pct.toFixed(1)}% more fuel than maximum-range cruise.`, action: 'Lower the cost index when fuel or carbon prices rise or when the schedule has slack; raise it only when time-related costs (crew, maintenance by the hour, missed connections) justify the burn.', basis: 'Minimum of (fuel flow + cost index) / speed' });
    return out;
  },
};

/** Energy-state climb: for each energy height pick the altitude that maximises `merit`; returns the path with time and fuel integrals. */
function energyPath(Ps, ff, merit, E0, E1, nE, hCand, ok) {
  const Es = N.linspace(E0, E1, nE), path = [];
  for (const E of Es) { let b = null; for (const h of hCand) { if (h > E) break; const V = Math.sqrt(2 * G0 * (E - h)); if (!ok(h, V)) continue; const ps = Ps(h, V); if (ps <= 0) continue; const v = merit(ps, ff(h, V)); if (!b || v > b.v) b = { v, h, V, ps, ff: ff(h, V) }; } path.push(b ? { E, ...b } : { E, h: NaN, V: NaN, ps: NaN, ff: NaN }); }
  let t = 0, fu = 0, done = true; for (let j = 1; j < path.length; j++) { const a = path[j - 1], b = path[j]; if (!(a.ps > 0) || !(b.ps > 0)) { done = false; break; } const dE = b.E - a.E; t += 0.5 * dE * (1 / a.ps + 1 / b.ps); fu += 0.5 * dE * (a.ff / a.ps + b.ff / b.ps); }
  return { path, time: t, fuel: fu, done };
}
const climbopt = {
  id: 'climbopt', title: 'Optimal climb by the energy-state method', fidelity: 'reduced-order',
  summary: 'Plots specific excess power over speed and altitude and finds the climb paths that reach the cruise energy in the least time or with the least fuel, compared with the constant-airspeed schedule.',
  equations: ['Aircraft energy equations', 'Optimal control equations', 'Dynamic programming equations', 'Three-degree-of-freedom trajectory equations'],
  applicable: fuelWing,
  inputs: [...AC, ...FLT.filter((f) => f.key !== 'range_km'),
    { key: 'climb_cas', label: 'Reference climb calibrated airspeed', unit: 'm/s', default: 150, min: 5, group: 'Mission' }, { key: 'climb_rating', label: 'Climb thrust / thrust available at altitude', unit: '-', default: 1, min: 0.5, max: 1, group: 'Mission' },
    { key: 'ps_floor', label: 'Specific excess power required at top of climb', unit: 'm/s', default: 0.5, min: 0.2, max: 5, group: 'Limits' },
    { key: 'mach_max', label: 'Maximum operating Mach', unit: '-', default: 0.82, min: 0.05, max: 0.95, group: 'Limits' }, { key: 'vmo_eas', label: 'Maximum operating speed (EAS)', unit: 'm/s', default: 180, min: 5, group: 'Limits' },
    { key: 'nGrid', label: 'Grid points per axis', unit: '', default: 40, min: 10, max: 160, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => ({ ...acDefaults(c, up, d), climb_cas: profile.defaults(c, up, d).climb_cas, alt_max: undefined, mach_max: c.aero.Mmo || 0.6, vmo_eas: c.aero.Vmo_ms || 1.3 * (c.mission.cruise_V_ms || c.flight.V_ms) }),
  run(i) {
    const p = i, m = Math.min(i.mtow, i.oew + i.payload + i.fuel_load), W = m * G0, ng = Math.round(i.nGrid), warnings = []; let hTop = i.cruise_alt;
    const Ps = (h, V) => { const a = isa(h, i.dISA); return ((i.climb_rating * thrustAvail(p, V, h, i.dISA) - dragOf(p, m, V, a.rho)) * V) / W; }, ff = (h, V) => fuelFlow(p, i.climb_rating * thrustAvail(p, V, h, i.dISA), V);
    const ok = (h, V) => { const a = isa(h, i.dISA); return V >= 1.2 * Math.sqrt((2 * W) / (a.rho * p.S * i.CLmax)) && V <= i.mach_max * a.a && V * Math.sqrt(a.sigma) <= i.vmo_eas; };
    const sched = (h) => Math.min(i.climb_cas / Math.sqrt(isa(h, i.dISA).sigma), i.cruise_V);
    while (hTop > 200 && !(Ps(hTop, sched(hTop)) >= i.ps_floor)) hTop -= Math.max(50, 0.01 * i.cruise_alt); // climb target: highest altitude with a usable rate on the schedule
    if (hTop < i.cruise_alt - 1) warnings.push(`At this mass less than ${i.ps_floor} m/s of specific excess power is left on the climb schedule above ${hTop.toFixed(0)} m; the climb is evaluated to that altitude instead of ${i.cruise_alt.toFixed(0)} m.`);
    const Vend = Ps(hTop, i.cruise_V) > 0.1 && ok(hTop, i.cruise_V) ? i.cruise_V : sched(hTop);
    if (Vend < i.cruise_V - 1e-9) warnings.push('The planned cruise speed cannot be reached at the top of climb at this mass (no excess power, or outside the speed limits); the climb ends at the schedule speed.');
    const hs = N.linspace(0, Math.max(i.cruise_alt * 1.1, 500), ng), Vmax = Math.max(i.cruise_V * 1.15, i.vmo_eas), Vs = N.linspace(0.5 * i.climb_cas, Vmax, ng), Z = hs.map((h) => Vs.map((V) => (ok(h, V) ? Math.max(0, Ps(h, V)) : 0)));
    const V0 = sched(0), E0 = (V0 * V0) / (2 * G0), E1 = hTop + Vend ** 2 / (2 * G0), hc = N.linspace(0, hTop, 4 * ng), nE = 3 * ng;
    const okS = (h, V) => ok(h, V) || Math.abs(V - sched(h)) < 1e-9, tmin = energyPath(Ps, ff, (ps) => ps, E0, E1, nE, hc, okS), fmin = energyPath(Ps, ff, (ps, f) => ps / f, E0, E1, nE, hc, okS);
    // reference: constant-CAS / cruise-TAS schedule to the top, then level acceleration to cruise speed; same end energy
    const hr = N.linspace(0, hTop, nE), ref = [...hr.map((h) => [h, sched(h)]), ...(Vend > sched(hTop) ? N.linspace(sched(hTop), Vend, 8).slice(1).map((V) => [hTop, V]) : [])]; let tr = 0, fr = 0, refOk = true;
    for (let j = 1; j < ref.length; j++) { const a = { h: ref[j - 1][0], V: ref[j - 1][1] }, b = { h: ref[j][0], V: ref[j][1] }, pa = Ps(a.h, a.V), pb = Ps(b.h, b.V); if (pa <= 0 || pb <= 0) { refOk = false; break; } const dE = b.h - a.h + (b.V ** 2 - a.V ** 2) / (2 * G0); tr += 0.5 * dE * (1 / pa + 1 / pb); fr += 0.5 * dE * (ff(a.h, a.V) / pa + ff(b.h, b.V) / pb); }
    if (!tmin.done) warnings.push('Specific excess power runs out before the cruise energy is reached: the cruise altitude and speed are not attainable at this mass and thrust rating.');
    if (!refOk) warnings.push('The reference speed schedule cannot reach the cruise altitude.');
    const fin = (pth) => pth.filter((q) => Number.isFinite(q.h));
    return {
      kpis: [
        { key: 'climb_time_min_s', label: 'Minimum time to cruise energy', value: tmin.done ? tmin.time : NaN, unit: 's', status: tmin.done ? 'ok' : 'bad' }, { key: 'climb_fuel_at_min_time_kg', label: 'Fuel on the minimum-time path', value: tmin.done ? tmin.fuel : NaN, unit: 'kg' },
        { key: 'climb_fuel_min_kg', label: 'Minimum fuel to cruise energy', value: fmin.done ? fmin.fuel : NaN, unit: 'kg' }, { key: 'climb_time_at_min_fuel_s', label: 'Time on the minimum-fuel path', value: fmin.done ? fmin.time : NaN, unit: 's' },
        { key: 'climb_time_ref_s', label: 'Time on the reference schedule', value: refOk ? tr : NaN, unit: 's' }, { key: 'climb_fuel_ref_kg', label: 'Fuel on the reference schedule', value: refOk ? fr : NaN, unit: 'kg' },
        { key: 'time_saving_pct', label: 'Time saved by the optimal path', value: refOk && tmin.done ? 100 * (1 - tmin.time / tr) : NaN, unit: '%' }, { key: 'fuel_saving_pct', label: 'Fuel saved by the optimal path', value: refOk && fmin.done ? 100 * (1 - fmin.fuel / fr) : NaN, unit: '%' },
        { key: 'ps_max_sl_ms', label: 'Maximum specific excess power at sea level', value: Math.max(...Z[0]), unit: 'm/s' }, { key: 'climb_top_m', label: 'Top of climb evaluated', value: hTop, unit: 'm', status: hTop >= i.cruise_alt - 1 ? 'ok' : 'warn' },
      ],
      plots: [
        { type: 'heat', title: 'Specific excess power with optimal climb paths', xlabel: 'True airspeed [m/s]', ylabel: 'Altitude [m]', zlabel: 'Ps [m/s]', x: Vs, y: hs, z: Z, contours: 12, overlay: [{ name: 'Minimum time', x: fin(tmin.path).map((q) => q.V), y: fin(tmin.path).map((q) => q.h) }, { name: 'Minimum fuel', x: fin(fmin.path).map((q) => q.V), y: fin(fmin.path).map((q) => q.h) }, { name: 'Reference schedule', x: ref.map((q) => q[1]), y: ref.map((q) => q[0]) }] },
        { type: 'line', title: 'Specific excess power along each path', xlabel: 'Energy height [m]', ylabel: 'Ps [m/s]', series: [{ name: 'Minimum time', x: fin(tmin.path).map((q) => q.E), y: fin(tmin.path).map((q) => q.ps) }, { name: 'Minimum fuel', x: fin(fmin.path).map((q) => q.E), y: fin(fmin.path).map((q) => q.ps) }] },
      ],
      warnings, models: ['Energy-state approximation: total energy height as the single state', 'Specific excess power from thrust lapse and the parabolic polar', 'Pointwise maximisation of Ps (minimum time) and Ps per unit fuel flow (minimum fuel)'],
      assumptions: ['Mass constant during the climb', 'Kinetic and potential energy can be exchanged instantly (zooms and dives along constant-energy lines cost nothing)', 'Lift equals weight', 'Speed limits: 1.2× stall, VMO and MMO'],
    };
  },
  convergence: { param: 'nGrid', label: 'Grid points per axis', levels: [10, 20, 40, 80], metric: 'climb_time_min_s' },
  verify() {
    const e = energyPath(() => 12, () => 2, (ps) => ps, 1000, 9000, 40, N.linspace(0, 8000, 50), () => true);
    const p = { ptype: 'turbojet', n_eng: 1, T0_N: 5e4, P0_W: 0, bpr: 0, S: 30, CD0: 0.02, k: 0.06, eta_prop: 0.8 }, V = 150, m = 8000, a = isa(3000), ps = ((0.9 * thrustAvail(p, V, 3000) - dragOf(p, m, V, a.rho)) * V) / (m * G0);
    const o = N.kv(climbopt.run({ ...Object.fromEntries(climbopt.inputs.map((f) => [f.key, f.default])) }));
    return [
      N.check('Constant Ps: time = ΔE / Ps', e.time, 8000 / 12, 1e-12, 'Energy-state integral'),
      N.check('Constant Ps and fuel flow: fuel = ṁ·ΔE / Ps', e.fuel, (2 * 8000) / 12, 1e-12, 'Energy-state integral'),
      N.check('Ps = V(T − D)/W at a point', ps, (V * (0.9 * 5e4 * Math.max(0.05, a.delta * (1 + 0.2 * (V / a.a) ** 2) ** 3.5 * (1 - 0.16 * Math.sqrt(V / a.a))) - (0.5 * a.rho * V * V * 30 * 0.02 + (0.06 * (m * G0) ** 2) / (0.5 * a.rho * V * V * 30)))) / (m * G0), 1e-12, 'Definition'),
      N.check('Optimal path is never slower than the reference schedule', Math.min(o.climb_time_min_s, o.climb_time_ref_s) / o.climb_time_min_s, 1, 5e-3, 'Optimality (grid tolerance)'),
    ];
  },
  recommend(res) {
    const o = res.outputs, out = [];
    if (!Number.isFinite(o.climb_time_min_s)) return [{ severity: 'critical', title: 'Cruise condition is not reachable', detail: 'Specific excess power falls to zero below the cruise energy.', action: 'Lower the cruise altitude or take-off mass, or raise the climb rating.', basis: 'Ps = 0 boundary' }];
    if (o.fuel_saving_pct > 2) out.push({ severity: 'advise', title: 'The climb schedule can be improved', detail: `The minimum-fuel path uses ${o.fuel_saving_pct.toFixed(1)}% less climb fuel than the reference schedule; minimum time is ${o.time_saving_pct.toFixed(1)}% quicker.`, action: 'Move the climb speed towards the path shown on the Ps map (typically a little faster low down, easing to best-climb speed aloft). Climb fuel saved is CO₂ saved on every departure.', basis: 'Energy-state optimal climb' });
    else out.push({ severity: 'info', title: 'The reference climb schedule is near-optimal', detail: `Within ${Math.max(0, o.fuel_saving_pct || 0).toFixed(1)}% of minimum fuel.`, action: 'No change needed; review again if mass or thrust rating changes.', basis: 'Energy-state optimal climb' });
    return out;
  },
};

const heli = {
  id: 'heli', title: 'Helicopter missions: offshore, search and rescue, external load', fidelity: 'reduced-order',
  summary: 'Builds the rotor power-required curve, then flies an out-and-back mission with hover segments and time on station, giving fuel, radius of action, point of no return and the payload–radius trade.',
  equations: ['Aircraft energy equations', 'Fuel consumption equations', 'Mission constraint equations', 'Flight path kinematic equations'],
  applicable: (c) => (c.rotor.R_m > 0 && c.wing.S_m2 === 0 && c.prop.type !== 'electric' ? true : 'This analysis is for fuel-burning rotorcraft without a wing.'),
  inputs: [
    { key: 'mission', label: 'Mission type', type: 'select', options: ['Offshore transport', 'Search and rescue', 'External load'], default: 'Offshore transport', group: 'Mission' },
    { key: 'mtow', label: 'Maximum take-off mass', unit: 'kg', default: 9980, min: 1, group: 'Mass' }, { key: 'oew', label: 'Operating empty mass', unit: 'kg', default: 5220, min: 1, group: 'Mass' }, { key: 'payload', label: 'Payload', unit: 'kg', default: 3600, min: 0, group: 'Mass' }, { key: 'fuel_cap', label: 'Fuel capacity', unit: 'kg', default: 1100, min: 0, group: 'Mass' },
    ...ROTOR_IN,
    { key: 'P_inst', label: 'Total installed power', unit: 'W', default: 2.82e6, min: 1, group: 'Propulsion' }, { key: 'bsfc', label: 'BSFC', unit: 'kg/W/s', default: 7.8e-8, min: 1e-8, group: 'Propulsion' }, { key: 'idle_frac', label: 'Ground idle fuel flow / maximum', unit: '-', default: 0.2, min: 0.02, max: 0.5, group: 'Propulsion' },
    { key: 'radius_km', label: 'Mission radius', unit: 'km', default: 250, min: 0, group: 'Mission' }, { key: 'V', label: 'Cruise speed', unit: 'm/s', default: 72, min: 5, max: 110, group: 'Mission' }, { key: 'alt', label: 'Cruise altitude', unit: 'm', default: 500, min: 0, max: 6000, group: 'Mission' }, { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -40, max: 40, group: 'Mission' },
    { key: 'wind_out', label: 'Headwind outbound (tailwind home)', unit: 'm/s', default: 0, min: -30, max: 30, group: 'Mission' }, { key: 'hover_min', label: 'Total hover time', unit: 'min', default: 10, min: 0, max: 240, group: 'Mission', help: 'Split evenly between departure, destination and return' },
    { key: 'ground_min', label: 'Rotors-running ground time', unit: 'min', default: 8, min: 0, max: 60, group: 'Mission' }, { key: 'payload_drop', label: 'Payload left at the destination', unit: 'kg', default: 0, min: 0, group: 'Mission', help: 'External load or passengers not carried back' },
    { key: 'reserve_min', label: 'Final reserve at best-endurance speed', unit: 'min', default: 30, min: 0, max: 90, group: 'Reserves' }, { key: 'contingency_pct', label: 'Contingency fuel', unit: '% of trip', default: 10, min: 0, max: 30, group: 'Reserves' },
    { key: 'fuel', label: 'Fuel', type: 'select', options: ['Jet A-1', 'Avgas 100LL', 'SAF (HEFA-SPK)'], default: 'Jet A-1', group: 'Propulsion' },
  ],
  defaults: (c, up, d) => ({ ...rotorDefaults(c, d), mission: /sar|search/i.test(c.mission.profile) ? 'Search and rescue' : /external|load/i.test(c.mission.profile) ? 'External load' : 'Offshore transport', mtow: c.mass.mtow_kg, oew: c.mass.oew_kg, payload: c.mass.payload_kg, fuel_cap: c.mass.fuel_kg, P_inst: d.P_total, bsfc: c.prop.bsfc_kg_Ws,
    radius_km: 0.3 * c.mission.range_km, V: c.mission.cruise_V_ms || c.flight.V_ms, alt: c.mission.cruise_alt_m, dISA: c.atm.dISA_K, hover_min: c.mission.hover_min || 6, reserve_min: Math.min(c.mission.reserve_min, 30), fuel: c.prop.fuel === 'Avgas 100LL' ? c.prop.fuel : 'Jet A-1' }),
  run(i) {
    const r = i, rho = isa(i.alt, i.dISA).rho, rho0 = isa(0, i.dISA).rho, fl = FLUIDS[i.fuel] || FLUIDS['Jet A-1'], warnings = [], sar = i.mission === 'Search and rescue';
    const fuelMax = Math.max(0, Math.min(i.fuel_cap, i.mtow - i.oew - i.payload)), tom = i.oew + i.payload + fuelMax, ffAt = (m, V, rh = rho) => i.bsfc * rotorPower(r, m, V, rh).P;
    /** Fly out, work and return for a given radius; returns phase list, with optional on-station time at best-endurance speed. */
    const fly = (radius, payload, fuel0, tStation = 0) => {
      const ph = []; let m = i.oew + payload + fuel0, t = 0; const seg = (name, ff, dt, dist = 0) => { const f = ff * dt; ph.push({ name, fuel: f, time: dt, dist }); m -= f; t += dt; };
      seg('Start and ground running', i.idle_frac * i.bsfc * i.P_inst, i.ground_min * 60); seg('Hover and departure', ffAt(m, 0, rho0), (i.hover_min * 60) / 3);
      const leg = (name, wind) => { const gs = Math.max(1, i.V - wind), dt = radius / gs, f1 = ffAt(m, i.V), f2 = ffAt(m - 0.5 * f1 * dt, i.V); seg(name, f2, dt, radius); };
      leg('Cruise out', i.wind_out); if (tStation > 0) seg('On station at best-endurance speed', ffAt(m, rotorSpeeds(r, m, rho).Vbe), tStation);
      seg('Hover at destination', ffAt(m, 0, rho0), (i.hover_min * 60) / 3); m -= Math.min(i.payload_drop, payload); leg('Cruise home', -i.wind_out); seg('Hover and landing', ffAt(m, 0, rho0), (i.hover_min * 60) / 3);
      const trip = N.sum(ph.map((q) => q.fuel)), sp = rotorSpeeds(r, m, rho), res = sp.Pbe * i.bsfc * i.reserve_min * 60 + (i.contingency_pct / 100) * trip;
      return { ph, trip, time: t, res, mEnd: m, margin: fuel0 - trip - res };
    };
    const base0 = fly(i.radius_km * 1e3, i.payload, fuelMax);
    // SAR: longest search time that still leaves the reserves
    let tSt = 0; if (sar && base0.margin > 0) tSt = N.brent((ts) => fly(i.radius_km * 1e3, i.payload, fuelMax, ts).margin, 0, 20 * 3600, 1);
    const b = fly(i.radius_km * 1e3, i.payload, fuelMax, tSt), sp = rotorSpeeds(r, tom, rho, i.wind_out), Ph = rotorPower(r, tom, 0, rho0).P, pm = (i.P_inst * isa(0, i.dISA).sigma ** 0.7 - Ph) / (i.P_inst * isa(0, i.dISA).sigma ** 0.7);
    const radMax = (pay) => { const f0 = Math.max(0, Math.min(i.fuel_cap, i.mtow - i.oew - pay)), g = (R) => fly(R, pay, f0).margin; return g(0) <= 0 ? 0 : g(2e6) > 0 ? 2e6 : N.brent(g, 0, 2e6, 10); };
    const pays = N.linspace(0, Math.max(i.payload, i.mtow - i.oew - 0.2 * i.fuel_cap), 12), rads = pays.map((q) => radMax(q) / 1e3), rMaxNow = radMax(i.payload);
    // radius of action and point of no return on usable fuel (after hover, ground and reserves), at cruise fuel flow
    const ffc = ffAt(tom - 0.5 * fuelMax, i.V), usable = Math.max(0, fuelMax - b.res - (b.ph[0].fuel + b.ph[1].fuel + b.ph[b.ph.length - 1].fuel)), E = usable / ffc, gso = Math.max(1, i.V - i.wind_out), gsh = Math.max(1, i.V + i.wind_out), roa = (E * gso * gsh) / (gso + gsh), tPnr = (E * gsh) / (gso + gsh);
    const feas = b.margin >= -1e-6 && pm > 0 && tom <= i.mtow + 1e-6;
    if (pm <= 0) warnings.push('Hover power required exceeds installed power at take-off mass: the mission cannot start out of ground effect at this mass and temperature.');
    if (b.margin < 0) warnings.push(`Fuel is ${(-b.margin).toFixed(0)} kg short of trip plus reserves for this radius; the largest radius with this payload is ${(rMaxNow / 1e3).toFixed(0)} km.`);
    if (fuelMax < i.fuel_cap - 1e-6) warnings.push(`Tanks cannot be filled: fuel is limited to ${fuelMax.toFixed(0)} kg by maximum take-off mass with this payload.`);
    if (i.V / i.vt > 0.4) warnings.push('Advance ratio above 0.4: retreating-blade stall and compressibility are not modelled and power is optimistic.');
    const block = b.trip, co2 = block * (fl.co2_per_kg ?? 3.16);
    return {
      kpis: [
        { key: 'block_fuel_kg', label: 'Mission fuel', value: block, unit: 'kg' }, { key: 'block_time_h', label: 'Mission time', value: b.time / 3600, unit: 'h' },
        { key: 'reserve_margin_kg', label: 'Fuel margin above reserves', value: b.margin, unit: 'kg', status: b.margin >= -1e-6 ? 'ok' : 'bad' }, { key: 'mission_feasible', label: 'Mission feasible', value: feas ? 1 : 0, unit: '', status: feas ? 'ok' : 'bad' },
        { key: 'co2_kg', label: 'CO₂ emitted', value: co2, unit: 'kg' }, { key: 'mission_energy_kWh', label: 'Fuel energy used', value: (block * fl.LHV) / 3.6e6, unit: 'kWh' },
        { key: 'time_on_station_min', label: 'Time on station', value: tSt / 60, unit: 'min', note: sar ? 'Longest search time leaving the reserves' : 'Search and rescue missions only' },
        { key: 'radius_max_km', label: 'Largest radius with this payload', value: rMaxNow / 1e3, unit: 'km' }, { key: 'radius_of_action_km', label: 'Radius of action with wind', value: roa / 1e3, unit: 'km' }, { key: 'pnr_time_min', label: 'Time to the point of no return', value: tPnr / 60, unit: 'min' },
        { key: 'V_be_ms', label: 'Best-endurance speed', value: sp.Vbe, unit: 'm/s' }, { key: 'V_br_ms', label: 'Best-range speed (with wind)', value: sp.Vbr, unit: 'm/s' },
        { key: 'hover_power_W', label: 'Hover power at take-off mass', value: Ph, unit: 'W' }, { key: 'hover_power_margin_pct', label: 'Hover power margin', value: 100 * pm, unit: '%', status: pm > 0.1 ? 'ok' : pm > 0 ? 'warn' : 'bad' },
        { key: 'cruise_fuel_flow_kgs', label: 'Cruise fuel flow', value: ffc, unit: 'kg/s' }, { key: 'hover_fuel_flow_kgs', label: 'Hover fuel flow', value: i.bsfc * Ph, unit: 'kg/s' },
      ],
      plots: [
        { type: 'line', title: 'Power required in level flight at take-off mass', xlabel: 'True airspeed [m/s]', ylabel: 'Power [kW]', series: [{ name: 'Total', x: sp.Vs, y: sp.P.map((v) => v / 1e3) }, { name: 'Induced', x: sp.Vs, y: sp.Vs.map((V) => rotorPower(r, tom, V, rho).Pi / r.eta / 1e3), style: 'dash' }, { name: 'Profile', x: sp.Vs, y: sp.Vs.map((V) => rotorPower(r, tom, V, rho).Po / r.eta / 1e3), style: 'dash' }, { name: 'Parasite', x: sp.Vs, y: sp.Vs.map((V) => rotorPower(r, tom, V, rho).Pp / r.eta / 1e3), style: 'dash' }], annotations: [{ x: sp.Vbe, label: 'Best endurance' }, { x: sp.Vbr, label: 'Best range' }, { y: (i.P_inst * isa(i.alt, i.dISA).sigma ** 0.7) / 1e3, label: 'Available' }] },
        { type: 'line', title: 'Payload–radius', xlabel: 'Radius [km]', ylabel: 'Payload [kg]', series: [{ name: 'Payload–radius boundary', x: rads, y: pays, style: 'line+points' }, { name: 'This mission', x: [i.radius_km], y: [i.payload], style: 'points' }] },
        { type: 'bar', title: 'Fuel by phase', ylabel: 'Fuel [kg]', categories: [...b.ph.map((q) => q.name), 'Reserves'], series: [{ name: 'Fuel', y: [...b.ph.map((q) => q.fuel), b.res] }] },
      ],
      tables: [{ title: 'Mission phases', columns: ['Phase', 'Fuel [kg]', 'Time [min]', 'Distance [km]'], rows: [...b.ph.map((q) => [q.name, q.fuel, q.time / 60, q.dist / 1e3]), ['Total', b.trip, b.time / 60, (2 * i.radius_km)]] }],
      outputs: { tom_kg: tom, fuel_loaded_kg: fuelMax },
      warnings, models: ['Momentum theory with Glauert forward-flight inflow', 'Profile power with advance-ratio growth (1 + 4.65μ²)', 'Parasite power from equivalent flat-plate area', 'Constant brake specific fuel consumption', 'Radius of action and point of no return with wind'],
      assumptions: ['Transmission efficiency covers tail rotor and accessories', 'Hover out of ground effect at sea-level density of the day', 'Fuel loaded to capacity or to the maximum take-off mass limit', 'Specific fuel consumption does not rise at part power (optimistic at low power)', 'No blade stall or compressibility limits'],
    };
  },
  verify() {
    const r = { R: 6, nR: 1, sig: 0.07, vt: 210, cd0: 0, kappa: 1, f: 0, eta: 1 }, m = 4000, rho = 1.2, T = m * G0, A = Math.PI * 36, h = rotorPower(r, m, 0, rho), f = rotorPower(r, m, 80, rho);
    const E = 2, V = 60, w = 10, roa = (E * (V - w) * (V + w)) / (2 * V);
    return [
      N.check('Hover limit: P = T^1.5/√(2ρA)', h.P, T ** 1.5 / Math.sqrt(2 * rho * A), 1e-9, 'Rankine–Froude momentum theory'),
      N.check('Glauert inflow satisfies vi·√(V² + vi²) = vh²', f.vi * Math.hypot(80, f.vi), T / (2 * rho * A), 1e-8, 'Glauert forward-flight momentum equation'),
      N.check('High-speed induced velocity → T/(2ρAV)', f.vi, T / (2 * rho * A * 80), 2e-3, 'Glauert high-speed limit'),
      N.check('Radius of action E·(V² − w²)/(2V)', roa, E * (V * V - w * w) / (2 * V), 1e-12, 'Out-and-back kinematics'),
    ];
  },
  calibration: { params: [{ key: 'f', min: 0.05, max: 10 }, { key: 'kappa', min: 1, max: 1.6 }, { key: 'bsfc', min: 4e-8, max: 2e-7 }], sweep: 'V', target: 'cruise_fuel_flow_kgs', note: 'Supply measured fuel flow at several cruise speeds at a known mass and altitude.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!o.mission_feasible) out.push({ severity: 'critical', title: 'The mission is not feasible', detail: res.warnings[0], action: o.hover_power_margin_pct <= 0 ? 'Reduce take-off mass or wait for cooler conditions; consider a rolling or in-ground-effect departure only where the procedure allows.' : `Reduce the radius to ${o.radius_max_km.toFixed(0)} km, offload payload for fuel, or plan a refuelling stop.`, basis: 'Fuel with reserves and hover power margin' });
    if (Math.abs(i.V - o.V_br_ms) > 0.12 * o.V_br_ms) out.push({ severity: 'advise', title: 'Cruise speed is away from best range', detail: `Planned ${i.V} m/s against best-range ${o.V_br_ms.toFixed(0)} m/s in this wind.`, action: 'Cruise near best-range speed on fuel-critical legs and at best-endurance speed while searching or holding: both cut fuel and CO₂.', basis: 'Power curve: maximum ground distance per unit energy' });
    if (o.time_on_station_min > 0) out.push({ severity: 'info', title: `${o.time_on_station_min.toFixed(0)} min on station at ${i.radius_km} km`, detail: 'After transit, hover and reserves.', action: 'Each 100 kg of payload traded for fuel adds roughly ' + ((100 / Math.max(o.cruise_fuel_flow_kgs * 0.8, 1e-9)) / 60).toFixed(0) + ' min of search time.', basis: 'Best-endurance fuel flow' });
    out.push({ severity: 'info', title: 'Know the point of no return', detail: `Radius of action ${o.radius_of_action_km.toFixed(0)} km; point of no return ${o.pnr_time_min.toFixed(0)} min after departure.`, action: 'Brief the PNR for every offshore sector and update it with the actual wind.', basis: 'Out-and-back fuel kinematics' });
    return out;
  },
};

const electric = {
  id: 'electric', title: 'Electric, eVTOL and multirotor mission energy', fidelity: 'reduced-order',
  summary: 'Adds up battery energy for hover, transition, climb, cruise and reserve, tracks state of charge through the mission, and shows how range and margin react to headwind and temperature.',
  equations: ['Aircraft energy equations', 'Mission constraint equations', 'Flight path kinematic equations', 'Atmospheric equations'],
  applicable: (c) => (c.prop.type === 'electric' ? true : 'This analysis is for battery-electric aircraft.'),
  inputs: [
    { key: 'mass', label: 'Flight mass', unit: 'kg', default: 2400, min: 0.05, group: 'Mass' }, { key: 'batt_kWh', label: 'Battery energy (nameplate)', unit: 'kWh', default: 160, min: 0.001, group: 'Battery' },
    { key: 'usable_frac', label: 'Usable share of nameplate energy', unit: '-', default: 0.9, min: 0.3, max: 1, group: 'Battery', help: 'State-of-charge window kept for cell life and end-of-life capacity' }, { key: 'pack_Wh_kg', label: 'Pack specific energy', unit: 'Wh/kg', default: 190, min: 30, max: 800, group: 'Battery' },
    { key: 'eta_elec', label: 'Battery-to-shaft efficiency', unit: '-', default: 0.9, min: 0.5, max: 0.99, group: 'Battery' }, { key: 'temp_C', label: 'Battery temperature at start', unit: '°C', default: 20, min: -30, max: 50, group: 'Battery' },
    { key: 'cold_derate', label: 'Capacity loss per K below 20 °C', unit: '1/K', default: 0.006, min: 0, max: 0.03, group: 'Battery', help: 'Illustrative; depends strongly on chemistry and discharge rate' },
    { key: 'has_wing', label: 'Wing-borne cruise', type: 'bool', default: true, group: 'Aerodynamics' }, { key: 'S', label: 'Wing area', unit: 'm²', default: 13.5, min: 0, group: 'Aerodynamics' }, { key: 'CD0', label: 'Zero-lift drag coefficient', unit: '-', default: 0.035, min: 0.001, group: 'Aerodynamics' }, { key: 'k', label: 'Induced drag factor', unit: '-', default: 0.04, min: 0.001, group: 'Aerodynamics' }, { key: 'eta_prop', label: 'Cruise propulsive efficiency', unit: '-', default: 0.83, min: 0.2, max: 0.95, group: 'Aerodynamics' },
    ...ROTOR_IN.map((f) => (f.key === 'R' ? { ...f, default: 1.4, min: 0 } : f.key === 'nR' ? { ...f, default: 6 } : f.key === 'sig' ? { ...f, default: 0.18 } : f.key === 'vt' ? { ...f, default: 161 } : f.key === 'cd0' ? { ...f, default: 0.012 } : f.key === 'f' ? { ...f, default: 0.9 } : f.key === 'eta' ? { ...f, default: 0.92, label: 'Motor-shaft to rotor efficiency' } : f)),
    { key: 'range_km', label: 'Stage length', unit: 'km', default: 120, min: 0, group: 'Mission' }, { key: 'V', label: 'Cruise true airspeed', unit: 'm/s', default: 67, min: 1, group: 'Mission' }, { key: 'alt', label: 'Cruise height above the pad', unit: 'm', default: 600, min: 0, max: 6000, group: 'Mission' },
    { key: 'elev', label: 'Pad elevation', unit: 'm', default: 0, min: -400, max: 4500, group: 'Mission' }, { key: 'dISA', label: 'ISA deviation', unit: 'K', default: 0, min: -40, max: 40, group: 'Mission' }, { key: 'headwind', label: 'Headwind', unit: 'm/s', default: 0, min: -30, max: 30, group: 'Mission' },
    { key: 'hover_min', label: 'Total hover time (take-off and landing)', unit: 'min', default: 2, min: 0, max: 60, group: 'Mission' }, { key: 'trans_s', label: 'Transition time (each way)', unit: 's', default: 30, min: 0, max: 180, group: 'Mission' }, { key: 'roc', label: 'Climb rate', unit: 'm/s', default: 4, min: 0.3, max: 15, group: 'Mission' },
    { key: 'reserve_min', label: 'Reserve time at cruise or loiter power', unit: 'min', default: 20, min: 0, max: 60, group: 'Reserves' }, { key: 'aux_W', label: 'Avionics, payload and thermal-management power', unit: 'W', default: 3000, min: 0, group: 'Mission' },
    { key: 'grid_g_kWh', label: 'Grid carbon intensity', unit: 'g CO₂/kWh', default: 0, min: 0, max: 1200, group: 'Reserves', help: '0 to report in-flight emissions only (zero). Enter the local grid value for charging emissions.' },
  ],
  defaults: (c, up, d) => { const rotorBorne = c.rotor.R_m > 0, wing = c.wing.S_m2 > 0; return { ...rotorDefaults(c, d), R: c.rotor.R_m || 0, mass: c.mass.mtow_kg, batt_kWh: c.systems.batt_kWh, pack_Wh_kg: 0.75 * (BATTERIES[c.systems.batt_chem]?.wh_kg || 250), eta_elec: N.clamp((up.electrical?.motor_eff ?? 0.94) * 0.96, 0.5, 0.99), temp_C: c.site.T_C ?? 15,
    has_wing: wing, S: c.wing.S_m2, CD0: up.cfd?.CD0 ?? c.aero.CD0, k: d.k_induced || undefined, eta_prop: c.prop.eta_prop, range_km: c.mission.range_km, V: c.mission.cruise_V_ms || c.flight.V_ms, alt: Math.max(0, c.mission.cruise_alt_m - (c.site.elev_m || 0)), elev: c.site.elev_m ?? 0, dISA: c.atm.dISA_K,
    hover_min: rotorBorne ? c.mission.hover_min : 0, trans_s: rotorBorne && wing ? 30 : 0, roc: wing ? 4 : 3, reserve_min: Math.min(c.mission.reserve_min, wing ? 30 : 5), aux_W: Math.max(5, 1.2 * c.mass.mtow_kg) }; },
  run(i) {
    const warnings = [], rot = i.R > 0, wing = i.has_wing && i.S > 0;
    const model = (wind, tempC, dT = i.dISA) => {
      const rho0 = isa(i.elev, dT).rho, rhoC = isa(i.elev + i.alt, dT).rho, m = i.mass, r = i, ph = [];
      const Pcr = (V, rho) => (wing ? (dragOf(i, m, V, rho) * V) / i.eta_prop : rot ? rotorPower(r, m, V, rho).P : 0) / i.eta_elec + i.aux_W, Ph = rot ? rotorPower(r, m, 0, rho0).P / i.eta_elec + i.aux_W : 0;
      const add = (name, P, t, dist = 0) => ph.push({ name, P, time: t, E: (P * t) / 3.6e6, dist });
      if (rot) add('Hover (take-off)', Ph, i.hover_min * 30); if (rot && wing) add('Transition out', 0.5 * (Ph + Pcr(i.V, rho0)) * 1.1, i.trans_s, 0.5 * i.V * i.trans_s);
      const tCl = i.alt / i.roc, Vc = wing ? i.V : Math.min(i.V, 0.7 * i.V + 1), Pcl = Pcr(Vc, 0.5 * (rho0 + rhoC)) + (m * G0 * i.roc) / ((wing ? i.eta_prop : i.eta) * i.eta_elec), gs = Math.max(0.5, i.V - wind);
      add('Climb', Pcl, tCl, Math.max(0, Vc - wind) * tCl); const tDe = i.alt / i.roc, Pde = Math.max(i.aux_W, Pcr(Vc, 0.5 * (rho0 + rhoC)) - (0.7 * m * G0 * i.roc) / i.eta_elec), dCD = ph.reduce((s, q) => s + q.dist, 0) + Math.max(0, Vc - wind) * tDe + (rot && wing ? 0.5 * i.V * i.trans_s : 0);
      const dCr = Math.max(0, i.range_km * 1e3 - dCD), PcrV = Pcr(i.V, rhoC); add('Cruise', PcrV, dCr / gs, dCr); add('Descent', Pde, tDe, Math.max(0, Vc - wind) * tDe);
      if (rot && wing) add('Transition in', 0.5 * (Ph + Pcr(i.V, rho0)) * 1.1, i.trans_s, 0.5 * i.V * i.trans_s); if (rot) add('Hover (landing)', Ph, i.hover_min * 30);
      // loiter power: minimum of the cruise power curve
      const Vs = N.linspace(Math.max(2, 0.3 * i.V), 1.3 * i.V, 40), Pl = wing || rot ? Math.min(...Vs.map((V) => Pcr(V, rhoC))) : PcrV, Pres = wing ? PcrV : Pl, Eres = (Pres * i.reserve_min * 60) / 3.6e6;
      const Euse = i.batt_kWh * i.usable_frac * N.clamp(1 - i.cold_derate * Math.max(0, 20 - tempC), 0.3, 1), Em = N.sum(ph.map((q) => q.E)), Efix = Em - ph.find((q) => q.name === 'Cruise').E;
      return { ph, Em, Eres, Euse, margin: Euse - Em - Eres, Ph, PcrV, Pl, gs, maxRange: Math.max(0, ((Euse - Eres - Efix) * 3.6e6 / PcrV) * gs + dCD) };
    };
    const b = model(i.headwind, i.temp_C), feasP = !rot || b.Ph > 0, feas = b.margin >= 0 && feasP, tBlock = N.sum(b.ph.map((q) => q.time));
    const soc = [1]; const tt = [0]; let e = 0, t = 0; for (const q of b.ph) { e += q.E; t += q.time; soc.push(1 - e / i.batt_kWh); tt.push(t / 60); }
    const winds = N.linspace(-10, 20, 16), temps = N.linspace(-20, 40, 13), mw = winds.map((w) => model(w, i.temp_C)), mt = temps.map((T) => model(i.headwind, T, i.dISA + (T - i.temp_C)));
    if (b.margin < 0) warnings.push(`Battery energy is ${(-b.margin).toFixed(2)} kWh short of mission plus reserve; the reachable stage length is ${(b.maxRange / 1e3).toFixed(1)} km.`);
    if (i.temp_C < 5) warnings.push('Cold battery: usable capacity is reduced by the illustrative derating factor; pre-heat the pack where possible.');
    if (!wing && !rot) warnings.push('Neither a wing nor a rotor is defined: no lift model is available and powers are zero.');
    const cRate = Math.max(...b.ph.map((q) => q.P)) / 1e3 / i.batt_kWh; if (cRate > 4) warnings.push(`Peak discharge is ${cRate.toFixed(1)} C: voltage sag and heating at this rate reduce usable energy beyond what is modelled.`);
    return {
      kpis: [
        { key: 'mission_energy_kWh', label: 'Mission energy from the battery', value: b.Em, unit: 'kWh' }, { key: 'block_time_h', label: 'Mission time', value: tBlock / 3600, unit: 'h' }, { key: 'block_fuel_kg', label: 'Fuel burned', value: 0, unit: 'kg', note: 'Battery-electric' },
        { key: 'reserve_energy_kWh', label: 'Reserve energy required', value: b.Eres, unit: 'kWh' }, { key: 'energy_margin_kWh', label: 'Energy margin above reserve', value: b.margin, unit: 'kWh', status: b.margin >= 0 ? 'ok' : 'bad' },
        { key: 'reserve_margin_kg', label: 'Margin as equivalent battery mass', value: (1e3 * b.margin) / i.pack_Wh_kg, unit: 'kg', status: b.margin >= 0 ? 'ok' : 'bad', note: 'Energy margin ÷ pack specific energy' },
        { key: 'mission_feasible', label: 'Mission feasible', value: feas ? 1 : 0, unit: '', status: feas ? 'ok' : 'bad' }, { key: 'soc_end', label: 'State of charge at landing', value: 1 - b.Em / i.batt_kWh, unit: '-', status: 1 - b.Em / i.batt_kWh > 0.2 ? 'ok' : 'warn' },
        { key: 'co2_kg', label: 'CO₂ (in flight, or charging if a grid intensity is set)', value: (i.grid_g_kWh * b.Em) / 0.92 / 1e3, unit: 'kg', note: 'Charging efficiency 92%' },
        { key: 'max_range_km', label: 'Reachable stage length with reserve', value: b.maxRange / 1e3, unit: 'km' }, { key: 'hover_power_W', label: 'Hover electrical power', value: b.Ph, unit: 'W' }, { key: 'cruise_power_W', label: 'Cruise electrical power', value: b.PcrV, unit: 'W' },
        { key: 'hover_endurance_min', label: 'Hover endurance on usable energy', value: b.Ph > 0 ? (b.Euse * 3.6e6) / b.Ph / 60 : 0, unit: 'min' }, { key: 'cruise_endurance_min', label: 'Endurance at minimum-power speed', value: b.Pl > 0 ? (b.Euse * 3.6e6) / b.Pl / 60 : 0, unit: 'min' },
        { key: 'energy_per_km_Wh', label: 'Energy per kilometre', value: i.range_km > 0 ? (1e3 * b.Em) / i.range_km : 0, unit: 'Wh/km' }, { key: 'peak_c_rate', label: 'Peak discharge rate', value: cRate, unit: 'C', status: cRate < 4 ? 'ok' : 'warn' },
      ],
      plots: [
        { type: 'line', title: 'State of charge through the mission', xlabel: 'Time [min]', ylabel: 'State of charge [-]', series: [{ name: 'State of charge', x: tt, y: soc, style: 'line+points' }], annotations: [{ y: 1 - i.usable_frac + b.Eres / i.batt_kWh, label: 'Reserve + unusable' }] },
        { type: 'bar', title: 'Energy by phase', ylabel: 'Energy [kWh]', categories: [...b.ph.map((q) => q.name), 'Reserve'], series: [{ name: 'Energy', y: [...b.ph.map((q) => q.E), b.Eres] }] },
        { type: 'line', title: 'Reachable stage length versus headwind', xlabel: 'Headwind [m/s]', ylabel: 'Stage length [km]', series: [{ name: 'With reserve', x: winds, y: mw.map((q) => q.maxRange / 1e3) }], annotations: [{ y: i.range_km, label: 'Planned stage' }] },
        { type: 'line', title: 'Energy margin versus temperature', xlabel: 'Battery and air temperature [°C]', ylabel: 'Energy margin [kWh]', series: [{ name: 'Margin above reserve', x: temps, y: mt.map((q) => q.margin) }], annotations: [{ y: 0, label: 'No margin' }] },
      ],
      tables: [{ title: 'Mission phases', columns: ['Phase', 'Power [kW]', 'Time [min]', 'Energy [kWh]', 'Distance [km]'], rows: b.ph.map((q) => [q.name, q.P / 1e3, q.time / 60, q.E, q.dist / 1e3]) }],
      outputs: { usable_energy_kWh: b.Euse },
      warnings, models: ['Momentum theory hover and Glauert forward-flight rotor power', 'Parabolic-polar wing-borne cruise', 'Constant battery-to-shaft efficiency', 'Linear cold-temperature capacity derating (illustrative)'],
      assumptions: ['Mass is constant', 'Transition power is 10% above the mean of hover and cruise power', 'Descent recovers no energy; it uses cruise power reduced by 70% of the potential-energy rate', 'Reserve at cruise power for winged aircraft and at minimum-power speed for rotor-borne aircraft', 'No voltage sag, Peukert or ageing model beyond the usable-energy fraction'],
    };
  },
  verify() {
    const i = Object.fromEntries(electric.inputs.map((f) => [f.key, f.default])), r = electric.run(i), o = N.kv(r), sum = N.sum(r.tables[0].rows.map((q) => q[3]));
    const j = { ...i, has_wing: false, S: 0, cd0: 0, kappa: 1, eta: 1, eta_elec: 1, aux_W: 0, elev: 0, dISA: 0, R: 1, nR: 4, mass: 100 }, h = N.kv(electric.run(j)), T = (100 * G0) / 4;
    return [
      N.check('Phase energies sum to the state-of-charge drop', sum, (1 - o.soc_end) * i.batt_kWh, 1e-12, 'Energy conservation'),
      N.check('Ideal hover power of four rotors', h.hover_power_W, (4 * T ** 1.5) / Math.sqrt(2 * isa(0).rho * Math.PI), 1e-9, 'Rankine–Froude momentum theory'),
      N.check('Hover endurance = usable energy / hover power', h.hover_endurance_min, (j.batt_kWh * j.usable_frac * 3.6e6) / h.hover_power_W / 60, 1e-12, 'Definition'),
    ];
  },
  calibration: { params: [{ key: 'eta_elec', min: 0.6, max: 0.99 }, { key: 'CD0', min: 0.01, max: 0.2 }, { key: 'f', min: 0.005, max: 5 }], sweep: 'V', target: 'cruise_power_W', note: 'Supply logged battery power at several cruise speeds.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (!o.mission_feasible) out.push({ severity: 'critical', title: 'Not enough battery energy for mission plus reserve', detail: `Short by ${(-o.energy_margin_kWh).toFixed(2)} kWh.`, action: `Shorten the stage to ${o.max_range_km.toFixed(1)} km, cut hover time (each minute costs ${(o.hover_power_W / 6e4).toFixed(2)} kWh), reduce payload, or fit a larger pack.`, basis: 'Energy balance with reserve' });
    else if (o.energy_margin_kWh < 0.1 * o.mission_energy_kWh) out.push({ severity: 'warn', title: 'Thin energy margin', detail: `${o.energy_margin_kWh.toFixed(2)} kWh above reserve.`, action: 'Check the headwind and temperature plots: a modest headwind or a cold pack removes this margin.', basis: 'Energy balance with reserve' });
    if (o.hover_power_W > 2.5 * o.cruise_power_W && i.hover_min > 0) out.push({ severity: 'advise', title: 'Hover dominates energy use', detail: `Hover draws ${(o.hover_power_W / Math.max(o.cruise_power_W, 1)).toFixed(1)}× cruise power.`, action: 'Shorten hover and transition by procedure design (direct approach paths, prompt transition). It is the cheapest energy saving available and also cuts noise exposure.', basis: 'Phase energy breakdown' });
    if (i.grid_g_kWh === 0) out.push({ severity: 'info', title: 'Charging emissions are not counted', detail: 'In-flight CO₂ is zero; life-cycle emissions depend on the electricity source.', action: 'Enter the grid carbon intensity to see charging emissions, and use low-carbon electricity to realise the benefit.', basis: 'Well-to-wake accounting' });
    return out;
  },
};

const stochastic = {
  id: 'stochastic', title: 'Stochastic fuel planning and daily operations', fidelity: 'numerical',
  summary: 'Samples wind, temperature, payload and delay to get the spread of mission fuel or energy, the chance of eating into reserves, and the discretionary fuel that covers a chosen share of days; also estimates daily utilisation from turnaround times.',
  equations: ['Fuel consumption equations', 'Mission constraint equations', 'Breguet range equation'],
  inputs: [
    { key: 'electric', label: 'Battery-electric aircraft', type: 'bool', default: false, group: 'Nominal mission' },
    { key: 'dist_km', label: 'Stage length', unit: 'km', default: 4500, min: 0.1, group: 'Nominal mission' }, { key: 'V', label: 'Cruise true airspeed', unit: 'm/s', default: 231, min: 1, group: 'Nominal mission' },
    { key: 'rate', label: 'Nominal cruise consumption', unit: 'kg/s or kW', default: 0.72, min: 0, group: 'Nominal mission', help: 'Fuel flow in kg/s, or battery power in kW for electric aircraft' },
    { key: 'fixed', label: 'Fixed fuel or energy (taxi, take-off, climb increment, approach)', unit: 'kg or kWh', default: 1500, min: 0, group: 'Nominal mission' },
    { key: 'hold_rate_frac', label: 'Holding consumption / cruise consumption', unit: '-', default: 0.8, min: 0.2, max: 3, group: 'Nominal mission' },
    { key: 'loaded', label: 'Fuel or usable energy loaded', unit: 'kg or kWh', default: 18800, min: 0, group: 'Nominal mission' }, { key: 'reserve', label: 'Required final reserve', unit: 'kg or kWh', default: 2400, min: 0, group: 'Nominal mission' },
    { key: 'mass', label: 'Nominal flight mass', unit: 'kg', default: 72000, min: 0.05, group: 'Nominal mission' }, { key: 'mass_elast', label: 'Consumption elasticity to mass', unit: '%/%', default: 0.7, min: 0, max: 2, group: 'Sensitivities', help: '≈ 2 × induced share of drag for wings (0.6–1.0); about 1.5 in hover' },
    { key: 'temp_sens', label: 'Consumption change per K', unit: '1/K', default: 0.001, min: -0.01, max: 0.02, group: 'Sensitivities' },
    { key: 'wind_sd', label: 'Headwind standard deviation', unit: 'm/s', default: 12, min: 0, max: 50, group: 'Uncertainty' }, { key: 'wind_mean', label: 'Mean headwind', unit: 'm/s', default: 0, min: -60, max: 60, group: 'Uncertainty', help: 'From live winds aloft when available' },
    { key: 'temp_sd', label: 'Temperature deviation standard deviation', unit: 'K', default: 5, min: 0, max: 25, group: 'Uncertainty' }, { key: 'payload_sd', label: 'Payload standard deviation', unit: 'kg', default: 800, min: 0, group: 'Uncertainty' },
    { key: 'delay_mean_min', label: 'Mean airborne delay (holding, vectors)', unit: 'min', default: 6, min: 0, max: 90, group: 'Uncertainty', help: 'Exponentially distributed' },
    { key: 'coverage', label: 'Share of flights to cover without touching reserves', unit: '-', default: 0.99, min: 0.5, max: 0.9999, group: 'Uncertainty' },
    { key: 'turn_min', label: 'Turnaround time', unit: 'min', default: 45, min: 1, max: 600, group: 'Daily operations' }, { key: 'day_h', label: 'Operating day', unit: 'h', default: 16, min: 1, max: 24, group: 'Daily operations' }, { key: 'days_yr', label: 'Operating days per year', unit: '', default: 340, min: 1, max: 366, group: 'Daily operations' },
    { key: 'ground_min', label: 'Taxi and ground time per flight', unit: 'min', default: 18, min: 0, max: 120, group: 'Daily operations' },
    { key: 'nSamples', label: 'Monte Carlo samples', unit: '', default: 4000, min: 200, max: 200000, step: 1, discrete: true, group: 'Numerics' }, { key: 'seed', label: 'Random seed', unit: '', default: 24, min: 1, max: 1e9, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => {
    const e = estCruise(c, d), dist = (c.wing.S_m2 === 0 && c.prop.type !== 'electric' ? 0.6 : 1) * (c.mission.range_km || 100), t = (dist * 1e3) / e.V, rate = e.elec ? e.P / 1e3 : (up.mission?.cruise_fuel_flow_kgs > 0 ? up.mission.cruise_fuel_flow_kgs : e.ff), per = e.elec ? 1 / 3600 : 1, cruise = rate * t * per, rot = c.wing.S_m2 === 0;
    const block = e.elec ? up.mission?.mission_energy_kWh : up.mission?.block_fuel_kg, fixed = block > cruise ? block - cruise : (rot ? 0.12 : 0.06) * cruise + (e.elec ? c.mission.hover_min / 60 * 3 * rate : 0), loaded = e.elec ? 0.9 * c.systems.batt_kWh : c.mass.fuel_kg, small = c.mass.mtow_kg < 5700;
    return { electric: e.elec, dist_km: dist, V: e.V, rate, fixed, loaded, reserve: Math.min(0.5 * loaded, rate * per * Math.min(c.mission.reserve_min, e.elec ? 10 : 45) * 60 * 0.8 + (e.elec || rot ? 0 : rate * (c.mission.alternate_km * 1e3) / e.V)), mass: c.mass.mtow_kg - 0.4 * c.mass.fuel_kg, mass_elast: rot ? 1.2 : 0.7,
      wind_sd: Math.min(12, 0.12 * e.V), payload_sd: 0.05 * c.mass.payload_kg, delay_mean_min: small ? 2 : 6, hold_rate_frac: rot ? 0.75 : 0.8, turn_min: small ? 20 : c.mass.mtow_kg > 40000 ? 45 : 30, day_h: c.meta.type === 'uav' ? 12 : small ? 10 : 16, days_yr: small ? 250 : 340, ground_min: small ? 6 : 18 };
  },
  run(i) {
    const u = N.rng(i.seed), n = Math.round(i.nSamples), per = i.electric ? 1 / 3600 : 1, unit = i.electric ? 'kWh' : 'kg', warnings = [];
    const burn = (w, dT, dP, delay) => { const rate = i.rate * (1 + (i.mass_elast * dP) / i.mass) * (1 + i.temp_sens * dT), gs = Math.max(0.1 * i.V, i.V - w); return i.fixed + rate * per * ((i.dist_km * 1e3) / gs) + i.hold_rate_frac * rate * per * delay; };
    const nominal = burn(i.wind_mean, 0, 0, 0), F = new Array(n), W = new Array(n);
    for (let s = 0; s < n; s++) { const w = i.wind_mean + i.wind_sd * N.randn(u), dT = i.temp_sd * N.randn(u), dP = i.payload_sd * N.randn(u), dl = -i.delay_mean_min * 60 * Math.log(1 - u() * (1 - 1e-12)); W[s] = w; F[s] = burn(w, dT, dP, dl); }
    const mean = N.mean(F), sd = N.std(F), q = N.quantile(F, i.coverage), avail = i.loaded - i.reserve, pInf = F.filter((f) => f > avail).length / n, extra = q - nominal, hist = N.histogram(F, 30), se = sd / Math.sqrt(n);
    const ns = [200, 400, 800, 1600, 3200, 6400].filter((v) => v <= n), run = ns.map((k) => N.quantile(F.slice(0, k), i.coverage));
    // daily operations: cycles that fit in the operating day
    const tBlock = (i.dist_km * 1e3) / Math.max(1, i.V - i.wind_mean) / 60 + i.ground_min, cyc = Math.max(0, Math.floor((i.day_h * 60 + i.turn_min) / (tBlock + i.turn_min))), fhDay = (cyc * (tBlock - i.ground_min)) / 60;
    if (avail <= 0) warnings.push('The required reserve is at least as large as the load: every flight infringes it.');
    if (pInf > 1 - i.coverage) warnings.push(`${(100 * pInf).toFixed(1)}% of sampled flights land with less than the required reserve; the target is ${(100 * (1 - i.coverage)).toFixed(2)}%.`);
    warnings.push('Uncertainty magnitudes are planning assumptions; replace them with route statistics from your own operation.');
    return {
      kpis: [
        { key: 'nominal_burn', label: `Nominal mission ${i.electric ? 'energy' : 'fuel'}`, value: nominal, unit }, { key: 'mean_burn', label: 'Mean over sampled days', value: mean, unit, note: `± ${(1.96 * se).toPrecision(2)} (95% sampling error)` },
        { key: 'sd_burn', label: 'Standard deviation', value: sd, unit }, { key: 'p95_burn', label: '95th percentile', value: N.quantile(F, 0.95), unit }, { key: 'cover_burn', label: `${(100 * i.coverage).toFixed(1)}th percentile`, value: q, unit },
        { key: 'discretionary', label: `Recommended extra ${i.electric ? 'energy' : 'fuel'} above nominal`, value: Math.max(0, extra), unit, note: 'Covers the chosen share of flights without using the final reserve' },
        { key: 'p_reserve_infringed', label: 'Probability of infringing the reserve with the current load', value: pInf, unit: '-', status: pInf <= 1 - i.coverage ? 'ok' : pInf < 0.05 ? 'warn' : 'bad' },
        { key: 'surplus_at_coverage', label: 'Load above the coverage requirement', value: avail - q, unit, status: avail >= q ? 'ok' : 'bad' },
        { key: 'wind_correlation', label: 'Correlation of burn with headwind', value: N.corr(W, F), unit: '-' },
        { key: 'cycles_per_day', label: 'Flights per operating day', value: cyc, unit: '' }, { key: 'fh_per_day', label: 'Flight hours per day', value: fhDay, unit: 'h' }, { key: 'util_fh_yr', label: 'Annual utilisation', value: fhDay * i.days_yr, unit: 'FH/yr' },
        { key: 'ground_share_pct', label: 'Share of the operating day on the ground', value: 100 * (1 - fhDay / i.day_h), unit: '%' },
      ],
      plots: [
        { type: 'bar', title: `Distribution of mission ${i.electric ? 'energy' : 'fuel'}`, ylabel: 'Flights [-]', categories: hist.centers.map((v) => v.toPrecision(4)), series: [{ name: 'Sampled flights', y: hist.counts }] },
        { type: 'line', title: 'Exceedance curve', xlabel: `Mission ${i.electric ? 'energy [kWh]' : 'fuel [kg]'}`, ylabel: 'Probability of exceeding [-]', ylog: true, series: [{ name: 'Sampled', x: N.range(60, (k) => N.quantile(F, k / 60 * 0.999)), y: N.range(60, (k) => Math.max(1e-4, 1 - (k / 60) * 0.999)) }], annotations: [{ x: avail, label: 'Load − reserve' }, { x: nominal, label: 'Nominal' }] },
        { type: 'line', title: 'Convergence of the coverage percentile', xlabel: 'Samples [-]', ylabel: `Percentile [${unit}]`, xlog: true, series: [{ name: 'Running estimate', x: ns, y: run, style: 'line+points' }] },
      ],
      warnings, models: ['Monte Carlo sampling (seeded) of wind, temperature, payload and delay', 'Mission burn = fixed part + consumption × air time + holding', 'Daily cycle count from block and turnaround time'],
      assumptions: ['Wind, temperature and payload are independent normal variables; airborne delay is exponential', 'Consumption scales linearly with mass and temperature deviations', 'Ground speed is floored at 10% of airspeed', 'Daily operations repeat one stage with a fixed turnaround'],
    };
  },
  convergence: { param: 'nSamples', label: 'Monte Carlo samples', levels: [500, 1000, 2000, 4000, 8000], metric: 'mean_burn', hOf: (n) => 1 / Math.sqrt(n) },
  verify() {
    const b = { ...Object.fromEntries(stochastic.inputs.map((f) => [f.key, f.default])), wind_sd: 0, temp_sd: 0, payload_sd: 0, wind_mean: 0 };
    const d = N.kv(stochastic.run({ ...b, delay_mean_min: 0 })), e = N.kv(stochastic.run({ ...b, delay_mean_min: 10, nSamples: 40000 })), nom = b.fixed + (b.rate * b.dist_km * 1e3) / b.V;
    return [
      N.check('Zero uncertainty reproduces the deterministic burn', d.mean_burn, nom, 1e-12, 'Degenerate distribution'),
      N.check('Zero uncertainty has zero spread', d.sd_burn / d.mean_burn + 1, 1, 1e-9, 'Degenerate distribution'),
      N.check('Mean with exponential delay = nominal + holding rate × mean delay', e.mean_burn, nom + 0.8 * b.rate * 600, 2e-3, 'Expectation of an exponential variable; sampling error ≈ 0.01%'),
      N.check('Standard deviation with exponential delay = holding rate × mean delay', e.sd_burn, 0.8 * b.rate * 600, 0.02, 'Exponential distribution: σ = mean'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], unit = i.electric ? 'kWh' : 'kg';
    if (o.surplus_at_coverage < 0) out.push({ severity: 'warn', title: 'Current load does not cover the chosen share of days', detail: `Short by ${(-o.surplus_at_coverage).toFixed(0)} ${unit} at ${(100 * i.coverage).toFixed(1)}% coverage; ${(100 * o.p_reserve_infringed).toFixed(1)}% of flights would use final reserve.`, action: i.electric ? 'Shorten the stage, raise the minimum departure state of charge, or accept a lower coverage with a firm diversion plan.' : 'Uplift the recommended discretionary fuel on this route, or plan an en-route alternate that allows a reduced contingency.', basis: 'Monte Carlo percentile of mission burn against load less reserve' });
    else if (!i.electric && o.surplus_at_coverage > 0.05 * i.loaded) out.push({ severity: 'advise', title: 'Statistical fuel planning allows a lower uplift', detail: `${o.surplus_at_coverage.toFixed(0)} kg more than the ${(100 * i.coverage).toFixed(1)}% requirement is carried.`, action: `Carrying it costs roughly 2.5–4% of its mass in fuel per flight hour. Trimming the uplift to nominal + ${o.discretionary.toFixed(0)} kg keeps the same protection and saves fuel and CO₂ every flight.`, basis: 'Cost of weight; statistical contingency fuel' });
    if (Math.abs(o.wind_correlation) > 0.7) out.push({ severity: 'info', title: 'Wind drives the spread', detail: `Correlation with headwind ${o.wind_correlation.toFixed(2)}.`, action: 'Use day-of-operation winds from the live-data connector: it removes most of the uncertainty and the extra fuel that goes with it.', basis: 'Sample correlation' });
    if (o.ground_share_pct > 60) out.push({ severity: 'advise', title: 'The aircraft spends most of the day on the ground', detail: `${o.cycles_per_day} flights and ${o.fh_per_day.toFixed(1)} flight hours in a ${i.day_h} h day.`, action: 'Shorter turnarounds or a longer operating day raise utilisation and spread fixed ownership cost over more hours (see Suite 26).', basis: 'Daily cycle model' });
    return out;
  },
};

export default {
  id: 'mission', n: 24,
  tagline: 'Can the aircraft fly this mission with legal reserves, how much fuel or energy and time does it take, and what would make it cheaper or more robust.',
  analyses: [profile, route, trade, climbopt, heli, electric, stochastic],
  consumes: [
    { from: 'cfd', keys: ['CD0', 'k_induced', 'CLmax'], why: 'Drag polar and lift limit' },
    { from: 'propulsion', keys: ['thrust_static_N', 'tsfc_kg_Ns'], why: 'Installed thrust and fuel consumption' },
    { from: 'propeller', keys: ['eta_prop'], why: 'Propulsive efficiency' },
    { from: 'performance', keys: ['fuel_flow_cruise_kgs'], why: 'Cruise fuel flow for route fuel' },
    { from: 'electrical', keys: ['motor_eff'], why: 'Drive-train efficiency for electric missions' },
  ],
  provides: [
    { key: 'block_fuel_kg', label: 'Block fuel', unit: 'kg' }, { key: 'block_time_h', label: 'Block time', unit: 'h' }, { key: 'mission_energy_kWh', label: 'Mission energy', unit: 'kWh' },
    { key: 'reserve_margin_kg', label: 'Reserve margin', unit: 'kg' }, { key: 'co2_kg', label: 'Mission CO₂', unit: 'kg' }, { key: 'mission_feasible', label: 'Mission feasible', unit: '' },
  ],
  handoff: [
    { model: 'Six-degree-of-freedom mission simulation', why: 'Mission fuel and time are governed by point-mass energy; attitude dynamics belong to Suites 4 and 16', tool: 'Suite 4 / 16, or a full flight simulator' },
    { model: 'Trajectory optimisation by Pontryagin’s principle or direct collocation with 4-D weather', why: 'Needs gridded wind and temperature fields and an optimal-control solver; the energy-state method and speed/altitude maps are provided instead', tool: 'Flight-planning systems, optimal-control software' },
    { model: 'Air-traffic network and fleet scheduling optimisation', why: 'Combinatorial problems over many aircraft, crews and slots', tool: 'Airline operations-research suites' },
    { model: 'Certified performance and flight-planning data', why: 'Operational planning must use the approved flight manual and operator fuel policy', tool: 'Manufacturer performance software, operator flight-planning system' },
    { model: 'Hybrid-electric energy management', why: 'Needs engine and battery maps with a power-split strategy; only fuel-only and battery-only missions are solved', tool: 'Suite 19 with a dedicated energy-management optimiser' },
  ],
};
