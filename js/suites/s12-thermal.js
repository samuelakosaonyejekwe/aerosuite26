// Suite 12 — Thermal Engineering and Heat Transfer.
// Transient multilayer conduction, 2-D board/plate conduction, aerodynamic heating, heat exchangers,
// lumped thermal networks for electronics and batteries, thermal stress and radiation exchange.

import * as N from '../core/numerics.js';
import { isa, kAir, sutherland, CP_AIR, PR_AIR, R_AIR, G0 } from '../core/atmosphere.js';
import { METALS, FLUIDS } from '../data/materials.js';

// ---- shared helpers --------------------------------------------------------------------------
const SIG = 5.670374419e-8; // Stefan–Boltzmann constant [W/m²/K⁴]
const num = (key, label, unit, def, min, max, group, help, x) => ({ key, label, unit, default: def, min, max, group, ...(help ? { help } : {}), ...x });
const sel = (key, label, options, def, group, help) => ({ key, label, type: 'select', options, default: def, group, ...(help ? { help } : {}) });
const kp = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const thin = (a, n = 300) => { if (a.length <= n) return a; const s = (a.length - 1) / (n - 1); return N.range(n, (j) => a[Math.round(j * s)]); };
const METAL_NAMES = Object.keys(METALS);
const metal = (name) => METALS[name] || METALS['Al 2024-T3'];
const siteT = (c) => (c.site?.T_C ?? 15) + 273.15;
/** Turbulent flat-plate local heat-transfer coefficient [W/m²/K] at distance x (Colburn analogy). */
const hTurb = (V, x, a) => { const k = kAir(a.T), Re = (a.rho * V * x) / a.mu; return (0.0296 * Re ** 0.8 * PR_AIR ** (1 / 3) * k) / x; };
const cruise = (c) => { const a = isa(c.mission.cruise_alt_m || c.atm.alt_m || 0, c.atm.dISA_K || 0), V = c.mission.cruise_V_ms || c.flight.V_ms || 50; return { a, V, M: V / a.a }; };
const recovery = (T, M, turbulent = true) => T * (1 + (turbulent ? PR_AIR ** (1 / 3) : Math.sqrt(PR_AIR)) * 0.2 * M * M);

// ---- 1. transient multilayer wall ------------------------------------------------------------
const SCHEMES = ['Crank–Nicolson', 'Backward Euler'];
/** Finite-volume θ-scheme for a 1-D multilayer wall with convection + radiation on both faces. */
function wallSolve(i) {
  const m1 = metal(i.mat1);
  const L = [{ t: i.t1_mm / 1e3, k: m1.k, rc: m1.rho * m1.cp, g: i.qgen, contact: true }, { t: i.t2_mm / 1e3, k: i.k2, rc: i.rho2 * i.cp2, g: 0 }, { t: i.t3_mm / 1e3, k: i.k3, rc: i.rho3 * i.cp3, g: 0 }].filter((l) => l.t > 0);
  if (!L.length) L.push({ t: 1e-3, k: m1.k, rc: m1.rho * m1.cp, g: i.qgen });
  const Lt = N.sum(L.map((l) => l.t)), nReq = Math.max(Math.round(i.nNodes), 2 * L.length), dx = [], k = [], rc = [], g = [], lay = [];
  L.forEach((l, j) => { const n = Math.max(2, Math.round((nReq * l.t) / Lt)); for (let q = 0; q < n; q++) { dx.push(l.t / n); k.push(l.k); rc.push(l.rc); g.push(l.g); lay.push(j); } });
  const n = dx.length, xc = []; let x0 = 0;
  for (let j = 0; j < n; j++) { xc.push(x0 + dx[j] / 2); x0 += dx[j]; }
  const G = N.range(n - 1, (j) => 1 / (dx[j] / (2 * k[j]) + dx[j + 1] / (2 * k[j + 1]) + (lay[j] !== lay[j + 1] && L[lay[j]].contact ? i.Rc : 0)));
  const C = dx.map((d, j) => rc[j] * d), th = i.scheme === SCHEMES[1] ? 1 : 0.5, nt = Math.max(1, Math.round(i.nSteps)), dt = i.t_end / nt, gen = N.sum(g.map((v, j) => v * dx[j]));
  // massless surface node: conduction to the first cell balances convection + linearised radiation + absorbed flux
  const side = (Tc, a, h, Tf, eps, Tr, q, Ts) => {
    const hr = eps * SIG * (Ts * Ts + Tr * Tr) * (Ts + Tr), b = h + hr;
    if (b < 1e-9) return { U: 0, Teq: 0, q0: q, Ts: Tc + q / a, flux: q };
    const Teq = (h * Tf + hr * Tr + q) / b, U = (a * b) / (a + b);
    return { U, Teq, q0: 0, Ts: (a * Tc + b * Teq) / (a + b), flux: U * (Teq - Tc) };
  };
  const aO = (2 * k[0]) / dx[0], aI = (2 * k[n - 1]) / dx[n - 1];
  const sO = (Tc, Ts) => side(Tc, aO, i.h_o, i.T_o, i.eps_o, i.T_sky, i.q_abs, Ts), sI = (Tc, Ts) => side(Tc, aI, i.h_i, i.T_i, i.eps_i, i.T_i, 0, Ts);
  let T = new Array(n).fill(i.T_init), bO = sO(T[0], i.T_init), bI = sI(T[n - 1], i.T_init);
  for (let q = 0; q < 3; q++) { bO = sO(T[0], bO.Ts); bI = sI(T[n - 1], bI.Ts); }
  const E0 = N.dot(C, T), time = [0], TsO = [bO.Ts], TsI = [bI.Ts], qO = [bO.flux], qI = [bI.flux], snaps = [{ t: 0, T: [bO.Ts, ...T, bI.Ts] }], snapAt = [0.02, 0.1, 0.3, 0.6, 1].map((f) => Math.max(1, Math.round(f * nt)));
  let Ein = 0, Tmax = Math.max(bO.Ts, bI.Ts, i.T_init), Tmin = Math.min(bO.Ts, bI.Ts, i.T_init), qmax = Math.abs(bO.flux);
  const sub = new Array(n), dia = new Array(n), sup = new Array(n), rhs = new Array(n);
  for (let s = 1; s <= nt; s++) {
    const Rold = T.map((Tj, j) => (j > 0 ? G[j - 1] * (T[j - 1] - Tj) : bO.flux) + (j < n - 1 ? G[j] * (T[j + 1] - Tj) : bI.flux) + g[j] * dx[j]);
    if (n === 1) Rold[0] = bO.flux + bI.flux + g[0] * dx[0];
    let cO = bO, cI = bI, Tn = T, fO = 0, fI = 0;
    for (let it = 0; it < 3; it++) {
      for (let j = 0; j < n; j++) {
        const gl = j > 0 ? G[j - 1] : 0, gr = j < n - 1 ? G[j] : 0;
        sub[j] = -th * gl; sup[j] = -th * gr; dia[j] = C[j] / dt + th * (gl + gr); rhs[j] = (C[j] / dt) * T[j] + (1 - th) * Rold[j] + th * g[j] * dx[j];
      }
      dia[0] += th * cO.U; rhs[0] += th * (cO.U * cO.Teq + cO.q0); dia[n - 1] += th * cI.U; rhs[n - 1] += th * (cI.U * cI.Teq + cI.q0);
      Tn = N.solveTridiag(sub, dia, sup, rhs);
      fO = cO.U * (cO.Teq - Tn[0]) + cO.q0; fI = cI.U * (cI.Teq - Tn[n - 1]) + cI.q0;
      cO = sO(Tn[0], cO.Ts); cI = sI(Tn[n - 1], cI.Ts);
    }
    Ein += dt * (th * (fO + fI + gen) + (1 - th) * (bO.flux + bI.flux + gen));
    T = Tn; bO = cO; bI = cI;
    time.push(s * dt); TsO.push(bO.Ts); TsI.push(bI.Ts); qO.push(bO.flux); qI.push(bI.flux);
    Tmax = Math.max(Tmax, bO.Ts, bI.Ts, N.amax(T)); Tmin = Math.min(Tmin, bO.Ts, bI.Ts, N.amin(T)); qmax = Math.max(qmax, Math.abs(bO.flux));
    if (snapAt.includes(s)) snaps.push({ t: s * dt, T: [bO.Ts, ...T, bI.Ts] });
  }
  const dE = N.dot(C, T) - E0, Rcond = N.sum(L.map((l) => l.t / l.k)) + (L.length > 1 && L[0].contact ? i.Rc : 0);
  return { x: [0, ...xc, Lt], xc, T, time, TsO, TsI, qO, qI, snaps, Tmax, Tmin, qmax, Lt, Rcond, Ctot: N.sum(C), eErr: Math.abs(dE - Ein) / Math.max(Math.abs(dE), Math.abs(Ein), 1e-9), n };
}

const wall = {
  id: 'wall', title: 'Transient conduction through a multilayer skin or wall', fidelity: 'numerical',
  summary: 'Temperature history through up to three layers (skin, insulation, liner) with convection, radiation and absorbed heat flux on both faces, for example a hot-soaked airframe cooling down in cruise.',
  equations: ['Fourier heat conduction equation', 'Transient heat diffusion equation', 'Newton’s law of cooling', 'Stefan–Boltzmann radiation law', 'Thermal resistance relations', 'Conservation of thermal energy equation'],
  inputs: [
    sel('mat1', 'Layer 1 (outer skin) metal', METAL_NAMES, 'Al 2024-T3', 'Layers'),
    num('t1_mm', 'Layer 1 thickness', 'mm', 2, 0, 200, 'Layers', 'Set 0 to omit the metal skin'),
    num('qgen', 'Internal heat generation in layer 1', 'W/m³', 0, 0, 1e9, 'Layers', 'e.g. embedded heater mat or ohmic heating'),
    num('Rc', 'Contact resistance between layers 1 and 2', 'm²K/W', 0, 0, 1, 'Layers', 'Typical 1e-4 (clamped metal joint) to 1e-3 (dry, light pressure)'),
    num('t2_mm', 'Layer 2 (insulation) thickness', 'mm', 50, 0, 500, 'Layers', 'Set 0 to omit'),
    num('k2', 'Layer 2 conductivity', 'W/m/K', 0.04, 0.005, 500, 'Layers', 'Glass-fibre blanket ≈ 0.035–0.045'),
    num('rho2', 'Layer 2 density', 'kg/m³', 10, 1, 20000, 'Layers'), num('cp2', 'Layer 2 specific heat', 'J/kg/K', 840, 100, 5000, 'Layers'),
    num('t3_mm', 'Layer 3 (liner) thickness', 'mm', 0, 0, 200, 'Layers', 'Set 0 to omit'),
    num('k3', 'Layer 3 conductivity', 'W/m/K', 0.3, 0.005, 500, 'Layers'), num('rho3', 'Layer 3 density', 'kg/m³', 1500, 1, 20000, 'Layers'), num('cp3', 'Layer 3 specific heat', 'J/kg/K', 1100, 100, 5000, 'Layers'),
    num('h_o', 'Outer convection coefficient', 'W/m²/K', 60, 0, 1e5, 'Outer boundary', 'Cruise boundary layer 30–150; still air 5–10'),
    num('T_o', 'Outer fluid (recovery) temperature', 'K', 240, 50, 3000, 'Outer boundary', 'Use the recovery temperature in flight'),
    num('eps_o', 'Outer emissivity', '-', 0.85, 0, 1, 'Outer boundary', 'Paint 0.85–0.95, bare polished aluminium 0.05–0.1'),
    num('T_sky', 'Outer radiative environment temperature', 'K', 230, 3, 3000, 'Outer boundary'),
    num('q_abs', 'Absorbed external heat flux', 'W/m²', 0, -1e6, 1e7, 'Outer boundary', 'Absorbed solar (α·G) or prescribed aerodynamic heating'),
    num('h_i', 'Inner convection coefficient', 'W/m²/K', 5, 0, 1e5, 'Inner boundary', 'Set 0 with zero emissivity for an adiabatic back face'),
    num('T_i', 'Inner air temperature', 'K', 295, 50, 3000, 'Inner boundary'), num('eps_i', 'Inner emissivity', '-', 0, 0, 1, 'Inner boundary'),
    num('T_init', 'Initial temperature (uniform)', 'K', 303, 50, 3000, 'Initial state', 'Filled from the site temperature'),
    num('T_limit', 'Allowable material temperature', 'K', 393, 100, 3000, 'Limits', 'Typical long-term limit ≈ 390–420 K for aluminium alloys, ≈ 350–450 K for epoxy composites'),
    num('t_end', 'Simulated time', 's', 1800, 1, 1e7, 'Numerics'),
    sel('scheme', 'Time scheme', SCHEMES, SCHEMES[0], 'Numerics', 'Crank–Nicolson is second order; Backward Euler is first order but never oscillates'),
    num('nNodes', 'Finite-volume cells through the wall', '', 40, 4, 2000, 'Numerics', '', { step: 1, discrete: true }),
    num('nSteps', 'Time steps', '', 400, 4, 100000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => {
    const { a, V, M } = cruise(c), cabin = (c.mission.pax || 0) > 0;
    return { mat1: METALS[c.struct.material] ? c.struct.material : undefined, t1_mm: c.struct.t_skin_mm, t2_mm: cabin ? 50 : 0, T_o: recovery(a.T, M), h_o: hTurb(V, Math.max(0.3, 0.3 * (c.fuselage.len_m || 1)), a), T_sky: a.T, T_init: siteT(c), T_i: cabin ? 295 : siteT(c), h_i: cabin ? 5 : 3 };
  },
  run(i) {
    const r = wallSolve(i), U = 1 / ((i.h_o > 0 ? 1 / i.h_o : Infinity) + r.Rcond + (i.h_i > 0 ? 1 / i.h_i : Infinity)), m1 = metal(i.mat1);
    const Bi = (Math.max(i.h_o, i.h_i) * r.Lt) / (r.Lt / r.Rcond || 1), margin = i.T_limit - r.Tmax, warnings = [];
    const tau = r.Ctot / Math.max(1e-12, (i.h_o + 4 * i.eps_o * SIG * i.T_o ** 3) + i.h_i), last = r.time.length - 1;
    if (r.eErr > 1e-6) warnings.push('Energy balance closure is poorer than 1e-6: increase the number of time steps (radiation is linearised within each step).');
    if (i.scheme === SCHEMES[0] && (i.t_end / i.nSteps) * (Math.max(i.h_o, i.h_i, 1) / (r.Ctot / r.n)) > 20) warnings.push('The time step is large compared with the surface-cell time constant: Crank–Nicolson may show decaying oscillations. Use more steps or Backward Euler.');
    if (r.Tmax > 0.6 * (m1.JC?.Tm || 900) && i.t1_mm > 0) warnings.push('Skin temperature approaches a range where conductivity, specific heat and strength vary strongly; constant room-temperature properties are used.');
    if (Math.abs(r.TsI[last] - r.TsI[Math.round(0.9 * last)]) > 0.02 * Math.abs(r.TsI[last] - i.T_init) + 0.05) warnings.push('The wall has not reached a steady state by the end of the simulated time.');
    return {
      kpis: [
        kp('T_max_K', 'Peak temperature in the wall', r.Tmax, 'K', margin > 20 ? 'ok' : margin > 0 ? 'warn' : 'bad', `Limit ${i.T_limit.toFixed(0)} K`),
        kp('T_min_K', 'Lowest temperature in the wall', r.Tmin, 'K'),
        kp('T_outer_K', 'Outer surface temperature at end', r.TsO[last], 'K'),
        kp('T_inner_K', 'Inner surface temperature at end', r.TsI[last], 'K'),
        kp('q_max_Wm2', 'Peak outer-surface heat flux', r.qmax, 'W/m²'),
        kp('q_through_Wm2', 'Heat flux into the inner space at end', -r.qI[last], 'W/m²', undefined, 'Positive = heat entering the cabin or bay'),
        kp('U_overall', 'Steady overall heat-transfer coefficient (convection only)', Number.isFinite(U) ? U : 0, 'W/m²/K'),
        kp('thermal_margin_K', 'Margin to the allowable temperature', margin, 'K', margin > 20 ? 'ok' : margin > 0 ? 'warn' : 'bad'),
        kp('time_constant_s', 'Lumped time constant C/(h_o+h_r+h_i)', tau, 's'),
        kp('biot', 'Biot number of the wall', Bi, '-', undefined, 'Below 0.1 the wall is nearly isothermal'),
        kp('energy_error', 'Energy balance closure error', r.eErr, '-', r.eErr < 1e-6 ? 'ok' : 'warn'),
      ],
      plots: [
        { type: 'line', title: 'Temperature profiles through the wall', xlabel: 'Distance from outer surface [mm]', ylabel: 'Temperature [K]', series: r.snaps.map((s) => ({ name: `t = ${s.t.toFixed(0)} s`, x: r.x.map((v) => v * 1e3), y: s.T })) },
        { type: 'line', title: 'Surface temperature histories', xlabel: 'Time [s]', ylabel: 'Temperature [K]', series: [{ name: 'Outer surface', x: thin(r.time), y: thin(r.TsO) }, { name: 'Inner surface', x: thin(r.time), y: thin(r.TsI) }], annotations: [{ y: i.T_limit, label: 'Allowable' }] },
        { type: 'line', title: 'Surface heat flux into the wall', xlabel: 'Time [s]', ylabel: 'Heat flux [W/m²]', series: [{ name: 'Outer face', x: thin(r.time), y: thin(r.qO) }, { name: 'Inner face', x: thin(r.time), y: thin(r.qI) }] },
      ],
      warnings, models: [`Finite-volume conduction, ${i.scheme} in time`, 'Harmonic-mean interface conductance with contact resistance', 'Radiation linearised by Picard iteration within each step'],
      assumptions: ['One-dimensional heat flow normal to the wall', 'Constant material properties', 'Grey diffuse surfaces exchanging with a single environment temperature per face', 'Boundary conditions constant over the simulated time'],
    };
  },
  convergence: { param: 'nNodes', label: 'Cells through the wall', levels: [8, 16, 32, 64], metric: 'T_inner_K' },
  verify() {
    const base = { mat1: 'Al 2024-T3', t1_mm: 0, qgen: 0, Rc: 0, t2_mm: 100, k2: 1, rho2: 1000, cp2: 1000, t3_mm: 0, k3: 1, rho3: 1, cp3: 1, h_o: 50, T_o: 300, eps_o: 0, T_sky: 300, q_abs: 0, h_i: 50, T_i: 300, eps_i: 0, T_init: 400, T_limit: 500, t_end: 1000, scheme: SCHEMES[0], nNodes: 80, nSteps: 400 };
    // (1) symmetric slab, half-thickness L, Bi = 2.5, Fo = 0.4: exact eigenfunction series for the mid-plane temperature
    const r = wallSolve(base), Lh = 0.05, Bi = (50 * Lh) / 1, Fo = (1e-6 * 1000) / Lh ** 2; let th = 0;
    for (let m = 0; m < 30; m++) { const z = N.brent((x) => x * Math.tan(x) - Bi, m * Math.PI + 1e-9, m * Math.PI + Math.PI / 2 - 1e-9, 1e-13); th += ((4 * Math.sin(z)) / (2 * z + Math.sin(2 * z))) * Math.exp(-z * z * Fo); }
    const mid = 0.5 * (r.T[39] + r.T[40]);
    // (2) steady three-layer wall with contact resistance: q = ΔT / ΣR
    const s = wallSolve({ ...base, t1_mm: 5, Rc: 2e-3, t2_mm: 20, k2: 0.5, t3_mm: 10, k3: 2, h_o: 100, T_o: 400, h_i: 10, T_i: 300, T_init: 350, t_end: 4e6, scheme: SCHEMES[1], nNodes: 30, nSteps: 200 });
    const qRef = 100 / (1 / 100 + 0.005 / 121 + 2e-3 + 0.02 / 0.5 + 0.01 / 2 + 1 / 10);
    // (3) radiative equilibrium of an insulated-back plate under absorbed flux: T = (q/εσ + Tsky⁴)^¼
    const e = wallSolve({ ...base, t2_mm: 2, h_o: 0, eps_o: 0.8, T_sky: 250, q_abs: 800, h_i: 0, T_init: 300, t_end: 2e6, scheme: SCHEMES[1], nNodes: 8, nSteps: 400 });
    return [
      N.check('Slab mid-plane θ at Bi = 2.5, Fo = 0.4', (mid - 300) / 100, th, 3e-3, 'Exact series Σ Cn·exp(−ζn²Fo), ζn·tan ζn = Bi (Incropera, ch. 5)'),
      N.check('Steady multilayer heat flux = ΔT/ΣR', s.qO[s.qO.length - 1], qRef, 1e-4, 'Series thermal resistances'),
      N.check('Radiative-equilibrium temperature', e.TsO[e.TsO.length - 1], (800 / (0.8 * SIG) + 250 ** 4) ** 0.25, 1e-5, 'Stefan–Boltzmann balance'),
      N.check('Discrete energy conservation', r.eErr, 0, 1e-9, 'Stored energy change equals integrated boundary heat'),
    ];
  },
  calibration: { params: [{ key: 'h_o', min: 1, max: 2000 }, { key: 'k2', min: 0.005, max: 5 }, { key: 'Rc', min: 0, max: 0.05 }], sweep: 't_end', target: 'T_inner_K', note: 'Supply thermocouple temperatures of the inner surface against time from a thermal-chamber or flight cold-soak test.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.thermal_margin_K < 0) out.push({ severity: 'critical', title: 'Wall exceeds its allowable temperature', detail: `Peak ${o.T_max_K.toFixed(0)} K against a limit of ${i.T_limit.toFixed(0)} K.`, action: 'Add insulation or a thermal barrier, raise surface emissivity, lower absorbed flux (lighter paint) or select a higher-temperature material (titanium, CMC).', basis: 'Material temperature limit' });
    if (o.q_through_Wm2 && Math.abs(o.q_through_Wm2) > 60 && i.t2_mm > 0) out.push({ severity: 'advise', title: 'Significant heat leak through the wall', detail: `${Math.abs(o.q_through_Wm2).toFixed(0)} W/m² ${o.q_through_Wm2 < 0 ? 'leaves' : 'enters'} the inner space; every kW the environmental control system must make up costs bleed or electrical power and therefore fuel.`, action: 'Increase insulation thickness or reduce thermal bridges; carry the heat leak to Suite 20 (cabin load budget).', basis: 'Steady heat balance' });
    if (i.T_i > 285 && o.T_inner_K < 283) out.push({ severity: 'advise', title: 'Cold inner surface: condensation risk', detail: `Inner surface settles at ${o.T_inner_K.toFixed(0)} K, likely below the cabin dew point.`, action: 'Check the dew point in Suite 20 and provide drainage or a vapour barrier; trapped moisture adds mass and promotes corrosion.', basis: 'Surface temperature versus dew point' });
    return out;
  },
};

// ---- 2. 2-D steady conduction on a board / plate ---------------------------------------------
const EDGES = ['All four edges held at sink temperature', 'Left/right edges held (card guides), others insulated', 'All edges insulated (surface convection only)'];
function plateSolve(i) {
  const nx = Math.max(4, Math.round(i.nx)), ny = Math.max(4, Math.round((nx * i.Ly) / i.Lx)), dx = i.Lx / nx, dy = i.Ly / ny, kt = (i.k * i.t_mm) / 1e3;
  const dirX = i.edges !== EDGES[2], dirY = i.edges === EDGES[0]; let hf = i.h * i.faces, clampH = false;
  if (!dirX && hf < 0.05) { hf = 0.05; clampH = true; }
  const xs = N.linspace(0, i.Lx, nx + 1), ys = N.linspace(0, i.Ly, ny + 1), cvx = (j, n, d) => [Math.max(0, (j - 0.5) * d), Math.min(n * d, (j + 0.5) * d)];
  const q = N.zeros(ny + 1, nx + 1), bg = i.P_bg / (i.Lx * i.Ly);
  for (let j = 0; j <= ny; j++) for (let ii = 0; ii <= nx; ii++) q[j][ii] = bg;
  for (const [P, fx, fy, s] of [[i.P1, i.x1, i.y1, i.s1], [i.P2, i.x2, i.y2, i.s2], [i.P3, i.x3, i.y3, i.s3]]) {
    if (!(P > 0)) continue;
    const hx = Math.min(s / 2, i.Lx / 2), hy = Math.min(s / 2, i.Ly / 2), cx = N.clamp(fx * i.Lx, hx, i.Lx - hx), cy = N.clamp(fy * i.Ly, hy, i.Ly - hy), qf = P / (4 * hx * hy);
    for (let j = 0; j <= ny; j++) { const [y0, y1] = cvx(j, ny, dy), oy = Math.min(y1, cy + hy) - Math.max(y0, cy - hy); if (oy <= 0) continue;
      for (let ii = 0; ii <= nx; ii++) { const [x0, x1] = cvx(ii, nx, dx), ox = Math.min(x1, cx + hx) - Math.max(x0, cx - hx); if (ox > 0) q[j][ii] += (qf * ox * oy) / ((x1 - x0) * (y1 - y0)); } }
  }
  const T = N.range(ny + 1, () => new Array(nx + 1).fill(dirX ? i.T_edge : i.T_air)), ax = kt / (dx * dx), ay = kt / (dy * dy), den = 2 * ax + 2 * ay + hf;
  const i0 = dirX ? 1 : 0, i1 = dirX ? nx - 1 : nx, j0 = dirY ? 1 : 0, j1 = dirY ? ny - 1 : ny, om = 2 / (1 + Math.sin(Math.PI / (Math.max(nx, ny) * (dirY ? 1 : 2))));
  const maxIt = 30000; let it = 0, dmax = Infinity, scale = 1;
  for (; it < maxIt && dmax > 1e-9 * scale; it++) {
    dmax = 0; let lo = Infinity, hi = -Infinity;
    for (let j = j0; j <= j1; j++) { const row = T[j], S = j > 0 ? T[j - 1] : T[1], Nn = j < ny ? T[j + 1] : T[ny - 1], qj = q[j];
      for (let ii = i0; ii <= i1; ii++) {
        const W = ii > 0 ? row[ii - 1] : row[1], E = ii < nx ? row[ii + 1] : row[nx - 1], d = om * ((ax * (W + E) + ay * (S[ii] + Nn[ii]) + hf * i.T_air + qj[ii]) / den - row[ii]);
        row[ii] += d; if (Math.abs(d) > dmax) dmax = Math.abs(d); if (row[ii] < lo) lo = row[ii]; if (row[ii] > hi) hi = row[ii];
      } }
    scale = Math.max(1, hi - lo);
  }
  // discrete heat balance over the free nodes
  let Pin = 0, Qc = 0, Qe = 0, Tmax = -Infinity, Tsum = 0, At = 0, jm = 0, im = 0;
  for (let j = j0; j <= j1; j++) for (let ii = i0; ii <= i1; ii++) {
    const [x0, x1] = cvx(ii, nx, dx), [y0, y1] = cvx(j, ny, dy), A = (x1 - x0) * (y1 - y0), Tp = T[j][ii];
    Pin += q[j][ii] * A; Qc += hf * (Tp - i.T_air) * A; Tsum += Tp * A; At += A;
    if (dirX && ii === 1) Qe += ax * (Tp - i.T_edge) * A; if (dirX && ii === nx - 1) Qe += ax * (Tp - i.T_edge) * A;
    if (dirY && j === 1) Qe += ay * (Tp - i.T_edge) * A; if (dirY && j === ny - 1) Qe += ay * (Tp - i.T_edge) * A;
    if (Tp > Tmax) { Tmax = Tp; jm = j; im = ii; }
  }
  let gmax = 0;
  for (let j = 1; j < ny; j++) for (let ii = 1; ii < nx; ii++) gmax = Math.max(gmax, Math.hypot((T[j][ii + 1] - T[j][ii - 1]) / (2 * dx), (T[j + 1][ii] - T[j - 1][ii]) / (2 * dy)));
  return { xs, ys, T, Tmax, Tmean: Tsum / At, Pin, Qc, Qe, it, converged: it < maxIt, res: Math.abs(Pin - Qc - Qe) / Math.max(Pin, 1e-12), qmax: i.k * gmax, xm: xs[im], ym: ys[jm], nx, ny, clampH, jm };
}

const plate2d = {
  id: 'plate2d', title: '2-D conduction in an avionics board or heat-spreader plate', fidelity: 'numerical',
  summary: 'Steady temperature map of a thin plate with component heat sources, surface convection and cooled or insulated edges, to locate hot spots and the share of heat leaving by each path.',
  equations: ['Fourier heat conduction equation', 'Newton’s law of cooling', 'Conservation of thermal energy equation', 'Fourier’s law'],
  inputs: [
    num('Lx', 'Plate length', 'm', 0.16, 0.005, 20, 'Geometry'), num('Ly', 'Plate width', 'm', 0.1, 0.005, 20, 'Geometry'),
    num('t_mm', 'Plate thickness', 'mm', 2.4, 0.05, 200, 'Geometry'),
    num('k', 'In-plane thermal conductivity', 'W/m/K', 60, 0.1, 2000, 'Material', 'FR-4 with copper planes ≈ 10–40 (effective); conduction-cooled card with thermal core ≈ 50–150; bare FR-4 ≈ 0.3'),
    sel('edges', 'Edge condition', EDGES, EDGES[1], 'Boundary conditions'),
    num('T_edge', 'Edge sink temperature', 'K', 318, 150, 1000, 'Boundary conditions', 'Card-guide or cold-wall temperature'),
    num('h', 'Surface convection coefficient (per face)', 'W/m²/K', 6, 0, 5000, 'Boundary conditions', 'Natural convection 3–8 at sea level (falls with √density at altitude); forced air 20–80'),
    num('faces', 'Cooled faces', '', 2, 0, 2, 'Boundary conditions', '', { step: 1, discrete: true }),
    num('T_air', 'Local air temperature', 'K', 318, 150, 1000, 'Boundary conditions'),
    num('P_bg', 'Uniformly distributed power', 'W', 3, 0, 1e6, 'Heat sources'),
    num('P1', 'Component 1 power', 'W', 5, 0, 1e6, 'Heat sources'), num('x1', 'Component 1 x / Lx', '-', 0.3, 0, 1, 'Heat sources'), num('y1', 'Component 1 y / Ly', '-', 0.5, 0, 1, 'Heat sources'), num('s1', 'Component 1 footprint side', 'm', 0.025, 0.001, 5, 'Heat sources'),
    num('P2', 'Component 2 power', 'W', 3, 0, 1e6, 'Heat sources'), num('x2', 'Component 2 x / Lx', '-', 0.7, 0, 1, 'Heat sources'), num('y2', 'Component 2 y / Ly', '-', 0.3, 0, 1, 'Heat sources'), num('s2', 'Component 2 footprint side', 'm', 0.02, 0.001, 5, 'Heat sources'),
    num('P3', 'Component 3 power', 'W', 2, 0, 1e6, 'Heat sources'), num('x3', 'Component 3 x / Lx', '-', 0.7, 0, 1, 'Heat sources'), num('y3', 'Component 3 y / Ly', '-', 0.75, 0, 1, 'Heat sources'), num('s3', 'Component 3 footprint side', 'm', 0.015, 0.001, 5, 'Heat sources'),
    num('T_limit', 'Allowable board temperature', 'K', 358, 200, 1500, 'Limits', '85 °C (358 K) is a common industrial component limit'),
    num('nx', 'Grid intervals along the length', '', 48, 8, 400, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => {
    // unpressurised bays: natural convection weakens with altitude (h ∝ ρ^½ for laminar free convection)
    const pressurised = (c.fuselage.cabin_dp_Pa || 0) > 0, sg = pressurised ? 1 : isa(c.mission.cruise_alt_m || c.atm.alt_m || 0).sigma, Tb = Math.max(siteT(c), 288) + 15;
    return { h: 6 * Math.sqrt(sg), T_air: Tb, T_edge: Tb };
  },
  run(i) {
    const r = plateSolve(i), C = (v) => v - 273.15, margin = i.T_limit - r.Tmax, P = i.P_bg + i.P1 + i.P2 + i.P3, warnings = [];
    if (!r.converged) warnings.push('The SOR iteration reached its limit before meeting the tolerance; the heat balance residual shows the remaining error.');
    if (r.clampH) warnings.push('With insulated edges and no surface convection the steady problem has no solution; a minimal coefficient of 0.05 W/m²/K was applied.');
    if ((i.h * i.t_mm) / 1e3 / i.k > 0.1) warnings.push('Through-thickness Biot number exceeds 0.1: the thin-plate (2-D) assumption is weak and component-side temperatures will be higher than shown.');
    return {
      kpis: [
        kp('T_max_K', 'Hot-spot temperature', r.Tmax, 'K', margin > 10 ? 'ok' : margin > 0 ? 'warn' : 'bad', `${C(r.Tmax).toFixed(1)} °C at x = ${(r.xm * 1e3).toFixed(0)} mm, y = ${(r.ym * 1e3).toFixed(0)} mm`),
        kp('T_mean_K', 'Area-mean temperature', r.Tmean, 'K'),
        kp('thermal_margin_K', 'Margin to the allowable temperature', margin, 'K', margin > 10 ? 'ok' : margin > 0 ? 'warn' : 'bad'),
        kp('P_total_W', 'Total dissipated power', P, 'W'),
        kp('Q_edge_W', 'Heat conducted to the edges', r.Qe, 'W'), kp('Q_conv_W', 'Heat convected from the faces', r.Qc, 'W'),
        kp('q_max_Wm2', 'Peak in-plane conduction flux', r.qmax, 'W/m²'),
        kp('R_hotspot_KW', 'Hot-spot rise per watt', (r.Tmax - (i.edges === EDGES[2] ? i.T_air : i.T_edge)) / Math.max(P, 1e-12), 'K/W'),
        kp('balance_residual', 'Discrete heat-balance residual', r.res, '-', r.res < 1e-5 ? 'ok' : 'warn'),
        kp('iterations', 'SOR iterations', r.it, ''),
      ],
      plots: [
        { type: 'heat', title: 'Plate temperature', xlabel: 'x [mm]', ylabel: 'y [mm]', zlabel: 'Temperature [°C]', x: r.xs.map((v) => v * 1e3), y: r.ys.map((v) => v * 1e3), z: r.T.map((row) => row.map(C)), contours: 12, equalAspect: true },
        { type: 'line', title: 'Temperature along the line through the hot spot', xlabel: 'x [mm]', ylabel: 'Temperature [°C]', series: [{ name: `y = ${(r.ym * 1e3).toFixed(0)} mm`, x: r.xs.map((v) => v * 1e3), y: r.T[r.jm].map(C) }], annotations: [{ y: C(i.T_limit), label: 'Allowable' }] },
        { type: 'bar', title: 'Heat paths', ylabel: 'Heat flow [W]', categories: ['Edges', 'Surface convection'], series: [{ name: 'Heat removed', y: [r.Qe, r.Qc] }] },
      ],
      warnings, models: ['Five-point finite-difference conduction with surface-loss term', 'Successive over-relaxation with near-optimal relaxation factor', 'Area-weighted (conservative) source deposition'],
      assumptions: ['Thin plate: temperature uniform through the thickness', 'Uniform effective in-plane conductivity', 'Uniform convection coefficient and air temperature; radiation folded into h', 'Component heat spread uniformly over its footprint'],
    };
  },
  convergence: { param: 'nx', label: 'Grid intervals along the length', levels: [12, 24, 48, 96], metric: 'T_mean_K' },
  verify() {
    const i = { Lx: 0.2, Ly: 0.2, t_mm: 2, k: 50, edges: EDGES[0], T_edge: 300, h: 0, faces: 2, T_air: 300, P_bg: 40, P1: 0, x1: 0.5, y1: 0.5, s1: 0.01, P2: 0, x2: 0.5, y2: 0.5, s2: 0.01, P3: 0, x3: 0.5, y3: 0.5, s3: 0.01, T_limit: 400, nx: 48 };
    const r = plateSolve(i); let s = 0;
    for (let m = 1; m < 200; m += 2) for (let n = 1; n < 200; n += 2) s += ((-1) ** ((m - 1) / 2 + (n - 1) / 2)) / (m * n * (m * m + n * n));
    const dT = ((16 * (40 / 0.04) * 0.04) / (50 * 0.002 * Math.PI ** 4)) * s;
    // fin limit: insulated edges in y, uniform source, convection only -> T = T_air + q''/h everywhere
    const f = plateSolve({ ...i, edges: EDGES[2], h: 10, nx: 16 });
    return [
      N.check('Centre rise of a uniformly heated square plate with cold edges', r.Tmax - 300, dT, 2e-3, 'Double Fourier series of the Poisson equation (0.0737·q″L²/kt)'),
      N.check('Uniform source balanced by convection: T − T_air = q″/(2h)', f.Tmax - 300, 40 / 0.04 / 20, 1e-6, 'Lumped surface balance'),
      N.check('Heat balance residual', r.res, 0, 1e-6, 'Source power equals edge conduction plus convection'),
    ];
  },
  calibration: { params: [{ key: 'k', min: 0.2, max: 400 }, { key: 'h', min: 0.5, max: 200 }], sweep: 'P1', target: 'T_max_K', note: 'Supply measured hot-spot temperature (thermocouple or infra-red) against component power from a bench or chamber test.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.thermal_margin_K < 0) out.push({ severity: 'critical', title: 'Board hot spot exceeds its limit', detail: `Hot spot ${(o.T_max_K - 273.15).toFixed(0)} °C against ${(i.T_limit - 273.15).toFixed(0)} °C.`, action: 'Move the hottest component towards a cooled edge, add copper planes or a metal core (higher k·t), or add forced-air or conduction cooling.', basis: 'Component temperature limit' });
    else if (o.thermal_margin_K < 10) out.push({ severity: 'warn', title: 'Thin thermal margin', detail: `Only ${o.thermal_margin_K.toFixed(1)} K below the limit.`, action: 'Component failure rates roughly double for each 10–15 K of temperature rise (rule of thumb): improving the margin extends life and reduces unscheduled removals.', basis: 'Arrhenius-type reliability trend' });
    if (o.Q_conv_W > 0.6 * o.P_total_W && i.h < 10) out.push({ severity: 'advise', title: 'Cooling relies on natural convection', detail: `${((100 * o.Q_conv_W) / o.P_total_W).toFixed(0)}% of the heat leaves by surface convection, which weakens with altitude in an unpressurised bay.`, action: 'Re-run with the coefficient reduced for the maximum bay altitude, or provide a conduction path to a cold wall.', basis: 'h ∝ ρ^½ for laminar free convection' });
    return out;
  },
};

// ---- 3. aerodynamic heating -----------------------------------------------------------------
const BL = ['turbulent', 'laminar'];
function aeroState(i, M, turbulent) {
  const a = isa(i.alt_m, i.dISA), V = M * a.a, T0 = a.T * (1 + 0.2 * M * M), Taw = recovery(a.T, M, turbulent);
  // Eckert reference-temperature method: incompressible flat-plate correlations with properties at T*
  const q = (Tw, x = i.x) => {
    const Ts = a.T + 0.5 * (Tw - a.T) + 0.22 * (Taw - a.T), rho = a.p / (R_AIR * Ts), mu = sutherland(Ts), Re = Math.max(1, (rho * V * x) / mu);
    const St = (turbulent ? 0.0296 * Re ** -0.2 : 0.332 * Re ** -0.5) * PR_AIR ** (-2 / 3), h = St * rho * V * CP_AIR;
    return { h, q: h * (Taw - Tw), Re, Ts, k: kAir(Ts), mu };
  };
  const bal = (Tw) => q(Tw).q + i.q_solar - i.eps * SIG * (Tw ** 4 - i.T_sky ** 4);
  const Teq = i.eps > 0 || V > 0 ? N.findRoot(bal, 20, Math.max(Taw, i.T_sky, 300) + 2500, 80, 1e-8) : NaN;
  return { a, V, T0, Taw, q, Teq, bal };
}
const aeroheat = {
  id: 'aeroheat', title: 'Aerodynamic heating and skin equilibrium temperature', fidelity: 'analytical',
  summary: 'Recovery temperature, boundary-layer heat flux by the Eckert reference-temperature method, leading-edge stagnation heating and the radiative-equilibrium skin temperature across the Mach range.',
  equations: ['Newton’s law of cooling', 'Stefan–Boltzmann radiation law', 'Conservation of thermal energy equation', 'Conjugate heat-transfer equations'],
  inputs: [
    num('M', 'Flight Mach number', '-', 0.78, 0.01, 8, 'Flow'), num('alt_m', 'Altitude', 'm', 10668, -500, 40000, 'Flow'), num('dISA', 'ISA deviation', 'K', 0, -60, 50, 'Flow'),
    sel('bl', 'Boundary-layer state', BL, BL[0], 'Flow', 'Turbulent is the conservative choice for heating'),
    num('x', 'Distance from the leading edge', 'm', 1, 0.005, 100, 'Geometry'),
    num('R_le', 'Leading-edge (nose) radius', 'm', 0.03, 0.0005, 5, 'Geometry'),
    num('T_wall', 'Wall temperature for heat-flux evaluation', 'K', 250, 50, 3000, 'Wall', 'Cold-wall flux is the design case for transient heating'),
    num('eps', 'Surface emissivity', '-', 0.85, 0, 1, 'Wall'), num('T_sky', 'Radiative sink temperature', 'K', 220, 3, 400, 'Wall', 'Roughly the ambient static temperature at altitude'),
    num('q_solar', 'Absorbed solar flux', 'W/m²', 0, 0, 1500, 'Wall'),
    num('T_limit', 'Allowable skin temperature', 'K', 393, 100, 3000, 'Limits', 'Aluminium ≈ 390–420 K, titanium ≈ 700 K, nickel alloys ≈ 1100 K (typical long-term figures)'),
  ],
  defaults: (c) => { const { a, M } = cruise(c); return { M: Math.max(0.02, c.mission.cruise_mach || M), alt_m: c.mission.cruise_alt_m || c.atm.alt_m, dISA: c.atm.dISA_K, x: Math.max(0.05, 0.25 * (c.fuselage.len_m || 1)), R_le: Math.max(0.002, 0.012 * (c.wing.S_m2 > 0 && c.wing.b_m > 0 ? c.wing.S_m2 / c.wing.b_m : c.rotor.chord_m || c.fuselage.dia_m || 0.5)), T_wall: a.T, T_sky: a.T }; },
  run(i) {
    const turb = i.bl === BL[0], s = aeroState(i, i.M, turb), w = s.q(i.T_wall), warnings = [];
    // stagnation line of the leading edge: cylinder in cross-flow, film properties; Sutton–Graves above Mach 3
    const Tf = 0.5 * (s.T0 + i.T_wall), rf = s.a.p / (R_AIR * Tf), ReD = (rf * s.V * 2 * i.R_le) / sutherland(Tf), hSt = (1.14 * Math.sqrt(ReD) * PR_AIR ** 0.4 * kAir(Tf)) / (2 * i.R_le);
    const hyper = i.M >= 3, qSt = hyper ? 1.7415e-4 * Math.sqrt(s.a.rho / i.R_le) * s.V ** 3 * Math.max(0, 1 - i.T_wall / s.T0) : hSt * (s.T0 - i.T_wall);
    const TeqSt = hyper ? N.findRoot((T) => 1.7415e-4 * Math.sqrt(s.a.rho / i.R_le) * s.V ** 3 * Math.max(0, 1 - T / s.T0) - i.eps * SIG * (T ** 4 - i.T_sky ** 4), 20, s.T0, 80) : N.findRoot((T) => hSt * (s.T0 - T) + i.q_solar - i.eps * SIG * (T ** 4 - i.T_sky ** 4), 20, s.T0 + 800, 80);
    const margin = i.T_limit - Math.max(s.Teq, TeqSt), ReTr = 5e5;
    if (turb && w.Re < ReTr) warnings.push(`Reynolds number ${w.Re.toExponential(2)} is below typical transition (≈5·10⁵): the turbulent estimate is conservative here.`);
    if (!turb && w.Re > 3e6) warnings.push('Reynolds number is far above typical transition: laminar flow is unlikely to persist, so the laminar heat flux is optimistic.');
    if (i.M > 1 && !hyper) warnings.push('Between Mach 1 and 3 the leading-edge estimate uses the incompressible cylinder correlation with total temperature and ignores the bow shock; treat it as indicative.');
    if (i.M > 5) warnings.push('Above about Mach 5 real-gas effects (dissociation, variable γ) are not modelled; total temperature is overestimated.');
    const Ms = N.linspace(0.05, Math.max(3, 1.3 * i.M), 60), sw = Ms.map((m) => aeroState(i, m, turb)), xs = N.logspace(Math.max(0.005, i.x / 50), i.x * 3, 50), lam = aeroState(i, i.M, false), tb = aeroState(i, i.M, true);
    return {
      kpis: [
        kp('T_recovery_K', 'Recovery (adiabatic-wall) temperature', s.Taw, 'K'),
        kp('T_total_K', 'Total (stagnation) temperature', s.T0, 'K'),
        kp('T_radeq_K', 'Radiative-equilibrium skin temperature', s.Teq, 'K', s.Teq < i.T_limit ? 'ok' : 'bad'),
        kp('T_max_K', 'Leading-edge equilibrium temperature', TeqSt, 'K', TeqSt < i.T_limit ? 'ok' : 'bad', `Limit ${i.T_limit.toFixed(0)} K`),
        kp('q_wall_Wm2', 'Convective heat flux at the stated wall temperature', w.q, 'W/m²'),
        kp('h_Wm2K', 'Heat-transfer coefficient', w.h, 'W/m²/K'),
        kp('q_stag_Wm2', 'Leading-edge stagnation heat flux', qSt, 'W/m²'),
        kp('q_max_Wm2', 'Largest heat flux (surface or leading edge)', Math.max(Math.abs(w.q), Math.abs(qSt)), 'W/m²'),
        kp('Re_x', 'Reference-temperature Reynolds number', w.Re, '-'),
        kp('T_ref_K', 'Eckert reference temperature', w.Ts, 'K'),
        kp('thermal_margin_K', 'Margin to the allowable skin temperature', margin, 'K', margin > 20 ? 'ok' : margin > 0 ? 'warn' : 'bad'),
      ],
      plots: [
        { type: 'line', title: 'Skin temperatures versus Mach number', xlabel: 'Mach number [-]', ylabel: 'Temperature [K]', series: [{ name: 'Total temperature', x: Ms, y: sw.map((v) => v.T0) }, { name: 'Recovery temperature', x: Ms, y: sw.map((v) => v.Taw) }, { name: 'Radiative equilibrium', x: Ms, y: sw.map((v) => v.Teq) }], annotations: [{ y: i.T_limit, label: 'Allowable' }, { x: i.M, label: 'Flight point' }] },
        { type: 'line', title: 'Heat flux along the surface', xlabel: 'Distance from leading edge [m]', ylabel: 'Heat flux [W/m²]', xlog: true, series: [{ name: 'Turbulent', x: xs, y: xs.map((x) => tb.q(i.T_wall, x).q) }, { name: 'Laminar', x: xs, y: xs.map((x) => lam.q(i.T_wall, x).q) }] },
      ],
      warnings, models: ['Recovery factor Pr^⅓ (turbulent) or Pr^½ (laminar)', 'Eckert reference-temperature flat-plate heat transfer', hyper ? 'Sutton–Graves stagnation-point heating (empirical, cold-wall with enthalpy-ratio correction)' : 'Cylinder stagnation-line correlation Nu = 1.14·Re^½·Pr^0.4'],
      assumptions: ['Calorically perfect air, γ = 1.4, Pr = 0.71', 'Zero pressure gradient flat plate; local edge conditions equal to free stream', 'Skin in radiative equilibrium with no conduction into the structure', 'No shock or interference heating'],
    };
  },
  verify() {
    const b = { M: 2, alt_m: 11000, dISA: 0, bl: BL[0], x: 1, R_le: 0.03, T_wall: 300, eps: 0.8, T_sky: 216.65, q_solar: 0, T_limit: 600 };
    const s = aeroState(b, 2, true), lo = aeroState({ ...b, alt_m: 0 }, 0.05, false), w = lo.q(288.15 + 0.01), Re = (1.225 * 0.05 * 340.294 * 1) / 1.7894e-5;
    return [
      N.check('Recovery temperature at Mach 2, r = Pr^⅓', s.Taw, 216.65 * (1 + 0.2 * 0.71 ** (1 / 3) * 4), 1e-9, 'T_aw = T(1 + r(γ−1)/2·M²)'),
      N.check('Low-speed laminar limit Nu_x = 0.332·Re^½·Pr^⅓', (w.h * 1) / w.k, 0.332 * Math.sqrt(Re) * 0.71 ** (1 / 3), 0.02, 'Pohlhausen flat-plate solution (property set consistent to 2%)'),
      N.check('Radiative-equilibrium balance residual', s.bal(s.Teq) / (0.8 * SIG * s.Teq ** 4), 0, 1e-6, 'q_conv = εσ(T⁴ − T_sky⁴)'),
    ];
  },
  calibration: { params: [{ key: 'eps', min: 0.02, max: 1 }, { key: 'x', min: 0.01, max: 50 }], sweep: 'M', target: 'T_radeq_K', note: 'Supply measured skin temperatures against Mach number from flight thermal survey or wind-tunnel thermography.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.thermal_margin_K < 0) out.push({ severity: 'critical', title: 'Equilibrium skin temperature exceeds the material limit', detail: `Up to ${Math.max(o.T_radeq_K, o.T_max_K).toFixed(0)} K against ${i.T_limit.toFixed(0)} K.`, action: 'Reduce Mach number or dwell time (run the transient wall analysis), blunt the leading edge, raise emissivity, or change to titanium / high-temperature composite.', basis: 'Radiative equilibrium with aerodynamic heating' });
    if (i.M < 1 && o.T_recovery_K < 273.15) out.push({ severity: 'info', title: 'Skin runs below freezing in cruise', detail: `Recovery temperature is ${(o.T_recovery_K - 273.15).toFixed(0)} °C.`, action: 'Use it as the fuel-tank and cabin-wall boundary temperature in Suite 20 and as the surface temperature for icing in Suite 13.', basis: 'Adiabatic-wall temperature' });
    return out;
  },
};

// ---- 4. heat exchanger ------------------------------------------------------------------------
const ARR = ['Counterflow', 'Parallel flow', 'Crossflow (both unmixed)'];
const HXF = {
  'Air': { cp: CP_AIR, mu: 1.9e-5, rho: null },
  'Engine oil (MIL-PRF-23699)': FLUIDS['MIL-PRF-23699 oil'], 'Jet A-1 fuel': FLUIDS['Jet A-1'], 'Hydraulic fluid (Skydrol LD-4)': FLUIDS['Skydrol LD-4'],
  'Water-glycol 50/50 (typical)': { rho: 1070, mu: 3.5e-3, cp: 3300, k: 0.4 },
};
const HXN = Object.keys(HXF);
/** Effectiveness ε(NTU, Cr) for the three arrangements. Crossflow uses the standard unmixed–unmixed approximation. */
export function hxEff(arr, NTU, Cr) {
  if (NTU <= 0) return 0;
  if (Cr < 1e-9) return 1 - Math.exp(-NTU);
  if (arr === ARR[0]) return Math.abs(1 - Cr) < 1e-9 ? NTU / (1 + NTU) : (1 - Math.exp(-NTU * (1 - Cr))) / (1 - Cr * Math.exp(-NTU * (1 - Cr)));
  if (arr === ARR[1]) return (1 - Math.exp(-NTU * (1 + Cr))) / (1 + Cr);
  return 1 - Math.exp((NTU ** 0.22 / Cr) * (Math.exp(-Cr * NTU ** 0.78) - 1));
}
function hxRate(i) {
  const fh = HXF[i.hot] || HXF[HXN[1]], fc = HXF[i.cold] || HXF.Air, Ch = i.m_hot * fh.cp, Cc = i.m_cold * fc.cp, Cmin = Math.min(Ch, Cc), Cr = Cmin / Math.max(Ch, Cc), UA = i.U * i.A, NTU = UA / Cmin;
  const eff = hxEff(i.arr, NTU, Cr), Q = eff * Cmin * (i.Th_in - i.Tc_in), Tho = i.Th_in - Q / Ch, Tco = i.Tc_in + Q / Cc;
  const d1 = i.arr === ARR[1] ? i.Th_in - i.Tc_in : i.Th_in - Tco, d2 = i.arr === ARR[1] ? Tho - Tco : Tho - i.Tc_in, lmtd = Math.abs(d1 - d2) < 1e-9 * Math.abs(d1) ? d1 : (d1 - d2) / Math.log(d1 / d2);
  return { fh, fc, Ch, Cc, Cmin, Cr, UA, NTU, eff, Q, Tho, Tco, lmtd, d1 };
}
const hx = {
  id: 'hx', title: 'Heat exchanger rating and sizing (ε–NTU)', fidelity: 'analytical',
  summary: 'Rates a fuel–oil, air–oil or air–liquid cooler: effectiveness, duty and outlet temperatures, the area needed for a required duty, and the tube-side pressure drop.',
  equations: ['Heat exchanger energy balance equations', 'Thermal resistance relations', 'Newton’s law of cooling', 'Conservation of thermal energy equation'],
  inputs: [
    sel('arr', 'Flow arrangement', ARR, ARR[0], 'Exchanger'),
    num('U', 'Overall heat-transfer coefficient', 'W/m²/K', 600, 1, 20000, 'Exchanger', 'Liquid–liquid 300–1500; air–liquid 40–150 (referred to the stated area); air–air 20–80'),
    num('A', 'Heat-transfer area', 'm²', 2, 0.001, 5000, 'Exchanger'),
    sel('hot', 'Hot fluid', HXN, HXN[1], 'Hot side'), num('m_hot', 'Hot mass flow', 'kg/s', 1, 1e-5, 1000, 'Hot side'), num('Th_in', 'Hot inlet temperature', 'K', 400, 150, 1500, 'Hot side'),
    sel('cold', 'Cold fluid', HXN, HXN[2], 'Cold side'), num('m_cold', 'Cold mass flow', 'kg/s', 0.6, 1e-5, 1000, 'Cold side'), num('Tc_in', 'Cold inlet temperature', 'K', 290, 150, 1500, 'Cold side'),
    num('Q_req', 'Required duty', 'W', 60000, 0, 1e9, 'Requirement', 'Heat that must be removed from the hot stream'),
    num('Tc_max', 'Maximum cold-side outlet temperature', 'K', 393, 200, 1500, 'Requirement', 'Jet fuel is usually kept below about 390–420 K at the engine to avoid coking (typical)'),
    num('n_tubes', 'Hot-side parallel passages', '', 60, 1, 100000, 'Hot-side pressure drop', '', { step: 1, discrete: true }),
    num('D_tube', 'Passage hydraulic diameter', 'm', 0.004, 0.0003, 0.5, 'Hot-side pressure drop'), num('L_tube', 'Passage length', 'm', 0.6, 0.01, 100, 'Hot-side pressure drop'),
  ],
  defaults: (c, up, d) => {
    const n = Math.max(1, c.prop.n_eng), turbine = ['turbofan', 'turbojet', 'turboprop', 'turboshaft'].includes(c.prop.type), elec = c.prop.type === 'electric', { a, V, M } = cruise(c);
    const Pref = d.P_total || d.T_total * Math.max(V, 60), Q = Math.max(20, (up.propulsion?.heat_rejection_W ?? (elec ? up.electrical?.elec_losses_W ?? 0.08 * Pref : (turbine ? 0.0025 : 0.08) * Pref)) / n);
    const fuel = turbine && c.mass.fuel_kg > 0, ff = (up.propulsion?.fuel_flow_kgs ?? up.performance?.fuel_flow_cruise_kgs ?? (d.T_total ? c.prop.tsfc_kg_Ns * d.T_total * 0.22 : c.prop.bsfc_kg_Ws * d.P_total * 0.7)) / n;
    const hot = elec ? HXN[4] : HXN[1], cold = fuel ? HXN[2] : HXN[0], Th = elec ? 328 : 400, Tc = fuel ? 290 : Math.min(recovery(a.T, M), 320), m_hot = Q / (HXF[hot].cp * (elec ? 8 : 30));
    const m_cold = fuel ? Math.max(ff, 1e-4) : Q / (CP_AIR * 0.45 * (Th - Tc)), Cmin = Math.min(m_hot * HXF[hot].cp, m_cold * HXF[cold].cp), U = fuel ? 600 : 90, D = fuel || !elec ? 0.004 : 0.006;
    // passages sized for about 1 m/s of liquid so that the default pressure drop is realistic
    return { hot, cold, Th_in: Th, Tc_in: Tc, m_hot, m_cold, Q_req: Q, U, A: (2 * Cmin) / U, arr: fuel ? ARR[0] : ARR[2], D_tube: D, n_tubes: Math.max(4, Math.ceil(m_hot / (HXF[hot].rho * 1.0 * Math.PI * D * D / 4))), L_tube: 0.4 };
  },
  run(i) {
    const r = hxRate(i), warnings = [], effReq = i.Q_req / Math.max(1e-9, r.Cmin * (i.Th_in - i.Tc_in));
    const effMax = i.arr === ARR[1] ? 1 / (1 + r.Cr) : 1;
    const ntuReq = effReq <= 0 ? 0 : effReq >= effMax * 0.9995 ? NaN : N.brent((n) => hxEff(i.arr, n, r.Cr) - effReq, 1e-9, 500, 1e-10), Areq = (ntuReq * r.Cmin) / i.U;
    // hot-side passage pressure drop (Darcy–Weisbach, laminar 64/Re or Blasius/Swamee–Jain smooth)
    const Tm = 0.5 * (i.Th_in + r.Tho), rho = r.fh.rho ?? 101325 / (R_AIR * Tm), v = i.m_hot / (rho * i.n_tubes * Math.PI * i.D_tube ** 2 / 4), Re = (rho * v * i.D_tube) / r.fh.mu;
    const f = Re < 2300 ? 64 / Math.max(Re, 1e-9) : 0.25 / Math.log10(5.74 / Re ** 0.9) ** 2, dp = f * (i.L_tube / i.D_tube) * 0.5 * rho * v * v, Ppump = (i.m_hot / rho) * dp;
    if (i.Th_in <= i.Tc_in) warnings.push('Hot inlet is not hotter than the cold inlet: no heat can be transferred in the intended direction.');
    if (Number.isNaN(ntuReq)) warnings.push(`The required duty needs an effectiveness of ${effReq.toFixed(2)}, beyond what this arrangement can reach with the available cold-side capacity rate: no finite area will do it.`);
    if (r.Tco > i.Tc_max) warnings.push(`Cold-side outlet ${r.Tco.toFixed(0)} K exceeds the stated limit of ${i.Tc_max.toFixed(0)} K.`);
    if (!r.fh.rho) warnings.push('Air density on the hot side is evaluated at 1 atm and the mean temperature for the pressure-drop estimate.');
    const ntus = N.linspace(0, Math.max(5, 1.5 * r.NTU), 60), xa = N.linspace(0, 1, 40), par = i.arr === ARR[1], m = i.U * i.A * (1 / r.Ch + (par ? 1 : -1) / r.Cc);
    const g = (x) => (Math.abs(m) < 1e-9 ? x : (1 - Math.exp(-m * x)) / m), Thx = xa.map((x) => i.Th_in - ((i.U * i.A * r.d1) / r.Ch) * g(x)), Tcx = xa.map((x, j) => (par ? i.Tc_in + ((i.U * i.A * r.d1) / r.Cc) * g(x) : r.Tco - ((i.U * i.A * r.d1) / r.Cc) * g(x)));
    const marginQ = i.Q_req > 0 ? r.Q / i.Q_req - 1 : Infinity, plots = [{ type: 'line', title: `Effectiveness versus NTU at Cr = ${r.Cr.toFixed(2)}`, xlabel: 'NTU [-]', ylabel: 'Effectiveness [-]', series: ARR.map((a) => ({ name: a, x: ntus, y: ntus.map((n) => hxEff(a, n, r.Cr)) })), annotations: [{ x: r.NTU, label: 'Design' }] }];
    if (i.arr !== ARR[2]) plots.push({ type: 'line', title: 'Temperature along the exchanger', xlabel: 'Fraction of area from hot inlet [-]', ylabel: 'Temperature [K]', series: [{ name: 'Hot stream', x: xa, y: Thx }, { name: 'Cold stream', x: xa, y: Tcx }] });
    return {
      kpis: [
        kp('hx_effectiveness', 'Effectiveness', r.eff, '-'),
        kp('Q_W', 'Heat duty', r.Q, 'W', marginQ >= 0 ? 'ok' : 'bad', `Required ${i.Q_req.toFixed(0)} W`),
        kp('duty_margin_pct', 'Duty margin', Number.isFinite(marginQ) ? 100 * marginQ : 0, '%', marginQ >= 0.1 ? 'ok' : marginQ >= 0 ? 'warn' : 'bad'),
        kp('T_hot_out_K', 'Hot outlet temperature', r.Tho, 'K'), kp('T_cold_out_K', 'Cold outlet temperature', r.Tco, 'K', r.Tco <= i.Tc_max ? 'ok' : 'warn'),
        kp('NTU', 'Number of transfer units', r.NTU, '-'), kp('Cr', 'Capacity-rate ratio', r.Cr, '-'),
        kp('LMTD_K', 'Log-mean temperature difference (counter/parallel basis)', r.lmtd, 'K'),
        kp('F_correction', 'LMTD correction factor', r.Q / Math.max(1e-12, r.UA * r.lmtd), '-'),
        kp('area_required_m2', 'Area for the required duty', Areq, 'm²'),
        kp('dp_hot_Pa', 'Hot-side pressure drop', dp, 'Pa'), kp('pump_power_W', 'Hot-side pumping power', Ppump, 'W'), kp('Re_tube', 'Passage Reynolds number', Re, '-'),
      ],
      plots, warnings,
      tables: [{ title: 'Stream summary', columns: ['Stream', 'Fluid', 'Capacity rate [W/K]', 'Inlet [K]', 'Outlet [K]'], rows: [['Hot', i.hot, r.Ch, i.Th_in, r.Tho], ['Cold', i.cold, r.Cc, i.Tc_in, r.Tco]] }],
      models: ['ε–NTU relations (counterflow, parallel flow; crossflow unmixed–unmixed approximation)', 'Log-mean temperature difference cross-check', 'Darcy–Weisbach passage pressure drop'],
      assumptions: ['Constant specific heats and overall coefficient', 'No heat loss to the surroundings, no phase change', 'Smooth passages; entrance, exit and header losses not included'],
    };
  },
  verify() {
    const b = { arr: ARR[0], U: 500, A: 2, hot: HXN[1], m_hot: 1, Th_in: 400, cold: HXN[1], m_cold: 1, Tc_in: 300, Q_req: 1, Tc_max: 500, n_tubes: 10, D_tube: 0.01, L_tube: 1 };
    const r1 = hxRate(b), r2 = hxRate({ ...b, m_cold: 2.5 }), r3 = hxRate({ ...b, arr: ARR[1], m_cold: 2.5 });
    return [
      N.check('Balanced counterflow ε = NTU/(1+NTU)', r1.eff, 0.5 / 1.5, 1e-9, 'Kays & London'),
      N.check('Counterflow: Q = UA·LMTD', r2.Q, r2.UA * r2.lmtd, 1e-9, 'LMTD method equivalence'),
      N.check('Parallel flow: Q = UA·LMTD', r3.Q, r3.UA * r3.lmtd, 1e-9, 'LMTD method equivalence'),
      N.check('Crossflow approximation tends to 1 − exp(−NTU) as Cr → 0', hxEff(ARR[2], 2, 1e-6), 1 - Math.exp(-2), 1e-4, 'Limit for a condensing/evaporating stream'),
    ];
  },
  calibration: { params: [{ key: 'U', min: 5, max: 10000 }], sweep: 'm_cold', target: 'Q_W', note: 'Supply measured duty against coolant flow from a heat-exchanger rig test to fit the overall coefficient.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.duty_margin_pct < 0) out.push({ severity: 'warn', title: 'Exchanger does not meet the required duty', detail: `Delivers ${(o.Q_W / 1e3).toFixed(1)} kW of ${(i.Q_req / 1e3).toFixed(1)} kW.`, action: Number.isFinite(o.area_required_m2) ? `Increase area to about ${o.area_required_m2.toFixed(2)} m², or raise the cold-side flow.` : 'Raise the cold-side flow or add a second sink (air–oil cooler): more area alone cannot deliver this duty.', basis: 'ε–NTU rating' });
    if (o.T_cold_out_K > i.Tc_max) out.push({ severity: 'warn', title: 'Cold stream overheats', detail: `Outlet ${o.T_cold_out_K.toFixed(0)} K against a limit of ${i.Tc_max.toFixed(0)} K.`, action: 'Recirculate fuel to tank, add an air–oil cooler in parallel, or reduce the heat load at low fuel flow (descent, idle).', basis: 'Fuel thermal stability limit (typical)' });
    if (i.cold === HXN[2] && o.duty_margin_pct >= 0) out.push({ severity: 'info', title: 'Waste heat is recovered into the fuel', detail: `${(o.Q_W / 1e3).toFixed(1)} kW pre-heats the fuel instead of being rejected to ram air.`, action: 'Using fuel as the heat sink avoids ram-air cooler drag and returns the heat to the cycle; keep it as the primary sink where fuel flow allows.', basis: 'Energy recovery' });
    if (o.NTU > 4) out.push({ severity: 'advise', title: 'Exchanger is oversized for its flows', detail: `NTU = ${o.NTU.toFixed(1)}: effectiveness gains are flat beyond NTU ≈ 3–4.`, action: 'Area (mass and volume) can be reduced with little loss of duty.', basis: 'ε–NTU curve saturation' });
    return out;
  },
};

// ---- 5. lumped thermal network -----------------------------------------------------------------
const MODES = ['Natural convection + radiation', 'Forced air + radiation', 'Liquid cold plate'];
/** Sink-to-ambient conductance [W/K] with density-dependent air properties at the bay pressure. */
function sinkG(i, Ts, p) {
  if (i.mode === MODES[2]) return { G: 1 / i.R_cp, h: 1 / (i.R_cp * i.A_sink) };
  const Tf = 0.5 * (Ts + i.T_amb), rho = p / (R_AIR * Tf), mu = sutherland(Tf), k = kAir(Tf); let Nu;
  if (i.mode === MODES[0]) { const Ra = (G0 * (Math.abs(Ts - i.T_amb) / Tf) * i.L_sink ** 3 * rho * rho * CP_AIR) / (mu * k); Nu = (0.825 + (0.387 * Ra ** (1 / 6)) / (1 + (0.492 / PR_AIR) ** (9 / 16)) ** (8 / 27)) ** 2; }
  else { const Re = (rho * i.v_air * i.L_sink) / mu; Nu = Re < 5e5 ? 0.664 * Math.sqrt(Re) * PR_AIR ** (1 / 3) : (0.037 * Re ** 0.8 - 871) * PR_AIR ** (1 / 3); }
  const h = (Nu * k) / i.L_sink, hr = i.eps * SIG * (Ts * Ts + i.T_amb * i.T_amb) * (Ts + i.T_amb);
  return { G: (h + hr) * i.A_sink, h, hr };
}
function netSteady(i, P, p) {
  const T3 = P <= 0 ? i.T_amb : N.brent((T) => sinkG(i, T, p).G * (T - i.T_amb) - P, i.T_amb, i.T_amb + 5000, 1e-9), T2 = T3 + P * i.R_cs;
  return [T2 + P * i.R_jc, T2, T3];
}
function netSolve(i) {
  const p = isa(i.alt_bay).p, nt = Math.max(2, Math.round(i.nSteps)), dt = i.t_end / nt, C = [i.C1, i.C2, i.C3], g12 = 1 / i.R_jc, g23 = 1 / i.R_cs;
  const Pt = (t) => i.P_W * (t < i.t_peak ? i.peak : 1);
  let T = [i.T_init, i.T_init, i.T_init], Tp = T; const time = [0], H = [T.slice()];
  for (let s = 1; s <= nt; s++) {
    const t = s * dt, bdf = s > 1, a = bdf ? 1.5 / dt : 1 / dt; let Tn = T;
    for (let it = 0; it < 3; it++) { // BDF2 (first step backward Euler), sink conductance lagged by Picard iteration
      const gs = sinkG(i, Tn[2], p).G;
      const A = [[a * C[0] + g12, -g12, 0], [-g12, a * C[1] + g12 + g23, -g23], [0, -g23, a * C[2] + g23 + gs]];
      const b = C.map((c, j) => (bdf ? (c * (4 * T[j] - Tp[j])) / (2 * dt) : (c * T[j]) / dt)); b[0] += Pt(t); b[2] += gs * i.T_amb;
      Tn = N.solve(A, b);
    }
    Tp = T; T = Tn; time.push(t); H.push(T.slice());
  }
  return { time, H, p };
}
const network = {
  id: 'network', title: 'Electronics and battery cooling: lumped thermal network', fidelity: 'reduced-order',
  summary: 'Junction (or cell), case and heat-sink temperatures of a powered unit through a resistance–capacitance network, with natural, forced-air or liquid cooling and the loss of air cooling at altitude.',
  equations: ['Thermal resistance relations', 'Conservation of thermal energy equation', 'Newton’s law of cooling', 'Stefan–Boltzmann radiation law', 'Transient heat diffusion equation'],
  inputs: [
    num('P_W', 'Steady heat dissipation', 'W', 150, 0, 1e7, 'Heat load', 'Filled from the electrical suite losses when available'),
    num('peak', 'Peak power factor', '-', 1.5, 0.1, 20, 'Heat load', 'Multiplier applied during the initial peak (e.g. take-off or hover)'), num('t_peak', 'Peak duration', 's', 120, 0, 1e5, 'Heat load'),
    num('R_jc', 'Junction/cell-to-case resistance', 'K/W', 0.05, 1e-7, 1000, 'Network'), num('R_cs', 'Case-to-sink (interface) resistance', 'K/W', 0.03, 1e-7, 1000, 'Network', 'Includes thermal interface material and contact resistance'),
    num('C1', 'Junction/cell heat capacity', 'J/K', 3000, 0, 1e9, 'Network'), num('C2', 'Case heat capacity', 'J/K', 6000, 0, 1e9, 'Network'), num('C3', 'Sink heat capacity', 'J/K', 9000, 0, 1e9, 'Network'),
    sel('mode', 'Cooling mode', MODES, MODES[0], 'Cooling'),
    num('A_sink', 'Sink wetted area', 'm²', 0.5, 1e-4, 1000, 'Cooling', 'Total finned surface exchanging with the air'), num('L_sink', 'Sink characteristic length (height)', 'm', 0.15, 0.005, 10, 'Cooling'),
    num('eps', 'Sink emissivity', '-', 0.85, 0, 1, 'Cooling', 'Black anodised ≈ 0.85'), num('v_air', 'Forced-air velocity', 'm/s', 4, 0.05, 100, 'Cooling'),
    num('R_cp', 'Cold-plate resistance to coolant', 'K/W', 0.05, 1e-6, 100, 'Cooling', 'Liquid mode only'),
    num('T_amb', 'Ambient air or coolant temperature', 'K', 313, 150, 600, 'Environment'), num('alt_bay', 'Pressure altitude of the equipment bay', 'm', 0, -500, 25000, 'Environment', 'Cabin altitude for pressurised bays'),
    num('T_init', 'Initial temperature', 'K', 303, 150, 600, 'Initial state'), num('T_limit', 'Allowable junction/cell temperature', 'K', 358, 200, 1000, 'Limits', 'Li-ion cells ≈ 318–333 K in operation; silicon junctions 398–423 K; equipment cases 343–358 K (typical)'),
    num('dT_cool', 'Allowed coolant temperature rise', 'K', 10, 0.5, 100, 'Cooling'),
    num('t_end', 'Simulated time', 's', 3600, 1, 1e7, 'Numerics'), num('nSteps', 'Time steps', '', 300, 10, 100000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c, up, d) => {
    const elec = c.prop.type === 'electric', P = Math.max(5, up.electrical?.elec_losses_W ?? (elec ? 0.06 * d.P_total : 0.04 * (c.systems.gen_kVA || 0.1) * 1000 * Math.max(1, c.prop.n_eng) + 20));
    const pressurised = (c.fuselage.cabin_dp_Pa || 0) > 0, big = P > 3000, Ta = Math.max(siteT(c), 288) + (big ? 5 : elec ? 10 : 15), lim = elec ? 333 : 358;
    return { P_W: P, R_jc: 4 / P, R_cs: 3 / P, R_cp: 8 / P, C1: 25 * P, C2: 30 * P, C3: 40 * P, mode: big ? MODES[2] : elec ? MODES[1] : MODES[0], A_sink: Math.max(0.02, P / (elec ? 400 : 250)), L_sink: N.clamp(0.02 * Math.sqrt(P), 0.03, 0.6), v_air: elec ? Math.max(4, 0.4 * (c.flight.V_ms || 10)) : 4, T_amb: Ta, T_init: siteT(c), alt_bay: pressurised ? c.systems.cabin_alt_m || 2400 : c.mission.cruise_alt_m || c.atm.alt_m || 0, T_limit: lim, t_peak: (c.mission.hover_min || 2) * 60, peak: elec ? 1.6 : 1.2 };
  },
  run(i) {
    const r = netSolve(i), ss = netSteady(i, i.P_W, r.p), pk = netSteady(i, i.P_W * i.peak, r.p), Tj = r.H.map((h) => h[0]), TjMax = N.amax(Tj), g = sinkG(i, ss[2], r.p), g0 = sinkG(i, ss[2], 101325);
    const margin = i.T_limit - Math.max(TjMax, ss[0]), cpCool = i.mode === MODES[2] ? 3300 : CP_AIR, mdot = i.P_W / (cpCool * i.dT_cool), warnings = [];
    const Rtot = i.R_jc + i.R_cs + 1 / g.G, tau = (i.C1 + i.C2 + i.C3) * Math.max(1 / g.G, 1e-12) + i.C1 * i.R_jc + (i.C1 + i.C2) * i.R_cs;
    if (i.mode !== MODES[2] && r.p < 30000) warnings.push('Bay pressure is below 30 kPa: air cooling is severely weakened and the continuum correlations lose accuracy; consider conduction or liquid cooling.');
    if (Math.abs(Tj[Tj.length - 1] - ss[0]) > 0.05 * Math.abs(ss[0] - i.T_init) + 0.5) warnings.push('The unit has not reached thermal steady state within the simulated time; the steady KPI shows where it is heading.');
    const alts = N.linspace(0, 14000, 29), plots = [{ type: 'line', title: 'Temperature histories', xlabel: 'Time [s]', ylabel: 'Temperature [K]', series: ['Junction / cell', 'Case', 'Sink'].map((name, j) => ({ name, x: thin(r.time), y: thin(r.H.map((h) => h[j])) })), annotations: [{ y: i.T_limit, label: 'Allowable' }] }];
    if (i.mode !== MODES[2]) plots.push({ type: 'line', title: 'Steady junction temperature versus bay altitude', xlabel: 'Bay pressure altitude [m]', ylabel: 'Temperature [K]', series: [{ name: 'Junction / cell (steady)', x: alts, y: alts.map((h) => netSteady(i, i.P_W, isa(h).p)[0]) }], annotations: [{ y: i.T_limit, label: 'Allowable' }, { x: i.alt_bay, label: 'Bay altitude' }] });
    return {
      kpis: [
        kp('T_max_K', 'Peak junction/cell temperature in the run', TjMax, 'K', TjMax < i.T_limit ? 'ok' : 'bad', `Limit ${i.T_limit.toFixed(0)} K`),
        kp('T_junction_ss_K', 'Steady junction/cell temperature', ss[0], 'K', ss[0] < i.T_limit ? 'ok' : 'bad'),
        kp('T_junction_peak_ss_K', 'Steady temperature if the peak power were sustained', pk[0], 'K', pk[0] < i.T_limit ? 'ok' : 'warn'),
        kp('T_sink_ss_K', 'Steady sink temperature', ss[2], 'K'),
        kp('thermal_margin_K', 'Margin to the allowable temperature', margin, 'K', margin > 10 ? 'ok' : margin > 0 ? 'warn' : 'bad'),
        kp('R_total_KW', 'Total junction-to-ambient resistance', Rtot, 'K/W'), kp('h_sink_Wm2K', 'Sink convection coefficient', g.h, 'W/m²/K'),
        kp('altitude_derate', 'Sink conductance relative to sea level', g.G / g0.G, '-', undefined, 'Loss of air cooling with density'),
        kp('time_constant_s', 'Dominant thermal time constant (Elmore estimate)', tau, 's'),
        kp('coolant_flow_kgs', 'Coolant mass flow for the allowed temperature rise', mdot, 'kg/s', undefined, i.mode === MODES[2] ? 'Water-glycol, cp ≈ 3300 J/kg/K' : 'Air'),
        kp('T_end_K', 'Junction/cell temperature at end of run', Tj[Tj.length - 1], 'K'),
      ],
      plots, warnings,
      models: ['Three-node resistance–capacitance network, BDF2 implicit integration', i.mode === MODES[0] ? 'Churchill–Chu vertical-plate natural convection (empirical)' : i.mode === MODES[1] ? 'Flat-plate forced convection, laminar/mixed (empirical)' : 'Prescribed cold-plate resistance', 'Linearised grey-body radiation from the sink'],
      assumptions: ['Each node is isothermal (lumped)', 'Air properties at film temperature and bay pressure', 'Constant ambient or coolant temperature', 'Single dissipating node; no temperature dependence of the heat load'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [25, 50, 100, 200, 400], metric: 'T_end_K' },
  verify() {
    const b = { P_W: 100, peak: 1, t_peak: 0, R_jc: 0.2, R_cs: 0.1, C1: 500, C2: 0, C3: 0, mode: MODES[2], A_sink: 1, L_sink: 0.1, eps: 0, v_air: 1, R_cp: 0.2, T_amb: 300, alt_bay: 0, T_init: 300, T_limit: 400, dT_cool: 10, t_end: 300, nSteps: 600 };
    const r = netSolve(b), ss = netSteady({ ...b, mode: MODES[1], eps: 0, v_air: 5 }, 100, 101325), g = sinkG({ ...b, mode: MODES[1], v_air: 5 }, ss[2], 101325);
    return [
      N.check('Single-capacity step response T = T∞ + P·R(1 − e^(−t/RC))', r.H[600][0] - 300, 100 * 0.5 * (1 - Math.exp(-300 / (0.5 * 500))), 1e-4, 'Lumped-capacitance solution'),
      N.check('Steady series chain ΔT = P·ΣR', ss[0] - 300, 100 * (0.3 + 1 / g.G), 1e-8, 'Thermal resistance network'),
    ];
  },
  calibration: { params: [{ key: 'R_cs', min: 1e-5, max: 10 }, { key: 'A_sink', min: 1e-3, max: 100 }, { key: 'C2', min: 1, max: 1e7 }], sweep: 't_end', target: 'T_end_K', note: 'Supply measured junction or cell temperature against time after a power step from a thermal-chamber test.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.thermal_margin_K < 0) out.push({ severity: 'critical', title: 'Unit overheats', detail: `Reaches ${Math.max(o.T_max_K, o.T_junction_ss_K).toFixed(0)} K against a limit of ${i.T_limit.toFixed(0)} K.`, action: i.mode === MODES[0] ? 'Add forced air or a cold plate, enlarge the sink, or derate the power at altitude.' : 'Lower the interface resistance, increase coolant flow or sink area, or reduce the duty cycle.', basis: 'Junction/cell temperature limit' });
    if (o.altitude_derate < 0.75 && i.mode !== MODES[2]) out.push({ severity: 'advise', title: 'Air cooling is degraded at bay altitude', detail: `Sink conductance is ${(100 * o.altitude_derate).toFixed(0)}% of its sea-level value.`, action: 'Qualify the unit at maximum bay altitude; a conduction path to structure or a liquid loop is insensitive to altitude.', basis: 'Density dependence of convection' });
    if (o.thermal_margin_K > 40) out.push({ severity: 'info', title: 'Generous thermal margin', detail: `${o.thermal_margin_K.toFixed(0)} K below the limit.`, action: 'Sink mass or fan power can be reduced; cooler electronics and cells also age more slowly, so trade mass against life.', basis: 'Mass and energy saving' });
    return out;
  },
};

// ---- 6. thermal stress ----------------------------------------------------------------------
const CONSTR = ['Uniaxial, fully restrained', 'Biaxial plate, fully restrained', 'Partially restrained (uniaxial)'];
/** Cycles to failure from the strain-life curve εa = (σf'/E)(2N)^b + εf'(2N)^c. */
function strainLife(m, ea) {
  const f = (lg) => (m.sf / m.E) * (2 * 10 ** lg) ** m.b + m.ef * (2 * 10 ** lg) ** m.c - ea;
  if (f(0) <= 0) return 1; if (f(12) >= 0) return 1e12; // run-out cap
  return 10 ** N.brent(f, 0, 12, 1e-10);
}
function stressCalc(i) {
  const m = metal(i.mat), m2 = metal(i.mat2), fr = i.constraint === CONSTR[2] ? i.restraint : 1, bi = i.constraint === CONSTR[1] ? 1 / (1 - m.nu) : 1;
  const sig = m.E * m.alpha * i.dT * fr * bi, sigGrad = (m.E * m.alpha * i.dT_grad) / (2 * (1 - m.nu));
  // bimaterial strip (Timoshenko 1925): curvature when free, interface force when held flat
  const t1 = i.t1_mm / 1e3, t2 = i.t2_mm / 1e3, hh = t1 + t2, mm = t1 / t2, nn = m.E / m2.E, da = m2.alpha - m.alpha;
  const kappa = (6 * da * i.dT * (1 + mm) ** 2) / (hh * (3 * (1 + mm) ** 2 + (1 + mm * nn) * (mm * mm + 1 / (mm * nn))));
  const F = (da * i.dT) / (1 / (m.E * t1) + 1 / (m2.E * t2)), s1 = F / t1, s2 = -F / t2;
  const sTot = Math.abs(sig) + Math.abs(sigGrad), ea = (i.Kt * sTot) / (2 * m.E), life = strainLife(m, ea);
  return { m, m2, sig, sigGrad, sTot, kappa, s1, s2, ea, life, free: m.alpha * i.dT * i.L, tip: 0.5 * kappa * i.L ** 2 };
}
const stress = {
  id: 'stress', title: 'Thermal stress, expansion mismatch and thermal fatigue', fidelity: 'analytical',
  summary: 'Stress from restrained thermal expansion, through-thickness temperature gradients and dissimilar-material joints, with a strain-life estimate of thermal cycles to crack initiation.',
  equations: ['Thermal expansion equations', 'Thermomechanical stress equations'],
  inputs: [
    sel('mat', 'Material', METAL_NAMES, 'Al 2024-T3', 'Material'),
    num('dT', 'Temperature change of the part', 'K', -70, -1500, 1500, 'Loading', 'Relative to the stress-free (assembly) temperature; ground-to-cruise is filled from the case'),
    sel('constraint', 'Restraint', CONSTR, CONSTR[2], 'Loading'),
    num('restraint', 'Restraint fraction', '-', 0.3, 0, 1, 'Loading', '0 = free to expand, 1 = rigidly held. Built-up airframe joints are typically 0.1–0.5'),
    num('dT_grad', 'Through-thickness temperature difference', 'K', 10, -1000, 1000, 'Loading', 'Hot face minus cold face, bending fully restrained'),
    num('L', 'Part length', 'm', 0.5, 0.001, 100, 'Geometry'),
    sel('mat2', 'Joined material (bimaterial)', METAL_NAMES, 'Ti-6Al-4V', 'Bimaterial joint'), num('t1_mm', 'Thickness of material 1', 'mm', 2, 0.01, 500, 'Bimaterial joint'), num('t2_mm', 'Thickness of material 2', 'mm', 2, 0.01, 500, 'Bimaterial joint'),
    num('Kt', 'Stress concentration factor', '-', 2.5, 1, 10, 'Fatigue', 'Fastener hole ≈ 2.5–3'),
    num('cycles_per_flight', 'Thermal cycles per flight', '', 1, 0.01, 1000, 'Fatigue'),
  ],
  defaults: (c) => ({ mat: METALS[c.struct.material] ? c.struct.material : undefined, mat2: c.struct.material === 'Ti-6Al-4V' ? 'Al 2024-T3' : 'Ti-6Al-4V', dT: Math.min(-5, recovery(cruise(c).a.T, cruise(c).M) - Math.max(siteT(c), 288.15)), t1_mm: c.struct.t_skin_mm, t2_mm: c.struct.t_skin_mm, L: N.clamp(0.05 * (c.wing.b_m || c.fuselage.len_m || 2), 0.05, 0.5) }),
  run(i) {
    const r = stressCalc(i), m = r.m, ms = m.Sy / Math.max(r.sTot, 1e-9) - 1, msB = Math.min(m.Sy / Math.max(Math.abs(r.s1), 1e-9), r.m2.Sy / Math.max(Math.abs(r.s2), 1e-9)) - 1, warnings = [];
    if (r.sTot > m.Sy) warnings.push('Elastic thermal stress exceeds yield: the part will yield or buckle and the linear result overstates the stress; the fatigue estimate (elastic strain) understates plastic strain.');
    if (r.sig < 0 && i.constraint !== CONSTR[2]) warnings.push('The restrained part is in compression: check thin skins for thermal buckling in Suite 2.');
    if (Math.abs(i.dT) > 150) warnings.push('Large temperature change: room-temperature modulus, expansion coefficient and strength are used; high-temperature knock-downs are not applied.');
    const dTs = N.linspace(0, Math.max(50, 1.5 * Math.abs(i.dT)), 30), Nf = N.logspace(10, 1e8, 50);
    return {
      kpis: [
        kp('sigma_thermal_Pa', 'Stress from restrained expansion', r.sig, 'Pa', undefined, r.sig < 0 ? 'Compressive' : 'Tensile'),
        kp('sigma_gradient_Pa', 'Surface stress from the temperature gradient', r.sigGrad, 'Pa'),
        kp('sigma_total_Pa', 'Combined thermal stress magnitude', r.sTot, 'Pa', ms > 0.5 ? 'ok' : ms > 0 ? 'warn' : 'bad'),
        kp('ms_yield', 'Margin of safety on yield', ms, '-', ms > 0.5 ? 'ok' : ms > 0 ? 'warn' : 'bad'),
        kp('free_expansion_m', 'Free thermal growth over the part length', r.free, 'm'),
        kp('bimaterial_stress1_Pa', 'Bimaterial joint stress, material 1 (held flat)', r.s1, 'Pa', msB > 0 ? 'ok' : 'bad'),
        kp('bimaterial_stress2_Pa', 'Bimaterial joint stress, material 2 (held flat)', r.s2, 'Pa'),
        kp('bimaterial_curvature', 'Free bimaterial curvature', r.kappa, '1/m'), kp('bimaterial_tip_m', 'Free bimaterial bow over the part length', r.tip, 'm'),
        kp('strain_amplitude', 'Local strain amplitude per thermal cycle', r.ea, '-'),
        kp('thermal_life_cycles', 'Thermal cycles to crack initiation', r.life, 'cycles', r.life > 1e5 ? 'ok' : r.life > 1e4 ? 'warn' : 'bad', r.life >= 1e12 ? 'Run-out: beyond 10¹² cycles on the strain-life curve' : ''),
        kp('thermal_life_flights', 'Equivalent flights', r.life / i.cycles_per_flight, 'flights'),
      ],
      plots: [
        { type: 'line', title: 'Thermal stress versus temperature change', xlabel: '|ΔT| [K]', ylabel: 'Stress [MPa]', series: [{ name: 'Uniaxial restrained', x: dTs, y: dTs.map((d) => (m.E * m.alpha * d) / 1e6) }, { name: 'Biaxial restrained', x: dTs, y: dTs.map((d) => (m.E * m.alpha * d) / (1 - m.nu) / 1e6) }, { name: 'As specified', x: dTs, y: dTs.map((d) => (Math.abs(r.sig) / Math.max(Math.abs(i.dT), 1e-9)) * d / 1e6) }], annotations: [{ y: m.Sy / 1e6, label: 'Yield' }, { x: Math.abs(i.dT), label: 'Design ΔT' }] },
        { type: 'line', title: 'Strain-life curve and operating point', xlabel: 'Cycles to initiation [-]', ylabel: 'Strain amplitude [-]', xlog: true, ylog: true, series: [{ name: i.mat, x: Nf, y: Nf.map((n) => (m.sf / m.E) * (2 * n) ** m.b + m.ef * (2 * n) ** m.c) }, { name: 'Operating point', x: [Math.min(r.life, 1e8)], y: [Math.max(r.ea, 1e-6)], style: 'points' }] },
      ],
      warnings, models: ['Restrained thermal expansion σ = EαΔT (÷(1−ν) biaxial)', 'Linear-gradient plate bending stress EαΔT/(2(1−ν))', 'Timoshenko bimaterial strip', 'Coffin–Manson–Basquin strain-life with elastic local strain Kt·σ/E'],
      assumptions: ['Linear elastic, temperature-independent properties (handbook room-temperature values, not design allowables)', 'Uniform temperature change plus a linear through-thickness gradient', 'No mean-stress or creep–fatigue interaction', 'Perfect bond in the bimaterial joint; edge peel stresses not evaluated'],
    };
  },
  verify() {
    const b = { mat: 'Steel 4340 (QT)', dT: 100, constraint: CONSTR[0], restraint: 1, dT_grad: 0, L: 1, mat2: 'Steel 4340 (QT)', t1_mm: 1, t2_mm: 1, Kt: 1, cycles_per_flight: 1 }, r = stressCalc(b), m = METALS['Steel 4340 (QT)'];
    // equal-thickness, equal-modulus strip with a fictitious expansion difference: κ = 3Δα·ΔT/(2h)
    const a2 = { ...METALS['Al 2024-T3'] }, tm = (da) => { const hh = 0.002; return (6 * da * 100 * 4) / (hh * (3 * 4 + 2 * 2)); };
    const bi = stressCalc({ ...b, mat: 'Al 2024-T3', mat2: 'Al 7075-T6', t1_mm: 2, t2_mm: 2 }), F = ((METALS['Al 7075-T6'].alpha - a2.alpha) * 100) / (1 / (a2.E * 0.002) + 1 / (METALS['Al 7075-T6'].E * 0.002));
    const Nl = strainLife(m, 0.004);
    return [
      N.check('Fully restrained bar σ = EαΔT', r.sig, 205e9 * 12.3e-6 * 100, 1e-12, 'Timoshenko & Goodier'),
      N.check('Timoshenko strip, m = n = 1: κ = 1.5·Δα·ΔT/h', tm(1e-6), (1.5 * 1e-6 * 100) / 0.002, 1e-12, 'Timoshenko (1925)'),
      N.check('Bimaterial joint force balance', bi.s1 * 0.002, F, 1e-12, 'Compatibility of two bonded bars'),
      N.check('Strain-life inversion round trip', (m.sf / m.E) * (2 * Nl) ** m.b + m.ef * (2 * Nl) ** m.c, 0.004, 1e-8, 'Coffin–Manson–Basquin'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.ms_yield < 0) out.push({ severity: 'critical', title: 'Thermal stress exceeds yield', detail: `${(o.sigma_total_Pa / 1e6).toFixed(0)} MPa combined thermal stress.`, action: 'Introduce expansion joints, slotted fastener holes or flexible mounts to lower the restraint, or match expansion coefficients.', basis: 'Yield criterion' });
    if (o.thermal_life_cycles < 1e5) out.push({ severity: 'warn', title: 'Thermal fatigue is a life driver', detail: `About ${o.thermal_life_cycles.toExponential(1)} cycles to initiation at this strain amplitude.`, action: 'Pass the thermal stress range to Suite 9 for combination with the mechanical spectrum; reducing restraint or Kt extends life and avoids early part replacement.', basis: 'Strain-life (Coffin–Manson–Basquin)' });
    if (Math.abs(o.bimaterial_stress1_Pa) > 0.3 * metal(i.mat).Sy) out.push({ severity: 'advise', title: 'Significant mismatch stress in the dissimilar joint', detail: `${(Math.abs(o.bimaterial_stress1_Pa) / 1e6).toFixed(0)} MPa in ${i.mat} joined to ${i.mat2}.`, action: 'Use a material pair with closer expansion coefficients (e.g. titanium against carbon composite), or a compliant interlayer.', basis: 'Expansion mismatch' });
    return out;
  },
};

// ---- 7. radiation exchange and ramp hot soak --------------------------------------------------
const SURR = ['Black surroundings at stated temperature', 'Re-radiating (adiabatic) shroud'];
/** View factor between coaxial parallel disks of radii r1, r2 at spacing L. */
export function diskF(r1, r2, L) { const R1 = r1 / L, R2 = r2 / L, S = 1 + (1 + R2 * R2) / (R1 * R1); return 0.5 * (S - Math.sqrt(S * S - 4 * (r2 / r1) ** 2)); }
function radSolve(i) {
  const A1 = Math.PI * i.r1 ** 2, A2 = Math.PI * i.r2 ** 2, F12 = diskF(i.r1, i.r2, i.gap), F21 = (F12 * A1) / A2, F13 = 1 - F12, F23 = 1 - F21, E = (T) => SIG * T ** 4;
  // radiosity network: unknown J1, J2 (and J3 for a re-radiating shroud); surface resistance (1−ε)/(εA)
  const g1 = (i.eps1 * A1) / Math.max(1e-12, 1 - i.eps1), g2 = (i.eps2 * A2) / Math.max(1e-12, 1 - i.eps2), rerad = i.surr === SURR[1];
  const M = [[g1 + A1 * F12 + A1 * F13, -A1 * F12, -A1 * F13], [-A2 * F21, g2 + A2 * F21 + A2 * F23, -A2 * F23], rerad ? [-A1 * F13, -A2 * F23, A1 * F13 + A2 * F23] : [0, 0, 1]];
  const J = N.solve(M, [g1 * E(i.T1), g2 * E(i.T2), rerad ? 0 : E(i.T3)]), q1 = g1 * (E(i.T1) - J[0]), q2 = g2 * (E(i.T2) - J[1]);
  return { A1, A2, F12, F21, J, q1, q2, q3: -(q1 + q2), q12: A1 * F12 * (J[0] - J[1]), T3: (J[2] / SIG) ** 0.25 };
}
/** Steady sunlit skin temperature: absorbed solar + convection + sky/ground radiation + conduction to the interior. */
function soak(i, alpha = i.alpha_s, wind = i.wind) {
  const h = 5.7 + 3.8 * wind; // Jürges/McAdams wind correlation (empirical)
  const f = (T) => alpha * i.G_solar + h * (i.T_air - T) + i.eps_ir * SIG * (i.T_skyr ** 4 - T ** 4) + i.U_in * (i.T_in - T);
  return { T: N.brent(f, 100, 800, 1e-9), h };
}
const radiation = {
  id: 'radiation', title: 'Radiation exchange and ramp hot-soak temperature', fidelity: 'analytical',
  summary: 'Net radiative heat exchange between two facing surfaces inside an enclosure by the radiosity method, and the sunlit skin temperature of a parked aircraft from the site weather.',
  equations: ['Stefan–Boltzmann radiation law', 'Radiative transfer equation', 'Conservation of thermal energy equation', 'Newton’s law of cooling'],
  inputs: [
    num('T1', 'Surface 1 temperature (hot)', 'K', 600, 3, 3000, 'Enclosure', 'e.g. engine casing, brake, exhaust duct'), num('eps1', 'Surface 1 emissivity', '-', 0.6, 0.01, 1, 'Enclosure'), num('r1', 'Surface 1 equivalent radius', 'm', 0.3, 0.001, 50, 'Enclosure'),
    num('T2', 'Surface 2 temperature (receiver)', 'K', 330, 3, 3000, 'Enclosure', 'e.g. nacelle skin, wheel-well structure'), num('eps2', 'Surface 2 emissivity', '-', 0.3, 0.01, 1, 'Enclosure', 'Polished heat shield ≈ 0.05–0.1; oxidised or painted metal 0.6–0.9'), num('r2', 'Surface 2 equivalent radius', 'm', 0.3, 0.001, 50, 'Enclosure'),
    num('gap', 'Separation', 'm', 0.1, 0.0005, 50, 'Enclosure'),
    sel('surr', 'Third surface', SURR, SURR[0], 'Enclosure'), num('T3', 'Surroundings temperature', 'K', 320, 3, 3000, 'Enclosure'),
    num('G_solar', 'Solar irradiance on the skin', 'W/m²', 1000, 0, 1400, 'Ramp hot soak', 'Clear-sky noon ≈ 900–1100 on a horizontal surface'),
    num('alpha_s', 'Solar absorptivity of the finish', '-', 0.3, 0.02, 1, 'Ramp hot soak', 'White paint 0.2–0.3, grey 0.5–0.7, dark or bare carbon 0.9'),
    num('eps_ir', 'Infra-red emissivity of the finish', '-', 0.9, 0.02, 1, 'Ramp hot soak'),
    num('T_air', 'Ambient air temperature', 'K', 303, 200, 340, 'Ramp hot soak', 'Filled from the site weather'), num('wind', 'Wind speed', 'm/s', 1, 0, 40, 'Ramp hot soak', 'Filled from the site weather'),
    num('T_skyr', 'Effective sky temperature', 'K', 283, 150, 340, 'Ramp hot soak', 'Clear sky is roughly 15–25 K below air temperature; overcast is close to air temperature'),
    num('U_in', 'Conductance from skin to interior', 'W/m²/K', 1, 0, 500, 'Ramp hot soak'), num('T_in', 'Interior air temperature', 'K', 303, 200, 400, 'Ramp hot soak'),
    num('T_limit', 'Allowable skin temperature', 'K', 353, 250, 1000, 'Limits', 'Wet-conditioned epoxy composites are often limited to about 345–365 K (typical)'),
  ],
  defaults: (c) => { const Ta = siteT(c), comp = !METALS[c.struct.material]; return { T_air: Ta, wind: c.site?.wind_ms ?? 1, T_skyr: Ta - 20, T_in: Ta + 5, T3: Ta + 15, T2: Ta + 30, alpha_s: comp ? 0.5 : 0.3 }; },
  run(i) {
    const r = radSolve(i), s = soak(i), hr = r.q12 / Math.max(1e-9, r.A1 * (i.T1 - i.T2)), margin = i.T_limit - s.T, warnings = [];
    if (i.gap < 0.02 * Math.min(i.r1, i.r2)) warnings.push('Very small gap: conduction and convection across the gap (not modelled here) may exceed radiation.');
    const al = N.linspace(0.1, 1, 30), ws = N.linspace(0, 15, 30), e2 = N.linspace(0.03, 1, 30);
    return {
      kpis: [
        kp('F12', 'View factor surface 1 → 2', r.F12, '-'),
        kp('q12_W', 'Net radiative exchange 1 → 2', r.q12, 'W'), kp('q1_W', 'Net heat leaving surface 1', r.q1, 'W'), kp('q2_W', 'Net heat leaving surface 2', r.q2, 'W', undefined, 'Negative = surface 2 absorbs heat'),
        kp('q_max_Wm2', 'Radiative flux absorbed by surface 2', Math.abs(r.q2) / r.A2, 'W/m²'),
        kp('h_rad_Wm2K', 'Equivalent radiation coefficient', hr, 'W/m²/K'),
        kp('T_shroud_K', i.surr === SURR[1] ? 'Re-radiating shroud temperature' : 'Surroundings temperature', r.T3, 'K'),
        kp('T_soak_K', 'Sunlit skin temperature on the ramp', s.T, 'K', margin > 10 ? 'ok' : margin > 0 ? 'warn' : 'bad', `${(s.T - 273.15).toFixed(0)} °C`),
        kp('T_max_K', 'Hot-soak skin temperature', s.T, 'K'),
        kp('soak_rise_K', 'Skin temperature above ambient', s.T - i.T_air, 'K'),
        kp('thermal_margin_K', 'Margin to the allowable skin temperature', margin, 'K', margin > 10 ? 'ok' : margin > 0 ? 'warn' : 'bad'),
        kp('h_wind_Wm2K', 'External convection coefficient', s.h, 'W/m²/K'),
      ],
      plots: [
        { type: 'line', title: 'Hot-soak skin temperature versus solar absorptivity', xlabel: 'Solar absorptivity [-]', ylabel: 'Skin temperature [K]', series: [{ name: `Wind ${i.wind.toFixed(1)} m/s`, x: al, y: al.map((a) => soak(i, a).T) }, { name: 'Calm', x: al, y: al.map((a) => soak(i, a, 0).T) }], annotations: [{ y: i.T_limit, label: 'Allowable' }, { x: i.alpha_s, label: 'Finish' }] },
        { type: 'line', title: 'Hot-soak skin temperature versus wind speed', xlabel: 'Wind speed [m/s]', ylabel: 'Skin temperature [K]', series: [{ name: 'Sunlit skin', x: ws, y: ws.map((w) => soak(i, i.alpha_s, w).T) }] },
        { type: 'line', title: 'Heat radiated to surface 2 versus its emissivity', xlabel: 'Emissivity of surface 2 [-]', ylabel: 'Heat absorbed [W]', series: [{ name: 'Absorbed by surface 2', x: e2, y: e2.map((e) => -radSolve({ ...i, eps2: e }).q2) }] },
      ],
      warnings, models: ['Radiosity network for a three-surface grey diffuse enclosure', 'Exact view factor for coaxial parallel disks', 'Jürges/McAdams wind convection h = 5.7 + 3.8·V (empirical)', 'Steady solar/convective/radiative skin balance'],
      assumptions: ['Grey, diffuse, opaque isothermal surfaces; non-participating gas', 'Facing surfaces idealised as coaxial disks of equivalent area', 'Hot soak is steady with a single effective sky temperature', 'No ground-reflected solar radiation'],
    };
  },
  verify() {
    const b = { T1: 800, eps1: 0.7, r1: 100, T2: 400, eps2: 0.4, r2: 100, gap: 0.01, surr: SURR[0], T3: 300, G_solar: 1000, alpha_s: 0.5, eps_ir: 0.9, T_air: 300, wind: 0, T_skyr: 300, U_in: 0, T_in: 300, T_limit: 400 };
    const r = radSolve(b), bb = radSolve({ ...b, r1: 1, r2: 1, gap: 1, eps1: 1 - 1e-12, eps2: 1 - 1e-12 }), s = soak(b);
    return [
      N.check('View factor of equal coaxial disks at r = L', diskF(1, 1, 1), 0.5 * (3 - Math.sqrt(5)), 1e-12, 'Closed form (Howell catalogue C-41): 0.38197'),
      N.check('Close parallel plates q = σA(T1⁴−T2⁴)/(1/ε1+1/ε2−1)', r.q12, (SIG * r.A1 * (800 ** 4 - 400 ** 4)) / (1 / 0.7 + 1 / 0.4 - 1), 1e-3, 'Infinite parallel-plate limit'),
      N.check('Black disks: q12 = A1·F12·σ(T1⁴−T2⁴)', bb.q12, Math.PI * diskF(1, 1, 1) * SIG * (800 ** 4 - 400 ** 4), 1e-8, 'Black-body exchange'),
      N.check('Enclosure energy conservation', (r.q1 + r.q2 + r.q3) / r.q1, 0, 1e-12, 'Σq = 0'),
      N.check('Hot-soak balance residual', (0.5 * 1000 + 5.7 * (300 - s.T) + 0.9 * SIG * (300 ** 4 - s.T ** 4)) / 500, 0, 1e-8, 'Steady skin energy balance'),
    ];
  },
  calibration: { params: [{ key: 'alpha_s', min: 0.05, max: 1 }, { key: 'T_skyr', min: 200, max: 320 }], sweep: 'G_solar', target: 'T_soak_K', note: 'Supply measured skin temperatures against measured solar irradiance from a ramp soak test.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.thermal_margin_K < 0) out.push({ severity: 'warn', title: 'Ramp hot soak exceeds the skin temperature limit', detail: `Sunlit skin reaches ${(o.T_soak_K - 273.15).toFixed(0)} °C.`, action: 'Use a lower-absorptivity finish, restrict dark liveries on composite structure, or apply a hot-day structural temperature allowance.', basis: 'Steady solar balance' });
    if (i.alpha_s > 0.5) out.push({ severity: 'advise', title: 'Dark finish raises the cooling demand', detail: `The skin runs ${o.soak_rise_K.toFixed(0)} K above ambient; a white finish would cut this substantially (see chart).`, action: 'A light finish lowers ground cooling energy (APU fuel or ground power) and cabin pull-down time.', basis: 'Absorbed solar flux α·G' });
    if (-o.q2_W > 0 && i.eps2 > 0.2) out.push({ severity: 'info', title: 'A low-emissivity shield would cut radiant heating', detail: `Surface 2 absorbs ${(-o.q2_W).toFixed(0)} W at ε = ${i.eps2.toFixed(2)}.`, action: 'A polished or gold-coated heat shield (ε ≈ 0.05–0.1) reduces the absorbed heat roughly in proportion; see the emissivity chart.', basis: 'Radiosity network' });
    return out;
  },
};

export default {
  id: 'thermal', n: 12,
  tagline: 'How hot or cold each part gets, how much heat must be moved, and whether materials, electronics and coolers have margin.',
  analyses: [wall, plate2d, aeroheat, hx, network, stress, radiation],
  consumes: [
    { from: 'propulsion', keys: ['heat_rejection_W', 'fuel_flow_kgs'], why: 'Oil-cooler heat load and fuel flow available as heat sink' },
    { from: 'performance', keys: ['fuel_flow_cruise_kgs'], why: 'Cruise fuel flow for the fuel–oil cooler when the propulsion suite has not run' },
    { from: 'electrical', keys: ['elec_losses_W'], why: 'Electrical losses to be removed by the cooling system' },
  ],
  provides: [
    { key: 'T_max_K', label: 'Peak component temperature', unit: 'K' }, { key: 'q_max_Wm2', label: 'Peak heat flux', unit: 'W/m²' }, { key: 'T_recovery_K', label: 'Recovery temperature', unit: 'K' },
    { key: 'hx_effectiveness', label: 'Heat-exchanger effectiveness', unit: '-' }, { key: 'thermal_margin_K', label: 'Thermal margin', unit: 'K' },
  ],
  handoff: [
    { model: '3-D conjugate heat transfer (coupled CFD–thermal)', why: 'Needs a resolved 3-D flow and solid mesh of the real geometry; here convection enters through correlations and user coefficients', tool: 'CHT-capable CFD solver (finite volume) with imported CAD' },
    { model: '3-D finite-element thermal and thermo-structural models', why: 'Only 1-D wall and 2-D thin-plate conduction are solved natively; thermal stress uses closed-form restraint models', tool: 'General-purpose FE solver with thermal–structural coupling' },
    { model: 'Turbine blade internal and film cooling', why: 'Requires engine gas-path data, coolant passage geometry and film-effectiveness correlations specific to the blade; a steady resistance estimate can be made with the wall analysis using gas-side and coolant-side coefficients', tool: 'Engine secondary-air and blade cooling design codes, CHT CFD' },
    { model: 'Participating-media and spectral radiation (radiative transfer equation in gases)', why: 'Grey diffuse surface exchange only; no gas absorption, emission or scattering', tool: 'Discrete-ordinates or Monte Carlo ray-tracing radiation solver' },
    { model: 'Battery electrochemical–thermal coupling and thermal runaway propagation', why: 'Heat generation is prescribed; cell electrochemistry is handled at equivalent-circuit level in Suite 19', tool: 'Electrochemical battery models and abuse-test data' },
    { model: 'Engine gas-path–solid heat transfer coupling', why: 'Engine cycle temperatures come from Suite 7; component metal temperatures need detailed secondary-air and cooling models', tool: 'Engine thermal models / whole-engine FE' },
  ],
};
