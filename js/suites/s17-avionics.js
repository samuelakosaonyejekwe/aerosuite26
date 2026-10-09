// Suite 17 — Avionics, Navigation and Sensor Systems.
// Inertial navigation error growth, GNSS/INS Kalman filtering, satellite geometry and integrity, radar range and
// Doppler, radio link budgets and data-bus loading, air-data computation, attitude estimation and fault detection.

import * as N from '../core/numerics.js';
import { isa, G0, R_AIR, GAMMA, P0, A0, pressureAltitude, casFromMach } from '../core/atmosphere.js';

// ---- shared helpers --------------------------------------------------------------------------
const num = (key, label, unit, def, min, max, group, help, x) => ({ key, label, unit, default: def, min, max, group, ...(help ? { help } : {}), ...x });
const sel = (key, label, options, def, group, help) => ({ key, label, type: 'select', options, default: def, group, ...(help ? { help } : {}) });
const kp = (key, label, value, unit, status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const thin = (a, n = 300) => { if (a.length <= n) return a; const s = (a.length - 1) / (n - 1); return N.range(n, (j) => a[Math.round(j * s)]); };
const RE = 6371000, C0 = 299792458, KB = 1.380649e-23, D2R = Math.PI / 180, DPH = D2R / 3600, SQH = 60; // °/h → rad/s; √h → √s
const dB = (x) => 10 * Math.log10(x), lin = (d) => 10 ** (d / 10);
/** Complementary error function with fractional error below 1.2e-7 everywhere (Chebyshev fit, Numerical Recipes erfcc). */
export function erfc(x) {
  const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}
const flightOf = (c) => { const alt = c.mission.cruise_alt_m || c.atm.alt_m || 0, V = c.mission.cruise_V_ms || c.flight.V_ms || 30; return { alt, V, M: V / isa(alt, c.atm.dISA_K || 0).a }; };
const sym = (P) => { const n = P.length; for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) P[a][b] = P[b][a] = 0.5 * (P[a][b] + P[b][a]); return P; };
/** Kalman measurement update. Returns the updated state, covariance, innovation and its covariance. */
export function kfUpdate(x, P, z, H, R) {
  const Ht = N.transpose(H), S = N.madd(N.matmul(N.matmul(H, P), Ht), R), K = N.matmul(N.matmul(P, Ht), N.inv(S)), nu = z.map((v, a) => v - N.dot(H[a], x));
  const xn = x.map((v, a) => v + N.dot(K[a], nu)), KH = N.matmul(K, H), P1 = N.matmul(N.eye(x.length).map((r, a) => r.map((v, b) => v - KH[a][b])), P);
  return { x: xn, P: sym(P1), nu, S };
}

// ---- 1. inertial navigation error growth -------------------------------------------------------
// Typical order-of-magnitude sensor grades: [gyro bias °/h, accelerometer bias µg, angle random walk °/√h, velocity random walk m/s/√h]
const GRADE = { 'Consumer MEMS': [30, 5000, 0.5, 0.1], 'Tactical': [1, 1000, 0.1, 0.03], 'Navigation': [0.01, 50, 0.002, 0.003] }, GRADES = [...Object.keys(GRADE), 'Custom'];
const gradeOf = (i) => { const g = GRADE[i.grade] || [i.bg_dph, i.ba_ug, i.arw, i.vrw]; return { bg: g[0] * DPH, ba: g[1] * 1e-6 * G0, arw: (g[2] * D2R) / SQH, vrw: g[3] / SQH }; };
/** Covariance of the single-channel Schuler error model with states [δp, δv, tilt, accel bias, gyro bias]. */
function insCov(s, i, T, nt, part = 'all') {
  const on = (k) => part === 'all' || part === k, n = 5, P0m = N.zeros(n);
  if (on('init')) { P0m[0][0] = i.p0 ** 2; P0m[1][1] = i.v0 ** 2; P0m[2][2] = (i.tilt0_mrad * 1e-3) ** 2; }
  if (on('ba')) P0m[3][3] = s.ba ** 2; if (on('bg')) P0m[4][4] = s.bg ** 2;
  const qv = on('rw') ? s.vrw ** 2 : 0, qf = on('rw') ? s.arw ** 2 : 0, F = [[0, 1, 0, 0, 0], [0, 0, -G0, 1, 0], [0, 1 / RE, 0, 0, 1], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]];
  const f = (t, y) => { const out = new Array(25).fill(0); for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) { let v = 0; for (let k = 0; k < n; k++) v += F[a][k] * y[k * n + b] + y[a * n + k] * F[b][k]; out[a * n + b] = v; } out[6] += qv; out[12] += qf; return out; };
  const r = N.rk4(f, 0, P0m.flat(), T, nt);
  return { t: r.t, sp: r.y.map((y) => Math.sqrt(Math.max(0, y[0]))), sv: r.y.map((y) => Math.sqrt(Math.max(0, y[6]))), sf: r.y.map((y) => Math.sqrt(Math.max(0, y[12]))) };
}
const ins = {
  id: 'ins', title: 'Inertial navigation error growth (Schuler channel)', fidelity: 'reduced-order',
  summary: 'How position error grows when navigating on inertial sensors alone: gyro and accelerometer bias, random walk and initial misalignment propagated through the Schuler-tuned error dynamics for consumer, tactical and navigation-grade sensors.',
  equations: ['Inertial navigation equations', 'Strapdown navigation equations', 'Navigation error propagation equations'],
  inputs: [
    sel('grade', 'Sensor grade', GRADES, 'Tactical', 'Sensors', 'Preset values are typical orders of magnitude, not a specific product'),
    num('bg_dph', 'Gyro bias (custom)', '°/h', 1, 0, 1000, 'Sensors'), num('ba_ug', 'Accelerometer bias (custom)', 'µg', 1000, 0, 1e5, 'Sensors'), num('arw', 'Angle random walk (custom)', '°/√h', 0.1, 0, 10, 'Sensors'), num('vrw', 'Velocity random walk (custom)', 'm/s/√h', 0.03, 0, 5, 'Sensors'),
    num('tilt0_mrad', 'Initial alignment (tilt) error', 'mrad', 0.5, 0, 100, 'Initial state', 'Gyro-compass alignment ≈ 0.05–1 mrad; coarse levelling of MEMS ≈ 5–20 mrad'), num('p0', 'Initial position error', 'm', 3, 0, 1e4, 'Initial state'), num('v0', 'Initial velocity error', 'm/s', 0.05, 0, 50, 'Initial state'),
    num('t_end', 'Unaided navigation time', 's', 7200, 10, 86400, 'Numerics'), num('nSteps', 'Time steps', '', 600, 20, 100000, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => ({ grade: c.mass.mtow_kg > 5700 ? 'Navigation' : c.mass.mtow_kg > 600 ? 'Tactical' : 'Consumer MEMS', tilt0_mrad: c.mass.mtow_kg > 5700 ? 0.1 : c.mass.mtow_kg > 600 ? 1 : 10, t_end: N.clamp(((c.mission.range_km || 50) * 1e3) / Math.max(5, c.mission.cruise_V_ms || c.flight.V_ms || 20), 600, 6 * 3600) }),
  run(i) {
    const s = gradeOf(i), nt = Math.round(i.nSteps), r = insCov(s, i, i.t_end, nt), last = r.t.length - 1, ws = Math.sqrt(G0 / RE), Ts = (2 * Math.PI) / ws, at = (T) => { const q = insCov(s, i, T, Math.max(200, Math.round((nt * T) / i.t_end))); return q.sp[q.sp.length - 1]; };
    const d1h = at(3600), d60 = at(60), parts = ['bg', 'ba', 'rw', 'init'].map((k) => insCov(s, i, i.t_end, nt, k)), tm = r.t.map((v) => v / 60), warnings = [];
    const gr = Object.keys(GRADE).map((gname) => insCov(gradeOf({ grade: gname }), { ...i, tilt0_mrad: { 'Consumer MEMS': 10, 'Tactical': 1, 'Navigation': 0.1 }[gname] }, i.t_end, nt));
    if (i.t_end > 4 * 3600) warnings.push('Beyond a few hours the 24-hour (Earth-rate) and Foucault modes, heading error coupling and vertical-channel instability matter; the single-channel model then understates the error.');
    if (N.amax(r.sf) > 0.1) warnings.push('Tilt error grows beyond about 0.1 rad: the small-angle linear error model is no longer valid and the figures only indicate that unaided navigation is impossible for this long.');
    if (d1h > 1e5) warnings.push('Unaided drift exceeds 100 km in the first hour: this sensor grade is only usable with continuous external aiding.');
    return {
      kpis: [
        kp('nav_drift_m_h', 'Position error after 1 h unaided (1σ)', d1h, 'm', d1h < 1852 * 2 ? 'ok' : 'warn', `${(d1h / 1852).toFixed(2)} NM/h`), kp('pos_sigma_end_m', 'Position error at the end (1σ)', r.sp[last], 'm'),
        kp('pos_60s_m', 'Position error after a 60 s aiding outage (1σ)', d60, 'm'), kp('vel_sigma_end_ms', 'Velocity error at the end (1σ)', r.sv[last], 'm/s'), kp('tilt_sigma_end_mrad', 'Tilt error at the end (1σ)', 1e3 * r.sf[last], 'mrad'),
        kp('schuler_period_min', 'Schuler period', Ts / 60, 'min'), kp('gyro_bias_dph', 'Gyro bias used', s.bg / DPH, '°/h'), kp('accel_bias_ug', 'Accelerometer bias used', (s.ba / G0) * 1e6, 'µg'),
        kp('drift_gyro_m_h', 'Mean drift rate from gyro bias (b·R)', s.bg * RE * 3600, 'm/h'), kp('pos_peak_accel_m', 'Peak position error from accelerometer bias (2b/ωs²)', (2 * s.ba) / (ws * ws), 'm'),
      ],
      plots: [
        { type: 'line', title: 'Unaided position error by sensor grade', xlabel: 'Time [min]', ylabel: 'Position error 1σ [m]', ylog: true, series: Object.keys(GRADE).map((gname, k) => ({ name: gname, x: thin(tm).slice(1), y: thin(gr[k].sp).slice(1) })) },
        { type: 'line', title: 'Error contributions for the selected sensor', xlabel: 'Time [min]', ylabel: 'Position error 1σ [m]', series: [{ name: 'Total', x: thin(tm), y: thin(r.sp) }, { name: 'Gyro bias', x: thin(tm), y: thin(parts[0].sp) }, { name: 'Accelerometer bias', x: thin(tm), y: thin(parts[1].sp) }, { name: 'Random walk', x: thin(tm), y: thin(parts[2].sp) }, { name: 'Initial errors', x: thin(tm), y: thin(parts[3].sp) }], annotations: [{ x: Ts / 60, label: 'Schuler period' }] },
        { type: 'line', title: 'Velocity error', xlabel: 'Time [min]', ylabel: 'Velocity error 1σ [m/s]', series: [{ name: 'Velocity', x: thin(tm), y: thin(r.sv) }] },
      ],
      warnings, models: ['Single-axis Schuler-tuned INS error model with bias states', 'Linear covariance propagation (Lyapunov equation, RK4)', 'Typical sensor-grade error budgets'],
      assumptions: ['One horizontal channel on a non-rotating spherical Earth; channels uncoupled', 'Constant (turn-on) biases with the stated 1σ values; white sensor noise', 'No scale-factor, misalignment or g-sensitivity errors, no vertical channel', 'Stationary or benign trajectory', 'Sensor-grade presets and alignment errors are typical orders of magnitude, not a specific product'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps', levels: [50, 100, 200, 400, 800], metric: 'pos_sigma_end_m' },
  verify() {
    const z = { grade: 'Custom', bg_dph: 0.01, ba_ug: 100, arw: 0.002, vrw: 0, tilt0_mrad: 0, p0: 0, v0: 0, t_end: 5000, nSteps: 2000 }, s = gradeOf(z), ws = Math.sqrt(G0 / RE), T = 5000;
    const g = insCov(s, z, T, 2000, 'bg'), a = insCov(s, z, T, 2000, 'ba'), w = insCov(s, z, 60, 400, 'rw'), t0 = insCov(s, { ...z, tilt0_mrad: 1 }, T, 2000, 'init');
    return [
      N.check('Gyro-bias response b·R·(t − sin(ωs·t)/ωs)', g.sp[2000], s.bg * RE * (T - Math.sin(ws * T) / ws), 1e-8, 'Britting, Inertial Navigation Systems Analysis'),
      N.check('Accelerometer-bias response (b/ωs²)(1 − cos ωs·t)', a.sp[2000], (s.ba / (ws * ws)) * (1 - Math.cos(ws * T)), 1e-5, 'Schuler oscillation (near a null of the response; RK4 with 2000 steps)'),
      N.check('Initial tilt response φ0·R·(1 − cos ωs·t)', t0.sp[2000], 1e-3 * RE * (1 - Math.cos(ws * T)), 1e-5, 'Schuler oscillation (near a null of the response)'),
      N.check('Schuler period 2π·sqrt(R/g) = 84.4 min', (2 * Math.PI) / ws / 60, 84.4, 1e-3, 'Schuler (1923)'),
      N.check('Short-term angle-random-walk growth g·N·sqrt(t⁵/20)', w.sp[400], G0 * s.arw * Math.sqrt(60 ** 5 / 20), 5e-3, 'Triple integral of white noise (Schuler feedback negligible at 60 s)'),
    ];
  },
  calibration: { params: [{ key: 'bg_dph', min: 0, max: 100 }, { key: 'ba_ug', min: 0, max: 20000 }, { key: 'tilt0_mrad', min: 0, max: 50 }], sweep: 't_end', target: 'pos_sigma_end_m', note: 'Supply measured free-inertial position error against time from a static or flight navigation trial (set the grade to Custom).' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.nav_drift_m_h > 2 * 1852) out.push({ severity: 'advise', title: 'Inertial-only navigation drifts quickly', detail: `${(o.nav_drift_m_h / 1852).toFixed(1)} NM in the first hour; ${o.pos_60s_m.toFixed(1)} m after a one-minute outage.`, action: 'Rely on GNSS/INS integration (next analysis) and define the longest tolerable aiding outage; add air-data, magnetometer or vision aiding for unmanned aircraft.', basis: 'Free-inertial error growth' });
    else out.push({ severity: 'info', title: 'Navigation-grade performance', detail: `${(o.nav_drift_m_h / 1852).toFixed(2)} NM/h unaided.`, action: 'Sufficient as an independent means of navigation during GNSS loss; a tactical-grade unit with tight GNSS coupling is lighter, cheaper and uses less power if that independence is not required.', basis: 'Typical long-range inertial performance of about 1–2 NM/h' });
    return out;
  },
};

// ---- 2. GNSS/INS extended Kalman filter ----------------------------------------------------------
function fusionRun(i, seed, noise = 1, aid = true) {
  const u = N.rng(seed), gn = () => noise * N.randn(u), dt = 1 / i.imu_hz, nt = Math.round(i.t_end * i.imu_hz), s = { bg: i.bg_dph * DPH, ba: i.ba_mg * 1e-3 * G0, arw: (i.arw * D2R) / SQH, vrw: i.vrw / SQH }, wt = i.turn_dps * D2R;
  const tTurn = wt > 0 ? Math.PI / wt : 0, cyc = 2 * (i.t_leg + tTurn), wOf = (t) => { const p = t % cyc; return p < i.t_leg ? 0 : p < i.t_leg + tTurn ? wt : p < 2 * i.t_leg + tTurn ? 0 : wt; };
  // truth, and true sensor biases drawn from their 1σ so that the filter statistics can be checked
  let pt = [0, 0], psi = 0; const bT = [s.ba * gn(), s.ba * gn(), s.bg * gn()], sp0 = i.gnss_sigma, sv0 = 0.5, sh0 = 2 * D2R;
  let x = [sp0 * gn(), sp0 * gn(), i.V + sv0 * gn(), sv0 * gn(), sh0 * gn(), 0, 0, 0], P = N.zeros(8); [sp0, sp0, sv0, sv0, sh0, s.ba, s.ba, s.bg].forEach((v, k) => (P[k][k] = Math.max(v * v, 1e-12)));
  const H = [[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0]], R = [[i.gnss_sigma ** 2, 0], [0, i.gnss_sigma ** 2]], every = Math.max(1, Math.round(i.imu_hz / i.gnss_hz)), qb = 1e-10;
  const T = [], tr = [], es = [], err = [], sig = [], nees = [], hErr = []; let nu2 = 0, nUp = 0;
  for (let k = 0; k < nt; k++) {
    const t = k * dt, w = wOf(t), c0 = Math.cos(psi), s0 = Math.sin(psi), p1 = psi + w * dt;
    pt = w === 0 ? [pt[0] + i.V * c0 * dt, pt[1] + i.V * s0 * dt] : [pt[0] + (i.V / w) * (Math.sin(p1) - s0), pt[1] + (i.V / w) * (c0 - Math.cos(p1))]; psi = p1;
    // IMU sample: body specific force [0, V·ω] and yaw rate, with bias and white noise
    const fm = [bT[0] + (s.vrw / Math.sqrt(dt)) * gn(), i.V * w + bT[1] + (s.vrw / Math.sqrt(dt)) * gn()], wm = w + bT[2] + (s.arw / Math.sqrt(dt)) * gn();
    // strapdown mechanisation (mid-point heading)
    const wh = wm - x[7], fb = [fm[0] - x[5], fm[1] - x[6]], pm = x[4] + 0.5 * wh * dt, cm = Math.cos(pm), sm = Math.sin(pm), an = [cm * fb[0] - sm * fb[1], sm * fb[0] + cm * fb[1]];
    x = [x[0] + x[2] * dt + 0.5 * an[0] * dt * dt, x[1] + x[3] * dt + 0.5 * an[1] * dt * dt, x[2] + an[0] * dt, x[3] + an[1] * dt, x[4] + wh * dt, x[5], x[6], x[7]];
    // covariance propagation with the linearised transition matrix Φ = I + F·dt
    const Phi = N.eye(8); Phi[0][2] = Phi[1][3] = dt; Phi[2][4] = (-sm * fb[0] - cm * fb[1]) * dt; Phi[3][4] = (cm * fb[0] - sm * fb[1]) * dt; Phi[2][5] = -cm * dt; Phi[2][6] = sm * dt; Phi[3][5] = -sm * dt; Phi[3][6] = -cm * dt; Phi[4][7] = -dt;
    P = N.matmul(N.matmul(Phi, P), N.transpose(Phi)); P[2][2] += s.vrw ** 2 * dt; P[3][3] += s.vrw ** 2 * dt; P[4][4] += s.arw ** 2 * dt; P[5][5] += qb * dt; P[6][6] += qb * dt; P[7][7] += qb * 1e-4 * dt; sym(P);
    const tn = t + dt, out = tn >= i.t_out && tn < i.t_out + i.d_out;
    if ((k + 1) % every === 0) {
      if (aid && !out) { const z = [pt[0] + i.gnss_sigma * gn(), pt[1] + i.gnss_sigma * gn()], up = kfUpdate(x, P, z, H, R); x = up.x; P = up.P; const Si = N.inv(up.S); nu2 += up.nu[0] * (Si[0][0] * up.nu[0] + Si[0][1] * up.nu[1]) + up.nu[1] * (Si[1][0] * up.nu[0] + Si[1][1] * up.nu[1]); nUp++; }
      const e = [x[0] - pt[0], x[1] - pt[1]], dtm = P[0][0] * P[1][1] - P[0][1] * P[1][0];
      T.push(tn); tr.push(pt.slice()); es.push([x[0], x[1]]); err.push(Math.hypot(e[0], e[1])); sig.push(Math.sqrt(P[0][0] + P[1][1])); hErr.push(x[4] - psi);
      nees.push((P[1][1] * e[0] * e[0] - 2 * P[0][1] * e[0] * e[1] + P[0][0] * e[1] * e[1]) / dtm);
    }
  }
  return { T, tr, es, err, sig, nees, hErr, nis: nUp ? nu2 / nUp : NaN, bErr: [x[5] - bT[0], x[6] - bT[1], x[7] - bT[2]], dist: i.V * i.t_end };
}
const fusion = {
  id: 'fusion', title: 'GNSS/INS integration with an extended Kalman filter', fidelity: 'numerical',
  summary: 'Simulates a loosely coupled GNSS/inertial navigator on a racetrack trajectory with seeded sensor noise and a GNSS outage: estimated versus true track, error against the filter’s own 3σ bound, and statistical consistency checks.',
  equations: ['Kalman filter equations', 'Extended Kalman filter equations', 'Strapdown navigation equations', 'Measurement observation equations', 'Coordinate transformation equations', 'Navigation error propagation equations'],
  inputs: [
    num('V', 'Ground speed', 'm/s', 70, 1, 400, 'Trajectory'), num('t_leg', 'Straight leg duration', 's', 60, 5, 3600, 'Trajectory'), num('turn_dps', 'Turn rate', '°/s', 3, 0, 30, 'Trajectory', '3 °/s is a standard-rate turn'),
    num('bg_dph', 'Gyro bias (1σ)', '°/h', 10, 0, 1000, 'Inertial sensors'), num('ba_mg', 'Accelerometer bias (1σ)', 'mg', 2, 0, 100, 'Inertial sensors'), num('arw', 'Angle random walk', '°/√h', 0.3, 0, 10, 'Inertial sensors'), num('vrw', 'Velocity random walk', 'm/s/√h', 0.1, 0, 5, 'Inertial sensors'),
    num('gnss_sigma', 'GNSS horizontal position noise per axis (1σ)', 'm', 2.5, 0.01, 100, 'GNSS', 'From the satellite-geometry analysis when available'), num('gnss_hz', 'GNSS update rate', 'Hz', 1, 0.1, 20, 'GNSS'),
    num('t_out', 'GNSS outage start', 's', 120, 0, 1e5, 'GNSS'), num('d_out', 'GNSS outage duration', 's', 40, 0, 3600, 'GNSS', 'Jamming, masking or antenna shadowing in a turn'),
    num('imu_hz', 'Inertial sample rate', 'Hz', 10, 2, 400, 'Numerics'), num('t_end', 'Simulated time', 's', 240, 20, 3600, 'Numerics'),
    num('nMC', 'Monte Carlo runs', '', 4, 1, 200, 'Numerics', '', { step: 1, discrete: true }), num('seed', 'Random seed', '', 11, 1, 1e6, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c, up) => { const V = flightOf(c).V, small = c.mass.mtow_kg < 600, big = c.mass.mtow_kg > 5700; return { V, turn_dps: N.clamp((G0 * Math.tan(25 * D2R)) / V / D2R, 1.5, 20), bg_dph: big ? 0.05 : small ? 30 : 3, ba_mg: big ? 0.1 : small ? 5 : 1, arw: big ? 0.005 : small ? 0.5 : 0.1, vrw: big ? 0.005 : small ? 0.1 : 0.03, gnss_sigma: up.avionics?.hpe_sigma_axis_m }; },
  run(i) {
    const nmc = Math.round(i.nMC), runs = N.range(nmc, (k) => fusionRun(i, Math.round(i.seed) + 101 * k)), r0 = runs[0], T = r0.T, inOut = (t) => t >= i.t_out && t < i.t_out + i.d_out, warm = Math.min(20, 0.2 * i.t_end), warnings = [];
    const rmsAt = T.map((t, k) => Math.sqrt(N.mean(runs.map((r) => r.err[k] ** 2)))), neesAt = T.map((t, k) => N.mean(runs.map((r) => r.nees[k]))), aided = T.map((t, k) => (!inOut(t) && t > warm && !(t >= i.t_out && t < i.t_out + i.d_out + 5) ? k : -1)).filter((k) => k >= 0);
    const posErr = Math.sqrt(N.mean(aided.map((k) => rmsAt[k] ** 2))), kEnd = T.findIndex((t) => t >= i.t_out + i.d_out - 1e-9) - 1, hasOut = i.d_out > 0 && kEnd > 0 && kEnd < T.length, outErr = hasOut ? rmsAt[kEnd] : 0, out3 = hasOut ? 3 * N.mean(runs.map((r) => r.sig[kEnd])) : 0;
    const neesM = N.mean(aided.map((k) => neesAt[k])), nis = N.mean(runs.map((r) => r.nis)), hRms = Math.sqrt(N.mean(runs.flatMap((r) => aided.map((k) => r.hErr[k] ** 2)))) / D2R, inB = N.mean(runs.map((r) => r.err.filter((e, k) => e <= 3 * r.sig[k]).length / r.err.length));
    const free = fusionRun(i, Math.round(i.seed), 1, false);
    if (neesM > 3.5) warnings.push('Normalised estimation error is well above its expected value of 2: the filter is over-confident (linearisation error or process noise set too low).');
    if (neesM < 0.8) warnings.push('Normalised estimation error is well below 2: the filter is conservative (covariance larger than the actual error).');
    if (hasOut && outErr > 10 * i.gnss_sigma) warnings.push(`Position error grows to ${outErr.toFixed(0)} m (RMS) by the end of the outage: this inertial grade bridges only short GNSS gaps.`);
    if (nmc < 3) warnings.push('Few Monte Carlo runs: the consistency statistics are noisy.');
    return {
      kpis: [
        kp('pos_err_m', 'Horizontal position error with GNSS (RMS)', posErr, 'm'), kp('pos_err_outage_m', 'Position error at the end of the outage (RMS)', outErr, 'm', outErr < 10 * i.gnss_sigma ? 'ok' : 'warn'),
        kp('sigma3_outage_m', 'Filter 3σ bound at the end of the outage', out3, 'm'), kp('heading_rms_deg', 'Heading error (RMS)', hRms, '°'),
        kp('nees_mean', 'Mean normalised position error squared (NEES)', neesM, '-', neesM > 0.8 && neesM < 3.5 ? 'ok' : 'warn', 'Expected value 2 for a consistent filter'), kp('nis_mean', 'Mean normalised innovation squared (NIS)', nis, '-', nis > 1.2 && nis < 3 ? 'ok' : 'warn', 'Expected value 2'),
        kp('within_3sigma', 'Share of samples inside the 3σ bound', inB, '-', inB > 0.95 ? 'ok' : 'warn'), kp('free_inertial_end_m', 'Error at the end without any GNSS (same sensors)', free.err[free.err.length - 1], 'm'),
        kp('gyro_bias_err_dph', 'Gyro-bias estimation error at the end (RMS)', Math.sqrt(N.mean(runs.map((r) => r.bErr[2] ** 2))) / DPH, '°/h'), kp('improvement', 'GNSS noise / fused position error', (i.gnss_sigma * Math.SQRT2) / posErr, '-'),
      ],
      plots: [
        { type: 'line', title: 'True and estimated track (first run)', xlabel: 'East [m]', ylabel: 'North [m]', equalAspect: true, series: [{ name: 'True', x: thin(r0.tr.map((p) => p[1])), y: thin(r0.tr.map((p) => p[0])) }, { name: 'Estimated', x: thin(r0.es.map((p) => p[1])), y: thin(r0.es.map((p) => p[0])), style: 'dash' }] },
        { type: 'line', title: 'Position error and filter 3σ bound', xlabel: 'Time [s]', ylabel: 'Horizontal error [m]', series: [{ name: 'Error (first run)', x: thin(T), y: thin(r0.err) }, { name: 'RMS over the runs', x: thin(T), y: thin(rmsAt) }, { name: '3σ bound (first run)', x: thin(T), y: thin(r0.sig).map((v) => 3 * v), style: 'dash' }], annotations: hasOut ? [{ x: i.t_out, label: 'GNSS lost' }, { x: i.t_out + i.d_out, label: 'GNSS back' }] : [] },
        { type: 'line', title: 'Filter consistency (NEES averaged over the runs)', xlabel: 'Time [s]', ylabel: 'NEES [-]', series: [{ name: 'NEES', x: thin(T), y: thin(neesAt) }], annotations: [{ y: 2, label: 'Expected' }] },
      ],
      warnings, models: ['Planar strapdown mechanisation (heading, velocity, position)', 'Eight-state extended Kalman filter: position, velocity, heading, accelerometer and gyro biases', 'Loosely coupled GNSS position updates', 'Seeded Monte Carlo with biases drawn from their 1σ values'],
      assumptions: ['Horizontal plane only; level flight with coordinated turns', 'White GNSS position errors (real errors are time-correlated)', 'Constant sensor biases; no scale-factor or misalignment errors', 'Outage removes all GNSS measurements simultaneously', 'Inertial sensor error defaults are typical for the aircraft class, not a specific unit'],
    };
  },
  verify() {
    // scalar random walk with direct measurement: steady prior variance P = (Q + sqrt(Q² + 4QR))/2
    let P = [[10]]; const Q = 0.3, R = 2; for (let k = 0; k < 300; k++) P = [[kfUpdate([0], P, [0], [[1]], [[R]]).P[0][0] + Q]];
    const b = { V: 50, t_leg: 30, turn_dps: 6, bg_dph: 0, ba_mg: 0, arw: 0, vrw: 0, gnss_sigma: 2, gnss_hz: 1, t_out: 1e9, d_out: 0, imu_hz: 50, t_end: 120, nMC: 1, seed: 3 }, clean = fusionRun(b, 3, 0, false);
    const mc = N.range(24, (k) => fusionRun({ ...b, bg_dph: 10, ba_mg: 2, arw: 0.3, vrw: 0.1, imu_hz: 10 }, 500 + 7 * k)), nees = N.mean(mc.map((r) => N.mean(r.nees.slice(20)))), nis = N.mean(mc.map((r) => r.nis));
    return [
      N.check('Scalar Kalman filter steady-state variance', P[0][0], (Q + Math.sqrt(Q * Q + 4 * Q * R)) / 2, 1e-10, 'Algebraic Riccati equation'),
      N.check('Noise-free strapdown mechanisation reproduces the trajectory', clean.err[clean.err.length - 1] / clean.dist, 0, 1e-4, 'Position error / distance flown after two turns'),
      N.check('Filter consistency: mean NEES (2 degrees of freedom)', nees, 2, 0.25, 'Bar-Shalom consistency test, 24 Monte Carlo runs'),
      N.check('Filter consistency: mean NIS', nis, 2, 0.15, 'Innovation whiteness / magnitude test'),
    ];
  },
  calibration: { params: [{ key: 'bg_dph', min: 0, max: 500 }, { key: 'ba_mg', min: 0, max: 50 }, { key: 'vrw', min: 0, max: 2 }], sweep: 'd_out', target: 'pos_err_outage_m', note: 'Supply measured position error at the end of GNSS outages of different lengths (against a reference trajectory) to fit the inertial error parameters.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.pos_err_outage_m > 10 * i.gnss_sigma) out.push({ severity: 'warn', title: 'Outage bridging is weak', detail: `${o.pos_err_outage_m.toFixed(0)} m RMS after ${i.d_out.toFixed(0)} s without GNSS.`, action: 'Use a better inertial grade, add velocity aiding (air data, Doppler, optical flow) or tightly coupled integration that keeps using fewer than four satellites.', basis: 'Monte Carlo error at the end of the outage' });
    if (o.nees_mean > 3.5 || o.within_3sigma < 0.9) out.push({ severity: 'warn', title: 'Filter is over-confident', detail: `NEES ${o.nees_mean.toFixed(1)}; ${(100 * o.within_3sigma).toFixed(0)}% of samples inside 3σ.`, action: 'Increase process noise or model the neglected error sources; an optimistic covariance defeats integrity monitoring downstream.', basis: 'NEES consistency test' });
    out.push({ severity: 'info', title: 'Fusion gain', detail: `Fused position error ${o.pos_err_m.toFixed(2)} m RMS against ${(i.gnss_sigma * Math.SQRT2).toFixed(2)} m for GNSS alone; heading observable to ${o.heading_rms_deg.toFixed(2)}° from the manoeuvres.`, action: 'Accurate navigation permits tighter route keeping and continuous-descent approaches, which save fuel and noise; pass pos_err_m to Suite 16 (guidance) and Suite 24 (mission).', basis: 'Kalman filter error statistics' });
    return out;
  },
};

// ---- 3. GNSS geometry, accuracy and integrity -----------------------------------------------------
const MU_E = 3.986004418e14, W_E = 7.2921159e-5, A_GPS = 26559.7e3;
/** Dilution-of-precision factors from line-of-sight unit vectors in east-north-up axes. */
export function dop(units) {
  if (units.length < 4) return { g: Infinity, p: Infinity, h: Infinity, v: Infinity, t: Infinity };
  const H = units.map((u) => [u[0], u[1], u[2], 1]); let Q;
  try { Q = N.inv(N.matmul(N.transpose(H), H)); } catch { return { g: Infinity, p: Infinity, h: Infinity, v: Infinity, t: Infinity }; }
  const q = (k) => Math.max(0, Q[k][k]);
  return { g: Math.sqrt(q(0) + q(1) + q(2) + q(3)), p: Math.sqrt(q(0) + q(1) + q(2)), h: Math.sqrt(q(0) + q(1)), v: Math.sqrt(q(2)), t: Math.sqrt(q(3)) };
}
/** Visible satellites of a nominal Walker 24/6/1, 55° constellation (approximate almanac) at time t [s]. */
function skyAt(i, t) {
  const la = i.lat * D2R, lo = i.lon * D2R, rs = RE + i.alt_m, n = Math.sqrt(MU_E / A_GPS ** 3), inc = 55 * D2R, cl = Math.cos(la), sl = Math.sin(la), co = Math.cos(lo), so = Math.sin(lo);
  const site = [rs * cl * co, rs * cl * so, rs * sl], e = [-so, co, 0], nn = [-sl * co, -sl * so, cl], up = [cl * co, cl * so, sl], vis = [], np = Math.round(i.n_planes), ns = Math.round(i.n_per_plane);
  for (let p = 0; p < np; p++) for (let s = 0; s < ns; s++) {
    const Om = (2 * Math.PI * p) / np - W_E * t, uu = (2 * Math.PI * s) / ns + (2 * Math.PI * p) / (np * ns) + n * t, cu = Math.cos(uu), su = Math.sin(uu), cO = Math.cos(Om), sO = Math.sin(Om);
    const r = [A_GPS * (cu * cO - su * Math.cos(inc) * sO), A_GPS * (cu * sO + su * Math.cos(inc) * cO), A_GPS * su * Math.sin(inc)], d = [r[0] - site[0], r[1] - site[1], r[2] - site[2]], rg = N.norm(d), le = [N.dot(d, e) / rg, N.dot(d, nn) / rg, N.dot(d, up) / rg], el = Math.asin(le[2]) / D2R;
    if (el >= i.mask_deg) vis.push({ u: le, el, az: (Math.atan2(le[0], le[1]) / D2R + 360) % 360 });
  }
  return vis;
}
const uere = (i) => { const iono = i.dual ? 0.1 : i.s_iono * (1 + 0.15 * i.kp_index), parts = [iono, i.s_tropo, i.s_clock, i.s_mp, i.s_noise]; return { parts, total: Math.hypot(...parts) }; };
const gnss = {
  id: 'gnss', title: 'GNSS geometry, position accuracy and integrity availability', fidelity: 'reduced-order',
  summary: 'Satellite visibility and dilution of precision over a day at the operating site from a nominal 24-satellite constellation, the resulting position error budget, and how often enough satellites are in view for receiver integrity monitoring.',
  equations: ['Coordinate transformation equations', 'Measurement observation equations', 'Navigation error propagation equations'],
  inputs: [
    num('lat', 'Site latitude', '°', 45, -90, 90, 'Site', 'Filled from the case site'), num('lon', 'Site longitude', '°', 0, -180, 180, 'Site'), num('alt_m', 'Receiver altitude', 'm', 0, -400, 20000, 'Site'),
    num('mask_deg', 'Elevation mask angle', '°', 5, 0, 45, 'Site', '5° en route; 10–25° in terrain, urban canyons or banked flight'),
    num('n_planes', 'Orbital planes', '', 6, 3, 8, 'Constellation', '', { step: 1, discrete: true }), num('n_per_plane', 'Satellites per plane', '', 4, 2, 8, 'Constellation', '', { step: 1, discrete: true }),
    sel('dual', 'Receiver', ['Single frequency', 'Dual frequency (ionosphere-free)'], 'Single frequency', 'Error budget'),
    num('s_iono', 'Ionospheric range error, quiet conditions (1σ)', 'm', 4, 0, 30, 'Error budget', 'Single-frequency residual after the broadcast model (typical)'), num('kp_index', 'Geomagnetic Kp index', '-', 2, 0, 9, 'Error budget', 'From live space-weather data when available; raises the ionospheric term (heuristic scaling)'),
    num('s_tropo', 'Tropospheric residual (1σ)', 'm', 0.5, 0, 5, 'Error budget'), num('s_clock', 'Satellite clock and ephemeris (1σ)', 'm', 1, 0, 10, 'Error budget'), num('s_mp', 'Multipath (1σ)', 'm', 1, 0, 20, 'Error budget'), num('s_noise', 'Receiver noise (1σ)', 'm', 0.5, 0, 10, 'Error budget'),
    num('pdop_max', 'Largest acceptable PDOP', '-', 6, 1.5, 20, 'Requirement'), num('hal_m', 'Horizontal alert limit', 'm', 556, 5, 10000, 'Requirement', '556 m (0.3 NM) for non-precision approach, 40 m for approaches with vertical guidance'),
    num('t0_h', 'Time of the sky plot', 'h', 0, 0, 24, 'Numerics'), num('nTimes', 'Time samples over 24 h', '', 144, 12, 1440, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => ({ lat: c.site?.lat ?? undefined, lon: c.site?.lon ?? undefined, alt_m: c.site?.elev_m ?? 0, mask_deg: c.meta.type === 'aeroplane' ? 5 : 10 }),
  run(i0) {
    const i = { ...i0, dual: i0.dual !== 'Single frequency' }, nt = Math.round(i.nTimes), ts = N.range(nt, (k) => (86164 * k) / nt), sk = ts.map((t) => skyAt(i, t)), dp = sk.map((v) => dop(v.map((s) => s.u))), cap = (v) => (Number.isFinite(v) ? Math.min(v, 50) : 50);
    const pd = dp.map((d) => cap(d.p)), hd = dp.map((d) => cap(d.h)), vd = dp.map((d) => cap(d.v)), nv = sk.map((v) => v.length), ue = uere(i), now = skyAt(i, i.t0_h * 3600), warnings = [];
    const avail = pd.filter((v) => v <= i.pdop_max).length / nt, fd = nv.filter((v, k) => v >= 5 && pd[k] <= i.pdop_max).length / nt, fde = nv.filter((v, k) => v >= 6 && pd[k] <= i.pdop_max).length / nt, h95 = 2 * ue.total * N.mean(hd), v95 = 1.96 * ue.total * N.mean(vd), hpl = 6 * ue.total * N.quantile(hd, 0.95);
    if (avail < 0.999) warnings.push(`PDOP exceeds ${i.pdop_max} for ${(100 * (1 - avail)).toFixed(1)}% of the day at this mask angle: plan around the outage windows or add a second constellation.`);
    if (fd < 0.99) warnings.push('Fewer than five usable satellites at times: receiver autonomous integrity monitoring (fault detection) is not continuously available.');
    if (i.kp_index >= 5 && !i.dual) warnings.push('Geomagnetic storm conditions (Kp ≥ 5): single-frequency ionospheric errors and scintillation can be much larger than this budget; augmentation or dual-frequency is advisable.');
    if (Math.abs(i.lat) > 70) warnings.push('At high latitude satellites stay low on the horizon: vertical accuracy degrades (high VDOP).');
    const th = ts.map((t) => t / 3600);
    return {
      kpis: [
        kp('pdop_mean', 'Mean PDOP', N.mean(pd), '-', N.mean(pd) < 3 ? 'ok' : 'warn'), kp('pdop_95', '95th-percentile PDOP', N.quantile(pd, 0.95), '-'), kp('hdop_mean', 'Mean HDOP', N.mean(hd), '-'), kp('vdop_mean', 'Mean VDOP', N.mean(vd), '-'),
        kp('n_visible_mean', 'Mean satellites in view', N.mean(nv), ''), kp('n_visible_min', 'Fewest satellites in view', N.amin(nv), '', N.amin(nv) >= 5 ? 'ok' : 'warn'),
        kp('uere_m', 'User-equivalent range error (1σ)', ue.total, 'm'), kp('hpe_95_m', 'Horizontal position error (95%)', h95, 'm', h95 < i.hal_m ? 'ok' : 'bad'), kp('vpe_95_m', 'Vertical position error (95%)', v95, 'm'),
        kp('hpe_sigma_axis_m', 'Horizontal error per axis (1σ)', (ue.total * N.mean(hd)) / Math.SQRT2, 'm'), kp('pdop_availability', `Share of the day with PDOP ≤ ${i.pdop_max}`, avail, '-', avail > 0.999 ? 'ok' : 'warn'),
        kp('raim_fd_availability', 'Fault-detection availability (≥ 5 satellites)', fd, '-', fd > 0.99 ? 'ok' : 'warn', 'Heuristic count-based criterion'), kp('raim_fde_availability', 'Fault detection and exclusion availability (≥ 6 satellites)', fde, '-'),
        kp('hpl_indicative_m', 'Indicative horizontal protection level', hpl, 'm', hpl < i.hal_m ? 'ok' : 'warn', '≈ 6·UERE·HDOP(95%), a rough stand-in for a RAIM computation'),
      ],
      plots: [
        { type: 'line', title: 'Dilution of precision over a day', xlabel: 'Time [h]', ylabel: 'DOP [-]', series: [{ name: 'PDOP', x: th, y: pd }, { name: 'HDOP', x: th, y: hd }, { name: 'VDOP', x: th, y: vd }], annotations: [{ y: i.pdop_max, label: 'PDOP limit' }] },
        { type: 'line', title: 'Satellites above the mask angle', xlabel: 'Time [h]', ylabel: 'Satellites in view [-]', series: [{ name: 'Visible', x: th, y: nv, style: 'step' }], annotations: [{ y: 5, label: 'Fault detection' }] },
        { type: 'polar', title: `Sky plot at ${i.t0_h.toFixed(1)} h (${now.length} satellites)`, rlabel: 'Zenith angle [°]', series: [{ name: 'Satellites', theta_deg: now.map((s) => s.az), r: now.map((s) => 90 - s.el) }] },
        { type: 'bar', title: 'Range error budget', ylabel: 'Error 1σ [m]', categories: ['Ionosphere', 'Troposphere', 'Clock / ephemeris', 'Multipath', 'Receiver noise', 'Total (RSS)'], series: [{ name: 'Range error', y: [...ue.parts, ue.total] }] },
      ],
      warnings, models: ['Nominal Walker 24/6/1 constellation at 55° inclination on circular orbits (approximate almanac)', 'Spherical rotating Earth, east-north-up line-of-sight geometry', 'DOP from (HᵀH)⁻¹', 'Root-sum-square user-equivalent range error; position error = UERE × DOP'],
      assumptions: ['Nominal constellation: real constellations have more satellites and uneven slots, so real DOP is usually better', 'No terrain or airframe masking beyond the elevation mask', 'Equal, uncorrelated range errors on all satellites', 'Kp scaling of the ionospheric term is a heuristic, not an ionospheric model; no augmentation (SBAS/GBAS) credited', 'Range-error terms are typical single-frequency values; the dual-frequency option only removes the ionospheric term (the higher noise and multipath of the ionosphere-free combination are not added)', 'Integrity availability is a satellite-count screen with an indicative protection level, not a RAIM computation'],
    };
  },
  convergence: { param: 'nTimes', label: 'Time samples over the day', levels: [36, 72, 144, 288], metric: 'pdop_mean' },
  verify() {
    const c = [0, 120, 240].map((a) => [Math.cos(a * D2R), Math.sin(a * D2R), 0]), d = dop([[0, 0, 1], ...c]), b = { lat: 40, lon: -75, alt_m: 0, mask_deg: 0, n_planes: 6, n_per_plane: 4 }, v = skyAt(b, 1000);
    return [
      N.check('GDOP of one zenith and three horizon satellites', d.g, Math.sqrt(3), 1e-10, 'Closed form of (HᵀH)⁻¹ for the symmetric geometry'),
      N.check('HDOP of the same geometry', d.h, Math.sqrt(4 / 3), 1e-10, 'Closed form'), N.check('VDOP of the same geometry', d.v, Math.sqrt(4 / 3), 1e-10, 'Closed form'),
      N.check('Orbital period at the GPS semi-major axis', (2 * Math.PI * Math.sqrt(A_GPS ** 3 / MU_E)) / 43082, 1, 5e-4, 'Half a sidereal day (Kepler’s third law)'),
      N.check('Line-of-sight vectors are unit vectors', N.norm(v[0].u), 1, 1e-12, 'Geometry'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.raim_fd_availability < 0.99 || o.pdop_availability < 0.999) out.push({ severity: 'warn', title: 'GNSS integrity is not continuously available', detail: `Fault detection ${(100 * o.raim_fd_availability).toFixed(1)}% of the day; fewest satellites ${o.n_visible_min}.`, action: 'Use a multi-constellation receiver, lower the antenna mask (placement clear of the fin and rotor), or plan operations with a RAIM prediction; keep an inertial or radio-navigation backup.', basis: 'Satellite count and PDOP criteria' });
    if (i.dual === 'Single frequency' && o.uere_m > 3) out.push({ severity: 'advise', title: 'Ionosphere dominates the error budget', detail: `UERE ${o.uere_m.toFixed(1)} m at Kp ${i.kp_index}.`, action: 'A dual-frequency or SBAS-augmented receiver removes most of it and enables approaches with vertical guidance, reducing diversions and go-arounds (fuel, emissions).', basis: 'Range error budget' });
    out.push({ severity: 'info', title: 'Position accuracy', detail: `${o.hpe_95_m.toFixed(1)} m horizontal and ${o.vpe_95_m.toFixed(1)} m vertical (95%).`, action: 'Feed the per-axis 1σ value to the GNSS/INS analysis as measurement noise.', basis: 'UERE × DOP' });
    return out;
  },
};

// ---- 4. radar ---------------------------------------------------------------------------------
/** Albersheim's formula: single-pulse SNR [dB] for detection probability Pd and false-alarm probability Pfa with n non-coherently integrated pulses (non-fluctuating target). */
export function albersheim(Pd, Pfa, n = 1) { const A = Math.log(0.62 / Pfa), B = Math.log(Pd / (1 - Pd)); return -5 * Math.log10(n) + (6.2 + 4.54 / Math.sqrt(n + 0.44)) * Math.log10(A + 0.12 * A * B + 1.7 * B); }
function radarCalc(i) {
  const lam = C0 / i.f_Hz, G = lin(i.G_dB), B = 1 / i.tau, Nn = KB * 290 * B * lin(i.F_dB), gam = i.rain_k * i.rain_mmh ** i.rain_a; // one-way rain attenuation [dB/km]
  const snr = (R, rcs = i.rcs) => dB((i.P_t * G * G * lam * lam * rcs) / ((4 * Math.PI) ** 3 * R ** 4 * Nn * lin(i.L_dB))) - (2 * gam * R) / 1000, req = albersheim(i.Pd, i.Pfa, Math.max(1, Math.round(i.n_pulses)));
  const rmax = (rcs) => (snr(10, rcs) < req ? 10 : snr(5e7, rcs) > req ? 5e7 : N.brent((R) => snr(R, rcs) - req, 10, 5e7, 1e-3)), Rm = rmax(i.rcs), hz = Math.sqrt(2 * (4 / 3) * RE) * (Math.sqrt(Math.max(0, i.h_radar)) + Math.sqrt(Math.max(0, i.h_target)));
  return { lam, G, B, Nn, gam, snr, req, rmax, Rm, hz, Ru: C0 / (2 * i.prf), vb: (lam * i.prf) / 2, fd: (2 * i.V_close) / lam, dR: (C0 * i.tau) / 2, bw: 70 * lam / Math.max(1e-9, Math.sqrt((G * lam * lam) / (0.6 * Math.PI ** 2))) };
}
const radar = {
  id: 'radar', title: 'Radar detection range, Doppler and radio altimeter', fidelity: 'analytical',
  summary: 'Maximum detection range from the radar equation for a given target, with pulse integration and rain attenuation; Doppler shift, blind speeds and range ambiguity; and the resolution of an FMCW radio altimeter.',
  equations: ['Radar range equations', 'Doppler shift equations', 'Antenna radiation equations', 'Electromagnetic wave equations'],
  inputs: [
    num('P_t', 'Peak transmit power', 'W', 150, 0.001, 1e7, 'Radar', 'Solid-state weather radar ≈ 40–150 W; small sense-and-avoid radar ≈ 1–10 W'), num('G_dB', 'Antenna gain', 'dBi', 34, 0, 60, 'Radar'), num('f_Hz', 'Carrier frequency', 'Hz', 9.375e9, 1e8, 1e11, 'Radar', 'X band 9.3–9.5 GHz'),
    num('tau', 'Transmitted pulse length', 's', 1e-6, 1e-9, 1e-2, 'Radar', 'Sets the pulse energy P·τ and the matched-filter bandwidth 1/τ. With pulse compression enter the uncompressed length: detection is then right and the range resolution is finer than shown by the compression ratio'), num('prf', 'Pulse repetition frequency', 'Hz', 1500, 10, 1e6, 'Radar'), num('n_pulses', 'Pulses integrated per dwell', '', 16, 1, 10000, 'Radar', '', { step: 1, discrete: true }),
    num('F_dB', 'Receiver noise figure', 'dB', 4, 0, 20, 'Radar'), num('L_dB', 'System losses', 'dB', 6, 0, 30, 'Radar'),
    num('rcs', 'Target radar cross-section', 'm²', 5, 1e-4, 1e5, 'Target', 'Typical orders: small drone 0.01–0.1, light aircraft 1–3, airliner 20–100 (aspect dependent)'), num('Pd', 'Required detection probability', '-', 0.9, 0.1, 0.999, 'Target'), num('Pfa', 'False-alarm probability', '-', 1e-6, 1e-12, 1e-2, 'Target'),
    num('V_close', 'Closing speed', 'm/s', 250, 0, 2000, 'Target'), num('h_radar', 'Radar altitude', 'm', 10668, 0, 30000, 'Geometry'), num('h_target', 'Target altitude', 'm', 3000, 0, 30000, 'Geometry'),
    num('rain_mmh', 'Rain rate along the path', 'mm/h', 0, 0, 150, 'Weather'), num('rain_k', 'Rain attenuation coefficient k', 'dB/km', 0.012, 0, 1, 'Weather', 'γ = k·R^α one-way; ≈ 0.012 and 1.26 near 10 GHz, ≈ 0.0007 and 1.1 near 3 GHz (ITU-R P.838 type values)'), num('rain_a', 'Rain attenuation exponent α', '-', 1.26, 0.5, 2, 'Weather'),
    num('B_alt', 'Radio-altimeter sweep bandwidth', 'Hz', 150e6, 1e6, 2e9, 'Radio altimeter', '4.2–4.4 GHz band altimeters sweep about 100–150 MHz'), num('T_sweep', 'Radio-altimeter sweep time', 's', 0.01, 1e-5, 1, 'Radio altimeter'), num('h_agl', 'Height above ground', 'm', 300, 0.5, 2500, 'Radio altimeter'),
  ],
  defaults: (c) => { const f = flightOf(c), small = c.mass.mtow_kg < 600; return { h_radar: f.alt, h_target: Math.max(0, f.alt - 1500), V_close: 2 * f.V, P_t: small ? 5 : c.mass.mtow_kg < 5700 ? 40 : 150, G_dB: small ? 24 : c.mass.mtow_kg < 5700 ? 28 : 34, f_Hz: small ? 24e9 : 9.375e9, rcs: small ? 0.5 : 5, rain_k: small ? 0.1 : 0.012, rain_a: small ? 1.06 : 1.26, h_agl: Math.min(300, Math.max(20, f.alt)) }; },
  run(i) {
    const r = radarCalc(i), use = Math.min(r.Rm, r.hz), warnings = [], Rs = N.logspace(Math.max(50, r.Rm / 50), Math.max(200, r.Rm * 3), 60), rc = N.logspace(1e-3, 1e3, 40), dRa = C0 / (2 * i.B_alt), fb = (2 * i.h_agl * i.B_alt) / (C0 * i.T_sweep);
    if (r.Rm > r.Ru) warnings.push(`Detection range exceeds the unambiguous range c/(2·PRF) = ${(r.Ru / 1e3).toFixed(0)} km: second-time-around echoes will be ambiguous; lower the PRF or use PRF staggering.`);
    if (r.Rm > r.hz) warnings.push(`The radar horizon (${(r.hz / 1e3).toFixed(0)} km) limits detection before the power budget does.`);
    if (r.gam > 0 && 2 * r.gam * r.Rm / 1000 > 6) warnings.push('Two-way rain attenuation exceeds 6 dB at the detection range: storms behind heavy rain may be hidden (attenuation shadow).');
    return {
      kpis: [
        kp('radar_range_m', 'Usable detection range', use, 'm', undefined, `${(use / 1852).toFixed(1)} NM`), kp('range_power_m', 'Detection range from the radar equation', r.Rm, 'm'), kp('radar_horizon_m', 'Radar horizon (4/3 Earth)', r.hz, 'm'),
        kp('snr_required_dB', 'Single-pulse SNR required', r.req, 'dB', undefined, `Pd ${i.Pd}, Pfa ${i.Pfa}, ${Math.round(i.n_pulses)} pulses`), kp('snr_at_10km_dB', 'Single-pulse SNR at 10 km', r.snr(1e4), 'dB'),
        kp('range_unambiguous_m', 'Unambiguous range', r.Ru, 'm', r.Rm <= r.Ru ? 'ok' : 'warn'), kp('range_resolution_m', 'Range resolution', r.dR, 'm', undefined, 'c·τ/2 for an unmodulated pulse'), kp('doppler_Hz', 'Doppler shift at the closing speed', r.fd, 'Hz'),
        kp('blind_speed_ms', 'First blind speed', r.vb, 'm/s', undefined, i.V_close <= r.vb ? 'Doppler unambiguous at the closing speed' : 'Closing speed is above it: Doppler is ambiguous (normal for a low-PRF mode; use staggered or higher PRF for velocity measurement)'), kp('wavelength_m', 'Wavelength', r.lam, 'm'), kp('beamwidth_deg', 'Antenna beamwidth, 70·λ/D', r.bw, '°', undefined, 'Aperture diameter from the gain at 60% efficiency'), kp('rain_atten_dBkm', 'One-way rain attenuation', r.gam, 'dB/km'),
        kp('ralt_resolution_m', 'Radio-altimeter range resolution c/(2B)', dRa, 'm'), kp('ralt_beat_Hz', 'Radio-altimeter beat frequency', fb, 'Hz'),
      ],
      plots: [
        { type: 'line', title: 'Signal-to-noise ratio versus range', xlabel: 'Range [km]', ylabel: 'Single-pulse SNR [dB]', xlog: true, series: [{ name: `Target ${i.rcs} m²`, x: Rs.map((v) => v / 1e3), y: Rs.map((R) => r.snr(R)) }], annotations: [{ y: r.req, label: 'Required' }, { x: r.Rm / 1e3, label: 'Detection range' }] },
        { type: 'line', title: 'Detection range versus target cross-section', xlabel: 'Radar cross-section [m²]', ylabel: 'Detection range [km]', xlog: true, ylog: true, series: [{ name: 'Radar equation', x: rc, y: rc.map((s) => r.rmax(s) / 1e3) }], annotations: [{ x: i.rcs, label: 'Target' }] },
      ],
      warnings, models: ['Monostatic pulse radar range equation with thermal noise kT₀BF', 'Albersheim detection formula with non-coherent integration (empirical fit, non-fluctuating target)', 'Power-law rain attenuation γ = k·R^α', 'Doppler, blind-speed and ambiguity relations; FMCW altimeter relations'],
      assumptions: ['Non-fluctuating (Swerling 0) point target in free space; no clutter, multipath or jamming', 'Matched filter with bandwidth 1/τ; all other losses in the stated loss term', 'Uniform rain along the whole path', 'Horizon from 4/3 effective Earth radius, smooth Earth', 'Radar, target cross-section and rain-attenuation defaults are typical orders of magnitude, not a specific radar; the closing speed defaults to twice the cruise speed (head-on traffic at the same speed)'],
    };
  },
  verify() {
    const b = { P_t: 1000, G_dB: 30, f_Hz: 1e10, tau: 1e-6, prf: 1000, n_pulses: 1, F_dB: 3, L_dB: 3, rcs: 1, Pd: 0.9, Pfa: 1e-6, V_close: 100, h_radar: 1000, h_target: 0, rain_mmh: 0, rain_k: 0.012, rain_a: 1.26, B_alt: 150e6, T_sweep: 0.01, h_agl: 300 }, r = radarCalc(b), lam = C0 / 1e10;
    const hand = ((1000 * 1e6 * lam * lam * 1) / ((4 * Math.PI) ** 3 * KB * 290 * 1e6 * lin(3) * lin(3) * lin(r.req))) ** 0.25;
    return [
      N.check('Detection range from the radar equation', r.Rm, hand, 1e-6, 'R⁴ = P·G²·λ²·σ/((4π)³·kT₀BF·L·SNR)'), N.check('Range scales with the fourth root of cross-section', r.rmax(16) / r.rmax(1), 2, 1e-6, 'Radar equation'),
      N.check('Required SNR for Pd = 0.9, Pfa = 1e-6, one pulse', r.req, 13.2, 0.012, 'Skolnik detection curves (≈13.2 dB); Albersheim is accurate to about 0.2 dB'),
      N.check('First blind speed λ·PRF/2', r.vb, (lam * 1000) / 2, 1e-12, 'Doppler ambiguity'), N.check('Radar horizon 4.12·√h km', r.hz, 4121.6 * Math.sqrt(1000), 1e-4, '4/3-Earth geometry'),
    ];
  },
  calibration: { params: [{ key: 'L_dB', min: 0, max: 25 }, { key: 'F_dB', min: 0.5, max: 15 }], sweep: 'rcs', target: 'range_power_m', note: 'Supply measured detection range against calibrated target cross-section from range trials.' },
  recommend(res, i) {
    const o = res.outputs, out = [], tWarn = o.radar_range_m / Math.max(1, i.V_close);
    out.push({ severity: tWarn < 25 ? 'warn' : 'info', title: 'Warning time against the target', detail: `Detection at ${(o.radar_range_m / 1e3).toFixed(1)} km gives ${tWarn.toFixed(0)} s at a ${i.V_close.toFixed(0)} m/s closing speed.`, action: tWarn < 25 ? 'Below the roughly 25–40 s usually wanted for traffic avoidance: raise power-aperture, integrate more pulses, or fuse with cooperative surveillance (ADS-B).' : 'Adequate for avoidance manoeuvres; weather radar range also allows early re-routing around storms, which saves fuel compared with late deviations.', basis: 'Range / closing speed; 25–40 s is a typical avoidance allowance, not a sourced requirement' });
    if (o.range_power_m > o.range_unambiguous_m) out.push({ severity: 'advise', title: 'Range ambiguity', detail: `Unambiguous range ${(o.range_unambiguous_m / 1e3).toFixed(0)} km is shorter than the detection range.`, action: 'Lower the PRF for long-range search and use a medium/high PRF mode for Doppler.', basis: 'c/(2·PRF)' });
    return out;
  },
};

// ---- 5. communication link and data-bus loading --------------------------------------------------
const MOD = ['BPSK', 'QPSK'], A429 = ['High speed (100 kbit/s)', 'Low speed (12.5 kbit/s)'];
const ber = (ebn0) => 0.5 * erfc(Math.sqrt(Math.max(ebn0, 0)));
function linkCalc(i, d = i.dist) {
  const lam = C0 / i.f_Hz, fspl = 20 * Math.log10((4 * Math.PI * d) / lam), prx = dB(i.P_t) + 30 + i.Gt_dB + i.Gr_dB - fspl - i.L_dB, n0 = -174 + i.NF_dB, ebn0 = prx - n0 - dB(i.Rb);
  const req = dB(N.brent((x) => Math.log10(ber(x)) - Math.log10(i.ber_req), 0.01, 200, 1e-12)) + i.L_impl;
  return { lam, fspl, prx, n0, ebn0, req, margin: ebn0 - req - i.fade_dB, ber: ber(lin(ebn0 - i.L_impl)) };
}
/** Rate-monotonic analysis of a task set [[period, wcet], ...]: utilisation, Liu–Layland bound and exact response times. */
export function rmAnalysis(tasks) {
  const ts = tasks.slice().sort((a, b) => a[0] - b[0]), n = ts.length, U = N.sum(ts.map((t) => t[1] / t[0])), bound = n * (2 ** (1 / n) - 1), resp = [];
  ts.forEach((t, k) => { let R = t[1]; for (let it = 0; it < 500; it++) { let Rn = t[1]; for (let j = 0; j < k; j++) Rn += Math.ceil(R / ts[j][0] - 1e-12) * ts[j][1]; if (Math.abs(Rn - R) < 1e-12 || Rn > 50 * t[0]) { R = Rn; break; } R = Rn; } resp.push(R); });
  return { ts, U, bound, resp, ok: resp.every((R, k) => R <= ts[k][0] + 1e-12) };
}
const parseTasks = (s) => String(s || '').split(/[,;]+/).map((p) => p.trim().split(':').map(Number)).filter((p) => p.length === 2 && p[0] > 0 && p[1] > 0);
const link = {
  id: 'link', title: 'Radio link budget, radio horizon and data-bus loading', fidelity: 'analytical',
  summary: 'Margin of a command, telemetry or voice/data radio link from transmit power to bit-error rate, the line-of-sight range, and the loading and timing margins of the on-board data buses and processor.',
  equations: ['Communication link budget equations', 'Antenna radiation equations', 'Electromagnetic wave equations'],
  inputs: [
    num('f_Hz', 'Carrier frequency', 'Hz', 127e6, 1e6, 1e11, 'Radio link', 'VHF comm 118–137 MHz; L band ≈ 1 GHz; unmanned C2 often 900 MHz, 2.4 GHz or 5 GHz (C band)'), num('P_t', 'Transmit power', 'W', 10, 1e-4, 1e4, 'Radio link'),
    num('Gt_dB', 'Transmit antenna gain', 'dBi', 2, -10, 50, 'Radio link'), num('Gr_dB', 'Receive antenna gain', 'dBi', 6, -10, 50, 'Radio link', 'Whip 2–3 dBi, patch or sector 8–14 dBi, tracking dish 20–30 dBi (typical)'), num('L_dB', 'Cable, pointing, polarisation and atmospheric losses', 'dB', 6, 0, 40, 'Radio link'),
    num('NF_dB', 'Receiver noise figure', 'dB', 5, 0, 20, 'Radio link'), num('Rb', 'Data rate', 'bit/s', 31500, 10, 1e9, 'Radio link', 'VHF data link mode 2 carries 31.5 kbit/s; unmanned C2 and video links 0.1–10 Mbit/s'), sel('mod', 'Modulation', MOD, MOD[1], 'Radio link', 'Coherent BPSK and QPSK have the same bit-error rate against Eb/N0'),
    num('ber_req', 'Required bit-error rate', '-', 1e-6, 1e-12, 1e-2, 'Radio link'), num('L_impl', 'Implementation loss', 'dB', 2, 0, 10, 'Radio link'), num('fade_dB', 'Fade margin required', 'dB', 10, 0, 40, 'Radio link', 'Multipath and airframe shadowing: 10–20 dB is common for air-to-ground links'),
    num('dist', 'Link distance', 'm', 100000, 10, 5e6, 'Geometry'), num('h_gnd', 'Ground antenna height', 'm', 10, 0, 5000, 'Geometry'), num('h_air', 'Aircraft altitude above the ground antenna site', 'm', 3000, 0, 30000, 'Geometry'),
    num('n50', 'ARINC 429: labels at 50 Hz', '', 12, 0, 500, 'Data bus', '', { step: 1, discrete: true }), num('n10', 'ARINC 429: labels at 10 Hz', '', 40, 0, 500, 'Data bus', '', { step: 1, discrete: true }), num('n1', 'ARINC 429: labels at 1 Hz', '', 60, 0, 500, 'Data bus', '', { step: 1, discrete: true }),
    sel('a429', 'ARINC 429 speed', A429, A429[0], 'Data bus'),
    num('n_vl', 'Switched-Ethernet virtual links on the port', '', 60, 1, 5000, 'Data bus', '', { step: 1, discrete: true }), num('frame_B', 'Frame size', 'bytes', 300, 64, 1518, 'Data bus'), num('bag_ms', 'Bandwidth allocation gap', 'ms', 8, 1, 128, 'Data bus'), num('link_bps', 'Port speed', 'bit/s', 100e6, 1e6, 1e10, 'Data bus'),
    { key: 'tasks', label: 'Processor tasks as period:execution time [ms]', type: 'text', default: '5:1, 20:4, 50:8, 200:30', group: 'Processor', help: 'Comma-separated list, e.g. 10:2, 40:8. Rate-monotonic priorities are assumed.' },
  ],
  defaults: (c) => { const f = flightOf(c), uav = c.meta.type === 'uav', big = c.mass.mtow_kg > 5700, hz = Math.sqrt(2 * (4 / 3) * RE) * (Math.sqrt(10) + Math.sqrt(Math.max(f.alt, 1))), dist = uav ? Math.min(0.6 * hz, Math.max(2000, (c.mission.range_km || 10) * 500)) : 0.6 * hz; return { h_air: f.alt, dist, f_Hz: uav ? 2.4e9 : 127e6, P_t: uav ? (c.mass.mtow_kg < 20 ? 0.5 : 2) : big ? 25 : 10, Rb: uav ? 1e6 : 31500, Gr_dB: uav ? (dist > 20e3 ? 24 : 12) : 3, n50: big ? 12 : 4, n10: big ? 40 : 12, n1: big ? 60 : 20, n_vl: big ? 60 : 8 }; },
  run(i) {
    const r = linkCalc(i), hz = Math.sqrt(2 * (4 / 3) * RE) * (Math.sqrt(i.h_gnd) + Math.sqrt(i.h_air)), d0 = i.dist * 10 ** (r.margin / 20), use = Math.min(d0, hz), warnings = [];
    // ARINC 429: 32-bit word + 4-bit gap; switched Ethernet: one frame per virtual link every BAG, FIFO worst case behind every other link
    const bps = i.a429 === A429[0] ? 100e3 : 12.5e3, u429 = ((i.n50 * 50 + i.n10 * 10 + i.n1) * 36) / bps, tf = ((i.frame_B + 20) * 8) / i.link_bps, uEth = (i.n_vl * tf) / (i.bag_ms / 1e3), wc = i.n_vl * tf;
    const tasks = parseTasks(i.tasks), rm = tasks.length ? rmAnalysis(tasks) : { U: 0, bound: 1, ok: true, resp: [], ts: [] }, slack = rm.ts.length ? N.amin(rm.ts.map((t, k) => 1 - rm.resp[k] / t[0])) : 1;
    if (r.margin < 0) warnings.push('Link margin is negative at this distance with the required fade margin: the bit-error-rate requirement is not met.');
    if (i.dist > hz) warnings.push(`The link distance exceeds the radio horizon (${(hz / 1e3).toFixed(0)} km): there is no line of sight; use a relay, satellite link or higher altitude.`);
    if (u429 > 1) warnings.push('The ARINC 429 label set does not fit on one bus at this speed: split it across buses or reduce rates.');
    if (uEth > 1) warnings.push('The switched-Ethernet port is over-subscribed: frames will be dropped.');
    if (!tasks.length) warnings.push('No valid tasks were parsed from the task list (use period:time pairs separated by commas).');
    if (!rm.ok) warnings.push('At least one task misses its deadline under rate-monotonic scheduling (exact response-time analysis).');
    const ds = N.logspace(Math.max(100, i.dist / 100), Math.max(i.dist * 5, hz * 1.2), 60), eb = N.linspace(0, 14, 57);
    return {
      kpis: [
        kp('link_margin_dB', 'Link margin after the fade allowance', r.margin, 'dB', r.margin >= 3 ? 'ok' : r.margin >= 0 ? 'warn' : 'bad'), kp('ebn0_dB', 'Received Eb/N0', r.ebn0, 'dB'), kp('ebn0_req_dB', 'Eb/N0 required (incl. implementation loss)', r.req, 'dB'),
        kp('ber', 'Bit-error rate without fading', r.ber, '-'), kp('p_rx_dBm', 'Received power', r.prx, 'dBm'), kp('fspl_dB', 'Free-space path loss', r.fspl, 'dB'),
        kp('radio_horizon_m', 'Radio horizon (4/3 Earth)', hz, 'm', i.dist <= hz ? 'ok' : 'bad'), kp('range_zero_margin_m', 'Distance at zero margin', d0, 'm'), kp('link_range_m', 'Usable link range', use, 'm', undefined, d0 > hz ? 'Limited by line of sight' : 'Limited by the power budget'),
        kp('a429_utilisation', 'ARINC 429 bus utilisation', u429, '-', u429 < 0.7 ? 'ok' : u429 <= 1 ? 'warn' : 'bad', 'Growth margin of 30% is customary'), kp('eth_utilisation', 'Switched-Ethernet port utilisation', uEth, '-', uEth < 0.5 ? 'ok' : uEth <= 1 ? 'warn' : 'bad'),
        kp('eth_latency_wc_s', 'Worst-case port queuing latency', wc, 's', wc < i.bag_ms / 1e3 ? 'ok' : 'warn', 'One frame from every virtual link ahead in the queue'),
        kp('cpu_utilisation', 'Processor utilisation', rm.U, '-', rm.U <= rm.bound ? 'ok' : rm.ok ? 'warn' : 'bad', `Liu–Layland bound ${rm.bound.toFixed(3)} for ${rm.ts.length} tasks`), kp('cpu_schedulable', 'All deadlines met (1 = yes)', rm.ok ? 1 : 0, '-', rm.ok ? 'ok' : 'bad'), kp('cpu_slack_min', 'Smallest deadline slack', slack, '-', slack > 0.2 ? 'ok' : slack >= 0 ? 'warn' : 'bad'),
      ],
      plots: [
        { type: 'line', title: 'Link margin versus distance', xlabel: 'Distance [km]', ylabel: 'Margin [dB]', xlog: true, series: [{ name: 'Margin after fade allowance', x: ds.map((v) => v / 1e3), y: ds.map((d) => linkCalc(i, d).margin) }], annotations: [{ y: 0, label: 'Limit' }, { x: hz / 1e3, label: 'Radio horizon' }, { x: i.dist / 1e3, label: 'Link distance' }] },
        { type: 'line', title: 'Bit-error rate of coherent BPSK/QPSK', xlabel: 'Eb/N0 [dB]', ylabel: 'Bit-error rate [-]', ylog: true, series: [{ name: 'BER = ½·erfc(√(Eb/N0))', x: eb, y: eb.map((e) => ber(lin(e))) }], annotations: [{ y: i.ber_req, label: 'Required' }] },
        { type: 'bar', title: 'Bus and processor loading', ylabel: 'Utilisation [-]', categories: ['ARINC 429 bus', 'Ethernet port', 'Processor'], series: [{ name: 'Utilisation', y: [u429, uEth, rm.U] }] },
      ],
      tables: rm.ts.length ? [{ title: 'Rate-monotonic response times', columns: ['Period [ms]', 'Execution time [ms]', 'Worst-case response [ms]', 'Deadline met'], rows: rm.ts.map((t, k) => [t[0], t[1], rm.resp[k], rm.resp[k] <= t[0] ? 'yes' : 'no']) }] : [],
      warnings, models: ['Friis free-space link budget with thermal noise −174 dBm/Hz + noise figure', 'Coherent BPSK/QPSK bit-error rate ½·erfc(√(Eb/N0))', '4/3-Earth radio horizon', 'ARINC 429 word timing (32 bits + 4-bit gap)', 'Switched-Ethernet (AFDX-style) bandwidth-allocation-gap loading with FIFO worst case', 'Rate-monotonic utilisation bound and exact response-time analysis'],
      assumptions: ['Free-space propagation with a lumped fade margin: no terrain diffraction, ground-reflection nulls or rain cells', 'Uncoded modulation (forward error correction would add coding gain)', 'Generic message sets: equal frame size and gap for all virtual links, one switch port', 'Independent periodic tasks with deadlines equal to periods and no blocking', 'Radio, antenna, bus-load and task-set defaults are generic illustrative values; the default link distance is 60% of the radio horizon (or the mission radius of a small unmanned aircraft)'],
    };
  },
  verify() {
    const b = { f_Hz: 1e9, P_t: 1, Gt_dB: 0, Gr_dB: 0, L_dB: 0, NF_dB: 0, Rb: 1e6, mod: MOD[0], ber_req: 1e-5, L_impl: 0, fade_dB: 0, dist: 1000, h_gnd: 0, h_air: 1000 }, r = linkCalc(b), rm = rmAnalysis([[50, 12], [40, 10], [30, 10]]);
    return [
      N.check('Free-space path loss at 1 km, 1 GHz', r.fspl, 92.45, 1e-4, '20·log10(4πd/λ) = 32.45 + 20·log f[MHz] + 20·log d[km]'), N.check('erfc(1)', erfc(1), 0.15729920705, 2e-7, 'Tabulated value'),
      N.check('Eb/N0 required for BER 1e-5 (BPSK)', r.req, 9.59, 1e-3, 'Proakis, Digital Communications'), N.check('Received power: +30 dBm − FSPL', r.prx, 30 - r.fspl, 1e-12, 'Friis equation'),
      N.check('Liu–Layland bound for three tasks', rm.bound, 3 * (2 ** (1 / 3) - 1), 1e-12, 'Liu & Layland (1973)'), N.check('Exact response time of the lowest-priority task', rm.resp[2], 52, 1e-12, 'Response-time analysis by hand: 12 + 2·10 + 2·10 = 52 > 50'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.link_margin_dB < 3) out.push({ severity: o.link_margin_dB < 0 ? 'critical' : 'warn', title: 'Radio link margin is thin', detail: `${o.link_margin_dB.toFixed(1)} dB at ${(i.dist / 1e3).toFixed(0)} km; usable range about ${(o.link_range_m / 1e3).toFixed(0)} km.`, action: 'Raise antenna gain (directional or tracking ground antenna) before raising power, lower the data rate, or add forward error correction (typically 4–7 dB of coding gain).', basis: 'Link budget' });
    else out.push({ severity: 'info', title: 'Link closes with margin', detail: `${o.link_margin_dB.toFixed(1)} dB spare; range limited by ${o.range_zero_margin_m > o.radio_horizon_m ? 'line of sight' : 'power'} at ${(o.link_range_m / 1e3).toFixed(0)} km.`, action: 'Surplus margin can be traded for lower transmit power (less heat and electrical load, Suite 19) or a higher data rate.', basis: 'Link budget' });
    if (o.a429_utilisation > 0.7 || o.eth_utilisation > 0.5 || o.cpu_slack_min < 0.2) out.push({ severity: o.cpu_schedulable ? 'advise' : 'critical', title: 'Avionics network or processor is heavily loaded', detail: `ARINC 429 ${(100 * o.a429_utilisation).toFixed(0)}%, Ethernet port ${(100 * o.eth_utilisation).toFixed(0)}%, processor ${(100 * o.cpu_utilisation).toFixed(0)}% with ${(100 * o.cpu_slack_min).toFixed(0)}% minimum slack.`, action: 'Keep at least 30% growth margin: split buses, reduce update rates of slow parameters, or move functions to another module.', basis: 'Bus timing and rate-monotonic schedulability' });
    return out;
  },
};

// ---- 6. air-data system --------------------------------------------------------------------------
/** Air-data computation from static pressure, impact pressure and total temperature (subsonic pitot relation). */
export function airData(ps, qc, tat, rec = 1) {
  const M = Math.sqrt(5 * ((qc / ps + 1) ** (2 / 7) - 1)), sat = tat / (1 + 0.2 * rec * M * M);
  return { M, alt: pressureAltitude(ps), cas: A0 * Math.sqrt(5 * ((qc / P0 + 1) ** (2 / 7) - 1)), sat, tas: M * Math.sqrt(GAMMA * R_AIR * sat) };
}
function airMeas(i, es = i.e_static, ep = i.e_pitot, et = i.e_tat, cp = i.cp_static) {
  const a = isa(i.alt_m, i.dISA), qc = a.p * ((1 + 0.2 * i.M ** 2) ** 3.5 - 1), tat = a.T * (1 + 0.2 * i.rec * i.M ** 2), ps = a.p + es + cp * qc, pt = a.p + qc + ep;
  return { a, qc, truth: { M: i.M, alt: pressureAltitude(a.p), cas: casFromMach(i.M, i.alt_m, i.dISA), sat: a.T, tas: i.M * a.a }, m: airData(ps, Math.max(1e-9, pt - ps), tat + et, i.rec) };
}
const airdata = {
  id: 'airdata', title: 'Air-data system: airspeed, Mach and altitude with sensor errors', fidelity: 'analytical',
  summary: 'Computes indicated altitude, calibrated and true airspeed, Mach number and air temperature from pitot, static and temperature sensors, and shows how transducer errors and static-source position error corrupt them.',
  equations: ['Measurement observation equations', 'Navigation error propagation equations', 'Coordinate transformation equations'],
  inputs: [
    num('alt_m', 'True pressure altitude', 'm', 10668, -400, 20000, 'Flight condition'), num('M', 'True Mach number', '-', 0.78, 0.02, 0.99, 'Flight condition'), num('dISA', 'ISA deviation', 'K', 0, -60, 50, 'Flight condition'),
    num('e_static', 'Static transducer error', 'Pa', 15, -2000, 2000, 'Sensor errors', 'Good air-data modules ≈ 10–30 Pa; low-cost sensors 50–200 Pa'), num('cp_static', 'Static-source position error (ΔCp)', '-', 0.005, -0.2, 0.2, 'Sensor errors', 'Static pressure error as a fraction of impact pressure, after calibration'),
    num('e_pitot', 'Pitot (total pressure) error', 'Pa', 15, -2000, 2000, 'Sensor errors'), num('e_tat', 'Total-temperature probe error', 'K', 0.5, -10, 10, 'Sensor errors'), num('rec', 'Temperature probe recovery factor', '-', 0.98, 0.6, 1, 'Sensor errors'),
    num('s_noise', 'Random pressure noise per sensor (1σ)', 'Pa', 5, 0, 500, 'Sensor errors'), num('alt_tol', 'Altitude-keeping error allowed for the air-data system', 'm', 25, 1, 300, 'Requirement'),
    num('nMC', 'Monte Carlo samples', '', 400, 20, 20000, 'Numerics', '', { step: 1, discrete: true }), num('seed', 'Random seed', '', 5, 1, 1e6, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => { const f = flightOf(c), q = c.mass.mtow_kg > 5700 ? 1 : c.mass.mtow_kg > 600 ? 3 : 8; return { alt_m: f.alt, M: N.clamp(f.M, 0.02, 0.95), dISA: c.atm.dISA_K, e_static: 15 * q, e_pitot: 15 * q, s_noise: 5 * q, cp_static: 0.005 * Math.min(q, 4) }; },
  run(i) {
    const r = airMeas(i), t = r.truth, m = r.m, u = N.rng(Math.round(i.seed)), n = Math.round(i.nMC), da = [], dc = [], warnings = [];
    for (let k = 0; k < n; k++) { const q = airMeas(i, i.e_static + i.s_noise * N.randn(u), i.e_pitot + i.s_noise * N.randn(u)).m; da.push(q.alt - t.alt); dc.push(q.cas - t.cas); }
    const dhdp = -1 / (r.a.p / (R_AIR * (r.a.T - i.dISA)) * G0), eAlt = m.alt - t.alt, hs = N.linspace(0, 15000, 31), cps = N.linspace(-0.03, 0.03, 31);
    if (Math.abs(eAlt) > i.alt_tol) warnings.push(`Altitude error ${eAlt.toFixed(0)} m exceeds the stated allowance: reduced vertical separation operations need tighter static-source calibration.`);
    if (r.qc < 20 * Math.max(Math.abs(i.e_pitot), i.s_noise, 1)) warnings.push('Impact pressure is small compared with the pressure errors: airspeed is poorly measured at this low speed (use a more sensitive differential sensor).');
    if (i.M > 0.95) warnings.push('Near and above Mach 1 the Rayleigh supersonic pitot relation is needed; only the subsonic relation is implemented.');
    return {
      kpis: [
        kp('alt_indicated_m', 'Indicated pressure altitude', m.alt, 'm'), kp('alt_error_m', 'Altitude error', eAlt, 'm', Math.abs(eAlt) <= i.alt_tol ? 'ok' : 'warn', `Allowance ±${i.alt_tol} m`),
        kp('cas_ms', 'Calibrated airspeed (computed)', m.cas, 'm/s'), kp('cas_error_ms', 'Calibrated airspeed error', m.cas - t.cas, 'm/s'), kp('mach_computed', 'Mach number (computed)', m.M, '-'), kp('mach_error', 'Mach error', m.M - t.M, '-'),
        kp('tas_ms', 'True airspeed (computed)', m.tas, 'm/s'), kp('tas_error_ms', 'True airspeed error', m.tas - t.tas, 'm/s'), kp('sat_error_K', 'Static air temperature error', m.sat - t.sat, 'K'),
        kp('dh_dp_m_Pa', 'Altitude sensitivity to static pressure', dhdp, 'm/Pa'), kp('alt_noise_sigma_m', 'Altitude noise (1σ)', N.std(da), 'm'), kp('cas_noise_sigma_ms', 'Airspeed noise (1σ)', N.std(dc), 'm/s'), kp('qc_Pa', 'Impact pressure', r.qc, 'Pa'),
      ],
      plots: [
        { type: 'line', title: 'Altitude error from the static-pressure error versus altitude', xlabel: 'Altitude [m]', ylabel: 'Altitude error [m]', series: [{ name: 'Transducer + position error', x: hs, y: hs.map((h) => { const q = airMeas({ ...i, alt_m: h }); return q.m.alt - q.truth.alt; }) }], annotations: [{ x: i.alt_m, label: 'Flight altitude' }] },
        { type: 'line', title: 'Mach error versus static-source position error', xlabel: 'Position error ΔCp [-]', ylabel: 'Mach error [-]', series: [{ name: 'Mach', x: cps, y: cps.map((c) => airMeas(i, i.e_static, i.e_pitot, i.e_tat, c).m.M - t.M) }], annotations: [{ x: i.cp_static, label: 'As stated' }, { y: 0, label: '' }] },
        { type: 'line', title: 'Airspeed error versus static-source position error', xlabel: 'Position error ΔCp [-]', ylabel: 'Calibrated airspeed error [m/s]', series: [{ name: 'CAS', x: cps, y: cps.map((c) => { const q = airMeas(i, i.e_static, i.e_pitot, i.e_tat, c).m; return q.cas - t.cas; }) }] },
      ],
      warnings, models: ['Standard-atmosphere pressure altitude', 'Subsonic compressible pitot relation for Mach and calibrated airspeed', 'Total-temperature recovery relation', 'Seeded Monte Carlo of pressure-sensor noise'],
      assumptions: ['Subsonic flight; calorically perfect air', 'Position error proportional to impact pressure (single ΔCp); no angle-of-attack or sideslip dependence', 'No pneumatic lag, probe icing or blockage', 'Probe recovery factor known exactly', 'Sensor error defaults are typical for the equipment class, not a specific air-data unit'],
    };
  },
  verify() {
    const b = { alt_m: 9000, M: 0.7, dISA: 10, e_static: 0, cp_static: 0, e_pitot: 0, e_tat: 0, rec: 0.98, s_noise: 0, alt_tol: 25, nMC: 20, seed: 1 }, r = airMeas(b), a = isa(9000, 0), d = (pressureAltitude(a.p + 5) - pressureAltitude(a.p - 5)) / 10;
    return [
      N.check('Error-free sensors recover Mach', r.m.M, 0.7, 1e-12, 'Inverse of the isentropic pitot relation'), N.check('Error-free sensors recover pressure altitude', r.m.alt, 9000, 1e-9, 'ISA inversion'),
      N.check('Error-free sensors recover true airspeed', r.m.tas, 0.7 * isa(9000, 10).a, 1e-12, 'TAS = M·sqrt(γRT) with the recovery relation'), N.check('Calibrated airspeed agrees with the atmosphere module', r.m.cas, casFromMach(0.7, 9000, 10), 1e-12, 'Same pitot relation at sea-level reference'),
      N.check('Altitude sensitivity dh/dp = −1/(ρg)', d, -1 / (a.rho * G0), 1e-5, 'Hydrostatic equation'),
    ];
  },
  calibration: { params: [{ key: 'cp_static', min: -0.1, max: 0.1 }, { key: 'e_static', min: -500, max: 500 }], sweep: 'M', target: 'alt_error_m', note: 'Supply measured altitude error against Mach number from a trailing-cone or pacer-aircraft static-source calibration.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (Math.abs(o.alt_error_m) > i.alt_tol) out.push({ severity: 'warn', title: 'Altimetry error exceeds the allowance', detail: `${o.alt_error_m.toFixed(0)} m at this condition (${(o.dh_dp_m_Pa).toFixed(2)} m per Pa of static error).`, action: 'Recalibrate the static-source error correction (position error grows with Mach and altitude) or relocate the static ports; apply the correction in the air-data computer.', basis: 'Altimetry system error budget' });
    out.push({ severity: 'info', title: 'Air-data accuracy feeds efficiency', detail: `Mach error ${(o.mach_error * 1e3).toFixed(1)}·10⁻³, true airspeed error ${o.tas_error_ms.toFixed(2)} m/s.`, action: 'Flying a few thousandths of Mach or a few hundred feet away from the optimum costs fuel on every flight: accurate air data lets the flight management system hold the economic speed and level (Suite 24).', basis: 'Sensitivity of cruise fuel to speed and altitude errors' });
    return out;
  },
};

// ---- 7. attitude estimation and fault detection ---------------------------------------------------
const qmul = (a, b) => [a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3], a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2], a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1], a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]];
/** Propagate a body-to-navigation quaternion through a rotation vector ω·dt (exact for constant rate). */
export function qstep(q, w, dt) { const th = N.norm(w) * dt; if (th < 1e-14) return q; const s = Math.sin(th / 2) / (th / dt), r = qmul(q, [Math.cos(th / 2), w[0] * s, w[1] * s, w[2] * s]), nr = Math.hypot(...r); return r.map((v) => v / nr); }
const qEuler = (q) => [Math.atan2(2 * (q[0] * q[1] + q[2] * q[3]), 1 - 2 * (q[1] * q[1] + q[2] * q[2])), Math.asin(N.clamp(2 * (q[0] * q[2] - q[3] * q[1]), -1, 1)), Math.atan2(2 * (q[0] * q[3] + q[1] * q[2]), 1 - 2 * (q[2] * q[2] + q[3] * q[3]))];
const eulerQ = (ph, th, ps) => { const c1 = Math.cos(ph / 2), s1 = Math.sin(ph / 2), c2 = Math.cos(th / 2), s2 = Math.sin(th / 2), c3 = Math.cos(ps / 2), s3 = Math.sin(ps / 2); return [c1 * c2 * c3 + s1 * s2 * s3, s1 * c2 * c3 - c1 * s2 * s3, c1 * s2 * c3 + s1 * c2 * s3, c1 * c2 * s3 - s1 * s2 * c3]; };
/** Third row of the direction-cosine matrix: the navigation down axis expressed in body axes. */
const downB = (q) => [2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[2] * q[3] + q[0] * q[1]), 1 - 2 * (q[1] * q[1] + q[2] * q[2])];
const wrap = (a) => Math.atan2(Math.sin(a), Math.cos(a)), TB = 1; // TB [s]: correlation time of the gyro-bias allowance in the heading filter
function attSim(i, noise = 1, init = null) {
  const u = N.rng(Math.round(i.seed)), gn = () => noise * N.randn(u), dt = 1 / i.rate_hz, nt = Math.round(i.t_end * i.rate_hz), A = i.amp_deg * D2R, w1 = 2 * Math.PI * i.f_man, w2 = 0.7 * w1, yr = i.yaw_dps * D2R;
  const bg = [i.bg_dps, -0.6 * i.bg_dps, 0.8 * i.bg_dps].map((v) => v * D2R * noise), sg = i.sg_dps * D2R, sa = i.sa_g, sm = i.sm_deg * D2R, tau = i.tau_c, al = tau / (tau + dt), thr = N.normInv(1 - i.pfa / 2) ** 2;
  let qg = init || eulerQ(0, 0.5 * A * Math.sin(1), 0), qm = qg.slice(), bm = [0, 0, 0], ce = init ? [0, 0] : [0, 0.5 * A * Math.sin(1)], hk = 0, Pk = (5 * D2R) ** 2, hn = 0, Pn = Pk, tDet = NaN, fa = 0;
  const T = [], tru = [], eg = [], em = [], ec = [], hT = [], hF = [], hN = [], hP = [], nis = [];
  for (let k = 0; k < nt; k++) {
    const t = (k + 0.5) * dt, ph = A * Math.sin(w1 * t), th = 0.5 * A * Math.sin(w2 * t + 1), ps = yr * t, dph = A * w1 * Math.cos(w1 * t), dth = 0.5 * A * w2 * Math.cos(w2 * t + 1);
    // body rates from the Euler kinematics; accelerometer senses −g plus a lateral manoeuvre disturbance
    const w = [dph - yr * Math.sin(th), dth * Math.cos(ph) + yr * Math.cos(th) * Math.sin(ph), -dth * Math.sin(ph) + yr * Math.cos(th) * Math.cos(ph)], wm = w.map((v, a) => v + bg[a] + sg * gn());
    const te = (k + 1) * dt, phe = A * Math.sin(w1 * te), the = 0.5 * A * Math.sin(w2 * te + 1), pse = yr * te, dist = i.a_dist_g * Math.sin(2 * Math.PI * 0.05 * te) * noise;
    const f = [Math.sin(the) + sa * gn(), -Math.sin(phe) * Math.cos(the) + dist + sa * gn(), -Math.cos(phe) * Math.cos(the) + sa * gn()], fn = Math.hypot(...f), mag = pse + sm * gn() + (te >= i.t_fault ? i.fault_deg * D2R * noise : 0);
    // (a) gyro integration only
    qg = qstep(qg, wm, dt);
    // (b) Mahony filter: gravity-direction error drives proportional–integral correction; magnetometer corrects yaw
    const dn = downB(qm), a = f.map((v) => -v / fn), e = [a[1] * dn[2] - a[2] * dn[1], a[2] * dn[0] - a[0] * dn[2], a[0] * dn[1] - a[1] * dn[0]], ey = wrap(mag - qEuler(qm)[2]) * (te >= i.t_fault && i.fdi_on && !Number.isNaN(tDet) ? 0 : 1);
    const wc = wm.map((v, b) => v - bm[b] + i.kp_m * e[b] + i.kp_yaw * ey * dn[b]); bm = bm.map((v, b) => v - i.ki_m * e[b] * dt); qm = qstep(qm, wc, dt);
    // (c) Euler-angle complementary filter (roll and pitch)
    const er = [wm[0] + Math.tan(ce[1]) * (wm[1] * Math.sin(ce[0]) + wm[2] * Math.cos(ce[0])), wm[1] * Math.cos(ce[0]) - wm[2] * Math.sin(ce[0])];
    ce = [al * (ce[0] + er[0] * dt) + (1 - al) * Math.atan2(-f[1], -f[2]), al * (ce[1] + er[1] * dt) + (1 - al) * Math.atan2(f[0], Math.hypot(f[1], f[2]))];
    // heading Kalman filters with and without the innovation χ² monitor (magnetometer step fault injected)
    const em2 = qEuler(qm), hd = ((wm[1] * Math.sin(em2[0]) + wm[2] * Math.cos(em2[0])) / Math.cos(em2[1])) * dt, Rm = sm * sm, Q = (sg * dt) ** 2 + (i.bg_dps * D2R) ** 2 * TB * dt; // white gyro noise per sample integrates to (σg·dt)²; the unmodelled constant bias is allowed for as a random walk that matches its drift after TB
    hk += hd; Pk += Q; hn += hd; Pn += Q; const nu = wrap(mag - hk), S = Pk + Rm, d2 = (nu * nu) / S, bad = d2 > thr;
    if (bad && te < i.t_fault) fa++; if (bad && te >= i.t_fault && Number.isNaN(tDet)) tDet = te - i.t_fault;
    if (!(bad && i.fdi_on)) { const K = Pk / S; hk += K * nu; Pk *= 1 - K; } { const K = Pn / (Pn + Rm); hn += K * wrap(mag - hn); Pn *= 1 - K; }
    const ag = qEuler(qg);
    T.push(te); tru.push([phe, the]); eg.push(Math.hypot(wrap(ag[0] - phe), ag[1] - the)); em.push(Math.hypot(wrap(em2[0] - phe), em2[1] - the)); ec.push(Math.hypot(wrap(ce[0] - phe), ce[1] - the)); hT.push(pse); hF.push(wrap(hk - pse)); hN.push(wrap(hn - pse)); hP.push(Pn); nis.push(d2);
  }
  const k0 = Math.round(0.2 * nt), rms = (a) => Math.sqrt(N.mean(a.slice(k0).map((v) => v * v))) / D2R, kf = T.findIndex((t) => t >= i.t_fault), post = (a) => (kf >= 0 && kf < nt - 2 ? N.amax(a.slice(kf).map(Math.abs)) / D2R : 0);
  return { T, tru, eg, em, ec, hF, hN, hP, nis, thr, rg: rms(eg), rm: rms(em), rc: rms(ec), tDet, fa, hFmax: post(hF), hNmax: post(hN), qm, bm, bg };
}
const attitude = {
  id: 'attitude', title: 'Attitude estimation and sensor fault detection', fidelity: 'numerical',
  summary: 'Compares gyro-only integration, a complementary filter and a quaternion (Mahony) filter on simulated inertial data during manoeuvres, and tests an innovation chi-square monitor against an injected magnetometer fault.',
  equations: ['Quaternion propagation equations', 'Coordinate transformation equations', 'Kalman filter equations', 'Measurement observation equations', 'Strapdown navigation equations'],
  inputs: [
    num('amp_deg', 'Roll manoeuvre amplitude', '°', 20, 0, 60, 'Motion', 'Pitch amplitude is half of this'), num('f_man', 'Manoeuvre frequency', 'Hz', 0.1, 0.005, 5, 'Motion'), num('yaw_dps', 'Steady turn rate', '°/s', 1.5, -30, 30, 'Motion'),
    num('a_dist_g', 'Unmodelled lateral acceleration', 'g', 0.03, 0, 1, 'Motion', 'Skidding turns, gusts and vibration corrupt the gravity reference'),
    num('bg_dps', 'Gyro bias', '°/s', 0.1, 0, 5, 'Sensors', 'Consumer MEMS ≈ 0.05–0.5 °/s after warm-up; tactical ≈ 0.0003 °/s'), num('sg_dps', 'Gyro noise per sample (1σ)', '°/s', 0.1, 0, 5, 'Sensors'), num('sa_g', 'Accelerometer noise per sample (1σ)', 'g', 0.01, 0, 0.5, 'Sensors'), num('sm_deg', 'Magnetometer heading noise (1σ)', '°', 1.5, 0.01, 20, 'Sensors'),
    num('tau_c', 'Complementary-filter time constant', 's', 3, 0.05, 100, 'Filters'), num('kp_m', 'Mahony proportional gain', '1/s', 0.4, 0, 50, 'Filters'), num('ki_m', 'Mahony integral gain', '1/s²', 0.03, 0, 20, 'Filters'), num('kp_yaw', 'Heading correction gain', '1/s', 0.5, 0, 20, 'Filters'),
    num('t_fault', 'Magnetometer fault time', 's', 60, 0, 1e5, 'Fault detection'), num('fault_deg', 'Magnetometer fault size (step)', '°', 25, 0, 180, 'Fault detection', 'Hard-iron shift, e.g. a payload or current change'), num('pfa', 'False-alarm probability per test', '-', 0.001, 1e-8, 0.2, 'Fault detection'),
    { key: 'fdi_on', label: 'Reject measurements that fail the monitor', type: 'bool', default: true, group: 'Fault detection' },
    num('rate_hz', 'Sample rate', 'Hz', 100, 5, 2000, 'Numerics'), num('t_end', 'Simulated time', 's', 120, 5, 3600, 'Numerics'), num('seed', 'Random seed', '', 21, 1, 1e6, 'Numerics', '', { step: 1, discrete: true }),
  ],
  defaults: (c) => { const g = c.mass.mtow_kg > 5700 ? 0.02 : c.mass.mtow_kg > 600 ? 0.3 : 1; return { bg_dps: 0.1 * g, sg_dps: 0.1 * Math.max(g, 0.1), sa_g: 0.01 * Math.max(g, 0.2), amp_deg: c.meta.type === 'uav' ? 25 : 20, yaw_dps: N.clamp((G0 * Math.tan(15 * D2R)) / flightOf(c).V / D2R, 0.5, 6) }; },
  run(i) {
    const r = attSim(i), dg = (a) => thin(a).map((v) => v / D2R), warnings = [], best = Math.min(r.rm, r.rc), faulted = i.t_fault < i.t_end && i.fault_deg > 0;
    if (Number.isNaN(r.tDet) && faulted) warnings.push('The injected magnetometer fault was not detected: it is too small relative to the heading uncertainty, or the filter absorbed it. Lower the threshold or add an independent heading source.');
    if (r.fa > 0.02 * r.T.length) warnings.push('The monitor raises frequent false alarms before the fault: the filter noise model is optimistic.');
    if (best > 2) warnings.push('Attitude error exceeds 2° RMS: the gravity reference is disturbed by manoeuvre accelerations; add airspeed/GNSS-velocity compensation of the accelerometer.');
    return {
      kpis: [
        kp('att_rms_mahony_deg', 'Roll/pitch error, quaternion (Mahony) filter (RMS)', r.rm, '°', r.rm < 1 ? 'ok' : r.rm < 2.5 ? 'warn' : 'bad'), kp('att_rms_comp_deg', 'Roll/pitch error, complementary filter (RMS)', r.rc, '°', r.rc < 1 ? 'ok' : r.rc < 2.5 ? 'warn' : 'bad'),
        kp('att_rms_gyro_deg', 'Roll/pitch error, gyro integration only (RMS)', r.rg, '°'), kp('gyro_drift_end_deg', 'Gyro-only error at the end', r.eg[r.eg.length - 1] / D2R, '°'),
        kp('bias_est_dps', 'Gyro bias estimated by the Mahony filter (x axis)', r.bm[0] / D2R, '°/s', undefined, `True ${(r.bg[0] / D2R).toFixed(3)} °/s`),
        kp('fault_detect_s', 'Time to detect the magnetometer fault', Number.isNaN(r.tDet) ? i.t_end : r.tDet, 's', !faulted ? undefined : Number.isNaN(r.tDet) ? 'bad' : r.tDet < 2 ? 'ok' : 'warn', !faulted ? 'No fault is injected within the simulated time' : Number.isNaN(r.tDet) ? 'Not detected' : ''),
        kp('false_alarms', 'False alarms before the fault', r.fa, ''), kp('chi2_threshold', 'Chi-square threshold (1 degree of freedom)', r.thr, '-'),
        kp('heading_err_fdi_deg', 'Largest heading error after the fault, with the monitor', r.hFmax, '°', r.hFmax < 5 ? 'ok' : 'warn'), kp('heading_err_nofdi_deg', 'Largest heading error after the fault, without the monitor', r.hNmax, '°'),
      ],
      plots: [
        { type: 'line', title: 'Roll/pitch estimation error', xlabel: 'Time [s]', ylabel: 'Error [°]', series: [{ name: 'Gyro integration only', x: thin(r.T), y: dg(r.eg) }, { name: 'Complementary filter', x: thin(r.T), y: dg(r.ec) }, { name: 'Quaternion (Mahony) filter', x: thin(r.T), y: dg(r.em) }] },
        { type: 'line', title: 'Heading error around the magnetometer fault', xlabel: 'Time [s]', ylabel: 'Heading error [°]', series: [{ name: i.fdi_on ? 'With innovation monitor' : 'Monitor flags only (measurement still used)', x: thin(r.T), y: dg(r.hF) }, { name: 'No monitor', x: thin(r.T), y: dg(r.hN) }], annotations: [{ x: i.t_fault, label: 'Fault' }] },
        { type: 'line', title: 'Normalised innovation squared of the heading filter', xlabel: 'Time [s]', ylabel: 'ν²/S [-]', ylog: true, series: [{ name: 'Test statistic', x: thin(r.T), y: thin(r.nis).map((v) => Math.max(v, 1e-6)) }], annotations: [{ y: r.thr, label: 'Threshold' }] },
      ],
      warnings, models: ['Quaternion strapdown attitude propagation (exact rotation-vector update)', 'Mahony explicit complementary filter on SO(3) with gyro-bias estimation', 'Euler-angle complementary filter (gyro high-pass, accelerometer low-pass)', 'Scalar heading Kalman filter with innovation chi-square fault monitor', 'Seeded sensor simulation with constant gyro bias and white noise'],
      assumptions: ['Accelerometer measures gravity plus a prescribed lateral disturbance (no centripetal compensation)', 'Magnetometer provides tilt-compensated heading with white noise', 'Constant gyro biases; no scale-factor, misalignment or vibration rectification errors', 'Single simulated realisation (change the seed to see the scatter)', 'Heading filter process noise: integrated white gyro noise plus a random-walk allowance for the unestimated gyro bias (1 s correlation time)', 'Sensor noise, bias and filter-gain defaults are typical values, not a specific unit'],
    };
  },
  verify() {
    const w = [0.3, -0.2, 0.5], q = N.range(1000, () => 0).reduce((a) => qstep(a, w, 0.01), [1, 0, 0, 0]), ang = 2 * Math.acos(N.clamp(q[0], -1, 1)), tot = N.norm(w) * 10;
    const b = { amp_deg: 0, f_man: 0.1, yaw_dps: 0, a_dist_g: 0, bg_dps: 0.5, sg_dps: 0, sa_g: 0, sm_deg: 1, tau_c: 2, kp_m: 2, ki_m: 0.5, kp_yaw: 1, t_fault: 1e9, fault_deg: 0, pfa: 0.05, fdi_on: true, rate_hz: 100, t_end: 120, seed: 2 };
    const c = attSim(b), cv = attSim({ ...b, bg_dps: 0, t_end: 40 }, 0, eulerQ(20 * D2R, -10 * D2R, 0)), e = qEuler(cv.qm);
    // heading filter with a bias-free noisy gyro: the actual error variance must equal the variance the filter reports
    const hc = attSim({ ...b, bg_dps: 0, sg_dps: 2, sm_deg: 0.5, t_end: 300, seed: 9 }), k0 = 2000, hvar = N.mean(hc.hN.slice(k0).map((v) => v * v)) / N.mean(hc.hP.slice(k0));
    return [
      N.check('Quaternion propagation: rotation angle = |ω|·t', Math.abs(ang - 2 * Math.PI * Math.round(tot / (2 * Math.PI))), Math.abs(tot - 2 * Math.PI * Math.round(tot / (2 * Math.PI))), 1e-9, 'Exact rotation-vector integration at constant rate'),
      N.check('Complementary filter steady error with a constant gyro bias = b·τ', c.ec[c.ec.length - 1] / D2R, Math.hypot(0.5 * 2, 0.3 * 2), 2e-3, 'First-order blend: error = bias × time constant'),
      N.check('Mahony filter removes a constant gyro bias (integral action)', c.em[c.em.length - 1] / D2R, 0, 1e-3, 'Zero steady-state attitude error'),
      N.check('Mahony filter converges from a 20° initial error', Math.hypot(e[0], e[1]) / D2R, 0, 1e-3, 'Almost-global convergence of the explicit complementary filter'),
      N.check('Chi-square threshold for 5% false alarms, 1 degree of freedom', c.thr, 3.8415, 1e-3, 'χ² tables'),
      N.check('Heading Kalman filter: actual error variance / reported variance', hvar, 1, 0.15, 'Consistency of a correctly tuned Kalman filter (about 1000 independent samples)'),
    ];
  },
  calibration: { params: [{ key: 'kp_m', min: 0.05, max: 20 }, { key: 'tau_c', min: 0.1, max: 30 }], sweep: 'a_dist_g', target: 'att_rms_mahony_deg', note: 'Supply measured attitude error (against a reference attitude system or motion table) for different manoeuvre disturbance levels to tune the filter gains.' },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (Math.min(o.att_rms_mahony_deg, o.att_rms_comp_deg) > 1) out.push({ severity: 'warn', title: 'Attitude accuracy is marginal', detail: `Best filter error ${Math.min(o.att_rms_mahony_deg, o.att_rms_comp_deg).toFixed(2)}° RMS.`, action: 'Lower the accelerometer correction gain during manoeuvres (gain scheduling on |a| ≠ g), compensate centripetal acceleration with airspeed or GNSS velocity, or use an EKF with a manoeuvre model.', basis: 'Simulated estimation error' });
    if (o.fault_detect_s >= i.t_end && i.fault_deg > 0 && i.t_fault < i.t_end) out.push({ severity: 'warn', title: 'Magnetometer fault goes undetected', detail: `A ${i.fault_deg}° heading step was absorbed without an alarm.`, action: 'Tighten the monitor threshold, cross-check heading against GNSS track in steady flight, and treat the event in the Suite 22 failure analysis.', basis: 'Innovation chi-square test' });
    else if (i.fault_deg > 0 && i.t_fault < i.t_end) out.push({ severity: 'info', title: 'Fault detection works for this case', detail: `Detected ${o.fault_detect_s.toFixed(2)} s after onset; heading error held to ${o.heading_err_fdi_deg.toFixed(1)}° instead of ${o.heading_err_nofdi_deg.toFixed(1)}°.`, action: 'Analytical redundancy of this kind can replace a duplicated sensor, saving mass, power and cost, provided coverage is demonstrated for all credible faults.', basis: 'Innovation monitoring' });
    return out;
  },
};

export default {
  id: 'avionics', n: 17,
  tagline: 'How accurately the aircraft knows where it is and how it is oriented, how far its radar and radios reach, and whether sensors, buses and processors have margin.',
  analyses: [ins, gnss, fusion, radar, link, airdata, attitude],
  consumes: [],
  provides: [
    { key: 'pos_err_m', label: 'Navigation position error', unit: 'm' }, { key: 'nav_drift_m_h', label: 'Inertial drift', unit: 'm/h' }, { key: 'radar_range_m', label: 'Radar detection range', unit: 'm' }, { key: 'link_margin_dB', label: 'Radio link margin', unit: 'dB' },
  ],
  handoff: [
    { model: 'Electromagnetic field solvers for installed antenna patterns, radome effects and radar cross-section', why: 'Needs full-wave or asymptotic solution of Maxwell’s equations on the airframe geometry; antennas are represented here by gain values', tool: 'Method-of-moments / FDTD / physical-optics electromagnetic solvers' },
    { model: 'Unscented Kalman filter and particle filter (nonlinear Bayesian estimation)', why: 'Not implemented natively; the extended Kalman filter and its consistency tests cover the mildly nonlinear GNSS/INS case shown', tool: 'Estimation toolboxes; the EKF analysis here as the baseline for comparison' },
    { model: 'Full three-dimensional strapdown INS with Earth rate, transport rate, vertical channel and tightly coupled GNSS (pseudorange/carrier phase)', why: 'The native models are a single Schuler channel and a planar loosely coupled filter', tool: 'Navigation simulation toolkits with GNSS signal simulators' },
    { model: 'Vision–inertial odometry and radar–inertial navigation', why: 'Require image or radar signal simulation and feature tracking', tool: 'Robotics perception frameworks with sensor simulators' },
    { model: 'Receiver autonomous integrity monitoring with protection levels (weighted least-squares RAIM, ARAIM, SBAS/GBAS)', why: 'Only a satellite-count and DOP heuristic with an indicative protection level is given', tool: 'Certified GNSS availability prediction tools' },
    { model: 'Radar in clutter and fluctuating targets (Swerling 1–4), terrain masking and propagation ducting', why: 'Free-space, non-fluctuating detection only', tool: 'Radar system simulators with terrain and clutter models' },
    { model: 'Network calculus / trajectory approach for certified AFDX worst-case delays; protocol-level and software integration testing', why: 'A single-port FIFO bound and generic message sets are used', tool: 'Avionics network analysis tools and hardware-in-the-loop integration rigs' },
  ],
};
