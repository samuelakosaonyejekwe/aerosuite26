// Explicit nonlinear dynamics finite-element kernel in three dimensions.
// Elements: 8-node hexahedron with one-point quadrature and Flanagan–Belytschko hourglass control; 4-node
// co-rotational thin shell (Belytschko–Lin–Tsay: one in-plane point, Gauss integration of plane-stress plasticity through
// the thickness, so plate bending is carried by the stress field and not by the stabilisation); 2-node co-rotational beam with axial, torsional and bi-axial bending elastic–plastic resultants (plastic hinges at
// the ends); 2-node nonlinear spring/damper (crush element). Lumped mass, central-difference time integration
// with an automatic Courant time step and optional mass scaling. J2 plasticity with Johnson–Cook hardening,
// rate and thermal terms and a Cowper–Symonds option, failure by equivalent plastic strain with erosion.
// Penalty / crush-spring contact with a rigid plane (Coulomb friction, or a water-like pressure law) and with a
// rigid sphere or cylinder impactor. Energy accounting by category with a global balance.
// Pure computation on typed arrays: runs unchanged in a browser, a Web Worker and Node.

const GAM = [1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1];
// Gauss–Legendre points and weights on [−1, 1] for 1 to 5 points (through-thickness integration of the shell)
const GQ = [[[0], [2]], [[-0.5773502691896257, 0.5773502691896257], [1, 1]], [[-0.7745966692414834, 0, 0.7745966692414834], [5 / 9, 8 / 9, 5 / 9]],
  [[-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526], [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]],
  [[-0.906179845938664, -0.5384693101056831, 0, 0.5384693101056831, 0.906179845938664], [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891]]];
const PERM = [0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 0, 5, 6, 7, 4, 2, 3, 0, 1, 6, 7, 4, 5, 3, 0, 1, 2, 7, 4, 5, 6, 4, 7, 6, 5, 0, 3, 2, 1, 5, 4, 7, 6, 1, 0, 3, 2, 6, 5, 4, 7, 2, 1, 0, 3, 7, 6, 5, 4, 3, 2, 1, 0];

/** Flow stress [Pa] of material m at equivalent plastic strain ep, strain rate [1/s] and temperature [K]. */
export function flowStress(m, ep, rate, T) {
  let s = m.A + (m.B ? m.B * Math.pow(ep > 0 ? ep : 0, m.n) : 0);
  if (m.C && rate > m.eps0) s *= 1 + m.C * Math.log(rate / m.eps0);
  if (m.csD) s *= 1 + Math.pow(Math.max(rate, 0) / m.csD, 1 / m.csQ);
  if (m.thermal && T > m.T0) s *= Math.max(0, 1 - Math.pow(Math.min(1, (T - m.T0) / (m.Tm - m.T0)), m.m));
  return s;
}

/**
 * Volume-weighted gradient operator of the 8-node hexahedron: Bm[0..7] = ∫∂N/∂x dV, Bm[8..15] = ∫∂N/∂y dV,
 * Bm[16..23] = ∫∂N/∂z dV (exact for any straight-edged hexahedron). X holds x[0..7], y[8..15], z[16..23].
 * Returns the element volume.
 */
export function hexGradient(X, Bm) {
  for (let I = 0; I < 8; I++) {
    const p = 8 * I, n2 = PERM[p + 1], n3 = PERM[p + 2], n4 = PERM[p + 3], n5 = PERM[p + 4], n6 = PERM[p + 5], n8 = PERM[p + 7];
    for (let c = 0; c < 3; c++) {
      const a = 8 * ((c + 1) % 3), b = 8 * ((c + 2) % 3);
      const b2 = X[b + n2], b3 = X[b + n3], b4 = X[b + n4], b5 = X[b + n5], b6 = X[b + n6], b8 = X[b + n8];
      Bm[8 * c + I] = (X[a + n2] * (b6 - b3 - (b4 - b5)) + X[a + n3] * (b2 - b4) + X[a + n4] * (b3 - b8 - (b5 - b2)) + X[a + n5] * (b8 - b6 - (b2 - b4)) + X[a + n6] * (b5 - b2) + X[a + n8] * (b4 - b5)) / 12;
    }
  }
  let V = 0; for (let I = 0; I < 8; I++) V += X[I] * Bm[I];
  return V;
}

/**
 * Create an explicit FE model. Build it with node/material/hex/shell/beam/spring/mass/fix/load/plane/impactor, call
 * init(), then step() or run(tEnd, onStep). Options: gravity [3], cfl, hourglass (stiffness coefficient),
 * hourglassVisc, shellHourglass (stiffness coefficient of the shell stabilisation), bulkLinear, bulkQuad, penalty (contact stiffness as a fraction of m/Δt²), dtMax, massScaleDt
 * (hexahedra whose stable step is below this value receive added mass, reported in addedMass), damping
 * (mass-proportional, 1/s). The step of beams, springs and crush contacts follows from the assembled nodal masses;
 * nodal rotary inertia is raised where rotation would otherwise control it (reported in addedInertia).
 */
export function createFE(opt = {}) {
  const gv = opt.gravity || [0, 0, 0], cfl = opt.cfl ?? 0.9, kHG = opt.hourglass ?? 0.05, qHG = opt.hourglassVisc ?? 0.02, bq1 = opt.bulkLinear ?? 0.06, bq2 = opt.bulkQuad ?? 1.5, pf = opt.penalty ?? 0.1, rSH = opt.shellHourglass ?? 0.03;
  const bld = { xyz: [], hex: [], hexM: [], shell: [], shM: [], shT: [], shQ: [], beam: [], sec: [], ref: [], spr: [], pm: [], fix: [], load: [], v0: [] };
  const mats = [], planes = [], imps = [];
  const fe = { t: 0, dt: 0, nStep: 0, damping: opt.damping || 0, mats, planes, impactors: imps, addedMass: 0, addedInertia: 0, erodedHex: 0, erodedShell: 0, erodedBeam: 0 };
  const E = fe.energy = { kinetic: 0, internal: 0, plastic: 0, hourglass: 0, viscous: 0, contact: 0, friction: 0, crush: 0, external: 0, damping: 0, initial: 0, total: 0, error: 0, peak: 0 };

  fe.node = (x, y, z) => { bld.xyz.push(x, y, z); return bld.xyz.length / 3 - 1; };
  /** Material: { E, nu, rho, A (yield), B, n, C, eps0, m, Tm, T0, cp, beta, csD, csQ, efail }. B = 0 is perfectly plastic; A = Infinity is elastic. */
  fe.material = (p) => {
    const m = { nu: 0.3, B: 0, n: 1, C: 0, eps0: 1, m: 1, Tm: 0, T0: 293, cp: 900, beta: 0.9, csD: 0, csQ: 5, efail: 0, ...p };
    m.G = m.E / (2 * (1 + m.nu)); m.K = m.E / (3 * (1 - 2 * m.nu)); m.cd = Math.sqrt((m.K + (4 / 3) * m.G) / m.rho); m.cs = Math.sqrt(Math.max(m.K + (4 / 3) * m.G, (opt.dtBulk ?? 1) * 3 * m.K) / m.rho); // m.cs bounds the highest element frequency (volumetric mode of a single element)
    m.thermal = m.Tm > m.T0; // thermal softening and adiabatic heating are active only when a melting temperature is given
    mats.push(m); return mats.length - 1;
  };
  fe.hex = (n, mat) => { bld.hex.push(...n); bld.hexM.push(mat); return bld.hexM.length - 1; };
  /**
   * Thin shell: four nodes in order around the element (the normal follows the right-hand rule), material, thickness and
   * the number of Gauss points through the thickness (1–5; 2 is exact for elastic bending, 5 resolves spreading plasticity).
   * Uses the nodal rotations. A point fails at the material's efail; the element is eroded when all its points have failed.
   */
  fe.shell = (n, mat, t, nq = 5) => { bld.shell.push(...n); bld.shM.push(mat); bld.shT.push(t); bld.shQ.push(Math.min(5, Math.max(1, Math.round(nq)))); return bld.shM.length - 1; };
  /**
   * Beam section: { E, G, rho, A, Iy, Iz, J, Np, Mpy, Mpz, hard, cap, thetaSoft, residual, thetaFail, epsFail }. Hinge capacity is
   * Mp·min(cap, 1 + hard·θp) up to the plastic rotation thetaSoft, then Mp·residual; the element is eroded at thetaFail.
   * ref: vector fixing the local y axis (Iy resists rotation about it).
   */
  fe.beam = (n1, n2, sec, ref = [0, 0, 1]) => { bld.beam.push(n1, n2); bld.sec.push({ Np: Infinity, Mpy: Infinity, Mpz: Infinity, hard: 0, cap: 1, thetaSoft: Infinity, residual: 1, thetaFail: Infinity, epsFail: Infinity, ...sec }); bld.ref.push(ref); return bld.sec.length - 1; };
  /** Spring: law { k, Fyc, Fyt, dmax, kb, c, kLat }; dir = fixed unit direction from n1 to n2 (compression shortens along it) or null for the current line. */
  fe.spring = (n1, n2, law, dir = null) => { bld.spr.push({ n1, n2, dir, law: { Fyc: Infinity, Fyt: Infinity, dmax: Infinity, kb: 0, c: 0, kLat: 0, ...law }, dp: 0, d: 0, F: 0, wp: 0, e: 0, bottomed: false, dMaxSeen: 0, dpMax: 0 }); return bld.spr.length - 1; };
  fe.mass = (node, m, inertia = 0) => { bld.pm.push([node, m, inertia]); };
  fe.fix = (node, mask = [1, 1, 1, 1, 1, 1]) => { bld.fix.push([node, mask]); };
  fe.load = (node, f) => { bld.load.push([node, f]); };
  fe.velocity = (node, v) => { bld.v0.push([node, v]); };
  /**
   * Rigid plane: { p0, n, mu, fluid: { rho, Cd } }. Contact points are added with contact(node, { off, law, area }):
   * off is a body-fixed offset from the node (rotates with it), law an optional crush law { k, Fy, dmax, kb, failAtMax }, mu an optional friction override.
   */
  fe.plane = (p) => { const pl = { p0: [0, 0, 0], n: [0, 0, 1], mu: 0, fluid: null, ...p, pts: [], force: 0 }; pl.contact = (node, c = {}) => { pl.pts.push({ node, off: c.off || null, law: c.law || null, area: c.area || 0, mu: c.mu, name: c.name || '', dp: 0, wp: 0, failed: false, dMax: 0, F: 0, tFail: NaN }); return pl.pts.length - 1; }; planes.push(pl); return pl; };
  /** Rigid impactor: { c [3], v [3], R, mass, axis (unit vector → cylinder), nodes, fcap (force cap per node, N) or fcaps (one cap per listed node), free [3] (0/1 mask) }. An infinite mass gives a prescribed velocity. */
  fe.impactor = (p) => { const im = { axis: null, fcap: Infinity, free: [1, 1, 1], ...p, c: p.c.slice(), v: p.v.slice(), f: [0, 0, 0], force: 0 }; imps.push(im); return im; };

  fe.init = () => {
    const nn = fe.nn = bld.xyz.length / 3, nh = fe.nh = bld.hexM.length, nb = fe.nb = bld.sec.length;
    const x = fe.x = Float64Array.from(bld.xyz), x0 = fe.x0 = Float64Array.from(bld.xyz), v = fe.v = new Float64Array(3 * nn), a = fe.a = new Float64Array(3 * nn), f = new Float64Array(3 * nn), fc = new Float64Array(3 * nn), fl = new Float64Array(3 * nn), fh = new Float64Array(3 * nn), fq = new Float64Array(3 * nn);
    const ms = fe.mass_ = new Float64Array(nn), fixed = new Uint8Array(6 * nn), ns = fe.ns = bld.shM.length, hasRot = nb > 0 || ns > 0 || planes.some((p) => p.pts.some((c) => c.off));
    const w = fe.w = new Float64Array(hasRot ? 3 * nn : 0), mo = new Float64Array(hasRot ? 3 * nn : 0), mc = new Float64Array(hasRot ? 3 * nn : 0), In = fe.inertia = new Float64Array(hasRot ? nn : 0), R = fe.R = new Float64Array(hasRot ? 9 * nn : 0);
    for (let i = 0; i < (hasRot ? nn : 0); i++) R[9 * i] = R[9 * i + 4] = R[9 * i + 8] = 1;
    const hx = fe.hexN = Int32Array.from(bld.hex), hm = Int32Array.from(bld.hexM), sig = fe.sig = new Float64Array(6 * nh), ep = fe.ep = new Float64Array(nh), Te = fe.temp = new Float64Array(nh).fill(293), alive = fe.hexAlive = new Uint8Array(nh).fill(1), hq = new Float64Array(12 * nh), hrho = new Float64Array(nh), hvol = fe.hexVol = new Float64Array(nh), sys = new Float64Array(nh); // sys: rate-free flow stress at the current plastic strain and temperature
    const X = new Float64Array(24), U = new Float64Array(24), Bm = new Float64Array(24), gm = new Float64Array(8);
    let dtE = Infinity; const target = opt.massScaleDt || 0, kT = new Float64Array(nn), kR = new Float64Array(nn); // nodal stiffness sums (Gershgorin bound on the highest frequency)
    // hexahedra: mass and stable step
    for (let e = 0; e < nh; e++) {
      const m = mats[hm[e]];
      for (let I = 0; I < 8; I++) { const k = 3 * hx[8 * e + I]; X[I] = x[k]; X[8 + I] = x[k + 1]; X[16 + I] = x[k + 2]; }
      const V = hexGradient(X, Bm); if (!(V > 0)) throw new Error('explicitfe: hexahedron ' + e + ' has non-positive volume (check node ordering)');
      let gx = 0, gy = 0, gz = 0, gxy = 0, gxz = 0, gyz = 0; for (let I = 0; I < 8; I++) { gx += Bm[I] * Bm[I]; gy += Bm[8 + I] * Bm[8 + I]; gz += Bm[16 + I] * Bm[16 + I]; gxy += Bm[I] * Bm[8 + I]; gxz += Bm[I] * Bm[16 + I]; gyz += Bm[8 + I] * Bm[16 + I]; }
      const l = V / Math.sqrt(2 * Math.max(gx + Math.abs(gxy) + Math.abs(gxz), gy + Math.abs(gxy) + Math.abs(gyz), gz + Math.abs(gxz) + Math.abs(gyz)));
      let rho = m.rho, dt = (cfl * l) / m.cs;
      if (target > dt) { const s = (target / dt) ** 2; fe.addedMass += (s - 1) * rho * V; rho *= s; dt = target; }
      hrho[e] = rho; hvol[e] = V; Te[e] = m.T0; sys[e] = m.A; if (dt < dtE) dtE = dt;
      for (let I = 0; I < 8; I++) ms[hx[8 * e + I]] += (rho * V) / 8;
    }
    // shells: mass, rotary inertia and stable step. The rotary inertia is the larger of the physical t²/12 and A/8 per unit mass,
    // so that the rotational modes do not control the step (the excess is reported in addedInertia).
    const sn = fe.shellN = Int32Array.from(bld.shell), sm = Int32Array.from(bld.shM), sT = Float64Array.from(bld.shT), sQ = Int32Array.from(bld.shQ), sO = new Int32Array(ns + 1); for (let e = 0; e < ns; e++) sO[e + 1] = sO[e] + sQ[e];
    const sSig = fe.shellSig = new Float64Array(3 * sO[ns]), sEp = fe.shellEp = new Float64Array(sO[ns]), sTe = fe.shellTemp = new Float64Array(sO[ns]), sYs = new Float64Array(sO[ns]), sDead = fe.shellPointFailed = new Uint8Array(sO[ns]), sAl = fe.shellAlive = new Uint8Array(ns).fill(1), sRes = new Float64Array(7 * ns), sA0 = fe.shellArea = new Float64Array(ns); // sRes: transverse shear resultants (2) and hourglass forces (5)
    fe.shellPoints = sO;
    for (let e = 0; e < ns; e++) {
      const m = mats[sm[e]], t = sT[e], q = [0, 1, 2, 3].map((I) => 3 * sn[4 * e + I]), d = (a, b, k) => x[q[a] + k] - x[q[b] + k], side = (a, b) => Math.hypot(d(a, b, 0), d(a, b, 1), d(a, b, 2));
      const A = 0.5 * Math.hypot(d(2, 0, 1) * d(3, 1, 2) - d(2, 0, 2) * d(3, 1, 1), d(2, 0, 2) * d(3, 1, 0) - d(2, 0, 0) * d(3, 1, 2), d(2, 0, 0) * d(3, 1, 1) - d(2, 0, 1) * d(3, 1, 0));
      if (!(A > 0)) throw new Error('explicitfe: shell ' + e + ' has no area');
      const dt = (cfl * A) / Math.max(side(1, 0), side(2, 1), side(3, 2), side(0, 3)) / Math.sqrt(m.E / (m.rho * (1 - m.nu * m.nu))), mn = (m.rho * t * A) / 4, al = Math.max((t * t) / 12, A / 8);
      if (dt < dtE) dtE = dt; sA0[e] = A;
      for (let I = 0; I < 4; I++) { ms[sn[4 * e + I]] += mn; In[sn[4 * e + I]] += mn * al; fe.addedInertia += mn * (al - (t * t) / 12); }
      for (let p = sO[e]; p < sO[e + 1]; p++) { sTe[p] = m.T0; sYs[p] = m.A; }
    }
    // beams
    const bn = fe.beamN = Int32Array.from(bld.beam), sec = fe.sec = bld.sec, bE = new Float64Array(9 * nb), bL = fe.beamL = new Float64Array(nb), bd = new Float64Array(6 * nb), bF = fe.beamF = new Float64Array(6 * nb), bk = fe.beamKappa = new Float64Array(2 * nb), bep = fe.beamEps = new Float64Array(nb), bAl = fe.beamAlive = new Uint8Array(nb).fill(1), bEi = new Float64Array(nb), bT = fe.beamFailT = new Float64Array(nb).fill(NaN), bY = fe.beamYieldT = new Float64Array(nb).fill(NaN);
    for (let e = 0; e < nb; e++) {
      const s = sec[e], i1 = 3 * bn[2 * e], i2 = 3 * bn[2 * e + 1], dx = x[i2] - x[i1], dy = x[i2 + 1] - x[i1 + 1], dz = x[i2 + 2] - x[i1 + 2], L = Math.hypot(dx, dy, dz), e1 = [dx / L, dy / L, dz / L];
      let r = bld.ref[e], d = r[0] * e1[0] + r[1] * e1[1] + r[2] * e1[2]; if (Math.abs(d) > 0.999 * Math.hypot(r[0], r[1], r[2])) { r = Math.abs(e1[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]; d = r[0] * e1[0] + r[1] * e1[1] + r[2] * e1[2]; }
      let e2 = [r[0] - d * e1[0], r[1] - d * e1[1], r[2] - d * e1[2]]; const n2 = Math.hypot(e2[0], e2[1], e2[2]); e2 = e2.map((q) => q / n2);
      const e3 = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
      for (let k = 0; k < 3; k++) { bE[9 * e + k] = e1[k]; bE[9 * e + 3 + k] = e2[k]; bE[9 * e + 6 + k] = e3[k]; }
      bL[e] = L; s.G = s.G || s.E / 2.6; s.J = s.J || s.Iy + s.Iz;
      const m = s.rho * s.A * L, Imax = Math.max(s.Iy, s.Iz), r2 = Imax / s.A, kt = Math.max((2 * s.E * s.A) / L, (36 * s.E * Imax) / L ** 3), kr = Math.max((18 * s.E * Imax) / L, (2 * s.G * s.J) / L);
      // rotary inertia per end from the element's own mass; raised below where a larger step would otherwise be controlled by rotation
      const Ir = m * Math.max((L * L) / 24, 4.5 * r2, (0.5 * s.G * s.J) / (s.E * s.A));
      for (const nd of [bn[2 * e], bn[2 * e + 1]]) { ms[nd] += m / 2; In[nd] += Ir; kT[nd] += kt; kR[nd] += kr; }
    }
    for (const [nd, m, I] of bld.pm) { ms[nd] += m; if (hasRot && I) In[nd] += I; }
    for (const [nd, mask] of bld.fix) for (let k = 0; k < 6; k++) if (mask[k]) fixed[6 * nd + k] = 1;
    for (const [nd, q] of bld.load) for (let k = 0; k < 3; k++) fl[3 * nd + k] += q[k];
    for (const [nd, q] of bld.v0) { if (nd < 0) for (let i = 0; i < nn; i++) { v[3 * i] = q[0]; v[3 * i + 1] = q[1]; v[3 * i + 2] = q[2]; } else { v[3 * nd] = q[0]; v[3 * nd + 1] = q[1]; v[3 * nd + 2] = q[2]; } }
    for (let i = 0; i < nn; i++) { if (!(ms[i] > 0)) { if (fixed[6 * i] && fixed[6 * i + 1] && fixed[6 * i + 2]) ms[i] = 1; else throw new Error('explicitfe: node ' + i + ' has no mass'); } for (let k = 0; k < 3; k++) if (fixed[6 * i + k]) v[3 * i + k] = 0; }
    // springs and crush contacts
    const spr = fe.springs = bld.spr;
    for (const s of spr) {
      const i1 = 3 * s.n1, i2 = 3 * s.n2; s.r0 = [x[i2] - x[i1], x[i2 + 1] - x[i1 + 1], x[i2 + 2] - x[i1 + 2]]; s.L0 = Math.hypot(...s.r0);
      const free = (n) => !(fixed[6 * n] && fixed[6 * n + 1] && fixed[6 * n + 2]), mi = (free(s.n1) ? 1 / ms[s.n1] : 0) + (free(s.n2) ? 1 / ms[s.n2] : 0), km = Math.max(s.law.k, s.law.kb, s.law.kLat);
      if (mi > 0 && km > 0) { const om = Math.sqrt(km * mi), z = (s.law.c * mi) / (2 * om), dt = (cfl * 2 * (Math.sqrt(1 + z * z) - z)) / om; if (dt < dtE) dtE = dt; kT[s.n1] += 2 * km; kT[s.n2] += 2 * km; }
    }
    for (const pl of planes) for (const c of pl.pts) if (c.law) kT[c.node] += Math.max(c.law.k, c.law.kb || 0) * (c.off ? 2.25 : 1);
    // highest translational frequency from the assembled nodal masses; rotary inertia is raised where needed so that rotation does not control the step
    let w2 = 0; for (let i = 0; i < nn; i++) if (kT[i] > 0 && !(fixed[6 * i] && fixed[6 * i + 1] && fixed[6 * i + 2])) w2 = Math.max(w2, kT[i] / ms[i]);
    if (w2 > 0) { const dt = (cfl * 2) / Math.sqrt(w2); if (dt < dtE) dtE = dt; const wr = 4 / ((opt.dtMax ? Math.min(dtE, opt.dtMax) : dtE) / cfl) ** 2; for (let i = 0; i < In.length; i++) if (kR[i] > 0) { const need = kR[i] / wr; if (need > In[i]) { fe.addedInertia += need - In[i]; In[i] = need; } } }
    if (opt.dtMax && dtE > opt.dtMax) dtE = opt.dtMax;
    if (!Number.isFinite(dtE)) dtE = opt.dtMax || 1e-4;
    fe.dt = dtE; fe.totalMass = 0; for (let i = 0; i < nn; i++) fe.totalMass += ms[i];
    let dtPrev = dtE, dtHex = Infinity, dtShell = Infinity;
    const keRot = () => { let k = 0; for (let i = 0; i < In.length; i++) k += 0.5 * In[i] * (w[3 * i] ** 2 + w[3 * i + 1] ** 2 + w[3 * i + 2] ** 2); return k; };
    const keAll = () => { let k = 0; for (let i = 0; i < nn; i++) k += 0.5 * ms[i] * (v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2); for (const im of imps) if (Number.isFinite(im.mass)) k += 0.5 * im.mass * (im.v[0] ** 2 + im.v[1] ** 2 + im.v[2] ** 2); return k + keRot(); };
    E.initial = E.kinetic = E.peak = keAll();

    // ---- internal forces -------------------------------------------------------------------------
    function hexForces(dts) {
      let dtMin = Infinity, eInt = 0, ePl = 0;
      for (let e = 0; e < nh; e++) {
        if (!alive[e]) continue;
        const m = mats[hm[e]], e8 = 8 * e, s6 = 6 * e;
        for (let I = 0; I < 8; I++) { const k = 3 * hx[e8 + I]; X[I] = x[k]; X[8 + I] = x[k + 1]; X[16 + I] = x[k + 2]; U[I] = v[k]; U[8 + I] = v[k + 1]; U[16 + I] = v[k + 2]; }
        const V = hexGradient(X, Bm);
        if (!(V > 0.02 * hvol[e])) { alive[e] = 0; fe.erodedHex++; continue; } // inverted or collapsed: erode
        const iV = 1 / V; let l11 = 0, l12 = 0, l13 = 0, l21 = 0, l22 = 0, l23 = 0, l31 = 0, l32 = 0, l33 = 0, gx = 0, gy = 0, gz = 0, gxy = 0, gxz = 0, gyz = 0;
        for (let I = 0; I < 8; I++) {
          const bx = Bm[I], by = Bm[8 + I], bz = Bm[16 + I], ux = U[I], uy = U[8 + I], uz = U[16 + I];
          l11 += ux * bx; l12 += ux * by; l13 += ux * bz; l21 += uy * bx; l22 += uy * by; l23 += uy * bz; l31 += uz * bx; l32 += uz * by; l33 += uz * bz;
          gx += bx * bx; gy += by * by; gz += bz * bz; gxy += bx * by; gxz += bx * bz; gyz += by * bz;
        }
        const d11 = l11 * iV, d22 = l22 * iV, d33 = l33 * iV, d12 = 0.5 * (l12 + l21) * iV, d13 = 0.5 * (l13 + l31) * iV, d23 = 0.5 * (l23 + l32) * iV, w12 = 0.5 * (l12 - l21) * iV, w13 = 0.5 * (l13 - l31) * iV, w23 = 0.5 * (l23 - l32) * iV, tr = d11 + d22 + d33;
        let s11 = sig[s6], s22 = sig[s6 + 1], s33 = sig[s6 + 2], s12 = sig[s6 + 3], s13 = sig[s6 + 4], s23 = sig[s6 + 5];
        const o11 = s11, o22 = s22, o33 = s33, o12 = s12, o13 = s13, o23 = s23;
        // Jaumann rate: σ ← σ + Δt (W σ − σ W)
        if (dts > 0) {
          const r11 = 2 * (w12 * s12 + w13 * s13), r22 = 2 * (-w12 * s12 + w23 * s23), r33 = -2 * (w13 * s13 + w23 * s23);
          const r12 = w12 * (s22 - s11) + w13 * s23 + w23 * s13, r13 = w13 * (s33 - s11) + w12 * s23 - w23 * s12, r23 = w23 * (s33 - s22) - w12 * s13 - w13 * s12;
          s11 += dts * r11; s22 += dts * r22; s33 += dts * r33; s12 += dts * r12; s13 += dts * r13; s23 += dts * r23;
        }
        let p = (s11 + s22 + s33) / 3 + m.K * tr * dts; const G2 = 2 * m.G * dts, t3 = tr / 3;
        let q11 = s11 - (s11 + s22 + s33) / 3 + G2 * (d11 - t3), q22 = s22 - (s11 + s22 + s33) / 3 + G2 * (d22 - t3), q33, q12 = s12 + G2 * d12, q13 = s13 + G2 * d13, q23 = s23 + G2 * d23;
        q33 = -q11 - q22;
        const qe = Math.sqrt(1.5 * (q11 * q11 + q22 * q22 + q33 * q33 + 2 * (q12 * q12 + q13 * q13 + q23 * q23)));
        if (qe > sys[e]) { // strain rate can only raise the flow stress above the cached rate-free value
          const rate = Math.sqrt((2 / 3) * ((d11 - t3) ** 2 + (d22 - t3) ** 2 + (d33 - t3) ** 2 + 2 * (d12 * d12 + d13 * d13 + d23 * d23))), e0 = ep[e], T = Te[e], sy0 = flowStress(m, e0, rate, T);
          if (qe > sy0) {
            const G3 = 3 * m.G; let lo = 0, hi = (qe - sy0) / G3, glo = qe - sy0, ghi = qe - G3 * hi - flowStress(m, e0 + hi, rate, T), de = hi, sy = sy0;
            if (ghi < 0) for (let k = 0; k < 20; k++) {
              de = k % 3 === 2 ? 0.5 * (lo + hi) : lo + ((hi - lo) * glo) / (glo - ghi); sy = flowStress(m, e0 + de, rate, T);
              const g = qe - G3 * de - sy; if (g > 0) { lo = de; glo = g; } else { hi = de; ghi = g; }
              if (Math.abs(g) < 1e-6 * qe) break;
            } else sy = qe - G3 * hi;
            const sc = Math.max(0, qe - G3 * de) / qe; q11 *= sc; q22 *= sc; q33 *= sc; q12 *= sc; q13 *= sc; q23 *= sc;
            ep[e] = e0 + de; ePl += sc * qe * de * V; if (m.thermal) Te[e] = T + (m.beta * sc * qe * de) / (m.rho * m.cp);
            sys[e] = flowStress(m, ep[e], 0, Te[e]);
            if (m.efail > 0 && ep[e] >= m.efail) { alive[e] = 0; fe.erodedHex++; }
          }
        }
        s11 = q11 + p; s22 = q22 + p; s33 = q33 + p; s12 = q12; s13 = q13; s23 = q23;
        sig[s6] = s11; sig[s6 + 1] = s22; sig[s6 + 2] = s33; sig[s6 + 3] = s12; sig[s6 + 4] = s13; sig[s6 + 5] = s23;
        // bulk viscosity and stable time step
        const Gmax = Math.max(gx + Math.abs(gxy) + Math.abs(gxz), gy + Math.abs(gxy) + Math.abs(gyz), gz + Math.abs(gxz) + Math.abs(gyz)), l = V / Math.sqrt(2 * Gmax), rho = (hrho[e] * hvol[e]) * iV, sr = Math.sqrt(m.rho / hrho[e]), c = m.cd * sr, cs = m.cs * sr;
        let qv = 0, Qd = 0; if (tr < 0) { Qd = bq1 * c + bq2 * l * -tr; qv = rho * l * -tr * Qd; }
        const dt = (cfl * l) / (Qd + Math.sqrt(Qd * Qd + cs * cs)); if (dt < dtMin) dtMin = dt;
        eInt += dts * V * (0.5 * ((o11 + s11) * d11 + (o22 + s22) * d22 + (o33 + s33) * d33) + (o12 + s12) * d12 + (o13 + s13) * d13 + (o23 + s23) * d23);
        if (!alive[e]) continue;
        for (let I = 0; I < 8; I++) { const k = 3 * hx[e8 + I], bx = Bm[I], by = Bm[8 + I], bz = Bm[16 + I]; f[k] -= s11 * bx + s12 * by + s13 * bz; f[k + 1] -= s12 * bx + s22 * by + s23 * bz; f[k + 2] -= s13 * bx + s23 * by + s33 * bz; }
        if (qv > 0) for (let I = 0; I < 8; I++) { const k = 3 * hx[e8 + I]; fq[k] += qv * Bm[I]; fq[k + 1] += qv * Bm[8 + I]; fq[k + 2] += qv * Bm[16 + I]; }
        // hourglass control: stiffness form plus a small viscous part
        // The stiffness scales with the harmonic mean of the eigenvalues of G = B·Bᵀ (the smallest sectional dimension), which
        // gives a thin element roughly its plate-bending stiffness; the hourglass force is capped at the plastic-bending level.
        const I2 = gx * gy + gx * gz + gy * gz - gxy * gxy - gxz * gxz - gyz * gyz, dG = gx * (gy * gz - gyz * gyz) - gxy * (gxy * gz - gyz * gxz) + gxz * (gxy * gyz - gy * gxz), lh = I2 > 0 ? Math.max(0, (3 * dG) / I2) : 0;
        const ks = 0.5 * kHG * (m.K + (4 / 3) * m.G) * lh * iV, cvs = (qHG * rho * c * Math.cbrt(V * V)) / 4, qp = (sys[e] * Math.sqrt((4 * lh) / 3)) / 8;
        for (let al = 0; al < 4; al++) {
          let hx_ = 0, hy = 0, hz = 0; const g8 = 8 * al;
          for (let I = 0; I < 8; I++) { const g = GAM[g8 + I]; hx_ += X[I] * g; hy += X[8 + I] * g; hz += X[16 + I] * g; }
          let qx = 0, qy = 0, qz = 0;
          for (let I = 0; I < 8; I++) { const g = GAM[g8 + I] - (Bm[I] * hx_ + Bm[8 + I] * hy + Bm[16 + I] * hz) * iV; gm[I] = g; qx += U[I] * g; qy += U[8 + I] * g; qz += U[16 + I] * g; }
          const h = 12 * e + 3 * al; let nx = hq[h] + dts * ks * qx, ny = hq[h + 1] + dts * ks * qy, nz = hq[h + 2] + dts * ks * qz; const qn = Math.sqrt(nx * nx + ny * ny + nz * nz);
          if (qn > qp) { const r = qp / qn; nx *= r; ny *= r; nz *= r; }
          const fx = nx + cvs * qx, fy = ny + cvs * qy, fz = nz + cvs * qz;
          hq[h] = nx; hq[h + 1] = ny; hq[h + 2] = nz;
          for (let I = 0; I < 8; I++) { const k = 3 * hx[e8 + I], g = gm[I]; fh[k] -= fx * g; fh[k + 1] -= fy * g; fh[k + 2] -= fz * g; }
        }
      }
      dtHex = dtMin; E.internal += eInt; E.plastic += ePl;
    }
    const sX = new Float64Array(4), sY = new Float64Array(4), sB1 = new Float64Array(4), sB2 = new Float64Array(4), sGm = new Float64Array(4), sV = new Float64Array(20), sI = new Int32Array(4), SH = [1, -1, 1, -1];
    /**
     * Shell forces. Velocity strains in the element's co-rotational frame (e3 normal to the two diagonals) at the single
     * in-plane point: membrane d = B·v, curvature rates from the nodal rotation rates (κx = B1·θy, κy = −B2·θx,
     * κxy = B2·θy − B1·θx) and transverse shear γxz = B1·vz + mean θy, γyz = B2·vz − mean θx. Plane-stress J2 plasticity is
     * integrated at each thickness point; the transverse shear resultant is elastic (factor 5/6) up to the shear yield of the
     * section. The five zero-energy modes of the one-point element carry a small elastic stabilisation.
     */
    function shellForces(dts) {
      let dtMin = Infinity, eInt = 0, ePl = 0, eHg = 0;
      for (let e = 0; e < ns; e++) {
        if (!sAl[e]) continue;
        const m = mats[sm[e]], t = sT[e], nq = sQ[e], e4 = 4 * e; for (let I = 0; I < 4; I++) sI[I] = 3 * sn[e4 + I];
        const n1 = sI[0], n2 = sI[1], n3 = sI[2], n4 = sI[3], ax = x[n3] - x[n1], ay = x[n3 + 1] - x[n1 + 1], az = x[n3 + 2] - x[n1 + 2], bx = x[n4] - x[n2], by = x[n4 + 1] - x[n2 + 1], bz = x[n4 + 2] - x[n2 + 2];
        let e3x = ay * bz - az * by, e3y = az * bx - ax * bz, e3z = ax * by - ay * bx; const n3l = Math.hypot(e3x, e3y, e3z), A = 0.5 * n3l;
        if (!(A > 0.02 * sA0[e])) { sAl[e] = 0; fe.erodedShell++; continue; } // collapsed: erode
        e3x /= n3l; e3y /= n3l; e3z /= n3l;
        let e1x = x[n2] + x[n3] - x[n1] - x[n4], e1y = x[n2 + 1] + x[n3 + 1] - x[n1 + 1] - x[n4 + 1], e1z = x[n2 + 2] + x[n3 + 2] - x[n1 + 2] - x[n4 + 2]; const d3 = e1x * e3x + e1y * e3y + e1z * e3z; e1x -= d3 * e3x; e1y -= d3 * e3y; e1z -= d3 * e3z;
        const n1l = Math.hypot(e1x, e1y, e1z); e1x /= n1l; e1y /= n1l; e1z /= n1l; const e2x = e3y * e1z - e3z * e1y, e2y = e3z * e1x - e3x * e1z, e2z = e3x * e1y - e3y * e1x;
        for (let I = 0; I < 4; I++) {
          const k = sI[I], px = x[k] - x[n1], py = x[k + 1] - x[n1 + 1], pz = x[k + 2] - x[n1 + 2], o = 5 * I; sX[I] = px * e1x + py * e1y + pz * e1z; sY[I] = px * e2x + py * e2y + pz * e2z;
          sV[o] = v[k] * e1x + v[k + 1] * e1y + v[k + 2] * e1z; sV[o + 1] = v[k] * e2x + v[k + 1] * e2y + v[k + 2] * e2z; sV[o + 2] = v[k] * e3x + v[k + 1] * e3y + v[k + 2] * e3z;
          sV[o + 3] = w[k] * e1x + w[k + 1] * e1y + w[k + 2] * e1z; sV[o + 4] = w[k] * e2x + w[k + 1] * e2y + w[k + 2] * e2z;
        }
        const i2A = 1 / (2 * A); sB1[0] = (sY[1] - sY[3]) * i2A; sB1[1] = (sY[2] - sY[0]) * i2A; sB1[2] = -sB1[0]; sB1[3] = -sB1[1]; sB2[0] = (sX[3] - sX[1]) * i2A; sB2[1] = (sX[0] - sX[2]) * i2A; sB2[2] = -sB2[0]; sB2[3] = -sB2[1];
        const hx_ = sX[0] - sX[1] + sX[2] - sX[3], hy = sY[0] - sY[1] + sY[2] - sY[3];
        let dxm = 0, dym = 0, gm = 0, kx = 0, ky = 0, kxy = 0, gxz = 0, gyz = 0, BB = 0, qmx = 0, qmy = 0, qw = 0, qbx = 0, qby = 0;
        for (let I = 0; I < 4; I++) {
          const o = 5 * I, b1 = sB1[I], b2 = sB2[I], vx = sV[o], vy = sV[o + 1], vz = sV[o + 2], tx = sV[o + 3], ty = sV[o + 4], g = SH[I] - hx_ * b1 - hy * b2; sGm[I] = g;
          dxm += b1 * vx; dym += b2 * vy; gm += b2 * vx + b1 * vy; kx += b1 * ty; ky -= b2 * tx; kxy += b2 * ty - b1 * tx; gxz += b1 * vz + 0.25 * ty; gyz += b2 * vz - 0.25 * tx;
          BB += b1 * b1 + b2 * b2; qmx += g * vx; qmy += g * vy; qw += g * vz; qbx += g * tx; qby += g * ty;
        }
        // through-thickness integration of the plane-stress response
        const nu = m.nu, C1 = m.E / (1 - nu * nu), G = m.G, kp = m.E / (3 * (1 - nu)), kd = 2 * G, ht = 0.5 * t, gp = GQ[nq - 1], p0 = sO[e];
        let Nx = 0, Ny = 0, Nxy = 0, Mx = 0, My = 0, Mxy = 0, live = 0, ysum = 0;
        for (let p = 0; p < nq; p++) {
          const ip = p0 + p; if (sDead[ip]) continue;
          const z = ht * gp[0][p], wt = ht * gp[1][p], dx = dxm + z * kx, dy = dym + z * ky, gxy = gm + z * kxy, s3 = 3 * ip, ox = sSig[s3], oy = sSig[s3 + 1], oxy = sSig[s3 + 2];
          let sx = ox + dts * C1 * (dx + nu * dy), sy = oy + dts * C1 * (dy + nu * dx), txy = oxy + dts * G * gxy;
          const pm = 0.5 * (sx + sy), dm = 0.5 * (sx - sy), q2 = pm * pm + 3 * (dm * dm + txy * txy);
          if (q2 > sYs[ip] * sYs[ip]) { // strain rate can only raise the flow stress above the cached rate-free value
            const qe = Math.sqrt(q2), rate = Math.sqrt((2 / 3) * (dx * dx + dy * dy + (dx + dy) ** 2 + 0.5 * gxy * gxy)), e0 = sEp[ip], T = sTe[ip], sy0 = flowStress(m, e0, rate, T);
            if (qe > sy0) {
              // plane-stress radial return: the mean part p = (σx + σy)/2 and the deviatoric parts shrink by 1/(1 + λ·E/(3(1 − ν))) and
              // 1/(1 + 2Gλ); the plastic multiplier λ follows from q(λ) = σy(εp + ⅔·λ·q(λ)), bracketed between 0 and the perfectly plastic bound
              let lo = 0, glo = qe - sy0, hi = (qe / sy0 - 1) / Math.min(kp, kd), ghi = NaN, lam = hi, fp = 1, fd = 1, qn = qe, de = 0;
              for (let k = 0; k < 40; k++) {
                fp = 1 / (1 + lam * kp); fd = 1 / (1 + lam * kd); qn = Math.sqrt(pm * pm * fp * fp + 3 * (dm * dm + txy * txy) * fd * fd); de = (2 / 3) * lam * qn;
                const g = qn - flowStress(m, e0 + de, rate, T); if (Math.abs(g) < 1e-10 * qe || (k === 0 && g > 0)) break;
                if (g > 0) { lo = lam; glo = g; } else { hi = lam; ghi = g; }
                lam = k % 3 === 2 ? 0.5 * (lo + hi) : lo + ((hi - lo) * glo) / (glo - ghi);
              }
              sx = pm * fp + dm * fd; sy = pm * fp - dm * fd; txy *= fd; sEp[ip] = e0 + de; ePl += A * wt * qn * de;
              if (m.thermal) sTe[ip] = T + (m.beta * qn * de) / (m.rho * m.cp);
              sYs[ip] = flowStress(m, sEp[ip], 0, sTe[ip]);
              if (m.efail > 0 && sEp[ip] >= m.efail) { sDead[ip] = 1; sx = sy = txy = 0; }
            }
          }
          sSig[s3] = sx; sSig[s3 + 1] = sy; sSig[s3 + 2] = txy;
          eInt += dts * A * wt * (0.5 * ((ox + sx) * dx + (oy + sy) * dy) + 0.5 * (oxy + txy) * gxy);
          if (sDead[ip]) continue;
          live++; ysum += sYs[ip]; Nx += wt * sx; Ny += wt * sy; Nxy += wt * txy; Mx += wt * z * sx; My += wt * z * sy; Mxy += wt * z * txy;
        }
        if (!live) { sAl[e] = 0; fe.erodedShell++; continue; }
        // transverse shear resultants and stabilisation forces, reduced with the failed share of the thickness
        const r7 = 7 * e, fa = live / nq, kS = (5 / 6) * G * t * fa, oQx = sRes[r7], oQy = sRes[r7 + 1], qc = (fa * t * ysum) / (live * Math.sqrt(3));
        let Qx = oQx + dts * kS * gxz, Qy = oQy + dts * kS * gyz; const qn = Math.hypot(Qx, Qy); if (qn > qc) { Qx *= qc / qn; Qy *= qc / qn; }
        sRes[r7] = Qx; sRes[r7 + 1] = Qy; eInt += dts * A * 0.5 * ((oQx + Qx) * gxz + (oQy + Qy) * gyz);
        const cm = (rSH * fa * m.E * t * A * BB) / 8, cw = (rSH * fa * (5 / 6) * G * t ** 3 * A * BB) / 12, cb = (rSH * fa * m.E * t ** 3 * A * BB) / 192;
        const Hmx = sRes[r7 + 2] + dts * cm * qmx, Hmy = sRes[r7 + 3] + dts * cm * qmy, Hw = sRes[r7 + 4] + dts * cw * qw, Hbx = sRes[r7 + 5] + dts * cb * qbx, Hby = sRes[r7 + 6] + dts * cb * qby;
        eHg += dts * 0.5 * ((sRes[r7 + 2] + Hmx) * qmx + (sRes[r7 + 3] + Hmy) * qmy + (sRes[r7 + 4] + Hw) * qw + (sRes[r7 + 5] + Hbx) * qbx + (sRes[r7 + 6] + Hby) * qby);
        sRes[r7 + 2] = Hmx; sRes[r7 + 3] = Hmy; sRes[r7 + 4] = Hw; sRes[r7 + 5] = Hbx; sRes[r7 + 6] = Hby;
        const s12 = Math.hypot(sX[1], sY[1]), s23 = Math.hypot(sX[2] - sX[1], sY[2] - sY[1]), s34 = Math.hypot(sX[3] - sX[2], sY[3] - sY[2]), s41 = Math.hypot(sX[3], sY[3]), dt = (cfl * A) / Math.max(s12, s23, s34, s41) / Math.sqrt(C1 / m.rho); if (dt < dtMin) dtMin = dt;
        for (let I = 0; I < 4; I++) {
          const k = sI[I], b1 = sB1[I], b2 = sB2[I], g = sGm[I];
          const fx = A * (b1 * Nx + b2 * Nxy) + g * Hmx, fy = A * (b2 * Ny + b1 * Nxy) + g * Hmy, fz = A * (b1 * Qx + b2 * Qy) + g * Hw, mx = A * (-b2 * My - b1 * Mxy - 0.25 * Qy) + g * Hbx, my = A * (b1 * Mx + b2 * Mxy + 0.25 * Qx) + g * Hby;
          f[k] -= fx * e1x + fy * e2x + fz * e3x; f[k + 1] -= fx * e1y + fy * e2y + fz * e3y; f[k + 2] -= fx * e1z + fy * e2z + fz * e3z;
          mo[k] -= mx * e1x + my * e2x; mo[k + 1] -= mx * e1y + my * e2y; mo[k + 2] -= mx * e1z + my * e2z;
        }
      }
      dtShell = dtMin; E.internal += eInt; E.plastic += ePl; E.hourglass += eHg;
    }
    const dd = new Float64Array(6), Fo = new Float64Array(6), tA = new Float64Array(9), tB = new Float64Array(9);
    const hp = new Float64Array(4); // plastic rotation increments of the two hinges: [end 1 about y, about z, end 2 about y, about z]
    /**
     * Radial return of one beam end onto its moment interaction surface |(My/Mpy, Mz/Mpz)| = allow, from the moments (ay, az)
     * the end would carry without plastic rotation of its own. Stores that rotation in hp[q], hp[q + 1] and returns its
     * work-equivalent magnitude (the plastic rotation itself in uniaxial bending).
     */
    function hinge(ay, az, s, ky, kz, allow, q) {
      const ry = ay / s.Mpy, rz = az / s.Mpz, mm = Math.hypot(ry, rz);
      if (!(mm > allow)) { hp[q] = hp[q + 1] = 0; return 0; }
      const u = 1 - allow / mm; hp[q] = ry ? (u * ay) / (4 * ky) : 0; hp[q + 1] = rz ? (u * az) / (4 * kz) : 0;
      return (hp[q] * ry + hp[q + 1] * rz) / mm;
    }
    function beamForces(dts) {
      for (let e = 0; e < nb; e++) {
        if (!bAl[e]) continue;
        const s = sec[e], n1 = bn[2 * e], n2 = bn[2 * e + 1], i1 = 3 * n1, i2 = 3 * n2, L0 = bL[e], b9 = 9 * e;
        const dx = x[i2] - x[i1], dy = x[i2 + 1] - x[i1 + 1], dz = x[i2 + 2] - x[i1 + 2], L = Math.hypot(dx, dy, dz), e1x = dx / L, e1y = dy / L, e1z = dz / L;
        // nodal triads T = R·E0 (columns t1, t2, t3)
        const r1 = 9 * n1, r2 = 9 * n2;
        for (let c = 0; c < 3; c++) { const ax = bE[b9 + 3 * c], ay = bE[b9 + 3 * c + 1], az = bE[b9 + 3 * c + 2], q = 3 * c;
          tA[q] = R[r1] * ax + R[r1 + 1] * ay + R[r1 + 2] * az; tA[q + 1] = R[r1 + 3] * ax + R[r1 + 4] * ay + R[r1 + 5] * az; tA[q + 2] = R[r1 + 6] * ax + R[r1 + 7] * ay + R[r1 + 8] * az;
          tB[q] = R[r2] * ax + R[r2 + 1] * ay + R[r2 + 2] * az; tB[q + 1] = R[r2 + 3] * ax + R[r2 + 4] * ay + R[r2 + 5] * az; tB[q + 2] = R[r2 + 6] * ax + R[r2 + 7] * ay + R[r2 + 8] * az; }
        let mx = tA[3] + tB[3], my = tA[4] + tB[4], mz = tA[5] + tB[5]; const dm = mx * e1x + my * e1y + mz * e1z; mx -= dm * e1x; my -= dm * e1y; mz -= dm * e1z;
        const nm = Math.hypot(mx, my, mz) || 1, e2x = mx / nm, e2y = my / nm, e2z = mz / nm, e3x = e1y * e2z - e1z * e2y, e3y = e1z * e2x - e1x * e2z, e3z = e1x * e2y - e1y * e2x;
        // deformation increments in rate form, work-conjugate to the nodal forces and evaluated in the mid-step
        // element frame: axial stretch, twist and the two bending rotations of each end relative to the chord
        const d6 = 6 * e, vx = v[i2] - v[i1], vy = v[i2 + 1] - v[i1 + 1], vz = v[i2 + 2] - v[i1 + 2], hx_ = dx - 0.5 * dts * vx, hy = dy - 0.5 * dts * vy, hz = dz - 0.5 * dts * vz, Lm = Math.hypot(hx_, hy, hz), f1x = hx_ / Lm, f1y = hy / Lm, f1z = hz / Lm;
        let f2x = e2x, f2y = e2y, f2z = e2z; const df = f2x * f1x + f2y * f1y + f2z * f1z; f2x -= df * f1x; f2y -= df * f1y; f2z -= df * f1z; const nf = Math.hypot(f2x, f2y, f2z) || 1; f2x /= nf; f2y /= nf; f2z /= nf;
        const f3x = f1y * f2z - f1z * f2y, f3y = f1z * f2x - f1x * f2z, f3z = f1x * f2y - f1y * f2x, sy = ((f3x * vx + f3y * vy + f3z * vz) / Lm), sz = ((f2x * vx + f2y * vy + f2z * vz) / Lm);
        dd[0] = L - L0 - bd[d6]; bd[d6] = L - L0;
        dd[1] = dts * (f1x * (w[i2] - w[i1]) + f1y * (w[i2 + 1] - w[i1 + 1]) + f1z * (w[i2 + 2] - w[i1 + 2]));
        dd[2] = dts * (f2x * w[i1] + f2y * w[i1 + 1] + f2z * w[i1 + 2] + sy); dd[3] = dts * (f3x * w[i1] + f3y * w[i1 + 1] + f3z * w[i1 + 2] - sz);
        dd[4] = dts * (f2x * w[i2] + f2y * w[i2 + 1] + f2z * w[i2 + 2] + sy); dd[5] = dts * (f3x * w[i2] + f3y * w[i2 + 1] + f3z * w[i2 + 2] - sz);
        for (let k = 0; k < 6; k++) Fo[k] = bF[d6 + k];
        const EA = (s.E * s.A) / L0, GJ = (s.G * s.J) / L0, ky = (s.E * s.Iy) / L0, kz = (s.E * s.Iz) / L0;
        let Nf = Fo[0] + EA * dd[0], Tq = Fo[1] + GJ * dd[1], M1y = Fo[2] + ky * (4 * dd[2] + 2 * dd[4]), M2y = Fo[4] + ky * (2 * dd[2] + 4 * dd[4]), M1z = Fo[3] + kz * (4 * dd[3] + 2 * dd[5]), M2z = Fo[5] + kz * (2 * dd[3] + 4 * dd[5]);
        // plastic resultants: axial squash load and a hinge at each end with axial-force interaction. A plastic rotation at one
        // end relieves the other end through the carry-over term of the elastic element (M1 = 4k·θ1 + 2k·θ2), so the two radial
        // returns are coupled; they are solved together by block Gauss–Seidel (contraction 1/4 per sweep).
        if (s.Np < Infinity || s.Mpy < Infinity || s.Mpz < Infinity) {
          const NpH = s.Np;
          if (Math.abs(Nf) > NpH) { if (!(bY[e] >= 0)) bY[e] = fe.t; bep[e] += (Math.abs(Nf) - NpH) / (s.E * s.A); Nf = Math.sign(Nf) * NpH; }
          const ax = Math.max(0.02, 1 - (Nf / NpH) ** 2), k1 = bk[2 * e], k2 = bk[2 * e + 1], al1 = ax * (k1 >= s.thetaSoft ? s.residual : Math.min(s.cap, 1 + s.hard * k1)), al2 = ax * (k2 >= s.thetaSoft ? s.residual : Math.min(s.cap, 1 + s.hard * k2));
          let g1 = 0, g2 = 0; hp[0] = hp[1] = hp[2] = hp[3] = 0;
          for (let it = 0; it < 60; it++) {
            const o0 = hp[0], o1 = hp[1], o2 = hp[2], o3 = hp[3];
            g1 = hinge(M1y - 2 * ky * o2, M1z - 2 * kz * o3, s, ky, kz, al1, 0); g2 = hinge(M2y - 2 * ky * hp[0], M2z - 2 * kz * hp[1], s, ky, kz, al2, 2);
            if (Math.abs(hp[0] - o0) + Math.abs(hp[1] - o1) + Math.abs(hp[2] - o2) + Math.abs(hp[3] - o3) <= 1e-14 * (Math.abs(hp[0]) + Math.abs(hp[1]) + Math.abs(hp[2]) + Math.abs(hp[3]))) break;
          }
          if (g1 > 0 || g2 > 0) {
            if (!(bY[e] >= 0)) bY[e] = fe.t; bk[2 * e] += g1; bk[2 * e + 1] += g2;
            M1y -= ky * (4 * hp[0] + 2 * hp[2]); M2y -= ky * (2 * hp[0] + 4 * hp[2]); M1z -= kz * (4 * hp[1] + 2 * hp[3]); M2z -= kz * (2 * hp[1] + 4 * hp[3]);
          }
        }
        bF[d6] = Nf; bF[d6 + 1] = Tq; bF[d6 + 2] = M1y; bF[d6 + 3] = M1z; bF[d6 + 4] = M2y; bF[d6 + 5] = M2z;
        const dE = 0.5 * ((Fo[0] + Nf) * dd[0] + (Fo[1] + Tq) * dd[1] + (Fo[2] + M1y) * dd[2] + (Fo[3] + M1z) * dd[3] + (Fo[4] + M2y) * dd[4] + (Fo[5] + M2z) * dd[5]);
        bEi[e] += dE; E.internal += dE;
        const el = (Nf * Nf) / (2 * EA) + (Tq * Tq) / (2 * GJ) + (M1y * M1y - M1y * M2y + M2y * M2y) / (6 * ky) + (M1z * M1z - M1z * M2z + M2z * M2z) / (6 * kz);
        s.wp = Math.max(s.wp || 0, bEi[e] - el);
        if (bk[2 * e] >= s.thetaFail || bk[2 * e + 1] >= s.thetaFail || bep[e] >= s.epsFail) { bAl[e] = 0; fe.erodedBeam++; bT[e] = fe.t; s.wp = bEi[e]; continue; }
        const qy = (M1y + M2y) / L, qz = (M1z + M2z) / L, fx = Nf * e1x + qy * e3x - qz * e2x, fy = Nf * e1y + qy * e3y - qz * e2y, fz = Nf * e1z + qy * e3z - qz * e2z;
        f[i2] -= fx; f[i2 + 1] -= fy; f[i2 + 2] -= fz; f[i1] += fx; f[i1 + 1] += fy; f[i1 + 2] += fz;
        mo[i1] -= -Tq * e1x + M1y * e2x + M1z * e3x; mo[i1 + 1] -= -Tq * e1y + M1y * e2y + M1z * e3y; mo[i1 + 2] -= -Tq * e1z + M1y * e2z + M1z * e3z;
        mo[i2] -= Tq * e1x + M2y * e2x + M2z * e3x; mo[i2 + 1] -= Tq * e1y + M2y * e2y + M2z * e3y; mo[i2 + 2] -= Tq * e1z + M2y * e2z + M2z * e3z;
      }
    }
    function springForces() {
      for (const s of spr) {
        const i1 = 3 * s.n1, i2 = 3 * s.n2, lw = s.law, rx = x[i2] - x[i1], ry = x[i2 + 1] - x[i1 + 1], rz = x[i2 + 2] - x[i1 + 2];
        let ux, uy, uz, d;
        if (s.dir) { ux = s.dir[0]; uy = s.dir[1]; uz = s.dir[2]; d = -((rx - s.r0[0]) * ux + (ry - s.r0[1]) * uy + (rz - s.r0[2]) * uz); }
        else { const L = Math.hypot(rx, ry, rz) || 1e-300; ux = rx / L; uy = ry / L; uz = rz / L; d = s.L0 - L; }
        const rate = -((v[i2] - v[i1]) * ux + (v[i2 + 1] - v[i1 + 1]) * uy + (v[i2 + 2] - v[i1 + 2]) * uz);
        let F = lw.k * (d - s.dp);
        if (F > lw.Fyc) { const ndp = d - lw.Fyc / lw.k; s.wp += lw.Fyc * (ndp - s.dp); E.crush += lw.Fyc * (ndp - s.dp); s.dp = ndp; F = lw.Fyc; }
        else if (F < -lw.Fyt) { const ndp = d + lw.Fyt / lw.k; s.wp += lw.Fyt * (s.dp - ndp); E.crush += lw.Fyt * (s.dp - ndp); s.dp = ndp; F = -lw.Fyt; }
        if (d > lw.dmax) { F += lw.kb * (d - lw.dmax); s.bottomed = true; }
        F += lw.c * rate; if (d > s.dMaxSeen) s.dMaxSeen = d; if (s.dp > s.dpMax) s.dpMax = s.dp;
        E.internal += 0.5 * (F + s.F) * (d - s.d); s.d = d; s.F = F;
        f[i2] += F * ux; f[i2 + 1] += F * uy; f[i2 + 2] += F * uz; f[i1] -= F * ux; f[i1 + 1] -= F * uy; f[i1 + 2] -= F * uz;
        if (s.dir && lw.kLat) {
          const qx = rx - s.r0[0], qy = ry - s.r0[1], qz = rz - s.r0[2], al = qx * ux + qy * uy + qz * uz, px = qx - al * ux, py = qy - al * uy, pz = qz - al * uz, en = 0.5 * lw.kLat * (px * px + py * py + pz * pz);
          E.internal += en - s.e; s.e = en;
          f[i2] -= lw.kLat * px; f[i2 + 1] -= lw.kLat * py; f[i2 + 2] -= lw.kLat * pz; f[i1] += lw.kLat * px; f[i1 + 1] += lw.kLat * py; f[i1 + 2] += lw.kLat * pz;
        }
      }
    }
    function contactForces(dt) {
      for (const pl of planes) {
        const n = pl.n, p0 = pl.p0; pl.force = 0;
        for (const c of pl.pts) {
          c.F = 0; if (c.failed) continue;
          const i = 3 * c.node; let ox = 0, oy = 0, oz = 0, vx = v[i], vy = v[i + 1], vz = v[i + 2];
          if (c.off) { const r = 9 * c.node, o = c.off; ox = R[r] * o[0] + R[r + 1] * o[1] + R[r + 2] * o[2]; oy = R[r + 3] * o[0] + R[r + 4] * o[1] + R[r + 5] * o[2]; oz = R[r + 6] * o[0] + R[r + 7] * o[1] + R[r + 8] * o[2]; vx += w[i + 1] * oz - w[i + 2] * oy; vy += w[i + 2] * ox - w[i] * oz; vz += w[i] * oy - w[i + 1] * ox; }
          const d = -((x[i] + ox - p0[0]) * n[0] + (x[i + 1] + oy - p0[1]) * n[1] + (x[i + 2] + oz - p0[2]) * n[2]);
          if (!(d > 0)) continue;
          const vn = vx * n[0] + vy * n[1] + vz * n[2]; let F;
          if (pl.fluid) F = c.area * ((vn < 0 ? 0.5 * pl.fluid.rho * pl.fluid.Cd * vn * vn : 0) + pl.fluid.rho * 9.80665 * d);
          else if (c.law) {
            const lw = c.law; F = lw.k * (d - c.dp);
            if (F > lw.Fy) { const ndp = d - lw.Fy / lw.k; c.wp += lw.Fy * (ndp - c.dp); E.crush += lw.Fy * (ndp - c.dp); c.dp = ndp; F = lw.Fy; }
            if (F < 0) F = 0;
            if (d > lw.dmax) { if (lw.failAtMax) { c.failed = true; c.tFail = fe.t; c.dMax = d; continue; } F += (lw.kb || 0) * (d - lw.dmax); }
          } else F = ((pf * ms[c.node]) / (dt * dt)) * d;
          if (d > c.dMax) c.dMax = d;
          let fx = F * n[0], fy = F * n[1], fz = F * n[2];
          const mu = c.mu ?? pl.mu;
          if (mu > 0) {
            const tx = vx - vn * n[0], ty = vy - vn * n[1], tz = vz - vn * n[2], vt = Math.hypot(tx, ty, tz);
            if (vt > 1e-12) { const ft = Math.min(mu * F, (ms[c.node] * vt) / dt); fx -= (ft * tx) / vt; fy -= (ft * ty) / vt; fz -= (ft * tz) / vt; E.friction += ft * vt * dt; }
          }
          c.F = F; pl.force += F; fc[i] += fx; fc[i + 1] += fy; fc[i + 2] += fz;
          if (c.off) { mc[i] += oy * fz - oz * fy; mc[i + 1] += oz * fx - ox * fz; mc[i + 2] += ox * fy - oy * fx; }
        }
      }
      for (const im of imps) {
        im.f[0] = im.f[1] = im.f[2] = 0; im.force = 0; const ax = im.axis, c = im.c, fin = Number.isFinite(im.mass);
        for (let q = 0; q < im.nodes.length; q++) {
          const nd = im.nodes[q], i = 3 * nd; let rx = x[i] - c[0], ry = x[i + 1] - c[1], rz = x[i + 2] - c[2];
          if (ax) { const al = rx * ax[0] + ry * ax[1] + rz * ax[2]; rx -= al * ax[0]; ry -= al * ax[1]; rz -= al * ax[2]; }
          const d2 = rx * rx + ry * ry + rz * rz; if (d2 >= im.R * im.R) continue;
          const dist = Math.sqrt(d2) || 1e-300, pen = im.R - dist, me = fin ? (ms[nd] * im.mass) / (ms[nd] + im.mass) : ms[nd], F = Math.min(((pf * me) / (dt * dt)) * pen, im.fcaps ? im.fcaps[q] : im.fcap), fx = (F * rx) / dist, fy = (F * ry) / dist, fz = (F * rz) / dist;
          fc[i] += fx; fc[i + 1] += fy; fc[i + 2] += fz; im.f[0] -= fx; im.f[1] -= fy; im.f[2] -= fz; im.force += F;
        }
      }
    }

    /** Advance one time step. */
    fe.step = () => {
      const dts = fe.nStep ? dtPrev : 0; // the configuration x_n was reached with dtPrev and the half-step velocities v
      f.fill(0); fc.fill(0); if (nh) { fh.fill(0); fq.fill(0); } if (hasRot) { mo.fill(0); mc.fill(0); }
      if (nh) hexForces(dts); if (ns) shellForces(dts); if (nb) beamForces(dts); if (spr.length) springForces();
      let dt = Math.min(fe.dt, dtHex, dtShell); if (opt.dtMax && dt > opt.dtMax) dt = opt.dtMax; if (fe.nStep && dt > 1.05 * dtPrev) dt = 1.05 * dtPrev;
      contactForces(dt);
      const dta = fe.nStep ? 0.5 * (dt + dtPrev) : dt, al = fe.damping; let wExt = 0, wCon = 0, wDmp = 0, wHg = 0, wQ = 0, ke = 0;
      for (let i = 0; i < nn; i++) {
        const m = ms[i], i3 = 3 * i;
        for (let k = 0; k < 3; k++) {
          const j = i3 + k; if (fixed[6 * i + k]) { a[j] = 0; v[j] = 0; continue; }
          const fe_ = m * gv[k] + fl[j], vo = v[j], fd = -al * m * vo, acc = (f[j] + fh[j] + fq[j] + fc[j] + fe_ + fd) / m, vn = vo + dta * acc, vm = 0.5 * (vo + vn);
          a[j] = acc; v[j] = vn; wExt += fe_ * vm; wCon -= fc[j] * vm; wHg -= fh[j] * vm; wQ -= fq[j] * vm; wDmp -= fd * vm; ke += 0.5 * m * vn * vn;
        }
      }
      if (hasRot) for (let i = 0; i < nn; i++) {
        const I = In[i]; if (!(I > 0)) continue; const i3 = 3 * i;
        for (let k = 0; k < 3; k++) { const j = i3 + k; if (fixed[6 * i + 3 + k]) { w[j] = 0; continue; } const wo = w[j], md = -al * I * wo, wn = wo + (dta * (mo[j] + mc[j] + md)) / I, wm = 0.5 * (wo + wn); w[j] = wn; wCon -= mc[j] * wm; wDmp -= md * wm; ke += 0.5 * I * wn * wn; }
        // rotate the nodal triad by the exponential map of ω·Δt
        const wx = w[i3], wy = w[i3 + 1], wz = w[i3 + 2], wn2 = wx * wx + wy * wy + wz * wz;
        if (wn2 > 0) {
          const wn = Math.sqrt(wn2), th = wn * dt, s = Math.sin(th), c1 = 1 - Math.cos(th), kx = wx / wn, ky = wy / wn, kz = wz / wn, r = 9 * i;
          const q00 = 1 - c1 * (ky * ky + kz * kz), q01 = -s * kz + c1 * kx * ky, q02 = s * ky + c1 * kx * kz, q10 = s * kz + c1 * kx * ky, q11 = 1 - c1 * (kx * kx + kz * kz), q12 = -s * kx + c1 * ky * kz, q20 = -s * ky + c1 * kx * kz, q21 = s * kx + c1 * ky * kz, q22 = 1 - c1 * (kx * kx + ky * ky);
          for (let cc = 0; cc < 3; cc++) { const b0 = R[r + cc], b1 = R[r + 3 + cc], b2 = R[r + 6 + cc]; R[r + cc] = q00 * b0 + q01 * b1 + q02 * b2; R[r + 3 + cc] = q10 * b0 + q11 * b1 + q12 * b2; R[r + 6 + cc] = q20 * b0 + q21 * b1 + q22 * b2; }
        }
      }
      for (const im of imps) {
        if (!Number.isFinite(im.mass)) { for (let k = 0; k < 3; k++) { wCon -= im.f[k] * im.v[k]; wExt -= im.f[k] * im.v[k]; im.c[k] += dt * im.v[k]; } continue; }
        for (let k = 0; k < 3; k++) { if (!im.free[k]) { im.v[k] = 0; continue; } const vo = im.v[k], vn = vo + (dta * im.f[k]) / im.mass; wCon -= im.f[k] * 0.5 * (vo + vn); im.v[k] = vn; ke += 0.5 * im.mass * vn * vn; im.c[k] += dt * vn; }
      }
      for (let j = 0; j < 3 * nn; j++) x[j] += dt * v[j];
      E.external += dta * wExt; E.contact += dta * wCon; E.damping += dta * wDmp; E.hourglass += dta * wHg; E.viscous += dta * wQ; E.kinetic = ke;
      E.total = ke + E.internal + E.viscous + E.hourglass + E.contact + E.damping - E.external;
      const ref = Math.max(E.initial, E.peak = Math.max(E.peak, ke, E.internal, Math.abs(E.external)));
      E.error = ref > 0 ? Math.abs(E.total - E.initial) / ref : 0;
      dtPrev = dt; fe.dtLast = dt; fe.t += dt; fe.nStep++;
      return dt;
    };
    /** March to tEnd; onStep(fe) is called after every step and may return true to stop. */
    fe.run = (tEnd, onStep, maxSteps = 5e6) => { while (fe.t < tEnd && fe.nStep < maxSteps) { fe.step(); if (!Number.isFinite(E.kinetic)) throw new Error('explicitfe: solution became unstable'); if (onStep && onStep(fe) === true) break; } return fe; };
    /** Total linear momentum [3] of the nodes and free impactors. */
    fe.momentum = () => { const p = [0, 0, 0]; for (let i = 0; i < nn; i++) for (let k = 0; k < 3; k++) p[k] += ms[i] * v[3 * i + k]; for (const im of imps) if (Number.isFinite(im.mass)) for (let k = 0; k < 3; k++) p[k] += im.mass * im.v[k]; return p; };
    /** Evaluate internal nodal forces for the current state without advancing (dts = strain increment time). */
    fe.internalForces = (dts = 0) => { f.fill(0); fh.fill(0); fq.fill(0); if (hasRot) mo.fill(0); if (nh) hexForces(dts); if (ns) shellForces(dts); if (nb) beamForces(dts); if (spr.length) springForces(); for (let j = 0; j < f.length; j++) f[j] += fh[j] + fq[j]; return f; };
    return fe;
  };
  return fe;
}
