// Three-dimensional incompressible Navier–Stokes solver: DNS, LES and RANS on a Cartesian staggered (MAC) grid with an
// immersed boundary. Pure JavaScript on flat typed arrays; no DOM and no dependencies, so it runs unchanged in a
// browser tab, a Web Worker and Node.
//
// METHOD
//   Grid          nx × ny × nz cells, uniform in x and z, uniform or stretched in y (pass `yFaces`). Velocities live on
//                 the cell faces, pressure and eddy viscosity at the cell centres. Two ghost layers on every side.
//   Momentum      conservative (divergence) form, second-order central differences; the advected value can be blended
//                 with QUICK or first-order upwind (`scheme`, `blend`). On a uniform periodic grid the central scheme
//                 conserves kinetic energy exactly in the inviscid limit.
//   Time          fractional step (incremental-pressure projection: the provisional velocity carries −∇pⁿ, so the
//                 immersed surface does not leak at steady state). `time: 'rk3'` is the three-stage low-storage Runge–Kutta scheme of Wray
//                 with a projection per stage; `time: 'ab2'` is variable-step Adams–Bashforth with one projection.
//                 `diffusion: 'implicit'` adds an approximately factored implicit viscous operator in delta form
//                 (I − θΔt δx νδx)(I − θΔt δy νδy)(I − θΔt δz νδz) Δu = Δu_explicit, which removes the viscous step limit
//                 (steady states are unchanged; the time accuracy of the viscous term drops to first order).
//   Pressure      exact discrete projection. The Poisson equation is diagonalised by real Fourier transforms in periodic
//                 directions and cosine transforms in bounded ones (x and z), with a tridiagonal solve in y, or a third
//                 Fourier transform when y is periodic. Power-of-two sizes in x and z use an FFT; other sizes use a dense
//                 O(n²) transform, which is fine up to n ≈ 48 and slow beyond.
//   Turbulence    'dns' (none) · 'les' (Smagorinsky with van Driest damping, or WALE) · 'rans' (Spalart–Allmaras
//                 one-equation model with wall distance, or a mixing-length model). `wallModel: true` replaces the
//                 resolved wall shear by the log-law (equilibrium) value on domain walls and around immersed bodies.
//   Boundaries    per side: 'periodic', 'wall' (no slip, optional wall velocity = moving lid), 'slip', 'inflow'
//                 (uniform vector or a profile function) and 'outflow' (convective, with global mass correction).
//   Bodies        direct-forcing immersed boundary. Without a wall model the surface is no-slip: velocity points inside
//                 the body are held at rest and those within half a cell outside are damped. With `wallModel` the
//                 surface is a slip wall carrying the log-law shear stress: only the wall-normal velocity is removed
//                 near the surface, so an under-resolved boundary layer is not forced to separate. The body is given as a
//                 signed-distance function, a solid-mask function, or a triangle surface that is voxelised by ray
//                 casting along the three axes. The momentum removed by the forcing is the hydrodynamic force.
//                 The wall is located to within about one cell: this is first-order accurate at the surface.
//
// API
//   const s = createSolver({
//     n: [nx, ny, nz], L: [Lx, Ly, Lz], origin: [x0, y0, z0], yFaces: number[ny + 1] (optional, stretched y),
//     nu,                                            // kinematic viscosity (density is 1: forces are per unit density)
//     bc: { x: 'periodic' | [lo, hi], y: ..., z: ... },  // side: 'wall' | 'slip' | 'inflow' | 'outflow' | { type, velocity: [u, v, w] | (x, y, z) => [u, v, w] }
//     model: 'dns' | 'les' | 'rans', sgs: 'smagorinsky' | 'wale', Cs: 0.17 (0.1 with walls), Cw: 0.5 (0.325 with walls), vanDriest: true,
//     rans: 'sa' | 'mixing-length', mixingLengthMax, nuTildeInflow, nuTildeInit: (x, y, z) => ν̃, wallModel: false,
//     scheme: 'central' | 'quick' | 'upwind', blend: 0..1, time: 'rk3' | 'ab2', diffusion: 'explicit' | 'implicit',
//     dt (fixed) or cfl (adaptive), dtMax,
//     forcing: { gradient: [fx, fy, fz] } | { bulk: U, dir: 0 },   // body force, or constant mass flux
//     body: { sdf: (x, y, z) => d } | { mask: (x, y, z) => bool } | { positions: number[], triangles: number[] },
//     init: [u, v, w] | (x, y, z) => [u, v, w], perturb: { amplitude, seed, shape: (x, y, z) => factor },
//     stats: { start: time, fields: false }, historyEvery: 0, residual: false,
//   });
//   s.step()                         advance one time step
//   s.run(n, onProgress, every)      advance n steps; onProgress(fraction, solver) every `every` steps; return false to stop
//   s.diagnostics()                  { t, step, dt, cfl, ke, dissipation, enstrophy, sgsDissipation, divMax, divL2, uMax,
//                                      force: [Fx, Fy, Fz], bulk, residual, wallShear: number[6], yPlus: { mean, max }, nutMax }
//   s.monitor()                      cheap per-step values: { wallShear, wallYPlus, uTau, yPlus, nutMax, sgsDissipation, force, residual, gradient }
//   s.history                        { t, ke, dissipation, enstrophy, sgs, div, residual, fx, fy, fz, bulk, dt, gradient } (plain arrays)
//   s.statistics()                   time- and plane-averaged profiles along y: { y, U, V, W, uu, vv, ww, uv, nut, samples }
//   s.slice(axis, index, quantity, { max })   plane of cell-centred values { x, y, z[row][col] } (plain arrays);
//                                    quantity: 'u' 'v' 'w' 'speed' 'p' 'vorticity' 'wx' 'wy' 'wz' 'q' 'nut' 'chi' 'dist' 'umean' 'vmean' 'wmean' 'pmean'
//   s.field(quantity)                Float32Array of nx·ny·nz cell-centred values, index i + nx·(j + ny·k)
//   s.outline(axis, index)           body outline in that plane as a broken polyline { x, y } (NaN separates segments)
//   s.surfacePressure(max)           pressure samples in the first fluid cells around the body { x, y, z, p }
//   s.setBody(body)                  replace the immersed body
//   s.fields                         { u, v, w, p, nuEff, nuTilde, chi, dist } typed arrays with ghosts; s.index(i, j, k)
//   s.grid                           { nx, ny, nz, x, y, z (centres), xf, yf, zf (faces), dx, dz, dy[] }
//
// Storage is about 110–150 bytes per cell (128³ needs roughly 300 MB).

const G = 2;
const WALL = 1, SLIP = 2, INFLOW = 3, OUTFLOW = 4;
const TYPE = { periodic: 0, wall: WALL, slip: SLIP, inflow: INFLOW, outflow: OUTFLOW };
const KAPPA = 0.41, B_LOG = 5.2, Y_LAM = 11.05;
const SLIP_DEPTH = 3; // cells below the surface that still slide in slip-wall mode

function mulberry(seed) {
  let a = seed >>> 0;
  return () => { a = (a + 0x6d2b79f5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
}

/** Friction velocity from the speed u at wall distance d: u⁺ = y⁺ below y⁺ = 11.05, u⁺ = ln(y⁺)/κ + 5.2 above. */
export function wallFriction(u, d, nu) {
  if (!(u > 0) || !(d > 0)) return 0;
  const Re = (u * d) / nu;
  if (Re <= Y_LAM * Y_LAM) return Math.sqrt((nu * u) / d);
  let ut = u / (Math.log(Re) / KAPPA + B_LOG) * 1.6;
  for (let it = 0; it < 8; it++) {
    const lg = Math.log((d * ut) / nu) / KAPPA + B_LOG, f = ut * lg - u, df = lg + 1 / KAPPA, un = ut - f / df;
    if (!(un > 0)) { ut *= 0.5; continue; }
    if (Math.abs(un - ut) < 1e-10 * ut) { ut = un; break; }
    ut = un;
  }
  return ut;
}

// ---- real transforms that diagonalise the second-difference operator ---------------------------
/** In-place radix-2 complex FFT with precomputed tables. */
function fftRadix2(n) {
  const rev = new Uint32Array(n), cs = new Float64Array(n >> 1), sn = new Float64Array(n >> 1);
  let bits = 0; while (1 << bits < n) bits++;
  for (let i = 0; i < n; i++) { let r = 0; for (let b = 0; b < bits; b++) if (i & (1 << b)) r |= 1 << (bits - 1 - b); rev[i] = r; }
  for (let t = 0; t < n >> 1; t++) { cs[t] = Math.cos((2 * Math.PI * t) / n); sn[t] = Math.sin((2 * Math.PI * t) / n); }
  return (re, im, inverse) => {
    for (let i = 0; i < n; i++) { const j = rev[i]; if (j > i) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
    for (let len = 2; len <= n; len <<= 1) {
      const half = len >> 1, stp = n / len;
      for (let i = 0; i < n; i += len) for (let k = 0, t = 0; k < half; k++, t += stp) {
        const wr = cs[t], wi = inverse ? sn[t] : -sn[t], a = i + k, b = a + half, vr = re[b] * wr - im[b] * wi, vi = re[b] * wi + im[b] * wr;
        re[b] = re[a] - vr; im[b] = im[a] - vi; re[a] += vr; im[a] += vi;
      }
    }
  };
}
/**
 * Transform plan for one direction. periodic: real Fourier basis; otherwise the cosine (DCT-II) basis of a cell-centred
 * Neumann problem. fwd/inv act in place on a Float64Array of length n; lam[m] is the eigenvalue of the unit-spacing
 * second difference for coefficient m.
 */
function makePlan(n, periodic) {
  const lam = new Float64Array(n), pow2 = n >= 4 && (n & (n - 1)) === 0;
  for (let m = 0; m < n; m++) lam[m] = periodic ? -4 * Math.sin((Math.PI * Math.min(m, n - m)) / n) ** 2 : -4 * Math.sin((Math.PI * m) / (2 * n)) ** 2;
  if (pow2) {
    const fft = fftRadix2(n), re = new Float64Array(n), im = new Float64Array(n), h = n >> 1;
    if (periodic) return {
      n, lam,
      fwd(x) { for (let j = 0; j < n; j++) { re[j] = x[j]; im[j] = 0; } fft(re, im, false); x[0] = re[0]; x[h] = re[h]; for (let k = 1; k < h; k++) { x[k] = re[k]; x[n - k] = im[k]; } },
      inv(x) { re[0] = x[0]; im[0] = 0; re[h] = x[h]; im[h] = 0; for (let k = 1; k < h; k++) { re[k] = x[k]; im[k] = x[n - k]; re[n - k] = x[k]; im[n - k] = -x[n - k]; } fft(re, im, true); const s = 1 / n; for (let j = 0; j < n; j++) x[j] = re[j] * s; },
    };
    const c = new Float64Array(n), s = new Float64Array(n);
    for (let k = 0; k < n; k++) { c[k] = Math.cos((Math.PI * k) / (2 * n)); s[k] = Math.sin((Math.PI * k) / (2 * n)); }
    return {
      n, lam,
      fwd(x) { for (let j = 0; j < h; j++) { re[j] = x[2 * j]; re[n - 1 - j] = x[2 * j + 1]; } for (let j = 0; j < n; j++) im[j] = 0; fft(re, im, false); for (let k = 0; k < n; k++) x[k] = re[k] * c[k] + im[k] * s[k]; },
      inv(x) { re[0] = x[0]; im[0] = 0; for (let k = 1; k < n; k++) { const a = x[k], b = x[n - k]; re[k] = a * c[k] + b * s[k]; im[k] = a * s[k] - b * c[k]; } fft(re, im, true); const q = 1 / n; for (let j = 0; j < n; j++) im[j] = re[j] * q; for (let j = 0; j < h; j++) { x[2 * j] = im[j]; x[2 * j + 1] = im[n - 1 - j]; } },
    };
  }
  // dense orthonormal eigenvector matrix (rows are basis vectors), same coefficient ordering as the FFT path
  const Q = new Float64Array(n * n), t = new Float64Array(n);
  for (let m = 0; m < n; m++) for (let j = 0; j < n; j++) {
    let v;
    if (periodic) { const k = Math.min(m, n - m), sinRow = m > n - m; v = m === 0 ? Math.sqrt(1 / n) : 2 * k === n ? Math.sqrt(1 / n) * (j % 2 ? -1 : 1) : Math.sqrt(2 / n) * (sinRow ? Math.sin((2 * Math.PI * k * j) / n) : Math.cos((2 * Math.PI * k * j) / n)); }
    else v = m === 0 ? Math.sqrt(1 / n) : Math.sqrt(2 / n) * Math.cos((Math.PI * m * (j + 0.5)) / n);
    Q[m * n + j] = v;
  }
  return {
    n, lam,
    fwd(x) { for (let m = 0; m < n; m++) { let a = 0; const r = m * n; for (let j = 0; j < n; j++) a += Q[r + j] * x[j]; t[m] = a; } for (let m = 0; m < n; m++) x[m] = t[m]; },
    inv(x) { t.fill(0); for (let m = 0; m < n; m++) { const a = x[m], r = m * n; if (a !== 0) for (let j = 0; j < n; j++) t[j] += Q[r + j] * a; } for (let j = 0; j < n; j++) x[j] = t[j]; },
  };
}

/** Godunov update of |∇d| = 1 from the smallest upwind neighbours a, b, c with spacings ha, hb, hc. */
function eikonal(a, b, c, ha, hb, hc) {
  let t;
  if (a > b) { t = a; a = b; b = t; t = ha; ha = hb; hb = t; }
  if (b > c) { t = b; b = c; c = t; t = hb; hb = hc; hc = t; }
  if (a > b) { t = a; a = b; b = t; t = ha; ha = hb; hb = t; }
  let d = a + ha;
  if (!(d > b)) return d;
  const ia = 1 / (ha * ha), ib = 1 / (hb * hb);
  let A = ia + ib, Bq = a * ia + b * ib, C = a * a * ia + b * b * ib - 1;
  d = (Bq + Math.sqrt(Math.max(Bq * Bq - A * C, 0))) / A;
  if (!(d > c)) return d;
  const ic = 1 / (hc * hc);
  A += ic; Bq += c * ic; C += c * c * ic;
  return (Bq + Math.sqrt(Math.max(Bq * Bq - A * C, 0))) / A;
}

export function createSolver(cfg) {
  const nx = Math.round(cfg.n[0]), ny = Math.round(cfg.n[1]), nz = Math.round(cfg.n[2]);
  if (!(nx >= 4 && ny >= 4 && nz >= 2)) throw new Error('cfd3d: the grid needs at least 4 × 4 × 2 cells');
  const Ldom = cfg.L, org = cfg.origin || [0, 0, 0], nu = cfg.nu;
  const Nx = nx + 2 * G, Ny = ny + 2 * G, Nz = nz + 2 * G, sx = 1, sy = Nx, sz = Nx * Ny, NT = Nx * Ny * Nz;
  const nn = [nx, ny, nz], NN = [Nx, Ny, Nz], st = [sx, sy, sz];
  const model = cfg.model || 'dns', sgs = cfg.sgs || 'smagorinsky', ransModel = cfg.rans || 'sa';
  const variable = model !== 'dns', useSA = model === 'rans' && ransModel === 'sa';
  const scheme = cfg.scheme || 'central', blend = scheme === 'central' ? 0 : Math.min(Math.max(cfg.blend ?? 1, 0), 1), quick = scheme === 'quick';
  const timeScheme = cfg.time || 'rk3', implicit = cfg.diffusion === 'implicit', wallModel = !!cfg.wallModel;

  // ---- boundary conditions ----
  const per = [false, false, false], sides = [];
  for (let a = 0; a < 3; a++) {
    const spec = (cfg.bc || {})['xyz'[a]] ?? 'periodic';
    if (spec === 'periodic' || (Array.isArray(spec) && spec[0] === 'periodic')) { per[a] = true; sides.push({ type: 0 }, { type: 0 }); continue; }
    for (let h = 0; h < 2; h++) {
      const sp = Array.isArray(spec) ? spec[h] : spec, o = typeof sp === 'string' ? { type: sp } : sp, type = TYPE[o.type];
      if (!type) throw new Error(`cfd3d: unknown boundary type "${o.type}"`);
      sides.push({ type, velocity: o.velocity || [0, 0, 0], val: null, tauw: 0, yplus: 0 });
    }
  }
  if (per[1] && cfg.yFaces) throw new Error('cfd3d: a periodic y direction needs a uniform grid');

  // ---- metrics (arrays include ghosts) ----
  const xc = [], xf = [], dcw = [], idc = [], idf = [], i2c = [];
  for (let a = 0; a < 3; a++) {
    const n = nn[a], N = NN[a], w = new Float64Array(N), c = new Float64Array(N), f = new Float64Array(N);
    for (let m = 0; m < n; m++) w[m + G] = a === 1 && cfg.yFaces ? cfg.yFaces[m + 1] - cfg.yFaces[m] : Ldom[a] / n;
    for (let q = 0; q < G; q++) { w[G - 1 - q] = per[a] ? w[G + n - 1 - q] : w[G + q]; w[G + n + q] = per[a] ? w[G + q] : w[G + n - 1 - q]; }
    const x0 = a === 1 && cfg.yFaces ? cfg.yFaces[0] : org[a];
    f[G - 1] = x0; for (let m = G; m < N; m++) f[m] = f[m - 1] + w[m]; for (let m = G - 2; m >= 0; m--) f[m] = f[m + 1] - w[m + 1];
    for (let m = 0; m < N; m++) c[m] = f[m] - 0.5 * w[m];
    const ic = new Float64Array(N), iff = new Float64Array(N), i2 = new Float64Array(N);
    for (let m = 0; m < N; m++) { ic[m] = 1 / w[m]; iff[m] = m < N - 1 ? 1 / (c[m + 1] - c[m]) : 1 / w[m]; i2[m] = m > 0 && m < N - 1 ? 1 / (c[m + 1] - c[m - 1]) : 0.5 / w[m]; }
    xc.push(c); xf.push(f); dcw.push(w); idc.push(ic); idf.push(iff); i2c.push(i2);
  }
  const dx = Ldom[0] / nx, dz = Ldom[2] / nz, idx = 1 / dx, idz = 1 / dz;
  const Ly = xf[1][G + ny - 1] - xf[1][G - 1], volTot = Ldom[0] * Ly * Ldom[2];
  let dyMin = Infinity; for (let j = G; j < G + ny; j++) dyMin = Math.min(dyMin, dcw[1][j]);
  const hMin = Math.min(dx, dyMin, dz);

  // ---- fields ----
  const u = new Float64Array(NT), v = new Float64Array(NT), w = new Float64Array(NT), p = new Float64Array(NT), phi = new Float64Array(NT);
  const vel = [u, v, w];
  let Hn = [new Float64Array(NT), new Float64Array(NT), new Float64Array(NT)], Ho = [new Float64Array(NT), new Float64Array(NT), new Float64Array(NT)];
  const ne = variable ? new Float64Array(NT).fill(nu) : null;
  const nt = useSA ? new Float64Array(NT) : null, ntRhs = useSA ? new Float64Array(NT) : null, vortS = useSA ? new Float32Array(NT) : null;
  const dist = new Float32Array(NT).fill(1e30);
  let chi = null, chiC = null, sd = null, band = null, bandF = null, drag = null, hasBody = false;
  const slipBody = wallModel && variable;
  const prev = cfg.residual ? [new Float32Array(NT), new Float32Array(NT), new Float32Array(NT)] : null;
  const ntIn = cfg.nuTildeInflow ?? 3 * nu;

  // ---- boundary value planes ----
  const oth = (a) => [a === 0 ? 1 : 0, a === 2 ? 1 : 2];
  for (let a = 0; a < 3; a++) if (!per[a]) for (let h = 0; h < 2; h++) {
    const s = sides[2 * a + h], [b, c] = oth(a), Nb = NN[b], Nc = NN[c], xa = h ? xf[a][G + nn[a] - 1] : xf[a][G - 1], pos = [0, 0, 0];
    s.val = [new Float64Array(Nb * Nc), new Float64Array(Nb * Nc), new Float64Array(Nb * Nc)];
    for (let pc = 0; pc < Nc; pc++) for (let pb = 0; pb < Nb; pb++) {
      pos[a] = xa; pos[b] = xc[b][pb]; pos[c] = xc[c][pc];
      const vv = typeof s.velocity === 'function' ? s.velocity(pos[0], pos[1], pos[2]) : s.velocity, q = pb + Nb * pc;
      for (let d = 0; d < 3; d++) s.val[d][q] = s.type === SLIP ? 0 : d === a && s.type === WALL ? 0 : vv[d] || 0;
    }
  }

  /** Fill both ghost layers of a field. kind: 0–2 velocity component, 3 zero-gradient scalar, 4 effective viscosity, 5 SA variable. */
  function fillGhosts(f, kind) {
    for (let a = 0; a < 3; a++) {
      const [b, c] = oth(a), sa = st[a], sb = st[b], sc = st[c], n = nn[a], Nb = NN[b], Nc = NN[c];
      if (per[a]) {
        for (let pc = 0; pc < Nc; pc++) for (let pb = 0; pb < Nb; pb++) {
          const I0 = pb * sb + pc * sc;
          f[I0 + (G - 1) * sa] = f[I0 + (G + n - 1) * sa]; f[I0 + (G - 2) * sa] = f[I0 + (G + n - 2) * sa];
          f[I0 + (G + n) * sa] = f[I0 + G * sa]; f[I0 + (G + n + 1) * sa] = f[I0 + (G + 1) * sa];
        }
        continue;
      }
      for (let h = 0; h < 2; h++) {
        const s = sides[2 * a + h], ty = s.type;
        if (kind === a) {
          const val = s.val[a];
          if (h === 0) for (let pc = 0; pc < Nc; pc++) for (let pb = 0; pb < Nb; pb++) { const fB = pb * sb + pc * sc + (G - 1) * sa, bv = val[pb + Nb * pc]; f[fB] = bv; f[fB - sa] = 2 * bv - f[fB + sa]; }
          else for (let pc = 0; pc < Nc; pc++) for (let pb = 0; pb < Nb; pb++) { const fB = pb * sb + pc * sc + (G + n - 1) * sa, bv = val[pb + Nb * pc]; f[fB] = bv; f[fB + sa] = 2 * bv - f[fB - sa]; f[fB + 2 * sa] = 2 * bv - f[fB - 2 * sa]; }
          continue;
        }
        const dirichlet = kind < 3 ? ty === WALL || ty === INFLOW : kind === 4 ? ty === WALL : kind === 5 ? ty === WALL || ty === INFLOW : false;
        const val = kind < 3 ? s.val[kind] : null, cst = kind === 4 ? nu : kind === 5 ? (ty === INFLOW ? ntIn : 0) : 0;
        const c0 = h ? (G + n - 1) * sa : G * sa, o = h ? sa : -sa;
        for (let pc = 0; pc < Nc; pc++) for (let pb = 0; pb < Nb; pb++) {
          const I = pb * sb + pc * sc + c0;
          if (dirichlet) { const bv = 2 * (val ? val[pb + Nb * pc] : cst); f[I + o] = bv - f[I]; f[I + 2 * o] = bv - f[I - o]; }
          else { f[I + o] = f[I]; f[I + 2 * o] = f[I - o]; }
        }
      }
    }
  }
  const fillVel = () => { fillGhosts(u, 0); fillGhosts(v, 1); fillGhosts(w, 2); };

  // ---- pressure Poisson solver ----
  const planX = makePlan(nx, per[0]), planZ = makePlan(nz, per[2]), planY = per[1] ? makePlan(ny, true) : null;
  const W = new Float64Array(nx * ny * nz), nm = nx * nz, lamXZ = new Float64Array(nm);
  for (let k = 0; k < nz; k++) for (let i = 0; i < nx; i++) lamXZ[i + nx * k] = planX.lam[i] * idx * idx + planZ.lam[k] * idz * idz;
  const ay = new Float64Array(ny), cy = new Float64Array(ny);
  let invB = null;
  if (!per[1]) {
    invB = new Float64Array(nm * ny);
    for (let j = 0; j < ny; j++) { ay[j] = j > 0 ? idc[1][j + G] * idf[1][j + G - 1] : 0; cy[j] = j < ny - 1 ? idc[1][j + G] * idf[1][j + G] : 0; }
    for (let m = 0; m < nm; m++) {
      let Bp = 0;
      for (let j = 0; j < ny; j++) { let Bj = -(ay[j] + cy[j]) + lamXZ[m]; if (j > 0) Bj -= ay[j] * cy[j - 1] * Bp; Bp = 1 / Bj; invB[j * nm + m] = Bp; }
    }
    invB[(ny - 1) * nm] = 0; // the constant mode is singular: pin it
  }
  const lineBuf = [new Float64Array(nx), new Float64Array(ny), new Float64Array(nz)];
  function lines(plan, buf, n, stride, count1, stride1, count2, stride2, inverse) {
    for (let q = 0; q < count2; q++) for (let r = 0; r < count1; r++) {
      const o = q * stride2 + r * stride1;
      for (let m = 0; m < n; m++) buf[m] = W[o + m * stride];
      if (inverse) plan.inv(buf); else plan.fwd(buf);
      for (let m = 0; m < n; m++) W[o + m * stride] = buf[m];
    }
  }
  /** Solve the discrete Poisson equation in place: W (layout i + nx·(k + nz·j)) holds the right-hand side, then the solution. */
  function poisson() {
    lines(planX, lineBuf[0], nx, 1, nz, nx, ny, nm, false);
    lines(planZ, lineBuf[2], nz, nx, nx, 1, ny, nm, false);
    if (per[1]) {
      lines(planY, lineBuf[1], ny, nm, nx, 1, nz, nx, false);
      const idy2 = idc[1][G] * idc[1][G];
      for (let j = 0; j < ny; j++) { const ly = planY.lam[j] * idy2, o = j * nm; for (let m = 0; m < nm; m++) { const l = lamXZ[m] + ly; W[o + m] = l < -1e-14 * idy2 ? W[o + m] / l : 0; } }
      lines(planY, lineBuf[1], ny, nm, nx, 1, nz, nx, true);
    } else {
      for (let j = 1; j < ny; j++) { const a = ay[j], o = j * nm, om = o - nm; for (let m = 0; m < nm; m++) W[o + m] -= a * invB[om + m] * W[om + m]; }
      { const o = (ny - 1) * nm; for (let m = 0; m < nm; m++) W[o + m] *= invB[o + m]; }
      for (let j = ny - 2; j >= 0; j--) { const c = cy[j], o = j * nm, op = o + nm; for (let m = 0; m < nm; m++) W[o + m] = (W[o + m] - c * W[op + m]) * invB[o + m]; }
    }
    lines(planZ, lineBuf[2], nz, nx, nx, 1, ny, nm, true);
    lines(planX, lineBuf[0], nx, 1, nz, nx, ny, nm, true);
  }
  function project(dtk) {
    const idt = 1 / dtk, idyA = idc[1];
    for (let j = 0; j < ny; j++) { const idy = idyA[j + G]; for (let k = 0; k < nz; k++) { let I = G + Nx * (j + G + Ny * (k + G)), o = nx * (k + nz * j); for (let i = 0; i < nx; i++, I++, o++) W[o] = ((u[I] - u[I - 1]) * idx + (v[I] - v[I - sy]) * idy + (w[I] - w[I - sz]) * idz) * idt; } }
    poisson();
    for (let j = 0; j < ny; j++) for (let k = 0; k < nz; k++) { let I = G + Nx * (j + G + Ny * (k + G)), o = nx * (k + nz * j); for (let i = 0; i < nx; i++, I++, o++) phi[I] = W[o]; }
    fillGhosts(phi, 3);
    gradientStep(phi, dtk);
    for (let I = 0; I < NT; I++) p[I] += phi[I];
  }
  /** u ← u − dtk·∇f on every face that carries an unknown. */
  function gradientStep(f, dtk) {
    const idfy = idf[1], gx = dtk * idx, gz = dtk * idz;
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) {
      const gy = dtk * idfy[j]; let I = G + Nx * (j + Ny * k);
      for (let i = 0; i < nx; i++, I++) { const p0 = f[I]; u[I] -= gx * (f[I + 1] - p0); v[I] -= gy * (f[I + sy] - p0); w[I] -= gz * (f[I + sz] - p0); }
    }
  }

  // ---- immersed body ----
  const coordOf = (a, m) => xc[a][m];
  const lowerBound = (arr, lo, hi, val) => { while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m] < val) lo = m + 1; else hi = m; } return lo; };
  /** Extend |d| from frozen interface cells to the whole grid by fast sweeping, keeping the sign of `inside`. */
  function sweepDistance(d, frozen) {
    const wy = dcw[1];
    for (let pass = 0; pass < 2; pass++) for (let s = 0; s < 8; s++) {
      const di = s & 1 ? -1 : 1, dj = s & 2 ? -1 : 1, dk = s & 4 ? -1 : 1;
      for (let kk = 0; kk < nz; kk++) { const k = G + (dk > 0 ? kk : nz - 1 - kk); for (let jj = 0; jj < ny; jj++) { const j = G + (dj > 0 ? jj : ny - 1 - jj), hy = wy[j]; for (let ii = 0; ii < nx; ii++) {
        const i = G + (di > 0 ? ii : nx - 1 - ii), I = i + Nx * (j + Ny * k);
        if (frozen[I]) continue;
        const a = Math.min(i > G ? d[I - 1] : 1e30, i < G + nx - 1 ? d[I + 1] : 1e30), b = Math.min(j > G ? d[I - sy] : 1e30, j < G + ny - 1 ? d[I + sy] : 1e30), c = Math.min(k > G ? d[I - sz] : 1e30, k < G + nz - 1 ? d[I + sz] : 1e30);
        if (a > 1e29 && b > 1e29 && c > 1e29) continue;
        const val = eikonal(a, b, c, dx, hy, dz);
        if (val < d[I]) d[I] = val;
      } } }
    }
  }
  function voxelize(body) {
    const P = body.positions, flat = typeof P[0] === 'number', T = body.triangles, nV = flat ? P.length / 3 : P.length;
    const tri = T ? (typeof T[0] === 'number' ? T : T.flat()) : null, nT = tri ? Math.floor(tri.length / 3) : Math.floor(nV / 3);
    const px = (i, a) => (flat ? P[3 * i + a] : P[i][a]);
    const votes = new Uint8Array(NT), dax = [new Float32Array(NT).fill(1e30), new Float32Array(NT).fill(1e30), new Float32Array(NT).fill(1e30)];
    for (let a = 0; a < 3; a++) {
      const b = (a + 1) % 3, c = (a + 2) % 3, nb = nn[b], nc = nn[c], na = nn[a], sa = st[a], sb = st[b], sc = st[c], ca = xc[a], cb = xc[b], cc = xc[c];
      const eb = 3.711e-7 * dcw[b][G], ec = 7.193e-7 * dcw[c][G], rayOf = [], tOf = [];
      for (let t = 0; t < nT; t++) {
        const i0 = tri ? tri[3 * t] : 3 * t, i1 = tri ? tri[3 * t + 1] : 3 * t + 1, i2 = tri ? tri[3 * t + 2] : 3 * t + 2;
        const b0 = px(i0, b), b1 = px(i1, b), b2 = px(i2, b), c0 = px(i0, c), c1 = px(i1, c), c2 = px(i2, c), a0 = px(i0, a), a1 = px(i1, a), a2 = px(i2, a);
        const ar = (b1 - b0) * (c2 - c0) - (b2 - b0) * (c1 - c0);
        if (!(Math.abs(ar) > 1e-300)) continue;
        const pb0 = lowerBound(cb, G, G + nb, Math.min(b0, b1, b2) - eb), pb1 = lowerBound(cb, G, G + nb, Math.max(b0, b1, b2) - eb);
        const pc0 = lowerBound(cc, G, G + nc, Math.min(c0, c1, c2) - ec), pc1 = lowerBound(cc, G, G + nc, Math.max(c0, c1, c2) - ec);
        for (let pc = pc0; pc < pc1; pc++) for (let pb = pb0; pb < pb1; pb++) {
          const rb = cb[pb] + eb, rc = cc[pc] + ec, w1 = ((rb - b0) * (c2 - c0) - (b2 - b0) * (rc - c0)) / ar, w2 = ((b1 - b0) * (rc - c0) - (rb - b0) * (c1 - c0)) / ar, w0 = 1 - w1 - w2;
          if (w0 >= 0 && w1 >= 0 && w2 >= 0) { rayOf.push(pb - G + nb * (pc - G)); tOf.push(w0 * a0 + w1 * a1 + w2 * a2); }
        }
      }
      const nRay = nb * nc, start = new Int32Array(nRay + 1), nH = rayOf.length;
      for (let q = 0; q < nH; q++) start[rayOf[q] + 1]++;
      for (let r = 0; r < nRay; r++) start[r + 1] += start[r];
      const fillAt = start.slice(0, nRay), ts = new Float64Array(nH);
      for (let q = 0; q < nH; q++) ts[fillAt[rayOf[q]]++] = tOf[q];
      for (let r = 0; r < nRay; r++) {
        const s0 = start[r], s1 = start[r + 1]; if (s1 === s0) continue;
        for (let q = s0 + 1; q < s1; q++) { const tv = ts[q]; let m = q - 1; while (m >= s0 && ts[m] > tv) { ts[m + 1] = ts[m]; m--; } ts[m + 1] = tv; }
        const I0 = ((r % nb) + G) * sb + (Math.floor(r / nb) + G) * sc;
        for (let q = s0; q + 1 < s1; q += 2) { const m0 = lowerBound(ca, G, G + na, ts[q]), m1 = lowerBound(ca, G, G + na, ts[q + 1]); for (let m = m0; m < m1; m++) votes[I0 + m * sa]++; }
        for (let q = s0; q < s1; q++) {
          const m1 = lowerBound(ca, G, G + na, ts[q]), d = dax[a];
          if (m1 < G + na) { const I = I0 + m1 * sa, dd = ca[m1] - ts[q]; if (dd < d[I]) d[I] = dd; }
          if (m1 > G) { const I = I0 + (m1 - 1) * sa, dd = ts[q] - ca[m1 - 1]; if (dd < d[I]) d[I] = dd; }
        }
      }
    }
    const frozen = new Uint8Array(NT), d = new Float32Array(NT).fill(1e30);
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) for (let i = G; i < G + nx; i++) {
      const I = i + Nx * (j + Ny * k), a = dax[0][I], b = dax[1][I], c = dax[2][I];
      const ha = a < 1.01 * dx, hb = b < 1.01 * dcw[1][j], hc = c < 1.01 * dz;
      if (ha || hb || hc) { d[I] = 1 / Math.sqrt((ha ? 1 / (a * a + 1e-30) : 0) + (hb ? 1 / (b * b + 1e-30) : 0) + (hc ? 1 / (c * c + 1e-30) : 0)); frozen[I] = 1; }
    }
    sweepDistance(d, frozen);
    for (let I = 0; I < NT; I++) if (votes[I] >= 2) d[I] = -d[I];
    return d;
  }
  function maskToDistance(inside) {
    const frozen = new Uint8Array(NT), d = new Float32Array(NT).fill(1e30);
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) for (let i = G; i < G + nx; i++) {
      const I = i + Nx * (j + Ny * k), f = inside[I]; let h = 1e30;
      if ((i > G && inside[I - 1] !== f) || (i < G + nx - 1 && inside[I + 1] !== f)) h = Math.min(h, dx);
      if ((j > G && inside[I - sy] !== f) || (j < G + ny - 1 && inside[I + sy] !== f)) h = Math.min(h, dcw[1][j]);
      if ((k > G && inside[I - sz] !== f) || (k < G + nz - 1 && inside[I + sz] !== f)) h = Math.min(h, dz);
      if (h < 1e29) { d[I] = 0.5 * h; frozen[I] = 1; }
    }
    sweepDistance(d, frozen);
    for (let I = 0; I < NT; I++) if (inside[I]) d[I] = -d[I];
    return d;
  }
  function wallDistance() {
    for (let k = 0; k < Nz; k++) for (let j = 0; j < Ny; j++) for (let i = 0; i < Nx; i++) {
      const I = i + Nx * (j + Ny * k), m = [i, j, k]; let d = 1e30;
      for (let a = 0; a < 3; a++) if (!per[a]) {
        if (sides[2 * a].type === WALL) d = Math.min(d, Math.abs(xc[a][m[a]] - xf[a][G - 1]));
        if (sides[2 * a + 1].type === WALL) d = Math.min(d, Math.abs(xf[a][G + nn[a] - 1] - xc[a][m[a]]));
      }
      if (hasBody) d = Math.min(d, Math.max(sd[I], 0));
      dist[I] = d;
    }
  }
  function setBody(body) {
    if (!body) { hasBody = false; chi = chiC = sd = band = bandF = drag = null; wallDistance(); return; }
    hasBody = true;
    sd = new Float32Array(NT); chiC = new Float32Array(NT); chi = [new Float32Array(NT), new Float32Array(NT), new Float32Array(NT)]; drag = slipBody ? new Float32Array(NT) : null;
    const dF = [new Float32Array(NT), new Float32Array(NT), new Float32Array(NT)]; // signed distance at the velocity points
    if (typeof body.sdf === 'function') {
      for (let k = 0; k < Nz; k++) for (let j = 0; j < Ny; j++) for (let i = 0; i < Nx; i++) {
        const I = i + Nx * (j + Ny * k), x = xc[0][i], y = xc[1][j], z = xc[2][k];
        sd[I] = body.sdf(x, y, z); dF[0][I] = body.sdf(xf[0][i], y, z); dF[1][I] = body.sdf(x, xf[1][j], z); dF[2][I] = body.sdf(x, y, xf[2][k]);
      }
    } else {
      let d;
      if (typeof body.mask === 'function') { const inside = new Uint8Array(NT); for (let k = 0; k < Nz; k++) for (let j = 0; j < Ny; j++) for (let i = 0; i < Nx; i++) inside[i + Nx * (j + Ny * k)] = body.mask(xc[0][i], xc[1][j], xc[2][k]) ? 1 : 0; d = maskToDistance(inside); }
      else if (body.positions) d = voxelize(body);
      else throw new Error('cfd3d: body needs sdf, mask or positions/triangles');
      sd.set(d); fillGhosts(sd, 3);
      for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) for (let i = G; i < G + nx; i++) { const I = i + Nx * (j + Ny * k); dF[0][I] = 0.5 * (sd[I] + sd[I + 1]); dF[1][I] = 0.5 * (sd[I] + sd[I + sy]); dF[2][I] = 0.5 * (sd[I] + sd[I + sz]); }
      for (let a = 0; a < 3; a++) fillGhosts(dF[a], 3);
    }
    // No-slip mode: a velocity point is held at rest anywhere inside the body and damped within half a cell outside.
    // Slip mode (wall model on): the wall-normal velocity component is removed from one cell outside the surface to
    // SLIP_DEPTH cells inside it, the log-law wall stress acts on the first cell outside, and only points deeper than
    // that are held at rest. The fluid just inside the surface therefore slides with the outer flow instead of forming a
    // dead layer that the stair-stepped grid would mix into the near-wall flow (which separates it on a coarse grid).
    const frac = (d, h) => { if (slipBody) return d < -SLIP_DEPTH * h ? 1 : 0; if (d < 0) return 1; const c = 0.5 - d / h; return c < 0 ? 0 : c; };
    bandF = slipBody ? [] : null;
    for (let a = 0; a < 3; a++) {
      const b = a === 0 ? 1 : 0, c = a === 2 ? 1 : 2, sa = st[a], sb = st[b], sc = st[c], ids = [], dat = [], D = dF[a], ch = chi[a];
      for (let k = 0; k < Nz; k++) for (let j = 0; j < Ny; j++) { const h = Math.min(dx, dcw[1][j], dz); for (let i = 0; i < Nx; i++) { const I = i + Nx * (j + Ny * k); ch[I] = frac(D[I], h); } }
      if (!slipBody) continue;
      for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) for (let i = G; i < G + nx; i++) {
        const I = i + Nx * (j + Ny * k), d = D[I], h = Math.min(dx, dcw[1][j], dz);
        if (!(d >= -SLIP_DEPTH * h && d < h)) continue;
        const m = [i, j, k], g = [0, 0, 0];
        if (typeof body.sdf === 'function') { // exact normal by differencing the distance function at the velocity point
          const P = [xc[0][i], xc[1][j], xc[2][k]], e = 1e-3 * h; P[a] = xf[a][m[a]];
          for (let q = 0; q < 3; q++) { P[q] += e; const dp = body.sdf(P[0], P[1], P[2]); P[q] -= 2 * e; g[q] = dp - body.sdf(P[0], P[1], P[2]); P[q] += e; }
        } else {
          g[a] = (sd[I + sa] - sd[I]) * idf[a][m[a]];
          g[b] = 0.5 * (sd[I + sb] - sd[I - sb] + sd[I + sa + sb] - sd[I + sa - sb]) * i2c[b][m[b]];
          g[c] = 0.5 * (sd[I + sc] - sd[I - sc] + sd[I + sa + sc] - sd[I + sa - sc]) * i2c[c][m[c]];
        }
        const gn = Math.hypot(g[0], g[1], g[2]); if (!(gn > 1e-9)) continue;
        ids.push(I); dat.push(d <= 0.5 * h ? 1 : 2 * (1 - d / h), g[0] / gn, g[1] / gn, g[2] / gn, dx * dz * (a === 1 ? 1 / idf[1][j] : dcw[1][j]), d >= 0 ? 1 : 0);
      }
      bandF.push({ I: Int32Array.from(ids), d: Float32Array.from(dat) });
    }
    let nSolid = 0; const list = [];
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) for (let i = G; i < G + nx; i++) {
      const I = i + Nx * (j + Ny * k), h = Math.min(dx, dcw[1][j], dz);
      chiC[I] = sd[I] < 0 ? 1 : Math.max(0.5 - sd[I] / h, 0); if (sd[I] < 0) nSolid++;
      if (Math.abs(sd[I]) < 1.5 * h) list.push(I);
    }
    fillGhosts(chiC, 3);
    band = Int32Array.from(list); solver.solidCells = nSolid;
    wallDistance();
  }
  /** Slip-wall immersed boundary: remove the wall-normal velocity near the surface and apply the wall-model friction outside it. */
  function slipForcing(dtk) {
    for (let a = 0; a < 3; a++) {
      const b = a === 0 ? 1 : 0, c = a === 2 ? 1 : 2, ua = vel[a], ub = vel[b], uc = vel[c], sa = st[a], sb = st[b], sc = st[c], II = bandF[a].I, D = bandF[a].d; let sum = 0;
      for (let q = 0, o = 0; q < II.length; q++, o += 6) {
        const I = II[q], u0 = ua[I], na = D[o + 1 + a], un = u0 * na + 0.25 * (ub[I] + ub[I - sb] + ub[I + sa] + ub[I + sa - sb]) * D[o + 1 + b] + 0.25 * (uc[I] + uc[I - sc] + uc[I + sa] + uc[I + sa - sc]) * D[o + 1 + c];
        const fr = D[o + 5] * dtk * 0.5 * (drag[I] + drag[I + sa]), r = D[o] * un * na + ((u0 - un * na) * fr) / (1 + fr);
        ua[I] = u0 - r; sum += r * D[o + 4];
      }
      impulse[a] += sum;
    }
  }

  // ---- momentum right-hand side ----
  function rhsComp(a, H) {
    const b = a === 0 ? 1 : 0, c = a === 2 ? 1 : 2, ua = vel[a], ub = vel[b], uc = vel[c], sa = st[a], sb = st[b], sc = st[c];
    const wA = dcw[a], icA = idc[a], ifA = idf[a], icB = idc[b], ifB = idf[b], icC = idc[c], ifC = idf[c], q8 = 0.125 * blend, hb = 0.5 * blend, hs = 0.15 * blend, strY = !!cfg.yFaces, qa = quick && !(strY && a === 1), qb = quick && !(strY && b === 1), qc = quick && !(strY && c === 1);
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) {
      let I = G + Nx * (j + Ny * k);
      for (let i = G; i < G + nx; i++, I++) {
        const pa = a === 0 ? i : a === 1 ? j : k, pb = b === 0 ? i : b === 1 ? j : k, pc = c === 0 ? i : c === 1 ? j : k;
        const ida = ifA[pa], idb = icB[pb], idcc = icC[pc], w0 = 0.5 * wA[pa] * ida, w1 = 0.5 * wA[pa + 1] * ida;
        const u0 = ua[I], um = ua[I - sa], up = ua[I + sa], ubp = ua[I + sb], ubm = ua[I - sb], ucp = ua[I + sc], ucm = ua[I - sc];
        const acm = 0.5 * (um + u0), acp = 0.5 * (u0 + up);
        const vbp = w0 * ub[I] + w1 * ub[I + sa], vbm = w0 * ub[I - sb] + w1 * ub[I + sa - sb];
        const vcp = w0 * uc[I] + w1 * uc[I + sa], vcm = w0 * uc[I - sc] + w1 * uc[I + sa - sc];
        let Fap = acp * acp, Fam = acm * acm, Fbp = vbp * 0.5 * (u0 + ubp), Fbm = vbm * 0.5 * (ubm + u0), Fcp = vcp * 0.5 * (u0 + ucp), Fcm = vcm * 0.5 * (ucm + u0);
        if (blend > 0) {
          // QUICK curvature terms assume uniform spacing; in a stretched direction a consistent upwind blend is used instead
          if (qa) { Fap -= q8 * acp * (acp > 0 ? um - 2 * u0 + up : ua[I + 2 * sa] - 2 * up + u0); Fam -= q8 * acm * (acm > 0 ? ua[I - 2 * sa] - 2 * um + u0 : up - 2 * u0 + um); }
          else { const e = a === 1 && quick ? hs : hb; Fap -= e * Math.abs(acp) * (up - u0); Fam -= e * Math.abs(acm) * (u0 - um); }
          if (qb) { Fbp -= q8 * vbp * (vbp > 0 ? ubm - 2 * u0 + ubp : ua[I + 2 * sb] - 2 * ubp + u0); Fbm -= q8 * vbm * (vbm > 0 ? ua[I - 2 * sb] - 2 * ubm + u0 : ubp - 2 * u0 + ubm); }
          else { const e = b === 1 && quick ? hs : hb; Fbp -= e * Math.abs(vbp) * (ubp - u0); Fbm -= e * Math.abs(vbm) * (u0 - ubm); }
          if (qc) { Fcp -= q8 * vcp * (vcp > 0 ? ucm - 2 * u0 + ucp : ua[I + 2 * sc] - 2 * ucp + u0); Fcm -= q8 * vcm * (vcm > 0 ? ua[I - 2 * sc] - 2 * ucm + u0 : ucp - 2 * u0 + ucm); }
          else { const e = c === 1 && quick ? hs : hb; Fcp -= e * Math.abs(vcp) * (ucp - u0); Fcm -= e * Math.abs(vcm) * (u0 - ucm); }
        }
        let diff;
        if (variable) {
          const n0 = ne[I], n1 = ne[I + sa], nab = n0 + n1;
          const tap = 2 * n1 * (up - u0) * icA[pa + 1], tam = 2 * n0 * (u0 - um) * icA[pa];
          const tbp = 0.25 * (nab + ne[I + sb] + ne[I + sa + sb]) * ((ubp - u0) * ifB[pb] + (ub[I + sa] - ub[I]) * ida);
          const tbm = 0.25 * (nab + ne[I - sb] + ne[I + sa - sb]) * ((u0 - ubm) * ifB[pb - 1] + (ub[I + sa - sb] - ub[I - sb]) * ida);
          const tcp = 0.25 * (nab + ne[I + sc] + ne[I + sa + sc]) * ((ucp - u0) * ifC[pc] + (uc[I + sa] - uc[I]) * ida);
          const tcm = 0.25 * (nab + ne[I - sc] + ne[I + sa - sc]) * ((u0 - ucm) * ifC[pc - 1] + (uc[I + sa - sc] - uc[I - sc]) * ida);
          diff = (tap - tam) * ida + (tbp - tbm) * idb + (tcp - tcm) * idcc;
        } else diff = nu * (((up - u0) * icA[pa + 1] - (u0 - um) * icA[pa]) * ida + ((ubp - u0) * ifB[pb] - (u0 - ubm) * ifB[pb - 1]) * idb + ((ucp - u0) * ifC[pc] - (u0 - ucm) * ifC[pc - 1]) * idcc);
        H[I] = diff - (Fap - Fam) * ida - (Fbp - Fbm) * idb - (Fcp - Fcm) * idcc;
      }
    }
  }

  // ---- implicit viscous operator (approximate factorisation, delta form) ----
  const maxN = Math.max(nx, ny, nz), tA = new Float64Array(maxN), tB = new Float64Array(maxN), tC = new Float64Array(maxN), tD = new Float64Array(maxN), tE = new Float64Array(maxN), tF = new Float64Array(maxN);
  function implicitDir(q, comp, d, fac) {
    const [b, c] = oth(d), sd_ = st[d], sb = st[b], sc = st[c], n = nn[d], normal = d === comp, sC = st[comp];
    const m = per[d] || !normal ? n : n - 1, icD = idc[d], ifD = idf[d], lo = sides[2 * d].type, hi = sides[2 * d + 1].type;
    const dirLo = lo === WALL || lo === INFLOW, dirHi = hi === WALL || hi === INFLOW;
    for (let pc = G; pc < G + nn[c]; pc++) for (let pb = G; pb < G + nn[b]; pb++) {
      const I0 = pb * sb + pc * sc + G * sd_;
      for (let r = 0; r < m; r++) {
        const I = I0 + r * sd_, pd = r + G, nuP = fac * (variable ? 0.5 * (ne[I] + ne[I + sC]) : nu);
        const cm = normal ? nuP * icD[pd] * ifD[pd] : nuP * ifD[pd - 1] * icD[pd], cp = normal ? nuP * icD[pd + 1] * ifD[pd] : nuP * ifD[pd] * icD[pd];
        tA[r] = -cm; tC[r] = -cp; tB[r] = 1 + cm + cp; tD[r] = q[I];
      }
      if (per[d]) { // cyclic system by Sherman–Morrison
        const al = tC[m - 1], be = tA[0], gam = -tB[0];
        tB[0] -= gam; tB[m - 1] -= (al * be) / gam;
        for (let r = 0; r < m; r++) tE[r] = 0; tE[0] = gam; tE[m - 1] = al;
        for (let r = 1; r < m; r++) { const f = tA[r] / tB[r - 1]; tB[r] -= f * tC[r - 1]; tD[r] -= f * tD[r - 1]; tE[r] -= f * tE[r - 1]; }
        tD[m - 1] /= tB[m - 1]; tE[m - 1] /= tB[m - 1];
        for (let r = m - 2; r >= 0; r--) { tD[r] = (tD[r] - tC[r] * tD[r + 1]) / tB[r]; tE[r] = (tE[r] - tC[r] * tE[r + 1]) / tB[r]; }
        const f = (tD[0] + (be * tD[m - 1]) / gam) / (1 + tE[0] + (be * tE[m - 1]) / gam);
        for (let r = 0; r < m; r++) q[I0 + r * sd_] = tD[r] - f * tE[r];
      } else {
        if (!normal) { tB[0] += dirLo ? -tA[0] : tA[0]; tB[m - 1] += dirHi ? -tC[m - 1] : tC[m - 1]; }
        for (let r = 1; r < m; r++) { const f = tA[r] / tB[r - 1]; tB[r] -= f * tC[r - 1]; tD[r] -= f * tD[r - 1]; }
        tD[m - 1] /= tB[m - 1];
        for (let r = m - 2; r >= 0; r--) tD[r] = (tD[r] - tC[r] * tD[r + 1]) / tB[r];
        for (let r = 0; r < m; r++) q[I0 + r * sd_] = tD[r];
      }
    }
  }

  // ---- turbulence closures ----
  const CS = cfg.Cs ?? (sides.some((s) => s.type === WALL) ? 0.1 : 0.17), CW = cfg.Cw ?? (sides.some((s) => s.type === WALL) ? 0.325 : 0.5), vanDriest = cfg.vanDriest !== false, lMax = cfg.mixingLengthMax ?? 1e30;
  const cb1 = 0.1355, sig = 2 / 3, cb2 = 0.622, cv1 = 7.1, cw1 = cb1 / (KAPPA * KAPPA) + (1 + cb2) / sig, cw2 = 0.3, cw3 = 2, cv13 = cv1 ** 3, cw36 = cw3 ** 6;
  let sgsDiss = 0, utauRef = 0, nutMax = 0, ypMean = 0, ypMax = 0;
  /**
   * One pass over the cell centres with the velocity-gradient tensor.
   * mode 1 Smagorinsky → nuEff, 2 WALE → nuEff, 3 mixing length → nuEff, 4 |ω| → out, 5 Q → out, 6 ωx, 7 ωy, 8 ωz → out.
   */
  function gradPass(mode, out) {
    const i2x = i2c[0][G], i2z = i2c[2][G]; let sD = 0, nmx = 0;
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) {
      const idy = idc[1][j], i2y = i2c[1][j], wyj = dcw[1][j], del2 = Math.cbrt(dx * wyj * dz) ** 2; let I = G + Nx * (j + Ny * k);
      for (let i = G; i < G + nx; i++, I++) {
        const g00 = (u[I] - u[I - 1]) * idx, g11 = (v[I] - v[I - sy]) * idy, g22 = (w[I] - w[I - sz]) * idz;
        const g01 = 0.5 * (u[I + sy] + u[I - 1 + sy] - u[I - sy] - u[I - 1 - sy]) * i2y, g02 = 0.5 * (u[I + sz] + u[I - 1 + sz] - u[I - sz] - u[I - 1 - sz]) * i2z;
        const g10 = 0.5 * (v[I + 1] + v[I + 1 - sy] - v[I - 1] - v[I - 1 - sy]) * i2x, g12 = 0.5 * (v[I + sz] + v[I + sz - sy] - v[I - sz] - v[I - sz - sy]) * i2z;
        const g20 = 0.5 * (w[I + 1] + w[I + 1 - sz] - w[I - 1] - w[I - 1 - sz]) * i2x, g21 = 0.5 * (w[I + sy] + w[I + sy - sz] - w[I - sy] - w[I - sy - sz]) * i2y;
        if (mode >= 4) {
          const ox = g21 - g12, oy = g02 - g20, oz = g10 - g01;
          out[I] = mode === 4 ? Math.sqrt(ox * ox + oy * oy + oz * oz) : mode === 5 ? -0.5 * (g00 * g00 + g11 * g11 + g22 * g22) - (g01 * g10 + g02 * g20 + g12 * g21) : mode === 6 ? ox : mode === 7 ? oy : oz;
          continue;
        }
        const s01 = 0.5 * (g01 + g10), s02 = 0.5 * (g02 + g20), s12 = 0.5 * (g12 + g21), SS = g00 * g00 + g11 * g11 + g22 * g22 + 2 * (s01 * s01 + s02 * s02 + s12 * s12);
        let nt_;
        if (mode === 1) {
          const Sm = Math.sqrt(2 * SS); let D = 1;
          if (vanDriest && dist[I] < 1e29) { const yp = utauRef > 0 ? (dist[I] * utauRef) / nu : dist[I] * Math.sqrt(Sm / nu); D = 1 - Math.exp(-yp / 25); }
          nt_ = CS * CS * del2 * D * D * Sm;
        } else if (mode === 2) {
          const a00 = g00 * g00 + g01 * g10 + g02 * g20, a01 = g00 * g01 + g01 * g11 + g02 * g21, a02 = g00 * g02 + g01 * g12 + g02 * g22;
          const a10 = g10 * g00 + g11 * g10 + g12 * g20, a11 = g10 * g01 + g11 * g11 + g12 * g21, a12 = g10 * g02 + g11 * g12 + g12 * g22;
          const a20 = g20 * g00 + g21 * g10 + g22 * g20, a21 = g20 * g01 + g21 * g11 + g22 * g21, a22 = g20 * g02 + g21 * g12 + g22 * g22;
          const tr = (a00 + a11 + a22) / 3, d00 = a00 - tr, d11 = a11 - tr, d22 = a22 - tr, d01 = 0.5 * (a01 + a10), d02 = 0.5 * (a02 + a20), d12 = 0.5 * (a12 + a21);
          const DD = d00 * d00 + d11 * d11 + d22 * d22 + 2 * (d01 * d01 + d02 * d02 + d12 * d12);
          nt_ = DD > 0 ? (CW * CW * del2 * DD * Math.sqrt(DD)) / (SS * SS * Math.sqrt(SS) + DD * Math.sqrt(Math.sqrt(DD)) + 1e-300) : 0;
        } else {
          const Sm = Math.sqrt(2 * SS), d = dist[I]; let l = Math.min(KAPPA * d, lMax);
          if (d < 1e29) { const yp = utauRef > 0 ? (d * utauRef) / nu : d * Math.sqrt(Sm / nu); l *= 1 - Math.exp(-yp / 26); }
          nt_ = l < 1e29 ? l * l * Sm : 0;
        }
        if (hasBody) nt_ *= 1 - chiC[I];
        if (nt_ > nmx) nmx = nt_;
        sD += 2 * nt_ * SS * wyj;
        ne[I] = nu + nt_;
      }
    }
    if (mode < 4) { sgsDiss = sD / (nx * nz * Ly); nutMax = nmx; }
  }
  /** Advance the Spalart–Allmaras variable by dt (first-order upwind convection, point-implicit destruction). */
  function stepSA(dt) {
    gradPass(4, vortS);
    const k2 = KAPPA * KAPPA; let nmx = 0;
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) {
      const idy = idc[1][j], ifp = idf[1][j], ifm = idf[1][j - 1], i2y = i2c[1][j]; let I = G + Nx * (j + Ny * k);
      for (let i = G; i < G + nx; i++, I++) {
        const n0 = nt[I], xm = nt[I - 1], xp = nt[I + 1], ym = nt[I - sy], yp = nt[I + sy], zm = nt[I - sz], zp = nt[I + sz];
        const uc = 0.5 * (u[I] + u[I - 1]), vc = 0.5 * (v[I] + v[I - sy]), wc = 0.5 * (w[I] + w[I - sz]);
        const conv = (uc > 0 ? uc * (n0 - xm) : uc * (xp - n0)) * idx + (vc > 0 ? vc * (n0 - ym) * ifm : vc * (yp - n0) * ifp) + (wc > 0 ? wc * (n0 - zm) : wc * (zp - n0)) * idz;
        const dif = ((nu + 0.5 * (n0 + xp)) * (xp - n0) - (nu + 0.5 * (n0 + xm)) * (n0 - xm)) * idx * idx + ((nu + 0.5 * (n0 + yp)) * (yp - n0) * ifp - (nu + 0.5 * (n0 + ym)) * (n0 - ym) * ifm) * idy + ((nu + 0.5 * (n0 + zp)) * (zp - n0) - (nu + 0.5 * (n0 + zm)) * (n0 - zm)) * idz * idz;
        const gx = 0.5 * (xp - xm) * idx, gy = (yp - ym) * i2y, gz = 0.5 * (zp - zm) * idz;
        const d = Math.max(dist[I], 1e-12), S = vortS[I], X = n0 / nu, X3 = X * X * X, fv1 = X3 / (X3 + cv13), fv2 = 1 - X / (1 + X * fv1);
        let prod = 0, dest = 0;
        if (d < 1e29) {
          const kd2 = k2 * d * d, St = Math.max(S + (n0 * fv2) / kd2, 0.3 * S), r = St > 0 ? Math.min(n0 / (St * kd2), 10) : 10, g = r + cw2 * (r ** 6 - r), fw = g * ((1 + cw36) / (g ** 6 + cw36)) ** (1 / 6);
          prod = cb1 * St * n0; dest = (cw1 * fw * n0) / (d * d);
        } else prod = cb1 * S * n0;
        let nn_ = (n0 + dt * (prod - conv + (dif + cb2 * (gx * gx + gy * gy + gz * gz)) / sig)) / (1 + dt * dest);
        if (!(nn_ > 0)) nn_ = 0;
        if (hasBody) nn_ *= 1 - chiC[I];
        ntRhs[I] = nn_;
      }
    }
    if (implicit && !per[1]) { // wall-normal diffusion of the update treated implicitly (delta form)
      const fac = (1.5 * dt) / sig, dLo = sides[2].type === WALL || sides[2].type === INFLOW, dHi = sides[3].type === WALL || sides[3].type === INFLOW;
      for (let k = G; k < G + nz; k++) for (let i = G; i < G + nx; i++) {
        const I0 = i + Nx * (G + Ny * k);
        for (let r = 0; r < ny; r++) {
          const I = I0 + r * sy, j = r + G, cm = fac * (nu + 0.5 * (nt[I] + nt[I - sy])) * idf[1][j - 1] * idc[1][j], cp = fac * (nu + 0.5 * (nt[I] + nt[I + sy])) * idf[1][j] * idc[1][j];
          tA[r] = -cm; tC[r] = -cp; tB[r] = 1 + cm + cp; tD[r] = ntRhs[I] - nt[I];
        }
        tB[0] += dLo ? -tA[0] : tA[0]; tB[ny - 1] += dHi ? -tC[ny - 1] : tC[ny - 1];
        for (let r = 1; r < ny; r++) { const f = tA[r] / tB[r - 1]; tB[r] -= f * tC[r - 1]; tD[r] -= f * tD[r - 1]; }
        tD[ny - 1] /= tB[ny - 1];
        for (let r = ny - 2; r >= 0; r--) tD[r] = (tD[r] - tC[r] * tD[r + 1]) / tB[r];
        for (let r = 0; r < ny; r++) { const I = I0 + r * sy, vn = nt[I] + tD[r]; ntRhs[I] = vn > 0 ? vn : 0; }
      }
    }
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) { let I = G + Nx * (j + Ny * k); for (let i = G; i < G + nx; i++, I++) { const n0 = ntRhs[I], X = n0 / nu, X3 = X * X * X, t = (n0 * X3) / (X3 + cv13); nt[I] = n0; ne[I] = nu + t; if (t > nmx) nmx = t; } }
    fillGhosts(nt, 5); nutMax = nmx; sgsDiss = 0;
  }
  /** Mean wall shear on domain walls and, with the wall model, the equivalent wall-edge viscosity. */
  function wallTreatment() {
    let tmax = 0;
    for (let a = 0; a < 3; a++) if (!per[a]) for (let h = 0; h < 2; h++) {
      const s = sides[2 * a + h]; if (s.type !== WALL) continue;
      const [b, c] = oth(a), sa = st[a], sb = st[b], sc = st[c], ub = vel[b], uc = vel[c], pa = h ? G + nn[a] - 1 : G, o = h ? sa : -sa, Nb = NN[b];
      let tsum = 0, area = 0, ypS = 0;
      for (let pc = G; pc < G + nn[c]; pc++) for (let pb = G; pb < G + nn[b]; pb++) {
        const I = pb * sb + pc * sc + pa * sa, q = pb + Nb * pc, d = 0.5 * dcw[a][pa], ar = dcw[b][pb] * dcw[c][pc];
        const tb = 0.5 * (ub[I] + ub[I - sb]) - s.val[b][q], tc = 0.5 * (uc[I] + uc[I - sc]) - s.val[c][q], ut = Math.hypot(tb, tc);
        let tau = (nu * ut) / d;
        if (wallModel && variable) {
          const uf = wallFriction(ut, d, nu); tau = uf * uf;
          const nw = ut > 1e-12 ? Math.max((tau * d) / ut, nu) : nu;
          ne[I + o] = 2 * nw - ne[I]; ne[I + 2 * o] = ne[I + o];
        }
        tsum += tau * ar; area += ar; ypS += (Math.sqrt(tau) * d) / nu * ar;
      }
      s.tauw = tsum / area; s.yplus = ypS / area; if (s.tauw > tmax) tmax = s.tauw;
    }
    utauRef = Math.sqrt(tmax);
  }
  /** Wall-layer treatment around the immersed body: y⁺ of the first fluid cells and, optionally, the log-law viscosity. */
  function bodyWallTreatment() {
    if (!hasBody) { ypMean = ypMax = 0; return; }
    let s = 0, n = 0, mx = 0; const i2x = i2c[0][G], i2z = i2c[2][G];
    for (let q = 0; q < band.length; q++) {
      const I = band[q], d0 = sd[I]; if (!(d0 > 0)) continue;
      const j = Math.floor(I / Nx) % Ny, h = Math.min(dx, dcw[1][j], dz), ux = 0.5 * (u[I] + u[I - 1]), uy = 0.5 * (v[I] + v[I - sy]), uz = 0.5 * (w[I] + w[I - sz]);
      if (!slipBody) { // resolved wall: friction velocity from the local velocity gradient
        const d = Math.max(d0, 0.25 * h), yp = (Math.sqrt(((variable ? ne[I] : nu) * Math.hypot(ux, uy, uz)) / d) * d) / nu;
        if (d0 < h) { s += yp; n++; if (yp > mx) mx = yp; }
        continue;
      }
      // wall model: log-law friction velocity from the tangential speed, stored as the drag rate τ_w/(u_t·h) of the near-wall layer
      let gx = (sd[I + 1] - sd[I - 1]) * i2x, gy = (sd[I + sy] - sd[I - sy]) * i2c[1][j], gz = (sd[I + sz] - sd[I - sz]) * i2z; const gn = Math.hypot(gx, gy, gz) || 1; gx /= gn; gy /= gn; gz /= gn;
      const un = ux * gx + uy * gy + uz * gz, ut = Math.sqrt(Math.max(ux * ux + uy * uy + uz * uz - un * un, 0)), d = Math.max(d0, 0.5 * h), uf = wallFriction(ut, d, nu), yp = (uf * d) / nu;
      if (d0 < h) { s += yp; n++; if (yp > mx) mx = yp; }
      drag[I] = ut > 1e-9 ? (uf * uf) / (ut * h) : 0;
    }
    ypMean = n ? s / n : 0; ypMax = mx;
  }
  function turbulence(dt) {
    if (!variable) { if (sides.some((s) => s.type === WALL)) wallTreatment(); bodyWallTreatment(); return; }
    if (useSA) stepSA(dt); else gradPass(model === 'les' ? (sgs === 'wale' ? 2 : 1) : 3, null);
    bodyWallTreatment();
    fillGhosts(ne, 4);
    wallTreatment();
  }

  // ---- time stepping ----
  const cflTarget = cfg.cfl ?? (timeScheme === 'rk3' ? 1.0 : 0.3), vnLimit = timeScheme === 'rk3' ? 0.5 : 0.2, adaptive = !(cfg.dt > 0);
  const fGrad = cfg.forcing?.gradient || [0, 0, 0], bulkTarget = cfg.forcing?.bulk, bulkDir = cfg.forcing?.dir ?? 0;
  let dt = cfg.dt > 0 ? cfg.dt : 0, dtOld = 0, t = 0, nStep = 0, residual = NaN, impliedGradient = fGrad[bulkDir] || 0;
  const force = [0, 0, 0], impulse = [0, 0, 0];
  function stableDt() {
    let cm = 1e-300, dm = 0; const idfy = idf[1], idyA = idc[1], sxz = idx * idx + idz * idz;
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) {
      const iy = idfy[j], dy2 = sxz + idyA[j] * idyA[j]; let I = G + Nx * (j + Ny * k);
      for (let i = G; i < G + nx; i++, I++) {
        const c = Math.abs(u[I]) * idx + Math.abs(v[I]) * iy + Math.abs(w[I]) * idz; if (c > cm) cm = c;
        if (variable) { const d = implicit ? 0 : ne[I] * dy2, e = useSA ? 1.5 * (nu + nt[I]) * (implicit ? sxz : dy2) : 0; if (d > dm) dm = d; if (e > dm) dm = e; }
      }
    }
    if (!variable && !implicit) dm = nu * (sxz + 1 / (dyMin * dyMin));
    let d = cflTarget / cm;
    if (dm > 0) d = Math.min(d, vnLimit / dm);
    if (cfg.dtMax > 0) d = Math.min(d, cfg.dtMax);
    if (dtOld > 0) d = Math.min(d, 1.3 * dtOld);
    return d;
  }
  function outflowUpdate(dtk) {
    let qin = 0, qout = 0, aout = 0;
    for (let a = 0; a < 3; a++) if (!per[a]) for (let h = 0; h < 2; h++) {
      const s = sides[2 * a + h], [b, c] = oth(a), sa = st[a], sb = st[b], sc = st[c], Nb = NN[b], ua = vel[a], val = s.val[a], sg = h ? 1 : -1;
      const fB = (h ? G + nn[a] - 1 : G - 1) * sa, inn = h ? -sa : sa, ic = idc[a][h ? G + nn[a] - 1 : G];
      if (s.type === OUTFLOW) {
        let un = 0, ar = 0;
        for (let pc = G; pc < G + nn[c]; pc++) for (let pb = G; pb < G + nn[b]; pb++) { const A = dcw[b][pb] * dcw[c][pc]; un += sg * ua[pb * sb + pc * sc + fB] * A; ar += A; }
        const Uc = Math.max(un / ar, 0);
        for (let pc = G; pc < G + nn[c]; pc++) for (let pb = G; pb < G + nn[b]; pb++) {
          const I = pb * sb + pc * sc + fB, q = pb + Nb * pc, A = dcw[b][pb] * dcw[c][pc];
          val[q] = ua[I] - dtk * Uc * (ua[I] - ua[I + inn]) * ic;
          qout += sg * val[q] * A; aout += A;
        }
        s._sg = sg;
      } else for (let pc = G; pc < G + nn[c]; pc++) for (let pb = G; pb < G + nn[b]; pb++) qin -= sg * val[pb + Nb * pc] * dcw[b][pb] * dcw[c][pc];
    }
    if (!(aout > 0)) return;
    const corr = (qin - qout) / aout;
    for (let a = 0; a < 3; a++) if (!per[a]) for (let h = 0; h < 2; h++) { const s = sides[2 * a + h]; if (s.type !== OUTFLOW) continue; const val = s.val[a], sg = h ? 1 : -1; for (let q = 0; q < val.length; q++) val[q] += sg * corr; }
  }
  const hasOutflow = sides.some((s) => s.type === OUTFLOW);
  function substep(g, z, dtk, dtFull) {
    for (let a = 0; a < 3; a++) rhsComp(a, Hn[a]);
    for (let a = 0; a < 3; a++) {
      const q = vel[a], h1 = Hn[a], h0 = Ho[a], fa = (g + z) * dtFull * (fGrad[a] || 0), c1 = g * dtFull, c0 = z * dtFull;
      if (implicit) {
        // Δu from the explicit update, smoothed by the factored implicit operator, then added
        for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) { let I = G + Nx * (j + Ny * k); for (let i = G; i < G + nx; i++, I++) h0[I] = c1 * h1[I] + c0 * h0[I] + fa; }
        const fac = 1.5 * dtFull;
        for (let d = 0; d < 3; d++) implicitDir(h0, a, d, fac);
        for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) { let I = G + Nx * (j + Ny * k); for (let i = G; i < G + nx; i++, I++) q[I] += h0[I]; }
      } else for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) { let I = G + Nx * (j + Ny * k); for (let i = G; i < G + nx; i++, I++) q[I] += c1 * h1[I] + c0 * h0[I] + fa; }
    }
    { const s = Hn; Hn = Ho; Ho = s; }
    gradientStep(p, dtk); // incremental-pressure form: the projection below only adds the pressure change
    if (hasBody) for (let a = 0; a < 3; a++) {
      const q = vel[a], ch = chi[a]; let sum = 0;
      for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) {
        const wj = a === 1 ? 1 / idf[1][j] : dcw[1][j]; let I = G + Nx * (j + Ny * k), sr = 0;
        for (let i = G; i < G + nx; i++, I++) { const c = ch[I]; if (c > 0) { const r = c * q[I]; sr += r; q[I] -= r; } }
        sum += sr * wj;
      }
      impulse[a] += sum * dx * dz;
    }
    if (hasBody && slipBody) slipForcing(dtk);
    if (hasOutflow) outflowUpdate(dtk);
    fillVel();
    project(dtk);
    fillVel();
  }
  function bulkVelocity(a) {
    const q = vel[a]; let s = 0;
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) { const wj = a === 1 ? 1 / idf[1][j] : dcw[1][j]; let I = G + Nx * (j + Ny * k), r = 0; for (let i = G; i < G + nx; i++, I++) r += q[I]; s += r * wj; }
    return s / (nx * nz * Ly);
  }
  function step() {
    turbulence(dt > 0 ? dt : 0);
    if (adaptive) { dtOld = dt; dt = stableDt(); }
    if (prev) for (let a = 0; a < 3; a++) prev[a].set(vel[a]);
    impulse[0] = impulse[1] = impulse[2] = 0;
    if (timeScheme === 'rk3') {
      substep(8 / 15, 0, (8 / 15) * dt, dt); substep(5 / 12, -17 / 60, (2 / 15) * dt, dt); substep(3 / 4, -5 / 12, (1 / 3) * dt, dt);
    } else {
      const r = nStep === 0 || !(dtOld > 0) ? 0 : dt / dtOld;
      substep(1 + 0.5 * r, -0.5 * r, dt, dt);
      if (!adaptive) dtOld = dt;
    }
    for (let a = 0; a < 3; a++) force[a] = impulse[a] / dt;
    if (bulkTarget !== undefined) {
      const q = vel[bulkDir], dU = bulkTarget - bulkVelocity(bulkDir);
      for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) { let I = G + Nx * (j + Ny * k); for (let i = G; i < G + nx; i++, I++) q[I] += dU; }
      impliedGradient = dU / dt; fillGhosts(q, bulkDir);
    }
    t += dt; nStep++;
    if (prev) {
      let s = 0;
      for (let a = 0; a < 3; a++) { const q = vel[a], o = prev[a]; for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) { let I = G + Nx * (j + Ny * k); for (let i = G; i < G + nx; i++, I++) { const d = q[I] - o[I]; s += d * d; } } }
      residual = Math.sqrt(s / (3 * nx * ny * nz)) / dt;
    }
    if (stats && t >= stats.start) accumulate();
    if (histEvery > 0 && nStep % histEvery === 0) record();
    solver.time = t; solver.stepCount = nStep; solver.dt = dt;
  }
  function run(n, onProgress, every = 10) {
    for (let s = 0; s < n; s++) {
      step();
      if (!Number.isFinite(u[G + Nx * (G + Ny * G)])) { solver.diverged = true; break; }
      if (onProgress && (s + 1) % every === 0 && onProgress((s + 1) / n, solver) === false) break;
    }
    return solver;
  }

  // ---- diagnostics ----
  function measure() {
    let ke = 0, dmax = 0, d2 = 0, cfl = 0, gr = 0, en = 0, um = 0; const idfy = idf[1], idyA = idc[1];
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) {
      const wy = dcw[1][j], iy = idfy[j], idy = idyA[j], wf = 1 / iy; let I = G + Nx * (j + Ny * k);
      for (let i = G; i < G + nx; i++, I++) {
        const u0 = u[I], v0 = v[I], w0 = w[I];
        ke += (u0 * u0 + w0 * w0) * wy + v0 * v0 * wf;
        const dv = (u0 - u[I - 1]) * idx + (v0 - v[I - sy]) * idy + (w0 - w[I - sz]) * idz, ad = Math.abs(dv); if (ad > dmax) dmax = ad; d2 += dv * dv;
        const c = Math.abs(u0) * idx + Math.abs(v0) * iy + Math.abs(w0) * idz; if (c > cfl) cfl = c;
        const sp = u0 * u0 + v0 * v0 + w0 * w0; if (sp > um) um = sp;
        const uxx = (u0 - u[I - 1]) * idx, uy = (u[I + sy] - u0) * iy, uz = (u[I + sz] - u0) * idz;
        const vx = (v[I + 1] - v0) * idx, vyy = (v0 - v[I - sy]) * idy, vz = (v[I + sz] - v0) * idz;
        const wx = (w[I + 1] - w0) * idx, wyv = (w[I + sy] - w0) * iy, wzz = (w0 - w[I - sz]) * idz;
        gr += (uxx * uxx + uz * uz + vyy * vyy + wx * wx + wzz * wzz) * wy + (uy * uy + vx * vx + vz * vz + wyv * wyv) * wf;
        const ox = wyv - vz, oy = uz - wx, oz = vx - uy;
        en += (oy * oy) * wy + (ox * ox + oz * oz) * wf;
      }
    }
    const nrm = 1 / (nx * nz * Ly);
    return { ke: 0.5 * ke * nrm, dissipation: nu * gr * nrm, enstrophy: 0.5 * en * nrm, divMax: dmax, divL2: Math.sqrt(d2 / (nx * ny * nz)), cfl: cfl * dt, uMax: Math.sqrt(um) };
  }
  const histEvery = cfg.historyEvery | 0;
  const history = { t: [], ke: [], dissipation: [], enstrophy: [], sgs: [], div: [], residual: [], fx: [], fy: [], fz: [], bulk: [], dt: [], gradient: [] };
  function record() {
    const m = measure();
    history.t.push(t); history.ke.push(m.ke); history.dissipation.push(m.dissipation); history.enstrophy.push(m.enstrophy); history.sgs.push(sgsDiss); history.div.push(m.divMax);
    history.residual.push(residual); history.fx.push(force[0]); history.fy.push(force[1]); history.fz.push(force[2]); history.bulk.push(bulkVelocity(bulkDir)); history.dt.push(dt); history.gradient.push(impliedGradient);
  }
  function diagnostics() {
    const m = measure();
    return { t, step: nStep, dt, ...m, sgsDissipation: sgsDiss, force: force.slice(), bulk: bulkVelocity(bulkDir), residual, gradient: impliedGradient, wallShear: sides.map((s) => s.tauw || 0), wallYPlus: sides.map((s) => s.yplus || 0), yPlus: { mean: ypMean, max: ypMax }, nutMax, uTau: utauRef };
  }

  // ---- statistics ----
  const stats = cfg.stats ? { start: cfg.stats.start ?? 0, n: 0, U: new Float64Array(ny), V: new Float64Array(ny), W: new Float64Array(ny), uu: new Float64Array(ny), vv: new Float64Array(ny), ww: new Float64Array(ny), uv: new Float64Array(ny), nut: new Float64Array(ny), f: cfg.stats.fields ? [new Float32Array(nx * ny * nz), new Float32Array(nx * ny * nz), new Float32Array(nx * ny * nz), new Float32Array(nx * ny * nz)] : null } : null;
  function accumulate() {
    const inv = 1 / (nx * nz), f = stats.f;
    for (let j = G; j < G + ny; j++) {
      let a = 0, b = 0, c = 0, aa = 0, bb = 0, cc = 0, ab = 0, e = 0;
      for (let k = G; k < G + nz; k++) { let I = G + Nx * (j + Ny * k), o = nx * (j - G + ny * (k - G)); for (let i = G; i < G + nx; i++, I++, o++) {
        const uc = 0.5 * (u[I] + u[I - 1]), vc = 0.5 * (v[I] + v[I - sy]), wc = 0.5 * (w[I] + w[I - sz]);
        a += uc; b += vc; c += wc; aa += uc * uc; bb += vc * vc; cc += wc * wc; ab += uc * vc; if (variable) e += ne[I] - nu;
        if (f) { f[0][o] += uc; f[1][o] += vc; f[2][o] += wc; f[3][o] += p[I]; }
      } }
      const q = j - G; stats.U[q] += a * inv; stats.V[q] += b * inv; stats.W[q] += c * inv; stats.uu[q] += aa * inv; stats.vv[q] += bb * inv; stats.ww[q] += cc * inv; stats.uv[q] += ab * inv; stats.nut[q] += e * inv;
    }
    stats.n++;
  }
  function statistics() {
    if (!stats || !stats.n) return null;
    const m = 1 / stats.n, arr = (A) => Array.from(A, (x) => x * m), U = arr(stats.U), V = arr(stats.V), Wm = arr(stats.W);
    return { y: Array.from(xc[1].subarray(G, G + ny)), U, V, W: Wm, uu: arr(stats.uu).map((x, j) => x - U[j] * U[j]), vv: arr(stats.vv).map((x, j) => x - V[j] * V[j]), ww: arr(stats.ww).map((x, j) => x - Wm[j] * Wm[j]), uv: arr(stats.uv).map((x, j) => x - U[j] * V[j]), nut: arr(stats.nut), samples: stats.n };
  }

  // ---- field access ----
  let scratch = null, scratchKey = '';
  function cellFn(q) {
    const derived = { vorticity: 4, q: 5, wx: 6, wy: 7, wz: 8 }[q];
    if (derived) { const key = q + ':' + nStep; if (!scratch) scratch = new Float32Array(NT); if (scratchKey !== key) { gradPass(derived, scratch); scratchKey = key; } const s = scratch; return (I) => s[I]; }
    const mi = { umean: 0, vmean: 1, wmean: 2, pmean: 3 }[q];
    if (mi !== undefined) { if (!stats?.f || !stats.n) throw new Error('cfd3d: mean fields were not accumulated'); const f = stats.f[mi], m = 1 / stats.n; return (I, i, j, k) => f[i - G + nx * (j - G + ny * (k - G))] * m; }
    switch (q) {
      case 'u': return (I) => 0.5 * (u[I] + u[I - 1]);
      case 'v': return (I) => 0.5 * (v[I] + v[I - sy]);
      case 'w': return (I) => 0.5 * (w[I] + w[I - sz]);
      case 'speed': return (I) => Math.hypot(0.5 * (u[I] + u[I - 1]), 0.5 * (v[I] + v[I - sy]), 0.5 * (w[I] + w[I - sz]));
      case 'p': return (I) => p[I];
      case 'nut': return (I) => (variable ? (ne[I] - nu) / nu : 0);
      case 'nutilde': return (I) => (nt ? nt[I] / nu : 0);
      case 'chi': return (I) => (hasBody ? chiC[I] : 0);
      case 'sd': return (I) => (hasBody ? sd[I] : 1e30);
      case 'dist': return (I) => dist[I];
      default: throw new Error(`cfd3d: unknown quantity "${q}"`);
    }
  }
  const axisOf = (a) => (typeof a === 'number' ? a : 'xyz'.indexOf(a));
  const planeAxes = (a) => (a === 2 ? [0, 1] : a === 1 ? [0, 2] : [2, 1]);
  function slice(axis, index, quantity, opt = {}) {
    const a = axisOf(axis), [h, vv] = planeAxes(a), max = opt.max || 80, f = cellFn(quantity);
    const sh = Math.max(1, Math.ceil(nn[h] / max)), sv = Math.max(1, Math.ceil(nn[vv] / max)), pa = G + Math.min(Math.max(Math.round(index), 0), nn[a] - 1);
    const xs = [], ys = [], z = [], m = [0, 0, 0]; m[a] = pa;
    for (let ph = G + (sh >> 1); ph < G + nn[h]; ph += sh) xs.push(xc[h][ph]);
    for (let pv = G + (sv >> 1); pv < G + nn[vv]; pv += sv) {
      ys.push(xc[vv][pv]); const row = []; m[vv] = pv;
      for (let ph = G + (sh >> 1); ph < G + nn[h]; ph += sh) { m[h] = ph; row.push(f(m[0] + Nx * (m[1] + Ny * m[2]), m[0], m[1], m[2])); }
      z.push(row);
    }
    return { x: xs, y: ys, z, axis: a, position: xc[a][pa] };
  }
  function field(quantity) {
    const f = cellFn(quantity), out = new Float32Array(nx * ny * nz);
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) for (let i = G; i < G + nx; i++) out[i - G + nx * (j - G + ny * (k - G))] = f(i + Nx * (j + Ny * k), i, j, k);
    return out;
  }
  function outline(axis, index) {
    const X = [], Y = []; if (!hasBody) return { x: X, y: Y };
    const a = axisOf(axis), [h, vv] = planeAxes(a), pa = G + Math.min(Math.max(Math.round(index), 0), nn[a] - 1), sh = st[h], sv = st[vv];
    const cross = (x0, y0, d0, x1, y1, d1) => { const t = d0 / (d0 - d1); return [x0 + t * (x1 - x0), y0 + t * (y1 - y0)]; };
    for (let pv = G; pv < G + nn[vv] - 1; pv++) for (let ph = G; ph < G + nn[h] - 1; ph++) {
      const I = pa * st[a] + ph * sh + pv * sv, d00 = sd[I], d10 = sd[I + sh], d01 = sd[I + sv], d11 = sd[I + sh + sv], x0 = xc[h][ph], x1 = xc[h][ph + 1], y0 = xc[vv][pv], y1 = xc[vv][pv + 1], pts = [];
      if (d00 < 0 !== d10 < 0) pts.push(cross(x0, y0, d00, x1, y0, d10));
      if (d10 < 0 !== d11 < 0) pts.push(cross(x1, y0, d10, x1, y1, d11));
      if (d01 < 0 !== d11 < 0) pts.push(cross(x0, y1, d01, x1, y1, d11));
      if (d00 < 0 !== d01 < 0) pts.push(cross(x0, y0, d00, x0, y1, d01));
      for (let q = 0; q + 1 < pts.length; q += 2) { X.push(pts[q][0], pts[q + 1][0], NaN); Y.push(pts[q][1], pts[q + 1][1], NaN); }
    }
    return { x: X, y: Y };
  }
  function surfacePressure(max = 2000) {
    const out = { x: [], y: [], z: [], p: [] }; if (!hasBody) return out;
    const cand = []; for (let q = 0; q < band.length; q++) { const I = band[q], j = Math.floor(I / Nx) % Ny; if (sd[I] > 0 && sd[I] < Math.min(dx, dcw[1][j], dz)) cand.push(I); }
    const stp = Math.max(1, Math.ceil(cand.length / max));
    for (let q = 0; q < cand.length; q += stp) { const I = cand[q], i = I % Nx, j = Math.floor(I / Nx) % Ny, k = Math.floor(I / (Nx * Ny)); out.x.push(xc[0][i]); out.y.push(xc[1][j]); out.z.push(xc[2][k]); out.p.push(p[I]); }
    return out;
  }

  const solver = {
    grid: { nx, ny, nz, x: Array.from(xc[0].subarray(G, G + nx)), y: Array.from(xc[1].subarray(G, G + ny)), z: Array.from(xc[2].subarray(G, G + nz)), xf: Array.from(xf[0].subarray(G - 1, G + nx)), yf: Array.from(xf[1].subarray(G - 1, G + ny)), zf: Array.from(xf[2].subarray(G - 1, G + nz)), dx, dz, dy: Array.from(dcw[1].subarray(G, G + ny)), ghost: G },
    fields: { u, v, w, p, nuEff: ne, nuTilde: nt, get chi() { return chiC; }, get sd() { return sd; }, dist },
    index: (i, j, k) => i + G + Nx * (j + G + Ny * (k + G)),
    time: 0, stepCount: 0, dt: 0, diverged: false, solidCells: 0,
    step, run, diagnostics, history,
    /** Cheap per-step monitor (no field pass): wall shear and y⁺ per side, body y⁺, eddy-viscosity peak, SGS dissipation, body force. */
    monitor: () => ({ wallShear: sides.map((s) => s.tauw || 0), wallYPlus: sides.map((s) => s.yplus || 0), uTau: utauRef, yPlus: { mean: ypMean, max: ypMax }, nutMax, sgsDissipation: sgsDiss, force: force.slice(), residual, gradient: impliedGradient }), statistics, slice, field, outline, surfacePressure, setBody,
    /** Discrete Poisson solve on a compact right-hand side (layout i + nx·(k + nz·j)); exposed for verification. */
    poisson: (rhs) => { W.set(rhs); poisson(); return Float64Array.from(W); },
  };

  // ---- initial state ----
  if (cfg.body) setBody(cfg.body); else wallDistance();
  {
    const ini = cfg.init, rng = cfg.perturb ? (cfg.perturb.rng || mulberry(cfg.perturb.seed ?? 12345)) : null, amp = cfg.perturb?.amplitude || 0, shape = cfg.perturb?.shape;
    for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) for (let i = G; i < G + nx; i++) {
      const I = i + Nx * (j + Ny * k), x = xc[0][i], y = xc[1][j], z = xc[2][k];
      if (typeof ini === 'function') { u[I] = ini(xf[0][i], y, z)[0]; v[I] = ini(x, xf[1][j], z)[1]; w[I] = ini(x, y, xf[2][k])[2]; }
      else if (ini) { u[I] = ini[0]; v[I] = ini[1]; w[I] = ini[2]; }
      if (rng) { const s = amp * (shape ? shape(x, y, z) : 1); u[I] += s * (2 * rng() - 1); v[I] += s * (2 * rng() - 1); w[I] += s * (2 * rng() - 1); }
    }
    if (hasBody) for (let a = 0; a < 3; a++) { const q = vel[a], ch = chi[a]; for (let I = 0; I < NT; I++) q[I] *= 1 - ch[I]; }
    for (let a = 0; a < 3; a++) if (!per[a]) for (let h = 0; h < 2; h++) { // outflow planes start from the initial field
      const s = sides[2 * a + h]; if (s.type !== OUTFLOW) continue;
      const [b, c] = oth(a), fB = (h ? G + nn[a] - 1 : G - 1) * st[a], inn = h ? -st[a] : st[a];
      for (let pc = G; pc < G + nn[c]; pc++) for (let pb = G; pb < G + nn[b]; pb++) s.val[a][pb + NN[b] * pc] = vel[a][pb * st[b] + pc * st[c] + fB + inn];
    }
    if (hasOutflow) outflowUpdate(0);
    if (nt) { nt.fill(ntIn); if (typeof cfg.nuTildeInit === 'function') for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) for (let i = G; i < G + nx; i++) nt[i + Nx * (j + Ny * k)] = Math.max(cfg.nuTildeInit(xc[0][i], xc[1][j], xc[2][k]), 0); if (hasBody) for (let I = 0; I < NT; I++) nt[I] *= 1 - chiC[I]; fillGhosts(nt, 5); }
    fillVel(); project(1); fillVel(); p.fill(0);
    if (bulkTarget !== undefined) { const q = vel[bulkDir], dU = bulkTarget - bulkVelocity(bulkDir); for (let k = G; k < G + nz; k++) for (let j = G; j < G + ny; j++) { let I = G + Nx * (j + Ny * k); for (let i = G; i < G + nx; i++, I++) q[I] += dU; } fillGhosts(q, bulkDir); }
    if (histEvery > 0) record();
  }
  return solver;
}
