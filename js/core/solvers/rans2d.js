// Two-dimensional body-fitted RANS solver for airfoil sections and flat plates.
// Pure JavaScript on flat typed arrays; no DOM and no dependencies, so it runs unchanged in a browser tab, a Web Worker
// and Node. The iteration allocates no arrays (only the force record and three history entries per step).
//
// METHOD
//   Equations     steady incompressible Reynolds-averaged Navier–Stokes by artificial compressibility
//                 (∂p/∂τ + β ∇·u = 0), with the Spalart–Allmaras one-equation model in its standard form without the
//                 trip term (fully turbulent, 'SA-noft2'), including the negative-ν̃ continuation of Allmaras, Johnson &
//                 Spalart (2012). Unit chord, unit free-stream speed, unit density: ν = 1/Re, Cp = 2p.
//   Grid          structured, cell-centred finite volume. `cGrid` wraps a C-grid with a wake cut around any closed
//                 section (algebraic: wall-normal lines blended into rays to a circular far field 20–30 chords away,
//                 two-zone wall-normal stretching from a first-cell height, tangential smoothing of the outer field);
//                 `plateGrid` and `boxGrid` make the rectangular grids used by the verification cases.
//   Convection    Roe flux-difference splitting for the artificial-compressibility system with third-order upwind-biased
//                 reconstruction along grid lines (κ scheme, κ = 1/2 = QUICK) whose weights account for the grid
//                 stretching; no limiter (the flow is smooth). Spalart–Allmaras convection is first-order upwind on
//                 the same mass flux. The artificial-compressibility parameter β follows the largest speed in the field.
//   Diffusion     full viscous stress on curvilinear cells: face gradients are the average of the Green–Gauss cell
//                 gradients, corrected along the line of centres by the compact two-point difference.
//   Time march    implicit backward Euler in delta form with local time steps: block (3 × 3) tridiagonal solves along
//                 the wall-normal grid lines (the two halves of a wake line are solved as one line through the cut),
//                 symmetric Gauss–Seidel between lines repeated until the update settles, first-order Jacobians
//                 (defect correction). The CFL number follows the residual while it is high (switched evolution
//                 relaxation) and a step that makes the residual jump is taken back. The turbulence equation is solved
//                 the same way after the mean flow (loosely coupled, under-relaxed).
//   Boundaries    no-slip wall with the wall pressure extrapolated from the first two cells; optional log-law wall
//                 function (first cell at y⁺ ≈ 30–200, ν̃ = κ u_τ y there); symmetry; slip wall; far field through the
//                 Roe flux with fixed pressure and free velocity on outflow, and the point-vortex correction for lift.
//
// ACCURACY (NACA 0012, Re = 6 million, against the NASA Turbulence Modeling Resource Spalart–Allmaras results)
//   wall functions, 56 × 22 section grid (1 760 cells): drag +7% at α = 0°, +12% at 10°; lift −2.5% at 10°
//   wall-resolved, 96 × 48 (6 912 cells):               drag +2% at 0°, +5% at 10°, +15% at 15°; lift −1% and −4%
//   wall-resolved, 128 × 64 (11 264 cells):             drag +1% at 0°, +3% at 10°, +9% at 15°; lift −1% and −3%
//   wall-resolved, 192 × 80 (20 480 cells):             drag +0.5% at 0°, +1% at 10°, +5% at 15°; lift −0.4% and −2%
//
// API
//   const sec = resampleSection(X, Y, n)                    any closed section polyline → n cells clustered at the nose and tail
//   const g = cGrid({ X, Y, nWake, nj, h1, radius })        X, Y: trailing edge → lower → leading edge → upper → trailing edge, chord 1
//   plateGrid({ ni, nj, h1 }), boxGrid({ xf, yf, bottom, distort })   rectangular grids of the verification cases
//   g.quality                                               { minOrthogonality_deg, maxAspectRatio, maxGrowth, minArea, cells }
//   const s = createRans2d(g, { Re, alpha, turbulent, wallFunction, clInit, kappa, beta, cflMax, farState, init, slipWall })
//   s.iterate(n, tol)            up to n implicit iterations (stops at residual < tol × initial); returns the residual drop in orders
//   s.forces()                   { cl, cd, cdp, cdf, cm }   (cm about the quarter chord, nose-up positive)
//   s.surface()                  wall distributions { x, y, cp, cf, yplus, utau, side } in grid order
//   s.profile(x, upper)          boundary-layer profile on the wall-normal grid line nearest to x: { yplus, uplus, y, u, nut, utau }
//   s.field(name)                cell values 'p' 'u' 'v' 'speed' 'nut' (ν_t/ν) 'nutilde' 'dist' as Float64Array(ni·nj), index i·nj + j
//   s.setAlpha(a)                change the angle of attack keeping the solution (warm start of a polar)
//   s.history                    { res, cl, cd } plain arrays, one entry per iteration
//   wallSpacing(Re, yplus)       first-cell height that puts the first cell centre at the target y⁺

const KAPPA = 0.41, B_LOG = 5.0, CB1 = 0.1355, SIG = 2 / 3, CB2 = 0.622, CW2 = 0.3, CW3 = 2, CV1 = 7.1, CV2 = 0.7, CV3 = 0.9, CN1 = 16;
const CW1 = CB1 / (KAPPA * KAPPA) + (1 + CB2) / SIG, CW36 = CW3 ** 6, CV13 = CV1 ** 3;
const FAR = 4, WALL = 1, SYM = 2, CUT = 3; // 0 marks an interior face
export const BOUNDARY = { FAR, WALL, SYM, CUT };

/** Height of the first cell whose centre sits at the target y⁺ (flat-plate estimate cf = 0.026 Re^(-1/7)); chord = 1. */
export function wallSpacing(Re, yplus = 1) {
  const cf = 0.026 / Math.max(Re, 1e3) ** (1 / 7);
  return (2 * yplus) / (Math.max(Re, 1e3) * Math.sqrt(cf / 2));
}

const geoSum = (r, m) => (Math.abs(r - 1) < 1e-12 ? m : (r ** m - 1) / (r - 1));
function bisect(f, lo, hi) { for (let k = 0; k < 80; k++) { const m = 0.5 * (lo + hi); if (f(m) > 0) hi = m; else lo = m; } return 0.5 * (lo + hi); }
/**
 * Wall-normal node distances 0 … D for n cells starting with height h1: geometric growth rBL through the boundary layer
 * (thickness about `delta`), then accelerating growth to reach D. Returns { s: Float64Array(n + 1), growth, nBL }.
 */
export function wallStretch(h1, D, n, rBL = 1.2, delta = 0.06) {
  const s = new Float64Array(n + 1);
  if (h1 * n >= D) { for (let j = 0; j <= n; j++) s[j] = (D * j) / n; return { s, growth: 1, nBL: n }; }
  const rg = bisect((r) => h1 * geoSum(r, n) - D, 1, 1e3);
  let r1 = rg, nb = n, q = rg;
  const ratio = (j, r2) => { if (j < nb) return r1; const u = Math.min(1, (j - nb + 1) / 4); return r1 + (r2 - r1) * u * u * (3 - 2 * u); }; // outer zone: ratio r2 reached over four cells
  if (rg > rBL) {
    const nbMax = Math.max(2, Math.round(0.62 * n)), dl = Math.min(delta, 0.5 * D);
    r1 = rBL; nb = Math.ceil(Math.log(1 + (dl * (rBL - 1)) / h1) / Math.log(rBL));
    if (nb > nbMax) { nb = nbMax; r1 = Math.min(rg, bisect((r) => h1 * geoSum(r, nb) - dl, rBL, rg)); }
    const total = (qq) => { let h = h1, t = h1; for (let j = 1; j < n; j++) { h *= ratio(j, qq); t += h; } return t; };
    q = total(r1) >= D ? r1 : bisect((qq) => total(qq) - D, r1, 200);
  }
  let h = h1; s[1] = h1;
  for (let j = 1; j < n; j++) { h *= ratio(j, q); s[j + 1] = s[j] + h; }
  const k = D / s[n]; for (let j = 1; j <= n; j++) s[j] *= k; s[n] = D;
  return { s, growth: r1, nBL: nb };
}

/** Two-sided (Vinokur tanh) distribution of n cells on [0, 1] with first spacing d0 and last spacing d1. */
export function twoSided(n, d0, d1) {
  const u = new Float64Array(n + 1), B = 1 / (n * Math.sqrt(d0 * d1)), A = Math.sqrt(d1 / d0);
  if (!(B > 1.0005)) { for (let k = 0; k <= n; k++) u[k] = k / n; return u; }
  const dl = bisect((d) => Math.sinh(d) / d - B, 1e-6, 40);
  for (let k = 0; k <= n; k++) { const t = 0.5 * (1 + Math.tanh(dl * (k / n - 0.5)) / Math.tanh(dl / 2)); u[k] = t / (A + (1 - A) * t); }
  u[0] = 0; u[n] = 1;
  return u;
}
/**
 * Resample a closed section polyline (trailing edge → lower → leading edge → upper → trailing edge, chord 1) to n cells,
 * clustered in arc length at the leading edge (spacing dLE, default about a third of the nose radius) and at the trailing edge.
 */
export function resampleSection(X, Y, n, { dLE = 0, dTE = 0 } = {}) {
  const m = X.length - 1, t = new Float64Array(m + 1); let kLE = 0, dmax = -1;
  for (let k = 1; k <= m; k++) t[k] = t[k - 1] + Math.hypot(X[k] - X[k - 1], Y[k] - Y[k - 1]);
  for (let k = 0; k <= m; k++) { const d = Math.hypot(X[k] - X[0], Y[k] - Y[0]); if (d > dmax) { dmax = d; kLE = k; } }
  const cr = (A, k, u) => { // Catmull–Rom on the non-uniform arc-length parameter
    const k0 = Math.max(k - 1, 0), k3 = Math.min(k + 2, m), h = t[k + 1] - t[k];
    const m1 = k > 0 ? (A[k + 1] - A[k0]) / (t[k + 1] - t[k0]) : (A[1] - A[0]) / h, m2 = k + 1 < m ? (A[k3] - A[k]) / (t[k3] - t[k]) : (A[m] - A[m - 1]) / h;
    const u2 = u * u, u3 = u2 * u;
    return (2 * u3 - 3 * u2 + 1) * A[k] + (u3 - 2 * u2 + u) * h * m1 + (-2 * u3 + 3 * u2) * A[k + 1] + (u3 - u2) * h * m2;
  };
  let kk = 0;
  const at = (sv) => { if (t[kk] > sv) kk = 0; while (kk < m - 1 && t[kk + 1] < sv) kk++; const u = Math.min(1, Math.max(0, (sv - t[kk]) / (t[kk + 1] - t[kk] || 1))); return [cr(X, kk, u), cr(Y, kk, u)]; };
  // nose radius from the circle through the leading edge and the points about 1% of the arc either side
  const pa = at(t[kLE] - 0.01), pb = at(t[kLE] + 0.01), ax = pa[0] - X[kLE], ay = pa[1] - Y[kLE], bx = pb[0] - X[kLE], by = pb[1] - Y[kLE], cr2 = Math.abs(ax * by - ay * bx);
  const rLE = cr2 > 1e-14 ? (Math.hypot(ax, ay) * Math.hypot(bx, by) * Math.hypot(ax - bx, ay - by)) / (2 * cr2) : 0.02;
  const h = n >> 1, Xo = new Array(n + 1), Yo = new Array(n + 1), sL = t[kLE], sU = t[m] - t[kLE];
  const d0 = dLE > 0 ? dLE : Math.min(0.13 / h, 0.3 * rLE), d1 = dTE > 0 ? dTE : 0.3 / h, uL = twoSided(h, d0 / sL, d1 / sL), uU = twoSided(n - h, d0 / sU, d1 / sU);
  for (let k = 0; k <= n; k++) { const sv = k <= h ? sL * (1 - uL[h - k]) : sL + sU * uU[k - h], p = at(sv); Xo[k] = p[0]; Yo[k] = p[1]; }
  Xo[0] = X[0]; Yo[0] = Y[0]; Xo[n] = X[0]; Yo[n] = Y[0]; Xo[h] = X[kLE]; Yo[h] = Y[kLE];
  return { X: Xo, Y: Yo, rLE, dLE: d0, dTE: d1 };
}

function minArea(x, y, ni, nj) {
  let mn = Infinity; const N1 = nj + 1;
  for (let i = 0; i < ni; i++) for (let j = 0; j < nj; j++) {
    const a = i * N1 + j, b = a + N1, c = b + 1, d = a + 1, ar = 0.5 * ((x[c] - x[a]) * (y[d] - y[b]) - (x[d] - x[b]) * (y[c] - y[a]));
    if (ar < mn) mn = ar;
  }
  return mn;
}
/** Orthogonality, aspect ratio and growth of a structured grid. */
export function gridQuality(g) {
  const { ni, nj, x, y } = g, N1 = nj + 1; let skewMax = 0, skewWall = 0, arMax = 0, grMax = 1, aMin = Infinity;
  for (let i = 0; i < ni; i++) {
    let hPrev = 0;
    for (let j = 0; j < nj; j++) {
      const a = i * N1 + j, b = a + N1, c = b + 1, d = a + 1;
      const ix = 0.5 * (x[b] + x[c] - x[a] - x[d]), iy = 0.5 * (y[b] + y[c] - y[a] - y[d]), jx = 0.5 * (x[d] + x[c] - x[a] - x[b]), jy = 0.5 * (y[d] + y[c] - y[a] - y[b]);
      const li = Math.hypot(ix, iy), lj = Math.hypot(jx, jy), sk = Math.abs(90 - (Math.acos(Math.max(-1, Math.min(1, (ix * jx + iy * jy) / (li * lj || 1)))) * 180) / Math.PI);
      const ar = 0.5 * ((x[c] - x[a]) * (y[d] - y[b]) - (x[d] - x[b]) * (y[c] - y[a]));
      if (sk > skewMax) skewMax = sk; if (j < 4 && g.bottom[i] === WALL && sk > skewWall) skewWall = sk;
      if (Math.max(li / lj, lj / li) > arMax) arMax = Math.max(li / lj, lj / li);
      if (j > 0 && hPrev > 0) { const gr = Math.max(lj / hPrev, hPrev / lj); if (gr > grMax) grMax = gr; }
      hPrev = lj; if (ar < aMin) aMin = ar;
    }
  }
  return { minOrthogonality_deg: 90 - skewMax, minOrthogonalityWall_deg: 90 - skewWall, maxAspectRatio: arMax, maxGrowth: grMax, minArea: aMin, cells: ni * nj };
}

/**
 * C-grid around a closed section. X, Y: nS + 1 surface nodes from the trailing edge along the lower surface to the leading
 * edge and back along the upper surface (first node repeated at the end); chord 1. Grid index i runs from the lower wake
 * exit around the section to the upper wake exit, j from the wall (or the wake cut) to the far field.
 */
export function cGrid({ X, Y, nWake = 24, nj = 40, h1 = 1e-5, radius = 30, wakeLength = 0, rBL = 1.2, delta = 0.06, smooth = 2 } = {}) {
  const nS = X.length - 1, nW = Math.max(2, Math.round(nWake)), ni = nS + 2 * nW, N1 = nj + 1, R = Math.max(5, radius), Lw = wakeLength > 0 ? wakeLength : R;
  const xT = X[0], yT = Y[0], t = new Float64Array(nS + 1);
  for (let k = 1; k <= nS; k++) t[k] = t[k - 1] + Math.hypot(X[k] - X[k - 1], Y[k] - Y[k - 1]);
  const arc = t[nS]; for (let k = 0; k <= nS; k++) t[k] /= arc;
  let nx = new Float64Array(nS + 1), ny = new Float64Array(nS + 1);
  for (let k = 1; k < nS; k++) { const tx = X[k + 1] - X[k - 1], ty = Y[k + 1] - Y[k - 1], l = Math.hypot(tx, ty) || 1; nx[k] = -ty / l; ny[k] = tx / l; }
  { const l0 = Math.hypot(X[1] - X[0], Y[1] - Y[0]), l1 = Math.hypot(X[nS] - X[nS - 1], Y[nS] - Y[nS - 1]); nx[0] = -(Y[1] - Y[0]) / l0; ny[0] = (X[1] - X[0]) / l0; nx[nS] = -(Y[nS] - Y[nS - 1]) / l1; ny[nS] = (X[nS] - X[nS - 1]) / l1; }
  const phiL = Math.atan2(nx[0], -ny[0]), phiU = Math.atan2(nx[nS], ny[nS]), LB = 2, LPHI = 1.5; // lean of the trailing-edge lines, carried into the wake and decaying
  for (let p = 0; p < 4; p++) {
    const ax = nx.slice(), ay = ny.slice();
    for (let k = 1; k < nS; k++) { const u = 0.25 * nx[k - 1] + 0.5 * nx[k] + 0.25 * nx[k + 1], v = 0.25 * ny[k - 1] + 0.5 * ny[k] + 0.25 * ny[k + 1], l = Math.hypot(u, v) || 1; ax[k] = u / l; ay[k] = v / l; }
    nx = ax; ny = ay;
  }
  const thRay = new Float64Array(nS + 1);
  for (let k = 0; k <= nS; k++) {
    const ax = X[k] - xT, ay = Y[k] - yT, b = ax * nx[k] + ay * ny[k], sv = -b + Math.sqrt(Math.max(0, b * b - (ax * ax + ay * ay - R * R)));
    let th = Math.atan2(ay + sv * ny[k], ax + sv * nx[k]); const ref = -Math.PI / 2 - Math.PI * t[k];
    th += 2 * Math.PI * Math.round((ref - th) / (2 * Math.PI)); thRay[k] = th;
  }
  const dTE = Math.max(Math.hypot(X[1] - X[0], Y[1] - Y[0]), Math.hypot(X[nS] - X[nS - 1], Y[nS] - Y[nS - 1])), xw = wallStretch(Math.min(dTE, Lw / nW), Lw, nW, 1.25, 1).s;
  const wall = wallStretch(h1, R, nj, rBL, delta);
  const sw = new Float64Array(nj + 1);
  const build = (w) => {
    const x = new Float64Array((ni + 1) * N1), y = new Float64Array((ni + 1) * N1), th = new Float64Array(nS + 1), dec = new Float64Array(nS + 1); let sum = 0;
    for (let k = 1; k <= nS; k++) { const as = Math.PI * (t[k] - t[k - 1]); dec[k] = Math.max((1 - w) * (thRay[k - 1] - thRay[k]) + w * as, 0.25 * as); sum += dec[k]; }
    th[0] = -Math.PI / 2; for (let k = 1; k <= nS; k++) th[k] = th[k - 1] - (dec[k] * Math.PI) / sum;
    for (let k = 0; k <= nS; k++) {
      const ox = xT + R * Math.cos(th[k]), oy = yT + R * Math.sin(th[k]), dx = ox - X[k], dy = oy - Y[k], D = Math.hypot(dx, dy), ex = dx / D, ey = dy / D, o = (nW + k) * N1;
      for (let j = 0; j <= nj; j++) { const d = (wall.s[j] / R) * D, f = 1 - Math.exp(-d / LB); x[o + j] = X[k] + d * ((1 - f) * nx[k] + f * ex); y[o + j] = Y[k] + d * ((1 - f) * ny[k] + f * ey); }
    }
    for (let m = 1; m <= nW; m++) {
      const dh = Math.min((h1 * xw[m]) / 0.02, 0.004), a = (nW - m) * N1, b = (nW + nS + m) * N1; // the cells on the wake cut thicken downstream
      for (let j = 0; j <= nj; j++) sw[j] = wall.s[j] + dh * j * Math.exp(-wall.s[j]);
      const dk = Math.exp(-xw[m] / LPHI), sl = Math.sin(phiL * dk), cl = Math.cos(phiL * dk), su = Math.sin(phiU * dk), cu = Math.cos(phiU * dk);
      for (let j = 0; j <= nj; j++) { const d = sw[j], f = 1 - Math.exp(-d / LB); x[a + j] = xT + xw[m] + d * (1 - f) * sl; y[a + j] = yT - d * ((1 - f) * cl + f); x[b + j] = xT + xw[m] + d * (1 - f) * su; y[b + j] = yT + d * ((1 - f) * cu + f); }
    }
    return { x, y };
  };
  let grid = null, blend = 0;
  for (const w of [0, 0.25, 0.5, 0.75, 1]) { const g = build(w); if (minArea(g.x, g.y, ni, nj) > 0) { grid = g; blend = w; break; } }
  if (!grid) throw new Error('The section grid could not be generated without inverted cells: use fewer wall-normal cells, a larger first-cell height or more surface cells.');
  // tangential smoothing of the outer field (elliptic along the wrap-around direction; the near-wall layers are untouched)
  for (let p = 0; p < smooth; p++) {
    const x = grid.x.slice(), y = grid.y.slice();
    for (let j = 1; j < nj; j++) {
      const om = 0.5 * Math.max(0, Math.min(1, (wall.s[j] - 0.5) / 2)); if (!(om > 0)) continue;
      for (let i = 1; i < ni; i++) { const q = i * N1 + j; x[q] += om * (0.5 * (grid.x[q - N1] + grid.x[q + N1]) - grid.x[q]); y[q] += om * (0.5 * (grid.y[q - N1] + grid.y[q + N1]) - grid.y[q]); }
    }
    if (minArea(x, y, ni, nj) > 0) { grid.x = x; grid.y = y; } else break;
  }
  const bottom = new Uint8Array(ni); for (let i = 0; i < ni; i++) bottom[i] = i < nW || i >= nW + nS ? CUT : WALL;
  const g = { kind: 'c', ni, nj, x: grid.x, y: grid.y, bottom, nWake: nW, nSurf: nS, radius: R, h1, growth: wall.growth, nBL: wall.nBL, blend };
  g.quality = gridQuality(g);
  return g;
}

/** Rectangular grid from face coordinates xf (ni + 1) and yf (nj + 1); bottom(i) gives the boundary type of each bottom face; distort moves interior nodes. */
export function boxGrid({ xf, yf, bottom = () => FAR, distort = null }) {
  const ni = xf.length - 1, nj = yf.length - 1, N1 = nj + 1, x = new Float64Array((ni + 1) * N1), y = new Float64Array((ni + 1) * N1), bt = new Uint8Array(ni);
  for (let i = 0; i <= ni; i++) for (let j = 0; j <= nj; j++) {
    let px = xf[i], py = yf[j];
    if (distort && i > 0 && i < ni && j > 0 && j < nj) { const d = distort(px, py, i, j); px += d[0]; py += d[1]; }
    x[i * N1 + j] = px; y[i * N1 + j] = py;
  }
  for (let i = 0; i < ni; i++) bt[i] = bottom(i, 0.5 * (xf[i] + xf[i + 1]));
  const g = { kind: 'box', ni, nj, x, y, bottom: bt };
  g.quality = gridQuality(g);
  return g;
}
/** Flat plate of unit length starting at x = 0 with a symmetry strip upstream; first cell height h1, clustered at the leading edge. */
export function plateGrid({ ni = 64, nj = 40, h1 = 1e-4, upstream = 0.25, height = 0.6, rBL = 1.2, delta = 0.03, dxLE = 0.004 } = {}) {
  const nUp = Math.max(4, Math.round(ni * 0.2)), nP = ni - nUp, a = wallStretch(dxLE, 1, nP, 1.15, 1).s, b = wallStretch(dxLE, upstream, nUp, 1.3, upstream).s, xf = [];
  for (let k = nUp; k >= 1; k--) xf.push(-b[k]);
  for (let k = 0; k <= nP; k++) xf.push(a[k]);
  const g = boxGrid({ xf, yf: Array.from(wallStretch(h1, height, nj, rBL, delta).s), bottom: (i) => (i < nUp ? SYM : WALL) });
  g.kind = 'plate'; g.nUp = nUp; g.h1 = h1;
  return g;
}

function inv3(a, o, b, p) {
  const a0 = a[o], a1 = a[o + 1], a2 = a[o + 2], a3 = a[o + 3], a4 = a[o + 4], a5 = a[o + 5], a6 = a[o + 6], a7 = a[o + 7], a8 = a[o + 8];
  const c0 = a4 * a8 - a5 * a7, c1 = a5 * a6 - a3 * a8, c2 = a3 * a7 - a4 * a6, id = 1 / (a0 * c0 + a1 * c1 + a2 * c2);
  b[p] = c0 * id; b[p + 1] = (a2 * a7 - a1 * a8) * id; b[p + 2] = (a1 * a5 - a2 * a4) * id;
  b[p + 3] = c1 * id; b[p + 4] = (a0 * a8 - a2 * a6) * id; b[p + 5] = (a2 * a3 - a0 * a5) * id;
  b[p + 6] = c2 * id; b[p + 7] = (a1 * a6 - a0 * a7) * id; b[p + 8] = (a0 * a4 - a1 * a3) * id;
}

/** Create a solver on a grid from cGrid / plateGrid / boxGrid. */
export function createRans2d(g, o = {}) {
  const { ni, nj, x: gx, y: gy } = g, N1 = nj + 1, nc = ni * nj;
  const Re = o.Re ?? 1e6, nu = 1 / Re, turb = o.turbulent !== false, wf = !!o.wallFunction && turb, kap = o.firstOrder ? -2 : o.kappa ?? 0.5, betaFix = o.beta > 0, betaK = o.betaFactor ?? 1;
  const betaNow = new Float64Array([betaFix ? o.beta : 1]);
  const slip = !!o.slipWall, ntInit = o.nuTildeInit ?? 0.01, ntInf = (o.nuTildeInf ?? 3) * nu, vortex = o.vortex ?? g.kind === 'c', xRef = o.xRef ?? 0.25, farFn = o.farState || null;
  let alpha = o.alpha ?? 0, ca = Math.cos(alpha), sa = Math.sin(alpha), gammaFar = 0;

  // ---- geometry: cells, faces, ghosts ------------------------------------------------------------
  const xc = new Float64Array(nc), yc = new Float64Array(nc), vol = new Float64Array(nc);
  for (let i = 0; i < ni; i++) for (let j = 0; j < nj; j++) {
    const a = i * N1 + j, b = a + N1, c = b + 1, d = a + 1, q = i * nj + j;
    xc[q] = 0.25 * (gx[a] + gx[b] + gx[c] + gx[d]); yc[q] = 0.25 * (gy[a] + gy[b] + gy[c] + gy[d]);
    vol[q] = 0.5 * ((gx[c] - gx[a]) * (gy[d] - gy[b]) - (gx[d] - gx[b]) * (gy[c] - gy[a]));
  }
  let nCut = 0, nBnd = 2 * nj + ni; for (let i = 0; i < ni; i++) { if (g.bottom[i] === CUT) nCut++; else nBnd++; }
  const nf = (ni - 1) * nj + ni * (nj - 1) + nCut / 2 + nBnd, nt = nc + nBnd;
  const fL = new Int32Array(nf), fR = new Int32Array(nf), fLL = new Int32Array(nf).fill(-1), fRR = new Int32Array(nf).fill(-1), fT = new Uint8Array(nf);
  const fnx = new Float64Array(nf), fny = new Float64Array(nf), fdx = new Float64Array(nf), fdy = new Float64Array(nf), fw = new Float64Array(nf), fg = new Float64Array(nf), fxm = new Float64Array(nf), fym = new Float64Array(nf);
  const cellS = new Int32Array(nc).fill(-1), cellN = new Int32Array(nc).fill(-1), cellW = new Int32Array(nc).fill(-1), cellE = new Int32Array(nc).fill(-1);
  const C = (i, j) => i * nj + j;
  let f = 0, gh = nc;
  const geom = (k, a, b, sx, sy, L, Rr, type) => { // face from node a to node b, normal (sx, sy)·|ab| pointing from L to R
    const tx = gx[b] - gx[a], ty = gy[b] - gy[a]; fnx[k] = sx * ty; fny[k] = -sx * tx; if (sy) { fnx[k] = -sy * ty; fny[k] = sy * tx; }
    fxm[k] = 0.5 * (gx[a] + gx[b]); fym[k] = 0.5 * (gy[a] + gy[b]); fL[k] = L; fR[k] = Rr; fT[k] = type;
    if (Rr < nc) { fdx[k] = xc[Rr] - xc[L]; fdy[k] = yc[Rr] - yc[L]; const dl = Math.hypot(fxm[k] - xc[L], fym[k] - yc[L]), dr = Math.hypot(fxm[k] - xc[Rr], fym[k] - yc[Rr]); fw[k] = dr / (dl + dr); }
    else { // ghost at the mirror image of the cell centre in the face
      const s2 = fnx[k] * fnx[k] + fny[k] * fny[k], dn = ((fxm[k] - xc[L]) * fnx[k] + (fym[k] - yc[L]) * fny[k]) / s2, m = type === FAR ? 1 : 2; // a far-field ghost holds the state at the face itself
      fdx[k] = m * dn * fnx[k]; fdy[k] = m * dn * fny[k]; fw[k] = type === FAR ? 0 : 0.5;
    }
    fg[k] = (fnx[k] * fnx[k] + fny[k] * fny[k]) / (fdx[k] * fnx[k] + fdy[k] * fny[k]);
  };
  const ghostOfS = new Int32Array(ni).fill(-1), ghostOfW = new Int32Array(nj), ghostOfE = new Int32Array(nj), ghostOfN = new Int32Array(ni);
  // boundary faces first (their ghosts are needed for the reconstruction stencils)
  for (let j = 0; j < nj; j++) { geom(f, j, j + 1, -1, 0, C(0, j), gh, FAR); cellW[C(0, j)] = f; ghostOfW[j] = gh++; f++; }
  for (let j = 0; j < nj; j++) { geom(f, ni * N1 + j, ni * N1 + j + 1, 1, 0, C(ni - 1, j), gh, FAR); cellE[C(ni - 1, j)] = f; ghostOfE[j] = gh++; f++; }
  for (let i = 0; i < ni; i++) { geom(f, i * N1 + nj, (i + 1) * N1 + nj, 0, 1, C(i, nj - 1), gh, FAR); cellN[C(i, nj - 1)] = f; ghostOfN[i] = gh++; f++; }
  for (let i = 0; i < ni; i++) if (g.bottom[i] !== CUT) { geom(f, i * N1, (i + 1) * N1, 0, -1, C(i, 0), gh, g.bottom[i]); cellS[C(i, 0)] = f; ghostOfS[i] = gh++; f++; }
  const nB = f;
  const below = (i) => (g.bottom[i] === CUT ? C(ni - 1 - i, 0) : ghostOfS[i]);
  for (let i = 0; i < ni; i++) if (g.bottom[i] === CUT && i > ni - 1 - i) { // wake cut: one face per pair of cells
    const m = ni - 1 - i; geom(f, i * N1, (i + 1) * N1, 0, 1, C(m, 0), C(i, 0), 0); fLL[f] = nj > 1 ? C(m, 1) : -1; fRR[f] = nj > 1 ? C(i, 1) : -1; cellS[C(i, 0)] = f; cellS[C(m, 0)] = f; f++;
  }
  for (let i = 1; i < ni; i++) for (let j = 0; j < nj; j++) {
    geom(f, i * N1 + j, i * N1 + j + 1, 1, 0, C(i - 1, j), C(i, j), 0); fLL[f] = i >= 2 ? C(i - 2, j) : ghostOfW[j]; fRR[f] = i + 1 < ni ? C(i + 1, j) : ghostOfE[j]; cellE[C(i - 1, j)] = f; cellW[C(i, j)] = f; f++;
  }
  for (let i = 0; i < ni; i++) for (let j = 1; j < nj; j++) {
    geom(f, i * N1 + j, (i + 1) * N1 + j, 0, 1, C(i, j - 1), C(i, j), 0); fLL[f] = j >= 2 ? C(i, j - 2) : below(i); fRR[f] = j + 1 < nj ? C(i, j + 1) : ghostOfN[i]; cellN[C(i, j - 1)] = f; cellS[C(i, j)] = f; f++;
  }
  // reconstruction weights that are exact for a linear field on a stretched grid (they reduce to the classical κ scheme on a uniform one)
  const mL1 = new Float64Array(nf), mL2 = new Float64Array(nf), mR1 = new Float64Array(nf), mR2 = new Float64Array(nf);
  {
    const xa = new Float64Array(nt), ya = new Float64Array(nt); xa.set(xc); ya.set(yc);
    for (let k = 0; k < nB; k++) { xa[fR[k]] = xc[fL[k]] + fdx[k]; ya[fR[k]] = yc[fL[k]] + fdy[k]; }
    const dd = (a, b) => Math.hypot(xa[a] - xa[b], ya[a] - ya[b]) || 1e-300, weighted = o.weightedMuscl !== false && o.weightedMuscl !== 0;
    for (let k = nB; k < nf; k++) {
      const L = fL[k], Rr = fR[k], a = fLL[k], b = fRR[k], sL = Math.hypot(fxm[k] - xc[L], fym[k] - yc[L]), sR = Math.hypot(fxm[k] - xc[Rr], fym[k] - yc[Rr]), dLR = dd(L, Rr);
      if (weighted) { if (a >= 0) { mL1[k] = (0.5 * (1 - kap) * sL) / dd(a, L); mL2[k] = (0.5 * (1 + kap) * sL) / dLR; } if (b >= 0) { mR1[k] = (0.5 * (1 - kap) * sR) / dd(b, Rr); mR2[k] = (0.5 * (1 + kap) * sR) / dLR; } }
      else { mL1[k] = mR1[k] = 0.25 * (1 - kap); mL2[k] = mR2[k] = 0.25 * (1 + kap); }
    }
  }
  // wall faces and wall distance
  const wallF = []; for (let k = 0; k < nB; k++) if (fT[k] === WALL) wallF.push(k);
  const nWall = wallF.length, dist = new Float64Array(nc).fill(1e30), wx0 = new Float64Array(nWall), wy0 = new Float64Array(nWall), wx1 = new Float64Array(nWall), wy1 = new Float64Array(nWall);
  wallF.forEach((k, m) => { const i = (fL[k] / nj) | 0; wx0[m] = gx[i * N1]; wy0[m] = gy[i * N1]; wx1[m] = gx[(i + 1) * N1]; wy1[m] = gy[(i + 1) * N1]; });
  if (nWall) for (let q = 0; q < nc; q++) {
    let d2 = 1e60; const px = xc[q], py = yc[q];
    for (let m = 0; m < nWall; m++) { const ex = wx1[m] - wx0[m], ey = wy1[m] - wy0[m], ux = px - wx0[m], uy = py - wy0[m]; let tt = (ux * ex + uy * ey) / (ex * ex + ey * ey); tt = tt < 0 ? 0 : tt > 1 ? 1 : tt; const bx = ux - tt * ex, by = uy - tt * ey, dd = bx * bx + by * by; if (dd < d2) d2 = dd; }
    dist[q] = Math.sqrt(d2);
  }
  // wall pressure: extrapolated from the first two cells (the normal pressure gradient of a curved wall matters when the first cell is tall)
  const pwc = new Float64Array(nWall);
  if (nj > 1) wallF.forEach((k, m) => { const L = fL[k]; pwc[m] = ((o.wallPressure ?? 0.6) * dist[L]) / Math.max(dist[L + 1] - dist[L], 1e-300); });
  // implicit lines: wall-normal grid lines; the two halves of a wake line are joined through the cut
  const lineStart = [0], lineCell = [], lineOf = new Int32Array(nc);
  for (let i = 0; i < ni; i++) {
    if (g.bottom[i] === CUT) { const m = ni - 1 - i; if (i > m) continue; for (let j = nj - 1; j >= 0; j--) lineCell.push(C(i, j)); for (let j = 0; j < nj; j++) lineCell.push(C(m, j)); }
    else for (let j = 0; j < nj; j++) lineCell.push(C(i, j));
    lineStart.push(lineCell.length);
  }
  const nLines = lineStart.length - 1, upB = new Int32Array(nc).fill(-1), loB = new Int32Array(nc).fill(-1), xB = new Int32Array(2 * nc).fill(-1), xNb = new Int32Array(2 * nc).fill(-1);
  {
    const facesOf = (q) => [cellS[q], cellN[q], cellW[q], cellE[q]], other = (k, q) => (fL[k] === q ? fR[k] : fL[k]);
    for (let l = 0; l < nLines; l++) for (let p = lineStart[l]; p < lineStart[l + 1]; p++) lineOf[lineCell[p]] = l;
    for (let l = 0; l < nLines; l++) for (let p = lineStart[l]; p < lineStart[l + 1]; p++) {
      const q = lineCell[p], prev = p > lineStart[l] ? lineCell[p - 1] : -1, next = p + 1 < lineStart[l + 1] ? lineCell[p + 1] : -1; let nx = 0;
      for (const k of facesOf(q)) {
        if (k < 0 || fT[k] !== 0) continue; const ob = other(k, q), side = 2 * k + (fL[k] === q ? 0 : 1);
        if (ob === next) upB[p] = side; else if (ob === prev) loB[p] = side; else { xB[2 * p + nx] = side; xNb[2 * p + nx] = ob; nx++; }
      }
    }
  }

  // ---- state -------------------------------------------------------------------------------------
  const P = new Float64Array(nt), U = new Float64Array(nt), V = new Float64Array(nt), NT = new Float64Array(nt), MUT = new Float64Array(nt);
  const gux = new Float64Array(nc), guy = new Float64Array(nc), gvx = new Float64Array(nc), gvy = new Float64Array(nc);
  const R0 = new Float64Array(nc), R1 = new Float64Array(nc), R2 = new Float64Array(nc), RT = new Float64Array(nc), lam = new Float64Array(nc), mdot = new Float64Array(nf);
  const D = new Float64Array(9 * nc), OFF = new Float64Array(18 * nf), DS = new Float64Array(nc), SOFF = new Float64Array(2 * nf);
  const BI = new Float64Array(9 * nc), MM = new Float64Array(9 * nc), dQ = new Float64Array(3 * nc), rr = new Float64Array(3 * nc), tb = new Float64Array(9), sB = new Float64Array(nc), sM = new Float64Array(nc), dT = new Float64Array(nc), rs = new Float64Array(nc);
  const farP = new Float64Array(nB), farU = new Float64Array(nB), farV = new Float64Array(nB), utau = new Float64Array(nWall), muW = new Float64Array(nWall).fill(nu);
  const history = { res: [], cl: [], cd: [] };
  let iter = 0, cfl = o.cflStart ?? 10, res0 = 0, resNow = 0, resT = 0, diverged = false;
  const cflMax = o.cflMax ?? 1000, cflGrow = o.cflGrow ?? 1.25, cflMin = o.cflStart ?? 10, dMax = o.maxChange ?? 0.4, serDrop = o.serDrop ?? 3e-3, rejectAt = o.rejectAt ?? 2.5, sP = new Float64Array(nt), sU = new Float64Array(nt), sV = new Float64Array(nt), sN = new Float64Array(nt); const sw2 = new Float64Array(2); let resPrev = 0, rejected = 0, rejections = 0, nSweeps = 0;

  function setFar() {
    const k0 = gammaFar / (2 * Math.PI);
    for (let k = 0; k < nB; k++) {
      if (fT[k] !== FAR) continue;
      if (farFn) { const s = farFn(fxm[k], fym[k]); farP[k] = s[0]; farU[k] = s[1]; farV[k] = s[2]; continue; }
      const rx = fxm[k] - xRef, ry = fym[k], r2 = rx * rx + ry * ry || 1, u = ca + (k0 * ry) / r2, v = sa - (k0 * rx) / r2;
      farU[k] = u; farV[k] = v; farP[k] = 0.5 * (1 - u * u - v * v);
    }
  }
  function initialise() {
    if (o.clInit && vortex) gammaFar = 0.5 * o.clInit;
    setFar();
    const dl = o.dampInit ?? (turb ? 0.002 : 0.02);
    for (let q = 0; q < nc; q++) {
      let u = ca, v = sa, p = 0;
      if (o.init) { const s = o.init(xc[q], yc[q]); p = s[0]; u = s[1]; v = s[2]; }
      else {
        if (o.clInit) { const k0 = o.clInit / (4 * Math.PI), rx = xc[q] - xRef, ry = yc[q], r2 = rx * rx + ry * ry + 0.25; u += (k0 * ry) / r2; v -= (k0 * rx) / r2; p = 0.5 * (1 - u * u - v * v); }
        if (nWall && !o.noDamp) { const s = 1 - Math.exp(-dist[q] / dl); u *= s; v *= s; }
      }
      P[q] = p; U[q] = u; V[q] = v; NT[q] = turb ? ntInf + (nWall && !slip ? ntInit * dist[q] * Math.exp(-dist[q] / 0.02) : 0) : 0; // a turbulent layer from the start avoids a laminar-separation transient
    }
    ghosts();
  }
  function ghosts() {
    for (let k = 0; k < nB; k++) {
      const L = fL[k], G = fR[k], ty = fT[k];
      if (ty === WALL && !slip) { P[G] = P[L]; U[G] = -U[L]; V[G] = -V[L]; NT[G] = -NT[L]; MUT[G] = -MUT[L]; }
      else if (ty !== FAR) { const s2 = fnx[k] * fnx[k] + fny[k] * fny[k], un = (2 * (U[L] * fnx[k] + V[L] * fny[k])) / s2; P[G] = P[L]; U[G] = U[L] - un * fnx[k]; V[G] = V[L] - un * fny[k]; NT[G] = NT[L]; MUT[G] = MUT[L]; }
      else { const out = U[L] * fnx[k] + V[L] * fny[k] > 0; P[G] = farP[k]; if (out && !farFn) { U[G] = U[L]; V[G] = V[L]; } else { U[G] = farU[k]; V[G] = farV[k]; } // outflow: pressure fixed, velocity leaves freely
        NT[G] = out ? NT[L] : turb ? ntInf : 0; MUT[G] = out ? MUT[L] : MUT[G]; }
    }
  }
  const fv1 = (chi) => { const c3 = chi * chi * chi; return c3 / (c3 + CV13); };
  function eddy() {
    if (!turb) return;
    for (let q = 0; q < nc; q++) { const nt2 = NT[q]; MUT[q] = nt2 > 0 ? nt2 * fv1(nt2 / nu) : 0; }
    const mi = ntInf * fv1(ntInf / nu); for (let k = 0; k < nB; k++) if (fT[k] === FAR) MUT[fR[k]] = mi;
  }
  /** Friction velocity from the log law (with the viscous sublayer below y⁺ = 11.06) for speed u at distance d. */
  function logLaw(u, d) {
    const Rey = (u * d) / nu; if (!(Rey > 0)) return 0;
    if (Rey <= 122.3) return Math.sqrt((nu * u) / d);
    let ut = u / (Math.log(Rey) / KAPPA + B_LOG) * 1.5;
    for (let it = 0; it < 8; it++) { const lg = Math.log((d * ut) / nu) / KAPPA + B_LOG, un = ut - (ut * lg - u) / (lg + 1 / KAPPA); if (!(un > 0)) { ut *= 0.5; continue; } if (Math.abs(un - ut) < 1e-12 * ut) { ut = un; break; } ut = un; }
    return ut;
  }
  function wallModel() {
    for (let m = 0; m < nWall; m++) {
      const k = wallF[m], L = fL[k], s = Math.sqrt(fnx[k] * fnx[k] + fny[k] * fny[k]), tx = -fny[k] / s, ty = fnx[k] / s, ut = U[L] * tx + V[L] * ty, d = 0.5 * Math.hypot(fdx[k], fdy[k]);
      if (slip) { utau[m] = 0; muW[m] = 0; }
      else if (wf) { const us = logLaw(Math.abs(ut), d); utau[m] = us; muW[m] = Math.abs(ut) > 1e-12 ? Math.max(nu, (us * us * d) / Math.abs(ut)) : nu; if (us * d > 11.06 * nu) NT[L] = KAPPA * us * d; }
      else { utau[m] = Math.sqrt((nu * Math.abs(ut)) / d); muW[m] = nu; }
    }
  }

  // ---- residual ----------------------------------------------------------------------------------
  function gradients() {
    gux.fill(0); guy.fill(0); gvx.fill(0); gvy.fill(0);
    for (let k = 0; k < nf; k++) {
      const L = fL[k], Rr = fR[k], w = fw[k], uf = w * U[L] + (1 - w) * U[Rr], vf = w * V[L] + (1 - w) * V[Rr], ax = fnx[k], ay = fny[k];
      gux[L] += uf * ax; guy[L] += uf * ay; gvx[L] += vf * ax; gvy[L] += vf * ay;
      if (Rr < nc) { gux[Rr] -= uf * ax; guy[Rr] -= uf * ay; gvx[Rr] -= vf * ax; gvy[Rr] -= vf * ay; }
    }
    for (let q = 0; q < nc; q++) { const iv = 1 / vol[q]; gux[q] *= iv; guy[q] *= iv; gvx[q] *= iv; gvy[q] *= iv; }
  }
  function residual(jac) {
    const beta = betaNow[0];
    ghosts(); gradients();
    R0.fill(0); R1.fill(0); R2.fill(0); lam.fill(0);
    if (jac) D.fill(0);
    const hi = kap > -1.5;
    let mw = 0;
    for (let k = 0; k < nf; k++) {
      const L = fL[k], Rr = fR[k], ty = fT[k], ax = fnx[k], ay = fny[k], s2 = ax * ax + ay * ay, inner = Rr < nc;
      let f0, f1, f2, kv = 0;
      if (ty === WALL || ty === SYM) {
        const pw = ty === WALL ? P[L] - pwc[mw] * (P[L + 1] - P[L]) : P[L]; f0 = 0; f1 = pw * ax; f2 = pw * ay; mdot[k] = 0;
        lam[L] += Math.sqrt(beta * s2);
        if (jac) { const o9 = 9 * L; D[o9 + 3] += ax; D[o9 + 6] += ay; }
        if (ty === SYM || slip) { R1[L] += f1; R2[L] += f2; if (ty === WALL) mw++; continue; }
      } else {
        let pl = P[L], ul = U[L], vl = V[L], pr = P[Rr], ur = U[Rr], vr = V[Rr];
        if (hi && ty === 0) {
          const a = fLL[k], b = fRR[k];
          if (a >= 0) { const k1 = mL1[k], k2 = mL2[k]; pl += k1 * (pl - P[a]) + k2 * (pr - pl); ul += k1 * (ul - U[a]) + k2 * (ur - ul); vl += k1 * (vl - V[a]) + k2 * (vr - vl); }
          if (b >= 0) { const k1 = mR1[k], k2 = mR2[k], p0 = P[Rr], u0 = U[Rr], v0 = V[Rr], pL0 = P[L], uL0 = U[L], vL0 = V[L]; pr = p0 - k1 * (P[b] - p0) - k2 * (p0 - pL0); ur = u0 - k1 * (U[b] - u0) - k2 * (u0 - uL0); vr = v0 - k1 * (V[b] - v0) - k2 * (v0 - vL0); }
        }
        const thl = ul * ax + vl * ay, thr = ur * ax + vr * ay, uh = 0.5 * (ul + ur), vh = 0.5 * (vl + vr), th = 0.5 * (thl + thr), c = Math.sqrt(th * th + beta * s2), at = Math.abs(th);
        const dp = pr - pl, du = ur - ul, dv = vr - vl, t1 = du * ax + dv * ay, w0 = -th * dp + beta * t1, w1 = ax * dp + uh * t1, w2 = ay * dp + vh * t1, t2 = w1 * ax + w2 * ay;
        const e1 = th / c, e2 = (c - at) / (c * c);
        f0 = 0.5 * (beta * (thl + thr) - (at * dp + e1 * w0 + e2 * (-th * w0 + beta * t2)));
        f1 = 0.5 * (ul * thl + ur * thr + (pl + pr) * ax - (at * du + e1 * w1 + e2 * (ax * w0 + uh * t2)));
        f2 = 0.5 * (vl * thl + vr * thr + (pl + pr) * ay - (at * dv + e1 * w2 + e2 * (ay * w0 + vh * t2)));
        mdot[k] = f0 / beta;
        const sr = at + c; lam[L] += sr; if (inner) lam[Rr] += sr;
        if (jac) {
          // first-order Jacobians A± = (Â ± |Â|)/2 at the cell states, |Â| = |θ| I + (θ/c) B + ((c − |θ|)/c²) B², B = Â − θ I
          const u1 = 0.5 * (U[L] + U[Rr]), v1 = 0.5 * (V[L] + V[Rr]), t0 = u1 * ax + v1 * ay, cc = Math.sqrt(t0 * t0 + beta * s2), a0 = Math.abs(t0), g1 = t0 / cc, g2 = (cc - a0) / (cc * cc);
          const b0 = -t0, b1 = beta * ax, b2 = beta * ay, b3 = ax, b4 = u1 * ax, b5 = u1 * ay, b6 = ay, b7 = v1 * ax, b8 = v1 * ay;
          const q0 = b0 * b0 + b1 * b3 + b2 * b6, q1 = b0 * b1 + b1 * b4 + b2 * b7, q2 = b0 * b2 + b1 * b5 + b2 * b8, q3 = b3 * b0 + b4 * b3 + b5 * b6, q4 = b3 * b1 + b4 * b4 + b5 * b7, q5 = b3 * b2 + b4 * b5 + b5 * b8, q6 = b6 * b0 + b7 * b3 + b8 * b6, q7 = b6 * b1 + b7 * b4 + b8 * b7, q8 = b6 * b2 + b7 * b5 + b8 * b8;
          const m0 = a0 + g1 * b0 + g2 * q0, m1 = g1 * b1 + g2 * q1, m2 = g1 * b2 + g2 * q2, m3 = g1 * b3 + g2 * q3, m4 = a0 + g1 * b4 + g2 * q4, m5 = g1 * b5 + g2 * q5, m6 = g1 * b6 + g2 * q6, m7 = g1 * b7 + g2 * q7, m8 = a0 + g1 * b8 + g2 * q8;
          const h0 = b0 + t0, h4 = b4 + t0, h8 = b8 + t0, oL = 9 * L;
          D[oL] += 0.5 * (h0 + m0); D[oL + 1] += 0.5 * (b1 + m1); D[oL + 2] += 0.5 * (b2 + m2); D[oL + 3] += 0.5 * (b3 + m3); D[oL + 4] += 0.5 * (h4 + m4); D[oL + 5] += 0.5 * (b5 + m5); D[oL + 6] += 0.5 * (b6 + m6); D[oL + 7] += 0.5 * (b7 + m7); D[oL + 8] += 0.5 * (h8 + m8);
          if (inner) {
            const oR = 9 * Rr, oa = 18 * k, ob = oa + 9;
            D[oR] -= 0.5 * (h0 - m0); D[oR + 1] -= 0.5 * (b1 - m1); D[oR + 2] -= 0.5 * (b2 - m2); D[oR + 3] -= 0.5 * (b3 - m3); D[oR + 4] -= 0.5 * (h4 - m4); D[oR + 5] -= 0.5 * (b5 - m5); D[oR + 6] -= 0.5 * (b6 - m6); D[oR + 7] -= 0.5 * (b7 - m7); D[oR + 8] -= 0.5 * (h8 - m8);
            OFF[oa] = 0.5 * (h0 - m0); OFF[oa + 1] = 0.5 * (b1 - m1); OFF[oa + 2] = 0.5 * (b2 - m2); OFF[oa + 3] = 0.5 * (b3 - m3); OFF[oa + 4] = 0.5 * (h4 - m4); OFF[oa + 5] = 0.5 * (b5 - m5); OFF[oa + 6] = 0.5 * (b6 - m6); OFF[oa + 7] = 0.5 * (b7 - m7); OFF[oa + 8] = 0.5 * (h8 - m8);
            OFF[ob] = -0.5 * (h0 + m0); OFF[ob + 1] = -0.5 * (b1 + m1); OFF[ob + 2] = -0.5 * (b2 + m2); OFF[ob + 3] = -0.5 * (b3 + m3); OFF[ob + 4] = -0.5 * (h4 + m4); OFF[ob + 5] = -0.5 * (b5 + m5); OFF[ob + 6] = -0.5 * (b6 + m6); OFF[ob + 7] = -0.5 * (b7 + m7); OFF[ob + 8] = -0.5 * (h8 + m8);
          }
        }
      }
      // viscous stress: averaged cell gradients corrected along the line of centres
      const dx = fdx[k], dy = fdy[k], id2 = 1 / (dx * dx + dy * dy), w = fw[k];
      let ux, uy, vx, vy, mu;
      if (inner) { ux = w * gux[L] + (1 - w) * gux[Rr]; uy = w * guy[L] + (1 - w) * guy[Rr]; vx = w * gvx[L] + (1 - w) * gvx[Rr]; vy = w * gvy[L] + (1 - w) * gvy[Rr]; mu = nu + w * MUT[L] + (1 - w) * MUT[Rr]; }
      else { ux = gux[L]; uy = guy[L]; vx = gvx[L]; vy = gvy[L]; mu = ty === FAR ? nu + MUT[L] : muW[mw++]; }
      const cu = (U[Rr] - U[L] - ux * dx - uy * dy) * id2, cv = (V[Rr] - V[L] - vx * dx - vy * dy) * id2;
      ux += cu * dx; uy += cu * dy; vx += cv * dx; vy += cv * dy;
      const sxy = uy + vx; f1 -= mu * (2 * ux * ax + sxy * ay); f2 -= mu * (sxy * ax + 2 * vy * ay);
      kv = mu * fg[k];
      R0[L] += f0; R1[L] += f1; R2[L] += f2; lam[L] += 2 * kv;
      if (inner) { R0[Rr] -= f0; R1[Rr] -= f1; R2[Rr] -= f2; lam[Rr] += 2 * kv; }
      if (jac) {
        if (inner) { const oL = 9 * L, oR = 9 * Rr, oa = 18 * k, ob = oa + 9; D[oL + 4] += kv; D[oL + 8] += kv; D[oR + 4] += kv; D[oR + 8] += kv; OFF[oa + 4] -= kv; OFF[oa + 8] -= kv; OFF[ob + 4] -= kv; OFF[ob + 8] -= kv; }
        else { const oL = 9 * L, k2v = ty === FAR ? kv : 2 * kv; D[oL + 4] += k2v; D[oL + 8] += k2v; }
      }
    }
    let s = 0; for (let q = 0; q < nc; q++) { const a = R0[q] / (beta * vol[q]), b = R1[q] / vol[q], c = R2[q] / vol[q]; s += (a * a + b * b + c * c) * vol[q]; }
    return Math.sqrt(s);
  }
  function residualSA(jac) {
    RT.fill(0); if (jac) DS.fill(0);
    let mw = 0;
    for (let k = 0; k < nf; k++) {
      const L = fL[k], Rr = fR[k], ty = fT[k], inner = Rr < nc;
      if (ty === SYM) continue;
      const nl = NT[L], nr = NT[Rr], md = mdot[k];
      if (md > 0) { if (inner) { RT[Rr] -= md * (nl - nr); if (jac) { DS[Rr] += md; SOFF[2 * k + 1] = -md; SOFF[2 * k] = 0; } } }
      else { RT[L] += md * (nr - nl); if (jac) { DS[L] -= md; if (inner) { SOFF[2 * k] = md; SOFF[2 * k + 1] = 0; } } }
      if (ty === FAR) continue;
      const wall = ty === WALL; if (wall && (wf || slip)) { mw++; continue; }
      const xl = nl / nu, xr = nr / nu, el = nl >= 0 ? nu + nl : nu + nl * ((CN1 + xl * xl * xl) / (CN1 - xl * xl * xl)), er = wall ? 2 * nu - el : nr >= 0 ? nu + nr : nu + nr * ((CN1 + xr * xr * xr) / (CN1 - xr * xr * xr));
      const w = fw[k], ef = w * el + (1 - w) * er, gg = fg[k] / SIG, cL = ((1 + CB2) * ef - CB2 * el) * gg, cR = ((1 + CB2) * ef - CB2 * er) * gg, dn = nr - nl;
      RT[L] -= cL * dn; if (inner) RT[Rr] += cR * dn;
      if (jac) { const pL = cL > 0 ? cL : 0, pR = cR > 0 ? cR : 0; if (inner) { DS[L] += pL; DS[Rr] += pR; SOFF[2 * k] -= pL; SOFF[2 * k + 1] -= pR; } else DS[L] += 2 * pL; }
    }
    let s = 0;
    for (let q = 0; q < nc; q++) {
      const nt2 = NT[q], d = dist[q], om = Math.abs(gvx[q] - guy[q]), kd2 = KAPPA * KAPPA * d * d; let src, dsrc;
      if (nt2 >= 0) {
        const chi = nt2 / nu, c3 = chi * chi * chi, f1 = c3 / (c3 + CV13), f2 = 1 - chi / (1 + chi * f1), sb = (nt2 * f2) / kd2;
        const st = sb >= -CV2 * om ? om + sb : om + (om * (CV2 * CV2 * om + CV3 * sb)) / ((CV3 - 2 * CV2) * om - sb);
        let r = st > 1e-30 ? nt2 / (st * kd2) : 10; if (r > 10 || r < 0) r = 10;
        const r2 = r * r, gq = r + CW2 * (r2 * r2 * r2 - r), g2 = gq * gq, fwq = gq * Math.cbrt(Math.sqrt((1 + CW36) / (g2 * g2 * g2 + CW36))), ds = (CW1 * fwq * nt2) / (d * d);
        src = CB1 * st * nt2 - ds * nt2; dsrc = 2 * ds;
      } else { src = CB1 * om * nt2 + (CW1 * nt2 * nt2) / (d * d); dsrc = CB1 * om; }
      RT[q] -= src * vol[q]; if (jac) DS[q] += dsrc * vol[q];
      const a = RT[q] / vol[q]; s += a * a * vol[q];
    }
    if (wf) for (let m = 0; m < nWall; m++) { const L = fL[wallF[m]]; if (utau[m] * 0.5 * Math.hypot(fdx[wallF[m]], fdy[wallF[m]]) > 11.06 * nu) { RT[L] = 0; DS[L] = 1e30; } }
    return Math.sqrt(s) * Re;
  }

  // ---- implicit line solves ----------------------------------------------------------------------
  function factor() {
    const ic = 1 / cfl;
    for (let l = 0; l < nLines; l++) {
      const a = lineStart[l], b = lineStart[l + 1];
      for (let p = a; p < b; p++) {
        const q = lineCell[p], o9 = 9 * q, t = lam[q] * ic, pp = 9 * p;
        for (let m = 0; m < 9; m++) tb[m] = D[o9 + m];
        tb[0] += t; tb[4] += t; tb[8] += t;
        if (p > a) { // M = A · inv(B'_{p-1});  B' = B − M · C_{p-1}
          const A = 9 * loB[p], Bp = pp - 9, Cu = 9 * upB[p - 1];
          for (let r = 0; r < 3; r++) {
            const x0 = OFF[A + 3 * r], x1 = OFF[A + 3 * r + 1], x2 = OFF[A + 3 * r + 2];
            const y0 = x0 * BI[Bp] + x1 * BI[Bp + 3] + x2 * BI[Bp + 6], y1 = x0 * BI[Bp + 1] + x1 * BI[Bp + 4] + x2 * BI[Bp + 7], y2 = x0 * BI[Bp + 2] + x1 * BI[Bp + 5] + x2 * BI[Bp + 8];
            MM[pp + 3 * r] = y0; MM[pp + 3 * r + 1] = y1; MM[pp + 3 * r + 2] = y2;
            tb[3 * r] -= y0 * OFF[Cu] + y1 * OFF[Cu + 3] + y2 * OFF[Cu + 6]; tb[3 * r + 1] -= y0 * OFF[Cu + 1] + y1 * OFF[Cu + 4] + y2 * OFF[Cu + 7]; tb[3 * r + 2] -= y0 * OFF[Cu + 2] + y1 * OFF[Cu + 5] + y2 * OFF[Cu + 8];
          }
        }
        inv3(tb, 0, BI, pp);
      }
    }
  }
  function solveLine(l) {
    const a = lineStart[l], b = lineStart[l + 1];
    for (let p = a; p < b; p++) {
      const q = lineCell[p]; let r0 = -R0[q], r1 = -R1[q], r2 = -R2[q];
      for (let e = 0; e < 2; e++) { const sd = xB[2 * p + e]; if (sd < 0) continue; const o9 = 9 * sd, n3 = 3 * xNb[2 * p + e], d0 = dQ[n3], d1 = dQ[n3 + 1], d2 = dQ[n3 + 2]; r0 -= OFF[o9] * d0 + OFF[o9 + 1] * d1 + OFF[o9 + 2] * d2; r1 -= OFF[o9 + 3] * d0 + OFF[o9 + 4] * d1 + OFF[o9 + 5] * d2; r2 -= OFF[o9 + 6] * d0 + OFF[o9 + 7] * d1 + OFF[o9 + 8] * d2; }
      if (p > a) { const pp = 9 * p, s0 = rr[3 * p - 3], s1 = rr[3 * p - 2], s2 = rr[3 * p - 1]; r0 -= MM[pp] * s0 + MM[pp + 1] * s1 + MM[pp + 2] * s2; r1 -= MM[pp + 3] * s0 + MM[pp + 4] * s1 + MM[pp + 5] * s2; r2 -= MM[pp + 6] * s0 + MM[pp + 7] * s1 + MM[pp + 8] * s2; }
      rr[3 * p] = r0; rr[3 * p + 1] = r1; rr[3 * p + 2] = r2;
    }
    let x0 = 0, x1 = 0, x2 = 0, ch = 0, nm = 0;
    for (let p = b - 1; p >= a; p--) {
      let r0 = rr[3 * p], r1 = rr[3 * p + 1], r2 = rr[3 * p + 2];
      if (p < b - 1) { const o9 = 9 * upB[p]; r0 -= OFF[o9] * x0 + OFF[o9 + 1] * x1 + OFF[o9 + 2] * x2; r1 -= OFF[o9 + 3] * x0 + OFF[o9 + 4] * x1 + OFF[o9 + 5] * x2; r2 -= OFF[o9 + 6] * x0 + OFF[o9 + 7] * x1 + OFF[o9 + 8] * x2; }
      const pp = 9 * p; x0 = BI[pp] * r0 + BI[pp + 1] * r1 + BI[pp + 2] * r2; x1 = BI[pp + 3] * r0 + BI[pp + 4] * r1 + BI[pp + 5] * r2; x2 = BI[pp + 6] * r0 + BI[pp + 7] * r1 + BI[pp + 8] * r2;
      const q3 = 3 * lineCell[p], e0 = x0 - dQ[q3], e1 = x1 - dQ[q3 + 1], e2 = x2 - dQ[q3 + 2]; ch += e0 * e0 + e1 * e1 + e2 * e2; nm += x0 * x0 + x1 * x1 + x2 * x2;
      dQ[q3] = x0; dQ[q3 + 1] = x1; dQ[q3 + 2] = x2;
    }
    sw2[0] += ch; sw2[1] += nm; // typed-array accumulators: a closure variable would allocate on every write
  }
  function factorSA() {
    const ic = 1 / cfl;
    for (let l = 0; l < nLines; l++) {
      const a = lineStart[l], b = lineStart[l + 1];
      for (let p = a; p < b; p++) {
        const q = lineCell[p]; let d = DS[q] + lam[q] * ic;
        if (p > a) { sM[p] = SOFF[loB[p]] * sB[p - 1]; d -= sM[p] * SOFF[upB[p - 1]]; }
        sB[p] = 1 / d;
      }
    }
  }
  function solveLineSA(l) {
    const a = lineStart[l], b = lineStart[l + 1];
    for (let p = a; p < b; p++) {
      const q = lineCell[p]; let r = -RT[q];
      for (let e = 0; e < 2; e++) { const sd = xB[2 * p + e]; if (sd >= 0) r -= SOFF[sd] * dT[xNb[2 * p + e]]; }
      if (p > a) r -= sM[p] * rs[p - 1];
      rs[p] = r;
    }
    let xn = 0;
    for (let p = b - 1; p >= a; p--) { let r = rs[p]; if (p < b - 1) r -= SOFF[upB[p]] * xn; xn = sB[p] * r; dT[lineCell[p]] = xn; }
  }

  // ---- iteration ---------------------------------------------------------------------------------
  const relaxF = o.relaxFlow ?? 1, relaxT = o.relaxTurb ?? 0.35, sweeps = o.sweeps ?? 2, sweepsMax = Math.max(sweeps, o.sweepsMax ?? 12), sweepTol = o.sweepTol ?? 0.1, vortexEvery = o.vortexEvery ?? 5;
  function step() {
    if (!betaFix && iter % 5 === 0) { // artificial-compressibility parameter: at least the largest speed squared in the field (suction peak)
      let u2 = 1; for (let q = 0; q < nc; q++) { const w = U[q] * U[q] + V[q] * V[q]; if (w > u2) u2 = w; }
      const bt = Math.min(Math.max(1, betaK * u2), 25); betaNow[0] = iter === 0 ? bt : Math.max(bt, 0.9 * betaNow[0]);
    }
    eddy(); if (nWall) wallModel();
    resNow = residual(true);
    if (iter > 0 && !(resNow <= rejectAt * resPrev) && cfl > 1.01 * cflMin && rejected < 4) {
      // the last update made the residual jump: take it back and repeat it with a smaller pseudo-time step
      P.set(sP); U.set(sU); V.set(sV); NT.set(sN); cfl = Math.max(cflMin, 0.25 * cfl); rejected++; rejections++;
      eddy(); if (nWall) wallModel(); resNow = residual(true);
    } else rejected = 0;
    if (!Number.isFinite(resNow)) { diverged = true; return; }
    if (iter === 0) res0 = resNow || 1e-300;
    sP.set(P); sU.set(U); sV.set(V); sN.set(NT);
    factor(); dQ.fill(0);
    // symmetric line Gauss–Seidel until the update stops changing; a system too stiff for the sweeps lowers the CFL number
    let stiff = true;
    for (let sw = 0; sw < sweepsMax; sw++) {
      sw2[0] = 0; sw2[1] = 0;
      if (sw % 2 === 0) for (let l = 0; l < nLines; l++) solveLine(l); else for (let l = nLines - 1; l >= 0; l--) solveLine(l);
      nSweeps++;
      if (sw + 1 >= sweeps && sw2[0] <= sweepTol * sweepTol * sw2[1]) { stiff = false; break; }
    }
    let mx = 0; for (let q = 0; q < nc; q++) { const a = Math.abs(dQ[3 * q + 1]), b = Math.abs(dQ[3 * q + 2]), c = Math.abs(dQ[3 * q]); if (a > mx) mx = a; if (b > mx) mx = b; if (c > mx) mx = c; }
    const om = mx > dMax ? dMax / mx : 1, of = om * relaxF;
    for (let q = 0; q < nc; q++) { P[q] += of * dQ[3 * q]; U[q] += of * dQ[3 * q + 1]; V[q] += of * dQ[3 * q + 2]; }
    if (turb) {
      ghosts(); gradients();
      // mass fluxes of the updated field are close enough to the stored ones for the first-order upwind turbulence convection
      resT = residualSA(true); factorSA(); dT.fill(0);
      for (let sw = 0; sw < 4; sw++) { if (sw % 2 === 0) for (let l = 0; l < nLines; l++) solveLineSA(l); else for (let l = nLines - 1; l >= 0; l--) solveLineSA(l); }
      for (let q = 0; q < nc; q++) { const v = NT[q], w = v + relaxT * dT[q]; NT[q] = v > 0 && w < 0.1 * v ? 0.1 * v : w; }
    }
    iter++;
    // switched evolution relaxation while the residual is still high: the CFL number follows the residual, so a rough
    // transient is marched with small steps; once the residual has dropped the CFL number simply grows to its maximum
    const ratio = resPrev > 0 && resNow > 0 ? (resPrev / resNow) * 1.1 : cflGrow, early = resNow > serDrop * res0; resPrev = resNow;
    cfl = om < 1 || stiff ? Math.max(cflMin, 0.5 * cfl) : Math.max(cflMin, Math.min(cflMax, cfl * (early ? Math.max(0.5, Math.min(cflGrow, ratio)) : cflGrow)));
    const F = forces(); history.res.push(resNow / res0); history.cl.push(F.cl); history.cd.push(F.cd);
    if (vortex && !farFn && iter % vortexEvery === 0) { gammaFar = 0.5 * F.cl; setFar(); }
  }
  function iterate(n = 1, tol = 0) {
    for (let k = 0; k < n && !diverged; k++) { step(); if (tol > 0 && resNow < tol * res0) break; }
    return Math.log10(res0 / Math.max(resNow, 1e-300));
  }
  function forces() {
    let fxp = 0, fyp = 0, fxv = 0, fyv = 0, mz = 0;
    for (let m = 0; m < nWall; m++) {
      const k = wallF[m], L = fL[k], pw = P[L] - pwc[m] * (P[L + 1] - P[L]), px = pw * fnx[k], py = pw * fny[k], kv = 2 * muW[m] * fg[k], s2 = fnx[k] * fnx[k] + fny[k] * fny[k], un = (U[L] * fnx[k] + V[L] * fny[k]) / s2, vx = kv * (U[L] - un * fnx[k]), vy = kv * (V[L] - un * fny[k]);
      fxp += px; fyp += py; fxv += vx; fyv += vy; mz += (fxm[k] - xRef) * (py + vy) - fym[k] * (px + vx);
    }
    const cdp = 2 * (fxp * ca + fyp * sa), cdf = 2 * (fxv * ca + fyv * sa);
    return { cl: 2 * ((fyp + fyv) * ca - (fxp + fxv) * sa), cd: cdp + cdf, cdp, cdf, cm: -2 * mz };
  }
  function surface() {
    const out = { x: [], y: [], cp: [], cf: [], yplus: [], utau: [], side: [] };
    for (let m = 0; m < nWall; m++) {
      const k = wallF[m], L = fL[k], s = Math.sqrt(fnx[k] * fnx[k] + fny[k] * fny[k]), tx = -fny[k] / s, ty = fnx[k] / s, ut = U[L] * tx + V[L] * ty, d = 0.5 * Math.hypot(fdx[k], fdy[k]), tw = (muW[m] * ut) / d;
      // the wall tangent (tx, ty) points in the direction of increasing i; on the lower surface of a C-grid that is upstream
      const dir = g.kind === 'c' ? (tx >= 0 ? 1 : -1) : 1, upper = g.kind === 'c' ? m >= nWall / 2 : true;
      out.x.push(fxm[k]); out.y.push(fym[k]); out.cp.push(2 * (P[L] - pwc[m] * (P[L + 1] - P[L]))); out.cf.push(2 * tw * dir); out.yplus.push((Math.sqrt(Math.abs(tw)) * d) / nu); out.utau.push(Math.sqrt(Math.abs(tw))); out.side.push(upper ? 1 : -1);
    }
    return out;
  }
  function profile(xs, upper = true) {
    let best = -1, bd = 1e30;
    for (let m = 0; m < nWall; m++) { if (g.kind === 'c' && (m >= nWall / 2) !== upper) continue; const d = Math.abs(fxm[wallF[m]] - xs); if (d < bd) { bd = d; best = m; } }
    if (best < 0) return null;
    const k = wallF[best], i = (fL[k] / nj) | 0, s = Math.sqrt(fnx[k] * fnx[k] + fny[k] * fny[k]), tx = -fny[k] / s, ty = fnx[k] / s, d1 = 0.5 * Math.hypot(fdx[k], fdy[k]), q0 = C(i, 0), ut0 = U[q0] * tx + V[q0] * ty, us = Math.sqrt((muW[best] * Math.abs(ut0)) / d1), sg = ut0 < 0 ? -1 : 1;
    const out = { x: fxm[k], utau: us, y: [], u: [], yplus: [], uplus: [], nut: [] };
    for (let j = 0; j < nj; j++) { const q = C(i, j), d = dist[q], u = sg * (U[q] * tx + V[q] * ty); out.y.push(d); out.u.push(u); out.yplus.push((d * us) / nu); out.uplus.push(us > 0 ? u / us : 0); out.nut.push(MUT[q] / nu); }
    return out;
  }
  function field(name) {
    const out = new Float64Array(nc);
    for (let q = 0; q < nc; q++) out[q] = name === 'p' ? P[q] : name === 'u' ? U[q] : name === 'v' ? V[q] : name === 'nut' ? MUT[q] / nu : name === 'nutilde' ? NT[q] / nu : name === 'dist' ? dist[q] : Math.hypot(U[q], V[q]);
    return out;
  }
  initialise();
  return {
    grid: g, xc, yc, vol, dist, history, iterate, step, forces, surface, profile, field, wallSpacing: nWall ? 0.5 * Math.hypot(fdx[wallF[0]], fdy[wallF[0]]) : 0,
    residualOnly: () => { eddy(); if (nWall) wallModel(); return residual(false); },
    setAlpha(a) { alpha = a; ca = Math.cos(a); sa = Math.sin(a); setFar(); cfl = Math.max(cflMin, Math.min(cfl, o.cflRestart ?? 1000)); resPrev = 0; res0 = 0; iter = 0; },
    get iterations() { return history.res.length; }, get residualDrop() { return Math.log10(res0 / Math.max(resNow, 1e-300)); }, get residual() { return resNow; }, get residualTurb() { return resT; },
    get cfl() { return cfl; }, get sweepCount() { return nSweeps; }, get beta() { return betaNow[0]; }, get rejections() { return rejections; }, get diverged() { return diverged; }, get alpha() { return alpha; }, state: { P, U, V, NT, MUT },
  };
}
