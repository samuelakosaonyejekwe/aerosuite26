// Two-dimensional incompressible Navier–Stokes kernel with an immersed boundary, for fluid–structure coupling.
// Uniform staggered (MAC) grid, typed arrays. Fractional-step (projection) method:
//   1. nsPredict  — explicit advection and diffusion: fourth-order central advection and viscous terms by second-order
//                   Adams–Bashforth, plus the fourth-difference upwind dissipation (third-order upwind / Kawamura–Kuwahara
//                   family) by forward Euler;
//   2. nsForce    — direct-forcing immersed boundary: the velocity inside a smoothed body mask is driven to the rigid-body
//                   velocity; the momentum exchanged is the hydrodynamic load (may be applied several times per step:
//                   force–project–force gives the load including the pressure reaction of the current step);
//   3. nsProject  — exact pressure projection by a discrete cosine transform across the stream and a tridiagonal solve along it.
// Boundary conditions: uniform inflow on the left, convective outflow with zero pressure on the right, free-slip walls
// at the bottom and top. Flow is in +x. All quantities are in consistent (SI or non-dimensional) units.

const G = 2; // ghost layers

/**
 * Create a flow state. o: { nx, ny (even; a power of two uses the fast transform), h, nu, U, rho = 1, upwind = 1 }.
 * upwind scales the fourth-difference dissipation of the advection scheme (1 = third-order upwind, 3 = Kawamura–Kuwahara).
 */
export function nsCreate(o) {
  const nx = Math.round(o.nx), ny = Math.round(o.ny), su = nx + 1 + 2 * G, sv = nx + 2 * G, ru = ny + 2 * G, rv = ny + 1 + 2 * G;
  const s = {
    nx, ny, h: o.h, nu: o.nu, U: o.U, rho: o.rho ?? 1, kk: o.upwind ?? 1, su, sv, t: 0, nStep: 0,
    u: new Float64Array(su * ru).fill(o.U), v: new Float64Array(sv * rv), p: new Float64Array(nx * ny), // p: pressure summed over the projections of the current step
    hu: new Float64Array(su * ru), hv: new Float64Array(sv * rv), hu0: new Float64Array(su * ru), hv0: new Float64Array(sv * rv),
    bu: new Float64Array(su * ru), bv: new Float64Array(sv * rv),
    cs: null, tm: new Float64Array(nx * ny), w1: new Float64Array(nx * ny), w2: new Float64Array(nx * ny),
  };
  if ((ny & (ny - 1)) === 0 && ny >= 8) {
    // power-of-two heights: fast cosine transform through a complex FFT of the same length (Makhoul's reordering)
    const f = { rev: new Uint16Array(ny), wr: new Float64Array(ny / 2), wi: new Float64Array(ny / 2), cr: new Float64Array(ny), sr: new Float64Array(ny), re: new Float64Array(ny), im: new Float64Array(ny) };
    for (let i = 1, j = 0; i < ny; i++) { let bit = ny >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; f.rev[i] = j; }
    for (let k = 0; k < ny / 2; k++) { f.wr[k] = Math.cos((2 * Math.PI * k) / ny); f.wi[k] = -Math.sin((2 * Math.PI * k) / ny); }
    for (let k = 0; k < ny; k++) { f.cr[k] = Math.cos((Math.PI * k) / (2 * ny)); f.sr[k] = Math.sin((Math.PI * k) / (2 * ny)); }
    s.fft = f;
  } else {
    // otherwise a matrix transform folded about mid-height (even and odd modes use the symmetric and antisymmetric halves)
    if (ny % 2) throw new Error('nsCreate: the number of cells across the stream must be even');
    const hf = ny >> 1; s.hf = hf; s.cs = new Float64Array(ny * hf); s.fe = new Float64Array(hf); s.fo = new Float64Array(hf);
    for (let k = 0; k < ny; k++) for (let j = 0; j < hf; j++) s.cs[k * hf + j] = Math.cos((Math.PI * k * (j + 0.5)) / ny);
  }
  // Thomas factors of ψ(i−1) + ψ(i+1) − (2 + L_k)ψ(i) = D with a Neumann inlet and a Dirichlet (zero pressure) outlet face
  for (let k = 0; k < ny; k++) {
    const L = 2 - 2 * Math.cos((Math.PI * k) / ny); let c = 0;
    for (let i = 0; i < nx; i++) { const b = -(2 + L) + (i === 0 ? 1 : 0) - (i === nx - 1 ? 1 : 0); c = 1 / (b - c); s.tm[i * ny + k] = c; }
  }
  return s;
}

/** In-place radix-2 FFT of f.re, f.im (already in bit-reversed order); sgn = 1 forward, −1 inverse (unscaled). */
function fftCore(f, n, sgn) {
  const re = f.re, im = f.im, wr = f.wr, wi = f.wi;
  for (let i = 0; i < n; i += 2) { const xr = re[i + 1], xi = im[i + 1]; re[i + 1] = re[i] - xr; im[i + 1] = im[i] - xi; re[i] += xr; im[i] += xi; }
  for (let len = 4; len <= n; len <<= 1) {
    const half = len >> 1, st = n / len;
    for (let k = 0; k < half; k++) {
      const c = wr[k * st], d = sgn * wi[k * st];
      for (let a = k; a < n; a += len) {
        const b = a + half, br = re[b], bi = im[b], xr = br * c - bi * d, xi = br * d + bi * c, ar = re[a], ai = im[a];
        re[b] = ar - xr; im[b] = ai - xi; re[a] = ar + xr; im[a] = ai + xi;
      }
    }
  }
}
/** X[k] = Σ x[j]·cos(πk(j+½)/n) of src[o..o+n) into dst[o..o+n). */
function dctFwd(s, src, dst, o) {
  const n = s.ny, f = s.fft;
  if (!f) {
    const cs = s.cs, hf = s.hf, fe = s.fe, fo = s.fo;
    for (let j = 0; j < hf; j++) { const a = src[o + j], b = src[o + n - 1 - j]; fe[j] = a + b; fo[j] = a - b; }
    for (let k = 0; k < n; k++) { const c = k * hf, w = k & 1 ? fo : fe; let a = 0, b = 0, j = 0; for (; j + 1 < hf; j += 2) { a += cs[c + j] * w[j]; b += cs[c + j + 1] * w[j + 1]; } if (j < hf) a += cs[c + j] * w[j]; dst[o + k] = a + b; }
    return;
  }
  const { re, im, rev, cr, sr } = f, h = n >> 1;
  for (let m = 0; m < h; m++) { re[rev[m]] = src[o + 2 * m]; re[rev[n - 1 - m]] = src[o + 2 * m + 1]; }
  im.fill(0); fftCore(f, n, 1);
  for (let k = 0; k < n; k++) dst[o + k] = cr[k] * re[k] + sr[k] * im[k];
}
/** Inverse of dctFwd. */
function dctInv(s, src, dst, o) {
  const n = s.ny, f = s.fft;
  if (!f) {
    const cs = s.cs, hf = s.hf, wk = 2 / n;
    for (let j = 0; j < hf; j++) {
      let e = 0.5 * src[o], d = 0;
      for (let k = 2; k < n; k += 2) e += src[o + k] * cs[k * hf + j];
      for (let k = 1; k < n; k += 2) d += src[o + k] * cs[k * hf + j];
      e *= wk; d *= wk; dst[o + j] = e + d; dst[o + n - 1 - j] = e - d;
    }
    return;
  }
  const { re, im, rev, cr, sr } = f, h = n >> 1, sc = 1 / n;
  re[0] = src[o] * sc; im[0] = 0;
  for (let k = 1; k < n; k++) { const a = src[o + k] * sc, b = -src[o + n - k] * sc, r = rev[k]; re[r] = cr[k] * a - sr[k] * b; im[r] = cr[k] * b + sr[k] * a; }
  fftCore(f, n, -1);
  for (let m = 0; m < h; m++) { dst[o + 2 * m] = re[m]; dst[o + 2 * m + 1] = re[n - 1 - m]; }
}

function ghosts(s) {
  const { nx, ny, su, sv, u, v, U } = s;
  for (let j = 0; j < ny; j++) {
    const r = (j + G) * su;
    u[r + G] = U; u[r + G - 1] = U; u[r + G - 2] = U;
    u[r + G + nx + 1] = u[r + G + nx]; u[r + G + nx + 2] = u[r + G + nx];
    const q = (j + G) * sv;
    v[q + G - 1] = -v[q + G]; v[q + G - 2] = -v[q + G + 1];
    v[q + G + nx] = v[q + G + nx - 1]; v[q + G + nx + 1] = v[q + G + nx - 1];
  }
  for (let i = 0; i < su; i++) {
    u[i + su] = u[i + 2 * su]; u[i] = u[i + 3 * su];
    u[i + (ny + G) * su] = u[i + (ny + G - 1) * su]; u[i + (ny + G + 1) * su] = u[i + (ny + G - 2) * su];
  }
  for (let i = 0; i < sv; i++) {
    v[i + G * sv] = 0; v[i + (ny + G) * sv] = 0;
    v[i + sv] = -v[i + 3 * sv]; v[i] = -v[i + 4 * sv];
    v[i + (ny + G + 1) * sv] = -v[i + (ny + G - 1) * sv]; v[i + (ny + G + 2) * sv] = -v[i + (ny + G - 2) * sv];
  }
}

/** Advection–diffusion predictor (no pressure, no body). Advances s.u, s.v to provisional values and s.t by dt. */
export function nsPredict(s, dt) {
  const { nx, ny, su, sv, u, v, hu, hv, hu0, hv0, h, nu, bu: du, bv: dv } = s, c12 = 1 / (12 * h), dif = nu / (h * h), ck = (s.kk * dt) / (12 * h);
  ghosts(s);
  for (let j = 0; j < ny; j++) {
    for (let i = 1; i < nx; i++) {
      const k = i + G + (j + G) * su, q = i + G + (j + G) * sv, uc = u[k], vc = 0.25 * (v[q - 1] + v[q] + v[q - 1 + sv] + v[q + sv]);
      const xm2 = u[k - 2], xm1 = u[k - 1], xp1 = u[k + 1], xp2 = u[k + 2], ym2 = u[k - 2 * su], ym1 = u[k - su], yp1 = u[k + su], yp2 = u[k + 2 * su];
      hu[k] = -c12 * (uc * (-xp2 + 8 * (xp1 - xm1) + xm2) + vc * (-yp2 + 8 * (yp1 - ym1) + ym2)) + dif * (xp1 + xm1 + yp1 + ym1 - 4 * uc);
      du[k] = -ck * (Math.abs(uc) * (xp2 - 4 * (xp1 + xm1) + 6 * uc + xm2) + Math.abs(vc) * (yp2 - 4 * (yp1 + ym1) + 6 * uc + ym2));
    }
  }
  for (let j = 1; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const k = i + G + (j + G) * sv, q = i + G + (j + G) * su, vc = v[k], uc = 0.25 * (u[q - su] + u[q + 1 - su] + u[q] + u[q + 1]);
      const xm2 = v[k - 2], xm1 = v[k - 1], xp1 = v[k + 1], xp2 = v[k + 2], ym2 = v[k - 2 * sv], ym1 = v[k - sv], yp1 = v[k + sv], yp2 = v[k + 2 * sv];
      hv[k] = -c12 * (uc * (-xp2 + 8 * (xp1 - xm1) + xm2) + vc * (-yp2 + 8 * (yp1 - ym1) + ym2)) + dif * (xp1 + xm1 + yp1 + ym1 - 4 * vc);
      dv[k] = -ck * (Math.abs(uc) * (xp2 - 4 * (xp1 + xm1) + 6 * vc + xm2) + Math.abs(vc) * (yp2 - 4 * (yp1 + ym1) + 6 * vc + ym2));
    }
  }
  const a = s.nStep === 0 ? dt : 1.5 * dt, b = s.nStep === 0 ? 0 : -0.5 * dt, cfl = (s.U * dt) / h;
  for (let j = 0; j < ny; j++) {
    const r = G + (j + G) * su;
    u[r + nx] -= cfl * (u[r + nx] - u[r + nx - 1]); // convective outflow
    for (let i = 1; i < nx; i++) { const k = r + i; u[k] += a * hu[k] + b * hu0[k] + du[k]; }
  }
  for (let j = 1; j < ny; j++) { const r = G + (j + G) * sv; for (let i = 0; i < nx; i++) { const k = r + i; v[k] += a * hv[k] + b * hv0[k] + dv[k]; } }
  s.hu = hu0; s.hu0 = hu; s.hv = hv0; s.hv0 = hv;
  s.p.fill(0); s.t += dt; s.nStep++;
}

/** Signed distance from point (px, py) to the body surface (negative inside). */
function sdist(b, px, py) {
  const dx = px - b.x, dy = py - b.y;
  if (b.type === 'circle') return Math.sqrt(dx * dx + dy * dy) - b.R;
  const xi = dx * b.ca + dy * b.sa - b.xm, eta = -dx * b.sa + dy * b.ca, e = Math.abs(xi) - b.half;
  return (e > 0 ? Math.sqrt(e * e + eta * eta) : Math.abs(eta)) - 0.5 * b.thick;
}

/**
 * Describe a rigid body. Circle: { type: 'circle', R }. Rounded flat plate: { type: 'plate', c, thick, xm } with xm the
 * position of mid-chord along the chord measured from the reference point (pivot), positive downstream.
 * State (mutable): x, y of the reference point, ang (counter-clockwise), vx, vy, om.
 */
export function nsBody(o) {
  const b = { type: o.type, R: o.R ?? 0, c: o.c ?? 0, thick: o.thick ?? 0, xm: o.xm ?? 0, x: o.x, y: o.y, ang: o.ang ?? 0, vx: 0, vy: 0, om: 0, ca: 1, sa: 0, half: 0, ext: 0 };
  b.half = Math.max(0, 0.5 * (b.c - b.thick)); b.ext = b.type === 'circle' ? b.R : Math.abs(b.xm) + 0.5 * b.c;
  return b;
}

/**
 * One direct-forcing pass towards the body's current rigid velocity. Returns the load on the body from this pass
 * { Fx, Fy, Mz } (moment about the reference point, counter-clockwise) and the momentum of the fluid inside the mask
 * { Px, Py, Lz, area }. The hydrodynamic load is Σ passes of F plus the rate of change of the interior momentum.
 */
export function nsForce(s, b, dt) {
  const { nx, ny, su, sv, u, v, h, rho } = s, m = rho * h * h, out = { Fx: 0, Fy: 0, Mz: 0, Px: 0, Py: 0, Lz: 0, area: 0 };
  b.ca = Math.cos(b.ang); b.sa = Math.sin(b.ang);
  const i0 = Math.max(1, Math.floor((b.x - b.ext) / h) - 2), i1 = Math.min(nx - 1, Math.ceil((b.x + b.ext) / h) + 2), j0 = Math.max(1, Math.floor((b.y - b.ext) / h) - 2), j1 = Math.min(ny - 1, Math.ceil((b.y + b.ext) / h) + 2);
  let sx = 0, sy = 0, mz = 0, px = 0, py = 0, lz = 0, ar = 0;
  for (let j = j0; j < j1; j++) for (let i = i0; i <= i1; i++) {
    const X = i * h, Y = (j + 0.5) * h, d = sdist(b, X, Y); if (d >= 0.5 * h) continue;
    const chi = d <= -0.5 * h ? 1 : 0.5 - d / h, k = i + G + (j + G) * su, ry = Y - b.y, ub = b.vx - b.om * ry, du = chi * (ub - u[k]);
    u[k] += du; sx += du; mz -= ry * du; px += chi * ub; lz -= ry * chi * ub; ar += 0.5 * chi;
  }
  for (let j = j0; j <= j1; j++) for (let i = i0; i < i1; i++) {
    const X = (i + 0.5) * h, Y = j * h, d = sdist(b, X, Y); if (d >= 0.5 * h) continue;
    const chi = d <= -0.5 * h ? 1 : 0.5 - d / h, k = i + G + (j + G) * sv, rx = X - b.x, vb = b.vy + b.om * rx, dv = chi * (vb - v[k]);
    v[k] += dv; sy += dv; mz += rx * dv; py += chi * vb; lz += rx * chi * vb; ar += 0.5 * chi;
  }
  out.Fx = (-m * sx) / dt; out.Fy = (-m * sy) / dt; out.Mz = (-m * mz) / dt; out.Px = m * px; out.Py = m * py; out.Lz = m * lz; out.area = ar * h * h;
  return out;
}

/** Momentum { Px, Py, Lz, area } that fluid moving rigidly with the body inside its mask would have (no change to the field). */
export function nsInterior(s, b) {
  const { nx, ny, h, rho } = s, m = rho * h * h; b.ca = Math.cos(b.ang); b.sa = Math.sin(b.ang);
  const i0 = Math.max(1, Math.floor((b.x - b.ext) / h) - 2), i1 = Math.min(nx - 1, Math.ceil((b.x + b.ext) / h) + 2), j0 = Math.max(1, Math.floor((b.y - b.ext) / h) - 2), j1 = Math.min(ny - 1, Math.ceil((b.y + b.ext) / h) + 2);
  let px = 0, py = 0, lz = 0, ar = 0;
  for (let j = j0; j < j1; j++) for (let i = i0; i <= i1; i++) { const X = i * h, Y = (j + 0.5) * h, d = sdist(b, X, Y); if (d >= 0.5 * h) continue; const chi = d <= -0.5 * h ? 1 : 0.5 - d / h, ry = Y - b.y, ub = b.vx - b.om * ry; px += chi * ub; lz -= ry * chi * ub; ar += 0.5 * chi; }
  for (let j = j0; j <= j1; j++) for (let i = i0; i < i1; i++) { const X = (i + 0.5) * h, Y = j * h, d = sdist(b, X, Y); if (d >= 0.5 * h) continue; const chi = d <= -0.5 * h ? 1 : 0.5 - d / h, rx = X - b.x, vb = b.vy + b.om * rx; py += chi * vb; lz += rx * chi * vb; ar += 0.5 * chi; }
  return { Px: m * px, Py: m * py, Lz: m * lz, area: ar * h * h };
}

/** Pressure projection: makes the velocity discretely divergence-free and accumulates the pressure of this step in s.p. */
export function nsProject(s, dt) {
  const { nx, ny, su, sv, u, v, tm, w1, w2, h, p } = s;
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) { const k = i + G + (j + G) * su, q = i + G + (j + G) * sv; w1[i * ny + j] = u[k + 1] - u[k] + v[q + sv] - v[q]; }
  for (let i = 0; i < nx; i++) dctFwd(s, w1, w2, i * ny);
  for (let k = 0; k < ny; k++) w2[k] *= tm[k];
  for (let i = 1; i < nx; i++) { const r = i * ny; for (let k = 0; k < ny; k++) w2[r + k] = (w2[r + k] - w2[r - ny + k]) * tm[r + k]; }
  for (let i = nx - 2; i >= 0; i--) { const r = i * ny; for (let k = 0; k < ny; k++) w2[r + k] -= tm[r + k] * w2[r + ny + k]; }
  for (let i = 0; i < nx; i++) dctInv(s, w2, w1, i * ny);
  const pf = (s.rho * h) / dt;
  for (let j = 0; j < ny; j++) {
    const r = G + (j + G) * su;
    for (let i = 1; i < nx; i++) u[r + i] -= w1[i * ny + j] - w1[(i - 1) * ny + j];
    u[r + nx] += 2 * w1[(nx - 1) * ny + j];
  }
  for (let j = 1; j < ny; j++) { const r = G + (j + G) * sv; for (let i = 0; i < nx; i++) v[r + i] -= w1[i * ny + j] - w1[i * ny + j - 1]; }
  for (let k = 0; k < nx * ny; k++) p[k] += pf * w1[k];
}

/** Copy the velocity field into numbered storage slot k (allocated on first use). */
export function nsStore(s, k) {
  const sl = s.slots || (s.slots = []); if (!sl[k]) sl[k] = { u: new Float64Array(s.u.length), v: new Float64Array(s.v.length) };
  sl[k].u.set(s.u); sl[k].v.set(s.v);
}
/** Restore the velocity field from slot k. */
export function nsRecall(s, k) { s.u.set(s.slots[k].u); s.v.set(s.slots[k].v); }
/** Turn slot k into the difference slot k − current field (a stored sensitivity field). */
export function nsDiff(s, k) { const u = s.u, v = s.v, su = s.slots[k].u, sv = s.slots[k].v; for (let i = 0; i < u.length; i++) su[i] -= u[i]; for (let i = 0; i < v.length; i++) sv[i] -= v[i]; }
/** Linear combination in place: field ← c0·field + Σ c[d]·slot[ks[d]]. Used to superpose sub-solutions that are affine in the body motion. */
export function nsCombine(s, c0, ks, c) {
  const u = s.u, v = s.v;
  if (c0 !== 1) { for (let k = 0; k < u.length; k++) u[k] *= c0; for (let k = 0; k < v.length; k++) v[k] *= c0; }
  for (let d = 0; d < ks.length; d++) { const a = c[d], su = s.slots[ks[d]].u, sv = s.slots[ks[d]].v; if (a === 0) continue; for (let k = 0; k < u.length; k++) u[k] += a * su[k]; for (let k = 0; k < v.length; k++) v[k] += a * sv[k]; }
}

/** Largest absolute discrete divergence (per unit time) — zero to round-off after a projection. */
export function nsDivergence(s) {
  const { nx, ny, su, sv, u, v, h } = s; let m = 0;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const k = i + G + (j + G) * su, q = i + G + (j + G) * sv, d = Math.abs(u[k + 1] - u[k] + v[q + sv] - v[q]) / h; if (d > m) m = d; }
  return m;
}

/** Largest velocity magnitude component (for the Courant number). */
export function nsMaxSpeed(s) {
  const { nx, ny, su, sv, u, v } = s; let m = 0;
  for (let j = 0; j < ny; j++) for (let i = 0; i <= nx; i++) { const a = Math.abs(u[i + G + (j + G) * su]); if (a > m) m = a; }
  for (let j = 0; j <= ny; j++) for (let i = 0; i < nx; i++) { const a = Math.abs(v[i + G + (j + G) * sv]); if (a > m) m = a; }
  return m;
}

/** Vorticity at cell corners, sampled every `skip` cells: { x[], y[], z[ny][nx] } as plain arrays. */
export function nsVorticity(s, skip = 1) {
  const { nx, ny, su, sv, u, v, h } = s, x = [], y = [], z = [];
  for (let i = 1; i < nx; i += skip) x.push(i * h);
  for (let j = 1; j < ny; j += skip) {
    y.push(j * h); const row = [];
    for (let i = 1; i < nx; i += skip) { const k = i + G + (j + G) * su, q = i + G + (j + G) * sv; row.push((v[q] - v[q - 1] - u[k] + u[k - su]) / h); }
    z.push(row);
  }
  return { x, y, z };
}

/** Cross-stream velocity v at a point (nearest v-face), e.g. as a wake probe. */
export function nsProbeV(s, x, y) {
  const i = Math.min(s.nx - 1, Math.max(0, Math.floor(x / s.h))), j = Math.min(s.ny, Math.max(0, Math.round(y / s.h)));
  return s.v[i + G + (j + G) * s.sv];
}
