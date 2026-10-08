// Shared numerical kernel used by every suite: linear algebra, eigen-solvers, ODE
// integration, root finding, optimisation, statistics, sampling and FFT.
// Matrices are plain arrays of row arrays. Everything is dependency-free and pure.

export const PI = Math.PI;
export const deg = (r) => (r * 180) / PI;
export const rad = (d) => (d * PI) / 180;
export const clamp = (x, a, b) => (x < a ? a : x > b ? b : x);
export const lerp = (a, b, t) => a + (b - a) * t;
export const sq = (x) => x * x;
export const sign = (x) => (x < 0 ? -1 : 1);

// ---------- vectors ----------
export function linspace(a, b, n) {
  if (n < 2) return [a];
  const out = new Array(n), h = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out[i] = a + h * i;
  return out;
}
/** n values geometrically spaced from a to b (end VALUES, not exponents). */
export const logspace = (a, b, n) => linspace(Math.log10(a), Math.log10(b), n).map((v) => 10 ** v);
export const range = (n, f = (i) => i) => Array.from({ length: n }, (_, i) => f(i));
export const sum = (a) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s; };
export const mean = (a) => (a.length ? sum(a) / a.length : NaN);
export function variance(a) {
  if (a.length < 2) return 0;
  const m = mean(a); let s = 0;
  for (const v of a) s += (v - m) * (v - m);
  return s / (a.length - 1);
}
export const std = (a) => Math.sqrt(variance(a));
export const amin = (a) => { let m = Infinity; for (const v of a) if (v < m) m = v; return m; };
export const amax = (a) => { let m = -Infinity; for (const v of a) if (v > m) m = v; return m; };
export const argmax = (a) => { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] > a[k]) k = i; return k; };
export const argmin = (a) => { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] < a[k]) k = i; return k; };
export const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
export const norm = (a) => Math.sqrt(dot(a, a));
export const vadd = (a, b, s = 1) => a.map((v, i) => v + s * b[i]);
export const vscale = (a, s) => a.map((v) => v * s);
export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

// ---------- interpolation ----------
/** Piecewise-linear interpolation, clamped at the ends. xs must be ascending. */
export function interp1(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] > x) hi = m; else lo = m; }
  const t = (x - xs[lo]) / (xs[hi] - xs[lo]);
  return ys[lo] + t * (ys[hi] - ys[lo]);
}
/** Bilinear interpolation on a rectilinear table Z[iy][ix], clamped. */
export function interp2(xs, ys, Z, x, y) {
  const col = ys.map((_, j) => interp1(xs, Z[j], x));
  return interp1(ys, col, y);
}
/** Natural cubic spline; returns an evaluator function. */
export function spline(xs, ys) {
  const n = xs.length, h = [], al = [], l = [1], mu = [0], z = [0];
  const c = new Array(n).fill(0), b = [], d = [];
  for (let i = 0; i < n - 1; i++) h[i] = xs[i + 1] - xs[i];
  for (let i = 1; i < n - 1; i++) al[i] = (3 / h[i]) * (ys[i + 1] - ys[i]) - (3 / h[i - 1]) * (ys[i] - ys[i - 1]);
  for (let i = 1; i < n - 1; i++) {
    l[i] = 2 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
    mu[i] = h[i] / l[i];
    z[i] = (al[i] - h[i - 1] * z[i - 1]) / l[i];
  }
  for (let j = n - 2; j >= 0; j--) {
    c[j] = (z[j] || 0) - (mu[j] || 0) * c[j + 1];
    b[j] = (ys[j + 1] - ys[j]) / h[j] - (h[j] * (c[j + 1] + 2 * c[j])) / 3;
    d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
  }
  return (x) => {
    let i = 0, hi = n - 1;
    if (x <= xs[0]) i = 0; else if (x >= xs[n - 1]) i = n - 2;
    else { while (hi - i > 1) { const m = (i + hi) >> 1; if (xs[m] > x) hi = m; else i = m; } }
    const dx = x - xs[i];
    return ys[i] + dx * (b[i] + dx * (c[i] + dx * d[i]));
  };
}

// ---------- root finding ----------
/** Brent's method on a bracketing interval [a,b]. Throws if the root is not bracketed. */
export function brent(f, a, b, tol = 1e-10, maxIter = 200) {
  let fa = f(a), fb = f(b);
  if (fa === 0) return a;
  if (fb === 0) return b;
  if (fa * fb > 0) throw new Error('brent: root not bracketed');
  let c = a, fc = fa, d = b - a, e = d;
  for (let i = 0; i < maxIter; i++) {
    if (fb * fc > 0) { c = a; fc = fa; d = e = b - a; }
    if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
    const tol1 = 2 * Number.EPSILON * Math.abs(b) + 0.5 * tol, xm = 0.5 * (c - b);
    if (Math.abs(xm) <= tol1 || fb === 0) return b;
    if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
      const s = fb / fa; let p, q;
      if (a === c) { p = 2 * xm * s; q = 1 - s; }
      else { const qq = fa / fc, r = fb / fc; p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1)); q = (qq - 1) * (r - 1) * (s - 1); }
      if (p > 0) q = -q;
      p = Math.abs(p);
      if (2 * p < Math.min(3 * xm * q - Math.abs(tol1 * q), Math.abs(e * q))) { e = d; d = p / q; } else { d = xm; e = d; }
    } else { d = xm; e = d; }
    a = b; fa = fb;
    b += Math.abs(d) > tol1 ? d : tol1 * sign(xm);
    fb = f(b);
  }
  return b;
}
/** Scan [a,b] in n steps for a sign change, then refine with Brent. Returns NaN when none is found. */
export function findRoot(f, a, b, n = 60, tol = 1e-9) {
  let x0 = a, f0 = f(a);
  for (let i = 1; i <= n; i++) {
    const x1 = a + ((b - a) * i) / n, f1 = f(x1);
    if (Number.isFinite(f0) && Number.isFinite(f1) && f0 * f1 <= 0) return brent(f, x0, x1, tol);
    x0 = x1; f0 = f1;
  }
  return NaN;
}
/** Damped Newton for systems F(x)=0 with a finite-difference Jacobian. */
export function fsolve(F, x0, { tol = 1e-9, maxIter = 60 } = {}) {
  let x = x0.slice(), f = F(x), nf = norm(f);
  const n = x.length;
  for (let it = 0; it < maxIter; it++) {
    if (nf < tol) return { x, converged: true, iterations: it, residual: nf };
    const J = zeros(f.length, n);
    for (let j = 0; j < n; j++) {
      const h = 1e-7 * Math.max(1, Math.abs(x[j])), xp = x.slice(); xp[j] += h;
      const fp = F(xp);
      for (let i = 0; i < f.length; i++) J[i][j] = (fp[i] - f[i]) / h;
    }
    let dx;
    try { dx = f.length === n ? solve(J, f.map((v) => -v)) : lstsq(J, f.map((v) => -v)); } catch { break; }
    let lam = 1, ok = false;
    for (let k = 0; k < 25; k++) {
      const xn = vadd(x, dx, lam), fn = F(xn), nn = norm(fn);
      if (Number.isFinite(nn) && nn < nf) { x = xn; f = fn; nf = nn; ok = true; break; }
      lam *= 0.5;
    }
    if (!ok) break;
  }
  return { x, converged: nf < tol * 100, iterations: maxIter, residual: nf };
}

// ---------- quadrature ----------
export function trapz(x, y) { let s = 0; for (let i = 1; i < x.length; i++) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]); return s; }
export function simpson(f, a, b, n = 200) {
  if (n % 2) n++;
  const h = (b - a) / n; let s = f(a) + f(b);
  for (let i = 1; i < n; i++) s += f(a + i * h) * (i % 2 ? 4 : 2);
  return (s * h) / 3;
}
export function cumtrapz(x, y) { const o = [0]; for (let i = 1; i < x.length; i++) o[i] = o[i - 1] + 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]); return o; }

// ---------- dense linear algebra ----------
export const zeros = (n, m = n) => Array.from({ length: n }, () => new Array(m).fill(0));
export const eye = (n) => { const I = zeros(n); for (let i = 0; i < n; i++) I[i][i] = 1; return I; };
export const transpose = (A) => A[0].map((_, j) => A.map((r) => r[j]));
export const matvec = (A, x) => A.map((r) => dot(r, x));
export function matmul(A, B) {
  const n = A.length, m = B[0].length, k = B.length, C = zeros(n, m);
  for (let i = 0; i < n; i++) for (let p = 0; p < k; p++) { const a = A[i][p]; if (a !== 0) for (let j = 0; j < m; j++) C[i][j] += a * B[p][j]; }
  return C;
}
export const madd = (A, B, s = 1) => A.map((r, i) => r.map((v, j) => v + s * B[i][j]));
export const mscale = (A, s) => A.map((r) => r.map((v) => v * s));
/** LU factorisation with partial pivoting. */
export function lu(A) {
  const n = A.length, M = A.map((r) => r.slice()), piv = range(n); let sgn = 1;
  for (let k = 0; k < n; k++) {
    let p = k, mx = Math.abs(M[k][k]);
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > mx) { mx = Math.abs(M[i][k]); p = i; }
    if (mx < 1e-300) throw new Error('Singular matrix');
    if (p !== k) { [M[p], M[k]] = [M[k], M[p]]; [piv[p], piv[k]] = [piv[k], piv[p]]; sgn = -sgn; }
    for (let i = k + 1; i < n; i++) {
      const f = (M[i][k] /= M[k][k]);
      if (f !== 0) for (let j = k + 1; j < n; j++) M[i][j] -= f * M[k][j];
    }
  }
  return { M, piv, sgn };
}
export function luSolve({ M, piv }, b) {
  const n = M.length, x = new Array(n);
  for (let i = 0; i < n; i++) { let s = b[piv[i]]; for (let j = 0; j < i; j++) s -= M[i][j] * x[j]; x[i] = s; }
  for (let i = n - 1; i >= 0; i--) { let s = x[i]; for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j]; x[i] = s / M[i][i]; }
  return x;
}
export const solve = (A, b) => luSolve(lu(A), b);
export function inv(A) {
  const n = A.length, f = lu(A), cols = range(n, (j) => luSolve(f, range(n, (i) => (i === j ? 1 : 0))));
  return transpose(cols);
}
export function det(A) { try { const { M, sgn } = lu(A); let d = sgn; for (let i = 0; i < M.length; i++) d *= M[i][i]; return d; } catch { return 0; } }
/** Least squares via the normal equations with a tiny Tikhonov term for robustness. */
export function lstsq(A, b, ridge = 0) {
  const At = transpose(A), N = matmul(At, A);
  for (let i = 0; i < N.length; i++) N[i][i] += ridge + 1e-14 * (Math.abs(N[i][i]) + 1e-300);
  return solve(N, matvec(At, b));
}
/** Thomas algorithm: a=sub, b=diag, c=super, d=rhs. */
export function solveTridiag(a, b, c, d) {
  const n = d.length, cp = new Array(n), dp = new Array(n), x = new Array(n);
  cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0];
  for (let i = 1; i < n; i++) { const m = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / m; dp[i] = (d[i] - a[i] * dp[i - 1]) / m; }
  x[n - 1] = dp[n - 1];
  for (let i = n - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
  return x;
}
export function cholesky(A) {
  const n = A.length, L = zeros(n);
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
    let s = A[i][j];
    for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
    if (i === j) { if (s <= 0) throw new Error('Matrix not positive definite'); L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j];
  }
  return L;
}
/** Polynomial least-squares fit; returns coefficients c0 + c1 x + ... */
export function polyfit(x, y, order) { return lstsq(x.map((v) => range(order + 1, (k) => v ** k)), y); }
export const polyval = (c, x) => { let s = 0; for (let k = c.length - 1; k >= 0; k--) s = s * x + c[k]; return s; };

// ---------- eigenvalue problems ----------
/** Symmetric eigenproblem by cyclic Jacobi. Returns ascending values and vectors[k] (k-th eigenvector). */
export function eigSym(Ain) {
  const n = Ain.length, A = Ain.map((r) => r.slice()), V = eye(n);
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0, diag = 0;
    for (let i = 0; i < n; i++) { diag += A[i][i] * A[i][i]; for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j]; }
    if (off <= 1e-28 * (diag + off) || off < 1e-300) break; // relative, so it works at any matrix scale
    for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(A[p][q]) < 1e-300) continue;
      const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
      const t = sign(th) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq; }
      for (let k = 0; k < n; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk; }
      for (let k = 0; k < n; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
    }
  }
  const idx = range(n).sort((a, b) => A[a][a] - A[b][b]);
  return { values: idx.map((i) => A[i][i]), vectors: idx.map((i) => V.map((r) => r[i])) };
}
/** Generalised symmetric problem K x = lambda M x (M positive definite). Vectors are M-orthonormal. */
export function eigGenSym(K, M) {
  const n = K.length, L = cholesky(M), Li = zeros(n);
  for (let j = 0; j < n; j++) {
    Li[j][j] = 1 / L[j][j];
    for (let i = j + 1; i < n; i++) { let s = 0; for (let k = j; k < i; k++) s -= L[i][k] * Li[k][j]; Li[i][j] = s / L[i][i]; }
  }
  const C = matmul(matmul(Li, K), transpose(Li));
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) C[i][j] = C[j][i] = 0.5 * (C[i][j] + C[j][i]);
  const e = eigSym(C), LiT = transpose(Li);
  return { values: e.values, vectors: e.vectors.map((v) => matvec(LiT, v)) };
}

// complex scalars as [re, im]
export const C = {
  add: (a, b) => [a[0] + b[0], a[1] + b[1]],
  sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
  mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
  div: (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
  abs: (a) => Math.hypot(a[0], a[1]),
  arg: (a) => Math.atan2(a[1], a[0]),
  conj: (a) => [a[0], -a[1]],
  scale: (a, s) => [a[0] * s, a[1] * s],
  exp: (a) => { const e = Math.exp(a[0]); return [e * Math.cos(a[1]), e * Math.sin(a[1])]; },
  sqrt: (a) => { const r = Math.hypot(a[0], a[1]), re = Math.sqrt((r + a[0]) / 2), im = Math.sqrt(Math.max(0, (r - a[0]) / 2)); return [re, a[1] < 0 ? -im : im]; },
  from: (x) => (Array.isArray(x) ? x : [x, 0]),
};
/** Complex linear solve. A is a matrix of [re,im] (or real numbers), b a vector of the same. */
export function csolve(Ain, bin) {
  const n = Ain.length, A = Ain.map((r) => r.map(C.from)), b = bin.map(C.from);
  for (let k = 0; k < n; k++) {
    let p = k, mx = C.abs(A[k][k]);
    for (let i = k + 1; i < n; i++) { const v = C.abs(A[i][k]); if (v > mx) { mx = v; p = i; } }
    if (mx < 1e-300) throw new Error('Singular complex matrix');
    [A[p], A[k]] = [A[k], A[p]]; [b[p], b[k]] = [b[k], b[p]];
    for (let i = k + 1; i < n; i++) {
      const f = C.div(A[i][k], A[k][k]);
      for (let j = k; j < n; j++) A[i][j] = C.sub(A[i][j], C.mul(f, A[k][j]));
      b[i] = C.sub(b[i], C.mul(f, b[k]));
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let j = i + 1; j < n; j++) s = C.sub(s, C.mul(A[i][j], x[j])); x[i] = C.div(s, A[i][i]); }
  return x;
}
/**
 * Eigenvalues of a general (real or complex) square matrix by shifted QR on the
 * upper-Hessenberg form with deflation. Entries may be numbers or [re,im].
 * Returns [[re,im], ...] sorted by imaginary part then real part.
 */
export function eig(Ain) {
  let n = Ain.length;
  const H = Ain.map((r) => r.map((v) => C.from(v).slice()));
  // Hessenberg reduction by Gaussian similarity transforms with pivoting
  for (let m = 1; m < n - 1; m++) {
    let p = m, mx = 0;
    for (let i = m; i < n; i++) { const v = C.abs(H[i][m - 1]); if (v > mx) { mx = v; p = i; } }
    if (mx === 0) continue;
    if (p !== m) { [H[p], H[m]] = [H[m], H[p]]; for (let i = 0; i < n; i++) [H[i][p], H[i][m]] = [H[i][m], H[i][p]]; }
    for (let i = m + 1; i < n; i++) {
      const f = C.div(H[i][m - 1], H[m][m - 1]);
      if (f[0] === 0 && f[1] === 0) continue;
      for (let j = 0; j < n; j++) H[i][j] = C.sub(H[i][j], C.mul(f, H[m][j]));
      for (let j = 0; j < n; j++) H[j][m] = C.add(H[j][m], C.mul(f, H[j][i]));
    }
  }
  const vals = [];
  let iter = 0;
  while (n > 0) {
    if (n === 1) { vals.push(H[0][0]); break; }
    // deflation check
    let l = n - 1;
    for (; l > 0; l--) {
      const s = C.abs(H[l - 1][l - 1]) + C.abs(H[l][l]) || 1;
      if (C.abs(H[l][l - 1]) < 1e-14 * s) { H[l][l - 1] = [0, 0]; break; }
    }
    if (l === n - 1) { vals.push(H[n - 1][n - 1]); n--; iter = 0; continue; }
    if (++iter > 500) throw new Error('eig: QR iteration did not converge');
    // Wilkinson shift from trailing 2x2 (exceptional shift every 10 iterations)
    const a = H[n - 2][n - 2], b = H[n - 2][n - 1], c = H[n - 1][n - 2], d = H[n - 1][n - 1];
    let mu;
    if (iter % 11 === 10) mu = [C.abs(c) + C.abs(H[n - 2][n - 3] || [0, 0]), 0];
    else {
      const tr = C.add(a, d), dt = C.sub(C.mul(a, d), C.mul(b, c));
      const disc = C.sqrt(C.sub(C.mul(tr, tr), C.scale(dt, 4)));
      const l1 = C.scale(C.add(tr, disc), 0.5), l2 = C.scale(C.sub(tr, disc), 0.5);
      mu = C.abs(C.sub(l1, d)) < C.abs(C.sub(l2, d)) ? l1 : l2;
    }
    // QR step on active block l..n-1 using Givens rotations
    for (let i = l; i < n; i++) H[i][i] = C.sub(H[i][i], mu);
    const rot = [];
    for (let k = l; k < n - 1; k++) {
      const x = H[k][k], y = H[k + 1][k], r = Math.hypot(C.abs(x), C.abs(y));
      if (r === 0) { rot.push(null); continue; }
      const cs = C.scale(x, 1 / r), sn = C.scale(y, 1 / r); // G = [conj(cs) conj(sn); -sn cs]
      rot.push([cs, sn]);
      for (let j = k; j < n; j++) {
        const u = H[k][j], v = H[k + 1][j];
        H[k][j] = C.add(C.mul(C.conj(cs), u), C.mul(C.conj(sn), v));
        H[k + 1][j] = C.sub(C.mul(cs, v), C.mul(sn, u));
      }
    }
    for (let k = l; k < n - 1; k++) {
      const g = rot[k - l]; if (!g) continue;
      const [cs, sn] = g;
      for (let i = l; i <= Math.min(k + 2, n - 1); i++) {
        const u = H[i][k], v = H[i][k + 1];
        H[i][k] = C.add(C.mul(u, cs), C.mul(v, sn));
        H[i][k + 1] = C.sub(C.mul(v, C.conj(cs)), C.mul(u, C.conj(sn)));
      }
    }
    for (let i = l; i < n; i++) H[i][i] = C.add(H[i][i], mu);
  }
  return vals.map((v) => [v[0], Math.abs(v[1]) < 1e-12 * (1 + Math.abs(v[0])) ? 0 : v[1]]).sort((p, q) => p[1] - q[1] || p[0] - q[0]);
}
/** Roots of c0 + c1 x + ... + cn x^n via the companion matrix. Returns [[re,im],...]. */
export function polyRoots(c) {
  let n = c.length - 1;
  while (n > 0 && c[n] === 0) n--;
  if (n < 1) return [];
  const A = zeros(n);
  for (let i = 1; i < n; i++) A[i][i - 1] = 1;
  for (let i = 0; i < n; i++) A[i][n - 1] = -c[i] / c[n];
  return eig(A);
}

// ---------- control helpers ----------
/** Continuous Lyapunov equation A X + X A' + Q = 0 (Kronecker solve; small n). */
export function lyap(A, Q) {
  const n = A.length, N = n * n, K = zeros(N), b = new Array(N);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const r = i * n + j; b[r] = -Q[i][j];
    for (let k = 0; k < n; k++) { K[r][k * n + j] += A[i][k]; K[r][i * n + k] += A[j][k]; }
  }
  const x = solve(K, b);
  return range(n, (i) => range(n, (j) => x[i * n + j]));
}
/** Continuous algebraic Riccati equation A'P + PA - P B R^-1 B' P + Q = 0 via the matrix sign function. */
export function care(A, B, Q, R) {
  const n = A.length, G = matmul(matmul(B, inv(R)), transpose(B)), At = transpose(A);
  let Z = zeros(2 * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { Z[i][j] = A[i][j]; Z[i][j + n] = -G[i][j]; Z[i + n][j] = -Q[i][j]; Z[i + n][j + n] = -At[i][j]; }
  for (let it = 0; it < 100; it++) {
    const Zi = inv(Z), c = Math.abs(det(Z)) ** (-1 / (2 * n)) || 1;
    const Zn = madd(mscale(Z, c), mscale(Zi, 1 / c)).map((r) => r.map((v) => v / 2));
    let diff = 0, nz = 0;
    for (let i = 0; i < 2 * n; i++) for (let j = 0; j < 2 * n; j++) { diff += (Zn[i][j] - Z[i][j]) ** 2; nz += Zn[i][j] ** 2; }
    Z = Zn;
    if (diff < 1e-22 * nz) break;
  }
  // [W12; W22 + I] P = -[W11 + I; W21]
  const M = zeros(2 * n, n), Nn = zeros(2 * n, n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    M[i][j] = Z[i][j + n]; M[i + n][j] = Z[i + n][j + n] + (i === j ? 1 : 0);
    Nn[i][j] = -(Z[i][j] + (i === j ? 1 : 0)); Nn[i + n][j] = -Z[i + n][j];
  }
  const Mt = transpose(M), MtM = matmul(Mt, M), rhs = matmul(Mt, Nn), f = lu(MtM);
  const P = transpose(range(n, (j) => luSolve(f, rhs.map((r) => r[j]))));
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) P[i][j] = P[j][i] = 0.5 * (P[i][j] + P[j][i]);
  return P;
}
/** LQR gain K = R^-1 B' P. */
export function lqr(A, B, Q, R) { const P = care(A, B, Q, R); return { K: matmul(matmul(inv(R), transpose(B)), P), P }; }

// ---------- ODE integration ----------
/** Fixed-step classical RK4. f(t,y)->dy. Returns {t:[], y:[state arrays]}. */
export function rk4(f, t0, y0, t1, n) {
  const h = (t1 - t0) / n, T = [t0], Y = [y0.slice()];
  let y = y0.slice(), t = t0;
  for (let i = 0; i < n; i++) {
    const k1 = f(t, y), k2 = f(t + h / 2, vadd(y, k1, h / 2)), k3 = f(t + h / 2, vadd(y, k2, h / 2)), k4 = f(t + h, vadd(y, k3, h));
    y = y.map((v, j) => v + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
    t = t0 + (i + 1) * h; T.push(t); Y.push(y.slice());
  }
  return { t: T, y: Y };
}
/**
 * Adaptive Dormand–Prince RK45. opts: rtol, atol, hmax, h0, maxSteps, stop(t,y)->bool (terminates when true).
 */
export function rk45(f, t0, y0, t1, opts = {}) {
  const { rtol = 1e-6, atol = 1e-9, hmax = Infinity, maxSteps = 200000, stop } = opts;
  const c = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
  const a = [[], [1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9], [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729], [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656], [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]];
  const b5 = a[6], e = [71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40];
  const n = y0.length, T = [t0], Y = [y0.slice()];
  let t = t0, y = y0.slice(), h = opts.h0 || Math.min(hmax, (t1 - t0) / 100), k = [f(t, y)], steps = 0;
  while (t < t1 - 1e-14 * Math.abs(t1) && steps++ < maxSteps) {
    if (t + h > t1) h = t1 - t;
    for (let s = 1; s < 7; s++) {
      const ys = y.slice();
      for (let j = 0; j < s; j++) { const as = a[s][j]; if (as) for (let i = 0; i < n; i++) ys[i] += h * as * k[j][i]; }
      k[s] = f(t + c[s] * h, ys);
      if (s === 6) k.yn = ys;
    }
    let err = 0;
    for (let i = 0; i < n; i++) {
      let ei = 0; for (let s = 0; s < 7; s++) ei += e[s] * k[s][i];
      const sc = atol + rtol * Math.max(Math.abs(y[i]), Math.abs(k.yn[i]));
      err += ((h * ei) / sc) ** 2;
    }
    err = Math.sqrt(err / n);
    if (err <= 1 || h < 1e-13) {
      t += h; y = k.yn; T.push(t); Y.push(y.slice());
      const kl = k[6]; k = [kl];
      if (stop && stop(t, y)) break;
    }
    h = Math.min(hmax, h * clamp(0.9 * (err || 1e-10) ** -0.2, 0.2, 5));
    if (!Number.isFinite(h)) break;
  }
  return { t: T, y: Y, steps };
}
/** Backward Euler with Newton iterations for stiff systems. */
export function beuler(f, t0, y0, t1, n) {
  const h = (t1 - t0) / n, T = [t0], Y = [y0.slice()];
  let y = y0.slice();
  for (let i = 0; i < n; i++) {
    const tn = t0 + (i + 1) * h, yo = y;
    y = fsolve((z) => z.map((v, j) => v - yo[j] - h * f(tn, z)[j]), yo, { tol: 1e-10, maxIter: 30 }).x;
    T.push(tn); Y.push(y.slice());
  }
  return { t: T, y: Y };
}
/** Newmark-beta (average acceleration) for M a + C v + K u = F(t). Returns {t,u,v,a}. */
export function newmark(M, Cm, K, F, u0, v0, dt, nSteps, beta = 0.25, gamma = 0.5) {
  const n = u0.length, Minv = lu(M);
  let u = u0.slice(), v = v0.slice();
  let a = luSolve(Minv, vadd(vadd(F(0), matvec(Cm, v), -1), matvec(K, u), -1));
  const Keff = lu(K.map((r, i) => r.map((k, j) => k + (gamma / (beta * dt)) * Cm[i][j] + M[i][j] / (beta * dt * dt))));
  const T = [0], U = [u.slice()], V = [v.slice()], A = [a.slice()];
  for (let s = 1; s <= nSteps; s++) {
    const t = s * dt, f = F(t), rhs = new Array(n);
    const mu = u.map((ui, i) => ui / (beta * dt * dt) + v[i] / (beta * dt) + (1 / (2 * beta) - 1) * a[i]);
    const cu = u.map((ui, i) => (gamma / (beta * dt)) * ui + (gamma / beta - 1) * v[i] + dt * (gamma / (2 * beta) - 1) * a[i]);
    const Mm = matvec(M, mu), Cc = matvec(Cm, cu);
    for (let i = 0; i < n; i++) rhs[i] = f[i] + Mm[i] + Cc[i];
    const un = luSolve(Keff, rhs);
    const an = un.map((x, i) => (x - u[i]) / (beta * dt * dt) - v[i] / (beta * dt) - (1 / (2 * beta) - 1) * a[i]);
    const vn = v.map((x, i) => x + dt * ((1 - gamma) * a[i] + gamma * an[i]));
    u = un; v = vn; a = an; T.push(t); U.push(u.slice()); V.push(v.slice()); A.push(a.slice());
  }
  return { t: T, u: U, v: V, a: A };
}

// ---------- optimisation ----------
export function goldenSection(f, a, b, tol = 1e-8) {
  const g = (Math.sqrt(5) - 1) / 2; let c = b - g * (b - a), d = a + g * (b - a), fc = f(c), fd = f(d);
  while (Math.abs(b - a) > tol * (1 + Math.abs(a) + Math.abs(b))) {
    if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = f(c); } else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = f(d); }
  }
  return 0.5 * (a + b);
}
/** Nelder–Mead simplex minimiser. Returns {x, f, iterations, history}. */
export function nelderMead(f, x0, { tol = 1e-9, maxIter = 2000, step = 0.1 } = {}) {
  const n = x0.length;
  let S = [x0.slice()];
  for (let i = 0; i < n; i++) { const x = x0.slice(); x[i] += x[i] !== 0 ? step * Math.abs(x[i]) : step; S.push(x); }
  let F = S.map(f); const history = [];
  for (let it = 0; it < maxIter; it++) {
    const idx = range(n + 1).sort((a, b) => F[a] - F[b]); S = idx.map((i) => S[i]); F = idx.map((i) => F[i]);
    history.push(F[0]);
    if (Math.abs(F[n] - F[0]) <= tol * (Math.abs(F[0]) + tol)) return { x: S[0], f: F[0], iterations: it, history };
    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) c[j] += S[i][j] / n;
    const xr = c.map((v, j) => v + (v - S[n][j])), fr = f(xr);
    if (fr < F[0]) { const xe = c.map((v, j) => v + 2 * (v - S[n][j])), fe = f(xe); if (fe < fr) { S[n] = xe; F[n] = fe; } else { S[n] = xr; F[n] = fr; } }
    else if (fr < F[n - 1]) { S[n] = xr; F[n] = fr; }
    else {
      const xc = fr < F[n] ? c.map((v, j) => v + 0.5 * (xr[j] - v)) : c.map((v, j) => v + 0.5 * (S[n][j] - v)), fc = f(xc);
      if (fc < Math.min(fr, F[n])) { S[n] = xc; F[n] = fc; }
      else for (let i = 1; i <= n; i++) { S[i] = S[i].map((v, j) => S[0][j] + 0.5 * (v - S[0][j])); F[i] = f(S[i]); }
    }
  }
  const k = argmin(F);
  return { x: S[k], f: F[k], iterations: maxIter, history };
}
/** Levenberg–Marquardt for nonlinear least squares. resid(p)->array. Returns {p, cost, cov, iterations}. */
export function levenbergMarquardt(resid, p0, { maxIter = 100, tol = 1e-10 } = {}) {
  let p = p0.slice(), r = resid(p), cost = dot(r, r), lam = 1e-3, J; const n = p.length;
  const jac = (pp, rr) => { const Jm = zeros(rr.length, n); for (let j = 0; j < n; j++) { const h = 1e-6 * Math.max(1e-3, Math.abs(pp[j])), q = pp.slice(); q[j] += h; const rq = resid(q); for (let i = 0; i < rr.length; i++) Jm[i][j] = (rq[i] - rr[i]) / h; } return Jm; };
  let it = 0;
  for (; it < maxIter; it++) {
    J = jac(p, r);
    const Jt = transpose(J), H = matmul(Jt, J), g = matvec(Jt, r);
    let improved = false;
    for (let k = 0; k < 20; k++) {
      const A = H.map((row, i) => row.map((v, j) => (i === j ? v * (1 + lam) + 1e-300 : v)));
      let dp; try { dp = solve(A, g.map((v) => -v)); } catch { lam *= 10; continue; }
      const pn = vadd(p, dp), rn = resid(pn), cn = dot(rn, rn);
      if (Number.isFinite(cn) && cn < cost) { const rel = (cost - cn) / (cost + 1e-300); p = pn; r = rn; cost = cn; lam = Math.max(lam / 5, 1e-12); improved = true; if (rel < tol) k = 99; break; }
      lam *= 5;
    }
    if (!improved) break;
    if (cost < 1e-28) break;
  }
  let cov = null;
  try { J = jac(p, r); const dof = Math.max(1, r.length - n); cov = mscale(inv(matmul(transpose(J), J)), cost / dof); } catch { /* singular */ }
  return { p, cost, cov, iterations: it };
}

// ---------- statistics & sampling ----------
/** Deterministic seeded uniform generator (mulberry32). */
export function rng(seed = 12345) {
  let a = seed >>> 0;
  return () => { a = (a + 0x6d2b79f5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
}
export function erf(x) {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return x >= 0 ? y : -y;
}
export const normCdf = (x) => 0.5 * (1 + erf(x / Math.SQRT2));
/** Inverse standard normal CDF (Acklam's rational approximation). */
export function normInv(p) {
  if (p <= 0) return -Infinity; if (p >= 1) return Infinity;
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239];
  const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572];
  const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
  const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
  let q, r;
  if (p < 0.02425) { q = Math.sqrt(-2 * Math.log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  if (p > 1 - 0.02425) { q = Math.sqrt(-2 * Math.log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  q = p - 0.5; r = q * q;
  return ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q) / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}
export const randn = (u) => normInv(clamp(u(), 1e-12, 1 - 1e-12));
export function gammaln(x) {
  const g = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.001208650973866179, -0.000005395239384953];
  let y = x, t = x + 5.5; t -= (x + 0.5) * Math.log(t); let s = 1.000000000190015;
  for (let j = 0; j < 6; j++) s += g[j] / ++y;
  return -t + Math.log((2.5066282746310005 * s) / x);
}
export const gamma = (x) => Math.exp(gammaln(x));
export function quantile(arr, q) {
  const a = arr.slice().sort((x, y) => x - y), pos = (a.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
  return a[lo] + (a[hi] - a[lo]) * (pos - lo);
}
export function histogram(arr, bins = 20) {
  const lo = amin(arr), hi = amax(arr), w = (hi - lo) / bins || 1, counts = new Array(bins).fill(0);
  for (const v of arr) counts[Math.min(bins - 1, Math.floor((v - lo) / w))]++;
  return { centers: range(bins, (i) => lo + (i + 0.5) * w), counts, width: w };
}
export function corr(a, b) {
  const ma = mean(a), mb = mean(b); let sab = 0, sa = 0, sb = 0;
  for (let i = 0; i < a.length; i++) { sab += (a[i] - ma) * (b[i] - mb); sa += (a[i] - ma) ** 2; sb += (b[i] - mb) ** 2; }
  return sa && sb ? sab / Math.sqrt(sa * sb) : 0;
}
/** Latin hypercube sample: n points in d dimensions on [0,1). */
export function lhs(n, d, u = rng()) {
  const out = range(n, () => new Array(d));
  for (let j = 0; j < d; j++) {
    const perm = range(n);
    for (let i = n - 1; i > 0; i--) { const k = Math.floor(u() * (i + 1)); [perm[i], perm[k]] = [perm[k], perm[i]]; }
    for (let i = 0; i < n; i++) out[i][j] = (perm[i] + u()) / n;
  }
  return out;
}
/** Map a uniform variate to a distribution {type:'normal'|'uniform'|'lognormal'|'triangular'|'weibull', ...}. */
export function sampleDist(d, p) {
  switch (d.type) {
    case 'uniform': return d.min + (d.max - d.min) * p;
    case 'lognormal': return Math.exp(Math.log(d.mean) - 0.5 * Math.log(1 + (d.sd / d.mean) ** 2) + Math.sqrt(Math.log(1 + (d.sd / d.mean) ** 2)) * normInv(p));
    case 'triangular': { const f = (d.mode - d.min) / (d.max - d.min); return p < f ? d.min + Math.sqrt(p * (d.max - d.min) * (d.mode - d.min)) : d.max - Math.sqrt((1 - p) * (d.max - d.min) * (d.max - d.mode)); }
    case 'weibull': return d.scale * (-Math.log(1 - p)) ** (1 / d.shape);
    default: return d.mean + d.sd * normInv(p);
  }
}

// ---------- signal processing ----------
/** In-place radix-2 FFT. Lengths must be a power of two. */
export function fft(re, im) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) { let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; } }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (-2 * PI) / len, wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      for (let k = 0; k < len / 2; k++) {
        const ur = re[i + k], ui = im[i + k], vr = re[i + k + len / 2] * cr - im[i + k + len / 2] * ci, vi = re[i + k + len / 2] * ci + im[i + k + len / 2] * cr;
        re[i + k] = ur + vr; im[i + k] = ui + vi; re[i + k + len / 2] = ur - vr; im[i + k + len / 2] = ui - vi;
        const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
      }
    }
  }
}
/** Single-sided amplitude spectrum of a uniformly sampled signal. Returns {f, amp}. */
export function spectrum(y, dt) {
  let n = 1; while (n * 2 <= y.length) n *= 2;
  const m = mean(y.slice(0, n)), re = y.slice(0, n).map((v, i) => (v - m) * (0.5 - 0.5 * Math.cos((2 * PI * i) / (n - 1)))), im = new Array(n).fill(0);
  fft(re, im);
  const f = [], amp = [];
  for (let k = 0; k < n / 2; k++) { f.push(k / (n * dt)); amp.push((4 * Math.hypot(re[k], im[k])) / n); }
  return { f, amp };
}

// ---------- verification helpers ----------
/**
 * Richardson extrapolation and Grid Convergence Index (Roache) from three solutions.
 * h: representative cell sizes fine->coarse [h1,h2,h3]; f: matching solutions [f1,f2,f3]. Fs is the safety factor.
 */
export function gci(h, f, Fs = 1.25) {
  const [h1, h2, h3] = h, [f1, f2, f3] = f, r21 = h2 / h1, r32 = h3 / h2, e21 = f2 - f1, e32 = f3 - f2;
  const out = { r21, r32, e21, e32, p: NaN, fExact: NaN, gciFine: NaN, gciCoarse: NaN, asymptotic: NaN, monotonic: e21 * e32 > 0 };
  if (e21 === 0 || e32 === 0) { out.p = Infinity; out.fExact = f1; out.gciFine = 0; out.gciCoarse = 0; out.asymptotic = 1; return out; }
  const s = sign(e32 / e21); let p = Math.abs(Math.log(Math.abs(e32 / e21))) / Math.log(r21);
  for (let i = 0; i < 60; i++) { const q = Math.log((r21 ** p - s) / (r32 ** p - s)); const pn = Math.abs(Math.log(Math.abs(e32 / e21)) + q) / Math.log(r21); if (!Number.isFinite(pn)) break; if (Math.abs(pn - p) < 1e-10) { p = pn; break; } p = pn; }
  out.p = p;
  out.fExact = f1 + (f1 - f2) / (r21 ** p - 1);
  out.gciFine = (Fs * Math.abs(e21 / (f1 || 1e-300))) / (r21 ** p - 1);
  out.gciCoarse = (Fs * Math.abs(e32 / (f2 || 1e-300))) / (r32 ** p - 1);
  out.asymptotic = out.gciCoarse / (r21 ** p * out.gciFine);
  return out;
}
/** Error metrics between prediction and observation arrays. */
export function errorMetrics(pred, obs) {
  const n = Math.min(pred.length, obs.length); let se = 0, ae = 0, b = 0, ape = 0, m = 0;
  for (let i = 0; i < n; i++) { const e = pred[i] - obs[i]; se += e * e; ae += Math.abs(e); b += e; if (obs[i] !== 0) { ape += Math.abs(e / obs[i]); m++; } }
  const om = mean(obs.slice(0, n)); let st = 0; for (let i = 0; i < n; i++) st += (obs[i] - om) ** 2;
  return { n, rmse: Math.sqrt(se / n), mae: ae / n, bias: b / n, mape: m ? (100 * ape) / m : NaN, r2: st ? 1 - se / st : NaN };
}
/** Helper for suite self-tests. */
export function check(name, actual, expected, tol, ref = '') {
  const err = expected === 0 ? Math.abs(actual) : Math.abs((actual - expected) / expected);
  return { name, actual, expected, tol, error: err, pass: Number.isFinite(actual) && err <= tol, ref };
}
/** Flatten a result's KPI list into {key: value}; handy inside verify() where run() is called directly. */
export const kv = (res) => Object.fromEntries((res.kpis || []).map((k) => [k.key, k.value]));
