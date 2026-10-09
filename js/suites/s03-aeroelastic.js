// Suite 3 — Aeroelasticity and Fluid–Structure Interaction.
// One unsteady strip-theory kernel (Theodorsen in the frequency domain, Wagner/Küssner indicial functions in time)
// drives a typical section, a Rayleigh–Ritz bending–torsion wing and a rotating pitch–flap blade: flutter by the
// p-k method, divergence, aileron reversal, load redistribution, limit cycles, gust response and panel flutter.
// Two coupled time-domain solvers complete the suite: an unsteady vortex lattice strongly coupled to a finite-element beam
// wing, and a two-dimensional Navier–Stokes immersed-boundary solver coupled to a spring-mounted cylinder or plate.

import * as N from '../core/numerics.js';
import { isa, G0, RHO0 } from '../core/atmosphere.js';
import { METALS } from '../data/materials.js';
import { nsCreate, nsBody, nsPredict, nsForce, nsProject, nsStore, nsRecall, nsCombine, nsDiff, nsInterior, nsVorticity, nsDivergence, nsMaxSpeed } from '../core/solvers/nsib2d.js';

const TAU = 2 * Math.PI, fin = Number.isFinite, C = N.C;

// ---- Bessel functions and the Theodorsen / Sears functions ------------------------------------
/** J0, J1, Y0, Y1 by ascending series (x ≤ 12) or Hankel asymptotic expansions. */
export function bessel(x) {
  if (x <= 12) {
    const q = (x * x) / 4, lg = Math.log(x / 2) + 0.5772156649015329;
    let t = 1, j0 = 1, j1 = 1, u = 1, y0 = 0, y1 = 1, H = 0; // u = (−q)^m/(m!(m+1)!)
    for (let m = 1; m < 60; m++) {
      t *= -q / (m * m); j0 += t; H += 1 / m; y0 -= t * H;
      u *= -q / (m * (m + 1)); j1 += u; y1 += u * (2 * H + 1 / (m + 1));
      if (Math.abs(t) < 1e-17 && Math.abs(u) < 1e-17) break;
    }
    const J1 = (x / 2) * j1;
    return { j0, j1: J1, y0: (2 / Math.PI) * (lg * j0 + y0), y1: -2 / (Math.PI * x) + (2 / Math.PI) * lg * J1 - (x / (2 * Math.PI)) * (y1 - 1) - (x / (2 * Math.PI)) * 1 };
  }
  const out = {}, s = Math.sqrt(2 / (Math.PI * x));
  for (const nu of [0, 1]) {
    const mu = 4 * nu * nu; let P = 1, Q = 0, a = 1;
    for (let k = 1; k < 30; k++) { const an = (a * (mu - (2 * k - 1) ** 2)) / (k * 8 * x); if (Math.abs(an) > Math.abs(a) && k > 2) break; a = an; if (k % 2) Q += (k % 4 === 1 ? 1 : -1) * a; else P += (k % 4 === 2 ? -1 : 1) * a; if (Math.abs(a) < 1e-16) break; }
    const chi = x - (nu / 2 + 0.25) * Math.PI;
    out['j' + nu] = s * (P * Math.cos(chi) - Q * Math.sin(chi)); out['y' + nu] = s * (P * Math.sin(chi) + Q * Math.cos(chi));
  }
  return out;
}
/** Theodorsen's function C(k) = H1⁽²⁾/(H1⁽²⁾ + i·H0⁽²⁾) as [F, G]. */
export function theodorsen(k) {
  if (!(k > 1e-9)) return [1, 0];
  if (k > 400) return [0.5, -1 / (8 * k)];
  const b = bessel(k), h1 = [b.j1, -b.y1], h0 = [b.j0, -b.y0];
  return C.div(h1, C.add(h1, C.mul([0, 1], h0)));
}
/** Sears' function S(k) = [J0 − i·J1]·C(k) + i·J1 (gust referenced to mid-chord). */
export function sears(k) { if (!(k > 1e-9)) return [1, 0]; const b = bessel(Math.min(k, 400)); return C.add(C.mul([b.j0, -b.j1], theodorsen(k)), [0, b.j1]); }
const JONES = { A: [0.165, 0.335], b: [0.0455, 0.3] }, KUSSNER = { A: [0.5, 0.5], b: [0.13, 1] };

/** Complex eigenvalue solver on typed arrays (shifted QR on the Hessenberg form), equivalent to the kernel's eig(). */
function ceig(A) {
  const n = A.length, hr = new Float64Array(n * n), hi = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { const v = A[i][j]; if (typeof v === 'number') hr[i * n + j] = v; else { hr[i * n + j] = v[0]; hi[i * n + j] = v[1]; } }
  return ceigFlat(n, hr, hi, new Float64Array(4 * n));
}
/** Eigenvalues [[re, im], ...] of the complex n×n matrix held row-major in hr/hi (destroyed); wk is scratch of length 4n. */
function ceigFlat(n, hr, hi, wk) {
  const cr = wk.subarray(0, n), ci = wk.subarray(n, 2 * n), sr = wk.subarray(2 * n, 3 * n), si = wk.subarray(3 * n, 4 * n);
  const ab = (x, y) => Math.sqrt(x * x + y * y), swap = (a, b) => { let t = hr[a]; hr[a] = hr[b]; hr[b] = t; t = hi[a]; hi[a] = hi[b]; hi[b] = t; };
  for (let m = 1; m < n - 1; m++) {
    let p = m, mx = 0;
    for (let i = m; i < n; i++) { const v = ab(hr[i * n + m - 1], hi[i * n + m - 1]); if (v > mx) { mx = v; p = i; } }
    if (mx === 0) continue;
    if (p !== m) { for (let j = 0; j < n; j++) swap(p * n + j, m * n + j); for (let j = 0; j < n; j++) swap(j * n + p, j * n + m); }
    const pr = hr[m * n + m - 1], pi = hi[m * n + m - 1], pd = pr * pr + pi * pi;
    for (let i = m + 1; i < n; i++) {
      const ar = hr[i * n + m - 1], ai = hi[i * n + m - 1]; if (ar === 0 && ai === 0) continue;
      const fr = (ar * pr + ai * pi) / pd, fi = (ai * pr - ar * pi) / pd;
      for (let j = 0; j < n; j++) { const br = hr[m * n + j], bi = hi[m * n + j]; hr[i * n + j] -= fr * br - fi * bi; hi[i * n + j] -= fr * bi + fi * br; }
      for (let j = 0; j < n; j++) { const br = hr[j * n + i], bi = hi[j * n + i]; hr[j * n + m] += fr * br - fi * bi; hi[j * n + m] += fr * bi + fi * br; }
    }
  }
  const vals = []; let nn = n, iter = 0;
  while (nn > 0) {
    if (nn === 1) { vals.push([hr[0], hi[0]]); break; }
    let l = nn - 1;
    for (; l > 0; l--) { const d = ab(hr[(l - 1) * n + l - 1], hi[(l - 1) * n + l - 1]) + ab(hr[l * n + l], hi[l * n + l]) || 1; if (ab(hr[l * n + l - 1], hi[l * n + l - 1]) < 1e-14 * d) { hr[l * n + l - 1] = 0; hi[l * n + l - 1] = 0; break; } }
    if (l === nn - 1) { vals.push([hr[l * n + l], hi[l * n + l]]); nn--; iter = 0; continue; }
    if (++iter > 500) throw new Error('ceig: QR iteration did not converge');
    const q = (nn - 2) * n + nn - 2, ar = hr[q], ai = hi[q], br = hr[q + 1], bi = hi[q + 1], er = hr[q + n], ei = hi[q + n], dr = hr[q + n + 1], di = hi[q + n + 1];
    let mr, mi;
    if (iter % 11 === 10) { mr = ab(er, ei) + (nn > 2 ? ab(hr[q - 1], hi[q - 1]) : 0); mi = 0; }
    else {
      const tr = ar + dr, ti = ai + di, xr = tr * tr - ti * ti - 4 * (ar * dr - ai * di - br * er + bi * ei), xi = 2 * tr * ti - 4 * (ar * di + ai * dr - br * ei - bi * er);
      const r = ab(xr, xi), zr = Math.sqrt((r + xr) / 2), zi = (xi < 0 ? -1 : 1) * Math.sqrt(Math.max(0, (r - xr) / 2));
      const l1r = (tr + zr) / 2, l1i = (ti + zi) / 2, l2r = (tr - zr) / 2, l2i = (ti - zi) / 2;
      if (ab(l1r - dr, l1i - di) < ab(l2r - dr, l2i - di)) { mr = l1r; mi = l1i; } else { mr = l2r; mi = l2i; }
    }
    for (let i = l; i < nn; i++) { hr[i * n + i] -= mr; hi[i * n + i] -= mi; }
    for (let k = l; k < nn - 1; k++) {
      const xr = hr[k * n + k], xi = hi[k * n + k], yr = hr[(k + 1) * n + k], yi = hi[(k + 1) * n + k], r = Math.sqrt(xr * xr + xi * xi + yr * yr + yi * yi);
      if (r === 0) { cr[k] = NaN; continue; }
      const c1 = xr / r, c2 = xi / r, s1 = yr / r, s2 = yi / r; cr[k] = c1; ci[k] = c2; sr[k] = s1; si[k] = s2;
      for (let j = k; j < nn; j++) {
        const a = k * n + j, b = a + n, ur = hr[a], ui = hi[a], vr = hr[b], vi = hi[b];
        hr[a] = c1 * ur + c2 * ui + s1 * vr + s2 * vi; hi[a] = c1 * ui - c2 * ur + s1 * vi - s2 * vr;
        hr[b] = c1 * vr - c2 * vi - s1 * ur + s2 * ui; hi[b] = c1 * vi + c2 * vr - s1 * ui - s2 * ur;
      }
    }
    for (let k = l; k < nn - 1; k++) {
      if (Number.isNaN(cr[k])) continue;
      const c1 = cr[k], c2 = ci[k], s1 = sr[k], s2 = si[k], top = Math.min(k + 2, nn - 1);
      for (let i = l; i <= top; i++) {
        const a = i * n + k, ur = hr[a], ui = hi[a], vr = hr[a + 1], vi = hi[a + 1];
        hr[a] = ur * c1 - ui * c2 + vr * s1 - vi * s2; hi[a] = ur * c2 + ui * c1 + vr * s2 + vi * s1;
        hr[a + 1] = vr * c1 + vi * c2 - ur * s1 - ui * s2; hi[a + 1] = vi * c1 - vr * c2 - ui * s1 + ur * s2;
      }
    }
    for (let i = l; i < nn; i++) { hr[i * n + i] += mr; hi[i * n + i] += mi; }
  }
  return vals;
}

// ---- unsteady strip-theory kernel -------------------------------------------------------------
// A model is { n, M, C, K, strips }. Each strip: { b, a, dy, cl, U, ph[n], pa[n], sw[n] } with plunge h (down) = Σ ph·q,
// pitch α (nose-up) = Σ pa·q about the elastic axis at a·b aft of mid-chord, and sw the sweep incidence term tanΛ·dh/dy.
/** Per-strip vectors: downwash at ¾-chord w = Wd·q̇ + Wk·q, and the circulatory force column Fc (lift at ¼-chord). */
function stripVecs(m, rho) {
  return m.strips.map((s) => ({
    Wd: s.ph.map((h, j) => h + s.b * (0.5 - s.a) * s.pa[j]), Wk: s.pa.map((p, j) => s.U * (p + (s.sw ? s.sw[j] : 0))),
    Fc: s.ph.map((h, j) => s.cl * rho * s.U * s.b * s.dy * (h - s.b * (s.a + 0.5) * s.pa[j])), s,
  }));
}
/** Non-circulatory (apparent-mass) matrices: real added mass and damping. */
function apparent(m, rho) {
  const n = m.n, Ma = N.zeros(n), Ca = N.zeros(n);
  for (const s of m.strips) {
    const g = Math.PI * rho * s.b * s.b * s.dy;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      Ma[i][j] += g * (s.ph[i] * s.ph[j] - s.b * s.a * (s.ph[i] * s.pa[j] + s.pa[i] * s.ph[j]) + s.b * s.b * (0.125 + s.a * s.a) * s.pa[i] * s.pa[j]);
      Ca[i][j] += g * s.U * (s.ph[i] * s.pa[j] + s.b * (0.5 - s.a) * s.pa[i] * s.pa[j]);
    }
  }
  return { Ma, Ca };
}
/**
 * Complex aerodynamic + structural damping and stiffness at the oscillation frequency w (quasi → C(k) = 1), flat row-major.
 * The four matrices live in buffers owned by pre and are overwritten by the next call.
 */
function aeroCK(m, pre, w, quasi) {
  const nn = m.n * m.n, { Cr, Ci, Kr, Ki, C0, K0, st, D, Kk } = pre;
  Cr.set(C0); Kr.set(K0); Ci.fill(0); Ki.fill(0);
  for (let s = 0, o = 0; s < st.length; s++, o += nn) {
    let F = 1, G = 0; if (!quasi) { const ck = theodorsen((w * st[s].b) / Math.max(st[s].U, 1e-9)); F = ck[0]; G = ck[1]; }
    for (let q = 0; q < nn; q++) { const d = D[o + q], k = Kk[o + q]; Cr[q] += F * d; Ci[q] += G * d; Kr[q] += F * k; Ki[q] += G * k; }
  }
  return pre;
}
/** Roots of [s²(M+Ma) + s(C+Ca+ΣC(k)·Fc·Wdᵀ) + K + ΣC(k)·Fc·Wkᵀ] with the aerodynamics evaluated at the oscillation frequency w. */
function pRoots(m, w, pre, quasi) {
  const n = m.n, n2 = 2 * n, a = aeroCK(m, pre, w, quasi), hr = pre.hr, hi = pre.hi, Mi = pre.Mif;
  hr.fill(0); hi.fill(0);
  for (let i = 0; i < n; i++) {
    hr[i * n2 + i + n] = 1;
    for (let j = 0; j < n; j++) { let kr = 0, ki = 0, c1 = 0, c2 = 0; for (let l = 0; l < n; l++) { const mi = Mi[i * n + l], q = l * n + j; kr += mi * a.Kr[q]; ki += mi * a.Ki[q]; c1 += mi * a.Cr[q]; c2 += mi * a.Ci[q]; } const r = (i + n) * n2 + j; hr[r] = -kr; hi[r] = -ki; hr[r + n] = -c1; hi[r + n] = -c2; }
  }
  return ceigFlat(n2, hr, hi, pre.wk);
}
const cdist = (a, b) => Math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2), cmod = (a) => Math.sqrt(a[0] * a[0] + a[1] * a[1]);
/** Everything that does not depend on frequency at one condition: apparent mass, rank-one strip matrices (flat) and work buffers. */
function prep(m, rho) {
  const n = m.n, nn = n * n, st = m.strips, ns = st.length, Mt = m.M.map((r) => r.slice()), D = new Float64Array(ns * nn), Kk = new Float64Array(ns * nn), C0 = new Float64Array(nn), K0 = new Float64Array(nn), wd = new Float64Array(n), wk = new Float64Array(n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { C0[i * n + j] = m.C[i][j]; K0[i * n + j] = m.K[i][j]; }
  for (let s = 0; s < ns; s++) {
    const t = st[s], b = t.b, a = t.a, ph = t.ph, pa = t.pa, sw = t.sw, g = Math.PI * rho * b * b * t.dy, fc = t.cl * rho * t.U * b * t.dy, o = s * nn;
    for (let j = 0; j < n; j++) { wd[j] = ph[j] + b * (0.5 - a) * pa[j]; wk[j] = t.U * (pa[j] + (sw ? sw[j] : 0)); }
    for (let i = 0; i < n; i++) {
      const f = fc * (ph[i] - b * (a + 0.5) * pa[i]), Mr = Mt[i];
      for (let j = 0; j < n; j++) {
        Mr[j] += g * (ph[i] * ph[j] - b * a * (ph[i] * pa[j] + pa[i] * ph[j]) + b * b * (0.125 + a * a) * pa[i] * pa[j]);
        C0[i * n + j] += g * t.U * (ph[i] * pa[j] + b * (0.5 - a) * pa[i] * pa[j]);
        D[o + i * n + j] = f * wd[j]; Kk[o + i * n + j] = f * wk[j];
      }
    }
  }
  const Mi = N.inv(Mt), Mif = new Float64Array(nn); for (let i = 0; i < n; i++) Mif.set(Mi[i], i * n);
  return { st, Mt, Mi, Mif, D, Kk, C0, K0, Cr: new Float64Array(nn), Ci: new Float64Array(nn), Kr: new Float64Array(nn), Ki: new Float64Array(nn), re: new Float64Array(nn), im: new Float64Array(nn), hr: new Float64Array(4 * nn), hi: new Float64Array(4 * nn), wk: new Float64Array(8 * n) };
}
/** Determinant of the complex n×n matrix held in re/im (row-major, overwritten) by Gaussian elimination with pivoting, as [re, im]. */
function cdet(n, re, im) {
  let pr = 1, pi = 0;
  for (let k = 0; k < n; k++) {
    let p = k, mx = -1; for (let i = k; i < n; i++) { const v = re[i * n + k] ** 2 + im[i * n + k] ** 2; if (v > mx) { mx = v; p = i; } }
    if (!(mx > 0)) return [0, 0];
    if (p !== k) { for (let j = 0; j < n; j++) { let t = re[p * n + j]; re[p * n + j] = re[k * n + j]; re[k * n + j] = t; t = im[p * n + j]; im[p * n + j] = im[k * n + j]; im[k * n + j] = t; } pr = -pr; pi = -pi; }
    const ar = re[k * n + k], ai = im[k * n + k], tr = pr * ar - pi * ai; pi = pr * ai + pi * ar; pr = tr;
    for (let i = k + 1; i < n; i++) { const fr = (re[i * n + k] * ar + im[i * n + k] * ai) / mx, fi = (im[i * n + k] * ar - re[i * n + k] * ai) / mx; if (fr === 0 && fi === 0) continue; for (let j = k + 1; j < n; j++) { const br = re[k * n + j], bi = im[k * n + j]; re[i * n + j] -= fr * br - fi * bi; im[i * n + j] -= fr * bi + fi * br; } }
  }
  return [pr, pi];
}
/** Flutter determinant det[s²(M+Ma) + s·C + K] at the complex s = [σ, ω] for aerodynamic matrices a = aeroCK(...). */
function fdet(m, pre, a, s) {
  const n = m.n, re = pre.re, im = pre.im, s2r = s[0] * s[0] - s[1] * s[1], s2i = 2 * s[0] * s[1];
  for (let i = 0; i < n; i++) for (let j = 0, q = i * n; j < n; j++, q++) { const mt = pre.Mt[i][j]; re[q] = s2r * mt + s[0] * a.Cr[q] - s[1] * a.Ci[q] + a.Kr[q]; im[q] = s2i * mt + s[0] * a.Ci[q] + s[1] * a.Cr[q] + a.Ki[q]; }
  return cdet(n, re, im);
}
/**
 * p-k iteration from the root guess [σ, ω]. One full eigen-solution at the guessed frequency picks the nearest root; the
 * aerodynamics are then re-evaluated at that root's frequency and the root is polished on the flutter determinant by the
 * secant method (the determinant is a polynomial in s for frozen aerodynamics), until the frequency stops moving.
 */
function pkIter(m, pre, guess, quasi) {
  let lam = null, bd = Infinity, w = Math.abs(guess[1]);
  for (const x of pRoots(m, w, pre, quasi)) { if (x[1] < -1e-9 * (1 + cmod(x))) continue; const d = cdist(x, guess); if (d < bd) { bd = d; lam = x; } }
  if (!lam) return guess;
  lam = [lam[0], Math.max(lam[1], 0)];
  for (let it = 0; it < 20 && !quasi; it++) {
    if (Math.abs(lam[1] - w) < 1e-8 * (1 + lam[1])) break;
    w = it < 8 ? lam[1] : 0.5 * (w + lam[1]); // plain fixed point first; under-relaxed if it has not settled (heavily damped roots)
    const a = aeroCK(m, pre, w, false), sc = 1e-5 * (1 + cmod(lam)); let s0 = lam, f0 = fdet(m, pre, a, s0), s1 = [lam[0] + sc, lam[1] + sc], f1 = fdet(m, pre, a, s1), ok = false;
    for (let k = 0; k < 12; k++) {
      const dr = f1[0] - f0[0], di = f1[1] - f0[1], dd = dr * dr + di * di; if (!(dd > 0)) { ok = cmod(f1) === 0; break; }
      const qr = (f1[0] * dr + f1[1] * di) / dd, qi = (f1[1] * dr - f1[0] * di) / dd, hr = s1[0] - s0[0], hi = s1[1] - s0[1], s2 = [s1[0] - (qr * hr - qi * hi), s1[1] - (qr * hi + qi * hr)];
      s0 = s1; f0 = f1; s1 = s2; f1 = fdet(m, pre, a, s1);
      if (!fin(s1[0]) || !fin(s1[1])) break;
      if (cdist(s1, s0) < 1e-11 * (1 + cmod(s1))) { ok = true; break; }
    }
    // the polish must stay on the root it started from; otherwise fall back to the eigen-solution at this frequency
    if (ok && s1[1] > -1e-9 * (1 + cmod(s1)) && cdist(s1, lam) < 0.2 * (cmod(lam) + 1e-12)) lam = [s1[0], Math.max(s1[1], 0)];
    else { let best = null; bd = Infinity; for (const x of pRoots(m, w, pre, false)) { if (x[1] < -1e-9 * (1 + cmod(x))) continue; const d = cdist(x, lam); if (d < bd) { bd = d; best = x; } } if (!best) break; lam = [best[0], Math.max(best[1], 0)]; }
  }
  return lam;
}
/** p-k iteration for one mode at one condition, started from the root guess [σ, ω]. */
const pkPoint = (m, rho, guess, quasi) => pkIter(m, prep(m, rho), guess, quasi);
/**
 * p-k roots of every tracked mode at one condition, continued from the seeds (the roots at the previous condition). Nearest-root
 * continuation alone lets two trackers settle on the same root where two frequencies approach each other, and the other root —
 * often the one that flutters — is then lost for the rest of the sweep. Whenever two trackers coincide, the one that moved
 * farther is restarted from every other root of the eigen-solution until it converges on a root nobody else holds.
 */
function pkModes(m, rho, seeds, quasi) {
  const pre = prep(m, rho), out = seeds.map((g) => pkIter(m, pre, g, quasi)), same = (a, b) => cdist(a, b) < 1e-5 * (1 + cmod(a));
  for (let k = 1; k < out.length; k++) for (let j = 0; j < k; j++) {
    if (!same(out[k], out[j])) continue;
    const mv = cdist(out[k], seeds[k]) >= cdist(out[j], seeds[j]) ? k : j, held = (r) => out.some((o, q) => q !== mv && same(o, r));
    const cands = [...pRoots(m, Math.abs(out[mv][1]), pre, quasi), ...pRoots(m, Math.abs(seeds[mv][1]), pre, quasi)].filter((x) => x[1] >= -1e-9 * (1 + cmod(x))).sort((a, b) => cdist(a, seeds[mv]) - cdist(b, seeds[mv]));
    for (const c of cands) { const r = pkIter(m, pre, c, quasi); if (!held(r)) { out[mv] = r; break; } }
  }
  return out;
}
/**
 * Number of unstable roots of the flutter determinant D(s) = det[s²(M+Ma) + s·C(k) + K(k)] from its phase along s = iω
 * (argument principle): the phase rises by (n − Z)·π between ω = 0 and ∞ when Z roots lie in the right half-plane. It needs
 * neither root tracking nor the p-k iteration, so it is an independent check on the p-k sweep. Returns { Z, osc } with osc
 * the number of unstable roots beyond one static (divergence) root, or NaN values when the phase cannot be resolved.
 */
function nyquistCount(m, rho, quasi) {
  const n = m.n, pre = prep(m, rho), Mt = pre.Mt, k0 = Float64Array.from(aeroCK(m, pre, 0, true).Kr);
  let wHi = 0, wLo = Infinity; for (let i = 0; i < n; i++) { const kk = (Math.abs(m.K[i][i]) + Math.abs(k0[i * n + i])) / Mt[i][i], ks = Math.abs(m.K[i][i]) / Mt[i][i]; if (kk > wHi) wHi = kk; if (ks > 0 && ks < wLo) wLo = ks; }
  wHi = 100 * Math.sqrt(wHi); wLo = fin(wLo) ? 1e-4 * Math.sqrt(wLo) : 1e-8 * wHi;
  if (!(wHi > 0) || !(wLo > 0)) return { Z: NaN, osc: NaN };
  const arg = (w) => { const d = fdet(m, pre, aeroCK(m, pre, w, quasi), [0, w]); return d[0] === 0 && d[1] === 0 ? NaN : Math.atan2(d[1], d[0]); };
  const wrap = (d) => d - TAU * Math.round(d / TAU); let ok = true;
  const seg = (w0, a0, w1, a1, depth) => { const d = wrap(a1 - a0); if (Math.abs(d) < 0.6) return d; if (depth > 44) { if (Math.abs(d) > 2.4) ok = false; return d; } const wm = Math.sqrt(w0 * w1), am = arg(wm); return seg(w0, a0, wm, am, depth + 1) + seg(wm, am, w1, a1, depth + 1); };
  const nP = Math.ceil(12 * Math.log10(wHi / wLo)), a0 = arg(wLo); let tot = 0, wp = wLo, ap = a0;
  for (let k = 1; k <= nP; k++) { const w = wLo * (wHi / wLo) ** (k / nP), a = arg(w); tot += seg(wp, ap, w, a, 0); wp = w; ap = a; }
  const z = n - tot / Math.PI, Z = Math.round(z);
  if (!ok || !fin(z) || Math.abs(z - Z) > 0.3) return { Z: NaN, osc: NaN };
  return { Z, osc: Z - (Math.abs(a0) > Math.PI / 2 ? 1 : 0) };
}
/** Static aeroelastic stiffness per unit (reference speed)²: K̂ = Σ Fc·Wkᵀ evaluated by the caller's model. */
function aeroStiff(m, rho) { const n = m.n, Kh = N.zeros(n); for (const v of stripVecs(m, rho)) for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Kh[i][j] += v.Fc[i] * v.Wk[j]; return Kh; }
/**
 * p-k sweep over a parameter grid P (airspeed or rotor speed). mk(P) builds the model. Every mode in w0 is continued through
 * the sweep with pkModes; an interval is halved (up to 8 sub-steps) where a lightly damped root moves quickly or two roots
 * approach each other, so that branches are not confused near a frequency coalescence. Flutter is the lowest P at which ANY
 * branch acquires positive damping (refined by bisection), cross-checked with the determinant phase count; divergence is the
 * first zero of det(K + K_aero). Options: quasi (C(k) = 1), ptol (relative tolerance on the flutter value of P), check (false skips the phase count). Returns the root tracks on the
 * grid, { flutter: { P, w, mode }, div, check }.
 */
function pkSweep(mk, rho, Ps, w0, { quasi = false, ptol = 1e-7, check: doCheck = true } = {}) {
  const nm = w0.length, tracks = w0.map(() => []), tol = (r) => 1e-9 * (1 + cmod(r)), unst = (s) => s.findIndex((r) => r[0] > tol(r) && r[1] > 1e-6);
  let flutter = null;
  const light = (r) => r[0] > -0.6 * cmod(r) && r[1] > 1e-6;
  const rough = (s0, s1) => {
    for (let k = 0; k < nm; k++) { if (!light(s1[k])) continue; if (cdist(s0[k], s1[k]) > 0.25 * cmod(s1[k])) return true; for (let j = 0; j < k; j++) if (light(s1[j]) && cdist(s1[k], s1[j]) < 0.2 * cmod(s1[k]) && cdist(s0[k], s1[k]) > 0.04 * cmod(s1[k])) return true; }
    return false;
  };
  const cross = (P0, s0, P1, s1) => {
    // bracketed false position (Illinois) on the damping of the crossing branch, every branch being continued from the stable side
    let lo = P0, hi = P1, g = s0, k = unst(s1), root = s1[k], flo = g[k][0] - tol(g[k]), fhi = root[0] - tol(root);
    for (let it = 0; it < 60 && hi - lo > ptol * hi && fhi > 1e-7 * cmod(root); it++) {
      let mid = it % 3 !== 2 && flo < 0 && fhi > 0 ? lo - (flo * (hi - lo)) / (fhi - flo) : 0.5 * (lo + hi); mid = N.clamp(mid, lo + 0.02 * (hi - lo), hi - 0.02 * (hi - lo)); // every third step is a plain bisection
      const sm = pkModes(mk(mid), rho, g, quasi), km = unst(sm);
      if (km >= 0) { hi = mid; root = sm[km]; fhi = root[0] - tol(root); if (km !== k) { k = km; flo = g[k][0] - tol(g[k]); } else flo *= 0.5; }
      else { lo = mid; g = sm; flo = g[k][0] - tol(g[k]); fhi *= 0.5; }
    }
    flutter = { P: hi, w: root[1], mode: k };
  };
  const advance = (P0, s0, P1, depth) => {
    const s1 = pkModes(mk(P1), rho, s0, quasi);
    if (depth < 3 && rough(s0, s1)) { const Pm = 0.5 * (P0 + P1); return advance(Pm, advance(P0, s0, Pm, depth + 1), P1, depth + 1); }
    if (!flutter && unst(s1) >= 0) { if (unst(s0) < 0) cross(P0, s0, P1, s1); }
    return s1;
  };
  let state = pkModes(mk(Ps[0]), rho, w0.map((w) => [0, w]), quasi); state.forEach((l, k) => tracks[k].push(l));
  if (unst(state) >= 0) { const k = unst(state); flutter = { P: Ps[0], w: state[k][1], mode: k, atStart: true }; }
  for (let j = 1; j < Ps.length; j++) { state = advance(Ps[j - 1], state, Ps[j], 0); state.forEach((l, k) => tracks[k].push(l)); }
  // independent check: the determinant phase count must show no oscillatory instability just below the p-k flutter speed
  // (or at the top of the range when the sweep found none); if it does, a branch was missed and the lower boundary is located with it
  const check = { ok: true, Z: NaN }, top = Ps[Ps.length - 1], Pc = flutter ? 0.98 * flutter.P : top, oscAt = (P) => nyquistCount(mk(P), rho, quasi).osc;
  if (doCheck && !(flutter && flutter.atStart) && Pc > Ps[0]) {
    const z = oscAt(Pc); check.Z = z;
    if (z >= 2) {
      check.ok = false; let j = 0; while (j < Ps.length - 1 && Ps[j] < Pc && !(oscAt(Ps[j]) >= 2)) j++;
      let hi = Math.min(Ps[j], Pc), lo = j > 0 ? Ps[j - 1] : 0.5 * Ps[0];
      if (j === 0 || !(oscAt(lo) >= 2)) for (let it = 0; it < 14; it++) { const mid = 0.5 * (lo + hi); if (oscAt(mid) >= 2) hi = mid; else lo = mid; }
      const sm = pkModes(mk(hi), rho, tracks.map((t) => t[Math.max(0, j - 1)]), quasi), k = unst(sm);
      flutter = { P: hi, w: k >= 0 ? sm[k][1] : NaN, mode: k, byCount: true };
    }
  }
  const det = (P) => { const m = mk(P); return N.det(N.madd(m.K, aeroStiff(m, rho))) / N.det(m.K); };
  let div = NaN, Pp = 1e-3 * Ps[0], d0 = det(Pp);
  for (let j = 0; j < Ps.length; j++) { const d1 = det(Ps[j]); if (d0 > 0 && d1 <= 0) { div = N.brent(det, Pp, Ps[j], 1e-7 * Ps[j]); break; } d0 = d1; Pp = Ps[j]; }
  return { tracks, flutter, div, check };
}
/** Time-domain state matrix with two Jones lag states per strip: x = [q, q̇, y1, y2]. Also returns the layout for forcing. */
function stateSpace(m, rho) {
  const n = m.n, ap = apparent(m, rho), sv = stripVecs(m, rho), Mi = N.inv(N.madd(m.M, ap.Ma)), ns = sv.length, nx = 2 * n + 2 * ns, A = N.zeros(nx), a0 = 1 - JONES.A[0] - JONES.A[1];
  const Kq = m.K.map((r) => r.slice()), Cq = N.madd(m.C, ap.Ca), Fy = N.zeros(n, 2 * ns);
  sv.forEach((v, s) => {
    for (let i = 0; i < n; i++) { for (let j = 0; j < n; j++) { Kq[i][j] += a0 * v.Fc[i] * v.Wk[j]; Cq[i][j] += a0 * v.Fc[i] * v.Wd[j]; } Fy[i][2 * s] = Fy[i][2 * s + 1] = v.Fc[i]; }
    for (let l = 0; l < 2; l++) { const r = 2 * n + 2 * s + l, be = (JONES.b[l] * v.s.U) / v.s.b; A[r][r] = -be; for (let j = 0; j < n; j++) { A[r][j] = be * JONES.A[l] * v.Wk[j]; A[r][n + j] = be * JONES.A[l] * v.Wd[j]; } }
  });
  const MK = N.matmul(Mi, Kq), MC = N.matmul(Mi, Cq), MF = N.matmul(Mi, Fy);
  for (let i = 0; i < n; i++) { A[i][n + i] = 1; for (let j = 0; j < n; j++) { A[n + i][j] = -MK[i][j]; A[n + i][n + j] = -MC[i][j]; } for (let j = 0; j < 2 * ns; j++) A[n + i][2 * n + j] = -MF[i][j]; }
  return { A, Mi, n, nx };
}
const growth = (m, rho) => N.amax(ceig(stateSpace(m, rho).A).map((l) => l[0]));

// ---- model builders -----------------------------------------------------------------------------
/** Typical section per unit span. i: b, a, x_alpha, r_alpha, mu·πρb² = m, wh, wa, zeta, cl. */
function sectionModel(i, m, U) {
  const S = m * i.b * i.x_alpha, Ia = m * (i.b * i.r_alpha) ** 2;
  return { n: 2, M: [[m, S], [S, Ia]], C: [[2 * i.zeta * m * i.wh, 0], [0, 2 * i.zeta * Ia * i.wa]], K: [[m * i.wh ** 2, 0], [0, Ia * i.wa ** 2]], strips: [{ b: i.b, a: i.a, dy: 1, cl: i.cl, U, ph: [1, 0], pa: [0, 1] }] };
}
/** Rayleigh–Ritz cantilever wing: bending η², η³…, torsion η, η²… along the swept elastic axis; strip aerodynamics normal to it. */
function wingRitz(i, nb, nt, ns) {
  const L = i.semi_span, n = nb + nt, sw = N.rad(i.sweep_deg || 0), cs = Math.cos(sw), tn = Math.tan(sw), a = 2 * i.x_ea - 1, xa = 2 * (i.x_cg - i.x_ea);
  const ch = (e) => i.c_root * (1 - (1 - i.taper) * e) * cs, st = (e) => (ch(e) / (i.c_root * cs)) ** i.stiff_exp, m0 = (3 * i.m_semi) / (L * (1 + i.taper + i.taper ** 2)), mm = (e) => m0 * (ch(e) / (i.c_root * cs)) ** 2;
  const ph = (k, e) => (k < nb ? e ** (k + 2) : 0), dph = (k, e) => (k < nb ? ((k + 2) * e ** (k + 1)) / L : 0), d2 = (k, e) => (k < nb ? ((k + 2) * (k + 1) * e ** k) / (L * L) : 0);
  const pa = (k, e) => (k >= nb ? e ** (k - nb + 1) : 0), dpa = (k, e) => (k >= nb ? ((k - nb + 1) * e ** (k - nb)) / L : 0);
  const M = N.zeros(n), K = N.zeros(n), I = (f) => L * N.simpson(f, 0, 1, 80);
  for (let p = 0; p < n; p++) for (let q = p; q < n; q++) {
    M[p][q] = M[q][p] = I((e) => { const b = ch(e) / 2, m = mm(e); return m * (ph(p, e) * ph(q, e) + xa * b * (ph(p, e) * pa(q, e) + pa(p, e) * ph(q, e)) + (2 * i.r_alpha * b) ** 2 * pa(p, e) * pa(q, e)); });
    K[p][q] = K[q][p] = I((e) => st(e) * (i.EI_root * d2(p, e) * d2(q, e) + i.GJ_root * dpa(p, e) * dpa(q, e)));
  }
  const eg = N.eigGenSym(N.mscale(K, 1 / K[0][0]), N.mscale(M, 1 / M[0][0])), sc = K[0][0] / M[0][0], wn = eg.values.map((v) => Math.sqrt(Math.max(v * sc, 0)));
  // modal damping C = M·Φ·diag(2ζω)·Φᵀ·M with M-orthonormal Φ
  const Cm = N.zeros(n);
  if (i.zeta > 0) eg.vectors.forEach((v, k) => { const u = N.matvec(M, v).map((x) => x / Math.sqrt(M[0][0])); for (let p = 0; p < n; p++) for (let q = 0; q < n; q++) Cm[p][q] += 2 * i.zeta * wn[k] * u[p] * u[q]; });
  const geo = N.range(ns, (s) => { const e = (s + 0.5) / ns; return { e, y: e * L * cs, b: ch(e) / 2, a, dy: L / ns, cl: i.cl, ph: N.range(n, (k) => ph(k, e)), pa: N.range(n, (k) => pa(k, e)), sw: N.range(n, (k) => tn * dph(k, e)) }; });
  // torsion-dominated modes: classify each normal mode by its torsional share of kinetic energy
  const tors = eg.vectors.map((v) => { let t = 0, tot = 0; for (let p = 0; p < n; p++) for (let q = 0; q < n; q++) { const e = v[p] * M[p][q] * v[q]; tot += e; if (p >= nb && q >= nb) t += e; } return N.clamp(t / tot, 0, 1); });
  return { n, M, C: Cm, K, wn, tors, cs, at: (V) => ({ n, M, C: Cm, K, strips: geo.map((g) => ({ ...g, U: V * cs })) }) };
}
/** Rigid pitch–flap blade in the rotating frame: q = [β (flap up), θ (pitch nose-up)]. */
function bladeModel(i, Om, ns = 10) {
  const R = i.R, e = i.e_hinge * R, l = R - e, mp = i.m_blade / l, b = i.chord / 2, a = 2 * i.x_pa - 1, xI = (i.x_cg - i.x_pa) * i.chord;
  const Ib = (mp * l ** 3) / 3, Ix = (xI * mp * l * l) / 2, If = mp * l * (i.r_alpha * i.chord) ** 2, nu2 = 1 + (1.5 * e) / l + (i.f_flap_nr * TAU) ** 2 / (Om * Om || 1e-12), wt = i.f_tors * TAU;
  const M = [[Ib, -Ix], [-Ix, If]], K = [[Ib * nu2 * Om * Om, -Ix * Om * Om], [-Ix * Om * Om, If * (Om * Om + wt * wt)]];
  const strips = N.range(ns, (s) => { const r = e + ((s + 0.5) / ns) * l; return { b, a, dy: l / ns, cl: i.cl, U: Math.max(Om * r, 1e-6), ph: [-(r - e), 0], pa: [0, 1] }; });
  return { n: 2, M, C: [[0, 0], [0, 2 * i.zeta * If * Math.sqrt(Om * Om + wt * wt)]], K, strips, Ib, If, wn: [Math.sqrt(nu2) * Om, Math.sqrt(Om * Om + wt * wt)] };
}

// ---- defaults from the shared case --------------------------------------------------------------
const helmbold = (AR, sw = 0) => (AR > 0 ? (TAU * AR) / (2 + Math.sqrt(4 + (AR * AR) / Math.cos(sw) ** 2)) : TAU);
function wingDefaults(c, up, d) {
  if (!(c.wing.S_m2 > 0)) return {};
  const mat = METALS[c.struct.material] || METALS['Al 2024-T3'], cr = d.c_root, w = c.struct.box_chord_frac * cr, h = c.struct.box_height_frac * c.wing.tc * cr, ts = c.struct.t_skin_mm / 1e3, tw = c.struct.t_spar_mm / 1e3;
  const EI = mat.E * (2 * 1.6 * w * ts * (h / 2) ** 2 + (2 * tw * h ** 3) / 12), GJ = (mat.G * 4 * (w * h) ** 2) / ((2 * w) / ts + (2 * h) / tw), sw = N.rad(c.wing.sweep_deg);
  const Vd = up.performance?.V_d_eas ?? 1.25 * (c.aero.Vmo_ms || c.flight.V_ms);
  return {
    semi_span: c.wing.b_m / 2 / Math.cos(sw), c_root: cr, taper: c.wing.taper, sweep_deg: c.wing.sweep_deg, EI_root: up.fea?.EI_root_Nm2 ?? EI, GJ_root: up.fea?.GJ_root_Nm2 ?? GJ,
    m_semi: (up.fea?.wing_struct_mass_kg > 0 ? up.fea.wing_struct_mass_kg / 2 + 0.02 * c.mass.mtow_kg : 0.06 * c.mass.mtow_kg) + 0.3 * c.mass.fuel_kg,
    cl: up.cfd?.CLa_per_rad ?? helmbold(d.AR, sw), zeta: c.struct.zeta, alt_m: c.atm.alt_m, V: c.flight.V_ms, V_d_eas: Vd, mass_kg: c.mass.mtow_kg, aileron_frac: Math.min(0.35, 4 * c.controls.Sa_S),
  };
}
const WING_INPUTS = [
  { key: 'semi_span', label: 'Structural semi-span', unit: 'm', default: 17, min: 0.1, group: 'Geometry', help: 'Root to tip along the elastic axis' },
  { key: 'c_root', label: 'Root chord (streamwise)', unit: 'm', default: 5.8, min: 0.02, group: 'Geometry' },
  { key: 'taper', label: 'Taper ratio', unit: '-', default: 0.3, min: 0.05, max: 1, group: 'Geometry' },
  { key: 'sweep_deg', label: 'Elastic-axis sweep (aft +)', unit: 'deg', default: 0, min: -45, max: 60, group: 'Geometry' },
  { key: 'x_ea', label: 'Elastic axis position / chord', unit: '-', default: 0.4, min: 0.1, max: 0.7, group: 'Section', help: '0.35–0.45 for a two-spar box' },
  { key: 'x_cg', label: 'Section mass centre / chord', unit: '-', default: 0.45, min: 0.1, max: 0.8, group: 'Section', help: 'Aft of the elastic axis is destabilising; mass balance moves it forward' },
  { key: 'r_alpha', label: 'Radius of gyration about the elastic axis / chord', unit: '-', default: 0.25, min: 0.1, max: 0.5, group: 'Section' },
  { key: 'EI_root', label: 'Root bending stiffness EI', unit: 'N·m²', default: 1.5e8, min: 1e-3, group: 'Structure', help: 'From the structures suite when available' },
  { key: 'GJ_root', label: 'Root torsional stiffness GJ', unit: 'N·m²', default: 1.2e8, min: 1e-3, group: 'Structure' },
  { key: 'stiff_exp', label: 'Stiffness taper exponent', unit: '-', default: 3, min: 0, max: 5, group: 'Structure', help: 'EI, GJ ∝ chord^p' },
  { key: 'm_semi', label: 'Semi-wing mass including fuel', unit: 'kg', default: 10000, min: 1e-4, group: 'Structure', help: 'Distributed ∝ chord²' },
  { key: 'zeta', label: 'Structural damping ratio', unit: '-', default: 0.02, min: 0, max: 0.2, group: 'Structure', help: '0.01–0.03 metallic; flutter clearance is often shown with zero and with g = 0.03' },
  { key: 'cl', label: 'Strip lift-curve slope', unit: '1/rad', default: 5, min: 1, max: 7, group: 'Aerodynamics', help: 'Wing lift-curve slope applied to every strip (finite-span correction of 2π)' },
  { key: 'alt_m', label: 'Altitude', unit: 'm', default: 0, min: -500, max: 25000, group: 'Flight' },
  { key: 'V_d_eas', label: 'Design dive speed VD (EAS)', unit: 'm/s', default: 225, min: 1, group: 'Flight', help: 'From the performance suite when available' },
];
const wingPick = (o) => Object.fromEntries(WING_INPUTS.map((f) => [f.key, o[f.key]]));
const hasWing = (c) => (c.wing.S_m2 > 0 ? true : 'This analysis needs a wing; use the rotor-blade analyses for rotorcraft and multirotors.');
const hasRotor = (c) => (c.rotor.R_m > 0 || c.prop.prop_dia_m > 0 ? true : 'This analysis needs a rotor or propeller.');
const margin = (VfEas, Vd) => VfEas / (1.15 * Vd) - 1;
/** Plain-language notes on how the p-k sweep located its flutter point. */
function pkNotes(r, warnings, what = 'speed') {
  if (r.flutter?.byCount) warnings.push(`The determinant phase count found an oscillatory instability below the lowest crossing of the tracked p-k branches; the flutter ${what} reported is the one located by the count. Refine the ${what} grid and inspect the damping diagram.`);
  if (r.flutter?.atStart) warnings.push(`A branch is already unstable at the lowest ${what} of the sweep: the flutter ${what} lies at or below the value reported.`);
}
const PK_MODEL = 'p-k roots of every mode continued together (no two trackers on one root, grid refined where branches approach), flutter taken as the first branch to lose damping and cross-checked by the phase count of the flutter determinant (argument principle)';
const MARGIN_NOTE = 'Flutter EAS / (1.15 · VD) − 1: positive means flutter-free to 1.15 VD';

// ---- 1. typical-section flutter (p-k) and limit cycles ------------------------------------------
const SEC_INPUTS = [
  { key: 'b', label: 'Semi-chord', unit: 'm', default: 1, min: 0.005, group: 'Section', help: 'Half the chord at the representative station (about 75% semi-span)' },
  { key: 'a', label: 'Elastic axis aft of mid-chord / semi-chord', unit: '-', default: -0.2, min: -0.9, max: 0.9, group: 'Section', help: '−0.2 is 40% chord' },
  { key: 'x_alpha', label: 'Mass centre aft of the elastic axis / semi-chord', unit: '-', default: 0.1, min: -0.5, max: 0.8, group: 'Section' },
  { key: 'r_alpha', label: 'Radius of gyration about the elastic axis / semi-chord', unit: '-', default: 0.5, min: 0.1, max: 1.5, group: 'Section' },
  { key: 'm', label: 'Mass per unit span', unit: 'kg/m', default: 100, min: 1e-6, group: 'Section' },
  { key: 'f_h', label: 'Uncoupled bending (plunge) frequency', unit: 'Hz', default: 2, min: 1e-3, group: 'Structure', help: 'First bending from the vibration suite when available' },
  { key: 'f_a', label: 'Uncoupled torsion (pitch) frequency', unit: 'Hz', default: 8, min: 1e-3, group: 'Structure' },
  { key: 'zeta', label: 'Structural damping ratio', unit: '-', default: 0, min: 0, max: 0.2, group: 'Structure' },
  { key: 'cl', label: 'Lift-curve slope', unit: '1/rad', default: 6.2832, min: 1, max: 7, group: 'Aerodynamics' },
  { key: 'alt_m', label: 'Altitude', unit: 'm', default: 0, min: -500, max: 25000, group: 'Flight' },
  { key: 'V_d_eas', label: 'Design dive speed VD (EAS)', unit: 'm/s', default: 225, min: 1, group: 'Flight' },
  { key: 'kappa', label: 'Cubic torsional stiffening', unit: '1/rad²', default: 20, min: 0, max: 1e4, group: 'Non-linearity', help: 'Restoring moment Kα·(α + κ·α³); 0 disables the limit-cycle study' },
  { key: 'V_lco_frac', label: 'Time-history speed / flutter speed', unit: '-', default: 1.05, min: 0.2, max: 2, group: 'Non-linearity' },
  { key: 'alpha0_deg', label: 'Initial pitch disturbance', unit: 'deg', default: 1, min: 0.01, max: 20, group: 'Initial conditions' },
  { key: 'nSpeeds', label: 'Speed points', unit: '', default: 40, min: 10, max: 400, step: 1, discrete: true, group: 'Numerics' },
  { key: 'nSteps', label: 'Time steps per torsion period', unit: '', default: 40, min: 16, max: 400, step: 1, discrete: true, group: 'Numerics' },
];
const secI = (i) => ({ ...i, wh: TAU * i.f_h, wa: TAU * i.f_a });
function sectionFlutter(i, rho, nSp) {
  const s = secI(i), Uref = s.b * s.wa, Ud = Math.sqrt((2 * s.m * (s.b * s.r_alpha) ** 2 * s.wa ** 2) / (rho * 2 * s.b * s.cl * s.b * (0.5 + s.a) || 1e-12));
  const top = (fin(Ud) && s.a > -0.5 ? Math.min(1.3 * Ud, 12 * Uref) : 8 * Uref), Us = N.linspace(top / nSp, top, nSp);
  return { ...pkSweep((U) => sectionModel(s, s.m, U), rho, Us, [s.wh, s.wa]), Us, Ud: s.a > -0.5 ? Ud : NaN, s };
}
/** RK4 time march of the Jones state-space model with a cubic torsional spring; returns histories and the settled pitch amplitude. */
function lcoRun(s, rho, U, kappa, alpha0, nPer, nSteps) {
  const ss = stateSpace(sectionModel(s, s.m, U), rho), A = ss.A, nx = ss.nx, Ka = s.m * (s.b * s.r_alpha) ** 2 * s.wa ** 2, g0 = ss.Mi[0][1] * Ka * kappa, g1 = ss.Mi[1][1] * Ka * kappa;
  const f = (x) => { const d = N.matvec(A, x); const c3 = x[1] ** 3; d[2] -= g0 * c3; d[3] -= g1 * c3; return d; };
  const T = TAU / s.wa, dt = T / nSteps, n = Math.round(nPer * nSteps); let x = new Array(nx).fill(0); x[1] = alpha0;
  const t = [0], al = [alpha0], hh = [0]; let amp = 0, blown = false;
  for (let k = 1; k <= n; k++) {
    const k1 = f(x), k2 = f(N.vadd(x, k1, dt / 2)), k3 = f(N.vadd(x, k2, dt / 2)), k4 = f(N.vadd(x, k3, dt));
    x = x.map((v, j) => v + (dt / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
    if (!fin(x[1]) || Math.abs(x[1]) > 10) { blown = true; break; }
    t.push(k * dt); al.push(x[1]); hh.push(x[0]); if (k > 0.8 * n) amp = Math.max(amp, Math.abs(x[1]));
  }
  return { t, al, hh, amp: blown ? Infinity : amp, blown };
}
const section = {
  id: 'section', title: 'Typical-section flutter (p-k, Theodorsen) and limit cycles', fidelity: 'reduced-order',
  summary: 'Pitch–plunge flutter of a representative wing section with exact Theodorsen aerodynamics by the p-k method, the section divergence speed, and time-domain response with a hardening torsional spring.',
  equations: ['Theodorsen unsteady airfoil theory', 'Wagner indicial lift formulation', 'Flutter eigenvalue equations', 'Divergence equilibrium equations', 'Coupled aerodynamic–structural equations of motion', 'Euler–Lagrange equations'],
  applicable: hasWing, inputs: SEC_INPUTS,
  defaults: (c, up, d) => {
    const w = wingDefaults(c, up, d), wi = { ...Object.fromEntries(WING_INPUTS.map((f) => [f.key, f.default])), ...w }, r = wingRitz(wi, 3, 2, 4), kb = r.tors.findIndex((t) => t < 0.5), kt = r.tors.findIndex((t) => t >= 0.5);
    const ch = wi.c_root * (1 - (1 - wi.taper) * 0.75) * Math.cos(N.rad(wi.sweep_deg)), m0 = (3 * wi.m_semi) / (wi.semi_span * (1 + wi.taper + wi.taper ** 2));
    return { b: ch / 2, a: 2 * wi.x_ea - 1, x_alpha: 2 * (wi.x_cg - wi.x_ea), r_alpha: 2 * wi.r_alpha, m: m0 * (1 - (1 - wi.taper) * 0.75) ** 2, f_h: up.vibration?.f1_Hz ?? (kb >= 0 ? r.wn[kb] / TAU : undefined), f_a: up.vibration?.f_torsion_Hz ?? (kt >= 0 ? r.wn[kt] / TAU : undefined), cl: w.cl, alt_m: w.alt_m, V_d_eas: w.V_d_eas };
  },
  run(i) {
    const at = isa(i.alt_m), rho = at.rho, r = sectionFlutter(i, rho, Math.round(i.nSpeeds)), s = r.s, Uf = r.flutter?.P ?? NaN, wf = r.flutter?.w ?? NaN, sg = Math.sqrt(at.sigma), mu = s.m / (Math.PI * rho * s.b ** 2);
    const names = ['Bending branch', 'Torsion branch'], warnings = [];
    // time histories and limit-cycle amplitude sweep above the linear flutter speed
    const U1 = fin(Uf) ? i.V_lco_frac * Uf : r.Us[r.Us.length - 1] * 0.5, a0 = N.rad(i.alpha0_deg), th = lcoRun(s, rho, U1, i.kappa, a0, 60, Math.round(i.nSteps));
    const fr = [1.02, 1.05, 1.1, 1.15, 1.2, 1.3], lco = fin(Uf) && i.kappa > 0 ? fr.map((x) => lcoRun(s, rho, x * Uf, i.kappa, a0, 60, Math.round(i.nSteps)).amp) : [];
    const Me = Uf / at.a;
    if (!fin(Uf)) warnings.push('No flutter crossing was found up to the highest speed analysed.');
    pkNotes(r, warnings);
    if (Me > 0.7) warnings.push(`Flutter Mach number is ${Me.toFixed(2)}: incompressible Theodorsen theory is not valid in the transonic range, where the flutter boundary dips.`);
    if (fin(r.div) && (!fin(Uf) || r.div < Uf)) warnings.push('Divergence occurs before flutter for this section.');
    if (i.kappa > 0 && th.blown) warnings.push('The time history diverged: the cubic stiffening is too weak to bound the motion at this speed.');
    const m = margin(Uf * sg, i.V_d_eas), ds = (a) => { const st = Math.max(1, Math.ceil(a.length / 400)); return a.filter((_, k) => k % st === 0); };
    return {
      kpis: [
        { key: 'V_flutter_section_ms', label: 'Section flutter speed (TAS)', value: Uf, unit: 'm/s' },
        { key: 'V_flutter_section_eas', label: 'Section flutter speed (EAS)', value: Uf * sg, unit: 'm/s' },
        { key: 'f_flutter_section_Hz', label: 'Flutter frequency', value: wf / TAU, unit: 'Hz' },
        { key: 'k_flutter', label: 'Reduced frequency at flutter', value: (wf * s.b) / Uf, unit: '-' },
        { key: 'U_flutter_reduced', label: 'Reduced flutter speed U/(b·ωα)', value: Uf / (s.b * s.wa), unit: '-' },
        { key: 'V_div_section_ms', label: 'Section divergence speed (TAS)', value: r.div, unit: 'm/s' },
        { key: 'section_margin', label: 'Section flutter margin', value: m, unit: '-', status: !fin(Uf) || m > 0 ? 'ok' : 'bad', note: MARGIN_NOTE },
        { key: 'mass_ratio', label: 'Mass ratio μ = m/(πρb²)', value: mu, unit: '-' },
        { key: 'freq_ratio', label: 'Frequency ratio ωh/ωα', value: s.wh / s.wa, unit: '-', status: s.wh / s.wa < 0.8 ? 'ok' : 'warn', note: 'Flutter speed is lowest as the ratio approaches 1' },
        { key: 'lco_amp_deg', label: 'Pitch amplitude at the time-history speed', value: N.deg(th.amp), unit: 'deg', note: 'Settled amplitude over the last fifth of the record' },
      ].filter((k) => fin(k.value) || k.key === 'V_flutter_section_ms'),
      plots: [
        { type: 'line', title: 'V–g diagram', xlabel: 'True airspeed [m/s]', ylabel: 'Damping g = 2σ/ω [-]', series: r.tracks.map((tk, k) => ({ name: names[k], x: r.Us, y: tk.map((l) => (l[1] > 1e-6 ? N.clamp((2 * l[0]) / l[1], -2, 2) : NaN)) })), annotations: [{ y: 0, label: 'Flutter boundary' }, ...(fin(Uf) ? [{ x: Uf, label: 'Flutter' }] : [])] },
        { type: 'line', title: 'V–f diagram', xlabel: 'True airspeed [m/s]', ylabel: 'Frequency [Hz]', series: r.tracks.map((tk, k) => ({ name: names[k], x: r.Us, y: tk.map((l) => l[1] / TAU) })) },
        { type: 'line', title: `Pitch response at ${U1.toFixed(1)} m/s`, xlabel: 'Time [s]', ylabel: 'Pitch angle [deg]', series: [{ name: i.kappa > 0 ? 'Cubic torsional spring' : 'Linear', x: ds(th.t), y: ds(th.al).map(N.deg) }] },
        ...(lco.length ? [{ type: 'line', title: 'Limit-cycle amplitude above the flutter speed', xlabel: 'Speed / linear flutter speed [-]', ylabel: 'Pitch amplitude [deg]', series: [{ name: 'Time-domain', x: [1, ...fr], y: [0, ...lco.map((v) => (fin(v) ? N.deg(v) : NaN))], style: 'line+points' }, { name: 'Describing function', x: [1, ...fr], y: [0, ...fr.map((x) => N.deg(Math.sqrt(Math.max(0, describing(s, rho, x * Uf) / (0.75 * i.kappa)))))], style: 'dash' }] }] : []),
      ],
      warnings, models: ['Linear flutter model: two-degree-of-freedom typical section', 'p–k flutter model with exact Theodorsen function (Bessel functions)', PK_MODEL, 'Non-linear time-domain flutter model: Jones two-lag approximation of the Wagner function, RK4', 'Limit-cycle oscillation model: cubic hardening spring'],
      assumptions: ['Two-dimensional incompressible potential flow; no thickness, viscosity or shocks', 'One representative section stands for the wing', 'p-k roots are exact at the flutter boundary and approximate away from it', 'No control-surface degree of freedom'],
    };
  },
  convergence: { param: 'nSpeeds', label: 'Speed points in the p-k sweep', levels: [10, 20, 40, 80], metric: 'V_flutter_section_ms' },
  calibration: { params: [{ key: 'f_a', min: 0.01, max: 1e3 }, { key: 'x_alpha', min: -0.3, max: 0.6 }, { key: 'zeta', min: 0, max: 0.1 }], sweep: 'alt_m', target: 'V_flutter_section_ms', note: 'Wind-tunnel or flight flutter-test speeds at several densities; frequencies from a ground vibration test' },
  verify() {
    const c = [[0.1, 0.8319, -0.1723], [0.2, 0.7276, -0.1886], [0.5, 0.5979, -0.1507], [1, 0.5394, -0.1003]], be = bessel(1);
    // benchmark section: a = −1/5, xα = 1/10, μ = 20, rα² = 6/25, ωh/ωα = 2/5
    const rho = 1.225, b = 0.5, i = { b, a: -0.2, x_alpha: 0.1, r_alpha: Math.sqrt(0.24), m: 20 * Math.PI * rho * b * b, f_h: 0.4 * 10 / TAU, f_a: 10 / TAU, zeta: 0, cl: TAU }, r = sectionFlutter(i, rho, 40), s = r.s;
    // Jones state-space flutter speed by bisection on eigenvalue growth
    const Uj = N.brent((U) => growth(sectionModel(s, s.m, U), rho), 0.6 * r.flutter.P, 1.3 * r.flutter.P, 1e-6);
    const amp = lcoRun(s, rho, 1.1 * r.flutter.P, 50, 0.01, 150, 40).amp, df = Math.sqrt(describing(s, rho, 1.1 * r.flutter.P) / (0.75 * 50));
    return [
      N.check('Bessel J0(1)', be.j0, 0.7651976866, 1e-9, 'Abramowitz & Stegun'), N.check('Bessel Y1(1)', be.y1, -0.7812128213, 1e-9, 'Abramowitz & Stegun'),
      N.check('Bessel Wronskian at x = 20 (asymptotic branch)', bessel(20).j1 * bessel(20).y0 - bessel(20).j0 * bessel(20).y1, 2 / (Math.PI * 20), 1e-9, 'J1·Y0 − J0·Y1 = 2/(πx)'),
      ...c.map(([k, F]) => N.check(`Theodorsen F(${k})`, theodorsen(k)[0], F, 2e-4, 'Theodorsen (1935), NACA Report 496 tables')),
      ...c.map(([k, , G]) => N.check(`Theodorsen G(${k})`, theodorsen(k)[1], G, 6e-4, 'Theodorsen (1935), NACA Report 496 tables')),
      N.check('Flutter root is neutrally stable at the p-k flutter speed', pkPoint(sectionModel(s, s.m, r.flutter.P), rho, [0, r.flutter.w])[0] / r.flutter.w, 0, 1e-5, 'Definition of the flutter boundary'),
      N.check('Section divergence speed b·ωθ·r·sqrt(μ/(1+2a))', r.div / (b * 10), Math.sqrt(0.24 * 20 / 0.6), 1e-5, 'Typical-section divergence, exact'),
      N.check('Jones two-lag time-domain model reproduces the Theodorsen flutter speed', Uj, r.flutter.P, 0.02, 'R. T. Jones (1940) approximation of the Wagner function'),
      N.check('Limit-cycle amplitude agrees with the describing-function estimate', amp, df, 0.1, 'First-harmonic balance: Kα,eq = Kα(1 + ¾κA²)'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (fin(o.V_flutter_section_ms) && o.section_margin < 0) out.push({ severity: 'critical', title: 'Section flutter inside 1.15 VD', detail: `Flutter at ${o.V_flutter_section_eas.toFixed(0)} m/s EAS against 1.15·VD = ${(1.15 * i.V_d_eas).toFixed(0)} m/s.`, action: 'Raise torsional stiffness, move the section mass centre forward (mass balance), or separate bending and torsion frequencies; confirm with the multi-mode wing analysis.', basis: 'Flutter clearance to 1.15 VD (CS/FAR 25.629 principle)' });
    if (o.freq_ratio > 0.8) out.push({ severity: 'warn', title: 'Bending and torsion frequencies nearly coincide', detail: `ωh/ωα = ${o.freq_ratio.toFixed(2)}.`, action: 'Increase GJ or reduce the pitch inertia to lift the torsion frequency.', basis: 'Frequency coalescence mechanism of classical flutter' });
    if (o.lco_amp_deg > 0.5 && i.V_lco_frac > 1) out.push({ severity: 'advise', title: 'Limit-cycle oscillation beyond the linear boundary', detail: `Pitch amplitude ${o.lco_amp_deg.toFixed(1)}° at ${i.V_lco_frac} × flutter speed.`, action: 'A bounded limit cycle still consumes fatigue life (Suite 9); treat the linear flutter speed as the limit.', basis: 'Non-linear time-domain response' });
    return out;
  },
};
/** Describing-function estimate: the ¾κA² stiffening that puts the linear flutter speed at U (returns ¾κA², ≥ 0). */
function describing(s, rho, U) {
  const g = (x) => growth(sectionModel({ ...s, wa: s.wa * Math.sqrt(1 + x) }, s.m, U), rho);
  if (g(0) <= 0) return 0;
  let hi = 0.5; for (let k = 0; k < 8 && g(hi) > 0; k++) hi *= 2;
  return g(hi) > 0 ? NaN : N.brent(g, 0, hi, 1e-7);
}

// ---- 2. wing static aeroelasticity ----------------------------------------------------------
const wingStatic = {
  id: 'static', title: 'Wing divergence, aileron reversal and load redistribution', fidelity: 'reduced-order',
  summary: 'Static aeroelastic behaviour of the flexible wing by strip theory on assumed bending and torsion modes: divergence speed, aileron reversal speed, roll effectiveness and the inboard shift of lift.',
  equations: ['Divergence equilibrium equations', 'Control reversal relations', 'Coupled aerodynamic–structural equations of motion', 'Modal superposition equations'],
  applicable: hasWing,
  inputs: [...WING_INPUTS,
    { key: 'V', label: 'Flight speed for the load distribution (TAS)', unit: 'm/s', default: 230, min: 1, group: 'Flight' },
    { key: 'mass_kg', label: 'Aircraft mass', unit: 'kg', default: 70000, min: 0.1, group: 'Flight' },
    { key: 'ail_in', label: 'Aileron inboard end / semi-span', unit: '-', default: 0.7, min: 0.05, max: 0.95, group: 'Aileron' },
    { key: 'ail_out', label: 'Aileron outboard end / semi-span', unit: '-', default: 0.95, min: 0.1, max: 1, group: 'Aileron' },
    { key: 'aileron_frac', label: 'Aileron chord / wing chord', unit: '-', default: 0.25, min: 0.05, max: 0.5, group: 'Aileron' },
    { key: 'k_visc', label: 'Aileron effectiveness factor', unit: '-', default: 0.8, min: 0.3, max: 1, group: 'Aileron', help: 'Viscous and gap loss applied to the thin-aerofoil flap derivatives' },
    { key: 'nModes', label: 'Assumed modes per motion', unit: '', default: 4, min: 1, max: 6, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => wingDefaults(c, up, d),
  run(i) {
    const at = isa(i.alt_m), rho = at.rho, nm = Math.round(i.nModes), ns = 40, W = wingRitz(i, nm, nm, ns), n = W.n, sg = Math.sqrt(at.sigma);
    const m1 = W.at(1), Kh = aeroStiff(m1, rho), sv1 = stripVecs(m1, rho), Ki = N.inv(W.K); // K̂ per unit V²
    // divergence: (K + V²K̂)q = 0 → eigenvalues ν of K⁻¹K̂, V_D = sqrt(−1/ν) for the most negative real ν
    const nus = ceig(N.matmul(Ki, Kh)).filter((l) => Math.abs(l[1]) < 1e-9 * Math.abs(l[0]) && l[0] < 0).map((l) => l[0]), Vd0 = nus.length ? Math.sqrt(-1 / N.amin(nus)) : Infinity, Vdiv = Vd0 > 20 * at.a ? Infinity : Vd0; // beyond 20× the speed of sound is no divergence in any practical sense
    // thin-aerofoil flap derivatives
    const thh = Math.acos(2 * i.aileron_frac - 1), cld = i.k_visc * 2 * (Math.PI - thh + Math.sin(thh)) * (i.cl / TAU), cmd = -i.k_visc * 0.5 * Math.sin(thh) * (1 - Math.cos(thh));
    const inA = (s) => s.e >= i.ail_in && s.e <= i.ail_out;
    /** Elastic solution at speed V for load case: rigid incidence alpha(s) [rad] and aileron angle delta. Returns strip lifts. */
    const solve = (V, alpha, delta) => {
      const rhs = new Array(n).fill(0), Uc = V * W.cs;
      sv1.forEach((v) => { const s = v.s, q2b = rho * Uc * Uc * s.b * s.dy, al = alpha(s), La = q2b * cld * delta * (inA(s) ? 1 : 0), Ma = q2b * 2 * s.b * cmd * delta * (inA(s) ? 1 : 0);
        for (let k = 0; k < n; k++) rhs[k] += -(v.Fc[k] * V * Uc * al) - La * (s.ph[k] - s.b * (s.a + 0.5) * s.pa[k]) + Ma * s.pa[k]; });
      let q; try { q = N.solve(N.madd(W.K, Kh, V * V), rhs); } catch { q = new Array(n).fill(NaN); }
      const lift = sv1.map((v) => { const s = v.s, f0 = s.cl * rho * Uc * s.b * s.dy; return f0 * (V * N.dot(v.Wk, q) + Uc * alpha(s)) + rho * Uc * Uc * s.b * s.dy * cld * delta * (inA(s) ? 1 : 0); });
      return { q, lift, roll: N.sum(lift.map((l, k) => l * sv1[k].s.y)), L: N.sum(lift), twist: N.sum(q.slice(nm)), tipDefl: -N.sum(q.slice(0, nm)) };
    };
    const rigidRoll = (V) => N.sum(sv1.map((v) => (inA(v.s) ? rho * (V * W.cs) ** 2 * v.s.b * v.s.dy * cld * v.s.y : 0)));
    const eff = (V) => solve(V, () => 0, 1).roll / rigidRoll(V), top = Math.min(fin(Vdiv) ? 0.98 * Vdiv : Infinity, 3 * i.V_d_eas / sg, 1.2 * at.a * 3);
    const Vs = N.linspace(top / 40, top, 40), E = Vs.map(eff); let Vrev = Infinity;
    for (let k = 1; k < Vs.length; k++) if (E[k - 1] > 0 && E[k] <= 0) { Vrev = N.brent(eff, Vs[k - 1], Vs[k], 1e-6 * top); break; }
    // lift-curve-slope ratio and spanwise load at the flight point
    const Vf = Math.min(i.V, 0.95 * top), flex = solve(Vf, () => 1, 0), rig = sv1.map((v) => v.s.cl * rho * (Vf * W.cs) ** 2 * v.s.b * v.s.dy), Lr = N.sum(rig), clRatio = flex.L / Lr;
    const aTrim = (i.mass_kg * G0) / 2 / flex.L, eta = sv1.map((v) => v.s.e), ycp = (l) => N.sum(l.map((x, k) => x * sv1[k].s.y)) / N.sum(l), liftEff = Vs.map((V) => solve(V, () => 1, 0).L / N.sum(sv1.map((v) => v.s.cl * rho * (V * W.cs) ** 2 * v.s.b * v.s.dy)));
    const effD = eff(Math.min(i.V_d_eas / sg, 0.97 * top)), warnings = [], Md = Vdiv / at.a;
    if (!fin(Vdiv)) warnings.push('No divergence: wash-out from aft sweep (or an elastic axis ahead of the aerodynamic centre) keeps the wing statically stable at every speed.');
    if (!fin(Vrev)) warnings.push('No aileron reversal was found below the divergence speed or the top of the speed range.');
    if ((fin(Md) && Md > 0.7) || Vrev / at.a > 0.7) warnings.push('Critical speeds are in the compressible range: strip theory with a constant lift-curve slope is optimistic there; correct for Mach number.');
    if (i.V > 0.95 * top) warnings.push('The flight speed is at or beyond the static stability limit; the load distribution is shown at a reduced speed.');
    return {
      kpis: [
        { key: 'V_divergence_ms', label: 'Divergence speed (TAS)', value: Vdiv, unit: 'm/s', status: Vdiv * sg > 1.15 * i.V_d_eas ? 'ok' : 'bad' },
        { key: 'V_reversal_ms', label: 'Aileron reversal speed (TAS)', value: Vrev, unit: 'm/s', status: Vrev * sg > 1.15 * i.V_d_eas ? 'ok' : 'bad' },
        { key: 'V_divergence_eas', label: 'Divergence speed (EAS)', value: Vdiv * sg, unit: 'm/s' }, { key: 'V_reversal_eas', label: 'Aileron reversal speed (EAS)', value: Vrev * sg, unit: 'm/s' },
        { key: 'roll_eff_VD', label: 'Aileron roll effectiveness at VD', value: effD, unit: '-', status: effD > 0.5 ? 'ok' : effD > 0 ? 'warn' : 'bad', note: 'Flexible / rigid rolling moment' },
        { key: 'lift_eff', label: 'Flexible / rigid lift-curve slope at the flight speed', value: clRatio, unit: '-' },
        { key: 'tip_twist_deg', label: 'Elastic tip twist in 1 g flight', value: N.deg(flex.twist * aTrim), unit: 'deg', note: 'Nose-up positive' },
        { key: 'tip_defl_m', label: 'Elastic tip deflection in 1 g flight', value: flex.tipDefl * aTrim, unit: 'm' },
        { key: 'cp_shift_pct', label: 'Spanwise centre-of-pressure shift', value: (100 * (ycp(flex.lift) - ycp(rig))) / (i.semi_span * W.cs), unit: '% semi-span', note: 'Negative = inboard, which relieves root bending' },
      ],
      plots: [
        { type: 'line', title: 'Static aeroelastic effectiveness', xlabel: 'True airspeed [m/s]', ylabel: 'Flexible / rigid [-]', series: [{ name: 'Aileron rolling moment', x: Vs, y: E }, { name: 'Lift-curve slope', x: Vs, y: liftEff.map((v) => N.clamp(v, -3, 5)) }], annotations: [{ y: 0, label: 'Reversal' }, { x: i.V_d_eas / sg, label: 'VD' }] },
        { type: 'line', title: 'Spanwise lift at 1 g', xlabel: 'Span fraction [-]', ylabel: 'Lift per unit span [kN/m]', series: [{ name: 'Flexible', x: eta, y: flex.lift.map((l, k) => (l * aTrim) / sv1[k].s.dy / 1e3) }, { name: 'Rigid (same total lift)', x: eta, y: rig.map((l, k) => ((l * i.mass_kg * G0) / 2 / Lr) / sv1[k].s.dy / 1e3), style: 'dash' }] },
        { type: 'line', title: 'Elastic twist and bending at 1 g', xlabel: 'Span fraction [-]', ylabel: 'Twist [deg]', series: [{ name: 'Elastic twist', x: eta, y: sv1.map((v) => N.deg(N.dot(v.s.pa, flex.q) * aTrim)) }] },
      ],
      warnings, models: ['Strip theory with a uniform lift-curve slope', 'Rayleigh–Ritz bending–torsion wing (polynomial assumed modes)', 'Swept-wing bending-slope incidence coupling', 'Thin-aerofoil flap lift and moment derivatives'],
      assumptions: ['Steady incompressible aerodynamics, no spanwise induction beyond the lift-curve-slope correction', 'Wing clamped at the root; straight elastic axis', 'Antisymmetric aileron case evaluated as rolling moment at zero roll rate'],
    };
  },
  convergence: { param: 'nModes', label: 'Assumed modes per motion', levels: [1, 2, 3, 4, 5], metric: 'V_reversal_ms' },
  calibration: { params: [{ key: 'GJ_root', min: 1, max: 1e11 }, { key: 'x_ea', min: 0.2, max: 0.6 }], sweep: 'V', target: 'lift_eff', note: 'Measured flexible-to-rigid lift or roll effectiveness from wind-tunnel or flight loads tests' },
  verify() {
    const b = { semi_span: 6, c_root: 1.5, taper: 1, sweep_deg: 0, x_ea: 0.4, x_cg: 0.4, r_alpha: 0.25, EI_root: 5e6, GJ_root: 8e5, stiff_exp: 0, m_semi: 300, zeta: 0, cl: 5.5, alt_m: 0, V_d_eas: 100, V: 50, mass_kg: 1000, ail_in: 0.6, ail_out: 0.95, aileron_frac: 0.25, k_visc: 1, nModes: 4 };
    const r = N.kv(wingStatic.run(b)), qD = (Math.PI ** 2 * 8e5) / (4 * 36 * 1.5 * (0.15 * 1.5) * 5.5), stiff = N.kv(wingStatic.run({ ...b, GJ_root: 8e11, EI_root: 5e12 }));
    return [
      N.check('Uniform straight wing divergence q = π²GJ/(4L²·c·e·a)', r.V_divergence_ms, Math.sqrt((2 * qD) / isa(0).rho), 1e-4, 'Bisplinghoff, Ashley & Halfman, Aeroelasticity; strip theory exact'),
      N.check('Rigid limit: roll effectiveness → 1', stiff.roll_eff_VD, 1, 1e-4, 'Limit of infinite stiffness'),
      N.check('Rigid limit: lift effectiveness → 1', stiff.lift_eff, 1, 1e-4, 'Limit of infinite stiffness'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], lim = 1.15 * i.V_d_eas;
    if (o.V_divergence_eas < lim) out.push({ severity: 'critical', title: 'Torsional divergence inside 1.15 VD', detail: `Divergence at ${o.V_divergence_eas.toFixed(0)} m/s EAS against ${lim.toFixed(0)} m/s required.`, action: 'Increase GJ, move the elastic axis forward towards the quarter chord, or use aft sweep / bend–twist-coupled laminates (Suite 21).', basis: 'Divergence clearance to 1.15 VD' });
    if (o.V_reversal_eas < lim) out.push({ severity: 'critical', title: 'Aileron reversal inside 1.15 VD', detail: `Reversal at ${o.V_reversal_eas.toFixed(0)} m/s EAS.`, action: 'Stiffen the wing in torsion, move the ailerons inboard, or use spoilers for high-speed roll control.', basis: 'Control reversal clearance' });
    else if (o.roll_eff_VD < 0.5) out.push({ severity: 'warn', title: 'Low roll effectiveness at VD', detail: `Ailerons retain ${(100 * o.roll_eff_VD).toFixed(0)}% of rigid effectiveness.`, action: 'Check the roll-rate requirement at high speed in Suite 4 with this factor applied to Clδa.', basis: 'Static aeroelastic control effectiveness' });
    if (o.cp_shift_pct < -1) out.push({ severity: 'info', title: 'Aeroelastic load relief', detail: `Lift moves ${(-o.cp_shift_pct).toFixed(1)}% of semi-span inboard at 1 g.`, action: 'Pass the flexible load distribution to Suite 2: the root bending relief can be traded for structural mass and fuel burn.', basis: 'Flexible lift redistribution' });
    return out;
  },
};

// ---- 3. multi-mode wing flutter ---------------------------------------------------------------
/** p-k sweep of the Ritz wing W (built here when not supplied) from Vtop/nSp to Vtop; every normal mode is tracked. */
function wingFlutterSolve(i, rho, nm, nSp, Vtop, W = wingRitz(i, nm, nm, 14), opt) {
  const track = W.wn.map((_, k) => k), Vs = N.linspace(Vtop / nSp, Vtop, nSp);
  return { W, Vs, track, ...pkSweep(W.at, rho, Vs, W.wn, opt) };
}
const wingFlutter = {
  id: 'wingflutter', title: 'Wing bending–torsion flutter (multi-mode p-k)', fidelity: 'reduced-order',
  summary: 'Flutter of the cantilever wing with assumed bending and torsion modes and Theodorsen strip aerodynamics: damping and frequency against speed, flutter speed, frequency and margin to 1.15 VD across altitude.',
  equations: ['Flutter eigenvalue equations', 'Theodorsen unsteady airfoil theory', 'Modal superposition equations', 'Linearised structural dynamic equations', 'Coupled aerodynamic–structural equations of motion'],
  applicable: hasWing,
  inputs: [...WING_INPUTS,
    { key: 'nModes', label: 'Assumed modes per motion', unit: '', default: 3, min: 1, max: 5, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nSpeeds', label: 'Speed points', unit: '', default: 30, min: 10, max: 300, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c, up, d) => wingDefaults(c, up, d),
  run(i) {
    const at = isa(i.alt_m), sg = Math.sqrt(at.sigma), nm = Math.round(i.nModes), top0 = Math.max(2.2 * i.V_d_eas / sg, 1.5 * at.a);
    let r = wingFlutterSolve(i, at.rho, nm, Math.round(i.nSpeeds), top0);
    if (!r.flutter) r = wingFlutterSolve(i, at.rho, nm, Math.round(i.nSpeeds), 3 * top0);
    const Vf = r.flutter?.P ?? Infinity, wf = r.flutter?.w ?? NaN, mg = margin(Vf * sg, i.V_d_eas), Vdiv = fin(r.div) ? r.div : Infinity, W = r.W, warnings = [];
    // flutter boundary against altitude (coarser sweep)
    const alts = [0, 3000, 6000, 9000, 12000].filter((h) => h <= Math.max(i.alt_m, 3000) + 3001), W2 = nm > 2 ? wingRitz(i, 2, 2, 14) : r.W, env = alts.map((h) => { const a = isa(h), x = wingFlutterSolve(i, a.rho, Math.min(nm, 2), 12, Math.max(3 * i.V_d_eas / Math.sqrt(a.sigma), 2 * a.a), W2, { ptol: 1e-4, check: false }); return (x.flutter?.P ?? NaN) * Math.sqrt(a.sigma); });
    const kb = W.tors.findIndex((t) => t < 0.5), kt = W.tors.findIndex((t) => t >= 0.5), envOk = alts.map((_, k) => k).filter((k) => fin(env[k]));
    if (envOk.length < alts.length) warnings.push(envOk.length ? 'No flutter was found inside the coarse altitude sweep at some altitudes; those points are left off the flutter margin diagram.' : 'No flutter was found inside the coarse altitude sweep, so the flutter margin diagram shows only the VD lines: the boundary lies above three times VD.');
    if (!fin(Vf)) warnings.push('No flutter crossing was found in the speed range analysed (three times the default range).');
    pkNotes(r, warnings);
    const Mreq = (1.15 * i.V_d_eas) / sg / at.a, Mf = Vf / at.a; // Mach number at the clearance speed 1.15 VD and at flutter
    if (fin(Vf) && Mf > 0.7) warnings.push(Mreq > 0.6 ? `Flutter Mach number ${Mf.toFixed(2)} is beyond the validity of incompressible strip theory and the clearance speed 1.15 VD is itself at Mach ${Mreq.toFixed(2)}: the transonic dip can lower the true flutter speed substantially.` : `The flutter speed corresponds to Mach ${Mf.toFixed(2)}, where incompressible theory no longer holds, but the clearance speed 1.15 VD is only Mach ${Mreq.toFixed(2)}: the theory is valid over the whole envelope, in which no flutter is found, and the flutter speed itself is a nominal value.`);
    if (Math.abs(i.sweep_deg) > 15) warnings.push('Swept wing: only the bending-slope incidence term of swept strip theory is retained; expect larger uncertainty.');
    const lbl = (k) => `Mode ${k + 1} (${W.tors[r.track[k]] >= 0.5 ? 'torsion' : 'bending'}, ${(W.wn[r.track[k]] / TAU).toPrecision(3)} Hz)`;
    return {
      kpis: [
        { key: 'V_flutter_ms', label: 'Flutter speed (TAS)', value: Vf, unit: 'm/s' },
        { key: 'V_flutter_eas', label: 'Flutter speed (EAS)', value: Vf * sg, unit: 'm/s' },
        { key: 'f_flutter_Hz', label: 'Flutter frequency', value: wf / TAU, unit: 'Hz' },
        { key: 'flutter_margin', label: 'Flutter margin', value: fin(mg) ? mg : 9.99, unit: '-', status: mg > 0 ? 'ok' : 'bad', note: MARGIN_NOTE + (fin(mg) ? '' : '; no flutter found, reported as the cap 9.99') },
        { key: 'M_flutter', label: 'Flutter Mach number (incompressible theory)', value: Mf, unit: '-', status: Mf < 0.7 || Mreq <= 0.6 ? 'ok' : 'warn', note: `Clearance speed 1.15 VD is Mach ${Mreq.toFixed(2)}; compressibility matters when both exceed about 0.6–0.7` },
        { key: 'M_clearance', label: 'Mach number at 1.15 VD', value: Mreq, unit: '-' },
        { key: 'V_div_dyn_ms', label: 'Divergence speed from the same model', value: Vdiv, unit: 'm/s' },
        { key: 'f_bend1_Hz', label: 'First bending mode in vacuum', value: W.wn[kb] / TAU, unit: 'Hz' },
        { key: 'f_tors1_Hz', label: 'First torsion mode in vacuum', value: (W.wn[kt] ?? NaN) / TAU, unit: 'Hz' },
        { key: 'flutter_mode', label: 'Branch that goes unstable', value: r.flutter ? r.flutter.mode + 1 : 0, unit: '-' },
      ].filter((k) => fin(k.value) || ['V_flutter_ms', 'f_flutter_Hz'].includes(k.key)),
      plots: [
        { type: 'line', title: 'V–g diagram', xlabel: 'True airspeed [m/s]', ylabel: 'Damping g = 2σ/ω [-]', series: r.tracks.slice(0, 6).map((tk, k) => ({ name: lbl(k), x: r.Vs, y: tk.map((l) => (l[1] > 1e-6 ? N.clamp((2 * l[0]) / l[1], -1.5, 1.5) : NaN)) })), annotations: [{ y: 0, label: 'Flutter boundary' }, { x: (1.15 * i.V_d_eas) / sg, label: '1.15 VD' }] },
        { type: 'line', title: 'V–f diagram', xlabel: 'True airspeed [m/s]', ylabel: 'Frequency [Hz]', series: r.tracks.slice(0, 6).map((tk, k) => ({ name: lbl(k), x: r.Vs, y: tk.map((l) => l[1] / TAU) })) },
        { type: 'line', title: 'Flutter margin diagram', xlabel: 'Flutter speed [m/s EAS]', ylabel: 'Altitude [m]', series: [...(envOk.length ? [{ name: 'Flutter boundary', x: envOk.map((k) => env[k]), y: envOk.map((k) => alts[k]), style: 'line+points' }] : []), { name: '1.15 VD', x: alts.map(() => 1.15 * i.V_d_eas), y: alts, style: 'dash' }, { name: 'VD', x: alts.map(() => i.V_d_eas), y: alts, style: 'dash' }] },
      ],
      tables: [{ title: 'In-vacuum normal modes', columns: ['Mode', 'Frequency [Hz]', 'Torsional energy share'], rows: W.wn.slice(0, 6).map((w, k) => [k + 1, w / TAU, W.tors[k]]) }],
      warnings, models: ['p–k flutter model with Theodorsen strip aerodynamics at the local reduced frequency', PK_MODEL, 'Finite-span lift-curve-slope correction of the circulatory terms (modified strip theory)', 'Rayleigh–Ritz structural dynamics model with inertial bending–torsion coupling'],
      assumptions: ['Incompressible, attached flow; no tip-loss or spanwise aerodynamic coupling', 'Cantilever wing without engines, stores or control surfaces', 'Constant true-airspeed sweep at fixed density (matched-point iteration on Mach is not performed)'],
    };
  },
  convergence: { param: 'nModes', label: 'Assumed modes per motion', levels: [1, 2, 3, 4], metric: 'V_flutter_ms' },
  calibration: { params: [{ key: 'GJ_root', min: 1, max: 1e11 }, { key: 'EI_root', min: 1, max: 1e12 }, { key: 'x_cg', min: 0.2, max: 0.7 }, { key: 'zeta', min: 0, max: 0.1 }], sweep: 'alt_m', target: 'V_flutter_ms', note: 'Flutter-model wind-tunnel boundaries or flight flutter-test damping trends; update stiffness from the ground vibration test first' },
  verify() {
    // Goland wing: uniform cantilever, exact strip-theory flutter 137.2 m/s at 70.7 rad/s
    const g = { semi_span: 6.096, c_root: 1.8288, taper: 1, sweep_deg: 0, x_ea: 0.33, x_cg: 0.43, r_alpha: Math.sqrt(8.64 / 35.71) / 1.8288, EI_root: 9.77e6, GJ_root: 0.987e6, stiff_exp: 0, m_semi: 35.71 * 6.096, zeta: 0, cl: TAU };
    const r = wingFlutterSolve(g, 1.225, 3, 30, 220), k = Math.sqrt(9.77e6 / (35.71 * 6.096 ** 4));
    // exact coalescence flutter of two undamped modes joined by a circulatory (non-symmetric) stiffness P·[[0, a], [−b, 0]]:
    // s⁴ + (ω1² + ω2²)s² + ω1²ω2² + a·b·P² = 0 has a double root at P = (ω2² − ω1²)/(2·sqrt(a·b)), ω² = (ω1² + ω2²)/2
    const w1 = 3, w2 = 5, ca = 2, cb = 0.5, co = pkSweep((P) => ({ n: 2, M: [[1, 0], [0, 1]], C: [[0, 0], [0, 0]], K: [[w1 * w1, ca * P], [-cb * P, w2 * w2]], strips: [] }), 1, N.linspace(0.5, 12, 20), [w1, w2]);
    // a tapered wing whose second-bending and torsion branches approach each other well below the flutter speed (the case in
    // which nearest-root tracking alone loses the fluttering branch): p-k against the tracking-free Jones state-space eigenvalues
    const tw = { semi_span: 13.54, c_root: 2.837, taper: 0.59, sweep_deg: 3, x_ea: 0.4, x_cg: 0.45, r_alpha: 0.25, EI_root: 2.669e7, GJ_root: 1.897e7, stiff_exp: 3, m_semi: 2880, zeta: 0.02, cl: 5.316 }, rt = 0.6597;
    const pt = wingFlutterSolve(tw, rt, 3, 30, 480), Wt = pt.W, osc = (V) => N.amax(ceig(stateSpace(Wt.at(V), rt).A).filter((l) => Math.abs(l[1]) > 1).map((l) => l[0])), Vss = N.findRoot(osc, 100, 480, 19, 1e-4);
    const below = nyquistCount(Wt.at(0.95 * pt.flutter.P), rt), above = nyquistCount(Wt.at(1.05 * pt.flutter.P), rt);
    return [
      N.check('Two-mode coalescence flutter: P = (ω2² − ω1²)/(2·sqrt(a·b))', co.flutter.P, (w2 * w2 - w1 * w1) / (2 * Math.sqrt(ca * cb)), 1e-5, 'Exact double root of the biquadratic characteristic equation'),
      N.check('Two-mode coalescence flutter frequency sqrt((ω1² + ω2²)/2)', co.flutter.w, Math.sqrt((w1 * w1 + w2 * w2) / 2), 1e-3, 'Exact double root of the biquadratic characteristic equation'),
      N.check('Tapered wing with approaching branches: p-k flutter speed equals the state-space (Jones) flutter speed', pt.flutter.P, Vss, 0.03, 'Tracking-free eigenvalues of the two-lag Wagner state-space model of the same wing; Jones’ approximation of C(k) is within about 2%'),
      N.check('Every p-k branch is distinct at the top of that sweep (no two trackers on one root)', Math.min(...pt.tracks.flatMap((a, p) => pt.tracks.slice(0, p).map((b) => cdist(a[29], b[29]) / (1 + cmod(a[29]))))) > 1e-4 ? 1 : 0, 1, 0, 'A linear system of n modes has n distinct p-k roots here'),
      N.check('Determinant phase count: no unstable root 5% below the p-k flutter speed', below.Z, 0, 0, 'Argument principle on the flutter determinant along s = iω'),
      N.check('Determinant phase count: one unstable pair 5% above the p-k flutter speed', above.Z, 2, 0, 'Argument principle on the flutter determinant along s = iω'),
      N.check('Goland wing flutter speed 137.2 m/s', r.flutter.P, 137.2, 0.02, 'Goland (1945); Goland & Luke (1948) exact solution'),
      N.check('Goland wing flutter frequency 70.7 rad/s', r.flutter.w, 70.7, 0.02, 'Goland & Luke (1948)'),
      N.check('Ritz torsion frequency of the uniform wing', wingRitz({ ...g, x_cg: 0.33 }, 3, 3, 4).wn[1], (Math.PI / (2 * 6.096)) * Math.sqrt(0.987e6 / 8.64), 2e-3, 'Fixed–free torsion rod, exact'),
      N.check('Ritz first bending frequency of the uniform wing', wingRitz({ ...g, x_cg: 0.33 }, 3, 3, 4).wn[0], 1.875104 ** 2 * k, 2e-3, 'Euler–Bernoulli cantilever, exact'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.flutter_margin < 0) out.push({ severity: 'critical', title: 'Flutter inside 1.15 VD', detail: `Flutter at ${o.V_flutter_eas.toFixed(0)} m/s EAS and ${fin(o.f_flutter_Hz) ? o.f_flutter_Hz.toFixed(1) : '—'} Hz; margin ${(100 * o.flutter_margin).toFixed(0)}%.`, action: 'Raise GJ (thicker skins, ±45° plies), move mass forward of the elastic axis, or add tip/leading-edge mass balance. Re-run Suites 2 and 10 with the changed structure.', basis: 'Flutter clearance to 1.15 VD' });
    else if (o.flutter_margin < 0.15) out.push({ severity: 'warn', title: 'Small flutter margin', detail: `Margin over 1.15 VD is ${(100 * o.flutter_margin).toFixed(0)}%.`, action: 'Strip theory carries roughly ±15% uncertainty: confirm with a doublet-lattice analysis and a ground vibration test before relying on this margin.', basis: 'Model-form uncertainty of strip-theory flutter prediction' });
    if (o.M_flutter > 0.7 && o.M_clearance > 0.6) out.push({ severity: 'advise', title: 'Transonic flutter assessment required', detail: `Predicted flutter Mach ${o.M_flutter.toFixed(2)}; clearance speed 1.15 VD at Mach ${o.M_clearance.toFixed(2)}.`, action: 'Use compressible unsteady aerodynamics (doublet lattice with transonic correction or CFD-based aeroelasticity).', basis: 'Validity limit of incompressible theory' });
    return out;
  },
};

// ---- 4. gust response -----------------------------------------------------------------------
const gust = {
  id: 'gust', title: 'Discrete and continuous gust response (Küssner, Sears)', fidelity: 'reduced-order',
  summary: 'Heave response of the aircraft to a 1-cosine gust with Wagner and Küssner lift build-up, tuned over gust length, and the load-factor response to von Kármán turbulence.',
  equations: ['Küssner gust-response formulation', 'Wagner indicial lift formulation', 'Theodorsen unsteady airfoil theory', 'Coupled aerodynamic–structural equations of motion'],
  applicable: hasWing,
  inputs: [
    { key: 'mass_kg', label: 'Aircraft mass', unit: 'kg', default: 70000, min: 0.1, group: 'Aircraft' },
    { key: 'S', label: 'Wing area', unit: 'm²', default: 122, min: 0.01, group: 'Aircraft' },
    { key: 'mac', label: 'Mean aerodynamic chord', unit: 'm', default: 4.2, min: 0.02, group: 'Aircraft' },
    { key: 'cl', label: 'Lift-curve slope', unit: '1/rad', default: 5, min: 1, max: 7, group: 'Aircraft' },
    { key: 'V', label: 'True airspeed', unit: 'm/s', default: 230, min: 1, group: 'Flight' },
    { key: 'alt_m', label: 'Altitude', unit: 'm', default: 10000, min: -500, max: 25000, group: 'Flight' },
    { key: 'U_ds', label: 'Design gust velocity (EAS)', unit: 'm/s', default: 15.24, min: 0, max: 40, group: 'Gust', help: 'Derived gust velocity for the classical 12.5-chord discrete gust and the Pratt formula: 15.24 m/s (50 ft/s) at VC. The current transport rule instead uses a reference velocity of 17.07 m/s (56 ft/s) EAS at sea level, falling to 13.41 m/s (44 ft/s) at 15 000 ft, with gradient distances from 9 to 107 m (30–350 ft) searched for the worst case' },
    { key: 'H_chords', label: 'Gust gradient distance', unit: 'chords', default: 12.5, min: 1, max: 200, group: 'Gust', help: 'Distance to peak velocity; 12.5 chords in the classical discrete-gust criterion' },
    { key: 'sigma_w', label: 'Turbulence RMS velocity (TAS)', unit: 'm/s', default: 1.5, min: 0, max: 20, group: 'Turbulence' },
    { key: 'L_turb', label: 'Turbulence scale length', unit: 'm', default: 762, min: 10, max: 5000, group: 'Turbulence', help: '762 m (2500 ft) is the usual von Kármán scale at altitude' },
    { key: 'n_limit', label: 'Limit load factor', unit: 'g', default: 2.5, min: 1.2, max: 9, group: 'Limits' },
    { key: 'nSteps', label: 'Time steps per gust length', unit: '', default: 200, min: 40, max: 4000, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c, up, d) => ({ mass_kg: c.mass.mtow_kg, S: c.wing.S_m2 || undefined, mac: d.mac || undefined, cl: up.cfd?.CLa_per_rad ?? helmbold(d.AR, N.rad(c.wing.sweep_deg)), V: c.flight.V_ms, alt_m: c.atm.alt_m, sigma_w: Math.max(0.5, c.atm.turb_intensity * c.flight.V_ms), n_limit: c.aero.n_pos }),
  run(i) {
    const at = isa(i.alt_m), rho = at.rho, b = i.mac / 2, q = 0.5 * rho * i.V * i.V, kL = 0.5 * rho * i.V * i.S * i.cl, ma = (Math.PI * rho * b * b * i.S) / i.mac, mt = i.mass_kg + ma, Ude = i.U_ds / Math.sqrt(at.sigma);
    const a0 = 1 - JONES.A[0] - JONES.A[1], be = (c) => (c * i.V) / b;
    /** 1-cosine gust of gradient H [m]; unsteady = Wagner + Küssner, otherwise quasi-steady. Returns Δn history. */
    const sim = (H, unsteady = true, n0 = Math.round(i.nSteps)) => {
      const n = Math.max(n0, Math.ceil((1.5 * H) / b)), Tg = (2 * H) / i.V, tEnd = 3 * Tg, dt = Tg / n, wg = (t) => (t < Tg ? 0.5 * Ude * (1 - Math.cos((TAU * t) / Tg)) : 0);
      const f = (t, x) => { const w = wg(t), ge = unsteady ? x[3] + x[4] : w, ve = unsteady ? a0 * x[0] + x[1] + x[2] : x[0];
        return [(kL * (ge - ve)) / mt, be(JONES.b[0]) * (JONES.A[0] * x[0] - x[1]), be(JONES.b[1]) * (JONES.A[1] * x[0] - x[2]), be(KUSSNER.b[0]) * (KUSSNER.A[0] * w - x[3]), be(KUSSNER.b[1]) * (KUSSNER.A[1] * w - x[4])]; };
      const r = N.rk4(f, 0, [0, 0, 0, 0, 0], tEnd, 3 * n), dn = r.y.map((x, k) => f(r.t[k], x)[0] / G0);
      return { t: r.t, dn, peak: N.amax(dn), wg: r.t.map(wg) };
    };
    const H0 = i.H_chords * i.mac, u = sim(H0), qs = sim(H0, false), Hs = N.logspace(2, 60, 14).map((x) => x * i.mac), pk = Hs.map((H) => sim(H, true, 80).peak), kw = N.argmax(pk);
    const mug = (2 * i.mass_kg) / (rho * i.mac * i.cl * i.S), Kg = (0.88 * mug) / (5.3 + mug), dnSharp = (kL * Ude) / (i.mass_kg * G0), dnPratt = Kg * dnSharp;
    // continuous turbulence: |Δn/w_g|² with Sears and Theodorsen, von Kármán spectrum
    const om = N.logspace(1e-3, 1e3, 500).map((x) => (x * i.V) / i.L_turb), H2 = om.map((w) => { const k = (w * b) / i.V, v = C.div(C.scale(sears(k), kL), C.add([0, w * mt], C.scale(theodorsen(k), kL))); return ((w * C.abs(v)) / G0) ** 2; });
    const vk = (w) => { const x = (1.339 * i.L_turb * w) / i.V; return ((i.sigma_w ** 2 * i.L_turb) / (Math.PI * i.V)) * (1 + (8 / 3) * x * x) / (1 + x * x) ** (11 / 6); };
    const Pn = om.map((w, k) => H2[k] * vk(w)), sn = Math.sqrt(N.trapz(om, Pn)), Abar = i.sigma_w ? sn / i.sigma_w : 0, N0 = Math.sqrt(N.trapz(om, om.map((w, k) => w * w * Pn[k])) / (sn * sn || 1)) / TAU;
    const nPeak = 1 + u.peak, exceed = sn > 0 ? N0 * Math.exp(-((i.n_limit - 1) ** 2) / (2 * sn * sn)) : 0, warnings = [], st = Math.ceil(u.t.length / 300), ds = (a) => a.filter((_, k) => k % st === 0);
    const gustAng = Math.atan2(Ude, i.V);
    if (nPeak > i.n_limit) warnings.push('The discrete gust load factor exceeds the manoeuvre limit: the structure is gust-critical at this condition.');
    if (gustAng > 0.2) warnings.push(`The gust changes the incidence by ${N.deg(gustAng).toFixed(0)}°: the wing stalls before the linear lift assumed here is reached, so the load factor is an upper bound. The design gust velocity is large for an aircraft this slow.`);
    if (i.V / at.a > 0.75) warnings.push('Compressibility is not included in the indicial functions; use Mach-dependent lift build-up above about Mach 0.7.');
    return {
      kpis: [
        { key: 'dn_gust', label: 'Peak incremental load factor (unsteady)', value: u.peak, unit: 'g' }, { key: 'gust_angle_deg', label: 'Incidence change at the gust peak', value: N.deg(gustAng), unit: 'deg', status: gustAng > 0.2 ? 'warn' : 'ok', note: 'Linear lift holds to roughly 10–12°' },
        { key: 'n_gust_peak', label: 'Peak load factor', value: nPeak, unit: 'g', status: nPeak <= i.n_limit ? 'ok' : 'warn', note: `Limit ${i.n_limit} g` },
        { key: 'dn_quasi', label: 'Quasi-steady 1-cosine result', value: qs.peak, unit: 'g' }, { key: 'dn_sharp', label: 'Sharp-edged rigid result', value: dnSharp, unit: 'g' },
        { key: 'dn_pratt', label: 'Pratt formula result', value: dnPratt, unit: 'g', note: `Kg = ${Kg.toFixed(3)}` },
        { key: 'alleviation', label: 'Gust alleviation factor (this model)', value: u.peak / dnSharp, unit: '-' },
        { key: 'H_tuned_chords', label: 'Most severe gradient distance', value: Hs[kw] / i.mac, unit: 'chords' }, { key: 'dn_tuned', label: 'Load increment at the tuned gust', value: pk[kw], unit: 'g' },
        { key: 'mass_param', label: 'Mass parameter μg', value: mug, unit: '-' },
        { key: 'sigma_n', label: 'RMS load factor in turbulence', value: sn, unit: 'g' }, { key: 'A_bar', label: 'Ā = σn/σw', value: Abar, unit: 'g/(m/s)' }, { key: 'N0_Hz', label: 'Characteristic frequency N0', value: N0, unit: 'Hz' },
        { key: 'exceed_per_h', label: 'Limit-load exceedances per hour in this turbulence', value: 3600 * exceed, unit: '1/h', note: 'Rice: N0·exp(−Δn²/2σ²)' },
      ],
      plots: [
        { type: 'line', title: '1-cosine gust response', xlabel: 'Time [s]', ylabel: 'Incremental load factor [g]', series: [{ name: 'Wagner + Küssner', x: ds(u.t), y: ds(u.dn) }, { name: 'Quasi-steady', x: ds(qs.t), y: ds(qs.dn), style: 'dash' }] },
        { type: 'line', title: 'Gust velocity profile', xlabel: 'Time [s]', ylabel: 'Gust velocity [m/s TAS]', series: [{ name: '1 − cos gust', x: ds(u.t), y: ds(u.wg) }] },
        { type: 'line', title: 'Gust tuning', xlabel: 'Gradient distance [chords]', ylabel: 'Peak incremental load factor [g]', xlog: true, series: [{ name: 'Unsteady model', x: Hs.map((h) => h / i.mac), y: pk, style: 'line+points' }], annotations: [{ y: dnPratt, label: 'Pratt' }] },
        { type: 'line', title: 'Load-factor spectrum in von Kármán turbulence', xlabel: 'Frequency [Hz]', ylabel: 'PSD of load factor [g²/(rad/s)]', xlog: true, ylog: true, series: [{ name: 'Response', x: om.map((w) => w / TAU), y: Pn.map((v) => Math.max(v, 1e-30)) }] },
      ],
      warnings, models: ['Gust-response model: rigid heave with Wagner (Jones) and Küssner indicial functions', 'Sears function for harmonic gusts', 'von Kármán vertical turbulence spectrum', 'Pratt gust formula for comparison'],
      assumptions: ['Heave only: pitch response and wing flexibility are not included (pitch generally relieves, flexibility can amplify wing loads)', 'Unswept gust penetration: the whole span meets the gust at once', 'Incompressible indicial functions'],
    };
  },
  convergence: { param: 'nSteps', label: 'Time steps per gust length', levels: [25, 50, 100, 200], metric: 'dn_gust' },
  verify() {
    // very heavy aircraft: no heave relief, so the load history is the gust convolved with the Küssner function
    const b = { mass_kg: 1e12, S: 20, mac: 2, cl: 5, V: 80, alt_m: 0, U_ds: 10, H_chords: 12.5, sigma_w: 1, L_turb: 762, n_limit: 2.5, nSteps: 400 }, r = N.kv(gust.run(b)), Tg = 50 / 80, s = (t) => (80 * t) / 1;
    const psi = (x) => 1 - 0.5 * Math.exp(-0.13 * x) - 0.5 * Math.exp(-x), duh = (tt) => N.simpson((t) => 0.5 * ((TAU / Tg) * Math.sin((TAU * t) / Tg)) * psi(s(tt - t)), 0, Math.min(tt, Tg), 400);
    let pk = 0; for (const tt of N.linspace(0.2 * Tg, 1.5 * Tg, 300)) pk = Math.max(pk, duh(tt));
    const cj = (k) => C.sub([1, 0], C.add(C.div([0, 0.165 * k], [0.0455, k]), C.div([0, 0.335 * k], [0.3, k])));
    return [
      N.check('Rigid 1-cosine load equals the Duhamel integral with the Küssner function', r.dn_gust / r.dn_sharp, pk, 2e-3, 'Küssner (1936); Bisplinghoff, Ashley & Halfman'),
      N.check('Sears function S(0) = 1', sears(1e-12)[0], 1, 1e-12, 'Sears (1941)'),
      N.check('Sears function modulus tends to 1/sqrt(2πk)', C.abs(sears(40)), 1 / Math.sqrt(TAU * 40), 5e-3, 'Sears (1941) asymptote'),
      N.check('Jones approximation of C(k) at k = 0.2', cj(0.2)[0], theodorsen(0.2)[0], 0.02, 'R. T. Jones (1940)'),
      N.check('von Kármán spectrum integrates to σ²', N.simpson((x) => { const w = Math.exp(x), y = 1.339 * w; return (w / Math.PI) * (1 + (8 / 3) * y * y) / (1 + y * y) ** (11 / 6); }, -12, 12, 4000), 1, 2e-3, 'Definition of the spectrum'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.n_gust_peak > i.n_limit) out.push({ severity: o.gust_angle_deg > 11.5 ? 'advise' : 'warn', title: o.gust_angle_deg > 11.5 ? 'Gust load factor above the limit, but beyond the linear range' : 'Gust-critical condition', detail: `Peak ${o.n_gust_peak.toFixed(2)} g against the ${i.n_limit} g manoeuvre limit${o.gust_angle_deg > 11.5 ? `; the gust incidence of ${o.gust_angle_deg.toFixed(0)}° exceeds the stall margin, so stall caps the real load` : ''}.`, action: 'Carry the gust load factor to Suites 2 and 9; consider gust-load alleviation in Suite 16 or a higher wing loading.', basis: 'Discrete tuned-gust response' });
    if (o.dn_tuned > 1.05 * o.dn_gust) out.push({ severity: 'advise', title: 'A different gust length is more severe', detail: `Tuned gust at ${o.H_tuned_chords.toFixed(0)} chords gives ${o.dn_tuned.toFixed(2)} g against ${o.dn_gust.toFixed(2)} g.`, action: 'Use the tuned-gust value for loads; modern requirements sweep the gradient distance.', basis: 'Tuned discrete gust' });
    return out;
  },
};

// ---- 5. supersonic panel flutter ------------------------------------------------------------
/** Critical dynamic-pressure parameter λ = 2q·a³/(β·D) of a simply supported 2-D panel (Galerkin, nm sine modes). rx = in-plane compression / Euler buckling load. */
function panelLambda(nm, rx) {
  const P = Math.PI, Km = (lam) => N.range(nm, (p) => N.range(nm, (q) => (p === q ? ((p + 1) * P) ** 4 - rx * P * P * ((p + 1) * P) ** 2 : (p + q) % 2 ? (lam * 4 * (p + 1) * (q + 1)) / ((p + 1) ** 2 - (q + 1) ** 2) : 0)));
  const cplx = (lam) => N.amax(ceig(Km(lam)).map((l) => Math.abs(l[1])));
  let lo = 0, hi = 50; while (cplx(hi) < 1e-6 && hi < 1e5) { lo = hi; hi *= 1.5; }
  for (let k = 0; k < 50; k++) { const mid = 0.5 * (lo + hi); if (cplx(mid) > 1e-6) hi = mid; else lo = mid; }
  return { lam: 0.5 * (lo + hi), Km };
}
const panel = {
  id: 'panel', title: 'Supersonic panel flutter (piston theory)', fidelity: 'reduced-order',
  summary: 'Flutter boundary of a flat skin panel in supersonic flow by first-order piston theory and Galerkin sine modes: critical dynamic pressure, minimum gauge and frequency coalescence.',
  equations: ['Piston theory for suitable supersonic applications', 'Flutter eigenvalue equations', 'Modal superposition equations'],
  applicable: hasWing,
  inputs: [
    { key: 'a', label: 'Panel length (streamwise)', unit: 'm', default: 0.4, min: 0.01, group: 'Panel' },
    { key: 't_mm', label: 'Panel thickness', unit: 'mm', default: 1.6, min: 0.05, group: 'Panel' },
    { key: 'E', label: 'Young\'s modulus', unit: 'Pa', default: 73.1e9, min: 1e8, group: 'Panel' }, { key: 'nu', label: 'Poisson ratio', unit: '-', default: 0.33, min: 0, max: 0.49, group: 'Panel' },
    { key: 'rho_m', label: 'Panel density', unit: 'kg/m³', default: 2780, min: 100, group: 'Panel' },
    { key: 'rx', label: 'In-plane compression / buckling load', unit: '-', default: 0, min: -3, max: 0.99, group: 'Panel', help: 'Thermal or mechanical compression lowers the flutter boundary; negative is tension' },
    { key: 'Mach', label: 'Mach number', unit: '-', default: 2, min: 1.3, max: 5, group: 'Flow', help: 'Piston theory needs roughly M > 1.6' },
    { key: 'alt_m', label: 'Altitude', unit: 'm', default: 11000, min: 0, max: 40000, group: 'Flow' },
    { key: 'nModes', label: 'Sine modes', unit: '', default: 6, min: 2, max: 12, step: 1, discrete: true, group: 'Numerics' },
  ],
  defaults: (c) => { const m = METALS[c.struct.material] || METALS['Al 2024-T3']; return { t_mm: Math.min(c.struct.t_skin_mm, 2), E: m.E, nu: m.nu, rho_m: m.rho, a: N.clamp(0.08 * (c.wing.S_m2 / (c.wing.b_m || 1)), 0.05, 0.6), Mach: Math.max(2, c.aero.Mmo || 0) }; },
  run(i, ctx) {
    const at = isa(i.alt_m), q = 0.7 * at.p * i.Mach ** 2, be = Math.sqrt(i.Mach ** 2 - 1), t = i.t_mm / 1e3, D = (i.E * t ** 3) / (12 * (1 - i.nu ** 2)), nm = Math.round(i.nModes), r = panelLambda(nm, i.rx);
    const lam = (2 * q * i.a ** 3) / (be * D), qcr = (r.lam * be * D) / (2 * i.a ** 3), tReq = ((2 * q * i.a ** 3 * 12 * (1 - i.nu ** 2)) / (r.lam * be * i.E)) ** (1 / 3), w0 = Math.sqrt(D / (i.rho_m * t * i.a ** 4));
    const ls = N.linspace(0, 1.15 * r.lam, 60), f = [[], []]; for (const l of ls) { const e = ceig(r.Km(l)).map((x) => x[0]).sort((x, y) => x - y); f[0].push((Math.sqrt(Math.max(e[0], 0)) * w0) / TAU); f[1].push((Math.sqrt(Math.max(e[1], 0)) * w0) / TAU); }
    const Mmo = ctx?.case?.aero?.Mmo || 0, warnings = [];
    if (Mmo < 1) warnings.push(`This vehicle is subsonic (MMO ${Mmo || 'not set'}): panel flutter is not a design case for it. The result is a reference calculation at the Mach number entered.`);
    if (i.Mach < 1.6) warnings.push('Below about Mach 1.6 first-order piston theory is inaccurate and single-mode (low supersonic) panel flutter can occur; use linearised potential flow.');
    if (i.rx > 0.9) warnings.push('The panel is close to buckling: post-buckled panels need a non-linear (von Kármán plate) analysis.');
    return {
      kpis: [
        { key: 'lambda_cr', label: 'Critical dynamic-pressure parameter λcr', value: r.lam, unit: '-', note: 'λ = 2q·a³/(β·D)' },
        { key: 'lambda', label: 'λ at this flight condition', value: lam, unit: '-' },
        { key: 'panel_margin', label: 'Dynamic-pressure margin qcr/q − 1', value: qcr / q - 1, unit: '-', status: qcr / q > 1.5 ? 'ok' : qcr > q ? 'warn' : 'bad' },
        { key: 'q_cr_Pa', label: 'Critical dynamic pressure', value: qcr, unit: 'Pa' }, { key: 'q_Pa', label: 'Flight dynamic pressure', value: q, unit: 'Pa' },
        { key: 't_min_mm', label: 'Minimum flutter-free thickness', value: tReq * 1e3, unit: 'mm' },
        { key: 'f_panel1_Hz', label: 'First panel frequency in vacuum', value: (Math.PI ** 2 * Math.sqrt(1 - i.rx) * w0) / TAU, unit: 'Hz' },
        { key: 'f_flutter_panel_Hz', label: 'Coalescence (flutter) frequency', value: f[0][Math.round(59 / 1.15)], unit: 'Hz' },
      ],
      plots: [{ type: 'line', title: 'Frequency coalescence', xlabel: 'Dynamic-pressure parameter λ [-]', ylabel: 'Frequency [Hz]', series: [{ name: 'Mode 1', x: ls, y: f[0] }, { name: 'Mode 2', x: ls, y: f[1] }], annotations: [{ x: r.lam, label: 'Flutter' }, { x: lam, label: 'Flight' }] },
        { type: 'line', title: 'Minimum thickness against Mach number', xlabel: 'Mach number [-]', ylabel: 'Thickness [mm]', series: [{ name: 'Flutter-free gauge', x: N.linspace(1.5, 4, 26), y: N.linspace(1.5, 4, 26).map((M) => 1e3 * ((2 * 0.7 * at.p * M * M * i.a ** 3 * 12 * (1 - i.nu ** 2)) / (r.lam * Math.sqrt(M * M - 1) * i.E)) ** (1 / 3)) }], annotations: [{ y: i.t_mm, label: 'Panel' }] }],
      warnings, models: ['First-order piston theory (static aerodynamic term, Ackeret limit)', 'Galerkin solution with simply supported beam modes', 'Frequency-coalescence flutter criterion'],
      assumptions: ['Two-dimensional (infinite-width) flat panel, simply supported, flow over one side', 'Aerodynamic and structural damping neglected (slightly conservative)', 'Linear plate: no curvature, cavity or post-buckling effects'],
    };
  },
  convergence: { param: 'nModes', label: 'Sine modes', levels: [2, 4, 6, 8], metric: 'lambda_cr' },
  verify() {
    return [N.check('Simply supported 2-D panel λcr = 343.36', panelLambda(8, 0).lam, 343.36, 5e-4, 'Hedgepeth (1957); Dowell, Aeroelasticity of Plates and Shells'),
      N.check('Two-mode Galerkin λcr = 45π⁴/16', panelLambda(2, 0).lam, (3 * 15 * Math.PI ** 4) / 16, 1e-6, 'Coalescence of (π⁴, 16π⁴) with coupling ±8λ/3: λ = 45π⁴/16')];
  },
  recommend(res, i) {
    const o = res.outputs;
    return o.panel_margin < 0.5 ? [{ severity: o.panel_margin < 0 ? 'critical' : 'warn', title: 'Insufficient panel flutter margin', detail: `qcr/q − 1 = ${(100 * o.panel_margin).toFixed(0)}%; flutter-free gauge is ${o.t_min_mm.toFixed(2)} mm against ${i.t_mm} mm.`, action: 'Thicken the skin, shorten the panel with an extra stiffener or frame, or orient stiffeners streamwise; include thermal compression at the hot condition.', basis: 'Piston-theory panel flutter boundary' }] : [];
  },
};

// ---- 6 & 7. rotor / propeller blade pitch–flap stability ------------------------------------
const BLADE_INPUTS = [
  { key: 'R', label: 'Blade radius', unit: 'm', default: 8, min: 0.02, group: 'Blade' },
  { key: 'chord', label: 'Blade chord', unit: 'm', default: 0.5, min: 0.002, group: 'Blade' },
  { key: 'm_blade', label: 'Blade mass', unit: 'kg', default: 110, min: 1e-5, group: 'Blade' },
  { key: 'e_hinge', label: 'Flap hinge offset / R', unit: '-', default: 0.05, min: 0, max: 0.4, group: 'Blade' },
  { key: 'f_flap_nr', label: 'Non-rotating flap frequency', unit: 'Hz', default: 0, min: 0, group: 'Blade', help: '0 for an articulated blade; the first cantilever frequency for hingeless rotors and propellers' },
  { key: 'f_tors', label: 'Non-rotating pitch/torsion frequency', unit: 'Hz', default: 20, min: 0.01, group: 'Blade', help: 'Set by control-system stiffness and blade torsional stiffness; typically 3–8 /rev on helicopters' },
  { key: 'x_pa', label: 'Pitch axis / chord', unit: '-', default: 0.25, min: 0.1, max: 0.5, group: 'Section' },
  { key: 'x_cg', label: 'Section mass centre / chord', unit: '-', default: 0.27, min: 0.1, max: 0.6, group: 'Section', help: 'Blades are mass-balanced to about 25% chord; aft positions are destabilising' },
  { key: 'r_alpha', label: 'Pitch radius of gyration / chord', unit: '-', default: 0.25, min: 0.1, max: 0.5, group: 'Section' },
  { key: 'zeta', label: 'Pitch structural damping ratio', unit: '-', default: 0.01, min: 0, max: 0.3, group: 'Section' },
  { key: 'cl', label: 'Section lift-curve slope', unit: '1/rad', default: 5.73, min: 1, max: 7, group: 'Aerodynamics' },
  { key: 'rpm', label: 'Operating speed', unit: 'rpm', default: 258, min: 1, group: 'Operation' },
  { key: 'alt_m', label: 'Altitude', unit: 'm', default: 0, min: -500, max: 15000, group: 'Operation' },
];
function bladeDefaults(c) {
  const r = c.rotor, heli = c.meta.type === 'helicopter', rot = r.R_m > 0, R = rot ? r.R_m : c.prop.prop_dia_m / 2;
  if (!(R > 0)) return {};
  const chord = rot && r.chord_m > 0 ? r.chord_m : 0.14 * R, t = 0.12 * chord, m = rot && r.blade_mass_kg > 0 ? r.blade_mass_kg : 163 * chord * chord * R, rpm = (rot ? r.rpm : c.prop.rpm) || 2000, e = heli ? r.hinge_offset : 0.12, l = R * (1 - e);
  const fnr = heli ? 0 : (0.5596 * Math.sqrt((70e9 * 0.011 * chord * t ** 3) / ((m / l) * l ** 4))); // first cantilever mode of the section-shape stiffness estimate
  return { R, chord, m_blade: m, e_hinge: e, f_flap_nr: fnr, f_tors: heli ? (4.5 * rpm) / 60 : Math.max((4 * rpm) / 60, 6 * fnr), cl: (rot && r.cla) || 5.73, rpm, alt_m: c.atm.alt_m };
}
const rotor = {
  id: 'rotor', title: 'Rotor blade pitch–flap flutter and divergence', fidelity: 'reduced-order',
  summary: 'Coupled flapping and pitching of a rotating blade with centrifugal coupling and unsteady strip aerodynamics, solved by the p-k method against rotor speed: flutter and divergence speeds and modal damping.',
  equations: ['Flutter eigenvalue equations', 'Theodorsen unsteady airfoil theory', 'Divergence equilibrium equations', 'Coupled aerodynamic–structural equations of motion', 'Euler–Lagrange equations'],
  applicable: hasRotor,
  inputs: [...BLADE_INPUTS, { key: 'rpm_max_frac', label: 'Top speed / operating', unit: '-', default: 1.5, min: 1.05, max: 4, group: 'Operation' }, { key: 'nSpeeds', label: 'Rotor speed points', unit: '', default: 30, min: 8, max: 300, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c) => bladeDefaults(c),
  run(i, ctx) {
    const rho = isa(i.alt_m).rho, Om0 = (i.rpm * TAU) / 60, top = Om0 * i.rpm_max_frac, nS = Math.round(i.nSpeeds), Os = N.linspace(top / nS, top, nS), m0 = bladeModel(i, Os[0]);
    const r = pkSweep((Om) => bladeModel(i, Om), rho, Os, m0.wn), rpm = Os.map((o) => (o * 60) / TAU), Of = r.flutter?.P ?? Infinity, Od = fin(r.div) ? r.div : Infinity, lim = Math.min(Of, Od);
    const op = bladeModel(i, Om0), lock = (rho * i.cl * i.chord * i.R ** 4) / op.Ib, rt = pkModes(op, rho, op.wn.map((w) => [0, w])), mg = lim / (1.2 * Om0) - 1, warnings = [];
    if (i.x_cg > i.x_pa + 1e-9 && !fin(lim)) warnings.push('The mass centre is aft of the pitch axis but no instability was found in this speed range.');
    pkNotes(r, warnings, 'rotor speed');
    if (fin(lim) && mg < 0) warnings.push(`Pitch–flap ${Of < Od ? 'flutter' : 'divergence'} at ${((lim * 60) / TAU).toFixed(0)} rpm is inside 120% of operating speed.`);
    if ((top * i.R) / isa(i.alt_m).a > 0.9) warnings.push('Tip Mach number exceeds 0.9 at the top of the range: compressibility and stall flutter are outside this model.');
    const pure = ctx?.case && !(ctx.case.wing.S_m2 > 0);
    return {
      kpis: [
        { key: 'rotor_flutter_rpm', label: 'Pitch–flap flutter speed', value: (Of * 60) / TAU, unit: 'rpm' },
        { key: 'rotor_div_rpm', label: 'Pitch divergence speed', value: (Od * 60) / TAU, unit: 'rpm' },
        { key: 'rotor_speed_margin', label: 'Stability margin on rotor speed', value: fin(mg) ? mg : 9.99, unit: '-', status: mg > 0 ? 'ok' : 'bad', note: 'Lowest instability speed / (1.2 · operating speed) − 1; 9.99 when none was found' },
        { key: 'flap_damping', label: 'Flap mode damping ratio at operating speed', value: -rt[0][0] / Math.hypot(rt[0][0], rt[0][1]), unit: '-' },
        { key: 'pitch_damping', label: 'Pitch mode damping ratio at operating speed', value: -rt[1][0] / Math.hypot(rt[1][0], rt[1][1]), unit: '-', status: rt[1][0] < 0 ? 'ok' : 'bad' },
        { key: 'flap_freq_rev', label: 'Flap frequency', value: rt[0][1] / Om0, unit: '/rev' }, { key: 'pitch_freq_rev', label: 'Pitch frequency', value: rt[1][1] / Om0, unit: '/rev' },
        { key: 'lock_number', label: 'Lock number', value: lock, unit: '-' },
        { key: 'cg_offset_pct', label: 'Mass centre aft of the pitch axis', value: 100 * (i.x_cg - i.x_pa), unit: '% chord', status: i.x_cg <= i.x_pa + 0.02 ? 'ok' : 'warn' },
      ].filter((k) => fin(k.value)),
      plots: [
        { type: 'line', title: 'Modal damping against rotor speed', xlabel: 'Rotor speed [rpm]', ylabel: 'Damping ratio [-]', series: r.tracks.map((tk, k) => ({ name: ['Flap branch', 'Pitch branch'][k], x: rpm, y: tk.map((l) => N.clamp(-l[0] / (Math.hypot(l[0], l[1]) || 1), -1, 1)) })), annotations: [{ y: 0, label: 'Stability boundary' }, { x: i.rpm, label: 'Operating' }] },
        { type: 'line', title: 'Modal frequency against rotor speed', xlabel: 'Rotor speed [rpm]', ylabel: 'Frequency [Hz]', series: [...r.tracks.map((tk, k) => ({ name: ['Flap branch', 'Pitch branch'][k], x: rpm, y: tk.map((l) => l[1] / TAU) })), { name: '1/rev', x: rpm, y: rpm.map((v) => v / 60), style: 'dash' }] },
      ],
      outputs: pure && fin(Of) ? { f_flutter_Hz: r.flutter.w / TAU } : {},
      warnings, models: ['Rotor aeroelastic model: rigid flap and rigid pitch with centrifugal and propeller-moment coupling', 'p–k method with Theodorsen strip aerodynamics at the local reduced frequency', PK_MODEL, 'Hover, zero inflow perturbation'],
      assumptions: ['Rigid blade: elastic bending and torsion modes, lag motion and pitch–lag/pitch–flap kinematic coupling (δ3) are not included', 'Shed wake of preceding blades and revolutions (Loewy lift deficiency) neglected', 'Hover only; forward-flight periodic coefficients are not modelled', 'Uniform blade mass and chord'],
    };
  },
  convergence: { param: 'nSpeeds', label: 'Rotor speed points', levels: [8, 15, 30, 60], metric: 'pitch_damping' },
  calibration: { params: [{ key: 'f_tors', min: 0.1, max: 1e3 }, { key: 'x_cg', min: 0.15, max: 0.5 }], sweep: 'rpm', target: 'pitch_damping', note: 'Measured pitch-mode damping from whirl-tower shake tests against rotor speed' },
  verify() {
    // pitch locked out and mass-balanced: quasi-steady flap damping is the Lock-number result γ/16 (per rev)
    const b = { R: 6, chord: 0.05, m_blade: 60, e_hinge: 0, f_flap_nr: 0, f_tors: 2000, x_pa: 0.25, x_cg: 0.25, r_alpha: 0.25, zeta: 0, cl: 5.7 }, Om = 30, m = bladeModel(b, Om, 60), rho = 1.2, l = pkPoint(m, rho, [0, Om], true);
    const Ib = 60 * 36 / 3, Ia = N.sum(m.strips.map((s) => Math.PI * rho * s.b * s.b * s.dy * s.ph[0] ** 2)), gam = (rho * 5.7 * 0.05 * 6 ** 4) / (Ib + Ia);
    const d = { ...b, chord: 0.3, m_blade: 40, R: 5, f_tors: 3, x_pa: 0.4, x_cg: 0.4 }, If = 40 * (0.25 * 0.3) ** 2, Od = Math.sqrt((If * (TAU * 3) ** 2) / ((rho * 5.7 * 0.3 * 0.045 * 5 ** 3) / 6 - If));
    const dn = N.brent((o) => { const mm = bladeModel(d, o, 80); return N.det(N.madd(mm.K, aeroStiff(mm, rho))); }, 0.5 * Od, 1.5 * Od, 1e-9);
    return [N.check('Flap damping σ = −γΩ/16', l[0], (-gam * Om) / 16, 1e-3, 'Rigid flapping blade in hover; Johnson, Helicopter Theory'),
      N.check('Static pitch divergence speed of a mass-balanced blade', dn, Od, 1e-3, 'I_f(Ω² + ω_θ²) = ρ·a·c·e·Ω²·R³/6 with e the pitch axis aft of the quarter chord'),
      N.check('Flap frequency ω² + σ² = Ω²', Math.hypot(l[0], l[1]), Om, 1e-3, 'Centrally hinged blade')];
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (o.rotor_speed_margin < 0) out.push({ severity: 'critical', title: 'Pitch–flap instability inside the rotor speed range', detail: `Margin ${(100 * o.rotor_speed_margin).toFixed(0)}% on 120% rotor speed with the mass centre ${o.cg_offset_pct.toFixed(1)}% chord aft of the pitch axis.`, action: 'Add leading-edge balance mass to bring the section mass centre to or ahead of the quarter chord, and stiffen the pitch-control path.', basis: 'Pitch–flap flutter and divergence criteria' });
    else if (o.cg_offset_pct > 2) out.push({ severity: 'advise', title: 'Blade is not mass-balanced', detail: `Mass centre is ${o.cg_offset_pct.toFixed(1)}% chord aft of the pitch axis.`, action: 'Use the stability-map analysis to see how much control stiffness this offset needs, including a failed or soft pitch link.', basis: 'Chordwise mass balance rule' });
    return out;
  },
};
const rotorMap = {
  id: 'rotormap', title: 'Blade stability map: mass balance against control stiffness', fidelity: 'reduced-order',
  summary: 'Time-domain (Wagner) eigenvalue stability of the pitch–flap blade over chordwise mass-centre position and pitch frequency, with the critical aft mass-centre limit and a disturbed-pitch time history.',
  equations: ['Wagner indicial lift formulation', 'Flutter eigenvalue equations', 'Divergence equilibrium equations', 'Coupled aerodynamic–structural equations of motion'],
  applicable: hasRotor,
  inputs: [...BLADE_INPUTS, { key: 'overspeed', label: 'Assessment speed / operating', unit: '-', default: 1.2, min: 1, max: 2, group: 'Operation' }, { key: 'theta0_deg', label: 'Initial pitch disturbance', unit: 'deg', default: 1, min: 0.01, max: 10, group: 'Initial conditions' },
    { key: 'nGrid', label: 'Map resolution', unit: '', default: 14, min: 6, max: 40, step: 1, discrete: true, group: 'Numerics' }],
  defaults: (c) => bladeDefaults(c),
  run(i) {
    const rho = isa(i.alt_m).rho, Om = ((i.rpm * TAU) / 60) * i.overspeed, ng = Math.round(i.nGrid), g = (xc, ft) => growth(bladeModel({ ...i, x_cg: xc, f_tors: ft }, Om, 6), rho);
    const xs = N.linspace(i.x_pa - 0.05, i.x_pa + 0.25, ng), fr = N.linspace(1.5, 10, ng), z = fr.map((f) => xs.map((x) => g(x, (f * Om) / TAU) / Om));
    // critical mass-centre position at the actual pitch frequency
    const gx = (x) => g(x, i.f_tors); let xcr = Infinity;
    if (gx(i.x_pa + 0.3) > 0) xcr = gx(i.x_pa - 0.05) > 0 ? -Infinity : N.brent(gx, i.x_pa - 0.05, i.x_pa + 0.3, 1e-6);
    const now = gx(i.x_cg), cgM = 100 * (xcr - i.x_cg), ss = stateSpace(bladeModel(i, Om, 6), rho), x0 = new Array(ss.nx).fill(0); x0[1] = N.rad(i.theta0_deg);
    const T = TAU / Om, sim = N.rk4((t, x) => N.matvec(ss.A, x), 0, x0, 12 * T, 1200), st = 4, warnings = [];
    if (now > 0) warnings.push(`The blade is unstable at ${(100 * i.overspeed).toFixed(0)}% speed (amplitude doubles in ${(Math.LN2 / now).toFixed(3)} s).`);
    if (!fin(xcr)) warnings.push(xcr > 0 ? 'No instability up to 30% chord aft of the pitch axis at this control stiffness.' : 'Unstable even with the mass centre ahead of the pitch axis: pitch stiffness is too low.');
    return {
      kpis: [
        { key: 'cg_crit_pct', label: 'Critical mass-centre position', value: 100 * xcr, unit: '% chord' },
        { key: 'rotor_cg_margin_pct_c', label: 'Mass-balance margin', value: fin(cgM) ? cgM : cgM > 0 ? 30 : -30, unit: '% chord', status: cgM > 2 ? 'ok' : cgM > 0 ? 'warn' : 'bad', note: 'Critical minus actual mass-centre position; capped at ±30 when no boundary lies in the range' },
        { key: 'growth_per_rev', label: 'Least stable root, real part / Ω', value: now / Om, unit: '-', status: now < 0 ? 'ok' : 'bad' },
        { key: 'pitch_freq_rev_nr', label: 'Non-rotating pitch frequency at this speed', value: (i.f_tors * TAU) / Om, unit: '/rev' },
      ].filter((k) => fin(k.value)),
      plots: [
        { type: 'heat', title: 'Stability map: least stable root', xlabel: 'Mass centre [chord fraction]', ylabel: 'Non-rotating pitch frequency [/rev]', zlabel: 'Real part / Ω [-]', x: xs, y: fr, z: z.map((r) => r.map((v) => N.clamp(v, -0.3, 0.3))), contours: 10, diverging: true, overlay: [{ name: 'This blade', x: [i.x_cg], y: [N.clamp((i.f_tors * TAU) / Om, 1.5, 10)] }] },
        { type: 'line', title: 'Response to a pitch disturbance', xlabel: 'Rotor revolutions [-]', ylabel: 'Angle [deg]', series: [{ name: 'Pitch', x: sim.t.filter((_, k) => k % st === 0).map((t) => t / T), y: sim.y.filter((_, k) => k % st === 0).map((x) => N.clamp(N.deg(x[1]), -90, 90)) }, { name: 'Flap', x: sim.t.filter((_, k) => k % st === 0).map((t) => t / T), y: sim.y.filter((_, k) => k % st === 0).map((x) => N.clamp(N.deg(x[0]), -90, 90)) }] },
      ],
      warnings, models: ['Non-linear-ready time-domain aeroelastic model: Jones two-lag Wagner states per blade strip', 'Eigenvalue stability map (covers flutter and divergence together)'],
      assumptions: ['Same rigid pitch–flap blade model as the rotor flutter analysis', 'Hover, no returning-wake effects', 'Map is drawn at the assessment (overspeed) condition'],
    };
  },
  convergence: { param: 'nGrid', label: 'Map resolution', levels: [6, 10, 14, 20], metric: 'growth_per_rev' },
  verify() {
    // the Jones time-domain eigenvalues must put the stability boundary where the Theodorsen p-k solution does
    const b = { R: 5, chord: 0.3, m_blade: 40, e_hinge: 0, f_flap_nr: 0, f_tors: 12, x_pa: 0.25, x_cg: 0.35, r_alpha: 0.25, zeta: 0, cl: 5.7 }, rho = 1.2, Os = N.linspace(2, 80, 40);
    const pk = pkSweep((o) => bladeModel(b, o), rho, Os, bladeModel(b, 2).wn), lim = Math.min(pk.flutter?.P ?? Infinity, fin(pk.div) ? pk.div : Infinity);
    const td = N.findRoot((o) => growth(bladeModel(b, o), rho), 2, 80, 78, 1e-7);
    return [N.check('Time-domain stability boundary matches the p-k boundary', td, lim, 0.03, 'Jones (1940) two-lag approximation against exact Theodorsen aerodynamics')];
  },
  recommend(res) {
    const o = res.outputs;
    return o.rotor_cg_margin_pct_c < 2 ? [{ severity: o.rotor_cg_margin_pct_c < 0 ? 'critical' : 'warn', title: 'Little or no mass-balance margin', detail: `Margin is ${o.rotor_cg_margin_pct_c.toFixed(1)}% chord at the assessment speed.`, action: 'Move the section mass centre forward with balance weights or raise the pitch-control stiffness; repeat with degraded control stiffness.', basis: 'Eigenvalue stability map of the pitch–flap blade' }] : [];
  },
};

// ---- 8. unsteady vortex lattice coupled to a finite-element beam wing ---------------------------
/** Vertical velocity at (px, py) from a unit vortex ring with corners x[4], y[4] lying in the plane z = 0. */
function ringW(px, py, x, y) {
  let w = 0;
  for (let k = 0; k < 4; k++) {
    const m = (k + 1) & 3, ax = px - x[k], ay = py - y[k], bx = px - x[m], by = py - y[m], cr = ax * by - ay * bx, ra = Math.sqrt(ax * ax + ay * ay), rb = Math.sqrt(bx * bx + by * by);
    if (cr * cr < 1e-20 * ra * ra * rb * rb || ra * rb === 0) continue;
    const dx = x[m] - x[k], dy = y[m] - y[k];
    w += ((dx * ax + dy * ay) / ra - (dx * bx + dy * by) / rb) / cr;
  }
  return w / (4 * Math.PI);
}
/**
 * Semi-span vortex-ring lattice of the swept, tapered wing with its mirror image, and a flat wake of nw rows of length dxw
 * convected with the free stream (the last row is extended far downstream). Panel index k = j·nc + i (i chordwise).
 */
function uvlmLattice(i, nc, ns, nw, dxw) {
  const sw = N.rad(i.sweep_deg || 0), ys = i.semi_span * Math.cos(sw), tn = Math.tan(sw), ch = (y) => i.c_root * (1 - (1 - i.taper) * (y / ys)), P = (xi, y) => i.x_ea * i.c_root + y * tn + (xi - i.x_ea) * ch(y);
  const n = nc * ns, yj = N.linspace(0, ys, ns + 1), rings = [], g = { n, nc, ns, nw, dxw, ys, xc: new Float64Array(n), yc: new Float64Array(n), xq: new Float64Array(n), xm: new Float64Array(n), dy: new Float64Array(n), area: new Float64Array(n), chord: new Float64Array(ns) };
  const rear = (ic, y) => (ic < nc - 1 ? P((ic + 1.25) / nc, y) : P(1, y) + 0.25 * dxw);
  for (let j = 0; j < ns; j++) for (let ic = 0; ic < nc; ic++) {
    const k = j * nc + ic, y0 = yj[j], y1 = yj[j + 1], ym = 0.5 * (y0 + y1), f = (ic + 0.25) / nc;
    rings.push([[P(f, y0), P(f, y1), rear(ic, y1), rear(ic, y0)], [y0, y1, y1, y0]]);
    g.xc[k] = P((ic + 0.75) / nc, ym); g.yc[k] = ym; g.xq[k] = P(f, ym); g.xm[k] = P((ic + 0.5) / nc, ym); g.dy[k] = y1 - y0; g.area[k] = (ch(ym) / nc) * (y1 - y0); g.chord[j] = ch(ym);
  }
  const A = N.zeros(n), Aw = new Float64Array(n * nw * ns), far = 200 * ys + 50 * i.c_root;
  for (let k = 0; k < n; k++) {
    for (let l = 0; l < n; l++) A[k][l] = ringW(g.xc[k], g.yc[k], rings[l][0], rings[l][1]) + ringW(g.xc[k], -g.yc[k], rings[l][0], rings[l][1]);
    for (let r = 0; r < nw; r++) for (let j = 0; j < ns; j++) {
      const y0 = yj[j], y1 = yj[j + 1], t0 = rear(nc - 1, y0), t1 = rear(nc - 1, y1), b = r === nw - 1 ? far : (r + 1) * dxw, x = [t0 + r * dxw, t1 + r * dxw, t1 + b, t0 + b], y = [y0, y1, y1, y0];
      Aw[k * nw * ns + r * ns + j] = ringW(g.xc[k], g.yc[k], x, y) + ringW(g.xc[k], -g.yc[k], x, y);
    }
  }
  // steady matrix: every wake row carries the trailing-edge circulation of its column
  const As = A.map((r) => r.slice());
  for (let k = 0; k < n; k++) for (let j = 0; j < ns; j++) { let a = 0; for (let r = 0; r < nw; r++) a += Aw[k * nw * ns + r * ns + j]; As[k][j * nc + nc - 1] += a; }
  const flat = (M) => { const F = new Float64Array(n * n); for (let k = 0; k < n; k++) F.set(M[k], k * n); return F; };
  return { ...g, Ai: flat(N.inv(A)), Asi: flat(N.inv(As)), Aw, S_semi: N.sum(g.area) };
}
/**
 * Finite-element beam along the elastic axis: Hermite bending + linear torsion, three freedoms per node (w up, dw/ds, twist
 * nose-up), consistent mass with the inertial bending–torsion coupling of the section mass offset. Clamped at the root.
 */
function beamFE(i, ne) {
  const L = i.semi_span, sw = N.rad(i.sweep_deg || 0), cs = Math.cos(sw), le = L / ne, nd = 3 * (ne + 1), K = N.zeros(nd), M = N.zeros(nd), m0 = (3 * i.m_semi) / (L * (1 + i.taper + i.taper ** 2));
  const gp = [0.5 - Math.sqrt(0.15), 0.5, 0.5 + Math.sqrt(0.15)], gw = [5 / 18, 8 / 18, 5 / 18], EI = [], GJ = [];
  for (let e = 0; e < ne; e++) {
    const r = 1 - (1 - i.taper) * ((e + 0.5) / ne), chn = i.c_root * r * cs, ei = i.EI_root * r ** i.stiff_exp, gj = i.GJ_root * r ** i.stiff_exp, m = m0 * r * r, Ia = m * (i.r_alpha * chn) ** 2, S = m * (i.x_cg - i.x_ea) * chn;
    EI.push(ei); GJ.push(gj);
    const b = [3 * e, 3 * e + 1, 3 * e + 3, 3 * e + 4], t = [3 * e + 2, 3 * e + 5], l = le;
    const kb = [[12, 6 * l, -12, 6 * l], [6 * l, 4 * l * l, -6 * l, 2 * l * l], [-12, -6 * l, 12, -6 * l], [6 * l, 2 * l * l, -6 * l, 4 * l * l]], mb = [[156, 22 * l, 54, -13 * l], [22 * l, 4 * l * l, 13 * l, -3 * l * l], [54, 13 * l, 156, -22 * l], [-13 * l, -3 * l * l, -22 * l, 4 * l * l]];
    for (let p = 0; p < 4; p++) for (let q = 0; q < 4; q++) { K[b[p]][b[q]] += (ei / l ** 3) * kb[p][q]; M[b[p]][b[q]] += ((m * l) / 420) * mb[p][q]; }
    for (let p = 0; p < 2; p++) for (let q = 0; q < 2; q++) { K[t[p]][t[q]] += (gj / l) * (p === q ? 1 : -1); M[t[p]][t[q]] += ((Ia * l) / 6) * (p === q ? 2 : 1); }
    // z of the mass centre = w − d·φ: coupling −S ∫ N_w N_φ ds by three-point Gauss quadrature
    for (let g = 0; g < 3; g++) { const x = gp[g], Nw = hermite(x, l), Nt = [1 - x, x]; for (let p = 0; p < 4; p++) for (let q = 0; q < 2; q++) { const v = -S * gw[g] * l * Nw[p] * Nt[q]; M[b[p]][t[q]] += v; M[t[q]][b[p]] += v; } }
  }
  const free = N.range(nd - 3, (k) => k + 3), Kf = free.map((p) => free.map((q) => K[p][q])), Mf = free.map((p) => free.map((q) => M[p][q]));
  return { L, ne, le, nd, cs, sn: Math.sin(sw), Kf, Mf, EI, GJ, x0: i.x_ea * i.c_root };
}
const hermite = (x, l) => [1 - 3 * x * x + 2 * x ** 3, l * (x - 2 * x * x + x ** 3), 3 * x * x - 2 * x ** 3, l * (x ** 3 - x * x)];
const dHermite = (x, l) => [(-6 * x + 6 * x * x) / l, 1 - 4 * x + 3 * x * x, (6 * x - 6 * x * x) / l, 3 * x * x - 2 * x];
/** Interpolation rows over the free beam freedoms for a planform point: vertical displacement z and streamwise slope dz/dx. */
function beamRows(fe, x, y) {
  const s = N.clamp((x - fe.x0) * fe.sn + y * fe.cs, 0, fe.L), d = (x - fe.x0) * fe.cs - y * fe.sn, e = Math.min(fe.ne - 1, Math.floor(s / fe.le)), xi = s / fe.le - e, Nw = hermite(xi, fe.le), dN = dHermite(xi, fe.le);
  const z = new Float64Array(fe.nd), sx = new Float64Array(fe.nd), b = [3 * e, 3 * e + 1, 3 * e + 3, 3 * e + 4];
  for (let p = 0; p < 4; p++) { z[b[p]] = Nw[p]; sx[b[p]] = fe.sn * dN[p]; }
  z[3 * e + 2] = -d * (1 - xi); z[3 * e + 5] = -d * xi;
  sx[3 * e + 2] = -fe.cs * (1 - xi) + (d * fe.sn) / fe.le; sx[3 * e + 5] = -fe.cs * xi - (d * fe.sn) / fe.le;
  return { z: z.subarray(3), sx: sx.subarray(3), zFull: z };
}
/** Coupled model: lattice, beam modes and the conservative transfer operators projected on the retained modes. */
function fsiModel(i, o) {
  const nc = o.nc, ns = o.ns, sw = N.rad(i.sweep_deg || 0), ys = i.semi_span * Math.cos(sw), cref = (i.c_root * (1 + i.taper)) / 2, dxw = o.dtau * cref, lat = uvlmLattice(i, nc, ns, o.nw, dxw), fe = beamFE(i, o.ne), nf = fe.nd - 3;
  const sK = fe.Kf[0][0], sM = fe.Mf[0][0], eg = N.eigGenSym(N.mscale(fe.Kf, 1 / sK), N.mscale(fe.Mf, 1 / sM)), nm = Math.min(o.nm, nf), wn = [], Phi = [];
  for (let m = 0; m < nm; m++) { const v = eg.vectors[m], mv = N.matvec(fe.Mf, v), g = Math.sqrt(N.dot(v, mv)); Phi.push(v.map((x) => x / g)); wn.push(Math.sqrt(Math.max((eg.values[m] * sK) / sM, 0))); }
  const n = lat.n, Zc = new Float64Array(n * nm), Xc = new Float64Array(n * nm), Zq = new Float64Array(n * nm), Zm = new Float64Array(n * nm), rowsQ = [];
  for (let k = 0; k < n; k++) {
    const rc = beamRows(fe, lat.xc[k], lat.yc[k]), rq = beamRows(fe, lat.xq[k], lat.yc[k]), rm = beamRows(fe, lat.xm[k], lat.yc[k]); rowsQ.push(rq);
    for (let m = 0; m < nm; m++) { let a = 0, b = 0, c = 0, d = 0; const p = Phi[m]; for (let q = 0; q < nf; q++) { a += rc.z[q] * p[q]; b += rc.sx[q] * p[q]; c += rq.z[q] * p[q]; d += rm.z[q] * p[q]; } Zc[k * nm + m] = a; Xc[k * nm + m] = b; Zq[k * nm + m] = c; Zm[k * nm + m] = d; }
  }
  // output rows on the modes: tip deflection, tip twist and root bending moment EI·w''(0)
  const tipW = Phi.map((p) => p[nf - 3]), tipT = Phi.map((p) => p[nf - 1]), l = fe.le, rootM = Phi.map((p) => (fe.EI[0] * (6 * p[0] - 2 * l * p[1])) / (l * l));
  const tors = Phi.map((p) => { let t = 0, tot = 0; for (let a = 0; a < nf; a++) for (let b = 0; b < nf; b++) { const e = p[a] * fe.Mf[a][b] * p[b]; tot += e; if (a % 3 === 2 && b % 3 === 2) t += e; } return N.clamp(t / tot, 0, 1); });
  return { lat, fe, nm, nf, wn, Phi, Zc, Xc, Zq, Zm, rowsQ, tipW, tipT, rootM, tors, cref, dxw, ys, zeta: i.zeta, inp: i, opt: o, coarse: {} };
}
/**
 * Strongly coupled time march. Each step: convect the wake, then sub-iterate lattice solution ↔ Newmark structure with
 * Aitken relaxation of the modal accelerations until the interface residual is below tol.
 * o: { rho, U, nSteps, alpha0, gust(x, t) → upwash, v0[], rigid, maxIt, tol, hist }.
 */
function fsiMarch(md, o) {
  const { lat, nm } = md, { n, nc, ns, nw } = lat, nwk = nw * ns, U = o.U, rho = o.rho, dt = md.dxw / U, rigid = o.rigid || nm === 0, maxIt = rigid ? 1 : o.maxIt, zeta = md.zeta;
  const gam = new Float64Array(n), gamN = new Float64Array(n), gw = new Float64Array(nwk), base = new Float64Array(n), rhs = new Float64Array(n), q = new Float64Array(nm), v = new Float64Array(nm), a = new Float64Array(nm), ak = new Float64Array(nm), vk = new Float64Array(nm), qk = new Float64Array(nm), Q = new Float64Array(nm), r = new Float64Array(nm), rp = new Float64Array(nm);
  if (o.v0) v.set(o.v0);
  const al0 = o.alpha0 || 0, gust = o.gust || null;
  const modes = o.hist === 'modes', H = { t: [], tip: [], twist: [], lift: [], rootM: [], E: modes ? new Float64Array(o.nSteps) : [], q: modes ? new Float64Array(o.nSteps * nm) : null }; let its = 0, wA = 0, wS = 0, omg = 1, lift = 0, resMax = 0, liftPk = 0;
  const Zc = md.Zc, Xc = md.Xc, Zq = md.Zq, Zm = md.Zm, Ai = lat.Ai, Aw = lat.Aw, dyp = lat.dy, area = lat.area, wn = md.wn, xcol = lat.xc, rdt = rho / dt, rU = rho * U;
  /** Normal wash → circulation → panel loads → generalised forces. Returns the total lift; acc also accumulates the interface work. */
  const loads = (acc) => {
    for (let k = 0, off = 0; k < n; k++, off += nm) { let b = base[k]; for (let m = 0; m < nm; m++) b += Zc[off + m] * vk[m] + U * Xc[off + m] * qk[m]; rhs[k] = b; }
    for (let k = 0, off = 0; k < n; k++, off += n) { let g0 = 0, g1 = 0, l = 0; for (; l + 1 < n; l += 2) { g0 += Ai[off + l] * rhs[l]; g1 += Ai[off + l + 1] * rhs[l + 1]; } if (l < n) g0 += Ai[off + l] * rhs[l]; gam[k] = g0 + g1; }
    for (let m = 0; m < nm; m++) Q[m] = 0;
    let L = 0, w = 0;
    for (let k = 0, off = 0; k < n; k++, off += nm) {
      const fs = rU * (gam[k] - (k % nc ? gam[k - 1] : 0)) * dyp[k], fu = rdt * area[k] * (gam[k] - gamN[k]); L += fs + fu;
      if (acc) { let zq = 0, zm = 0; for (let m = 0; m < nm; m++) { const a1 = Zq[off + m], a2 = Zm[off + m]; Q[m] += fs * a1 + fu * a2; zq += a1 * vk[m]; zm += a2 * vk[m]; } w += fs * zq + fu * zm; }
      else for (let m = 0; m < nm; m++) Q[m] += fs * Zq[off + m] + fu * Zm[off + m];
    }
    if (acc) { wA += w * dt; let sw = 0; for (let m = 0; m < nm; m++) sw += Q[m] * vk[m]; wS += sw * dt; }
    return L;
  };
  for (let st = 1; st <= o.nSteps; st++) {
    const t = st * dt;
    gw.copyWithin(ns, 0, nwk - ns); for (let j = 0; j < ns; j++) gw[j] = gamN[j * nc + nc - 1];
    for (let k = 0, off = 0; k < n; k++, off += nwk) { let b0 = 0, b1 = 0, l = 0; for (; l + 1 < nwk; l += 2) { b0 += Aw[off + l] * gw[l]; b1 += Aw[off + l + 1] * gw[l + 1]; } if (l < nwk) b0 += Aw[off + l] * gw[l]; base[k] = -U * al0 - (gust ? gust(xcol[k], t) : 0) - b0 - b1; }
    ak.set(a); let it = 0, res = 0;
    for (; it < maxIt; it++) {
      for (let m = 0; m < nm; m++) { vk[m] = v[m] + 0.5 * dt * (a[m] + ak[m]); qk[m] = q[m] + dt * v[m] + 0.25 * dt * dt * (a[m] + ak[m]); }
      lift = loads(false);
      if (rigid) break;
      let rr = 0, aa = 0, num = 0, den = 0;
      for (let m = 0; m < nm; m++) {
        const w = wn[m], an = (Q[m] - 2 * zeta * w * (v[m] + 0.5 * dt * a[m]) - w * w * (q[m] + dt * v[m] + 0.25 * dt * dt * a[m])) / (1 + zeta * w * dt + 0.25 * w * w * dt * dt);
        r[m] = an - ak[m]; rr += r[m] * r[m]; aa += an * an; const dr = r[m] - rp[m]; num += rp[m] * dr; den += dr * dr;
      }
      omg = it === 0 || den === 0 ? Math.min(omg, 1) : N.clamp((-omg * num) / den, 0.05, 1.5);
      for (let m = 0; m < nm; m++) { ak[m] += omg * r[m]; rp[m] = r[m]; }
      res = Math.sqrt(rr / (aa + 1e-300));
      if (res < o.tol) { it++; break; }
    }
    its += it; if (res > resMax) resMax = res;
    if (!rigid) { for (let m = 0; m < nm; m++) { vk[m] = v[m] + 0.5 * dt * (a[m] + ak[m]); qk[m] = q[m] + dt * v[m] + 0.25 * dt * dt * (a[m] + ak[m]); } lift = loads(true); for (let m = 0; m < nm; m++) { q[m] = qk[m]; v[m] = vk[m]; a[m] = ak[m]; } }
    gamN.set(gam); if (Math.abs(lift) > Math.abs(liftPk)) liftPk = lift;
    if (modes) { let E = 0; for (let m = 0; m < nm; m++) { E += 0.5 * (v[m] * v[m] + md.wn[m] ** 2 * q[m] * q[m]); H.q[(st - 1) * nm + m] = q[m]; } H.E[st - 1] = E; }
    else if (o.hist) {
      let tip = 0, tw = 0, rm = 0, E = 0; for (let m = 0; m < nm; m++) { tip += md.tipW[m] * q[m]; tw += md.tipT[m] * q[m]; rm += md.rootM[m] * q[m]; E += 0.5 * (v[m] * v[m] + md.wn[m] ** 2 * q[m] * q[m]); }
      H.t.push(t); H.tip.push(tip); H.twist.push(tw); H.lift.push(lift); H.rootM.push(rm); H.E.push(E);
    }
  }
  return { H, dt, its: its / o.nSteps, resMax, wA, wS, q: Array.from(q), liftPk, gam: Array.from(gam) };
}
/** Static aeroelastic equilibrium at incidence alpha by partitioned fixed-point iteration with Aitken relaxation. */
function fsiStatic(md, rho, U, alpha, rigid) {
  const { lat, nm } = md, { n, nc } = lat, q = new Array(nm).fill(0), rhs = new Float64Array(n), gam = new Float64Array(n); let rp = null, omg = 1, it = 0, res = Infinity, L = 0, Fk = [];
  for (; it < 200; it++) {
    for (let k = 0; k < n; k++) { let b = -U * alpha; for (let m = 0; m < nm; m++) b += U * md.Xc[k * nm + m] * q[m]; rhs[k] = b; }
    for (let k = 0; k < n; k++) { let g = 0; for (let l = 0; l < n; l++) g += lat.Asi[k * n + l] * rhs[l]; gam[k] = g; }
    const Q = new Array(nm).fill(0); L = 0; Fk = [];
    for (let k = 0; k < n; k++) { const f = rho * U * (gam[k] - (k % nc ? gam[k - 1] : 0)) * lat.dy[k]; L += f; Fk.push(f); for (let m = 0; m < nm; m++) Q[m] += f * md.Zq[k * nm + m]; }
    if (rigid || nm === 0) { res = 0; break; }
    const r = Q.map((x, m) => x / md.wn[m] ** 2 - q[m]);
    if (rp) { const dr = r.map((x, m) => x - rp[m]), den = N.dot(dr, dr); omg = den > 0 ? N.clamp((-omg * N.dot(rp, dr)) / den, 0.02, 1.5) : omg; }
    for (let m = 0; m < nm; m++) q[m] += omg * r[m]; rp = r;
    res = N.norm(r) / (N.norm(q) + 1e-300);
    if (!(res < 1e12)) break;
    if (res < 1e-10) { it++; break; }
  }
  return { q, L, Fk, gam: Array.from(gam), it, converged: res < 1e-8, tip: N.dot(md.tipW, q), twist: N.dot(md.tipT, q), rootM: N.dot(md.rootM, q) };
}
/**
 * Free-decay run after a modal velocity kick: exponential growth rate σ of the energy envelope and the oscillation frequency.
 * When the record would need more than o.maxSteps steps, the time step (and with it the wake row length) is doubled as needed.
 */
function fsiGrowth(md0, rho, U, nPer, o) {
  const T1 = TAU / md0.wn[0]; let mult = 1; while ((nPer * T1 * U) / (md0.dxw * mult) > o.maxSteps && mult < 64) mult *= 2;
  const md = mult === 1 ? md0 : (md0.coarse[mult] = md0.coarse[mult] || { ...md0, dxw: md0.dxw * mult, lat: uvlmLattice(md0.inp, md0.opt.nc, md0.opt.ns, md0.opt.nw, md0.dxw * mult) });
  const dt = md.dxw / U, nm = md.nm, nSteps = Math.max(40, Math.round((nPer * T1) / dt)), v0 = md.wn.map((w, m) => (m < 4 ? 1 : 0));
  const r = fsiMarch(md, { rho, U, nSteps, v0, maxIt: o.maxIt, tol: o.tol, hist: 'modes' }), E = r.H.E, Q = r.H.q, k0 = Math.floor(0.35 * nSteps);
  if (!fin(E[nSteps - 1])) return { sigma: Infinity, w: NaN, its: r.its, nSteps, dt, blown: true };
  // least-squares slope of ln E over the last 65% of the record
  let sx = 0, sy = 0, sxx = 0, sxy = 0, n = 0; for (let k = k0; k < nSteps; k++) { const x = (k + 1) * dt, y = Math.log(Math.max(E[k], 1e-300)); sx += x; sy += y; sxx += x * x; sxy += x * y; n++; }
  const sigma = (n * sxy - sx * sy) / (n * sxx - sx * sx) / 2;
  // frequency of the mode carrying most energy at the end of the record, from its zero crossings
  let mb = 0, eb = -1; for (let m = 0; m < nm; m++) { let e = 0; for (let k = Math.floor(0.8 * nSteps); k < nSteps; k++) e += (md.wn[m] * Q[k * nm + m]) ** 2; if (e > eb) { eb = e; mb = m; } }
  let first = 0, last = 0, nz = 0; for (let k = k0 + 1; k < nSteps; k++) { const y0 = Q[(k - 1) * nm + mb], y1 = Q[k * nm + mb]; if (y0 < 0 && y1 >= 0) { const t = (k + (-y0) / (y1 - y0)) * dt; if (!nz) first = t; last = t; nz++; } }
  return { sigma, w: nz > 1 ? (TAU * (nz - 1)) / (last - first) : 0, its: r.its, nSteps, dt, mode: mb, blown: false };
}
const UVLM_BASE = { semi_span: 6, c_root: 1, taper: 1, sweep_deg: 0, x_ea: 0.4, x_cg: 0.45, r_alpha: 0.25, EI_root: 2e5, GJ_root: 1e5, stiff_exp: 0, m_semi: 60, zeta: 0.02, cl: 5.5, alt_m: 0, V_d_eas: 60, V: 40, mass_kg: 600, n_load: 1, V_flutter_pk: 0, f_flutter_pk: 0, w_gust: 5, H_chords: 12.5, nChord: 4, nSpan: 8, nWake: 16, dtau: 0.25, nElem: 8, nModes: 6, nSpeeds: 5, nPeriods: 4, stepsPerPeriod: 48, maxIt: 10, tol: 1e-6 };
const fsiOpts = (i) => ({ nc: Math.max(1, Math.round(i.nChord)), ns: Math.max(2, Math.round(i.nSpan)), nw: Math.max(4, Math.round(i.nWake)), dtau: i.dtau, ne: Math.max(2, Math.round(i.nElem)), nm: Math.max(1, Math.round(i.nModes)) });
const uvlmfsi = {
  id: 'uvlmfsi', title: 'Coupled unsteady vortex-lattice / beam wing (time-domain fluid–structure interaction)', fidelity: 'numerical',
  summary: 'An unsteady vortex-ring lattice with a shed wake is strongly coupled, step by step, to a finite-element bending–torsion beam of the wing. It gives the static aeroelastic shape in flight, the response to a 1-cosine gust and the flutter speed found from growth or decay of free oscillations, each compared with the strip-theory analyses of this suite.',
  equations: ['Unsteady aerodynamic integral equations', 'Coupled aerodynamic–structural equations of motion', 'Linearised structural dynamic equations', 'Wagner indicial lift formulation', 'Küssner gust-response formulation', 'Modal superposition equations', 'Interface displacement compatibility equations', 'Interface traction equilibrium equations', 'Strongly coupled and loosely coupled time-marching formulations'],
  applicable: hasWing,
  inputs: [...WING_INPUTS,
    { key: 'V', label: 'Flight speed (TAS)', unit: 'm/s', default: 230, min: 1, group: 'Flight' },
    { key: 'mass_kg', label: 'Aircraft mass', unit: 'kg', default: 70000, min: 0.1, group: 'Flight' },
    { key: 'n_load', label: 'Load factor for the static shape', unit: 'g', default: 1, min: 0.1, max: 9, group: 'Flight' },
    { key: 'V_flutter_pk', label: 'p-k flutter speed for comparison (TAS)', unit: 'm/s', default: 0, min: 0, group: 'Flight', help: 'From the multi-mode wing flutter analysis of this suite when it has been run; 0 recomputes it here with the same settings (3 modes per motion, 30 speeds)' },
    { key: 'f_flutter_pk', label: 'p-k flutter frequency for comparison', unit: 'Hz', default: 0, min: 0, group: 'Flight' },
    { key: 'w_gust', label: 'Gust velocity (TAS)', unit: 'm/s', default: 15.24, min: 0, max: 60, group: 'Gust', help: 'Peak of the 1-cosine gust; 15.24 m/s (50 ft/s) is the classical derived gust velocity' },
    { key: 'H_chords', label: 'Gust gradient distance', unit: 'chords', default: 12.5, min: 1, max: 200, group: 'Gust', help: 'Distance to the gust peak in mean chords' },
    { key: 'nChord', label: 'Chordwise panels', unit: '', default: 4, min: 1, max: 16, step: 1, discrete: true, group: 'Numerics', help: 'Four resolve lift and pitching moment to a few per cent; use 8 or more for final values' },
    { key: 'nSpan', label: 'Spanwise panels (semi-span)', unit: '', default: 8, min: 2, max: 40, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nWake', label: 'Wake rows kept', unit: '', default: 12, min: 4, max: 400, step: 1, discrete: true, group: 'Numerics', help: 'Rows of shed vortex rings; the last row extends far downstream' },
    { key: 'dtau', label: 'Time step U·Δt / mean chord', unit: '-', default: 0.25, min: 0.02, max: 1, group: 'Numerics', help: 'Also the length of one wake row. About one panel chord (1 / chordwise panels) is the usual choice' },
    { key: 'nElem', label: 'Beam finite elements', unit: '', default: 8, min: 2, max: 40, step: 1, discrete: true, group: 'Numerics' },
    { key: 'nModes', label: 'Beam modes retained', unit: '', default: 6, min: 1, max: 30, step: 1, discrete: true, group: 'Numerics', help: 'The finite-element model is reduced to its lowest natural modes before time integration' },
    { key: 'nSpeeds', label: 'Speeds in the flutter search', unit: '', default: 5, min: 0, max: 40, step: 1, discrete: true, group: 'Numerics', help: '0 skips the time-domain flutter search' },
    { key: 'nPeriods', label: 'First-bending periods simulated per speed', unit: '-', default: 4, min: 2, max: 40, group: 'Numerics' },
    { key: 'stepsPerPeriod', label: 'Largest number of steps per first-bending period in the flutter search', unit: '', default: 36, min: 20, max: 2000, step: 1, discrete: true, group: 'Numerics', help: 'When the time step above would need more, it is doubled (with the wake row length) for the flutter runs; raise this for a finer flutter result' },
    { key: 'maxIt', label: 'Maximum coupling sub-iterations per step', unit: '', default: 10, min: 1, max: 50, step: 1, discrete: true, group: 'Numerics', help: '1 gives a loosely coupled (staggered) scheme' },
    { key: 'tol', label: 'Coupling tolerance', unit: '-', default: 1e-5, min: 1e-12, max: 1e-2, group: 'Numerics', help: 'Relative change of the structural accelerations between sub-iterations; 1e-6 or tighter for final values' },
  ],
  defaults: (c, up, d) => ({ ...wingDefaults(c, up, d), n_load: 1, V_flutter_pk: up.aeroelastic?.V_flutter_ms > 0 ? up.aeroelastic.V_flutter_ms : undefined, f_flutter_pk: up.aeroelastic?.f_flutter_Hz > 0 ? up.aeroelastic.f_flutter_Hz : undefined }),
  run(i, ctx) {
    const at = isa(i.alt_m), rho = at.rho, sg = Math.sqrt(at.sigma), o = fsiOpts(i), md = fsiModel(i, o), warnings = [], W = (i.mass_kg * G0 * i.n_load) / 2, cref = md.cref;
    ctx?.progress?.(0.05, 'Static aeroelastic equilibrium');
    // --- static shape at the flight condition (linear: solve for unit incidence, scale to the required lift)
    const fl = fsiStatic(md, rho, i.V, 1, false), rg = fsiStatic(md, rho, i.V, 1, true), qS = 0.5 * rho * i.V ** 2 * md.lat.S_semi, CLa_r = rg.L / qS, CLa_f = fl.L / qS, aTrim = fl.converged && fl.L > 0 ? W / fl.L : NaN;
    if (!fl.converged) warnings.push('The static coupling iteration did not converge: the flight speed is at or above the static divergence speed of this model.');
    // conservation of the load transfer: virtual work on both sides of the interface and the transferred resultant
    const fe = md.fe, fn = new Float64Array(fe.nd); fl.Fk.forEach((f, k) => { const z = md.rowsQ[k].zFull; for (let p = 0; p < fe.nd; p++) fn[p] += f * z[p]; });
    let sumNode = 0; for (let p = 0; p < fe.nd; p += 3) sumNode += fn[p];
    const du = md.Phi[0].map((_, p) => N.sum(md.Phi.map((ph, m) => ph[p] * (m + 1) * 0.37))), wAero = N.sum(fl.Fk.map((f, k) => f * N.dot(Array.from(md.rowsQ[k].z), du))), wStr = N.dot(Array.from(fn.subarray(3)), du), vwErr = Math.abs(wAero - wStr) / (Math.abs(wAero) + 1e-300), resErr = Math.abs(sumNode - fl.L) / (Math.abs(fl.L) + 1e-300);
    // reference: strip-theory static analysis of this suite at the same condition
    let ref = {}; try { ref = N.kv(wingStatic.run({ ...wingPick(i), V: i.V, mass_kg: i.mass_kg * i.n_load, ail_in: 0.7, ail_out: 0.95, aileron_frac: 0.25, k_visc: 0.8, nModes: 4 })); } catch { ref = {}; }
    ctx?.progress?.(0.2, 'Gust response');
    // --- 1-cosine gust, flexible and rigid
    const H = i.H_chords * cref, x0 = -0.5 * cref, gust = (x, t) => { const s = i.V * t - (x - x0); return s > 0 && s < 2 * H ? 0.5 * i.w_gust * (1 - Math.cos((Math.PI * s) / H)) : 0; }, dt = md.dxw / i.V;
    const nG = Math.min(4000, Math.ceil((2 * H + 10 * cref + 2 * i.c_root + md.ys * Math.abs(Math.tan(N.rad(i.sweep_deg)))) / md.dxw)), gf = fsiMarch(md, { rho, U: i.V, nSteps: nG, gust, maxIt: o.nm ? Math.round(i.maxIt) : 1, tol: i.tol, hist: true }), gr = fsiMarch(md, { rho, U: i.V, nSteps: nG, gust, rigid: true, hist: true });
    const pk = (a) => a.reduce((m, x) => (Math.abs(x) > Math.abs(m) ? x : m), 0), tipPk = pk(gf.H.tip), rmPk = pk(gf.H.rootM), dnF = (2 * pk(gf.H.lift)) / (i.mass_kg * G0), dnR = (2 * pk(gr.H.lift)) / (i.mass_kg * G0), rmStatic = fl.rootM * aTrim;
    const blown = !fin(tipPk) || Math.abs(tipPk) > 50 * i.semi_span;
    if (blown) warnings.push('The gust response grows without bound: the flight speed is above the flutter or divergence speed of the coupled model.');
    if (gf.resMax > 100 * i.tol && Math.round(i.maxIt) > 1) warnings.push(`The coupling sub-iterations hit the limit of ${Math.round(i.maxIt)} with a residual of ${gf.resMax.toExponential(1)}; raise the limit.`);
    // --- flutter onset from time-domain growth or decay
    const nV = Math.round(i.nSpeeds), top0 = Math.max(2.2 * i.V_d_eas / sg, 1.5 * at.a); let pkr = null; if (!(i.V_flutter_pk > 0)) { try { pkr = wingFlutterSolve(i, rho, 3, 30, top0); if (!pkr.flutter) pkr = wingFlutterSolve(i, rho, 3, 30, 3 * top0); } catch { pkr = null; } }
    const VfPk = i.V_flutter_pk > 0 ? i.V_flutter_pk : pkr?.flutter?.P ?? NaN, wfPk = i.V_flutter_pk > 0 ? (i.f_flutter_pk > 0 ? TAU * i.f_flutter_pk : NaN) : pkr?.flutter?.w ?? NaN, gro = { maxIt: Math.round(i.maxIt), tol: i.tol, maxSteps: Math.round(i.stepsPerPeriod) * i.nPeriods }, Vs = [], sig = [], frq = [], dec = [];
    let Vf = Infinity, wf = NaN, itsF = 0, nRun = 0, stepsPer = 0, dtF = NaN;
    if (nV >= 2 && md.nm >= 2) {
      const hi = fin(VfPk) ? Math.min(1.6 * VfPk, 2 * top0) : top0, lo = 0.25 * hi, run1 = (V) => { const g = fsiGrowth(md, rho, V, i.nPeriods, gro); itsF += g.its; nRun++; stepsPer = Math.max(stepsPer, g.nSteps); return g; };
      let prev = null, br = null;
      for (let k = 0; k < nV; k++) {
        const V = lo + ((hi - lo) * k) / (nV - 1), g = run1(V); Vs.push(V); sig.push(g.sigma); frq.push(g.w / TAU); dec.push(g.w > 0 && fin(g.sigma) ? (-TAU * g.sigma) / g.w : NaN);
        if (!br && prev && prev.g.sigma <= 0 && g.sigma > 0) br = [prev, { V, g }];
        prev = { V, g }; ctx?.progress?.(0.3 + (0.6 * (k + 1)) / (nV + 3), 'Flutter search');
        if (br) { if (k < nV - 1 && g.blown) break; }
      }
      if (br) {
        let [a, b] = br;
        for (let it = 0; it < 2; it++) { const sb = fin(b.g.sigma) ? b.g.sigma : 10 * Math.abs(a.g.sigma) + 1, Vm = N.clamp(a.V - (a.g.sigma * (b.V - a.V)) / (sb - a.g.sigma), a.V + 0.1 * (b.V - a.V), b.V - 0.1 * (b.V - a.V)), g = run1(Vm); if (g.sigma > 0) b = { V: Vm, g }; else a = { V: Vm, g }; }
        const sb = fin(b.g.sigma) ? b.g.sigma : NaN; Vf = fin(sb) ? a.V - (a.g.sigma * (b.V - a.V)) / (sb - a.g.sigma) : 0.5 * (a.V + b.V); wf = b.g.w > 0 ? b.g.w : a.g.w; dtF = b.g.dt;
      } else if (sig[0] > 0) { Vf = Vs[0]; wf = TAU * frq[0]; dtF = md.dxw / Vf; warnings.push('The coupled model is already unstable at the lowest speed of the search; the flutter speed lies below it.'); }
      else warnings.push('No instability was found in the time-domain search range (up to ' + hi.toFixed(0) + ' m/s).');
    } else warnings.push('The time-domain flutter search was skipped (it needs at least two speeds and two beam modes).');
    const isDiv = fin(Vf) && !(wf > 0.2 * md.wn[0]);
    if (isDiv) warnings.push('The instability found in the time domain is non-oscillatory: it is static divergence rather than flutter.');
    const stepsPerCycle = fin(wf) && wf > 0 && fin(Vf) ? TAU / (wf * dtF) : NaN, kf = fin(Vf) && wf > 0 ? (wf * cref) / (2 * Vf) : NaN;
    if (fin(stepsPerCycle) && stepsPerCycle < 16) warnings.push(`Only ${stepsPerCycle.toFixed(0)} time steps per flutter cycle: reduce the time step for a firmer flutter speed.`);
    if (i.V / at.a > 0.6 || (fin(Vf) && Vf / at.a > 0.6)) warnings.push('Mach number above about 0.6: the vortex lattice is incompressible, so lift is under-predicted and transonic effects (the flutter dip) are absent.');
    warnings.push(`Default resolution is deliberately coarse for interactive use (${o.nc} × ${o.ns} panels, ${md.nm} modes, ${o.nw} wake rows): expect a few per cent discretisation error. For final values use 8 or more chordwise and 16 or more spanwise panels, 24 or more wake rows, 96 steps per period and a coupling tolerance of 1e-6, and run the convergence study.`);
    const mgn = margin(Vf * sg, i.V_d_eas), tms = gf.H.t, th = (a) => { const st = Math.max(1, Math.ceil(a.length / 300)); return a.filter((_, k) => k % st === 0); }, eta = N.range(o.ns, (j) => (j + 0.5) / o.ns);
    const spanLift = (r, sc) => N.range(o.ns, (j) => { let s = 0; for (let ic = 0; ic < o.nc; ic++) s += r.Fk[j * o.nc + ic]; return (s * sc) / md.lat.dy[j * o.nc] / 1e3; });
    return {
      kpis: [
        { key: 'tip_defl_uvlm_m', label: `Elastic tip deflection at ${i.n_load} g (coupled lattice–beam)`, value: fl.tip * aTrim, unit: 'm', note: `Strip-theory static analysis of this suite: ${fin(ref.tip_defl_m) ? ref.tip_defl_m.toPrecision(3) : '—'} m` },
        { key: 'tip_twist_uvlm_deg', label: 'Elastic tip twist', value: N.deg(fl.twist * aTrim), unit: 'deg', note: `Nose-up positive; strip theory: ${fin(ref.tip_twist_deg) ? ref.tip_twist_deg.toPrecision(3) : '—'} deg` },
        { key: 'lift_eff_uvlm', label: 'Flexible / rigid lift-curve slope', value: CLa_f / CLa_r, unit: '-', note: `Strip theory: ${fin(ref.lift_eff) ? ref.lift_eff.toPrecision(3) : '—'}` },
        { key: 'CLa_rigid_uvlm', label: 'Rigid wing lift-curve slope (vortex lattice)', value: CLa_r, unit: '1/rad', note: `Strip value supplied as input: ${i.cl.toPrecision(3)}` },
        { key: 'alpha_trim_deg', label: 'Incidence for the required lift (flexible)', value: N.deg(aTrim), unit: 'deg' },
        { key: 'root_bm_static_Nm', label: 'Root bending moment in steady flight', value: rmStatic, unit: 'N·m' },
        { key: 'gust_tip_defl_m', label: 'Peak incremental tip deflection in the gust', value: tipPk, unit: 'm' },
        { key: 'gust_root_bm_Nm', label: 'Peak incremental root bending moment in the gust', value: rmPk, unit: 'N·m' },
        { key: 'gust_dn_flex', label: 'Peak incremental load factor, flexible wing', value: dnF, unit: 'g', note: 'Wing lift increment over aircraft weight; root clamped (no heave or pitch relief)' },
        { key: 'gust_dn_rigid', label: 'Peak incremental load factor, rigid wing', value: dnR, unit: 'g' },
        { key: 'V_flutter_td_ms', label: isDiv ? 'Instability (divergence) speed from the time domain (TAS)' : 'Flutter speed from time-domain growth (TAS)', value: Vf, unit: 'm/s', status: !fin(Vf) || mgn > 0 ? 'ok' : 'bad', note: `p-k strip theory: ${fin(VfPk) ? VfPk.toFixed(1) + ' m/s' : 'no flutter found'}` },
        { key: 'f_flutter_td_Hz', label: 'Flutter frequency from the time domain', value: wf / TAU, unit: 'Hz', note: `p-k strip theory: ${fin(wfPk) ? (wfPk / TAU).toPrecision(3) + ' Hz' : '—'}` },
        { key: 'V_flutter_pk_ms', label: 'Flutter speed, p-k strip theory (reference)', value: VfPk, unit: 'm/s' },
        { key: 'flutter_ratio_td_pk', label: 'Time-domain / p-k flutter speed', value: Vf / VfPk, unit: '-', note: 'Within about 10% of 1 when both models follow the same mechanism (the lattice adds tip loss and spanwise induction, strip theory uses one lift-curve slope for every strip); a larger difference points to a different instability or to insufficient resolution' },
        { key: 'flutter_margin_td', label: 'Flutter margin from the time domain', value: fin(mgn) ? mgn : 9.99, unit: '-', status: mgn > 0 ? 'ok' : 'bad', note: MARGIN_NOTE + (fin(mgn) ? '' : '; no instability found, reported as the cap 9.99') },
        { key: 'k_flutter_td', label: 'Reduced frequency at flutter', value: kf, unit: '-' },
        { key: 'subiter_mean', label: 'Mean coupling sub-iterations per step (gust run)', value: gf.its, unit: '-' },
        { key: 'coupling_residual', label: 'Largest interface residual left after sub-iteration', value: gf.resMax, unit: '-', status: gf.resMax <= 100 * i.tol ? 'ok' : 'warn' },
        { key: 'virtual_work_err', label: 'Virtual-work mismatch of the load transfer', value: vwErr, unit: '-', status: vwErr < 1e-9 ? 'ok' : 'bad', note: 'Aerodynamic panel forces × interpolated displacements against nodal forces × nodal displacements' },
        { key: 'work_balance_err', label: 'Work balance across the interface in the gust run', value: Math.abs(gf.wA - gf.wS) / (Math.abs(gf.wA) + 1e-300), unit: '-', note: 'Work done by the air loads on the lattice against work received by the structure' },
        { key: 'force_transfer_err', label: 'Resultant-force mismatch of the load transfer', value: resErr, unit: '-', status: resErr < 1e-9 ? 'ok' : 'bad' },
        { key: 'dt_s', label: 'Time step at the flight speed', value: dt, unit: 's' },
      ].filter((k) => fin(k.value) || ['V_flutter_td_ms', 'tip_defl_uvlm_m'].includes(k.key)),
      plots: [
        { type: 'line', title: 'Gust response: tip deflection', xlabel: 'Time [s]', ylabel: 'Incremental tip deflection [m]', series: [{ name: 'Coupled lattice–beam', x: th(tms), y: th(gf.H.tip) }] },
        { type: 'line', title: 'Gust response: wing lift increment', xlabel: 'Time [s]', ylabel: 'Semi-wing lift increment [kN]', series: [{ name: 'Flexible', x: th(tms), y: th(gf.H.lift).map((v) => v / 1e3) }, { name: 'Rigid', x: th(gr.H.t), y: th(gr.H.lift).map((v) => v / 1e3), style: 'dash' }] },
        ...(Vs.length ? [{ type: 'line', title: 'Time-domain damping against speed', xlabel: 'True airspeed [m/s]', ylabel: 'Logarithmic decrement [-]', series: [{ name: 'Least-damped motion', x: Vs, y: dec.map((v) => (fin(v) ? N.clamp(v, -3, 3) : NaN)), style: 'line+points' }], annotations: [{ y: 0, label: 'Neutral' }, ...(fin(Vf) ? [{ x: Vf, label: 'Time domain' }] : []), ...(fin(VfPk) ? [{ x: VfPk, label: 'p-k' }] : [])] },
          { type: 'line', title: 'Time-domain growth rate against speed', xlabel: 'True airspeed [m/s]', ylabel: 'Growth rate σ [1/s]', series: [{ name: 'Energy-envelope fit', x: Vs, y: sig.map((v) => (fin(v) ? v : NaN)), style: 'line+points' }], annotations: [{ y: 0, label: 'Neutral' }] }] : []),
        { type: 'line', title: `Spanwise lift at ${i.n_load} g`, xlabel: 'Span fraction [-]', ylabel: 'Lift per unit span [kN/m]', series: [{ name: 'Flexible', x: eta, y: spanLift(fl, fin(aTrim) ? aTrim : 0) }, { name: 'Rigid (same total lift)', x: eta, y: spanLift(rg, W / rg.L), style: 'dash' }] },
      ],
      tables: [{ title: 'Beam modes retained (in vacuum)', columns: ['Mode', 'Frequency [Hz]', 'Torsional energy share'], rows: md.wn.map((w, k) => [k + 1, w / TAU, md.tors[k]]) }],
      outputs: { flutter_runs: nRun, flutter_steps_per_run: stepsPer },
      warnings,
      models: [`Unsteady vortex-ring lattice, ${o.nc} × ${o.ns} panels on the semi-span with a mirror image, Kutta condition by shedding the trailing-edge circulation into a ${o.nw}-row prescribed wake`, 'Unsteady Bernoulli panel loads: Kutta–Joukowski force on the bound segments plus ρ·∂Γ/∂t on the panel', `Euler–Bernoulli bending and St Venant torsion beam, ${o.ne} finite elements reduced to ${md.nm} natural modes, Newmark average-acceleration integration`, 'Partitioned strong coupling: sub-iterations with Aitken dynamic relaxation each time step', 'Conservative transfer: displacements interpolated with the beam shape functions (rigid sections normal to the elastic axis), loads by the transposed operator', 'Nonlinear time-domain flutter model: flutter onset from the exponential fit of the energy envelope of free oscillations at each speed', 'Gust-response model: 1-cosine gust convected over the lattice', 'Finite element structural dynamics model'],
      assumptions: ['Incompressible potential flow, small displacements: the lattice and the flat, free-stream-convected wake stay in the undeformed wing plane and the motion enters through the normal-wash boundary condition (no wake roll-up)', 'Wing clamped at the root; rigid-body heave and pitch of the aircraft are not included', 'No thickness, viscosity, stall or shock effects', 'The wake is truncated after the stated rows, the last one extending far downstream'],
    };
  },
  convergence: { param: 'nSpan', label: 'Spanwise panels (semi-span)', levels: [4, 6, 8, 12, 16], metric: 'tip_defl_uvlm_m' },
  calibration: { params: [{ key: 'GJ_root', min: 1, max: 1e11 }, { key: 'EI_root', min: 1, max: 1e12 }, { key: 'x_cg', min: 0.2, max: 0.7 }, { key: 'zeta', min: 0, max: 0.1 }], sweep: 'V', target: 'tip_defl_uvlm_m', note: 'Measured wing deflections in flight or in the wind tunnel at several speeds; flutter-test damping trends constrain the mass and stiffness distribution' },
  verify() {
    const out = [];
    // (a) beam finite elements against exact cantilever results
    const bi = { ...UVLM_BASE, x_cg: 0.4 }, fe = beamFE(bi, 8), nf = fe.nd - 3, f = new Array(nf).fill(0); f[nf - 3] = 100;
    const u = N.solve(fe.Kf, f), md0 = fsiModel(bi, { nc: 2, ns: 4, nw: 6, dtau: 0.5, ne: 8, nm: 6 }), mL = 60 / 6, Ia = mL * 0.25 ** 2, kb = md0.tors.findIndex((t) => t < 0.5), kt = md0.tors.findIndex((t) => t >= 0.5);
    out.push(N.check('Beam FE: tip deflection under a tip load = P·L³/(3·EI)', u[nf - 3], (100 * 216) / (3 * 2e5), 1e-9, 'Euler–Bernoulli cantilever (Hermite elements are exact for nodal loads)'));
    out.push(N.check('Beam FE: first bending frequency 1.8751²·√(EI/(m·L⁴))', md0.wn[kb], 1.875104 ** 2 * Math.sqrt(2e5 / (mL * 6 ** 4)), 1e-4, 'Euler–Bernoulli cantilever, exact'));
    out.push(N.check('Beam FE: first torsion frequency (π/2L)·√(GJ/Iα)', md0.wn[kt], (Math.PI / 12) * Math.sqrt(1e5 / Ia), 3e-3, 'Fixed–free torsion rod, exact (linear elements, 8 along the span)'));
    // (b) steady lattice against Prandtl lifting-line theory (Glauert series) for a rectangular wing of aspect ratio 16
    const AR = 16, wi = { ...UVLM_BASE, semi_span: 8, c_root: 1 }, lat = uvlmLattice(wi, 4, 16, 8, 0.25), md1 = { lat, nm: 0, Xc: [], Zq: [], wn: [], tipW: [], tipT: [], rootM: [] }, st = fsiStatic(md1, 1, 1, 1, true), CLa = st.L / (0.5 * 8);
    const nt = 24, th = N.range(nt, (k) => ((k + 0.5) * Math.PI) / (2 * nt)), mu = TAU / (4 * AR), Am = N.solve(th.map((t) => N.range(nt, (k) => { const n = 2 * k + 1; return Math.sin(n * t) * (1 + (mu * n) / Math.sin(t)); })), th.map(() => mu));
    out.push(N.check('Steady lattice lift-curve slope, rectangular wing AR 16, against lifting-line theory', CLa, Math.PI * AR * Am[0], 0.03, 'Prandtl lifting line (Glauert sine series, 24 terms); 4 × 16 panels, lifting-surface result lies 1–3% below'));
    // (c) Wagner's indicial lift: step change of incidence on a wing of aspect ratio 400 (two-dimensional limit)
    const w2 = { ...UVLM_BASE, semi_span: 200, c_root: 1 }, l2 = uvlmLattice(w2, 8, 4, 160, 0.125), m2 = { lat: l2, nm: 0, dxw: 0.125, Zc: [], Xc: [], Zq: [], Zm: [], wn: [], tipW: [], tipT: [], rootM: [], zeta: 0 }, wg = fsiMarch(m2, { rho: 1, U: 1, nSteps: 96, alpha0: 0.01, rigid: true, hist: true }), L2 = Math.PI * 0.01 * 200;
    const jones = (s) => 1 - JONES.A[0] * Math.exp(-JONES.b[0] * s) - JONES.A[1] * Math.exp(-JONES.b[1] * s);
    for (const s of [2, 6, 16]) { const k = Math.round(s / 2 / 0.125) - 1; out.push(N.check(`Wagner indicial lift at s = ${s} semi-chords`, wg.H.lift[k] / L2, jones(s), 0.03, 'Wagner (1925) function in R. T. Jones’ two-exponential form (itself within 1%); 8 chordwise panels, U·Δt = c/8')); }
    // (d) conservative transfer and coupling on a flexible wing
    const r = N.kv(uvlmfsi.run({ ...UVLM_BASE, nSpeeds: 0 }));
    out.push(N.check('Virtual work is conserved by the load / displacement transfer', r.virtual_work_err, 0, 1e-10, 'Loads transferred with the transpose of the displacement interpolation'));
    out.push(N.check('Resultant force is conserved by the load transfer', r.force_transfer_err, 0, 1e-10, 'Partition of unity of the beam shape functions'));
    out.push(N.check('Work done by the air loads equals work received by the structure (gust run)', r.work_balance_err, 0, 1e-9, 'Discrete energy consistency of the coupled step'));
    // (e) time-domain flutter of a slender wing against p-k strip theory of the same wing (strip theory is exact as AR → ∞)
    const g = { semi_span: 18.288, c_root: 1.8288, taper: 1, sweep_deg: 0, x_ea: 0.33, x_cg: 0.43, r_alpha: Math.sqrt(8.64 / 35.71) / 1.8288, EI_root: 9.77e6 * 81, GJ_root: 0.987e6 * 9, stiff_exp: 0, m_semi: 35.71 * 18.288, zeta: 0, cl: TAU, alt_m: 0, V_d_eas: 100 };
    const pkf = wingFlutterSolve(g, 1.225, 3, 30, 400), mdg = fsiModel(g, { nc: 4, ns: 10, nw: 20, dtau: 0.25, ne: 8, nm: 6 }), gr = (V) => fsiGrowth(mdg, 1.225, V, 5, { maxIt: 10, tol: 1e-8, maxSteps: 400 }).sigma;
    let Vtd = NaN; try { Vtd = N.brent(gr, 0.85 * pkf.flutter.P, 1.3 * pkf.flutter.P, 1, 8); } catch { Vtd = NaN; }
    out.push(N.check('Time-domain flutter speed of a slender uniform wing (AR 20) against p-k strip theory', Vtd, pkf.flutter.P, 0.1, 'Goland-section wing lengthened threefold with frequencies preserved; the lattice adds tip loss, which raises the flutter speed by a few per cent'));
    // (f) the same comparison on a tapered wing whose second-bending and torsion branches approach each other below the flutter
    // speed: the time march has no root tracking, so agreement shows the p-k sweep kept the fluttering branch
    const tw = { semi_span: 13.54, c_root: 2.837, taper: 0.59, sweep_deg: 3, x_ea: 0.4, x_cg: 0.45, r_alpha: 0.25, EI_root: 2.669e7, GJ_root: 1.897e7, stiff_exp: 3, m_semi: 2880, zeta: 0.02, cl: 5.316, alt_m: 6100, V_d_eas: 160 }, rt = 0.6597;
    const pkt = wingFlutterSolve(tw, rt, 3, 30, 480), mdt = fsiModel(tw, { nc: 4, ns: 10, nw: 20, dtau: 0.25, ne: 8, nm: 6 }), grt = (V) => fsiGrowth(mdt, rt, V, 5, { maxIt: 10, tol: 1e-8, maxSteps: 400 }).sigma;
    let Vtt = NaN; try { Vtt = N.brent(grt, 0.75 * pkt.flutter.P, 1.3 * pkt.flutter.P, 1, 8); } catch { Vtt = NaN; }
    out.push(N.check('Time-domain flutter speed of a tapered wing (AR 11) with approaching branches against p-k strip theory', Vtt, pkt.flutter.P, 0.1, 'Strip theory with the finite-span lift-curve slope against the lattice: the two differ in spanwise load distribution, hence the 10% tolerance'));
    return out;
  },
  recommend(res, i) {
    const o = res.outputs, out = [];
    if (fin(o.V_flutter_td_ms) && o.flutter_margin_td < 0) out.push({ severity: 'critical', title: 'Time-domain instability inside 1.15 VD', detail: `The coupled lattice–beam model becomes unstable at ${o.V_flutter_td_ms.toFixed(0)} m/s TAS${fin(o.f_flutter_td_Hz) ? ` (${o.f_flutter_td_Hz.toFixed(1)} Hz)` : ''}.`, action: 'Raise torsional stiffness, move section mass forward of the elastic axis, or add tip mass balance; then repeat with finer panels and time step.', basis: 'Flutter clearance to 1.15 VD' });
    if (fin(o.flutter_ratio_td_pk) && Math.abs(o.flutter_ratio_td_pk - 1) > 0.25) out.push({ severity: 'advise', title: 'Time-domain and p-k flutter speeds differ by more than 25%', detail: `Ratio ${o.flutter_ratio_td_pk.toFixed(2)}.`, action: 'Check panel, mode and time-step convergence first (the two agree within 10% on the verification wings). The p-k sweep continues every mode and is cross-checked by a tracking-free determinant count, and the time-domain run has no mode tracking at all, so a remaining gap is a model-form difference: for low aspect ratio or strongly swept wings the lattice result is the more credible of the two. Use the lower speed for the margin.', basis: 'Model-form difference between strip theory and lifting-surface aerodynamics' });
    if (o.gust_dn_flex > 0 && o.gust_dn_rigid > 0) out.push({ severity: 'info', title: 'Gust loads with flexibility', detail: `Peak incremental load factor ${o.gust_dn_flex.toFixed(2)} g flexible against ${o.gust_dn_rigid.toFixed(2)} g rigid; root bending increment ${(o.gust_root_bm_Nm / 1e3).toFixed(0)} kN·m on ${(o.root_bm_static_Nm / 1e3).toFixed(0)} kN·m in steady flight.`, action: 'Use the flexible root bending increment for the gust load case in Suite 2; dynamic overshoot or relief here translates directly into wing structural mass and fuel burn.', basis: 'Coupled time-domain gust response' });
    if (o.coupling_residual > 100 * i.tol) out.push({ severity: 'warn', title: 'Coupling not converged within the sub-iteration limit', detail: `Residual ${o.coupling_residual.toExponential(1)}.`, action: 'Raise the maximum number of sub-iterations or reduce the time step.', basis: 'Strong-coupling convergence criterion' });
    return out;
  },
};

// ---- 9. Navier–Stokes / structure coupling (immersed boundary) -----------------------------------
const IB_BODIES = ['Cylinder (vortex-induced vibration)', 'Flat plate (pitch–plunge)'];
/** Amplitude, mean and zero-crossing frequency of a uniformly sampled signal over [k0, k1). */
function oscStats(y, dt, k0, k1) {
  let mn = Infinity, mx = -Infinity, m = 0; for (let k = k0; k < k1; k++) { const v = y[k]; if (v < mn) mn = v; if (v > mx) mx = v; m += v; } m /= Math.max(1, k1 - k0);
  // frequency from all sign changes about the mean (half-periods), so that short records still give an estimate
  let first = NaN, last = NaN, nz = 0, ss = 0; for (let k = k0; k < k1; k++) { ss += (y[k] - m) ** 2; if (k > k0 && (y[k - 1] - m) * (y[k] - m) < 0) { const t = (k - 1 + (m - y[k - 1]) / (y[k] - y[k - 1])) * dt; if (nz === 0) first = t; last = t; nz++; } }
  return { amp: 0.5 * (mx - mn), mean: m, rms: Math.sqrt(ss / Math.max(1, k1 - k0)), f: nz > 1 ? (nz - 1) / (2 * (last - first)) : 0, max: mx, min: mn };
}
/**
 * Coupled immersed-boundary Navier–Stokes / rigid-body simulation in non-dimensional units (U = ρ = 1; length unit = cylinder
 * diameter or plate semi-chord). o.fixed holds the body (after a short transverse kick that triggers shedding); o.still
 * removes the stream (body released in fluid at rest, time still measured in the same units). o.load(vNew, vOld, dt, t) replaces
 * the fluid sub-solve by a prescribed load [Fx, Fy, Mz] on the body, the mean over the step that ends at time t with the body
 * velocities vNew (verification of the coupling algebra).
 */
function ibRun(i, ctx, o = {}) {
  const plate = i.body === IB_BODIES[1], cells = i.cells, h = 1 / cells, nx = Math.max(16, Math.round(i.len * cells)), ny = 2 * Math.max(4, Math.round((i.height * cells) / 2)), Lre = plate ? 2 : 1, nu = Lre / i.Re;
  const s = nsCreate({ nx, ny, h, nu, U: o.still ? 0 : 1, rho: 1, upwind: i.upwind }), y0 = 0.5 * ny * h, x0 = o.still ? 0.5 * nx * h : Math.min(i.x_body, 0.5 * nx * h), thick = plate ? Math.max(2 * i.tc, 2 * h) : 0;
  const body = plate ? nsBody({ type: 'plate', c: 2, thick, xm: -i.a_ea, x: x0, y: y0 }) : nsBody({ type: 'circle', R: 0.5, x: x0, y: y0 });
  const dt = Math.min((i.cfl * h) / 1, (0.2 * h * h) / nu), nSw = Math.max(1, Math.round(i.n_sweep)), Urs = nSw === 1 ? [i.Ur_min] : N.linspace(i.Ur_min, i.Ur_max, nSw);
  const stages = Urs.map((Ur, k) => ({ Ur, n: Math.max(20, Math.round((i.t_stage + (k === 0 ? i.t_start : 0)) / dt)) })), nTot = N.sum(stages.map((g) => g.n));
  // structure: q = [y] or [y, x] for the cylinder; q = [h (down), α (nose-up)] for the plate.
  // Time integration: the fluid sub-step returns the momentum it exchanged with the body over the step, i.e. the step-mean load
  // H. It is paired with the structural equation in the same time-centred (impulse) form, M·ā + C·(v + ½Δt·ā) + K·(q + ½Δt·v + ¼Δt²·ā) = H,
  // with the step-mean acceleration ā as the unknown and v ← v + Δt·ā, q ← q + Δt·v + ½Δt²·ā (the average-acceleration rule).
  // Body momentum then changes by exactly H·Δt, the work of the fluid is H·Δq, and an added-mass reaction H = −m_a·ā is carried
  // in full by the interface operator. In the arrays below a, a1, aj, an, ap are step-mean accelerations.
  const fixed = !!o.fixed, nd = fixed ? 0 : plate ? 2 : i.inline ? 2 : 1, M = [[0, 0], [0, 0]], Kd = [0, 0], Cd = [0, 0], Mv = [0, 0];
  const setStage = (Ur) => {
    if (plate) { const m = i.mu * Math.PI, wa = 1 / Ur, wh = i.w_ratio * wa, Ia = m * i.r_alpha ** 2; M[0][0] = m; M[0][1] = M[1][0] = m * i.x_alpha; M[1][1] = Ia; Kd[0] = m * wh * wh; Kd[1] = Ia * wa * wa; Cd[0] = 2 * i.zeta * m * wh; Cd[1] = 2 * i.zeta * Ia * wa; Mv[0] = Math.PI; Mv[1] = Math.PI * (0.125 + i.a_ea ** 2); }
    else { const m = (i.m_star * Math.PI) / 4, wn = TAU / Ur; M[0][0] = M[1][1] = m; Kd[0] = Kd[1] = m * wn * wn; Cd[0] = Cd[1] = 2 * i.zeta * m * wn; Mv[0] = Mv[1] = Math.PI / 4; }
  };
  const q = [0, 0], v = [0, 0], a = [0, 0], a1 = [0, 0], aj = [0, 0], an = [0, 0], ap = [0, 0], vk = [0, 0], qk = [0, 0], cc = [0, 0], dH = [null, null];
  if (!fixed) { if (plate) q[1] = N.rad(i.alpha0_deg); else v[0] = 0.2; }
  const pose = (qq) => { if (plate) { body.y = y0 - qq[0]; body.ang = -qq[1]; } else { body.y = y0 + qq[0]; body.x = x0 + (nd > 1 ? qq[1] : 0); } };
  const vel = (vv) => { if (plate) { body.vy = -vv[0]; body.om = -vv[1]; } else { body.vy = vv[0]; body.vx = nd > 1 ? vv[1] : 0; } };
  const nNewton = fixed ? 0 : Math.max(0, Math.round(i.subiter)), T = [], Y = [], X = [], CL = [], CDh = [], CM = [], marks = [], slots = nd > 1 ? [1, 2] : [1], del = 1;
  // the fluid inside the mask starts with the free-stream momentum, not the body's: that is the reference for the interior-inertia term
  pose(q); body.vx = s.U; body.vy = 0; body.om = 0; const P0 = nsInterior(s, body); body.vx = 0;
  let Pp = [0, 0, 0, P0.Px, P0.Py, P0.Lz], its = 0, resMax = 0, resSum = 0, cflMax = 0, divMax = 0, step = 0, amMax = 0, tNow = 0;
  /** One fluid sub-solve from the stored predictor field for step-mean body accelerations acc: force – project – force. Returns [Hx, Hy, Hz, Px, Py, Lz]. */
  const evalAt = (acc) => {
    nsRecall(s, 0); for (let d = 0; d < nd; d++) vk[d] = v[d] + dt * acc[d]; vel(vk);
    if (o.load) return [...o.load(vk, v, dt, tNow), 0, 0, 0];
    const f0 = nsForce(s, body, dt); nsProject(s, dt); const f1 = nsForce(s, body, dt), pp = Pp;
    return [f0.Fx + f1.Fx + (f1.Px - pp[3]) / dt, f0.Fy + f1.Fy + (f1.Py - pp[4]) / dt, f0.Mz + f1.Mz + (f1.Lz - pp[5]) / dt + body.vx * f1.Py - body.vy * f1.Px, f1.Px, f1.Py, f1.Lz];
  };
  const gen = (H, d) => (plate ? (d === 0 ? -H[1] : -H[2]) : d === 0 ? H[1] : H[0]);
  for (let sg = 0; sg < stages.length; sg++) {
    setStage(stages[sg].Ur); marks.push(step);
    for (let k = 0; k < stages[sg].n; k++, step++) {
      const t = tNow = (step + 1) * dt; let H = null, res = 0;
      nsPredict(s, dt);
      if (fixed || nNewton === 0) {
        // fixed body, or loosely coupled: one forcing pass with the load of the previous pressure field and a lagged added-mass term
        if (fixed) { const vy = t < 4 ? 0.3 * Math.sin((Math.PI * t) / 2) : 0; body.vy = vy; body.y = t < 4 ? body.y + vy * dt : y0; }
        else { for (let d = 0; d < nd; d++) { qk[d] = q[d] + dt * v[d] + 0.5 * dt * dt * a[d]; vk[d] = v[d] + dt * a[d]; } pose(qk); vel(vk); } // motion predicted with the mean acceleration of the previous step
        const f = nsForce(s, body, dt), pp = Pp; nsProject(s, dt);
        H = o.load ? [...o.load(vk, v, dt, t), 0, 0, 0] : [f.Fx + (f.Px - pp[3]) / dt, f.Fy + (f.Py - pp[4]) / dt, f.Mz + (f.Lz - pp[5]) / dt + body.vx * f.Py - body.vy * f.Px, f.Px, f.Py, f.Lz];
        if (!fixed) {
          const b0 = gen(H, 0) + Mv[0] * a[0] - Cd[0] * v[0] - Kd[0] * (q[0] + 0.5 * dt * v[0]), A00 = M[0][0] + Mv[0] + 0.5 * dt * Cd[0] + 0.25 * dt * dt * Kd[0];
          if (nd === 1) aj[0] = b0 / A00;
          else { const b1 = gen(H, 1) + Mv[1] * a[1] - Cd[1] * v[1] - Kd[1] * (q[1] + 0.5 * dt * v[1]), A11 = M[1][1] + Mv[1] + 0.5 * dt * Cd[1] + 0.25 * dt * dt * Kd[1], A01 = M[0][1], det = A00 * A11 - A01 * A01; aj[0] = (b0 * A11 - A01 * b1) / det; aj[1] = (A00 * b1 - A01 * b0) / det; }
          res = dt * Math.hypot(aj[0] - a[0], nd > 1 ? aj[1] - a[1] : 0); its++;
        }
      } else {
        // strong coupling: for a frozen pose the fluid sub-step is affine in the body's step-mean accelerations, so nd + 1 sub-solves give the
        // load, its discrete added-mass operator J and the matching sensitivity fields; the interface condition is then solved
        // directly (a Newton step) and the field corrected by superposition. With one iteration per step J and the sensitivity
        // fields are refreshed every fourth step (quasi-Newton); with more they are rebuilt in every iteration.
        nsStore(s, 0); const ext = k >= 2 ? 1 : 0; aj[0] = a[0] + ext * (a[0] - a1[0]); aj[1] = a[1] + ext * (a[1] - a1[1]);
        for (let j = 0; j < nNewton; j++) {
          for (let d = 0; d < nd; d++) qk[d] = q[d] + dt * v[d] + 0.5 * dt * dt * aj[d];
          pose(qk); const fresh = nNewton > 1 || k % 4 === 0;
          if (fresh) for (let d = 0; d < nd; d++) { ap[0] = aj[0]; ap[1] = aj[1]; ap[d] += del; dH[d] = evalAt(ap); nsStore(s, d + 1); }
          const Hb = evalAt(aj);
          if (fresh) for (let d = 0; d < nd; d++) { nsDiff(s, d + 1); for (let m = 0; m < 6; m++) dH[d][m] -= Hb[m]; }
          const J00 = gen(dH[0], 0) / del, b0 = gen(Hb, 0) - Cd[0] * v[0] - Kd[0] * (q[0] + 0.5 * dt * v[0]), A00 = M[0][0] + 0.5 * dt * Cd[0] + 0.25 * dt * dt * Kd[0] - J00;
          if (nd === 1) { an[0] = (b0 - J00 * aj[0]) / A00; an[1] = 0; amMax = Math.max(amMax, -J00 / M[0][0]); }
          else {
            const J10 = gen(dH[0], 1) / del, J01 = gen(dH[1], 0) / del, J11 = gen(dH[1], 1) / del;
            const b1 = gen(Hb, 1) - Cd[1] * v[1] - Kd[1] * (q[1] + 0.5 * dt * v[1]), A11 = M[1][1] + 0.5 * dt * Cd[1] + 0.25 * dt * dt * Kd[1] - J11, A01 = M[0][1] - J01, A10 = M[1][0] - J10;
            const r0 = b0 - J00 * aj[0] - J01 * aj[1], r1 = b1 - J10 * aj[0] - J11 * aj[1], det = A00 * A11 - A01 * A10;
            an[0] = (r0 * A11 - A01 * r1) / det; an[1] = (A00 * r1 - A10 * r0) / det; amMax = Math.max(amMax, -J00 / M[0][0], -J11 / M[1][1]);
          }
          for (let d = 0; d < nd; d++) cc[d] = (an[d] - aj[d]) / del;
          nsCombine(s, 1, slots, cc);
          H = Hb.map((x, m) => x + (nd > 1 ? cc[0] * dH[0][m] + cc[1] * dH[1][m] : cc[0] * dH[0][m]));
          // pose inconsistency left by this iteration: structural position minus the position the fluid was solved with
          res = 0.5 * dt * dt * Math.hypot(an[0] - aj[0], an[1] - aj[1]); aj[0] = an[0]; aj[1] = an[1]; its++;
          if (res < i.tol) break;
        }
        nsProject(s, dt);
      }
      Pp = H;
      if (!fixed) { for (let d = 0; d < nd; d++) { q[d] += dt * v[d] + 0.5 * dt * dt * aj[d]; v[d] += dt * aj[d]; a1[d] = a[d]; a[d] = aj[d]; } if (res > resMax) resMax = res; resSum += res; }
      const Hy = H[1], Hx = H[0], Hz = H[2];
      if (!fin(Hy) || !fin(q[0]) || Math.abs(q[0]) > 0.45 * ny * h) throw new Error('The coupled flow solution diverged or the body left the domain: reduce the Courant number, raise the sub-iteration count or enlarge the domain.');
      T.push(t); Y.push(plate ? -q[0] : fixed ? body.y - y0 : q[0]); X.push(plate ? q[1] : nd > 1 ? q[1] : 0); CL.push(Hy / (0.5 * Lre)); CDh.push(Hx / (0.5 * Lre)); CM.push(-Hz / (0.5 * Lre * Lre));
      if (step % 40 === 0) { cflMax = Math.max(cflMax, (nsMaxSpeed(s) * dt) / h); divMax = Math.max(divMax, nsDivergence(s)); ctx?.progress?.(step / nTot, 'Coupled flow–structure time stepping'); }
    }
  }
  const snap = nsVorticity(s, Math.max(1, Math.ceil(nx / 160)));
  // lift and drag for the statistics are box-filtered over about one convective time to remove the grid-crossing noise of the moving interface
  const nb = Math.max(1, Math.round(0.6 / dt)), box = (y) => { const c = [0]; for (let k = 0; k < y.length; k++) c.push(c[k] + y[k]); return y.map((_, k) => { const a = Math.max(0, k - nb), b = Math.min(y.length, k + nb + 1); return (c[b] - c[a]) / (b - a); }); }, CLf = box(CL), CDf = box(CDh);
  const st = stages.map((g, k) => { const k0 = marks[k] + Math.floor((k === 0 ? 0.5 : 0.25) * g.n), k1 = marks[k] + g.n; return { Ur: g.Ur, y: oscStats(Y, dt, k0, k1), x: oscStats(X, dt, k0, k1), cl: oscStats(CLf, dt, k0, k1), cd: oscStats(CDf, dt, k0, k1) }; });
  return { plate, nx, ny, h, dt, nTot, T, Y, X, CL, CD: CDh, CM, st, marks, snap, its: its / nTot, resMax, resMean: resSum / nTot, cflMax, divMax, body: { x: body.x, y: body.y, ang: body.ang, thick }, y0, x0, nNewton, nd, Lre, amMax };
}
const IB_BASE = { body: IB_BODIES[0], Re: 150, L_ref: 0.05, U_ms: 10, rho: 1.225, m_star: 4, zeta: 0.005, inline: false, mu: 4, a_ea: 0, x_alpha: 0.1, r_alpha: 0.5, w_ratio: 0.5, alpha0_deg: 4, tc: 0.12, Ur_min: 3, Ur_max: 7, n_sweep: 3, cells: 5, len: 12.8, height: 6.4, x_body: 3.2, cfl: 0.3, t_start: 10, t_stage: 12, subiter: 1, tol: 1e-5, upwind: 1.5 };
const cfdfsi = {
  id: 'cfdfsi', title: 'Navier–Stokes flow coupled to a moving structure (vortex-induced vibration, plate flutter)', fidelity: 'numerical',
  summary: 'A two-dimensional incompressible Navier–Stokes solver with an immersed boundary is strongly coupled to a spring-mounted body: a circular cylinder free to vibrate in the vortex street it sheds (lock-in), or a flat plate free to pitch and plunge. The reduced velocity is stepped through a range to find where the motion grows; vorticity field, motion and force histories and spectra are reported.',
  equations: ['Partitioned CFD–CSD coupling equations', 'Interface displacement compatibility equations', 'Interface traction equilibrium equations', 'Strongly coupled and loosely coupled time-marching formulations', 'Coupled aerodynamic–structural equations of motion'],
  inputs: [
    { key: 'body', label: 'Body', type: 'select', options: IB_BODIES, default: IB_BODIES[0], group: 'Body' },
    { key: 'Re', label: 'Reynolds number (diameter or chord)', unit: '-', default: 150, min: 20, max: 2000, group: 'Flow', help: 'Laminar two-dimensional shedding: 50–200 for a cylinder. The real flight Reynolds number is far higher and cannot be resolved here' },
    { key: 'L_ref', label: 'Diameter or chord', unit: 'm', default: 0.05, min: 1e-4, group: 'Scales', help: 'Only used to express frequencies, stiffness and forces in physical units' },
    { key: 'U_ms', label: 'Flow speed', unit: 'm/s', default: 10, min: 0.01, group: 'Scales' },
    { key: 'rho', label: 'Fluid density', unit: 'kg/m³', default: 1.225, min: 1e-3, group: 'Scales' },
    { key: 'm_star', label: 'Cylinder mass ratio m / (ρ·π·D²/4)', unit: '-', default: 4, min: 0.3, max: 1000, group: 'Cylinder', help: 'About 1–10 in water, hundreds to thousands for metal tubes in air. Low values need more coupling sub-iterations' },
    { key: 'inline', label: 'Cylinder also free in the stream direction', type: 'bool', default: false, group: 'Cylinder' },
    { key: 'mu', label: 'Plate mass ratio m / (π·ρ·b²)', unit: '-', default: 4, min: 0.5, max: 1000, group: 'Plate', help: 'b is the semi-chord' },
    { key: 'a_ea', label: 'Plate elastic axis aft of mid-chord / semi-chord', unit: '-', default: 0, min: -0.8, max: 0.8, group: 'Plate' },
    { key: 'x_alpha', label: 'Plate mass centre aft of the elastic axis / semi-chord', unit: '-', default: 0.1, min: -0.5, max: 0.8, group: 'Plate' },
    { key: 'r_alpha', label: 'Plate radius of gyration / semi-chord', unit: '-', default: 0.5, min: 0.1, max: 1.5, group: 'Plate' },
    { key: 'w_ratio', label: 'Plunge / pitch natural frequency ratio', unit: '-', default: 0.5, min: 0.05, max: 3, group: 'Plate' },
    { key: 'alpha0_deg', label: 'Plate initial pitch disturbance', unit: 'deg', default: 4, min: 0.1, max: 30, group: 'Plate' },
    { key: 'tc', label: 'Plate thickness / chord', unit: '-', default: 0.12, min: 0.01, max: 0.4, group: 'Plate', help: 'Rounded-edge plate; at least two grid cells are always used' },
    { key: 'zeta', label: 'Structural damping ratio', unit: '-', default: 0.005, min: 0, max: 0.5, group: 'Structure' },
    { key: 'Ur_min', label: 'Reduced velocity, first stage', unit: '-', default: 3, min: 0.3, max: 40, group: 'Structure', help: 'Cylinder: U/(f_n·D), lock-in is expected around 4–8. Plate: U/(b·ω_α), of order 1 for a light plate' },
    { key: 'Ur_max', label: 'Reduced velocity, last stage', unit: '-', default: 7, min: 0.3, max: 40, group: 'Structure' },
    { key: 'n_sweep', label: 'Reduced-velocity stages', unit: '', default: 3, min: 1, max: 24, step: 1, discrete: true, group: 'Structure', help: 'The spring stiffness is stepped during one continuous simulation, as in a wind-tunnel speed sweep' },
    { key: 'cells', label: 'Grid cells per diameter / per semi-chord', unit: '', default: 5, min: 3, max: 40, step: 1, discrete: true, group: 'Numerics', help: 'The default of 5 is a demonstration grid; 16–32 are needed for quantitative forces and amplitudes' },
    { key: 'len', label: 'Domain length / reference length', unit: '-', default: 9.6, min: 6, max: 60, group: 'Numerics', help: 'The default is short to keep the run interactive; 12.8 or more keeps the outflow boundary out of the near wake' },
    { key: 'height', label: 'Domain height / reference length', unit: '-', default: 6.4, min: 3, max: 40, group: 'Numerics', help: 'Free-slip walls; a cell count that is a power of two uses the fast transform' },
    { key: 'x_body', label: 'Body position from the inlet / reference length', unit: '-', default: 3.2, min: 1.5, max: 20, group: 'Numerics' },
    { key: 'cfl', label: 'Courant number U·Δt/h', unit: '-', default: 0.3, min: 0.05, max: 0.45, group: 'Numerics' },
    { key: 't_start', label: 'Start-up time before the first stage is evaluated', unit: 'L/U', default: 10, min: 2, max: 400, group: 'Numerics' },
    { key: 't_stage', label: 'Time per stage', unit: 'L/U', default: 12, min: 4, max: 2000, group: 'Numerics', help: 'Amplitudes need tens of cycles to settle; the default shows the trend only' },
    { key: 'subiter', label: 'Strong-coupling iterations per step', unit: '', default: 1, min: 0, max: 6, step: 1, discrete: true, group: 'Numerics', help: '0 = loosely coupled (staggered, lagged added mass). 1 = implicit interface solve with the discrete added-mass operator of the fluid sub-step, refreshed every fourth step. 2 or more = operator rebuilt every iteration and the body position iterated as well' },
    { key: 'tol', label: 'Coupling tolerance on the body position', unit: 'L', default: 1e-5, min: 1e-10, max: 0.01, group: 'Numerics', help: 'Extra coupling iterations stop when the structural position and the position used by the fluid differ by less than this many reference lengths' },
    { key: 'upwind', label: 'Upwind dissipation factor', unit: '-', default: 1.5, min: 0.5, max: 3, group: 'Numerics', help: '1 = third-order upwind, 3 = Kawamura–Kuwahara' },
  ],
  defaults(c, up, d) {
    const wing = c.wing.S_m2 > 0, at = isa(c.atm.alt_m, c.atm.dISA_K);
    // aeroplanes: a wing-section-like plate at the cruise speed; rotorcraft and multirotors: a landing-gear or boom tube
    return wing ? { body: IB_BODIES[1], L_ref: d.mac || undefined, U_ms: c.flight.V_ms, rho: at.rho, Re: 300, Ur_min: 0.7, Ur_max: 1.7, n_sweep: 2, t_start: 4, t_stage: 12, zeta: c.struct.zeta }
      : { body: IB_BODIES[0], L_ref: N.clamp(0.04 * (c.fuselage.dia_m || 1), 0.005, 0.2), U_ms: Math.max(1, c.flight.V_ms), rho: at.rho };
  },
  run(i, ctx) {
    const R = ibRun(i, ctx), plate = R.plate, warnings = [], Ls = plate ? i.L_ref / 2 : i.L_ref, tS = Ls / i.U_ms, st = R.st, fS = i.rho * i.U_ms ** 2 * Ls;
    const amps = st.map((g) => (plate ? g.x.amp : g.y.amp)), kp = N.argmax(amps), pkS = st[kp], thin = (arr, k0 = 0) => { const n = arr.length - k0, sk = Math.max(1, Math.ceil(n / 380)), o = []; for (let k = k0; k < arr.length; k += sk) o.push(arr[k]); return o; };
    const k0 = R.marks[kp], k1 = kp + 1 < R.marks.length ? R.marks[kp + 1] : R.nTot, seg = (arr) => arr.slice(k0 + Math.floor(0.3 * (k1 - k0)), k1), spY = N.spectrum(seg(plate ? R.X : R.Y), R.dt), spL = N.spectrum(seg(R.CL), R.dt), fcut = (sp) => { const n = Math.min(sp.f.length, Math.max(8, sp.f.findIndex((f) => f > (plate ? 0.5 : 0.6)))); return { f: sp.f.slice(0, n), a: sp.amp.slice(0, n) }; }, sy = fcut(spY), sl = fcut(spL);
    const nuAir = 1.5e-5, ReReal = (i.U_ms * i.L_ref) / nuAir, cycles = (plate ? pkS.x.f : pkS.y.f) * i.t_stage;
    warnings.push(`Demonstration resolution: ${i.cells} cells per ${plate ? 'semi-chord' : 'diameter'} on a ${R.nx} × ${R.ny} grid, domain height ${(R.ny * R.h).toFixed(1)} reference lengths with free-slip walls (blockage ${(100 * (plate ? 2 * Math.sin(Math.max(pkS.x.amp, 0.05)) : 1) / (R.ny * R.h)).toFixed(0)}%). Forces and amplitudes are indicative; use 16 or more cells per reference length, a taller and longer domain (length 12.8 or more) and longer stages for quantitative work.`);
    warnings.push(`The simulated Reynolds number is ${i.Re}; at the stated size and speed the real value would be about ${ReReal.toExponential(1)} in sea-level air. Laminar two-dimensional results do not carry over to turbulent, three-dimensional flow.`);
    if (i.cells < 12) warnings.push('At this resolution the smeared interface adds numerical damping of a few per cent of critical and over-states the added mass (free-decay test in still fluid), so response amplitudes are under-predicted and the response frequency is somewhat low.');
    if (cycles < 8) warnings.push(`Each stage covers only about ${cycles.toFixed(1)} oscillation cycles: amplitudes have not fully settled. Lengthen the time per stage.`);
    if (R.nNewton === 0) warnings.push('Loosely coupled run: the fluid load lags the body motion by one step. This is adequate for heavy bodies only; use at least one strong-coupling iteration when the mass ratio is below about 5.');
    if (R.nNewton > 0 && R.resMax > 5e-3) warnings.push(`The body position used by the fluid differs from the structural position by up to ${R.resMax.toExponential(1)} reference lengths within a step; add a strong-coupling iteration or reduce the time step.`);
    if (R.cflMax > 0.8) warnings.push(`Local Courant number reached ${R.cflMax.toFixed(2)}; reduce the Courant number input.`);
    const common = [
      { key: 'subiter_mean', label: 'Mean coupling iterations per step', value: R.its, unit: '-', note: R.nNewton > 1 ? `${R.nd + 1} fluid sub-solves per iteration` : R.nNewton ? `1 fluid sub-solve per step plus ${R.nd} more every fourth step` : 'staggered' },
      { key: 'coupling_residual', label: R.nNewton ? 'Largest body-position mismatch between fluid and structure in a step' : 'Largest body-velocity change not yet seen by the fluid in a step', value: R.resMax, unit: R.nNewton ? 'L' : 'U', status: R.resMax <= 5e-3 ? 'ok' : 'warn', note: R.nNewton ? 'Velocity and load are matched exactly for the position used; this is what remains' : 'Staggered scheme' },
      ...(R.nNewton ? [{ key: 'added_mass_ratio', label: 'Largest implicit added mass / structural mass', value: R.amMax, unit: '-', note: 'Added mass (or inertia) seen by the coupling Jacobian of the force–project–force sub-step and carried implicitly in full; the part of the fluid reaction that one sub-step does not yet develop follows one step later' }] : []),
      { key: 'cfl_max', label: 'Largest local Courant number', value: R.cflMax, unit: '-', status: R.cflMax <= 0.8 ? 'ok' : 'warn' },
      { key: 'div_max', label: 'Largest velocity divergence after projection', value: R.divMax, unit: 'U/L', status: R.divMax < 1e-8 ? 'ok' : 'bad' },
      { key: 'n_cells', label: 'Grid cells', value: R.nx * R.ny, unit: '-' }, { key: 'n_steps', label: 'Time steps', value: R.nTot, unit: '-' },
    ];
    const vort = { type: 'heat', title: 'Vorticity field at the end of the run', xlabel: `x / ${plate ? 'semi-chord' : 'diameter'} [-]`, ylabel: `y / ${plate ? 'semi-chord' : 'diameter'} [-]`, zlabel: 'Vorticity ω·L/U [-]', x: R.snap.x, y: R.snap.y, z: R.snap.z.map((r) => r.map((v) => N.clamp(v, -4, 4))), diverging: true, equalAspect: true, contours: 12,
      overlay: [plate ? (() => { const ca = Math.cos(R.body.ang), sa = Math.sin(R.body.ang), xm = -i.a_ea; return { name: 'Plate', x: [R.body.x + (xm - 1) * ca, R.body.x + (xm + 1) * ca], y: [R.body.y + (xm - 1) * sa, R.body.y + (xm + 1) * sa] }; })() : { name: 'Cylinder', x: N.range(25, (k) => R.body.x + 0.5 * Math.cos((TAU * k) / 24)), y: N.range(25, (k) => R.body.y + 0.5 * Math.sin((TAU * k) / 24)) }] };
    const tt = thin(R.T), ann = R.marks.slice(1).map((m, k) => ({ x: R.T[m], label: `Ur ${st[k + 1].Ur.toPrecision(3)}` }));
    if (!plate) {
      const fn = 1 / pkS.Ur, lock = pkS.y.amp > 0.1 && Math.abs(pkS.cl.f / Math.max(pkS.y.f, 1e-9) - 1) < 0.1, St0 = st[0].cl.f, fnHz = fn / tS, mPhys = (i.m_star * Math.PI * i.rho * i.L_ref ** 2) / 4;
      if (amps[kp] > 0.1 && (kp === 0 || kp === st.length - 1) && st.length > 1) warnings.push('The largest amplitude is at the end of the reduced-velocity range: widen the range to bracket the lock-in peak.');
      return {
        kpis: [
          { key: 'amp_response', label: 'Peak cross-flow amplitude A/D', value: pkS.y.amp, unit: '-', status: pkS.y.amp < 0.1 ? 'ok' : pkS.y.amp < 0.3 ? 'warn' : 'bad', note: 'Half the peak-to-peak displacement over the evaluation window of the stage; above about 0.1 D indicates lock-in' },
          { key: 'Ur_peak', label: 'Reduced velocity of the largest response', value: pkS.Ur, unit: '-' },
          { key: 'f_ratio', label: 'Oscillation frequency / natural frequency at the peak', value: pkS.y.f / fn, unit: '-', note: 'Close to 1 in lock-in for heavy cylinders; added mass shifts it for light ones' },
          { key: 'St_shed', label: 'Shedding frequency f·D/U in the first stage', value: St0, unit: '-', note: 'About 0.16–0.19 for a fixed cylinder at Re 100–150 (Strouhal number)' },
          { key: 'lock_in', label: 'Lock-in at the peak (1 = yes)', value: lock ? 1 : 0, unit: '-', status: lock ? 'warn' : 'ok', note: 'Lift and motion at the same frequency with A/D above 0.1' },
          { key: 'CL_rms', label: 'RMS lift coefficient at the peak', value: pkS.cl.rms, unit: '-', note: 'Fluctuating part, box-filtered over about one D/U to remove interface grid noise' }, { key: 'CD_mean', label: 'Mean drag coefficient at the peak', value: pkS.cd.mean, unit: '-', note: 'Over-predicted on coarse grids with wall blockage' },
          { key: 'f_shed_Hz', label: 'Shedding frequency at the stated size and speed', value: St0 / tS, unit: 'Hz' },
          { key: 'f_n_peak_Hz', label: 'Natural frequency that gives the largest response', value: fnHz, unit: 'Hz' },
          { key: 'amp_peak_m', label: 'Peak amplitude at the stated size', value: pkS.y.amp * i.L_ref, unit: 'm' },
          { key: 'lift_amp_Npm', label: 'Lift amplitude per unit length at the peak', value: pkS.cl.amp * 0.5 * fS, unit: 'N/m' },
          { key: 'scruton', label: 'Mass-damping (Scruton) number π²·m*·ζ', value: Math.PI ** 2 * i.m_star * i.zeta, unit: '-', note: `Mass per unit length ${mPhys.toPrecision(3)} kg/m` },
          ...common,
        ],
        plots: [
          vort,
          { type: 'line', title: 'Cylinder displacement history', xlabel: 'Time t·U/D [-]', ylabel: 'Displacement / D [-]', series: [{ name: 'Cross-flow y/D', x: tt, y: thin(R.Y) }, ...(R.nd > 1 ? [{ name: 'In-line x/D', x: tt, y: thin(R.X) }] : [])], annotations: ann },
          { type: 'line', title: 'Force coefficient histories', xlabel: 'Time t·U/D [-]', ylabel: 'Force coefficient [-]', series: [{ name: 'Lift', x: tt, y: thin(R.CL).map((v) => N.clamp(v, -6, 6)) }, { name: 'Drag', x: tt, y: thin(R.CD).map((v) => N.clamp(v, -6, 6)) }] },
          { type: 'line', title: 'Lock-in: amplitude against reduced velocity', xlabel: 'Reduced velocity U/(f_n·D) [-]', ylabel: 'Amplitude A/D [-]', series: [{ name: 'Cross-flow amplitude', x: st.map((g) => g.Ur), y: st.map((g) => g.y.amp), style: 'line+points' }] },
          { type: 'line', title: 'Frequencies against reduced velocity', xlabel: 'Reduced velocity U/(f_n·D) [-]', ylabel: 'Frequency f·D/U [-]', series: [{ name: 'Body motion', x: st.map((g) => g.Ur), y: st.map((g) => g.y.f), style: 'line+points' }, { name: 'Lift (shedding)', x: st.map((g) => g.Ur), y: st.map((g) => g.cl.f), style: 'points' }, { name: 'Natural frequency 1/Ur', x: st.map((g) => g.Ur), y: st.map((g) => 1 / g.Ur), style: 'dash' }] },
          { type: 'line', title: 'Spectra at the stage of largest response', xlabel: 'Frequency f·D/U [-]', ylabel: 'Amplitude [-]', series: [{ name: 'Displacement y/D', x: sy.f, y: sy.a }, { name: 'Lift coefficient', x: sl.f, y: sl.a }] },
        ],
        tables: [{ title: 'Response by stage', columns: ['Reduced velocity', 'A/D', 'f_motion·D/U', 'f_lift·D/U', 'CL rms', 'CD mean'], rows: st.map((g) => [g.Ur, g.y.amp, g.y.f, g.cl.f, g.cl.rms, g.cd.mean]) }],
        warnings,
        models: IB_MODELS(R, i), assumptions: IB_ASSUME,
      };
    }
    // plate: growth of the pitch oscillation and the stage at which it becomes self-sustained
    const a0 = N.rad(i.alpha0_deg), onset = st.find((g) => g.x.amp > 2 * a0), kRed = pkS.x.f > 0 ? pkS.x.f * TAU : NaN, sec = { b: 1, a: i.a_ea, x_alpha: i.x_alpha, r_alpha: i.r_alpha, m: i.mu * Math.PI, f_h: i.w_ratio / TAU, f_a: 1 / TAU, zeta: i.zeta, cl: TAU };
    let UrPk = NaN, UrDiv = NaN; try { const p = sectionFlutter(sec, 1, 40); UrPk = p.flutter?.P ?? NaN; UrDiv = p.div; } catch { UrPk = NaN; }
    const bTrue = i.L_ref / 2, waHz = (g) => i.U_ms / (g.Ur * bTrue) / TAU;
    if (!onset) warnings.push('The pitch oscillation did not grow to twice the initial disturbance in any stage: the plate is stable over this reduced-velocity range within the simulated time.');
    if (pkS.x.amp > 0.35) warnings.push('Pitch amplitude exceeds about 20°: the flow is separated (stall flutter / limit cycle) and the linearised inertial coupling of the two-degree-of-freedom section is only approximate.');
    return {
      kpis: [
        { key: 'amp_response', label: 'Largest pitch amplitude', value: pkS.x.amp, unit: 'rad', status: pkS.x.amp < 2 * a0 ? 'ok' : pkS.x.amp < 0.35 ? 'warn' : 'bad', note: `Initial disturbance ${a0.toPrecision(2)} rad; more than twice that is counted as self-excited` },
        { key: 'pitch_amp_deg', label: 'Largest pitch amplitude', value: N.deg(pkS.x.amp), unit: 'deg' },
        { key: 'plunge_amp_b', label: 'Plunge amplitude at that stage / semi-chord', value: pkS.y.amp, unit: '-' },
        { key: 'Ur_onset', label: 'Reduced velocity U/(b·ω_α) at which the motion becomes self-excited', value: onset ? onset.Ur : Infinity, unit: '-', note: 'First stage with pitch amplitude above twice the disturbance (resolution: the stage spacing)' },
        { key: 'Ur_flutter_pk', label: 'Flutter reduced velocity, inviscid p-k (Theodorsen) for the same section', value: UrPk, unit: '-', note: 'Thin-aerofoil potential flow; the low-Reynolds-number viscous flow has a lower lift slope and separates' },
        { key: 'Ur_div_pk', label: 'Divergence reduced velocity, inviscid', value: UrDiv, unit: '-' },
        { key: 'k_reduced', label: 'Reduced frequency ω·b/U of the pitch motion', value: kRed, unit: '-' },
        { key: 'f_ratio', label: 'Pitch oscillation frequency / pitch natural frequency', value: kRed * pkS.Ur, unit: '-' },
        { key: 'CL_rms', label: 'RMS lift coefficient at the largest response', value: pkS.cl.rms, unit: '-' }, { key: 'CD_mean', label: 'Mean drag coefficient at the largest response', value: pkS.cd.mean, unit: '-' },
        { key: 'f_alpha_onset_Hz', label: 'Pitch natural frequency at onset for the stated size and speed', value: onset ? waHz(onset) : NaN, unit: 'Hz' },
        { key: 'f_osc_Hz', label: 'Oscillation frequency at the stated size and speed', value: pkS.x.f > 0 ? pkS.x.f / (bTrue / i.U_ms) : NaN, unit: 'Hz' },
        ...common,
      ].filter((k) => fin(k.value) || k.key === 'Ur_onset'),
      plots: [
        vort,
        { type: 'line', title: 'Plate pitch history', xlabel: 'Time t·U/b [-]', ylabel: 'Pitch angle [deg]', series: [{ name: 'Pitch (nose-up)', x: tt, y: thin(R.X).map(N.deg) }], annotations: ann },
        { type: 'line', title: 'Plate plunge history', xlabel: 'Time t·U/b [-]', ylabel: 'Displacement / semi-chord [-]', series: [{ name: 'Elastic-axis height', x: tt, y: thin(R.Y) }], annotations: ann },
        { type: 'line', title: 'Force coefficient histories', xlabel: 'Time t·U/b [-]', ylabel: 'Coefficient [-]', series: [{ name: 'Lift', x: tt, y: thin(R.CL).map((v) => N.clamp(v, -6, 6)) }, { name: 'Drag', x: tt, y: thin(R.CD).map((v) => N.clamp(v, -6, 6)) }, { name: 'Moment about the elastic axis (nose-up)', x: tt, y: thin(R.CM).map((v) => N.clamp(v, -6, 6)) }] },
        { type: 'line', title: 'Pitch amplitude against reduced velocity', xlabel: 'Reduced velocity U/(b·ω_α) [-]', ylabel: 'Pitch amplitude [deg]', series: [{ name: 'Navier–Stokes coupled', x: st.map((g) => g.Ur), y: st.map((g) => N.deg(g.x.amp)), style: 'line+points' }], annotations: [{ y: N.deg(2 * a0), label: 'Twice the disturbance' }, ...(fin(UrPk) ? [{ x: UrPk, label: 'Inviscid flutter' }] : [])] },
        { type: 'line', title: 'Spectra at the stage of largest response', xlabel: 'Frequency f·b/U [-]', ylabel: 'Amplitude [-]', series: [{ name: 'Pitch angle [rad]', x: sy.f, y: sy.a }, { name: 'Lift coefficient', x: sl.f, y: sl.a }] },
      ],
      tables: [{ title: 'Response by stage', columns: ['Reduced velocity', 'Pitch amplitude [deg]', 'Plunge amplitude / b', 'f_pitch·b/U', 'CL rms'], rows: st.map((g) => [g.Ur, N.deg(g.x.amp), g.y.amp, g.x.f, g.cl.rms]) }],
      warnings,
      models: IB_MODELS(R, i), assumptions: [...IB_ASSUME, 'Plate: rigid section on linear plunge and pitch springs with linearised inertial coupling; rounded edges and finite thickness of at least two grid cells'],
    };
  },
  convergence: { param: 'cells', label: 'Grid cells per reference length', levels: [4, 5, 6], metric: 'CD_mean' },
  calibration: { params: [{ key: 'zeta', min: 0, max: 0.2 }, { key: 'm_star', min: 0.3, max: 1000 }], sweep: 'Ur_min', target: 'amp_response', note: 'Measured vibration amplitudes against reduced velocity from a water-channel or wind-tunnel test calibrate the structural damping and effective mass ratio.' },
  verify() {
    // fixed cylinder: Strouhal number, exact mass conservation, and the potential-flow added mass through the coupling terms
    const b = { ...IB_BASE, Re: 150, cells: 5, len: 16, height: 6.4, x_body: 4, n_sweep: 1, t_start: 20, t_stage: 22 }, R = ibRun(b, null, { fixed: true }), St = R.st[0].cl.f;
    // heavy cylinder released on a spring in the wake-free start: with a very large mass ratio the body ignores the fluid and oscillates at f_n
    const hv = ibRun({ ...IB_BASE, m_star: 400, zeta: 0, n_sweep: 1, Ur_min: 5, t_start: 4, t_stage: 16 }, null), fh = hv.st[0].y.f;
    // light cylinder (m* = 1) released with a velocity kick in fluid at rest: the strongly coupled period shows the added mass
    const sf = ibRun({ ...IB_BASE, m_star: 1, zeta: 0, Re: 2000, cells: 8, len: 8, height: 8, n_sweep: 1, Ur_min: 5, t_start: 2, t_stage: 13 }, null, { still: true }), ex = [];
    for (let k = 1; k < sf.Y.length - 1 && ex.length < 4; k++) if ((sf.Y[k] - sf.Y[k - 1]) * (sf.Y[k + 1] - sf.Y[k]) < 0) ex.push(sf.T[k]);
    const fStill = ex.length > 2 ? (ex.length - 1) / (2 * (ex[ex.length - 1] - ex[0])) : NaN;
    // coupling algebra alone: the flow sub-solve is replaced by a prescribed added-mass reaction H = −m_a·Δv/Δt, so the coupled system is
    // a linear oscillator of mass m + m_a whose average-acceleration solution obeys an exact linear recurrence
    const cy = { ...IB_BASE, m_star: 1, zeta: 0.02, cells: 3, len: 6, height: 3, n_sweep: 1, Ur_min: 5, t_start: 2, t_stage: 20, cfl: 0.05, tol: 1e-13 }, mc = Math.PI / 4, ma = 1.5 * mc, wn0 = TAU / cy.Ur_min;
    const am = (vn, vo, dt) => [0, (-ma * (vn[0] - vo[0])) / dt, 0], c1 = ibRun(cy, null, { still: true, load: am }), c3 = ibRun({ ...cy, subiter: 3 }, null, { still: true, load: am }), Yc = c1.Y, hs = c1.dt / 2;
    const w2 = (mc * wn0 * wn0) / (mc + ma), cz = (2 * cy.zeta * mc * wn0) / (mc + ma), den = 1 + cz * hs + w2 * hs * hs, trA = (2 * (2 + cz * hs)) / den - 2, detA = (1 - cz * hs + w2 * hs * hs) / den, wd = Math.sqrt(w2 - 0.25 * cz * cz), yA = 0.2 / wd;
    let rec1 = 0, dev = 0, d13 = 0; for (let k = 1; k + 1 < Yc.length; k++) rec1 = Math.max(rec1, Math.abs(Yc[k + 1] - trA * Yc[k] + detA * Yc[k - 1]) / yA);
    for (let k = 0; k < Yc.length; k++) { dev = Math.max(dev, Math.abs(Yc[k] - yA * Math.exp(-0.5 * cz * c1.T[k]) * Math.sin(wd * c1.T[k])) / yA); d13 = Math.max(d13, Math.abs(c3.Y[k] - Yc[k]) / yA); }
    // plate: two degrees of freedom with inertial coupling and a full added-mass matrix [[π, −πa], [−πa, π(1/8 + a²)]]; the free response is a
    // sum of two modes, so y(k+2) + y(k−2) − 2(c1 + c2)·(y(k+1) + y(k−1)) + (2 + 4·c1·c2)·y(k) = 0 with c = cos of the discrete modal phase step
    const pb = { ...IB_BASE, body: IB_BODIES[1], mu: 2, a_ea: -0.2, x_alpha: 0.2, r_alpha: 0.5, w_ratio: 0.6, zeta: 0, alpha0_deg: 4, cells: 3, len: 8, height: 6, n_sweep: 1, Ur_min: 1, t_start: 2, t_stage: 20, cfl: 0.1, tol: 1e-13 };
    const A2 = [[Math.PI, -Math.PI * pb.a_ea], [-Math.PI * pb.a_ea, Math.PI * (0.125 + pb.a_ea ** 2)]], amP = (vn, vo, dt) => [0, (A2[0][0] * (vn[0] - vo[0]) + A2[0][1] * (vn[1] - vo[1])) / dt, (A2[1][0] * (vn[0] - vo[0]) + A2[1][1] * (vn[1] - vo[1])) / dt];
    const p1 = ibRun(pb, null, { still: true, load: amP }), mp = pb.mu * Math.PI, Ip = mp * pb.r_alpha ** 2, m11 = mp + A2[0][0], m12 = mp * pb.x_alpha + A2[0][1], m22 = Ip + A2[1][1], k1 = mp * pb.w_ratio ** 2, k2 = Ip;
    const qa = m11 * m22 - m12 * m12, qb = -(k1 * m22 + k2 * m11), disc = Math.sqrt(qb * qb - 4 * qa * k1 * k2), cosd = [(-qb - disc) / (2 * qa), (-qb + disc) / (2 * qa)].map((l) => { const x = 0.25 * l * p1.dt * p1.dt; return (1 - x) / (1 + x); });
    let rec2 = 0; for (const sig of [p1.X, p1.Y]) { const sc = N.amax(sig.map(Math.abs)); for (let k = 2; k + 2 < sig.length; k++) rec2 = Math.max(rec2, Math.abs(sig[k + 2] + sig[k - 2] - 2 * (cosd[0] + cosd[1]) * (sig[k + 1] + sig[k - 1]) + (2 + 4 * cosd[0] * cosd[1]) * sig[k]) / sc); }
    return [
      N.check('Coupling with a prescribed added-mass reaction: damped oscillator of mass m + m_a, exact recurrence of the average-acceleration rule', rec1, 0, 1e-10, 'y(k+1) − tr(A)·y(k) + det(A)·y(k−1) = 0 with A = (I − ½Δt·B)⁻¹(I + ½Δt·B), B the state matrix of (m + m_a)ÿ + cẏ + ky = 0; any staggering between fluid load and structural equation would show as numerical damping'),
      N.check('Coupling with a prescribed added-mass reaction: response against the analytic damped oscillator with the added-mass frequency shift', dev, 0, 1e-3, 'y = (v0/ω_d)·e^(−ζωt)·sin(ω_d·t), ω² = k/(m + m_a) with m_a = 1.5 m; the phase error of the second-order rule, (ωΔt)²/12·ωt, is about 3e-4 at the end of the run'),
      N.check('Coupling iteration: one Newton step solves the interface equation (iterations per step with a limit of three)', c3.its, 2, 1e-12, 'The sub-step is affine in the body acceleration, so Newton with its exact operator converges in one step; the second iteration only confirms a zero residual'),
      N.check('Coupling iteration: one and three iterations per step give the same response', d13, 0, 1e-10, 'Exactness of the single Newton step'),
      N.check('Coupling residual after convergence', c3.resMax, 0, 1e-13, 'Body position used by the flow minus the structural position'),
      N.check('Two-degree-of-freedom plate with a prescribed full added-mass matrix: exact two-mode recurrence', rec2, 0, 1e-9, 'Generalised eigenvalues of (M + M_a, K) mapped through the average-acceleration rule, cos(ω_d·Δt) = (1 − ¼ω²Δt²)/(1 + ¼ω²Δt²); verifies the 2 × 2 interface solve with off-diagonal structural and added mass'),
      N.check('Light cylinder in still fluid: frequency f_n·√(m/(m + m_a)) with m_a = ρπD²/4', fStill, 0.2 * Math.sqrt(0.5), 0.12, 'Potential-flow added mass of a circular cylinder; with 8 cells per diameter the smeared interface, the 12.5% confinement and the unresolved Stokes layer raise the effective added mass to about 1.3, so the frequency is 7–8% low (it approaches the exact value as the grid is refined)'),
      N.check('Fixed-cylinder Strouhal number at Re = 150', St, 0.183, 0.1, 'Williamson’s St–Re relation (0.183); deliberately coarse grid of 5 cells per diameter with 16% blockage, hence the 10% tolerance'),
      N.check('Velocity field is divergence-free after projection', R.divMax, 0, 1e-9, 'Exact discrete projection (cosine transform + tridiagonal solve)'),
      N.check('Very heavy spring-mounted cylinder oscillates at its natural frequency', fh, 1 / 5, 0.03, 'Limit of infinite mass ratio: fluid loads are negligible, Newmark average acceleration with about 80 steps per period'),
      N.check('Mean drag of the fixed cylinder is positive and of order one', R.st[0].cd.mean > 0.8 && R.st[0].cd.mean < 3 ? 1 : 0, 1, 1e-12, 'Sanity bound: CD ≈ 1.3 at Re 150 in an unbounded stream; blockage and the diffuse interface raise it on this grid'),
    ];
  },
  recommend(res, i) {
    const o = res.outputs, out = [], plate = i.body === IB_BODIES[1];
    if (!plate && o.lock_in > 0) out.push({ severity: o.amp_response > 0.3 ? 'warn' : 'advise', title: 'Vortex-induced vibration with lock-in', detail: `Amplitude ${o.amp_response.toFixed(2)} D at reduced velocity ${o.Ur_peak.toFixed(1)}; shedding at about ${o.f_shed_Hz.toFixed(0)} Hz for the stated size and speed.`, action: `Keep the natural frequency of the tube or strut away from the shedding frequency over the flight envelope (reduced velocity below about 3 or above about 10), raise the mass-damping number (now ${o.scruton.toFixed(2)}), or break up the shedding with a fairing, strakes or a splitter plate.`, basis: 'Lock-in when the shedding frequency approaches the structural frequency; fatigue from sustained oscillation' });
    if (!plate && !(o.lock_in > 0)) out.push({ severity: 'info', title: 'No lock-in in the simulated range', detail: `Largest amplitude ${o.amp_response.toFixed(3)} D.`, action: 'Widen the reduced-velocity range and lengthen the stages before concluding: resonance builds over tens of cycles.', basis: 'Vortex-induced vibration screening' });
    if (plate && fin(o.Ur_onset)) out.push({ severity: 'advise', title: 'Self-excited pitch–plunge oscillation of the model plate', detail: `Onset at reduced velocity ${o.Ur_onset.toFixed(1)} with pitch amplitude up to ${o.pitch_amp_deg.toFixed(0)}°${fin(o.Ur_flutter_pk) ? `; inviscid flutter estimate ${o.Ur_flutter_pk.toFixed(1)}` : ''}.`, action: 'Raise the torsional frequency, move the mass centre forward of the elastic axis or add damping; confirm the boundary with the multi-mode wing flutter and coupled lattice–beam analyses at the real Reynolds number. The plate is a generic light section (mass ratio, frequency ratio and reduced-velocity range as entered), not this aircraft’s wing: the sweep is set up to cross its flutter boundary.', basis: 'Time-domain growth of the coupled Navier–Stokes / two-degree-of-freedom section; demonstration case, not a clearance of the aircraft' });
    if (plate && !fin(o.Ur_onset)) out.push({ severity: 'info', title: 'Plate stable in the simulated range', detail: `Pitch amplitude stayed at ${o.pitch_amp_deg.toFixed(1)}°.`, action: 'Extend the reduced-velocity range or the time per stage to find the onset.', basis: 'Time-domain response' });
    out.push({ severity: 'advise', title: 'Resolution of the flow solution', detail: `${o.n_cells} cells, ${o.n_steps} steps, ${i.cells} cells per reference length.`, action: 'For quantitative amplitudes run the convergence study and use at least 16 cells per reference length; for flight Reynolds numbers export the case from the High-fidelity bridge page to a body-fitted solver.', basis: 'Grid-convergence practice for immersed-boundary methods' });
    return out;
  },
};
const IB_MODELS = (R, i) => [`Incompressible Navier–Stokes on a ${R.nx} × ${R.ny} staggered grid, fractional-step projection with an exact transform-based pressure solve`, 'Third-order upwind-biased advection (fourth-order central part by Adams–Bashforth, dissipative part by forward Euler), explicit diffusion', 'Direct-forcing immersed boundary with a one-cell smoothed body mask; loads from the momentum exchanged plus the inertia of the fluid inside the mask', `Rigid body on linear springs and dampers, average-acceleration rule in time-centred (impulse) form, which pairs with the step-mean fluid load without a half-step lag; ${R.nNewton ? `strong coupling: force–project–force fluid sub-solves give the load and, by perturbing each degree of freedom, its discrete added-mass operator; the interface condition is solved implicitly with it and the flow corrected by superposition (stable for mass ratios below one); ${R.nNewton > 1 ? `operator rebuilt in each of up to ${R.nNewton} iterations per step` : 'operator refreshed every fourth step (quasi-Newton)'}` : 'loose (staggered) coupling with a lagged analytical added-mass term'}`, 'Stepped reduced-velocity sweep in one continuous simulation', 'CFD-based aeroelastic model; limit-cycle oscillation model'];
const IB_ASSUME = ['Two-dimensional, laminar, incompressible flow at low Reynolds number', 'Uniform inflow, convective outflow, free-slip walls at the top and bottom (blockage raises forces and shedding frequency)', 'First-order accurate interface treatment on a Cartesian grid: the body surface is smeared over one cell', 'Structure moves as a rigid body; springs are linear'];

export default {
  id: 'aeroelastic', n: 3,
  tagline: 'At what speed the structure and the airflow start feeding each other — flutter, divergence, control reversal — and how much margin remains.',
  analyses: [rotor, rotorMap, section, wingStatic, wingFlutter, gust, panel, uvlmfsi, cfdfsi],
  consumes: [
    { from: 'fea', keys: ['EI_root_Nm2', 'GJ_root_Nm2', 'wing_struct_mass_kg'], why: 'Wing stiffness and mass' },
    { from: 'vibration', keys: ['f1_Hz', 'f_torsion_Hz'], why: 'Bending and torsion frequencies of the typical section' },
    { from: 'cfd', keys: ['CLa_per_rad'], why: 'Lift-curve slope for the strip aerodynamics' },
    { from: 'performance', keys: ['V_d_eas'], why: 'Design dive speed for the flutter margin' },
  ],
  provides: [
    { key: 'V_flutter_ms', label: 'Flutter speed', unit: 'm/s' }, { key: 'f_flutter_Hz', label: 'Flutter frequency', unit: 'Hz' }, { key: 'V_divergence_ms', label: 'Divergence speed', unit: 'm/s' },
    { key: 'V_reversal_ms', label: 'Aileron reversal speed', unit: 'm/s' }, { key: 'flutter_margin', label: 'Flutter margin over 1.15 VD', unit: '-' },
  ],
  handoff: [
    { model: 'Three-dimensional body-fitted transonic CFD–CSD (Euler / RANS, Arbitrary Lagrangian–Eulerian moving meshes, monolithic coupling)', why: 'Coupled time-domain fluid–structure solutions are now native in two forms — an unsteady vortex lattice strongly coupled to a finite-element beam wing, and a two-dimensional incompressible Navier–Stokes immersed-boundary solver coupled to a spring-mounted body — but shock motion, the transonic flutter dip, separation at flight Reynolds number and full-aircraft geometry need a compressible body-fitted solver on millions of cells', tool: 'Coupled CFD–CSD solver; the app’s High-fidelity bridge page exports a ready-to-run case (geometry, structural model, flight condition) for it' },
    { model: 'Doublet-lattice unsteady aerodynamics on the full aircraft', why: 'Needs a frequency-domain lifting-surface model of wing, tail and fuselage with splining to a full structural finite-element model; the native lattice covers one cantilever wing in incompressible flow', tool: 'Aeroelastic FE solver with doublet-lattice (SOL 145-class); case export from the High-fidelity bridge page' },
    { model: 'Free-wake roll-up, large-deflection (geometrically nonlinear) aeroelasticity and free-flying flexible aircraft', why: 'The native lattice keeps a flat, free-stream-convected wake and small displacements on a clamped wing; very flexible wings and coupled rigid-body motion are not represented', tool: 'Nonlinear aeroelastic simulation framework (geometrically exact beam + free-wake vortex lattice)' },
    { model: 'Control-surface flutter, freeplay and aeroservoelasticity', why: 'Needs Theodorsen flap functions or panel aerodynamics, actuator impedance and the control laws', tool: 'Aeroservoelastic analysis package' },
    { model: 'Whirl flutter of propeller–nacelle and tiltrotor systems', why: 'Requires validated propeller aerodynamic derivatives and nacelle mount dynamics', tool: 'Whirl-flutter module of an aeroelastic or rotorcraft comprehensive code' },
    { model: 'Comprehensive rotor aeroelasticity (geometrically exact blades, flap–lag–torsion, forward flight, Loewy/free wake)', why: 'Periodic-coefficient multi-blade analysis with trim and wake models', tool: 'Rotorcraft comprehensive analysis code' },
    { model: 'Fluid–structure–thermal interaction', why: 'Aerothermal heating coupled with structural response needs CFD heat-flux and thermal FE', tool: 'Coupled aerothermoelastic toolchain' },
  ],
};
