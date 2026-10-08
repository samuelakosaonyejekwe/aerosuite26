// Suite 3 — Aeroelasticity and Fluid–Structure Interaction.
// One unsteady strip-theory kernel (Theodorsen in the frequency domain, Wagner/Küssner indicial functions in time)
// drives a typical section, a Rayleigh–Ritz bending–torsion wing and a rotating pitch–flap blade: flutter by the
// p-k method, divergence, aileron reversal, load redistribution, limit cycles, gust response and panel flutter.

import * as N from '../core/numerics.js';
import { isa, G0, RHO0 } from '../core/atmosphere.js';
import { METALS } from '../data/materials.js';

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

/** Allocation-free complex eigenvalue solver (shifted QR on the Hessenberg form), equivalent to the kernel's eig(). */
function ceig(A) {
  const n = A.length, hr = new Float64Array(n * n), hi = new Float64Array(n * n), cr = new Float64Array(n), ci = new Float64Array(n), sr = new Float64Array(n), si = new Float64Array(n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { const v = A[i][j]; if (typeof v === 'number') hr[i * n + j] = v; else { hr[i * n + j] = v[0]; hi[i * n + j] = v[1]; } }
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
/** Roots of [s²(M+Ma) + s(C+Ca+ΣC(k)·Fc·Wdᵀ) + K + ΣC(k)·Fc·Wkᵀ] at the oscillation frequency w (quasi → C(k) = 1). */
function pRoots(m, rho, w, pre, quasi) {
  const n = m.n, Ct = m.C.map((r, i) => r.map((v, j) => [v + pre.ap.Ca[i][j], 0])), Kt = m.K.map((r) => r.map((v) => [v, 0]));
  for (const v of pre.sv) {
    const ck = quasi ? [1, 0] : theodorsen((w * v.s.b) / Math.max(v.s.U, 1e-9));
    for (let i = 0; i < n; i++) { const f = v.Fc[i]; if (f) for (let j = 0; j < n; j++) { Ct[i][j][0] += ck[0] * f * v.Wd[j]; Ct[i][j][1] += ck[1] * f * v.Wd[j]; Kt[i][j][0] += ck[0] * f * v.Wk[j]; Kt[i][j][1] += ck[1] * f * v.Wk[j]; } }
  }
  const A = N.range(2 * n, () => new Array(2 * n).fill(0));
  for (let i = 0; i < n; i++) {
    A[i][i + n] = 1;
    for (let j = 0; j < n; j++) { let kr = 0, ki = 0, c1 = 0, c2 = 0; for (let l = 0; l < n; l++) { const mi = pre.Mi[i][l]; kr += mi * Kt[l][j][0]; ki += mi * Kt[l][j][1]; c1 += mi * Ct[l][j][0]; c2 += mi * Ct[l][j][1]; } A[i + n][j] = [-kr, -ki]; A[i + n][j + n] = [-c1, -c2]; }
  }
  return ceig(A);
}
const prep = (m, rho) => { const ap = apparent(m, rho); return { ap, sv: stripVecs(m, rho), Mi: N.inv(N.madd(m.M, ap.Ma)) }; };
/** p-k iteration for one mode at one condition, started from the root guess [σ, ω]. */
function pkPoint(m, rho, guess, quasi) {
  const pre = prep(m, rho); let lam = guess;
  for (let it = 0; it < 40; it++) {
    const r = pRoots(m, rho, Math.abs(lam[1]), pre, quasi); let best = null, bd = Infinity;
    for (const x of r) { if (x[1] < -1e-9) continue; const d = Math.hypot(x[0] - lam[0], x[1] - lam[1]); if (d < bd) { bd = d; best = x; } }
    if (!best) break;
    const dw = Math.abs(best[1] - lam[1]); lam = best;
    if (dw < 1e-7 * (1 + Math.abs(lam[1])) || quasi) break;
  }
  return lam;
}
/** Static aeroelastic stiffness per unit (reference speed)²: K̂ = Σ Fc·Wkᵀ evaluated by the caller's model. */
function aeroStiff(m, rho) { const n = m.n, Kh = N.zeros(n); for (const v of stripVecs(m, rho)) for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Kh[i][j] += v.Fc[i] * v.Wk[j]; return Kh; }
/**
 * p-k sweep over a parameter grid P (airspeed or rotor speed). mk(P) builds the model. Returns per-mode root tracks,
 * the first flutter crossing (refined by bisection) and the divergence point from det(K + K_aero) = 0.
 */
function pkSweep(mk, rho, Ps, w0, quasi) {
  const nm = w0.length, tracks = w0.map(() => []); let last = w0.map((w) => [0, w]);
  for (const P of Ps) { const m = mk(P); last = last.map((g) => pkPoint(m, rho, g, quasi)); last.forEach((l, k) => tracks[k].push(l)); }
  let flutter = null;
  for (let k = 0; k < nm; k++) for (let j = 1; j < Ps.length; j++) {
    const a = tracks[k][j - 1], b = tracks[k][j];
    if (a[0] <= 0 && b[0] > 0 && b[1] > 1e-6 && (!flutter || Ps[j] < flutter.hi)) {
      let lo = Ps[j - 1], hi = Ps[j], g = a, root = b;
      for (let it = 0; it < 30 && hi - lo > 1e-6 * hi; it++) { const mid = 0.5 * (lo + hi), r = pkPoint(mk(mid), rho, g, quasi); if (r[0] > 0) { hi = mid; root = r; } else { lo = mid; g = r; } }
      flutter = { P: 0.5 * (lo + hi), w: root[1], mode: k, hi: Ps[j] };
    }
  }
  const det = (P) => { const m = mk(P); return N.det(N.madd(m.K, aeroStiff(m, rho))) / N.det(m.K); };
  let div = NaN, d0 = det(Ps[0]);
  for (let j = 1; j < Ps.length; j++) { const d1 = det(Ps[j]); if (d0 > 0 && d1 <= 0) { div = N.brent(det, Ps[j - 1], Ps[j], 1e-7 * Ps[j]); break; } d0 = d1; }
  return { tracks, flutter, div };
}
/** Time-domain state matrix with two Jones lag states per strip: x = [q, q̇, y1, y2]. Also returns the layout for forcing. */
function stateSpace(m, rho) {
  const n = m.n, pre = prep(m, rho), ns = pre.sv.length, nx = 2 * n + 2 * ns, A = N.zeros(nx), a0 = 1 - JONES.A[0] - JONES.A[1];
  const Kq = m.K.map((r) => r.slice()), Cq = N.madd(m.C, pre.ap.Ca), Fy = N.zeros(n, 2 * ns);
  pre.sv.forEach((v, s) => {
    for (let i = 0; i < n; i++) { for (let j = 0; j < n; j++) { Kq[i][j] += a0 * v.Fc[i] * v.Wk[j]; Cq[i][j] += a0 * v.Fc[i] * v.Wd[j]; } Fy[i][2 * s] = Fy[i][2 * s + 1] = v.Fc[i]; }
    for (let l = 0; l < 2; l++) { const r = 2 * n + 2 * s + l, be = (JONES.b[l] * v.s.U) / v.s.b; A[r][r] = -be; for (let j = 0; j < n; j++) { A[r][j] = be * JONES.A[l] * v.Wk[j]; A[r][n + j] = be * JONES.A[l] * v.Wd[j]; } }
  });
  const MK = N.matmul(pre.Mi, Kq), MC = N.matmul(pre.Mi, Cq), MF = N.matmul(pre.Mi, Fy);
  for (let i = 0; i < n; i++) { A[i][n + i] = 1; for (let j = 0; j < n; j++) { A[n + i][j] = -MK[i][j]; A[n + i][n + j] = -MC[i][j]; } for (let j = 0; j < 2 * ns; j++) A[n + i][2 * n + j] = -MF[i][j]; }
  return { A, Mi: pre.Mi, n, nx };
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
  const tors = eg.vectors.map((v) => { let t = 0, tot = 0; for (let p = 0; p < n; p++) for (let q = 0; q < n; q++) { const e = v[p] * M[p][q] * v[q]; tot += e; if (p >= nb && q >= nb) t += e; } return t / tot; });
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
    return { b: ch / 2, a: 2 * wi.x_ea - 1, x_alpha: 2 * (wi.x_cg - wi.x_ea), r_alpha: 2 * wi.r_alpha, m: m0 * (1 - (1 - wi.taper) * 0.75) ** 2, f_h: up.vibration?.f1_Hz ?? r.wn[kb] / TAU, f_a: up.vibration?.f_torsion_Hz ?? r.wn[kt] / TAU, cl: w.cl, alt_m: w.alt_m, V_d_eas: w.V_d_eas };
  },
  run(i) {
    const at = isa(i.alt_m), rho = at.rho, r = sectionFlutter(i, rho, Math.round(i.nSpeeds)), s = r.s, Uf = r.flutter?.P ?? NaN, wf = r.flutter?.w ?? NaN, sg = Math.sqrt(at.sigma), mu = s.m / (Math.PI * rho * s.b ** 2);
    const names = ['Bending branch', 'Torsion branch'], warnings = [];
    // time histories and limit-cycle amplitude sweep above the linear flutter speed
    const U1 = fin(Uf) ? i.V_lco_frac * Uf : r.Us[r.Us.length - 1] * 0.5, a0 = N.rad(i.alpha0_deg), th = lcoRun(s, rho, U1, i.kappa, a0, 60, Math.round(i.nSteps));
    const fr = [1.02, 1.05, 1.1, 1.15, 1.2, 1.3], lco = fin(Uf) && i.kappa > 0 ? fr.map((x) => lcoRun(s, rho, x * Uf, i.kappa, a0, 60, Math.round(i.nSteps)).amp) : [];
    const Me = Uf / at.a;
    if (!fin(Uf)) warnings.push('No flutter crossing was found up to the highest speed analysed.');
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
      warnings, models: ['Linear flutter model: two-degree-of-freedom typical section', 'p–k flutter model with exact Theodorsen function (Bessel functions)', 'Non-linear time-domain flutter model: Jones two-lag approximation of the Wagner function, RK4', 'Limit-cycle oscillation model: cubic hardening spring'],
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
function wingFlutterSolve(i, rho, nm, nSp, Vtop) {
  const W = wingRitz(i, nm, nm, 14), track = W.wn.map((_, k) => k).slice(0, Math.min(4, W.n)), Vs = N.linspace(Vtop / nSp, Vtop, nSp);
  return { W, Vs, track, ...pkSweep(W.at, rho, Vs, track.map((k) => W.wn[k])) };
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
    const alts = [0, 3000, 6000, 9000, 12000].filter((h) => h <= Math.max(i.alt_m, 3000) + 3001), env = alts.map((h) => { const a = isa(h), x = wingFlutterSolve(i, a.rho, Math.min(nm, 2), 16, Math.max(3 * i.V_d_eas / Math.sqrt(a.sigma), 2 * a.a)); return (x.flutter?.P ?? NaN) * Math.sqrt(a.sigma); });
    const kb = W.tors.findIndex((t) => t < 0.5), kt = W.tors.findIndex((t) => t >= 0.5), envOk = alts.map((_, k) => k).filter((k) => fin(env[k]));
    if (envOk.length < alts.length) warnings.push(envOk.length ? 'No flutter was found inside the coarse altitude sweep at some altitudes; those points are left off the flutter margin diagram.' : 'No flutter was found inside the coarse altitude sweep, so the flutter margin diagram shows only the VD lines: the boundary lies above three times VD.');
    if (!fin(Vf)) warnings.push('No flutter crossing was found in the speed range analysed (three times the default range).');
    if (fin(Vf) && Vf / at.a > 0.7) warnings.push(`Flutter Mach number ${(Vf / at.a).toFixed(2)} is beyond the validity of incompressible strip theory: the transonic dip can lower the true flutter speed substantially.`);
    if (Math.abs(i.sweep_deg) > 15) warnings.push('Swept wing: only the bending-slope incidence term of swept strip theory is retained; expect larger uncertainty.');
    const lbl = (k) => `Mode ${k + 1} (${W.tors[r.track[k]] >= 0.5 ? 'torsion' : 'bending'}, ${(W.wn[r.track[k]] / TAU).toPrecision(3)} Hz)`;
    return {
      kpis: [
        { key: 'V_flutter_ms', label: 'Flutter speed (TAS)', value: Vf, unit: 'm/s' },
        { key: 'V_flutter_eas', label: 'Flutter speed (EAS)', value: Vf * sg, unit: 'm/s' },
        { key: 'f_flutter_Hz', label: 'Flutter frequency', value: wf / TAU, unit: 'Hz' },
        { key: 'flutter_margin', label: 'Flutter margin', value: fin(mg) ? mg : 9.99, unit: '-', status: mg > 0 ? 'ok' : 'bad', note: MARGIN_NOTE + (fin(mg) ? '' : '; no flutter found, reported as the cap 9.99') },
        { key: 'M_flutter', label: 'Flutter Mach number (incompressible theory)', value: Vf / at.a, unit: '-', status: Vf / at.a < 0.7 ? 'ok' : 'warn' },
        { key: 'V_div_dyn_ms', label: 'Divergence speed from the same model', value: Vdiv, unit: 'm/s' },
        { key: 'f_bend1_Hz', label: 'First bending mode in vacuum', value: W.wn[kb] / TAU, unit: 'Hz' },
        { key: 'f_tors1_Hz', label: 'First torsion mode in vacuum', value: (W.wn[kt] ?? NaN) / TAU, unit: 'Hz' },
        { key: 'flutter_mode', label: 'Branch that goes unstable', value: r.flutter ? r.flutter.mode + 1 : 0, unit: '-' },
      ].filter((k) => fin(k.value) || ['V_flutter_ms', 'f_flutter_Hz'].includes(k.key)),
      plots: [
        { type: 'line', title: 'V–g diagram', xlabel: 'True airspeed [m/s]', ylabel: 'Damping g = 2σ/ω [-]', series: r.tracks.map((tk, k) => ({ name: lbl(k), x: r.Vs, y: tk.map((l) => (l[1] > 1e-6 ? N.clamp((2 * l[0]) / l[1], -1.5, 1.5) : NaN)) })), annotations: [{ y: 0, label: 'Flutter boundary' }, { x: (1.15 * i.V_d_eas) / sg, label: '1.15 VD' }] },
        { type: 'line', title: 'V–f diagram', xlabel: 'True airspeed [m/s]', ylabel: 'Frequency [Hz]', series: r.tracks.map((tk, k) => ({ name: lbl(k), x: r.Vs, y: tk.map((l) => l[1] / TAU) })) },
        { type: 'line', title: 'Flutter margin diagram', xlabel: 'Flutter speed [m/s EAS]', ylabel: 'Altitude [m]', series: [...(envOk.length ? [{ name: 'Flutter boundary', x: envOk.map((k) => env[k]), y: envOk.map((k) => alts[k]), style: 'line+points' }] : []), { name: '1.15 VD', x: alts.map(() => 1.15 * i.V_d_eas), y: alts, style: 'dash' }, { name: 'VD', x: alts.map(() => i.V_d_eas), y: alts, style: 'dash' }] },
      ],
      tables: [{ title: 'In-vacuum normal modes', columns: ['Mode', 'Frequency [Hz]', 'Torsional energy share'], rows: W.wn.slice(0, 6).map((w, k) => [k + 1, w / TAU, W.tors[k]]) }],
      warnings, models: ['p–k flutter model with Theodorsen strip aerodynamics at the local reduced frequency', 'Finite-span lift-curve-slope correction of the circulatory terms (modified strip theory)', 'Rayleigh–Ritz structural dynamics model with inertial bending–torsion coupling'],
      assumptions: ['Incompressible, attached flow; no tip-loss or spanwise aerodynamic coupling', 'Cantilever wing without engines, stores or control surfaces', 'Constant true-airspeed sweep at fixed density (matched-point iteration on Mach is not performed)'],
    };
  },
  convergence: { param: 'nModes', label: 'Assumed modes per motion', levels: [1, 2, 3, 4], metric: 'V_flutter_ms' },
  calibration: { params: [{ key: 'GJ_root', min: 1, max: 1e11 }, { key: 'EI_root', min: 1, max: 1e12 }, { key: 'x_cg', min: 0.2, max: 0.7 }, { key: 'zeta', min: 0, max: 0.1 }], sweep: 'alt_m', target: 'V_flutter_ms', note: 'Flutter-model wind-tunnel boundaries or flight flutter-test damping trends; update stiffness from the ground vibration test first' },
  verify() {
    // Goland wing: uniform cantilever, exact strip-theory flutter 137.2 m/s at 70.7 rad/s
    const g = { semi_span: 6.096, c_root: 1.8288, taper: 1, sweep_deg: 0, x_ea: 0.33, x_cg: 0.43, r_alpha: Math.sqrt(8.64 / 35.71) / 1.8288, EI_root: 9.77e6, GJ_root: 0.987e6, stiff_exp: 0, m_semi: 35.71 * 6.096, zeta: 0, cl: TAU };
    const r = wingFlutterSolve(g, 1.225, 3, 30, 220), k = Math.sqrt(9.77e6 / (35.71 * 6.096 ** 4));
    return [
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
    if (o.M_flutter > 0.7) out.push({ severity: 'advise', title: 'Transonic flutter assessment required', detail: `Predicted flutter Mach ${o.M_flutter.toFixed(2)}.`, action: 'Use compressible unsteady aerodynamics (doublet lattice with transonic correction or CFD-based aeroelasticity).', basis: 'Validity limit of incompressible theory' });
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
    { key: 'U_ds', label: 'Design gust velocity (EAS)', unit: 'm/s', default: 15.24, min: 0, max: 40, group: 'Gust', help: 'Derived gust velocity; 15.24 m/s (50 ft/s) at VC is the classical value' },
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
    if (nPeak > i.n_limit) warnings.push('The discrete gust load factor exceeds the manoeuvre limit: the structure is gust-critical at this condition.');
    if (i.V / at.a > 0.75) warnings.push('Compressibility is not included in the indicial functions; use Mach-dependent lift build-up above about Mach 0.7.');
    return {
      kpis: [
        { key: 'dn_gust', label: 'Peak incremental load factor (unsteady)', value: u.peak, unit: 'g' },
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
    if (o.n_gust_peak > i.n_limit) out.push({ severity: 'warn', title: 'Gust-critical condition', detail: `Peak ${o.n_gust_peak.toFixed(2)} g against the ${i.n_limit} g manoeuvre limit.`, action: 'Carry the gust load factor to Suites 2 and 9; consider gust-load alleviation in Suite 16 or a higher wing loading.', basis: 'Discrete tuned-gust response' });
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
    const op = bladeModel(i, Om0), lock = (rho * i.cl * i.chord * i.R ** 4) / op.Ib, rt = [pkPoint(op, rho, [0, op.wn[0]]), pkPoint(op, rho, [0, op.wn[1]])], mg = lim / (1.2 * Om0) - 1, warnings = [];
    if (i.x_cg > i.x_pa + 1e-9 && !fin(lim)) warnings.push('The mass centre is aft of the pitch axis but no instability was found in this speed range.');
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
      warnings, models: ['Rotor aeroelastic model: rigid flap and rigid pitch with centrifugal and propeller-moment coupling', 'p–k method with Theodorsen strip aerodynamics at the local reduced frequency', 'Hover, zero inflow perturbation'],
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

export default {
  id: 'aeroelastic', n: 3,
  tagline: 'At what speed the structure and the airflow start feeding each other — flutter, divergence, control reversal — and how much margin remains.',
  analyses: [rotor, rotorMap, section, wingStatic, wingFlutter, gust, panel],
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
    { model: 'Doublet-lattice unsteady aerodynamics on the full aircraft', why: 'Needs a lifting-surface panel model with aerodynamic influence coefficients and splining to a full structural model', tool: 'Aeroelastic FE solver with doublet-lattice (SOL 145-class)' },
    { model: 'Transonic aeroelasticity and CFD-based aeroelastic models', why: 'Shock motion and the transonic dip need time-accurate Euler/RANS coupled to the structure', tool: 'Coupled CFD–CSD solver' },
    { model: 'Partitioned / monolithic CFD–CSD coupling with ALE mesh motion and conservative load transfer', why: 'Three-dimensional moving-mesh flow solution is beyond in-browser computation', tool: 'CFD and FE codes with a coupling library' },
    { model: 'Control-surface flutter, freeplay and aeroservoelasticity', why: 'Needs Theodorsen flap functions or panel aerodynamics, actuator impedance and the control laws', tool: 'Aeroservoelastic analysis package' },
    { model: 'Whirl flutter of propeller–nacelle and tiltrotor systems', why: 'Requires validated propeller aerodynamic derivatives and nacelle mount dynamics', tool: 'Whirl-flutter module of an aeroelastic or rotorcraft comprehensive code' },
    { model: 'Comprehensive rotor aeroelasticity (geometrically exact blades, flap–lag–torsion, forward flight, Loewy/free wake)', why: 'Periodic-coefficient multi-blade analysis with trim and wake models', tool: 'Rotorcraft comprehensive analysis code' },
    { model: 'Fluid–structure–thermal interaction', why: 'Aerothermal heating coupled with structural response needs CFD heat-flux and thermal FE', tool: 'Coupled aerothermoelastic toolchain' },
  ],
};
