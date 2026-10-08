// Lightweight WebGL viewer for imported geometry and meshes: shaded triangles, element edges, 1-D
// elements and point clouds, with orbit / pan / zoom on mouse and touch. No third-party code.

import { h, icon } from './dom.js';

const VS = `attribute vec3 p; attribute vec3 n; uniform mat4 mvp; uniform mat3 nm; uniform float ps; varying vec3 vn;
void main(){ vn = nm * n; gl_Position = mvp * vec4(p, 1.0); gl_PointSize = ps; }`;
const FS = `precision mediump float; varying vec3 vn; uniform vec3 col; uniform float lit;
void main(){ vec3 N = normalize(vn); float d = abs(dot(N, normalize(vec3(0.35, 0.55, 0.75)))); float k = mix(1.0, 0.28 + 0.72 * d, lit); gl_FragColor = vec4(col * k, 1.0); }`;
const MAX_TRIS = 1500000, MAX_PTS = 2000000;

const mul = (a, b) => { const o = new Float32Array(16); for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) { let s = 0; for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k]; o[c * 4 + r] = s; } return o; };
const persp = (fov, asp, n, f) => { const t = 1 / Math.tan(fov / 2), o = new Float32Array(16); o[0] = t / asp; o[5] = t; o[10] = (f + n) / (n - f); o[11] = -1; o[14] = (2 * f * n) / (n - f); return o; };

/** Create a viewer inside `host`. Returns { setModel(model), fit(), destroy() }. */
export function createViewer(host) {
  const canvas = h('canvas', { 'aria-label': '3-D view of the imported geometry. Drag to rotate, scroll or pinch to zoom, shift-drag or two-finger drag to pan.', tabIndex: 0 });
  const hud = h('div', { class: 'hud' }, 'Drag to rotate · scroll or pinch to zoom · shift-drag to pan');
  let edges = false;
  const edgeBtn = h('button', { class: 'icon-btn', title: 'Show or hide element edges', onclick: () => { edges = !edges; draw(); } }, icon('layers', 18));
  const wrap = h('div', { class: 'viewer' }, canvas, h('div', { class: 'tools' }, edgeBtn, h('button', { class: 'icon-btn', title: 'Fit to view', onclick: () => fit() }, icon('refresh', 18)), ...[['X', 0, 0], ['Y', Math.PI / 2, 0], ['Z', 0, Math.PI / 2 - 0.001]].map(([l, y, p]) => h('button', { class: 'icon-btn', title: `Look along ${l}`, onclick: () => { yaw = y; pitch = p; draw(); } }, l))), hud);
  host.append(wrap);
  const gl = canvas.getContext('webgl', { antialias: true, alpha: true, preserveDrawingBuffer: true });
  if (!gl) { hud.textContent = 'This device does not provide WebGL; the 3-D preview is unavailable, but all measurements still work.'; return { setModel() {}, fit() {}, destroy() { wrap.remove(); } }; }
  gl.getExtension('OES_element_index_uint');
  const sh = (t, s) => { const o = gl.createShader(t); gl.shaderSource(o, s); gl.compileShader(o); return o; };
  const prog = gl.createProgram(); gl.attachShader(prog, sh(gl.VERTEX_SHADER, VS)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, FS)); gl.linkProgram(prog); gl.useProgram(prog);
  const loc = { p: gl.getAttribLocation(prog, 'p'), n: gl.getAttribLocation(prog, 'n'), mvp: gl.getUniformLocation(prog, 'mvp'), nm: gl.getUniformLocation(prog, 'nm'), col: gl.getUniformLocation(prog, 'col'), lit: gl.getUniformLocation(prog, 'lit'), ps: gl.getUniformLocation(prog, 'ps') };
  const buf = { tri: gl.createBuffer(), triN: gl.createBuffer(), edge: gl.createBuffer(), line: gl.createBuffer(), pts: gl.createBuffer() };
  let n = { tri: 0, edge: 0, line: 0, pts: 0 }, centre = [0, 0, 0], radius = 1, yaw = 0.7, pitch = 0.45, dist = 3, pan = [0, 0];

  function setModel(m) {
    const P = m.positions, T = m.triangles || new Uint32Array(0), nv = P.length / 3;
    let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < P.length; i += 3) for (let k = 0; k < 3; k++) { const v = P[i + k]; if (v < lo[k]) lo[k] = v; if (v > hi[k]) hi[k] = v; }
    if (!nv) { lo = [0, 0, 0]; hi = [1, 1, 1]; }
    centre = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2]; radius = Math.max(1e-9, Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2);
    const s = 1 / radius, X = (i, k) => (P[i * 3 + k] - centre[k]) * s; // normalise to a unit sphere for float precision
    const nt = T.length / 3, step = Math.max(1, Math.ceil(nt / MAX_TRIS)), kept = Math.ceil(nt / step);
    const tp = new Float32Array(kept * 9), tn = new Float32Array(kept * 9), ep = new Float32Array(kept * 18);
    let o = 0;
    for (let t = 0; t < nt; t += step) {
      const a = T[t * 3], b = T[t * 3 + 1], c = T[t * 3 + 2]; if (a >= nv || b >= nv || c >= nv) continue;
      const ax = X(a, 0), ay = X(a, 1), az = X(a, 2), bx = X(b, 0), by = X(b, 1), bz = X(b, 2), cx = X(c, 0), cy = X(c, 1), cz = X(c, 2);
      const ux = bx - ax, uy = by - ay, uz = bz - az, vx = cx - ax, vy = cy - ay, vz = cz - az; let nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx; const l = Math.hypot(nx, ny, nz) || 1; nx /= l; ny /= l; nz /= l;
      tp.set([ax, ay, az, bx, by, bz, cx, cy, cz], o * 9); tn.set([nx, ny, nz, nx, ny, nz, nx, ny, nz], o * 9);
      ep.set([ax, ay, az, bx, by, bz, bx, by, bz, cx, cy, cz, cx, cy, cz, ax, ay, az], o * 18); o++;
    }
    n.tri = o * 3; n.edge = o * 6;
    gl.bindBuffer(gl.ARRAY_BUFFER, buf.tri); gl.bufferData(gl.ARRAY_BUFFER, tp.subarray(0, o * 9), gl.STATIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, buf.triN); gl.bufferData(gl.ARRAY_BUFFER, tn.subarray(0, o * 9), gl.STATIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, buf.edge); gl.bufferData(gl.ARRAY_BUFFER, ep.subarray(0, o * 18), gl.STATIC_DRAW);
    const L = m.lines || new Uint32Array(0), lp = new Float32Array(L.length * 3); let lc = 0;
    for (let i = 0; i < L.length; i++) { if (L[i] >= nv) continue; lp[lc * 3] = X(L[i], 0); lp[lc * 3 + 1] = X(L[i], 1); lp[lc * 3 + 2] = X(L[i], 2); lc++; }
    n.line = lc - (lc % 2); gl.bindBuffer(gl.ARRAY_BUFFER, buf.line); gl.bufferData(gl.ARRAY_BUFFER, lp.subarray(0, n.line * 3), gl.STATIC_DRAW);
    // show vertices as points when there is no surface (point clouds, CAD control points)
    n.pts = 0;
    if (!n.tri) { const ps = Math.max(1, Math.ceil(nv / MAX_PTS)), pp = new Float32Array(Math.ceil(nv / ps) * 3); let pc = 0; for (let i = 0; i < nv; i += ps) { pp[pc * 3] = X(i, 0); pp[pc * 3 + 1] = X(i, 1); pp[pc * 3 + 2] = X(i, 2); pc++; } n.pts = pc; gl.bindBuffer(gl.ARRAY_BUFFER, buf.pts); gl.bufferData(gl.ARRAY_BUFFER, pp, gl.STATIC_DRAW); }
    hud.textContent = `${nt.toLocaleString()} triangles${step > 1 ? ` (showing every ${step}${step === 2 ? 'nd' : 'th'})` : ''}${n.line ? ` · ${(n.line / 2).toLocaleString()} line elements` : ''}${n.pts ? ` · ${nv.toLocaleString()} points` : ''} · drag to rotate, scroll or pinch to zoom`;
    edges = nt > 0 && nt < 60000 && (m.kind || '').includes('mesh');
    fit();
  }
  function fit() { dist = 2.6; pan = [0, 0]; draw(); }

  function draw() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2), w = Math.max(1, Math.round(canvas.clientWidth * dpr)), hh = Math.max(1, Math.round(canvas.clientHeight * dpr));
    if (canvas.width !== w || canvas.height !== hh) { canvas.width = w; canvas.height = hh; }
    gl.viewport(0, 0, w, hh); gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT); gl.enable(gl.DEPTH_TEST);
    const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
    // model axes: x aft, y span, z up  ->  view: rotate about z (yaw) then tilt (pitch)
    const R = new Float32Array([cy, -sy * sp, sy * cp, 0, sy, cy * sp, -cy * cp, 0, 0, cp, sp, 0, 0, 0, 0, 1]);
    const V = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, pan[0], pan[1], -dist, 1]);
    const mvp = mul(persp(0.7, w / hh, Math.max(0.01, dist - 1.5), dist + 1.5), mul(V, R));
    gl.uniformMatrix4fv(loc.mvp, false, mvp); gl.uniformMatrix3fv(loc.nm, false, new Float32Array([R[0], R[1], R[2], R[4], R[5], R[6], R[8], R[9], R[10]]));
    const dark = document.documentElement.dataset.theme === 'dark', bind = (b, a) => { gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.enableVertexAttribArray(a); gl.vertexAttribPointer(a, 3, gl.FLOAT, false, 0, 0); };
    gl.uniform1f(loc.ps, 2 * dpr);
    if (n.tri) { bind(buf.tri, loc.p); bind(buf.triN, loc.n); gl.uniform1f(loc.lit, 1); gl.uniform3fv(loc.col, dark ? [0.42, 0.62, 0.92] : [0.33, 0.56, 0.88]); gl.enable(gl.POLYGON_OFFSET_FILL); gl.polygonOffset(1, 1); gl.drawArrays(gl.TRIANGLES, 0, n.tri); gl.disable(gl.POLYGON_OFFSET_FILL); }
    gl.disableVertexAttribArray(loc.n); gl.vertexAttrib3f(loc.n, 0, 0, 1); gl.uniform1f(loc.lit, 0);
    if (edges && n.edge) { bind(buf.edge, loc.p); gl.uniform3fv(loc.col, dark ? [0.07, 0.1, 0.16] : [0.1, 0.2, 0.38]); gl.drawArrays(gl.LINES, 0, n.edge); }
    if (n.line) { bind(buf.line, loc.p); gl.uniform3fv(loc.col, [0.92, 0.41, 0.2]); gl.drawArrays(gl.LINES, 0, n.line); }
    if (n.pts) { bind(buf.pts, loc.p); gl.uniform3fv(loc.col, dark ? [0.5, 0.72, 0.98] : [0.16, 0.44, 0.82]); gl.drawArrays(gl.POINTS, 0, n.pts); }
  }

  // ---- interaction ----
  const ptrs = new Map(); let lastPinch = 0;
  canvas.addEventListener('pointerdown', (e) => { canvas.setPointerCapture(e.pointerId); ptrs.set(e.pointerId, [e.clientX, e.clientY]); canvas.style.cursor = 'grabbing'; });
  canvas.addEventListener('pointermove', (e) => {
    const p = ptrs.get(e.pointerId); if (!p) return;
    const dx = e.clientX - p[0], dy = e.clientY - p[1]; ptrs.set(e.pointerId, [e.clientX, e.clientY]);
    if (ptrs.size === 2) { const [a, b] = [...ptrs.values()], d = Math.hypot(a[0] - b[0], a[1] - b[1]); if (lastPinch) dist = Math.max(0.3, Math.min(20, dist * (lastPinch / d))); lastPinch = d; pan[0] += (dx * dist) / (2 * canvas.clientHeight); pan[1] -= (dy * dist) / (2 * canvas.clientHeight); }
    else if (e.shiftKey || e.buttons === 4 || e.buttons === 2) { pan[0] += (dx * dist) / canvas.clientHeight; pan[1] -= (dy * dist) / canvas.clientHeight; }
    else { yaw += dx * 0.008; pitch = Math.max(-1.57, Math.min(1.57, pitch + dy * 0.008)); }
    draw();
  });
  const up = (e) => { ptrs.delete(e.pointerId); lastPinch = 0; canvas.style.cursor = 'grab'; };
  canvas.addEventListener('pointerup', up); canvas.addEventListener('pointercancel', up);
  canvas.addEventListener('wheel', (e) => { e.preventDefault(); dist = Math.max(0.3, Math.min(20, dist * Math.exp(e.deltaY * 0.0012))); draw(); }, { passive: false });
  canvas.addEventListener('contextmenu', (e) => e.preventDefault());
  canvas.addEventListener('keydown', (e) => { const k = { ArrowLeft: [-0.1, 0], ArrowRight: [0.1, 0], ArrowUp: [0, -0.1], ArrowDown: [0, 0.1] }[e.key]; if (k) { yaw += k[0]; pitch = Math.max(-1.57, Math.min(1.57, pitch + k[1])); draw(); e.preventDefault(); } if (e.key === '+' || e.key === '=') { dist *= 0.9; draw(); } if (e.key === '-') { dist *= 1.1; draw(); } });
  const ro = new ResizeObserver(() => draw()); ro.observe(canvas);
  const onTheme = () => draw(); window.addEventListener('themechange', onTheme);
  return { setModel, fit, destroy() { ro.disconnect(); window.removeEventListener('themechange', onTheme); gl.getExtension('WEBGL_lose_context')?.loseContext(); wrap.remove(); } };
}
