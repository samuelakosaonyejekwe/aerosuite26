// Canvas chart renderer for the suite plot specs: line, bar, heat (filled contours), tri (unstructured
// field), polar. Interactive by default: hover read-out, legend toggling, table view, CSV and PNG export.
// Colours come from CSS custom properties so light and dark themes are both first-class.

import { h, setKids, icon, num, toCsv, downloadText, downloadBlob } from './dom.js';

const css = (el, name) => getComputedStyle(el).getPropertyValue(name).trim();
const SEQ = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b'];
const hex = (c) => [parseInt(c.slice(1, 3), 16), parseInt(c.slice(3, 5), 16), parseInt(c.slice(5, 7), 16)];
const SEQ_RGB = SEQ.map(hex);
function ramp(stops, t) {
  t = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(t)), f = t - i, a = stops[i], b = stops[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}
function colourMap(dark, diverging, mid) {
  if (diverging) { const lo = hex('#2a78d6'), hi = hex('#e34948'), m = hex(mid); return (t) => (t < 0.5 ? ramp([lo, m], t * 2) : ramp([m, hi], (t - 0.5) * 2)); }
  const stops = dark ? SEQ_RGB.slice().reverse() : SEQ_RGB; // near-zero recedes towards the surface in both themes
  return (t) => ramp(stops, t);
}

export function niceTicks(lo, hi, n = 6) {
  if (!(hi > lo)) { const d = Math.abs(lo) || 1; lo -= d * 0.5; hi += d * 0.5; }
  const raw = (hi - lo) / n, mag = 10 ** Math.floor(Math.log10(raw)), r = raw / mag;
  const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag, t = [];
  for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) t.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return t;
}
function logTicks(lo, hi) { const t = []; for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) t.push(10 ** e); return t.filter((v) => v >= lo * 0.999 && v <= hi * 1.001); }
const tickLabel = (v, span) => { const a = Math.abs(v); if (v === 0) return '0'; if (a >= 1e5 || a < 1e-3) return v.toExponential(span && a / span > 50 ? 2 : 1).replace('e+', 'e'); const d = Math.max(0, Math.min(6, 2 - Math.floor(Math.log10(span || a)))); return String(Number(v.toFixed(d))); };
const finite = Number.isFinite;

function extent(arrs, positive = false) {
  let lo = Infinity, hi = -Infinity;
  for (const a of arrs) for (const v of a) if (finite(v) && (!positive || v > 0)) { if (v < lo) lo = v; if (v > hi) hi = v; }
  return lo === Infinity ? [0, 1] : [lo, hi];
}

/** Render a plot spec into `host`. Returns { destroy }. */
export function renderPlot(host, spec) {
  const hidden = new Set();
  let showTable = false, hover = null, geom = null;
  const canvas = h('canvas', { class: 'plot-canvas', role: 'img', 'aria-label': spec.title, tabIndex: 0 });
  const tip = h('div', { class: 'plot-tip', hidden: true });
  const legend = h('div', { class: 'plot-legend' });
  const tableHost = h('div', { class: 'plot-table', hidden: true });
  const wrap = h('figure', { class: `plot plot-${spec.type}` },
    h('figcaption', { class: 'plot-h' }, h('span', { class: 'plot-title' }, spec.title),
      h('span', { class: 'plot-tools' },
        h('button', { class: 'icon-btn', title: 'Show the data as a table', onclick: () => { showTable = !showTable; tableHost.hidden = !showTable; if (showTable) fillTable(); } }, icon('table', 16)),
        h('button', { class: 'icon-btn', title: 'Download data (CSV)', onclick: () => downloadText(slug(spec.title) + '.csv', toCsv(...tabular(spec)), 'text/csv') }, icon('download', 16)),
        h('button', { class: 'icon-btn', title: 'Download image (PNG)', onclick: () => canvas.toBlob((b) => b && downloadBlob(slug(spec.title) + '.png', b)) }, icon('image', 16)))),
    legend, h('div', { class: 'plot-area' }, canvas, tip), tableHost);
  host.append(wrap);

  const series = spec.series || [];
  const named = spec.type === 'line' || spec.type === 'bar' || spec.type === 'polar';
  if (named && series.length > 1) series.forEach((s, i) => {
    const b = h('button', { class: 'leg', type: 'button', title: 'Click to hide or show', onclick: () => { hidden.has(i) ? hidden.delete(i) : hidden.add(i); b.classList.toggle('off', hidden.has(i)); draw(); } },
      h('i', { class: `sw ${s.style === 'dash' ? 'dash' : ''}`, style: { '--c': `var(--series-${(i % 8) + 1})` } }), s.name);
    legend.append(b);
  });
  function fillTable() {
    const [cols, rows] = tabular(spec), max = 300;
    setKids(tableHost, h('table', { class: 'data' }, h('thead', null, h('tr', null, cols.map((c) => h('th', null, c)))), h('tbody', null, rows.slice(0, max).map((r) => h('tr', null, r.map((v) => h('td', null, typeof v === 'number' ? num(v, 5) : v)))))),
      rows.length > max ? h('p', { class: 'muted' }, `Showing ${max} of ${rows.length} rows. Download the CSV for all.`) : null);
  }

  function draw() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2.5), W = Math.max(240, wrap.clientWidth), small = W < 440;
    const H = spec.type === 'polar' ? Math.min(W, 380) : (spec.equalAspect && spec.type !== 'line') ? Math.min(Math.max(240, W * 0.62), 460) : Math.round(Math.min(360, Math.max(230, W * 0.56)));
    canvas.width = W * dpr; canvas.height = H * dpr; canvas.style.height = H + 'px';
    const g = canvas.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
    const dark = css(wrap, '--is-dark') === '1';
    const T = { ink: css(wrap, '--ink'), ink2: css(wrap, '--ink-2'), grid: css(wrap, '--grid'), surf: css(wrap, '--surface'), mid: dark ? '#383835' : '#f0efec', col: (i) => css(wrap, `--series-${(i % 8) + 1}`), dark, small };
    g.font = `${small ? 10.5 : 11.5}px ${css(wrap, '--font') || 'system-ui, sans-serif'}`; g.textBaseline = 'middle'; g.lineJoin = 'round'; g.lineCap = 'round';
    geom = ({ line: drawLine, bar: drawBar, heat: drawField, tri: drawField, polar: drawPolar }[spec.type] || drawLine)(g, W, H, T);
  }

  // ---- shared frame ---------------------------------------------------------------------
  function frame(g, W, H, T, xr, yr, { xlog, ylog, equal, right = 12, xcats } = {}) {
    const m = { l: 0, r: right, t: 10, b: 0 };
    const yt = ylog ? logTicks(yr[0], yr[1]) : niceTicks(yr[0], yr[1], Math.max(3, Math.round(H / 55)));
    if (!ylog && yt.length) { yr = [Math.min(yr[0], yt[0]), Math.max(yr[1], yt[yt.length - 1])]; }
    const ylab = yt.map((v) => tickLabel(v, yr[1] - yr[0]));
    m.l = Math.max(...ylab.map((s) => g.measureText(s).width), 14) + (spec.ylabel ? 30 : 12);
    m.b = (spec.xlabel ? 40 : 24) + (xcats && T.small ? 14 : 0);
    let pw = W - m.l - m.r, ph = H - m.t - m.b;
    if (equal) { // same scale on both axes
      const sx = pw / (xr[1] - xr[0]), sy = ph / (yr[1] - yr[0]), s = Math.min(sx, sy);
      const cx = (xr[0] + xr[1]) / 2, cy = (yr[0] + yr[1]) / 2; xr = [cx - pw / (2 * s), cx + pw / (2 * s)]; yr = [cy - ph / (2 * s), cy + ph / (2 * s)];
    }
    const xt = xcats ? [] : xlog ? logTicks(xr[0], xr[1]) : niceTicks(xr[0], xr[1], Math.max(3, Math.round(pw / 85)));
    const tx = xlog ? (v) => m.l + ((Math.log10(v) - Math.log10(xr[0])) / (Math.log10(xr[1]) - Math.log10(xr[0]))) * pw : (v) => m.l + ((v - xr[0]) / (xr[1] - xr[0])) * pw;
    const ty = ylog ? (v) => m.t + ph - ((Math.log10(v) - Math.log10(yr[0])) / (Math.log10(yr[1]) - Math.log10(yr[0]))) * ph : (v) => m.t + ph - ((v - yr[0]) / (yr[1] - yr[0])) * ph;
    const yts = equal ? niceTicks(yr[0], yr[1], Math.max(3, Math.round(ph / 55))) : yt;
    g.strokeStyle = T.grid; g.lineWidth = 1; g.fillStyle = T.ink2;
    g.textAlign = 'right';
    for (const v of yts) { const y = Math.round(ty(v)) + 0.5; if (y < m.t - 1 || y > m.t + ph + 1) continue; g.beginPath(); g.moveTo(m.l, y); g.lineTo(m.l + pw, y); g.stroke(); g.fillText(tickLabel(v, yr[1] - yr[0]), m.l - 7, y); }
    g.textAlign = 'center';
    for (const v of xt) { const x = Math.round(tx(v)) + 0.5; if (x < m.l - 1 || x > m.l + pw + 1) continue; g.beginPath(); g.moveTo(x, m.t); g.lineTo(x, m.t + ph); g.stroke(); g.fillText(tickLabel(v, xr[1] - xr[0]), x, m.t + ph + 12); }
    g.strokeStyle = T.ink2; g.globalAlpha = 0.55; g.strokeRect(m.l + 0.5, m.t + 0.5, pw, ph); g.globalAlpha = 1;
    g.fillStyle = T.ink2;
    if (spec.xlabel) g.fillText(spec.xlabel, m.l + pw / 2, H - 10);
    if (spec.ylabel) { g.save(); g.translate(12, m.t + ph / 2); g.rotate(-Math.PI / 2); g.fillText(spec.ylabel, 0, 0); g.restore(); }
    return { m, pw, ph, tx, ty, xr, yr };
  }
  const clip = (g, f) => { g.save(); g.beginPath(); g.rect(f.m.l, f.m.t, f.pw, f.ph); g.clip(); };

  // ---- line -------------------------------------------------------------------------------
  function drawLine(g, W, H, T) {
    const vis = series.map((s, i) => ({ s, i })).filter((o) => !hidden.has(o.i));
    let xr = extent(vis.map((o) => o.s.x), spec.xlog), yr = extent(vis.map((o) => o.s.y), spec.ylog);
    for (const a of spec.annotations || []) { if (finite(a.x) && a.x >= xr[0] - (xr[1] - xr[0]) && a.x <= xr[1] + (xr[1] - xr[0]) * 0.25) xr = [Math.min(xr[0], a.x), Math.max(xr[1], a.x)]; if (finite(a.y) && Math.abs(a.y) <= Math.max(Math.abs(yr[0]), Math.abs(yr[1])) * 1.5 + 1e-12) yr = [Math.min(yr[0], a.y), Math.max(yr[1], a.y)]; }
    const pad = (r, log) => (log ? [r[0] / 1.15, r[1] * 1.15] : r[0] === r[1] ? [r[0] - (Math.abs(r[0]) || 1) * 0.5, r[1] + (Math.abs(r[1]) || 1) * 0.5] : [r[0] - (r[1] - r[0]) * 0.04, r[1] + (r[1] - r[0]) * 0.04]);
    xr = pad(xr, spec.xlog); yr = pad(yr, spec.ylog);
    const f = frame(g, W, H, T, xr, yr, { xlog: spec.xlog, ylog: spec.ylog, equal: spec.equalAspect });
    clip(g, f);
    for (const a of spec.annotations || []) {
      g.strokeStyle = T.ink2; g.setLineDash([4, 4]); g.lineWidth = 1; g.beginPath(); g.fillStyle = T.ink2;
      if (finite(a.x)) { const x = f.tx(a.x); g.moveTo(x, f.m.t); g.lineTo(x, f.m.t + f.ph); g.stroke(); g.setLineDash([]); g.save(); g.translate(x - 5, f.m.t + 6); g.rotate(-Math.PI / 2); g.textAlign = 'right'; g.fillText(a.label || '', 0, 0); g.restore(); }
      else if (finite(a.y)) { const y = f.ty(a.y); g.moveTo(f.m.l, y); g.lineTo(f.m.l + f.pw, y); g.stroke(); g.setLineDash([]); g.textAlign = 'right'; g.fillText(a.label || '', f.m.l + f.pw - 6, y - 8); }
      g.setLineDash([]);
    }
    const pts = [];
    for (const { s, i } of vis) {
      const c = T.col(i), st = s.style || 'line';
      g.strokeStyle = c; g.fillStyle = c; g.lineWidth = 2; g.setLineDash(st === 'dash' ? [6, 5] : []);
      if (st !== 'points') {
        g.beginPath(); let pen = false, py = 0;
        for (let k = 0; k < s.x.length; k++) {
          const xv = s.x[k], yv = s.y[k];
          if (!finite(xv) || !finite(yv) || (spec.xlog && xv <= 0) || (spec.ylog && yv <= 0)) { pen = false; continue; }
          const x = f.tx(xv), y = f.ty(yv);
          if (!pen) g.moveTo(x, y); else if (st === 'step') { g.lineTo(x, py); g.lineTo(x, y); } else g.lineTo(x, y);
          pen = true; py = y;
        }
        g.stroke(); g.setLineDash([]);
      }
      const marks = st === 'points' || st === 'line+points' || s.x.length === 1;
      for (let k = 0; k < s.x.length; k++) {
        if (!finite(s.x[k]) || !finite(s.y[k])) continue;
        const x = f.tx(s.x[k]), y = f.ty(s.y[k]); pts.push({ x, y, i, k });
        if (marks) { g.beginPath(); g.arc(x, y, 4, 0, 7); g.fillStyle = c; g.fill(); g.strokeStyle = T.surf; g.lineWidth = 1.5; g.stroke(); }
      }
    }
    if (hover) { const p = hover; g.strokeStyle = T.ink2; g.lineWidth = 1; g.globalAlpha = 0.6; g.beginPath(); g.moveTo(p.x, f.m.t); g.lineTo(p.x, f.m.t + f.ph); g.moveTo(f.m.l, p.y); g.lineTo(f.m.l + f.pw, p.y); g.stroke(); g.globalAlpha = 1; g.beginPath(); g.arc(p.x, p.y, 5.5, 0, 7); g.fillStyle = T.col(p.i); g.fill(); g.strokeStyle = T.surf; g.lineWidth = 2; g.stroke(); }
    g.restore();
    return {
      pick(mx, my) {
        let best = null, bd = 40 * 40;
        for (const p of pts) { const d = (p.x - mx) ** 2 + (p.y - my) ** 2; if (d < bd) { bd = d; best = p; } }
        if (!best) return null;
        const s = series[best.i];
        return { ...best, rows: [[spec.xlabel || 'x', num(s.x[best.k], 5)], [spec.ylabel || 'y', num(s.y[best.k], 5)]], head: series.length > 1 ? s.name : null, sw: best.i };
      },
    };
  }

  // ---- bar --------------------------------------------------------------------------------
  function drawBar(g, W, H, T) {
    const cats = spec.categories || [], vis = series.map((s, i) => ({ s, i })).filter((o) => !hidden.has(o.i)), stacked = !!spec.stacked;
    let lo = 0, hi = 0;
    cats.forEach((_, k) => { if (stacked) { let p = 0, n = 0; for (const { s } of vis) { const v = s.y[k]; if (finite(v)) v >= 0 ? (p += v) : (n += v); } hi = Math.max(hi, p); lo = Math.min(lo, n); } else for (const { s } of vis) { const v = s.y[k]; if (finite(v)) { hi = Math.max(hi, v); lo = Math.min(lo, v); } } });
    if (hi === lo) hi = lo + 1;
    const f = frame(g, W, H, T, [0, cats.length], [lo < 0 ? lo * 1.06 : 0, hi > 0 ? hi * 1.06 : 0], { xcats: true });
    const bw = f.pw / cats.length, inner = Math.min(bw * 0.72, 90), marks = [];
    g.textAlign = 'center'; g.fillStyle = T.ink2;
    cats.forEach((c, k) => { const label = String(c), maxw = bw - 6; let t = label; while (t.length > 3 && g.measureText(t).width > maxw) t = t.slice(0, -2); g.fillText(t === label ? t : t + '…', f.m.l + bw * (k + 0.5), f.m.t + f.ph + 12); });
    const y0 = f.ty(0);
    cats.forEach((_, k) => {
      let accP = 0, accN = 0;
      vis.forEach(({ s, i }, j) => {
        const v = s.y[k]; if (!finite(v)) return;
        const w = stacked ? inner : Math.max(3, inner / vis.length - 2), x = f.m.l + bw * (k + 0.5) - inner / 2 + (stacked ? 0 : j * (inner / vis.length) + 1);
        const base = stacked ? (v >= 0 ? accP : accN) : 0, top = base + v;
        if (stacked) v >= 0 ? (accP = top) : (accN = top);
        const ya = f.ty(base), yb = f.ty(top), yy = Math.min(ya, yb), hh = Math.max(1, Math.abs(ya - yb) - (stacked ? 2 : 0));
        g.fillStyle = T.col(i); const r = Math.min(4, hh / 2, w / 2);
        g.beginPath(); if (g.roundRect) g.roundRect(x, yy, w, hh, stacked ? 1.5 : v >= 0 ? [r, r, 0, 0] : [0, 0, r, r]); else g.rect(x, yy, w, hh); g.fill();
        marks.push({ x, y: yy, w, h: hh, i, k, v });
      });
    });
    g.strokeStyle = T.ink2; g.lineWidth = 1; g.beginPath(); g.moveTo(f.m.l, Math.round(y0) + 0.5); g.lineTo(f.m.l + f.pw, Math.round(y0) + 0.5); g.stroke();
    if (hover) { g.strokeStyle = T.ink; g.lineWidth = 1.5; g.strokeRect(hover.bx - 1, hover.by - 1, hover.bw + 2, hover.bh + 2); }
    return { pick(mx, my) { const mk = marks.find((b) => mx >= b.x - 3 && mx <= b.x + b.w + 3 && my >= b.y - 4 && my <= b.y + b.h + 4); return mk ? { x: mk.x + mk.w / 2, y: mk.y, bx: mk.x, by: mk.y, bw: mk.w, bh: mk.h, head: String(cats[mk.k]), sw: mk.i, rows: [[series[mk.i].name || spec.ylabel || 'value', num(mk.v, 5)]] } : null; } };
  }

  // ---- heat (rectilinear) and tri (unstructured) scalar fields -----------------------------
  function drawField(g, W, H, T) {
    const isTri = spec.type === 'tri';
    const xs = isTri ? spec.nodes.map((p) => p[0]) : spec.x, ys = isTri ? spec.nodes.map((p) => p[1]) : spec.y;
    const vals = isTri ? spec.values : spec.z.flat();
    let [zlo, zhi] = extent([vals]);
    if (spec.diverging) { const a = Math.max(Math.abs(zlo), Math.abs(zhi)); zlo = -a; zhi = a; }
    if (zhi === zlo) zhi = zlo + 1;
    const cm = colourMap(T.dark, spec.diverging, T.mid), cbw = 54;
    const f = frame(g, W, H, T, extent([xs]), extent([ys]), { equal: spec.equalAspect, right: cbw + 14 });
    const col = (v) => { const c = cm((v - zlo) / (zhi - zlo)); return `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`; };
    clip(g, f);
    if (isTri) {
      for (const t of spec.tris) {
        const a = spec.nodes[t[0]], b = spec.nodes[t[1]], c = spec.nodes[t[2]], v = (spec.values[t[0]] + spec.values[t[1]] + spec.values[t[2]]) / 3;
        if (!finite(v)) continue;
        g.beginPath(); g.moveTo(f.tx(a[0]), f.ty(a[1])); g.lineTo(f.tx(b[0]), f.ty(b[1])); g.lineTo(f.tx(c[0]), f.ty(c[1])); g.closePath();
        g.fillStyle = col(v); g.fill(); g.strokeStyle = spec.edges ? T.surf : g.fillStyle; g.lineWidth = spec.edges ? 0.5 : 0.7; g.stroke();
      }
    } else {
      // sample the grid bilinearly into an off-screen image no larger than the plot area
      const nx = xs.length, ny = ys.length, iw = Math.max(2, Math.min(Math.round(f.pw), 420)), ih = Math.max(2, Math.min(Math.round(f.ph), 320));
      const x0 = Math.max(f.xr[0], Math.min(xs[0], xs[nx - 1])), x1 = Math.min(f.xr[1], Math.max(xs[0], xs[nx - 1])), y0 = Math.max(f.yr[0], Math.min(ys[0], ys[ny - 1])), y1 = Math.min(f.yr[1], Math.max(ys[0], ys[ny - 1]));
      const off = document.createElement('canvas'); off.width = iw; off.height = ih;
      const og = off.getContext('2d'), img = og.createImageData(iw, ih), find = (arr, v) => { let lo = 0, hi = arr.length - 1; const asc = arr[hi] >= arr[0]; while (hi - lo > 1) { const m = (lo + hi) >> 1; if ((arr[m] > v) === asc) hi = m; else lo = m; } return lo; };
      const ix = new Array(iw), fx = new Array(iw);
      for (let i = 0; i < iw; i++) { const xv = x0 + ((i + 0.5) / iw) * (x1 - x0), k = find(xs, xv); ix[i] = k; fx[i] = Math.max(0, Math.min(1, (xv - xs[k]) / (xs[k + 1] - xs[k] || 1))); }
      for (let j = 0; j < ih; j++) {
        const yv = y1 - ((j + 0.5) / ih) * (y1 - y0), k = find(ys, yv), fy = Math.max(0, Math.min(1, (yv - ys[k]) / (ys[k + 1] - ys[k] || 1))), r0 = spec.z[k], r1 = spec.z[k + 1] || r0;
        for (let i = 0; i < iw; i++) {
          const a = r0[ix[i]], b = r0[ix[i] + 1] ?? a, c = r1[ix[i]], d = r1[ix[i] + 1] ?? c, p = (j * iw + i) * 4;
          if (!finite(a) || !finite(b) || !finite(c) || !finite(d)) { img.data[p + 3] = 0; continue; }
          const v = (a * (1 - fx[i]) + b * fx[i]) * (1 - fy) + (c * (1 - fx[i]) + d * fx[i]) * fy, cc = cm((v - zlo) / (zhi - zlo));
          img.data[p] = cc[0]; img.data[p + 1] = cc[1]; img.data[p + 2] = cc[2]; img.data[p + 3] = 255;
        }
      }
      og.putImageData(img, 0, 0); g.imageSmoothingEnabled = true;
      g.drawImage(off, f.tx(x0), f.ty(y1), f.tx(x1) - f.tx(x0), f.ty(y0) - f.ty(y1));
      if (spec.contours) { // marching squares iso-lines
        g.strokeStyle = T.ink; g.globalAlpha = 0.35; g.lineWidth = 0.8; g.beginPath();
        const lv = niceTicks(zlo, zhi, spec.contours);
        for (const L of lv) for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
          const q = [spec.z[j][i], spec.z[j][i + 1], spec.z[j + 1][i + 1], spec.z[j + 1][i]];
          if (!q.every(finite)) continue;
          const P = [[xs[i], ys[j]], [xs[i + 1], ys[j]], [xs[i + 1], ys[j + 1]], [xs[i], ys[j + 1]]], cut = [];
          for (let e = 0; e < 4; e++) { const a = q[e], b = q[(e + 1) % 4]; if ((a < L) !== (b < L)) { const t = (L - a) / (b - a), A = P[e], B = P[(e + 1) % 4]; cut.push([A[0] + t * (B[0] - A[0]), A[1] + t * (B[1] - A[1])]); } }
          for (let c = 0; c + 1 < cut.length; c += 2) { g.moveTo(f.tx(cut[c][0]), f.ty(cut[c][1])); g.lineTo(f.tx(cut[c + 1][0]), f.ty(cut[c + 1][1])); }
        }
        g.stroke(); g.globalAlpha = 1;
      }
    }
    for (const o of spec.overlay || []) { g.strokeStyle = T.ink; g.lineWidth = 1.6; g.beginPath(); let pen = false; for (let k = 0; k < o.x.length; k++) { if (!finite(o.x[k]) || !finite(o.y[k])) { pen = false; continue; } pen ? g.lineTo(f.tx(o.x[k]), f.ty(o.y[k])) : g.moveTo(f.tx(o.x[k]), f.ty(o.y[k])); pen = true; } if (o.fill) { g.fillStyle = T.surf; g.fill(); } g.stroke(); }
    if (hover) { g.strokeStyle = T.ink; g.lineWidth = 1.5; g.beginPath(); g.arc(hover.x, hover.y, 5, 0, 7); g.stroke(); }
    g.restore();
    // colour bar
    const bx = W - cbw, bh = f.ph, by = f.m.t;
    for (let k = 0; k < bh; k++) { g.fillStyle = col(zhi - (k / bh) * (zhi - zlo)); g.fillRect(bx, by + k, 12, 1.5); }
    g.strokeStyle = T.ink2; g.globalAlpha = 0.5; g.strokeRect(bx + 0.5, by + 0.5, 12, bh); g.globalAlpha = 1;
    g.fillStyle = T.ink2; g.textAlign = 'left';
    for (const v of niceTicks(zlo, zhi, 5)) { if (v < zlo || v > zhi) continue; g.fillText(tickLabel(v, zhi - zlo), bx + 16, by + ((zhi - v) / (zhi - zlo)) * bh); }
    if (spec.zlabel) { g.save(); g.translate(W - 5, by + bh / 2); g.rotate(-Math.PI / 2); g.textAlign = 'center'; g.fillText(spec.zlabel, 0, 0); g.restore(); }
    const ux = (px) => f.xr[0] + ((px - f.m.l) / f.pw) * (f.xr[1] - f.xr[0]), uy = (py) => f.yr[0] + ((f.m.t + f.ph - py) / f.ph) * (f.yr[1] - f.yr[0]);
    return {
      pick(mx, my) {
        if (mx < f.m.l || mx > f.m.l + f.pw || my < f.m.t || my > f.m.t + f.ph) return null;
        const X = ux(mx), Y = uy(my); let v, px = mx, py = my;
        if (isTri) { let bd = Infinity, bi = -1; for (let k = 0; k < spec.nodes.length; k++) { const d = (f.tx(spec.nodes[k][0]) - mx) ** 2 + (f.ty(spec.nodes[k][1]) - my) ** 2; if (d < bd) { bd = d; bi = k; } } if (bd > 900) return null; v = spec.values[bi]; px = f.tx(spec.nodes[bi][0]); py = f.ty(spec.nodes[bi][1]); }
        else { const near = (arr, q) => { let b = 0; for (let k = 1; k < arr.length; k++) if (Math.abs(arr[k] - q) < Math.abs(arr[b] - q)) b = k; return b; }; const i = near(xs, X), j = near(ys, Y); v = spec.z[j][i]; if (!finite(v)) return null; }
        return { x: px, y: py, rows: [[spec.zlabel || 'value', num(v, 5)], [spec.xlabel || 'x', num(X, 4)], [spec.ylabel || 'y', num(Y, 4)]] };
      },
    };
  }

  // ---- polar --------------------------------------------------------------------------------
  function drawPolar(g, W, H, T) {
    const vis = series.map((s, i) => ({ s, i })).filter((o) => !hidden.has(o.i));
    let [lo, hi] = extent(vis.map((o) => o.s.r)); const ticks = niceTicks(lo, hi, 4); lo = Math.min(lo, ticks[0]); hi = Math.max(hi, ticks[ticks.length - 1]); if (lo > 0 && lo < 0.3 * hi) lo = 0;
    const cx = W / 2, cy = H / 2, R = Math.min(W, H) / 2 - 26, rr = (v) => ((v - lo) / (hi - lo || 1)) * R, pts = [];
    g.strokeStyle = T.grid; g.fillStyle = T.ink2; g.lineWidth = 1;
    for (const t of ticks) { if (t < lo) continue; g.beginPath(); g.arc(cx, cy, rr(t), 0, 7); g.stroke(); g.textAlign = 'left'; g.fillText(tickLabel(t, hi - lo), cx + 4, cy - rr(t) - 7); }
    for (let a = 0; a < 360; a += 30) { const th = (a * Math.PI) / 180; g.beginPath(); g.moveTo(cx, cy); g.lineTo(cx + R * Math.cos(th), cy - R * Math.sin(th)); g.stroke(); g.textAlign = 'center'; g.fillText(a + '°', cx + (R + 13) * Math.cos(th), cy - (R + 13) * Math.sin(th)); }
    for (const { s, i } of vis) {
      g.strokeStyle = T.col(i); g.lineWidth = 2; g.beginPath(); let pen = false;
      for (let k = 0; k < s.r.length; k++) { if (!finite(s.r[k])) { pen = false; continue; } const th = (s.theta_deg[k] * Math.PI) / 180, x = cx + rr(s.r[k]) * Math.cos(th), y = cy - rr(s.r[k]) * Math.sin(th); pen ? g.lineTo(x, y) : g.moveTo(x, y); pen = true; pts.push({ x, y, i, k }); }
      g.stroke();
    }
    if (spec.rlabel) { g.fillStyle = T.ink2; g.textAlign = 'left'; g.fillText(spec.rlabel, 6, H - 9); }
    if (hover) { g.beginPath(); g.arc(hover.x, hover.y, 5.5, 0, 7); g.fillStyle = T.col(hover.i); g.fill(); g.strokeStyle = T.surf; g.lineWidth = 2; g.stroke(); }
    return { pick(mx, my) { let best = null, bd = 1600; for (const p of pts) { const d = (p.x - mx) ** 2 + (p.y - my) ** 2; if (d < bd) { bd = d; best = p; } } if (!best) return null; const s = series[best.i]; return { ...best, head: series.length > 1 ? s.name : null, sw: best.i, rows: [['Angle', num(s.theta_deg[best.k], 4) + '°'], [spec.rlabel || 'r', num(s.r[best.k], 5)]] }; } };
  }

  // ---- interaction --------------------------------------------------------------------------
  function onMove(ev) {
    if (!geom) return;
    const r = canvas.getBoundingClientRect(), p = geom.pick(ev.clientX - r.left, ev.clientY - r.top);
    const changed = (p && (!hover || p.x !== hover.x || p.y !== hover.y)) || (!p && hover);
    hover = p;
    if (changed) draw();
    if (!p) { tip.hidden = true; return; }
    setKids(tip, p.head ? h('div', { class: 'tip-h' }, p.sw != null ? h('i', { class: 'sw', style: { '--c': `var(--series-${(p.sw % 8) + 1})` } }) : null, p.head) : null, ...p.rows.map(([k, v]) => h('div', { class: 'tip-r' }, h('span', null, k), h('b', null, v))));
    tip.hidden = false;
    const tw = tip.offsetWidth, th = tip.offsetHeight, x = p.x + 14 + tw > r.width ? p.x - tw - 14 : p.x + 14, y = Math.max(4, Math.min(r.height - th - 4, p.y - th / 2));
    tip.style.transform = `translate(${Math.max(2, x)}px, ${y}px)`;
  }
  const onLeave = () => { if (hover) { hover = null; draw(); } tip.hidden = true; };
  canvas.addEventListener('pointermove', onMove); canvas.addEventListener('pointerdown', onMove); canvas.addEventListener('pointerleave', onLeave);
  let raf = 0; const ro = new ResizeObserver(() => { cancelAnimationFrame(raf); raf = requestAnimationFrame(draw); }); ro.observe(wrap);
  const onTheme = () => draw(); window.addEventListener('themechange', onTheme);
  draw();
  return { destroy() { ro.disconnect(); window.removeEventListener('themechange', onTheme); wrap.remove(); } };
}

const slug = (s) => String(s).toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '').slice(0, 60) || 'chart';

/** Flatten any plot spec to [columns, rows] for the table view and CSV export. */
export function tabular(spec) {
  if (spec.type === 'bar') return [['Category', ...spec.series.map((s) => s.name)], spec.categories.map((c, k) => [c, ...spec.series.map((s) => s.y[k])])];
  if (spec.type === 'heat') return [[`${spec.ylabel || 'y'} \\ ${spec.xlabel || 'x'}`, ...spec.x], spec.y.map((y, j) => [y, ...spec.z[j]])];
  if (spec.type === 'tri') return [[spec.xlabel || 'x', spec.ylabel || 'y', spec.zlabel || 'value'], spec.nodes.map((p, k) => [p[0], p[1], spec.values[k]])];
  if (spec.type === 'polar') { const rows = []; for (const s of spec.series) s.r.forEach((r, k) => rows.push([s.name, s.theta_deg[k], r])); return [['Series', 'Angle [deg]', spec.rlabel || 'r'], rows]; }
  const rows = []; for (const s of spec.series) s.x.forEach((x, k) => rows.push([s.name, x, s.y[k]]));
  return [['Series', spec.xlabel || 'x', spec.ylabel || 'y'], rows];
}
