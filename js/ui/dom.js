// Tiny DOM toolkit: element builder, number formatting, icons, toasts and file helpers.
// All text goes through textContent, never innerHTML, so imported data can never inject markup.

export function h(tag, attrs, ...kids) {
  const el = tag === 'svg' || tag === 'path' || tag === 'circle' || tag === 'g' || tag === 'line' || tag === 'rect' || tag === 'polyline' || tag === 'text'
    ? document.createElementNS('http://www.w3.org/2000/svg', tag) : document.createElement(tag);
  if (attrs) for (const [k, v] of Object.entries(attrs)) {
    if (v === false || v == null) continue;
    if (k === 'class') el.setAttribute('class', v);
    else if (k === 'style' && typeof v === 'object') { for (const [p, q] of Object.entries(v)) p.startsWith('--') ? el.style.setProperty(p, q) : (el.style[p] = q); }
    else if (k.startsWith('on') && typeof v === 'function') el.addEventListener(k.slice(2), v);
    else if (k === 'dataset') Object.assign(el.dataset, v);
    else if (k in el && !(el instanceof SVGElement) && k !== 'list' && k !== 'form') { try { el[k] = v; } catch { el.setAttribute(k, v); } }
    else el.setAttribute(k, v === true ? '' : v);
  }
  add(el, kids);
  return el;
}
export function add(el, kids) {
  for (const k of kids.flat(Infinity)) {
    if (k == null || k === false) continue;
    el.append(k instanceof Node ? k : document.createTextNode(String(k)));
  }
  return el;
}
/** Replace an element's children (null/false entries are skipped, unlike Element.replaceChildren). */
export const setKids = (el, ...kids) => add(clear(el), kids);
export const clear = (el) => { while (el.firstChild) el.firstChild.remove(); return el; };
export const $ = (sel, root = document) => root.querySelector(sel);

/** Human-friendly number: 4 significant figures, thin-space thousands, scientific outside 1e-4..1e9. */
export function num(v, sig = 4) {
  if (typeof v !== 'number') return v == null ? '–' : String(v);
  if (Number.isNaN(v)) return 'n/a';
  if (!Number.isFinite(v)) return v > 0 ? '∞' : '−∞';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e9 || a < 1e-4) { const [m, e] = v.toExponential(sig - 1).split('e'); return `${trimZeros(m)}×10${sup(Number(e))}`.replace('-', '−'); }
  const digits = Math.max(0, sig - 1 - Math.floor(Math.log10(a)));
  let s = v.toFixed(Math.min(digits, 8));
  if (s.includes('.')) s = trimZeros(s);
  const [i, f] = s.split('.');
  const ii = i.replace(/\B(?=(\d{3})+(?!\d))/g, ' ');
  return (f ? `${ii}.${f}` : ii).replace('-', '−');
}
const trimZeros = (s) => (s.includes('.') ? s.replace(/0+$/, '').replace(/\.$/, '') : s);
const sup = (n) => String(n).split('').map((c) => ({ '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' }[c] ?? c)).join('');
export const withUnit = (v, unit) => (unit && unit !== '-' && unit !== '' ? `${num(v)} ${unit}` : num(v));
export const ago = (ts) => {
  if (!ts) return 'never';
  const s = (Date.now() - ts) / 1000;
  return s < 60 ? 'just now' : s < 3600 ? `${Math.floor(s / 60)} min ago` : s < 86400 ? `${Math.floor(s / 3600)} h ago` : `${Math.floor(s / 86400)} d ago`;
};

// Stroke icons on a 24-unit grid.
const ICONS = {
  back: 'M15 5l-7 7 7 7', fwd: 'M9 5l7 7-7 7', up: 'M5 15l7-7 7 7', down: 'M5 9l7 7 7-7', menu: 'M4 6h16M4 12h16M4 18h16', close: 'M6 6l12 12M18 6L6 18',
  home: 'M4 11l8-7 8 7v9h-5v-6H9v6H4z', play: 'M7 5l12 7-12 7z', stop: 'M7 7h10v10H7z', refresh: 'M20 12a8 8 0 1 1-2.3-5.7M20 4v5h-5', download: 'M12 4v11m-5-5l5 5 5-5M5 20h14',
  upload: 'M12 20V9m-5 5l5-5 5 5M5 4h14', search: 'M11 4a7 7 0 1 1 0 14 7 7 0 0 1 0-14zm5 12l4 4', sun: 'M12 8a4 4 0 1 1 0 8 4 4 0 0 1 0-8zM12 2v3M12 19v3M2 12h3M19 12h3M5 5l2 2M17 17l2 2M19 5l-2 2M7 17l-2 2',
  moon: 'M20 14A8 8 0 1 1 10 4a6.5 6.5 0 0 0 10 10z', info: 'M12 3a9 9 0 1 1 0 18 9 9 0 0 1 0-18zm0 8v6m0-9v.5', warn: 'M12 4l9 16H3zM12 10v5m0 2v.5', check: 'M5 12l5 5 9-10',
  link: 'M10 14a4 4 0 0 0 6 0l3-3a4 4 0 0 0-6-6l-1 1M14 10a4 4 0 0 0-6 0l-3 3a4 4 0 0 0 6 6l1-1', table: 'M4 5h16v14H4zM4 10h16M4 15h16M10 5v14', image: 'M4 5h16v14H4zM4 16l5-5 4 4 3-3 4 4M9 9v.5',
  globe: 'M12 3a9 9 0 1 1 0 18 9 9 0 0 1 0-18zM3 12h18M12 3c3 3 3 15 0 18M12 3c-3 3-3 15 0 18', cube: 'M12 3l8 4.5v9L12 21l-8-4.5v-9zM4 7.5l8 4.5 8-4.5M12 12v9', sliders: 'M4 7h10M18 7h2M4 17h4M12 17h8M14 5v4M8 15v4',
  layers: 'M12 4l9 5-9 5-9-5zM3 14l9 5 9-5', bulb: 'M9 18h6M10 21h4M12 3a6 6 0 0 0-4 10.5c.7.7 1 1.5 1 2.5h6c0-1 .3-1.8 1-2.5A6 6 0 0 0 12 3z', doc: 'M7 3h7l4 4v14H7zM14 3v4h4M10 12h5M10 16h5',
  graph: 'M6 6a2 2 0 1 1 0 4 2 2 0 0 1 0-4zM18 4a2 2 0 1 1 0 4 2 2 0 0 1 0-4zM17 16a2 2 0 1 1 0 4 2 2 0 0 1 0-4zM8 8.5l8-2M7.5 9.8l8 6.6', install: 'M12 3v10m-4-4l4 4 4-4M5 15v4h14v-4', pin: 'M12 21s7-6.2 7-11a7 7 0 1 0-14 0c0 4.8 7 11 7 11zM12 8a2 2 0 1 1 0 4 2 2 0 0 1 0-4z',
  flow: 'M3 8c4-4 8 4 12 0s4-2 6-1M3 14c4-4 8 4 12 0s4-2 6-1M3 19c5-2 9 1 18-2', truss: 'M3 18h18M3 18l4.5-9 4.5 9 4.5-9 4.5 9M7.5 9h9', flutter: 'M3 12c2-6 4-6 6 0s4 6 6 0 4-6 6 0', axes: 'M5 19V5M5 19h14M5 19l8-8M5 5l-2 2M5 5l2 2M19 19l-2-2M19 19l-2 2',
  envelope: 'M4 19V8c3-5 9-5 13 0l3 5v6z', rotor: 'M12 12L4 8M12 12l8-4M12 12v9M12 10.5a1.5 1.5 0 1 1 0 3 1.5 1.5 0 0 1 0-3z', engine: 'M4 9h4l2-2h5l3 2h2v6h-2l-3 2h-5l-2-2H4zM10 9v6M14 9v6', prop: 'M12 12c-2-6 2-9 0-9s-2 3 0 9zm0 0c6 2 9-2 9 0s-3 2-9 0zm0 0c2 6-2 9 0 9s2-3 0-9zm0 0c-6-2-9 2-9 0s3-2 9 0z',
  crack: 'M4 4h16v16H4zM12 4l-2 5 3 3-2 4 1 4', wave: 'M3 12h3l2-7 4 14 3-10 2 3h4', sound: 'M4 10v4h3l5 4V6L7 10zM16 9c1.5 1.5 1.5 4.5 0 6M18.5 6.5c3 3 3 8 0 11', heat: 'M12 3c1 4 5 5 5 10a5 5 0 0 1-10 0c0-2 1-3 2-4 0 2 1 3 2 3 0-3-1-5 1-9z',
  ice: 'M12 3v18M4.2 7.5l15.6 9M19.8 7.5l-15.6 9M9 4.5l3 2 3-2M9 19.5l3-2 3 2', gear: 'M12 3v9M8 12h8l-1 4H9zM7.5 19.5a2 2 0 1 1 4 0 2 2 0 0 1-4 0zM12.5 19.5a2 2 0 1 1 4 0 2 2 0 0 1-4 0z', impact: 'M3 20h18M12 4v10m-4-4l4 4 4-4M5 17l-2-2M19 17l2-2',
  loop: 'M5 8h10a4 4 0 0 1 0 8H7M5 8l3-3M5 8l3 3M7 16l3-3M7 16l3 3', radar: 'M12 12L19 5M12 3a9 9 0 1 0 9 9M12 7a5 5 0 1 0 5 5M12 11.2a.8.8 0 1 1 0 1.6.8.8 0 0 1 0-1.6z', piston: 'M4 9h10v6H4zM14 12h6M8 9v6M20 9v6',
  bolt: 'M13 3L5 13h6l-1 8 8-10h-6z', drop: 'M12 3c4 5 6 8 6 11a6 6 0 0 1-12 0c0-3 2-6 6-11z', shield: 'M12 3l8 3v6c0 5-3.5 8-8 9-4.5-1-8-4-8-9V6zM8.5 12l2.5 2.5 4.5-5', pareto: 'M4 4v16h16M7 8c1 5 5 8 11 9M7 8v.5M11 13v.5M18 17v.5',
  route: 'M6 19a2 2 0 1 1 0-4 2 2 0 0 1 0 4zM18 9a2 2 0 1 1 0-4 2 2 0 0 1 0 4zM8 17h6a3 3 0 0 0 0-6h-4a3 3 0 0 1 0-6h6', coin: 'M12 4a8 8 0 1 1 0 16 8 8 0 0 1 0-16zM14.5 9.5c-.5-1-1.5-1.5-2.5-1.5-1.5 0-2.5.8-2.5 2s1 1.7 2.5 2 2.5.8 2.5 2-1 2-2.5 2c-1.2 0-2.2-.6-2.7-1.6M12 6.5V8m0 8v1.5',
  jet: 'M3 13l7-1 5-7h2l-2 7 5-.5 1.5-2H23l-1 3 1 3h-1.5L20 13.5l-5-.5 2 7h-2l-5-7-7-1z', heli: 'M3 6h18M12 6v3M6 13c0-2 2-4 6-4s6 2 6 4-2 3-6 3H4l-1-4zM9 16v3h6v-3M6 19h12', ga: 'M3 12h18M12 5v14M8 19h8M10 5h4', uav: 'M12 5v14M3 10h18M9 19h6', quad: 'M6 6a2.5 2.5 0 1 1 0 .1zM18 6a2.5 2.5 0 1 1 0 .1zM6 18a2.5 2.5 0 1 1 0 .1zM18 18a2.5 2.5 0 1 1 0 .1zM8 8l8 8M16 8l-8 8', evtol: 'M4 7h6M14 7h6M7 7v4M17 7v4M3 13h18l-2 4H5z',
  leaf: 'M5 19c0-9 5-14 15-14 0 10-5 15-14 14zM5 19l8-8', clock: 'M12 4a8 8 0 1 1 0 16 8 8 0 0 1 0-16zM12 8v4l3 2', star: 'M12 4l2.4 5 5.6.7-4.1 3.8 1.1 5.5-5-2.8-5 2.8 1.1-5.5L4 9.7 9.6 9z', external: 'M14 5h5v5M19 5l-8 8M11 7H6v11h11v-5',
};
export function icon(name, size = 20) {
  return h('svg', { viewBox: '0 0 24 24', width: size, height: size, fill: 'none', stroke: 'currentColor', 'stroke-width': 1.7, 'stroke-linecap': 'round', 'stroke-linejoin': 'round', 'aria-hidden': 'true', class: 'ico' },
    h('path', { d: ICONS[name] || ICONS.info }));
}

let toastHost;
export function toast(msg, kind = 'info', ms = 4200) {
  if (!toastHost) { toastHost = h('div', { class: 'toasts', role: 'status', 'aria-live': 'polite' }); document.body.append(toastHost); }
  const t = h('div', { class: `toast ${kind}` }, icon(kind === 'bad' ? 'warn' : kind === 'ok' ? 'check' : 'info', 18), h('span', null, msg));
  toastHost.append(t);
  setTimeout(() => { t.classList.add('out'); setTimeout(() => t.remove(), 300); }, ms);
}

export function downloadText(name, text, type = 'text/plain') { downloadBlob(name, new Blob([text], { type })); }
export function downloadBlob(name, blob) {
  const a = h('a', { href: URL.createObjectURL(blob), download: name });
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 5000);
}
export const csvEscape = (v) => { const s = v == null ? '' : String(v); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
export const toCsv = (cols, rows) => [cols, ...rows].map((r) => r.map(csvEscape).join(',')).join('\n');

export function pickFiles({ multiple = false, accept = '' } = {}) {
  return new Promise((res) => { const i = h('input', { type: 'file', multiple, accept, style: { display: 'none' } }); i.addEventListener('change', () => { res([...i.files]); i.remove(); }); document.body.append(i); i.click(); });
}
export const debounce = (fn, ms = 250) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };

/** Section with a heading; optional collapsible body. */
export function card(title, body, { collapsible = false, open = true, actions = null, cls = '' } = {}) {
  if (!collapsible) return h('section', { class: `card ${cls}` }, title ? h('header', { class: 'card-h' }, h('h3', null, title), actions) : null, h('div', { class: 'card-b' }, body));
  return h('details', { class: `card ${cls}`, open }, h('summary', { class: 'card-h' }, h('h3', null, title), actions), h('div', { class: 'card-b' }, body));
}
export const badge = (text, kind = '') => h('span', { class: `badge ${kind}` }, text);
export const btn = (label, onclick, { ic, kind = '', title, disabled } = {}) => h('button', { class: `btn ${kind}`, type: 'button', onclick, title: title || null, disabled: disabled || null }, ic ? icon(ic, 18) : null, label ? h('span', null, label) : null);
export const empty = (ic, title, text, action) => h('div', { class: 'empty' }, icon(ic, 34), h('strong', null, title), text ? h('p', null, text) : null, action);
