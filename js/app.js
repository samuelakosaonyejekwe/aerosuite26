// Application shell: navigation, routing, theme, install prompt, offline service worker, command palette.

import { h, icon, clear, toast, add } from './ui/dom.js';
import { SUITES, GROUPS } from './core/registry.js';
import { state, on, setSetting } from './core/store.js';
import { startLive, status as liveStatus } from './core/live.js';

const VIEWS = {
  home: () => import('./ui/views/home.js'), case: () => import('./ui/views/case.js'), geometry: () => import('./ui/views/geometry.js'),
  suite: () => import('./ui/views/suite.js'), integrated: () => import('./ui/views/integrated.js'), decisions: () => import('./ui/views/decisions.js'),
  live: () => import('./ui/views/livehub.js'), reports: () => import('./ui/views/reports.js'), about: () => import('./ui/views/about.js'),
};
const MAIN = [
  { path: 'home', label: 'Overview', ic: 'home' }, { path: 'case', label: 'Case & input portal', ic: 'sliders' }, { path: 'geometry', label: 'Geometry & mesh', ic: 'cube' },
  { path: 'integrated', label: 'Integrated run', ic: 'graph' }, { path: 'decisions', label: 'Decision support', ic: 'bulb' }, { path: 'live', label: 'Live data', ic: 'globe' },
  { path: 'reports', label: 'Reports & assurance', ic: 'doc' }, { path: 'about', label: 'Install, offline & about', ic: 'install' },
];
/** Linear page order used by the previous/next arrows at the foot of every page. */
export const SEQUENCE = [
  { hash: '#/home', label: 'Overview' }, { hash: '#/case', label: 'Case & input portal' }, { hash: '#/geometry', label: 'Geometry & mesh' },
  ...SUITES.map((s) => ({ hash: `#/suite/${s.id}`, label: `${s.n}. ${s.short}` })),
  { hash: '#/integrated', label: 'Integrated run' }, { hash: '#/decisions', label: 'Decision support' }, { hash: '#/live', label: 'Live data' }, { hash: '#/reports', label: 'Reports & assurance' }, { hash: '#/about', label: 'Install, offline & about' },
];

const app = document.getElementById('app');
const page = h('main', { class: 'page', id: 'page', tabIndex: -1 });
const crumb = h('div', { class: 'crumb' });
const backBtn = h('button', { class: 'arrow', title: 'Back (Alt+←)', 'aria-label': 'Go back', onclick: () => history.back() }, icon('back', 20));
const fwdBtn = h('button', { class: 'arrow', title: 'Forward (Alt+→)', 'aria-label': 'Go forward', onclick: () => history.forward() }, icon('fwd', 20));
const livePill = h('a', { class: 'pill', href: '#/live', title: 'Live data status' }, h('i', { class: 'dot' }), h('span', { class: 'txt' }, 'Live data'));
const installBtn = h('button', { class: 'btn primary sm', hidden: true, onclick: () => installApp() }, icon('install', 16), h('span', { class: 'hide-sm' }, 'Install'));
const themeBtn = h('button', { class: 'icon-btn', title: 'Switch light / dark theme', 'aria-label': 'Switch theme', onclick: toggleTheme });
const sideNav = h('nav', { class: 'nav', 'aria-label': 'Main' });
const bottom = h('nav', { class: 'bottom', 'aria-label': 'Quick' });

function buildShell() {
  const side = h('aside', { class: 'side' },
    h('a', { class: 'brand', href: '#/home' }, h('span', { class: 'logo' }, icon('jet', 22)), h('span', null, h('b', null, 'AeroSuite 26'), h('small', null, 'Aircraft engineering simulation'))),
    sideNav);
  const top = h('header', { class: 'top' },
    h('button', { class: 'icon-btn menu-btn', 'aria-label': 'Open menu', onclick: () => app.classList.toggle('open') }, icon('menu')),
    backBtn, fwdBtn, crumb,
    h('button', { class: 'icon-btn', title: 'Search suites, analyses and pages (Ctrl+K or /)', 'aria-label': 'Search', onclick: openPalette }, icon('search')),
    livePill, installBtn, themeBtn);
  clear(app);
  add(app, [side, h('div', { class: 'scrim', onclick: () => app.classList.remove('open') }), h('div', { class: 'main' }, top, page), bottom]);
  renderNav(); paintTheme(); paintLive();
}

function renderNav() {
  const cur = location.hash || '#/home', links = [];
  links.push(...MAIN.map((m) => h('a', { href: `#/${m.path}`, class: cur.startsWith(`#/${m.path}`) ? 'on' : '' }, icon(m.ic, 18), h('span', { class: 't' }, m.label))));
  for (const g of GROUPS) {
    links.push(h('div', { class: 'nav-h' }, g));
    for (const s of SUITES.filter((x) => x.group === g)) {
      const runs = Object.entries(state.runs).filter(([k]) => k.startsWith(s.id + '.')).map(([, v]) => v.status);
      const st = runs.includes('bad') ? 'bad' : runs.includes('warn') ? 'warn' : runs.length ? 'ok' : '';
      links.push(h('a', { href: `#/suite/${s.id}`, class: cur.startsWith(`#/suite/${s.id}`) ? 'on' : '', title: s.title }, h('span', { class: 'n' }, s.n), h('span', { class: 't' }, s.short), h('i', { class: `dot ${st}`, title: st ? `Last run: ${st}` : 'Not run yet' })));
    }
  }
  clear(sideNav); add(sideNav, links);
  clear(bottom);
  add(bottom, [['home', 'Overview', 'home'], ['case', 'Case', 'sliders'], ['suite/' + lastSuite(), 'Suites', 'layers'], ['integrated', 'Run all', 'graph'], ['decisions', 'Advice', 'bulb']].map(([p, l, ic]) => h('a', { href: `#/${p}`, class: cur.startsWith(`#/${p.split('/')[0]}`) ? 'on' : '' }, icon(ic, 21), l)));
}
const lastSuite = () => { try { return sessionStorage.getItem('lastSuite') || 'cfd'; } catch { return 'cfd'; } };

// ---- router -----------------------------------------------------------------------------
let cleanup = null, navToken = 0;
const stack = []; let ptr = -1;
function trackHistory(hash) {
  if (stack[ptr] === hash) return;
  if (ptr > 0 && stack[ptr - 1] === hash) ptr--;
  else if (ptr < stack.length - 1 && stack[ptr + 1] === hash) ptr++;
  else { stack.length = ptr + 1; stack.push(hash); ptr++; }
  backBtn.disabled = ptr <= 0 && history.length <= 1; fwdBtn.disabled = ptr >= stack.length - 1;
}
export function pager() {
  const base = (location.hash || '#/home').split('/').slice(0, location.hash.startsWith('#/suite') ? 3 : 2).join('/');
  const i = SEQUENCE.findIndex((s) => s.hash === base); if (i < 0) return null;
  const prev = SEQUENCE[i - 1], next = SEQUENCE[i + 1];
  return h('nav', { class: 'pager', 'aria-label': 'Previous and next page' },
    prev ? h('a', { href: prev.hash, rel: 'prev' }, h('span', { class: 'arrow' }, icon('back', 18)), h('span', null, h('small', null, 'Previous'), prev.label)) : h('span'),
    next ? h('a', { class: 'next', href: next.hash, rel: 'next' }, h('span', null, h('small', null, 'Next'), next.label), h('span', { class: 'arrow' }, icon('fwd', 18))) : null);
}
async function route() {
  const hash = location.hash || '#/home', parts = hash.replace(/^#\/?/, '').split('/').map(decodeURIComponent), name = VIEWS[parts[0]] ? parts[0] : 'home', token = ++navToken;
  trackHistory(hash); app.classList.remove('open');
  if (name === 'suite' && parts[1]) { try { sessionStorage.setItem('lastSuite', parts[1]); } catch { /* ignore */ } }
  try { cleanup?.(); } catch (e) { console.error(e); } cleanup = null;
  renderNav();
  try {
    const mod = await VIEWS[name]();
    if (token !== navToken) return;
    clear(page);
    const res = await mod.render(page, parts.slice(1), { setCrumb: (t) => { crumb.textContent = t; document.title = `${t} · AeroSuite 26`; } });
    if (token !== navToken) return;
    cleanup = typeof res === 'function' ? res : null;
    const pg = pager(); if (pg) page.append(pg);
    if (!history.state?.keepScroll) window.scrollTo(0, 0);
  } catch (e) {
    console.error(e); clear(page);
    add(page, [h('div', { class: 'note bad' }, icon('warn'), h('div', null, h('b', null, 'This page could not be shown. '), e.message || String(e), navigator.onLine === false ? ' You are offline and this part of the app has not been stored on the device yet. Open it once while online.' : ''))]);
  }
}
/** Change the hash without adding to history (used for tabs inside a page). */
export function replaceHash(hash) { history.replaceState({ keepScroll: true }, '', hash); trackHistory(hash); renderNav(); }

// ---- theme --------------------------------------------------------------------------------
function paintTheme() { const dark = document.documentElement.dataset.theme === 'dark'; clear(themeBtn).append(icon(dark ? 'sun' : 'moon')); document.querySelector('meta[name=theme-color]')?.setAttribute('content', dark ? '#0e1116' : '#1f6fd1'); }
function toggleTheme() { const t = document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark'; document.documentElement.dataset.theme = t; setSetting('theme', t); paintTheme(); window.dispatchEvent(new Event('themechange')); }

// ---- live pill ----------------------------------------------------------------------------
function paintLive() {
  const on = navigator.onLine !== false, vals = Object.values(liveStatus), ok = vals.filter((v) => v.ok).length, failed = vals.filter((v) => v.ok === false).length;
  livePill.className = `pill ${on ? 'live' : 'off'}`;
  livePill.lastChild.textContent = on ? (vals.length ? `Live · ${ok} feed${ok === 1 ? '' : 's'}${failed ? `, ${failed} unavailable` : ''}` : 'Live data') : 'Offline · using saved data';
}

// ---- install ------------------------------------------------------------------------------
let deferredPrompt = null;
export const installState = () => ({ canPrompt: !!deferredPrompt, standalone: matchMedia('(display-mode: standalone)').matches || navigator.standalone === true, ios: /iphone|ipad|ipod/i.test(navigator.userAgent) || (navigator.platform === 'MacIntel' && navigator.maxTouchPoints > 1) });
export async function installApp() {
  if (deferredPrompt) { deferredPrompt.prompt(); const r = await deferredPrompt.userChoice.catch(() => null); deferredPrompt = null; installBtn.hidden = true; if (r?.outcome === 'accepted') toast('Installed. AeroSuite 26 now opens from your home screen or app list, with or without a connection.', 'ok'); return; }
  location.hash = '#/about';
}
window.addEventListener('beforeinstallprompt', (e) => { e.preventDefault(); deferredPrompt = e; installBtn.hidden = false; });
window.addEventListener('appinstalled', () => { deferredPrompt = null; installBtn.hidden = true; });

// ---- service worker -----------------------------------------------------------------------
async function registerSW() {
  if (!('serviceWorker' in navigator) || location.protocol === 'file:') return;
  try {
    const reg = await navigator.serviceWorker.register('sw.js');
    const offer = (w) => { if (!w) return; const t = h('button', { class: 'btn sm', onclick: () => w.postMessage({ type: 'skip-waiting' }) }, 'Reload'); toast('A new version is ready.', 'info', 12000); document.querySelector('.toasts .toast:last-child')?.append(t); };
    if (reg.waiting && navigator.serviceWorker.controller) offer(reg.waiting);
    reg.addEventListener('updatefound', () => { const w = reg.installing; w?.addEventListener('statechange', () => { if (w.state === 'installed' && navigator.serviceWorker.controller) offer(w); }); });
    let reloaded = false; navigator.serviceWorker.addEventListener('controllerchange', () => { if (!reloaded) { reloaded = true; location.reload(); } });
    setInterval(() => reg.update().catch(() => {}), 30 * 60e3);                     // look for app updates while open
    document.addEventListener('visibilitychange', () => { if (document.visibilityState === 'visible') reg.update().catch(() => {}); });
    // background refresh of live data when the installed app is closed (supported browsers only)
    try { const perm = await navigator.permissions.query({ name: 'periodic-background-sync' }); if (perm.state === 'granted' && reg.periodicSync) await reg.periodicSync.register('refresh-live', { minInterval: 6 * 3600e3 }); } catch { /* not supported */ }
  } catch (e) { console.warn('Service worker registration failed', e); }
}

// ---- command palette ----------------------------------------------------------------------
let paletteItems = null;
async function openPalette() {
  if (document.querySelector('.palette')) return;
  if (!paletteItems) {
    paletteItems = [...MAIN.map((m) => ({ label: m.label, hint: 'Page', hash: `#/${m.path}`, ic: m.ic })), ...SUITES.map((s) => ({ label: `${s.n}. ${s.title}`, hint: s.group, hash: `#/suite/${s.id}`, ic: s.icon }))];
    import('./core/jobs.js').then(async ({ execute }) => { // enrich with every analysis, lazily
      for (const s of SUITES) { try { const d = await execute({ kind: 'describe', suite: s.id, case: state.case, up: state.up }); for (const a of d.analyses) paletteItems.push({ label: a.title, hint: `${s.n}. ${s.short}`, hash: `#/suite/${s.id}/${a.id}`, ic: s.icon, extra: a.summary }); } catch { /* suite not available */ } }
    });
  }
  const input = h('input', { type: 'search', placeholder: 'Search suites, analyses and pages…', 'aria-label': 'Search', autocomplete: 'off' }), list = h('ul', { role: 'listbox' });
  const el = h('div', { class: 'palette', onclick: (e) => { if (e.target === el) close(); } }, h('div', { class: 'box' }, input, list));
  let sel = 0, shown = [];
  const close = () => el.remove();
  const go = (it) => { close(); location.hash = it.hash; };
  const paint = () => {
    const q = input.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
    shown = paletteItems.filter((it) => q.every((w) => (it.label + ' ' + it.hint + ' ' + (it.extra || '')).toLowerCase().includes(w))).slice(0, 40); sel = Math.min(sel, Math.max(0, shown.length - 1));
    clear(list); add(list, shown.length ? shown.map((it, i) => h('li', { class: i === sel ? 'on' : '', role: 'option', onclick: () => go(it) }, icon(it.ic, 18), h('span', null, it.label), h('small', null, it.hint))) : [h('li', null, 'Nothing matches. Try “flutter”, “take-off”, “NPV”…')]);
    list.querySelector('.on')?.scrollIntoView({ block: 'nearest' });
  };
  input.addEventListener('input', () => { sel = 0; paint(); });
  input.addEventListener('keydown', (e) => { if (e.key === 'ArrowDown') { sel = Math.min(shown.length - 1, sel + 1); paint(); e.preventDefault(); } else if (e.key === 'ArrowUp') { sel = Math.max(0, sel - 1); paint(); e.preventDefault(); } else if (e.key === 'Enter' && shown[sel]) go(shown[sel]); else if (e.key === 'Escape') close(); });
  document.body.append(el); paint(); input.focus();
}

// ---- boot ---------------------------------------------------------------------------------
window.addEventListener('keydown', (e) => {
  const typing = /^(INPUT|TEXTAREA|SELECT)$/.test(document.activeElement?.tagName || '');
  if ((e.key === 'k' && (e.ctrlKey || e.metaKey)) || (e.key === '/' && !typing)) { e.preventDefault(); openPalette(); }
  if (e.altKey && e.key === 'ArrowLeft') history.back();
  if (e.altKey && e.key === 'ArrowRight') history.forward();
});
buildShell();
window.addEventListener('hashchange', route);
window.addEventListener('online', () => { paintLive(); toast('Back online. Refreshing live data.', 'ok'); });
window.addEventListener('offline', () => { paintLive(); toast('You are offline. Everything keeps working with the data already saved on this device.', 'info'); });
on('live', paintLive); on('run', renderNav);
route();
registerSW();
startLive();
