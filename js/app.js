// Application shell: navigation, routing, theme, install prompt, offline service worker, command palette.

import { h, icon, clear, toast, add } from './ui/dom.js';
import { SUITES, GROUPS } from './core/registry.js';
import { state, on, setSetting } from './core/store.js';
import { startLive, status as liveStatus } from './core/live.js';
import { t, LANGS, lang, applyLang } from './ui/i18n.js';
import { unitSystem } from './ui/units.js';
import { pageGuide, toggleGuide } from './ui/guide.js';

const VIEWS = {
  home: () => import('./ui/views/home.js'), case: () => import('./ui/views/case.js'), geometry: () => import('./ui/views/geometry.js'),
  suite: () => import('./ui/views/suite.js'), integrated: () => import('./ui/views/integrated.js'), decisions: () => import('./ui/views/decisions.js'),
  live: () => import('./ui/views/livehub.js'), bridge: () => import('./ui/views/bridge.js'), share: () => import('./ui/views/share.js'), reports: () => import('./ui/views/reports.js'), about: () => import('./ui/views/about.js'),
};
const MAIN = [
  { path: 'home', label: 'Overview', ic: 'home' }, { path: 'case', label: 'Case & input portal', ic: 'sliders' }, { path: 'geometry', label: 'Geometry & mesh', ic: 'cube' },
  { path: 'integrated', label: 'Integrated run', ic: 'graph' }, { path: 'bridge', label: 'High-fidelity bridge', ic: 'external' }, { path: 'decisions', label: 'Decision support', ic: 'bulb' }, { path: 'live', label: 'Live data', ic: 'globe' },
  { path: 'reports', label: 'Reports & assurance', ic: 'doc' }, { path: 'about', label: 'Install, offline & about', ic: 'install' },
];
/** Linear page order used by the previous/next arrows at the foot of every page. */
export const SEQUENCE = [
  { hash: '#/home', label: 'Overview' }, { hash: '#/case', label: 'Case & input portal' }, { hash: '#/geometry', label: 'Geometry & mesh' },
  ...SUITES.map((s) => ({ hash: `#/suite/${s.id}`, label: `${s.d}. ${s.short}` })),
  { hash: '#/integrated', label: 'Integrated run' }, { hash: '#/bridge', label: 'High-fidelity bridge' }, { hash: '#/decisions', label: 'Decision support' }, { hash: '#/live', label: 'Live data' }, { hash: '#/reports', label: 'Reports & assurance' }, { hash: '#/about', label: 'Install, offline & about' },
];

const app = document.getElementById('app');
const page = h('main', { class: 'page', id: 'page', tabIndex: -1 });
const crumb = h('div', { class: 'crumb' });
const backBtn = h('button', { class: 'arrow', title: 'Back (Alt+←)', 'aria-label': 'Go back', onclick: () => history.back() }, icon('back', 20));
const fwdBtn = h('button', { class: 'arrow', title: 'Forward (Alt+→)', 'aria-label': 'Go forward', onclick: () => history.forward() }, icon('fwd', 20));
const livePill = h('a', { class: 'pill', href: '#/live', title: 'Live data status' }, h('i', { class: 'dot' }), h('span', { class: 'txt' }, 'Live data'));
const installBtn = h('button', { class: 'btn primary sm', hidden: true, onclick: () => installApp() }, icon('install', 16), h('span', { class: 'hide-sm' }, t('Install')));
const themeBtn = h('button', { class: 'icon-btn', title: 'Switch light / dark theme', 'aria-label': 'Switch theme', onclick: toggleTheme });
const sideNav = h('nav', { class: 'nav', 'aria-label': 'Main' });
const bottom = h('nav', { class: 'bottom', 'aria-label': 'Quick' });
const versionLabel = h('span', null, 'Version …');
const updateBtn = h('button', { class: 'btn sm', title: 'Look for a newer version of AeroSuite 26 now', onclick: () => checkForUpdate(true) }, icon('refresh', 16), h('span', null, 'Check for updates'));
const sideFoot = h('div', { class: 'side-foot' }, versionLabel, updateBtn);

function buildShell() {
  const side = h('aside', { class: 'side' },
    h('a', { class: 'brand', href: '#/home' }, h('span', { class: 'logo' }, icon('jet', 22)), h('span', null, h('b', null, 'AeroSuite 26'), h('small', null, 'Aircraft engineering simulation'))),
    sideNav, sideFoot);
  const top = h('header', { class: 'top' },
    h('button', { class: 'icon-btn menu-btn', 'aria-label': 'Open menu', onclick: () => app.classList.toggle('open') }, icon('menu')),
    backBtn, fwdBtn, crumb,
    h('button', { class: 'icon-btn', title: 'How this page works', 'aria-label': 'Help for this page', onclick: toggleGuide }, h('b', { style: { fontSize: '1.05rem' } }, '?')),
    h('button', { class: 'icon-btn', title: 'Search suites, analyses and pages (Ctrl+K or /)', 'aria-label': 'Search', onclick: openPalette }, icon('search')),
    livePill, installBtn,
    h('button', { class: 'pill hide-sm', title: 'Switch between SI units and aviation units (ft, kt, lb, nmi)', onclick: () => { setSetting('units', unitSystem() === 'si' ? 'aviation' : 'si'); buildShell(); route(); } }, unitSystem() === 'si' ? 'SI' : 'ft · kt · lb'),
    h('select', { class: 'pill', 'aria-label': 'Interface language', title: 'Interface language', onchange: (e) => { setSetting('lang', e.target.value); applyLang(); buildShell(); route(); } }, Object.entries(LANGS).map(([k, v]) => h('option', { value: k, selected: k === lang() }, v))),
    themeBtn);
  clear(app);
  add(app, [side, h('div', { class: 'scrim', onclick: () => app.classList.remove('open') }), h('div', { class: 'main' }, top, page), bottom]);
  renderNav(); paintTheme(); paintLive();
}

function renderNav() {
  const cur = location.hash || '#/home', links = [];
  links.push(...MAIN.map((m) => h('a', { href: `#/${m.path}`, class: cur.startsWith(`#/${m.path}`) ? 'on' : '' }, icon(m.ic, 18), h('span', { class: 't' }, t(m.label)))));
  for (const g of GROUPS) {
    links.push(h('div', { class: 'nav-h' }, g));
    for (const s of SUITES.filter((x) => x.group === g)) {
      const runs = Object.entries(state.runs).filter(([k]) => k.startsWith(s.id + '.')).map(([, v]) => v.status);
      const st = runs.includes('bad') ? 'bad' : runs.includes('warn') ? 'warn' : runs.length ? 'ok' : '';
      links.push(h('a', { href: `#/suite/${s.id}`, class: cur.startsWith(`#/suite/${s.id}`) ? 'on' : '', title: s.title }, h('span', { class: 'n' }, s.d), h('span', { class: 't' }, s.short), h('i', { class: `dot ${st}`, title: st ? `Last run: ${st}` : 'Not run yet' })));
    }
  }
  const keep = sideNav.scrollTop;                       // rebuilding the list must not move the menu
  clear(sideNav); add(sideNav, links);
  sideNav.scrollTop = keep;
  const cur_el = sideNav.querySelector('a.on'); if (cur_el) { const r = cur_el.getBoundingClientRect(), b = sideNav.getBoundingClientRect(); if (r.top < b.top || r.bottom > b.bottom) cur_el.scrollIntoView({ block: 'nearest' }); }
  clear(bottom);
  add(bottom, [['home', 'Overview', 'home'], ['case', 'Case', 'sliders'], ['suite/' + lastSuite(), 'Suites', 'layers'], ['integrated', 'Run all', 'graph'], ['decisions', 'Advice', 'bulb']].map(([p, l, ic]) => h('a', { href: `#/${p}`, class: cur.startsWith(`#/${p.split('/')[0]}`) ? 'on' : '' }, icon(ic, 21), t(l))));
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
    prev ? h('a', { href: prev.hash, rel: 'prev' }, h('span', { class: 'arrow' }, icon('back', 18)), h('span', null, h('small', null, t('Previous')), t(prev.label))) : h('span'),
    next ? h('a', { class: 'next', href: next.hash, rel: 'next' }, h('span', null, h('small', null, t('Next')), t(next.label)), h('span', { class: 'arrow' }, icon('fwd', 18))) : null);
}
async function route() {
  const hash = location.hash || '#/home', parts = hash.replace(/^#\/?/, '').split('/').map(decodeURIComponent), name = VIEWS[parts[0]] ? parts[0] : 'home', token = ++navToken;
  trackHistory(hash); app.classList.remove('open');
  if (name === 'suite' && parts[1]) { try { sessionStorage.setItem('lastSuite', parts[1]); } catch { /* ignore */ } }
  try { cleanup?.(); } catch (e) { console.error(e); } cleanup = null; delete page.dataset.route;
  renderNav();
  try {
    const mod = await VIEWS[name]();
    if (token !== navToken) return;
    clear(page);
    const res = await mod.render(page, parts.slice(1), { setCrumb: (t) => { crumb.textContent = t; document.title = `${t} · AeroSuite 26`; } });
    if (token !== navToken) return;
    cleanup = typeof res === 'function' ? res : null;
    paintGuide(parts);
    page.dataset.route = hash; // marks the page as fully rendered for this address (used by automated tests)
    const pg = pager(); if (pg) page.append(pg);
    if (!history.state?.keepScroll) window.scrollTo(0, 0);
  } catch (e) {
    console.error(e); clear(page);
    add(page, [h('div', { class: 'note bad' }, icon('warn'), h('div', null, h('b', null, 'This page could not be shown. '), e.message || String(e), navigator.onLine === false ? ' You are offline and this part of the app has not been stored on the device yet. Open it once while online.' : ''))]);
  }
}
function paintGuide(parts) { page.querySelector('details.guide')?.remove(); const g = pageGuide(parts); if (g) page.prepend(g); }
/** Change the hash without adding to history (used for tabs inside a page). */
export function replaceHash(hash) { history.replaceState({ keepScroll: true }, '', hash); trackHistory(hash); renderNav(); paintGuide(hash.replace(/^#\/?/, '').split('/')); }

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

// ---- service worker and updates ----------------------------------------------------------------
let swReg = null, runningVersion = null;
/** Ask the active service worker which build it is serving. */
function askVersion() {
  return new Promise((res) => { const c = navigator.serviceWorker?.controller; if (!c) return res(null); const on = (e) => { if (e.data?.type === 'version') { navigator.serviceWorker.removeEventListener('message', on); res(e.data.version); } }; navigator.serviceWorker.addEventListener('message', on); c.postMessage({ type: 'version' }); setTimeout(() => res(null), 1500); });
}
const applyUpdate = async () => {
  // switch to the waiting build; if none is waiting yet, fetch it first; as a last resort clear the stored copy and reload
  try { await swReg?.update(); } catch { /* offline */ }
  const w = swReg?.waiting || swReg?.installing;
  if (w) { if (w.state === 'installed') w.postMessage({ type: 'skip-waiting' }); else w.addEventListener('statechange', () => { if (w.state === 'installed') w.postMessage({ type: 'skip-waiting' }); }); }
  setTimeout(async () => { try { for (const k of await caches.keys()) if (k.startsWith('aerosuite-app-')) await caches.delete(k); await swReg?.unregister(); } catch { /* ignore */ } location.reload(); }, 12000);
};
/** Compare the running build with the published one; show the prompt if a newer one exists. */
export async function checkForUpdate(manual = false) {
  if (location.protocol === 'file:') { if (manual) toast('This is the single-file copy. Download a fresh copy from the website to update.', 'info'); return; }
  if (manual) { updateBtn.disabled = true; updateBtn.lastChild.textContent = 'Checking…'; }
  try {
    try { await swReg?.update(); } catch { /* offline */ }
    const latest = await fetch('version.json', { cache: 'no-store' }).then((r) => r.json()).catch(() => null);
    runningVersion = (await askVersion()) || runningVersion;
    versionLabel.textContent = runningVersion ? `Version ${runningVersion}` : latest?.version ? `Version ${latest.version}` : 'Version unknown';
    const newer = (swReg?.waiting && navigator.serviceWorker.controller) || (latest?.version && runningVersion && latest.version !== runningVersion && runningVersion !== 'dev');
    if (newer) showUpdatePrompt(applyUpdate, latest?.version);
    else if (manual) toast(navigator.onLine === false ? 'You are offline. Updates are checked when you reconnect.' : 'You have the latest version.', navigator.onLine === false ? 'info' : 'ok');
  } finally { if (manual) { updateBtn.disabled = false; updateBtn.lastChild.textContent = 'Check for updates'; } }
}
async function registerSW() {
  if (!('serviceWorker' in navigator) || location.protocol === 'file:') { versionLabel.textContent = location.protocol === 'file:' ? 'Single-file copy' : 'Version (no offline engine)'; return; }
  try {
    swReg = await navigator.serviceWorker.register('sw.js');
    const hadController = !!navigator.serviceWorker.controller; let reloading = false;
    swReg.addEventListener('updatefound', () => { const w = swReg.installing; w?.addEventListener('statechange', () => { if (w.state === 'installed' && navigator.serviceWorker.controller) checkForUpdate(); }); });
    navigator.serviceWorker.addEventListener('controllerchange', () => { if (hadController && !reloading) { reloading = true; location.reload(); } });
    await navigator.serviceWorker.ready; setTimeout(() => checkForUpdate(), 1200);
    setInterval(() => checkForUpdate(), 20 * 60e3);                                  // look for a newer version while open
    document.addEventListener('visibilitychange', () => { if (document.visibilityState === 'visible') checkForUpdate(); });
    // background refresh of live data when the installed app is closed (supported browsers only)
    try { const perm = await navigator.permissions.query({ name: 'periodic-background-sync' }); if (perm.state === 'granted' && swReg.periodicSync) await swReg.periodicSync.register('refresh-live', { minInterval: 6 * 3600e3 }); } catch { /* not supported */ }
  } catch (e) { console.warn('Service worker registration failed', e); }
}

/** Persistent banner inviting the person to switch to the newly published version. */
function showUpdatePrompt(apply, version) {
  if (document.querySelector('.update-bar')) return;
  updateBtn.classList.add('primary'); updateBtn.lastChild.textContent = 'Update available';
  const bar = h('div', { class: 'update-bar', role: 'alert' }, icon('refresh', 20),
    h('span', null, h('b', null, 'A new version of AeroSuite 26 is ready. '), version ? `(${version}) ` : '', 'Update to get the latest fixes and features. Your case and results are kept.'),
    h('button', { class: 'btn primary', onclick: (e) => { e.target.disabled = true; e.target.textContent = 'Updating…'; apply(); } }, 'Update now'),
    h('button', { class: 'btn ghost', onclick: () => bar.remove(), title: 'Keep working; use “Check for updates” in the menu whenever you are ready' }, 'Later'));
  document.body.append(bar);
}

// ---- command palette ----------------------------------------------------------------------
let paletteItems = null;
async function openPalette() {
  if (document.querySelector('.palette')) return;
  if (!paletteItems) {
    paletteItems = [...MAIN.map((m) => ({ label: m.label, hint: 'Page', hash: `#/${m.path}`, ic: m.ic })), ...SUITES.map((s) => ({ label: `${s.d}. ${s.title}`, hint: s.group, hash: `#/suite/${s.id}`, ic: s.icon }))];
    import('./core/jobs.js').then(async ({ execute }) => { // enrich with every analysis, lazily
      for (const s of SUITES) { try { const d = await execute({ kind: 'describe', suite: s.id, case: state.case, up: state.up }); for (const a of d.analyses) paletteItems.push({ label: a.title, hint: `${s.d}. ${s.short}`, hash: `#/suite/${s.id}/${a.id}`, ic: s.icon, extra: a.summary }); } catch { /* suite not available */ } }
    });
  }
  const input = h('input', { type: 'search', placeholder: t('Search suites, analyses and pages…'), 'aria-label': 'Search', autocomplete: 'off' }), list = h('ul', { role: 'listbox' });
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
applyLang();
buildShell();
window.addEventListener('hashchange', route);
window.addEventListener('online', () => { paintLive(); toast('Back online. Refreshing live data.', 'ok'); });
window.addEventListener('offline', () => { paintLive(); toast('You are offline. Everything keeps working with the data already saved on this device.', 'info'); });
on('live', paintLive); on('run', renderNav);
route();
registerSW();
startLive();
