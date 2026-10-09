// Install, offline & about: how to put the app on any device, how it stays available, what it can
// and cannot do, and where its data comes from.

import { h, icon, clear, add, num, btn, badge, card, toast } from '../dom.js';
import { installApp, installState } from '../../app.js';
import { SUITES } from '../../core/registry.js';
import { state } from '../../core/store.js';
import { runJob, usingWorker } from '../../core/runner.js';
import { CONNECTORS } from '../../core/live.js';
import { licenceSection } from './livehub.js';

export async function render(root, _p, { setCrumb }) {
  setCrumb('Install, offline & about');
  const host = h('div'); root.append(host);
  const ins = installState(), ua = navigator.userAgent, android = /android/i.test(ua), mac = /macintosh/i.test(ua) && !ins.ios, firefox = /firefox/i.test(ua), safari = /safari/i.test(ua) && !/chrome|chromium|edg|crios|fxios/i.test(ua);
  const licHost = h('div', { id: 'licences' }), offHost = h('dl', { class: 'kv' }), mirHost = h('div', { class: 'stack' }), capHost = h('div', null, h('div', { class: 'row muted' }, h('i', { class: 'spin' }), 'Reading the 26 suites…'));

  const steps = ins.standalone ? [h('div', { class: 'note ok' }, icon('check'), h('div', null, 'AeroSuite 26 is installed and running as an app on this device.'))]
    : ins.ios ? [h('ol', null, h('li', null, 'Open this page in ', h('b', null, 'Safari'), ' (on iOS 16.4 or later Chrome and Edge also work).'), h('li', null, 'Tap the ', h('b', null, 'Share'), ' button (the square with an upward arrow).'), h('li', null, 'Choose ', h('b', null, 'Add to Home Screen'), ', then ', h('b', null, 'Add'), '.'), h('li', null, 'Open AeroSuite 26 from the home screen. It now works in aeroplane mode.'))]
    : ins.canPrompt ? [h('p', null, 'This browser can install the app directly.'), btn('Install AeroSuite 26', installApp, { ic: 'install', kind: 'primary big' })]
    : android ? [h('ol', null, h('li', null, 'Open the browser menu (⋮).'), h('li', null, 'Tap ', h('b', null, 'Install app'), ' or ', h('b', null, 'Add to Home screen'), '.'), h('li', null, 'Confirm. The app appears in your launcher and works offline.'))]
    : safari && mac ? [h('ol', null, h('li', null, 'In Safari choose ', h('b', null, 'File → Add to Dock'), ' (macOS Sonoma or later).'), h('li', null, 'Open AeroSuite 26 from the Dock or Launchpad.'))]
    : firefox ? [h('p', null, 'Firefox on desktop does not install web apps, but the app is still stored for offline use after this first visit: bookmark this page and it will open without a connection. For an app window, open this address in Chrome, Edge or Safari, or download the single-file copy below.')]
    : [h('ol', null, h('li', null, 'In Chrome or Edge, click the install icon at the right of the address bar, or open the menu and choose ', h('b', null, 'Install AeroSuite 26'), '.'), h('li', null, 'The app opens in its own window and is added to your Start menu, Dock or app list.'))];

  host.append(
    h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Any device, with or without a connection'), h('h1', null, 'Install, offline & about'), h('p', null, 'AeroSuite 26 is a web application that installs like a native one on Windows, macOS, Linux, ChromeOS, Android, iPhone and iPad. After the first visit all solvers, pages and reference data are stored on the device, so it starts instantly and keeps working in aeroplane mode.'))),
    h('div', { class: 'grid g2' },
      card('Install on this device', h('div', { class: 'stack' }, ...steps,
        h('div', { class: 'row' }, h('a', { class: 'btn', href: 'standalone.html', download: 'AeroSuite26.html' }, icon('download', 18), 'Download single-file copy'), h('span', { class: 'muted small' }, 'One HTML file with the whole app inside. Keep it on a USB stick or shared drive and open it in any browser — no installation, no server, no connection.')))),
      card('Offline readiness', h('div', { class: 'stack' }, offHost, h('div', { class: 'row' }, btn('Store everything for offline use now', precache, { ic: 'download' }), btn('Check for an update', async () => { const r = await navigator.serviceWorker?.getRegistration(); if (!r) { toast('Updates are checked automatically when the app is served from the web.', 'info'); return; } await r.update(); toast(r.waiting || r.installing ? 'An update is being installed.' : 'You have the latest version.', 'ok'); }, { ic: 'refresh', kind: 'ghost' }))))),
    h('div', { class: 'gap' }),
    card('Always available: mirrors', h('div', { class: 'stack' }, h('p', { class: 'muted small' }, 'The same build is published to several independent hosts. If one provider has an outage, open another address — and once installed, the app does not need any of them to run.'), mirHost)),
    h('div', { class: 'gap' }),
    card('What each suite computes on this device', capHost),
    h('div', { class: 'gap' }),
    h('div', { class: 'grid g2' },
      card('Capabilities and limits', h('div', { class: 'stack small' },
        h('p', null, h('b', null, 'What it is. '), 'An integrated analysis environment covering 26 engineering disciplines for aeroplanes, helicopters, rotorcraft, UAVs and electric VTOL aircraft, built on classical governing equations, reduced-order and numerical solvers that run in real time in a browser, coupled through one shared aircraft case, with mesh-convergence, uncertainty, calibration and validation tools on every analysis.'),
        h('p', null, h('b', null, 'What it is not. '), 'It is not a CAD system, and it does not run three-dimensional RANS, LES or DNS, full-aircraft nonlinear finite-element crash simulations, or coupled high-fidelity CFD–structure solutions inside the browser. Where the specification names such models, each suite lists them under “Models handed off” with the class of external tool to use, and accepts their results back for calibration and validation.'),
        h('p', null, h('b', null, 'Credibility. '), 'Verification benchmarks show the equations are solved correctly. They do not show a model represents your aircraft: that needs validation against independent test data within the intended operating envelope. Default material properties, failure rates, cost coefficients and correlations are typical or illustrative values, labelled as such, to be replaced with sourced data. Nothing here is airworthiness evidence on its own.'))),
      card('Privacy, security and data sources', h('div', { class: 'stack small' },
        h('p', null, h('b', null, 'Your data stays with you. '), 'Cases, geometry, results and reports are stored only in this browser on this device. Nothing is uploaded. Imported files are parsed locally and treated as untrusted data.'),
        h('p', null, h('b', null, 'Locked down. '), 'A strict content-security policy allows scripts only from the app itself and connections only to the listed public data providers. There are no trackers, adverts, accounts or third-party scripts.'),
        h('p', null, h('b', null, 'Live data providers: '), Object.values(CONNECTORS).map((c, i) => [i ? ' · ' : '', h('a', { href: c.home, target: '_blank', rel: 'noopener noreferrer' }, c.provider)]), '. Terms, licences and required attributions are listed under “Data sources and licences” below.'),
        h('p', null, h('b', null, 'Licence. '), 'AeroSuite 26 is proprietary software: Copyright © 2026 Samuel Akosa Onyejekwe. All rights reserved (', h('a', { href: 'LICENSE', target: '_blank', rel: 'noopener' }, 'LICENSE'), '). It includes third-party components under their own licences, among them the OpenCASCADE geometry kernel under the GNU LGPL 2.1 (', h('a', { href: 'THIRD_PARTY_NOTICES.md', target: '_blank', rel: 'noopener' }, 'THIRD_PARTY_NOTICES.md'), ').'),
        h('p', null, h('b', null, 'Shortcuts: '), h('span', { class: 'mono' }, 'Ctrl/⌘ K'), ' or ', h('span', { class: 'mono' }, '/'), ' search · ', h('span', { class: 'mono' }, 'Alt ←'), ' back · ', h('span', { class: 'mono' }, 'Alt →'), ' forward.')))),
    h('div', { class: 'gap' }), licHost);
  licenceSection({ pack: true }).then((el) => licHost.append(el)).catch(() => {});

  async function paintOffline() {
    const reg = await navigator.serviceWorker?.getRegistration().catch(() => null); let files = 0, bytes = null;
    try { for (const n of await caches.keys()) files += (await (await caches.open(n)).keys()).length; } catch { /* no cache api */ }
    try { const e = await navigator.storage?.estimate(); bytes = e?.usage; } catch { /* unsupported */ }
    let persisted = false; try { persisted = await navigator.storage?.persisted(); } catch { /* unsupported */ }
    clear(offHost); add(offHost, [
      h('dt', null, 'Offline engine'), h('dd', null, location.protocol === 'file:' ? badge('Single-file copy: fully offline', 'ok') : reg?.active ? badge('Active', 'ok') : badge('Starting — reload once', 'warn')),
      h('dt', null, 'Files stored'), h('dd', null, files ? `${files}` : 'none yet'), h('dt', null, 'Storage used'), h('dd', null, bytes != null ? `${num(bytes / 1e6, 3)} MB` : 'unknown'),
      h('dt', null, 'Protected from clean-up'), h('dd', null, persisted ? badge('Yes', 'ok') : h('span', null, 'No ', btn('Protect', async () => { const ok = await navigator.storage?.persist?.(); toast(ok ? 'Storage is now protected from automatic clean-up.' : 'The browser declined; installing the app usually grants this.', ok ? 'ok' : 'info'); paintOffline(); }, { kind: 'sm ghost' }))),
      h('dt', null, 'Connection'), h('dd', null, navigator.onLine === false ? badge('Offline', 'warn') : badge('Online', 'ok')), h('dt', null, 'Solver thread'), h('dd', null, usingWorker() ? 'Background worker' : 'Main thread'),
    ]);
  }
  async function precache() {
    const reg = await navigator.serviceWorker?.getRegistration().catch(() => null);
    if (!reg?.active) { toast(location.protocol === 'file:' ? 'The single-file copy is already complete.' : 'The offline engine is still starting. Reload the page once, then try again.', 'info'); return; }
    reg.active.postMessage({ type: 'precache-all' }); toast('Storing every suite and data file for offline use…', 'info'); setTimeout(paintOffline, 4000);
  }
  navigator.serviceWorker?.addEventListener('message', (e) => { if (e.data?.type === 'precache-done') { toast(`Ready for offline use: ${e.data.count} files stored.`, 'ok'); paintOffline(); } });

  async function paintMirrors() {
    let list = []; let ver = null;
    try { list = (await (await fetch('mirrors.json', { cache: 'no-store' })).json()).mirrors || []; } catch { /* offline or standalone */ }
    try { ver = await (await fetch('version.json', { cache: 'no-store' })).json(); } catch { /* offline */ }
    clear(mirHost);
    if (ver) mirHost.append(h('p', { class: 'small' }, 'This copy: version ', h('b', null, ver.version), ` built ${ver.built}.`));
    if (!list.length) { mirHost.append(h('p', { class: 'muted small' }, 'No mirror addresses are configured in this build yet (mirrors.json). The deployment guide in the repository explains how to publish to several hosts in one step.')); return; }
    for (const m of list) {
      const st = h('span', { class: 'badge' }, 'checking…'), here = location.href.startsWith(m.url);
      mirHost.append(h('div', { class: 'row' }, h('a', { href: m.url, target: '_blank', rel: 'noopener noreferrer' }, m.name), h('span', { class: 'muted small grow' }, m.url), here ? badge('You are here', 'accent') : null, st));
      fetch(m.url.replace(/\/?$/, '/') + 'version.json', { cache: 'no-store', mode: 'cors' }).then((r) => (r.ok ? r.json() : Promise.reject())).then((v) => { st.textContent = `Up · v${v.version}`; st.className = 'badge ok'; }).catch(() => { st.textContent = navigator.onLine === false ? 'Offline' : 'Not reachable'; st.className = 'badge warn'; });
    }
  }
  async function paintCaps() {
    const rows = [];
    for (const s of SUITES) { try { const d = await runJob({ kind: 'describe', suite: s.id, case: state.case, up: state.up }); const f = (k) => d.analyses.filter((a) => a.fidelity === k).length; rows.push(h('tr', null, h('td', null, h('a', { href: `#/suite/${s.id}` }, `${s.d}. ${s.short}`)), h('td', { class: 'num' }, d.analyses.length), h('td', { class: 'num' }, f('numerical')), h('td', { class: 'num' }, f('reduced-order')), h('td', { class: 'num' }, f('analytical')), h('td', { class: 'num' }, d.analyses.filter((a) => a.convergence).length), h('td', { class: 'num' }, d.handoff.length), h('td', null, d.analyses.map((a) => a.title).join(' · ')))); } catch { rows.push(h('tr', null, h('td', null, `${s.d}. ${s.short}`), h('td', { colSpan: 7 }, 'not available in this build'))); } }
    clear(capHost); capHost.append(h('div', { class: 'table-wrap' }, h('table', { class: 'data' }, h('thead', null, h('tr', null, ['Suite', 'Analyses', 'Numerical solvers', 'Reduced-order', 'Closed-form', 'With mesh study', 'Handed off', 'Analyses included'].map((x, i) => h('th', { class: i > 0 && i < 7 ? 'num' : '' }, x)))), h('tbody', null, rows))));
  }
  paintOffline(); paintMirrors(); paintCaps();
}
