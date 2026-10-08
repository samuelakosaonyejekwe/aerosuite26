// Offline engine. Stores the whole application on the device, serves it instantly from that copy,
// refreshes it in the background when a newer build is published, and keeps the last answer from each
// live-data provider so the app still has data in aeroplane mode.

/* BUILD:BEGIN */
const VERSION = "1.0.0+7225dec049";
const ASSETS = ["./","assets/icon-180.png","assets/icon-192.png","assets/icon-512.png","assets/icon-maskable-512.png","assets/icon.svg","css/app.css","index.html","js/app.js","js/boot-theme.js","js/core/atmosphere.js","js/core/case.js","js/core/geometry/analysis.js","js/core/geometry/formats.js","js/core/geometry/index.js","js/core/geometry/parsers-cad.js","js/core/geometry/parsers-mesh.js","js/core/geometry/parsers-surface.js","js/core/geometry/parsers-util.js","js/core/geometry/parsers.js","js/core/geometry/section.js","js/core/jobs.js","js/core/live.js","js/core/numerics.js","js/core/registry.js","js/core/runner.js","js/core/spec.js","js/core/store.js","js/core/studies.js","js/core/worker.js","js/data/materials.js","js/data/spec.json","js/suites/s01-cfd.js","js/suites/s02-fea.js","js/suites/s03-aeroelastic.js","js/suites/s04-flightdyn.js","js/suites/s05-performance.js","js/suites/s06-rotorcraft.js","js/suites/s07-propulsion.js","js/suites/s08-propeller.js","js/suites/s09-fatigue.js","js/suites/s10-vibration.js","js/suites/s11-acoustics.js","js/suites/s12-thermal.js","js/suites/s13-icing.js","js/suites/s14-gear.js","js/suites/s15-crash.js","js/suites/s16-control.js","js/suites/s17-avionics.js","js/suites/s18-hydmech.js","js/suites/s19-electrical.js","js/suites/s20-fuelecs.js","js/suites/s21-composites.js","js/suites/s22-safety.js","js/suites/s23-mdao.js","js/suites/s24-mission.js","js/suites/s25-vvuq.js","js/suites/s26-economics.js","js/ui/dom.js","js/ui/plots.js","js/ui/results.js","js/ui/viewer.js","js/ui/views/about.js","js/ui/views/case.js","js/ui/views/decisions.js","js/ui/views/geometry.js","js/ui/views/home.js","js/ui/views/integrated.js","js/ui/views/livehub.js","js/ui/views/reports.js","js/ui/views/suite.js","manifest.webmanifest","mirrors.json","standalone.html"];
/* BUILD:END */

const APP_CACHE = `aerosuite-app-${VERSION}`, LIVE_CACHE = 'aerosuite-live-v1';
const LIVE_HOSTS = ['api.open-meteo.com', 'archive-api.open-meteo.com', 'air-quality-api.open-meteo.com', 'marine-api.open-meteo.com', 'geocoding-api.open-meteo.com', 'services.swpc.noaa.gov', 'api.frankfurter.dev', 'api.worldbank.org', 'raw.githubusercontent.com', 'api.carbonintensity.org.uk', 'api.openalex.org', 'api.github.com', 'overpass-api.de', 'overpass.kumi.systems', 'overpass.private.coffee', 'maps.mail.ru'];

async function precache(all) {
  const cache = await caches.open(APP_CACHE), list = all ? ASSETS : ASSETS.filter((a) => !/^js\/suites\/|^js\/data\/|standalone\.html$/.test(a));
  let count = 0;
  // fetch individually so one missing file never blocks installation
  await Promise.all(list.map(async (url) => { try { if (all && (await cache.match(url))) { count++; return; } const r = await fetch(url, { cache: 'reload' }); if (r.ok) { await cache.put(url, r); count++; } } catch { /* offline during install: picked up later */ } }));
  return count;
}

self.addEventListener('install', (e) => { e.waitUntil(precache(false)); });
self.addEventListener('activate', (e) => {
  e.waitUntil((async () => {
    for (const k of await caches.keys()) if (k.startsWith('aerosuite-app-') && k !== APP_CACHE) await caches.delete(k);
    await self.clients.claim();
    precache(true).then((count) => self.clients.matchAll().then((cs) => cs.forEach((c) => c.postMessage({ type: 'precache-done', count })))); // everything else, without delaying start-up
  })());
});
self.addEventListener('message', (e) => {
  if (e.data?.type === 'skip-waiting') self.skipWaiting();
  if (e.data?.type === 'precache-all') e.waitUntil(precache(true).then((count) => e.source?.postMessage({ type: 'precache-done', count })));
});

self.addEventListener('fetch', (e) => {
  const req = e.request; if (req.method !== 'GET') return;
  const url = new URL(req.url);
  if (url.origin === self.location.origin) { e.respondWith(appAsset(req, url, e)); return; }
  if (LIVE_HOSTS.includes(url.hostname)) e.respondWith(liveData(req));
});

/** App files: answer from the device immediately, refresh the stored copy in the background. */
async function appAsset(req, url, e) {
  const cache = await caches.open(APP_CACHE), nav = req.mode === 'navigate';
  const key = nav ? 'index.html' : req;
  if (/\/(version|mirrors)\.json$/.test(url.pathname)) { try { return await fetch(req); } catch { return (await cache.match(req, { ignoreSearch: true })) || new Response('{}', { headers: { 'Content-Type': 'application/json' } }); } }
  const hit = await cache.match(key, { ignoreSearch: true });
  const refresh = fetch(nav ? 'index.html' : req, { cache: 'no-cache' }).then((r) => { if (r.ok && r.type === 'basic') cache.put(key, r.clone()); return r; });
  if (hit) { e.waitUntil(refresh.catch(() => {})); return hit; }
  try { return await refresh; } catch { return nav ? ((await cache.match('index.html')) || offlinePage()) : new Response('', { status: 504, statusText: 'Offline and not stored yet' }); }
}
const offlinePage = () => new Response('<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>AeroSuite 26</title><body style="font:16px system-ui;padding:2rem;max-width:40rem;margin:auto"><h1>Offline</h1><p>AeroSuite 26 has not been stored on this device yet. Connect once and open it again; after that it works without a connection.</p>', { headers: { 'Content-Type': 'text/html; charset=utf-8' } });

/** Live data: always try the network first; if unreachable, answer with the last stored reply, labelled. */
async function liveData(req) {
  const cache = await caches.open(LIVE_CACHE);
  try {
    const r = await fetch(req);
    if (r.ok) { const copy = r.clone(), h = new Headers(copy.headers); h.set('x-sw-stored-at', String(Date.now())); cache.put(req.url, new Response(await copy.blob(), { status: 200, headers: h })).then(() => trim(cache)); }
    return r;
  } catch (err) {
    const hit = await cache.match(req.url);
    if (!hit) throw err;
    const h = new Headers(hit.headers); h.set('x-sw-cached-at', hit.headers.get('x-sw-stored-at') || '0');
    return new Response(await hit.blob(), { status: 200, headers: h });
  }
}
async function trim(cache, max = 120) { const keys = await cache.keys(); for (let i = 0; i < keys.length - max; i++) await cache.delete(keys[i]); }

// Background refresh of stored live-data replies while the installed app is closed (where supported).
self.addEventListener('periodicsync', (e) => {
  if (e.tag !== 'refresh-live') return;
  e.waitUntil((async () => {
    const cache = await caches.open(LIVE_CACHE);
    for (const k of (await cache.keys()).slice(-40)) { try { const r = await fetch(k.url, { cache: 'no-store' }); if (r.ok) { const h = new Headers(r.headers); h.set('x-sw-stored-at', String(Date.now())); await cache.put(k.url, new Response(await r.blob(), { status: 200, headers: h })); } } catch { /* still offline */ } }
    (await self.clients.matchAll()).forEach((c) => c.postMessage({ type: 'periodic-refresh' }));
  })());
});
