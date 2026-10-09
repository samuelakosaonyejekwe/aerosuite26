// Release build: `node tools/build.mjs`
//  1. stamps the service worker with a content-derived version and the full offline asset list,
//  2. writes version.json,
//  3. sets the content-security policy's connect-src from the deployment configuration (config.json): configured
//     mirror origins and operator addresses are added, and with "commercial": true the hosts a commercial
//     deployment may not use are taken out, so the browser itself refuses them,
//  4. bundles the whole application into one self-contained file, standalone.html,
//  5. lists LICENSE, THIRD_PARTY_NOTICES.md and config.json in the offline set, so the notices ship with the site.
// The app itself has no runtime dependencies; esbuild is used here only to produce the single-file copy.

import { readFile, writeFile, readdir } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { join, relative, resolve, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = resolve(fileURLToPath(new URL('..', import.meta.url)));
const read = (p) => readFile(join(root, p), 'utf8');
async function walk(dir, out = []) { for (const e of await readdir(join(root, dir), { withFileTypes: true })) { const p = join(dir, e.name); if (e.isDirectory()) await walk(p, out); else out.push(p.split(sep).join('/')); } return out; }

const pkg = JSON.parse(await read('package.json'));
// js/ includes the bundled airport database (js/data/airports); data/snapshot.json is the cloud snapshot of live feeds
// The licence notices and the deployment configuration are part of the site: LICENSE, THIRD_PARTY_NOTICES.md (and the
// vendored licence texts under js/vendor/, picked up by the walk) must reach every recipient of the application.
const NOTICES = ['LICENSE', 'THIRD_PARTY_NOTICES.md', 'COMPLIANCE.md', 'config.json'];
for (const f of NOTICES) await readFile(join(root, f)).catch(() => { throw new Error(`${f} is missing: it must ship with the site`); });
for (const f of ['js/vendor/occt/LICENSE.occt-import-js.txt', 'js/vendor/occt/LICENSE.occt.txt', 'js/vendor/h5wasm/LICENSE.txt', 'js/vendor/lazperf/LICENSE.txt', 'js/data/licences.json']) await readFile(join(root, f)).catch(() => { throw new Error(`${f} is missing: third-party licence texts must be deployed`); });
const gridFiles = await walk('data/grid').catch(() => []); // forecast grids: rebuilt with the snapshot
const assets = ['index.html', 'manifest.webmanifest', 'mirrors.json', 'data/snapshot.json', ...gridFiles, ...NOTICES, ...(await walk('css')), ...(await walk('js')), ...(await walk('assets'))].filter((f) => !f.endsWith('.map')).sort();
// Files that change on a schedule or per deployment without the application changing: stored offline, but left out of
// the version hash so that a fresh snapshot or a configuration change does not make every installed copy announce an update.
const VOLATILE = new Set(['data/snapshot.json', 'config.json', ...gridFiles]);

// ---- 3. CSP connect-src (before hashing, since it changes index.html)
// Every host a connector can use, with the licence-registry id that governs it (js/data/licences.json). Hosts of a
// source that a commercial deployment may not use are left out of the policy when config.json says "commercial": true
// (unless the operator listed the id under "accept"); customer hosts are listed only when they can be used.
const LIVE_HOSTS = [
  ['open-meteo-free', ['api.open-meteo.com', 'geocoding-api.open-meteo.com', 'archive-api.open-meteo.com', 'air-quality-api.open-meteo.com', 'marine-api.open-meteo.com']],
  ['open-meteo-customer', ['customer-api.open-meteo.com', 'customer-geocoding-api.open-meteo.com', 'customer-archive-api.open-meteo.com', 'customer-air-quality-api.open-meteo.com', 'customer-marine-api.open-meteo.com']],
  ['met-norway', ['api.met.no']], ['nws', ['api.weather.gov']],
  ['overpass-public', ['overpass-api.de', 'overpass.private.coffee', 'maps.mail.ru']],
  ['noaa-swpc', ['services.swpc.noaa.gov']], ['ecb-statistics', ['data-api.ecb.europa.eu']], ['frankfurter', ['api.frankfurter.dev']], ['world-bank', ['api.worldbank.org']], ['uk-ets-price', ['www.gov.uk']],
  ['oil-prices-dataset', ['raw.githubusercontent.com']], ['gb-carbon-intensity', ['api.carbonintensity.org.uk']], ['openalex', ['api.openalex.org']],
  ['github-user-bridge', ['api.github.com', '*.blob.core.windows.net']], // solver bridge: the user's own token; artifact downloads are redirected to blob storage
];
const { cleanConfig } = await import('../js/core/live.js');
let rawConfig = {}; try { rawConfig = JSON.parse(await read('config.json')); } catch (e) { throw new Error(`config.json is not valid JSON: ${e.message}`); }
const cfg = cleanConfig(rawConfig), registry = JSON.parse(await read('js/data/licences.json'));
const hostAllowed = (id) => {
  if (id === 'github-user-bridge') return true;
  if (id === 'open-meteo-customer') return !cfg.commercial || !!cfg.openMeteoApiKey; // listed in a non-commercial build so that adding a key later needs no rebuild
  if (!cfg.commercial || cfg.accept.includes(id)) return true;
  return registry.sources.find((x) => x.id === id)?.class === 'commercial-ok';
};
const mirrors = JSON.parse(await read('mirrors.json')).mirrors || [];
let html = await read('index.html');
const origins = [...LIVE_HOSTS.filter(([id]) => hostAllowed(id)).flatMap(([, hs]) => hs.map((x) => `https://${x}`)), ...new Set(mirrors.map((m) => new URL(m.url).origin))];
for (const u of [rawConfig.carbonPriceUrl, rawConfig.overpassUrl]) { try { const o = new URL(u).origin; if (/^https:/.test(o) && !origins.includes(o)) origins.push(o); } catch { /* empty, or a path on this site */ } }
if (!/connect-src 'self'[^;"]*/.test(html)) throw new Error('index.html: connect-src not found in the content-security policy');
html = html.replace(/connect-src 'self'[^;"]*/, `connect-src 'self' ${origins.join(' ')}`);
await writeFile(join(root, 'index.html'), html);
console.log(`content-security policy: ${origins.length} hosts (${cfg.commercial ? 'commercial' : 'non-commercial'} deployment)`);

// ---- essential data for the single-file copy -------------------------------------------------------
// standalone.html is opened from disk and cannot read the data files next to the site. It therefore carries a compact
// subset inline: large and medium airports with their runways (OurAirports), cities of 100 000 inhabitants or more
// (GeoNames), the design-temperature grid (NCEP-DOE Reanalysis 2) and the latest cloud snapshot (carbon panel, exchange
// rates, fuel prices, policy rates). js/core/airports.js, places.js and live.js fall back to it when there is no site.
async function essentialData() {
  const files = {}, dir = 'js/data/airports/', idx = JSON.parse(await read(dir + 'index.json')), tiles = {}, codes = []; let nAir = 0;
  for (const id of Object.keys(idx.tiles)) {
    const out = [], cs = new Set(); let keep = false, head = null, n = 0;
    for (const l of (await read(`${dir}t${id}.txt`)).split('\n')) {
      if (!l) continue; if (l[0] === '#') { keep = /\|[LM]$/.test(l); head = l; continue; } if (!keep) continue;
      if (head) { out.push(head); head = null; } out.push(l); n++; const f = l.split('|'); for (const c of [f[0], f[1] === '1' ? '' : f[1], f[2]]) if (c) cs.add(c);
    }
    if (n) { files[`airports/t${id}.txt`] = out.join('\n') + '\n'; tiles[id] = n; codes.push(`${id}:${[...cs].join(' ')}`); nAir += n; }
  }
  const names = []; let on = false; for (const l of (await read(dir + 'names.txt')).split('\n')) { if (l[0] === '#') { on = l === '#L' || l === '#M'; if (on) names.push(l); continue; } if (on && l) names.push(l); }
  files['airports/codes.txt'] = codes.join('\n') + '\n'; files['airports/names.txt'] = names.join('\n') + '\n';
  files['airports/index.json'] = JSON.stringify({ ...idx, tiles, counts: { ...idx.counts, airports: nAir }, subset: 'large and medium airports only (single-file copy)' });
  const big = []; let head = null, nCity = 0; for (const l of (await read('js/data/places/cities.txt')).split('\n')) { if (!l) continue; if (l[0] === '#') { head = l; continue; } if (Number(l.split('|')[5]) >= 100) { if (head) { big.push(head); head = null; } big.push(l); nCity++; } }
  const pm = JSON.parse(await read('js/data/places/index.json')); delete pm.small;
  files['places/cities.txt'] = big.join('\n') + '\n'; files['places/index.json'] = JSON.stringify({ ...pm, counts: { places: nCity }, subset: 'cities of 100 000 inhabitants or more (single-file copy)' });
  files['climate/index.json'] = await read('js/data/climate/index.json'); files['climate/grid.bin'] = { b64: (await readFile(join(root, 'js/data/climate/grid.bin'))).toString('base64') };
  let snapshot = null; try { snapshot = JSON.parse(await read('data/snapshot.json')); } catch { /* none yet */ }
  return { built: new Date().toISOString(), carries: { airports: nAir, cities: nCity, climate: true, snapshot: snapshot?.generated || null }, files, snapshot };
}

// ---- 4. standalone single file
let standaloneOk = false;
try {
  const esbuild = await import('esbuild');
  const out = await esbuild.build({ entryPoints: [join(root, 'js/app.js')], bundle: true, format: 'iife', minify: true, write: false, target: ['es2020'], logLevel: 'silent', define: { 'import.meta.url': '""' }, legalComments: 'none' });
  const essential = await essentialData(); console.log(`single-file copy carries ${essential.carries.airports} airports, ${essential.carries.cities} cities, the climate grid and the snapshot (${(JSON.stringify(essential).length / 1e6).toFixed(2)} MB)`);
  const js = out.outputFiles[0].text, css = await read('css/app.css'), spec = (await read('js/data/spec.json')).replace(/</g, '\\u003c'), boot = await read('js/boot-theme.js'), inline = (o) => JSON.stringify(o).replace(/</g, '\\u003c');
  const icon = 'data:image/svg+xml;base64,' + Buffer.from(await read('assets/icon.svg')).toString('base64');
  const csp = html.match(/Content-Security-Policy" content="([^"]+)"/)[1].replace("script-src 'self'", "script-src 'unsafe-inline'").replace("style-src 'self'", "style-src 'unsafe-inline'");
  const page = `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>AeroSuite 26 — integrated aircraft engineering simulation</title>
<meta http-equiv="Content-Security-Policy" content="${csp}"><meta name="referrer" content="no-referrer"><meta name="theme-color" content="#1f6fd1"><meta name="color-scheme" content="light dark">
<link rel="icon" href="${icon}"><style>${css}</style><script>${boot}</script></head>
<body><div id="app" class="app"></div>
<script>globalThis.__AEROSUITE_STANDALONE__=true;globalThis.__AEROSUITE_MIRRORS__=${JSON.stringify(mirrors.map((m) => m.url)).replace(/</g, '\\u003c')};globalThis.__AEROSUITE_SPEC__=${spec};globalThis.__AEROSUITE_CONFIG__=${inline(rawConfig)};globalThis.__AEROSUITE_LICENCES__=${inline(registry)};globalThis.__AEROSUITE_GRID_FACTORS__=${inline(JSON.parse(await read('js/data/grid-factors.json')))};globalThis.__AEROSUITE_TOOLS__=${inline(JSON.parse(await read('js/data/tools.json')))};globalThis.__AEROSUITE_ESSENTIAL__=${inline(essential)};</script>
<script>${js.replace(/<\/script/gi, '<\\/script')}</script></body></html>
`;
  await writeFile(join(root, 'standalone.html'), page);
  standaloneOk = true;
  console.log(`standalone.html  ${(page.length / 1e6).toFixed(2)} MB`);
} catch (e) { console.warn('standalone.html was not rebuilt (' + (e.code === 'ERR_MODULE_NOT_FOUND' ? 'run `npm install` first to get esbuild' : e.message) + ')'); }

// ---- 1. + 2. version, asset list
const list = [...assets, 'standalone.html'];
const hash = createHash('sha256');
for (const f of list) { if (VOLATILE.has(f)) continue; try { hash.update(f); hash.update(await readFile(join(root, f))); } catch { /* optional file */ } }
const version = `${pkg.version}+${hash.digest('hex').slice(0, 10)}`, built = new Date().toISOString();
let sw = await read('sw.js');
sw = sw.replace(/\/\* BUILD:BEGIN \*\/[\s\S]*?\/\* BUILD:END \*\//, `/* BUILD:BEGIN */\nconst VERSION = ${JSON.stringify(version)};\nconst ASSETS = ${JSON.stringify(['./', ...list])};\n/* BUILD:END */`);
await writeFile(join(root, 'sw.js'), sw);
await writeFile(join(root, 'version.json'), JSON.stringify({ name: 'AeroSuite 26', version, built, files: list.length }, null, 1) + '\n');
console.log(`version ${version}  ·  ${list.length} files in the offline set${standaloneOk ? '' : '  ·  (standalone not rebuilt)'}`);
