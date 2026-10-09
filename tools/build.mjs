// Release build: `node tools/build.mjs`
//  1. stamps the service worker with a content-derived version and the full offline asset list,
//  2. writes version.json,
//  3. adds configured mirror origins to the content-security policy (for mirror health checks),
//  4. bundles the whole application into one self-contained file, standalone.html.
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
const assets = ['index.html', 'manifest.webmanifest', 'mirrors.json', 'data/snapshot.json', ...(await walk('css')), ...(await walk('js')), ...(await walk('assets'))].filter((f) => !f.endsWith('.map')).sort();
// Files that change on a schedule without the application changing: stored offline, but left out of the version hash
// so that a fresh snapshot does not make every installed copy announce an update.
const VOLATILE = new Set(['data/snapshot.json']);

// ---- 3. CSP mirror origins (before hashing, since it changes index.html)
const mirrors = JSON.parse(await read('mirrors.json')).mirrors || [];
let html = await read('index.html');
const origins = [...new Set(mirrors.map((m) => new URL(m.url).origin))];
html = html.replace(/(connect-src 'self')((?: https:\/\/[^\s";]+)*)/, (_, a, rest) => { const keep = rest.trim().split(/\s+/).filter(Boolean); for (const o of origins) if (!keep.includes(o)) keep.push(o); return `${a} ${keep.join(' ')}`; });
await writeFile(join(root, 'index.html'), html);

// ---- 4. standalone single file
let standaloneOk = false;
try {
  const esbuild = await import('esbuild');
  const out = await esbuild.build({ entryPoints: [join(root, 'js/app.js')], bundle: true, format: 'iife', minify: true, write: false, target: ['es2020'], logLevel: 'silent', define: { 'import.meta.url': '""' }, legalComments: 'none' });
  const js = out.outputFiles[0].text, css = await read('css/app.css'), spec = (await read('js/data/spec.json')).replace(/</g, '\\u003c'), boot = await read('js/boot-theme.js');
  const icon = 'data:image/svg+xml;base64,' + Buffer.from(await read('assets/icon.svg')).toString('base64');
  const csp = html.match(/Content-Security-Policy" content="([^"]+)"/)[1].replace("script-src 'self'", "script-src 'unsafe-inline'").replace("style-src 'self'", "style-src 'unsafe-inline'");
  const page = `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>AeroSuite 26 — integrated aircraft engineering simulation</title>
<meta http-equiv="Content-Security-Policy" content="${csp}"><meta name="referrer" content="no-referrer"><meta name="theme-color" content="#1f6fd1"><meta name="color-scheme" content="light dark">
<link rel="icon" href="${icon}"><style>${css}</style><script>${boot}</script></head>
<body><div id="app" class="app"></div>
<script>globalThis.__AEROSUITE_STANDALONE__=true;globalThis.__AEROSUITE_MIRRORS__=${JSON.stringify(mirrors.map((m) => m.url)).replace(/</g, '\\u003c')};globalThis.__AEROSUITE_SPEC__=${spec};</script>
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
