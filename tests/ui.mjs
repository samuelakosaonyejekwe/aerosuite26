// Browser smoke test: `node tests/ui.mjs [--shots <dir>] [--mobile]`
// Serves the app, walks every page in headless Chromium, runs analyses and studies through the real
// UI, imports a generated geometry file, and fails on any console error or uncaught exception.

import { spawn } from 'node:child_process';
import { mkdirSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright-core';
import { SUITES } from '../js/core/registry.js';

const args = process.argv.slice(2), shotDir = args.includes('--shots') ? args[args.indexOf('--shots') + 1] : null, mobile = args.includes('--mobile');
const only = args.includes('--only') ? args[args.indexOf('--only') + 1].split(',') : null;
const port = 8000 + Math.floor(Math.random() * 900), base = `http://localhost:${port}/`;
const server = spawn(process.execPath, [fileURLToPath(new URL('../tools/serve.mjs', import.meta.url)), String(port)], { stdio: 'ignore' });
await new Promise((r) => setTimeout(r, 700));
if (shotDir) mkdirSync(shotDir, { recursive: true });

const browser = await chromium.launch();
const ctx = await browser.newContext(mobile ? { viewport: { width: 390, height: 800 }, deviceScaleFactor: 2, isMobile: true, hasTouch: true } : { viewport: { width: 1360, height: 900 } });
const page = await ctx.newPage();
const errors = []; let step = 'boot';
page.on('console', (m) => { if (m.type() === 'error' && !/Failed to load resource|net::ERR|status of 4\d\d|status of 5\d\d/.test(m.text())) errors.push(`[${step}] console: ${m.text()}`); });
page.on('pageerror', (e) => errors.push(`[${step}] exception: ${e.message}`));
const shot = async (name) => { if (shotDir) await page.screenshot({ path: join(shotDir, `${mobile ? 'm-' : ''}${name}.png`), fullPage: false }); };
const go = async (hash, name) => { step = name; await page.goto(base + hash); await page.waitForSelector(`#page[data-route="${hash}"]`, { timeout: 60000 }); await page.waitForTimeout(150); };
const want = (name) => !only || only.includes(name);
const fail = (m) => errors.push(`[${step}] ${m}`);

try {
  if (want('home')) { await go('#/home', 'home'); if (!(await page.locator('.suite-card').count())) fail('no suite cards'); if ((await page.locator('.suite-card').count()) !== 26) fail('expected 26 suite cards'); await shot('home'); }
  if (want('case')) {
    await go('#/case', 'case'); await page.locator('.suite-card', { hasText: 'Regional twin turboprop' }).click(); await page.waitForTimeout(300); await shot('case');
    await page.locator('.tab', { hasText: 'Location' }).click(); await page.waitForTimeout(300); await shot('case-site');
    await page.locator('.tab', { hasText: 'Import' }).click(); await page.waitForTimeout(200);
    await go('#/home', 'case-reset'); await go('#/case', 'case-reset'); await page.locator('.suite-card', { hasText: 'Narrow-body' }).click();
  }
  if (want('suites')) for (const s of SUITES) {
    await go(`#/suite/${s.id}`, `suite ${s.id}`);
    if (await page.locator('.note.bad').count()) { fail(`suite page error: ${await page.locator('.note.bad').first().innerText()}`); continue; }
    const chips = await page.locator('.chip:not(.na)').count(); if (!chips) { if (!(await page.locator('.note.warn', { hasText: 'Not applicable' }).count())) fail('no applicable analyses and no explanation shown'); continue; } // e.g. rotor suites on a jet airliner
    await page.locator('.btn.primary.big').click();
    try { await page.waitForSelector('.kpis .kpi', { timeout: 20000 }); } catch { fail('run produced no KPIs'); }
    await page.waitForTimeout(250);
    if (!(await page.locator('.plot canvas').count()) && !(await page.locator('table.data').count())) fail('no plots or tables rendered');
    if (s.id === 'performance' || s.id === 'cfd' || s.id === 'fea' || s.id === 'economics') await shot(`suite-${s.id}`);
  }
  if (want('studies')) {
    await go('#/suite/performance/field/mesh', 'mesh tab'); await page.locator('.btn.primary', { hasText: 'Run convergence study' }).click(); await page.waitForSelector('.kpis .kpi', { timeout: 30000 }); await page.waitForTimeout(300); await shot('mesh');
    await go('#/suite/performance/range/studies', 'studies tab'); await page.locator('.btn.primary', { hasText: 'Rank inputs' }).click(); await page.waitForSelector('.plot canvas', { timeout: 30000 });
    await page.locator('.btn.primary', { hasText: 'Run sweep' }).click(); await page.locator('.btn.primary', { hasText: 'Propagate uncertainty' }).click(); for (let k = 0; k < 160 && (await page.locator('.plot canvas').count()) < 4; k++) await page.waitForTimeout(250); if ((await page.locator('.plot canvas').count()) < 4) fail('studies did not render four charts'); await shot('studies');
    await go('#/suite/performance/point/vv', 'vv tab'); await page.locator('.btn.primary', { hasText: 'Run verification' }).click(); await page.waitForSelector('table.data', { timeout: 20000 });
    await page.locator('textarea.inp').fill('60000, 150\n70000, 162\n78000, 171'); await page.locator('.btn.primary', { hasText: 'Validate' }).click(); await page.waitForSelector('.plot canvas', { timeout: 20000 }); await shot('vv');
    await go('#/suite/performance/point/spec', 'spec tab'); await page.waitForSelector('.term', { timeout: 10000 }); await shot('spec');
    await go('#/suite/performance/point/live', 'live tab'); await page.waitForTimeout(600);
  }
  if (want('geometry')) {
    await go('#/geometry', 'geometry');
    // a closed NACA-like extruded wing as ASCII STL
    const N = 40, span = 5, prof = []; for (let i = 0; i <= N; i++) { const x = 0.5 * (1 - Math.cos((Math.PI * i) / N)), t = 0.6 * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1036 * x ** 4); prof.push([x, t]); }
    const loop = [...prof, ...prof.slice(1, -1).reverse().map(([x, t]) => [x, -t])], tri = [], f = (a, b, c) => tri.push(`facet normal 0 0 0\nouter loop\nvertex ${a.join(' ')}\nvertex ${b.join(' ')}\nvertex ${c.join(' ')}\nendloop\nendfacet`);
    for (let i = 0; i < loop.length; i++) { const p = loop[i], q = loop[(i + 1) % loop.length]; f([p[0], 0, p[1]], [q[0], 0, q[1]], [q[0], span, q[1]]); f([p[0], 0, p[1]], [q[0], span, q[1]], [p[0], span, p[1]]); }
    for (let i = 1; i < loop.length - 1; i++) { f([loop[0][0], 0, loop[0][1]], [loop[i + 1][0], 0, loop[i + 1][1]], [loop[i][0], 0, loop[i][1]]); f([loop[0][0], span, loop[0][1]], [loop[i][0], span, loop[i][1]], [loop[i + 1][0], span, loop[i + 1][1]]); }
    const file = join(tmpdir(), `aerosuite-wing-${process.pid}.stl`); writeFileSync(file, `solid wing\n${tri.join('\n')}\nendsolid wing\n`);
    const [chooser] = await Promise.all([page.waitForEvent('filechooser'), page.locator('.drop').click()]); await chooser.setFiles(file);
    await page.waitForSelector('.viewer canvas', { timeout: 20000 }); await page.waitForTimeout(500);
    await page.locator('.btn', { hasText: 'Confirm' }).click(); await page.waitForTimeout(500); await shot('geometry');
    await page.locator('.btn.primary', { hasText: 'Cut section' }).click(); await page.waitForSelector('.plot canvas', { timeout: 10000 });
    await page.locator('.btn', { hasText: 'with mesh-sensitivity study' }).click(); await page.waitForSelector('.plot-tri canvas', { timeout: 30000 }); await page.waitForTimeout(300);
    await page.locator('.plot-tri').scrollIntoViewIfNeeded(); await shot('geometry-section');
  }
  if (want('cad')) { // exact CAD through the embedded geometry kernel
    await go('#/home', 'cad'); await go('#/geometry', 'cad');
    const [ch] = await Promise.all([page.waitForEvent('filechooser'), page.locator('.drop').click()]); await ch.setFiles(fileURLToPath(new URL('./fixtures/box.stp', import.meta.url)));
    let txt = ''; for (let k = 0; k < 480; k++) { txt = await page.locator('#page').innerText(); if (/box\.stp/.test(txt) && /12 triangles/.test(txt)) break; await page.waitForTimeout(250); }
    if (!(await page.locator('.viewer canvas').count())) fail('STEP file produced no 3-D view');
    if (!/Read in full/.test(txt)) fail('STEP not reported as read in full');
    if (!/12 triangles/.test(txt)) fail('STEP faces were not tessellated in the browser (expected 12 triangles for a box)'); await shot('geometry-step');
  }
  if (want('integrated')) {
    await go('#/integrated', 'integrated'); await page.waitForSelector('.flowmap .nd', { timeout: 30000 }); await page.locator('.btn.primary.big').click();
    await page.waitForSelector('table.data', { timeout: 180000 }); await page.waitForTimeout(400); await shot('integrated');
    const failed = await page.locator('td', { hasText: /^Failed:/ }).count(); if (failed) fail(`${failed} analyses failed in the integrated run: ${await page.locator('td', { hasText: /^Failed:/ }).first().innerText()}`);
    await go('#/decisions', 'decisions'); await page.waitForTimeout(400); if (!(await page.locator('.rec').count())) fail('no recommendations after integrated run'); await shot('decisions');
  }
  if (want('pages')) {
    await go('#/live', 'live'); await page.waitForTimeout(500); await shot('live');
    await go('#/reports', 'reports'); await page.locator('.btn.primary', { hasText: 'Verify all suites' }).click(); await page.waitForSelector('table.data', { timeout: 120000 }); const failedV = await page.locator('td', { hasText: /^FAIL$/ }).count(); if (failedV) fail(`${failedV} verification benchmarks fail in the browser`); await shot('reports');
    await page.locator('.btn', { hasText: 'Build printable report' }).click(); await page.waitForTimeout(1500);
    await go('#/about', 'about'); await page.waitForSelector('table.data', { timeout: 60000 }); await shot('about');
    // navigation arrows and palette
    await go('#/home', 'nav'); await page.locator('.pager a.next').click(); await page.waitForTimeout(300); if (!page.url().includes('#/case')) fail('pager next did not go to the case page');
    await page.locator('.top .arrow').first().click(); await page.waitForTimeout(300); if (!page.url().includes('#/home')) fail('back arrow did not return home');
    await page.keyboard.press('Control+k'); await page.locator('.palette input').fill('flutter'); await page.waitForTimeout(200); if (!(await page.locator('.palette li').count())) fail('palette empty'); await page.keyboard.press('Escape');
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth - window.innerWidth); if (overflow > 2) fail(`page overflows horizontally by ${overflow}px`);
  }
} catch (e) { errors.push(`[${step}] test aborted: ${e.message.split('\n')[0]}`); await shot('failure'); }

await browser.close(); server.kill();
console.log(errors.length ? errors.join('\n') : 'UI smoke test passed');
console.log(`${errors.length} problem(s)`);
process.exit(errors.length ? 1 : 0);
