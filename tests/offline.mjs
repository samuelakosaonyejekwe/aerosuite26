// Offline and single-file checks: `node tests/offline.mjs`
//  1. serve the built app, let the service worker store it, kill the network, reload, run an analysis;
//  2. open standalone.html straight from disk and run an analysis.
import { spawn } from 'node:child_process';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { chromium } from 'playwright-core';

const root = fileURLToPath(new URL('..', import.meta.url)), port = 8200 + Math.floor(Math.random() * 500), errors = [];
const server = spawn(process.execPath, [root + 'tools/serve.mjs', String(port)], { stdio: 'ignore' });
await new Promise((r) => setTimeout(r, 700));
const browser = await chromium.launch();
const watch = (page, tag) => { page.on('pageerror', (e) => errors.push(`[${tag}] exception: ${e.message}`)); page.on('console', (m) => { if (m.type() === 'error' && !/Failed to load resource|net::ERR|ERR_INTERNET_DISCONNECTED|Failed to fetch/.test(m.text())) errors.push(`[${tag}] console: ${m.text().slice(0, 200)}`); }); };
const runSuite = async (page, tag) => { await page.waitForSelector('.btn.primary.big', { timeout: 20000 }); await page.locator('.btn.primary.big').click(); try { await page.waitForSelector('.kpis .kpi', { timeout: 30000 }); } catch { errors.push(`[${tag}] analysis produced no results`); } };

try {
  // ---- 1. installed/offline behaviour
  const ctx = await browser.newContext(), page = await ctx.newPage(); watch(page, 'offline');
  await page.goto(`http://localhost:${port}/`); await page.waitForSelector('.suite-card');
  let files = 0;
  for (let k = 0; k < 120 && files < 60; k++) { await page.waitForTimeout(500); files = await page.evaluate(async () => { let n = 0; for (const c of await caches.keys()) if (c.startsWith('aerosuite-app-')) n += (await (await caches.open(c)).keys()).length; return n; }); }
  if (files < 60) errors.push(`[offline] only ${files} files were stored by the service worker`);
  await ctx.setOffline(true); server.kill();
  await page.reload(); await page.waitForSelector('.suite-card', { timeout: 20000 });
  await page.goto(`http://localhost:${port}/#/suite/economics`); await runSuite(page, 'offline');
  await page.goto(`http://localhost:${port}/#/suite/cfd/airfoil/spec`); try { await page.waitForSelector('.term', { timeout: 15000 }); } catch { errors.push('[offline] specification not available offline'); }
  await page.goto(`http://localhost:${port}/#/geometry`); try { await page.waitForSelector('.drop', { timeout: 15000 }); } catch { errors.push('[offline] geometry page not available offline'); }
  console.log(`offline: ${files} files stored, app reloaded and solved with the network and server off`);
  await ctx.close();
  // ---- 2. single-file copy from disk
  const ctx2 = await browser.newContext(), p2 = await ctx2.newPage(); watch(p2, 'standalone');
  await p2.goto(pathToFileURL(root + 'standalone.html').href); await p2.waitForSelector('.suite-card', { timeout: 20000 });
  await p2.goto(pathToFileURL(root + 'standalone.html').href + '#/suite/performance'); await runSuite(p2, 'standalone');
  await p2.goto(pathToFileURL(root + 'standalone.html').href + '#/suite/fea/beam/spec'); try { await p2.waitForSelector('.term', { timeout: 15000 }); } catch { errors.push('[standalone] specification missing'); }
  console.log('standalone: opened from disk and solved');
} catch (e) { errors.push('aborted: ' + e.message.split('\n')[0]); }
await browser.close(); server.kill();
console.log(errors.length ? errors.join('\n') : 'Offline checks passed'); process.exit(errors.length ? 1 : 0);
