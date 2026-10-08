// Update-prompt check: `node tests/update.mjs`
// Serves a copy of the built app, loads it, publishes a "new version" by re-stamping the service worker,
// and checks that the "Update now" prompt appears, that clicking it switches version, and that the
// prompt does not appear on a first visit.
import { spawn } from 'node:child_process';
import { cpSync, mkdtempSync, readFileSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright-core';

const root = fileURLToPath(new URL('..', import.meta.url)), dir = mkdtempSync(join(tmpdir(), 'aerosuite-upd-')), port = 8900 + Math.floor(Math.random() * 90), errors = [];
for (const f of ['index.html', 'sw.js', 'manifest.webmanifest', 'version.json', 'mirrors.json', 'css', 'js', 'assets', 'tools']) cpSync(root + f, join(dir, f), { recursive: true });
const server = spawn(process.execPath, [join(dir, 'tools/serve.mjs'), String(port)], { stdio: 'ignore' });
await new Promise((r) => setTimeout(r, 700));
const browser = await chromium.launch(), page = await (await browser.newContext()).newPage();
page.on('pageerror', (e) => errors.push('exception: ' + e.message));
try {
  await page.goto(`http://localhost:${port}/`); await page.waitForSelector('.suite-card');
  for (let k = 0; k < 60 && !(await page.evaluate(() => !!navigator.serviceWorker.controller)); k++) await page.waitForTimeout(250);
  await page.waitForTimeout(1500);
  if (await page.locator('.update-bar').count()) errors.push('update prompt shown on a first visit');
  // publish a new version
  const sw = readFileSync(join(dir, 'sw.js'), 'utf8'), v2 = 'test-' + Date.now();
  writeFileSync(join(dir, 'sw.js'), sw.replace(/const VERSION = [^;]+;/, `const VERSION = ${JSON.stringify(v2)};`));
  await page.evaluate(() => navigator.serviceWorker.getRegistration().then((r) => r.update()));
  try { await page.waitForSelector('.update-bar', { timeout: 30000 }); } catch { errors.push('no update prompt after a new version was published'); }
  if (!errors.length) {
    await Promise.all([page.waitForEvent('load', { timeout: 30000 }), page.locator('.update-bar .btn.primary').click()]);
    await page.waitForSelector('.suite-card'); await page.waitForTimeout(800);
    const caches2 = await page.evaluate(() => caches.keys());
    if (!caches2.some((c) => c.includes(v2))) errors.push('new version did not take over after clicking Update now: ' + caches2.join(','));
    if (await page.locator('.update-bar').count()) errors.push('prompt still shown after updating');
  }
} catch (e) { errors.push('aborted: ' + e.message.split('\n')[0]); }
await browser.close(); server.kill(); rmSync(dir, { recursive: true, force: true });
console.log(errors.length ? errors.join('\n') : 'Update prompt check passed'); process.exit(errors.length ? 1 : 0);
