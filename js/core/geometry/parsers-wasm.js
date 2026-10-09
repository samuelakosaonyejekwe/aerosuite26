// Lazy loaders for the vendored WebAssembly kernels (js/vendor/…). Nothing here is fetched or compiled until a
// file that needs it is opened. Each loader resolves URLs relative to this module, supplies the .wasm bytes
// itself (fetch in a browser or worker, the file system under Node) and caches the instantiated kernel.
// A page that uses them must allow WebAssembly compilation: CSP `script-src 'self' 'wasm-unsafe-eval'`.

const IS_NODE = typeof process === 'object' && !!process?.versions?.node && typeof window === 'undefined' && typeof importScripts === 'undefined';
const url = (rel) => new URL(rel, import.meta.url);
const cache = {};
/** Node built-ins are imported by computed name so that browser bundlers do not try to resolve them. */
export const nodeImport = (name) => import('node:' + name);

async function bytesOf(rel) {
  const u = url(rel);
  if (IS_NODE && u.protocol === 'file:') { const fs = await nodeImport('fs/promises'); return new Uint8Array(await fs.readFile(u)); }
  const r = await fetch(u);
  if (!r.ok) throw new Error(`could not fetch ${u.pathname.split('/').pop()} (HTTP ${r.status})`);
  return new Uint8Array(await r.arrayBuffer());
}
/** Run `make` once; a failed load is not cached, so a later attempt (e.g. back online) can succeed. */
function once(key, make) {
  if (!cache[key]) cache[key] = make().catch((e) => { cache[key] = null; throw new Error(`${key} kernel could not be loaded: ${e?.message || e}`); });
  return cache[key];
}

/** OpenCASCADE importer (occt-import-js): { ReadStepFile, ReadIgesFile, ReadBrepFile }. */
export const loadOcct = () => once('OpenCASCADE', async () => {
  const [mod, wasmBinary] = await Promise.all([import(url('../../vendor/occt/occt-import-js.js').href), bytesOf('../../vendor/occt/occt-import-js.wasm')]);
  return mod.default({ wasmBinary, print() {}, printErr() {} });
});

/** h5wasm (HDF5 C library): the high-level module namespace with File, FS … */
export const loadH5 = () => once('HDF5', async () => {
  const h5 = await import(url('../../vendor/h5wasm/hdf5_hl.js').href);
  const M = await h5.ready;
  if (typeof M.activate_throwing_error_handler === 'function') M.activate_throwing_error_handler();   // library errors become exceptions instead of console output
  return h5;
});

/** laz-perf (LASzip decoder): emscripten module with LASZip, _malloc, _free, HEAPU8. */
export const loadLazPerf = () => once('LASzip', async () => {
  const [mod, wasmBinary] = await Promise.all([import(url('../../vendor/lazperf/laz-perf.js').href), bytesOf('../../vendor/lazperf/laz-perf.wasm')]);
  if (!IS_NODE) return mod.default({ wasmBinary, print() {}, printErr() {} });
  // The emscripten glue expects CommonJS globals under Node; lend it `require` and `__dirname` while it starts.
  const { createRequire } = await nodeImport('module'), { fileURLToPath } = await nodeImport('url');
  const had = { r: Object.getOwnPropertyDescriptor(globalThis, 'require'), d: Object.getOwnPropertyDescriptor(globalThis, '__dirname') };
  globalThis.require = createRequire(import.meta.url); globalThis.__dirname = fileURLToPath(url('../../vendor/lazperf/'));
  try { return await mod.default({ wasmBinary, print() {}, printErr() {} }); }
  finally { for (const [k, d] of [['require', had.r], ['__dirname', had.d]]) { if (d) Object.defineProperty(globalThis, k, d); else delete globalThis[k]; } }
});

/** Forget an instantiated kernel after it has thrown from inside WebAssembly, so the next use starts from a clean instance. */
export function dropKernel(key) { cache[key] = null; }
