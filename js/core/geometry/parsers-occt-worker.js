// Worker entry for the OpenCASCADE kernel. Exact-CAD translation runs here, off the caller's thread, so that a
// file which drives the kernel into a very long (or endless) computation can be stopped by terminating the worker.
// Works as a browser module worker and as a Node worker thread.

import { loadOcct, nodeImport } from './parsers-wasm.js';

const node = typeof process === 'object' && !!process?.versions?.node && typeof self === 'undefined';
const port = node ? (await nodeImport('worker_threads')).parentPort : self;

async function run({ id, kind, bytes, params }) {
  try {
    const occt = await loadOcct();
    const res = kind === 'step' ? occt.ReadStepFile(bytes, params) : kind === 'iges' ? occt.ReadIgesFile(bytes, params) : occt.ReadBrepFile(bytes, params);
    if (!res || !res.success) return { reply: { id, ok: false, error: 'the OpenCASCADE kernel could not translate this file' }, transfer: [] };
    const transfer = [], meshes = (res.meshes || []).map((m) => {
      const position = Float64Array.from(m.attributes?.position?.array || []), index = Uint32Array.from(m.index?.array || []);
      transfer.push(position.buffer, index.buffer);
      return { name: m.name ?? '', color: m.color ?? null, brep_faces: m.brep_faces || [], position, index };
    });
    return { reply: { id, ok: true, root: res.root, meshes }, transfer };
  } catch (e) {
    return { reply: { id, ok: false, fatal: true, error: typeof e === 'number' || /exception/i.test(String(e?.message)) ? 'the OpenCASCADE kernel raised an internal exception on this file' : String(e?.message || e) }, transfer: [] };
  }
}
const onMessage = async (msg) => { const { reply, transfer } = await run(msg); port.postMessage(reply, transfer); };
if (node) port.on('message', onMessage); else port.onmessage = (e) => onMessage(e.data);
