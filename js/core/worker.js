// Web Worker entry: keeps heavy computation off the UI thread.
import { execute } from './jobs.js';

self.onmessage = async (ev) => {
  const { id, job } = ev.data;
  let last = 0;
  const progress = (f, msg) => { const now = Date.now(); if (now - last > 60 || f >= 1) { last = now; self.postMessage({ id, progress: f, msg }); } };
  try { self.postMessage({ id, ok: true, result: await execute(job, progress) }); }
  catch (e) { self.postMessage({ id, ok: false, error: e?.message || String(e) }); }
};
