// Runs jobs in a module Web Worker when available, otherwise on the main thread (e.g. the
// single-file standalone build opened from disk). Supports progress reporting and cancellation.

import { execute } from './jobs.js';

let worker = null, seq = 0, workerBroken = typeof Worker === 'undefined' || globalThis.__AEROSUITE_STANDALONE__ === true;
const pending = new Map();

function spawn() {
  try {
    worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = (ev) => {
      const { id, progress, msg, ok, result, error } = ev.data, p = pending.get(id);
      if (!p) return;
      if (progress !== undefined) { p.onProgress?.(progress, msg); return; }
      pending.delete(id);
      ok ? p.resolve(result) : p.reject(new Error(error));
    };
    worker.onerror = () => { // module workers unsupported or script failed to load: fall back for everything
      workerBroken = true; worker = null;
      for (const [id, p] of pending) { pending.delete(id); execute(p.job, p.onProgress).then(p.resolve, p.reject); }
    };
  } catch { workerBroken = true; worker = null; }
}

export function runJob(job, onProgress) {
  if (workerBroken) return yieldThen(() => execute(job, onProgress));
  if (!worker) spawn();
  if (!worker) return yieldThen(() => execute(job, onProgress));
  const id = ++seq;
  return new Promise((resolve, reject) => { pending.set(id, { resolve, reject, onProgress, job }); worker.postMessage({ id, job }); });
}
const yieldThen = (fn) => new Promise((res, rej) => setTimeout(() => fn().then(res, rej), 20));

/** Abort everything in flight by restarting the worker. */
export function cancelAll() {
  if (worker) { worker.terminate(); worker = null; }
  for (const [id, p] of pending) { pending.delete(id); p.reject(new Error('Cancelled')); }
}
export const usingWorker = () => !workerBroken;
