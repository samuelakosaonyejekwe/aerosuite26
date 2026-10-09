// Runs jobs in a module Web Worker when available, otherwise on the main thread (e.g. the
// single-file standalone build opened from disk). Supports progress reporting and cancellation.

import { execute } from './jobs.js';
import { customMaterials } from './matlib.js';

const withMaterials = (job) => (job.materials ? job : { ...job, materials: customMaterials() });

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
  job = withMaterials(job);
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

// ---- parallel pool: independent jobs on several cores at once ---------------------------------
/**
 * Run jobs concurrently on up to `size` workers. onProgress(index, fraction, message) reports per job.
 * Falls back to sequential execution where module workers are unavailable. Resolves to results in order.
 */
export async function runJobsParallel(jobs, onProgress = () => {}, size = Math.max(1, Math.min(6, (globalThis.navigator?.hardwareConcurrency || 4) - 1))) {
  jobs = jobs.map(withMaterials);
  if (workerBroken || jobs.length < 2 || size < 2) { const out = []; for (let i = 0; i < jobs.length; i++) out.push(await runJob(jobs[i], (f, m) => onProgress(i, f, m))); return out; }
  const results = new Array(jobs.length); let next = 0, failed = null;
  const lane = async () => {
    let w; try { w = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' }); } catch { w = null; }
    while (next < jobs.length && !failed) {
      const i = next++;
      try {
        results[i] = w ? await new Promise((res, rej) => { w.onmessage = (ev) => { const d = ev.data; if (d.progress !== undefined) onProgress(i, d.progress, d.msg); else d.ok ? res(d.result) : rej(new Error(d.error)); }; w.onerror = (e) => rej(new Error(e.message || 'Worker failed')); w.postMessage({ id: i, job: jobs[i] }); })
          : await execute(jobs[i], (f, m) => onProgress(i, f, m));
      } catch (e) { failed = e; }
    }
    w?.terminate(); pool.delete(w);
  };
  const lanes = []; for (let k = 0; k < Math.min(size, jobs.length); k++) lanes.push(lane());
  await Promise.all(lanes);
  if (failed) throw failed;
  return results;
}
const pool = new Set();
