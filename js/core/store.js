// Project state and persistence. Small state lives in localStorage; bulky items (results, uploaded
// geometry, cached live data) live in IndexedDB. Everything degrades to in-memory if storage is blocked.

import { defaultCase, normaliseCase } from './case.js';

const LS = 'aerosuite26.';
const mem = new Map();
const ls = {
  get(k, d) { try { const v = localStorage.getItem(LS + k); return v == null ? d : JSON.parse(v); } catch { return mem.has(k) ? mem.get(k) : d; } },
  set(k, v) { mem.set(k, v); try { localStorage.setItem(LS + k, JSON.stringify(v)); } catch { /* quota or private mode */ } },
};

// ---- IndexedDB key-value -----------------------------------------------------------------
let dbp;
function db() {
  if (!dbp) dbp = new Promise((res) => {
    try {
      const rq = indexedDB.open('aerosuite26', 1);
      rq.onupgradeneeded = () => rq.result.createObjectStore('kv');
      rq.onsuccess = () => res(rq.result); rq.onerror = () => res(null); rq.onblocked = () => res(null);
    } catch { res(null); }
  });
  return dbp;
}
const idbMem = new Map();
export const idb = {
  async get(k) { const d = await db(); if (!d) return idbMem.get(k); return new Promise((res) => { const r = d.transaction('kv').objectStore('kv').get(k); r.onsuccess = () => res(r.result); r.onerror = () => res(undefined); }); },
  async set(k, v) { const d = await db(); if (!d) { idbMem.set(k, v); return; } return new Promise((res) => { const t = d.transaction('kv', 'readwrite'); t.objectStore('kv').put(v, k); t.oncomplete = () => res(); t.onerror = () => res(); t.onabort = () => res(); }); },
  async del(k) { const d = await db(); if (!d) { idbMem.delete(k); return; } return new Promise((res) => { const t = d.transaction('kv', 'readwrite'); t.objectStore('kv').delete(k); t.oncomplete = () => res(); t.onerror = () => res(); }); },
  async keys(prefix = '') { const d = await db(); if (!d) return [...idbMem.keys()].filter((k) => k.startsWith(prefix)); return new Promise((res) => { const r = d.transaction('kv').objectStore('kv').getAllKeys(); r.onsuccess = () => res(r.result.filter((k) => String(k).startsWith(prefix))); r.onerror = () => res([]); }); },
};

// ---- events ------------------------------------------------------------------------------
const subs = new Map();
export function on(evt, fn) { if (!subs.has(evt)) subs.set(evt, new Set()); subs.get(evt).add(fn); return () => subs.get(evt)?.delete(fn); }
export function emit(evt, data) { for (const fn of subs.get(evt) || []) { try { fn(data); } catch (e) { console.error(e); } } }

// ---- state -------------------------------------------------------------------------------
export const state = {
  case: normaliseCase(ls.get('case', null) || defaultCase()),
  preset: ls.get('preset', 'narrowbody'),
  overrides: ls.get('overrides', {}),     // { 'suite.analysis': { key: value } } — user-entered inputs
  up: ls.get('up', {}),                   // published outputs per suite (the coupling bus)
  runs: ls.get('runs', {}),               // { 'suite.analysis': { ts, ms, status, nWarn } }
  recs: ls.get('recs', {}),               // { 'suite.analysis': [recommendation] }
  settings: ls.get('settings', { theme: 'auto', autoApplyLive: true, autoRefreshLive: true, autoRun: false, currency: 'USD', fuelCrack: 0.2 }),
  results: {},                            // full results, lazily loaded from IndexedDB
  geometry: [],                           // imported models (session); summaries persisted
};
const saveCase = () => ls.set('case', state.case);

export function setCase(c, preset) { state.case = normaliseCase(c); if (preset !== undefined) { state.preset = preset; ls.set('preset', preset); } saveCase(); emit('case', state.case); }
export function patchCase(section, key, value) { state.case[section][key] = value; saveCase(); emit('case', state.case); }
export function patchCaseMany(section, obj, silent = false) { Object.assign(state.case[section], obj); saveCase(); if (!silent) emit('case', state.case); }
export function setSetting(k, v) { state.settings[k] = v; ls.set('settings', state.settings); emit('settings', state.settings); }

const key = (s, a) => `${s}.${a}`;
export const getOverrides = (s, a) => state.overrides[key(s, a)] || {};
export function setOverride(s, a, k, v) { const o = (state.overrides[key(s, a)] ||= {}); if (v === undefined) delete o[k]; else o[k] = v; ls.set('overrides', state.overrides); }
export function clearOverrides(s, a) { delete state.overrides[key(s, a)]; ls.set('overrides', state.overrides); }

/** Store a finished run: full result to IndexedDB, summary + published outputs to localStorage. */
export async function saveRun(s, a, payload) {
  const k = key(s, a);
  state.results[k] = payload;
  const publish = Object.fromEntries(Object.entries(payload.res.outputs || {}).filter(([, v]) => typeof v === 'number' || (Array.isArray(v) && JSON.stringify(v).length < 6000)));
  state.up[s] = { ...(state.up[s] || {}), ...publish };
  const bad = payload.res.kpis.some((x) => x.status === 'bad'), warn = payload.res.kpis.some((x) => x.status === 'warn') || payload.res.warnings.length > 0;
  state.runs[k] = { ts: Date.now(), ms: payload.res.elapsed_ms, status: bad ? 'bad' : warn ? 'warn' : 'ok', nWarn: payload.res.warnings.length };
  state.recs[k] = payload.recs || [];
  ls.set('up', state.up); ls.set('runs', state.runs); ls.set('recs', state.recs);
  idb.set('res.' + k, payload);
  emit('run', { suite: s, analysis: a });
}
export async function loadRun(s, a) { const k = key(s, a); if (!state.results[k]) { const p = await idb.get('res.' + k); if (p) state.results[k] = p; } return state.results[k]; }
export async function clearRuns() {
  state.results = {}; state.up = {}; state.runs = {}; state.recs = {};
  ls.set('up', {}); ls.set('runs', {}); ls.set('recs', {});
  for (const k of await idb.keys('res.')) await idb.del(k);
  emit('run', {});
}

/** Whole-project export/import (case, inputs, published outputs, run summaries, recommendations). */
export function exportProject() {
  return { app: 'AeroSuite 26', schema: 1, exported: new Date().toISOString(), preset: state.preset, case: state.case, overrides: state.overrides, up: state.up, runs: state.runs, recs: state.recs, geometry: state.geometry.map((g) => g.summary).filter(Boolean) };
}
export function importProject(p) {
  if (!p || typeof p !== 'object' || !p.case) throw new Error('This file is not an AeroSuite project or case file.');
  state.overrides = p.overrides && typeof p.overrides === 'object' ? p.overrides : {}; ls.set('overrides', state.overrides);
  state.up = p.up && typeof p.up === 'object' ? p.up : {}; state.runs = {}; state.recs = p.recs && typeof p.recs === 'object' ? p.recs : {}; state.results = {};
  ls.set('up', state.up); ls.set('runs', state.runs); ls.set('recs', state.recs);
  setCase(p.case, p.preset || 'custom');
  emit('run', {});
}
export { ls };
