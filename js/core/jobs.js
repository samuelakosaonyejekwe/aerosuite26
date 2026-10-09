// Job execution shared by the Web Worker and the main-thread fallback. A job is plain data in,
// plain data out, so it can cross the worker boundary.

import { loadSuite, executionOrder, SUITES } from './registry.js';
import * as S from './studies.js';
import { registerCustom, extendOptions } from './matlib.js';

// Solvers that march a 3-D or coupled field in time. They are optional in the integrated run because each takes seconds.
export const HEAVY = { cfd: ['dns3d', 'les3d', 'rans3d', 'bluff3d', 'cavity'], crash: ['barrel3d', 'aircraft3d', 'impact3d'], aeroelastic: ['uvlmfsi', 'cfdfsi'] };

export async function execute(job, progress = () => {}) {
  if (job.materials) registerCustom(job.materials); // the user's own materials, selectable like built-in ones
  const def = job.suite ? await loadSuite(job.suite) : null;
  const an = def && job.analysis ? def.analyses.find((a) => a.id === job.analysis) : null;
  const ctx = job.case ? S.makeCtx(job.case, job.up || {}, progress) : null;
  const inputs = () => S.resolveInputs(an, ctx, job.overrides || {});
  switch (job.kind) {
    case 'describe': { // serialisable description of a suite for the UI
      const c = S.makeCtx(job.case, job.up || {});
      return {
        id: def.id, n: def.n, tagline: def.tagline, consumes: def.consumes || [], provides: def.provides || [], handoff: def.handoff || [],
        analyses: def.analyses.map((a) => {
          const ok = S.applicable(a, c), r = S.resolveInputs(a, c, (job.overridesBy || {})[a.id] || {});
          return { id: a.id, title: a.title, summary: a.summary, fidelity: a.fidelity, equations: a.equations || [], inputs: a.inputs.map((f) => (f.type === 'select' ? { ...f, options: extendOptions(f.options) } : f)), values: r.inp, linked: r.linked, applicable: ok,
            convergence: a.convergence ? { param: a.convergence.param, label: a.convergence.label, levels: a.convergence.levels, metric: a.convergence.metric } : null,
            hasVerify: !!a.verify, validation: a.validation || [], calibration: a.calibration || null };
        }),
      };
    }
    case 'run': {
      const { inp, linked } = inputs(), res = await S.runAnalysis(an, inp, ctx);
      let recs = [];
      try { recs = (an.recommend ? an.recommend(res, inp, ctx) : []) || []; } catch { recs = []; }
      return { res, inp, linked, recs };
    }
    case 'convergence': return S.convergenceStudy(an, inputs().inp, ctx, job.opts || {});
    case 'sensitivity': return S.sensitivity(an, inputs().inp, ctx, job.opts.metric, job.opts.keys, job.opts.pct);
    case 'sweep': return S.sweep(an, inputs().inp, ctx, job.opts.key, job.opts.values);
    case 'uq': return S.uqStudy(an, inputs().inp, ctx, job.opts);
    case 'calibrate': return S.calibrate(an, inputs().inp, ctx, job.opts);
    case 'validate': return S.validate(an, inputs().inp, ctx, job.opts);
    case 'verify': return S.verifySuite(def);
    case 'verifyAll': {
      const out = [];
      for (let i = 0; i < SUITES.length; i++) { progress(i / SUITES.length, SUITES[i].short); const d = await loadSuite(SUITES[i].id); out.push({ suite: d.id, n: d.n, checks: S.verifySuite(d) }); }
      return out;
    }
    case 'integrated': { // run every applicable analysis of every (selected) suite in dependency order
      const defs = []; for (const m of SUITES) if (!job.only || job.only.includes(m.id)) defs.push(await loadSuite(m.id)); // a single-suite job loads only that suite
      const order = executionOrder(defs);
      const up = JSON.parse(JSON.stringify(job.up || {})), log = [], total = order.reduce((n, id) => n + defs.find((d) => d.id === id).analyses.length, 0); let done = 0;
      for (const id of order) {
        const d = defs.find((x) => x.id === id); up[id] = up[id] || {};
        for (const a of d.analyses) {
          progress(done++ / total, a.title);
          const c = S.makeCtx(job.case, up), ok = S.applicable(a, c);
          if (ok !== true) { log.push({ suite: id, analysis: a.id, title: a.title, skipped: ok }); continue; }
          if (!job.includeHeavy && HEAVY[id]?.includes(a.id)) { log.push({ suite: id, analysis: a.id, title: a.title, skipped: 'Heavy 3-D solver, left out for speed: tick “Include the heavy 3-D solvers” or run it from its suite' }); continue; }
          try {
            const { inp, linked } = S.resolveInputs(a, c, (job.overridesAll || {})[`${id}.${a.id}`] || {}), res = await S.runAnalysis(a, inp, c);
            let recs = []; try { recs = (a.recommend ? a.recommend(res, inp, c) : []) || []; } catch { recs = []; }
            Object.assign(up[id], Object.fromEntries(Object.entries(res.outputs).filter(([, v]) => typeof v === 'number' || Array.isArray(v))));
            log.push({ suite: id, analysis: a.id, title: a.title, payload: { res, inp, linked, recs } });
          } catch (e) { log.push({ suite: id, analysis: a.id, title: a.title, error: e.message }); }
        }
      }
      progress(1, 'done');
      return { order, log };
    }
    default: throw new Error('Unknown job kind ' + job.kind);
  }
}
