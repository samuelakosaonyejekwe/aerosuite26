// The engineering specification (equations, models, initial/boundary conditions, calibration,
// verification, validation, inputs, outputs and geometry requirements for each suite), loaded on demand.

let cache = null;
export async function loadSpec() {
  if (cache) return cache;
  if (globalThis.__AEROSUITE_SPEC__) return (cache = globalThis.__AEROSUITE_SPEC__);
  const r = await fetch(new URL('../data/spec.json', import.meta.url));
  if (!r.ok) throw new Error('Specification data is not available');
  return (cache = await r.json());
}
/** Split a specification sentence list ("... include A, B and C.") into individual terms. */
export function terms(text) {
  if (!text) return [];
  let t = text.replace(/^[^:]*?(?:include|includes|are|shall include|cover)\s+/i, '');
  const stop = t.search(/\.\s+[A-Z]/); if (stop > 0) t = t.slice(0, stop);
  return t.replace(/\.$/, '').split(/,\s*(?:and\s+)?|;\s*|\s+and\s+(?=[a-zA-Z-]+\s+(?:equations?|models?|relations?|formulations?|theorem|law|criteri\w+)\b)/).map((s) => s.trim().replace(/^and\s+/, '')).filter((s) => s.length > 2 && s.length < 110);
}
