// Rows of the compliance pack, derived from the licence registry (js/data/licences.json). Shared by the generator
// (tools/compliance.mjs → COMPLIANCE.md) and the in-app printable view on the About page.

const ROUTE = { api: 'contacted from the user’s browser', 'web-data': 'read by the scheduled snapshot job (server-side)', 'bundled-dataset': 'shipped with the app', citations: 'shipped with the app (references)', 'data-licence': 'data licence', operator: 'configured by the operator', 'user-service': 'contacted with the user’s own account' };
/** Is the source used in the given mode ('commercial' = no keys, nothing accepted; 'nonCommercial' = default configuration)? */
export function usedIn(s, mode) {
  if (s.used === false) return false;
  if (mode === 'nonCommercial') return s.class !== 'commercial-with-key' || s.kind === 'user-service' || s.id === 'overpass-public';
  return s.class === 'commercial-ok' && s.usedInCommercial !== false || s.kind === 'user-service';
}
/** One row per source contacted or shipped in the mode. */
export function complianceRows(reg, mode) {
  return reg.sources.filter((s) => usedIn(s, mode)).map((s) => ({ id: s.id, name: s.name, use: `${ROUTE[s.kind] || s.kind}: ${s.usedFor}`, licence: s.licence, obligations: (s.obligations || []).join('; '), shownAt: s.shownAt || '', attribution: s.attribution || '', quote: s.quote ? `“${s.quote}”` : '—', url: s.url || '—', retrieved: s.retrieved }));
}
/** Items that still need a human decision, each as a yes/no question. */
export const complianceQuestions = (reg) => reg.sources.filter((s) => s.question && s.used !== false).map((s) => ({ id: s.id, name: s.name, question: s.question }));
