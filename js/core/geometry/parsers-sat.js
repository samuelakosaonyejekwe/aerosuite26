// ACIS reader: SAT (text) and SAB (binary) save files → entity records → boundary-representation graph →
// in-memory STEP AP214 for the OpenCASCADE kernel.
// Record layout follows the published SAT save-file description (header, "$n" entity pointers, "#" terminators,
// "{ subtype }" blocks, "@n" length-prefixed strings from ACIS 7) and, for SAB, the tagged binary encoding as
// documented in open-source readers. Analytic geometry is translated exactly; spline curves and surfaces are
// taken from the B-spline (bs3) data the file stores for them, which is exact for "exactcur"/"exactsur" and an
// approximation within the file's fit tolerance for procedural (blend, offset, sweep, loft …) definitions.

import { str, view } from './parsers-util.js';
import { V, StepWriter } from './parsers-brep.js';

const isPtr = (t) => t !== null && typeof t === 'object' && 'p' in t;
const isStr = (t) => t !== null && typeof t === 'object' && 's' in t;
const isRev = (t) => t === true || t === 'reversed' || t === 'reverse_v' || t === 'reverse_u';

// ---------- SAT text ----------
function parseSATText(b) {
  const src = str(b, 0, b.length), n = src.length; let p = 0;
  const line = () => { const e = src.indexOf('\n', p), s = src.slice(p, e < 0 ? n : e); p = e < 0 ? n : e + 1; return s.replace(/\r$/, ''); };
  const l1 = line().trim().split(/\s+/).map(Number);
  // ACIS 1.x saves have a one-line header; from 2.0 a product line and a units/tolerance line follow
  const old = l1[0] < 200, l2 = old ? '' : line(), l3 = old ? [] : line().trim().split(/\s+/).map(Number);
  if (!(l1[0] >= 100 && l1[0] < 1e6)) throw new Error('the first line is not an ACIS save-file version record');
  const header = { version: l1[0], records: l1[1] ?? 0, bodies: l1[2] ?? 0, flags: l1[3] ?? 0, product: null, acisVersion: null, date: null, unitsMM: null, resabs: null, resnor: null };
  // header strings: "<len> text" (the '@' prefix was introduced for strings inside records, not here)
  const strs = []; for (let i = 0, k = 0; i < l2.length && k < 6; k++) { const m = /^\s*@?(\d+)\s/.exec(l2.slice(i)); if (!m) break; strs.push(l2.substr(i + m[0].length, +m[1])); i += m[0].length + +m[1]; }
  [header.product, header.acisVersion, header.date] = [strs[0]?.trim() ?? null, strs[1]?.trim() ?? null, strs[2]?.trim() ?? null];
  if (l3.length >= 2 && l3[0] > 0) { header.unitsMM = l3[0]; header.resabs = l3[1]; header.resnor = l3[2] ?? null; }
  // newer saves add a fourth header line ("T @52 <key>" or a lone flag) before the first record
  { const e = src.indexOf('\n', p), l4 = src.slice(p, e < 0 ? n : e); if (/^\s*[TF](\s+@\d+\s+\S+)?\s*$/.test(l4)) p = e < 0 ? n : e + 1; }
  const records = new Map(); let seq = 0, rec = null, depth = 0;
  const re = /\S+/g; re.lastIndex = p;
  for (let m; (m = re.exec(src));) {
    const w = m[0];
    if (!rec) {
      if (/^End-of-(ACIS|ASM)-data$/.test(w)) break;
      let index = seq, type = w;
      if (/^-\d+$/.test(w)) { index = +w.slice(1); const m2 = re.exec(src); if (!m2) break; type = m2[0]; }
      rec = { index, type, t: [] }; seq = index + 1; depth = 0; continue;
    }
    if (w === '#' && depth === 0) { records.set(rec.index, rec); rec = null; continue; }
    if (w === '{') depth++; else if (w === '}') depth--;
    if (w[0] === '$' && /^\$-?\d+$/.test(w)) rec.t.push({ p: +w.slice(1) });
    else if (w[0] === '@' && /^@\d+$/.test(w)) { const len = +w.slice(1), s0 = re.lastIndex + 1; rec.t.push({ s: src.substr(s0, len) }); re.lastIndex = s0 + len; }
    else { const v = +w; rec.t.push(v === v && /^[-+.\d]/.test(w) ? v : w); }
    if (rec.t.length > 5e7) throw new Error('an ACIS record is implausibly long (missing "#" terminator)');
  }
  return { header, records, encoding: 'SAT (text)' };
}

// ---------- SAB binary ----------
function parseSAB(b) {
  const dv = view(b), n = b.length; let p = b[1] === 0x53 /* ASM BinaryFile4 */ ? 15 : 15;
  if (str(b, 0, 15) !== 'ACIS BinaryFile' && str(b, 0, 15) !== 'ASM BinaryFile4') throw new Error('not an ACIS binary (SAB) file');
  const need = (k) => { if (p + k > n) throw new Error('the SAB stream ends inside a record (truncated file)'); };
  const i32 = () => { need(4); const v = dv.getInt32(p, true); p += 4; return v; }, f64 = () => { need(8); const v = dv.getFloat64(p, true); p += 8; return v; }, chars = (k) => { need(k); const s = str(b, p, p + k); p += k; return s; };
  const header = { version: i32(), records: i32(), bodies: i32(), flags: i32(), product: null, acisVersion: null, date: null, unitsMM: null, resabs: null, resnor: null };
  const tagStr = () => { need(1); if (b[p] !== 0x07) return null; p++; need(1); return chars(b[p++]); }, tagDbl = () => { need(1); if (b[p] !== 0x06) return null; p++; return f64(); };
  header.product = tagStr(); header.acisVersion = tagStr(); header.date = tagStr(); header.unitsMM = tagDbl(); header.resabs = tagDbl(); header.resnor = tagDbl();
  const records = new Map(); let seq = 0, rec = null, name = [];
  while (p < n) {
    const tag = b[p++];
    switch (tag) {
      case 0x04: if (!rec) { i32(); break; } rec.t.push(i32()); break;
      case 0x06: case 0x17: { const v = f64(); if (rec) rec.t.push(v); break; }
      case 0x07: { need(1); const s = chars(b[p++]); if (rec) rec.t.push({ s }); break; }
      case 0x08: { need(2); const k = dv.getUint16(p, true); p += 2; const s = chars(k); if (rec) rec.t.push({ s }); break; }
      case 0x09: case 0x12: { const k = i32(); if (!(k >= 0 && k <= n)) throw new Error('SAB string length is not plausible'); const s = chars(k); if (rec) rec.t.push({ s }); break; }
      case 0x0a: if (rec) rec.t.push(true); break;
      case 0x0b: if (rec) rec.t.push(false); break;
      case 0x0c: { const v = i32(); if (rec) rec.t.push({ p: v }); break; }
      case 0x0e: need(1); name.push(chars(b[p++])); break;
      case 0x0d: {
        need(1); name.push(chars(b[p++])); const type = name.join('-'); name = [];
        if (rec) { rec.t.push(type); break; }                        // an identifier inside a record (subtype name, keyword)
        if (/^End-of-(ACIS|ASM)-data$/.test(type)) { p = n; break; }
        rec = { index: seq++, type, t: [] }; break;
      }
      case 0x0f: if (rec) rec.t.push('{'); break;
      case 0x10: if (rec) rec.t.push('}'); break;
      case 0x11: if (rec) { records.set(rec.index, rec); rec = null; } break;
      case 0x13: case 0x14: { const x = f64(), y = f64(), z = f64(); if (rec) rec.t.push(x, y, z); break; }
      case 0x15: { const v = i32(); if (rec) rec.t.push({ e: v }); break; }
      case 0x16: { const u = f64(), v = f64(); if (rec) rec.t.push(u, v); break; }
      default: throw new Error(`unknown SAB tag 0x${tag.toString(16)} at offset ${p - 1}`);
    }
    if (rec && rec.t.length > 5e7) throw new Error('an ACIS record is implausibly long');
  }
  return { header, records, encoding: 'SAB (binary)' };
}

/** Parse an ACIS save file (text or binary) into { header, records: Map(index → { type, t: tokens, at }) }. */
export function parseACIS(b) {
  const sab = b[0] === 0x41 && (str(b, 0, 15) === 'ACIS BinaryFile' || str(b, 0, 15) === 'ASM BinaryFile4'), out = sab ? parseSAB(b) : parseSATText(b), v = out.header.version;
  // fields shared by every entity: attribute pointer; from ACIS 7 also one or two integers (tag/id) and a history pointer
  for (const r of out.records.values()) {
    const t = r.t; let i = 0;
    if (isPtr(t[0])) { r.attrib = t[0].p; i = 1; if (v >= 700) { let k = i; while (k < t.length && typeof t[k] === 'number') k++; if (k > i && isPtr(t[k])) i = k + 1; } }
    r.at = i;
  }
  return out;
}

/** Topology records: the entity pointers after the common fields (stray integers between them are version padding), and the first keyword that follows. */
function topo(r) { const p = []; let word = null; for (let i = r.at; i < r.t.length; i++) { const t = r.t[i]; if (isPtr(t)) p.push(t.p); else if (typeof t === 'string' || typeof t === 'boolean') { word = t; break; } } return { p, word }; }
/** Sequential reader over the tokens of one record. */
class Cur {
  constructor(r, i = r.at) { this.t = r.t; this.i = i; }
  next() { return this.t[this.i++]; }
  peek() { return this.t[this.i]; }
  ptr() { const t = this.t[this.i]; if (isPtr(t)) { this.i++; return t.p; } return -1; }
  num() { const t = this.t[this.i++]; if (typeof t !== 'number') throw new Error('expected a number in an ACIS record'); return t; }
  vec() { return [this.num(), this.num(), this.num()]; }
  /** One interval bound: "I" (infinite) or "F value". */
  bound() { const t = this.t[this.i++]; if (t === 'F' || t === true) return this.num(); return null; }
  word() { return this.t[this.i++]; }
}

/** B-spline data that follows a "nubs" / "nurbs" keyword. Returns a curve or surface description, or null. */
function bs3(t, i) {
  const rational = t[i] === 'nurbs'; let k = i + 1;
  const d1 = t[k++]; if (typeof d1 !== 'number') return null;
  const pairs = (count) => { const kn = [], mu = []; for (let q = 0; q < count; q++) { const a = t[k++], m = t[k++]; if (typeof a !== 'number' || typeof m !== 'number') return null; kn.push(a); mu.push(m); } return { kn, mu }; };
  const pts = (count) => { const out = []; for (let q = 0; q < count; q++) { const x = t[k++], y = t[k++], z = t[k++], w = rational ? t[k++] : 1; if (![x, y, z, w].every((v) => typeof v === 'number')) return null; out.push([x, y, z, w]); } return out; };
  // ACIS omits one multiplicity at each end of a knot vector; STEP wants the full clamped form
  const full = (mu) => mu.map((m, q) => m + (q === 0 || q === mu.length - 1 ? 1 : 0));
  if (typeof t[k] !== 'number') {                            // curve: degree form nknots
    k++; const nk = t[k++]; if (!(nk >= 2 && nk < 1e6)) return null;
    const kv = pairs(nk); if (!kv) return null; const n = kv.mu.reduce((s, m) => s + m, 0) - d1 + 1; if (!(n >= 2)) return null;
    const P = pts(n); if (!P) return null;
    return { kind: 'curve', degree: d1, knots: kv.kn, mults: full(kv.mu), poles: P.map((q) => q.slice(0, 3)), weights: rational ? P.map((q) => q[3]) : null, end: k };
  }
  const d2 = t[k++];
  while (k < t.length && !(typeof t[k] === 'number' && typeof t[k + 1] === 'number')) k++;      // rational spec, forms and singularities
  const nu = t[k++], nv = t[k++]; if (!(nu >= 2 && nv >= 2 && nu < 1e6 && nv < 1e6)) return null;
  const ku = pairs(nu), kvv = ku && pairs(nv); if (!ku || !kvv) return null;
  const cu = ku.mu.reduce((s, m) => s + m, 0) - d1 + 1, cv = kvv.mu.reduce((s, m) => s + m, 0) - d2 + 1; if (!(cu >= 2 && cv >= 2) || cu * cv > 4e6) return null;
  const P = pts(cu * cv); if (!P) return null;
  const poles = [], weights = []; for (let a = 0; a < cu; a++) { const row = [], wr = []; for (let c = 0; c < cv; c++) { const q = P[c * cu + a]; row.push(q.slice(0, 3)); wr.push(q[3]); } poles.push(row); weights.push(wr); }   // file order: u fastest
  return { kind: 'surface', udeg: d1, vdeg: d2, uknots: ku.kn, umults: full(ku.mu), vknots: kvv.kn, vmults: full(kvv.mu), poles, weights: rational ? weights : null, end: k };
}

/**
 * Translate parsed ACIS records to STEP. Returns { step, stats }.
 */
export function acisToStep(parsed) {
  const R = parsed.records, ver = parsed.header.version, W = new StepWriter();
  const stats = { bodies: 0, faces: 0, facesTranslated: 0, edges: 0, skippedFaces: {}, approximatedCurves: {}, approximatedSurfaces: {}, bodyNames: [] };
  const skip = (why) => { stats.skippedFaces[why] = (stats.skippedFaces[why] || 0) + 1; }, bump = (o, k) => { o[k] = (o[k] || 0) + 1; };
  const base = (type) => type.replace(/^t(?=(edge|coedge|vertex)$)/, '');
  // subtype blocks are numbered in file order so that "{ ref n }" can point back at them
  const subtypes = [];
  for (const r of R.values()) { if (!r.type.endsWith('-curve') && !r.type.endsWith('-surface') && r.type !== 'pcurve') continue; const st = []; for (let i = 0; i < r.t.length; i++) { if (r.t[i] === '{') { if (r.t[i + 1] === 'ref') { st.push(null); continue; } const e = { rec: r, i }; subtypes.push(e); st.push(e); } else if (r.t[i] === '}') { const e = st.pop(); if (e) e.end = i; } } }
  /** B-spline data of the subtype block opening at token i of record r (following refs), with its subtype name. */
  const splineData = (r, i, want, depth = 0) => {
    const t = r.t; if (t[i] !== '{' || depth > 8) return null;
    if (t[i + 1] === 'ref') { const s = subtypes[t[i + 2]]; return s ? splineData(s.rec, s.i, want, depth + 1) : null; }
    const name = String(t[i + 1]); let lvl = 0, found = null;
    for (let k = i; k < t.length; k++) {
      if (t[k] === '{') lvl++; else if (t[k] === '}') { if (--lvl === 0) break; }
      else if (lvl === 1 && (t[k] === 'nubs' || t[k] === 'nurbs')) { const d = bs3(t, k); if (d) { if (d.kind === want) found = d; k = d.end - 1; } }
    }
    return found ? { ...found, name, exact: /^exact/.test(name) } : { fail: name };
  };

  const emitBody = (body) => {
    const bp = topo(body).p, lump0 = bp[0] ?? -1, tr = R.get(bp[2]);
    let P = (p) => p, D = (d) => d, sc = 1, mirror = false;
    if (tr && tr.type === 'transform') {
      // the 13 numbers (rotation rows, translation, scale) are the last ones before the flag words
      // (old saves write the three flags as integers after them)
      // some writers store the whole transform of a binary file as one text literal: split it back into tokens
      for (let i = tr.t.length - 1; i >= tr.at; i--) { const x = tr.t[i]; if (x && typeof x === 'object' && typeof x.s === 'string' && /^\s*[-+0-9.]/.test(x.s) && /\s/.test(x.s.trim())) tr.t.splice(i, 1, ...x.s.trim().split(/\s+/).map((w) => (/^[-+]?(\d|\.\d)/.test(w) ? +w : w))); }
      const nums = tr.t.slice(tr.at).filter((x) => typeof x === 'number'), words = tr.t.some((x) => typeof x === 'string' || typeof x === 'boolean');      // flags as words (text) or logicals (binary)
      const q = words ? nums.slice(-13) : nums.slice(-16, -3), m = q.slice(0, 9), tv = q.slice(9, 12), s = q[12]; if (q.length < 13) throw new Error('malformed transform');
      // ACIS transforms act on row vectors: p' = (p · M) · scale + t
      const col = (j) => [m[j], m[3 + j], m[6 + j]]; sc = s || 1;
      P = (p) => [V.dot(p, col(0)) * sc + tv[0], V.dot(p, col(1)) * sc + tv[1], V.dot(p, col(2)) * sc + tv[2]]; D = (d) => [V.dot(d, col(0)), V.dot(d, col(1)), V.dot(d, col(2))];
      mirror = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]) < 0;
    }
    const cache = new Map(), memo = (key, make) => { if (!cache.has(key)) cache.set(key, make()); return cache.get(key); };
    const vertexId = (vi) => memo('v' + vi, () => { const v = R.get(vi); if (!v || base(v.type) !== 'vertex') return null; const ptrs = v.t.slice(v.at).filter(isPtr), pt = R.get(ptrs[ptrs.length - 1]?.p); if (!pt || pt.type !== 'point') return null; const p = new Cur(pt).vec(); return { id: W.vertex(P(p)), p }; });
    const curveId = (ci) => memo('c' + ci, () => {
      const c = R.get(ci); if (!c) return { fail: 'missing curve' }; const k = new Cur(c);
      try {
        if (c.type === 'straight-curve') { const root = k.vec(), dir = k.vec(); return { id: W.line(P(root), D(dir)), rev: false }; }
        if (c.type === 'ellipse-curve') { const cen = k.vec(), nrm = k.vec(), maj = k.vec(), ratio = k.num(), a = V.len(maj) * sc; if (!(a > 0)) return { fail: 'ellipse-curve (zero radius)' }; return { id: Math.abs(ratio - 1) < 1e-12 ? W.circle(P(cen), D(nrm), D(maj), a) : ratio < 1 ? W.ellipse(P(cen), D(nrm), D(maj), a, a * ratio) : W.ellipse(P(cen), D(nrm), D(V.cross(nrm, maj)), a * ratio, a), rev: false }; }
        if (c.type === 'intcurve-curve') {
          const rev = isRev(k.word()), d = splineData(c, k.i, 'curve'); if (!d || d.fail) return { fail: `intcurve-curve (${d?.fail ?? 'no subtype data'})` };
          if (!d.exact) bump(stats.approximatedCurves, d.name);
          return { id: W.bspline({ degree: d.degree, poles: d.poles.map(P), knots: d.knots, mults: d.mults, weights: d.weights }), rev };
        }
      } catch (e) { return { fail: `${c.type} (malformed record)` }; }
      return { fail: c.type };
    });
    const edgeId = (ei) => memo('e' + ei, () => {
      const e = R.get(ei); if (!e || base(e.type) !== 'edge') return { fail: 'missing edge' };
      const tp = topo(e), sv = tp.p[0], ev = tp.p[1], ci = tp.p[3], rev = isRev(tp.word);
      const s = vertexId(sv), t = vertexId(ev), c = curveId(ci); if (c.fail) return { fail: c.fail }; if (!s || !t) return { fail: 'edge without vertices' };
      stats.edges++; return { id: W.edge(s.id, t.id, c.id, rev === c.rev) };
    });
    const surfaceId = (si) => memo('s' + si, () => {
      const s = R.get(si); if (!s) return { fail: 'face without a surface' }; const k = new Cur(s);
      try {
        switch (s.type) {
          case 'plane-surface': { const root = k.vec(), nrm = k.vec(), u = k.vec(); return { id: W.plane(P(root), D(nrm), D(u)) }; }
          case 'cone-surface': {
            const cen = k.vec(), axis = k.vec(), maj = k.vec(), ratio = k.num(); k.bound(); k.bound(); const sin = k.num(), cos = k.num(), r = V.len(maj) * sc;
            if (Math.abs(ratio - 1) > 1e-9) return { fail: 'cone-surface (elliptical)' };
            if (Math.abs(sin) < 1e-12) return { id: W.cylinder(P(cen), D(axis), D(maj), r), flip: cos < 0 };
            const tn = sin / cos;                              // radius grows by tan(half angle) per unit along +axis
            return { id: W.cone(P(cen), D(tn > 0 ? axis : V.mul(axis, -1)), D(maj), r, Math.atan(Math.abs(tn))), flip: cos < 0 };
          }
          case 'sphere-surface': { const cen = k.vec(), rad = k.num(), u = k.vec(), pole = k.vec(); return { id: W.sphere(P(cen), D(pole), D(u), Math.abs(rad) * sc), flip: rad < 0, pole: V.add(cen, V.mul(V.unit(pole), Math.abs(rad))) }; }
          case 'torus-surface': { const cen = k.vec(), nrm = k.vec(), major = k.num(), minor = k.num(), u = k.vec(); if (!(Math.abs(major) > Math.abs(minor))) return { fail: 'torus-surface (self-intersecting form)' }; return { id: W.torus(P(cen), D(nrm), D(u), Math.abs(major) * sc, Math.abs(minor) * sc), flip: minor < 0, pole: V.add(cen, V.mul(V.unit(u), Math.abs(major) + Math.abs(minor))) }; }
          case 'spline-surface': {
            const rev = isRev(k.word()), d = splineData(s, k.i, 'surface'); if (!d || d.fail) return { fail: `spline-surface (${d?.fail ?? 'no subtype data'})` };
            if (!d.exact) bump(stats.approximatedSurfaces, d.name);
            return { id: W.bsurface({ udeg: d.udeg, vdeg: d.vdeg, poles: d.poles.map((row) => row.map(P)), uknots: d.uknots, umults: d.umults, vknots: d.vknots, vmults: d.vmults, weights: d.weights }), flip: rev };
          }
          default: return { fail: s.type };
        }
      } catch (e) { return { fail: `${s.type} (malformed record)` }; }
    });

    const shells = [], faceColors = [];
    for (let li = lump0, gl = 0; li >= 0 && gl < 1e5; gl++) {
      const lump = R.get(li); if (!lump || lump.type !== 'lump') break; const lp0 = topo(lump).p, nextLump = lp0[0] ?? -1;
      for (let si = lp0[1] ?? -1, gs = 0; si >= 0 && gs < 1e5; gs++) {
        const sh = R.get(si); if (!sh || sh.type !== 'shell') break; const sp = topo(sh).p, nextShell = sp[0] ?? -1, shell = { faces: [], closed: false };
        for (let fi = sp[2] ?? -1, gf = 0; fi >= 0 && gf < 1e7; gf++) {
          const face = R.get(fi); if (!face || face.type !== 'face') break; stats.faces++;
          const ft = topo(face), nextFace = ft.p[0] ?? -1, loop0 = ft.p[1] ?? -1, surf = surfaceId(ft.p[4] ?? -1), frev = isRev(ft.word); fi = nextFace;
          if (surf.fail) { skip(surf.fail); continue; }
          const bounds = []; let bad = null;
          for (let lpi = loop0, g = 0; lpi >= 0 && g < 1e5 && !bad; g++) {
            const lp = R.get(lpi); if (!lp || lp.type !== 'loop') break; const lt = topo(lp).p, nextLoop = lt[0] ?? -1, first = lt[1] ?? -1, oriented = [];
            for (let ce = first, q = 0; ce >= 0 && q < 1e5; q++) {
              const co = R.get(ce); if (!co || base(co.type) !== 'coedge') break; const ct = topo(co), nx = ct.p[0] ?? -1, ei = ct.p[3] ?? -1, crev = isRev(ct.word);
              if (ei >= 0) { const e = edgeId(ei); if (e.fail) { bad = `edge on ${e.fail}`; break; } oriented.push(W.oriented(e.id, !crev)); }
              ce = nx; if (ce === first) break;
            }
            if (oriented.length) bounds.push({ loop: W.loop(oriented) });
            lpi = nextLoop;
          }
          if (bad) { skip(bad); continue; }
          if (!bounds.length) { if (!surf.pole) { skip('unbounded face'); continue; } bounds.push({ loop: W.vertexLoop(W.vertex(P(surf.pole))) }); }
          let same = !frev; if (surf.flip) same = !same; if (mirror) same = !same;
          const fid = W.face(bounds, surf.id, same), col = colourOf(face); if (col) faceColors.push([fid, col]);
          shell.faces.push(fid); stats.facesTranslated++;
        }
        shells.push(shell); si = nextShell;
      }
      li = nextLump;
    }
    const name = nameOf(body) || `body ${stats.bodies + 1}`;
    W.body(name, shells, { color: colourOf(body), faceColors }); stats.bodies++; stats.bodyNames.push(name);
  };
  /** Walk an entity's attribute chain. */
  const attribs = (r) => { const out = []; for (let ai = r.attrib, g = 0; ai >= 0 && g < 1000; g++) { const a = R.get(ai); if (!a) break; out.push(a); const c = new Cur(a); ai = c.ptr(); } return out; };
  const colourOf = (r) => { for (const a of attribs(r)) { if (/^rgb_color-st-attrib$/.test(a.type)) { const v = a.t.slice(a.at).filter((x) => typeof x === 'number').slice(-3); if (v.length === 3 && v.every((x) => x >= 0 && x <= 1)) return v; } if (/truecolor-adesk-attrib$/.test(a.type)) { const v = a.t.slice(a.at).filter((x) => typeof x === 'number').pop(); if (Number.isInteger(v)) return [((v >> 16) & 255) / 255, ((v >> 8) & 255) / 255, (v & 255) / 255]; } } return null; };
  const nameOf = (r) => { let fallback = null; for (const a of attribs(r)) { if (!/name_attrib/.test(a.type)) continue; const s = a.t.filter(isStr).map((x) => x.s); if (s.length >= 2 && /name/i.test(s[0]) && s[1].trim()) return s[1].trim(); if (s.length) fallback ??= s[s.length - 1].trim(); } return fallback || null; };

  for (const r of R.values()) if (r.type === 'body') { try { emitBody(r); } catch (e) { if (typeof process === 'object' && process.env?.SAT_DEBUG) console.error(e); skip(`body ${r.index} (malformed topology records)`); } }
  const unit = parsed.header.unitsMM;
  void ver;
  return { step: stats.facesTranslated ? W.finish({ unit: '.MILLI.', tol: Math.max(1e-7, Math.min(1e-3, parsed.header.resabs || 1e-6)), source: 'ACIS' }) : null, stats, unitsMM: unit };
}
