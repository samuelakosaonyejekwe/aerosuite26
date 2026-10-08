// Readers for exact-CAD exchange files (STEP, IGES, Parasolid text) and the header-only inspectors for
// formats that are recognised but not decoded (ACIS, JT, HDF5/NetCDF containers, LAZ, E57, GeoTIFF …).
// None of these tessellates exact B-rep geometry: what is returned is metadata, a census, and whatever
// explicit points / polylines / facets the file itself carries.

import { Scanner, str, nums, toNum, view, LIMITS } from './parsers-util.js';
import { lasHeader } from './parsers-surface.js';

// ---------- STEP (ISO 10303-21) ----------
/** Parse a STEP argument list (the text between the outer parentheses) into nested arrays. Iterative, no recursion. */
function stepArgs(s) {
  const root = [], stack = [root], n = s.length; let i = 0;
  while (i < n) {
    const c = s[i], top = stack[stack.length - 1];
    if (c === "'") {
      let j = i + 1, out = '';
      for (;;) { const k = s.indexOf("'", j); if (k < 0) { out += s.slice(j); j = n; break; } out += s.slice(j, k); if (s[k + 1] === "'") { out += "'"; j = k + 2; } else { j = k + 1; break; } }
      top.push(out); i = j;
    } else if (c === '(') { const a = []; top.push(a); stack.push(a); i++; }
    else if (c === ')') { if (stack.length > 1) stack.pop(); i++; }
    else if (c === ',' || c <= ' ') i++;
    else { let j = i; while (j < n && s[j] > ' ' && s[j] !== ',' && s[j] !== '(' && s[j] !== ')' && s[j] !== "'") j++; const t = s.slice(i, j), v = +t; top.push(t !== '' && v === v && /^[-+.\d]/.test(t) ? v : t); i = j; }
  }
  return root;
}
/** Decode the \X2\…\X0\ and \X\hh escapes used in STEP strings. */
const stepStr = (s) => (typeof s === 'string' && s.includes('\\') ? s.replace(/\\X2\\([0-9A-F]+)\\X0\\/g, (m, h) => h.match(/.{4}/g).map((x) => String.fromCharCode(parseInt(x, 16))).join('')).replace(/\\X\\([0-9A-F]{2})/g, (m, h) => String.fromCharCode(parseInt(h, 16))) : s);
const STEP_UNIT_WORDS = [[/^(INCH|IN)$/i, 'in'], [/^(FOOT|FEET|FT)$/i, 'ft'], [/^(MILLIMETRE|MILLIMETER|MM)$/i, 'mm'], [/^(CENTIMETRE|CENTIMETER|CM)$/i, 'cm'], [/^(METRE|METER|M)$/i, 'm']];
const STEP_TESS = new Set(['COORDINATES_LIST', 'TRIANGULATED_FACE', 'COMPLEX_TRIANGULATED_FACE', 'TRIANGULATED_SURFACE_SET', 'COMPLEX_TRIANGULATED_SURFACE_SET']);

export function readSTEP(b, ctx) {
  const n = b.length, headEnd = Math.min(n, 1 << 18), head = str(b, 0, headEnd).replace(/\/\*[\s\S]*?\*\//g, '');
  const meta = { schema: [], ap: null, description: null, fileName: null, timeStamp: null, author: null, organization: null, preprocessor: null, originatingSystem: null, products: [], assemblyUsages: 0, entityCount: 0, census: {}, lengthUnit: null, cartesianPoints: 0, tessellated: { coordinateLists: 0, faces: 0, triangles: 0 }, polylines: 0 };
  let m = /FILE_SCHEMA\s*\(([\s\S]*?)\)\s*;/.exec(head); if (m) meta.schema = stepArgs(m[1]).flat(3).filter((x) => typeof x === 'string');
  m = /FILE_DESCRIPTION\s*\(([\s\S]*?)\)\s*;/.exec(head); if (m) { const a = stepArgs(m[1]); meta.description = (Array.isArray(a[0]) ? a[0].join(' ') : String(a[0] ?? '')).trim() || null; }
  m = /FILE_NAME\s*\(([\s\S]*?)\)\s*;\s*FILE_SCHEMA/.exec(head) || /FILE_NAME\s*\(([\s\S]*?)\)\s*;/.exec(head);
  if (m) { const a = stepArgs(m[1]), one = (v) => stepStr(Array.isArray(v) ? v.filter(Boolean).join('; ') : v) || null; meta.fileName = one(a[0]); meta.timeStamp = one(a[1]); meta.author = one(a[2]); meta.organization = one(a[3]); meta.preprocessor = one(a[4]); meta.originatingSystem = one(a[5]); }
  const sch = meta.schema.join(' ').toUpperCase();
  meta.ap = /AP242|MANAGED_MODEL_BASED/.test(sch) ? 'AP242' : /AUTOMOTIVE_DESIGN|AP214/.test(sch) ? 'AP214' : /CONFIG_CONTROL_DESIGN|CONFIGURATION_CONTROL|AP203/.test(sch) ? 'AP203' : null;

  // DATA section: scan instance by instance ("#id = NAME(args);") over the bytes; only entities that are
  // used are decoded to strings, everything else is just counted.
  let p = 0;   // start after the HEADER section's ENDSEC
  for (let i = 0; i + 6 < Math.min(n, 1 << 20); i++) if (b[i] === 69 && b[i + 1] === 78 && b[i + 2] === 68 && b[i + 3] === 83 && b[i + 4] === 69 && b[i + 5] === 67) { p = i + 6; break; }
  const census = new Map(), pts = [], ptIndex = new Map(), polylines = [], coordLists = new Map(), faces = [], lengthUnits = new Map(), unitCtx = [], products = new Set();
  const count = (k) => census.set(k, (census.get(k) || 0) + 1);
  while (p < n) {
    while (p < n && b[p] !== 35) { if (b[p] === 69 && b[p + 1] === 78 && b[p + 2] === 68 && b[p + 3] === 83 && b[p + 4] === 69 && b[p + 5] === 67) { p = n; break; } p++; }   // stop at ENDSEC
    if (p >= n) break;
    let q = p + 1, id = 0; while (q < n && b[q] >= 48 && b[q] <= 57) { id = id * 10 + b[q] - 48; q++; }
    while (q < n && b[q] <= 32) q++;
    if (b[q] !== 61) { p = q + 1; continue; }
    q++; while (q < n && b[q] <= 32) q++;
    const nameStart = q; while (q < n && ((b[q] >= 65 && b[q] <= 90) || (b[q] >= 48 && b[q] <= 57) || b[q] === 95)) q++;
    const name = str(b, nameStart, q);
    let end = q, inStr = false; for (; end < n; end++) { const c = b[end]; if (c === 39) inStr = !inStr; else if (c === 59 && !inStr) break; }
    meta.entityCount++;
    if (name === 'CARTESIAN_POINT') {
      count(name);
      const s = str(b, q, end), o = s.lastIndexOf('('), c = s.indexOf(')', o), v = o < 0 ? [] : nums(s.slice(o + 1, c < 0 ? undefined : c));
      if (pts.length < 3 * LIMITS.verts) { ptIndex.set(id, pts.length / 3); pts.push(v[0] ?? 0, v[1] ?? 0, v[2] ?? 0); }
    } else if (name === '') {                               // complex entity: ( A(...) B(...) … )
      const s = str(b, q, end);
      for (const mm of s.matchAll(/([A-Z_][A-Z0-9_]*)\s*\(/g)) count(mm[1]);
      if (s.includes('LENGTH_UNIT')) {
        const si = /SI_UNIT\s*\(\s*(\.\w+\.|\$)\s*,\s*\.METRE\.\s*\)/.exec(s), cb = /CONVERSION_BASED_UNIT\s*\(\s*'([^']*)'/.exec(s);
        if (cb) { const w = cb[1].trim(), hit = STEP_UNIT_WORDS.find(([re]) => re.test(w)); lengthUnits.set(id, { unit: hit ? hit[1] : null, label: w }); }
        else if (si) { const pre = si[1].replace(/\./g, ''); lengthUnits.set(id, { unit: { $: 'm', MILLI: 'mm', CENTI: 'cm' }[pre] ?? null, label: pre === '$' ? 'METRE' : `${pre}METRE` }); }
      }
      if (s.includes('GLOBAL_UNIT_ASSIGNED_CONTEXT')) { const g = /GLOBAL_UNIT_ASSIGNED_CONTEXT\s*\(\s*\(([^)]*)\)/.exec(s); if (g) unitCtx.push(...(g[1].match(/#\d+/g) || []).map((r) => +r.slice(1))); }
    } else {
      count(name);
      if (name === 'PRODUCT') { const a = stepArgs(str(b, q, end)); const nm = stepStr(a[0]?.[1] || a[0]?.[0]); if (nm && products.size < 500) products.add(String(nm)); }
      else if (name === 'NEXT_ASSEMBLY_USAGE_OCCURRENCE') meta.assemblyUsages++;
      else if (name === 'POLYLINE') { const s = str(b, q, end); polylines.push((s.match(/#\d+/g) || []).map((r) => +r.slice(1))); }
      else if (name === 'GLOBAL_UNIT_ASSIGNED_CONTEXT') unitCtx.push(...(str(b, q, end).match(/#\d+/g) || []).map((r) => +r.slice(1)));
      else if (STEP_TESS.has(name)) {
        const a = stepArgs(str(b, q, end))[0] || [];
        if (name === 'COORDINATES_LIST') coordLists.set(id, (a[2] || []).filter(Array.isArray));
        else { const set = name.endsWith('SURFACE_SET'), ref = typeof a[1] === 'string' ? +a[1].slice(1) : NaN; faces.push({ ref, pn: Array.isArray(a[set ? 4 : 5]) ? a[set ? 4 : 5] : [], data: a.slice(set ? 5 : 6).filter(Array.isArray), complex: name.startsWith('COMPLEX') }); }
      }
    }
    p = end + 1;
  }
  meta.census = Object.fromEntries([...census].sort((x, y) => y[1] - x[1]).slice(0, 400));
  meta.products = [...products]; meta.cartesianPoints = pts.length / 3; meta.polylines = polylines.length;
  if (!meta.entityCount) throw Object.assign(new Error('no entity instances were found in the DATA section'), { meta });

  // length unit: the one referenced by the global unit context, else the only one present
  let cand = unitCtx.map((r) => lengthUnits.get(r)).filter(Boolean); if (!cand.length) cand = [...lengthUnits.values()];
  const distinct = [...new Set(cand.map((u) => u.label))]; let units = { length: null, source: null };
  if (distinct.length === 1) { meta.lengthUnit = distinct[0]; if (cand[0].unit) units = { length: cand[0].unit, source: 'file' }; else ctx.warn(`The STEP length unit "${distinct[0]}" has no platform equivalent; assign the units manually.`); }
  else if (distinct.length > 1) { meta.lengthUnit = distinct.join(' / '); ctx.warn(`The STEP file declares several length units (${distinct.join(', ')}); the units must be confirmed manually.`); }

  // AP242 tessellated geometry, when the file carries it
  const tp = [], tri = [], base = new Map();
  for (const [id, list] of coordLists) { base.set(id, tp.length / 3); for (const c of list) tp.push(+c[0] || 0, +c[1] || 0, +c[2] || 0); }
  meta.tessellated.coordinateLists = coordLists.size;
  for (const f of faces) {
    const b0 = base.get(f.ref), cl = coordLists.get(f.ref); if (b0 === undefined) continue;
    const ix = (v) => { const k = f.pn.length ? f.pn[v - 1] : v; return k >= 1 && k <= cl.length ? b0 + k - 1 : -1; };
    const T = (a, c, d) => { a = ix(a); c = ix(c); d = ix(d); if (a >= 0 && c >= 0 && d >= 0) tri.push(a, c, d); };
    meta.tessellated.faces++;
    if (!f.complex) { for (const t of f.data[0] || []) if (Array.isArray(t)) T(t[0], t[1], t[2]); }
    else {
      for (const s of f.data[0] || []) if (Array.isArray(s)) for (let i = 2; i < s.length; i++) (i % 2 ? T(s[i - 1], s[i - 2], s[i]) : T(s[i - 2], s[i - 1], s[i]));   // strips
      for (const s of f.data[1] || []) if (Array.isArray(s)) for (let i = 2; i < s.length; i++) T(s[0], s[i - 1], s[i]);                                               // fans
    }
  }
  meta.tessellated.triangles = tri.length / 3;
  const brep = (census.get('ADVANCED_FACE') || 0) + (census.get('FACE_SURFACE') || 0) + (census.get('MANIFOLD_SOLID_BREP') || 0) + (census.get('SHELL_BASED_SURFACE_MODEL') || 0);
  if (tri.length) {
    if (brep) ctx.warn(`The file carries both exact B-rep geometry (${brep} face/solid entities) and an AP242 tessellation; the tessellation is shown, the exact faces are not evaluated.`);
    return { positions: Float64Array.from(tp), triangles: Uint32Array.from(tri), kind: 'surface', units, meta };
  }
  const lines = [];
  for (const pl of polylines) for (let i = 1; i < pl.length; i++) { const a = ptIndex.get(pl[i - 1]), c = ptIndex.get(pl[i]); if (a !== undefined && c !== undefined) lines.push(a, c); }
  ctx.warn(brep
    ? `Exact B-rep geometry (${census.get('ADVANCED_FACE') || 0} ADVANCED_FACE, ${census.get('MANIFOLD_SOLID_BREP') || 0} MANIFOLD_SOLID_BREP) is NOT tessellated here. The ${pts.length / 3} CARTESIAN_POINTs (vertices, control points, placement origins) are shown as a point cloud; extents are indicative only because control points and placement origins may lie off the surface. No area, volume or section can be computed from this file.`
    : 'No tessellated faces were found; the CARTESIAN_POINTs are shown as a point cloud.');
  return { positions: Float64Array.from(pts), lines: Uint32Array.from(lines), kind: 'cad-brep', units, meta };
}

// ---------- IGES ----------
const IGES_UNITS = { 1: ['in', 'inch'], 2: ['mm', 'millimetre'], 4: ['ft', 'foot'], 5: [null, 'mile'], 6: ['m', 'metre'], 7: [null, 'kilometre'], 8: [null, 'mil'], 9: [null, 'micron'], 10: ['cm', 'centimetre'], 11: [null, 'microinch'] };
const IGES_NAMES = { 100: 'circular arc', 102: 'composite curve', 104: 'conic arc', 106: 'copious data', 108: 'plane', 110: 'line', 112: 'parametric spline curve', 114: 'parametric spline surface', 116: 'point', 118: 'ruled surface', 120: 'surface of revolution', 122: 'tabulated cylinder', 124: 'transformation matrix', 126: 'rational B-spline curve', 128: 'rational B-spline surface', 130: 'offset curve', 140: 'offset surface', 141: 'boundary', 142: 'curve on parametric surface', 143: 'bounded surface', 144: 'trimmed parametric surface', 186: 'manifold solid B-rep object', 190: 'plane surface', 192: 'right circular cylindrical surface', 194: 'right circular conical surface', 196: 'spherical surface', 198: 'toroidal surface', 212: 'general note', 308: 'subfigure definition', 314: 'colour definition', 402: 'associativity instance', 406: 'property', 408: 'singular subfigure instance', 502: 'vertex list', 504: 'edge list', 508: 'loop', 510: 'face', 514: 'shell' };
/** Split the IGES global section into its fields (Hollerith strings "nHtext" honoured). */
function igesGlobal(g) {
  let i = 0, pd = ',', rd = ';';
  let m = /^1H(.)/.exec(g); if (m) { pd = m[1]; i = 3; } if (g[i] === pd) i++;
  m = /^1H(.)/.exec(g.slice(i)); if (m) { rd = m[1]; i += 3; } if (g[i] === pd) i++;
  const out = [pd, rd];
  for (let k = 0; k < 40 && i < g.length; k++) {
    while (g[i] === ' ') i++;
    const h = /^(\d+)H/.exec(g.slice(i, i + 12)); let v;
    if (h) { v = g.substr(i + h[0].length, +h[1]); i += h[0].length + +h[1]; } else { let j = i; while (j < g.length && g[j] !== pd && g[j] !== rd) j++; v = g.slice(i, j).trim(); i = j; }
    out.push(v);
    while (i < g.length && g[i] !== pd && g[i] !== rd) i++;
    if (g[i] === rd || i >= g.length) break;
    i++;
  }
  return { f: out, pd, rd };
}
export function readIGES(b, ctx) {
  const sc = new Scanner(b); let start = '', glob = ''; const dir = new Map(), pdata = new Map(), census = new Map(), WANT = new Set([116, 110, 106, 126, 128]);
  let dPrev = null, withXf = 0, lines = 0, flavour = null;
  for (let ln; (ln = sc.line()) !== null;) {
    if (ln.length < 73) continue;
    const sec = ln[72]; lines++;
    if (sec === 'S') start += ln.slice(0, 72).trimEnd() + ' ';
    else if (sec === 'G') glob += ln.slice(0, 72);
    else if (sec === 'D') {
      if (dPrev === null) dPrev = ln;
      else {
        const type = parseInt(dPrev.slice(0, 8), 10), seq = parseInt(dPrev.slice(73, 80), 10), xf = parseInt(dPrev.slice(48, 56), 10) || 0, status = dPrev.slice(64, 72).padStart(8, '0'), form = parseInt(ln.slice(32, 40), 10) || 0;
        census.set(type, (census.get(type) || 0) + 1);
        if (WANT.has(type)) dir.set(seq, { type, form, xf, use: status.slice(4, 6) });
        dPrev = null;
      }
    } else if (sec === 'P') { const de = parseInt(ln.slice(64, 72), 10), d = dir.get(de); if (d) pdata.set(de, (pdata.get(de) || '') + ln.slice(0, 64)); }
    else if (sec === 'C' || sec === 'B') flavour = sec === 'C' ? 'compressed' : 'binary';
  }
  if (flavour) throw new Error(`${flavour} IGES is not read; export fixed-format ASCII IGES`);
  if (!glob) throw new Error('no IGES global section was found');
  const { f, pd, rd } = igesGlobal(glob), flag = parseInt(f[13], 10), u = IGES_UNITS[flag] ?? (flag === 3 ? [null, f[14]] : [null, null]);
  const uname = String(f[14] || '').toUpperCase(), length = u[0] ?? ({ IN: 'in', INCH: 'in', MM: 'mm', FT: 'ft', M: 'm', CM: 'cm' }[uname] || null);
  const meta = { start: start.trim().slice(0, 400), productId: f[2] || null, fileName: f[3] || null, originatingSystem: f[4] || null, preprocessor: f[5] || null, modelScale: toNum(f[12] || '') || null, unitsFlag: flag === flag ? flag : null, unitsName: f[14] || u[1] || null, created: f[17] || null, minResolution: toNum(f[18] || '') || null, author: f[20] || null, organization: f[21] || null, igesVersionFlag: parseInt(f[22], 10) || null, lineCount: lines, entityCount: 0, census: {}, read: { points: 0, lines: 0, copious: 0, curveControlPoints: 0, surfaceControlPoints: 0 }, skippedParametricCurves: 0 };
  for (const [t, c] of [...census].sort((x, y) => y[1] - x[1])) { meta.census[`${t} ${IGES_NAMES[t] ?? 'entity'}`] = c; meta.entityCount += c; }
  const xyz = [], seg = [];
  const add = (x, y, z) => { xyz.push(x, y, z); return xyz.length / 3 - 1; };
  for (const [de, raw] of pdata) {
    const d = dir.get(de), e = raw.indexOf(rd), v = (e < 0 ? raw : raw.slice(0, e)).split(pd).map((s) => toNum(s.trim()));
    if (d.use === '05') { meta.skippedParametricCurves++; continue; }      // 2-D curves in a surface's parameter space
    if (d.xf) withXf++;
    if (d.type === 116 && v.length >= 4) { add(v[1], v[2], v[3]); meta.read.points++; }
    else if (d.type === 110 && v.length >= 7) { seg.push(add(v[1], v[2], v[3]), add(v[4], v[5], v[6])); meta.read.lines++; }
    else if (d.type === 106) {
      const ip = v[1], np = v[2]; if (!(np >= 1 && np <= v.length)) continue;
      const step = ip === 1 ? 2 : ip === 2 ? 3 : 6, o = ip === 1 ? 4 : 3; let prev = -1;
      for (let k = 0; k < np && o + k * step + (ip === 1 ? 1 : 2) < v.length; k++) { const i = o + k * step, id = ip === 1 ? add(v[i], v[i + 1], v[3]) : add(v[i], v[i + 1], v[i + 2]); if (prev >= 0 && d.form !== 1 && d.form !== 2 && d.form !== 3) seg.push(prev, id); prev = id; }
      meta.read.copious++;
    } else if (d.type === 126) {
      const K = v[1], Mm = v[2], A = K + 1 + Mm, o = 7 + (A + 1) + (K + 1); if (!(K >= 0 && Mm >= 0) || o + 3 * (K + 1) > v.length) continue;
      for (let k = 0; k <= K; k++) add(v[o + 3 * k], v[o + 3 * k + 1], v[o + 3 * k + 2]); meta.read.curveControlPoints += K + 1;
    } else if (d.type === 128) {
      const K1 = v[1], K2 = v[2], M1 = v[3], M2 = v[4], C = (K1 + 1) * (K2 + 1), o = 10 + (K1 + M1 + 2) + (K2 + M2 + 2) + C; if (!(K1 >= 0 && K2 >= 0 && M1 >= 0 && M2 >= 0) || o + 3 * C > v.length) continue;
      for (let k = 0; k < C; k++) add(v[o + 3 * k], v[o + 3 * k + 1], v[o + 3 * k + 2]); meta.read.surfaceControlPoints += C;
    }
  }
  if (flag === flag && !length) ctx.warn(`The IGES unit "${meta.unitsName ?? flag}" has no platform equivalent; assign the units manually.`);
  if (withXf) ctx.warn(`${withXf} of the entities read reference a transformation matrix (entity 124); the transformations are NOT applied.`);
  const surf = [114, 118, 120, 122, 128, 140, 143, 144, 186, 510].reduce((s, t) => s + (census.get(t) || 0), 0);
  ctx.warn(`${surf ? `Surfaces and solids (${surf} surface/face entities) are NOT tessellated. ` : ''}Points, lines, copious data and B-spline control points are shown as a point cloud / polylines; extents are indicative only (control points may lie off the geometry).`);
  return { positions: Float64Array.from(xyz), lines: Uint32Array.from(seg), kind: 'cad-brep', units: { length, source: length ? 'file' : null }, meta };
}

// ---------- Parasolid text transmit file ----------
export function readXT(b, ctx) {
  const head = str(b, 0, Math.min(b.length, 16384)), meta = { header: {}, modellerVersion: null, schema: null, format: 'text', bytes: b.length, entityCensus: null };
  const hEnd = head.indexOf('**END_OF_HEADER');
  for (const m of (hEnd < 0 ? head : head.slice(0, hEnd)).replace(/\r?\n/g, '').matchAll(/([A-Z_]+)=([^;]*);/g)) meta.header[m[1]] = m[2];
  let m = /TRANSMIT FILE[^\n]*?version\s+(\d+)/i.exec(head); if (m) meta.modellerVersion = m[1];
  m = /\bSCH_\w+/.exec(head); if (m) meta.schema = m[0];
  ctx.warn('Parasolid transmit file: the header was read, but the topology/geometry stream is schema-driven and is NOT decoded — no entity census, extents or tessellation are available.');
  return { kind: 'metadata-only', meta };
}

// ---------- header-only inspectors ----------
function sat(b, ctx) {
  const meta = { encoding: 'SAT (text)' }; let units = { length: null, source: null };
  const sc = new Scanner(b), l1 = nums(sc.line() ?? ''), l2 = sc.line() ?? '', l3 = nums(sc.line() ?? '');
  meta.saveVersion = l1[0] ?? null; meta.records = l1[1] ?? null; meta.bodies = l1[2] ?? null;
  const strs = []; for (let i = 0, k = 0; i < l2.length && k < 6; k++) { const m = /^\s*@?(\d+)\s/.exec(l2.slice(i)); if (!m) break; strs.push(l2.substr(i + m[0].length, +m[1])); i += m[0].length + +m[1]; }
  [meta.product, meta.acisVersion, meta.date] = [strs[0] ?? null, strs[1] ?? null, strs[2] ?? null];
  if (l3.length >= 2 && l3[0] > 0) { meta.millimetresPerUnit = l3[0]; meta.resabs = l3[1]; const u = { 1: 'mm', 10: 'cm', 1000: 'm', 25.4: 'in', 304.8: 'ft' }[l3[0]]; if (u) units = { length: u, source: 'file' }; }
  if (b.length <= 64e6) { const c = {}; for (const m of str(b, 0, b.length).matchAll(/^(?:-\d+ )?([a-z][\w-]*) \$/gm)) c[m[1]] = (c[m[1]] || 0) + 1; meta.census = c; }
  return { kind: 'metadata-only', meta, units };
}
function tiff(b) {
  const le = b[0] === 0x49, dv = view(b), big = dv.getUint16(2, le) === 43, meta = { byteOrder: le ? 'little-endian' : 'big-endian', bigTiff: big, geoTiff: false };
  if (big) return { kind: 'metadata-only', meta };
  const ifd = dv.getUint32(4, le); if (ifd + 2 > b.length) return { kind: 'metadata-only', meta };
  const n = Math.min(dv.getUint16(ifd, le), 4096), SZ = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 12: 8 };
  const vals = (type, cnt, at, max = 64) => {
    const s = SZ[type]; if (!s) return []; let o = cnt * s <= 4 ? at : dv.getUint32(at, le); const out = [];
    for (let i = 0; i < Math.min(cnt, max) && o + s <= b.length; i++, o += s) out.push(type === 3 ? dv.getUint16(o, le) : type === 4 ? dv.getUint32(o, le) : type === 12 ? dv.getFloat64(o, le) : dv.getUint8(o));
    return out;
  };
  for (let i = 0; i < n && ifd + 2 + 12 * i + 12 <= b.length; i++) {
    const e = ifd + 2 + 12 * i, tag = dv.getUint16(e, le), type = dv.getUint16(e + 2, le), cnt = dv.getUint32(e + 4, le), v = vals(type, cnt, e + 8);
    if (tag === 256) meta.width = v[0]; else if (tag === 257) meta.height = v[0]; else if (tag === 258) meta.bitsPerSample = v; else if (tag === 259) meta.compression = v[0]; else if (tag === 277) meta.samplesPerPixel = v[0]; else if (tag === 339) meta.sampleFormat = v[0];
    else if (tag === 33550) { meta.pixelScale = v.slice(0, 3); meta.geoTiff = true; } else if (tag === 33922) { meta.tiePoints = v.slice(0, 6); meta.geoTiff = true; }
    else if (tag === 34735) {
      meta.geoTiff = true; const keys = {}; const NAMES = { 1024: 'modelType', 1025: 'rasterType', 2048: 'geographicCRS', 3072: 'projectedCRS', 3076: 'linearUnits', 4096: 'verticalCRS', 4099: 'verticalUnits' };
      for (let k = 4; k + 3 < v.length; k += 4) if (v[k + 1] === 0 && NAMES[v[k]]) keys[NAMES[v[k]]] = v[k + 3];
      meta.geoKeys = keys;
    }
  }
  return { kind: 'metadata-only', meta };
}
/** Inspect a recognised-but-undecoded file: report what its header states and nothing more. */
export function readMetadataOnly(id) {
  return (b, ctx) => {
    const dv = view(b), head = str(b, 0, Math.min(b.length, 512));
    if (id === 'acis') return head.startsWith('ACIS BinaryFile') ? { kind: 'metadata-only', meta: { encoding: 'SAB (binary)' } } : sat(b, ctx);
    if (id === 'geotiff') return tiff(b);
    if (id === 'laz') { const h = lasHeader(b); return { kind: 'metadata-only', meta: { ...h, compressed: true } }; }
    const meta = {};
    if (['cgns', 'exodus', 'med', 'fluent-h5', 'hdf5'].includes(id)) {
      if (b[0] === 0x89) { meta.container = 'HDF5'; meta.superblockVersion = b[8]; }
      else if (head.startsWith('CDF')) { meta.container = 'NetCDF classic'; meta.netcdfVersion = b[3]; if (b.length >= 8) meta.records = dv.getUint32(4, false); }
      else if (head.startsWith('@(#)ADF')) { meta.container = 'ADF'; meta.adfVersion = /ADF Database Version\s+(\S+)/.exec(head)?.[1] ?? null; }
      else meta.container = 'unconfirmed (extension only)';
    } else if (id === 'jt') meta.version = /^Version\s+([\d.]+)\s+JT/.exec(head)?.[1] ?? null;
    else if (id === 'e57' && b.length >= 48) { meta.version = `${dv.getUint32(8, true)}.${dv.getUint32(12, true)}`; meta.fileLength = Number(dv.getBigUint64(16, true)); meta.xmlOffset = Number(dv.getBigUint64(24, true)); meta.xmlLength = Number(dv.getBigUint64(32, true)); meta.pageSize = Number(dv.getBigUint64(40, true)); }
    else if (id === 'tecplot-bin') { meta.magic = head.slice(0, 8).replace(/[^\x20-\x7e]/g, ''); meta.flavour = head.startsWith('#!SZPLT') ? 'SZL (.szplt)' : 'classic binary (.plt)'; }
    else if (id === 'xb') meta.signature = 'Parasolid binary transmit';
    else if (id === 'nativecad') meta.container = b[0] === 0xd0 && b[1] === 0xcf ? 'OLE compound document' : b[0] === 0x50 && b[1] === 0x4b ? 'ZIP package' : head.startsWith('V5_CFV2') ? 'CATIA V5' : 'proprietary';
    return { kind: 'metadata-only', meta };
  };
}
