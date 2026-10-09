// AutoCAD DXF reader (ASCII and binary), written from the published DXF reference. A DXF file is a flat sequence of
// (group code, value) pairs organised in sections. Read here: HEADER ($ACADVER, $INSUNITS), TABLES (layers), BLOCKS,
// ENTITIES and ACDSDATA. Faces come from 3DFACE, SOLID/TRACE, MESH and polyface / polygon-mesh POLYLINEs; wires from
// LINE, LWPOLYLINE, POLYLINE, CIRCLE, ARC, ELLIPSE and SPLINE; INSERTs are expanded with their block transforms.
// 3DSOLID / BODY / REGION / SURFACE entities carry ACIS data - obfuscated SAT text in the entity (R2000-R2010) or
// binary SAB in the ACDSDATA section (R2013+) - which is returned for the ACIS translator.

import { str, view, Mesh, guard, LIMITS, M4 } from './parsers-util.js';

const BIN_SENTINEL = 'AutoCAD Binary DXF\r\n\u001a\u0000';
const ACIS_TYPES = new Set(['3DSOLID', 'BODY', 'REGION', 'SURFACE', 'PLANESURFACE', 'EXTRUDEDSURFACE', 'REVOLVEDSURFACE', 'SWEPTSURFACE', 'LOFTEDSURFACE', 'NURBSSURFACE']);
const INSUNITS = { 0: null, 1: 'in', 2: 'ft', 4: 'mm', 5: 'cm', 6: 'm' };
const INSUNIT_NAMES = { 3: 'miles', 7: 'kilometres', 8: 'microinches', 9: 'mils', 10: 'yards', 11: 'angstroms', 12: 'nanometres', 13: 'microns', 14: 'decimetres', 15: 'decametres', 16: 'hectometres', 17: 'gigametres', 18: 'astronomical units', 19: 'light years', 20: 'parsecs', 21: 'US survey feet' };

/** Kind of value a binary-DXF group code carries. */
function binKind(c) {
  if (c <= 9 || (c >= 100 && c <= 105) || (c >= 300 && c <= 309) || (c >= 320 && c <= 369) || (c >= 390 && c <= 399) || (c >= 410 && c <= 419) || (c >= 430 && c <= 439) || (c >= 470 && c <= 481) || c === 999 || (c >= 1000 && c <= 1003) || (c >= 1005 && c <= 1009)) return 's';
  if ((c >= 10 && c <= 59) || (c >= 110 && c <= 149) || (c >= 210 && c <= 239) || (c >= 460 && c <= 469) || (c >= 1010 && c <= 1059)) return 'f';
  if ((c >= 60 && c <= 79) || (c >= 170 && c <= 179) || (c >= 270 && c <= 289) || (c >= 370 && c <= 389) || (c >= 400 && c <= 409) || (c >= 1060 && c <= 1070)) return 'h';
  if ((c >= 90 && c <= 99) || (c >= 420 && c <= 429) || (c >= 440 && c <= 459) || c === 1071) return 'i';
  if (c >= 160 && c <= 169) return 'q';
  if (c >= 290 && c <= 299) return 'b';
  if ((c >= 310 && c <= 319) || c === 1004) return 'x';
  return null;
}

/** Split the file into parallel arrays of group codes and values (strings; numbers and byte chunks in binary DXF). */
function dxfPairs(b) {
  const codes = [], vals = [];
  if (str(b, 0, 22) === BIN_SENTINEL) {
    const dv = view(b), n = b.length, wide = b[23] === 0; let p = 22;
    while (p < n) {
      let c;
      if (wide) { if (p + 2 > n) break; c = dv.getUint16(p, true); p += 2; } else { c = b[p++]; if (c === 255) { if (p + 2 > n) break; c = dv.getUint16(p, true); p += 2; } }
      const k = binKind(c); if (!k) throw new Error(`binary DXF: unknown group code ${c} at offset ${p}`);
      let v;
      if (k === 's') { let e = p; while (e < n && b[e] !== 0) e++; v = str(b, p, e); p = e + 1; }
      else if (k === 'f') { if (p + 8 > n) break; v = dv.getFloat64(p, true); p += 8; }
      else if (k === 'h') { if (p + 2 > n) break; v = dv.getInt16(p, true); p += 2; }
      else if (k === 'i') { if (p + 4 > n) break; v = dv.getInt32(p, true); p += 4; }
      else if (k === 'q') { if (p + 8 > n) break; v = Number(dv.getBigInt64(p, true)); p += 8; }
      else if (k === 'b') { v = b[p++]; }
      else { const len = b[p++]; if (p + len > n) break; v = b.subarray(p, p + len); p += len; }
      codes.push(c); vals.push(v); guard(codes.length, 'DXF group', 2e8);
    }
    return { codes, vals, binary: true };
  }
  // ASCII: a group code line followed by a value line
  const s = new TextDecoder('utf-8').decode(b[0] === 0xef && b[1] === 0xbb && b[2] === 0xbf ? b.subarray(3) : b), n = s.length; let p = 0, bad = 0;
  const line = () => { let e = s.indexOf('\n', p); if (e < 0) e = n; let q = e; if (q > p && s.charCodeAt(q - 1) === 13) q--; const t = s.slice(p, q); p = e + 1; return t; };
  while (p < n) {
    const cl = line(); if (p >= n && !cl.trim()) break;
    const c = parseInt(cl, 10);
    if (!(c >= -5 && c < 1100)) { if (++bad > 50) throw new Error('the file is not a sequence of DXF group codes and values'); continue; }
    if (p > n) break;
    codes.push(c); vals.push(line()); guard(codes.length, 'DXF group', 2e8);
  }
  if (bad && codes.length < 4) throw new Error('the file is not a sequence of DXF group codes and values');
  return { codes, vals, binary: false };
}

const num = (v) => (typeof v === 'number' ? v : +String(v).trim() || 0);
const txt = (v) => (typeof v === 'string' ? v.trim() : String(v ?? ''));

/** Arbitrary-axis algorithm: the object coordinate system of an extrusion direction, as a matrix. */
function ocs(nx, ny, nz) {
  const l = Math.hypot(nx, ny, nz); if (!(l > 0) || (Math.abs(nx) < 1e-12 && Math.abs(ny) < 1e-12 && nz > 0)) return null;
  nx /= l; ny /= l; nz /= l;
  let ax, ay, az; if (Math.abs(nx) < 1 / 64 && Math.abs(ny) < 1 / 64) { ax = nz; ay = 0; az = -nx; } else { ax = -ny; ay = nx; az = 0; }
  const al = Math.hypot(ax, ay, az); ax /= al; ay /= al; az /= al;
  return [ax, ay, az, 0, ny * az - nz * ay, nz * ax - nx * az, nx * ay - ny * ax, 0, nx, ny, nz, 0, 0, 0, 0, 1];
}
/** The obfuscated SAT text of DXF R2000-R2010: every character but the blank is mirrored about 79.5, and "^ " stands for "A". */
export function dxfDecodeSat(line) {
  let o = '';
  for (let i = 0; i < line.length; i++) { const c = line.charCodeAt(i); if (c === 32) o += ' '; else if (c === 94 && line.charCodeAt(i + 1) === 32) { o += 'A'; i++; } else o += c > 32 && c < 127 ? String.fromCharCode(159 - c) : line[i]; }
  return o;
}
/** Inverse of dxfDecodeSat (used to build test fixtures and nothing else in the library). */
export function dxfEncodeSat(line) {
  let o = '';
  for (let i = 0; i < line.length; i++) { const c = line.charCodeAt(i); if (c === 32) o += ' '; else if (c === 65) o += '^ '; else o += c > 32 && c < 127 ? String.fromCharCode(159 - c) : line[i]; }
  return o;
}

/**
 * Read a DXF file. Returns { part, acis, meta, warnings }: `part` is the mesh part (faces, wires, layer groups),
 * `acis` the ACIS solids still to be translated: [{ type, layer, handle, sat | sab, matrix }].
 */
export function readDXFNative(b, opts = {}) {
  const { codes, vals, binary } = dxfPairs(b), N = codes.length, warnings = [];
  const meta = { encoding: binary ? 'binary' : 'ASCII', version: null, insunits: null, layers: [], blocks: 0, entityCensus: {}, skippedEntities: {}, acisSolids: 0 };
  // ---- sections ----
  const sections = new Map();
  for (let i = 0; i < N; i++) if (codes[i] === 0 && txt(vals[i]) === 'SECTION' && codes[i + 1] === 2) { const name = txt(vals[i + 1]).toUpperCase(); let e = i + 2; while (e < N && !(codes[e] === 0 && txt(vals[e]) === 'ENDSEC')) e++; sections.set(name, [i + 2, e]); i = e; }
  if (!sections.size) {
    if (!codes.some((c, i) => c === 0 && /^(LINE|3DFACE|POLYLINE|LWPOLYLINE|CIRCLE|ARC|INSERT)$/.test(txt(vals[i])))) throw new Error('no DXF sections or entities were found');
    sections.set('ENTITIES', [0, N]);                         // a bare entity list without section markers
  }
  const records = (s, e) => { const out = []; let i = s; while (i < e && codes[i] !== 0) i++; while (i < e) { let j = i + 1; while (j < e && codes[j] !== 0) j++; out.push([i, j]); i = j; } return out; };
  const first = (r, code, dflt = undefined) => { for (let i = r[0] + 1; i < r[1]; i++) if (codes[i] === code) return vals[i]; return dflt; };
  const typeOf = (r) => txt(vals[r[0]]).toUpperCase();
  // ---- header ----
  let ext = null;
  if (sections.has('HEADER')) {
    const [s, e] = sections.get('HEADER');
    for (let i = s; i < e; i++) {
      if (codes[i] !== 9) continue; const name = txt(vals[i]);
      if (name === '$ACADVER') meta.version = txt(vals[i + 1]);
      else if (name === '$INSUNITS') meta.insunits = num(vals[i + 1]);
      else if (name === '$MEASUREMENT') meta.measurement = num(vals[i + 1]);
      else if (name === '$EXTMIN' || name === '$EXTMAX') { (ext ??= {})[name] = [num(vals[i + 1]), num(vals[i + 2]), num(vals[i + 3])]; }
    }
  }
  const RELEASE = { AC1006: 'R10', AC1009: 'R11/R12', AC1012: 'R13', AC1014: 'R14', AC1015: '2000', AC1018: '2004', AC1021: '2007', AC1024: '2010', AC1027: '2013', AC1032: '2018' };
  meta.release = RELEASE[meta.version] ?? null;
  // ---- layers ----
  if (sections.has('TABLES')) { const [s, e] = sections.get('TABLES'); for (const r of records(s, e)) if (typeOf(r) === 'LAYER') { const nm = first(r, 2); if (nm !== undefined && meta.layers.length < 5000) meta.layers.push(txt(nm)); } }
  // ---- blocks ----
  const blocks = new Map();
  if (sections.has('BLOCKS')) {
    const [s, e] = sections.get('BLOCKS'); let cur = null;
    for (const r of records(s, e)) {
      const t = typeOf(r);
      if (t === 'BLOCK') { cur = { name: txt(first(r, 2, '')), base: [num(first(r, 10, 0)), num(first(r, 20, 0)), num(first(r, 30, 0))], ents: [] }; blocks.set(cur.name.toUpperCase(), cur); }
      else if (t === 'ENDBLK') cur = null;
      else if (cur) cur.ents.push(r);
    }
    meta.blocks = blocks.size;
  }
  // ---- binary ACIS data of R2013+ (ACDSDATA): owner handle → SAB bytes ----
  const sabOf = new Map();
  if (sections.has('ACDSDATA')) {
    const [s, e] = sections.get('ACDSDATA');
    for (const r of records(s, e)) {
      if (typeOf(r) !== 'ACDSRECORD') continue;
      let handle = null, chunks = [], size = 0, asm = false;
      for (let i = r[0] + 1; i < r[1]; i++) {
        const c = codes[i], v = vals[i];
        if (c === 320) handle = txt(v).toUpperCase();
        else if (c === 2) asm = txt(v) === 'ASM_Data';
        else if (c === 310 && asm) { let u; if (typeof v === 'string') { const h = v.trim(); u = new Uint8Array(h.length >> 1); for (let k = 0; k < u.length; k++) u[k] = parseInt(h.substr(2 * k, 2), 16); } else u = v; chunks.push(u); size += u.length; guard(size, 'ACIS byte', 400e6); }
      }
      if (handle && size) { const u = new Uint8Array(size); let o = 0; for (const c of chunks) { u.set(c, o); o += c.length; } sabOf.set(handle, u); }
    }
  }
  // ---- entities ----
  const M = new Mesh(), acis = [], census = meta.entityCensus, skipped = meta.skippedEntities;
  const layerGroup = (name) => M.group('layer', name, name);
  const put = (m, x, y, z) => { const q = m ? M4.apply(m, x, y, z) : [x, y, z]; return M.node(q[0], q[1], q[2]); };
  const polyline = (pts, m, g, closed) => { if (pts.length < 2) return; guard(M.nNodes + pts.length, 'DXF vertex', LIMITS.verts); const ids = pts.map((q) => put(m, q[0], q[1], q[2] || 0)); for (let i = 0; i + 1 < ids.length; i++) M.elem('line2', [ids[i], ids[i + 1]], g); if (closed && ids.length > 2) M.elem('line2', [ids[ids.length - 1], ids[0]], g); };
  const arcPts = (cx, cy, z, r, a0, a1) => { let sw = a1 - a0; while (sw <= 1e-12) sw += 2 * Math.PI; const n = Math.max(2, Math.min(720, Math.ceil(sw / (Math.PI / 36)))), out = []; for (let i = 0; i <= n; i++) { const a = a0 + (sw * i) / n; out.push([cx + r * Math.cos(a), cy + r * Math.sin(a), z]); } return out; };
  /** 2-D vertices with bulges → points along the arcs. */
  const bulged = (v, closed) => {
    const out = [], n = v.length;
    for (let i = 0; i < n; i++) {
      const a = v[i], bq = v[(i + 1) % n]; out.push([a[0], a[1], a[2]]);
      if (!a[3] || (i === n - 1 && !closed)) continue;
      const bu = a[3], dx = bq[0] - a[0], dy = bq[1] - a[1], c = Math.hypot(dx, dy); if (!(c > 0)) continue;
      const k = (1 - bu * bu) / (4 * bu), cx = (a[0] + bq[0]) / 2 - dy * k, cy = (a[1] + bq[1]) / 2 + dx * k, r = Math.hypot(a[0] - cx, a[1] - cy), a0 = Math.atan2(a[1] - cy, a[0] - cx), th = 4 * Math.atan(bu), steps = Math.max(1, Math.min(360, Math.ceil(Math.abs(th) / (Math.PI / 36))));
      for (let s = 1; s < steps; s++) { const an = a0 + (th * s) / steps; out.push([cx + r * Math.cos(an), cy + r * Math.sin(an), a[2]]); }
    }
    return out;
  };
  const extrusion = (r) => ocs(num(first(r, 210, 0)), num(first(r, 220, 0)), num(first(r, 230, 1)));
  const withOcs = (m, o) => (o ? (m ? M4.mul(m, o) : o) : m);
  const bump = (o, k) => { o[k] = (o[k] || 0) + 1; };
  let inserts = 0;
  const emitList = (list, m, depth, layer0) => {
    for (let q = 0; q < list.length; q++) {
      const r = list[q], t = typeOf(r); if (t === 'SEQEND' || t === 'VERTEX' || t === 'ATTRIB') continue;
      if (num(first(r, 67, 0)) === 1) continue;               // paper space
      let layer = txt(first(r, 8, '0')); if (layer === '0' && layer0) layer = layer0;      // entities on layer 0 of a block take the layer of the INSERT
      bump(census, t);
      const g = () => layerGroup(layer);
      const P = (k) => [num(first(r, 10 + k, 0)), num(first(r, 20 + k, 0)), num(first(r, 30 + k, 0))];
      switch (t) {
        case '3DFACE': case 'SOLID': case 'TRACE': {
          const c = [P(0), P(1), P(2), first(r, 13) === undefined ? P(2) : P(3)], mm = t === '3DFACE' ? m : withOcs(m, extrusion(r));
          if (t !== '3DFACE') { const sw = c[2]; c[2] = c[3]; c[3] = sw; }             // SOLID / TRACE list their corners in a Z pattern
          const same = (a, d) => a[0] === d[0] && a[1] === d[1] && a[2] === d[2], ids = c.map((p) => put(mm, p[0], p[1], p[2]));
          if (same(c[2], c[3])) M.elem('tri3', [ids[0], ids[1], ids[2]], g()); else M.elem('quad4', ids, g());
          break;
        }
        case 'LINE': polyline([P(0), P(1)], m, g(), false); break;
        case 'CIRCLE': case 'ARC': {
          const c = P(0), rad = num(first(r, 40, 0)), a0 = t === 'ARC' ? (num(first(r, 50, 0)) * Math.PI) / 180 : 0, a1 = t === 'ARC' ? (num(first(r, 51, 360)) * Math.PI) / 180 : 2 * Math.PI;
          if (rad > 0) polyline(arcPts(c[0], c[1], c[2], rad, a0, a1), withOcs(m, extrusion(r)), g(), false);
          break;
        }
        case 'ELLIPSE': {
          const c = P(0), a = P(1), ratio = num(first(r, 40, 1)), u0 = num(first(r, 41, 0)), u1raw = num(first(r, 42, 2 * Math.PI)), nrm = [num(first(r, 210, 0)), num(first(r, 220, 0)), num(first(r, 230, 1))];
          const bx = [nrm[1] * a[2] - nrm[2] * a[1], nrm[2] * a[0] - nrm[0] * a[2], nrm[0] * a[1] - nrm[1] * a[0]].map((x) => x * ratio); let sw = u1raw - u0; while (sw <= 1e-12) sw += 2 * Math.PI;
          const n = Math.max(4, Math.ceil(sw / (Math.PI / 36))), pts = []; for (let i = 0; i <= n; i++) { const u = u0 + (sw * i) / n, cu = Math.cos(u), su = Math.sin(u); pts.push([c[0] + a[0] * cu + bx[0] * su, c[1] + a[1] * cu + bx[1] * su, c[2] + a[2] * cu + bx[2] * su]); }
          polyline(pts, m, g(), false); break;
        }
        case 'LWPOLYLINE': {
          const v = [], elev = num(first(r, 38, 0)); for (let i = r[0] + 1; i < r[1]; i++) { if (codes[i] === 10) v.push([num(vals[i]), 0, elev, 0]); else if (codes[i] === 20 && v.length) v[v.length - 1][1] = num(vals[i]); else if (codes[i] === 42 && v.length) v[v.length - 1][3] = num(vals[i]); }
          const closed = (num(first(r, 70, 0)) & 1) === 1; polyline(bulged(v, closed), withOcs(m, extrusion(r)), g(), closed); break;
        }
        case 'POLYLINE': {
          const flags = num(first(r, 70, 0)), verts = []; let k = q + 1;
          for (; k < list.length && typeOf(list[k]) === 'VERTEX'; k++) verts.push(list[k]);
          if (flags & 64) {                                    // polyface mesh: coordinate vertices, then face records with 1-based indices
            const ids = [], gg = g();
            for (const vr of verts) {
              const vf = num(first(vr, 70, 0));
              if ((vf & 192) === 192) ids.push(put(m, num(first(vr, 10, 0)), num(first(vr, 20, 0)), num(first(vr, 30, 0))));
              else if (vf & 128) { const f = [71, 72, 73, 74].map((c) => Math.abs(num(first(vr, c, 0)))).filter((x) => x > 0).map((x) => ids[x - 1]); if (f.length === 3) M.elem('tri3', f, gg); else if (f.length === 4) M.elem('quad4', f, gg); else if (f.length === 2) M.elem('line2', f, gg); }
            }
          } else if (flags & 16) {                             // polygon mesh, M × N vertices
            const mN = num(first(r, 71, 0)), nN = num(first(r, 72, 0)), gg = g();
            if (mN >= 2 && nN >= 2 && verts.length >= mN * nN) {
              const ids = verts.slice(0, mN * nN).map((vr) => put(m, num(first(vr, 10, 0)), num(first(vr, 20, 0)), num(first(vr, 30, 0)))), mc = flags & 1 ? mN : mN - 1, nc = flags & 32 ? nN : nN - 1;
              for (let i = 0; i < mc; i++) for (let j = 0; j < nc; j++) M.elem('quad4', [ids[i * nN + j], ids[((i + 1) % mN) * nN + j], ids[((i + 1) % mN) * nN + ((j + 1) % nN)], ids[i * nN + ((j + 1) % nN)]], gg);
            } else warnings.push('A polygon-mesh POLYLINE has fewer vertices than its M × N counts state and was skipped.');
          } else {
            const is3d = (flags & 8) !== 0, elev = num(first(r, 30, 0)), v = verts.map((vr) => [num(first(vr, 10, 0)), num(first(vr, 20, 0)), is3d ? num(first(vr, 30, 0)) : elev, is3d ? 0 : num(first(vr, 42, 0))]), closed = (flags & 1) === 1;
            polyline(is3d ? v : bulged(v, closed), is3d ? m : withOcs(m, extrusion(r)), g(), closed);
          }
          q = k - 1; break;
        }
        case 'MESH': {
          const pts = [], faces = []; let mode = 0, lvl = 0;
          for (let i = r[0] + 1; i < r[1]; i++) {
            const c = codes[i], v = vals[i];
            if (c === 91) lvl = num(v); else if (c === 92) mode = 1; else if (c === 93) mode = 2; else if (c === 94 || c === 95) mode = 3;
            else if (c === 10 && mode === 1) pts.push([num(v), 0, 0]); else if (c === 20 && mode === 1 && pts.length) pts[pts.length - 1][1] = num(v); else if (c === 30 && mode === 1 && pts.length) pts[pts.length - 1][2] = num(v);
            else if (c === 90 && mode === 2) faces.push(num(v));
          }
          const base = M.nNodes, gg = g(); for (const p of pts) put(m, p[0], p[1], p[2]);
          for (let i = 0; i < faces.length;) { const n = faces[i++]; if (!(n >= 3 && i + n <= faces.length)) break; const f = faces.slice(i, i + n); i += n; if (f.some((x) => !(x >= 0 && x < pts.length))) continue; if (n === 3) M.elem('tri3', f.map((x) => base + x), gg); else if (n === 4) M.elem('quad4', f.map((x) => base + x), gg); else for (let k = 1; k + 1 < n; k++) M.elem('tri3', [base + f[0], base + f[k], base + f[k + 1]], gg); }
          if (lvl > 0) warnings.push(`A MESH entity has subdivision level ${lvl}; its control mesh is shown, not the smoothed surface.`);
          break;
        }
        case 'SPLINE': {
          const ctrl = [], fit = [], knots = []; const deg = num(first(r, 71, 3));
          for (let i = r[0] + 1; i < r[1]; i++) { const c = codes[i], v = num(vals[i]); if (c === 10) ctrl.push([v, 0, 0, 1]); else if (c === 20 && ctrl.length) ctrl[ctrl.length - 1][1] = v; else if (c === 30 && ctrl.length) ctrl[ctrl.length - 1][2] = v; else if (c === 41 && ctrl.length) ctrl[ctrl.length - 1][3] = v || 1; else if (c === 11) fit.push([v, 0, 0]); else if (c === 21 && fit.length) fit[fit.length - 1][1] = v; else if (c === 31 && fit.length) fit[fit.length - 1][2] = v; else if (c === 40) knots.push(v); }
          let pts = null;
          if (ctrl.length > deg && knots.length === ctrl.length + deg + 1 && deg >= 1 && deg <= 12) {
            const n = ctrl.length, u0 = knots[deg], u1 = knots[n], steps = Math.min(2000, Math.max(8, 8 * (n - deg))); pts = [];
            for (let s = 0; s <= steps; s++) {
              const u = Math.min(u1, u0 + ((u1 - u0) * s) / steps); let k = deg; while (k < n - 1 && knots[k + 1] <= u) k++;
              const d = []; for (let j = 0; j <= deg; j++) { const c = ctrl[k - deg + j]; d.push([c[0] * c[3], c[1] * c[3], c[2] * c[3], c[3]]); }
              for (let rr = 1; rr <= deg; rr++) for (let j = deg; j >= rr; j--) { const i0 = k - deg + j, den = knots[i0 + deg - rr + 1] - knots[i0], a = den > 0 ? (u - knots[i0]) / den : 0; for (let c = 0; c < 4; c++) d[j][c] = (1 - a) * d[j - 1][c] + a * d[j][c]; }
              const w = d[deg][3] || 1; pts.push([d[deg][0] / w, d[deg][1] / w, d[deg][2] / w]);
            }
          } else if (fit.length >= 2) { pts = fit; warnings.push('A SPLINE without usable control data is drawn as the polyline through its fit points.'); }
          else if (ctrl.length >= 2) { pts = ctrl; warnings.push('A SPLINE with inconsistent knots is drawn as its control polygon.'); }
          if (pts) polyline(pts, m, g(), false);
          break;
        }
        case 'INSERT': {
          const blk = blocks.get(txt(first(r, 2, '')).toUpperCase());
          if (!blk) { bump(skipped, 'INSERT of an undefined block'); break; }
          if (depth >= 16) { bump(skipped, 'INSERT nested deeper than 16 levels'); break; }
          if (++inserts > 200000) { bump(skipped, 'INSERT beyond the expansion limit'); break; }
          const ins = P(0), sx = num(first(r, 41, 1)) || 1, sy = num(first(r, 42, 1)) || 1, sz = num(first(r, 43, 1)) || 1, rot = (num(first(r, 50, 0)) * Math.PI) / 180, cr = Math.cos(rot), sr = Math.sin(rot);
          const cols = Math.max(1, num(first(r, 70, 1))), rows = Math.max(1, num(first(r, 71, 1))), dc = num(first(r, 44, 0)), dr = num(first(r, 45, 0)), o = extrusion(r);
          for (let ci = 0; ci < Math.min(cols, 1000); ci++) for (let ri = 0; ri < Math.min(rows, 1000); ri++) {
            // p' = OCS · ( T(ins) · R(rot) · ( T(col, row offset) · S · (p − base) ) )
            const ox = ci * dc, oy = ri * dr, tx = ins[0] + cr * ox - sr * oy, ty = ins[1] + sr * ox + cr * oy;
            let loc = [cr * sx, sr * sx, 0, 0, -sr * sy, cr * sy, 0, 0, 0, 0, sz, 0, tx, ty, ins[2], 1];
            loc = M4.mul(loc, [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -blk.base[0], -blk.base[1], -blk.base[2], 1]);
            emitList(blk.ents, withOcs(m, o) ? M4.mul(withOcs(m, o), loc) : loc, depth + 1, layer);
          }
          break;
        }
        default:
          if (ACIS_TYPES.has(t)) {
            const handle = txt(first(r, 5, '')).toUpperCase(); let sat = null; const sab = sabOf.get(handle) || null;
            if (!sab) { const lines = []; for (let i = r[0] + 1; i < r[1]; i++) { if (codes[i] === 1) lines.push(String(vals[i])); else if (codes[i] === 3 && lines.length) lines[lines.length - 1] += String(vals[i]); } if (lines.length) sat = lines.map(dxfDecodeSat).join('\n') + '\n'; }
            if (sab || sat) { if (acis.length < 5000) acis.push({ type: t, layer, handle, sat, sab, matrix: m }); meta.acisSolids++; } else bump(skipped, `${t} without ACIS data`);
          } else bump(skipped, t);
      }
    }
  };
  if (sections.has('ENTITIES')) { const [s, e] = sections.get('ENTITIES'); emitList(records(s, e), null, 0, null); }
  const part = M.nNodes ? M.result({}) : { positions: new Float64Array(0), elements: [], groups: [] };
  const sk = Object.entries(skipped).filter(([k]) => !/^(TEXT|MTEXT|DIMENSION|HATCH|VIEWPORT|ATTDEF|LEADER|MLEADER|POINT|IMAGE|TOLERANCE|ACAD_PROXY_ENTITY)$/.test(k));
  const annot = Object.entries(skipped).filter(([k]) => /^(TEXT|MTEXT|DIMENSION|HATCH|VIEWPORT|ATTDEF|LEADER|MLEADER|POINT|IMAGE|TOLERANCE|ACAD_PROXY_ENTITY)$/.test(k)).reduce((s, [, n]) => s + n, 0);
  if (annot) warnings.push(`${annot} annotation entity(ies) (text, dimensions, hatches, points, viewports) carry no model geometry and were skipped.`);
  if (sk.length) warnings.push(`Entities not read: ${sk.map(([k, n]) => `${n} × ${k}`).join(', ')}.`);
  let length = null;
  if (meta.insunits !== null && meta.insunits !== undefined) { length = INSUNITS[meta.insunits] ?? null; if (!length && meta.insunits !== 0) warnings.push(`The drawing unit ($INSUNITS = ${meta.insunits}, ${INSUNIT_NAMES[meta.insunits] || 'unknown'}) has no platform equivalent; assign the units manually.`); }
  if (opts.extents !== false && ext && ext.$EXTMIN && ext.$EXTMAX) meta.headerExtents = { min: ext.$EXTMIN, max: ext.$EXTMAX };
  return { part, acis, meta, warnings, units: { length, source: length ? 'file' : null } };
}
