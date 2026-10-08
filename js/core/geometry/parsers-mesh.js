// Readers for analysis meshes: VTK legacy + XML, Gmsh, SU2, Nastran, Abaqus, UNV, ANSYS CDB, Fluent,
// Tecplot, Plot3D and OpenFOAM polyMesh. Each reader takes (bytes, ctx) and returns a partial model.

import { Scanner, Mesh, ET, guard, str, text, nums, toNum, xmlParse, xmlAll, xmlFirst, xmlChild, structuredBlock, LIMITS } from './parsers-util.js';

// VTK cell-type codes (also used by SU2).
const VTK_TYPES = { 3: 'line2', 5: 'tri3', 9: 'quad4', 10: 'tet4', 12: 'hex8', 13: 'wedge6', 14: 'pyr5', 21: 'line3', 22: 'tri6', 23: 'quad8', 24: 'tet10', 25: 'hex20', 26: 'wedge15', 27: 'pyr13', 28: 'quad9', 29: 'hex27' };
/** Collapse an 8-node brick with repeated nodes to the element it really is. */
function degen8(n) {
  const [a, b, c, d, e, f, g, h] = n;
  if (c === d && e === f && f === g && g === h) return ['tet4', [a, b, c, e]];
  if (e === f && f === g && g === h) return ['pyr5', [a, b, c, d, e]];
  if (c === d && g === h) return ['wedge6', [a, b, c, e, f, g]];
  return ['hex8', n];
}
/** Add one VTK cell (0-based node indices). Returns false for an unsupported type. */
function vtkCell(M, type, n, g = -1) {
  const t = VTK_TYPES[type];
  if (t) { if (n.length < ET[t].n) return false; M.elem(t, n, g); return true; }
  switch (type) {
    case 1: case 2: return true;                                                                    // vertices carry no connectivity
    case 4: for (let i = 1; i < n.length; i++) M.elem('line2', [n[i - 1], n[i]], g); return true;      // polyline
    case 6: for (let i = 2; i < n.length; i++) M.elem('tri3', i % 2 ? [n[i - 1], n[i - 2], n[i]] : [n[i - 2], n[i - 1], n[i]], g); return true; // strip
    case 7: for (let i = 2; i < n.length; i++) M.elem('tri3', [n[0], n[i - 1], n[i]], g); return true; // polygon (fan)
    case 8: if (n.length < 4) return false; M.elem('quad4', [n[0], n[1], n[3], n[2]], g); return true; // pixel
    case 11: if (n.length < 8) return false; M.elem('hex8', [n[0], n[1], n[3], n[2], n[4], n[5], n[7], n[6]], g); return true; // voxel
    default: return false;
  }
}

// ---------- VTK legacy ASCII ----------
export function readVTK(b, ctx) {
  const sc = new Scanner(b), version = (sc.line() ?? '').replace(/^#\s*/, '').trim(), title = (sc.line() ?? '').trim(), enc = (sc.line() ?? '').trim().toUpperCase();
  const meta = { version, title, encoding: enc, dataset: null, pointData: [], cellData: [] };
  if (enc !== 'ASCII') { const m = /DATASET\s+(\w+)/.exec(str(b, sc.p, sc.p + 200)); meta.dataset = m ? m[1] : null; throw Object.assign(new Error('VTK legacy BINARY encoding is not read; re-save as ASCII'), { meta }); }
  const M = new Mesh(); let dims = null, cells = null, unsupported = 0, where = null;
  const readCells = (a, c) => {
    guard(a, 'VTK cell'); guard(c, 'VTK connectivity', LIMITS.elems * 8);
    const off = [], conn = [];
    if (sc.peek() === 'OFFSETS') {                       // VTK 5.x layout: offsets array then connectivity array
      sc.token(); sc.token(); for (let i = 0; i < a; i++) off.push(sc.int('VTK offset'));
      sc.token(); sc.token(); for (let i = 0; i < c; i++) conn.push(sc.int('VTK connectivity'));
    } else {
      for (let i = 0; i < a; i++) { const k = sc.int('VTK cell size'); if (!(k >= 0 && k <= 1e6)) throw new Error('VTK cell size is not plausible'); off.push(conn.length); for (let j = 0; j < k; j++) conn.push(sc.int('VTK connectivity')); }
      off.push(conn.length);
    }
    return { off, conn };
  };
  const emit = (cl, typeOf) => { for (let i = 0; i + 1 < cl.off.length; i++) if (!vtkCell(M, typeOf(i), cl.conn.slice(cl.off[i], cl.off[i + 1]))) unsupported++; };
  for (let t; (t = sc.token()) !== null;) {
    switch (t.toUpperCase()) {
      case 'DATASET': meta.dataset = sc.token(); if (['STRUCTURED_POINTS', 'RECTILINEAR_GRID', 'FIELD'].includes(meta.dataset)) throw Object.assign(new Error(`VTK dataset type ${meta.dataset} is not read`), { meta }); break;
      case 'DIMENSIONS': dims = [sc.int(), sc.int(), sc.int()]; break;
      case 'POINTS': { const n = guard(sc.int('VTK point count'), 'VTK point', LIMITS.verts); sc.token(); for (let i = 0; i < n; i++) M.node(sc.float(), sc.float(), sc.float()); break; }
      case 'POLYGONS': emit(readCells(sc.int(), sc.int()), () => 7); break;
      case 'TRIANGLE_STRIPS': emit(readCells(sc.int(), sc.int()), () => 6); break;
      case 'LINES': emit(readCells(sc.int(), sc.int()), () => 4); break;
      case 'VERTICES': readCells(sc.int(), sc.int()); break;
      case 'CELLS': cells = readCells(sc.int(), sc.int()); break;
      case 'CELL_TYPES': { const n = guard(sc.int(), 'VTK cell type'), ty = []; for (let i = 0; i < n; i++) ty.push(sc.int('VTK cell type')); if (cells) emit(cells, (i) => ty[i]); cells = null; break; }
      case 'POINT_DATA': where = meta.pointData; sc.token(); break;
      case 'CELL_DATA': where = meta.cellData; sc.token(); break;
      case 'SCALARS': case 'VECTORS': case 'NORMALS': case 'TENSORS': case 'TEXTURE_COORDINATES': case 'COLOR_SCALARS': { const name = sc.token(); if (where && name) where.push(name); break; }
      case 'FIELD': {
        sc.token(); const na = sc.int('VTK field array count');
        for (let k = 0; k < na && k < 10000; k++) { const name = sc.token(), nc = sc.int(), nt = sc.int(); sc.token(); if (where && name) where.push(name); guard(nc * nt, 'VTK field value', LIMITS.elems * 8); for (let i = 0; i < nc * nt; i++) if (sc.token() === null) break; }
        break;
      }
      default:
    }
  }
  if (meta.dataset === 'STRUCTURED_GRID' && dims) { if (dims[0] * dims[1] * dims[2] !== M.nNodes) throw new Error('STRUCTURED_GRID DIMENSIONS do not match the point count'); structuredBlock(M, dims[0], dims[1], dims[2], 0, -1); meta.dimensions = dims; }
  if (unsupported) ctx.warn(`${unsupported} VTK cell(s) of unsupported type (e.g. polyhedra) were skipped.`);
  return M.result({ meta });
}

// ---------- VTK XML (.vtu / .vtp, inline ASCII only) ----------
export function readVTKXML(b, ctx) {
  const root = xmlFirst(xmlParse(text(b)), 'VTKFile');
  if (!root) throw new Error('no <VTKFile> root element');
  const type = root.attrs.type, pieces = xmlAll(root, 'Piece');
  const meta = { type, version: root.attrs.version ?? null, byteOrder: root.attrs.byte_order ?? null, compressor: root.attrs.compressor ?? null, pieces: pieces.map((p) => ({ ...p.attrs })), dataArrays: [...new Set(xmlAll(root, 'DataArray').map((d) => d.attrs.Name).filter(Boolean))] };
  if (type !== 'UnstructuredGrid' && type !== 'PolyData') throw Object.assign(new Error(`VTK XML type ${type} is recognised but not read (only UnstructuredGrid and PolyData)`), { meta });
  const M = new Mesh(); let unsupported = 0;
  const arr = (parent, name) => {
    const d = parent && (name ? xmlAll(parent, 'DataArray').find((x) => x.attrs.Name === name) : xmlFirst(parent, 'DataArray'));
    if (!d) return null;
    if ((d.attrs.format ?? 'ascii') !== 'ascii') throw Object.assign(new Error(`DataArray "${d.attrs.Name ?? ''}" is stored as ${d.attrs.format}; only inline ASCII DataArrays are read`), { meta });
    return nums(d.text);
  };
  for (const pc of pieces) {
    const base = M.nNodes, P = arr(xmlChild(pc, 'Points')); if (!P) continue;
    guard(base + P.length / 3, 'VTK point', LIMITS.verts);
    for (let i = 0; i + 2 < P.length; i += 3) M.node(P[i], P[i + 1], P[i + 2]);
    const block = (el, typeOf) => {
      const conn = arr(el, 'connectivity'), off = arr(el, 'offsets'); if (!conn || !off) return;
      const ty = typeOf === null ? arr(el, 'types') : null;
      for (let i = 0, s = 0; i < off.length; i++) { const e = off[i]; if (!(e >= s && e <= conn.length)) break; const n = conn.slice(s, e).map((v) => v + base); if (!vtkCell(M, ty ? ty[i] : typeOf, n)) unsupported++; s = e; }
    };
    if (type === 'UnstructuredGrid') block(xmlChild(pc, 'Cells'), null);
    else { block(xmlChild(pc, 'Polys'), 7); block(xmlChild(pc, 'Strips'), 6); block(xmlChild(pc, 'Lines'), 4); }
  }
  if (unsupported) ctx.warn(`${unsupported} VTK cell(s) of unsupported type were skipped.`);
  return M.result({ meta });
}

// ---------- Gmsh MSH 2.2 / 4.1 ASCII ----------
const GMSH = { 1: 'line2', 2: 'tri3', 3: 'quad4', 4: 'tet4', 5: 'hex8', 6: 'wedge6', 7: 'pyr5', 8: 'line3', 9: 'tri6', 10: 'quad9', 11: 'tet10', 12: 'hex27', 13: 'wedge18', 14: 'pyr14', 16: 'quad8', 17: 'hex20', 18: 'wedge15', 19: 'pyr13' };
export function readGmsh(b, ctx) {
  const sc = new Scanner(b), M = new Mesh(), names = new Map(), entPhys = new Map(), meta = { version: null, physicalNames: 0, pointElements: 0 };
  let major = 0, maxDim = 0;
  const phys = (dim, tag) => { if (!tag) return -1; const g = M.group('phys', dim + ':' + tag, names.get(dim + ':' + tag) ?? `physical ${dim}D #${tag}`); M.groups[g].dim = dim; M.groups[g].tag = tag; return g; };
  const add = (type, ids, dim, tag) => {
    if (type === 15) { meta.pointElements++; return; }
    const t = GMSH[type]; if (!t) throw new Error(`Gmsh element type ${type} is not supported (high-order or polyhedral)`);
    maxDim = Math.max(maxDim, ET[t].dim); M.elemIds(t, ids, phys(dim ?? ET[t].dim, tag));
  };
  const nodesOf = (type) => (type === 15 ? 1 : ET[GMSH[type]]?.n ?? (() => { throw new Error(`Gmsh element type ${type} is not supported (high-order or polyhedral)`); })());
  for (let ln; (ln = sc.line()) !== null;) {
    ln = ln.trim(); if (ln[0] !== '$' || ln.startsWith('$End')) continue;
    const sec = ln.slice(1);
    if (sec === 'MeshFormat') {
      const v = sc.token(), ft = sc.int('Gmsh file type'); meta.version = v; major = parseFloat(v);
      if (ft !== 0) throw Object.assign(new Error('binary Gmsh MSH is not read; re-save as ASCII (msh22 or msh41)'), { meta });
      if (!(major >= 2 && major < 3) && !(major >= 4.1 && major < 5)) throw Object.assign(new Error(`Gmsh MSH version ${v} is not read (supported: 2.x and 4.1 ASCII)`), { meta });
    } else if (sec === 'NOD' || sec === 'ELM') throw Object.assign(new Error('Gmsh MSH version 1 is not read; re-save as version 2.2 or 4.1'), { meta });
    else if (sec === 'PhysicalNames') {
      const n = guard(sc.int('physical name count'), 'physical name', 1e6); sc.line();
      for (let i = 0; i < n; i++) { const m = /^\s*(\d+)\s+(-?\d+)\s+"?(.*?)"?\s*$/.exec(sc.line() ?? ''); if (m) names.set(m[1] + ':' + m[2], m[3]); }
      meta.physicalNames = names.size;
    } else if (sec === 'Entities' && major >= 4) {
      const cnt = [sc.int(), sc.int(), sc.int(), sc.int()];
      for (let dim = 0; dim < 4; dim++) for (let i = 0; i < guard(cnt[dim], 'Gmsh entity', 1e7); i++) {
        const tag = sc.int('entity tag'); for (let k = 0; k < (dim === 0 ? 3 : 6); k++) sc.float('entity bounds');
        const np = sc.int('entity physical count'); for (let k = 0; k < np; k++) { const p = Math.abs(sc.int()); if (k === 0) entPhys.set(dim + ':' + tag, p); }
        if (dim > 0) { const nb = sc.int('entity boundary count'); for (let k = 0; k < nb; k++) sc.int(); }
      }
      meta.entities = { points: cnt[0], curves: cnt[1], surfaces: cnt[2], volumes: cnt[3] };
    } else if (sec === 'Nodes') {
      if (major < 4) { const n = guard(sc.int('node count'), 'Gmsh node', LIMITS.verts); for (let i = 0; i < n; i++) { const id = sc.int('node id'); M.node(sc.float(), sc.float(), sc.float(), id); } }
      else {
        const nb = guard(sc.int('node block count'), 'Gmsh node block'), total = guard(sc.int('node count'), 'Gmsh node', LIMITS.verts); sc.int(); sc.int();
        for (let k = 0, seen = 0; k < nb; k++) {
          const dim = sc.int(); sc.int(); const par = sc.int(), n = sc.int('nodes in block'); if (!(n >= 0) || (seen += n) > total) throw new Error('Gmsh node blocks exceed the declared node count');
          const tags = new Array(n); for (let i = 0; i < n; i++) tags[i] = sc.int('node tag');
          for (let i = 0; i < n; i++) { M.node(sc.float(), sc.float(), sc.float(), tags[i]); if (par) for (let j = 0; j < dim; j++) sc.float(); }
        }
      }
    } else if (sec === 'Elements') {
      if (major < 4) {
        const n = guard(sc.int('element count'), 'Gmsh element');
        for (let i = 0; i < n; i++) {
          sc.int('element id'); const type = sc.int('element type'), nt = sc.int('tag count'); let tag = 0;
          if (!(nt >= 0 && nt < 100)) throw new Error('Gmsh element tag count is not plausible');
          for (let k = 0; k < nt; k++) { const v = sc.int(); if (k === 0) tag = v; }
          const nn = nodesOf(type), ids = new Array(nn); for (let k = 0; k < nn; k++) ids[k] = sc.int('element node');
          add(type, ids, null, tag);
        }
      } else {
        const nb = guard(sc.int('element block count'), 'Gmsh element block'), total = guard(sc.int('element count'), 'Gmsh element'); sc.int(); sc.int();
        for (let k = 0, seen = 0; k < nb; k++) {
          const dim = sc.int(), etag = sc.int(), type = sc.int('element type'), n = sc.int('elements in block'), nn = nodesOf(type);
          if (!(n >= 0) || (seen += n) > total) throw new Error('Gmsh element blocks exceed the declared element count');
          for (let i = 0; i < n; i++) { sc.int('element tag'); const ids = new Array(nn); for (let j = 0; j < nn; j++) ids[j] = sc.int('element node'); add(type, ids, dim, entPhys.get(dim + ':' + etag) ?? 0); }
        }
      }
    } else if (sec === 'Periodic' || sec === 'NodeData' || sec === 'ElementData' || sec === 'ElementNodeData') (meta.otherSections ??= []).push(sec);
    for (let l2; (l2 = sc.line()) !== null;) if (l2.trimStart().startsWith('$End')) break;
  }
  if (!meta.version) throw new Error('no $MeshFormat section');
  for (const g of M.groups) g.kind = g.dim === 3 ? 'zone' : g.dim === maxDim ? 'surface' : 'boundary';
  return M.result({ meta });
}

// ---------- SU2 ----------
export function readSU2(b, ctx) {
  const sc = new Scanner(b), M = new Mesh(), meta = { dimension: null, markers: [] };
  let g = -1, zones = 1;
  const next = () => { for (let ln; (ln = sc.line()) !== null;) { const h = ln.indexOf('%'); if (h >= 0) ln = ln.slice(0, h); ln = ln.trim(); if (ln) return ln; } return null; };
  const elems = (n, grp) => {
    for (let i = 0; i < n; i++) {
      const ln = next(); if (ln === null) throw new Error('SU2 file ends inside an element list');
      const p = ln.split(/\s+/), t = VTK_TYPES[+p[0]];
      if (!t) throw new Error(`SU2 element type ${p[0]} is not recognised`);
      const nn = ET[t].n, ids = new Array(nn); for (let k = 0; k < nn; k++) ids[k] = parseInt(p[k + 1], 10);
      M.elem(t, ids, grp);
    }
  };
  for (let ln; (ln = next()) !== null;) {
    const m = /^([A-Z_0-9]+)\s*=\s*(.*)$/i.exec(ln); if (!m) continue;
    const key = m[1].toUpperCase(), val = m[2].trim();
    if (key === 'NZONE') zones = parseInt(val, 10);
    else if (key === 'IZONE' && parseInt(val, 10) > 1) { ctx.warn(`Multi-zone SU2 file (${zones} zones): only the first zone was read.`); break; }
    else if (key === 'NDIME') meta.dimension = parseInt(val, 10);
    else if (key === 'NELEM') elems(guard(parseInt(val, 10), 'SU2 element'), -1);
    else if (key === 'NPOIN') {
      const n = guard(parseInt(val, 10), 'SU2 point', LIMITS.verts), nd = meta.dimension === 2 ? 2 : 3;
      for (let i = 0; i < n; i++) { const l2 = next(); if (l2 === null) throw new Error('SU2 file ends inside the point list'); const p = l2.split(/\s+/); M.node(toNum(p[0]), toNum(p[1] ?? ''), nd === 3 ? toNum(p[2] ?? '') : 0); }
    } else if (key === 'MARKER_TAG') { g = M.group('boundary', val, val); meta.markers.push(val); }
    else if (key === 'MARKER_ELEMS') elems(guard(parseInt(val, 10), 'SU2 marker element'), g);
  }
  if (meta.dimension === null) throw new Error('no NDIME keyword');
  return M.result({ meta, kind: meta.dimension === 2 ? 'surface-mesh' : undefined });
}

// ---------- Nastran bulk data ----------
/** Nastran real: 1.5+3, 1.-3, 1.5D3 and ordinary forms; blank = 0. */
function nasNum(s) {
  if (!s) return 0;
  const v = +s; if (v === v) return v;
  const m = /^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$/.exec(s);
  return m ? +(m[1] + 'e' + m[2]) : +s.replace(/[dD]/, 'e');
}
const NAS_SOLID = { CTETRA: { 4: 'tet4', 10: 'tet10' }, CHEXA: { 8: 'hex8', 20: 'hex20' }, CPENTA: { 6: 'wedge6', 15: 'wedge15' }, CPYRAM: { 5: 'pyr5', 13: 'pyr13' } };
export function readNastran(b, ctx) {
  const sc = new Scanner(b), M = new Mesh(), census = {}, props = new Map(), mats = new Map(), cords = [];
  const meta = { cards: census, properties: [], materials: [], coordinateSystems: cords, gridsInLocalSystems: 0, includes: [] };
  const head = str(b, 0, Math.min(b.length, 1 << 20));
  let inBulk = !/^\s*BEGIN\s+BULK/im.test(head);   // a deck without BEGIN BULK is taken as bulk data throughout
  const I = (s) => parseInt(s, 10);
  const card = (f) => {
    const name = f[0]; census[name] = (census[name] || 0) + 1;
    const nodesFrom = (k, n) => { const a = []; for (let i = k; i < f.length && a.length < n; i++) if (f[i]) a.push(I(f[i])); return a; };
    switch (name) {
      case 'GRID': M.node(nasNum(f[3]), nasNum(f[4]), nasNum(f[5]), I(f[1])); if (f[2] && I(f[2]) !== 0) meta.gridsInLocalSystems++; break;
      case 'CTRIA3': case 'CTRIAR': M.later('tri3', [I(f[3]), I(f[4]), I(f[5])], M.group('property', I(f[2]))); break;
      case 'CQUAD4': case 'CQUADR': case 'CSHEAR': M.later('quad4', [I(f[3]), I(f[4]), I(f[5]), I(f[6])], M.group('property', I(f[2]))); break;
      case 'CTRIA6': M.later('tri6', nodesFrom(3, 6).length === 6 ? nodesFrom(3, 6) : [-1], M.group('property', I(f[2]))); break;
      case 'CQUAD8': M.later('quad8', nodesFrom(3, 8).length === 8 ? nodesFrom(3, 8) : [-1], M.group('property', I(f[2]))); break;
      case 'CTETRA': case 'CHEXA': case 'CPENTA': case 'CPYRAM': { const ids = nodesFrom(3, 20), t = NAS_SOLID[name][ids.length]; if (t) M.later(t, ids, M.group('property', I(f[2]))); else M.badRefs++; break; }
      case 'CBAR': case 'CBEAM': case 'CROD': case 'CTUBE': M.later('line2', [I(f[3]), I(f[4])], M.group('property', I(f[2]))); break;
      case 'CONROD': M.later('line2', [I(f[2]), I(f[3])], -1); break;
      case 'PSHELL': props.set(I(f[1]), { card: name, mid: I(f[2]), t: nasNum(f[3]) }); break;
      case 'PCOMP': case 'PCOMPG': props.set(I(f[1]), { card: name }); break;
      case 'PSOLID': case 'PBAR': case 'PBARL': case 'PBEAM': case 'PBEAML': case 'PROD': case 'PTUBE': case 'PSHEAR': props.set(I(f[1]), { card: name, mid: I(f[2]) }); break;
      case 'MAT1': mats.set(I(f[1]), { card: name, E: f[2] ? nasNum(f[2]) : null, G: f[3] ? nasNum(f[3]) : null, nu: f[4] ? nasNum(f[4]) : null, rho: f[5] ? nasNum(f[5]) : null }); break;
      case 'MAT2': case 'MAT8': case 'MAT9': case 'MAT10': mats.set(I(f[1]), { card: name }); break;
      case 'CORD2R': case 'CORD2C': case 'CORD2S': case 'CORD1R': case 'CORD1C': case 'CORD1S': cords.push({ card: name, id: I(f[1]) }); break;
      default:
    }
  };
  let cur = null;
  const flush = () => { if (cur) { card(cur); cur = null; } };
  for (let ln; (ln = sc.line()) !== null;) {
    if (!inBulk) { if (/^\s*BEGIN\s+BULK/i.test(ln)) inBulk = true; continue; }
    if (!ln.trim() || ln[0] === '$') continue;
    const dollar = ln.indexOf('$'); if (dollar > 0) ln = ln.slice(0, dollar);
    if (/^ENDDATA/i.test(ln)) break;
    if (/^INCLUDE/i.test(ln)) { meta.includes.push(ln.slice(7).trim().replace(/['"]/g, '')); continue; }
    let f, cont;
    if (ln.includes(',')) {                               // free field
      const p = ln.split(',').map((s) => s.trim()); cont = p[0] === '' || p[0][0] === '+' || (p[0][0] === '*' && cur);
      f = [p[0], ...p.slice(1, 9)];
    } else {                                              // fixed field: 8-character fields, or 16-character when the name carries '*'
      if (ln.includes('\t')) { let o = ''; for (const ch of ln) o += ch === '\t' ? ' '.repeat(8 - (o.length % 8)) : ch; ln = o; }
      const f0 = ln.slice(0, 8), large = f0.includes('*'); f = [f0.trim()];
      if (large) for (let k = 0; k < 4; k++) f.push(ln.slice(8 + 16 * k, 24 + 16 * k).trim());
      else for (let k = 1; k < 9; k++) f.push(ln.slice(8 * k, 8 * k + 8).trim());
      cont = f[0] === '' || f[0][0] === '+' || f[0][0] === '*';
    }
    if (cont && cur) { for (let k = 1; k < f.length; k++) cur.push(f[k]); continue; }
    if (cont) continue;
    flush();
    f[0] = f[0].replace('*', '').trim().toUpperCase();
    if (!/^[A-Z][A-Z0-9]*$/.test(f[0])) continue;
    cur = f;
  }
  flush();
  if (!M.nNodes && !M.pending) throw new Error('no GRID or element cards were found');
  const out = M.result({ meta });
  // property groups carry a readable name; material groups count the elements that reference them through a property
  const matCount = new Map();
  for (const g of M.groups) { const p = props.get(g.tag); g.name = p ? `${p.card} ${g.tag}${p.t ? ` (t = ${p.t})` : ''}${p.mid ? ` → MAT ${p.mid}` : ''}` : `PID ${g.tag}`; if (p?.mid) matCount.set(p.mid, (matCount.get(p.mid) || 0) + g.count); }
  for (const [id, m] of mats) M.groups.push({ id: M.groups.length, name: `${m.card} ${id}`, kind: 'material', count: matCount.get(id) || 0, tag: id });
  meta.properties = [...props].map(([id, p]) => ({ id, ...p })); meta.materials = [...mats].map(([id, m]) => ({ id, ...m }));
  if (meta.gridsInLocalSystems) ctx.warn(`${meta.gridsInLocalSystems} GRID point(s) are defined in a non-basic coordinate system (CP ≠ 0); the transformation is NOT applied, so their positions are wrong until the model is exported in the basic system.`);
  if (meta.includes.length) ctx.warn(`The deck INCLUDEs ${meta.includes.length} other file(s) (${meta.includes.slice(0, 3).join(', ')}…) which were not read.`);
  return out;
}

// ---------- Abaqus input deck ----------
function abqType(t) {
  let m;
  if ((m = /^C3D(\d+)/.exec(t))) return { 4: 'tet4', 10: 'tet10', 8: 'hex8', 20: 'hex20', 6: 'wedge6', 15: 'wedge15', 5: 'pyr5' }[m[1]];
  if ((m = /^SC(\d)/.exec(t))) return { 6: 'wedge6', 8: 'hex8' }[m[1]];
  if ((m = /^(?:STRI|S|M3D|R3D|DS|SFM3D|CPS|CPE|CPEG|CAX|DC2D|AC2D)(\d)/.exec(t))) return { 3: 'tri3', 4: 'quad4', 6: 'tri6', 8: 'quad8', 9: 'quad9' }[m[1]];
  if ((m = /^(?:B|PIPE|FRAME)\d([123])/.exec(t))) return m[1] === '2' ? 'line3' : 'line2';
  if ((m = /^T\dD([23])/.exec(t))) return m[1] === '3' ? 'line3' : 'line2';
  return null;
}
export function readAbaqus(b, ctx) {
  const sc = new Scanner(b), M = new Mesh();
  const meta = { heading: null, parts: [], instances: [], nsets: [], elsets: [], materials: [], sections: [], elementTypes: {}, unsupportedElementTypes: {} };
  const global = { name: '', ids: new Map(), n0: 0 }; let part = global;
  const parts = new Map(); let mode = null, etype = null, eg = -1, pend = [], setRec = null, gen = false, inst = null, headLines = 0;
  const params = (ln) => { const o = {}; for (const p of ln.split(',').slice(1)) { const e = p.indexOf('='); if (e < 0) o[p.trim().toUpperCase()] = true; else o[p.slice(0, e).trim().toUpperCase()] = p.slice(e + 1).trim().replace(/^"|"$/g, ''); } return o; };
  for (let ln; (ln = sc.line()) !== null;) {
    if (ln.startsWith('**') || !ln.trim()) continue;
    if (ln[0] === '*') {
      const kw = ln.split(',')[0].slice(1).trim().toUpperCase().replace(/\s+/g, ' '), P = params(ln);
      mode = null; pend = [];
      if (kw === 'HEADING') { mode = 'heading'; headLines = 0; }
      else if (kw === 'PART') { part = { name: P.NAME || `part ${parts.size + 1}`, ids: new Map(), n0: M.nNodes, n1: M.nNodes }; parts.set(part.name.toUpperCase(), part); meta.parts.push(part.name); }
      else if (kw === 'END PART') { part.n1 = M.nNodes; part = global; }
      else if (kw === 'INSTANCE') { inst = { name: P.NAME, part: P.PART, data: [] }; meta.instances.push(inst); mode = 'instance'; const src = parts.get(String(P.PART).toUpperCase()); part = src && !src.used ? src : { name: P.NAME || '', ids: new Map(), n0: M.nNodes }; }
      else if (kw === 'END INSTANCE') { part = global; inst = null; }
      else if (kw === 'NODE') mode = 'node';
      else if (kw === 'ELEMENT') {
        const tn = String(P.TYPE || '').toUpperCase(); etype = abqType(tn);
        if (etype) meta.elementTypes[tn] = meta.elementTypes[tn] || 0; else meta.unsupportedElementTypes[tn] = 0;
        etype = etype ? { t: etype, tn } : { t: null, tn };
        const gname = P.ELSET ? (part.name ? `${part.name}.${P.ELSET}` : P.ELSET) : part.name || null;
        eg = gname ? M.group('component', gname, gname) : -1; mode = 'element';
      } else if (kw === 'NSET' || kw === 'ELSET') { setRec = { name: P[kw] || '(unnamed)', count: 0 }; (kw === 'NSET' ? meta.nsets : meta.elsets).push(setRec); gen = !!P.GENERATE; mode = 'set'; }
      else if (kw === 'MATERIAL') meta.materials.push(P.NAME || '(unnamed)');
      else if (/SECTION$/.test(kw)) meta.sections.push({ type: kw, elset: P.ELSET ?? null, material: P.MATERIAL ?? null });
      continue;
    }
    if (mode === 'heading') { if (headLines++ === 0) meta.heading = ln.trim(); }
    else if (mode === 'node') { const p = ln.split(','); const id = parseInt(p[0], 10); if (id === id) part.ids.set(id, M.node(toNum((p[1] ?? '').trim()), toNum((p[2] ?? '').trim()) || 0, toNum((p[3] ?? '').trim()) || 0)); }
    else if (mode === 'element') {
      if (!etype.t) { meta.unsupportedElementTypes[etype.tn]++; continue; }
      for (const s of ln.split(',')) { const t = s.trim(); if (t) pend.push(parseInt(t, 10)); }
      const need = ET[etype.t].n + 1;
      while (pend.length >= need) { const ids = pend.splice(0, need), n = new Array(need - 1); for (let k = 1; k < need; k++) n[k - 1] = part.ids.get(ids[k]) ?? -1; M.elem(etype.t, n, eg); meta.elementTypes[etype.tn]++; }
    } else if (mode === 'set') { const v = ln.split(',').map((s) => s.trim()).filter(Boolean); if (gen) { const a = +v[0], e = +v[1], s = +v[2] || 1; if (e >= a && s > 0) setRec.count += Math.floor((e - a) / s) + 1; } else setRec.count += v.length; }
    else if (mode === 'instance' && inst) inst.data.push(nums(ln));
  }
  // apply the first instance's placement to its part (translation, then rotation about the axis a→b)
  const X = M.xyz;
  for (const it of meta.instances) {
    const p = parts.get(String(it.part).toUpperCase()); if (!p) continue;
    if (p.used) { ctx.warn(`Part "${p.name}" is instanced more than once; only its first instance ("${p.used}") is materialised.`); continue; }
    p.used = it.name || true;
    const t = it.data[0] && it.data[0].length >= 3 ? it.data[0] : null, r = it.data[1] && it.data[1].length >= 7 ? it.data[1] : null;
    if (!t && !r) continue;
    let ax = null, c = 1, s = 0;
    if (r) { const d = [r[3] - r[0], r[4] - r[1], r[5] - r[2]], L = Math.hypot(...d); if (L > 0) { ax = d.map((v) => v / L); c = Math.cos((r[6] * Math.PI) / 180); s = Math.sin((r[6] * Math.PI) / 180); } }
    for (let i = p.n0; i < p.n1; i++) {
      let x = X[3 * i] + (t ? t[0] : 0), y = X[3 * i + 1] + (t ? t[1] : 0), z = X[3 * i + 2] + (t ? t[2] : 0);
      if (ax) {                                           // Rodrigues rotation about the axis through r[0..2]
        const px = x - r[0], py = y - r[1], pz = z - r[2], dt = ax[0] * px + ax[1] * py + ax[2] * pz;
        const cx = ax[1] * pz - ax[2] * py, cy = ax[2] * px - ax[0] * pz, cz = ax[0] * py - ax[1] * px;
        x = r[0] + px * c + cx * s + ax[0] * dt * (1 - c); y = r[1] + py * c + cy * s + ax[1] * dt * (1 - c); z = r[2] + pz * c + cz * s + ax[2] * dt * (1 - c);
      }
      X[3 * i] = x; X[3 * i + 1] = y; X[3 * i + 2] = z;
    }
  }
  meta.instances = meta.instances.map((it) => ({ name: it.name, part: it.part, translation: it.data[0] ?? null, rotation: it.data[1] ?? null }));
  const un = Object.entries(meta.unsupportedElementTypes).filter(([, n]) => n > 0);
  if (un.length) ctx.warn(`Element types not mapped and skipped: ${un.map(([t, n]) => `${t} (${n})`).join(', ')}.`);
  if (!M.nNodes) throw new Error('no *NODE data were found');
  const out = M.result({ meta });
  for (const name of meta.materials) M.groups.push({ id: M.groups.length, name, kind: 'material', count: 0 });
  return out;
}

// ---------- I-deas universal file ----------
const UNV_UNITS = { 1: ['m', 'SI: metre (newton)'], 2: ['ft', 'BG: foot (pound-force)'], 3: ['m', 'MG: metre (kilogram-force)'], 4: ['ft', 'BA: foot (poundal)'], 5: ['mm', 'MM: millimetre (milli-newton)'], 6: ['cm', 'CM: centimetre (centi-newton)'], 7: ['in', 'IN: inch (pound-force)'], 8: ['mm', 'GM: millimetre (kilogram-force)'], 9: [null, 'US: user defined'], 10: ['mm', 'MN: millimetre (newton)'] };
const UNV_FE = { 11: 'line2', 21: 'line2', 22: 'line2', 23: 'line2', 24: 'line3', 41: 'tri3', 51: 'tri3', 61: 'tri3', 74: 'tri3', 81: 'tri3', 91: 'tri3', 42: 'tri6', 52: 'tri6', 62: 'tri6', 72: 'tri6', 82: 'tri6', 92: 'tri6', 44: 'quad4', 54: 'quad4', 64: 'quad4', 71: 'quad4', 84: 'quad4', 94: 'quad4', 45: 'quad8', 55: 'quad8', 65: 'quad8', 75: 'quad8', 85: 'quad8', 95: 'quad8', 111: 'tet4', 118: 'tet10', 112: 'wedge6', 113: 'wedge15', 115: 'hex8', 116: 'hex20' };
// UNV lists parabolic elements corner, mid, corner, mid …; reorder so the corner nodes come first.
const UNV_ORDER = { tri6: [0, 2, 4, 1, 3, 5], quad8: [0, 2, 4, 6, 1, 3, 5, 7], line3: [0, 2, 1], tet10: [0, 2, 4, 9, 1, 3, 5, 6, 7, 8], wedge15: [0, 2, 4, 9, 11, 13, 1, 3, 5, 6, 7, 8, 10, 12, 14], hex20: [0, 2, 4, 6, 12, 14, 16, 18, 1, 3, 5, 7, 8, 9, 10, 11, 13, 15, 17, 19] };
export function readUNV(b, ctx) {
  const sc = new Scanner(b), M = new Mesh(), meta = { datasets: {}, unitSystem: null, unsupportedElements: {} };
  let units = { length: null, source: null };
  const eid = new Map(), groupsRaw = [];
  const isDelim = (ln) => ln.trim() === '-1';
  for (let ln; (ln = sc.line()) !== null;) {
    if (!isDelim(ln)) continue;
    const ds = sc.line(); if (ds === null) break;
    const id = parseInt(ds, 10); if (id !== id) continue;
    meta.datasets[id] = (meta.datasets[id] || 0) + 1;
    let ended = false;   // every dataset is consumed up to its closing -1 so the delimiter pairing never slips
    const body = () => { if (ended) return null; const l2 = sc.line(); if (l2 === null || isDelim(l2)) { ended = true; return null; } return l2; };
    if (id === 164) {
      const l2 = body(); if (l2 !== null) { const code = parseInt(l2, 10), u = UNV_UNITS[code]; meta.unitSystem = u ? u[1] : l2.trim(); if (u && u[0]) units = { length: u[0], source: 'file' }; if (code !== 1 && u && u[0]) ctx.warn(`UNV dataset 164 declares "${u[1]}". Some writers store coordinates in SI regardless of this record — check the extents against a known dimension.`); }
      for (let l3; (l3 = body()) !== null;);
    } else if (id === 2411) {
      for (let l1; (l1 = body()) !== null;) { const l2 = body(); if (l2 === null) break; const c = nums(l2); M.node(c[0], c[1] ?? 0, c[2] ?? 0, parseInt(l1, 10)); }
    } else if (id === 2412) {
      for (let l1; (l1 = body()) !== null;) {
        const r = nums(l1), fe = r[1], nn = r[5]; if (!(nn > 0 && nn <= 64)) break;
        if (fe >= 11 && fe <= 32) { if (body() === null) break; }       // beam elements carry an extra orientation record
        const ids = []; let ok = true;
        while (ids.length < nn) { const l2 = body(); if (l2 === null) { ok = false; break; } ids.push(...nums(l2)); }
        if (!ok) break;
        const t = UNV_FE[fe];
        if (!t || ET[t].n !== nn) { meta.unsupportedElements[fe] = (meta.unsupportedElements[fe] || 0) + 1; continue; }
        const ord = UNV_ORDER[t]; eid.set(r[0], [t, ord ? ord.map((k) => ids[k]) : ids.slice(0, nn)]);
      }
    } else if (id === 2467 || id === 2477 || id === 2452 || id === 2435) {
      for (let l1; (l1 = body()) !== null;) {
        const r = nums(l1), n = r[7]; const name = (body() ?? '').trim(); if (!(n >= 0)) break;
        const g = { name: name || `group ${r[0]}`, elems: [], nodes: 0 };
        for (let got = 0; got < n;) { const l3 = body(); if (l3 === null) break; const v = nums(l3); for (let k = 0; k + 1 < v.length && got < n; k += id === 2467 || id === 2477 ? 4 : 2, got++) { if (v[k] === 8) g.elems.push(v[k + 1]); else if (v[k] === 7) g.nodes++; } }
        groupsRaw.push(g);
      }
    }
    while (body() !== null);
  }
  const eg = new Map();
  for (const g of groupsRaw) { if (!g.elems.length) { M.groups.push({ id: M.groups.length, name: g.name, kind: 'component', count: 0, nodes: g.nodes }); continue; } const gi = M.group('component', g.name, g.name); for (const e of g.elems) if (!eg.has(e)) eg.set(e, gi); }
  for (const [id, [t, ids]] of eid) M.later(t, ids, eg.get(id) ?? -1);
  if (!M.nNodes) throw new Error('no dataset 2411 (nodes) was found');
  const un = Object.entries(meta.unsupportedElements);
  if (un.length) ctx.warn(`Elements with unmapped FE descriptors were skipped: ${un.map(([t, n]) => `${t} (${n})`).join(', ')}.`);
  return M.result({ meta, units });
}

// ---------- ANSYS CDB ----------
const ANSYS_SHELL = new Set([28, 41, 43, 63, 93, 131, 132, 157, 163, 181, 281, 42, 82, 182, 183, 55, 77, 152, 153, 154, 208, 209]);
const ANSYS_LINE = new Set([1, 3, 4, 8, 10, 23, 24, 44, 54, 180, 188, 189, 288, 289, 33, 14, 11]);
const ANSYS_UNITS = { SI: 'm', MKS: 'm', CGS: 'cm', MPA: 'mm', BFT: 'ft', BIN: 'in' };
export function readCDB(b, ctx) {
  const sc = new Scanner(b), M = new Mesh(), et = new Map(), meta = { elementTypes: {}, components: [], units: null, unmapped: 0 };
  let units = { length: null, source: null };
  const fixedInts = (ln, w) => { const o = []; for (let p = 0; p + 1 <= ln.length; p += w) { const s = ln.slice(p, p + w).trim(); if (!s) break; o.push(parseInt(s, 10)); } return o; };
  for (let ln; (ln = sc.line()) !== null;) {
    const up = ln.trimStart().toUpperCase();
    if (up.startsWith('ET,')) { const p = up.split(','); const id = parseInt(p[1], 10), num = parseInt(p[2], 10) || parseInt(String(p[2]).replace(/\D+/g, ''), 10); et.set(id, num); meta.elementTypes[id] = num; }
    else if (up.startsWith('/UNITS')) { const u = up.split(',')[1]?.trim(); meta.units = u ?? null; if (ANSYS_UNITS[u]) units = { length: ANSYS_UNITS[u], source: 'file' }; }
    else if (up.startsWith('CMBLOCK')) { const p = ln.split(','); meta.components.push({ name: (p[1] || '').trim(), kind: (p[2] || '').trim(), entries: parseInt(p[3], 10) || 0 }); }
    else if (up.startsWith('NBLOCK')) {
      const m = /\((\d+)i(\d+)\s*,\s*(\d+)e(\d+)/i.exec(sc.line() ?? ''); if (!m) throw new Error('NBLOCK format line not understood');
      const ni = +m[1], wi = +m[2], wf = +m[4], c0 = ni * wi;
      for (let l2; (l2 = sc.line()) !== null;) {
        const id = parseInt(l2.slice(0, wi), 10); if (!(id > 0) || /^\s*N\s*,/i.test(l2)) break;
        M.node(toNum(l2.slice(c0, c0 + wf).trim()) || 0, toNum(l2.slice(c0 + wf, c0 + 2 * wf).trim()) || 0, toNum(l2.slice(c0 + 2 * wf, c0 + 3 * wf).trim()) || 0, id);
      }
    } else if (up.startsWith('EBLOCK')) {
      if (!/SOLID/.test(up)) { ctx.warn('An EBLOCK without the SOLID key (non-solid layout) was skipped.'); continue; }
      const m = /\((\d+)i(\d+)/i.exec(sc.line() ?? ''); if (!m) throw new Error('EBLOCK format line not understood');
      const w = +m[2];
      for (let l2; (l2 = sc.line()) !== null;) {
        const f = fixedInts(l2, w); if (!f.length || f[0] === -1) break;
        if (f.length < 11) continue;
        const nn = f[8]; let n = f.slice(11);
        while (n.length < nn) { const l3 = sc.line(); if (l3 === null) break; n = n.concat(fixedInts(l3, w)); }
        if (n.length < nn || !(nn > 0)) continue;
        n.length = nn;
        const num = et.get(f[1]), g = M.group('material', f[0], `MAT ${f[0]}`);
        let t = null, ids = n;
        if (num !== undefined && ANSYS_LINE.has(num)) { t = 'line2'; ids = n.slice(0, 2); }
        else if (nn === 8 && num !== undefined && ANSYS_SHELL.has(num)) t = 'quad8';
        else if (nn === 8) [t, ids] = degen8(n);
        else if (nn === 20) t = 'hex20';
        else if (nn === 10) t = 'tet10';
        else if (nn === 4 && num === 285) t = 'tet4';
        else if (nn === 4) { if (n[2] === n[3]) { t = 'tri3'; ids = n.slice(0, 3); } else t = 'quad4'; }
        else if (nn === 3 && (num === undefined || ANSYS_SHELL.has(num))) t = 'tri3';
        else if (nn === 6) t = 'tri6';
        else if (nn === 2 || nn === 3) { t = 'line2'; ids = n.slice(0, 2); }
        if (t) M.later(t, ids, g); else meta.unmapped++;
      }
    }
  }
  if (!M.nNodes) throw new Error('no NBLOCK node data were found');
  if (meta.unmapped) ctx.warn(`${meta.unmapped} element(s) with an unmapped node count/type were skipped.`);
  return M.result({ meta, units });
}

// ---------- Fluent mesh (ASCII) ----------
const FLUENT_BC = { 2: 'interior', 3: 'wall', 4: 'pressure-inlet', 5: 'pressure-outlet', 7: 'symmetry', 8: 'periodic-shadow', 9: 'pressure-far-field', 10: 'velocity-inlet', 12: 'periodic', 14: 'fan/porous-jump', 20: 'mass-flow-inlet', 24: 'interface', 31: 'parent', 36: 'outflow', 37: 'axis' };
export function readFluent(b, ctx) {
  const n = b.length, M = new Mesh(), meta = { dimension: null, nodes: 0, faces: 0, cells: 0, zones: [], header: null, volumeCellsReconstructed: false };
  const names = new Map(), faceZones = [], cellZones = []; let p = 0, pos = null, binary = false;
  // end of the list that opens at byte `open` ('('), ignoring parentheses inside strings
  const close = (open) => { let d = 0, q = false; for (let i = open; i < n; i++) { const c = b[i]; if (c === 34) q = !q; else if (!q) { if (c === 40) d++; else if (c === 41 && --d === 0) return i; } } return -1; };
  const hex = (s) => parseInt(s, 16);
  while (p < n) {
    while (p < n && b[p] !== 40) p++;
    if (p >= n) break;
    let q = p + 1; while (q < n && b[q] <= 32) q++;
    let e = q; while (e < n && b[e] > 32 && b[e] !== 40 && b[e] !== 41) e++;
    const idx = parseInt(str(b, q, e), 10);
    if (idx !== idx) { const c = close(p); p = c < 0 ? n : c + 1; continue; }
    if (idx >= 2000) { binary = true; break; }                              // binary sections cannot be skipped safely
    if (idx === 10 || idx === 12 || idx === 13) {
      let h0 = e; while (h0 < n && b[h0] !== 40) h0++;
      const h1 = close(h0); if (h1 < 0) break;
      const h = str(b, h0 + 1, h1).trim().split(/\s+/).map(hex), [zone, first, last, type, kind] = h;
      let body0 = h1 + 1; while (body0 < n && b[body0] <= 32) body0++;
      let body1 = -1;
      if (b[body0] === 40) { body1 = body0; while (body1 < n && b[body1] !== 41) body1++; }
      const count = last - first + 1;
      if (zone === 0) { if (idx === 10) { meta.nodes = last; pos = new Float64Array(3 * guard(last, 'Fluent node', LIMITS.verts)); if (kind) meta.dimension ??= kind; } else if (idx === 13) meta.faces = last; else meta.cells = last; }
      else if (idx === 10 && body1 > 0) {
        const nd = kind || meta.dimension || 3, s = new Scanner(b.subarray(body0 + 1, body1));
        if (!pos || 3 * last > pos.length) { const np = new Float64Array(3 * guard(last, 'Fluent node', LIMITS.verts)); if (pos) np.set(pos); pos = np; meta.nodes = Math.max(meta.nodes, last); }
        for (let i = first - 1; i < last; i++) { pos[3 * i] = s.float('Fluent coordinate'); pos[3 * i + 1] = s.float('Fluent coordinate'); pos[3 * i + 2] = nd === 3 ? s.float('Fluent coordinate') : 0; }
      } else if (idx === 13 && body1 > 0) {
        const z = { id: zone, bc: type, faces: count, faceType: kind }; faceZones.push(z);
        if (type !== 2) {                                                     // boundary zone: keep its faces as a surface
          const s = new Scanner(b.subarray(body0 + 1, body1)); z.list = [];
          for (let i = 0; i < guard(count, 'Fluent face'); i++) { const k = kind === 0 || kind === 5 ? hex(s.token() ?? '') : kind; if (!(k >= 2 && k <= 1000)) throw new Error('Fluent face record is not plausible'); const f = new Array(k); for (let j = 0; j < k; j++) f[j] = hex(s.token() ?? '') - 1; s.token(); s.token(); z.list.push(f); }
        }
      } else if (idx === 12) cellZones.push({ id: zone, cells: count, active: type, cellType: kind });
      const endAll = close(p); p = endAll < 0 ? n : endAll + 1;
      continue;
    }
    const c = close(p); if (c < 0) break;
    const inner = str(b, e, Math.min(c, e + 400)).trim();
    if (idx === 2) meta.dimension = parseInt(inner, 10);
    else if (idx === 1 || idx === 0) { if (idx === 1) meta.header = inner.replace(/"/g, ''); }
    else if (idx === 39 || idx === 45) { const m = /\(\s*(\d+)\s+(\S+)\s+([^\s()]+)/.exec(inner); if (m) names.set(+m[1], { type: m[2], name: m[3] }); }
    p = c + 1;
  }
  if (binary) ctx.warn('The file contains binary sections, which are not read; write the mesh as ASCII.');
  if (!pos) throw new Error('no ASCII node section (10 …) was found');
  const nv = pos.length / 3; M.xyz = Array.from(pos);
  for (const z of faceZones) {
    const nm = names.get(z.id), bcName = FLUENT_BC[z.bc] ?? `type ${z.bc}`;
    meta.zones.push({ id: z.id, kind: 'faces', type: nm?.type ?? bcName, name: nm?.name ?? null, count: z.faces });
    if (!z.list) continue;
    const g = M.group('boundary', z.id, nm?.name ?? `zone ${z.id} (${bcName})`);
    for (const f of z.list) {
      if (f.some((v) => !(v >= 0 && v < nv))) { M.badRefs++; continue; }
      if (f.length === 2) M.elem('line2', f, g); else if (f.length === 3) M.elem('tri3', f, g); else if (f.length === 4) M.elem('quad4', f, g); else for (let j = 2; j < f.length; j++) M.elem('tri3', [f[0], f[j - 1], f[j]], g);
    }
  }
  for (const z of cellZones) { const nm = names.get(z.id); meta.zones.push({ id: z.id, kind: 'cells', type: nm?.type ?? 'fluid', name: nm?.name ?? null, count: z.cells }); M.groups.push({ id: M.groups.length, name: nm?.name ?? `cell zone ${z.id}`, kind: 'zone', count: z.cells, tag: z.id }); }
  ctx.warn(`Fluent mesh: ${meta.nodes} nodes, ${meta.faces} faces, ${meta.cells} cells declared. Boundary face zones are shown as a surface; cell connectivity is NOT reconstructed, so this model carries no volume elements. Face orientation is as stored.`);
  return M.result({ meta, kind: M.blocks.size ? 'surface-mesh' : 'pointcloud' });
}

// ---------- Tecplot ASCII ----------
const TEC_ET = { TRIANGLE: 'tri3', QUADRILATERAL: 'quad4', TETRAHEDRON: 'tet4', BRICK: 'hex8', LINESEG: 'line2' };
export function readTecplot(b, ctx) {
  const sc = new Scanner(b, { sep: ',' }), M = new Mesh(), meta = { title: null, variables: [], zones: [] };
  let pre = '', zone = null, vars = null, ignored = 0;
  const coordCols = () => {
    if (!vars) { const m = /VARIABLES\s*=\s*([\s\S]*)$/i.exec(pre); let v = []; if (m) { const q = m[1].match(/"[^"]*"/g); v = q ? q.map((s) => s.slice(1, -1)) : m[1].split(/[\s,]+/).filter(Boolean); } vars = v; meta.variables = v; const t = /TITLE\s*=\s*"([^"]*)"/i.exec(pre); meta.title = t ? t[1] : null; }
    const find = (c) => vars.findIndex((v) => new RegExp(`^${c}($|[^a-z])|^coordinate${c}$`, 'i').test(v.trim()));
    let ix = find('x'), iy = find('y'), iz = find('z');
    if (ix < 0 || iy < 0) { ix = 0; iy = 1; iz = vars.length === 0 || vars.length > 2 ? 2 : -1; }
    return { nv: vars.length || 3, ix, iy, iz };
  };
  const data = (hdr) => {
    const kv = {}; for (const m of hdr.replace(/^\s*ZONE\b/i, '').matchAll(/([A-Za-z]\w*)\s*=\s*("[^"]*"|\([^)]*\)|[^\s,]+)/g)) kv[m[1].toUpperCase()] = m[2].replace(/^"|"$/g, '');
    const { nv, ix, iy, iz } = coordCols(), fmt = (kv.F || '').toUpperCase(), zt = (kv.ZONETYPE || '').toUpperCase();
    const fe = zt.startsWith('FE') || fmt.startsWith('FE'), block = (kv.DATAPACKING || '').toUpperCase() === 'BLOCK' || fmt.endsWith('BLOCK');
    const etName = zt.startsWith('FE') ? zt.slice(2) : (kv.ET || 'TRIANGLE').toUpperCase(), et = TEC_ET[etName];
    if (fe && !et) throw new Error(`Tecplot zone type ${zt || kv.ET} is not read`);
    if (kv.VARLOCATION && /CELLCENTERED/i.test(kv.VARLOCATION)) throw new Error('Tecplot zone with cell-centred variables is not read');
    if (kv.VARSHARELIST || kv.CONNECTIVITYSHAREZONE) throw new Error('Tecplot zones that share variables/connectivity are not read');
    const I = parseInt(kv.I ?? 1, 10), J = parseInt(kv.J ?? 1, 10), K = parseInt(kv.K ?? 1, 10);
    const n = guard(fe ? parseInt(kv.N ?? kv.NODES, 10) : I * J * K, 'Tecplot node', LIMITS.verts), ne = fe ? guard(parseInt(kv.E ?? kv.ELEMENTS, 10), 'Tecplot element') : 0;
    const base = M.nNodes, name = kv.T ?? `zone ${meta.zones.length + 1}`, g = M.group('zone', M.groups.length, name);
    if (block) { const cols = []; for (let v = 0; v < nv; v++) { const keep = v === ix || v === iy || v === iz, a = keep ? new Float64Array(n) : null; for (let i = 0; i < n; i++) { const x = sc.float('Tecplot value'); if (a) a[i] = x; } cols.push(a); } for (let i = 0; i < n; i++) M.node(cols[ix][i], cols[iy][i], iz >= 0 ? cols[iz][i] : 0); }
    else for (let i = 0; i < n; i++) { let x = 0, y = 0, z = 0; for (let v = 0; v < nv; v++) { const q = sc.float('Tecplot value'); if (v === ix) x = q; else if (v === iy) y = q; else if (v === iz) z = q; } M.node(x, y, z); }
    if (fe) {
      const nn = ET[et].n;
      for (let e = 0; e < ne; e++) {
        const ids = new Array(nn); for (let k = 0; k < nn; k++) { const v = sc.int('Tecplot connectivity'); ids[k] = v >= 1 && v <= n ? base + v - 1 : -1; }
        if (et === 'hex8') { const [t, c] = degen8(ids); M.elem(t, c, g); } else if (et === 'quad4' && ids[2] === ids[3]) M.elem('tri3', ids.slice(0, 3), g); else M.elem(et, ids, g);
      }
    } else structuredBlock(M, I, J, K, base, g);
    meta.zones.push({ name, type: fe ? 'FE' + etName : 'ORDERED', nodes: n, elements: fe ? ne : null, ijk: fe ? null : [I, J, K], packing: block ? 'BLOCK' : 'POINT' });
  };
  for (;;) {
    const save = sc.p, ln = sc.line(); if (ln === null) break;
    const t = ln.trim(); if (!t || t[0] === '#') continue;
    if (/^[-+.\d]/.test(t)) { if (!zone) { ignored++; continue; } sc.p = save; const h = zone; zone = null; data(h); continue; }
    if (/^ZONE\b/i.test(t)) zone = t; else if (zone !== null) zone += ' ' + t; else pre += ' ' + t;
  }
  if (!M.nNodes) throw new Error('no zone data were found');
  if (ignored) ctx.warn(`${ignored} numeric line(s) outside any zone (text/geometry records) were ignored.`);
  if (meta.variables.length > 3) ctx.warn(`Field variables (${meta.variables.slice(0, 8).join(', ')}${meta.variables.length > 8 ? '…' : ''}) are listed in the metadata but not loaded.`);
  return M.result({ meta });
}

// ---------- Plot3D ASCII (formatted, whole, optionally multi-block) ----------
export function readPlot3D(b, ctx) {
  const sc = new Scanner(b, { sep: ',' }), first = nums(sc.line() ?? '');
  let dims = [], two = false;
  if (first.length === 1) { const nb = guard(first[0], 'Plot3D block', 100000); for (let k = 0; k < nb; k++) dims.push([sc.int('block dimension'), sc.int('block dimension'), sc.int('block dimension')]); }
  else if (first.length === 3) dims = [first];
  else if (first.length === 2) { dims = [[first[0], first[1], 1]]; two = true; }
  else throw new Error('the first line is not a Plot3D block count or block dimension record');
  let total = 0; for (const d of dims) { if (!d.every((v) => Number.isInteger(v) && v >= 1)) throw new Error('Plot3D block dimensions are not positive integers'); total += d[0] * d[1] * d[2]; }
  guard(total, 'Plot3D point', LIMITS.verts);
  const M = new Mesh();
  dims.forEach((d, k) => {
    const n = d[0] * d[1] * d[2], base = M.nNodes, c = [new Float64Array(n), new Float64Array(n), new Float64Array(n)];
    for (let a = 0; a < (two ? 2 : 3); a++) for (let i = 0; i < n; i++) c[a][i] = sc.float('Plot3D coordinate');
    for (let i = 0; i < n; i++) M.node(c[0][i], c[1][i], c[2][i]);
    structuredBlock(M, d[0], d[1], d[2], base, M.group('zone', k, `block ${k + 1} (${d.join('×')})`));
  });
  const iblank = sc.peek() !== null;
  if (iblank) ctx.warn('Extra values follow the coordinates (IBLANK array or further data); they were ignored.');
  return M.result({ meta: { blocks: dims.map((d) => ({ ni: d[0], nj: d[1], nk: d[2] })), planar: two, iblankPresent: iblank, structured: true } });
}

// ---------- OpenFOAM polyMesh (points + faces [+ boundary, owner]) ----------
export function readOpenFOAM(b, ctx) {
  const files = {};
  for (const f of [{ name: ctx.name, bytes: b }, ...(ctx.opts.companions || [])]) {
    if (!(f?.bytes instanceof Uint8Array)) continue;
    const head = str(f.bytes, 0, Math.min(f.bytes.length, 2000)), obj = /\bobject\s+(\w+)\s*;/.exec(head)?.[1] ?? String(f.name).split(/[\\/]/).pop();
    const fmt = /\bformat\s+(\w+)\s*;/.exec(head)?.[1] ?? 'ascii', at = head.indexOf('FoamFile'), end = at < 0 ? 0 : head.indexOf('}', at) + 1;
    if (['points', 'faces', 'boundary', 'owner', 'neighbour'].includes(obj)) files[obj] = { bytes: f.bytes, fmt, start: Math.max(0, end), cls: /\bclass\s+(\w+)\s*;/.exec(head)?.[1] ?? null, note: /\bnote\s+"([^"]*)"/.exec(head)?.[1] ?? null };
  }
  const meta = { filesSupplied: Object.keys(files), patches: [], cells: null, cellsReconstructed: false };
  for (const k of Object.keys(files)) if (files[k].fmt !== 'ascii') throw Object.assign(new Error(`OpenFOAM file "${k}" is written in ${files[k].fmt} format; only ascii is read (foamFormatConvert)`), { meta });
  if (!files.points || !files.faces) throw Object.assign(new Error(`OpenFOAM polyMesh needs both "points" and "faces" (supplied: ${meta.filesSupplied.join(', ') || 'none'}); pass the missing files in opts.companions`), { meta });
  const M = new Mesh(), open = (f) => new Scanner(f.bytes, { pos: f.start, comment: '/', sep: '()' });
  let sc = open(files.points); const np = guard(sc.int('point count'), 'OpenFOAM point', LIMITS.verts);
  for (let i = 0; i < np; i++) M.node(sc.float('point'), sc.float('point'), sc.float('point'));
  sc = open(files.faces); const nf = guard(sc.int('face count'), 'OpenFOAM face'), faces = new Array(nf);
  for (let i = 0; i < nf; i++) { const k = sc.int('face size'); if (!(k >= 3 && k <= 10000)) throw new Error('OpenFOAM face size is not plausible'); const f = new Array(k); for (let j = 0; j < k; j++) f[j] = sc.int('face vertex'); faces[i] = f; }
  const emit = (s, e, g) => { for (let i = s; i < e && i < nf; i++) { const f = faces[i]; if (f.length === 3) M.elem('tri3', f, g); else if (f.length === 4) M.elem('quad4', f, g); else for (let j = 2; j < f.length; j++) M.elem('tri3', [f[0], f[j - 1], f[j]], g); } };
  if (files.boundary) {
    const src = text(files.boundary.bytes.subarray(files.boundary.start)).replace(/\/\/[^\n]*/g, '');
    for (const m of src.matchAll(/([A-Za-z_][\w.-]*)\s*\{([^{}]*)\}/g)) {
      const nF = +(/\bnFaces\s+(\d+)/.exec(m[2])?.[1] ?? NaN), sF = +(/\bstartFace\s+(\d+)/.exec(m[2])?.[1] ?? NaN), type = /\btype\s+(\w+)/.exec(m[2])?.[1] ?? 'patch';
      if (!(nF >= 0 && sF >= 0)) continue;
      meta.patches.push({ name: m[1], type, nFaces: nF, startFace: sF });
      emit(sF, sF + nF, M.group('boundary', m[1], `${m[1]} (${type})`));
    }
  } else { emit(0, nf, -1); ctx.warn('No "boundary" file was supplied: ALL faces (internal ones included) are shown and no patch names are available.'); }
  if (files.owner) { const m = /nCells:\s*(\d+)/.exec(files.owner.note ?? ''); meta.cells = m ? +m[1] : null; }
  ctx.warn('OpenFOAM polyMesh: boundary patches are shown as a surface; polyhedral cells are NOT reconstructed, so this model carries no volume elements.');
  return M.result({ meta, kind: 'surface-mesh' });
}
