// High-fidelity bridge tests: `node tests/bridge.mjs`
// Builds every deck for several presets and checks structural validity, round-trips the ZIP code,
// parses representative solver outputs (trimmed from real SU2 8.5.0, OpenFOAM v2412 and CalculiX 2.23
// runs of generated cases) and exercises the GitHub client against a mock.
//
// Optional end-to-end checks against the real solvers (never downloaded by this script):
//   BRIDGE_SU2=/path/to/SU2_CFD  BRIDGE_CCX=/path/to/ccx  BRIDGE_GMSH=/path/to/gmsh
//   BRIDGE_FOAM_BASHRC=/path/to/openfoam/etc/bashrc  [BRIDGE_SCRATCH=/dir]  node tests/bridge.mjs
// Exits non-zero when any check fails.

import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { PRESETS } from '../js/core/case.js';
import * as D from '../js/core/bridge/decks.js';
import * as R from '../js/core/bridge/results.js';
import * as Z from '../js/core/bridge/zip.js';
import * as GH from '../js/core/bridge/remote.js';
import * as G from '../js/core/geometry/index.js';

let pass = 0, fail = 0, section = '';
const failures = [];
const ok = (name, cond, detail = '') => { if (cond) pass++; else { fail++; failures.push(`[${section}] ${name}${detail ? ' — ' + detail : ''}`); } };
const near = (name, a, b, tol = 1e-9) => ok(name, Number.isFinite(a) && Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `got ${a}, expected ${b} (tol ${tol})`);
const throws = (name, fn, re) => { try { fn(); ok(name, false, 'did not throw'); } catch (e) { ok(name, !re || re.test(e.message), `threw "${e.message}"`); } };
const sec = (s) => { section = s; };
const enc = (s) => new TextEncoder().encode(s);
const file = (b, p) => b.files.find((f) => f.path === p)?.text;
const USE = ['narrowbody', 'turboprop', 'ga', 'fixedUav', 'evtol', 'helicopter'];

// ---------------------------------------------------------------------------------------------
// validators
// ---------------------------------------------------------------------------------------------
/** Parse and check an SU2 native mesh. Returns { ok, why, nDim, nElem, nPoin, markers, minMeasure }. */
function checkSu2(text) {
  const L = text.split('\n'); let i = 0; const next = (key) => { for (; i < L.length; i++) if (L[i].startsWith(key)) return Number(L[i++].slice(key.length)); return NaN; };
  const nDim = next('NDIME='), nElem = next('NELEM='), elems = [];
  for (let e = 0; e < nElem; e++) elems.push(L[i++].trim().split(/\s+/).map(Number));
  const nPoin = next('NPOIN='), P = [];
  for (let p = 0; p < nPoin; p++) P.push(L[i++].trim().split(/\s+/).map(Number).slice(0, nDim));
  const nMark = next('NMARK='), markers = {};
  for (let m = 0; m < nMark; m++) { const tag = L[i++].replace('MARKER_TAG=', '').trim(), n = Number(L[i++].replace('MARKER_ELEMS=', '')), els = []; for (let k = 0; k < n; k++) els.push(L[i++].trim().split(/\s+/).map(Number)); markers[tag] = els; }
  const bad = (why) => ({ ok: false, why, nDim, nElem, nPoin, markers });
  if (![2, 3].includes(nDim) || !(nElem > 0) || !(nPoin > 0) || elems.length !== nElem || P.length !== nPoin) return bad('header counts do not match the content');
  if (P.some((p) => p.length !== nDim || p.some((v) => !Number.isFinite(v)))) return bad('non-finite or short coordinates');
  const used = new Uint8Array(nPoin), faces = new Map(), key = (ids) => ids.slice().sort((a, b) => a - b).join(',');
  let minMeasure = Infinity;
  const sub = (a, b) => a.map((v, k) => v - b[k]), cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]], dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  for (const el of elems) {
    const type = el[0], ids = el.slice(1, type === 9 ? 5 : 9);
    if ((nDim === 2 && type !== 9) || (nDim === 3 && type !== 12)) return bad(`unexpected element type ${type}`);
    if (ids.some((q) => !(q >= 0 && q < nPoin) || !Number.isInteger(q))) return bad('element refers to a node that does not exist');
    if (new Set(ids).size !== ids.length) return bad('degenerate element (repeated node)');
    for (const q of ids) used[q] = 1;
    if (nDim === 2) {
      let a = 0; for (let k = 0; k < 4; k++) { const p = P[ids[k]], q = P[ids[(k + 1) % 4]]; a += p[0] * q[1] - q[0] * p[1]; } a /= 2; if (a < minMeasure) minMeasure = a;
      for (let k = 0; k < 4; k++) { const kk = key([ids[k], ids[(k + 1) % 4]]); faces.set(kk, (faces.get(kk) || 0) + 1); }
    } else {
      // corner Jacobians of the VTK-ordered hexahedron must all be positive
      for (const [c, x, y, z] of [[0, 1, 3, 4], [1, 2, 0, 5], [2, 3, 1, 6], [3, 0, 2, 7], [4, 7, 5, 0], [5, 4, 6, 1], [6, 5, 7, 2], [7, 6, 4, 3]]) { const j = dot(cross(sub(P[ids[x]], P[ids[c]]), sub(P[ids[y]], P[ids[c]])), sub(P[ids[z]], P[ids[c]])); if (j < minMeasure) minMeasure = j; }
      for (const fc of [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]) { const kk = key(fc.map((q) => ids[q])); faces.set(kk, (faces.get(kk) || 0) + 1); }
    }
  }
  if (!(minMeasure > 0)) return bad(`inverted or zero-size cell (smallest measure ${minMeasure})`);
  if (used.some((u) => !u)) return bad('node not used by any element');
  if ([...faces.values()].some((n) => n > 2)) return bad('a face is shared by more than two cells');
  const boundary = new Set([...faces].filter(([, n]) => n === 1).map(([k]) => k)), marked = new Set();
  for (const [tag, els] of Object.entries(markers)) for (const el of els) {
    if (el[0] !== (nDim === 2 ? 3 : 9)) return bad(`marker ${tag}: wrong boundary element type ${el[0]}`);
    const kk = key(el.slice(1)); if (!boundary.has(kk)) return bad(`marker ${tag}: element ${kk} is not on the boundary of the mesh`); if (marked.has(kk)) return bad(`marker ${tag}: boundary element ${kk} listed twice`); marked.add(kk);
  }
  if (marked.size !== boundary.size) return bad(`${boundary.size - marked.size} boundary faces carry no marker`);
  return { ok: true, why: '', nDim, nElem, nPoin, markers, minMeasure };
}
const cfgMap = (text) => { const m = new Map(), dup = []; for (const l of text.split('\n')) { if (!l.trim() || l.startsWith('%')) continue; const k = l.split('=')[0].trim(), v = l.slice(l.indexOf('=') + 1).split('%')[0].trim(); if (m.has(k)) dup.push(k); m.set(k, v); } return { m, dup }; };
const tuple = (v) => (v || '').replace(/[()]/g, '').split(',').map((x) => x.trim()).filter(Boolean);

/** Braces / parentheses balance of an OpenFOAM dictionary, ignoring comments and strings. */
function foamBalanced(text) {
  const s = text.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/.*$/gm, '').replace(/"[^"]*"/g, '""'); const st = [];
  for (const ch of s) { if (ch === '{' || ch === '(') st.push(ch); else if (ch === '}') { if (st.pop() !== '{') return false; } else if (ch === ')') { if (st.pop() !== '(') return false; } }
  return st.length === 0;
}
/** Closedness of an ASCII STL: every directed edge is matched by exactly one opposite edge; volume positive. */
function checkStl(text) {
  const v = [...text.matchAll(/vertex\s+(\S+)\s+(\S+)\s+(\S+)/g)].map((m) => [Number(m[1]), Number(m[2]), Number(m[3])]), ids = new Map(), id = (p) => { const k = p.map((x) => x.toPrecision(9)).join(','); if (!ids.has(k)) ids.set(k, ids.size); return ids.get(k); };
  const edges = new Map(); let vol = 0, degenerate = 0;
  for (let t = 0; t < v.length; t += 3) {
    const a = v[t], b = v[t + 1], c = v[t + 2], q = [id(a), id(b), id(c)];
    if (new Set(q).size < 3) { degenerate++; continue; }
    vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
    for (let k = 0; k < 3; k++) { const e = `${q[k]}>${q[(k + 1) % 3]}`; edges.set(e, (edges.get(e) || 0) + 1); }
  }
  let open = 0; for (const [e, n] of edges) { const [a, b] = e.split('>'); if (n !== 1 || edges.get(`${b}>${a}`) !== 1) open++; }
  return { tris: v.length / 3, open, vol, degenerate };
}
/** Parse a CalculiX / Abaqus deck (keywords, nodes, elements, sets) and check every reference. */
function checkInp(text, include = {}) {
  const NODES = { S3: 3, S4: 4, S6: 6, S8R: 8, C3D4: 4, C3D10: 10, C3D8: 8, C3D20R: 20, C3D6: 6, C3D15: 15 };
  const nodes = new Map(), elems = new Map(), elsets = new Map(), nsets = new Map(), kw = [], sections = [], cloads = [], bounds = []; let mode = null, par = {}, pend = [], err = '';
  const lines = []; for (const l of text.split('\n')) { const m = /^\*INCLUDE,\s*INPUT=(\S+)/i.exec(l); if (m) { if (include[m[1]] == null) return { err: `include ${m[1]} missing` }; lines.push(...include[m[1]].split('\n')); } else lines.push(l); }
  for (const raw of lines) {
    const l = raw.trim(); if (!l || l.startsWith('**')) continue;
    if (l.startsWith('*')) { const parts = l.split(',').map((x) => x.trim()); mode = parts[0].toUpperCase(); kw.push(mode); par = {}; for (const p of parts.slice(1)) { const [k, v] = p.split('='); par[k.trim().toUpperCase()] = (v ?? '').trim(); } pend = []; if (mode === '*SHELL SECTION' || mode === '*SOLID SECTION') sections.push({ ...par, kind: mode }); if (mode === '*ELEMENT' && par.ELSET && !elsets.has(par.ELSET)) elsets.set(par.ELSET, []); continue; }
    const t = l.split(',').map((x) => x.trim()).filter((x) => x !== '');
    if (mode === '*NODE') { const id = Number(t[0]); if (nodes.has(id)) err ||= `duplicate node ${id}`; if (t.length !== 4 || t.slice(1).some((x) => !Number.isFinite(Number(x)))) err ||= `bad node line "${l}"`; nodes.set(id, t.slice(1).map(Number)); }
    else if (mode === '*ELEMENT') { pend.push(...t.map(Number)); const need = NODES[par.TYPE]; if (!need) { err ||= `unknown element type ${par.TYPE}`; continue; } if (pend.length >= need + 1) { const [id, ...ns] = pend; if (ns.length !== need) err ||= `element ${id} has ${ns.length} nodes, ${par.TYPE} needs ${need}`; if (elems.has(id)) err ||= `duplicate element ${id}`; elems.set(id, { type: par.TYPE, nodes: ns }); elsets.get(par.ELSET)?.push(id); pend = []; } }
    else if (mode === '*NSET') { if (!nsets.has(par.NSET)) nsets.set(par.NSET, []); nsets.get(par.NSET).push(...t.map(Number)); }
    else if (mode === '*ELSET') { if (!elsets.has(par.ELSET)) elsets.set(par.ELSET, []); elsets.get(par.ELSET).push(...t.map(Number)); }
    else if (mode === '*CLOAD') cloads.push([Number(t[0]), Number(t[1]), Number(t[2])]);
    else if (mode === '*BOUNDARY') bounds.push(t);
  }
  if (par.NSET && mode === '*NODE') { /* NALL */ }
  for (const [id, e] of elems) { if (new Set(e.nodes).size !== e.nodes.length) err ||= `element ${id} repeats a node`; for (const n of e.nodes) if (!nodes.has(n)) err ||= `element ${id} refers to missing node ${n}`; }
  for (const [nm, ids] of nsets) for (const n of ids) if (!nodes.has(n)) err ||= `node set ${nm} refers to missing node ${n}`;
  for (const [nm, ids] of elsets) for (const n of ids) if (!elems.has(n)) err ||= `element set ${nm} refers to missing element ${n}`;
  for (const s of sections) { if (!elsets.has(s.ELSET)) err ||= `section refers to missing element set ${s.ELSET}`; if (s.MATERIAL !== 'MAT') err ||= 'section without material'; }
  for (const [n, dof, F] of cloads) if (!nodes.has(n) || !(dof >= 1 && dof <= 6) || !Number.isFinite(F)) err ||= `bad load on node ${n}`;
  for (const b of bounds) if (!nsets.has(b[0]) && !nodes.has(Number(b[0]))) err ||= `boundary refers to unknown set ${b[0]}`;
  const usedN = new Set(); for (const e of elems.values()) for (const n of e.nodes) usedN.add(n);
  const sectioned = new Set(sections.flatMap((s) => elsets.get(s.ELSET) || []));
  if (usedN.size !== nodes.size) err ||= `${nodes.size - usedN.size} nodes are not used by any element`;
  if (sectioned.size !== elems.size) err ||= `${elems.size - sectioned.size} elements have no section`;
  return { err, nodes, elems, elsets, nsets, kw, cloads, loadZ: cloads.filter((c) => c[1] === 3).reduce((a, c) => a + c[2], 0) };
}

// ---------------------------------------------------------------------------------------------
sec('hash and geometry primitives');
ok('SHA-256 of "abc"', D.sha256('abc') === 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad');
ok('SHA-256 of the empty string', D.sha256('') === 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855');
ok('SHA-256 across a block boundary', D.sha256('a'.repeat(119)) === D.sha256(enc('a'.repeat(119))) && D.sha256('a'.repeat(64)) === 'ffe054fe7ae0cb6dc65c3af9b61d5209f439851db43d0ba5997337df154668eb');
ok('canonical JSON ignores key order', D.caseHash({ a: 1, b: { c: 2, d: [1, 2] } }) === D.caseHash({ b: { d: [1, 2], c: 2 }, a: 1 }));
near('geometric ratio reproduces the total', (() => { const r = D.geometricRatio(1e-6, 80); return (1e-6 * (r ** 80 - 1)) / (r - 1); })(), 1, 1e-9);
{ const a = D.parseNaca('2412'), p = D.nacaAt(a, 0.3); near('NACA 2412 thickness at 30 % chord', p.upper[1] - p.lower[1], 0.12 * 1.0006, 2e-3); ok('non-NACA code falls back to a symmetric section of the case thickness', D.parseNaca('SC(2)-0714', 0.14).code === '0014' && !D.parseNaca('x', 0.1).exact); }
{ const lp = D.nacaLoop(D.parseNaca('4412'), 50); ok('section loop is closed at a sharp trailing edge', lp.length === 101 && lp[0][0] === 1 && lp[100][0] === 1 && lp[0][1] === 0 && lp[50][0] === 0); }
near('wall spacing for y+ = 1 at Re 6e6 is a few micro-chords', D.wallSpacing(6e6, 1, 1), 4.5e-6, 0.4);

// ---------------------------------------------------------------------------------------------
sec('SU2 grids');
for (const naca of ['0012', '2412', '4412', '4318', '4415', '6409', '0006', '2424']) for (const type of ['C', 'O']) {
  const g = D.airfoilGrid({ naca, type, chord: 1.7, nSurf: 48, nWake: 24, growth: 1.2, firstCell: 8e-6, radius: 30, wake: 20 }), c = checkSu2(D.su2MeshText(g));
  ok(`${type}-grid NACA ${naca}: valid, positive cells, boundary fully marked`, c.ok, c.why);
  ok(`${type}-grid NACA ${naca}: counts`, c.nElem === g.info.cells && c.nPoin === g.info.points && c.markers.airfoil?.length === 96 && c.markers.farfield?.length > 0);
  near(`${type}-grid NACA ${naca}: first node height`, g.info.firstCell, 8e-6, 0.03);
}
{ const g = D.extrudeGrid(D.airfoilGrid({ naca: '2412', type: 'C', nSurf: 24, nWake: 12, firstCell: 1e-3, radius: 20 }), 2, 4), c = checkSu2(D.su2MeshText(g));
  ok('extruded C-grid: valid hexahedra and four markers', c.ok && c.nDim === 3 && Object.keys(c.markers).join() === 'airfoil,farfield,sym_root,sym_tip', c.why);
  const m = await G.importFile('extruded.su2', enc(D.su2MeshText(g)));
  ok('extruded grid is read back by the app\'s own SU2 importer', m.elements.some((e) => e.type === 'hex8' && e.count === g.info.cells) && m.groups.map((x) => x.name).includes('sym_tip'));
}

// ---------------------------------------------------------------------------------------------
sec('SU2 cases');
const su2Variants = [{ analysis: 'rans-sa' }, { analysis: 'rans-sst', regime: 'compressible' }, { analysis: 'euler', regime: 'incompressible' }, { analysis: 'euler', regime: 'compressible', mesh: 'O-grid', mg_levels: 3 }, { analysis: 'rans-sa', regime: 'incompressible', mesh: 'O-grid' }, { analysis: 'rans-sa', dim: '3-D extruded' }];
for (const preset of USE) for (const v of su2Variants) {
  const c = PRESETS[preset].data(), b = D.buildCase(c, { solver: 'su2', settings: { ...v, resolution: 'coarse' }, appVersion: 'test' }), tag = `${preset} ${JSON.stringify(v)}`;
  const mesh = checkSu2(file(b, 'mesh.su2')), { m, dup } = cfgMap(file(b, 'config.cfg')), s = b.settings, three = v.dim === '3-D extruded';
  ok(`${tag}: mesh valid`, mesh.ok && mesh.nDim === (three ? 3 : 2), mesh.why);
  ok(`${tag}: no duplicate options`, dup.length === 0, dup.join());
  const names = Object.keys(mesh.markers), wallKey = v.analysis === 'euler' ? 'MARKER_EULER' : 'MARKER_HEATFLUX', walls = tuple(m.get(wallKey)).filter((x) => Number.isNaN(Number(x)));
  const cfgMarkers = [...walls, ...tuple(m.get('MARKER_FAR')), ...tuple(m.get('MARKER_SYM'))];
  ok(`${tag}: every mesh marker has a boundary condition and vice versa`, names.every((n) => cfgMarkers.includes(n)) && cfgMarkers.every((n) => names.includes(n)) && new Set(cfgMarkers).size === cfgMarkers.length, `${names} vs ${cfgMarkers}`);
  ok(`${tag}: monitored and plotted markers exist`, [...tuple(m.get('MARKER_MONITORING')), ...tuple(m.get('MARKER_PLOTTING'))].every((n) => names.includes(n)));
  const inc = s.regime === 'incompressible', expectSolver = inc ? (v.analysis === 'euler' ? 'INC_EULER' : 'INC_RANS') : v.analysis === 'euler' ? 'EULER' : 'RANS';
  ok(`${tag}: solver ${expectSolver}`, m.get('SOLVER') === expectSolver && (v.analysis === 'euler' ? !m.has('KIND_TURB_MODEL') : m.get('KIND_TURB_MODEL') === (v.analysis === 'rans-sst' ? 'SST' : 'SA')));
  for (const k of ['MESH_FILENAME', 'MESH_FORMAT', 'ITER', 'CFL_NUMBER', 'REF_AREA', 'REF_LENGTH', 'CONV_NUM_METHOD_FLOW', 'LINEAR_SOLVER', 'CONV_FIELD', 'HISTORY_OUTPUT', 'OUTPUT_FILES', 'MGLEVEL', 'TIME_DISCRE_FLOW']) ok(`${tag}: has ${k}`, m.has(k) && m.get(k) !== '');
  if (inc) { const vel = tuple(m.get('INC_VELOCITY_INIT')).map(Number), fsx = D.flowState(c); near(`${tag}: velocity magnitude from the case`, Math.hypot(...vel), fsx.V, 1e-6); near(`${tag}: density from ISA`, Number(m.get('INC_DENSITY_INIT')), fsx.rho, 1e-4); ok(`${tag}: lift direction is ${three ? 'z' : 'y'}`, three ? vel[1] === 0 : vel[2] === 0); if (v.analysis !== 'euler') near(`${tag}: viscosity gives the case Reynolds number`, (fsx.rho * fsx.V * s.chord_m) / Number(m.get('MU_CONSTANT')), fsx.reynolds, 2e-3); }
  else { const fsx = D.flowState(c); near(`${tag}: Mach from speed and ISA speed of sound`, Number(m.get('MACH_NUMBER')), fsx.mach, 1e-4); near(`${tag}: angle of attack`, Number(m.get('AOA')), c.flight.alpha_deg, 1e-9); near(`${tag}: temperature`, Number(m.get('FREESTREAM_TEMPERATURE')), fsx.T, 1e-5); if (v.analysis !== 'euler') near(`${tag}: Reynolds number`, Number(m.get('REYNOLDS_NUMBER')), fsx.reynolds, 1e-4); }
  near(`${tag}: reference area`, Number(m.get('REF_AREA')), three ? s.chord_m * s.span_m : s.chord_m, 1e-6);
  if (v.analysis !== 'euler' && !three) near(`${tag}: first cell matches the y+ target`, b.meshInfo.firstCell, D.wallSpacing(s.reynolds, s.chord_m, 1), 0.03);
  ok(`${tag}: run script, readme, case and manifest present`, ['run.sh', 'README.txt', 'case.json', 'manifest.json'].every((p) => file(b, p) != null) && b.files.find((f) => f.path === 'run.sh').exec === true);
  const man = JSON.parse(file(b, 'manifest.json'));
  ok(`${tag}: manifest traces every deck by checksum`, man.caseHash === D.caseHash(c) && man.appVersion === 'test' && man.decks.length === b.files.length - 1 && man.decks.every((d) => d.sha256 === D.sha256(file(b, d.path))) && man.expectedOutputs.includes('history.csv'));
}
throws('imported-mesh route without a mesh is refused in plain words', () => D.buildCase(PRESETS.ga.data(), { solver: 'su2', settings: { mesh: 'imported' } }), /volume mesh/);
throws('out-of-range setting is refused', () => D.buildCase(PRESETS.ga.data(), { solver: 'su2', settings: { mach: -1 } }), /Mach number/);
{ // imported volume mesh → SU2 deck
  const src = D.su2MeshText(D.extrudeGrid(D.airfoilGrid({ naca: '0012', type: 'O', nSurf: 16, firstCell: 5e-3, radius: 10 }), 1, 2)), model = await G.importFile('wing.su2', enc(src));
  const def = D.defaultSettings(PRESETS.turboprop.data(), 'su2', { model });
  ok('imported mesh: roles guessed from boundary names', def.mesh === 'imported' && def.wall_markers === 'airfoil' && def.far_markers === 'farfield' && def.sym_markers === 'sym_root, sym_tip', JSON.stringify([def.wall_markers, def.far_markers, def.sym_markers]));
  const b = D.buildCase(PRESETS.turboprop.data(), { solver: 'su2', model, settings: { analysis: 'euler' } }), chk = checkSu2(file(b, 'mesh.su2')), { m } = cfgMap(file(b, 'config.cfg'));
  ok('imported mesh: exported mesh is valid and fully marked', chk.ok, chk.why);
  ok('imported mesh: configuration uses its markers and the wing reference area', m.get('MARKER_EULER') === '( airfoil )' && m.get('MARKER_SYM') === '( sym_root, sym_tip )' && Number(m.get('REF_AREA')) === 61);
  throws('imported mesh: an unassigned boundary is reported', () => D.buildCase(PRESETS.turboprop.data(), { solver: 'su2', model, settings: { analysis: 'euler', sym_markers: '' } }), /sym_root/);
}

// ---------------------------------------------------------------------------------------------
sec('OpenFOAM cases');
const FOAM_REQ = { 'system/controlDict': ['application', 'endTime', 'deltaT', 'writeInterval', 'functions', 'forceCoeffs', 'liftDir', 'dragDir', 'magUInf', 'lRef', 'Aref', 'rhoInf'], 'system/fvSchemes': ['ddtSchemes', 'gradSchemes', 'divSchemes', 'laplacianSchemes', 'interpolationSchemes', 'snGradSchemes', 'wallDist'], 'system/fvSolution': ['solvers', 'tolerance', 'relTol'], 'system/blockMeshDict': ['vertices', 'blocks', 'boundary', 'freestream'], 'system/snappyHexMeshDict': ['castellatedMesh', 'geometry', 'body.stl', 'triSurfaceMesh', 'refinementSurfaces', 'locationInMesh', 'snapControls', 'addLayersControls', 'meshQualityControls', 'mergeTolerance'], 'constant/transportProperties': ['transportModel', 'nu'], 'constant/turbulenceProperties': ['simulationType'], '0/U': ['dimensions', 'internalField', 'boundaryField', 'freestream', 'noSlip'], '0/p': ['dimensions', 'internalField', 'boundaryField'], '0/nut': ['nutUSpaldingWallFunction'], '0/nuTilda': ['boundaryField'], '0/k': ['kqRWallFunction'], '0/omega': ['omegaWallFunction'] };
for (const preset of USE) for (const analysis of ['simpleFoam', 'pimpleFoam', 'rhoSimpleFoam']) for (const turbulence of analysis === 'simpleFoam' ? ['kOmegaSST', 'SpalartAllmaras'] : ['kOmegaSST']) {
  const c = PRESETS[preset].data(), b = D.buildCase(c, { solver: 'openfoam', settings: { analysis, turbulence } }), tag = `${preset} ${analysis} ${turbulence}`, s = b.settings;
  for (const [p, keys] of Object.entries(FOAM_REQ)) { const t = file(b, p); ok(`${tag}: ${p} exists with FoamFile header`, !!t && /FoamFile\s*\{[^}]*object\s+\w+;/.test(t)); ok(`${tag}: ${p} balanced`, !!t && foamBalanced(t)); for (const k of keys) ok(`${tag}: ${p} has ${k}`, !!t && t.includes(k)); }
  const cd = file(b, 'system/controlDict'), tp = file(b, 'constant/turbulenceProperties'), fs2 = file(b, 'system/fvSolution');
  ok(`${tag}: application`, new RegExp(`application\\s+${analysis};`).test(cd));
  ok(`${tag}: turbulence model`, analysis === 'pimpleFoam' ? /simulationType LES;/.test(tp) && /LESModel\s+WALE;/.test(tp) && /delta\s+cubeRootVol;/.test(tp) : new RegExp(`RASModel\\s+${turbulence};`).test(tp));
  ok(`${tag}: pressure-velocity coupling`, analysis === 'pimpleFoam' ? /PIMPLE\s*\{/.test(fs2) && /pFinal/.test(fs2) && /ddtSchemes\s*\{\s*default\s+backward;/.test(file(b, 'system/fvSchemes')) : /SIMPLE\s*\{/.test(fs2) && /steadyState/.test(file(b, 'system/fvSchemes')));
  if (analysis === 'rhoSimpleFoam') { ok(`${tag}: thermophysical model and temperature field`, /hePsiThermo/.test(file(b, 'constant/thermophysicalProperties') || '') && /dimensions\s+\[0 0 0 1 0 0 0\]/.test(file(b, '0/T') || '') && foamBalanced(file(b, '0/alphat') || '(')); ok(`${tag}: temperature limiter present and balanced`, /limitTemperature/.test(file(b, 'system/fvOptions') || '') && foamBalanced(file(b, 'system/fvOptions'))); near(`${tag}: absolute pressure from ISA`, Number(/internalField\s+uniform\s+(\S+);/.exec(file(b, '0/p'))[1]), D.flowState(c).p, 1e-5); }
  else ok(`${tag}: kinematic pressure`, /dimensions\s+\[0 2 -2 0 0 0 0\]/.test(file(b, '0/p')) && file(b, 'constant/thermophysicalProperties') == null);
  const U = /internalField\s+uniform\s+\(([^)]+)\)/.exec(file(b, '0/U'))[1].split(/\s+/).map(Number);
  near(`${tag}: speed`, Math.hypot(...U), c.flight.V_ms, 1e-6); near(`${tag}: angle of attack in the x–z plane`, Math.atan2(U[2], U[0]) * 180 / Math.PI, c.flight.alpha_deg, 1e-6);
  const lift = /liftDir\s+\(([^)]+)\)/.exec(cd)[1].split(/\s+/).map(Number), drag = /dragDir\s+\(([^)]+)\)/.exec(cd)[1].split(/\s+/).map(Number);
  near(`${tag}: lift and drag directions are orthonormal and drag follows the flow`, Math.abs(lift[0] * drag[0] + lift[2] * drag[2]) + Math.abs(Math.hypot(...lift) - 1) + Math.abs(drag[0] * U[2] - drag[2] * U[0]) / c.flight.V_ms, 0, 1e-6);
  near(`${tag}: reference area`, Number(/Aref\s+(\S+);/.exec(cd)[1]), s.ref_area_m2, 1e-6);
  const stl = checkStl(file(b, 'constant/triSurface/body.stl'));
  ok(`${tag}: surface is closed, outward and non-degenerate (${stl.tris} facets)`, stl.tris > 200 && stl.open === 0 && stl.vol > 0 && stl.degenerate === 0, JSON.stringify(stl));
  const bm = file(b, 'system/blockMeshDict'), verts = [...bm.matchAll(/^\s{4}\(([-\d.e+ ]+)\)$/gm)].map((x) => x[1].split(/\s+/).map(Number)), cells = /\(0 1 2 3 4 5 6 7\) \((\d+) (\d+) (\d+)\)/.exec(bm).slice(1).map(Number), mi = b.meshInfo, d0 = mi.coreCell_m;
  const grad = [...bm.matchAll(/^\s{8}\(\(([^)]+)\) \(([^)]+)\) \(([^)]+)\)\)$/gm)].map((m) => m.slice(1).map((t) => t.split(/\s+/).map(Number)));
  ok(`${tag}: background block has 8 vertices spanning the stated domain`, verts.length === 8 && [0, 1, 2].every((k) => { const tol = 1e-6 * (mi.domain.max[k] - mi.domain.min[k]); return Math.abs(verts[0][k] - mi.domain.min[k]) < tol && Math.abs(verts[6][k] - mi.domain.max[k]) < tol; }) && cells.join() === mi.background.join());
  ok(`${tag}: three-part grading per axis, fractions sum to one, uniform core`, grad.length === 3 && grad.every((g) => Math.abs(g[0][0] + g[1][0] + g[2][0] - 1) < 1e-6 && Math.abs(g[0][1] + g[1][1] + g[2][1] - 1) < 1e-6 && g[1][2] === 1 && g[0][2] < 1 && g[2][2] > 1));
  ok(`${tag}: core cells are cubes of the stated size`, [0, 1, 2].every((k) => { const L = mi.domain.max[k] - mi.domain.min[k], n = grad[k][1][1] * cells[k], len = grad[k][1][0] * L; return Math.abs(n - Math.round(n)) < 1e-5 && Math.abs(len / Math.round(n) - d0) < 1e-5 * d0 && Math.abs(len - (mi.core.max[k] - mi.core.min[k])) < 1e-5 * len; }));
  ok(`${tag}: stretched cells start at the core size on both sides`, [0, 1, 2].every((k) => { const L = mi.domain.max[k] - mi.domain.min[k]; return [0, 2].every((q) => { const len = grad[k][q][0] * L, n = Math.round(grad[k][q][1] * cells[k]), E = q === 0 ? 1 / grad[k][q][2] : grad[k][q][2], r = E ** (1 / (n - 1)), first = (len * (r - 1)) / (r ** n - 1); return Math.abs(first - d0) < 0.02 * d0 && r < 1.31; }); }));
  const feat = s.geometry === 'wing' ? D.flowState(c).d.c_tip : c.fuselage.dia_m;
  ok(`${tag}: core cells are small enough for snappyHexMesh to find the body (≤ half its smallest planform dimension)`, d0 <= 0.5 * feat * 1.0001 && mi.surfaceLevel >= 1 && Math.abs(mi.surfaceCell_m * 2 ** mi.surfaceLevel - d0) < 1e-9 * d0, `${d0} vs ${feat}`);
  const loc = /locationInMesh\s+\(([^)]+)\)/.exec(file(b, 'system/snappyHexMeshDict'))[1].split(/\s+/).map(Number), pts = [...file(b, 'constant/triSurface/body.stl').matchAll(/vertex\s+(\S+)\s+(\S+)\s+(\S+)/g)].map((m) => m.slice(1).map(Number));
  const bmin = [0, 1, 2].map((k) => Math.min(...pts.map((p) => p[k]))), bmax = [0, 1, 2].map((k) => Math.max(...pts.map((p) => p[k])));
  ok(`${tag}: mesh seed point is in the first core cell, outside the body and off the cell faces`, loc.every((v, k) => { const o = (v - mi.core.min[k]) / d0; return o > 0.3 && o < 0.5; }) && loc.every((v, k) => v < bmin[k] - d0));
  ok(`${tag}: body lies inside the uniform core with at least two cells to spare`, [0, 1, 2].every((k) => bmin[k] - mi.core.min[k] >= 2 * d0 && mi.core.max[k] - bmax[k] >= 2 * d0));
  ok(`${tag}: surface refinement levels match the stated level`, new RegExp(`level \\(${Math.max(1, mi.surfaceLevel - 1)} ${mi.surfaceLevel}\\);`).test(file(b, 'system/snappyHexMeshDict')));
  ok(`${tag}: run script meshes then solves`, /blockMesh[\s\S]*snappyHexMesh -overwrite[\s\S]*checkMesh/.test(file(b, 'run.sh')) && file(b, 'run.sh').includes(analysis));
}
{ // imported surface → OpenFOAM
  const stlSrc = (() => { const tris = D.fuselageTriangles(PRESETS.helicopter.data()); let s = 'solid h\n'; for (const t of tris) s += `facet normal 0 0 0\n outer loop\n${t.map((p) => `  vertex ${p.join(' ')}`).join('\n')}\n endloop\nendfacet\n`; return s + 'endsolid h\n'; })();
  const model = await G.importFile('pod.stl', enc(stlSrc)), def = D.defaultSettings(PRESETS.helicopter.data(), 'openfoam', { model }), b = D.buildCase(PRESETS.helicopter.data(), { solver: 'openfoam', model });
  const st = checkStl(file(b, 'constant/triSurface/body.stl'));
  ok('imported surface is used as the body', def.geometry === 'imported' && st.open === 0 && st.vol > 0 && /^solid body/.test(file(b, 'constant/triSurface/body.stl')));
  throws('wing route without a wing is refused', () => D.buildCase(PRESETS.helicopter.data(), { solver: 'openfoam', settings: { geometry: 'wing' } }), /no wing/);
}

// ---------------------------------------------------------------------------------------------
sec('CalculiX / Abaqus cases');
for (const preset of USE) for (const element of ['S8R', 'S4']) for (const analysis of ['all', 'static', 'buckle', 'modal', 'dynamic', 'explicit']) {
  const c = PRESETS[preset].data(), b = D.buildCase(c, { solver: 'calculix', settings: { analysis, element, n_span: 12, n_chord: 3 } }), tag = `${preset} ${element} ${analysis}`, s = b.settings, model = file(b, 'model.inp');
  const jobs = analysis === 'all' ? ['static', 'buckle', 'modal'] : [analysis], mi = b.meshInfo;
  ok(`${tag}: one deck per analysis`, jobs.every((j) => file(b, `${j}.inp`) != null) && b.files.filter((f) => /\.inp$/.test(f.path)).length === jobs.length + 1);
  for (const j of jobs) {
    const r = checkInp(file(b, `${j}.inp`), { 'model.inp': model });
    ok(`${tag}/${j}: every node, element, set, section and load reference resolves`, !r.err, r.err);
    ok(`${tag}/${j}: element type and counts`, [...r.elems.values()].every((e) => e.type === element) && r.elems.size === mi.elements && r.nodes.size === mi.nodes);
    ok(`${tag}/${j}: sets`, ['SKIN_U', 'SKIN_L', 'SPAR_F', 'SPAR_R', 'RIBS'].every((n) => r.elsets.get(n)?.length > 0) && r.nsets.get('ROOT').length > 0 && r.nsets.get('TIP').length > 0 && r.nsets.get('TIPREF').length === 1);
    ok(`${tag}/${j}: root nodes lie on y = 0 and tip nodes at the half span`, r.nsets.get('ROOT').every((n) => Math.abs(r.nodes.get(n)[1]) < 1e-9) && r.nsets.get('TIP').every((n) => Math.abs(r.nodes.get(n)[1] - mi.half) < 1e-9 * mi.half));
    const kws = r.kw.join(' '), step = { static: '*STATIC', buckle: '*BUCKLE', modal: '*FREQUENCY', dynamic: '*DYNAMIC', explicit: '*DYNAMIC' }[j];
    ok(`${tag}/${j}: material, clamp and one complete step`, kws.includes('*MATERIAL *ELASTIC') && kws.includes('*DENSITY') && kws.includes('*BOUNDARY') && r.kw.filter((k) => k === '*STEP').length === 1 && r.kw.indexOf(step) > r.kw.indexOf('*STEP') && r.kw[r.kw.length - 1] === '*END STEP');
    if (j !== 'modal') { near(`${tag}/${j}: applied load equals the lift minus the root half-cell`, r.loadZ, mi.appliedLoad_N, 1e-6); near(`${tag}/${j}: mesh load + root share = load factor × 1 g lift`, mi.appliedLoad_N + mi.rootReaction_N, s.load_factor * s.lift_1g_N, 1e-9); ok(`${tag}/${j}: load is most of the lift`, mi.appliedLoad_N / mi.halfLift_N > 0.9 && mi.appliedLoad_N / mi.halfLift_N < 1); }
    else ok(`${tag}/${j}: no load in the frequency step`, r.cloads.length === 0);
    if (j === 'explicit') ok(`${tag}/${j}: explicit keyword and pulse amplitude`, /\*DYNAMIC, EXPLICIT/.test(file(b, 'explicit.inp')) && /\*AMPLITUDE, NAME=GUST/.test(file(b, 'explicit.inp')) && /\*CLOAD, AMPLITUDE=GUST/.test(file(b, 'explicit.inp')));
    if (j === 'dynamic') ok(`${tag}/${j}: implicit direct time stepping`, /\*DYNAMIC, DIRECT/.test(file(b, 'dynamic.inp')) && !/EXPLICIT/.test(file(b, 'dynamic.inp')));
  }
  ok(`${tag}: material constants from the materials table`, model.includes(String(Number((c.struct.material === 'Al 7075-T6' ? 71.7e9 : 73.1e9).toPrecision(8)))));
  ok(`${tag}: run script runs each job`, jobs.every((j) => file(b, 'run.sh').includes(j)));
}
{ // quarter-chord resultant and box dimensions
  const c = PRESETS.ga.data(), m = D.wingBoxMesh(c, D.defaultSettings(c, 'calculix')), d = m.info;
  let F = 0, Mx = 0; for (const [n, f] of m.loads) { const eta = m.nodes[n][1] / d.half, ch = d.c_root + (d.c_tip - d.c_root) * eta, xle = 0.25 * d.c_root - 0.25 * ch; F += f; Mx += f * ((m.nodes[n][0] - xle) / ch); }
  near('wing box: load resultant acts at 25 % chord', Mx / F, 0.25, 1e-9);
  near('wing box: root cell is the native rectangular cell', d.box_w_root / d.c_root, c.struct.box_chord_frac, 1e-12); near('wing box: depth', d.box_h_root / d.c_root, c.struct.box_height_frac * c.wing.tc, 1e-12);
  // imported structural mesh (the generated S4 deck read back through the app's Abaqus importer)
  const b4 = D.buildCase(c, { solver: 'calculix', settings: { analysis: 'static', element: 'S4', n_span: 8, n_chord: 2 } }), model = await G.importFile('box.inp', enc(file(b4, 'model.inp')));
  ok('generated deck is read back by the app\'s own Abaqus importer', model.elements.some((e) => e.type === 'quad4' && e.count === b4.meshInfo.elements));
  const bi = D.buildCase(c, { solver: 'calculix', model, settings: { analysis: 'all', source: 'imported' } }), r = checkInp(file(bi, 'static.inp'), { 'model.inp': file(bi, 'model.inp') });
  ok('imported structural mesh: valid deck', !r.err && r.elems.size === b4.meshInfo.elements && [...r.elems.values()].every((e) => e.type === 'S4'), r.err);
  near('imported structural mesh: total load', r.loadZ, bi.settings.load_factor * bi.settings.lift_1g_N, 1e-6);
  ok('imported structural mesh: clamped at the low end of its longest axis', bi.meshInfo.clampAxis === 'y' && r.nsets.get('ROOT').every((n) => r.nodes.get(n)[1] < 0.06));
  throws('imported route without a mesh is refused', () => D.buildCase(c, { solver: 'calculix', settings: { source: 'imported' } }), /structural mesh/);
}

// ---------------------------------------------------------------------------------------------
sec('Gmsh script');
for (const preset of ['narrowbody', 'turboprop', 'ga', 'fixedUav', 'evtol']) {
  const geo = D.gmshGeo(PRESETS[preset].data());
  ok(`${preset}: balanced braces and required statements`, foamBalanced(geo.replace(/\[\]/g, '')) && ['SetFactory("OpenCASCADE")', 'ThruSections', 'BooleanDifference', 'Physical Surface("wall")', 'Physical Surface("symmetry")', 'Physical Surface("farfield")', 'Physical Volume("fluid")', 'Background Field'].every((k) => geo.includes(k)));
  const pts = [...geo.matchAll(/^Point\((\d+)\) = \{([^}]+)\};/gm)], splines = [...geo.matchAll(/^Spline\(\d+\) = \{([^}]+)\};/gm)].flatMap((m) => m[1].split(',').map(Number));
  ok(`${preset}: three sections, every spline point defined`, pts.length === 3 * 60 && splines.every((p) => p >= 1 && p <= pts.length) && pts.every((m) => m[2].split(',').every((x) => Number.isFinite(Number(x)))));
  ok(`${preset}: included in the CFD archives`, D.buildCase(PRESETS[preset].data(), { solver: 'su2', settings: { resolution: 'coarse' } }).files.some((f) => f.path === 'geometry/aircraft.geo'));
}
throws('Gmsh script needs a wing', () => D.gmshGeo(PRESETS.multirotor.data()), /wing/);

// ---------------------------------------------------------------------------------------------
sec('ZIP');
{
  const big = 'The quick brown fox jumps over the lazy dog. '.repeat(4000), bin = Uint8Array.from({ length: 5000 }, (_, i) => (i * 7919) & 255);
  const files = [{ path: 'config.cfg', text: 'A= 1\n' }, { path: 'system/controlDict', text: big }, { path: 'data/raw.bin', bytes: bin }, { path: 'run.sh', text: '#!/bin/sh\n', exec: true }, { path: 'empty.txt', text: '' }, { path: 'ünï/cödé.txt', text: 'π' }];
  near('CRC-32 of "123456789"', Z.crc32(enc('123456789')), 0xcbf43926, 0);
  for (const deflate of [true, false]) {
    const z = await Z.zipWrite(files, { deflate, date: new Date(2026, 0, 2, 3, 4, 6) }), back = await Z.zipRead(z);
    ok(`round trip (${deflate ? 'deflate' : 'stored'}): names`, back.map((f) => f.name).join('|') === files.map((f) => f.path).join('|'));
    ok(`round trip (${deflate ? 'deflate' : 'stored'}): contents`, back.every((f, i) => Buffer.compare(Buffer.from(f.bytes), Buffer.from(files[i].bytes ?? enc(files[i].text))) === 0));
    ok(`${deflate ? 'deflate shrinks' : 'stored keeps'} repetitive text`, deflate ? z.length < big.length / 10 : z.length > big.length);
    const tmp = path.join(os.tmpdir(), `bridge-zip-${process.pid}-${deflate ? 'd' : 's'}.zip`);
    try { fs.writeFileSync(tmp, z); const r = spawnSync('unzip', ['-t', tmp], { encoding: 'utf8' }); if (!r.error) ok(`archive passes "unzip -t" (${deflate ? 'deflate' : 'stored'})`, r.status === 0 && /No errors detected/.test(r.stdout), r.stdout.slice(-200)); } finally { fs.rmSync(tmp, { force: true }); }
    const broken = z.slice(); broken[30 + 'config.cfg'.length + 1] ^= 0xff;                // a byte of the first entry's data
    let caught = ''; try { await Z.zipRead(broken); } catch (e) { caught = e.message; } ok(`corruption is detected (${deflate ? 'deflate' : 'stored'})`, /checksum|damaged/.test(caught), caught);
  }
  for (const p of ['../evil', '/abs/path', 'C:\\x', 'a/../../b']) { let msg = ''; try { await Z.zipWrite([{ path: p, text: 'x' }]); } catch (e) { msg = e.message; } ok(`unsafe path "${p}" refused`, /Unsafe path/.test(msg)); }
  let msg = ''; try { await Z.zipRead(enc('not a zip at all')); } catch (e) { msg = e.message; } ok('non-archive rejected in plain words', /not a ZIP/.test(msg));
  near('base64 round trip', Z.base64Decode(Z.base64Encode(bin)).reduce((a, v, i) => a + Math.abs(v - bin[i]), 0), 0, 0);
  const b = D.buildCase(PRESETS.ga.data(), { solver: 'calculix', settings: { n_span: 8 } }), z = await Z.zipWrite(b.files.map((f) => ({ path: f.path, text: f.text, exec: f.exec }))), back = await Z.zipRead(z);
  ok('a whole case survives the archive byte for byte', back.length === b.files.length && back.every((f) => Z.toText(f.bytes) === file(b, f.name)));
}

// ---------------------------------------------------------------------------------------------
sec('result readers');
// --- SU2 history.csv (rows 0, 1, 250, 498, 499 of a real incompressible RANS-SA run of the generated NACA 2412 C-grid, SU2 8.5.0)
const HIST = `"Time_Iter","Outer_Iter","Inner_Iter",     "rms[P]"     ,     "rms[U]"     ,     "rms[V]"     ,     "rms[nu]"    ,    "RefForce"    ,       "CD"       ,       "CL"       ,       "CSF"      ,       "CMx"      ,       "CMy"      ,       "CMz"      ,       "CFx"      ,       "CFy"      ,       "CFz"      ,      "CEff"
          0,           0,           0,       -3.47617665,      -15.84173657,      -17.17550729,      -8.735906538,          0.744625,      0.6180573699,      0.3667392498,                 0,                 0,                 0,      0.2169122995,      0.6048818508,      0.3880857335,                 0,       0.593374123
          0,           0,           1,      -3.607091525,      -3.484616405,      -3.360415274,      -9.059064677,          0.744625,      0.4869051969,       0.302087767,                 0,                 0,                 0,      0.1759010231,      0.4760658763,      0.3188964897,                 0,      0.6204241995
          0,           0,         250,        -5.5061534,      -5.101510517,      -5.919688848,      -8.887835482,          0.744625,     0.01264482455,      0.4503833289,                 0,                 0,                 0,     0.05089349713,   -0.003081029937,      0.4505502656,                 0,       35.61799747
          0,           0,         498,      -3.905870208,      -4.054123434,      -5.214058435,      -8.524827761,          0.744625,     0.01198948931,       0.438108938,                 0,                 0,                 0,     0.04838946303,   -0.003307595894,       0.438260481,                 0,       36.54108415
          0,           0,         499,      -3.943648509,      -4.056276731,      -5.212554529,       -8.52358346,          0.744625,     0.01198774439,       0.438057899,                 0,                 0,                 0,     0.04837788044,   -0.003307558526,      0.4382094122,                 0,       36.54214545
`;
{
  const r = R.parseSu2History(HIST), v = r.values;
  near('SU2 history: final CL', v.CL, 0.438057899, 1e-12); near('SU2 history: final CD', v.CD, 0.01198774439, 1e-12); near('SU2 history: L/D', v.LD, 36.5421, 1e-5);
  near('SU2 history: 2-D pitching moment is reported nose-up positive', v.Cm, -0.04837788044, 1e-12); near('SU2 history: residual drop', v.residual_drop_decades, 0.467471859, 1e-8);
  ok('SU2 history: KPIs, three charts and a table in app shapes', r.kpis.every((k) => typeof k.key === 'string' && typeof k.label === 'string' && typeof k.value === 'number') && r.plots.length === 3 && r.plots.every((p) => p.type === 'line' && p.title && p.xlabel && p.ylabel && p.series.every((s) => s.x.length === s.y.length)) && r.tables[0].rows.some((x) => x[0] === 'CEff'));
  ok('SU2 history: unconverged run is flagged', r.kpis.find((k) => k.key === 'residual_drop_decades').status === 'warn' && r.warnings.length === 1);
  const long = HIST.split('\n')[0] + '\n' + Array.from({ length: 3000 }, (_, i) => `0,0,${i},${-3 - i / 500},-4,-4,-9,1,${0.01 + 1 / (i + 1)},${0.5 - 1 / (i + 1)},0,0,0,0.05,0,0.5,0,40`).join('\n');
  const rl = R.parseSu2History(long); ok('SU2 history: long histories are thinned to 400 points and a converged run passes', rl.plots[0].series[0].x.length === 400 && rl.plots[0].series[0].x[399] === 2999 && rl.kpis[0].status === 'ok' && rl.warnings.length === 0);
  throws('SU2 history: wrong file is refused', () => R.parseSu2History('a,b\n1,2\n'), /history/);
  const sr = R.parseSu2SurfaceCsv('"PointID","x","y","Pressure","Velocity_x","Velocity_y","Nu_Tilde"\n2368, 1.48925, 0.0, 0.1186, 0, 0, 0\n2442, 1.4883, -7.0e-05, 0.1150, 0, 0, 0\n9, 0.0, 0.0, 0.5, 0, 0, 0\n');
  ok('SU2 surface CSV: pressure plotted against x/c, with a hint to add the VTU', sr.plots[0].series[0].x[0] === 1 && sr.plots[0].series[0].x[2] === 0 && /vtu/.test(sr.warnings[0]));
}
// --- VTU: ASCII and raw appended binary (the layout SU2 writes)
{
  const pts = [[1, 0, 0], [0.5, 0.06, 0], [0, 0, 0], [0.5, -0.06, 0]], cp = [0.9, -0.7, 1.0, -0.2], cf = [[0.001, 0, 0], [0.004, 0, 0], [0, 0, 0], [0.003, 0, 0]], yp = [0.2, 0.9, 0.1, 0.8], conn = [0, 1, 1, 2, 2, 3, 3, 0], offs = [2, 4, 6, 8], types = [3, 3, 3, 3];
  const ascii = `<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">\n<UnstructuredGrid>\n<Piece NumberOfPoints="4" NumberOfCells="4">\n<Points>\n<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n${pts.flat().join(' ')}\n</DataArray>\n</Points>\n<Cells>\n<DataArray type="Int32" Name="connectivity" format="ascii">${conn.join(' ')}</DataArray>\n<DataArray type="Int32" Name="offsets" format="ascii">${offs.join(' ')}</DataArray>\n<DataArray type="UInt8" Name="types" format="ascii">${types.join(' ')}</DataArray>\n</Cells>\n<PointData>\n<DataArray type="Float32" Name="Pressure_Coefficient" NumberOfComponents="1" format="ascii">${cp.join(' ')}</DataArray>\n<DataArray type="Float32" Name="Skin_Friction_Coefficient" NumberOfComponents="3" format="ascii">${cf.flat().join(' ')}</DataArray>\n<DataArray type="Float32" Name="Y_Plus" NumberOfComponents="1" format="ascii">${yp.join(' ')}</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n`;
  const blocks = [[Float32Array.from(pts.flat())], [Int32Array.from(conn)], [Int32Array.from(offs)], [Uint8Array.from(types)], [Float32Array.from(cp)], [Float32Array.from(cf.flat())], [Float32Array.from(yp)]]; let off = 0; const offsets = [], parts = [];
  for (const [a] of blocks) { offsets.push(off); const hd = new DataView(new ArrayBuffer(8)); hd.setBigUint64(0, BigInt(a.byteLength), true); parts.push(new Uint8Array(hd.buffer), new Uint8Array(a.buffer)); off += 8 + a.byteLength; }
  const head = `<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n<UnstructuredGrid>\n<Piece NumberOfPoints="4" NumberOfCells="4">\n<Points>\n<DataArray type="Float32" Name="" NumberOfComponents= "3" offset="${offsets[0]}" format="appended"/>\n</Points>\n<Cells>\n<DataArray type="Int32" Name="connectivity" NumberOfComponents= "1" offset="${offsets[1]}" format="appended"/>\n<DataArray type="Int32" Name="offsets" NumberOfComponents= "1" offset="${offsets[2]}" format="appended"/>\n<DataArray type="UInt8" Name="types" NumberOfComponents= "1" offset="${offsets[3]}" format="appended"/>\n</Cells>\n<PointData>\n<DataArray type="Float32" Name="Pressure_Coefficient" NumberOfComponents= "1" offset="${offsets[4]}" format="appended"/>\n<DataArray type="Float32" Name="Skin_Friction_Coefficient" NumberOfComponents= "3" offset="${offsets[5]}" format="appended"/>\n<DataArray type="Float32" Name="Y_Plus" NumberOfComponents= "1" offset="${offsets[6]}" format="appended"/>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n<AppendedData encoding="raw">\n_`;
  const raw = Buffer.concat([Buffer.from(head, 'latin1'), ...parts.map((p) => Buffer.from(p)), Buffer.from('\n</AppendedData>\n</VTKFile>\n')]);
  for (const [label, input] of [['ASCII', enc(ascii)], ['raw appended', new Uint8Array(raw)]]) {
    const v = R.parseVtu(input), s = R.su2SurfaceFromVtu(v);
    ok(`VTU ${label}: points, cells and three fields`, v.nPoints === 4 && v.conn.length === 8 && v.types[0] === 3 && Object.keys(v.pointData).length === 3 && v.pointData.Skin_Friction_Coefficient.nc === 3);
    near(`VTU ${label}: minimum pressure coefficient`, s.values.cp_min, -0.7, 1e-6); near(`VTU ${label}: largest y+`, s.values.yplus_max, 0.9, 1e-6);
    ok(`VTU ${label}: wall walked in order from the trailing edge, −Cp plotted`, s.plots.length === 3 && s.plots[0].series[0].x.join() === '1,0.5,0,0.5' && Math.abs(s.plots[0].series[0].y[1] - 0.7) < 1e-6 && s.kpis.find((k) => k.key === 'yplus_max').status === 'ok');
  }
  throws('VTU: compressed files are refused in plain words', () => R.parseVtu(ascii.replace('byte_order', 'compressor="vtkZLibDataCompressor" byte_order')), /compressed/);
  // 2-D volume file → contour
  const g = D.airfoilGrid({ naca: '0012', type: 'O', nSurf: 12, firstCell: 0.01, radius: 6 }), np = g.nodes.length / 2, nq = g.quads.length / 4;
  const vol = { nPoints: np, points: Float64Array.from({ length: 3 * np }, (_, i) => (i % 3 === 2 ? 0 : g.nodes[2 * Math.floor(i / 3) + (i % 3)])), conn: Float64Array.from(g.quads), offsets: Float64Array.from({ length: nq }, (_, i) => 4 * (i + 1)), types: new Float64Array(nq).fill(9), pointData: { Mach: { nc: 1, data: Float64Array.from({ length: np }, (_, i) => 0.3 + 0.001 * i) } } };
  const fl = R.su2FieldFromVtu(vol, { chord: 1, xle: 0 });
  ok('VTU volume: near-field Mach contour as a triangle plot', fl.plots.length === 1 && fl.plots[0].type === 'tri' && fl.plots[0].tris.length > 20 && fl.plots[0].tris.every((t) => t.every((i) => i >= 0 && i < fl.plots[0].nodes.length)) && fl.plots[0].values.length === fl.plots[0].nodes.length && fl.values.mach_max > 0.3);
}
// --- OpenFOAM v2412 (trimmed from a real simpleFoam run of the generated wing case)
const COEFF = `# Force and moment coefficients
# dragDir       : (9.99390827e-01 0.00000000e+00 3.48994969e-02)
# magUInf       : 6.20000000e+01
# lRef          : 1.48925000e+00
# Aref          : 1.62000000e+01
#
# Time          \tCd              \tCd(f)           \tCd(r)           \tCl              \tCl(f)           \tCl(r)           \tCmPitch         \tCmRoll          \tCmYaw           \tCs              \tCs(f)           \tCs(r)
1               \t2.83474912e-02\t1.40162152e-02\t1.43312760e-02\t2.45939843e-02\t-7.72216256e-03\t3.23161468e-02\t-2.00191547e-02\t-1.57530439e-04\t-1.00791774e-04\t-3.51267234e-07\t-1.00967408e-04\t1.00616141e-04
2               \t4.51106220e-02\t2.24255612e-02\t2.26850607e-02\t5.43382166e-02\t-6.47540836e-03\t6.08136249e-02\t-3.36445166e-02\t-1.29749737e-04\t-1.26384552e-04\t-4.78490078e-07\t-1.26623797e-04\t1.26145307e-04
79              \t3.01342962e-02\t1.56313045e-02\t1.45029916e-02\t1.68481656e-01\t4.02441938e-02\t1.28237463e-01\t-4.39966344e-02\t5.64156469e-04\t6.80844115e-05\t-4.48825020e-06\t6.58402864e-05\t-7.03285366e-05
80              \t3.01345443e-02\t1.56324247e-02\t1.45021196e-02\t1.68487364e-01\t4.02446017e-02\t1.28242762e-01\t-4.39990801e-02\t5.65152556e-04\t6.80776345e-05\t-4.50075943e-06\t6.58272548e-05\t-7.03280142e-05
`;
const SINFO = `# Solver information
# Time          \tU_solver        \tUx_initial      \tUx_final        \tUx_iters        \tUy_initial      \tUy_final        \tUy_iters        \tUz_initial      \tUz_final        \tUz_iters        \tU_converged     \tk_solver        \tk_initial       \tk_final         \tk_iters         \tk_converged     \tomega_solver    \tomega_initial   \tomega_final     \tomega_iters     \tomega_converged \tp_solver        \tp_initial       \tp_final         \tp_iters         \tp_converged
1               \tsmoothSolver\t9.99968720e-01\t7.35557910e-02\t8\t1.00000000e+00\t2.11952440e-02\t8\t9.99985510e-01\t6.76813280e-02\t8\tfalse\tsmoothSolver\t1.00000000e+00\t1.05555430e-05\t2\ttrue\tsmoothSolver\t9.99948250e-01\t1.02768840e-06\t2\ttrue\tGAMG\t1.00000000e+00\t2.56365880e-02\t5\ttrue
80              \tsmoothSolver\t3.03155130e-04\t1.01519460e-05\t4\t9.85614300e-04\t1.56161330e-05\t4\t3.69729730e-04\t5.05097980e-06\t4\tfalse\tsmoothSolver\t1.78959720e-08\t1.35157230e-10\t2\ttrue\tsmoothSolver\t9.94567800e-09\t9.94567800e-09\t0\ttrue\tGAMG\t1.44469900e-05\t2.17740360e-07\t2\ttrue
`;
const FORCE = `# Force         \n# CofR          : (3.72312500e-01 0.00000000e+00 9.64136840e-02)\n#\n# Time          \ttotal_x total_y total_z\tpressure_x pressure_y pressure_z\tviscous_x viscous_y viscous_z\n80               7.29443464e+02 -1.35461339e-01 5.09960217e+03 5.74673156e+02 -1.22251279e-01 5.09306761e+03 1.54770308e+02 -1.32100607e-02 6.53456055e+00\n`;
const YPLUS = `# y+ ()         \n# Time          \tpatch           \tmin             \tmax             \taverage         \n27              \tbody\t5.60658344e+02\t9.67984722e+03\t4.93542468e+03\n54              \tbody\t5.58069784e+02\t9.67689010e+03\t4.93738436e+03\n`;
const FOAMLOG = `Starting time loop\n\nTime = 1\n\nsmoothSolver:  Solving for Ux, Initial residual = 0.99996872, Final residual = 0.073555791, No Iterations 8\nGAMG:  Solving for p, Initial residual = 1, Final residual = 0.025636588, No Iterations 5\nGAMG:  Solving for p, Initial residual = 0.2, Final residual = 0.01, No Iterations 3\nExecutionTime = 1.2 s  ClockTime = 3 s\n\nTime = 2\n\nsmoothSolver:  Solving for Ux, Initial residual = 0.10496277, Final residual = 0.0084006783, No Iterations 4\nGAMG:  Solving for p, Initial residual = 0.029827393, Final residual = 0.0011220367, No Iterations 3\nEnd\n`;
const CHECK = `Mesh stats \n    points:           53048\n    cells:            43207\n\nChecking geometry...\n    Mesh non-orthogonality Max: 59.784146 average: 9.6348559\n ***Max skewness = 6.0301124, 15 highly skew faces detected which may impair the quality of the results\nFailed 1 mesh checks.\n`;
{
  const r = R.parseFoamCoeffs(COEFF);
  near('OpenFOAM coefficients: Cl', r.values.CL, 1.68487364e-01, 1e-12); near('OpenFOAM coefficients: Cd', r.values.CD, 3.01345443e-02, 1e-12); near('OpenFOAM coefficients: pitching moment', r.values.Cm, -4.39990801e-02, 1e-12);
  ok('OpenFOAM coefficients: charts and table', r.plots.length === 2 && r.tables[0].rows.length === 12 && r.tables[0].rows[0][0] === 'Cd');
  const old = R.parseFoamCoeffs('# Time        Cm            Cd            Cl            Cl(f)         Cl(r)\n1 0.01 0.03 0.4 0.2 0.2\n2 0.011 0.031 0.41 0.2 0.21\n');
  near('OpenFOAM coefficients: older forceCoeffs.dat column order', old.values.CL, 0.41, 1e-12); near('…and its moment column', old.values.Cm, 0.011, 1e-12);
  const si = R.parseFoamSolverInfo(SINFO); ok('OpenFOAM solver information: residual chart on a log axis and final pressure residual', si.plots[0].ylog === true && si.plots[0].series.map((s) => s.name).join() === 'Ux,Uy,Uz,k,omega,p' && Math.abs(si.values.p_residual - 1.444699e-5) < 1e-12);
  const fo = R.parseFoamForces(FORCE); near('OpenFOAM forces: total z force', fo.values.force_z_N, 5099.60217, 1e-9);
  const yp = R.parseFoamYPlus(YPLUS); near('OpenFOAM y+: maximum of the last write', yp.values.yplus_max, 9676.8901, 1e-9); ok('OpenFOAM y+: outside the wall-function range is flagged', yp.kpis[0].status === 'warn');
  const lg = R.parseFoamLog(FOAMLOG); ok('OpenFOAM log: first residual of each equation per step', lg.values.iterations === 2 && lg.plots[0].series.find((s) => s.name === 'p').y.join() === '1,0.029827393');
  const cm = R.parseFoamCheckMesh(CHECK); ok('OpenFOAM checkMesh: cell count and failed check', cm.values.mesh_cells === 43207 && /1 failed/.test(cm.warnings[0]));
}
// --- CalculiX 2.23 (.dat trimmed from real runs of the generated wing box; .frd built here in the real fixed-width layout)
const DAT_STATIC = `\n                        S T E P       1\n\n\n                                INCREMENT     1\n\n\n displacements (vx,vy,vz) for set TIP and time  0.1000000E+01\n\n       321 -6.649267E-04  6.347190E-03  4.216487E-01\n       322 -6.657403E-04  6.342864E-03  4.226339E-01\n       336 -4.071725E-11  5.380887E-12  4.167110E-01\n\n total force (fx,fy,fz) for set ROOT and time  0.1000000E+01\n\n       -1.294059E-07 -1.035109E-04 -2.006634E+04\n`;
const DAT_MODAL = `\n                        S T E P       1\n\n\n     E I G E N V A L U E   O U T P U T\n\n MODE NO    EIGENVALUE                       FREQUENCY   \n                                     REAL PART            IMAGINARY PART\n                           (RAD/TIME)      (CYCLES/TIME     (RAD/TIME)\n\n      1   0.2238654E+04   0.4731441E+02   0.7530323E+01   0.0000000E+00\n      2   0.9446279E+04   0.9719197E+02   0.1546858E+02   0.0000000E+00\n      3   0.1157100E+05   0.1075686E+03   0.1712007E+02   0.0000000E+00\n\n     P A R T I C I P A T I O N   F A C T O R S\n\nMODE NO.   X-COMPONENT     Y-COMPONENT     Z-COMPONENT     X-ROTATION      Y-ROTATION      Z-ROTATION\n\n      1  -0.2104221E-07   0.7245494E-09  -0.4169990E+01  -0.1645744E+02   0.2522212E+01   0.5941823E-08\n      2  -0.5029533E-04   0.2946415E-02  -0.1525297E-05  -0.6650130E-06   0.9698065E-06   0.1870866E-02\n`;
const DAT_BUCKLE = `\n                        S T E P       1\n\n\n     B U C K L I N G   F A C T O R   O U T P U T\n\n MODE NO       BUCKLING\n                FACTOR\n\n      1   0.9373322E+00\n      2   0.9931912E+00\n      3   0.1019650E+01\n      4   0.1098173E+01\n`;
const DAT_DYN = [0.001, 0.002, 0.003, 0.004].map((t, i) => `\n displacements (vx,vy,vz) for set TIPREF and time  ${t.toExponential(7).toUpperCase()}\n\n        51  2.310198E-11  2.314787E-12  ${(0.01 * Math.sin(i + 1)).toExponential(6).toUpperCase()}\n`).join('');
const e12 = (v) => { const s = v.toExponential(5).toUpperCase().replace(/E([+-])(\d)$/, 'E$10$2'); return s.padStart(12); };
function frd(nodes, blocks) {
  const L = ['    1C', '    1UUSER', `    2C${String(nodes.length).padStart(30)}${' '.repeat(37)}1`];
  nodes.forEach((p, i) => L.push(` -1${String(i + 1).padStart(10)}${p.map(e12).join('')}`)); L.push(' -3');
  blocks.forEach((b, k) => {
    L.push(`    1PSTEP${String(k + 1).padStart(26)}${String(1).padStart(12)}${String(1).padStart(12)}`);
    if (b.mode) L.push(`    1PMODE${String(b.mode).padStart(26)}`);
    L.push(`  100CL  ${100 + k + 1} ${String(b.value)}${String(nodes.length).padStart(12)}                     ${b.ictype}    ${k + 1}${b.ictype === 2 ? 'MODAL' : '     '}      1`);
    L.push(` -4  ${b.name.padEnd(8)}${String(b.name === 'STRESS' ? 6 : 4).padStart(5)}    1`, ...Array.from({ length: b.name === 'STRESS' ? 6 : 4 }, () => ' -5  X           1    2    1    0'));
    b.data.forEach((row, i) => L.push(` -1${String(i + 1).padStart(10)}${row.map(e12).join('')}`)); L.push(' -3');
  });
  return L.join('\n') + '\n 9999\n';
}
{
  const s = R.parseCcxDat(DAT_STATIC); near('CalculiX .dat: mean tip deflection', s.values.tip_deflection_m, (0.4216487 + 0.4226339 + 0.416711) / 3, 1e-9); near('CalculiX .dat: root reaction', s.values.root_reaction_z_N, -20066.34, 1e-9);
  const m = R.parseCcxDat(DAT_MODAL); ok('CalculiX .dat: eigenfrequencies only (participation factors not mistaken for modes)', m.values.frequencies_Hz.join() === '7.530323,15.46858,17.12007' && m.tables[0].rows.length === 3 && m.plots[0].type === 'bar');
  const b = R.parseCcxDat(DAT_BUCKLE); near('CalculiX .dat: first buckling factor', b.values.buckling_load_factor, 0.9373322, 1e-12); ok('CalculiX .dat: a factor below one is "not acceptable"', b.kpis[0].status === 'bad' && b.tables[0].rows.length === 4);
  const d = R.parseCcxDat(DAT_DYN); ok('CalculiX .dat: time history of the tip node', d.plots[0].series[0].x.join() === '0.001,0.002,0.003,0.004' && Math.abs(d.values.tip_deflection_dynamic_peak_m - 0.01 * Math.sin(2)) < 1e-8);
  // synthetic cantilever strip: 21 stations × 3 chordwise nodes
  const nodes = [], NS = 21; for (let k = 0; k < NS; k++) for (const x of [0, 0.5, 1]) nodes.push([x, (5 * k) / (NS - 1), 0]);
  const eta = (p) => p[1] / 5, bend1 = nodes.map((p) => [0, 0, eta(p) ** 2]), bend2 = nodes.map((p) => [0, 0, eta(p) ** 2 * (1 - 1.4 * eta(p))]), tors = nodes.map((p) => [0, 0, (p[0] - 0.5) * eta(p)]), lag = nodes.map((p) => [eta(p) ** 2, 0, 0]), local = nodes.map((p) => [0, 0, (p[0] === 0.5 ? 1 : -0.5) * Math.sin(8 * Math.PI * eta(p))]);
  const stat = nodes.map((p) => [1e-4, 2e-4, 0.4 * eta(p) ** 2]), stress = nodes.map((p) => [1e8 * (1 - eta(p)), 0, 0, 0, 0, 0]);
  const fs1 = R.parseCcxFrd(frd(nodes, [{ name: 'DISP', ictype: 0, value: 1, data: stat }, { name: 'STRESS', ictype: 0, value: 1, data: stress }]));
  near('CalculiX .frd: largest vertical displacement', fs1.values.tip_deflection_max_m, 0.4, 1e-5); near('CalculiX .frd: peak von Mises stress', fs1.values.sigma_max_Pa, 1e8, 1e-5);
  ok('CalculiX .frd: deflected shape along the span', fs1.plots[0].series[0].x.length === 40 && fs1.plots[0].xlabel.includes('y'));
  const fm = R.parseCcxFrd(frd(nodes, [{ name: 'DISP', ictype: 2, mode: 1, value: 7.53, data: bend1 }, { name: 'DISP', ictype: 2, mode: 2, value: 15.4, data: local }, { name: 'DISP', ictype: 2, mode: 3, value: 21.5, data: lag }, { name: 'DISP', ictype: 2, mode: 4, value: 38.2, data: tors }, { name: 'DISP', ictype: 2, mode: 5, value: 46.9, data: bend2 }]));
  ok('CalculiX .frd: modes classified from their shapes', fm.values.mode_kinds.join() === 'bending,local,inplane,torsion,bending', fm.values.mode_kinds.join());
  ok('CalculiX .frd: bending and torsion frequencies skip local and in-plane modes', fm.values.f1_Hz === 7.53 && fm.values.f2_Hz === 46.9 && fm.values.f3_Hz === undefined && fm.values.f_torsion_Hz === 38.2);
  const fb = R.parseCcxFrd(frd(nodes, [{ name: 'DISP', ictype: 4, value: 0, data: stat }, { name: 'DISP', ictype: 4, value: 0.937, data: local }])); near('CalculiX .frd: buckling factor (reference state skipped)', fb.values.buckling_load_factor, 0.937, 1e-9);
  throws('CalculiX .frd: wrong file refused', () => R.parseCcxFrd('hello\n'), /frd/);

  // dispatcher, merge and comparison
  const c = PRESETS.ga.data(), built = D.buildCase(c, { solver: 'calculix', settings: { n_span: 8 } });
  const all = R.readResults([{ name: 'out/static.dat', bytes: enc(DAT_STATIC) }, { name: 'out/modal.dat', bytes: enc(DAT_MODAL) }, { name: 'out/buckle.dat', bytes: enc(DAT_BUCKLE) }, { name: 'out/static.frd', bytes: enc(frd(nodes, [{ name: 'DISP', ictype: 0, value: 1, data: stat }, { name: 'STRESS', ictype: 0, value: 1, data: stress }])) }, { name: 'out/modal.frd', bytes: enc(frd(nodes, [{ name: 'DISP', ictype: 2, mode: 1, value: 7.53, data: bend1 }, { name: 'DISP', ictype: 2, mode: 2, value: 15.4, data: local }])) }, { name: 'out/manifest.json', bytes: enc(file(built, 'manifest.json')) }, { name: 'out/notes.txt', bytes: enc('hi') }]);
  ok('dispatcher: solver, manifest and every file accounted for', all.solver === 'calculix' && all.manifest.caseHash === D.caseHash(c) && all.sources.length === 7 && all.sources.filter((s) => s.ok).length === 6 && all.sources.find((s) => s.name === 'out/notes.txt').ok === false);
  ok('dispatcher: shape-based frequencies replace solver order (the local mode is not "second bending")', all.values.f1_Hz === 7.53 && all.values.f2_Hz === undefined && !all.kpis.some((k) => /^mode\d_Hz$/.test(k.key)));
  ok('dispatcher: result is plain data', JSON.stringify(structuredClone({ k: all.kpis, p: all.plots, t: all.tables })) === JSON.stringify({ k: all.kpis, p: all.plots, t: all.tables }));
  const up = { fea: { tip_deflection_m: 0.4, sigma_max_Pa: 0.8e8, buckling_factor: 1.2 }, vibration: { f1_Hz: 7.0, f2_Hz: 44 } }, rows = R.compare(all.values, up, { manifest: all.manifest }), row = (k) => rows.find((x) => x.key === k);
  near('comparison: tip deflection difference', row('tip_deflection_m').pct, (100 * (all.values.tip_deflection_m - 0.4)) / 0.4, 1e-12); near('comparison: calibration ratio for the first frequency', row('f1_Hz').ratio, 7.53 / 7, 1e-12);
  ok('comparison: only quantities present on both sides get numbers; missing native values are explained', row('sigma_max_Pa').diff === 1e8 - 0.8e8 || Math.abs(row('sigma_max_Pa').diff - 0.2e8) < 1 && !row('f2_Hz') && R.compare({ CL: 0.5 }, {}).find((x) => x.key === 'CL').native === null && /Run the native/.test(R.compare({ CL: 0.5 }, {})[0].note));
  const su = R.readResults([{ name: 'history.csv', bytes: enc(HIST) }, { name: 'postProcessing/forceCoeffs1/0/coefficient.dat', bytes: enc(COEFF) }]);
  ok('dispatcher: files are recognised by content and path', su.sources.map((s) => s.kind).join() === 'su2-history,foam-coeffs' && su.values.CL === 0.438057899);
  const rec = R.hifiRecord(all, { ts: 5 }); ok('data-bus record: numbers plus provenance', rec.tip_deflection_m > 0.41 && rec.solver === 'calculix' && rec.analysis === 'all' && rec.caseHash === D.caseHash(c) && rec.ts === 5 && !('frequencies_Hz' in rec));
  const cp = R.comparePlot(rows); ok('comparison chart is a bar spec', cp.type === 'bar' && cp.categories.length === cp.series[0].y.length && cp.categories.length >= 3);
  const two = R.compare({ CL: 0.44, CD: 0.012 }, { cfd: { CL: 0.4, CD: 0.03 } }, { manifest: D.buildCase(c, { solver: 'su2', settings: { resolution: 'coarse' } }).manifest }); ok('comparison: a 2-D section result is labelled as such', /Section \(2-D\)/.test(two[0].note));
}

// ---------------------------------------------------------------------------------------------
sec('GitHub client (mock network)');
{
  const calls = [], state = { branch: false, file: null, runs: [], dispatched: null, polls: 0 }, zipArt = await Z.zipWrite([{ path: 'history.csv', text: HIST }, { path: 'manifest.json', text: '{}' }]);
  const json = (status, body) => ({ ok: status >= 200 && status < 300, status, json: async () => body, arrayBuffer: async () => new ArrayBuffer(0) });
  const fetchImpl = async (url, init) => {
    calls.push({ url, method: init.method, auth: init.headers.Authorization, body: init.body ? JSON.parse(init.body) : null });
    const p = url.replace('https://api.github.com/repos/jane/aero', '');
    if (p === '') return json(200, { default_branch: 'main', private: false, full_name: 'jane/aero', html_url: 'https://github.com/jane/aero' });
    if (p === '/actions/workflows/solve.yml') return json(200, { id: 7, state: 'active' });
    if (p === '/git/ref/heads/aerosuite-cases') return state.branch ? json(200, { object: { sha: 'b' } }) : json(404, { message: 'Not Found' });
    if (p === '/git/ref/heads/main') return json(200, { object: { sha: 'abc123' } });
    if (p === '/git/refs' && init.method === 'POST') { state.branch = true; return json(201, {}); }
    if (p.startsWith('/contents/cases/') && init.method === 'GET') return json(404, { message: 'Not Found' });
    if (p.startsWith('/contents/cases/') && init.method === 'PUT') { state.file = JSON.parse(init.body); return json(201, { content: { path: 'cases/x/case.zip', sha: 's' } }); }
    if (p === '/actions/workflows/solve.yml/dispatches') { state.dispatched = JSON.parse(init.body); state.runs.push({ id: 99, display_title: `Solve su2 · x ${state.dispatched.inputs.tag}`, status: 'queued', html_url: 'https://github.com/jane/aero/actions/runs/99' }); return { ok: true, status: 204, json: async () => { throw new Error('no body'); } }; }
    if (p.startsWith('/actions/workflows/solve.yml/runs')) return json(200, { workflow_runs: [{ id: 1, display_title: 'Solve su2 · other r000' }, ...state.runs] });
    if (p === '/actions/runs/99') { state.polls++; return json(200, { id: 99, status: state.polls >= 2 ? 'completed' : 'in_progress', conclusion: state.polls >= 2 ? 'success' : null, html_url: 'https://github.com/jane/aero/actions/runs/99' }); }
    if (p.startsWith('/actions/runs/99/jobs')) return json(200, { jobs: [{ steps: [{ name: 'Run the solver', status: 'in_progress' }] }] });
    if (p.startsWith('/actions/runs/99/artifacts')) return json(200, { artifacts: [{ id: 5, name: 'results-su2-x', size_in_bytes: zipArt.length, expired: false }, { id: 6, name: 'old', size_in_bytes: 1, expired: true }] });
    if (p === '/actions/artifacts/5/zip') return { ok: true, status: 200, arrayBuffer: async () => zipArt.buffer.slice(zipArt.byteOffset, zipArt.byteOffset + zipArt.byteLength), json: async () => ({}) };
    return json(500, { message: `unexpected ${init.method} ${p}` });
  };
  ok('repository parsing accepts owner/name, URLs and git remotes', ['jane/aero', 'https://github.com/jane/aero', 'https://github.com/jane/aero.git', 'git@github.com:jane/aero.git', ' jane/aero/ '].every((s) => { const r = GH.parseRepo(s); return r.owner === 'jane' && r.repo === 'aero'; }));
  throws('repository parsing refuses anything else', () => GH.parseRepo('https://evil.example/jane/aero'), /owner\/name/); throws('…including path tricks', () => GH.parseRepo('jane/aero/../../x'), /owner\/name/);
  throws('a missing token is refused before any request', () => GH.createGithub({ token: '', owner: 'jane', repo: 'aero', fetchImpl }), /token/);
  const TOKEN = 'github_pat_' + 'A1b2C3d4E5'.repeat(8), gh = GH.createGithub({ token: TOKEN, owner: 'jane', repo: 'aero', fetchImpl }), noSleep = async () => {};
  // small case → inline
  const small = await Z.zipWrite([{ path: 'run.sh', text: 'echo hi\n' }]), status = [];
  const sub = await GH.submitCase(gh, { zipBytes: small, name: 'GA case/su2 #1', solver: 'su2', threads: 2, onStatus: (m) => status.push(m) });
  ok('small case travels inline, dispatched on the default branch', sub.route === 'inline' && state.dispatched.ref === 'main' && state.dispatched.inputs.case_b64 === Z.base64Encode(small) && state.dispatched.inputs.solver === 'su2' && state.dispatched.inputs.threads === '2' && !('branch' in state.dispatched.inputs) && state.file === null);
  ok('case name is made safe for the workflow', /^[A-Za-z0-9._-]+$/.test(state.dispatched.inputs.case) && state.dispatched.inputs.case === 'GA-case-su2-1');
  const run = await GH.followRun(gh, sub, { onStatus: (m) => status.push(m), sleepFn: noSleep });
  ok('run is found by its tag and followed to completion', run.id === 99 && run.conclusion === 'success' && status.some((m) => /Run the solver/.test(m)) && status.some((m) => /Queued/.test(m)));
  const arts = await GH.fetchArtifacts(gh, 99), inner = await Z.zipRead(arts[0].bytes);
  ok('artifact archive is downloaded (expired ones skipped) and unpacked', arts.length === 1 && inner.map((f) => f.name).join() === 'history.csv,manifest.json' && R.readResults(inner).values.CL === 0.438057899);
  // large case → branch upload
  state.dispatched = null; const bigCase = await Z.zipWrite([{ path: 'mesh.su2', bytes: Uint8Array.from({ length: 90000 }, (_, i) => (i * 2654435761) >>> 24) }], { deflate: false });
  const sub2 = await GH.submitCase(gh, { zipBytes: bigCase, name: 'big', solver: 'calculix' });
  ok('large case is committed to a separate branch, never the default one', sub2.route === 'branch' && sub2.branch === 'aerosuite-cases' && state.branch && state.file.branch === 'aerosuite-cases' && state.file.content === Z.base64Encode(bigCase) && state.dispatched.inputs.branch === 'aerosuite-cases' && !('case_b64' in state.dispatched.inputs) && state.dispatched.ref === 'main');
  ok('upload path and commit message carry no attribution', calls.some((c) => c.method === 'PUT' && c.url.endsWith('/contents/cases/big/case.zip')) && state.file.message === 'Add solver case big' && !('committer' in state.file) && !('author' in state.file));
  ok('every request went to api.github.com with the token; nothing else was contacted', calls.length > 10 && calls.every((c) => c.url.startsWith('https://api.github.com/repos/jane/aero') && c.auth === `Bearer ${TOKEN}`));
  let refused = ''; try { await GH.createGithub({ token: TOKEN, owner: 'jane', repo: 'aero', fetchImpl: async (u) => { refused = u; return json(200, {}); } }).downloadArtifact('1/../../../../evil'); } catch { /* fine either way */ } ok('artifact ids cannot redirect the token elsewhere', refused === '' || refused.startsWith('https://api.github.com/'));
  const deny = GH.createGithub({ token: TOKEN, owner: 'jane', repo: 'aero', fetchImpl: async () => json(403, { message: 'Resource not accessible by personal access token' }) }); let msg = ''; try { await deny.dispatch('main', {}); } catch (e) { msg = e.message; } ok('a permission error names the permission to grant', /Actions \(read and write\)/.test(msg) && /Resource not accessible/.test(msg));
  const offline = GH.createGithub({ token: TOKEN, owner: 'jane', repo: 'aero', fetchImpl: async () => { throw new TypeError('Failed to fetch'); } }); msg = ''; try { await offline.repoInfo(); } catch (e) { msg = e.message; } ok('a network failure is explained', /did not go through/.test(msg));
}

// ---------------------------------------------------------------------------------------------
sec('workflow file');
{
  const y = fs.readFileSync(new URL('../.github/workflows/solve.yml', import.meta.url), 'utf8');
  ok('manually dispatched only, with the documented inputs', /^on:\n  workflow_dispatch:/m.test(y) && !/^\s+(push|pull_request|schedule):/m.test(y) && ['solver:', 'case:', 'branch:', 'case_b64:', 'tag:', 'threads:'].every((k) => y.includes(`      ${k}`)) && /options: \[su2, calculix, openfoam\]/.test(y));
  ok('read-only: no write permission, no commit, no push', /^permissions:\n  contents: read$/m.test(y) && !/contents: write|git (push|commit)|GITHUB_TOKEN|secrets\./.test(y) && /persist-credentials: false/.test(y));
  ok('solvers come from their official sources', y.includes('https://github.com/su2code/SU2/releases/download/') && y.includes('calculix-ccx') && y.includes('https://dl.openfoam.com/add-debian-repo.sh'));
  ok('inputs reach the shell only through environment variables', !/run: \|[\s\S]*?\$\{\{ inputs\./.test(y.split('steps:')[1].replace(/with:[\s\S]*?(?=\n      - name:|$)/g, '').replace(/if: \$\{\{[^\n]*\n/g, '')));
  ok('results are uploaded as an artifact even when the solver fails', /actions\/upload-artifact@v4/.test(y) && /if: \$\{\{ always\(\) \}\}/.test(y));
  ok('the run name carries the tag the app searches for', /^run-name: .*\$\{\{ inputs\.tag \}\}/m.test(y));
  const tabs = y.split('\n').filter((l) => /\t/.test(l)); ok('no tab characters (YAML)', tabs.length === 0);
}

// ---------------------------------------------------------------------------------------------
// optional: the real solvers
// ---------------------------------------------------------------------------------------------
const scratch = process.env.BRIDGE_SCRATCH || fs.mkdtempSync(path.join(os.tmpdir(), 'bridge-real-'));
const writeCase = (b, dir) => { fs.rmSync(dir, { recursive: true, force: true }); for (const f of b.files) { const p = path.join(dir, f.path); fs.mkdirSync(path.dirname(p), { recursive: true }); fs.writeFileSync(p, f.text); } };
const slurp = (dir, names) => names.filter((n) => fs.existsSync(path.join(dir, n))).map((n) => ({ name: n, bytes: new Uint8Array(fs.readFileSync(path.join(dir, n))) }));
const real = [];
if (process.env.BRIDGE_SU2) {
  sec('real SU2');
  const dir = path.join(scratch, 'su2-euler'), b = D.buildCase(PRESETS.ga.data(), { solver: 'su2', settings: { analysis: 'euler', regime: 'compressible', mach: 0.5, mesh: 'O-grid', resolution: 'coarse', iter: 250 } }); writeCase(b, dir);
  const r = spawnSync(process.env.BRIDGE_SU2, ['config.cfg'], { cwd: dir, encoding: 'utf8', timeout: 1500e3, maxBuffer: 64e6, env: { ...process.env, OMP_NUM_THREADS: '2', OMP_WAIT_POLICY: 'passive' } });
  ok('SU2_CFD accepts the generated configuration and mesh and exits successfully', r.status === 0 && /Exit Success \(SU2_CFD\)/.test(r.stdout), (r.stdout || '').slice(-600) + (r.stderr || '').slice(-300));
  if (r.status === 0) { const res = R.readResults(slurp(dir, ['history.csv', 'surface_flow.vtu', 'flow.vtu', 'manifest.json'])); ok('real history.csv, surface_flow.vtu and flow.vtu are parsed', res.values.iterations === 250 && Number.isFinite(res.values.CL) && Number.isFinite(res.values.cp_min) && res.plots.some((p) => p.type === 'tri')); ok('inviscid NACA 2412 at 2° and Mach 0.5: lift near thin-aerofoil theory with compressibility, drag near zero', res.values.CL > 0.4 && res.values.CL < 0.62 && Math.abs(res.values.CD) < 0.01, `CL ${res.values.CL}, CD ${res.values.CD}`); real.push(`SU2 Euler O-grid: ${res.values.iterations} iterations, CL ${res.values.CL.toFixed(4)}, CD ${res.values.CD.toFixed(5)}`); }
}
if (process.env.BRIDGE_CCX) {
  sec('real CalculiX');
  const dir = path.join(scratch, 'ccx-box'), c = PRESETS.ga.data(), b = D.buildCase(c, { solver: 'calculix', settings: { analysis: 'all', n_span: 16, n_chord: 4 } }); writeCase(b, dir);
  for (const job of ['static', 'buckle', 'modal']) { const r = spawnSync(process.env.BRIDGE_CCX, ['-i', job], { cwd: dir, encoding: 'utf8', timeout: 1500e3, maxBuffer: 64e6, env: { ...process.env, OMP_NUM_THREADS: '2' } }); ok(`ccx runs ${job}.inp to completion`, r.status === 0 && /Job finished/.test(r.stdout) && !/\*ERROR/.test(r.stdout), (r.stdout || '').slice(-500)); }
  const res = R.readResults(slurp(dir, ['static.dat', 'static.frd', 'buckle.dat', 'buckle.frd', 'modal.dat', 'modal.frd', 'manifest.json'])), v = res.values;
  near('root reaction balances the applied load', -v.root_reaction_z_N, b.meshInfo.appliedLoad_N, 1e-4);
  ok('tip deflects upward by a plausible fraction of the span', v.tip_deflection_m > 0 && v.tip_deflection_m / b.meshInfo.half > 0.005 && v.tip_deflection_m / b.meshInfo.half < 0.3, `${v.tip_deflection_m}`);
  ok('first mode is vertical bending; stress, buckling factor and frequency are finite', v.mode_kinds?.[0] === 'bending' && v.f1_Hz > 0.5 && v.f1_Hz < 100 && v.sigma_max_Pa > 1e6 && v.buckling_load_factor > 0, JSON.stringify({ f1: v.f1_Hz, s: v.sigma_max_Pa, b: v.buckling_load_factor, k: v.mode_kinds }));
  real.push(`CalculiX wing box: tip ${v.tip_deflection_m.toFixed(4)} m, f1 ${v.f1_Hz.toFixed(3)} Hz, buckling factor ${v.buckling_load_factor.toFixed(3)}, peak stress ${(v.sigma_max_Pa / 1e6).toFixed(1)} MPa`);
}
if (process.env.BRIDGE_GMSH) {
  sec('real Gmsh');
  const dir = path.join(scratch, 'gmsh'); fs.mkdirSync(dir, { recursive: true }); fs.writeFileSync(path.join(dir, 'aircraft.geo'), D.gmshGeo(PRESETS.ga.data()) + '\nhWall = 0.25;\n');
  const r = spawnSync(process.env.BRIDGE_GMSH, ['aircraft.geo', '-2', '-format', 'su2', '-o', 'surface.su2', '-v', '3'], { cwd: dir, encoding: 'utf8', timeout: 900e3, maxBuffer: 64e6 }), out = fs.existsSync(path.join(dir, 'surface.su2')) ? fs.readFileSync(path.join(dir, 'surface.su2'), 'utf8') : '';
  ok('gmsh builds the geometry, names the boundaries and meshes the surfaces', r.status === 0 && !/Error/.test(r.stdout + r.stderr) && ['wall', 'symmetry', 'farfield'].every((m) => new RegExp(`MARKER_TAG= ${m}\\b`).test(out)), (r.stdout + r.stderr).slice(-600));
  if (r.status === 0) real.push('Gmsh: parametric half-model surface-meshed with wall / symmetry / farfield groups');
}
if (process.env.BRIDGE_FOAM_BASHRC) {
  sec('real OpenFOAM');
  const dir = path.join(scratch, 'foam'), b = D.buildCase(PRESETS.fixedUav.data(), { solver: 'openfoam', settings: { analysis: 'simpleFoam', cells_per_chord: 10, layers: 0, iterations: 20 } }); writeCase(b, dir);
  const r = spawnSync('bash', ['-c', `set +u; source "${process.env.BRIDGE_FOAM_BASHRC}" >/dev/null 2>&1; bash run.sh`], { cwd: dir, encoding: 'utf8', timeout: 2400e3, maxBuffer: 256e6 });
  const files = []; const walk = (d) => { for (const e of fs.readdirSync(d, { withFileTypes: true })) { const p = path.join(d, e.name); if (e.isDirectory()) { if (d !== dir || e.name === 'postProcessing') walk(p); } else if (/^log\.|\.dat$|manifest\.json/.test(e.name)) files.push({ name: path.relative(dir, p), bytes: new Uint8Array(fs.readFileSync(p)) }); } }; walk(dir);
  const res = R.readResults(files);
  ok('blockMesh, snappyHexMesh and simpleFoam run the generated case', r.status === 0 && res.values.iterations === 20 && Number.isFinite(res.values.CL) && res.values.mesh_cells > 1000, (r.stdout || '').slice(-500));
  if (r.status === 0) real.push(`OpenFOAM simpleFoam: ${res.values.mesh_cells} cells, ${res.values.iterations} iterations, CL ${res.values.CL.toFixed(4)}, CD ${res.values.CD.toFixed(4)}`);
}
if (!process.env.BRIDGE_SCRATCH) fs.rmSync(scratch, { recursive: true, force: true });

// ---------- summary ----------
for (const r of real) console.log('  real solver: ' + r);
console.log(`bridge: ${pass} passed, ${fail} failed${real.length ? '' : ' (real-solver checks skipped: set BRIDGE_SU2 / BRIDGE_CCX / BRIDGE_GMSH / BRIDGE_FOAM_BASHRC to run them)'}`);
if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }

