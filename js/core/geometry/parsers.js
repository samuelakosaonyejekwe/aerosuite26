// importFile(): detect the format from the content, run the matching reader, validate what it produced and
// return a Model with a provenance log. Files are untrusted: a recognised file that cannot be read never
// throws - it comes back as a metadata-only Model whose warnings say why and what conversion route to use.
// Nothing is repaired on import and units are never guessed.

import { FORMATS, detectFormat, formatById } from './formats.js';
import { ET, displayTriangles, displayLines } from './parsers-util.js';
import * as S from './parsers-surface.js';
import * as Me from './parsers-mesh.js';
import * as C from './parsers-cad.js';
import * as K from './parsers-kernel.js';
import * as H from './parsers-hdf.js';
import * as R from './parsers-raster.js';

const READERS = {
  stl: S.readSTL, obj: S.readOBJ, ply: S.readPLY, off: S.readOFF, gltf: S.readGLTF, '3mf': S.read3MF, amf: S.readAMF, json: S.readJSON,
  dae: S.readDAE, vrml: S.readVRML, x3d: S.readX3D, xyz: S.readXYZ, pcd: S.readPCD, las: S.readLAS,
  vtk: Me.readVTK, vtkxml: Me.readVTKXML, gmsh: Me.readGmsh, su2: Me.readSU2, nastran: Me.readNastran, abaqus: Me.readAbaqus, unv: Me.readUNV,
  cdb: Me.readCDB, fluent: Me.readFluent, tecplot: Me.readTecplot, plot3d: Me.readPlot3D, openfoam: Me.readOpenFOAM,
  step: K.readSTEPKernel, iges: K.readIGESKernel, brep: K.readBREP, xt: C.readXT, laz: K.readLAZ,
  cgns: H.readHDF5Any, med: H.readHDF5Any, 'fluent-h5': H.readHDF5Any, exodus: H.readHDF5Any, hdf5: H.readHDF5Any, geotiff: R.readGeoTIFF, e57: R.readE57,
};
const STRUCTURAL = new Set(['nastran', 'abaqus', 'cdb', 'unv']);

/** Drop elements whose connectivity points outside the vertex list; returns the number removed. */
function validateElements(elements, nVerts) {
  let dropped = 0;
  for (let b = elements.length - 1; b >= 0; b--) {
    const el = elements[b], info = ET[el.type];
    if (!info || el.nodesPer !== info.n || !(el.conn instanceof Uint32Array) || el.conn.length !== el.count * el.nodesPer) { dropped += el.count || 0; elements.splice(b, 1); continue; }
    let bad = 0; const np = el.nodesPer, c = el.conn;
    for (let e = 0; e < el.count; e++) for (let k = 0; k < np; k++) if (c[e * np + k] >= nVerts) { bad++; break; }
    if (!bad) continue;
    const conn = new Uint32Array((el.count - bad) * np), grp = el.group ? new Int32Array(el.count - bad) : null; let w = 0;
    for (let e = 0; e < el.count; e++) { let ok = true; for (let k = 0; k < np; k++) if (c[e * np + k] >= nVerts) { ok = false; break; } if (!ok) continue; for (let k = 0; k < np; k++) conn[w * np + k] = c[e * np + k]; if (grp) grp[w] = el.group[e]; w++; }
    el.conn = conn; el.count = w; if (grp) el.group = grp; dropped += bad;
    if (!w) elements.splice(b, 1);
  }
  return dropped;
}
function validIndices(arr, per, nVerts) {
  let bad = 0; const n = Math.floor(arr.length / per);
  for (let e = 0; e < n; e++) for (let k = 0; k < per; k++) if (arr[e * per + k] >= nVerts) { bad++; break; }
  if (!bad && arr.length === n * per) return { arr, bad: 0 };
  const out = new Uint32Array((n - bad) * per); let w = 0;
  for (let e = 0; e < n; e++) { let ok = true; for (let k = 0; k < per; k++) if (arr[e * per + k] >= nVerts) { ok = false; break; } if (ok) { for (let k = 0; k < per; k++) out[w * per + k] = arr[e * per + k]; w++; } }
  return { arr: out, bad };
}

/**
 * Import one file. opts:
 *   format      force a format id from FORMATS instead of detecting it
 *   companions  [{ name, bytes }] extra files (OpenFOAM points/faces/boundary/owner, external glTF .bin buffers)
 *   maxPoints   cap for subsampled point clouds (LAS / LAZ), default 200 000
 *   maxGrid     largest number of height-field samples per side for GeoTIFF rasters, default 300
 *   linearDeflection, angularDeflection, linearDeflectionType   tessellation accuracy for STEP / IGES / BREP
 *                (default 0.001 of the bounding box and 0.5 rad; 'absolute_value' makes linearDeflection a length)
 *   kernelTimeout  milliseconds after which a STEP / IGES / BREP translation is stopped (default 180 000)
 *   wasm        false to skip the WebAssembly kernels (text-level results only)
 * Throws only for empty input or content that matches no known format.
 */
export async function importFile(fileName, bytes, opts = {}) {
  if (bytes instanceof ArrayBuffer) bytes = new Uint8Array(bytes);
  else if (ArrayBuffer.isView(bytes) && !(bytes instanceof Uint8Array)) bytes = new Uint8Array(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  if (!(bytes instanceof Uint8Array) || bytes.length === 0) throw new Error(`"${fileName}" is empty or was not supplied as bytes.`);
  const det = opts.format && formatById(opts.format) ? { id: opts.format, confidence: 1, note: 'format chosen by the user' } : detectFormat(fileName, bytes);
  const fmt = formatById(det.id);
  if (!fmt) throw new Error(`"${fileName}" was not recognised as any supported geometry or mesh format (${det.note}). Supported formats: ${FORMATS.map((f) => f.ext[0]).join(', ')}.`);
  const model = {
    name: String(fileName || 'model').split(/[\\/]/).pop(), format: fmt.id, formatName: fmt.name, support: fmt.support, kind: 'metadata-only',
    positions: new Float64Array(0), triangles: new Uint32Array(0), lines: new Uint32Array(0), elements: [], groups: [],
    units: { length: null, source: null }, meta: {}, log: [], warnings: [],
  };
  const log = (step, detail) => model.log.push({ step, detail });
  log('source', `${model.name}, ${bytes.length} bytes`);
  log('detect', `${fmt.name} [${fmt.support}] — confidence ${det.confidence.toFixed(2)}: ${det.note}`);
  if (det.confidence < 0.5) model.warnings.push(`Format identified with low confidence (${det.note}).`);
  const ctx = { name: model.name, opts, warn: (m) => { if (model.warnings.length < 200) model.warnings.push(m); } };
  const reader = READERS[fmt.id] || C.readMetadataOnly(fmt.id);
  let part = null;
  try { part = await reader(bytes, ctx); }
  catch (e) {
    // RangeError / TypeError from a malformed file are reported the same way as the readers' own errors
    const why = e instanceof RangeError ? 'the data ends before the header says it should, or an offset points outside the file' : e?.message || String(e);
    model.warnings.push(`${fmt.name}: the file could not be read — ${why}.`);
    if (e && e.meta && typeof e.meta === 'object') model.meta = e.meta;
    log('parse', `reader failed: ${why}`);
  }
  if (part) {
    if (part.formatId && part.formatId !== fmt.id && formatById(part.formatId)) { const f2 = formatById(part.formatId); model.format = f2.id; model.formatName = f2.name; model.support = f2.support; log('schema', `content identified as ${f2.name}`); }
    if (part.support) model.support = part.support;
    if (part.positions instanceof Float64Array && part.positions.length % 3 === 0) model.positions = part.positions;
    const nv = model.positions.length / 3;
    model.meta = part.meta || {};
    if (part.units && part.units.length) model.units = { length: part.units.length, source: part.units.source || 'file' };
    for (const w of part.warnings || []) ctx.warn(w);
    model.groups = (part.groups || []).map((g, i) => ({ ...g, id: i }));
    model.elements = (part.elements || []).filter((e) => e.count > 0);
    let dropped = validateElements(model.elements, nv);
    if (part.triangles) { const v = validIndices(part.triangles, 3, nv); model.triangles = v.arr; dropped += v.bad; }
    else model.triangles = displayTriangles(model.elements, nv);
    const ln = validIndices(part.lines && part.lines.length ? part.lines : displayLines(model.elements), 2, nv); model.lines = ln.arr; dropped += ln.bad;
    if (dropped) ctx.warn(`${dropped} element(s)/facet(s) referenced vertices beyond the vertex list and were dropped.`);
    let nonFinite = 0; const P = model.positions; for (let i = 0; i < P.length; i++) if (!Number.isFinite(P[i])) nonFinite++;
    if (nonFinite) ctx.warn(`${nonFinite} coordinate value(s) are not finite numbers (NaN/Inf); measurements of this model are unreliable until they are removed.`);
    const dimMax = model.elements.reduce((d, e) => Math.max(d, ET[e.type].dim), 0);
    model.kind = part.kind && part.kind !== 'metadata-only' ? part.kind
      : part.kind === 'metadata-only' ? 'metadata-only'
        : dimMax === 3 ? (STRUCTURAL.has(model.format) ? 'structural-mesh' : 'volume-mesh')
          : dimMax >= 1 ? (STRUCTURAL.has(fmt.id) ? 'structural-mesh' : fmt.category === 'mesh' ? 'surface-mesh' : 'surface')
            : model.triangles.length ? 'surface' : nv ? 'pointcloud' : 'metadata-only';
    if (model.kind !== 'metadata-only' && !nv) model.kind = 'metadata-only';
    const census = model.elements.map((e) => `${e.count} ${e.type}`).join(', ');
    log('parse', `${fmt.id} reader: ${nv} vertices, ${model.triangles.length / 3} display triangles, ${model.lines.length / 2} line segments${census ? `, elements: ${census}` : ''}, ${model.groups.length} group(s); kind = ${model.kind}`);
  }
  log('units', model.units.length ? `length unit stated by the file: ${model.units.length}` : 'the file does not state a length unit — to be confirmed by the user');
  if (model.kind === 'metadata-only') {
    const f3 = formatById(model.format) || fmt;
    model.warnings.push(`${f3.name} is ${f3.support === 'metadata' ? 'recognised but not decoded here' : 'readable here in general, but this file yielded no geometry'}: ${f3.limits}${f3.pathway ? ` Conversion pathway: ${f3.pathway}` : ''}`);
    log('pathway', f3.pathway || 'see the format capability profile');
  }
  return model;
}
