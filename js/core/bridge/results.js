// High-fidelity bridge — result readers. Parses what the external solvers write (SU2 history and surface
// files, OpenFOAM postProcessing tables and logs, CalculiX .dat and .frd) into the app's own
// { kpis, plots, tables } shapes, and sets the values beside the native suite outputs.
// Pure computation: no DOM, no network. Runs unchanged in the browser and in Node.

const dec = new TextDecoder();
const asText = (v) => (typeof v === 'string' ? v : dec.decode(v));
const last = (a) => a[a.length - 1];
const finite = (v) => typeof v === 'number' && Number.isFinite(v);
const base = (name) => String(name).replace(/\\/g, '/').split('/').pop();
/** Keep at most n points of a history (always including the last one). */
export function thin(arr, n = 400) {
  if (arr.length <= n) return arr.slice();
  const out = [], step = (arr.length - 1) / (n - 1);
  for (let i = 0; i < n; i++) out.push(arr[Math.round(i * step)]);
  return out;
}
const kpi = (key, label, value, unit = '-', status, note) => ({ key, label, value, unit, ...(status ? { status } : {}), ...(note ? { note } : {}) });
const empty = () => ({ kpis: [], plots: [], tables: [], values: {}, warnings: [] });

// ---------------------------------------------------------------------------------------------
// Delimited tables
// ---------------------------------------------------------------------------------------------
/** CSV with optional quoted header names; numeric cells become numbers. */
export function parseCsv(text) {
  const lines = asText(text).split(/\r?\n/).filter((l) => l.trim() !== '');
  if (!lines.length) return { columns: [], rows: [] };
  const columns = lines[0].split(',').map((s) => s.trim().replace(/^"|"$/g, '').trim()), rows = [];
  for (let i = 1; i < lines.length; i++) { const c = lines[i].split(','); if (c.length < 2) continue; rows.push(c.slice(0, columns.length).map((s) => { const t = s.trim(); if (t === '') return NaN; const v = Number(t); return Number.isNaN(v) && !/^-?nan$/i.test(t) ? t : v; })); }
  return { columns, rows };
}
const col = (t, name) => { const j = t.columns.indexOf(name); return j < 0 ? null : t.rows.map((r) => r[j]); };
const colLike = (t, re) => { const j = t.columns.findIndex((c) => re.test(c)); return j < 0 ? null : t.rows.map((r) => r[j]); };

// ---------------------------------------------------------------------------------------------
// SU2
// ---------------------------------------------------------------------------------------------
/** SU2 history.csv: convergence of residuals and force coefficients. */
export function parseSu2History(text) {
  const t = parseCsv(text), out = empty();
  if (!t.rows.length || !t.columns.includes('CL')) throw new Error('This does not look like an SU2 history file (no CL column).');
  const it = col(t, 'Inner_Iter') || t.rows.map((_, i) => i), CL = col(t, 'CL'), CD = col(t, 'CD'), n = t.rows.length;
  // nose-up positive: in 3-D that is SU2's CMy; in 2-D SU2's CMz is positive nose-down (x aft, y up), so its sign is flipped
  const threeD = t.columns.some((c) => /^rms\[(RhoW|W)\]$/.test(c)), cmRaw = col(t, threeD ? 'CMy' : 'CMz'), CM = cmRaw ? (threeD ? cmRaw : cmRaw.map((v) => -v)) : null;
  const res = t.columns.filter((c) => /^rms\[/.test(c)), lead = col(t, res[0] || '');
  const cl = last(CL), cd = last(CD), tail = Math.max(2, Math.round(0.1 * n)), clTail = CL.slice(-tail), drift = Math.max(...clTail) - Math.min(...clTail);
  const drop = lead ? lead[0] - last(lead) : NaN, bad = !finite(cl) || !finite(cd), settled = drift < 2e-3 * Math.max(0.05, Math.abs(cl));
  out.values = { CL: cl, CD: cd, LD: cd !== 0 ? cl / cd : NaN, iterations: n, residual_drop_decades: drop, CL_drift: drift };
  if (CM) out.values.Cm = last(CM);
  out.kpis.push(kpi('CL', 'Lift coefficient', cl, '-', bad ? 'bad' : settled ? 'ok' : 'warn', settled ? 'Steady over the last tenth of the run' : `Still moving: changed by ${drift.toPrecision(2)} over the last tenth of the run`),
    kpi('CD', 'Drag coefficient', cd), kpi('LD', 'Lift-to-drag ratio', out.values.LD));
  if (CM) out.kpis.push(kpi('Cm', 'Pitching-moment coefficient (about 25 % chord, nose-up positive)', last(CM)));
  out.kpis.push(kpi('iterations', 'Iterations run', n, ''), kpi('residual_drop_decades', 'Residual reduction', drop, 'decades', finite(drop) ? (drop >= 3 ? 'ok' : 'warn') : 'bad', 'At least 3 orders of magnitude is the usual minimum'));
  if (bad) out.warnings.push('The run ended with non-finite coefficients: the solver diverged. Lower the CFL number or check the mesh.');
  else if (!settled || !(drop >= 3)) out.warnings.push('The solution is not fully converged yet: the residual has dropped less than three decades or the lift is still changing. Run more iterations before using the numbers.');
  const idx = thin(it.map((_, i) => i)), x = idx.map((i) => it[i]);
  if (res.length) out.plots.push({ type: 'line', title: 'Residual history', xlabel: 'Iteration', ylabel: 'log10 of RMS residual', series: res.slice(0, 6).map((c) => { const v = col(t, c); return { name: c.replace(/^rms\[|\]$/g, ''), x, y: idx.map((i) => v[i]) }; }) });
  const from = Math.min(n - 1, Math.round(0.1 * n)), idx2 = thin(it.map((_, i) => i).slice(from)), x2 = idx2.map((i) => it[i]);
  out.plots.push({ type: 'line', title: 'Lift coefficient convergence', xlabel: 'Iteration', ylabel: 'CL [-]', series: [{ name: 'CL', x: x2, y: idx2.map((i) => CL[i]) }] },
    { type: 'line', title: 'Drag coefficient convergence', xlabel: 'Iteration', ylabel: 'CD [-]', series: [{ name: 'CD', x: x2, y: idx2.map((i) => CD[i]) }] });
  out.tables.push({ title: 'Final coefficients (SU2)', columns: ['Quantity', 'Value'], rows: t.columns.filter((c) => /^C[A-Z]/.test(c) || c === 'RefForce').map((c) => [c, last(col(t, c))]) });
  return out;
}

/** SU2 surface_flow.csv: the solution on the wall. Pressure coefficient is plotted when the file carries it. */
export function parseSu2SurfaceCsv(text) {
  const t = parseCsv(text), out = empty(), x = col(t, 'x');
  if (!x) throw new Error('This does not look like an SU2 surface file (no x column).');
  const x0 = Math.min(...x), c = Math.max(...x) - x0 || 1, xc = x.map((v) => (v - x0) / c);
  let name = 'Pressure_Coefficient', y = col(t, name);
  if (!y) { name = 'Pressure'; y = col(t, name); }
  if (!y) { const r = col(t, 'Density'), mx = col(t, 'Momentum_x'), my = col(t, 'Momentum_y'), mz = col(t, 'Momentum_z'), E = col(t, 'Energy'); if (r && mx && my && E) { name = 'Pressure (from the conservative variables, γ = 1.4)'; y = r.map((rr, i) => 0.4 * (E[i] - (0.5 * (mx[i] ** 2 + my[i] ** 2 + (mz ? mz[i] ** 2 : 0))) / rr)); } }
  if (y) {
    const cp = name === 'Pressure_Coefficient';
    out.plots.push({ type: 'line', title: cp ? 'Surface pressure coefficient' : 'Surface pressure (solver units)', xlabel: 'x / c [-]', ylabel: cp ? '−Cp [-]' : name, series: [{ name: cp ? '−Cp' : 'p', x: xc, y: cp ? y.map((v) => -v) : y, style: 'points' }] });
    if (cp) { out.values.cp_min = Math.min(...y); out.kpis.push(kpi('cp_min', 'Minimum pressure coefficient', out.values.cp_min)); }
  }
  out.kpis.push(kpi('surface_points', 'Surface points read', x.length, ''));
  if (!col(t, 'Pressure_Coefficient')) out.warnings.push('The surface CSV holds the raw solution only. Drop surface_flow.vtu as well to get pressure coefficient, skin friction and y+.');
  return out;
}

// --- VTK XML unstructured grid (.vtu): inline ASCII or raw appended binary, as SU2 writes it -----------
const VT = { Float32: [4, 'getFloat32'], Float64: [8, 'getFloat64'], Int8: [1, 'getInt8'], UInt8: [1, 'getUint8'], Int16: [2, 'getInt16'], UInt16: [2, 'getUint16'], Int32: [4, 'getInt32'], UInt32: [4, 'getUint32'], Int64: [8, 'getBigInt64'], UInt64: [8, 'getBigUint64'] };
function indexOfAscii(b, s, from = 0) { const n = s.length; outer: for (let i = from; i <= b.length - n; i++) { for (let k = 0; k < n; k++) if (b[i + k] !== s.charCodeAt(k)) continue outer; return i; } return -1; }
/** Read points, cells and point data of a .vtu file. Compressed and base64 files are rejected with a clear message. */
export function parseVtu(input, { maxBytes = 400e6 } = {}) {
  const b = typeof input === 'string' ? new TextEncoder().encode(input) : input instanceof Uint8Array ? input : new Uint8Array(input);
  if (b.length > maxBytes) throw new Error('This VTU file is too large to read in the browser.');
  const ap = indexOfAscii(b, '<AppendedData'), head = new TextDecoder('latin1').decode(ap < 0 ? b : b.subarray(0, ap));
  const root = /<VTKFile\b([^>]*)>/.exec(head);
  if (!root || !/UnstructuredGrid/.test(root[1])) throw new Error('This is not a VTK unstructured-grid (.vtu) file.');
  if (/compressor\s*=/.test(root[1])) throw new Error('This VTU file is compressed, which is not read here. Write it uncompressed (SU2 does by default).');
  if (/byte_order\s*=\s*"BigEndian"/.test(root[1])) throw new Error('Big-endian VTU files are not read here.');
  const hdr = /header_type\s*=\s*"UInt64"/.test(root[1]) ? 8 : 4, piece = /<Piece\b([^>]*)>/.exec(head), nP = Number(/NumberOfPoints\s*=\s*"(\d+)"/.exec(piece?.[1] || '')?.[1] || 0);
  let dataStart = -1;
  if (ap >= 0) { const tagEnd = b.indexOf(0x3e, ap); if (!/encoding\s*=\s*"raw"/.test(new TextDecoder('latin1').decode(b.subarray(ap, tagEnd + 1)))) throw new Error('Only raw appended data is read; this VTU file uses base64.'); dataStart = b.indexOf(0x5f, tagEnd) + 1; }
  const dv = new DataView(b.buffer, b.byteOffset, b.byteLength), section = (tag) => { const a = head.indexOf(`<${tag}`), e = head.indexOf(`</${tag}>`); return a < 0 ? [-1, -1] : [a, e < 0 ? head.length : e]; };
  const S = { Points: section('Points'), Cells: section('Cells'), PointData: section('PointData') }, out = { nPoints: nP, points: null, conn: null, offsets: null, types: null, pointData: {} };
  const re = /<DataArray\b([^>]*?)(\/>|>([\s\S]*?)<\/DataArray>)/g; let m;
  while ((m = re.exec(head))) {
    const at = {}; for (const a of m[1].matchAll(/([\w:]+)\s*=\s*"([^"]*)"/g)) at[a[1]] = a[2];
    const where = Object.keys(S).find((k) => m.index > S[k][0] && m.index < S[k][1]); if (!where) continue;
    const ty = VT[at.type]; if (!ty) continue;
    let data;
    if (at.format === 'appended') {
      if (dataStart <= 0) throw new Error('The VTU file refers to appended data that is missing.');
      const o = dataStart + Number(at.offset || 0), nb = hdr === 8 ? Number(dv.getBigUint64(o, true)) : dv.getUint32(o, true), n = nb / ty[0];
      if (o + hdr + nb > b.length) throw new Error('The VTU file is truncated.');
      data = new Float64Array(n); const big = ty[1].includes('Big');
      for (let i = 0; i < n; i++) data[i] = big ? Number(dv[ty[1]](o + hdr + i * ty[0], true)) : ty[0] === 1 ? dv[ty[1]](o + hdr + i) : dv[ty[1]](o + hdr + i * ty[0], true);
    } else if (!at.format || at.format === 'ascii') data = Float64Array.from((m[3] || '').trim().split(/\s+/).filter(Boolean).map(Number));
    else throw new Error('Inline base64 VTU data is not read here. Write ASCII or raw appended data.');
    if (where === 'Points') out.points = data;
    else if (where === 'Cells') { if (at.Name === 'connectivity') out.conn = data; else if (at.Name === 'offsets') out.offsets = data; else if (at.Name === 'types') out.types = data; }
    else out.pointData[at.Name || 'unnamed'] = { nc: Number(at.NumberOfComponents || 1), data };
  }
  if (!out.points) throw new Error('The VTU file has no point coordinates.');
  if (!out.nPoints) out.nPoints = out.points.length / 3;
  return out;
}
const comp = (fld, k) => { if (!fld) return null; const n = fld.data.length / fld.nc, o = new Array(n); for (let i = 0; i < n; i++) o[i] = fld.data[i * fld.nc + k]; return o; };

/** SU2 surface_flow.vtu: pressure coefficient, skin friction and y+ along the wall. */
export function su2SurfaceFromVtu(v) {
  const out = empty(), n = v.nPoints, X = [], Y = [], Z = [];
  for (let i = 0; i < n; i++) { X.push(v.points[3 * i]); Y.push(v.points[3 * i + 1]); Z.push(v.points[3 * i + 2]); }
  const planar = Z.every((z) => Math.abs(z) < 1e-12);
  // 2-D: walk the wall in connectivity order starting from the trailing edge; 3-D: keep the points of the station nearest mid-span
  let order = [];
  if (planar && v.conn && v.conn.length === 2 * (v.offsets?.length || 0)) {
    const nb = Array.from({ length: n }, () => []); for (let e = 0; e < v.conn.length; e += 2) { nb[v.conn[e]].push(v.conn[e + 1]); nb[v.conn[e + 1]].push(v.conn[e]); }
    let cur = 0; for (let i = 1; i < n; i++) if (X[i] > X[cur]) cur = i;
    let prev = -1; for (let k = 0; k < n; k++) { order.push(cur); const nx = nb[cur].find((q) => q !== prev && !(k > 1 && q === order[0])); if (nx === undefined) break; prev = cur; cur = nx; }
  } else if (!planar) { const ys = [...new Set(Y.map((y) => Number(y.toPrecision(9))))].sort((a, b) => a - b), ym = ys[Math.floor(ys.length / 2)]; order = X.map((_, i) => i).filter((i) => Math.abs(Y[i] - ym) < 1e-7 * (1 + Math.abs(ym))).sort((a, b) => X[a] - X[b]); }
  if (order.length < 3) order = X.map((_, i) => i).sort((a, b) => X[a] - X[b]);
  const x0 = Math.min(...X), c = Math.max(...X) - x0 || 1, xc = order.map((i) => (X[i] - x0) / c), style = planar && order.length === n ? 'line' : 'points';
  out.values.chord_m = c; out.values.x_le_m = x0;
  const cp = comp(v.pointData.Pressure_Coefficient, 0), cf = comp(v.pointData.Skin_Friction_Coefficient, 0), yp = comp(v.pointData.Y_Plus, 0);
  if (cp) { out.values.cp_min = Math.min(...cp); out.kpis.push(kpi('cp_min', 'Minimum pressure coefficient', out.values.cp_min)); out.plots.push({ type: 'line', title: 'Surface pressure coefficient', xlabel: 'x / c [-]', ylabel: '−Cp [-]', series: [{ name: '−Cp', x: xc, y: order.map((i) => -cp[i]), style }] }); }
  if (cf && cf.some((q) => q !== 0)) out.plots.push({ type: 'line', title: 'Skin-friction coefficient (streamwise)', xlabel: 'x / c [-]', ylabel: 'Cf,x [-]', series: [{ name: 'Cf,x', x: xc, y: order.map((i) => cf[i]), style }], annotations: [{ y: 0, label: 'Separation below this line' }] });
  if (yp && yp.some((q) => q > 0)) {
    const mx = Math.max(...yp), mean = yp.reduce((a, q) => a + q, 0) / yp.length; out.values.yplus_max = mx; out.values.yplus_mean = mean;
    out.kpis.push(kpi('yplus_max', 'Largest y+ on the wall', mx, '-', mx <= 5 ? 'ok' : 'warn', 'About 1 (at most 5) is needed when the boundary layer is resolved without wall functions'), kpi('yplus_mean', 'Mean y+ on the wall', mean));
    out.plots.push({ type: 'line', title: 'Wall y+', xlabel: 'x / c [-]', ylabel: 'y+ [-]', series: [{ name: 'y+', x: xc, y: order.map((i) => yp[i]), style }] });
  }
  if (!cp && !yp) out.warnings.push('The surface VTU file has no pressure-coefficient or y+ field.');
  return out;
}
/** SU2 flow.vtu (2-D): contour of Mach number or pressure coefficient in a window around the section. */
export function su2FieldFromVtu(v, { chord = null, xle = 0, maxTris = 24000 } = {}) {
  const out = empty(), P = v.points, n = v.nPoints;
  for (let i = 2; i < P.length; i += 3) if (Math.abs(P[i]) > 1e-12) { out.warnings.push('The volume file is three-dimensional; contour plots are drawn for 2-D section runs only. Open it in ParaView for 3-D views.'); return out; }
  if (!v.conn || !v.offsets || !v.types) return out;
  const names = ['Mach', 'Pressure_Coefficient', 'Pressure'], name = names.find((k) => v.pointData[k]); if (!name) return out;
  const val = comp(v.pointData[name], 0), c = chord || 1;
  let win = { x0: xle - 0.6 * c, x1: xle + 1.9 * c, y: 0.9 * c }, tris = [];
  for (let attempt = 0; attempt < 6; attempt++) {
    tris = []; let s = 0;
    for (let e = 0; e < v.offsets.length; e++) {
      const en = v.offsets[e], ty = v.types[e], ids = []; for (let k = s; k < en; k++) ids.push(v.conn[k]); s = en;
      if (ty !== 9 && ty !== 5) continue;
      if (!ids.every((i) => P[3 * i] >= win.x0 && P[3 * i] <= win.x1 && Math.abs(P[3 * i + 1]) <= win.y)) continue;
      tris.push([ids[0], ids[1], ids[2]]); if (ty === 9) tris.push([ids[0], ids[2], ids[3]]);
    }
    if (tris.length <= maxTris) break;
    win = { x0: xle - 0.6 * c * 0.8 ** (attempt + 1), x1: xle + c + 0.9 * c * 0.8 ** (attempt + 1), y: 0.9 * c * 0.8 ** (attempt + 1) };
  }
  if (!tris.length || tris.length > maxTris) return out;
  const map = new Map(), nodes = [], values = [];
  const id = (i) => { let q = map.get(i); if (q === undefined) { q = nodes.length; map.set(i, q); nodes.push([P[3 * i], P[3 * i + 1]]); values.push(val[i]); } return q; };
  out.plots.push({ type: 'tri', title: `${name.replace(/_/g, ' ')} around the section`, xlabel: 'x [m]', ylabel: 'y [m]', zlabel: name === 'Mach' ? 'Mach [-]' : name === 'Pressure' ? 'p (solver units)' : 'Cp [-]', equalAspect: true, edges: false, diverging: name === 'Pressure_Coefficient', nodes, tris: tris.map((t) => t.map(id)), values });
  if (name === 'Mach') { let mx = -Infinity; for (let i = 0; i < n; i++) if (val[i] > mx) mx = val[i]; out.values.mach_max = mx; out.kpis.push(kpi('mach_max', 'Highest local Mach number', mx, '-', mx > 1 ? 'warn' : 'ok', mx > 1 ? 'Supersonic region present: expect a shock and wave drag' : 'Subsonic everywhere')); }
  return out;
}

// ---------------------------------------------------------------------------------------------
// OpenFOAM
// ---------------------------------------------------------------------------------------------
/** A postProcessing .dat table: '#' header lines (the last one names the columns), bracketed vectors flattened. */
export function parseFoamDat(text) {
  const lines = asText(text).split(/\r?\n/), rows = []; let columns = [], meta = {};
  for (const l of lines) {
    const s = l.trim(); if (!s) continue;
    if (s.startsWith('#')) { const body = s.slice(1).trim(), kv = /^([\w ()]+?)\s*:\s*(.+)$/.exec(body); if (kv) meta[kv[1].trim()] = kv[2].trim(); else if (body) columns = body.split(/\s+/).filter(Boolean); continue; }
    const r = s.replace(/[()]/g, ' ').trim().split(/\s+/).map((x) => { const v = Number(x); return Number.isNaN(v) && !/nan/i.test(x) ? x : v; });
    rows.push(r);
  }
  return { columns, rows, meta };
}
const tailMean = (a, frac = 0.2) => { const k = Math.max(1, Math.round(frac * a.length)), s = a.slice(-k); return s.reduce((p, q) => p + q, 0) / s.length; };

/** forceCoeffs function object (coefficient.dat in current releases, forceCoeffs.dat in older ones). */
export function parseFoamCoeffs(text) {
  const t = parseFoamDat(text), out = empty(), time = col(t, 'Time'), Cl = col(t, 'Cl'), Cd = col(t, 'Cd'), Cm = col(t, 'CmPitch') || col(t, 'Cm');
  if (!time || !Cl || !Cd || !t.rows.length) throw new Error('This does not look like an OpenFOAM force-coefficient file (Time, Cd and Cl columns expected).');
  const n = time.length, cl = last(Cl), cd = last(Cd), tail = Math.max(2, Math.round(0.1 * n)), sl = Cl.slice(-tail), drift = Math.max(...sl) - Math.min(...sl), settled = drift < 2e-3 * Math.max(0.05, Math.abs(cl));
  out.values = { CL: cl, CD: cd, LD: cd !== 0 ? cl / cd : NaN, iterations: n, CL_mean_last20pct: tailMean(Cl), CD_mean_last20pct: tailMean(Cd), CL_drift: drift };
  if (Cm) out.values.Cm = last(Cm);
  out.kpis.push(kpi('CL', 'Lift coefficient', cl, '-', !finite(cl) ? 'bad' : settled ? 'ok' : 'warn', settled ? 'Steady over the last tenth of the run' : `Still moving: changed by ${drift.toPrecision(2)} over the last tenth of the run`),
    kpi('CD', 'Drag coefficient', cd), kpi('LD', 'Lift-to-drag ratio', out.values.LD));
  if (Cm) out.kpis.push(kpi('Cm', 'Pitching-moment coefficient', last(Cm)));
  out.kpis.push(kpi('CL_mean_last20pct', 'Mean lift coefficient over the last fifth', out.values.CL_mean_last20pct), kpi('CD_mean_last20pct', 'Mean drag coefficient over the last fifth', out.values.CD_mean_last20pct), kpi('iterations', 'Steps written', n, ''));
  if (!settled) out.warnings.push('The force coefficients are still changing. For a steady run, continue it; for a time-accurate run, use the mean over the last part.');
  const from = Math.min(n - 1, Math.round(0.1 * n)), idx = thin(time.map((_, i) => i).slice(from)), x = idx.map((i) => time[i]);
  out.plots.push({ type: 'line', title: 'Lift coefficient history', xlabel: 'Iteration or time [s]', ylabel: 'CL [-]', series: [{ name: 'CL', x, y: idx.map((i) => Cl[i]) }] }, { type: 'line', title: 'Drag coefficient history', xlabel: 'Iteration or time [s]', ylabel: 'CD [-]', series: [{ name: 'CD', x, y: idx.map((i) => Cd[i]) }] });
  out.tables.push({ title: 'Final coefficients (OpenFOAM)', columns: ['Quantity', 'Value'], rows: t.columns.slice(1).map((c, j) => [c, last(t.rows)[j + 1]]) });
  return out;
}
/** forces function object (force.dat): final force components. */
export function parseFoamForces(text) {
  const t = parseFoamDat(text), out = empty(); if (!t.rows.length) return out;
  const r = last(t.rows), named = t.columns.length === r.length ? t.columns : ['Time', 'total_x', 'total_y', 'total_z', 'pressure_x', 'pressure_y', 'pressure_z', 'viscous_x', 'viscous_y', 'viscous_z'].slice(0, r.length);
  out.tables.push({ title: 'Final forces on the body (OpenFOAM)', columns: ['Component', 'Force [N]'], rows: named.slice(1).map((c, j) => [c, r[j + 1]]) });
  const g = (k) => r[named.indexOf(k)]; if (finite(g('total_x'))) { out.values.force_x_N = g('total_x'); out.values.force_y_N = g('total_y'); out.values.force_z_N = g('total_z'); }
  return out;
}
/** solverInfo function object: initial residual of each equation per step. */
export function parseFoamSolverInfo(text) {
  const t = parseFoamDat(text), out = empty(), time = col(t, 'Time'); if (!time || !t.rows.length) return out;
  const init = t.columns.filter((c) => /_initial$/.test(c)).slice(0, 6), idx = thin(time.map((_, i) => i)), x = idx.map((i) => time[i]);
  if (!init.length) return out;
  out.plots.push({ type: 'line', title: 'Residual history', xlabel: 'Iteration or time [s]', ylabel: 'Initial residual [-]', ylog: true, series: init.map((c) => { const v = col(t, c); return { name: c.replace(/_initial$/, ''), x, y: idx.map((i) => (v[i] > 0 ? v[i] : NaN)) }; }) });
  const p = col(t, 'p_initial'); if (p) { out.values.p_residual = last(p); out.kpis.push(kpi('p_residual', 'Final pressure residual', last(p), '-', last(p) < 1e-3 ? 'ok' : 'warn', 'Below 1e-3 is a common minimum for a steady run')); }
  return out;
}
/** yPlus function object. */
export function parseFoamYPlus(text) {
  const t = parseFoamDat(text), out = empty(); if (!t.rows.length) return out;
  const r = last(t.rows), nums = r.filter(finite); if (nums.length < 4) return out;
  const [mn, mx, av] = nums.slice(-3); out.values.yplus_max = mx; out.values.yplus_mean = av;
  out.kpis.push(kpi('yplus_max', 'Largest y+ on the wall', mx, '-', mx <= 300 ? 'ok' : 'warn', 'Wall functions are valid for roughly 30 to 300'), kpi('yplus_mean', 'Mean y+ on the wall', av), kpi('yplus_min', 'Smallest y+ on the wall', mn));
  return out;
}
/** Solver log (log.simpleFoam …): first initial residual of every equation in each step, and the mesh size if printed. */
export function parseFoamLog(text) {
  const out = empty(), src = asText(text), series = new Map(), times = []; let cur = null, seen = null;
  for (const l of src.split(/\r?\n/)) {
    const tm = /^Time = ([-+0-9.eE]+)/.exec(l); if (tm) { cur = Number(tm[1]); times.push(cur); seen = new Set(); continue; }
    if (cur == null) continue;
    const m = /Solving for (\w+), Initial residual = ([-+0-9.eE]+)/.exec(l); if (!m || seen.has(m[1])) continue;
    seen.add(m[1]); if (!series.has(m[1])) series.set(m[1], []); series.get(m[1]).push([cur, Number(m[2])]);
  }
  if (!times.length) return out;
  out.values.iterations = times.length; out.kpis.push(kpi('iterations', 'Steps in the log', times.length, ''));
  const names = [...series.keys()].slice(0, 6);
  if (names.length) out.plots.push({ type: 'line', title: 'Residual history (from the solver log)', xlabel: 'Iteration or time [s]', ylabel: 'Initial residual [-]', ylog: true, series: names.map((k) => { const a = thin(series.get(k)); return { name: k, x: a.map((q) => q[0]), y: a.map((q) => (q[1] > 0 ? q[1] : NaN)) }; }) });
  if (/FOAM FATAL/.test(src)) out.warnings.push('The log contains a FOAM FATAL ERROR: the run stopped early. Read the end of the log for the reason.');
  if (/Floating point exception|nan/i.test(src.slice(-4000)) && !/End\s*$/.test(src.trim())) out.warnings.push('The run appears to have diverged (floating-point exception at the end of the log).');
  return out;
}
/** checkMesh log: cell count and whether the mesh passed. */
export function parseFoamCheckMesh(text) {
  const out = empty(), src = asText(text), cells = /^\s*cells:\s+(\d+)/m.exec(src), fail = /Failed (\d+) mesh checks/.exec(src), non = /Max non-orthogonality.*?=\s*([\d.]+)/i.exec(src) || /non-orthogonality Max:\s*([\d.]+)/i.exec(src);
  if (cells) { out.values.mesh_cells = Number(cells[1]); out.kpis.push(kpi('mesh_cells', 'Mesh cells', Number(cells[1]), '')); }
  if (non) out.kpis.push(kpi('mesh_nonortho_max', 'Largest non-orthogonality', Number(non[1]), 'deg', Number(non[1]) <= 70 ? 'ok' : 'warn', 'Above about 70° the solution loses accuracy'));
  if (fail) out.warnings.push(`checkMesh reported ${fail[1]} failed mesh check(s). The solution may still run, but inspect the mesh before trusting forces.`);
  return out;
}

// ---------------------------------------------------------------------------------------------
// CalculiX
// ---------------------------------------------------------------------------------------------
/** CalculiX .dat: printed displacements, reaction totals, eigenfrequencies and buckling factors. */
export function parseCcxDat(text) {
  const out = empty(), L = asText(text).split(/\r?\n/), nums = (l) => l.trim().split(/\s+/).map(Number), disp = [], freqs = [], buck = []; let react = null;
  for (let i = 0; i < L.length; i++) {
    const l = L[i]; let m;
    if ((m = /displacements \(vx,vy,vz\) for set (\S+) and time\s+(\S+)/.exec(l))) {
      const rows = []; let j = i + 1; while (j < L.length && !L[j].trim()) j++;
      for (; j < L.length && L[j].trim(); j++) { const r = nums(L[j]); if (r.length < 4 || r.some(Number.isNaN)) break; rows.push(r); }
      disp.push({ set: m[1], time: Number(m[2]), rows }); i = j - 1;
    } else if ((m = /total force \(fx,fy,fz\) for set (\S+) and time/.exec(l))) { let j = i + 1; while (j < L.length && !L[j].trim()) j++; const r = nums(L[j] || ''); if (r.length >= 3) react = { set: m[1], f: r.slice(0, 3) }; i = j; }
    else if (/E I G E N V A L U E\s+O U T P U T/.test(l)) { for (let j = i + 1; j < L.length; j++) { const r = nums(L[j]); if (r.length >= 4 && r.every(finite) && Number.isInteger(r[0])) freqs.push({ mode: r[0], eigenvalue: r[1], omega: r[2], f: r[3] }); else if (freqs.length && L[j].trim() === '') break; else if (/[A-Z] [A-Z] [A-Z]/.test(L[j]) && freqs.length) break; } }
    else if (/B U C K L I N G\s+F A C T O R\s+O U T P U T/.test(l)) { for (let j = i + 1; j < L.length; j++) { const r = nums(L[j]); if (r.length === 2 && r.every(finite)) buck.push({ mode: r[0], factor: r[1] }); else if (buck.length && L[j].trim() === '') break; } }
  }
  const times = [...new Set(disp.map((d) => d.time))];
  if (disp.length && times.length <= 2) {
    const d = disp.find((q) => /^TIP$/i.test(q.set)) || disp[0], uz = d.rows.map((r) => r[3]), mean = uz.reduce((a, b) => a + b, 0) / uz.length, peak = uz.reduce((a, b) => (Math.abs(b) > Math.abs(a) ? b : a), 0);
    out.values.tip_deflection_m = mean; out.values.tip_deflection_peak_m = peak; out.values.tip_twist_spread_m = Math.max(...uz) - Math.min(...uz);
    out.kpis.push(kpi('tip_deflection_m', `Tip deflection (mean of node set ${d.set})`, mean, 'm'), kpi('tip_deflection_peak_m', 'Largest tip-node deflection', peak, 'm'), kpi('tip_twist_spread_m', 'Vertical spread across the tip section (twist)', out.values.tip_twist_spread_m, 'm'));
  } else if (disp.length) {
    const t = disp.map((d) => d.time), uz = disp.map((d) => d.rows[0]?.[3] ?? NaN), idx = thin(t.map((_, i) => i)), peak = uz.reduce((a, b) => (Math.abs(b) > Math.abs(a) ? b : a), 0);
    out.values.tip_deflection_dynamic_peak_m = peak; out.kpis.push(kpi('tip_deflection_dynamic_peak_m', 'Peak dynamic tip deflection', peak, 'm'), kpi('dynamic_steps', 'Time points printed', t.length, ''));
    out.plots.push({ type: 'line', title: 'Tip response to the load pulse', xlabel: 'Time [s]', ylabel: 'Tip deflection [m]', series: [{ name: `Node set ${disp[0].set}`, x: idx.map((i) => t[i]), y: idx.map((i) => uz[i]) }] });
  }
  if (react) { out.values.root_reaction_z_N = react.f[2]; out.kpis.push(kpi('root_reaction_z_N', `Vertical reaction at ${react.set}`, react.f[2], 'N', undefined, 'Equal and opposite to the load applied to the mesh')); }
  if (freqs.length) {
    freqs.forEach((q, i) => { if (i < 3) out.values[`f${i + 1}_Hz`] = q.f; });
    out.values.frequencies_Hz = freqs.map((q) => q.f);
    out.kpis.push(...freqs.slice(0, 4).map((q, i) => kpi(`mode${i + 1}_Hz`, `Natural frequency ${i + 1}`, q.f, 'Hz')));
    out.tables.push({ title: 'Natural frequencies (CalculiX)', columns: ['Mode', 'Frequency [Hz]', 'Angular frequency [rad/s]', 'Eigenvalue [rad²/s²]'], rows: freqs.map((q) => [q.mode, q.f, q.omega, q.eigenvalue]) });
    out.plots.push({ type: 'bar', title: 'Natural frequencies', ylabel: 'Frequency [Hz]', categories: freqs.map((q) => `Mode ${q.mode}`), series: [{ name: 'CalculiX', y: freqs.map((q) => q.f) }] });
  }
  if (buck.length) {
    const pos = buck.filter((q) => q.factor > 0), first = (pos[0] || buck[0]).factor;
    out.values.buckling_load_factor = first; out.values.buckling_factors = buck.map((q) => q.factor);
    out.kpis.push(kpi('buckling_load_factor', 'First buckling factor (critical load ÷ applied load)', first, '-', first >= 1.5 ? 'ok' : first >= 1 ? 'warn' : 'bad', 'Below 1 the structure buckles before the applied load is reached'));
    out.tables.push({ title: 'Buckling factors (CalculiX)', columns: ['Mode', 'Critical load ÷ applied load'], rows: buck.map((q) => [q.mode, q.factor]) });
  }
  if (!disp.length && !freqs.length && !buck.length && !react) out.warnings.push('No printed results were found in this .dat file.');
  return out;
}

/**
 * CalculiX .frd (ASCII): nodal displacements and stresses of the static step, and the shape of every mode.
 * Modes are classified from their shapes (vertical bending, in-plane bending, torsion or local panel) so that
 * like is compared with like.
 */
export function parseCcxFrd(text) {
  const out = empty(), L = asText(text).split(/\r?\n/), ids = new Map(), xyz = [];
  const fx = (l, k) => Number(l.slice(13 + 12 * k, 25 + 12 * k));
  let i = 0;
  for (; i < L.length; i++) if (/^\s{4}2C/.test(L[i])) { for (i++; i < L.length && L[i].startsWith(' -1'); i++) { ids.set(Number(L[i].slice(3, 13)), xyz.length / 3); xyz.push(fx(L[i], 0), fx(L[i], 1), fx(L[i], 2)); } break; }
  const n = xyz.length / 3; if (!n) throw new Error('This does not look like a CalculiX .frd results file (no node block).');
  const mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity]; for (let q = 0; q < n; q++) for (let k = 0; k < 3; k++) { const v = xyz[3 * q + k]; if (v < mn[k]) mn[k] = v; if (v > mx[k]) mx[k] = v; }
  const size = mx.map((v, k) => v - mn[k]), sp = size.indexOf(Math.max(...size)), up = 2, ch = [0, 1, 2].find((k) => k !== sp && k !== up) ?? 0, NB = 40;
  const bin = new Int32Array(n), cnt = new Float64Array(NB), xm = new Float64Array(NB), sxx = new Float64Array(NB);
  for (let q = 0; q < n; q++) { const b = Math.min(NB - 1, Math.floor(((xyz[3 * q + sp] - mn[sp]) / (size[sp] || 1)) * NB)); bin[q] = b; cnt[b]++; xm[b] += xyz[3 * q + ch]; }
  for (let b = 0; b < NB; b++) if (cnt[b]) xm[b] /= cnt[b];
  for (let q = 0; q < n; q++) sxx[bin[q]] += (xyz[3 * q + ch] - xm[bin[q]]) ** 2;
  const station = Array.from({ length: NB }, (_, b) => mn[sp] + ((b + 0.5) / NB) * size[sp]), modes = [], buck = [];
  let staticDisp = null, vmMax = NaN, vmNode = -1;
  for (; i < L.length; i++) {
    if (!/^\s{2}100CL/.test(L[i])) continue;
    const tk = L[i].trim().split(/\s+/), value = Number(tk[2]), ictype = Number(tk[4]), nameLine = L[i + 1] || '', name = nameLine.startsWith(' -4') ? nameLine.slice(5, 13).trim() : '';
    let j = i + 2; while (j < L.length && L[j].startsWith(' -5')) j++;
    if (name === 'DISP') {
      const suz = new Float64Array(NB), sux = new Float64Array(NB), su2 = new Float64Array(NB), sxu = new Float64Array(NB); let umax = 0, zmax = 0;
      for (; j < L.length && L[j].startsWith(' -1'); j++) { const q = ids.get(Number(L[j].slice(3, 13))); if (q === undefined) continue; const u = [fx(L[j], 0), fx(L[j], 1), fx(L[j], 2)], b = bin[q], m2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]; suz[b] += u[up]; sux[b] += u[ch]; su2[b] += m2; sxu[b] += (xyz[3 * q + ch] - xm[b]) * u[up]; if (m2 > umax) umax = m2; if (Math.abs(u[up]) > Math.abs(zmax)) zmax = u[up]; }
      let eT = 0, eB = 0, eL = 0, eR = 0; for (let b = 0; b < NB; b++) if (cnt[b]) { eT += su2[b]; eB += suz[b] ** 2 / cnt[b]; eL += sux[b] ** 2 / cnt[b]; if (sxx[b] > 0) eR += sxu[b] ** 2 / sxx[b]; }
      const fr = eT > 0 ? { bending: eB / eT, inplane: eL / eT, torsion: eR / eT } : { bending: 0, inplane: 0, torsion: 0 }, best = Object.entries(fr).sort((a, b) => b[1] - a[1])[0];
      const kind = best[1] >= 0.5 ? best[0] : 'local', mean = Array.from(suz, (v, b) => (cnt[b] ? v / cnt[b] : NaN));
      const rec = { value, kind, fractions: fr, shape: mean, umax: Math.sqrt(umax), zmax };
      if (ictype === 2) modes.push(rec); else if (ictype === 4) { if (value !== 0) buck.push(rec); } else if (!staticDisp || ictype === 0) staticDisp = rec;
    } else if (name === 'STRESS') {
      for (; j < L.length && L[j].startsWith(' -1'); j++) { const s = [0, 1, 2, 3, 4, 5].map((k) => fx(L[j], k)), vm = Math.sqrt(0.5 * ((s[0] - s[1]) ** 2 + (s[1] - s[2]) ** 2 + (s[2] - s[0]) ** 2) + 3 * (s[3] ** 2 + s[4] ** 2 + s[5] ** 2)); if (!(vm <= vmMax)) { vmMax = vm; vmNode = Number(L[j].slice(3, 13)); } }
    }
    i = j - 1;
  }
  const axis = 'xyz'[sp];
  if (staticDisp && !modes.length && !buck.length) {
    out.values.tip_deflection_max_m = staticDisp.zmax; out.values.displacement_max_m = staticDisp.umax;
    out.kpis.push(kpi('tip_deflection_max_m', 'Largest vertical displacement', staticDisp.zmax, 'm'), kpi('displacement_max_m', 'Largest displacement magnitude', staticDisp.umax, 'm'));
    out.plots.push({ type: 'line', title: 'Deflected shape along the span', xlabel: `Spanwise position ${axis} [m]`, ylabel: 'Vertical deflection [m]', series: [{ name: 'Section mean', x: station, y: staticDisp.shape }] });
  }
  if (finite(vmMax)) { out.values.sigma_max_Pa = vmMax; out.kpis.push(kpi('sigma_max_Pa', 'Peak von Mises stress', vmMax, 'Pa', undefined, `At node ${vmNode}; nodal values extrapolated from the integration points`)); }
  const LABEL = { bending: 'vertical bending', inplane: 'in-plane (fore-aft) bending', torsion: 'torsion', local: 'local panel mode' };
  if (modes.length) {
    const bend = modes.filter((m) => m.kind === 'bending'), tors = modes.find((m) => m.kind === 'torsion');
    bend.slice(0, 3).forEach((m, k) => { out.values[`f${k + 1}_Hz`] = m.value; out.kpis.push(kpi(`f${k + 1}_Hz`, `Vertical bending frequency ${k + 1}`, m.value, 'Hz')); });
    if (tors) { out.values.f_torsion_Hz = tors.value; out.kpis.push(kpi('f_torsion_Hz', 'First torsion frequency', tors.value, 'Hz')); }
    out.values.frequencies_Hz = modes.map((m) => m.value); out.values.mode_kinds = modes.map((m) => m.kind);
    out.tables.push({ title: 'Modes identified from their shapes', columns: ['Mode', 'Frequency [Hz]', 'Character', 'Bending share', 'In-plane share', 'Torsion share'], rows: modes.map((m, k) => [k + 1, m.value, LABEL[m.kind], m.fractions.bending, m.fractions.inplane, m.fractions.torsion]) });
    const shown = (bend.length ? bend : modes).slice(0, 4);
    out.plots.push({ type: 'line', title: 'Mode shapes along the span (normalised)', xlabel: `Spanwise position ${axis} [m]`, ylabel: 'Vertical displacement, normalised [-]', series: shown.map((m) => { const pk = Math.max(...m.shape.filter(finite).map(Math.abs)) || 1, sg = Math.sign(last(m.shape.filter(finite)) || 1); return { name: `${m.value.toPrecision(4)} Hz`, x: station, y: m.shape.map((v) => (sg * v) / pk) }; }) });
    if (!bend.length) out.warnings.push('None of the computed modes is a clean vertical bending mode; ask for more modes or check the model.');
    if (modes.some((m) => m.kind === 'local')) out.warnings.push('Some modes are local skin-panel modes. They do not exist in a beam model, so only the bending and torsion modes are compared with the native suites.');
  }
  if (buck.length) { const pos = buck.filter((b) => b.value > 0), first = (pos[0] || buck[0]).value; out.values.buckling_load_factor = first; out.tables.push({ title: 'Buckling modes', columns: ['Mode', 'Critical load ÷ applied load', 'Character'], rows: buck.map((b, k) => [k + 1, b.value, LABEL[b.kind]]) }); }
  out.kpis.push(kpi('frd_nodes', 'Result nodes read', n, ''));
  return out;
}

// ---------------------------------------------------------------------------------------------
// Dispatcher: a set of dropped files → one merged result
// ---------------------------------------------------------------------------------------------
/** What a file is, from its name and (when the name is not enough) its first bytes. */
export function detectOutput(name, bytes) {
  const b = base(name).toLowerCase(), p = String(name).replace(/\\/g, '/').toLowerCase(), headTxt = () => dec.decode(bytes.subarray(0, Math.min(bytes.length, 3000)));
  if (b === 'manifest.json') return 'manifest';
  if (/\.vtu$/.test(b)) return /surface/.test(b) ? 'su2-surface-vtu' : 'su2-volume-vtu';
  if (/\.frd$/.test(b)) return 'ccx-frd';
  if (/\.csv$/.test(b)) { const h = headTxt(); if (/"?Inner_Iter"?|"?Time_Iter"?/.test(h)) return 'su2-history'; if (/"?PointID"?/.test(h)) return 'su2-surface-csv'; return null; }
  if (/\.dat$/.test(b)) {
    if (/^(coefficient|forcecoeffs)/.test(b)) return 'foam-coeffs'; if (/^forces?(_\d+)?\.dat$/.test(b)) return 'foam-forces'; if (/^solverinfo/.test(b)) return 'foam-solverinfo'; if (/^yplus/.test(b)) return 'foam-yplus';
    const h = headTxt();
    if (/^#\s*Force coefficients/m.test(h) || /#\s*Time\s+.*\bCd\b/.test(h)) return 'foam-coeffs'; if (/#\s*Solver information/i.test(h)) return 'foam-solverinfo';
    if (/S T E P|E I G E N V A L U E|B U C K L I N G|displacements \(vx,vy,vz\)/.test(dec.decode(bytes.subarray(0, Math.min(bytes.length, 20000))))) return 'ccx-dat';
    if (/^#/.test(h.trim())) return /moment\.dat$/.test(b) ? null : 'foam-forces';
    return null;
  }
  if (/^log\.checkmesh/.test(b)) return 'foam-checkmesh';
  if (/^log\.(simplefoam|pimplefoam|rhosimplefoam|potentialfoam)/.test(b)) return 'foam-log';
  if (/postprocessing\//.test(p)) return null;
  return null;
}
const ORDER = ['manifest', 'su2-history', 'su2-surface-vtu', 'su2-surface-csv', 'su2-volume-vtu', 'foam-coeffs', 'foam-forces', 'foam-solverinfo', 'foam-log', 'foam-yplus', 'foam-checkmesh', 'ccx-dat', 'ccx-frd'];
const SOLVER_OF = (kind) => (kind.startsWith('su2') ? 'su2' : kind.startsWith('foam') ? 'openfoam' : kind.startsWith('ccx') ? 'calculix' : null);

/**
 * Read a set of solver output files.
 * @param {{name: string, bytes: Uint8Array}[]} files
 * @param {{manifest?: object}} [opts] the manifest of the case that produced them, when known
 * @returns {{solver, kpis, plots, tables, values, warnings, sources: {name, kind, ok, note}[], manifest}}
 */
export function readResults(files, { manifest = null } = {}) {
  const out = { solver: null, ...empty(), sources: [], manifest };
  const tagged = files.map((f) => ({ ...f, kind: detectOutput(f.name, f.bytes) })).filter((f) => { if (!f.kind) { out.sources.push({ name: f.name, kind: null, ok: false, note: 'Not a result file this page reads' }); return false; } return true; });
  tagged.sort((a, b) => ORDER.indexOf(a.kind) - ORDER.indexOf(b.kind));
  const seenKpi = new Set(), hasSolverInfo = tagged.some((f) => f.kind === 'foam-solverinfo'), hasSurfVtu = tagged.some((f) => f.kind === 'su2-surface-vtu'), ctx = { chord: null, xle: 0 };
  for (const fl of tagged) {
    try {
      let r = null;
      if (fl.kind === 'manifest') { const m = JSON.parse(asText(fl.bytes)); if (m && m.app === 'AeroSuite 26' && m.solver) { out.manifest = m; ctx.chord = m.settings?.chord_m || null; } out.sources.push({ name: fl.name, kind: fl.kind, ok: true, note: 'Case manifest' }); continue; }
      if (fl.kind === 'su2-history') r = parseSu2History(fl.bytes);
      else if (fl.kind === 'su2-surface-csv') { if (hasSurfVtu) { out.sources.push({ name: fl.name, kind: fl.kind, ok: true, note: 'Skipped: the surface VTU file carries the same points with more fields' }); continue; } r = parseSu2SurfaceCsv(fl.bytes); }
      else if (fl.kind === 'su2-surface-vtu') { r = su2SurfaceFromVtu(parseVtu(fl.bytes)); if (!ctx.chord && r.values.chord_m) ctx.chord = r.values.chord_m; ctx.xle = r.values.x_le_m ?? 0; delete r.values.chord_m; delete r.values.x_le_m; }
      else if (fl.kind === 'su2-volume-vtu') r = su2FieldFromVtu(parseVtu(fl.bytes), ctx);
      else if (fl.kind === 'foam-coeffs') r = parseFoamCoeffs(fl.bytes);
      else if (fl.kind === 'foam-forces') r = parseFoamForces(fl.bytes);
      else if (fl.kind === 'foam-solverinfo') r = parseFoamSolverInfo(fl.bytes);
      else if (fl.kind === 'foam-yplus') r = parseFoamYPlus(fl.bytes);
      else if (fl.kind === 'foam-checkmesh') r = parseFoamCheckMesh(fl.bytes);
      else if (fl.kind === 'foam-log') { r = parseFoamLog(fl.bytes); if (hasSolverInfo) r.plots = []; }
      else if (fl.kind === 'ccx-dat') r = parseCcxDat(fl.bytes);
      else if (fl.kind === 'ccx-frd') r = parseCcxFrd(fl.bytes);
      if (!r) continue;
      out.solver ||= SOLVER_OF(fl.kind);
      for (const k of r.kpis) { if (seenKpi.has(k.key)) continue; seenKpi.add(k.key); out.kpis.push(k); }
      // shape-based bending frequencies from the .frd replace the raw mode order of the .dat
      if (fl.kind === 'ccx-frd' && r.values.mode_kinds) for (const k of ['f1_Hz', 'f2_Hz', 'f3_Hz', 'f_torsion_Hz']) delete out.values[k];
      for (const [k, v] of Object.entries(r.values)) if (fl.kind === 'ccx-frd' || out.values[k] === undefined) out.values[k] = v;
      out.plots.push(...r.plots); out.tables.push(...r.tables); out.warnings.push(...r.warnings);
      out.sources.push({ name: fl.name, kind: fl.kind, ok: true, note: `${r.kpis.length} values, ${r.plots.length} charts` });
    } catch (e) { out.sources.push({ name: fl.name, kind: fl.kind, ok: false, note: e.message }); out.warnings.push(`${base(fl.name)}: ${e.message}`); }
  }
  // the .frd classifies modes; when it is present, drop the unclassified f1..f3 KPIs that came from the .dat order
  if (tagged.some((f) => f.kind === 'ccx-frd') && out.values.mode_kinds) out.kpis = out.kpis.filter((k) => !/^mode\d_Hz$/.test(k.key));
  else if (out.values.frequencies_Hz && !out.values.mode_kinds) out.warnings.push('Frequencies are listed in solver order. Drop the .frd file too so that bending, torsion and local panel modes can be told apart.');
  out.warnings = [...new Set(out.warnings)];
  return out;
}

// ---------------------------------------------------------------------------------------------
// Comparison with the native suites, and the record stored on the data bus
// ---------------------------------------------------------------------------------------------
const MAP = [
  { key: 'CL', label: 'Lift coefficient', unit: '-', suite: 'cfd', native: 'CL' },
  { key: 'CD', label: 'Drag coefficient', unit: '-', suite: 'cfd', native: 'CD' },
  { key: 'LD', label: 'Lift-to-drag ratio', unit: '-', suite: 'cfd', native: (u) => (finite(u.CL) && u.CD ? u.CL / u.CD : undefined) },
  { key: 'Cm', label: 'Pitching-moment coefficient', unit: '-', suite: 'cfd', native: 'Cm_ac' },
  { key: 'cp_min', label: 'Minimum pressure coefficient', unit: '-', suite: 'cfd', native: 'cp_min' },
  { key: 'tip_deflection_m', label: 'Tip deflection', unit: 'm', suite: 'fea', native: 'tip_deflection_m' },
  { key: 'sigma_max_Pa', label: 'Peak von Mises stress', unit: 'Pa', suite: 'fea', native: 'sigma_max_Pa' },
  { key: 'buckling_load_factor', label: 'Buckling factor', unit: '-', suite: 'fea', native: 'buckling_factor' },
  { key: 'f1_Hz', label: 'First bending frequency', unit: 'Hz', suite: 'vibration', native: 'f1_Hz' },
  { key: 'f2_Hz', label: 'Second bending frequency', unit: 'Hz', suite: 'vibration', native: 'f2_Hz' },
  { key: 'f3_Hz', label: 'Third bending frequency', unit: 'Hz', suite: 'vibration', native: 'f3_Hz' },
  { key: 'f_torsion_Hz', label: 'First torsion frequency', unit: 'Hz', suite: 'vibration', native: 'f_torsion_Hz' },
];
/**
 * Set high-fidelity values beside the native suite outputs on the data bus.
 * Each row: { key, label, unit, hifi, native, suite, diff, pct, ratio, note }. `ratio` (high-fidelity ÷ native) is the
 * factor that would calibrate the native model to the high-fidelity answer.
 */
export function compare(values, up = {}, { manifest = null } = {}) {
  const rows = [], twoD = manifest?.solver === 'su2' && manifest?.mesh?.source === 'generated' && !manifest?.mesh?.span, sf = manifest?.settings?.load_factor;
  for (const m of MAP) {
    const hifi = values[m.key]; if (!finite(hifi)) continue;
    const u = up[m.suite] || {}, native = typeof m.native === 'function' ? m.native(u) : u[m.native], has = finite(native);
    let note = has ? '' : `Run the native ${m.suite} suite to compare`;
    if (has && twoD && ['CL', 'CD', 'LD', 'Cm'].includes(m.key)) note = 'Section (2-D) value beside a whole-aircraft native value: no induced drag or finite-span lift loss in the section result';
    if (has && m.key === 'buckling_load_factor') note = 'Native factor is for a single skin panel and already divides by the ultimate safety factor; the shell model finds the first buckle anywhere under the applied load';
    if (has && m.key === 'tip_deflection_m' && finite(sf)) note = `Both at ${sf} g when the native run used the same load factor`;
    if (has && m.key === 'sigma_max_Pa') note = 'Shell peak includes local stress concentrations at ribs and the root that a beam model averages out';
    rows.push({ key: m.key, label: m.label, unit: m.unit, suite: m.suite, hifi, native: has ? native : null, diff: has ? hifi - native : null, pct: has && native !== 0 ? (100 * (hifi - native)) / Math.abs(native) : null, ratio: has && native !== 0 ? hifi / native : null, note });
  }
  return rows;
}
/** Plain-data record of a high-fidelity result for the data bus (state.up.hifi). */
export function hifiRecord(result, { ts = Date.now() } = {}) {
  const rec = {}; for (const [k, v] of Object.entries(result.values)) if (finite(v)) rec[k] = v;
  const m = result.manifest;
  return { ...rec, solver: result.solver || m?.solver || '', analysis: m?.analysis || '', caseHash: m?.caseHash || '', caseName: m?.caseName || '', ts };
}
/** Comparison rows as a bar chart of percentage differences. */
export function comparePlot(rows) {
  const r = rows.filter((x) => finite(x.pct)); if (!r.length) return null;
  return { type: 'bar', title: 'High-fidelity result relative to the native suites', ylabel: 'Difference [%]', categories: r.map((x) => x.label), series: [{ name: 'High-fidelity − native', y: r.map((x) => x.pct) }] };
}
