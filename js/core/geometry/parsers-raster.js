// GeoTIFF terrain rasters → height-field surface, and the E57 header/XML inspector.
// The TIFF decoder is native: classic TIFF and BigTIFF, strips or tiles, uncompressed / LZW / Deflate / PackBits,
// 8/16/32-bit integer and 32/64-bit float samples, horizontal and floating-point predictors.

import { guard, str, view, LIMITS } from './parsers-util.js';

const TYPE_SIZE = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8, 16: 8, 17: 8, 18: 8 };
const COMPRESSION = { 1: 'none', 5: 'LZW', 8: 'Deflate', 32946: 'Deflate', 32773: 'PackBits', 6: 'JPEG (old)', 7: 'JPEG', 34712: 'JPEG 2000', 50000: 'ZSTD', 50001: 'WebP', 34887: 'LERC', 34925: 'LZMA' };

/** TIFF LZW (MSB-first codes, "early change"). Decodes at most `expected` bytes. */
function lzw(src, expected) {
  const out = new Uint8Array(expected), prefix = new Int32Array(4096), suffix = new Uint8Array(4096), firstCh = new Uint8Array(4096), stack = new Uint8Array(4097);
  for (let i = 0; i < 256; i++) { suffix[i] = i; firstCh[i] = i; }
  let op = 0, bits = 9, next = 258, buf = 0, cnt = 0, ip = 0, prev = -1;
  const emit = (code) => { let len = 0, c = code; while (c >= 258) { if (len > 4095) throw new Error('corrupt LZW data'); stack[len++] = suffix[c]; c = prefix[c]; } stack[len++] = c; while (len && op < expected) out[op++] = stack[--len]; };
  while (op < expected) {
    while (cnt < bits && ip < src.length) { buf = ((buf << 8) | src[ip++]) & 0xffffff; cnt += 8; }
    if (cnt < bits) break;
    const code = (buf >> (cnt - bits)) & ((1 << bits) - 1); cnt -= bits;
    if (code === 257) break;
    if (code === 256) { bits = 9; next = 258; prev = -1; continue; }
    if (prev < 0) { if (code > 255) throw new Error('corrupt LZW data'); out[op++] = code; prev = code; continue; }
    if (code < next) { emit(code); if (next < 4096) { prefix[next] = prev; suffix[next] = firstCh[code]; firstCh[next] = firstCh[prev]; next++; } }
    else if (code === next && next < 4096) { prefix[next] = prev; suffix[next] = firstCh[prev]; firstCh[next] = firstCh[prev]; next++; emit(code); }
    else throw new Error('corrupt LZW data');
    prev = code;
    if (next + 1 >= 1 << bits && bits < 12) bits++;
  }
  return out;
}
function packBits(src, expected) {
  const out = new Uint8Array(expected); let op = 0, ip = 0;
  while (ip < src.length && op < expected) { const n = src[ip++] << 24 >> 24; if (n >= 0) { for (let k = 0; k <= n && ip < src.length && op < expected; k++) out[op++] = src[ip++]; } else if (n !== -128) { const v = src[ip++]; for (let k = 0; k < 1 - n && op < expected; k++) out[op++] = v; } }
  return out;
}
async function inflateZlib(data, expected) {
  if (typeof DecompressionStream === 'undefined') throw new Error('this runtime has no DecompressionStream; Deflate-compressed TIFF cannot be read');
  const ds = new DecompressionStream('deflate'), w = ds.writable.getWriter(); w.write(data).catch(() => {}); w.close().catch(() => {});
  const r = ds.readable.getReader(), out = new Uint8Array(expected); let o = 0;
  for (;;) { const { done, value } = await r.read(); if (done) break; const k = Math.min(value.length, expected - o); out.set(value.subarray(0, k), o); o += k; if (o >= expected) { r.cancel().catch(() => {}); break; } }
  return out;
}

/** Read the first image directory of a TIFF / BigTIFF file into { tag: number[] | string }. */
function readIFD(b) {
  const dv = view(b), le = b[0] === 0x49, magic = dv.getUint16(2, le), big = magic === 43;
  if (magic !== 42 && !big) throw new Error('not a TIFF file');
  const u64 = (p) => Number(dv.getBigUint64(p, le)), ifd = big ? u64(8) : dv.getUint32(4, le);
  if (!(ifd > 0) || ifd + (big ? 8 : 2) > b.length) throw new Error('TIFF image directory lies outside the file (truncated)');
  const n = Math.min(big ? u64(ifd) : dv.getUint16(ifd, le), 4096), esz = big ? 20 : 12, base = ifd + (big ? 8 : 2), tags = {};
  for (let i = 0; i < n && base + esz * (i + 1) <= b.length; i++) {
    const e = base + esz * i, tag = dv.getUint16(e, le), type = dv.getUint16(e + 2, le), cnt = big ? u64(e + 4) : dv.getUint32(e + 4, le), s = TYPE_SIZE[type]; if (!s || !(cnt >= 0)) continue;
    const inl = big ? 8 : 4, vp = e + (big ? 12 : 8); let o = cnt * s <= inl ? vp : big ? u64(vp) : dv.getUint32(vp, le);
    const m = Math.min(cnt, 4e6); if (o + m * s > b.length) continue;
    if (type === 2) { tags[tag] = str(b, o, o + m).replace(/\0+$/, ''); continue; }
    const out = new Array(m);
    for (let k = 0; k < m; k++, o += s) out[k] = type === 3 ? dv.getUint16(o, le) : type === 4 ? dv.getUint32(o, le) : type === 12 ? dv.getFloat64(o, le) : type === 11 ? dv.getFloat32(o, le) : type === 8 ? dv.getInt16(o, le) : type === 9 ? dv.getInt32(o, le) : type === 16 || type === 18 ? u64(o) : type === 17 ? Number(dv.getBigInt64(o, le)) : type === 5 ? dv.getUint32(o, le) / (dv.getUint32(o + 4, le) || 1) : type === 6 ? dv.getInt8(o) : dv.getUint8(o);
    tags[tag] = out;
  }
  return { tags, le, big };
}
const GEOKEY = { 1024: 'modelType', 1025: 'rasterType', 2048: 'geographicCRS', 2052: 'geogLinearUnits', 2054: 'geogAngularUnits', 3072: 'projectedCRS', 3076: 'linearUnits', 4096: 'verticalCRS', 4099: 'verticalUnits' };
const GEOCITE = { 1026: 'citation', 2049: 'geographicCitation', 3073: 'projectedCitation', 4097: 'verticalCitation' };

/**
 * GeoTIFF elevation raster → height-field surface. The raster is sampled on a regular stride so that neither side
 * exceeds opts.maxGrid samples (default 300); georeferencing is kept in meta and applied to x/y.
 */
export async function readGeoTIFF(b, ctx) {
  const { tags, le, big } = readIFD(b), one = (t, d) => (Array.isArray(tags[t]) && tags[t].length ? tags[t][0] : d);
  const W = one(256, 0), H = one(257, 0), bits = one(258, 1), comp = one(259, 1), S = one(277, 1), planar = one(284, 1), pred = one(317, 1), fmt = one(339, 1), photo = one(262, 1);
  const meta = { byteOrder: le ? 'little-endian' : 'big-endian', bigTiff: big, width: W, height: H, bitsPerSample: bits, samplesPerPixel: S, sampleFormat: { 1: 'unsigned integer', 2: 'signed integer', 3: 'floating point' }[fmt] ?? `code ${fmt}`, compression: COMPRESSION[comp] ?? `code ${comp}`, predictor: pred, tiled: !!tags[322], geoTiff: false };
  // georeferencing
  const scale = tags[33550], tie = tags[33922], xf = tags[34264], gk = tags[34735];
  if (scale) { meta.pixelScale = scale.slice(0, 3); meta.geoTiff = true; } if (tie) { meta.tiePoints = tie.slice(0, 6); meta.geoTiff = true; } if (xf && xf.length >= 16) { meta.modelTransformation = xf.slice(0, 16); meta.geoTiff = true; }
  if (gk) { meta.geoTiff = true; const keys = {}; for (let k = 4; k + 3 < gk.length; k += 4) { if (gk[k + 1] === 0 && GEOKEY[gk[k]]) keys[GEOKEY[gk[k]]] = gk[k + 3]; else if (gk[k + 1] === 34737 && GEOCITE[gk[k]] && typeof tags[34737] === 'string') keys[GEOCITE[gk[k]]] = tags[34737].slice(gk[k + 3], gk[k + 3] + gk[k + 2]).replace(/\|$/, ''); } meta.geoKeys = keys; }
  const nodata = typeof tags[42113] === 'string' && tags[42113].trim() !== '' ? Number(tags[42113]) : null; meta.noData = nodata;
  const fail = (why) => Object.assign(new Error(why), { meta });
  if (!(W > 0 && H > 0)) throw fail('TIFF image size is missing');
  if (S >= 3 || photo === 2 || photo === 3 || photo === 6) throw fail(`the raster is imagery (${S} samples per pixel, photometric code ${photo}), not a single-band elevation grid`);
  const Bs = bits / 8;
  if (![1, 2, 4, 8].includes(Bs) || (fmt === 3 && Bs < 4) || (fmt !== 3 && Bs === 8) || ![1, 2, 3].includes(fmt)) throw fail(`${bits}-bit ${meta.sampleFormat} samples are not read (supported: 8/16/32-bit integers, 32/64-bit floats)`);
  if (![1, 5, 8, 32946, 32773].includes(comp)) throw fail(`TIFF compression "${meta.compression}" is not decoded (supported: none, LZW, Deflate, PackBits)`);
  if (![1, 2, 3].includes(pred) || (pred === 3 && fmt !== 3)) throw fail(`TIFF predictor ${pred} is not supported`);
  const tiled = !!tags[322], bw = tiled ? one(322, 0) : W, bh = tiled ? one(323, 0) : Math.min(one(278, H), H), offs = tags[tiled ? 324 : 273], lens = tags[tiled ? 325 : 279];
  if (!offs || !(bw > 0 && bh > 0)) throw fail('TIFF strip/tile layout tags are missing');
  const spp = planar === 2 ? 1 : S, rowBytes = bw * spp * Bs, blockBytes = guard(rowBytes * bh, 'TIFF block byte', 1 << 30);
  const maxGrid = Math.max(2, Math.min(2000, Math.floor(ctx.opts.maxGrid ?? 300))), step = Math.max(1, Math.ceil((Math.max(W, H) - 1) / (maxGrid - 1))), nx = Math.floor((W - 1) / step) + 1, ny = Math.floor((H - 1) / step) + 1;
  const Zg = new Float64Array(guard(nx * ny, 'height-field sample', LIMITS.verts)).fill(NaN), across = tiled ? Math.ceil(W / bw) : 1, down = Math.ceil(H / bh), dv0 = new DataView(new ArrayBuffer(8));
  let missing = 0;
  for (let by = 0; by < down; by++) for (let bx = 0; bx < across; bx++) {
    const r0 = by * bh, c0 = bx * bw, r1 = Math.min(H, r0 + bh), c1 = Math.min(W, c0 + bw), fr = Math.ceil(r0 / step) * step, fc = Math.ceil(c0 / step) * step;
    if (fr >= r1 || fc >= c1) continue;                    // no sampled pixel falls in this block: it is never decoded
    const k = by * across + bx, o = offs[k], l = lens ? lens[k] : blockBytes;
    if (!(o > 0) || !(l > 0) || o + l > b.length) { missing++; continue; }
    const rawIn = b.subarray(o, o + l);
    let raw = comp === 1 ? rawIn : comp === 5 ? lzw(rawIn, blockBytes) : comp === 32773 ? packBits(rawIn, blockBytes) : await inflateZlib(rawIn, blockBytes);
    if (raw.length < rowBytes) { missing++; continue; }
    if (pred !== 1 && comp === 1) raw = raw.slice();
    const rows = Math.min(bh, Math.floor(raw.length / rowBytes)), dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength), tmp = pred === 3 ? new Uint8Array(Bs) : null;
    for (let r = fr; r < r1 && r - r0 < rows; r += step) {
      const ro = (r - r0) * rowBytes;
      if (pred === 2) { if (Bs === 1) for (let i = spp; i < bw * spp; i++) raw[ro + i] = (raw[ro + i] + raw[ro + i - spp]) & 255; else if (Bs === 2) for (let i = spp; i < bw * spp; i++) dv.setUint16(ro + 2 * i, (dv.getUint16(ro + 2 * i, le) + dv.getUint16(ro + 2 * (i - spp), le)) & 65535, le); else for (let i = spp; i < bw * spp; i++) dv.setUint32(ro + 4 * i, (dv.getUint32(ro + 4 * i, le) + dv.getUint32(ro + 4 * (i - spp), le)) >>> 0, le); }
      else if (pred === 3) for (let i = spp; i < rowBytes; i++) raw[ro + i] = (raw[ro + i] + raw[ro + i - spp]) & 255;
      for (let c = fc; c < c1; c += step) {
        const i = (c - c0) * spp; let v;
        if (pred === 3) { for (let q = 0; q < Bs; q++) tmp[q] = raw[ro + q * bw * spp + i]; for (let q = 0; q < Bs; q++) dv0.setUint8(q, tmp[q]); v = Bs === 4 ? dv0.getFloat32(0, false) : dv0.getFloat64(0, false); }
        else { const p = ro + i * Bs; v = fmt === 3 ? (Bs === 4 ? dv.getFloat32(p, le) : dv.getFloat64(p, le)) : fmt === 2 ? (Bs === 1 ? dv.getInt8(p) : Bs === 2 ? dv.getInt16(p, le) : dv.getInt32(p, le)) : Bs === 1 ? raw[p] : Bs === 2 ? dv.getUint16(p, le) : dv.getUint32(p, le); }
        Zg[(r / step) * nx + c / step] = v;
      }
    }
  }
  if (missing) ctx.warn(`${missing} strip(s)/tile(s) were missing or truncated; the corresponding area is left empty.`);
  // pixel → model coordinates
  const area = meta.geoKeys?.rasterType !== 2, half = area ? 0.5 : 0; let toXY, georef = true;
  if (meta.modelTransformation) { const m = meta.modelTransformation; toXY = (c, r) => [m[0] * (c + half) + m[1] * (r + half) + m[3], m[4] * (c + half) + m[5] * (r + half) + m[7]]; }
  else if (scale && tie && tie.length >= 6 && scale[0] > 0 && scale[1] > 0) toXY = (c, r) => [tie[3] + (c + half - tie[0]) * scale[0], tie[4] - (r + half - tie[1]) * scale[1]];
  else { georef = false; toXY = (c, r) => [c, -r]; }
  const pos = [], idx = new Int32Array(nx * ny).fill(-1), tri = []; let zmin = Infinity, zmax = -Infinity, holes = 0;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const z = Zg[j * nx + i]; if (!Number.isFinite(z) || (nodata !== null && z === nodata) || z <= -1e30) { holes++; continue; }
    const [x, y] = toXY(i * step, j * step); idx[j * nx + i] = pos.length / 3; pos.push(x, y, z); if (z < zmin) zmin = z; if (z > zmax) zmax = z;
  }
  for (let j = 0; j + 1 < ny; j++) for (let i = 0; i + 1 < nx; i++) { const a = idx[j * nx + i], c = idx[j * nx + i + 1], d = idx[(j + 1) * nx + i], e = idx[(j + 1) * nx + i + 1]; if (a < 0 || c < 0 || d < 0 || e < 0) continue; tri.push(d, e, c, d, c, a); }
  Object.assign(meta, { grid: { nx, ny, stride: step }, samplesUsed: pos.length / 3, noDataSamples: holes, zRange: pos.length ? [zmin, zmax] : null, georeferenced: georef });
  if (!pos.length) throw fail('the raster holds no valid elevation samples');
  if (step > 1) ctx.warn(`Raster ${W} × ${H} sampled every ${step} pixel(s) to a ${nx} × ${ny} height field (raise opts.maxGrid, up to 2000, for more detail).`);
  ctx.warn('Single-band raster interpreted as elevation: sample values are used as z. Check that the band really is height and in what vertical unit.');
  let units = { length: null, source: null }; const mt = meta.geoKeys?.modelType, lu = meta.geoKeys?.linearUnits, vu = meta.geoKeys?.verticalUnits;
  if (!georef) ctx.warn('The file carries no usable georeferencing (ModelPixelScale + ModelTiepoint or ModelTransformation): x/y are pixel indices, so horizontal scale is NOT physical.');
  else if (mt === 2) ctx.warn('Geographic coordinate system: x/y are in degrees of longitude/latitude while z is a length, so areas, slopes and volumes of this surface are not physical. Reproject to a projected CRS (gdalwarp) first.');
  else if (mt === 1 && (lu === 9001 || lu === 9002) && (vu === undefined || vu === lu)) units = { length: lu === 9001 ? 'm' : 'ft', source: 'file' };
  else if (mt === 1 && lu !== undefined) ctx.warn(`Projected CRS with linear unit code ${lu}${vu !== undefined && vu !== lu ? ` and vertical unit code ${vu}` : ''}: the units must be confirmed manually.`);
  return { positions: Float64Array.from(pos), triangles: Uint32Array.from(tri), kind: 'surface', units, meta };
}

/** E57: binary header plus the XML section (scan names, point counts, bounds). Point records are not decoded. */
export function readE57(b) {
  const meta = {}; if (b.length < 48) return { kind: 'metadata-only', meta };
  const dv = view(b), big = (p) => Number(dv.getBigUint64(p, true));
  Object.assign(meta, { version: `${dv.getUint32(8, true)}.${dv.getUint32(12, true)}`, fileLength: big(16), xmlOffset: big(24), xmlLength: big(32), pageSize: big(40) });
  const page = meta.pageSize;
  if (page >= 64 && page <= 1 << 20 && meta.xmlLength > 0 && meta.xmlLength < 64e6 && meta.xmlOffset < b.length) {
    // logical bytes skip the 4-byte checksum that ends every physical page
    const out = new Uint8Array(meta.xmlLength); let p = meta.xmlOffset, o = 0;
    while (o < out.length && p < b.length) { const end = (Math.floor(p / page) + 1) * page - 4, k = Math.min(end - p, out.length - o, b.length - p); if (k <= 0) { p = end + 4; continue; } out.set(b.subarray(p, p + k), o); o += k; p = end + 4; }
    const xml = str(out, 0, o), scans = [];
    for (const m of xml.matchAll(/<points\b[^>]*\brecordCount="(\d+)"/g)) scans.push({ points: +m[1] });
    const names = [...xml.matchAll(/<name\b[^>]*>(?:<!\[CDATA\[)?([^<\]]*)/g)].map((m) => m[1].trim()).filter(Boolean);
    scans.forEach((s, i) => { if (names[i]) s.name = names[i]; });
    const bound = (k) => [...xml.matchAll(new RegExp(`<${k}\\b[^>]*>([-+.\\deE]+)<`, 'g'))].map((m) => +m[1]);
    const lo = ['xMinimum', 'yMinimum', 'zMinimum'].map((k) => bound(k)), hi = ['xMaximum', 'yMaximum', 'zMaximum'].map((k) => bound(k));
    if (lo.every((a) => a.length) && hi.every((a) => a.length)) meta.cartesianBounds = { min: lo.map((a) => Math.min(...a)), max: hi.map((a) => Math.max(...a)) };
    meta.scans = scans; meta.totalPoints = scans.reduce((s, q) => s + q.points, 0); meta.xmlRead = o === out.length;
    meta.fields = [...new Set([...xml.matchAll(/<(cartesian[XYZ]|spherical\w+|intensity|color\w+|rowIndex|columnIndex)\b/g)].map((m) => m[1]))];
  }
  return { kind: 'metadata-only', meta };
}
