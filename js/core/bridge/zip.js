// Minimal ZIP archive writer and reader (PKZIP 2.0: stored or raw-deflate entries, no ZIP64, no
// encryption). Works in browsers and in Node 20+; compression uses the platform's
// CompressionStream / DecompressionStream ('deflate-raw') and falls back to stored entries when
// that is unavailable. No third-party code.

const enc = new TextEncoder(), dec = new TextDecoder();

let TABLE = null;
/** CRC-32 (IEEE 802.3) of a byte array. */
export function crc32(bytes) {
  if (!TABLE) { TABLE = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; TABLE[n] = c >>> 0; } }
  let c = 0xffffffff;
  for (let i = 0; i < bytes.length; i++) c = TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

export const toBytes = (v) => (v instanceof Uint8Array ? v : v instanceof ArrayBuffer ? new Uint8Array(v) : enc.encode(String(v)));
export const toText = (bytes) => dec.decode(bytes);

async function pipe(bytes, stream) {
  const w = stream.writable.getWriter();
  const writing = w.write(bytes).then(() => w.close());
  const chunks = []; let n = 0;
  const r = stream.readable.getReader();
  for (;;) { const { done, value } = await r.read(); if (done) break; chunks.push(value); n += value.length; }
  await writing;
  const out = new Uint8Array(n); let o = 0;
  for (const c of chunks) { out.set(c, o); o += c.length; }
  return out;
}
const canDeflate = () => typeof CompressionStream === 'function';
const canInflate = () => typeof DecompressionStream === 'function';
export const deflateRaw = (bytes) => pipe(bytes, new CompressionStream('deflate-raw'));
export const inflateRaw = (bytes) => pipe(bytes, new DecompressionStream('deflate-raw'));

/** Reject absolute paths, drive letters and parent-directory steps; normalise separators. */
export function safePath(p) {
  const s = String(p).replace(/\\/g, '/').replace(/^\.\//, '');
  if (!s || s.startsWith('/') || /^[A-Za-z]:/.test(s) || s.split('/').some((x) => x === '..') || /[\0\r\n]/.test(s)) throw new Error(`Unsafe path in archive: "${p}"`);
  return s;
}

function dosTime(d) {
  const y = Math.max(1980, d.getFullYear());
  return { time: (d.getHours() << 11) | (d.getMinutes() << 5) | (d.getSeconds() >> 1), date: ((y - 1980) << 9) | ((d.getMonth() + 1) << 5) | d.getDate() };
}

/**
 * Build a ZIP archive.
 * @param {{path: string, text?: string, bytes?: Uint8Array, exec?: boolean}[]} files
 * @param {{deflate?: boolean, date?: Date}} [opts] deflate defaults to true when the platform supports it.
 * @returns {Promise<Uint8Array>}
 */
export async function zipWrite(files, { deflate = true, date = new Date() } = {}) {
  if (files.length > 65535) throw new Error('Too many files for a plain ZIP archive (limit 65 535).');
  const { time, date: dd } = dosTime(date), parts = [], central = [], seen = new Set();
  let offset = 0;
  for (const f of files) {
    const path = safePath(f.path);
    if (seen.has(path)) throw new Error(`Duplicate path in archive: ${path}`);
    seen.add(path);
    const name = enc.encode(path), raw = toBytes(f.bytes ?? f.text ?? ''), crc = crc32(raw);
    let data = raw, method = 0;
    if (deflate && canDeflate() && raw.length > 64) { const z = await deflateRaw(raw); if (z.length < raw.length) { data = z; method = 8; } }
    if (raw.length > 0xfffffffe || data.length > 0xfffffffe || offset > 0xfffffffe) throw new Error('Archive exceeds the 4 GB limit of a plain ZIP file.');
    const lh = new DataView(new ArrayBuffer(30));
    lh.setUint32(0, 0x04034b50, true); lh.setUint16(4, 20, true); lh.setUint16(6, 0x0800, true); lh.setUint16(8, method, true);
    lh.setUint16(10, time, true); lh.setUint16(12, dd, true); lh.setUint32(14, crc, true); lh.setUint32(18, data.length, true); lh.setUint32(22, raw.length, true);
    lh.setUint16(26, name.length, true); lh.setUint16(28, 0, true);
    parts.push(new Uint8Array(lh.buffer), name, data);
    const ch = new DataView(new ArrayBuffer(46));
    ch.setUint32(0, 0x02014b50, true); ch.setUint16(4, (3 << 8) | 20, true); ch.setUint16(6, 20, true); ch.setUint16(8, 0x0800, true); ch.setUint16(10, method, true);
    ch.setUint16(12, time, true); ch.setUint16(14, dd, true); ch.setUint32(16, crc, true); ch.setUint32(20, data.length, true); ch.setUint32(24, raw.length, true);
    ch.setUint16(28, name.length, true); ch.setUint32(38, ((0o100000 | (f.exec ? 0o755 : 0o644)) << 16) >>> 0, true); ch.setUint32(42, offset, true);
    central.push(new Uint8Array(ch.buffer), name);
    offset += 30 + name.length + data.length;
  }
  const cdSize = central.reduce((s, p) => s + p.length, 0), end = new DataView(new ArrayBuffer(22));
  end.setUint32(0, 0x06054b50, true); end.setUint16(8, files.length, true); end.setUint16(10, files.length, true); end.setUint32(12, cdSize, true); end.setUint32(16, offset, true);
  const all = [...parts, ...central, new Uint8Array(end.buffer)], out = new Uint8Array(offset + cdSize + 22);
  let o = 0; for (const p of all) { out.set(p, o); o += p.length; }
  return out;
}

/**
 * Read a ZIP archive. Directory entries are skipped; every entry's CRC is checked.
 * @param {Uint8Array|ArrayBuffer} input
 * @param {{maxBytes?: number, filter?: (name: string, size: number) => boolean}} [opts]
 * @returns {Promise<{name: string, bytes: Uint8Array}[]>}
 */
export async function zipRead(input, { maxBytes = 800e6, filter = null } = {}) {
  const b = toBytes(input), v = new DataView(b.buffer, b.byteOffset, b.byteLength);
  let e = -1;
  for (let i = b.length - 22; i >= Math.max(0, b.length - 22 - 65535); i--) if (v.getUint32(i, true) === 0x06054b50) { e = i; break; }
  if (e < 0) throw new Error('This is not a ZIP archive (no end-of-central-directory record).');
  const n = v.getUint16(e + 10, true), cdOff = v.getUint32(e + 16, true);
  if (n === 0xffff || cdOff === 0xffffffff) throw new Error('ZIP64 archives are not supported. Re-pack as a plain ZIP below 4 GB.');
  const out = []; let p = cdOff, total = 0;
  for (let k = 0; k < n; k++) {
    if (p + 46 > b.length || v.getUint32(p, true) !== 0x02014b50) throw new Error('The ZIP central directory is damaged.');
    const flags = v.getUint16(p + 8, true), method = v.getUint16(p + 10, true), crc = v.getUint32(p + 16, true), csize = v.getUint32(p + 20, true), usize = v.getUint32(p + 24, true);
    const nl = v.getUint16(p + 28, true), xl = v.getUint16(p + 30, true), cl = v.getUint16(p + 32, true), lho = v.getUint32(p + 42, true);
    const name = dec.decode(b.subarray(p + 46, p + 46 + nl));
    p += 46 + nl + xl + cl;
    if (name.endsWith('/')) continue;
    if (flags & 1) throw new Error(`Entry "${name}" is encrypted, which is not supported.`);
    if (filter && !filter(name, usize)) continue;
    total += usize; if (total > maxBytes) throw new Error('The archive expands to more data than this device should hold in memory.');
    if (lho + 30 > b.length || v.getUint32(lho, true) !== 0x04034b50) throw new Error(`Entry "${name}" has a damaged header.`);
    const start = lho + 30 + v.getUint16(lho + 26, true) + v.getUint16(lho + 28, true), comp = b.subarray(start, start + csize);
    if (comp.length !== csize) throw new Error(`Entry "${name}" is truncated.`);
    let bytes;
    if (method === 0) bytes = comp.slice();
    else if (method === 8) { if (!canInflate()) throw new Error('This browser cannot decompress ZIP entries (DecompressionStream is missing).'); bytes = await inflateRaw(comp); }
    else throw new Error(`Entry "${name}" uses compression method ${method}; only stored and deflate are supported.`);
    if (bytes.length !== usize || crc32(bytes) !== crc) throw new Error(`Entry "${name}" failed its checksum.`);
    out.push({ name: safePath(name), bytes });
  }
  return out;
}

/** Base64 of a byte array (chunked; safe for multi-megabyte inputs). */
export function base64Encode(bytes) {
  if (typeof Buffer !== 'undefined') return Buffer.from(bytes.buffer, bytes.byteOffset, bytes.byteLength).toString('base64');
  let s = ''; for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000));
  return btoa(s);
}
export function base64Decode(str) {
  if (typeof Buffer !== 'undefined') return new Uint8Array(Buffer.from(str, 'base64'));
  const s = atob(str), out = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
  return out;
}
