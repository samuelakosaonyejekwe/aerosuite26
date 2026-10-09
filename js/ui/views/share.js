// Opens a case shared as a link: #/share/<compressed case>. Nothing is uploaded; the case travels in the URL.

import { h, icon, btn, toast } from '../dom.js';
import { setCase, state } from '../../core/store.js';

const b64 = { enc: (u8) => btoa(String.fromCharCode(...u8)).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, ''), dec: (s) => Uint8Array.from(atob(s.replace(/-/g, '+').replace(/_/g, '/')), (c) => c.charCodeAt(0)) };
async function pump(stream, bytes) { const w = stream.writable.getWriter(); w.write(bytes); w.close(); return new Uint8Array(await new Response(stream.readable).arrayBuffer()); }

/** Build a share link for the current case (site weather tables and imported shapes are left out to keep it short). */
export async function shareLink() {
  const c = JSON.parse(JSON.stringify(state.case)); if (c.site) c.site.winds_aloft = []; delete c.shape; if (c.econ) delete c.econ.fx_rates;
  const raw = new TextEncoder().encode(JSON.stringify(c));
  const packed = typeof CompressionStream === 'function' ? 'z' + b64.enc(await pump(new CompressionStream('deflate-raw'), raw)) : 'j' + b64.enc(raw);
  return `${location.origin}${location.pathname}#/share/${packed}`;
}
export async function copyShareLink() {
  const url = await shareLink();
  try { if (navigator.share && /android|iphone|ipad/i.test(navigator.userAgent)) { await navigator.share({ title: `AeroSuite 26 — ${state.case.meta.name}`, url }); return; } await navigator.clipboard.writeText(url); toast('Link copied. Anyone who opens it gets this exact case.', 'ok'); }
  catch { prompt('Copy this link:', url); }
}

export async function render(root, [data], { setCrumb }) {
  setCrumb('Shared case');
  const host = h('div'); root.append(host);
  try {
    if (!data || data.length > 200000) throw new Error('The link is empty or too long.');
    const bytes = b64.dec(data.slice(1)), raw = data[0] === 'z' ? await pump(new DecompressionStream('deflate-raw'), bytes) : bytes, c = JSON.parse(new TextDecoder().decode(raw));
    if (!c || typeof c !== 'object' || !c.meta || !c.mass) throw new Error('The link does not contain an aircraft case.');
    host.append(h('div', { class: 'note' }, icon('link'), h('div', null, h('b', null, `Shared case: ${String(c.meta.name || 'Unnamed').slice(0, 80)}. `), 'Opening it replaces the case you are working on (your stored results are kept until you re-run). ', btn('Open this case', () => { setCase(c, 'custom'); toast('Shared case opened.', 'ok'); location.hash = '#/case'; }, { kind: 'primary sm' }), ' ', btn('Cancel', () => (location.hash = '#/home'), { kind: 'sm ghost' }))));
  } catch (e) { host.append(h('div', { class: 'note bad' }, icon('warn'), h('div', null, `This share link could not be read: ${e.message}`))); }
}
