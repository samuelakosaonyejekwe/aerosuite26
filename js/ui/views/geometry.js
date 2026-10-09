// Geometry & mesh workbench: import from other CAD/CAE tools, recognise, inspect, heal, set units,
// measure, section, mesh sections, and pass quantities to the shared case. Not a CAD system.

import { h, setKids, icon, clear, add, num, btn, badge, card, empty, toast, pickFiles, downloadText } from '../dom.js';
import { state, patchCase, setCase, idb, emit } from '../../core/store.js';
import { createViewer } from '../viewer.js';
import { renderPlot } from '../plots.js';
import { dataTable } from '../results.js';
import { METALS } from '../../data/materials.js';
import { loadSpec } from '../../core/spec.js';
import { SUITES } from '../../core/registry.js';

const UNIT_M = { m: 1, mm: 1e-3, cm: 1e-2, in: 0.0254, ft: 0.3048 };
const SUPPORT = { native: ['Read in full', 'ok'], partial: ['Documented subset', 'warn'], metadata: ['Recognised — needs conversion', 'bad'] };
let G = null; const geo = async () => (G ||= await import('../../core/geometry/index.js'));
let selected = null;

export async function render(root, _p, { setCrumb }) {
  setCrumb('Geometry & mesh');
  const host = h('div'); root.append(host);
  const lib = await geo(), listHost = h('div', { class: 'stack' }), panel = h('div', { class: 'stack' }), busy = h('div', { class: 'row muted', hidden: true }, h('i', { class: 'spin' }), h('span', null, 'Reading file…'));
  let viewer = null; const plots = [];
  const killPlots = () => plots.splice(0).forEach((p) => p.destroy());

  async function take(files) {
    const companions = []; for (const f of files) if (/^(points|faces|boundary|owner|neighbour)$/i.test(f.name)) companions.push({ name: f.name, bytes: new Uint8Array(await f.arrayBuffer()) });
    for (const f of files) {
      if (companions.length && /^(faces|boundary|owner|neighbour)$/i.test(f.name)) continue;
      busy.hidden = false; busy.lastChild.textContent = `Reading ${f.name} (${num(f.size / 1e6, 3)} MB)…`;
      await new Promise((r) => setTimeout(r, 30));
      try {
        if (f.size > 600e6) throw new Error('Files above 600 MB are beyond what a browser tab can hold in memory. Decimate or split the model first.');
        const bytes = new Uint8Array(await f.arrayBuffer()), model = await lib.importFile(f.name, bytes, { companions });
        const entry = register(model, f.size);
        if (f.size < 60e6) idb.set('geo.' + entry.id, { name: f.name, bytes }).catch(() => {});
        selected = entry.id;
        toast(`${f.name}: ${model.formatName}, ${SUPPORT[model.support]?.[0].toLowerCase() || model.support}.`, model.support === 'metadata' ? 'bad' : 'ok');
      } catch (e) { toast(`${f.name} could not be imported: ${e.message}`, 'bad', 8000); }
    }
    busy.hidden = true; paintList(); paintPanel();
  }
  function register(model, size) {
    let analysis = null, quality = null;
    try { analysis = lib.analyse(model); } catch (e) { model.warnings.push(`Analysis failed: ${e.message}`); }
    try { quality = lib.meshQuality(model); } catch { quality = null; }
    const entry = { id: `${Date.now().toString(36)}${Math.random().toString(36).slice(2, 6)}`, model, analysis, quality, size, label: analysis?.classification?.label || model.kind };
    entry.summary = summarise(entry); state.geometry.push(entry); emit('geometry'); return entry;
  }
  const summarise = (e) => ({ name: e.model.name, format: e.model.formatName, support: e.model.support, kind: e.model.kind, units: e.model.units?.length || null, bbox: e.analysis?.bbox?.size || null, nVerts: e.analysis?.nVerts, nTris: e.analysis?.nTris, area: e.analysis?.area, volume: e.analysis?.volume, watertight: e.analysis?.watertight, classification: e.label, log: e.model.log });
  const refresh = (e, model, note) => { e.model = model; try { e.analysis = lib.analyse(model); e.quality = lib.meshQuality(model); } catch { /* keep previous */ } e.summary = summarise(e); if (note) toast(note, 'ok'); paintList(); paintPanel(); };

  function paintList() {
    clear(listHost);
    if (!state.geometry.length) { listHost.append(h('p', { class: 'muted small' }, 'Nothing imported in this session yet.')); return; }
    add(listHost, state.geometry.map((e) => h('div', { class: 'row', style: { padding: '6px 8px', borderRadius: '9px', background: e.id === selected ? 'var(--accent-soft)' : '', cursor: 'pointer' }, onclick: () => { selected = e.id; paintList(); paintPanel(); } },
      icon('cube', 18), h('span', { class: 'grow', style: { minWidth: 0 } }, h('b', null, e.model.name), h('div', { class: 'muted small' }, `${e.model.formatName} · ${e.label}`)), badge(SUPPORT[e.model.support]?.[0] || e.model.support, SUPPORT[e.model.support]?.[1]),
      h('button', { class: 'icon-btn', title: 'Remove from this session', onclick: (ev) => { ev.stopPropagation(); state.geometry.splice(state.geometry.indexOf(e), 1); idb.del('geo.' + e.id); if (selected === e.id) selected = state.geometry[0]?.id || null; paintList(); paintPanel(); } }, icon('close', 16)))));
  }

  function paintPanel() {
    killPlots(); viewer?.destroy(); viewer = null; clear(panel);
    const e = state.geometry.find((x) => x.id === selected);
    if (!e) { panel.append(empty('cube', 'Import a geometry or mesh file', 'Drop a file above. Exact CAD (STEP, IGES), tessellated surfaces (STL, OBJ, PLY, glTF, 3MF), CFD and FE meshes (Gmsh, SU2, VTK, Nastran, Abaqus, UNV, CGNS…) and point clouds are recognised. Each format is read to a documented level, shown in the capability table below.')); return; }
    const m = e.model, a = e.analysis, q = e.quality, sup = SUPPORT[m.support] || [m.support, ''];
    const fmt = lib.FORMATS.find((f) => f.id === m.format);
    // --- header notes
    if (m.support !== 'native' && fmt) panel.append(h('div', { class: `note ${m.support === 'metadata' ? 'bad' : 'warn'}` }, icon('info'), h('div', null, h('b', null, `${m.formatName}: ${sup[0].toLowerCase()}. `), fmt.limits ? `Not read: ${fmt.limits}. ` : '', fmt.pathway ? h('span', null, h('b', null, 'Conversion route: '), fmt.pathway) : null)));
    if (m.warnings?.length) panel.append(h('div', { class: 'note warn' }, icon('warn'), h('div', null, h('ul', { style: { margin: 0, paddingLeft: '18px' } }, m.warnings.slice(0, 8).map((w) => h('li', null, w))))));
    // --- units: never assumed
    const hasGeom = m.positions?.length > 0;
    if (hasGeom && !m.units?.length) {
      const sel = h('select', { class: 'inp', style: { width: 'auto' } }, Object.keys(UNIT_M).map((u) => h('option', { value: u }, u === 'm' ? 'metres' : u === 'mm' ? 'millimetres' : u === 'cm' ? 'centimetres' : u === 'in' ? 'inches' : 'feet')));
      const guess = a?.bbox ? Math.max(...a.bbox.size) : 0; if (guess > 400) sel.value = 'mm';
      panel.append(h('div', { class: 'note warn' }, icon('warn'), h('div', null, h('b', null, 'This file does not state its length unit. '), `The largest extent is ${num(guess)} in file units. What are those units? `, sel, ' ', btn('Confirm', () => { const u = sel.value; const mm = lib.transform(m, { scale: UNIT_M[u] }); mm.units = { length: 'm', source: 'user' }; mm.log.push({ step: 'Units confirmed by user', detail: `File coordinates declared as ${u}; scaled by ${UNIT_M[u]} to metres.` }); refresh(e, mm, 'Units set. The model is now in metres.'); }, { kind: 'sm primary' }), ' All measurements below are in file units until you confirm.')));
    } else if (hasGeom && m.units.length !== 'm') panel.append(h('div', { class: 'note' }, icon('info'), h('div', null, `The file states its length unit as ${m.units.length}. `, btn('Convert to metres', () => { const mm = lib.transform(m, { scale: UNIT_M[m.units.length] || 1 }); mm.units = { length: 'm', source: 'converted' }; refresh(e, mm, 'Converted to metres.'); }, { kind: 'sm primary' }))));
    const U = m.units?.length || 'file units';
    // --- viewer + summary
    const vHost = h('div');
    panel.append(h('div', { class: 'grid g2' }, vHost,
      card('What was recognised', h('div', { class: 'stack' },
        h('dl', { class: 'kv' }, h('dt', null, 'File'), h('dd', null, `${m.name} · ${num(e.size / 1e6, 3)} MB`), h('dt', null, 'Format'), h('dd', null, m.formatName, ' ', badge(sup[0], sup[1])), h('dt', null, 'Content'), h('dd', null, m.kind),
          a ? [h('dt', null, 'Looks like'), h('dd', null, classSelect(e), a.classification ? h('div', { class: 'muted small' }, `${Math.round(100 * (a.classification.confidence || 0))}% confidence — ${a.classification.why || ''}`) : null),
            h('dt', null, 'Vertices / triangles'), h('dd', null, `${num(a.nVerts)} / ${num(a.nTris)}`), h('dt', null, 'Extents x · y · z'), h('dd', null, `${num(a.bbox.size[0])} · ${num(a.bbox.size[1])} · ${num(a.bbox.size[2])} ${U}`),
            a.area ? [h('dt', null, 'Surface area'), h('dd', null, `${num(a.area)} ${U}²`)] : null, h('dt', null, 'Enclosed volume'), h('dd', null, a.volume != null ? `${num(a.volume)} ${U}³` : 'not a closed surface'),
            a.nTris ? [h('dt', null, 'Watertight'), h('dd', null, a.watertight ? badge('Yes', 'ok') : badge(`No — ${num(a.boundaryEdges)} open edges`, 'warn')), h('dt', null, 'Separate bodies'), h('dd', null, num(a.components))] : null] : null,
          m.elements?.length ? [h('dt', null, 'Elements'), h('dd', null, m.elements.map((x) => `${num(x.count)} ${x.type}`).join(', '))] : null,
          m.groups?.length ? [h('dt', null, 'Groups / zones'), h('dd', null, m.groups.slice(0, 14).map((g) => badge(`${g.name}${g.count ? ` (${g.count})` : ''}`)), m.groups.length > 14 ? ` +${m.groups.length - 14} more` : '')] : null),
        metaBlock(m.meta)))));
    if (hasGeom) { try { viewer = createViewer(vHost); viewer.setModel(m); } catch (err) { vHost.append(h('div', { class: 'note warn' }, `3-D preview unavailable: ${err.message}`)); } } else vHost.append(empty('info', 'No displayable geometry', 'Only header information could be read from this file. Follow the conversion route above to bring in the shape itself.'));
    if (!hasGeom) { panel.append(provenance(m)); return; }

    // --- quality
    if (q?.perType?.length) {
      const pg = h('div', { class: 'plots' }), gradeKind = { good: 'ok', acceptable: 'ok', poor: 'warn', unusable: 'bad' }[q.grade] || '';
      panel.append(card('Mesh and tessellation quality', h('div', { class: 'stack' },
        h('div', { class: `note ${gradeKind || 'warn'}` }, icon(gradeKind === 'ok' ? 'check' : 'warn'), h('div', null, h('b', null, `Overall: ${q.grade}. `), q.issues?.length ? q.issues.join(' ') : 'No solver-relevant defects found.')),
        dataTable({ title: '', columns: ['Element type', 'Count', 'Worst aspect ratio', 'Mean aspect ratio', 'Worst skewness', 'Smallest angle [°]', 'Largest angle [°]', 'Inverted'], rows: q.perType.map((t) => [t.type, t.count, t.aspect?.max, t.aspect?.mean, t.skew?.max, t.minAngle?.min, t.maxAngle?.max, t.negJacobian]) }), pg), { collapsible: true }));
      const t0 = q.perType[0];
      if (t0?.skew?.hist) plots.push(renderPlot(pg, { type: 'bar', title: `Skewness distribution (${t0.type})`, ylabel: 'Elements', categories: t0.skew.hist.centers.map((c) => num(c, 2)), series: [{ name: 'Elements', y: t0.skew.hist.counts }] }));
      if (t0?.aspect?.hist) plots.push(renderPlot(pg, { type: 'bar', title: `Aspect-ratio distribution (${t0.type})`, ylabel: 'Elements', categories: t0.aspect.hist.centers.map((c) => num(c, 3)), series: [{ name: 'Elements', y: t0.aspect.hist.counts }] }));
    }
    // --- healing
    if (a?.nTris) {
      const tol = h('input', { class: 'inp', type: 'number', step: 'any', min: 0, placeholder: 'automatic', style: { width: '130px' }, 'aria-label': 'Weld tolerance' }), o = { deg: h('input', { type: 'checkbox', checked: true }), ori: h('input', { type: 'checkbox', checked: true }), unref: h('input', { type: 'checkbox', checked: true }) }, logHost = h('div');
      panel.append(card('Inspect and heal', h('div', { class: 'stack' },
        h('div', { class: 'kpis' }, kp('Open (boundary) edges', a.boundaryEdges, a.boundaryEdges ? 'warn' : 'ok'), kp('Non-manifold edges', a.nonManifoldEdges, a.nonManifoldEdges ? 'bad' : 'ok'), kp('Duplicate vertices', a.duplicateVerts, a.duplicateVerts ? 'warn' : 'ok'), kp('Degenerate triangles', a.degenerateTris, a.degenerateTris ? 'warn' : 'ok'), kp('Orientation conflicts', a.orientationConflicts, a.orientationConflicts ? 'warn' : 'ok')),
        h('p', { class: 'muted small' }, 'Healing never runs on its own. Each operation below is explicit, reports exactly what it changed and how far any vertex moved, and is written to the provenance record. Engineering dimensions are not altered beyond the weld tolerance you set.'),
        h('div', { class: 'row small' }, h('label', { class: 'row' }, 'Weld vertices closer than', tol, U), h('label', { class: 'row' }, o.deg, 'Remove zero-area triangles'), h('label', { class: 'row' }, o.ori, 'Make facet orientation consistent'), h('label', { class: 'row' }, o.unref, 'Drop unused vertices')),
        h('div', { class: 'row' }, btn('Heal', () => { try { const r = lib.heal(m, { weldTol: tol.value === '' ? null : Number(tol.value), removeDegenerate: o.deg.checked, fixOrientation: o.ori.checked, removeUnreferenced: o.unref.checked }); clear(logHost); logHost.append(dataTable({ title: 'Healing log', columns: ['Step', 'Result'], rows: r.log.map((l) => [l.step, l.detail]) })); if (r.changed) refresh(e, r.model, 'Healing applied and logged.'); else toast('Nothing needed changing.', 'info'); } catch (err) { toast(`Healing failed: ${err.message}`, 'bad'); } }, { ic: 'check', kind: 'primary' }),
          btn('Swap Y and Z axes', () => refresh(e, lib.transform(m, { swapYZ: true }), 'Axes swapped (logged).'), { kind: 'ghost' }), h('span', { class: 'muted small' }, 'Suites assume x aft, y spanwise, z up.')), logHost), { collapsible: true, open: !a.watertight }));
    }
    // --- measure and use in the case
    const metals = Object.keys(METALS), matSel = h('select', { class: 'inp', style: { width: 'auto' } }, metals.map((k) => h('option', { value: k, selected: k === state.case.struct.material }, k))), massHost = h('div');
    const mass = () => { try { const r = lib.massProperties(m, METALS[matSel.value].rho); clear(massHost); massHost.append(r.closed ? h('dl', { class: 'kv' }, h('dt', null, 'Solid volume'), h('dd', null, `${num(r.volume)} ${U}³`), h('dt', null, 'Mass if solid'), h('dd', null, `${num(r.mass)} kg`), h('dt', null, 'Centre of gravity'), h('dd', null, r.cg.map((v) => num(v)).join(', ')), h('dt', null, 'Inertia Ixx, Iyy, Izz'), h('dd', null, `${num(r.inertia[0][0])}, ${num(r.inertia[1][1])}, ${num(r.inertia[2][2])} kg·m²`)) : h('p', { class: 'muted small' }, 'Mass properties need a closed surface. Heal the model first, or import a watertight tessellation.')); } catch (err) { setKids(massHost, h('p', { class: 'muted small' }, err.message)); } };
    const d = a?.bbox?.size || [0, 0, 0], metric = m.units?.length === 'm';
    const apply = (label, sec, key, val) => btn(`${label}: ${num(val)} m`, () => { patchCase(sec, key, Number(val.toPrecision(6))); toast(`${label} written to the case.`, 'ok'); }, { kind: 'sm', ic: 'link', disabled: !metric || !(val > 0), title: metric ? 'Write this dimension into the shared case' : 'Confirm units first' });
    panel.append(card('Measure and send to the case', h('div', { class: 'stack' },
      h('p', { class: 'muted small' }, metric ? 'Dimensions taken from the model can be written straight into the shared case, where every suite picks them up.' : 'Confirm the length unit above to enable these.'),
      h('div', { class: 'row' }, apply('Wing span', 'wing', 'b_m', d[1]), apply('Fuselage length', 'fuselage', 'len_m', d[0]), apply('Fuselage diameter', 'fuselage', 'dia_m', Math.max(d[2], 0)), apply('Rotor radius', 'rotor', 'R_m', Math.max(d[0], d[1]) / 2), apply('Propeller diameter', 'prop', 'prop_dia_m', Math.max(d[1], d[2]))),
      a?.nTris && lib.simplifyForSolver ? h('div', { class: 'row' }, btn('Fly this shape in the 3-D flow solver', () => { try { const sh = lib.simplifyForSolver(m, 3000); state.case.shape = { ...sh, name: m.name, units: 'm' }; setCase(state.case, state.preset); toast('Shape sent. Open Suite 1 → 3-D RANS to run it.', 'ok'); location.hash = '#/suite/cfd/rans3d'; } catch (err) { toast(`Could not prepare the shape: ${err.message}`, 'bad'); } }, { ic: 'flow', kind: 'primary sm', disabled: !metric, title: metric ? 'Simplify the surface and hand it to the immersed-boundary Navier–Stokes solver' : 'Confirm units first' }), state.case.shape ? btn('Remove the shape from the case', () => { delete state.case.shape; setCase(state.case, state.preset); toast('Shape removed; the solver returns to the parametric wing.', 'info'); paintPanel(); }, { kind: 'sm ghost' }) : null) : null,
      h('div', { class: 'row small' }, 'Material for mass properties', matSel, btn('Compute mass properties', mass, { kind: 'sm' })), massHost), { collapsible: true }));

    // --- sectioning, meshing and section properties
    if (a?.nTris) {
      const axis = h('select', { class: 'inp', style: { width: 'auto' } }, ['y', 'x', 'z'].map((x) => h('option', { value: x }, `${x} = constant`))), posPct = h('input', { type: 'range', min: 1, max: 99, value: 30, style: { flex: 1, minWidth: '140px' }, 'aria-label': 'Section position' }), posLab = h('span', { class: 'small mono' }), secHost = h('div', { class: 'stack' });
      const value = () => { const k = 'xyz'.indexOf(axis.value); return a.bbox.min[k] + (Number(posPct.value) / 100) * a.bbox.size[k]; };
      const lab = () => (posLab.textContent = `${axis.value} = ${num(value())} ${U}`); lab(); posPct.addEventListener('input', lab); axis.addEventListener('change', lab);
      const cut = () => {
        killPlots(); clear(secHost);
        try {
          const loops = lib.slice(m, { axis: axis.value, value: value() }).filter((l) => l.length > 2);
          if (!loops.length) { secHost.append(h('p', { class: 'muted' }, 'The plane does not cut the model here.')); return; }
          const main = loops.slice().sort((p, r) => r.length - p.length)[0], sm = lib.sectionMetrics(main), pg = h('div', { class: 'plots' });
          secHost.append(h('div', { class: 'kpis' }, kp('Section loops', loops.length), kp(`Chord / length [${U}]`, sm.chord), kp(`Thickness [${U}]`, sm.thickness), kp('Thickness ratio', sm.tc), kp(`Enclosed area [${U}²]`, sm.area), kp(`Perimeter [${U}]`, sm.perimeter)), pg);
          plots.push(renderPlot(pg, { type: 'line', title: 'Section outline', xlabel: `In-plane u [${U}]`, ylabel: `In-plane v [${U}]`, equalAspect: true, series: loops.slice(0, 6).map((l, i) => ({ name: `Loop ${i + 1}`, x: [...l.map((p) => p[0]), l[0][0]], y: [...l.map((p) => p[1]), l[0][1]] })) }));
          const hIn = h('input', { class: 'inp', type: 'number', step: 'any', value: Number((sm.chord / 40).toPrecision(2)), style: { width: '110px' }, 'aria-label': 'Target element size' }), propHost = h('div', { class: 'stack' });
          const solve = (study) => { clear(propHost); try {
            const hh = Number(hIn.value); if (!(hh > 0)) throw new Error('Enter a positive element size.');
            const mesh = lib.meshSection(loops, { h: hh }), pr = lib.sectionProperties(mesh), g2 = h('div', { class: 'plots' });
            propHost.append(h('div', { class: 'kpis' }, kp('Triangles', mesh.quality.nTris), kp('Smallest angle [°]', mesh.quality.minAngle, mesh.quality.minAngle > 20 ? 'ok' : 'warn'), kp(`Area A [${U}²]`, pr.A), kp(`Ixx [${U}⁴]`, pr.Ixx), kp(`Iyy [${U}⁴]`, pr.Iyy), kp(`Torsion constant J [${U}⁴]`, pr.J)), pr.note ? h('p', { class: 'muted small' }, pr.note) : null, g2);
            plots.push(renderPlot(g2, { type: 'tri', title: 'Section mesh (distance from centroid)', xlabel: `u [${U}]`, ylabel: `v [${U}]`, zlabel: `r [${U}]`, equalAspect: true, edges: true, nodes: mesh.nodes, tris: mesh.tris, values: mesh.nodes.map((p) => Math.hypot(p[0] - pr.cx, p[1] - pr.cy)) }));
            if (study) { const s = lib.sectionConvergence(loops, [hh / 2, hh, hh * 2]); propHost.append(h('div', { class: 'note' }, icon('info'), h('div', null, h('b', null, 'Mesh sensitivity of J: '), `observed order ${num(s.gci?.p, 3)}, grid convergence index ${num(100 * (s.gci?.gciFine ?? NaN), 3)}% at the finest level, extrapolated J = ${num(s.gci?.fExact)}.`)), dataTable({ title: 'Refinement levels', columns: [`h [${U}]`, 'Triangles', `J [${U}⁴]`, `A [${U}²]`], rows: s.rows.map((r) => [r.h, r.nTris, r.J, r.A]) })); }
            if (metric) { const E = METALS[matSel.value]; propHost.append(h('p', { class: 'small' }, `As a solid ${matSel.value} section: bending stiffness E·Ixx = ${num(E.E * pr.Ixx)} N·m², torsional stiffness G·J = ${num(E.G * pr.J)} N·m², mass per length = ${num(E.rho * pr.A)} kg/m. `), h('div', { class: 'row' }, btn(`Use thickness ratio ${num(sm.tc, 3)} for the wing`, () => { patchCase('wing', 'tc', Number(sm.tc.toPrecision(4))); toast('Wing thickness ratio written to the case.', 'ok'); }, { kind: 'sm', ic: 'link', disabled: !(sm.tc > 0 && sm.tc < 0.5) }))); }
          } catch (err) { propHost.append(h('div', { class: 'note bad' }, icon('warn'), h('div', null, `Section meshing failed: ${err.message}`))); } };
          secHost.append(h('div', { class: 'row small' }, h('label', { class: 'row' }, 'Element size', hIn, U), btn('Mesh and compute section properties', () => solve(false), { kind: 'primary sm', ic: 'play' }), btn('…with mesh-sensitivity study', () => solve(true), { kind: 'sm' })), propHost);
        } catch (err) { secHost.append(h('div', { class: 'note bad' }, icon('warn'), h('div', null, `Sectioning failed: ${err.message}`))); }
      };
      panel.append(card('Cut a section, mesh it and get section properties', h('div', { class: 'stack' }, h('p', { class: 'muted small' }, 'A plane cut gives the section outline (for a wing or blade: chord and thickness ratio). Meshing the outline with triangles and solving the Saint-Venant torsion problem by finite elements gives area, second moments and the torsion constant, with a mesh-sensitivity check.'), h('div', { class: 'row' }, axis, posPct, posLab, btn('Cut section', cut, { ic: 'play', kind: 'primary' })), secHost), { collapsible: true }));
    }
    // --- export + provenance
    panel.append(card('Export for other tools', h('div', { class: 'row' }, [['stl', 'STL'], ['obj', 'OBJ'], ['vtk', 'VTK'], ['msh', 'Gmsh MSH'], ['su2', 'SU2'], ['json', 'JSON']].map(([id, label]) => btn(label, () => { try { downloadText(`${m.name.replace(/\.[^.]+$/, '')}.${id}`, lib.exportModel(m, id)); } catch (err) { toast(`Export failed: ${err.message}`, 'bad'); } }, { kind: 'sm', ic: 'download' }))), { collapsible: true, open: false }), provenance(m));
  }
  const kp = (label, value, st = '') => h('div', { class: `kpi ${st}` }, h('span', { class: 'l' }, label), h('div', { class: 'v' }, num(value)));
  const provenance = (m) => card('Provenance record', dataTable({ title: '', columns: ['#', 'Step', 'Detail'], rows: (m.log || []).map((l, i) => [i + 1, l.step, typeof l.detail === 'string' ? l.detail : JSON.stringify(l.detail)]) }), { collapsible: true, open: false });
  const CLASSES = ['complete aircraft', 'wing / lifting surface', 'rotor or propeller blade', 'fuselage / body of revolution', 'landing gear', 'engine / nacelle', 'structural mesh', 'CFD volume mesh', 'point cloud / scan', 'generic part'];
  function classSelect(e) { const cur = e.label; return h('select', { class: 'inp', style: { width: 'auto' }, title: 'Correct the classification if it is wrong', onchange: (ev) => { e.label = ev.target.value; e.model.log.push({ step: 'Classification set by user', detail: ev.target.value }); e.summary = summarise(e); paintList(); } }, [...new Set([cur, ...CLASSES])].map((c) => h('option', { value: c, selected: c === cur }, c))); }
  function metaBlock(meta) { const ent = Object.entries(meta || {}).filter(([, v]) => v != null && v !== '' && (typeof v !== 'object' || Object.keys(v).length)); if (!ent.length) return null; return h('details', null, h('summary', { class: 'small', style: { cursor: 'pointer' } }, 'File header and metadata'), h('dl', { class: 'kv small', style: { marginTop: '6px' } }, ent.slice(0, 40).map(([k, v]) => [h('dt', null, k), h('dd', null, typeof v === 'object' ? JSON.stringify(v).slice(0, 400) : String(v).slice(0, 400))]))); }

  const drop = h('div', { class: 'drop', tabIndex: 0, role: 'button', onclick: async () => take(await pickFiles({ multiple: true })), onkeydown: async (ev) => { if (ev.key === 'Enter') take(await pickFiles({ multiple: true })); }, ondragover: (ev) => { ev.preventDefault(); drop.classList.add('over'); }, ondragleave: () => drop.classList.remove('over'), ondrop: (ev) => { ev.preventDefault(); drop.classList.remove('over'); take([...ev.dataTransfer.files]); } },
    icon('upload', 32), h('b', null, 'Drop geometry or mesh files here, or click to choose'), h('span', { class: 'muted small' }, 'STEP · IGES · Parasolid · STL · OBJ · PLY · 3MF · glTF · Gmsh · SU2 · VTK · CGNS · Nastran · Abaqus · ANSYS · UNV · Tecplot · Plot3D · OpenFOAM · LAS · and more'), busy);
  const fmtHost = h('div'), suiteHost = h('div');
  host.append(
    h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Step 2 · import, never draw'), h('h1', null, 'Geometry & mesh'), h('p', null, 'Bring in geometry and meshes made in your CAD, meshing or scanning tools. The workbench tells you exactly what it understood, checks quality, heals defects on request with a full record, confirms units instead of guessing, and turns shapes into the numbers the suites need.'))),
    h('div', { class: 'grid', style: { gridTemplateColumns: 'minmax(0, 2fr) minmax(240px, 1fr)' } }, drop, card('Imported in this session', listHost)),
    h('div', { class: 'gap' }), panel, h('div', { class: 'gap' }),
    card('Format capability profiles', fmtHost, { collapsible: true, open: false }), h('div', { class: 'gap' }),
    card('What each suite expects', suiteHost, { collapsible: true, open: false }));
  fmtHost.append(h('p', { class: 'muted small' }, 'A file extension proves nothing, so every reader checks the content itself and is tested individually. “Read in full” formats are parsed completely; “documented subset” formats are read to the level stated; “recognised” formats are identified from their signature and need the conversion route shown.'),
    dataTable({ title: '', columns: ['Format', 'Extensions', 'Type', 'Support', 'What is read', 'What is preserved', 'Not read', 'Conversion route'], rows: lib.FORMATS.map((f) => [f.name, f.ext.map((x) => '.' + x).join(' '), f.category, SUPPORT[f.support]?.[0] || f.support, f.reads || '', f.preserves || '', f.limits || '', f.pathway || '']) }, { max: 200 }));
  loadSpec().then((sp) => suiteHost.append(h('div', { class: 'stack spec' }, sp.platform.geometry.filter((p) => p.includes(':')).slice(0, 10).map((p) => h('p', null, h('b', null, p.split(':')[0] + ': '), p.slice(p.indexOf(':') + 1))), h('h4', null, 'By suite'), SUITES.map((s) => h('details', null, h('summary', { style: { cursor: 'pointer' } }, `${s.d}. ${s.short}`), h('p', { style: { marginTop: '6px' } }, sp.suites[s.n].geometry)))))).catch(() => {});

  // restore files stored on this device from an earlier session
  if (!state.geometry.length) { try { for (const k of await idb.keys('geo.')) { const r = await idb.get(k); if (!r?.bytes) continue; try { const e = register(await lib.importFile(r.name, r.bytes), r.bytes.length); idb.del(k); idb.set('geo.' + e.id, r); } catch { idb.del(k); } } } catch { /* storage unavailable */ } }
  if (!selected || !state.geometry.some((x) => x.id === selected)) selected = state.geometry[0]?.id || null;
  paintList(); paintPanel();
  return () => { killPlots(); viewer?.destroy(); };
}
