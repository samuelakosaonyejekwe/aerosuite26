// Parasolid XT (transmit) reader: text (.x_t) and binary (.x_b) node streams → boundary-representation graph →
// in-memory STEP AP214 for the OpenCASCADE kernel. Written from the published "Parasolid XT Format Reference":
// the file is a schema-driven sequence of nodes; files from Parasolid V14 on embed their schema as edits against
// the V13 base schema (SCH_13006), which is tabulated below. Lengths are metres (Parasolid's modelling unit).

import { str } from './parsers-util.js';
import { V, StepWriter, expandKnots, nurbsCurvePoint, nurbsSurfacePoint } from './parsers-brep.js';

// ---- base schema: "<type> <NAME> field:code[:n] …"; code = d int, n short, u byte, c char, l logical, w wide char,
// f double, p pointer, v vector, i interval, b box, h intersection point; n = fixed array length, * = variable length.
const GEO = 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c';
const BASE_TEXT = `
10 ASSEMBLY highest_node_id:d attributes_groups:p attribute_chains:p list:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p type:u sub_instance:p
11 INSTANCE node_id:d attributes_groups:p type:u part:p transform:p assembly:p next_in_part:p prev_in_part:p next_of_part:p prev_of_part:p
12 BODY highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p
13 SHELL node_id:d attributes_groups:p body:p next:p face:p edge:p vertex:p region:p front_face:p
14 FACE node_id:d attributes_groups:p tolerance:f next:p previous:p loop:p shell:p surface:p sense:c next_on_surface:p previous_on_surface:p next_front:p previous_front:p front_shell:p
15 LOOP node_id:d attributes_groups:p fin:p face:p next:p
16 EDGE node_id:d attributes_groups:p tolerance:f fin:p previous:p next:p curve:p next_on_curve:p previous_on_curve:p owner:p
17 FIN attributes_groups:p loop:p forward:p backward:p vertex:p other:p edge:p curve:p next_at_vx:p sense:c
18 VERTEX node_id:d attributes_groups:p fin:p previous:p next:p point:p tolerance:f owner:p
19 REGION node_id:d attributes_groups:p body:p next:p previous:p shell:p type:c
29 POINT node_id:d attributes_groups:p owner:p next:p previous:p pvec:v
30 LINE ${GEO} pvec:v direction:v
31 CIRCLE ${GEO} centre:v normal:v x_axis:v radius:f
32 ELLIPSE ${GEO} centre:v normal:v x_axis:v major_radius:f minor_radius:f
38 INTERSECTION ${GEO} surface:p:2 chart:p start:p end:p
40 CHART base_parameter:f base_scale:f chart_count:d chordal_error:f angular_error:f parameter_error:f:2 hvec:h:*
41 LIMIT type:c hvec:h:*
45 BSPLINE_VERTICES vertices:f:*
50 PLANE ${GEO} pvec:v normal:v x_axis:v
51 CYLINDER ${GEO} pvec:v axis:v radius:f x_axis:v
52 CONE ${GEO} pvec:v axis:v radius:f sin_half_angle:f cos_half_angle:f x_axis:v
53 SPHERE ${GEO} centre:v radius:f axis:v x_axis:v
54 TORUS ${GEO} centre:v axis:v major_radius:f minor_radius:f x_axis:v
56 BLENDED_EDGE ${GEO} blend_type:c surface:p:2 spine:p range:f:2 thumb_weight:f:2 boundary:p:2 start:p end:p
59 BLEND_BOUND ${GEO} boundary:n blend:p
60 OFFSET_SURF ${GEO} check:c true_offset:l surface:p offset:f scale:f
67 SWEPT_SURF ${GEO} section:p sweep:v scale:f
68 SPUN_SURF ${GEO} profile:p base:v axis:v start:v end:v start_param:f end_param:f x_axis:v scale:f
70 LIST node_id:d owner:p next:p previous:p list_type:d list_length:d block_length:d size_of_entry:d list_block:p finger_block:p finger_index:d notransmit:l
74 POINTER_LIS_BLOCK n_entries:d next_block:p entries:p:*
79 ATT_DEF_ID string:c:*
80 ATTRIB_DEF next:p identifier:p type_id:d actions:u:8 field_names:p legal_owners:l:14 fields:u:*
81 ATTRIBUTE node_id:d definition:p owner:p next:p previous:p next_of_type:p previous_of_type:p fields:p:*
82 INT_VALUES values:d:*
83 REAL_VALUES values:f:*
84 CHAR_VALUES values:c:*
85 POINT_VALUES values:v:*
86 VECTOR_VALUES values:v:*
87 AXIS_VALUES values:v:*
88 TAG_VALUES values:d:*
89 DIRECTION_VALUES values:v:*
90 GROUP node_id:d attributes_groups:p owner:p next:p previous:p type:u first_member:p
91 MEMBER_OF_GROUP dummy_node_id:d owning_group:p owner:p next:p previous:p next_member:p previous_member:p
98 UNICODE_VALUES values:w:*
99 FIELD_NAMES names:p:*
100 TRANSFORM node_id:d owner:p next:p previous:p rotation_matrix:f:9 translation_vector:v scale:f flag:d perspective_vector:v
101 WORLD assembly:p attribute:p body:p transform:p surface:p curve:p point:p alive:l attrib_def:p highest_id:d current_id:d
102 KEY string:c:*
120 PE_SURF ${GEO} type:c data:p tf:p internal_geom:p:*
121 INT_PE_DATA geom_type:d real_array:p int_array:p
122 EXT_PE_DATA key:p real_array:p int_array:p
124 B_SURFACE ${GEO} nurbs:p data:p
125 SURFACE_DATA original_uint:i original_vint:i extended_uint:i extended_vint:i self_int:u original_u_start:c original_u_end:c original_v_start:c original_v_end:c extended_u_start:c extended_u_end:c extended_v_start:c extended_v_end:c analytic_form_type:c swept_form_type:c spun_form_type:c blend_form_type:c analytic_form:p swept_form:p spun_form:p blend_form:p
126 NURBS_SURF u_periodic:l v_periodic:l u_degree:n v_degree:n n_u_vertices:d n_v_vertices:d u_knot_type:u v_knot_type:u n_u_knots:d n_v_knots:d rational:l u_closed:l v_closed:l surface_form:u vertex_dim:n bspline_vertices:p u_knot_mult:p v_knot_mult:p u_knots:p v_knots:p
127 KNOT_MULT mult:n:*
128 KNOT_SET knots:f:*
130 PE_CURVE ${GEO} type:c data:p tf:p internal_geom:p:*
133 TRIMMED_CURVE ${GEO} basis_curve:p point_1:v point_2:v parm_1:f parm_2:f
134 B_CURVE ${GEO} nurbs:p data:p
135 CURVE_DATA self_int:u analytic_form:p
136 NURBS_CURVE degree:n n_vertices:d vertex_dim:n n_knots:d knot_type:u periodic:l closed:l rational:l curve_form:u bspline_vertices:p knot_mult:p knots:p
137 SP_CURVE ${GEO} surface:p b_curve:p original:p tolerance_to_original:f
141 GEOMETRIC_OWNER owner:p next:p previous:p shared_geometry:p`;
/** The tabulated V13 base schema: Map(node type → { name, fields: [{ name, code, n }] }). */
export function xtBaseSchema() { return baseSchema(); }
function baseSchema() {
  const m = new Map();
  for (const ln of BASE_TEXT.trim().split('\n')) { const t = ln.trim().split(/\s+/); m.set(+t[0], { name: t[1], fields: t.slice(2).map((f) => { const [name, code, n] = f.split(':'); return { name, code, n: n === '*' ? -1 : n ? +n : 0 }; }) }); }
  return m;
}
/** Adjustments for schemas older than the V13 base, as far as they could be checked against real files. */
function legacySchema(S, version) {
  if (version < 11000) { const fn = S.get(17); fn.fields = fn.fields.filter((f) => f.name !== 'attributes_groups'); }
  if (version < 10000) { for (const [t, names] of [[40, ['chart_count']], [136, ['n_vertices', 'n_knots']], [126, ['n_u_vertices', 'n_v_vertices', 'n_u_knots', 'n_v_knots']]]) { const d = S.get(t); d.fields = d.fields.map((f) => (names.includes(f.name) ? { ...f, code: 'n' } : f)); } }      // counts are shorts in V9 binary files (seen for CHART and NURBS_CURVE; assumed alike for NURBS_SURF)
  if (version < 10000) { const bd = S.get(12); bd.fields = bd.fields.filter((f) => f.name !== 'nom_geom_state'); }
  if (version < 13000) { const ad = S.get(80); ad.fields = ad.fields.filter((f) => f.name !== 'field_names').map((f) => (f.name === 'legal_owners' ? { ...f, n: 13 } : f)); }
  if (version < 13000) { for (const t of [125, 135]) { const d = S.get(t); d.fields = d.fields.filter((f) => !/_form$|_form_type$/.test(f.name)); } }
  return S;
}

// ---- learned schemas ----
// Files that do not embed their schema need the layouts of their schema version. The vendor's schema files are not
// public, but every file that does embed its schema states the layouts of the node types it uses. The table below
// holds every layout that differs between versions, as observed in real files: "<type> <NAME> <versions seen> <fields>"
// ("=" stands for the V13 base layout). For a version that is not listed, the layout is inferred (see inferLayouts).
const LEARNED_TEXT = `
12 BODY 35102,36001 highest_node_id:d attributes_groups:p attribute_chains:p lattice:p surface:p curve:p point:p mesh:p polyline:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_lattice:p boundary_surface:p boundary_curve:p boundary_point:p boundary_mesh:p boundary_polyline:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d mesh_offset_data:p
12 BODY 25001 highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d
12 BODY 31001,31100 highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p mesh:p polyline:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p boundary_mesh:p boundary_polyline:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d mesh_offset_data:p
12 BODY 26105,28002 highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p boundary_mesh:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d
12 BODY 19008,20000 highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p
12 BODY 16100 =
19 REGION 35102,36001 node_id:d attributes_groups:p body:p next:p previous:p shell:p frame:p type:c owner:p
19 REGION 25001,26105,28002,31001,31100 node_id:d attributes_groups:p body:p next:p previous:p shell:p type:c owner:p
19 REGION 16100,19008,20000 =
38 INTERSECTION 19008 =
38 INTERSECTION 36001 node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c surface:p:2 chart:p start:p end:p intersection_data:p
41 LIMIT 19008 =
41 LIMIT 36001 type:c term_use:c hvec:h:*
70 LIST 16100,20000,25001,26105,28002,31001,31100,35102,36001 node_id:d list_type:u notransmit:l owner:p next:p previous:p list_length:d block_length:d finger_index:d finger_block:p list_block:p
74 POINTER_LIS_BLOCK 20000,25001,26105,28002,31001,31100,35102,36001 n_entries:d index_map_offset:d next_block:p entries:p:*
74 POINTER_LIS_BLOCK 16100 =
80 ATTRIB_DEF 35102 next:p identifier:p type_id:d actions:u:8 field_names:p legal_owners:l:16 fields:u:*
80 ATTRIB_DEF 16100,20000,25001,26105,28002,31001,31100 =
80 ATTRIB_DEF 36001 next:p identifier:p type_id:d actions:u:8 field_names:p legal_owners:l:17 fields:u:*
204 INTERSECTION_DATA 36001 uv_type:u values:f:*`;
// node types whose layout equalled the V13 base in every embedded schema seen (versions 16100 to 36001): taken as fixed
const STABLE_TYPES = new Set([13, 14, 15, 16, 17, 18, 29, 30, 31, 32, 40, 45, 50, 51, 52, 53, 54, 56, 59, 79, 81, 82, 83, 84, 87, 90, 91, 98, 99, 124, 125, 126, 127, 128, 133, 134, 135, 136, 137, 141]);
let learnedCache = null;
function learnedLayouts() {
  if (learnedCache) return learnedCache;
  const base = baseSchema(), m = new Map();
  for (const ln of LEARNED_TEXT.trim().split('\n')) {
    const t = ln.trim().split(/\s+/), type = +t[0], vers = t[2].split(',').map(Number);
    const fields = t[3] === '=' ? base.get(type).fields : t.slice(3).map((f) => { const [name, code, n] = f.split(':'); return { name, code, n: n === '*' ? -1 : n ? +n : 0 }; });
    if (!m.has(type)) m.set(type, []); m.get(type).push({ name: t[1], fields, vers, base: t[3] === '=' });
  }
  return (learnedCache = m);
}
/** Schema versions whose layouts were learned from files that embed them. */
export function xtLearnedVersions() { return [...new Set([...learnedLayouts().values()].flat().flatMap((c) => c.vers))].sort((a, b) => a - b); }
/** Candidate layouts of a node type for a schema version, nearest observed version first. */
function candidateLayouts(type, version, generic = null) {
  const base = legacySchema(baseSchema(), Math.min(version, 13006)).get(type), learned = learnedLayouts().get(type) || [], out = [];
  const dist = (vers) => Math.min(...vers.map((v) => Math.abs(v - version) + (v > version ? 0.5 : 0)));
  for (const c of learned) out.push({ name: c.name, fields: c.fields, d: dist(c.base ? [...c.vers, 13006] : c.vers), from: c.vers });
  if (base && (version < 13006 || !learned.some((c) => c.base))) out.push({ name: base.name, fields: base.fields, d: version < 13006 ? -1 : learned.length ? dist([13006]) : 0, from: [13006] });
  if (version < 13006) { const sg = (c) => c.fields.map((f) => f.name + f.code + f.n).join(); for (const v of [13006, 12000, 10000]) { if (v <= version) continue; const alt = legacySchema(baseSchema(), v).get(type); if (alt && !out.some((c) => sg(c) === sg(alt))) out.push({ name: alt.name, fields: alt.fields, d: 5e5 + v, from: [v] }); } }     // old schemas written by newer modellers follow later layouts
  if (type === 80 && out.length) {      // the legal-owner table grows with the number of owner classes: also try the lengths in between
    const tpl = out.slice().sort((x, y) => x.d - y.d)[0], have = new Set(out.map((c) => c.fields.find((f) => f.name === 'legal_owners')?.n));
    for (const t2 of version < 13006 ? out.filter((c) => c.d < 6e5) : [tpl]) for (let n = 13; n <= 22; n++) { const fields = t2.fields.map((f) => (f.name === 'legal_owners' ? { ...f, n } : f)), k = fields.map((f) => f.name + f.code + f.n).join(); if (!out.some((c) => c.fields.map((f) => f.name + f.code + f.n).join() === k)) out.push({ name: tpl.name, fields, d: 1e6 + n, from: [], rank: 1 }); }
  }
  for (const c of out) c.rank = c.from.length && Math.min(...c.from, c.base !== false && c.from.includes(13006) ? 13006 : Infinity) <= version && Math.max(...c.from) >= version ? 0 : 1;
  out.sort((x, y) => x.d - y.d);
  if (!generic || (version > 13006 && STABLE_TYPES.has(type))) return out;
  // Layouts nobody has published or embedded: edits of the nearest known layout by one inserted block of 1-4 fields, or,
  // for a node type that is not known at all, a run of anonymous scalar fields with an optional trailing variable array.
  const codes = generic.text ? ['d', 'c', 'l'] : ['p', 'd', 'u', 'f'], sigs = new Set(out.map((c) => c.fields.map((f) => f.name + f.code + f.n).join()));
  const push = (name, fields, d) => { const k = fields.map((f) => f.name + f.code + f.n).join(); if (!sigs.has(k)) { sigs.add(k); out.push({ name, fields, d, from: [], guessed: true, rank: 2 }); } };
  if (out.length) {
    let d = 2e6;
    for (const tpl of out.slice()) {
      const n = tpl.fields.length, at = tpl.fields.findIndex((f) => f.n === -1), last = at >= 0 ? at : n;
      for (let k = 1; k <= 4; k++) for (const code of codes) for (let i = last; i >= 0; i--) { const extra = []; for (let j = 0; j < k; j++) extra.push({ name: `unknown_${i}_${j}`, code, n: 0 }); push(tpl.name, [...tpl.fields.slice(0, i), ...extra, ...tpl.fields.slice(i)], d++); }
      for (let i = n - 1; i >= 0; i--) if (tpl.fields[i].n !== -1) push(tpl.name, [...tpl.fields.slice(0, i), ...tpl.fields.slice(i + 1)], d++);
      if (version <= 13006 && !generic.text) for (let i = 0; i < n; i++) if (tpl.fields[i].code === 'd' || tpl.fields[i].code === 'n') push(tpl.name, tpl.fields.map((f, q) => (q === i ? { ...f, code: f.code === 'd' ? 'n' : 'd' } : f)), d++);   // integer width differences in old binary files
    }
  } else {
    const code = generic.text ? 'd' : 'p'; let d = 3e6;
    for (let k = 0; k <= 24; k++) { const f = []; for (let j = 0; j < k; j++) f.push({ name: `unknown_${j}`, code, n: 0 }); push(`UNKNOWN_${type}`, [...f, { name: 'unknown_list', code, n: -1 }], d++); push(`UNKNOWN_${type}`, f, d++); }
  }
  return out;
}
const CURVE_T = (t) => (t >= 30 && t <= 39) || (t >= 130 && t <= 139), SURF_T = (t) => (t >= 50 && t <= 69) || (t >= 120 && t <= 124);
const TOPO_PTR = { face: [14], loop: [15], fin: [17], shell: [13], region: [19], body: [12], vertex: [18], edge: [16], point: [29], forward: [17], backward: [17], other: [17], next_at_vx: [17], front_face: [14], next_on_surface: [14], previous_on_surface: [14], next_on_curve: [16], previous_on_curve: [16] };
/** Structural check of a parsed node graph: every topological pointer must land on a node of the permitted type. Returns the number of violations. */
function graphViolations(nodes) {
  let bad = 0;
  for (const nd of nodes.values()) {
    if (nd.type < 12 || nd.type > 19) continue;
    for (const [k, v] of Object.entries(nd.f)) {
      if (typeof v !== 'number' || v <= 0 || k === 'node_id') continue;
      let okT = null;
      if (k === 'next' || k === 'previous') okT = (t) => t === nd.type;
      else if (k === 'surface') okT = SURF_T; else if (k === 'curve') okT = CURVE_T;
      else if (TOPO_PTR[k]) okT = (t) => TOPO_PTR[k].includes(t);
      else continue;
      const tg = nodes.get(v); if (!tg && nd.type === 12 && (k === 'next' || k === 'previous')) continue;      // bodies are chained through their partition, which need not be transmitted
      if (!tg || !okT(tg.type)) bad++;
    }
  }
  return bad;
}
const DEBUG_XT = typeof process === 'object' && !!process.env?.XT_DEBUG;
/**
 * How doubtful a parsed graph is: [node types nobody documents, nodes that no pointer refers to]. In a correctly decoded
 * stream practically every node is referenced by another one, so a layout that turns pointers into plain numbers shows
 * up as orphaned nodes.
 */
function graphDoubt(nodes) {
  const ref = new Set(), kinds = new Set();
  for (const nd of nodes.values()) {
    if (/^UNKNOWN_/.test(nd.kind)) kinds.add(nd.kind);
    for (const fd of nd.layout?.fields || []) { if (fd.code !== 'p' && !fd.name.startsWith('unknown_')) continue; const v = nd.f[fd.name]; if (Array.isArray(v)) { for (const x of v) if (x > 0) ref.add(x); } else if (v > 0) ref.add(v); }
  }
  let orphans = 0; for (const i of nodes.keys()) if (!ref.has(i)) orphans++;
  return [kinds.size, orphans];
}
const emptyValue = (v) => v === null || v === undefined || v === 0 || v === false || v === '' || (Array.isArray(v) && v.every(emptyValue));
// housekeeping fields that the translator never reads; two layouts that differ only in these are the same for our purposes
const BOOKKEEPING = new Set(['attribute_chains', 'list', 'key', 'highest_node_id', 'index_map_offset', 'index_map', 'node_id_index_map', 'schema_embedding_map', 'child', 'lowest_node_id', 'previous']);
/** Two parses of the same stream are equivalent when every node agrees on the fields both layouts name (housekeeping fields aside) and the other fields are empty. */
function graphDifferences(A, B) {
  const diff = new Set(); if (A.size !== B.size) return diff.add('node count');
  for (const [i, a] of A) {
    const b = B.get(i); if (!b || b.type !== a.type) return diff.add('node types');
    if (/^UNKNOWN_/.test(a.kind)) continue;                // node types nobody documents are not interpreted, only stepped over
    for (const k of new Set([...Object.keys(a.f), ...Object.keys(b.f)])) { if (k.startsWith('unknown_') || BOOKKEEPING.has(k)) continue; const x = a.f[k], y = b.f[k]; if (x === undefined || y === undefined) { if (!emptyValue(x ?? y)) diff.add(`${a.kind}.${k}`); } else if (x !== y && (typeof x !== 'object' || JSON.stringify(x) !== JSON.stringify(y))) diff.add(`${a.kind}.${k}`); }
    if (diff.size > 8) break;
  }
  return diff;
}
// differences that do not change geometry: the fields are blanked and the fact is reported
const SOFT_DIFFERENCES = new Set(['ASSEMBLY.attributes_groups']);
/**
 * Layouts for a file whose schema is not embedded: a depth-first search over the candidate layouts of each node type in
 * order of first use, nearest observed version first. A combination is a solution when the whole stream parses to its
 * terminator and every topological pointer lands on a node of the right type. The first solution is accepted only if
 * replacing the layout of any single node type by another candidate yields no different valid graph.
 */
function inferLayouts(run, version, text, budget = 60000) {
  let generic = null, cands = new Map(); const cand = (t) => { let c = cands.get(t); if (!c) { c = candidateLayouts(t, version, generic); cands.set(t, c); } return c; };
  let attempts = 0, firstError = null, solution = null, choice = null; const t0 = Date.now();
  for (let phase = 0; phase < 2 && !solution; phase++) {
    if (phase === 1) { generic = { text }; cands = new Map(); }       // second pass: also try layouts nobody has documented
    // depth-first search; the stream is re-read only from the first node of the type whose layout was changed
    const state = {}, order = [], marks = new Map(); choice = new Map(); let resume = null;
    const pick = (t, p, n) => { const c = cand(t); if (t < 10) throw new Error(`invalid node type ${t}`); if (!c.length) throw new Error(`node type ${t} is not in the documented schema nor in any schema learned from files that embed theirs`); if (c[0].name.startsWith('UNKNOWN_') && !marks.has(t) && order.filter((u) => cand(u)[0].name.startsWith('UNKNOWN_')).length >= 2) throw new Error('too many undocumented node types (wrong field layout)'); if (!marks.has(t)) { marks.set(t, { p, n }); order.push(t); } return c[choice.get(t) || 0]; };
    for (; attempts < 200000 && Date.now() - t0 < (version <= 13006 ? Math.min(20000, budget) : phase ? budget : budget / 3); attempts++) {
      let err = null;
      try { const r = run(pick, state, resume), bad = graphViolations(r.nodes); if (bad) err = new Error(`${bad} topological pointer(s) land on nodes of the wrong type`); else { solution = r; break; } } catch (e) { err = e; }
      firstError ||= err;
      if (DEBUG_XT) { const nn = +(/\((\d+) nodes read/.exec(err.message)?.[1] ?? -1); if (nn > (inferLayouts.best ?? -1)) { inferLayouts.best = nn; console.error('best', nn, attempts, err.message.slice(0, 160), JSON.stringify([...choice])); } }
      let moved = false;
      while (order.length) { const t = order[order.length - 1], k = (choice.get(t) || 0) + 1; if (k < cand(t).length) { choice.set(t, k); resume = marks.get(t); moved = true; break; } choice.delete(t); marks.delete(t); order.pop(); }
      if (!moved) break;
    }
    if (solution) solution = { result: solution, order: order.slice(), marks, state };
  }
  if (!solution) throw Object.assign(new Error(`the file uses schema ${version} without embedding it, and no combination of the known node layouts (documented V13 base, and layouts learned from schema versions ${xtLearnedVersions().join(', ')}) or of single-block edits of them parses it consistently (${attempts + 1} attempt(s)); with the nearest layouts: ${firstError?.message}`), { meta: firstError?.meta });
  // detach the accepted graph from the search state, then compare it with every single-type substitution: a less doubtful
  // graph replaces it, an equally doubtful different one makes the file ambiguous
  solution.result = { ...solution.result, nodes: new Map(solution.result.nodes) };
  const worse = (a, b) => a[0] - b[0] || a[1] - b[1], soft = new Set();
  let doubt = graphDoubt(solution.result.nodes);
  for (let round = 0, changed = true; changed && round < 6; round++) {
    changed = false;
    for (const t of solution.order.slice().reverse()) {     // latest first, so that the shared parse state can be rewound to each type's first node
      const k0 = choice.get(t) || 0, c0 = cand(t)[k0];
      for (let k = 0; k < cand(t).length && !changed; k++) {
        if (k === k0) continue;
        const alt = new Map(choice), order = []; alt.set(t, k); let r = null;
        const c1 = cand(t)[k], sig = (c) => c.fields.map((f) => f.code + f.n).join();
        if (sig(c1) === sig(c0)) {                           // same byte layout, other field names: compare the nodes of this type only
          if ((c1.rank ?? 2) > (c0.rank ?? 2)) continue;
          const A = new Map(), B = new Map();
          for (const nd of solution.result.nodes.values()) if (nd.type === t) { const vals = c0.fields.map((f) => nd.f[f.name]), f2 = {}; c1.fields.forEach((f, q) => { f2[f.name] = vals[q]; }); A.set(nd.index, nd); B.set(nd.index, { ...nd, f: f2 }); }
          const df = graphDifferences(A, B); if (!df.size) continue;
          if ([...df].every((x) => SOFT_DIFFERENCES.has(x))) { for (const x of df) soft.add(x); continue; }
          throw Object.assign(new Error(`the file uses schema ${version} without embedding it, and its node layouts cannot be determined uniquely: two different layouts of ${c0.name} nodes both parse the whole file consistently`), { meta: solution.result.meta });
        }
        const fast = solution.marks && solution.marks.has(t); if (fast) order.push(...solution.order.slice(0, solution.order.indexOf(t)));
        try { r = run((u) => { const c = cand(u); if (!c.length) throw new Error('unknown node type'); if (!order.includes(u)) order.push(u); return c[alt.get(u) || 0]; }, fast ? solution.state : {}, fast ? solution.marks.get(t) : null); if (graphViolations(r.nodes)) r = null; } catch { r = null; }
        if (!r) continue;
        const d2 = graphDoubt(r.nodes), cmp = worse(d2, doubt);
        if (cmp > 0) continue;
        if (cmp < 0) { solution = { result: { ...r, nodes: new Map(r.nodes) }, order }; choice = alt; doubt = d2; changed = true; break; }
        if ((cand(t)[k].rank ?? 2) > (c0.rank ?? 2)) continue;                 // weaker evidence than the accepted layout
        const df = graphDifferences(solution.result.nodes, r.nodes); if (!df.size) continue;
        if ([...df].every((x) => SOFT_DIFFERENCES.has(x))) { for (const x of df) soft.add(x); continue; }
        throw Object.assign(new Error(`the file uses schema ${version} without embedding it, and its node layouts cannot be determined uniquely: two different layouts of ${c0.name} nodes both parse the whole file consistently`), { meta: solution.result.meta });
      }
      if (changed) break;
    }
  }
  const inferred = [], guessed = [];
  for (const t of solution.order) { const c0 = cand(t)[choice.get(t) || 0]; if (c0.guessed) guessed.push(`${c0.name} (${c0.fields.filter((f) => /^unknown_/.test(f.name)).length} field(s) of unknown meaning)`); else if (c0.from.length && c0.from[0] !== 13006) inferred.push(`${c0.name} as in schema ${c0.from.join('/')}`); }
  for (const nd of solution.result.nodes.values()) { delete nd.layout; for (const x of soft) { const [kind, field] = x.split('.'); if (nd.kind === kind) nd.f[field] = 0; } }
  solution.result.meta.inferredSchema = { version, attempts: attempts + 1, layouts: inferred, guessedLayouts: guessed, unreferencedNodes: doubt[1], undecidedFields: [...soft] };
  return solution.result;
}

// ---- stream readers ----
class TextIn {
  constructor(s) { this.s = s; this.p = 0; }
  eof() { return this.p >= this.s.length; }
  ch() { return this.s[this.p++]; }
  peek() { return this.s[this.p]; }
  num() {
    const s = this.s; if (s[this.p] === '?') { this.p++; return null; }
    let e = this.p; while (e < s.length && s.charCodeAt(e) !== 32) e++;
    if (e === this.p || e - this.p > 40) throw new Error(`malformed number at offset ${this.p} of the XT stream`);
    const v = +s.slice(this.p, e); if (v !== v) throw new Error(`malformed number "${s.slice(this.p, Math.min(e, this.p + 16))}" at offset ${this.p} of the XT stream`);
    this.p = e + 1; return v;
  }
  short() { return this.num(); } int() { return this.num(); } byte() { return this.num(); } dbl() { return this.num(); } ptr() { return this.num() || 0; } posInt() { return this.num(); }
  logical() { const c = this.ch(); if (c !== 'T' && c !== 'F') throw new Error(`expected a logical at offset ${this.p - 1} of the XT stream`); return c === 'T'; }
  chars(n) { const v = this.s.substr(this.p, n); this.p += n; return v; }
  sstr() { const n = this.num(); if (!(n >= 0 && n < 4096)) throw new Error('bad schema string length'); return this.chars(n); }
  vec() { if (this.s[this.p] === '?') { this.p++; return null; } return [this.num(), this.num(), this.num()]; }
  flagByte() { return this.num(); }
}
class BinIn {
  constructor(b, p, le) { this.b = b; this.dv = new DataView(b.buffer, b.byteOffset, b.byteLength); this.p = p; this.le = le; }
  need(n) { if (this.p + n > this.b.length) throw new Error('the XT binary stream ends inside a node (truncated file)'); }
  eof() { return this.p >= this.b.length; }
  byte() { this.need(1); return this.b[this.p++]; }
  ch() { return String.fromCharCode(this.byte()); }
  peek() { return String.fromCharCode(this.b[this.p]); }
  short() { this.need(2); const v = this.dv.getInt16(this.p, this.le); this.p += 2; return v; }
  int() { this.need(4); const v = this.dv.getInt32(this.p, this.le); this.p += 4; return v === -32764 ? null : v; }
  dbl() { this.need(8); const v = this.dv.getFloat64(this.p, this.le); this.p += 8; return v === -3.14158e13 ? null : v; }
  ptr() { let q = 0, r = this.short(); if (r < 0) { q = this.short(); r = -r; } return q * 32767 + r - 1; }
  posInt() { return this.ptr(); }
  logical() { return this.byte() !== 0; }
  chars(n) { this.need(n); const v = str(this.b, this.p, this.p + n); this.p += n; return v; }
  sstr() { return this.chars(this.byte()); }
  vec() { const v = [this.dbl(), this.dbl(), this.dbl()]; return v[0] === null || v[1] === null || v[2] === null ? null : v; }
  flagByte() { return this.byte(); }
}

/** Parse the textual file header and the node stream. Returns { header, meta, nodes: Map(index → node) }. */
export function parseXT(b) {
  const cache = {};
  try { return parseStream(b, null, cache); }
  catch (e) {
    // no embedded schema: newer than the documented base, or an old file that does not fit the documented layouts
    const m = e.meta, version = e.needLayouts || (m && !m.embeddedSchema && m.schemaVersion && m.nodeCount > 0 ? m.schemaVersion : 0);
    if (!version) throw e;
    // a stream without its terminator is truncated: no layout can make it parse, so do not search
    const tail = m.encoding === 'text' ? /1 [01]\s*$/.test(str(b, Math.max(0, b.length - 64), b.length).replace(/[\r\n]/g, '')) : b.length > 8 && ((b[b.length - 4] === 0 && b[b.length - 3] === 1) || (b[b.length - 4] === 1 && b[b.length - 3] === 0) || (b[b.length - 8] | b[b.length - 7]) <= 1);
    if (!tail) throw Object.assign(new Error(e.needLayouts ? `the file uses schema ${version} without embedding it and ends without the stream terminator (truncated), so its node layouts cannot be inferred` : e.message), { meta: m });
    const budget = Math.min(90000, 5000 + b.length / 50);     // search time grows with the file: 5 s for small files, 45 s at 2 MB
    try { return inferLayouts((pick, state, resume) => parseStream(b, pick, cache, state, resume), version, m.encoding === 'text', budget); }
    catch (e2) { if (e.needLayouts) throw e2; throw e; }
  }
}
function parseStream(b, pick, cache, state = null, resume = null) {
  const headLen = Math.min(b.length, 1 << 16), head = cache.head || str(b, 0, headLen), hEnd = head.indexOf('**END_OF_HEADER'), header = {};
  let start = 0; cache.head = head;
  if (hEnd >= 0) { for (const m of head.slice(0, hEnd).replace(/\r?\n/g, '').matchAll(/([A-Z_0-9]+)=([^;]*);/g)) header[m[1]] = m[2].replace(/\^_/g, ' ').replace(/\^;/g, ';'); const nl = head.indexOf('\n', hEnd); start = nl < 0 ? headLen : nl + 1; }
  // byte offset of `start` (the header is ASCII, so character and byte offsets agree)
  const c0 = b[start], meta = { header, encoding: null, modellerVersion: null, schema: null, schemaVersion: null, embeddedSchema: false, userFieldSize: 0, nodeCount: 0 };
  let io, text = false;
  if (c0 === 0x54 /* T */) {
    text = true; meta.encoding = 'text'; io = new TextIn((cache.text ||= str(b, start, b.length).replace(/[\r\n]/g, ''))); io.p = 1;
    const mv = io.chars(io.num()), sv = io.chars(io.num()); meta.modeller = mv.trim(); meta.schema = sv;
  } else if (c0 === 0x50 && b[start + 1] === 0x53 /* PS\0\0: neutral binary, big-endian */) {
    meta.encoding = 'neutral binary'; io = new BinIn(b, start + 4, false);
    const mv = io.chars(io.short()), sl = io.int(); if (!(sl > 0 && sl < 200)) throw new Error('XT binary prefix is not understood (schema name length)'); meta.modeller = mv.trim(); meta.schema = io.chars(sl);
  } else if (c0 === 0x42 /* B: bare binary, taken as little-endian IEEE */) {
    meta.encoding = 'bare binary (little-endian assumed)'; io = new BinIn(b, start + 1, true);
    const ml = io.int(); if (!(ml > 0 && ml < 400)) throw new Error('XT bare-binary prefix is not understood'); const mv = io.chars(ml), sl = io.int(); if (!(sl > 0 && sl < 200)) throw new Error('XT bare-binary prefix is not understood'); meta.modeller = mv.trim(); meta.schema = io.chars(sl);
  } else throw new Error('no XT node stream was found after the file header (expected T, PS or B)');
  meta.modellerVersion = /version\s+(\d+)/.exec(meta.modeller)?.[1] ?? null;
  const sm = /^SCH_(\d+)_(\d+)(?:_(\d+))?/.exec(meta.schema) || []; meta.schemaVersion = +sm[2] || null; meta.embeddedSchema = !!sm[3]; meta.baseSchema = sm[3] ? +sm[3] : null;
  let maxTypes = 0;
  if (meta.embeddedSchema) maxTypes = text ? io.num() : io.short();
  meta.userFieldSize = (text ? io.num() : io.int()) || 0;
  if (!(meta.userFieldSize >= 0 && meta.userFieldSize <= 16)) throw new Error('XT user-field size is not plausible (the prefix was not decoded correctly)');
  if (!meta.embeddedSchema && meta.schemaVersion > 13006 && !pick) throw Object.assign(new Error('layouts needed'), { needLayouts: meta.schemaVersion, meta });
  const S = pick ? null : meta.embeddedSchema ? baseSchema() : legacySchema(baseSchema(), meta.schemaVersion || 13006), resolved = new Map(); let nodes = new Map(), seq = null;
  if (state) { if (resume && state.nodes) { nodes = state.nodes; seq = state.seq; while (seq.length > resume.n) nodes.delete(seq.pop()); } else { state.nodes = nodes; seq = state.seq = []; } }
  const readField = () => { const name = io.sstr(), cls = io.short(), n = io.posInt(); let code = 'p'; if (!cls) code = io.sstr(); let xmt = true; if (n === 1) xmt = text ? (io.peek() === 'T' || io.peek() === 'F' ? io.logical() : io.num() !== 0) : io.logical(); return { name, code, n: n === 1 ? -1 : n, cls, xmt }; };
  const schemaFor = (type, p0) => {
    if (pick) return pick(type, p0, nodes.size);
    if (!meta.embeddedSchema) { const d = S.get(type); if (!d) throw new Error(`node type ${type} is not in the documented schema`); return d; }
    let d = resolved.get(type); if (d) return d;
    const flag = io.flagByte(), base = S.get(type);
    if (flag === 255) { if (!base) throw new Error(`node type ${type} is declared equal to the base schema but is not in it`); d = base; }
    else {
      const next = io.peek(), isEdit = 'CDIAZ'.includes(next) && !!base;
      if (!isEdit) { const name = io.sstr(); io.sstr(); const fields = []; for (let i = 0; i < flag; i++) fields.push(readField()); d = { name, fields }; }
      else {
        const fields = []; let bi = 0;
        for (let guardN = 0; guardN < 2000; guardN++) { const op = io.ch(); if (op === 'Z') break; if (op === 'C') { if (bi >= base.fields.length) throw new Error(`schema edit for node type ${type} copies past the base definition`); fields.push(base.fields[bi++]); } else if (op === 'D') bi++; else if (op === 'I' || op === 'A') fields.push(readField()); else throw new Error(`unknown schema edit "${op}" for node type ${type}`); }
        d = { name: base.name, fields };
      }
    }
    d = { name: d.name, fields: d.fields.filter((f) => f.xmt !== false || f.n === -1) }; resolved.set(type, d); return d;
  };
  const value = (code) => {
    switch (code) {
      case 'd': case 'n': case 'w': case 'u': { const v = code === 'd' ? io.int() : code === 'u' ? io.byte() : io.short(); if (strict && v !== null && (!Number.isInteger(v) || (code === 'u' && (v < 0 || v > 255)))) throw new Error('non-integer value in an integer field (wrong field layout)'); return v; }
      case 'f': { const v = io.dbl(); if (strict && v !== null && !(Math.abs(v) < 1e30 && (v === 0 || Math.abs(v) > 1e-60 || cur > 19))) throw new Error('implausible real value (wrong field layout)'); return v; }
      case 'p': { const v = io.ptr(); if (strict && !(v >= 0 && v < 1e8 && Number.isInteger(v))) throw new Error('invalid pointer (wrong field layout)'); return v; }
      case 'c': { const v = io.ch(); if (strict && !(v >= ' ' && v <= '~')) throw new Error('non-printable character (wrong field layout)'); return v; }
      case 'l': { if (strict && !text) { const v = io.byte(); if (v > 1) throw new Error('invalid logical (wrong field layout)'); return v !== 0; } return io.logical(); } case 'v': case 'h': return io.vec(); case 'i': return [io.dbl(), io.dbl()]; case 'b': return [io.dbl(), io.dbl(), io.dbl(), io.dbl(), io.dbl(), io.dbl()];
      case 't': return io.int();
      default: throw new Error(`unknown field type "${code}" in the XT schema`);
    }
  };
  const strict = !!pick; let last = 'the stream prefix', cur = 0;
  if (resume) io.p = resume.p;
  try {
  for (let count = 0; count < 5e7; count++) {
    if (io.eof()) throw new Error('the XT node stream ends without its terminator (truncated file)');
    const p0 = io.p, type = io.short();
    if (type === 1) { const z = io.ptr(); if (z !== 0 && z !== 1) throw new Error('the XT stream lost synchronisation (a node was decoded with the wrong field layout)'); if (strict && (text ? io.s.slice(io.p).trim().length > 0 : io.b.length - io.p > 16)) throw new Error('a terminator was decoded before the end of the stream (wrong field layout)'); break; }   // terminator: 1 followed by index 0
    if (!(type >= 2 && type < 1000)) throw new Error(`invalid node type ${type} in the XT stream (offset ${io.p})`);
    cur = type; const d = schemaFor(type, p0), variable = d.fields.some((f) => f.n === -1), vlen = variable ? io.int() : 0;
    if (variable && !(vlen >= 0 && vlen <= 5e7)) throw new Error(`implausible variable length ${vlen} for a ${d.name} node`);
    const index = io.ptr(), f = {}; if (strict && (!(index > 0) || nodes.has(index))) throw new Error('invalid or repeated node index (wrong field layout)');
    for (const fd of d.fields) {
      if (fd.xmt === false) continue;
      if (fd.n === 0) f[fd.name] = value(fd.code);
      else { const k = fd.n === -1 ? vlen : fd.n; if (fd.code === 'c') { f[fd.name] = io.chars(k); if (strict && /[^\x20-\x7e]/.test(f[fd.name])) throw new Error('non-printable text (wrong field layout)'); } else { const a = new Array(k); for (let i = 0; i < k; i++) a[i] = value(fd.code); f[fd.name] = a; } }
    }
    for (let i = 0; i < meta.userFieldSize; i++) io.int();
    nodes.set(index, strict ? { type, kind: d.name, index, f, layout: d } : { type, kind: d.name, index, f }); if (seq) seq.push(index); last = `${d.name} node ${index}`;
  }
  } catch (e) { throw Object.assign(new Error(`${e.message} — in a node of type ${cur} after ${last} (${nodes.size} nodes read, schema ${meta.schema})`), { meta: { ...meta, nodeCount: nodes.size } }); }
  meta.nodeCount = nodes.size; meta.maxNodeTypes = maxTypes || null;
  return { header, meta, nodes, layouts: resolved };
}

// ---- geometry helpers on parsed nodes ----
const SURF_TYPES = new Set([50, 51, 52, 53, 54, 56, 59, 60, 67, 68, 120, 124]);
function nurbsOf(N, nurbs, isSurface) {
  const f = nurbs?.f; if (!f) return null;
  const verts = N.get(f.bspline_vertices)?.f.vertices, dim = f.vertex_dim; if (!verts || !(dim >= 2)) return null;
  const split = (count, at) => { const out = []; for (let i = 0; i < count; i++) out.push(verts.slice(at + i * dim, at + (i + 1) * dim)); return out; };
  if (!isSurface) {
    const knots = N.get(f.knots)?.f.knots, mults = N.get(f.knot_mult)?.f.mult; if (!knots || !mults) return null;
    return { degree: f.degree, poles: split(f.n_vertices, 0), knots, mults, U: expandKnots(knots, mults), rational: !!f.rational, dim, closed: !!f.closed, periodic: !!f.periodic };
  }
  const uk = N.get(f.u_knots)?.f.knots, um = N.get(f.u_knot_mult)?.f.mult, vk = N.get(f.v_knots)?.f.knots, vm = N.get(f.v_knot_mult)?.f.mult; if (!uk || !um || !vk || !vm) return null;
  const poles = []; for (let i = 0; i < f.n_u_vertices; i++) poles.push(split(f.n_v_vertices, i * f.n_v_vertices * dim));   // v varies fastest
  return { udeg: f.u_degree, vdeg: f.v_degree, poles, uknots: uk, umults: um, vknots: vk, vmults: vm, U: expandKnots(uk, um), Vk: expandKnots(vk, vm), rational: !!f.rational, dim, uclosed: !!f.u_closed, vclosed: !!f.v_closed };
}
/** Natural parametrisation of a curve node as { at(t), range: [t0, t1] } where it can be evaluated here. */
function curveEval(N, c) {
  const f = c.f;
  switch (c.kind) {
    case 'LINE': return { at: (t) => V.add(f.pvec, V.mul(f.direction, t)), range: null };
    case 'CIRCLE': case 'ELLIPSE': { const y = V.cross(f.normal, f.x_axis), a = f.radius ?? f.major_radius, bq = f.radius ?? f.minor_radius; return { at: (t) => V.add(f.centre, V.add(V.mul(f.x_axis, a * Math.cos(t)), V.mul(y, bq * Math.sin(t)))), range: [0, 2 * Math.PI], periodic: true }; }
    case 'B_CURVE': { const nb = nurbsOf(N, N.get(f.nurbs), false); if (!nb) return null; return { at: (t) => nurbsCurvePoint(nb, t), range: [nb.U[nb.degree], nb.U[nb.U.length - 1 - nb.degree]], nurbs: nb }; }
    default: return null;
  }
}
/** Natural parametrisation of a surface node as a function (u, v) → point, where it can be evaluated here. */
function surfaceEval(N, s, depth = 0) {
  const f = s.f; if (depth > 4) return null;
  const frame = (axis, x) => ({ a: axis, x, y: V.cross(axis, x) });
  switch (s.kind) {
    case 'PLANE': { const { x, y } = frame(f.normal, f.x_axis); return (u, v) => V.add(f.pvec, V.add(V.mul(x, u), V.mul(y, v))); }
    case 'CYLINDER': { const { a, x, y } = frame(f.axis, f.x_axis); return (u, v) => V.add(f.pvec, V.add(V.add(V.mul(x, f.radius * Math.cos(u)), V.mul(y, f.radius * Math.sin(u))), V.mul(a, v))); }
    case 'CONE': { const { a, x, y } = frame(f.axis, f.x_axis), tn = f.sin_half_angle / f.cos_half_angle; return (u, v) => { const r = f.radius + v * tn; return V.add(V.add(f.pvec, V.mul(a, v)), V.add(V.mul(x, r * Math.cos(u)), V.mul(y, r * Math.sin(u)))); }; }
    case 'SPHERE': { const { a, x, y } = frame(f.axis, f.x_axis); return (u, v) => V.add(f.centre, V.add(V.mul(V.add(V.mul(x, Math.cos(u)), V.mul(y, Math.sin(u))), f.radius * Math.cos(v)), V.mul(a, f.radius * Math.sin(v)))); }
    case 'TORUS': { const { a, x, y } = frame(f.axis, f.x_axis); return (u, v) => V.add(f.centre, V.add(V.mul(V.add(V.mul(x, Math.cos(u)), V.mul(y, Math.sin(u))), f.major_radius + f.minor_radius * Math.cos(v)), V.mul(a, f.minor_radius * Math.sin(v)))); }
    case 'B_SURFACE': { const nb = nurbsOf(N, N.get(f.nurbs), true); return nb ? (u, v) => nurbsSurfacePoint(nb, u, v) : null; }
    case 'SWEPT_SURF': { const c = N.get(f.section), ce = c && curveEval(N, c); return ce ? (u, v) => V.add(ce.at(u), V.mul(f.sweep, v)) : null; }
    case 'SPUN_SURF': { const c = N.get(f.profile), ce = c && curveEval(N, c); if (!ce) return null; return (u, v) => { const p = ce.at(u), z = V.add(f.base, V.mul(f.axis, V.dot(V.sub(p, f.base), f.axis))), d = V.sub(p, z); return V.add(z, V.add(V.mul(d, Math.cos(v)), V.mul(V.cross(f.axis, d), Math.sin(v)))); }; }
    case 'OFFSET_SURF': {
      const base = N.get(f.surface), be = base && surfaceEval(N, base, depth + 1); if (!be) return null; const sgn = base.f.sense === '-' ? -1 : 1;
      return (u, v) => { const h = 1e-6, p = be(u, v), n = V.unit(V.cross(V.sub(be(u + h, v), be(u - h, v)), V.sub(be(u, v + h), be(u, v - h)))); return V.add(p, V.mul(n, sgn * f.offset)); };
    }
    default: return null;
  }
}

const TYSA_COLOUR = 'SDL/TYSA_COLOUR', TYSA_NAME = 'SDL/TYSA_NAME';
/** Attribute values attached to a node: { definitionName: [value-node fields…] }. */
function attributesOf(N, node) {
  const out = {}; let a = N.get(node.f.attributes_groups);
  for (let g = 0; a && g < 1000; g++) { if (a.kind === 'ATTRIBUTE') { const def = N.get(a.f.definition), id = def && N.get(def.f.identifier)?.f.string; if (id && !(id in out)) out[id] = (a.f.fields || []).map((p) => N.get(p)?.f.values); } a = N.get(a.f.next); }
  return out;
}
const colourOf = (at) => { const v = at[TYSA_COLOUR]?.[0]; return Array.isArray(v) && v.length >= 3 && v.every((x) => x >= 0 && x <= 1) ? [v[0], v[1], v[2]] : null; };
const nameOf = (at) => { for (const k of [TYSA_NAME, 'SDL/TYSA_UNAME']) { const v = at[k]?.[0]; if (typeof v === 'string' && v.trim()) return v.trim(); if (Array.isArray(v) && v.length) { const s = String.fromCharCode(...v.filter((c) => c > 0)).trim(); if (s) return s; } } return null; };

/**
 * Translate a parsed XT node graph to STEP. Returns { step, stats } where stats counts what was translated and
 * what had to be skipped (by surface / curve type).
 */
/**
 * `ext` (optional) lets a container format drive the translation: { writer, stats, placements: [{ xf, name }] } emits the
 * file's root part once per placement into a shared writer and leaves finishing the STEP text to the caller.
 */
export function xtToStep(parsed, ext = null) {
  const N = parsed.nodes, W = ext?.writer || new StepWriter(), stats = ext?.stats || { bodies: 0, faces: 0, facesTranslated: 0, edges: 0, skippedFaces: {}, approximatedCurves: {}, instances: 0, wireOnlyBodies: 0, bodyNames: [] };
  const skip = (why) => { stats.skippedFaces[why] = (stats.skippedFaces[why] || 0) + 1; }, approx = (why) => { stats.approximatedCurves[why] = (stats.approximatedCurves[why] || 0) + 1; };
  const facesOfBody = new Map();
  for (const n of N.values()) if (n.kind === 'FACE') { const sh = N.get(n.f.shell), rg = sh && N.get(sh.f.region), body = rg ? rg.f.body : sh?.f.body; if (body) (facesOfBody.get(body) || facesOfBody.set(body, []).get(body)).push(n); }

  const emitBody = (body, xf, label) => {
    const faces = facesOfBody.get(body.index) || []; if (!faces.length) { stats.wireOnlyBodies++; return; }
    const P = xf ? (p) => V.mul(V.add([V.dot(xf.r[0], p), V.dot(xf.r[1], p), V.dot(xf.r[2], p)], xf.t), xf.s) : (p) => p;
    const D = xf ? (d) => [V.dot(xf.r[0], d), V.dot(xf.r[1], d), V.dot(xf.r[2], d)] : (d) => d, sc = xf ? xf.s : 1;
    const cache = new Map(), memo = (key, make) => { if (!cache.has(key)) cache.set(key, make()); return cache.get(key); };
    const vertexId = (vi) => memo('v' + vi, () => { const v = N.get(vi), pt = v && N.get(v.f.point); return pt && V.finite(pt.f.pvec) ? { id: W.vertex(P(pt.f.pvec)), p: pt.f.pvec } : null; });

    // sampled 3-D points of a curve that has no direct STEP equivalent (model coordinates), or null
    const sample = (c, t0 = null, t1 = null) => {
      const f = c.f;
      if (c.kind === 'INTERSECTION') { const h = N.get(f.chart)?.f.hvec; return h && h.length >= 2 && h.every(V.finite) ? h : null; }
      if (c.kind === 'SP_CURVE') {
        const s = N.get(f.surface), se = s && surfaceEval(N, s), bc = N.get(f.b_curve), nb = bc && nurbsOf(N, N.get(bc.f.nurbs), false); if (!se || !nb) return null;
        let a = nb.U[nb.degree], e = nb.U[nb.U.length - 1 - nb.degree]; if (t0 !== null && t1 !== null) { a = t0; e = t1; }
        const spans = Math.max(1, nb.poles.length - nb.degree), n = Math.min(400, Math.max(16, spans * 6)), out = [];
        for (let i = 0; i <= n; i++) { const uv = nurbsCurvePoint(nb, a + ((e - a) * i) / n), p = se(uv[0], uv[1]); if (!V.finite(p)) return null; out.push(p); }
        return out;
      }
      return null;
    };
    /** STEP curve for a curve node: { id, forward, ring } in the direction "start vertex → end vertex" of the edge. */
    const edgeCurve = (c, startP, endP, closed) => {
      const f = c.f;
      if (c.kind === 'TRIMMED_CURVE') {
        const basis = N.get(f.basis_curve); if (!basis) return null;
        if (basis.kind === 'SP_CURVE' || basis.kind === 'INTERSECTION') { const pts = sample(basis, basis.kind === 'SP_CURVE' ? f.parm_1 : null, basis.kind === 'SP_CURVE' ? f.parm_2 : null); return pts ? polyEdge(pts, startP, endP, closed, basis.kind, basis.kind === 'SP_CURVE') : null; }
        return edgeCurve(basis, startP, endP, closed);
      }
      const fwd = f.sense !== '-';
      switch (c.kind) {
        case 'LINE': return { id: memo('c' + c.index, () => W.line(P(f.pvec), D(f.direction))), forward: fwd };
        case 'CIRCLE': return { id: memo('c' + c.index, () => W.circle(P(f.centre), D(f.normal), D(f.x_axis), f.radius * sc)), forward: fwd, ring: V.add(f.centre, V.mul(f.x_axis, f.radius)) };
        case 'ELLIPSE': return { id: memo('c' + c.index, () => W.ellipse(P(f.centre), D(f.normal), D(f.x_axis), f.major_radius * sc, f.minor_radius * sc)), forward: fwd, ring: V.add(f.centre, V.mul(f.x_axis, f.major_radius)) };
        case 'B_CURVE': {
          const nb = nurbsOf(N, N.get(f.nurbs), false); if (!nb || nb.dim < 3) return null;
          const sum = nb.mults.reduce((s, m) => s + m, 0); if (sum !== nb.poles.length + nb.degree + 1) { const pts = []; const [a, e] = [nb.U[nb.degree], nb.U[nb.U.length - 1 - nb.degree]]; for (let i = 0; i <= 200; i++) pts.push(nurbsCurvePoint(nb, a + ((e - a) * i) / 200)); return polyEdge(pts, startP, endP, closed, 'periodic B_CURVE', true, fwd); }
          const id = memo('c' + c.index, () => W.bspline({ degree: nb.degree, poles: nb.poles.map((q) => P(nb.rational ? q.slice(0, 3).map((x) => x / q[nb.dim - 1]) : q.slice(0, 3))), knots: nb.knots, mults: nb.mults, weights: nb.rational ? nb.poles.map((q) => q[nb.dim - 1]) : null }));
          return { id, forward: fwd, ring: nurbsCurvePoint(nb, nb.U[nb.degree]).slice(0, 3) };
        }
        case 'INTERSECTION': case 'SP_CURVE': { const pts = sample(c); return pts ? polyEdge(pts, startP, endP, closed, c.kind, false, fwd) : null; }
        default: return null;
      }
    };
    /** Polyline stand-in for a curve, cut to the edge and directed start → end. */
    const polyEdge = (pts, startP, endP, closed, why, trimmed, fwd = true) => {
      approx(why);
      let q = pts;
      if (!closed && startP && endP) {
        const near = (p) => { let k = 0, best = Infinity; for (let i = 0; i < pts.length; i++) { const d = V.dist(pts[i], p); if (d < best) { best = d; k = i; } } return k; };
        let i0 = near(startP), i1 = near(endP);
        if (i0 === i1) q = [startP, endP];
        else { const mid = i0 < i1 ? pts.slice(i0 + 1, i1) : pts.slice(i1 + 1, i0).reverse(); q = [startP, ...mid.filter((p) => V.dist(p, startP) > 1e-12 && V.dist(p, endP) > 1e-12), endP]; }
      } else if (!fwd) q = pts.slice().reverse();
      return { id: W.polyline(q.map(P)), forward: true, ring: q[0], poly: true };
    };
    const edgeId = (ei) => memo('e' + ei, () => {
      const e = N.get(ei); if (!e) return null;
      const fins = []; for (let fi = e.f.fin, g = 0; fi && g < 64; g++) { const fn = N.get(fi); if (!fn || fins.includes(fn)) break; fins.push(fn); fi = fn.f.other; }
      // the forward vertex of a '+' fin is the edge end; of a '-' fin the edge start; the tail of a fin is its predecessor's forward vertex
      let sv = 0, ev = 0;
      for (const fn of fins) { const tail = N.get(fn.f.backward)?.f.vertex || 0; if (fn.f.sense === '+') { ev ||= fn.f.vertex; sv ||= fn.f.loop ? tail : 0; } else { sv ||= fn.f.vertex; ev ||= fn.f.loop ? tail : 0; } }
      const c = N.get(e.f.curve) || fins.map((fn) => N.get(fn.f.curve)).find(Boolean); if (!c) return null;
      const s = sv ? vertexId(sv) : null, t = ev ? vertexId(ev) : null, closed = !s || !t || sv === ev;
      const ec = edgeCurve(c, s?.p, t?.p, closed); if (!ec) return { fail: c.kind === 'TRIMMED_CURVE' ? `TRIMMED_CURVE of ${N.get(c.f.basis_curve)?.kind}` : c.kind };
      let v1 = s?.id, v2 = t?.id;
      if (!v1 || !v2) { if (!ec.ring) return { fail: `ring edge on ${c.kind}` }; v1 = v2 = W.vertex(P(ec.ring)); }
      stats.edges++;
      return { id: W.edge(v1, v2, ec.id, ec.forward) };
    });
    const surfaceId = (s) => memo('s' + s.index, () => {
      const f = s.f;
      switch (s.kind) {
        case 'PLANE': return { id: W.plane(P(f.pvec), D(f.normal), D(f.x_axis)) };
        case 'CYLINDER': return { id: W.cylinder(P(f.pvec), D(f.axis), D(f.x_axis), f.radius * sc) };
        // Real files show the cone widening along +axis from `radius` at pvec, with an outward natural normal (the sign
        // of v in the reference's formula does not match what Parasolid writes), so the STEP cone uses the same axis.
        case 'CONE': { const ang = Math.atan2(f.sin_half_angle, f.cos_half_angle); if (!(ang > 1e-9) || !(f.radius >= 0)) return { fail: 'CONE (degenerate angle)' }; return { id: W.cone(P(f.pvec), D(f.axis), D(f.x_axis), f.radius * sc, ang) }; }
        case 'SPHERE': return { id: W.sphere(P(f.centre), D(f.axis), D(f.x_axis), f.radius * sc), pole: V.add(f.centre, V.mul(f.axis, f.radius)) };
        case 'TORUS': if (!(f.major_radius > f.minor_radius && f.minor_radius > 0)) return { fail: 'TORUS (apple/lemon form)' }; return { id: W.torus(P(f.centre), D(f.axis), D(f.x_axis), f.major_radius * sc, f.minor_radius * sc), pole: V.add(f.centre, V.mul(f.x_axis, f.major_radius + f.minor_radius)) };
        case 'B_SURFACE': {
          const nb = nurbsOf(N, N.get(f.nurbs), true); if (!nb || nb.dim < 3) return { fail: 'B_SURFACE (no NURBS data)' };
          const su = nb.umults.reduce((a, m) => a + m, 0), sv = nb.vmults.reduce((a, m) => a + m, 0); if (su !== nb.poles.length + nb.udeg + 1 || sv !== nb.poles[0].length + nb.vdeg + 1) return { fail: 'B_SURFACE (periodic knot form)' };
          const w = nb.rational ? nb.poles.map((row) => row.map((q) => q[nb.dim - 1])) : null;
          return { id: W.bsurface({ udeg: nb.udeg, vdeg: nb.vdeg, poles: nb.poles.map((row) => row.map((q) => P(nb.rational ? q.slice(0, 3).map((x) => x / q[nb.dim - 1]) : q.slice(0, 3)))), uknots: nb.uknots, umults: nb.umults, vknots: nb.vknots, vmults: nb.vmults, weights: w }) };
        }
        case 'SWEPT_SURF': { const c = N.get(f.section), ec = c && ['LINE', 'CIRCLE', 'ELLIPSE', 'B_CURVE'].includes(c.kind) ? edgeCurve(c, null, null, true) : null; return ec && !ec.poly ? { id: W.extrusion(ec.id, D(f.sweep)), flip: c.f.sense === '-' } : { fail: `SWEPT_SURF of ${c?.kind}` }; }
        case 'SPUN_SURF': { const c = N.get(f.profile), ec = c && ['LINE', 'CIRCLE', 'ELLIPSE', 'B_CURVE'].includes(c.kind) ? edgeCurve(c, null, null, true) : null; return ec && !ec.poly ? { id: W.revolution(ec.id, P(f.base), D(f.axis)), flip: c.f.sense === '-' } : { fail: `SPUN_SURF of ${c?.kind}` }; }
        case 'OFFSET_SURF': { const base = N.get(f.surface), bs = base && SURF_TYPES.has(base.type) ? surfaceId(base) : null; if (!bs || bs.fail) return { fail: `OFFSET_SURF of ${base?.kind}` }; const sgn = (base.f.sense === '-' ? -1 : 1) * (bs.flip ? -1 : 1); return { id: W.offset(bs.id, sgn * f.offset * sc), flip: bs.flip }; }
        default: return { fail: s.kind };
      }
    });

    const shells = new Map(), faceColors = [];
    for (const face of faces) {
      stats.faces++;
      const s = N.get(face.f.surface); if (!s) { skip('face without a surface'); continue; }
      const sid = surfaceId(s); if (sid.fail) { skip(sid.fail); continue; }
      const bounds = []; let bad = null;
      for (let li = face.f.loop, g = 0; li && g < 1e5 && !bad; g++) {
        const lp = N.get(li); if (!lp) break;
        const ring = []; for (let fi = lp.f.fin, k = 0; fi && k < 1e5; k++) { const fn = N.get(fi); if (!fn) break; ring.push(fn); fi = fn.f.forward; if (fi === lp.f.fin) break; }
        const oriented = [];
        for (const fn of ring) { if (!fn.f.edge) continue; const e = edgeId(fn.f.edge); if (!e || e.fail) { bad = `edge on ${e?.fail || 'missing curve'}`; break; } oriented.push(W.oriented(e.id, fn.f.sense === '+')); }
        if (bad) break;
        if (oriented.length) bounds.push({ loop: W.loop(oriented) });
        else { const vtx = ring[0] && ring[0].f.vertex ? vertexId(ring[0].f.vertex) : null; if (vtx) bounds.push({ loop: W.vertexLoop(vtx.id) }); }
        li = lp.f.next;
      }
      if (bad) { skip(bad); continue; }
      if (!bounds.length) { if (!sid.pole) { skip(`unbounded ${s.kind} face`); continue; } bounds.push({ loop: W.vertexLoop(W.vertex(P(sid.pole))) }); }
      // the face normal is parallel to the natural surface normal when face and surface senses agree
      let same = (face.f.sense === '+') === (s.f.sense !== '-'); if (sid.flip) same = !same; if (xf && xf.mirror) same = !same;
      const at = attributesOf(N, face), fid = W.face(bounds, sid.id, same, nameOf(at) || ''), col = colourOf(at); if (col) faceColors.push([fid, col]);
      const key = face.f.shell || 0; (shells.get(key) || shells.set(key, { faces: [], closed: false }).get(key)).faces.push(fid);
      stats.facesTranslated++;
    }
    const at = attributesOf(N, body), name = label || nameOf(at) || `body ${stats.bodies + 1}`;
    W.body(name, [...shells.values()], { color: colourOf(at), faceColors }); stats.bodies++; stats.bodyNames.push(name);
  };

  // walk the root: a body, an assembly of instances, or a list of parts
  const seenAsm = new Set();
  // (xf = { r: 3×3 rows, t, s, mirror } acts as x' = (r·x + t)·s)
  const compose = (a, b) => { if (!a) return b; if (!b) return a; const r = [0, 1, 2].map((i) => [0, 1, 2].map((j) => a.r[i][0] * b.r[0][j] + a.r[i][1] * b.r[1][j] + a.r[i][2] * b.r[2][j])); const t = V.add(V.mul([V.dot(a.r[0], b.t), V.dot(a.r[1], b.t), V.dot(a.r[2], b.t)], b.s), a.t); return { r, t: V.mul(t, 1 / (b.s || 1)), s: a.s * b.s, mirror: a.mirror !== b.mirror }; };
  const xfOf = (ti) => { const t = N.get(ti); if (!t || t.kind !== 'TRANSFORM') return null; const m = t.f.rotation_matrix || [1, 0, 0, 0, 1, 0, 0, 0, 1], r = [[m[0], m[1], m[2]], [m[3], m[4], m[5]], [m[6], m[7], m[8]]], det = V.dot(r[0], V.cross(r[1], r[2])); return { r, t: t.f.translation_vector || [0, 0, 0], s: t.f.scale || 1, mirror: det < 0 }; };
  const walkPart = (pi, xf, depth, label) => {
    const p = N.get(pi); if (!p || depth > 32) return;
    if (p.kind === 'BODY') emitBody(p, xf, label);
    else if (p.kind === 'ASSEMBLY') {
      if (seenAsm.has(pi) && depth > 16) return; seenAsm.add(pi);
      for (let ii = p.f.sub_instance, g = 0; ii && g < 1e5; g++) { const inst = N.get(ii); if (!inst) break; stats.instances++; walkPart(inst.f.part, compose(xf, xfOf(inst.f.transform)), depth + 1, nameOf(attributesOf(N, inst))); ii = inst.f.next_in_part; }
    } else if (p.kind === 'POINTER_LIS_BLOCK') for (const e of p.f.entries || []) if (e) walkPart(e, xf, depth + 1, null);
  };
  const before = stats.bodies;
  for (const pl of ext?.placements || [{ xf: null, name: null }]) {
    const b0 = stats.bodies; walkPart(1, pl.xf, 0, pl.name);
    if (stats.bodies === b0) for (const n of N.values()) if (n.kind === 'BODY') emitBody(n, pl.xf, pl.name);      // root not understood: take every body as it stands
  }
  if (ext) return { step: null, stats, added: stats.bodies - before };
  return { step: stats.facesTranslated ? W.finish({ unit: '$', tol: 1e-6, source: 'Parasolid XT' }) : null, stats };
}
