// Shared back end for the kernel-format translators (Parasolid XT, ACIS SAT): a small vector toolkit, NURBS
// evaluation, and a writer that turns a boundary-representation graph into an in-memory STEP AP214 file.
// The translators parse the proprietary file into faces / loops / edges / curves / surfaces, describe them
// through StepWriter, and the vendored OpenCASCADE kernel then tessellates that STEP text like any other.

export const V = {
  add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]], sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], mul: (a, s) => [a[0] * s, a[1] * s, a[2] * s],
  dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2], cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
  len: (a) => Math.hypot(a[0], a[1], a[2]), dist: (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]),
  unit(a) { const l = Math.hypot(a[0], a[1], a[2]); return l > 0 ? [a[0] / l, a[1] / l, a[2] / l] : [0, 0, 1]; },
  /** Some unit vector perpendicular to n. */
  perp(n) { const a = Math.abs(n[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]; return V.unit(V.cross(V.cross(n, a), n)); },
  finite: (a) => Array.isArray(a) && a.length >= 3 && Number.isFinite(a[0]) && Number.isFinite(a[1]) && Number.isFinite(a[2]),
};

// ---------- NURBS evaluation ----------
/** Expand distinct knots + multiplicities into a flat knot vector. */
export function expandKnots(knots, mults) { const out = []; for (let i = 0; i < knots.length; i++) for (let k = 0; k < (mults[i] || 1); k++) out.push(knots[i]); return out; }
function findSpan(U, p, n, t) { if (t >= U[n]) return n - 1; if (t <= U[p]) return p; let lo = p, hi = n; while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (t < U[mid]) hi = mid; else lo = mid; } return lo; }
function basis(U, p, span, t) {
  const N = new Array(p + 1).fill(0), L = new Array(p + 1), R = new Array(p + 1); N[0] = 1;
  for (let j = 1; j <= p; j++) { L[j] = t - U[span + 1 - j]; R[j] = U[span + j] - t; let saved = 0; for (let r = 0; r < j; r++) { const d = R[r + 1] + L[j - r], tmp = d !== 0 ? N[r] / d : 0; N[r] = saved + R[r + 1] * tmp; saved = L[j - r] * tmp; } N[j] = saved; }
  return N;
}
/** Point on a B-spline curve. poles: array of dim-vectors (homogeneous, last = weight, when rational). */
export function nurbsCurvePoint({ degree: p, U, poles, rational }, t) {
  const n = poles.length, dim = poles[0].length, s = findSpan(U, p, n, t), N = basis(U, p, s, t), out = new Array(dim).fill(0);
  for (let j = 0; j <= p; j++) { const P = poles[Math.min(n - 1, Math.max(0, s - p + j))]; for (let d = 0; d < dim; d++) out[d] += N[j] * P[d]; }
  if (rational) { const w = out[dim - 1] || 1; return out.slice(0, dim - 1).map((v) => v / w); }
  return out;
}
/** Point on a B-spline surface. poles[i][j] with i along u. */
export function nurbsSurfacePoint({ udeg, vdeg, U, Vk, poles, rational }, u, v) {
  const nu = poles.length, nv = poles[0].length, dim = poles[0][0].length, su = findSpan(U, udeg, nu, u), sv = findSpan(Vk, vdeg, nv, v), Nu = basis(U, udeg, su, u), Nv = basis(Vk, vdeg, sv, v), out = new Array(dim).fill(0);
  for (let i = 0; i <= udeg; i++) { const row = poles[Math.min(nu - 1, Math.max(0, su - udeg + i))]; for (let j = 0; j <= vdeg; j++) { const P = row[Math.min(nv - 1, Math.max(0, sv - vdeg + j))], w = Nu[i] * Nv[j]; for (let d = 0; d < dim; d++) out[d] += w * P[d]; } }
  if (rational) { const w = out[dim - 1] || 1; return out.slice(0, dim - 1).map((x) => x / w); }
  return out;
}

// ---------- STEP AP214 writer ----------
const num = (v) => { if (!Number.isFinite(v)) v = 0; if (Number.isInteger(v) && Math.abs(v) < 1e15) return v + '.'; const s = String(v); return s.includes('e') ? s.replace(/^(-?\d+)e/, '$1.e').toUpperCase() : s; };
const quote = (s) => "'" + String(s ?? '').replace(/[^\x20-\x7e]/g, '_').replace(/'/g, "''").replace(/\\/g, '/') + "'";
const bool = (b) => (b ? '.T.' : '.F.');

/** Incremental writer for a STEP file made of shell-based surface models (one product per body). */
export class StepWriter {
  constructor() { this.L = []; this.id = 10; this.bodies = []; this.cache = new Map(); }
  add(s) { this.L.push(`#${++this.id}=${s};`); return `#${this.id}`; }
  point(p) { return this.add(`CARTESIAN_POINT('',(${num(p[0])},${num(p[1])},${num(p[2])}))`); }
  dir(d) { const u = V.unit(d); return this.add(`DIRECTION('',(${num(u[0])},${num(u[1])},${num(u[2])}))`); }
  /** Right-handed placement: origin, z axis and a reference direction (made perpendicular to z). */
  axis(o, z, x) { const zz = V.unit(z); let xx = x && V.finite(x) ? V.sub(x, V.mul(zz, V.dot(x, zz))) : null; if (!xx || V.len(xx) < 1e-9) xx = V.perp(zz); return this.add(`AXIS2_PLACEMENT_3D('',${this.point(o)},${this.dir(zz)},${this.dir(xx)})`); }
  vertex(p) { return this.add(`VERTEX_POINT('',${this.point(p)})`); }
  // curves
  line(p, d) { return this.add(`LINE('',${this.point(p)},${this.add(`VECTOR('',${this.dir(d)},1.)`)})`); }
  circle(c, n, x, r) { return this.add(`CIRCLE('',${this.axis(c, n, x)},${num(r)})`); }
  ellipse(c, n, x, a, b) { return this.add(`ELLIPSE('',${this.axis(c, n, x)},${num(a)},${num(b)})`); }
  /** B-spline curve from distinct knots + multiplicities; weights optional. */
  bspline({ degree, poles, knots, mults, weights = null, closed = false }) {
    const pts = `(${poles.map((p) => this.point(p)).join(',')})`, kn = `(${mults.join(',')}),(${knots.map(num).join(',')}),.UNSPECIFIED.`;
    if (!weights) return this.add(`B_SPLINE_CURVE_WITH_KNOTS('',${degree},${pts},.UNSPECIFIED.,${bool(closed)},.F.,${kn})`);
    return this.add(`( BOUNDED_CURVE() B_SPLINE_CURVE(${degree},${pts},.UNSPECIFIED.,${bool(closed)},.F.) B_SPLINE_CURVE_WITH_KNOTS(${kn}) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((${weights.map(num).join(',')})) REPRESENTATION_ITEM('') )`);
  }
  /** Degree-1 B-spline through the given points (chord-length knots). */
  polyline(points) {
    const pts = points.filter((p, i) => i === 0 || V.dist(p, points[i - 1]) > 0); if (pts.length < 2) pts.push(V.add(pts[0], [1e-9, 0, 0]));
    const knots = [0]; for (let i = 1; i < pts.length; i++) knots.push(knots[i - 1] + V.dist(pts[i], pts[i - 1]));
    return this.bspline({ degree: 1, poles: pts, knots, mults: knots.map((_, i) => (i === 0 || i === knots.length - 1 ? 2 : 1)) });
  }
  // surfaces
  plane(o, n, x) { return this.add(`PLANE('',${this.axis(o, n, x)})`); }
  cylinder(o, a, x, r) { return this.add(`CYLINDRICAL_SURFACE('',${this.axis(o, a, x)},${num(r)})`); }
  /** Cone whose radius grows along +a from r at o, with half angle `ang` (rad). */
  cone(o, a, x, r, ang) { return this.add(`CONICAL_SURFACE('',${this.axis(o, a, x)},${num(r)},${num(ang)})`); }
  sphere(c, a, x, r) { return this.add(`SPHERICAL_SURFACE('',${this.axis(c, a, x)},${num(r)})`); }
  torus(c, a, x, R, r) { return this.add(`TOROIDAL_SURFACE('',${this.axis(c, a, x)},${num(R)},${num(r)})`); }
  /** B-spline surface; poles[i][j] with i along u. */
  bsurface({ udeg, vdeg, poles, uknots, umults, vknots, vmults, weights = null, uclosed = false, vclosed = false }) {
    const pts = `(${poles.map((row) => `(${row.map((p) => this.point(p)).join(',')})`).join(',')})`, kn = `(${umults.join(',')}),(${vmults.join(',')}),(${uknots.map(num).join(',')}),(${vknots.map(num).join(',')}),.UNSPECIFIED.`;
    if (!weights) return this.add(`B_SPLINE_SURFACE_WITH_KNOTS('',${udeg},${vdeg},${pts},.UNSPECIFIED.,${bool(uclosed)},${bool(vclosed)},.F.,${kn})`);
    return this.add(`( BOUNDED_SURFACE() B_SPLINE_SURFACE(${udeg},${vdeg},${pts},.UNSPECIFIED.,${bool(uclosed)},${bool(vclosed)},.F.) B_SPLINE_SURFACE_WITH_KNOTS(${kn}) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE((${weights.map((row) => `(${row.map(num).join(',')})`).join(',')})) REPRESENTATION_ITEM('') SURFACE() )`);
  }
  extrusion(curve, d) { return this.add(`SURFACE_OF_LINEAR_EXTRUSION('',${curve},${this.add(`VECTOR('',${this.dir(d)},1.)`)})`); }
  revolution(curve, p, a) { return this.add(`SURFACE_OF_REVOLUTION('',${curve},${this.add(`AXIS1_PLACEMENT('',${this.point(p)},${this.dir(a)})`)})`); }
  offset(surface, d) { return this.add(`OFFSET_SURFACE('',${surface},${num(d)},.F.)`); }
  // topology
  edge(v1, v2, curve, same) { return this.add(`EDGE_CURVE('',${v1},${v2},${curve},${bool(same)})`); }
  oriented(edge, forward) { return this.add(`ORIENTED_EDGE('',*,*,${edge},${bool(forward)})`); }
  loop(orientedEdges) { return this.add(`EDGE_LOOP('',(${orientedEdges.join(',')}))`); }
  vertexLoop(v) { return this.add(`VERTEX_LOOP('',${v})`); }
  /** bounds: [{ loop, outer }]; `same` = face normal agrees with the surface normal. */
  face(bounds, surface, same, name = '') { return this.add(`ADVANCED_FACE(${quote(name)},(${bounds.map((b) => this.add(`${b.outer ? 'FACE_OUTER_BOUND' : 'FACE_BOUND'}('',${b.loop},.T.)`)).join(',')}),${surface},${bool(same)})`); }
  /** Register one body: name, its shells (arrays of face ids, closed flag), optional colour and per-face colours. */
  body(name, shells, { color = null, faceColors = [] } = {}) { this.bodies.push({ name, shells: shells.filter((s) => s.faces.length), color, faceColors }); }
  style(item, rgb) { const c = this.add(`COLOUR_RGB('',${num(rgb[0])},${num(rgb[1])},${num(rgb[2])})`); this.add(`STYLED_ITEM('',(${this.add(`PRESENTATION_STYLE_ASSIGNMENT((${this.add(`SURFACE_STYLE_USAGE(.BOTH.,${this.add(`SURFACE_SIDE_STYLE('',(${this.add(`SURFACE_STYLE_FILL_AREA(${this.add(`FILL_AREA_STYLE('',(${this.add(`FILL_AREA_STYLE_COLOUR('',${c})`)}))`)})`)}))`)})`)}))`)}),${item})`); }
  /** Finish the file. unit: SI prefix enum for the length unit ('$' = metre, '.MILLI.' …); tol: geometric uncertainty in that unit. */
  finish({ unit = '$', tol = 1e-6, source = 'translated' } = {}) {
    const lu = this.add(`( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(${unit},.METRE.) )`), au = this.add('( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )'), su = this.add('( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() )');
    const unc = this.add(`UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(${num(tol)}),${lu},'distance_accuracy_value','')`);
    const ctx = this.add(`( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((${unc})) GLOBAL_UNIT_ASSIGNED_CONTEXT((${lu},${au},${su})) REPRESENTATION_CONTEXT('','') )`);
    const app = this.add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')"); this.add(`APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2010,${app})`);
    let n = 0;
    for (const b of this.bodies) {
      if (!b.shells.length) continue;
      const shells = b.shells.map((s) => this.add(`${s.closed ? 'CLOSED_SHELL' : 'OPEN_SHELL'}('',(${s.faces.join(',')}))`)), model = this.add(`SHELL_BASED_SURFACE_MODEL(${quote(b.name)},(${shells.join(',')}))`);
      const rep = this.add(`MANIFOLD_SURFACE_SHAPE_REPRESENTATION(${quote(b.name)},(${model}),${ctx})`);
      const prod = this.add(`PRODUCT(${quote(b.name || `body ${n + 1}`)},${quote(b.name || `body ${n + 1}`)},'',(${this.add(`PRODUCT_CONTEXT('',${app},'mechanical')`)}))`); this.add(`PRODUCT_RELATED_PRODUCT_CATEGORY('part','',(${prod}))`);
      const pd = this.add(`PRODUCT_DEFINITION('design','',${this.add(`PRODUCT_DEFINITION_FORMATION('','',${prod})`)},${this.add(`PRODUCT_DEFINITION_CONTEXT('part definition',${app},'design')`)})`);
      this.add(`SHAPE_DEFINITION_REPRESENTATION(${this.add(`PRODUCT_DEFINITION_SHAPE('','',${pd})`)},${rep})`);
      if (b.color) this.style(model, b.color);
      for (const [face, rgb] of b.faceColors) this.style(face, rgb);
      n++;
    }
    return `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((${quote(source)}),'2;1');\nFILE_NAME('translated.stp','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n${this.L.join('\n')}\nENDSEC;\nEND-ISO-10303-21;\n`;
  }
}
