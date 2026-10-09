// Format registry and content-based detection for the geometry/mesh import layer.
// Every format the platform names has a capability profile here: what is read, what metadata survives,
// what is NOT read and - for formats that are not read natively - the validated conversion pathway.
// A filename extension alone proves nothing: detectFormat() sniffs the content and only falls back to the
// extension (with a low confidence) when the bytes carry no recognisable signature.

import { str, looksText, view } from './parsers-util.js';

const VIA_TESS = 'Export a tessellation (STL, OBJ, glTF or 3MF) from the originating CAD system, or mesh the geometry with Gmsh/Salome and import the .msh / .unv / .vtk result.';
const VIA_STEP = 'Export STEP AP242 (with tessellation if the CAD system offers it) or a tessellated STL/OBJ/glTF from the originating CAD system; for analysis meshes, mesh in Gmsh/Salome and import .msh.';
const VIA_MESH = 'Convert with the originating tool, meshio or Gmsh to Gmsh .msh (2.2/4.1 ASCII), VTK legacy ASCII, SU2, Nastran or UNV and import that file.';
const VIA_CLOUD = 'Convert with CloudCompare or PDAL to uncompressed LAS, PLY, PCD (ASCII/binary) or XYZ.';
const F = (id, name, ext, category, support, reads, preserves, limits, pathway = '') => ({ id, name, ext, category, support, reads, preserves, limits, pathway });

export const FORMATS = [
  // ---- tessellated surfaces (native)
  F('stl', 'STL (stereolithography)', ['stl'], 'tessellated', 'native', 'Triangles from ASCII and binary files (binary files whose header starts with "solid" are recognised by their size).', 'Solid names (ASCII) as groups; binary header text.', 'No units, no colour/attribute bytes, no shared vertices in the file (coincident vertices are identified at analysis time, welded only by an explicit heal).'),
  F('obj', 'Wavefront OBJ', ['obj'], 'tessellated', 'native', 'Vertices, faces (polygons fan-triangulated, negative indices), polylines.', 'Object/group names as groups; material and mtllib names in the metadata.', 'Normals, texture coordinates, free-form curves/surfaces and material files are not read. No units.'),
  F('ply', 'PLY (Stanford polygon)', ['ply'], 'tessellated', 'native', 'ASCII, binary little-endian and binary big-endian; vertex x/y/z, polygon faces (fan-triangulated) or bare point clouds.', 'Header comments, element/property census.', 'Colours, normals and custom properties are skipped. No units.'),
  F('off', 'OFF (object file format)', ['off'], 'tessellated', 'native', 'Vertices and polygon faces (fan-triangulated); COFF/NOFF variants read for geometry.', 'Counts only.', 'Colours and normals skipped; 4OFF/nOFF higher-dimensional variants not read. No units.'),
  F('gltf', 'glTF 2.0 / GLB', ['gltf', 'glb'], 'tessellated', 'native', 'Triangle primitives (lists, strips, fans) with POSITION and indices, scene-node transforms applied; .gltf with embedded base64 buffers (or companion .bin files) and binary .glb.', 'Node/mesh names as groups; generator and version; units are metres by specification; Y-up axis noted.', 'Draco/meshopt-compressed primitives, sparse accessors, skins, morph targets, materials and animations are not read.'),
  F('3mf', '3MF (3D manufacturing format)', ['3mf'], 'tessellated', 'native', 'ZIP container; mesh objects (vertices/triangles), components and build-item transforms.', 'Object names as groups; the model unit attribute becomes the model units.', 'Materials, colours, slices, beam lattices and production-extension multi-file models are not read.'),
  F('amf', 'AMF (additive manufacturing format)', ['amf'], 'tessellated', 'native', 'XML, plain or zipped: object meshes, volumes as groups.', 'Volume/object names; the unit attribute becomes the model units.', 'Curved triangles, constellations (instances), materials and textures are not read.'),
  F('json', 'Platform geometry JSON', ['json'], 'tessellated', 'native', 'The JSON written by exportModel(model, "json"): positions, triangles, lines, element blocks, groups.', 'Units, metadata, provenance log and warnings.', 'Only this platform\'s own schema.'),
  // ---- tessellated (partial)
  F('dae', 'COLLADA', ['dae'], 'tessellated', 'partial', 'Geometry meshes with <triangles>, <polylist> or <polygons> primitives (POSITION source).', 'Geometry names as groups; unit and up-axis from <asset>.', 'Scene-graph transforms and instancing are NOT applied (geometry is read in its local frame); no materials, skins or splines.', 'If the scene uses node transforms, bake them on export or export glTF/OBJ instead.'),
  F('vrml', 'VRML 2.0 (.wrl)', ['wrl', 'vrml'], 'tessellated', 'partial', 'IndexedFaceSet nodes with inline Coordinate points.', 'None.', 'Transform nodes, DEF/USE reuse, PROTOs, primitives (Box, Sphere…) and elevation grids are not evaluated.', 'Export OBJ, STL or glTF from the authoring tool when the scene uses transforms or prototypes.'),
  F('x3d', 'X3D (XML encoding)', ['x3d'], 'tessellated', 'partial', 'IndexedFaceSet and IndexedTriangleSet with a child Coordinate node.', 'DEF names as groups.', 'Transforms, USE references, primitives and the ClassicVRML/binary encodings are not evaluated.', 'Export OBJ, STL or glTF from the authoring tool when the scene uses transforms.'),
  // ---- point clouds
  F('xyz', 'XYZ / delimited point table', ['xyz', 'txt', 'csv', 'asc', 'pts', 'dat'], 'pointcloud', 'native', 'One point per line: x y z (space, comma, semicolon or tab separated); two-column files are read as planar coordinates (e.g. aerofoil sections) joined into a polyline.', 'Column count; skipped header lines.', 'Extra columns (intensity, colour, normals) are ignored. No units.'),
  F('pcd', 'PCD (Point Cloud Data)', ['pcd'], 'pointcloud', 'native', 'ASCII and uncompressed binary data; x/y/z fields.', 'Header fields, viewpoint, organised width/height.', 'binary_compressed (LZF) data is not decoded; non-finite points are dropped and counted.', VIA_CLOUD),
  F('las', 'LAS 1.x (uncompressed lidar)', ['las'], 'pointcloud', 'partial', 'Public header (version, scale, offset, bounding box, point count) and X/Y/Z of point records, subsampled by stride to opts.maxPoints (default 200 000).', 'Header metadata, system/software identifiers, VLR count.', 'Classification, intensity, colour, waveforms and the CRS records are not interpreted; units come from the CRS and are therefore left unset.', VIA_CLOUD),
  F('laz', 'LAZ (compressed LAS)', ['laz'], 'pointcloud', 'native', 'Public header and X/Y/Z of all point formats, decompressed by the vendored LASzip decoder (laz-perf) and subsampled by stride to opts.maxPoints (default 200 000).', 'Header metadata, system/software identifiers.', 'Classification, intensity, colour and the CRS records are not interpreted, so units are left unset. Every record must be decompressed in sequence, so very large files take time even when subsampled.', 'If the decoder cannot be loaded, decompress with LASzip/PDAL to .las, or convert to PLY/PCD/XYZ.'),
  F('e57', 'ASTM E57', ['e57'], 'pointcloud', 'metadata', 'File header and the XML section: scan names, point counts, field names and Cartesian bounds.', 'Scan census and bounds in the metadata.', 'Point records (bit-packed CompressedVector packets) are not decoded, so no points are shown.', VIA_CLOUD),
  // ---- meshes (native)
  F('vtk', 'VTK legacy (ASCII)', ['vtk'], 'mesh', 'native', 'POLYDATA (polygons, strips, lines, vertices), UNSTRUCTURED_GRID (linear and quadratic cell types) and STRUCTURED_GRID, including the 5.x OFFSETS/CONNECTIVITY layout.', 'Title, dataset type, names of point/cell data arrays.', 'BINARY legacy files, STRUCTURED_POINTS/RECTILINEAR_GRID and polyhedral cells are not read; data arrays are listed, not loaded.', VIA_MESH),
  F('gmsh', 'Gmsh MSH (2.2 and 4.1 ASCII)', ['msh', 'gmsh'], 'mesh', 'native', 'Nodes and elements (lines, triangles, quads, tets, hexes, prisms, pyramids and their quadratic forms).', 'Physical names and tags as groups (zones / boundaries), entity-to-physical mapping (4.1).', 'Binary MSH, MSH 1.0/4.0, periodic links, partitions, post-processing views and high-order (>2) elements are not read.', 'Re-save from Gmsh as "Version 2 ASCII" or "Version 4 ASCII" (-format msh22 / msh41).'),
  F('su2', 'SU2 native mesh', ['su2'], 'mesh', 'native', '2-D and 3-D element connectivity (VTK type codes), points and boundary markers.', 'MARKER_TAG names as boundary groups.', 'Multi-zone files (NZONE > 1) read the first zone only; periodic/FFD blocks ignored.'),
  F('nastran', 'Nastran bulk data', ['bdf', 'nas', 'dat', 'fem'], 'mesh', 'native', 'GRID, CTRIA3/6, CQUAD4/8, CTETRA, CHEXA, CPENTA, CPYRAM, CBAR/CBEAM/CROD/CONROD in free-field, 8-character and 16-character fixed-field form with continuations.', 'Property ids (PSHELL/PSOLID/PBAR…) and material ids as groups; shell thickness, MAT1 constants, card census and coordinate-system ids in the metadata.', 'Non-basic GRID coordinate systems (CP ≠ 0) are reported but NOT applied; INCLUDE files, loads, constraints, RBE/MPC and superelements are counted only. No units in the format.'),
  F('abaqus', 'Abaqus input deck', ['inp'], 'mesh', 'native', '*NODE, *ELEMENT (continuum, shell, membrane, beam/truss types mapped to linear/quadratic elements), parts with the first *INSTANCE translation/rotation applied.', '*ELSET names on elements as groups; *NSET/*ELSET, *MATERIAL and section definitions in the metadata.', 'Repeated instances of one part are not replicated; *INCLUDE, constraints, contact, loads and steps are not read. No units in the format.'),
  F('unv', 'I-deas universal file', ['unv'], 'mesh', 'native', 'Dataset 2411 (nodes), 2412 (elements: beams, shells, solids, linear and parabolic), 164 (units), 2467/2477 (groups).', 'Unit system; permanent group names; dataset census.', 'Results datasets, coordinate systems (2420) and physical property tables are listed only.'),
  // ---- meshes (partial)
  F('vtkxml', 'VTK XML (.vtu/.vtp…)', ['vtu', 'vtp', 'vts', 'vtr', 'vti', 'vtm', 'pvtu', 'pvtp'], 'mesh', 'partial', '.vtu UnstructuredGrid and .vtp PolyData pieces whose DataArrays are inline ASCII.', 'Piece counts and data-array names.', 'Binary (base64) and appended/compressed DataArrays, structured/rectilinear/image/multiblock files and parallel wrappers are recognised but not decoded.', 'Re-save from ParaView/VTK as ASCII (.vtu/.vtp) or as legacy VTK ASCII.'),
  F('cdb', 'ANSYS CDB archive', ['cdb'], 'mesh', 'partial', 'NBLOCK nodes and solid-format EBLOCK elements (SOLID185/186/187/285, SHELL181/281, BEAM188 and similar, with degenerate shapes).', 'ET element-type table, material ids as groups, /UNITS, component (CMBLOCK) census.', 'Non-solid EBLOCK layout, real constants/sections, loads and coordinate systems are not read.', VIA_MESH),
  F('fluent', 'Fluent / ANSYS mesh (.msh, ASCII)', ['msh', 'cas'], 'mesh', 'partial', 'Node sections and the zone census; boundary face zones (tri/quad/polygon faces) are read as a surface.', 'Zone ids, types and names as groups; node/face/cell totals.', 'Cell connectivity is NOT reconstructed (no volume elements); binary sections and compressed .msh.gz / .cas files are not read.', 'For volume cells, convert with meshio or the originating mesher to Gmsh .msh, SU2 or VTK ASCII and import that file.'),
  F('tecplot', 'Tecplot ASCII data', ['dat', 'tec', 'plt'], 'mesh', 'partial', 'FE zones (FETRIANGLE, FEQUADRILATERAL, FETETRAHEDRON, FEBRICK, FELINESEG) and ordered I/J/K zones, POINT packing (BLOCK packing when all variables are nodal).', 'Zone titles as groups; variable names; title.', 'Cell-centred variables, polyhedral/polygonal zones, shared connectivity and face-neighbour data are not read; field variables are listed, not loaded.', 'Write the zone with DATAPACKING=POINT, or convert with Tecplot/meshio to VTK ASCII.'),
  F('plot3d', 'Plot3D structured grid (ASCII)', ['p3d', 'xyz', 'x', 'g', 'grd', 'fmt'], 'mesh', 'partial', 'Formatted single- and multi-block grids (whole format, 3-D; 2-D when the header has two dimensions) → hexahedra / quads.', 'Block dimensions; blocks as zones.', 'Unformatted/binary files, IBLANK interpretation and solution (.q) files are not read.', 'Write the grid as formatted (ASCII) whole-format Plot3D, or convert to VTK ASCII / Gmsh .msh.'),
  F('plot3dq', 'Plot3D solution file', ['q', 'f', 'fun'], 'mesh', 'metadata', 'Recognised by extension only.', 'None.', 'Solution/function files contain no coordinates.', 'Import the matching grid file (.x/.xyz/.p3d, ASCII) for geometry.'),
  F('openfoam', 'OpenFOAM polyMesh', ['foam'], 'mesh', 'partial', 'ASCII `points` and `faces` (plus `boundary` and `owner` when supplied as companion files): boundary patches as surface groups.', 'Patch names and types; cell count from `owner`.', 'Polyhedral cells are NOT reconstructed as volume elements; binary and compressed (.gz) field files are not read.', 'Supply constant/polyMesh/points, faces and boundary together, or run foamToVTK / foamMeshToFluent and import the result.'),
  // ---- meshes (metadata only)
  F('cgns', 'CGNS (HDF5)', ['cgns'], 'mesh', 'partial', 'HDF5-flavour CGNS through the vendored h5wasm kernel: bases, zones, GridCoordinates, structured zones (→ hexahedra / quads), Elements_t sections of fixed type and MIXED, NGON_n/NFACE_n boundary faces.', 'Section names, ZoneBC names / BC types / families as groups (when given on faces or elements), DimensionalUnits, CGNS version, zone and family census.', 'Legacy ADF-container files are recognised but NOT decoded; polyhedral cells are shown by their boundary faces only; vertex-located and structured-zone BCs, zone connectivity, flow solutions and linked files are listed or ignored, not applied; cylindrical coordinates are not converted.', 'Convert ADF files with `cgnsconvert -h in.cgns out.cgns`; for polyhedral meshes export tetra/hexa/prism cells or convert with the originating mesher to Gmsh .msh / SU2.'),
  F('exodus', 'Exodus II (NetCDF classic or NetCDF-4)', ['e', 'exo', 'ex2', 'exii', 'g', 'gen', 'nc'], 'mesh', 'partial', 'Coordinates, element blocks (HEX, TET, WEDGE, PYRAMID, QUAD/SHELL, TRI, BAR/BEAM families, linear and quadratic), side sets as boundary faces, node-set census. Classic NetCDF (CDF-1/2/5) is parsed natively; NetCDF-4 files go through the vendored h5wasm kernel.', 'Block ids and names as zones, side-set and node-set names as boundary groups, title and versions.', 'Results variables, element attributes, maps, polyhedral (NSIDED/NFACED) and superelement blocks are not read; side sets on shells in 3-D meshes are not resolved.', VIA_MESH),
  F('med', 'MED (Salome)', ['med', 'rmed'], 'mesh', 'partial', 'Unstructured MED 3.x/4.x meshes through the vendored h5wasm kernel: node coordinates and linear/quadratic segments, triangles, quadrangles, tetrahedra, hexahedra, pentahedra and pyramids.', 'Families resolved to their group names as groups; mesh name, description, axis unit, MED version.', 'Structured MED grids, polygons/polyhedra, fields (results) and node families are not read. Volume cells are renumbered to the platform orientation convention (logged); only the first computation step of each mesh is read.', 'Export UNV or Gmsh .msh from Salome, or convert with meshio.'),
  F('fluent-h5', 'Fluent CFF mesh/case (.msh.h5 / .cas.h5)', ['h5'], 'mesh', 'partial', 'Node coordinates and the boundary faces of the face zones (through the vendored h5wasm kernel), read from /meshes/*/nodes and /meshes/*/faces.', 'Face-zone and cell-zone names and ids as groups; node/face/cell totals.', 'Cell connectivity is NOT reconstructed (no volume elements); settings and results are ignored. The layout follows the published CFF structure and has been verified against synthetic files only.', 'For volume cells write an ASCII .msh or convert with the originating mesher / meshio to Gmsh .msh, SU2 or VTK ASCII.'),
  F('hdf5', 'HDF5 container (generic)', ['h5', 'hdf5', 'hdf', 'cgns', 'med'], 'mesh', 'metadata', 'The schema is identified from the content (CGNS, MED, Exodus/NetCDF-4, Fluent CFF are then read by their own readers); for anything else the object tree is listed.', 'Top of the object tree in the metadata.', 'HDF5 files that follow none of the mesh schemas above cannot be interpreted as geometry.', VIA_MESH),
  F('tecplot-bin', 'Tecplot binary (.plt / .szplt)', ['plt', 'szplt'], 'mesh', 'metadata', 'Magic number and file version.', 'None.', 'Binary zones are not decoded.', 'Write the data as ASCII .dat from Tecplot (File ▸ Write Data File ▸ ASCII, point format), or convert with meshio to VTK ASCII.'),
  // ---- CAD (exact geometry through the OpenCASCADE kernel)
  F('step', 'STEP (ISO 10303-21: AP203 / AP214 / AP242)', ['step', 'stp', 'p21'], 'cad', 'native', 'Exact B-rep solids, shells and faces tessellated by the vendored OpenCASCADE kernel (deflection set by opts.linearDeflection / opts.angularDeflection), with assembly placements applied; plus, from the text parser, the header, schema/AP, length unit, product names, entity census, and AP242 tessellated faces.', 'One group per body with its product/instance name and colour; per-face triangle ranges and colours (meta.brepFaces, elements[0].face); schema, originating system, author, units.', 'PMI, annotations, validation properties and the assembly tree itself (bodies are flattened into world coordinates) are not kept; wireframe-only and compressed .stpZ files give no faces; tessellation is an approximation whose accuracy is the chosen deflection.', 'If the OpenCASCADE kernel cannot be loaded (offline copy without js/vendor/occt, or WebAssembly blocked by the content-security policy) only the text-level content is shown; then export a tessellation (STL/OBJ/glTF) or a Gmsh mesh from the CAD system.'),
  F('iges', 'IGES 5.x', ['iges', 'igs'], 'cad', 'native', 'Surfaces, trimmed/bounded faces and solids tessellated by the vendored OpenCASCADE kernel; plus, from the text parser, the global section (units, product, system) and the entity-type census.', 'One group per body/face set with its name and colour; units, originating system, author, census.', 'IGES surfaces are usually not sewn, so the result is often an open, unstitched set of faces (heal() welds coincident edges only); annotation and drafting entities are ignored; compressed/binary IGES is not read.', 'If the OpenCASCADE kernel cannot be loaded (offline copy without js/vendor/occt, or WebAssembly blocked by the content-security policy) only the text-level content is shown; then export a tessellation (STL/OBJ/glTF) or a Gmsh mesh from the CAD system.'),
  F('brep', 'OpenCASCADE BREP (.brep)', ['brep', 'brp'], 'cad', 'native', 'Exact shapes in the OpenCASCADE text format, tessellated by the vendored OpenCASCADE kernel.', 'One group per body; per-face triangle ranges.', 'The format carries no units, names or colours; binary BREP (.bbrep) is not read.', 'If the OpenCASCADE kernel cannot be loaded (offline copy without js/vendor/occt, or WebAssembly blocked by the content-security policy) only the text-level content is shown; then export a tessellation (STL/OBJ/glTF) or a Gmsh mesh from the CAD system.'),
  F('xt', 'Parasolid text (.x_t)', ['x_t', 'xmt_txt'], 'cad', 'metadata', 'Transmit-file header fields (application, date, key), modeller version and schema.', 'Header metadata.', 'No open-source reader exists for the Parasolid kernel format: topology and geometry are not decoded, so no entity census, extents or tessellation are produced.', VIA_STEP),
  // ---- CAD (metadata only)
  F('xb', 'Parasolid binary (.x_b)', ['x_b', 'xmt_bin'], 'cad', 'metadata', 'Signature only.', 'None.', 'No open-source reader exists for the Parasolid kernel format; binary transmit data is not decoded.', VIA_STEP),
  F('acis', 'ACIS SAT / SAB', ['sat', 'sab'], 'cad', 'metadata', 'SAT text header (version, product, date, units scale) and record-type census; SAB by signature.', 'Header metadata; units when the header scale is a standard value.', 'No open-source reader evaluates ACIS geometry: surfaces and solids are not tessellated.', VIA_STEP),
  F('jt', 'JT (Siemens)', ['jt'], 'cad', 'metadata', 'Version string from the file header.', 'None.', 'No open-source reader exists for JT compressed LOD / XT B-rep segments; nothing beyond the header is decoded.', VIA_STEP),
  F('nativecad', 'Native CAD (CATIA, SolidWorks, NX, Creo, Inventor, Fusion, Rhino, FreeCAD…)', ['catpart', 'catproduct', 'sldprt', 'sldasm', 'prt', 'asm', 'ipt', 'iam', 'par', 'psm', 'f3d', '3dm', 'fcstd', 'dwg', 'dxf', 'ifc', 'skp', 'blend', 'fbx', '3ds'], 'cad', 'metadata', 'Recognised by extension and container signature only.', 'None.', 'Proprietary native containers have no open reader and are not decoded.', VIA_TESS),
  // ---- geospatial
  F('geotiff', 'GeoTIFF terrain raster', ['tif', 'tiff', 'gtiff'], 'geospatial', 'partial', 'Single-band elevation rasters in classic TIFF or BigTIFF: strips or tiles, uncompressed / LZW / Deflate / PackBits, 8/16/32-bit integer and 32/64-bit float samples, horizontal and floating-point predictors → a height-field surface sampled to at most opts.maxGrid (default 300) points per side.', 'ModelPixelScale, ModelTiepoint, ModelTransformation, GeoKeys (model/raster type, CRS codes, units, citations), GDAL no-data value; metres/feet as units when a projected CRS states them.', 'JPEG, JPEG 2000, ZSTD, LERC, LZMA and WebP compression, multi-band imagery, palette images, overviews and masks are not read; coordinates are not reprojected (a geographic CRS leaves x/y in degrees).', 'Convert with GDAL: gdal_translate -co COMPRESS=DEFLATE (or LZW) in.tif out.tif, and gdalwarp to a projected CRS in metres.'),
];

const byId = Object.fromEntries(FORMATS.map((f) => [f.id, f]));
export const formatById = (id) => byId[id] ?? null;

const startsWith = (b, sig, off = 0) => { if (b.length < off + sig.length) return false; for (let i = 0; i < sig.length; i++) if (b[off + i] !== (typeof sig === 'string' ? sig.charCodeAt(i) : sig[i])) return false; return true; };
const HDF5 = [0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a];
const PLOT3D_EXT = ['p3d', 'xyz', 'x', 'g', 'grd', 'fmt'];
const intTok = (t) => /^[+-]?\d+$/.test(t);
const numTok = (t) => t !== '' && Number.isFinite(+t.replace(/[dD]/, 'e'));

/** Identify by content. Returns { id, c, note } or null. */
function sniff(b, ext, base) {
  const n = b.length;
  if (startsWith(b, HDF5) || startsWith(b, HDF5, 512)) {
    const id = ext === 'cgns' ? 'cgns' : ext === 'med' || ext === 'rmed' ? 'med' : /\.(msh|cas|dat)\.h5$/.test(base) ? 'fluent-h5' : ['e', 'exo', 'ex2', 'exii', 'g', 'gen'].includes(ext) ? 'exodus' : 'hdf5';
    return { id, c: 0.9, note: 'HDF5 signature (\\x89HDF)' + (id === 'hdf5' ? '; the kind of HDF5 content cannot be told from the signature' : `; taken as ${byId[id].name} from the extension`) };
  }
  if (startsWith(b, 'CDF') && n > 3 && [1, 2, 5].includes(b[3])) return { id: 'exodus', c: 0.85, note: `NetCDF classic signature (CDF version ${b[3]}); Exodus II files use this container` };
  if (startsWith(b, '@(#)ADF Database')) return { id: 'cgns', c: 0.9, note: 'ADF database signature (legacy CGNS container)' };
  if (startsWith(b, 'glTF')) return { id: 'gltf', c: 0.99, note: 'GLB magic number' };
  if (startsWith(b, [0x50, 0x4b, 0x03, 0x04])) {
    const head = str(b, 0, Math.min(n, 65536));
    if (ext === '3mf' || /3D\/[^\0]*\.model/i.test(head)) return { id: '3mf', c: 0.95, note: 'ZIP container with a 3MF model part' };
    if (ext === 'amf') return { id: 'amf', c: 0.85, note: 'ZIP container with .amf extension (zipped AMF)' };
    if (byId.nativecad.ext.includes(ext)) return { id: 'nativecad', c: 0.7, note: 'ZIP-based native CAD container' };
    return null;
  }
  if (startsWith(b, 'LASF')) {
    const fmt = n > 104 ? b[104] : 0, laz = ext === 'laz' || (fmt & 0xc0) !== 0;
    return { id: laz ? 'laz' : 'las', c: 0.97, note: laz ? 'LASF signature with compressed (LASzip) point records' : 'LASF signature' };
  }
  if (startsWith(b, 'ASTM-E57')) return { id: 'e57', c: 0.99, note: 'ASTM-E57 signature' };
  if (startsWith(b, [0x49, 0x49, 0x2a, 0x00]) || startsWith(b, [0x4d, 0x4d, 0x00, 0x2a]) || startsWith(b, [0x49, 0x49, 0x2b, 0x00]) || startsWith(b, [0x4d, 0x4d, 0x00, 0x2b])) return { id: 'geotiff', c: 0.9, note: 'TIFF signature' };
  if (startsWith(b, '#!TDV') || startsWith(b, '#!SZPLT')) return { id: 'tecplot-bin', c: 0.99, note: 'Tecplot binary magic number' };
  if (startsWith(b, 'ply\n') || startsWith(b, 'ply\r')) return { id: 'ply', c: 0.99, note: 'PLY magic number' };
  if (startsWith(b, 'ACIS BinaryFile')) return { id: 'acis', c: 0.99, note: 'ACIS binary (SAB) signature' };
  if (startsWith(b, 'Version ') && /^Version \d+\.\d+ JT/.test(str(b, 0, 40))) return { id: 'jt', c: 0.99, note: 'JT version header' };
  if (startsWith(b, 'PS\0\0')) return { id: 'xb', c: 0.9, note: 'Parasolid binary transmit signature' };
  if (startsWith(b, [0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1])) return { id: 'nativecad', c: 0.7, note: 'OLE compound document (used by several native CAD formats)' };
  if (startsWith(b, 'V5_CFV2')) return { id: 'nativecad', c: 0.9, note: 'CATIA V5 container signature' };

  const head = str(b, b[0] === 0xef ? 3 : 0, Math.min(n, 8192)), lead = head.trimStart(), low = lead.toLowerCase();
  // STL: a binary file is identified by its size (84 + 50·n) even when the 80-byte header starts with "solid".
  if (n >= 84) {
    const nt = view(b).getUint32(80, true), asciiLike = low.startsWith('solid') && /\bfacet\b|\bendsolid\b/.test(low);
    if (84 + 50 * nt === n && !(asciiLike && nt === 0)) return { id: 'stl', c: asciiLike ? 0.8 : 0.95, note: low.startsWith('solid') ? 'binary STL (size matches 84 + 50·n) whose header begins with "solid"' : 'binary STL (size matches 84 + 50·n)' };
  }
  if (!looksText(b)) {
    if (ext === 'stl' && n >= 84) return { id: 'stl', c: 0.5, note: 'binary data with .stl extension, but the size does not match the triangle count' };
    if (ext === 'x_b' || ext === 'xmt_bin') return { id: 'xb', c: 0.6, note: 'binary data with Parasolid binary extension' };
    if (ext === 'sab') return { id: 'acis', c: 0.6, note: 'binary data with .sab extension' };
    if (PLOT3D_EXT.includes(ext) || ext === 'q') return { id: ext === 'q' ? 'plot3dq' : 'plot3d', c: 0.4, note: 'binary data with a Plot3D extension (unformatted Plot3D is not read)' };
    return null;
  }
  if (low.startsWith('iso-10303-21')) return { id: 'step', c: 0.99, note: 'ISO-10303-21 header' };
  if (/^(DBRep_DrawableShape|CASCADE Topology V\d)/.test(lead)) return { id: 'brep', c: 0.98, note: 'OpenCASCADE BREP header' };
  if (lead.startsWith('**ABCDEFGHIJKLMNOPQRSTUVWXYZ') || /\*\*PARASOLID/.test(head)) return { id: /FORMAT=binary/i.test(head) ? 'xb' : 'xt', c: 0.95, note: 'Parasolid transmit-file header' };
  if (head.length >= 80 && 'SCB'.includes(head[72]) && /^.{72}[SCB]\s*\d+\s*(\r?\n|$)/.test(head)) return { id: 'iges', c: 0.95, note: 'IGES fixed 80-column records (section letter in column 73)' };
  if (lead.startsWith('# vtk DataFile')) return { id: 'vtk', c: 0.99, note: 'VTK legacy header' };
  if (/<VTKFile\b/.test(head)) return { id: 'vtkxml', c: 0.99, note: 'VTK XML root element' };
  if (lead.startsWith('$MeshFormat') || lead.startsWith('$NOD') || lead.startsWith('$PhysicalNames')) return { id: 'gmsh', c: 0.99, note: 'Gmsh section header' };
  if (/<COLLADA\b/.test(head)) return { id: 'dae', c: 0.99, note: 'COLLADA root element' };
  if (/<X3D\b/.test(head)) return { id: 'x3d', c: 0.99, note: 'X3D root element' };
  if (lead.startsWith('#VRML') || lead.startsWith('#X3D')) return { id: 'vrml', c: 0.95, note: 'VRML header' };
  if (/<amf\b/i.test(head)) return { id: 'amf', c: 0.99, note: 'AMF root element' };
  if (/<model\b[^>]*3dmanufacturing/i.test(head)) return { id: '3mf', c: 0.9, note: 'bare 3MF model XML (not zipped)' };
  if (/\bFoamFile\b/.test(head)) return { id: 'openfoam', c: 0.97, note: 'OpenFOAM FoamFile header' };
  if (/^# \.PCD|^VERSION\s+[.\d]+/m.test(head) && /^FIELDS\s/m.test(head)) return { id: 'pcd', c: 0.98, note: 'PCD header' };
  if (lead[0] === '{') {
    if (/"aerosuiteGeometry"/.test(head)) return { id: 'json', c: 0.99, note: 'platform geometry JSON' };
    if (/"asset"\s*:/.test(str(b, 0, Math.min(n, 1 << 20))) || ext === 'gltf') return { id: 'gltf', c: 0.9, note: 'glTF JSON (asset object)' };
    return null;
  }
  if (low.startsWith('solid')) {
    if (/\bfacet\b|\bendsolid\b/.test(low)) return { id: 'stl', c: 0.97, note: 'ASCII STL keywords' };
    if (ext === 'stl') return { id: 'stl', c: 0.6, note: 'starts with "solid" but no facets were seen in the first block' };
  }
  if (/^(ST)?C?N?4?n?OFF\b/.test(lead)) return { id: 'off', c: 0.97, note: 'OFF keyword' };
  if (/^-1\s*\r?\n\s*\d+\s*\r?\n/.test(lead)) return { id: 'unv', c: 0.95, note: 'universal-file dataset delimiter (-1) followed by a dataset number' };
  if (/^\s*NDIME\s*=/m.test(head)) return { id: 'su2', c: 0.98, note: 'SU2 NDIME keyword' };
  if (lead[0] === '(' && /\(\s*(2\s+[23]\s*\)|10\s*\(|12\s*\(|13\s*\(|0\s*"|1\s*")/.test(head)) return { id: 'fluent', c: 0.9, note: 'Fluent mesh section headers' };
  if (/^\s*(NBLOCK|EBLOCK|\/PREP7|\/COM,\s*ANSYS|\/NOPR|ET,\s*\d+\s*,)/im.test(head)) return { id: 'cdb', c: 0.9, note: 'ANSYS archive commands' };
  if (/^\*(HEADING|NODE|PART|ELEMENT|ASSEMBLY)\b/im.test(head)) return { id: 'abaqus', c: 0.95, note: 'Abaqus keyword lines' };
  if (/^(BEGIN BULK|CEND|SOL\s+\w+|GRID\*?[\s,]|CQUAD4|CTRIA3|CTETRA|CHEXA|ENDDATA|PSHELL|MAT1)/im.test(head)) return { id: 'nastran', c: 0.9, note: 'Nastran bulk-data cards' };
  if (/^\s*(VARIABLES|ZONE)\b/im.test(head)) return { id: 'tecplot', c: 0.9, note: 'Tecplot ASCII keywords (VARIABLES / ZONE)' };
  if (/^v\s+[-+\d.]/m.test(head) && (/^f\s+-?\d/m.test(str(b, 0, Math.min(n, 1 << 20))) || ext === 'obj')) return { id: 'obj', c: 0.9, note: 'OBJ vertex/face statements' };
  if (ext === 'obj' && /^(#|v|vn|vt|o|g|mtllib)\s/m.test(head)) return { id: 'obj', c: 0.6, note: 'OBJ-like statements' };
  if (/^\d+ \d+ \d+ \d+\s*\r?\n/.test(lead) && (ext === 'sat' || /ACIS|SAT/i.test(head.slice(0, 400)))) return { id: 'acis', c: 0.9, note: 'ACIS SAT text header' };

  // Numeric tables: Plot3D (integer block header followed by 3·ni·nj·nk values) versus plain point columns.
  const lines = head.split(/\r?\n/).map((l) => l.trim()).filter((l) => l && !/^[#%!/]/.test(l)).slice(0, 12);
  if (lines.length >= 2) {
    const tk = lines.map((l) => l.split(/[\s,;]+/).filter(Boolean));
    const allInt = (t) => t.length > 0 && t.every(intTok);
    let p3d = false, cells = 0, hdrTok = 0;
    if (tk[0].length === 1 && allInt(tk[0]) && +tk[0][0] >= 1 && +tk[0][0] <= 1000) {
      const nb = +tk[0][0], dimsTok = [];
      for (let i = 1; i < tk.length && dimsTok.length < 3 * nb && allInt(tk[i]); i++) dimsTok.push(...tk[i]);
      if (dimsTok.length >= 3 * nb) { p3d = true; hdrTok = 1 + 3 * nb; for (let k = 0; k < nb; k++) cells += dimsTok[3 * k] * dimsTok[3 * k + 1] * dimsTok[3 * k + 2]; }
    } else if (tk[0].length === 3 && allInt(tk[0]) && tk[0].every((t) => +t >= 1) && tk[0].some((t) => +t > 1)) { p3d = true; hdrTok = 3; cells = tk[0][0] * tk[0][1] * tk[0][2]; }
    if (p3d) {
      // confirm with the token count (exact for small files, estimated from the head for large ones)
      const sample = str(b, 0, Math.min(n, 1 << 20)), tokens = (sample.match(/\S+/g) || []).length * (n / Math.min(n, 1 << 20)) - hdrTok;
      const ratio = tokens / (3 * cells);
      if (ratio > 0.97 && ratio < 1.36) return { id: 'plot3d', c: PLOT3D_EXT.includes(ext) ? 0.9 : 0.75, note: 'integer block-dimension header followed by 3·ni·nj·nk coordinate values (Plot3D formatted grid)' };
    }
    const numeric = tk.filter((t) => t.length >= 2 && t.every(numTok));
    if (numeric.length >= Math.max(2, lines.length - 2)) {
      const cols = numeric[numeric.length - 1].length;
      return { id: 'xyz', c: cols >= 3 ? 0.75 : 0.6, note: cols >= 3 ? `plain numeric table with ${cols} columns (point coordinates)` : 'two-column numeric table (planar coordinates, e.g. an aerofoil section)' };
    }
  }
  return null;
}

/**
 * Identify a file from its name and content. Content signatures decide; the extension is used to choose
 * between containers that share a signature and as a low-confidence fallback.
 * Returns { id, confidence: 0..1, note }; id is 'unknown' when nothing matches.
 */
export function detectFormat(fileName, bytes) {
  const base = String(fileName || '').split(/[\\/]/).pop().toLowerCase();
  const ext = base.includes('.') ? base.split('.').pop() : '';
  if (!(bytes instanceof Uint8Array) || bytes.length === 0) return { id: 'unknown', confidence: 0, note: 'empty file' };
  const extFormats = FORMATS.filter((f) => f.ext.includes(ext));
  let s = null;
  try { s = sniff(bytes, ext, base); } catch { s = null; }
  if (s) {
    const agrees = byId[s.id].ext.includes(ext) || !ext;
    let note = s.note;
    if (!agrees && extFormats.length) note += `; the .${ext} extension usually means ${extFormats.map((f) => f.name).join(' or ')}, but the content decides`;
    else if (!agrees && ext) note += `; unusual extension .${ext}`;
    return { id: s.id, confidence: Math.min(1, agrees ? s.c : s.c * 0.9), note };
  }
  if (['points', 'faces', 'boundary', 'owner', 'neighbour'].includes(base)) return { id: 'openfoam', confidence: 0.35, note: 'file name matches an OpenFOAM polyMesh file, but no FoamFile header was found' };
  if (extFormats.length) {
    const f = extFormats.find((x) => x.support === 'metadata' && x.id !== 'hdf5') && !looksText(bytes) ? extFormats.find((x) => x.support === 'metadata' && x.id !== 'hdf5') : extFormats[0];
    return { id: f.id, confidence: 0.3, note: `content not recognised; guessed from the .${ext} extension only (${extFormats.map((x) => x.name).join(' / ')})` };
  }
  return { id: 'unknown', confidence: 0, note: ext ? `unrecognised content and unknown extension .${ext}` : 'unrecognised content' };
}
