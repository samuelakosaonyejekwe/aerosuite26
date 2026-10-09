# Third-party notices — AeroSuite 26

AeroSuite 26 is proprietary software (see `LICENSE`). It is distributed together with the third-party
components and datasets listed here. Each stays under its own licence; the licence texts are shipped with the
site at the paths given. The machine-readable form of this list, with the page and sentence each conclusion
rests on, is `js/data/licences.json`; the app shows it under *Live data → Data sources and licences* and on
the *Install, offline & about* page. Licence terms were read on 2026-10-09.

## Software components shipped with the application

All of them are loaded as **separate files at run time**, only when a file of the matching format is opened.
None is compiled into the application's own code, and none is part of `standalone.html`.

| Component | Version | Licence | Files | Licence text |
|---|---|---|---|---|
| occt-import-js | 0.0.23 | GNU LGPL 2.1 | `js/vendor/occt/occt-import-js.js`, `js/vendor/occt/occt-import-js.wasm` | `js/vendor/occt/LICENSE.occt-import-js.txt` |
| Open CASCADE Technology (compiled into `occt-import-js.wasm`) | submodule commit `d2abb6d844231cb8f29be6894440874a4700e4a5` of occt-import-js 0.0.23 | GNU LGPL 2.1 with the Open CASCADE exception 1.0 | `js/vendor/occt/occt-import-js.wasm` | `js/vendor/occt/LICENSE.occt.txt` |
| h5wasm, with the HDF5 library | as vendored on 2026-10-09 (the files do not record a version number) | NIST software statement; HDF5 BSD-style licence (The HDF Group, University of Illinois) | `js/vendor/h5wasm/hdf5_hl.js`, `js/vendor/h5wasm/hdf5_util.js` | `js/vendor/h5wasm/LICENSE.txt` |
| laz-perf | 0.0.7 | Apache License 2.0 | `js/vendor/lazperf/laz-perf.js`, `js/vendor/lazperf/laz-perf.wasm` | `js/vendor/lazperf/LICENSE.txt` |

### Modifications

- **`js/vendor/occt/occt-import-js.js`** (changed 2026-10-09; upstream file `dist/occt-import-js.js`, 0.0.23):
  two additions so that the file loads as an ES module — a module-scope declaration `var process;` near the
  top (so that the loader, not the glue script, supplies the `.wasm` bytes in every runtime) and a final line
  `export default occtimportjs;`. Nothing else is changed. `occt-import-js.wasm` is the unmodified upstream
  build.
- **`js/vendor/lazperf/laz-perf.js`** (changed 2026-10-09; upstream file `lib/laz-perf.js`, 0.0.7): a final
  line `export default createLazPerf;` was added. Nothing else is changed. `laz-perf.wasm` is unmodified.
- **h5wasm**: no modification.

Each modified file states its change in its first lines.

### OpenCASCADE kernel: LGPL 2.1

This application uses the library occt-import-js, which contains Open CASCADE Technology. Both are free
software under the GNU Lesser General Public License, version 2.1; Open CASCADE Technology additionally
carries the *Open CASCADE exception (version 1.0)*, reproduced here because it is not part of the licence
file:

> The object code (i.e. not a source) form of a "work that uses the Library" can incorporate material from a
> header file that is part of the Library. As a special exception to the GNU Lesser General Public License
> version 2.1, you may distribute such object code incorporating material from header files provided with the
> Open CASCADE Technology libraries (including code of CDL generic classes) under terms of your choice,
> provided that you give prominent notice in supporting documentation to this code that it makes use of or is
> based on facilities provided by the Open CASCADE Technology software.

AeroSuite 26 makes use of facilities provided by the Open CASCADE Technology software. The application's own
code contains no part of the library: it calls three functions the kernel exports (`ReadStepFile`,
`ReadIgesFile`, `ReadBrepFile`) through `js/core/geometry/parsers-wasm.js`.

**Source code.** The corresponding source of the library as shipped is:

- occt-import-js 0.0.23 — <https://github.com/kovacsv/occt-import-js/tree/0.0.23> (with its `CMakeLists.txt` and `tools/` build files);
- Open CASCADE Technology — <https://git.dev.opencascade.org/repos/occt.git>, the commit named above
  (the `occt` submodule of that tag);
- the modified glue script — `js/vendor/occt/occt-import-js.js` is itself the source form, with the two
  additions described above.

Whoever distributes AeroSuite 26 must make this source available to recipients from the same place as the
application, or accompany it with a written offer, valid for at least three years, to supply it (LGPL 2.1,
sections 4 and 6). A link to the upstream repositories alone does not discharge that duty if they disappear.

**Replacing the library.** The kernel is two files, `js/vendor/occt/occt-import-js.js` and
`js/vendor/occt/occt-import-js.wasm`, fetched when a STEP, IGES or BREP file is opened. To use a modified
version:

1. build occt-import-js from the source above, with your changes, to obtain `occt-import-js.js` and
   `occt-import-js.wasm`;
2. replace the two files under `js/vendor/occt/` in your copy of the site with the files you built, keeping
   the names. The glue script must have a default export that is the module factory (add the line
   `export default occtimportjs;` and the `var process;` declaration described above);
3. reload the application (if a service worker holds the old files, choose *Check for an update* or clear the
   site data).

The interface the application relies on is the one the library documents: the module factory, called with
`{ wasmBinary }`, resolves to an object with `ReadStepFile`, `ReadIgesFile` and `ReadBrepFile`, each taking
the file content as a `Uint8Array` and a triangulation-parameter object.
The proprietary licence of the application permits this replacement and the reverse engineering needed to
debug it (see `LICENSE`).

### h5wasm and HDF5

h5wasm was developed at the National Institute of Standards and Technology (NIST), which is acknowledged as
its source; its notice and the HDF5 copyright notice, conditions and disclaimer are kept in full in
`js/vendor/h5wasm/LICENSE.txt`.

### laz-perf

Apache License 2.0; the licence text is in `js/vendor/lazperf/LICENSE.txt`. The vendored copy came without a
NOTICE file.

## Development tools (not distributed with the application)

| Tool | Version | Licence | Use |
|---|---|---|---|
| esbuild | 0.28.x | MIT | builds `standalone.html` (`tools/build.mjs`) |
| playwright-core | 1.49.1 | Apache License 2.0 | browser tests |

They are installed by `npm install` for development only and are not part of the deployed site. esbuild's
output (`standalone.html`) contains only the application's own code.

## Datasets bundled with the application

| Dataset | Files | Licence or terms | Notice |
|---|---|---|---|
| OurAirports airports and runways | `js/data/airports/` | Public domain | Airports and runways: OurAirports (public domain). The data comes with no guarantee of accuracy or fitness for use. |
| ICAO Aircraft Engine Emissions Databank, issue 32, as published by EASA (reduced extract) | `js/data/ref/engine-emissions.json`, `.js` | EASA copyright notice | © European Union Aviation Safety Agency. Reproduction is authorised, provided the source is acknowledged. |
| EASA certification noise levels (TCDSN databases, reduced extract) | `js/data/ref/cert-noise.json`, `.js` | EASA copyright notice | © European Union Aviation Safety Agency. Reproduction is authorised, provided the source is acknowledged. |
| FAA Aircraft Characteristics Database (extract) | `js/data/ref/aircraft-characteristics.json` | Work of the United States Government (17 U.S.C. § 105) | Source: U.S. Federal Aviation Administration. |
| Source registry with short quotations | `js/data/sources.json` | US Government works (eCFR, FAA Advisory Circulars, NASA reports, MIL handbooks): not subject to copyright in the United States. Other publishers: short quotations with the source named, each publisher's copyright reserved | Each entry names its document, section and address. |

## Live data sources

The providers contacted at run time, their terms, the attribution each requires and whether a commercial
deployment uses them are listed in `js/data/licences.json` and shown in the app. Required attribution
statements:

- Weather data by Open-Meteo.com (CC BY 4.0) — <https://open-meteo.com/>
- Weather: based on data from MET Norway (NLOD 2.0 / CC BY 4.0) — <https://api.met.no/doc/License>
- Weather observation: US National Weather Service (NOAA), public domain
- Site climate statistics: NASA Langley Research Center POWER project (only where a deployment enables it)
- © OpenStreetMap contributors, available under the Open Database License — <https://www.openstreetmap.org/copyright>
- Space weather: NOAA Space Weather Prediction Center (public domain)
- Source: ECB statistics (euro foreign exchange reference rates and key interest rates); rates per US dollar are own calculations
- Source: World Bank, World Development Indicators (CC BY 4.0)
- Source: U.S. Energy Information Administration
- Source: Board of Governors of the Federal Reserve System, H.15; Federal Reserve Bank of New York (effective federal funds rate)
- UK ETS carbon price: UK ETS Authority, GOV.UK. Contains public sector information licensed under the Open Government Licence v3.0.
- EU allowance auction price: EEX (European Energy Exchange) — non-commercial deployments only
- Carbon intensity: National Energy System Operator (NESO) Carbon Intensity API (CC BY 4.0)
- Literature index: OpenAlex (CC0)

None of these providers endorses AeroSuite 26.

## Where to get the source code of the LGPL components

The complete corresponding source code of occt-import-js 0.0.23 and of the Open CASCADE Technology revision it was built
from (commit d2abb6d844231cb8f29be6894440874a4700e4a5), together with the modified loader script, is published with this
application at:

https://github.com/samuelakosaonyejekwe/aerosuite26/releases/tag/third-party-sources

The loader script `js/vendor/occt/occt-import-js.js` was modified on 2026-10-09 (a module-scope `process` declaration and a
default export were added). The library is shipped as separate run-time files and may be replaced by a rebuilt version.
