# AeroSuite 26

An integrated, offline-first web application for aircraft engineering analysis: **26 connected suites**
covering aeroplanes, helicopters, rotorcraft, UAVs and electric VTOL aircraft, from aerodynamics and
structures to mission, safety, optimisation, verification and economics — with live data, mesh and
uncertainty studies on every analysis, and decision support.

It is a static site with **no runtime dependencies and no server**. All solvers run in the user's browser
(in a Web Worker), on any computer, Android phone, iPhone or iPad, and keep working in aeroplane mode once
the app has been opened or installed.

## The 26 suites

| # | Suite | # | Suite |
|---|---|---|---|
| 1 | Aerodynamics and CFD | 14 | Landing gear, ground dynamics and impact |
| 2 | Structural mechanics and FEA | 15 | Crashworthiness and impact mechanics |
| 3 | Aeroelasticity and FSI | 16 | Flight control systems and autopilot |
| 4 | Flight dynamics, stability and control | 17 | Avionics, navigation and sensor systems |
| 5 | Aircraft performance and flight envelope | 18 | Hydraulic, pneumatic and mechanical systems |
| 6 | Rotorcraft aerodynamics and aeromechanics | 19 | Electrical and hybrid-electric power systems |
| 7 | Propulsion and engine performance | 20 | Fuel systems and environmental control |
| 8 | Propeller and rotor performance | 21 | Materials and composite mechanics |
| 9 | Fatigue, fracture and damage tolerance | 22 | Safety, reliability and failure analysis |
| 10 | Vibration, modal and rotor dynamics | 23 | Multidisciplinary design analysis and optimisation |
| 11 | Aeroacoustics and noise prediction | 24 | Mission and operational simulation |
| 12 | Thermal engineering and heat transfer | 25 | Verification, validation and uncertainty quantification |
| 13 | Aircraft icing and ice protection | 26 | Economics and techno-economic analysis |

Each suite holds several analyses. Every analysis has: inputs with units, initial and boundary conditions;
defaults linked to the shared case and to upstream suites; a solver; KPIs, charts, contours and tables;
stated models, assumptions and validity warnings; recommendations; and — through the common study engine —
mesh/step convergence with Richardson extrapolation and GCI, sensitivity ranking, parameter sweeps,
uncertainty propagation, calibration and validation against user data, and code-verification benchmarks.

## Run it

```bash
node tools/serve.mjs 8080        # then open http://localhost:8080/
```

Any static file server works. There is nothing to compile for development.

```bash
npm install                      # dev tools only: esbuild (single-file build) and playwright-core (UI test)
npm test                         # all suites on all aircraft presets + verification benchmarks + geometry readers
npm run test:full                # the same plus convergence studies and the coupled chain, verbose
npm run test:ui                  # headless-browser walk through every page (needs a Chromium for Playwright)
node tests/licences.mjs          # licence registry, commercial-mode host policy, attribution and notices (no network)
npm run build                    # stamp the offline asset list/version and rebuild standalone.html
```

## How it is organised

```
index.html, css/, sw.js, manifest.webmanifest     app shell, styles, offline engine, install manifest
js/app.js                                         navigation, routing, theme, install, search
js/core/numerics.js, atmosphere.js                shared numerical kernel and standard atmosphere
js/core/case.js                                   the shared aircraft case and presets
js/core/registry.js, studies.js, jobs.js          suite index, generic study engine, job execution
js/core/runner.js, worker.js                      Web Worker execution with main-thread fallback
js/core/live.js                                   live-data connectors with on-device caching
config.json                                       deployment configuration (commercial mode, keys, providers)
js/data/licences.json, COMPLIANCE.md              audited registry of data sources, terms and components; compliance pack
js/core/gridwx.js, places.js                      forecast-grid reader; bundled places, climate and terrain
tools/grid.mjs, grib2.mjs, bunzip2.mjs            forecast grids from NOAA and DWD model output (snapshot job)
LICENSE, THIRD_PARTY_NOTICES.md                   proprietary licence; third-party components and notices
js/core/geometry/                                 geometry and mesh importers, healing, quality, sectioning
js/suites/sNN-*.js                                the 26 suites (pure computation, no DOM)
js/ui/                                            charts, 3-D viewer, pages
js/data/                                          materials database, engineering specification
tests/                                            Node and browser test harnesses
docs/SUITE_CONTRACT.md                            the contract every suite module follows
```

The suites are connected through a shared case (entered once) and a data bus of published outputs
(`docs/SUITE_CONTRACT.md` lists the keys). *Integrated run* executes all suites in dependency order.

## Live data

The same features are available to every deployment — commercial or not, with or without API keys — because
each one rests on a **keyless source whose licence explicitly allows commercial reuse**. Services limited to
non-commercial use (the Open-Meteo free API, GitHub search) are only an optional refinement on top when
`"commercial": false`. Data reaches the app by three routes:

1. **Bundled with the app (no network).**
   - Airports and runways: [OurAirports](https://ourairports.com/data/) (public domain), `js/data/airports/`,
     built by `node tools/fetch-data.mjs`.
   - Place names: GeoNames (CC BY 4.0), `js/data/places/`, built by `node tools/fetch-data.mjs --places`:
     34 000 cities above 15 000 inhabitants (1.2 MB) and, read only when a search needs them, 35 600 places of
     5 000 to 15 000 inhabitants in one file per first letter (`small/`, 1.4 MB). Place search finds cities,
     small towns and airports offline. (GeoNames' next list, places above 1 000, is about 5 MB: not bundled.)
   - Terrain: NOAA ETOPO 2022 (CC0), mean elevation of 0.5° cells, `js/data/terrain/`, built by
     `node tools/fetch-data.mjs --terrain` — route profiles and last-resort site elevation.
   - Design temperatures: hot day (99th percentile of daily maximum), cold day (1st percentile of daily
     minimum) and strong wind from ten years (2015–2024; wind 2022–2024) of the NCEP-DOE Reanalysis 2 (NOAA
     PSL, public domain) on a 2.5° grid, `js/data/climate/`, built once by `node tools/fetch-climate.mjs`;
     corrected to the site elevation with 6.5 K/km. A coastal reanalysis cell is partly sea and far
     windier than the land beside it, so the strong-wind field for land sites replaces sea cells by nearby
     land values (land–sea mask of the reanalysis); the unmodified field serves sites at sea. Over
     mid-latitude land the reanalysis itself is windier than finer ones: the strong-wind value is a
     conservative regional figure.
   - National grid emission factors for 208 economies: Ember yearly electricity data (CC BY 4.0),
     `js/data/grid-factors.json`, built by `node tools/fetch-grid-factors.mjs`.
   - Open-source tools per suite: a curated list, `js/data/tools.json`.
   - Read through `js/core/airports.js` and `js/core/places.js`.
2. **Cloud snapshot** (`.github/workflows/snapshot.yml` runs `node tools/snapshot.mjs` every three hours and
   republishes the site; nothing is committed). Fetched server-side from the primary publishers:
   - `data/grid/` — **global forecast grids** (`tools/grid.mjs`): NOAA GFS weather (surface fields at 1°
     and wind, temperature and height at eleven pressure levels at 2.5° for the three valid times nearest to
     now; 2.5° and 5° for the two later ones; freezing level, gust, visibility, cloud, precipitation), DWD
     GWAM sea state (2° cells), NOAA GEFS-Aerosols PM2.5/PM10/dust/optical depth (2.5°). Each holds valid
     times about 24 h ahead; the app interpolates in space and time (`js/core/gridwx.js`), loads only the
     two slices around the present, and keeps the files for offline use. Files are difference-coded and
     gzipped (about 2.4 MB in all; the browser's own gzip decompression unpacks them). The GRIB2 and bzip2 readers are in `tools/grib2.mjs` and
     `tools/bunzip2.mjs` (no dependencies).
   - `data/snapshot.json` — exchange rates (ECB), Brent and U.S. Gulf Coast jet fuel (US EIA), policy rates
     (Federal Reserve Board H.15, New York Fed as fall-back; ECB), space weather, World Bank world
     aggregates, and the **carbon panel**: EU ETS (the European Commission's price of CBAM certificates —
     the average EU allowance auction clearing price of the period, CC BY 4.0), UK ETS (UK ETS Authority
     annual determination, OGL v3.0) and the California–Québec joint auction (Gouvernement du Québec, CC BY
     4.0). The case takes the market of the site's country (the EU ETS price elsewhere) unless the user picks
     another on the Live data page; the operator can add its own price (`carbonPriceUrl`).
3. **Straight from each user's browser**, where it adds something: a point forecast from MET Norway or the
   latest US NWS observation (both commercial-use sources; their gaps are filled from the NOAA grid), NOAA
   SWPC, the ECB Data Portal, World Bank, OpenAlex; and, in non-commercial deployments only, Open-Meteo
   (finer weather, ERA5 year at the site, gases, smaller places), the open oil-prices dataset, the public
   Overpass servers, GitHub search and the GB Carbon Intensity API.

Each feed records its source, address, fetch time, success and the licence-registry ids of the sources it was
built from; a failed feed keeps its last good value. FRED, the EEX auction report (unless the operator holds
EEX's approval and lists `eex-auction` under `accept`) and NASA POWER are not used in any mode. No free,
keyless price series for sustainable aviation fuel was found, so SAF stays a stated multiple of jet fuel.

Feeds are cached on the device with a time-to-live and reused offline. `tests/data.mjs` checks the airport
database, the snapshot schema and the unit conversions; `tests/licences.mjs` checks the licence registry, the
commercial-mode host policy and **feature parity** between the modes.

## Data sources, licences and commercial use

Every external data source, bundled dataset and third-party component is recorded in
**`js/data/licences.json`**: what it is used for, its licence or terms, the page and the sentence the
conclusion rests on (read on the date given), its class, its obligations, the attribution text and where it is
shown, whether a commercial deployment uses it and what replaces it when it does not. The app shows this
registry under *Live data → Data sources and licences* and on the *Install, offline & about* page, and prints
the required attribution next to the data it belongs to. `node tools/compliance.mjs` writes the same record
as **`COMPLIANCE.md`** (per mode: every source contacted or shipped, licence, obligation, where the
attribution is shown, the sentence relied on, and the questions still open); the About page shows it as a
printable *Compliance pack*. The registry is a compliance record, not legal advice.

| Class | Meaning | Sources |
|---|---|---|
| `commercial-ok` | terms explicitly allow commercial use; the foundation of both modes | NOAA GFS, NOAA GEFS-Aerosols, DWD wave model, NCEP-DOE Reanalysis 2, NOAA ETOPO, GeoNames, OurAirports, MET Norway, US NWS, NOAA SWPC, ECB statistics, Frankfurter, World Bank WDI, US EIA, Federal Reserve Board, New York Fed, European Commission (CBAM certificate price), UK ETS Authority, Gouvernement du Québec, Ember, OpenAlex, OpenStreetMap data (ODbL), EASA-published databases, FAA database, US Government citations |
| `commercial-ok`, kept out of commercial mode by policy | licence allows it, but a term needs the operator's own acceptance | GB Carbon Intensity API (indemnity clause): list it under `accept` to use it |
| `commercial-with-key` | only with a paid key, own server or own licence | Open-Meteo customer API, an Overpass server of the operator, an operator-supplied carbon price, the solver bridge (the user's own GitHub token) |
| `unclear` | no explicit grant found; never used in commercial mode | GitHub repository search, the GitHub-hosted oil dataset (NASA POWER: removed) |
| `not-allowed` | excluded without the provider's consent | Open-Meteo free API (non-commercial refinement only), EEX auction report (off in every mode), FRED (removed) |

### Deployment configuration (`config.json`)

`config.json` at the site root is read when the app starts (and by `tools/snapshot.mjs` and
`tools/build.mjs`). It is stored for offline use; if it is missing, the defaults below apply. A device
that has once read `"commercial": true` stays in commercial mode while offline. Run `npm run build`
after changing it: the build writes the content-security policy from it.

| Key | Default | Meaning |
|---|---|---|
| `commercial` | `false` | `true` when the app is sold, monetised or offered as part of a commercial product or service. Switches every feed to the compliant set: a source is contacted only if the registry classes it `commercial-ok`, or the key or address it needs is configured below, or its id is listed under `accept`. The snapshot job then skips the EEX report and the GitHub-hosted oil dataset, the build removes the non-commercial hosts from the content-security policy, and values cached under non-commercial terms are not reused. |
| `openMeteoApiKey` | `""` | Key of an Open-Meteo API subscription. When set, weather, winds aloft, reanalysis, air quality, sea state, place search and terrain use the customer hosts (`customer-api.open-meteo.com` and siblings) with `&apikey=…`. The reanalysis feed needs the Professional plan. The key is sent from each user's browser and is therefore visible to users. |
| `weatherProvider` | `"auto"` | `"auto"`, `"open-meteo"`, `"met-norway"`, `"nws"`, `"noaa-grid"` or `"none"`. `auto` is Open-Meteo, MET Norway, NOAA grid in a non-commercial deployment; in a commercial one Open-Meteo (only with a key), MET Norway, US NWS, NOAA grid. A named provider is tried first and the others stay as fall-backs; `open-meteo` without a key is ignored in commercial mode. MET Norway and NWS replies are completed from the NOAA grid. |
| `carbonPriceUrl` | `""` | Address (https, or a path on the site) of a JSON document with a carbon price the operator is licensed to publish: `{ "date": "2026-10-08", "price": 85.07, "currency": "EUR", "market": "EU ETS (EUA Dec-26)", "source": "…", "attribution": "…" }`. When set it is added to the carbon panel and used for the case in every mode (unless the user picks another market). The other host must send CORS headers; the build adds it to the content-security policy. |
| `attribution` | `true` | Show the providers' attribution lines next to the data. Forced to `true` in commercial mode, where the notices are licence conditions; the registry tables are always shown. |
| `overpassUrl` | `""` | Address of an Overpass API server the operator runs or rents. Needed for the optional OpenStreetMap detail in commercial mode, because the public servers ask commercial users to use their own. |
| `openAlexApiKey` | `""` | Optional OpenAlex API key for the literature lists (the keyless allowance is small). Visible to users. |
| `accept` | `[]` | Registry ids the operator has cleared or licensed itself, for example `["gb-carbon-intensity"]` after accepting that API's terms, or `["eex-auction"]` with EEX's written approval. Each listed source is then used in commercial mode as in non-commercial mode. `["nasa-geos-cf"]` adds NO₂ and ozone from NASA's GEOS-CF forecast (snapshot job and app) in any mode. |

### Provider of each feature in each mode

| Feature | Commercial (no keys) | Non-commercial | Difference |
|---|---|---|---|
| Surface weather | MET Norway point forecast (US NWS observation as second choice), completed from the NOAA GFS grid; the grid alone when neither answers or offline | Open-Meteo; then the same chain | Open-Meteo is a finer model blend |
| Winds and temperatures aloft, freezing level, gusts, visibility | NOAA GFS grid (11 levels at 2.5°; surface 1°) | Open-Meteo; NOAA GFS grid as fall-back | the grid is coarser than Open-Meteo's models |
| Design temperatures | NCEP-DOE Reanalysis 2 statistics, 2015–2024, bundled (offline) | Open-Meteo ERA5 at the site, last 12 months; bundled statistics as fall-back | bundled values are multi-year but regional (2.5° grid) |
| Sea state | DWD wave model grid (2° cells) | Open-Meteo marine; DWD grid as fall-back | coarser near coasts |
| Air quality and dust | NOAA GEFS-Aerosols grid: PM2.5, PM10, dust, optical depth | Open-Meteo / CAMS; NOAA grid as fall-back | **NO₂ and ozone**: no keyless source with an explicit commercial-reuse grant publishes them worldwide (CAMS needs an account key; NOAA's global model has no gases). They are shown with Open-Meteo (free tier when non-commercial, or `openMeteoApiKey`), or from NASA GEOS-CF once the operator lists `nasa-geos-cf` under `accept` — NASA gives no explicit reuse grant for that experimental forecast, so it is off by default; the card says so |
| Place search | bundled GeoNames places (5 000 inhabitants and more) and OurAirports airports (offline) | the same, plus smaller places from the Open-Meteo geocoder | places under 5 000 inhabitants |
| Elevation | airport within 10 km, else city within 15 km, else ETOPO 0.5° terrain (all bundled); MET Norway model height with its forecast | Open-Meteo point elevation; bundled data as fall-back | point value versus nearest feature / cell mean |
| Aerodromes and runways | OurAirports (bundled) | the same; optional extra detail from OpenStreetMap through the public Overpass servers | nothing needed for runway data; OSM detail needs `overpassUrl` in commercial mode |
| Carbon price | official panel: EU ETS, UK ETS, California–Québec; operator's own price if configured | the same | none |
| Exchange rates, Brent, jet fuel, policy rates, space weather, inflation, literature | the same primary sources in both modes | the same | none (Brent is read by the browser from an open dataset when non-commercial, from the snapshot otherwise) |
| Grid emission factor | national average of the site's country (Ember, bundled) | the same; in Great Britain the current half-hour from the NESO API | live GB value only |
| Open-source tools per suite | curated list (bundled) | the same, plus GitHub search results | extra search results |

With `openMeteoApiKey` a commercial deployment uses Open-Meteo's customer hosts first, exactly as the
non-commercial column, and loses nothing. `tests/licences.mjs` asserts the parity above for Lagos, London,
Denver and Sydney with every online provider unreachable.

Weather providers return one normalised object (`js/core/live.js`, `WEATHER`): time, elevation, temperature,
humidity, station and sea-level pressure, wind, gust, precipitation, cloud, weather code, visibility,
freezing level and winds aloft. Fields a point provider does not have are taken from the NOAA grid and named
in `filled`; the location and Live data pages say so. MET Norway requires each client to identify itself,
which a page does through its `Origin` header: the single-file copy opened from disk therefore skips it.

**Single-file copy.** `standalone.html` cannot read the data files beside the site, so `tools/build.mjs` embeds
an essential subset (about 1.2 MB): large and medium airports with their runways, cities of 100 000
inhabitants or more, the design-temperature grid, the grid emission factors, the tool lists and the cloud
snapshot of the build (carbon panel, exchange rates, fuel prices, policy rates). Smaller airports and
places, terrain and the forecast grids need the served app; the Live data page of the single-file copy
says what it carries.

### Before selling: operator checklist

1. Set `"commercial": true` in `config.json`, run `npm run build`, deploy, and run `node tests/licences.mjs`.
   No key or contract is needed for any feature.
2. Keep the scheduled snapshot job running (it publishes the forecast grids and the carbon panel). Serve the
   app from its own domain; above low traffic put a caching proxy in front of `api.met.no`, as its terms ask —
   or set `"weatherProvider": "noaa-grid"` to use the grid only.
3. Optional refinements: an Open-Meteo subscription (`openMeteoApiKey`) for the finer weather, the site's own
   reanalysis year, NO₂/ozone and small places; a contracted carbon market price at `carbonPriceUrl`; an own
   Overpass server (`overpassUrl`).
4. Ship `LICENSE`, `THIRD_PARTY_NOTICES.md`, `COMPLIANCE.md` and the licence texts under `js/vendor/` with
   every copy (the build and the publish workflow do), and offer the source of the LGPL-licensed OpenCASCADE
   kernel as `THIRD_PARTY_NOTICES.md` describes.
5. Have counsel answer the questions listed at the top of `COMPLIANCE.md`.

## Licence

The application is proprietary: Copyright (c) 2026 Samuel Akosa Onyejekwe. All rights reserved. See
`LICENSE`. Third-party components and datasets remain under their own licences: see
`THIRD_PARTY_NOTICES.md` and `js/data/licences.json`.

## Geometry and meshes

The app imports; it does not draw. `js/core/geometry/formats.js` holds a capability profile for every
format: read in full, read to a documented subset, or recognised by signature with a stated conversion
route. Units are never guessed, healing is explicit and logged, and each reader is tested individually
(`tests/geometry.mjs`).

## Install and offline

The service worker stores the application on first visit and serves it from the device afterwards.
Browsers offer installation (desktop Chrome/Edge, Android, iOS “Add to Home Screen”, macOS “Add to Dock”).
`standalone.html` is the whole application in one file that opens from disk with no server.

## Deploy, with redundancy

The site is plain static files, so the same commit can be published to several independent hosts:

- **GitHub Pages** — `.github/workflows/deploy.yml` tests, builds and deploys on every push to `main`
  (repository Settings → Pages → Source: GitHub Actions).
- **Cloudflare Pages / Netlify / Vercel** — connect the repository; no build command, publish directory `.`
  (`_headers`, `netlify.toml` and `vercel.json` are included).
- **GitLab Pages / Codeberg Pages** — push the repository to a second remote (`.gitlab-ci.yml` included).

List the resulting addresses in `mirrors.json`, run `npm run build`, and commit: the *Install, offline &
about* page then shows each mirror with a live health check. Once hosts have deployed a build they keep
serving it even if the source repository host is unavailable, and installed copies need no host at all.

## Engineering credibility

Results state their model, assumptions and validity limits and distinguish calculated values, empirical
estimates and assumed inputs. Verification benchmarks demonstrate that equations are solved correctly; they
do not demonstrate that a model represents a particular aircraft. Each analysis states the resolution it
runs at and what that resolution supports: where a browser-sized grid gives a trend rather than a design
value, the result says so and the high-fidelity bridge writes the same case for SU2, OpenFOAM or CalculiX.
Default material properties, failure rates, cost coefficients and correlations carry their source beside
the input; values marked as estimated or assumed are to be replaced with programme data. Nothing produced here is, on its own, airworthiness evidence.
