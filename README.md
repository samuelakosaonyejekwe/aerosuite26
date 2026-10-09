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

Data reaches the app by three independent routes, so freshness never depends on a particular machine:

1. **Bundled with the app (no network).** Airports and runways worldwide come from the public-domain
   [OurAirports](https://ourairports.com/data/) datasets, packed by `node tools/fetch-data.mjs` into
   10° × 10° tiles under `js/data/airports/` (dataset date in `js/data/airports/index.json`) and read
   through `js/core/airports.js` (`nearest`, `findByCode`, `search`). Re-run the tool to update. The
   location page uses this database first; OpenStreetMap (© OpenStreetMap contributors, ODbL, through
   the public Overpass servers) is an optional enrichment on request.
2. **Straight from each user's browser** to public, CORS-enabled providers: Open-Meteo (weather, winds
   aloft, reanalysis, air quality, sea state, geocoding, terrain), NOAA SWPC (space weather),
   Frankfurter/ECB (exchange rates), World Bank (inflation, lending rates), the open oil-prices dataset
   of EIA spot series (Brent), the GB Carbon Intensity API, OpenAlex (literature) and GitHub
   (open-source tools).
3. **Cloud snapshot.** `.github/workflows/snapshot.yml` runs `node tools/snapshot.mjs` every three
   hours on GitHub's runners and republishes the site with a fresh `data/snapshot.json` (artifact
   deploy; nothing is committed). It holds every location-independent feed, so the first screen after
   opening is current even on browsers without background sync, and it is the only route to providers
   that send no CORS headers:
   - jet fuel — U.S. Gulf Coast kerosene-type jet fuel spot price, US EIA (public domain), read from the
     FRED CSV export of series `DJFUELUSGULF` with the EIA history page as fall-back; USD/gal converted
     to USD/kg at 0.804 kg/L;
   - carbon — clearing price of the latest EU ETS allowance (EUA) primary auction, from the public
     auction report of EEX, the EU common auction platform, in EUR/t CO₂, converted to USD with the
     ECB reference rate;
   - policy rates — effective federal funds rate (Federal Reserve, FRED `DFF`) and ECB deposit
     facility rate (ECB Data Portal).
   Each feed records its source, address, fetch time and success; a failed feed keeps its last good
   value. Other hosts read the newest snapshot from the addresses in `mirrors.json`. No free, keyless
   price series for sustainable aviation fuel was found, so SAF stays a stated multiple of jet fuel.

Feeds are cached on the device with a time-to-live and reused offline. `tests/data.mjs` checks the
airport database, the snapshot schema and the unit conversions (`--browser` adds the offline runway
look-up in a real browser).

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
do not demonstrate that a model represents a particular aircraft. Models that cannot be solved credibly in
a browser (three-dimensional RANS/LES/DNS, full-aircraft explicit crash FE, coupled high-fidelity CFD–CSD
and similar) are listed per suite as handed off to external solvers. Default material properties, failure
rates, cost coefficients and correlations are typical or illustrative values to be replaced with sourced
data. Nothing produced here is, on its own, airworthiness evidence.
