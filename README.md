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
js/data/licences.json                             audited registry of data sources, terms and components
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

Data reaches the app by three independent routes, so freshness never depends on a particular machine:

1. **Bundled with the app (no network).** Airports and runways worldwide come from the public-domain
   [OurAirports](https://ourairports.com/data/) datasets, packed by `node tools/fetch-data.mjs` into
   10° × 10° tiles under `js/data/airports/` (dataset date in `js/data/airports/index.json`) and read
   through `js/core/airports.js` (`nearest`, `findByCode`, `search`). Re-run the tool to update. The
   location page uses this database first; OpenStreetMap (© OpenStreetMap contributors, ODbL) is an
   optional enrichment on request.
2. **Straight from each user's browser** to public, CORS-enabled providers: weather (Open-Meteo, MET
   Norway, US National Weather Service — see the provider table below), NOAA SWPC (space weather), the
   European Central Bank Data Portal (euro reference exchange rates; Frankfurter as fall-back), World Bank
   (inflation, lending rates), the GB Carbon Intensity API and OpenAlex (literature); in non-commercial
   deployments also the open oil-prices dataset (Brent), the public Overpass servers and GitHub
   repository search.
3. **Cloud snapshot.** `.github/workflows/snapshot.yml` runs `node tools/snapshot.mjs` every three
   hours on GitHub's runners and republishes the site with a fresh `data/snapshot.json` (artifact
   deploy; nothing is committed). It holds every location-independent feed, so the first screen after
   opening is current even on browsers without background sync, and it is the only route to providers
   that send no CORS headers. Every feed goes to the primary publisher:
   - Brent and jet fuel — Europe Brent spot price and U.S. Gulf Coast kerosene-type jet fuel spot price,
     from the US EIA history pages (US Government, public domain); USD/gal converted to USD/kg at
     0.804 kg/L;
   - policy rates — effective federal funds rate from the Federal Reserve Board's H.15 data download
     (New York Fed Markets Data API as fall-back) and the ECB deposit facility rate (ECB Data Portal);
   - carbon — non-commercial deployments: clearing price of the latest EU ETS allowance (EUA) primary
     auction from the EEX auction report, in EUR/t CO₂; commercial deployments: the UK ETS carbon price
     determined by the UK ETS Authority for the scheme year (GOV.UK, Open Government Licence v3.0), in
     GBP/t CO₂e — a UK figure fixed once a year, labelled as such. Both are converted to USD with the ECB
     reference rate.
   Each feed records its source, address, fetch time, success and the licence-registry ids of the sources
   it was built from; a failed feed keeps its last good value. Other hosts read the newest snapshot from
   the addresses in `mirrors.json`. FRED is not used: its terms exclude commercial use and caching, and
   the series are available from their publishers. No free, keyless price series for sustainable
   aviation fuel was found, so SAF stays a stated multiple of jet fuel.

Feeds are cached on the device with a time-to-live and reused offline. `tests/data.mjs` checks the
airport database, the snapshot schema and the unit conversions (`--browser` adds the offline runway
look-up in a real browser); `tests/licences.mjs` checks the licence registry and the commercial mode.

## Data sources, licences and commercial use

Every external data source, bundled dataset and third-party component is recorded in
**`js/data/licences.json`**: what it is used for, its licence or terms, the page and the sentence the
conclusion rests on (read on the date given), its class, its obligations, the attribution text, whether a
commercial deployment uses it and what replaces it when it does not. The app shows this registry under
*Live data → Data sources and licences* and on the *Install, offline & about* page, and prints the
required attribution next to the data it belongs to. The registry is a compliance record, not legal
advice.

| Class | Meaning | Sources |
|---|---|---|
| `commercial-ok` | terms explicitly allow commercial use | MET Norway, US NWS, OurAirports, OpenStreetMap data (ODbL), NOAA SWPC, ECB statistics, Frankfurter, World Bank WDI, US EIA, Federal Reserve Board, New York Fed, UK ETS Authority (GOV.UK), GB Carbon Intensity API, OpenAlex, EASA-published databases, FAA database, US Government citations |
| `commercial-with-key` | only with a paid key, own server or own licence | Open-Meteo customer API, an Overpass server of the operator, an operator-supplied carbon price, the solver bridge (the user's own GitHub token) |
| `unclear` | no explicit grant found; off in commercial mode | NASA POWER, GitHub repository search, the GitHub-hosted oil dataset, short quotations from non-government publishers in `js/data/sources.json` |
| `not-allowed` | excluded without the provider's consent; off in commercial mode | Open-Meteo free API, EEX auction report, FRED (removed in every mode) |

### Deployment configuration (`config.json`)

`config.json` at the site root is read when the app starts (and by `tools/snapshot.mjs` and
`tools/build.mjs`). It is stored for offline use; if it is missing, the defaults below apply. A device
that has once read `"commercial": true` stays in commercial mode while offline. Run `npm run build`
after changing it: the build writes the content-security policy from it.

| Key | Default | Meaning |
|---|---|---|
| `commercial` | `false` | `true` when the app is sold, monetised or offered as part of a commercial product or service. Switches every feed to the compliant set: a source is contacted only if the registry classes it `commercial-ok`, or the key or address it needs is configured below, or its id is listed under `accept`. The snapshot job then skips the EEX report and the GitHub-hosted oil dataset, the build removes the non-commercial hosts from the content-security policy, and values cached under non-commercial terms are not reused. |
| `openMeteoApiKey` | `""` | Key of an Open-Meteo API subscription. When set, weather, winds aloft, reanalysis, air quality, sea state, place search and terrain use the customer hosts (`customer-api.open-meteo.com` and siblings) with `&apikey=…`. The reanalysis feed needs the Professional plan. The key is sent from each user's browser and is therefore visible to users. |
| `weatherProvider` | `"auto"` | `"auto"`, `"open-meteo"`, `"met-norway"`, `"nws"` or `"none"`. `auto` is Open-Meteo then MET Norway in a non-commercial deployment; in a commercial one Open-Meteo (only with a key), then MET Norway, then the US NWS. A named provider is tried first and the others stay as fall-backs; `open-meteo` without a key is ignored in commercial mode. |
| `carbonPriceUrl` | `""` | Address (https, or a path on the site) of a JSON document with a carbon price the operator is licensed to publish: `{ "date": "2026-10-08", "price": 85.07, "currency": "EUR", "market": "EU ETS (EUA Dec-26)", "source": "…", "attribution": "…" }`. When set it replaces the built-in carbon feed in every mode. The other host must send CORS headers; the build adds it to the content-security policy. |
| `attribution` | `true` | Show the providers' attribution lines next to the data. Forced to `true` in commercial mode, where the notices are licence conditions; the registry tables are always shown. |
| `overpassUrl` | `""` | Address of an Overpass API server the operator runs or rents. Needed for the optional OpenStreetMap detail in commercial mode, because the public servers ask commercial users to use their own. |
| `openAlexApiKey` | `""` | Optional OpenAlex API key for the literature lists (the keyless allowance is small). Visible to users. |
| `accept` | `[]` | Registry ids the operator has cleared or licensed itself, for example `["nasa-power"]` after its own legal review, or `["eex-auction"]` with an EEX data licence. Each listed source is then used in commercial mode as in non-commercial mode. |

### What changes in commercial mode

| Feed | Non-commercial | Commercial, no key | Commercial, with key or address |
|---|---|---|---|
| Surface weather | Open-Meteo (MET Norway as fall-back) | MET Norway Locationforecast, then US NWS (United States only) | Open-Meteo customer API, same fall-backs |
| Winds and temperatures aloft, freezing level, visibility | Open-Meteo | not supplied: still air and ISA aloft, freezing level from the standard lapse rate, shown in a note (NWS reports visibility) | Open-Meteo customer API |
| Design temperatures (12 months) | Open-Meteo ERA5 (NASA POWER as fall-back) | off (the case value stays); NASA POWER if listed under `accept` | Open-Meteo customer API (Professional plan) |
| Air quality, sea state | Open-Meteo | off | Open-Meteo customer API |
| Place search | Open-Meteo geocoder + bundled airports | bundled airports (name, town, ICAO/IATA) and coordinates | Open-Meteo customer geocoder + bundled airports |
| Terrain elevation | Open-Meteo | MET Norway model height at the site; airport elevation from the bundled database | Open-Meteo customer API |
| Aerodromes and runways | OurAirports (bundled), optional OpenStreetMap through the public Overpass servers | OurAirports (bundled) | plus OpenStreetMap through `overpassUrl` |
| Exchange rates | ECB Data Portal, Frankfurter as fall-back | same | same |
| Brent | EIA series via the open dataset (browser) and EIA directly (snapshot) | EIA directly, through the snapshot | same |
| Jet fuel, policy rates | EIA; Federal Reserve Board / New York Fed; ECB | same | same |
| Carbon price | EEX EU allowance auction (latest price, attributed) | UK ETS Authority determination (GOV.UK, OGL v3.0), labelled as a UK annual figure | the operator's own price (`carbonPriceUrl`) |
| Space weather, World Bank indicators, GB grid intensity, literature (OpenAlex) | unchanged | unchanged | unchanged |
| Open-source tools list | GitHub repository search | a link that opens the search on github.com | same |

Weather providers return one normalised object (`js/core/live.js`, `WEATHER`): time, elevation,
temperature, humidity, station and sea-level pressure, wind, gust, precipitation, cloud, visibility,
freezing level and winds aloft. MET Norway supplies no station pressure (derived from sea-level pressure
and elevation), visibility, freezing level or winds aloft, and gusts only in the Nordic area; the US NWS
supplies the latest station observation without freezing level or winds aloft. Missing fields are named
in `missing`, and the location and Live data pages show the note. MET Norway requires each client to
identify itself, which a page does through its `Origin` header: the single-file copy opened from disk
therefore skips MET Norway.

### Before selling: operator checklist

1. Set `"commercial": true` in `config.json`, run `npm run build`, deploy, and run `node tests/licences.mjs`.
2. Decide on weather: either accept MET Norway / NWS coverage (identify the site by serving it from its own
   domain; above low traffic put a caching proxy in front of `api.met.no`, as its terms require), or buy an
   Open-Meteo subscription and set `openMeteoApiKey`.
3. Carbon price: accept the UK ETS annual figure, or contract a market-data source and publish it at
   `carbonPriceUrl`.
4. Ship `LICENSE`, `THIRD_PARTY_NOTICES.md` and the licence texts under `js/vendor/` with every copy (the
   build and the publish workflow do), and offer the source of the LGPL-licensed OpenCASCADE kernel as
   `THIRD_PARTY_NOTICES.md` describes.
5. Have counsel review the items marked “legal review advised” or `needsLegalDecision` in the registry.

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
do not demonstrate that a model represents a particular aircraft. Models that cannot be solved credibly in
a browser (three-dimensional RANS/LES/DNS, full-aircraft explicit crash FE, coupled high-fidelity CFD–CSD
and similar) are listed per suite as handed off to external solvers. Default material properties, failure
rates, cost coefficients and correlations are typical or illustrative values to be replaced with sourced
data. Nothing produced here is, on its own, airworthiness evidence.
