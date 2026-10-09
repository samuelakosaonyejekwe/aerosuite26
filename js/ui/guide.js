// "How this page works": a short plain-language guide at the top of every page. It opens by itself the
// first time a page is visited, remembers when it is closed, and can be reopened with the ? button.

import { h, icon } from './dom.js';
import { ls } from '../core/store.js';

const G = {
  home: { what: 'Your starting point: where the study stands and what to do next.', steps: ['Follow the four numbered steps from left to right.', 'The three cards show the aircraft, the operating site and how many suites have results.', 'Every suite is one tap away in the lists below; a coloured dot shows its last outcome.'], tip: 'Press / or Ctrl+K anywhere to jump straight to any suite or analysis.' },
  case: { what: 'The one place where the aircraft and its operating site are described. All 26 suites read from here.', steps: ['Pick the representative aircraft closest to yours, then change any value.', 'Open “Location & live conditions”, search for the airport or place, and choose a runway.', 'Or import a file of values under “Import & export”. Nothing needs to be typed twice.'], tip: 'The consistency checks warn about values that do not fit together before you run anything.' },
  geometry: { what: 'Bring in shapes and meshes made in other tools. This page reads, checks and measures them; it does not draw.', steps: ['Drop a file. The page says what format it really is and how completely it was read.', 'Confirm the length unit if the file does not state one.', 'Inspect quality, heal defects if needed, then cut a section or send dimensions to the case.'], tip: 'Every repair and unit change is written to the provenance record at the bottom.' },
  integrated: { what: 'Runs every suite in the right order so each one receives the results it depends on.', steps: ['Tap a numbered circle to see what that suite receives and what it passes on.', 'Press “Run all suites”. Two passes let later results flow back upstream.', 'Read the outcome table, then open Decision support.'], tip: 'Analyses that do not apply to this aircraft type are skipped and listed, not hidden.' },
  bridge: { what: 'Prepares complete input files for full-scale open-source solvers and reads their results back, for work too large for a browser.', steps: ['Choose the solver and the kind of analysis.', 'Download the ready-to-run case and run it on your computer, a cluster or free cloud runners.', 'Drop the result files back here to compare them with the built-in suites.'], tip: 'Use this when you need body-fitted meshes with millions of cells or certification-grade detail.' },
  decisions: { what: 'Every finding from every suite in one list, most urgent first, with what to do about it.', steps: ['Start with “Act now”, then “Needs attention”.', 'Each card names the suite it came from: tap it to see the numbers behind it.', 'The scorecard and sustainability panel summarise the whole aircraft.'], tip: 'Each recommendation states the rule or criterion it rests on.' },
  live: { what: 'The outside data the app uses, where each comes from and how fresh it is.', steps: ['Green means fresh; amber means a saved value is in use until you are online.', 'Press Refresh on any card, or “Refresh everything”.', 'Switch off “Write live values into the case” to keep your own numbers.'], tip: 'Your browser fetches these directly from the providers; nothing passes through a server of this app.' },
  reports: { what: 'Turn the study into something others can check and reproduce.', steps: ['“Save project” stores the case, every input and all results in one file.', '“Build printable report” assembles all results with their assumptions; print or save as PDF.', '“Verify all suites” re-runs every built-in benchmark against exact solutions on this device.'], tip: '“Share case as a link” copies an address that opens this exact aircraft anywhere.' },
  about: { what: 'How to install the app, use it without a connection, and what it can and cannot do.', steps: ['Follow the install steps shown for this device.', '“Store everything for offline use now” makes aeroplane mode safe immediately.', 'The table lists what each suite computes on the device and what it hands to external solvers.'], tip: 'Installed copies update themselves the next time they are opened online.' },
  share: { what: 'Someone sent you an aircraft case as a link.', steps: ['Check the name shown.', 'Press “Open this case” to load it, or Cancel to keep yours.'], tip: 'The case travelled inside the link itself; nothing was uploaded.' },
};
const TAB = {
  run: ['Set the inputs and solve.', ['Blue fields are filled from your case, live data or an earlier suite; change any of them and they turn amber.', 'Press the large Run button. Results appear on the right: headline numbers, what they mean, charts and tables.', 'Hover or tap a chart to read values; use its icons to see the data table or download it.']],
  mesh: ['Check that the answer does not depend on how finely the problem is divided.', ['The solver is repeated at several resolutions.', 'You get the observed order of accuracy, an error estimate and the coarsest resolution that meets your target.', 'Press “Use this resolution” to adopt it for this analysis.']],
  studies: ['Explore how the answer responds to its inputs.', ['Ranking shows which inputs matter most.', 'A sweep varies one input across a range.', 'Uncertainty propagation gives the likely spread of the result when inputs are not known exactly.']],
  vv: ['Evidence that the result can be trusted.', ['Verification compares the solver with exact mathematical solutions.', 'Validation compares the model with measurements you paste in.', 'Calibration adjusts chosen parameters to fit measurements; keep separate data for validation.']],
  spec: ['What this suite is required to cover and how each analysis does it.', ['Green terms are equations and models an analysis here implements.', 'The hand-off table lists models that need an external solver, and why.', 'Further down: initial and boundary conditions, inputs and outputs, accepted geometry formats and data exchanged with other suites.']],
  live: ['Current outside knowledge for this discipline.', ['Recent or most-cited papers on the left.', 'Maintained open-source solvers and tools on the right.', 'Both are refreshed daily and kept for offline reading.']],
};

/** Build the guide strip for the current route. `parts` is the hash split on '/'. */
export function pageGuide(parts) {
  const name = parts[0] || 'home', isSuite = name === 'suite', tab = isSuite ? parts[3] || 'run' : null;
  const g = isSuite ? { what: TAB[tab]?.[0], steps: TAB[tab]?.[1], tip: 'Choose a different analysis with the rounded buttons above the tabs; greyed-out ones do not apply to this aircraft type.' } : G[name];
  if (!g?.what) return null;
  const key = `guide.${isSuite ? 'suite-' + tab : name}`, open = ls.get(key, true);
  const det = h('details', { class: 'guide', open, ontoggle: () => ls.set(key, det.open) },
    h('summary', null, icon('info', 18), h('b', null, 'How this page works'), h('span', { class: 'muted' }, g.what)),
    h('ol', null, g.steps.map((s) => h('li', null, s))), g.tip ? h('p', { class: 'muted small' }, h('b', null, 'Tip: '), g.tip) : null);
  return det;
}
/** Open or close the guide on the current page (the ? button). */
export function toggleGuide() { const d = document.querySelector('details.guide'); if (d) { d.open = !d.open; if (d.open) d.scrollIntoView({ block: 'nearest', behavior: 'smooth' }); } }
