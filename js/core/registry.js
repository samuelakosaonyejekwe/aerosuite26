// Static index of the 26 suites. Suite modules are loaded on demand so the shell starts instantly.

export const GROUPS = ['Aerodynamics', 'Structures', 'Dynamics & Control', 'Propulsion & Power', 'Systems', 'Integration & Assurance', 'Economics'];

export const SUITES = [
  { n: 1, id: 'cfd', file: 's01-cfd', title: 'Aerodynamics and Computational Fluid Dynamics (CFD)', short: 'Aerodynamics & CFD', group: 'Aerodynamics', icon: 'flow' },
  { n: 2, id: 'fea', file: 's02-fea', title: 'Structural Mechanics and Finite Element Analysis (FEA)', short: 'Structures & FEA', group: 'Structures', icon: 'truss' },
  { n: 3, id: 'aeroelastic', file: 's03-aeroelastic', title: 'Aeroelasticity and Fluid–Structure Interaction (FSI)', short: 'Aeroelasticity & FSI', group: 'Structures', icon: 'flutter' },
  { n: 4, id: 'flightdyn', file: 's04-flightdyn', title: 'Flight Dynamics, Stability and Control', short: 'Flight Dynamics', group: 'Dynamics & Control', icon: 'axes' },
  { n: 5, id: 'performance', file: 's05-performance', title: 'Aircraft Performance and Flight Envelope', short: 'Performance & Envelope', group: 'Aerodynamics', icon: 'envelope' },
  { n: 6, id: 'rotorcraft', file: 's06-rotorcraft', title: 'Rotorcraft Aerodynamics and Aeromechanics', short: 'Rotorcraft Aeromechanics', group: 'Aerodynamics', icon: 'rotor' },
  { n: 7, id: 'propulsion', file: 's07-propulsion', title: 'Propulsion and Engine Performance', short: 'Propulsion & Engines', group: 'Propulsion & Power', icon: 'engine' },
  { n: 8, id: 'propeller', file: 's08-propeller', title: 'Propeller and Rotor Performance', short: 'Propeller & Rotor', group: 'Propulsion & Power', icon: 'prop' },
  { n: 9, id: 'fatigue', file: 's09-fatigue', title: 'Fatigue, Fracture and Damage Tolerance', short: 'Fatigue & Fracture', group: 'Structures', icon: 'crack' },
  { n: 10, id: 'vibration', file: 's10-vibration', title: 'Vibration, Modal and Rotor Dynamics', short: 'Vibration & Modal', group: 'Structures', icon: 'wave' },
  { n: 11, id: 'acoustics', file: 's11-acoustics', title: 'Aeroacoustics and Noise Prediction', short: 'Aeroacoustics & Noise', group: 'Aerodynamics', icon: 'sound' },
  { n: 12, id: 'thermal', file: 's12-thermal', title: 'Thermal Engineering and Heat Transfer', short: 'Thermal & Heat Transfer', group: 'Systems', icon: 'heat' },
  { n: 13, id: 'icing', file: 's13-icing', title: 'Aircraft Icing and Ice Protection', short: 'Icing & Ice Protection', group: 'Systems', icon: 'ice' },
  { n: 14, id: 'gear', file: 's14-gear', title: 'Landing Gear, Ground Dynamics and Impact', short: 'Landing Gear & Ground', group: 'Dynamics & Control', icon: 'gear' },
  { n: 15, id: 'crash', file: 's15-crash', title: 'Crashworthiness and Impact Mechanics', short: 'Crashworthiness & Impact', group: 'Structures', icon: 'impact' },
  { n: 16, id: 'control', file: 's16-control', title: 'Flight Control Systems and Autopilot', short: 'Flight Control & Autopilot', group: 'Dynamics & Control', icon: 'loop' },
  { n: 17, id: 'avionics', file: 's17-avionics', title: 'Avionics, Navigation and Sensor Systems', short: 'Avionics & Navigation', group: 'Dynamics & Control', icon: 'radar' },
  { n: 18, id: 'hydmech', file: 's18-hydmech', title: 'Hydraulic, Pneumatic and Mechanical Systems', short: 'Hydraulic & Mechanical', group: 'Systems', icon: 'piston' },
  { n: 19, id: 'electrical', file: 's19-electrical', title: 'Electrical and Hybrid-Electric Power Systems', short: 'Electrical & Hybrid', group: 'Propulsion & Power', icon: 'bolt' },
  { n: 20, id: 'fuelecs', file: 's20-fuelecs', title: 'Fuel Systems and Environmental Control', short: 'Fuel & ECS', group: 'Systems', icon: 'drop' },
  { n: 21, id: 'composites', file: 's21-composites', title: 'Materials and Composite Mechanics', short: 'Materials & Composites', group: 'Structures', icon: 'layers' },
  { n: 22, id: 'safety', file: 's22-safety', title: 'Aircraft Safety, Reliability and Failure Analysis', short: 'Safety & Reliability', group: 'Integration & Assurance', icon: 'shield' },
  { n: 23, id: 'mdao', file: 's23-mdao', title: 'Multidisciplinary Design Analysis and Optimisation (MDAO)', short: 'MDAO', group: 'Integration & Assurance', icon: 'pareto' },
  { n: 24, id: 'mission', file: 's24-mission', title: 'Aircraft Mission and Operational Simulation', short: 'Mission & Operations', group: 'Integration & Assurance', icon: 'route' },
  { n: 25, id: 'vvuq', file: 's25-vvuq', title: 'Numerical Verification, Experimental Validation and Uncertainty Quantification', short: 'V&V and UQ', group: 'Integration & Assurance', icon: 'check' },
  { n: 26, id: 'economics', file: 's26-economics', title: 'Aircraft Economics and Techno-Economic Analysis', short: 'Economics', group: 'Economics', icon: 'coin' },
];

export const suiteMeta = (id) => SUITES.find((s) => s.id === id);

// Static import map (rather than a computed path) so bundlers can follow it for the standalone build.
const LOADERS = {
  cfd: () => import('../suites/s01-cfd.js'), fea: () => import('../suites/s02-fea.js'), aeroelastic: () => import('../suites/s03-aeroelastic.js'),
  flightdyn: () => import('../suites/s04-flightdyn.js'), performance: () => import('../suites/s05-performance.js'), rotorcraft: () => import('../suites/s06-rotorcraft.js'),
  propulsion: () => import('../suites/s07-propulsion.js'), propeller: () => import('../suites/s08-propeller.js'), fatigue: () => import('../suites/s09-fatigue.js'),
  vibration: () => import('../suites/s10-vibration.js'), acoustics: () => import('../suites/s11-acoustics.js'), thermal: () => import('../suites/s12-thermal.js'),
  icing: () => import('../suites/s13-icing.js'), gear: () => import('../suites/s14-gear.js'), crash: () => import('../suites/s15-crash.js'),
  control: () => import('../suites/s16-control.js'), avionics: () => import('../suites/s17-avionics.js'), hydmech: () => import('../suites/s18-hydmech.js'),
  electrical: () => import('../suites/s19-electrical.js'), fuelecs: () => import('../suites/s20-fuelecs.js'), composites: () => import('../suites/s21-composites.js'),
  safety: () => import('../suites/s22-safety.js'), mdao: () => import('../suites/s23-mdao.js'), mission: () => import('../suites/s24-mission.js'),
  vvuq: () => import('../suites/s25-vvuq.js'), economics: () => import('../suites/s26-economics.js'),
};
const cache = new Map();
export async function loadSuite(id) {
  if (!cache.has(id)) cache.set(id, LOADERS[id]().then((m) => m.default));
  return cache.get(id);
}

/** Execution order for the integrated run: upstream suites first (Kahn's algorithm, ties by suite number). */
export function executionOrder(defs) {
  const ids = defs.map((d) => d.id), indeg = Object.fromEntries(ids.map((i) => [i, 0])), out = Object.fromEntries(ids.map((i) => [i, []]));
  for (const d of defs) for (const c of d.consumes || []) if (c.from in indeg && c.from !== d.id) { out[c.from].push(d.id); indeg[d.id]++; }
  const num = Object.fromEntries(defs.map((d) => [d.id, d.n])), order = [], ready = ids.filter((i) => !indeg[i]);
  while (order.length < ids.length) {
    if (!ready.length) ready.push(ids.filter((i) => !order.includes(i)).sort((a, b) => indeg[a] - indeg[b] || num[a] - num[b])[0]); // break a cycle
    ready.sort((a, b) => num[a] - num[b]);
    const i = ready.shift();
    if (order.includes(i)) continue;
    order.push(i);
    for (const j of out[i]) if (--indeg[j] === 0) ready.push(j);
  }
  return order;
}
