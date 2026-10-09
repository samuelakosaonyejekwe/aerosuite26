// Display units. Everything is computed and stored in SI; this layer only changes what people see and
// type. "Aviation" shows the units used on flight decks and in most operational documents.

import { state } from '../core/store.js';

// SI unit string -> [display unit, factor, offset]
const AVIATION = {
  m: ['ft', 3.280839895, 0], km: ['nmi', 0.5399568035, 0], 'm/s': ['kt', 1.943844492, 0], kg: ['lb', 2.204622622, 0], N: ['lbf', 0.2248089431, 0],
  Pa: ['psi', 1.450377377e-4, 0], 'N/m²': ['lbf/ft²', 0.02088543423, 0], 'm²': ['ft²', 10.76391042, 0], 'm³': ['ft³', 35.31466672, 0], W: ['hp', 1 / 745.6998716, 0],
  K: ['°F', 1.8, -459.67], '°C': ['°F', 1.8, 32], 'kg/s': ['lb/h', 7936.641439, 0], 'N·m': ['lbf·ft', 0.7375621493, 0], 'Nm': ['lbf·ft', 0.7375621493, 0],
  mm: ['in', 0.03937007874, 0], 'kg/m³': ['lb/ft³', 0.06242796058, 0], J: ['ft·lbf', 0.7375621493, 0], 'm/s²': ['ft/s²', 3.280839895, 0], 'kg m²': ['slug·ft²', 0.7375621493, 0],
  'USD/kg': ['USD/lb', 1 / 2.204622622, 0], 'km²': ['nmi²', 0.2915533496, 0],
};
export const unitSystem = () => state.settings.units || 'si';
const entry = (unit) => (unitSystem() === 'aviation' ? AVIATION[unit] : null);
/** Convert an SI value for display. Returns [value, unitLabel]. */
export function toDisplay(v, unit) { const e = entry(unit); return e && typeof v === 'number' ? [v * e[1] + e[2], e[0]] : [v, unit]; }
/** Convert a typed display value back to SI. */
export function fromDisplay(v, unit) { const e = entry(unit); return e ? (v - e[2]) / e[1] : v; }
export const displayUnit = (unit) => entry(unit)?.[0] ?? unit;
