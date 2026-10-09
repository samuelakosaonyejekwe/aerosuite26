# Suite module contract

Every one of the 26 suites is one ES module in `js/suites/sNN-<id>.js` with a single default export.
Suite modules are **pure computation**: no DOM, no `window`, no `fetch`, no timers, no file access,
no third-party imports. They run unchanged in the browser, in a Web Worker and in Node (the test harness).

Allowed imports only:

```js
import * as N from '../core/numerics.js';      // linear algebra, eig, ODE, roots, optimisation, stats, FFT, gci, check, kv
import { isa, G0, ... } from '../core/atmosphere.js';
import { METALS, PLIES, FLUIDS, BATTERIES } from '../data/materials.js';
import { createSolver } from '../core/solvers/cfd3d.js';   // shared heavy kernels: cfd3d, explicitfe, nsib2d
import { ENGINES } from '../data/ref/engine-emissions.js'; // sourced reference datasets (compact JS wrappers)
```

`js/suites/s05-performance.js` is the reference implementation. Read it first.

## Suite object

```js
export default {
  id: 'performance', n: 5,                 // must match js/core/registry.js
  tagline: 'One plain-language sentence saying what the suite answers.',
  analyses: [ ...analysis objects... ],    // 4–7 per suite
  consumes: [{ from: 'cfd', keys: ['CD0'], why: 'Drag polar' }],   // upstream suites actually read in defaults()
  provides: [{ key: 'range_km', label: 'Range', unit: 'km' }],     // headline outputs other suites may read
  handoff:  [{ model: 'LES / DNS', why: 'why it is not solved natively', tool: 'recommended external solver class' }],
};
```

`handoff` is the honesty list: models named in the specification that are **not** solved natively in the
browser (e.g. 3-D RANS/LES/DNS, explicit crash FE, full CFD-CSD). Never pretend; list them there.

## Analysis object

```js
{
  id: 'field',                              // unique inside the suite
  title: 'Take-off and landing field performance',
  summary: 'One or two sentences a non-specialist can follow.',
  fidelity: 'analytical' | 'reduced-order' | 'numerical',
  equations: ['Takeoff ground-roll equations', ...],   // names of governing equations/models from the spec that this analysis really implements
  applicable: (c, d) => true | 'reason string',         // optional; c = case, d = derived(c)
  inputs: [ { key, label, unit, default, min, max, step, group, help,
              type: 'number' (default) | 'select' | 'bool' | 'text', options: [...], discrete: true } ],
  defaults: (c, up, d) => ({ key: value, ... }),        // pull values from the shared case and upstream suites
  run: (inp, ctx) => result,                             // sync or async; MUST NOT mutate inp/ctx
  convergence: { param: 'nSteps', label: 'Time steps', levels: [25, 50, 100, 200], metric: 'kpi_key', hOf: (n) => 1 / n },  // when the analysis has a mesh / step / sample / panel count
  verify: () => [ N.check(name, actual, expected, relTol, 'reference') ],   // code verification against analytical / exact solutions
  validation: [ { name, source, inputs: {...}, sweep: { key, values: [...] }, target: 'kpi_key', observed: [...], tol_pct: 10 } ],
  calibration: { params: [{ key, min, max }], sweep: 'input_key', target: 'kpi_key', note: 'what measured data to supply' },
  recommend: (res, inp, ctx) => [ { severity: 'info'|'advise'|'warn'|'critical', title, detail, action, basis } ],
}
```

Rules for inputs:

- Every input has a `label`, a `default` and (if numeric) sensible `min`/`max`. `unit` uses SI text (`'m/s'`, `'Pa'`, `'-'`).
- `group` clusters fields in the form (e.g. 'Geometry', 'Flow', 'Material', 'Numerics', 'Limits').
- `help` is a short plain-language hint, with typical values where useful.
- `defaults(c, up, d)` returns only keys it can fill. Return `undefined` for a key to keep the declared default
  (never return `0`/`NaN` for a missing quantity). `up.<suiteId>?.<key>` holds upstream outputs and may be absent.
- Mesh / resolution inputs go in group `'Numerics'` and are marked `discrete: true` when integer.
- The declared defaults alone must give a sensible run (they are used when no case is loaded).

`ctx` = `{ case, d, up, atm, progress(fraction, message) }` where `d = derived(case)` (see `js/core/case.js`)
and `atm = isa(case.atm.alt_m, case.atm.dISA_K)`. Call `ctx.progress?.()` in long loops.

## Result object

```js
return {
  kpis:   [{ key: 'tofl_m', label: 'Take-off distance', value: 1749, unit: 'm', status: 'ok'|'warn'|'bad', note: '' }],
  plots:  [ ...plot specs... ],
  tables: [{ title, columns: ['A', 'B'], rows: [[1, 2], ...] }],
  outputs: { extra_key: number | array },   // optional extra published data; KPI values are published automatically under their key
  warnings: ['validity-limit or failure notes in plain language'],
  models: ['models actually used in this run'],
  assumptions: ['key assumptions'],
};
```

- Results must be plain data (structured-cloneable): numbers, strings, arrays, plain objects. No functions, no typed arrays, no class instances.
- KPI `value` is always a number (NaN/Infinity allowed only when genuinely undefined, and then add a warning).
- KPI keys are `snake_case` with a unit suffix where helpful and are stable: they are the coupling interface.
- `status` expresses engineering acceptability against a stated criterion (put the criterion in `note`).
- A default run must finish in well under 1.5 s in Node; convergence levels must keep the finest level under ~3 s.
- Give at least one plot per analysis, usually 2–4. Every chart needs a title and axis labels with units.

### Plot specs

```js
{ type: 'line',  title, xlabel, ylabel, xlog, ylog, equalAspect,
  series: [{ name, x: [], y: [], style: 'line'|'dash'|'points'|'line+points'|'step' }],
  annotations: [{ x: 12, label: 'Stall' }, { y: 0, label: 'Limit' }] }
{ type: 'bar',   title, ylabel, categories: [], stacked: bool, series: [{ name, y: [] }] }
{ type: 'heat',  title, xlabel, ylabel, zlabel, x: [nx], y: [ny], z: [ny][nx], contours: 10, equalAspect,
  diverging: bool, overlay: [{ name, x: [], y: [] }] }                       // filled contour map on a rectilinear grid
{ type: 'tri',   title, xlabel, ylabel, zlabel, nodes: [[x, y], ...], tris: [[i, j, k], ...], values: [per node],
  equalAspect, edges: bool, diverging: bool, overlay: [...] }               // contour on an unstructured/deformed mesh
{ type: 'polar', title, rlabel, series: [{ name, theta_deg: [], r: [] }] }   // directivity etc.
```

Use NaN in a line series to break the line. Keep series to ≤ 6 per chart and ≤ ~400 points per series
(down-sample long histories). Never put two different units on one chart; make two charts.

## Verification, validation, calibration — keep them honest

- `verify()` compares the implementation with **exact/analytical** results (closed forms, manufactured solutions,
  conservation checks, limiting cases). Every analysis with non-trivial numerics should have at least one check.
  Call `run()` directly and read KPIs with `N.kv(result)`. Tolerances must be justified by the method's accuracy.
- `validation` entries must use **only** reference values you are certain of: textbook-exact results or universally
  documented benchmark values (e.g. thin-airfoil lift slope 2π, Sod shock-tube exact solution, Euler buckling,
  Blasius, Theodorsen). State the source. **Do not invent or approximate experimental data from memory.** If you are
  not certain of the numbers, leave `validation` out — users upload their own test data in the app.
- `calibration` declares which inputs are legitimately uncertain model parameters and which swept input/target pair
  measured data would constrain. The fitting itself is done by the core.

## Recommendations

`recommend()` turns results into decision support: what the numbers mean, what to do next, and why. Each item
names its `basis` (the criterion, rule or regulation paragraph). Include sustainability levers where relevant
(fuel, CO₂, NOx, noise, energy, material use, life extension). Do not cite regulation numbers unless sure of them.

## Published output keys (coupling interface)

These keys must be published (as KPI keys or in `outputs`) by the owning suite, by at least one analysis that is
applicable to the vehicle. Downstream suites read them as `up.<suite>?.<key>` in `defaults()` and must fall back
to the case when absent.

| Suite | Keys |
|---|---|
| cfd | `CL`, `CD`, `CD0`, `k_induced`, `CLa_per_rad`, `CLmax`, `LD_max`, `Cm_ac`, `e_oswald`, `cp_min`, `x_cp_frac` |
| fea | `sigma_max_Pa`, `tip_deflection_m`, `margin_of_safety`, `buckling_factor`, `EI_root_Nm2`, `GJ_root_Nm2`, `wing_struct_mass_kg` |
| aeroelastic | `V_flutter_ms`, `f_flutter_Hz`, `V_divergence_ms`, `V_reversal_ms`, `flutter_margin` |
| flightdyn | `static_margin`, `sp_omega_rads`, `sp_zeta`, `ph_period_s`, `ph_zeta`, `dr_omega_rads`, `dr_zeta`, `roll_tau_s`, `spiral_T2_s`, `trim_alpha_deg`, `trim_de_deg` |
| performance | `V_stall_ms`, `tofl_m`, `ldg_dist_m`, `roc_max_ms`, `ceiling_m`, `V_max_ms`, `range_km`, `endurance_h`, `LD_max`, `fuel_flow_cruise_kgs`, `n_limit`, `V_d_eas`, `brake_energy_J` |
| rotorcraft | `hover_power_W`, `FM`, `CT`, `v_induced_ms`, `coning_deg`, `a1_deg`, `b1_deg`, `collective_deg`, `P_cruise_W`, `hub_vib_N`, `autorotation_index` |
| propulsion | `thrust_N`, `thrust_static_N`, `tsfc_kg_Ns`, `P_shaft_W`, `eta_thermal`, `eta_overall`, `fuel_flow_kgs`, `T4_K`, `EINOx_g_kg`, `co2_kg_s`, `heat_rejection_W` |
| propeller | `prop_thrust_N`, `prop_power_W`, `eta_prop`, `CT_prop`, `CP_prop`, `J_adv`, `FM_prop`, `tip_mach` |
| fatigue | `life_cycles`, `life_fh`, `crit_crack_m`, `inspection_interval_fh`, `damage_per_flight` |
| vibration | `f1_Hz`, `f2_Hz`, `f3_Hz`, `f_torsion_Hz`, `crit_speed_rpm`, `gr_margin` |
| acoustics | `OASPL_dB`, `SPL_peak_dBA`, `bpf_Hz`, `cabin_SPL_dBA`, `footprint_km2` |
| thermal | `T_max_K`, `q_max_Wm2`, `T_recovery_K`, `hx_effectiveness`, `thermal_margin_K` |
| icing | `ice_thickness_mm`, `ice_mass_kg_m`, `dCLmax_pct`, `dCD_pct`, `antiice_power_W`, `freezing_fraction` |
| gear | `gear_load_N`, `gear_load_factor`, `stroke_used_m`, `stop_dist_m`, `brake_temp_K` |
| crash | `peak_g`, `crush_m`, `absorbed_J`, `bird_force_N`, `HIC` |
| control | `gm_dB`, `pm_deg`, `settling_s`, `overshoot_pct`, `bandwidth_rads` |
| avionics | `pos_err_m`, `nav_drift_m_h`, `radar_range_m`, `link_margin_dB` |
| hydmech | `pump_power_W`, `actuator_force_N`, `actuator_rate_ms`, `line_dp_Pa`, `gearbox_eff`, `bearing_L10_h` |
| electrical | `batt_soc_end`, `batt_temp_K`, `bus_V_min`, `elec_losses_W`, `endurance_elec_min`, `motor_eff`, `energy_used_kWh` |
| fuelecs | `fuel_cg_shift_m`, `slosh_freq_Hz`, `fuel_pump_power_W`, `cabin_alt_m`, `ecs_cooling_W`, `bleed_kgs` |
| composites | `Ex_Pa`, `Ey_Pa`, `Gxy_Pa`, `fpf_load_Npm`, `min_RF`, `laminate_t_m` |
| safety | `p_catastrophic_per_fh`, `system_reliability`, `mtbf_h`, `availability`, `dispatch_reliability` |
| mdao | `opt_mtow_kg`, `opt_AR`, `opt_S_m2`, `opt_fuel_kg`, `opt_improvement_pct` |
| mission | `block_fuel_kg`, `block_time_h`, `mission_energy_kWh`, `reserve_margin_kg`, `co2_kg`, `mission_feasible` |
| vvuq | `gci_fine_pct`, `observed_order`, `validation_rmse` |
| economics | `doc_usd_fh`, `cask_usd`, `npv_usd`, `irr`, `payback_yr`, `lcc_usd`, `co2_cost_usd_fh`, `breakeven_load_factor` |

## Vehicle classes

The test harness runs every analysis on seven presets (jet airliner, turboprop, piston single, helicopter,
fixed-wing UAV, electric quadrotor, eVTOL). Wing area is 0 for the helicopter and quadrotor; fuel mass is 0 for
electric aircraft; rotor radius is 0 for pure aeroplanes. Use `applicable` to exclude an analysis from a vehicle
it cannot represent, and make `defaults()` robust (no division by zero, no NaN) for every preset it does accept.

## Testing

```
node tests/run.mjs <suiteId> [<suiteId> ...]     # contract, all presets, verification, convergence
node tests/run.mjs                               # everything plus the coupled chain
```

A suite is done when its run reports `0 failures` and its numbers are physically plausible for each preset.
