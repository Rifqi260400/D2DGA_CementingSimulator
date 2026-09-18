# Solver inventory for the UI

Read before writing any UI code, per the build spec §2. Authority order: **solver
code > mockup > spec**. Every row below was read from the code at `9f6b9b2`, not
from the mockup.

Mockup read from `https://claude.ai/artifact/2CRc9Wdrvg4EPMkT6jNyJB`, four
artboards: `Main` (Case setup), `Run`, `Results`, `Validation`.

---

## A. Inputs

`d2dga/config.py` is already a frozen-dataclass tree (`Config` → `WellConfig`,
`StandoffConfig`, `SyntheticWallConfig`, `FluidsConfig`, `GridConfig`) and is
the CLI's only source of defaults. It is a usable single source of truth, with
one gap (A-13).

### A.1 Present in the solver

| # | Parameter | Type / unit | Default | Range | Where | Mockup control |
|---|---|---|---|---|---|---|
| A-1 | `well.casing_od_m` | m | `0.17780` (7 in) | > 0, < hole | `WellConfig` | "Casing OD (mm)" — **mockup shows 127.0, stale** |
| A-2 | `well.gauge_hole_diameter_m` | m | `0.26416` (10.4 in) | > casing OD | `WellConfig` | "Gauge hole (in)" ✓ |
| A-3 | `well.casing_shoe_m` | m | `194.0` | 0 … TD | `WellConfig` | "Open hole (m)" (derived: TD − shoe = 196) |
| A-4 | `well.total_depth_m` | m | `390.0` | > shoe | `WellConfig` | — (not shown) |
| A-5 | `well.washout_top_m` / `_bottom_m` / `washout_max_diameter_m` | m | 195 / 217 / 0.5842 | — | `WellConfig` | — (only via wall preset) |
| A-6 | `well.inclination_rad` | rad | `0.0` | **0 only** on the tabulated path | `WellConfig` | — |
| A-7 | `standoff.eccentricity_in_gauge_hole` | – | `0.3` | 0 … <1 | `StandoffConfig` | "ê (mm)" — mockup uses offset in mm, solver uses gauge eccentricity; convertible |
| A-8 | `synthetic_wall.amplitude_m` | m | `0.0381` (1.5 in) | feasibility-checked | `SyntheticWallConfig` | "Amplitude A (in)" ✓ |
| A-9 | `synthetic_wall.wavelength_m` | m | `20.0` | > 0 | `SyntheticWallConfig` | "Wavelength L (m)" ✓ |
| A-10 | `synthetic_wall.mode` | `symmetric` \| `enlargement` | `symmetric` | — | `SyntheticWallConfig` | — (not shown; matters, see GEO-04) |
| A-11 | `synthetic_wall.amplitude_is_diameter` | bool | `True` | — | `SyntheticWallConfig` | — (not shown; GEO-03) |
| A-12 | `fluids.mud_density`, `mud_viscosity` | kg/m³, Pa·s | 998, 1e-3 | > 0 | `FluidsConfig` | mud ρ̂ ✓ |
| A-13 | `fluids.cement_density`, `_consistency`, `_power_law_index`, `_yield_stress` | kg/m³, Pa·sⁿ, –, Pa | 1200, 0.6, 0.4, 1.4 | > 0; 0 < n ≤ 1 | `FluidsConfig` | cement row ✓ |
| A-14 | `fluids.mean_velocities_m_s` | m/s | `(0.05, 0.2, 0.5)` | > 0 | `FluidsConfig` | "Flow rate Q̂" — solver takes **ŵ₀**, not Q̂ (both accepted by `Scaling`) |
| A-15 | `grid.n_phi`, `grid.n_xi` | cells | 20, 400 | ≥ 2 | `GridConfig` | "Nφ", "Nξ" ✓ |
| A-16 | `Simulation(cfl=)` | – | `0.5` | `0 < cfl ≤ 1` (enforced) | `TransportSolver.__init__` | "CFL" ✓ |
| A-17 | `Simulation(inflow=)` | `no_axial_gradient` \| `uniform` | `no_axial_gradient` | — | `Simulation` | "Inflow BC" ✓ |
| A-18 | `Simulation(inflow_concentration=)` | – | `1.0` | `[0, 1]` (enforced) | `Simulation` | — |
| A-19 | `Simulation(picard_tol=, picard_max_iter=, picard_warm_start=)` | –, int, bool | 1e-8, 100, True | > 0 | `Simulation` | — → **Advanced** |
| A-20 | `Simulation(mobility_floor=)` | – | `1e-12` | > 0 | `Simulation` | — → **Advanced** (NUM-13) |
| A-21 | `TransportSolver(wavespeed=)` | `interval` \| `endpoints` | `interval` | — | `TransportSolver` | — → **Advanced, with a warning**: `endpoints` is NUM-26's non-monotone form, kept only to demonstrate failure |
| A-22 | `TransportSolver(inlet_flux=)` | `upwind` \| `llf` | `upwind` | — | `TransportSolver` | — → same treatment (NUM-29) |
| A-23 | `run(t_end=)` via `--volumes` | annulus volumes | `1.2` | > 0 | `kgep1_run.py` | "Pumped volume (m³)" — convert with `Geometry.annulus_volume_direct()` |
| A-24 | wall choice | `synthetic` \| `caliper` (CLI); 4 classes exist | `synthetic` | — | `geometry.py` | "Wall profile" — solver has **Uniform, Sinusoidal, Washout, CaliperLog** ✓ all four |
| A-25 | closure path | implicit | chosen by fluid type | — | `NewtonianClosures` / `TabulatedClosures` | "Closures" — **solver picks this from the fluids, it is not a free choice**; a Newtonian pair cannot use the table meaningfully and an HB pair cannot use the analytic form |
| A-26 | table axes `--n-c`, `--n-h` | int | 31, 9 | ≥ 2 | `kgep1_run.py` | — → **Advanced** (NUM-03 is still open) |

### A.2 Mockup controls with **no** solver parameter

| Mockup field | Verdict |
|---|---|
| **Mud τ̂_Y, κ̂, n** (shown as 4.79 Pa, 0.020, 0.70) | **`FluidsConfig` has no such fields.** `as_fluids()` hard-codes the mud as Newtonian: `HerschelBulkleyFluid(..., mud_viscosity, 1.0, 0.0)`. See gap **G-1**. |
| **Offset model "From centralizers"** | `EccentricityModel` is an ABC with two concrete subclasses, both constant-offset. No centralizer model exists. `BLOCKED.md` BLK-1. Show the option **disabled with the BLK-1 reason**. |
| **Reynolds number Re** | `Scaling.reynolds_zf22` exists, but the D2DGA model is **non-inertial — Re never enters the solution**. ZF22's own cases 9/10 differ only in Re and give identical answers. Display it labelled "reported, not used by the model", or drop it. |
| **"Predicted regime, ZF23 map (b > 80)"** | **No such predictor exists, and ZF23 has no such map.** ZF23 (3.6) is a *post-hoc* classification computed from a finished field. See gap **G-5**. |
| **"Pause"** button (Run monitor) | `Simulation.run` has no pause. It has `on_step` and periodic checkpointing. Stop = terminate the process; the checkpoint makes it resumable. See **G-6**. |
| **"Snapshot" button / snapshot slider "6 of 9"** | **No snapshots are saved anywhere.** See **G-2** — the largest gap. |
| **"Compare with / Overlay"** (Results) | No mechanism; needs two runs on disk, which needs G-2 first. |
| **η_N "narrow quartile"** | `narrow_side_profile` returns the narrow-side *column*, not a quartile-weighted scalar. No `eta_N` exists. See **G-4**. |
| **Δw_f** (Results, one of the four classification numbers) | `dispersion_metrics` returns `sigma_plus, sigma_minus, area_plus, area_minus`. **`Δw_f` is not computed.** The four numbers the solver has are σ_w+, σ_w−, \|w̄+\|, \|w̄−\|. |
| **"Flux through slice Q ± 3e−13"** (live invariant) | `StreamFunctionSolver.axial_flux(psi)` exists but is **not called during a run** — only in tests. Cheap to add via `on_step`. |
| **"10 / 10 WITHIN TOLERANCE"** (Validation) | **Reality is 9/10.** ZF22 case 1 does not reproduce (`BLOCKED.md` BLK-6). Must not be shown as 10/10. |
| **Blocked: "Casing outer diameter", "Caliper log"** | Both stale: GEO-01 is **resolved** (7 in, user-selected) and the caliper log **exists** (`data/K-GEP-01_2024-03-30.las`). The real blocked list is `BLOCKED.md`'s 8 items, led by the **centralizer record** (130× sensitivity). |
| **"Stage 2 · Newtonian, sinusoidal wall"** | No stage concept in the solver. Drop. |

---

## B. Outputs

### B.1 In-process, returned by `Simulation.run()` → `RunResult`

| Field | Shape / dtype | When |
|---|---|---|
| `concentration` | `(n_phi, n_xi)` float64 | end |
| `stream_function` | `(n_phi+1, n_xi+1)` float64 | end |
| `times`, `efficiency`, `outlet_concentration` | `(n_steps+1,)` float64 | **every step** (since R-1) |
| `reports` | `list[StepReport]`, one per `record_every` | incremental |
| `breakthrough_time`, `conservation_error` | float | end |
| `efficiency_at(t)`, `breakthrough_at(threshold)` | float | end |

`StepReport` fields: `n, t, dt, Q, mass, efficiency, boundary_flux, c_min, c_max, picard_iterations, static_cells`.

### B.2 Post-processing, `d2dga/postprocess.py` — **call these, never reimplement**

| Function | Returns |
|---|---|
| `cell_volume(geometry)` | `(n_phi, n_xi)` volume weights |
| `displacement_efficiency(geometry, c)` | scalar η_E |
| `axial_profile(geometry, c)` | `(n_xi,)` volume-weighted c̄(ξ) |
| `narrow_side_profile(geometry, c)` | `(n_xi,)` — the narrow column `c[-1, :]` |
| `residual_fraction(geometry, c, threshold=0.5)` | scalar |
| `front_speed_profile(geometry, c, t, n_bins=100)` | `(levels, w_f)` |
| `dispersion_metrics(levels, w_f)` → `DispersionMetrics` | `sigma_plus, sigma_minus, area_plus, area_minus, n_bins`, `.is_dispersive` (ZF23 3.6: σ_w+ > 0.08 **and** \|w̄+\| > 0.05) |
| `zf23_metrics(geometry, c, t, n_bins)` | the above in one call |

### B.3 Geometry, for the validity envelope

`Geometry.narrow_gap_parameter(xi)`, `.e(xi)`, `.r_a(xi)`, `.H(phi, xi)`,
`.depth(xi)`, `.wall_gradient(xi)`, `.annulus_volume_direct()`, `.summary()`.
`CaliperLogWall.stats` carries the log QC dict.

### B.4 Derived dimensionless groups, `Scaling` — for the live setup panel

`m`, `B`, `Fr_star`, `delta_rho`, `buoyancy_number` (= b), `w0_hat`,
`flow_rate_hat`, `gamma_dot_0`, `mu_e_hat`, `tau_0_hat`, `time_hat`,
`length_hat`, `reynolds_zf22`, `.summary()`.

### B.5 Closure-table diagnostics (added during remediation)

`ClosureTable.range_report()`, `.assumption_report()`, `.angle_sensitivity`,
`.n_range_points_out`, `.worst_excursion`, `.tuned_r`, `.n_failed`.

### B.6 On disk after a completed run

**Updated 2026-09-18 (G-2 / R-2).** A run now leaves three files with a common
stem, plus the closure table:

| File | Written by | Contents |
|---|---|---|
| `output/kgep1_ckpt_<...>.npz` | `Simulation.run` | `c` (final), `t`, `n`, `mass_running`, `cons_err`, `t_br`, `times`, `effs`, `outlet`, picard counters, **`tag`** (settings fingerprint) |
| `output/kgep1_ckpt_<...>.config.json` | `runio.write_config` | whole `Config` tree, CLI args, `physical` subset, command line, git SHA + dirty flag, Python/numpy/scipy versions |
| `output/kgep1_ckpt_<...>.metrics.json` | `runio.write_metrics` | `eta_E`, `t_br` at three thresholds **with their half-order error bars**, narrow-side minimum, residual fraction, volume balance, ZF23 metrics (all four + `is_dispersive` + sample time), static cells, picard counters, outlet import, closure-table range report, angle-sensitivity report, the `Scaling` groups, and the geometry ranges |
| `output/closure_table_<hash>.npz` | `kgep1_run.make_table` | the table, keyed by fluid pair + axes |

Still **not** persisted: `ψ`, and intermediate snapshots of `c̄`. So the Results
screen can show the final field and the full scalar history, and cannot show a
snapshot slider or a similarity collapse.

<details><summary>What it was before (kept, because it is why R-2 existed)</summary>

### The old situation

| File | Written by | Contents |
|---|---|---|
| `output/kgep1_ckpt_<...>.npz` | `Simulation.run` | `c` (final only), `t`, `n`, `mass_running`, `cons_err`, `t_br`, `times`, `effs`, `outlet`, `n_picard_unconverged`, `worst_picard_residual` |
| `output/closure_table_<hash>.npz` | `kgep1_run.make_table` | the table, keyed by fluid pair + axes |
| `output/*.log` | shell redirect of stdout | free text |

**Nothing else.** `scripts/kgep1_run.py` wrote **no machine-readable result
file** — it printed. `output/kgep1_results.md` was written by hand, not by the
script. There was **no saved config, no ψ, no snapshots, no metrics** — and no
way to tell two runs apart, which is R-2.

</details>

---

## C. Live signals

| Signal | Available? | How, today |
|---|---|---|
| step `n`, `t`, `dt`, `Q` | ✅ | `StepReport` → `on_step(rep, c, psi)` callback |
| `c_min`, `c_max` | ✅ | `StepReport` |
| volume error | ✅ | running `cons_err` in `run`; `StepReport.mass` + `boundary_flux` |
| efficiency, outlet c̄ | ✅ | `StepReport`, and the full arrays every step |
| Picard iterations | ✅ | `StepReport.picard_iterations` |
| static cells (NUM-13) | ✅ | `StepReport.static_cells` |
| closure-table range | ✅ | `closures.range_report()` — **live state**, readable at any moment |
| outlet import (A-1) | ✅ | `transport.outlet_import` — live accumulator |
| angle sensitivity (B-1) | ✅ | `table.assumption_report()` — fixed after build |
| **c̄ field, ψ field** | ✅ | passed to `on_step` |
| axial flux through a slice | ⚠️ | `elliptic.axial_flux(psi)` exists, **not called during a run** |
| front position wide/narrow | ⚠️ | derivable from `c` in `on_step` via `postprocess` |
| **anything outside the process** | ❌ | only stdout and the 180 s checkpoint |

**Consequence:** every live signal the mockup wants already exists *inside* the
process. Nothing needs to be added to the numerics. What is missing is a way to
get them *out*.

---

## D. Gaps and proposed resolutions

Ordered by how much they block the build. **G-1, G-2 and G-3 need your decision
before I write UI code** (spec §9).

| ID | Gap | Proposal |
|---|---|---|
| **G-1** | **Not approved — left as is.** `FluidsConfig` cannot express a non-Newtonian mud. The mockup asks for mud τ̂_Y, κ̂, n; `as_fluids()` hard-codes `(mud_viscosity, n=1, τ_Y=0)`. | **Add `mud_consistency`, `mud_power_law_index`, `mud_yield_stress` to `FluidsConfig`, defaulting to the current Newtonian water** (1e-3, 1.0, 0.0) so every existing number is bit-identical. ~8 lines, one dataclass + `as_fluids`. **Needs your OK — it touches `config.py`.** Without it the UI cannot offer the fluid the well actually has. ⚠️ Note this walks straight into finding **B-1**: with a yield stress in both fluids the closure table's θ = 0 slice is O(1) wrong. The UI must surface `assumption_report()` prominently on that path. |
| **G-2** | ~~A completed run persists nothing but a checkpoint.~~ **DONE 2026-09-18, and it turned out to hide a correctness hazard, not just a bookkeeping one.** | `d2dga/runio.py` now writes `<checkpoint>.config.json` (whole `Config` tree + CLI args + git SHA + library versions) and `<checkpoint>.metrics.json` (every number the run prints). The checkpoint path is unchanged, so runs recorded earlier still resume. **The hazard:** `Simulation.run` resumed on a matching path and grid *shape* only, and fluid properties are in neither — so editing `mud_density` and relaunching the same command silently continued a field computed under the old physics. `Simulation.run(checkpoint_tag=)` now refuses across a settings change, and `runio.check_resume` names the setting that moved. Logged as **R-2** in `docs/remediation_log.md`, registered as REM-10. 9 tests; verified end to end on the real runner with `--cfl` changed. |
| **G-3** | Live signals cannot leave the process. | **No solver change needed.** `on_step` is already a public hook. The UI's launcher passes an `on_step` that atomically rewrites `status.json` every ~2 s (temp file + `os.replace`, the pattern `run`'s checkpointing already uses); Streamlit polls it. Runs stay plain subprocesses, so **a UI run and a CLI run are the same process invocation** — which is what makes §8's byte-identical requirement testable. I'll state this choice in the UI README. |
| **G-4** | No `eta_N` (narrow-side efficiency) scalar. | Add `narrow_side_efficiency(geometry, c, fraction=0.25)` to `postprocess.py` — volume-weighted η over the narrowest quartile of φ. It belongs in the solver, not the UI (spec §5). ~6 lines. Low risk; I'd fold it into G-2's OK. |
| **G-5** | Mockup's "predicted regime from the ZF23 map (b > 80)" does not exist and ZF23 has no such map. | **Replace with the real a-priori predictor I added during remediation:** BF25 §3.3's Muskat criterion, `tests/benchmarks/bf25_muskat.py`. It takes exactly the setup-screen inputs (𝓘₁, 𝓘₂, q₀, 𝓘₃, b) and returns Δw(0⁺) — finger penetrates or is absorbed — validated against BF25's own `M₃^min = 3/2`. That is a genuine prediction from the config, which is what the panel is for. I'd move `bf25_muskat.py` from `tests/benchmarks/` into `d2dga/` so the UI imports it without importing tests. |
| **G-6** | No pause. | Drop the Pause button. Keep **Stop** (terminate; the checkpoint makes it resumable) and say so in the UI. Do not fake a pause. |
| **G-7** | Validation screen has no machine-readable source. `docs/gate_status.md` is prose; the suite prints. | The remediation already used a pytest plugin that emits `{exit, counts, failures}` JSON. **Commit it as `tests/run_gates.py`** and have the Validation screen read its output plus `output/zf22_table3.md` and `BLOCKED.md`. Smallest possible change; no test is modified. |
| **G-8** | Mockup's Validation numbers are wrong (10/10; blocked = casing OD + caliper). | Show the truth: **9/10**, case 1 failing, with a link to BLK-6; and the real `BLOCKED.md` list. |
| **G-9** | `Config` has no field for CFL, inflow BC, Picard settings, wall choice, or volumes — those live as CLI args and `Simulation` kwargs. **Still open.** | G-2 solved *provenance* (a run records what produced it) but not *definition*: the UI still has nowhere single to read CFL's default, range and choices from. **Cheapest fix, and smaller than the `RunConfig` dataclass I first proposed: extract `kgep1_run.py`'s `ArgumentParser` into a `build_parser()` function and have the UI introspect it** for names, defaults, `choices` and help text. That is genuinely one source of truth with no restatement, and it is a ~3-line move inside a script rather than a change to `config.py`. **Needs your OK before the Case-setup screen is built** — until then that screen would have to restate defaults, which §5 forbids. |
| **G-10** | Re is reported but unused by the model. | Show it greyed with "not used — D2DGA is non-inertial", rather than dropping it: its absence is informative (ZF22 cases 9/10 test exactly this). |

## E. Proposed deviations from the mockup

| Deviation | Reason |
|---|---|
| Snapshot slider is driven by a **requested snapshot schedule** set at launch, not a fixed "9 of 9" | Snapshots do not exist yet (G-2); their number is a config choice |
| "Closures" becomes a **read-only derived field** | The solver picks the path from the fluid types; it is not a free choice (A-25) |
| Fluid preset "Bittleston 2002 Example 1" replaced by the actual default (Materials 2025 Table 1) with its FLU-01 caveat shown | The config's fluids are not Bittleston's, and the caveat is load-bearing |
| A fifth screen, **Geometry & envelope**, for the A–L sweep and caliper QC | `scripts/phase1_sweep.py` and `CaliperLogWall.stats` produce a comparison matrix that will not fit in Results (spec §4 invites this) |
| Casing OD default 177.8 mm, not 127.0 | GEO-01, user-selected 2026-09-16 |
| Advanced section for A-19…A-22, A-26 | `wavespeed="endpoints"` and `inlet_flux="llf"` exist only to demonstrate known-broken behaviour and must be labelled as such |

## F. What I will **not** show

No panel gets a number the solver did not produce. Specifically blank-with-a-reason
until the underlying gap closes: snapshots and the similarity-collapse chart
(needs G-2), η_N (needs G-4), Δw_f (not computed by `dispersion_metrics` at all —
I will show the four metrics the solver *does* compute and label the axis
accordingly), centralizer offset model (BLK-1), and the caliper comparison
overlay (no completed caliper run exists, before or after remediation).
