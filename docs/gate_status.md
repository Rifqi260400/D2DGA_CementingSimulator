# Gate status

One row per M-gate. **A gate that has passed may never silently regress**
(remediation brief §4): if a fix breaks a previously green gate, the fix is
wrong even if it resolves its own finding.

`Last green at` is the commit of the run that produced the state, and every
entry in this file was produced by an **actual run in this session** — none is
carried over from the build.

> **The build spec is not in this repository.** `D2DGA_BUILD_SPEC.md` does not
> exist as a file; it was supplied in the conversation that produced the code.
> The gate list below is therefore reconstructed from the `M<n>-T<k>` markers in
> `tests/` and `docs/assumptions.md`, which is the only surviving enumeration.
> If a gate existed in the spec but was never given a marker, it is not listed
> here and its absence is invisible — recorded as a limitation, not a pass.

**Verdict of the run this table is built from:** `python tests/run_gates.py`
-> **332 passed, 0 failed, 0 errors** (exit 0), 2026-09-23. Before the remediation the suite was
223 passed / 0 failed; everything since has been added, and none removed or
skipped. No previously green gate has regressed.

Since 2026-09-19 the verdict is also written as data, to
`output/gates.json`, by `tests/run_gates.py` — a pytest **plugin**, so it
changes no collection, selection or assertion. The Validation screen reads that
file and **refuses to fall back to this document**, because a page that quotes
prose cannot tell you the prose is stale.

Of those, 8 are `tests/test_runio.py`, added with **R-2**: a checkpoint could be
resumed under different physics, because `Simulation.run` checked only the grid
shape and the checkpoint filename carries six of ~20 settings.

`tests/test_m9_picard_axis.py` (M9-T1..T4) gates **A-3b**, which is the one
finding in this project that a passing suite actively hid. A-3 — my own
remediation — anchored the closure table's velocity axis at zero with a single
node and recorded that it "changes nothing". It put the elliptic Picard
iteration at its cap on **37.5 % of steps**: 6 iterations became 100, residual
1.5e-09 became 3.4e-05. Every gate passed throughout the five days it was
live. It surfaced only when the production case was re-run.

⚠️ **Read that module's docstring before trusting these four gates.** They do
not reproduce the failure. Two cheap tests that would have were tried and both
failed to discriminate — a reduced table gives 52 iterations on both axes, a
synthetic flat front gives 100 on both. What M9 locks is the axis and its
presence in the run fingerprint (which was absent, so R-2's guard could not see
an axis change). The failure itself is caught at **run** level, where the cap
fraction is now printed with a verdict and recorded in `metrics.json`.

`tests/test_m8_conservation.py` (M8-T1..T7) gates **BCF25 §IV**, the
verification this stage of the project rests on: there is no CFD comparison,
so the model is checked against the ten published ZF22 cases and against the
volume ledger that paper applies to the same scheme. Measured **2.04 × 10⁻¹⁵**
against their reported ~10⁻¹⁵, with the total-flux imbalance **exactly 0.0**.

Three of those gates are about the *definition* rather than the magnitude,
because a conservation figure is easy to make look good: M8-T1 fixes the
normalisation (by the annulus volume, BCF25's, not by the volume present);
M8-T3 requires fluid 1's ledger to be built from the closure-free total flux
rather than read off fluid 2's; and M8-T7 injects an unbooked volume and
requires the ledger to report it at the size injected — a conservation check
that cannot fail is decoration. ⚠️ With `K = 2` the two fluid curves are near
mirror images and that is **not** corroboration; see REM-16 and
`docs/assumptions.md`.

`tests/test_ui_contract.py` (UI-1..17) gates the seam between the runner and
the interface, which is where a display can quietly stop describing the thing
it claims to describe. Three of its gates exist because of defects that were
invisible in the source and showed up only in the rendered page or in a
relaunch:

* **UI-12/13/14** — the Validation screen's first ZF22 parser dropped case 1,
  the only case that does not reproduce, and displayed **9 / 9**. It
  manufactured a clean bill of health by discarding its own counter-evidence.
* **UI-6/7** — `RunResult.steps`. Callers read `reports[-1].n`, which is short
  by up to `record_every - 1` on every run and raises `IndexError` on a
  relaunch of a finished one.
* **UI-15/16** — a run with non-default fluids completed, wrote its files, and
  was **invisible** to the interface, because the checkpoint suffix that stops
  such runs colliding was missing from the discovery pattern.
* **UI-17** — the "reproduce from the CLI" line omitted the fluid arguments, so
  a run on a non-default mud was shown a command reproducing a different run.

The 31 before those are `tests/test_m0_config.py` (M0-C1..C11), added with
**G-1 / REM-11**: the mud may now be Herschel-Bulkley. Two of them are the ones
to read. **M0-C9** asserts that a pair with both fluids Newtonian is
angle-*exact* (0.00%, the one case where the θ = 0 table is not an
approximation). **M0-C10** locks the corrected account of **B-1**: holding the
mud strictly Newtonian and moving only its viscosity, the angle sensitivity runs
0.53% → 11.4% → 93.1% as `m` goes 0.0063 → 0.032 → 0.127, so what makes the
shipped pair mild is the viscosity *ratio*, not the mud's rheology — and the
hazard therefore predates G-1 rather than being created by it.

⚠️ **Read those two gates with their axis caveat.** They probe a table whose
velocity axis stops at `|ū| = 10`. On the production axis (to 3000, which
FLU-06 forces) the shipped pair measures **9.8% — SIGNIFICANT, not benign**.
That had never been measured, because the production runs loaded a cached table
and a loaded table reports `not measured`. No gate changes and no computed
number moves, but `AUDIT_REPORT.md` B-1's "benign for the K-GEP-1 pair" is not
supported at the resolution the runs use. See `docs/assumptions.md` FLU-07.

| Gate | Test(s) | Status | Last green at |
|---|---|---|---|
| **M0-T1** | `test_m0_geometry.py` — 2 test(s) |  **green** | `987f8ea` |
| **M0-T2** | `test_m0_geometry.py` — 1 test(s) |  **green** | `987f8ea` |
| **M0-T3** | `test_m0_geometry.py` — 3 test(s) |  **green** | `987f8ea` |
| **M0-T4** | `test_m0_geometry.py` — 2 test(s) |  **green** | `987f8ea` |
| **M0-T6** | `test_m0_geometry.py` — 3 test(s) |  **green** | `987f8ea` |
| **M1-T1** | `test_m1_scaling.py` — 4 test(s) |  **green** | `987f8ea` |
| **M1-T2** | `test_m1_scaling.py` — 6 test(s) |  **green** | `987f8ea` |
| **M2-T1** | `test_m2_closures.py` — 2 test(s) |  **green** | `987f8ea` |
| **M2-T2** | `test_m2_closures.py` — 2 test(s) |  **green** | `987f8ea` |
| **M2-T3** | `test_m2_closures.py` — 2 test(s) |  **green** | `987f8ea` |
| **M2-T4** | `test_m2_closures.py` — 2 test(s) |  **green** | `987f8ea` |
| **M2-T5** | `test_m2_closures.py` — 2 test(s) |  **green** | `987f8ea` |
| **M3-T1** | `test_m3_tables.py` — 4 test(s) |  **green** | `987f8ea` |
| **M3-T2** | `test_m3_tables.py` — 3 test(s) |  **green** | `987f8ea` |
| **M3-T3** | `test_m3_tables.py` — 3 test(s) |  **green** | `987f8ea` |
| **M4-T1** | `test_m4_elliptic.py` — 2 test(s) |  **green** | `987f8ea` |
| **M4-T2** | `test_m4_elliptic.py` — 3 test(s) |  **green** | `987f8ea` |
| **M4-T3** | `test_m4_elliptic.py` — 5 test(s) |  **green** | `987f8ea` |
| **M4-T4** | `test_m4_elliptic.py` — 2 test(s) |  **green** | `987f8ea` |
| **M5-T1** | `test_m5_transport.py` — 2 test(s) |  **green** | `987f8ea` |
| **M5-T2** | `test_m5_transport.py` — 2 test(s) |  **green** | `987f8ea` |
| **M5-T3** | `test_m5_transport.py` — 1 test(s) |  **green** | `987f8ea` |
| **M5-T4** | `test_m5_transport.py` — 2 test(s) |  **green** | `987f8ea` |
| **M6-T1** | `test_m6_simulation.py` — 2 test(s) |  **green** | `987f8ea` |
| **M6-T2** | `test_m6_simulation.py` — 2 test(s) |  **green** | `987f8ea` |
| **M6-T3** | `test_m6_simulation.py` — 1 test(s) |  **green** | `987f8ea` |
| **M6-T4** | `test_m6_simulation.py` — 1 test(s) |  **green** | `987f8ea` |
| **M7-T1** | `test_m7_postprocess.py` — 2 test(s) |  **green** | `987f8ea` |
| **BENCH-01/02/03/05** | `test_m4_elliptic.py`, `test_m4b_twodga.py` — PF04 analytic steady states |  **green** | `987f8ea` |
| **BENCH-06** | `tests/benchmarks/kinematic_wave.py` + `test_m6_simulation.py` |  **green** | `987f8ea` |
| **BENCH-07/09** | ZF22 Table 3, via `scripts/zf22_table3.py` |  **green** | `987f8ea` |
| **BENCH-08** | ZF23 Table 3, `test_m7_postprocess.py` |  **green** | `987f8ea` |
| **BENCH-10** | `test_bench10_muskat.py` — BF25 §3.3, **added this session** |  **green** | `987f8ea` |

## Gates that assert nothing

| Gate | Why |
|---|---|
| **M0-T5** | Referenced in `docs/derivation.md` and implemented as a figure-producing script (`scripts/caliper_report.py`). **No test carries the name**, so it cannot pass and cannot regress. Without `D2DGA_BUILD_SPEC.md` (see `BLOCKED.md` BLK-5) there is no statement of what it was supposed to assert, so no test could be written for it here. |

## Not covered by the suite, and run separately

| Item | Status | When |
|---|---|---|
| **ZF22 Table 3**, all ten cases, 20 x 200, CFL 0.5 | **9 of 10 reproduce**; case 1 does not (`BLOCKED.md` BLK-6). Re-run in full this session and reproduces `output/zf22_table3.md` exactly. | 2026-09-18 |
| **K-GEP-1 production run** | **NOT re-run** (~3.3 h). Its numbers in `output/kgep1_results.md` are carried over from commit `c47c42f`; the argument that they do not move is in `REMEDIATION_SUMMARY.md` section 7 item 2, and it is an argument, not a run. | 2026-09-17 |
| **Measured-caliper run** | Never completed, before or after (~16 h). Smoke-tested only. | — |
