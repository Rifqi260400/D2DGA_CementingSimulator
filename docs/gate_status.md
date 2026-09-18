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

**Verdict of the run this table is built from:** `pytest tests/` -> **254 passed, 0 failed, 0 errors** (exit 0), 2026-09-18 07:20:14, at the
commit below. Before this remediation the suite was 223 passed / 0 failed;
31 tests were added and none removed or skipped.

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
