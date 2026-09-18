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

| Gate | Test(s) | Status | Last green at |
|---|---|---|---|
| **M0-T1** | `test_m0_geometry.py` — 2 test(s) | | |
| **M0-T2** | `test_m0_geometry.py` — 1 test(s) | | |
| **M0-T3** | `test_m0_geometry.py` — 3 test(s) | | |
| **M0-T4** | `test_m0_geometry.py` — 2 test(s) | | |
| **M0-T6** | `test_m0_geometry.py` — 3 test(s) | | |
| **M1-T1** | `test_m1_scaling.py` — 4 test(s) | | |
| **M1-T2** | `test_m1_scaling.py` — 6 test(s) | | |
| **M2-T1** | `test_m2_closures.py` — 2 test(s) | | |
| **M2-T2** | `test_m2_closures.py` — 2 test(s) | | |
| **M2-T3** | `test_m2_closures.py` — 2 test(s) | | |
| **M2-T4** | `test_m2_closures.py` — 2 test(s) | | |
| **M2-T5** | `test_m2_closures.py` — 2 test(s) | | |
| **M3-T1** | `test_m3_tables.py` — 4 test(s) | | |
| **M3-T2** | `test_m3_tables.py` — 3 test(s) | | |
| **M3-T3** | `test_m3_tables.py` — 3 test(s) | | |
| **M4-T1** | `test_m4_elliptic.py` — 2 test(s) | | |
| **M4-T2** | `test_m4_elliptic.py` — 3 test(s) | | |
| **M4-T3** | `test_m4_elliptic.py` — 5 test(s) | | |
| **M4-T4** | `test_m4_elliptic.py` — 2 test(s) | | |
| **M5-T1** | `test_m5_transport.py` — 2 test(s) | | |
| **M5-T2** | `test_m5_transport.py` — 2 test(s) | | |
| **M5-T3** | `test_m5_transport.py` — 1 test(s) | | |
| **M5-T4** | `test_m5_transport.py` — 2 test(s) | | |
| **M6-T1** | `test_m6_simulation.py` — 2 test(s) | | |
| **M6-T2** | `test_m6_simulation.py` — 2 test(s) | | |
| **M6-T3** | `test_m6_simulation.py` — 1 test(s) | | |
| **M6-T4** | `test_m6_simulation.py` — 1 test(s) | | |
| **M7-T1** | `test_m7_postprocess.py` — 2 test(s) | | |
| **BENCH-01/02/03/05** | `test_m4_elliptic.py`, `test_m4b_twodga.py` — PF04 analytic steady states | | |
| **BENCH-06** | `tests/benchmarks/kinematic_wave.py` + `test_m6_simulation.py` | | |
| **BENCH-07/09** | ZF22 Table 3, via `scripts/zf22_table3.py` | | |
| **BENCH-08** | ZF23 Table 3, `test_m7_postprocess.py` | | |
| **BENCH-10** | `test_bench10_muskat.py` — BF25 §3.3, **added this session** | | |
