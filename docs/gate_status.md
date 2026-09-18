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

| Gate | What it asserts | Status | Last green at |
|---|---|---|---|
| M0-T1 | (see tests) | pending run | — |
| M0-T2 | (see tests) | pending run | — |
| M0-T3 | (see tests) | pending run | — |
| M0-T4 | (see tests) | pending run | — |
| M0-T6 | (see tests) | pending run | — |
| M1-T1 | (see tests) | pending run | — |
| M1-T2 | (see tests) | pending run | — |
| M2-T1 | (see tests) | pending run | — |
| M2-T2 | (see tests) | pending run | — |
| M2-T3 | (see tests) | pending run | — |
| M2-T4 | (see tests) | pending run | — |
| M2-T5 | (see tests) | pending run | — |
| M3-T1 | (see tests) | pending run | — |
| M3-T2 | (see tests) | pending run | — |
| M3-T3 | (see tests) | pending run | — |
| M4-T1 | (see tests) | pending run | — |
| M4-T2 | (see tests) | pending run | — |
| M4-T3 | (see tests) | pending run | — |
| M4-T4 | (see tests) | pending run | — |
| M5-T1 | (see tests) | pending run | — |
| M5-T2 | (see tests) | pending run | — |
| M5-T3 | (see tests) | pending run | — |
| M5-T4 | (see tests) | pending run | — |
| M6-T1 | (see tests) | pending run | — |
| M6-T2 | (see tests) | pending run | — |
| M6-T3 | (see tests) | pending run | — |
| M6-T4 | (see tests) | pending run | — |
| M7-T1 | (see tests) | pending run | — |
