# Remediation summary

Work done against `AUDIT_REPORT.md`, from commit `be4cd99` to `987f8ea`.
Full suite at the end: **254 passed, 0 failed, 0 errors**, run 2026-09-18. Full detail in `docs/remediation_log.md`; unresolved items in
`BLOCKED.md`.

**No tolerance, test, or gate was weakened.** Two tests had a *claim* corrected
after the claim turned out to be wrong — both are written up with the
derivation, and both ended up asserting *more* than before, not less. They are
listed explicitly in §7.

---

## 1. Findings resolved, with the root cause of each

| ID | Root cause in one line |
|---|---|
| **A-3** | The closure-table range guard was a `UserWarning` and nothing else, so `-W ignore::UserWarning` erased the only evidence of extrapolation; and the `umag` axis started at 0.02 while `\|ū\| ≥ 0`. Record is now unsuppressible state on the table; axis anchored at 0. |
| **B-2** | `derivative_bounds` — which supplies every LLF wavespeed on the tabulated path — never range-checked at all, so it was silent *by construction*. Now records through the same path. |
| **A-2** | Not a bug: `t_br` is a level set in a diffusive tail, so it is `O(Δξ^{1/2})` where `η_E` is `O(Δξ)`. The *claim* that `t_br@0.01` is a converged number was wrong. Achievable accuracy derived, gate added, error bar now printed. |
| **A-4.1** | "Half the suite is saturated" was an impression; now a measurement — five of ten ZF22 cases move `η_E` by < 0.01 under a **doubling** of `b`. |
| **A-4.3** | The whole Herschel–Bulkley path had no external benchmark because every benchmark in the repo is Newtonian. Added **BENCH-10** from BF25 §3.3. |
| **A-5** | `BENCH-09` established mesh convergence on `t_br` only and refined both directions together. Refining one at a time shows neither converges and they cancel along the diagonal it used. |
| **B-1** | The table stores the `θ = 0` slice and the docstring wrongly called that a *reduction*. It is an approximation, and for a yield-stress pair it is `O(1)`. Build-time measurement + report added. **Upgraded MINOR → MAJOR.** |
| **B-3** | Stale prose in two places said the opposite of NUM-26 — it survived the NUM-26 fix. Corrected; behaviour pinned by a test. |
| **B-4** | Duplicate `CONV-03` and `NUM-02` register rows whose content contradicted the authoritative ones. Settled from the primary sources, both marked superseded in place. |
| **B-5** | Two records of the same run disagreed; the register rows predated the NUM-26/NUM-29 fixes. Settled by re-running all ten ZF22 cases. |
| **B-6** | `record_every=0` raised `ZeroDivisionError`. |
| **Q-1** | Confirmed and bounded; affects no reported number, but the ZF23 tail statistics are computed where the interpolant is worst. Caveated. |
| **Q-2** | BF25 (2.23) prints two forms joined by "=" that are not equal. The first is correct. Now asserted in code. |
| **R-1** *(new, found during remediation)* | `breakthrough_at` interpolated the *recorded* history, so a reported physical quantity depended on `record_every` — 31 % on ZF22 case 1. The scalar history is now recorded every step. |

## 2. Findings not reproduced, and what was observed instead

**A-1 — "the outflow face imports displacing fluid; the unfixed mirror of
NUM-29" — NOT REPRODUCED as a defect.**

The evidence reproduces: the `ξ = Z` face does run backwards, and the audit's
measurements are correct. The interpretation does not. For these flux functions
the upper concave envelope of `f = q₀ + b𝓘₃` on `[0,1]` is the straight chord,
so **the intermediate concentrations that make `f` negative never occur in the
exact solution** — the spurious flux is manufactured by the first-order
smearing. Measured against the exact entropy solution:

| `b_ZF22` | outlet import, `n_xi` = 50 → 400 | rate |
|---|---|---|
| 1000 | 0.0363 → 0.0181 → 0.0091 | **1.00** |
| 100 | 0.00375 → 0.00188 → 0.00094 → 0.00047 | **1.00** |
| 10 | identically 0 (`f ≥ 0` everywhere) | — |

It is `O(Δξ)` and vanishes. The net error in `η_E` is ~100× smaller than the
gross boundary traffic the audit quoted (the audit itself flagged medium
confidence there). The zero-gradient ghost is also the **only** closure that
preserves a uniform state exactly — injecting a `c̄ = 0` ghost breaks
`test_uniform_states_are_preserved_exactly`. **The boundary condition was left
unchanged.** What was missing was observability, so the volume is now
accumulated and printed.

## 3. Items in `BLOCKED.md`, and what each needs from you

| ID | Needs |
|---|---|
| **BLK-1** | **Centralizer record for K-GEP-1.** The largest sensitivity in the model: `e_gauge` 0.1 → 0.5 collapses narrow-side velocity fraction by **130×**. No efficiency number means much until this is fixed. |
| **BLK-2** | Pump schedule. (Less critical: FLU-06 measures the displacement as *not* pump-rate controlled.) |
| **BLK-3** | **Real mud properties.** The source paper's "drilling fluid" is water. `η_E = 0.994` is an upper bound, not a prediction. |
| **BLK-4** | Casing ID/weight/grade, to fix the caliper's +0.24 – 0.43 in bias. It propagates as `H³` in every flux. |
| **BLK-5** | **`D2DGA_BUILD_SPEC.md` itself** — it is not in the repository or in git history. The gate list is reconstructed from surviving `M<n>-T<k>` markers, so a gate the spec required but never marked is invisible. One gap is already visible: **M0-T5 has no test at all**. |
| **BLK-6** | ZF22 case 1. Not fixed; see §7. |
| **BLK-7** | BF25 (3.13)'s stability range is degenerate as printed; resolved reading implemented, alternative behind a flag, difference quantified. |
| **BLK-8** | Whether an extrapolated derivative bound ever fell *below* the true slope — not measured; cannot recur after the A-3 axis fix. |

## 4. Places where a published paper is wrong

Stated prominently, as required. Each is derived, not asserted.

### 4.1 BF25 (2.27) is wrong three ways — and two other papers confirm two of them

The printed form is
`𝓘₃ = c̄²(1−c̄)³[4mc̄ + 3(1−c̄²)] / (12[mc̄³+1−c̄³])`. The correct form is

```
𝓘₃ = − c̄²(1−c̄)³[4mc̄ + 3(1−c̄)] / (12 √m [mc̄³ + 1 − c̄³])
```

* **Numerator.** `3(1−c̄²)` should be `3(1−c̄)`. **ZF22 (4.26) and BCF25 (23)
  both print `3(1−c̄)`** — two independent papers, so this is confirmed without
  reference to this code.
* **Missing `√m`** in the denominator, although the neighbouring (2.24) and
  (2.25) both carry their `m^{±½}`.
* **Sign.** `𝓘₃ < 0` on `0 < c̄ < 1`. Confirmed three ways: integrating the
  velocity profile of (2.9)–(2.10) directly (which also reproduces (2.13),
  (2.24), (2.25), (2.26) exactly); the numerical AL solver; and **ZF22 (4.25)
  read with ZF22's own conventions** — at its (4.24) ZF22 writes
  `ρ = ρ₂ + (1−c̄)Δρ`, i.e. `Δρ = ρ₁ − ρ₂`, so `Δρ/F² = −b_ZF22` — which makes
  the two papers agree to **3.9e-16** over 144 states.

Conversions, verified to 1e-13: `BCF25 (23) = −√m · MASTER`,
`ZF22 (4.26) = −(6/√m) · MASTER`.

### 4.2 BF25 (2.23) prints two expressions joined by "=" that are not equal

They differ by `−c̄³√m (c̄−1)²(c̄+2) / (6(mc̄³ − c̄³ + 1))`, strictly negative on
`0 < c̄ < 1`. The **first** is correct, established without using (2.23) at all.
This is a *second, separate* defect in the same equation; the project's own
`CONV-03` recorded only (2.27).

### 4.3 BF25 (3.13)'s stability range is degenerate as printed

"For an unstable regime, `Δw(c₀) > 0` for all `c₀ ∈ [0,1]`" can never classify
anything: at `c₀ → 1`, `w_f → q₀'(1) + b𝓘₃'(1) = 0` for **every** pair and every
`b`, while `w_finger → 1`, so `Δw(1⁻) → +1` unconditionally. All ten ZF22 cases
come out "partial penetration", including case 4 at `b = 1000`, which displaces
as a near-piston. Resolved per brief §5b — see `BLOCKED.md` BLK-7.

### 4.4 Already recorded by the build, re-confirmed here

BCF25 (30)/(31)/(34)/(35) omit the mesh spacings in the buoyancy flux (NUM-15);
BCF25 (18)/(20) carry the azimuthal `Ψ`-differences with the wrong sign,
inconsistent with its own (11)/(14) (NUM-16); BF25 (2.18) carries a `1/r_a` on
the buoyancy term that the pressure elimination does not produce (CONV-04);
BF25 (2.21)/BCF25 (5) print the buoyancy flux so `r_a` cancels, contradicting
BF25 (2.20) (CONV-10); ZF22 Table 1 case 9's `Q̂₀` is a typo (CONV-08).

## 5. Places where the build spec misread a paper

`D2DGA_BUILD_SPEC.md` **is not in this repository** (BLK-5), so this section can
only relay what `docs/assumptions.md` attributes to it, not verify it:
GEO-03 (radius/diameter conflation), GEO-04 (a wall-amplitude sweep half of
which is geometrically infeasible), GEO-07 (`δ/π ≈ 0.34`, unattainable since
`δ/π < 1/π` for any `r̂_i > 0`), CONV-07 (ZF22 Table 2's `Re` treated as exact
when it holds rounded design targets). **None of these was independently
re-checked in this remediation**, because the document they refer to is not
available.

## 6. Gate status, before and after

Every "after" state below was produced by an **actual run in this session**.
`docs/gate_status.md` carries the per-gate table.

| | Before (audit, commit `be4cd99`) | After |
|---|---|---|
| Test suite | 223 passed, 0 failed | **254 passed, 0 failed, 0 errors** (exit 0) — 31 tests added, none removed or skipped |
| M0 | pass (M0-T5 has no test) | pass (M0-T5 **still** has no test — BLK-5) |
| M1 | pass | pass |
| M2 | pass | pass, + Q-2 and the cross-paper CONV-03 checks |
| M3 | pass with reservations | pass; range guard unsuppressible (A-3), `derivative_bounds` guarded (B-2), angle sensitivity measured per table (B-1) |
| M4 | pass | pass |
| M5 | pass | pass, + the A-1 first-order lock and the A-2 achievable-accuracy gate |
| M6 | conditional | conditional — unchanged: case 1 still does not reproduce (BLK-6), and R-1 fixed |
| M7 | pass | pass |
| **BENCH-10** | did not exist | **new** — first external analytic check on the Herschel–Bulkley path |
| ZF22 Table 3 | 9/10 | **9/10, re-run in full this session**, reproducing `output/zf22_table3.md` exactly |
| PF04 analytic | pass | pass, re-run this session |
| K-GEP-1 production | `η_E = 0.9946` at 16 × 80 | **not re-run** — see §7 |

**No previously green gate regressed.** The only failure during the work was a
test I had just written (B-1) whose bound was calibrated on the wrong fluid
pair; correcting it is §7.

## 7. What I am least confident about

In order.

1. **ZF22 case 1 is still not reproduced, and the stop condition is not met.**
   The brief requires all ten cases. Nine reproduce; case 1 does not, and I did
   not fix it. What I did establish is that it is probably not fixable as
   stated: `η_E` has no mesh limit in either direction (increments *grow*
   azimuthally, stay flat axially, and point opposite ways), and BF25's own
   Muskat criterion classifies it as the only unstable case of the ten. I
   believe a fixed-time `η_E` is not a well-posed target for this case — but
   that is an argument, and ZF22 published a number, so I may be wrong. The one
   untested structural difference is azimuthal periodicity vs the half-annulus
   symmetry (NUM-05/Q6), which I did not try because it is a new capability.

2. **The K-GEP-1 production result was not re-run.** At ~3.3 h per run it did
   not fit alongside the rest of this work. I assert that its numbers do not
   move, on the basis that (a) A-3's extrapolation changed the closures by
   nothing at six significant figures, (b) A-1's boundary condition was left
   unchanged, (c) R-1 did not affect it since it already ran at
   `record_every = 1`, and (d) B-1's angle sensitivity at that table's operating
   range (`|ū| ≤ 75`) measures 0.4 %. **That is reasoning, not a run.** The
   `output/kgep1_results.md` numbers are carried over, and the new diagnostics
   (range report, outlet import, angle sensitivity, `t_br` error bars) will only
   appear when it is next run.

3. **B-1 may still be under-measured.** I found it is `O(1)` for a yield-stress
   pair only because a test I calibrated on K-GEP-1 failed on a different pair.
   The probe now runs per table, but it samples three concentrations and three
   speeds at one `g_b`. A pair worse than the one I happened to test would not
   surprise me.

4. **Two tests had a claim corrected mid-work** (brief §1 requires these be
   called out):
   * `test_bench10_leading_edge_criterion_reproduces_bf25_m_equals_1_5` asserted
     `Δw(0⁺) = m − 3/2` for all `m` and failed at `m = 3, 5`. That was the test
     being wrong: above `M₃^min = 3/2` the leading wave is a *shock*, not the
     characteristic, so the equality does not hold. The corrected test asserts
     **more**: the dispersive/shock transition at `3/2` exactly, the equality
     below it, and the sign flip for every `m`.
   * `test_b1_angle_assumption_is_bounded` asserted bounds calibrated on
     K-GEP-1 and failed on the yield-stress pair. The corrected test asserts the
     measured `O(1)` deviations and the non-monotonicity in `|ū|`.

   Neither is a loosened tolerance — both now assert strictly more than the
   version that failed. But both are cases where I wrote the claim before I had
   measured it, and I would not have caught either without the suite.

5. **`M0-T5` asserts nothing.** It exists only as a figure-producing script. I
   could not write a test for it because, without the build spec (BLK-5), there
   is no statement of what it was supposed to assert.

6. **BLK-8 is unmeasured.** I established that `derivative_bounds` extrapolated
   silently; I did not establish that no extrapolated bound ever fell *below*
   the true slope, which is the direction that breaks monotonicity. The indirect
   evidence is strong and the situation cannot recur after the axis fix, but
   indirect is not measured.
