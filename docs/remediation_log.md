# Remediation log

One entry per finding in `AUDIT_REPORT.md`. Each records: whether the finding
reproduced, the ground truth **from the papers** (section and equation), the
root cause, the fix, the test that would have caught it, and the gate re-run.

**No tolerance, test, or gate was weakened in this work.** Where a quantity
turned out not to be convergent at the accuracy previously implied, the
achievable accuracy is *derived* below and the claim is corrected — the
tolerance is not.

> **`D2DGA_BUILD_SPEC.md` is not in this repository.** It was supplied in the
> conversation that produced the code and never committed. Every reference to it
> below is therefore to `docs/assumptions.md`'s citations of it, not to the
> document. Consequences are recorded in `BLOCKED.md`.

---

## A-1 — outflow face imports displacing fluid — **NOT REPRODUCED as a defect**

**Audit claim.** The `ξ = Z` face carries a negative flux and so imports fluid 2
from the ghost region; measured −10.1 × Q instantaneously and 3.5 % of the job
entering that way on a 10 × 50 mesh; read as "the unfixed mirror of NUM-29".

**Evidence reproduced, interpretation refuted.** The negative flux is real and
the measurement stands. The interpretation does not.

*Ground truth.* BF25 §3.1 (3.8) plus Oleinik: for the Riemann problem
`c̄ = 1` below, `c̄ = 0` above (so `u_l > u_r`), the entropy solution is the
**upper concave envelope** of `f(c̄) = q₀ + b𝓘₃` on `[0,1]`. For the ZF22 and
K-GEP-1 flux functions that envelope is the straight chord from `(0,0)` to
`(1,1)` — the audit itself measured `max f(c̄)/c̄ = 1.00000` for K-GEP-1. **The
intermediate concentrations that make `f` negative therefore never occur in the
exact solution.** The spurious flux is manufactured entirely by the first-order
smearing of the front, so it must vanish under refinement, and it does.

Measured on a concentric annulus (`tests/test_m5_transport.py`, and the sweep in
this session), against the exact entropy solution built from the envelope:

| `b_ZF22` | `n_xi` | outlet import (volumes) | rate | `η_E` error vs exact | rate |
|---|---|---|---|---|---|
| 1000 | 50 / 100 / 200 | 0.0363 / 0.0181 / 0.0091 | **1.00, 1.00** | 2.2e-3 / 1.1e-3 / 5.9e-4 | 0.95, 0.95 |
| 100 | 50 / 100 / 200 / 400 | 0.00375 / 0.00188 / 0.00094 / 0.00047 | **1.00, 1.00, 1.01** | 6.9e-3 / 3.7e-3 / 2.0e-3 / 1.1e-3 | 0.90, 0.90, 0.90 |
| 10 | any | **0** | — | 1.7e-2 … 3.0e-3 | 0.81–0.84 |

At `b_ZF22 = 10` the flux function is non-negative on all of `[0,1]`, no face can
run backwards, and the import is identically zero — which pins the mechanism to
the sign of `f`, not to the boundary.

**Why the boundary condition is left unchanged.** Zero-gradient replication is
the only closure that preserves a uniform state exactly: `F_out(c*, c*)` must
equal the interior flux or a constant field evolves. Injecting a `c̄ = 0` ghost —
the obvious "no fluid 2 above" alternative — breaks both
`test_uniform_states_are_preserved_exactly` and the new A-1 lock. Verified by
doing it.

**What was actually missing: observability.** An artefact nobody measures becomes
permanent (the lesson of A-3). `TransportSolver.outlet_import` now accumulates
the volume and `scripts/kgep1_run.py` prints it against capacity, with the
expectation ("halves when `n_xi` doubles") stated in the output.

**Tests.** `test_a1_outlet_import_vanishes_at_first_order` (a genuine boundary
defect plateaus instead of halving — verified to fail under an injected `c̄ = 0`
ghost) and `test_a1_outlet_import_is_zero_when_the_flux_never_goes_negative`.

**Correction to the audit:** the audit's "3.5 % of the job" is gross boundary
traffic, not a bias on `η_E`. The error in `η_E` is ~100× smaller and converges.
The audit flagged its own confidence here as "medium on how much of the 3.5 %
survives"; the answer is ~1 %.

---

## A-3 — closure-table range guard suppressed — **REPRODUCED, fixed at the root**

**Reproduced exactly.** Re-running the K-GEP-1 configuration with warnings
enabled: 76 `ClosureTableRangeWarning`s in the first 400 of 85 474 steps, on the
**low** end of the `umag` axis (requests to 1.08e-3 against a floor of 0.02).
The production logs contain zero warning lines.

**Two root causes, both fixed.**

1. *The guard was only a warning.* A `-W ignore::UserWarning`, a
   `warnings.simplefilter`, or a pytest config erases it, and then the fact that
   a result rests on linear **extrapolation** leaves no trace. The record is now
   state on the `ClosureTable` (`n_range_queries`, `n_range_points_out`,
   `worst_excursion`, `range_report()`), updated unconditionally on every query.
   `scripts/kgep1_run.py` prints it whether or not it is clean.
2. *The axis floor was wrong.* `|ū| ≥ 0` always, but the `umag` axis started at
   0.02. It is now anchored at 0, so an under-range query is impossible.

**Measured impact of the extrapolation that did occur: none that is visible.**
The gap-scale solve at `umag = 0`, `1e-3` and `0.02` agrees to **six significant
figures** on the K-GEP-1 pair (`𝓘₁ = 53.9673`, `q₀ = 0.148652`,
`𝓘₃ = −0.266212` at `c̄ = 0.1`, identical at all three). Re-measured on the real
configuration, 315 points of 25 899 552 point-queries fell outside, worst
excursion 0.0189. **No previously reported number moves.**

**Tests.** `test_a3_range_record_survives_warning_suppression` runs the query
under `warnings.simplefilter("ignore")`, asserts no warning was raised **and**
that the record is intact. Fails against the old code with `AttributeError`.

---

## B-2 — `derivative_bounds` had no range check — **REPRODUCED, fixed with A-3**

`ClosureTable.derivative_bounds` supplies every LLF wavespeed on the tabulated
path (29 probe nodes per face, every step) and extrapolates from the same axes
as `__call__`, but never checked them — silent **by construction**, not by a
flag. It now records through the same `_record_range`.

**Test.** `test_b2_derivative_bounds_records_out_of_range_too`. Fails against the
old code with `AttributeError`.

**Not measured, and still open:** whether an extrapolated derivative bound ever
fell *below* the true slope, which is the direction that would break
monotonicity. The audit said the same. Carried to `BLOCKED.md`.

---

## B-3 — two places stated the opposite of NUM-26 — **REPRODUCED, corrected**

`d2dga/transport.py`'s "Monotonicity" docstring and `docs/assumptions.md`
NUM-17 both read *"endpoint evaluation is sufficient — no interval maximum is
needed"*, contradicting NUM-26 in the same two files, which measures a
maximum-principle violation of `c̄ = −0.0089` under endpoint-only wavespeeds.
The prose was stale — it survived the NUM-26 fix — and it is exactly the
justification a maintainer would cite to flip `wavespeed` back to `"endpoints"`,
which is still selectable.

Both texts corrected. **Test:** `test_b3_interval_wavespeed_is_the_default` pins
the behaviour the prose contradicted; the existing
`test_num26_endpoint_wavespeeds_break_monotonicity` remains the substantive
guard.

---

## B-4 — duplicate register IDs with contradictory content — **REPRODUCED, settled from the papers**

`docs/assumptions.md` carried two `CONV-03` rows and two `NUM-02` rows. The
stale `CONV-03` printed the master `𝓘₃` **without its leading minus** and gave
the conversions as `BCF25 (23) = MASTER·√m`, `ZF22 (4.26) = MASTER·6/√m` — also
without their minus signs. Followed literally it reintroduces the sign error
that CONV-03 exists to record.

**Ground truth, derived rather than chosen (brief §3.3).**

*The numerator.* ZF22 (4.26) prints `c̄²(1−c̄)³[4mc̄+3(1−c̄)] / (2m[mc̄³+1−c̄³])`
and BCF25 (23) prints `c̄²(1−c̄)³[4mc̄+3(1−c̄)] / (12[mc̄³+1−c̄³])`. **Both have
`3(1−c̄)`.** BF25 (2.27)'s `3(1−c̄²)` is therefore a transcription error
confirmed by two independent papers, not merely by this code's derivation.

*The sign.* Checked on the physical quantity, not on anyone's printed `𝓘₃`.
BF25 (2.21) ξ-component and ZF22 (4.25) ξ-component must give the same buoyant
flux. Read with **ZF22's own conventions, both stated in ZF22's own text** —
at (4.24) ZF22 writes `ρ = ρ₂ + (1−c̄)Δρ` against `ρ = (1−c̄)ρ₁ + c̄ρ₂`, so
`Δρ = ρ₁ − ρ₂`, the *same* ordering as BF25, giving `Δρ/F² = −b_ZF22`; and
`η₂ = μ̂₂/μ̂₁ = 1/m` — the two papers agree to **3.9e-16 over 144 states**
spanning `m ∈ [0.2, 5]`, `b_ZF22 ∈ [−50, 1000]`, `c̄`, `H`. Reading ZF22's
`Δρ` with the opposite sign flips the result, which is the mistake the stale row
encoded.

*The conversions*, verified to 1e-13 over `c̄ ∈ [0.02, 0.98]` and
`m ∈ {0.006, 0.2, 1, 5, 160}`:

```
BCF25 (23)  = -sqrt(m)     * MASTER
ZF22 (4.26) = -(6/sqrt(m)) * MASTER
```

— the signs the authoritative row already carried. Both stale rows are marked
**SUPERSEDED** in place rather than deleted, so the contradiction stays visible.

**Tests.** `test_b4_conv03_buoyant_flux_agrees_with_zf22_4_25` and
`test_b4_conv03_conversion_factors_carry_their_minus_signs`. Verified to **fail**
when `script_I3` is replaced by the form the stale row describes.

---

## B-6 — `record_every=0` raised `ZeroDivisionError` — **REPRODUCED, fixed**

`d2dga/simulation.py` `if n % record_every == 0`. Zero now means endpoints only,
which is what it reads as. **Test:** `test_b6_record_every_zero_means_endpoints_only`,
which also asserts the resulting field is bit-identical to a fully recorded run.

---

## A-4.3 — the Herschel–Bulkley path had no published benchmark — **REPRODUCED, closed**

Every benchmark in the repository was Newtonian: ZF22 Table 3, ZF23 Table 3 and
PF04 all are, so the M3 closure table, the Picard elliptic and the NUM-13
mobility floor — the path the K-GEP-1 result runs on — were tested against no
external analytic statement.

**BF25 §3.3 is the one closed-form result in the six papers that applies to an
arbitrary pair.** Its (3.12),

```
w_finger(c0) = I1(1)/I1(c0) + b I1(1) [ I2(c0)/I1(c0) + c0 - 1 ]
```

is built from `𝓘₁` and `𝓘₂` alone. Implemented verbatim in
`tests/benchmarks/bf25_muskat.py`; gate **BENCH-10** in
`tests/test_bench10_muskat.py`.

Three results, none of them fitted:

1. **`𝓘₂(1) = 0`** — BF25 uses it to reach (3.12) and it was **not** among the
   endpoint identities this repository checked. Exactly zero analytically and
   below 1e-9 on a tabulated HB pair.
2. **The b = 0 reduction recovers BF25's own `M₃^min = 1.5`.** At `b = 0`,
   `w_finger(0⁺) = 𝓘₁(1)/𝓘₁(0) = m` and the front's leading wave travels at
   `q₀'(0) = 3/2`, so `Δw(0⁺) = m − 3/2` and the finger penetrates exactly when
   `m > 3/2`. BF25 §3.2 gives 3/2 as `M₃^min`, attributing it to Lajeunesse
   et al. (1999) from a completely different argument. Reproduced to 2.3e-9 at
   `m = 1.5`. This is what pins the reading of (3.13) to the paper.
3. **K-GEP-1 is Muskat-stable on the tabulated HB path**, `Δw(0⁺) < −10`.

---

## A-5 — ZF22 case 1 is not mesh-converged in `η_E` — **REPRODUCED, and explained**

Reproduced exactly: `η_E` 0.3657 (20 × 200) → 0.4124 (40 × 400), **+12.8 %**,
while `t_br` moves only −1.3 %. `BENCH-09`'s claim "(i) it is mesh-converged"
rests on `t_br` alone and does not hold for `η_E`, which is the quantity the
failure is stated in and which `assumptions.md` itself calls "the cleanest
comparison … threshold-free, unlike `t_br`".

**New: which mesh direction drives it.** BENCH-09 refined both together, so it
could not say. Refining one at a time (this session):

**Azimuthal refinement at fixed `n_xi = 200`:**

| `n_phi` | `η_E` | increment |
|---|---|---|
| 20 | 0.3657 | — |
| 40 | 0.3939 | **+0.0282** |
| 80 | 0.4773 | **+0.0834** |

**Axial refinement at fixed `n_phi = 20`:**

| `n_xi` | `η_E` | increment |
|---|---|---|
| 200 | 0.3657 | — |
| 400 | 0.3508 | −0.0149 |
| 800 | 0.3340 | −0.0168 |

**Neither direction converges, and they push opposite ways.** Azimuthally the
increments *grow* — they nearly triple on the second halving of `Δφ`. Axially
they are flat (−0.0149, −0.0168) where a first-order scheme should halve them.
A converging sequence has shrinking increments; neither of these does.

**So the "mesh convergence" BENCH-09 observed at 20 × 200 → 40 × 400 was a
cancellation of two individually non-convergent trends**, one upward and one
downward, that happened to nearly offset in `t_br`. That is the strongest
result of this remediation: `η_E` for ZF22 case 1 has **no mesh limit**, and it
was only ever checked along the one diagonal where the two errors cancel.

This confirms the BENCH-10 reading independently of BF25: an analytically
Muskat-unstable base state produces a finger whose width is set by the mesh, so
`η_E` at a fixed time is not a converged functional of it. `BENCH-09`'s claim
"(i) it is mesh-converged" is **false**.

⚠️ The `t_br` column of that study is **discarded**: it was run with
`record_every = 200`, and finding **R-1** below shows `breakthrough_at` depended
on exactly that setting. `η_E` is computed from the final concentration field
and is unaffected, and `η_E` is the quantity A-5 is about.

**New: an analytic reason, from the paper.** BF25 §3.3's Muskat criterion
(BENCH-10 above) classifies **case 1 as the only Muskat-unstable case of the
ten**: `Δw(0⁺) = +0.97` against `−1.96 … −162.9` for the other nine, and case 1
is the only case with `b < 0`. A Muskat-unstable front has no mesh-independent
finger width, so `η_E` at a fixed time is **not a converged functional of the
mesh** — which is precisely the behaviour measured. That is an analytic
statement from BF25, with no solver involved, and it is the piece BENCH-09 was
missing.

It does **not** close the gap to ZF22's 0.66, and it is not offered as doing so.
See `BLOCKED.md`.

---

## A-2 — `t_br` is mesh-limited — **REPRODUCED; the claim is corrected, no tolerance is touched**

Reproduced exactly: against the exact rarefaction, the `c̄ = 0.01` front runs
+17.5 % fast at `n_xi = 80` (the K-GEP-1 production mesh), +11.0 % at 200 (the
ZF22 mesh), +7.5 % at 400.

**Derivation of the achievable accuracy** (brief §1 requires this before any
claim about accuracy changes). The LLF update has modified equation

```
c_t + f(c)_ξ = ∂_ξ[ D_num ∂_ξ c ] + O(Δξ²),   D_num = ½ α Δξ (1 − ν)
```

with `ν` the Courant number. A front spreads diffusively to width
`W ≈ 2√(D_num t) = O(√Δξ)`. Hence two different rates for two different kinds
of quantity:

* a **level set in the diffusive tail** — `t_br` at threshold 0.01 — is
  displaced by `O(W)`, so its error is `O(Δξ^{1/2})`: **half order**;
* an **integral of `c̄`** — `η_E` — is insensitive to the tail width at leading
  order because the spreading is conservative: **first order**.

Both measured, and they agree with the derivation:

| quantity | measured order |
|---|---|
| `t_br` at `c̄ = 0.01` (exact rarefaction) | **0.45 – 0.55** |
| `L1` error of the whole profile | 0.6 – 0.8 |
| `η_E` vs the exact entropy solution (A-1 sweep, `b_ZF22` = 100 / 1000) | **0.90 / 0.95** |

**Consequence, which is the actual finding.** Halving the error in `t_br@0.01`
costs **four times** the cells. At `n_xi = 80` that error is ~18 %; reaching 1 %
would need ~300× the cells — not affordable, now or later. **`t_br@0.01` is
therefore not an economically converged quantity and must not be quoted to three
decimals.** `t_br@0.5` (~3 % at `n_xi = 80`) and `η_E` are.

**What changed.** Nothing was relaxed. `scripts/kgep1_run.py` now prints each
`t_br` with its half-order discretisation bracket and says so in the output.
`NUM-21`'s three-threshold reporting was already right; what was missing was the
error bar.

**Test.** `test_a2_threshold_front_is_half_order_and_integrals_are_first_order`
asserts the measured order is in `(0.35, 0.75)` — i.e. genuinely half order and
*not* first — that it is strictly worse than the `L1` order of the same solution,
and that the error at `n_xi = 100` exceeds 10 %, so the size cannot silently
drift.

**Still open, and recorded rather than closed:** *no completed mesh-refinement
run of K-GEP-1 exists.* Both quoted runs are 16 × 80; the only finer attempt
(16 × 120) reached `t/Z = 0.025` before the container was reclaimed. At ~3.3 h
per run at 16 × 80, a 16 × 160 run is ~13 h — beyond one container window. The
half-order result above bounds the `t_br` error without needing that run; it does
**not** bound the 2-D coupling error, which only a K-GEP-1 refinement can.

---

## B-1 — the closure table assumes `ū ∥ G̃_b` — **REPRODUCED, and it is worse than the audit found: upgraded MINOR → MAJOR**

`d2dga/gapscale/tables.py`'s module docstring claimed the angle axis "collapses"
for a vertical well because `G̃_b` is purely axial. **That is wrong**: `G̃_b`
having one component does not collapse `θ = ∠(ū, G̃_b)`, because `ū` is still a
2-vector and `v̄ ≠ 0` whenever the front is tilted — i.e. in every case of
interest. BF25 (2.8)–(2.10) make the gap-scale stress a 2-vector and
`η_k = η_k(|τ_k|)`, so all four closures depend on that angle.

What the table stores is the `θ = 0` slice; that is an **approximation**, not a
reduction. Corrected in the docstring, with the measured size (K-GEP-1 pair,
`g_b = 27.9`, `H = 1`, deviation of the `θ = 90°` solve from the stored `θ = 0`):

| `\|ū\|` | `𝓘₁` | `𝓘₂` | `q₀` | `𝓘₃` |
|---|---|---|---|---|
| 0.5 | +0.0 % | +0.0 % | +0.0 % | −0.0 % |
| 5 | +0.1 % | +0.1 % | +0.0 % | −0.3 % |
| 50 | +0.6 % | +1.4 % | +0.1 % | **−2.6 %** |

The production run reached `|ū| = 75`, so ~3 % on `𝓘₃` is the operating error
**for that pair**.

### The audit under-measured this. It is O(1) for a general pair.

The audit measured only the K-GEP-1 pair, which is a benign corner in two
independent ways: its displaced fluid is **water** (Newtonian, no yield stress),
and `b·𝓘₁ ≈ 1500` makes buoyancy dominate the stress so the angle barely enters.
Re-measured on a pair with a **yield stress in both fluids**
(`κ = 0.45/0.25`, `n = 0.6/0.8`, `τ_Y = 0.30/0.12`, `g_b = 25`), the `θ = 90°`
closures differ from the stored `θ = 0` ones by:

| `c̄` | `\|ū\|` | `𝓘₁` | `𝓘₂` | `q₀` | `𝓘₃` |
|---|---|---|---|---|---|
| 0.25 | 5 | +4.9 % | +40.7 % | +8.7 % | −17.2 % |
| 0.75 | 5 | **+42.5 %** | **+134.0 %** | +7.1 % | **−93.9 %** |
| 0.25 | 50 | +3.6 % | +41.7 % | +8.0 % | **−106.5 %** |
| 0.75 | 500 | −0.2 % | +1.4 % | +0.2 % | −0.0 % |

**These are O(1) errors on the closures, not corrections.** The mechanism is the
yield stress: whether material yields at all is set by `|τ|`, and `|τ|` depends
on the angle, so the unyielded fraction — and hence `𝓘₁`, an integral of the
fluidity — moves sharply. The dependence is also **non-monotone in `|ū|`**: it
peaks where the pressure-driven and buoyancy-driven stress scales are
comparable and falls again once `|ū|` dominates, so a single "small `|ū|`"
caveat would not have caught it.

**Why no gate caught it: every benchmark in this repository is Newtonian,** and
for a Newtonian pair there is *no* angle dependence at all — `η` is constant, so
rotating `ū` relative to `G̃_b` changes the stress direction and nothing else.
Measured: angle sensitivity `< 1e-6` for a Newtonian pair against `> 0.5` for
the yield-stress pair. That is asserted as a test, because it is the reason the
assumption looked safe for twenty years of published cases.

**Remediation (diagnostic, not a fix).** A fifth axis is the fix and is a new
capability, which the brief forbids. Instead, `ClosureTable.build()` now
**measures the angle sensitivity for its own fluid pair** (18 extra
augmented-Lagrangian solves against 4 743 for the table itself) and
`assumption_report()` states it with a verdict — "benign" below 5 %, otherwise
"LARGE — the θ = 0 table is not valid for this pair". `scripts/kgep1_run.py`
prints it unconditionally, the same pattern as the A-3 range record. A table
that is not valid for its pair now says so on every run instead of being
silently wrong.

**K-GEP-1's reported numbers are not affected** — the 2.6 % figure was measured
directly on that pair — but the table design is not generally valid, and that is
now recorded where it cannot be missed.

**Tests.** `test_b1_angle_assumption_is_bounded` (measures the O(1) deviations
and the non-monotonicity), `test_b1_build_measures_and_reports_the_angle_sensitivity`
(the verdict must say "not valid for this pair"), and
`test_b1_newtonian_pair_has_no_angle_sensitivity` (the reason it was missed).
All three fail against the old code.

---

## Q-2 — BF25 (2.23) prints two forms that disagree — **CONFIRMED; a published error**

BF25 (2.23) gives `𝓘₃` twice, joined by "=". In the unit channel they differ by

```
first − second = −c̄³√m (c̄−1)² (c̄+2) / (6 (m c̄³ − c̄³ + 1))
```

which is strictly negative on `0 < c̄ < 1`. One of the two is wrong.

**The first is correct**, established without using (2.23) at all: integrating
the velocity profile of BF25 (2.9)–(2.10) directly reproduces (2.13), (2.24),
(2.25) and (2.26) exactly and yields the first form. `docs/assumptions.md`
CONV-03 recorded that **(2.27)** is misprinted; it did not record that **(2.23)
itself is internally inconsistent**, which is a second, separate defect in the
same equation and a trap for anyone re-deriving from the paper.

**Test.** `test_q2_bf25_2_23_two_printed_forms_disagree` asserts the gap
symbolically, asserts the velocity-profile derivation selects the first form,
and asserts the code implements it.

---

## Q-1 — the table's first `c̄`-interval misrepresents the flux near the tip — **REPRODUCED, bounded, not fixed**

Confirmed: at `H = 1` on the K-GEP-1 pair the table's linear interpolation over
`c̄ ∈ [0, 1/30]` gives a constant `f/c̄ = −31.4` where a direct solve gives
`−0.74` at `c̄ = 0.002` rising to `−21.6` at `c̄ = 0.022`.

**It does not affect any reported number.** The leading wave speed
`max_c f(c̄)/c̄` is attained at `c̄ = 1` and equals **1.00000 from both the table
and a direct solve**, so the "single shock at the mean speed" reading of
`output/kgep1_results.md` stands. Away from the tip the interpolation is good:
the worst error in `f = q₀ + b𝓘₃` over all `c̄`-midpoints is **0.93–0.95 % of
the flux span** at `H ∈ {0.3, 1.0, 1.75}`.

**What is affected, and is now caveated rather than fixed:** the ZF23 dispersion
statistics `σ_{w+r}` and `|w̄_r^+|` are computed from the small-`c̄` tail, which is
exactly where the interpolant is worst. Both are far from their (3.6)
thresholds (0.0200 and 0.0031 against 0.08 and 0.05), so the *classification*
"not dispersive" is robust; the *numbers* are not established to better than the
tip error, which is not quantified. `NUM-03` (closure-table resolution) remains
**pending** and this is part of why.

---

## R-1 — `breakthrough_at` depended on the logging frequency — **found during remediation, not in the audit**

`RunResult.breakthrough_at` and `efficiency_at` interpolate the **recorded**
history, and `record_every` sub-sampled that history. So a reported physical
quantity depended on an output setting. Measured on ZF22 case 1 while re-running
it here: `t_br@0.01` = **0.0409** at `record_every = 200` against **0.0537** at
`record_every = 1` — a **31 %** difference on a headline number, from a logging
option.

The hazard was already known in the same function: `run`'s internal `t_br`
detector carries the comment *"taking the recorded one would make t_br depend on
the output frequency"* and interpolates between **steps**. The fix was applied
there and never to the accessor.

**Fix.** The scalar history (`times`, `effs`, `outlet`) is appended on every step
regardless of `record_every`; ~2 MB over a 10⁵-step run. `record_every` now
governs only the expensive part — the `StepReport` list and the `on_step`
callback, which sees the whole concentration field.

**No previously reported number moves.** `scripts/zf22_table3.py` and
`scripts/kgep1_run.py` both already ran with `record_every = 1`, and the K-GEP-1
checkpoint holds a full-length history (85 475 samples for 85 474 steps). What
was affected is any *diagnostic* run using `record_every > 1` — including the
case-1 mesh study above.

**Test.** `test_r1_breakthrough_does_not_depend_on_the_logging_frequency`
asserts bit-identical `t_br` and efficiency across `record_every ∈ {1, 7, 50, 0}`.
Verified to fail against the previous commit.

---

## B-5 — two records of the same run disagreed — **REPRODUCED, settled by re-running**

`docs/assumptions.md` recorded ZF22 case 1 as `η_E = 0.375`, `t_br = 0.059`;
`output/zf22_table3.md` recorded 0.366 and 0.054, both at 20 × 200, CFL 0.5.

**Re-run in this session, all ten cases:** case 1 gives `t_br@0.01 = 0.054`,
`η_E = 0.366`, and **every one of the ten reproduces the values already in
`output/zf22_table3.md` exactly**. So the results file is current and the
register rows were stale — written from a superseded run and never updated
after the NUM-26/NUM-29 fixes. Both register locations are corrected and now
point at the results file as the single record.

This also confirms that nothing in this remediation moved any ZF22 number:
`scripts/zf22_table3.py` already used `record_every = 1`, so R-1 did not touch
it, and A-1's boundary condition was deliberately left unchanged.

---

## R-2 — a checkpoint could be resumed under different physics — **found during the UI inventory**

Not in the audit. Surfaced while answering the UI build spec's §5 requirement
that "every run writes the exact config that produced it".

**The hazard.** `Simulation.run` resumed whenever the checkpoint path existed
and the grid **shape** matched:

```python
if checkpoint_path and os.path.exists(checkpoint_path):
    saved = z["c"]
    if saved.shape != c.shape:
        raise ValueError(...)
    c = saved                     # <- and nothing else was checked
```

The checkpoint filename carries six settings
(`wall`, `inflow`, `n_phi`, `n_xi`, `w0`, `volumes`); a run has about twenty.
**Fluid properties, eccentricity, CFL, casing OD and the table resolution are in
neither the filename nor the grid shape.** So editing `mud_density` in
`config.py` and relaunching the same command produced the same path and the same
shape, and the run **silently continued a concentration field computed under the
old physics**. The only trace was a line reading `resumed from ...`.

`--cfl` is the sharpest case, because it changes the answer and is not in the
filename at all.

**Fix, in two layers.**

1. `d2dga/runio.py` (new, ~200 lines, no numerics): writes the whole `Config`
   tree, the CLI arguments and the environment (git SHA, library versions) as
   `<checkpoint>.config.json`, and the post-processed results as
   `<checkpoint>.metrics.json`. `fingerprint()` hashes the physics-and-mesh
   subset — deliberately not the command line, paths or machine, so re-running
   the same physics elsewhere still resumes. `check_resume()` returns
   `(level, message)` with level `None` / `"warn"` / `"block"` and names the
   settings that moved, because "hashes differ" is not actionable.
2. `Simulation.run(..., checkpoint_tag=...)` stores the tag in the checkpoint and
   **refuses** to resume across a change. This is the backstop for a caller that
   forgets to pre-check.

A **legacy** checkpoint (written before tags existed) warns and proceeds rather
than blocking: refusing would break relaunching a run that finished earlier, and
all three existing K-GEP-1 checkpoints are complete.

**Verified end to end on the real runner.** Same command, `--cfl` changed from
0.5 to 0.25 — identical filename, identical grid:

```
REFUSING TO RESUME
  output/kgep1_ckpt_synthetic_no_axial_gradient_6x20_w0.2_v0.02.npz was written
  with DIFFERENT settings, so resuming it would continue a field computed under
  other physics:
    args.cfl: 0.5 -> 0.25
  Pass --no-resume to start over, or use a different output path.
```

and with the same CFL it resumes as before (`resumed from ...: step 969`).

**Also fixed while testing it:** `describe_mismatch` reported
`mean_velocities_m_s: [0.05, 0.2, 0.5] -> (0.05, 0.2, 0.5)` as a change, because
a tuple comes back from JSON as a list. The fingerprint never had that problem —
both serialise identically — but noise in the one message that has to be trusted
is worse than noise anywhere else. Both sides are normalised through JSON now.

**Tests.** `tests/test_runio.py`, 9 tests. The two that matter —
`test_simulation_refuses_to_resume_across_a_settings_change` and
`test_a_legacy_checkpoint_still_resumes_with_a_warning` — **fail against the
previous commit** with `KeyError: 'tag is not a file in the archive'`.

**No previously reported number moves.** The guard only refuses; it changes no
computation, and the three existing checkpoints are unaffected (legacy path).

---

## G-1 — a non-Newtonian mud, and a correction to B-1's cause — **REM-11**

**Requested by the user on 2026-09-19** (UI gap G-1, raised in `ui/INVENTORY.md`
and held until approved because it touches `config.py`). The mockup's setup
screen offers the mud a τ̂_Y, κ̂ and n; `FluidsConfig` could not express any of
them — `as_fluids()` wrote `n = 1` and `τ_Y = 0` as literals.

### What changed

`FluidsConfig` gains `mud_consistency`, `mud_power_law_index` and
`mud_yield_stress`. The field `mud_viscosity` is **gone**, replaced by a property
that returns κ̂ for a Newtonian mud and **raises** otherwise: a Herschel–Bulkley
fluid has no single viscosity, and returning κ̂ (Pa·s^n) under a name that means
Pa·s is the unit error this change exists to prevent. Nothing else read that
field. `from_record` migrates a pre-G-1 payload (`mud_viscosity` → consistency,
exact, since those runs were Newtonian by construction) and **raises on any key
it does not know**, so a future field cannot be silently dropped and the
reconstruction labelled reproduced.

Both fluids are now validated at construction: ρ̂ > 0, κ̂ > 0, n ∈ (0, 1],
τ̂_Y ≥ 0, every pumped velocity > 0. `n > 1` is *refused*, not extrapolated —
BF25 (2.8)–(2.10) and B02 (10) derive the gap-scale closures for shear-thinning
Herschel–Bulkley only, and nothing here is validated above 1.

**No previously reported number moves.** The default pair is bit-identical,
`repr` included — asserted by `test_m0_c1_default_fluids_are_bit_identical_to_the_pre_g1_code`,
which reproduces the old literal construction. `repr` matters because
`kgep1_run.make_table` hashes the scaled fluids into the closure table's cache
key; a changed label would have invalidated every cached table.

### The measurement that changed a claim

G-1 was held back because a yield-stress mud puts **both** fluids in the
yield-stress class, which is where B-1 is O(1). Measuring it produced something
else. Holding the mud **strictly Newtonian** (n = 1, τ̂_Y = 0) and moving only
its viscosity against the unchanged K-GEP-1 cement:

| κ̂₁ | m | B-1 worst |
|---|---|---|
| 1 mPa·s (shipped) | 0.0063 | **0.53%** |
| 2 mPa·s | 0.013 | 2.0% |
| 5 mPa·s | 0.032 | **11.4%** |
| 10 mPa·s | 0.063 | 35.7% |
| 20 mPa·s | 0.127 | **93.1%** |
| 50 mPa·s | 0.317 | 168.6% |

And a pair in which *both* fluids are Newtonian measures **exactly 0.00%** at
m = 0.0017 and m = 0.033.

`d2dga/gapscale/tables.py` explained the benign K-GEP-1 result by saying *"the
displaced fluid is Newtonian"*, and **`AUDIT_REPORT.md` B-1 repeats that
wording**. It is the wrong attribution. What makes the shipped pair benign is
that **m ≈ 0.006** — the mud is ~160× thinner, so the cement's yielded
structure, the only θ-sensitive part, barely feels the mud layer. The mud's own
rheology is a second-order lever on the same quantity.

Two consequences, and the second is the one that decided the user's fallback
instruction ("if the three fields don't work, go back to Newtonian-only with a
disclaimer on screen"):

1. The docstring is corrected in place, with the sweep printed in it. The audit
   report is **not** edited — that is forbidden — so the correction lives here
   and in `docs/assumptions.md` FLU-07.
2. **G-1 did not create this hazard and reverting G-1 would not close it.** A
   Newtonian mud at an ordinary 5 mPa·s — expressible before G-1 existed, by
   editing one field that was already there — is already at 11%. Removing the
   three fields would hide the hazard behind a field that remains. So the fields
   stay and the **disclaimer is keyed to the measured value**, not to whether the
   mud is Newtonian.

Read the figure as an **upper bound** on the angle error: it compares θ = 0
against θ = 90°, and a run whose front stays nearly axial never visits the worst
angle. It bounds what the stored slice can be wrong by; it does not say what
η_E is wrong by. That needs the fifth axis (NUM-14), not this probe.

### The second correction, and it is the sharper one

The sweep above uses a probe whose velocity axis stops at `|ū| = 10`. **The
sensitivity grows with `|ū|`,** and the production table's axis reaches **3000** —
it has to, because FLU-06 drives `|ū|` past 370. Rebuilt on the production axes
(31 × 5 × 17, umag to 3000, 141 s) the **shipped K-GEP-1 pair** measures

| `\|ū\|` | 0 | 8.3 | 3000 |
|---|---|---|---|
| B-1 | 0.0% | 0.4% | **9.8%** |

worst **9.8%**, which this method's own grading calls **SIGNIFICANT, not
benign**. It had never been measured on those axes before: the production runs
**loaded a cached table**, and a loaded table reports `not measured`. So the
audit's own diagnostic was present and silent on the very run it was written for.

**No computed number moves** — this is a diagnostic, and no gate changes. What
changes is how K-GEP-1 numbers may be quoted: `AUDIT_REPORT.md` B-1's "benign for
the K-GEP-1 pair" is not supported at the resolution the runs actually use. The
9.8% must travel with them, and it peaks at large `|ū|` — which is exactly where
the front is most tilted, i.e. where θ is furthest from 0. The two errors are
correlated, not independent.

A by-product worth recording: the same production-axis build on the
Herschel–Bulkley mud (1400 / 0.020 / 0.70 / 4.79) took **954 s and left 12 of
935 gap solves unconverged**, against 0 for the water pair. A non-Newtonian mud
is not only harder to justify, it is measurably harder to solve.

### The disclaimer

Both the CLI and the Results screen now state, before any number:

* **provenance** — `FluidsConfig.departures_from_validated()` names every fluid
  field that differs from the pair every published figure was computed with, as
  `name: was -> is`. A screen can therefore say *which* field moved, not merely
  that something did.
* **measurement** — the run's own `assumption_report()`, read from
  `metrics.json`. A run that recorded none gets a **warning, not a green tick**,
  and the number is not guessed. The two are kept separate because only the
  second licenses a result.

### Verification

`tests/test_m0_config.py`, 11 gates (M0-C1..C11): bit-identity, the refusing
`mud_viscosity` property, validation of both fluids over 8 × 2 bad inputs, the
derived rheology label, `from_record` round-trip and legacy migration, the
unknown-key refusal, the angle-exact Newtonian pair, the m-sweep locking the
corrected claim, and the named departures. M0-C10's thresholds are one-sided and
far from the measured values (benign < 5% vs 0.53%; not-benign > 50% vs 93%) so
the gate tests the claim rather than the solver's tolerance.

---

## UI screens 1, 2 and 4 — what building them found

The interface was built in the build spec's order. Screen 3 (Results) shipped
first because it could be verified against data that already existed; the other
three write, launch or judge, and each turned up something the code did not
show.

### Four defects, none of them visible in the source

**The Validation screen manufactured a clean bill of health.** Its first ZF22
parser filtered rows with `startswith("| ")` and then sliced `[2:]`.
`output/zf22_table3.md` writes its separator as `|---|---|`, with no space, so
the filter removed the separator and the slice removed **case 1** — the one
case of ten that does not reproduce. The screen displayed **9 / 9**, which is
exactly the error in the mockup I had corrected G-8 for, arrived at
independently by discarding the counter-evidence. Parsing moved into
`ui/reports.py` where it is testable without a browser; the separator is now
found by its content; UI-12/13/14 lock all three properties.

**A run with non-default fluids was invisible.** The `_f<hash>` suffix that
stops two fluid pairs writing to one checkpoint was missing from
`runs.discover`'s pattern, so a run launched from the interface completed,
wrote its config and metrics, and never appeared in the list. Nothing failed.
Found by trying to open a run the interface had just produced. UI-15/16.

**"Reproduce from the CLI" did not.** `Run.command()` listed six flags by hand
and omitted `--mud-*`, so a run on a non-default mud was shown a command that
reproduces a *different* run, under that heading. It is now built by the
launcher from the recorded payload, and UI-17 closes the loop by parsing the
displayed line back and requiring the same `Config`.

**Relaunching a finished run crashed.** `res.reports[-1].n` raises `IndexError`
when a resume finds the checkpoint already at `t_end`. Pressing Run twice is
free, so a launcher hits this immediately. `RunResult.steps` (REM-13) fixes it
and also fixes a quieter error: with `record_every > 1` the last *report* is
not the last *step*, so every printed step count was short.

### And one I nearly shipped

The Results screen threw `NameError: name 'json' is not defined` on the
"Recorded configuration" panel — a path that had never executed before,
because until today no run had a recorded payload. I had screenshotted that
page and read the first sixty lines of its text, which did not reach the
error. The fix is not the import; it is that every screen is now swept with
Playwright in several states and the captured text is searched for
`Traceback`, `NameError` and the rest. Reading the top of a page is not
looking at it.

### What was added to the solver, and why none of it is a UI calculation

`--status-json` and the `preparing` stub; `Simulation.conservation_error`;
`RunResult.steps`; `postprocess.narrow_side_efficiency` and
`tbr_relative_error`. REM-12..15. The rule the spec sets — "if a display needs
data the solver does not emit, add an emitting hook, not a calculation in the
UI layer" — is what forced each of them; the alternative for the volume error
in particular was a second copy of an invariant that could disagree with the
run's own report.

---

## The verification this stage rests on — BCF25 §IV — **REM-16**

**User direction, 2026-09-19:** at this stage validation is not a comparison
against CFD; it is verification against mass conservation and its percentage
error, as paper 11 (BCF25) does it.

That is a sharper instruction than it looks, because BCF25 specify the
computation exactly, and this codebase was not doing it.

### What BCF25 do

> "for each fluid, we compute the total volume in the annulus in two ways.
> First, we compute Vol₁(k), by integrating cⁿ_{k,i+1/2,j+1/2} numerically over
> the interior volume of the annulus, i.e. we multiply by H rₐ Δφ Δξ and sum
> over the interior cells. Second, we compute Vol₂(k) at each timestep,
> starting with the initial volumes in the annulus and adding/subtracting the
> inflow/outflow of each fluid at the start and finish of the annulus, at each
> timestep. We normalize these quantities with the total annulus volume and
> subtracted to give a relative volumetric error."

Their Figs. 7 and 15 plot that against time, one curve per fluid, across five
cases and both the 2DGA and D2DGA models, and report it stays at **~10⁻¹⁵**.

### What this codebase was doing, and why it was not comparable

`Simulation.run` computed `max |Vol₁ − Vol₂| / |Vol₂|`. Two differences:

1. **The normalisation.** Dividing by the volume *present* rather than by the
   annulus volume. `Vol₂ → 0` at the start of a run, so that denominator
   inflates the error exactly where the absolute error is smallest. It is a
   *stricter* measure, not a wrong one — but it is not BCF25's, and a number
   computed one way cannot be set beside a number computed the other.
2. **One scalar, not a series.** A running maximum cannot be plotted as their
   Fig. 7, and cannot show whether the error grows, drifts or spikes.

Both are now produced. The old number is kept and labelled, because it is the
stricter statement and every earlier result in this repository quotes it.

### What the two-fluid version proves, and what it does not

This is the part worth being careful about, because the figure invites a wrong
reading.

BCF25 evolve **K = 3** concentrations independently, so their three curves are
three separate ledgers and "the concentrations sum to 1" is a genuine result of
their scheme. Here **K = 2** and only `c̄₂` is evolved: `c̄₁ ≡ 1 − c̄₂`
pointwise. The two-fluid closure then makes fluid 1's face flux *identically*
the total volumetric face flux minus fluid 2's — `q_{1,0} = 1 − q₀`,
`𝓘_{1,3} = −𝓘₃`, and the LLF dissipation term `−½a(c_N − c_S)` changes sign
with `c` and cancels exactly. So

    err₁ + err₂ = −(total volume in − total volume out) / V

identically, and fluid 1's curve is **not independent evidence**. The only part
that is independent is the total-flux imbalance, which asks whether the
elliptic solve delivered the same `Q` at both ends. It is reported as its own
curve rather than being folded into a second fluid curve that would look like
corroboration. The identity is asserted as a gate (M8-T4) to a residual of
~9 ε, which is cancellation in forming `V − Vol₁`, not slack.

The total volumetric flux is taken from the stream function alone —
`0.5 (Ψ[-1,j] − Ψ[0,j])`, because the advective flux with `q₀ = 1` telescopes —
so it is closure-free and measures what the elliptic solve produced, not what
the transport scheme did with it. Checked against the face flux with the
annulus full of fluid 2, where the two must coincide: they agree to **0.0**
(M8-T3).

### Measured

| | |
|---|---|
| Run | steps `N` | worst \|err\| | √N·ε | ratio |
|---|---|---|---|---|
| 8 × 40 Newtonian test pair | 150 | 2.04 × 10⁻¹⁵ | 2.7 × 10⁻¹⁵ | **0.75** |
| 8 × 24 **K-GEP-1 fluids**, 0.35 volumes | 15 359 | 7.82 × 10⁻¹⁴ | 2.7 × 10⁻¹⁴ | **2.84** |
| 16 × 80 production (older measure) | 85 474 | 2.93 × 10⁻¹⁴ | 6.5 × 10⁻¹⁴ | **0.45** |

Total-flux imbalance: **exactly 0.0** in every case. BCF25 report ~10⁻¹⁵.

**The raw exponent is not the comparison, and quoting only the 2 × 10⁻¹⁵ from
the small case would be picking the flattering one.** Error accumulates with
the step count, so the scale-free reading is the ratio to √N·ε. On that
reading the K-GEP-1 pair sits at 2.8 — the highest of the three, and the
excess over a pure random walk is where FLU-06 shows up: with `b·𝓘₁ ≈ 1500`
each step's boundary flux is orders larger than the net volume change it
produces, so there is more cancellation per step than a random-walk estimate
assumes. It is still roundoff, not bookkeeping; a bookkeeping error would be
orders away, not a factor of three. The Validation screen prints the ratio
rather than the exponent for exactly this reason, and its pass band is
20 × √N·ε rather than a flat 10⁻¹⁴, which would mark a correct long run as
failing.

Read a longer run against roundoff, not against 10⁻¹⁵ flat: √N·ε is
6.5 × 10⁻¹⁴ at N = 85 474 steps, which is why the gate is set at 10⁻¹⁴ with a
second assertion tying it to √N·ε rather than to a bare constant.

The K-GEP-1 production run recorded **2.93 × 10⁻¹⁴** over 85 474 steps — but
under the *older* normalisation, because it predates this ledger. That is
**0.45 √N·ε**, i.e. below the random-walk roundoff scale for its length; and
since the older measure is an upper bound on the BCF25-normalised one (it
divides by a volume that is never larger than the annulus), the BCF25 number
for that run is at or below 2.93 × 10⁻¹⁴ too. It is arithmetic, not
bookkeeping. The run has **not** been repeated to record the ledger properly —
it is ~3.3 h — so the Validation screen shows that run as having no ledger
rather than presenting the old number as if it were BCF25's.

### Where it appears

Printed by `scripts/kgep1_run.py`, stored as a series in the checkpoint (and
restored across a resume — M8-T5, because a ledger that restarts at a
checkpoint would show a clean run that was not clean), summarised in
`metrics.json`, shown live on the run monitor, and plotted on the Validation
screen as BCF25's Figs. 7 and 15 with their 10⁻¹⁵ line drawn on the axes. That
screen now opens by saying, in those words, that this stage does not validate
against CFD and that verification and validation are different questions.

---

## A-3b — my own A-3 fix broke the elliptic solve — **REM-17**

**Found by re-running the production case**, which is the only reason it was
found at all. Every gate passed throughout.

### What A-3 did, and what it actually changed

A-3 (remediation, 2026-09-18) anchored the closure table's velocity axis at
zero, because `|ū|` reaches ~6 × 10⁻³ at stagnation points and queries below
the old axis minimum of 0.02 were being served by linear extrapolation. It
anchored with **one** node — `[0, 0.02]` — and the comment written at the time
said it "changes nothing: the gap solve at umag = 0, 1e-3 and 0.02 agrees to
six significant figures".

That agreement is exactly the problem. Because those values agree, the
`[0, 0.02]` segment of the interpolant is nearly **flat**, while the next
segment `[0.02, 0.1]` is not. The kink between them breaks the elliptic
**Picard** iteration, which updates `|ū|` from the previous `Ψ` and
re-evaluates the mobility: it oscillates across the kink instead of
contracting.

### Measured, same field, only the axis changed

16 × 80, K-GEP-1 fluids, the concentration field at step 1906 of the re-run:

| axis | nodes | Picard | residual | under-range queries |
|---|---|---|---|---|
| `[0.02] + geom` (pre-A-3) | 16 | **6** | 1.54e-09 | **36 per solve** |
| `[0, 0.02] + geom` (A-3) | 17 | **100, at the cap** | 3.37e-05 | none |
| `[0, .005, .01, .02] + geom` | 19 | **6** | 1.54e-09 | **none** |

The third fixes both problems. Along a fresh 2982-step trajectory the fraction
of steps hitting the cap is **37.5 % → 1.5 %**; on the live production run,
**37.5 % → 2.8 %**, and the wall-clock estimate fell from **8.9 h to 2.7 h**,
because the solver had been spending 100 iterations per step.

**It is not 0 %.** The pre-A-3 production run reported "0 steps hit the
iteration cap, worst residual 0.00e+00" over 85 474 steps. The refined axis
leaves 1.5–2.8 %. I am not claiming the axis is now optimal — only that it is
13–25× better and removes the under-range queries A-3 existed to remove.

### The second defect: R-2's guard could not see this

The axis is hard-coded in `kgep1_run`; it is in neither `Config` nor the CLI
arguments, so the **settings fingerprint did not cover it**. A checkpoint
computed with the broken axis would have resumed silently under the fixed one —
the precise failure R-2 was built to prevent, reached by a route R-2 did not
cover. *The guard only ever sees what it is handed.* The runner now hands it
the axis (`closure_umag_axis` in `_PHYSICAL_ARG_KEYS`), and the step-1906
checkpoint is archived rather than resumed, because it is now refused by tag
rather than by anyone remembering to refuse it.

### Why there is no direct gate, said plainly

Two attempts at a cheap unit test that reproduces the Picard failure both
**failed to discriminate**, and they are recorded in `tests/test_m9_picard_axis.py`
rather than quietly dropped:

* a reduced table (11 c-nodes, 3 H-nodes) gives 52 iterations on **both** axes —
  too coarse everywhere for this kink to dominate;
* a synthetic flat front at 16 × 80 gives 100 on **both** axes — harsher than
  the real flow, so it cannot tell them apart.

Reproducing it needs a production-resolution table (~10 min to build) and a
field from a real trajectory. That is not a unit test. M9-T1..T4 therefore gate
the axis itself, as a regression lock on the measured boundary, and its presence
in the fingerprint, which is a genuine functional property that was genuinely
absent. The failure is caught at **run** level instead: the cap count existed
already and nothing read it, so it is now stated as a fraction with a verdict,
recorded in `metrics.json`, and shown on the monitor and the Results screen.

A gate that cannot fail would have been worse than admitting the gap.

### What this says about the audit

A-3 was my own remediation, written with a claim ("changes nothing") that I did
not test. It passed every gate for five days. The thing that caught it was
running the case the numbers are quoted from — which is what the user asked for,
and which the audit had recorded as outstanding without acting on it.

