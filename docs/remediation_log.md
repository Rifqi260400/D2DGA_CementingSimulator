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

| `n_phi` | `n_xi` | `t_br@0.01` | `η_E` |
|---|---|---|---|
| 20 | 200 | 0.0409 | 0.3657 |
| 40 | 200 | 0.0456 | **0.3939** |
| 80 | 200 | *(see `gate_status.md`)* | |
| 20 | 400 | | |
| 20 | 800 | | |

`η_E` moves **+7.7 % on azimuthal refinement alone**, at fixed `n_xi`.

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

## B-1 — the closure table assumes `ū ∥ G̃_b` — **REPRODUCED, prose corrected, bound measured**

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

The production run reached `\|ū\| = 75`, so ~3 % on `𝓘₃` is the operating error.
`NUM-14` already refuses `β ≠ 0`, where a second angle appears; this is about
`β = 0`, where the approximation is live and was undocumented.

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
