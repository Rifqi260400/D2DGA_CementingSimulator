# Independent audit — D2DGA cementing displacement solver

Scope: the whole repository as of `c47c42f`, audited against the six source
papers in `docs/*.pdf` (text extractions in the session scratchpad), not
against the build spec and not against the builder's own record.

**No file other than this one has been modified. Nothing has been fixed.**

Measurement scripts referenced below as `scratchpad/audit/*.py` were written to
the session scratchpad, **not** to the repository, so that the no-modification
constraint held. They are not committed and will not survive this container; each
finding therefore quotes its numbers in full rather than relying on them.

---

## 0. Disclosure — this audit is not independent

I am the same model that built this solver. I do not have a second, uncontaminated
reading of these papers: where I "confirm" an earlier decision I may simply be
repeating the reasoning that produced it. Read this report with that weighting:

* **Confirmations of the builder's own conclusions are weak evidence.** Where I
  could, I re-derived by a *different route* than the one recorded in
  `docs/derivation.md` (for example, integrating the velocity profile of
  BF25 (2.9)–(2.10) directly rather than manipulating (2.23)) so that a shared
  algebraic slip would show up as a disagreement. Those are flagged.
* **New findings are worth more**, because producing them required contradicting
  the builder. Four of the findings below (A-1, A-2, A-3, B-1) are things the
  build process did not record.
* **One of my own initial suspicions turned out to be wrong** and is written up
  as such in §4, rather than deleted, because a reader deserves to know which
  way the errors ran.

Phase A (blind re-derivation and measurement) was completed before
`docs/assumptions.md` or `docs/derivation.md` were opened. Phase B compared.

---

## 1. Findings

Severity: **CRITICAL** = a reported number is wrong; **MAJOR** = a reported
number is materially uncertain or unsupported; **MINOR** = real but small or
documentation-only; **QUESTION** = unresolved, needs a decision or a measurement.

### There are no CRITICAL findings.

I looked specifically for a wrong sign, a wrong scaling conversion, a broken
conservation law, a violated maximum principle and a mis-stated closure, and
found none. The five conventions most likely to be wrong (§4) are all right,
two of them by a derivation route the builder did not use. I am not
manufacturing a headline finding in their place.

---

### A-1 — MAJOR — The outflow boundary imports displacing fluid. It is the unfixed mirror of NUM-29.

**Location** — `d2dga/transport.py:265-271` (`_pad_xi`), `:368-376` (outlet
column replication), `:399` (`Xi = adv_x + buoy_x - ...`). Registered as NUM-20
in `docs/assumptions.md`.

**What the papers say** — BF25 (2.20)–(2.21) define the fluid-2 flux as
`q = r_a(∫_0^{y_i} v dy, ∫_0^{y_i} w dy)`, whose axial part is
`r_a H w̄ q₀ + b̃ r_a H³𝓘₃` with `𝓘₃ < 0`. None of the six papers states an
outflow condition; BCF25 §III does not either.

**What the code does** — `ξ = Z` is closed by replicating the last interior
column into a ghost cell. The LLF dissipation on that face then vanishes
identically (equal states) and the face carries the full model flux, buoyant
term included. Because `𝓘₃ < 0` and `b > 0`, that flux is **negative** over most
of `0 < c̄ < 1`: the face runs backwards and fluid 2 enters the annulus through
the outlet, out of a ghost region that does not exist.

**Why it matters** — `η_E` is the mass in the annulus divided by its capacity, so
every unit that enters through the top is added directly to the headline number.
It also delays `t_br` measured at the outlet. The existing volume-balance
diagnostic in `scripts/kgep1_run.py` cannot see it: after breakthrough it lumps
all outlet traffic into one signed number.

**Evidence** — measured, `scratchpad/audit/a5_outlet.py`, 10 × 50 mesh, CFL 0.5:

| case | `b_BF25` | 1-D flux `f = q₀ + b𝓘₃` | most negative `Σ_φ Ξ[:,-1]` | fluid 2 entering through the outlet |
|---|---|---|---|---|
| ZF22 case 4 | +447.2 | spans [−5.996, +1.000] | **−10.07** (pumped rate = 1) at `t/Z = 0.911` | **3.26 volume-units = 3.5 % of the 1.2-volume job** |
| ZF22 case 1 | −22.36 | spans [0, 1] — never negative | 0 | 0 |

For the K-GEP-1 production configuration the same quantity is bounded at
**−4.61 × Q** (`scratchpad/audit/a5_outlet_kgep1.py`, worst uniform-`c̄` outlet
state, `c̄ = 0.36`); an order-of-magnitude volume estimate over one front
crossing of the outlet cell is ~5.6 % of an annulus volume. That was **not**
measured over a real run — see §6.

Note the structure of the error: it is zero whenever `f(c̄) ≥ 0` for all `c̄`
(case 1, and every weakly buoyant case), and grows with `|b|H³|𝓘₃|`. That is
the same scaling as NUM-29, which was found and fixed at the inlet on exactly
this reasoning ("at a Dirichlet inflow there is no second cell and it becomes a
pure SOURCE"). The inlet received a physically argued treatment; the outlet did
not.

**Confidence** — high that the mechanism and the measured magnitudes are real.
Medium on how much of the 3.5 % survives to bias `η_E`, because some of the
imported fluid leaves again later; I did not separate the two.

---

### A-2 — MAJOR — `t_br` is mesh-limited on the meshes actually used, and the K-GEP-1 result has no mesh-refinement evidence at all.

**Location** — `d2dga/simulation.py:189` (`breakthrough_threshold`, default 0.01),
`scripts/kgep1_run.py` (`--n-xi` default 150; the two completed runs used 80),
`output/kgep1_results.md`, `output/zf22_table3.md`.

**What the papers say** — ZF22 §5 defines `t_br` as when the displacing fluid
"first exits the annulus". BF25 §3.1–3.2 make that ambiguous by construction:
`q₀'(0) = 1.5`, so the D2DGA tip is asymptotically thin and any positive
threshold is a choice. The code records this honestly as NUM-21 and reports
three thresholds.

**What the code does** — the scheme is first-order LLF. On the exact reduction
where the answer is known in closed form, the front at the quoted 0.01
threshold runs systematically ahead of the true characteristic.

**Why it matters** — `t_br @0.01` is the leading number in both result files.
The K-GEP-1 headline "breakthrough at 0.99 of the piston time" rests on it.

**Evidence** — concentric annulus, `m = 1`, `b = 0`. The exact entropy solution
is a pure rarefaction, `c̄(ξ,t) = √(1 − (ξ/t)/1.5)` for `ξ/t ≤ 1.5`, derived here
from `f(c̄) = 1.5c̄ − 0.5c̄³` (concave, `f'` decreasing) — no repository code was
consulted. `scratchpad/audit/a3_r2.py`:

| `n_xi` | L1 error | front speed at `c̄=0.01` (exact 1.4998) | at `c̄=0.1` (1.4850) | at `c̄=0.5` (1.1250) |
|---|---|---|---|---|
| **80** (K-GEP-1 production) | 2.62e-2 | 1.7625 → **+17.5 % fast** | 1.5875 → +6.9 % | 1.0875 → −3.3 % |
| **150** (script default) | 1.71e-2 | 1.6867 → **+12.5 %** | 1.5400 → +3.7 % | 1.1000 → −2.2 % |
| **200** (ZF22 table) | 1.38e-2 | 1.6650 → **+11.0 %** | 1.5350 → +3.4 % | 1.1050 → −1.8 % |
| 400 | 8.42e-3 | 1.6125 → +7.5 % | 1.5125 → +1.9 % | 1.1125 → −1.1 % |

L1 converges at ~0.7 order; the 0.01-threshold front speed converges far more
slowly than that and is still 7.5 % fast at `n_xi = 400`. Azimuthal spread is
machine zero (4e-15), and `c̄ ∈ [0, 1]` throughout, so this is purely the
first-order smearing and not a coupling artefact.

Separately: **there is no completed mesh-refinement run of K-GEP-1.** Both
quoted runs are 16 × 80 (`output/kgep1_bc_b02.log`, `kgep1_bc_uniform.log`). The
only finer attempt, 16 × 120, reached `t/Z = 0.025` before the container was
reclaimed (`output/kgep1_synthetic.log`). The cell aspect ratio at 16 × 80 is
`dξ/dφ = 7.05/0.0625 = 113`; `docs/assumptions.md` NUM-07 flags ~28 as a
concern and is still marked **Pending**. `M6-T3` exists as a test but runs a
Newtonian pair at 12 × 75 … 48 × 300 and deliberately measures **volume**, not
front position — "so it is not contaminated by how much the first-order scheme
smears the front". That is a reasonable design for a conservation gate and it
does not cover the quantity at issue here.

**Confidence** — high. The exact solution is closed-form and the trend is
monotone over four meshes.

---

### A-3 — MAJOR — The closure-table range guard was suppressed in the production runs, and those runs did go out of range.

**Location** — `d2dga/gapscale/tables.py:269-283` (`_check_range`), `:257`
(`__call__(..., warn_out_of_range=True)`); `ClosureTableRangeWarning` subclasses
`UserWarning`.

**What the papers say** — nothing; this is an implementation guard.

**What the code does** — `_check_range` emits a `ClosureTableRangeWarning` and
then **extrapolates anyway** (`RegularGridInterpolator(..., bounds_error=False,
fill_value=None)`). The warning is therefore the only signal that a result rests
on extrapolated closures. `docs/assumptions.md` FLU-06 explicitly credits it:
*"the closure table's M3-T3 range guard is what caught it, by name, rather than
silently extrapolating."*

The two completed K-GEP-1 runs' logs (`output/kgep1_bc_b02.log`,
`output/kgep1_bc_uniform.log`, 227 lines between them) contain **zero** warning
lines. That is proof of suppression, not of cleanliness: under Python's *default*
filter every occurrence would have printed, because the message text differs each
time (it names the cell and the value) so `default` never de-duplicates it. At
the rate measured below that is of order 10^4 warning lines that are absent. The
runs were launched under `python -W ignore::UserWarning`, which is exactly what
produces this; the launch command is not recorded anywhere in the repository, so
treat the mechanism as inferred and the suppression itself as established.

**Evidence** — re-running the identical configuration with warnings enabled
(`scratchpad/audit/a5_range.py`, 400 of the run's 85 474 steps):

```
axis c     requested [0, 1]                 table [0, 1]
axis H     requested [0.2619, 1.7596]       table [0.2370, 1.9340]
axis umag  requested [0.00108, 75.296]      table [0.02, 3000]   UNDER by 0.0189
axis gb    requested [27.9021, 27.9021]     table [27.9021, 27.9021]

ClosureTableRangeWarnings raised: 76
```

So the runs did query out of range, on the **low** end of the velocity axis —
the opposite end from the one FLU-06 discusses and widened the axis for.

**Measured impact is small.** Comparing the extrapolated table against a direct
gap-scale solve at `umag` down to 0.001 (`scratchpad/audit/a5_extrap_error.py`):
`q₀` 0.00 %, `𝓘₁` ≤ 0.01 %, `𝓘₃` ≤ 0.10 % at `c̄ ∈ {0.1, 0.5, 0.9}`. The closures
are flat in `|ū|` below 0.02 for this fluid pair, so the extrapolation happens to
be almost exact.

**Why it matters anyway** — the numbers are fine; the *control* is not. A run
that silently ignores its own correctness guard cannot be quoted as having
passed it, and the register does quote it. The next fluid pair or the caliper
wall may not be as forgiving.

**Confidence** — high on the suppression and on the 76 warnings; high on the
≤0.1 % impact for this fluid pair specifically.

---

### A-4 — MAJOR — The validation evidence is thinner than "nine of ten within 2.5 %" reads, and the entire Herschel–Bulkley path is validated by no published case at all.

**Location** — `output/zf22_table3.md`, `docs/assumptions.md` M6-T1 section,
`tests/benchmarks/zf22_cases.py`.

**What the papers say** — ZF22 Table 3 (verified against the paper, p. 27) lists
ten cases. Five of them (2, 4, 5, 9, 10) have `b_ZF22 ≥ 100`; ZF22's own text
says of them *"the effect of eccentricity is not very apparent … since the
buoyancy number is large enough to compete against the geometry influence"*, and
their published `η_E` all sit in 0.95–1.00.

**What the code does** — reproduces `η_E` to within 2.5 % on nine cases and
fails case 1 by −45 % (0.366 against 0.66; `t_br` 0.054 against 0.44).

**Why it matters** — three distinct weaknesses, none of which is visible from the
headline:

1. **Saturation.** Half the suite sits where any model with a stabilising
   buoyancy term returns ~0.95–1.00. Agreement there does not discriminate. The
   genuinely discriminating cases are 3, 7, 8 (`b = 10`, dispersive) — and those
   *do* pass, including the `m`-ordering between 7 and 8 (0.897/0.857 against
   0.90/0.84), which is real evidence.
2. **The one case that tests adverse buoyancy in 2-D is the one that fails.**
   Case 1 is the only `b < 0` case in any of the six papers. The 1-D reduction
   agrees with ZF22 (leading-wave `t_br` 0.391 against 0.44), so the closures are
   not at fault; it is the 2-D coupling that diverges. `docs/assumptions.md`
   BENCH-09 records the failure and two refuted hypotheses, which is the right
   posture — but it means the **sign and magnitude of the buoyancy term in the
   2-D solver are validated only in the direction that makes displacement
   better**, never in the direction that makes it worse.
3. **Every ZF22 case is Newtonian.** The whole Herschel–Bulkley path — the
   augmented-Lagrangian gap solve, the `M3` closure table and its interpolation,
   the Picard elliptic iteration, the `NUM-13` mobility floor, and the
   parallel-`ū`/`G̃_b` assumption of B-1 below — is exercised by **no published
   benchmark in this repository**. BENCH-08 (ZF23 Table 3) is also Newtonian.
   The K-GEP-1 run is therefore the first and only use of that path, and it is
   the one being quoted.

**Evidence** — ZF22 Table 3 read from the paper; `output/zf22_table3.md`;
`tests/benchmarks/zf22_cases.py` `CASES` (all built with
`HerschelBulkleyFluid.newtonian`); `d2dga/elliptic.py:220` (`TabulatedClosures`
is reached only from `scripts/kgep1_run.py`).

**Confidence** — high. This is a coverage statement, not a defect claim.

---

### A-5 — MAJOR — ZF22 case 1 is **not** mesh-converged in `η_E`, contrary to `BENCH-09`, and it moves toward ZF22 on refinement.

**Location** — `docs/assumptions.md` BENCH-09; `output/zf22_table3.md` row 1.

**What the record says** — BENCH-09 states, as the first of five established
facts about the case-1 failure: *"(i) it is mesh-converged — `t_br` moves from
0.0534 to 0.0530 between 20×200 and 40×400 while `Δφ, Δξ` halve, so it is not a
resolution artefact"*. That claim is supported by `t_br` alone.

**What I measured** — I re-ran the case at both meshes, and at `b = 0` as a
control (`scratchpad/audit/a4_case1.py`):

| mesh | `b_BF25` | `t_br @0.01` | `t_br @0.1` | `η_E` | steps | conservation |
|---|---|---|---|---|---|---|
| 20 × 200 | 0 | 0.3687 | 0.4049 | 0.7999 | 896 | 4.7e-15 |
| 40 × 400 | 0 | 0.3796 | 0.4086 | 0.8056 | 1 811 | 2.3e-15 |
| 20 × 200 | −22.36 | 0.0537 | 0.0540 | **0.3657** | 8 531 | 1.1e-14 |
| 40 × 400 | −22.36 | 0.0530 | 0.0532 | **0.4124** | 18 801 | 3.9e-15 |

`t_br` is indeed stable (−1.3 %). **`η_E` is not: it moves +12.8 %, from 0.366 to
0.412**, and it moves *toward* ZF22's 0.66. The `b = 0` control moves only
+0.7 %, so this is specific to the adverse-buoyancy case, not general mesh error.

**Why it matters** — `η_E` is the quantity the failure is stated in (−45 %), and
it is the quantity `docs/assumptions.md` itself calls "the cleanest comparison —
it is threshold-free, unlike `t_br`". Establishing mesh convergence on `t_br` and
then reporting the failure in `η_E` does not close the loop. At +12.8 % per mesh
halving and no sign of a plateau, the mesh-refinement hypothesis for case 1 is
**not** excluded — it is the one hypothesis BENCH-09 believes it has excluded,
and on this evidence it has not. Two further halvings at this rate would not
reach 0.66, so refinement is unlikely to be the whole story, but it is clearly
part of it and the current record says it is none of it.

**Also new, and it narrows the search** — the same run isolates where the
discrepancy is *not*:

* At `b = 0`, `e = 0.8`, the solver gives `t_br = 0.369–0.380` against the
  analytic wide-side channel prediction `1/(1.5 × 1.6531) = 0.4033` derived here
  from `w̄(φ) = H²⟨H⟩/⟨H³⟩` and `q₀'(0) = 1.5`. The ~7 % shortfall is exactly the
  A-2 tip bias at these meshes. **Eccentric channelling is right.**
* At `b = −22.36`, 1-D theory times the same channel factor predicts
  `1/(2.557 × 1.6531) = 0.2366`. The solver gives **0.053 — 4.5× faster**.

So the failure is neither in the closures (1-D theory agrees with ZF22), nor in
the elliptic solver (§2, machine precision), nor in the eccentric channelling
(above), nor in the bottom-hole condition (BENCH-09, refuted), nor in inlet
buoyant exchange (BENCH-09, refuted). It is specifically in **how the 2-D
solution responds to an adverse buoyancy source**, and it is worth 4.5× in the
front speed.

**Confidence** — high on the numbers (two meshes, conservation 1e-14, `b = 0`
control). Medium on the interpretation, because I did not run a third mesh.

---

### B-1 — MINOR (structural) — The closure table assumes `ū` is parallel to `G̃_b`; it never is, once the front tilts.

**Location** — `d2dga/gapscale/tables.py:167`
(`solve_fixed_mean_velocity(c, [0.0, umag], [0.0, gb*H])`),
`d2dga/elliptic.py:245-262` (`_umag`, `q0_I3`).

**What the papers say** — BF25 (2.8)–(2.10) make the gap-scale stress a
**2-vector**: `τ₂ = [−G + (1−c̄)G_b] y`, `τ₁ = τ_i − [G + c̄G_b](y − y_i)`, and for
a Herschel–Bulkley fluid `η_k = η_k(|τ_k|)`. `|τ|` depends on the angle between
`G` and `G_b`, so the four closures do too. BF25 figure 4 is precisely an
illustration of the direction change.

**What the code does** — tabulates on `(c̄, H, |ū|, g_b)` and reconstructs the
gap-scale state with both `ū` and `G̃_b` purely axial. `TabulatedClosures.__init__`
documents the angle problem for `β ≠ 0` and refuses that case
(`NotImplementedError`, NUM-14) — but the same angle appears at `β = 0` whenever
`v̄ ≠ 0`, which is exactly when the displacement is interesting.

**Why it matters** — it is an uncontrolled approximation on the only path the
production result uses, and it is invisible to every Newtonian benchmark
(constant `η`, no angle dependence).

**Evidence** — measured directly with the AL solver, K-GEP-1 pair, `g_b = 27.9`,
`H = 1` (`scratchpad/audit/a2_angle.py`); `θ` is the angle of `ū` from the axial
direction, `θ = 0` is what the table stores:

| `c̄` | `\|ū\|` | dev. of `𝓘₁` at θ=90° | `𝓘₂` | `q₀` | `𝓘₃` |
|---|---|---|---|---|---|
| 0.25 | 50 | +0.0 % | +0.0 % | +0.0 % | −0.1 % |
| 0.50 | 50 | +0.1 % | +0.3 % | +0.0 % | −0.7 % |
| 0.75 | 50 | +0.6 % | +1.4 % | +0.1 % | **−2.6 %** |

The error grows with `|ū|/(g_b·𝓘₁)`. The K-GEP-1 run reached `|ū| = 75`, so
~2.6 % on `𝓘₃` is the right order — small. **But the table's `umag` axis runs to
3000**, and I did not measure the deviation at the top of that axis, so the
approximation is bounded only over the range the one completed run happened to
visit.

**Confidence** — high that the approximation exists and is undocumented for
`β = 0`; medium on its size, because it is bounded only where measured.

---

### B-2 — MINOR — `derivative_bounds` has no range check at all.

**Location** — `d2dga/gapscale/tables.py:245-252`. Compare `__call__` at `:257`,
which does check.

Every LLF wavespeed on the tabulated path goes through `derivative_bounds`
(`d2dga/elliptic.py:264-270`), at 29 probe nodes per face. It builds its own
`RegularGridInterpolator` with `fill_value=None` and never calls `_check_range`,
so out-of-range queries there are silent **by construction**, not merely by
command-line flag. Given A-3's measurement that the production run does go
out of range on `umag`, this path was extrapolating silently in every one of
those 85 474 steps. The consequence is a wrong LLF coefficient, which affects
`Δt` and the amount of artificial dissipation — not conservation, and not the
maximum principle (the bound can only be too large or too small; too small is
the dangerous direction and is not detectable here).

**Confidence** — high on the code fact; **I did not measure** whether the
extrapolated bound ever fell below the true slope.

---

### B-3 — MINOR — Two places in the repository state the opposite of NUM-26.

**Location** — `d2dga/transport.py:154-158` (the "Monotonicity" docstring
section) and `docs/assumptions.md` NUM-17, "Sensitivity tested?" column.

Both say, of the LLF wavespeeds:

> "endpoint evaluation is sufficient — no interval maximum is needed"

NUM-26, forty lines earlier in the same docstring and two rows later in the same
table, says it is **not** sufficient, and reports a measured maximum-principle
violation (`c̄ = −0.0089` within two steps at `β = π/4`, `b = 20`, `e = 0.5`).
The implementation follows NUM-26 (`wavespeed="interval"` is the default). The
stale text is a live hazard: it is exactly the justification a maintainer would
need to switch the default back to `"endpoints"`, which is retained in the code.

**Confidence** — certain; both texts are quoted above.

---

### B-4 — MINOR — The assumption register has duplicate IDs whose contents contradict each other.

**Location** — `docs/assumptions.md`.

* **`CONV-03` appears twice.** The row under *Convention conversions* gives the
  master form with a **leading minus** and says BF25 (2.27) is "wrong three
  times over", sign included; the conversions read `BCF25 (23) = −MASTER·√m`.
  The row under *Scaling and gap-scale solver* gives the master **without** the
  minus, says (2.27) is "wrong twice", and reads `BCF25 (23) = MASTER·√m`.
  The code implements the **negative** form (`d2dga/gapscale/newtonian.py`), and
  §4 below confirms that is correct — so the second row is stale and, taken at
  face value, would reintroduce the sign bug that the first row describes as the
  most consequential error found in the whole build.
* **`NUM-02` appears twice**, with `ρ = r` probed per fluid pair and `tol = 1e-9`
  in one row and `r = ρ = 1`, `tol = 1e-10`, `N_y = 400` in the other.

The register's own header says "One row per modelling choice".

**Confidence** — certain.

---

### B-5 — MINOR — Two records of the same run disagree.

`docs/assumptions.md` (M6-T1 table and BENCH-09) records ZF22 case 1 as
`η_E = 0.375`, `t_br = 0.059`. `output/zf22_table3.md` records `η_E = 0.366`,
`t_br @0.01 = 0.054`. Both state mesh 20 × 200, CFL 0.5. One of them is from a
superseded run and was not updated. The discrepancy is small but it is in the
one case that is under active investigation, where a 2.5 % shift in `η_E` is
the size of the effects being chased.

---

### B-6 — MINOR — `Simulation.run(record_every=0)` raises `ZeroDivisionError`.

`d2dga/simulation.py:290`, `if n % record_every == 0 or t >= t_end:`. Zero is the
natural way to ask for "no history", and the docstring does not exclude it.
Found by hitting it (`scratchpad/audit/a3_reductions.py`). Cosmetic; no result
is affected.

---

### Q-1 — QUESTION — The closure table's first `c̄`-interval badly misrepresents the flux near the tip.

At `H = 1` with the K-GEP-1 pair, the table's linear interpolation on
`c̄ ∈ [0, 1/30]` gives a *constant* `f(c̄)/c̄ = −31.41`, where a direct gap-scale
solve gives `−0.74` at `c̄ = 0.002` rising to `−21.6` at `c̄ = 0.022`
(`scratchpad/audit/a6_tipspeed.py`).

This does **not** affect K-GEP-1's headline result: for that flux function the
leading wave speed `max_c f(c̄)/c̄` is attained at `c̄ = 1` and equals **1.00000
from both the table and the direct solve**, so the "single shock at the mean
speed, `t_br → 1`" claim in `output/kgep1_results.md` is independently confirmed.
Away from the tip the interpolation is good: the worst error in `f = q₀ + b𝓘₃`
over all `c̄`-midpoints is **0.93–0.95 % of the flux span** at `H ∈ {0.3, 1.0,
1.75}` (`scratchpad/audit/a6_interp.py`).

But every small-`c̄` statistic is computed in exactly the region where the
interpolant is worst, and `output/kgep1_results.md` quotes two of them —
`σ_{w+r} = 0.0200` and `|w̄_r^+| = 0.0031`, used for the ZF23 (3.6) "not
dispersive" classification (thresholds 0.08 and 0.05, verified against ZF23
p. 1277). Those are far from their thresholds, so the *classification* is
robust; the *numbers* are not established to better than the tip error, which I
did not quantify. `NUM-03` (closure-table resolution) is still marked **pending**
in the register.

---

### Q-2 — QUESTION — BF25 (2.23)'s two printed forms disagree with each other, and nothing records it.

BF25 prints `𝓘₃` twice in (2.23), joined by "=". Integrating both for the
Newtonian pair:

```
first form  - second form  =  c³√m(−c³ + 3c − 2) / (6(mc³ − c³ + 1))
                           =  −c³√m(c−1)²(c+2) / (6(…))   ≤ 0,  not identically 0
```

The code implements the **first** form verbatim (`d2dga/gapscale/closures.py:86`)
and §4 below shows by an independent route that the first form is the correct
one. So the choice made is right. But `docs/assumptions.md` CONV-03 records only
that **(2.27)** is misprinted; it does not record that **(2.23) itself is
internally inconsistent**, which is a second, separate defect in the same source
equation and a trap for anyone re-deriving from the paper.

---

## 2. What I re-derived and found CORRECT

These are confirmations. Per §0 they are weak evidence where I reused the
builder's reasoning, and stronger where I did not. The route is stated for each.

| Item | Verdict | Route |
|---|---|---|
| **CONV-01** — `∇_aΨ = (2Hw̄, −2Hv̄)` | correct | Read verbatim at BF25 line 585 of the extraction. Independent of the builder. |
| **CONV-02** — `b_BF25 = √m · b_ZF22` | correct | ZF22 §2.3 (`τ̂₀ = μ̂₁ŵ₀/d̂`, `m = μ̂₁/μ̂₂`, `b = (ρ̂₂−ρ̂₁)ĝd̂²/(μ̂₁ŵ₀)`) against BF25 (2.4) (`Fr* = √(τ̂₀/ρ̂₁ĝd̂*)`) and (2.12). Also checked that ZF22's `F` and BF25's `Fr*` use the same length (`δ₀r̂*_a = d̂`), which was the factor most likely to be lost. |
| **CONV-03** — `𝓘₃ < 0`, and the exact Newtonian form | **correct, by a route the builder did not use** | Integrated the velocity profile directly: BF25 (2.9)–(2.10) → `du/dy = τ/η` → no-slip at `y=H`, continuity at `y=y_i` → `∫_0^{y_i} u dy`, then matched coefficients against `ūq₀ + 𝓘₃g_b`. This reproduced **(2.13) exactly, (2.24) exactly, (2.25) exactly, (2.26) exactly**, and gave `𝓘₃` equal to the code's `script_I3` **identically in sympy** — and to neither sign of (2.27) as printed. The builder derived from (2.23); I never used (2.23). |
| **CONV-04** — `b_s` carries **no** `r_a` | **correct, independently** | Eliminated `p` from BF25 (2.11)–(2.13) by `∂_ξ(∂_φp) = ∂_φ(∂_ξp)` in plain `(φ,ξ)` coordinates. Result: `∂_φ[a_φ∂_φΨ] + ∂_ξ[a_ξ∂_ξΨ] + ∂_φ[b_s cosβ] + ∂_ξ[b_s r_a sinβ sin πφ] = 0` with `a_φ = 1/(2r_aH³𝓘₁)`, `a_ξ = r_a/(2H³𝓘₁)`, `b_s = (1/Fr*²)[1 − Δρ(c̄ + 𝓘₂/𝓘₁)]` — **exactly** `d2dga/elliptic.py:11-14`. Both printed forms of BF25 (2.18) carry a `1/r_a` on the buoyancy part that this derivation says does not belong; they differ from each other as well. |
| **CONV-10** — one `r_a` survives in the transport buoyancy flux | **correct, independently** | Same velocity-profile derivation: `q_ξ = r_a H w̄ q₀ + b̃ r_a H³𝓘₃ cos β` and `q_φ/r_a = H v̄ q₀ − b̃ H³𝓘₃ sin β sin πφ`, matching `transport.py` (T2)/(T3) term for term. Note the asymmetry is *predicted*, not fitted: `r_a` appears on the axial buoyant flux and not the azimuthal one. |
| Elliptic equation under caliper geometry | **exact — no residual** | See §4: my own first attempt found a spurious `O(dr_a/dξ)` residual. It does not exist. |
| Elliptic solver, single fluid | **machine precision** | Derived `w̄(φ) = H²⟨H⟩/⟨H³⟩` from `Hw̄ = H³𝓘₁G_ξ` with `G_ξ` uniform in `φ`, then compared: max relative error **1.3e-14 … 6.0e-14** for `e ∈ {0, 0.4, 0.8}`, `n_φ ∈ {20, 40}`, with `max\|v̄\| ≤ 9e-15` (`scratchpad/audit/a3_reductions.py`). This rules the elliptic solver out as a cause of the case-1 failure. |
| Maximum principle | **structural, not clipped** | `grep -rn "np.clip\|np.maximum\|np.minimum" d2dga/` finds no clipping of `c̄`. The update `c_new = c − λ·div/Hra` is a pure flux update. Measured `c̄ ∈ [0, 1]` on every run I made and every run in `output/`. The only `maximum` on a physical quantity is the documented NUM-13 mobility floor. |
| Conservation | **1e-14 – 1e-15** | Reported on every run in `output/`; reproduced in every run I made during this audit. |
| Repository test suite | **223 passed, 0 failed** | `pytest tests/`, run during this audit. See §6 on why that is not reassurance. |
| Endpoint identities | **exact** | `q₀(0) = 0`, `q₀(1) − 1 = 0`, `𝓘₃(0) = 𝓘₃(1) = 0` to **0.0e+00** for `m ∈ {0.006, 0.2, 1, 5, 160}`, and on the K-GEP-1 table over **all** other axes. Table `q₀` monotone in `c̄`; table `𝓘₃ ≤ −1.1e-4 < 0` strictly on `0 < c̄ < 1` (`scratchpad/audit/a6_endpoints.py`). |
| K-GEP-1's "single shock at the mean speed" | **confirmed** | `max_c f(c̄)/c̄ = 1.00000` from both the table and a direct gap-scale solve. |
| ZF23 (3.6) thresholds | **correct** | `σ_{w+r} > 0.08` **and** `\|w̄_r^+\| > 0.05`, read at ZF23 p. 1277 of the extraction, matches `postprocess.DispersionMetrics.is_dispersive`. |
| ZF22 Table 3 transcription | **correct** | All ten rows of `tests/benchmarks/zf22_cases.py` checked against ZF22 p. 27. |

---

## 3. Verdict per M-gate

| Gate | Verdict | Basis |
|---|---|---|
| **M0** geometry / validity envelope | **Pass** | `H(φ,ξ)` and `δ/π` re-derived and checked; `δ/π < 1/π` bound is arithmetic and the register reports the spec's 0.34 as unattainable, correctly. The caliper envelope is honestly reported as strained. |
| **M1** scaling | **Pass** | CONV-01, CONV-02 verified against both papers independently. `Fr*`, `b`, `Δρ` orderings check out. |
| **M2** closures | **Pass, with the strongest evidence in the report** | `𝓘₁, 𝓘₂, q₀, 𝓘₃` all re-derived from the velocity profile and matching to sympy-zero. CONV-03's sign confirmed by a new route. |
| **M3** closure table | **Pass with reservations** | Endpoint identities exact, monotone, sign-definite, interpolation ≤1 % of the flux span. But: NUM-03 still "pending"; range guard disabled in production (A-3); `derivative_bounds` unguarded (B-2); tip region poorly resolved (Q-1). |
| **M4** elliptic | **Pass** | Reproduces the analytic axial profile to 1e-14 at `e = 0.8`; the equation itself re-derived as exact for `H(φ,ξ)`, `r_a(ξ)`, `β(ξ)`. |
| **M5** transport | **Pass** | Conservation to 1e-14, maximum principle holds and is structural, `r_a` placement independently confirmed, converges to the exact rarefaction. NUM-15/16/26/29 are all real defects in the printed scheme, correctly identified. Docstring contains one stale self-contradiction (B-3). |
| **M6** coupled loop | **Conditional** | Nine of ten ZF22 cases reproduced; case 1 is an open failure in the only configuration that tests adverse buoyancy (A-4), and its stated mesh convergence does not hold for `η_E` (A-5). Mesh independence is asserted on a Newtonian pair at mild aspect ratio and not at the production configuration (A-2). The outflow condition imports fluid (A-1). |
| **M7** post-processing | **Pass** | Volume-weighted averaging is right; ZF23 thresholds correct; `σ_{w−r}` bin dependence honestly flagged as ~15 %. |
| **K-GEP-1 production result** | **Do not quote `η_E` to four figures** | `η_E = 0.9946` is reported to 1e-4. A-1 bounds an unmeasured upward bias at a few percent, A-2 shows the mesh is not converged for front-position quantities and no refinement run exists, A-4 notes the HB path has no published validation. The *qualitative* conclusion — a near-piston displacement, dominated by a strongly favourable buoyancy — is well supported by the 1-D theory (`max f/c̄ = 1` exactly) and is not in doubt. The register's own FLU-01…06 caveats remain the largest uncertainty of all, and they are correctly stated. |

---

## 4. Where I was wrong

I opened a finding that the elliptic equation inherits an `O(dr_a/dξ)`
inconsistency from BF25, on the grounds that BF25 solves `∇_a·[r_a X] = 0` where
the exact condition is `∇_a·[X] = 0`, and measured `max|dr_a/dξ| = 0.52` on the
caliper wall to argue it was not small.

**That was wrong.** The single-valuedness condition is
`∂_ξ(∂_φ p) − ∂_φ(∂_ξ p) = 0` in plain `(φ, ξ)` coordinates, and `∂_φ p = r_a·(∇_a p)_φ`
carries the `r_a` exactly. Carrying it through gives the code's operator with no
residual whatever, for arbitrary `H(φ,ξ)`, `r_a(ξ)` and `β(ξ)`. `docs/derivation.md`
says so, and it is right. I had used a rotation identity that is not valid in
this metric.

The measurement stands as a fact about the geometry — `max|dr_a/dξ| = 0.516`
resolved and `0.093` on the 16 × 80 mesh for the caliper wall, against `0.0094`
and `0.0092` for the synthetic one — and it is relevant to the *upstream*
Hele-Shaw reduction, which is where `docs/derivation.md` correctly locates the
slow-variation assumption. It is not an elliptic-solver finding.

I also briefly over-read the closure-table interpolation error as 87 % on `𝓘₃`.
In absolute terms it is 0.9 % of the flux span; the 87 % was a relative error at
`c̄ = 0.017` where `𝓘₃ ≈ −0.003`. Q-1 states it correctly.

---

## 5. The three things most likely to be wrong

1. **The outflow boundary (A-1).** It is the same class of error as NUM-29,
   which cost the build an entire wrong result before it was caught, and it has
   the same signature — a flux term that scales with `b H³𝓘₃` and is therefore
   invisible on the published benchmarks but enormous on K-GEP-1. It is the one
   finding here that can move `η_E` in the direction that flatters the answer,
   and it was measured at 3.5 % of a job on a case the repository already runs.
2. **ZF22 case 1 (A-4, A-5, BENCH-09).** Four hypotheses have now been refuted
   and a fifth — mesh convergence — turns out not to have been excluded at all
   (A-5: `η_E` moves +12.8 % on one refinement, toward ZF22). The failure is in
   the 2-D response to an adverse buoyancy source, and it is worth 4.5× in the
   front speed against 1-D theory times the analytic channel factor. Until it is
   explained,
   the possibility that the 2-D buoyancy coupling is wrong in a way that happens
   to be benign for `b > 0` cannot be dismissed — and every quoted result has
   `b > 0`. I would put the next effort here, specifically on NUM-05/Q6 (the
   `φ`-symmetry assumption), which is the one structural choice BENCH-09 has not
   yet tested and the only one that plausibly differs between this code and
   ZF22's D2DGA in an unstable configuration.
3. **The Herschel–Bulkley path as a whole (A-4.3, B-1, B-2, Q-1, A-3).** Five
   separate findings in this report land on it, and not one of them is caught by
   any benchmark, because every benchmark in the repository is Newtonian. The
   parallel-`ū`/`G̃_b` assumption in particular is an approximation that is
   documented for `β ≠ 0` and undocumented for `β = 0`, where it is live.

---

## 6. What could not be checked

Stated as prominently as the findings, because several of these are larger than
anything above.

* **The measured-caliper result does not exist.** `--wall caliper` has been
  smoke-tested only. Its `δ/π = 0.168` and `|dr̂_o/dξ̂| = 0.163` are an order of
  magnitude beyond every validated case, so the two conventions that were
  derived rather than transcribed *precisely because caliper geometry can
  distinguish them* — CONV-04 and CONV-10 — remain untested by any run. I could
  not run it: ~16 h at 16 × 80.
* **A-1's effect on the K-GEP-1 numbers was bounded, not measured.** I measured
  the instantaneous bound (−4.61 Q) and estimated a per-crossing volume (~5.6 %),
  but measuring the real cumulative import needs a full 3.3 h run.
* **No mesh-refinement run of K-GEP-1 was performed**, by me or by the build. The
  only attempt died at 2.5 %.
* **B-1 was bounded only to `|ū| = 50`**, against a table axis reaching 3000. If
  the caliper run drives `|ū|` an order of magnitude higher — and its `H³` is 10×
  larger, so it will — the 2.6 % becomes an unknown.
* **B-2 was not measured.** I established that `derivative_bounds` cannot detect
  out-of-range queries; I did not establish whether an extrapolated bound ever
  fell *below* the true slope, which is the direction that would break
  monotonicity.
* **The test suite passes and that is not reassurance.** `pytest tests/` is
  **223 passed, 0 failed**, run during this audit. None of A-1 (the outflow
  face), A-2 (front position vs an exact solution), A-5 (case-1 `η_E` mesh
  sensitivity) or B-1 (the parallel-`ū` assumption) is covered by any test, and
  a green suite alongside this report is the point: the gaps are in what is
  asserted, not in whether the assertions hold.
* **The `1/Fr*²` term** (CONV-06) is dropped for `β = 0` on the grounds that
  `∇_a·f ≡ 0`. I verified that algebra
  (`∇_a·f = sin πφ[sin β dr_a/dξ + r_a cos β dβ/dξ]`) and it is exact. For any
  inclined caliper well it must be retained; nothing has been run in that regime.
* **The Materials 2025 fluid properties** were taken on trust from the register's
  transcription; I did not re-read that paper's Table 1.
* **`n_y = 200` for the production table** (against `N_y = 400` in the NUM-02
  register row) was not convergence-tested by me.
* **Nothing in this audit tests the physics against experiment.** ZF22/ZF23
  compare against 3-D simulation and experiment; this repository compares against
  ZF22's D2DGA column. Agreement with ZF22's D2DGA is agreement with a model, not
  with a well.
