# Blocked

Two kinds of item, handled differently (remediation brief §5).

* **§5a — missing physical data.** Cannot be derived, computed or reasoned out.
  A provisional value is in use so that work continues; each is a single
  clearly-marked constant so it can be changed in one place.
* **§5b — irreducible ambiguity.** The papers do not determine the answer, the
  derivation was attempted, and at least one numerical experiment was run to
  discriminate. The most defensible option is implemented, the alternative is
  behind a flag, and the difference is quantified.

---

## §5a — Missing physical data

### BLK-1 — Centralizer record for K-GEP-1 (dominates the answer)

**What is needed.** Centralizer type, spacing and placement over the 194–390 m
open hole, or a measured or computed standoff profile.

**Why it cannot be substituted.** `docs/assumptions.md` GEO-02 already measures
the sensitivity and it is the largest in the model: sweeping `e_gauge` from 0.1
to 0.5 collapses the far-field narrow-side velocity fraction
`(1−e)²/(1+1.5e²)` from **0.643 to 0.005 — a factor of 130** on the synthetic
wall, and 0.798 → 0.182 on the washout wall. Narrow-side displacement is the
whole question in primary cementing, so no efficiency number is meaningful until
this is fixed.

**Provisional value in use.** `e_gauge = 0.3` (≈70 % standoff), imposed as a
constant physical offset. Single constant:
`d2dga/config.py` → `GeometryConfig.eccentricity_gauge`.

**Status.** Unchanged by this remediation. Recorded here because it was in
`assumptions.md` only, where it is one row among sixty.

### BLK-2 — Pump schedule for K-GEP-1

**What is needed.** The actual pump rate(s) and their duration for the cement
job.

**Provisional value in use.** `ŵ₀ ∈ {0.05, 0.2, 0.5} m/s`, read from the inlet
velocities of *Materials* 2025, 18, 3098 — a different, much smaller annulus.
This is a reading, not a measurement (`assumptions.md` FLU-04). Single constant:
`d2dga/config.py` → `FluidsConfig.velocities`.

**Note.** FLU-06 measures that this displacement is **not pump-rate controlled**
(`b·𝓘₁ ≈ 1500`), so the schedule matters less here than BLK-1 or BLK-3.

### BLK-3 — Real drilling-mud properties

**What is needed.** Density and rheology of the mud actually in the hole.

**Why it matters more than any numerical finding in this report.** The source
paper's "drilling fluid" is **water** (998 kg/m³, 1 mPa·s, Newtonian), giving
`m ≈ 0.006` — the cement 160× more viscous than what it displaces — and a
202 kg/m³ *favourable* density difference. Both are strongly stabilising. A real
mud is 1100–1600 kg/m³ with a yield stress: `m` would rise by two orders of
magnitude and `Δρ̂` would shrink or reverse. `η_E = 0.994` is an **upper bound**,
not a prediction. `assumptions.md` FLU-01, FLU-06.

**Provisional values in use.** `d2dga/config.py` → `FluidsConfig` (mud
998 kg/m³ / 1e-3 Pa·s; cement 1200 kg/m³, κ = 0.6, n = 0.4, τ_Y = 1.4).
Re-affirmed by the user on 2026-09-17.

### BLK-4 — Casing ID / weight / grade, for the caliper bias correction

**What is needed.** The 10″ casing's exact ID, to fix the caliper tool's bias.

**Effect.** `assumptions.md` CAL-01: the cased section reads 10.432 in ± 0.007,
which is the casing ID, not rock. Against a nominal 10.00 in ID the tool carries
+0.43 in (+4.3 %); against 10.19 in, +0.24 in. True open-hole gauge is therefore
≈10.07–10.12 in rather than the logged 10.51 in. The bias propagates to `H`
cubed in every flux.

**Provisional handling.** The log is used as recorded, with the bias bracketed
rather than applied. Single constant: `min_diameter_m` in
`scripts/kgep1_run.py`'s `make_geometry`.

### BLK-5 — `D2DGA_BUILD_SPEC.md` does not exist in the repository

**What is needed.** The build spec itself, if the gate list is to be
authoritative.

**Why this is a blocker and not a nit.** The remediation brief names the spec as
a source and requires `docs/gate_status.md` to carry "every M-gate from the build
spec". The file is not in the repository and is not in git history — it was
supplied in the conversation that produced the code. `docs/gate_status.md` is
therefore reconstructed from the `M<n>-T<k>` markers that survive in `tests/`
and `docs/assumptions.md`. **A gate that the spec required but that was never
given a marker is invisible to this reconstruction, and its absence cannot be
detected from inside the repository.**

One gap is already visible from the reconstruction itself: **M0-T5** is
referenced in `docs/derivation.md` and implemented as a *script*
(`scripts/caliper_report.py`, which produces a figure) but has **no test**, so it
asserts nothing and cannot regress. Recorded, not fixed: without the spec there
is no statement of what M0-T5 was supposed to assert.

---

## §5b — Irreducible ambiguity

### BLK-6 — ZF22 case 1 (`BENCH-09` / audit A-4, A-5): still not reproduced

**Correct behaviour.** ZF22 Table 3 case 1 (`e = 0.8`, `m = 0.2`,
`b_ZF22 = −50`): `t_br = 0.44`, `η_E = 0.66`. This code gives `t_br ≈ 0.05`,
`η_E = 0.366` (20 × 200) rising to 0.412 (40 × 400).

**What has been ruled out** — four hypotheses, each by measurement:

1. *The gap-scale closures.* BF25 §3.1's 1-D leading wave for these parameters
   gives `t_br = 0.391` against ZF22's 0.44. The closures agree with ZF22.
2. *The elliptic solver.* It reproduces the analytic axial profile
   `w̄(φ) = H²⟨H⟩/⟨H³⟩` to **1.3e-14 … 6.0e-14** at `e = 0.8` (audit §2).
3. *Eccentric channelling.* At `b = 0, e = 0.8` the solver gives
   `t_br = 0.369–0.380` against the analytic wide-side prediction
   `1/(1.5 × 1.6531) = 0.4033`; the ~7 % shortfall is exactly the A-2 tip bias.
   Channelling is right.
4. *The bottom-hole condition and inlet buoyant exchange* — both refuted in
   `assumptions.md` BENCH-09.

**What this session added.**

5. *Mesh convergence was never established for `η_E`, and it does not hold*
   (A-5). At fixed `n_xi = 200`, azimuthal refinement gives
   `η_E` = 0.3657 → 0.3939 → **0.4773** for `n_phi` = 20 → 40 → 80. **The
   increments grow — +0.028 then +0.083 — so the sequence is diverging, not
   converging.** BENCH-09's claim "(i) it is mesh-converged" is correct for
   `t_br` and false for `η_E`, in the one direction it never refined.
6. *An analytic reason, from the paper.* BF25 §3.3's Muskat criterion
   (gate **BENCH-10**, added this session) classifies **case 1 as the only
   Muskat-unstable case of the ten**: `Δw(0⁺) = +0.97` against `−1.96 … −162.9`
   for the other nine. A Muskat-unstable front has **no mesh-independent finger
   width**, so `η_E` at a fixed time is not a converged functional of the mesh.
   The measured azimuthal mesh-dependence is what an unstable base state
   produces.

**Best hypothesis, and what I now doubt.** The disagreement is not a localised
coding error — every component it could live in has been pinned independently to
machine precision. It is more likely that ZF22's D2DGA computation of case 1 and
this one are *both* resolving an unstable flow and are not obliged to agree on a
fixed-time functional of it. What I doubt most is the implicit assumption, in
both the build and the audit, that `η_E` at `t = 1.2` is a well-posed target for
this case. ZF22 themselves single case 1 out (*"we verify again via the poor
efficiency of case 1 that `b < 0` should be strictly avoided"*) and report it to
two significant figures.

**Evidence that would distinguish the remaining possibilities**, in order of
value:

* ~~A φ-refinement sequence at fixed `n_xi`.~~ **Done this session, and it
  answers the question:** `η_E` climbs with growing increments through
  `n_phi` = 20, 40, 80 (0.3657, 0.3939, 0.4773). The functional is not
  converged and **no agreement with a fixed-time `η_E` should be expected** for
  this case. What remains open is only whether ZF22's own computation was
  equally mesh-dependent, which their paper does not report.
* Running the full annulus with azimuthal periodicity instead of the
  half-annulus symmetry (`assumptions.md` NUM-05/Q6). ZF22 §5 notes periodicity
  permits azimuthal asymmetry. This is the one structural choice never tested
  and the only one that plausibly differs from ZF22 in an unstable
  configuration. **Not attempted here: it is a new capability, and the
  remediation brief §8 forbids adding features.**
* ZF22's own mesh and time-integration for case 1, which the paper does not give.

**Not fixed. K-GEP-1 has `b > 0` and is Muskat-stable by a wide margin
(`Δw(0⁺) < −10`), so this does not affect the well being modelled** — but the
stop condition "ZF22 Table 3 reproduces for all ten cases" is **not met**, and
that is stated plainly rather than worked around.

### BLK-7 — BF25 (3.13)'s stability range is degenerate as printed

**The ambiguity.** BF25 §3.3 says "For an unstable regime, `Δw(c₀) > 0` for all
`c₀ ∈ [0,1]`". Taken literally it can never classify anything: at `c₀ → 1`,
`w_f → q₀'(1) + b𝓘₃'(1) = 0` for **every** pair and every `b`, while
`w_finger → 1`, so `Δw(1⁻) → +1` unconditionally. "Stable" is unreachable and
"unstable" unfalsifiable at that end. All ten ZF22 cases come out "partial
penetration", including case 4 at `b = 1000`, which displaces as a near-piston.

**What was done (brief §5b).** Primary source read directly (BF25 §3.1 (3.8)
and §3.3 (3.10)–(3.13)); derivation attempted and it closes, which is how the
degeneracy was found; numerical experiment run on all ten ZF22 cases plus a
`b = 0` sweep in `m`.

**Resolution implemented.** The discriminating quantity is `Δw` at the **leading
edge** (`c₀ → 0`), which is where a finger must outrun the front to penetrate it.
This reading is pinned to the paper rather than to preference: at `b = 0` it
reduces to `Δw(0⁺) = m − 3/2`, so the finger penetrates exactly when `m > 3/2` —
and BF25 §3.2 gives **3/2 as `M₃^min`**, attributing it to Lajeunesse et al.
(1999) from an entirely different argument. Reproduced to 2.3e-9.

**The alternative is behind a flag and the difference is quantified.**
`classify(..., c0_max=1.0)` is BF25 verbatim and returns "partial penetration"
for all ten cases; `c0_max = 0.9` separates them into stable (2, 3, 4, 5, 6, 9,
10) and not-stable (1, 7, 8). Both behaviours are asserted in
`test_bench10_bf25_3_13_literal_range_is_degenerate`, so the choice is visible.

### BLK-8 — Whether an extrapolated derivative bound ever fell below the true slope

**What is unresolved.** B-2 established that `ClosureTable.derivative_bounds`
extrapolated silently and now records. It did **not** establish whether any
extrapolated bound came out *below* the true `|dq₀/dc̄|` or `|d𝓘₃/dc̄|`, which is
the direction that would break LLF monotonicity.

**Why it is not closed here.** Answering it needs a direct gap-scale solve at
every out-of-range query of a full run and a comparison of slopes — ~10⁴ extra AL
solves inside an 85 000-step run. The indirect evidence is strong (`c̄ ∈ [0,1]`
held on every run, and the out-of-range queries were 315 of 25 899 552, all on
`umag` below 0.02 where the closures are flat to six significant figures), but
indirect is not the same as measured, and after the A-3 fix the situation cannot
recur: the `umag` axis is now anchored at 0.

**What would settle it.** A run-mode that, on each recorded out-of-range query,
performs the direct solve and asserts `bound ≥ |slope|`.
