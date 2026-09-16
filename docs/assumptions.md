# Assumption register

One row per modelling choice not fully pinned down by the six source papers.
Updated as the build proceeds. IDs are referenced from code comments.

Status key — **Conflict**: sources disagree, a choice was forced.
**Silent**: sources agree but are silent. **Pending**: awaiting user data.

---

## Convention conversions (build spec §1.1)

| ID | Choice | Value used | Source status | Justification | Sensitivity tested? |
|---|---|---|---|---|---|
| CONV-01 | Stream-function definition | `∇aΨ = (2Hw̄, −2Hv̄)` | **Conflict.** B02 (39) and PF04 (5) use `∂Ψ/∂φ = r_a H w̄` — no factor 2. BF25 §2.1, ZF22 (2.2), ZF23 (2.7) carry the 2. | Spec §1.1 mandates BF25. Everything taken from B02/PF04 is rescaled at the point of use. Getting this wrong puts every velocity out by 2. | Enforced by M4-T4 (flux through any ξ-slice = Q) |
| CONV-02 | Buoyancy scalar | BF25 (2.12), `b = −Δρ̂/(r_a Fr*²)`, `Δρ̂ = ρ̂₁ − ρ̂₂`, `b > 0` for denser-displacing-lighter upward | **Conflict.** B02/PF04 use a Stokes number `St*`; ZF22 uses `b = Δρ̂ĝd̂²/(µ̂₁ŵ₀)` with the *opposite* density difference sign; BCF25 (24) uses a geometric-mean viscosity. | Spec §1.1 mandates BF25. ZF22 conversion is exercised at M1-T2. | M1-T2 (pending) |
| CONV-03 | `I₃` prefactor | BF25 (2.23) | **Conflict.** ZF22 (4.26) has denominator `2m[·]` with prefactor `Δρ H³/(6η₂)`; BCF25 (27) has `12[·]` with prefactor `b H³ r_a`. Neither is wrong — the difference sits in the prefactor. | Spec §1.1 mandates BF25 and requires the others to be shown to reduce to it. | M2-T3 (pending) — will report the conversion factor |
| CONV-04 | `r_a` in the elliptic buoyancy scalar | Derived `r_a`-explicit form: `b_vec = (1/Fr*²)[1 − Δρ̂(c̄ + 𝓘₂/𝓘₁)]·(r_a cos β, r_a sin β sin πφ)` | **Conflict, newly identified.** The symbol `b` in BF25 (2.12) and in BF25 (2.18) differ by a factor `r_a`. See `docs/derivation.md` Finding 2. | At `r_a = 1` (every published case) the two agree, so the paper cannot distinguish them. Under a caliper wall `r_a` varies and they diverge. Derived form used; reduction to (2.18) at `r_a = 1` asserted as a unit test. | Derivation is symbolic and exact; unit test pending M4 |
| CONV-05 | `Ψ(1,ξ,t) = Q(t)` from B02 (68) | Rescaled into CONV-01 convention before use | **Conflict** (consequence of CONV-01) | Same reason as CONV-01. | M4-T4 (pending) |
| CONV-06 | `1/Fr*²` term in `b` | Dropped for K-GEP-1 | **Silent→resolved.** BF25 (2.18) drops it citing `∇a·f = 0`. | `docs/derivation.md` Finding 1 shows `∇a·f = [r_a cos β dβ/dξ + sin β dr_a/dξ] sin πφ`, which is *identically* zero at β = 0 however `r_a` varies — a stronger justification than BF25's. **Must be retained for any inclined caliper well.** | **Yes — M0-T4 passes**, and the companion test confirms the diagnostic detects a non-zero value at β = π/4 |

---

## Geometry (build spec §3, §4)

| ID | Choice | Value used | Source status | Justification | Sensitivity tested? |
|---|---|---|---|---|---|
| GEO-01 | Casing OD | **177.8 mm (7 in)** | **Pending→resolved by user 2026-09-16.** Spec offered 127 mm (project notes) vs 7 in / 5½ in (earlier analysis). | User-selected. ⚠️ **Unresolved tension:** 7 in gives δ/π = 0.062 in gauge hole, but the spec states an expectation of ~0.12, which is matched only by 127 mm (0.112). The spec's two δ/π figures are mutually inconsistent — see GEO-07. Worth re-confirming against the well file before the production run. | Not yet — all three were costed at decision time (0.062 / 0.098 / 0.112 gauge) |
| GEO-02 | Casing eccentricity | `e_gauge = 0.3` (≈70% standoff), imposed as a **constant physical offset** `ê = 12.95 mm` | **Pending.** No K-GEP-1 centralizer record supplied. | Routine field value. Constant-`ê` is the spec-mandated default (§4) and is the physically correct behaviour: `e = ê/(2d̂)` falls inside a washout because the casing does not move when the hole enlarges. Model is pluggable (`EccentricityModel` ABC). | **No — needs a sweep.** Also controls the synthetic-wall feasibility limit (GEO-04) |
| GEO-03 | Synthetic amplitude basis | **Diameter** (`amplitude_is_diameter=True`) | **Conflict.** Spec §3.5 writes `gauge + A·sin(2πξ̂/L)` for the *radius*, but quotes "gauge hole 10.4 in" (a diameter) and "amplitude 1.5 in" side by side. | Reading both as diameters keeps them on the same footing. Radius reading leaves only a 0.2 in minimum gap at A = 1.5 in with 7 in casing. Flag exposed for the literal reading. | Switchable; not yet swept |
| GEO-04 | Synthetic wall mode | `symmetric` default, `enlargement` available | **Silent.** Spec writes a symmetric sinusoid but also demands "strictly positive annular gap everywhere" *and* a sweep to A = 6.0 in — mutually impossible. | **Quantified:** with 7 in casing at `e_gauge = 0.3` the symmetric form admits A < 2.38 in (`e < 1` limit) and A < 3.40 in (positive-gap limit). The spec's sweep {0.5, 1.5, **3.0**, **6.0**} in is **half infeasible**: A = 3.0 gives e = 2.55, A = 6.0 gives a negative gap. One-sided `enlargement` mode admits all four and is the physically honest model (holes wash out, they do not wash in). | **Yes — `test_m0_synthetic_wall_feasibility_boundary`** |
| GEO-05 | Caliper interpolation | Linear, clamped at ends | **Silent.** | Higher-order interpolation would manufacture detail the log does not resolve. Digitisation uncertainty ±0.3–0.5 in, features under ~0.5 m unresolved, possible clipping above 24 in. | Not applicable until the real log arrives |
| GEO-06 | Washout shape | Gaussian, centre ξ̂ = 184 m, σ = 5.0 m | **Silent.** Spec gives only "195–217 m reaching ≈23 in". | ±2σ spans depth 196–216 m, matching the stated interval. Smooth and C^∞, so it stresses the slow-variation assumption only through its gradient — which is what we want to measure — not through a discontinuity. | Not yet |
| GEO-07 | Washout δ/π expectation | Computed, **not** taken from spec | **Spec error.** Spec §3.4 states δ/π "rising to ~0.34 in the washout". | δ/π = (r̂_o−r̂_i)/[π(r̂_o+r̂_i)] is bounded above by **1/π ≈ 0.3183** for any r̂_i > 0, so 0.34 is unattainable by that formula. Computed value for a 23 in washout: **0.170** (7 in casing). Reported rather than reproduced. | n/a — arithmetic fact |

---

## Numerical (build spec §7)

| ID | Choice | Value used | Source status | Justification | Sensitivity tested? |
|---|---|---|---|---|---|
| NUM-01 (Q2) | CFL number | *pending* | BCF25 says "< 1", never gives a value | To be set at M5/M6; confirm via M6-T4 at 0.5 and 0.25 | Pending |
| NUM-02 (Q3) | AL parameters `r`, `ρ`, tolerance | *pending* | PF04 gives only `0 < ρ < r(1+√5)/2`; BF25 Figs 18–19 show convergence but not values | To be tuned at M2; BF25 notes sensitivity to the initial guess — try continuation from the previous timestep | Pending |
| NUM-03 (Q4) | Closure table resolution | *pending* | Silent | Refine until M3-T1 passes with margin | Pending |
| NUM-04 (Q5) | Bottom-hole inflow condition | Both B02 variants implemented, selectable | B02 §3.2 calls its own choice uncertain | B02's justification (entry effects negligible in long annuli) is **weaker here**: their well is 1000 m, K-GEP-1's open hole is 196 m | Pending — will report the difference |
| NUM-05 (Q6) | φ symmetry vs periodicity | Symmetry (default) | All published 2D work imposes symmetry | ZF22 §5 notes periodicity permits azimuthal asymmetry and is a small change. Hook to be left. | Pending |
| NUM-06 (Q7) | Initial interface at t = 0 | *pending* | Silent | Compare sharp step at ξ = 0 vs one-cell smoothed step; confirm the difference decays | Pending |
| NUM-07 | Grid aspect ratio | `n_phi=20, n_xi=400` (BCF25's) | BCF25 §IV | ⚠️ On K-GEP-1, Z ≈ 564, so `dξ = 1.41` against `dφ = 0.05` — a cell aspect ratio of ~28. BCF25's own case had Z much smaller. Flagged for M6-T3. | Pending — M6-T3 |
