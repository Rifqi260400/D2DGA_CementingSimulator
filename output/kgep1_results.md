# K-GEP-1 displacement — D2DGA result

Synthetic wall (A = 1.5 in, L = 20 m, enlargement), 16 x 80, CFL 0.5, run to
1.2 pumped volumes.  Fluids: Table 1 of *Materials* 2025, 18, 3098 — cement
slurry (Herschel–Bulkley, 1200 kg/m³, κ = 0.6 Pa·sⁿ, n = 0.4, τ_Y = 1.4 Pa)
displacing the paper's drilling fluid (Newtonian water, 998 kg/m³, 1 mPa·s),
at ŵ₀ = 0.2 m/s.  Dimensionless: `m = 0.00638`, `B = 11.94`, `b = 27.90`,
`Fr* = 0.0852`, `δ/π = 0.0825`, `Z = 564.1`.

| | **B02 (70)** | **uniform inflow** |
|---|---|---|
| steps / wall time | 85 474 / 3.3 h | 85 565 / 3.3 h |
| volume conservation | 2.9e-14 | 1.8e-14 |
| `t_br` @ 0.01 / 0.1 / 0.5 | 0.953 / 0.986 / 0.988 | 0.959 / 0.991 / 0.993 |
| **`η_E` at 1.2 volumes** | **0.9946** | **0.9945** |
| narrow-side minimum `c̄` | 0.9494 | 0.9486 |
| residual (`c̄ < 0.5`) | 0.0000 | 0.0000 |
| ZF23 `σ_w+r` / `\|w_r+\|` | 0.0200 / 0.0031 | 0.0209 / 0.0018 |
| ZF23 (3.6) classification | not dispersive | not dispersive |
| static cells (NUM-13) | 0 | 0 |
| Picard steps at the cap | 0 | 0 |

## What it says

**A near-perfect displacement.** Breakthrough at 0.99 of the piston time, 99.5%
of the annulus displaced, the narrow side left at 95% cement, and no cell
anywhere below `c̄ = 0.5`.  By ZF23 (3.6a,b) the front is not dispersive.

**It is the predicted answer, not a surprise.** BF25 §3.1's planar theory gives
the flux function `q₀ + b𝓘₃` reaching **−25** in mid-range — far below the chord
from (0,0) to (1,1) — so the entropy solution is a single shock at the mean
speed and `t_br → 1`. Measured 0.988–0.993.  The 2-D solver and the 1-D theory
agree without being fitted to each other.

**The bottom-hole condition does not matter here.** `η_E` differs by 0.0001 and
`t_br` by 0.005 between B02 (70) and uniform inflow.  NUM-04b's measured 0.005-volume
leak under B02 (70) is boundary bookkeeping, not a change in the displacement —
and the contrast with ZF22 case 1, where the two conditions give *completely*
different answers, locates exactly when the choice matters: when the buoyancy is
adverse and the front is unstable.

## ⚠️ What it does not say

**This is a water displacement, not a mud displacement.**  The source paper's
"drilling fluid" is water.  That gives `m = 0.0064` — the cement is 160× more
viscous than the fluid it displaces — and a 202 kg/m³ *favourable* density
difference.  Both are strongly stabilising, and `b·𝓘₁ ≈ 1500` means the
gap-scale buoyancy exceeds the pumping by three orders of magnitude (FLU-06).
A displacement this clean is the arithmetic consequence of those inputs.

A real drilling mud is 1100–1600 kg/m³ and carries a yield stress: `m` would rise
by two orders of magnitude, `Δρ̂` would shrink or reverse, and the mud channel
that `narrow-side minimum = 0.95` says is absent here is precisely what a real
mud leaves behind.  **`η_E = 0.994` is an upper bound on what K-GEP-1 would
achieve, not a prediction of it.**  See assumptions.md FLU-01..06.

Secondary caveats: the pump rate is a reading of that paper's inlet velocities,
not a K-GEP-1 measurement (FLU-04); the wall is synthetic, and the measured
caliper is a harder problem still (`δ/π` 0.168 against 0.083, `r_a` varying 94%
against 19%) which has been smoke-tested but not run to completion.
