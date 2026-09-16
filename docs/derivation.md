# Section 3.3 — does the pressure elimination survive caliper geometry?

**Verdict: yes, exactly. No residual terms. Proceed.**

Reproduce with `python3 tests/derivation_check.py`.

## Question

BF25 eliminates the pressure between the two components of (2.13) to reach the
elliptic equation (2.16), `∇a·[S + b] = 0`. That algebra was written for a
uniform annulus: `r_a = 1`, `H = H(φ)`. The caliper modification makes
`H = H(φ,ξ)` and `r_a = r_a(ξ)`. The build spec requires us to confirm the step
still holds and to identify the order of any residual terms.

## Method

`tests/derivation_check.py` carries `r_a(ξ)`, `H(φ,ξ)`, `β(ξ)`, `c̄(φ,ξ)`,
`𝓘₁(φ,ξ)`, `𝓘₂(φ,ξ)` as unevaluated sympy functions — nothing is assumed
constant, nothing is linearised. It starts from BF25 (2.13), substitutes the
stream function (BF25 factor-of-2 convention), solves each momentum component
for its pressure-gradient component, and imposes single-valuedness of `p` via
`∂ξ(∂φp) − ∂φ(∂ξp) = 0`.

## Result

    residual  =  [∂ξ(p_φ) − ∂φ(p_ξ)]  −  ∇a·[S + b_vec]  ≡  0

Identically zero in sympy, with no simplifying substitutions.

**The cross-differentiation is exact algebra.** It is not an approximation and
it does not generate O(δ) terms. Two structural facts make it work:

1. `r_a` depends on ξ only, never on φ. So `r_a` commutes with `∂/∂φ`, and
   `(1/r_a)∂φ[S_φ]` collapses to `∂φ[Ψ_φ/(2 r_a I₁)]` with no leftover.
2. The buoyant-mobility ratio is gap-width-free:
   `I₂/(H·I₁) = H⁴𝓘₂/(H·H³𝓘₁) = 𝓘₂/𝓘₁`, confirmed symbolically. This is the
   BF25 Appendix A4 rescaling doing its job, and it is what lets one closure
   table serve all depths (build spec §3.2).

### Where the slow-variation assumption actually lives

Not here. It lives **upstream**, in the reduction of Navier–Stokes to the 2D
shear flow BF25 (2.2)–(2.3) / B02 §2.1, where terms of O(δ/π) are dropped and
the annulus is taken as locally uniform in ξ. Once (2.13) is granted, everything
downstream is exact for arbitrary `H(φ,ξ)` and `r_a(ξ)`.

The practical consequence: the caliper modification cannot be validated by
checking the elliptic equation. Its validity is entirely controlled by whether
the *upstream* reduction still holds — i.e. by δ/π and by `dr̂_o/dξ̂`. That is
precisely what the M0-T5 validity envelope and the §3.5 A–L sweep are for. This
is the same footing B02 stood on when it simulated a washout in §4.2.3.

## Finding 1 — `∇a·f` does not vanish in general, but does vanish for K-GEP-1

BF25 (2.18) drops the `1/Fr*²` term on the grounds that `∇a·f = 0`. Symbolically,
for variable geometry:

    ∇a·f = [ r_a cos β · dβ/dξ  +  sin β · dr_a/dξ ] · sin πφ

Two separate ways this can be non-zero under a caliper wall: a building
inclination (`dβ/dξ ≠ 0`) and, new to this work, an axially varying mean radius
(`dr_a/dξ ≠ 0`) in an inclined well.

**At β = 0 — K-GEP-1, vertical — both terms carry a `sin β` or a `dβ/dξ` factor
and `∇a·f ≡ 0` exactly, however violently `r_a` varies.** So for the production
case the term is legitimately dropped, and it is dropped for a stronger reason
than in BF25: not "small", but exactly zero.

For any inclined well with a caliper wall the term must be retained. The code
keeps `β` as a field rather than a constant and evaluates `∇a·f` numerically at
M0-T4 rather than assuming it away.

## Finding 2 — the symbol `b` means two different things in BF25

This one is invisible in the published work and bites only under caliper geometry.

BF25 (2.12) defines the scalar

    b = −Δρ̂ / (r_a Fr*²)                                    ... (2.12)

The `1/r_a` there is **correct for `G_b`**: it cancels against the `r_a` carried
inside `f = (r_a cos β, r_a sin β sin πφ)`, leaving `G_b` free of `r_a`, as a
modified-pressure-gradient must be.

But the scalar multiplying `f` in the *elliptic* equation is, from this
derivation,

    b_s = (1/Fr*²) [ 1 − Δρ̂ ( c̄ + 𝓘₂/𝓘₁ ) ]

which is **free of `r_a`**. Relative to (2.12) that is a factor `r_a`:

    b_s = 1/Fr*²  +  ( r_a · b_(2.12) ) · ( c̄ + 𝓘₂/𝓘₁ )

whereas BF25 (2.18) as printed reads `1/Fr*² + b·(c̄ + 𝓘₂/𝓘₁)` with the (2.12)
`b`. **At `r_a = 1` — every published case — the two agree and the ambiguity is
undetectable.** With a caliper wall `r_a` varies along the well and they do not.

The same latent factor appears between BF25 (2.21) and BCF25 (5) for the
buoyancy flux; both are self-consistent once `b` is read as the `r_a`-free
scalar paired with an `f` that carries `r_a`.

**Decision:** this code uses the derived `r_a`-explicit forms throughout and
never the printed (2.18). The full vector is

    b_vec = (1/Fr*²) [ 1 − Δρ̂ ( c̄ + 𝓘₂/𝓘₁ ) ] · ( r_a cos β , r_a sin β sin πφ )

Reduction to BF25 (2.18) at `r_a = 1` is asserted as a unit test. Logged as
CONV-04 in `docs/assumptions.md`.
