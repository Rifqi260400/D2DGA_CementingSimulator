"""
The original 2DGA closure of B02 / PF04.

Two reasons this belongs in the codebase even though D2DGA is the target model:

  1. PF04's analytic steady states (38) and (40)-(43) are 2DGA results, so
     reproducing them -- test gates M4-T2 and M4-T3 -- requires this closure.
  2. The whole point of the thesis is that D2DGA improves on 2DGA, and
     ZF22/ZF23 compare the two throughout.  Without 2DGA there is no baseline.

The closure, B02 (57) / PF04 (8), defines chi implicitly from the areal flux:

    |grad Psi| = H^(m+2) / (kappa^m (m+2))
                 * chi^(m+1) / (chi + tau_Y/H)^2
                 * [ chi + (m+2) tau_Y / ((m+1) H) ]

with m = 1/n the INVERSE power-law index, and the modified pressure gradient
is then G = chi + tau_Y/H, B02 (58).

This is algebraically identical to B02 (53), the plane-Poiseuille law already
used to validate the augmented Lagrangian solver at M2-T2 -- substituting
G = chi + tau_Y/H into (53) reproduces (8) term for term.  So the closure
arrives pre-verified rather than on trust.

CONVENTION WARNING.  B02 and PF04 define the stream function WITHOUT the factor
of 2 that BF25/ZF22/ZF23 carry (CONV-01): B02 (39) has dPsi/dphi = r_a H w_bar
where BF25 has 2 H r_a w_bar.  Every |grad Psi| passed to `chi` must therefore
be the BF25 value HALVED.  `chi_from_bf25_gradient` does that conversion in one
place so it cannot be forgotten at a call site.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq


def areal_flux(chi, H, kappa, m, tau_y):
    """PF04 (8): |grad Psi| as a function of chi.  Monotone increasing in chi."""
    G = chi + tau_y / H
    return (H ** (m + 2) / (kappa ** m * (m + 2))
            * chi ** (m + 1) / G ** 2
            * (chi + (m + 2) * tau_y / ((m + 1) * H)))


def chi(grad_psi, H, kappa, n, tau_y, xtol=1e-14):
    """
    Invert PF04 (8) for chi, given the PF04-convention |grad Psi|.

    All of grad_psi, H, kappa, n and tau_y broadcast together, so the
    rheological parameters may vary cell by cell -- which they do once the two
    fluids are mixed by concentration.

    chi <= 0 means the fluid has not yielded and nothing flows; we return 0,
    which makes the modified pressure gradient G = tau_Y/H, the yield value.
    """
    g, Hs, ks, ns, ts = np.broadcast_arrays(
        *[np.asarray(x, dtype=float) for x in (grad_psi, H, kappa, n, tau_y)])
    out = np.zeros(g.shape, dtype=float)
    flat = [a.ravel() for a in (g, Hs, ks, ns, ts)]
    of = out.ravel()
    for i in range(of.size):
        gi, Hi, ki, ni, ti = (a[i] for a in flat)
        if gi <= 0.0:
            continue
        mi = 1.0 / ni
        hi = 1.0
        while areal_flux(hi, Hi, ki, mi, ti) < gi:
            hi *= 2.0
            if hi > 1e12:
                break
        of[i] = brentq(lambda c: areal_flux(c, Hi, ki, mi, ti) - gi,
                       1e-300, hi, xtol=xtol, rtol=8.9e-16)
    return out


def chi_from_bf25_gradient(grad_psi_bf25, H, kappa, n, tau_y):
    """CONV-01 bridge: BF25's grad Psi is twice B02/PF04's."""
    return chi(0.5 * np.asarray(grad_psi_bf25, dtype=float), H, kappa, n, tau_y)


def dchi_dgrad(grad_psi, H, kappa, n, tau_y, rel=1e-6):
    """chi'(|grad Psi|), needed for PF04 (41)'s alpha_k."""
    g = np.asarray(grad_psi, dtype=float)
    h = np.maximum(rel * np.abs(g), 1e-10)
    return (chi(g + h, H, kappa, n, tau_y) - chi(g - h, H, kappa, n, tau_y)) / (2 * h)


def effective_script_I1(grad_psi_bf25, H, kappa, n, tau_y):
    """
    2DGA expressed in the same interface the D2DGA solver uses.

    D2DGA has S_phi = Psi_phi / (2 I1) with I1 = H^3 script_I1.
    2DGA  has S_phi = (chi + tau_Y/H) Psi_phi / |grad_a Psi|.

    Matching the two gives an EFFECTIVE mobility

        script_I1_eff = |grad_a Psi| / ( 2 H^3 (chi + tau_Y/H) )

    so the same elliptic assembly serves both models.  It is velocity
    dependent, so the problem is nonlinear and needs Picard iteration --
    unlike Newtonian D2DGA, which ZF23 Section 2.2.2 notes is linear.

    The 2DGA buoyancy has no layering term, i.e. script_I2 = 0 identically;
    that single difference is what separates the two models' b vectors.
    """
    g = np.asarray(grad_psi_bf25, dtype=float)
    c = chi_from_bf25_gradient(g, H, kappa, n, tau_y)
    G = c + tau_y / np.asarray(H, dtype=float)
    return np.where(g > 0, g / (2.0 * np.asarray(H) ** 3 * np.maximum(G, 1e-300)), 0.0)


# ---------------------------------------------------------------------------
# PF04 (40)-(43): mildly eccentric steady state
# ---------------------------------------------------------------------------
def P(chi_val, tau_y, m):
    """PF04's strictly positive function, printed immediately after (43)."""
    return (((m + 1) ** 2 * chi_val ** 2
             + (m + 2) * (2 * m + 1) * chi_val * tau_y
             + (m + 1) * (m + 2) * tau_y ** 2)
            / (chi_val * ((m + 1) * chi_val + (m + 2) * tau_y)))


def alpha(chi_1, dchi_1, tau_y):
    """PF04 (41):  alpha^2 = pi^2 chi'(1,1) / (chi(1,1) + tau_Y) > 0."""
    return np.sqrt(np.pi ** 2 * dchi_1 / (chi_1 + tau_y))


def steady_interface_eccentric(phi, e, b, beta, L, fluids):
    """
    PF04 (40).  Interface shape for a mildly eccentric annulus, to O(e).

    `fluids` is ((chi1, dchi1, tau1Y, m1), (chi2, dchi2, tau2Y, m2)) with chi
    and chi' evaluated at (|grad Psi|, H) = (1, 1), i.e. linearised about the
    concentric geometry, exactly as PF04 states below (40).

    The denominator is the SAME for both terms: [chi_k + tau_kY]^2_1 + b cos(beta),
    where [.]^2_1 is PF04's jump, fluid 2 minus fluid 1.  The first term is (38).
    """
    (c1, d1, t1, m1), (c2, d2, t2, m2) = fluids
    G1, G2 = c1 + t1, c2 + t2
    denom = (G2 - G1) + b * np.cos(beta)

    first = -(1.0 / np.pi) * b * np.sin(beta) * np.cos(np.pi * phi) / denom

    a1, a2 = alpha(c1, d1, t1), alpha(c2, d2, t2)
    total = (P(c1, t1, m1) * G1 * a1 * np.tanh(a1 * L)
             + P(c2, t2, m2) * G2 * a2 * np.tanh(a2 * L))
    second = -(e / np.pi ** 2) * total * np.cos(np.pi * phi) / denom
    return first + second


def moving_frame_stream_function(phi, z, e, L, fluids):
    """
    PF04 (42)-(43), the O(e) moving-frame stream function.

        Phi ~ e P(chi1, tau1Y, m1) (1 - cosh a1(L+z)/cosh a1 L) sin(pi phi)/pi,  z<0
        Phi ~ e P(chi2, tau2Y, m2) (1 - cosh a2(L-z)/cosh a2 L) sin(pi phi)/pi,  z>0

    Both branches vanish at z = 0 -- cosh(a L)/cosh(a L) = 1 -- so Phi is
    continuous across the interface, and the interface is the streamline
    Phi = 0.  That is the consistency check the two expressions have to pass,
    and they do.

    Returned in PF04's convention.  Multiply by 2 for BF25 (CONV-01).
    """
    (c1, d1, t1, m1), (c2, d2, t2, m2) = fluids
    a1, a2 = alpha(c1, d1, t1), alpha(c2, d2, t2)
    phi = np.asarray(phi, dtype=float)
    z = np.asarray(z, dtype=float)
    shape = np.sin(np.pi * phi) / np.pi
    below = e * P(c1, t1, m1) * (1 - np.cosh(a1 * (L + z)) / np.cosh(a1 * L)) * shape
    above = e * P(c2, t2, m2) * (1 - np.cosh(a2 * (L - z)) / np.cosh(a2 * L)) * shape
    return np.where(z < 0, below, above)
