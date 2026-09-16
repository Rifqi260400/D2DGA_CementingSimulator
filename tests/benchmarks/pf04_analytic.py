"""
PF04 analytic steady-state benchmarks, re-derived.

The display equations in the supplied PF04 PDF do not extract (the text layer
returns them as scrambled fragments and no page renderer is available), so
rather than transcribe (38) and (40)-(43) by eye we DERIVE the concentric
steady state from the model equations and check the structure against the
fragments that do survive.

Derivation (see tests/test_m4_elliptic.py::test_m4_t2_derivation_matches_pf04):
concentric annulus, sharp interface xi = g(phi) + t travelling at unit speed,
pure fluid on each side.  Far from the interface each fluid is fully developed,
so with BF25 (2.13) and I2 = 0 for a single fluid,

    v_bar = 0  =>  G_phi = 0  =>  dp/dphi = rho_k sin(beta) sin(pi phi) / Fr^2
    w_bar = 1  =>  G_xi  = 1/I1_k  =>  dp/dxi = -1/I1_k - rho_k cos(beta)/Fr^2

Pressure continuity along the interface then gives

    g(phi) = -(1/pi) b sin(beta) cos(pi phi) / ( J + b cos(beta) ),
    J = [1/I1_k]^2_1 = 1/I1_2 - 1/I1_1

which is exactly PF04 (38) once its jump notation [.]^2_1 is recognised, with
the identification 1/I1_k <-> chi_k(1) + tau_k,Y.  Both sides are the modified
pressure gradient that drives unit areal flow through a unit half-gap.

Note b here is the r_a-FREE buoyancy scalar of CONV-04, b = (rho2-rho1)/Fr*^2,
positive when the denser fluid displaces upward -- the same sign as PF04 (39).
"""

import numpy as np


def concentric_steady_interface(phi, b, beta, inv_I1_displaced, inv_I1_displacing):
    """
    PF04 (38).  Interface offset g(phi) about its mean, in scaled axial units.

    At beta = 0 this is identically zero: a vertical concentric displacement has
    no azimuthal buoyancy component, so the steady interface is flat.
    """
    J = inv_I1_displacing - inv_I1_displaced
    denom = J + b * np.cos(beta)
    return -(1.0 / np.pi) * b * np.sin(beta) * np.cos(np.pi * phi) / denom


def concentric_steady_interface_slope(phi, b, beta, inv_I1_displaced,
                                      inv_I1_displacing):
    """dg/dphi, needed for the kinematic condition PF04 (15)."""
    J = inv_I1_displacing - inv_I1_displaced
    denom = J + b * np.cos(beta)
    return b * np.sin(beta) * np.sin(np.pi * phi) / denom


def newtonian_inverse_mobilities(m):
    """
    1/I1 for each pure fluid, D2DGA Newtonian.

    BF25 (2.24) at c = 0 gives I1 = m^(-1/2)/3 (displaced fluid alone) and at
    c = 1 gives I1 = m^(1/2)/3 (displacing fluid alone).
    """
    return 3.0 * np.sqrt(m), 3.0 / np.sqrt(m)


# ---------------------------------------------------------------------------
# PF04's own 2DGA closure, needed to evaluate (38) with PF04's parameters
# ---------------------------------------------------------------------------
def chi_at_unit_flux(kappa, m, tau_y):
    """
    Invert PF04 (8) / B02 (57) at |grad Psi| = 1, H = 1:

        1 = chi^(m+1) [chi + (m+2) tau_Y/(m+1)]
            / ( kappa^m (m+2) (chi + tau_Y)^2 )

    chi + tau_Y is then the modified pressure gradient carrying unit areal flow
    through a unit half-gap -- PF04's counterpart of the D2DGA 1/I1.  Note m
    here is PF04's INVERSE power-law index, m = 1/n.
    """
    from scipy.optimize import brentq
    def residual(c):
        return (c ** (m + 1) * (c + (m + 2) * tau_y / (m + 1))
                / (kappa ** m * (m + 2) * (c + tau_y) ** 2) - 1.0)
    return brentq(residual, 1e-12, 1e4, xtol=1e-14, rtol=1e-15)


# PF04 Figs 3-5 parameter set (stated in their captions)
PF04_FIG_PARAMS = dict(e=0.0, St=0.1, rho1=1.0, rho2=0.9,
                       tau1Y=0.9, tau2Y=0.7, kappa1=0.5, kappa2=0.4,
                       m1=1.0, m2=1.2)


def pf04_figure_jump_and_b():
    """(J, b) for the PF04 Fig. 3-5 parameter set, in PF04's own convention."""
    p = PF04_FIG_PARAMS
    g1 = chi_at_unit_flux(p["kappa1"], p["m1"], p["tau1Y"]) + p["tau1Y"]
    g2 = chi_at_unit_flux(p["kappa2"], p["m2"], p["tau2Y"]) + p["tau2Y"]
    return g2 - g1, (p["rho2"] - p["rho1"]) / p["St"]
