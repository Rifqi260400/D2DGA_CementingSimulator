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
