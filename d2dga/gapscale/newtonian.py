"""
Analytic gap-scale closures for a Newtonian fluid pair.

Two roles:
  1. the production path for Newtonian runs -- no AL solve needed, BF25 gives
     these in closed form;
  2. the verification target for the numerical AL solver (test M2-T3).

BF25 Section 2.2, with eta1 = m^(1/2), eta2 = m^(-1/2):

    I1(c, m) = [ m^(1/2) c^3 + m^(-1/2)(1 - c^3) ] / 3               (2.24)
    I2(c, m) = [ 2 m^(1/2) c^3 (1-c) + m^(-1/2) c (1-c)^2 (1+2c) ]/6 (2.25)
    q0(c, m) = c [ m c^2 + 1.5 (1 - c^2) ] / [ m c^3 + 1 - c^3 ]     (2.26)

(2.24)-(2.26) are reproduced exactly by integrating BF25 (2.14), (2.15),
(2.22) -- verified symbolically in tests/test_m2_closures.py.


I3 -- CONV-03, and a correction to BF25 (2.27)
----------------------------------------------
BF25 (2.27) as printed reads

    I3 = c^2 (1-c)^3 [ 4 m c + 3 (1 - c^2) ] / ( 12 [ m c^3 + 1 - c^3 ] )

This does not follow from BF25's own definition (2.23), and it is wrong three
times over -- in the numerator, in the denominator, and in the SIGN.
Integrating (2.23) with eta1 = m^(1/2), eta2 = m^(-1/2) gives, exactly,

    I3 = - c^2 (1-c)^3 [ 4 m c + 3 (1 - c) ]
           / ( 12 m^(1/2) [ m c^3 + 1 - c^3 ] )                  <-- MASTER

Differences from the printed (2.27): the numerator is 3(1-c), not 3(1-c^2);
the denominator carries an m^(1/2) that is missing in print (and (2.24), (2.25)
DO carry their m^(+/-1/2), so the omission is inconsistent with the paper's own
neighbouring equations); and the whole thing is negative.  The derivation is
asserted symbolically, for general m and c, in tests/test_m2_closures.py.


The sign is not cosmetic -- it is the model's central physical statement
-------------------------------------------------------------------------
I3 enters the transport flux as `q0 + b I3` (BF25 2.21, 3.7), so with I3 < 0 a
FAVOURABLE density difference (b > 0, heavy fluid displacing upwards) subtracts
from the dispersive flux, and a large enough b drives the flux function below
the chord from (0,0) to (1,1).  For the Riemann problem c = 1 below, c = 0
above, that turns the front from a spreading fan into a single shock travelling
at exactly the mean speed.  Three independent confirmations that this, and not
the opposite sign, is BF25's:

  1. (2.23) integrated symbolically.  Exact, for all m and c -- and (2.14) and
     (2.15) integrate to the printed (2.24) and (2.25) in the same computation,
     so the layer convention and the algebra are pinned by the same check.
  2. BF25 Section 3.2 in words: "the shock velocity approaches the mean
     velocity (w_f -> 1), and the front becomes more stable when b >> 1 ...
     There is less dispersion for larger b."  With the opposite sign, larger b
     would give a flux hump many times the mean flux and a violent
     counter-current exchange -- the exact opposite.
  3. ZF22 Table 3 case 1, the one case with b < 0 (rho_2 < rho_1): it has by
     far the worst displacement of the ten (t_br = 0.44, eta_E = 0.66).  With
     I3 < 0 a negative b flips the sign of the buoyancy flux and produces
     exactly that instability.

BF25 Section 3.2 also states "the effect of the buoyancy number and the
viscosity ratio on q0'(c) + b I3'(c) is negligible at c = 0, and it maintains a
constant value of 1.5".  MASTER satisfies this identically -- I3 ~ c^2 near
zero, so I3'(0) = 0 for every b and m -- which is a check on the numerator's
leading power rather than on the sign.

Conversions to the other papers.  The other papers print the MAGNITUDE with
their own sign conventions; these helpers reproduce what each paper prints, so
the relation to MASTER carries the minus:

    BCF25 (23) form  =  -MASTER * m^(1/2)
    ZF22  (4.26)     =  -MASTER * 6 / m^(1/2)

ZF22 writes the same physics as `+ (Delta_rho H^3 / 6 eta2) I3` with
Delta_rho < 0 in the favourable case, so the two sign conventions cancel and
the physics -- favourable buoyancy suppresses forward dispersion -- agrees.
"""

from __future__ import annotations

import numpy as np


def _denominator(c, m):
    # Grouped as m c^3 + (1 - c^3), NOT (m c^3 + 1) - c^3.  The latter loses the
    # exact cancellation at c = 1 -- 0.3 + 1.0 - 1.0 = 0.30000000000000004 --
    # and q0(1) then comes out 1 - 2e-16 instead of exactly 1.  M2-T4 requires
    # q0(1) == 1 exactly, and the transport scheme's maximum principle leans on
    # it, so the grouping is load-bearing rather than cosmetic.
    return m * c ** 3 + (1.0 - c ** 3)


def script_I1(c, m):
    """Mean mobility.  BF25 (2.24)."""
    c = np.asarray(c, dtype=float)
    return (np.sqrt(m) * c ** 3 + (1.0 - c ** 3) / np.sqrt(m)) / 3.0


def script_I2(c, m):
    """Buoyant mobility.  BF25 (2.25).  Vanishes at c = 0 and c = 1, because
    buoyancy arises only when two fluids are present in the gap."""
    c = np.asarray(c, dtype=float)
    return (2.0 * np.sqrt(m) * c ** 3 * (1.0 - c)
            + c * (1.0 - c) ** 2 * (1.0 + 2.0 * c) / np.sqrt(m)) / 6.0


def q0(c, m):
    """Isotropic flux.  BF25 (2.26).  q0(0) = 0 and q0(1) = 1 exactly."""
    c = np.asarray(c, dtype=float)
    return c * (m * c ** 2 + 1.5 * (1.0 - c ** 2)) / _denominator(c, m)


def script_I3(c, m):
    """
    Buoyant flux distribution -- MASTER form, see module docstring.

    NEGATIVE on 0 < c < 1, with I3(0) = I3(1) = 0.  The sign is BF25 (2.23)'s
    and is what makes favourable buoyancy suppress dispersion rather than
    drive a counter-current exchange.
    """
    c = np.asarray(c, dtype=float)
    return -(c ** 2 * (1.0 - c) ** 3 * (4.0 * m * c + 3.0 * (1.0 - c))
             / (12.0 * np.sqrt(m) * _denominator(c, m)))


# -- derivatives, needed for LLF wavespeeds (BCF25 32-35) and for the
#    kinematic-wave front speed (BF25 3.8) -------------------------------
#
# These are ANALYTIC, not finite-differenced.  The LLF monotonicity argument
# needs the wavespeed to bound |dq0/dc| from ABOVE at the neighbouring cell
# values; a finite difference can fall below the true slope by O(h^2) and
# would turn the discrete maximum principle into an approximate one.  Both
# functions are rational in c, so exact differentiation is cheap.
#
#   q0 = N / D,   N = (m - 3/2) c^3 + (3/2) c ,       D = (m - 1) c^3 + 1
#   I3 = P / (12 sqrt(m) D),
#                 P = c^2 (1-c)^3 [ (4m - 3) c + 3 ]

def dq0_dc(c, m):
    """d q0 / d c.  Analytic; equals 3/2 at c = 0 for m = 1, the Poiseuille
    centreline speed that sets the leading edge of a dispersive front."""
    c = np.asarray(c, dtype=float)
    D = _denominator(c, m)
    N = (m - 1.5) * c ** 3 + 1.5 * c
    dN = 3.0 * (m - 1.5) * c ** 2 + 1.5
    dD = 3.0 * (m - 1.0) * c ** 2
    return (dN * D - N * dD) / D ** 2


def dscript_I3_dc(c, m):
    """d I3 / d c.  Analytic, MASTER form."""
    c = np.asarray(c, dtype=float)
    D = _denominator(c, m)
    dD = 3.0 * (m - 1.0) * c ** 2
    lin = (4.0 * m - 3.0) * c + 3.0
    P = c ** 2 * (1.0 - c) ** 3 * lin
    dP = (2.0 * c * (1.0 - c) ** 3 * lin
          - 3.0 * c ** 2 * (1.0 - c) ** 2 * lin
          + c ** 2 * (1.0 - c) ** 3 * (4.0 * m - 3.0))
    return -(dP * D - P * dD) / (12.0 * np.sqrt(m) * D ** 2)


# -- other papers' forms, kept ONLY so the conversions can be tested -------
def script_I3_bcf25_form(c, m):
    """BCF25 (23) as printed.  Equals -MASTER * sqrt(m).  Not used in the
    solver."""
    c = np.asarray(c, dtype=float)
    return (c ** 2 * (1.0 - c) ** 3 * (4.0 * m * c + 3.0 * (1.0 - c))
            / (12.0 * _denominator(c, m)))


def script_I3_zf22_form(c, m):
    """ZF22 (4.26) as printed.  Equals -MASTER * 6 / sqrt(m).  Not used in the
    solver."""
    c = np.asarray(c, dtype=float)
    return (c ** 2 * (1.0 - c) ** 3 * (4.0 * m * c + 3.0 * (1.0 - c))
            / (2.0 * m * _denominator(c, m)))


def script_I3_bf25_as_printed(c, m):
    """BF25 (2.27) verbatim.  Retained to document the discrepancy; NOT used."""
    c = np.asarray(c, dtype=float)
    return (c ** 2 * (1.0 - c) ** 3 * (4.0 * m * c + 3.0 * (1.0 - c ** 2))
            / (12.0 * _denominator(c, m)))


def all_closures(c, m):
    """(I1, I2, q0, I3) at concentration c and viscosity ratio m."""
    return script_I1(c, m), script_I2(c, m), q0(c, m), script_I3(c, m)
