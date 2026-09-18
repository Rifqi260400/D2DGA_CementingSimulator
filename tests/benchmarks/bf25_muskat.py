"""
BF25 Section 3.3 -- Muskat stability analysis of the dispersive front.

This is the only quantitative, closed-form result in any of the six papers that
applies to a HERSCHEL-BULKLEY pair, which is why it is here: the ZF22 and ZF23
benchmarks are Newtonian without exception, so before this file nothing in the
repository tested the M3 closure table against an external analytic statement
(audit finding A-4.3).

BF25's construction.  The base flow is the dispersing planar front (3.7) with
gap-averaged concentration c_bar_0.  The axial pressure gradient through it is
set by c_bar_0, BF25 (3.10):

    dp/dxi = -rho(c_bar_0)/Fr*^2 - [1 + b I2(c_bar_0)] / I1(c_bar_0)

A thin finger of pure fluid 2 advancing in the (phi, xi) plane feels that SAME
pressure gradient, BF25 (3.11), because the finger is thin and the pressure
transmits across it.  Evaluating (3.2) inside the finger -- where c_bar = 1 --
then gives the finger velocity, BF25 (3.12):

    w_finger(c_bar_0) = I1(1)/I1(c_bar_0)
                        + b I1(1) [ I2(c_bar_0)/I1(c_bar_0) + c_bar_0 - 1 ]

using I2(1) = 0.  The first term is a mobility ratio; the second is buoyancy.

Stability, BF25 (3.13):  Delta_w(c_bar_0) = w_finger(c_bar_0) - w_f(c_bar_0),
with w_f the front speed of the dispersing base state.  BF25: "For an unstable
regime, Delta_w(c_0) > 0 for all c_0 in [0,1]" -- the finger then advances
through the front at every concentration.  The three regimes BF25 draws are
therefore

    Delta_w > 0 everywhere            -> unstable
    Delta_w < 0 everywhere            -> stable
    Delta_w changes sign              -> partial penetration

Nothing here is fitted: I1 and I2 come from whatever closure provider is passed
in, so the same code runs on the analytic Newtonian forms and on an M3 table.
"""

from __future__ import annotations

import numpy as np


def finger_velocity(script_I1, script_I2, c0, b):
    """
    BF25 (3.12), verbatim.

    `script_I1`, `script_I2` are callables of c_bar returning the H-free
    closures.  `c0` may be an array.  I2(1) = 0 is used by BF25 in deriving
    (3.12); `check_I2_at_one` asserts it rather than assuming it.
    """
    c0 = np.asarray(c0, dtype=float)
    I1_1 = float(np.asarray(script_I1(1.0)).reshape(()))
    I1_0 = np.asarray(script_I1(c0), dtype=float)
    I2_0 = np.asarray(script_I2(c0), dtype=float)
    return I1_1 / I1_0 + b * I1_1 * (I2_0 / I1_0 + c0 - 1.0)


def check_I2_at_one(script_I2, tol=1e-10) -> float:
    """BF25 uses I2(1) = 0 to reach (3.12).  Returns |I2(1)|, for assertion."""
    return abs(float(np.asarray(script_I2(1.0)).reshape(())))


def front_speed(q0, script_I3, c0, b):
    """
    The dispersing front speed w_f(c_bar_0) of BF25 Section 3.1: the speed of
    the wave carrying concentration c_bar_0 in the entropy solution of the
    planar problem (3.7), i.e. the local slope of the UPPER CONCAVE ENVELOPE of
    f = q0 + b I3 on [0, 1].

    Using the envelope rather than f' is what BF25 means by "w_f(c_bar_0) has
    been computed to include any segments of c_bar_0 over which there is a
    shock front".  Inside a shock the envelope is straight and every
    concentration it spans travels at the one shock speed.
    """
    c0 = np.atleast_1d(np.asarray(c0, dtype=float))
    cg = np.linspace(0.0, 1.0, 200001)
    fg = np.asarray(q0(cg), dtype=float) + b * np.asarray(script_I3(cg),
                                                          dtype=float)
    V = _upper_concave_hull(cg, fg)
    cs, fs = V[:, 0], V[:, 1]
    slope = np.diff(fs) / np.diff(cs)
    idx = np.clip(np.searchsorted(cs, c0, side="right") - 1, 0, len(slope) - 1)
    return slope[idx]


def _upper_concave_hull(c, f):
    hull = []
    for i in range(len(c)):
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (y2 - y1) * (c[i] - x1) > (f[i] - y1) * (x2 - x1):
                break
            hull.pop()
        hull.append((c[i], f[i]))
    return np.array(hull)


def classify(script_I1, script_I2, q0, script_I3, b, n=401, c0_max=1.0):
    """
    BF25 (3.13).  Returns (regime, dw, c0), regime in
    {"unstable", "stable", "partial penetration"}.

    `c0_max = 1.0` is BF25 verbatim -- "for an unstable regime, Delta_w(c_0) > 0
    for all c_0 in [0,1]".  KNOWN PROBLEM WITH THE LITERAL READING, quantified
    in docs/remediation_log.md under Q-3: at c_0 -> 1 the base state has
    w_f -> q0'(1) + b I3'(1) = 0 + 0 = 0 for EVERY pair and every b, while
    w_finger -> I1(1)/I1(1) + b I1(1)[0 + 1 - 1] = 1.  So Delta_w(1-) -> +1
    unconditionally, which makes "stable" unreachable and "unstable"
    unfalsifiable at that end: every fluid pair is classified "partial
    penetration", including the ten ZF22 cases, which cannot be what BF25's
    three-regime figure shows.  The c_0 -> 1 limit is also physically empty --
    c_bar_0 = 1 is an annulus already full of fluid 2, with no front to finger
    into.

    `c0_max < 1` restricts the comparison to concentrations where a front
    actually exists.  The classification is then stable over a wide range of
    c0_max because Delta_w is monotone in the cases tested, and the
    DISCRIMINATING quantity is the sign of Delta_w at SMALL c_0 -- the leading
    edge of the front, which is what a finger must outrun to penetrate.
    `dw_at_leading_edge` returns it directly, and it is the number this
    repository uses.
    """
    c0 = np.linspace(1.0 / n, c0_max - (1.0 - c0_max) / n if c0_max < 1.0
                     else 1.0 - 1.0 / n, n)
    dw = finger_velocity(script_I1, script_I2, c0, b) \
        - front_speed(q0, script_I3, c0, b)
    if np.all(dw > 0.0):
        return "unstable", dw, c0
    if np.all(dw < 0.0):
        return "stable", dw, c0
    return "partial penetration", dw, c0


def dw_at_leading_edge(script_I1, script_I2, q0, script_I3, b, c0=1e-3):
    """
    Delta_w at the LEADING EDGE of the dispersing front (c_bar_0 -> 0).

    This is the half of BF25 (3.13) that is both well posed and discriminating.
    A finger of pure fluid 2 penetrates the front only if it outruns the front's
    own leading wave, so Delta_w(0+) > 0 is the condition for the finger to
    break through rather than be absorbed.

    For a Newtonian pair at b = 0 it reduces to a clean classical statement:
    w_finger(0+) = I1(1)/I1(0) = m and w_f(0+) = q0'(0) = 3/2, so
    Delta_w(0+) = m - 3/2 -- the finger penetrates exactly when m > 3/2, which
    is BF25 Section 3.2's own M_3^min = 1.5 attributed to Lajeunesse et al.
    (1999).  That coincidence is asserted as a test, and it is what pins this
    reading of (3.13) to the paper rather than to preference.
    """
    return float(finger_velocity(script_I1, script_I2, c0, b)
                 - front_speed(q0, script_I3, c0, b)[0])
