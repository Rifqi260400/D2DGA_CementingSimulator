"""
BF25 Section 3.1 -- kinematic wave analysis of the 1-D flux f(c) = q0 + b I3.

The planar reduction (BF25 3.7) is

    dc/dt + d/dxi [ q0(c) + b I3(c) ] = 0,

a scalar conservation law.  For the displacement Riemann problem -- c = 1 below
(small xi), c = 0 above -- the entropy solution is read off the UPPER CONCAVE
ENVELOPE of f on [0, 1], because the left state exceeds the right one.  Where
the envelope touches f the solution is a rarefaction (BF25's "dispersion");
each straight segment of the envelope is a shock whose speed is the segment's
slope, which is BF25 (3.9), the Rankine-Hugoniot condition.

This is an analytic reference, not part of the solver: it gives the exact
similarity solution c(xi/t) that the 2-D code must reproduce in its planar
limit, and it reproduces BF25's four regimes (dispersion, shock, spike, static
wall layer) from the flux function alone.
"""

import numpy as np

from d2dga.gapscale import newtonian as _newt


def flux(c, m, b):
    """BF25 (3.7)'s flux function for a Newtonian pair."""
    return _newt.q0(c, m) + b * _newt.script_I3(c, m)


def upper_concave_envelope(x, y):
    """
    Upper concave envelope of the sampled function y(x), x increasing.

    Monotone-chain upper hull.  Returns the indices of the hull vertices; the
    envelope is the piecewise-linear interpolant through them, and between two
    consecutive vertices that are NOT adjacent in x the solution is a shock.
    """
    hull = []
    for i in range(len(x)):
        while len(hull) >= 2:
            a, bb = hull[-2], hull[-1]
            # drop the middle point if it lies on or below the chord a->i
            cross = ((x[bb] - x[a]) * (y[i] - y[a])
                     - (y[bb] - y[a]) * (x[i] - x[a]))
            if cross >= 0:
                hull.pop()
            else:
                break
        hull.append(i)
    return np.array(hull)


def riemann_structure(m, b, n=40001):
    """
    Decompose the c = 1 -> c = 0 Riemann solution.

    Returns a dict with
        vertices     : the envelope's vertex concentrations, 0 ... 1
        shocks       : list of (c_lo, c_hi, speed) for each straight segment
                       that spans more than one sample
        leading_speed: the envelope's slope leaving c = 0, i.e. the speed of the
                       fastest part of the front
        main_shock   : the (c_lo, c_hi, speed) of the segment spanning the
                       largest concentration interval, or None
    """
    c = np.linspace(0.0, 1.0, n)
    f = flux(c, m, b)
    hull = upper_concave_envelope(c, f)
    verts = c[hull]

    shocks = []
    for k in range(len(hull) - 1):
        i, j = hull[k], hull[k + 1]
        if j - i > 1:
            speed = (f[j] - f[i]) / (c[j] - c[i])
            shocks.append((c[i], c[j], float(speed)))

    i, j = hull[0], hull[1]
    leading = float((f[j] - f[i]) / (c[j] - c[i]))
    main = max(shocks, key=lambda s: s[1] - s[0]) if shocks else None
    return dict(c=c, f=f, vertices=verts, shocks=shocks,
                leading_speed=leading, main_shock=main)


def similarity_solution(zeta, m, b, n=40001):
    """
    c(xi/t) for the Riemann problem, evaluated at the similarity variable zeta.

    Built from the envelope: the wave speed at concentration c is the envelope
    slope there, which decreases monotonically in c, so inverting it gives c as
    a function of zeta.
    """
    st = riemann_structure(m, b, n)
    c, f = st["c"], st["f"]
    hull = upper_concave_envelope(c, f)
    # envelope value and slope on the hull
    g = np.interp(c, c[hull], f[hull])
    speed = np.gradient(g, c)
    # speed is non-increasing in c; invert
    zeta = np.atleast_1d(np.asarray(zeta, dtype=float))
    out = np.interp(-zeta, -speed[::-1][::-1], c, left=1.0, right=0.0)
    order = np.argsort(-speed)
    out = np.interp(zeta, speed[order][::-1], c[order][::-1], left=1.0, right=0.0)
    return out
