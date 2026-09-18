"""
Exact 1-D entropy solution for a concentric annulus.
Exact 1-D entropy solution for a concentric annulus, from the upper concave
envelope of f(c) = q0 + b I3 (Oleinik; BF25 Section 3.1).

Riemann problem u_l = 1 (below), u_r = 0 (above) -> u_l > u_r -> upper concave
envelope g of f on [0,1].  g' is decreasing; the similarity solution is
u(s) = sup{ c : g'(c) >= s }.  Then, in a uniform annulus,

    eta_E(t) = (1/Z) int_0^Z u(xi/t) dxi = (t/Z) int_0^{Z/t} u(s) ds.

Added in remediation: A-1 needed an exact reference to show that the outlet
import converges away, and A-4.1 needs one to measure how much each published
case actually constrains the model.
"""

import numpy as np


def concave_hull(c, f):
    """Upper concave envelope vertices of (c, f), c increasing."""
    hull = []
    for i in range(len(c)):
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            # drop the middle point if it is BELOW the chord x1->current
            if (y2 - y1) * (c[i] - x1) > (f[i] - y1) * (x2 - x1):
                break
            hull.pop()
        hull.append((c[i], f[i]))
    return np.array(hull)


def similarity(c_grid, f_grid):
    """Return (s_nodes, u_nodes) with s increasing and u decreasing."""
    V = concave_hull(c_grid, f_grid)
    cs, fs = V[:, 0], V[:, 1]
    sig = np.diff(fs) / np.diff(cs)             # decreasing
    # walk the hull from c = 1 (smallest slope) down to c = 0 (largest slope)
    s_nodes, u_nodes = [], []
    for k in range(len(sig) - 1, -1, -1):
        s_nodes += [sig[k], sig[k]]
        u_nodes += [cs[k + 1], cs[k]]
    return np.array(s_nodes), np.array(u_nodes)


def u_of_s(s, s_nodes, u_nodes):
    s = np.atleast_1d(np.asarray(s, float))
    out = np.empty_like(s)
    out[s <= s_nodes[0]] = 1.0
    out[s >= s_nodes[-1]] = 0.0
    mid = (s > s_nodes[0]) & (s < s_nodes[-1])
    out[mid] = np.interp(s[mid], s_nodes, u_nodes)
    return out


def exact_eta(t, Z, s_nodes, u_nodes, n=400001):
    """(1/Z) int_0^Z u(xi/t) dxi."""
    xi = np.linspace(0.0, Z, n)
    return float(np.trapezoid(u_of_s(xi / t, s_nodes, u_nodes), xi) / Z)


def exact_tbr(threshold, Z, s_nodes, u_nodes):
    """t at which u(Z/t) first exceeds `threshold`."""
    lead = s_nodes[-1]
    lo, hi = Z / lead, 50.0 * Z
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if float(u_of_s(Z / mid, s_nodes, u_nodes)[0]) > threshold:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)
