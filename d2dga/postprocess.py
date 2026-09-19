"""
M7 -- post-processing: the metrics the displacement is actually judged by.

Two families, from the two papers that define them.

ZF22 Section 6 -- how much and how soon
    t_br   = t_hat_br w0_hat / L_annu_hat, the dimensionless breakthrough time.
             In BF25 scaling the piston time is exactly Z, so t_br = t / Z.
    eta_E  = displaced fraction of the annulus volume at t = 1.2 Z.
Both are on RunResult (`breakthrough_at`, `efficiency_at`), because they need
the time history rather than a single field.

ZF23 Section 3.2 -- what SHAPE the front has
    w_f(c_bar)   front speed of each concentration level, from the converged
                 profile; for a planar flow this is BF25 (3.8).
    w_r          = w_f - 1, split into its positive and negative parts:
                 w_r+ is dispersion running AHEAD of the mean front, w_r- is
                 fluid left BEHIND it -- residual wall layers and narrow-side
                 mud channels, which is what actually ruins a cement job.
    sigma_w+r, sigma_w-r   ZF23 (3.2), (3.3)
    |w_r+|, |w_r-|         ZF23 (3.4), (3.5)
    classification         ZF23 (3.6): dispersive iff sigma_w+r > 0.08 AND
                           |w_r+| > 0.05.

A note on sigma_w-r, because it is easy to over-read.  ZF23 computes these by
"segmenting c_bar and summing over a finite number of discrete c_bar values",
and the result is NOT independent of that segmentation for the negative part:
w_f -> 0 as c_bar -> 1 (q0'(1) = 0 for every m), so w_r -> -1 there and the
last bins dominate the sum.  Refining the bins therefore moves sigma_w-r and
does not converge.  The bin count is exposed rather than hidden, and the
integral measures |w_r+|, |w_r-| -- which ZF23 introduced precisely because
"the standard deviations are influenced by large deviations" -- do converge.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ===========================================================================
# fields
# ===========================================================================
def cell_volume(geometry):
    """H r_a at cell centres times the cell area: the volume weight."""
    g = geometry.grid
    H = geometry.H(g.phi_centres, g.xi_centres)
    ra = geometry.r_a(g.xi_centres)[None, :]
    return H * ra * g.dphi * g.dxi


def displacement_efficiency(geometry, c):
    """Displaced fraction of the (half) annulus volume."""
    w = cell_volume(geometry)
    return float(np.sum(w * np.asarray(c, dtype=float)) / np.sum(w))


def axial_profile(geometry, c):
    """
    Volume-weighted azimuthal average of c_bar at each xi.

    Volume-weighted, not plain: H varies by a factor of 9 across phi at
    e = 0.8, so a plain mean would count the narrow side -- where the mud
    actually survives -- as heavily as the wide side.
    """
    g = geometry.grid
    H = geometry.H(g.phi_centres, g.xi_centres)
    ra = geometry.r_a(g.xi_centres)[None, :]
    w = H * ra
    return np.sum(w * np.asarray(c, dtype=float), axis=0) / np.sum(w, axis=0)


def narrow_side_profile(geometry, c):
    """c_bar on the narrow side (phi = 1) at each xi -- where a mud channel
    survives if one survives anywhere."""
    return np.asarray(c, dtype=float)[-1, :]


def narrow_side_efficiency(geometry, c, fraction=0.25):
    """Displacement efficiency over the NARROW quarter of the annulus.

    `displacement_efficiency` is a volume-weighted mean over the whole
    annulus, and on an eccentric well it is dominated by the wide side -- the
    part that was never in doubt.  At e = 0.3 the wide half carries most of the
    volume, so a job that leaves a continuous mud channel on the narrow side
    can still report eta_E above 0.99.  K-GEP-1 does exactly that: 0.9946
    overall against a narrow-side minimum of 0.9494.

    So this is the same integral restricted to the narrowest `fraction` of phi.
    It is not a different measure of a different thing; it is the measure that
    answers the question primary cementing is about.  `narrow_side_profile`
    gives the single narrowest column, which is the worst case rather than an
    average -- both are worth having, and they answer different questions.

    phi runs wide (0) to narrow (1), so the narrow cells are at the end.  At
    least one column is always included, so the function is defined on any
    mesh.
    """
    g = geometry.grid
    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"fraction must be in (0, 1], got {fraction!r}")
    k = max(1, int(round(fraction * g.n_phi)))
    w = cell_volume(geometry)[-k:, :]
    return float(np.sum(w * np.asarray(c, dtype=float)[-k:, :]) / np.sum(w))


def residual_fraction(geometry, c, threshold=0.5):
    """Volume fraction of the annulus still holding mostly displaced fluid."""
    w = cell_volume(geometry)
    mask = np.asarray(c, dtype=float) < threshold
    return float(np.sum(w[mask]) / np.sum(w))


# ===========================================================================
# ZF23 Section 3.2
# ===========================================================================
@dataclass(frozen=True)
class DispersionMetrics:
    sigma_plus: float
    sigma_minus: float
    area_plus: float
    area_minus: float
    n_bins: int

    @property
    def is_dispersive(self) -> bool:
        """ZF23 (3.6a,b)."""
        return self.sigma_plus > 0.08 and self.area_plus > 0.05


def front_speed_profile(geometry, c, t, n_bins=100):
    """
    w_f(c_bar) at bin midpoints, from the axial profile at time t.

    The front is treated as a similarity solution, which is what ZF23 means by
    "the converged profiles": each concentration level sits at xi = w_f t, so
    inverting the profile gives w_f directly.  Requires the profile to have
    left the initial condition behind, i.e. a late time.
    """
    if t <= 0:
        raise ValueError("t must be positive")
    prof = axial_profile(geometry, c)
    xi = geometry.grid.xi_centres
    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    levels = 0.5 * (edges[1:] + edges[:-1])
    # prof decreases with xi; np.interp needs an increasing x, so reverse
    xi_of_c = np.interp(levels, prof[::-1], xi[::-1])
    return levels, xi_of_c / t


def dispersion_metrics(levels, w_f):
    """ZF23 (3.2)-(3.5) and the (3.6) classification."""
    wr = np.asarray(w_f, dtype=float) - 1.0
    levels = np.asarray(levels, dtype=float)
    plus, minus = wr[wr > 0], wr[wr < 0]
    sigma_plus = float(np.sqrt(np.mean(plus ** 2))) if plus.size else 0.0
    sigma_minus = float(np.sqrt(np.mean(minus ** 2))) if minus.size else 0.0
    area_plus = float(np.trapezoid(np.clip(wr, 0.0, None), levels))
    area_minus = float(-np.trapezoid(np.clip(wr, None, 0.0), levels))
    return DispersionMetrics(sigma_plus, sigma_minus, area_plus, area_minus,
                             len(levels))


def zf23_metrics(geometry, c, t, n_bins=100):
    """Convenience: profile then metrics."""
    levels, w_f = front_speed_profile(geometry, c, t, n_bins)
    return dispersion_metrics(levels, w_f)


# --------------------------------------------------------------------------
# Discretisation error of a threshold breakthrough time (A-2)
# --------------------------------------------------------------------------
# A threshold front is a LEVEL SET of a front the first-order scheme smears, so
# its error is O(dxi^(1/2)) -- half order -- while eta_E, an integral, is
# O(dxi).  Measured against the exact rarefaction by
# `test_a2_threshold_front_is_half_order_and_integrals_are_first_order` at
# n_xi = 80.
#
# These constants lived in `scripts/kgep1_run.py` AND in the UI, which is two
# copies of a measurement: the day the test is re-run at a finer mesh, one of
# them silently becomes a lie.  They live here, with the other post-processing,
# and both callers ask.
TBR_RELATIVE_ERROR_AT_80 = {0.01: 0.175, 0.1: 0.069, 0.5: 0.033}
TBR_ERROR_REFERENCE_N_XI = 80


def tbr_relative_error(threshold: float, n_xi: int) -> float:
    """Relative error to attach to a breakthrough time at this mesh.

    Half-order scaling from the measured reference, so halving the error costs
    four times the cells.  Raises for a threshold that was never measured
    rather than interpolating between two measurements -- the whole point of
    the number is that it is measured.
    """
    try:
        ref = TBR_RELATIVE_ERROR_AT_80[threshold]
    except KeyError:
        raise KeyError(
            f"no measured t_br error for threshold {threshold!r}; measured "
            f"thresholds are {sorted(TBR_RELATIVE_ERROR_AT_80)}") from None
    return ref * (TBR_ERROR_REFERENCE_N_XI / float(n_xi)) ** 0.5
