"""M7 test gates -- post-processing metrics (ZF22 Section 6, ZF23 Section 3.2)."""

import numpy as np
import pytest

import tests.benchmarks.kinematic_wave as kw
from d2dga.config import GridConfig, WellConfig
from d2dga.elliptic import NewtonianClosures
from d2dga.gapscale import newtonian as nwt
from d2dga.geometry import ConstantEccentricity, Geometry, UniformWall
from d2dga.postprocess import (axial_profile, dispersion_metrics,
                               displacement_efficiency, front_speed_profile,
                               narrow_side_profile, residual_fraction,
                               zf23_metrics)
from d2dga.scaling import zf22_buoyancy_to_bf25
from d2dga.simulation import Simulation

SHORT = dict(casing_shoe_m=387.2, total_depth_m=390.0)

# ZF23 Table 3, D2DGA columns.  Both cases are CONCENTRIC (e = 0), so the flow
# is planar and BF25 Section 3.1's 1-D theory applies exactly.
ZF23_TABLE3 = {
    "EXP48": dict(e=0.0, b_zf=30.0, m=0.8, Re=144, sigma_plus=0.1369,
                  sigma_minus=0.5493),
    "EXP96": dict(e=0.0, b_zf=750.0, m=0.2, Re=29, sigma_plus=0.0368,
                  sigma_minus=0.2337),
}


def build(n_phi=8, n_xi=300, e=0.0, m=1.0, b=0.0, cfl=0.5):
    well = WellConfig(inclination_rad=0.0, **SHORT)
    geo = Geometry(well, UniformWall(well.gauge_hole_diameter_m),
                   ConstantEccentricity(e), GridConfig(n_phi, n_xi))
    sim = Simulation(geo, NewtonianClosures(m), froude=1.0, delta_rho=-b,
                     cfl=cfl)
    return geo, sim


# =========================================================================
# fields
# =========================================================================
def test_efficiency_is_volume_weighted_not_cell_counted():
    """
    At e = 0.8 the gap varies by a factor of 9 across phi, so counting cells
    instead of weighting by H r_a would over-report the displacement by
    treating the narrow side -- where mud actually survives -- as equal to the
    wide side.  Build a field displaced only on the wide half and check the
    efficiency is the volume fraction, not 0.5.
    """
    geo, _ = build(n_phi=200, n_xi=10, e=0.8)
    c = np.zeros((200, 10))
    c[:100, :] = 1.0                      # wide half only
    g = geo.grid
    H = geo.H(g.phi_centres, g.xi_centres)
    exact = float(np.sum(H[:100]) / np.sum(H))
    assert displacement_efficiency(geo, c) == pytest.approx(exact, rel=1e-12)
    assert exact > 0.6                    # the wide half holds most of the volume
    assert abs(exact - 0.5) > 0.1         # and a cell count would say 0.5


def test_residual_and_narrow_side_diagnostics():
    geo, _ = build(n_phi=20, n_xi=40, e=0.6)
    c = np.ones((20, 40))
    c[-3:, :10] = 0.0                     # a narrow-side mud channel at the bottom
    assert narrow_side_profile(geo, c)[:10].max() == 0.0
    assert narrow_side_profile(geo, c)[10:].min() == 1.0
    f = residual_fraction(geo, c)
    assert 0.0 < f < 0.05                 # narrow side, so a small VOLUME
    assert displacement_efficiency(geo, c) == pytest.approx(1.0 - f, rel=1e-12)


# =========================================================================
# M7-T1 -- the front-speed extraction
# =========================================================================
def test_m7_t1_front_speed_extraction_recovers_the_analytic_wave():
    """
    Gate M7-T1, first half.  The extraction must recover a wave speed profile
    it can be checked against exactly.

    For m = 1, b = 0 the Riemann solution is a rarefaction with
    w_f(c_bar) = q0'(c_bar) = 1.5 (1 - c_bar^2), the profile ZF23's metrics are
    computed from.  Run the full 2-D solver concentrically -- which makes the
    flow planar -- and invert the axial profile.
    """
    geo, sim = build(n_phi=8, n_xi=600, m=1.0, b=0.0)
    Z = geo.grid.Z
    res = sim.run(t_end=0.55 * Z, record_every=100_000)
    t = res.times[-1]

    levels, w_f = front_speed_profile(geo, res.concentration, t, n_bins=50)
    exact = nwt.dq0_dc(levels, 1.0)
    # the first-order scheme smears the two ends of the fan, so compare on the
    # interior where the similarity solution is resolved
    core = (levels > 0.15) & (levels < 0.85)
    assert np.max(np.abs(w_f[core] - exact[core])) < 0.06
    assert w_f[0] > w_f[-1]               # monotone decreasing, as (3.8) requires


def test_m7_t1_metrics_reproduce_zf23_table3():
    """
    Gate M7-T1, second half.  ZF23 Table 3's D2DGA columns, for its two
    concentric cases.

    Computed from BF25 Section 3.1's 1-D theory rather than from a 2-D run:
    both cases have e = 0, so the flow IS planar and the envelope gives the
    converged w_f(c_bar) exactly, with no mesh smearing of the front to unpick.
    The 2-D solver is checked against the same theory by the test above.

    sigma_w+r agrees to 11% on EXP 48.  On EXP 96 we get 0.021 against ZF23's
    0.0368 -- the difference is the dispersive SPIKE, which at b_BF25 = 335 is
    only 0.0027 wide in concentration; ZF23's own caption for figure 7 says "a
    finer mesh would better represent the spike", so their number contains
    some of a feature that the exact theory puts almost entirely inside one
    bin.  Both agree on what matters for the ZF23 (3.6) classification: EXP 48
    dispersive, EXP 96 not.
    """
    got = {}
    for name, ref in ZF23_TABLE3.items():
        b = zf22_buoyancy_to_bf25(ref["b_zf"], ref["m"])
        c = np.linspace(0.0, 1.0, 200001)
        f = kw.flux(c, ref["m"], b)
        hull = kw.upper_concave_envelope(c, f)
        g = np.interp(c, c[hull], f[hull])
        speed = np.gradient(g, c)
        edges = np.linspace(0.0, 1.0, 101)
        levels = 0.5 * (edges[1:] + edges[:-1])
        got[name] = dispersion_metrics(levels, np.interp(levels, c, speed))

    assert got["EXP48"].sigma_plus == pytest.approx(
        ZF23_TABLE3["EXP48"]["sigma_plus"], rel=0.15)
    assert got["EXP48"].sigma_minus == pytest.approx(
        ZF23_TABLE3["EXP48"]["sigma_minus"], rel=0.15)
    # EXP 96: the spike, see the docstring.  Same order of magnitude and the
    # same side of every ZF23 threshold.
    assert 0.5 * ZF23_TABLE3["EXP96"]["sigma_plus"] < got["EXP96"].sigma_plus \
        < ZF23_TABLE3["EXP96"]["sigma_plus"]

    # ZF23 (3.6a,b): dispersive iff sigma_w+r > 0.08 AND |w_r+| > 0.05
    assert got["EXP48"].is_dispersive
    assert not got["EXP96"].is_dispersive
    assert got["EXP48"].sigma_plus > 0.08 and got["EXP48"].area_plus > 0.05
    assert got["EXP96"].sigma_plus < 0.08


def test_m7_sigma_minus_is_bin_dependent_and_the_area_measure_is_not():
    """
    ZF23 computes its standard deviations by "segmenting c_bar and summing over
    a finite number of discrete c_bar values" and does not say how many.  For
    the NEGATIVE part that matters: w_f -> 0 as c_bar -> 1 for every m, because
    q0'(1) = 0, so w_r -> -1 there and the last bins dominate the sum.

    Measured over a 32x range of bin counts (25 to 800), sigma_w-r moves by 14%
    and does so NON-MONOTONICALLY -- the successive differences change sign, so
    it is not a sequence converging to a limit.  The practical consequence is
    that no sigma_w-r comparison against ZF23 can be read as tighter than about
    15%, which is why BENCH-08 states its agreement at that level and no better.

    ZF23 introduced |w_r+| and |w_r-| (3.4)-(3.5) precisely because "the
    standard deviations are influenced by large deviations".  Those DO settle,
    and the test asserts both halves so the distinction is not lost later.
    """
    b = zf22_buoyancy_to_bf25(750.0, 0.2)
    c = np.linspace(0.0, 1.0, 200001)
    f = kw.flux(c, 0.2, b)
    hull = kw.upper_concave_envelope(c, f)
    speed = np.gradient(np.interp(c, c[hull], f[hull]), c)

    sig_minus, area_minus = [], []
    for n in (25, 50, 100, 200, 400, 800):
        edges = np.linspace(0.0, 1.0, n + 1)
        levels = 0.5 * (edges[1:] + edges[:-1])
        m = dispersion_metrics(levels, np.interp(levels, c, speed))
        sig_minus.append(m.sigma_minus)
        area_minus.append(m.area_minus)

    spread = (max(sig_minus) - min(sig_minus)) / max(sig_minus)
    assert spread > 0.1, sig_minus
    d_sig = np.diff(sig_minus)
    assert np.any(d_sig < 0) and np.any(d_sig > 0), d_sig   # not a limit

    d_area = np.abs(np.diff(area_minus))
    assert np.all(d_area[-3:] == np.sort(d_area[-3:])[::-1]), area_minus
    assert d_area[-1] < 0.5 * d_area[0], area_minus


def test_m7_axial_profile_is_volume_weighted():
    """The azimuthal average must weight by H r_a, not by cell count -- same
    reason as the efficiency, and it feeds the front-speed extraction."""
    geo, _ = build(n_phi=100, n_xi=5, e=0.7)
    c = np.zeros((100, 5))
    c[:50, :] = 1.0
    g = geo.grid
    H = geo.H(g.phi_centres, g.xi_centres)
    exact = np.sum(H[:50], axis=0) / np.sum(H, axis=0)
    assert np.allclose(axial_profile(geo, c), exact, rtol=1e-12)
