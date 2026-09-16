"""M0 test gates -- geometry."""

import numpy as np
import pytest

from d2dga.config import Config, GridConfig, WellConfig
from d2dga.geometry import (  # noqa: F401
    max_symmetric_amplitude_m, synthetic_wall_feasible,
    CaliperLogWall, ConstantEccentricity, ConstantOffsetEccentricity,
    Geometry, GeometryError, SinusoidalWall, UniformWall, WashoutWall,
    build_geometry,
)

CFG = Config()


def _uniform_geometry(e=0.3, n_phi=64, n_xi=128, diameter=None):
    well = CFG.well
    d = diameter if diameter is not None else well.gauge_hole_diameter_m
    return Geometry(well, UniformWall(d), ConstantEccentricity(e),
                    GridConfig(n_phi=n_phi, n_xi=n_xi))


# ---------------------------------------------------------------- M0-T1 ---
def test_m0_t1_uniform_wall_reduces_to_published_form():
    """UniformWall must reproduce H(phi) = 1 + e cos(pi phi) to machine
    precision -- the specialisation used by ZF22/ZF23/BF25/BCF25.  This proves
    the geometry module adds nothing spurious in the uniform limit."""
    for e in (0.0, 0.2, 0.5, 0.8):
        g = _uniform_geometry(e=e)
        phi = g.grid.phi_centres
        xi = g.grid.xi_centres
        H = g.H(phi, xi)
        expected = (1.0 + e * np.cos(np.pi * phi))[:, None] * np.ones(xi.size)[None, :]
        assert np.allclose(H, expected, rtol=0, atol=1e-14), f"e={e}"
        # r_a == 1 and delta == delta_star identically in the uniform limit
        assert np.allclose(g.r_a(xi), 1.0, atol=1e-14)
        assert np.allclose(g.delta(xi), g.delta_star, atol=1e-14)


def test_m0_t1_uniform_wall_H_is_xi_independent():
    g = _uniform_geometry(e=0.4)
    H = g.H(g.grid.phi_centres, g.grid.xi_centres)
    assert np.allclose(H, H[:, :1], atol=1e-14)


# ---------------------------------------------------------------- M0-T2 ---
@pytest.mark.parametrize("wall_name", ["uniform", "sinusoidal", "washout"])
def test_m0_t2_annulus_volume(wall_name):
    """Volume from integrating H*r_a must match pi(r_o^2 - r_i^2) to < 0.1%."""
    well = CFG.well
    walls = {
        "uniform": UniformWall(well.gauge_hole_diameter_m),
        "sinusoidal": SinusoidalWall(well.gauge_hole_diameter_m,
                                     CFG.synthetic_wall.amplitude_m,
                                     CFG.synthetic_wall.wavelength_m),
        "washout": WashoutWall(well.gauge_hole_diameter_m,
                               well.washout_max_diameter_m,
                               centre_xi_hat_m=185.0, sigma_m=5.0),
    }
    d_gauge = 0.5 * (0.5 * well.gauge_hole_diameter_m - well.casing_outer_radius_m)
    g = Geometry(well, walls[wall_name],
                 ConstantOffsetEccentricity(0.3, d_gauge), GridConfig())
    direct = g.annulus_volume_direct()
    modelled = g.annulus_volume_from_H()
    rel = abs(modelled - direct) / direct
    assert rel < 1e-3, f"{wall_name}: {rel:.3e} (direct {direct:.6f}, H {modelled:.6f})"


# ---------------------------------------------------------------- M0-T3 ---
@pytest.mark.parametrize("mode", ["symmetric", "enlargement"])
def test_m0_t3_positivity_and_finiteness(mode):
    well = CFG.well
    wall = SinusoidalWall(well.gauge_hole_diameter_m, CFG.synthetic_wall.amplitude_m,
                          CFG.synthetic_wall.wavelength_m, mode=mode)
    d_gauge = 0.5 * (0.5 * well.gauge_hole_diameter_m - well.casing_outer_radius_m)
    g = Geometry(well, wall, ConstantOffsetEccentricity(0.3, d_gauge), GridConfig())
    phi = np.linspace(0, 1, 201)
    xi = np.linspace(0, g.Z, 2001)
    H = g.H(phi, xi)
    assert np.all(np.isfinite(H))
    assert np.all(H > 0), f"min H = {H.min()}"
    assert np.all(g.e(xi) < 1.0)
    assert np.all(np.isfinite(g.narrow_gap_parameter(xi)))


def test_m0_t3_impossible_geometry_raises():
    """Guards must fire rather than silently producing a negative gap."""
    well = CFG.well
    with pytest.raises(GeometryError, match="non-positive"):
        Geometry(well, UniformWall(0.9 * well.casing_od_m),
                 ConstantEccentricity(0.3), GridConfig())
    with pytest.raises(GeometryError, match="contacts"):
        Geometry(well, UniformWall(well.gauge_hole_diameter_m),
                 ConstantEccentricity(1.05), GridConfig())


def test_m0_t3_caliper_rejects_nan():
    with pytest.raises(GeometryError, match="non-finite"):
        CaliperLogWall.from_arrays([0.0, 1.0, 2.0], [0.26, np.nan, 0.26])


# ---------------------------------------------------------------- M0-T4 ---
def test_m0_t4_div_a_f_vanishes_for_vertical_well():
    """At beta = 0 the term BF25 (2.18) drops is EXACTLY zero, not merely
    small -- and it stays zero under the most violent caliper variation."""
    well = CFG.well
    # A = 3 in in ONE-SIDED mode: the most violent feasible wall for this
    # casing.  The symmetric form at this amplitude pinches the hole to
    # 7.4 in against 7 in casing and drives e = 2.55 -- see
    # test_m0_synthetic_wall_feasibility_boundary.
    wall = SinusoidalWall(well.gauge_hole_diameter_m, 3.0 * 0.0254, 5.0,
                          mode="enlargement")
    d_gauge = 0.5 * (0.5 * well.gauge_hole_diameter_m - well.casing_outer_radius_m)
    g = Geometry(well, wall, ConstantOffsetEccentricity(0.3, d_gauge), GridConfig())
    assert g.well.inclination_rad == 0.0
    div = g.div_a_f(g.grid.phi_centres, g.grid.xi_centres)
    f_phi, f_xi = g.buoyancy_direction(g.grid.phi_centres, g.grid.xi_centres)
    scale = np.abs(np.hypot(f_phi, f_xi)).max()
    assert np.max(np.abs(div)) < 1e-12 * max(scale, 1.0)


def test_m0_t4_div_a_f_is_detected_when_inclined():
    """The test above must not be vacuous: with beta != 0 AND a varying wall
    the divergence is genuinely non-zero, and the machinery sees it."""
    well = WellConfig(inclination_rad=np.pi / 4)
    wall = SinusoidalWall(well.gauge_hole_diameter_m, 3.0 * 0.0254, 5.0,
                          mode="enlargement")
    d_gauge = 0.5 * (0.5 * well.gauge_hole_diameter_m - well.casing_outer_radius_m)
    g = Geometry(well, wall, ConstantOffsetEccentricity(0.3, d_gauge), GridConfig())
    div = g.div_a_f(g.grid.phi_centres, g.grid.xi_centres)
    # closed form, docs/derivation.md Finding 1, with dbeta/dxi = 0
    h = 1e-6
    dra = (g.r_a(g.grid.xi_centres + h) - g.r_a(g.grid.xi_centres - h)) / (2 * h)
    closed = (np.sin(well.inclination_rad) * dra)[None, :] * \
             np.sin(np.pi * g.grid.phi_centres)[:, None]
    assert np.max(np.abs(div)) > 1e-6
    assert np.allclose(div, closed, rtol=1e-5, atol=1e-10)


# ---------------------------------------------------------------- M0-T6 ---
def test_m0_t6_sinusoidal_smooth_positive_and_bounded():
    g = build_geometry(CFG)
    xi = np.linspace(0, g.Z, 4001)
    H = g.H(np.linspace(0, 1, 101), xi)
    assert np.all(H > 0)
    # C1: no jumps in the wall gradient
    grad = g.wall_gradient(xi)
    assert np.all(np.isfinite(grad))
    assert np.max(np.abs(np.diff(grad))) < 1e-2
    dp = g.narrow_gap_parameter(xi)
    assert dp.max() < 1.0 / np.pi


def test_m0_t6_wall_swap_touches_only_geometry():
    """Swapping the wall implementation must change nothing but the profile
    object: same class, same interface, same downstream field shapes."""
    well = CFG.well
    d_gauge = 0.5 * (0.5 * well.gauge_hole_diameter_m - well.casing_outer_radius_m)
    ecc = ConstantOffsetEccentricity(0.3, d_gauge)
    walls = [
        UniformWall(well.gauge_hole_diameter_m),
        SinusoidalWall(well.gauge_hole_diameter_m, 0.0381, 20.0),
        WashoutWall(well.gauge_hole_diameter_m, well.washout_max_diameter_m, 185.0, 5.0),
        CaliperLogWall.from_arrays(np.linspace(0, 196, 50),
                                   np.full(50, well.gauge_hole_diameter_m)),
    ]
    shapes = set()
    for w in walls:
        g = Geometry(well, w, ecc, GridConfig())
        H = g.H(g.grid.phi_centres, g.grid.xi_centres)
        shapes.add(H.shape)
        assert np.all(np.isfinite(H)) and np.all(H > 0)
    assert len(shapes) == 1


def test_m0_t6_caliper_log_reproduces_uniform_wall():
    """A caliper log at constant gauge must be indistinguishable from
    UniformWall -- the interface contract the real log will rely on."""
    well = CFG.well
    xi_hat = np.linspace(0, well.open_hole_length_m, 197)
    cal = CaliperLogWall.from_arrays(xi_hat, np.full_like(xi_hat,
                                                          well.gauge_hole_diameter_m))
    uni = UniformWall(well.gauge_hole_diameter_m)
    probe = np.linspace(0, well.open_hole_length_m, 1000)
    assert np.allclose(cal.outer_radius(probe), uni.outer_radius(probe), atol=1e-15)


def test_m0_caliper_parsers_are_explicit_stubs():
    for fn in (CaliperLogWall.from_las, CaliperLogWall.from_csv,
               CaliperLogWall.from_two_column):
        with pytest.raises(NotImplementedError, match="stubbed"):
            fn("nonexistent.file")


# ---------------------------------------------------- coordinate mapping ---
def test_coordinate_round_trip():
    g = build_geometry(CFG)
    xi = np.linspace(0, g.Z, 101)
    assert np.allclose(g.xi_of(g.xi_hat(xi)), xi, atol=1e-12)
    assert np.isclose(g.depth(0.0), CFG.well.total_depth_m)
    assert np.isclose(g.depth(g.Z), CFG.well.casing_shoe_m)


def test_eccentricity_falls_inside_washout():
    """Physical requirement from build spec Section 4: a constant casing offset
    means e must FALL where the hole enlarges."""
    well = CFG.well
    wall = WashoutWall(well.gauge_hole_diameter_m, well.washout_max_diameter_m,
                       centre_xi_hat_m=185.0, sigma_m=5.0)
    d_gauge = 0.5 * (0.5 * well.gauge_hole_diameter_m - well.casing_outer_radius_m)
    g = Geometry(well, wall, ConstantOffsetEccentricity(0.3, d_gauge), GridConfig())
    xi_peak = g.xi_of(185.0)
    xi_far = g.xi_of(20.0)
    assert g.e(np.array([xi_peak]))[0] < g.e(np.array([xi_far]))[0]
    assert np.isclose(g.e(np.array([xi_far]))[0], 0.3, rtol=1e-6)


# ------------------------------------------------- synthetic-wall feasibility ---
def test_m0_synthetic_wall_feasibility_boundary():
    """
    The build spec's A-sweep {0.5, 1.5, 3.0, 6.0} in is only half-admissible
    with 7 in casing in the literal symmetric form.  Documented, not dodged.

    A = 3.0 in  -> hole pinches to 7.4 in, d_hat = 2.54 mm against a fixed
                   12.95 mm offset, so e = 2.55 >= 1.
    A = 6.0 in  -> hole pinches to 4.4 in, smaller than the 7 in casing.

    In 'enlargement' mode the wall never drops below gauge and all four are fine.
    """
    cfg = Config()
    IN = 0.0254
    expected = {0.5: True, 1.5: True, 3.0: False, 6.0: False}
    for A, ok_expected in expected.items():
        ok, _ = synthetic_wall_feasible(cfg, A * IN, mode="symmetric")
        assert ok is ok_expected, f"A={A} in"
        ok_enl, _ = synthetic_wall_feasible(cfg, A * IN, mode="enlargement")
        assert ok_enl is True

    a_max = max_symmetric_amplitude_m(cfg)
    assert 2.3 * IN < a_max < 2.4 * IN
    ok_below, _ = synthetic_wall_feasible(cfg, a_max * 0.99, mode="symmetric")
    ok_above, _ = synthetic_wall_feasible(cfg, a_max * 1.01, mode="symmetric")
    assert ok_below and not ok_above
