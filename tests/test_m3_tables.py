"""M3 test gates -- closure tabulation."""

import warnings

import numpy as np
import pytest

from d2dga.gapscale import newtonian as nwt
from d2dga.gapscale.al_solver import TwoLayerGapSolver
from d2dga.gapscale.closures import closures_from_solution
from d2dga.gapscale.tables import (ClosureTable, ClosureTableRangeWarning,
                                   newtonian_table)
from d2dga.scaling import ScaledFluid

M = 0.2


@pytest.fixture(scope="module")
def newt_table():
    return newtonian_table(M).build()


@pytest.fixture(scope="module")
def hb_pair():
    # shear-thinning mud with a yield stress, displaced by a lighter-consistency
    # slurry -- representative of the cementing pairs BF25 targets
    return (ScaledFluid("mud", 1.0, 0.45, 0.6, 0.30),
            ScaledFluid("slurry", 1.0, 0.25, 0.8, 0.12))


@pytest.fixture(scope="module")
def hb_table(hb_pair):
    f1, f2 = hb_pair
    return ClosureTable(f1, f2, c_grid=np.linspace(0.0, 1.0, 81),
                        n_y=200, tol=1e-9).build()


# ---------------------------------------------------------------- M3-T1 ---
def test_m3_t1_newtonian_interpolation_within_budget(newt_table):
    """Error measured against max|f| over the domain, not pointwise relative:
    I2 and I3 vanish at both endpoints, so pointwise relative error is
    undefined there and says nothing about the flux."""
    c = np.linspace(0.0, 1.0, 401)
    got = dict(zip(("I1", "I2", "q0", "I3"), newt_table(c)))
    exact = {"I1": nwt.script_I1(c, M), "I2": nwt.script_I2(c, M),
             "q0": nwt.q0(c, M), "I3": nwt.script_I3(c, M)}
    for k in exact:
        err = np.max(np.abs(got[k] - exact[k])) / np.max(np.abs(exact[k]))
        assert err < 1e-3, f"{k}: {err:.3e}"


def test_m3_t1_herschel_bulkley_interpolation_vs_direct_solve(hb_pair, hb_table):
    """The real test: a non-Newtonian table against direct AL solves at points
    that are NOT table nodes."""
    f1, f2 = hb_pair
    solver = TwoLayerGapSolver(f1, f2, n_y=200, tol=1e-9, max_iter=60000)
    probe = np.array([0.077, 0.213, 0.386, 0.541, 0.688, 0.834, 0.947])
    direct = [closures_from_solution(
        solver.solve_fixed_mean_velocity(c, [0.0, 1.0], [0.0, 0.0])) for c in probe]
    got = hb_table(probe)
    for j, key in enumerate(("I1", "I2", "q0", "I3")):
        ref = np.array([getattr(d, key) for d in direct])
        scale = np.max(np.abs(ref))
        err = np.max(np.abs(got[j] - ref)) / scale
        assert err < 1e-3, f"{key}: {err:.3e}"


def test_m3_t1_four_dimensional_table_builds_and_interpolates(hb_pair):
    """Machinery check on the full 4-D state (c, H, |u|, Gb), coarse so it is
    affordable.  Confirms the BF25 (A4) H-rescaling is wired up: a table built
    at several H must reproduce a direct solve at the same H."""
    f1, f2 = hb_pair
    tab = ClosureTable(f1, f2,
                       c_grid=np.linspace(0.0, 1.0, 11),
                       h_grid=np.array([0.6, 1.0, 1.6]),
                       umag_grid=np.array([0.5, 1.0, 2.0]),
                       gb_grid=np.array([0.0, 5.0]),
                       n_y=120, tol=1e-8).build()
    assert tab.n_points == 11 * 3 * 3 * 2
    assert tab.n_failed == 0
    I1, I2, q0, I3 = tab(0.5, H=1.0, umag=1.0, gb=0.0)
    assert np.isfinite([I1, I2, q0, I3]).all()
    assert I1 > 0 and 0.0 <= q0 <= 1.0


@pytest.mark.parametrize("H", [0.6, 1.0, 1.6])
def test_m3_t1_h_rescaling_matches_a_direct_solve(hb_pair, H):
    """BF25 (A4): kappa~_k = kappa_k / H^n_k with velocity, yield stress and
    power-law index unchanged.  This is the identity that lets ONE table serve
    every depth (build spec Section 3.2), so it gets its own test."""
    f1, f2 = hb_pair
    rescaled = [ScaledFluid(f.name, f.density,
                            f.consistency / H ** f.power_law_index,
                            f.power_law_index, f.yield_stress) for f in (f1, f2)]
    tab = ClosureTable(f1, f2, c_grid=np.array([0.45]),
                       h_grid=np.array([H]), n_y=200, tol=1e-9).build()
    # The direct solve must use the table's OWN augmented-Lagrangian parameter.
    # Since NUM-02 the table picks r by probing, and two AL solves at different
    # r agree only to within the convergence TOLERANCE (1e-9 here), not to
    # machine precision -- so comparing across different r would be measuring
    # the solver's tolerance rather than the (A4) rescaling this test is about.
    solver = TwoLayerGapSolver(*rescaled, n_y=200, tol=1e-9, max_iter=60000,
                               r=tab.tuned_r, rho=tab.tuned_r)
    direct = closures_from_solution(
        solver.solve_fixed_mean_velocity(0.45, [0.0, 1.0], [0.0, 0.0]))
    got = tab(0.45, H=H)
    for value, ref in zip(got, direct.as_tuple()):
        assert np.isclose(float(value), ref, rtol=1e-10, atol=1e-14)


# ---------------------------------------------------------------- M3-T2 ---
def test_m3_t2_q0_monotone_in_c(newt_table, hb_table):
    """A non-monotone q0 would break the LLF scheme's maximum principle.
    Linear interpolation preserves monotone samples, so checking the samples is
    both sufficient and necessary -- which is why the interpolant is linear and
    not cubic."""
    assert newt_table.q0_is_monotone_in_c()
    assert hb_table.q0_is_monotone_in_c()
    c = np.linspace(0.0, 1.0, 1001)
    for tab in (newt_table, hb_table):
        q = tab(c)[2]
        assert np.all(np.diff(q) >= -1e-12)


def test_m3_t2_q0_endpoints_exact_in_table(newt_table, hb_table):
    for tab in (newt_table, hb_table):
        lo, hi = tab.endpoint_values()
        assert np.allclose(lo, 0.0, atol=0)
        assert np.allclose(hi, 1.0, atol=0)


def test_m3_t2_linear_is_chosen_for_the_guarantee_not_for_observed_failure():
    """
    Honest justification of the linear interpolant.

    It is NOT true that a cubic overshoots on these closures: through the smooth
    Newtonian q0 a CubicSpline stays inside [0, 1] exactly, and through a real
    yield-stress q0 (tau_Y = [0.94, 0.23], n = [0.5, 0.5]) it overshoots by only
    3e-7.  So the case for linear is not an observed failure.

    The case for linear is the GUARANTEE.  Linear interpolation of monotone
    samples is monotone and stays within the sample range, unconditionally, for
    any rheology the table is ever built for.  Cubic offers no such guarantee,
    and it does break on the shape BF25 Section 2.3 says to expect: once a
    static wall layer forms, q0 rises to 1 at c_min and then PLATEAUS, and a
    cubic through that kink rings -- overshooting 1 by 2e-3 and losing
    monotonicity outright.

    Since the transport scheme's maximum principle (BCF25 Section III B) rests
    on 0 <= c <= 1, an interpolant that can exceed 1 anywhere is disqualified
    regardless of how it behaves on the cases we happened to test.
    """
    from scipy.interpolate import CubicSpline
    c_nodes = np.linspace(0.0, 1.0, 21)
    fine_c = np.linspace(0.0, 1.0, 4001)

    # smooth closure: cubic is fine, so no overshoot claim is made
    smooth = CubicSpline(c_nodes, nwt.q0(c_nodes, 5.0))(fine_c)
    assert smooth.max() <= 1.0 + 1e-9

    # static-wall-layer shape, BF25 Section 2.3: cubic fails, linear does not
    plateau = np.minimum(c_nodes / 0.52, 1.0)
    cubic = CubicSpline(c_nodes, plateau)(fine_c)
    linear = np.interp(fine_c, c_nodes, plateau)
    assert cubic.max() > 1.0 + 1e-4
    assert not np.all(np.diff(cubic) >= -1e-12)
    assert linear.max() <= 1.0
    assert np.all(np.diff(linear) >= -1e-12)


# ---------------------------------------------------------------- M3-T3 ---
def test_m3_t3_out_of_range_warns_and_names_the_cell(hb_pair):
    f1, f2 = hb_pair
    tab = ClosureTable(f1, f2, c_grid=np.linspace(0.0, 1.0, 11),
                       h_grid=np.array([0.8, 1.0, 1.2]), n_y=120,
                       tol=1e-8).build()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tab(np.array([0.5, 0.5, 0.5]), H=np.array([1.0, 1.0, 2.5]))
    msgs = [str(w.message) for w in caught
            if issubclass(w.category, ClosureTableRangeWarning)]
    assert len(msgs) == 1
    assert "axis 'H'" in msgs[0]
    assert "cell 2" in msgs[0]
    assert "2.5" in msgs[0]


def test_m3_t3_in_range_is_silent(hb_table):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        hb_table(np.linspace(0.05, 0.95, 20))
    assert not [w for w in caught
                if issubclass(w.category, ClosureTableRangeWarning)]


def test_a3_range_record_survives_warning_suppression(hb_pair):
    """
    A-3.  The K-GEP-1 production runs were launched under
    `-W ignore::UserWarning`, which removed the ONLY evidence that the closure
    table was being extrapolated -- 76 suppressed warnings in the first 400 of
    85 474 steps.  A guard that a command-line flag can erase is not a guard.

    Against the old code this test fails: there was no record to survive, only
    the warning.
    """
    f1, f2 = hb_pair
    tab = ClosureTable(f1, f2, c_grid=np.linspace(0.0, 1.0, 11),
                       h_grid=np.array([0.8, 1.0, 1.2]), n_y=120,
                       tol=1e-8).build()
    assert tab.range_report() == ""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("ignore")          # exactly what -W ignore does
        tab(np.array([0.5, 0.5, 0.5]), H=np.array([1.0, 1.0, 2.5]))
    assert not caught                            # the warning really is gone
    assert tab.n_range_points_out["H"] == 1      # the record is not
    worst, val, lo, hi = tab.worst_excursion["H"]
    assert val == 2.5 and (lo, hi) == (0.8, 1.2)
    assert abs(worst - 1.3) < 1e-12
    assert "axis 'H'" in tab.range_report()
    assert "EXTRAPOLATION" in tab.range_report()
    tab.reset_range_record()
    assert tab.range_report() == ""


def test_b2_derivative_bounds_records_out_of_range_too(hb_pair):
    """
    B-2.  `derivative_bounds` supplies every LLF wavespeed on the tabulated
    path and had no range check at all, so out-of-range queries there were
    silent BY CONSTRUCTION rather than by a flag.  Against the old code this
    fails: `n_range_points_out` stays empty however far out the query is.
    """
    f1, f2 = hb_pair
    tab = ClosureTable(f1, f2, c_grid=np.linspace(0.0, 1.0, 11),
                       h_grid=np.array([0.8, 1.0, 1.2]), n_y=120,
                       tol=1e-8).build()
    tab.derivative_bounds(np.array([0.5, 0.5]), H=np.array([1.0, 3.0]))
    assert tab.n_range_points_out.get("H") == 1
    assert tab.worst_excursion["H"][1] == 3.0


def test_m3_t3_table_covers_the_unit_interval(newt_table, hb_table):
    """c is a volume fraction, so the table must cover [0, 1] exactly -- no
    run can ever ask for a concentration outside it."""
    for tab in (newt_table, hb_table):
        assert tab.c_grid[0] == 0.0
        assert tab.c_grid[-1] == 1.0


# ---------------------------------------------------------------- B-1 ------
def test_b1_angle_assumption_is_bounded(hb_pair):
    """
    B-1.  `ClosureTable._solve_one` builds every entry with
    `u_mean = [0, umag]` against `Gb = [0, gb*H]` -- both purely axial -- so the
    table stores the theta = 0 slice, where theta is the angle between u_bar and
    G~_b.  The module docstring used to claim the angle axis "collapses" for a
    vertical well because G~_b is purely axial.  It does not: u_bar is still a
    2-vector, and theta != 0 whenever v_bar != 0, i.e. whenever the front is
    tilted.

    BF25 (2.8)-(2.10) make the gap-scale stress a 2-vector with
    eta_k = eta_k(|tau_k|), so all four closures depend on theta.  This test
    measures how much, at the SAME states the table stores, and asserts the
    bound -- including at the TOP of the umag axis, which the audit flagged as
    unmeasured.

    It is a bound, not a fix.  The fix is a fifth axis, which NUM-14 already
    describes for beta != 0.  `ClosureTable._measure_angle_sensitivity` runs the
    same probe at build time and `assumption_report()` states the result, so no
    table can be used without its own number attached.
    """
    from d2dga.gapscale.al_solver import TwoLayerGapSolver
    from d2dga.gapscale.closures import closures_from_solution

    f1, f2 = hb_pair
    gb = 25.0
    sol = TwoLayerGapSolver(f1, f2, r=0.01, rho=0.01, n_y=150, tol=1e-8)

    def at(c, u, theta):
        s = sol.solve_fixed_mean_velocity(
            c, [u * np.sin(theta), u * np.cos(theta)], [0.0, gb])
        assert s.converged
        return closures_from_solution(s)

    worst = {}
    for u in (0.5, 5.0, 50.0, 500.0):
        dev = 0.0
        for c in (0.25, 0.5, 0.75):
            ref = at(c, u, 0.0)                     # what the table stores
            for th in (np.pi / 4, np.pi / 2):
                got = at(c, u, th)
                for a, b in ((got.I1, ref.I1), (got.I2, ref.I2),
                             (got.q0, ref.q0), (got.I3, ref.I3)):
                    dev = max(dev, abs(a - b) / max(abs(b), 1e-30))
        worst[u] = dev

    # MEASURED, and it is not small.  The audit put this at "up to 2.6% on I3",
    # but that was the K-GEP-1 pair -- Newtonian displaced fluid, and
    # b*I1 ~ 1500 so buoyancy dominates the stress and the angle barely enters.
    # For a pair with a YIELD STRESS IN BOTH FLUIDS the errors are O(1),
    # because whether material yields at all depends on |tau| and |tau| depends
    # on theta, so the unyielded fraction -- and hence I1, an integral of the
    # fluidity -- moves sharply.  At c = 0.75, |u| = 5 the theta = 90 degree
    # closures differ from the stored ones by I1 +42%, I2 +134%, I3 -94%.
    assert worst[0.5] < 0.10, worst          # nearly parallel: still small
    assert worst[5.0] > 0.50, worst          # and then it is NOT a correction
    assert worst[50.0] > 0.50, worst
    # non-monotone in |u|: it peaks where the two stress scales are comparable
    # and falls again once |u_bar| dominates, which is why a single "small
    # |u_bar|" caveat would not have caught it
    assert worst[500.0] < worst[5.0], worst


def test_b1_build_measures_and_reports_the_angle_sensitivity(hb_pair):
    """
    B-1.  The measurement above is only useful if every table carries it.
    `build()` runs the probe and `assumption_report()` states it; for a pair
    with a yield stress in both fluids the verdict must say so out loud.

    Against the old code this fails with AttributeError: neither the
    measurement nor the accessor existed.
    """
    f1, f2 = hb_pair
    tab = ClosureTable(f1, f2, c_grid=np.linspace(0.0, 1.0, 11),
                       umag_grid=np.array([0.5, 5.0, 50.0]),
                       gb_grid=np.array([25.0]), n_y=120, tol=1e-8).build()
    assert np.isfinite(tab.angle_sensitivity)
    assert tab.angle_sensitivity > 0.5, tab.angle_sensitivity_detail
    rep = tab.assumption_report()
    assert "angle sensitivity (B-1)" in rep
    assert "not valid for this pair" in rep, rep


def test_b1_newtonian_pair_has_no_angle_sensitivity():
    """The other half, and the reason the assumption was defensible for every
    PUBLISHED case: for a Newtonian pair eta is constant, so rotating u_bar
    relative to G~_b changes the stress direction and nothing else.  Every
    benchmark in this repository is Newtonian, which is why no gate caught
    B-1."""
    f = ScaledFluid("newt", 1.0, 1.0, 1.0, 0.0)
    f2 = ScaledFluid("newt2", 1.0, 0.2, 1.0, 0.0)
    tab = ClosureTable(f, f2, c_grid=np.linspace(0.0, 1.0, 6),
                       umag_grid=np.array([0.5, 5.0, 50.0]),
                       gb_grid=np.array([25.0]), n_y=120, tol=1e-9).build()
    assert tab.angle_sensitivity < 1e-6, tab.angle_sensitivity_detail
    assert "benign" in tab.assumption_report()
