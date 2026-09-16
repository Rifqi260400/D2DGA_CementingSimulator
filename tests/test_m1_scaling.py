"""M1 test gates -- scaling."""

import numpy as np
import pytest

from d2dga.scaling import (HerschelBulkleyFluid, Scaling, bf25_buoyancy_to_zf22,
                           zf22_buoyancy_to_bf25)
from tests.benchmarks.zf22_cases import CASES, build_scaling


def _sample_scaling():
    return Scaling(
        HerschelBulkleyFluid("mud", 1200.0, 0.8, 0.6, 3.0),
        HerschelBulkleyFluid("cement", 1800.0, 1.4, 0.75, 6.0),
        r_a_hat_star=0.1106, delta_star=0.19, mean_velocity=0.35,
    )


# ---------------------------------------------------------------- M1-T1 ---
@pytest.mark.parametrize("kind", ["length", "radial", "velocity", "time", "stress",
                                  "pressure", "density", "flow_rate"])
def test_m1_t1_round_trip_identity(kind):
    """dimensional -> dimensionless -> dimensional must be the identity to
    machine precision, for every variable."""
    s = _sample_scaling()
    fwd = getattr(s, f"to_dimensionless_{kind}")
    bwd = getattr(s, f"to_dimensional_{kind}")
    values = np.array([1e-6, 1e-3, 0.5, 1.0, 7.3, 1e3, 1e6])
    assert np.allclose(bwd(fwd(values)), values, rtol=1e-15, atol=0)
    assert np.allclose(fwd(bwd(values)), values, rtol=1e-15, atol=0)


def test_m1_t1_flow_rate_velocity_consistency():
    """Q = A w0 must hold whichever of the two is supplied."""
    a = Scaling(HerschelBulkleyFluid.newtonian("a", 1000, 0.004),
                HerschelBulkleyFluid.newtonian("b", 1200, 0.02),
                r_a_hat_star=0.019845, delta_star=0.12018, mean_velocity=0.032)
    b = Scaling(HerschelBulkleyFluid.newtonian("a", 1000, 0.004),
                HerschelBulkleyFluid.newtonian("b", 1200, 0.02),
                r_a_hat_star=0.019845, delta_star=0.12018, flow_rate=a.flow_rate_hat)
    assert np.isclose(a.w0_hat, b.w0_hat, rtol=1e-15)


def test_m1_t1_bf25_scaled_consistency_relations():
    """BF25 (2.32): kappa1 = m^0.5/(1+B) and kappa2 = m^-0.5/(1+B).  These are
    algebraic identities of the scaling, so any drift means a factor is wrong."""
    s = _sample_scaling()
    assert np.isclose(s.scaled_fluid1.consistency,
                      np.sqrt(s.m) / (1.0 + s.B), rtol=1e-13)
    assert np.isclose(s.scaled_fluid2.consistency,
                      1.0 / np.sqrt(s.m) / (1.0 + s.B), rtol=1e-13)
    # scaled yield stresses must land in [0, 1), which is the point of using
    # max(tauY) in the stress scale (BF25 after 2.33)
    for f in (s.scaled_fluid1, s.scaled_fluid2):
        assert 0.0 <= f.yield_stress < 1.0


def test_m1_t1_newtonian_limit_recovers_bf25_section_2_2():
    """n=1, tau_Y=0 must give B=0, eta1 = m^0.5, eta2 = m^-0.5 (BF25 Section 2.2)."""
    s = Scaling(HerschelBulkleyFluid.newtonian("a", 1000, 0.004),
                HerschelBulkleyFluid.newtonian("b", 1200, 0.02),
                r_a_hat_star=0.019845, delta_star=0.12018, mean_velocity=0.032)
    assert s.B == 0.0
    assert np.isclose(s.m, 0.004 / 0.02)
    assert np.isclose(s.scaled_fluid1.consistency, np.sqrt(s.m), rtol=1e-14)
    assert np.isclose(s.scaled_fluid2.consistency, 1 / np.sqrt(s.m), rtol=1e-14)
    assert np.isclose(s.mu_e_hat, np.sqrt(0.004 * 0.02), rtol=1e-14)


# ---------------------------------------------------------------- M1-T2 ---
@pytest.mark.parametrize("case", CASES, ids=lambda c: f"case{c.case}")
def test_m1_t2_zf22_viscosity_ratio(case):
    """m is convention-independent for Newtonian pairs and must match exactly."""
    s = build_scaling(case)
    assert np.isclose(s.m, case.m, rtol=1e-12), f"{s.m} vs {case.m}"


@pytest.mark.parametrize("case", CASES, ids=lambda c: f"case{c.case}")
def test_m1_t2_zf22_buoyancy_number(case):
    """
    THE convention test.  Compute b in BF25 convention from Table 1, convert
    back to ZF22 convention, compare with Table 2.

    Tolerance 1% absolute-on-magnitude: Table 2 prints b to 1-2 significant
    figures (-50, 10, 100, 1000), so exact agreement is not available.
    """
    s = build_scaling(case)
    b_zf22 = bf25_buoyancy_to_zf22(s.buoyancy_number, s.m)
    assert np.isclose(b_zf22, case.b, rtol=0.02), (
        f"case {case.case}: recovered {b_zf22:.4f}, table {case.b}")


@pytest.mark.parametrize("case", CASES, ids=lambda c: f"case{c.case}")
def test_m1_t2_zf22_reynolds_is_a_rounded_design_value(case):
    """
    Re does NOT reproduce Table 2 exactly, and the discrepancy is systematic.

    With d_hat = 2.385 mm (stated explicitly in ZF23 Section 2.1) and the
    printed w0, every case comes out at 0.954 x the tabulated Re: 19.05 vs 20,
    95.4 vs 100, 954 vs 1000.  A constant ratio across three decades is a
    rounding convention, not an error in our scaling -- and it cannot be a
    d_hat mismatch, because the b column (which carries d_hat^2) reproduces to
    better than 1% with the same d_hat.

    Conclusion: Table 2's Re column holds rounded DESIGN targets.  Recorded as
    assumption CONV-07.  Re is a diagnostic only here -- the D2DGA model is
    non-inertial, so this affects nothing downstream.
    """
    s = build_scaling(case)
    ratio = s.reynolds_zf22 / case.Re
    assert np.isclose(ratio, 0.954, atol=0.005), (
        f"case {case.case}: Re {s.reynolds_zf22:.3f} vs table {case.Re}, "
        f"ratio {ratio:.4f}")


def test_m1_t2_conversion_is_an_involution():
    for m in (0.2, 0.5, 1.0, 2.0, 5.0):
        for b in (-50.0, 0.0, 10.0, 1000.0):
            assert np.isclose(bf25_buoyancy_to_zf22(zf22_buoyancy_to_bf25(b, m), m),
                              b, rtol=1e-14, atol=1e-14)


def test_m1_t2_conversion_factor_matches_bf25_own_translation():
    """
    Independent check of CONV-02, using BF25's own words rather than our algebra.

    BF25 Section 3.2: "we have b = 3 m^0.5 / U and m = M", where U is
    Lajeunesse's parameter U = 3 mu1 w0 / (h^2 g (rho1 - rho2)).  Since
    ZF22's b = (rho2-rho1) g d^2/(mu1 w0) with h = d, we have U = 3/b_ZF22,
    hence b_BF25 = 3 m^0.5 / (3/b_ZF22) = m^0.5 b_ZF22.
    """
    for case in CASES:
        s = build_scaling(case)
        mu1, mu2 = case.mu1, case.mu2
        b_zf22 = (case.rho2 - case.rho1) * 9.81 * (2.385e-3) ** 2 / (mu1 * case.w0)
        U = 3.0 / b_zf22 if b_zf22 != 0 else np.inf
        b_via_lajeunesse = 3.0 * np.sqrt(mu1 / mu2) / U
        assert np.isclose(s.buoyancy_number, b_via_lajeunesse, rtol=1e-9)


def test_m1_t2_case9_flow_rate_typo_is_documented():
    """Table 1 case 9 prints Q = 2.38e-4, inconsistent with its own w0 = 0.04.
    We drive from w0; confirm the implied Q is 2.38e-5, one decade smaller."""
    s = build_scaling(CASES[8])
    assert np.isclose(s.flow_rate_hat, 2.38e-5, rtol=0.02)


def test_buoyancy_sign_convention():
    """b > 0 for the favourable case (denser displacing lighter), b < 0 when
    the displacing fluid is lighter.  ZF22 case 1 is the density-unstable one."""
    assert build_scaling(CASES[0]).buoyancy_number < 0    # rho2 = 885 < 1000
    assert build_scaling(CASES[1]).buoyancy_number > 0    # rho2 = 1229 > 1000
