"""
BENCH-10 -- BF25 Section 3.3's Muskat stability analysis.

Added in remediation of audit finding A-4.3: *every* benchmark in this
repository was Newtonian (ZF22 Table 3, ZF23 Table 3, PF04), so the M3 closure
table -- and therefore the whole Herschel-Bulkley path, which is the path the
K-GEP-1 result uses -- was tested against no external analytic statement at all.

BF25 (3.12) is the one closed-form result in the six papers that applies to an
arbitrary fluid pair.  It is built from I1 and I2 only, so it exercises exactly
the two closures the HB path had no independent check on.
"""

import numpy as np
import pytest

from d2dga.gapscale import newtonian as nwt
from d2dga.scaling import zf22_buoyancy_to_bf25
from d2dga.muskat import (check_I2_at_one, classify,
                                          dw_at_leading_edge, finger_velocity,
                                          front_speed)
from tests.benchmarks.zf22_cases import CASES


def _newt(m):
    return (lambda c: nwt.script_I1(c, m), lambda c: nwt.script_I2(c, m),
            lambda c: nwt.q0(c, m), lambda c: nwt.script_I3(c, m))


# ---------------------------------------------------------------- I2(1) ---
@pytest.mark.parametrize("m", [0.006, 0.2, 1.0, 5.0, 160.0])
def test_bench10_I2_at_one_is_exactly_zero_analytic(m):
    """BF25 uses I2(1) = 0 to reach (3.12) -- "Note that we have used the fact
    that I2(1) = 0 in the above".  It was not among the endpoint identities this
    repository checked.  It is structural, not approximate: I2 = (1-c)A + cE
    with E = int_c^1 y(1-y)/eta1 dy, and both terms vanish at c = 1."""
    assert check_I2_at_one(lambda c: nwt.script_I2(c, m)) == 0.0


def test_bench10_I2_at_one_is_zero_on_a_herschel_bulkley_table(hb_table_fixture):
    """The same identity on a tabulated HB pair, which is where it could fail:
    the gap solve is iterative and the c = 1 node is a degenerate single-layer
    problem."""
    tab = hb_table_fixture
    I1, I2, q0, I3 = tab(1.0, H=1.0, umag=1.0)
    assert abs(float(I2)) < 1e-9, float(I2)
    assert float(I1) > 0.0
    assert abs(float(q0) - 1.0) < 1e-12
    assert abs(float(I3)) < 1e-9


@pytest.fixture(scope="module")
def hb_table_fixture():
    from d2dga.gapscale.tables import ClosureTable
    from d2dga.scaling import ScaledFluid
    f1 = ScaledFluid("mud", 1.0, 0.45, 0.6, 0.30)
    f2 = ScaledFluid("slurry", 1.0, 0.25, 0.8, 0.12)
    return ClosureTable(f1, f2, c_grid=np.linspace(0.0, 1.0, 41),
                        n_y=200, tol=1e-9).build()


# ------------------------------------------------- the paper's own threshold ---
@pytest.mark.parametrize("m", [0.2, 1.0, 1.4, 1.5, 1.6, 3.0, 5.0])
def test_bench10_leading_edge_criterion_reproduces_bf25_m_equals_1_5(m):
    """
    The check that pins this reading of BF25 (3.13) to the paper instead of to
    preference.

    At b = 0 the finger velocity (3.12) collapses to the mobility ratio
    I1(1)/I1(c0) -> m as c0 -> 0, so

        Delta_w(0+) = m - w_f(0+)

    and the leading wave speed w_f(0+) is where BF25's own threshold lives.
    For m <= 3/2 the flux q0 is concave, the front DISPERSES, and the leading
    wave is the characteristic q0'(0) = 3/2 exactly.  For m > 3/2 it is not:
    the upper concave envelope leaves the origin along a chord, the leading
    wave is a SHOCK, and w_f(0+) = max_c q0(c)/c > 3/2.  That transition at
    m = 3/2 is BF25 Section 3.2's `M_3^min`, attributed there to Lajeunesse
    et al. (1999) from an entirely different argument.

    So the criterion reproduces the paper in three separate ways, all asserted
    here: w_f(0+) = 3/2 exactly below the threshold and strictly above it after;
    Delta_w(0+) = m - 3/2 exactly below it; and the SIGN of Delta_w(0+) flips at
    m = 3/2 for every m.  Recovering 3/2 from (3.12) is independent evidence
    that I1 and the q0 derivative are both right.

    (An earlier version of this test asserted Delta_w = m - 3/2 for all m and
    failed at m = 3 and m = 5.  That was the test being wrong, not the code:
    above the threshold the leading wave is a shock and the equality does not
    hold.  The claim is corrected, not the tolerance.)
    """
    I1, I2, q0, I3 = _newt(m)
    wf0 = float(front_speed(q0, I3, [1e-6], 0.0)[0])
    if m <= 1.5:
        assert abs(wf0 - 1.5) < 1e-5, (m, wf0)          # dispersive
    else:
        assert wf0 > 1.5 + 1e-6, (m, wf0)               # frontal shock
        assert abs(wf0 - float(np.max(q0(np.linspace(1e-9, 1, 200001))
                                      / np.linspace(1e-9, 1, 200001)))) < 1e-4

    dw = dw_at_leading_edge(I1, I2, q0, I3, 0.0, c0=1e-6)
    if m <= 1.5:
        assert abs(dw - (m - 1.5)) < 2e-4, (dw, m - 1.5)
    assert (dw > 0) == (m > 1.5), (m, dw)


# ------------------------------------------------------ ZF22, all ten cases ---
def test_bench10_case_1_is_the_only_muskat_unstable_zf22_case():
    """
    The analytic statement BENCH-09 was missing.

    ZF22 case 1 is the only one of the ten where BF25's finger outruns the
    front's leading wave: Delta_w(0+) = +0.97 against -1.96 ... -162.9 for the
    other nine.  So the base state of case 1 is Muskat-UNSTABLE and the other
    nine are not -- an analytic result, from the paper, with no solver involved.

    That matters for how case 1's disagreement is read (audit A-5): a
    Muskat-unstable front has no mesh-independent finger width, so eta_E at a
    fixed time is not a converged functional of the mesh, and the +12.8% drift
    in eta_E between 20x200 and 40x400 is what an unstable base state produces.
    It is not evidence that the closures or the elliptic solve are wrong -- both
    are pinned to machine precision by other gates.
    """
    dw = {}
    for cs in CASES:
        b = zf22_buoyancy_to_bf25(cs.b, cs.m)
        dw[cs.case] = dw_at_leading_edge(*_newt(cs.m), b)
    assert dw[1] > 0.0, dw
    for k, v in dw.items():
        if k != 1:
            assert v < 0.0, (k, v, dw)
    # and case 1 is the only one with b < 0, which is the mechanism
    assert [cs.case for cs in CASES if cs.b < 0] == [1]


def test_bench10_bf25_3_13_literal_range_is_degenerate():
    """
    Q-3, recorded rather than papered over.

    BF25 (3.13) says "for an unstable regime, Delta_w(c_0) > 0 for all
    c_0 in [0,1]".  Taken literally that can never classify anything, because at
    c_0 -> 1

        w_f -> q0'(1) + b I3'(1) = 0 + 0 = 0     (both derivatives vanish)
        w_finger -> I1(1)/I1(1) + b I1(1)[0 + 1 - 1] = 1

    for EVERY pair and every b, so Delta_w(1-) -> +1 unconditionally: "stable"
    is unreachable and "unstable" is unfalsifiable at that end.  Every one of
    the ten ZF22 cases comes out "partial penetration" under the literal
    reading, including case 4 at b = 1000, which displaces as a near-piston.

    This test asserts the degeneracy so that the restricted reading used above
    is visibly a choice, with a measured consequence, and not a silent one.
    """
    for cs in CASES:
        b = zf22_buoyancy_to_bf25(cs.b, cs.m)
        regime_literal, dw, _ = classify(*_newt(cs.m), b, c0_max=1.0)
        assert regime_literal == "partial penetration", (cs.case, regime_literal)
        assert dw[-1] > 0.9, (cs.case, dw[-1])       # the degenerate end
    # restricted to where a front exists, the ten cases separate
    regimes = {}
    for cs in CASES:
        b = zf22_buoyancy_to_bf25(cs.b, cs.m)
        regimes[cs.case], _, _ = classify(*_newt(cs.m), b, c0_max=0.9)
    assert regimes[1] != "stable"
    assert {regimes[k] for k in (2, 3, 4, 5, 6, 9, 10)} == {"stable"}


# ------------------------------------------------------- the HB path, K-GEP-1 ---
def test_bench10_kgep1_is_muskat_stable_on_the_tabulated_hb_path():
    """
    A-4.3.  The first test in this repository that puts the K-GEP-1 Herschel-
    Bulkley closure table against an external analytic criterion.

    BF25 Section 3.1 states the closures are evaluated "imposing the velocity
    (v_bar, w_bar) = (0, 1)", so H = 1 and |u_bar| = 1 are the right axes.
    """
    import sys
    sys.path.insert(0, "scripts")
    from d2dga.config import Config, GridConfig
    from d2dga.scaling import Scaling
    import kgep1_run as K

    cfg = Config(grid=GridConfig(16, 80))
    geo = K.make_geometry(cfg, "synthetic")
    mud, cement = cfg.fluids.as_fluids()
    sc = Scaling(mud, cement, r_a_hat_star=geo.r_a_hat_star,
                 delta_star=geo.delta_star, mean_velocity=0.2)
    tab, _, _ = K.make_table(sc, geo, 31, 9)
    gb = sc.buoyancy_number

    def I1(c):
        return tab(c, H=1.0, umag=1.0, gb=gb)[0]

    def I2(c):
        return tab(c, H=1.0, umag=1.0, gb=gb)[1]

    def q0(c):
        return tab(c, H=1.0, umag=1.0, gb=gb)[2]

    def I3(c):
        return tab(c, H=1.0, umag=1.0, gb=gb)[3]

    assert check_I2_at_one(I2) < 1e-9
    dw = dw_at_leading_edge(I1, I2, q0, I3, gb)
    assert dw < 0.0, dw          # the finger is absorbed: Muskat-stable
    # and it is stable by a wide margin, not marginally
    assert dw < -10.0, dw
    regime, dwv, _ = classify(I1, I2, q0, I3, gb, n=101, c0_max=0.9)
    assert regime == "stable", (regime, dwv.min(), dwv.max())


# --------------------------------------------------- A-4.1, made quantitative ---
def test_a41_five_of_the_ten_zf22_cases_cannot_discriminate():
    """
    A-4.1.  The audit observed that half the ZF22 suite is "buoyancy-saturated",
    so agreement there does not discriminate.  That was an impression; this
    measures it.

    For each case, the EXACT 1-D entropy solution (upper concave envelope of
    f = q0 + b I3, evaluated at t = 1.2 volumes) gives eta_E as a function of b
    with no solver involved.  Its logarithmic sensitivity |d eta_E / d ln b| is
    how much the published eta_E constrains the model:

        case    b_ZF22    eta_1D    |d eta / d ln b|
          1       -50     0.5593        0.382
          8        10     0.8640        0.038
          6       100     0.9413        0.028
          7        10     0.9047        0.019
          3        10     0.9233        0.016
          2,5,9,10 100    0.9607        0.014
          4      1000     0.9840        0.007

    So for cases 2, 4, 5, 9 and 10 a **doubling of b** would move eta_E by less
    than 0.01 -- comfortably inside the +-0.025 agreement the reproduction
    claims.  Those five cases are consistent with almost any buoyancy treatment
    and cannot be cited as evidence for this one.

    And the sharp end of it: case 1 is TEN TIMES more discriminating than the
    next best case and fifty times more than the weakest -- and case 1 is the
    one that does not reproduce (BLOCKED.md BLK-6).
    """
    from tests.benchmarks.entropy_solution import exact_eta, similarity

    cg = np.linspace(0.0, 1.0, 200001)

    def sensitivity(case):
        def eta(bz):
            b = zf22_buoyancy_to_bf25(bz, case.m)
            sn, un = similarity(cg, nwt.q0(cg, case.m)
                                + b * nwt.script_I3(cg, case.m))
            return exact_eta(1.2, 1.0, sn, un)
        return abs(eta(case.b * 1.1) - eta(case.b / 1.1)) / (2 * np.log(1.1))

    sens = {cs.case: sensitivity(cs) for cs in CASES}

    # the five saturated cases: a doubling of b moves eta_E by < 0.01
    for k in (2, 4, 5, 9, 10):
        assert sens[k] * np.log(2.0) < 0.010, (k, sens[k])
    # case 1 is an order of magnitude more discriminating than anything else
    others = max(v for k, v in sens.items() if k != 1)
    assert sens[1] > 10.0 * others, (sens[1], others, sens)
