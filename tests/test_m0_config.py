"""M0 test gates -- the configuration dataclasses.

These exist because of UI gap G-1: `FluidsConfig` gained a non-Newtonian mud on
2026-09-19, and the whole point of that change is that it must be *possible*
without being *silent*.  Three things are locked here.

  1. The default pair is unchanged to the bit, so every number in
     AUDIT_REPORT.md, docs/gate_status.md and output/kgep1_results.md still
     describes the fluids the code builds today.
  2. An incoherent fluid cannot be constructed at all.
  3. A yield-stress mud makes finding B-1 O(1), and the table's own
     `assumption_report()` must SAY so rather than reporting it as benign.

Point 3 is the load-bearing one.  It asserts that the hazard is surfaced, not
that the physics is fine -- it is not fine, and the test's job is to guarantee
nobody can read a number off that path without the warning attached.
"""

import numpy as np
import pytest

from d2dga.config import Config, FluidsConfig
from d2dga.scaling import HerschelBulkleyFluid

# A plausible field mud: 11.7 ppg, and the Herschel-Bulkley fit quoted for a
# bentonite/polymer system.  Used wherever a test needs a mud that is NOT water.
FIELD_MUD = dict(mud_density=1400.0, mud_consistency=0.020,
                 mud_power_law_index=0.70, mud_yield_stress=4.79)


# ---------------------------------------------------------------- M0-C1 ---
def test_m0_c1_default_fluids_are_bit_identical_to_the_pre_g1_code():
    """Regression lock on the G-1 change.

    Before 2026-09-19 `as_fluids` read

        HerschelBulkleyFluid("drilling fluid (water)", mud_density,
                             mud_viscosity, 1.0, 0.0)

    with n and tau_Y written as literals.  That construction is reproduced here
    verbatim.  Equality is not enough on its own: `repr` is compared too,
    because `scripts/kgep1_run.make_table` hashes the *repr* of the scaled
    fluids into the closure table's cache key, so a changed name or a float that
    merely compares equal would silently invalidate every cached table and every
    recorded run tag.
    """
    mud, cement = FluidsConfig().as_fluids()
    pre_g1_mud = HerschelBulkleyFluid("drilling fluid (water)", 998.0, 1.0e-3,
                                      1.0, 0.0)
    pre_g1_cement = HerschelBulkleyFluid("cement slurry", 1200.0, 0.6, 0.4, 1.4)
    assert mud == pre_g1_mud
    assert cement == pre_g1_cement
    assert repr(mud) == repr(pre_g1_mud)
    assert repr(cement) == repr(pre_g1_cement)
    assert mud.is_newtonian and not cement.is_newtonian


# ---------------------------------------------------------------- M0-C2 ---
def test_m0_c2_mud_viscosity_refuses_when_there_is_no_single_viscosity():
    """`mud_viscosity` used to be a stored field.  As a property it must return
    kappa_hat for a Newtonian mud and REFUSE otherwise -- returning kappa_hat in
    Pa s^n under a name that means Pa s is exactly the unit error this change
    was meant to make impossible."""
    assert FluidsConfig().mud_viscosity == 1.0e-3
    assert FluidsConfig(mud_consistency=2.5e-3).mud_viscosity == 2.5e-3

    with pytest.raises(ValueError, match="no single viscosity"):
        FluidsConfig(**FIELD_MUD).mud_viscosity
    # a yield stress ALONE (Bingham) is enough; so is n alone
    with pytest.raises(ValueError, match="no single viscosity"):
        FluidsConfig(mud_yield_stress=4.79).mud_viscosity
    with pytest.raises(ValueError, match="no single viscosity"):
        FluidsConfig(mud_power_law_index=0.7, mud_consistency=0.02).mud_viscosity


# ---------------------------------------------------------------- M0-C3 ---
@pytest.mark.parametrize("who", ["mud", "cement"])
@pytest.mark.parametrize("suffix,bad,match", [
    ("_density", 0.0, "must be > 0"),
    ("_density", -1200.0, "must be > 0"),
    ("_consistency", 0.0, r"must be > 0 Pa s\^n"),
    ("_consistency", -0.6, r"must be > 0 Pa s\^n"),
    ("_power_law_index", 0.0, r"must be in \(0, 1\]"),
    ("_power_law_index", 1.4, r"must be in \(0, 1\]"),
    ("_power_law_index", -0.4, r"must be in \(0, 1\]"),
    ("_yield_stress", -1e-12, "must be >= 0 Pa"),
])
def test_m0_c3_incoherent_fluids_are_refused_at_construction(who, suffix, bad,
                                                             match):
    """Symmetric in the two fluids on purpose: the cement fields had no
    validation before, and a UI that can edit the mud can edit the cement.

    n > 1 is refused rather than extrapolated.  It is a real class of fluid, but
    BF25 (2.8)-(2.10) and B02 (10) derive the gap-scale closures for
    shear-thinning Herschel-Bulkley, n in (0, 1]; nothing in this codebase is
    validated above 1, so accepting it would return numbers with no authority.
    """
    with pytest.raises(ValueError, match=match):
        FluidsConfig(**{who + suffix: bad})


# ---------------------------------------------------------------- M0-C4 ---
def test_m0_c4_mean_velocities_must_be_a_non_empty_set_of_positive_speeds():
    with pytest.raises(ValueError, match="must not be empty"):
        FluidsConfig(mean_velocities_m_s=())
    with pytest.raises(ValueError, match="must be > 0 m/s"):
        FluidsConfig(mean_velocities_m_s=(0.05, 0.0, 0.5))
    with pytest.raises(ValueError, match="must be > 0 m/s"):
        FluidsConfig(mean_velocities_m_s=(-0.2,))


# ---------------------------------------------------------------- M0-C5 ---
@pytest.mark.parametrize("kwargs,expected", [
    ({}, "drilling fluid (water)"),
    ({"mud_density": 1400.0}, "drilling fluid (Newtonian)"),
    ({"mud_consistency": 2.0e-3}, "drilling fluid (Newtonian)"),
    ({"mud_yield_stress": 4.79}, "drilling fluid (Bingham)"),
    ({"mud_power_law_index": 0.7, "mud_consistency": 0.02},
     "drilling fluid (power-law)"),
    (FIELD_MUD, "drilling fluid (Herschel-Bulkley)"),
])
def test_m0_c5_mud_name_states_the_rheology_it_actually_has(kwargs, expected):
    """The label is derived, never asserted.  It matters twice: it goes into the
    closure table's cache key, and it is what the UI prints next to the fluid.
    A mud at 1400 kg/m3 must not still be called water."""
    assert FluidsConfig(**kwargs).mud_name == expected
    assert FluidsConfig(**kwargs).as_fluids()[0].name == expected


# ---------------------------------------------------------------- M0-C6 ---
def test_m0_c6_from_record_round_trips_a_payload_written_today():
    from d2dga import runio
    for fluids in (FluidsConfig(), FluidsConfig(**FIELD_MUD)):
        payload = runio.config_payload(Config(fluids=fluids), {"w0": 0.2})
        # through JSON, because that is how the file is actually read back
        import json
        recorded = json.loads(json.dumps(payload))["config"]["fluids"]
        assert FluidsConfig.from_record(recorded) == fluids


# ---------------------------------------------------------------- M0-C7 ---
def test_m0_c7_from_record_migrates_a_pre_g1_payload():
    """Runs recorded between R-2 and G-1 carry `mud_viscosity` and no rheology
    fields.  They were Newtonian by construction, so the migration is exact --
    and must produce a config equal to the default, not merely similar."""
    legacy = {
        "mud_density": 998.0, "mud_viscosity": 1.0e-3,
        "cement_density": 1200.0, "cement_consistency": 0.6,
        "cement_power_law_index": 0.4, "cement_yield_stress": 1.4,
        "mean_velocities_m_s": [0.05, 0.2, 0.5],
    }
    assert FluidsConfig.from_record(legacy) == FluidsConfig()
    # and the fluids it builds are the ones that run produced
    assert FluidsConfig.from_record(legacy).as_fluids() == FluidsConfig().as_fluids()


# ---------------------------------------------------------------- M0-C8 ---
def test_m0_c8_from_record_refuses_a_key_it_does_not_understand():
    """A future field dropped silently would reconstruct a DIFFERENT run and
    label it reproduced.  The reader is told to check out the recorded git SHA
    instead."""
    with pytest.raises(ValueError, match="unknown key"):
        FluidsConfig.from_record({"mud_gel_strength_pa": 7.0})
    with pytest.raises(ValueError, match="git_sha"):
        FluidsConfig.from_record({"mud_slip_coefficient": 0.1})


# ---------------------------------------------------------------- M0-C9 ---
def _tiny_table(fluids, mean_velocity=0.2):
    """One 5x1x3x1 closure table for a fluid pair, built (not loaded) so that
    the B-1 angle probe runs.  Small on purpose: the probe samples three |u| at
    three c, and its verdict does not depend on the table's own resolution."""
    from d2dga.config import GridConfig
    from d2dga.gapscale.tables import ClosureTable
    from d2dga.geometry import ConstantEccentricity, Geometry, UniformWall
    from d2dga.scaling import Scaling

    well = Config().well
    geo = Geometry(well, UniformWall(well.gauge_hole_diameter_m),
                   ConstantEccentricity(0.3), GridConfig(8, 20))
    mud, cement = fluids.as_fluids()
    sc = Scaling(mud, cement, r_a_hat_star=geo.r_a_hat_star,
                 delta_star=geo.delta_star, mean_velocity=mean_velocity)
    tab = ClosureTable(sc.scaled_fluid1, sc.scaled_fluid2,
                       c_grid=np.linspace(0.0, 1.0, 5),
                       h_grid=np.array([1.0]),
                       umag_grid=np.array([0.0, 1.0, 10.0]),
                       gb_grid=np.array([sc.buoyancy_number]),
                       n_y=120, tol=1e-8).build()
    assert tab.n_failed == 0, f"{tab.n_failed} gap solves failed to converge"
    return sc, tab


def test_m0_c9_a_newtonian_pair_has_no_angle_dependence_at_all():
    """The one case where the theta = 0 table is EXACT.

    With eta constant in both fluids, rotating tau changes its direction and
    nothing else, so all four closures are invariant.  Asserted at two very
    different viscosity ratios, because the next gate shows the ratio is what
    governs the non-Newtonian case -- if it leaked into this one, the claim
    would be about something else.
    """
    for kappa1 in (1.0e-3, 2.0e-2):
        _, tab = _tiny_table(FluidsConfig(
            mud_consistency=kappa1, cement_consistency=0.6,
            cement_power_law_index=1.0, cement_yield_stress=0.0))
        assert tab.angle_sensitivity == pytest.approx(0.0, abs=1e-12), (
            f"kappa1={kappa1}: a Newtonian pair must be angle-exact, got "
            f"{tab.angle_sensitivity:.3e}\n{tab.assumption_report()}")


# --------------------------------------------------------------------- M0-C10 ---
def test_m0_c10_b1_is_governed_by_the_viscosity_ratio_not_by_the_muds_rheology():
    """The correction recorded on 2026-09-19, locked as a gate.

    `tables.py` used to explain the benign K-GEP-1 measurement by saying "the
    displaced fluid is Newtonian", and AUDIT_REPORT.md B-1 repeats it.  That
    attribution is wrong.  Holding the mud strictly Newtonian (n = 1,
    tau_Y = 0) and moving only its viscosity against the unchanged K-GEP-1
    cement gives 0.53% at 1 mPa s, 11.4% at 5 mPa s and 93.1% at 20 mPa s.

    So the shipped pair is benign because m ~ 0.006, not because the mud is
    Newtonian -- and G-1's non-Newtonian fields did not create this hazard, they
    only made it easier to reach.  An ordinary 5 mPa s Newtonian mud, expressible
    before G-1 existed, is already outside the benign band.

    All rows share the same probe axes, so the sweep is apples-to-apples; the
    ABSOLUTE numbers are not portable, because the sensitivity grows with |u_bar|
    and `_tiny_table`'s axis stops at 10.  On the production axes (umag to 3000)
    the same shipped pair measures 9.8% -- SIGNIFICANT, not benign -- which is
    recorded in `tables.py` and FLU-07 and is why nothing here asserts that the
    shipped pair is benign in general.  What this gate locks is the DEPENDENCE:
    monotone in m, small at m ~ 0.006 and O(1) by m ~ 0.13, with the mud held
    strictly Newtonian throughout.

    Thresholds are one-sided and far from the measured values (< 5% vs 0.53%;
    > 50% vs 93%) so the gate tests the claim, not the solver's tolerance.
    """
    devs = {}
    for kappa1 in (1.0e-3, 5.0e-3, 2.0e-2):
        sc, tab = _tiny_table(FluidsConfig(mud_consistency=kappa1))
        assert sc.scaled_fluid1.yield_stress == 0.0
        assert sc.scaled_fluid1.power_law_index == 1.0
        devs[kappa1] = tab.angle_sensitivity
        if kappa1 == 1.0e-3:
            shipped_report = tab.assumption_report()

    assert devs[1.0e-3] < 0.05, (
        f"on this probe axis the shipped pair measured {devs[1.0e-3]:.2%}; it "
        f"was 0.53%.  Note this is NOT a claim that the pair is benign on the "
        f"production axis, where it measures 9.8% -- see FLU-07.")
    assert "benign" in shipped_report
    assert devs[1.0e-3] < devs[5.0e-3] < devs[2.0e-2], (
        f"B-1 must grow monotonically with the viscosity ratio: {devs}")
    assert devs[2.0e-2] > 0.5, (
        f"a NEWTONIAN mud at 20 mPa s measured only {devs[2.0e-2]:.2%}; the "
        f"correction in tables.py claims 93%.  Re-measure and rewrite it.")


# --------------------------------------------------------------------- M0-C11 ---
def test_m0_c11_departures_from_the_validated_pair_are_named_not_merely_flagged():
    """What the UI's refusal banner is built on.  A screen must be able to say
    WHICH field moved, because "not the validated pair" alone does not tell a
    reader whether to distrust the efficiency or the breakthrough time."""
    assert FluidsConfig().is_validated_pair
    assert FluidsConfig().departures_from_validated() == []

    f = FluidsConfig(**FIELD_MUD)
    assert not f.is_validated_pair
    assert f.departures_from_validated() == [
        "mud_density: 998 -> 1400",
        "mud_consistency: 0.001 -> 0.02",
        "mud_power_law_index: 1 -> 0.7",
        "mud_yield_stress: 0 -> 4.79",
    ]
    # a cement change must be named too -- the mud is not the only way out
    assert FluidsConfig(cement_density=1900.0).departures_from_validated() == [
        "cement_density: 1200 -> 1900"]
