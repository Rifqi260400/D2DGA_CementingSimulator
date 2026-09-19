"""M8 gates -- BCF25 Section IV, the conservation verification.

At this stage of the project this is the verification the model rests on:
there is no CFD comparison, so what is checked is (a) the ten published ZF22
cases and (b) this ledger, which is the test BCF25 apply to the same scheme on
the same equations, reporting ~1e-15.

The gates here pin the DEFINITION, not just the magnitude. A conservation
figure is easy to make look good -- normalise by something large, take the
error at a convenient time, or let one fluid's ledger be computed from the
other's so it corroborates itself. Each of those is closed below.
"""

from __future__ import annotations

import numpy as np
import pytest

from d2dga.config import Config, GridConfig
from d2dga.conservation import ConservationLedger, total_axial_flux
from d2dga.elliptic import NewtonianClosures
from d2dga.geometry import build_geometry
from d2dga.simulation import Simulation


def _sim(n_phi=8, n_xi=40, m=0.5, delta_rho=0.0):
    geo = build_geometry(Config(grid=GridConfig(n_phi, n_xi)))
    return geo, Simulation(geo, NewtonianClosures(m=m), froude=1.0,
                           delta_rho=delta_rho)


# ---------------------------------------------------------------- M8-T1 ---
def test_m8_t1_the_ledger_is_bcf25s_definition_not_a_rescaling_of_it():
    """Vol1 − Vol2, normalised by the TOTAL ANNULUS VOLUME.

    The codebase's older `cons_err` divides by the volume *present*, which
    tends to zero at the start of a run and so inflates the early error. That
    is a stricter measure, not a wrong one -- but it is not the quantity BCF25
    report, and only one of the two can be put beside their 1e-15. This test
    fixes which is which by constructing both from the same run.
    """
    geo, sim = _sim()
    res = sim.run(t_end=0.3 * geo.grid.Z)
    led = sim.ledger

    assert led.capacity == pytest.approx(sim._capacity)
    assert len(led.err2) == len(led.times) == res.steps + 1
    assert led.err1[0] == led.err2[0] == 0.0     # t = 0 agrees by construction

    # the two normalisations differ, and in the direction predicted: dividing
    # by the volume present cannot be smaller than dividing by the whole
    # annulus, because the annulus is never fuller than full
    assert res.conservation_error >= max(abs(e) for e in led.err2) * 0.999


# ---------------------------------------------------------------- M8-T2 ---
def test_m8_t2_the_error_stays_at_bcf25s_order():
    """BCF25: "consistently O(1e-15), comparable to the standard Matlab
    precision", over five cases and both models.

    The bound here is 1e-14 rather than 1e-15 because roundoff accumulates
    with the number of steps and these meshes take more of them than a
    tolerance of exactly 1e-15 would survive; √N·eps is the scale to read it
    against. It is still two orders below anything that could be called a
    bookkeeping defect, and the margin is measured rather than assumed -- the
    assertion message prints both.
    """
    for n_xi, drho in ((40, 0.0), (40, -0.2), (80, -0.2)):
        geo, sim = _sim(n_xi=n_xi, delta_rho=drho)
        sim.run(t_end=0.3 * geo.grid.Z)
        led = sim.ledger
        n = len(led.err2) - 1
        roundoff = n ** 0.5 * np.finfo(float).eps
        assert led.worst < 1e-14, (
            f"n_xi={n_xi} drho={drho}: worst {led.worst:.2e} over {n} steps; "
            f"sqrt(N)*eps = {roundoff:.2e}")
        assert led.worst < max(200 * roundoff, 1e-15), (
            f"n_xi={n_xi} drho={drho}: {led.worst:.2e} is far above the "
            f"roundoff scale {roundoff:.2e}, so it is bookkeeping, not "
            f"arithmetic")


# ---------------------------------------------------------------- M8-T3 ---
def test_m8_t3_the_total_flux_is_closure_free_and_exact():
    """Fluid 1's ledger must not be read off fluid 2's.

    It is built from the TOTAL volumetric boundary flux, which comes from the
    stream function alone -- the advective flux with q0 = 1 telescopes. The
    check that this is the right quantity: with the annulus full of fluid 2,
    the fluid-2 face flux must equal it exactly.
    """
    geo, sim = _sim()
    for c in (sim.initial_condition(), np.ones((geo.grid.n_phi, geo.grid.n_xi))):
        psi, _ = sim.solve_stream_function(c, sim.Q(0.0))
        for j in (0, -1):
            assert total_axial_flux(psi, j) == pytest.approx(1.0, rel=1e-12)

    full = np.ones((geo.grid.n_phi, geo.grid.n_xi))
    psi, _ = sim.solve_stream_function(full, sim.Q(0.0))
    _, Xi, _, _ = sim.transport.fluxes(full, psi)
    for j in (0, -1):
        assert float(np.sum(Xi[:, j])) == pytest.approx(
            total_axial_flux(psi, j), abs=1e-14)


# ---------------------------------------------------------------- M8-T4 ---
def test_m8_t4_the_elliptic_solve_delivers_the_same_flux_at_both_ends():
    """The independent content of the two-fluid ledger.

    With K = 2, fluid 1's error is algebraically -fluid 2's PLUS the total-flux
    imbalance, so the imbalance is the only part that is not self-corroborating.
    """
    geo, sim = _sim(delta_rho=-0.2)
    sim.run(t_end=0.3 * geo.grid.Z)
    assert sim.ledger.total_flux_imbalance < 1e-14

    # and the algebra above is real, not asserted: err1 + err2 + imbalance is
    # identically zero.
    #
    #   err1 + err2 = [(V - Vol1_2) - (V - init) - S_tot + S_2
    #                  + Vol1_2 - init - S_2] / V = -S_tot / V
    #
    # In exact arithmetic that is zero.  In floating point it is not: err1 is
    # formed from V - Vol1_2 with both O(V), so the residual is cancellation
    # at the scale of eps, and the bound has to be written in units of eps
    # rather than as a constant that happens to fit.  Measured here at ~9 eps;
    # 100 eps leaves room for a longer run and still fails a real bookkeeping
    # error by three orders.
    e1 = np.asarray(sim.ledger.err1)
    e2 = np.asarray(sim.ledger.err2)
    im = np.asarray(sim.ledger.imbalance)
    residual = float(np.abs(e1 + e2 + im).max())
    eps = np.finfo(float).eps
    assert residual < 100 * eps, (
        f"the identity err1 + err2 = -imbalance held only to "
        f"{residual / eps:.0f} eps; it is exact in exact arithmetic, so "
        f"anything much past roundoff means the two ledgers are not measuring "
        f"the same thing")


# ---------------------------------------------------------------- M8-T5 ---
def test_m8_t5_the_ledger_is_not_reset_by_a_resume(tmp_path):
    """A conservation figure that restarts at a checkpoint would show a clean
    run that was not clean.

    Note what is NOT asserted: an equal step count. Stopping at 0.15 Z
    truncates the last step to land exactly on it, so the resumed run takes
    one step more than an uninterrupted one -- an inherent consequence of
    checkpointing at a time that is not step-aligned, not a ledger defect.
    What must hold is that the series covers the WHOLE run from t = 0 and that
    the error is the same size, which is what a reader takes from the figure.
    """
    geo, _ = _sim()
    ckpt = str(tmp_path / "c.npz")

    _, whole = _sim()
    whole.run(t_end=0.3 * geo.grid.Z)

    _, a = _sim()
    a.run(t_end=0.15 * geo.grid.Z, checkpoint_path=ckpt)
    partial_len = len(a.ledger.err2)
    _, b = _sim()
    b.run(t_end=0.3 * geo.grid.Z, checkpoint_path=ckpt)

    assert not b.ledger.resumed_without_history
    assert b.ledger.times[0] == 0.0               # from the start, not the resume
    assert len(b.ledger.err2) > partial_len       # and it grew, not restarted
    assert abs(len(b.ledger.err2) - len(whole.ledger.err2)) <= 1
    assert b.ledger.times[-1] == pytest.approx(whole.ledger.times[-1], rel=1e-12)
    # same order, and both far below the gate in M8-T2
    assert b.ledger.worst == pytest.approx(whole.ledger.worst, rel=0.5)
    assert b.ledger.worst < 1e-14


def test_m8_t6_a_pre_ledger_checkpoint_is_reported_as_partial(tmp_path):
    """Resuming a checkpoint written before the ledger existed gives a series
    covering only the steps since.  It must SAY so, not present a partial run
    as a whole one."""
    import numpy as np

    geo, a = _sim()
    ckpt = str(tmp_path / "c.npz")
    a.run(t_end=0.15 * geo.grid.Z, checkpoint_path=ckpt)

    # strip the ledger arrays, as an older writer would have left it
    with np.load(ckpt) as z:
        keep = {k: z[k] for k in z.files if not k.startswith("ledger_")}
    np.savez_compressed(ckpt, **keep)

    _, b = _sim()
    b.run(t_end=0.3 * geo.grid.Z, checkpoint_path=ckpt)
    assert b.ledger.resumed_without_history
    assert "PARTIAL" in b.ledger.report()


# ---------------------------------------------------------------- M8-T7 ---
def test_m8_t7_the_ledger_catches_a_flux_that_is_not_accounted_for():
    """A conservation check that cannot fail is decoration.

    Injecting a volume of fluid 2 that no boundary flux reports must show up
    in the ledger at the size injected, divided by the annulus volume.
    """
    # A real step has flux_total = 0 NET: the annulus is full and
    # incompressible, so whatever volume enters at one end leaves at the
    # other.  10 units of fluid 2 in therefore means 10 units of fluid 1 out.
    led = ConservationLedger(capacity=100.0, vol1_2_initial=0.0)
    led.start(0.0, 0.0)
    led.record(1.0, vol1_2=10.0, flux_2=10.0, flux_total=0.0)
    assert led.worst == pytest.approx(0.0, abs=1e-15)

    led.record(2.0, vol1_2=20.5, flux_2=10.0, flux_total=0.0)   # 0.5 unbooked
    assert led.err2[-1] == pytest.approx(0.5 / 100.0)
    assert led.err1[-1] == pytest.approx(-0.5 / 100.0)
    assert led.worst == pytest.approx(0.005)

    # and a total flux that does not balance shows up as the imbalance, which
    # is the part fluid 1's ledger adds over fluid 2's
    led2 = ConservationLedger(capacity=100.0, vol1_2_initial=0.0)
    led2.start(0.0, 0.0)
    led2.record(1.0, vol1_2=10.0, flux_2=10.0, flux_total=2.0)
    assert led2.total_flux_imbalance == pytest.approx(0.02)
    assert led2.err2[-1] == pytest.approx(0.0, abs=1e-15)
    assert led2.err1[-1] == pytest.approx(-0.02)
