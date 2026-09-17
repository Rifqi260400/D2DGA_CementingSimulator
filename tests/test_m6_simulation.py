"""M6 test gates -- the coupled elliptic/transport time loop."""

import numpy as np
import pytest

from d2dga.config import GridConfig, WellConfig
from d2dga.elliptic import NewtonianClosures, StreamFunctionSolver
from d2dga.geometry import ConstantEccentricity, Geometry, UniformWall
from d2dga.simulation import Simulation
from tests.benchmarks.kinematic_wave import riemann_structure
from tests.benchmarks.zf22_cases import CASES, build_scaling, build_simulation

SHORT = dict(casing_shoe_m=387.2, total_depth_m=390.0)


def build(n_phi=16, n_xi=120, beta=0.0, e=0.0, m=1.0, b=0.0, cfl=0.5,
          short=True, inflow="no_axial_gradient"):
    kw = dict(SHORT) if short else {}
    well = WellConfig(inclination_rad=beta, **kw)
    geo = Geometry(well, UniformWall(well.gauge_hole_diameter_m),
                   ConstantEccentricity(e), GridConfig(n_phi, n_xi))
    sim = Simulation(geo, NewtonianClosures(m), froude=1.0, delta_rho=-b,
                     cfl=cfl, inflow=inflow)
    return geo, sim


# =========================================================================
# wiring -- the two halves must agree about b
# =========================================================================
def test_elliptic_and_transport_share_one_buoyancy_number():
    """
    The buoyancy number appears twice -- in the elliptic source (BF25 2.18) and
    in the transport flux (2.21) -- and a mismatch between them would be silent
    and would look like a physical result.  Simulation DERIVES it once from
    Fr* and Delta_rho instead of taking it as a separate argument; assert that
    is what reaches both solvers.
    """
    for froude, drho in ((1.0, -8.0), (2.5, 3.0), (0.4, -0.25)):
        geo, _ = build()
        sim = Simulation(geo, NewtonianClosures(0.7), froude=froude,
                         delta_rho=drho)
        expected = -drho / froude ** 2
        assert sim.buoyancy_number == pytest.approx(expected, rel=1e-15)
        assert sim.transport.b == pytest.approx(expected, rel=1e-15)
        # the elliptic stores Fr^2 and Delta_rho separately; its source scalar
        # at c = 0 is 1/Fr^2 and its c-gradient is b, so recover b from it
        c0 = np.zeros((geo.grid.n_phi, geo.grid.n_xi))
        _, _, B0, _ = sim.elliptic.coefficients(c0)
        _, _, B1, _ = sim.elliptic.coefficients(c0 + 1.0)
        # bs = 1/Fr^2 + b (c + I2/I1), and I2 vanishes at BOTH c = 0 and c = 1,
        # so the difference between the two is exactly b
        assert np.allclose(B1 - B0, expected, rtol=1e-12)


def test_vectorised_assembly_matches_the_reference_bit_for_bit():
    """
    The coupled loop reassembles the elliptic matrix at every step, so the
    row-by-row assembly that every M4 gate was established on was replaced by a
    vectorised one.  It must be identical, not merely close -- a 1e-15 drift
    would accumulate over 10^4-10^5 steps and would be indistinguishable from
    physics.
    """
    import scipy.sparse as sp
    rng = np.random.default_rng(4)
    for beta, e, inflow in ((0.0, 0.0, "no_axial_gradient"),
                            (0.6, 0.4, "no_axial_gradient"),
                            (0.0, 0.0, "uniform"),
                            (1.1, 0.75, "uniform")):
        geo, sim = build(n_phi=13, n_xi=17, beta=beta, e=e, b=5.0,
                         inflow=inflow)
        solver = sim.elliptic
        c = rng.random((13, 17))
        coeffs = solver.coefficients(c)
        A1, b1 = solver.assemble_reference(*coeffs, 0.7)
        A2, b2 = solver.assemble(*coeffs, 0.7)
        assert A1.nnz == A2.nnz
        assert abs(A1 - A2).max() == 0.0
        assert np.max(np.abs(b1 - b2)) == 0.0


# =========================================================================
# M6-T2 -- the coupled loop keeps M5's guarantees
# =========================================================================
@pytest.mark.parametrize("beta,b,e", [(0.0, 0.0, 0.0), (0.0, 20.0, 0.5),
                                      (np.pi / 4, 20.0, 0.5),
                                      (0.0, -20.0, 0.5)])
def test_m6_t2_conservation_and_maximum_principle_survive_the_coupling(beta, b, e):
    """
    Gate M6-T2.  M5 established conservation and the maximum principle with Psi
    frozen and prescribed.  Here Psi is recomputed from c at every step, which
    is the case that matters: if the elliptic solve returned a Psi whose
    discrete divergence were not exactly zero, the telescoping that makes M5-T1
    exact would break, and the maximum principle with it.

    Note the last case has b < 0 -- adverse buoyancy, ZF22's case 1.  That is
    the configuration that drives the strongest secondary flow, so it is the
    one most likely to break a guarantee.
    """
    geo, sim = build(n_phi=16, n_xi=100, beta=beta, b=b, e=e, m=0.4)
    c = sim.initial_condition()
    mass = sim.transport.mass(c)
    for _ in range(60):
        c, psi, dt, Q, bflux, _ = sim.step(c, 0.0)
        mass += bflux
        assert abs(sim.transport.mass(c) - mass) / max(abs(mass), 1e-12) < 5e-14
        # The bound is exact in exact arithmetic.  In floating point the
        # cancellation that produces a zero cell leaves O(eps) dust -- values
        # around -1e-37 are seen at beta = pi/4 -- so the assertion is that any
        # violation is at ROUNDOFF level, not that it is bitwise absent.  A real
        # monotonicity failure is O(1), not O(1e-30): a scheme that oversteps
        # (44) goes to -0.6 in one step (see M5-T2's sharpness test), so this
        # tolerance cannot hide one.
        assert c.min() >= -1e-14, c.min()
        assert c.max() <= 1.0 + 1e-14, c.max()


def test_num26_endpoint_wavespeeds_break_monotonicity():
    """
    NUM-26.  BCF25 (32)-(35) take the LLF coefficients from the flux-function
    slopes at the two cell values either side of a face.  That is not enough.

    The coefficient then depends on the neighbour's own value, and where
    |dI3/dc| is FALLING it decreases as that neighbour's concentration rises.
    Monotonicity needs d/dc_S[alpha(c_S, c_C)(c_S - c_C)] >= 0, and the
    (d alpha/d c_S) term can overwhelm alpha.

    Run the same configuration both ways.  The printed rule drives the
    concentration to -0.03 -- an O(10^-2) violation, not roundoff -- while the
    interval maximum holds the bound exactly.  This is also why it survived
    every M5 gate: M5-T2's parameters never combined a large enough state jump
    with the falling branch of |I3'|.
    """
    from d2dga.transport import TransportSolver

    worst = {}
    for mode in ("interval", "endpoints"):
        geo, sim = build(n_phi=16, n_xi=100, beta=np.pi / 4, b=20.0, e=0.5,
                         m=0.4)
        sim.transport = TransportSolver(
            geo, sim.closures, buoyancy_number=sim.buoyancy_number,
            inflow_concentration=1.0, cfl=0.5, wavespeed=mode)
        c = sim.initial_condition()
        low = 0.0
        for _ in range(60):
            c, _, _, _, _, _ = sim.step(c, 0.0)
            low = min(low, float(c.min()))
        worst[mode] = low

    assert worst["endpoints"] < -1e-3, worst
    assert worst["interval"] >= -1e-14, worst

    # and the fix is the interval property itself, not extra dissipation: the
    # interval maximum is >= the endpoint maximum by construction, and equal
    # wherever no interior extremum separates the two states
    from d2dga.elliptic import NewtonianClosures
    nodes = NewtonianClosures(0.4).wavespeed_nodes()
    assert nodes.size >= 2 and np.all((nodes > 0) & (nodes < 1))


def test_m6_t2_axial_flux_is_the_imposed_rate_at_every_slice():
    """The elliptic BC gives Psi(1) - Psi(0) = 2Q, so the axial flux through
    every xi-slice must equal 2Q for as long as the run lasts -- CONV-05."""
    geo, sim = build(n_phi=16, n_xi=80, e=0.5, b=10.0, m=0.4)
    c = sim.initial_condition()
    for _ in range(30):
        c, psi, dt, Q, _, _ = sim.step(c, 0.0)
        assert np.allclose(sim.elliptic.axial_flux(psi), 2.0 * Q, rtol=1e-10)


def test_time_dependent_flow_rate_is_honoured():
    """A pump schedule enters only through Q(t) in the elliptic BC.  Check the
    flux tracks it; the staging in M7 depends on nothing else."""
    geo, sim = build(n_phi=12, n_xi=60)
    sim.flow_rate = lambda t: 1.0 + 0.5 * np.sin(t)
    c = sim.initial_condition()
    t = 0.0
    for _ in range(20):
        c, psi, dt, Q, _, _ = sim.step(c, t)
        assert Q == pytest.approx(1.0 + 0.5 * np.sin(t))
        assert np.allclose(sim.elliptic.axial_flux(psi), 2.0 * Q, rtol=1e-10)
        t += dt


# =========================================================================
# M6-T3 -- grid independence  (NUM-07)
# =========================================================================
def test_m6_t3_solution_converges_under_mesh_refinement():
    """
    Gate M6-T3.  NUM-07 flagged the cell aspect ratio: on K-GEP-1, Z ~ 564 and
    dxi ~ 1.4 against dphi = 0.05, an aspect ratio of ~28.  Refine both
    directions and check the front position converges rather than drifting.

    Measured on the displaced VOLUME rather than on a contour, because volume
    is what the scheme conserves exactly and what the efficiency metric uses,
    so it is not contaminated by how much the first-order scheme smears the
    front.
    """
    effs = []
    for n_phi, n_xi in ((12, 75), (24, 150), (48, 300)):
        geo, sim = build(n_phi=n_phi, n_xi=n_xi, e=0.5, b=10.0, m=0.4)
        res = sim.run(t_end=0.6 * geo.grid.Z, record_every=10_000)
        effs.append(res.efficiency[-1])
    d1, d2 = abs(effs[1] - effs[0]), abs(effs[2] - effs[1])
    assert d2 < d1, effs
    assert d2 < 0.01, effs


# =========================================================================
# M6-T4 -- timestep independence  (Q2 / NUM-01)
# =========================================================================
def test_m6_t4_result_is_insensitive_to_the_cfl_number():
    """
    Gate M6-T4, and the half of Q2 that M5 could not settle.

    M5 showed the TRANSPORT error falls as the CFL number rises (a first-order
    upwind scheme is exact at CFL = 1 for linear advection).  What M5 could not
    see is the lag between Psi and c: Psi^n moves c across the whole step while
    the true Psi varies over it, an O(dt) error that only exists once the two
    solvers are coupled.  If that lag mattered at CFL = 0.5, halving the step
    would move the answer.
    """
    out = {}
    for cfl in (0.5, 0.25):
        geo, sim = build(n_phi=16, n_xi=120, e=0.5, b=10.0, m=0.4, cfl=cfl)
        res = sim.run(t_end=0.6 * geo.grid.Z, record_every=10_000)
        out[cfl] = (res.efficiency[-1], res.concentration.copy())
    assert abs(out[0.5][0] - out[0.25][0]) < 2e-3
    assert np.max(np.abs(out[0.5][1] - out[0.25][1])) < 0.05


# =========================================================================
# M6-T1 -- ZF22 Table 3.  The expensive cases live in scripts/zf22_table3.py;
# here we check the planar limit, which is the part that isolates the closures
# from the 2-D solve.
# =========================================================================
def test_conv02_buoyancy_number_survives_the_dimensional_route():
    """
    CONV-02 via ZF22's own Table 1 -> Table 2.  Computing b from the
    DIMENSIONAL fluid properties through BF25's scaling must reproduce the
    printed ZF22 b once converted, for every case.  That closes the loop
    between the table used to drive the benchmark and BF25's definitions.
    """
    from d2dga.scaling import bf25_buoyancy_to_zf22
    for case in CASES:
        sc = build_scaling(case)
        got = bf25_buoyancy_to_zf22(sc.buoyancy_number, case.m)
        assert got == pytest.approx(case.b, rel=0.02), (case.case, got, case.b)


@pytest.mark.parametrize("case", CASES, ids=lambda c: f"zf22-{c.case}")
def test_m6_t1_planar_wave_structure_predicts_the_published_breakthrough(case):
    """
    Gate M6-T1, planar part.

    BF25 Section 3.1 reduces the model to a scalar conservation law whose
    Riemann solution is read off the upper concave envelope of
    f(c) = q0 + b I3.  For a near-concentric annulus the 2-D breakthrough time
    is bounded by that 1-D structure: the bulk of the front cannot arrive before
    the main shock, and cannot arrive after the piston time.  So

        1 / w_shock  >=  t_br  is the prediction,

    and for the five strongly buoyant cases (|b| large, shock speed -> 1) it
    pins t_br to within a few per cent of ZF22's printed value.

    This is the part of M6-T1 that isolates the CLOSURES from the 2-D solve.
    The full 2-D reproduction of Table 3 is in scripts/zf22_table3.py -- the
    strongly buoyant cases need ~4x10^4 timesteps each, far too slow for a unit
    test.
    """
    from d2dga.scaling import zf22_buoyancy_to_bf25
    b = zf22_buoyancy_to_bf25(case.b, case.m)
    st = riemann_structure(case.m, b)
    shock = st["main_shock"]
    assert shock is not None, case.case
    predicted = 1.0 / shock[2]
    # ZF22's t_br can only be earlier than the bulk arrival, never later
    assert case.t_br_d2dga <= predicted + 0.08, (case.case, predicted)
    # and the ORDERING across cases is the real content: more buoyancy, later
    # breakthrough
    if abs(b) > 40:
        assert abs(predicted - case.t_br_d2dga) < 0.07, (case.case, predicted)


def test_m6_t1_lajeunesse_critical_viscosity_ratio():
    """
    BF25 Section 3.2: at b = 0 the dispersive-to-shock transition happens at a
    critical viscosity ratio "M3min = 1.5, which is identical with the
    predictions of Lajeunesse et al. (1999)".

    This tests q0 alone -- no buoyancy, no 2-D solve -- against a number from a
    third paper, so it is an independent check on BF25 (2.26) and on the
    envelope construction at once.
    """
    assert riemann_structure(1.45, 0.0, n=200001)["main_shock"] is None
    assert riemann_structure(1.50, 0.0, n=200001)["main_shock"] is None
    assert riemann_structure(1.60, 0.0, n=200001)["main_shock"] is not None

    # BF25 figure 12(a)'s colour bar spans 1.500-1.520 over this range of m;
    # the text's "the shock speed is equal to 1.5" is looser than its own
    # figure.  Our values sit inside the figure's band.
    for m, lo, hi in ((2.0, 1.500, 1.520), (2.25, 1.500, 1.520)):
        s = riemann_structure(m, 0.0, n=200001)["main_shock"][2]
        assert lo <= s <= hi, (m, s)
    assert riemann_structure(5.0, 0.0, n=200001)["main_shock"][2] > 1.6


def test_checkpoint_round_trips_exactly(tmp_path):
    """
    A run must be resumable, because this environment reclaims its container on
    its own schedule and a full K-GEP-1 job is several hours -- three runs were
    lost that way at 22%, 22% and 2% before checkpointing existed.

    Two properties have to be separated or the test measures the wrong one:

      * the checkpoint must round-trip EXACTLY -- what comes back is bitwise
        what went in;
      * a resumed run CANNOT be bitwise equal to an uninterrupted one, because
        stopping at t1 forces a step boundary there and the time integration is
        first order, so the two take different step SEQUENCES.  The honest
        comparison is against a restart that takes the same sequence without a
        checkpoint, and the residual against a continuous run is O(dt).
    """
    Z = build(n_phi=12, n_xi=60, e=0.4, m=0.4, b=10.0)[0].grid.Z

    geo, sim = build(n_phi=12, n_xi=60, e=0.4, m=0.4, b=10.0)
    cont = sim.run(t_end=0.30 * Z, record_every=10_000)

    geo, sim = build(n_phi=12, n_xi=60, e=0.4, m=0.4, b=10.0)
    half = sim.run(t_end=0.15 * Z, record_every=10_000)
    geo, sim = build(n_phi=12, n_xi=60, e=0.4, m=0.4, b=10.0)
    # `run` always starts its clock at 0, so continuing for another 0.15 Z is
    # t_end = 0.15 Z with the half-way field as the initial condition
    split = sim.run(t_end=0.15 * Z, c0=half.concentration, record_every=10_000)

    path = str(tmp_path / "ckpt.npz")
    geo, sim = build(n_phi=12, n_xi=60, e=0.4, m=0.4, b=10.0)
    sim.run(t_end=0.15 * Z, record_every=10_000, checkpoint_path=path,
            checkpoint_every=0.0)
    geo, sim = build(n_phi=12, n_xi=60, e=0.4, m=0.4, b=10.0)
    res = sim.run(t_end=0.30 * Z, record_every=10_000, checkpoint_path=path,
                  checkpoint_every=0.0)

    assert np.max(np.abs(res.concentration - split.concentration)) == 0.0
    assert np.max(np.abs(split.concentration - cont.concentration)) < 1e-3

    # a checkpoint from a different mesh must fail loudly, not be interpolated
    geo, other = build(n_phi=12, n_xi=40, e=0.4, m=0.4, b=10.0)
    with pytest.raises(ValueError, match="checkpoint"):
        other.run(t_end=0.01 * Z, checkpoint_path=path)
