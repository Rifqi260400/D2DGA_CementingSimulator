"""M2 test gates -- gap-scale closure solver."""

import numpy as np
import pytest
import sympy as sp

import tests.benchmarks.kinematic_wave as kw

from d2dga.gapscale import newtonian as nwt
from d2dga.gapscale.al_solver import TwoLayerGapSolver
from d2dga.gapscale.closures import closures_from_solution, closures_newtonian_analytic
from d2dga.scaling import ScaledFluid

N_Y = 400          # chosen at M2-T5
TOL = 1e-10


def newtonian_pair(m):
    """BF25 Section 2.2: eta1 = m^(1/2), eta2 = m^(-1/2)."""
    return (ScaledFluid("displaced", 1.0, np.sqrt(m), 1.0, 0.0),
            ScaledFluid("displacing", 1.0, 1.0 / np.sqrt(m), 1.0, 0.0))


def F_b02(G, H, kappa, n, tau_y):
    """B02 (53): plane-Poiseuille pressure-drop / flow-rate law for a
    Herschel-Bulkley fluid.  H |s_bar| = F(G), with m = 1/n.  Independent of
    BF25, so it is a genuine external check on the AL solver."""
    m = 1.0 / n
    if H * G <= tau_y:
        return 0.0
    return ((H * G - tau_y) * ((m + 1) * H * G + tau_y) / (G ** 2 * (m + 1) * (m + 2))
            * ((H * G - tau_y) / kappa) ** m)


# ---------------------------------------------------------------- M2-T1 ---
@pytest.mark.parametrize("m", [0.2, 0.5, 1.0, 2.0, 5.0])
@pytest.mark.parametrize("c", [0.2, 0.5, 0.8])
def test_m2_t1_newtonian_isodense_matches_closed_form(m, c):
    """AL solve + quadrature must reproduce BF25 (2.24)-(2.26) and the
    corrected I3."""
    f1, f2 = newtonian_pair(m)
    solver = TwoLayerGapSolver(f1, f2, n_y=N_Y, tol=TOL)
    sol = solver.solve_fixed_mean_velocity(c, [0.0, 1.0], [0.0, 0.0])
    assert sol.converged
    num = closures_from_solution(sol)
    exact = closures_newtonian_analytic(c, m)
    assert np.isclose(num.I1, exact.I1, rtol=2e-5)
    assert np.isclose(num.I2, exact.I2, rtol=2e-5)
    assert np.isclose(num.q0, exact.q0, rtol=2e-5)
    assert np.isclose(num.I3, exact.I3, rtol=5e-4, atol=1e-9)


def test_m2_t1_analytic_forms_follow_from_the_integral_definitions():
    """BF25 (2.24)-(2.26) are not independent statements: they are what
    (2.14), (2.15), (2.22) give when eta1 = m^0.5 and eta2 = m^-0.5.
    Verified symbolically so no arithmetic slip can hide in the closed forms."""
    c, m, y = sp.symbols('c m y', positive=True)
    e1, e2 = sp.sqrt(m), 1 / sp.sqrt(m)
    A = sp.integrate(y ** 2 / e2, (y, 0, c))
    C = sp.integrate(y ** 2 / e1, (y, c, 1))
    Bq = sp.integrate(y / e1, (y, c, 1))
    D = sp.integrate((1 - y) / e1, (y, c, 1))
    E = sp.integrate(y * (1 - y) / e1, (y, c, 1))
    I1, I2 = A + C, (1 - c) * A + c * E
    q0 = (A + c * Bq) / I1
    I3 = ((A + c * Bq) * I2 - ((1 - c) * A + c ** 2 * D) * I1) / I1

    assert sp.simplify(I1 - (sp.sqrt(m) * c**3 + (1 - c**3) / sp.sqrt(m)) / 3) == 0
    assert sp.simplify(I2 - (2 * sp.sqrt(m) * c**3 * (1 - c)
                             + c * (1 - c)**2 * (1 + 2 * c) / sp.sqrt(m)) / 6) == 0
    assert sp.simplify(q0 - c * (m * c**2 + sp.Rational(3, 2) * (1 - c**2))
                       / (m * c**3 + (1 - c**3))) == 0
    master = -(c**2 * (1 - c)**3 * (4 * m * c + 3 * (1 - c))
               / (12 * sp.sqrt(m) * (m * c**3 + (1 - c**3))))
    assert sp.simplify(I3 - master) == 0
    # the same computation pins the layer convention: I1 and I2 above came out
    # equal to the PRINTED (2.24) and (2.25), so the sign of I3 is not a free
    # choice about which fluid sits at the wall -- it is fixed by the same
    # integrals that reproduce the paper's own two neighbouring formulae.
    assert master.subs({c: sp.Rational(1, 2), m: 1}) < 0


# ---------------------------------------------------------------- M2-T2 ---
@pytest.mark.parametrize("kappa,n,tau_y,G", [
    (1.0, 1.0, 0.0, 3.0),      # Newtonian -> plane Poiseuille
    (1.0, 1.0, 0.3, 3.0),      # Bingham
    (0.8, 0.6, 0.0, 4.0),      # power law, the case that exposed NUM-08
    (0.8, 0.6, 0.25, 4.0),     # Herschel-Bulkley
    (0.5, 0.4, 0.5, 6.0),      # strongly shear-thinning, high yield
    (1.2, 0.75, 0.45, 5.0),
])
def test_m2_t2_single_fluid_recovers_b02_poiseuille(kappa, n, tau_y, G):
    """One fluid filling the channel must reproduce B02 (53) exactly, and the
    central plug must occupy y < tau_Y/G."""
    f = ScaledFluid("single", 1.0, kappa, n, tau_y)
    solver = TwoLayerGapSolver(f, f, n_y=800, tol=1e-12, max_iter=200000)
    sol = solver.solve_fixed_pressure_gradient(0.0, [0.0, G], [0.0, 0.0])
    assert sol.converged
    assert np.isclose(sol.mean_velocity[1], F_b02(G, 1.0, kappa, n, tau_y), rtol=1e-5)
    assert np.isclose(sol.unyielded_fraction, tau_y / G, atol=2.0 / 800)


def test_m2_t2_newtonian_centreline_is_1p5_times_mean():
    """Plane Poiseuille: peak velocity = 1.5 x mean.  This is the number that
    sets the leading-edge dispersion speed in ZF23 (test M5-T4)."""
    f = ScaledFluid("single", 1.0, 1.0, 1.0, 0.0)
    solver = TwoLayerGapSolver(f, f, n_y=800, tol=1e-12)
    sol = solver.solve_fixed_mean_velocity(1.0, [0.0, 1.0], [0.0, 0.0])
    assert sol.converged
    assert np.isclose(sol.u[0, 1] / sol.mean_velocity[1], 1.5, rtol=1e-5)


# ---------------------------------------------------------------- M2-T3 ---
def test_m2_t3_cross_paper_reduction_and_conversion_factors():
    """
    THE cross-paper gate.  Settles CONV-03.

    BF25 (2.27) AS PRINTED does not follow from BF25's own (2.23).  Integrating
    (2.23) with eta1 = m^0.5, eta2 = m^-0.5 gives

        I3 = - c^2 (1-c)^3 [4 m c + 3(1-c)] / (12 m^(1/2) [m c^3 + 1 - c^3])

    differing from print in THREE places: 3(1-c) not 3(1-c^2); an m^(1/2) that
    is absent in (2.27) although (2.24)-(2.25) both carry their m^(+-1/2); and
    the overall SIGN.  The first two are transcription slips.  The third is the
    one that changes results -- see test_conv03_sign_is_what_makes_buoyancy
    _stabilising and newtonian.py's module docstring.

    Conversion factors, the deliverable this gate asks for.  The other papers
    print the same magnitude under their own sign conventions, so the factors
    carry the minus:
        BCF25 (23) form = -MASTER * m^(1/2)
        ZF22  (4.26)    = -MASTER * 6 / m^(1/2)
    """
    c, m = sp.symbols('c m', positive=True)
    master = -(c**2 * (1 - c)**3 * (4 * m * c + 3 * (1 - c))
               / (12 * sp.sqrt(m) * (m * c**3 + (1 - c**3))))
    bcf25 = c**2 * (1 - c)**3 * (4 * m * c + 3 * (1 - c)) / (12 * (m * c**3 + (1 - c**3)))
    zf22 = c**2 * (1 - c)**3 * (4 * m * c + 3 * (1 - c)) / (2 * m * (m * c**3 + (1 - c**3)))
    printed = (c**2 * (1 - c)**3 * (4 * m * c + 3 * (1 - c**2))
               / (12 * (m * c**3 + (1 - c**3))))

    assert sp.simplify(bcf25 / master + sp.sqrt(m)) == 0
    assert sp.simplify(zf22 / master + 6 / sp.sqrt(m)) == 0
    assert sp.simplify(printed - master) != 0          # the discrepancy is real

    # numeric confirmation through the module functions
    cc = np.linspace(0.05, 0.95, 19)
    for mm in (0.2, 1.0, 5.0):
        assert np.allclose(nwt.script_I3_bcf25_form(cc, mm),
                           -nwt.script_I3(cc, mm) * np.sqrt(mm), rtol=1e-12)
        assert np.allclose(nwt.script_I3_zf22_form(cc, mm),
                           -nwt.script_I3(cc, mm) * 6 / np.sqrt(mm), rtol=1e-12)


def test_conv03_sign_is_what_makes_buoyancy_stabilising():
    """
    The sign of I3 is the model's central physical claim, so test the claim and
    not just the algebra.

    BF25 Section 3.2, verbatim: "the shock velocity approaches the mean velocity
    (w_f -> 1), and the front becomes more stable when b >> 1 ... There is less
    dispersion for larger b", and "buoyancy completely suppresses the spike
    regime as b increases from 100 to 1000".

    The object those sentences describe is the SHOCK, not the tip.  The tip
    speed is 1.5 for every b (BF25 says so in the same paragraph: q0'(0) = 1.5
    and I3'(0) = 0), so it cannot be what "approaches the mean velocity".  The
    shock comes out of the upper concave envelope of f(c) = q0 + b I3 on [0, 1]
    -- see tests/benchmarks/kinematic_wave.py -- and its speed is the envelope
    segment's slope, BF25 (3.9).
    """
    for m in (0.2, 1.0, 5.0):
        speeds, spikes = [], []
        for b in (10.0, 100.0, 1000.0):
            shock = kw.riemann_structure(m, b)["main_shock"]
            assert shock is not None, (m, b)
            speeds.append(shock[2])
            spikes.append(shock[0])        # the spike occupies 0 < c < c_lo
        assert speeds == sorted(speeds, reverse=True), (m, speeds)
        assert abs(speeds[-1] - 1.0) < 0.05, (m, speeds)        # -> mean speed
        assert spikes == sorted(spikes, reverse=True), (m, spikes)
        assert spikes[-1] < 0.1 * spikes[0], (m, spikes)        # suppressed

        # reversing the sign turns the same statement inside out: the flux
        # grows a hump many times the mean flux and the leading wave runs away
        assert kw.riemann_structure(m, -1000.0)["leading_speed"] > 20.0

    # BF25 Section 3.2 again: "the effect of the buoyancy number and the
    # viscosity ratio on q0' + b I3' is negligible at c = 0, and it maintains a
    # constant value of 1.5".  I3 ~ c^2 near zero, so this is exact.
    for m in (0.2, 1.0, 5.0):
        for b in (0.0, 100.0, 1000.0):
            speed = (nwt.dq0_dc(np.array(0.0), m)
                     + b * nwt.dscript_I3_dc(np.array(0.0), m))
            assert abs(float(speed) - 1.5) < 1e-12


@pytest.mark.parametrize("m", [0.2, 1.0, 5.0])
def test_m2_t3_al_solver_confirms_master_not_printed(m):
    """The numerical solver is an independent third route.  It must agree with
    MASTER and disagree with BF25 (2.27) as printed."""
    f1, f2 = newtonian_pair(m)
    solver = TwoLayerGapSolver(f1, f2, n_y=N_Y, tol=TOL)
    for c in (0.3, 0.6):
        sol = solver.solve_fixed_mean_velocity(c, [0.0, 1.0], [0.0, 0.0])
        num = closures_from_solution(sol).I3
        assert np.isclose(num, nwt.script_I3(c, m), rtol=5e-4, atol=1e-10)
        if not np.isclose(m, 1.0):
            assert not np.isclose(num, nwt.script_I3_bf25_as_printed(c, m), rtol=0.05)


# ---------------------------------------------------------------- M2-T4 ---
@pytest.mark.parametrize("kappa1,n1,ty1,kappa2,n2,ty2", [
    (1.0, 1.0, 0.0, 1.0, 1.0, 0.0),
    (2.0, 0.6, 0.3, 0.5, 0.8, 0.1),
    (0.5, 0.4, 0.9, 1.5, 0.7, 0.2),
])
def test_m2_t4_q0_endpoints_are_exact(kappa1, n1, ty1, kappa2, n2, ty2):
    """q0(0) = 0 and q0(1) = 1 EXACTLY, for every rheology.  These are the
    minimum requirement BCF25 places on any closure: they are what make
    sum_k q_k0 = 1 hold, and hence the transport scheme conservative."""
    f1 = ScaledFluid("d", 1.0, kappa1, n1, ty1)
    f2 = ScaledFluid("g", 1.0, kappa2, n2, ty2)
    solver = TwoLayerGapSolver(f1, f2, n_y=200, tol=1e-9, max_iter=60000)
    for c, expected in ((0.0, 0.0), (1.0, 1.0)):
        sol = solver.solve_fixed_mean_velocity(c, [0.0, 1.0], [0.0, 0.0])
        assert closures_from_solution(sol).q0 == expected


def test_m2_t4_analytic_q0_endpoints_are_exact():
    for m in (0.2, 0.3, 1.0, 5.0, 17.3):
        assert nwt.q0(0.0, m) == 0.0
        assert nwt.q0(1.0, m) == 1.0
        assert nwt.script_I2(0.0, m) == 0.0
        assert nwt.script_I2(1.0, m) == 0.0
        assert nwt.script_I3(0.0, m) == 0.0
        assert nwt.script_I3(1.0, m) == 0.0
        assert np.all(nwt.script_I3(np.linspace(0.05, 0.95, 19), m) < 0.0)


# ---------------------------------------------------------------- M2-T5 ---
def test_m2_t5_second_order_convergence_in_n_y():
    """
    BF25 A.2.3 reports second-order integration, so tolerances below 1/N_y^2
    buy nothing.  Confirm the order, then justify N_y.

    This gate previously FAILED at first order.  Cause: BF25 (A18) fixes the
    constant of integration with -r q(0) + lam~(0), and evaluating those by
    extrapolating cell-centred data is O(dy) wrong for du/dy ~ y^(1/n).  Both
    fields are odd about the channel centre, so the bracket is exactly zero.
    See al_solver.py and assumptions.md NUM-08.
    """
    kappa, n, tau_y, G = 0.8, 0.6, 0.0, 4.0     # the worst case: fractional
    exact = F_b02(G, 1.0, kappa, n, tau_y)      # power, no plug to smooth it
    f = ScaledFluid("pl", 1.0, kappa, n, tau_y)
    errors = []
    for n_y in (100, 200, 400, 800):
        solver = TwoLayerGapSolver(f, f, n_y=n_y, tol=1e-12, max_iter=200000)
        sol = solver.solve_fixed_pressure_gradient(0.0, [0.0, G], [0.0, 0.0])
        errors.append(abs(sol.mean_velocity[1] - exact) / exact)
    for coarse, fine in zip(errors[:-1], errors[1:]):
        assert 3.5 < coarse / fine < 4.5, f"order ratio {coarse/fine:.2f}, want ~4"
    # N_Y = 400 is the production choice: ~2.5e-6, comfortably inside M3-T1's
    # 0.1% interpolation budget, at a quarter the cost of 800.
    assert errors[2] < 1e-5


def test_m2_t5_tolerance_floor_is_the_mesh_not_the_iteration():
    """Below ~1/N_y^2 the tolerance stops mattering, exactly as BF25 says.
    Tightening from 1e-8 to 1e-14 must not move the answer."""
    kappa, n, tau_y, G = 0.8, 0.6, 0.0, 4.0
    f = ScaledFluid("pl", 1.0, kappa, n, tau_y)
    results = []
    for tol in (1e-8, 1e-11, 1e-14):
        solver = TwoLayerGapSolver(f, f, n_y=400, tol=tol, max_iter=200000)
        results.append(solver.solve_fixed_pressure_gradient(
            0.0, [0.0, G], [0.0, 0.0]).mean_velocity[1])
    # 1e-7 is the right bar, not 1e-9: the iteration floor sits at ~6e-9
    # relative, which is already 400x below the 2.5e-6 MESH error at N_y = 400.
    # That gap is exactly the claim -- the mesh, not the tolerance, is the
    # binding constraint.
    assert np.allclose(results, results[0], rtol=1e-7)


# ------------------------------------------------- two-layer specifics ----
def test_buoyancy_rotates_the_stress_direction_across_the_interface():
    """BF25 Fig. 4: with Gb != 0 the shear-stress direction JUMPS at the
    interface, because buoyancy acts oppositely in the two layers.  This is why
    the gap problem must be solved as a 2-vector and cannot be reduced to a
    scalar per direction."""
    f1 = ScaledFluid("d", 1.0, 1.0, 1.0, 0.0)
    f2 = ScaledFluid("g", 1.0, 1.0, 1.0, 0.0)
    solver = TwoLayerGapSolver(f1, f2, n_y=400, tol=1e-11)

    def stress_gradient_directions(sol):
        kappa, n, tau_y = solver._props(sol.in_fluid2)
        g = np.linalg.norm(sol.dudy, axis=1)
        safe = np.where(g > 0, g, 1.0)
        tau = (kappa * g ** n + tau_y)[:, None] * sol.dudy / safe[:, None]
        out = {}
        for key, mask in (("inner", sol.in_fluid2), ("outer", ~sol.in_fluid2)):
            d = np.diff(tau[mask], axis=0)
            out[key] = d / np.linalg.norm(d, axis=1)[:, None]
        return out

    # BF25 (2.9)-(2.10): d tau/dy is a CONSTANT vector in each layer,
    #     layer 2: -G + (1-c) Gb        layer 1: -(G + c Gb)
    # so in the (tau_phi, tau_xi) plane the stress traces two straight segments.
    # They are not collinear -- they differ by Gb -- which is the kink in
    # BF25 Fig. 4(a).  Note the segments' DIRECTIONS differ; the stress
    # direction itself rotates through the outer layer, because that segment
    # does not pass through the origin.
    d = stress_gradient_directions(
        solver.solve_fixed_mean_velocity(0.5, [0.0, 1.0], [1.0, 1.0]))
    assert np.max(np.ptp(d["inner"], axis=0)) < 1e-9
    assert np.max(np.ptp(d["outer"], axis=0)) < 1e-9
    assert np.linalg.norm(d["inner"][0] - d["outer"][0]) > 1e-3

    # with Gb = 0 the two constants coincide and the kink disappears
    d0 = stress_gradient_directions(
        solver.solve_fixed_mean_velocity(0.5, [0.0, 1.0], [0.0, 0.0]))
    assert np.linalg.norm(d0["inner"][0] - d0["outer"][0]) < 1e-9


def test_interface_lands_exactly_on_a_cell_face():
    f1 = ScaledFluid("d", 1.0, 1.0, 1.0, 0.0)
    f2 = ScaledFluid("g", 1.0, 2.0, 1.0, 0.0)
    solver = TwoLayerGapSolver(f1, f2, n_y=100)
    for c in (0.01, 0.13, 0.5, 0.77, 0.99):
        sol = solver.solve_fixed_mean_velocity(c, [0.0, 1.0], [0.0, 0.0])
        assert np.min(np.abs(sol.y_edges - c)) < 1e-14


# ------------------------------------------------- NUM-02 / Q3 ------------
def test_num02_al_parameter_sets_the_rate_and_never_the_answer():
    """
    NUM-02 (Q3).  PF04 gives only 0 < rho < r(1+sqrt5)/2 and BF25 A.2.3 reports
    "rho = r = 1"; neither says how to choose r.  Three things have to hold for
    `TwoLayerGapSolver.tuned` to be a legitimate optimisation rather than a
    tuning knob that moves results:

      1. the converged closures must be INDEPENDENT of r -- to within the
         convergence tolerance, which is what that phrase actually means, not
         to machine precision;
      2. the iteration count must really depend on it, or there is nothing to
         tune;
      3. the best r must be FLUID-PAIR dependent, or a new fixed default would
         do and probing would be waste.

    Point 3 is the one that settles the design, and it is easy to get wrong:
    for the K-GEP-1 pair r = 1 costs 3434 iterations against 56 at r = 0.01, a
    factor of 61 -- but for the pair below r = 1 is already near optimal and
    r = 0.01 is ~20x WORSE.  A single new default would therefore have traded
    one badly scaled case for another.
    """
    from d2dga.config import FluidsConfig
    from d2dga.geometry import build_geometry
    from d2dga.config import Config
    from d2dga.scaling import Scaling

    cfg = Config()
    geo = build_geometry(cfg)
    mud, cement = cfg.fluids.as_fluids()
    sc = Scaling(mud, cement, r_a_hat_star=geo.r_a_hat_star,
                 delta_star=geo.delta_star, mean_velocity=0.2)
    f1, f2, gb = sc.scaled_fluid1, sc.scaled_fluid2, sc.buoyancy_number
    assert gb > 20.0                       # strongly buoyant, where r matters

    got = {}
    for r in (1.0, 0.1, 0.01):
        solver = TwoLayerGapSolver(f1, f2, n_y=200, r=r, rho=r, tol=1e-9,
                                   max_iter=200000)
        sol = solver.solve_fixed_mean_velocity(0.5, [0.0, 1.0], [0.0, gb])
        assert sol.converged, r
        got[r] = (closures_from_solution(sol), sol.iterations)

    # (1) r-independence, to the convergence tolerance
    ref = got[1.0][0]
    for r, (cl, _) in got.items():
        assert cl.I1 == pytest.approx(ref.I1, rel=1e-7), r
        assert cl.I3 == pytest.approx(ref.I3, rel=1e-6, abs=1e-12), r

    # (2) the cost really moves -- and in this direction for this pair
    assert got[1.0][1] > 20 * got[0.01][1], {k: v[1] for k, v in got.items()}

    # (3) the SAME r is far from optimal for a different pair, so a fixed
    # default cannot serve both
    g1, g2 = newtonian_pair(0.2)
    g2 = ScaledFluid("displacing", 1.0, g2.consistency, 0.5, 0.8)
    other = {}
    for r in (1.0, 0.01):
        solver = TwoLayerGapSolver(g1, g2, n_y=200, r=r, rho=r, tol=1e-9,
                                   max_iter=200000)
        sol = solver.solve_fixed_mean_velocity(0.5, [0.0, 1.0], [0.0, 25.0])
        assert sol.converged, r
        other[r] = sol.iterations
    assert other[0.01] > 5 * other[1.0], other

    # (4) probing picks a rung at least as good as the default, on both pairs
    for (a, bb, g, worst) in ((f1, f2, gb, got[1.0][1]),
                              (g1, g2, 25.0, other[1.0])):
        tuned = TwoLayerGapSolver.tuned(a, bb, gb_probe=g, n_y=200, tol=1e-9,
                                        max_iter=200000)
        sol = tuned.solve_fixed_mean_velocity(0.5, [0.0, 1.0], [0.0, g])
        assert sol.converged and sol.iterations <= worst


# =========================================================================
# B-4 / CONV-03 -- the cross-paper check that settles the duplicate register
# rows, done on the PHYSICAL FLUX rather than on anyone's printed I3.
# =========================================================================
def test_b4_conv03_buoyant_flux_agrees_with_zf22_4_25():
    """
    B-4.  `docs/assumptions.md` carried two CONV-03 rows that disagreed on the
    SIGN of the conversions to ZF22 (4.26) and BCF25 (23).  This settles it
    from the primary sources, on the quantity that actually enters the model.

    BF25 (2.21), xi-component, r_a = 1, beta = 0:

        q_buoy = b~ H^3 I3^BF25(c, m) ,        b~ = sqrt(m) b_ZF22   (CONV-02)

    ZF22 (4.25), same component:

        q_buoy = [Delta_rho / F^2] (H^3 / 6 eta_2) I3^ZF22(c, m)

    read with ZF22's OWN conventions, both stated in its own text: at (4.24)
    ZF22 writes `rho = rho_2 + (1 - c) Delta_rho` against `rho = (1-c)rho_1 +
    c rho_2`, i.e. **Delta_rho = rho_1 - rho_2**, the same ordering as BF25 --
    so Delta_rho/F^2 = -b_ZF22, NOT +b_ZF22; and eta_2 = mu_2/mu_1 = 1/m.

    With those readings the two papers agree to machine precision, and they
    agree on the NEGATIVE I3 that this code uses.  Reading ZF22's Delta_rho
    with the opposite sign flips the result, which is precisely the mistake the
    stale register row encoded.
    """
    from d2dga.scaling import zf22_buoyancy_to_bf25

    def I3_zf22(c, m):                      # ZF22 (4.26), verbatim
        return (c ** 2 * (1 - c) ** 3 * (4 * m * c + 3 * (1 - c))
                / (2 * m * (m * c ** 3 + 1 - c ** 3)))

    worst = 0.0
    for m in (0.2, 0.5, 2.0, 5.0):
        for b_zf in (-50.0, 10.0, 100.0, 1000.0):
            for c in (0.2, 0.5, 0.8):
                for H in (0.4, 1.0, 1.7):
                    ours = zf22_buoyancy_to_bf25(b_zf, m) * H ** 3 \
                        * nwt.script_I3(c, m)
                    theirs = -b_zf * (m * H ** 3 / 6.0) * I3_zf22(c, m)
                    worst = max(worst, abs(ours - theirs)
                                / max(abs(theirs), 1e-300))
    assert worst < 1e-13, worst


def test_b4_conv03_conversion_factors_carry_their_minus_signs():
    """
    B-4.  The two stale-row conversions, asserted with their signs.  Against the
    register row that omitted them this fails by a factor of -1.

    Both published forms are printed POSITIVE on 0 < c < 1 while the derived
    master is negative, so both conversions carry a minus:

        BCF25 (23)  = -sqrt(m)   * MASTER
        ZF22 (4.26) = -(6/sqrt(m)) * MASTER

    Note also that BCF25 (23) and BF25 (2.27) are the same printed formula
    except for the numerator -- BCF25 has 3(1-c), BF25 has 3(1-c^2) -- so ZF22
    and BCF25 between them confirm, independently of this code, that BF25
    (2.27)'s 3(1-c^2) is a transcription error.
    """
    def I3_bcf25(c, m):                     # BCF25 (23), verbatim
        return (c ** 2 * (1 - c) ** 3 * (4 * m * c + 3 * (1 - c))
                / (12 * (m * c ** 3 + 1 - c ** 3)))

    def I3_zf22(c, m):                      # ZF22 (4.26), verbatim
        return (c ** 2 * (1 - c) ** 3 * (4 * m * c + 3 * (1 - c))
                / (2 * m * (m * c ** 3 + 1 - c ** 3)))

    c = np.linspace(0.02, 0.98, 97)
    for m in (0.006, 0.2, 1.0, 5.0, 160.0):
        master = nwt.script_I3(c, m)
        assert np.all(master < 0.0)
        assert np.allclose(I3_bcf25(c, m), -np.sqrt(m) * master, rtol=1e-13)
        assert np.allclose(I3_zf22(c, m), -(6 / np.sqrt(m)) * master, rtol=1e-13)
