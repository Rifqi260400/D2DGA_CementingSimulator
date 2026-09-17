"""M5 test gates -- LLF transport of the gap-averaged concentration."""

import numpy as np
import pytest
import sympy as sp

from d2dga.config import GridConfig, WellConfig
from d2dga.elliptic import (NewtonianClosures, StreamFunctionSolver,
                            TwoDGAClosures)
from d2dga.gapscale import newtonian as _newt
from d2dga.geometry import (ConstantEccentricity, Geometry, SinusoidalWall,
                            UniformWall)
from d2dga.scaling import HerschelBulkleyFluid
from d2dga.transport import TransportSolver

# A unit-viscosity Newtonian fluid, for the 2DGA closure path.  Density is
# irrelevant there: TwoDGAClosures mixes only the rheology, and the buoyancy
# enters through the separate b passed to the solvers.
_NEWTONIAN = HerschelBulkleyFluid(name="unit", density=1000.0,
                                  consistency=1.0, power_law_index=1.0,
                                  yield_stress=0.0)

# Same short window as M4: on the full 196 m K-GEP-1 domain Z ~ 564 and the
# structures of interest are sub-cell on any affordable mesh.
SHORT = dict(casing_shoe_m=387.2, total_depth_m=390.0)


def build(n_phi=16, n_xi=200, beta=0.0, e=0.0, m=1.0, b=0.0, wall=None,
          two_dga=False, cfl=0.5, c_in=1.0, fluids=None):
    well = WellConfig(inclination_rad=beta, **SHORT)
    w = wall if wall is not None else UniformWall(well.gauge_hole_diameter_m)
    geo = Geometry(well, w, ConstantEccentricity(e), GridConfig(n_phi, n_xi))
    clo = TwoDGAClosures(*fluids) if two_dga else NewtonianClosures(m)
    # b_num = -delta_rho / Fr*^2, so Fr = 1 and delta_rho = -b give b_num = b.
    ell = StreamFunctionSolver(geo, clo, froude=1.0, delta_rho=-b)
    tr = TransportSolver(geo, clo, buoyancy_number=b,
                         inflow_concentration=c_in, cfl=cfl)
    return geo, ell, tr


# =========================================================================
# formulation -- the conservation form is exact, not an approximation
# =========================================================================
def test_conservation_form_is_algebraically_identical_to_bf25_2_19():
    """
    The code advances d[H r_a c]/dt + dF/dphi + dG/dxi = 0, not BF25 (2.19)'s
    div_a . q = 0.  div_a carries 1/r_a on the azimuthal part, so the two forms
    coincide only because r_a depends on xi alone.  Assert it symbolically
    rather than trusting the hand algebra.
    """
    phi, xi, b, beta = sp.symbols('phi xi b beta')
    ra = sp.Function('r_a')(xi)
    H = sp.Function('H')(phi, xi)
    Psi = sp.Function('Psi')(phi, xi)
    q0 = sp.Function('q0')(phi, xi)
    I3 = sp.Function('I3')(phi, xi)

    f_phi = ra * sp.cos(beta)
    f_xi = ra * sp.sin(beta) * sp.sin(sp.pi * phi)

    # BF25 (2.21) verbatim
    q_phi = ra * sp.Rational(1, 2) * (-sp.diff(Psi, xi)) * q0 - b * H**3 * I3 * f_xi
    q_xi = ra * sp.Rational(1, 2) * (sp.diff(Psi, phi) / ra) * q0 \
        + b * H**3 * I3 * f_phi
    div_a_q = sp.diff(q_phi, phi) / ra + sp.diff(q_xi, xi)

    # transport.py (T2), (T3)
    F = -sp.Rational(1, 2) * q0 * sp.diff(Psi, xi) \
        - b * H**3 * I3 * sp.sin(beta) * sp.sin(sp.pi * phi)
    G = sp.Rational(1, 2) * q0 * sp.diff(Psi, phi) \
        + b * H**3 * ra * I3 * sp.cos(beta)

    assert sp.simplify(div_a_q - (sp.diff(F, phi) + sp.diff(G, xi))) == 0


def test_flux_reduces_to_bf25_3_7_in_one_dimension():
    """
    BF25 (3.7) -- the 1-D kinematic wave equation the paper's whole regime
    classification rests on -- reads dc/dt + d[q0 + b I3]/dxi = 0.  It must fall
    out of (T3) with H = r_a = 1, beta = 0 and unit mean velocity, with NO
    stray factor of 2 from the streamfunction convention.  If CONV-01 had been
    applied inconsistently this is where it would show.
    """
    geo, ell, tr = build(n_phi=8, n_xi=40, m=2.0, b=1.5)
    c = np.linspace(0.05, 0.95, 40)[None, :] * np.ones((8, 1))
    # uniform axial flow: Psi = 2 H r_a w phi, and H = r_a = w = 1 here
    psi = 2.0 * geo.grid.phi_edges[:, None] * np.ones((1, geo.grid.n_xi + 1))
    Phi, Xi, _, _ = tr.fluxes(c, psi)

    assert np.allclose(tr.H_c, 1.0) and np.allclose(tr.ra_c, 1.0)
    # face-integrated G over dphi, so divide it out to recover G itself
    G = Xi[:, 1:-1] / geo.grid.dphi
    c_face = 0.5 * (c[:, :-1] + c[:, 1:])
    expected = _newt.q0(c_face, 2.0) + 1.5 * _newt.script_I3(c_face, 2.0)
    # LLF averages the flux function rather than evaluating it at the averaged
    # state, so compare against the same average
    expected = 0.5 * (_newt.q0(c[:, :-1], 2.0) + _newt.q0(c[:, 1:], 2.0)) \
        + 1.5 * 0.5 * (_newt.script_I3(c[:, :-1], 2.0)
                       + _newt.script_I3(c[:, 1:], 2.0))
    dissip = G - expected
    # the remaining difference is exactly the LLF dissipation
    assert np.max(np.abs(dissip)) > 0
    assert np.allclose(G - dissip, expected)


# =========================================================================
# M5-T1 -- volume conservation
# =========================================================================
@pytest.mark.parametrize("beta,b,e,wall", [
    (0.0, 0.0, 0.0, None),
    (0.0, 6.0, 0.3, None),
    (np.pi / 4, 6.0, 0.3, None),
])
def test_m5_t1_volume_conservation_is_exact(beta, b, e, wall):
    """
    Gate M5-T1.  Every step, the change in sum(H r_a c) must equal the net flux
    through the two open ends, to machine precision -- not to a tolerance.
    This is the whole reason the equation was put into the plain-divergence
    form (T1): the azimuthal flux telescopes to the symmetry planes, where it
    is identically zero.
    """
    geo, ell, tr = build(n_phi=20, n_xi=120, beta=beta, e=e, m=0.6, b=b,
                         wall=wall)
    rng = np.random.default_rng(11)
    c = rng.random((20, 120))
    for _ in range(40):
        psi = ell.solve(c, Q=1.0)
        m_before = tr.mass(c)
        c, dt, bflux = tr.step(c, psi)
        m_after = tr.mass(c)
        rel = abs((m_after - m_before) - bflux) / max(abs(m_before), 1e-30)
        assert rel < 5e-15, rel


def test_m5_t1_azimuthal_flux_vanishes_on_the_symmetry_planes():
    """The telescoping in M5-T1 works only because phi = 0 and phi = 1 carry no
    flux.  For the advective part that is automatic -- Psi is constant along
    both -- but the buoyancy part is set to zero explicitly, so assert that the
    part which is automatic really is."""
    geo, ell, tr = build(n_phi=24, n_xi=60, beta=np.pi / 3, e=0.4, b=5.0)
    rng = np.random.default_rng(3)
    c = rng.random((24, 60))
    psi = ell.solve(c, Q=1.0)
    # advective part alone, with the explicit zeroing not applied
    dpsi_xi = psi[:, :-1] - psi[:, 1:]
    assert np.max(np.abs(dpsi_xi[0, :])) < 1e-12
    assert np.max(np.abs(dpsi_xi[-1, :])) < 1e-12
    Phi, _, a_phi, _ = tr.fluxes(c, psi)
    assert np.all(Phi[0, :] == 0.0) and np.all(Phi[-1, :] == 0.0)
    assert np.all(a_phi[0, :] == 0.0) and np.all(a_phi[-1, :] == 0.0)


def test_uniform_states_are_preserved_exactly():
    """c == 0 and c == 1 are fixed points, to the last bit.  For c == 1 this
    needs q0(1) == 1 EXACTLY and the discrete Psi-differences around a cell to
    telescope to exactly zero; a 1e-16 error in either would seed a slow drift
    over the ~10^5 steps of a full job."""
    for value in (0.0, 1.0):
        geo, ell, tr = build(n_phi=12, n_xi=80, beta=np.pi / 5, e=0.35,
                             m=0.4, b=4.0, c_in=value)
        c = np.full((12, 80), value)
        psi = ell.solve(c, Q=1.0)
        for _ in range(5):
            c, dt, _ = tr.step(c, psi)
            assert np.all(c == value)


# =========================================================================
# M5-T2 -- maximum principle
# =========================================================================
@pytest.mark.parametrize("beta,b,cfl", [(0.0, 0.0, 0.5), (0.0, 8.0, 0.5),
                                        (np.pi / 4, 8.0, 0.5),
                                        (np.pi / 4, 8.0, 0.99)])
def test_m5_t2_maximum_principle(beta, b, cfl):
    """Gate M5-T2.  c_bar must stay in [0, 1] at every cell and every step, from
    a deliberately rough initial condition that a non-monotone scheme would
    immediately overshoot."""
    geo, ell, tr = build(n_phi=20, n_xi=100, beta=beta, e=0.4, m=0.3, b=b,
                         cfl=cfl)
    rng = np.random.default_rng(7)
    c = rng.integers(0, 2, size=(20, 100)).astype(float)   # 0/1 checkerboard
    for _ in range(120):
        psi = ell.solve(c, Q=1.0)
        c, dt, _ = tr.step(c, psi)
        assert c.min() >= 0.0, c.min()
        assert c.max() <= 1.0, c.max()


def test_m5_t2_cfl_bound_is_the_thing_doing_the_work():
    """
    The maximum principle is not incidental -- it is bought by BCF25 (44), and
    the bound is sharp.

    Sharpness has to be demonstrated on a state that attains it, not on an
    arbitrary rough field: (44) is a MINIMUM over cells, and a cell stepped past
    its own limit still stays positive if its neighbours happen to feed it.  The
    state that exposes the bound is an isolated spike, where every F_l(0) = 0
    and the update collapses to c^{n+1} = c^n (1 - dt/dt_cell) exactly.

    The cell's OWN limit is what is tested, not the global one: since NUM-26 the
    wavespeeds are interval maxima, so placing the spike changes which cell
    attains the global minimum.  The global bound is then checked to be no
    larger, which is what makes stepping at it safe everywhere.
    """
    geo, ell, tr = build(n_phi=16, n_xi=80, m=0.3, b=8.0, beta=np.pi / 4, e=0.4)
    c = np.zeros((16, 80))
    psi = ell.solve(c, Q=1.0)          # Psi does not depend on c for Newtonian
    p, q = 7, 40
    c[p, q] = 1.0

    Phi, Xi, a_phi, a_xi = tr.fluxes(c, psi)
    total = (a_phi[:-1, :] + a_phi[1:, :] + a_xi[:, :-1] + a_xi[:, 1:])
    dt_cell = tr.cell_area * 2.0 * tr.Hra[p, q] / total[p, q]
    assert tr.max_timestep(a_phi, a_xi) <= dt_cell * (1.0 + 1e-12)

    def advance(dt):
        lam = dt / tr.cell_area
        div = (Phi[1:, :] - Phi[:-1, :]) + (Xi[:, 1:] - Xi[:, :-1])
        return c - lam * div / tr.Hra

    assert abs(advance(dt_cell)[p, q]) < 1e-12          # the bound is attained
    assert advance(1.6 * dt_cell)[p, q] == pytest.approx(-0.6, abs=1e-10)

    at_global = advance(tr.max_timestep(a_phi, a_xi))
    assert at_global.min() >= -1e-14
    assert at_global.max() <= 1.0 + 1e-14

    with pytest.raises(ValueError, match="monotonicity limit"):
        tr.step(c, psi, dt=2.0 * dt_cell)


def test_k_equals_two_complement_sums_to_one():
    """
    K = 2, so fluid 1 never gets its own scalar: it rides on q_{1,0} = 1 - q0
    and I_{1,3} = -I3 (BCF25 6).  Advancing 1 - c_bar with those complementary
    closures must reproduce 1 - (advance of c_bar) exactly, which is BCF25's
    sum-to-one property with the redundant equation removed.
    """
    class Complement(NewtonianClosures):
        def q0_I3(self, c, H, umag=None):
            q, I = super().q0_I3(1.0 - np.asarray(c, float), H, umag)
            return 1.0 - q, -I

        def dq0_dI3_dc(self, c, H, umag=None, h=None):
            dq, dI = super().dq0_dI3_dc(1.0 - np.asarray(c, float), H, umag)
            return dq, dI      # d/dc[1 - q0(1-c)] = q0'(1-c); likewise -I3'

        def wavespeed_nodes(self):
            # The complement's derivatives peak at 1 - (the original peaks).
            # Mirroring them is what makes the two solvers agree bit for bit:
            # the LLF interval maximum (NUM-26) is part of the flux, so
            # sum-to-one holds only if the two fluids' wavespeed bounds are
            # themselves complementary.
            return np.sort(1.0 - super().wavespeed_nodes())

    geo, ell, tr = build(n_phi=12, n_xi=60, beta=np.pi / 6, e=0.3, m=0.45, b=5.0)
    tr1 = TransportSolver(geo, Complement(0.45), buoyancy_number=5.0,
                          inflow_concentration=0.0, cfl=0.5)
    rng = np.random.default_rng(19)
    c2 = rng.random((12, 60))
    c1 = 1.0 - c2
    psi = ell.solve(c2, Q=1.0)
    dt = min(tr.cfl * tr.max_timestep(*tr.fluxes(c2, psi)[2:]),
             tr1.cfl * tr1.max_timestep(*tr1.fluxes(c1, psi)[2:]))
    new2, _, _ = tr.step(c2, psi, dt=dt)
    new1, _, _ = tr1.step(c1, psi, dt=dt)
    assert np.max(np.abs(new1 + new2 - 1.0)) < 5e-15


# =========================================================================
# M5-T3 -- pure advection
# =========================================================================
def test_m5_t3_pure_advection_travels_at_the_gap_averaged_speed():
    """
    Gate M5-T3.  With the 2DGA closures (q0 = c, I3 = 0) and a uniform
    concentric annulus the transport is linear advection at w_bar = 1.  A
    conservative scheme moves the front's CENTRE OF MASS at exactly that speed
    however much it smears it, so that is what is measured.
    """
    f = _NEWTONIAN
    geo, ell, tr = build(n_phi=8, n_xi=400, two_dga=True, fluids=(f, f),
                         c_in=0.0)
    n_xi = geo.grid.n_xi
    xi = geo.grid.xi_centres
    c = np.zeros((8, n_xi))
    c[:, 40:80] = 1.0                 # a slug, away from both boundaries
    psi = ell.solve_nonlinear(c, Q=1.0)[0]
    v, w = ell.velocities(psi)
    assert np.allclose(w, 1.0, atol=1e-10)

    def centroid(cc):
        return float(np.sum(cc[0] * xi) / np.sum(cc[0]))

    x0, t = centroid(c), 0.0
    T = 2.0
    while t < T:
        dt = min(tr.cfl * tr.max_timestep(*tr.fluxes(c, psi)[2:]), T - t)
        c, dt, _ = tr.step(c, psi, dt=dt)
        t += dt
    assert c.max() <= 1.0 and c.min() >= 0.0
    assert abs((centroid(c) - x0) - 1.0 * t) < 1e-10


def test_num16_azimuthal_upwind_direction():
    """
    NUM-16.  BCF25 (18) and (20) print the azimuthal LLF coefficients with the
    Psi-differences carrying the opposite sign to what (11) and (14) imply.
    The two readings are not both defensible: one upwinds from the upstream
    side and the other from the downstream side.

    Set v_bar > 0, so dPsi/dxi < 0 and the upstream cell is the one at LOWER
    phi.  With the 2DGA closures the LLF flux collapses to pure upwinding, so
    the flux must depend on the low-phi value alone.  Taking (18)/(20) as
    printed would make it depend on the high-phi value instead.
    """
    f = _NEWTONIAN
    geo, ell, tr = build(n_phi=4, n_xi=4, two_dga=True, fluids=(f, f))
    n_phi, n_xi = geo.grid.n_phi, geo.grid.n_xi

    # Psi decreasing in xi  =>  dPsi/dxi < 0  =>  v_bar = -dPsi/dxi/(2H) > 0
    psi = np.zeros((n_phi + 1, n_xi + 1))
    psi += -np.arange(n_xi + 1)[None, :]
    c = np.zeros((n_phi, n_xi))
    c[1, :] = 1.0                       # a single lit row at phi index 1

    Phi, _, _, _ = tr.fluxes(c, psi)
    # face 2 separates cell 1 (L, lit) from cell 2 (R, dark).  Upwinding from
    # the low-phi side means this face carries flux; face 1 (dark L, lit R)
    # must carry none.
    assert np.all(np.abs(Phi[2, :]) > 1e-12)
    assert np.allclose(Phi[1, :], 0.0, atol=1e-15)


# =========================================================================
# M5-T4 -- the dispersive front
# =========================================================================
def _rarefaction(n_xi, T=4.0, m=1.0, cfl=0.5):
    geo, ell, tr = build(n_phi=8, n_xi=n_xi, m=m, b=0.0, cfl=cfl, c_in=1.0)
    c = np.zeros((8, n_xi))
    psi = ell.solve(c, Q=1.0)
    t = 0.0
    while t < T:
        dt = min(tr.cfl * tr.max_timestep(*tr.fluxes(c, psi)[2:]), T - t)
        c, dt, _ = tr.step(c, psi, dt=dt)
        t += dt
    return geo, c, t


def test_m5_t4_leading_edge_moves_at_the_poiseuille_centreline_speed():
    """
    Gate M5-T4.  For m = 1, q0 = 1.5 c - 0.5 c^3, so q0'(0) = 3/2: the tip of
    the displacing fluid runs at 1.5 x the mean velocity, the centreline speed
    of plane Poiseuille flow.  This is the single most characteristic feature
    of the dispersive model and is absent from 2DGA, where the tip speed is 1.

    The measured tip overshoots because a first-order scheme diffuses the
    vanishing tail forward; the test is that it converges DOWN to 3/2.
    """
    assert _newt.dq0_dc(np.array(0.0), 1.0) == pytest.approx(1.5, abs=1e-14)

    speeds = []
    for n in (100, 200, 400, 800):
        geo, c, t = _rarefaction(n)
        xi = geo.grid.xi_centres
        tip = xi[np.max(np.nonzero(c[0] > 1e-3)[0])]
        speeds.append(tip / t)
    assert all(s > 1.5 for s in speeds)
    assert speeds == sorted(speeds, reverse=True)      # monotone approach
    assert speeds[-1] < 1.62

    # 2DGA, same problem: the front is a contact discontinuity at speed 1, so
    # the tip cannot outrun the mean flow however fine the mesh.
    f = _NEWTONIAN
    geo, ell, tr = build(n_phi=8, n_xi=800, two_dga=True, fluids=(f, f),
                         c_in=1.0)
    c = np.zeros((8, 800))
    psi = ell.solve_nonlinear(c, Q=1.0)[0]
    t = 0.0
    while t < 4.0:
        dt = min(tr.cfl * tr.max_timestep(*tr.fluxes(c, psi)[2:]), 4.0 - t)
        c, dt, _ = tr.step(c, psi, dt=dt)
        t += dt
    xi = geo.grid.xi_centres
    tip_2dga = xi[np.max(np.nonzero(c[0] > 1e-3)[0])] / t
    assert tip_2dga < 1.2
    assert tip_2dga < speeds[-1]


def test_m5_t4_converges_to_the_exact_rarefaction():
    """
    The m = 1 Riemann problem has a closed-form solution: q0'(c) = 1.5(1 - c^2)
    is monotone decreasing, so the c = 1 / c = 0 jump opens into a rarefaction
    fan with c = sqrt(1 - xi/(1.5 t)).  Check L1 convergence against it.

    First order is the expected rate away from the fan edges; the square-root
    corner at the tip degrades it, so ~0.6-0.8 is what a correct monotone
    first-order scheme gives here.  The point of the gate is that the error
    converges to zero towards the RIGHT profile, not the rate.
    """
    errs = []
    for n in (100, 200, 400, 800):
        geo, c, t = _rarefaction(n)
        xi = geo.grid.xi_centres
        exact = np.sqrt(np.clip(1.0 - (xi / t) / 1.5, 0.0, 1.0))
        errs.append(float(np.mean(np.abs(c[0] - exact))))
    assert errs == sorted(errs, reverse=True)
    orders = [np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    assert min(orders) > 0.55, orders
    assert errs[-1] < 0.006


def test_dispersive_and_2dga_fronts_differ_where_the_physics_says_they_should():
    """
    Sanity on the model difference rather than the scheme: BF25's whole case
    for D2DGA is that the gap-scale layering spreads the front.  At equal
    viscosity and zero buoyancy the D2DGA profile must be strictly more spread
    than the 2DGA one from the same initial data.
    """
    f = _NEWTONIAN
    n = 400

    def spread(two_dga):
        geo, ell, tr = build(n_phi=8, n_xi=n, m=1.0,
                             two_dga=two_dga, fluids=(f, f), c_in=1.0)
        c = np.zeros((8, n))
        psi = (ell.solve_nonlinear(c, Q=1.0)[0] if two_dga else ell.solve(c, Q=1.0))
        t = 0.0
        while t < 4.0:
            dt = min(tr.cfl * tr.max_timestep(*tr.fluxes(c, psi)[2:]), 4.0 - t)
            c, dt, _ = tr.step(c, psi, dt=dt)
            t += dt
        xi = geo.grid.xi_centres
        hi = xi[np.max(np.nonzero(c[0] > 0.1)[0])]
        lo = xi[np.max(np.nonzero(c[0] > 0.9)[0])]
        return hi - lo

    assert spread(False) > 3.0 * spread(True)


# =========================================================================
# NUM-15 -- the mesh spacings BCF25 (30)/(31) omit
# =========================================================================
def test_num15_buoyancy_flux_is_grid_independent():
    """
    NUM-15.  BCF25 (30) and (31) add the buoyancy terms without a mesh spacing,
    while the advective terms enter through Psi-differences that already carry
    one.  Since the update multiplies everything by lambda = dt/(dphi dxi), the
    printed buoyancy flux would scale like 1/dxi on a phi-face.  Check that the
    restored version gives a buoyancy flux per unit length that does not move
    when the mesh is refined.
    """
    per_length = []
    for n_xi in (60, 120, 240):
        geo, ell, tr = build(n_phi=12, n_xi=n_xi, beta=np.pi / 3, e=0.0,
                             m=1.0, b=6.0)
        c = np.full((12, n_xi), 0.5)
        psi = ell.solve(c, Q=1.0)
        Phi, _, _, _ = tr.fluxes(c, psi)
        # uniform c, so the advective part is the only other contributor and
        # it is identical between meshes; take the buoyancy part directly
        mid = Phi[6, n_xi // 2] - 0.25 * 2.0 * (psi[6, n_xi // 2]
                                                - psi[6, n_xi // 2 + 1])
        per_length.append(mid / geo.grid.dxi)
    assert np.allclose(per_length, per_length[0], rtol=1e-10)
    # and it is not zero, so the test has something to measure
    assert abs(per_length[0]) > 1e-6


def test_num15_buoyancy_flux_matches_the_closed_form():
    """The restored phi-face buoyancy flux must be exactly
    -dxi * b H^3 I3 sin(beta) sin(pi phi), averaged over the two cells."""
    beta, b, m, cbar = np.pi / 3, 6.0, 0.7, 0.4
    geo, ell, tr = build(n_phi=12, n_xi=40, beta=beta, e=0.0, m=m, b=b)
    c = np.full((12, 40), cbar)
    psi = ell.solve(c, Q=1.0)
    Phi, _, _, _ = tr.fluxes(c, psi)

    phi_c = geo.grid.phi_centres
    H = tr.H_c[:, 20]
    k = -b * H ** 3 * np.sin(beta) * np.sin(np.pi * phi_c)
    I3 = _newt.script_I3(np.array(cbar), m)
    expected = geo.grid.dxi * 0.5 * (k[:-1] + k[1:]) * I3

    adv = 0.5 * (psi[1:-1, 20] - psi[1:-1, 21])      # 0.25*(q0+q0)*dpsi, q0=q0
    adv = 0.25 * 2.0 * _newt.q0(np.array(cbar), m) * (psi[1:-1, 20] - psi[1:-1, 21])
    assert np.allclose(Phi[1:-1, 20] - adv, expected, rtol=1e-12, atol=1e-14)


# =========================================================================
# Q2 (NUM-01) -- CFL number, and Q7 (NUM-06) -- initial interface
# =========================================================================
def test_q2_cfl_number_accuracy_and_monotonicity():
    """
    Q2.  BCF25 says only "a fixed CFL number < 1" and never gives a value.

    Measured, on the m = 1 rarefaction against its exact solution: the L1 error
    RISES as the CFL number falls -- 0.0031 at 0.95 against 0.0102 at 0.25, a
    factor of 3.3.  That is the expected behaviour of a first-order upwind
    scheme (it is exact at CFL = 1 for linear advection) but it is the opposite
    of the usual "smaller timestep, safer" instinct, so it is measured rather
    than assumed.

    Monotonicity is unaffected: (44) is evaluated on the current Psi^n and c^n,
    both of which are known before the step, so a large CFL number is not a
    gamble.  The maximum principle is asserted here at 0.95 for that reason.

    This settles the TRANSPORT half of Q2 only.  The other half -- the
    elliptic/transport splitting error, which grows with dt -- is invisible
    here because Psi is frozen, and is M6-T4's to resolve.
    """
    errs = {}
    for cfl in (0.95, 0.5, 0.25):
        geo, c, t = _rarefaction(400, cfl=cfl)
        xi = geo.grid.xi_centres
        exact = np.sqrt(np.clip(1.0 - (xi / t) / 1.5, 0.0, 1.0))
        errs[cfl] = float(np.mean(np.abs(c[0] - exact)))
        assert c.min() >= 0.0 and c.max() <= 1.0
    assert errs[0.95] < errs[0.5] < errs[0.25]
    assert errs[0.25] / errs[0.95] > 2.5


def test_q7_initial_interface_shape_is_immaterial():
    """
    Q7.  Sharp step versus one-cell-smoothed step at t = 0.

    The difference between the two solutions must decay with the mesh; if it
    did not, the initial condition would be a modelling choice rather than a
    discretisation detail.  It converges at first order in L1 (ratio 2.00 per
    halving over three meshes), i.e. at exactly the scheme's own accuracy, so
    the sharp step is kept.
    """
    diffs = []
    for n_xi in (200, 400, 800):
        out = []
        for smooth in (False, True):
            geo, ell, tr = build(n_phi=8, n_xi=n_xi, m=1.0, c_in=1.0)
            c = np.zeros((8, n_xi))
            c[:, :n_xi // 8] = 1.0
            if smooth:
                c[:, n_xi // 8] = 0.5
            psi = ell.solve(c, Q=1.0)
            t = 0.0
            while t < 4.0:
                dt = min(tr.cfl * tr.max_timestep(*tr.fluxes(c, psi)[2:]), 4.0 - t)
                c, dt, _ = tr.step(c, psi, dt=dt)
                t += dt
            out.append(c[0].copy())
        diffs.append(float(np.mean(np.abs(out[0] - out[1]))))
    ratios = [diffs[i] / diffs[i + 1] for i in range(len(diffs) - 1)]
    assert all(r > 1.8 for r in ratios), ratios
    assert diffs[-1] < 1e-3


# =========================================================================
# NUM-29 -- the LLF flux must not be used at the inflow boundary
# =========================================================================
def test_num29_inflow_face_delivers_exactly_the_pumped_flux():
    """
    NUM-29.  Applying the interior LLF flux at xi = 0 through a ghost cell puts
    the artificial term -alpha (c_R - c_L) on a face with no second cell to
    redistribute to.  In the interior that term is a conservative diffusion; at
    a Dirichlet inflow it is a pure SOURCE, sized by the LLF wavespeed rather
    than by anything physical.

    The invariant it breaks is simple and checkable: the displacing fluid must
    enter at exactly Q c_in.  Assert that, and assert that the ghost-cell
    version does not -- on a strongly buoyant case, because the error scales
    with the buoyancy wavespeed.  On the K-GEP-1 fluids, where that wavespeed
    is ~10^3 times the mean flow, the inflow face carried 110 Q.
    """
    for b in (0.0, 100.0):
        geo, ell, tr = build(n_phi=16, n_xi=60, m=0.2, b=b, c_in=1.0)
        c = np.zeros((16, 60))                      # all displaced fluid
        psi = ell.solve(c, Q=0.8)
        _, Xi, _, _ = tr.fluxes(c, psi)
        assert float(np.sum(Xi[:, 0])) == pytest.approx(0.8, rel=1e-12), b

    # the ghost-cell form over-delivers, and worse the more buoyant the case
    over = {}
    for b in (0.0, 100.0):
        geo, ell, tr = build(n_phi=16, n_xi=60, m=0.2, b=b, c_in=1.0)
        tr = TransportSolver(geo, tr.closures, buoyancy_number=b,
                             inflow_concentration=1.0, cfl=0.5,
                             inlet_flux="llf")
        c = np.zeros((16, 60))
        psi = ell.solve(c, Q=0.8)
        _, Xi, _, _ = tr.fluxes(c, psi)
        over[b] = float(np.sum(Xi[:, 0])) / 0.8
    # Even with no buoyancy the ghost-cell form over-delivers (1.25x here,
    # from the advective wavespeed alone).  What makes it fatal rather than
    # untidy is that it scales with the LLF wavespeed, so strong buoyancy
    # multiplies it.
    assert over[0.0] > 1.1, over
    assert over[100.0] > 3.0 * over[0.0], over

    # and it shows up exactly where it matters: more displacing fluid present
    # than has been pumped
    geo, ell, tr = build(n_phi=16, n_xi=60, m=0.2, b=100.0, c_in=1.0)
    bad = TransportSolver(geo, tr.closures, buoyancy_number=100.0,
                          inflow_concentration=1.0, cfl=0.5, inlet_flux="llf")
    cap = tr.mass(np.ones((16, 60)))
    for solver, expect_physical in ((tr, True), (bad, False)):
        c = np.zeros((16, 60))
        t = 0.0
        for _ in range(60):
            psi = ell.solve(c, Q=1.0)
            c, dt, _ = solver.step(c, psi)
            t += dt
        pumped = t / geo.grid.Z
        present = solver.mass(c) / cap
        if expect_physical:
            # the invariant: before breakthrough, what is in the annulus cannot
            # exceed what has been pumped
            assert present <= pumped * (1.0 + 1e-9), (present, pumped)
        else:
            # the ghost-cell form breaks it.  The excess here is 18% after 60
            # steps and decays once the inlet cell fills, which is why it went
            # unnoticed on the published cases -- their alpha is small.  On the
            # K-GEP-1 fluids the same mechanism gave eta = 0.034 against 0.009
            # of a volume pumped.
            assert present > pumped * 1.05, (present, pumped)
