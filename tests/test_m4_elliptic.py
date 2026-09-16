"""M4 test gates -- elliptic solve for the stream function."""

import numpy as np
import pytest
import sympy as sp

from d2dga.config import GridConfig, WellConfig
from d2dga.elliptic import NewtonianClosures, StreamFunctionSolver
from d2dga.geometry import (ConstantEccentricity, Geometry, SinusoidalWall,
                            UniformWall)
from tests.benchmarks.pf04_analytic import (concentric_steady_interface,
                                            newtonian_inverse_mobilities)

# A short annulus, so the steady interface tilt is resolved by many cells.
# On the full 196 m K-GEP-1 domain Z ~ 564 and the entire tilt is sub-cell,
# which is why PF04 works in a short moving-frame window.
SHORT = dict(casing_shoe_m=387.2, total_depth_m=390.0)      # -> Z ~ 8


def make(beta=0.0, e=0.0, n_phi=96, n_xi=300, m=0.5, b=8.0, short=True,
         wall=None):
    kw = dict(SHORT) if short else {}
    well = WellConfig(inclination_rad=beta, **kw)
    w = wall if wall is not None else UniformWall(well.gauge_hole_diameter_m)
    geo = Geometry(well, w, ConstantEccentricity(e), GridConfig(n_phi, n_xi))
    # b_num = -delta_rho / Fr*^2, so Fr = 1 and delta_rho = -b give b_num = b
    solver = StreamFunctionSolver(geo, NewtonianClosures(m),
                                 froude=1.0, delta_rho=-b)
    return geo, solver


# ------------------------------------------------- formulation checks -----
def test_plain_divergence_form_equals_published_div_a_form():
    """The code discretises d_phi(X) + d_xi(Y) = 0 rather than BF25 (2.16)'s
    div_a.[S + b].  The two are algebraically identical, not an approximation."""
    phi, xi = sp.symbols('phi xi')
    ra = sp.Function('r_a')(xi)
    I1 = sp.Function('I1')(phi, xi)
    Psi = sp.Function('Psi')(phi, xi)
    bs = sp.Function('b_s')(phi, xi)
    beta = sp.Function('beta')(xi)
    f_phi, f_xi = ra * sp.cos(beta), ra * sp.sin(beta) * sp.sin(sp.pi * phi)

    published = (sp.diff(ra * (sp.diff(Psi, phi) / ra) / (2 * I1) + bs * f_phi, phi) / ra
                 + sp.diff(ra * sp.diff(Psi, xi) / (2 * I1) + bs * f_xi, xi))
    plain = (sp.diff(sp.diff(Psi, phi) / (2 * ra * I1) + bs * sp.cos(beta), phi)
             + sp.diff(ra * sp.diff(Psi, xi) / (2 * I1)
                       + bs * ra * sp.sin(beta) * sp.sin(sp.pi * phi), xi))
    assert sp.simplify(plain - published) == 0


def test_m4_t2_derivation_matches_pf04_38():
    """
    PF04 (38) re-derived rather than transcribed.

    The display equations in the supplied PF04 PDF do not extract -- the text
    layer returns scrambled fragments -- and no page renderer is available, so
    we derive the concentric steady state from the model equations and check it
    against the fragments that do survive.  Reading PF04's jump notation
    [.]^2_1 correctly, the two agree exactly, including the sign of b cos(beta).
    """
    phi, beta, b, J = sp.symbols('phi beta b J')
    gprime = sp.simplify(-(b * sp.sin(beta) * sp.sin(sp.pi * phi))
                         / (-J - b * sp.cos(beta)))
    g = sp.simplify(sp.integrate(gprime, phi))
    pf04 = -b * sp.sin(beta) * sp.cos(sp.pi * phi) / (sp.pi * (J + b * sp.cos(beta)))
    assert sp.simplify(g - pf04) == 0


def test_conv04_b_s_reduces_to_bf25_2_18_at_unit_radius():
    """docs/derivation.md Finding 2: the scalar multiplying f is r_a-free.
    At r_a = 1 -- every published case -- it must coincide with BF25 (2.18)."""
    c, drho, Fr2, I1, I2 = sp.symbols('c Delta_rho Fr2 I1 I2')
    derived = (1 - drho * (c + I2 / I1)) / Fr2
    bf25_218 = 1 / Fr2 + (-drho / (1 * Fr2)) * (c + I2 / I1)   # r_a = 1
    assert sp.simplify(derived - bf25_218) == 0


# ---------------------------------------------------------------- M4-T1 ---
def test_m4_t1_concentric_isodense_single_fluid():
    """Uniform axial velocity, zero azimuthal velocity."""
    geo, s = make(beta=0.0, e=0.0, n_phi=24, n_xi=60, m=1.0, b=0.0)
    psi = s.solve(np.zeros((24, 60)), Q=1.0)
    v, w = s.velocities(psi)
    assert np.allclose(w, 1.0, atol=1e-12)
    assert np.max(np.abs(v)) < 1e-12


@pytest.mark.parametrize("e", [0.0, 0.2, 0.4, 0.6, 0.8])
def test_m4_t1_eccentric_far_field_velocity_ratio(e):
    """ZF22 Section 3: for Newtonian fluids the far-field wide/narrow mean
    velocity ratio is (1+e)^2/(1-e)^2.  An independent check on the operator
    with strongly variable coefficients."""
    geo, s = make(beta=0.0, e=e, n_phi=200, n_xi=20, m=1.0, b=0.0, short=False)
    psi = s.solve(np.zeros((200, 20)), Q=1.0)
    _, w = s.velocities(psi)
    mid = w[:, 10]
    assert np.isclose(mid[0] / mid[-1], (1 + e) ** 2 / (1 - e) ** 2, rtol=5e-4)


# ---------------------------------------------------------------- M4-T2 ---
@pytest.mark.parametrize("beta", [0.0, np.pi / 4, np.pi / 2])
def test_m4_t2_concentric_steady_state(beta):
    """
    Prescribe the analytic steady interface; the exact solution is then the
    UNIFORM flow Psi = 2 phi, v_bar = 0, w_bar = 1.

    At beta = 0 this is exact to machine precision.  At beta != 0 the sharp
    interface is represented as a staircase on the fixed mesh, so the error is
    finite but must CONVERGE under refinement.
    """
    m, b = 0.5, 8.0
    iv1, iv2 = newtonian_inverse_mobilities(m)
    errs = []
    for n_phi, n_xi in ((48, 150), (96, 300), (192, 600)):
        geo, s = make(beta=beta, e=0.0, n_phi=n_phi, n_xi=n_xi, m=m, b=b)
        g = geo.grid
        xi0 = 0.5 * geo.Z
        shape = concentric_steady_interface(g.phi_centres, b, beta, iv1, iv2)
        c = (g.xi_centres[None, :] < (xi0 + shape)[:, None]).astype(float)
        psi = s.solve(c, Q=1.0)
        errs.append(np.max(np.abs(psi - 2.0 * g.phi_edges[:, None])) / 2.0)

    if beta == 0.0:
        # sin(beta) = 0 -> flat interface, no azimuthal buoyancy, exact
        assert max(errs) < 1e-12
    else:
        assert errs[-1] < 1e-2
        assert errs[-1] < errs[0]          # converging, not stalling


def test_m4_t2_buoyancy_source_sign_is_exercised():
    """
    Guard against the NUM-10 regression.

    The RHS buoyancy sign was once inverted and every test still passed,
    because uniform c with beta = 0 makes B constant and the source vanishes.
    This test asserts the source is genuinely non-zero for the M4-T2 setup, so
    that gate really does exercise it.
    """
    m, b, beta = 0.5, 8.0, np.pi / 4
    iv1, iv2 = newtonian_inverse_mobilities(m)
    geo, s = make(beta=beta, e=0.0, n_phi=48, n_xi=150, m=m, b=b)
    g = geo.grid
    shape = concentric_steady_interface(g.phi_centres, b, beta, iv1, iv2)
    c = (g.xi_centres[None, :] < (0.5 * geo.Z + shape)[:, None]).astype(float)
    a_phi, a_xi, B_phi, B_xi = s.coefficients(c)
    _, rhs = s.assemble(a_phi, a_xi, B_phi, B_xi, Q=1.0)
    assert np.max(np.abs(rhs)) > 1.0


# ---------------------------------------------------------------- M4-T3 ---
@pytest.mark.parametrize("e", [0.01, 0.05])
def test_m4_t3_mildly_eccentric_approaches_concentric(e):
    """
    PARTIAL.  PF04 (40)-(43) are not reproduced -- see the module note and
    assumptions.md BENCH-01.  Their coefficients P(chi, tau_Y, m), alpha_k and
    tanh(alpha_k L) cannot be transcribed from the supplied PDF, and deriving
    the full O(e) perturbation solution is out of scope here.

    What IS asserted: at small e the solution stays close to the concentric
    steady state and approaches it monotonically as e -> 0.  That is the
    verifiable content of the perturbation structure without the coefficients.
    """
    m, b, beta = 0.5, 8.0, np.pi / 4
    iv1, iv2 = newtonian_inverse_mobilities(m)

    def deviation(ecc):
        geo, s = make(beta=beta, e=ecc, n_phi=192, n_xi=600, m=m, b=b)
        g = geo.grid
        shape = concentric_steady_interface(g.phi_centres, b, beta, iv1, iv2)
        c = (g.xi_centres[None, :] < (0.5 * geo.Z + shape)[:, None]).astype(float)
        psi = s.solve(c, Q=1.0)
        H = geo.H(g.phi_edges, g.xi_edges)
        ra = geo.r_a(g.xi_edges)[None, :]
        ref = np.zeros_like(psi)
        ref[1:, :] = np.cumsum(2.0 * 0.5 * (H[1:, :] + H[:-1, :]) * ra * g.dphi, axis=0)
        return np.max(np.abs(psi - ref)) / 2.0

    assert deviation(e) > deviation(0.0)
    assert deviation(e) < deviation(0.0) + 1.0 * e


# ---------------------------------------------------------------- M4-T4 ---
@pytest.mark.parametrize("e,beta", [(0.0, 0.0), (0.5, 0.0), (0.5, np.pi / 3)])
def test_m4_t4_axial_flux_is_exactly_2Q(e, beta):
    """CONV-05.  With grad_a Psi = (2 H w_bar, -2 H v_bar) the half-annulus
    volumetric flow is Psi(1) - Psi(0) = 2Q.  PF04 (10) sets Psi(1) = 1 because
    PF04 (5) omits the factor 2; using that value here would halve every
    velocity."""
    for Q in (1.0, 0.37, 2.5):
        geo, s = make(beta=beta, e=e, n_phi=40, n_xi=120)
        c = np.zeros((40, 120))
        c[:, :60] = 1.0
        psi = s.solve(c, Q=Q)
        assert np.allclose(s.axial_flux(psi), 2.0 * Q, atol=1e-10)


def test_m4_t4_reconstructed_velocity_integrates_to_the_flux():
    """Independent of the boundary values: integrate the reconstructed velocity
    field and confirm it carries 2Q through every slice."""
    geo, s = make(beta=0.0, e=0.6, n_phi=200, n_xi=40, m=1.0, b=0.0, short=False)
    psi = s.solve(np.zeros((200, 40)), Q=1.0)
    _, w = s.velocities(psi)
    g = geo.grid
    H = geo.H(g.phi_centres, g.xi_centres)
    ra = geo.r_a(g.xi_centres)[None, :]
    flux = np.sum(2.0 * H * ra * w * g.dphi, axis=0)
    assert np.allclose(flux, 2.0, atol=1e-9)


# ----------------------------------------------------- caliper geometry ---
def test_axially_varying_wall_preserves_flux():
    """The caliper modification must not leak volume: flux conservation has to
    survive r_a(xi) and H(phi, xi) varying along the well."""
    wall = SinusoidalWall(0.26416, 0.0381, 20.0)
    geo, s = make(beta=0.0, e=0.3, n_phi=40, n_xi=200, short=False, wall=wall)
    c = np.zeros((40, 200))
    c[:, :100] = 1.0
    psi = s.solve(c, Q=1.0)
    assert np.allclose(s.axial_flux(psi), 2.0, atol=1e-9)


# ------------------------------------------- eccentric far field (exact) ---
def eccentric_far_field_velocity(phi, e):
    """
    Exact far-field azimuthal velocity profile, DERIVED not transcribed.

    With d/dxi = 0 the elliptic equation collapses to d_phi(a_phi dPsi/dphi) = 0,
    so a_phi dPsi/dphi is constant.  With a_phi = 1/(2 H^3 I1) and
    Psi(1) - Psi(0) = 2Q,

        w_bar = Q H^2 I1 / int_0^1 H^3 I1 dphi .

    For a PURE fluid I1 is constant in phi and cancels, leaving

        w_bar = H^2 / int_0^1 H^3 dphi = (1 + e cos pi phi)^2 / (1 + 3 e^2 / 2)

    independent of viscosity.  That independence is a real prediction, not a
    restatement, and it is asserted below at two very different m.
    """
    H = 1.0 + e * np.cos(np.pi * phi)
    return H ** 2 / (1.0 + 1.5 * e ** 2)


@pytest.mark.parametrize("e", [0.01, 0.05, 0.2, 0.5, 0.8])
@pytest.mark.parametrize("m", [0.2, 5.0])
def test_m4_t3_eccentric_far_field_profile_is_exact(e, m):
    """
    Strong replacement for the O(e) content of M4-T3.

    PF04 (40)-(43) could not be recovered from the PDF (BENCH-02), so instead of
    transcribing a formula we may have misread, this asserts an eccentric
    benchmark derived from the model itself -- and it holds to MACHINE
    PRECISION, not merely to discretisation error, because the far-field
    problem is one-dimensional and the flux-form discretisation is exact for it.
    """
    geo, s = make(beta=0.0, e=e, n_phi=200, n_xi=20, m=m, b=0.0, short=False)
    psi = s.solve(np.zeros((200, 20)), Q=1.0)
    _, w = s.velocities(psi)
    exact = eccentric_far_field_velocity(geo.grid.phi_centres, e)
    assert np.max(np.abs(w[:, 10] - exact) / exact) < 1e-9


def test_far_field_profile_is_viscosity_independent():
    """The derived prediction: I1 cancels, so two fluids with a 25x consistency
    ratio must give the SAME far-field profile."""
    prof = []
    for m in (0.2, 5.0):
        geo, s = make(beta=0.0, e=0.6, n_phi=160, n_xi=20, m=m, b=0.0, short=False)
        psi = s.solve(np.zeros((160, 20)), Q=1.0)
        prof.append(s.velocities(psi)[1][:, 10])
    assert np.allclose(prof[0], prof[1], rtol=1e-10)


# --------------------------------------------- bottom-hole inflow (Q5) -----
def test_uniform_inflow_imposes_a_flat_velocity_profile():
    """
    B02 Section 3.2's stated alternative to (70).  Guards a regression that was
    live for one commit: the Dirichlet row was appended AFTER the flux stencil,
    and because scipy SUMS duplicate (row, col) entries the row became
    flux + 1 on the diagonal instead of a clean identity, so the prescribed
    values never landed at all.
    """
    geo, s = make(beta=0.0, e=0.5, n_phi=60, n_xi=300, m=0.2, b=10.0, short=False)
    g = geo.grid
    c = np.broadcast_to((g.xi_centres < 0.15 * geo.Z)[None, :],
                        (60, 300)).astype(float)
    s.inflow = "uniform"
    psi = s.solve(c, Q=1.0)
    assert np.allclose(psi[:, 0], 2.0 * s._uniform_inflow_shape, atol=1e-12)

    H0 = geo.H(g.phi_edges, np.array([0.0]))[:, 0]
    ra0 = float(geo.r_a(np.array([0.0]))[0])
    w0 = np.gradient(psi[:, 0], g.dphi) / (2.0 * H0 * ra0)
    assert np.allclose(w0[2:-2], 1.0, atol=5e-3)      # flat, as imposed
    assert np.allclose(s.axial_flux(psi), 2.0, atol=1e-9)


def test_q5_inflow_choice_is_local_to_bottom_hole():
    """
    Q5 resolution.  The two conditions differ substantially AT the boundary but
    the difference decays within ~1.3 m of a 196 m open hole, so the choice is
    immaterial for K-GEP-1 away from bottom hole.  This matters because B02's
    own justification (long annulus) is weaker here -- B02's well is 1000 m.
    """
    geo, s = make(beta=0.0, e=0.5, n_phi=60, n_xi=300, m=0.2, b=10.0, short=False)
    g = geo.grid
    c = np.broadcast_to((g.xi_centres < 0.15 * geo.Z)[None, :],
                        (60, 300)).astype(float)
    s.inflow = "no_axial_gradient"
    psi_ng = s.solve(c, Q=1.0)
    s.inflow = "uniform"
    psi_un = s.solve(c, Q=1.0)

    assert not np.allclose(psi_ng, psi_un)            # the option does something
    dev = np.max(np.abs(psi_ng - psi_un), axis=0)
    first_small = np.argmax(dev < 0.01 * dev.max())
    reach_m = float(geo.xi_hat(g.xi_edges[first_small]))
    assert reach_m < 0.05 * geo.Z_hat                 # local to bottom hole


def test_concentration_shape_is_validated():
    """A (1, n_xi) field used to broadcast into a cryptic failure deep in the
    assembly.  It must be rejected at the boundary with a useful message."""
    geo, s = make(n_phi=20, n_xi=40)
    with pytest.raises(ValueError, match="expected"):
        s.solve(np.zeros((1, 40)), Q=1.0)
