"""
2DGA closure and PF04 (40)-(43).

The equations were unreadable in the supplied PDF (BENCH-02) until the user
provided page images; they are now transcribed and implemented in
`d2dga.gapscale.twodga`.  This file verifies what can be verified.
"""

import numpy as np
import pytest

from d2dga.config import GridConfig, WellConfig
from d2dga.elliptic import StreamFunctionSolver, TwoDGAClosures
from d2dga.gapscale import twodga as T
from d2dga.geometry import ConstantEccentricity, Geometry, UniformWall
from d2dga.scaling import ScaledFluid

# PF04 Figs 3-5 parameter set, from their captions:
# e=0, St*=0.1, rho1=1, rho2=0.9, tau1Y=0.9, tau2Y=0.7,
# kappa1=0.5, kappa2=0.4, m1=1.0, m2=1.2   (m is the INVERSE power index)
PF04_K = [(0.5, 1.0, 0.9), (0.4, 1.2, 0.7)]
PF04_B = (0.9 - 1.0) / 0.1          # = -1.0, density UNSTABLE
SHORT = dict(casing_shoe_m=388.6, total_depth_m=390.0)      # -> Z ~ 4


def pf04_fluids():
    out = []
    for kappa, m, ty in PF04_K:
        n = 1.0 / m
        out.append((float(T.chi(1.0, 1.0, kappa, n, ty)),
                    float(T.dchi_dgrad(1.0, 1.0, kappa, n, ty)), ty, m))
    return out


def pf04_scaled_fluids():
    """
    PF04's fluids, mapped into BF25 ordering.  CONV-09: PF04 Section 2.2 numbers
    fluid 1 as the DISPLACING (lower) one, BF25 Section 2.3 numbers fluid 1 as
    the DISPLACED one.  So PF04's (kappa1, m1, tau1Y) = (0.5, 1.0, 0.9) is the
    DISPLACING fluid here, and (0.4, 1.2, 0.7) the displaced one.

    Returned as (displaced, displacing), which is what StreamFunctionSolver
    expects.
    """
    return (ScaledFluid("displaced", 0.9, 0.4, 1.0 / 1.2, 0.7),
            ScaledFluid("displacing", 1.0, 0.5, 1.0, 0.9))


# BF25 delta_rho = (rho_displaced - rho_displacing)/rho_displaced, and
# b_num = -delta_rho/Fr^2.  The displacing fluid is the denser one here, so the
# configuration is density STABLE and b_num > 0 -- even though PF04 writes the
# same case as b = -1 in its own sign convention.
PF04_DRHO = (0.9 - 1.0) / 0.9
PF04_FR = np.sqrt(abs(PF04_DRHO))


# ------------------------------------------------- the 2DGA closure itself ---
def test_pf04_8_is_algebraically_b02_53():
    """PF04 (8) and B02 (53) are the same law written two ways: substituting
    G = chi + tau_Y/H into (53) reproduces (8).  Since the AL solver was already
    validated against (53) at M2-T2, the 2DGA closure arrives pre-verified."""
    for kappa, m, tau_y in PF04_K:
        for chi_val in (0.3, 1.0, 3.0):
            H = 1.0
            pf04 = T.areal_flux(chi_val, H, kappa, m, tau_y)
            G = chi_val + tau_y / H
            b02 = ((H * G - tau_y) * ((m + 1) * H * G + tau_y)
                   / (G ** 2 * (m + 1) * (m + 2)) * ((H * G - tau_y) / kappa) ** m)
            assert np.isclose(pf04, b02, rtol=1e-12)


def test_chi_inversion_round_trips():
    for kappa, m, tau_y in PF04_K:
        n = 1.0 / m
        for grad in (0.1, 1.0, 5.0):
            c = float(T.chi(grad, 1.0, kappa, n, tau_y))
            assert np.isclose(T.areal_flux(c, 1.0, kappa, m, tau_y), grad, rtol=1e-10)


def test_conv01_bridge_halves_the_bf25_gradient():
    """B02/PF04 omit the factor 2 in the stream function (CONV-01), so every
    |grad Psi| handed to chi must be halved.  Centralised so it cannot be
    forgotten at a call site."""
    kappa, m, tau_y = PF04_K[0]
    n = 1.0 / m
    assert np.isclose(float(T.chi_from_bf25_gradient(2.0, 1.0, kappa, n, tau_y)),
                      float(T.chi(1.0, 1.0, kappa, n, tau_y)), rtol=1e-12)


# --------------------------------------------------- PF04 (40)-(43) itself ---
def test_pf04_40_reduces_to_38_at_zero_eccentricity():
    fl = pf04_fluids()
    phi = np.linspace(0, 1, 11)
    got = T.steady_interface_eccentric(phi, 0.0, PF04_B, np.pi / 2, 2.0, fl)
    J = (fl[1][0] + fl[1][2]) - (fl[0][0] + fl[0][2])
    want = -(1 / np.pi) * PF04_B * np.sin(np.pi / 2) * np.cos(np.pi * phi) / (J + 0.0)
    assert np.allclose(got, want, rtol=1e-12)
    assert np.isclose(abs(got[0]), 0.411, atol=0.002)     # PF04 Fig. 5 shows ~0.42


def test_pf04_42_43_are_continuous_and_vanish_at_the_interface():
    """Both branches contain 1 - cosh(a L)/cosh(a L) at z = 0, so Phi = 0 there.
    That is the consistency condition the pair has to satisfy -- the interface
    is the streamline Phi = 0 -- and it is what makes the two expressions a
    single continuous field rather than two unrelated ones."""
    fl = pf04_fluids()
    phi = np.linspace(0.05, 0.95, 9)
    for eps in (-1e-10, 1e-10):
        assert np.allclose(T.moving_frame_stream_function(phi, eps, 0.05, 2.0, fl),
                           0.0, atol=1e-8)


def test_pf04_42_43_satisfy_the_b02_end_conditions():
    """B02 (69)-(70) impose dPsi/dxi = 0 at both ends.  d/dz cosh(a(L+z)) is
    a sinh(a(L+z)), which vanishes at z = -L; likewise at z = +L.  So PF04's
    perturbation solution is consistent with the same boundary conditions the
    solver uses."""
    fl = pf04_fluids()
    L, e, h = 2.0, 0.05, 1e-6
    for z0, sign in ((-L, +1), (L, -1)):
        d = (T.moving_frame_stream_function(0.5, z0 + sign * h, e, L, fl)
             - T.moving_frame_stream_function(0.5, z0, e, L, fl)) / h
        assert abs(float(d)) < 1e-3


def test_pf04_41_alpha_and_P_are_positive():
    """PF04 states alpha_k^2 > 0 and calls P 'the strictly positive function'.
    Both are structural: a complex alpha or negative P would make (42)-(43)
    meaningless."""
    for (c, d, ty, m) in pf04_fluids():
        assert T.alpha(c, d, ty) > 0
        assert T.P(c, ty, m) > 0


def test_pf04_perturbation_breaks_down_well_before_e_equals_0p2():
    """The O(e) term overtakes the leading term and flips the interface sign by
    e = 0.2, so 'mildly eccentric' is a real restriction, not boilerplate."""
    fl = pf04_fluids()
    g = [float(T.steady_interface_eccentric(np.array([0.0]), e, PF04_B, np.pi/2, 2.0, fl)[0])
         for e in (0.0, 0.05, 0.2)]
    assert g[0] < 0 and g[2] > 0


# ------------------------------------------------- the 2DGA solver itself ---
def test_2dga_solver_reproduces_the_concentric_steady_state():
    """
    M4-T2 for the 2DGA model, with PF04's OWN Herschel-Bulkley fluids and
    beta = pi/2.  Prescribing the (38) interface must give the uniform
    moving-frame flow, Phi = Psi - 2 phi = 0, and the residual must converge.
    """
    f1, f2 = pf04_scaled_fluids()
    fl = pf04_fluids()
    errs = []
    for n_phi, n_xi in ((32, 100), (64, 200)):
        well = WellConfig(inclination_rad=np.pi / 2, **SHORT)
        geo = Geometry(well, UniformWall(well.gauge_hole_diameter_m),
                       ConstantEccentricity(0.0), GridConfig(n_phi, n_xi))
        s = StreamFunctionSolver(geo, TwoDGAClosures(f1, f2),
                                 froude=PF04_FR, delta_rho=PF04_DRHO)
        g = geo.grid
        shape = T.steady_interface_eccentric(g.phi_centres, 0.0, PF04_B,
                                             np.pi / 2, 0.5 * geo.Z, fl)
        c = (g.xi_centres[None, :] < (0.5 * geo.Z + shape)[:, None]).astype(float)
        psi, _, _ = s.solve_nonlinear(c, Q=1.0, tol=1e-10, max_iter=200)
        errs.append(np.abs(psi - 2.0 * g.phi_edges[:, None]).max() / 2)
        assert np.allclose(s.axial_flux(psi), 2.0, atol=1e-9)
    assert errs[-1] < 1e-2
    assert errs[-1] < errs[0]


def test_2dga_buoyancy_has_no_layering_term():
    """The single structural difference from D2DGA: 2DGA has script_I2 = 0, so
    its b vector carries only the mean-density gradient.  That is exactly the
    simplification BF25 Section 2.1 removes."""
    f1, f2 = pf04_scaled_fluids()
    prov = TwoDGAClosures(f1, f2)
    c = np.linspace(0, 1, 7)[:, None] * np.ones((1, 3))
    I1, I2 = prov.script_I1_I2(c, np.ones_like(c), np.ones_like(c))
    assert np.all(I2 == 0.0)
    assert np.all(I1 > 0)
    q0, I3 = prov.q0_I3(c, np.ones_like(c), np.ones_like(c))
    assert np.allclose(q0, c)          # advects with the gap-averaged velocity
    assert np.all(I3 == 0.0)           # and has no buoyancy flux


# ------------------------------------------------------------------ M4-T3 ---
def _phi_error(e, n_phi, n_xi):
    """Solve with the (40) interface prescribed and compare the moving-frame
    stream function against (42)/(43).  Returns (far-field rel, frontal rel)."""
    f1, f2 = pf04_scaled_fluids()
    fl = pf04_fluids()
    well = WellConfig(inclination_rad=0.0, **SHORT)
    geo = Geometry(well, UniformWall(well.gauge_hole_diameter_m),
                   ConstantEccentricity(e), GridConfig(n_phi, n_xi))
    s = StreamFunctionSolver(geo, TwoDGAClosures(f1, f2),
                             froude=PF04_FR, delta_rho=PF04_DRHO)
    g = geo.grid
    L = 0.5 * geo.Z
    shape = T.steady_interface_eccentric(g.phi_centres, e, PF04_B, 0.0, L, fl)
    c = (g.xi_centres[None, :] < (L + shape)[:, None]).astype(float)
    psi, _, _ = s.solve_nonlinear(c, Q=1.0, tol=1e-10, max_iter=250)

    pe = g.phi_edges
    z = g.xi_edges - L
    phi_num = psi - 2.0 * (pe + (e / np.pi) * np.sin(np.pi * pe))[:, None]
    phi_an = 2.0 * T.moving_frame_stream_function(pe[:, None], z[None, :], e, L, fl)
    scale = np.abs(phi_an).max()
    d = np.abs(phi_num - phi_an)
    return (d[:, np.abs(z) > 0.6 * L].max() / scale,
            d[:, np.abs(z) < 0.15 * L].max() / scale)


@pytest.mark.parametrize("e", [0.02, 0.05])
def test_m4_t3_far_field_matches_pf04_to_its_own_truncation(e):
    """
    M4-T3, resolved.

    beta = 0, which is the regime PF04's perturbation is derived for: ZF22
    Section 4.2 notes it holds "when the annulus is mildly eccentric, e << 1,
    AND when the concentric steady-state profile is approximately perpendicular
    to the annulus axis".  At beta = pi/2 that profile is tilted +-0.41, so the
    formula does not apply there at all -- my first attempt tested it outside
    its own regime (BENCH-05).

    (40)-(43) are O(e) asymptotic, so their own truncation error is O(e)
    relative.  The far-field agreement tracks exactly that: 1.10% at e = 0.02
    and 2.70% at e = 0.05, a ratio of 2.45 against an e-ratio of 2.5.  Scaling
    linearly in e is the signature of the FORMULA's truncation, not of an error
    in the solver -- a solver bug would not know about e.
    """
    far, _ = _phi_error(e, 64, 200)
    assert far < 0.9 * e


def test_m4_t3_far_field_error_is_linear_in_e():
    """The discriminating test.  If the residual were a solver error it would
    not scale with the perturbation parameter; scaling linearly identifies it
    as PF04's own O(e^2) truncation."""
    far2, _ = _phi_error(0.02, 64, 200)
    far5, _ = _phi_error(0.05, 64, 200)
    assert 1.8 < far5 / far2 < 3.2          # e-ratio is 2.5


def test_m4_t3_frontal_residual_converges_under_refinement():
    """
    The remaining discrepancy sits at the front, where the sharp interface is
    represented as a staircase on the fixed mesh.  It CONVERGES once the O(e)
    interface tilt is actually resolved -- the same NUM-12 lesson as M4-T2, and
    the reason it matters more here is that the O(e) tilt is small, so it needs
    a fine axial mesh before it spans enough cells to mean anything.
    """
    _, coarse = _phi_error(0.05, 64, 200)
    _, fine = _phi_error(0.05, 96, 300)
    assert fine < coarse


def test_conv09_fluid_numbering_is_the_one_that_bites():
    """
    Regression guard for CONV-09.

    PF04 Section 2.2 numbers fluid 1 as the DISPLACING (lower) one; BF25
    Section 2.3 numbers fluid 1 as the DISPLACED one.  Swapping them changes
    which P and alpha sit on which side of z = 0, and the far-field amplitudes
    of the two fluids differ by ~15% here -- so getting it wrong produces a
    mismatch that does NOT converge under refinement, which is exactly how it
    was found.
    """
    fl = pf04_fluids()
    amp = [T.P(c, ty, m) for (c, _, ty, m) in fl]
    assert abs(amp[0] - amp[1]) / max(amp) > 0.1     # the two sides differ
    displaced, displacing = pf04_scaled_fluids()
    # PF04's fluid 1 (kappa=0.5, tauY=0.9) must be our DISPLACING fluid
    assert displacing.consistency == 0.5 and displacing.yield_stress == 0.9
    assert displaced.consistency == 0.4 and displaced.yield_stress == 0.7
    # and the displacing fluid is the denser one, so the case is density stable
    assert displacing.density > displaced.density
    assert -PF04_DRHO / PF04_FR ** 2 > 0             # b_num > 0 in BF25 signs
