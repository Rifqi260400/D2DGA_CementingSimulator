"""
M4 -- elliptic solve for the stream function.

The published form is BF25 (2.16), div_a . [S + b] = 0, with
div_a = ((1/r_a) d/dphi, d/dxi).  That is awkward to discretise because it is
not a plain divergence.  The derivation in docs/derivation.md produces, en
route, an equivalent form that IS:

    d/dphi [ a_phi dPsi/dphi + B_phi ] + d/dxi [ a_xi dPsi/dxi + B_xi ] = 0

    a_phi = 1 / (2 r_a H^3 I1)          B_phi = b_s cos(beta)
    a_xi  = r_a / (2 H^3 I1)            B_xi  = b_s r_a sin(beta) sin(pi phi)

    b_s   = (1/Fr*^2) [ 1 - Delta_rho ( c + I2/I1 ) ]

where I1, I2 are the H-free script-I closures.  The equivalence is exact, not
asymptotic, and is asserted symbolically in tests/test_m4_elliptic.py.

Note b_s carries NO r_a -- see docs/derivation.md Finding 2 (CONV-04).  BF25
(2.12) and (2.18) use the symbol `b` for two quantities differing by r_a; the
distinction is invisible at r_a = 1, which is every published case, but not
under a caliper wall where r_a varies by 65% through the washout.

Since a_phi, a_xi > 0 the operator is symmetric negative-definite, so for
Newtonian fluids -- where I1 does not depend on the solution -- this is a single
sparse direct solve.  ZF23 Section 2.2.2 states this explicitly: the D2DGA
stream-function equation is LINEAR for Newtonian fluids.  No augmented
Lagrangian is needed or used there.

Mesh (BCF25 Fig. 2): Psi at cell CORNERS, concentration at cell CENTRES.
Coefficients therefore live on the two staggered face families, and the
concentration is averaged onto them.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .gapscale import newtonian as _newt
from .gapscale import twodga as _twodga
from .geometry import Geometry


def _extrema_of(*funcs, n_scan: int = 4001):
    """
    Interior local maxima of |f| for each f, refined to machine precision.

    Used to make the LLF interval maximum EXACT for closed-form closures rather
    than a sampled bound.  A coarse scan locates each peak's bracketing triple;
    a bounded scalar minimisation of -f then refines it.
    """
    from scipy.optimize import minimize_scalar

    grid = np.linspace(0.0, 1.0, n_scan)
    out = []
    for f in funcs:
        v = np.asarray(f(grid), dtype=float)
        peaks = np.nonzero((v[1:-1] > v[:-2]) & (v[1:-1] > v[2:]))[0] + 1
        for k in peaks:
            r = minimize_scalar(lambda c: -float(f(np.array(c))),
                                bounds=(grid[k - 1], grid[k + 1]),
                                method="bounded",
                                options={"xatol": 1e-13})
            out.append(float(r.x))
    return np.array(sorted(out))


# ===========================================================================
# Closure providers
# ===========================================================================
class ClosureProvider(ABC):
    """Supplies the H-free script-I closures at a given state."""

    @abstractmethod
    def script_I1_I2(self, c, H, umag):
        ...

    @abstractmethod
    def q0_I3(self, c, H, umag):
        ...

    def dq0_dI3_dc(self, c, H, umag, h: float = 1e-6):
        """
        d q0/d c and d I3/d c, needed for the LLF wavespeeds BCF25 (32)-(35).

        The default is a central difference clipped to [0, 1], with a small
        safety factor.  The factor matters: LLF monotonicity requires the
        wavespeed to bound the true slope from ABOVE at the neighbouring cell
        values, and a difference quotient can sit O(h^2) BELOW it.  Providers
        with closed-form closures override this with exact derivatives and no
        safety factor.  assumptions.md NUM-17.
        """
        c = np.asarray(c, dtype=float)
        lo = np.clip(c - h, 0.0, 1.0)
        hi = np.clip(c + h, 0.0, 1.0)
        span = np.where(hi > lo, hi - lo, 1.0)
        q_hi, I_hi = self.q0_I3(hi, H, umag)
        q_lo, I_lo = self.q0_I3(lo, H, umag)
        safety = 1.0 + 1e-6
        return (safety * (q_hi - q_lo) / span,
                safety * (I_hi - I_lo) / span)

    def wavespeed_nodes(self):
        """
        Interior concentrations at which |dq0/dc| or |dI3/dc| can peak.

        The LLF coefficient must be the maximum of the flux-function slope over
        the WHOLE interval between the two states at a face, not just at its
        two endpoints -- see transport.py's NUM-26 note and BCF25 (32)-(35),
        which evaluate the endpoints only.  The endpoints are always probed; the
        transport solver additionally probes every node returned here that lies
        inside the interval, so this list must contain every interior extremum
        of both derivatives for the bound to be exact.

        The default is a uniform sample, which is a bound rather than a
        guarantee; providers with closed-form closures return their actual
        extrema.
        """
        return np.linspace(0.0, 1.0, 65)[1:-1]

    @property
    def is_linear(self) -> bool:
        """True when the closures do not depend on the velocity, so the
        elliptic problem is linear and one direct solve suffices."""
        return False


class NewtonianClosures(ClosureProvider):
    """Analytic BF25 (2.24)-(2.26) plus the corrected I3.  Velocity-independent,
    hence linear."""

    def __init__(self, m: float):
        self.m = float(m)

    def script_I1_I2(self, c, H, umag=None):
        return _newt.script_I1(c, self.m), _newt.script_I2(c, self.m)

    def q0_I3(self, c, H, umag=None):
        return _newt.q0(c, self.m), _newt.script_I3(c, self.m)

    def dq0_dI3_dc(self, c, H, umag=None, h=None):
        return _newt.dq0_dc(c, self.m), _newt.dscript_I3_dc(c, self.m)

    def wavespeed_nodes(self):
        """The exact interior extrema of |q0'| and |I3'| at this m, located
        once and cached.  Both are smooth rationals with at most a couple of
        interior peaks -- |q0'| has none for m <= 1.5 and one above it, |I3'|
        has two -- so the interval maximum this buys is exact, not a sample."""
        if getattr(self, "_nodes", None) is None:
            self._nodes = _extrema_of(
                lambda c: np.abs(_newt.dq0_dc(c, self.m)),
                lambda c: np.abs(_newt.dscript_I3_dc(c, self.m)))
        return self._nodes

    @property
    def is_linear(self):
        return True


class TwoDGAClosures(ClosureProvider):
    """
    The ORIGINAL B02/PF04 closure, for the 2DGA baseline.

    Two uses: reproducing PF04's analytic steady states (M4-T2, M4-T3), and
    giving the thesis the 2DGA-vs-D2DGA comparison that ZF22/ZF23 run
    throughout -- without it there is no baseline to improve on.

    Expressed through the same script_I1 interface as D2DGA via the effective
    mobility in `twodga.effective_script_I1`, so one elliptic assembly serves
    both models.  The ONLY other difference is the buoyancy: 2DGA has no
    layering term, so script_I2 = 0 identically.  That single change is what
    separates the two models' b vectors, which makes the comparison clean.

    Rheological parameters are mixed LINEARLY in concentration, as B02
    Section 2.2.1 does -- B02 notes the mixture laws are "used for simplicity,
    with no particular physical justification".
    """

    def __init__(self, fluid1, fluid2):
        self.f1, self.f2 = fluid1, fluid2

    def _mix(self, c):
        c = np.asarray(c, dtype=float)
        kappa = (1 - c) * self.f1.consistency + c * self.f2.consistency
        n = (1 - c) * self.f1.power_law_index + c * self.f2.power_law_index
        tau_y = (1 - c) * self.f1.yield_stress + c * self.f2.yield_stress
        return kappa, n, tau_y

    def script_I1_I2(self, c, H, umag):
        if umag is None:
            umag = np.ones_like(np.asarray(c, dtype=float))
        kappa, n, tau_y = self._mix(c)
        grad = 2.0 * np.asarray(H, dtype=float) * np.asarray(umag, dtype=float)
        I1 = _twodga.effective_script_I1(grad, H, kappa, n, tau_y)
        return I1, np.zeros_like(I1)

    def q0_I3(self, c, H, umag):
        """2DGA advects with the gap-averaged velocity and has no buoyancy
        flux: q0 = c and I3 = 0.  That is precisely the simplification D2DGA
        exists to remove (BF25 Section 2.1)."""
        c = np.asarray(c, dtype=float)
        return c, np.zeros_like(c)

    def dq0_dI3_dc(self, c, H, umag=None, h=None):
        """q0 = c and I3 = 0 identically, so the wavespeeds are exact: the
        transport reduces to BCF25 (11)-(21), the 2DGA scheme."""
        c = np.asarray(c, dtype=float)
        return np.ones_like(c), np.zeros_like(c)

    def wavespeed_nodes(self):
        """Both derivatives are constant, so there is nothing between the
        endpoints that the endpoints do not already bound."""
        return np.empty(0)


class TabulatedClosures(ClosureProvider):
    """Herschel-Bulkley path: interpolate a pre-built M3 table.  Velocity
    dependent, so the elliptic problem is nonlinear and is driven to a fixed
    point by Picard iteration."""

    def __init__(self, table, gb: float = 0.0, inclination_rad: float = 0.0):
        # The gap-scale problem depends on the buoyancy vector G~_b, which is
        # two-dimensional in general:  G_b = -b (f_xi, -f_phi)  with
        # f = (r_a cos beta, r_a sin beta sin(pi phi)).  For a VERTICAL well
        # f_xi = 0, so G_b is purely axial and a single scalar axis suffices --
        # which is the case tabulated here, and is K-GEP-1.
        #
        # For an inclined well G_b acquires a phi-dependent azimuthal component
        # and the table would need an extra axis plus the angle between u_bar
        # and G~_b.  Rather than silently return the wrong closure, refuse.
        # assumptions.md NUM-14.
        if abs(float(inclination_rad)) > 1e-12:
            raise NotImplementedError(
                "TabulatedClosures assumes a vertical well (beta = 0), where the "
                "gap-scale buoyancy vector is purely axial and one scalar gb axis "
                "is enough. For beta != 0 the table needs an azimuthal buoyancy "
                "axis and an angle axis; see assumptions.md NUM-14.")
        self.table = table
        self.gb = float(gb)

    @staticmethod
    def _umag(umag, c):
        """The first elliptic solve of a Picard sequence has no velocity yet
        and passes None.  NewtonianClosures ignores it and TwoDGAClosures
        substitutes 1; a tabulated table must do the same, or np.asarray(None)
        puts NaN on the interpolation axis and the elliptic matrix comes back
        exactly singular."""
        if umag is None:
            return np.ones_like(np.asarray(c, dtype=float))
        return umag

    def script_I1_I2(self, c, H, umag=None):
        I1, I2, _, _ = self.table(c, H=H, umag=self._umag(umag, c), gb=self.gb)
        return I1, I2

    def q0_I3(self, c, H, umag=None):
        _, _, q0, I3 = self.table(c, H=H, umag=self._umag(umag, c), gb=self.gb)
        return q0, I3

    def dq0_dI3_dc(self, c, H, umag=None, h=None):
        """Upper bounds on the two c-derivatives, straight from the table's own
        piecewise-linear structure rather than by finite-differencing it.  See
        ClosureTable._build_derivative_bounds: a bound is what the LLF
        wavespeeds need, and this costs one interpolation instead of eight."""
        return self.table.derivative_bounds(c, H=H, umag=self._umag(umag, c),
                                            gb=self.gb)

    def wavespeed_nodes(self):
        """The table's own c nodes.  Between them the interpolant is linear, so
        its derivative is piecewise constant and the endpoints of each segment
        bound it -- probing the nodes therefore makes the interval maximum
        exact for the interpolated closure that the solver actually uses."""
        return np.asarray(self.table.c_grid, dtype=float)[1:-1]


# ===========================================================================
# Solver
# ===========================================================================
class StreamFunctionSolver:
    """
    Assembles and solves the 5-point conservative discretisation.

    Boundary conditions, B02 Section 3.2, converted to BF25's factor-of-2
    stream-function convention (CONV-01, CONV-05):

        Psi(0, xi)  = 0                       wide side
        Psi(1, xi)  = 2 Q(t)                  narrow side
        dPsi/dxi    = 0    at xi = 0 and Z    B02 (69), (70)

    The factor 2 on the narrow-side value is NOT cosmetic.  With
    grad_a Psi = (2 H w_bar, -2 H v_bar) we have dPsi/dphi = 2 H r_a w_bar, so
    the half-annulus volumetric flow is int_0^1 dPsi/dphi dphi = Psi(1) - Psi(0).
    A uniform concentric annulus at unit mean velocity has H = r_a = w_bar = 1
    and integrates to 2.  PF04 (10) sets Psi(1) = 1 because PF04 (5) omits the
    factor 2.  Using PF04's value with BF25's convention halves every velocity.

    The bottom-hole condition is selectable.  B02 itself calls its choice
    uncertain, and its justification -- entry effects negligible over a long
    annulus -- is weaker here: B02's well is 1000 m, K-GEP-1's open hole is
    196 m.  See assumptions.md NUM-04.
    """

    def __init__(self, geometry: Geometry, closures: ClosureProvider,
                 froude: float, delta_rho: float,
                 inflow: str = "no_axial_gradient",
                 mobility_floor: float = 1e-12):
        self.geom = geometry
        self.closures = closures
        self.Fr2 = float(froude) ** 2
        self.delta_rho = float(delta_rho)
        if inflow not in ("no_axial_gradient", "uniform"):
            raise ValueError(f"unknown inflow condition {inflow!r}")
        self.inflow = inflow
        self.mobility_floor = float(mobility_floor)
        self.n_static_cells = 0

        g = geometry.grid
        self.n_phi, self.n_xi = g.n_phi, g.n_xi
        self.dphi, self.dxi = g.dphi, g.dxi

        # --- geometry on the two staggered face families -------------------
        # a_phi faces: (phi cell-centre, xi node)   shape (n_phi,   n_xi+1)
        # a_xi  faces: (phi node,        xi centre) shape (n_phi+1, n_xi)
        self.H_pf = geometry.H(g.phi_centres, g.xi_edges)
        self.H_xf = geometry.H(g.phi_edges, g.xi_centres)
        self.ra_pf = geometry.r_a(g.xi_edges)[None, :]
        self.ra_xf = geometry.r_a(g.xi_centres)[None, :]
        self.beta = geometry.well.inclination_rad
        self.sin_pi_phi_xf = np.sin(np.pi * g.phi_edges)[:, None]

        # Uniform-inflow profile, B02 Section 3.2's stated alternative to (70):
        # impose a uniform inflow VELOCITY rather than zero axial gradient.
        # w_bar = const at xi = 0 means dPsi/dphi = 2 H r_a there, so
        #     Psi(phi, 0) = 2Q * int_0^phi H r_a dphi' / int_0^1 H r_a dphi'.
        H0 = geometry.H(g.phi_edges, np.array([0.0]))[:, 0]
        ra0 = float(geometry.r_a(np.array([0.0]))[0])
        seg = 0.5 * (H0[1:] + H0[:-1]) * ra0 * g.dphi
        cum = np.concatenate([[0.0], np.cumsum(seg)])
        self._uniform_inflow_shape = cum / cum[-1]      # 0 at phi=0, 1 at phi=1

    # -- concentration onto the face families -------------------------------
    def _c_on_faces(self, c):
        """c is cell-centred, shape (n_phi, n_xi).  Average onto each face
        family, replicating at the domain edges."""
        c = np.asarray(c, dtype=float)
        expected = (self.n_phi, self.n_xi)
        if c.shape != expected:
            raise ValueError(
                f"concentration field has shape {c.shape}, expected {expected} "
                f"(cell centres: n_phi x n_xi). A common slip is building it as "
                f"xi_centres[None, :] < ..., which gives (1, n_xi) and then "
                f"broadcasts into a cryptic failure deeper in the assembly.")
        pad_x = np.concatenate([c[:, :1], c, c[:, -1:]], axis=1)
        c_pf = 0.5 * (pad_x[:, :-1] + pad_x[:, 1:])          # (n_phi, n_xi+1)
        pad_p = np.concatenate([c[:1, :], c, c[-1:, :]], axis=0)
        c_xf = 0.5 * (pad_p[:-1, :] + pad_p[1:, :])          # (n_phi+1, n_xi)
        return c_pf, c_xf

    # -- coefficients --------------------------------------------------------
    def coefficients(self, c, umag_pf=None, umag_xf=None):
        c_pf, c_xf = self._c_on_faces(c)
        I1_pf, I2_pf = self.closures.script_I1_I2(c_pf, self.H_pf, umag_pf)
        I1_xf, I2_xf = self.closures.script_I1_I2(c_xf, self.H_xf, umag_xf)

        # Static cells.  For a yield-stress fluid I1 -> 0 where the material is
        # unyielded.  In this formulation a = 1/(2 r_a H^3 I1) -> infinity, which
        # is the CORRECT limit -- it drives grad Psi -> 0 so nothing flows -- but
        # it cannot be represented in floating point, and I2/I1 becomes 0/0.
        #
        # We therefore floor the mobility.  This is a REGULARISATION, and PF04
        # Section 5 is explicit that regularisation gives an overly optimistic
        # picture of mud removal because it lets nominally static fluid creep,
        # with the error accumulating over the long timescale of a cementing
        # job.  Capturing true unyielded zones needs the augmented Lagrangian
        # treatment of PF04 Section 3.1 applied at the 2-D level.
        #
        # For Newtonian fluids I1 > 0 everywhere and this never triggers, which
        # is why `n_static_cells` is reported: a non-zero count on a Newtonian
        # run means something is wrong, and on a yield-stress run it marks
        # exactly where the regularisation is load-bearing.
        # assumptions.md NUM-13.
        floor = self.mobility_floor
        self.n_static_cells = int(np.sum(I1_pf < floor) + np.sum(I1_xf < floor))
        I1_pf = np.maximum(I1_pf, floor)
        I1_xf = np.maximum(I1_xf, floor)

        a_phi = 1.0 / (2.0 * self.ra_pf * self.H_pf ** 3 * I1_pf)
        a_xi = self.ra_xf / (2.0 * self.H_xf ** 3 * I1_xf)

        bs_pf = (1.0 - self.delta_rho * (c_pf + I2_pf / I1_pf)) / self.Fr2
        bs_xf = (1.0 - self.delta_rho * (c_xf + I2_xf / I1_xf)) / self.Fr2
        B_phi = bs_pf * np.cos(self.beta)
        B_xi = bs_xf * self.ra_xf * np.sin(self.beta) * self.sin_pi_phi_xf
        return a_phi, a_xi, B_phi, B_xi

    # -- assembly -------------------------------------------------------------
    def _index(self, i, j):
        """Unknown ordering: interior phi nodes i = 1..n_phi-1, all xi nodes."""
        return (i - 1) * (self.n_xi + 1) + j

    def assemble_reference(self, a_phi, a_xi, B_phi, B_xi, Q):
        """
        Row-by-row assembly.  Kept as the REFERENCE: it is the version that was
        written first, reviewed line by line against the discretisation, and
        that every M4 gate was established on.  `assemble` is a vectorised
        rewrite that must reproduce it bit for bit -- asserted in
        tests/test_m6_simulation.py -- because the coupled loop reassembles the
        matrix at every one of ~10^3-10^5 timesteps and a Python loop over
        n_phi x n_xi cells dominates the runtime.
        """
        n_i, n_j = self.n_phi - 1, self.n_xi + 1
        n = n_i * n_j
        dphi2, dxi2 = self.dphi ** 2, self.dxi ** 2
        rows, cols, vals = [], [], []
        rhs = np.zeros(n)
        psi_narrow = 2.0 * Q                # CONV-05, see class docstring

        for i in range(1, self.n_phi):
            for j in range(self.n_xi + 1):
                k = self._index(i, j)

                # --- bottom-hole Dirichlet, uniform-inflow option -----------
                # B02 Section 3.2 offers this as the alternative to (70):
                # impose a uniform inflow VELOCITY instead of zero axial
                # gradient.  B02 calls its own choice uncertain, and its
                # justification -- entry effects negligible over a long annulus
                # -- is weaker here: B02's well is 1000 m, K-GEP-1's open hole
                # is 196 m.  assumptions.md NUM-04.
                #
                # This MUST be the first thing done for the row.  scipy sums
                # duplicate (row, col) entries, so appending an identity after
                # the flux stencil gives flux + 1 on the diagonal rather than a
                # clean Dirichlet row, and the prescribed value never lands.
                if self.inflow == "uniform" and j == 0:
                    rows.append(k)
                    cols.append(k)
                    vals.append(1.0)
                    rhs[k] = 2.0 * Q * self._uniform_inflow_shape[i]
                    continue

                diag = 0.0

                # --- phi fluxes (Dirichlet at i = 0 and i = n_phi) ----------
                ae, aw = a_phi[i, j], a_phi[i - 1, j]
                diag -= (ae + aw) / dphi2
                if i + 1 == self.n_phi:
                    rhs[k] -= ae * psi_narrow / dphi2
                else:
                    rows.append(k); cols.append(self._index(i + 1, j)); vals.append(ae / dphi2)
                if i - 1 != 0:
                    rows.append(k); cols.append(self._index(i - 1, j)); vals.append(aw / dphi2)
                # (Psi = 0 on the wide side contributes nothing to the rhs)

                # --- xi fluxes (homogeneous Neumann via ghost nodes) --------
                an = a_xi[i, j] if j < self.n_xi else 0.0
                asf = a_xi[i, j - 1] if j > 0 else 0.0
                diag -= (an + asf) / dxi2
                if j + 1 <= self.n_xi:
                    rows.append(k); cols.append(self._index(i, j + 1)); vals.append(an / dxi2)
                if j - 1 >= 0:
                    if self.inflow == "uniform" and j == 1:
                        # the j = 0 row is Dirichlet, so its value is known:
                        # move it to the right-hand side
                        rhs[k] -= asf * 2.0 * Q * self._uniform_inflow_shape[i] / dxi2
                    else:
                        rows.append(k); cols.append(self._index(i, j - 1)); vals.append(asf / dxi2)

                rows.append(k); cols.append(k); vals.append(diag)

                # --- source ------------------------------------------------
                # The equation is
                #     d_phi(a_phi dPsi/dphi) + d_xi(a_xi dPsi/dxi)
                #         = -d_phi(B_phi) - d_xi(B_xi),
                # so the buoyancy divergence enters the right-hand side with a
                # MINUS.  Note that every configuration with uniform c and
                # beta = 0 has B constant and hence zero source, which is why
                # this sign is invisible in the single-fluid and flat-interface
                # checks and only shows up on an inclined two-fluid steady
                # state.  See assumptions.md NUM-10.
                rhs[k] -= (B_phi[i, j] - B_phi[i - 1, j]) / self.dphi
                bn = B_xi[i, j] if j < self.n_xi else B_xi[i, self.n_xi - 1]
                bs_ = B_xi[i, j - 1] if j > 0 else B_xi[i, 0]
                rhs[k] -= (bn - bs_) / self.dxi

        A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
        return A, rhs

    def assemble(self, a_phi, a_xi, B_phi, B_xi, Q):
        """Vectorised equivalent of `assemble_reference`; see its docstring."""
        n_phi, n_xi = self.n_phi, self.n_xi
        n_i, n_j = n_phi - 1, n_xi + 1
        dphi2, dxi2 = self.dphi ** 2, self.dxi ** 2
        psi_narrow = 2.0 * Q
        uniform = self.inflow == "uniform"

        idx = (np.arange(n_i)[:, None] * n_j + np.arange(n_j)[None, :])
        ii = np.arange(1, n_phi)[:, None] * np.ones((1, n_j), dtype=int)
        jj = np.ones((n_i, 1), dtype=int) * np.arange(n_j)[None, :]

        ae = a_phi[1:n_phi, :]                       # face between i and i+1
        aw = a_phi[0:n_phi - 1, :]                   # face between i-1 and i
        an = np.concatenate([a_xi[1:n_phi, :], np.zeros((n_i, 1))], axis=1)
        asf = np.concatenate([np.zeros((n_i, 1)), a_xi[1:n_phi, :]], axis=1)

        # ---- right-hand side, in the SAME order as assemble_reference ------
        rhs = np.zeros((n_i, n_j))
        rhs[-1, :] -= ae[-1, :] * psi_narrow / dphi2        # narrow-side Dirichlet
        interior = jj != 0 if uniform else np.ones_like(jj, dtype=bool)
        dirich = (2.0 * Q * self._uniform_inflow_shape[1:n_phi]
                  if uniform else None)
        if uniform:
            # the j = 1 rows lose their south NEIGHBOUR to the right-hand side
            # exact same association as assemble_reference, so the two
            # agree bit for bit rather than to a tolerance
            rhs[:, 1] -= (asf[:, 1] * 2.0 * Q
                          * self._uniform_inflow_shape[1:n_phi] / dxi2)
        rhs -= (B_phi[1:n_phi, :] - B_phi[0:n_phi - 1, :]) / self.dphi
        bn = np.concatenate([B_xi[1:n_phi, :], B_xi[1:n_phi, n_xi - 1:n_xi]], axis=1)
        bs_ = np.concatenate([B_xi[1:n_phi, 0:1], B_xi[1:n_phi, :]], axis=1)
        rhs -= (bn - bs_) / self.dxi
        if uniform:
            rhs[:, 0] = dirich

        rows, cols, vals = [], [], []

        def add(mask, col, val):
            m = mask & interior
            rows.append(idx[m]); cols.append(col[m]); vals.append(val[m])

        add(np.ones_like(jj, dtype=bool), idx,
            -(ae + aw) / dphi2 - (an + asf) / dxi2)
        add(ii + 1 != n_phi, idx + n_j, ae / dphi2)
        add(ii - 1 != 0, idx - n_j, aw / dphi2)
        add(jj + 1 <= n_xi, idx + 1, an / dxi2)
        add(jj - 1 >= 0 if not uniform else jj - 1 >= 1, idx - 1, asf / dxi2)

        if uniform:
            rows.append(idx[:, 0]); cols.append(idx[:, 0])
            vals.append(np.ones(n_i))

        A = sp.csr_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n_i * n_j, n_i * n_j))
        return A, rhs.ravel()

    # -- solve -----------------------------------------------------------------
    def solve(self, c, Q=1.0, umag_pf=None, umag_xf=None):
        a_phi, a_xi, B_phi, B_xi = self.coefficients(c, umag_pf, umag_xf)
        A, rhs = self.assemble(a_phi, a_xi, B_phi, B_xi, Q)
        x = spla.spsolve(A.tocsc(), rhs)
        psi = np.zeros((self.n_phi + 1, self.n_xi + 1))
        psi[-1, :] = 2.0 * Q
        psi[1:-1, :] = x.reshape(self.n_phi - 1, self.n_xi + 1)
        return psi

    def solve_nonlinear(self, c, Q=1.0, tol=1e-10, max_iter=100, relax=1.0,
                        psi0=None, relax_min=0.05):
        """
        Picard iteration for velocity-dependent (non-Newtonian) closures.

        `psi0` seeds the iteration.  In a time loop the previous step's Psi is
        an excellent seed -- c moves by at most one CFL-limited step -- and
        BF25 A.2.3 suggests exactly this continuation.  Measured on the K-GEP-1
        Herschel-Bulkley pair: from a cold start the residual falls only
        linearly and 40 iterations reach 1.8e-5; warm-started it converges in a
        handful.  Over the ~10^5 steps of a full run that is the difference
        between hours and days.
        """
        if self.closures.is_linear:
            return self.solve(c, Q), 1, 0.0
        psi = self.solve(c, Q) if psi0 is None else np.asarray(psi0, dtype=float)
        err = np.inf
        prev = np.inf
        for it in range(1, max_iter + 1):
            u_pf, u_xf = self.speed_on_faces(psi)
            new = self.solve(c, Q, umag_pf=u_pf, umag_xf=u_xf)
            if relax != 1.0:
                new = relax * new + (1.0 - relax) * psi
            err = float(np.max(np.abs(new - psi)) / max(1.0, np.max(np.abs(new))))
            psi = new
            if err < tol:
                return psi, it, err
            # Adaptive under-relaxation.  The Picard gain grows with the
            # buoyancy source, so on a strongly buoyant run the iteration count
            # creeps up as the front develops -- measured on K-GEP-1 it went
            # 4 -> 26 over the first 15% of the run and was heading for the
            # iteration cap, which would have killed the run after hours.
            # Halving the step whenever the residual stops falling bounds it.
            if err >= prev:
                relax = max(0.5 * relax, relax_min)
            prev = err
        return psi, max_iter, err

    # -- reconstruction ---------------------------------------------------------
    def velocities(self, psi):
        """
        Gap-averaged velocities at CELL CENTRES, from
            dPsi/dphi = 2 H r_a w_bar ,   dPsi/dxi = -2 H v_bar .
        """
        g = self.geom.grid
        dpsi_dphi = 0.5 * ((psi[1:, :-1] - psi[:-1, :-1])
                           + (psi[1:, 1:] - psi[:-1, 1:])) / self.dphi
        dpsi_dxi = 0.5 * ((psi[:-1, 1:] - psi[:-1, :-1])
                          + (psi[1:, 1:] - psi[1:, :-1])) / self.dxi
        H_c = self.geom.H(g.phi_centres, g.xi_centres)
        ra_c = self.geom.r_a(g.xi_centres)[None, :]
        w = dpsi_dphi / (2.0 * H_c * ra_c)
        v = -dpsi_dxi / (2.0 * H_c)
        return v, w

    def speed_on_faces(self, psi):
        """|u_bar| = |grad_a Psi| / (2H) on each face family, for the
        velocity-dependent closures."""
        dpdphi_pf = (psi[1:, :] - psi[:-1, :]) / self.dphi           # (n_phi, n_xi+1)
        dpdxi_pf = np.gradient(0.5 * (psi[1:, :] + psi[:-1, :]), self.dxi, axis=1)
        mag_pf = np.hypot(dpdphi_pf / self.ra_pf, dpdxi_pf) / (2.0 * self.H_pf)

        dpdxi_xf = (psi[:, 1:] - psi[:, :-1]) / self.dxi             # (n_phi+1, n_xi)
        dpdphi_xf = np.gradient(0.5 * (psi[:, 1:] + psi[:, :-1]), self.dphi, axis=0)
        mag_xf = np.hypot(dpdphi_xf / self.ra_xf, dpdxi_xf) / (2.0 * self.H_xf)
        return mag_pf, mag_xf

    def axial_flux(self, psi):
        """int_0^1 2 H r_a w_bar dphi at every xi node.  Must equal 2Q."""
        return psi[-1, :] - psi[0, :]
