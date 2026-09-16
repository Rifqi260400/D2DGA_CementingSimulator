"""
M5 -- transport of the gap-averaged concentration.

Solves BF25 (2.19),

    d[H r_a c_bar]/dt + div_a . q = 0 ,    div_a = ((1/r_a) d/dphi, d/dxi)

with the D2DGA volumetric flux BF25 (2.21) = BCF25 (5),

    q = r_a H (v_bar, w_bar) q0(c_bar) + b H^3 r_a I3(c_bar) (-sin beta sin pi phi,
                                                              cos beta)

using the first-order local Lax-Friedrichs (LLF) monotone finite-volume scheme
of BCF25 Section III.


Why this is written as a plain conservation law
-----------------------------------------------
div_a is NOT a plain divergence -- it carries 1/r_a on the azimuthal part --
so at first sight the equation is not in conservation form.  But r_a = r_a(xi)
only, so

    (1/r_a) d(q_phi)/dphi = d(q_phi / r_a)/dphi

exactly.  Writing F = q_phi / r_a and G = q_xi gives

    d[H r_a c]/dt + dF/dphi + dG/dxi = 0                              (T1)
    F = -(1/2) q0 dPsi/dxi  -  b H^3 I3 sin beta sin(pi phi)          (T2)
    G = +(1/2) q0 dPsi/dphi +  b H^3 r_a I3 cos beta                  (T3)

which IS a plain conservation law for the density m = H r_a c_bar.  A finite
volume scheme on m therefore conserves sum(H r_a c_bar) to machine precision,
which is gate M5-T1.  (T2), (T3) use BF25's streamfunction convention
grad_a Psi = (2 H w_bar, -2 H v_bar), i.e. dPsi/dxi = -2 H v_bar and
dPsi/dphi = 2 H r_a w_bar -- CONV-01.  Substituting those recovers
F = r_a H v_bar q0 / r_a and G = r_a H w_bar q0 + ..., i.e. BCF25 (5) exactly.

K = 2 throughout (build spec Section 4.1): fluid 1 is the displaced mud, fluid 2
the displacing cement, c_bar = c_bar_2 is a scalar.  The other fluid follows
from q_{1,0} = 1 - q0 and I_{1,3} = -I3, so BCF25's sum-to-one bookkeeping is
satisfied identically rather than enforced.


Two departures from BCF25 as printed
------------------------------------
NUM-15 (mesh spacings in the buoyancy flux).  BCF25 (30) and (31) add the
buoyancy terms to h_phi and h_xi with no mesh spacing, while the advective
terms enter through DIFFERENCES of Psi, which already carry the spacing
(Psi_{i,j} - Psi_{i,j+1} ~ dPsi/dxi * dxi).  Since the update multiplies
everything by the single factor lambda = dt/(dphi dxi), the printed buoyancy
terms are short by dxi on a phi-face and by dphi on a xi-face; they would make
the buoyancy flux grid-dependent and the scheme inconsistent.  The spacings are
restored here.  The same omission propagates to (34) and (35), and so to the
CFL condition (44), where it is restored identically -- consistently, so the
monotonicity argument is unchanged.

NUM-16 (sign of the printed coefficients (18) and (20)).  BCF25 (11)/(14) and
(18)/(20) are mutually inconsistent: the Psi-differences in the azimuthal
coefficients a_E, a_W appear with the opposite sign to what (11) and (14)
imply.  The axial pair (17)/(19) is consistent with (12)/(14), which pins the
convention: taking (18)/(20) as printed upwinds the azimuthal flux from the
DOWNSTREAM side.  (11)/(14) are the correct pair and are what is implemented.
The check is in tests/test_m5_transport.py::test_num16_azimuthal_upwind_direction.


Monotonicity
------------
The scheme is monotone in the strict sense: c^{n+1}_P is a non-decreasing
function of every entry of the five-point stencil at time n.  The argument
(BCF25 Section III D) needs three ingredients, all of which are structural
here rather than assumed:

  * q0(0) = 0, q0(1) = 1, I3(0) = I3(1) = 0 exactly.  The grouping of the
    Newtonian denominator (see gapscale/newtonian.py) is what makes q0(1) = 1
    exact rather than 1 - 2e-16.
  * the LLF coefficients bound the flux-function slopes evaluated AT THE
    NEIGHBOURING CELL VALUES.  Because each neighbour enters the update through
    its own value only, endpoint evaluation is sufficient -- no interval
    maximum is needed.
  * the buoyancy coefficient [b H^3 ...] is evaluated at CELL CENTRES, not on
    the face.  That is what makes the centre cell's own flux contributions
    cancel between its two opposite faces, leaving F_C with no flux term at
    all.  Evaluating it on the face instead would break the cancellation.

Together with F(0,0,0,0,0) = 0 and F(1,1,1,1,1) = 1 -- the latter because the
discrete Psi-differences telescope to zero EXACTLY around any cell -- this
gives the discrete maximum principle 0 <= c_bar <= 1 (gate M5-T2).
"""

from __future__ import annotations

import numpy as np

from .elliptic import ClosureProvider
from .geometry import Geometry


class TransportSolver:
    """
    First-order LLF finite-volume advance of the cell-centred concentration.

    Parameters
    ----------
    geometry : Geometry
    closures : ClosureProvider          -- supplies q0, I3 and their c-derivatives
    buoyancy_number : float             -- BF25 b = -Delta_rho / Fr*^2 (CONV-04:
                                           the r_a-FREE scalar; the r_a that
                                           belongs to the axial flux is carried
                                           explicitly in (T3))
    inflow_concentration : float        -- c_bar of the fluid entering at xi = 0
    cfl : float                         -- multiplier on the BCF25 (44) bound
    """

    def __init__(self, geometry: Geometry, closures: ClosureProvider,
                 buoyancy_number: float = 0.0,
                 inflow_concentration: float = 1.0,
                 cfl: float = 0.5):
        if not 0.0 < cfl <= 1.0:
            raise ValueError("cfl must lie in (0, 1]")
        if not 0.0 <= inflow_concentration <= 1.0:
            raise ValueError("inflow_concentration must lie in [0, 1]")

        self.geom = geometry
        self.closures = closures
        self.b = float(buoyancy_number)
        self.c_in = float(inflow_concentration)
        self.cfl = float(cfl)

        g = geometry.grid
        self.n_phi, self.n_xi = g.n_phi, g.n_xi
        self.dphi, self.dxi = g.dphi, g.dxi

        # ---- cell-centred geometry ----------------------------------------
        self.H_c = geometry.H(g.phi_centres, g.xi_centres)        # (n_phi, n_xi)
        self.ra_c = np.broadcast_to(geometry.r_a(g.xi_centres)[None, :],
                                    self.H_c.shape).copy()
        self.Hra = self.H_c * self.ra_c
        beta = geometry.well.inclination_rad
        sin_pi_phi = np.sin(np.pi * g.phi_centres)[:, None]

        # Coefficients multiplying I3 in (T2) and (T3), at cell centres.
        self.k_phi = -self.b * self.H_c ** 3 * np.sin(beta) * sin_pi_phi
        self.k_xi = self.b * self.H_c ** 3 * self.ra_c * np.cos(beta)

        self.cell_area = self.dphi * self.dxi

    # ------------------------------------------------------------------
    # state at cell centres
    # ------------------------------------------------------------------
    def _state(self, c, umag):
        """q0, I3 and dq0/dc, dI3/dc at every cell centre."""
        c = np.asarray(c, dtype=float)
        if c.shape != (self.n_phi, self.n_xi):
            raise ValueError(
                f"concentration has shape {c.shape}, expected "
                f"{(self.n_phi, self.n_xi)} (cell centres)")
        q0, I3 = self.closures.q0_I3(c, self.H_c, umag)
        dq0, dI3 = self.closures.dq0_dI3_dc(c, self.H_c, umag)
        return (np.broadcast_to(q0, c.shape), np.broadcast_to(I3, c.shape),
                np.broadcast_to(dq0, c.shape), np.broadcast_to(dI3, c.shape))

    # -- ghost padding in xi only: inflow at xi = 0, outflow at xi = Z ----
    def _pad_xi(self, a, inflow_value=None):
        """Replicate the edge column, except that the inflow column may be
        overridden (Dirichlet inflow state)."""
        left = np.full((a.shape[0], 1), inflow_value) if inflow_value is not None \
            else a[:, :1]
        return np.concatenate([left, a, a[:, -1:]], axis=1)

    # ------------------------------------------------------------------
    # face fluxes
    # ------------------------------------------------------------------
    def fluxes(self, c, psi, umag=None):
        """
        Face-INTEGRATED LLF fluxes.

        Returns
        -------
        Phi : (n_phi + 1, n_xi)   integral of F over each phi-face
        Xi  : (n_phi, n_xi + 1)   integral of G over each xi-face
        a_phi, a_xi : the LLF coefficients on the same faces, for the CFL bound
        """
        q0, I3, dq0, dI3 = self._state(c, umag)
        c = np.asarray(c, dtype=float)

        # ---------------- phi-faces ------------------------------------
        # Face i lies between cell i-1 (L) and cell i (R).  Faces 0 and n_phi
        # are the symmetry planes and carry no flux; pad by replication so the
        # arrays line up, then zero them.
        def pad_phi(a):
            return np.concatenate([a[:1, :], a, a[-1:, :]], axis=0)

        cP, q0P, I3P, dq0P, dI3P = map(pad_phi, (c, q0, I3, dq0, dI3))
        kphiP = pad_phi(self.k_phi)

        dpsi_xi = psi[:, :-1] - psi[:, 1:]              # (n_phi+1, n_xi)

        cL, cR = cP[:-1, :], cP[1:, :]
        adv = 0.25 * (q0P[:-1, :] + q0P[1:, :]) * dpsi_xi
        buoy = 0.5 * self.dxi * (kphiP[:-1, :] * I3P[:-1, :]
                                 + kphiP[1:, :] * I3P[1:, :])

        a_q = 0.5 * np.abs(dpsi_xi) * np.maximum(np.abs(dq0P[:-1, :]),
                                                 np.abs(dq0P[1:, :]))
        a_b = self.dxi * np.maximum(np.abs(kphiP[:-1, :]), np.abs(kphiP[1:, :])) \
            * np.maximum(np.abs(dI3P[:-1, :]), np.abs(dI3P[1:, :]))
        a_phi = a_q + a_b

        Phi = adv + buoy - 0.5 * a_phi * (cR - cL)
        # symmetry planes: no flux, and no dissipation either (so they drop out
        # of the CFL sum as well).  The advective part is already exactly zero
        # there because Psi is constant along both phi boundaries.
        Phi[0, :] = 0.0
        Phi[-1, :] = 0.0
        a_phi = a_phi.copy()
        a_phi[0, :] = 0.0
        a_phi[-1, :] = 0.0

        # ---------------- xi-faces -------------------------------------
        # Face j lies between cell j-1 (S) and cell j (N).
        cX = self._pad_xi(c, self.c_in)
        q0_in, I3_in = self.closures.q0_I3(
            np.full((self.n_phi, 1), self.c_in), self.H_c[:, :1],
            None if umag is None else np.asarray(umag)[:, :1])
        dq0_in, dI3_in = self.closures.dq0_dI3_dc(
            np.full((self.n_phi, 1), self.c_in), self.H_c[:, :1],
            None if umag is None else np.asarray(umag)[:, :1])
        q0X = np.concatenate([np.broadcast_to(q0_in, (self.n_phi, 1)),
                              q0, q0[:, -1:]], axis=1)
        I3X = np.concatenate([np.broadcast_to(I3_in, (self.n_phi, 1)),
                              I3, I3[:, -1:]], axis=1)
        dq0X = np.concatenate([np.broadcast_to(dq0_in, (self.n_phi, 1)),
                               dq0, dq0[:, -1:]], axis=1)
        dI3X = np.concatenate([np.broadcast_to(dI3_in, (self.n_phi, 1)),
                               dI3, dI3[:, -1:]], axis=1)
        kxiX = self._pad_xi(self.k_xi)

        dpsi_phi = psi[1:, :] - psi[:-1, :]             # (n_phi, n_xi+1)

        cS, cN = cX[:, :-1], cX[:, 1:]
        adv_x = 0.25 * (q0X[:, :-1] + q0X[:, 1:]) * dpsi_phi
        buoy_x = 0.5 * self.dphi * (kxiX[:, :-1] * I3X[:, :-1]
                                    + kxiX[:, 1:] * I3X[:, 1:])

        b_q = 0.5 * np.abs(dpsi_phi) * np.maximum(np.abs(dq0X[:, :-1]),
                                                  np.abs(dq0X[:, 1:]))
        b_b = self.dphi * np.maximum(np.abs(kxiX[:, :-1]), np.abs(kxiX[:, 1:])) \
            * np.maximum(np.abs(dI3X[:, :-1]), np.abs(dI3X[:, 1:]))
        a_xi = b_q + b_b

        Xi = adv_x + buoy_x - 0.5 * a_xi * (cN - cS)
        return Phi, Xi, a_phi, a_xi

    # ------------------------------------------------------------------
    # timestep
    # ------------------------------------------------------------------
    def max_timestep(self, a_phi, a_xi):
        """BCF25 (44), with the mesh spacings restored (NUM-15).  Returns the
        LARGEST monotonicity-preserving step; the caller applies the CFL
        number."""
        total = (a_phi[:-1, :] + a_phi[1:, :]
                 + a_xi[:, :-1] + a_xi[:, 1:])
        with np.errstate(divide="ignore"):
            k = np.where(total > 0.0, 2.0 * self.Hra / np.maximum(total, 1e-300),
                         np.inf)
        return float(self.cell_area * np.min(k))

    # ------------------------------------------------------------------
    # one step
    # ------------------------------------------------------------------
    def step(self, c, psi, dt=None, umag=None, max_dt=None):
        """
        Advance one timestep.  If dt is None the CFL-limited step is used,
        capped at max_dt when one is given (so a run can land exactly on a
        requested output time without a second flux evaluation).

        Returns (c_new, dt, boundary_flux) where boundary_flux is the net
        volume of fluid 2 that entered through xi = 0 minus that which left
        through xi = Z, so that

            mass(c_new) - mass(c) == boundary_flux

        holds to machine precision.  That identity is gate M5-T1.
        """
        Phi, Xi, a_phi, a_xi = self.fluxes(c, psi, umag)
        dt_max = self.max_timestep(a_phi, a_xi)
        if dt is None:
            dt = self.cfl * dt_max
            if max_dt is not None:
                dt = min(dt, float(max_dt))
        elif dt > dt_max * (1.0 + 1e-12):
            raise ValueError(
                f"dt = {dt:.6g} exceeds the monotonicity limit {dt_max:.6g} "
                f"(BCF25 44); the maximum principle is not guaranteed.")

        lam = dt / self.cell_area
        div = (Phi[1:, :] - Phi[:-1, :]) + (Xi[:, 1:] - Xi[:, :-1])
        c_new = np.asarray(c, dtype=float) - lam * div / self.Hra

        influx = dt * float(np.sum(Xi[:, 0]))
        outflux = dt * float(np.sum(Xi[:, -1]))
        return c_new, dt, influx - outflux

    # ------------------------------------------------------------------
    def mass(self, c):
        """Total volume of fluid 2 in the domain (half annulus, scaled)."""
        return float(np.sum(self.Hra * np.asarray(c, dtype=float))
                     * self.cell_area)
