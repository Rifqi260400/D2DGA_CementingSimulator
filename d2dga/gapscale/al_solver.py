"""
M2 -- gap-scale two-layer flow, solved by augmented Lagrangian.

BF25 Appendix A.  The unit-channel problem: fluid 2 (displacing) occupies the
centre of the gap, y in [0, c), fluid 1 (displaced) the wall layer, y in (c, 1].
Symmetry at y = 0, no slip at y = 1, stress and velocity continuous at y = c.

    0 = d tau2/dy + G - (1-c) Gb ,   y in [0, c)      (A1)
    0 = d tau1/dy + G +   c   Gb ,   y in (c, 1]      (A2)

    tau_k = (kappa_k |du/dy|^(n_k - 1) + tau_Y,k |du/dy|^-1) du/dy,
                                          |tau_k| >  tau_Y,k
    du/dy = 0,                            |tau_k| <= tau_Y,k        (A3)

u is a TWO-COMPONENT vector, (phi, xi).  This is not decoration: the buoyancy
vector acts in opposite directions in the two layers (BF25 after A12), so the
shear-stress direction rotates across the gap and the components do not
decouple.  BF25 Fig. 4(b) shows the resulting kink at the interface.

The yield stress makes the minimisation (A12) non-differentiable, so we relax
du/dy -> q, enforce it with a multiplier, and augment.  The saddle-point problem
(A14)-(A15) is driven by an Uzawa iteration.

Two closure scenarios, both from BF25 A.2:
  * `solve_fixed_pressure_gradient`  (A.2.1) -- G given, find u.
  * `solve_fixed_mean_velocity`      (A.2.2) -- u-bar given (i.e. grad_a Psi
    from the 2-D model), find G along with u.  This is the one the D2DGA
    closures need.

Everything here is in TILDE variables: the unit channel, BF25 (A4)

    y~ = y/H ,  kappa~_k = kappa_k / H^n_k ,  G~ = H G ,  Gb~ = H Gb

with velocity, yield stresses and power-law indices unchanged.  That rescaling
is what lets ONE closure table serve every depth (build spec Section 3.2).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..scaling import ScaledFluid


@dataclass
class GapSolution:
    """Converged (or not) state of the unit-channel problem."""
    c: float
    y_edges: np.ndarray       # (N+1,)
    y_centres: np.ndarray     # (N,)
    dy: np.ndarray            # (N,)
    in_fluid2: np.ndarray     # (N,) bool, True where y < c
    u: np.ndarray             # (N+1, 2) velocity at edges
    dudy: np.ndarray          # (N, 2) at cell centres
    fluidity: np.ndarray      # (N,) 1/eta at centres; 0 where unyielded
    G: np.ndarray             # (2,) modified pressure gradient
    Gb: np.ndarray            # (2,) buoyancy vector
    mean_velocity: np.ndarray  # (2,)
    converged: bool
    iterations: int
    residual: float

    @property
    def unyielded_fraction(self) -> float:
        return float(np.sum(self.dy[self.fluidity == 0.0]))


class TwoLayerGapSolver:
    """
    Augmented Lagrangian solver for the unit-channel two-layer problem.

    Numerical parameters r and rho:  PF04 (27) requires 0 < rho < r(1+sqrt5)/2
    for convergence; BF25 A.2.3 reports using rho = r = 1 and neither paper
    gives a tolerance.  We default to rho = r = 1, which satisfies the bound
    with a wide margin (1 < 1.618).  See assumptions.md NUM-02.
    """

    def __init__(self, fluid1: ScaledFluid, fluid2: ScaledFluid,
                 n_y: int = 200, r: float = 1.0, rho: float = 1.0,
                 tol: float = 1e-10, max_iter: int = 20000,
                 bisection_steps: int = 60):
        self.f1 = fluid1
        self.f2 = fluid2
        self.n_y = int(n_y)
        self.r = float(r)
        self.rho = float(rho)
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.bisection_steps = int(bisection_steps)
        if not (0.0 < self.rho < self.r * (1.0 + np.sqrt(5.0)) / 2.0):
            raise ValueError("rho, r violate the PF04 (27) convergence bound")
        # norm exponent, BF25 A.2.3
        self.p = 1.0 + min(fluid1.power_law_index, fluid2.power_law_index)

    # -- grid ---------------------------------------------------------------
    def _grid(self, c):
        """Two uniform sub-grids meeting exactly at the interface y = c, so the
        jump in fluid properties always sits on a cell face."""
        if c <= 0.0:
            y_e = np.linspace(0.0, 1.0, self.n_y + 1)
            n2 = 0
        elif c >= 1.0:
            y_e = np.linspace(0.0, 1.0, self.n_y + 1)
            n2 = self.n_y
        else:
            n2 = int(np.clip(round(self.n_y * c), 1, self.n_y - 1))
            n1 = self.n_y - n2
            y_e = np.concatenate([np.linspace(0.0, c, n2 + 1),
                                  np.linspace(c, 1.0, n1 + 1)[1:]])
        y_c = 0.5 * (y_e[:-1] + y_e[1:])
        dy = np.diff(y_e)
        in_f2 = np.arange(y_c.size) < n2
        return y_e, y_c, dy, in_f2

    # -- per-cell fluid properties -----------------------------------------
    def _props(self, in_f2):
        kappa = np.where(in_f2, self.f2.consistency, self.f1.consistency)
        n = np.where(in_f2, self.f2.power_law_index, self.f1.power_law_index)
        tau_y = np.where(in_f2, self.f2.yield_stress, self.f1.yield_stress)
        return kappa, n, tau_y

    # -- q update, BF25 (A20)/(A29) ----------------------------------------
    def _update_q(self, m_vec, kappa, n, tau_y):
        """
        q = 0                      if |m| <= tau_Y
        q = theta m, with  kappa theta^n |m|^n + r theta |m| = |m| - tau_Y

        The left side is strictly increasing in theta from 0, so the root is
        unique and bracketed by [0, (|m|-tau_Y)/(r|m|)]: at that upper end the
        r-term alone already equals the right side.  Vectorised bisection --
        slower than Newton but it cannot fail, and M2 is a correctness gate.
        """
        mag = np.linalg.norm(m_vec, axis=1)
        yielded = mag > tau_y
        theta = np.zeros_like(mag)
        if not np.any(yielded):
            return np.zeros_like(m_vec), yielded

        mg = mag[yielded]
        k, nn, ty = kappa[yielded], n[yielded], tau_y[yielded]
        rhs = mg - ty
        lo = np.zeros_like(mg)
        hi = rhs / (self.r * mg)
        for _ in range(self.bisection_steps):
            mid = 0.5 * (lo + hi)
            f = k * (mid * mg) ** nn + self.r * mid * mg - rhs
            too_big = f > 0.0
            hi = np.where(too_big, mid, hi)
            lo = np.where(too_big, lo, mid)
        theta[yielded] = 0.5 * (lo + hi)
        return theta[:, None] * m_vec, yielded

    # -- helpers ------------------------------------------------------------
    @staticmethod
    def _integrate_from_wall(dudy, dy):
        """u(1) = 0, u_i = u_{i+1} - dy_i * (du/dy)_i."""
        n = dudy.shape[0]
        u = np.zeros((n + 1, dudy.shape[1]))
        for i in range(n - 1, -1, -1):
            u[i] = u[i + 1] - dy[i] * dudy[i]
        return u

    @staticmethod
    def _mean(u, dy):
        return np.sum(0.5 * (u[:-1] + u[1:]) * dy[:, None], axis=0)

    # NOTE on the constant of integration, BF25 (A18) / (A26).
    #
    # Those equations read
    #     r du/dy = r q(y) - lam~(y) - r q(0) + lam~(0),
    # i.e. the bracket -r q(0) + lam~(0) is the constant that imposes symmetry,
    # du/dy = 0 at y = 0.  With q and lam~ stored at CELL CENTRES neither is
    # available at y = 0, and the obvious move -- extrapolate from the first two
    # centres -- is wrong, in a way that quietly costs an order of accuracy.
    #
    # Both q (which approximates du/dy) and lam~ (which approximates tau - lam0)
    # are ODD about the channel centre, because u is even there.  An odd function
    # vanishes at 0, so the correct value of the bracket is exactly ZERO.  Linear
    # extrapolation does not know that: applied to du/dy ~ y^(1/n) it returns
    # -0.49 dy^(1/n) instead of 0, and the fixed point then absorbs the mismatch
    # as a CONSTANT offset in lam~, which shifts the whole stress field.
    #
    # For a power-law fluid the damage is O(dy), not O(dy^(1/n)): solving
    #     1.5 du/dy(dy/2) = 0.5 du/dy(3dy/2)   with   tau = -G y + lam~
    # gives lam~ = -0.5717 G dy for n = 0.6.  Measured: -0.01143 at G = 4,
    # dy = 0.005; predicted: -0.01143.  That single constant was responsible for
    # a 1.6e-3 error in u_bar and for first-order mesh convergence, against the
    # second order BF25 A.2.3 reports.
    #
    # Setting the bracket to zero restores tau = lam0 exactly at the fixed point
    # (hence lam~ -> 0) and recovers second order.  See assumptions.md NUM-08.

    def _lambda0_fixed_velocity(self, y_c, c, Gb):
        """BF25 (A23): buoyancy-balancing multiplier, velocity-driven case."""
        lam = np.empty((y_c.size, 2))
        inner = y_c < c
        lam[inner] = (1.0 - c) * Gb[None, :] * y_c[inner, None]
        lam[~inner] = (1.0 - y_c[~inner, None]) * Gb[None, :] * c
        return lam

    def _lambda0_fixed_gradient(self, y_c, c, G, Gb):
        """BF25 (A17): as above but also absorbing the pressure gradient."""
        return self._lambda0_fixed_velocity(y_c, c, Gb) - G[None, :] * y_c[:, None]

    def _fluidity(self, dudy, kappa, n, tau_y):
        """1/eta at each cell.  Zero where unyielded, which is exactly what the
        closure integrals need -- no infinities ever appear."""
        g = np.linalg.norm(dudy, axis=1)
        denom = kappa * np.power(g, n, where=g > 0, out=np.zeros_like(g)) + tau_y
        out = np.zeros_like(g)
        ok = (g > 0) & (denom > 0)
        out[ok] = g[ok] / denom[ok]
        return out

    # -- main solves --------------------------------------------------------
    def solve_fixed_mean_velocity(self, c, u_mean, Gb, u_init=None):
        """BF25 A.2.2.  u_mean and Gb are 2-vectors."""
        c = float(c)
        u_mean = np.asarray(u_mean, dtype=float).reshape(2)
        Gb = np.asarray(Gb, dtype=float).reshape(2)
        y_e, y_c, dy, in_f2 = self._grid(c)
        kappa, n, tau_y = self._props(in_f2)
        lam0 = self._lambda0_fixed_velocity(y_c, c, Gb)

        # unit-pressure-gradient Poiseuille reference, BF25 (A22)
        dudy_P = -y_c[:, None] / self.r * np.ones((1, 2))
        u_P = self._integrate_from_wall(dudy_P, dy)
        mean_P = self._mean(u_P, dy)

        q = np.zeros((y_c.size, 2)) if u_init is None else u_init
        lam_t = np.zeros((y_c.size, 2))
        G = np.zeros(2)
        residual, converged, it = np.inf, False, 0

        for it in range(1, self.max_iter + 1):
            dudy_I = q - lam_t / self.r          # symmetry constant == 0, see above
            u_I = self._integrate_from_wall(dudy_I, dy)
            mean_I = self._mean(u_I, dy)

            alpha = (u_mean - mean_I) / mean_P      # BF25 (A28)
            G = alpha
            dudy = dudy_I + alpha[None, :] * dudy_P
            u = u_I + alpha[None, :] * u_P

            m_vec = lam0 + lam_t + self.r * dudy
            q_new, _ = self._update_q(m_vec, kappa, n, tau_y)
            lam_t = lam_t + self.rho * (dudy - q_new)

            # Two residuals.  Stagnation ||q^{k+1}-q^k|| alone can report success
            # while the AL constraint du/dy = q is still violated, so we also
            # require feasibility.  PF04 monitors the former; the latter is what
            # the method is actually solving.
            stagnation = float(np.sum(np.abs(q_new - q) ** self.p
                                      * dy[:, None]) ** (1.0 / self.p))
            feasibility = float(np.max(np.abs(dudy - q_new)))
            residual = max(stagnation, feasibility)
            q = q_new
            if residual < self.tol:
                converged = True
                break

        fluidity = self._fluidity(dudy, kappa, n, tau_y)
        return GapSolution(c, y_e, y_c, dy, in_f2, u, dudy, fluidity,
                           G, Gb, self._mean(u, dy), converged, it, residual)

    def solve_fixed_pressure_gradient(self, c, G, Gb):
        """BF25 A.2.1.  Used for verification and for inverting the closure."""
        c = float(c)
        G = np.asarray(G, dtype=float).reshape(2)
        Gb = np.asarray(Gb, dtype=float).reshape(2)
        y_e, y_c, dy, in_f2 = self._grid(c)
        kappa, n, tau_y = self._props(in_f2)
        lam0 = self._lambda0_fixed_gradient(y_c, c, G, Gb)

        q = np.zeros((y_c.size, 2))
        lam_t = np.zeros((y_c.size, 2))
        residual, converged, it = np.inf, False, 0
        dudy = np.zeros((y_c.size, 2))

        for it in range(1, self.max_iter + 1):
            dudy = q - lam_t / self.r            # symmetry constant == 0, see above
            u = self._integrate_from_wall(dudy, dy)
            m_vec = lam0 + lam_t + self.r * dudy
            q_new, _ = self._update_q(m_vec, kappa, n, tau_y)
            lam_t = lam_t + self.rho * (dudy - q_new)
            stagnation = float(np.sum(np.abs(q_new - q) ** self.p
                                      * dy[:, None]) ** (1.0 / self.p))
            feasibility = float(np.max(np.abs(dudy - q_new)))
            residual = max(stagnation, feasibility)
            q = q_new
            if residual < self.tol:
                converged = True
                break

        u = self._integrate_from_wall(dudy, dy)
        fluidity = self._fluidity(dudy, kappa, n, tau_y)
        return GapSolution(c, y_e, y_c, dy, in_f2, u, dudy, fluidity,
                           G, Gb, self._mean(u, dy), converged, it, residual)
