"""
M3 -- closure tabulation.

BF25's own conclusion (Section 4): computing the gap-scale closures per cell
per timestep over a wellbore domain is "impractical without, e.g. tabulating
each function".  So this is a first-class module, not an optimisation.

What has to be tabulated, and what does not
-------------------------------------------
For a NEWTONIAN pair the gap problem is linear: eta is constant in each layer,
so I1, I2, q0, I3 depend on c alone (through m).  H drops out by the BF25 (A4)
rescaling, and the velocity magnitude and buoyancy drop out by linearity.  A
1-D table in c is exact, and in fact the analytic forms in `newtonian.py` make
even that unnecessary.

For HERSCHEL-BULKLEY the effective viscosity depends on the local shear rate,
so the closures pick up a dependence on the STATE of the cell.  In tilde
variables (BF25 A4) that state is

    c          gap-averaged concentration of the displacing fluid
    H          local half-gap        -> kappa~_k = kappa_k / H^n_k
    |u_bar|    local mean speed      = |grad_a Psi| / (2H)
    G~_b       buoyancy vector       = H G_b

so the table is 4-D once the rheology of the pair is fixed, which it is for any
one simulation.  We exploit the rotational structure: the gap problem depends on
the two vectors u_bar and G~_b only through their magnitudes and the angle
between them, and for a VERTICAL well (beta = 0) G~_b has a single component
aligned with the axis, so the angle axis collapses.  K-GEP-1 is vertical.

Interpolation
-------------
LINEAR, deliberately.  A cubic or spline interpolant can overshoot, and a
non-monotone q0 would break the LLF transport scheme's maximum principle
(BCF25 Section III B leans on monotone fluxes).  Linear interpolation of
monotone samples is monotone; that is asserted at M3-T2.
"""

from __future__ import annotations

import itertools
import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from ..scaling import ScaledFluid
from .al_solver import TwoLayerGapSolver
from .closures import Closures, closures_from_solution


class ClosureTableRangeWarning(UserWarning):
    """Raised at runtime when a cell asks for a state outside the tabulated box."""


@dataclass
class ClosureTable:
    """
    Pre-computed (I1, I2, q0, I3) over a regular grid in the gap-scale state.

    Axes are given as 1-D arrays.  `c_grid` is mandatory; the other three are
    optional and any that is omitted (or length-1) is treated as a fixed value,
    which is the Newtonian situation.
    """

    fluid1: ScaledFluid
    fluid2: ScaledFluid
    c_grid: np.ndarray
    h_grid: np.ndarray | None = None
    umag_grid: np.ndarray | None = None
    gb_grid: np.ndarray | None = None
    n_y: int = 400
    tol: float = 1e-10
    max_iter: int = 60000

    _values: dict = field(default_factory=dict, init=False, repr=False)
    _interp: dict = field(default_factory=dict, init=False, repr=False)
    _axes: list = field(default_factory=list, init=False, repr=False)
    n_failed: int = field(default=0, init=False)

    # -- axis bookkeeping ---------------------------------------------------
    def _axis_list(self):
        axes = [("c", np.atleast_1d(np.asarray(self.c_grid, float)))]
        for name, g, default in (("H", self.h_grid, 1.0),
                                 ("umag", self.umag_grid, 1.0),
                                 ("gb", self.gb_grid, 0.0)):
            axes.append((name, np.atleast_1d(np.asarray(
                default if g is None else g, float))))
        return axes

    @property
    def shape(self):
        return tuple(len(a) for _, a in self._axis_list())

    @property
    def n_points(self):
        return int(np.prod(self.shape))

    # -- build ---------------------------------------------------------------
    def build(self, verbose: bool = False) -> "ClosureTable":
        self._axes = self._axis_list()
        shape = self.shape
        out = {k: np.zeros(shape) for k in ("I1", "I2", "q0", "I3")}
        solver = TwoLayerGapSolver(self.fluid1, self.fluid2, n_y=self.n_y,
                                   tol=self.tol, max_iter=self.max_iter)
        self.n_failed = 0

        grids = [a for _, a in self._axes]
        for idx in itertools.product(*[range(n) for n in shape]):
            c, H, umag, gb = (grids[k][idx[k]] for k in range(4))
            cl = self._solve_one(solver, c, H, umag, gb)
            for k in out:
                out[k][idx] = getattr(cl, k)
        if verbose:
            print(f"built {self.n_points} points, {self.n_failed} unconverged")

        self._values = out
        pts = tuple(a for _, a in self._axes)
        self._interp = {
            k: RegularGridInterpolator(pts, v, method="linear",
                                       bounds_error=False, fill_value=None)
            for k, v in out.items()
        }
        return self

    def _solve_one(self, solver, c, H, umag, gb) -> Closures:
        """
        One gap-scale solve, in tilde variables.

        BF25 (A4): kappa~_k = kappa_k / H^n_k, G~_b = H G_b, with velocity,
        yield stress and power-law index unchanged.  We therefore rescale the
        CONSISTENCIES by H and leave everything else alone -- which is exactly
        why one table serves every depth.
        """
        if H != 1.0:
            f1 = ScaledFluid(self.fluid1.name,
                             self.fluid1.density,
                             self.fluid1.consistency / H ** self.fluid1.power_law_index,
                             self.fluid1.power_law_index, self.fluid1.yield_stress)
            f2 = ScaledFluid(self.fluid2.name,
                             self.fluid2.density,
                             self.fluid2.consistency / H ** self.fluid2.power_law_index,
                             self.fluid2.power_law_index, self.fluid2.yield_stress)
            solver = TwoLayerGapSolver(f1, f2, n_y=self.n_y, tol=self.tol,
                                       max_iter=self.max_iter)
        sol = solver.solve_fixed_mean_velocity(c, [0.0, umag], [0.0, gb * H])
        if not sol.converged:
            self.n_failed += 1
        return closures_from_solution(sol)

    # -- evaluate -------------------------------------------------------------
    def __call__(self, c, H=1.0, umag=1.0, gb=0.0, warn_out_of_range=True):
        """Interpolate all four closures.  Broadcasting over array inputs."""
        if not self._interp:
            raise RuntimeError("call build() first")
        c, H, umag, gb = np.broadcast_arrays(*[np.asarray(x, float)
                                               for x in (c, H, umag, gb)])
        pts = np.stack([c.ravel(), H.ravel(), umag.ravel(), gb.ravel()], axis=-1)
        if warn_out_of_range:
            self._check_range(pts)
        vals = {k: f(pts).reshape(c.shape) for k, f in self._interp.items()}
        return vals["I1"], vals["I2"], vals["q0"], vals["I3"]

    def _check_range(self, pts):
        """M3-T3: name the offending cell rather than silently extrapolating."""
        for j, (name, axis) in enumerate(self._axes):
            if axis.size == 1:
                continue
            lo, hi = axis[0], axis[-1]
            bad = (pts[:, j] < lo - 1e-12) | (pts[:, j] > hi + 1e-12)
            if np.any(bad):
                k = int(np.argmax(bad))
                warnings.warn(
                    f"closure table: axis '{name}' out of range at flat cell "
                    f"{k}: requested {pts[k, j]:.6g}, table covers "
                    f"[{lo:.6g}, {hi:.6g}]. Linear extrapolation is being used; "
                    f"widen the table or check the solution.",
                    ClosureTableRangeWarning, stacklevel=3)

    # -- diagnostics ----------------------------------------------------------
    def q0_is_monotone_in_c(self) -> bool:
        """M3-T2.  Linear interpolation preserves monotonicity of its samples,
        so checking the samples is sufficient -- and necessary, because a
        non-monotone q0 would break the transport scheme's maximum principle."""
        q = self._values["q0"]
        return bool(np.all(np.diff(q, axis=0) >= -1e-12))

    def endpoint_values(self):
        """q0 must be 0 at c = 0 and 1 at c = 1 for every other axis value."""
        q = self._values["q0"]
        return q[0].ravel(), q[-1].ravel()


def newtonian_table(m, n_c=161, **kw) -> ClosureTable:
    """
    Convenience 1-D table for a Newtonian pair.

    Exact by construction: the Newtonian gap problem is linear, so the closures
    have no H, |u| or Gb dependence at all and one axis suffices.  Used to test
    the interpolation machinery against `newtonian.py` analytics without paying
    for a 4-D build.

    n_c = 161 UNIFORM is the resolution choice (assumptions.md NUM-03), from a
    convergence study at m = 0.2 measuring error against max|f| over the domain:

        n_c    I1        I2        q0        I3
         41    3.7e-4    9.3e-4    7.8e-4    5.2e-3
         81    9.0e-5    2.3e-4    2.3e-4    1.4e-3
        161    2.3e-5    6.0e-5    6.3e-5    3.5e-4   <- chosen
        321    5.9e-6    1.5e-5    1.5e-5    8.5e-5

    I3 is the binding constraint: it vanishes at both endpoints and behaves like
    c^2 near zero, which linear interpolation resolves worst.  161 leaves a 3x
    margin on M3-T1's 0.1% budget.

    Chebyshev clustering was tried and REJECTED: it improves q0 (2.4e-4 vs
    7.8e-4 at n_c=41) but degrades I2 by more than 2x, because I2 peaks in the
    interior where clustering thins the grid.  Uniform wins on the binding term.

    Note on the error norm: pointwise RELATIVE error is meaningless for I2 and
    I3, which vanish at c = 0 and c = 1.  Error is measured against max|f| over
    the domain, which is the scale at which these enter the flux additively.
    """
    f1 = ScaledFluid("displaced", 1.0, np.sqrt(m), 1.0, 0.0)
    f2 = ScaledFluid("displacing", 1.0, 1.0 / np.sqrt(m), 1.0, 0.0)
    return ClosureTable(f1, f2, c_grid=np.linspace(0.0, 1.0, n_c), **kw)
