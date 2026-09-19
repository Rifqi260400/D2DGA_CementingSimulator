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
the two vectors u_bar and G~_b only through their magnitudes and the ANGLE
between them.

B-1.  An earlier version of this paragraph continued "and for a VERTICAL well
(beta = 0) G~_b has a single component aligned with the axis, so the angle axis
collapses".  That is wrong.  G~_b being purely axial does not collapse the
angle, because u_bar is still a 2-vector: theta = angle(u_bar, G~_b) is non-zero
whenever v_bar != 0, which is whenever the front is tilted -- i.e. in every case
of interest.  What the table actually stores is the theta = 0 slice
(`_solve_one` passes u_mean = [0, umag] against G~_b = [0, gb*H]), and that is
an APPROXIMATION, not a reduction.

Its size is set by |u_bar| / (gb * I1), i.e. by how far the pressure-driven
stress is from negligible against the buoyancy-driven one.  Measured on the
K-GEP-1 pair (gb = 27.9, H = 1) as the deviation of the theta = 90 degree solve
from the stored theta = 0 one -- see test_b1_angle_assumption_is_bounded:

    |u_bar|      I1       I2       q0       I3
        0.5   +0.0%    +0.0%    +0.0%    -0.0%
        5     +0.1%    +0.1%    +0.0%    -0.3%
       50     +0.6%    +1.4%    +0.1%    -2.6%

The production run reached |u_bar| = 75, so ~3% on I3 is the operating error.
It grows with |u_bar| and the axis runs to 3000, so the bound is NOT uniform
over the table -- the test measures the top of the axis too and states it.
NUM-14 already refuses beta != 0, where a second angle appears; this note is
about beta = 0, where the approximation is live and was undocumented.

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
    # None = choose r by probing at build time (NUM-02).  Give a number to pin
    # it, e.g. to reproduce an older table exactly.
    r: float | None = None

    tuned_r: float = field(default=float("nan"), init=False)
    _values: dict = field(default_factory=dict, init=False, repr=False)
    _interp: dict = field(default_factory=dict, init=False, repr=False)
    _axes: list = field(default_factory=list, init=False, repr=False)
    n_failed: int = field(default=0, init=False)
    # A-3/B-2.  The range guard used to be a UserWarning and nothing else, so a
    # `-W ignore::UserWarning` on the command line removed the only evidence
    # that a result rested on extrapolated closures -- which is exactly what
    # happened to the K-GEP-1 production runs.  The warning is kept for
    # interactive use, but the RECORD is now state on the table: it is updated
    # on every query, it cannot be switched off, and `range_report()` is
    # printed by the run scripts.  `derivative_bounds` records too; it did not
    # even warn before.
    angle_sensitivity: float = field(default=float("nan"), init=False)
    angle_sensitivity_detail: dict = field(default_factory=dict, init=False,
                                           repr=False)
    n_range_queries: int = field(default=0, init=False)
    n_range_points_out: dict = field(default_factory=dict, init=False)
    worst_excursion: dict = field(default_factory=dict, init=False, repr=False)

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
        gb_axis = self._axes[3][1]
        gb_probe = float(gb_axis[np.argmax(np.abs(gb_axis))])
        common = dict(n_y=self.n_y, tol=self.tol, max_iter=self.max_iter)
        if self.r is None:
            # Probe at the LARGEST |gb| on the table's own axis: that is the
            # stiffest point, and r matters only where buoyancy is strong.
            solver = TwoLayerGapSolver.tuned(self.fluid1, self.fluid2,
                                             gb_probe=gb_probe, **common)
            self.tuned_r = solver.r
        else:
            solver = TwoLayerGapSolver(self.fluid1, self.fluid2, r=self.r,
                                       rho=self.r, **common)
            self.tuned_r = self.r
        self.n_failed = 0

        grids = [a for _, a in self._axes]
        for idx in itertools.product(*[range(n) for n in shape]):
            c, H, umag, gb = (grids[k][idx[k]] for k in range(4))
            cl = self._solve_one(solver, c, H, umag, gb)
            for k in out:
                out[k][idx] = getattr(cl, k)
        self._measure_angle_sensitivity(solver)
        if verbose:
            print(f"built {self.n_points} points, {self.n_failed} unconverged, "
                  f"r = {self.tuned_r:g}")
            print(self.assumption_report())

        self._values = out
        self._dinterp = None
        pts = tuple(a for _, a in self._axes)
        self._interp = {
            k: RegularGridInterpolator(pts, v, method="linear",
                                       bounds_error=False, fill_value=None)
            for k, v in out.items()
        }
        return self

    def _measure_angle_sensitivity(self, solver):
        """
        B-1.  Every entry of this table is built with u_mean and G~_b PARALLEL
        (`_solve_one` passes u_mean = [0, umag] against Gb = [0, gb*H]).  The
        real flow has theta = angle(u_bar, G~_b) != 0 whenever v_bar != 0, i.e.
        whenever the front is tilted.  BF25 (2.8)-(2.10) make the gap-scale
        stress a 2-vector with eta_k = eta_k(|tau_k|), so all four closures
        depend on theta.

        For a pair in which BOTH fluids are Newtonian there is no dependence at
        all -- eta is constant, so only the direction of |tau| changes and the
        closures do not.  Measured: exactly 0.00% at m = 0.0017 and m = 0.033.
        Once EITHER fluid has a yield stress the dependence is not small:
        whether material yields at all is set by |tau|, and |tau| depends on
        theta, so the unyielded fraction -- and therefore I1, an integral of the
        fluidity -- moves sharply.  On a pair with a yield stress in both fluids
        (kappa 0.45/0.25, n 0.6/0.8, tau_Y 0.30/0.12, gb = 25) the theta = 90
        degree closures differ from the stored theta = 0 ones by up to

            I1 +42%,  I2 +134%,  q0 +8.7%,  I3 -107%

        CORRECTION, 2026-09-19.  This docstring used to explain the benign
        K-GEP-1 result by saying "the displaced fluid is Newtonian".  That
        attribution is wrong, and AUDIT_REPORT.md B-1 repeats it.  A sweep
        holding the mud NEWTONIAN (n = 1, tau_Y = 0) and moving only its
        viscosity against the unchanged K-GEP-1 cement gives

            kappa1 = 1 mPa s   m = 0.0063   ->   0.53%   (the shipped default)
            kappa1 = 2 mPa s   m = 0.013    ->   2.0%
            kappa1 = 5 mPa s   m = 0.032    ->  11.4%
            kappa1 = 10 mPa s  m = 0.063    ->  35.7%
            kappa1 = 20 mPa s  m = 0.127    ->  93.1%
            kappa1 = 50 mPa s  m = 0.317    -> 168.6%

        So what governs the size is m, not whether the mud is Newtonian: at
        m ~ 0.006 the mud is ~160x thinner, so the cement's yielded structure --
        the only theta-sensitive part -- barely feels the mud layer.  A perfectly
        Newtonian mud at an ordinary 5 mPa s already puts B-1 at 11%.

        SECOND CORRECTION, same day, and it is the sharper one.  The figures
        above come from a probe whose velocity axis stops at |u_bar| = 10.  The
        sensitivity GROWS with |u_bar|, and the production table's axis reaches
        3000 -- it has to, because FLU-06 drives |u_bar| past 370.  Rebuilt on
        the production axes (31 x 5 x 17, umag to 3000) the SHIPPED K-GEP-1 pair
        measures

            |u| = 0: 0.0%    |u| = 8.3: 0.4%    |u| = 3000: 9.8%

        -- worst 9.8%, which this method's own grading calls SIGNIFICANT, not
        benign.  It had never been measured on those axes before: the production
        runs LOADED a cached table, and a loaded table reports "not measured".
        No computed number moves -- this is a diagnostic -- but "benign for the
        K-GEP-1 pair", as AUDIT_REPORT.md B-1 puts it, is not supported at the
        resolution the runs actually use.  Quote K-GEP-1 numbers with 9.8%
        attached, and note that |u_bar| is largest exactly where the front is
        most tilted, which is where the angle error also peaks.

        Read the number as an upper bound on the angle error, not as the error
        in a run: it compares theta = 0 against theta = 90 degrees, and a run
        whose front stays nearly axial never visits the worst angle.  It bounds
        what the stored slice can be wrong by; it does not say what eta_E is
        wrong by.  Deciding that needs the fifth axis (NUM-14), not this probe.

        This is a diagnostic, not a fix.  The fix is a fifth axis (the angle),
        which NUM-14 already describes for the beta != 0 case.
        """
        probe_c = (0.25, 0.5, 0.75)
        gb_axis = self._axes[3][1]
        gb = float(gb_axis[np.argmax(np.abs(gb_axis))])
        u_axis = self._axes[2][1]
        probe_u = sorted({float(u_axis[0]), float(np.median(u_axis)),
                          float(u_axis[-1])})
        worst, detail = 0.0, {}
        for u in probe_u:
            local = 0.0
            for c in probe_c:
                ref = closures_from_solution(solver.solve_fixed_mean_velocity(
                    c, [0.0, u], [0.0, gb]))
                rot = closures_from_solution(solver.solve_fixed_mean_velocity(
                    c, [u, 0.0], [0.0, gb]))
                for a, b in ((rot.I1, ref.I1), (rot.I2, ref.I2),
                             (rot.q0, ref.q0), (rot.I3, ref.I3)):
                    local = max(local, abs(a - b) / max(abs(b), 1e-30))
            detail[u] = local
            worst = max(worst, local)
        self.angle_sensitivity = worst
        self.angle_sensitivity_detail = detail

    def assumption_report(self) -> str:
        """B-1.  One line stating how far the theta = 0 assumption is from the
        truth for THIS fluid pair, so no result can be quoted without it."""
        if not np.isfinite(self.angle_sensitivity):
            return "angle sensitivity (B-1): not measured (table was loaded, " \
                   "not built)"
        det = sorted(self.angle_sensitivity_detail.items())
        d = "  ".join(f"|u|={u:g}: {v:.1%}" for u, v in det)
        where = max(det, key=lambda kv: kv[1])[0] if det else float("nan")
        if self.angle_sensitivity < 0.05:
            verdict = "benign"
        elif self.angle_sensitivity < 0.5:
            verdict = (f"SIGNIFICANT, worst at |u_bar| = {where:g} -- valid only "
                       f"where |u_bar| stays well below that; check the run's "
                       f"actual range")
        else:
            verdict = "LARGE -- the theta = 0 table is not valid for this pair"
        return (f"angle sensitivity (B-1): worst "
                f"{self.angle_sensitivity:.1%} [{verdict}]\n  {d}")

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
                                       max_iter=self.max_iter,
                                       r=solver.r, rho=solver.rho)
        sol = solver.solve_fixed_mean_velocity(c, [0.0, umag], [0.0, gb * H])
        if not sol.converged:
            self.n_failed += 1
        return closures_from_solution(sol)

    # -- persistence ----------------------------------------------------------
    def save(self, path):
        """Cache a built table.  Building is the expensive part of a
        Herschel-Bulkley run -- thousands of augmented-Lagrangian solves -- and
        the result depends only on the two fluids and the axes, so it is worth
        keeping between runs."""
        if not self._values:
            raise RuntimeError("build() first")
        np.savez_compressed(
            path, tuned_r=self.tuned_r, n_failed=self.n_failed,
            angle_sensitivity=self.angle_sensitivity,
            angle_u=np.array(sorted(self.angle_sensitivity_detail)),
            angle_dev=np.array([self.angle_sensitivity_detail[k] for k in
                                sorted(self.angle_sensitivity_detail)]),
            **{f"axis_{n}": a for n, a in self._axes},
            **{f"val_{k}": v for k, v in self._values.items()})
        return path

    def load(self, path):
        """Restore a cached table.  The axes are checked against this object's
        own, so a stale cache is a loud failure rather than a silent wrong
        answer."""
        z = np.load(path)
        axes = self._axis_list()
        for name, axis in axes:
            cached = z[f"axis_{name}"]
            if cached.shape != axis.shape or not np.allclose(cached, axis):
                raise ValueError(
                    f"cached table at {path} has a different '{name}' axis; "
                    f"delete it or change the cache key")
        self._axes = axes
        self._values = {k: z[f"val_{k}"] for k in ("I1", "I2", "q0", "I3")}
        self.tuned_r = float(z["tuned_r"])
        self.n_failed = int(z["n_failed"])
        if "angle_sensitivity" in z.files:
            self.angle_sensitivity = float(z["angle_sensitivity"])
            self.angle_sensitivity_detail = dict(zip(
                [float(x) for x in z["angle_u"]],
                [float(x) for x in z["angle_dev"]]))
        self._dinterp = None
        pts = tuple(a for _, a in self._axes)
        self._interp = {
            k: RegularGridInterpolator(pts, v, method="linear",
                                       bounds_error=False, fill_value=None)
            for k, v in self._values.items()
        }
        return self

    # -- c-derivative bounds --------------------------------------------------
    def _build_derivative_bounds(self):
        """
        Interpolators for an UPPER BOUND on |dq0/dc| and |dI3/dc|.

        The table is linear in c, so the interpolant's derivative is the
        segment slope -- piecewise constant and discontinuous at the nodes.
        At node k we store max(|slope_{k-1}|, |slope_k|), the larger of the two
        adjacent slopes.  Linearly interpolating THAT still bounds the true
        slope everywhere: inside segment k both endpoint values are >= |slope_k|
        by construction, so any convex combination of them is too.  A bound is
        exactly what the LLF monotonicity argument needs (NUM-26), and this one
        costs a single interpolation instead of the two four-quantity calls a
        finite difference of q0_I3 would take -- which measured at 244 ms per
        transport step against 48 ms for the whole elliptic Picard sequence.
        """
        c_axis = self._axes[0][1]
        if c_axis.size < 2:
            self._dinterp = None
            return
        dc = np.diff(c_axis).reshape((-1,) + (1,) * 3)
        out = {}
        for key in ("q0", "I3"):
            v = self._values[key]
            slope = np.abs(np.diff(v, axis=0) / dc)          # (n_c-1, ...)
            bound = np.empty_like(v)
            bound[0] = slope[0]
            bound[-1] = slope[-1]
            bound[1:-1] = np.maximum(slope[:-1], slope[1:])
            out[key] = bound
        pts = tuple(a for _, a in self._axes)
        self._dinterp = {
            k: RegularGridInterpolator(pts, v, method="linear",
                                       bounds_error=False, fill_value=None)
            for k, v in out.items()
        }

    def derivative_bounds(self, c, H=1.0, umag=1.0, gb=0.0):
        """Upper bounds on |dq0/dc| and |dI3/dc|; see `_build_derivative_bounds`."""
        if getattr(self, "_dinterp", None) is None:
            self._build_derivative_bounds()
        c, H, umag, gb = np.broadcast_arrays(*[np.asarray(x, float)
                                               for x in (c, H, umag, gb)])
        pts = np.stack([c.ravel(), H.ravel(), umag.ravel(), gb.ravel()], axis=-1)
        # B-2: this path had no range check at all, so out-of-range queries for
        # the LLF wavespeeds were silent BY CONSTRUCTION, not merely by a
        # command-line flag.  It extrapolates from the same axes as __call__ and
        # is called on every face of every step, so it is recorded identically.
        self._record_range(pts, warn=False)
        return (np.abs(self._dinterp["q0"](pts)).reshape(c.shape),
                np.abs(self._dinterp["I3"](pts)).reshape(c.shape))

    # -- evaluate -------------------------------------------------------------
    def __call__(self, c, H=1.0, umag=1.0, gb=0.0, warn_out_of_range=True):
        """Interpolate all four closures.  Broadcasting over array inputs."""
        if not self._interp:
            raise RuntimeError("call build() first")
        c, H, umag, gb = np.broadcast_arrays(*[np.asarray(x, float)
                                               for x in (c, H, umag, gb)])
        pts = np.stack([c.ravel(), H.ravel(), umag.ravel(), gb.ravel()], axis=-1)
        # The record is unconditional; `warn_out_of_range` now governs only
        # whether a warning is ALSO emitted (A-3).
        self._record_range(pts, warn=warn_out_of_range)
        vals = {k: f(pts).reshape(c.shape) for k, f in self._interp.items()}
        return vals["I1"], vals["I2"], vals["q0"], vals["I3"]

    def _record_range(self, pts, warn=True):
        """
        M3-T3, A-3, B-2.  Record -- always -- every query that falls outside the
        tabulated box, and optionally warn as well.

        The counters are the durable evidence.  A warning can be silenced by a
        `-W` flag, a `warnings.simplefilter`, or a pytest configuration, and if
        that happens the fact that a result rests on linear EXTRAPOLATION of the
        closures leaves no trace at all.  These do not go away.
        """
        self.n_range_queries += int(pts.shape[0])
        for j, (name, axis) in enumerate(self._axes):
            if axis.size == 1:
                continue
            lo, hi = float(axis[0]), float(axis[-1])
            below, above = lo - pts[:, j], pts[:, j] - hi
            worst = float(np.max(np.maximum(below, above)))
            if worst <= 1e-12:
                continue
            bad = (below > 1e-12) | (above > 1e-12)
            n_bad = int(np.count_nonzero(bad))
            self.n_range_points_out[name] = \
                self.n_range_points_out.get(name, 0) + n_bad
            k = int(np.argmax(np.maximum(below, above)))
            prev = self.worst_excursion.get(name)
            if prev is None or worst > prev[0]:
                self.worst_excursion[name] = (worst, float(pts[k, j]), lo, hi)
            if warn:
                warnings.warn(
                    f"closure table: axis '{name}' out of range at flat cell "
                    f"{k}: requested {pts[k, j]:.6g}, table covers "
                    f"[{lo:.6g}, {hi:.6g}]. Linear extrapolation is being used; "
                    f"widen the table or check the solution.",
                    ClosureTableRangeWarning, stacklevel=3)

    def reset_range_record(self):
        """Forget the accumulated out-of-range record (e.g. between runs)."""
        self.n_range_queries = 0
        self.n_range_points_out = {}
        self.worst_excursion = {}

    def range_report(self) -> str:
        """One line per offending axis; the empty string when the table was
        never queried outside its own box.  Printed by the run scripts so that
        no result can be quoted without it."""
        if not self.n_range_points_out:
            return ""
        out = [f"closure table queried OUTSIDE its range "
               f"({self.n_range_queries} point-queries total):"]
        for name, n in sorted(self.n_range_points_out.items()):
            worst, val, lo, hi = self.worst_excursion[name]
            out.append(f"  axis '{name}': {n} points out, worst request "
                       f"{val:.6g} against [{lo:.6g}, {hi:.6g}] "
                       f"(excursion {worst:.4g}); LINEAR EXTRAPOLATION was used")
        return "\n".join(out)

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
