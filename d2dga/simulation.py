"""
M6 -- the coupled time loop.

Each step does exactly two things, in this order:

  1. solve the elliptic problem for Psi^n from c^n and the imposed flow rate
     Q(t^n)  (M4);
  2. advance c^n -> c^{n+1} explicitly with that frozen Psi^n  (M5).

This is not operator splitting in the usual sense.  BF25 (2.16) is a
CONSTRAINT, not an evolution equation: there is no time derivative in it, and
Psi is determined instantaneously by c and Q.  So the only time discretisation
in the whole model is the explicit Euler step in the transport equation, and
the scheme is first-order in time for that reason alone.  The consequence that
matters is that the LLF monotonicity bound BCF25 (44) is evaluated on Psi^n and
c^n -- both already known when the step is taken -- so the discrete maximum
principle carries over to the coupled loop unchanged.  That is asserted
directly rather than assumed (M6-T2).

What DOES enter here and not in M4 or M5 is the lag between Psi and c: Psi^n is
used to move c across the whole of [t^n, t^{n+1}], while the true Psi varies
over that interval.  That error is O(dt) and is the reason the CFL number is
not simply pushed to 1 even though M5's evidence favours large steps
(assumptions.md NUM-01).  M6-T4 measures it.

Fluid numbering follows BF25 Section 2.3 throughout -- fluid 1 displaced (mud),
fluid 2 displacing (cement), c_bar = c_bar_2.  NOT PF04's, which is reversed;
see gapscale/twodga.py and assumptions.md CONV-09.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .elliptic import ClosureProvider, StreamFunctionSolver
from .geometry import Geometry
from .transport import TransportSolver


@dataclass
class StepReport:
    """One timestep's diagnostics."""
    n: int
    t: float
    dt: float
    Q: float
    mass: float
    efficiency: float
    boundary_flux: float
    c_min: float
    c_max: float
    picard_iterations: int
    static_cells: int


@dataclass
class RunResult:
    times: np.ndarray
    efficiency: np.ndarray
    outlet_concentration: np.ndarray        # max over phi at xi = Z, per step
    concentration: np.ndarray               # final field
    stream_function: np.ndarray             # final Psi
    reports: list = field(default_factory=list)
    breakthrough_time: float = np.nan
    conservation_error: float = 0.0

    def efficiency_at(self, t):
        """Linearly interpolated displacement efficiency at time t."""
        return float(np.interp(t, self.times, self.efficiency))


class Simulation:
    """
    Parameters
    ----------
    geometry, closures : as for the M4/M5 solvers
    froude, delta_rho  : BF25 (2.4) Fr* and Delta_rho.  The buoyancy number
                         b = -Delta_rho / Fr*^2 (CONV-04) is DERIVED from them
                         rather than passed separately, so the elliptic source
                         and the transport flux cannot disagree about it -- a
                         mismatch there would be silent and would look like a
                         physical result.
    flow_rate          : float, or a callable t -> Q(t) for a pump schedule
    """

    def __init__(self, geometry: Geometry, closures: ClosureProvider,
                 froude: float, delta_rho: float,
                 flow_rate=1.0,
                 cfl: float = 0.5,
                 inflow_concentration: float = 1.0,
                 inflow: str = "no_axial_gradient",
                 mobility_floor: float = 1e-12,
                 picard_tol: float = 1e-10,
                 picard_max_iter: int = 100):
        self.geom = geometry
        self.closures = closures
        self.buoyancy_number = -float(delta_rho) / float(froude) ** 2
        self.elliptic = StreamFunctionSolver(
            geometry, closures, froude=froude, delta_rho=delta_rho,
            inflow=inflow, mobility_floor=mobility_floor)
        self.transport = TransportSolver(
            geometry, closures, buoyancy_number=self.buoyancy_number,
            inflow_concentration=inflow_concentration, cfl=cfl)
        self.flow_rate = flow_rate
        self.picard_tol = picard_tol
        self.picard_max_iter = picard_max_iter

        # total pore volume of the (half) annulus, for the efficiency metric
        self._capacity = self.transport.mass(np.ones(
            (geometry.grid.n_phi, geometry.grid.n_xi)))

    # ------------------------------------------------------------------
    def Q(self, t):
        return float(self.flow_rate(t)) if callable(self.flow_rate) \
            else float(self.flow_rate)

    def initial_condition(self, value=0.0):
        """Uniform field; 0.0 = annulus full of the displaced fluid."""
        g = self.geom.grid
        return np.full((g.n_phi, g.n_xi), float(value))

    def solve_stream_function(self, c, Q):
        if self.closures.is_linear:
            return self.elliptic.solve(c, Q), 1
        psi, iters, err = self.elliptic.solve_nonlinear(
            c, Q, tol=self.picard_tol, max_iter=self.picard_max_iter)
        if err > self.picard_tol and iters >= self.picard_max_iter:
            raise RuntimeError(
                f"Picard iteration did not converge: {iters} iterations, "
                f"residual {err:.3e} > {self.picard_tol:.3e}")
        return psi, iters

    # ------------------------------------------------------------------
    def step(self, c, t, dt=None, max_dt=None):
        """One coupled step.  Returns (c_new, psi, dt, Q, boundary_flux, iters)."""
        Q = self.Q(t)
        psi, iters = self.solve_stream_function(c, Q)
        v, w = self.elliptic.velocities(psi)
        umag = np.hypot(v, w)
        c_new, dt, bflux = self.transport.step(c, psi, dt=dt, umag=umag,
                                               max_dt=max_dt)
        return c_new, psi, dt, Q, bflux, iters

    # ------------------------------------------------------------------
    def run(self, t_end, c0=None, max_steps=2_000_000,
            breakthrough_threshold=0.01, on_step=None, record_every=1):
        """
        Advance to t_end.

        breakthrough_threshold : ZF22 defines breakthrough as the time the
            displacing fluid "first exits the annulus", which is exact for a
            sharp front and ambiguous for a dispersive one -- with q0'(0) = 3/2
            the D2DGA tip is asymptotically thin, so ANY positive threshold is a
            choice.  It is exposed here and its sensitivity is reported rather
            than buried.  assumptions.md NUM-21.
        """
        c = self.initial_condition() if c0 is None else np.array(c0, dtype=float)
        t, n = 0.0, 0
        times, effs, outlet, reports = [0.0], [self.transport.mass(c) / self._capacity], \
            [float(np.max(c[:, -1]))], []
        mass_running = self.transport.mass(c)
        cons_err = 0.0
        t_br = np.nan
        psi = None

        prev_t, prev_out = 0.0, outlet[0]
        while t < t_end and n < max_steps:
            c_new, psi, dt, Q, bflux, iters = self.step(
                c, t, max_dt=t_end - t)
            mass_running += bflux
            prev_t, prev_out = t, float(np.max(c[:, -1]))
            c = c_new
            t += dt
            n += 1

            m = self.transport.mass(c)
            cons_err = max(cons_err,
                           abs(m - mass_running) / max(abs(mass_running), 1e-30))
            out = float(np.max(c[:, -1]))
            if np.isnan(t_br) and out >= breakthrough_threshold:
                # linear interpolation between the bracketing steps.  Uses the
                # previous STEP, not the previous recorded sample: with
                # record_every > 1 those differ, and taking the recorded one
                # would make t_br depend on the output frequency.
                if out > prev_out:
                    frac = (breakthrough_threshold - prev_out) / (out - prev_out)
                    t_br = prev_t + frac * (t - prev_t)
                else:
                    t_br = t

            if n % record_every == 0 or t >= t_end:
                times.append(t)
                effs.append(m / self._capacity)
                outlet.append(out)
                rep = StepReport(n=n, t=t, dt=dt, Q=Q, mass=m,
                                 efficiency=m / self._capacity,
                                 boundary_flux=bflux,
                                 c_min=float(c.min()), c_max=float(c.max()),
                                 picard_iterations=iters,
                                 static_cells=self.elliptic.n_static_cells)
                reports.append(rep)
                if on_step is not None:
                    on_step(rep, c, psi)

        return RunResult(times=np.array(times), efficiency=np.array(effs),
                         outlet_concentration=np.array(outlet),
                         concentration=c, stream_function=psi,
                         reports=reports, breakthrough_time=t_br,
                         conservation_error=cons_err)
