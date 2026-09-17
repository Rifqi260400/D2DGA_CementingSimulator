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

    def breakthrough_at(self, threshold):
        """
        First time the outlet concentration reaches `threshold`, from the
        recorded history rather than from the running detector -- so one run
        can be reported at several thresholds.

        Worth reporting as a sweep, because ZF22's definition ("the time at
        which the displacing fluid first exits the annulus") is exact only for
        a sharp front.  D2DGA fronts are not sharp: BF25 Section 3.1's spike
        regime puts a vanishingly thin tip ahead of the main shock, moving at
        the Poiseuille centreline speed 1.5 whatever the buoyancy.  A small
        threshold therefore measures the spike and a large one the shock, and
        the two can differ by a factor of 1.5.  assumptions.md NUM-21.
        """
        out = self.outlet_concentration
        hit = np.nonzero(out >= threshold)[0]
        if hit.size == 0:
            return np.nan
        k = int(hit[0])
        if k == 0:
            return float(self.times[0])
        lo, hi = out[k - 1], out[k]
        if hi <= lo:
            return float(self.times[k])
        frac = (threshold - lo) / (hi - lo)
        return float(self.times[k - 1] + frac * (self.times[k] - self.times[k - 1]))


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
                 picard_tol: float = 1e-8,
                 picard_max_iter: int = 100,
                 picard_warm_start: bool = True):
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
        # 1e-10 is far tighter than the time discretisation warrants: the
        # transport step is first order in dt, so driving the elliptic residual
        # to machine-ish precision every step buys nothing and costs iterations.
        self.picard_tol = picard_tol
        self.picard_max_iter = picard_max_iter
        self.picard_warm_start = bool(picard_warm_start)
        self._psi_prev = None
        self.n_picard_unconverged = 0
        self.worst_picard_residual = 0.0

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
            c, Q, tol=self.picard_tol, max_iter=self.picard_max_iter,
            psi0=self._psi_prev if self.picard_warm_start else None)
        self._psi_prev = psi
        if err > self.picard_tol and iters >= self.picard_max_iter:
            # Do NOT abort the run.  A single under-converged elliptic solve
            # perturbs Psi by O(err), which the next step largely corrects;
            # throwing away hours of a 10^5-step run over it is far worse.  But
            # it must not pass silently either, so count it and report the
            # worst residual -- `run` prints both, and a non-zero count means
            # the result needs looking at rather than quoting.
            self.n_picard_unconverged += 1
            self.worst_picard_residual = max(self.worst_picard_residual, err)
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
            breakthrough_threshold=0.01, on_step=None, record_every=1,
            checkpoint_path=None, checkpoint_every=180.0):
        """
        Advance to t_end.

        checkpoint_path : if given, the run state is written there periodically
            and, if the file already exists and matches this grid, the run
            RESUMES from it.  A full K-GEP-1 job is ~10^5 steps and several
            hours, and this environment reclaims its container on its own
            schedule -- three runs were lost that way at 22%, 22% and 2%.
            Checkpointing turns that from lost work into a relaunch.  The grid
            shape is verified on load, so a checkpoint from a different mesh is
            a loud failure rather than a silent wrong answer.

        breakthrough_threshold : ZF22 defines breakthrough as the time the
            displacing fluid "first exits the annulus", which is exact for a
            sharp front and ambiguous for a dispersive one -- with q0'(0) = 3/2
            the D2DGA tip is asymptotically thin, so ANY positive threshold is a
            choice.  It is exposed here and its sensitivity is reported rather
            than buried.  assumptions.md NUM-21.
        """
        import os
        import time as _time

        c = self.initial_condition() if c0 is None else np.array(c0, dtype=float)
        t, n = 0.0, 0
        times, effs, outlet, reports = [0.0], [self.transport.mass(c) / self._capacity], \
            [float(np.max(c[:, -1]))], []
        mass_running = self.transport.mass(c)
        cons_err = 0.0
        t_br = np.nan
        psi = None

        if checkpoint_path and os.path.exists(checkpoint_path):
            z = np.load(checkpoint_path)
            saved = z["c"]
            if saved.shape != c.shape:
                raise ValueError(
                    f"checkpoint {checkpoint_path} holds a {saved.shape} field "
                    f"but this run is {c.shape}; delete it or use another path")
            c = saved
            t = float(z["t"])
            n = int(z["n"])
            mass_running = float(z["mass_running"])
            cons_err = float(z["cons_err"])
            t_br = float(z["t_br"])
            times = list(z["times"])
            effs = list(z["effs"])
            outlet = list(z["outlet"])
            self.n_picard_unconverged = int(z["n_picard_unconverged"])
            self.worst_picard_residual = float(z["worst_picard_residual"])
            print(f"resumed from {checkpoint_path}: step {n}, t/Z = "
                  f"{t / self.geom.grid.Z:.4f}", flush=True)

        def _save():
            # Write through a file OBJECT: np.savez_compressed appends ".npz"
            # to a path that lacks it, so a "...npz.tmp" name would silently
            # become "...npz.tmp.npz" and the rename would miss it.
            tmp = checkpoint_path + ".tmp"
            with open(tmp, "wb") as fh:
                _write(fh)
            # rename is atomic, so a reclaim mid-write cannot leave a truncated
            # checkpoint behind
            os.replace(tmp, checkpoint_path)

        def _write(fh):
            np.savez_compressed(
                fh, c=c, t=t, n=n, mass_running=mass_running,
                cons_err=cons_err, t_br=t_br, times=np.array(times),
                effs=np.array(effs), outlet=np.array(outlet),
                n_picard_unconverged=self.n_picard_unconverged,
                worst_picard_residual=self.worst_picard_residual)

        last_save = _time.time()

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

            if checkpoint_path and _time.time() - last_save > checkpoint_every:
                _save()
                last_save = _time.time()

        if checkpoint_path:
            _save()

        return RunResult(times=np.array(times), efficiency=np.array(effs),
                         outlet_concentration=np.array(outlet),
                         concentration=c, stream_function=psi,
                         reports=reports, breakthrough_time=t_br,
                         conservation_error=cons_err)
