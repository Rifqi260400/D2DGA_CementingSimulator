"""
BCF25 Section IV -- the conservation verification, implemented as printed.

At this stage of the project this IS the verification: there is no CFD
comparison, and the model is checked against (a) the published ZF22 cases and
(b) this volume ledger, which is the test BCF25 uses for exactly the same
scheme on exactly the same equations.

BCF25's construction, quoted from the paper
-------------------------------------------
    "for each fluid, we compute the total volume in the annulus in two ways.
     First, we compute Vol1(k), by integrating c^n_{k,i+1/2,j+1/2}
     numerically over the interior volume of the annulus, i.e. we multiply by
     H r_a dphi dxi and sum over the interior cells.  Second, we compute
     Vol2(k) at each timestep, starting with the initial volumes in the
     annulus and adding/subtracting the inflow/outflow of each fluid at the
     start and finish of the annulus, at each timestep.  We normalize these
     quantities with the total annulus volume and subtracted to give a
     relative volumetric error."

So, per fluid k and per timestep,

    err_k = [ Vol1(k) - Vol2(k) ] / V_total.

Their figures 7 and 15 plot err_k against time, one curve per fluid, and their
finding is that it stays at O(1e-15) -- "comparable to the standard Matlab
precision" -- across five cases, two models and both Newtonian and
non-Newtonian fluids.

Two things this implementation has to be honest about
-----------------------------------------------------
**1. The normalisation was different here before.**  `Simulation.run` divided
by `Vol2` rather than by the total annulus volume.  Early in a run `Vol2 -> 0`,
so that denominator inflates the error exactly where the absolute error is
smallest and least interesting; late in a run the two agree, because
`Vol2 -> V_total`.  Neither is wrong, but only one is BCF25's, and comparing a
number against their 1e-15 requires computing the same number.  Both are
reported.

**2. With K = 2, "the concentrations sum to 1" is structural, not a check.**
BCF25 evolve K = 3 concentrations independently, so their per-fluid curves are
three genuinely separate ledgers and their sum-to-one property is a real
result.  Here only c_2 is evolved and c_1 = 1 - c_2 pointwise by definition.
The two-fluid closure then makes fluid 1's face flux *identically* the total
volumetric face flux minus fluid 2's -- including the LLF dissipation term,
which changes sign with c and cancels exactly.  So err_1 is not independent
evidence; it differs from -err_2 only by

    (total volume in) - (total volume out),

which is the elliptic solver's business and IS an independent check.  That
quantity is reported separately rather than being hidden inside a second curve
that looks like corroboration.  `total_flux_imbalance` is the number to read.

The total volumetric flux through an axial face is exact and closure-free:
summing the advective flux with q0 = 1 over phi telescopes the stream
function, giving 0.5 * (Psi[-1, j] - Psi[0, j]).  Verified against the
face flux with the annulus full of fluid 2, where the two must coincide: they
agree to 0.0 (test M8-T3).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def total_axial_flux(psi, j):
    """Total volumetric flux through axial face `j`, from the stream function.

    Closure-free: the advective flux of ALL fluids through a face is
    0.25 (q0_S + q0_N) dPsi with q0 = 1, and summing over phi telescopes.  It
    therefore measures what the elliptic solve delivered, independently of
    anything the transport scheme did with it.
    """
    psi = np.asarray(psi, dtype=float)
    return 0.5 * float(psi[-1, j] - psi[0, j])


@dataclass
class ConservationLedger:
    """BCF25 Section IV, accumulated over a run.

    `capacity` is the total annulus volume, `mass(ones)`.  Feed it one
    `record` per timestep; it keeps the two volumes for each fluid and the
    relative errors, so the series can be plotted as their Figs. 7 and 15.
    """

    capacity: float
    vol1_2_initial: float = 0.0

    times: list = field(default_factory=list)
    err1: list = field(default_factory=list)      # fluid 1 (displaced)
    err2: list = field(default_factory=list)      # fluid 2 (displacing)
    imbalance: list = field(default_factory=list)  # (in - out) of TOTAL volume

    _vol2_2: float = 0.0        # ledger volume of fluid 2
    _vol2_1: float = 0.0        # ledger volume of fluid 1
    _net_total: float = 0.0     # cumulative (total in - total out)
    # True when the run resumed from a checkpoint written before this ledger
    # existed: the series then covers only the steps taken since, and any
    # figure or verdict drawn from it describes part of a run.
    resumed_without_history: bool = False

    def __post_init__(self):
        self._vol2_2 = float(self.vol1_2_initial)
        self._vol2_1 = float(self.capacity) - float(self.vol1_2_initial)

    def start(self, t, vol1_2):
        """Record t = 0, where both ledgers agree by construction."""
        self.times.append(float(t))
        self.err1.append(0.0)
        self.err2.append(0.0)
        self.imbalance.append(0.0)

    def record(self, t, vol1_2, flux_2, flux_total):
        """One timestep.

        vol1_2     : integral of c_2 over the interior, AFTER the step
        flux_2     : net volume of fluid 2 that crossed the boundaries
        flux_total : net volume of ALL fluid that crossed them
        """
        self._vol2_2 += float(flux_2)
        self._vol2_1 += float(flux_total) - float(flux_2)
        self._net_total += float(flux_total)
        v = float(self.capacity)
        vol1_1 = v - float(vol1_2)
        self.times.append(float(t))
        self.err2.append((float(vol1_2) - self._vol2_2) / v)
        self.err1.append((vol1_1 - self._vol2_1) / v)
        self.imbalance.append(self._net_total / v)

    # -- summary ---------------------------------------------------------
    @property
    def worst(self) -> float:
        """Largest |err| over the run, over both fluids -- the single number
        to compare against BCF25's ~1e-15."""
        if not self.err2:
            return 0.0
        return max(max(abs(e) for e in self.err1),
                   max(abs(e) for e in self.err2))

    @property
    def total_flux_imbalance(self) -> float:
        """Largest |(total in) - (total out)| / V over the run.

        The only part of fluid 1's ledger that is NOT algebraically implied by
        fluid 2's, and therefore the independent content of the two-fluid
        version of this check: it asks whether the elliptic solve delivered
        the same volumetric flux at both ends.
        """
        return max((abs(x) for x in self.imbalance), default=0.0)

    def report(self) -> str:
        if not self.err2:
            return "conservation ledger (BCF25 IV): no steps recorded"
        n = len(self.err2) - 1
        partial = (" PARTIAL -- resumed from a checkpoint written before this "
                   "ledger existed, so it covers only the steps taken since."
                   if self.resumed_without_history else "")
        return (f"conservation ledger (BCF25 IV, normalised by the annulus "
                f"volume): worst |err| = {self.worst:.2e} over {n} steps "
                f"[fluid 1 {max(abs(e) for e in self.err1):.2e}, "
                f"fluid 2 {max(abs(e) for e in self.err2):.2e}]; "
                f"total-flux imbalance {self.total_flux_imbalance:.2e}.  "
                f"BCF25 report ~1e-15 over their five cases." + partial)

    def as_dict(self) -> dict:
        return {"times": list(self.times), "err1": list(self.err1),
                "err2": list(self.err2), "imbalance": list(self.imbalance),
                "worst": self.worst,
                "total_flux_imbalance": self.total_flux_imbalance,
                "capacity": float(self.capacity),
                "partial": bool(self.resumed_without_history)}
