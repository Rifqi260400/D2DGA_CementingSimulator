"""
M2 -- closure functionals from a gap-scale solution.

BF25 (2.14), (2.15), (2.22), (2.23), evaluated in the unit channel (H = 1) so
that what comes out is the H-free script-I family.  H is reinstated by the
caller via the BF25 (A4) rescaling, which is the whole point of build spec
Section 3.2: one table, every depth.

Every integrand is written in terms of the FLUIDITY 1/eta rather than eta, so
unyielded material contributes exactly zero instead of dividing by infinity.
That is not a numerical dodge -- it is what BF25 Section 2.3 describes when it
says the displacing fluid's contribution to the mean mobility "is therefore
zero since eta2 -> infinity and gammadot -> 0".

The five primitive integrals (H = 1, interface at c):

    A  = int_0^c  y^2 / eta2  dy          displacing-layer second moment
    C  = int_c^1  y^2 / eta1  dy          wall-layer second moment
    Bq = int_c^1  y   / eta1  dy
    D  = int_c^1  (1-y) / eta1  dy
    E  = int_c^1  y(1-y) / eta1  dy

from which

    I1 = A + C                                            (2.14)
    I2 = (1-c) A + c E                                    (2.15)
    q0 = (A + c Bq) / I1                                  (2.22)
    I3 = [ ((1-c) A + c^2 D) I1 - (A + c Bq) I2 ] / I1    (2.23)

The bracket ORDER in I3 is the one that reproduces BF25's own Lajeunesse
translation (Section 3.2); the opposite order gives the right magnitude with
the wrong sign.  See newtonian.py for the full CONV-03 write-up.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .al_solver import GapSolution


@dataclass(frozen=True)
class Closures:
    """The four functionals the 2-D model needs, at one (c, rheology, G) state."""
    c: float
    I1: float      # script_I1, mean mobility
    I2: float      # script_I2, buoyant mobility
    q0: float      # isotropic flux
    I3: float      # buoyant flux distribution

    def as_tuple(self):
        return self.I1, self.I2, self.q0, self.I3


def _cell_integral(integrand, dy):
    """Midpoint rule on the cell-centred grid.  Second order in 1/N_y, which is
    what sets the useful tolerance floor (test M2-T5)."""
    return float(np.sum(integrand * dy))


def closures_from_solution(sol: GapSolution) -> Closures:
    """Evaluate BF25 (2.14), (2.15), (2.22), (2.23) from a converged gap solve."""
    c = sol.c
    y = sol.y_centres
    dy = sol.dy
    fl = sol.fluidity          # 1/eta at cell centres
    f2 = sol.in_fluid2
    f1 = ~f2

    A = _cell_integral((y ** 2 * fl)[f2], dy[f2]) if np.any(f2) else 0.0
    C = _cell_integral((y ** 2 * fl)[f1], dy[f1]) if np.any(f1) else 0.0
    Bq = _cell_integral((y * fl)[f1], dy[f1]) if np.any(f1) else 0.0
    D = _cell_integral(((1.0 - y) * fl)[f1], dy[f1]) if np.any(f1) else 0.0
    E = _cell_integral((y * (1.0 - y) * fl)[f1], dy[f1]) if np.any(f1) else 0.0

    I1 = A + C
    if I1 <= 0.0:
        # Everything unyielded: no flow anywhere.  I1 = 0 makes S singular, so
        # the caller must treat this as a static cell rather than divide by it.
        return Closures(c, 0.0, 0.0, 0.0, 0.0)

    I2 = (1.0 - c) * A + c * E
    q0 = (A + c * Bq) / I1
    I3 = (((1.0 - c) * A + c ** 2 * D) * I1 - (A + c * Bq) * I2) / I1
    return Closures(c, I1, I2, q0, I3)


def closures_newtonian_analytic(c, m) -> Closures:
    """Convenience wrapper on the closed forms, same return type."""
    from . import newtonian as _n
    return Closures(float(c), float(_n.script_I1(c, m)), float(_n.script_I2(c, m)),
                    float(_n.q0(c, m)), float(_n.script_I3(c, m)))
