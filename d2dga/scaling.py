"""
M1 -- non-dimensionalisation.  BF25 convention ONLY.

Every scaling in the codebase is defined here and nowhere else.  Anything
imported from another paper is converted at the boundary, with the conversion
written down (see `zf22_buoyancy_to_bf25` and docs/assumptions.md CONV-02).

BF25 Section 2 relations implemented:

    gamma_dot_0 = w0_hat / d_star_hat                              (2.29)
    mu_e_hat    = [kappa1 gd0^(n1-1) * kappa2 gd0^(n2-1)]^(1/2)    (2.29)
    tau_0_hat   = mu_e_hat gamma_dot_0 + max(tauY1, tauY2)         (2.30)
    kappa_k     = kappa_k_hat gd0^n_k / tau_0_hat                  (2.32)
                = m^(+/-1/2) / (1 + B)
    tau_Y_k     = tauY_k_hat / tau_0_hat                           (2.32)
    m           = kappa1_hat gd0^(n1-n2) / kappa2_hat              (2.32)
    B           = max(tauY1, tauY2) / (mu_e_hat gamma_dot_0)       (2.33)
    Fr_star     = sqrt( tau_0_hat / (rho1_hat g d_star_hat) )      (2.4)

Buoyancy.  Per docs/derivation.md Finding 2 (CONV-04) the scalar that multiplies
f in the ELLIPTIC equation is r_a-free:

    b_num = -Delta_rho / Fr_star^2 = (rho2_hat - rho1_hat) g d_star_hat / tau_0_hat

BF25 (2.12) carries an extra 1/r_a which belongs to G_b, not here.  We store the
r_a-free number and multiply by r_a where the geometry demands it.

Index convention throughout: fluid 1 = displaced (in situ mud),
fluid 2 = displacing (cement slurry).  K = 2 (build spec Section 4.1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

G_ACCEL = 9.81  # m/s^2


# ===========================================================================
# Fluids
# ===========================================================================
@dataclass(frozen=True)
class HerschelBulkleyFluid:
    """Dimensional Herschel-Bulkley fluid: tau = tau_Y + kappa gammadot^n."""

    name: str
    density: float          # kg/m^3
    consistency: float      # Pa s^n
    power_law_index: float  # -
    yield_stress: float = 0.0   # Pa

    @property
    def is_newtonian(self) -> bool:
        return self.power_law_index == 1.0 and self.yield_stress == 0.0

    def apparent_viscosity_at(self, gamma_dot: float) -> float:
        """kappa gd^(n-1) + tau_Y/gd -- the effective viscosity at a shear rate."""
        visc = self.consistency * gamma_dot ** (self.power_law_index - 1.0)
        if self.yield_stress:
            visc += self.yield_stress / gamma_dot
        return visc

    @classmethod
    def newtonian(cls, name, density, viscosity):
        return cls(name, density, viscosity, 1.0, 0.0)


@dataclass(frozen=True)
class ScaledFluid:
    """Dimensionless counterpart. BF25 (2.31)-(2.32)."""
    name: str
    density: float          # scaled with rho1_hat
    consistency: float      # kappa_k
    power_law_index: float  # n_k, unchanged by scaling
    yield_stress: float     # tau_Y,k, in [0, 1)


# ===========================================================================
# Scaling
# ===========================================================================
class Scaling:
    """
    Holds every scale factor and every dimensionless group for one simulation.

    Construct from the dimensional problem; read dimensionless quantities off
    the properties.  `to_dimensional_*` / `to_dimensionless_*` are exact
    inverses (test M1-T1).
    """

    def __init__(self,
                 fluid_displaced: HerschelBulkleyFluid,
                 fluid_displacing: HerschelBulkleyFluid,
                 r_a_hat_star: float,
                 delta_star: float,
                 flow_rate: float | None = None,
                 mean_velocity: float | None = None,
                 g: float = G_ACCEL):
        self.fluid1 = fluid_displaced
        self.fluid2 = fluid_displacing
        self.r_a_hat_star = float(r_a_hat_star)
        self.delta_star = float(delta_star)
        self.g = float(g)

        # d_star_hat := delta_star * r_a_hat_star.
        # This is forced, not chosen: B02 (21) defines H so that H = 1
        # corresponds to a half-gap of delta_star * r_a_hat_star.  Using an
        # independently averaged d_hat would make H = 1 mean something else.
        # ASSUMPTION SCA-01.
        self.d_hat_star = self.delta_star * self.r_a_hat_star

        # Annulus reference area, B02: A_star = 4 pi delta_star r_a_star^2
        #                                    = 4 pi d_star_hat r_a_star
        self.area_hat = 4.0 * np.pi * self.d_hat_star * self.r_a_hat_star

        if (flow_rate is None) == (mean_velocity is None):
            raise ValueError("give exactly one of flow_rate, mean_velocity")
        if mean_velocity is None:
            self.flow_rate_hat = float(flow_rate)
            self.w0_hat = self.flow_rate_hat / self.area_hat
        else:
            self.w0_hat = float(mean_velocity)
            self.flow_rate_hat = self.w0_hat * self.area_hat

        # ---- BF25 (2.29): shear-rate and effective-viscosity scales --------
        self.gamma_dot_0 = self.w0_hat / self.d_hat_star
        mu1 = self.fluid1.consistency * self.gamma_dot_0 ** (self.fluid1.power_law_index - 1.0)
        mu2 = self.fluid2.consistency * self.gamma_dot_0 ** (self.fluid2.power_law_index - 1.0)
        self.mu_e_hat = float(np.sqrt(mu1 * mu2))

        # ---- BF25 (2.33), (2.30): Bingham number and stress scale ----------
        self.tau_y_max_hat = max(self.fluid1.yield_stress, self.fluid2.yield_stress)
        self.B = self.tau_y_max_hat / (self.mu_e_hat * self.gamma_dot_0)
        self.tau_0_hat = self.mu_e_hat * self.gamma_dot_0 + self.tau_y_max_hat

        # ---- BF25 (2.32): generalised viscosity ratio ----------------------
        self.m = (self.fluid1.consistency
                  * self.gamma_dot_0 ** (self.fluid1.power_law_index
                                         - self.fluid2.power_law_index)
                  / self.fluid2.consistency)

        # ---- derived scales -------------------------------------------------
        self.length_hat = np.pi * self.r_a_hat_star        # axial / azimuthal
        self.time_hat = self.length_hat / self.w0_hat
        self.velocity_hat = self.w0_hat
        self.pressure_hat = self.tau_0_hat / self.d_hat_star * self.length_hat

        # ---- BF25 (2.4): Froude number --------------------------------------
        self.Fr_star = float(np.sqrt(self.tau_0_hat
                                     / (self.fluid1.density * self.g * self.d_hat_star)))

    # -- dimensionless densities ------------------------------------------
    @property
    def rho1(self) -> float:
        return 1.0

    @property
    def rho2(self) -> float:
        return self.fluid2.density / self.fluid1.density

    @property
    def delta_rho(self) -> float:
        """Delta_rho = (rho1_hat - rho2_hat)/rho1_hat.  BF25 (2.7).
        Positive when the DISPLACED fluid is the heavier one."""
        return (self.fluid1.density - self.fluid2.density) / self.fluid1.density

    @property
    def buoyancy_number(self) -> float:
        """
        b_num = -Delta_rho / Fr*^2, the r_a-FREE scalar (CONV-04).

        Sign: positive for the normal cementing case, denser fluid displacing
        lighter upward.  BF25 signs b so that b > 0 is the favourable case.
        """
        return -self.delta_rho / self.Fr_star ** 2

    # -- scaled fluids, BF25 (2.32) ---------------------------------------
    def _scale_fluid(self, fluid: HerschelBulkleyFluid) -> ScaledFluid:
        return ScaledFluid(
            name=fluid.name,
            density=fluid.density / self.fluid1.density,
            consistency=(fluid.consistency
                         * self.gamma_dot_0 ** fluid.power_law_index / self.tau_0_hat),
            power_law_index=fluid.power_law_index,
            yield_stress=fluid.yield_stress / self.tau_0_hat,
        )

    @property
    def scaled_fluid1(self) -> ScaledFluid:
        return self._scale_fluid(self.fluid1)

    @property
    def scaled_fluid2(self) -> ScaledFluid:
        return self._scale_fluid(self.fluid2)

    # -- round trip, M1-T1 --------------------------------------------------
    def to_dimensionless_length(self, x):   return np.asarray(x) / self.length_hat
    def to_dimensional_length(self, x):     return np.asarray(x) * self.length_hat
    def to_dimensionless_radial(self, y):   return np.asarray(y) / self.d_hat_star
    def to_dimensional_radial(self, y):     return np.asarray(y) * self.d_hat_star
    def to_dimensionless_velocity(self, u): return np.asarray(u) / self.velocity_hat
    def to_dimensional_velocity(self, u):   return np.asarray(u) * self.velocity_hat
    def to_dimensionless_time(self, t):     return np.asarray(t) / self.time_hat
    def to_dimensional_time(self, t):       return np.asarray(t) * self.time_hat
    def to_dimensionless_stress(self, s):   return np.asarray(s) / self.tau_0_hat
    def to_dimensional_stress(self, s):     return np.asarray(s) * self.tau_0_hat
    def to_dimensionless_pressure(self, p): return np.asarray(p) / self.pressure_hat
    def to_dimensional_pressure(self, p):   return np.asarray(p) * self.pressure_hat
    def to_dimensionless_density(self, r):  return np.asarray(r) / self.fluid1.density
    def to_dimensional_density(self, r):    return np.asarray(r) * self.fluid1.density
    def to_dimensionless_flow_rate(self, q): return np.asarray(q) / self.flow_rate_hat
    def to_dimensional_flow_rate(self, q):   return np.asarray(q) * self.flow_rate_hat

    # -- reference numbers not used by the model, but reported --------------
    @property
    def reynolds_zf22(self) -> float:
        """Re = rho1 w0 d_hat / mu1, ZF22/ZF23 definition.  Diagnostic only:
        the D2DGA model is non-inertial, so Re never enters the solution.  It
        is carried because ZF22 Tables 1-3 are indexed by it."""
        mu1 = self.fluid1.consistency * self.gamma_dot_0 ** (
            self.fluid1.power_law_index - 1.0)
        return self.fluid1.density * self.w0_hat * self.d_hat_star / mu1

    def summary(self) -> str:
        return "\n".join([
            f"w0_hat        {self.w0_hat:.6g} m/s     Q_hat {self.flow_rate_hat:.6g} m^3/s",
            f"d_star_hat    {self.d_hat_star*1e3:.4f} mm   r_a_star {self.r_a_hat_star*1e3:.4f} mm",
            f"gamma_dot_0   {self.gamma_dot_0:.6g} 1/s",
            f"mu_e_hat      {self.mu_e_hat:.6g} Pa s",
            f"tau_0_hat     {self.tau_0_hat:.6g} Pa",
            f"B             {self.B:.6g}",
            f"m             {self.m:.6g}",
            f"Fr_star       {self.Fr_star:.6g}",
            f"Delta_rho     {self.delta_rho:.6g}",
            f"b (r_a-free)  {self.buoyancy_number:.6g}",
            f"time_hat      {self.time_hat:.6g} s",
        ])


# ===========================================================================
# Cross-paper conversions (build spec Section 1.1)
# ===========================================================================
def zf22_buoyancy_to_bf25(b_zf22: float, m: float) -> float:
    """
    ZF22 buoyancy number -> BF25 buoyancy number.   CONV-02.

        b_ZF22  = (rho2 - rho1) g d^2 / (mu1 w0)          ZF22 Section 2.3
        b_BF25  = (rho2 - rho1) g d^2 / ((mu1 mu2)^0.5 w0) BF25 (2.12) + (2.4)

    The two differ only in which viscosity scales the stress: ZF22 uses the
    displaced fluid, BF25 the geometric mean.  Hence

        b_BF25 = b_ZF22 * mu1 / (mu1 mu2)^0.5 = b_ZF22 * m^(1/2)

    Independent confirmation: BF25 Section 3.2 states its own translation to
    Lajeunesse's parameters as "b = 3 m^0.5 / U", and U = 3 / b_ZF22.  Composing
    gives b_BF25 = m^0.5 b_ZF22, the same factor.  Two routes, one answer.
    """
    return b_zf22 * np.sqrt(m)


def bf25_buoyancy_to_zf22(b_bf25: float, m: float) -> float:
    return b_bf25 / np.sqrt(m)
