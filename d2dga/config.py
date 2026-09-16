"""
Single home for every physical and numerical constant.

Convention: BF25 Section 2 throughout.  Nothing in this file is expressed in
another paper's scaling.  Dimensional quantities carry a `_m`, `_pa`, `_s`, ...
suffix; anything without one is already dimensionless.

Constants whose value is NOT pinned down by a source are marked
``# ASSUMPTION <ID>`` and carry a matching row in docs/assumptions.md.
"""

from dataclasses import dataclass, field

IN_TO_M = 0.0254


# --------------------------------------------------------------------------
# K-GEP-1 well
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class WellConfig:
    """Geothermal well K-GEP-1."""

    # ---- casing -----------------------------------------------------------
    # ASSUMPTION GEO-01.  User-selected 2026-09-16 from three candidates
    # (127 mm / 139.7 mm / 177.8 mm).  NOTE the tension recorded in
    # docs/assumptions.md: 177.8 mm yields delta/pi = 0.062 in gauge hole,
    # whereas the build spec states an expectation of ~0.12.  The spec's own
    # two delta/pi figures are mutually inconsistent; see assumptions.md.
    casing_od_m: float = 7.0 * IN_TO_M            # 0.17780 m

    # ---- hole -------------------------------------------------------------
    gauge_hole_diameter_m: float = 10.4 * IN_TO_M  # 0.26416 m
    washout_max_diameter_m: float = 23.0 * IN_TO_M

    # ---- depths (m below surface) ----------------------------------------
    casing_shoe_m: float = 194.0
    total_depth_m: float = 390.0
    washout_top_m: float = 195.0
    washout_bottom_m: float = 217.0

    # ---- inclination ------------------------------------------------------
    # Vertical well.  Kept as a field, NOT hard-coded: the buoyancy vector f
    # depends on beta, and docs/derivation.md Finding 1 shows that div_a.f
    # vanishes identically only because beta == 0.  An inclined caliper well
    # would need the 1/Fr*^2 term retained.
    inclination_rad: float = 0.0

    @property
    def open_hole_length_m(self) -> float:
        return self.total_depth_m - self.casing_shoe_m   # 196.0 m

    @property
    def casing_outer_radius_m(self) -> float:
        return 0.5 * self.casing_od_m

    def depth_of(self, xi_hat_m):
        """Map axial coordinate (B02 convention: measured UP from bottom hole)
        to true vertical depth below surface."""
        return self.total_depth_m - xi_hat_m


# --------------------------------------------------------------------------
# Eccentricity / standoff
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class StandoffConfig:
    """
    Casing offset model.

    e(xi) = e_hat(xi) / (2 d_hat(xi)),  B02 (16).

    Centralizer positions for K-GEP-1 are not available, so the default is a
    CONSTANT physical offset e_hat.  That offset is specified indirectly, by
    naming the eccentricity the casing would have in GAUGE hole; the code then
    back-computes e_hat = 2 * e_gauge * d_hat_gauge.  Because d_hat grows inside
    a washout while e_hat does not, e automatically FALLS in enlargements, which
    is the behaviour the build spec (Section 4) asks for.
    """

    # ASSUMPTION GEO-02 -- PENDING USER CONFIRMATION.
    # 0.3 corresponds to ~70% standoff, a routine field value, but no K-GEP-1
    # centralizer record has been supplied.  Sensitivity not yet tested.
    eccentricity_in_gauge_hole: float = 0.3


# --------------------------------------------------------------------------
# Synthetic wall (build spec Section 3.5) -- used until the caliper log lands
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class SyntheticWallConfig:
    amplitude_m: float = 1.5 * IN_TO_M
    wavelength_m: float = 20.0

    # ASSUMPTION GEO-03.  The spec writes the profile as `gauge + A sin(2 pi
    # xi/L)` for the RADIUS, but quotes both "gauge hole 10.4 in" and
    # "amplitude 1.5 in" -- and 10.4 in is unambiguously a DIAMETER.  Reading A
    # as a diameter amplitude keeps the two figures on the same footing and is
    # the default here.  Set amplitude_is_diameter=False for the literal
    # radius reading.  See assumptions.md for the feasibility consequences:
    # with 7 in casing the radius reading leaves only a 0.2 in minimum gap.
    amplitude_is_diameter: bool = True

    # 'symmetric'   : gauge + A sin(...)          -- literal spec form
    # 'enlargement' : gauge + A (1 - cos(...))/2  -- one-sided, >= gauge always
    # Real boreholes wash out but do not wash in, and the spec's A-sweep to
    # 6.0 in is geometrically impossible in the symmetric form.  See
    # assumptions.md GEO-04.
    mode: str = "symmetric"


# --------------------------------------------------------------------------
# Fluids
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class FluidsConfig:
    """
    Mud and cement properties.

    SOURCE, user-selected 2026-09-16: Table 1 of "Numerical Analysis of Cement
    Placement Into Drilling Fluid in Oilwell Applications", Materials 2025, 18,
    3098.  Taken verbatim; only the properties are taken, not that paper's
    geometry or flow direction.

    Three things about this source have to travel with any result built on it.
    See docs/assumptions.md FLU-01..03.

      1. Its "drilling fluid" is WATER -- 998 kg/m3, 1 mPa s, Newtonian.  A real
         drilling mud is 1100-1600 kg/m3 and carries a yield stress.  The
         consequence is m ~ 0.003-0.011 on K-GEP-1: the displacing fluid comes
         out 100-360 times more viscous than the displaced one, which is an
         unusually FAVOURABLE viscosity ratio.  Any displacement efficiency
         computed from this pair is optimistic.
      2. Its cement is denser than its "mud" (1200 vs 998), so b > 0 -- a
         favourable, stabilising density difference -- and b lands at 26-29 on
         K-GEP-1, comparable to ZF22's strongly buoyant cases 2/5/9.
      3. Its flow is DOWNWARD ("the cement slurry inlet velocity, which is
         downward").  Ours, and B02/PF04/ZF22/BF25 without exception, is upward
         annular displacement.  Only the properties cross over.

    The cement is Herschel-Bulkley with a large Bingham number (B = 6-32 over
    the paper's velocity range), so runs with this pair go through the M3
    tabulated-closure path and the Picard elliptic, NOT the fast Newtonian one
    -- and NUM-13's mobility-floor regularisation becomes load-bearing.
    """

    # Table 1, displaced fluid ("Drilling Fluid", Newtonian)
    mud_density: float = 998.0            # kg/m^3
    mud_viscosity: float = 1.0e-3         # Pa s

    # Table 1, displacing fluid ("Cement Slurry", Herschel-Bulkley)
    cement_density: float = 1200.0        # kg/m^3
    cement_consistency: float = 0.6       # Pa s^n
    cement_power_law_index: float = 0.4   # -
    cement_yield_stress: float = 1.4      # Pa

    # The paper's three inlet velocities.  Used here as the ANNULAR mean
    # velocity w0_hat, which is what BF25's scaling takes; the paper quotes
    # them as inlet velocities in its own much smaller annulus, so this is a
    # reading, not a measurement.  ASSUMPTION FLU-04.
    mean_velocities_m_s: tuple = (0.05, 0.2, 0.5)

    def as_fluids(self):
        """(displaced, displacing) in BF25's index convention."""
        from .scaling import HerschelBulkleyFluid
        return (HerschelBulkleyFluid("drilling fluid (water)", self.mud_density,
                                     self.mud_viscosity, 1.0, 0.0),
                HerschelBulkleyFluid("cement slurry", self.cement_density,
                                     self.cement_consistency,
                                     self.cement_power_law_index,
                                     self.cement_yield_stress))


# --------------------------------------------------------------------------
# Numerical grid
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class GridConfig:
    """
    Staggered mesh, BCF25 Fig. 2: Psi at cell corners, concentration at cell
    centres, velocities at face midpoints.

    n_phi / n_xi count CELLS.  BCF25 used 400 axial x 20 azimuthal over 200 m.
    """
    n_phi: int = 20
    n_xi: int = 400


@dataclass(frozen=True)
class Config:
    well: WellConfig = field(default_factory=WellConfig)
    standoff: StandoffConfig = field(default_factory=StandoffConfig)
    synthetic_wall: SyntheticWallConfig = field(default_factory=SyntheticWallConfig)
    fluids: FluidsConfig = field(default_factory=FluidsConfig)
    grid: GridConfig = field(default_factory=GridConfig)
