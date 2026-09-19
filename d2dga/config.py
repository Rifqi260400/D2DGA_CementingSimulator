"""
Single home for every physical and numerical constant.

Convention: BF25 Section 2 throughout.  Nothing in this file is expressed in
another paper's scaling.  Dimensional quantities carry a `_m`, `_pa`, `_s`, ...
suffix; anything without one is already dimensionless.

Constants whose value is NOT pinned down by a source are marked
``# ASSUMPTION <ID>`` and carry a matching row in docs/assumptions.md.
"""

from dataclasses import dataclass, field, fields

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
         computed from this pair is optimistic.  This was raised with the user
         and re-affirmed on 2026-09-17: the mud is to be treated as water-like.
         Results stand, with the caveat attached -- they describe a WATER
         displacement and bound from above what a real mud would give.
         Since 2026-09-19 the DEFAULT is still that water, but it is no longer
         the only expressible mud: `mud_consistency`, `mud_power_law_index` and
         `mud_yield_stress` make a real Herschel-Bulkley mud constructible.
         Read the B-1 warning on those fields before using one.
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

    # Table 1, displaced fluid ("Drilling Fluid", Newtonian by default).
    #
    # The three rheology fields are NOT a restatement of one viscosity: a
    # Herschel-Bulkley fluid has no single viscosity, so `mud_consistency`
    # carries kappa_hat in Pa s^n and coincides with the Newtonian viscosity
    # only at n = 1.  The field `mud_viscosity` therefore no longer exists;
    # `mud_viscosity` survives as a PROPERTY that returns kappa_hat for a
    # Newtonian mud and refuses for any other, so no caller can read a number
    # in the wrong units.  Defaults reproduce Table 1's water exactly
    # (1.0e-3, n = 1, tau_Y = 0), so every result computed before this field
    # existed is bit-identical -- verified by test M0-T5.
    #
    # WARNING, and it is the reason this was not added earlier: a mud with a
    # yield stress puts BOTH fluids in the yield-stress class, which is the
    # regime where finding B-1 is O(1) rather than a correction.  Every entry
    # of the closure table is built with u_bar parallel to G~_b (theta = 0);
    # for a pair with tau_Y on both sides the theta = 90 degree closures differ by
    # up to I2 +134%.  `ClosureTable.assumption_report()` measures this per
    # table and must be read before quoting any number from such a run.  See
    # docs/assumptions.md FLU-05 and AUDIT_REPORT.md B-1.
    mud_density: float = 998.0            # kg/m^3
    mud_consistency: float = 1.0e-3       # Pa s^n  (= viscosity when n = 1)
    mud_power_law_index: float = 1.0      # -
    mud_yield_stress: float = 0.0         # Pa

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

    # ---- validation ------------------------------------------------------
    # A dataclass with silent defaults is exactly how a UI ends up displaying a
    # plausible number for an input nobody supplied, which the UI build spec
    # forbids.  These checks make an incoherent fluid impossible to construct,
    # so the failure surfaces at the point of entry rather than as a strange
    # closure table 40 minutes later.
    def __post_init__(self):
        for who, rho, kap, n, ty in (
                ("mud", self.mud_density, self.mud_consistency,
                 self.mud_power_law_index, self.mud_yield_stress),
                ("cement", self.cement_density, self.cement_consistency,
                 self.cement_power_law_index, self.cement_yield_stress)):
            if not rho > 0.0:
                raise ValueError(f"{who}_density must be > 0, got {rho!r}")
            if not kap > 0.0:
                # kappa = 0 is not a thin fluid, it is no fluid: BF25 (2.29)
                # takes the geometric mean of the two consistencies, so a zero
                # collapses mu_e_hat and every dimensionless group with it.
                raise ValueError(f"{who}_consistency must be > 0 Pa s^n, "
                                 f"got {kap!r}")
            if not 0.0 < n <= 1.0:
                # BF25 (2.8)-(2.10), B02 (10): the gap-scale closure theory is
                # written for shear-thinning Herschel-Bulkley, n in (0, 1].
                # n > 1 is a real fluid but not a fluid this solver's closures
                # were derived or validated for, so it is refused rather than
                # silently extrapolated.
                raise ValueError(f"{who}_power_law_index must be in (0, 1] "
                                 f"(shear-thinning Herschel-Bulkley, "
                                 f"BF25 2.8-2.10), got {n!r}")
            if ty < 0.0:
                raise ValueError(f"{who}_yield_stress must be >= 0 Pa, "
                                 f"got {ty!r}")
        if not self.mean_velocities_m_s:
            raise ValueError("mean_velocities_m_s must not be empty")
        for w in self.mean_velocities_m_s:
            if not w > 0.0:
                raise ValueError("every mean velocity must be > 0 m/s, got "
                                 f"{self.mean_velocities_m_s!r}")

    # ---- derived ---------------------------------------------------------
    @property
    def mud_is_newtonian(self) -> bool:
        return self.mud_power_law_index == 1.0 and self.mud_yield_stress == 0.0

    @property
    def cement_is_newtonian(self) -> bool:
        return (self.cement_power_law_index == 1.0
                and self.cement_yield_stress == 0.0)

    @property
    def mud_viscosity(self) -> float:
        """The mud's Newtonian viscosity in Pa s.

        Refuses for a non-Newtonian mud instead of returning kappa_hat, whose
        units are Pa s^n.  This used to be a stored field; it is a property so
        that a mud given a power-law index cannot also carry a stale viscosity
        that a caller -- or a recorded config.json -- would read as the truth.
        """
        if not self.mud_is_newtonian:
            raise ValueError(
                "the mud is not Newtonian (n = "
                f"{self.mud_power_law_index}, tau_Y = {self.mud_yield_stress} "
                "Pa), so it has no single viscosity.  Use mud_consistency "
                "(Pa s^n) with mud_power_law_index, or "
                "as_fluids()[0].apparent_viscosity_at(gamma_dot).")
        return self.mud_consistency

    @property
    def mud_name(self) -> str:
        """Label for the displaced fluid, derived rather than asserted.

        Load-bearing beyond cosmetics: the name is a field of `ScaledFluid` and
        `scripts/kgep1_run.make_table` hashes the scaled fluids into the closure
        table's cache key.  The exact string for the documented Table 1 pair is
        preserved so tables cached before these fields existed still hit.
        """
        if self.mud_is_newtonian:
            if (self.mud_density, self.mud_consistency) == (998.0, 1.0e-3):
                return "drilling fluid (water)"   # Table 1 verbatim; see above
            return "drilling fluid (Newtonian)"
        if self.mud_power_law_index == 1.0:
            return "drilling fluid (Bingham)"
        if self.mud_yield_stress == 0.0:
            return "drilling fluid (power-law)"
        return "drilling fluid (Herschel-Bulkley)"

    # The pair every published number in this repository was computed with.
    # "Validated" means only that: it is the pair the gates were run against.
    # It does NOT mean the theta = 0 closure assumption is benign for it -- on
    # the production velocity axis that pair measures B-1 = 9.8%, which the
    # grading in gapscale/tables.py calls SIGNIFICANT.  Named here, next to the
    # fields, so a screen can say "this is not that pair" without restating any
    # number, and so that saying so is never mistaken for a clean bill.
    VALIDATED = (998.0, 1.0e-3, 1.0, 0.0, 1200.0, 0.6, 0.4, 1.4)

    @property
    def is_validated_pair(self) -> bool:
        return (self.mud_density, self.mud_consistency,
                self.mud_power_law_index, self.mud_yield_stress,
                self.cement_density, self.cement_consistency,
                self.cement_power_law_index,
                self.cement_yield_stress) == self.VALIDATED

    def departures_from_validated(self) -> list[str]:
        """Which fluid fields differ from the validated pair, as `name: a -> b`.

        Provenance, not physics: it compares this config's own fields against a
        recorded tuple.  A screen uses it to say what changed; it does NOT
        decide whether the result is trustworthy.  Only the per-table
        measurement -- `ClosureTable.assumption_report()` -- can do that, and
        the caller must show it.
        """
        names = ("mud_density", "mud_consistency", "mud_power_law_index",
                 "mud_yield_stress", "cement_density", "cement_consistency",
                 "cement_power_law_index", "cement_yield_stress")
        got = (self.mud_density, self.mud_consistency,
               self.mud_power_law_index, self.mud_yield_stress,
               self.cement_density, self.cement_consistency,
               self.cement_power_law_index, self.cement_yield_stress)
        return [f"{n}: {ref:g} -> {g:g}"
                for n, ref, g in zip(names, self.VALIDATED, got) if ref != g]

    @classmethod
    def from_record(cls, d: dict):
        """Rebuild from a `runio` config.json payload, migrating old schemas.

        Two shifts have to be absorbed: `mean_velocities_m_s` comes back from
        JSON as a list, and payloads written before the non-Newtonian mud fields
        existed carry `mud_viscosity`, which was Newtonian by construction.  Any
        OTHER unknown key raises, because a key this class does not understand
        means the record and the code have diverged in a way the reader must be
        told about -- dropping it would turn that into a silently wrong config.
        """
        d = dict(d)
        if "mud_viscosity" in d:
            # Pre-G-1 record: n = 1 and tau_Y = 0 were hard-coded in
            # as_fluids(), so viscosity IS the consistency for these runs.
            d.setdefault("mud_consistency", d.pop("mud_viscosity"))
            d.pop("mud_viscosity", None)
        if "mean_velocities_m_s" in d:
            d["mean_velocities_m_s"] = tuple(d["mean_velocities_m_s"])
        known = {f.name for f in fields(cls)}
        unknown = sorted(set(d) - known)
        if unknown:
            raise ValueError(
                f"recorded fluids config carries unknown key(s) {unknown}; "
                "this code cannot reproduce that run.  Check out the commit "
                "named in the payload's environment.git_sha instead of "
                "guessing.")
        return cls(**d)

    def as_fluids(self):
        """(displaced, displacing) in BF25's index convention."""
        from .scaling import HerschelBulkleyFluid
        return (HerschelBulkleyFluid(self.mud_name, self.mud_density,
                                     self.mud_consistency,
                                     self.mud_power_law_index,
                                     self.mud_yield_stress),
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
