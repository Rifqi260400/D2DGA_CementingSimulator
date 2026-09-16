"""
M0 -- wellbore geometry.

Turns a wall profile plus a casing programme into the fields the rest of the
solver needs: r_a(xi), delta(xi), e(xi), H(phi, xi), and the validity
diagnostic delta/pi.

All geometric definitions follow B02 Section 2.1:

    r_hat_a(xi_hat) = (r_o + r_i)/2                             B02 (6)
    d_hat(xi_hat)   = (r_o - r_i)/2                             B02 (6)
    r_hat_a_star    = (1/Z_hat) * int r_hat_a  d xi_hat         B02 (13)
    xi              = (xi_hat - xi_hat_bh) / (pi r_hat_a_star)  B02 (14)
    delta(xi_hat)   = d_hat / r_hat_a                           B02 (15)
    delta_star      = (1/Z_hat) * int delta   d xi_hat          B02 (15)
    r_a(xi)         = r_hat_a / r_hat_a_star                    B02 (16)
    e(xi)           = e_hat / (2 d_hat)                         B02 (16)
    H(phi, xi)      = delta(xi) r_a(xi) [1 + e(xi) cos(pi phi)] / delta_star
                                                                B02 (21)

B02 (21) is the *only* place in the published literature where H is allowed to
vary with xi; ZF22 / ZF23 / BF25 / BCF25 all specialise to r_a = 1 and
H = 1 + e cos(pi phi).  Setting a UniformWall here must reproduce that
specialisation to machine precision -- that is test M0-T1.

Array convention: fields indexed [phi, xi].  xi_hat is measured UPWARD from
bottom hole, so depth_below_surface = total_depth - xi_hat.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from .config import (IN_TO_M, Config, GridConfig, StandoffConfig,
                     SyntheticWallConfig, WellConfig)


class GeometryError(ValueError):
    """Raised when a wall profile produces a physically impossible annulus."""


# ===========================================================================
# Wall profiles
# ===========================================================================
class WallProfile(ABC):
    """Outer (formation) wall of the annulus."""

    @abstractmethod
    def outer_radius(self, xi_hat: np.ndarray) -> np.ndarray:
        """Outer radius in METRES at axial position(s) xi_hat, in metres,
        measured upward from bottom hole."""

    def describe(self) -> str:
        return self.__class__.__name__


@dataclass(frozen=True)
class UniformWall(WallProfile):
    """Constant gauge hole.  The reduction case: must recover the published
    uniform-annulus geometry exactly (M0-T1)."""

    diameter_m: float

    def outer_radius(self, xi_hat):
        xi_hat = np.asarray(xi_hat, dtype=float)
        return np.full_like(xi_hat, 0.5 * self.diameter_m)

    def describe(self):
        return f"UniformWall(D={self.diameter_m:.5f} m)"


@dataclass(frozen=True)
class SinusoidalWall(WallProfile):
    """
    Gauge hole with a controlled periodic perturbation.

    This is the workhorse for the A-L sweep of build spec Section 3.5.  Because
    amplitude and wavelength are free, it gives a controlled study the real
    caliper cannot: it says in advance which parts of the measured log the
    reduced model can be trusted on.

    mode='symmetric'   : r_o = r_gauge + A_r sin(2 pi xi_hat / L)
                         the literal form written in the spec.
    mode='enlargement' : r_o = r_gauge + A_r (1 - cos(2 pi xi_hat / L)) / 2
                         one-sided, never below gauge.  Real holes wash out but
                         do not wash in, and the spec's sweep to A = 6 in is
                         geometrically impossible in the symmetric form with
                         7 in casing.  See assumptions.md GEO-04.
    """

    gauge_diameter_m: float
    amplitude_m: float
    wavelength_m: float
    amplitude_is_diameter: bool = True
    mode: str = "symmetric"

    @property
    def _amplitude_radius_m(self) -> float:
        return 0.5 * self.amplitude_m if self.amplitude_is_diameter else self.amplitude_m

    def outer_radius(self, xi_hat):
        xi_hat = np.asarray(xi_hat, dtype=float)
        k = 2.0 * np.pi / self.wavelength_m
        a = self._amplitude_radius_m
        if self.mode == "symmetric":
            perturbation = a * np.sin(k * xi_hat)
        elif self.mode == "enlargement":
            perturbation = a * 0.5 * (1.0 - np.cos(k * xi_hat))
        else:
            raise ValueError(f"unknown mode {self.mode!r}")
        return 0.5 * self.gauge_diameter_m + perturbation

    def describe(self):
        basis = "diam" if self.amplitude_is_diameter else "rad"
        return (f"SinusoidalWall(D={self.gauge_diameter_m:.5f} m, "
                f"A={self.amplitude_m:.5f} m [{basis}], L={self.wavelength_m:.1f} m, "
                f"{self.mode})")


@dataclass(frozen=True)
class WashoutWall(WallProfile):
    """Gauge hole with a single smooth Gaussian enlargement -- the K-GEP-1
    195-217 m feature.  Smooth and C-infinity, so it never violates the
    slow-variation assumption through a discontinuity; only through its
    gradient, which is what we want to measure."""

    gauge_diameter_m: float
    peak_diameter_m: float
    centre_xi_hat_m: float
    sigma_m: float

    def outer_radius(self, xi_hat):
        xi_hat = np.asarray(xi_hat, dtype=float)
        bulge = 0.5 * (self.peak_diameter_m - self.gauge_diameter_m)
        z = (xi_hat - self.centre_xi_hat_m) / self.sigma_m
        return 0.5 * self.gauge_diameter_m + bulge * np.exp(-0.5 * z * z)

    def describe(self):
        return (f"WashoutWall(D={self.gauge_diameter_m:.5f}->"
                f"{self.peak_diameter_m:.5f} m at xi={self.centre_xi_hat_m:.1f} m)")


class CaliperLogWall(WallProfile):
    """
    Measured caliper log.

    The K-GEP-1 log has not been supplied yet, so the FILE PARSERS ARE STUBS.
    `from_arrays` is fully implemented and tested, so once the log is read into
    two arrays this class works; the parsers only have to reach that point.

    Deliberately no format is assumed.  LAS, CSV and digitised two-column
    depth/diameter all dispatch to their own reader, each of which raises until
    an actual file has been inspected.

    Swapping SinusoidalWall -> CaliperLogWall must touch no module but this one.
    That is asserted by test M0-T6.
    """

    def __init__(self, xi_hat_m: np.ndarray, outer_diameter_m: np.ndarray,
                 source: str = "<arrays>"):
        xi_hat_m = np.asarray(xi_hat_m, dtype=float)
        outer_diameter_m = np.asarray(outer_diameter_m, dtype=float)
        if xi_hat_m.ndim != 1 or xi_hat_m.shape != outer_diameter_m.shape:
            raise ValueError("xi_hat_m and outer_diameter_m must be 1-D, same length")
        order = np.argsort(xi_hat_m)
        self._xi = xi_hat_m[order]
        self._r = 0.5 * outer_diameter_m[order]
        if not np.all(np.isfinite(self._r)):
            raise GeometryError("caliper log contains non-finite radii; "
                                "gaps must be filled before construction")
        self.source = source

    @classmethod
    def from_arrays(cls, xi_hat_m, outer_diameter_m, source="<arrays>"):
        return cls(xi_hat_m, outer_diameter_m, source)

    # -- parsers: stubs until the real file is in hand ----------------------
    @classmethod
    def from_las(cls, path, curve=None, total_depth_m=None):
        raise NotImplementedError(
            "LAS parser stubbed: the K-GEP-1 log has not been supplied. "
            "Needs the caliper curve mnemonic, its unit (in vs mm), the depth "
            "reference, and the null value. Route to from_arrays once known.")

    @classmethod
    def from_csv(cls, path, depth_col=None, diameter_col=None, total_depth_m=None):
        raise NotImplementedError(
            "CSV parser stubbed: column names, units and depth convention "
            "unknown until the file is inspected. Route to from_arrays.")

    @classmethod
    def from_two_column(cls, path, total_depth_m=None):
        raise NotImplementedError(
            "Digitised two-column parser stubbed: delimiter, header presence, "
            "units and depth direction unknown. Route to from_arrays.")

    def outer_radius(self, xi_hat):
        xi_hat = np.asarray(xi_hat, dtype=float)
        # Linear interpolation, clamped at the ends.  A higher-order scheme
        # would manufacture detail the log does not resolve (features under
        # ~0.5 m are not resolved -- assumptions.md GEO-05).
        return np.interp(xi_hat, self._xi, self._r)

    def describe(self):
        return f"CaliperLogWall(n={self._xi.size}, source={self.source})"


# ===========================================================================
# Eccentricity models
# ===========================================================================
class EccentricityModel(ABC):
    """Supplies the casing offset e_hat(xi_hat) in metres."""

    @abstractmethod
    def offset(self, xi_hat: np.ndarray, d_hat: np.ndarray) -> np.ndarray:
        ...

    def describe(self) -> str:
        return self.__class__.__name__


@dataclass(frozen=True)
class ConstantOffsetEccentricity(EccentricityModel):
    """
    Constant physical offset e_hat, the default in the absence of centralizer
    records.  e_hat is pinned by naming the eccentricity the casing would have
    in GAUGE hole.

    Consequence, and it is the physically right one: inside a washout d_hat
    grows while e_hat does not, so e = e_hat/(2 d_hat) FALLS.  The casing does
    not move just because the hole got bigger.
    """

    eccentricity_in_gauge_hole: float
    d_hat_gauge_m: float

    @property
    def offset_m(self) -> float:
        return 2.0 * self.eccentricity_in_gauge_hole * self.d_hat_gauge_m

    def offset(self, xi_hat, d_hat):
        return np.full_like(np.asarray(d_hat, dtype=float), self.offset_m)

    def describe(self):
        return (f"ConstantOffsetEccentricity(e_gauge="
                f"{self.eccentricity_in_gauge_hole:.3f}, "
                f"e_hat={self.offset_m*1e3:.2f} mm)")


@dataclass(frozen=True)
class ConstantEccentricity(EccentricityModel):
    """Fixed e regardless of gap -- i.e. the casing slides outward as the hole
    enlarges.  Not physical for a real wellbore, but it is what every published
    uniform-annulus case assumes, so it is the model needed to reproduce them."""

    value: float

    def offset(self, xi_hat, d_hat):
        return 2.0 * self.value * np.asarray(d_hat, dtype=float)

    def describe(self):
        return f"ConstantEccentricity(e={self.value:.3f})"


# ===========================================================================
# Grid
# ===========================================================================
class Grid:
    """
    Staggered rectangular mesh, BCF25 Fig. 2.

    Concentration lives at cell centres, Psi at cell corners, velocities on
    faces.  phi in [0, 1] spans HALF the annulus (0 = wide side, 1 = narrow),
    the symmetry that all published 2D work assumes (build spec Q6).
    """

    def __init__(self, n_phi: int, n_xi: int, Z: float):
        self.n_phi = int(n_phi)
        self.n_xi = int(n_xi)
        self.Z = float(Z)

        self.phi_edges = np.linspace(0.0, 1.0, self.n_phi + 1)
        self.xi_edges = np.linspace(0.0, self.Z, self.n_xi + 1)
        self.phi_centres = 0.5 * (self.phi_edges[:-1] + self.phi_edges[1:])
        self.xi_centres = 0.5 * (self.xi_edges[:-1] + self.xi_edges[1:])

        self.dphi = 1.0 / self.n_phi
        self.dxi = self.Z / self.n_xi

    def __repr__(self):
        return (f"Grid(n_phi={self.n_phi}, n_xi={self.n_xi}, Z={self.Z:.3f}, "
                f"dphi={self.dphi:.4g}, dxi={self.dxi:.4g})")


# ===========================================================================
# Geometry
# ===========================================================================
class Geometry:
    """
    Assembled annular geometry on a grid.

    Construction order matters: r_hat_a_star and delta_star are integrals over
    the whole open-hole interval (B02 13, 15), so they must be formed before
    the scaled coordinate xi and hence before the grid exists.  We therefore
    quadrature them on a dense reference grid in xi_hat first.
    """

    def __init__(self,
                 well: WellConfig,
                 wall: WallProfile,
                 eccentricity: EccentricityModel,
                 grid_cfg: GridConfig,
                 n_reference: int = 20001):
        self.well = well
        self.wall = wall
        self.ecc_model = eccentricity

        self.r_i_hat = well.casing_outer_radius_m
        self.Z_hat = well.open_hole_length_m

        # ---- reference quadrature for the two global averages -------------
        xi_ref = np.linspace(0.0, self.Z_hat, n_reference)
        r_o_ref = np.asarray(wall.outer_radius(xi_ref), dtype=float)
        self._check_wall(xi_ref, r_o_ref)

        r_a_hat_ref = 0.5 * (r_o_ref + self.r_i_hat)
        d_hat_ref = 0.5 * (r_o_ref - self.r_i_hat)
        delta_ref = d_hat_ref / r_a_hat_ref

        self.r_a_hat_star = float(np.trapezoid(r_a_hat_ref, xi_ref) / self.Z_hat)
        self.delta_star = float(np.trapezoid(delta_ref, xi_ref) / self.Z_hat)

        # ---- scaled axial extent, B02 (14) --------------------------------
        self.Z = self.Z_hat / (np.pi * self.r_a_hat_star)
        self.grid = Grid(grid_cfg.n_phi, grid_cfg.n_xi, self.Z)

        self._validate_eccentricity(xi_ref, d_hat_ref)

    # -- coordinate mapping -------------------------------------------------
    def xi_hat(self, xi):
        """Scaled axial coordinate -> metres above bottom hole.  B02 (14)."""
        return np.asarray(xi, dtype=float) * np.pi * self.r_a_hat_star

    def xi_of(self, xi_hat):
        return np.asarray(xi_hat, dtype=float) / (np.pi * self.r_a_hat_star)

    def depth(self, xi):
        """Scaled axial coordinate -> true depth below surface, metres."""
        return self.well.depth_of(self.xi_hat(xi))

    # -- dimensional fields at scaled xi -----------------------------------
    def r_o_hat(self, xi):
        return np.asarray(self.wall.outer_radius(self.xi_hat(xi)), dtype=float)

    def r_a_hat(self, xi):
        return 0.5 * (self.r_o_hat(xi) + self.r_i_hat)

    def d_hat(self, xi):
        return 0.5 * (self.r_o_hat(xi) - self.r_i_hat)

    # -- dimensionless fields, B02 (15)-(16) -------------------------------
    def delta(self, xi):
        return self.d_hat(xi) / self.r_a_hat(xi)

    def r_a(self, xi):
        return self.r_a_hat(xi) / self.r_a_hat_star

    def e(self, xi):
        d = self.d_hat(xi)
        e_hat = self.ecc_model.offset(self.xi_hat(xi), d)
        return e_hat / (2.0 * d)

    # -- gap half-width, B02 (21) ------------------------------------------
    def H(self, phi, xi):
        """
        H(phi, xi) = delta(xi) r_a(xi) [1 + e(xi) cos(pi phi)] / delta_star

        Returns shape (len(phi), len(xi)) when both are 1-D.
        """
        phi = np.atleast_1d(np.asarray(phi, dtype=float))
        xi = np.atleast_1d(np.asarray(xi, dtype=float))
        axial = (self.delta(xi) * self.r_a(xi) / self.delta_star)[None, :]
        azimuthal = 1.0 + self.e(xi)[None, :] * np.cos(np.pi * phi)[:, None]
        return axial * azimuthal

    # -- validity diagnostic, build spec 3.4 -------------------------------
    def narrow_gap_parameter(self, xi):
        """
        delta/pi = (r_o - r_i) / [pi (r_o + r_i)].

        Bounded above by 1/pi ~ 0.3183 for any r_i > 0: the numerator can never
        exceed the bracket.  Frigaard's validated experiments sit at ~0.038.
        """
        r_o = self.r_o_hat(xi)
        return (r_o - self.r_i_hat) / (np.pi * (r_o + self.r_i_hat))

    def wall_gradient(self, xi):
        """d r_o_hat / d xi_hat, dimensionless.  The other half of the
        slow-axial-variation question (build spec Q8)."""
        xh = self.xi_hat(xi)
        step = max(1e-4, 1e-6 * self.Z_hat)
        fwd = np.asarray(self.wall.outer_radius(xh + step), dtype=float)
        bwd = np.asarray(self.wall.outer_radius(xh - step), dtype=float)
        return (fwd - bwd) / (2.0 * step)

    def div_a_f(self, phi, xi):
        """
        div_a . f, evaluated numerically.  Build spec M0-T4.

        BF25 (2.18) drops the 1/Fr*^2 term on the grounds that this vanishes.
        docs/derivation.md Finding 1 gives the closed form under variable
        geometry:

            div_a . f = [ r_a cos(beta) dbeta/dxi + sin(beta) dr_a/dxi ] sin(pi phi)

        Both terms carry a sin(beta) or a dbeta/dxi, so at beta = 0 this is
        identically zero NO MATTER how violently r_a varies.  For an inclined
        well with a caliper wall it is not, and the term must be kept.

        We evaluate by finite difference rather than substituting the closed
        form, so that the test actually exercises the geometry rather than
        restating the algebra.
        """
        phi = np.atleast_1d(np.asarray(phi, dtype=float))
        xi = np.atleast_1d(np.asarray(xi, dtype=float))
        beta = self.well.inclination_rad

        # f_phi = r_a cos(beta) has no phi dependence -> d f_phi / d phi = 0.
        d_fphi_dphi = np.zeros((phi.size, xi.size))

        h = max(1e-6, 1e-7 * self.Z)
        g = lambda z: self.r_a(z) * np.sin(beta)
        d_ra_sinbeta_dxi = (g(xi + h) - g(xi - h)) / (2.0 * h)
        d_fxi_dxi = np.sin(np.pi * phi)[:, None] * d_ra_sinbeta_dxi[None, :]

        return d_fphi_dphi / self.r_a(xi)[None, :] + d_fxi_dxi

    def buoyancy_direction(self, phi, xi):
        """f = (f_phi, f_xi) = (r_a cos beta, r_a sin beta sin(pi phi)).

        Note the crossed labelling inherited from BF25 (2.11)-(2.12): the
        phi-component carries cos(beta) (the axial gravity direction) because
        the buoyancy enters the momentum balance rotated, as (f_xi, -f_phi).
        """
        phi = np.atleast_1d(np.asarray(phi, dtype=float))
        xi = np.atleast_1d(np.asarray(xi, dtype=float))
        beta = self.well.inclination_rad
        ra = self.r_a(xi)[None, :]
        f_phi = np.broadcast_to(ra * np.cos(beta), (phi.size, xi.size)).copy()
        f_xi = ra * np.sin(beta) * np.sin(np.pi * phi)[:, None]
        return f_phi, f_xi

    # -- volume, for M0-T2 --------------------------------------------------
    def annulus_volume_direct(self, n=200001):
        """
        Exact annular volume, integral of pi (r_o^2 - r_i^2) d xi_hat.

        The cross-sectional area between two non-intersecting circles is
        independent of their offset, so this holds for any eccentricity.
        """
        xh = np.linspace(0.0, self.Z_hat, n)
        r_o = np.asarray(self.wall.outer_radius(xh), dtype=float)
        return float(np.trapezoid(np.pi * (r_o**2 - self.r_i_hat**2), xh))

    def annulus_volume_from_H(self, n_phi=4001, n_xi=4001):
        """
        Volume as the gap-averaged model sees it.

            V = 4 pi^2 delta_star r_hat_a_star^3  int int H r_a  d phi d xi

        Derivation: the physical half-annulus cross-section is
        int_0^pi 2 d_hat (1 + e cos theta) r_hat_a d theta = 2 pi d_hat r_hat_a,
        so the full section is 4 pi d_hat r_hat_a -- which equals
        pi(r_o^2 - r_i^2) identically.  Pushing that through the B02 scaling
        gives the prefactor above.  Agreement with annulus_volume_direct is
        therefore a real test of the scaling chain, not a tautology.
        """
        phi = np.linspace(0.0, 1.0, n_phi)
        xi = np.linspace(0.0, self.Z, n_xi)
        integrand = self.H(phi, xi) * self.r_a(xi)[None, :]
        inner = np.trapezoid(integrand, phi, axis=0)
        total = float(np.trapezoid(inner, xi))
        return 4.0 * np.pi**2 * self.delta_star * self.r_a_hat_star**3 * total

    # -- guards -------------------------------------------------------------
    def _check_wall(self, xi_ref, r_o_ref):
        if not np.all(np.isfinite(r_o_ref)):
            raise GeometryError("wall profile returned non-finite radii")
        gap = r_o_ref - self.r_i_hat
        if np.any(gap <= 0.0):
            bad = int(np.argmin(gap))
            raise GeometryError(
                f"annular gap non-positive: r_o - r_i = {gap[bad]*1e3:.3f} mm at "
                f"xi_hat = {xi_ref[bad]:.2f} m (depth "
                f"{self.well.depth_of(xi_ref[bad]):.2f} m). "
                f"Casing OD {self.well.casing_od_m*1e3:.1f} mm is larger than the "
                f"hole there.")

    def _validate_eccentricity(self, xi_ref, d_hat_ref):
        e_hat = self.ecc_model.offset(xi_ref, d_hat_ref)
        e = e_hat / (2.0 * d_hat_ref)
        if np.any(e >= 1.0):
            bad = int(np.argmax(e))
            raise GeometryError(
                f"eccentricity e = {e[bad]:.4f} >= 1 at xi_hat = {xi_ref[bad]:.2f} m "
                f"(depth {self.well.depth_of(xi_ref[bad]):.2f} m): casing contacts "
                f"the wall. B02 assumption (i) requires e < 1.")
        self.e_max = float(np.max(e))
        self.e_min = float(np.min(e))

    def summary(self) -> str:
        xi = self.grid.xi_centres
        dp = self.narrow_gap_parameter(xi)
        return "\n".join([
            f"wall              : {self.wall.describe()}",
            f"eccentricity      : {self.ecc_model.describe()}",
            f"casing r_i        : {self.r_i_hat*1e3:.3f} mm",
            f"r_hat_a_star      : {self.r_a_hat_star*1e3:.3f} mm",
            f"delta_star        : {self.delta_star:.6f}",
            f"Z (scaled)        : {self.Z:.3f}",
            f"grid              : {self.grid!r}",
            f"e range           : {self.e_min:.4f} .. {self.e_max:.4f}",
            f"delta/pi range    : {dp.min():.4f} .. {dp.max():.4f}  "
            f"(ceiling 1/pi = {1/np.pi:.4f}; ZF23 experiments 0.038)",
        ])


# ===========================================================================
# Convenience builders
# ===========================================================================
def build_geometry(cfg: Config, wall: WallProfile | None = None,
                   eccentricity: EccentricityModel | None = None) -> Geometry:
    """Assemble a Geometry from config, defaulting to the synthetic wall of
    build spec Section 3.5."""
    if wall is None:
        sw: SyntheticWallConfig = cfg.synthetic_wall
        wall = SinusoidalWall(
            gauge_diameter_m=cfg.well.gauge_hole_diameter_m,
            amplitude_m=sw.amplitude_m,
            wavelength_m=sw.wavelength_m,
            amplitude_is_diameter=sw.amplitude_is_diameter,
            mode=sw.mode,
        )
    if eccentricity is None:
        st: StandoffConfig = cfg.standoff
        d_gauge = 0.5 * (0.5 * cfg.well.gauge_hole_diameter_m
                         - cfg.well.casing_outer_radius_m)
        eccentricity = ConstantOffsetEccentricity(
            eccentricity_in_gauge_hole=st.eccentricity_in_gauge_hole,
            d_hat_gauge_m=d_gauge,
        )
    return Geometry(cfg.well, wall, eccentricity, cfg.grid)


# ===========================================================================
# Feasibility of a synthetic wall (build spec Section 3.5 A-L sweep)
# ===========================================================================
def synthetic_wall_feasible(cfg: Config, amplitude_m: float, mode: str = "symmetric",
                            amplitude_is_diameter: bool = True,
                            eccentricity_in_gauge_hole: float | None = None):
    """
    Can this (amplitude, mode) actually be built with the configured casing?

    Returns (ok: bool, reason: str).

    Two independent limits bite, in this order:

      1. positive gap   : min(r_o) > r_i
      2. e < 1          : min(d_hat) > e_hat/2, with e_hat fixed by the gauge-hole
                          eccentricity (ConstantOffsetEccentricity)

    Limit 2 is the tighter one and is a direct consequence of holding the casing
    offset fixed: pinching the hole while the casing stays put drives e up.
    In 'enlargement' mode the wall never goes below gauge, so neither limit can
    be violated and every amplitude is admissible.
    """
    well = cfg.well
    r_i = well.casing_outer_radius_m
    D_g = well.gauge_hole_diameter_m
    e_g = (cfg.standoff.eccentricity_in_gauge_hole
           if eccentricity_in_gauge_hole is None else eccentricity_in_gauge_hole)
    a_r = 0.5 * amplitude_m if amplitude_is_diameter else amplitude_m

    if mode == "enlargement":
        return True, "enlargement mode never goes below gauge"

    r_o_min = 0.5 * D_g - a_r
    gap = r_o_min - r_i
    if gap <= 0.0:
        return False, (f"negative gap: min hole {2*r_o_min/IN_TO_M:.2f} in vs casing "
                       f"{D_g and well.casing_od_m/IN_TO_M:.2f} in")
    d_min = 0.5 * gap
    d_gauge = 0.5 * (0.5 * D_g - r_i)
    e_hat_half = e_g * d_gauge
    e_max = e_hat_half / d_min
    if e_max >= 1.0:
        return False, (f"e = {e_max:.2f} >= 1 at the pinch "
                       f"(d_hat = {d_min*1e3:.2f} mm vs offset/2 = "
                       f"{e_hat_half*1e3:.2f} mm)")
    return True, f"ok (max e = {e_max:.3f}, min gap = {gap/IN_TO_M:.3f} in)"


def max_symmetric_amplitude_m(cfg: Config, eccentricity_in_gauge_hole=None):
    """Largest symmetric diameter amplitude admitting e < 1 everywhere."""
    well = cfg.well
    r_i = well.casing_outer_radius_m
    D_g = well.gauge_hole_diameter_m
    e_g = (cfg.standoff.eccentricity_in_gauge_hole
           if eccentricity_in_gauge_hole is None else eccentricity_in_gauge_hole)
    d_gauge = 0.5 * (0.5 * D_g - r_i)
    return D_g - 2.0 * (r_i + 2.0 * e_g * d_gauge)
