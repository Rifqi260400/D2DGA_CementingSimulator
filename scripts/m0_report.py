"""
M0 report: validity envelope and geometry gate summary.

Produces output/validity_envelope.png -- build spec Section 3.4 calls this a
primary thesis output, not a diagnostic: it defines where the reduced model is
expected to degrade.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from d2dga.config import IN_TO_M, Config, GridConfig
from d2dga.geometry import (ConstantOffsetEccentricity, Geometry, SinusoidalWall,
                            UniformWall, WashoutWall, synthetic_wall_feasible)

OUT = pathlib.Path(__file__).resolve().parents[1] / "output"
OUT.mkdir(exist_ok=True)

cfg = Config()
well = cfg.well
d_gauge = 0.5 * (0.5 * well.gauge_hole_diameter_m - well.casing_outer_radius_m)
ecc = ConstantOffsetEccentricity(cfg.standoff.eccentricity_in_gauge_hole, d_gauge)

# K-GEP-1 washout 195-217 m depth -> xi_hat 173-195 m above bottom hole.
# Gaussian centred at 184 m; sigma 5.0 m puts +-2 sigma at 174-194 m, i.e.
# depth 196-216 m.  ASSUMPTION GEO-06.
cases = {
    "UniformWall (gauge 10.4 in)":
        UniformWall(well.gauge_hole_diameter_m),
    "SinusoidalWall (A=1.5 in, L=20 m)":
        SinusoidalWall(well.gauge_hole_diameter_m, cfg.synthetic_wall.amplitude_m,
                       cfg.synthetic_wall.wavelength_m),
    "WashoutWall (23 in at 195-217 m)":
        WashoutWall(well.gauge_hole_diameter_m, well.washout_max_diameter_m,
                    centre_xi_hat_m=184.0, sigma_m=5.0),
}

fig, axes = plt.subplots(1, 3, figsize=(15, 6.5), sharey=True)
print("=" * 74)
print("M0 REPORT -- geometry and validity envelope")
print("=" * 74)
print(f"casing OD {well.casing_od_m*1e3:.1f} mm ({well.casing_od_m/IN_TO_M:.1f} in), "
      f"gauge hole {well.gauge_hole_diameter_m/IN_TO_M:.1f} in, "
      f"open hole {well.open_hole_length_m:.0f} m\n")

for ax, (name, wall) in zip(axes, cases.items()):
    g = Geometry(well, wall, ecc, GridConfig())
    xi = np.linspace(0.0, g.Z, 3000)
    depth = g.depth(xi)
    dp = g.narrow_gap_parameter(xi)
    hole_in = 2.0 * g.r_o_hat(xi) / IN_TO_M

    ax.plot(dp, depth, lw=1.6, color="#1f4e79", label=r"$\delta/\pi$")
    ax.axvline(0.038, color="#c0392b", ls="--", lw=1.2,
               label=r"ZF23 experiments (0.038)")
    ax.axvline(1 / np.pi, color="#7f8c8d", ls=":", lw=1.4,
               label=r"ceiling $1/\pi$")
    ax.invert_yaxis()
    ax.set_xlim(0, 0.34)
    ax.set_title(name, fontsize=9)
    ax.set_xlabel(r"$\delta/\pi$")
    ax.grid(alpha=0.3)

    twin = ax.twiny()
    twin.plot(hole_in, depth, color="#27ae60", lw=0.9, alpha=0.75)
    twin.set_xlabel("hole diameter (in)", color="#27ae60", fontsize=8)
    twin.tick_params(axis="x", labelcolor="#27ae60", labelsize=7)

    print(f"{name}")
    print(f"   delta/pi   {dp.min():.4f} .. {dp.max():.4f}   "
          f"(x{dp.max()/0.038:.1f} the validated experimental value)")
    print(f"   e          {g.e(xi).min():.4f} .. {g.e(xi).max():.4f}")
    print(f"   r_a        {g.r_a(xi).min():.4f} .. {g.r_a(xi).max():.4f}")
    print(f"   |dr_o/dxi| max {np.abs(g.wall_gradient(xi)).max():.5f}")
    vd, vh = g.annulus_volume_direct(), g.annulus_volume_from_H()
    print(f"   volume     direct {vd:.6f} m^3, via H {vh:.6f} m^3, "
          f"rel err {abs(vh-vd)/vd:.2e}")
    print()

axes[0].set_ylabel("depth below surface (m)")
axes[0].legend(loc="lower right", fontsize=7)
fig.suptitle("K-GEP-1 narrow-gap validity envelope  "
             r"$\delta/\pi = (\hat r_o-\hat r_i)/[\pi(\hat r_o+\hat r_i)]$",
             fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "validity_envelope.png", dpi=150)
print(f"wrote {OUT/'validity_envelope.png'}")

print("\n" + "-" * 74)
print("synthetic-wall feasibility, spec sweep, symmetric mode, 7 in casing")
print("-" * 74)
for A in (0.5, 1.5, 3.0, 6.0):
    for mode in ("symmetric", "enlargement"):
        ok, why = synthetic_wall_feasible(cfg, A * IN_TO_M, mode=mode)
        print(f"  A={A:4.1f} in  {mode:12s}  {'OK  ' if ok else 'FAIL'}  {why}")
