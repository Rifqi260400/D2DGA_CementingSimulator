"""Validity envelope from the measured K-GEP-1 caliper (M0-T5, real data)."""
import sys, pathlib, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from d2dga.config import IN_TO_M, Config, GridConfig
from d2dga.geometry import (CaliperLogWall, ConstantOffsetEccentricity, Geometry,
                            SinusoidalWall)

OUT = pathlib.Path(__file__).resolve().parents[1] / "output"; OUT.mkdir(exist_ok=True)
cfg = Config(); well = cfg.well
d_gauge = 0.5 * (0.5 * well.gauge_hole_diameter_m - well.casing_outer_radius_m)
ecc = ConstantOffsetEccentricity(cfg.standoff.eccentricity_in_gauge_hole, d_gauge)

rep = {}
wall = CaliperLogWall.from_las("data/K-GEP-01_2024-03-30.las", total_depth_m=390.0,
                               top_m=194.0, bottom_m=386.7,
                               min_diameter_m=7.2 * IN_TO_M, resample_m=0.25, report=rep)
synth = SinusoidalWall(well.gauge_hole_diameter_m, cfg.synthetic_wall.amplitude_m,
                       cfg.synthetic_wall.wavelength_m)

fig, axes = plt.subplots(1, 4, figsize=(17, 8), sharey=True)
for ax, (name, w) in zip(axes[:2], [("measured caliper", wall), ("synthetic (A=1.5in, L=20m)", synth)]):
    g = Geometry(well, w, ecc, GridConfig(20, 400))
    xi = np.linspace(0, g.Z, 4000); depth = g.depth(xi)
    ax.plot(2 * g.r_o_hat(xi) / IN_TO_M, depth, lw=0.7, color="#27ae60")
    ax.axvline(well.casing_od_m / IN_TO_M, color="k", ls="--", lw=1,
               label=f"casing OD {well.casing_od_m/IN_TO_M:.0f} in")
    ax.axvline(10.43, color="#7f8c8d", ls=":", lw=1, label="gauge 10.43 in")
    ax.set_xlabel("hole diameter (in)"); ax.set_title(name, fontsize=9)
    ax.grid(alpha=0.3); ax.invert_yaxis(); ax.legend(fontsize=7)

g_real = Geometry(well, wall, ecc, GridConfig(20, 400))
g_syn = Geometry(well, synth, ecc, GridConfig(20, 400))
xi = np.linspace(0, g_real.Z, 4000)
axes[2].plot(g_real.narrow_gap_parameter(xi), g_real.depth(xi), lw=0.7,
             color="#1f4e79", label="measured")
xs = np.linspace(0, g_syn.Z, 4000)
axes[2].plot(g_syn.narrow_gap_parameter(xs), g_syn.depth(xs), lw=0.9,
             color="#95a5a6", label="synthetic")
axes[2].axvline(0.038, color="#c0392b", ls="--", lw=1.2, label="ZF23 validated 0.038")
axes[2].axvline(1/np.pi, color="#7f8c8d", ls=":", lw=1.2, label=r"ceiling $1/\pi$")
axes[2].set_xlabel(r"$\delta/\pi$"); axes[2].set_title("validity envelope", fontsize=9)
axes[2].grid(alpha=0.3); axes[2].legend(fontsize=7); axes[2].set_xlim(0, 0.34)

axes[3].plot(np.abs(g_real.wall_gradient(xi)), g_real.depth(xi), lw=0.7,
             color="#8e44ad", label="measured")
axes[3].plot(np.abs(g_syn.wall_gradient(xs)), g_syn.depth(xs), lw=0.9,
             color="#95a5a6", label="synthetic")
axes[3].set_xlabel(r"$|d\hat r_o/d\hat\xi|$"); axes[3].set_xscale("log")
axes[3].set_title("slow-variation stress", fontsize=9)
axes[3].grid(alpha=0.3); axes[3].legend(fontsize=7)
axes[0].set_ylabel("depth below surface (m)")
fig.suptitle("K-GEP-1: measured caliper vs synthetic wall", fontsize=12)
fig.tight_layout(); fig.savefig(OUT / "caliper_validity_envelope.png", dpi=150)

dp = g_real.narrow_gap_parameter(xi); gr = np.abs(g_real.wall_gradient(xi))
print("MEASURED CALIPER, K-GEP-1, casing OD {:.0f} in".format(well.casing_od_m/IN_TO_M))
print("  samples used {} of {} raw ({} null, {} outside interval)".format(
    rep['n_used'], rep['n_raw'], rep['n_null'], rep['n_raw']-rep['n_null']-rep['n_in_interval']))
if rep.get('n_below_min'):
    print("  {} samples rejected below the minimum, {:.2f}-{:.2f} m".format(
        rep['n_below_min'], *rep['below_min_span_m']))
else:
    print("  no sub-minimum samples inside the interval (the 386.76-389.95 m "
          "tool artefact is already excluded by bottom_m=386.7)")
print("  hole diameter  {:.2f} .. {:.2f} in".format(*[x/IN_TO_M for x in rep["diameter_range_m"]]))
print("  delta/pi       min {:.4f}  median {:.4f}  p95 {:.4f}  max {:.4f}"
      .format(dp.min(), np.median(dp), np.percentile(dp,95), dp.max()))
print("     fraction above ZF23's 0.038: {:.1f}%".format(100*(dp>0.038).mean()))
print("  |dr_o/dxi|     median {:.5f}  p95 {:.5f}  max {:.5f}".format(
    np.median(gr), np.percentile(gr,95), gr.max()))
print("  e range        {:.3f} .. {:.3f}".format(g_real.e(xi).min(), g_real.e(xi).max()))
vd, vh = g_real.annulus_volume_direct(), g_real.annulus_volume_from_H()
print("  annulus volume {:.4f} m^3 (via H {:.4f}, rel err {:.2e})".format(vd, vh, abs(vh-vd)/vd))
print("\nSYNTHETIC for comparison")
gs = np.abs(g_syn.wall_gradient(xs)); ds = g_syn.narrow_gap_parameter(xs)
print("  delta/pi max {:.4f}   |dr_o/dxi| max {:.5f}".format(ds.max(), gs.max()))
print("\nwrote", OUT / "caliper_validity_envelope.png")
