"""
K-GEP-1 displacement run with the fluid pair from the Materials 2025 paper.

Run:  PYTHONPATH=. python scripts/kgep1_run.py [--wall synthetic|caliper]
                                               [--w0 0.2] [--n-xi 150]

Read docs/assumptions.md "Fluids (K-GEP-1)" before quoting any number from
this.  Three things in particular travel with every result here:

  FLU-01  the source's "drilling fluid" is WATER, so m ~ 0.006 and the
          viscosity ratio is unusually favourable -- efficiencies are
          optimistic;
  FLU-04  the pump rate is a reading of that paper's inlet velocities, not a
          K-GEP-1 measurement;
  FLU-05  the pair is Herschel-Bulkley, so NUM-13's mobility-floor
          regularisation is live -- `n_static_cells` is reported on every run
          and must be looked at, not assumed zero.

The closure table is cached to output/, keyed by the fluid pair and the
buoyancy number, because building it is the expensive part.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
import time

import numpy as np

sys.path.insert(0, ".")

from d2dga.config import Config, GridConfig, IN_TO_M          # noqa: E402
from d2dga.elliptic import TabulatedClosures                   # noqa: E402
from d2dga.gapscale.tables import ClosureTable                 # noqa: E402
from d2dga.geometry import CaliperLogWall, build_geometry      # noqa: E402
from d2dga.postprocess import (displacement_efficiency,        # noqa: E402
                               narrow_side_profile, residual_fraction,
                               zf23_metrics)
from d2dga.scaling import Scaling                              # noqa: E402
from d2dga.simulation import Simulation                        # noqa: E402

THRESHOLDS = (0.01, 0.1, 0.5)


def make_geometry(cfg, wall_kind):
    if wall_kind == "synthetic":
        return build_geometry(cfg)
    wall = CaliperLogWall.from_las(
        "data/K-GEP-01_2024-03-30.las",
        total_depth_m=cfg.well.total_depth_m,
        top_m=cfg.well.casing_shoe_m, bottom_m=386.7,
        min_diameter_m=7.2 * IN_TO_M, resample_m=0.25)
    return build_geometry(cfg, wall=wall)


def make_table(sc, geo, n_c, n_h, cache_dir="output"):
    g = geo.grid
    H = geo.H(g.phi_centres, g.xi_centres)
    h_lo, h_hi = 0.9 * float(H.min()), 1.1 * float(H.max())
    # The velocity axis has to reach far beyond O(1).  With b*I1 ~ 1500 on this
    # fluid pair (FLU-06) any azimuthal tilt of the front drives |u| into the
    # hundreds, so a linear axis to 4 -- ample for the published cases -- is off
    # the end within the first timesteps.  Geometric, to 3000.
    umag_grid = np.concatenate([[0.02], np.geomspace(0.1, 3000.0, 15)])
    # the key covers EVERY axis: `load` verifies them and raises on a mismatch,
    # but a key that ignored one would turn that loud signal into a routine
    # failure the first time an axis was tuned
    key = hashlib.sha1(
        repr((sc.scaled_fluid1, sc.scaled_fluid2, sc.buoyancy_number,
              n_c, n_h, round(h_lo, 6), round(h_hi, 6),
              np.round(umag_grid, 6).tolist())).encode()).hexdigest()[:12]
    path = os.path.join(cache_dir, f"closure_table_{key}.npz")

    tab = ClosureTable(sc.scaled_fluid1, sc.scaled_fluid2,
                       c_grid=np.linspace(0.0, 1.0, n_c),
                       h_grid=np.linspace(h_lo, h_hi, n_h),
                       umag_grid=umag_grid,
                       gb_grid=np.array([sc.buoyancy_number]),
                       n_y=200, tol=1e-9)
    if os.path.exists(path):
        print(f"reusing cached closure table {path}", flush=True)
        return tab.load(path), path, 0.0
    t0 = time.time()
    tab.build(verbose=True)
    tab.save(path)
    return tab, path, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wall", default="synthetic",
                    choices=("synthetic", "caliper"))
    ap.add_argument("--w0", type=float, default=0.2)
    ap.add_argument("--n-phi", type=int, default=16)
    ap.add_argument("--n-xi", type=int, default=150)
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--volumes", type=float, default=1.2)
    ap.add_argument("--n-c", type=int, default=31)
    ap.add_argument("--n-h", type=int, default=9)
    args = ap.parse_args()

    cfg = Config(grid=GridConfig(args.n_phi, args.n_xi))
    geo = make_geometry(cfg, args.wall)
    mud, cement = cfg.fluids.as_fluids()
    sc = Scaling(mud, cement, r_a_hat_star=geo.r_a_hat_star,
                 delta_star=geo.delta_star, mean_velocity=args.w0)
    g = geo.grid
    H = geo.H(g.phi_centres, g.xi_centres)
    print(f"wall={args.wall}  Z={g.Z:.1f}  dxi={g.dxi:.3f}  "
          f"H in [{H.min():.3f}, {H.max():.3f}]  "
          f"r_a in [{geo.r_a(g.xi_centres).min():.3f}, "
          f"{geo.r_a(g.xi_centres).max():.3f}]", flush=True)
    print(f"w0={args.w0} m/s  m={sc.m:.5g}  B={sc.B:.4g}  "
          f"b={sc.buoyancy_number:.4g}  Fr*={sc.Fr_star:.4g}  "
          f"delta/pi={float(np.max(geo.narrow_gap_parameter(g.xi_centres))):.4f}",
          flush=True)

    tab, path, secs = make_table(sc, geo, args.n_c, args.n_h)
    print(f"closure table {tab.shape}, r={tab.tuned_r:g}, "
          f"{tab.n_failed} unconverged, {secs:.0f} s", flush=True)

    sim = Simulation(geo, TabulatedClosures(tab, gb=sc.buoyancy_number),
                     froude=sc.Fr_star, delta_rho=sc.delta_rho, cfl=args.cfl)
    t_end = args.volumes * g.Z
    static = [0]
    last = [time.time()]
    # see scripts/phase1_sweep.py: the ZF23 metrics only mean anything while the
    # whole front is still inside the domain
    snap = {"c": None, "t": None}
    _SAMPLE_COL = int(0.8 * geo.grid.n_xi)

    def on_step(rep, c, psi):
        static[0] = max(static[0], rep.static_cells)
        # Sample while the tip is still well inside the domain -- at 80% of it.
        # Waiting until the tip reaches the OUTLET clips the w_r+ tail exactly
        # where it is largest, because that tail lives at small c_bar, right at
        # the front of the profile.
        if float(np.max(c[:, _SAMPLE_COL])) < 1e-3:
            snap["c"], snap["t"] = c.copy(), rep.t
        if time.time() - last[0] > 120:
            last[0] = time.time()
            print(f"   t/Z={rep.t / g.Z:.3f}  eta={rep.efficiency:.4f}  "
                  f"dt={rep.dt:.4g}  step {rep.n}  picard={rep.picard_iterations}"
                  f"  static={rep.static_cells}", flush=True)

    t0 = time.time()
    res = sim.run(t_end=t_end, record_every=1, on_step=on_step)
    wall = time.time() - t0

    print(f"\n{res.reports[-1].n} steps, {wall:.0f} s, "
          f"conservation {res.conservation_error:.2e}, "
          f"c in [{res.concentration.min():.4f}, {res.concentration.max():.4f}]")
    print("t_br  " + "  ".join(
        f"@{th}={res.breakthrough_at(th) / g.Z:.4f}" for th in THRESHOLDS))
    print(f"eta_E (at {args.volumes} volumes) = "
          f"{displacement_efficiency(geo, res.concentration):.4f}")
    print(f"narrow-side minimum c_bar = "
          f"{float(np.min(narrow_side_profile(geo, res.concentration))):.4f}")
    print(f"residual (c_bar < 0.5) volume fraction = "
          f"{residual_fraction(geo, res.concentration):.4f}")
    if snap["c"] is not None and snap["t"] > 0:
        m = zf23_metrics(geo, snap["c"], snap["t"])
        print(f"ZF23 (at t/Z={snap['t'] / g.Z:.3f}, before breakthrough)  "
              f"sigma_w+r={m.sigma_plus:.4f}  |w_r+|={m.area_plus:.4f}  "
              f"dispersive={m.is_dispersive}")
    else:
        print("ZF23: no pre-breakthrough state captured")
    print(f"NUM-13 mobility floor: max static cells over the run = {static[0]}")
    if static[0]:
        print("   ^ non-zero: the regularisation is load-bearing here, and "
              "PF04 section 5 warns it flatters mud removal.")


if __name__ == "__main__":
    main()
