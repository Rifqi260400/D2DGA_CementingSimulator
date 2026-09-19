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
from d2dga.postprocess import (cell_volume, displacement_efficiency,  # noqa: E402
                               narrow_side_profile, residual_fraction,
                               zf23_metrics)
from d2dga import runio                                        # noqa: E402
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
    # A-3.  The axis used to start at 0.02, but |u_bar| reaches 1e-3 at
    # stagnation points, so 315 of 25.9M closure queries in the production run
    # fell BELOW the axis and were served by linear extrapolation.  |u_bar| >= 0
    # always, so anchoring the axis at 0 makes an under-range query impossible.
    # It costs one extra plane (6% of the build) and changes nothing: the gap
    # solve at umag = 0, 1e-3 and 0.02 agrees to six significant figures on this
    # fluid pair, which is why the extrapolation was harmless here and why the
    # DETECTION (tables.range_report) mattered more than the extrapolation.
    umag_grid = np.concatenate([[0.0, 0.02], np.geomspace(0.1, 3000.0, 15)])
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


def build_parser() -> argparse.ArgumentParser:
    """
    The run parameters, defined once.

    Extracted from `main` so that the UI can INTROSPECT it -- names, defaults,
    types, `choices`, help text -- instead of restating them in a second place
    where they would drift.  `ui/` reads this; nothing about a parameter is
    written down twice.  Units and physical ranges that argparse cannot express
    live in `metadata()` below, next to it rather than in the UI layer.
    """
    ap = argparse.ArgumentParser(
        description="D2DGA cement displacement on K-GEP-1.")
    ap.add_argument("--wall", default="synthetic",
                    choices=("synthetic", "caliper"),
                    help="wall profile: the synthetic sinusoid, or the measured "
                         "caliper log (delta/pi 0.168 against 0.083, ~7x the cost)")
    ap.add_argument("--w0", type=float, default=0.2,
                    help="annular mean velocity w0_hat, m/s")
    ap.add_argument("--n-phi", type=int, default=16,
                    help="azimuthal cells over the HALF annulus")
    ap.add_argument("--n-xi", type=int, default=150,
                    help="axial cells; t_br error is O(dxi^1/2), see A-2")
    ap.add_argument("--cfl", type=float, default=0.5,
                    help="multiplier on the BCF25 (44) monotonicity bound")
    ap.add_argument("--volumes", type=float, default=1.2,
                    help="pumped volume, in annulus volumes")
    ap.add_argument("--no-resume", action="store_true",
                    help="delete any existing checkpoint and start over")
    ap.add_argument("--n-c", type=int, default=31,
                    help="closure-table concentration nodes (NUM-03 open)")
    ap.add_argument("--n-h", type=int, default=9,
                    help="closure-table half-gap nodes")
    ap.add_argument("--inflow", default="no_axial_gradient",
                    choices=("no_axial_gradient", "uniform"),
                    help="bottom-hole condition, assumptions.md NUM-04.  B02 "
                         "(70) permits backflow through the shoe; a real shoe "
                         "does not, which is what 'uniform' forbids.")
    return ap


# Units and admissible ranges, which argparse has no field for.  Kept here, next
# to the parser, so the UI reads both from one module.  `lo`/`hi` are inclusive
# bounds the solver itself enforces or that are physically meaningless outside.
PARAM_META = {
    "w0":      {"unit": "m/s",      "lo": 1e-4, "hi": 5.0},
    "n_phi":   {"unit": "cells",    "lo": 2,    "hi": 256},
    "n_xi":    {"unit": "cells",    "lo": 2,    "hi": 2000},
    "cfl":     {"unit": "-",        "lo": 1e-3, "hi": 1.0},
    "volumes": {"unit": "volumes",  "lo": 1e-3, "hi": 10.0},
    "n_c":     {"unit": "nodes",    "lo": 2,    "hi": 401},
    "n_h":     {"unit": "nodes",    "lo": 1,    "hi": 65},
}


def main():
    args = build_parser().parse_args()

    cfg = Config(grid=GridConfig(args.n_phi, args.n_xi))
    geo = make_geometry(cfg, args.wall)
    if args.no_resume:
        for f in os.listdir("output"):
            if f.startswith("kgep1_ckpt_"):
                os.remove(os.path.join("output", f))
    mud, cement = cfg.fluids.as_fluids()
    # Say it here, before any number is printed, and in the same words the UI
    # uses: every figure in AUDIT_REPORT.md / docs/gate_status.md was computed
    # with the Materials 2025 Table 1 pair, and a changed fluid invalidates the
    # comparison rather than merely shifting it.  The B-1 measurement for THIS
    # pair is printed later, from the table that was actually built.
    print(f"fluids: {mud.name} -> {cement.name}")
    _dep = cfg.fluids.departures_from_validated()
    if _dep:
        print("WARNING: not the validated fluid pair -- " + "; ".join(_dep))
        print("         the figures in AUDIT_REPORT.md, docs/gate_status.md and "
              "output/kgep1_results.md do NOT describe this run.")
        print("         read the angle-sensitivity (B-1) line below before "
              "quoting anything: it is benign for the shipped pair only.")
    sc = Scaling(mud, cement, r_a_hat_star=geo.r_a_hat_star,
                 delta_star=geo.delta_star, mean_velocity=args.w0)
    g = geo.grid
    H = geo.H(g.phi_centres, g.xi_centres)
    print(f"inflow={args.inflow}")
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
                     froude=sc.Fr_star, delta_rho=sc.delta_rho, cfl=args.cfl,
                     inflow=args.inflow)
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
    ckpt = os.path.join("output", f"kgep1_ckpt_{args.wall}_{args.inflow}_"
                                  f"{args.n_phi}x{args.n_xi}_w{args.w0}_"
                                  f"v{args.volumes}.npz")
    # Provenance and the resume guard.  Before this, the checkpoint filename
    # carried six of ~20 settings, so editing a fluid property and relaunching
    # the same command SILENTLY resumed a field computed under the old physics
    # -- `run` checked only the grid shape.  The payload is written beside the
    # checkpoint and its fingerprint is stored inside it.
    payload = runio.config_payload(cfg, vars(args))
    level, message = runio.check_resume(ckpt, payload)
    if level == "block":
        raise SystemExit("REFUSING TO RESUME\n  " + message)
    if level == "warn":
        print("WARNING: " + message, flush=True)
    tag = runio.fingerprint(payload)
    print(f"settings fingerprint {tag}", flush=True)
    runio.write_config(ckpt, payload)

    res = sim.run(t_end=t_end, record_every=1, on_step=on_step,
                  checkpoint_path=ckpt, checkpoint_tag=tag)
    wall = time.time() - t0

    print(f"\n{res.reports[-1].n} steps, {wall:.0f} s, "
          f"conservation {res.conservation_error:.2e}, "
          f"c in [{res.concentration.min():.4f}, {res.concentration.max():.4f}]")
    # A-2.  t_br is a LEVEL SET of a front the first-order scheme smears, so its
    # error is O(dxi^(1/2)) -- half order -- while eta_E, an integral, is
    # O(dxi).  Measured against the exact rarefaction
    # (test_a2_threshold_front_is_half_order_and_integrals_are_first_order):
    # ~18% at the 0.01 threshold on n_xi = 80, ~3% at the 0.5 threshold.  The
    # figures are printed with that attached rather than to four decimals,
    # because halving the 0.01 error costs four times the cells.
    _TBR_REL_ERR = {0.01: 0.175, 0.1: 0.069, 0.5: 0.033}      # at n_xi = 80
    scale = (80.0 / g.n_xi) ** 0.5                            # half order
    print("t_br  " + "  ".join(
        f"@{th}={res.breakthrough_at(th) / g.Z:.4f}"
        f"(+-{_TBR_REL_ERR[th] * scale * 100:.0f}%)" for th in THRESHOLDS))
    print("      the bracket is the DISCRETISATION error of a threshold front, "
          "O(dxi^1/2); eta_E below is an integral and is O(dxi)")
    eta = displacement_efficiency(geo, res.concentration)
    print(f"eta_E (at {args.volumes} volumes) = {eta:.4f}")

    # Volume balance.  Before breakthrough what is in the annulus cannot exceed
    # what was pumped; after it, it cannot exceed it either.  This is the check
    # that caught NUM-29 -- the inflow face was delivering 110 Q -- and it is
    # cheap, so it is printed on every run rather than done by hand when
    # something already looks wrong.
    capacity = float(np.sum(cell_volume(geo)))
    pumped = args.volumes * g.Z / capacity          # Q = 1, so volume = t
    tb = res.breakthrough_at(0.01) / g.Z
    note = ("expected negative: displacing fluid has been leaving through the "
            f"outlet since t/Z = {tb:.3f}" if np.isfinite(tb) and tb < args.volumes
            else "no breakthrough, so this must be <= 0 up to the startup "
                 "transient")
    print(f"volume balance: pumped {pumped:.4f}, present {eta:.4f}, "
          f"difference {eta - pumped:+.5f} "
          f"({100 * (eta - pumped) / args.volumes:+.2f}% of the job) -- {note}")
    # A-3.  The closure table's range guard used to be a UserWarning and nothing
    # more, and these runs were launched under `-W ignore::UserWarning`, so the
    # evidence that the closures were being EXTRAPOLATED was destroyed.  The
    # record is now state on the table and is printed unconditionally here.
    # A-1.  Gross volume of displacing fluid that entered through the OUTLET
    # while the smeared front was crossing it.  It is an O(dxi) artefact and
    # must fall roughly in half when n_xi doubles; if it does not, the outflow
    # treatment has a real defect.  Quoted against the annulus capacity.
    print(f"outlet import (A-1 artefact) = "
          f"{sim.transport.outlet_import / capacity:.5f} volumes "
          f"({100 * sim.transport.outlet_import / capacity / args.volumes:+.3f}% "
          f"of the job); expect it to halve when n_xi doubles", flush=True)
    # B-1.  The closure table stores the theta = 0 slice (u_bar parallel to
    # G~_b).  For a Newtonian pair that is exact; for a yield-stress pair it can
    # be an O(1) error.  The table measures it for ITS OWN fluid pair at build
    # time, and it is printed here so no number can be quoted without it.
    print(tab.assumption_report(), flush=True)
    rr = sim.closures.range_report()
    print("closure-table range: " + (rr if rr else "all queries inside the table"),
          flush=True)
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
    print(f"Picard: {sim.n_picard_unconverged} steps hit the iteration cap, "
          f"worst residual {sim.worst_picard_residual:.2e} "
          f"(tolerance {sim.picard_tol:.0e})")
    if static[0]:
        print("   ^ non-zero: the regularisation is load-bearing here, and "
              "PF04 section 5 warns it flatters mud removal.")

    # Everything printed above, written as JSON beside the checkpoint.  Same
    # values, same postprocess calls -- nothing is recomputed differently for
    # the file, so the file and the log cannot disagree.
    zf = None
    if snap["c"] is not None and snap["t"] > 0:
        _m = zf23_metrics(geo, snap["c"], snap["t"])
        zf = {"sampled_at_t_over_Z": snap["t"] / g.Z,
              "sigma_plus": _m.sigma_plus, "sigma_minus": _m.sigma_minus,
              "area_plus": _m.area_plus, "area_minus": _m.area_minus,
              "n_bins": _m.n_bins, "is_dispersive": bool(_m.is_dispersive)}
    runio.write_metrics(ckpt, {
        "schema": 1,
        "settings_fingerprint": tag,
        "checkpoint": os.path.basename(ckpt),
        "steps": res.reports[-1].n,
        "wall_clock_s": wall,
        "conservation_error": res.conservation_error,
        "c_min": float(res.concentration.min()),
        "c_max": float(res.concentration.max()),
        "eta_E": eta,
        "t_br_over_Z": {str(th): (res.breakthrough_at(th) / g.Z) for th in THRESHOLDS},
        "t_br_relative_error": {"0.01": 0.175 * (80.0 / g.n_xi) ** 0.5,
                                "0.1": 0.069 * (80.0 / g.n_xi) ** 0.5,
                                "0.5": 0.033 * (80.0 / g.n_xi) ** 0.5},
        "narrow_side_min": float(np.min(narrow_side_profile(geo, res.concentration))),
        "residual_fraction": residual_fraction(geo, res.concentration),
        "volume_pumped": args.volumes * g.Z / capacity,
        "annulus_capacity": capacity,
        "zf23": zf,
        "static_cells_max": static[0],
        "picard_unconverged": sim.n_picard_unconverged,
        "worst_picard_residual": sim.worst_picard_residual,
        "outlet_import_volumes": sim.transport.outlet_import / capacity,
        "closure_table_range_report": sim.closures.range_report(),
        "closure_table_assumption_report": tab.assumption_report(),
        "scaling": {"m": sc.m, "B": sc.B, "b": sc.buoyancy_number,
                    "Fr_star": sc.Fr_star, "delta_rho": sc.delta_rho,
                    "w0_hat": sc.w0_hat, "Z": g.Z},
        "geometry": {"delta_over_pi_min": float(np.min(geo.narrow_gap_parameter(g.xi_centres))),
                     "delta_over_pi_max": float(np.max(geo.narrow_gap_parameter(g.xi_centres))),
                     "e_min": float(np.min(geo.e(g.xi_centres))),
                     "e_max": float(np.max(geo.e(g.xi_centres))),
                     "H_min": float(H.min()), "H_max": float(H.max())},
    })
    print(f"wrote {runio.config_path_for(ckpt)}")
    print(f"wrote {runio.metrics_path_for(ckpt)}")


if __name__ == "__main__":
    main()
