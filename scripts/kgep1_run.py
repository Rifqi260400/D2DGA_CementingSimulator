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
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, ".")

from d2dga.config import (Config, FluidsConfig, GridConfig,   # noqa: E402
                          IN_TO_M, StandoffConfig)
from d2dga.elliptic import TabulatedClosures                   # noqa: E402
from d2dga.gapscale.tables import ClosureTable                 # noqa: E402
from d2dga.geometry import CaliperLogWall, build_geometry      # noqa: E402
from d2dga.postprocess import (cell_volume, displacement_efficiency,  # noqa: E402
                               narrow_side_efficiency, narrow_side_profile,
                               residual_fraction, tbr_relative_error,
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


def closure_umag_axis():
    """The velocity axis of the closure table.

    A function rather than a literal inside `make_table` because it is a
    PHYSICS INPUT that was invisible to the run fingerprint: it lives neither
    in `Config` nor in the CLI arguments, so changing it left the settings tag
    unchanged and a checkpoint computed with the old axis would have resumed
    silently under a new one -- precisely the hazard R-2 exists to stop.  It is
    now recorded in the config payload (see `main`) and therefore covered.
    """
    return np.concatenate([[0.0, 0.005, 0.01, 0.02],
                           np.geomspace(0.1, 3000.0, 15)])


def make_table(sc, geo, n_c, n_h, cache_dir="output"):
    g = geo.grid
    H = geo.H(g.phi_centres, g.xi_centres)
    h_lo, h_hi = 0.9 * float(H.min()), 1.1 * float(H.max())
    # The velocity axis has to reach far beyond O(1).  With b*I1 ~ 1500 on this
    # fluid pair (FLU-06) any azimuthal tilt of the front drives |u| into the
    # hundreds, so a linear axis to 4 -- ample for the published cases -- is off
    # the end within the first timesteps.  Geometric, to 3000.
    # A-3.  The axis used to start at 0.02, but |u_bar| reaches ~6e-3 at
    # stagnation points, so closure queries in the production run fell BELOW the
    # axis and were served by linear extrapolation.  |u_bar| >= 0 always, so
    # anchoring the axis at 0 makes an under-range query impossible.
    #
    # A-3b, 2026-09-23.  Anchoring it with a SINGLE node -- [0, 0.02] -- was a
    # defect, and the note that used to stand here ("changes nothing") was
    # wrong.  The gap solve at umag = 0 and 0.02 agrees to six significant
    # figures, so that segment is almost flat, while the next one, [0.02, 0.1],
    # is not.  The kink between them breaks the elliptic PICARD iteration,
    # which updates |u_bar| from the previous Psi and re-evaluates the mobility:
    # it oscillates across the kink instead of contracting.  Measured on the
    # K-GEP-1 pair at 16 x 80, same field, same everything but this axis:
    #
    #   [0.02] + geom              16 nodes   picard   6   residual 1.54e-09
    #   [0, 0.02] + geom           17 nodes   picard 100   residual 3.37e-05
    #   [0, .005, .01, .02] + geom 19 nodes   picard   6   residual 1.54e-09
    #
    # The third also leaves NO query outside the table, so it fixes A-3 without
    # the cost.  The first converges but extrapolates 36 queries per solve.
    #
    # This was invisible for five days because no test ran the Picard iteration
    # on a production-like field and the production case was not re-run; the
    # previous 16 x 80 run reported "0 steps hit the iteration cap, worst
    # residual 0.00e+00", the run on the [0, 0.02] axis hit it on 37.5% of
    # steps.  Gate M9-T1 now locks it.
    umag_grid = closure_umag_axis()
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
    ap.add_argument("--status-json", default=None, metavar="PATH",
                    help="rewrite this file with live progress every ~1.5 s "
                         "(atomic), so a monitor outside the process can read "
                         "it.  Purely additive: it is written from the existing "
                         "on_step hook and changes no computation.")

    # Fluids and standoff.  These used to live only in config.py, which meant a
    # UI could not launch a run on anything but the defaults without editing a
    # source file -- and then the run would not be reproducible from the command
    # line, which the build spec requires.  Defaults are READ FROM the
    # dataclasses, so nothing is restated here.
    fl = ap.add_argument_group(
        "fluids", "Herschel-Bulkley: tau = tau_Y + kappa gammadot^n.  Defaults "
                  "are the Materials 2025 Table 1 pair -- the only pair the "
                  "gates were run against.  Change one and the run is no longer "
                  "described by anything in AUDIT_REPORT.md; the runner says so.")
    for flag, dest, unit, helptext in FLUID_ARGS:
        fl.add_argument(flag, type=float, default=getattr(_FLUID_DEFAULTS, dest),
                        help=f"{helptext} [{unit}]")
    ap.add_argument("--eccentricity", type=float,
                    default=StandoffConfig().eccentricity_in_gauge_hole,
                    help="casing eccentricity IN GAUGE HOLE, B02 (16); the "
                         "physical offset is held constant, so e falls inside a "
                         "washout.  GEO-02: no K-GEP-1 centralizer record "
                         "exists, so this is an assumption, not a measurement.")
    return ap


_FLUID_DEFAULTS = FluidsConfig()

# (flag, dest, unit, help).  One row per fluid field, so the parser, the UI form
# and the reproducibility check cannot disagree about what exists.
FLUID_ARGS = (
    ("--mud-density", "mud_density", "kg/m^3",
     "displaced fluid density rho1_hat"),
    ("--mud-consistency", "mud_consistency", "Pa s^n",
     "displaced fluid consistency kappa1_hat (= viscosity when n = 1)"),
    ("--mud-power-law-index", "mud_power_law_index", "-",
     "displaced fluid n1, in (0, 1]; n > 1 is refused, not extrapolated"),
    ("--mud-yield-stress", "mud_yield_stress", "Pa",
     "displaced fluid tauY1_hat"),
    ("--cement-density", "cement_density", "kg/m^3",
     "displacing fluid density rho2_hat"),
    ("--cement-consistency", "cement_consistency", "Pa s^n",
     "displacing fluid consistency kappa2_hat"),
    ("--cement-power-law-index", "cement_power_law_index", "-",
     "displacing fluid n2, in (0, 1]"),
    ("--cement-yield-stress", "cement_yield_stress", "Pa",
     "displacing fluid tauY2_hat"),
)


def fluids_from_args(args) -> FluidsConfig:
    """`FluidsConfig` from the parsed arguments, validated by the dataclass.

    The validation is `FluidsConfig.__post_init__`'s, not a second copy: an
    incoherent pair raises here exactly as it would if written into config.py.
    """
    return FluidsConfig(**{dest: getattr(args, dest)
                           for _, dest, _, _ in FLUID_ARGS})


def config_from_args(args) -> Config:
    """The whole `Config` a set of arguments describes."""
    return Config(grid=GridConfig(args.n_phi, args.n_xi),
                  fluids=fluids_from_args(args),
                  standoff=StandoffConfig(args.eccentricity))


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
    # e < 1 is the geometric limit (the casing touches the wall at e = 1/2 in
    # B02's convention once d_hat is accounted for); 0.95 keeps H > 0.
    "eccentricity": {"unit": "-", "lo": 0.0, "hi": 0.95},
}


def _physics_suffix(cfg) -> str:
    """`""` for the shipped fluids and standoff, `_f<6 hex>` otherwise."""
    default = Config(grid=cfg.grid)
    if (cfg.fluids, cfg.standoff) == (default.fluids, default.standoff):
        return ""
    key = repr((cfg.fluids, cfg.standoff)).encode()
    return "_f" + hashlib.sha1(key).hexdigest()[:6]


def write_status_stub(path, phase, note=""):
    """A status document for the part of a run that happens BEFORE stepping.

    The closure table is built first and can take a quarter of an hour; the
    step loop -- and therefore `StatusWriter` -- does not exist yet.  Without
    this, a monitor watching a freshly launched run sees no file at all and
    cannot distinguish "building the table" from "failed to start", which are
    very different things to a person deciding whether to wait.
    """
    if not path:
        return
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump({"schema": 1, "state": "preparing", "phase": phase,
                   "note": note, "pid": os.getpid(),
                   "command": " ".join(sys.argv),
                   "updated_at": time.time()}, f)
    os.replace(tmp, path)


class StatusWriter:
    """Rewrites a small JSON file with live progress, atomically.

    Why a file and not a socket or a queue: a run is already a plain subprocess
    and must stay one, because the build spec requires a UI run and a CLI run to
    be the same invocation.  A file is the only channel that costs the solver
    nothing and works whether the reader is a UI, a tail, or nobody.

    Atomic by temp-file plus `os.replace`, the pattern `Simulation.run` already
    uses for checkpoints, so a reader can never see half a document.  Throttled,
    so a 10^5-step run does not spend its time serialising.  Every field here is
    read off the objects the run already maintains; nothing is computed for the
    file.
    """

    def __init__(self, path, sim, geo, capacity, t_end, command, every=1.5):
        self.path, self.sim, self.geo = path, sim, geo
        self.capacity, self.t_end = capacity, t_end
        self.command, self.every = command, every
        self.t0, self.last = time.time(), 0.0
        self.static_max = 0

    def _write(self, doc):
        tmp = self.path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(doc, f)
        os.replace(tmp, self.path)

    def update(self, rep, force=False):
        self.static_max = max(self.static_max, rep.static_cells)
        now = time.time()
        if not force and now - self.last < self.every:
            return
        self.last = now
        Z = self.geo.grid.Z
        frac = min(1.0, rep.t / self.t_end) if self.t_end > 0 else 0.0
        self._write({
            "schema": 1, "state": "running", "pid": os.getpid(),
            "command": self.command,
            "started_at": self.t0, "updated_at": now,
            "elapsed_s": now - self.t0,
            "fraction": frac,
            "eta_wall_s": ((now - self.t0) * (1.0 - frac) / frac
                           if frac > 1e-6 else None),
            "step": rep.n, "t": rep.t, "t_over_Z": rep.t / Z, "dt": rep.dt,
            "efficiency": rep.efficiency,
            "c_min": rep.c_min, "c_max": rep.c_max,
            "picard_iterations": rep.picard_iterations,
            "static_cells": rep.static_cells, "static_cells_max": self.static_max,
            # the invariants the run must be watched by, not just its progress
            "conservation_error": self.sim.conservation_error,
            "outlet_import_volumes": self.sim.transport.outlet_import / self.capacity,
            "picard_unconverged": self.sim.n_picard_unconverged,
            "worst_picard_residual": self.sim.worst_picard_residual,
            "closure_table_range_report": self.sim.closures.range_report(),
            # BCF25 IV, live.  The invariant this stage's verification rests
            # on should not first be visible when the run ends.
            "bcf25_worst": self.sim.ledger.worst,
            "bcf25_total_flux_imbalance": self.sim.ledger.total_flux_imbalance,
        })

    def finish(self, state, message=""):
        try:
            doc = json.load(open(self.path)) if os.path.exists(self.path) else {}
        except (OSError, ValueError):
            doc = {}
        doc.update({"state": state, "message": message,
                    "updated_at": time.time(),
                    "elapsed_s": time.time() - self.t0})
        self._write(doc)


def main():
    args = build_parser().parse_args()

    write_status_stub(args.status_json, "geometry", "building the wall profile")
    cfg = config_from_args(args)
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

    write_status_stub(args.status_json, "closure table",
                      f"{args.n_c} x {args.n_h} x 17 gap-scale solves; this is "
                      f"the long part, and it is cached")
    tab, path, secs = make_table(sc, geo, args.n_c, args.n_h)
    print(f"closure table {tab.shape}, r={tab.tuned_r:g}, "
          f"{tab.n_failed} unconverged, {secs:.0f} s", flush=True)

    sim = Simulation(geo, TabulatedClosures(tab, gb=sc.buoyancy_number),
                     froude=sc.Fr_star, delta_rho=sc.delta_rho, cfl=args.cfl,
                     inflow=args.inflow)
    t_end = args.volumes * g.Z
    static = [0]
    last = [time.time()]
    capacity = float(np.sum(cell_volume(geo)))
    status = None
    if args.status_json:
        status = StatusWriter(args.status_json, sim, geo, capacity, t_end,
                              " ".join(sys.argv))
        write_status_stub(args.status_json, "stepping",
                          "elliptic solve + transport step per timestep")
    # see scripts/phase1_sweep.py: the ZF23 metrics only mean anything while the
    # whole front is still inside the domain
    snap = {"c": None, "t": None}
    _SAMPLE_COL = int(0.8 * geo.grid.n_xi)

    def on_step(rep, c, psi):
        static[0] = max(static[0], rep.static_cells)
        if status is not None:
            status.update(rep)
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
    # The filename carried six settings and none of them a fluid property, so
    # two runs differing only in the mud would write to the SAME path.  Since
    # R-2 that is refused rather than silently resumed -- correct, but it also
    # made the second run impossible.  A short hash of the non-mesh physics is
    # appended when it is not the shipped default, so the default keeps its
    # historical name (old checkpoints still resolve) and anything else gets its
    # own file.  The hash is a disambiguator, not the guard: the guard is the
    # settings tag stored inside the checkpoint.
    phys = _physics_suffix(cfg)
    ckpt = os.path.join("output", f"kgep1_ckpt_{args.wall}_{args.inflow}_"
                                  f"{args.n_phi}x{args.n_xi}_w{args.w0}_"
                                  f"v{args.volumes}{phys}.npz")
    # Provenance and the resume guard.  Before this, the checkpoint filename
    # carried six of ~20 settings, so editing a fluid property and relaunching
    # the same command SILENTLY resumed a field computed under the old physics
    # -- `run` checked only the grid shape.  The payload is written beside the
    # checkpoint and its fingerprint is stored inside it.
    # A-3b.  The closure table's velocity axis is a physics input that lives in
    # neither Config nor the arguments.  Recording it here puts it inside the
    # fingerprint, so a checkpoint cannot be resumed across a change to it.
    _args = dict(vars(args))
    _args["closure_umag_axis"] = [float(x) for x in closure_umag_axis()]
    payload = runio.config_payload(cfg, _args)
    level, message = runio.check_resume(ckpt, payload)
    if level == "block":
        raise SystemExit("REFUSING TO RESUME\n  " + message)
    if level == "warn":
        print("WARNING: " + message, flush=True)
    tag = runio.fingerprint(payload)
    print(f"settings fingerprint {tag}", flush=True)
    runio.write_config(ckpt, payload)

    try:
        res = sim.run(t_end=t_end, record_every=1, on_step=on_step,
                      checkpoint_path=ckpt, checkpoint_tag=tag)
    except BaseException as exc:
        # A monitor that only ever sees "running" cannot tell a crashed run from
        # a slow one, and a run killed at 90% is the expensive case.  Record the
        # reason, then re-raise unchanged -- this handler must not swallow it.
        if status is not None:
            status.finish("failed", f"{type(exc).__name__}: {exc}")
        raise
    wall = time.time() - t0

    if not res.reports:
        print("\nnothing to do: the checkpoint is already at or past "
              f"{args.volumes} volumes (step {res.steps}, "
              f"t/Z = {res.times[-1] / g.Z:.4f}).  The results below are read "
              "back from it, not recomputed.", flush=True)
    print(f"\n{res.steps} steps, {wall:.0f} s, "
          f"conservation {res.conservation_error:.2e}, "
          f"c in [{res.concentration.min():.4f}, {res.concentration.max():.4f}]")
    # A-2.  t_br is a LEVEL SET of a front the first-order scheme smears, so its
    # error is O(dxi^(1/2)) -- half order -- while eta_E, an integral, is
    # O(dxi).  Measured against the exact rarefaction
    # (test_a2_threshold_front_is_half_order_and_integrals_are_first_order):
    # ~18% at the 0.01 threshold on n_xi = 80, ~3% at the 0.5 threshold.  The
    # figures are printed with that attached rather than to four decimals,
    # because halving the 0.01 error costs four times the cells.
    print("t_br  " + "  ".join(
        f"@{th}={res.breakthrough_at(th) / g.Z:.4f}"
        f"(+-{tbr_relative_error(th, g.n_xi) * 100:.0f}%)" for th in THRESHOLDS))
    print("      the bracket is the DISCRETISATION error of a threshold front, "
          "O(dxi^1/2); eta_E below is an integral and is O(dxi)")
    eta = displacement_efficiency(geo, res.concentration)
    print(f"eta_E (at {args.volumes} volumes) = {eta:.4f}")

    # Volume balance.  Before breakthrough what is in the annulus cannot exceed
    # what was pumped; after it, it cannot exceed it either.  This is the check
    # that caught NUM-29 -- the inflow face was delivering 110 Q -- and it is
    # cheap, so it is printed on every run rather than done by hand when
    # something already looks wrong.
    pumped = args.volumes * g.Z / capacity          # Q = 1, so volume = t
    tb = res.breakthrough_at(0.01) / g.Z
    note = ("expected negative: displacing fluid has been leaving through the "
            f"outlet since t/Z = {tb:.3f}" if np.isfinite(tb) and tb < args.volumes
            else "no breakthrough, so this must be <= 0 up to the startup "
                 "transient")
    print(f"volume balance: pumped {pumped:.4f}, present {eta:.4f}, "
          f"difference {eta - pumped:+.5f} "
          f"({100 * (eta - pumped) / args.volumes:+.2f}% of the job) -- {note}")
    # BCF25 Section IV, the verification this stage of the project rests on:
    # there is no CFD comparison, so the model is checked against the published
    # ZF22 cases and against this ledger -- the same test BCF25 apply to the
    # same scheme.  Reported with THEIR normalisation (by the annulus volume)
    # so the number can be put beside their ~1e-15, and next to the older
    # `cons_err`, which divides by the running volume instead and is therefore
    # a stricter early-time measure of the same thing.
    print(sim.ledger.report(), flush=True)
    print(f"  (the run's own cons_err, normalised by the volume present "
          f"rather than by the annulus, is {res.conservation_error:.2e})")
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
    # eta_E is dominated by the wide side, which was never in doubt.  The same
    # integral over the narrow quarter is the number the job is actually about.
    print(f"eta_N (narrow quarter) = "
          f"{narrow_side_efficiency(geo, res.concentration):.4f}")
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
    # A-3b.  This line existed and nothing read it.  A closure-table axis
    # change put 37.5% of steps at the cap for five days without anyone
    # noticing, because the only evidence was a count at the end of a log.
    # It is now stated as a fraction, with a verdict, and recorded in
    # metrics.json where the Validation screen reads it.
    _cap_frac = sim.n_picard_unconverged / max(res.steps, 1)
    print(f"Picard: {sim.n_picard_unconverged} of {res.steps} steps hit the "
          f"iteration cap ({_cap_frac:.1%}), worst residual "
          f"{sim.worst_picard_residual:.2e} (tolerance {sim.picard_tol:.0e})")
    if _cap_frac > 0.05:
        print("   ^^ WARNING: the elliptic solve is NOT converging on a large "
              "fraction of steps.  Psi is perturbed by O(residual) each time. "
              "This result needs looking at, not quoting.  A-3b was exactly "
              "this, caused by the closure table's velocity axis.")
    elif sim.n_picard_unconverged:
        print("   ^^ non-zero but small; the next step largely corrects an "
              "under-converged solve, and the pre-A-3 production run reported "
              "0 over 85474 steps, so a rise here is worth tracing.")
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
        "steps": res.steps,
        "wall_clock_s": wall,
        "conservation_error": res.conservation_error,
        "conservation_bcf25": {
            k: v for k, v in sim.ledger.as_dict().items()
            # the full series goes in the checkpoint, not in this file: a
            # 10^5-step run would make metrics.json tens of megabytes and
            # nothing reads it from here
            if k not in ("times", "err1", "err2", "imbalance")},
        "c_min": float(res.concentration.min()),
        "c_max": float(res.concentration.max()),
        "eta_E": eta,
        "t_br_over_Z": {str(th): (res.breakthrough_at(th) / g.Z) for th in THRESHOLDS},
        "t_br_relative_error": {str(th): tbr_relative_error(th, g.n_xi)
                                for th in THRESHOLDS},
        "eta_N": narrow_side_efficiency(geo, res.concentration),
        "eta_N_fraction": 0.25,
        "narrow_side_min": float(np.min(narrow_side_profile(geo, res.concentration))),
        "residual_fraction": residual_fraction(geo, res.concentration),
        "volume_pumped": args.volumes * g.Z / capacity,
        "annulus_capacity": capacity,
        "zf23": zf,
        "static_cells_max": static[0],
        "picard_unconverged": sim.n_picard_unconverged,
        "picard_unconverged_fraction": sim.n_picard_unconverged / max(res.steps, 1),
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
    if status is not None:
        status.finish("done", os.path.basename(ckpt))


if __name__ == "__main__":
    main()
