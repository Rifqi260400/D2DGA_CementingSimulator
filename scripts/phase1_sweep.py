"""
Phase 1 -- the full 196 m K-GEP-1 open hole on the synthetic wall, plus the
A-L sweep of build spec Section 3.5.

Run:  PYTHONPATH=. python scripts/phase1_sweep.py

Why the sweep is parameterised by (m, b) and not by dimensional rheology:
K-GEP-1's mud and cement rheology and its pump schedule have not been supplied,
and inventing them would put invented numbers into the thesis' headline result.
The A-L study is a MODEL-BEHAVIOUR study -- it asks which wall geometries the
reduced model can be trusted on -- and its natural inputs are the dimensionless
groups the model actually depends on.  Those are supplied directly.  When the
real fluids arrive, `d2dga.scaling.Scaling` turns them into the same two
numbers and nothing else changes.

What each case reports:
  * the VALIDITY diagnostics first -- delta/pi and the wall gradient.  The
    Hele-Shaw reduction assumes both are small, so a case that scores well on
    displacement while failing these is not a result, it is a warning.
  * displacement efficiency at 1.2 pumped volumes and breakthrough time at
    three outlet thresholds (NUM-21).
  * the narrow-side residual, which is where a mud channel survives if one
    survives at all, and the ZF23 (3.2)-(3.6) dispersion classification.
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

sys.path.insert(0, ".")

from d2dga.config import (Config, GridConfig, StandoffConfig,  # noqa: E402
                          SyntheticWallConfig, WellConfig)
from d2dga.elliptic import NewtonianClosures                    # noqa: E402
from d2dga.geometry import (ConstantOffsetEccentricity, Geometry,  # noqa: E402
                            SinusoidalWall, synthetic_wall_feasible)
from d2dga.postprocess import (displacement_efficiency,          # noqa: E402
                               narrow_side_profile, residual_fraction,
                               zf23_metrics)
from d2dga.simulation import Simulation                          # noqa: E402

IN = 0.0254
THRESHOLDS = (0.01, 0.1, 0.5)


def build(cfg, amplitude_m, wavelength_m, mode, n_phi, n_xi):
    wall = SinusoidalWall(
        gauge_diameter_m=cfg.well.gauge_hole_diameter_m,
        amplitude_m=amplitude_m, wavelength_m=wavelength_m,
        amplitude_is_diameter=cfg.synthetic_wall.amplitude_is_diameter,
        mode=mode)
    d_gauge = 0.5 * (0.5 * cfg.well.gauge_hole_diameter_m
                     - cfg.well.casing_outer_radius_m)
    ecc = ConstantOffsetEccentricity(
        eccentricity_in_gauge_hole=cfg.standoff.eccentricity_in_gauge_hole,
        d_hat_gauge_m=d_gauge)
    return Geometry(cfg.well, wall, ecc, GridConfig(n_phi, n_xi))


def run_one(cfg, amplitude_m, wavelength_m, mode, m, b, n_phi, n_xi, cfl,
            volumes):
    ok, reason = synthetic_wall_feasible(
        cfg, amplitude_m, mode=mode,
        amplitude_is_diameter=cfg.synthetic_wall.amplitude_is_diameter)
    if not ok:
        return dict(skipped=reason)

    geo = build(cfg, amplitude_m, wavelength_m, mode, n_phi, n_xi)
    xi = geo.grid.xi_centres
    delta_pi = float(np.max(geo.narrow_gap_parameter(xi)))
    grad = float(np.max(np.abs(geo.wall_gradient(xi))))

    sim = Simulation(geo, NewtonianClosures(m), froude=1.0, delta_rho=-b,
                     cfl=cfl, inflow_concentration=1.0)
    Z = geo.grid.Z
    t0 = time.time()

    # ZF23's dispersion metrics are read off w_f(c_bar), which is recovered by
    # inverting the axial profile.  That only works while the whole front is
    # still INSIDE the domain: once the fast part has left through the outlet,
    # the profile is truncated, w_f <= 1 everywhere and sigma_w+r comes out
    # identically zero -- which is an artefact of the sampling time, not a
    # non-dispersive flow.  So snapshot the last state before the tip reaches
    # the outlet and compute the metrics from that.
    snap = {"c": None, "t": None}
    _SAMPLE_COL = int(0.8 * geo.grid.n_xi)

    def on_step(rep, c, psi):
        # Sample while the tip is still well inside the domain -- at 80% of it.
        # Waiting until the tip reaches the OUTLET clips the w_r+ tail exactly
        # where it is largest, because that tail lives at small c_bar, right at
        # the front of the profile.
        if float(np.max(c[:, _SAMPLE_COL])) < 1e-3:
            snap["c"], snap["t"] = c.copy(), rep.t

    res = sim.run(t_end=volumes * Z, record_every=1, on_step=on_step)
    return dict(
        skipped=None, Z=Z, delta_pi=delta_pi, wall_gradient=grad,
        e_min=float(np.min(geo.e(xi))), e_max=float(np.max(geo.e(xi))),
        steps=res.steps, wall=time.time() - t0,
        eta=displacement_efficiency(geo, res.concentration),
        t_br={th: res.breakthrough_at(th) / Z for th in THRESHOLDS},
        narrow_min=float(np.min(narrow_side_profile(geo, res.concentration))),
        residual=residual_fraction(geo, res.concentration),
        zf23=(zf23_metrics(geo, snap["c"], snap["t"])
              if snap["c"] is not None and snap["t"] > 0 else None),
        zf23_time=None if snap["t"] is None else snap["t"] / Z,
        cons=res.conservation_error,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-phi", type=int, default=20)
    ap.add_argument("--n-xi", type=int, default=400)
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--volumes", type=float, default=1.2)
    ap.add_argument("--m", type=float, default=0.5,
                    help="viscosity ratio eta1/eta2 (displaced/displacing)")
    ap.add_argument("--b", type=float, default=10.0,
                    help="BF25 buoyancy number -Delta_rho / Fr*^2")
    ap.add_argument("--mode", default="enlargement",
                    choices=("symmetric", "enlargement"))
    ap.add_argument("--amplitudes-in", type=float, nargs="*",
                    default=[0.0, 1.5, 3.0, 6.0])
    ap.add_argument("--wavelengths-m", type=float, nargs="*",
                    default=[5.0, 20.0, 60.0])
    ap.add_argument("--out", default="output/phase1_sweep.md")
    args = ap.parse_args()

    cfg = Config(well=WellConfig(), standoff=StandoffConfig(),
                 synthetic_wall=SyntheticWallConfig(mode=args.mode),
                 grid=GridConfig(args.n_phi, args.n_xi))

    rows = []
    for A_in in args.amplitudes_in:
        for L in args.wavelengths_m:
            r = run_one(cfg, A_in * IN, L, args.mode, args.m, args.b,
                        args.n_phi, args.n_xi, args.cfl, args.volumes)
            r["A_in"], r["L"] = A_in, L
            rows.append(r)
            if r["skipped"]:
                print(f"A={A_in:4.1f} in  L={L:5.1f} m  SKIPPED: {r['skipped']}")
                continue
            print(f"A={A_in:4.1f} in  L={L:5.1f} m  "
                  f"delta/pi={r['delta_pi']:.4f}  |dr/dxi|={r['wall_gradient']:.4f}  "
                  f"e=[{r['e_min']:.2f},{r['e_max']:.2f}]")
            print(f"    eta_E={r['eta']:.4f}  t_br "
                  + "  ".join(f"@{th}={r['t_br'][th]:.3f}" for th in THRESHOLDS)
                  + f"  narrow_min={r['narrow_min']:.3f}  residual={r['residual']:.4f}")
            z = r["zf23"]
            zs = ("n/a" if z is None else
                  f"sigma+={z.sigma_plus:.4f} |w+|={z.area_plus:.4f} "
                  f"dispersive={z.is_dispersive} @t/Z={r['zf23_time']:.3f}")
            print(f"    ZF23 {zs}   "
                  f"{r['steps']} steps, {r['wall']:.0f}s, cons {r['cons']:.1e}")

    lines = [
        "# Phase 1 -- K-GEP-1 synthetic-wall A-L sweep",
        "",
        f"Newtonian pair, `m = {args.m}`, `b = {args.b}`; wall mode "
        f"`{args.mode}`; mesh `{args.n_phi} x {args.n_xi}`; CFL {args.cfl}; "
        f"run to {args.volumes} pumped volumes.",
        "",
        "`delta/pi` and `|dr_o/dxi|` are the Hele-Shaw validity diagnostics and",
        "come FIRST on purpose: a case that displaces well while failing them is",
        "a warning, not a result. `t_br` is quoted at three outlet thresholds",
        "(assumptions.md NUM-21).  The ZF23 dispersion metrics are evaluated while",
        "the tip is still at 80% of the domain: after breakthrough the axial",
        "profile is truncated and sigma_w+r collapses to zero, and sampling right",
        "at the outlet clips the w_r+ tail exactly where it is largest.",
        "",
        "| A (in) | L (m) | delta/pi | \\|dr/dxi\\| | e range | eta_E | t_br@0.01 | @0.1 | @0.5 | narrow min | residual | sigma_w+r | dispersive |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        if r["skipped"]:
            lines.append(f"| {r['A_in']} | {r['L']} | — | — | — | — | — | — | — "
                         f"| — | — | — | skipped: {r['skipped']} |")
            continue
        t = r["t_br"]
        lines.append(
            f"| {r['A_in']} | {r['L']} | {r['delta_pi']:.4f} | "
            f"{r['wall_gradient']:.4f} | {r['e_min']:.2f}–{r['e_max']:.2f} | "
            f"{r['eta']:.4f} | {t[0.01]:.3f} | {t[0.1]:.3f} | {t[0.5]:.3f} | "
            f"{r['narrow_min']:.3f} | {r['residual']:.4f} | "
            + ("— | — |" if r["zf23"] is None else
               f"{r['zf23'].sigma_plus:.4f} | {r['zf23'].is_dispersive} |"))
    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten {args.out}")


if __name__ == "__main__":
    main()
