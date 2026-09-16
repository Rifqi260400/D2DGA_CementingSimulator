"""
M6-T1 -- reproduce ZF22 Table 3 with the coupled D2DGA solver.

Run:  PYTHONPATH=. python scripts/zf22_table3.py [--n-phi 20] [--n-xi 200]

ZF22's two metrics (their Section 6):
    t_br   = t_hat_br * w0_hat / L_annu_hat, with t_hat_br the time at which
             the displacing fluid first exits the annulus.  In our scaling the
             piston time is exactly Z, so t_br = t / Z.
    eta_E  = the displaced fraction of the annulus volume at the fixed time
             1.2 * L_annu_hat / w0_hat, i.e. at t = 1.2 Z.

t_br is reported at three outlet thresholds rather than one.  ZF22's wording is
exact only for a sharp front, and D2DGA fronts are not sharp: BF25 Section 3.1's
spike regime puts a vanishingly thin tip ahead of the main shock, travelling at
the Poiseuille centreline speed of 1.5 whatever the buoyancy number.  A small
threshold measures that spike, a large one the shock behind it, and for the
weakly buoyant cases the two differ by a factor of 1.5.  Reporting the bracket
is more honest than picking the threshold that flatters the comparison.

The 1-D kinematic-wave prediction (BF25 3.1, upper concave envelope of
q0 + b I3) is printed alongside, because it isolates the closures from the 2-D
solve: where the 2-D answer and the 1-D shock time disagree, the difference is
the eccentric channelling, not the gap-scale physics.
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

sys.path.insert(0, ".")

from d2dga.scaling import zf22_buoyancy_to_bf25            # noqa: E402
from tests.benchmarks.kinematic_wave import riemann_structure  # noqa: E402
from tests.benchmarks.zf22_cases import CASES, build_simulation  # noqa: E402

THRESHOLDS = (0.01, 0.1, 0.5)


def run_case(case, n_phi, n_xi, cfl):
    geo, sim, b = build_simulation(case, n_phi=n_phi, n_xi=n_xi, cfl=cfl)
    Z = geo.grid.Z
    t0 = time.time()
    res = sim.run(t_end=1.2 * Z, record_every=1)
    wall = time.time() - t0
    st = riemann_structure(case.m, b)
    shock = st["main_shock"]
    return dict(
        case=case,
        b=b,
        Z=Z,
        steps=res.reports[-1].n,
        wall=wall,
        t_br={th: res.breakthrough_at(th) / Z for th in THRESHOLDS},
        eta=res.efficiency_at(1.2 * Z),
        cons=res.conservation_error,
        c_min=float(res.concentration.min()),
        c_max=float(res.concentration.max()),
        lead=st["leading_speed"],
        shock=shock,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-phi", type=int, default=20)
    ap.add_argument("--n-xi", type=int, default=200)
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--cases", type=int, nargs="*", default=None)
    ap.add_argument("--out", default="output/zf22_table3.md")
    args = ap.parse_args()

    wanted = CASES if args.cases is None else [c for c in CASES
                                               if c.case in args.cases]
    rows = []
    for case in wanted:
        r = run_case(case, args.n_phi, args.n_xi, args.cfl)
        rows.append(r)
        th = r["t_br"]
        print(f"case {case.case:2d}  e={case.e}  m={case.m:<4} b_ZF22={case.b:<6} "
              f"b_BF25={r['b']:8.1f}  {r['steps']:6d} steps  {r['wall']:6.0f}s")
        print(f"      t_br  0.01/{th[0.01]:.3f}  0.1/{th[0.1]:.3f}  "
              f"0.5/{th[0.5]:.3f}   ZF22 {case.t_br_d2dga}  "
              f"(1-D shock {np.nan if r['shock'] is None else 1/r['shock'][2]:.3f})")
        print(f"      eta_E {r['eta']:.3f}   ZF22 {case.eta_e_d2dga}   "
              f"cons {r['cons']:.1e}   c in [{r['c_min']:.3f}, {r['c_max']:.3f}]")

    lines = [
        "# ZF22 Table 3 -- D2DGA reproduction (M6-T1)",
        "",
        f"Mesh `n_phi = {args.n_phi}`, `n_xi = {args.n_xi}`, CFL = {args.cfl}.",
        "",
        "`t_br` is reported at three outlet-concentration thresholds; see the",
        "module docstring and assumptions.md NUM-21 for why one number is not",
        "enough. `1-D shock` is the arrival time of the main shock of BF25",
        "(3.7)'s planar reduction, which isolates the gap-scale closures from",
        "the 2-D solve.",
        "",
        "| case | e | m | b (ZF22) | b (BF25) | t_br @0.01 | @0.1 | @0.5 | ZF22 t_br | 1-D shock | eta_E | ZF22 eta_E |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        c = r["case"]
        th = r["t_br"]
        sh = "--" if r["shock"] is None else f"{1 / r['shock'][2]:.3f}"
        lines.append(
            f"| {c.case} | {c.e} | {c.m} | {c.b} | {r['b']:.1f} | "
            f"{th[0.01]:.3f} | {th[0.1]:.3f} | {th[0.5]:.3f} | {c.t_br_d2dga} | "
            f"{sh} | {r['eta']:.3f} | {c.eta_e_d2dga} |")
    lines += ["", "Volume conservation (max relative drift over the run):", ""]
    for r in rows:
        lines.append(f"* case {r['case'].case}: {r['cons']:.2e}")
    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten {args.out}")


if __name__ == "__main__":
    main()
