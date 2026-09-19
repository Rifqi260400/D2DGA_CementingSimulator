"""Screen 4 — Validation.

What the model has been shown to reproduce, what it has not, and what nobody
has checked. All three, with equal prominence: a validation page that shows
only the passes is an advertisement.

Every number is read from a file something else produced.

  * `output/gates.json` — written by `tests/run_gates.py`, a pytest plugin.
    Absent means the gates have not been run in this checkout; the page says
    so and does **not** fall back to the prose in `docs/gate_status.md`.
  * `output/zf22_table3.md` — the ten-case comparison, written by
    `scripts/zf22_table3.py`. Parsed, not transcribed.
  * `BLOCKED.md` — the open questions.

The mockup showed 10/10 on ZF22. That is not true and the page says so: case 1
is not reproduced (BLK-6), and the table below shows the numbers side by side
rather than a tick.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time

import numpy as np
import streamlit as st

from ui import figures, runs
from ui.common import (ROOT, banner, page_chrome, panel, rows,
                       unavailable)
from ui.reports import read_blocked, read_gates, zf22_comparison

page_chrome("Validation")
st.title("Validation")

GATES = os.path.join(ROOT, "output", "gates.json")
ZF22 = os.path.join(ROOT, "output", "zf22_table3.md")


# --------------------------------------------------------------------------
# the gates
# --------------------------------------------------------------------------
banner("warn",
       "<strong>This stage does not validate against CFD.</strong> There is no "
       "resolved three-dimensional simulation to compare with, and none is "
       "claimed. What is on this page is <em>verification</em>: the scheme is "
       "checked against the ten published ZF22 cases, and against the volume "
       "ledger of <strong>BCF25 §IV</strong> — the same test that paper applies "
       "to the same scheme, on the same equations, reporting ~10<sup>−15</sup>. "
       "Verification asks whether the equations are solved correctly. "
       "Validation against a measurement or a resolved simulation is a "
       "different question, and this page does not answer it.")

st.markdown("### Conservation — BCF25 §IV")

st.markdown(
    '<span class="note">BCF25 compute each fluid\'s volume two ways: '
    '<strong>Vol₁</strong> by integrating c&#772;<sub>k</sub> over the interior '
    '(multiplying by H·r<sub>a</sub>·Δφ·Δξ and summing), and <strong>Vol₂</strong> '
    'by starting from the initial volume and adding the inflow and subtracting '
    'the outflow at each timestep. Both are normalised by the total annulus '
    'volume and differenced. Their Figs. 7 and 15 plot that against time, one '
    'curve per fluid, and report it stays at ~10<sup>−15</sup> across five '
    'cases and both models. This is the same computation on this '
    'solver.</span>', unsafe_allow_html=True)

_led_runs = [r for r in runs.discover("output") if r.ledger is not None]
if not _led_runs:
    unavailable("Conservation ledger",
                "No run on disk carries one. It is written into the checkpoint "
                "from 2026-09-19; earlier runs recorded only the single "
                "running maximum, under a different normalisation "
                "(by the volume present, not by the annulus), so their number "
                "is <em>not</em> comparable with BCF25's and is not shown here "
                "as though it were. Re-run to record it.")
else:
    pick = st.selectbox("Run", _led_runs, key="val_ledger",
                        format_func=lambda r: r.name.replace("_", " · "))
    led = pick.ledger
    geo_v = runs.build_geometry(pick)
    e1 = np.abs(np.asarray(led["err1"], dtype=float))
    e2 = np.abs(np.asarray(led["err2"], dtype=float))
    im = np.abs(np.asarray(led["imbalance"], dtype=float))
    worst = float(max(e1.max(), e2.max()))

    c1, c2 = st.columns([1.0, 1.9], gap="medium")
    with c1:
        panel("Worst relative error",
              f'<div class="big">{worst:.2e}</div>'
              + rows((("fluid 1 (displaced)", f"{e1.max():.2e}"),
                      ("fluid 2 (displacing)", f"{e2.max():.2e}"),
                      ("total flux in − out", f"{im.max():.2e}"),
                      ("steps", f"{len(e2) - 1:,}"),
                      ("BCF25 report", "~1e-15"))),
              "normalised by the total annulus volume — BCF25's "
              "normalisation, so the numbers are comparable")
        # The criterion is the ROUNDOFF SCALE, not a flat constant.  BCF25's
        # 1e-15 is a number for their runs; error accumulates with the step
        # count, so a 15,000-step run cannot be held to the same figure and
        # a comparison against 1e-15 flat would mark a correct long run as
        # failing.  sqrt(N)*eps is the random-walk estimate; the ratio to it
        # is printed, because that ratio -- not the raw exponent -- is what
        # says whether this is arithmetic or bookkeeping.
        n_steps = len(e2) - 1
        roundoff = max(n_steps, 1) ** 0.5 * float(np.finfo(float).eps)
        ratio = worst / roundoff if roundoff else float("inf")
        panel("Against the roundoff scale",
              rows((("steps N", f"{n_steps:,}"),
                    ("√N·ε", f"{roundoff:.2e}"),
                    ("worst / √N·ε", f"{ratio:.2f}"))),
              "ε = 2.22e-16. A random walk of N roundoff errors grows as "
              "√N·ε, so this ratio is the scale-free reading; the raw "
              "exponent is not comparable between runs of different length.")
        if ratio <= 20.0:
            banner("good",
                   f"<strong>Arithmetic, not bookkeeping.</strong> "
                   f"{worst:.1e} over {n_steps:,} steps is "
                   f"<strong>{ratio:.1f}×</strong> the random-walk roundoff "
                   f"scale. BCF25 report ~1e-15 for their cases; this run is "
                   f"longer, so the comparable statement is the ratio, not "
                   f"the exponent.")
        else:
            banner("bad",
                   f"<strong>{worst:.1e} is {ratio:.0f}× the roundoff "
                   f"scale</strong> for {n_steps:,} steps ({roundoff:.1e}). "
                   f"That is too far past arithmetic to be arithmetic: look "
                   f"at the flux bookkeeping.")
        if im.max() == 0.0:
            banner("good", "<strong>The total volumetric flux balances "
                           "exactly.</strong> In and out agree to the bit, so "
                           "the elliptic solve delivers the same Q at both "
                           "ends.")
    with c2:
        st.pyplot(figures.conservation(led["times"], led["err1"], led["err2"],
                                       led["imbalance"], geo_v.grid.Z),
                  use_container_width=True)
        banner("warn",
               "<strong>The two fluid curves are near mirror images, and that "
               "is not corroboration.</strong> BCF25 evolve K = 3 "
               "concentrations independently, so their three curves are three "
               "separate ledgers and summing to 1 is a real result. Here K = 2 "
               "and only c&#772;<sub>2</sub> is evolved: c&#772;<sub>1</sub> ≡ "
               "1 − c&#772;<sub>2</sub> pointwise, and the two-fluid closure "
               "makes fluid 1's face flux identically the total minus fluid "
               "2's — the LLF dissipation term changes sign with c and cancels "
               "exactly. So fluid 1's ledger differs from −fluid 2's only by "
               "the <em>total-flux imbalance</em>, which is the third curve and "
               "the independent content of this check.")

st.markdown("### Test gates")

gates = read_gates()

if gates is None:
    unavailable("Gate results",
                "No <code>output/gates.json</code>. Produce one with"
                "<br><code>PYTHONPATH=. python tests/run_gates.py</code>"
                "<br>This page will not substitute the prose in "
                "<code>docs/gate_status.md</code> for a run: that file is "
                "written by hand and cannot go stale loudly.")
    if st.button("Run the gates now (several minutes)", key="val_run"):
        with st.spinner("pytest…"):
            subprocess.run([sys.executable, "tests/run_gates.py"], cwd=ROOT,
                           env=dict(os.environ, PYTHONPATH=ROOT))
        st.rerun()
else:
    counts = gates.get("counts", {})
    passed, failed = counts.get("passed", 0), counts.get("failed", 0)
    skipped = counts.get("skipped", 0)
    env = gates.get("environment", {})
    age_h = (time.time() - gates.get("finished_at", 0)) / 3600.0

    g1, g2 = st.columns([1.0, 1.6], gap="medium")
    with g1:
        panel("Suite",
              f'<div class="big">{passed} / {passed + failed}</div>'
              + rows((("failed", str(failed)), ("skipped", str(skipped)),
                      ("exit code", str(gates.get("exit"))),
                      ("duration", f'{gates.get("duration_s", 0) / 60:.1f} min'),
                      ("git", (env.get("git_sha") or "?")[:8]
                       + (" dirty" if env.get("git_dirty") else "")))),
              f'recorded {age_h:.1f} h ago by <code>tests/run_gates.py</code>')
        if env.get("git_dirty"):
            banner("warn", "<strong>The tree was dirty when these gates ran.</strong> "
                           "They describe uncommitted code, so the commit named "
                           "above does not reproduce them.")
        if age_h > 24:
            banner("warn", f"<strong>These results are {age_h / 24:.0f} days "
                           f"old.</strong> They describe the code as it was "
                           f"then.")
    with g2:
        if failed:
            banner("bad", f"<strong>{failed} gate(s) failing.</strong>")
            for f in gates.get("failures", []):
                with st.expander(f["nodeid"]):
                    st.code(f["message"])
        else:
            banner("good", "<strong>Every gate in the suite passes.</strong> "
                           "That is a statement about the gates that exist — "
                           "see the open items below for what has no gate.")
        mods = gates.get("modules", {})
        if mods:
            st.markdown(
                '<div class="panel"><span class="lbl">By module</span>'
                + "".join(
                    f'<div class="kv"><span>{os.path.basename(k)}</span>'
                    f'<span>{v.get("passed", 0)} passed'
                    + (f', <strong>{v["failed"]} failed</strong>'
                       if v.get("failed") else "") + '</span></div>'
                    for k, v in sorted(mods.items()))
                + '</div>', unsafe_allow_html=True)


# --------------------------------------------------------------------------
# ZF22 Table 3 -- the external benchmark
# --------------------------------------------------------------------------
st.markdown("### ZF22 Table 3 — ten published cases")

z = zf22_comparison()
if not z["rows"]:
    unavailable("ZF22 comparison",
                "No <code>output/zf22_table3.md</code>. Produce it with"
                "<br><code>PYTHONPATH=. python scripts/zf22_table3.py</code>")
else:
    col = {name: i for i, name in enumerate(z["header"])}
    z1, z2 = st.columns([1.0, 2.4], gap="medium")
    with z1:
        panel("Reproduced",
              f'<div class="big">{z["n_ok"]} / {z["n_total"]}</div>',
              f'η<sub>E</sub> within {z["band"]:.0%} of ZF22\'s published '
              f'value. They quote it to two figures, so this band is the '
              f'resolution of the comparison, not a tolerance picked to pass — '
              f'at 10% the failing case misses by 45%, and no defensible band '
              f'rescues it.')
        if z["failing"]:
            banner("bad",
                   f'<strong>Case {", ".join(z["failing"])} is not '
                   f'reproduced.</strong> This is the repository\'s headline '
                   f'open item, <strong>BLK-6</strong>. It is the only '
                   f'Muskat-unstable case of the ten — BENCH-10 gives the '
                   f'analytic reason — and neither mesh direction converges '
                   f'for it, so BENCH-09\'s agreement on the diagonal was a '
                   f'cancellation, not a check. The mockup showed 10/10.')
        else:
            banner("good", "<strong>All ten reproduce.</strong>")
    with z2:
        st.markdown(
            '<div class="panel"><span class="lbl">Case by case</span>'
            '<table class="zf"><tr>'
            '<th>case</th><th>b (BF25)</th><th>m</th>'
            '<th>η<sub>E</sub></th><th>ZF22</th><th>Δ</th></tr>'
            + "".join(
                '<tr class="{cls}"><td>{case}</td><td>{b}</td><td>{m}</td>'
                '<td>{mine:.3f}</td><td>{theirs:.3f}</td>'
                '<td>{d:+.1%}</td></tr>'.format(
                    cls="" if r["ok"] else "miss", case=r["case"],
                    b=r["cells"][col["b (BF25)"]], m=r["cells"][col["m"]],
                    mine=r["eta"], theirs=r["zf22_eta"], d=r["delta"])
                for r in z["rows"])
            + '</table></div>', unsafe_allow_html=True)

    with st.expander("The full table, as the script wrote it"):
        st.markdown("\n".join(z["lines"]))


# --------------------------------------------------------------------------
# open items
# --------------------------------------------------------------------------
st.markdown("### Open items")
blocked = read_blocked()
if not blocked:
    banner("warn", "<code>BLOCKED.md</code> could not be read.")
else:
    missing = [b for b in blocked if "missing" in b["kind"]]
    ambiguous = [b for b in blocked if "ambiguity" in b["kind"]]
    o1, o2 = st.columns(2, gap="medium")
    with o1:
        st.markdown(
            '<div class="panel"><span class="lbl">Missing physical data</span>'
            + "".join(f'<div class="kv"><span>{b["id"]}</span>'
                      f'<span style="text-align:right">{b["title"]}</span></div>'
                      for b in missing)
            + '<span class="note">A provisional value is in use for each so '
              'that work continues. None of them is a measurement, and BLK-1 '
              '— the centralizer record — moves the narrow-side velocity '
              'fraction by a factor of 130 on its own.</span></div>',
            unsafe_allow_html=True)
    with o2:
        st.markdown(
            '<div class="panel"><span class="lbl">Irreducible ambiguity</span>'
            + "".join(f'<div class="kv"><span>{b["id"]}</span>'
                      f'<span style="text-align:right">{b["title"]}</span></div>'
                      for b in ambiguous)
            + '<span class="note">The papers do not determine these. The most '
              'defensible option is implemented, the alternative is behind a '
              'flag, and the difference is quantified in '
              '<code>BLOCKED.md</code>.</span></div>',
            unsafe_allow_html=True)

st.markdown("### What has no gate at all")
banner("warn",
       "<strong>Read this next to the pass count.</strong> The suite tests what "
       "it tests. There is no K-GEP-1 <em>measurement</em> to validate against "
       "— no caliper-matched cement bond log, no returns record — so nothing "
       "on this page says the model predicts this well correctly. It says the "
       "model reproduces published benchmarks, conserves what it should, and "
       "fails where BLK-6 says it fails. The B-1 angle assumption is measured "
       "per run rather than gated, because its size depends on the fluid pair; "
       "on the production axes the shipped pair measures 9.8%, which the "
       "solver's own grading calls SIGNIFICANT.")
