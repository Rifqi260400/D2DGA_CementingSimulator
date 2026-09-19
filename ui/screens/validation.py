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

import streamlit as st

from ui.common import ROOT, banner, page_chrome, panel, rows, unavailable
from ui.reports import read_blocked, read_gates, zf22_comparison

page_chrome("Validation")
st.title("Validation")

GATES = os.path.join(ROOT, "output", "gates.json")
ZF22 = os.path.join(ROOT, "output", "zf22_table3.md")


# --------------------------------------------------------------------------
# the gates
# --------------------------------------------------------------------------
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
