"""Presentation helpers shared by every screen.

Nothing here computes physics or post-processing.  If a screen needs a number,
it calls `d2dga`; this module only decides how it looks on the page.

One Streamlit quirk shapes the whole file: **Streamlit closes every markdown
block it is handed.**  Writing `st.markdown('<div class="panel">')` and the
closing tag in a second call produces an empty white box and orphaned content,
which is how the first version of the Results screen shipped a row of blank
panels.  So a panel is always emitted in ONE call, and these helpers return
HTML strings rather than writing them.
"""

from __future__ import annotations

import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st                                 # noqa: E402

from ui import theme                                   # noqa: E402


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def kv(label, value) -> str:
    """One label/value row, as HTML."""
    return f'<div class="kv"><span>{label}</span><span>{value}</span></div>'


def rows(pairs) -> str:
    return "".join(kv(a, b) for a, b in pairs)


def panel(title: str, body: str, note: str = "") -> None:
    st.markdown(
        f'<div class="panel"><span class="lbl">{title}</span>{body}'
        + (f'<span class="note">{note}</span>' if note else "")
        + '</div>', unsafe_allow_html=True)


def unavailable(title, reason) -> None:
    """For a quantity the run does not contain.

    Distinct from a warning: this is an ABSENCE, and it always says why, so the
    reader can tell "the solver does not produce this" from "this run happens
    not to have it".  Never replaced by a plausible value.
    """
    st.markdown(f'<div class="unavail"><strong>{title} — unavailable</strong>'
                f'<br>{reason}</div>', unsafe_allow_html=True)


def banner(kind: str, html: str) -> None:
    """`good` | `warn` | `bad` — increasing strength.

    `warn`: outside a validated range, read with care.
    `bad`: the numbers refer to a different problem than the one on screen.
    """
    st.markdown(f'<div class="{kind}">{html}</div>', unsafe_allow_html=True)


def download_figure(fig, name) -> None:
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    st.download_button("Export figure (PNG)", buf.getvalue(), file_name=name,
                       mime="image/png", key=f"dl_{name}")


def page_chrome(_unused: str = "") -> None:
    """The stylesheet and the sidebar masthead, on every screen.

    The masthead names the WELL, not the screen: `st.navigation` already shows
    which screen is open, and repeating it was the kind of duplication that
    makes a reader check which one is authoritative.
    """
    st.markdown(theme.CSS, unsafe_allow_html=True)
    st.sidebar.markdown(
        f'<div style="font-weight:600;font-size:15px">D2DGA</div>'
        f'<div style="font-size:12px;color:{theme.MUTED}">K-GEP-1</div>',
        unsafe_allow_html=True)
    st.sidebar.markdown("---")


def fluid_disclaimer(fluids, assumption_report: str | None,
                     recorded: bool = True, measured_on: str = "") -> None:
    """The B-1 / validated-pair disclaimer, identical on every screen.

    Two statements, deliberately not merged.

    PROVENANCE -- is this the fluid pair every published figure was computed
    with?  `FluidsConfig.departures_from_validated()` answers from the config's
    own fields and names each one, because "not the validated pair" alone does
    not tell a reader whether to distrust the efficiency or the breakthrough
    time.

    MEASUREMENT -- how far the closure table's theta = 0 slice is from the truth
    for THIS pair.  Only the table knows, and only when it is BUILT rather than
    loaded, so it is passed in from the run's metrics and shown as a warning --
    never a green tick, never a guess -- when it is absent.

    Only the second licenses a result, and the two can disagree: the validated
    pair is not a benign pair.  On the production velocity axis it measures
    9.8%, which the solver's own grading calls SIGNIFICANT.
    """
    dep = fluids.departures_from_validated()
    if dep:
        banner("bad",
               '<strong>Not the validated fluid pair.</strong> Every figure in '
               '<code>AUDIT_REPORT.md</code>, <code>docs/gate_status.md</code> '
               'and <code>output/kgep1_results.md</code> was computed with the '
               'Materials 2025 Table 1 pair. This differs in: <code>'
               + '</code>, <code>'.join(dep) +
               '</code>. Those figures do not describe it.')
    if assumption_report:
        worst = assumption_report.splitlines()[0]
        kind = ("good" if "[benign]" in worst
                else "bad" if "LARGE" in worst else "warn")
        banner(kind,
               '<strong>Closure-table angle assumption (B-1), measured for this '
               f'pair{measured_on}:</strong><br><code>{worst}</code><br>'
               'Every table entry is built with u&#772; parallel to '
               'G&#771;<sub>b</sub> (θ = 0); the real flow has θ ≠ 0 wherever '
               'the front is tilted. The figure is an <em>upper bound</em> on '
               'the angle error — θ = 0 against θ = 90° — not the error in '
               'η<sub>E</sub>. The fix is a fifth table axis (NUM-14), not a '
               'correction factor.')
    else:
        banner("warn",
               '<strong>Closure-table angle assumption (B-1): not measured'
               + ('' if recorded else ' for this run') +
               '.</strong> The θ = 0 error is a property of the fluid pair '
               '<em>and of the table\'s velocity axis</em>, and is only known '
               'when the table is <em>built</em>. It is not guessed here. '
               'Rebuilt on the production axes the shipped water/cement pair '
               'measures <strong>9.8% — SIGNIFICANT, not benign</strong>, '
               'peaking at |u&#772;| = 3000 where the front is most tilted. '
               'That had never been measured on those axes before, because the '
               'production runs loaded a cached table and a loaded table '
               'reports "not measured".')
