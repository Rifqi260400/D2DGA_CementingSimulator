"""Screen 1 — Case setup and launch.

Everything settable, and everything the settings imply, on one page.

What this screen must not do, and the mechanism that stops it:

  * **Restate a parameter.** It does not know a single default. `launch.parameters()`
    introspects `kgep1_run.build_parser()`; adding an argument to the runner makes
    it appear here with no edit.
  * **Eyeball the flow classification.** The stability verdict is BF25 §3.3's
    Muskat criterion, computed from this case's own closures -- a real
    prediction, not a lookup. It costs a small closure-table build, so it is
    behind an explicit control that says what it costs, rather than being
    silently skipped or silently slow.
  * **Default a missing input.** The items in `BLOCKED.md` are shown as blocked,
    with the provisional value named as provisional.
  * **Launch something a terminal could not.** `launch.argv_from` builds the
    command, the command is shown before it runs, and the run is a subprocess.
"""

from __future__ import annotations

import numpy as np
import streamlit as st

from d2dga import muskat
from d2dga.config import Config, FluidsConfig, GridConfig, StandoffConfig
from d2dga.gapscale.tables import ClosureTable
from d2dga.scaling import Scaling
from ui import figures, launch, theme
from ui.common import (ROOT, banner, fluid_disclaimer, kv, page_chrome,
                       panel, rows)
from ui.reports import read_blocked

page_chrome("Case setup")
st.title("Case setup")

PARAMS = {p["dest"]: p for p in launch.parameters()}
GROUPS = (
    ("Case", ("wall", "inflow", "w0", "volumes", "eccentricity")),
    ("Mesh and numerics", ("n_phi", "n_xi", "cfl", "n_c", "n_h")),
    ("Displaced fluid (mud)", ("mud_density", "mud_consistency",
                               "mud_power_law_index", "mud_yield_stress")),
    ("Displacing fluid (cement)", ("cement_density", "cement_consistency",
                                   "cement_power_law_index",
                                   "cement_yield_stress")),
)


def _widget(p):
    """One input, entirely described by the parser."""
    label = p["dest"].replace("_", " ")
    if p["unit"]:
        label += f'  [{p["unit"]}]'
    key = "setup_" + p["dest"]
    if p["choices"]:
        return st.selectbox(label, p["choices"],
                            index=p["choices"].index(p["default"]),
                            key=key, help=p["help"])
    if p["type"] is int:
        return st.number_input(label, value=int(p["default"]), step=1,
                               min_value=int(p["lo"]) if p["lo"] else None,
                               max_value=int(p["hi"]) if p["hi"] else None,
                               key=key, help=p["help"])
    return st.number_input(label, value=float(p["default"]), format="%g",
                           min_value=float(p["lo"]) if p["lo"] is not None else None,
                           max_value=float(p["hi"]) if p["hi"] is not None else None,
                           key=key, help=p["help"])


# --------------------------------------------------------------------------
# the form
# --------------------------------------------------------------------------
values = {}
cols = st.columns(4, gap="medium")
for col, (title, dests) in zip(cols, GROUPS):
    with col:
        st.markdown(f'<span class="lbl">{title}</span>', unsafe_allow_html=True)
        for d in dests:
            values[d] = _widget(PARAMS[d])

st.markdown("---")


# --------------------------------------------------------------------------
# what the settings imply -- derived live, never typed in
# --------------------------------------------------------------------------
@st.cache_resource(show_spinner=False)
def _case(sig):
    """Geometry and scaling for a set of values. Cheap: no closure table."""
    vals = dict(sig)
    import sys, os
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import kgep1_run as runner
    cfg = Config(grid=GridConfig(int(vals["n_phi"]), int(vals["n_xi"])),
                 fluids=FluidsConfig(**{k: vals[k] for k in (
                     "mud_density", "mud_consistency", "mud_power_law_index",
                     "mud_yield_stress", "cement_density", "cement_consistency",
                     "cement_power_law_index", "cement_yield_stress")}),
                 standoff=StandoffConfig(float(vals["eccentricity"])))
    geo = runner.make_geometry(cfg, vals["wall"])
    mud, cement = cfg.fluids.as_fluids()
    sc = Scaling(mud, cement, r_a_hat_star=geo.r_a_hat_star,
                 delta_star=geo.delta_star, mean_velocity=float(vals["w0"]))
    return cfg, geo, sc


sig = tuple(sorted(values.items()))
try:
    cfg, geo, sc = _case(sig)
except (ValueError, Exception) as exc:            # noqa: BLE001 - shown, not hidden
    banner("bad", f"<strong>This case cannot be built.</strong> <code>{exc}</code>")
    st.stop()

g = geo.grid
dpi = np.asarray(geo.narrow_gap_parameter(g.xi_centres), dtype=float)
outside = float(np.mean(dpi > theme.DELTA_PI_VALIDATED))

left, right = st.columns([1.0, 1.0], gap="medium")
with left:
    panel("Dimensionless groups",
          rows((("b", f"{sc.buoyancy_number:.4g}"), ("m", f"{sc.m:.4g}"),
                ("B", f"{sc.B:.4g}"), ("Fr*", f"{sc.Fr_star:.4g}"),
                ("Δρ&#770;", f"{sc.delta_rho:.4g}"),
                ("Z", f"{g.Z:.1f}"), ("Δξ", f"{g.dxi:.3f}"))),
          "computed from these inputs by <code>d2dga.scaling.Scaling</code>, "
          "the same object the run uses")
    panel("Reynolds number",
          kv("Re", '<span style="color:#9A9A92">not used</span>'),
          "D2DGA is non-inertial: BF25 (2.16) is a constraint with no inertia "
          "term, so Re enters nothing. Shown because its absence is "
          "informative — ZF22 cases 9 and 10 exist to test exactly that.")
with right:
    fig = figures.envelope_only(geo)
    st.pyplot(fig, use_container_width=True)
    kind = "warn" if outside else "good"
    banner(kind,
           f"<strong>{outside:.0%} of the interval is above ZF23's validated "
           f"δ/π = {theme.DELTA_PI_VALIDATED}</strong>, peaking at "
           f"<code>{dpi.max():.4f}</code> "
           f"({dpi.max() / theme.DELTA_PI_VALIDATED:.1f}×). Above it the "
           f"Hele-Shaw reduction is a prediction outside the regime it was "
           f"demonstrated in.")

fluid_disclaimer(cfg.fluids, None)


# --------------------------------------------------------------------------
# a-priori flow classification -- BF25 (3.13), computed
# --------------------------------------------------------------------------
st.markdown("### Predicted regime, before running")

st.markdown(
    '<span class="note">BF25 §3.3\'s Muskat analysis is the only closed-form '
    'prediction in the six papers that applies to a Herschel–Bulkley pair, and '
    'it takes exactly these inputs. A thin finger of pure cement feels the same '
    'axial pressure gradient as the dispersing front (3.11); comparing its '
    'speed (3.12) with the front\'s own leading wave gives Δw(0⁺). Positive: '
    'the finger outruns the front and penetrates. Negative: it is absorbed. '
    'Validated against BF25\'s own M₃<sup>min</sup> = 3/2 for a Newtonian pair '
    '— <code>tests/test_bench10_muskat.py</code>.</span>',
    unsafe_allow_html=True)

want = st.checkbox(
    "Compute it (builds a small closure table for this pair — seconds for a "
    "Newtonian mud, up to a minute with two yield stresses)", value=False,
    key="setup_muskat")


@st.cache_resource(show_spinner="Solving the gap-scale problem…")
def _muskat(fluid_sig, b):
    cfg2, geo2, sc2 = _case(fluid_sig)
    tab = ClosureTable(sc2.scaled_fluid1, sc2.scaled_fluid2,
                       c_grid=np.linspace(0.0, 1.0, 21),
                       h_grid=np.array([1.0]),
                       umag_grid=np.array([1.0]),
                       gb_grid=np.array([b]),
                       n_y=160, tol=1e-9).build()
    # BF25 §3.1 evaluates the closures "imposing (v_bar, w_bar) = (0, 1)", so
    # H = 1 and |u_bar| = 1 are the axes the analysis is written on.
    def at(i):
        return lambda c: tab(c, H=1.0, umag=1.0, gb=b)[i]
    I1, I2, q0, I3 = (at(i) for i in range(4))
    dw0 = muskat.dw_at_leading_edge(I1, I2, q0, I3, b)
    regime, dw, c0 = muskat.classify(I1, I2, q0, I3, b, n=101, c0_max=0.9)
    return dw0, regime, dw, c0, tab.assumption_report(), tab.n_failed


if want:
    dw0, regime, dw, c0, report, n_failed = _muskat(sig, sc.buoyancy_number)
    m1, m2 = st.columns([1.0, 1.4], gap="medium")
    with m1:
        verdict = ("finger PENETRATES the front" if dw0 > 0
                   else "finger is ABSORBED by the front")
        panel("BF25 (3.13)",
              f'<div class="big">{dw0:+.4g}</div>'
              + rows((("Δw(0⁺)", f"{dw0:+.4g}"),
                      ("verdict", verdict),
                      ("regime over c₀ ∈ (0, 0.9]", regime),
                      ("unconverged gap solves", str(n_failed)))),
              "Δw(0⁺) is the discriminating half of (3.13); the c₀ → 1 limit is "
              "degenerate as printed — BLK-7.")
        if n_failed:
            banner("warn", f"<strong>{n_failed} of the gap-scale solves did not "
                           f"converge.</strong> The closures at those "
                           f"concentrations are the solver's last iterate, so "
                           f"this verdict is weaker than it looks.")
    with m2:
        st.pyplot(figures.muskat(c0, dw), use_container_width=True)
    fluid_disclaimer(
        cfg.fluids, report,
        measured_on=" <em>on this verdict's own probe table</em>")
    banner("warn",
           "<strong>That percentage is not the run's.</strong> The table just "
           "built for this verdict has a single velocity node at |u&#772;| = 1, "
           "because BF25 §3.1 evaluates the closures at "
           "(v&#772;, w&#772;) = (0, 1). The table a <em>run</em> builds has a "
           "velocity axis reaching 3000 — it must, because b·𝓘₁ ≈ 1500 on this "
           "pair drives |u&#772;| past 370 — and the angle sensitivity grows "
           "with |u&#772;|. On the shipped pair the run's table measures "
           "<strong>9.8%</strong> where this one measures a fraction of a "
           "percent. Read this one as applying to the verdict above it, and "
           "nothing else.")
else:
    st.markdown('<div class="unavail"><strong>Not computed</strong><br>'
                'It needs the gap-scale closures for this pair, which means '
                'solving the augmented-Lagrangian problem at 21 concentrations. '
                'Nothing is shown in its place.</div>', unsafe_allow_html=True)


# --------------------------------------------------------------------------
# what is blocked -- shown as blocked, not silently defaulted
# --------------------------------------------------------------------------
st.markdown("### Inputs that are still missing")
blocked = read_blocked()
if blocked:
    st.markdown(
        '<span class="note">From <code>BLOCKED.md</code>, read at page load. '
        'A provisional value is in use for each so that work continues; none of '
        'them is a measurement.</span>', unsafe_allow_html=True)
    bc = st.columns(min(3, len(blocked)), gap="medium")
    for i, item in enumerate(blocked):
        with bc[i % len(bc)]:
            panel(item["id"], f'<div style="font-size:12px">{item["title"]}</div>',
                  item["kind"])
else:
    banner("warn", "<code>BLOCKED.md</code> could not be read, so this page "
                   "cannot tell you what is missing. That is itself a gap.")


# --------------------------------------------------------------------------
# launch
# --------------------------------------------------------------------------
st.markdown("### Launch")
argv = launch.argv_from(values)
st.code(launch.shell_command(argv), language="bash")
st.markdown('<span class="note">Every parameter is written out, including the '
            'ones left at their default, so the line keeps reproducing this run '
            'after the defaults move. Pressing Run executes exactly this, as an '
            'ordinary subprocess.</span>', unsafe_allow_html=True)

c1, c2 = st.columns([1, 4])
with c1:
    if st.button("Run", type="primary", key="setup_launch"):
        handle = launch.launch(values)
        st.session_state["active_run"] = handle
        st.success(f"Started pid {handle['pid']}. Open **Run monitor**.")
with c2:
    st.markdown('<span class="note">The run detaches, so it survives a browser '
                'refresh and a restart of this server. There is no Pause: the '
                'solver has no such state and a fake one would misdescribe the '
                'button. Stop is safe — the run checkpoints every 180 s and the '
                'same command resumes it.</span>', unsafe_allow_html=True)
