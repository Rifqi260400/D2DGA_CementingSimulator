"""
D2DGA solver UI.

Build order (spec section 7) step 1: the **Results viewer, read-only**.  It
points at a run that already exists on disk.  No launching, no live monitoring,
no mocked numbers.

Run it with:   PYTHONPATH=. streamlit run ui/app.py

Design rules this file obeys and must keep obeying:
  * every number traces to a solver output -- see `ui/INVENTORY.md`;
  * the delta/pi validity envelope is drawn on the SAME depth axis as the
    concentration field, aligned row for row (spec section 3);
  * a quantity the run did not record is shown as UNAVAILABLE with the reason,
    never substituted;
  * physics and post-processing are imported from `d2dga`, never reimplemented.
"""

from __future__ import annotations

import io
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                    # noqa: E402
import streamlit as st                                # noqa: E402

from d2dga.postprocess import (displacement_efficiency,   # noqa: E402
                               narrow_side_profile, residual_fraction)
from d2dga.simulation import RunResult                 # noqa: E402
from ui import figures, runs, theme                    # noqa: E402

st.set_page_config(page_title="D2DGA — K-GEP-1", layout="wide",
                   initial_sidebar_state="expanded")
st.markdown(theme.CSS, unsafe_allow_html=True)


# --------------------------------------------------------------------------
# small helpers -- presentation only
# --------------------------------------------------------------------------
def kv(label, value) -> str:
    """One label/value row.  Returns HTML; panels are emitted in ONE markdown
    call, because Streamlit closes every block it is given -- an opening
    `<div class="panel">` on its own renders as an empty white box."""
    return f'<div class="kv"><span>{label}</span><span>{value}</span></div>'


def panel(title: str, body: str, note: str = "") -> None:
    st.markdown(
        f'<div class="panel"><span class="lbl">{title}</span>{body}'
        + (f'<span class="note">{note}</span>' if note else "")
        + '</div>', unsafe_allow_html=True)


def unavailable(title, reason):
    st.markdown(f'<div class="unavail"><strong>{title} — unavailable</strong>'
                f'<br>{reason}</div>', unsafe_allow_html=True)


def download_figure(fig, name):
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    st.download_button("Export figure (PNG)", buf.getvalue(), file_name=name,
                       mime="image/png", key=f"dl_{name}")


# Half-order discretisation error of a threshold front, measured against the
# exact rarefaction at n_xi = 80 (A-2).  Scaled here as O(dxi^(1/2)); the same
# constants live in scripts/kgep1_run.py, which prints them on every run.
_TBR_ERR_AT_80 = {0.01: 0.175, 0.1: 0.069, 0.5: 0.033}


# --------------------------------------------------------------------------
# sidebar: pick a run
# --------------------------------------------------------------------------
st.sidebar.markdown("### D2DGA")
st.sidebar.markdown('<span class="mono" style="font-size:12px;color:#5C5C55">'
                    'K-GEP-1</span>', unsafe_allow_html=True)
st.sidebar.markdown("---")

output_dir = st.sidebar.text_input("Output directory", value="output")
all_runs = runs.discover(output_dir)

if not all_runs:
    st.title("Results")
    unavailable("No runs found",
                f"No <code>kgep1_ckpt_*.npz</code> in <code>{output_dir}</code>. "
                f"Produce one with<br><code>PYTHONPATH=. python "
                f"scripts/kgep1_run.py</code>")
    st.stop()

labels = {r.name: r for r in all_runs}
chosen = st.sidebar.radio("Run", list(labels), index=0,
                          format_func=lambda s: s.replace("_", " · "))
run = labels[chosen]

st.sidebar.markdown("---")
st.sidebar.markdown('<span class="lbl">Provenance</span>', unsafe_allow_html=True)
if run.config_source == runs.RECORDED:
    st.sidebar.markdown(
        f'<div class="good">Settings <strong>recorded</strong> with the run.'
        f'<br>tag <code>{run.settings_tag}</code>'
        + (f'<br>git <code>{run.git_sha}</code>'
           + (' <em>dirty</em>' if run.git_dirty else '') if run.git_sha else '')
        + '</div>', unsafe_allow_html=True)
else:
    st.sidebar.markdown(
        '<div class="warn"><strong>Settings reconstructed, not recorded.</strong> '
        'This run predates the provenance writer, so only the six settings in '
        'its filename are known. Everything else — fluids, eccentricity, CFL, '
        'table resolution — is taken from the <em>current</em> '
        '<code>config.py</code> and is right only if that file has not changed '
        'since. Re-run it to get a recorded config.</div>',
        unsafe_allow_html=True)


# --------------------------------------------------------------------------
# rebuild what the run used
# --------------------------------------------------------------------------
@st.cache_resource(show_spinner="Rebuilding geometry…")
def _geometry(checkpoint: str, _mtime: float):
    r = runs.load(checkpoint)
    geo = runs.build_geometry(r)
    return geo, runs.build_scaling(r, geo)


geo, sc = _geometry(run.checkpoint, run.mtime)
g = geo.grid
c = run.concentration
times, effs, outlet = run.history

# `breakthrough_at` is the solver's own interpolation.  Reconstructing a
# RunResult around the recorded history calls it rather than copying it.
result = RunResult(times=np.asarray(times), efficiency=np.asarray(effs),
                   outlet_concentration=np.asarray(outlet), concentration=c,
                   stream_function=np.zeros((g.n_phi + 1, g.n_xi + 1)),
                   reports=[], breakthrough_time=float("nan"),
                   conservation_error=run.conservation_error)

eta = displacement_efficiency(geo, c)


# --------------------------------------------------------------------------
# header
# --------------------------------------------------------------------------
left, right = st.columns([3, 2])
with left:
    st.title("Results")
    st.markdown(
        f'<span class="note">{run.steps:,} steps &middot; final state at '
        f'{run.args["volumes"]:g} pumped volumes &middot; '
        f'{run.args["wall"]} wall &middot; {run.args["inflow"]}</span>',
        unsafe_allow_html=True)
with right:
    st.markdown('<span class="lbl">Reproduce from the CLI</span>',
                unsafe_allow_html=True)
    st.code(run.command(), language="bash")

st.markdown("")


# --------------------------------------------------------------------------
# the mandated pairing: field beside its validity envelope
# --------------------------------------------------------------------------
col_field, col_num = st.columns([1.05, 1.0], gap="medium")

with col_field:
    span = float(c.max() - c.min())
    stretch = False
    if span < 0.25:
        stretch = st.checkbox(
            f"Stretch the colour scale to this run's range "
            f"({c.min():.4f} – {c.max():.4f})",
            value=True,
            help="The absolute 0–1 ramp is the honest default, but this run is "
                 "almost fully displaced, so on that scale the field is a flat "
                 "block and the narrow-side deficit is invisible. Stretching "
                 "makes the structure visible; the colour bar states the range "
                 "either way, so no reading is implied that the data does not "
                 "support.")
    fig = figures.field_with_envelope(geo, c, stretch=stretch)
    st.pyplot(fig, use_container_width=True)
    dpi = np.asarray(geo.narrow_gap_parameter(g.xi_centres), dtype=float)
    st.markdown(
        f'<div class="warn">The right-hand panel is the narrow-gap parameter on '
        f'the same depth axis. <strong>{np.mean(dpi > theme.DELTA_PI_VALIDATED):.0%} '
        f'of this interval is above ZF23\'s validated 0.038</strong>, peaking at '
        f'<code>{dpi.max():.4f}</code> — {dpi.max() / theme.DELTA_PI_VALIDATED:.1f}× '
        f'it. Where the bar is amber the Hele-Shaw reduction is being used '
        f'outside the regime it was demonstrated in.</div>',
        unsafe_allow_html=True)
    download_figure(fig, f"{run.name}_field.png")

with col_num:
    panel("Displacement efficiency",
          f'<div class="big">{eta:.4f}</div>',
          f'volume-weighted, whole interval, at {run.args["volumes"]:g} volumes')

    scale = (80.0 / g.n_xi) ** 0.5
    rows = ""
    for th, rel in _TBR_ERR_AT_80.items():
        tb = result.breakthrough_at(th) / g.Z
        shown = "—" if not np.isfinite(tb) else f"{tb:.3f}"
        rows += kv(f"c&#772; &gt; {th}",
                   f'{shown} <span style="color:{theme.AMBER_INK}">'
                   f'±{rel * scale:.0%}</span>')
    panel("Breakthrough t<sub>br</sub>/Z", rows,
          "A threshold front is a level set in a diffusive tail, so its error "
          "is O(Δξ<sup>1/2</sup>) — half order — while η<sub>E</sub> is an "
          "integral and is O(Δξ). Halving the 0.01 error costs four times the "
          "cells.")

    panel("Invariants and residuals",
          kv("Volume conservation", f"{run.conservation_error:.2e}")
          + kv("c&#772; range", f"{c.min():.4f} – {c.max():.4f}")
          + kv("Narrow-side minimum",
               f"{float(np.min(narrow_side_profile(geo, c))):.4f}")
          + kv("Residual c&#772; &lt; 0.5", f"{residual_fraction(geo, c):.4f}")
          + kv("Picard at cap", f"{run.picard_unconverged}"))

    panel("Dimensionless groups",
          "".join(kv(lab, val) for lab, val in (
              ("b", f"{sc.buoyancy_number:.4g}"), ("m", f"{sc.m:.4g}"),
              ("B", f"{sc.B:.4g}"), ("Fr*", f"{sc.Fr_star:.4g}"),
              ("Δρ&#770;", f"{sc.delta_rho:.4g}"), ("Z", f"{g.Z:.1f}"))))


# --------------------------------------------------------------------------
# history and profiles
# --------------------------------------------------------------------------
st.markdown("### History")
h1, h2 = st.columns(2, gap="medium")
with h1:
    fig_h = figures.history(times, effs, outlet, g.Z,
                           t_br=result.breakthrough_at(0.01))
    st.pyplot(fig_h, use_container_width=True)
    st.markdown(f'<span class="note">{len(times):,} samples — recorded every '
                f'step, so t<sub>br</sub> does not depend on the logging '
                f'frequency (R-1).</span>', unsafe_allow_html=True)
    download_figure(fig_h, f"{run.name}_history.png")
with h2:
    fig_p = figures.profiles(geo, c)
    st.pyplot(fig_p, use_container_width=True)
    st.markdown('<span class="note">Both from <code>postprocess</code>: the '
                'volume-weighted azimuthal mean, and the narrow-side column '
                'where mud survives longest.</span>', unsafe_allow_html=True)
    download_figure(fig_p, f"{run.name}_profiles.png")


# --------------------------------------------------------------------------
# front classification -- from the run's own metrics, never recomputed here
# --------------------------------------------------------------------------
st.markdown("### Front classification")
zf = (run.metrics or {}).get("zf23")
if zf:
    z1, z2 = st.columns([2, 1], gap="medium")
    with z1:
        panel("ZF23 (3.6)",
              kv("σ<sub>w+r</sub>", f'{zf["sigma_plus"]:.4f} &nbsp;/&nbsp; 0.08')
              + kv("|w&#772;<sub>r</sub>&#8314;|",
                   f'{zf["area_plus"]:.4f} &nbsp;/&nbsp; 0.05')
              + kv("|w&#772;<sub>r</sub>&#8315;|", f'{zf["area_minus"]:.4f}')
              + kv("sampled at t/Z", f'{zf["sampled_at_t_over_Z"]:.3f}'),
              "ZF23 (3.6) requires <strong>both</strong> of the first two above "
              "threshold. σ<sub>w−r</sub> is deliberately omitted: it moves "
              "14 % non-monotonically over a 32× range of bin counts, so it is "
              "not a converged quantity and no conclusion rests on it (NUM-25).")
    with z2:
        verdict = "Dispersive" if zf["is_dispersive"] else "Not dispersive"
        css = "warn" if zf["is_dispersive"] else "good"
        st.markdown(f'<div class="{css}"><span class="lbl">Verdict</span><br>'
                    f'<strong style="font-size:15px">{verdict}</strong></div>',
                    unsafe_allow_html=True)
else:
    unavailable(
        "ZF23 dispersion metrics",
        "They are computed <em>during</em> a run, from a pre-breakthrough "
        "snapshot while the whole front is still inside the domain. Only the "
        "<strong>final</strong> field is saved, and by then the front has left, "
        "so recomputing them here would produce a number with no meaning. This "
        "run has no <code>.metrics.json</code>; re-run it to get one.")


# --------------------------------------------------------------------------
# what this run does not contain
# --------------------------------------------------------------------------
st.markdown("### Not available for this run")
n1, n2, n3 = st.columns(3, gap="medium")
with n1:
    unavailable("Time snapshots", "A run saves only the final field, so there "
                "is no time slider and no similarity collapse.")
with n2:
    unavailable("Stream function Ψ", "Not persisted, so velocity and streamline "
                "plots cannot be drawn from a finished run.")
with n3:
    unavailable("η<sub>N</sub>, Δw<sub>f</sub>",
                "Neither exists in <code>postprocess.py</code>. They would have "
                "to be added to the solver, not computed here.")

with st.expander("Recorded configuration"):
    if run.payload:
        st.json(run.payload, expanded=False)
        st.download_button("Download config.json",
                           json.dumps(run.payload, indent=1),
                           file_name=f"{run.name}.config.json",
                           mime="application/json")
    else:
        st.markdown(
            '<div class="warn">No configuration was recorded with this run. '
            'What is shown elsewhere on this page is reconstructed from the '
            'filename plus the current <code>config.py</code>.</div>',
            unsafe_allow_html=True)
