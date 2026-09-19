"""Screen 3 — Results.

Read-only. It points at a run that already exists on disk and shows what that
run recorded — nothing more. A quantity the run did not record is marked
unavailable with the reason, never substituted.

Design rules this file obeys and must keep obeying:
  * every number traces to a solver output — see `ui/INVENTORY.md`;
  * the δ/π validity envelope is drawn on the SAME depth axis as the
    concentration field, aligned row for row, in ONE figure with a shared y —
    so the alignment is structural and cannot drift (spec §3);
  * physics and post-processing are imported from `d2dga`, never reimplemented;
    even `breakthrough_at` is called on a reconstructed `RunResult` rather than
    copied.
"""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np                                    # noqa: E402
import streamlit as st                                # noqa: E402

from d2dga.postprocess import (displacement_efficiency,   # noqa: E402
                               narrow_side_efficiency,
                               narrow_side_profile, residual_fraction,
                               tbr_relative_error)
from d2dga.simulation import RunResult                 # noqa: E402
from ui import figures, runs, theme                    # noqa: E402
from ui.common import (download_figure, fluid_disclaimer, kv,  # noqa: E402
                       page_chrome, panel, unavailable)

page_chrome("Results")


# --------------------------------------------------------------------------
# sidebar: pick a run
# --------------------------------------------------------------------------
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
# fluid pair and the B-1 disclaimer
#
# Placed here, above every number, because it governs how all of them may be
# read.  Two separate statements, deliberately not merged:
#
#   * PROVENANCE -- is this the fluid pair every published number in the repo
#     was computed with?  `FluidsConfig.departures_from_validated` answers that
#     from the config's own fields.
#   * MEASUREMENT -- how far the closure table's theta = 0 slice is from the
#     truth for THIS pair.  Only the table knows, and only if it was built
#     rather than loaded, so it is read from the run's metrics.json and shown as
#     unavailable when the run did not record it.  It is never estimated here.
#
# The second is what actually licenses a result, which is why a run with no
# measurement gets a warning and not a green tick.
# --------------------------------------------------------------------------
_cfg_fluids = runs.build_config(run).fluids
_ang = (run.metrics or {}).get("closure_table_assumption_report")

fl_a, fl_b = st.columns([1.0, 1.0], gap="medium")
with fl_a:
    mud_f, cem_f = _cfg_fluids.as_fluids()
    panel("Fluid pair",
          kv("Displaced", f"{mud_f.name}")
          + kv("&nbsp;&nbsp;ρ&#770;, κ&#770;, n, τ&#770;<sub>Y</sub>",
               f"{mud_f.density:g} &middot; {mud_f.consistency:g} &middot; "
               f"{mud_f.power_law_index:g} &middot; {mud_f.yield_stress:g}")
          + kv("Displacing", f"{cem_f.name}")
          + kv("&nbsp;&nbsp;ρ&#770;, κ&#770;, n, τ&#770;<sub>Y</sub>",
               f"{cem_f.density:g} &middot; {cem_f.consistency:g} &middot; "
               f"{cem_f.power_law_index:g} &middot; {cem_f.yield_stress:g}"),
          "kg/m³ &middot; Pa s<sup>n</sup> &middot; – &middot; Pa"
          + ("" if run.config_source == runs.RECORDED
             else " &mdash; from the current config.py, not recorded"))
with fl_b:
    fluid_disclaimer(_cfg_fluids, _ang,
                     recorded=run.config_source == runs.RECORDED)

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
    eta_n = narrow_side_efficiency(geo, c)
    panel("Displacement efficiency",
          f'<div class="big">{eta:.4f}</div>'
          + kv("η<sub>N</sub>, narrow quarter", f"{eta_n:.4f}")
          + kv("narrowest column, worst",
               f"{float(np.min(narrow_side_profile(geo, c))):.4f}"),
          f'volume-weighted, at {run.args["volumes"]:g} volumes. η<sub>E</sub> '
          f'is dominated by the wide side, which was never in doubt; '
          f'η<sub>N</sub> is the same integral over the narrow quarter, which '
          f'is what the job is about.')

    rows = ""
    any_br = False
    for th in (0.01, 0.1, 0.5):
        tb = result.breakthrough_at(th) / g.Z
        if np.isfinite(tb):
            any_br = True
            # the error constants are measured, and live with the other
            # post-processing so this screen and the runner cannot disagree
            value = (f'{tb:.3f} <span style="color:{theme.AMBER_INK}">'
                     f'±{tbr_relative_error(th, g.n_xi):.0%}</span>')
        else:
            # no error bar on a value that does not exist: an em-dash with a
            # tolerance beside it reads as a number that was measured badly
            value = '<span style="color:#9A9A92">not reached</span>'
        rows += kv(f"c&#772; &gt; {th}", value)
    if not any_br:
        rows += ('<span class="note">The outlet concentration never reached '
                 'these thresholds, so this run has no breakthrough time. '
                 'Nothing is extrapolated.</span>')
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
    unavailable("Δw<sub>f</sub>, the measured front-speed difference",
                "BF25 (3.13)'s Δw can be <em>predicted</em> before a run — it "
                "is on <strong>Case setup</strong> — but measuring it from a "
                "finished run needs the finger and the front tracked "
                "separately in time, and a run keeps only its final field. "
                "η<sub>N</sub> was in this panel until 2026-09-19 and is now "
                "computed, beside η<sub>E</sub>.")

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
