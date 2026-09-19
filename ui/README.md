# D2DGA solver UI

    pip install streamlit matplotlib
    PYTHONPATH=. streamlit run ui/app.py

Read `ui/INVENTORY.md` first: it maps every control and panel to a real solver
input or output, and lists what the mockup asked for that the solver does not
have.

## What exists today

**Screen 1 of 4 — Results, read-only.** It points at a run that already exists
in `output/` and shows what that run recorded. No launching, no live monitoring.
That is build order step 1 of the spec, and it is deliberate: it forces the
output inventory to be right before anything depends on it.

## Choices, and why

**Streamlit.** Pure Python, no build step, one researcher on one machine. The
live-monitoring screen will poll a status file written by the run's own
`on_step` callback, so it needs no server and no change to the solver's
numerics.

**Runs are plain subprocesses.** A run launched from the UI will be the same
process invocation as the CLI one, which is what makes "byte-identical to the
CLI" testable rather than asserted. The Results screen prints the command that
reproduces what you are looking at.

**Parameters are read, never restated.** `scripts/kgep1_run.build_parser()` is
the single definition of names, defaults, types, choices and help; `PARAM_META`
beside it adds the units and ranges argparse cannot express. The UI introspects
both. Nothing about a parameter is written down twice.

**Post-processing is imported.** `displacement_efficiency`, `axial_profile`,
`narrow_side_profile`, `residual_fraction`, `zf23_metrics` and
`RunResult.breakthrough_at` are called from `d2dga`, never reimplemented here.
`ui/` holds arrangement and nothing else.

**One colour system.** `ui/theme.py` defines the ground, the single teal accent
and the concentration ramp (tan mud → teal cement) once. No other file may
declare a second ramp.

## Rules this interface keeps

* **The δ/π validity envelope is drawn on the same depth axis as the
  concentration field, in one figure with a shared y-axis** — so the alignment
  is structural, not eyeballed. For K-GEP-1 the whole interval is above ZF23's
  validated 0.038, and you can see where the front is and where the model is
  untrustworthy at once.
* **The concentration field is unwrapped**, azimuth across with the wide side
  left, depth down, as every source paper plots it.
* **Nothing is mocked.** A quantity a run did not record is shown as
  *unavailable* with the reason. There are no placeholder numbers anywhere.
* **Provenance is stated, not assumed.** A run carries either a recorded config
  (`RECORDED`) or one reconstructed from its filename plus the current
  `config.py` (`RECONSTRUCTED`). The sidebar says which, and warns about the
  second.

## Honest limits of the Results screen

A finished run persists the **final** field and the scalar history — not
intermediate snapshots and not Ψ. So:

| Panel | Status |
|---|---|
| Time slider, similarity collapse | **Unavailable** — no snapshots are saved |
| Velocity / streamline plots | **Unavailable** — Ψ is not persisted |
| η_N, Δw_f | **Unavailable** — neither exists in `postprocess.py` |
| ZF23 dispersion metrics | Shown **only from `metrics.json`**. They are computed during a run from a pre-breakthrough snapshot; by the final state the front has left the domain, so recomputing them here would produce a meaningless number. A run predating `runio` has none. |

The colour scale defaults to the absolute 0–1 ramp. For a nearly complete
displacement that is a flat block — true, and useless — so a checkbox stretches
it to the run's own range. The colour bar states which range is in use either
way.

## The fluid-pair disclaimer

Above every number the screen states two things, kept separate on purpose.

**Provenance** — whether these are the fluids every published figure was
computed with. `FluidsConfig.departures_from_validated()` names each field that
differs, as `name: was -> is`, so a reader learns *which* property moved rather
than only that something did. A run on the shipped pair gets no banner; a run on
anything else gets a red one saying the repository's figures do not describe it.

**Measurement** — how far the closure table's `θ = 0` slice is from the truth
for this pair (finding B-1). Only the table knows, and only when it is *built*
rather than loaded, so the screen reads the run's own `assumption_report()` out
of `metrics.json`. A run that recorded none gets a **warning, not a green
tick**, and the number is never estimated.

The two are not merged because only the second licenses a result, and they can
disagree: the validated pair is not a benign pair. Measured on the production
axes it reaches **9.8% — SIGNIFICANT by the solver's own grading** — a figure
that had never been taken before 2026-09-19, because the production runs loaded
a cached table and a loaded table reports `not measured`. Two things follow, and
both are on screen: the size is governed by the viscosity ratio `m` rather than
by whether the mud is Newtonian (a Newtonian mud at an ordinary 5 mPa·s already
measures 11%), and the figure is an **upper bound** on the angle error —
`θ = 0` against `θ = 90°` — not the error in `η_E`. See `docs/assumptions.md`
FLU-07 and `docs/remediation_log.md` REM-11.

Since 2026-09-19 the mud can be Herschel–Bulkley (UI gap G-1, approved by the
user): `mud_consistency`, `mud_power_law_index`, `mud_yield_stress`. The default
is unchanged to the bit. `mud_viscosity` is now a property that **refuses** for
a non-Newtonian mud instead of returning κ̂ under a name that means Pa·s.

## Deviations from the mockup

Listed in full, with reasons, in `ui/INVENTORY.md` §E, and the mockup itself has
been revised to match what the solver produces. The short version: the snapshot
slider, the similarity collapse, η_N,
Δw_f, Pause, Snapshot and Compare/Overlay were removed because the solver has no
data behind them (the mud Herschel–Bulkley fields were on that list until
2026-09-19 and are now real — see above); the casing OD, the pumped-volume unit, the flow-rate control,
the validation counts and the blocked list were corrected; and the ZF23 "regime
map" was replaced by BF25 §3.3's Muskat criterion, which exists and is
computable from the setup inputs.

## Not built yet

Screens 2–4 (Case setup and launch, Live run monitor, Validation), in that
order. The live monitor is last on purpose: it is the only one that needs the
solver to emit progress, and what is worth emitting is clearer once the other
two exist.
