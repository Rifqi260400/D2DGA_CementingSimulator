# The D2DGA interface

```bash
pip install streamlit matplotlib
PYTHONPATH=. streamlit run ui/app.py
```

Four screens: **Case setup**, **Run monitor**, **Results**, **Validation**.

```
ui/
  app.py          entry point; st.navigation over the four screens
  common.py       presentation helpers (panels, banners, the B-1 disclaimer)
  theme.py        the visual system, defined once
  figures.py      matplotlib, arranging quantities that d2dga computed
  runs.py         finding and loading finished runs
  reports.py      parsing what the repo writes: gates.json, zf22_table3.md, BLOCKED.md
  launch.py       building a command, running it as a subprocess, reading its status
  screens/        setup.py  monitor.py  results.py  validation.py
```

## The five rules, and what enforces each

**One source for every parameter.** The setup form knows no defaults.
`launch.parameters()` introspects `scripts/kgep1_run.build_parser()` for names,
types, `choices` and help; `PARAM_META` and `FLUID_ARGS` add the units and
ranges argparse has no field for. Add an argument to the runner and it appears
in the form. *Gate: UI-1 asserts the form offers exactly the parser's
arguments, minus the three the launcher owns.*

**A UI run is a CLI run.** The launcher builds an argument vector and executes
the runner as an ordinary detached subprocess. The command is shown before it
runs and stored with the result. Every parameter is written out, including the
defaults, so the line keeps reproducing the run after the defaults move.
*Gate: UI-2 parses the launcher's argv back with `build_parser()` and requires
the same `Config`.*

**The numerics are never changed to suit a display.** Two things the interface
needed did not exist, and both were added as *emissions*, never as
calculations in this layer: `--status-json` (written from the existing
`on_step`) and `Simulation.conservation_error` (the running volume error, which
was a local variable). Recomputing that invariant here would have created a
second copy that could disagree with the report. *REM-12; gates UI-5, UI-9.*

**Post-processing is called, not reimplemented.** Every derived quantity comes
from `d2dga.postprocess`, `d2dga.scaling` or `d2dga.muskat`. Even
`breakthrough_at` is called on a reconstructed `RunResult` rather than copied.
The t_br error constants moved into `postprocess` because they were written
down twice.

**Degrade honestly.** A quantity a run did not record is shown as
*unavailable*, with the reason. No plausible substitute, anywhere. During the
phase before stepping the live panels read **—**, not `nan` and not `0`: those
numbers do not exist yet, and `nan` would say something was computed and came
out undefined.

## The disclaimer that travels with every number

Two statements, deliberately not merged.

**Provenance** — is this the fluid pair every published figure was computed
with? `FluidsConfig.departures_from_validated()` names each field that differs,
as `name: was -> is`, because "not the validated pair" alone does not tell a
reader whether to distrust the efficiency or the breakthrough time.

**Measurement** — how far the closure table's `θ = 0` slice is from the truth
for this pair (finding B-1). Only the table knows, and only when it is *built*
rather than loaded, so the screen reads the run's own `assumption_report()` out
of `metrics.json`. A run that recorded none gets a **warning, not a green
tick**, and the number is never estimated.

Only the second licenses a result, and the two can disagree: **the validated
pair is not a benign pair.** On the production velocity axis it measures
**9.8%**, which the solver's own grading calls SIGNIFICANT. The size is
governed by the viscosity ratio `m`, not by whether the mud is Newtonian — a
Newtonian mud at an ordinary 5 mPa·s already measures 11% — and the figure is
an **upper bound** on the angle error (θ = 0 against θ = 90°), not the error in
η_E. See `docs/assumptions.md` FLU-07.

## What each screen shows, and what it refuses to

| Screen | Shows | Refuses |
|---|---|---|
| **Case setup** | Every parameter; b, m, B, Fr*, Δρ̂, Z, Δξ from `Scaling`; the δ/π envelope; BF25 §3.3's Muskat verdict computed from this case's own closures; the `BLOCKED.md` items as blocked; the exact command | Re is greyed as "not used" rather than dropped — its absence is informative. The Muskat verdict is behind a control that states its cost, because it builds a real closure table; nothing is shown if it is not computed |
| **Run monitor** | Progress, and next to it the **invariants** — volume error, c̄ range, outlet import, Picard at cap, static cells — because those are what decide whether a run is worth finishing | No live concentration field: writing it every two seconds would make the monitor a bottleneck on the run. No Pause: the solver has no paused state |
| **Results** | η_E, η_N, the narrowest column, t_br at three thresholds with measured error bars, invariants, dimensionless groups, the field beside its δ/π envelope on a shared axis, ZF23 metrics from `metrics.json` | Snapshots, Ψ, and the measured Δw_f — the run keeps only its final field |
| **Validation** | Opens by saying this stage does **not** validate against CFD. Then **BCF25 §IV's conservation ledger** plotted as their Figs. 7/15 with their 10⁻¹⁵ line on the axes; `output/gates.json` with per-module counts; ZF22 Table 3 **9 / 10**, case 1 at −44.5 % linked to BLK-6; the open items; and what has no gate at all | It will not fall back to the prose in `docs/gate_status.md` when the gates have not been run. It will not show a run's older conservation number in the BCF25 panel: that one is normalised differently and is not comparable |

## Two defects that only the rendered page revealed

Both are worth recording because neither was visible in the source.

**Empty panels.** Streamlit closes every markdown block it is given, so an
opening `<div class="panel">` in its own call renders as a blank white box. A
panel is now always emitted in one call.

**9 / 9.** The Validation screen's first ZF22 parser filtered rows on
`startswith("| ")` and then sliced `[2:]`. The table file writes its separator
as `|---|---|`, with no space — so the filter removed the separator and the
slice removed **case 1**, the one case that does not reproduce. The screen
displayed a clean 9 / 9, exactly the mockup's error, arrived at independently
by discarding its own counter-evidence. Parsing moved into `ui/reports.py`
where it can be tested without a browser, the separator is found by content,
and UI-12/13/14 lock it.

## The conservation panel, and what it does not prove

BCF25 §IV computes each fluid's volume two ways — integrated over the
interior, and accumulated from the boundary fluxes — normalises both by the
total annulus volume, and differences them. Their Figs. 7 and 15 plot that
against time and report **~10⁻¹⁵**. The Validation screen runs the same
computation on this solver and draws the same figure, with their level marked.

Two things the panel says out loud rather than letting the picture imply:

* **The two fluid curves being mirror images is not corroboration.** BCF25
  evolve K = 3 concentrations independently. Here K = 2 and only `c̄₂` is
  evolved, so `c̄₁ ≡ 1 − c̄₂` and fluid 1's flux is identically the total minus
  fluid 2's — the LLF dissipation changes sign with `c` and cancels. The
  independent content is the third curve, the **total-flux imbalance**, which
  asks whether the elliptic solve delivered the same Q at both ends. Measured:
  exactly 0.0.
* **A longer run sits higher, and that is roundoff.** √N·ε is 6.5 × 10⁻¹⁴ at
  85 000 steps, so the panel reads the measured value against that rather than
  against 10⁻¹⁵ flat, and says which it is doing.

Runs recorded before 2026-09-19 have no ledger. They are shown as having none
— their older conservation number divides by the volume present rather than
the annulus, and presenting it in this panel would be comparing two different
quantities against BCF25's one.

## Honest limits

* No intermediate snapshots and no Ψ are persisted, so there is no time
  slider, no similarity collapse and no streamline plot.
* ZF23 dispersion metrics come only from `metrics.json`: they are computed
  during a run from a pre-breakthrough state, and by the final state the front
  has left the domain, so recomputing them from what is stored would produce a
  meaningless number. Runs predating `runio` have none.
* The colour scale defaults to the absolute 0–1 ramp. For a nearly complete
  displacement that is a flat block — true, and useless — so a checkbox
  stretches it to the run's range. The colour bar states which is in use.
* Nothing here validates the model against a K-GEP-1 *measurement*, because
  there is none. The Validation screen says so in those words.
