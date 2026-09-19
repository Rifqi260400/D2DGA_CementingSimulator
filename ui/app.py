"""D2DGA solver UI — entry point.

    PYTHONPATH=. streamlit run ui/app.py

Four screens, in the build spec's order: set a case up and launch it, watch it,
read the result, and see what the model has and has not been shown to
reproduce.

The rules every screen keeps, and where each is enforced:

  * **One source for every parameter.** `ui/launch.py` introspects
    `scripts/kgep1_run.build_parser()`; nothing about a parameter is written
    down twice. Add an argument to the runner and it appears in the form.
  * **A UI run is a CLI run.** The launcher builds an argument vector and
    executes the runner as a subprocess. The command is shown before it runs
    and recorded with the result.
  * **The solver's numerics are never changed to suit a display.** The two
    things the UI needed that did not exist -- live progress and the running
    volume error -- were added as emissions (`--status-json`,
    `Simulation.conservation_error`), not as calculations in this layer.
  * **Post-processing is called, not reimplemented.** Every derived quantity
    comes from `d2dga.postprocess`, `d2dga.scaling` or `d2dga.muskat`.
  * **Degrade honestly.** A quantity a run did not record is shown as
    unavailable, with the reason. No plausible substitute, anywhere.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st                                 # noqa: E402

st.set_page_config(page_title="D2DGA — K-GEP-1", layout="wide",
                   initial_sidebar_state="expanded")

_HERE = os.path.dirname(os.path.abspath(__file__))
nav = st.navigation([
    st.Page(os.path.join(_HERE, "screens", "setup.py"),
            title="Case setup", icon=":material/tune:", default=True),
    st.Page(os.path.join(_HERE, "screens", "monitor.py"),
            title="Run monitor", icon=":material/monitoring:"),
    st.Page(os.path.join(_HERE, "screens", "results.py"),
            title="Results", icon=":material/insights:"),
    st.Page(os.path.join(_HERE, "screens", "validation.py"),
            title="Validation", icon=":material/verified:"),
])
nav.run()
