"""Screen 2 — Live run monitor.

Reads `status.json`, which the runner rewrites atomically from its existing
`on_step` hook. Nothing here talks to the solver, and nothing the solver does
depends on this page being open.

Two things the build spec asks for shape it.

**Invariants stay visible.** Progress is the least useful thing on the page. A
run that is 80% done with a volume error of 1e-3, a growing count of
under-converged Picard solves, or closure queries falling off the table is a
run to kill, and the only way to know is to be shown those numbers while it is
still running. They sit next to the progress bar, not behind a tab.

**A long run must not block the interface.** The child is detached and the page
is a reader, so refreshing, closing the browser or restarting this server
changes nothing about the run.

Honest gaps, stated on the page rather than filled: there are no intermediate
concentration fields — the solver does not emit them and this screen will not
invent one — and the phase before stepping (building the closure table, which
can be the longest part) reports a phase name and no percentage, because none
exists.
"""

from __future__ import annotations

import streamlit as st

from ui import launch
from ui.common import banner, page_chrome, panel, rows, unavailable

page_chrome("Run monitor")
st.title("Run monitor")

handles = launch.handles()
active = st.session_state.get("active_run")
if active and active["tag"] not in {h["tag"] for h in handles}:
    handles.insert(0, active)

if not handles:
    unavailable("No runs launched from this interface",
                "Start one on <strong>Case setup</strong>. A run started from a "
                "terminal is monitored the same way if it was given "
                "<code>--status-json PATH</code>; nothing else about it "
                "differs.")
    st.stop()

labels = {f'{h["tag"]}  ·  pid {h["pid"]}': h for h in handles}
chosen = st.sidebar.radio("Run", list(labels), index=0)
h = labels[chosen]

auto = st.sidebar.checkbox("Refresh every 3 s", value=True)
st.sidebar.markdown("---")

st.code(h["command"], language="bash")


@st.fragment(run_every=3 if auto else None)
def live(h=h):
    """Only this block reruns on the timer.

    The first version slept three seconds and called `st.rerun()`, which
    reruns the WHOLE page: every widget is dead for those three seconds and
    Streamlit shows the app as permanently busy.  A fragment reruns itself and
    leaves the rest of the page -- the run picker, the Stop button -- live.
    That is the build spec's "a long run must not block the interface", and a
    monitor that blocks is worse than no monitor.
    """
    _live_body(h)


def _live_body(h):
    status = launch.read_status(h["status"])
    running = launch.alive(h["pid"])

# --------------------------------------------------------------------------
# state
# --------------------------------------------------------------------------
    if status is None:
        if running:
            banner("warn", "<strong>Starting.</strong> The process is alive but has "
                           "not written a status document yet. Nothing is being "
                           "estimated for it.")
        else:
            banner("bad", "<strong>The process is gone and wrote no status.</strong> "
                          "It failed before it could report. The log below is all "
                          "there is.")
        st.text(launch.tail(h["log"], 40) or "(the log is empty)")
        return

    state = status.get("state", "unknown")
    if state == "preparing":
        banner("warn",
               f'<strong>Preparing — {status.get("phase", "?")}.</strong> '
               f'{status.get("note", "")}<br>No percentage is shown because none '
               f'exists: the closure table is built as a fixed set of gap-scale '
               f'solves before the first timestep, and it is cached, so a rebuild '
               f'happens only for a pair that has not been run before.')
    elif state == "running" and not running:
        banner("bad", "<strong>The process has died</strong> since its last status "
                      "write. The figures below are the last ones it recorded, not "
                      "the final state.")
    elif state == "failed":
        banner("bad", f'<strong>Failed.</strong> <code>{status.get("message", "")}'
                      f'</code>')
    elif state == "done":
        banner("good", f'<strong>Finished.</strong> Wrote '
                       f'<code>{status.get("message", "")}</code> — open it on '
                       f'<strong>Results</strong>.')

    # --------------------------------------------------------------------------
    # progress and invariants, side by side
    # --------------------------------------------------------------------------
    if state == "preparing":
        st.markdown('<span class="note">The panels below fill in once the first '
                    'timestep is taken. Until then they read &mdash;, not 0 and '
                    'not nan: the quantities do not exist yet.</span>',
                    unsafe_allow_html=True)
    if "fraction" in status:
        frac = float(status["fraction"])
        st.progress(min(1.0, max(0.0, frac)),
                    text=f'{frac:.1%} of the pumped volume  ·  step '
                         f'{status.get("step", 0):,}')

    def num(key, fmt="{:.4g}", scale=1.0):
        """A recorded value, or an em-dash.

        Not `nan`.  Before the first timestep these quantities do not EXIST -- the
        run is still building its closure table -- and printing `nan` says the
        opposite: that something was computed and came out undefined.  The
        distinction matters most in exactly this window, because a reader watching
        a slow start is deciding whether the run is broken.
        """
        v = status.get(key)
        if v is None:
            return "&mdash;"
        try:
            return fmt.format(float(v) * scale)
        except (TypeError, ValueError):
            return str(v)


    left, mid, right = st.columns(3, gap="medium")
    with left:
        eta_s = status.get("eta_wall_s")
        panel("Progress",
              rows((("t/Z", num("t_over_Z", "{:.4f}")),
                    ("step", num("step", "{:,.0f}")),
                    ("dt", num("dt")),
                    ("elapsed", num("elapsed_s", "{:.1f} min", 1 / 60)),
                    ("remaining",
                     f'{eta_s / 60:.0f} min' if eta_s else "&mdash;"))),
              "remaining is elapsed × (1 − fraction) / fraction — a linear "
              "extrapolation of the steps so far, not a model of the run")
    with mid:
        panel("Invariants",
              rows((("volume conservation", num("conservation_error", "{:.2e}")),
                    ("c&#772; range",
                     num("c_min", "{:.4f}") + " – " + num("c_max", "{:.4f}")),
                    ("outlet import",
                     num("outlet_import_volumes", "{:.5f} vol")),
                    ("Picard at cap", num("picard_unconverged", "{:,.0f}")),
                    ("worst Picard residual",
                     num("worst_picard_residual", "{:.2e}")))),
              "the worst value so far, not the current step — these are the "
              "numbers that decide whether the run is worth finishing")
    with right:
        panel("Regularisation and closures",
              rows((("static cells now", num("static_cells", "{:,.0f}")),
                    ("static cells worst", num("static_cells_max", "{:,.0f}")),
                    ("efficiency so far", num("efficiency", "{:.4f}")))),
              "NUM-13's mobility floor is load-bearing when static cells are "
              "non-zero, and PF04 §5 warns it flatters mud removal")
        rr = status.get("closure_table_range_report")
        if rr:
            banner("warn", f"<strong>Closure queries outside the table.</strong>"
                           f"<br><code>{rr}</code>")

    cmin, cmax = status.get("c_min"), status.get("c_max")
    if cmin is not None and (cmin < -1e-9 or cmax > 1 + 1e-9):
        banner("bad", f"<strong>The maximum principle is violated:</strong> "
                      f"c&#772; ∈ [{cmin:.6f}, {cmax:.6f}]. The LLF scheme is "
                      f"monotone under BCF25 (44), so this is a defect, not noise.")

    # --------------------------------------------------------------------------
    # what cannot be shown live
    # --------------------------------------------------------------------------
    st.markdown("### Field")
    unavailable("Live concentration field",
                "The runner's <code>on_step</code> receives the field, but writing "
                "it every couple of seconds would make the status document "
                "megabytes and the monitor a bottleneck on the run. Nothing is "
                "drawn in its place. The final field is on "
                "<strong>Results</strong>; the checkpoint is written every 180 s, "
                "so it is at most that far behind.")

    st.markdown("### Log")
    st.text(launch.tail(h["log"], 30) or "(the log is empty)")


live()

c1, c2 = st.columns([1, 5])
with c1:
    if launch.alive(h["pid"]) and st.button("Stop", key="monitor_stop"):
        launch.stop(h["pid"])
        st.warning("Sent SIGTERM. The last checkpoint is intact; the same "
                   "command resumes from it.")
with c2:
    st.markdown('<span class="note">There is no Pause. The solver has no '
                'paused state, and a button that claimed otherwise would '
                'misdescribe what it does. Stopping is not destructive: the '
                'run checkpoints every 180 s.</span>', unsafe_allow_html=True)
