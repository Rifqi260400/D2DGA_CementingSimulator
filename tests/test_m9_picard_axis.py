"""M9 gates -- A-3b: the closure table's velocity axis and the Picard iteration.

What happened
-------------
A-3 anchored the table's velocity axis at zero, so that a query below the axis
became impossible.  It was anchored with a SINGLE node -- `[0, 0.02]` -- and the
note written at the time said it "changes nothing".  That was wrong.  The gap
solve at `umag = 0` and `0.02` agrees to six significant figures, so that
segment is nearly flat while the next one, `[0.02, 0.1]`, is not.  The kink
between them breaks the elliptic Picard iteration, which updates `|u_bar|` from
the previous `Psi` and re-evaluates the mobility: it oscillates across the kink
instead of contracting.

Measured on the K-GEP-1 pair at 16 x 80, same field, only the axis changed:

    [0.02] + geom               16 nodes   picard   6   residual 1.54e-09
    [0, 0.02] + geom            17 nodes   picard 100   residual 3.37e-05
    [0, .005, .01, .02] + geom  19 nodes   picard   6   residual 1.54e-09

and along a fresh 2982-step trajectory, steps hitting the iteration cap:
37.5% on the broken axis, 1.5% on the refined one, against 0% reported by the
pre-A-3 production run.

Why there is no direct gate here, and what these gates do instead
-----------------------------------------------------------------
Two attempts at a cheap unit test that reproduces the failure both FAILED to
discriminate, and they are recorded here rather than quietly dropped:

  * a reduced table (11 c-nodes, 3 H-nodes) gives 52 iterations on BOTH axes --
    too coarse everywhere for this particular kink to dominate;
  * a synthetic flat front at 16 x 80 gives 100 on BOTH axes -- a harsher test
    than the real flow, so it cannot tell them apart.

Reproducing it needs a production-resolution table (~10 min to build) and a
field from a real trajectory.  That is not a unit test.  So what is gated here
is (T1) the axis a run will actually use, as a regression lock on the measured
boundary, and (T2) that the axis is inside the run fingerprint, which is a
genuine functional property and was genuinely absent.  The failure itself is
caught at RUN level: `kgep1_run` reports the cap count prominently and the
monitor shows it live.  Saying this plainly is the point -- a gate that cannot
fail would be worse than admitting the gap.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

sys.path.insert(0, "scripts")

import kgep1_run as runner                                    # noqa: E402
from d2dga import runio                                       # noqa: E402
from d2dga.config import Config, GridConfig                   # noqa: E402


# ---------------------------------------------------------------- M9-T1 ---
def test_m9_t1_the_velocity_axis_resolves_the_region_below_the_first_decade():
    """Regression lock on a MEASURED boundary, not a derived criterion.

    Two nodes below 0.1 (`[0, 0.02]`) measurably breaks the Picard iteration on
    the production case; four (`[0, .005, .01, .02]`) measurably does not.  The
    bound is four because that is what was measured, and this docstring says so
    rather than dressing it up as theory.

    It also re-asserts A-3 itself: the axis starts at exactly zero, so an
    under-range query is impossible -- `|u_bar| >= 0` always.
    """
    axis = np.asarray(runner.closure_umag_axis(), dtype=float)

    assert axis[0] == 0.0, "A-3: anchored at zero, or under-range queries return"
    assert np.all(np.diff(axis) > 0.0), "the axis must be strictly increasing"

    below = axis[axis < 0.1]
    assert len(below) >= 4, (
        f"only {len(below)} nodes below 0.1: {below}.  Two of them "
        f"([0, 0.02]) put the elliptic Picard iteration at its cap on 37.5% of "
        f"steps; four put it at 1.5%.  See A-3b.")

    # and no segment below the first decade may be a large multiple of the one
    # before it -- that ratio is what makes a kink out of an added node
    seg = np.diff(axis[axis <= 0.1])
    assert np.max(seg[1:] / seg[:-1]) <= 10.0, (
        f"segment length ratios below 0.1: {seg[1:] / seg[:-1]}")


def test_m9_t1b_the_axis_reaches_far_enough_for_this_fluid_pair():
    """FLU-06: b*I1 ~ 1500 drives |u_bar| past 370, so an axis that stops at
    O(1) -- ample for the published cases -- is off the end within the first
    timesteps."""
    axis = np.asarray(runner.closure_umag_axis(), dtype=float)
    assert axis[-1] >= 3000.0, axis[-1]


# ---------------------------------------------------------------- M9-T2 ---
def test_m9_t2_the_axis_is_inside_the_run_fingerprint():
    """The gap R-2 could not see.

    R-2 refuses to resume a checkpoint whose settings changed, by comparing a
    fingerprint of the physics.  The closure table's velocity axis IS physics --
    changing it changes the answer, as A-3b measured -- but it lives in neither
    `Config` nor the CLI arguments, so the fingerprint did not cover it and a
    checkpoint computed with the broken axis would have resumed silently under
    the fixed one.  The guard only ever sees what it is handed; the runner now
    hands it the axis.
    """
    assert "closure_umag_axis" in runio._PHYSICAL_ARG_KEYS

    args = dict(vars(runner.build_parser().parse_args(["--n-xi", "80"])))
    args["closure_umag_axis"] = [float(x) for x in runner.closure_umag_axis()]
    cfg = Config(grid=GridConfig(16, 80))
    good = runio.fingerprint(runio.config_payload(cfg, args))

    broken = dict(args)
    broken["closure_umag_axis"] = [0.0, 0.02] + args["closure_umag_axis"][4:]
    assert runio.fingerprint(runio.config_payload(cfg, broken)) != good

    # and a mismatch must be REPORTED by name, not just refused
    msg = "; ".join(runio.describe_mismatch(
        runio.config_payload(cfg, args)["physical"],
        runio.config_payload(cfg, broken)["physical"], ""))
    assert "closure_umag_axis" in msg, msg


# ---------------------------------------------------------------- M9-T3 ---
def test_m9_t3_the_runner_records_the_axis_it_used():
    """A recorded run must say which axis produced it, or the fingerprint is
    the only evidence and a reader cannot see what moved."""
    args = dict(vars(runner.build_parser().parse_args([])))
    args["closure_umag_axis"] = [float(x) for x in runner.closure_umag_axis()]
    payload = runio.config_payload(Config(), args)
    recorded = payload["physical"]["args"]["closure_umag_axis"]
    assert len(recorded) == len(runner.closure_umag_axis())
    assert recorded[0] == 0.0 and recorded[-1] >= 3000.0


# ---------------------------------------------------------------- M9-T4 ---
@pytest.mark.parametrize("axis,ok", [
    ([0.0, 0.005, 0.01, 0.02, 0.1, 3000.0], True),
    ([0.0, 0.02, 0.1, 3000.0], False),          # the broken A-3 axis
    ([0.02, 0.1, 3000.0], False),               # the pre-A-3 axis: no zero node
    ([0.0, 0.001, 0.05, 0.1, 3000.0], False),   # a 50x jump below 0.1
])
def test_m9_t4_the_axis_check_would_reject_each_axis_this_project_has_used(axis, ok):
    """The lock has to be able to fail, or it is decoration.  Each historical
    axis is put through it, including the two that were actually shipped."""
    a = np.asarray(axis, dtype=float)
    below = a[a < 0.1]
    seg = np.diff(a[a <= 0.1])
    passes = bool(a[0] == 0.0 and len(below) >= 4
                  and (len(seg) < 2 or np.max(seg[1:] / seg[:-1]) <= 10.0))
    assert passes == ok, (axis, passes)
