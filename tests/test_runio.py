"""
Run provenance and the resume guard, `d2dga/runio.py`.

The hazard these tests pin: `Simulation.run` resumed whenever the checkpoint
path existed and the grid SHAPE matched.  Fluid properties are in neither the
path nor the shape, so editing `mud_density` and relaunching the same command
silently continued a concentration field computed under the old physics.
"""

import json

import numpy as np
import pytest

from d2dga import runio
from d2dga.config import Config, FluidsConfig, GridConfig
from tests.test_m6_simulation import build


# ------------------------------------------------------- the payload ------
def test_config_payload_round_trips_through_json(tmp_path):
    cfg = Config(grid=GridConfig(8, 40))
    payload = runio.config_payload(cfg, {"wall": "synthetic", "w0": 0.2})
    ckpt = str(tmp_path / "run.npz")
    assert runio.read_config(ckpt) is None          # nothing written yet
    runio.write_config(ckpt, payload)
    back = runio.read_config(ckpt)
    assert back["config"]["grid"]["n_xi"] == 40
    assert back["config"]["fluids"]["mud_density"] == 998.0
    assert back["args"]["w0"] == 0.2
    assert back["environment"]["numpy"] == np.__version__
    assert runio.fingerprint(back) == runio.fingerprint(payload)


def test_fingerprint_is_stable_and_sensitive_to_physics_only():
    """The tag must ignore where a run writes and which machine ran it, and
    must not ignore anything that changes the answer."""
    cfg = Config(grid=GridConfig(8, 40))
    args = {"wall": "synthetic", "w0": 0.2, "n_phi": 8, "n_xi": 40,
            "cfl": 0.5, "volumes": 1.2, "inflow": "no_axial_gradient"}
    base = runio.fingerprint(runio.config_payload(cfg, args))

    # same physics, different invocation -> same tag
    same = runio.config_payload(cfg, dict(args, no_resume=True),
                                command=["somewhere/else.py", "--quiet"])
    same["environment"]["git_sha"] = "0" * 40
    assert runio.fingerprint(same) == base

    # each of these changes the answer -> must change the tag
    assert runio.fingerprint(runio.config_payload(
        Config(grid=GridConfig(8, 40), fluids=FluidsConfig(mud_density=1400.0)),
        args)) != base
    assert runio.fingerprint(runio.config_payload(
        Config(grid=GridConfig(8, 80)), args)) != base
    assert runio.fingerprint(runio.config_payload(
        cfg, dict(args, cfl=0.25))) != base
    assert runio.fingerprint(runio.config_payload(
        cfg, dict(args, inflow="uniform"))) != base


def test_describe_mismatch_names_the_setting_that_moved():
    """A refusal saying 'hashes differ' is not actionable."""
    a = runio.config_payload(Config(), {"w0": 0.2})["physical"]
    b = runio.config_payload(
        Config(fluids=FluidsConfig(mud_density=1400.0)), {"w0": 0.5})["physical"]
    diffs = runio.describe_mismatch(a, b)
    joined = " | ".join(diffs)
    assert "config.fluids.mud_density: 998.0 -> 1400.0" in joined
    assert "args.w0: 0.2 -> 0.5" in joined


# ------------------------------------------------------- the guard --------
def test_check_resume_levels(tmp_path):
    cfg = Config(grid=GridConfig(8, 40))
    args = {"wall": "synthetic", "w0": 0.2, "n_phi": 8, "n_xi": 40}
    payload = runio.config_payload(cfg, args)
    ckpt = str(tmp_path / "run.npz")

    # no checkpoint at all
    assert runio.check_resume(ckpt, payload)[0] is None

    # a legacy checkpoint: no recorded config -> warn, do not block
    np.savez_compressed(ckpt, c=np.zeros((8, 40)))
    level, msg = runio.check_resume(ckpt, payload)
    assert level == "warn" and "LEGACY" in msg

    # matching config -> silent
    runio.write_config(ckpt, payload)
    assert runio.check_resume(ckpt, payload)[0] is None

    # changed fluids -> block, and say which setting
    other = runio.config_payload(
        Config(grid=GridConfig(8, 40),
               fluids=FluidsConfig(mud_density=1400.0)), args)
    level, msg = runio.check_resume(ckpt, other)
    assert level == "block"
    assert "mud_density" in msg and "998.0 -> 1400.0" in msg


def test_simulation_refuses_to_resume_across_a_settings_change(tmp_path):
    """
    The scenario, end to end, through `Simulation.run` itself -- the backstop
    for when a caller forgets to pre-check.

    Against the old code this test FAILS: `run` took no `checkpoint_tag`, so the
    second call resumed the first run's field with different settings and
    returned quietly.
    """
    ckpt = str(tmp_path / "ckpt.npz")
    geo, sim = build(n_phi=8, n_xi=40, e=0.3, b=5.0, m=0.5)
    sim.run(t_end=0.10 * geo.grid.Z, record_every=0,
            checkpoint_path=ckpt, checkpoint_every=0.0, checkpoint_tag="tag-A")
    assert "tag" in np.load(ckpt).files

    # same grid shape, so the OLD shape check passes; different settings tag
    geo2, sim2 = build(n_phi=8, n_xi=40, e=0.3, b=5.0, m=0.5)
    with pytest.raises(ValueError, match="written with settings tag"):
        sim2.run(t_end=0.20 * geo2.grid.Z, record_every=0,
                 checkpoint_path=ckpt, checkpoint_tag="tag-B")

    # the matching tag still resumes
    geo3, sim3 = build(n_phi=8, n_xi=40, e=0.3, b=5.0, m=0.5)
    res = sim3.run(t_end=0.20 * geo3.grid.Z, record_every=0,
                   checkpoint_path=ckpt, checkpoint_tag="tag-A")
    assert res.reports[-1].t >= 0.20 * geo3.grid.Z


def test_a_legacy_checkpoint_still_resumes_with_a_warning(tmp_path, capsys):
    """Refusing a tagless checkpoint would break relaunching a run finished
    before this module existed, so it proceeds noisily instead."""
    ckpt = str(tmp_path / "legacy.npz")
    geo, sim = build(n_phi=8, n_xi=40, e=0.3, b=5.0, m=0.5)
    sim.run(t_end=0.10 * geo.grid.Z, record_every=0,
            checkpoint_path=ckpt, checkpoint_every=0.0)      # no tag: legacy
    assert str(np.load(ckpt)["tag"]) == ""

    geo2, sim2 = build(n_phi=8, n_xi=40, e=0.3, b=5.0, m=0.5)
    sim2.run(t_end=0.20 * geo2.grid.Z, record_every=0,
             checkpoint_path=ckpt, checkpoint_tag="tag-A")
    assert "carries no settings tag" in capsys.readouterr().out


def test_metrics_round_trip(tmp_path):
    ckpt = str(tmp_path / "run.npz")
    assert runio.read_metrics(ckpt) is None
    runio.write_metrics(ckpt, {"eta_E": 0.9946, "zf23": {"is_dispersive": False}})
    got = runio.read_metrics(ckpt)
    assert got["eta_E"] == 0.9946
    assert got["zf23"]["is_dispersive"] is False
    assert json.loads(open(runio.metrics_path_for(ckpt)).read())["eta_E"] == 0.9946


def test_describe_mismatch_ignores_tuple_versus_list_from_json_round_trip():
    """A tuple written to JSON comes back a list.  The fingerprint never cared
    -- both serialise identically -- but the diff reported
    `mean_velocities_m_s: [0.05, 0.2, 0.5] -> (0.05, 0.2, 0.5)` as a change,
    which is noise in exactly the message that has to be trusted."""
    live = runio.config_payload(Config(), {"cfl": 0.5})["physical"]
    stored = json.loads(json.dumps(live))           # what read_config returns
    assert runio.describe_mismatch(stored, live) == []
    # and a real change still shows through the same path
    changed = runio.config_payload(Config(), {"cfl": 0.25})["physical"]
    assert runio.describe_mismatch(stored, changed) == ["args.cfl: 0.5 -> 0.25"]
