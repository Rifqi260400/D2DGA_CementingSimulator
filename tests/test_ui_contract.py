"""The contract between the runner and the interface.

The UI makes three claims that are easy to assert and easy to quietly break:

  1. **It restates no parameter.** Everything it offers comes from
     `build_parser()`; if the two ever disagree, the form is editing something
     the runner ignores.
  2. **A run launched from the interface is a run launched from a terminal.**
     The argument vector the launcher builds must parse back into the same
     `Config`.  This is the whole basis of "reproducible from the CLI", and
     without a test it is a sentence in a README.
  3. **It never fills a gap with a plausible value.** The mechanisms that
     make that possible -- the emitted invariant, the status document, the
     measured t_br error -- are tested here rather than trusted.

Nothing in this file touches the physics; it tests the seams.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, "scripts")

import kgep1_run as runner                                     # noqa: E402
from d2dga.config import Config, FluidsConfig, GridConfig      # noqa: E402
from d2dga.elliptic import NewtonianClosures                   # noqa: E402
from d2dga.geometry import build_geometry                      # noqa: E402
from d2dga.postprocess import (TBR_RELATIVE_ERROR_AT_80,       # noqa: E402
                               tbr_relative_error)
from d2dga.simulation import Simulation                        # noqa: E402
from ui import launch                                          # noqa: E402
from ui.reports import (parse_markdown_table, read_blocked,    # noqa: E402
                        zf22_comparison)


# ------------------------------------------------------------------ UI-1 ---
def test_ui1_the_form_offers_exactly_what_the_runner_accepts():
    """No parameter is written down twice, and none is missing.

    `parameters()` is introspection, so the only way it can drift is if the
    launcher starts filtering something it should not.  The three it owns --
    help, the status path, the resume switch -- are named explicitly here, so
    adding a fourth exclusion is a deliberate act that fails this test first.
    """
    offered = {p["dest"] for p in launch.parameters()}
    parser = {a.dest for a in runner.build_parser()._actions}
    assert offered == parser - {"help", "status_json", "no_resume"}
    # and every fluid field of the dataclass is reachable from the form,
    # otherwise the UI can express a case the config cannot
    fluid_fields = {d for _, d, _, _ in runner.FLUID_ARGS}
    assert fluid_fields <= offered
    assert fluid_fields == {f for f in vars(FluidsConfig()).keys()
                            if f.startswith(("mud_", "cement_"))}


# ------------------------------------------------------------------ UI-2 ---
@pytest.mark.parametrize("edits", [
    {},
    {"w0": 0.5, "n_phi": 12, "inflow": "uniform"},
    {"mud_density": 1400.0, "mud_consistency": 0.02,
     "mud_power_law_index": 0.7, "mud_yield_stress": 4.79},
    {"eccentricity": 0.55, "cement_density": 1900.0, "wall": "caliper"},
])
def test_ui2_a_launched_command_parses_back_into_the_same_case(edits):
    """The claim "a UI run is a CLI run", made testable.

    The launcher builds an argv; the runner's own parser reads it; the `Config`
    that comes out must equal the one the form describes.  Round-tripping
    through `build_parser()` rather than comparing strings is the point: it is
    the same parsing the real process does.
    """
    values = {p["dest"]: p["default"] for p in launch.parameters()}
    values.update(edits)
    argv = launch.argv_from(values)
    # drop the interpreter and the script path -- the parser sees only flags
    parsed = runner.build_parser().parse_args(argv[2:])
    assert runner.config_from_args(parsed) == runner.config_from_args(
        runner.build_parser().parse_args(
            [x for d, v in values.items()
             for x in (next(p["flag"] for p in launch.parameters()
                            if p["dest"] == d), str(v))]))
    for dest, want in values.items():
        got = getattr(parsed, dest)
        assert got == type(got)(want), dest


def test_ui2b_the_command_names_every_parameter_not_just_the_changed_ones():
    """A command that relies on defaults stops reproducing when they move."""
    line = launch.shell_command(launch.argv_from({}))
    for p in launch.parameters():
        assert p["flag"] in line, p["flag"]
    assert line.startswith("PYTHONPATH=. python scripts/kgep1_run.py")


# ------------------------------------------------------------------ UI-3 ---
def test_ui3_fluid_arguments_build_the_dataclass_and_inherit_its_validation():
    ap = runner.build_parser()
    assert runner.fluids_from_args(ap.parse_args([])) == FluidsConfig()
    got = runner.fluids_from_args(ap.parse_args(
        ["--mud-yield-stress", "4.79", "--mud-consistency", "0.02",
         "--mud-power-law-index", "0.7"]))
    assert got.mud_yield_stress == 4.79 and got.mud_name.endswith("(Herschel-Bulkley)")
    # the dataclass's own refusal, not a second copy in the script
    with pytest.raises(ValueError, match=r"must be in \(0, 1\]"):
        runner.fluids_from_args(ap.parse_args(["--mud-power-law-index", "1.4"]))


# ------------------------------------------------------------------ UI-4 ---
def test_ui4_checkpoints_of_different_physics_do_not_collide():
    """Two cases differing only in a fluid wrote to the SAME file.

    Since R-2 that is refused rather than silently resumed -- correct, but it
    also made the second case impossible to run.  The suffix separates them.
    The DEFAULT must keep its historical name, or every checkpoint recorded
    before today stops resolving.
    """
    base = Config(grid=GridConfig(16, 80))
    assert runner._physics_suffix(base) == ""
    a = Config(grid=GridConfig(16, 80),
               fluids=FluidsConfig(mud_density=1400.0))
    b = Config(grid=GridConfig(16, 80),
               fluids=FluidsConfig(mud_density=1401.0))
    assert runner._physics_suffix(a).startswith("_f")
    assert runner._physics_suffix(a) != runner._physics_suffix(b)
    # and it depends only on the physics, not on the mesh
    assert runner._physics_suffix(a) == runner._physics_suffix(
        Config(grid=GridConfig(4, 4), fluids=FluidsConfig(mud_density=1400.0)))


# ------------------------------------------------------------------ UI-5 ---
def _tiny_sim():
    geo = build_geometry(Config(grid=GridConfig(8, 20)))
    return geo, Simulation(geo, NewtonianClosures(m=0.5), froude=1.0,
                           delta_rho=0.0)


def test_ui5_the_running_volume_error_is_emitted_not_recomputed():
    """`Simulation.conservation_error` must equal the returned one at every
    step, or a monitor and a report would disagree about the same run."""
    geo, sim = _tiny_sim()
    seen = []
    res = sim.run(t_end=0.05 * geo.grid.Z,
                  on_step=lambda r, c, p: seen.append(sim.conservation_error))
    assert seen, "on_step never fired"
    assert seen == sorted(seen), "the worst-so-far must be monotone"
    assert seen[-1] == res.conservation_error == sim.conservation_error


def test_ui6_run_result_reports_its_step_count_even_when_it_took_no_steps(tmp_path):
    """Regression: relaunching a FINISHED run raised IndexError.

    Callers read `reports[-1].n`, which is empty when a resume finds the
    checkpoint already at t_end -- and a launcher does that routinely, because
    pressing Run twice is free.  It was also silently short by up to
    `record_every - 1` on every other run.
    """
    geo, sim = _tiny_sim()
    ckpt = str(tmp_path / "c.npz")
    t_end = 0.05 * geo.grid.Z
    first = sim.run(t_end=t_end, checkpoint_path=ckpt)
    assert first.steps == first.reports[-1].n > 0

    _, sim2 = _tiny_sim()
    again = sim2.run(t_end=t_end, checkpoint_path=ckpt)
    assert again.reports == []                 # nothing left to do
    assert again.steps == first.steps          # and it still knows how far it got
    assert np.isfinite(again.conservation_error)


def test_ui7_step_count_is_not_the_last_recorded_report(tmp_path):
    """With `record_every > 1` the last REPORT is not the last STEP."""
    geo, sim = _tiny_sim()
    res = sim.run(t_end=0.05 * geo.grid.Z, record_every=7)
    assert res.steps >= res.reports[-1].n
    assert res.steps == len(res.times) - 1


# ------------------------------------------------------------------ UI-8 ---
def test_ui8_tbr_error_is_half_order_and_refuses_an_unmeasured_threshold():
    """The constants are a MEASUREMENT (A-2), so they live in one place and
    an unmeasured threshold is refused rather than interpolated."""
    for th, ref in TBR_RELATIVE_ERROR_AT_80.items():
        assert tbr_relative_error(th, 80) == ref
        # half order: four times the cells halves the error
        assert tbr_relative_error(th, 320) == pytest.approx(ref / 2.0)
    with pytest.raises(KeyError, match="no measured t_br error"):
        tbr_relative_error(0.25, 80)


# ------------------------------------------------------------------ UI-9 ---
def test_ui9_the_status_document_is_atomic_and_carries_the_invariants(tmp_path):
    """A monitor must never read half a document, and must be shown the
    invariants rather than only the progress."""
    geo, sim = _tiny_sim()
    path = str(tmp_path / "status.json")
    runner.write_status_stub(path, "closure table", "note")
    doc = json.load(open(path))
    assert doc["state"] == "preparing" and doc["phase"] == "closure table"
    assert not os.path.exists(path + ".tmp")

    writer = runner.StatusWriter(path, sim, geo, capacity=1.0,
                                 t_end=0.05 * geo.grid.Z, command="x", every=0.0)
    sim.run(t_end=0.05 * geo.grid.Z,
            on_step=lambda r, c, p: writer.update(r, force=True))
    doc = json.load(open(path))
    assert doc["state"] == "running"
    for key in ("conservation_error", "outlet_import_volumes",
                "picard_unconverged", "worst_picard_residual", "c_min", "c_max",
                "static_cells_max", "fraction", "step"):
        assert key in doc, key
    assert doc["conservation_error"] == sim.conservation_error
    assert 0.0 <= doc["fraction"] <= 1.0

    writer.finish("done", "ckpt.npz")
    doc = json.load(open(path))
    assert doc["state"] == "done" and doc["message"] == "ckpt.npz"
    # finishing must not discard what the run recorded
    assert "conservation_error" in doc


# ----------------------------------------------------------------- UI-10 ---
def test_ui10_the_blocked_list_is_parsed_from_the_document():
    """The screens must not carry their own copy of what is unresolved."""
    items = read_blocked()
    ids = [b["id"] for b in items]
    assert ids == sorted(ids, key=lambda s: int(s.split("-")[1]))
    assert "BLK-6" in ids            # ZF22 case 1, the headline open item
    kinds = {b["kind"] for b in items}
    assert kinds == {"missing physical data", "irreducible ambiguity"}
    assert read_blocked("/nonexistent/BLOCKED.md") == []


def test_ui11_the_muskat_predictor_is_importable_without_the_test_tree():
    """G-5: a screen must not import from `tests/`."""
    import d2dga.muskat as m
    assert hasattr(m, "dw_at_leading_edge") and hasattr(m, "classify")
    assert "tests" not in m.__file__.split(os.sep)


# ----------------------------------------------------------------- UI-12 ---
def test_ui12_the_zf22_table_keeps_every_row_including_the_failing_one():
    """The Validation screen's counter-evidence must survive parsing.

    The first version filtered rows with `startswith("| ")` and then dropped
    two more by position.  `output/zf22_table3.md` writes its separator as
    `|---|---|`, with no space, so the filter removed the separator and the
    slice removed CASE 1 -- the one case that does not reproduce.  The screen
    then displayed **9 / 9**: a clean bill of health manufactured by discarding
    the evidence against it.  Nothing about that was visible in the code.
    """
    z = zf22_comparison()
    assert z["n_total"] == 10, "all ten published cases must be present"
    assert [r["case"] for r in z["rows"]] == [str(i) for i in range(1, 11)]
    assert z["n_ok"] == 9 and z["failing"] == ["1"], (z["n_ok"], z["failing"])
    case1 = z["rows"][0]
    assert not case1["ok"] and case1["delta"] < -0.4, case1
    # and no band that could be called defensible rescues it
    assert abs(case1["eta"] - case1["zf22_eta"]) > 0.25 * abs(case1["zf22_eta"])


@pytest.mark.parametrize("sep", ["|---|---|", "| --- | --- |", "|:--|--:|"])
def test_ui13_the_separator_is_found_by_content_not_position(sep):
    """Every markdown separator spelling, since the bug was one spelling."""
    header, rows = parse_markdown_table(
        ["# title", "", "| a | b |", sep, "| 1 | 2 |", "| 3 | 4 |", "text"])
    assert header == ["a", "b"]
    assert rows == [["1", "2"], ["3", "4"]]


def test_ui14_an_unreadable_row_counts_as_not_reproduced(tmp_path):
    """A row whose numbers cannot be parsed must not be silently a pass."""
    f = tmp_path / "t.md"
    f.write_text("| case | eta_E | ZF22 eta_E |\n|---|---|---|\n"
                 "| 1 | 0.90 | 0.91 |\n| 2 | n/a | 0.80 |\n")
    z = zf22_comparison(str(f))
    assert z["n_total"] == 2 and z["n_ok"] == 1 and z["failing"] == ["2"]


# ----------------------------------------------------------------- UI-15 ---
def test_ui15_a_run_with_non_default_physics_is_discoverable(tmp_path):
    """The suffix that stops checkpoints colliding must not hide them.

    `kgep1_run` appends `_f<hash>` when the fluids or standoff are not the
    shipped defaults (UI-4).  `runs.discover` matched only the old name, so a
    run on a real mud completed, wrote its config and metrics -- and simply did
    not appear in the interface.  Nothing failed; the row was absent.  Caught
    by trying to open a run the interface had just produced.
    """
    from ui import runs as uiruns

    for name in ("kgep1_ckpt_synthetic_no_axial_gradient_8x24_w0.2_v0.35.npz",
                 "kgep1_ckpt_synthetic_no_axial_gradient_8x24_w0.2_v0.35_fedbaa3.npz",
                 "kgep1_ckpt_caliper_uniform_16x80_w0.5_v1.2_f0a1b2c.npz"):
        (tmp_path / name).write_bytes(b"")
    # things that must NOT be picked up: a cached closure table, an unrelated
    # file, and a suffix the runner cannot produce (it writes exactly 6 hex)
    (tmp_path / "not_a_run.npz").write_bytes(b"")
    (tmp_path / "closure_table_abc123.npz").write_bytes(b"")
    (tmp_path / "kgep1_ckpt_synthetic_uniform_8x24_w0.2_v0.35_f0a1b2c3.npz"
     ).write_bytes(b"")

    found = {r.name for r in uiruns.discover(str(tmp_path))}
    assert "synthetic_no_axial_gradient_8x24_w0.2_v0.35" in found
    assert "synthetic_no_axial_gradient_8x24_w0.2_v0.35_fedbaa3" in found
    assert "caliper_uniform_16x80_w0.5_v1.2_f0a1b2c" in found
    assert len(found) == 3, found
    # the parsed arguments must be unaffected by the suffix
    by_name = {r.name: r for r in uiruns.discover(str(tmp_path))}
    a = by_name["synthetic_no_axial_gradient_8x24_w0.2_v0.35"].args
    b = by_name["synthetic_no_axial_gradient_8x24_w0.2_v0.35_fedbaa3"].args
    for k in ("wall", "inflow", "n_phi", "n_xi", "w0", "volumes"):
        assert a[k] == b[k], k


def test_ui16_every_checkpoint_the_runner_can_write_is_matched():
    """Generated from the runner's own naming, not from examples.

    The two must not be able to drift: the interface's pattern is the only
    thing standing between a finished run and a run nobody can open.
    """
    from ui.runs import _NAME

    for fluids, standoff in ((FluidsConfig(), None),
                             (FluidsConfig(mud_density=1400.0), None)):
        cfg = Config(grid=GridConfig(16, 80), fluids=fluids)
        suffix = runner._physics_suffix(cfg)
        for wall in ("synthetic", "caliper"):
            for inflow in ("no_axial_gradient", "uniform"):
                name = (f"kgep1_ckpt_{wall}_{inflow}_16x80_w0.2_v1.2"
                        f"{suffix}.npz")
                assert _NAME.search(name), name


# ----------------------------------------------------------------- UI-17 ---
def test_ui17_the_reproduce_line_actually_reproduces_the_run(tmp_path):
    """"Reproduce from the CLI" has to be true, not decorative.

    The first version listed six flags by hand.  A run on a non-default mud
    was therefore shown a command that reproduces a DIFFERENT run, under that
    exact heading -- worse than showing nothing, because a reader would copy
    it.  The command is now built by the launcher from the recorded payload,
    so this test can close the loop: parse the displayed line back and require
    the same `Config`.
    """
    import shlex

    from d2dga import runio
    from ui import runs as uiruns

    fluids = FluidsConfig(mud_density=1400.0, mud_consistency=0.02,
                          mud_power_law_index=0.7, mud_yield_stress=4.79)
    cfg = Config(grid=GridConfig(8, 24), fluids=fluids)
    args = {"wall": "synthetic", "inflow": "uniform", "n_phi": 8, "n_xi": 24,
            "w0": 0.5, "volumes": 0.35, "cfl": 0.25, "n_c": 9, "n_h": 3}
    ckpt = str(tmp_path / "kgep1_ckpt_synthetic_uniform_8x24_w0.5_v0.35_fabcdef.npz")
    np.savez(ckpt, c=np.zeros((8, 24)))
    runio.write_config(ckpt, runio.config_payload(cfg, args))

    run = uiruns.discover(str(tmp_path))[0]
    assert run.config_source == uiruns.RECORDED

    line = run.command()
    tokens = shlex.split(line)
    assert tokens[0] == "PYTHONPATH=."
    parsed = runner.build_parser().parse_args(tokens[3:])
    assert runner.config_from_args(parsed) == cfg
    for key, want in args.items():
        assert getattr(parsed, key) == want, key
