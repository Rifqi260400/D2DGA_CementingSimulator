"""Starting a run from the interface, and watching it.

The one rule that shapes this module: **a run launched here and a run launched
from a terminal must be the same invocation.** So the UI does not import the
solver and call it in-process; it builds an argument vector for
`scripts/kgep1_run.py` and executes it as an ordinary subprocess. The command
shown on screen is the command that ran -- not a reconstruction of it -- which
is what makes "reproducible from the CLI" checkable rather than claimed.

Two consequences follow.

* Parameters are never restated here. `build_parser()` is introspected for
  names, defaults, types and `choices`; `PARAM_META` and `FLUID_ARGS` supply
  the units and ranges argparse has no field for. Adding an argument to the
  runner makes it appear in the form with no change to the UI.
* Progress leaves the process through `--status-json`, an atomic file the
  runner rewrites from its existing `on_step` hook. Nothing about the numerics
  is touched to make the monitor work.

The child is started in its own session, so it survives a Streamlit rerun, a
browser refresh and the server being restarted. A run outliving its launcher is
the normal case: a K-GEP-1 job is hours.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
if os.path.join(ROOT, "scripts") not in sys.path:
    sys.path.insert(0, os.path.join(ROOT, "scripts"))

import kgep1_run as runner                             # noqa: E402

RUNS_DIR = os.path.join("output", "ui_runs")


# --------------------------------------------------------------------------
# the parameter surface, read from the runner
# --------------------------------------------------------------------------
def parameters() -> list[dict]:
    """Every settable run parameter, as the runner defines it.

    One dict per argument: `dest`, `flag`, `type`, `default`, `choices`,
    `help`, and `unit`/`lo`/`hi` where `PARAM_META` or `FLUID_ARGS` supply them.
    Arguments that are not part of the physical case -- `--status-json`,
    `--no-resume` -- are left out; the launcher sets those itself.
    """
    fluid_meta = {dest: unit for _, dest, unit, _ in runner.FLUID_ARGS}
    skip = {"help", "status_json", "no_resume"}
    out = []
    for a in runner.build_parser()._actions:
        if a.dest in skip:
            continue
        meta = dict(runner.PARAM_META.get(a.dest, {}))
        if a.dest in fluid_meta:
            meta.setdefault("unit", fluid_meta[a.dest])
        out.append({
            "dest": a.dest,
            "flag": a.option_strings[-1],
            "type": a.type,
            "default": a.default,
            "choices": list(a.choices) if a.choices else None,
            "help": a.help or "",
            "unit": meta.get("unit"),
            "lo": meta.get("lo"),
            "hi": meta.get("hi"),
        })
    return out


def argv_from(values: dict, status_json: str | None = None) -> list[str]:
    """The full command for a case.

    EVERY parameter is written out, including the ones left at their default.
    A shorter command that relies on defaults reproduces the run only while
    those defaults hold; this one keeps working after `config.py` moves, which
    is the situation the reproducibility requirement exists for.
    """
    argv = [sys.executable, "scripts/kgep1_run.py"]
    for p in parameters():
        argv += [p["flag"], str(values.get(p["dest"], p["default"]))]
    if status_json:
        argv += ["--status-json", status_json]
    return argv


def shell_command(argv: list[str]) -> str:
    """The same command as a copy-pasteable line."""
    body = " ".join(["python" if a == sys.executable else a for a in argv])
    return "PYTHONPATH=. " + body


# --------------------------------------------------------------------------
# launching and watching
# --------------------------------------------------------------------------
def launch(values: dict, tag: str | None = None) -> dict:
    """Start a run detached. Returns the handle the monitor reads."""
    os.makedirs(RUNS_DIR, exist_ok=True)
    tag = tag or time.strftime("%Y%m%d-%H%M%S")
    status_path = os.path.join(RUNS_DIR, f"{tag}.status.json")
    log_path = os.path.join(RUNS_DIR, f"{tag}.log")
    argv = argv_from(values, status_json=status_path)

    env = dict(os.environ, PYTHONPATH=ROOT)
    log = open(log_path, "w")
    proc = subprocess.Popen(argv, cwd=ROOT, env=env, stdout=log,
                            stderr=subprocess.STDOUT, start_new_session=True)
    handle = {"tag": tag, "pid": proc.pid, "status": status_path,
              "log": log_path, "command": shell_command(argv),
              "started_at": time.time()}
    with open(os.path.join(RUNS_DIR, f"{tag}.handle.json"), "w") as f:
        json.dump(handle, f)
    return handle


def handles() -> list[dict]:
    """Every run this interface has started, newest first."""
    if not os.path.isdir(RUNS_DIR):
        return []
    out = []
    for name in os.listdir(RUNS_DIR):
        if name.endswith(".handle.json"):
            try:
                with open(os.path.join(RUNS_DIR, name)) as f:
                    out.append(json.load(f))
            except (OSError, ValueError):
                continue
    return sorted(out, key=lambda h: h.get("started_at", 0), reverse=True)


def read_status(path: str) -> dict | None:
    """The child's last atomic write, or None if it has not written yet.

    A partial read cannot happen -- the runner writes a temp file and renames --
    but the file may legitimately not exist for the first minutes of a run,
    while the closure table builds. The caller must show that as "starting",
    not as a failure.
    """
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def stop(pid: int) -> None:
    """Terminate a run.

    There is no pause: the solver has no such state, and a fake one would be a
    lie about what the button does. Stopping is safe rather than destructive --
    the run checkpoints every 180 s, so relaunching the same command resumes
    from the last checkpoint.
    """
    try:
        os.killpg(os.getpgid(pid), signal.SIGTERM)
    except OSError:
        pass


def tail(path: str, n: int = 40) -> str:
    try:
        with open(path) as f:
            return "".join(f.readlines()[-n:])
    except OSError:
        return ""
