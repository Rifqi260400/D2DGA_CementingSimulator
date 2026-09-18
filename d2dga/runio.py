"""
Run provenance -- the settings that produced a result, written beside it.

Why this module exists
----------------------
A finished run used to leave exactly one structured artefact: the checkpoint
`.npz`.  Its FILENAME carries six settings (wall, inflow, n_phi, n_xi, w0,
volumes); the run has about twenty.  Two consequences, and the first is a
correctness hazard, not a bookkeeping nuisance:

  * `Simulation.run` resumes whenever the checkpoint path exists and the grid
    SHAPE matches.  Fluid properties are in neither the path nor the shape, so
    editing `mud_density` and relaunching the same command **silently continued
    a concentration field computed with the old fluids**.  The only trace was a
    line reading `resumed from ...`.
  * a figure could not be traced back to the settings that made it, which a
    thesis needs.

What it does
------------
`write_config` puts the whole `Config` tree, the CLI arguments and the
environment (git SHA, library versions) in a JSON file beside the checkpoint.
`fingerprint` reduces the physics-and-mesh part of that to a short hash, which
`Simulation.run` stores inside the checkpoint and refuses to resume across --
see `checkpoint_tag` there.  `describe_mismatch` turns a refusal into a list of
the settings that actually changed, because "hashes differ" is not actionable.

It touches no numerics.  Nothing here is imported by the solver's physics
modules; the run SCRIPTS call it.

Layout, chosen for backward compatibility: the checkpoint path is unchanged, so
runs recorded before this module existed still resume.  The companions sit
beside it with the same stem.

    output/kgep1_ckpt_<...>.npz              the checkpoint, as before
    output/kgep1_ckpt_<...>.config.json      every setting + environment
    output/kgep1_ckpt_<...>.metrics.json     the post-processed results

A checkpoint with no `.config.json` is a legacy run: readable, with its
settings recorded as unknown rather than guessed.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from dataclasses import asdict

# Keys of the CLI-argument block that change the ANSWER, and so must not differ
# between a checkpoint and the run resuming it.  Everything else about a run --
# where it writes, how often it checkpoints, how loud it is -- may differ freely.
_PHYSICAL_ARG_KEYS = (
    "wall", "w0", "n_phi", "n_xi", "cfl", "volumes", "inflow", "n_c", "n_h",
)


# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
def _stem(checkpoint_path: str) -> str:
    return checkpoint_path[:-4] if checkpoint_path.endswith(".npz") else checkpoint_path


def config_path_for(checkpoint_path: str) -> str:
    return _stem(checkpoint_path) + ".config.json"


def metrics_path_for(checkpoint_path: str) -> str:
    return _stem(checkpoint_path) + ".metrics.json"


# --------------------------------------------------------------------------
# the payload
# --------------------------------------------------------------------------
def environment() -> dict:
    """Git SHA and library versions -- the other half of reproducibility."""
    def _git(*args):
        try:
            return subprocess.run(("git",) + args, capture_output=True, text=True,
                                  timeout=10).stdout.strip() or None
        except (OSError, subprocess.SubprocessError):
            return None

    import numpy
    import scipy
    return {
        "git_sha": _git("rev-parse", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
    }


def config_payload(cfg, args: dict, command: list[str] | None = None) -> dict:
    """
    The full record of a run.

    `cfg` is a `d2dga.config.Config`; `args` the run script's own arguments.
    `physical` is the subset that determines the answer and is what
    `fingerprint` hashes -- kept as its own block so the distinction is visible
    in the file rather than living only in this module.
    """
    conf = asdict(cfg)
    return {
        "schema": 1,
        "config": conf,
        "args": dict(args),
        "physical": {
            "config": conf,
            "args": {k: args[k] for k in _PHYSICAL_ARG_KEYS if k in args},
        },
        "command": list(command) if command is not None else list(sys.argv),
        "environment": environment(),
    }


def fingerprint(payload: dict) -> str:
    """
    Short hash of the physics-and-mesh part of a payload.

    Canonical JSON (sorted keys, no whitespace) so the hash depends on the
    values and not on dict ordering.  Deliberately does NOT cover the command
    line, the environment, or output paths: re-running the same physics from a
    different directory, or after a docs commit, must still resume.
    """
    physical = payload.get("physical", payload)
    blob = json.dumps(physical, sort_keys=True, separators=(",", ":"),
                      default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def _normalise(obj):
    """Round-trip through JSON so a tuple read back as a list does not read as a
    change.  The FINGERPRINT never had this problem -- a tuple and a list
    serialise identically -- but the human-readable diff did, and reported
    `mean_velocities_m_s: [0.05, 0.2, 0.5] -> (0.05, 0.2, 0.5)` as a
    difference."""
    return json.loads(json.dumps(obj, default=str))


def describe_mismatch(old: dict, new: dict, prefix: str = "") -> list[str]:
    """
    Flat list of the leaf settings that differ, as `path: old -> new`.

    This is the difference between a usable refusal and an unusable one: the
    hash says two runs differ, this says the mud density changed.
    """
    if not prefix:
        old, new = _normalise(old), _normalise(new)
    out = []
    for key in sorted(set(old) | set(new)):
        a, b = old.get(key, "<absent>"), new.get(key, "<absent>")
        path = f"{prefix}{key}"
        if isinstance(a, dict) and isinstance(b, dict):
            out.extend(describe_mismatch(a, b, prefix=path + "."))
        elif a != b:
            out.append(f"{path}: {a!r} -> {b!r}")
    return out


# --------------------------------------------------------------------------
# read / write
# --------------------------------------------------------------------------
def _atomic_write_json(path: str, obj: dict) -> str:
    """Temp file then `os.replace`, the same pattern `Simulation.run` uses for
    its checkpoint: a reclaim mid-write cannot leave a truncated record."""
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=False, default=str)
        fh.write("\n")
    os.replace(tmp, path)
    return path


def write_config(checkpoint_path: str, payload: dict) -> str:
    return _atomic_write_json(config_path_for(checkpoint_path), payload)


def read_config(checkpoint_path: str) -> dict | None:
    """The payload, or None for a run recorded before this module existed."""
    path = config_path_for(checkpoint_path)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def write_metrics(checkpoint_path: str, metrics: dict) -> str:
    return _atomic_write_json(metrics_path_for(checkpoint_path), metrics)


def read_metrics(checkpoint_path: str) -> dict | None:
    path = metrics_path_for(checkpoint_path)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


# --------------------------------------------------------------------------
# the guard
# --------------------------------------------------------------------------
def check_resume(checkpoint_path: str, payload: dict):
    """
    Called BEFORE `Simulation.run`, to fail early and readably.

    Returns `(level, message)` with level in `{None, "warn", "block"}`:

      None    no checkpoint, or one whose settings match -- resume freely.
      "warn"  a LEGACY checkpoint with no recorded config.  Its settings cannot
              be verified, but refusing would break relaunching a run finished
              before this module existed, so it proceeds noisily.
      "block" the recorded settings differ.  Resuming would continue a field
              computed under other physics.

    `Simulation.run` enforces the same thing through `checkpoint_tag`; this
    exists so the user is told *which* setting moved instead of two hashes.
    """
    if not os.path.exists(checkpoint_path):
        return None, ""
    old = read_config(checkpoint_path)
    if old is None:
        return "warn", (
            f"{checkpoint_path} is a LEGACY checkpoint with no recorded config, "
            f"so this resume cannot be verified against the current settings. "
            f"Pass --no-resume to start clean."
        )
    if fingerprint(old) == fingerprint(payload):
        return None, ""
    diffs = describe_mismatch(old.get("physical", {}), payload.get("physical", {}))
    listed = "\n".join("    " + d for d in diffs[:20])
    more = f"\n    ... and {len(diffs) - 20} more" if len(diffs) > 20 else ""
    return "block", (
        f"{checkpoint_path} was written with DIFFERENT settings, so resuming it "
        f"would continue a field computed under other physics:\n"
        f"{listed}{more}\n"
        f"  Pass --no-resume to start over, or use a different output path."
    )
