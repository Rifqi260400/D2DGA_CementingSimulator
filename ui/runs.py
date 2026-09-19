"""
Finding and loading completed runs.

Everything here reads; nothing computes physics.  Post-processing is
`d2dga.postprocess`, called, never reimplemented (build spec section 5).

What a run leaves behind, and what that costs the viewer
-------------------------------------------------------
Since R-2 a run writes three files with a common stem:

    output/kgep1_ckpt_<...>.npz           checkpoint: FINAL field + scalar history
    output/kgep1_ckpt_<...>.config.json   every setting, git SHA, versions
    output/kgep1_ckpt_<...>.metrics.json  the post-processed results

There are no intermediate SNAPSHOTS and no stream function, so a viewer can show
the final field and the whole scalar history and cannot show a time slider or a
similarity collapse.  Those panels are marked unavailable rather than faked.

A run recorded before R-2 has no `.config.json`.  Its settings are then
RECONSTRUCTED from the filename plus the current `config.py` defaults, which is
correct only if `config.py` has not changed since.  `Run.config_source` says
which of the two you are looking at, and the UI must show it.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field

import numpy as np

# `output/kgep1_ckpt_<wall>_<inflow>_<nphi>x<nxi>_w<w0>_v<volumes>.npz`
_NAME = re.compile(
    r"kgep1_ckpt_(?P<wall>synthetic|caliper)_"
    r"(?P<inflow>no_axial_gradient|uniform)_"
    r"(?P<n_phi>\d+)x(?P<n_xi>\d+)_w(?P<w0>[\d.]+)_v(?P<volumes>[\d.]+)\.npz$")

RECORDED = "recorded"
RECONSTRUCTED = "reconstructed"


@dataclass
class Run:
    """One completed (or checkpointed) run on disk."""

    checkpoint: str
    args: dict                       # wall, inflow, n_phi, n_xi, w0, volumes, ...
    config_source: str               # RECORDED | RECONSTRUCTED
    payload: dict | None = None      # the .config.json, when recorded
    metrics: dict | None = None      # the .metrics.json, when written
    _npz: dict = field(default_factory=dict, repr=False)

    # -- identity ----------------------------------------------------------
    @property
    def name(self) -> str:
        return os.path.basename(self.checkpoint)[len("kgep1_ckpt_"):-len(".npz")]

    @property
    def mtime(self) -> float:
        return os.path.getmtime(self.checkpoint)

    @property
    def git_sha(self) -> str | None:
        env = (self.payload or {}).get("environment") or {}
        sha = env.get("git_sha")
        return sha[:8] if sha else None

    @property
    def git_dirty(self) -> bool:
        return bool(((self.payload or {}).get("environment") or {}).get("git_dirty"))

    def command(self) -> str:
        """The CLI invocation that reproduces this run.

        Built from the recorded arguments when they exist.  For a reconstructed
        run only the six settings in the filename are known, so the command is
        complete only in so far as everything else was left at its default --
        which is exactly what `config_source` warns about.
        """
        a = self.args
        parts = ["PYTHONPATH=.", "python", "scripts/kgep1_run.py",
                 f"--wall {a['wall']}", f"--w0 {a['w0']}",
                 f"--n-phi {a['n_phi']}", f"--n-xi {a['n_xi']}",
                 f"--volumes {a['volumes']}", f"--inflow {a['inflow']}"]
        for key, flag in (("cfl", "--cfl"), ("n_c", "--n-c"), ("n_h", "--n-h")):
            if key in a:
                parts.append(f"{flag} {a[key]}")
        return " ".join(parts)

    # -- arrays ------------------------------------------------------------
    def _load(self):
        if not self._npz:
            with np.load(self.checkpoint) as z:
                self._npz = {k: z[k] for k in z.files}
        return self._npz

    @property
    def concentration(self) -> np.ndarray:
        """Final gap-averaged concentration, shape (n_phi, n_xi)."""
        return self._load()["c"]

    @property
    def t(self) -> float:
        return float(self._load()["t"])

    @property
    def steps(self) -> int:
        return int(self._load()["n"])

    @property
    def history(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """`(times, efficiency, outlet_concentration)`, one sample per step."""
        d = self._load()
        return d["times"], d["effs"], d["outlet"]

    @property
    def conservation_error(self) -> float:
        return float(self._load()["cons_err"])

    @property
    def settings_tag(self) -> str | None:
        d = self._load()
        return str(d["tag"]) if "tag" in d and str(d["tag"]) else None

    @property
    def picard_unconverged(self) -> int:
        return int(self._load().get("n_picard_unconverged", 0))

    # -- things a run does NOT contain, stated rather than guessed ----------
    @property
    def has_snapshots(self) -> bool:
        return False        # only the final field is written

    @property
    def has_stream_function(self) -> bool:
        return False


def discover(output_dir: str = "output") -> list[Run]:
    """Every checkpoint in `output_dir`, newest first."""
    runs = []
    if not os.path.isdir(output_dir):
        return runs
    for fn in sorted(os.listdir(output_dir)):
        m = _NAME.match(fn)
        if not m:
            continue
        runs.append(load(os.path.join(output_dir, fn), _match=m))
    runs.sort(key=lambda r: r.mtime, reverse=True)
    return runs


def load(checkpoint: str, _match=None) -> Run:
    from d2dga import runio

    payload = runio.read_config(checkpoint)
    metrics = runio.read_metrics(checkpoint)
    if payload is not None:
        return Run(checkpoint=checkpoint, args=dict(payload.get("args") or {}),
                   config_source=RECORDED, payload=payload, metrics=metrics)

    m = _match or _NAME.match(os.path.basename(checkpoint))
    if m is None:
        raise ValueError(f"{checkpoint} is not a recognised run checkpoint and "
                         f"has no recorded config, so its settings are unknown")
    g = m.groupdict()
    args = {"wall": g["wall"], "inflow": g["inflow"],
            "n_phi": int(g["n_phi"]), "n_xi": int(g["n_xi"]),
            "w0": float(g["w0"]), "volumes": float(g["volumes"])}
    return Run(checkpoint=checkpoint, args=args, config_source=RECONSTRUCTED,
               payload=None, metrics=metrics)


def build_config(run: Run):
    """
    Rebuild the `Config` a run used.

    From the recorded payload when there is one.  Otherwise from the current
    `config.py` defaults with the grid from the filename -- correct only if
    `config.py` has not changed since, which `Run.config_source` flags.
    """
    from d2dga.config import (Config, FluidsConfig, GridConfig, StandoffConfig,
                              SyntheticWallConfig, WellConfig)

    grid = GridConfig(int(run.args["n_phi"]), int(run.args["n_xi"]))
    if run.config_source == RECONSTRUCTED or not run.payload:
        return Config(grid=grid)
    c = run.payload["config"]
    return Config(well=WellConfig(**c["well"]),
                  standoff=StandoffConfig(**c["standoff"]),
                  synthetic_wall=SyntheticWallConfig(**c["synthetic_wall"]),
                  fluids=FluidsConfig.from_record(c["fluids"]),
                  grid=grid)


def build_geometry(run: Run):
    """The `Geometry` a run used, so `postprocess` can be called on its field.

    Delegates the wall construction to `scripts/kgep1_run.make_geometry`, which
    is the same code the run itself used -- the caliper path in particular has
    QC arguments that must not be restated here.
    """
    import sys
    if "scripts" not in sys.path:
        sys.path.insert(0, "scripts")
    import kgep1_run as runner
    return runner.make_geometry(build_config(run), run.args["wall"])


def build_scaling(run: Run, geometry):
    """The `Scaling` a run used: b, m, B, Fr*, and the dimensional scales."""
    from d2dga.scaling import Scaling
    cfg = build_config(run)
    mud, cement = cfg.fluids.as_fluids()
    return Scaling(mud, cement, r_a_hat_star=geometry.r_a_hat_star,
                   delta_star=geometry.delta_star,
                   mean_velocity=float(run.args["w0"]))
