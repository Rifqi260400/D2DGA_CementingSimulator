"""Run the test suite and write a machine-readable verdict.

Why this exists (UI gap G-7): the Validation screen has to show whether the
gates pass, and `docs/gate_status.md` is prose written by hand. Prose cannot go
stale loudly. pytest's own summary is a line of stdout that is easy to lose --
during the remediation it was lost twice on detached runs, and "I believe it
passed" is not a gate.

This is a pytest PLUGIN, not a second test runner: it adds an observer and
changes nothing about collection, selection or assertion. No test is modified,
skipped or reordered by it.

    PYTHONPATH=. python tests/run_gates.py [output.json]

The default output is `output/gates.json`, which the Validation screen reads.
A screen that finds no such file says the gates have not been run here -- it
does not fall back to the prose.
"""

from __future__ import annotations

import json
import pathlib
import platform
import subprocess
import sys
import time

import pytest

DEFAULT_OUT = pathlib.Path("output/gates.json")


def _git(*args):
    try:
        return subprocess.run(("git",) + args, capture_output=True, text=True,
                              timeout=10).stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


class Verdict:
    """Records the outcome of every test, and writes it once at the end."""

    def __init__(self, out: pathlib.Path):
        self.out = out
        self.counts: dict[str, int] = {}
        self.failures: list[dict] = []
        self.modules: dict[str, dict] = {}
        self.t0 = time.time()

    def pytest_runtest_logreport(self, report):
        # `call` is the test body; a non-passing `setup` is a failure or a skip
        # that never reached the body, and both matter.
        if report.when != "call" and not (report.when == "setup"
                                          and report.outcome != "passed"):
            return
        self.counts[report.outcome] = self.counts.get(report.outcome, 0) + 1
        module = report.nodeid.split("::", 1)[0]
        m = self.modules.setdefault(module, {})
        m[report.outcome] = m.get(report.outcome, 0) + 1
        if report.outcome == "failed":
            self.failures.append({
                "nodeid": report.nodeid,
                "message": str(report.longrepr)[-2000:],
            })

    def pytest_sessionfinish(self, session, exitstatus):
        self.out.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.out.with_suffix(self.out.suffix + ".tmp")
        tmp.write_text(json.dumps({
            "schema": 1,
            "exit": int(exitstatus),
            "counts": self.counts,
            "modules": self.modules,
            "failures": self.failures,
            "duration_s": time.time() - self.t0,
            "finished_at": time.time(),
            "environment": {
                "git_sha": _git("rev-parse", "HEAD"),
                "git_dirty": bool(_git("status", "--porcelain")),
                "python": sys.version.split()[0],
                "platform": platform.platform(),
            },
        }, indent=2))
        tmp.replace(self.out)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    out = pathlib.Path(argv.pop(0)) if argv and not argv[0].startswith("-") \
        else DEFAULT_OUT
    verdict = Verdict(out)
    code = pytest.main(["-q", *argv], plugins=[verdict])
    print(f"\nwrote {out}")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
