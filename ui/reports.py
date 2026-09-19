"""Reading the artefacts other parts of the repository write.

Kept out of the screens and out of `ui/common.py` (which is presentation) for
one reason: this is the code that decides what the Validation screen CLAIMS,
so it has to be testable without a browser. The first version of the ZF22
parser lived inside the screen, dropped the first data row because the
markdown separator was not where it assumed, and reported **9 / 9 reproduced**
on a table whose first row is the one case that is not reproduced. A screen
that quietly discards its own counter-evidence is worse than no screen, and
the defect was invisible until the page was looked at.

So: parse here, assert in `tests/test_ui_contract.py`, display there.
"""

from __future__ import annotations

import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ZF22 quotes eta_E to two significant figures, so agreement can only be
# claimed to about that resolution.  This band is the resolution of the
# COMPARISON, not a tolerance chosen so that cases pass -- at 10% case 1 misses
# by 45%, and no defensible band rescues it.
ZF22_BAND = 0.10


def read_blocked(path=None) -> list[dict]:
    """The items in `BLOCKED.md`, as `{id, title, kind}`.

    Parsed from the file rather than copied into the UI, so a screen cannot
    keep showing an item after it is resolved, or miss one that is added.
    `kind` is the section -- missing physical data, or irreducible ambiguity --
    which are different things to a reader deciding whether to trust a number.
    """
    path = path or os.path.join(ROOT, "BLOCKED.md")
    kind, out = "", []
    try:
        with open(path) as f:
            for line in f:
                if line.startswith("## "):
                    kind = ("missing physical data" if "5a" in line
                            else "irreducible ambiguity" if "5b" in line else "")
                elif line.startswith("### BLK-"):
                    head = line[4:].strip()
                    ident, _, title = head.partition(" — ")
                    out.append({"id": ident.strip(),
                                "title": (title or head).strip(),
                                "kind": kind})
    except OSError:
        return []
    return out


def parse_markdown_table(lines) -> tuple[list[str], list[list[str]]]:
    """`(header, rows)` from the first pipe table in `lines`.

    The separator row is identified by its CONTENT -- cells of dashes and
    colons -- not by its position. Assuming position is what dropped a data
    row: `output/zf22_table3.md` writes the separator with no space after the
    pipe, so a "starts with '| '" filter removes the separator from the list
    and then `[2:]` removes the first case as well.
    """
    def cells(line):
        return [c.strip() for c in line.strip().strip("|").split("|")]

    def is_separator(cs):
        return bool(cs) and all(set(c) <= set("-: ") and "-" in c for c in cs)

    table = [cells(ln) for ln in lines if ln.lstrip().startswith("|")]
    if not table:
        return [], []
    header, rest = table[0], table[1:]
    return header, [r for r in rest if not is_separator(r)]


def zf22_comparison(path=None) -> dict:
    """The ten published cases, each with a verdict.

    Returns `{header, rows, n_total, n_ok, failing, band}` where a row is
    `{cells, case, eta, zf22_eta, delta, ok}`. A row whose numbers cannot be
    read is `ok = False` and counted -- an unreadable row is not a pass.
    """
    path = path or os.path.join(ROOT, "output", "zf22_table3.md")
    try:
        with open(path) as f:
            lines = [ln.rstrip("\n") for ln in f]
    except OSError:
        return {"header": [], "rows": [], "n_total": 0, "n_ok": 0,
                "failing": [], "band": ZF22_BAND, "lines": []}

    header, body = parse_markdown_table(lines)
    col = {name: i for i, name in enumerate(header)}

    def num(cells, name):
        try:
            return float(cells[col[name]])
        except (KeyError, ValueError, IndexError):
            return float("nan")

    rows = []
    for cells in body:
        eta, ref = num(cells, "eta_E"), num(cells, "ZF22 eta_E")
        delta = (eta - ref) / ref if ref not in (0.0,) and ref == ref else float("nan")
        ok = bool(abs(eta - ref) <= ZF22_BAND * abs(ref)) if eta == eta and ref == ref \
            else False
        rows.append({"cells": cells, "case": cells[0] if cells else "?",
                     "eta": eta, "zf22_eta": ref, "delta": delta, "ok": ok})
    return {"header": header, "rows": rows, "n_total": len(rows),
            "n_ok": sum(r["ok"] for r in rows),
            "failing": [r["case"] for r in rows if not r["ok"]],
            "band": ZF22_BAND, "lines": lines}


def read_gates(path=None) -> dict | None:
    """`output/gates.json`, or None. Never falls back to the prose."""
    path = path or os.path.join(ROOT, "output", "gates.json")
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None
