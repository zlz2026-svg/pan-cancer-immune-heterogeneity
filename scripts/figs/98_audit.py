#!/usr/bin/env python
"""98_audit.py -- geometry + typography audit for the six manuscript figures.

Why this exists next to 99_qa.py
--------------------------------
99_qa.py owns the *journal contract* (physical size, embedded font floor, white
background, colour-vision separation, the three germ-layer words).  This script owns
the *rendered geometry*: it reads the final PDF and checks, per glyph box, that
nothing overlaps, nothing leaves the page, no string is drawn twice, and no glyph
run falls below the 5 pt floor.

Relationship to nature-figure's auditors
----------------------------------------
`audit_pdf_text.py` is used directly -- it is a content-stream scan and agrees with
this script.

`audit_figure_collisions.py` is NOT used as a gate.  It groups glyphs by PDF text
operator, so on cairo_pdf output one line of tightly spaced labels collapses into a
single wide box: on Fig2 it reported 77 failures whose small-box overlap was 98-100 %
between *different* labels drawn at different positions.  A span-level check with
PyMuPDF finds 0 real overlaps and 0 glyph runs outside the page on the same file, so
the operator-level grouping is the false-positive source, not the figure.  We keep the
span-level test (exact glyph boxes) and record the operator-level verdict for
reference only.

Exit code is non-zero on any failure, so run_all.R stops the chain.
"""

import itertools
import json
import os
import subprocess
import sys

import pymupdf

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIGDIR = os.path.join(ROOT, "figures")
QADIR = os.path.join(FIGDIR, "qa")
SKILL = r"C:\Users\zlz\.workbuddy\skills\nature-figure\scripts"

# name -> (width_mm, height_mm)
REGISTRY = {
    "Fig1_design": (183, 178),
    "Fig2_heterogeneity": (183, 216),
    "Fig3_attribution": (183, 152),
    "Fig4_reversal_validation": (183, 100),
    "FigS2_robustness": (183, 110),
    "Supplementary_Figure_S1": (183, 95),
}

SIZE_FLOOR_PT = 4.99      # journal floor is 5 pt; leave a hair for float rounding
SIZE_EDGE_MM = 1.0        # allowed bleed outside the declared page box
OVERLAP_LIMIT = 0.05      # share of the smaller box
DUP_LIMIT_PT = 0.5        # same string within this distance = drawn twice

# Known, registered technical debt.  A finding that matches gets reported as a note
# instead of a blocker, and the entry carries the batch that removes it.  Delete the
# entry when that batch lands -- an entry outliving its fix would hide a regression.
#
# Empty since batch B2: the Fig3 heatmap became full width with horizontal row
# labels, so the 43 % label collision it was registered for is gone.  Do not add an
# entry here without naming the batch that will remove it.
KNOWN = {}

failures = []
notes = []
lines = []


def ok(msg):
    lines.append(("ok  ", msg))


def bad(msg):
    lines.append(("FAIL", msg))
    failures.append(msg)


def warn(msg):
    lines.append(("warn", msg))
    notes.append(msg)


def spans_of(page):
    out = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for s in line["spans"]:
                t = s["text"].strip()
                if t:
                    out.append((t, s["bbox"], s["size"]))
    return out


def overlap_ratio(a, b):
    x0, y0 = max(a[0], b[0]), max(a[1], b[1])
    x1, y1 = min(a[2], b[2]), min(a[3], b[3])
    if x1 <= x0 or y1 <= y0:
        return 0.0
    inter = (x1 - x0) * (y1 - y0)
    area_a = (a[2] - a[0]) * (a[3] - a[1])
    area_b = (b[2] - b[0]) * (b[3] - b[1])
    return inter / min(area_a, area_b)


def audit(name, path, w_mm, h_mm):
    doc = pymupdf.open(path)
    page = doc[0]
    pt = 25.4 / 72.0
    got_w, got_h = page.rect.width * pt, page.rect.height * pt

    if abs(got_w - w_mm) > 1.0 or abs(got_h - h_mm) > 1.0:
        bad(f"{name}: page {got_w:.1f}x{got_h:.1f} mm, expected {w_mm}x{h_mm}")
    else:
        ok(f"{name}: page {got_w:.1f}x{got_h:.1f} mm")

    sp = spans_of(page)
    if not sp:
        bad(f"{name}: no extracted text -- glyphs may be outlined")
        doc.close()
        return {"name": name, "spans": 0}

    # 1. glyph runs below the journal floor
    small = [t for t, _, s in sp if s < SIZE_FLOOR_PT]
    smear = sorted({round(s, 2) for _, _, s in sp})
    if small:
        bad(f"{name}: {len(small)} glyph run(s) below {SIZE_FLOOR_PT} pt, "
            f"e.g. {small[:3]} (levels: {smear})")
    else:
        ok(f"{name}: every glyph run >= {min(smear):.1f} pt (levels: {smear})")

    # 2. text outside the page box
    m = SIZE_EDGE_MM / pt
    out = [(t, b) for t, b, _ in sp
           if b[0] < -m or b[1] < -m or b[2] > page.rect.width + m
           or b[3] > page.rect.height + m]
    if out:
        bad(f"{name}: {len(out)} glyph run(s) outside the page, e.g. {out[0][0]!r}")
    else:
        ok(f"{name}: all {len(sp)} glyph runs inside the page")

    # 3. real glyph-box overlaps
    hits = []
    for (t1, b1, _), (t2, b2, _) in itertools.combinations(sp, 2):
        r = overlap_ratio(b1, b2)
        if r > OVERLAP_LIMIT:
            hits.append((round(r, 3), t1, t2))
    if hits:
        worst = max(hits)
        msg = (f"{name}: {len(hits)} overlapping glyph pair(s), worst {worst[0]:.0%} "
               f"between {worst[1]!r} and {worst[2]!r}")
        if name in KNOWN:
            warn(f"{msg} -- KNOWN: {KNOWN[name]}")
        else:
            bad(msg)
    else:
        ok(f"{name}: no overlapping glyph boxes among {len(sp)} runs")

    # 4. the same string drawn twice at the same place
    dup = []
    seen = {}
    for t, b, _ in sp:
        key = (t, round(b[0] / DUP_LIMIT_PT), round(b[1] / DUP_LIMIT_PT))
        if key in seen:
            dup.append(t)
        seen[key] = True
    if dup:
        bad(f"{name}: {len(dup)} string(s) drawn twice, e.g. {dup[:3]}")
    else:
        ok(f"{name}: no duplicated strings")

    # 5. colourbar-free reference: how much of the page is ink
    doc.close()
    return {"name": name, "spans": len(sp), "min_pt": min(smear),
            "outside": len(out), "overlaps": len(hits), "dupes": len(dup)}


def skill_pdf_text(name, path):
    """Content-stream font audit from nature-figure (reliable on cairo output)."""
    script = os.path.join(SKILL, "audit_pdf_text.py")
    if not os.path.exists(script):
        warn(f"{name}: nature-figure audit_pdf_text.py not found -- skipped")
        return None
    r = subprocess.run([sys.executable, script, "--min-pt", str(SIZE_FLOOR_PT), path],
                       capture_output=True, text=True)
    verdict = "PASS" if "verdict: PASS" in r.stdout else "FAIL"
    if verdict == "FAIL":
        bad(f"{name}: audit_pdf_text verdict FAIL")
    else:
        ok(f"{name}: audit_pdf_text PASS (stream-level)")
    return verdict


def skill_collisions_reference(name, path):
    """Operator-level collision audit, recorded but not gated (see module docstring)."""
    script = os.path.join(SKILL, "audit_figure_collisions.py")
    if not os.path.exists(script):
        return None
    out = os.path.join(QADIR, f"{name}.collisions.json")
    r = subprocess.run([sys.executable, script, "--json-out", out, path],
                       capture_output=True, text=True)
    verdict = "PASS" if "verdict: PASS" in r.stdout else "REVIEW"
    n_fail = ""
    for line in r.stdout.splitlines():
        if line.startswith("summary:"):
            n_fail = line.replace("summary:", "").strip()
    warn(f"{name}: operator-level collision audit = {verdict} ({n_fail}); "
         f"span-level check above is authoritative")
    return verdict


def main():
    os.makedirs(QADIR, exist_ok=True)
    reports = []
    for name, (w, h) in REGISTRY.items():
        p = os.path.join(FIGDIR, f"{name}.pdf")
        if not os.path.exists(p):
            bad(f"{name}: missing PDF")
            continue
        lines.append(("----", f"{name}"))
        reports.append(audit(name, p, w, h))
        skill_pdf_text(name, p)
        skill_collisions_reference(name, p)

    with open(os.path.join(QADIR, "98_audit.json"), "w", encoding="utf-8") as fh:
        json.dump({"figures": reports, "failures": failures, "notes": notes},
                  fh, indent=2)

    print("=" * 74)
    for tag, msg in lines:
        print(f"  {tag:4s} {msg}")
    print("=" * 74)
    if failures:
        print(f"FAILED ({len(failures)} problem(s))")
        return 1
    print(f"ALL GEOMETRY CHECKS PASSED ({len(reports)} figures, "
          f"{sum(r['spans'] for r in reports)} glyph runs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
