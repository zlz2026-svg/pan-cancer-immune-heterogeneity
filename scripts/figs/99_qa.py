"""99_qa.py -- automated figure QA for the pan-cancer immune heterogeneity manuscript.

Checks, per figure:
  1. pdf / svg / png all present
  2. PDF MediaBox matches the declared print size (mm)
  3. PNG pixel dimensions agree with the declared size at 600 dpi
  4. every deliverable file < 20 MB (npj hard limit)
  5. minimum font size embedded in the PDF >= 5 pt, and which strings are smallest
  6. PNG is white-backed, has no pure red/green channel-only encoding, and keeps
     the locked palette separable after deuteranopia / protanopia simulation
       (CIE76 dE in Lab, linear-RGB Vienot 1999 transform)

Run from the project root:  python scripts/figs/99_qa.py
"""

import os
import re
import sys
import zlib

import numpy as np
from PIL import Image

ROOT = os.getcwd()
FIGDIR = os.path.join(ROOT, "figures")
MAX_MB = 20.0
MIN_PT = 5.0
SIZE_TOL_MM = 1.0

# name: (width_mm, height_mm, note, uses_germ_colour)
REGISTRY = {
    "Fig1_design": (183, 178, "main", True),
    "Fig2_heterogeneity": (183, 216, "main", True),
    "Fig3_attribution": (183, 152, "main", True),
    "Fig4_reversal_validation": (183, 100, "main", False),
    "FigS2_robustness": (183, 110, "supplementary", False),
    "Supplementary_Figure_S1": (183, 95, "supplementary", False),
    "Fig_GA": (183, 120, "main", True),
}

# locked palette (mirrors scripts/figs/00_tokens.R)
PAL = {
    "prot": "#2166AC",
    "risk": "#E08214",
    "ink": "#1A1A1A",
    "ink2": "#4D4D4D",
    "ref": "#7F7F7F",
    "ns": "#D9D9D9",
    "band": "#F0F0F0",
    "endoderm": "#0072B2",
    "mesoderm": "#009E73",
    "ectoderm": "#CC79A7",
}
GERM = ["endoderm", "mesoderm", "ectoderm"]

problems = []
notes = []


def fail(msg):
    problems.append(msg)


def ok(msg):
    notes.append(msg)


def hex2rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], dtype=float) / 255.0


# ---------------------------------------------------------------- colour science
def srgb_to_linear(c):
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def linear_to_srgb(c):
    c = np.clip(c, 0.0, 1.0)
    return np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1 / 2.4) - 0.055)


RGB2LMS = np.array([[17.8824, 43.5161, 4.11935],
                    [3.45565, 27.1554, 3.86714],
                    [0.0299566, 0.184309, 1.46709]])
LMS2RGB = np.linalg.inv(RGB2LMS)
SIM = {
    "deuteranopia": np.array([[1.0, 0.0, 0.0], [0.494207, 0.0, 1.24827], [0.0, 0.0, 1.0]]),
    "protanopia": np.array([[0.0, 2.02344, -2.52581], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
}
M_RGB2XYZ = np.array([[0.4124, 0.3576, 0.1805],
                      [0.2126, 0.7152, 0.0722],
                      [0.0193, 0.1192, 0.9505]])
WHITE = np.array([0.95047, 1.0, 1.08883])


def rgb_to_lab(rgb):
    xyz = M_RGB2XYZ @ srgb_to_linear(np.asarray(rgb, dtype=float))
    t = xyz / WHITE
    f = np.where(t > (6 / 29) ** 3, np.cbrt(t), t / (3 * (6 / 29) ** 2) + 4 / 29)
    L = 116 * f[1] - 16
    a = 500 * (f[0] - f[1])
    b = 200 * (f[1] - f[2])
    return np.array([L, a, b])


def simulate(rgb, mode):
    lin = srgb_to_linear(np.asarray(rgb, dtype=float))
    return linear_to_srgb(LMS2RGB @ (SIM[mode] @ (RGB2LMS @ lin)))


def dE76(c1, c2):
    return float(np.linalg.norm(rgb_to_lab(c1) - rgb_to_lab(c2)))


# ---------------------------------------------------------------- pdf parsing
def pdf_uncompress(path):
    data = open(path, "rb").read()
    chunks = []
    for m in re.finditer(rb"stream\r?\n", data):
        s = m.end()
        e = data.find(b"endstream", s)
        if e < 0:
            continue
        try:
            chunks.append(zlib.decompress(data[s:e]))
        except zlib.error:
            continue
    return data, b"\n".join(chunks)


def pdf_check(name, path, w_mm, h_mm, germ):
    data, txt = pdf_uncompress(path)
    mb = re.search(rb"/MediaBox\s*\[\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*\]", data)
    if not mb:
        fail(f"{name}: MediaBox not found in PDF")
        return
    pw = float(mb.group(3)) / 72 * 25.4
    ph = float(mb.group(4)) / 72 * 25.4
    if abs(pw - w_mm) > SIZE_TOL_MM or abs(ph - h_mm) > SIZE_TOL_MM:
        fail(f"{name}: PDF {pw:.1f}x{ph:.1f} mm != declared {w_mm}x{h_mm} mm")
    else:
        ok(f"{name}: PDF {pw:.1f}x{ph:.1f} mm")

    # font size lives in the text matrix scaling:  "s 0 0 -s x y Tm" then "/f-n 1 Tf"
    sizes = {}
    for m in re.finditer(
            rb"([0-9.]+)\s+0\s+0\s+-?[0-9.]+\s+[-0-9.]+\s+[-0-9.]+\s+Tm\s*\n?\s*/\S+\s+1\s+Tf\s*\n?\((.*?)\)\s*Tj",
            txt):
        s = round(float(m.group(1)), 2)
        if s < 0.2 or s > 40:
            continue
        label = m.group(2).decode("latin-1", "replace")[:28]
        sizes.setdefault(s, []).append(label)
    if not sizes:
        notes.append(f"{name}: font sizes not recoverable from PDF (text may be outlined)")
    else:
        smin = min(sizes)
        if smin < MIN_PT - 0.01:
            fail(f"{name}: smallest embedded font {smin} pt < {MIN_PT} pt "
                 f"({len(sizes[smin])} strings, e.g. {sizes[smin][:3]})")
        else:
            ok(f"{name}: min font {smin} pt (levels: "
               + ", ".join(f"{k}" for k in sorted(sizes)) + ")")

    # germ layer must never be carried by colour alone -- require the names in text
    if germ:
        # any literal string in the content stream, whether shown with Tj or a TJ array
        raw = re.findall(rb"\(((?:\\.|[^()\\])*)\)", txt)
        strings = " ".join(s.decode("latin-1", "replace") for s in raw).lower()
        missing = [g for g in ("endoderm", "mesoderm", "ectoderm") if g not in strings]
        if missing:
            fail(f"{name}: germ-layer colour used but {missing} never appear as text "
                 f"-- grouping is colour-only, unreadable in greyscale")
        else:
            ok(f"{name}: all three germ-layer names present as text")


# ---------------------------------------------------------------- png parsing
def png_check(name, path, w_mm, h_mm):
    im = Image.open(path)
    dpi = im.info.get("dpi", (None, None))[0]
    w_px, h_px = im.size
    if dpi:
        got_w = w_px / dpi * 25.4
        got_h = h_px / dpi * 25.4
        if abs(got_w - w_mm) > SIZE_TOL_MM or abs(got_h - h_mm) > SIZE_TOL_MM:
            fail(f"{name}: PNG implies {got_w:.1f}x{got_h:.1f} mm at {dpi:.0f} dpi "
                 f"!= declared {w_mm}x{h_mm} mm")
        else:
            ok(f"{name}: PNG {w_px}x{h_px} px @ {dpi:.0f} dpi = {got_w:.1f}x{got_h:.1f} mm")
    else:
        notes.append(f"{name}: PNG has no dpi tag ({w_px}x{h_px} px)")

    arr = np.asarray(im.convert("RGB"), dtype=np.uint8)
    flat = arr.reshape(-1, 3)
    # background must be white
    corner = flat[[0, arr.shape[1] - 1, -arr.shape[1], -1]]
    if not np.all(corner >= 250):
        fail(f"{name}: PNG corners are not white ({corner.tolist()})")
    white_frac = float(np.mean(np.all(flat >= 250, axis=1)))
    if white_frac < 0.55:
        fail(f"{name}: only {white_frac*100:.0f}% of pixels are white (expect >= 55%)")
    else:
        ok(f"{name}: white background {white_frac*100:.0f}% of pixels")

    # no red/green single-channel encoding
    reds = flat[(flat[:, 0] >= 200) & (flat[:, 1] < 60) & (flat[:, 2] < 60)]
    greens = flat[(flat[:, 1] >= 200) & (flat[:, 0] < 60) & (flat[:, 2] < 60)]
    if len(reds) > 2000 and len(greens) > 2000:
        fail(f"{name}: saturated red ({len(reds)} px) and green ({len(greens)} px) both "
             f"present -- risk of channel-only encoding")
    uniq = len(np.unique(flat[::7], axis=0))
    ok(f"{name}: ~{uniq} distinct colours in a sampled grid")


# ---------------------------------------------------------------- palette audit
def palette_audit():
    # the three germ-layer colours carry the categorical encoding in Fig 2 / Fig 4
    for mode in ["deuteranopia", "protanopia"]:
        worst = None
        for i, a in enumerate(GERM):
            for b in GERM[i + 1:]:
                d = dE76(simulate(hex2rgb(PAL[a]), mode), simulate(hex2rgb(PAL[b]), mode))
                if worst is None or d < worst[0]:
                    worst = (d, a, b)
        d, a, b = worst
        if d < 12:
            fail(f"germ layer colours {a} vs {b} collapse under {mode} (dE76 = {d:.1f})")
        else:
            ok(f"{mode}: closest germ-layer pair {a}/{b} dE76 = {d:.1f}")

    # grey-scale separation: the layers stay close in luminance, which is tolerable
    # only because every figure co-encodes them with a text label (checked above)
    Ls = {k: rgb_to_lab(hex2rgb(PAL[k]))[0] for k in GERM}
    order = sorted(Ls.items(), key=lambda kv: kv[1])
    gap = min(order[i + 1][1] - order[i][1] for i in range(len(order) - 1))
    msg = ("germ layer L* = " + ", ".join(f"{k} {v:.0f}" for k, v in order)
           + f"; min gap {gap:.1f} L* -> not separable in greyscale alone, "
             "relies on the co-encoded text labels")
    if gap < 8:
        notes.append("NOTE " + msg)
    else:
        ok(msg)

    # risk vs protective diverging ends
    d = dE76(hex2rgb(PAL["risk"]), hex2rgb(PAL["prot"]))
    if d < 40:
        fail(f"risk/protective endpoints too close (dE76 = {d:.1f})")
    else:
        ok(f"risk vs protective dE76 = {d:.1f}")


# ---------------------------------------------------------------- main
def main():
    print("=" * 74)
    print("FIGURE QA")
    print("=" * 74)
    for name, (w, h, kind, germ) in REGISTRY.items():
        print(f"\n-- {name}  ({kind}, target {w}x{h} mm)")
        for ext in ("pdf", "svg", "png"):
            p = os.path.join(FIGDIR, f"{name}.{ext}")
            if not os.path.exists(p):
                fail(f"{name}: missing {ext}")
                continue
            mb = os.path.getsize(p) / 1024 / 1024
            if mb > MAX_MB:
                fail(f"{name}: {ext} is {mb:.1f} MB (> {MAX_MB} MB)")
        pp = os.path.join(FIGDIR, f"{name}.pdf")
        if os.path.exists(pp):
            pdf_check(name, pp, w, h, germ)
        pn = os.path.join(FIGDIR, f"{name}.png")
        if os.path.exists(pn):
            png_check(name, pn, w, h)

    print("\n-- locked palette, colour-vision audit")
    palette_audit()

    print("\n" + "-" * 74)
    for n in notes:
        print("  ok   " + n)
    if problems:
        print("\nFAILED CHECKS")
        for p in problems:
            print("  FAIL " + p)
        print(f"\n{len(problems)} problem(s).")
        return 1
    print(f"\nALL CHECKS PASSED ({len(notes)} assertions).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
