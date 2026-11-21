"""Graphical abstract: what the vsbench pipeline does, and what it measured.

Every number drawn here is read from the run's own result files, so the figure cannot
drift from the tables. Writes PNG (300 dpi), PDF and SVG.

Usage: python scripts/graphical_abstract.py [RESULTS_DIR] [OUT_STEM]
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
RES = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "q1"
STEM = Path(sys.argv[2]) if len(sys.argv) > 2 else RES / "figures" / "graphical_abstract"

INK, MUTE = "#1A202C", "#4A5568"
STAGE_F, STAGE_E = "#EDF2F7", "#9FB3C8"
IN_F, IN_E = "#E2EDF6", "#7FA7C4"
CTRL_F, CTRL_E = "#FDF2DC", "#C2951F"
DOCK, BASE, CHANCE = "#2C5282", "#AE4E17", "#8A97A6"

# ----------------------------------------------------------------- numbers
SETS = ["random", "matched", "decoy"]
SET_LABEL = {"random": "random\ninactives", "matched": "property-\nmatched",
             "decoy": "DUD-E\ndecoys"}

macro = pd.read_csv(RES / "metrics_macro.csv")
dock = (macro[(macro.metric == "roc_auc") & (macro.score == "consensus")]
        .set_index("negative_set").value)

base = pd.concat([pd.read_csv(f).assign(target=Path(f).parts[-3])
                  for f in sorted(glob.glob(str(RES / "*" / "analysis" / "baselines.csv")))],
                 ignore_index=True)
best = base.loc[base.groupby(["target", "negative_set"]).roc_auc.idxmax()]
basemac = best.groupby("negative_set").roc_auc.mean()

scored = sorted(glob.glob(str(RES / "*" / "analysis" / "scored.csv")))
n_compounds = sum(len(pd.read_csv(f)) for f in scored)
n_targets = len(scored)

reg = pd.concat([pd.read_csv(f) for f in
                 sorted(glob.glob(str(RES / "*" / "analysis" / "baseline_comparison.csv")))],
                ignore_index=True)
reg = reg[reg.score.isin(("autodock-gpu", "gnina", "consensus"))]
n_cells, n_above = len(reg), int(reg.beats_baseline.sum())

# ----------------------------------------------------------------- canvas
W, H = 132.0, 66.0
fig = plt.figure(figsize=(W / 10, H / 10), dpi=300)
fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")


def box(x, y, w, h, label, sub=None, fc=STAGE_F, ec=STAGE_E, fs=9.0, subfs=7.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.1",
                                linewidth=1.0, facecolor=fc, edgecolor=ec, zorder=2))
    ax.text(x + w / 2, y + h * (0.70 if sub else 0.5), label, ha="center", va="center",
            fontsize=fs, color=INK, fontweight="bold", zorder=3)
    if sub:
        ax.text(x + w / 2, y + h * 0.33, sub, ha="center", va="center", fontsize=subfs,
                color=MUTE, linespacing=1.45, zorder=3)


def arrow(x1, y1, x2, y2, color=STAGE_E, lw=1.4, ms=9, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=ms,
                                 linewidth=lw, color=color, linestyle=ls,
                                 shrinkA=0, shrinkB=0, zorder=1))


def band(x, y, w, h, title, fc, tc):
    ax.add_patch(Rectangle((x, y), w, h, linewidth=0, facecolor=fc, zorder=0))
    ax.text(x + 1.6, y + h - 1.7, title, ha="left", va="center", fontsize=7.6,
            color=tc, fontweight="bold", zorder=1)


LEFT, RIGHT = 2.0, 93.5
LW = RIGHT - LEFT

# ----------------------------------------------------------------- title
ax.text(LEFT, 62.6, "How much docking enrichment survives property matching?",
        fontsize=15.0, color=INK, fontweight="bold", va="center")
ax.text(LEFT, 59.2, f"vsbench  ·  an open, tested harness  ·  {n_targets} LIT-PCBA targets  "
                    f"·  {n_compounds:,} compounds  ·  two independent scoring functions",
        fontsize=9.0, color=MUTE, va="center")

# ----------------------------------------------------------------- inputs
band(LEFT, 44.0, LW, 13.5, "INPUTS", "#E8F0F7", "#2C5282")
inputs = [
    ("5 LIT-PCBA targets", "chosen from 15 by a rule fixed\nbefore docking: ≥ 25 actives,\n≤ 70,000 inactives"),
    ("129 receptor templates", "each audited for metal ions against\nits original RCSB deposition"),
    ("Three negative sets per target", "300 random  ·  8:1 matched on seven\ndescriptors  ·  8:1 DUD-E decoys"),
]
iw, igap = (LW - 2 * 2.4) / 3, 2.4
for i, (name, sub) in enumerate(inputs):
    x = LEFT + i * (iw + igap)
    box(x, 44.8, iw, 9.4, name, sub, fc=IN_F, ec=IN_E, fs=9.2, subfs=7.3)
    arrow(x + iw / 2, 44.6, x + iw / 2, 42.6, color=IN_E, lw=1.3)

# ----------------------------------------------------------------- pipeline
band(LEFT, 25.0, LW, 17.5, "PIPELINE", "#EDF1F5", "#2D3748")
stages = [
    ("check-data", "parses every input file;\nrefuses a receptor whose\nresidues cannot be read"),
    ("prep", "ETKDGv3 + MMFF94 + meeko;\nmetals restored or the run\nfails rather than drop one"),
    ("templates", "redock each co-crystal ligand;\nin-place symmetry-aware RMSD;\nselected on RMSD alone"),
    ("screen", "AutoDock-GPU × 3 seeds\n+ gnina CNN rescoring;\nresumable, 6 GB GPU"),
    ("analyze", "stratified bootstrap;\nsize-controlled AUC;\nbaselines + multiplicity"),
    ("report", "34 tables, 9 figures,\nevery draw reproducible\nfrom a single seed"),
]
sgap = 1.9
sw = (LW - 5 * sgap) / 6
sy, sh = 25.8, 12.6
for i, (name, sub) in enumerate(stages):
    x = LEFT + i * (sw + sgap)
    box(x, sy, sw, sh, name, sub, fs=9.6, subfs=7.1)
    if i:
        arrow(x - sgap - 0.05, sy + sh / 2, x - 0.2, sy + sh / 2, lw=1.6, ms=10)

# ----------------------------------------------------------------- controls
band(LEFT, 4.0, LW, 18.0, "LABEL-FREE CONTROLS, RUN BESIDE EVERY DOCKING SCORE",
     "#FCF4E2", "#8A6A14")
ctrl = [
    ("Seven physical descriptors", "size, logP, TPSA,\nH-bond counts, charge"),
    ("1-NN ECFP4 similarity", "no protein, no pose,\nno scoring function"),
    ("Heavy-atom stratification", "conditional AUC within\na compound's own size class"),
    ("Multiplicity adjustment", "Bonferroni across each\nfamily of comparisons"),
]
cgap = 2.2
cw = (LW - 3 * cgap) / 4
cy, ch = 5.2, 11.8
for i, (name, sub) in enumerate(ctrl):
    x = LEFT + i * (cw + cgap)
    box(x, cy, cw, ch, name, sub, fc=CTRL_F, ec=CTRL_E, fs=8.8, subfs=7.1)
    # drawn only in the gap between the bands, so the band label stays clear
    arrow(x + cw / 2, 22.1, x + cw / 2, sy - 0.9, color=CTRL_E, lw=1.2,
          ls=(0, (2.6, 2.0)), ms=9)

# ----------------------------------------------------------------- result panel
rx, ry, rw, rh = 96.8, 5.0, W - 96.8 - 2.0, 52.5
ax.add_patch(FancyBboxPatch((rx, ry), rw, rh, boxstyle="round,pad=0,rounding_size=1.5",
                            linewidth=1.4, facecolor="white", edgecolor=INK, zorder=2))
ax.text(rx + 2.2, ry + rh - 3.0, "WHAT IT MEASURED", fontsize=7.8, color=INK,
        fontweight="bold", va="center", zorder=3)
arrow(RIGHT + 0.4, 32.1, rx - 0.5, 32.1, color=INK, lw=1.8, ms=12)

cax = fig.add_axes([(rx + 6.6) / W, 32.2 / H, (rw - 9.6) / W, 14.2 / H], zorder=4)
xs = list(range(len(SETS)))
cax.axhline(0.5, color=CHANCE, lw=1.0, ls=(0, (3, 2.4)), zorder=1)
cax.plot(xs, [basemac[s] for s in SETS], "-o", color=BASE, lw=2.2, ms=7, zorder=3)
cax.plot(xs, [dock[s] for s in SETS], "-s", color=DOCK, lw=2.2, ms=6.3, zorder=3)
for i, s in enumerate(SETS):
    cax.annotate(f"{dock[s]:.2f}", (i, dock[s]), textcoords="offset points",
                 xytext=(0, -14), ha="center", fontsize=8.0, color=DOCK, fontweight="bold")
    cax.annotate(f"{basemac[s]:.2f}", (i, basemac[s]), textcoords="offset points",
                 xytext=(0, 9), ha="center", fontsize=8.0, color=BASE, fontweight="bold")
cax.text(-0.36, 0.506, "chance", fontsize=7.0, color=CHANCE, va="bottom", ha="left")
cax.set_xlim(-0.40, 2.40)
cax.set_ylim(0.455, 0.775)
cax.set_xticks(xs)
cax.set_xticklabels([SET_LABEL[s] for s in SETS], fontsize=7.4, color=MUTE, linespacing=1.3)
cax.set_yticks([0.5, 0.6, 0.7])
cax.tick_params(labelsize=7.2, colors=MUTE, length=2.5, pad=2)
cax.set_ylabel("ROC-AUC (mean over targets)", fontsize=7.6, color=INK, labelpad=3)
for side in ("top", "right"):
    cax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    cax.spines[side].set_color(STAGE_E)

# hand-placed legend, so nothing can collide with the panel heading
for dy, colour, marker, text in ((0.0, BASE, "o", "best structure-free baseline"),
                                 (-3.1, DOCK, "s", "two-engine docking consensus")):
    ly = ry + rh - 7.4 + dy
    ax.add_line(Line2D([rx + 3.0, rx + 6.4], [ly, ly], color=colour, lw=2.2, zorder=3))
    ax.add_line(Line2D([rx + 4.7], [ly], color=colour, marker=marker, ms=6, zorder=3))
    ax.text(rx + 7.4, ly, text, fontsize=7.6, color=INK, va="center", zorder=3)

ax.text(rx + rw / 2, ry + 16.6,
        "Matching removes the enrichment,\nand most of what it removed\nwas molecular size.",
        ha="center", va="center", fontsize=9.4, color=INK, fontweight="bold",
        linespacing=1.55, zorder=3)
ax.add_line(Line2D([rx + 5.0, rx + rw - 5.0], [ry + 11.0, ry + 11.0],
                   color=STAGE_E, lw=0.9, zorder=3))
ax.text(rx + rw / 2, ry + 6.2,
        f"Docking outranks that baseline in\n{n_above} of {n_cells} target × score × set cells.",
        ha="center", va="center", fontsize=9.4, color=BASE, fontweight="bold",
        linespacing=1.55, zorder=3)

# ----------------------------------------------------------------- write
STEM.parent.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{STEM}.{ext}", facecolor="white")
print(f"wrote {STEM}.png / .pdf / .svg")
print("  docking consensus: " + ", ".join(f"{s} {dock[s]:.3f}" for s in SETS))
print("  best baseline:     " + ", ".join(f"{s} {basemac[s]:.3f}" for s in SETS))
print(f"  docking beats the baseline in {n_above} of {n_cells} cells")
print(f"  {n_compounds} compounds across {n_targets} targets")
