#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fig S3: Single-cell validation of immune-signal dilution
a: 单细胞恶性比例 vs 免疫比例 (负相关, 稀释)
b: bulk Prolifer 签名 vs 单细胞 Cycling 比例 (正相关, bulk-单细胞增殖一致)
输出: results/sc_dilution/FigS3_singlecell_dilution.{pdf,png}
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import os

cmp = pd.read_csv("results/sc_dilution/sc_bulk_comparison.csv")
OUT = "results/sc_dilution"

# 配色（与正文一致：蓝/橙，色盲友好——避免红绿，拟合线用蓝）
SCATTER = "#2c3e50"      # 深蓝灰
FIT = "#1f77b4"          # 蓝（色盲友好，替代红）
TXT = "#1a1a1a"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
})

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))

# --- panel a: 恶性比例 vs 免疫比例 ---
ax = axes[0]
x, y = cmp["p_malignant"], cmp["p_immune"]
r, p = spearmanr(x, y)
ax.scatter(x, y, s=34, color=SCATTER, edgecolor="white", linewidth=0.5, zorder=3)
# 拟合线
z = np.polyfit(x, y, 1)
xs = np.linspace(x.min(), x.max(), 50)
ax.plot(xs, np.polyval(z, xs), color=FIT, lw=1.2, zorder=2)
ax.set_xlabel("Malignant-cell fraction (single-cell)")
ax.set_ylabel("Immune-cell fraction (single-cell)")
ax.set_xlim(-0.05, 0.85)
ax.set_ylim(-0.05, 1.0)
ax.text(0.03, 0.95, f"Spearman r = {r:.2f}\np = {p:.1e}",
        transform=ax.transAxes, va="top", ha="left", fontsize=8.5)
ax.set_title("a", loc="left", fontweight="bold")

# --- panel b: bulk Prolifer vs Cycling ---
ax = axes[1]
x, y = cmp["bulk_prolif"], cmp["p_cycling"]
r, p = spearmanr(x, y)
ax.scatter(x, y, s=34, color=SCATTER, edgecolor="white", linewidth=0.5, zorder=3)
z = np.polyfit(x, y, 1)
xs = np.linspace(x.min(), x.max(), 50)
ax.plot(xs, np.polyval(z, xs), color=FIT, lw=1.2, zorder=2)
ax.set_xlabel("Bulk proliferation score (per cancer)")
ax.set_ylabel("Cycling-cell fraction (single-cell)")
ax.set_xlim(2.0, 6.5)
ax.set_ylim(-0.02, 0.20)
ax.text(0.03, 0.95, f"Spearman r = {r:.2f}\np = {p:.3f}",
        transform=ax.transAxes, va="top", ha="left", fontsize=8.5)
ax.set_title("b", loc="left", fontweight="bold")

fig.tight_layout()
for ext in ["pdf", "png"]:
    fp = os.path.join(OUT, f"FigS3_singlecell_dilution.{ext}")
    fig.savefig(fp, dpi=300 if ext == "png" else None, bbox_inches="tight")
    print(f"已写出 {fp}")

print(f"\n面板 a: 恶性 vs 免疫, r={spearmanr(cmp['p_malignant'], cmp['p_immune'])[0]:.3f}")
print(f"面板 b: Prolifer vs Cycling, r={spearmanr(cmp['bulk_prolif'], cmp['p_cycling'])[0]:.3f}")
print("完成")
