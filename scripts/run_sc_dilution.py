#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A5 单细胞稀释机制验证 - Step 2: 提取 per-癌种细胞比例
数据: data/atlas_celltype_annotated_481k.h5ad (481,738 cells, 30 cancer types, 19 cell types)
只读 anndata backed 模式, 不载入 X 矩阵。
只保留 sample_type == 'Tumour' 的细胞。
输出: results/sc_dilution/per_cancer_cell_fraction.csv
"""
import anndata as ad
import pandas as pd
import numpy as np
import os

H5AD = "data/atlas_celltype_annotated_481k.h5ad"
OUTDIR = "results/sc_dilution"
os.makedirs(OUTDIR, exist_ok=True)

print("读取 h5ad (backed)...")
adata = ad.read_h5ad(H5AD, backed="r")
obs = adata.obs

print(f"总细胞: {adata.n_obs}, 总癌种: {obs['cancer_type'].nunique()}")

# 只保留 Tumour
tumor = obs[obs["sample_type"] == "Tumour"].copy()
print(f"Tumour 细胞: {len(tumor)}")

# 细胞类型归类
IMMUNE = {
    "CD4+ T cell", "CD8+ T cell", "Macrophage", "B cell", "NK cell",
    "Regulatory T cell", "Plasma cell", "Mast cell", "Neutrophil/Granulocyte",
    "Monocyte",
}
STROMAL = {"Fibroblast/CAF", "Endothelial", "Pericyte/vSMC"}
NEURONAL = {"Neuronal (glioma)", "Neuronal (neuroblastoma)"}
MALIGNANT = {"Epithelial/Malignant"}
CYCLING = {"Cycling cell"}
NOISE = {"Low-quality/Doublet", "Erythroid"}

# 免疫亚型细分（用于与 7 个签名对应）
SUBTYPE = {
    "cd8": {"CD8+ T cell"},
    "cd4": {"CD4+ T cell"},
    "b": {"B cell", "Plasma cell"},
    "myeloid": {"Macrophage", "Monocyte", "Mast cell"},
    "neutrophil": {"Neutrophil/Granulocyte"},
    "nk": {"NK cell"},
    "treg": {"Regulatory T cell"},
}

rows = []
for ct, sub in tumor.groupby("cancer_type", observed=True):
    n_total = len(sub)
    def frac(cells):
        return sub["cell_type"].isin(cells).sum() / n_total
    row = {
        "cancer_type_sc": ct,
        "n_tumour": n_total,
        "p_malignant": frac(MALIGNANT),
        "p_cycling": frac(CYCLING),
        "p_immune": frac(IMMUNE),
        "p_stromal": frac(STROMAL),
        "p_neuronal": frac(NEURONAL),
        "p_noise": frac(NOISE),
    }
    for k, cells in SUBTYPE.items():
        row[f"p_{k}"] = frac(cells)
    rows.append(row)

df = pd.DataFrame(rows).sort_values("p_malignant", ascending=False)
out = os.path.join(OUTDIR, "per_cancer_cell_fraction.csv")
df.to_csv(out, index=False)
print(f"\n已写出 {out}, 癌种数 {len(df)}")

# 打印摘要 + Cycling 谱系线索
print("\n=== per-癌种比例摘要 (按恶性比例降序) ===")
print(df[["cancer_type_sc", "n_tumour", "p_malignant", "p_cycling", "p_immune", "p_stromal"]].round(4).to_string(index=False))

# Cycling 与恶性比例的关系（判断 Cycling 是否主要反映恶性增殖）
from scipy.stats import spearmanr
r, p = spearmanr(df["p_malignant"], df["p_cycling"])
print(f"\nCycling vs Malignant 比例 Spearman: r={r:.3f}, p={p:.4g}")
r2, p2 = spearmanr(df["p_immune"], df["p_cycling"])
print(f"Cycling vs Immune 比例 Spearman: r={r2:.3f}, p={p2:.4g}")

adata.file.close()
print("\n完成")
