#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A5 稀释检验链 - Step 3: 癌种映射 + 三重 Spearman 检验
检验 1: 单细胞 Cycling/恶性比例 vs bulk Prolifer 签名 (prolif_abs) —— 应正相关
检验 2: 单细胞恶性比例 vs 免疫比例 —— 应负相关 (稀释)
检验 3: bulk estimate_abs (stromal+immune) vs 单细胞 (免疫+基质) 比例 —— 应正相关
输出: results/sc_dilution/sc_bulk_comparison.csv + sc_dilution_tests.txt
"""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import os

OUTDIR = "results/sc_dilution"
sc = pd.read_csv(os.path.join(OUTDIR, "per_cancer_cell_fraction.csv"))
mods = pd.read_csv("results/specB/moderators_abs.csv")
per = pd.read_csv("results/specB/per_cancer_three.csv")

# bulk 只取 B2 版的 per-cancer 数据（用于后续如需要），moderators 用 prolif_abs/estimate_abs
mods.columns = [c.strip('"') for c in mods.columns]
mods["cancer"] = mods["cancer"].str.strip('"')

# 单细胞癌种 -> TCGA 癌种映射
# 一对一
one2one = {
    "Ovarian Carcinoma": "OV",
    "Head and Neck Squamous Cell Carcinoma": "HNSC",
    "Lung Adenocarcinoma": "LUAD",
    "Lung Squamous Cell Carcinoma": "LUSC",
    "Pancreatic Ductal Adenocarcinoma": "PAAD",
    "Bladder Cancer": "BLCA",
    "Breast Cancer": "BRCA",
    "Cervical Cancer": "CESC",
    "Clear Cell Renal Cell Carcinoma": "KIRC",
    "Papillary Renal Cell Carcinoma": "KIRP",
    "Hepatocellular Cancer": "LIHC",
    "Prostate Cancer": "PRAD",
    "Distal Cholangiocarcinoma": "CHOL",
    "Papillary Thyroid Carcinoma": "THCA",
}

# 合并映射: 单细胞癌种 -> (TCGA 癌种列表, 标签)
merged = {
    "Glioma": (["GBM", "LGG"], "Glioma"),
    "Colorectal Cancer": (["COAD", "READ"], "Colorectal"),
}

# Esophageal: 两个单细胞癌种合并 -> ESCA
eso_sc = ["Esophageal Adenocarcinoma", "Esophageal Squamous Cell Carcinoma"]

# 组装映射后的数据
rows = []
for sc_ct, tcga in one2one.items():
    r = sc[sc["cancer_type_sc"] == sc_ct]
    if len(r) == 0:
        continue
    m = mods[mods["cancer"] == tcga]
    if len(m) == 0:
        continue
    rows.append({
        "label": tcga, "sc_cancer": sc_ct,
        "n_tumour": r["n_tumour"].iloc[0],
        "p_malignant": r["p_malignant"].iloc[0],
        "p_cycling": r["p_cycling"].iloc[0],
        "p_immune": r["p_immune"].iloc[0],
        "p_stromal": r["p_stromal"].iloc[0],
        "p_cd8": r["p_cd8"].iloc[0],
        "p_myeloid": r["p_myeloid"].iloc[0],
        "bulk_prolif": m["prolif_abs"].iloc[0],
        "bulk_estimate": m["estimate_abs"].iloc[0],
    })

# 合并癌种 (Glioma, Colorectal): 单细胞已是合并, bulk 取 TCGA 两癌种均值
for sc_ct, (tcgas, label) in merged.items():
    r = sc[sc["cancer_type_sc"] == sc_ct]
    if len(r) == 0:
        continue
    m = mods[mods["cancer"].isin(tcgas)]
    if len(m) != len(tcgas):
        continue
    rows.append({
        "label": label, "sc_cancer": sc_ct,
        "n_tumour": r["n_tumour"].iloc[0],
        "p_malignant": r["p_malignant"].iloc[0],
        "p_cycling": r["p_cycling"].iloc[0],
        "p_immune": r["p_immune"].iloc[0],
        "p_stromal": r["p_stromal"].iloc[0],
        "p_cd8": r["p_cd8"].iloc[0],
        "p_myeloid": r["p_myeloid"].iloc[0],
        "bulk_prolif": m["prolif_abs"].mean(),
        "bulk_estimate": m["estimate_abs"].mean(),
    })

# Esophageal: 两个单细胞癌种按细胞数加权合并 -> ESCA
eso = sc[sc["cancer_type_sc"].isin(eso_sc)]
if len(eso) == 2:
    n = eso["n_tumour"].sum()
    m = mods[mods["cancer"] == "ESCA"]
    if len(m) == 1:
        rows.append({
            "label": "ESCA", "sc_cancer": "Esophageal (merged)",
            "n_tumour": n,
            "p_malignant": (eso["p_malignant"] * eso["n_tumour"]).sum() / n,
            "p_cycling": (eso["p_cycling"] * eso["n_tumour"]).sum() / n,
            "p_immune": (eso["p_immune"] * eso["n_tumour"]).sum() / n,
            "p_stromal": (eso["p_stromal"] * eso["n_tumour"]).sum() / n,
            "p_cd8": (eso["p_cd8"] * eso["n_tumour"]).sum() / n,
            "p_myeloid": (eso["p_myeloid"] * eso["n_tumour"]).sum() / n,
            "bulk_prolif": m["prolif_abs"].iloc[0],
            "bulk_estimate": m["estimate_abs"].iloc[0],
        })

cmp = pd.DataFrame(rows)
cmp["p_immune_stromal"] = cmp["p_immune"] + cmp["p_stromal"]
cmp.to_csv(os.path.join(OUTDIR, "sc_bulk_comparison.csv"), index=False)
print(f"映射后癌种数: {len(cmp)}")

# 检验 2 (纯单细胞, 全 30 癌种 + 映射子集)
lines = []
lines.append("=" * 60)
lines.append("A5 稀释检验链结果")
lines.append("=" * 60)

# 检验 1: Cycling vs bulk Prolifer
r, p = spearmanr(cmp["p_cycling"], cmp["bulk_prolif"])
lines.append(f"\n[检验1a] 单细胞 Cycling 比例 vs bulk Prolifer 签名: Spearman r={r:.3f}, p={p:.4g}  (k={len(cmp)})")

r, p = spearmanr(cmp["p_malignant"], cmp["bulk_prolif"])
lines.append(f"[检验1b] 单细胞 恶性比例 vs bulk Prolifer 签名: Spearman r={r:.3f}, p={p:.4g}")

# 检验 2: 恶性 vs 免疫 (纯单细胞)
r_all, p_all = spearmanr(sc["p_malignant"], sc["p_immune"])
lines.append(f"\n[检验2a] 恶性比例 vs 免疫比例 (全 30 单细胞癌种): Spearman r={r_all:.3f}, p={p_all:.4g}")
r_sub, p_sub = spearmanr(cmp["p_malignant"], cmp["p_immune"])
lines.append(f"[检验2b] 恶性比例 vs 免疫比例 (映射 17 癌种): Spearman r={r_sub:.3f}, p={p_sub:.4g}")

# 检验 3: bulk estimate vs 单细胞免疫+基质
r, p = spearmanr(cmp["bulk_estimate"], cmp["p_immune_stromal"])
lines.append(f"\n[检验3a] bulk estimate_abs vs 单细胞(免疫+基质)比例: Spearman r={r:.3f}, p={p:.4g}")
r, p = spearmanr(cmp["bulk_estimate"], cmp["p_immune"])
lines.append(f"[检验3b] bulk estimate_abs vs 单细胞免疫比例: Spearman r={r:.3f}, p={p:.4g}")

# 附加: 恶性比例 vs 免疫比例 的 Pearson (敏感性)
from scipy.stats import pearsonr
r_pe, p_pe = pearsonr(cmp["p_malignant"], cmp["p_immune"])
lines.append(f"\n[附加] 恶性比例 vs 免疫比例 Pearson (映射 17): r={r_pe:.3f}, p={p_pe:.4g}")

lines.append("\n" + "=" * 60)
lines.append("映射后数据 (per 癌种)")
lines.append(cmp[["label", "n_tumour", "p_malignant", "p_cycling", "p_immune", "p_cd8", "bulk_prolif", "bulk_estimate"]].round(3).to_string(index=False))

out_txt = os.path.join(OUTDIR, "sc_dilution_tests.txt")
with open(out_txt, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print("\n".join(lines))
print(f"\n已写出 {out_txt}")
