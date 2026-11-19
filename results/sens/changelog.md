# 敏感性分析 changelog（2026-11-20）

> 项目：泛癌免疫异质性归因 | 触发：审稿意见 M2/M6/M7 补做 + 自检发现 stage 正则 bug

## 0. 重大发现：主分析 stage 协变量从未生效（正则 bug）

- **bug**：`run_specB.R` 的 stage_grp 正则 `grepl("^Stage I|^Stage II", ...)` 有前缀误匹配——"Stage I" 是 "Stage II/III/IV" 的共同前缀，导致所有 stage 都被归为 "I-II"，`stage_grp` 恒为单水平，从未被 Cox 纳入。
- **证据**：S2（age+sex）与主分析（宣称 age+sex+stage）结果逐位相同；covariate_combo.csv 中 31 癌种 stage_ok 全 FALSE。
- **修复**：改为按序匹配 `^Stage III|^Stage IV → III-IV`，再 `^Stage II → I-II`，再 `^Stage I → I-II`（验证 14 种 stage 取值全对）。
- **旧版归档**：`results/specB_stagebug/`（stage 未生效，等价 age+sex）。

## 1. 主分析重跑（stage 修复后，B2）

**pooled HR（age+sex+stage）vs 旧（=age+sex）：**

| 签名 | 旧 pooled HR (p) | 新 pooled HR (p) |
|---|---|---|
| CD8 | 0.975 (0.544) | 0.963 (0.349) |
| TLS | 0.933 (0.096) | **0.914 (0.013) 变显著保护** |
| CAF | 1.218 (<0.001) | 1.183 (<0.001) |
| TAN | 1.087 (0.010) | 1.084 (0.014) |
| cDC2 | 0.899 (0.003) | 0.909 (0.005) |
| IFNg | 1.032 (0.518) | 1.022 (0.639) |
| Prolifer | 1.265 (0.001) | 1.256 (<0.001) |

**归因（stage 修复后，KH + 置换 + FDR）：**

| 签名 | 增殖 R² (KH p) | TAN 相关 |
|---|---|---|
| CD8 | 37.5% (0.016, perm 0.009) | — |
| TLS | 21.4% (0.048, perm 0.0505) | — |
| CAF | 27.6% (0.020, perm 0.0245) | — |
| IFNg | 35.6% (0.010, perm 0.004) | — |
| Prolifer | 73.7% (<0.001, 自指) | — |
| TAN | — | estimate_abs 36.5% (perm 0.013)；germ 30.3% (KH p=0.040，perm 不显著 0.08–0.44) |

**核心结论不变**：增殖主导（CD8/CAF/IFNg 仍显著）、TAN 非肿瘤成分 36.5% 主导、方向翻转稳。TLS 调整 stage 后变显著保护（唯一实质变化）。

## 2. M2 替代增殖度量（Hallmark E2F/G2M 去重 Prolifer8 → 320 基因，r=0.957）

| 签名 | prolif_abs R² | prolif_alt R² (perm_p) | 判定 |
|---|---|---|---|
| CD8 | 37.5% | 16.3% (0.0465) | ✅ 仍显著（边缘） |
| TLS | 21.4% | 14.2% (0.0615) | ❌ 不显著 |
| CAF | 27.6% | 9.9% (0.0735) | ❌ 不显著 |
| IFNg | 35.6% | 19.5% (0.0240) | ✅ 仍显著 |
| Prolifer | 73.7% | 67.5% (0.0005) | 自指 |

**结论**：增殖归因对 **CD8、IFNg 是真机制**（替代度量仍显著）；对 **CAF、TLS 部分依赖基因集选择 + stage**（替代度量不显著）。Discussion 须如实写"增殖主导在 CD8/IFNg 稳健，CAF/TLS 较弱"。

## 3. S2 统一协变量（age+sex 敏感性）

- stage 调整影响小：pooled HR 方向全不变，CAF 增殖 R² 38.2%→27.6%（降 10%）、IFNg 44.8%→35.6%（降 9%）、CD8 41.1%→37.5%。
- 结论：stage 与增殖相关，调整后增殖解释力部分被 stage 吸收，但方向/显著性不变。

## 4. S4 胚层敏感性

- 剔 SARC：TAN germ 30.3%→21.6%（p=0.078），其余不变。
- **合并 GBM+LGG 为胶质瘤组（4 水平）：TAN germ R² 飙到 85.6%（p=0.002）、CD8 germ 55.0%（p=0.050）** → 印证"胚层效应主要是胶质瘤效应"（审稿 M7 的担忧成立）。
- 各胚层组合并效应：CAF 三组均显著风险（ectoderm 最强 1.280）；TAN ectoderm 1.258（p=0.022）显著、endoderm 1.076、mesoderm 0.983。

## 5. 最终定稿口径（写正文用）

1. **增殖主导（分签名）**：CD8 37.5%、IFNg 35.6% 稳健（替代度量下仍显著）；CAF 27.6%、TLS 21.4% 较弱（替代度量不显著，诚实标注）。
2. **TAN 非肿瘤成分 36.5%**（perm p=0.013）主导；germ 30.3% 置换不显著。
3. **方向翻转**：GBM CD8 1.51 vs OV 0.84。
4. **"胚层效应"= 胶质瘤效应**：合并 GBM+LGG 后 germ R² 85.6%。
5. **TLS 调整 stage 后显著保护**（0.914, p=0.013）。
6. 外部验证 GLASS 7/7 同向 6/7 显著（不受 stage 影响）。

## 6. 输出文件

- 主分析：`results/specB/`（已重跑，stage 修复版）
- 旧版：`results/specB_stagebug/`（stage 未生效）
- 敏感性：`results/sens/`（m2_prolif_alt.txt、s2_age_sex.txt、s4_germ_sensitivity.txt、covariate_combo.csv、prolif_alt.csv、m2_gene_selection_log.txt）
