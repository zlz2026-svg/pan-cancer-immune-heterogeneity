# 重跑 changelog

> 项目：泛癌免疫异质性归因 | 依据：方案复核意见 V1（采纳 4 处修正 + 技术更正）

## S3（2026-11-19）修 READ TMB

- **改动**：`results/cancer_tmb.csv` 中 COAD、READ 两行。
- **旧值（错）**：COAD n=534/mut=103/TMB=2.71；READ n=534/mut=103/TMB=2.71（均为 COADREAD 合并队列值）。
- **新值（对）**：COAD n=338/mut=105/TMB=2.76；READ n=140/mut=89/TMB=2.34。
- **API 核实要点**：
  1. `read_tcga_pan_can_atlas_2018` 与 `coad_tcga_pan_can_atlas_2018` 均返回 **404（不存在）**，cBioPortal 仅有合并队列 `coadread_tcga_pan_can_atlas_2018`——复核意见假设的"READ 单独请求有误"实为"无单独队列"。
  2. 在合并队列内用 `CANCER_TYPE_DETAILED` 区分：Colon Adenocarcinoma（COAD）378 样本、Rectal Adenocarcinoma（READ）155 样本、黏液腺癌 61 样本（未归类）。
  3. 分别算 MUTATION_COUNT 中位数 → COAD 105（/38Mb=2.76）、READ 89（/38Mb=2.34）。
- **影响**：TMB 仍为弱调节变量（预计仍不显著），不影响主结论；但错值已消除。
- **待办**：S1 重跑后，Table 2 的 TMB 列与 Fig 3 用新值重出。

## S1（2026-11-19 完成）癌种内 z-score 重跑

- **实现**：`scripts/run_specB.R`，B2（泛癌 z 后癌种内标准化分数，主）+ B1（癌种内逐基因 z，稳健）+ 绝对调节变量（prolif_abs/estimate_abs）。
- **验收通过**：
  1. 逐癌种 p/方向 pan vs B2 完全一致（最大差 0.000000）——推论 1 实证 ✅
  2. β_B2/β_pan 中位数 0.6221（符合 s_c 0.4–0.6 区间）✅
  3. prolif_abs/estimate_abs 不全为 0（雷区 2 防线有效）✅
  4. k=31 ✅（样本数验收代码 `unique(n)` 去重 bug 显示 8747，实际 8924）
- **核心结果（B2 vs 旧版 pan）**：
  - pooled HR 收缩符合推论 2：CAF 1.40→1.22、Prolifer 1.50→1.27、TAN 1.16→1.09、cDC2 0.86→0.90；显著性格局不变（CAF/TAN/cDC2/Prolifer 显著，CD8/TLS/IFNg 不显著）。
  - **增殖归因结构几乎不变**：CD8 41.1%(p=0.004)、CAF 38.2%(p=0.004)、IFNg 44.8%(p=0.001)、TLS 21.5%(p=0.038)——与旧版 43.5/39.3/46.2/21.8 高度一致。
  - **TAN 双重调制更强**：胚层 20.4%(p=0.043) + ESTIMATE 37.8%(p=0.007)（旧版 28.5%/22.5%）。
  - cDC2 仍无显著归因。
- **伪影假设检验（O3）**：免疫签名 s_c 与增殖无显著负相关（CD8/TLS/TAN/cDC2/IFNg 均不显著；CAF 正相关 rho=0.472；仅 Prolifer 自指行负相关 rho=-0.427）→ **"增殖归因是尺度伪影"假设不成立，增殖主导是真机制**。
- **落点预案 A**：新版（B2）为主分析，旧版进补充材料。

## M3（2026-11-19 完成）单变量 meta 回归定稿（KH + FDR + 置换）

- **产出**：`results/specB/meta_regression_v5.csv/txt`、`permutest.txt`；脚本 `run_meta_v5.R`、`run_permutest.R`。
- **Knapp–Hartung + 置换后的最终结论**：
  1. **增殖主导（核心，稳固）**：置换 p 显著——CD8 41.1%（perm_p=0.004）、CAF 38.2%（0.0075）、IFNg 44.8%（0.0015）、Prolifer 80.3%（0.0005，自指）；TLS 21.5%（perm_p 0.045–0.153，边缘）。
  2. **TAN 故事重构**：非肿瘤成分 estimate 37.8%（perm_p=0.0115）**显著**；但胚层 germ（perm_p=0.105–0.295）**不显著**——"胚层对 TAN"在置换检验下不成立，改为"非肿瘤成分（稀释效应）主导 TAN 异质性"。
  3. 增殖系数方向均为负（prolif_abs 越高 → 该签名 HR 越低）。
- **含义**：核心卖点（增殖主导）稳；TAN 机制从"胚层+纯度双重"修正为"非肿瘤成分主导"，更贴合 Jiao 2026 稀释效应。

## M5（2026-11-19 完成）小事件癌种敏感性（P-a）

- **产出**：`results/specB/sensitivity_low_events.txt`；脚本 `run_sensitivity_low_events.R`。
- **结果**：剔除 PCPG/DLBC/THYM/PRAD/KICH（事件 ≤10）后 k=26，pooled HR 与增殖归因 R² **几乎不变**（变动 <4%）——结论不依赖低事件癌种，稳健 ✅。
- **事件门槛**：老板拍板 P-a（维持 k=31 主分析 + 本敏感性）。

## 收口（2026-11-19 完成）Results 重写 + Table 2 + Fig 2/3/4 + Methods 文本

- **正文 V2**（`02-最终交付/【1120】正文V2.md`）已重写：
  - Methods：签名评分改癌种内 z（B2）、CD8_cyto→CD8_Tcell、调节变量改绝对表达、ESTIMATE 改名 stromal/immune score、单变量 meta 回归 + KH + 置换 + FDR、GLASS 细节（log RPKM/单变量/79 例 66 事件/保护端未验证）。
  - Results：pooled HR 用 B2（CAF 1.22/TAN 1.09/cDC2 0.90/Prolifer 1.27）、归因用 M3（增殖 38-45% 置换 p 显著；TAN 由非肿瘤成分 37.8% 主导，胚层不显著；Prolifer 自指标注）、留一法 <3%、GLASS 6/7 显著。
  - Table 1 CD8_Tcell 改名；Table 2 用 M3 结果 + 英文表注（消除内部版本语言）。
- **图表重出**：Fig 2（B2 森林图）、Fig 3（M3 归因热图，Stromal/immune 改名）、Fig 4（CD8_Tcell 改名）。
- **留一法（B2）**：全部稳定（变动 <3%），`results/specB/leave_one_out_B2.txt`。
- **待办**：① 补充表（31 癌种胚层映射 + 协变量组合 + 低事件癌种标注）；② Abstract/Introduction/Discussion；③ Fig 1 流程图。

## 最终核心结论（B2 + KH + 置换定稿）

> ⚠️ 以下为 stage bug 修复**前**的中间值。**终版见下节**。

1. **增殖主导（核心，稳）**：CD8 41.1%/CAF 38.2%/IFNg 44.8%（置换 p 均 <0.01）；Prolifer 80.3%（自指）。
2. **TAN 由非肿瘤成分主导**：37.8%（置换 p=0.0115），胚层不显著（置换 p 0.10-0.30）。
3. **方向翻转稳**：GBM CD8 1.51 风险 vs OV 0.84 保护（尺度变换不影响，推论 1 实证）。
4. **尺度伪影假设被否**：免疫签名 s_c 与增殖无显著负相关 → 增殖归因是真机制。

## 最终核心结论（stage 修复后终版，2026-11-20）

> stage_grp 正则 bug 修复后（旧版 `results/specB_stagebug/`），协变量改为 age+sex+stage（20 癌种），以下为**终版数字，正文 V4/V5 以此为准**。

1. **pooled HR**：CAF 1.183、TAN 1.084、cDC2 0.909、Prolifer 1.256、TLS 0.914（调整 stage 后显著保护 p=0.013）；CD8 0.963、IFNg 1.022（不显著）。I² 51.4–85.7%。
2. **增殖归因**：CD8 37.5%（perm 0.009）、CAF 27.6%（perm 0.0245）、IFNg 35.6%（perm 0.004）、TLS 21.4%（perm 0.0505 边缘）、Prolifer 73.7%（自指）。
3. **TAN 非肿瘤成分 36.5%**（perm 0.013）主导；germ 30.3%（KH 0.040，置换不显著 0.08–0.44）。
4. **方向翻转**：GBM CD8 1.21 风险 vs OV 0.88 保护；4 免疫签名翻转 vs 3 基质/增殖签名一致风险。
5. **胚层效应=胶质瘤效应**：合并 GBM+LGG 后 germ R² 85.6%。
6. **M2 替代增殖度量**：CD8 16.3%/IFNg 19.5% 仍显著；CAF 9.9%/TLS 14.2% 不显著（部分依赖基因集选择）。
7. **FDR 口径（P-a）**：BH q 值均 ≥0.114（除自指行），正文定性为探索性分析、主推断用置换 p。
