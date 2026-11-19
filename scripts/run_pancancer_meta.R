# 泛癌免疫标志物跨癌种效应异质性归因 — 核心分析脚本
# 用法：cd 到项目根目录后 Rscript scripts/run_pancancer_meta.R
library(survival)
library(metafor)

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

# ============ 1. 读表达矩阵 ============
cat("[1/6] Reading expression matrix...\n")
expr_raw <- read.delim(gzfile("data/tcga_RSEM_gene_tpm.gz"),
                       header = TRUE, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
# Toil 格式：行 = gene（第1列），列 = sample。若第1列名是 "sample" 则基因在列，需转置。
if (colnames(expr_raw)[1] == "sample" || rownames(expr_raw)[1] == "1") {
  expr_raw <- as.data.frame(t(expr_raw), stringsAsFactors = FALSE)
}
expr <- as.matrix(expr_raw)
cat("  Expression dim (gene x sample):", nrow(expr), "x", ncol(expr), "\n")

# ============ 2. 读生存数据 ============
cat("[2/6] Reading survival data...\n")
surv <- read.delim("data/survival_pancanatlas.tsv", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
colnames(surv)[1] <- "sample"
cat("  Survival rows:", nrow(surv), "\n")

# ============ 3. 免疫签名基因集（公认标志物） ============
signatures <- list(
  CD8_cyto  = c("CD8A","CD8B","CD3D","CD3E","GZMA","GZMB","PRF1","IFNG","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX"),
  TLS       = c("CXCL13","CCL19","CCL21","CD79A","MS4A1","LTB","TNFRSF17"),
  CAF       = c("COL1A1","COL1A2","COL3A1","ACTA2","FAP","PDGFRB","POSTN","FN1","MMP2"),
  TAN       = c("CSF3R","S100A8","S100A9","FCGR3B","CXCR2","OLR1"),
  cDC2      = c("CD1C","FCER1A","CLEC10A","CD1E"),
  IFNg_path = c("IFNG","CXCL9","CXCL10","CXCL11","GBP1","STAT1","IRF1","GBP2"),
  Prolifer  = c("MKI67","PCNA","TOP2A","BIRC5","CCNB1","MCM2","AURKA","CDK1")  # meta-PCNA 增殖
)

# ============ 4. 胚层映射（Hoadley 2018 简化版） ============
endoderm <- c("LIHC","COAD","READ","STAD","ESCA","PAAD","CHOL","LUAD","LUSC","THCA","CESC","BLCA","PRAD","THYM")
mesoderm <- c("KIRC","KIRP","KICH","OV","UCS","SARC","MESO","LAML","DLBC","TGCT","ACC","UCEC")
ectoderm <- c("GBM","LGG","SKCM","HNSC","PCPG","UVM","BRCA")
germ_map <- function(abbr) {
  ifelse(abbr %in% endoderm, "endoderm",
  ifelse(abbr %in% ectoderm, "ectoderm", "mesoderm"))
}

# ============ 5. 签名评分（矩阵已是 log2 尺度 -> 跨样本 z-score -> 均值） ============
cat("[3/6] Reading symbol->ENSG mapping & computing signature scores...\n")
s2e <- read.delim("data/symbol2ensg.tsv", header = TRUE, stringsAsFactors = FALSE)
E <- expr                                   # Toil RSEM 已是 log2(TPM+0.001)，不再二次 log
rownames(E) <- sub("\\..*$", "", rownames(E))  # 去掉 ENSG 版本号后缀
tumor_cols <- grepl("-01[A-Z]?$", colnames(E))  # 只保留原发肿瘤样本（排除 -11 正常等）
E <- E[, tumor_cols, drop = FALSE]
cat("  Tumor samples kept:", ncol(E), "\n")
sig_ensg <- lapply(signatures, function(syms) s2e$ensg[s2e$symbol %in% syms])
sig_score <- function(ensg_ids) {
  g <- intersect(ensg_ids, rownames(E))
  if (length(g) < 2) stop("too few genes matched: ", paste(ensg_ids, collapse=","))
  sub <- E[g, , drop = FALSE]
  z <- t(scale(t(sub)))                      # 每个基因跨样本 z-score
  colMeans(z, na.rm = TRUE)                  # 签名 = 基因均值
}
cat("  Signature scores will be computed per column below.\n")

# ============ 6. 匹配样本 + 逐癌种 Cox ============
cat("[4/6] Matching samples & running per-cancer Cox...\n")
# 统一 sample ID：去 -01/-11 后缀
# 直接用原始样本 ID（都是 -01 后缀格式）匹配，不去 trim
score_df <- data.frame(sample = colnames(E), stringsAsFactors = FALSE)
for (nm in names(sig_ensg)) {
  score_df[[nm]] <- sig_score(sig_ensg[[nm]])
}
surv <- surv[grepl("-01[A-Z]?$", surv$sample), ]  # 只保留 -01 原发肿瘤记录

merged <- merge(score_df, surv, by = "sample", all = FALSE)
merged$age <- as.numeric(merged$age_at_initial_pathologic_diagnosis)
merged$OS.time <- as.numeric(merged$OS.time)
merged$OS <- as.numeric(merged$OS)
merged$stage <- as.character(merged$ajcc_pathologic_tumor_stage)
merged$stage_grp <- ifelse(grepl("^Stage I|^Stage II", merged$stage), "I-II",
                    ifelse(grepl("^Stage III|^Stage IV", merged$stage), "III-IV", NA))
merged$germ <- germ_map(merged[["cancer type abbreviation"]])
cat("  Merged samples:", nrow(merged), "\n")
cat("  [DIAG] cancer types:", length(unique(merged[["cancer type abbreviation"]])), "\n")
cat("  [DIAG] OS events:", sum(merged$OS == 1, na.rm = TRUE), "\n")
cat("  [DIAG] OS NA:", sum(is.na(merged$OS)), "\n")
cat("  [DIAG] score NA (CD8_cyto):", sum(is.na(merged$CD8_cyto)), "\n")
cat("  [DIAG] has cancer col:", "cancer type abbreviation" %in% colnames(merged), "\n")

sig_names <- names(signatures)
res_list <- list()
for (ca in unique(merged[["cancer type abbreviation"]])) {
  d <- merged[merged[["cancer type abbreviation"]] == ca, ]
  cat("  [DIAG]", ca, "n=", nrow(d), "events=", sum(d$OS == 1, na.rm = TRUE), "\n")
  if (nrow(d) < 30 || sum(d$OS == 1, na.rm = TRUE) < 5) next
  for (sig in sig_names) {
    d2 <- d[!is.na(d[[sig]]) & !is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
    if (nrow(d2) < 25) next
    d2$score <- as.numeric(d2[[sig]])       # 固定变量名，避免公式内索引
    mod_terms <- "score + age"              # 动态纳入有 2+ 水平的协变量
    if (length(unique(d2$gender)) >= 2) mod_terms <- paste(mod_terms, "+ gender")
    if (length(unique(d2$stage_grp[!is.na(d2$stage_grp)])) >= 2) mod_terms <- paste(mod_terms, "+ stage_grp")
    fit <- tryCatch(coxph(as.formula(paste("Surv(OS.time, OS) ~", mod_terms)), data = d2),
                    error = function(e) { cat("  [ERR]", ca, sig, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(fit)) next
    s <- summary(fit)
    res_list[[length(res_list)+1]] <- data.frame(
      cancer = ca, signature = sig, germ = germ_map(ca),
      n = nrow(d2),
      HR = s$coefficients[1, "exp(coef)"],
      beta = s$coefficients[1, "coef"],
      se = s$coefficients[1, "se(coef)"],
      p = s$coefficients[1, "Pr(>|z|)"],
      purity_proxy = mean(d2$CAF, na.rm = TRUE),      # 基质/纯度代理（后续换 CPE）
      prolif = mean(d2$Prolifer, na.rm = TRUE)
    )
  }
}
res <- do.call(rbind, res_list)
write.csv(res, file.path(out_dir, "per_cancer_cox.csv"), row.names = FALSE)
cat("  Per-cancer Cox done, rows:", nrow(res), "\n")

# ============ 7. meta 分析 + meta 回归（异质性归因） ============
cat("[5/6] Running meta-analysis & meta-regression...\n")
meta_res <- list()
for (sig in sig_names) {
  d <- res[res$signature == sig & !is.na(res$se) & res$se > 0, ]
  if (nrow(d) < 8) next
  # 随机效应 meta
  m <- rma(yi = d$beta, sei = d$se, method = "REML")
  # meta 回归：胚层 + 纯度 + 增殖
  m2 <- tryCatch(rma(yi = d$beta, sei = d$se, mods = ~ germ + prolif, method = "REML"),
                 error = function(e) NULL)
  meta_res[[sig]] <- list(
    k = nrow(d),
    pooled_HR = exp(m$b[1]), ci_lb = exp(m$ci.lb), ci_ub = exp(m$ci.ub),
    I2 = m$I2, tau2 = m$tau2, Qp = m$QEp,
    meta_reg = if (is.null(m2)) NULL else list(
      R2 = m2$R2, p_germ = ifelse("germmesoderm" %in% rownames(m2$b), NA, NA),
      coef = m2$b, pval = m2$pval
    )
  )
}
cat("[6/6] Saving summary...\n")
sink(file.path(out_dir, "meta_summary.txt"))
for (sig in names(meta_res)) {
  r <- meta_res[[sig]]
  cat("\n=== ", sig, " (k=", r$k, ") ===\n", sep="")
  cat(sprintf("  pooled HR = %.3f (95%%CI %.3f-%.3f), I2=%.1f%%, tau2=%.4f, Q p=%.4f\n",
      r$pooled_HR, r$ci_lb, r$ci_ub, r$I2, r$tau2, r$Qp))
  if (!is.null(r$meta_reg)) {
    cat(sprintf("  meta-reg R2=%.3f\n", r$meta_reg$R2))
  }
}
sink()
cat("DONE. Results in", out_dir, "\n")
