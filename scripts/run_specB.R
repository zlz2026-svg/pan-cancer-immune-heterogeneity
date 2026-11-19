# S1：癌种内 z-score 重跑（B2 主 + B1 稳健）+ 绝对调节变量 + 验收
# 用法：cd 到项目根目录后 Rscript scripts/run_specB.R
# 产出：results/specB/ 下 per_cancer_three.csv、pooled_comparison.csv、within_cancer_sd.csv、sd_prolif_cor.csv、verification.txt
library(survival)
library(metafor)

dir.create("results/specB", showWarnings = FALSE)

# ============ 1. 读数据 ============
cat("[1/6] Reading expression matrix...\n")
expr_raw <- read.delim(gzfile("data/tcga_RSEM_gene_tpm.gz"), header=TRUE, check.names=FALSE, row.names=1)
E <- as.matrix(expr_raw)
rownames(E) <- sub("\\..*$", "", rownames(E))
tumor_cols <- grepl("-01[A-Z]?$", colnames(E))
E <- E[, tumor_cols, drop = FALSE]
cat("  Tumor samples:", ncol(E), "\n")

surv <- read.delim("data/survival_pancanatlas.tsv", header=TRUE, check.names=FALSE)
colnames(surv)[1] <- "sample"
surv <- surv[grepl("-01[A-Z]?$", surv$sample), ]

s2e <- read.delim("data/symbol2ensg.tsv", header=TRUE, stringsAsFactors=FALSE)
egs <- read.delim("data/estimate_genes_ensg.tsv", header=TRUE, stringsAsFactors=FALSE)
stromal_ensg <- egs$ensg[egs$set == "stromal"]
immune_ensg <- egs$ensg[egs$set == "immune"]

# 签名基因集（与主分析一致）
signatures <- list(
  CD8_cyto  = c("CD8A","CD8B","CD3D","CD3E","GZMA","GZMB","PRF1","IFNG","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX"),
  TLS       = c("CXCL13","CCL19","CCL21","CD79A","MS4A1","LTB","TNFRSF17"),
  CAF       = c("COL1A1","COL1A2","COL3A1","ACTA2","FAP","PDGFRB","POSTN","FN1","MMP2"),
  TAN       = c("CSF3R","S100A8","S100A9","FCGR3B","CXCR2","OLR1"),
  cDC2      = c("CD1C","FCER1A","CLEC10A","CD1E"),
  IFNg_path = c("IFNG","CXCL9","CXCL10","CXCL11","GBP1","STAT1","IRF1","GBP2"),
  Prolifer  = c("MKI67","PCNA","TOP2A","BIRC5","CCNB1","MCM2","AURKA","CDK1")
)
sig_ensg <- lapply(signatures, function(syms) s2e$ensg[s2e$symbol %in% syms])

# ============ 2. 算泛癌 z 版签名分（S_pan，旧口径） ============
cat("[2/6] Computing pan-cancer z scores...\n")
pan_score <- function(ensg_ids) {
  g <- intersect(ensg_ids, rownames(E))
  sub <- E[g, , drop = FALSE]
  z <- t(scale(t(sub)))
  colMeans(z, na.rm = TRUE)
}
S_pan <- do.call(rbind, lapply(sig_ensg, pan_score))   # 7 签名 × 样本

# 样本-癌种映射
sample_to_cancer <- surv[["cancer type abbreviation"]][match(colnames(E), surv$sample)]
names(sample_to_cancer) <- colnames(E)

# ============ 3. 算 B2（癌种内标准化 S_pan）与 s_c ============
cat("[3/6] Computing B2 (within-cancer standardized) and B1 (within-cancer per-gene z)...\n")
sig_names <- rownames(S_pan)
cancers <- sort(unique(sample_to_cancer[!is.na(sample_to_cancer)]))

# s_c 表：31 癌种 × 7 签名（S_pan 的癌种内 SD）
s_c <- matrix(NA, nrow = length(cancers), ncol = length(sig_names),
              dimnames = list(cancers, sig_names))
S_B2 <- S_pan   # 复制，癌种内标准化
for (ca in cancers) {
  idx <- which(sample_to_cancer == ca)
  if (length(idx) < 3) next
  for (nm in sig_names) {
    vals <- S_pan[nm, idx]
    sd_c <- sd(vals, na.rm = TRUE)
    s_c[ca, nm] <- sd_c
    if (!is.na(sd_c) && sd_c > 1e-8) {
      S_B2[nm, idx] <- (vals - mean(vals, na.rm = TRUE)) / sd_c
    }
  }
}

# B1：癌种内逐基因 z + 均值
S_B1 <- matrix(NA, nrow = length(sig_names), ncol = ncol(E), dimnames = dimnames(S_pan))
for (ca in cancers) {
  idx <- which(sample_to_cancer == ca)
  if (length(idx) < 3) next
  for (nm in sig_names) {
    g <- intersect(sig_ensg[[nm]], rownames(E))
    sub <- E[g, idx, drop = FALSE]
    # 癌种内逐基因 z；基因方差过小(除零)则剔除该基因
    z <- apply(sub, 1, function(x) { s <- sd(x, na.rm=TRUE); if (is.na(s)||s<1e-8) rep(NA,length(x)) else (x-mean(x,na.rm=TRUE))/s })
    S_B1[nm, idx] <- rowMeans(z, na.rm = TRUE)
  }
}

# ============ 4. 绝对调节变量（雷区 2 防线） ============
cat("[4/6] Computing absolute-level moderators...\n")
prolif_abs <- sapply(cancers, function(ca) {
  idx <- which(sample_to_cancer == ca)
  g <- intersect(sig_ensg[["Prolifer"]], rownames(E))
  mean(E[g, idx], na.rm = TRUE)   # log2(TPM+0.001) 绝对均值
})
est_abs <- sapply(cancers, function(ca) {
  idx <- which(sample_to_cancer == ca)
  sg <- intersect(stromal_ensg, rownames(E)); ig <- intersect(immune_ensg, rownames(E))
  mean(E[sg, idx], na.rm = TRUE) + mean(E[ig, idx], na.rm = TRUE)
})
mod_abs <- data.frame(cancer = cancers, prolif_abs = prolif_abs, estimate_abs = est_abs, stringsAsFactors = FALSE)
write.csv(mod_abs, "results/specB/moderators_abs.csv", row.names = FALSE)
cat("  prolif_abs 全为0?", all(prolif_abs == 0), " | estimate_abs 全为0?", all(est_abs == 0), "\n")

# ============ 5. 逐癌种 Cox（三口径：S_pan / B2 / B1） ============
cat("[5/6] Running per-cancer Cox (3 score versions)...\n")
merged <- merge(data.frame(sample = colnames(E)), surv, by = "sample", all = FALSE)
merged$age <- as.numeric(merged$age_at_initial_pathologic_diagnosis)
merged$OS.time <- as.numeric(merged$OS.time)
merged$OS <- as.numeric(merged$OS)
merged$stage <- as.character(merged$ajcc_pathologic_tumor_stage)
merged$stage_grp <- ifelse(grepl("^Stage III", merged$stage) | grepl("^Stage IV", merged$stage), "III-IV",
                    ifelse(grepl("^Stage II", merged$stage), "I-II",
                    ifelse(grepl("^Stage I", merged$stage), "I-II", NA)))

run_cox <- function(score_mat, tag) {
  res_list <- list()
  for (ca in unique(merged[["cancer type abbreviation"]])) {
    d <- merged[merged[["cancer type abbreviation"]] == ca, ]
    if (nrow(d) < 30 || sum(d$OS == 1, na.rm = TRUE) < 5) next
    for (nm in sig_names) {
      scores <- score_mat[nm, match(d$sample, colnames(E))]
      d2 <- d[!is.na(scores) & !is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
      sc <- scores[!is.na(scores) & !is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0]
      if (length(sc) < 25) next
      d2$score <- as.numeric(sc)
      mod_terms <- "score + age"
      if (length(unique(d2$gender)) >= 2) mod_terms <- paste(mod_terms, "+ gender")
      if (length(unique(d2$stage_grp[!is.na(d2$stage_grp)])) >= 2) mod_terms <- paste(mod_terms, "+ stage_grp")
      fit <- tryCatch(coxph(as.formula(paste("Surv(OS.time, OS) ~", mod_terms)), data = d2),
                      error = function(e) NULL)
      if (is.null(fit)) next
      s <- summary(fit)
      res_list[[length(res_list)+1]] <- data.frame(
        cancer = ca, signature = nm, version = tag, n = nrow(d2),
        beta = s$coefficients[1, "coef"], se = s$coefficients[1, "se(coef)"],
        p = s$coefficients[1, "Pr(>|z|)"], stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, res_list)
}
res_pan <- run_cox(S_pan, "pan")
res_B2  <- run_cox(S_B2, "B2")
res_B1  <- run_cox(S_B1, "B1")
res3 <- rbind(res_pan, res_B2, res_B1)
write.csv(res3, "results/specB/per_cancer_three.csv", row.names = FALSE)
cat("  Per-cancer Cox rows:", nrow(res3), "\n")

# ============ 6. pooled 对比 + O1/O2/O3 + 验收 ============
cat("[6/6] Pooled comparison & verification...\n")
pooled <- function(d) {
  m <- rma(yi = beta, sei = se, data = d, method = "REML")
  c(HR = exp(as.numeric(m$b)), lo = exp(as.numeric(m$ci.lb)), hi = exp(as.numeric(m$ci.ub)),
    I2 = m$I2, tau2 = m$tau2, Qp = m$QEp)
}
o1 <- list()
for (nm in sig_names) {
  dp <- res_pan[res_pan$signature == nm & is.finite(res_pan$beta) & res_pan$se > 0, ]
  d2 <- res_B2[res_B2$signature == nm & is.finite(res_B2$beta) & res_B2$se > 0, ]
  d1 <- res_B1[res_B1$signature == nm & is.finite(res_B1$beta) & res_B1$se > 0, ]
  o1[[nm]] <- data.frame(signature = nm,
    pan_HR = pooled(dp)["HR"], B2_HR = pooled(d2)["HR"], B1_HR = pooled(d1)["HR"],
    pan_I2 = pooled(dp)["I2"], B2_I2 = pooled(d2)["I2"], B1_I2 = pooled(d1)["I2"])
}
O1 <- do.call(rbind, o1)
write.csv(O1, "results/specB/pooled_comparison.csv", row.names = FALSE)

# O2：s_c 表
write.csv(as.data.frame(s_c), "results/specB/within_cancer_sd.csv")

# O3：s_c 与增殖的相关
sd_prolif <- data.frame(signature = sig_names, pearson_r = NA, pearson_p = NA,
                        spearman_rho = NA, spearman_p = NA)
for (nm in sig_names) {
  x <- s_c[, nm]; y <- prolif_abs[rownames(s_c)]
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) >= 5) {
    ct <- cor.test(x[ok], y[ok], method = "pearson")
    cs <- cor.test(x[ok], y[ok], method = "spearman")
    sd_prolif[sd_prolif$signature == nm, ] <- c(nm, ct$estimate, ct$p.value, cs$estimate, cs$p.value)
  }
}
write.csv(sd_prolif, "results/specB/sd_prolif_cor.csv", row.names = FALSE)

# 验收
v <- c()
v <- c(v, sprintf("【验收1】逐癌种 p/方向一致：pan vs B2 的 p 最大差 = %.6f", 
                  max(abs(res_pan$p - res_B2$p), na.rm = TRUE)))
# 验收2：β_B2/β_pan ≈ s_c（抽查整体比值）
ratio <- res_B2$beta / res_pan$beta
v <- c(v, sprintf("【验收2】β_B2/β_pan 中位数 = %.4f（应接近各癌种 s_c）", median(ratio, na.rm = TRUE)))
v <- c(v, sprintf("【验收3】prolif_abs 全为0 = %s；estimate_abs 全为0 = %s", all(prolif_abs==0), all(est_abs==0)))
v <- c(v, sprintf("【验收4】k = %d（应31）；样本数 = %d（应8924）",
                  length(unique(res_B2$cancer)), sum(res_B2$n[!duplicated(res_B2$cancer)])))
writeLines(v, "results/specB/verification.txt")
cat(v, sep = "\n")
cat("\nDONE. Outputs in results/specB/\n")
