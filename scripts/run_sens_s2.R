# S2：统一协变量敏感性（只 age+sex，去 stage）+ 协变量组合清单导出
# 目的：检验 stage 纳入与否对 pooled HR 与归因 R² 的影响（M6）
# 用法：cd 到项目根目录后 Rscript scripts/run_sens_s2.R
# 产出：results/sens/s2_age_sex.txt、results/sens/covariate_combo.csv
library(survival)
library(metafor)

dir.create("results/sens", showWarnings = FALSE)

# ============ 1. 读数据（复用 run_specB.R） ============
cat("[1/4] Reading expression matrix...\n")
expr_raw <- read.delim(gzfile("data/tcga_RSEM_gene_tpm.gz"), header = TRUE, check.names = FALSE, row.names = 1)
E <- as.matrix(expr_raw)
rownames(E) <- sub("\\..*$", "", rownames(E))
tumor_cols <- grepl("-01[A-Z]?$", colnames(E))
E <- E[, tumor_cols, drop = FALSE]

surv <- read.delim("data/survival_pancanatlas.tsv", header = TRUE, check.names = FALSE)
colnames(surv)[1] <- "sample"
surv <- surv[grepl("-01[A-Z]?$", surv$sample), ]
sample_to_cancer <- surv[["cancer type abbreviation"]][match(colnames(E), surv$sample)]
names(sample_to_cancer) <- colnames(E)

s2e <- read.delim("data/symbol2ensg.tsv", header = TRUE, stringsAsFactors = FALSE)

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
sig_names <- names(signatures)

# ============ 2. 算 B2 分数（癌种内标准化） ============
cat("[2/4] Computing B2 (within-cancer standardized) scores...\n")
pan_score <- function(ensg_ids) {
  g <- intersect(ensg_ids, rownames(E)); sub <- E[g, , drop = FALSE]
  z <- t(scale(t(sub))); colMeans(z, na.rm = TRUE)
}
S_pan <- do.call(rbind, lapply(sig_ensg, pan_score))
cancers <- sort(unique(sample_to_cancer[!is.na(sample_to_cancer)]))
S_B2 <- S_pan
for (ca in cancers) {
  idx <- which(sample_to_cancer == ca)
  if (length(idx) < 3) next
  for (nm in sig_names) {
    vals <- S_pan[nm, idx]; sd_c <- sd(vals, na.rm = TRUE)
    if (!is.na(sd_c) && sd_c > 1e-8) S_B2[nm, idx] <- (vals - mean(vals, na.rm = TRUE)) / sd_c
  }
}

# ============ 3. 逐癌种 Cox（只 age+sex） + 协变量组合清单 ============
cat("[3/4] Running per-cancer Cox (age+sex only)...\n")
merged <- merge(data.frame(sample = colnames(E)), surv, by = "sample", all = FALSE)
merged$age <- as.numeric(merged$age_at_initial_pathologic_diagnosis)
merged$OS.time <- as.numeric(merged$OS.time)
merged$OS <- as.numeric(merged$OS)
merged$stage <- as.character(merged$ajcc_pathologic_tumor_stage)
merged$stage_grp <- ifelse(grepl("^Stage III", merged$stage) | grepl("^Stage IV", merged$stage), "III-IV",
                    ifelse(grepl("^Stage II", merged$stage), "I-II",
                    ifelse(grepl("^Stage I", merged$stage), "I-II", NA)))

cancers_in <- unique(merged[["cancer type abbreviation"]])
cov_combo <- list()
res_list <- list()
for (ca in cancers_in) {
  d <- merged[merged[["cancer type abbreviation"]] == ca, ]
  n_total <- nrow(d)
  n_events <- sum(d$OS == 1, na.rm = TRUE)
  if (n_total < 30 || n_events < 5) next
  gender_ok <- length(unique(d$gender)) >= 2
  stage_ok <- length(unique(d$stage_grp[!is.na(d$stage_grp)])) >= 2
  cov_combo[[length(cov_combo)+1]] <- data.frame(
    cancer = ca, n = n_total, events = n_events,
    gender_ok = gender_ok, stage_ok = stage_ok, stringsAsFactors = FALSE)
  for (nm in sig_names) {
    scores <- S_B2[nm, match(d$sample, colnames(E))]
    d2 <- d[!is.na(scores) & !is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
    sc <- scores[!is.na(scores) & !is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0]
    if (length(sc) < 25) next
    d2$score <- as.numeric(sc)
    mod_terms <- "score + age"
    if (gender_ok) mod_terms <- paste(mod_terms, "+ gender")
    fit <- tryCatch(coxph(as.formula(paste("Surv(OS.time, OS) ~", mod_terms)), data = d2),
                    error = function(e) NULL)
    if (is.null(fit)) next
    s <- summary(fit)
    res_list[[length(res_list)+1]] <- data.frame(
      cancer = ca, signature = nm, n = nrow(d2),
      beta = s$coefficients[1, "coef"], se = s$coefficients[1, "se(coef)"],
      p = s$coefficients[1, "Pr(>|z|)"], stringsAsFactors = FALSE)
  }
}
res_s2 <- do.call(rbind, res_list)
cov_df <- do.call(rbind, cov_combo)
write.csv(cov_df, "results/sens/covariate_combo.csv", row.names = FALSE)

# ============ 4. pooled + 归因（age+sex 版）对比主分析 ============
cat("[4/4] Pooled + meta-regression (age+sex) vs main...\n")
mods <- read.csv("results/specB/moderators_abs.csv")
tmb <- read.csv("results/cancer_tmb.csv")
res_s2 <- merge(res_s2, mods, by = "cancer", all.x = TRUE)
res_s2 <- merge(res_s2, tmb, by = "cancer", all.x = TRUE)

# 主分析（B2, age+sex+stage）读取已有值
res_main <- read.csv("results/specB/per_cancer_three.csv")
res_main <- res_main[res_main$version == "B2", ]
res_main <- merge(res_main, mods, by = "cancer", all.x = TRUE)
res_main <- merge(res_main, tmb, by = "cancer", all.x = TRUE)

out <- c("===== S2 统一协变量（age+sex）vs 主分析（age+sex+stage）=====", "")
for (nm in sig_names) {
  d2 <- res_s2[res_s2$signature == nm & is.finite(res_s2$beta) & res_s2$se > 0, ]
  dm <- res_main[res_main$signature == nm & is.finite(res_main$beta) & res_main$se > 0, ]
  m2 <- rma(yi = beta, sei = se, data = d2, method = "REML")
  mm <- rma(yi = beta, sei = se, data = dm, method = "REML")
  mp2 <- rma(yi = beta, sei = se, mods = ~ prolif_abs, data = d2, method = "REML", test = "knha")
  mpm <- rma(yi = beta, sei = se, mods = ~ prolif_abs, data = dm, method = "REML", test = "knha")
  out <- c(out, sprintf("%-10s pooled: age+sex HR=%.3f(p=%.3f) vs main HR=%.3f(p=%.3f) | prolif R2: age+sex %.1f%% vs main %.1f%%",
                        nm, exp(m2$b), m2$pval, exp(mm$b), mm$pval, mp2$R2, mpm$R2))
}
writeLines(out, "results/sens/s2_age_sex.txt")
cat(out, sep = "\n")
cat("\nDONE. results/sens/s2_age_sex.txt + covariate_combo.csv\n")
