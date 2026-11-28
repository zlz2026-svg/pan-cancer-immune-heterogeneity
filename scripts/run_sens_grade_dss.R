# Sensitivity analyses S1 (grade adjustment) + S2 (DSS endpoint)
# Usage: cd to project root; Rscript scripts/run_sens_grade_dss.R
# Outputs: results/sens_grade_dss/per_cancer_sens.csv, sens_meta_compare.txt
# Four versions: B2_repro (OS, no grade; reproduces main B2), B2_grade (OS + grade),
#                B2_dss (DSS, no grade), B2_dss_grade (DSS + grade)
library(survival)
library(metafor)

dir.create("results/sens_grade_dss", showWarnings = FALSE)

# ============ 1. Read data ============
cat("[1/4] Reading expression matrix...\n")
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

# ============ 2. Signature gene sets (same as run_specB.R) ============
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

# ============ 3. Pan-cancer z score + B2 within-cancer standardization ============
pan_score <- function(ensg_ids) {
  g <- intersect(ensg_ids, rownames(E))
  sub <- E[g, , drop = FALSE]
  z <- t(scale(t(sub)))
  colMeans(z, na.rm = TRUE)
}
S_pan <- do.call(rbind, lapply(sig_ensg, pan_score))
sig_names <- rownames(S_pan)
sample_to_cancer <- surv[["cancer type abbreviation"]][match(colnames(E), surv$sample)]
names(sample_to_cancer) <- colnames(E)
cancers <- sort(unique(sample_to_cancer[!is.na(sample_to_cancer)]))

S_B2 <- S_pan
for (ca in cancers) {
  idx <- which(sample_to_cancer == ca)
  if (length(idx) < 3) next
  for (nm in sig_names) {
    vals <- S_pan[nm, idx]
    sd_c <- sd(vals, na.rm = TRUE)
    if (!is.na(sd_c) && sd_c > 1e-8) {
      S_B2[nm, idx] <- (vals - mean(vals, na.rm = TRUE)) / sd_c
    }
  }
}

# ============ 4. Clinical covariates + grade + endpoints ============
merged <- merge(data.frame(sample = colnames(E)), surv, by = "sample", all = FALSE)
merged$age <- as.numeric(merged$age_at_initial_pathologic_diagnosis)
merged$OS.time <- as.numeric(merged$OS.time)
merged$OS <- as.numeric(merged$OS)
merged$DSS.time <- as.numeric(merged$DSS.time)
merged$DSS <- as.numeric(merged$DSS)
merged$stage <- as.character(merged$ajcc_pathologic_tumor_stage)
merged$stage_grp <- ifelse(grepl("^Stage III", merged$stage) | grepl("^Stage IV", merged$stage), "III-IV",
                    ifelse(grepl("^Stage II", merged$stage), "I-II",
                    ifelse(grepl("^Stage I", merged$stage), "I-II", NA)))

# grade: low = G1/G2/Low Grade, high = G3/G4/High Grade; GX/Unknown/blank = NA
g <- as.character(merged$histological_grade)
merged$grade_grp <- ifelse(g %in% c("G1","G2","Low Grade"), "low",
                    ifelse(g %in% c("G3","G4","High Grade"), "high", NA))

# ============ 5. Per-cancer Cox (4 versions) ============
run_cox <- function(score_mat, tag, endpoint, use_grade) {
  tcol <- paste0(endpoint, ".time")
  ecol <- endpoint
  res_list <- list()
  for (ca in unique(merged[["cancer type abbreviation"]])) {
    d <- merged[merged[["cancer type abbreviation"]] == ca, ]
    if (nrow(d) < 30 || sum(d[[ecol]] == 1, na.rm = TRUE) < 5) next
    for (nm in sig_names) {
      scores <- score_mat[nm, match(d$sample, colnames(E))]
      keep <- !is.na(scores) & !is.na(d[[ecol]]) & !is.na(d[[tcol]]) & d[[tcol]] > 0
      d2 <- d[keep, ]
      sc <- scores[keep]
      if (length(sc) < 25) next
      d2$score <- as.numeric(sc)
      mod_terms <- "score + age"
      if (length(unique(d2$gender)) >= 2) mod_terms <- paste(mod_terms, "+ gender")
      if (length(unique(d2$stage_grp[!is.na(d2$stage_grp)])) >= 2) mod_terms <- paste(mod_terms, "+ stage_grp")
      if (use_grade) {
        n_low  <- sum(d2$grade_grp == "low",  na.rm = TRUE)
        n_high <- sum(d2$grade_grp == "high", na.rm = TRUE)
        if (n_low >= 15 && n_high >= 15) mod_terms <- paste(mod_terms, "+ grade_grp")
      }
      fit <- tryCatch(coxph(as.formula(paste("Surv(", tcol, ",", ecol, ") ~", mod_terms)), data = d2),
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

cat("[2/4] B2 repro (OS, no grade)...\n")
res_repro <- run_cox(S_B2, "B2_repro", "OS", FALSE)
cat("[3/4] B2+grade (OS), B2 (DSS), B2+grade (DSS)...\n")
res_grade <- run_cox(S_B2, "B2_grade", "OS", TRUE)
res_dss   <- run_cox(S_B2, "B2_dss",   "DSS", FALSE)
res_dssg  <- run_cox(S_B2, "B2_dss_grade", "DSS", TRUE)
res_all <- rbind(res_repro, res_grade, res_dss, res_dssg)
write.csv(res_all, "results/sens_grade_dss/per_cancer_sens.csv", row.names = FALSE)
cat("  total rows:", nrow(res_all), "\n")

# ============ 6. Meta-regression comparison ============
mods <- read.csv("results/specB/moderators_abs.csv")
tmb <- read.csv("results/cancer_tmb.csv")
res <- merge(res_all, mods, by = "cancer", all.x = TRUE)
res <- merge(res, tmb, by = "cancer", all.x = TRUE)

endoderm <- c("LIHC","COAD","READ","STAD","ESCA","PAAD","CHOL","LUAD","LUSC","THCA","CESC","BLCA","PRAD","THYM")
mesoderm <- c("KIRC","KIRP","KICH","OV","UCS","SARC","MESO","LAML","DLBC","TGCT","ACC","UCEC")
ectoderm <- c("GBM","LGG","SKCM","HNSC","PCPG","UVM","BRCA")
germ <- function(a) ifelse(a %in% endoderm, "endoderm", ifelse(a %in% ectoderm, "ectoderm", "mesoderm"))
res$germ <- factor(germ(res$cancer), levels = c("endoderm","mesoderm","ectoderm"))

cat("[4/4] Meta-regression comparison...\n")
out <- c()
versions <- c("B2_repro","B2_grade","B2_dss","B2_dss_grade")
for (ver in versions) {
  dv <- res[res$version == ver & is.finite(res$beta) & res$se > 0, ]
  out <- c(out, sprintf("===== version %s (k=%d cancer types) =====", ver, length(unique(dv$cancer))))
  for (nm in sig_names) {
    d <- dv[dv$signature == nm, ]
    if (nrow(d) < 8) next
    pfit <- tryCatch(rma(yi = beta, sei = se, mods = ~ prolif_abs,   data = d, method = "REML", test = "knha"), error = function(e) NULL)
    efit <- tryCatch(rma(yi = beta, sei = se, mods = ~ estimate_abs, data = d, method = "REML", test = "knha"), error = function(e) NULL)
    pR2 <- if (is.null(pfit)) NA else pfit$R2
    pP  <- if (is.null(pfit)) NA else pfit$pval[2]
    eR2 <- if (is.null(efit)) NA else efit$R2
    eP  <- if (is.null(efit)) NA else efit$pval[2]
    out <- c(out, sprintf("  %-10s prolif R2=%5.1f%% p=%.4f | estimate R2=%5.1f%% p=%.4f",
        nm, pR2 * 100, pP, eR2 * 100, eP))
  }
}
writeLines(out, "results/sens_grade_dss/sens_meta_compare.txt")
cat(out, sep = "\n")
cat("\nDONE\n")
