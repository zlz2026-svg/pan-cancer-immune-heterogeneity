suppressPackageStartupMessages(library(survival))

expr_raw <- read.delim(gzfile("data/tcga_RSEM_gene_tpm.gz"), header = TRUE,
                       check.names = FALSE, row.names = 1)
E <- as.matrix(expr_raw)
rownames(E) <- sub("\\..*$", "", rownames(E))
E <- E[, grepl("-01[A-Z]?$", colnames(E)), drop = FALSE]

surv <- read.delim("data/survival_pancanatlas.tsv", header = TRUE, check.names = FALSE)
colnames(surv)[1] <- "sample"
surv <- surv[grepl("-01[A-Z]?$", surv$sample), ]

s2e <- read.delim("data/symbol2ensg.tsv", header = TRUE, stringsAsFactors = FALSE)
signatures <- list(
  CD8_cyto  = c("CD8A","CD8B","CD3D","CD3E","GZMA","GZMB","PRF1","IFNG","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX"),
  TLS       = c("CXCL13","CCL19","CCL21","CD79A","MS4A1","LTB","TNFRSF17"),
  CAF       = c("COL1A1","COL1A2","COL3A","ACTA2","FAP","PDGFRB","POSTN","FN1","MMP2"),
  TAN       = c("CSF3R","S100A8","S100A9","FCGR3B","CXCR2","OLR1"),
  cDC2      = c("CD1C","FCER1A","CLEC10A","CD1E"),
  IFNg_path = c("IFNG","CXCL9","CXCL10","CXCL11","GBP1","STAT1","IRF1","GBP2"),
  Prolifer  = c("MKI67","PCNA","TOP2A","BIRC5","CCNB1","MCM2","AURKA","CDK1"))
sig_ensg <- lapply(signatures, function(syms) s2e$ensg[s2e$symbol %in% syms])

pan_score <- function(ensg_ids) {
  g <- intersect(ensg_ids, rownames(E))
  colMeans(t(scale(t(E[g, , drop = FALSE]))), na.rm = TRUE)
}
S_pan <- do.call(rbind, lapply(sig_ensg, pan_score))
sample_to_cancer <- surv[["cancer type abbreviation"]][match(colnames(E), surv$sample)]
cancers <- sort(unique(sample_to_cancer[!is.na(sample_to_cancer)]))

S_B2 <- S_pan
for (ca in cancers) {
  idx <- which(sample_to_cancer == ca)
  for (nm in rownames(S_pan)) {
    vals <- S_pan[nm, idx]
    sdc <- sd(vals, na.rm = TRUE)
    if (!is.na(sdc) && sdc > 1e-8) S_B2[nm, idx] <- (vals - mean(vals, na.rm = TRUE)) / sdc
  }
}

merged <- merge(data.frame(sample = colnames(E)), surv, by = "sample", all = FALSE)
merged$age <- as.numeric(merged$age_at_initial_pathologic_diagnosis)
merged$OS.time <- as.numeric(merged$OS.time)
merged$OS <- as.numeric(merged$OS)
merged$stage <- as.character(merged$ajcc_pathologic_tumor_stage)
merged$stage_grp <- ifelse(grepl("^Stage III", merged$stage) | grepl("^Stage IV", merged$stage), "III-IV",
                    ifelse(grepl("^Stage II", merged$stage), "I-II",
                    ifelse(grepl("^Stage I", merged$stage), "I-II", NA)))

cov_viol_count <- c(age = 0, gender = 0, stage_grp = 0)
global_viol_models <- 0
cov_only_models <- 0
for (ca in unique(merged[["cancer type abbreviation"]])) {
  d <- merged[merged[["cancer type abbreviation"]] == ca, ]
  if (nrow(d) < 30 || sum(d$OS == 1, na.rm = TRUE) < 5) next
  for (nm in rownames(S_B2)) {
    scores <- S_B2[nm, match(d$sample, colnames(E))]
    d2 <- d[!is.na(scores) & !is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
    sc <- scores[!is.na(scores) & !is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0]
    if (length(sc) < 25) next
    d2$score <- as.numeric(sc)
    form <- "score + age"
    if (length(unique(d2$gender)) >= 2) form <- paste(form, "+ gender")
    if (length(unique(d2$stage_grp[!is.na(d2$stage_grp)])) >= 2) form <- paste(form, "+ stage_grp")
    fit <- tryCatch(coxph(as.formula(paste("Surv(OS.time, OS) ~", form)), data = d2),
                    error = function(e) NULL)
    if (is.null(fit)) next
    z <- tryCatch(cox.zph(fit), error = function(e) NULL)
    if (is.null(z)) next
    gp <- z$table["GLOBAL", "p"]
    if (gp >= 0.05) next
    global_viol_models <- global_viol_models + 1
    sp <- z$table["score", "p"]
    if (sp >= 0.05) cov_only_models <- cov_only_models + 1
    for (covname in names(cov_viol_count)) {
      if (covname %in% rownames(z$table) && z$table[covname, "p"] < 0.05) {
        cov_viol_count[covname] <- cov_viol_count[covname] + 1
      }
    }
  }
}
cat("global violations:", global_viol_models, "\n")
cat("covariate-only (score term ok):", cov_only_models, "\n")
cat("covariate rows violating across all global-violation models:\n")
print(cov_viol_count)
