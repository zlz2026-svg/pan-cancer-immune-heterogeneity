# A6: OV-AU ovarian cancer external validation (CD8 protective direction)
# Usage: cd project root; Rscript scripts/run_ovau_validation.R
# Outputs: results/ovau_validation_cox.csv, results/ovau_validation_summary.txt
library(survival)

cat("[1] Reading expression matrix...\n")
GX <- readRDS("data/ovau/ICGC_OV-AU_GX.rds")
cat("  GX:", nrow(GX), "genes x", ncol(GX), "samples | range:", paste(range(GX), collapse=" .. "), "\n")

sample_map <- read.delim("data/ovau/sample_OVAU.tsv", header = TRUE, stringsAsFactors = FALSE)
donor <- read.delim("data/ovau/donor_clean.tsv", header = TRUE, stringsAsFactors = FALSE)

# map sample -> donor
smp <- sample_map[match(colnames(GX), sample_map$icgc_sample_id), ]
map <- data.frame(sample = colnames(GX), donor = smp$icgc_donor_id, stringsAsFactors = FALSE)
cat("  mapped samples:", sum(!is.na(map$donor)), "/", ncol(GX),
    "| unique donors:", length(unique(map$donor[!is.na(map$donor)])), "\n")

# survival: donor-level (deceased = event, survival_time in days)
donor$event <- ifelse(donor$donor_vital_status == "deceased", 1, 0)
donor$time <- as.numeric(donor$donor_survival_time)
donor$age <- as.numeric(donor$donor_age_at_diagnosis)

map$event <- donor$event[match(map$donor, donor$icgc_donor_id)]
map$time <- donor$time[match(map$donor, donor$icgc_donor_id)]
map$age <- donor$age[match(map$donor, donor$icgc_donor_id)]

# keep samples with valid survival; one sample per donor (dedup)
map <- map[!is.na(map$event) & !is.na(map$time) & map$time > 0, ]
map <- map[!duplicated(map$donor), ]
cat("  final samples (one per donor):", nrow(map), "| events:", sum(map$event), "\n")

# signatures (names match main analysis)
signatures <- list(
  CD8_cyto  = c("CD8A","CD8B","CD3D","CD3E","GZMA","GZMB","PRF1","IFNG","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX"),
  TLS       = c("CXCL13","CCL19","CCL21","CD79A","MS4A1","LTB","TNFRSF17"),
  CAF       = c("COL1A1","COL1A2","COL3A1","ACTA2","FAP","PDGFRB","POSTN","FN1","MMP2"),
  TAN       = c("CSF3R","S100A8","S100A9","FCGR3B","CXCR2","OLR1"),
  cDC2      = c("CD1C","FCER1A","CLEC10A","CD1E"),
  IFNg_path = c("IFNG","CXCL9","CXCL10","CXCL11","GBP1","STAT1","IRF1","GBP2"),
  Prolifer  = c("MKI67","PCNA","TOP2A","BIRC5","CCNB1","MCM2","AURKA","CDK1")
)

# score: within-cohort per-gene z + mean (same as GLASS validation)
score_sig <- function(genes) {
  g <- genes[genes %in% rownames(GX)]
  sub <- GX[g, map$sample, drop = FALSE]
  z <- t(scale(t(sub)))
  colMeans(z, na.rm = TRUE)
}

# per-signature univariable Cox
res <- data.frame(signature = character(), n = numeric(), events = numeric(),
                  HR = numeric(), lo = numeric(), hi = numeric(), p = numeric(),
                  stringsAsFactors = FALSE)
for (nm in names(signatures)) {
  sc <- score_sig(signatures[[nm]])
  d <- data.frame(time = map$time, event = map$event, score = sc)
  d <- d[!is.na(d$score) & is.finite(d$score), ]
  if (nrow(d) < 30 || sum(d$event) < 5) next
  fit <- tryCatch(coxph(Surv(time, event) ~ score, data = d), error = function(e) NULL)
  if (is.null(fit)) next
  s <- summary(fit)
  ci <- exp(confint(fit))
  res <- rbind(res, data.frame(signature = nm, n = nrow(d), events = sum(d$event),
    HR = exp(coef(fit)), lo = ci[1], hi = ci[2], p = s$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE))
}
write.csv(res, "results/ovau_validation_cox.csv", row.names = FALSE)

# TCGA OV reference (B2)
main <- read.csv("results/specB/per_cancer_three.csv")
ov <- main[main$cancer == "OV" & main$version == "B2", c("signature","beta","p")]
ov$HR_TCGA_OV <- exp(ov$beta)
m <- merge(res, ov[, c("signature","HR_TCGA_OV")], by = "signature", all.x = TRUE)
m$direction_agree <- sign(log(m$HR)) == sign(log(m$HR_TCGA_OV))

out <- c()
out <- c(out, "=== OV-AU external validation (per-signature univariable Cox) ===")
for (i in seq_len(nrow(m))) {
  r <- m[i, ]
  out <- c(out, sprintf("%-10s n=%3d events=%2d  HR=%.2f (%.2f-%.2f) p=%.3f | TCGA-OV HR=%.2f | agree=%s",
    r$signature, r$n, r$events, r$HR, r$lo, r$hi, r$p, r$HR_TCGA_OV, r$direction_agree))
}
out <- c(out, sprintf("\nDirection agreement with TCGA OV: %d/%d",
  sum(m$direction_agree, na.rm = TRUE), nrow(m)))
writeLines(out, "results/ovau_validation_summary.txt")
cat(out, sep = "\n")
cat("\nDONE\n")
