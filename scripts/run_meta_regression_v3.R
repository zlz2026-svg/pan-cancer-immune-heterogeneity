# 泛癌免疫异质性归因 meta 回归 v3（含 TMB）
library(metafor)

res <- read.csv("results/per_cancer_cox.csv", stringsAsFactors = FALSE)
est <- read.csv("results/cancer_estimate_scores.csv", stringsAsFactors = FALSE)
tmb <- read.csv("results/cancer_tmb.csv", stringsAsFactors = FALSE)
res <- merge(res, est, by = "cancer", all.x = TRUE)
res <- merge(res, tmb, by = "cancer", all.x = TRUE)
res$germ <- factor(res$germ, levels = c("endoderm", "mesoderm", "ectoderm"))

out <- c()
for (sig in unique(res$signature)) {
  d <- res[res$signature == sig & !is.na(res$se) & res$se > 0, ]
  d <- d[is.finite(d$beta) & is.finite(d$se), ]
  if (nrow(d) < 8) next
  m <- rma(yi = beta, sei = se, data = d, method = "REML")
  m_germ <- rma(yi = beta, sei = se, mods = ~ germ, data = d, method = "REML")
  m_prol <- rma(yi = beta, sei = se, mods = ~ prolif, data = d, method = "REML")
  m_est  <- rma(yi = beta, sei = se, mods = ~ estimate, data = d, method = "REML")
  m_tmb  <- tryCatch(rma(yi = beta, sei = se, mods = ~ median_TMB, data = d, method = "REML"), error = function(e) NULL)
  m_full <- tryCatch(rma(yi = beta, sei = se, mods = ~ germ + prolif + estimate + median_TMB, data = d, method = "REML"), error = function(e) NULL)
  out <- c(out, sprintf(
    "=== %s (k=%d) ===\n  I2=%.1f%%, tau2=%.4f\n  germ R2=%.1f%% (p=%.3f) | prolif R2=%.1f%% (p=%.3f) | ESTIMATE R2=%.1f%% | TMB R2=%s (p=%s)\n  full R2=%s\n",
    sig, nrow(d), m$I2, m$tau2,
    m_germ$R2, m_germ$QMp, m_prol$R2, m_prol$QMp, m_est$R2,
    if (is.null(m_tmb)) "NA" else sprintf("%.1f%%", m_tmb$R2),
    if (is.null(m_tmb)) "NA" else sprintf("%.3f", m_tmb$QMp),
    if (is.null(m_full)) "NA" else sprintf("%.1f%%", m_full$R2)
  ))
}
writeLines(out, "results/meta_regression_v4.txt")
cat(out, sep = "")
cat("\nDONE\n")
