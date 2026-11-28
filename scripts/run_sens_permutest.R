# Permutation tests for the key sensitivity models (grade / DSS versions)
# Mirrors run_permutest.R but for results/sens_grade_dss/per_cancer_sens.csv
library(metafor)

res <- read.csv("results/sens_grade_dss/per_cancer_sens.csv")
mods <- read.csv("results/specB/moderators_abs.csv")
res <- merge(res, mods, by = "cancer", all.x = TRUE)

# key sensitivity models to permutation-test
targets <- list(
  c("CD8_cyto","prolif_abs","B2_grade"),
  c("CD8_cyto","prolif_abs","B2_dss"),
  c("IFNg_path","prolif_abs","B2_grade"),
  c("IFNg_path","prolif_abs","B2_dss"),
  c("CAF","prolif_abs","B2_grade"),
  c("TLS","prolif_abs","B2_grade"),
  c("TLS","prolif_abs","B2_dss"),
  c("TAN","estimate_abs","B2_grade"),
  c("TAN","estimate_abs","B2_dss")
)

out <- c()
for (t in targets) {
  nm <- t[1]; mn <- t[2]; ver <- t[3]
  d <- res[res$signature == nm & res$version == ver & is.finite(res$beta) & res$se > 0, ]
  if (nrow(d) < 8) { out <- c(out, sprintf("%s ~ %s (%s): skip (k=%d)", nm, mn, ver, nrow(d))); next }
  fit <- rma(yi = beta, sei = se, mods = as.formula(paste("~", mn)), data = d, method = "REML")
  pt <- tryCatch(permutest(fit, iter = 2000, exact = FALSE), error = function(e) NULL)
  if (is.null(pt)) { out <- c(out, sprintf("%s ~ %s (%s): permutest FAILED", nm, mn, ver)); next }
  # permutest pval: intercept + slope; take the slope (last)
  pslope <- tail(pt$pval, 1)
  out <- c(out, sprintf("%s ~ %s (%s, k=%d): perm_p=%.4f (z p=%.4f)", nm, mn, ver, nrow(d), pslope, tail(fit$pval, 1)))
}
writeLines(out, "results/sens_grade_dss/sens_permutest.txt")
cat(out, sep = "\n")
cat("\nDONE\n")
