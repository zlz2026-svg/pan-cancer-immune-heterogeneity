# Clean meta-regression comparison for sensitivity versions (no matrix read)
# Reads results/sens_grade_dss/per_cancer_sens.csv produced by run_sens_grade_dss.R
library(metafor)

res <- read.csv("results/sens_grade_dss/per_cancer_sens.csv")
mods <- read.csv("results/specB/moderators_abs.csv")
tmb <- read.csv("results/cancer_tmb.csv")
res <- merge(res, mods, by = "cancer", all.x = TRUE)
res <- merge(res, tmb, by = "cancer", all.x = TRUE)

endoderm <- c("LIHC","COAD","READ","STAD","ESCA","PAAD","CHOL","LUAD","LUSC","THCA","CESC","BLCA","PRAD","THYM")
mesoderm <- c("KIRC","KIRP","KICH","OV","UCS","SARC","MESO","LAML","DLBC","TGCT","ACC","UCEC")
ectoderm <- c("GBM","LGG","SKCM","HNSC","PCPG","UVM","BRCA")
germ <- function(a) ifelse(a %in% endoderm, "endoderm", ifelse(a %in% ectoderm, "ectoderm", "mesoderm"))
res$germ <- factor(germ(res$cancer), levels = c("endoderm","mesoderm","ectoderm"))

sigs <- c("CD8_cyto","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")
versions <- c("B2_repro","B2_grade","B2_dss","B2_dss_grade")
mod_names <- c("prolif_abs","estimate_abs","germ","median_TMB")

out <- c()
for (ver in versions) {
  dv <- res[res$version == ver & is.finite(res$beta) & res$se > 0, ]
  out <- c(out, sprintf("===== %s (k=%d) =====", ver, length(unique(dv$cancer))))
  for (nm in sigs) {
    d <- dv[dv$signature == nm, ]
    if (nrow(d) < 8) next
    for (mn in mod_names) {
      fit <- tryCatch(rma(yi = beta, sei = se, mods = as.formula(paste("~", mn)),
                          data = d, method = "REML", test = "knha"), error = function(e) NULL)
      if (is.null(fit)) next
      if (mn %in% c("prolif_abs","estimate_abs","median_TMB")) {
        pval <- fit$pval[2]
      } else {
        pval <- fit$QMp
      }
      out <- c(out, sprintf("  %-10s %-13s R2=%5.1f%% p=%.4f", nm, mn, fit$R2, pval))
    }
  }
}
writeLines(out, "results/sens_grade_dss/sens_meta_compare.txt")
cat(out, sep = "\n")
cat("\nDONE\n")
