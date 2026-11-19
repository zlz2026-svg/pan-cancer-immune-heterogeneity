# M2 归因部分（复用缓存的 prolif_alt，避免重读矩阵）
# 用法：cd 到项目根目录后 Rscript scripts/run_sens_m2_attribution.R
library(metafor)

res <- read.csv("results/specB/per_cancer_three.csv")
mods <- read.csv("results/specB/moderators_abs.csv")
tmb <- read.csv("results/cancer_tmb.csv")
pa <- read.csv("results/sens/prolif_alt.csv")
mods$prolif_alt <- pa$prolif_alt[match(mods$cancer, pa$cancer)]

endoderm <- c("LIHC","COAD","READ","STAD","ESCA","PAAD","CHOL","LUAD","LUSC","THCA","CESC","BLCA","PRAD","THYM")
mesoderm <- c("KIRC","KIRP","KICH","OV","UCS","SARC","MESO","LAML","DLBC","TGCT","ACC","UCEC")
ectoderm <- c("GBM","LGG","SKCM","HNSC","PCPG","UVM","BRCA")
germ <- function(a) ifelse(a %in% endoderm, "endoderm", ifelse(a %in% ectoderm, "ectoderm", "mesoderm"))
res$germ <- factor(germ(res$cancer), levels = c("endoderm","mesoderm","ectoderm"))
res <- merge(res, mods, by = "cancer", all.x = TRUE)
res <- merge(res, tmb, by = "cancer", all.x = TRUE)

sigs <- c("CD8_cyto","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")
out <- c()
out <- c(out, sprintf("prolif_abs 与 prolif_alt 的 Pearson r = %.4f",
                      cor(res$prolif_abs, res$prolif_alt, use = "complete.obs")), "")
for (nm in sigs) {
  d <- res[res$signature == nm & res$version == "B2" & is.finite(res$beta) & res$se > 0, ]
  if (nrow(d) < 8) next
  m_abs <- rma(yi = beta, sei = se, mods = ~ prolif_abs, data = d, method = "REML", test = "knha")
  m_alt <- rma(yi = beta, sei = se, mods = ~ prolif_alt, data = d, method = "REML", test = "knha")
  pt_alt <- tryCatch(permutest(rma(yi = beta, sei = se, mods = ~ prolif_alt, data = d, method = "REML"),
                               iter = 2000, exact = FALSE), error = function(e) NULL)
  perm_p <- if (is.null(pt_alt)) NA_real_ else as.numeric(tail(pt_alt$pval, 1))
  out <- c(out, sprintf("%-10s prolif_abs R2=%5.1f%% (KH p=%.3f) | prolif_alt R2=%5.1f%% (KH p=%.3f, perm_p=%s)",
                        nm, m_abs$R2, m_abs$QMp, m_alt$R2, m_alt$QMp,
                        if (is.na(perm_p)) "FAIL" else sprintf("%.4f", perm_p)))
}
writeLines(out, "results/sens/m2_prolif_alt.txt")
cat(out, sep = "\n")
cat("\nDONE\n")
