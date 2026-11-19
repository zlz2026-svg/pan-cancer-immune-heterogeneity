# 置换检验补充（对核心显著模型）
library(metafor)
res <- read.csv("results/specB/per_cancer_three.csv")
mods <- read.csv("results/specB/moderators_abs.csv")
endoderm <- c("LIHC","COAD","READ","STAD","ESCA","PAAD","CHOL","LUAD","LUSC","THCA","CESC","BLCA","PRAD","THYM")
mesoderm <- c("KIRC","KIRP","KICH","OV","UCS","SARC","MESO","LAML","DLBC","TGCT","ACC","UCEC")
ectoderm <- c("GBM","LGG","SKCM","HNSC","PCPG","UVM","BRCA")
germ <- function(a) ifelse(a %in% endoderm, "endoderm", ifelse(a %in% ectoderm, "ectoderm", "mesoderm"))
res$germ <- factor(germ(res$cancer), levels = c("endoderm","mesoderm","ectoderm"))
res <- merge(res, mods, by = "cancer", all.x = TRUE)

# 核心模型：prolif 对 CD8/CAF/IFNg/Prolifer/TLS + estimate 对 TAN
targets <- list(
  c("CD8_cyto","prolif_abs"), c("CAF","prolif_abs"),
  c("IFNg_path","prolif_abs"), c("Prolifer","prolif_abs"),
  c("TLS","prolif_abs"), c("TAN","estimate_abs"), c("TAN","germ"))
out <- c()
for (t in targets) {
  nm <- t[1]; mn <- t[2]
  d <- res[res$signature == nm & res$version == "B2" & is.finite(res$beta) & res$se > 0, ]
  fit <- rma(yi = beta, sei = se, mods = as.formula(paste("~", mn)), data = d, method = "REML")
  pt <- tryCatch(permutest(fit, iter = 2000, exact = FALSE), error = function(e) NULL)
  if (is.null(pt)) { out <- c(out, sprintf("%s ~ %s: permutest FAILED", nm, mn)); next }
  out <- c(out, sprintf("%s ~ %s: perm_p=%.4f (REML p=%.4f)", nm, mn, pt$pval, fit$pval))
}
writeLines(out, "results/specB/permutest.txt")
cat(out, sep = "\n")
