# B2 口径的 pooled p 值 + 归因结构预览（判断三档预案）
library(metafor)
res <- read.csv("results/specB/per_cancer_three.csv")
mods <- read.csv("results/specB/moderators_abs.csv")
tmb <- read.csv("results/cancer_tmb.csv")
endoderm <- c("LIHC","COAD","READ","STAD","ESCA","PAAD","CHOL","LUAD","LUSC","THCA","CESC","BLCA","PRAD","THYM")
mesoderm <- c("KIRC","KIRP","KICH","OV","UCS","SARC","MESO","LAML","DLBC","TGCT","ACC","UCEC")
ectoderm <- c("GBM","LGG","SKCM","HNSC","PCPG","UVM","BRCA")
germ <- function(a) ifelse(a %in% endoderm, "endoderm", ifelse(a %in% ectoderm, "ectoderm", "mesoderm"))
res$germ <- germ(res$cancer)
res <- merge(res, mods, by = "cancer", all.x = TRUE)
res <- merge(res, tmb, by = "cancer", all.x = TRUE)

sigs <- c("CD8_cyto","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")
out <- c()
for (nm in sigs) {
  d <- res[res$signature == nm & res$version == "B2" & is.finite(res$beta) & res$se > 0, ]
  m <- rma(yi = beta, sei = se, data = d, method = "REML")
  m_p <- rma(yi = beta, sei = se, mods = ~ prolif_abs, data = d, method = "REML")
  m_e <- rma(yi = beta, sei = se, mods = ~ estimate_abs, data = d, method = "REML")
  m_g <- rma(yi = beta, sei = se, mods = ~ germ, data = d, method = "REML")
  m_t <- tryCatch(rma(yi = beta, sei = se, mods = ~ median_TMB, data = d, method = "REML"), error = function(e) NULL)
  line <- sprintf("%s: pooled HR=%.3f (%.3f-%.3f) p=%.3f | prolif R2=%.1f%%(p=%.3f) germ R2=%.1f%%(p=%.3f) est R2=%.1f%%(p=%.3f) tmb R2=%s",
    nm, exp(as.numeric(m$b)), exp(as.numeric(m$ci.lb)), exp(as.numeric(m$ci.ub)), m$pval,
    m_p$R2, m_p$QMp, m_g$R2, m_g$QMp, m_e$R2, m_e$QMp,
    if (is.null(m_t)) "NA" else sprintf("%.1f%%(p=%.3f)", m_t$R2, m_t$QMp))
  out <- c(out, line)
}
writeLines(out, "results/specB/B2_attribution_preview.txt")
cat(out, sep = "\n")
