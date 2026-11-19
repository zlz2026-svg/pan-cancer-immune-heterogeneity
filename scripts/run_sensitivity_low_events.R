# M5：敏感性分析——剔除低事件癌种（P-a），对比 k=31 vs k=26
library(metafor)
res <- read.csv("results/specB/per_cancer_three.csv")
mods <- read.csv("results/specB/moderators_abs.csv")
endoderm <- c("LIHC","COAD","READ","STAD","ESCA","PAAD","CHOL","LUAD","LUSC","THCA","CESC","BLCA","PRAD","THYM")
mesoderm <- c("KIRC","KIRP","KICH","OV","UCS","SARC","MESO","LAML","DLBC","TGCT","ACC","UCEC")
ectoderm <- c("GBM","LGG","SKCM","HNSC","PCPG","UVM","BRCA")
germ <- function(a) ifelse(a %in% endoderm, "endoderm", ifelse(a %in% ectoderm, "ectoderm", "mesoderm"))
res$germ <- factor(germ(res$cancer), levels = c("endoderm","mesoderm","ectoderm"))
res <- merge(res, mods, by = "cancer", all.x = TRUE)

# 低事件癌种（事件数 ≤10，据 survival 数据）：PCPG 6、DLBC 9、THYM 8、PRAD 10、KICH 10
low_event <- c("PCPG","DLBC","THYM","PRAD","KICH")

sigs <- c("CD8_cyto","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")
out <- c("===== k=31（主分析，B2）=====")
for (nm in sigs) {
  d <- res[res$signature == nm & res$version == "B2" & is.finite(res$beta) & res$se > 0, ]
  m <- rma(yi = beta, sei = se, data = d, method = "REML")
  m_p <- rma(yi = beta, sei = se, mods = ~ prolif_abs, data = d, method = "REML", test = "knha")
  out <- c(out, sprintf("%s: k=%d pooled HR=%.3f p=%.3f | prolif R2=%.1f%% p=%.3f",
                        nm, nrow(d), exp(m$b), m$pval, m_p$R2, m_p$QMp))
}

out <- c(out, "", "===== k=26（剔除 PCPG/DLBC/THYM/PRAD/KICH，B2）=====")
for (nm in sigs) {
  d <- res[res$signature == nm & res$version == "B2" & is.finite(res$beta) & res$se > 0 & !(res$cancer %in% low_event), ]
  m <- rma(yi = beta, sei = se, data = d, method = "REML")
  m_p <- rma(yi = beta, sei = se, mods = ~ prolif_abs, data = d, method = "REML", test = "knha")
  out <- c(out, sprintf("%s: k=%d pooled HR=%.3f p=%.3f | prolif R2=%.1f%% p=%.3f",
                        nm, nrow(d), exp(m$b), m$pval, m_p$R2, m_p$QMp))
}
writeLines(out, "results/specB/sensitivity_low_events.txt")
cat(out, sep = "\n")
