# S4：胚层敏感性（剔 SARC + 合并 GBM+LGG）+ 各胚层组合并效应
# 用法：cd 到项目根目录后 Rscript scripts/run_sens_s4.R
# 产出：results/sens/s4_germ_sensitivity.txt
library(metafor)

dir.create("results/sens", showWarnings = FALSE)

res <- read.csv("results/specB/per_cancer_three.csv")
res <- res[res$version == "B2", ]
mods <- read.csv("results/specB/moderators_abs.csv")
tmb <- read.csv("results/cancer_tmb.csv")

endoderm <- c("LIHC","COAD","READ","STAD","ESCA","PAAD","CHOL","LUAD","LUSC","THCA","CESC","BLCA","PRAD","THYM")
mesoderm <- c("KIRC","KIRP","KICH","OV","UCS","SARC","MESO","LAML","DLBC","TGCT","ACC","UCEC")
ectoderm <- c("GBM","LGG","SKCM","HNSC","PCPG","UVM","BRCA")
germ <- function(a) ifelse(a %in% endoderm, "endoderm", ifelse(a %in% ectoderm, "ectoderm", "mesoderm"))
res$germ <- factor(germ(res$cancer), levels = c("endoderm","mesoderm","ectoderm"))
res <- merge(res, mods, by = "cancer", all.x = TRUE)
res <- merge(res, tmb, by = "cancer", all.x = TRUE)

sigs <- c("CD8_cyto","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")
out <- c("===== S4 胚层敏感性 =====", "")

# ---- 1. 各胚层组合并效应（组内 pooled HR + 组内 k）----
out <- c(out, "【1】各胚层组合并效应（pooled HR, 组内 k）")
for (nm in sigs) {
  d <- res[res$signature == nm & is.finite(res$beta) & res$se > 0, ]
  line <- sprintf("  %-10s", nm)
  for (g in c("endoderm","mesoderm","ectoderm")) {
    dg <- d[d$germ == g, ]
    if (nrow(dg) < 2) { line <- paste0(line, sprintf(" | %s: k=%d (不足)", g, nrow(dg))); next }
    m <- rma(yi = beta, sei = se, data = dg, method = "REML")
    line <- paste0(line, sprintf(" | %s: k=%d HR=%.3f(p=%.3f)", g, nrow(dg), exp(m$b), m$pval))
  }
  out <- c(out, line)
}

# ---- 2. 剔 SARC（k=30）重跑 germ 归因 ----
out <- c(out, "", "【2】剔 SARC 后 germ 归因（k=30）vs 主分析（k=31）")
for (nm in sigs) {
  d_full <- res[res$signature == nm & is.finite(res$beta) & res$se > 0, ]
  d_no <- d_full[d_full$cancer != "SARC", ]
  m_full <- rma(yi = beta, sei = se, mods = ~ germ, data = d_full, method = "REML", test = "knha")
  m_no <- rma(yi = beta, sei = se, mods = ~ germ, data = d_no, method = "REML", test = "knha")
  out <- c(out, sprintf("  %-10s germ R2: k=31 %.1f%%(p=%.3f) -> k=30(剔SARC) %.1f%%(p=%.3f)",
                        nm, m_full$R2, m_full$QMp, m_no$R2, m_no$QMp))
}

# ---- 3. 合并 GBM+LGG（胶质瘤单独成组，4 水平）重跑 germ 归因 ----
out <- c(out, "", "【3】合并 GBM+LGG 为胶质瘤组后 germ 归因（4 水平）vs 主分析")
res$germ4 <- as.character(res$germ)
res$germ4[res$cancer %in% c("GBM","LGG")] <- "glioma"
res$germ4 <- factor(res$germ4, levels = c("endoderm","mesoderm","ectoderm","glioma"))
for (nm in sigs) {
  d_full <- res[res$signature == nm & is.finite(res$beta) & res$se > 0, ]
  m_full <- rma(yi = beta, sei = se, mods = ~ germ, data = d_full, method = "REML", test = "knha")
  m4 <- rma(yi = beta, sei = se, mods = ~ germ4, data = d_full, method = "REML", test = "knha")
  out <- c(out, sprintf("  %-10s germ R2: 3水平 %.1f%%(p=%.3f) -> 4水平(GBM+LGG合并) %.1f%%(p=%.3f)",
                        nm, m_full$R2, m_full$QMp, m4$R2, m4$QMp))
}

writeLines(out, "results/sens/s4_germ_sensitivity.txt")
cat(out, sep = "\n")
cat("\nDONE. results/sens/s4_germ_sensitivity.txt\n")
