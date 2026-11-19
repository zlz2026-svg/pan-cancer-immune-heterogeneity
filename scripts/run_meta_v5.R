# M3：单变量 meta 回归定稿（Knapp-Hartung + FDR + 系数）
library(metafor)

res <- read.csv("results/specB/per_cancer_three.csv")
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
mod_names <- c("prolif_abs","estimate_abs","germ","median_TMB")

# 收集所有单变量 meta 回归的 p 值（供 FDR）
all_p <- c()
rows <- list()

for (nm in sigs) {
  d <- res[res$signature == nm & res$version == "B2" & is.finite(res$beta) & res$se > 0, ]
  if (nrow(d) < 8) next
  for (mn in mod_names) {
    fit <- tryCatch(rma(yi = beta, sei = se, mods = as.formula(paste("~", mn)),
                        data = d, method = "REML", test = "knha"), error = function(e) NULL)
    if (is.null(fit)) next
    # 提取系数（连续变量取第2个，germ 取各水平）
    b <- fit$b
    ci_lb <- fit$ci.lb; ci_ub <- fit$ci.ub
    p_vals <- fit$pval
    # 对连续 moderator，报告其系数；对 germ，报告 F 检验 p 即可
    if (mn %in% c("prolif_abs","estimate_abs","median_TMB")) {
      coef_line <- sprintf("beta=%.4f (%.4f, %.4f)", b[2], ci_lb[2], ci_ub[2])
      p_repr <- p_vals[2]
    } else {
      coef_line <- sprintf("endoderm=%.4f mesoderm=%.4f ectoderm=%.4f", b[1], b[2], b[3])
      p_repr <- fit$QMp   # 整体 F 检验 p
    }
    all_p <- c(all_p, p_repr)
    rows[[length(rows)+1]] <- data.frame(
      signature = nm, moderator = mn, R2 = fit$R2, tau2 = fit$tau2,
      beta = coef_line, p = p_repr, QMp = fit$QMp, stringsAsFactors = FALSE)
  }
}
resdf <- do.call(rbind, rows)

# FDR 校正（BH）
resdf$fdr <- p.adjust(resdf$p, method = "BH")

write.csv(resdf, "results/specB/meta_regression_v5.csv", row.names = FALSE)

# 文本摘要
out <- c()
for (nm in sigs) {
  dd <- resdf[resdf$signature == nm, ]
  out <- c(out, sprintf("=== %s (B2, k=%d) ===", nm,
        nrow(res[res$signature==nm & res$version=="B2" & is.finite(res$beta) & res$se>0, ])))
  for (i in seq_len(nrow(dd))) {
    out <- c(out, sprintf("  %-14s R2=%5.1f%%  %s  p=%.3f  FDR=%.3f",
        dd$moderator[i], dd$R2[i], dd$beta[i], dd$p[i], dd$fdr[i]))
  }
}
writeLines(out, "results/specB/meta_regression_v5.txt")
cat(out, sep = "\n")
cat("\nDONE\n")
