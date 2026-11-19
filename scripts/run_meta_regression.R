# 泛癌免疫标志物异质性归因 — meta 回归（读逐癌种 Cox 结果）
library(metafor)

res <- read.csv("results/per_cancer_cox.csv", stringsAsFactors = FALSE)
res$germ <- factor(res$germ, levels = c("endoderm", "mesoderm", "ectoderm"))

out <- c()
for (sig in unique(res$signature)) {
  d <- res[res$signature == sig & !is.na(res$se) & res$se > 0, ]
  d <- d[is.finite(d$beta) & is.finite(d$se), ]
  if (nrow(d) < 8) next
  # 随机效应 meta
  m <- rma(yi = beta, sei = se, data = d, method = "REML")
  # 单变量 meta 回归：胚层
  m_germ <- rma(yi = beta, sei = se, mods = ~ germ, data = d, method = "REML")
  # 单变量：增殖
  m_prolif <- rma(yi = beta, sei = se, mods = ~ prolif, data = d, method = "REML")
  # 多变量：胚层 + 增殖 + 纯度
  m_full <- tryCatch(rma(yi = beta, sei = se, mods = ~ germ + prolif + purity_proxy, data = d, method = "REML"),
                     error = function(e) NULL)
  out <- c(out, sprintf(
    "=== %s (k=%d) ===\n  pooled HR=%.3f (%.3f-%.3f), I2=%.1f%%, tau2=%.4f\n  germ meta-reg: R2=%.1f%%, QM p=%.4f\n  prolif meta-reg: R2=%.1f%%, QM p=%.4f\n  full model R2: %s\n",
    sig, nrow(d),
    exp(m$b[1]), exp(m$ci.lb), exp(m$ci.ub), m$I2, m$tau2,
    m_germ$R2, m_germ$QMp,
    m_prolif$R2, m_prolif$QMp,
    if (is.null(m_full)) "NA (failed)" else sprintf("%.1f%%", m_full$R2)
  ))
}
writeLines(out, "results/meta_regression_summary.txt")
cat(out, sep = "")
cat("\nDONE\n")
