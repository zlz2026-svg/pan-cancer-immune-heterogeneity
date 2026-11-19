# 留一法敏感性分析：逐个剔除癌种重算 pooled HR
library(metafor)

res <- read.csv("results/per_cancer_cox.csv", stringsAsFactors = FALSE)

out <- c()
for (sig in unique(res$signature)) {
  d <- res[res$signature == sig & !is.na(res$se) & res$se > 0, ]
  d <- d[is.finite(d$beta) & is.finite(d$se), ]
  if (nrow(d) < 8) next

  m_full <- rma(yi = beta, sei = se, data = d, method = "REML")
  full_hr <- exp(as.numeric(m_full$b))
  full_lo <- exp(as.numeric(m_full$ci.lb))
  full_hi <- exp(as.numeric(m_full$ci.ub))

  loo <- sapply(seq_len(nrow(d)), function(i) {
    m <- rma(yi = beta, sei = se, data = d[-i, ], method = "REML")
    exp(as.numeric(m$b))
  })
  # 影响最大的癌种：剔除后 pooled HR 偏离完整值最远
  idx <- which.max(abs(log(loo) - log(full_hr)))

  out <- c(out, sprintf(
    "=== %s (k=%d) ===\n  full pooled HR = %.3f (95%%CI %.3f-%.3f)\n  leave-one-out HR range = %.3f - %.3f\n  most influential: %s (HR without it = %.3f)\n",
    sig, nrow(d), full_hr, full_lo, full_hi, min(loo), max(loo), d$cancer[idx], loo[idx]
  ))
}
writeLines(out, "results/leave_one_out.txt")
cat(out, sep = "")
cat("\nDONE\n")
