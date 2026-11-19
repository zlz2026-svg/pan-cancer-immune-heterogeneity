library(metafor)
res <- read.csv("results/specB/per_cancer_three.csv")
sigs <- c("CD8_cyto","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")
out <- c()
for (nm in sigs) {
  d <- res[res$signature == nm & res$version == "B2" & is.finite(res$beta) & res$se > 0, ]
  m <- rma(yi = beta, sei = se, data = d, method = "REML")
  loo <- sapply(seq_len(nrow(d)), function(i) exp(as.numeric(rma(yi = beta, sei = se, data = d[-i, ], method = "REML")$b)))
  out <- c(out, sprintf("%s: full HR=%.3f, loo range=%.3f-%.3f", nm, exp(as.numeric(m$b)), min(loo), max(loo)))
}
writeLines(out, "results/specB/leave_one_out_B2.txt")
cat(out, sep = "\n")
