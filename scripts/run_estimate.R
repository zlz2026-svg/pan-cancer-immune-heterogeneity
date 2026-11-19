# 计算 ESTIMATE 分数（Stromal/Immune/ESTIMATE）并按癌种聚合
# 用法：cd 到项目根目录后 Rscript scripts/run_estimate.R
library(survival)

cat("[1/4] Reading expression matrix...\n")
expr_raw <- read.delim(gzfile("data/tcga_RSEM_gene_tpm.gz"), header=TRUE, check.names=FALSE, row.names=1)
E <- as.matrix(expr_raw)
rownames(E) <- sub("\\..*$", "", rownames(E))  # 去 ENSG 版本号
tumor_cols <- grepl("-01[A-Z]?$", colnames(E))
E <- E[, tumor_cols, drop = FALSE]
cat("  Tumor samples:", ncol(E), "\n")

cat("[2/4] Reading ESTIMATE gene mapping...\n")
gs <- read.delim("data/estimate_genes_ensg.tsv", header=TRUE, stringsAsFactors=FALSE)
stromal_ensg <- gs$ensg[gs$set == "stromal"]
immune_ensg <- gs$ensg[gs$set == "immune"]

cat("[3/4] Computing scores (z-score mean)...\n")
sig_score <- function(ensg_ids) {
  g <- intersect(ensg_ids, rownames(E))
  sub <- E[g, , drop = FALSE]
  z <- t(scale(t(sub)))
  colMeans(z, na.rm = TRUE)
}
stromal <- sig_score(stromal_ensg)
immune <- sig_score(immune_ensg)
estimate <- stromal + immune
cat("  stromal genes matched:", sum(stromal_ensg %in% rownames(E)), "/141\n")
cat("  immune genes matched:", sum(immune_ensg %in% rownames(E)), "/141\n")

cat("[4/4] Aggregating to cancer level...\n")
surv <- read.delim("data/survival_pancanatlas.tsv", header=TRUE, check.names=FALSE)
colnames(surv)[1] <- "sample"
surv <- surv[grepl("-01[A-Z]?$", surv$sample), ]
score_df <- data.frame(sample = colnames(E), stromal = stromal, immune = immune, estimate = estimate)
merged <- merge(score_df, surv, by = "sample", all = FALSE)
agg <- aggregate(cbind(stromal, immune, estimate) ~ merged[["cancer type abbreviation"]],
                 data = merged, FUN = mean)
colnames(agg)[1] <- "cancer"
write.csv(agg, "results/cancer_estimate_scores.csv", row.names = FALSE)
cat("  Cancer-level ESTIMATE scores saved (n cancer =", nrow(agg), ")\n")
print(agg[order(-agg$estimate), ][1:5, ])
cat("DONE\n")
