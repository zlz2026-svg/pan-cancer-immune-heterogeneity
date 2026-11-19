# M2：替代增殖度量敏感性（prolif_alt vs prolif_abs）
# 目的：检验"增殖归因"是否依赖 Prolifer 签名自指。prolif_alt 用 Hallmark E2F/G2M（剔除 Prolifer 8 基因）绝对表达均值。
# 用法：cd 到项目根目录后 Rscript scripts/run_sens_m2.R
# 产出：results/sens/m2_prolif_alt.txt
library(survival)
library(metafor)

dir.create("results/sens", showWarnings = FALSE)

# ============ 1. 读表达矩阵 ============
cat("[1/4] Reading expression matrix...\n")
expr_raw <- read.delim(gzfile("data/tcga_RSEM_gene_tpm.gz"), header = TRUE, check.names = FALSE, row.names = 1)
E <- as.matrix(expr_raw)
rownames(E) <- sub("\\..*$", "", rownames(E))
tumor_cols <- grepl("-01[A-Z]?$", colnames(E))
E <- E[, tumor_cols, drop = FALSE]

surv <- read.delim("data/survival_pancanatlas.tsv", header = TRUE, check.names = FALSE)
colnames(surv)[1] <- "sample"
surv <- surv[grepl("-01[A-Z]?$", surv$sample), ]
sample_to_cancer <- surv[["cancer type abbreviation"]][match(colnames(E), surv$sample)]
names(sample_to_cancer) <- colnames(E)

# ============ 2. 读替代增殖基因集，算 prolif_alt ============
cat("[2/4] Computing prolif_alt (Hallmark E2F/G2M minus Prolifer8)...\n")
alt <- read.delim("data/hallmark_prolif_alt_ensg.tsv", stringsAsFactors = FALSE)
alt_ensg <- unique(alt$ensg)
alt_ensg <- alt_ensg[alt_ensg %in% rownames(E)]
cat("  alternative proliferation genes mapped to matrix:", length(alt_ensg), "\n")

cancers <- sort(unique(sample_to_cancer[!is.na(sample_to_cancer)]))
prolif_alt <- sapply(cancers, function(ca) {
  idx <- which(sample_to_cancer == ca)
  mean(E[alt_ensg, idx], na.rm = TRUE)   # log2(TPM+0.001) 绝对均值
})
# 缓存 prolif_alt，避免归因出错时重读矩阵
write.csv(data.frame(cancer = names(prolif_alt), prolif_alt = unname(prolif_alt)),
          "results/sens/prolif_alt.csv", row.names = FALSE)

# ============ 3. 读现有结果 + 合并 ============
cat("[3/4] Merging prolif_alt into moderators...\n")
res <- read.csv("results/specB/per_cancer_three.csv")
mods <- read.csv("results/specB/moderators_abs.csv")
tmb <- read.csv("results/cancer_tmb.csv")

mods$prolif_alt <- prolif_alt[match(mods$cancer, names(prolif_alt))]

endoderm <- c("LIHC","COAD","READ","STAD","ESCA","PAAD","CHOL","LUAD","LUSC","THCA","CESC","BLCA","PRAD","THYM")
mesoderm <- c("KIRC","KIRP","KICH","OV","UCS","SARC","MESO","LAML","DLBC","TGCT","ACC","UCEC")
ectoderm <- c("GBM","LGG","SKCM","HNSC","PCPG","UVM","BRCA")
germ <- function(a) ifelse(a %in% endoderm, "endoderm", ifelse(a %in% ectoderm, "ectoderm", "mesoderm"))
res$germ <- factor(germ(res$cancer), levels = c("endoderm","mesoderm","ectoderm"))
res <- merge(res, mods, by = "cancer", all.x = TRUE)
res <- merge(res, tmb, by = "cancer", all.x = TRUE)

# 两个增殖度量的相关性（应正相关，非完全一致）
cor_abs_alt <- cor(res$prolif_abs, res$prolif_alt, use = "complete.obs")

# ============ 4. 归因对比 ============
cat("[4/4] Running meta-regression comparison...\n")
sigs <- c("CD8_cyto","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")
out <- c()
out <- c(out, sprintf("prolif_abs 与 prolif_alt 的 Pearson r = %.4f", cor_abs_alt), "")
for (nm in sigs) {
  d <- res[res$signature == nm & res$version == "B2" & is.finite(res$beta) & res$se > 0, ]
  if (nrow(d) < 8) next
  # prolif_abs（KH）
  m_abs <- rma(yi = beta, sei = se, mods = ~ prolif_abs, data = d, method = "REML", test = "knha")
  # prolif_alt（KH）
  m_alt <- rma(yi = beta, sei = se, mods = ~ prolif_alt, data = d, method = "REML", test = "knha")
  # prolif_alt 置换（取斜率系数 p，即 pval 最后一个元素）
  pt_alt <- tryCatch(permutest(rma(yi = beta, sei = se, mods = ~ prolif_alt, data = d, method = "REML"),
                               iter = 2000, exact = FALSE), error = function(e) NULL)
  perm_p <- if (is.null(pt_alt)) NA_real_ else as.numeric(tail(pt_alt$pval, 1))
  out <- c(out, sprintf("%-10s prolif_abs R2=%5.1f%% (KH p=%.3f) | prolif_alt R2=%5.1f%% (KH p=%.3f, perm_p=%s)",
                        nm, m_abs$R2, m_abs$QMp, m_alt$R2, m_alt$QMp,
                        if (is.na(perm_p)) "FAIL" else sprintf("%.4f", perm_p)))
}
writeLines(out, "results/sens/m2_prolif_alt.txt")
cat(out, sep = "\n")
cat("\nDONE. results/sens/m2_prolif_alt.txt\n")
