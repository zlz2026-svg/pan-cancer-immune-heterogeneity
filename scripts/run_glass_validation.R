# GLASS 胶质瘤独立队列 — 免疫签名 OS Cox 外部验证
# 与主分析同套签名基因集 + 同套 z-score 评分
library(survival)

d <- read.csv("results/glass_expr_clin.csv", check.names = FALSE, stringsAsFactors = FALSE)

signatures <- list(
  CD8_cyto  = c("CD8A","CD8B","CD3D","CD3E","GZMA","GZMB","PRF1","IFNG","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX"),
  TLS       = c("CXCL13","CCL19","CCL21","CD79A","MS4A1","LTB","TNFRSF17"),
  CAF       = c("COL1A1","COL1A2","COL3A1","ACTA2","FAP","PDGFRB","POSTN","FN1","MMP2"),
  TAN       = c("CSF3R","S100A8","S100A9","FCGR3B","CXCR2","OLR1"),
  cDC2      = c("CD1C","FCER1A","CLEC10A","CD1E"),
  IFNg_path = c("IFNG","CXCL9","CXCL10","CXCL11","GBP1","STAT1","IRF1","GBP2"),
  Prolifer  = c("MKI67","PCNA","TOP2A","BIRC5","CCNB1","MCM2","AURKA","CDK1")
)

# 签名评分：跨基因 z-score（跨样本）后取均值，与主分析一致
for (nm in names(signatures)) {
  genes <- intersect(signatures[[nm]], colnames(d))
  cat(sprintf("  %s: %d/%d 基因匹配\n", nm, length(genes), length(signatures[[nm]])))
  sub <- as.matrix(d[, genes, drop = FALSE])
  storage.mode(sub) <- "double"
  z <- scale(sub)                             # 每列(基因)跨样本 z-score，与主分析一致
  d[[nm]] <- rowMeans(z, na.rm = TRUE)        # 每个样本取基因均值
}

d$OS <- as.numeric(d$OS)
d$OS.time <- as.numeric(d$OS.time)
d2 <- d[!is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
cat(sprintf("有效样本: %d, 事件: %d\n", nrow(d2), sum(d2$OS == 1)))

# 单变量 Cox
res <- data.frame()
for (nm in names(signatures)) {
  fit <- tryCatch(coxph(Surv(OS.time, OS) ~ d2[[nm]], data = d2), error = function(e) NULL)
  if (!is.null(fit)) {
    s <- summary(fit)
    res <- rbind(res, data.frame(
      signature = nm,
      HR = s$conf.int[1],
      lower = s$conf.int[3],
      upper = s$conf.int[4],
      p = s$coefficients[5],
      stringsAsFactors = FALSE
    ))
  }
}
res$direction <- ifelse(res$HR > 1, "risk", "protective")
write.csv(res, "results/glass_validation_cox.csv", row.names = FALSE)
cat("\n===== GLASS 外部验证 Cox 结果 =====\n")
print(res, row.names = FALSE)
