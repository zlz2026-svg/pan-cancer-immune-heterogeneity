# 诊断 Cox 循环问题（假签名评分）
library(survival)
surv <- read.delim("data/survival_pancanatlas.tsv", header=TRUE, check.names=FALSE)
colnames(surv)[1] <- "sample"
surv <- surv[grepl("-01[A-Z]?$", surv$sample), ]

# 假签名评分（7个签名，随机数）
set.seed(1)
n <- nrow(surv)
scores <- data.frame(sample = surv$sample, CD8_cyto=rnorm(n), TLS=rnorm(n), CAF=rnorm(n), TAN=rnorm(n),
                     cDC2=rnorm(n), IFNg_path=rnorm(n), Prolifer=rnorm(n))
merged <- merge(scores, surv, by = "sample", all = FALSE)

merged$age <- as.numeric(merged$age_at_initial_pathologic_diagnosis)
merged$OS.time <- as.numeric(merged$OS.time)
merged$OS <- as.numeric(merged$OS)
merged$stage <- as.character(merged$ajcc_pathologic_tumor_stage)
merged$stage_grp <- ifelse(grepl("^Stage I|^Stage II", merged$stage), "I-II",
                    ifelse(grepl("^Stage III|^Stage IV", merged$stage), "III-IV", "NA"))

cat("OS.time NA:", sum(is.na(merged$OS.time)), " | OS.time<=0:", sum(merged$OS.time <= 0, na.rm=TRUE), "\n")

sig_names <- c("CD8_cyto","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")
res_list <- list()
for (ca in unique(merged[["cancer type abbreviation"]])) {
  d <- merged[merged[["cancer type abbreviation"]] == ca, ]
  if (nrow(d) < 30 || sum(d$OS == 1, na.rm = TRUE) < 5) next
  for (sig in sig_names) {
    d2 <- d[!is.na(d[[sig]]) & !is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
    if (nrow(d2) < 25) next
    d2$score <- as.numeric(d2[[sig]])
    fit <- tryCatch(coxph(Surv(OS.time, OS) ~ score + age + gender + stage_grp, data=d2),
                    error=function(e){cat("[ERR]",ca,sig,conditionMessage(e),"\n");NULL})
    if (is.null(fit)) next
    s <- summary(fit)
    res_list[[length(res_list)+1]] <- data.frame(cancer=ca, signature=sig, HR=s$coefficients[1,"exp(coef)"], beta=s$coefficients[1,"coef"], se=s$coefficients[1,"se(coef)"])
  }
}
cat("res_list length:", length(res_list), "\n")
if (length(res_list) > 0) {
  res <- do.call(rbind, res_list)
  cat("res rows:", nrow(res), "\n")
  print(head(res, 5))
}
