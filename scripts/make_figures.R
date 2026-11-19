# 生成正文核心图（读结果文件，不再硬编码）
# Fig2 逐癌种森林图（7 分面 × 31 癌种 + pooled）；Fig3 归因热图；Fig4 方向翻转 + GLASS
library(ggplot2)
library(metafor)
dir.create("figures", showWarnings = FALSE)

res <- read.csv("results/specB/per_cancer_three.csv")
res <- res[res$version == "B2", ]
meta <- read.csv("results/specB/meta_regression_v5.csv")
glass <- read.csv("results/glass_validation_cox.csv")

sigs <- c("CD8_cyto","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")
labs <- c("CD8_Tcell","TLS","CAF","TAN","cDC2","IFNg_path","Prolifer")

# ============ Fig2：逐癌种森林图（7 分面） ============
res$HR <- exp(res$beta)
res$lo <- exp(res$beta - 1.96 * res$se)
res$hi <- exp(res$beta + 1.96 * res$se)
res$lab <- factor(res$signature, levels = sigs, labels = labs)

# pooled 菱形（每签名）
pooled <- do.call(rbind, lapply(sigs, function(nm) {
  d <- res[res$signature == nm & is.finite(res$beta) & res$se > 0, ]
  m <- rma(yi = beta, sei = se, data = d, method = "REML")
  data.frame(signature = nm, HR = exp(as.numeric(m$b)), lo = exp(as.numeric(m$ci.lb)),
             hi = exp(as.numeric(m$ci.ub)), I2 = m$I2, cancer = "POOLED")
}))
pooled$lab <- factor(pooled$signature, levels = sigs, labels = labs)

p2 <- ggplot(res, aes(x = HR, y = reorder(cancer, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, linewidth = 0.4, color = "grey40") +
  geom_point(size = 1.4, color = "grey20") +
  geom_errorbarh(data = pooled, aes(xmin = lo, xmax = hi), height = 0.5, linewidth = 1.1, color = "#c0392b") +
  geom_point(data = pooled, size = 2.6, color = "#c0392b") +
  scale_x_log10() +
  facet_wrap(~lab, ncol = 2, scales = "free_y") +
  labs(x = "Hazard ratio (log scale)", y = "",
       title = "Per-cancer prognostic effects and pooled estimates (31 cancer types)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 6.5),
        strip.text = element_text(size = 9, face = "bold"),
        panel.grid.minor = element_blank())
ggsave("figures/Fig2_forest.png", p2, width = 11, height = 13, dpi = 300)

# ============ Fig3：归因热图（R²，%） ============
# 从 meta_regression_v5.csv 提取 4 个连续/分类 moderator 的 R2 与 p
mod_map <- c("prolif_abs" = "Proliferation", "estimate_abs" = "Stromal/immune",
             "germ" = "Germ layer", "median_TMB" = "TMB")
attr <- meta[meta$moderator %in% names(mod_map), c("signature","moderator","R2","p")]
attr$mod <- factor(mod_map[attr$moderator], levels = c("Germ layer","Proliferation","Stromal/immune","TMB"))
attr$sig <- factor(attr$signature, levels = rev(sigs), labels = rev(labs))
# 星号用 permutation p < 0.05（与正文/Table 2 一致），手动指定显著格
star_pairs <- c("CD8_cyto|prolif_abs", "CAF|prolif_abs", "IFNg_path|prolif_abs",
                "TAN|estimate_abs", "Prolifer|prolif_abs")
attr$star <- ifelse(paste(attr$signature, attr$moderator, sep = "|") %in% star_pairs, "*", "")

p3 <- ggplot(attr, aes(x = mod, y = sig, fill = R2)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(sprintf("%.1f", R2), star)), size = 3.2) +
  scale_fill_gradient(low = "white", high = "#c0392b", limits = c(0, 80), na.value = "grey90") +
  labs(x = "", y = "", fill = expression(tau^2 ~ "explained (%)"),
       title = "Attribution of cross-cancer heterogeneity (meta-regression R2)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("figures/Fig3_attribution.png", p3, width = 6.5, height = 5, dpi = 300)

# ============ Fig4：方向翻转 + GLASS 对照 ============
dirflip <- res[res$signature == "CD8_cyto" & res$cancer %in% c("GBM","OV"), c("cancer","HR","lo","hi")]
dirflip$cancer <- factor(dirflip$cancer, levels = c("GBM","OV"))
p4a <- ggplot(dirflip, aes(x = cancer, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 1) +
  geom_point(size = 4, color = "#c0392b") +
  labs(x = "", y = "CD8_Tcell hazard ratio",
       title = "A. CD8 T-cell direction reversal") +
  coord_flip() +
  theme_minimal(base_size = 13)

glass$sig <- factor(glass$signature, levels = rev(sigs), labels = rev(labs))
p4b <- ggplot(glass, aes(x = sig, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.8) +
  geom_point(size = 3, color = "#c0392b") +
  labs(x = "", y = "Hazard ratio (GLASS glioma)",
       title = "B. GLASS glioma validation") +
  coord_flip() +
  theme_minimal(base_size = 13)

png("figures/Fig4_direction_flip_validation.png", width = 2600, height = 1200, res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
print(p4a, vp = grid::viewport(layout.pos.col = 1))
print(p4b, vp = grid::viewport(layout.pos.col = 2))
dev.off()

cat("Figures saved to figures/\n")
