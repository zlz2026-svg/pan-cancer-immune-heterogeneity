# 01_data.R -- single data entry point for all figures. Fails loudly on drift.
# Run from the project root: source("scripts/figs/01_data.R")
# Figure scripts source 00_tokens.R first; sourcing it here keeps this file runnable
# on its own (run_all.R calls it as a standalone assertion pass).
if (!exists("SIG_MAP")) source("scripts/figs/00_tokens.R")

suppressPackageStartupMessages({
  library(metafor)
})

stop_here <- function(...) stop(paste0("ASSERTION FAILED: ", ...), call. = FALSE)

# ---- 1. per-cancer estimates (B2 = within-cancer z-scored scores) --------
raw <- read.csv("results/specB/per_cancer_three.csv", stringsAsFactors = FALSE)
d <- raw[raw$version == "B2", ]
if (nrow(d) != 217) stop_here("expected 217 B2 rows, got ", nrow(d))
if (length(unique(d$cancer)) != 31) stop_here("expected 31 cancer types")
d$signature <- unname(SIG_MAP[d$signature])
if (any(is.na(d$signature))) stop_here("unmapped signature label")
d$HR <- exp(d$beta)
d$lo <- exp(d$beta - 1.96 * d$se)
d$hi <- exp(d$beta + 1.96 * d$se)
d$sig <- d$p < 0.05
d$invse <- 1 / d$se

# ---- 2. cohort meta ------------------------------------------------------
cm <- read.csv("data/cancer_meta.csv", stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
names(cm)[1] <- "cancer"
if (nrow(cm) != 31) stop_here("cancer_meta must have 31 rows")
if (sum(cm$modelled_n) != 8924) stop_here("modelled n sum != 8924 (", sum(cm$modelled_n), ")")
if (sum(cm$eligible_n) != 9032) stop_here("eligible n sum != 9032 (", sum(cm$eligible_n), ")")
if (sum(cm$deaths) != 2680) stop_here("deaths sum != 2680 (", sum(cm$deaths), ")")
if (!all(c("endoderm", "mesoderm", "ectoderm") %in% cm$germ_layer)) stop_here("germ layer coding")
germ_n <- table(factor(cm$germ_layer, levels = GERM_ORDER))
if (!identical(as.integer(germ_n), c(14L, 10L, 7L))) stop_here("germ counts should be 14/10/7")
# modelled n in cancer_meta must equal n used in the Cox models
chk <- merge(cm[, c("cancer", "modelled_n")], unique(d[, c("cancer", "n")]), by = "cancer")
if (any(chk$modelled_n != chk$n)) stop_here("modelled_n does not match Cox model n")
cm$germ_layer <- factor(cm$germ_layer, levels = GERM_ORDER)

# ---- 3. pooled effects: recompute and reconcile with pooled_comparison ---
pooled_comparison <- read.csv("results/specB/pooled_comparison.csv", stringsAsFactors = FALSE)
pooled <- do.call(rbind, lapply(names(SIG_MAP), function(s) {
  dd <- d[d$signature == SIG_MAP[[s]] & is.finite(d$beta) & d$se > 0, ]
  m <- rma(yi = dd$beta, sei = dd$se, method = "REML")
  data.frame(signature = SIG_MAP[[s]], k = nrow(dd),
             HR = exp(as.numeric(m$b)), lo = exp(m$ci.lb), hi = exp(m$ci.ub),
             p = as.numeric(m$pval), I2 = as.numeric(m$I2), tau2 = as.numeric(m$tau2),
             stringsAsFactors = FALSE)
}))
pooled$signature <- factor(pooled$signature, levels = SIG_ORDER)
pooled <- pooled[order(pooled$signature), ]
csv_key <- c(CD8_Tcell = "CD8_cyto", TLS = "TLS", CAF = "CAF", TAN = "TAN",
             cDC2 = "cDC2", IFNg_path = "IFNg_path", Prolifer = "Prolifer")
for (i in seq_len(nrow(pooled))) {
  k <- csv_key[[as.character(pooled$signature[i])]]
  ref_hr <- pooled_comparison$B2_HR[pooled_comparison$signature == k]
  ref_i2 <- pooled_comparison$B2_I2[pooled_comparison$signature == k]
  if (abs(pooled$HR[i] - ref_hr) > 0.005) stop_here("pooled HR mismatch for ", k)
  if (abs(pooled$I2[i] - ref_i2) > 0.1) stop_here("I2 mismatch for ", k)
}
if (any(pooled$I2 < 51.4 - 0.1) || any(pooled$I2 > 85.7 + 0.1)) stop_here("I2 outside 51.4-85.7")

# ---- 4. significant-cell counts (Fig 2 bottom row) -----------------------
sig_tab <- do.call(rbind, lapply(SIG_ORDER, function(s) {
  g <- d[d$signature == s, ]
  data.frame(signature = s,
             prot = sum(g$sig & g$beta < 0),
             risk = sum(g$sig & g$beta > 0),
             n_sig = sum(g$sig), stringsAsFactors = FALSE)
}))
expected_prot <- c(5, 6, 2, 8, 0, 0, 0)
expected_risk <- c(4, 2, 4, 1, 8, 4, 12)
if (!isTRUE(all.equal(as.numeric(sig_tab$prot), expected_prot))) stop_here("protective counts changed: ", paste(sig_tab$prot, collapse = "/"))
if (!isTRUE(all.equal(as.numeric(sig_tab$risk), expected_risk))) stop_here("risk counts changed: ", paste(sig_tab$risk, collapse = "/"))
if (sum(sig_tab$n_sig) != 56) stop_here("total significant cells != 56")

# ---- 5. out-of-range cells (|log2 HR| > 1) ------------------------------
oor <- d[abs(log2(d$HR)) > 1, c("cancer", "signature", "HR")]
if (nrow(oor) != 10) stop_here("expected 10 cells with |log2 HR| > 1")

# ---- 6. direction-reversal cells quoted in the text ---------------------
rev_cells <- data.frame(
  cancer = c("GBM", "OV", "UVM", "CESC", "LGG", "SARC"),
  signature = c("CD8_Tcell", "CD8_Tcell", "TLS", "TLS", "IFNg_path", "IFNg_path"),
  stringsAsFactors = FALSE
)
rev_cells <- merge(rev_cells, d[, c("cancer", "signature", "HR", "lo", "hi", "p")],
                   by = c("cancer", "signature"), sort = FALSE)
rev_cells$expected <- c(1.21, 0.88, 1.69, 0.71, 1.74, 0.66)
if (any(abs(rev_cells$HR - rev_cells$expected) > 0.006)) {
  stop_here("reversal HR no longer matches manuscript: ",
            paste(rev_cells$cancer, sprintf("%.2f", rev_cells$HR), collapse = " "))
}

# ---- 7. meta-regression attribution -------------------------------------
mr <- read.csv("results/specB/meta_regression_v5.csv", stringsAsFactors = FALSE)
mod_lab <- c(germ = "Germ layer", prolif_abs = "Proliferation",
             estimate_abs = "Stromal\u2013immune", median_TMB = "TMB")
mr$signature <- unname(SIG_MAP[mr$signature])
mr$mod_lab <- unname(mod_lab[mr$moderator])
if (nrow(mr) != 28) stop_here("expected 28 meta-regression rows")
if (max(mr$R2) < 73.6 || max(mr$R2) > 73.8) stop_here("max R2 should be 73.7 (Prolifer, self-referential)")

# permutation p from permutest.txt (single run, fixed seed; one row per pair)
perm_lines <- readLines("results/specB/permutest.txt", warn = FALSE)
perm <- do.call(rbind, lapply(perm_lines, function(ln) {
  if (!grepl("perm_p=", ln)) return(NULL)
  sig <- trimws(sub(" ~ .*", "", ln))
  mod <- sub("^.* ~ ", "", ln); mod <- sub(":.*", "", mod)
  p <- as.numeric(sub(".*perm_p=", "", sub(" \\(REML.*", "", ln)))
  data.frame(signature = unname(SIG_MAP[sig]), moderator = mod, perm_p = p,
             stringsAsFactors = FALSE)
}))
mr <- merge(mr, perm, by = c("signature", "moderator"), all.x = TRUE)
mr$star_perm <- ifelse(!is.na(mr$perm_p) & mr$perm_p < 0.05, "*", "")
mr$star <- mr$star_perm
mr$star <- ifelse(mr$signature == "Prolifer" & mr$moderator == "prolif_abs", "\u2021", mr$star)
mr$star <- ifelse(mr$signature == "TAN" & mr$moderator == "germ" & mr$p < 0.05, "\u2020", mr$star)
star_tab <- mr[mr$star_perm == "*", c("signature", "moderator")]
expected_star <- data.frame(
  signature = c("CD8_Tcell", "TLS", "CAF", "IFNg_path", "TAN", "Prolifer"),
  moderator = c("prolif_abs", "prolif_abs", "prolif_abs", "prolif_abs", "estimate_abs", "prolif_abs"),
  stringsAsFactors = FALSE
)
if (!setequal(paste0(star_tab$signature, star_tab$moderator),
              paste0(expected_star$signature, expected_star$moderator))) {
  stop_here("permutation stars moved: ", paste(star_tab$signature, star_tab$moderator, collapse = "; "))
}

# ---- 8. moderators (sample-level axis for Fig 3 scatter) ----------------
mods <- read.csv("results/specB/moderators_abs.csv", stringsAsFactors = FALSE)
mods <- mods[mods$cancer %in% cm$cancer, ]
if (nrow(mods) != 31) stop_here("moderators after dropping TGCT should be 31 rows")
mods <- merge(mods, cm[, c("cancer", "germ_layer")], by = "cancer")

# ---- 9. proportional-hazards violations ---------------------------------
ph <- read.csv("results/specB/ph_check_review.csv", stringsAsFactors = FALSE)
ph$signature <- unname(SIG_MAP[ph$sig])
ph$ph_viol <- ph$score_zph_p < 0.05
if (sum(ph$ph_viol) != 21) stop_here("expected 21 signature-term PH violations, got ", sum(ph$ph_viol))

# ---- 10. external validation (GLASS) ------------------------------------
glass <- read.csv("results/glass_validation_cox.csv", stringsAsFactors = FALSE)
glass$signature <- unname(SIG_MAP[glass$signature])
glass$sig <- glass$p < 0.05
if (nrow(glass) != 7) stop_here("GLASS should have 7 signatures")
if (sum(glass$sig) != 5) stop_here("GLASS should have 5/7 significant")

# ---- 11. glioma-pooled estimates (cache once, never recompute at plot time)
dir.create("results/fig_inputs", showWarnings = FALSE)
glioma_file <- "results/fig_inputs/fig4c_glioma.csv"
if (!file.exists(glioma_file)) {
  gl <- do.call(rbind, lapply(SIG_ORDER, function(s) {
    dd <- d[d$signature == s & d$cancer %in% c("GBM", "LGG") & is.finite(d$beta) & d$se > 0, ]
    m <- rma(yi = dd$beta, sei = dd$se, method = "REML")
    data.frame(signature = s, HR = exp(as.numeric(m$b)), lo = exp(m$ci.lb),
               hi = exp(m$ci.ub), k = nrow(dd), stringsAsFactors = FALSE)
  }))
  write.csv(gl, glioma_file, row.names = FALSE)
}
glioma <- read.csv(glioma_file, stringsAsFactors = FALSE)
if (nrow(glioma) != 7) stop_here("glioma cache must have 7 rows")

write.csv(cm, "results/fig_inputs/cancer_meta.csv", row.names = FALSE)
write.csv(sig_tab, "results/fig_inputs/sig_counts.csv", row.names = FALSE)

message("01_data.R OK: 217 estimates, 31 cancers, 56 significant cells, 21 PH violations, ",
        "pooled reconciled, 7/7 GLASS direction, glioma cache ready.")
