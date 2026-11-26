# fig3.R -- attribution of between-cancer heterogeneity
#
# a: stepped R2 heatmap, transposed to 4 moderators x 7 signatures so that the
#    seven columns are the same seven columns as Fig 2 (cross-figure reading) and
#    no axis label has to be rotated;
# b-e: weighted meta-regression, 2 x 2 below the heatmap.
#
# Layout uses absolute millimetre units in grid.layout.  The previous version laid
# the two halves out in npc fractions and placed the panel letters at eyeballed npc
# constants, which is how the letter "e" ended up drawn twice (43 % overlap with
# the strip text).  Absolute units remove both failure modes.
source("scripts/figs/00_tokens.R")
source("scripts/figs/01_data.R")
source("scripts/figs/00_legend.R")

W <- 183; H <- 152
H_HEAT <- 60; H_PTS <- 80; H_LEG <- 12          # 60 + 80 + 12 = 152

MOD_LEVELS <- c("Germ layer", "Proliferation", "Stromal\u2013immune", "TMB")

# ---------------- panel a: mm heatmap ------------------------------------
X_GUT   <- 27                                    # left gutter for moderator names
X_RIGHT <- W - 5
CELL_W  <- (X_RIGHT - X_GUT) / length(SIG_ORDER)
CELL_H  <- 12.5
Y_HDR   <- 56.8                                  # signature names
Y_TOP   <- 52.5                                  # top edge of row 1

hm <- mr[mr$moderator %in% names(mod_lab), ]
hm$row <- as.integer(factor(hm$mod_lab, levels = MOD_LEVELS))
hm$col <- as.integer(factor(hm$signature, levels = SIG_ORDER))
stopifnot(nrow(hm) == 28, !any(is.na(hm$row)), !any(is.na(hm$col)))

hm$xmin <- X_GUT + (hm$col - 1) * CELL_W
hm$xmax <- hm$xmin + CELL_W
hm$ymax <- Y_TOP - (hm$row - 1) * CELL_H
hm$ymin <- hm$ymax - CELL_H
hm$fill <- PAL$r2_bins[r2_step(hm$R2)]
hm$txt  <- paste0(ifelse(hm$R2 < 0.05, "0", sprintf("%.1f", hm$R2)), hm$star)
# white on the two darkest steps only: #2171B5 is the first that carries white
# text at >= 4.5:1, so the threshold sits at R2 = 45, not at the mid-bin edge
hm$tcol <- ifelse(hm$R2 >= 45, "white", PAL$ink)

# the dagger cell gets a dashed inset outline instead of a hatch: a hatch needs
# clipping to the cell rectangle and the dagger glyph already carries the meaning
dag <- hm[hm$star == "\u2020", ]
stopifnot(nrow(dag) == 1)

hdr <- data.frame(x = X_GUT + (seq_along(SIG_ORDER) - 0.5) * CELL_W,
                  lab = SIG_ORDER, stringsAsFactors = FALSE)
rowlab <- data.frame(y = Y_TOP - (seq_along(MOD_LEVELS) - 0.5) * CELL_H,
                     lab = MOD_LEVELS, stringsAsFactors = FALSE)

p_heat <- mm_canvas(W, H_HEAT) +
  geom_rect(data = hm, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                           fill = fill), colour = "white", linewidth = lw(0.6)) +
  geom_rect(data = dag, aes(xmin = xmin + 0.5, xmax = xmax - 0.5,
                            ymin = ymin + 0.5, ymax = ymax - 0.5),
            fill = NA, colour = "white", linewidth = lw(0.6), linetype = "22") +
  scale_fill_identity() +
  geom_text(data = hm, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = txt),
            colour = hm$tcol, size = fsp(FS$label)) +
  geom_text(data = hdr, aes(x = x, y = Y_HDR, label = lab),
            size = fsp(FS$label), colour = PAL$ink) +
  geom_text(data = rowlab, aes(x = X_GUT - 3, y = y, label = lab),
            size = fsp(FS$label), colour = PAL$ink, hjust = 1) +
  annotate("segment", x = X_GUT, xend = X_RIGHT, y = Y_TOP + 1.6, yend = Y_TOP + 1.6,
           linewidth = lw(LW$edge), colour = PAL$ink)

# ---------------- panels b-e: weighted meta-regression --------------------
PANELS <- list(
  list(tag = "b", sig = "CD8_Tcell", mod = "prolif_abs",
       xlab = "", ylab = "Hazard ratio (log scale)"),
  list(tag = "c", sig = "IFNg_path", mod = "prolif_abs",
       xlab = "", ylab = ""),
  list(tag = "d", sig = "CAF", mod = "prolif_abs",
       xlab = "Tumour proliferation score", ylab = "Hazard ratio (log scale)"),
  list(tag = "e", sig = "TAN", mod = "estimate_abs",
       xlab = "Stromal\u2013immune content score", ylab = "")
)

mk_scatter <- function(p) {
  dd <- merge(d[d$signature == p$sig, c("cancer", "beta", "se", "HR")],
              mods[, c("cancer", "prolif_abs", "estimate_abs", "germ_layer")],
              by = "cancer")
  dd$x <- dd[[p$mod]]
  m  <- rma(yi = dd$beta, sei = dd$se, mods = ~ x, data = dd, method = "REML")
  gx <- seq(min(dd$x), max(dd$x), length.out = 100)
  pr <- predict(m, newmods = gx)
  ft <- data.frame(x = gx, y = exp(pr$pred), lo = exp(pr$ci.lb), hi = exp(pr$ci.ub))
  rr <- mr[mr$signature == p$sig & mr$moderator == p$mod, ]
  stopifnot(nrow(rr) == 1)
  pst <- ifelse(is.na(rr$perm_p), rr$p, rr$perm_p)

  ggplot(dd, aes(x = x, y = HR)) +
    geom_ribbon(data = ft, aes(x = x, ymin = lo, ymax = hi), inherit.aes = FALSE,
                fill = PAL$ns, alpha = 0.75) +
    geom_line(data = ft, aes(x = x, y = y), inherit.aes = FALSE,
              linewidth = lw(LW$main), colour = PAL$ink) +
    geom_hline(yintercept = 1, linetype = "22", linewidth = lw(LW$ref), colour = PAL$ref) +
    geom_point(aes(fill = germ_layer, shape = germ_layer), size = 2.5,
               colour = PAL$ink2, stroke = 0.28) +
    scale_y_log10(breaks = scales::breaks_log(n = 5),
                  labels = function(b) format(b, trim = TRUE, scientific = FALSE)) +
    scale_fill_manual(values = GERM_COL) +
    scale_shape_manual(values = GERM_SHAPE) +
    labs(x = p$xlab, y = p$ylab,
         title = p$sig,
         subtitle = sprintf("R\u00b2 = %.1f%%, perm p = %.4f", rr$R2, pst)) +
    theme_npj() +
    theme(legend.position = "none",
          # title left-aligned like the subtitle, so the two read as one block;
          # the top margin is the space the panel letter is drawn into
          plot.title = element_text(size = FS$label, face = "bold", colour = PAL$ink,
                                    hjust = 0, margin = margin(0, 0, 1.4, 0, unit = "mm")),
          plot.subtitle = element_text(size = FS$stat, colour = PAL$ink2, hjust = 0,
                                       margin = margin(0, 0, 1.2, 0, unit = "mm")),
          panel.grid.major = element_line(linewidth = lw(0.3), colour = "#F2F2F2"),
          plot.margin = margin(4.2, 4, 1, 1, unit = "mm"))
}

pts <- lapply(PANELS, mk_scatter)

# ---------------- legend strip -------------------------------------------
KEY <- mm_colour_key(x0 = 21, y = 8.4, block_w = 11.5, block_h = 4.0,
                     cols = PAL$r2_bins, labs = R2_BIN_LAB, col_lab = PAL$ink)

germ_sw <- do.call(rbind, lapply(seq_along(GERM_ORDER), function(i) {
  g <- GERM_ORDER[i]
  dd <- mm_swatch(100 + (i - 1) * 25, 8.4, 3.2, GERM_SHAPE[[g]], id = i)
  dd$fill <- unname(GERM_COL[[g]])
  dd
}))
germ_lab <- data.frame(x = 100 + (seq_along(GERM_ORDER) - 1) * 25 + 2.6, y = 8.4,
                       label = GERM_ORDER, stringsAsFactors = FALSE)

p_leg <- mm_canvas(W, H_LEG) +
  geom_rect(data = KEY$rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                  fill = fill), colour = "white", linewidth = lw(0.5)) +
  geom_polygon(data = germ_sw, aes(x = x, y = y, group = id, fill = fill),
               colour = PAL$ink2, linewidth = lw(0.3)) +
  scale_fill_identity() +
  geom_text(data = KEY$labs, aes(x = x, y = y, label = label),
            colour = KEY$labs$colour, size = fsp(FS$foot)) +
  geom_text(data = germ_lab, aes(x = x, y = y, label = label),
            size = fsp(FS$legend), colour = PAL$ink, hjust = 0) +
  annotate("text", x = 1.5, y = 8.4, label = "\u03c4\u00b2 explained (%)",
           size = fsp(FS$legend), colour = PAL$ink, hjust = 0, vjust = 0.5) +
  annotate("text", x = 1.5, y = 2.6,
           label = "\u2020 Knapp\u2013Hartung p < 0.05 without permutation significance      \u2021 self-referential (Prolifer regressed on its own score)",
           size = fsp(FS$foot), colour = PAL$ink2, hjust = 0, vjust = 0.5)

# ---------------- compose -------------------------------------------------
TAG_GP <- gpar(fontsize = FS$panel, fontface = "bold", fontfamily = "Arial")

draw_all <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    3, 1, heights = unit(c(H_HEAT, H_PTS, H_LEG), "mm"))))

  # Every print() gets an explicit vp: print.ggplot only skips its internal
  # grid.newpage() when vp is supplied, and a stray newpage would wipe the stack
  # pushed above ("cannot pop the top-level viewport").
  print(p_heat, vp = viewport(layout.pos.row = 1))
  grid.text("a", x = unit(1.6, "mm"), y = unit(1, "npc") - unit(0.8, "mm"),
            just = c("left", "top"), gp = TAG_GP,
            vp = viewport(layout.pos.row = 1))

  pushViewport(viewport(layout.pos.row = 2))
  pushViewport(viewport(layout = grid.layout(
    2, 2, widths = unit(c(W / 2, W / 2), "mm"),
    heights = unit(c(H_PTS / 2, H_PTS / 2), "mm"))))
  for (i in 1:4) {
    vp <- viewport(layout.pos.row = if (i <= 2) 1 else 2,
                   layout.pos.col = if (i %% 2 == 1) 1 else 2)
    print(pts[[i]], vp = vp)
    grid.text(c("b", "c", "d", "e")[i], x = unit(1.6, "mm"),
              y = unit(1, "npc") - unit(0.6, "mm"), just = c("left", "top"),
              gp = TAG_GP, vp = vp)
  }
  popViewport(2)

  print(p_leg, vp = viewport(layout.pos.row = 3))
  popViewport()
}

save_fig_fn(draw_all, "Fig3_attribution", W, H)

cat("Fig3 saved:", nrow(hm), "R2 cells +", length(pts), "scatter panels\n")
