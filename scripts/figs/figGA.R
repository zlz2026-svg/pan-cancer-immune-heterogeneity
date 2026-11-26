# figGA.R -- Graphical Abstract (183 x 120 mm, 3:2)
#
# Three-column narrative plus a bottom band, per figures/briefs/01_GraphicalAbstract.md:
#   Column 1 "Data"        -- 31 germ-layer squares (14 / 10 / 7 in three runs)
#   Column 2 "Heterogeneity" -- a two-tone bar with a direction reversal
#   Column 3 "Attribution"  -- five result cards
#   Bottom band             -- the one-line conclusion
#
# Coordinates are mm-canvas (y = 0 at the bottom), so the brief's top-down y
# values are flipped: the bottom band sits at the foot of the canvas.
source("scripts/figs/00_tokens.R")
source("scripts/figs/00_legend.R")

W <- 183; H <- 120

# ---- Column 1: 31 squares (upper left) -----------------------------------
sq <- 2.5; gap <- 1.2
grid_x0 <- 8; grid_y0 <- 110
n_endo <- 14; n_meso <- 10; n_ecto <- 7
types <- c(rep("endoderm", n_endo), rep("mesoderm", n_meso), rep("ectoderm", n_ecto))
sq_df <- do.call(rbind, lapply(seq_along(types), function(i) {
  r <- (i - 1) %/% 6L
  c <- (i - 1) %% 6L
  mm_rect(grid_x0 + c * (sq + gap), grid_y0 - r * (sq + gap) - sq,
          grid_x0 + c * (sq + gap) + sq, grid_y0 - r * (sq + gap))
}))
sq_df$fill <- GERM_COL[types]

c1_base <- grid_y0 - 5 * (sq + gap)
c1_stats <- data.frame(
  x = grid_x0 + 2.5 * (sq + gap),
  y = c1_base - c(5, 9, 13),
  label = c("31 cancer types", "8,924 tumors", "7 signatures"),
  size = FS$axis_title, stringsAsFactors = FALSE
)

# ---- Column 2: two-tone bar (middle) -------------------------------------
bar_x0 <- 63; bar_x1 <- 125; bar_y0 <- 65.5; bar_y1 <- 70.5
bar_mid <- bar_x0 + (bar_x1 - bar_x0) / 2
bar <- data.frame(
  xmin = c(bar_x0, bar_mid + 0.5), xmax = c(bar_mid - 0.5, bar_x1),
  ymin = bar_y0, ymax = bar_y1,
  fill = c(PAL$prot, PAL$risk), stringsAsFactors = FALSE
)
arrow_cy <- (bar_y0 + bar_y1) / 2
arrow_l <- data.frame(x = c(bar_x0, bar_x0 - 2.2, bar_x0),
                      y = c(arrow_cy + 1.6, arrow_cy, arrow_cy - 1.6))
arrow_r <- data.frame(x = c(bar_x1, bar_x1 + 2.2, bar_x1),
                      y = c(arrow_cy + 1.6, arrow_cy, arrow_cy - 1.6))
c2_lab <- data.frame(
  x = (bar_x0 + bar_x1) / 2,
  y = c(bar_y1 + 6.5, bar_y0 - 6.0),
  label = c("I\u00b2 = 51\u201386%", "Direction reversal"),
  size = c(FS$axis_title, FS$label), stringsAsFactors = FALSE
)

# ---- Column 3: five result cards (upper right) ---------------------------
cards <- data.frame(
  name = c("Proliferation", "Stromal\u2013immune content", "Germ layer",
           "Tumor mutational burden", "GLASS validation"),
  value = c("27.6\u201337.5%", "36.5%", "glioma-specific", "4.2% maximum", "7/7 consistent"),
  stringsAsFactors = FALSE
)
card_x0 <- 132; card_x1 <- 177; card_h <- 13; card_gap <- 2
card_y0 <- 108 - (seq_len(nrow(cards)) - 1) * (card_h + card_gap)
card_df <- do.call(rbind, lapply(seq_len(nrow(cards)), function(i) {
  mm_rect(card_x0, card_y0[i] - card_h, card_x1, card_y0[i])
}))
card_text <- do.call(rbind, lapply(seq_len(nrow(cards)), function(i) {
  rbind(
    data.frame(x = card_x0 + 2.5, y = card_y0[i] - 3.2, label = cards$name[i],
               colour = PAL$ink2, size = FS$label, stringsAsFactors = FALSE),
    data.frame(x = card_x0 + 2.5, y = card_y0[i] - 7.6, label = cards$value[i],
               colour = PAL$ink, size = FS$axis_title, stringsAsFactors = FALSE)
  )
}))

# ---- bottom band ----------------------------------------------------------
band <- mm_rect(5, 6, 178, 20)
band_lab <- data.frame(x = (5 + 178) / 2, y = 13,
                       label = "Immune signatures need proliferation adjustment",
                       size = FS$axis_title, stringsAsFactors = FALSE)

# ---- assemble -------------------------------------------------------------
p <- mm_canvas(W, H) +
  geom_rect(data = sq_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = sq_df$fill) +
  geom_text(data = c1_stats, aes(x = x, y = y, label = label),
            size = fsp(c1_stats$size), colour = PAL$ink) +
  geom_rect(data = bar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bar$fill) +
  geom_polygon(data = arrow_l, aes(x = x, y = y), fill = PAL$ink) +
  geom_polygon(data = arrow_r, aes(x = x, y = y), fill = PAL$ink) +
  geom_text(data = c2_lab, aes(x = x, y = y, label = label),
            size = fsp(c2_lab$size), colour = PAL$ink) +
  geom_rect(data = card_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = PAL$pale) +
  geom_text(data = card_text, aes(x = x, y = y, label = label, colour = colour),
            size = fsp(card_text$size), hjust = 0) +
  scale_colour_identity() +
  geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = PAL$pale) +
  geom_text(data = band_lab, aes(x = x, y = y, label = label),
            size = fsp(band_lab$size), colour = PAL$ink)

save_fig(p, "Fig_GA", W, H)
cat("Fig_GA saved.\n")
