# fig1.R -- study design and conceptual framework
#
# a: a three-stage workflow band instead of nine grey boxes on a serpentine arrow;
# b: cohort landscape, with the germ-layer grouping anchored by a colour band and
#    the deaths column given a header;
# c: the two candidate explanations as a fork, drawn with the same cell glyphs as
#    Supplementary Fig S1 so the paper uses one visual language for a "cell field".
source("scripts/figs/00_tokens.R")
source("scripts/figs/01_data.R")
source("scripts/figs/00_legend.R")
source("scripts/figs/00_cells.R")

W <- 183; H <- 178
H_A <- 38; H_ROW <- 120; H_LEG <- 10; GAP <- 2      # 38+2+120+2+10 = 172
W_B <- 92; W_C <- 87; W_GAP <- 4

# ---------------- panel a: three-stage workflow band ----------------------
STAGES <- list(
  list(name = "Cohort assembly",
       lines = c("33 cancer types",
                 "\u2212 LAML, \u2212 TGCT",
                 "31 cohorts, 8,924 tumors")),
  list(name = "Estimation",
       lines = c("7 signatures, z-scored",
                 "per-cancer Cox, 217 models",
                 "age, sex, stage")),
  list(name = "Synthesis",
       lines = c("random-effects meta, I\u00b2 51\u201386%",
                 "meta-regression, 4 moderators",
                 "GLASS replication, n = 79"))
)
CARD_W <- 54; CARD_H <- 22
card_x <- c(3, 64.5, 126)
Y_LINE <- 28; Y_NAME <- 33.5
Y_CARD_TOP <- 24.5
Y_LINES <- Y_CARD_TOP - c(5.0, 11.0, 17.0)

cards <- do.call(rbind, lapply(seq_along(STAGES), function(i) {
  data.frame(x0 = card_x[i], x1 = card_x[i] + CARD_W,
             xm = card_x[i] + CARD_W / 2, name = STAGES[[i]]$name,
             stringsAsFactors = FALSE)
}))
lines_df <- do.call(rbind, lapply(seq_along(STAGES), function(i) {
  data.frame(x = card_x[i] + CARD_W / 2, y = Y_LINES, lab = STAGES[[i]]$lines,
             stringsAsFactors = FALSE)
}))
chev <- do.call(rbind, lapply(1:2, function(i) {
  cx <- (card_x[i] + CARD_W + card_x[i + 1]) / 2
  data.frame(x = cx + c(-1.0, 1.3, -1.0), y = Y_LINE + c(1.3, 0, -1.3), id = i)
}))

pA <- mm_canvas(W, H_A) +
  geom_rect(data = cards, aes(xmin = x0, xmax = x1, ymin = Y_CARD_TOP - CARD_H,
                              ymax = Y_CARD_TOP),
            fill = PAL$pale, colour = NA) +
  annotate("segment", x = 3, xend = 177, y = Y_LINE, yend = Y_LINE,
           linewidth = lw(LW$edge), colour = PAL$ink2) +
  geom_polygon(data = chev, aes(x = x, y = y, group = id), fill = PAL$ink2) +
  geom_point(data = cards, aes(x = xm, y = Y_LINE), shape = 21, size = 1.9,
             fill = "white", colour = PAL$ink2, stroke = 0.5) +
  geom_text(data = cards, aes(x = xm, y = Y_NAME, label = name),
            size = fsp(FS$label), colour = PAL$ink, fontface = "bold") +
  geom_text(data = lines_df, aes(x = x, y = y, label = lab),
            size = fsp(FS$label), colour = PAL$ink2)

# ---------------- panel b: cohort landscape ------------------------------
mb <- cm
mb <- mb[order(as.integer(mb$germ_layer), -mb$modelled_n), ]
mb$y <- rev(seq_len(nrow(mb)))
mb$rowlab <- paste0(mb$cancer, ifelse(mb$deaths <= 10, "\u2020", ""))

X_DEATH <- 1250; X_BAND <- c(1370, 1530); X_GNAME <- 1730
gb1 <- do.call(rbind, lapply(GERM_ORDER, function(g) {
  yy <- mb$y[as.character(mb$germ_layer) == g]
  data.frame(germ = g, y0 = min(yy) - 0.45, y1 = max(yy) + 0.45,
             mid = (min(yy) + max(yy)) / 2, stringsAsFactors = FALSE)
}))
gb1$lab <- c(endoderm = "Endoderm", mesoderm = "Mesoderm",
             ectoderm = "Ectoderm")[gb1$germ]
gb1$fill <- unname(GERM_COL[gb1$germ])
stopifnot(nrow(gb1) == 3, !any(is.na(gb1$lab)))

pB <- ggplot(mb) +
  geom_rect(data = gb1, aes(xmin = X_BAND[1], xmax = X_BAND[2], ymin = y0, ymax = y1),
            fill = gb1$fill, colour = NA) +
  geom_segment(aes(x = eligible_n, xend = modelled_n, y = y, yend = y),
               linewidth = lw(LW$thin), colour = PAL$ns) +
  geom_point(aes(x = eligible_n, y = y), shape = 21, size = 1.45,
             fill = "white", colour = PAL$ink2, stroke = 0.42) +
  geom_point(aes(x = modelled_n, y = y), shape = 21, size = 1.45,
             fill = GERM_COL[as.character(mb$germ_layer)], colour = PAL$ink2,
             stroke = 0.28) +
  geom_text(aes(x = X_DEATH, y = y, label = deaths), size = fsp(FS$cell),
            colour = PAL$ink2, hjust = 1) +
  annotate("text", x = X_DEATH, y = 32.4, label = "deaths",
           size = fsp(FS$foot), colour = PAL$ink2, hjust = 1, fontface = "italic") +
  annotate("text", x = X_GNAME, y = gb1$mid, label = gb1$lab, angle = 90,
           size = fsp(FS$foot), colour = PAL$ink2) +
  scale_x_log10(breaks = c(30, 100, 300, 1000), expand = c(0, 0)) +
  coord_cartesian(xlim = c(30, 2100), clip = "off") +
  scale_y_continuous(breaks = mb$y, labels = mb$rowlab, expand = c(0.015, 0.015)) +
  labs(x = "Tumors (log scale)", y = NULL) +
  theme_npj() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = FS$cell, colour = PAL$ink, hjust = 1),
        plot.margin = margin(5, 2, 1, 1, unit = "mm"))

# ---------------- panel c: the two candidate explanations ----------------
# A "field" is a list of cells plus the polygons they render to; keeping the two
# apart lets assert_field() check the table rather than guess from the artwork.
mk_field <- function(xs, ys, types, d, seed) {
  cells <- data.frame(x = xs, y = ys, type = types, stringsAsFactors = FALSE)
  glyphs <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    cell_glyph(cells$x[i], cells$y[i], d, cells$type[i], id = i, seed = seed + i)
  }))
  list(cells = cells, glyphs = glyphs)
}

VIGN_D <- 5.4; VIGN_D2 <- 4.4
f1a <- mk_field(c(11, 22, 33, 16), c(104, 100, 103, 91),
                c("imm", "imm", "imm", "mal"), VIGN_D, 100)
f1b <- mk_field(c(54, 65, 76, 59), c(104, 100, 103, 91),
                c("mal", "mal", "mal", "imm"), VIGN_D, 200)
# branch 2 holds composition constant (three immune, three malignant) and varies
# only the arrangement -- the claim is about architecture, not about abundance
f2a <- mk_field(c(6, 12, 18, 7, 14, 22), c(58, 57, 58, 46, 45, 46),
                c("imm", "imm", "imm", "mal", "mal", "mal"), VIGN_D2, 300)
f2b <- mk_field(c(37, 46, 55, 40, 49, 52), c(58, 52, 46, 46, 58, 40),
                c("imm", "mal", "imm", "mal", "mal", "imm"), VIGN_D2, 400)
f2c <- mk_field(c(63, 65, 64, 73, 79, 83), c(58, 52, 46, 50, 44, 56),
                c("imm", "imm", "imm", "mal", "mal", "mal"), VIGN_D2, 500)

boxes <- data.frame(
  x0 = c(0, 45, 0, 30, 60), x1 = c(41, 86, 28, 58, 86),
  y0 = c(84, 84, 38, 38, 38), y1 = c(108, 108, 62, 62, 62)
)
glyphs <- rbind(f1a$glyphs, f1b$glyphs, f2a$glyphs, f2b$glyphs, f2c$glyphs)

hd  <- data.frame(x = c(0, 0), y = c(114, 67),
                  lab = c("Proliferation\u2013dilution", "Tissue of origin"),
                  stringsAsFactors = FALSE)
sub <- data.frame(x = c(20.5, 65.5), y = 80.5,
                  lab = c("Low proliferation", "High proliferation"),
                  stringsAsFactors = FALSE)
cap <- data.frame(x = c(20.5, 65.5), y = 76,
                  lab = c("Immune signal intact; HR far from 1",
                          "Immune signal diluted; HR closer to 1"),
                  stringsAsFactors = FALSE)
glab <- data.frame(x = c(14, 44, 73.5), y = 33,
                   lab = c("Endoderm", "Mesoderm", "Ectoderm"),
                   stringsAsFactors = FALSE)

pC <- mm_canvas(W_C, H_ROW) +
  geom_rect(data = boxes, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = NA, colour = PAL$ns, linewidth = lw(0.35)) +
  geom_polygon(data = glyphs, aes(x = x, y = y, group = group, fill = part),
               colour = NA) +
  scale_fill_manual(values = CELL_FILL, guide = "none") +
  geom_text(data = hd, aes(x = x, y = y, label = lab), size = fsp(FS$label),
            colour = PAL$ink, fontface = "bold", hjust = 0, vjust = 0) +
  geom_text(data = sub, aes(x = x, y = y, label = lab), size = fsp(FS$stat),
            colour = PAL$ink2) +
  geom_text(data = cap, aes(x = x, y = y, label = lab), size = fsp(FS$foot),
            colour = PAL$ink2) +
  geom_text(data = glab, aes(x = x, y = y, label = lab), size = fsp(FS$stat),
            colour = PAL$ink2) +
  annotate("text", x = 0, y = 23, hjust = 0, vjust = 0.5, size = fsp(FS$foot),
           colour = PAL$ink2,
           label = "Distinct immune architecture; effect not transferable") +
  annotate("rect", xmin = 0, xmax = W_C, ymin = 4, ymax = 15,
           fill = PAL$band, colour = NA) +
  annotate("text", x = W_C / 2, y = 9.5,
           label = "Both tested as univariate moderators",
           size = fsp(FS$label), colour = PAL$ink)

# ---------------- legend strip -------------------------------------------
leg_cells <- rbind(
  cell_glyph(5.0, 5.6, 4.2, "mal", id = 1, seed = 7),
  cell_glyph(52.0, 5.6, 4.2, "imm", id = 2, seed = 8)
)
p_leg <- mm_canvas(W, H_LEG) +
  geom_polygon(data = mm_disc(110.0, 6.2, 2.9), aes(x = x, y = y), fill = "white",
               colour = PAL$ink2, linewidth = lw(0.42)) +
  geom_polygon(data = mm_disc(110.0, 2.4, 2.9), aes(x = x, y = y), fill = PAL$ink2,
               colour = PAL$ink2, linewidth = lw(0.28)) +
  annotate("text", x = 113.0, y = 6.2, hjust = 0, size = fsp(FS$legend),
           colour = PAL$ink, label = "eligible, matched expression and survival") +
  annotate("text", x = 113.0, y = 2.4, hjust = 0, size = fsp(FS$legend),
           colour = PAL$ink, label = "modeled, entered the Cox models") +
  geom_polygon(data = leg_cells, aes(x = x, y = y, group = group, fill = part),
               colour = NA) +
  scale_fill_manual(values = CELL_FILL, guide = "none") +
  annotate("text", x = 8.0, y = 5.6, hjust = 0, size = fsp(FS$legend),
           colour = PAL$ink, label = "malignant cell") +
  annotate("text", x = 55.0, y = 5.6, hjust = 0, size = fsp(FS$legend),
           colour = PAL$ink, label = "immune cell") +
  annotate("text", x = 70.0, y = 5.6, hjust = 0, size = fsp(FS$foot),
           colour = PAL$ink2, label = NOTE_DEATHS) +
  annotate("text", x = 8.0, y = 1.6, hjust = 0, size = fsp(FS$foot),
           colour = PAL$ink2,
           label = "Panels b and c are schematic; the fields in c are not drawn to scale.")

# ---------------- compose -------------------------------------------------
TAG_GP <- gpar(fontsize = FS$panel, fontface = "bold", fontfamily = "Arial")

draw_all <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    5, 1, heights = unit(c(H_A, GAP, H_ROW, GAP, H_LEG), "mm"))))

  print(pA, vp = viewport(layout.pos.row = 1))
  grid.text("a", x = unit(1.6, "mm"), y = unit(1, "npc") - unit(0.6, "mm"),
            just = c("left", "top"), gp = TAG_GP, vp = viewport(layout.pos.row = 1))

  pushViewport(viewport(layout.pos.row = 3))
  pushViewport(viewport(layout = grid.layout(
    1, 3, widths = unit(c(W_B, W_GAP, W_C), "mm"), heights = unit(H_ROW, "mm"))))
  print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text("b", x = unit(1.6, "mm"), y = unit(1, "npc") - unit(0.6, "mm"),
            just = c("left", "top"), gp = TAG_GP,
            vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pC, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
  grid.text("c", x = unit(1.6, "mm"), y = unit(1, "npc") - unit(0.6, "mm"),
            just = c("left", "top"), gp = TAG_GP,
            vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
  popViewport(2)

  print(p_leg, vp = viewport(layout.pos.row = 5))
  popViewport()
}

save_fig_fn(draw_all, "Fig1_design", W, H)

# the four vignettes must not overlap: same diameter, minimum centre distance
for (fl in list(f1a, f1b)) assert_field(fl$cells, VIGN_D, sum(fl$cells$type == "mal"),
                                        sum(fl$cells$type == "imm"))
for (fl in list(f2a, f2b, f2c)) assert_field(fl$cells, VIGN_D2,
                                             sum(fl$cells$type == "mal"),
                                             sum(fl$cells$type == "imm"))
cat("Fig1 saved. Cohorts:", nrow(mb), " deaths sum:", sum(mb$deaths),
    " modelled sum:", sum(mb$modelled_n), "\n")
cat("Fig1 cell fields: d =", VIGN_D, "mm (branch 1),", VIGN_D2, "mm (branch 2)\n")
