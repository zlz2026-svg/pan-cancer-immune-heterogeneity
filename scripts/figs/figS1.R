# figS1.R -- Supplementary Figure S1: immune-signal dilution by proliferation
#
# Two matched fields, a signal bar under each, a legend strip and a footnote.  The
# fields are drawn with the shared cell glyphs (00_cells.R) rather than as flat
# discs: at this size a flat disc reads as a scatter point, not as a cell.
#
# Geometry is deliberately constrained and asserted, because the whole point of the
# panel is a controlled comparison:
#   * both fields are the same 70 x 50 mm area;
#   * every cell has exactly the same outer diameter (assert_diameter);
#   * no two cells touch (assert_field);
#   * composition is 6 + 18 on the left and 26 + 5 on the right.
# Malignant and immune cells differ in nucleus size/position and membrane ripple,
# never in size.
source("scripts/figs/00_tokens.R")
source("scripts/figs/00_legend.R")
source("scripts/figs/00_cells.R")

W <- 183; H <- 95
CELL_D <- 5.2
FIELD_W <- 70; FIELD_H <- 50
Y_FIELD0 <- 30
Y_TITLE <- 88
Y_BAR <- 22; BAR_W <- 66; BAR_H <- 3.2
Y_NOTE <- 16; Y_LEG <- 9; Y_FOOT <- 3

XL <- c(12, 12 + FIELD_W)
XR <- c(101, 101 + FIELD_W)

set.seed(20261120)
low  <- place_cells(6, 18, XL[1], Y_FIELD0, FIELD_W, FIELD_H, CELL_D, seed = 11)
high <- place_cells(26, 5, XR[1], Y_FIELD0, FIELD_W, FIELD_H, CELL_D, seed = 22)

assert_field(low,  CELL_D, 6, 18)
assert_field(high, CELL_D, 26, 5)

poly_low  <- render_field(low,  CELL_D, seed = 1000)
poly_high <- render_field(high, CELL_D, seed = 2000)
poly <- rbind(poly_low, poly_high)
assert_diameter(poly, CELL_D)

fields <- data.frame(x0 = c(XL[1], XR[1]), x1 = c(XL[2], XR[2]),
                     y0 = Y_FIELD0, y1 = Y_FIELD0 + FIELD_H)

bars <- data.frame(x0 = c(XL[1] + (FIELD_W - BAR_W) / 2, XR[1] + (FIELD_W - BAR_W) / 2),
                   fillfrac = c(0.75, 0.16))
bars$x1 <- bars$x0 + BAR_W

labs <- data.frame(x = c(mean(XL), mean(XR)),
                   title = c("Low-proliferation tumour", "High-proliferation tumour"),
                   note  = c("Immune signal intact; HR far from 1",
                             "Immune signal diluted; HR closer to 1"),
                   stringsAsFactors = FALSE)

leg_cells <- rbind(cell_glyph(10.0, Y_LEG, CELL_D, "mal", id = 1, seed = 31),
                   cell_glyph(52.0, Y_LEG, CELL_D, "imm", id = 2, seed = 32))

p <- mm_canvas(W, H) +
  # the field outline delimits the sampled area; without it the cells read as a
  # free-floating scatter rather than as two comparable fields
  geom_rect(data = fields, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = NA, colour = PAL$ns, linewidth = lw(0.35)) +
  geom_polygon(data = poly, aes(x = x, y = y, group = group, fill = part), colour = NA) +
  scale_fill_manual(values = CELL_FILL, guide = "none") +
  geom_rect(data = bars, aes(xmin = x0, xmax = x1, ymin = Y_BAR, ymax = Y_BAR + BAR_H),
            fill = PAL$pale, colour = NA) +
  geom_rect(data = bars, aes(xmin = x0, xmax = x0 + (x1 - x0) * fillfrac,
                             ymin = Y_BAR, ymax = Y_BAR + BAR_H),
            fill = CELL_STYLE$imm[["membrane"]], colour = NA) +
  geom_text(data = labs, aes(x = x, y = Y_TITLE, label = title),
            size = fsp(FS$label), colour = PAL$ink, fontface = "bold") +
  geom_text(data = labs, aes(x = x, y = Y_NOTE, label = note),
            size = fsp(FS$label), colour = PAL$ink2) +
  geom_polygon(data = leg_cells, aes(x = x, y = y, group = group, fill = part), colour = NA) +
  annotate("text", x = 13.1, y = Y_LEG, hjust = 0, size = fsp(FS$legend),
           colour = PAL$ink, label = "Malignant cell") +
  annotate("text", x = 55.1, y = Y_LEG, hjust = 0, size = fsp(FS$legend),
           colour = PAL$ink, label = "Immune cell") +
  annotate("text", x = W / 2, y = Y_FOOT,
           label = "Schematic; not to scale. Association, not causation.",
           size = fsp(FS$foot), colour = PAL$ink2)

save_fig(p, "Supplementary_Figure_S1", W, H)

cat("FigS1 saved:", nrow(low) + nrow(high), "cells (", nrow(low), "left:",
    sum(low$type == "mal"), "mal +", sum(low$type == "imm"), "imm;", nrow(high),
    "right:", sum(high$type == "mal"), "mal +", sum(high$type == "imm"), "imm )\n")
cat("FigS1 diameter check: max radius", round(max(sqrt((poly$x - poly$cx)^2 +
    (poly$y - poly$cy)^2)[poly$part %in% c("mal_membrane", "imm_membrane")]), 3),
    "mm, nominal", CELL_D / 2, "mm\n")
