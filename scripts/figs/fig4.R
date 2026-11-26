# fig4.R -- direction reversal and external validation
#
# a: three reversal pairs, full width, with the signature name carried on the
#    first row of each pair so the pairing no longer has to be inferred from the
#    caption;
# b: GLASS forest;
# c: TCGA glioma vs GLASS as a paired dumbbell.
#
# Why panel c is no longer a scatter with a line of identity: the two axes are not
# commensurate.  TCGA glioma HRs span 1.19-1.46 (0.09 log10 decades), GLASS HRs
# span 1.19-3.40 (0.46 decades).  A line through (1,1) is still *correct* -- it does
# mark equal hazard -- but with 5x more decades on y than on x its apparent slope is
# a property of the axis ranges, not of the data, so the perpendicular distance from
# the line (the thing a reader actually estimates) is meaningless.  Equalising the
# limits is not an option either: x would be squeezed into 14 % of the panel.
# A paired dumbbell states the same two facts -- same direction in 7/7, larger in
# GLASS -- without claiming a geometry the panel cannot support.
source("scripts/figs/00_tokens.R")
source("scripts/figs/01_data.R")
source("scripts/figs/00_legend.R")

W <- 183; H <- 100
H_A <- 45; H_B <- 43; H_LEG <- 10; GAP <- 2      # 45+2+43+10 = 100

f2 <- function(x, d = 2) sprintf(paste0("%.", d, "f"), x + 1e-9)

# ---------------- panel a: reversal pairs --------------------------------
SIG_LAB <- c(CD8_Tcell = "CD8 T cells", TLS = "TLS", IFNg_path = "IFN\u2013\u03b3 pathway")
pairs <- data.frame(
  signature = c("CD8_Tcell", "CD8_Tcell", "TLS", "TLS", "IFNg_path", "IFNg_path"),
  cancer    = c("GBM", "OV", "UVM", "CESC", "LGG", "SARC"),
  ypos      = c(6, 5, 4, 3, 2, 1),
  stringsAsFactors = FALSE
)
pa <- merge(pairs, d[, c("cancer", "signature", "n", "HR", "lo", "hi", "p")],
            by = c("cancer", "signature"), sort = FALSE)
pa <- pa[order(-pa$ypos), ]
stopifnot(nrow(pa) == 6)
pa$dir  <- ifelse(pa$HR < 1, "prot", "risk")
pa$lab1 <- sprintf("%s (%s\u2013%s), %s", f2(pa$HR), f2(pa$lo), f2(pa$hi),
                   ifelse(pa$p < 0.001, "p < 0.001", paste0("p = ", f2(pa$p, 3))))
# the signature name rides on the first row of each pair; no left-hand strip, no
# rotated text, no second plot to keep in register
pa$rowlab <- ifelse(duplicated(pa$signature),
                    sprintf("%s (n = %d)", pa$cancer, pa$n),
                    sprintf("%s\n%s (n = %d)",
                            unname(SIG_LAB[pa$signature]), pa$cancer, pa$n))

pA <- ggplot(pa, aes(y = ypos)) +
  annotate("segment", x = 1, xend = 1, y = 0.45, yend = 6.55,
           linewidth = lw(LW$ref), linetype = "22", colour = PAL$ref) +
  geom_segment(aes(x = lo, xend = hi, y = ypos, yend = ypos),
               linewidth = lw(LW$ci), colour = PAL$ink) +
  geom_point(aes(x = HR, y = ypos, fill = dir), shape = 21, size = 2.4,
             colour = PAL$ink, stroke = 0.32) +
  geom_text(aes(x = 2.78, y = ypos, label = lab1), size = fsp(FS$stat),
            colour = PAL$ink, hjust = 0) +
  scale_fill_manual(values = c(prot = PAL$prot, risk = PAL$risk), guide = "none") +
  scale_x_log10(breaks = c(0.6, 0.8, 1.0, 1.5, 2.0)) +
  scale_y_continuous(limits = c(0.4, 6.6), breaks = pa$ypos, labels = pa$rowlab,
                     expand = c(0, 0)) +
  labs(x = "Hazard ratio (log scale)", y = NULL) +
  coord_cartesian(xlim = c(0.55, 6.8)) +
  theme_npj() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = FS$label, colour = PAL$ink, hjust = 1,
                                   lineheight = 0.92),
        axis.title.x = element_text(size = FS$tick),
        plot.margin = margin(2, 2, 1, 1, unit = "mm"))

# ---------------- panel b: GLASS forest ----------------------------------
gb <- glass
gb <- gb[order(match(gb$signature, SIG_ORDER)), ]
gb$y <- rev(seq_len(nrow(gb)))
gb$lab1 <- sprintf("%s%s (%s\u2013%s)", f2(gb$HR), ifelse(gb$sig, "*", ""),
                   f2(gb$lower), f2(gb$upper))

pB <- ggplot(gb, aes(y = y)) +
  annotate("segment", x = 1, xend = 1, y = 0.4, yend = 7.6,
           linewidth = lw(LW$ref), linetype = "22", colour = PAL$ref) +
  geom_segment(aes(x = lower, xend = upper, y = y, yend = y),
               linewidth = lw(LW$ci), colour = PAL$ink) +
  geom_point(aes(x = HR, y = y, fill = sig), shape = 21, size = 2.2,
             colour = PAL$ink2, stroke = 0.30) +
  geom_text(aes(x = 5.4, y = y, label = lab1), size = fsp(FS$stat),
            colour = PAL$ink, hjust = 0) +
  scale_fill_manual(values = c(`TRUE` = PAL$risk, `FALSE` = PAL$ns), guide = "none") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 5)) +
  scale_y_continuous(limits = c(0.3, 7.7), breaks = gb$y, labels = gb$signature,
                     expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(xlim = c(0.95, 11.5)) +
  theme_npj() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = FS$label, colour = PAL$ink, hjust = 1),
        axis.title.x = element_text(size = FS$tick),
        plot.margin = margin(6, 1, 1, 1, unit = "mm"))

# ---------------- panel c: paired dumbbell -------------------------------
cc <- merge(glioma[, c("signature", "HR", "lo", "hi")],
            glass[, c("signature", "HR", "lower", "upper", "p")],
            by = "signature", suffixes = c(".tcga", ".glass"))
cc <- cc[order(match(cc$signature, SIG_ORDER)), ]
cc$y <- rev(seq_len(nrow(cc)))
cse <- (log(cc$upper) - log(cc$lower)) / (2 * 1.96)
cc$invse <- 1 / cse
cc$ratio <- cc$HR.glass / cc$HR.tcga
stopifnot(nrow(cc) == 7, sum(cc$ratio > 1) == 7)

pC <- ggplot(cc, aes(y = y)) +
  annotate("segment", x = 1, xend = 1, y = 0.4, yend = 7.6,
           linewidth = lw(LW$ref), linetype = "22", colour = PAL$ref) +
  geom_segment(aes(x = HR.tcga, xend = HR.glass, y = y, yend = y),
               linewidth = lw(LW$ci), colour = PAL$ns) +
  geom_point(aes(x = HR.tcga), shape = 21, size = 1.9, fill = "white",
             colour = PAL$ink2, stroke = 0.35) +
  geom_point(aes(x = HR.glass, size = invse), shape = 21, fill = PAL$risk,
             colour = PAL$ink, stroke = 0.30) +
  geom_text(aes(x = 3.95, y = y, label = sprintf("%.1f\u00d7", ratio)),
            size = fsp(FS$foot), colour = PAL$ink2, hjust = 0) +
  scale_size_area(max_size = 2.6, limits = c(0, max(cc$invse)), guide = "none") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3)) +
  scale_y_continuous(limits = c(0.3, 7.7), breaks = cc$y, labels = cc$signature,
                     expand = c(0, 0)) +
  labs(x = "Hazard ratio (log scale)", y = NULL) +
  coord_cartesian(xlim = c(0.95, 5.4)) +
  theme_npj() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = FS$label, colour = PAL$ink, hjust = 1),
        axis.title.x = element_text(size = FS$tick),
        plot.margin = margin(6, 1, 1, 1, unit = "mm"))

# ---------------- legend strip -------------------------------------------
# Three swatch columns, two rows.  Every mark that appears in a panel has exactly
# one legend entry and one meaning: an open circle is TCGA glioma (panel c) and is
# never used for "not significant", which is a pale grey disc in panel b.
sw <- function(cx, cy, fill, stroke_col = PAL$ink, sw = 0.32) {
  geom_polygon(data = mm_disc(cx, cy, 2.9), aes(x = x, y = y),
               fill = fill, colour = stroke_col, linewidth = lw(sw))
}
lg <- function(x, y, label, size = FS$legend, col = PAL$ink) {
  annotate("text", x = x, y = y, label = label, hjust = 0, size = fsp(size), colour = col)
}

p_leg <- mm_canvas(W, H_LEG) +
  sw(4.5, 6.2, PAL$prot) + sw(4.5, 2.4, PAL$risk) +
  sw(70.0, 6.2, PAL$ns, PAL$ink2) +
  sw(112.0, 6.2, "white", PAL$ink2, 0.40) + sw(112.0, 2.4, PAL$risk) +
  lg(7.4, 6.2, "protective (HR < 1)") +
  lg(7.4, 2.4, "risk (HR > 1)") +
  lg(72.9, 6.2, "not significant (b)") +
  lg(114.9, 6.2, "TCGA glioma") +
  lg(114.9, 2.4, "GLASS glioma") +
  lg(140.0, 2.4, "panels c; * p < 0.05 in GLASS", FS$foot, PAL$ink2)

# ---------------- compose -------------------------------------------------
TAG_GP <- gpar(fontsize = FS$panel, fontface = "bold", fontfamily = "Arial")

draw_all <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    4, 1, heights = unit(c(H_A, GAP, H_B, H_LEG), "mm"))))

  print(pA, vp = viewport(layout.pos.row = 1))
  grid.text("a", x = unit(1.6, "mm"), y = unit(1, "npc") - unit(0.4, "mm"),
            just = c("left", "top"), gp = TAG_GP, vp = viewport(layout.pos.row = 1))

  pushViewport(viewport(layout.pos.row = 3))
  pushViewport(viewport(layout = grid.layout(
    1, 2, widths = unit(c(96, 87), "mm"), heights = unit(H_B, "mm"))))
  for (i in 1:2) {
    vp <- viewport(layout.pos.row = 1, layout.pos.col = i)
    print(list(pB, pC)[[i]], vp = vp)
    grid.text(c("b", "c")[i], x = unit(1.6, "mm"),
              y = unit(1, "npc") - unit(0.4, "mm"), just = c("left", "top"),
              gp = TAG_GP, vp = vp)
  }
  popViewport(2)

  print(p_leg, vp = viewport(layout.pos.row = 4))
  popViewport()
}

save_fig_fn(draw_all, "Fig4_reversal_validation", W, H)

cat("Fig4 saved. Reversal pairs:\n")
print(data.frame(row = pa$rowlab, hr = pa$lab1), row.names = FALSE)
cat("\nGLASS / TCGA ratio by signature:\n")
print(data.frame(signature = cc$signature, tcga = round(cc$HR.tcga, 2),
                 glass = round(cc$HR.glass, 2), ratio = round(cc$ratio, 2)),
      row.names = FALSE)
