# figS2.R -- robustness of the attribution
# a: main vs alternative proliferation metric; b: with vs without stage adjustment;
# c: leave-one-out range of every pooled hazard ratio.
#
# The three data files are parsed and asserted exactly as before -- the parsing
# block below is unchanged, only the rendering is new.  Panels a and b used to be
# ggrepel slopegraphs that flung the seven signature names to the left edge, so
# reading one row meant tracking a line across 40 mm.  They are now dumbbells with
# the names as ordinary row labels.
source("scripts/figs/00_tokens.R")
source("scripts/figs/01_data.R")
source("scripts/figs/00_legend.R")

W <- 183; H <- 110
H_A <- 44; H_C <- 62; GAP <- 4                   # 44 + 4 + 62 = 110

# ---------------- parse sensitivity outputs ------------------------------
m2 <- readLines("results/sens/m2_prolif_alt.txt", warn = FALSE)
m2 <- m2[grepl("prolif_abs R2", m2)]
alt <- do.call(rbind, lapply(m2, function(ln) {
  s <- trimws(sub(" .*", "", ln))
  v <- regmatches(ln, gregexpr("[0-9]+\\.?[0-9]*%", ln))[[1]]
  data.frame(signature = unname(SIG_MAP[s]),
             main = as.numeric(sub("%", "", v[1])),
             alt  = as.numeric(sub("%", "", v[2])),
             perm = as.numeric(sub(".*perm_p=([0-9.]+).*", "\\1", ln)),
             stringsAsFactors = FALSE)
}))
alt$signature <- factor(alt$signature, levels = SIG_ORDER)
alt <- alt[order(alt$signature), ]
stopifnot(nrow(alt) == 7)
stopifnot(abs(alt$main[alt$signature == "CD8_Tcell"] - 37.5) < 0.05)
stopifnot(abs(alt$alt[alt$signature == "IFNg_path"] - 19.5) < 0.05)

s2 <- readLines("results/sens/s2_age_sex.txt", warn = FALSE)
s2 <- s2[grepl("prolif R2", s2)]
stg <- do.call(rbind, lapply(s2, function(ln) {
  s <- trimws(sub(" .*", "", ln))
  tail1 <- sub(".*prolif R2:", "", ln)
  v <- as.numeric(regmatches(tail1, gregexpr("[0-9]+\\.?[0-9]*", tail1))[[1]])
  vv <- tail(v, 2)
  data.frame(signature = unname(SIG_MAP[s]), nostage = vv[1], main = vv[2], stringsAsFactors = FALSE)
}))
stg$signature <- factor(stg$signature, levels = SIG_ORDER)
stg <- stg[order(stg$signature), ]
stopifnot(nrow(stg) == 7)
stopifnot(!any(is.na(stg$main)), !any(is.na(stg$nostage)))
stopifnot(abs(stg$main[stg$signature == "CAF"] - 27.6) < 0.05)
stopifnot(abs(stg$nostage[stg$signature == "CAF"] - 38.2) < 0.05)
stopifnot(abs(stg$nostage[stg$signature == "Prolifer"] - 80.3) < 0.05)

loo <- readLines("results/specB/leave_one_out_B2.txt", warn = FALSE)
loo <- loo[grepl("loo range", loo)]
lo <- do.call(rbind, lapply(loo, function(ln) {
  s <- trimws(sub(":.*", "", ln))
  rng <- as.numeric(strsplit(sub(".*loo range=", "", ln), "-")[[1]])
  data.frame(signature = unname(SIG_MAP[s]), lo = rng[1], hi = rng[2],
             full = as.numeric(sub(".*full HR=([0-9.]+).*", "\\1", ln)), stringsAsFactors = FALSE)
}))
lo$signature <- factor(lo$signature, levels = SIG_ORDER)
lo <- lo[order(lo$signature), ]
stopifnot(nrow(lo) == 7)

# ---------------- panels a and b: dumbbells ------------------------------
# Rows are ordered top to bottom in SIG_ORDER.
yrow <- rev(seq_len(7))

mk_dumbbell <- function(sig, left_v, right_v, still, xlab, ylim_v, breaks_v) {
  dq <- data.frame(signature = sig, y = yrow, l = left_v, r = right_v,
                   still = still, stringsAsFactors = FALSE)
  ggplot(dq) +
    geom_segment(aes(x = l, xend = r, y = y, yend = y),
                 linewidth = lw(LW$thin), colour = PAL$ns) +
    geom_point(aes(x = l, y = y), shape = 21, size = 2.3, fill = "white",
               colour = PAL$ink2, stroke = 0.40) +
    geom_point(aes(x = r, y = y, fill = still), shape = 21, size = 2.3,
               colour = PAL$ink, stroke = 0.28) +
    scale_fill_manual(values = c(`TRUE` = PAL$ink, `FALSE` = PAL$ring_ns), guide = "none") +
    scale_x_continuous(limits = c(0, ylim_v), breaks = breaks_v,
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0.4, 7.6), breaks = yrow, labels = sig,
                       expand = c(0, 0)) +
    labs(x = xlab, y = NULL) +
    theme_npj() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = FS$label, colour = PAL$ink, hjust = 1),
          axis.title.x = element_text(size = FS$tick),
          plot.margin = margin(8, 3, 1, 1, unit = "mm"))
}

alt$still <- alt$perm < 0.05
pA <- mk_dumbbell(alt$signature, alt$main, alt$alt, alt$still,
                  "Alternative metric: \u03c4\u00b2 explained (%)", 82, c(0, 20, 40, 60, 80)) +
  labs(subtitle = "open, main score; filled, alternative metric (dark: permutation p < 0.05)") +
  theme(plot.subtitle = element_text(size = FS$foot, colour = PAL$ink2, hjust = 0))

pB <- mk_dumbbell(stg$signature, stg$main, stg$nostage, rep(TRUE, 7),
                  "Age and sex only: \u03c4\u00b2 explained (%)", 92, c(0, 30, 60, 90)) +
  labs(subtitle = "open, adjusted for age, sex and stage; filled, age and sex only") +
  theme(plot.subtitle = element_text(size = FS$foot, colour = PAL$ink2, hjust = 0))

# ---------------- panel c: leave-one-out ---------------------------------
lo$dev <- pmax(abs(lo$lo - lo$full), abs(lo$hi - lo$full)) / lo$full * 100
lo$y <- rev(seq_len(nrow(lo)))

pC <- ggplot(lo, aes(y = y)) +
  geom_vline(xintercept = 1, linewidth = lw(LW$thin), colour = PAL$ns) +
  geom_segment(aes(x = lo, xend = hi, y = y, yend = y),
               linewidth = lw(LW$ci), colour = PAL$ink) +
  geom_point(aes(x = full), shape = 21, size = 2.3, fill = "white",
             colour = PAL$ink, stroke = 0.40) +
  geom_text(aes(x = 1.42, y = y, label = sprintf("%.1f%%", dev)),
            size = fsp(FS$foot), colour = PAL$ink2, hjust = 0) +
  scale_x_log10(breaks = c(0.9, 1.0, 1.1, 1.2, 1.3),
                expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0.4, 7.6), breaks = lo$y, labels = lo$signature,
                     expand = c(0, 0)) +
  labs(x = "Pooled HR, leave-one-out range (log scale)", y = NULL,
       subtitle = "open circle, full-sample estimate; right-hand column, largest departure from it (%)") +
  theme_npj() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = FS$label, colour = PAL$ink, hjust = 1),
        axis.title.x = element_text(size = FS$tick),
        plot.subtitle = element_text(size = FS$foot, colour = PAL$ink2, hjust = 0),
        plot.margin = margin(10, 3, 1, 1, unit = "mm"))

# ---------------- compose -------------------------------------------------
TAG_GP <- gpar(fontsize = FS$panel, fontface = "bold", fontfamily = "Arial")

draw_all <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    3, 1, heights = unit(c(H_A, GAP, H_C), "mm"))))

  pushViewport(viewport(layout.pos.row = 1))
  pushViewport(viewport(layout = grid.layout(
    1, 2, widths = unit(c(W / 2, W / 2), "mm"), heights = unit(H_A, "mm"))))
  for (i in 1:2) {
    vp <- viewport(layout.pos.row = 1, layout.pos.col = i)
    print(list(pA, pB)[[i]], vp = vp)
    grid.text(c("a", "b")[i], x = unit(1.6, "mm"),
              y = unit(1, "npc") - unit(0.6, "mm"), just = c("left", "top"),
              gp = TAG_GP, vp = vp)
  }
  popViewport(2)

  print(pC, vp = viewport(layout.pos.row = 3))
  grid.text("c", x = unit(1.6, "mm"), y = unit(1, "npc") - unit(0.6, "mm"),
            just = c("left", "top"), gp = TAG_GP, vp = viewport(layout.pos.row = 3))
  popViewport()
}

save_fig_fn(draw_all, "FigS2_robustness", W, H)

cat("FigS2 saved.\n")
print(alt[, c("signature", "main", "alt", "perm")], row.names = FALSE)
cat("max LOO departure:", sprintf("%.1f%%", max(lo$dev)), "\n")
