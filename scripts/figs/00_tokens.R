# 00_tokens.R -- single source of design tokens for all figures
# Target: npj Precision Oncology (double column 183 mm, single column 89 mm)
# All sizes given in pt and converted to ggplot units (mm) by fsp()/lw().
#
# v2 (2026-11-20): binned HR palette, unified stroke table, svglite vector output,
#                  token-level contrast assertions. Additive only -- PAL$neutral and
#                  scale_fill_hr() are untouched so the other five figures do not move.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# ---- colours -------------------------------------------------------------
PAL <- list(
  prot   = "#2166AC",   # protective (HR < 1)
  neutral = "#F7F7F7",  # HR = 1 -- used by the continuous scales, NOT by hr_bins
  risk   = "#E08214",   # risk (HR > 1)
  r2     = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
  ink    = "#1A1A1A",
  ink2   = "#4D4D4D",
  ref    = "#7F7F7F",
  ns     = "#D9D9D9",
  band   = "#F0F0F0",
  pale   = "#F2F2F2",   # neutral box fill (was hard-coded in fig1.R)
  ring_ns = "#BDBDBD"   # non-significant ring; L* = 78 vs ink L* = 10 -> greyscale safe
)

# Binned hazard-ratio palette: six steps from protective to risk.
# The two innermost steps are deliberately *tinted*, not near-white: with
# neutral = #F7F7F7 (L* 97.2) on white (L* 100) the 118 cells with
# |log2 HR| < 0.25 were invisible.  L* of the six steps: 42.5 66.3 79.0 81.5 63.2 47.7
# Adjacent CIE76 dE: 33.2 16.0 59.3 36.6 19.2 -- every step is separable.
PAL$hr_bins <- c("#2166AC", "#67A9CF", "#9ECAE1", "#FDBE85", "#E08214", "#B35806")

# Break points in log2(HR) and the labels printed on the legend.  A binned scale
# takes one label per breakpoint (n breaks -> n+1 colour steps), so these are the
# five cut-offs, not the six interval names.
# Cut-offs are 20 % / 25 % relative change in hazard, matching the manuscript's
# "close to unity" wording; they are NOT chosen from the data.
HR_BREAKS_L2 <- c(-0.58, -0.32, 0, 0.32, 1.0)   # log2 of 0.67, 0.80, 1.00, 1.25, 2.00
HR_BIN_LAB  <- c("0.67", "0.80", "1.00", "1.25", "2.00")

# Stepped R2 ramp for the Fig 3 heatmap.  The old continuous ramp started at
# #F7FBFF (L* 98.8) on white, so the five cells with R2 = 0 and the six below 5 %
# read as holes in the grid.  Six steps, monotonically darker, first one a mid
# grey so that "no variance explained" reads as an absence rather than a gap.
# L* = 86.7 83.0 74.0 61.5 48.5 32.1; adjacent CIE76 dE = 14.8 12.0 15.6 16.6 19.1
PAL$r2_bins <- c("#D9D9D9", "#B5D3E7", "#8CBCDB", "#5A9BC9", "#2878B5", "#084C8D")
R2_BREAKS   <- c(5, 15, 30, 45, 60)
R2_BIN_LAB  <- c("< 5", "5\u201315", "15\u201330", "30\u201345", "45\u201360", "\u2265 60")
# index 1..6 for a vector of R2 values; anything >= 60 falls in step 6
r2_step <- function(v) pmin(findInterval(v, R2_BREAKS) + 1L, length(PAL$r2_bins))

GERM_COL <- c(endoderm = "#0072B2", mesoderm = "#009E73", ectoderm = "#CC79A7")
GERM_SHAPE <- c(endoderm = 21, mesoderm = 24, ectoderm = 22)  # circle / triangle / square
GERM_ORDER <- c("endoderm", "mesoderm", "ectoderm")

# ---- font sizes (pt) -----------------------------------------------------
# Six steps, each with one job.  5 pt is the journal floor, never a convenience.
FS <- list(
  panel      = 8,    # bold lower-case panel letters
  axis_title = 7,    # axis titles and column-group headers
  tick       = 6,    # axis tick labels and row labels
  label      = 6,    # in-panel labels, strip text, legend text
  stat       = 5.5,  # statistics printed in panels (I2, R2, perm p)
  cell       = 5,    # numbers inside matrix cells
  foot       = 5,    # footnotes and legend body
  node       = 6.5   # schematic node titles
)
FS$legend <- FS$label
FS$legend_ttl <- FS$label
FS$ann <- FS$stat

# ---- line widths (pt) ----------------------------------------------------
LW <- list(main = 1.0, ci = 0.75, ref = 0.75, edge = 0.5, thin = 0.4)

# Disc stroke weights, one value per role (was five ad-hoc values across scripts).
STROKE <- c(disc = 0.30, ring_sig = 0.32, ring_ns = 0.22)

# Shared footnote: one symbol, one meaning, used by Fig1b and Fig2 alike.
NOTE_DEATHS <- "\u2020 ten or fewer death events"

# pt -> ggplot text size (mm)
fsp <- function(pt) pt / 2.845276
# pt -> ggplot linewidth (mm)
lw  <- function(pt) pt / 2.845276

# ---- signature display order --------------------------------------------
SIG_ORDER <- c("CD8_Tcell", "TLS", "IFNg_path", "cDC2", "CAF", "TAN", "Prolifer")
SIG_IMMUNE <- c("CD8_Tcell", "TLS", "IFNg_path", "cDC2")
SIG_STROMA <- c("CAF", "TAN", "Prolifer")
SIG_MAP <- c(CD8_cyto = "CD8_Tcell", TLS = "TLS", CAF = "CAF", TAN = "TAN",
             cDC2 = "cDC2", IFNg_path = "IFNg_path", Prolifer = "Prolifer")

# ---- theme ---------------------------------------------------------------
theme_npj <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text             = element_text(colour = PAL$ink),
      axis.line        = element_line(linewidth = lw(LW$edge), colour = PAL$ink),
      axis.ticks       = element_line(linewidth = lw(LW$edge), colour = PAL$ink),
      axis.ticks.length = unit(1.2, "mm"),
      axis.text        = element_text(size = FS$tick, colour = PAL$ink),
      axis.title       = element_text(size = FS$axis_title, colour = PAL$ink),
      plot.title       = element_blank(),
      legend.title     = element_text(size = FS$legend_ttl, colour = PAL$ink),
      legend.text      = element_text(size = FS$legend, colour = PAL$ink),
      legend.key.size  = unit(3.2, "mm"),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(1, 1, 1, 1, unit = "mm"),
      strip.background = element_blank(),
      strip.text       = element_text(size = FS$label, colour = PAL$ink)
    )
}

# Panel letters are placed by patchwork (plot_annotation(tag_levels = "a")),
# never by hand-computed npc constants.
TAG <- list(size = FS$panel, face = "bold", colour = PAL$ink)
theme_tag <- theme(plot.tag = element_text(size = TAG$size, face = TAG$face,
                                           colour = TAG$colour),
                   plot.tag.position = c(0, 1))

# ---- scales --------------------------------------------------------------
scale_fill_hr <- function(...) {
  scale_fill_gradient2(
    low = PAL$prot, mid = PAL$neutral, high = PAL$risk, midpoint = 0,
    limits = c(-1, 1), oob = scales::squish,
    guide = guide_colourbar(title = "HR", barwidth = unit(26, "mm"),
                            barheight = unit(2.6, "mm"), frame.colour = NA),
    ...
  )
}

# Binned version used by Fig2.  oob = squish keeps the extreme cells inside the
# top step; the ten cells beyond |log2 HR| > 1 are printed with their exact value,
# so no overflow glyph is needed.  Do NOT move the midpoint or use asinh/log-modulus
# here: that would claim HR = 0.9 is "more neutral" than HR = 1.0.
scale_fill_hr_bin <- function(..., guide = "none") {
  scale_fill_stepsn(
    colours = PAL$hr_bins,
    breaks  = HR_BREAKS_L2,
    labels  = HR_BIN_LAB,
    limits  = c(-1.02, 1.02),
    oob     = scales::squish,
    guide   = guide,
    ...
  )
}

# ---- export helpers ------------------------------------------------------
# Vector PDF (cairo, fonts embedded), vector SVG (svglite: live <text>, not outlined
# glyph paths), 600 dpi PNG.
save_fig <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  dir.create("figures", showWarnings = FALSE)
  w <- width_mm / 25.4; h <- height_mm / 25.4

  grDevices::cairo_pdf(file.path("figures", paste0(stem, ".pdf")),
                       width = w, height = h, family = "Arial", bg = "white")
  print(plot); grDevices::dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file.path("figures", paste0(stem, ".svg")),
                     width = w, height = h, bg = "white")
  } else {
    grDevices::svg(file.path("figures", paste0(stem, ".svg")), width = w, height = h,
                   family = "Arial", bg = "white", onefile = TRUE)
  }
  print(plot); grDevices::dev.off()

  grDevices::png(file.path("figures", paste0(stem, ".png")),
                 width = round(w * dpi), height = round(h * dpi), res = dpi,
                 type = "cairo", family = "Arial", bg = "white", units = "px")
  print(plot); grDevices::dev.off()

  invisible(stem)
}

# Same three-format export for figures built by a drawing function (grid layouts
# that print several viewports).  Keeps every figure on one SVG device, so all
# six deliverables carry live <text> instead of outlined glyph paths.
save_fig_fn <- function(draw, stem, width_mm, height_mm, dpi = 600) {
  dir.create("figures", showWarnings = FALSE)
  w <- width_mm / 25.4; h <- height_mm / 25.4

  grDevices::cairo_pdf(file.path("figures", paste0(stem, ".pdf")),
                       width = w, height = h, family = "Arial", bg = "white")
  draw(); grDevices::dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file.path("figures", paste0(stem, ".svg")),
                     width = w, height = h, bg = "white")
  } else {
    grDevices::svg(file.path("figures", paste0(stem, ".svg")), width = w, height = h,
                   family = "Arial", bg = "white", onefile = TRUE)
  }
  draw(); grDevices::dev.off()

  grDevices::png(file.path("figures", paste0(stem, ".png")),
                 width = round(w * dpi), height = round(h * dpi), res = dpi,
                 type = "cairo", family = "Arial", bg = "white", units = "px")
  draw(); grDevices::dev.off()

  invisible(stem)
}

# ---- small utilities -----------------------------------------------------
fmt_hr <- function(hr, lo, hi, dig = 2) {
  sprintf(paste0("%.", dig, "f (%.", dig, "f\u2013%.", dig, "f)"), hr, lo, hi)
}

# ---- token-level contrast guard -----------------------------------------
# Catch a drifting palette at the source rather than after rendering.
assert_tokens <- function() {
  if (!requireNamespace("farver", quietly = TRUE)) return(invisible(NULL))
  de76 <- function(a, b) {
    la <- farver::decode_colour(a, to = "lab")
    lb <- farver::decode_colour(b, to = "lab")
    sqrt(sum((la[1, ] - lb[1, ])^2))
  }
  # adjacent hazard bins must stay apart
  gaps <- vapply(seq_len(length(PAL$hr_bins) - 1),
                 function(i) de76(PAL$hr_bins[i], PAL$hr_bins[i + 1]), numeric(1))
  if (any(gaps < 6)) stop("hr_bins: adjacent dE76 < 6 -- ", paste(round(gaps, 1), collapse = " "))
  # germ layers must stay apart
  gp <- c(de76(GERM_COL[1], GERM_COL[2]), de76(GERM_COL[1], GERM_COL[3]),
          de76(GERM_COL[2], GERM_COL[3]))
  if (any(gp < 12)) stop("GERM_COL: pairwise dE76 < 12 -- ", paste(round(gp, 1), collapse = " "))
  # R2 steps must be separable too, including grey step 1 vs the palest blue
  r2g <- vapply(seq_len(length(PAL$r2_bins) - 1),
                function(i) de76(PAL$r2_bins[i], PAL$r2_bins[i + 1]), numeric(1))
  if (any(r2g < 8)) stop("r2_bins: adjacent dE76 < 8 -- ", paste(round(r2g, 1), collapse = " "))
  # significant vs non-significant ring must survive greyscale
  ldiff <- abs(farver::decode_colour(PAL$ink, to = "lab")[1, 1] -
               farver::decode_colour(PAL$ring_ns, to = "lab")[1, 1])
  if (ldiff < 40) stop("ring contrast: L* difference only ", round(ldiff, 1))
  invisible(list(bin_gaps = gaps, germ = gp, ring_L = ldiff))
}

assert_tokens()
