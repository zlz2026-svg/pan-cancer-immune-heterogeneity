# 00_legend.R -- millimetre-coordinate primitives for hand-drawn legends.
#
# Why this file exists: patchwork's guide collector is not usable here.  With
# ggplot2 4.0.4, guides = "collect" pushes every legend into one full-height
# right-hand column (shrinking the panels), guides = "keep" is ignored, and
# theme_void() + legend.position = "center" silently drops the legend entirely.
# Fig 2 already solved this by drawing its legend into a plain mm canvas, and
# these helpers generalise that solution.
#
# Convention: callers pass millimetre coordinates directly.  A panel drawn into a
# viewport of exactly W x H mm with plot.margin = 0 has 1 data unit = 1 mm.
# Every helper returns a bare data.frame; the caller adds the matching geom, so a
# legend is assembled with the same grammar as any other panel.

# ---- canvas --------------------------------------------------------------
mm_canvas <- function(W, H) {
  ggplot() +
    coord_cartesian(xlim = c(0, W), ylim = c(0, H), expand = FALSE, clip = "off") +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = "white", colour = NA))
}

# ---- shapes --------------------------------------------------------------
# Circle drawn as a polygon so the diameter is exact in millimetres.
mm_disc <- function(cx, cy, d, id = 1L, n = 72) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + d / 2 * cos(th), y = cy + d / 2 * sin(th), id = id)
}

# Filled bar from (x0, y0) to (x1, y1).
mm_rect <- function(x0, y0, x1, y1) {
  data.frame(xmin = x0, xmax = x1, ymin = y0, ymax = y1)
}

# Regular triangle / square centred on (cx, cy), sized so that the *bounding box*
# equals a disc of diameter d.  ggplot2 draws shapes 24 and 22 at the same bounding
# size as shape 21, so the legend has to match that, not match ink area -- an
# area-parity triangle would come out visibly smaller than the points it explains.
mm_tri <- function(cx, cy, d, id = 1L) {
  R <- d / 2
  a <- seq(pi / 2, pi / 2 + 2 * pi / 3, length.out = 3)
  data.frame(x = cx + R * cos(a), y = cy + R * sin(a), id = id)
}
mm_sq <- function(cx, cy, d, id = 1L) {
  h <- d / 2 * 0.90
  data.frame(x = cx + c(-h, h, h, -h), y = cy + c(-h, -h, h, h), id = id)
}

# Swatch dispatcher: one call per legend glyph, keyed on the germ-layer shapes so
# the legend and the panels can never drift apart.
mm_swatch <- function(cx, cy, d, shape = 21, id = 1L) {
  if (shape == 21) return(mm_disc(cx, cy, d, id))
  if (shape == 24) return(mm_tri(cx, cy, d, id))
  if (shape == 22) return(mm_sq(cx, cy, d, id))
  stop("mm_swatch: unsupported shape ", shape)
}

# ---- text ---------------------------------------------------------------
# Sizes are pt, converted exactly as geom_text wants them (mm units).
mm_label <- function(x, y, label, size = FS$legend, colour = PAL$ink,
                     hjust = 0, vjust = 0.5, fontface = "plain") {
  data.frame(x = x, y = y, label = label, colour = colour,
             hjust = hjust, vjust = vjust, fontface = fontface,
             stringsAsFactors = FALSE)
}

mm_text_layers <- function(df) {
  list(geom_text(data = df, aes(x = x, y = y, label = label, colour = colour,
                                hjust = hjust, vjust = vjust, fontface = fontface),
                 size = fsp(FS$legend)),
       scale_colour_identity())
}

# ---- composite keys -----------------------------------------------------
# A horizontal stepped colour key: n blocks of block_w mm, each labelled inside.
# Labels go inside the blocks because a 6-block R2 key with labels underneath
# needs 18 mm of height, which no bottom legend strip can spare.
mm_colour_key <- function(x0, y, block_w, block_h, cols, labs, col_lab) {
  n <- length(cols)
  rects <- do.call(rbind, lapply(seq_len(n), function(i) {
    mm_rect(x0 + (i - 1) * block_w, y - block_h / 2,
            x0 + i * block_w, y + block_h / 2)
  }))
  rects$fill <- cols
  # white text on the three darkest steps (L* < 75), ink on the pale ones
  lum <- if (requireNamespace("farver", quietly = TRUE))
    farver::decode_colour(cols, to = "lab")[, 1] else rep(90, n)
  labs_df <- data.frame(x = x0 + (seq_len(n) - 0.5) * block_w, y = y,
                        label = labs, colour = ifelse(lum < 75, "white", col_lab),
                        stringsAsFactors = FALSE)
  list(rects = rects, labs = labs_df)
}
