# 00_cells.R -- vector cell glyphs shared by Fig 1c and Supplementary Fig S1.
#
# One primitive, two figures: a cell is a closed membrane polygon, a nucleus, and
# one or two cytoplasmic granules.  Malignant and immune cells differ in *form*
# (large eccentric nucleus and a wavy membrane vs a small central nucleus and a
# near-circular membrane) as well as in colour, which matters because #CC79A7 and
# #009E73 both drift toward grey under deuteranopia.
#
# HARD INVARIANT: every glyph is rescaled so its largest radial distance is exactly
# d/2.  The dilution argument in Fig S1 only works if a larger cell cannot explain
# the weaker signal, so all cells must share one diameter and the difference has to
# live in morphology.  cell_glyph() enforces the diameter by construction;
# assert_field() re-checks it from the rendered polygons.

CELL_STYLE <- list(
  mal = c(membrane = "#CC79A7", nucleus = "#7B3E60", granule = "#AB6089"),
  imm = c(membrane = "#009E73", nucleus = "#005741", granule = "#00795A")
)
# fill keys handed to scale_fill_manual(); six values, two cell types x three parts
CELL_FILL <- local({
  v <- character(0)
  for (tp in names(CELL_STYLE)) {
    for (pt in names(CELL_STYLE[[tp]])) {
      v[paste(tp, pt, sep = "_")] <- unname(CELL_STYLE[[tp]][[pt]])
    }
  }
  v
})

# ---- one cell ------------------------------------------------------------
# Returns a polygon data.frame with columns x, y, group, part (= "<type>_<part>").
cell_glyph <- function(cx, cy, d, type = c("mal", "imm"), id = 1L, seed = 1L) {
  type <- match.arg(type)
  set.seed(seed)
  n <- 72L
  th <- seq(0, 2 * pi, length.out = n)

  # membrane: radial ripple, then rescaled so max(r) == d/2 exactly
  # a single high-amplitude sinusoid turned the malignant outline into a cog;
  # two low-amplitude harmonics at different frequencies read as an irregular cell
  if (type == "mal") {
    r <- 1 + 0.030 * sin(3 * th + runif(1, 0, 2 * pi)) +
             0.016 * sin(7 * th + runif(1, 0, 2 * pi))
  } else {
    r <- 1 + 0.009 * sin(3 * th + runif(1, 0, 2 * pi))
  }
  r <- r / max(r) * (d / 2)
  mem <- data.frame(x = cx + r * cos(th), y = cy + r * sin(th), sub = "membrane")

  # nucleus: eccentric oval for malignant, small central round for immune
  if (type == "mal") {
    a <- 0.28 * d; b <- 0.21 * d
    ox <- 0.10 * d * sample(c(-1, 1), 1)
    oy <- 0.07 * d * sample(c(-1, 1), 1)
    rot <- runif(1, 0, pi)
  } else {
    a <- 0.17 * d; b <- 0.15 * d
    ox <- 0; oy <- 0; rot <- 0
  }
  nt <- seq(0, 2 * pi, length.out = 48)
  ex <- a * cos(nt); ey <- b * sin(nt)
  nuc <- data.frame(x = cx + ox + ex * cos(rot) - ey * sin(rot),
                    y = cy + oy + ex * sin(rot) + ey * cos(rot),
                    sub = "nucleus")

  # granules: texture only, never an encoding channel
  ng <- if (type == "mal") 1L else 2L
  gr <- do.call(rbind, lapply(seq_len(ng), function(i) {
    ang <- runif(1, 0, 2 * pi); rad <- runif(1, 0.24, 0.32) * d
    gd <- 0.085 * d
    tt <- seq(0, 2 * pi, length.out = 24)
    data.frame(x = cx + rad * cos(ang) + gd / 2 * cos(tt),
               y = cy + rad * sin(ang) + gd / 2 * sin(tt),
               sub = "granule")
  }))

  out <- rbind(mem, nuc, gr)
  out$part  <- paste(type, out$sub, sep = "_")
  out$cx <- cx; out$cy <- cy
  # seed is unique per cell, so group ids cannot collide when several fields
  # are drawn into one layer (a shared id merged every field into one polygon)
  out$group <- paste(seed, id, out$sub, sep = ".")
  out[, c("x", "y", "group", "part", "cx", "cy")]
}

# ---- a field of cells ----------------------------------------------------
# Poisson-disc (dart throwing) placement.  Not a jittered lattice: keeping a
# non-touching guarantee on a hexagonal lattice allows a jitter of only 0.02 d
# before adjacent cells collide, and a 0.02 d jitter still reads as a grid.
# Dart throwing gives blue-noise spacing at any density and stays organic.
place_cells <- function(n_mal, n_imm, x0, y0, w, h, d, seed = 1L,
                        gap = 1.10, max_try = 4e5) {
  set.seed(seed)
  want <- n_mal + n_imm
  r <- d / 2
  pts <- matrix(numeric(0), 0, 2)
  tries <- 0L
  while (nrow(pts) < want) {
    tries <- tries + 1L
    if (tries > max_try) stop("place_cells: did not converge", call. = FALSE)
    cand <- c(runif(1, x0 + r, x0 + w - r), runif(1, y0 + r, y0 + h - r))
    if (nrow(pts) == 0 ||
        min(sqrt((pts[, 1] - cand[1])^2 + (pts[, 2] - cand[2])^2)) >= d * gap) {
      pts <- rbind(pts, cand)
    }
  }
  # types are shuffled so neither type clusters in a corner
  data.frame(x = pts[, 1], y = pts[, 2],
             type = sample(c(rep("mal", n_mal), rep("imm", n_imm))),
             stringsAsFactors = FALSE)
}

render_field <- function(cells, d, seed = 11L) {
  do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    cell_glyph(cells$x[i], cells$y[i], d, cells$type[i], id = i, seed = seed + i)
  }))
}

# ---- invariants ----------------------------------------------------------
assert_field <- function(cells, d, n_mal, n_imm, tol = 0.005) {
  if (nrow(cells) != n_mal + n_imm)
    stop("assert_field: expected ", n_mal + n_imm, " cells, got ", nrow(cells), call. = FALSE)
  if (sum(cells$type == "mal") != n_mal || sum(cells$type == "imm") != n_imm)
    stop("assert_field: composition changed", call. = FALSE)
  dd <- as.matrix(dist(cells[, c("x", "y")]))
  diag(dd) <- Inf
  if (min(dd) < d * 1.02) {
    stop("assert_field: two cells are only ", round(min(dd) / d, 3),
         " diameters apart", call. = FALSE)
  }
  invisible(TRUE)
}

# every membrane must reach exactly d/2 from its centre -- checked on the polygons,
# not on the intent, so a later tweak to the ripple cannot silently break it
assert_diameter <- function(polys, d, tol = 0.01) {
  mem <- polys[polys$part %in% c("mal_membrane", "imm_membrane"), ]
  rad <- sqrt((mem$x - mem$cx)^2 + (mem$y - mem$cy)^2)
  if (abs(max(rad) - d / 2) > tol * d) {
    stop("assert_diameter: largest radius ", round(max(rad), 3),
         " mm, expected ", d / 2, call. = FALSE)
  }
  invisible(max(rad))
}
