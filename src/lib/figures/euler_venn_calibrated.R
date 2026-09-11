# Calibrated Euler/Venn diagram (Figure 3d): ClinVar Benign/Pathogenic/
# Uncertain-significance ("VUS") variants overlapping with gnomAD population
# frequency data, plus a separate non-overlapping "Unobserved variants"
# circle for SNVs never seen in ClinVar or gnomAD.
#
# Renders at the *exact* final placed size (cairo_pdf, mm units, Arial),
# matching the convention in variant_bar_calibrated.R/
# bubble_scatter_calibrated.R -- solved backward from a fixed *total* canvas
# size (legend included, like the bar chart's "total figure" spec) rather
# than bubble_scatter_calibrated.R's fixed-plot-area/grow-outward approach.
#
# The circle *geometry* (Euler-diagram optimization -- finding circle radii/
# positions that best reproduce the requested subset sizes) is still done by
# the eulerr package, exactly as the original uncalibrated notebook did --
# that numerical fit isn't worth reimplementing, matching this session's
# convention of keeping an established algorithm (ggrepel's force-directed
# label placement in bubble_scatter_calibrated.R) while replacing only the
# *rendering* with our own mm-calibrated layout. This file takes over
# rendering because eulerr's own plot() isn't unit-calibrated and can't
# rotate/shrink individual labels to avoid collisions the way the published
# figure's labels are (each region's count is placed by hand here instead of
# eulerr's own labeling heuristic).
#
# Expects the caller to have already loaded: eulerr, dplyr, ggplot2,
# ggforce, extrafont (+ loadfonts(device = 'all')), and grid.

point_in_mm <- 0.3527778
FONT_FAMILY <- "Arial"
LABEL_PT <- 7
LABEL_PT_MIN <- 4
GAP_MM <- 1.2
LEGEND_LINEHEIGHT <- 0.8

# Matches the original (uncalibrated) notebook's fill colors/alpha exactly,
# so the overlap regions blend to the same colors the published figure used.
SET_COLORS <- c(B = "#1D7AAB", P = "#CA7682", G = "#A5A9BE", V = "#545454", D = "white")
FILL_ALPHA <- 0.6
SET_ORDER <- c("B", "P", "G", "V", "D")
# Draw G (which every other set overlaps) first so B/P/V paint on top of
# it -- alpha blending is order-dependent, and this matches the reference
# figure's overlap colors.
DRAW_ORDER <- c("G", "B", "P", "V", "D")
SET_DISPLAY <- c(
  B = "Benign", P = "Pathogenic", G = "gnomAD",
  V = "Uncertain\nsignificance", D = "Unobserved variants"
)
# Legend layout: a 2x2 grid (col1: B/P, col2: G/V) plus D as a third,
# vertically-centered column -- matches the reference figure.
LEGEND_GRID <- list(
  col1 = c("B", "P"),
  col2 = c("G", "V"),
  col3 = c("D")
)

text_width_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY, fontface = "plain") {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = fontface)
  vapply(label, function(l) convertWidth(grobWidth(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

text_height_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY, fontface = "plain", lineheight = 1.2) {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = fontface, lineheight = lineheight)
  vapply(label, function(l) convertHeight(grobHeight(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

# Widest single line of a (possibly multi-line, "\n"-joined) label -- used
# to size legend columns, since grobWidth() on a multi-line string already
# returns its widest line, but we sometimes need that per-candidate-string
# rather than per-grob.
widest_line_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY) {
  lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
  max(text_width_mm(lines, pt, family))
}

normalize <- function(v) v / sqrt(sum(v^2))

# Try LABEL_PT, then shrink (in whole points) until `label` fits within
# max_width_mm x max_height_mm, down to LABEL_PT_MIN.
fit_label_size <- function(label, max_width_mm, max_height_mm) {
  for (pt in seq(LABEL_PT, LABEL_PT_MIN, by = -0.5)) {
    w <- widest_line_mm(label, pt)
    h <- text_height_mm(label, pt)
    if (w <= max_width_mm && h <= max_height_mm) {
      return(pt)
    }
  }
  LABEL_PT_MIN
}

# gene_counts: a list with B, P, G, V (exclusive-of-gnomAD counts), D
# (Unreported_SNVs), and BG/PG/VG (each set's overlap with gnomAD) --
# already-collapsed integers, matching the euler() call the original
# notebook made.
plot_euler_venn_calibrated <- function(B, P, G, V, D, BG, PG, VG,
                                        total_width_mm = 3 * 25.4, total_height_mm = 2.4 * 25.4,
                                        layout_seed = 1) {
  measure_dev_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(measure_dev_file)
  on.exit({
    grDevices::dev.off()
    unlink(measure_dev_file)
  }, add = TRUE)

  counts <- c(B = B, P = P, G = G, V = V, D = D)

  # eulerr::euler()'s optimizer starts from a random initial layout and
  # isn't seeded internally, so calling it twice with *identical* input
  # produces two different (though equally valid -- same areas/overlaps)
  # circle arrangements. Fix a seed so the diagram's layout is reproducible
  # across reruns, restoring the caller's own RNG state afterward so this
  # doesn't perturb unrelated randomness elsewhere in the pipeline.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    },
    add = TRUE
  )
  set.seed(layout_seed)

  fit <- eulerr::euler(c(B = B, P = P, G = G, D = D, V = V, "B&G" = BG, "P&G" = PG, "V&G" = VG))
  ell <- fit$ellipses # h, k, a, b, phi -- a == b here (circles, not ellipses)

  circles <- lapply(SET_ORDER, function(s) list(h = ell[s, "h"], k = ell[s, "k"], r = ell[s, "a"]))
  names(circles) <- SET_ORDER

  # ---- Legend sizing (text measurements only; doesn't depend on the
  # diagram's own scale) ----
  swatch_mm <- text_height_mm("Mg", LABEL_PT) * 1.3
  legend_row_gap_mm <- GAP_MM * 0.6
  legend_col_gap_mm <- GAP_MM * 1.8

  col_width <- function(members) {
    labels <- SET_DISPLAY[members]
    swatch_mm + GAP_MM / 2 + max(vapply(labels, widest_line_mm, numeric(1), pt = LABEL_PT))
  }
  col1_w <- col_width(LEGEND_GRID$col1)
  col2_w <- col_width(LEGEND_GRID$col2)
  col3_w <- col_width(LEGEND_GRID$col3)
  legend_w <- col1_w + legend_col_gap_mm + col2_w + legend_col_gap_mm + col3_w

  row_height <- function(member) {
    max(swatch_mm, text_height_mm(SET_DISPLAY[[member]], LABEL_PT, lineheight = LEGEND_LINEHEIGHT))
  }
  row1_h <- max(row_height(LEGEND_GRID$col1[1]), row_height(LEGEND_GRID$col2[1]))
  row2_h <- max(row_height(LEGEND_GRID$col1[2]), row_height(LEGEND_GRID$col2[2]))
  legend_h <- row1_h + legend_row_gap_mm + row2_h

  # ---- Diagram scale: fit the 5 circles' bounding box into the remaining
  # canvas above the legend, preserving circle shape (one uniform scale for
  # both axes) ----
  top_pad_mm <- 1.5
  side_pad_mm <- 1.5
  diagram_to_legend_gap_mm <- GAP_MM * 1.5
  # text_height_mm()'s grobHeight-based measurement doesn't reliably reserve
  # room for descenders (e.g. the "g" in "Pathogenic"/"significance"), and
  # the legend's bottom row otherwise sits with its measured box flush
  # against y = 0 -- leaving zero tolerance and clipping any descender that
  # extends slightly further than measured. Reserve a small margin below it.
  bottom_pad_mm <- 1

  bbox <- list(
    xmin = min(vapply(circles, function(c) c$h - c$r, numeric(1))),
    xmax = max(vapply(circles, function(c) c$h + c$r, numeric(1))),
    ymin = min(vapply(circles, function(c) c$k - c$r, numeric(1))),
    ymax = max(vapply(circles, function(c) c$k + c$r, numeric(1)))
  )
  bbox_w <- bbox$xmax - bbox$xmin
  bbox_h <- bbox$ymax - bbox$ymin

  width_budget_mm <- total_width_mm - 2 * side_pad_mm
  height_budget_mm <- total_height_mm - top_pad_mm - diagram_to_legend_gap_mm - legend_h - bottom_pad_mm
  scale <- min(width_budget_mm / bbox_w, height_budget_mm / bbox_h)

  diagram_w_mm <- bbox_w * scale
  diagram_h_mm <- bbox_h * scale
  diagram_top_y_mm <- total_height_mm - top_pad_mm
  x_center_offset_mm <- (total_width_mm - diagram_w_mm) / 2
  # A uniform scale + translate (both axes scaled identically, so circles
  # stay circles): x grows rightward from the bbox's left edge; y grows
  # downward from the diagram's top edge, matching eulerr's own coordinate
  # convention, so no reflection is needed. Every downstream computation
  # (label offsets, angles) is re-derived from the *transformed* circle
  # centers, so this mapping's orientation doesn't need to match any
  # particular external convention -- it only has to be applied uniformly.
  to_mm <- function(pt) {
    list(
      x = (pt$h - bbox$xmin) * scale + x_center_offset_mm,
      y = diagram_top_y_mm - (bbox$ymax - pt$k) * scale
    )
  }
  circles_mm <- lapply(circles, function(c) {
    center <- to_mm(list(h = c$h, k = c$k))
    list(x0 = center$x, y0 = center$y, r = c$r * scale)
  })

  # ---- Label placement ----
  # Normalize an angle (degrees) to (-90, 90] so rotated text never renders
  # upside down.
  normalize_angle <- function(angle_deg) {
    if (angle_deg > 90) angle_deg <- angle_deg - 180
    if (angle_deg < -90) angle_deg <- angle_deg + 180
    angle_deg
  }

  # Overlap (X & G) labels sit at the midpoint of the lens's extent *along
  # the line connecting the two centers* (not the "radical axis" -- that's
  # the point equidistant in circle-power from both circles, which isn't
  # the same as the visual center of the lens when the two circles differ
  # a lot in size), rotated to align with the lens's long axis
  # (perpendicular to the center-to-center line, since the two circles'
  # intersection points are symmetric *across* that line).
  overlap_label <- function(set_name, overlap_count) {
    cx <- circles_mm[[set_name]]
    cg <- circles_mm$G
    d <- sqrt((cg$x0 - cx$x0)^2 + (cg$y0 - cx$y0)^2)
    dir_to_g <- c(cg$x0 - cx$x0, cg$y0 - cx$y0) / d
    # The lens, sliced along the center-to-center line and parametrized by
    # distance t from cx's own center (toward G), spans from where it
    # enters G (t = d - cg$r, or 0 if cx's center is already inside G) to
    # where it exits X (t = cx$r).
    lo <- max(0, d - cg$r)
    hi <- min(cx$r, d)
    t <- (lo + hi) / 2
    point <- c(cx$x0, cx$y0) + dir_to_g * t
    # The lens's half-width perpendicular to the center line at this point
    # is limited by whichever circle is narrower there.
    half_chord <- min(sqrt(max(cx$r^2 - t^2, 0)), sqrt(max(cg$r^2 - (d - t)^2, 0)))
    angle_deg <- normalize_angle(atan2(dir_to_g[2], dir_to_g[1]) * 180 / pi + 90)
    label <- format(overlap_count, big.mark = ",", trim = TRUE)
    pt <- fit_label_size(label, half_chord * 1.9, half_chord * 1.9)
    list(x = point[1], y = point[2], label = label, angle = angle_deg, pt = pt)
  }

  # Exclusive-region (X-only, X != G) labels sit inside X's own crescent,
  # offset from its center directly *away* from G along the line connecting
  # them (the simplest point guaranteed to avoid the overlap lens for a
  # circle that only overlaps one other set) -- and, since that's usually
  # the crescent's longest straight run, *tilted* to read along that same
  # direction so a longer/larger label can still fit than a horizontal one
  # would.
  exclusive_label <- function(set_name, count, away_from = "G") {
    cx <- circles_mm[[set_name]]
    away <- circles_mm[[away_from]]
    d <- sqrt((away$x0 - cx$x0)^2 + (away$y0 - cx$y0)^2)
    dir_away <- normalize(c(cx$x0 - away$x0, cx$y0 - away$y0))
    label <- format(count, big.mark = ",", trim = TRUE)
    angle_deg <- normalize_angle(atan2(dir_away[2], dir_away[1]) * 180 / pi)

    # Along the ray from cx's center directly away from `away_from`: the
    # point stays inside cx up to t = cx$r (its own far edge), and is
    # outside `away_from` once t >= away$r - d (0 if cx's center is already
    # outside it). For a small, heavily-overlapped circle the crescent is
    # actually *wider tangentially than radially* (it's a thin sliver
    # peeking out from under the bigger circle, curving sideways more than
    # it extends outward) -- so try the label both radially-tilted (reading
    # away from `away_from`) and tangentially-tilted (perpendicular to
    # that, i.e. angle_deg + 90), at a handful of positions along the ray,
    # and keep whichever combination lets it render largest.
    t_min <- max(0, away$r - d)
    t_max <- cx$r
    best <- NULL
    for (frac in c(0.3, 0.4, 0.5, 0.6, 0.7)) {
      t <- t_min + frac * (t_max - t_min)
      along <- 2 * min(t - t_min, t_max - t) * 0.9
      perp <- 2 * sqrt(max(cx$r^2 - t^2, 0)) * 0.85
      candidates <- list(
        list(angle = angle_deg, pt = fit_label_size(label, along, perp)),
        list(angle = normalize_angle(angle_deg + 90), pt = fit_label_size(label, perp, along))
      )
      for (cand in candidates) {
        if (is.null(best) || cand$pt > best$pt) {
          best <- list(t = t, angle = cand$angle, pt = cand$pt)
        }
      }
    }
    point <- c(cx$x0, cx$y0) + dir_away * best$t
    list(x = point[1], y = point[2], label = label, angle = best$angle, pt = best$pt)
  }

  center_label <- function(set_name, count, max_dim) {
    cx <- circles_mm[[set_name]]
    label <- format(count, big.mark = ",", trim = TRUE)
    pt <- fit_label_size(label, max_dim, max_dim)
    list(x = cx$x0, y = cx$y0, label = label, angle = 0, pt = pt)
  }

  labels <- list(
    exclusive_label("B", B),
    exclusive_label("P", P),
    exclusive_label("V", V),
    center_label("G", G, circles_mm$G$r * 0.9),
    center_label("D", D, circles_mm$D$r * 1.2),
    overlap_label("B", BG),
    overlap_label("P", PG),
    overlap_label("V", VG)
  )
  labels_df <- do.call(rbind, lapply(labels, as_tibble))

  # ---- Legend geometry (sits at the bottom of the canvas, y in [0, legend_h]) ----
  legend_left_x <- (total_width_mm - legend_w) / 2
  row1_y <- bottom_pad_mm + legend_h - row1_h / 2
  row2_y <- row1_y - row1_h / 2 - legend_row_gap_mm - row2_h / 2

  legend_col_x <- c(
    legend_left_x,
    legend_left_x + col1_w + legend_col_gap_mm,
    legend_left_x + col1_w + legend_col_gap_mm + col2_w + legend_col_gap_mm
  )

  legend_entries <- do.call(rbind, lapply(seq_along(SET_ORDER), function(i) {
    s <- SET_ORDER[i]
    if (s %in% LEGEND_GRID$col1) {
      col_x <- legend_col_x[1]
      y <- if (identical(LEGEND_GRID$col1[1], s)) row1_y else row2_y
    } else if (s %in% LEGEND_GRID$col2) {
      col_x <- legend_col_x[2]
      y <- if (identical(LEGEND_GRID$col2[1], s)) row1_y else row2_y
    } else {
      col_x <- legend_col_x[3]
      y <- (row1_y + row2_y) / 2
    }
    tibble(
      set = s, swatch_x = col_x, label_x = col_x + swatch_mm + GAP_MM / 2, y = y,
      label = SET_DISPLAY[[s]]
    )
  }))

  x_range <- c(0, total_width_mm)
  y_range <- c(0, total_height_mm)

  p <- ggplot()
  for (s in DRAW_ORDER) {
    cm <- circles_mm[[s]]
    fill <- if (s == "D") NA else scales::alpha(SET_COLORS[[s]], FILL_ALPHA)
    p <- p + ggforce::geom_circle(
      data = tibble(x0 = cm$x0, y0 = cm$y0, r = cm$r),
      aes(x0 = x0, y0 = y0, r = r), inherit.aes = FALSE,
      fill = fill, colour = "black", linewidth = 0.5 * point_in_mm
    )
  }
  p <- p +
    geom_text(
      data = labels_df, inherit.aes = FALSE,
      aes(x = x, y = y, label = label, angle = angle, size = pt * point_in_mm),
      family = FONT_FAMILY, colour = "black"
    ) +
    scale_size_identity() +
    # ---- Legend ----
    geom_tile(
      data = legend_entries, inherit.aes = FALSE,
      aes(x = swatch_x + swatch_mm / 2, y = y, fill = I(scales::alpha(SET_COLORS[set], FILL_ALPHA))),
      width = swatch_mm, height = swatch_mm, colour = "black", linewidth = 0.25
    ) +
    geom_text(
      data = legend_entries, inherit.aes = FALSE,
      aes(x = label_x, y = y, label = label),
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 0, vjust = 0.5,
      lineheight = LEGEND_LINEHEIGHT
    ) +
    coord_cartesian(xlim = x_range, ylim = y_range, expand = FALSE, clip = "off") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

  list(plot = p, width_mm = diff(x_range), height_mm = diff(y_range))
}
