# Shared calibrated stacked bar chart logic (Figure 3, "Total Variant effect
# measurements" per gene, stacked by assay type -- SGE/VAMP-seq/Other).
#
# Renders at the *exact* final placed size (cairo_pdf, mm units, Arial),
# matching the convention in bubble_scatter_calibrated.R/
# metrics_bar_calibrated.R -- axis titles/ticks and the legend are laid out
# directly in mm, solved backward from a *fixed total canvas size* (unlike
# bubble_scatter_calibrated.R's fixed-plot-area/grow-outward approach) since
# this chart's total dimensions, not just its plot area, were specified.
#
# Expects the caller to have already loaded: dplyr, ggplot2, extrafont (+
# loadfonts(device = 'all')), and grid.

point_in_mm <- 0.3527778
FONT_FAMILY <- "Arial"
LABEL_PT <- 7

GAP_MM <- 1.2
TICK_LEN_MM <- 0.8
AXIS_LWD <- 0.5 * point_in_mm
STRIPE_WIDTH_MM <- 0.6

ASSAY_COLORS <- c(Other = "grey", SGE = "#E57525", VAMP = "#09A3A3")
# Bottom-to-top stack order comes from alphabetically sorting the fill factor
# -- ggplot2's own default behavior for an unordered factor, confirmed
# against the original (uncalibrated) chart's actual rendering, matching
# this stack order exactly. A variant tested by more than one assay type for
# the same gene (see build_figure3c_variant_assay_combos in
# src/build_figure3_data.py -- e.g. LDLR has 16,413 variants tested by both
# an "Other" and a Vamp-seq dataset) gets its own "+"-joined combo category,
# stacked directly next to its own base categories and rendered as a
# diagonal two-color stripe rather than silently folded into one of them.
ASSAY_STACK_ORDER <- c("VAMP", "Other+VAMP", "SGE", "Other+SGE", "Other")
# Left-to-right legend order/labels match the original (preprint) figure's
# in-plot legend for the base categories; combo categories (only rendered
# when actually present in the data) follow as a second row.
ASSAY_LEGEND_ORDER <- c("VAMP", "SGE", "Other")
ASSAY_LEGEND_DISPLAY <- c(Other = "Other", SGE = "SGE", VAMP = "VAMP-seq")
ASSAY_COMBO_ORDER <- c("Other+VAMP", "Other+SGE")
ASSAY_COMBO_DISPLAY <- c("Other+VAMP" = "Other + VAMP-seq", "Other+SGE" = "Other + SGE")

text_width_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY, fontface = "plain") {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = fontface)
  vapply(label, function(l) convertWidth(grobWidth(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

text_height_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY, fontface = "plain") {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = fontface)
  vapply(label, function(l) convertHeight(grobHeight(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

# Clip a convex polygon (x/y vertex vectors, in order) against an
# axis-aligned box, via Sutherland-Hodgman restricted to the box's 4 half-
# planes (always applicable here since the box is axis-aligned and every
# polygon clipped against it -- a stripe parallelogram -- is convex). Used
# instead of a line-drawing approximation (a thick line's "butt" end caps
# aren't parallel to the box's own edges, which leaves small diamond-shaped
# notches right at a rect's top/bottom boundary) so a diagonal stripe fill
# tiles a rect with true, exactly-straight edges and no gaps.
clip_polygon_to_box <- function(x, y, xmin, xmax, ymin, ymax) {
  clip_edge <- function(x, y, inside_fn, intersect_fn) {
    n <- length(x)
    if (n == 0) {
      return(list(x = numeric(0), y = numeric(0)))
    }
    out_x <- numeric(0)
    out_y <- numeric(0)
    for (i in seq_len(n)) {
      cx <- x[i]
      cy <- y[i]
      prev_i <- if (i == 1) n else i - 1
      px <- x[prev_i]
      py <- y[prev_i]
      c_in <- inside_fn(cx, cy)
      p_in <- inside_fn(px, py)
      if (c_in) {
        if (!p_in) {
          pt <- intersect_fn(px, py, cx, cy)
          out_x <- c(out_x, pt$x)
          out_y <- c(out_y, pt$y)
        }
        out_x <- c(out_x, cx)
        out_y <- c(out_y, cy)
      } else if (p_in) {
        pt <- intersect_fn(px, py, cx, cy)
        out_x <- c(out_x, pt$x)
        out_y <- c(out_y, pt$y)
      }
    }
    list(x = out_x, y = out_y)
  }

  r <- list(x = x, y = y)
  r <- clip_edge(r$x, r$y, function(px, py) px >= xmin, function(x0, y0, x1, y1) {
    list(x = xmin, y = y0 + (xmin - x0) / (x1 - x0) * (y1 - y0))
  })
  r <- clip_edge(r$x, r$y, function(px, py) px <= xmax, function(x0, y0, x1, y1) {
    list(x = xmax, y = y0 + (xmax - x0) / (x1 - x0) * (y1 - y0))
  })
  r <- clip_edge(r$x, r$y, function(px, py) py >= ymin, function(x0, y0, x1, y1) {
    list(x = x0 + (ymin - y0) / (y1 - y0) * (x1 - x0), y = ymin)
  })
  r <- clip_edge(r$x, r$y, function(px, py) py <= ymax, function(x0, y0, x1, y1) {
    list(x = x0 + (ymax - y0) / (y1 - y0) * (x1 - x0), y = ymax)
  })
  r
}

# Diagonal (45-degree) two-color stripe fill for a rect, as a set of
# polygons (one per stripe, exactly clipped to the rect) rather than a
# pattern-fill package -- matching this codebase's convention of building
# every calibrated figure from first-party grid/ggplot2 primitives. Each
# stripe is a long thick parallelogram (centered on a 45-degree line,
# generously overextended in both directions) clipped to the rect via
# clip_polygon_to_box(); adjacent stripe centerlines are spaced
# stripe_width_mm*sqrt(2) apart in (x-y) so their perpendicular spacing
# exactly equals stripe_width_mm, tiling the rect edge-to-edge with no gaps
# and no line-cap artifacts at the rect's own edges. `group_prefix` keeps
# each call's polygons distinct when several stripe fills are combined into
# one geom_polygon() layer.
diagonal_stripe_polygons <- function(xmin, xmax, ymin, ymax, colors, stripe_width_mm = STRIPE_WIDTH_MM,
                                      group_prefix = "") {
  w <- xmax - xmin
  h <- ymax - ymin
  half <- stripe_width_mm / 2
  step <- stripe_width_mm * sqrt(2)
  diag_unit <- c(1, 1) / sqrt(2)
  perp_unit <- c(1, -1) / sqrt(2)
  # Each stripe's centerline must reach at least `w` in x and `h` in y from
  # its anchor point -- but moving `extend` units along a 45-degree
  # direction only covers extend/sqrt(2) on either axis, so extend has to be
  # (w + h) *after* that sqrt(2) reduction, not before it. Using plain
  # `w + h` here (rather than `(w + h) * sqrt(2)`) under-covers any rect
  # whose height is much larger than its width (or vice versa) -- e.g. a
  # tall stacked-bar segment -- silently truncating the fill well short of
  # the rect's far edge.
  extend <- (w + h) * sqrt(2) + step
  c_values <- seq(-h - step, w + step, by = step)

  polys <- lapply(seq_along(c_values), function(i) {
    c_val <- c_values[i]
    p1 <- c(c_val, 0) - diag_unit * extend
    p2 <- c(c_val, 0) + diag_unit * extend
    off <- perp_unit * half
    corners_x <- c(p1[1] + off[1], p2[1] + off[1], p2[1] - off[1], p1[1] - off[1])
    corners_y <- c(p1[2] + off[2], p2[2] + off[2], p2[2] - off[2], p1[2] - off[2])
    clipped <- clip_polygon_to_box(corners_x, corners_y, 0, w, 0, h)
    if (length(clipped$x) < 3) {
      return(NULL)
    }
    tibble(
      x = xmin + clipped$x, y = ymin + clipped$y,
      group = paste0(group_prefix, "_", i),
      colour = colors[[(i %% length(colors)) + 1]]
    )
  })
  do.call(rbind, polys)
}

# summary_df: one row per (Gene, assay_type) with n_variants (the stacked
# segment's own height) and Gene already ordered as a factor in the desired
# left-to-right bar order (descending total, matching the original chart).
plot_variant_bar_calibrated <- function(summary_df,
                                         total_width_mm = 4.25 * 25.4, total_height_mm = 3.25 * 25.4,
                                         bar_width_frac = 0.85,
                                         y_title = "Total unique variants") {
  measure_dev_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(measure_dev_file)
  on.exit({
    grDevices::dev.off()
    unlink(measure_dev_file)
  }, add = TRUE)

  genes <- levels(summary_df$Gene)
  n_genes <- length(genes)

  # ---- Stack geometry (in data units; converted to mm once chart_width_mm/
  # chart_height_mm are known below) ----
  totals <- summary_df %>%
    group_by(Gene) %>%
    summarise(total = sum(n_variants), .groups = "drop")
  y_data_max <- max(totals$total, na.rm = TRUE)
  AXIS_END_PAD <- 1.05
  y_max <- y_data_max * AXIS_END_PAD
  y_breaks <- Filter(function(b) b <= y_max, pretty(c(0, totals$total)))
  y_tick_labels <- format(y_breaks, big.mark = ",", trim = TRUE)

  bars <- summary_df %>%
    mutate(assay_type = factor(assay_type, levels = ASSAY_STACK_ORDER)) %>%
    arrange(Gene, assay_type) %>%
    group_by(Gene) %>%
    mutate(ymax_data = cumsum(n_variants), ymin_data = lag(ymax_data, default = 0)) %>%
    ungroup() %>%
    mutate(gene_index = as.integer(Gene) - 1, fill_color = ASSAY_COLORS[as.character(assay_type)])

  # ---- Text measurements that don't depend on chart size ----
  y_tick_w <- max(text_width_mm(y_tick_labels, LABEL_PT))
  y_title_h <- text_height_mm(y_title, LABEL_PT) # rotated: rendered width = text height
  x_title <- "Gene"
  x_title_h <- text_height_mm(x_title, LABEL_PT)
  # Rotated 90 degrees (matching the original angle = 90 gene labels), so
  # their rendered *width* (horizontal) becomes their vertical extent below
  # the axis, not their height.
  gene_label_w <- max(text_width_mm(genes, LABEL_PT, fontface = "italic"))

  # ---- Legend: a horizontal row floating inside the plot area (swatch +
  # label per assay, no title/box), matching the original (preprint)
  # figure's in-plot legend placement rather than an outside right-hand
  # column. A second row (only when present in the data) adds one
  # diagonally-striped swatch per multi-assay combo.
  legend_swatch_mm <- text_height_mm("Mg", LABEL_PT) * 1.3
  legend_item_gap_mm <- GAP_MM * 2
  legend_row_gap_mm <- GAP_MM

  legend_labels <- ASSAY_LEGEND_DISPLAY[ASSAY_LEGEND_ORDER]
  legend_label_w <- text_width_mm(legend_labels, LABEL_PT)
  legend_item_w <- legend_swatch_mm + GAP_MM / 2 + legend_label_w
  n_legend <- length(ASSAY_LEGEND_ORDER)

  present_combos <- intersect(ASSAY_COMBO_ORDER, as.character(unique(bars$assay_type)))
  n_combo_legend <- length(present_combos)
  combo_legend_labels <- ASSAY_COMBO_DISPLAY[present_combos]
  combo_legend_label_w <- text_width_mm(combo_legend_labels, LABEL_PT)
  combo_legend_item_w <- legend_swatch_mm + GAP_MM / 2 + combo_legend_label_w

  # ---- Margins built backward from the fixed total canvas ----
  left_margin_mm <- y_title_h + GAP_MM + y_tick_w + GAP_MM
  bottom_margin_mm <- GAP_MM + gene_label_w + GAP_MM + x_title_h
  right_margin_mm <- GAP_MM
  top_margin_mm <- 2

  chart_width_mm <- total_width_mm - left_margin_mm - right_margin_mm
  chart_height_mm <- total_height_mm - top_margin_mm - bottom_margin_mm

  mm_y <- function(v) v / y_max * chart_height_mm
  slot_w <- chart_width_mm / n_genes
  bar_w <- slot_w * bar_width_frac

  bars <- bars %>%
    mutate(
      x_center = (gene_index + 0.5) * slot_w,
      xmin = x_center - bar_w / 2, xmax = x_center + bar_w / 2,
      ymin = mm_y(ymin_data), ymax = mm_y(ymax_data)
    )

  # Segments whose assay_type is a "+"-joined combo (a variant tested by
  # more than one assay type) render as a diagonal two-color stripe instead
  # of a solid fill.
  bars_solid <- bars %>% filter(!grepl("\\+", as.character(assay_type)))
  bars_combo <- bars %>% filter(grepl("\\+", as.character(assay_type)))

  combo_stripe_polygons <- if (nrow(bars_combo) > 0) {
    do.call(rbind, lapply(seq_len(nrow(bars_combo)), function(i) {
      row <- bars_combo[i, ]
      parts <- strsplit(as.character(row$assay_type), "\\+")[[1]]
      diagonal_stripe_polygons(row$xmin, row$xmax, row$ymin, row$ymax, ASSAY_COLORS[parts], group_prefix = paste0("bar", i))
    }))
  } else {
    tibble(x = numeric(0), y = numeric(0), group = character(0), colour = character(0))
  }

  # ---- Axis ticks/titles (positioned outside the panel, at x=0/y=0) ----
  y_tick_y <- mm_y(y_breaks)
  y_tick_right_x <- 0 - GAP_MM
  y_tick_left_x <- y_tick_right_x - y_tick_w
  y_title_x <- y_tick_left_x - GAP_MM - y_title_h / 2

  gene_tick_x <- (seq_len(n_genes) - 0.5) * slot_w
  gene_tick_top_y <- 0 - GAP_MM
  x_title_y <- gene_tick_top_y - gene_label_w - GAP_MM - x_title_h / 2

  # ---- Legend geometry: floating inside the plot area in the "elbow"
  # above the shorter bars (right-of-center, clear of the tallest bars on
  # the left) -- a guide-based match to the original figure's in-plot
  # legend, not a pixel-exact reproduction (the original's exact spot is
  # specific to its own gene ordering/heights).
  legend_total_w <- sum(legend_item_w) + (n_legend - 1) * legend_item_gap_mm
  legend_right_x <- chart_width_mm * 0.95
  legend_left_x <- legend_right_x - legend_total_w
  legend_item_left <- legend_left_x + cumsum(c(0, legend_item_w[-n_legend] + legend_item_gap_mm))
  legend_row1_y <- chart_height_mm * 0.5

  legend_swatches <- tibble(
    label = legend_labels,
    swatch_x = legend_item_left + legend_swatch_mm / 2,
    label_x = legend_item_left + legend_swatch_mm + GAP_MM / 2,
    y = legend_row1_y,
    fill = ASSAY_COLORS[ASSAY_LEGEND_ORDER]
  )

  if (n_combo_legend > 0) {
    combo_total_w <- sum(combo_legend_item_w) + (n_combo_legend - 1) * legend_item_gap_mm
    combo_left_x <- chart_width_mm * 0.95 - combo_total_w
    combo_item_left <- combo_left_x + cumsum(c(0, head(combo_legend_item_w, -1) + legend_item_gap_mm))
    legend_row2_y <- legend_row1_y - legend_swatch_mm - legend_row_gap_mm

    combo_legend_swatches <- tibble(
      label = combo_legend_labels,
      combo = present_combos,
      swatch_xmin = combo_item_left,
      swatch_xmax = combo_item_left + legend_swatch_mm,
      label_x = combo_item_left + legend_swatch_mm + GAP_MM / 2,
      y = legend_row2_y
    )
    combo_legend_polygons <- do.call(rbind, lapply(seq_len(n_combo_legend), function(i) {
      row <- combo_legend_swatches[i, ]
      parts <- strsplit(row$combo, "\\+")[[1]]
      diagonal_stripe_polygons(
        row$swatch_xmin, row$swatch_xmax, row$y - legend_swatch_mm / 2, row$y + legend_swatch_mm / 2,
        ASSAY_COLORS[parts],
        stripe_width_mm = STRIPE_WIDTH_MM * 0.7,
        group_prefix = paste0("legend", i)
      )
    }))
  } else {
    combo_legend_swatches <- tibble(label = character(0), label_x = numeric(0), y = numeric(0))
    combo_legend_polygons <- tibble(x = numeric(0), y = numeric(0), group = character(0), colour = character(0))
  }

  x_range <- c(-left_margin_mm, chart_width_mm + right_margin_mm)
  y_range <- c(-bottom_margin_mm, chart_height_mm + top_margin_mm)

  p <- ggplot() +
    geom_rect(
      data = bars_solid, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill_color))
    ) +
    geom_polygon(
      data = combo_stripe_polygons, inherit.aes = FALSE,
      aes(x = x, y = y, group = group, fill = I(colour))
    ) +
    # Axis lines (bottom + left only, matching theme_classic()'s look).
    annotate("segment", x = 0, xend = 0, y = 0, yend = chart_height_mm, linewidth = AXIS_LWD, colour = "black") +
    annotate("segment", x = 0, xend = chart_width_mm, y = 0, yend = 0, linewidth = AXIS_LWD, colour = "black") +
    annotate("segment",
      x = -TICK_LEN_MM, xend = 0, y = y_tick_y, yend = y_tick_y,
      linewidth = AXIS_LWD, colour = "black"
    ) +
    annotate("text",
      x = y_tick_right_x, y = y_tick_y, label = y_tick_labels,
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 1
    ) +
    annotate("text",
      x = y_title_x, y = chart_height_mm / 2, label = y_title, angle = 90,
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black"
    ) +
    # Gene labels: rotated 90 degrees, hjust = 1 so they hang downward from
    # the tick (not upward into the bars), vjust = 0.5 to stay centered
    # under their own bar -- same convention the original chart used.
    annotate("text",
      x = gene_tick_x, y = gene_tick_top_y, label = genes, angle = 90, hjust = 1, vjust = 0.5,
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black", fontface = "italic"
    ) +
    annotate("text",
      x = chart_width_mm / 2, y = x_title_y, label = x_title,
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black"
    ) +
    # ---- Legend: horizontal row(s) floating inside the plot area ----
    geom_tile(
      data = legend_swatches, inherit.aes = FALSE,
      aes(x = swatch_x, y = y, fill = I(fill)),
      width = legend_swatch_mm, height = legend_swatch_mm
    ) +
    geom_text(
      data = legend_swatches, inherit.aes = FALSE,
      aes(x = label_x, y = y, label = label),
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 0, vjust = 0.5
    ) +
    # Multi-assay combo row (only non-empty when the data actually has one).
    geom_polygon(
      data = combo_legend_polygons, inherit.aes = FALSE,
      aes(x = x, y = y, group = group, fill = I(colour))
    ) +
    geom_text(
      data = combo_legend_swatches, inherit.aes = FALSE,
      aes(x = label_x, y = y, label = label),
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 0, vjust = 0.5
    ) +
    coord_cartesian(xlim = x_range, ylim = y_range, expand = FALSE, clip = "off") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

  list(plot = p, width_mm = diff(x_range), height_mm = diff(y_range))
}
