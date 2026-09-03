# Shared calibrated (real-print-size) grouped horizontal bar chart for Fig
# 5's ClinVar/ClinGen "controls metrics" panel (make_confusion_matrix_*'s
# own $metrics data frame, visualized via plot_clinvar_metrics() in
# Figure5_6.Rmd) -- reproduced here with exact physical dimensions,
# replicating the manually-edited PDF's abbreviated tick labels ("Ccd"/
# "Sens"/"Spec"/"MCC"/"Det"), axis titles, tick marks, and 2-row legend
# (none of which the original ggplot code renders -- it uses
# axis.title = element_blank() and legend.position = "none").
#
# Like confusion_matrix_calibrated.R, there's no data-driven geometry to
# extract from an existing ggplot -- bar positions are entirely fixed by
# this project's own dodge/width settings (position_dodge(width = 0.7),
# geom_col(width = 0.6), 3 predictors, 5 metrics). Confirmed once against
# ggplot_build() on the original plot_clinvar_metrics() code (not
# re-derived from ggplot2's internals at runtime, since the exact
# dodge/width-splitting formula isn't documented) and hardcoded here.
#
# Expects the caller to have already loaded: ggplot2, extrafont (+
# loadfonts(device = 'all')), and grid.

point_in_mm <- 0.3527778
TICK_PT <- 5
TITLE_PT <- 6
FONT_FAMILY <- "Arial"
# Same reused-gap convention as this project's other calibrated-chart
# libraries (see confusion_matrix_calibrated.R's own GAP_MM comment).
GAP_MM <- 1
# Small breathing room between a bar's own end and its value label --
# deliberately smaller than GAP_MM, matching confusion_matrix_calibrated.R's
# own AXIS_GAP_MM-vs-GAP_MM distinction.
VALUE_LABEL_GAP_MM <- 0.5
TICK_LEN_MM <- 0.8
AXIS_LWD <- 0.5 * point_in_mm

# Bottom-to-top metric order and abbreviated tick labels, matching
# plot_clinvar_metrics()'s own factor levels (c("Determinate_calls", "MCC",
# "Specificity", "Sensitivity", "Concordance")) and the manually-edited
# PDF's abbreviations.
CATEGORY_ORDER <- c("Determinate_calls", "MCC", "Specificity", "Sensitivity", "Concordance")
CATEGORY_DISPLAY <- c(
  Determinate_calls = "Det", MCC = "MCC", Specificity = "Spec",
  Sensitivity = "Sens", Concordance = "Ccd"
)
# Bottom-to-top order within each metric's 3-bar cluster -- matches
# ggplot2's own alphabetical dodge order for an unordered fill factor
# (confirmed via ggplot_build() on the original plot).
PREDICTOR_ORDER <- c("AlphaMissense", "MutPred2", "REVEL")
PREDICTOR_COLORS <- c(REVEL = "black", AlphaMissense = "darkgrey", MutPred2 = "whitesmoke")
# Legend reading order (row-major, 2 per row): matches the bar-cluster's
# own top-to-bottom order, not alphabetical -- reproduces the reference
# screenshot's manually-arranged legend.
LEGEND_ORDER <- c("REVEL", "MutPred2", "AlphaMissense")

# geom_col(position = position_dodge(width = 0.7), width = 0.6) with 3
# groups: dodge divides the 0.7-wide cluster into 3 evenly-spaced slots,
# and the bar's own width (0.6) is *also* divided by the group count for
# its rendered thickness -- confirmed via ggplot_build() on the original
# plot rather than re-derived from ggplot2's internals.
N_PREDICTORS <- length(PREDICTOR_ORDER)
DODGE_SLOT_UNITS <- 0.7 / N_PREDICTORS
BAR_HALF_THICKNESS_UNITS <- 0.6 / N_PREDICTORS / 2
# Default discrete-scale expansion (ggplot2's own expansion(add = 0.6)),
# the reason a stock discrete axis always shows whitespace above/below the
# outermost categories -- reproduced here since the calibrated chart should
# keep the current PNG's own bar padding.
DISCRETE_EXPAND_UNITS <- 0.6

text_width_mm <- function(label, pt, family = FONT_FAMILY) {
  gp <- gpar(fontfamily = family, fontsize = pt)
  vapply(label, function(l) convertWidth(grobWidth(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

text_height_mm <- function(label, pt, family = FONT_FAMILY) {
  gp <- gpar(fontfamily = family, fontsize = pt)
  vapply(label, function(l) convertHeight(grobHeight(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

# metrics_df: one row per predictor, with a "Predictor" column (values in
# PREDICTOR_ORDER) and one column per CATEGORY_ORDER entry -- exactly the
# make_confusion_matrix_clinvar_chek2()/make_confusion_matrix_clingen()
# $metrics data frame (after Predictor is set), i.e. what
# plot_clinvar_metrics() itself already takes. chart_width_mm/
# chart_height_mm size the 0-to-1-by-x / bottom-to-top-by-y core plot area
# only -- value labels, tick labels, axis titles, and the legend all extend
# beyond it.
make_metrics_bar_calibrated <- function(metrics_df,
                                         chart_width_mm = 26, chart_height_mm = 35,
                                         tick_pt = TICK_PT, title_pt = TITLE_PT) {
  measure_dev_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(measure_dev_file)
  on.exit({
    grDevices::dev.off()
    unlink(measure_dev_file)
  }, add = TRUE)

  n_categories <- length(CATEGORY_ORDER)
  y_unit_min <- 1 - DISCRETE_EXPAND_UNITS
  y_unit_max <- n_categories + DISCRETE_EXPAND_UNITS
  mm_y <- function(u) (u - y_unit_min) / (y_unit_max - y_unit_min) * chart_height_mm
  mm_x <- function(v) v * chart_width_mm

  # ---- Bars ----
  bars <- do.call(rbind, lapply(seq_along(CATEGORY_ORDER), function(i) {
    category <- CATEGORY_ORDER[i]
    do.call(rbind, lapply(seq_along(PREDICTOR_ORDER), function(j) {
      predictor <- PREDICTOR_ORDER[j]
      value <- metrics_df[[category]][metrics_df$Predictor == predictor]
      offset_units <- (j - (N_PREDICTORS + 1) / 2) * DODGE_SLOT_UNITS
      center_units <- i + offset_units
      data.frame(
        category = category, predictor = predictor, value = value,
        xmin = 0, xmax = mm_x(value),
        ymin = mm_y(center_units - BAR_HALF_THICKNESS_UNITS),
        ymax = mm_y(center_units + BAR_HALF_THICKNESS_UNITS),
        stringsAsFactors = FALSE
      )
    }))
  }))
  bars$fill_color <- PREDICTOR_COLORS[bars$predictor]
  bars$value_label <- sprintf("%.2f", bars$value)
  bars$label_x <- bars$xmax + VALUE_LABEL_GAP_MM
  bars$label_y <- (bars$ymin + bars$ymax) / 2
  bars$label_w <- text_width_mm(bars$value_label, tick_pt)

  # ---- Axis ticks/labels ----
  y_tick_labels <- unname(CATEGORY_DISPLAY[CATEGORY_ORDER])
  y_tick_y <- mm_y(seq_along(CATEGORY_ORDER))
  y_tick_w <- max(text_width_mm(y_tick_labels, tick_pt))

  x_tick_values <- c(0, 0.25, 0.5, 0.75, 1)
  x_tick_labels <- c("0", "0.25", "0.50", "0.75", "1")
  x_tick_x <- mm_x(x_tick_values)
  x_tick_h <- max(text_height_mm(x_tick_labels, tick_pt))

  y_title_h <- text_height_mm("Metric", title_pt) # rotated: rendered width = text height
  x_title_h <- text_height_mm("Value", title_pt)

  y_tick_x <- 0 - GAP_MM - y_tick_w / 2
  y_title_x <- 0 - GAP_MM - y_tick_w - GAP_MM - y_title_h / 2
  x_tick_y <- 0 - GAP_MM - x_tick_h / 2
  x_title_y <- 0 - GAP_MM - x_tick_h - GAP_MM - x_title_h / 2

  # ---- Legend: 2 rows (REVEL + MutPred2, then AlphaMissense alone),
  # matching the reference screenshot's own manual layout -- reading order
  # is each bar cluster's own top-to-bottom order (LEGEND_ORDER), wrapped 2
  # per row and left-aligned under the y-axis line.
  swatch_size <- text_height_mm("Mg", title_pt)
  legend_label_w <- text_width_mm(LEGEND_ORDER, title_pt)
  names(legend_label_w) <- LEGEND_ORDER
  item_w <- swatch_size + GAP_MM / 2 + legend_label_w
  row1_items <- LEGEND_ORDER[1:2]
  row2_items <- LEGEND_ORDER[3]
  row1_x <- setNames(c(0, item_w[[row1_items[1]]] + GAP_MM * 1.5), row1_items)
  row2_x <- setNames(0, row2_items)
  legend_row_h <- max(swatch_size, text_height_mm(LEGEND_ORDER, title_pt))
  legend_row_gap <- GAP_MM / 2

  legend_top <- x_title_y - x_title_h / 2 - GAP_MM
  legend_row1_y <- legend_top - legend_row_h / 2
  legend_row2_y <- legend_row1_y - legend_row_h - legend_row_gap
  legend_bottom <- legend_row2_y - legend_row_h / 2
  legend_right <- max(row1_x[[row1_items[2]]] + item_w[[row1_items[2]]], row2_x[[row2_items[1]]] + item_w[[row2_items[1]]])

  # ---- Canvas ----
  bars_right_edge <- max(bars$label_x + bars$label_w)
  x_range <- c(y_title_x - y_title_h / 2, max(bars_right_edge, legend_right)) + c(-GAP_MM, GAP_MM)
  y_range <- c(legend_bottom, chart_height_mm) + c(-GAP_MM, GAP_MM)

  # ---- Assemble ----
  legend_items <- do.call(rbind, lapply(LEGEND_ORDER, function(name) {
    is_row1 <- name %in% row1_items
    x <- if (is_row1) row1_x[[name]] else row2_x[[name]]
    y <- if (is_row1) legend_row1_y else legend_row2_y
    data.frame(
      predictor = name, swatch_x = x, swatch_y = y,
      label_x = x + swatch_size + GAP_MM / 2, label_y = y,
      stringsAsFactors = FALSE
    )
  }))
  legend_items$fill_color <- PREDICTOR_COLORS[legend_items$predictor]

  p <- ggplot() +
    geom_rect(
      data = bars, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
      colour = "black", linewidth = AXIS_LWD
    ) +
    scale_fill_identity() +
    geom_text(
      data = bars, inherit.aes = FALSE, aes(x = label_x, y = label_y, label = value_label),
      size = tick_pt * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 0
    ) +
    # Axis lines (bottom + left only, matching theme_classic()'s look).
    annotate("segment", x = 0, xend = 0, y = 0, yend = chart_height_mm, linewidth = AXIS_LWD, colour = "black") +
    annotate("segment", x = 0, xend = chart_width_mm, y = 0, yend = 0, linewidth = AXIS_LWD, colour = "black") +
    # Tick marks.
    annotate("segment",
      x = -TICK_LEN_MM, xend = 0, y = y_tick_y, yend = y_tick_y,
      linewidth = AXIS_LWD, colour = "black"
    ) +
    annotate("segment",
      x = x_tick_x, xend = x_tick_x, y = 0, yend = -TICK_LEN_MM,
      linewidth = AXIS_LWD, colour = "black"
    ) +
    # Tick labels.
    annotate("text",
      x = y_tick_x, y = y_tick_y, label = y_tick_labels,
      size = tick_pt * point_in_mm, family = FONT_FAMILY, colour = "black"
    ) +
    annotate("text",
      x = x_tick_x, y = x_tick_y, label = x_tick_labels,
      size = tick_pt * point_in_mm, family = FONT_FAMILY, colour = "black"
    ) +
    # Axis titles.
    annotate("text",
      x = y_title_x, y = chart_height_mm / 2, label = "Metric", angle = 90,
      size = title_pt * point_in_mm, family = FONT_FAMILY, colour = "black"
    ) +
    annotate("text",
      x = chart_width_mm / 2, y = x_title_y, label = "Value",
      size = title_pt * point_in_mm, family = FONT_FAMILY, colour = "black"
    ) +
    # Legend.
    geom_rect(
      data = legend_items, inherit.aes = FALSE,
      aes(
        xmin = swatch_x, xmax = swatch_x + swatch_size,
        ymin = swatch_y - swatch_size / 2, ymax = swatch_y + swatch_size / 2,
        fill = fill_color
      ),
      colour = "black", linewidth = AXIS_LWD
    ) +
    geom_text(
      data = legend_items, inherit.aes = FALSE, aes(x = label_x, y = label_y, label = predictor),
      size = title_pt * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 0
    ) +
    coord_cartesian(xlim = x_range, ylim = y_range, expand = FALSE, clip = "off") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

  list(plot = p, width_mm = diff(x_range), height_mm = diff(y_range))
}
