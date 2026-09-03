# Shared calibrated (real-print-size) 2x2 confusion-matrix chart logic --
# Fig 5's and Extended Data Fig 4/6's ClinVar/ClinGen-vs-predictor
# confusion matrices (make_confusion_matrix_clinvar_chek2()/
# make_confusion_matrix_clingen()'s own plot, reproduced here with exact
# physical dimensions and simplified axis/title text -- see the docstring
# above make_confusion_matrix_calibrated() below).
#
# Unlike sankey_calibrated.R, there's no data-driven geometry to extract
# from an existing ggplot -- a 2x2 grid's own layout is entirely fixed by
# matrix_width_mm/matrix_height_mm, so this works directly in millimeters
# as the plot's own data coordinate system throughout (no mm-per-data-unit
# conversion needed, unlike ggsankey's own arbitrary coordinate system).
#
# Expects the caller to have already loaded: ggplot2, extrafont (+
# loadfonts(device = 'all')), and grid.

point_in_mm <- 0.3527778
LABEL_PT <- 6
LABEL_MM <- LABEL_PT * point_in_mm
FONT_FAMILY <- "Arial"
# Gap between adjacent layout elements (matrix edge to tick labels, tick
# labels to axis label, matrix top to title, ...), and around the whole
# canvas -- both purposes reuse the same constant since a border stroke
# drawn centered on the canvas edge needs the same kind of small margin
# that separates any other two adjacent elements (see sankey_calibrated.R's
# own CANVAS_MARGIN_MM comment for the identical reasoning).
GAP_MM <- 1
# Small breathing room between the matrix's own bottom/left edge and the
# axis lines framing it (matching the reference screenshot, where the axis
# lines sit slightly clear of the colored cells rather than touching them
# directly) -- deliberately smaller than GAP_MM, which spaces the axis
# lines from the tick labels beyond them.
AXIS_GAP_MM <- 1
# grid's "lwd" is device-dependent but conventionally 1/96 inch per unit;
# convert our usual hairline (0.5pt, matching every other calibrated
# chart's own border convention) into that unit, for the axis line grob.
BORDER_LWD <- (0.5 * point_in_mm) / (25.4 / 96)

text_width_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY) {
  gp <- gpar(fontfamily = family, fontsize = pt)
  vapply(label, function(l) convertWidth(grobWidth(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

text_height_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY) {
  gp <- gpar(fontfamily = family, fontsize = pt)
  vapply(label, function(l) convertHeight(grobHeight(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

format_count <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# TP/TN/FP/FN: the 2x2 grid's own cell counts, matching
# make_confusion_matrix_clinvar_chek2()/make_confusion_matrix_clingen()'s
# own TP (x=Pathogenic, y=Pathogenic) / TN (x=Benign, y=Benign) / FP
# (x=Benign, y=Pathogenic) / FN (x=Pathogenic, y=Benign) convention exactly
# -- x is the ClinVar/ClinGen classification, y is "Exp + Pred
# classification". total_n is the *same* denominator those functions' own
# "Error"/discordance-rate figure already uses (nrow(df) as of their own
# return value, i.e. after any of their own internal filtering, such as
# ClinGen's VUS exclusion) -- pass the total_n field from that function's
# own returned list here rather than an independently-computed nrow(),
# since re-deriving it from a differently-filtered data frame would
# silently drift from the number the discordance-rate figure itself is
# built on. x_label names the ClinVar/ClinGen column ("ClinVar
# classification"/"ClinGen classification"); y_label defaults to "Exp +
# Pred classification" (dropping the predictor name the original plot's
# own y-axis carried -- this calibrated version relies on the surrounding
# figure/filename to convey which predictor it's showing instead, matching
# every other calibrated chart in this project).
make_confusion_matrix_calibrated <- function(TP, TN, FP, FN, total_n,
                                              matrix_width_mm, matrix_height_mm,
                                              x_label, y_label = "Exp + Pred classification",
                                              x_tick_labels = c("B/LB", "P/LP"),
                                              y_tick_labels = c("B/LB", "P/LP"),
                                              label_pt = LABEL_PT) {
  measure_dev_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(measure_dev_file)
  on.exit({
    grDevices::dev.off()
    unlink(measure_dev_file)
  }, add = TRUE)

  label_mm <- label_pt * point_in_mm

  discordance_rate <- (FP + FN) / total_n * 100
  title_text <- sprintf("Total = %s,\ndiscordance rate = %.2f%%", format_count(total_n), discordance_rate)

  # Cell layout: bottom-left origin, x increasing right (Benign column then
  # Pathogenic column), y increasing up (Benign row then Pathogenic row) --
  # matching the original plot's own alphabetical/factor-level ordering
  # (Benign before Pathogenic on both axes).
  half_w <- matrix_width_mm / 2
  half_h <- matrix_height_mm / 2
  cells <- data.frame(
    xmin = c(0, half_w, 0, half_w),
    xmax = c(half_w, matrix_width_mm, half_w, matrix_width_mm),
    ymin = c(half_h, half_h, 0, 0),
    ymax = c(matrix_height_mm, matrix_height_mm, half_h, half_h),
    count = c(FP, TP, TN, FN),
    is_diagonal = c(FALSE, TRUE, TRUE, FALSE)
  )
  cells$fill_color <- ifelse(
    cells$count == 0, "white",
    ifelse(cells$is_diagonal, "#D8BFD8", "#F2F2F2")
  )
  cells$count_label <- format_count(cells$count)
  cells$center_x <- (cells$xmin + cells$xmax) / 2
  cells$center_y <- (cells$ymin + cells$ymax) / 2

  # Layout: title above the matrix; y-axis label (rotated 90 degrees, so
  # its own *height* becomes the width it occupies) and y tick labels to
  # its left, in a column to the left of the matrix; x-axis label and x
  # tick labels below the matrix, same idea but unrotated.
  title_h_mm <- text_height_mm(title_text, pt = label_pt)
  # y_label = NULL omits both the text and the gap/column of space it would
  # otherwise reserve to the left of the y tick labels -- the formulas below
  # collapse to the correct tightened layout when both are 0.
  y_label_w_mm <- if (is.null(y_label)) 0 else text_height_mm(y_label, pt = label_pt)
  y_label_gap_mm <- if (is.null(y_label)) 0 else GAP_MM
  y_tick_w_mm <- max(text_width_mm(y_tick_labels, pt = label_pt))
  x_label_h_mm <- text_height_mm(x_label, pt = label_pt)
  x_tick_h_mm <- max(text_height_mm(x_tick_labels, pt = label_pt))
  # Both centered over the matrix's own width -- at matrix_width_mm large
  # relative to label_pt this is always narrower than the matrix itself, so
  # it was previously left out of x_range entirely; a small matrix_width_mm
  # (e.g. a 50%-scaled chart, same unchanged font size) can flip that,
  # letting the title/x_label overflow past the matrix's own left/right
  # edges -- accounted for explicitly below rather than assumed away.
  title_w_mm <- max(text_width_mm(strsplit(title_text, "\n")[[1]], pt = label_pt))
  x_label_w_mm <- text_width_mm(x_label, pt = label_pt)

  y_tick_x <- 0 - AXIS_GAP_MM - GAP_MM - y_tick_w_mm / 2
  y_label_x <- 0 - AXIS_GAP_MM - GAP_MM - y_tick_w_mm - y_label_gap_mm - y_label_w_mm / 2
  x_tick_y <- 0 - AXIS_GAP_MM - GAP_MM - x_tick_h_mm / 2
  x_label_y <- 0 - AXIS_GAP_MM - GAP_MM - x_tick_h_mm - GAP_MM - x_label_h_mm / 2
  title_y <- matrix_height_mm + GAP_MM + title_h_mm / 2
  matrix_center_x <- matrix_width_mm / 2

  x_range <- c(
    min(y_label_x - y_label_w_mm / 2, matrix_center_x - title_w_mm / 2, matrix_center_x - x_label_w_mm / 2),
    max(matrix_width_mm, matrix_center_x + title_w_mm / 2, matrix_center_x + x_label_w_mm / 2)
  ) + c(-GAP_MM, GAP_MM)
  y_range <- c(x_label_y - x_label_h_mm / 2, title_y + title_h_mm / 2) + c(-GAP_MM, GAP_MM)

  p <- ggplot() +
    geom_rect(
      data = cells, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
      colour = "white", linewidth = 0.5 * point_in_mm
    ) +
    scale_fill_identity() +
    # Bottom + left axis lines only (theme_classic()'s own look), no tick
    # marks -- just the two line segments framing the matrix, held off its
    # edge by AXIS_GAP_MM rather than touching the cells directly, and
    # sharing the (-AXIS_GAP_MM, -AXIS_GAP_MM) corner so they actually meet
    # at the origin instead of each stopping short at the other axis's own
    # unoffset position.
    annotate("segment", x = -AXIS_GAP_MM, xend = matrix_width_mm, y = -AXIS_GAP_MM, yend = -AXIS_GAP_MM,
              linewidth = 0.5 * point_in_mm, colour = "black") +
    annotate("segment", x = -AXIS_GAP_MM, xend = -AXIS_GAP_MM, y = -AXIS_GAP_MM, yend = matrix_height_mm,
              linewidth = 0.5 * point_in_mm, colour = "black") +
    geom_text(
      data = cells, inherit.aes = FALSE,
      aes(x = center_x, y = center_y, label = count_label),
      size = label_mm, family = FONT_FAMILY, colour = "black"
    ) +
    annotate("text", x = matrix_width_mm / 2, y = title_y, label = title_text,
              size = label_mm, family = FONT_FAMILY, colour = "black", lineheight = 1) +
    annotate("text", x = matrix_width_mm / 2, y = x_label_y, label = x_label,
              size = label_mm, family = FONT_FAMILY, colour = "black") +
    annotate("text", x = matrix_width_mm / 4, y = x_tick_y, label = x_tick_labels[1],
              size = label_mm, family = FONT_FAMILY, colour = "black") +
    annotate("text", x = matrix_width_mm * 3 / 4, y = x_tick_y, label = x_tick_labels[2],
              size = label_mm, family = FONT_FAMILY, colour = "black") +
    annotate("text", x = y_tick_x, y = matrix_height_mm / 4, label = y_tick_labels[1],
              size = label_mm, family = FONT_FAMILY, colour = "black") +
    annotate("text", x = y_tick_x, y = matrix_height_mm * 3 / 4, label = y_tick_labels[2],
              size = label_mm, family = FONT_FAMILY, colour = "black") +
    coord_cartesian(xlim = x_range, ylim = y_range, expand = FALSE, clip = "off") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

  if (!is.null(y_label)) {
    p <- p + annotate("text", x = y_label_x, y = matrix_height_mm / 2, label = y_label,
                       size = label_mm, family = FONT_FAMILY, colour = "black", angle = 90)
  }

  list(plot = p, width_mm = diff(x_range), height_mm = diff(y_range))
}
