# Shared calibrated bubble/scatter chart logic (Figure 3a -- clinical tests
# vs. possible SNVs, sized by ClinVar VUS count, colored by IGVF-produced
# status, with a red ring marking ACMG secondary-findings genes).
#
# Renders at the *exact* final placed size (cairo_pdf, mm units, Arial),
# matching the convention in metrics_bar_calibrated.R/
# confusion_matrix_calibrated.R/three_ring_donut.R -- the whole plot,
# including axis titles/ticks and the legend, is laid out directly in mm so
# that at 100% placement scale the labeled point size matches the requested
# point size. Point labels still use ggrepel::geom_text_repel() for their
# established force-directed overlap avoidance -- that algorithm isn't worth
# reimplementing, and its own padding/force parameters are given in real mm
# via grid::unit(), so it works correctly in this mm-space coordinate system.
#
# Expects the caller to have already loaded: dplyr, ggplot2, ggrepel,
# extrafont (+ loadfonts(device = 'all')), and grid.

point_in_mm <- 0.3527778
FONT_FAMILY <- "Arial"
LABEL_PT <- 7
# Slightly smaller than LABEL_PT (axis/legend text) -- kept as its own
# constant since the two are tuned independently.
GENE_LABEL_PT <- 6

GAP_MM <- 1.2
TICK_LEN_MM <- 0.8
AXIS_LWD <- 0.5 * point_in_mm

IGVF_COLORS <- c(
  "IGVF-produced" = "#176082",
  "IGVF-produced and community data" = "#7B9EAC",
  "Other" = "grey"
)
ACMG_STROKE_COLOR <- "red3"

text_width_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY, fontface = "plain") {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = fontface)
  vapply(label, function(l) convertWidth(grobWidth(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

text_height_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY, fontface = "plain") {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = fontface)
  vapply(label, function(l) convertHeight(grobHeight(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

# Matches ggplot2's scale_size_continuous() default area_pal(): the size
# *aesthetic* (point diameter in mm) is derived so that point AREA, not
# diameter, is linear in the data value -- reproduced manually here since the
# actual points and the manually-drawn size-legend circles both need the
# exact same value-to-mm mapping.
size_mm_for_value <- function(value, data_range, mm_range) {
  area_range <- mm_range^2
  area <- area_range[1] + (value - data_range[1]) / diff(data_range) * diff(area_range)
  sqrt(pmax(area, 0))
}

# gene_df: one row per gene, with gene_test_count (x), possible_SNVs (y),
# Uncertain.significance (VUS count -> point size), IGVF_produced (one of
# IGVF_COLORS' names -> fill color), and gene_symbol (ACMG-membership +
# point labels).
# size_legend_values are the two example VUS counts shown in the size
# legend (500/1,000 in the published figure) -- not derived from the data,
# since they're meant as fixed, easy-to-read reference points rather than
# data-driven breaks.
#
# chart_width_mm/chart_height_mm size the plot *area* only -- the rectangle
# from the origin to the far edges of the x/y axis lines -- matching how the
# published figure is described (3in x 2in there), narrowed by ~2 character
# widths at LABEL_PT (2 * width of "0"/"n" at 7pt Arial =~ 2.82mm) per a
# later request to trim the chart's overall width. Axis titles, tick labels,
# and the legend all extend beyond that rectangle, so the final returned
# width_mm/height_mm (what actually gets passed to ggsave) end up larger than
# chart_width_mm/chart_height_mm by however much those need.
plot_bubble_calibrated <- function(gene_df, acmg_sf_genes, label_genes,
                                    chart_width_mm = 3 * 25.4 - 2.82, chart_height_mm = 2 * 25.4,
                                    size_range_mm = c(1.3, 7.5) * 2 / 3,
                                    size_legend_values = c(500, 1000),
                                    top_margin_mm = 9, right_margin_mm = 11) {
  measure_dev_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(measure_dev_file)
  on.exit({
    grDevices::dev.off()
    unlink(measure_dev_file)
  }, add = TRUE)

  gene_df <- gene_df %>%
    mutate(
      is_acmg = gene_symbol %in% acmg_sf_genes,
      has_label = gene_symbol %in% label_genes,
      stroke_color = if_else(is_acmg, ACMG_STROKE_COLOR, "black"),
      stroke_width = if_else(is_acmg, 0.9, 0),
      fill_color = IGVF_COLORS[IGVF_produced]
    )

  size_data_range <- range(gene_df$Uncertain.significance, na.rm = TRUE)
  gene_df$point_size_mm <- size_mm_for_value(gene_df$Uncertain.significance, size_data_range, size_range_mm)

  # Axis line ends a little past the actual data max, not padded out to the
  # next round pretty() break (matching the original preprint figure's own
  # axes, which end between tick marks rather than exactly on one) -- so
  # pretty() breaks beyond that actual endpoint are dropped entirely rather
  # than drawn right at the axis's edge.
  AXIS_END_PAD <- 1.05
  x_data_max <- max(gene_df$gene_test_count, na.rm = TRUE)
  y_data_max <- max(gene_df$possible_SNVs, na.rm = TRUE)
  x_max <- x_data_max * AXIS_END_PAD
  y_max <- y_data_max * AXIS_END_PAD
  x_breaks <- Filter(function(b) b <= x_max, pretty(c(0, gene_df$gene_test_count)))
  y_breaks <- Filter(function(b) b <= y_max, pretty(c(0, gene_df$possible_SNVs)))

  # ---- Text measurements that don't depend on chart size ----
  x_tick_labels <- format(x_breaks, big.mark = ",", trim = TRUE)
  y_tick_labels <- format(y_breaks, big.mark = ",", trim = TRUE)
  x_tick_h <- max(text_height_mm(x_tick_labels, LABEL_PT))
  y_tick_w <- max(text_width_mm(y_tick_labels, LABEL_PT))
  x_title <- "Number of clinical tests"
  y_title <- "Number of possible SNVs"
  x_title_h <- text_height_mm(x_title, LABEL_PT)
  y_title_h <- text_height_mm(y_title, LABEL_PT) # rotated: rendered width = text height

  # ---- Legend: a left column (3 stacked, single-line rows -- the two
  # IGVF-status colors plus the ACMG ring, none word-wrapped) beside the
  # VUS-size group. ----
  LEFT_COLUMN_LABELS <- c(
    "IGVF-produced",
    "IGVF-produced and community data",
    "ACMG secondary finding gene"
  )
  SIZE_TITLE <- "Number of VUS (SNVs) in ClinVar"
  LEGEND_LINEHEIGHT <- 0.9

  swatch_mm <- text_height_mm("Mg", LABEL_PT)
  legend_line_h <- text_height_mm("Mg", LABEL_PT) * 1.05
  left_col_row_gap_mm <- GAP_MM * 0.6
  left_col_row_h <- max(swatch_mm, legend_line_h)
  left_col_h <- 3 * left_col_row_h + 2 * left_col_row_gap_mm

  left_col_label_w <- max(text_width_mm(LEFT_COLUMN_LABELS, LABEL_PT))
  left_col_w <- swatch_mm + GAP_MM / 2 + left_col_label_w
  size_title_w <- text_width_mm(SIZE_TITLE, LABEL_PT)
  size_circle_mm <- size_mm_for_value(size_legend_values, size_data_range, size_range_mm)
  size_value_labels <- format(size_legend_values, big.mark = ",", trim = TRUE)
  size_value_w <- text_width_mm(size_value_labels, LABEL_PT)

  # Chained left-to-right: circle 1, its own value label, then circle 2 and
  # its value label -- circle 2 must clear circle 1's *label text*, not just
  # circle 1 itself, or the two would overlap whenever the label is wider
  # than the small fixed gap between the circles.
  circle1_right <- size_circle_mm[1]
  label1_x <- circle1_right + GAP_MM / 2
  label1_right <- label1_x + size_value_w[1]
  circle2_left <- label1_right + GAP_MM
  circle2_center <- circle2_left + size_circle_mm[2] / 2
  label2_x <- circle2_left + size_circle_mm[2] + GAP_MM / 2
  size_circles_w <- label2_x + size_value_w[2]
  size_group_w <- max(size_title_w, size_circles_w)

  group_gap_mm <- GAP_MM * 1.5

  # ---- Margins built outward from the fixed chart_width_mm/chart_height_mm
  # (the plot area itself) -- the opposite direction from solving a panel
  # size backward out of a fixed total canvas: here the total canvas (this
  # function's returned width_mm/height_mm) is whatever results.
  # Left margin: distance from the y-axis line (x=0) to the axis title's own
  # left edge -- title, gap, tick labels (right-aligned against the axis,
  # hence hjust = 1 below), gap, axis line.
  left_margin_mm <- y_title_h + GAP_MM + y_tick_w + GAP_MM
  # bottom_margin_mm itself is computed further down, once the legend's
  # actual final geometry (legend_top_y, legend_circle_row_y) is known --
  # neither chart_width_mm/chart_height_mm nor mm_x/mm_y below depend on it,
  # so the ordering is safe.

  mm_x <- function(v) v / x_max * chart_width_mm
  mm_y <- function(v) v / y_max * chart_height_mm

  gene_df <- gene_df %>% mutate(x_mm = mm_x(gene_test_count), y_mm = mm_y(possible_SNVs))

  # ---- Axis ticks/titles (positioned outside the panel, at x=0/y=0) ----
  x_tick_x <- mm_x(x_breaks)
  y_tick_y <- mm_y(y_breaks)
  x_tick_label_y <- 0 - GAP_MM - x_tick_h / 2
  # y-tick labels are right-aligned (hjust = 1 below): the x position given
  # to annotate("text", hjust = 1, ...) is the text's own *right* edge, not
  # its center -- unlike every other hjust = 0.5 (default) label above/below.
  y_tick_right_x <- 0 - GAP_MM
  y_tick_left_x <- y_tick_right_x - y_tick_w
  x_title_y <- x_tick_label_y - x_tick_h / 2 - GAP_MM - x_title_h / 2
  y_title_x <- y_tick_left_x - GAP_MM - y_title_h / 2

  # ---- Legend item positions (row spans the chart's own width, left-aligned
  # under the y-axis line so it lines up with the panel like the reference).
  # The left column's 3 rows and the size-legend's own title are all
  # top-anchored at the same legend_top_y -- the column's rows via their own
  # row-center math below, the (single-line) size title via vjust = 1. ----
  legend_top_y <- x_title_y - x_title_h / 2 - GAP_MM * 1.5

  left_col_row1_y <- legend_top_y - left_col_row_h / 2
  left_col_row2_y <- left_col_row1_y - left_col_row_h - left_col_row_gap_mm
  left_col_row3_y <- left_col_row2_y - left_col_row_h - left_col_row_gap_mm
  left_col_row_y <- c(left_col_row1_y, left_col_row2_y, left_col_row3_y)

  # The size-legend's example circles sit half as far below the (1-line)
  # size title as the row's own text-block height would otherwise put them
  # -- i.e. half of the gap between the size title's actual bottom edge and
  # where the circles would land using the same GAP_MM/2 spacing a 2-line-
  # tall row would use below its own bottom edge (this figure's rows are
  # all single-line now, but the halved gap was tuned by eye against that
  # reference and is kept as-is rather than re-tuned).
  size_title_bottom <- legend_top_y - legend_line_h
  default_gap_below_size_title <- legend_line_h * LEGEND_LINEHEIGHT + GAP_MM / 2
  circles_gap <- default_gap_below_size_title / 2
  legend_circle_row_y <- size_title_bottom - circles_gap - max(size_circle_mm) / 2

  # Aligned with the y-axis tick labels' own left edge (not the axis line
  # itself, x = 0) so the legend block doesn't look indented relative to
  # them.
  left_col_x0 <- y_tick_left_x
  size_x0 <- left_col_x0 + left_col_w + group_gap_mm

  legend_swatches <- tibble(
    kind = c("igvf_all", "igvf_mixed"),
    xmin = left_col_x0, xmax = left_col_x0 + swatch_mm,
    ymin = left_col_row_y[1:2] - swatch_mm / 2, ymax = left_col_row_y[1:2] + swatch_mm / 2,
    fill = IGVF_COLORS[c("IGVF-produced", "IGVF-produced and community data")]
  )
  legend_labels <- tibble(
    x = left_col_x0 + swatch_mm + GAP_MM / 2,
    y = left_col_row_y,
    label = LEFT_COLUMN_LABELS
  )
  legend_acmg_ring <- tibble(
    x = left_col_x0 + swatch_mm / 2, y = left_col_row_y[3], size_mm = swatch_mm * 0.85
  )
  legend_size_title <- tibble(x = size_x0, y = legend_top_y, label = SIZE_TITLE)
  legend_size_circles <- tibble(
    x = size_x0 + c(size_circle_mm[1] / 2, circle2_center),
    y = legend_circle_row_y,
    size_mm = size_circle_mm,
    value_label = size_value_labels,
    label_x = size_x0 + c(label1_x, label2_x)
  )

  # Bottom margin: distance from the x-axis line (y=0) down to the legend
  # block's own bottom edge -- whichever sits lower, the left column's 3rd
  # row (e.g. "ACMG secondary finding gene", vjust = 0.5 centered) or the
  # size-legend's circles. text_height_mm()'s grobHeight-based measurement
  # doesn't reliably reserve room for descenders, so add a little extra
  # margin below rather than leaving the "g" in "finding"/"gene" sitting
  # with zero clearance above the canvas edge.
  bottom_pad_mm <- 1
  legend_bottom_edge <- min(legend_top_y - left_col_h, legend_circle_row_y - max(size_circle_mm) / 2)
  bottom_margin_mm <- (0 - legend_bottom_edge) + bottom_pad_mm

  # ---- Canvas -- width is chart_width_mm + right_margin_mm, unless the
  # legend (which sits below the chart, in the same horizontal band as the
  # right margin, so the two don't both need their own separate allowance)
  # needs more room than that, in which case the canvas grows to fit it
  # rather than clipping text.
  legend_right_edge <- size_x0 + max(size_title_w, size_circles_w)
  x_range <- c(-left_margin_mm, max(chart_width_mm + right_margin_mm, legend_right_edge))
  y_range <- c(-bottom_margin_mm, chart_height_mm + top_margin_mm)

  # A few labels sit in dense clusters where the shared default nudge
  # (rightward none, upward 6mm -- tuned for the common case) pulls them
  # into a neighbor's own label/callout line, making it hard to tell which
  # line belongs to which bubble. Each gets a one-off override, chosen to
  # send it toward the nearest open space instead:
  #  - BARD1 sits among KCNH2/JAG1/G6PD/F9; the shared upward nudge pulled
  #    it up into that cluster, crossing KCNH2's own line. Nudged right.
  #  - GCK's bubble sits under 1mm below G6PD's own bubble (nearly
  #    touching), so both labels competed for the same patch of space
  #    above. Nudged down instead, into the open space below.
  #  - TARDBP sits beside F9/RHO/DDX3X; nudged down-left, away from that
  #    cluster, into open space below-left of its own bubble.
  label_overrides <- tribble(
    ~gene_symbol, ~nudge_x_mm, ~nudge_y_mm,
    "BARD1", 6, 6,
    "GCK", 0, -6,
    "TARDBP", -6, -6
  )
  label_df <- gene_df %>%
    filter(has_label) %>%
    left_join(label_overrides, by = "gene_symbol") %>%
    mutate(
      nudge_x_mm = coalesce(nudge_x_mm, 0),
      nudge_y_mm = coalesce(nudge_y_mm, 6)
    )

  p <- ggplot() +
    # Axis lines (bottom + left only, matching theme_classic()'s look).
    annotate("segment", x = 0, xend = 0, y = 0, yend = chart_height_mm, linewidth = AXIS_LWD, colour = "black") +
    annotate("segment", x = 0, xend = chart_width_mm, y = 0, yend = 0, linewidth = AXIS_LWD, colour = "black") +
    annotate("segment",
      x = x_tick_x, xend = x_tick_x, y = 0, yend = -TICK_LEN_MM,
      linewidth = AXIS_LWD, colour = "black"
    ) +
    annotate("segment",
      x = -TICK_LEN_MM, xend = 0, y = y_tick_y, yend = y_tick_y,
      linewidth = AXIS_LWD, colour = "black"
    ) +
    annotate("text",
      x = x_tick_x, y = x_tick_label_y, label = x_tick_labels,
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black"
    ) +
    annotate("text",
      x = y_tick_right_x, y = y_tick_y, label = y_tick_labels,
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 1
    ) +
    annotate("text",
      x = chart_width_mm / 2, y = x_title_y, label = x_title,
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black"
    ) +
    annotate("text",
      x = y_title_x, y = chart_height_mm / 2, label = y_title, angle = 90,
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black"
    ) +
    geom_point(
      data = gene_df, aes(x = x_mm, y = y_mm, size = point_size_mm, fill = I(fill_color), colour = I(stroke_color)),
      shape = 21, stroke = gene_df$stroke_width
    ) +
    scale_size_identity() +
    geom_text_repel(
      data = label_df,
      aes(x = x_mm, y = y_mm, label = gene_symbol),
      size = GENE_LABEL_PT * point_in_mm, family = FONT_FAMILY, fontface = "italic",
      # point.padding turned out not to be an effective lever here -- it
      # only reserves space around each point's (x, y) coordinate, not its
      # actual rendered bubble radius, so raising it (tried up to 4.5mm)
      # produced no visible change. nudge_y instead displaces every label's
      # *starting* position before repulsion runs, which reliably forces a
      # visible gap -- and therefore a callout line, since min.segment.length
      # is 0 -- for virtually every label rather than only the ones that
      # would otherwise collide with something.
      nudge_x = label_df$nudge_x_mm, nudge_y = label_df$nudge_y_mm, box.padding = unit(0.6, "mm"), point.padding = unit(0.5, "mm"),
      force = 3, force_pull = 0.5, max.overlaps = Inf,
      segment.size = 0.15, segment.color = "grey40", min.segment.length = 0
    ) +
    # ---- Legend ----
    geom_rect(
      data = legend_swatches, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill))
    ) +
    geom_point(
      data = legend_acmg_ring, inherit.aes = FALSE,
      aes(x = x, y = y, size = size_mm), shape = 21, fill = NA, colour = ACMG_STROKE_COLOR, stroke = 0.9
    ) +
    geom_text(
      data = legend_labels, inherit.aes = FALSE, aes(x = x, y = y, label = label),
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 0, vjust = 0.5
    ) +
    geom_text(
      data = legend_size_title, inherit.aes = FALSE, aes(x = x, y = y, label = label),
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 0, vjust = 1,
      lineheight = LEGEND_LINEHEIGHT
    ) +
    geom_point(
      data = legend_size_circles, inherit.aes = FALSE,
      aes(x = x, y = y, size = size_mm), shape = 21, fill = NA, colour = "black", stroke = 0.7
    ) +
    geom_text(
      data = legend_size_circles, inherit.aes = FALSE, aes(x = label_x, y = y, label = value_label),
      size = LABEL_PT * point_in_mm, family = FONT_FAMILY, colour = "black", hjust = 0
    ) +
    coord_cartesian(xlim = x_range, ylim = y_range, expand = FALSE, clip = "off") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

  list(plot = p, width_mm = diff(x_range), height_mm = diff(y_range))
}
