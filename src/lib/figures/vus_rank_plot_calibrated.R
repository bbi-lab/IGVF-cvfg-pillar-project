# Calibrated rank-based dot plot (Extended Data Figure 8): per-gene %VUS
# remaining (0-5 combined-evidence points, i.e. still unresolved) after
# REVEL/AlphaMissense/MutPred2 gene-specific reclassification, one point per
# gene faceted by predictor (Overall/REVEL/AM/MP2), genes ranked by the
# Overall column.
#
# Unlike this session's other calibrated figures, this one is a standard
# ggplot2 construction (facet_grid + geom_point/geom_segment, theme_classic)
# -- nothing here needs a hand-rolled grid layout, since ggplot2's own gtable
# already sizes axis titles/tick labels/facet strips exactly from their
# theme() font sizes and fits them into whatever canvas ggsave() is given.
# The only real calibration problem is fitting ~40 gene rows legibly into a
# very short total height (about 4in here): text_height_mm() picks the
# largest y-axis label size that still lets all of them fit, rather than
# guessing a fixed point size that might overflow (too large) or waste
# space (too small) if the gene list's length changes. theme_classic()'s
# *default* spacing (strip padding, the gap above the x-axis title, etc.) is
# all sized off an 11pt base font, though, which wastes real room here at
# 6-7pt -- every element below explicitly overrides its margin instead of
# inheriting that default, freeing up that space for the panel (taller gene
# rows) rather than unused padding.
#
# This file's companion, vus_points_heatmap_calibrated.R, shares the same
# 40 genes and the same total_height_mm -- see its own header comment for
# how the two keep their y-axis tick marks lined up despite each ranking
# genes in a different order.
#
# Expects the caller to have already loaded: dplyr, ggplot2, extrafont (+
# loadfonts(device = 'all')), and grid.

point_in_mm <- 0.3527778
FONT_FAMILY <- "Arial"
LABEL_PT_MAX <- 7
LABEL_PT_MIN <- 4
# Explicit, small gaps -- replacing theme_classic()'s 11pt-base-font-sized
# defaults -- between chrome elements stacked below the panel. TITLE_GAP_MM
# in particular *is* "the distance between the x-axis and its title": kept
# as its own named constant (rather than folded into a font-derived margin)
# so it stays exactly this wide no matter how the title's own size changes.
STRIP_PAD_MM <- 0.6
TICK_LABEL_GAP_MM <- 0.8
TITLE_GAP_MM <- 1.3
BOTTOM_PAD_MM <- 0.5
TOP_PAD_MM <- 0.5
TITLE_LINEHEIGHT <- 0.9

text_height_mm <- function(label, pt, family = FONT_FAMILY, fontface = "plain", lineheight = 1) {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = fontface, lineheight = lineheight)
  vapply(label, function(l) convertHeight(grobHeight(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

text_width_mm <- function(label, pt, family = FONT_FAMILY, fontface = "plain") {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = fontface)
  vapply(label, function(l) convertWidth(grobWidth(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

# Greedy word-wrap of `text` to fit within max_width_mm at the given point
# size, returning a single "\n"-joined string -- the x-axis title here is
# long enough (in "% VUS remaining (within experimentally assessed
# regions)") that it doesn't fit on one line within this figure's narrow
# total width, and ggplot2 doesn't wrap axis titles on its own.
wrap_text_to_width <- function(text, pt, max_width_mm, family = FONT_FAMILY) {
  words <- strsplit(text, " ", fixed = TRUE)[[1]]
  lines <- character(0)
  current <- ""
  for (word in words) {
    candidate <- if (current == "") word else paste(current, word)
    if (current != "" && text_width_mm(candidate, pt, family) > max_width_mm) {
      lines <- c(lines, current)
      current <- word
    } else {
      current <- candidate
    }
  }
  lines <- c(lines, current)
  paste(lines, collapse = "\n")
}

# Chrome above the panel: top padding + the facet strip row.
top_chrome_mm <- function(strip_pt) {
  TOP_PAD_MM + text_height_mm("Mg", strip_pt, fontface = "bold") + 2 * STRIP_PAD_MM
}

# Chrome below the panel: axis ticks, the (angled) x tick label block, the
# gap this file preserves between the axis and its title (TITLE_GAP_MM),
# the title itself (measured as actually rendered -- title_text may be
# multi-line, "\n"-joined, and its height with TITLE_LINEHEIGHT applied
# depends on that line count in a way a naive single-line-height*n_lines
# estimate gets wrong, since lineheight only affects *inter*-line spacing),
# and bottom padding.
bottom_chrome_mm <- function(x_tick_pt, title_pt, title_text = "Mg") {
  tick_len_mm <- 0.8
  # 45-degree tick labels: their footprint below the axis is dominated by
  # the widest label's *width*, not its height, but a fixed allowance sized
  # off font metrics (rather than a magic constant) keeps this in step with
  # x_tick_pt if that's ever changed.
  x_tick_h <- text_height_mm("100", x_tick_pt) * 2.2 + TICK_LABEL_GAP_MM
  title_h <- text_height_mm(title_text, title_pt, lineheight = TITLE_LINEHEIGHT)
  tick_len_mm + x_tick_h + TITLE_GAP_MM + title_h + BOTTOM_PAD_MM
}

# plot_df: one row per (Gene, Predictor) with pct_vus_remaining (x) --
# Predictor already an ordered factor (Overall, REVEL, AM, MP2) and Gene
# already an ordered factor (ranked by the Overall column), matching
# Extended_data_figures.Rmd's own data prep.
#
# target_bottom_chrome_mm: when given (by the companion heatmap, to match
# this plot's own natural bottom chrome -- see that file), pads
# plot.margin's bottom by however much more this plot's *own* bottom chrome
# would otherwise need beyond that target, so the two plots' bottom margins
# -- and so their panels' bottom edges and all 40 tick marks -- line up
# exactly. Returns bottom_chrome_mm/top_chrome_mm so the companion heatmap
# (rendered second) can request this plot's own values as its target.
plot_vus_rank_calibrated <- function(plot_df, total_width_mm = 3.05 * 25.4, total_height_mm = 4.1 * 25.4,
                                      x_tick_pt = 6, title_pt = 7, strip_pt = 7,
                                      target_bottom_chrome_mm = NULL) {
  measure_dev_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(measure_dev_file)
  on.exit({
    grDevices::dev.off()
    unlink(measure_dev_file)
  }, add = TRUE)

  genes <- levels(plot_df$Gene)
  n_genes <- length(genes)

  # The x-axis title is centered under the *panel* area only, not the full
  # canvas -- it has to share the canvas with the y-axis gene labels, so its
  # wrap width depends on how wide those end up (gene_pt), while gene_pt's
  # own row-height budget depends on how many lines the title wraps to. Two
  # passes resolve that: an initial gene_pt from a 1-line assumption, then
  # the real wrap (now that the y-axis label width is known) determines the
  # actual line count for a final, corrected gene_pt.
  x_title <- "% VUS remaining (within experimentally assessed regions)"
  choose_gene_pt <- function(bottom_mm) {
    available_row_height_mm <- (total_height_mm - top_chrome_mm(strip_pt) - bottom_mm) / n_genes
    for (pt in seq(LABEL_PT_MAX, LABEL_PT_MIN, by = -0.5)) {
      if (text_height_mm("Ag", pt, fontface = "italic") <= available_row_height_mm) {
        return(pt)
      }
    }
    LABEL_PT_MIN
  }

  gene_pt <- choose_gene_pt(bottom_chrome_mm(x_tick_pt, title_pt, title_text = x_title))
  y_axis_label_w_mm <- max(text_width_mm(genes, gene_pt, fontface = "italic"))
  title_budget_mm <- total_width_mm - y_axis_label_w_mm - 4
  x_title_wrapped <- wrap_text_to_width(x_title, title_pt, title_budget_mm)

  own_bottom_chrome_mm <- bottom_chrome_mm(x_tick_pt, title_pt, x_title_wrapped)
  bottom_extra_mm <- if (is.null(target_bottom_chrome_mm)) {
    0
  } else {
    max(0, target_bottom_chrome_mm - own_bottom_chrome_mm)
  }
  gene_pt <- choose_gene_pt(own_bottom_chrome_mm + bottom_extra_mm)

  x_min <- floor(min(plot_df$pct_vus_remaining, na.rm = TRUE) / 5) * 5
  x_max <- ceiling(max(plot_df$pct_vus_remaining, na.rm = TRUE) / 5) * 5

  p <- ggplot(plot_df, aes(x = pct_vus_remaining, y = Gene)) +
    geom_segment(
      aes(x = x_min, xend = x_max, yend = Gene),
      linewidth = 0.35, color = "grey90"
    ) +
    geom_point(size = 1) +
    facet_grid(. ~ Predictor) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = seq(x_min, x_max, by = 25),
      expand = c(0.02, 0.02)
    ) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    # geom_point()'s rendered radius extends a little past its exact data
    # coordinate; with panels this narrow, the small default expansion
    # doesn't leave enough physical room for that at x = 0/100, so a point
    # right at either end was getting clipped to a half-circle by the panel
    # boundary (clip = "on" by default). Turning panel clipping off lets the
    # *point* draw its full circle just past the axis; nothing else in this
    # plot extends past the panel, so nothing else is affected.
    coord_cartesian(clip = "off") +
    labs(x = x_title_wrapped, y = NULL) +
    theme_classic(base_size = LABEL_PT_MAX, base_family = FONT_FAMILY) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold", size = strip_pt, family = FONT_FAMILY,
        margin = margin(STRIP_PAD_MM, 0, STRIP_PAD_MM, 0, unit = "mm")
      ),
      axis.text.y = element_text(size = gene_pt, family = FONT_FAMILY, face = "italic"),
      axis.text.x = element_text(
        angle = 45, hjust = 1, vjust = 1, size = x_tick_pt, family = FONT_FAMILY,
        margin = margin(t = TICK_LABEL_GAP_MM, unit = "mm")
      ),
      axis.title.x = element_text(
        size = title_pt, family = FONT_FAMILY, lineheight = TITLE_LINEHEIGHT,
        margin = margin(t = TITLE_GAP_MM, unit = "mm")
      ),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(0.8, "mm"),
      panel.spacing = unit(1, "mm"),
      # bottom_extra_mm (only ever nonzero if a *smaller* target_bottom_chrome_mm
      # than this plot's own natural bottom chrome were requested, which
      # doesn't happen in practice -- this plot has the taller of the two
      # bottom chromes) pads below the title, not the axis-to-title gap
      # itself, so TITLE_GAP_MM stays exactly what it says regardless.
      plot.margin = margin(TOP_PAD_MM, 1, BOTTOM_PAD_MM + bottom_extra_mm, 1, unit = "mm")
    )

  list(
    plot = p, width_mm = total_width_mm, height_mm = total_height_mm,
    top_chrome_mm = top_chrome_mm(strip_pt),
    bottom_chrome_mm = own_bottom_chrome_mm + bottom_extra_mm
  )
}
