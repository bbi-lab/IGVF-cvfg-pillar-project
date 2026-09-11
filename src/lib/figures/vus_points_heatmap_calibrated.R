# Calibrated heatmap (Extended Data Figure 8): per-gene distribution of
# combined-evidence points among ClinVar VUS, faceted by predictor (REVEL/
# AM/MP2), fill = proportion of that gene's VUS at each point value.
#
# Like vus_rank_plot_calibrated.R (the companion rank-based plot this
# shares its VUS data prep with), this is a standard ggplot2 construction
# (facet_grid + geom_tile, theme_classic) -- ggplot2's own gtable already
# sizes axis titles/tick labels/facet strips/legend from their theme() font
# sizes and fits them into whatever canvas ggsave() is given, so there's no
# hand-rolled grid layout here. The only real calibration problem is the
# same one as the rank plot: fitting ~40 gene rows legibly into a short
# total height, solved the same way (text_height_mm() picks the largest
# y-axis label size that still lets all of them fit). Like the rank plot,
# every element below explicitly overrides its margin instead of inheriting
# theme_classic()'s 11pt-base-font-sized default, to free that space up for
# the panel instead.
#
# Both this figure and the rank plot share the same 40 genes (albeit in a
# different order -- see below) and the same total_height_mm, and their
# y-axis tick marks are meant to line up despite that different order: pass
# this function's own bottom_chrome_mm (returned below) as the rank plot's
# target_bottom_chrome_mm -- or vice versa, whichever is rendered first --
# so both plots' panels end up the same height, with 40 evenly-spaced rows
# landing at the same y position in both.
#
# Expects the caller to have already loaded: dplyr, ggplot2, extrafont (+
# loadfonts(device = 'all')), and grid.

point_in_mm <- 0.3527778
FONT_FAMILY <- "Arial"
LABEL_PT_MAX <- 7
LABEL_PT_MIN <- 4
# Matches vus_rank_plot_calibrated.R's own constants of the same name --
# see that file for why these are named constants rather than left as
# theme_classic()'s (11pt-base-font-sized) defaults.
STRIP_PAD_MM <- 0.6
TICK_LABEL_GAP_MM <- 0.8
TITLE_GAP_MM <- 1.3
BOTTOM_PAD_MM <- 0.5
TOP_PAD_MM <- 0.5

text_height_mm <- function(label, pt, family = FONT_FAMILY, fontface = "plain") {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = fontface)
  vapply(label, function(l) convertHeight(grobHeight(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

# Chrome above the panel: top padding + the facet strip row -- same formula
# as vus_rank_plot_calibrated.R's, so the two plots' panels start at the
# same y position without needing to explicitly coordinate a target (as the
# bottom chrome does below).
top_chrome_mm <- function(strip_pt) {
  TOP_PAD_MM + text_height_mm("Mg", strip_pt, fontface = "bold") + 2 * STRIP_PAD_MM
}

# Chrome below the panel: axis ticks, the (angled) x tick label block, the
# axis-to-title gap, the title (measured as actually rendered -- see
# vus_rank_plot_calibrated.R's own bottom_chrome_mm(), whose 2-line wrapped
# title makes this matter more there), and bottom padding.
bottom_chrome_mm <- function(x_tick_pt, title_pt, title_text = "Total points") {
  tick_len_mm <- 0.8
  x_tick_h <- text_height_mm("100", x_tick_pt) * 2.2 + TICK_LABEL_GAP_MM
  title_h <- text_height_mm(title_text, title_pt)
  tick_len_mm + x_tick_h + TITLE_GAP_MM + title_h + BOTTOM_PAD_MM
}

# heat_df_all: one row per (Gene, Points, Predictor) with n/prop -- Predictor
# already an ordered factor (REVEL, AM, MP2), matching
# Extended_data_figures.Rmd's own data prep. Gene order isn't set explicitly
# there (a plain character column), so it renders in default ascending
# alphabetical order (bottom-to-top) -- unlike the companion rank plot's
# Overall-ranked order; kept as-is here to match the existing figure.
#
# target_bottom_chrome_mm: see vus_rank_plot_calibrated.R's own parameter of
# the same name -- this plot's "Total points" title is much shorter than
# the rank plot's, so it's normally the one padded to match (rather than
# the other way around). Returns bottom_chrome_mm/top_chrome_mm so the
# companion rank plot (if rendered second) can request this plot's own
# values as its target instead.
plot_vus_heatmap_calibrated <- function(heat_df_all, total_width_mm = 3.5 * 25.4, total_height_mm = 4.1 * 25.4,
                                         x_tick_pt = 6, title_pt = 7, strip_pt = 7, legend_pt = 6,
                                         target_bottom_chrome_mm = NULL) {
  measure_dev_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(measure_dev_file)
  on.exit({
    grDevices::dev.off()
    unlink(measure_dev_file)
  }, add = TRUE)

  genes <- unique(heat_df_all$Gene)
  n_genes <- length(genes)

  own_bottom_chrome_mm <- bottom_chrome_mm(x_tick_pt, title_pt)
  bottom_extra_mm <- if (is.null(target_bottom_chrome_mm)) {
    0
  } else {
    max(0, target_bottom_chrome_mm - own_bottom_chrome_mm)
  }

  available_row_height_mm <- (total_height_mm - top_chrome_mm(strip_pt) - own_bottom_chrome_mm - bottom_extra_mm) / n_genes
  gene_pt <- LABEL_PT_MIN
  for (pt in seq(LABEL_PT_MAX, LABEL_PT_MIN, by = -0.5)) {
    if (text_height_mm("Ag", pt, fontface = "italic") <= available_row_height_mm) {
      gene_pt <- pt
      break
    }
  }

  p <- ggplot(heat_df_all, aes(x = Points, y = Gene, fill = prop)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_grid(. ~ Predictor) +
    scale_fill_viridis_c(name = "Proportion") +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(x = "Total points", y = "Gene") +
    theme_classic(base_size = LABEL_PT_MAX, base_family = FONT_FAMILY) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, size = x_tick_pt, family = FONT_FAMILY,
        margin = margin(t = TICK_LABEL_GAP_MM, unit = "mm")
      ),
      axis.text.y = element_text(size = gene_pt, family = FONT_FAMILY, face = "italic"),
      axis.title.x = element_text(size = title_pt, family = FONT_FAMILY, margin = margin(t = TITLE_GAP_MM, unit = "mm")),
      axis.title.y = element_text(size = title_pt, family = FONT_FAMILY),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold", size = strip_pt, family = FONT_FAMILY,
        margin = margin(STRIP_PAD_MM, 0, STRIP_PAD_MM, 0, unit = "mm")
      ),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(0.8, "mm"),
      panel.spacing = unit(1, "mm"),
      legend.title = element_text(size = legend_pt, family = FONT_FAMILY),
      legend.text = element_text(size = legend_pt, family = FONT_FAMILY),
      legend.key.width = unit(3, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, -3, unit = "mm"),
      # bottom_extra_mm pads below the title, not the axis-to-title gap
      # itself, so TITLE_GAP_MM stays exactly what it says regardless.
      plot.margin = margin(TOP_PAD_MM, 1, BOTTOM_PAD_MM + bottom_extra_mm, 1, unit = "mm")
    )

  list(
    plot = p, width_mm = total_width_mm, height_mm = total_height_mm,
    top_chrome_mm = top_chrome_mm(strip_pt),
    bottom_chrome_mm = own_bottom_chrome_mm + bottom_extra_mm
  )
}
