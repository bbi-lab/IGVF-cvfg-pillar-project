# Shared calibrated sunburst/donut chart logic (Figure 3 -- disease -> assay
# type -> model system -> gene).
#
# Renders at the *exact* final placed size (cairo_pdf, mm units, Arial),
# matching the convention in three_ring_donut.R/confusion_matrix_calibrated.R/
# metrics_bar_calibrated.R -- at 100% placement scale in the assembled
# figure, text renders at true point size, with no rescaling needed.
#
# Unlike three_ring_donut.R's three FIXED categorical rings, this taxonomy
# has a variable number of branches at every level (assay type, model
# system, and especially gene), so the ring geometry is built bottom-up from
# a single flattened, ordered sequence of gene leaves (each leaf = 1 angular
# unit) rather than three_ring_donut.R's per-parent cumsum-and-rescale
# approach. The disease level itself is not drawn as its own ring (it would
# just repeat the same information as the fill color) -- it's used only to
# color every ring and to build the center legend, leaving the vacated
# center hole for that legend.
#
# The text-measurement/fit primitives (text_width_mm, polar_xy,
# arc_length_mm, local_room_mm, best_horizontal_angle, needs_reposition,
# wrap_two_lines) are the same techniques three_ring_donut.R uses, duplicated
# here rather than sourced from it -- each calibrated-figure lib is meant to
# be self-contained so multiple can be sourced into the same document without
# constant-name collisions (see three_ring_donut.R's own header comment).
#
# Expects the caller to have already loaded: dplyr, ggplot2, ggforce,
# extrafont (+ loadfonts(device = 'all')), and grid.

point_in_mm <- 0.3527778
FONT_FAMILY <- "Arial"

# geom_text()'s `size` aesthetic is in mm, converted to points internally via
# `size * .pt` where ggplot2's .pt = 72.27/25.4 -- NOT the same 72-points-
# per-inch convention point_in_mm uses for grid::gpar()-based measurement.
# Divide an actual point size by this before passing it to geom_text's size
# (whether as a fixed parameter or an aes()-mapped column) to get the
# intended point size on the page; measurement functions below (text_width_mm
# et al, which use grid::gpar(fontsize=)) take the point size directly and
# need no such conversion.
GGPLOT_PT_PER_MM <- 72.27 / 25.4
pt_to_gg_size <- function(pt) pt / GGPLOT_PT_PER_MM

# Base label size; auto-fit tries LABEL_PT_MAX down to LABEL_PT_MIN (see
# fit_label() below) rather than a single fixed size, so a handful of
# short/roomy labels end up larger and a handful of long/cramped ones
# (e.g. "Immortalized human cells" in a narrow wedge) end up smaller.
LABEL_PT <- 6
LABEL_PT_MAX <- 7
LABEL_PT_MIN <- 5

# Last-resort floors, tried only after LABEL_PT_MIN fails at every size and
# wrap combination -- for the handful of labels ("Direct protein function")
# that don't fit their wedge at any normal size. LAST_RESORT_WRAP_PT_MIN
# still wraps to two lines; LAST_RESORT_PT_MIN is single-line only, for the
# single narrowest wedge where even a wrapped 4.5pt doesn't fit.
LAST_RESORT_WRAP_PT_MIN <- 4.5
LAST_RESORT_PT_MIN <- 3

# wrap_two_lines()'s balanced (minimize length-difference) split picks
# "Direct protein" / "function" for this label; overridden here to split
# after the first word instead, matching how this specific label is meant
# to be wrapped when it's the tightest fit in the chart.
MANUAL_WRAP_OVERRIDES <- c("Direct protein function" = "Direct|protein function")

DISEASE_COLORS <- c(
  "Cancer" = "#F39B7FFF",
  "Cardio" = "#3C5488FF",
  "Metabolic" = "#00A087FF",
  "Other Rare Disease" = "#B09C85FF"
)
DISEASE_DISPLAY_NAMES <- c(
  "Cancer" = "Cancer",
  "Cardio" = "Cardiovascular",
  "Metabolic" = "Metabolic",
  "Other Rare Disease" = "Other rare\ndiseases"
)

###------------------------------------------------------------
### Text measurement / fit primitives (see three_ring_donut.R for the
### original derivation of each of these)
###------------------------------------------------------------

text_width_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY) {
  gp <- gpar(fontfamily = family, fontsize = pt, fontface = "plain")
  vapply(
    label,
    function(l) convertWidth(grobWidth(textGrob(l, gp = gp)), "mm", valueOnly = TRUE),
    numeric(1)
  )
}

# angle=0 is 3 o'clock (0 degrees from horizontal, on the right), increasing
# counterclockwise -- the standard math convention, matching how the
# published figure starts (Cancer/BAP1 at 0 degrees, sweeping
# counterclockwise from there). `theta`/`start`/`end` throughout this file
# (build_sunburst_rings() etc.) are plain data-sequence positions in this
# same convention. geom_arc_bar() itself uses a different, fixed convention
# (0 = 12 o'clock, increasing clockwise) that can't be changed -- see
# to_ggforce_angle() below, used only at the three geom_arc_bar() call sites,
# for the one-time conversion between the two. Wedge *widths* (end - start)
# are identical under both conventions (the conversion is a reflection, not
# a rescale), so every fit/measurement calculation elsewhere in this file
# that only cares about angular width is unaffected.
polar_xy <- function(r, angle) {
  tibble(x = r * cos(angle), y = r * sin(angle))
}

# Converts a [start, end] pair from this file's own (3-o'clock, CCW+)
# convention to geom_arc_bar()'s (12-o'clock, CW+) convention. Note the
# swap: reflecting an increasing range still needs start < end afterward.
to_ggforce_angle <- function(start, end) list(start = pi / 2 - end, end = pi / 2 - start)

arc_length_mm <- function(r, theta, mm_per_unit) r * theta * mm_per_unit
chord_mm <- function(r, theta, mm_per_unit) 2 * r * sin(theta / 2) * mm_per_unit

# Angular half-gap between a wrapped label's two side-by-side lines (see
# render_labels()) -- a function of font size and radius only, so it must be
# checked against the wedge's own angular width (theta) wherever it's used:
# on a narrow enough wedge this can otherwise exceed half the wedge's own
# span, pushing a line's position past its own boundary into a neighboring
# wedge (fit_label()'s wrap branch rejects any pt where that would happen).
line_gap_rad <- function(pt, r, mm_per_unit) (pt * point_in_mm * 1.1) / (r * mm_per_unit)

# Every label sits at its own wedge's bisector and runs radially (its
# baseline follows the line from center through that angle), rather than
# staying horizontal -- see plot_sunburst_calibrated()'s radial_angle_deg().
# That means the constraining dimensions for whether a label fits are:
# radially, the ring's own depth (r1 - r0); tangentially, the wedge's chord
# width at r_mid (how much angular room the font's line-height needs without
# touching the previous/next wedge's own label).
wrap_two_lines <- function(label) {
  if (label %in% names(MANUAL_WRAP_OVERRIDES)) {
    return(strsplit(MANUAL_WRAP_OVERRIDES[[label]], "\\|")[[1]])
  }
  words <- strsplit(label, " ")[[1]]
  if (length(words) < 2) return(c(label, NA_character_))
  best <- c(label, NA_character_)
  best_diff <- Inf
  for (i in seq_len(length(words) - 1)) {
    l1 <- paste(words[1:i], collapse = " ")
    l2 <- paste(words[(i + 1):length(words)], collapse = " ")
    d <- abs(nchar(l1) - nchar(l2))
    if (d < best_diff) {
      best_diff <- d
      best <- c(l1, l2)
    }
  }
  best
}

text_color_for_fill <- function(hex) {
  vapply(hex, function(h) {
    if (is.na(h)) return("black")
    rgb <- grDevices::col2rgb(h) / 255
    luminance <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
    if (luminance < 0.5) "white" else "black"
  }, character(1), USE.NAMES = FALSE)
}

###------------------------------------------------------------
### Ring geometry -- built bottom-up from a single ordered gene-leaf sequence
###------------------------------------------------------------

# gene_rows: one row per (Disease, Assay.Type, Model_system, Gene) MAVE
# dataset membership (duplicates across datasets collapse via distinct()
# below). Diseases/assay types/model systems keep first-appearance order
# (matching the original plotly sunburst's distinct()-preserves-order
# behavior); genes sort alphabetically within their (Disease, Assay.Type,
# Model_system) group. Ordering is nested one level at a time (disease, then
# assay type within disease, then model system within that pair) rather than
# a single flat distinct(Disease, Assay.Type, Model_system) pass -- the flat
# version only guarantees each triple is unique, not that every (Disease,
# Assay.Type) pair's model-system children stay contiguous in the final
# order, which the "parent span = min/max of its children's" step below
# depends on. If the source sheet interleaves the same assay type's rows
# with a different assay type's rows in between (a real occurrence in
# Supplementary_Data_3.xlsx), a flat ordering silently produces a parent span
# that swallows the unrelated assay type's wedge too.
# Clockwise starting position for each disease's arc, matching the published
# figure (Cancer starts at 12 o'clock) rather than the curation sheet's own
# row order.
DISEASE_ORDER <- c("Cancer", "Cardio", "Metabolic", "Other Rare Disease")

build_sunburst_rings <- function(gene_rows) {
  disease_order <- tibble(Disease = DISEASE_ORDER, disease_order = seq_along(DISEASE_ORDER))

  assay_order <- gene_rows %>%
    distinct(Disease, Assay.Type) %>%
    left_join(disease_order, by = "Disease") %>%
    group_by(Disease) %>%
    mutate(assay_order = row_number()) %>%
    ungroup()

  model_order <- gene_rows %>%
    distinct(Disease, Assay.Type, Model_system) %>%
    left_join(assay_order, by = c("Disease", "Assay.Type")) %>%
    group_by(Disease, Assay.Type) %>%
    mutate(model_order = row_number()) %>%
    ungroup()

  gene_level <- gene_rows %>%
    distinct(Disease, Assay.Type, Model_system, Gene) %>%
    left_join(model_order, by = c("Disease", "Assay.Type", "Model_system")) %>%
    arrange(disease_order, assay_order, model_order, Gene) %>%
    select(-disease_order, -assay_order, -model_order) %>%
    mutate(
      end = cumsum(rep(1, n())) / n() * 2 * pi,
      start = lag(end, default = 0),
      fill = DISEASE_COLORS[Disease]
    )

  # Each parent's [start, end] is exactly the min/max of its children's --
  # valid because gene_level's rows are grouped contiguously by construction
  # (sorted by group_order, which walks each (Disease, Assay.Type,
  # Model_system) triple as one contiguous block).
  model_level <- gene_level %>%
    group_by(Disease, Assay.Type, Model_system) %>%
    summarise(start = min(start), end = max(end), n = n(), .groups = "drop") %>%
    mutate(fill = DISEASE_COLORS[Disease])

  assay_level <- model_level %>%
    group_by(Disease, Assay.Type) %>%
    summarise(start = min(start), end = max(end), n = sum(n), .groups = "drop") %>%
    mutate(fill = DISEASE_COLORS[Disease])

  disease_level <- assay_level %>%
    group_by(Disease) %>%
    summarise(start = min(start), end = max(end), n = sum(n), .groups = "drop") %>%
    mutate(fill = DISEASE_COLORS[Disease])

  list(gene = gene_level, model = model_level, assay = assay_level, disease = disease_level)
}

###------------------------------------------------------------
### Rotation
###------------------------------------------------------------

# Our polar_xy() angle (theta) already is ggplot2's geom_text(angle=)
# convention (degrees, counterclockwise from the positive x-axis), so a
# label rotated by theta (in degrees) runs exactly along the radius at
# theta, reading outward. Anywhere that would put the label upside down
# (baseline pointing into the left half-plane) gets flipped 180 degrees so
# every label stays readable left-to-right, never inverted -- the same
# purpose geomtextpath's `upright = TRUE` serves for the curved labels in
# three_ring_donut.R, done manually here since these labels are straight.
radial_angle_deg <- function(theta) {
  deg <- theta * 180 / pi
  deg <- ((deg + 180) %% 360) - 180
  ifelse(deg > 90 | deg < -90, deg + 180, deg)
}

###------------------------------------------------------------
### Label fitting -- try LABEL_PT_MAX down to LABEL_PT_MIN, single line then
### two lines set side by side (tangentially), before giving up and using
### the smallest size regardless of fit.
###------------------------------------------------------------

fit_label <- function(label_text, r0, r1, start, end, mm_per_unit,
                       italic = FALSE, radial_fudge = 0.92, tangential_fudge = 0.85,
                       allow_wrap = TRUE) {
  theta <- end - start
  r_mid <- (r0 + r1) / 2
  mid_angle <- (start + end) / 2
  avail_radial_mm <- (r1 - r0) * mm_per_unit * radial_fudge
  # Measured at r0, the wedge's narrowest point (wedges fan out with
  # radius), not r_mid -- a label spans nearly the ring's full depth, so its
  # tangential footprint has to clear the tightest point along that span,
  # not just the middle.
  avail_tangential_mm <- chord_mm(r0, theta, mm_per_unit) * tangential_fudge
  fontface <- if (italic) "italic" else "plain"

  text_len_mm <- function(txt, pt) {
    gp <- gpar(fontfamily = FONT_FAMILY, fontsize = pt, fontface = fontface)
    convertWidth(grobWidth(textGrob(txt, gp = gp)), "mm", valueOnly = TRUE)
  }
  # A line's thickness (its extent in the tangential direction once
  # rotated radially) is well approximated by its font size -- cap height
  # plus a little for ascenders/descenders.
  line_thickness_mm <- function(pt) pt * point_in_mm * 1.3

  fits <- function(txt, pt, n_lines = 1) {
    text_len_mm(txt, pt) < avail_radial_mm && (line_thickness_mm(pt) * n_lines) < avail_tangential_mm
  }

  for (pt in seq(LABEL_PT_MAX, LABEL_PT_MIN, by = -0.5)) {
    if (fits(label_text, pt)) {
      return(list(lines = label_text, pt = pt, angle = mid_angle, r = r_mid, wrapped = FALSE))
    }
  }
  if (allow_wrap) {
    wrapped <- wrap_two_lines(label_text)
    if (!is.na(wrapped[2])) {
      for (pt in seq(LABEL_PT, LABEL_PT_MIN, by = -0.5)) {
        # The two lines' angular half-gap (each offset by half_gap from the
        # bisector) must stay well inside the wedge's own [start, end] --
        # otherwise they'd land in a neighboring wedge instead of this one.
        half_gap <- line_gap_rad(pt, r_mid, mm_per_unit) / 2
        gap_fits <- half_gap < theta / 2 * 0.8
        if (gap_fits &&
          max(text_len_mm(wrapped[1], pt), text_len_mm(wrapped[2], pt)) < avail_radial_mm &&
          line_thickness_mm(pt) * 2 < avail_tangential_mm) {
          return(list(lines = wrapped, pt = pt, angle = mid_angle, r = r_mid, wrapped = TRUE))
        }
      }
    }
  }
  # Last resort: nothing fit at any normal size (wrapped or not) down to
  # LABEL_PT_MIN, using the same safety margins every other label gets. Try
  # again down to LAST_RESORT_WRAP_PT_MIN/LAST_RESORT_PT_MIN with those
  # margins relaxed toward the actual measured text metrics (still real
  # grid::convertWidth()-based measurements, just less padded) -- appropriate
  # here since these are already the few labels we're deliberately squeezing
  # as tightly as will still read, not the general case -- before finally
  # giving up (render_labels() omits a still-no_fit label rather than force
  # it on).
  fits_loose <- function(txt, pt, n_lines = 1) {
    text_len_mm(txt, pt) < avail_radial_mm / 0.92 * 0.98 &&
      (pt * point_in_mm * 1.1 * n_lines) < avail_tangential_mm / tangential_fudge * 0.95
  }
  if (allow_wrap) {
    wrapped <- wrap_two_lines(label_text)
    if (!is.na(wrapped[2])) {
      for (pt in seq(LABEL_PT_MIN - 0.5, LAST_RESORT_PT_MIN, by = -0.5)) {
        half_gap <- line_gap_rad(pt, r_mid, mm_per_unit) / 2
        if (half_gap < theta / 2 * 0.9 &&
          max(text_len_mm(wrapped[1], pt), text_len_mm(wrapped[2], pt)) < avail_radial_mm / 0.92 * 0.98 &&
          (pt * point_in_mm * 1.1 * 2) < avail_tangential_mm / tangential_fudge * 0.95) {
          return(list(lines = wrapped, pt = pt, angle = mid_angle, r = r_mid, wrapped = TRUE))
        }
      }
    }
  }
  for (pt in seq(LAST_RESORT_WRAP_PT_MIN - 0.5, LAST_RESORT_PT_MIN, by = -0.5)) {
    if (fits_loose(label_text, pt)) {
      return(list(lines = label_text, pt = pt, angle = mid_angle, r = r_mid, wrapped = FALSE))
    }
  }
  # Truly nothing fits -- caller omits this label.
  list(lines = label_text, pt = LAST_RESORT_PT_MIN, angle = mid_angle, r = r_mid, wrapped = FALSE, no_fit = TRUE)
}

###------------------------------------------------------------
### Full plot
###------------------------------------------------------------

# diameter_mm defaults to 3.6in (91.44mm), the published figure's size.
# inner_hole_frac sizes the blank center (as a fraction of the outer radius)
# left for the disease-color legend -- kept small since the legend itself
# only needs a handful of short lines of text. ring_depth_weights splits the
# remaining radius unevenly rather than into equal thirds: since text runs
# radially, a ring's available label *length* is its own depth in mm, and
# assay/model-system labels ("Direct protein function", "immortalized human
# cells") run far longer than most gene symbols, so they get a larger share.
# Both defaults were tuned empirically against this project's actual curation
# data (minimizing labels that don't fit even wrapped at LABEL_PT_MIN, so get
# dropped by render_labels() below) -- the fit is a genuine trade-off between
# rings/hole size with real cliffs nearby (a change of ~0.01 can flip several
# labels from fitting to not), not a formula. Re-tune by hand if the
# underlying curation data changes enough to shift which wedges are
# narrowest.
plot_sunburst_calibrated <- function(gene_rows, diameter_mm = 91.44, inner_hole_frac = 0.25,
                                      ring_depth_weights = c(assay = 0.333, model = 0.40, gene = 0.267)) {
  outer_r_max <- 1
  mm_per_unit <- diameter_mm / (outer_r_max * 2)

  rings <- build_sunburst_rings(gene_rows)

  available_span <- outer_r_max - inner_hole_frac
  ring_depth_weights <- unname(ring_depth_weights / sum(ring_depth_weights))
  assay_r1 <- inner_hole_frac + available_span * ring_depth_weights[1]
  model_r1 <- assay_r1 + available_span * ring_depth_weights[2]
  radii <- list(
    assay = c(r0 = inner_hole_frac, r1 = assay_r1),
    model = c(r0 = assay_r1, r1 = model_r1),
    gene  = c(r0 = model_r1, r1 = outer_r_max)
  )

  assay_ring <- rings$assay %>% mutate(r0 = radii$assay["r0"], r1 = radii$assay["r1"])
  model_ring <- rings$model %>% mutate(r0 = radii$model["r0"], r1 = radii$model["r1"])
  gene_ring  <- rings$gene  %>% mutate(r0 = radii$gene["r0"],  r1 = radii$gene["r1"])

  # Model-system values in Supplementary_Data_3.xlsx are all lowercase
  # ("immortalized human cells", "yeast", ...); assay types are already
  # sentence case, so this is a no-op for them. Never applied to gene
  # symbols, which must stay in their actual HGNC casing (e.g. "BAP1").
  to_sentence_case <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

  fit_ring_labels <- function(ring_df, label_col, italic = FALSE, sentence_case = TRUE) {
    ring_df$.label_text <- ring_df[[label_col]]
    if (sentence_case) ring_df$.label_text <- to_sentence_case(ring_df$.label_text)
    fits <- Map(
      fit_label,
      ring_df$.label_text, ring_df$r0, ring_df$r1, ring_df$start, ring_df$end,
      MoreArgs = list(mm_per_unit = mm_per_unit, italic = italic)
    )
    ring_df$pt <- vapply(fits, `[[`, numeric(1), "pt")
    ring_df$gg_size <- pt_to_gg_size(ring_df$pt)
    ring_df$angle <- vapply(fits, `[[`, numeric(1), "angle")
    ring_df$r_label <- vapply(fits, `[[`, numeric(1), "r")
    ring_df$wrapped <- vapply(fits, `[[`, logical(1), "wrapped")
    ring_df$no_fit <- vapply(fits, function(f) isTRUE(f$no_fit), logical(1))
    ring_df$line1 <- vapply(fits, function(f) f$lines[1], character(1))
    ring_df$line2 <- vapply(fits, function(f) if (length(f$lines) > 1) f$lines[2] else NA_character_, character(1))
    ring_df$text_color <- text_color_for_fill(ring_df$fill)
    ring_df
  }

  assay_ring <- fit_ring_labels(assay_ring, "Assay.Type")
  model_ring <- fit_ring_labels(model_ring, "Model_system")
  gene_ring  <- fit_ring_labels(gene_ring, "Gene", italic = TRUE, sentence_case = FALSE)

  # A wrapped label's two lines sit side by side *tangentially* (offset in
  # angle, at the same radius) rather than stacked radially -- each line
  # still runs radially itself, just shifted along the arc from the other.
  # fit_label() already only selects `wrapped = TRUE` for a row when this
  # exact gap fits inside its own wedge, so no further clamping is needed
  # here -- see line_gap_rad()'s own definition for why it needs mm_per_unit.
  render_labels <- function(ring_df) {
    # Even a label flagged no_fit (nothing fit at any size/wrap combination
    # down to the last-resort floors) still renders, at its last-attempted
    # size/position -- may run past its own wedge's boundary into a
    # neighbor's, a deliberate trade-off for this handful of cases over
    # leaving them blank.
    single <- ring_df %>% filter(!wrapped)
    wrapped <- ring_df %>% filter(wrapped)

    out <- list()
    if (nrow(single) > 0) {
      xy <- polar_xy(single$r_label, single$angle)
      out$single <- bind_cols(single, xy) %>% mutate(rot_deg = radial_angle_deg(angle))
    }
    if (nrow(wrapped) > 0) {
      half_gap <- line_gap_rad(wrapped$pt, wrapped$r_label, mm_per_unit) / 2
      angle_a <- wrapped$angle - half_gap
      angle_b <- wrapped$angle + half_gap
      # Whichever of the two candidate positions sits higher on the page
      # (larger y = r*sin(angle)) is where line1 (the first word/phrase)
      # goes, so the two lines always read top-to-bottom in on-page order --
      # otherwise which offset direction is "first" is arbitrary and the
      # apparent line order flips depending on which side of the circle the
      # wedge is on.
      a_on_top <- sin(angle_a) >= sin(angle_b)
      line1_angle <- ifelse(a_on_top, angle_a, angle_b)
      line2_angle <- ifelse(a_on_top, angle_b, angle_a)
      xy1 <- polar_xy(wrapped$r_label, line1_angle)
      xy2 <- polar_xy(wrapped$r_label, line2_angle)
      out$line1 <- bind_cols(wrapped %>% mutate(.label = line1), xy1) %>%
        mutate(rot_deg = radial_angle_deg(line1_angle))
      out$line2 <- bind_cols(wrapped %>% mutate(.label = line2), xy2) %>%
        mutate(rot_deg = radial_angle_deg(line2_angle))
    }
    out
  }

  assay_labels <- render_labels(assay_ring)
  model_labels <- render_labels(model_ring)
  gene_labels <- render_labels(gene_ring)

  # geom_arc_bar() needs its own (12-o'clock, clockwise) angle convention --
  # see to_ggforce_angle()'s definition up top.
  to_arc_bar_data <- function(ring_df) {
    ggforce_angles <- to_ggforce_angle(ring_df$start, ring_df$end)
    ring_df %>% mutate(arc_start = ggforce_angles$start, arc_end = ggforce_angles$end)
  }
  assay_ring_arcs <- to_arc_bar_data(assay_ring)
  model_ring_arcs <- to_arc_bar_data(model_ring)
  gene_ring_arcs <- to_arc_bar_data(gene_ring)

  p <- ggplot() +
    geom_arc_bar(
      data = assay_ring_arcs,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = arc_start, end = arc_end, fill = I(fill)),
      linewidth = 0.5 * point_in_mm, color = "white"
    ) +
    geom_arc_bar(
      data = model_ring_arcs,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = arc_start, end = arc_end, fill = I(fill)),
      linewidth = 0.5 * point_in_mm, color = "white"
    ) +
    geom_arc_bar(
      data = gene_ring_arcs,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = arc_start, end = arc_end, fill = I(fill)),
      linewidth = 0.5 * point_in_mm, color = "white"
    )

  add_label_layers <- function(p, labels, italic = FALSE) {
    fontface <- if (italic) "italic" else "plain"
    if (!is.null(labels$single) && nrow(labels$single) > 0) {
      p <- p + geom_text(
        data = labels$single,
        aes(x = x, y = y, label = .label_text, size = gg_size, angle = rot_deg, color = I(text_color)),
        family = FONT_FAMILY, fontface = fontface
      )
    }
    if (!is.null(labels$line1) && nrow(labels$line1) > 0) {
      p <- p +
        geom_text(
          data = labels$line1,
          aes(x = x, y = y, label = .label, size = gg_size, angle = rot_deg, color = I(text_color)),
          family = FONT_FAMILY, fontface = fontface
        ) +
        geom_text(
          data = labels$line2,
          aes(x = x, y = y, label = .label, size = gg_size, angle = rot_deg, color = I(text_color)),
          family = FONT_FAMILY, fontface = fontface
        )
    }
    p
  }

  p <- p %>%
    add_label_layers(assay_labels) %>%
    add_label_layers(model_labels) %>%
    add_label_layers(gene_labels, italic = TRUE)

  # ---- center legend (disease colors) ----
  diseases <- rings$disease$Disease
  n_legend <- length(diseases)
  legend_labels <- DISEASE_DISPLAY_NAMES[diseases]
  swatch_pt <- LABEL_PT
  line_h_mm <- text_width_mm("Mg", swatch_pt) # ~line height proxy, matches metrics_bar_calibrated.R's convention
  line_h_mm <- max(line_h_mm, swatch_pt * point_in_mm * 1.4)
  # Swatches are spaced evenly regardless of label line count -- a wrapped
  # 2-line entry ("Other rare\ndiseases") uses LEGEND_LINEHEIGHT (a tight
  # multi-line spacing) to fit inside one uniform row rather than claiming
  # extra row height for itself.
  LEGEND_LINEHEIGHT <- 0.8
  row_h_mm <- line_h_mm * 1.35
  total_h_mm <- row_h_mm * n_legend
  swatch_size_mm <- swatch_pt * point_in_mm * 1.2
  legend_x0_mm <- -(inner_hole_frac * outer_r_max * mm_per_unit) * 0.6

  legend_df <- tibble(
    Disease = diseases,
    label = legend_labels,
    fill = DISEASE_COLORS[diseases],
    y_mm = total_h_mm / 2 - (row_h_mm * (seq_len(n_legend) - 0.5))
  ) %>%
    mutate(
      y = y_mm / mm_per_unit,
      swatch_x = legend_x0_mm / mm_per_unit,
      # Close to the swatch, not the more generous gap other calibrated
      # figures' legends use for a full-width legend -- this one has much
      # less horizontal room to work with, sitting inside the sunburst's
      # own center hole.
      label_x = (legend_x0_mm + swatch_size_mm * 1.15) / mm_per_unit
    )

  p <- p +
    geom_tile(
      data = legend_df,
      aes(x = swatch_x, y = y, fill = I(fill)),
      width = swatch_size_mm / mm_per_unit, height = swatch_size_mm / mm_per_unit
    ) +
    geom_text(
      data = legend_df,
      aes(x = label_x, y = y, label = label),
      hjust = 0, size = pt_to_gg_size(LABEL_PT), family = FONT_FAMILY, color = "black",
      lineheight = LEGEND_LINEHEIGHT
    )

  p <- p +
    scale_size_identity() +
    coord_fixed(xlim = c(-outer_r_max, outer_r_max), ylim = c(-outer_r_max, outer_r_max), expand = FALSE, clip = "off") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

  list(plot = p, width_mm = diameter_mm, height_mm = diameter_mm)
}
