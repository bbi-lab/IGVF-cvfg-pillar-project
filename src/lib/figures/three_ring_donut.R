# Shared three-ring donut chart logic (Fig 6c, Extended Data Fig 6/9).
#
# Source this from a donut-chart chunk after `colors_custom` is defined and
# before calling any of the make_three_ring_plot*() functions -- they look
# up `colors_custom` from the caller's environment at call time (matching
# the pre-existing make_three_ring_plot()/make_three_ring_plot_no_labels()
# convention), so any document that defines its own colors_custom can use
# these as-is. make_three_ring_plot_calibrated()'s `ring_diameter_mm`
# argument defaults to a `RING_DIAMETER_MM` the caller should define before
# invoking it (see Figure5_6.Rmd/Extended_data_figures.Rmd for the pattern).
#
# Expects the caller to have already loaded: dplyr, tidyr, ggplot2, ggforce,
# geomtextpath, extrafont (+ loadfonts(device = 'all')), and grid.
#
# df_subset passed to any of these functions is expected to be a subset of a
# `plot_summary_2`-shaped data frame: one row per (Total_Points_REVEL,
# point_source, conflict_status) combination, with columns `n` (count) and
# `total_variants` (the Total_Points_REVEL bin's total), matching the
# `cut()`-with-`complete()` pipeline in each calling document's own data-prep
# chunk (bin breaks/labels must match bin_display_labels below).

point_in_mm <- 0.3527778
LABEL_PT <- 5.5
LABEL_MM <- LABEL_PT * point_in_mm
# Smaller size for two-line-wrapped inner-ring labels, matching the original
# hand-placed figure's convention of shrinking text for this specific case.
WRAP_LABEL_PT <- 4.2
WRAP_LABEL_MM <- WRAP_LABEL_PT * point_in_mm
FONT_FAMILY <- "Arial"

# Count labels (middle/outer ring) get comma thousands separators; inner-ring
# bin names (e.g. "6 to 9") pass through unchanged since they aren't numeric.
format_count <- function(x) {
  if (is.numeric(x)) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  else as.character(x)
}

# Display text for inner-ring bin labels -- reads closer-to-zero value first,
# regardless of how the underlying Total_Points_REVEL bin factor is ordered
# (which also drives colors_custom's keys, left as-is). Uses the Unicode
# minus sign (U+2212), not a hyphen -- geomtextpath's glyph-curving silently
# drops a plain ASCII hyphen-minus regardless of font (reproduced with both
# Arial and sans), so a hyphen here would render as "1 to  3".
bin_display_labels <- c(
  "< -12" = "< −12",
  "-12 to -6" = "−6 to −12",
  "-5 to -4" = "−4 to −5",
  "-3 to -1" = "−1 to −3",
  "0" = "0",
  "1 to 3" = "1 to 3",
  "4 to 5" = "4 to 5",
  "6 to 9" = "6 to 9",
  "10 to 12" = "10 to 12",
  "> 12" = "> 12"
)

###------------------------------------------------------------
### Ring geometry
###------------------------------------------------------------

build_three_rings <- function(df_subset) {

  # ---- Ring 1: inner = Total_Points_REVEL bins ----
  inner_ring <- df_subset %>%
    distinct(Total_Points_REVEL, total_variants) %>%
    arrange(Total_Points_REVEL) %>%
    mutate(
      cum = cumsum(total_variants),
      start = lag(cum, default = 0),
      end = cum
    )

  total_sum <- sum(inner_ring$total_variants)

  inner_ring <- inner_ring %>%
    mutate(
      start = start / total_sum * 2*pi,
      end   = end   / total_sum * 2*pi,
      r0 = 0.3,
      r1 = 0.55
    )

  # ---- Ring 2: middle = point_source ----
  middle_ring <- df_subset %>%
    group_by(Total_Points_REVEL, point_source) %>%
    summarise(n_middle = sum(n), .groups = "drop") %>%
    left_join(inner_ring, by = "Total_Points_REVEL") %>%
    group_by(Total_Points_REVEL) %>%
    mutate(
      total_bin = sum(n_middle),
      local_start = lag(cumsum(n_middle), default = 0),
      local_end   = cumsum(n_middle),
      # Capture the bin's own start/end before overwriting them below -- the
      # end formula must divide the *original* bin span, not the already-
      # reassigned start (which it would otherwise pick up in the same
      # mutate(), inflating every category after the first and giving a
      # true-zero-count category a spurious nonzero width).
      bin_start = start,
      bin_end = end,
      start = bin_start + (local_start / total_bin) * (bin_end - bin_start),
      end   = bin_start + (local_end   / total_bin) * (bin_end - bin_start),
      r0 = 0.55,
      r1 = 0.85
    ) %>%
    ungroup() %>%
    select(-bin_start, -bin_end)

  # ---- Ring 3: outer = conflict_status ----
  outer_ring <- df_subset %>%
    group_by(Total_Points_REVEL, conflict_status) %>%
    summarise(n_outer = sum(n), .groups = "drop") %>%
    left_join(inner_ring, by = "Total_Points_REVEL") %>%
    group_by(Total_Points_REVEL) %>%
    mutate(
      total_bin = sum(n_outer),
      local_start = lag(cumsum(n_outer), default = 0),
      local_end   = cumsum(n_outer),
      bin_start = start,
      bin_end = end,
      start = bin_start + (local_start / total_bin) * (bin_end - bin_start),
      end   = bin_start + (local_end   / total_bin) * (bin_end - bin_start),
      r0 = 0.85,
      r1 = 1.15
    ) %>%
    ungroup() %>%
    select(-bin_start, -bin_end)

  list(inner = inner_ring, middle = middle_ring, outer = outer_ring)
}

###------------------------------------------------------------
### Original (arbitrary-size, hand-finished-in-Illustrator) plots
###------------------------------------------------------------

make_three_ring_plot <- function(df_subset, title = NULL) {

  rings <- build_three_rings(df_subset)
  inner_ring  <- rings$inner
  middle_ring <- rings$middle
  outer_ring  <- rings$outer

  inner_label <- inner_ring %>%
  mutate(
    mid = (start + end) / 2,
    x = (r0 + r1)/2 * cos(mid - pi/2),
    y = -(r0 + r1)/2 * sin(mid - pi/2),
    label_text = as.character(Total_Points_REVEL)
  )

  middle_label <- middle_ring %>%
  mutate(
    mid = (start + end)/2,
    x = (r0 + r1)/2 * cos(mid - pi/2),
    y = -(r0 + r1)/2 * sin(mid - pi/2),
    label_text = format_count(n_middle)
  )

  outer_label <- outer_ring %>%
  mutate(
    mid = (start + end)/2,
    x = (r0 + r1)/2 * cos(mid - pi/2),
    y = -(r0 + r1)/2 * sin(mid - pi/2),
    label_text = format_count(n_outer)
  )

  list(
  inner = inner_ring,
  middle = middle_ring,
  outer = outer_ring,
  inner_label = inner_label,
  middle_label = middle_label,
  outer_label = outer_label
)


  ggplot() +
    geom_arc_bar(
      data = inner_ring,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1,
          start = start, end = end, fill = Total_Points_REVEL)
    ) +
    geom_arc_bar(
      data = middle_ring,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1,
          start = start, end = end, fill = point_source)
    ) +
    geom_arc_bar(
      data = outer_ring,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1,
          start = start, end = end, fill = conflict_status)
    ) +
    geom_text(
    data = inner_label,
    aes(x = x, y = y, label = label_text),
    size = 3.8,
    fontface = "bold",
    color = "black"
  ) +
  geom_text(
    data = middle_label,
    aes(x = x, y = y, label = label_text),
    size = 3.2,
    fontface = "bold",
    color = "black"
  ) +
  geom_text(
    data = outer_label,
    aes(x = x, y = y, label = label_text),
    size = 2.7,
    fontface = "bold",
    color = "black"
  ) +
    coord_fixed() +
    theme_void() +
    scale_fill_manual(values = colors_custom) +
    labs(title = title) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

make_three_ring_plot_no_labels <- function(df_subset, title = NULL) {

  rings <- build_three_rings(df_subset)
  inner_ring  <- rings$inner
  middle_ring <- rings$middle
  outer_ring  <- rings$outer

  ggplot() +
    geom_arc_bar(
      data = inner_ring,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1,
          start = start, end = end, fill = Total_Points_REVEL)
    ) +
    geom_arc_bar(
      data = middle_ring,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1,
          start = start, end = end, fill = point_source)
    ) +
    geom_arc_bar(
      data = outer_ring,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1,
          start = start, end = end, fill = conflict_status)
    ) +
    coord_fixed() +
    theme_void() +
    scale_fill_manual(values = colors_custom) +
    labs(title = title) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

###------------------------------------------------------------
### Calibrated donut: real print size, baked-in labels
###------------------------------------------------------------
# Renders at the *exact* final placed size (cairo_pdf, mm units, Arial),
# matching Figure_6b.R's convention -- so at 100% placement scale in the
# assembled figure, text renders at true point size, instead of the
# make_three_ring_plot()/_no_labels() PNG workflow above (arbitrary size,
# default font, hand-finished label placement/sizing in Illustrator). This
# is an additional variant, not a replacement -- keep generating the
# lab/nolab versions above too.

# Device-independent-ish text width in mm, via grid (resolves Arial through
# fontconfig on the r-figures container's default cairo_pdf device).
text_width_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY) {
  gp <- gpar(fontfamily = family, fontsize = pt)
  vapply(
    label,
    function(l) convertWidth(grobWidth(textGrob(l, gp = gp)), "mm", valueOnly = TRUE),
    numeric(1)
  )
}

# angle=0 is 12 o'clock, increasing clockwise -- same convention as the
# inner_label/middle_label/outer_label x/y formulas above.
polar_xy <- function(r, angle) {
  tibble(x = r * cos(angle - pi/2), y = -r * sin(angle - pi/2))
}

arc_length_mm <- function(r, theta, mm_per_unit) r * theta * mm_per_unit
chord_mm <- function(r, theta, mm_per_unit) 2 * r * sin(theta / 2) * mm_per_unit

# For a straight (non-curved) label, the horizontal room actually available
# at its own bisector height is bounded by the wedge's two straight radial
# boundaries -- chord_mm (the straight-line distance between the arc's two
# endpoints) systematically overestimates that room whenever the bisector
# sits near the 90/270-degree sides (same effect best_horizontal_angle
# corrects for by repositioning): there, the boundaries are nearly vertical,
# so the true horizontal gap at the bisector's own height can be much
# smaller than the chord. A boundary line at angle phi crosses height y0 at
# x = y0*tan(phi), so the true gap is |y0|*|tan(end)-tan(start)|. Falls back
# to chord_mm if a boundary sits close enough to 90/270 to blow up tan().
local_room_mm <- function(r, mid, start, end, theta, mm_per_unit) {
  y0 <- r * cos(mid)
  room <- abs(y0) * abs(tan(end) - tan(start)) * mm_per_unit
  chord <- chord_mm(r, theta, mm_per_unit)
  ifelse(is.finite(room), pmin(room, chord), chord)
}

# For a horizontal label at radius r, angle theta (0/2pi = top, pi = bottom,
# pi/2 & 3pi/2 = sides), the room before hitting either straight radial
# boundary scales with |r * cos(theta)| -- maximal toward the top/bottom of
# the wedge, minimal toward the sides. Reposition within [start, end] toward
# whichever of {0, pi, 2pi} is reachable, keeping a margin off both edges so
# the label doesn't crowd the wedge boundary lines.
best_horizontal_angle <- function(start, end, margin_frac = 0.15) {
  theta <- end - start
  margin <- max(theta * margin_frac, 0.02)
  lo <- start + margin
  hi <- end - margin
  if (lo >= hi) return((start + end) / 2)
  candidates <- c(lo, hi)
  ks <- seq(floor(lo / pi), ceiling(hi / pi))
  interior <- ks * pi
  interior <- interior[interior >= lo & interior <= hi]
  candidates <- c(candidates, interior)
  candidates[which.max(abs(cos(candidates)))]
}

WIDE_WEDGE_THRESHOLD_RAD <- 80 * pi / 180

# Only reposition wide wedges (>80 degrees) whose bisector sits close to the
# 90/270-degree sides, where horizontal room is genuinely short -- not near
# the 0/180-degree top/bottom, where it's already generous (chord_mm can't
# tell these apart: it's the same number for a wedge centered at either
# position, even though the true room at the bisector differs a lot -- see
# best_horizontal_angle's |r*cos(theta)| reasoning). Narrow wedges, and wide
# ones already centered near the top/bottom, are left at their bisector.
needs_reposition <- function(start, end, cos_threshold = 0.5) {
  theta <- end - start
  if (theta <= WIDE_WEDGE_THRESHOLD_RAD) return(FALSE)
  mid <- (start + end) / 2
  abs(cos(mid)) < cos_threshold
}

# White text on dark/saturated fills (e.g. the "-12 to -6" navy, "Functional
# only" teal), black text otherwise -- generalizes across panels instead of
# hardcoding which categories are "dark".
text_color_for_fill <- function(fill_value) {
  hex <- colors_custom[fill_value]
  vapply(hex, function(h) {
    if (is.na(h)) return("black")
    rgb <- grDevices::col2rgb(h) / 255
    luminance <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
    if (luminance < 0.5) "white" else "black"
  }, character(1), USE.NAMES = FALSE)
}

# Balanced two-line word-wrap (not greedy) -- for "10 to 12" gives
# c("10", "to 12") rather than piling words onto the first line. Returns
# c(line1, NA) if the label has no space to break on.
wrap_two_lines <- function(label) {
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

make_three_ring_plot_calibrated <- function(df_subset, ring_diameter_mm = RING_DIAMETER_MM) {

  outer_r_max <- 1.15
  mm_per_unit <- ring_diameter_mm / (outer_r_max * 2)

  rings <- build_three_rings(df_subset)

  # ---- classify each label as "inside" (fits) or "outside" (leader line) ----
  # Curved (inner-ring) text follows the arc exactly, so arc length is a
  # precise measure of the space available; straight text at r_mid in a
  # tapering wedge only approximately fits within the chord width, so it
  # keeps a bigger safety margin.
  fit_fudge_straight <- 0.85
  fit_fudge_curved <- 0.95
  EXT_R <- outer_r_max + 0.15   # leader-line/label radius, same for every ring

  classify_ring <- function(ring_df, count_col, fill_col, curved = FALSE, relabel = FALSE) {
    fit_fudge <- if (curved) fit_fudge_curved else fit_fudge_straight
    df <- ring_df %>%
      mutate(theta = end - start) %>%
      filter(theta > 1e-9) %>%   # skip zero-size regions entirely -- nothing to label
      mutate(
        mid = (start + end) / 2,
        r_mid = (r0 + r1) / 2,
        label_text = format_count(.data[[count_col]]),
        label_text = if (relabel) unname(bin_display_labels[label_text]) else label_text,
        text_color = text_color_for_fill(as.character(.data[[fill_col]]))
      )
    # Straight labels reposition toward the top/bottom of their own wedge,
    # where horizontal room is greatest, for wide/crowded wedges (see
    # best_horizontal_angle/needs_reposition) -- computed here, before the
    # fit test, so the fit test itself measures room at the position the
    # label will actually render at, not just its raw (possibly much more
    # crowded) bisector. Curved inner-ring text follows the whole arc, so
    # neither repositioning nor this per-point room measure applies to it.
    if (!curved) {
      df <- df %>%
        rowwise() %>%
        mutate(
          label_angle = if (needs_reposition(start, end)) best_horizontal_angle(start, end) else mid
        ) %>%
        ungroup()
    } else {
      df$label_angle <- df$mid
    }
    df %>%
      mutate(
        avail_mm = if (curved) arc_length_mm(r_mid, theta, mm_per_unit)
                   else local_room_mm(r_mid, label_angle, start, end, theta, mm_per_unit),
        text_mm = text_width_mm(label_text),
        fits = text_mm < avail_mm * fit_fudge
      )
  }

  inner_ring  <- classify_ring(rings$inner,  "Total_Points_REVEL", "Total_Points_REVEL", curved = TRUE, relabel = TRUE)
  middle_ring <- classify_ring(rings$middle, "n_middle", "point_source")
  outer_ring  <- classify_ring(rings$outer,  "n_outer", "conflict_status")

  # For a curved label too wide for one line, try splitting it into two
  # shorter lines stacked within the ring's own radial band, set a size
  # smaller (WRAP_LABEL_PT) -- matches how the original hand-placed figure
  # handled "10 to 12" -- before falling back to an external leader line.
  fit_fudge_wrap <- 0.95
  inner_ring <- inner_ring %>%
    rowwise() %>%
    mutate(
      wrap_line1 = if (!fits) wrap_two_lines(label_text)[1] else NA_character_,
      wrap_line2 = if (!fits) wrap_two_lines(label_text)[2] else NA_character_,
      wrap_fits = !fits && !is.na(wrap_line2) &&
        max(text_width_mm(wrap_line1, pt = WRAP_LABEL_PT), text_width_mm(wrap_line2, pt = WRAP_LABEL_PT)) <
          avail_mm * fit_fudge_wrap
    ) %>%
    ungroup()

  # ---- inside placements ----
  # Straight (non-curved) labels reposition toward the top/bottom of their
  # own wedge, where horizontal room is greatest (see best_horizontal_angle) --
  # but only for wide, genuinely-crowded wedges (see needs_reposition);
  # curved inner-ring text follows the whole arc so this doesn't apply to it.
  inner_inside <- inner_ring %>% filter(fits)
  inner_inside_wrapped <- inner_ring %>% filter(!fits & wrap_fits)
  # label_angle (repositioned toward top/bottom for wide/crowded wedges, or
  # the raw bisector otherwise) was already computed in classify_ring, so the
  # fit test above measures room at the same position the label renders at.
  middle_inside_xy <- middle_ring %>% filter(fits) %>%
    bind_cols(polar_xy(.$r_mid, .$label_angle))
  outer_inside_xy <- outer_ring %>% filter(fits) %>%
    bind_cols(polar_xy(.$r_mid, .$label_angle))

  # ---- outside placements (any ring, fallback for inner too) ----
  # Non-outermost rings anchor their leader line at their own mid-radius,
  # not r1 -- starting exactly on a ring boundary reads as a continuation of
  # that boundary's wedge divisions in the ring(s) outside it. The outermost
  # ring has nothing beyond its own r1 to be confused with, so it keeps it.
  make_outside <- function(ring_df, anchor_at_mid = TRUE) {
    if (nrow(ring_df) == 0) return(ring_df)
    ring_df %>%
      mutate(r_anchor = if (anchor_at_mid) r_mid else r1) %>%
      bind_cols(polar_xy(.$r_anchor, .$mid) %>% rename(x_anchor = x, y_anchor = y)) %>%
      bind_cols(polar_xy(EXT_R, .$mid)) %>%
      mutate(
        hjust = case_when(
          x > 0.02 ~ 0,
          x < -0.02 ~ 1,
          TRUE ~ 0.5
        ),
        # Side labels (hjust 0/1) already clear the leader line -- the line
        # ends at the text's near edge, not its middle. Near-vertical labels
        # (hjust 0.5) need the same treatment on the vertical axis: anchor at
        # the near edge (bottom of text if above the ring, top if below) so
        # the line ends at the label instead of running into its middle.
        vjust = case_when(
          hjust != 0.5 ~ 0.5,
          y > 0 ~ 0,
          TRUE ~ 1
        )
      )
  }
  inner_outside  <- make_outside(inner_ring  %>% filter(!fits & !wrap_fits), anchor_at_mid = TRUE)
  middle_outside <- make_outside(middle_ring %>% filter(!fits), anchor_at_mid = TRUE)
  outer_outside  <- make_outside(outer_ring  %>% filter(!fits), anchor_at_mid = FALSE)
  all_outside <- bind_rows(inner_outside, middle_outside, outer_outside)

  # ---- canvas: keep the ring itself pixel-identical across panels, and make
  # sure it's big enough for every rendered label's actual footprint -- not
  # just outside/leader labels, but also inside labels that happen to sit
  # near the horizontal/vertical extreme of the ring (e.g. a wide count right
  # at the 9 o'clock position), which can otherwise get clipped by the
  # device's own canvas edge even though the ggplot panel itself has
  # clip = "off". hjust/vjust determine how a label's width/height is
  # distributed around its anchor (0 = anchor at that edge, so it all reads
  # out from there; 0.5 = centered; 1 = anchor at the other edge).
  half_h_r <- (LABEL_MM / 2) / mm_per_unit
  label_bboxes <- bind_rows(
    middle_inside_xy %>% transmute(x, y, w_r = text_mm / mm_per_unit, hjust = 0.5, vjust = 0.5),
    outer_inside_xy  %>% transmute(x, y, w_r = text_mm / mm_per_unit, hjust = 0.5, vjust = 0.5),
    all_outside      %>% transmute(x, y, w_r = text_mm / mm_per_unit, hjust, vjust)
  )
  if (nrow(label_bboxes) > 0) {
    label_bboxes <- label_bboxes %>%
      mutate(
        left   = x - w_r * hjust,
        right  = x + w_r * (1 - hjust),
        bottom = y - half_h_r * 2 * vjust,
        top    = y + half_h_r * 2 * (1 - vjust)
      )
    x_range <- c(min(-EXT_R, label_bboxes$left), max(EXT_R, label_bboxes$right))
    y_range <- c(min(-EXT_R, label_bboxes$bottom), max(EXT_R, label_bboxes$top))
  } else {
    x_range <- c(-EXT_R, EXT_R)
    y_range <- c(-EXT_R, EXT_R)
  }

  p <- ggplot() +
    geom_arc_bar(
      data = rings$inner,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end,
          fill = Total_Points_REVEL),
      linewidth = 0.5 * point_in_mm, color = "black"
    ) +
    geom_arc_bar(
      data = rings$middle,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end,
          fill = point_source),
      linewidth = 0.5 * point_in_mm, color = "black"
    ) +
    geom_arc_bar(
      data = rings$outer,
      aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end,
          fill = conflict_status),
      linewidth = 0.5 * point_in_mm, color = "black"
    )

  if (nrow(inner_inside) > 0) {
    # geomtextpath's vjust=0.5 centers the font's full ascent+descent box on
    # the path; our labels (digits/"to"/minus sign) never have descenders, so
    # the visible glyphs sit high (toward larger radius) within that box.
    # Nudge the path radius inward by a fixed physical amount to compensate --
    # proportional to font size, not ring size, since it's a font-metrics
    # effect, so it holds across panels/ring diameters.
    inward_nudge_mm <- 0.3
    arc_pts <- inner_inside %>%
      rowwise() %>%
      reframe(
        Total_Points_REVEL = Total_Points_REVEL,
        label_text = label_text,
        angle = seq(start, end, length.out = 40),
        r_mid = r_mid - inward_nudge_mm / mm_per_unit,
        text_color = text_color
      ) %>%
      bind_cols(polar_xy(.$r_mid, .$angle))
    p <- p + geom_textpath(
      data = arc_pts,
      aes(x = x, y = y, label = label_text, group = Total_Points_REVEL,
          color = text_color),
      size = LABEL_MM, family = FONT_FAMILY,
      upright = TRUE, hjust = 0.5, text_only = TRUE
    )
  }
  if (nrow(inner_inside_wrapped) > 0) {
    # Font-metrics nudge scales with font size (see the single-line case
    # above); line_gap is snug (1x font size, center-to-center) since these
    # are short fragments that read fine close together, and it directly
    # trades off against how much radial room is left within the ring band.
    inward_nudge_mm <- 0.3 * (WRAP_LABEL_PT / LABEL_PT)
    line_gap_mm <- WRAP_LABEL_MM * 1.0
    half_gap_r <- (line_gap_mm / 2) / mm_per_unit
    base_r <- inner_inside_wrapped$r_mid - inward_nudge_mm / mm_per_unit
    # Which physical radius (page-space, not path-space) reads as "line 1" (on
    # top) depends on which half of the circle the wedge sits in: in the
    # upper half, larger radius = higher on the page; in the lower half, it's
    # reversed. polar_xy's y = r*cos(angle) at the bisector gives that sign.
    page_y_sign <- sign(cos((inner_inside_wrapped$start + inner_inside_wrapped$end) / 2))
    build_wrapped_line <- function(line_col, r_sign) {
      inner_inside_wrapped %>%
        mutate(r_mid_line = base_r + r_sign * page_y_sign * half_gap_r) %>%
        rowwise() %>%
        reframe(
          Total_Points_REVEL = Total_Points_REVEL,
          label_text = .data[[line_col]],
          angle = seq(start, end, length.out = 40),
          r_mid = r_mid_line,
          text_color = text_color
        ) %>%
        bind_cols(polar_xy(.$r_mid, .$angle))
    }
    arc_pts_l1 <- build_wrapped_line("wrap_line1", +1)
    arc_pts_l2 <- build_wrapped_line("wrap_line2", -1)
    p <- p +
      geom_textpath(
        data = arc_pts_l1,
        aes(x = x, y = y, label = label_text, group = Total_Points_REVEL,
            color = text_color),
        size = WRAP_LABEL_MM, family = FONT_FAMILY,
        upright = TRUE, hjust = 0.5, text_only = TRUE
      ) +
      geom_textpath(
        data = arc_pts_l2,
        aes(x = x, y = y, label = label_text, group = Total_Points_REVEL,
            color = text_color),
        size = WRAP_LABEL_MM, family = FONT_FAMILY,
        upright = TRUE, hjust = 0.5, text_only = TRUE
      )
  }
  if (nrow(middle_inside_xy) > 0) {
    p <- p + geom_text(
      data = middle_inside_xy,
      aes(x = x, y = y, label = label_text, color = text_color),
      size = LABEL_MM, family = FONT_FAMILY
    )
  }
  if (nrow(outer_inside_xy) > 0) {
    p <- p + geom_text(
      data = outer_inside_xy,
      aes(x = x, y = y, label = label_text, color = text_color),
      size = LABEL_MM, family = FONT_FAMILY
    )
  }
  if (nrow(all_outside) > 0) {
    p <- p +
      geom_segment(
        data = all_outside,
        aes(x = x_anchor, y = y_anchor, xend = x, yend = y),
        color = "grey40", linewidth = 0.5 * point_in_mm
      ) +
      geom_text(
        data = all_outside,
        aes(x = x, y = y, label = label_text, hjust = hjust, vjust = vjust),
        size = LABEL_MM, family = FONT_FAMILY, color = "black"
      )
  }

  p <- p +
    coord_fixed(xlim = x_range, ylim = y_range, expand = FALSE, clip = "off") +
    theme_void() +
    scale_fill_manual(values = colors_custom) +
    scale_color_identity() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

  list(
    plot = p,
    width_mm = diff(x_range) * mm_per_unit,
    height_mm = diff(y_range) * mm_per_unit,
    ring_diameter_mm = 2 * outer_r_max * mm_per_unit,
    x_range = x_range, y_range = y_range, label_bboxes = label_bboxes,
    middle_inside_xy = middle_inside_xy
  )
}
