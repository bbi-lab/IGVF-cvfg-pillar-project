# Shared calibrated (real-print-size) Sankey chart logic -- Fig 5 (controls/
# ClinGen, name-labels-only) and Fig 6 (VUS/gnomAD/Unobserved, name+count),
# and any future supplementary sankeys built the same way.
#
# Source this from a sankey-chart chunk before calling make_sankey_calibrated().
# Unlike three_ring_donut.R, colors_custom is an explicit function argument
# here rather than looked up from the caller's environment -- both work, this
# is just the more explicit of the two conventions already used in this repo.
#
# Expects the caller to have already loaded: dplyr, tidyr, ggplot2, ggsankey,
# extrafont (+ loadfonts(device = 'all')), and grid.
#
# ggsankey's coordinate system (established empirically):
# - Node columns sit at integer x positions (1, 2, ...), each spanning
#   [x - NODE_HALF_WIDTH, x + NODE_HALF_WIDTH] (default half-width 0.05).
# - A node's height in y is proportional to its own count; ggsankey inserts
#   gaps between stacked nodes in a multi-node column so the column's total
#   span exceeds the raw sum of its nodes' heights.
# - Flows are drawn as smoothed ribbons between column edges, not needed for
#   our own label placement.

point_in_mm <- 0.3527778
LABEL_PT <- 7
LABEL_MM <- LABEL_PT * point_in_mm
FONT_FAMILY <- "Arial"
NODE_HALF_WIDTH <- 0.05
# Empty vertical space between the bottom of the name box and the top of
# the count box, matching the original hand-built chart's own snug (but
# still clearly visible) stacking. Since both boxes are now precisely-sized
# grobs (see make_box_grob()), this gap renders identically everywhere --
# no more renderer-dependent auto-sizing to leave extra margin for.
BOX_GAP_MM <- 0.7
# Empty space around the text within a label box, on all 4 sides -- used
# only to compute the box's own target rectangle (see make_box_grob()),
# never passed to any geom's own auto-sizing.
PAD_MM <- 0.25
# Tightened below ggplot2's own geom_label default (1.2) to keep a two-line
# wrapped label (e.g. "Likely\npathogenic") compact -- matching the original
# hand-built chart's snug line spacing. Must be passed to both geom_label()/
# geom_text() (the actual rendering) and text_height_mm() (so the box-height
# math matches what's actually drawn).
LINEHEIGHT <- 0.85

# A node/box border is stroked centered on its own nominal edge, so roughly
# half its width extends beyond that edge -- for whichever node/label/count
# box ends up flush against the canvas boundary (the common case, since the
# canvas is sized to the tightest fit around every box), that outer half
# would otherwise get clipped by the PDF page edge itself. This margin is
# added on all 4 sides of the canvas to leave room for it; comfortably more
# than the actual half-stroke-width (~0.09mm, from the 0.5pt hairline used
# for both node rects and box borders) rather than cutting it exactly.
CANVAS_MARGIN_MM <- 0.2

# Label/count box backgrounds are precisely-sized round-rect grobs (see
# GeomRoundBox/geom_round_box() below), not geom_label()'s own auto-sizing
# from text + label.padding: a PDF viewer can interpret that padding more
# generously than our own R/grid font metrics do, even when text POSITIONS
# match exactly between the two (confirmed: reported directly by a viewer
# screenshot showing a label box's rounded rectangle physically overlapping
# the count box below it, despite our computed coordinates showing a clean
# gap). A grob built at an exact physical size has no auto-sizing step left
# for any renderer to interpret differently.
BOX_RADIUS_MM <- 0.3
# grid's "lwd" is device-dependent but conventionally 1/96 inch per unit;
# convert our usual hairline (0.5pt, matching the donut charts' and node
# rects' own border) into that unit.
BOX_BORDER_LWD <- (0.5 * point_in_mm) / (25.4 / 96)

# A round-rect box at a data-space (x, y) CENTER but a true physical
# (width_mm, height_mm, radius) size. Deliberately not built via
# annotation_custom() (map a grob's own 0-1 npc square onto a target data
# rectangle): with this chart's wildly anisotropic data units (x spans ~1-2,
# y spans raw variant counts in the tens of thousands), that remapping
# distorts the round-rect's corner geometry into visibly square corners.
# coord$transform() converts only the box's CENTER point to npc (a plain
# translation, not a shape-distorting rescale); width/height/r stay in "mm"
# throughout, so the box's physical shape can't be affected by the panel's
# own data-to-npc scale factors.
GeomRoundBox <- ggplot2::ggproto("GeomRoundBox", ggplot2::Geom,
  required_aes = c("x", "y", "width_mm", "height_mm", "fill"),
  default_aes = ggplot2::aes(colour = "black"),
  draw_key = ggplot2::draw_key_rect,
  draw_panel = function(data, panel_params, coord) {
    coords <- coord$transform(data, panel_params)
    grobs <- lapply(seq_len(nrow(coords)), function(i) {
      grid::roundrectGrob(
        x = coords$x[i], y = coords$y[i],
        width = unit(coords$width_mm[i], "mm"), height = unit(coords$height_mm[i], "mm"),
        r = unit(BOX_RADIUS_MM, "mm"),
        gp = gpar(fill = coords$fill[i], col = coords$colour[i], lwd = BOX_BORDER_LWD)
      )
    })
    do.call(grid::grobTree, grobs)
  }
)
geom_round_box <- function(mapping = NULL, data = NULL, ...) {
  ggplot2::layer(geom = GeomRoundBox, mapping = mapping, data = data,
                 stat = "identity", position = "identity", inherit.aes = FALSE, ...)
}

format_count <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

text_width_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY) {
  gp <- gpar(fontfamily = family, fontsize = pt)
  vapply(label, function(l) convertWidth(grobWidth(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

text_height_mm <- function(label, pt = LABEL_PT, family = FONT_FAMILY, lineheight = LINEHEIGHT) {
  gp <- gpar(fontfamily = family, fontsize = pt, lineheight = lineheight)
  vapply(label, function(l) convertHeight(grobHeight(textGrob(l, gp = gp)), "mm", valueOnly = TRUE), numeric(1))
}

# Sentence case (first character capitalized, rest lowercase) -- for display
# only, never applied to the underlying category value used for color
# lookup/grouping. Node names outside gnomAD/preclass's hand-authored source
# labels (see make_sankey_calibrated's source_label) never contain "gnomAD",
# so no special-casing is needed here.
to_sentence_case <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))

# Balanced two-line word-wrap (not greedy) -- for "10 to 12" gives
# c("10", "to 12") rather than piling words onto the first line. Returns
# c(label, NA) if the label has no space to break on. Two-word labels (all
# of this project's source-node overrides) always split after the first
# word.
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

# Vectorized wrap_two_lines() for a whole label column -- single-word labels
# (no space to break on) pass through unchanged rather than gaining a
# trailing "\nNA".
wrap_labels_if_multiword <- function(labels) {
  vapply(labels, function(l) {
    w <- wrap_two_lines(l)
    if (is.na(w[2])) w[1] else paste(w, collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

# White text on dark/saturated fills (e.g. a navy "Benign"), black text
# otherwise (including mid-gray fallback fills) -- takes the already-resolved
# hex color directly, not a name lookup, since every box (label and count
# alike) now shares that color. Threshold sits below grey50's own luminance
# (~0.498) so a mid-gray fallback fill stays black, and above the point-value
# palette's own +6 ("#B85C6B", luminance ~0.4754) and -11 ("#1D7AAB",
# luminance ~0.3913) -- reproduced: at the old 0.45 threshold, +6 through
# +12 and -11/-12 rendered correctly, but +6 itself and -9/-10 (whose
# luminance, ~0.574/~0.496, both sit above +6's) still got black text
# despite their dark, saturated fills. Raising the threshold to 0.48 fixes
# +6 without disturbing any categorical chart's own text color (checked
# against every Pathogenic/Benign/etc. color already in use). Note this
# means white text now starts at -11, not -9 -- -9/-10's colors are
# measurably lighter than +6's by this same luminance formula, so no single
# threshold can treat +6 as "dark" while treating -9/-10 as "light".
text_color_for_fill <- function(hex) {
  vapply(hex, function(h) {
    rgb <- grDevices::col2rgb(h) / 255
    luminance <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
    if (luminance < 0.48) "white" else "black"
  }, character(1), USE.NAMES = FALSE)
}

# Flows can be genuinely semi-transparent (TRUE) or opaque using pre-blended
# colors (FALSE).
FLOW_TRANSPARENT <- TRUE
FLOW_ALPHA <- 0.7
blend_over_white <- function(hex, alpha) {
  rgb_vals <- grDevices::col2rgb(hex)
  blended <- alpha * rgb_vals + (1 - alpha) * 255
  grDevices::rgb(blended[1, ], blended[2, ], blended[3, ], maxColorValue = 255)
}

# source_col/points_col: unquoted column names (tidy-eval, like
# ggsankey::make_long()'s own arguments) -- source_col is the left-hand
# (stage 1) category for this chart (e.g. clinvar_sig_2025, or a literal
# constant column the caller adds beforehand, like gnomad_df$gnomad <-
# 'Gnomad'). levels_order must list every category appearing in EITHER
# column, in the chart's intended top-to-bottom node order. source_label, if
# given, overrides the *display* text of the source (stage 1) node only --
# word-wrapped to two lines and used verbatim (not sentence-cased, so it can
# carry its own styling, e.g. "gnomAD variant") -- without changing the
# underlying category value used for color lookup/grouping (Fig 6's charts,
# whose source is a single synthetic category, not a real one). Every other
# label is shown in sentence case. wrap_destination_labels wraps every
# multi-word label to two lines too (e.g. "Likely benign" -> "Likely
# \nbenign") -- for narrower charts; when source_label is also given, only
# non-source labels are wrapped this way (the source already got its own
# wrapped override). label_pt sets the font size (Fig 6 uses 7pt, Fig 5 uses
# 6pt); show_count draws a separate count box below each name box (Fig 6)
# or, when FALSE, just the name (Fig 5, which also has multiple distinct
# source categories rather than one dominant one -- see the flow-color
# comment below for why that's handled generically either way).
# center_label_nodes names (post sentence-case, e.g. "Likely pathogenic")
# any node whose name label should center on the node itself rather than
# block-center with its count, even though the node's own height would
# otherwise qualify it for the normal block-centered layout -- for a node
# tightly stacked next to another whose own label/count reaches into its
# space (see the is_small comment below for the reproduced case this was
# added for). label_overrides is a named vector keyed by a node's raw
# underlying value (e.g. c("0" = "No evidence")) whose display text should
# be replaced outright, verbatim -- for a numeric points chart, where "0"
# is really a stand-in for "no functional or predictive evidence at all"
# rather than a literal zero score. white_text_nodes names (matched the
# same way as center_label_nodes, post sentence-case/wrap) any node that
# should get white label text regardless of what text_color_for_fill()'s
# shared luminance threshold would otherwise choose for its fill color --
# for specific point-value colors a caller has judged too dark for black
# text but not dark enough to trip that global threshold, without raising
# the threshold itself (which would also affect every other node sharing
# a similar-luminance fill in every other chart, point-based or not,
# including this chart's own gray-fallback source node). label_box_height_scale
# multiplies the label (and count, if shown) box's own computed height --
# text height plus PAD_MM on top and bottom -- by this factor, leaving the
# box's width untouched; 1 (the default) reproduces the original unscaled
# height exactly, since (text_height_mm + 2*PAD_MM)/2*1 ==
# text_height_mm/2 + PAD_MM. clamp_label_right_edge caps every label's own
# right edge (text and background box together, as a unit -- shifting
# both left of the node's own center_x by the same amount, never shrinking
# either one) at its node's own xmax; every node in a ggsankey column
# shares the same xmax, so once every label in that column is clamped, a
# single shared x position just past that column edge is guaranteed clear
# of all of them regardless of any individual label's own width -- when a
# label already fits within its node's own width, this is a no-op and its
# default center_x positioning is unchanged. show_target_percent draws
# each non-source node's own share of
# the total source-side count (e.g. "12%", black text, same size as the
# node labels) just to the right of its (clamped) label box -- meant to
# be used together with clamp_label_right_edge, which is what makes a
# single shared percent-column position possible in the first place.
make_sankey_calibrated <- function(df, source_col, points_col, colors_custom, levels_order,
                                    width_mm, height_mm, source_label = NULL,
                                    wrap_destination_labels = FALSE,
                                    label_pt = LABEL_PT, show_count = TRUE,
                                    center_label_nodes = character(0),
                                    label_overrides = character(0),
                                    white_text_nodes = character(0),
                                    label_box_height_scale = 1,
                                    clamp_label_right_edge = FALSE,
                                    show_target_percent = FALSE) {
  # text_width_mm()/text_height_mm() (via grid::convertWidth/Height on a
  # textGrob) measure against whatever graphics device is currently active --
  # font metrics are device-dependent, and a knitr chunk device active at
  # call time isn't guaranteed to agree with cairo_pdf (what ggsave() below
  # actually renders with), which can silently under-measure text width and
  # produce a too-narrow canvas (reproduced: running this same function from
  # a plain Rscript vs. from inside an .Rmd knit gave different measured
  # widths for the same string). Open a dedicated cairo_pdf(NULL) device for
  # every measurement in this function, so results match the real render
  # regardless of whatever device the caller already has open.
  measure_dev_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(measure_dev_file)
  on.exit({
    grDevices::dev.off()
    unlink(measure_dev_file)
  }, add = TRUE)

  label_mm <- label_pt * point_in_mm

  source_col <- rlang::ensym(source_col)
  points_col <- rlang::ensym(points_col)

  summary_data <- df %>%
    group_by(!!source_col, !!points_col) %>%
    summarise(count = n(), .groups = "drop")
  expanded_df <- summary_data[rep(1:nrow(summary_data), summary_data$count), -3]
  rownames(expanded_df) <- NULL
  sankey_data <- expanded_df %>% ggsankey::make_long(!!source_col, !!points_col)

  sankey_data$node <- factor(sankey_data$node, levels = levels_order)
  sankey_data$next_node <- factor(sankey_data$next_node, levels = levels_order)
  # One scale_fill_identity() covers both the flow/node layers and the label
  # boxes below, since the label boxes' fill is already a resolved hex color
  # (can't mix a scale_fill_manual()-mapped layer with a literal-color layer
  # under a single fill scale). scale_fill_identity() has no na.value
  # fallback, so replicate one explicitly for any node not in colors_custom
  # (e.g. a chart's own literal source label, like "Gnomad") -- matching the
  # same #A0A0A0 gray already used for "Uncertain"/"No Classification", not
  # R's builtin "grey50" (#808080, visibly darker).
  resolve_fill <- function(node) {
    hex <- unname(colors_custom[as.character(node)])
    ifelse(is.na(hex), "#A0A0A0", hex)
  }
  # Build once with a temporary geom_sankey_label layer purely to extract
  # ggsankey's own node/flow geometry -- its internal stacking+gap algorithm
  # isn't something worth reimplementing. Fill is mapped to the plain `node`
  # factor here (not our own resolved hex colors): StatSankeyNode/
  # StatSankeyFlow's setup_data() internally does group_by_all() (grouping
  # by every mapped aes, fill included) then relies on the resulting row
  # order/grouping for its cumsum-based stacking -- if two DIFFERENT nodes
  # happen to resolve to the SAME literal fill color (e.g. two different
  # gray fallbacks unified to the same #A0A0A0), group_by_all() merges their
  # rows and corrupts the stacking (reproduced empirically: with a shared
  # gray, one node's flow crossed another's). `node` is guaranteed unique
  # per row identity, so it's the only safe thing to feed the stat; our own
  # colors are applied afterward, downstream of this extraction, never fed
  # back into a Stat computation.
  p_temp <- ggplot(sankey_data, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                                     fill = node, label = node)) +
    geom_sankey(flow.alpha = 0.7, node.color = "black", width = NODE_HALF_WIDTH * 2) +
    geom_sankey_label(size = 3)
  built <- ggplot_build(p_temp)
  node_geom <- built$data[[2]] %>%
    select(x, label, freq, xmin, xmax, ymin, ymax) %>%
    distinct() %>%
    mutate(
      fill_color = resolve_fill(label),
      # label_overrides is keyed by the raw underlying value (e.g. "0"),
      # looked up before sentence-casing, and used verbatim (not
      # sentence-cased) so the caller's own capitalization always wins --
      # for a numeric points chart's "0" node, which is really a stand-in
      # for "no functional or predictive evidence at all" rather than a
      # literal zero score, "No evidence" is a clearer label than a bare
      # digit. Doesn't touch fill_color's own lookup just above, which
      # still keys off the raw value.
      label = ifelse(
        as.character(label) %in% names(label_overrides),
        unname(label_overrides[as.character(label)]),
        to_sentence_case(as.character(label))
      ),
      center_x = (xmin + xmax) / 2,
      center_y = (ymin + ymax) / 2
    )
  is_source <- node_geom$x == min(node_geom$x)
  if (!is.null(source_label)) {
    node_geom$label[is_source] <- paste(wrap_two_lines(source_label), collapse = "\n")
  }
  if (wrap_destination_labels) {
    needs_wrap <- if (is.null(source_label)) rep(TRUE, nrow(node_geom)) else !is_source
    node_geom$label[needs_wrap] <- wrap_labels_if_multiword(node_geom$label[needs_wrap])
  }

  # Flow ribbon paths (shape) come straight from built$data[[1]] (their shape
  # doesn't depend on color), but their FILL is resolved per flow via each
  # flow's own SOURCE node -- not a single shared color -- since a chart's
  # source side can have several distinct categories (e.g. Fig 5's ClinVar/
  # ClinGen columns share the same category vocabulary on both sides, unlike
  # Fig 6's one dominant synthetic source). built$data[[1]]'s "label" column
  # (from p_temp's own label=node mapping) is constant within each flow's
  # `group` and equal to that flow's source category (confirmed empirically:
  # grouping the built flow data by `group` and counting distinct `label`
  # values always gives 1) -- reused here as a reliable group -> source-node
  # lookup, entirely downstream of the geometry extraction above so it can't
  # feed back into ggsankey's own stacking computation.
  flow_source <- built$data[[1]] %>%
    distinct(group, source_node = label) %>%
    mutate(flow_fill_color = if (FLOW_TRANSPARENT) resolve_fill(source_node)
                             else blend_over_white(resolve_fill(source_node), FLOW_ALPHA))
  flow_geom <- built$data[[1]] %>%
    select(x, y, group) %>%
    left_join(select(flow_source, group, flow_fill_color), by = "group")

  x_min <- min(node_geom$xmin); x_max <- max(node_geom$xmax)
  y_min <- min(node_geom$ymin); y_max <- max(node_geom$ymax)
  mm_per_x <- width_mm / (x_max - x_min)
  mm_per_y <- height_mm / (y_max - y_min)

  # Name label box height/width always computed from label_pt. With
  # show_count, stack the name label above a separate count box, both
  # sharing the node's own fill color -- same snug (1x font size) stacking
  # convention used for the donut charts' two-line wrap. For most nodes the
  # label+count pair is centered as a block on the node's own vertical
  # midpoint (label offset up by half_gap_y from center, count BOX_GAP_MM
  # below the label's own bottom edge). A node too short for that offset to
  # keep the label's own center within the node's own span instead centers
  # just the label on the node, with the count hung the normal distance
  # below it -- also forced for any node named in center_label_nodes,
  # regardless of its own height: this block-centering only accounts for a
  # node's own size, not a neighboring node's label/count reaching into its
  # space, so a node whose block-centered label collides with a tightly
  # adjacent node's own label/count (reproduced: Ext. Data Fig 7's AM/MP2
  # charts, "Likely pathogenic"'s label overlapping "Pathogenic"'s count)
  # needs the same override even though it isn't "small" on its own. Box
  # half-heights are computed per row (not a single shared constant) so a
  # two-line-wrapped label -- taller than the usual single-line box -- is
  # handled correctly too. Without show_count, there's no pair to center:
  # the (only) name box just centers on the node.
  if (show_count) {
    node_geom <- node_geom %>%
      mutate(
        count_text = format_count(freq),
        label_box_half_h_mm = (text_height_mm(label, pt = label_pt) + 2 * PAD_MM) / 2 * label_box_height_scale,
        count_box_half_h_mm = (text_height_mm(count_text, pt = label_pt) + 2 * PAD_MM) / 2 * label_box_height_scale,
        half_gap_y = ((label_box_half_h_mm + count_box_half_h_mm) / 2 + BOX_GAP_MM / 2) / mm_per_y,
        # gsub() undoes wrap_destination_labels' embedded newline (e.g.
        # "Likely\npathogenic") before matching against center_label_nodes,
        # which callers specify as the plain unwrapped name -- matching the
        # raw `label` column here would silently never fire for any
        # multi-word name once wrapped (reproduced: "Likely pathogenic" was
        # already wrapped to two lines by this point, so `label %in%
        # center_label_nodes` compared "Likely\npathogenic" against "Likely
        # pathogenic" and never matched).
        is_small = (ymax - ymin) * mm_per_y < (count_box_half_h_mm - label_box_half_h_mm + BOX_GAP_MM) |
          gsub("\n", " ", label, fixed = TRUE) %in% center_label_nodes,
        label_y = ifelse(is_small, center_y, center_y + half_gap_y),
        count_y = label_y - 2 * half_gap_y,
        label_w_mm = text_width_mm(label, pt = label_pt),
        count_w_mm = text_width_mm(count_text, pt = label_pt),
        text_color = ifelse(
          gsub("\n", " ", label, fixed = TRUE) %in% white_text_nodes,
          "white",
          text_color_for_fill(fill_color)
        )
      )
  } else {
    node_geom <- node_geom %>%
      mutate(
        label_box_half_h_mm = (text_height_mm(label, pt = label_pt) + 2 * PAD_MM) / 2 * label_box_height_scale,
        label_y = center_y,
        label_w_mm = text_width_mm(label, pt = label_pt),
        text_color = ifelse(
          gsub("\n", " ", label, fixed = TRUE) %in% white_text_nodes,
          "white",
          text_color_for_fill(fill_color)
        )
      )
  }

  # Canvas: node span is exactly width_mm/height_mm; grow only as needed so
  # every label/count box's *actual* footprint fits -- not just a single
  # global "widest box vs total width" heuristic, which ignores where a box
  # is actually centered (a label near the left edge can overflow even if
  # it's narrower than the full canvas).
  node_geom <- node_geom %>%
    mutate(
      # box_x is the label's own rendering x-center -- both its background
      # box (geom_round_box) and its text (geom_text) render at box_x, not
      # center_x directly, so clamp_label_right_edge (below) can shift the
      # whole label left as a unit without the text overflowing past its
      # own (now narrower-reaching) box. Defaults to center_x unchanged
      # when clamping is off or not needed for a given row.
      box_x = center_x,
      label_left   = box_x - (label_w_mm / 2 + PAD_MM) / mm_per_x,
      label_right  = box_x + (label_w_mm / 2 + PAD_MM) / mm_per_x,
      label_top    = label_y + label_box_half_h_mm / mm_per_y,
      label_bottom = label_y - label_box_half_h_mm / mm_per_y
    )
  if (clamp_label_right_edge) {
    node_geom <- node_geom %>%
      mutate(
        overflow_mm = pmax(0, label_right - xmax),
        box_x = box_x - overflow_mm,
        label_left = label_left - overflow_mm,
        label_right = label_right - overflow_mm
      ) %>%
      select(-overflow_mm)
  }
  if (show_target_percent) {
    # Each non-source node's own share of the total source-side count,
    # right-aligned to a single shared edge rather than left-aligned
    # individually right after each node's own (possibly
    # clamp_label_right_edge-shifted) label -- reproduced: left-aligning
    # each one right after its own label left a numeric node's percentage
    # sitting close enough to its thin node's own small ggsankey-drawn
    # connector-line stub to visually read as overlapping it, even though
    # a wide label like "Conflicting evidence" (which needs the most
    # clamping) had plenty of clearance. The shared edge is the widest of
    # every node's own *natural* left-aligned-with-a-BOX_GAP_MM-gap right
    # edge -- i.e. wherever the widest label's own percentage would have
    # landed anyway -- so every other (narrower) node's percentage gets
    # pushed out to that same generous clearance instead.
    total_source_freq <- sum(node_geom$freq[is_source])
    node_geom$pct_label <- ifelse(
      is_source, NA_character_,
      sprintf("%.0f%%", node_geom$freq / total_source_freq * 100)
    )
    natural_pct_left <- node_geom$label_right + BOX_GAP_MM / mm_per_x
    pct_w_mm <- rep(0, nrow(node_geom))
    pct_w_mm[!is_source] <- text_width_mm(node_geom$pct_label[!is_source], pt = label_pt)
    node_geom$pct_right <- max((natural_pct_left + pct_w_mm / mm_per_x)[!is_source])
  }
  if (show_count) {
    node_geom <- node_geom %>%
      mutate(
        count_left   = center_x - (count_w_mm / 2 + PAD_MM) / mm_per_x,
        count_right  = center_x + (count_w_mm / 2 + PAD_MM) / mm_per_x,
        count_top    = count_y + count_box_half_h_mm / mm_per_y,
        count_bottom = count_y - count_box_half_h_mm / mm_per_y
      )
  }
  x_range <- c(
    min(x_min, node_geom$label_left, if (show_count) node_geom$count_left) - CANVAS_MARGIN_MM / mm_per_x,
    max(x_max, node_geom$label_right, if (show_count) node_geom$count_right,
        if (show_target_percent) node_geom$pct_right) + CANVAS_MARGIN_MM / mm_per_x
  )
  y_range <- c(
    min(y_min, node_geom$label_bottom, if (show_count) node_geom$count_bottom) - CANVAS_MARGIN_MM / mm_per_y,
    max(y_max, node_geom$label_top, if (show_count) node_geom$count_top) + CANVAS_MARGIN_MM / mm_per_y
  )

  p <- ggplot() +
    # Flow ribbons via plain geom_polygon() reusing the bezier-smoothed
    # shapes already extracted above (flow_geom), colored per flow (see
    # flow_source above). alpha is only a live transparency group when
    # FLOW_TRANSPARENT (flow_fill_color is already opaque-equivalent
    # otherwise, so alpha=1 there is a no-op, not a second blend).
    geom_polygon(
      data = flow_geom, inherit.aes = FALSE,
      aes(x = x, y = y, group = group, fill = flow_fill_color),
      alpha = if (FLOW_TRANSPARENT) FLOW_ALPHA else 1
    ) +
    # Node rects via plain geom_rect() reusing node_geom's own xmin/xmax/
    # ymin/ymax (already extracted above for label positioning), rather than
    # a second StatSankeyNode layer -- simpler, and avoids recomputing the
    # same geometry twice. Hairline border matches the donut charts' own
    # geom_arc_bar(linewidth = 0.5 * point_in_mm) convention.
    geom_rect(
      data = node_geom, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color,
          linewidth = I(0.5 * point_in_mm)),
      colour = "black"
    ) +
    scale_fill_identity()

  # Label/count box backgrounds -- see GeomRoundBox/geom_round_box(). One row
  # per box (name, plus count when show_count), each carrying its own center
  # (data space) and true physical width/height (mm).
  boxes <- node_geom %>%
    transmute(x = box_x, y = label_y, width_mm = label_w_mm + 2 * PAD_MM,
              height_mm = 2 * label_box_half_h_mm, fill = fill_color)
  if (show_count) {
    boxes <- bind_rows(
      boxes,
      node_geom %>% transmute(x = center_x, y = count_y, width_mm = count_w_mm + 2 * PAD_MM,
                               height_mm = 2 * count_box_half_h_mm, fill = fill_color)
    )
  }

  p <- p +
    geom_round_box(data = boxes, aes(x = x, y = y, width_mm = width_mm, height_mm = height_mm, fill = fill)) +
    geom_text(
      data = node_geom, inherit.aes = FALSE,
      aes(x = box_x, y = label_y, label = label, color = text_color),
      size = label_mm, family = FONT_FAMILY, lineheight = LINEHEIGHT
    )
  if (show_count) {
    p <- p +
      geom_text(
        data = node_geom, inherit.aes = FALSE,
        aes(x = center_x, y = count_y, label = count_text, color = text_color),
        size = label_mm, family = FONT_FAMILY, lineheight = LINEHEIGHT
      )
  }
  if (show_target_percent) {
    p <- p +
      geom_text(
        data = node_geom[!is_source, ], inherit.aes = FALSE,
        aes(x = pct_right, y = label_y, label = pct_label),
        size = label_mm, family = FONT_FAMILY, colour = "black", hjust = 1
      )
  }

  p <- p +
    scale_color_identity() +
    coord_cartesian(xlim = x_range, ylim = y_range, expand = FALSE, clip = "off") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

  list(
    plot = p,
    width_mm = diff(x_range) * mm_per_x,
    height_mm = diff(y_range) * mm_per_y,
    node_geom = node_geom
  )
}
