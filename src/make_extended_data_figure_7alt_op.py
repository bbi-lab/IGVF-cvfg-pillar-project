#!/usr/bin/env python3
"""Build the OddsPath/Universal-calibration counterpart of Extended Data Figure
7alt: a calibrated (exact-print-size) heatmap of final variant classifications
across VUS, gnomAD, and unobserved variants.

Identical layout to `make_extended_data_figure_7alt.py` (9 rows: three
REVEL/AM/MP2 rows per variant-set group, with a gap between groups, x 5
classification columns P/LP/VUS/LB/B, plus a Total column; each cell shows
its count and row-wise percentage on a shared monochromatic 0-100% scale),
but reading `Supplementary_Data_6.xlsx`'s OddsPath-based
`{VUS,gnomAD,Unobserved}_{REVEL,AM,MP2}_OP` sheets and their
`Class_OP_{REVEL,AM,MP2}` classification column instead -- "functional
evidence from OddsPath likelihood ratios, plus predictive evidence
calibrated genome-wide only" (see
`notebooks/analysis/README_OddsPath_classifications.md`), rather than
Supplementary Data 5's ExCALIBR/gene-specific-with-fallback classification.
No figure in this repository currently visualizes these OP sheets, so this
script (and the missing-category handling below) is new rather than a port
of an existing chart.

The Unobserved group's MutPred2 sheet is spelled `Unobserved_mut_OP` (not
`Unobserved_MP2_OP`) in Supplementary_Data_6.xlsx -- a naming quirk already
present in `OddsPath_classifications.ipynb`'s own sheet dict, reproduced
here via SHEET_NAME_OVERRIDES rather than "fixed", so this script keeps
reading whatever that notebook actually writes. The classification column
itself is `Class_OP_MP2` regardless -- only the sheet name differs.

Unlike Supplementary Data 5's ExCALIBR/GeneSpecific classification, no
VUS/gnomAD/Unobserved row in Supplementary_Data_6.xlsx ever reaches the
"Pathogenic" bucket (genome-wide REVEL/AM/MP2 top out at PP3_Strong = +4,
and few if any variants in these already-uncertain categories also draw
the full +8 PS3_very_strong from OddsPath) -- so, unlike
`make_extended_data_figure_7alt.py`, a category absent from a sheet is
plotted as a genuine zero rather than treated as a schema error; only a
category outside the five expected labels raises.

--consequence-filter missense restricts every sheet to simplified_consequence
== "missense_variant" rows (same idea and column as
`make_extended_data_figure_7alt.py`'s own option) and defaults the output
path to extended_data_figure_7alt_op_missense/new_classification_heatmap_op_missense.pdf
instead of extended_data_figure_7alt_op/new_classification_heatmap_op.pdf, so it
doesn't overwrite the all-consequences run.
"""

from pathlib import Path

import click
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import colormaps
from matplotlib import colors as mcolors
from matplotlib.patches import Rectangle

DEFAULT_INPUT = Path("data/output/supplementary_data/Supplementary_Data_6.xlsx")

# "all" keeps every consequence type; any other key filters each sheet to
# simplified_consequence values in the given list before counting.
CONSEQUENCE_FILTERS = {
    "all": None,
    "missense": ["missense_variant"],
}


def default_output_path(consequence_filter):
    suffix = "" if consequence_filter == "all" else f"_{consequence_filter}"
    dirname = f"extended_data_figure_7alt_op{suffix}"
    filename = f"new_classification_heatmap_op{suffix}"
    return Path("data/output/figures") / dirname / f"{filename}.pdf"


GROUPS = ["VUS", "gnomAD", "Unobserved"]
PREDICTORS = ["REVEL", "AM", "MP2"]
CATEGORY_ORDER = ["Pathogenic", "Likely Pathogenic", "Uncertain", "Likely Benign", "Benign"]
CATEGORY_DISPLAY = {
    "Pathogenic": "P",
    "Likely Pathogenic": "LP",
    "Uncertain": "VUS",
    "Likely Benign": "LB",
    "Benign": "B",
}

# Supplementary_Data_6.xlsx sheet-naming quirk: every group_predictor sheet
# is named f"{group}_{predictor}_OP" except Unobserved's MutPred2 sheet,
# which is "Unobserved_mut_OP". The classification column is Class_OP_MP2
# in both cases -- only the sheet name is affected.
SHEET_NAME_OVERRIDES = {("Unobserved", "MP2"): "Unobserved_mut_OP"}

FONT_FAMILY = "Arial"
FONT_PT = 6
PCT_DECIMALS = 1
# Monochromatic (single-hue) sequential colormap -- a multi-hue map like
# viridis reads as several distinct colors across only 5-9 cells per row,
# which makes relative magnitude harder to judge at a glance than a single
# light-to-dark ramp.
CMAP_NAME = "Blues"
TOTAL_COL_FILL = "#F2F2F2"
LEGEND_TICKS = [0, 25, 50, 75, 100]
LEGEND_TITLE = "% of row total"

# Same reasoning as this project's R calibrated-chart libraries
# (sankey_calibrated.R/confusion_matrix_calibrated.R): a small, reused gap
# between adjacent layout elements, doubled between the three variant-set
# groups so they read as visually separated blocks rather than one
# undifferentiated 9-row grid.
PAD_MM = 1.2
LABEL_GAP_MM = 1.2
GROUP_GAP_MM = 2.4
MARGIN_MM = 1.2
LEGEND_GAP_MM = 2.4
LEGEND_BAR_H_MM = 2.5
LEGEND_TICK_GAP_MM = 0.8
LEGEND_TITLE_GAP_MM = 0.8


def sheet_name(group, predictor):
    return SHEET_NAME_OVERRIDES.get((group, predictor), f"{group}_{predictor}_OP")


def load_counts(input_path, allowed_consequences=None):
    """Return one row per (group, predictor) with each category's count, `total`
    (= the sheet's own row count after filtering), and `pct_<category>` columns.

    allowed_consequences, if given, restricts each sheet to rows whose
    simplified_consequence is in that list before counting.

    A category from CATEGORY_ORDER absent from a sheet is counted as 0
    (see module docstring: unlike Supplementary Data 5, "Pathogenic" is
    never observed in these OddsPath-based sheets). Raises ValueError only
    if a sheet's Class_OP_<predictor> column contains a category outside
    the five expected ones.
    """
    rows = []
    for group in GROUPS:
        for predictor in PREDICTORS:
            sheet = sheet_name(group, predictor)
            df = pd.read_excel(input_path, sheet_name=sheet)
            if allowed_consequences is not None:
                df = df[df["simplified_consequence"].isin(allowed_consequences)]
            class_col = f"Class_OP_{predictor}"
            counts = df[class_col].value_counts()

            unexpected = set(counts.index) - set(CATEGORY_ORDER)
            if unexpected:
                raise ValueError(f"{sheet}!{class_col} has unexpected categories: {sorted(unexpected)}")

            row = {"group": group, "predictor": predictor, "total": len(df)}
            row.update({category: int(counts.get(category, 0)) for category in CATEGORY_ORDER})
            rows.append(row)

    counts_df = pd.DataFrame(rows)
    for category in CATEGORY_ORDER:
        counts_df[f"pct_{category}"] = counts_df[category] / counts_df["total"] * 100
    return counts_df


def _measure_mm(strings, fontsize_pt=FONT_PT, family=FONT_FAMILY):
    """Return {string: (width_mm, height_mm)}, measuring each string's actual
    rendered extent (multi-line strings included) the same way this project's
    R calibrated-chart libraries measure against a throwaway device."""
    fig = plt.figure(dpi=600)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    sizes = {}
    for s in set(strings):
        text = fig.text(0.5, 0.5, s, fontsize=fontsize_pt, family=family, ha="center", va="center")
        fig.canvas.draw()
        bbox = text.get_window_extent(renderer=renderer)
        sizes[s] = (bbox.width / fig.dpi * 25.4, bbox.height / fig.dpi * 25.4)
        text.remove()
    plt.close(fig)
    return sizes


def _text_color_for(rgba):
    r, g, b = rgba[:3]
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "white" if luminance < 0.5 else "black"


def build_figure(counts_df, output_path):
    counts_df = counts_df.copy()
    for category in CATEGORY_ORDER:
        counts_df[f"text_{category}"] = [
            f"{count:,}\n({pct:.{PCT_DECIMALS}f}%)"
            for count, pct in zip(counts_df[category], counts_df[f"pct_{category}"])
        ]
    counts_df["text_total"] = counts_df["total"].apply(lambda v: f"{v:,}")

    header_labels = [CATEGORY_DISPLAY[category] for category in CATEGORY_ORDER]
    all_cell_texts = [text for category in CATEGORY_ORDER for text in counts_df[f"text_{category}"]]
    all_total_texts = list(counts_df["text_total"])

    legend_tick_labels = [f"{tick}%" for tick in LEGEND_TICKS]

    sizes = _measure_mm(
        set(all_cell_texts)
        | set(all_total_texts)
        | set(header_labels)
        | {"Total"}
        | set(GROUPS)
        | set(PREDICTORS)
        | set(legend_tick_labels)
        | {LEGEND_TITLE}
    )

    cell_w = max(sizes[t][0] for t in all_cell_texts + header_labels) + 2 * PAD_MM
    cell_h = max(sizes[t][1] for t in all_cell_texts) + 2 * PAD_MM
    header_h = max(sizes[h][1] for h in header_labels) + 2 * PAD_MM
    total_col_w = max(sizes[t][0] for t in all_total_texts + ["Total"]) + 2 * PAD_MM
    predictor_col_w = max(sizes[p][0] for p in PREDICTORS) + 2 * PAD_MM
    # Rotated 90 degrees, so its rendered width is the *unrotated* text height.
    group_col_w = max(sizes[g][1] for g in GROUPS) + 2 * PAD_MM

    n_categories = len(CATEGORY_ORDER)
    rows_per_group = len(PREDICTORS)
    grid_w = n_categories * cell_w
    grid_h = header_h + len(GROUPS) * rows_per_group * cell_h + (len(GROUPS) - 1) * GROUP_GAP_MM

    x_group_label = MARGIN_MM
    x_predictor_label = x_group_label + group_col_w
    x_grid = x_predictor_label + predictor_col_w + LABEL_GAP_MM
    x_categories = [x_grid + i * cell_w for i in range(n_categories)]
    x_total = x_grid + grid_w + LABEL_GAP_MM

    y_grid_top = MARGIN_MM + header_h
    grid_bottom = MARGIN_MM + grid_h

    legend_tick_label_h = max(sizes[label][1] for label in legend_tick_labels)
    legend_title_h = sizes[LEGEND_TITLE][1]
    legend_bar_top = grid_bottom + LEGEND_GAP_MM
    legend_bar_bottom = legend_bar_top + LEGEND_BAR_H_MM
    legend_tick_label_top = legend_bar_bottom + LEGEND_TICK_GAP_MM
    legend_title_top = legend_tick_label_top + legend_tick_label_h + LEGEND_TITLE_GAP_MM
    legend_bottom = legend_title_top + legend_title_h

    total_width = x_total + total_col_w + MARGIN_MM
    total_height = legend_bottom + MARGIN_MM

    cmap = colormaps[CMAP_NAME]
    norm = mcolors.Normalize(vmin=0, vmax=100)

    fig = plt.figure(figsize=(total_width / 25.4, total_height / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, total_width)
    ax.set_ylim(total_height, 0)
    ax.axis("off")

    text_kwargs = {"fontsize": FONT_PT, "family": FONT_FAMILY, "ha": "center", "va": "center"}

    for i, label in enumerate(header_labels):
        ax.text(x_categories[i] + cell_w / 2, MARGIN_MM + header_h / 2, label, **text_kwargs)
    ax.text(x_total + total_col_w / 2, MARGIN_MM + header_h / 2, "Total", **text_kwargs)

    for group_index, group in enumerate(GROUPS):
        group_top = y_grid_top + group_index * (rows_per_group * cell_h + GROUP_GAP_MM)
        group_height = rows_per_group * cell_h
        ax.text(
            x_group_label + group_col_w / 2,
            group_top + group_height / 2,
            group,
            rotation=90,
            **text_kwargs,
        )

        for predictor_index, predictor in enumerate(PREDICTORS):
            row = counts_df[(counts_df["group"] == group) & (counts_df["predictor"] == predictor)].iloc[0]
            row_top = group_top + predictor_index * cell_h
            row_center = row_top + cell_h / 2

            ax.text(x_predictor_label + predictor_col_w / 2, row_center, predictor, **text_kwargs)

            for i, category in enumerate(CATEGORY_ORDER):
                color = cmap(norm(row[f"pct_{category}"]))
                ax.add_patch(
                    Rectangle(
                        (x_categories[i], row_top), cell_w, cell_h, facecolor=color, edgecolor="white", linewidth=0.5
                    )
                )
                ax.text(
                    x_categories[i] + cell_w / 2,
                    row_center,
                    row[f"text_{category}"],
                    color=_text_color_for(color),
                    **text_kwargs,
                )

            ax.add_patch(
                Rectangle(
                    (x_total, row_top), total_col_w, cell_h, facecolor=TOTAL_COL_FILL, edgecolor="white", linewidth=0.5
                )
            )
            ax.text(x_total + total_col_w / 2, row_center, row["text_total"], **text_kwargs)

    # Legend: a horizontal 0-100% gradient bar under the heatmap grid (not the
    # label/Total columns, which aren't colored), shared by every cell since
    # they're all drawn on this same fixed 0-100% scale.
    gradient = np.linspace(0, 100, 256).reshape(1, -1)
    ax.imshow(
        gradient,
        extent=(x_grid, x_grid + grid_w, legend_bar_bottom, legend_bar_top),
        aspect="auto",
        cmap=cmap,
        norm=norm,
    )
    ax.add_patch(
        Rectangle(
            (x_grid, legend_bar_top),
            grid_w,
            LEGEND_BAR_H_MM,
            facecolor="none",
            edgecolor="black",
            linewidth=0.5,
        )
    )
    for tick, label in zip(LEGEND_TICKS, legend_tick_labels):
        tick_x = x_grid + tick / 100 * grid_w
        ax.text(tick_x, legend_tick_label_top + legend_tick_label_h / 2, label, **text_kwargs)
    ax.text(x_grid + grid_w / 2, legend_title_top + legend_title_h / 2, LEGEND_TITLE, **text_kwargs)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="pdf")
    plt.close(fig)


@click.command(help=__doc__)
@click.option(
    "--input", "input_path", default=DEFAULT_INPUT, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "--consequence-filter",
    "consequence_filter",
    default="all",
    type=click.Choice(sorted(CONSEQUENCE_FILTERS)),
    show_default=True,
)
@click.option(
    "--output",
    "output_path",
    default=None,
    type=click.Path(dir_okay=False, path_type=Path),
    help=f"Default: {default_output_path('all')}, or ..._<consequence-filter>/... if --consequence-filter isn't 'all'.",
)
def main(input_path, consequence_filter, output_path):
    plt.rcParams["pdf.fonttype"] = 42

    if output_path is None:
        output_path = default_output_path(consequence_filter)

    counts_df = load_counts(input_path, allowed_consequences=CONSEQUENCE_FILTERS[consequence_filter])
    build_figure(counts_df, output_path)
    click.echo(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
