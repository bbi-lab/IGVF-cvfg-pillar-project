# Figure 2i

# This script requires IGVFFI3804AVJR.csv.gz which is available from
# https://data.igvf.org/tabular-files/IGVFFI3804AVJR/

# Load libraries
library(tidyverse)
library(patchwork)
library(ggh4x)
library(extrafont)

# Font setup
loadfonts(device = 'all')

# Figure theme
point_in_mm = 0.3527778
nature_theme <- theme_linedraw() +
  theme(
    text = element_text(family = 'Arial', size = 7),
    line = element_line(linewidth = 1*point_in_mm),
    geom = element_geom(
      linewidth = 1*point_in_mm,
      borderwidth = 1*point_in_mm,
      pointsize = 3*point_in_mm
    ),
    #spacing = unit(3, 'points'),
    axis.ticks = element_line(linewidth = 0.5*point_in_mm),
    #axis.title = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, face = 'bold', color = 'black', margin = margin(3,3,3,3)),
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    ggh4x.facet.nestline = element_line(linewidth = 1*point_in_mm, color="black"),
    plot.tag = element_text(face = 'bold')
  )

# Load main table
or_df <- read_csv("IGVFFI3804AVJR.csv.gz")

# Filter to IGVF functional assays
func_df <- or_df %>%
  filter(
    str_ends(Dataset, '_IGVF'),
    Dataset != 'TSC2_IGVF', # We show TSC2 RapGAP separately in this figure
    Classifier == 'StandardizedClass',
    Classification %in% c('NORMAL', 'ABNORMAL'),
    `Cases with variants` > 0
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    Classification = factor(
      str_c('Functionally ', str_to_title(Classification)),
      levels = c('Functionally Normal', 'Functionally Abnormal')
    ),
    Gene = factor(
      case_when(
        Dataset == "TSC2_rapgap_IGVF" ~ "TSC2 RapGAP",
        .default = Gene
      ),
      levels = c('TSC2 RapGAP', 'BARD1', 'PALB2', 'RAD51D', 'XRCC2', 'CTCF', 'SFPQ')
    )
  )

# Calc limits for small plots
limits_df <- func_df %>% 
  filter(Gene != "TSC2 RapGAP") %>%
  summarise(
    OR_LI = min(OR_LI),
    OR_UI = max(OR_UI)
  )

# Make plot
make_func_plot <- function(plot_df, x_limits)
  ggplot(
    plot_df,
    aes(
      x=`Odds Ratio`,
      xmin=OR_LI,
      xmax=OR_UI,
      y=Classification
    )
  ) +
  facet_grid(cols = vars(Gene), scales = 'free_x') +
  geom_pointrange() +
  geom_errorbar(width = 0.2) +
  scale_x_log10(
    #labels = scales::label_number(accuracy = 0.1),
    minor_breaks = NULL,
    #minor_breaks = scales::minor_breaks_log(),
    guide = "axis_logticks"
  ) +
  geom_blank(aes(x = x_limits)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
  geom_vline(xintercept = 1, linetype = 'dashed')

fig2i_plot <- make_func_plot(func_df, deframe(limits_df))


# Show plot
print(fig2i_plot + nature_theme)


# Save plot
ggsave(
  'fig2i.pdf',
  fig2i_plot + nature_theme,
  width = 100,
  height = 30,
  units = 'mm',
  family = 'Arial',
  device = cairo_pdf
)
