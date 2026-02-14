# Figure 6b

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

# Broad gene-phenotype classes

gene_groups_df <- tribble(
  ~Gene, ~Disease, ~D2, ~D3,
  "BAP1", "Cancer", "Cancer", "Cancer",
  "BARD1", "Cancer", "Cancer", "Cancer",
  "BRCA1", "Cancer", "Cancer", "Cancer",
  "BRCA2", "Cancer", "Cancer", "Cancer",
  "CALM3", "Cardiovascular", "Cardio- vascular", "Cardiovascular",
  "CHEK2", "Cancer", "Cancer", "Cancer",
  "G6PD", "Metabolic", "Meta- bolic", "Meta- bolic",
  "GCK", "Metabolic", "Meta- bolic", "Meta- bolic",
  "KCNE1", "Cardiovascular", "Cardio- vascular", "Cardiovascular",
  "KCNH2", "Cardiovascular", "Cardio- vascular", "Cardiovascular",
  "KCNQ4", "Hearing loss", "Hearing loss", "Hearing loss", #Rare disease",
  "MSH2", "Cancer", "Cancer", "Cancer",
  "OTC", "Metabolic", "Meta- bolic", "Meta- bolic",
  "PALB2", "Cancer", "Cancer", "Cancer",
  "PTEN", "Cancer", "Cancer", "Cancer",
  "RAD51C", "Cancer", "Cancer", "Cancer",
  "RAD51D", "Cancer", "Cancer", "Cancer",
  "SCN5A", "Cardiovascular", "Cardio- vascular", "Cardiovascular",
  "TARDBP", "Rare disease", "Rare disease", "Rare disease",
  "TP53", "Cancer", "Cancer", "Cancer",
  "TSC2", "Cancer", "Cancer", "Cancer"
)


points_levels = c(
  "≤ -16", "≤ -15", "≤ -14", "≤ -13", "≤ -12", "≤ -11", "≤ -10", "≤ -9",
  "≤ -8", "≤ -7", "≤ -6", "≤ -5", "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  #"0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4", "≥ +5", "≥ +6", "≥ +7", "≥ +8",
  "≥ +9", "≥ +10", "≥ +11", "≥ +12", "≥ +13", "≥ +14", "≥ +15", "≥ +16"
)

combined_points_plot_df <-
  or_df %>%
  filter(
    Dataset == 'Total points',
    `Cases with variants` > 0,
    Classification %in% points_levels
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    significance = if_else(
      (LogOR_LI > 0) | (LogOR_UI < 0),
      'Significant at 95%',
      'Not significant at 95%'
    ),
    Classification = factor(
      Classification,
      levels = points_levels
    )
  )

combined_points_condensed_plot_df <-
  combined_points_plot_df %>%
  inner_join(
    combined_points_plot_df %>%
      group_by(Gene) %>%
      summarize(keep = any(significance == 'Significant at 95%')) %>%
      filter(keep) %>%
      select(Gene)
  ) %>%
  left_join(gene_groups_df)

combined_points_condensed_plot <- ggplot(
  combined_points_condensed_plot_df,
  aes(
    x = `Odds Ratio`,
    xmin = OR_LI,
    xmax = OR_UI,
    y = Classification,
    fill = significance
  )
) +
  facet_nested(
    cols = vars(D3, Gene),
    scales = 'free',
    labeller = as_labeller(
      function(s) s %>% str_replace_all('_', ' ') %>% str_wrap(width=7)
    )
  ) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  geom_pointrange(position = position_dodge(width=0.5), shape=21) +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks",
  ) +
  scale_fill_manual(
    values = c(
      'Not significant at 95%' = 'white',
      'Significant at 95%' = 'black'),
    breaks = c('Significant at 95%')
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_blank(aes(x = 6))

# Theme adjustments
cfig_theme <- nature_theme +
  theme(
    strip.text = element_text(
      size = 6,
      face = 'bold',
      color = 'black',
      margin = margin(0,0,0,0),
      family = "Arial"
    ),
    strip.clip = 'off',
    panel.spacing = unit(2, 'pt'),
    axis.title.x = element_text(margin = margin(0,0,0,0)),
    legend.position = 'bottom'
  )

fig6b_plot = combined_points_condensed_plot +
  cfig_theme +
  theme(axis.text.x = element_text(angle = 60))

print(fig6b_plot)

ggsave(
  'fig6b.pdf',
  fig6b_plot,
  width = 183, # Max 183
  height = 90,
  units = 'mm',
  device = cairo_pdf)
ggsave(
  'fig6b.svg',
  fig6b_plot,
  width = 183, # Max 183
  height = 75,
  units = 'mm')
