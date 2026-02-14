# Extended Data Figure 2

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

# Predictor plot component

# Order levels
vep_levels = c(
  "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  #"0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4"
)

# Select predictor results
predictor_plot_df <- or_df %>%
  filter(
    `Cases with variants` > 0,
    Classifier %in% c(
      'Gene-specific',
      'Genome-wide aggregation'
    ),
    Classification %in% vep_levels
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
      levels = vep_levels
    ),
    Classifier = Classifier %>% fct(
      levels = c(
        'Gene-specific',
        #'Domain aggregation',
        'Genome-wide aggregation'
      )
    )
    ,
    Dataset = fct(
      Dataset,
      levels = c('REVEL', 'AlphaMissense', 'MutPred2')
    )
  )

# Colors
color_mapping <- c(
  `Gene-specific` = "#CC6015", 
  `Domain aggregation` = "#D2A624",
  `Genome-wide aggregation` = "#544439"
)


# Plotting function
make_predictor_plot <- function(
    in_df,
    scales = 'free_x',
    legend_position = 'bottom'
)
  ggplot(
    in_df,
    aes(
      x=`Odds Ratio`,
      xmin=OR_LI,
      xmax=OR_UI,
      y=Classification,
      color = Classifier,
      shape = significance
    )
  ) +
  facet_grid(
    cols = vars(Gene),
    rows = vars(Dataset),
    scales = scales
  ) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  geom_pointrange(position = position_dodge(width=0.5), fill = 'white') +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks"
  ) +
  scale_shape_manual(
    values = c(
      'Not significant at 95%' = 21,
      'Significant at 95%' = 16),
    breaks = c('Significant at 95%'),
    guide = guide_legend(position = legend_position)
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_color_manual(
    values = color_mapping,
    guide = guide_legend(
      position = legend_position,
      override.aes = aes(shape = 21, fill = 'white')
    )
  )

# Make plot frame
figure_predictor_plot_df <- predictor_plot_df %>%
  filter(
    Gene %in% c('BRCA1', 'BRCA2', 'MSH2', 'TP53', 'TSC2'), # Have gene-specific calibration
  )

# Make plot for final figure
figure_predictor_plot <- figure_predictor_plot_df %>%
  make_predictor_plot(legend_position = 'top') +
  geom_blank(aes(x=deframe(
    predictor_plot_df %>%
      filter(Gene == 'MSH2') %>%
      summarize(max(OR_UI))
  )))


# Assay plot component

# Select assays for plot
assay_classification_levels = c(
  "NORMAL", "ABNORMAL",
  "≤ -8", "≤ -7", "≤ -6", "≤ -5", "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  "0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4", "≥ +5", "≥ +6", "≥ +7", "≥ +8"
)

assay_plot_df <-
  or_df %>%
  filter(
    Classifier %in% c('ExCALIBR', 'StandardizedClass'),
    `Cases with variants` > 0,
    Classification %in% assay_classification_levels
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
      levels = assay_classification_levels
    )
  )

condensed_assay_plot_df <-
  assay_plot_df %>%
  filter(
    Classifier == 'ExCALIBR',
    Dataset %in% c(
      'BAP1_Waters_2024',
      'BARD1_IGVF',
      #'BRCA1_Adamovich_2022_Cisplatin',
      #'BRCA1_Adamovich_2022_HDR',
      'BRCA1_Findlay_2018',
      'BRCA2_Hu_2024',
      #'BRCA2_Sahu_2025_HDR',
      'GCK_Gersing_2023_complementation',
      'KCNH2_Jiang_2022',
      'KCNQ4_Zheng_2022_current_homozygous',
      'MSH2_Jia_2021',
      'PALB2_IGVF',
      'RAD51C_Olvera-León_2024',
      'RAD51D_IGVF',
      'SCN5A_Ma_2024',
      'TP53_Fayer_2021_meta',
      #'TP53_Boettcher_2019'
      #'TP53_Fortuno_2021_Kato_meta',
      #'TP53_Giacomelli_2018_combined_score'
      'TSC2_IGVF'
    ),
    Classification != "0"
  ) %>%
  left_join(gene_groups_df)

# Limits for most panels
assay_plot_common_limits <- condensed_assay_plot_df %>% 
  filter(!(Gene %in% c("BRCA1", "MSH2", "KCNH2", "TSC2", "GCK"))) %>%
  summarise(
    OR_LI = min(OR_LI),
    OR_UI = max(OR_UI)
  ) %>% deframe()

# Build assay plot
condensed_assay_plot <- ggplot(
  condensed_assay_plot_df,
  aes(
    x=`Odds Ratio`,
    xmin=OR_LI,
    xmax=OR_UI,
    y=Classification,
    shape = significance#,
    #color = Classifier
  )
) +
  facet_nested_wrap(
    facets = vars(Disease, Gene),
    nrow = 3,
    scales = 'free_x',
  ) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  geom_pointrange(position = position_dodge(width=0.5), fill='white') +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks",
  ) +
  scale_shape_manual(
    values = c(
      'Not significant at 95%' = 21,
      'Significant at 95%' = 16),
    breaks = c('Significant at 95%'),
    guide = 'none'
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed')

# Compose figure together
fig_exd2 <- (condensed_assay_plot + nature_theme) +
  (figure_predictor_plot + nature_theme) +
  plot_annotation(tag_levels='a') +
  plot_layout(
    design = c(
      area(1,1,6,6),
      area(7,1,9,6)
    )
  )

# Show figure
print(fig_exd2)

# Save figure
ggsave(
  'Extended Data Figure 2.pdf',
  fig_exd2,
  width = 160, # Max 183
  height = 247,
  units = 'mm',
  device = cairo_pdf)
