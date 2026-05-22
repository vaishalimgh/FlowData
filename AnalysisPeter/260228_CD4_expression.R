# Peter van Galen, 260228
# Assess CD4 expression in different cell types

# Prerequisites
library(tidyverse)

# Clean environment & set working directory
rm(list = ls())
setwd("~/DropboxMGB/Sternum_BM/Sternum_BM_Flow/AnalysisPeter")

# Load data ---------------------------------------------------------------------------------------

# fmt: skip
progenitors_df <- read_csv("../AnalysisAdrienne/August 2023/Flow_Cell_Lists/20230816 All Progenitors.csv")
# fmt: skip
myeloid_df <- read_csv("../AnalysisAdrienne/August 2023/Flow_Cell_Lists/20230816 All Myeloid.csv")
# fmt: skip
lymphocytes_df <- read_csv("../AnalysisAdrienne/August 2023/Flow_Cell_Lists/20230816 All Lymphocytes.csv")

bind_rows(
  progenitors_df |> mutate(population = "1_Progenitors"),
  myeloid_df |> mutate(population = "2_Myeloid"),
  lymphocytes_df |> mutate(population = "3_Lymphocytes")
) |>
  ggplot(aes(
    x = `cFluor YG584-A___CD4`,
    fill = population
  )) +
  geom_density(alpha = 0.5, color = "black", linewidth = 0.3) +
  coord_cartesian(xlim = c(-10000, 200000)) +
  labs(
    x = "CD4 (cFluor YG584-A)",
    y = "Density",
    fill = "Population"
  ) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank())

ggsave("260228_CD4_expression.pdf", width = 6, height = 6)

# Note to self - for future analyses, 26 objects are found here:
rm(list = ls())
load("../AnalysisAdrienne/August 2023/Flow_Cell_Lists/20230817_FlowData.RData")

# From a quick assessment, it looks like data.seu is the most complete single object — it holds all 750,000 cells with 31 features (fluorescence intensities), the raw values in the counts layer, scaled data, UMAP/PCA reductions, and donor and cell type annotations in data.seu@meta.data.
