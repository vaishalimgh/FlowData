# Peter van Galen and Vaishali Kaushal, 260729
# Generate a summary table for the flow data on SBM samples

# Load necessary libraries
library(tidyverse)
library(readxl)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisVaishali"))

# Clear environment variables
rm(list = ls())


# LOAD DATA --------------------------------------------------------------------

# Load clinical variables
clinical_data <- read_excel(
  '../../AnalysisAdrienne/Counts and Fluor Level data/Flow_Cell_Lists/all_SBM_donor_clinicaldata.xlsx'
)

# Load counts data created by Adrienne
counts_data <- read_csv(
  '../../AnalysisAdrienne/Counts and Fluor Level data/Flow_Cell_Lists/20230814_Counts_Data.csv',
  col_types = cols(
    file = col_character(),
    .default = col_integer()
  )
)


# MERGE DATA -------------------------------------------------------------------

# Check that all the IDs match perfectly
counts_data_fcs <- substr(counts_data$file, 1, 7)
clinical_data_id <- paste0("SBM", clinical_data$`Record ID`)
identical(sort(counts_data_fcs), sort(clinical_data_id))

# Wrangle to enable merging
clinical_data <- clinical_data |>
  mutate(`Record ID` = paste0("SBM", `Record ID`))
counts_data <- counts_data |>
  mutate(file = substr(file, 1, 7)) |>
  rename(`Record ID` = file)

# Merge
merged_data <- left_join(clinical_data, counts_data)

# MODIFY COLUMN NAMES ----------------------------------------------------------

# fmt: skip
renaming_table <- tribble(~New, ~Old, ~Include,
  "CD45","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+|count",
  "HSCs","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34+/CD4-CD56-/CD20-CD123-/CD14-CD16-/CD11b-CD11c-/CD34+CD38-|count", # V included
  "Progenitors","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34+/CD4-CD56-/CD20-CD123-/CD14-CD16-/CD11b-CD11c-/CD34+CD38+|count", # V included
  "HSPCs","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34+|count",
  "Pro_B","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34+/CD4-CD56-/CD20-CD123-/CD14-CD16-/CD11b-CD11c-/CD34+CD38+/Pro-B|count",
  "Pre_Pro_B","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34+/CD4-CD56-/CD20-CD123-/CD14-CD16-/CD11b-CD11c-/CD34+CD38+/Pre-pro-B|count",
  "B.cells","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/B Cells|count",
  "Early.NK","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Early NK|count",
  "Mature.NK","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Mature NK|count",
  "Non_classical.monocyte","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Non-Classical Monocyte|count",
  "Classical.monocyte","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Classical Monocyte|count",
  "MDSC_like","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/CD14+/HLA-DR-/MDSC-like|count",
  "DCs","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells|count",
  "pDCs","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells/pDC|count",
  "cDCs","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells/cDC|count",
  "CD16pos_cDC","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells/cDC/CD16+ cDC|count",
  "CD16neg_cDC","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells/cDC/CD16- cDC|count",
  "ILC","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/CD14-/HLA-DR-/ILC|count",
  "CD8neg_NKT","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/NKT CD8-|count",
  "CD8pos_NKT","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/NKT CD8+|count",
  "T_cell","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell|count",
  "CD4_T","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell|count",
  "Tregs","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/Tregs|count",
  "Naive_CD4","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD4+/CD197+/Naive CD4+ T Cell|count",
  "CM_CD4","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD4+/CD197+/Central Memory CD4+ T Cell|count", # V included
  "Effector_CD4","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD4+/CD197-/Effector CD4+ T Cell|count",
  "PD1_CD4","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD279+/CD4+ T Cell|count", # V included
  "CD4_TPex","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD279+/CD4+ T Cell/CD4+ TPex|count",
  "CD8_T","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell|count",
  "Naive_CD8","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD197+/CD8+ T Cell/Naive CD8+ T Cell|count",
  "CM_CD8","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD197+/CD8+ T Cell/Central Memory CD8+ T Cell|count",
  "Effector_CD8","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD197-/CD8+ T Cell/Effector CD8+ T Cell|count",
  "PD1_CD8","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD279+/CD8+ T Cell|count", # V included
  "CD8_TPex","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD279+/CD8+ T Cell/CD8+ TPex|count",
  "gd_T","FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/gd T cell|count",
)

# Check that all the column names are present
renaming_table$Old %in% colnames(merged_data)

# Rename
renaming_vector <- deframe(renaming_table[c("New", "Old")])
merged_data_rename <- merged_data |> rename(all_of(renaming_vector))

# Select relevant columns
merged_data_subset <- merged_data_rename |>
  select(1:23, all_of(renaming_table$New))

# Save
write_csv(merged_data_subset, file = "../Merged_Flow_Data.csv")
