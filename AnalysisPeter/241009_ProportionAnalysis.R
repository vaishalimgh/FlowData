# Peter van Galen, 230504
# Analyze cell type proportions in relation to age and other clinical factors
# This is based on "/Sternum_BM/AnalysisAdrienne/Scripts/20230315 proportion analysis WIP.R"

# Prerequisites
library(readxl)
library(tidyverse)
library(ggforce)
library(ggpubr)
library(cowplot)

# Clean environment
rm(list=ls())

setwd("~/DropboxMGB/Projects/Sternum_BM/Sternum_BM_Flow/AnalysisPeter")

# Load data ---------------------------------------------------------------------------------------

# Read in the counts data
counts_data <- read_csv("~/DropboxMGB/Projects/Sternum_BM/Sternum_BM_Flow/AnalysisAdrienne/August 2023/Flow_Cell_Lists/20230814_Counts_Data.csv")

# Format the counts data
counts_data <- counts_data %>% mutate(Sample = gsub(".fcs", "", file), .before = 1) %>%
  select(!file)
colnames(counts_data) <- gsub("\\|count", "", colnames(counts_data))
colnames(counts_data)

# Read in and format the cohort data
cohort_data <- read_excel("~/DropboxMGB/Projects/Sternum_BM/Sternum_BM_Flow/AnalysisAdrienne/August 2023/Flow_Cell_Lists/all_SBM_donor_clinicaldata.xlsx") %>%
  mutate(Sample = paste0("SBM", `Record ID`), .before = 1)

# Add an atherosclerosis column to the cohort data # DOUBLE CHECK THIS MAKES SENSE
cohort_data$Atherosclerosis <- ifelse(cohort_data$`Primary pre-operative diagnosis - checkboxes (choice=Coronary artery disease)` == "Checked" | 
                                        cohort_data$`Peripheral vascular disease` == "Yes", yes = "Yes", no = "No")

# Merge
identical(sort(cohort_data$Sample), sort(counts_data$Sample))
merged_data <- inner_join(cohort_data, counts_data)
#data.frame(colnames(merged_data)) %>% view

# Need to check that this makes sense with the gating scheme
merged_data <- merged_data %>% mutate(all_myeloid =
                                        `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Non-Classical Monocyte` +
                                        `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Classical Monocyte` +
                                        `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells`,
                                      all_lymphoid = 
                                        `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/B Cells` +
                                        `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+` +
                                        `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/gd T cell` +
                                        `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Early NK` +
                                        `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Mature NK`)
merged_data <- merged_data %>%
  mutate(Myeloid = all_myeloid/`FlowCut-passed/Cells/Single Cells/Live Cells/CD45+`,
         Lymphoid = all_lymphoid/`FlowCut-passed/Cells/Single Cells/Live Cells/CD45+`) %>%
  mutate(myeloid_lymphoid_ratio = Myeloid/Lymphoid,
         `Naive CD4 T` =
           `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD4+/CD197+/Naive CD4+ T Cell` /
           `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+`,
         `Naive CD8 T` =
           `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD197+/CD8+ T Cell/Naive CD8+ T Cell` /
           `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+`)

merged_data %>% select(`Age at enrollment`, `Sex assigned at birth`, Myeloid, `Naive CD4 T`, `Naive CD8 T`) %>%
  pivot_longer(cols = -c(`Age at enrollment`, `Sex assigned at birth`), values_to = "Proportion") %>%
  ggplot(aes(x = `Age at enrollment`, y = Proportion, color = `Sex assigned at birth`)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = ..p.label..), method = "pearson", show.legend = FALSE) +
  scale_color_manual(values = c("#9370DB", "#DAA520")) +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(filename = "241009_AgeProportionsSex.pdf", width = 10, height = 3)



# Out of curiosity
merged_data <- merged_data %>% mutate(CD34 = `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34+` /
                                        `FlowCut-passed/Cells/Single Cells/Live Cells/CD45+`)
merged_data %>%
  ggplot(aes(x = `Age at enrollment`, y = CD34)) +
  geom_point()
cor.test(merged_data$CD34, merged_data$`Age at enrollment`) # nope


                                        

