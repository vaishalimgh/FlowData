# Vaishali Kaushal and Peter van Galen, 260709
# Generate a demographic table for the merged flow data, and save as PDF via Word.

# Load libraries
library(tidyverse)
library(flextable)
library(officer)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisVaishali"))

# Clear environment variables
rm(list = ls())

# Load flow data
df <- read.csv("../Merged_Flow_Data.csv") 

# -- Continuous variables ------------------------------------------------------

age_mean <- round(mean(df$Age.at.enrollment), 1)
age_sd   <- round(sd(df$Age.at.enrollment), 1)
bmi_mean <- round(mean(df$BMI), 1)
bmi_sd   <- round(sd(df$BMI), 1)

# -- Binary variables ----------------------------------------------------------

sex_female         <- sum(df$Sex.assigned.at.birth == "Female")
sex_male           <- sum(df$Sex.assigned.at.birth == "Male")
smoker_yes         <- sum(df$Smoking == "Yes")
smoker_no          <- sum(df$Smoking == "No")
diabetes_yes       <- sum(df$Diabetes == "Yes")
diabetes_no        <- sum(df$Diabetes == "No")
hypertension_yes   <- sum(df$Hypertension == "Yes")
hypertension_no    <- sum(df$Hypertension == "No")
hyperlipidemia_yes <- sum(df$Hyperlipidemia == "Yes")
hyperlipidemia_no  <- sum(df$Hyperlipidemia == "No")
hypo_yes           <- sum(df$Hypothyroidism == "Yes")
hypo_no            <- sum(df$Hypothyroidism == "No")
pvd_yes            <- sum(df$Peripheral.vascular.disease == "Yes")
pvd_no             <- sum(df$Peripheral.vascular.disease == "No")
stroke_yes         <- sum(df$Stroke == "Yes")
stroke_no          <- sum(df$Stroke == "No")
cancer_yes         <- sum(df$History.of.cancer == "Yes")
cancer_no          <- sum(df$History.of.cancer == "No")
autoimmune_yes     <- sum(df$Autoimmune.disease == "Yes")
autoimmune_no      <- sum(df$Autoimmune.disease == "No")
thrombosis_yes     <- sum(df$History.of.thrombosis.pulmonary.embolism == "Yes")
thrombosis_no      <- sum(df$History.of.thrombosis.pulmonary.embolism == "No")
cad_yes    <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Coronary.artery.disease. == "Checked")
cad_no     <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Coronary.artery.disease. == "Unchecked")
valve_yes  <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Valve.disease. == "Checked")
valve_no   <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Valve.disease. == "Unchecked")
hf_yes     <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Heart.failure. == "Checked")
hf_no      <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Heart.failure. == "Unchecked")
endo_yes   <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Endocarditis. == "Checked")
endo_no    <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Endocarditis. == "Unchecked")


# -- Build data frame ----------------------------------------------------------

table1 <- data.frame(
  num = as.character(c(1:17)),
  variable = c(
    "Sex Assigned at Birth",
    "Age at Enrollment",
    "BMI",
    "Smoking",
    "Diabetes",
    "Hypertension",
    "Hyperlipidemia",
    "Hypothyroidism",
    "Peripheral Vascular Disease",
    "Stroke",
    "History of Cancer",
    "Autoimmune Disease",
    "History of Thrombosis / Pulmonary Embolism",
    "Coronary Artery Disease",
    "Valve Disease",
    "Heart Failure",
    "Endocarditis"
  ),
  yes = c(
    paste0(sex_female, " (Female)"),
    paste0(age_mean, " \u00b1 ", age_sd, " years"),
    paste0(bmi_mean, " \u00b1 ", bmi_sd),
    smoker_yes, diabetes_yes, hypertension_yes, hyperlipidemia_yes,
    hypo_yes, pvd_yes, stroke_yes, cancer_yes, autoimmune_yes,
    thrombosis_yes, cad_yes, valve_yes, hf_yes, endo_yes
  ),
  no = c(
    paste0(sex_male, " (Male)"),
    "",   # will be merged
    "",   # will be merged
    smoker_no, diabetes_no, hypertension_no, hyperlipidemia_no,
    hypo_no, pvd_no, stroke_no, cancer_no, autoimmune_no,
    thrombosis_no, cad_no, valve_no, hf_no, endo_no
  )
)

# -- Build flextable -----------------------------------------------------------

ft <- flextable(table1) %>%
  
  # Column headers
  set_header_labels(
    num      = "Variable No.",
    variable = "Variable",
    yes      = "Number \u201cYes\u201d",
    no       = "Number \u201cNo\u201d"
  ) %>%
  
  # TRUE cell merge for Age and BMI rows (rows 2 and 3) across yes+no columns
  merge_at(i = 2, j = 3:4) %>%
  merge_at(i = 3, j = 3:4) %>%
  
  # Font — Times New Roman throughout
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 12, part = "all") %>%
  
  # Header styling: dark blue background, white bold text
  bg(bg = "#2E4057", part = "header") %>%
  color(color = "white", part = "header") %>%
  bold(part = "header") %>%
  
  # Alternating row shading
  bg(i = seq(1, 17, 2), bg = "#EAF2F8", part = "body") %>%
  bg(i = seq(2, 17, 2), bg = "#D6E4F0", part = "body") %>%
  
  # Center-align merged Age/BMI cells
  align(i = 2:3, j = 3, align = "center", part = "body") %>%
  
  # General alignment
  align(j = 1,   align = "center", part = "all") %>%
  align(j = 2,   align = "left",   part = "body") %>%
  align(j = 3:4, align = "center", part = "body") %>%
  align(j = 2:4, align = "center", part = "header") %>%
  
  # Column widths (inches)
  width(j = 1, width = 0.8) %>%
  width(j = 2, width = 3.2) %>%
  width(j = 3, width = 1.6) %>%
  width(j = 4, width = 1.6) %>%
  
  # Row height
  height_all(height = 0.3) %>%
  
  # Clean border styling
  border_outer(part = "all",    fp_border(color = "#2E4057", width = 1.5)) %>%
  border_inner_h(part = "body", fp_border(color = "#AAAAAA", width = 0.5)) %>%
  border_inner_v(part = "all",  fp_border(color = "#2E4057", width = 1))


# Save as Word document --------------------------------------------------------

# Save Word first
doc <- read_docx() %>%
  body_add_flextable(ft)
print(doc, target = "Demographic_Table.docx")
