library(flextable)
library(dplyr)
library(officer)

df <- read.csv("/Users/ritikajain/Desktop/Merged_Flow_Data 2.csv")

# ── Continuous variables ──────────────────────────────────────────────────────
age_mean <- round(mean(df$Age.at.enrollment, na.rm = TRUE), 1)
age_sd   <- round(sd(df$Age.at.enrollment,   na.rm = TRUE), 1)
bmi_mean <- round(mean(df$BMI, na.rm = TRUE), 1)
bmi_sd   <- round(sd(df$BMI,   na.rm = TRUE), 1)

# ── Binary variables ──────────────────────────────────────────────────────────
sex_female         <- sum(df$Sex.assigned.at.birth == "Female", na.rm = TRUE)
sex_male           <- sum(df$Sex.assigned.at.birth == "Male",   na.rm = TRUE)
smoker_yes         <- sum(df$Smoking == "Yes",              na.rm = TRUE)
smoker_no          <- sum(df$Smoking == "No",               na.rm = TRUE)
diabetes_yes       <- sum(df$Diabetes == "Yes",             na.rm = TRUE)
diabetes_no        <- sum(df$Diabetes == "No",              na.rm = TRUE)
hypertension_yes   <- sum(df$Hypertension == "Yes",         na.rm = TRUE)
hypertension_no    <- sum(df$Hypertension == "No",          na.rm = TRUE)
hyperlipidemia_yes <- sum(df$Hyperlipidemia == "Yes",       na.rm = TRUE)
hyperlipidemia_no  <- sum(df$Hyperlipidemia == "No",        na.rm = TRUE)
hypo_yes           <- sum(df$Hypothyroidism == "Yes",       na.rm = TRUE)
hypo_no            <- sum(df$Hypothyroidism == "No",        na.rm = TRUE)
pvd_yes            <- sum(df$Peripheral.vascular.disease == "Yes", na.rm = TRUE)
pvd_no             <- sum(df$Peripheral.vascular.disease == "No",  na.rm = TRUE)
stroke_yes         <- sum(df$Stroke == "Yes",               na.rm = TRUE)
stroke_no          <- sum(df$Stroke == "No",                na.rm = TRUE)
cancer_yes         <- sum(df$History.of.cancer == "Yes",    na.rm = TRUE)
cancer_no          <- sum(df$History.of.cancer == "No",     na.rm = TRUE)
autoimmune_yes     <- sum(df$Autoimmune.disease == "Yes",   na.rm = TRUE)
autoimmune_no      <- sum(df$Autoimmune.disease == "No",    na.rm = TRUE)
thrombosis_yes     <- sum(df$History.of.thrombosis.pulmonary.embolism == "Yes", na.rm = TRUE)
thrombosis_no      <- sum(df$History.of.thrombosis.pulmonary.embolism == "No",  na.rm = TRUE)
cad_yes    <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Coronary.artery.disease. == "Checked",   na.rm = TRUE)
cad_no     <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Coronary.artery.disease. == "Unchecked", na.rm = TRUE)
valve_yes  <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Valve.disease. == "Checked",   na.rm = TRUE)
valve_no   <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Valve.disease. == "Unchecked", na.rm = TRUE)
hf_yes     <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Heart.failure. == "Checked",   na.rm = TRUE)
hf_no      <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Heart.failure. == "Unchecked", na.rm = TRUE)
endo_yes   <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Endocarditis. == "Checked",   na.rm = TRUE)
endo_no    <- sum(df$Primary.pre.operative.diagnosis...checkboxes..choice.Endocarditis. == "Unchecked", na.rm = TRUE)

# ── Build data frame ──────────────────────────────────────────────────────────
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
  ),
  stringsAsFactors = FALSE
)

# ── Build flextable ───────────────────────────────────────────────────────────
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

# ── Save as PDF via Word ──────────────────────────────────────────────────────
doc <- read_docx() %>%
  body_add_flextable(ft)

# Save Word first, then convert to PDF
docx_path <- "/Users/ritikajain/Desktop/Table1_demographics.docx"
pdf_path  <- "/Users/ritikajain/Desktop/Table1_demographics.pdf"

print(doc, target = docx_path)
cat("Word file saved. Converting to PDF...\n")

# Convert to PDF using LibreOffice (built into Mac via R)
system(paste0(
  "/Applications/LibreOffice.app/Contents/MacOS/soffice ",
  "--headless --convert-to pdf --outdir ",
  "/Users/ritikajain/Desktop/ ",
  docx_path
))

cat("Done! PDF saved to Desktop as Table1_demographics.pdf\n")