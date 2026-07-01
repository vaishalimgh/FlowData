
#  Supplementary Figures: Correlation Plots
#  Clinical variables (Age, BMI) vs all immune cell populations
#  as a proportion of CD45+ cells

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(patchwork)


# --- SET SAVE LOCATION ---
setwd("/Users/ritikajain/Desktop/V_CorrelationPlot")

df <- read.csv("/Users/ritikajain/Desktop/Merged_Flow_Data 2.csv", check.names = FALSE)


# --- STEP 4: Define the CD45+ total column ---
cd45_col <- "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+|count"


# --- STEP 5: Define all cell types and their short display names ---
cell_cols <- c(
  "Pro-B"                    = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34+/CD4-CD56-/CD20-CD123-/CD14-CD16-/CD11b-CD11c-/CD34+CD38+/Pro-B|count",
  "Pre-pro-B"                = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34+/CD4-CD56-/CD20-CD123-/CD14-CD16-/CD11b-CD11c-/CD34+CD38+/Pre-pro-B|count",
  "B Cells"                  = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/B Cells|count",
  "Early NK"                 = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Early NK|count",
  "Mature NK"                = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Mature NK|count",
  "Non-Classical Monocyte"   = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Non-Classical Monocyte|count",
  "Classical Monocyte"       = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Classical Monocyte|count",
  "MDSC-like"                = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/CD14+/HLA-DR-/MDSC-like|count",
  "Dendritic Cells"          = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells|count",
  "pDC"                      = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells/pDC|count",
  "cDC"                      = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells/cDC|count",
  "CD16+ cDC"                = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells/cDC/CD16+ cDC|count",
  "CD16- cDC"                = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/Dendritic Cells/cDC/CD16- cDC|count",
  "ILC"                      = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3-/CD34-/CD19-/CD20-/CD14-/HLA-DR-/ILC|count",
  "NKT CD8-"                 = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/NKT CD8-|count",
  "NKT CD8+"                 = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/NKT CD8+|count",
  "Tregs"                    = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/Tregs|count",
  "Naive CD4+ T Cell"        = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD4+/CD197+/Naive CD4+ T Cell|count",
  "Central Memory CD4+ T Cell" = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD4+/CD197+/Central Memory CD4+ T Cell|count",
  "Effector CD4+ T Cell"     = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD4+/CD197-/Effector CD4+ T Cell|count",
  "CD4+ TPex"                = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD4+ T cell/CD279+/CD4+ T Cell/CD4+ TPex|count",
  "CD8+ TPex"                = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD279+/CD8+ T Cell/CD8+ TPex|count",
  "Central Memory CD8+ T Cell" = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD197+/CD8+ T Cell/Central Memory CD8+ T Cell|count",
  "Naive CD8+ T Cell"        = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD197+/CD8+ T Cell/Naive CD8+ T Cell|count",
  "Effector CD8+ T Cell"     = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD197-/CD8+ T Cell/Effector CD8+ T Cell|count",
  "gd T Cell"                = "FlowCut-passed/Cells/Single Cells/Live Cells/CD45+/CD3+/CD34-/TCRab+/T Cell/CD8+ T Cell/CD197-/CD8+ T Cell/Effector CD8+ T Cell|count"
)


# --- STEP 6: Calculate proportions ---
# For each cell type, divide its count by the CD45+ total
props <- df %>%
  select(`Record ID`, `Age at enrollment`, BMI, all_of(cd45_col), all_of(unname(cell_cols))) %>%
  filter(!is.na(`Age at enrollment`), !is.na(BMI), !!sym(cd45_col) > 0)

# Add proportion columns
for (short_name in names(cell_cols)) {
  full_col <- cell_cols[[short_name]]
  new_col  <- paste0("prop_", short_name)
  props[[new_col]] <- props[[full_col]] / props[[cd45_col]]
}


# --- STEP 7: Shared plot theme (grey background, clean look) ---
corr_theme <- theme_bw() +
  theme(
    panel.background   = element_rect(fill = "grey92", colour = NA),
    panel.grid.major   = element_line(colour = "white", linewidth = 0.4),
    panel.grid.minor   = element_line(colour = "white", linewidth = 0.2),
    panel.border       = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
    plot.title         = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title         = element_text(size = 8),
    axis.text          = element_text(size = 7),
    strip.background   = element_rect(fill = "grey85"),
    strip.text         = element_text(size = 8, face = "bold"),
    plot.margin        = margin(6, 6, 6, 6)
  )


# --- STEP 8: Function to make one scatter plot ---
make_corr_plot <- function(data, x_var, x_label, cell_name, x_breaks, x_limits) {
  
  prop_col <- paste0("prop_", cell_name)
  y_label  <- paste0(cell_name, " (% of CD45+)")
  
  plot_data <- data %>%
    select(x = all_of(x_var), y = all_of(prop_col)) %>%
    filter(!is.na(x), !is.na(y))
  
  # Calculate R2 and p-value manually
  fit     <- lm(y ~ x, data = plot_data)
  r2      <- round(summary(fit)$r.squared, 3)
  pval    <- summary(fit)$coefficients[2, 4]
  plab    <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3)))
  annot   <- paste0("R2 = ", r2, "\n", plab)
  
  ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(colour = "black", size = 1.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, colour = "steelblue", fill = "steelblue", alpha = 0.15, linewidth = 0.8) +
    annotate("text",
             x = -Inf, y = Inf,
             label = annot,
             hjust = -0.1, vjust = 1.3,
             size = 2.8) +
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    labs(
      title = cell_name,
      x     = x_label,
      y     = y_label
    ) +
    corr_theme
}


# --- STEP 9: Define axis settings ---

# Age: breaks at 40, 60, 80; limits slightly beyond data range
age_breaks <- c(40, 60, 80)
age_limits <- c(20, 90)   # adjust if your youngest/oldest patients fall outside this

# BMI: breaks at category boundaries (18.5, 25, 30) + label positions
# Tick marks at 18.5, 25, 30 show the category cut-offs on the x-axis
bmi_breaks <- c(15, 20, 25, 30, 35, 40, 45)
bmi_limits <- c(13, 47)# adjust to match your data range


# --- STEP 10: Generate all plots ---

# Get all cell names from our lookup table
all_cells <- names(cell_cols)

# ---- AGE PLOTS ----
age_plots <- lapply(all_cells, function(cell) {
  make_corr_plot(
    data     = props,
    x_var    = "Age at enrollment",
    x_label  = "Age at enrollment",
    cell_name = cell,
    x_breaks = age_breaks,
    x_limits = age_limits
  )
})
names(age_plots) <- all_cells

# ---- BMI PLOTS ----
bmi_plots <- lapply(all_cells, function(cell) {
  make_corr_plot(
    data     = props,
    x_var    = "BMI",
    x_label  = "BMI",
    cell_name = cell,
    x_breaks = bmi_breaks,
    x_limits = bmi_limits
  )
})
names(bmi_plots) <- all_cells

# --- Significant cell types ---
age_sig_cells  <- c("ILC", "NKT CD8+", "Naive CD4+ T Cell", "Naive CD8+ T Cell", "Non-Classical Monocyte")
bmi_sig_cells  <- c("Effector CD4+ T Cell")

age_nonsig_cells <- setdiff(all_cells, age_sig_cells)
bmi_nonsig_cells <- setdiff(all_cells, bmi_sig_cells)

# --- Split plot lists by significance ---
age_sig_plots    <- age_plots[age_sig_cells]
age_nonsig_plots <- age_plots[age_nonsig_cells]
bmi_sig_plots    <- bmi_plots[bmi_sig_cells]
bmi_nonsig_plots <- bmi_plots[bmi_nonsig_cells]


# --- STEP 11: Export as PDFs ---
# One PDF per clinical variable, with 3 plots per row

n_cells <- length(all_cells)
plots_per_row <- 3
n_rows <- ceiling(n_cells / plots_per_row)
page_height <- n_rows * 3.5 + 0.5   # inches

# -- Age PDF --
pdf("Supplementary_Correlations_Age.pdf",
    width = plots_per_row * 3.2,
    height = page_height)

wrap_plots(age_plots, ncol = plots_per_row) +
  plot_annotation(
    title    = "Supplementary Figure: Age at Enrollment vs Immune Cell Proportions",
    subtitle = "Each plot shows proportion of cell population relative to total CD45+ cells",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

dev.off()
message("Saved: Supplementary_Correlations_Age.pdf")

# -- BMI PDF --
pdf("Supplementary_Correlations_BMI.pdf",
    width = plots_per_row * 3.2,
    height = page_height)

wrap_plots(bmi_plots, ncol = plots_per_row) +
  plot_annotation(
    title    = "Supplementary Figure: BMI vs Immune Cell Proportions",
    subtitle = "Each plot shows proportion of cell population relative to total CD45+ cells\nBMI cut-offs: Underweight <18.5 | Normal 18.5-24.9 | Overweight 25-29.9 | Obese \u226530",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

dev.off()
message("Saved: Supplementary_Correlations_BMI.pdf")

message("All done! Check your working directory for the two PDF files.")

# --- Save significance-split PDFs ---

save_pdf <- function(plot_list, ncols, filename, title_text, subtitle_text) {
  n      <- length(plot_list)
  nrows  <- ceiling(n / ncols)
  pdf(filename, width = ncols * 3.2, height = nrows * 3.5 + 0.8)
  print(
    wrap_plots(plot_list, ncol = ncols) +
      plot_annotation(
        title    = title_text,
        subtitle = subtitle_text,
        theme    = theme(
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey40")
        )
      )
  )
  dev.off()
  message("Saved: ", filename)
}

save_pdf(age_sig_plots,    ncols = 2, "Sig_Age_Correlations.pdf",
         "Supplementary Figure: Significant Associations with Age",
         "Cell populations with p < 0.05 | Proportion of total CD45+ cells")

save_pdf(age_nonsig_plots, ncols = 3, "NonSig_Age_Correlations.pdf",
         "Supplementary Figure: Non-Significant Associations with Age",
         "Cell populations with p >= 0.05 | Proportion of total CD45+ cells")

save_pdf(bmi_sig_plots,    ncols = 2, "Sig_BMI_Correlations.pdf",
         "Supplementary Figure: Significant Associations with BMI",
         "Cell populations with p < 0.05 | Proportion of total CD45+ cells")

save_pdf(bmi_nonsig_plots, ncols = 3, "NonSig_BMI_Correlations.pdf",
         "Supplementary Figure: Non-Significant Associations with BMI",
         "Cell populations with p >= 0.05 | Proportion of total CD45+ cells")

message("All 6 PDFs saved!")