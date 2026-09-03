# Vaishali Kaushal and Peter van Galen, 260709
# Generate correlation plots of cell type proportions. vs continuous variables (age and BMI)

# Setup ------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(ggpubr)
library(patchwork)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisVaishali"))

# Clear environment variables
rm(list = ls())

# Load flow results (we still need to find out how to reproduce this)
df <- read.csv("../Merged_Flow_Data.csv", check.names = F)

# Define column names ----------------------------------------------------------

all_cells <- c( "CD45",
                "HSPCs",
                "Pro_B",
                "Pre_Pro_B",
                "B.cells", 
                "Early.NK",
                "Mature.NK",
                "Non_classical.monocyte",
                "Classical.monocyte",
                "MDSC_like",
                "pDCs",
                "cDCs",
                "ILC",
                "CD8neg_NKT",
                "CD8pos_NKT",
                "CD4_T",
                "Naive_CD4",
                "CM_CD4",
                "Effector_CD4",
                "PD1_CD4",
                "CD4_TPex",
                "Tregs",
                "CD8_T",
                "Naive_CD8",
                "CM_CD8",
                "Effector_CD8",
                "PD1_CD8",
                "CD8_TPex",
                "gd_T")



# Compute cell type proportions -----------------------------------------------

# Select relevant columns
props <- df %>%
  select(
    `Record ID`,
    `Age at enrollment`,
    BMI,
    all_of(all_cells)
  )

# Add cell proportions (of CD45+ cell count)
for (cell in all_cells) {
  new_col <- paste0("prop_", cell)
  props[[new_col]] <- props[[cell]] / props[["CD45"]]
}


# Plot theme -------------------------------------------------------------------
corr_theme <- theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    aspect.ratio = 1
  )


# Plotting function ------------------------------------------------------------

# Example values
# data <- props
# x_var <- "Age at enrollment"
# cell_name <- "Naive CD4+ T Cell"

make_corr_plot <- function(
  data,
  x_var,
  cell_name
) {
  prop_col <- paste0("prop_", cell_name)
  y_label <- paste0(cell_name, " (% of CD45+)")

  plot_data <- data %>%
    select(age = all_of(x_var), prop = all_of(prop_col))

  # Calculate statistics
  r_val <- round(cor(plot_data$age, plot_data$prop), 3)
  pval <- cor.test(plot_data$age, plot_data$prop)$p.value
  plab <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3)))
  annot <- paste0("r = ", r_val, "\n", plab)

  # Color y labe
  cell_name <- if (pval < 0.05) paste(cell_name, "(significant)") else cell_name
  title_col <- if (pval < 0.05) "forestgreen" else "black"

  # Plot
  ggplot(plot_data, aes(x = age, y = prop)) +
    geom_point(colour = "black", size = 1.5, alpha = 0.5) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      colour = "steelblue",
      fill = "steelblue",
      alpha = 0.15,
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = annot,
      hjust = -0.1,
      vjust = 1.3,
      size = 2.8
    ) +
    labs(
      title = cell_name,
      x = x_var,
      y = y_label
    ) +
    corr_theme +
    theme(plot.title = element_text(colour = title_col))
}


# Generate correlation plots ---------------------------------------------------

# Generate vector with cell type names to plot each one
all_celltypes <- setdiff(all_cells, "CD45")

age_plots <- lapply(all_celltypes, function(cell) {
  make_corr_plot(
    data = props,
    x_var = "Age at enrollment",
    cell_name = cell
  )
})
names(age_plots) <- all_celltypes

bmi_plots <- lapply(all_celltypes, function(cell) {
  make_corr_plot(
    data = props,
    x_var = "BMI",
    cell_name = cell
  )
})
names(bmi_plots) <- all_celltypes


# Split by significance --------------------------------------------------------

age_sig_celltypes <- c(
  "ILC",
  "CD8pos_NKT",
  "Naive_CD4",
  "Naive_CD8",
  "Non_classical.monocyte"
)
bmi_sig_celltypes <- c("Effector_CD4")

age_nonsig_celltypes <- setdiff(all_celltypes, age_sig_celltypes)
bmi_nonsig_celltypes <- setdiff(all_celltypes, bmi_sig_celltypes)

age_sig_plots <- age_plots[age_sig_celltypes]
age_nonsig_plots <- age_plots[age_nonsig_celltypes]
bmi_sig_plots <- bmi_plots[bmi_sig_celltypes]
bmi_nonsig_plots <- bmi_plots[bmi_nonsig_celltypes]


# Export PDFs ------------------------------------------------------------------

n_celltypes <- length(all_celltypes)
plots_per_row <- 6
n_rows <- ceiling(n_celltypes / plots_per_row)
page_height <- n_rows * 3.5 + 0.5

# Save all age correlations
pdf(
  "02.1_Correlations_Age.pdf",
  width = plots_per_row * 3.2,
  height = page_height
)

wrap_plots(age_plots, ncol = plots_per_row) +
  plot_annotation(
    title = "Age at Enrollment vs Immune Cell Proportions",
    subtitle = "Each plot shows proportion of cell population relative to total CD45+ cells"
  )

dev.off()

# Save all BMI correlations
pdf(
  "02.2_Correlations_BMI.pdf",
  width = plots_per_row * 3.2,
  height = page_height
)

wrap_plots(bmi_plots, ncol = plots_per_row) +
  plot_annotation(
    title = "BMI vs Immune Cell Proportions",
    subtitle = "Each plot shows proportion of cell population relative to total CD45+ cells"
  )

dev.off()

