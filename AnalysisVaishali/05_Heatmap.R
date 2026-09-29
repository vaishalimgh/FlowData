# Vaishali Kaushal, 260929
# Heatmap 

# Rows    = clinical variables (conditions)
# Columns = cell types
# Fill    = Log2 fold change (comparison1 vs comparison2, same groupings as the volcano script)
# Cell text = significance stars based on the (unadjusted) p-value, in 03_Volcano_Plots

# Load libraries
library(ggplot2)
library(reshape2) 

repo_root <- tryCatch(system("git rev-parse --show-toplevel", intern = TRUE, ignore.stderr = TRUE),
                      error = function(e) character(0))
if (length(repo_root) == 1 && nzchar(repo_root) && dir.exists(paste0(repo_root, "/AnalysisVaishali"))) {
  setwd(paste0(repo_root, "/AnalysisVaishali"))
}
rm(list = ls()[ls() != "repo_root"])

# Load data 
data_path <- if (file.exists("../Merged_Flow_Data.csv")) "../Merged_Flow_Data.csv" else "Merged_Flow_Data.csv"
data <- read.csv(data_path, check.names = TRUE)  # check.names=TRUE -> spaces become "."

# "Sex assigned at birth" -> "Sex.assigned.at.birth" 
if ("Sex.assigned.at.birth" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Sex.assigned.at.birth"] <- "Sex"
} else {
  colnames(data)[7] <- "Sex"
}
# "Age at enrollment" -> "Age.at.enrollment" 
if ("Age.at.enrollment" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Age.at.enrollment"] <- "Age"
}

# Long column names into short labels
all_cells <- c( "CD45",
                "HSPCs",
                "Pro_B",
                "Pre_Pro_B",
                "B_cells",
                "Early_NK",
                "Mature_NK",
                "Non_classical_monocyte",
                "Classical_monocyte",
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

# Define columns
cd45_col  <- "CD45"
cell_cols <- setdiff(all_cells, "CD45")

# Calculate proportions (cell count / CD45 count), 
prop_data <- data
for (col in cell_cols) {
  prop_data[, col] <- as.numeric(data[, col]) / as.numeric(data[, cd45_col])
}

# Add age and BMI (make binary using median as threshold),
prop_data$Age_Group <- ifelse(prop_data$Age >= median(prop_data$Age, na.rm = TRUE), "Yes", "No")
prop_data$BMI_Group <- ifelse(prop_data$BMI >= median(prop_data$BMI, na.rm = TRUE), "Yes", "No")

# Define conditions and their clean titles 
conditions <- c(
  "Sex",
  "Diabetes",
  "Hypertension",
  "Coronary.artery.disease",
  "Valve.disease",
  "Stroke",
  "Autoimmune.disease",
  "Peripheral.vascular.disease",
  "Atherosclerosis",
  "History.of.cancer",
  "Smoking",
  "Age_Group",
  "BMI_Group"
)

condition_titles <- c(
  "Sex",
  "Diabetes",
  "Hypertension",
  "Coronary Artery Disease",
  "Valve Disease",
  "Stroke",
  "Autoimmune Disease",
  "Peripheral Vascular Disease",
  "Atherosclerosis",
  "History of Cancer",
  "Smoking",
  "Age (High vs Low)",
  "BMI (High vs Low)"
)

# ---- Compute Log2FC + p-value for every (condition x cell type) 
# Uses the exact same group definitions as the volcano_plot()
#   comparison1 = "Yes"/"Female"/"Checked"   (numerator)
#   comparison2 = "No"/"Male"/"Unchecked"    (denominator)

compute_logfc_table <- function(data, cell_list, conditions, condition_titles) {
  
  results <- list()
  
  for (i in seq_along(conditions)) {
    condition <- conditions[i]
    title     <- condition_titles[i]
    
    comparison1 <- subset(data, data[, condition] == "Yes" | data[, condition] == "Female" | data[, condition] == "Checked")
    comparison2 <- subset(data, data[, condition] == "No"  | data[, condition] == "Male"   | data[, condition] == "Unchecked")
    
    n_cell <- length(cell_list)
    logFC_vec <- numeric(n_cell)
    pval_vec  <- numeric(n_cell)
    
    for (j in seq_len(n_cell)) {
      celltype <- cell_list[j]
      
      x <- as.numeric(comparison1[, celltype])
      y <- as.numeric(comparison2[, celltype])
      
      mean1 <- mean(x, na.rm = TRUE)
      mean2 <- mean(y, na.rm = TRUE)
      logFC_vec[j] <- log2(mean1 / mean2)
      
      ttest <- tryCatch(t.test(x, y), error = function(e) NULL)
      pval_vec[j] <- if (!is.null(ttest)) ttest$p.value else NA
    }
    
    adj_pval_vec <- p.adjust(pval_vec, method = "BH")
    
    results[[i]] <- data.frame(
      Condition   = title,
      Cell        = cell_list,
      LogFC       = logFC_vec,
      PValue      = pval_vec,
      AdjPValue   = adj_pval_vec,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, results)
}

heatmap_data <- compute_logfc_table(
  data             = prop_data,
  cell_list        = cell_cols,
  conditions       = conditions,
  condition_titles = condition_titles
)

# Cap extreme / infinite Log2FC values (e.g., when a group mean is 0) purely for color-scale
# purposes so a handful of +/-Inf values don't wash out the rest of the heatmap.
# Set FC_CAP <- NA to disable capping.
FC_CAP <- 5
heatmap_data$LogFC_capped <- pmin(pmax(heatmap_data$LogFC, -FC_CAP), FC_CAP)

# Significance stars based on the RAW p-value 
# Switch to AdjPValue below if we want to flag by adjusted significance instead.
heatmap_data$Sig <- with(heatmap_data, ifelse(is.na(PValue), "",
                                              ifelse(PValue < 0.001, "***",
                                                     ifelse(PValue < 0.01,  "**",
                                                            ifelse(PValue < 0.05,  "*", "")))))

# Keep axis order consistent with the vectors defined above rather than alphabetical
heatmap_data$Condition <- factor(heatmap_data$Condition, levels = rev(condition_titles))
heatmap_data$Cell      <- factor(heatmap_data$Cell, levels = cell_cols)

# ---- Plot heatmap 
max_abs_fc <- max(abs(heatmap_data$LogFC_capped[is.finite(heatmap_data$LogFC_capped)]), na.rm = TRUE)

p <- ggplot(heatmap_data, aes(x = Cell, y = Condition, fill = LogFC_capped)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = Sig), color = "black", size = 4, vjust = 0.75) +
  scale_fill_gradient2(
    low      = "#457B9D",
    mid      = "white",
    high     = "#E63946",
    midpoint = 0,
    limits   = c(-max_abs_fc, max_abs_fc),
    na.value = "grey85",
    name     = "Log2FC"
  ) +
  xlab("Cell Type") +
  ylab("Clinical Variable") +
  ggtitle("Cell Type Proportion Log2FC Across Clinical Variables") +
  theme_minimal() +
  theme(
    plot.title       = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.text.y      = element_text(size = 10),
    panel.grid       = element_blank(),
    legend.title     = element_text(face = "bold")
  )

# ---- Save outputs 
out_dir <- "05_Heatmap"
if (!dir.exists(out_dir)) dir.create(out_dir)

ggsave(file.path(out_dir, "LogFC_Heatmap.pdf"), p, device = "pdf", width = 14, height = 8)
ggsave(file.path(out_dir, "LogFC_Heatmap.png"), p, device = "png", width = 14, height = 8, dpi = 300)
write.csv(heatmap_data[, c("Condition", "Cell", "LogFC", "PValue", "AdjPValue")],
          file.path(out_dir, "LogFC_Heatmap_data.csv"), row.names = FALSE)

message("Saved: ", file.path(out_dir, "LogFC_Heatmap.pdf"))
message("Saved: ", file.path(out_dir, "LogFC_Heatmap.png"))
message("Saved: ", file.path(out_dir, "LogFC_Heatmap_data.csv"))

