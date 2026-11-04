# Vaishali Kaushal, 260929
# Heatmap 
# Rows    = clinical variables (conditions)
# Columns = cell types
# Fill    = Log2 fold change (comparison1 vs comparison2, same groupings as the volcano script)
# Cell text = significance stars based on the (unadjusted) p-value, in 03_Volcano_Plots

# Load libraries
library(ggplot2)
library(reshape2)
library(legendry)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisVaishali"))

# Clear environment variables
rm(list = ls())

# Load data
flow_df <- read.csv("../Merged_Flow_Data.csv", check.names = F)
all_cells <- colnames(flow_df)[c(25,26,29:31,33,36:39,41:59)]

# Define columns
cd45_col  <- "CD45"
cell_cols <- setdiff(all_cells, "CD45")

# Calculate proportions (cell count / CD45 count), 
prop_data <- flow_df
for (col in cell_cols) {
  prop_data[, col] <- as.numeric(flow_df[, col]) / as.numeric(flow_df[, cd45_col])
}

# Add age and BMI (make binary using median as threshold),
prop_data$Age_Group <- ifelse(prop_data$Age >= median(prop_data$Age, na.rm = TRUE), "Yes", "No")
prop_data$BMI_Group <- ifelse(prop_data$BMI >= median(prop_data$BMI, na.rm = TRUE), "Yes", "No")

# Define conditions
conditions <- c(
  "Sex",
  "Diabetes",
  "Hypertension",
  "Coronary artery disease",
  "Valve disease",
  "Stroke",
  "Autoimmune disease",
  "Peripheral vascular disease",
  "Atherosclerosis",
  "History of cancer",
  "Smoking",
  "Age_Group",
  "BMI_Group"
)

# ---- Compute Log2FC + p-value for every (condition x cell type) 
# Uses the exact same group definitions as the volcano_plot()
#   comparison1 = "Yes"/"Female"/"Checked"   (numerator)
#   comparison2 = "No"/"Male"/"Unchecked"    (denominator)

compute_logfc_table <- function(df, cell_list, conditions) {
  
  results <- list()
  
  for (i in seq_along(conditions)) {
    condition <- conditions[i]
    
    comparison1 <- subset(df, df[, condition] == "Yes" | df[, condition] == "Female" | df[, condition] == "Checked")
    comparison2 <- subset(df, df[, condition] == "No"  | df[, condition] == "Male"   | df[, condition] == "Unchecked")
    
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
      Condition   = condition,
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
  df         = prop_data,
  cell_list  = cell_cols,
  conditions = conditions
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

# Cluster conditions by similar Log2FC patterns across cell types
mat <- acast(heatmap_data, Condition ~ Cell, value.var = "LogFC_capped")
mat[is.na(mat)] <- 0
cond_order <- hclust(as.dist(1 - cor(t(mat))), method = "average")

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
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  # Orders rows by cond_order and draws the tree between the panel and the row labels
  scale_y_dendro(clust = as.hclust(rev(as.dendrogram(cond_order)))) +
  xlab("Cell Type") +
  ylab("Clinical Variable") +
  ggtitle("Cell Type Proportion Log2FC Across Clinical Variables") +
  theme_minimal() +
  theme(
    plot.title       = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.text.y      = element_text(size = 10),
axis.ticks.y.left = element_line(colour = "black"),
    panel.grid       = element_blank(),
    legend.title     = element_text(face = "bold")
  )

# ---- Save outputs
out_dir <- "05_Heatmap"
ggsave(file.path(out_dir, "LogFC_Heatmap.pdf"), p, device = "pdf", width = 10, height = 5.5)
ggsave(file.path(out_dir, "LogFC_Heatmap.png"), p, device = "png", width = 10, height = 5.5, dpi = 300)
write.csv(heatmap_data[, c("Condition", "Cell", "LogFC", "PValue", "AdjPValue")],
          file.path(out_dir, "LogFC_Heatmap_data.csv"), row.names = FALSE)
