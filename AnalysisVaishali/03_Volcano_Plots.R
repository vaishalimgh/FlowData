# Vaishali Kaushal and Peter van Galen, 260625
# Script to make volcano plots comparing cell type proportions across clinical conditions

# Load libraries
library(ggplot2)
library(ggrepel)

# Set working directory & clear environment
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisVaishali"))
rm(list = ls())

# Load data.  This file was likely created by merging files within 'Sternum_BM/Sternum_BM_Flow/AnalysisAdrienne/Counts and Fluor Level data/Flow_Cell_Lists' and we will rerun this script using the regenerated Merged_Flow_Data.csv
data <- read.csv("../Merged_Flow_Data.csv")
colnames(data)[7] <- "Sex"

# Long column names into short readable labels for the plots
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


# Define columns
cd45_col  <- "CD45"
cell_cols <- setdiff(all_cells, "CD45")

# Calculate proportions
prop_data <- data
for(col in cell_cols){
  prop_data[, col] <- as.numeric(data[, col]) / as.numeric(data[, cd45_col])
}
# Add age and BMI (make binary using median as threshold)
prop_data$Age_Group <- ifelse(prop_data$Age >= median(prop_data$Age, na.rm = TRUE), "Yes", "No")
prop_data$BMI_Group <- ifelse(prop_data$BMI >= median(prop_data$BMI, na.rm = TRUE), "Yes", "No")


# Define conditions and their clean titles
conditions <- c(
  "Sex",
  "Diabetes",
  "Hypertension",
  "Primary.pre.operative.diagnosis...checkboxes..choice.Coronary.artery.disease.",
  "Primary.pre.operative.diagnosis...checkboxes..choice.Valve.disease.",
  "Stroke",
  "Autoimmune.disease",
  "Peripheral.vascular.disease",
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
  "History of Cancer",
  "Smoking",       
  "Age (High vs Low)",
  "BMI (High vs Low)"
)

# Volcano plot function 
volcano_plot <- function(data, cell_list, conditions, condition_titles){
  
  for(i in 1:length(conditions)){
    condition <- conditions[i]
    title <- condition_titles[i]
    p_list <- list()
    
    comparison <- as.data.frame(matrix(NA, ncol = 3, nrow = length(cell_list)))
    colnames(comparison) <- c("Cell", "LogFC", "p-value")
    
    for(celltype in 1:length(cell_list)){
      comparison1 <- subset(data, data[,condition] == "Yes" | data[,condition] == "Female" | data[,condition] == "Checked")
      comparison2 <- subset(data, data[,condition] == "No"  | data[,condition] == "Male"   | data[,condition] == "Unchecked")
      
      x <- as.numeric(comparison1[, cell_list[celltype]])
      y <- as.numeric(comparison2[, cell_list[celltype]])
      
      mean1 <- mean(as.numeric(comparison1[, cell_list[celltype]]), na.rm = TRUE)
      mean2 <- mean(as.numeric(comparison2[, cell_list[celltype]]), na.rm = TRUE)
      logFC <- log2(mean1 / mean2)
      
      ttest <- tryCatch(
        t.test(x, y),
        error = function(e) NULL
      )
      p.value <- if(!is.null(ttest)) ttest$p.value else NA
      p_list  <- append(p_list, p.value)
      
      
      # Use short label instead of full column name
      comparison[celltype, 1] <- cell_list[celltype]
      comparison[celltype, 2] <- logFC
      comparison[celltype, 3] <- -log10(p.value)
    }
    
    # Adjusted p-values (BH method). We do not actually plot the adjusted value (as of 260819)!
    adj.p    <- p.adjust(p_list, method = "BH")
    adj.logp <- -log10(adj.p)
    comparison$adjusted <- adj.logp
    
    # Color coding (replace `p-value` with adjusted to color by adjusted P-value instead)
    for(i in 1:nrow(comparison)){
      if(!is.na(comparison$`p-value`[i]) & !is.na(comparison$LogFC[i])){
        if(comparison$`p-value`[i] > 1.301 & comparison$LogFC[i] > 0.3){
          comparison$color[i] <- "Increased"
        } else if(comparison$`p-value`[i] > 1.301 & comparison$LogFC[i] < -0.3){
          comparison$color[i] <- "Decreased"
        } else {
          comparison$color[i] <- "Not Significant"
        }
      } else {
        comparison$color[i] <- "Not Significant"
      }
    }
    
    # Only label significant points
    comparison$label <- ifelse(comparison$color != "Not Significant", comparison$Cell, "")
    
    # Plot
    p <- ggplot(data = comparison,
                aes(x = LogFC, y = `p-value`, label = label, colour = color)) + # Replace `p-value` with adjusted to color by adjusted P-value instead
      geom_point(size = 4, show.legend = TRUE) +
      geom_text_repel(size = 5,
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      point.padding = 0.3,
                      segment.color = "grey50",
                      segment.size = 0.3,
                      force = 2,
                      data = subset(comparison, label != ""),
                      show.legend = FALSE) +
      scale_colour_manual(values = c(
        "Increased"       = "#E63946",
        "Decreased"       = "#457B9D",
        "Not Significant" = "grey60"
      )) +
      guides(colour = guide_legend(override.aes = list(size = 4))) +
      xlab("Log2 Fold Change") +
      ylab("-log10 p-value") + # Note: we're plotting raw P-values
      ggtitle(title) +
      geom_hline(yintercept = 1.301, linetype = "dashed", colour = "grey40") +
      geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", colour = "grey40") +
      theme_bw() +
      theme(
        plot.title        = element_text(size = 22, hjust = 0.5, face = "bold"),
        axis.title        = element_text(size = 16, face = "bold"),
        text              = element_text(size = 14),
        aspect.ratio      = 1,
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        legend.title      = element_blank()
      )
    
    ggsave(paste0("03_Volcano_Plots/", title, "_volcano.pdf"), p, device = "pdf", width = 12, height = 10)
    message("Saved: ", title, "_volcano.pdf")
  }
}

# Save volcano plots
volcano_plot(
  data             = prop_data,
  cell_list        = cell_cols,
  conditions       = conditions,
  condition_titles = condition_titles
)

