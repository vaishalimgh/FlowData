# Vaishali Kaushal and Peter van Galen, 260625
# Script to make volcano plots comparing cell type proportions across clinical conditions

# Load libraries
library(ggplot2)
library(ggrepel)

# Set working directory & clear environment
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisVaishali"))
rm(list = ls())

# Load flow results & define cell types of interest
df <- read.csv("../Merged_Flow_Data.csv", check.names = F)
all_cells <- colnames(df)[c(25,26,29:31,33,36:39,41:59)]

# Define columns
cd45_col  <- "CD45"
cell_cols <- setdiff(all_cells, "CD45")

# Calculate proportions
prop_data <- df
for(col in cell_cols){
  prop_data[, col] <- as.numeric(df[, col]) / as.numeric(df[, cd45_col])
}
# Add age and BMI (make binary using median as threshold)
prop_data$Age_Group <- ifelse(prop_data$Age >= median(prop_data$Age, na.rm = TRUE), "Yes", "No")
prop_data$BMI_Group <- ifelse(prop_data$BMI >= median(prop_data$BMI, na.rm = TRUE), "Yes", "No")


# Define conditions and their clean titles
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

# Volcano plot function 
volcano_plot <- function(data, cell_list, conditions){
  
  for(i in 1:length(conditions)){
    condition <- conditions[i]
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
      ggtitle(condition) +
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
    
    ggsave(paste0("03.2_Volcano_Plots/", condition, "_volcano.pdf"), p, device = "pdf", width = 12, height = 10)
    message("Saved: ", condition, "_volcano.pdf")
  }
}

# Save volcano plots
volcano_plot(
  data             = prop_data,
  cell_list        = cell_cols,
  conditions       = conditions
)


# Saving the numerical results (csv) 

dir.create("03.1_Volcano_Results", showWarnings = FALSE)
dir.create("03.2_Volcano_Plots",   showWarnings = FALSE)

for(condition in conditions){
  g1 <- subset(prop_data, prop_data[,condition] %in% c("Yes", "Female", "Checked"))
  g2 <- subset(prop_data, prop_data[,condition] %in% c("No",  "Male",   "Unchecked"))
  
  res <- do.call(rbind, lapply(cell_cols, function(cell){
    x <- as.numeric(g1[, cell]); y <- as.numeric(g2[, cell])
    p <- tryCatch(t.test(x, y)$p.value, error = function(e) NA)
    data.frame(Condition   = condition,
               Cell        = cell,
               Mean_Group1 = mean(x, na.rm = TRUE),
               Mean_Group2 = mean(y, na.rm = TRUE),
               LogFC       = log2(mean(x, na.rm = TRUE) / mean(y, na.rm = TRUE)),
               P_value     = p)
  }))
  res$neg_log10_P    <- -log10(res$P_value)
  res$adj_P_BH       <- p.adjust(res$P_value, method = "BH")
  res$neg_log10_adjP <- -log10(res$adj_P_BH)
  res$Category <- ifelse(!is.na(res$neg_log10_P) & res$neg_log10_P > 1.301 & res$LogFC >  0.3, "Increased",
                         ifelse(!is.na(res$neg_log10_P) & res$neg_log10_P > 1.301 & res$LogFC < -0.3, "Decreased",
                                "Not Significant"))
  
  write.csv(res, paste0("03.1_Volcano_Results/", condition, "_volcano_results.csv"), row.names = FALSE)
  message("Saved: ", condition, "_volcano_results.csv")
}

