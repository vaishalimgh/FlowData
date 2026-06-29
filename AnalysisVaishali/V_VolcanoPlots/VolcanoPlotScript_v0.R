# Vaishali Kaushal and Peter van Galen, 260625
# Script to make volcano plots comparing cell type proportions across clinical conditions

# Library
library(ggplot2)
library(ggrepel)

# Set working directory & clear environment
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisVaishali"))
rm(list = ls())

# Load Data
data <- read.csv("../Merged_Flow_Data 2.csv")
colnames(data)[7] <- "Sex"

# Define columns
cd45_col <- colnames(data)[28]   # CD45+ total count column
cell_cols <- colnames(data)[29:ncol(data)]  # All immune cell columns

# Calculate proportions
prop_data <- data
for(col in cell_cols){
  prop_data[, col] <- as.numeric(data[, col]) / as.numeric(data[, cd45_col])
}

# Long column names into short readable labels for the plots
short_names <- c(
  "CD34+ Progenitors",
  "CD4-CD56- Progenitors",
  "CD20-CD123- Progenitors",
  "CD14-CD16- Progenitors",
  "CD11b-CD11c- Progenitors",
  "MLP",
  "MPP",
  "Pro-B",
  "Pre-pro-B",
  "CD3- CD34- Cells",
  "B Cells",
  "CD19-CD20- Cells",
  "Early NK",
  "Mature NK",
  "Non-Classical Monocyte",
  "Classical Monocyte",
  "CD14+ HLA-DR- Cells",
  "MDSC-like",
  "Dendritic Cells",
  "pDC",
  "cDC",
  "CD16+ cDC",
  "CD16- cDC",
  "CD14- HLA-DR- Cells",
  "ILC",
  "CD3+ T Cells",
  "TCRab+ T Cells",
  "NKT CD8-",
  "NKT CD8+",
  "T Cells",
  "CD4+ T Cell",
  "Tregs",
  "CD4+ CD197+",
  "Naive CD4+ T Cell",
  "Central Memory CD4+ T Cell",
  "CD4+ CD197-",
  "Effector CD4+ T Cell",
  "CD279+ CD4+ T Cell",
  "CD4+ TPex",
  "CD8+ T Cell",
  "CD279+ CD8+ T Cell",
  "CD8+ TPex",
  "CD197+ CD8+ T Cell",
  "Central Memory CD8+ T Cell",
  "Naive CD8+ T Cell",
  "CD197- CD8+ T Cell",
  "Effector CD8+ T Cell",
  "gd T Cell"
)

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
  "Smoking"
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
  "Smoking"
)
# Volcano plot function 
volcano_plot <- function(data, cell_list, short_labels, conditions, condition_titles, dir){
  setwd(dir)
  
  for(i in 1:length(conditions)){
    condition <- conditions[i]
    title <- condition_titles[i]
    p_list <- list()
    
    comparison <- as.data.frame(matrix(NA, ncol = 3, nrow = length(cell_list)))
    colnames(comparison) <- c("Cell", "LogFC", "p-value")
    
    for(celltype in 1:length(cell_list)){
      comparison1 <- subset(data, data[,condition] == "Yes" | data[,condition] == "Female" | data[,condition] == "Checked")
      comparison2 <- subset(data, data[,condition] == "No"  | data[,condition] == "Male"   | data[,condition] == "Unchecked")
      
      mean1 <- mean(as.numeric(comparison1[, cell_list[celltype]]), na.rm = TRUE)
      mean2 <- mean(as.numeric(comparison2[, cell_list[celltype]]), na.rm = TRUE)
      logFC <- log2(mean1 / mean2)
      
      ttest   <- t.test(as.numeric(comparison1[, cell_list[celltype]]),
                        as.numeric(comparison2[, cell_list[celltype]]))
      p.value <- ttest$p.value
      p_list  <- append(p_list, p.value)
      
      # Use short label instead of full column name
      comparison[celltype, 1] <- short_labels[celltype]
      comparison[celltype, 2] <- logFC
      comparison[celltype, 3] <- -log10(p.value)
    }
    
    # Adjusted p-values (BH method)
    adj.p    <- p.adjust(p_list, method = "BH")
    adj.logp <- -log10(adj.p)
    comparison$adjusted <- adj.logp
    
    # Color coding
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
    
    # Only label significant points (cleaner look)
    comparison$label <- ifelse(comparison$color != "Not Significant", comparison$Cell, "")
    
    # Plot
    p <- ggplot(data = comparison,
                aes(x = LogFC, y = `p-value`, label = label, colour = color)) +
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
      ylab("-log10 p-value") +
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
    
    ggsave(paste0(title, "_volcano_v0.pdf"), p, device = "pdf", width = 12, height = 10)
    message("Saved: ", title, "_volcano_v0.pdf")
  }
}

# UPDATE THIS PATH to your VolcanoPlots folder
volcano_plot(
  data            = prop_data,
  cell_list       = cell_cols,
  short_labels    = short_names,
  conditions      = conditions,
  condition_titles = condition_titles,
  dir             = "V_VolcanoPlots"
)
