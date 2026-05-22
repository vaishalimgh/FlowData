# Peter van Galen, 230504
# Analyze cell type proportions in relation to age and other clinical factors
# This is based on "/Sternum_BM/AnalysisAdrienne/Scripts/20230315 proportion analysis WIP.R"

# Prerequisites
library(readxl)
library(tidyverse)
library(ggforce)
library(cowplot)

# Clean environment
rm(list=ls())


# Load data ---------------------------------------------------------------------------------------

# Read in the counts data
counts_data <- as.data.frame(read_excel("~/DropboxMGB/Projects/Sternum_BM/AnalysisAdrienne/Previous Analysis/20230315_counts_data.xlsx", sheet = 2))

# Format the counts data
rownames(counts_data) <- counts_data$file
counts_data <- counts_data[,-1]
colnames(counts_data) <- gsub("\\|count", "", colnames(counts_data))
rownames(counts_data) <- gsub(".fcs", "", rownames(counts_data))
rownames(counts_data) <- gsub("_1", "", rownames(counts_data))
rownames(counts_data) <- gsub("-1", "", rownames(counts_data))

counts_data$Sample <- rownames(counts_data)

# Read in and format the cohort data
cohort_data <- as.data.frame(read_excel("~/DropboxMGB/Projects/Sternum_BM/AnalysisAdrienne/20230213 Cohort Clinical Data.xlsx"))
rownames(cohort_data) <- cohort_data$`Record ID`
cohort_data <- cohort_data[,-c(1,2)]

for(row in 1:nrow(cohort_data)){
  donor <- rownames(cohort_data)[row]
  new_name <- paste0("SBM", donor)
  rownames(cohort_data)[row] <- new_name
}

# Add an atherosclerosis column to the cohort data
cohort_data$Atherosclerosis <- ifelse(cohort_data$`Primary pre-operative diagnosis - checkboxes (choice=Coronary artery disease)` == "Checked" | 
                                        cohort_data$`Peripheral vascular disease` == "Yes", yes = "Yes", no = "No")

# Remove donor 1223 and re-save and re-load data
cohort_data$Sample <- rownames(cohort_data)
cohort_data[cohort_data$`Sample` == "SBM1223",] <- NA
cohort_data <- na.omit(cohort_data)


# QC ----------------------------------------------------------------------------------------------

# Make a plot of the CD45+ counts
ggplot(counts_data, aes(y = `CD45+`, x= Sample)) + 
  geom_bar(stat = "identity", width = 0.5, fill = "purple") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)) +
  theme(axis.text.y = element_text(size = 16)) +
  labs(y = "Count CD45+", x = "Sample ID") +
  theme(axis.title.x = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 16))# +
  #facet_zoom(ylim = c(0, 15000), zoom.data = ifelse(`CD45+` <= 15000, NA, FALSE), zoom.size = .25)

# Make a chart of the % viability
counts_data$`Percent viability` <- counts_data$`Live Cells`/counts_data$`Single Cells`

ggplot(counts_data, aes(y = `Percent viability`, x = Sample)) + 
  geom_bar(stat = "identity", width = 0.5, fill = "purple") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)) +
  theme(axis.text.y = element_text(size = 16)) +
  labs(y = "% viability", x = "Sample ID") +
  theme(axis.title.x = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 16))


# Myeloid/lymphoid lineage bias with age ----------------------------------------------------------

# Manually add in the cell types to assess
all_celltypes <- c("HSC","MLP", "MPP", "MEP", "B Cell Progenitors", "B Cells", "Early NK",
                   "Mature NK", "Non-Classical Monocyte", "Classical Monocyte", "MDSC-like",
                   "pDC", "cDC", "CD16+ cDC", "CD16- cDC", "ILC", "NKT CD8-", "NKT CD8+",
                   "CD4+ T cell", "Tregs", "Naive CD4+ T Cell", "Central Memory CD4+ T Cell",
                   "Effector CD4+ T Cell", "CD4+ TPex", "CD8+ T Cell", "CD8+ TPex", "Central Memory CD8+ T Cell",
                   "Naive CD8+ T Cell", "Effector CD8+ T Cell", "gd T cell")
setdiff(colnames(counts_data), all_celltypes) # <-- Should all gates be compared for cell type proportions?
all_progenitors <- c("HSC","MLP", "MPP", "MEP", "B Cell Progenitors")
#all_lymphocytes <- c("B Cells", "Early NK", "Mature NK", "NKT CD8-", "NKT CD8+", "CD4+ T cell", "Tregs", "Naive CD4+ T Cell", "Central Memory CD4+ T Cell",
#                     "Effector CD4+ T Cell", "CD4+ TPex", "CD8+ T Cell", "CD8+ TPex", "Central Memory CD8+ T Cell",
#                     "Naive CD8+ T Cell", "Effector CD8+ T Cell", "gd T cell") # <-- should remove overlap (so cells aren't counted double)
all_lymphocytes <- c("TCRab+", "gd T cell",  "B Cells", "Early NK", "Mature NK", "NKT CD8-", "NKT CD8+")
#all_tcells <- c("CD4+ T cell", "Tregs", "Naive CD4+ T Cell", "Central Memory CD4+ T Cell",
#                "Effector CD4+ T Cell", "CD4+ TPex", "CD8+ T Cell", "CD8+ TPex", "Central Memory CD8+ T Cell",
#                "Naive CD8+ T Cell", "Effector CD8+ T Cell", "gd T cell") # <-- should remove overlap (so cells aren't counted double)
all_tcells <- c("TCRab+", "gd T cell")
all_myeloid <- c("Non-Classical Monocyte", "Classical Monocyte", "MDSC-like",
                 "pDC", "cDC") # <-- removed cDC subsets (compared to Adrienne's script)

# Myeloid bias ----------------------------------
prop_myeloid <- data.frame(myeloid_of_cd45 = rowSums(counts_data[,all_myeloid]) / counts_data$`CD45+`)
# Check
all(rownames(cohort_data) == cohort_data$Sample)
# Merge & plot
prop_myeloid <- merge(prop_myeloid, cohort_data, by = "row.names")
p1 <- prop_myeloid %>%
  ggplot(aes(x = `Age at enrollment`, y = myeloid_of_cd45)) +
  geom_point() +
  geom_smooth(method=lm) +
  ylab("Proportion myeloid cells") +
  theme_bw() +
  theme(aspect.ratio = 0.5)
# Significance test
cor.test(prop_myeloid[,"Age at enrollment"], prop_myeloid[,"myeloid_of_cd45"])

# Relatedly
head(prop_myeloid)
prop_myeloid %>% filter(`Age at enrollment` <= 45) %>% .$myeloid_of_cd45 %>% mean*100
prop_myeloid %>% filter(`Age at enrollment` <= 45) %>% .$myeloid_of_cd45 %>% sd*100
prop_myeloid %>% filter(`Age at enrollment` >= 75) %>% .$myeloid_of_cd45 %>% mean*100
prop_myeloid %>% filter(`Age at enrollment` >= 75) %>% .$myeloid_of_cd45 %>% sd*100
t.test(filter(prop_myeloid, `Age at enrollment` <= 45)$myeloid_of_cd45,
       filter(prop_myeloid, `Age at enrollment` >= 75)$myeloid_of_cd45)
cor.test(prop_myeloid[,"BMI"], prop_myeloid[,"myeloid_of_cd45"])
t.test(filter(prop_myeloid, Diabetes == "Yes")$`myeloid_of_cd45`,
       filter(prop_myeloid, Diabetes == "No")$`myeloid_of_cd45`)
cohort_data %>%
  ggplot(aes(x = Diabetes, y = BMI)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1)

# Lymphoid bias ---------------------------------
prop_lymphoid <- data.frame(lymphoid_of_cd45 = rowSums(counts_data[,all_lymphocytes]) / counts_data$`CD45+`)
# Merge data
prop_lymphoid <- merge(prop_lymphoid, cohort_data, by = "row.names")
p2 <- prop_lymphoid %>%
  ggplot(aes(x = `Age at enrollment`, y = lymphoid_of_cd45)) +
  geom_point() +
  geom_smooth(method=lm) +
  ylab("Proportion lymphoid cells") +
  theme_bw() +
  theme(aspect.ratio = 0.5)
cor.test(prop_lymphoid[,"Age at enrollment"], prop_lymphoid[,"lymphoid_of_cd45"])

# Ratio ---------------------------------

# Calculate myeloid / lymphoid ratio & add age group factors
ratio_df <- merge(prop_myeloid, prop_lymphoid, by = c("Sample", "Age at enrollment")) %>%
  mutate(ratio = myeloid_of_cd45/lymphoid_of_cd45,
         age_group = case_when(`Age at enrollment` < 45 ~ "young",
                               `Age at enrollment` >= 45 & `Age at enrollment` <= 65 ~ "middle-age",
                               `Age at enrollment` > 65 ~ "aged")) %>%
  mutate(age_group = factor(age_group, levels = c("young", "middle-age", "aged")))

# Plot Myeloid/Lymphoid ratio
p3 <- ratio_df %>%
  ggplot(aes(x = `Age at enrollment`, y = ratio)) +
  geom_point() +
  geom_smooth(method=lm) +
  ylab("Myeloid/Lymphoid ratio") +
  theme_bw() +
  theme(aspect.ratio = 0.5)
cor.test(ratio_df[,"Age at enrollment"], ratio_df[,"ratio"])

# Box plot of three age groups
p4 <- ratio_df %>%
  ggplot(aes(x = age_group, y = ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  ylab("Myeloid/Lymphoid ratio") +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 1)

t.test(filter(ratio_df, age_group == "young")$ratio, filter(ratio_df, age_group == "aged")$ratio)
t.test(filter(ratio_df, age_group == "young")$ratio, filter(ratio_df, age_group == "middle-age")$ratio)
t.test(filter(ratio_df, age_group == "middle-age")$ratio, filter(ratio_df, age_group == "aged")$ratio)


pdf("230506_MyeloidBias.pdf")
plot_grid(plot_grid(p1, p2, p3, align = "v", ncol = 1), p4, ncol = 2, rel_heights = c(0.33, 0.33, 0.33))
dev.off()


# Progenitors -----------------------------------
# I tried a number of strategies here but did not find a convincing correlation (230520)
prop_prog <- data.frame(prog_of_cd45 = counts_data$`HSC` / counts_data$`CD45+`, row.names = rownames(counts_data))
# Check
all(rownames(cohort_data) == cohort_data$Sample)
# Merge & plot
prop_prog <- merge(prop_prog, cohort_data, by = "row.names")
prop_prog %>%
  ggplot(aes(x = `Age at enrollment`, y = prog_of_cd45, color = `Smoking`)) +
  geom_point() +
  geom_smooth(method=lm) +
  ylab("Proportion prog cells") +
  theme_bw() +
  theme(aspect.ratio = 0.5)
# Significance test
cor.test(prop_prog[,"Age at enrollment"], prop_prog[,"prog_of_cd45"])



# EVERYTHING BELOW IS ADRIENNE'S CODE

# Correlation function
# "propstring" is a string that describes which function you are taking the proportion of
# e.g. "As a proportion of all CD45+ cells"
cor_func <- function(celltype_list, data, dir, propstring){
  # Set the working directory and generate lists that will be populated with
  # significant cell types
  setwd(dir)
  sigcelltypes_age <- list()
  sigcelltypes_bmi <- list()
  
  # For each cell type in the list, run a correlation between age or BMI and the proportion
  # of that cell type
  for(celltype in celltype_list){
    cor_age <- cor.test(as.numeric(data[,celltype]), as.numeric(data$`Age at enrollment`))
    cor_bmi <- cor.test(as.numeric(data[,celltype]), as.numeric(data$BMI))
    
    # If the correlations have a p-value of less than 0.05, print the correlation info
    # Also add the cell type to the list of significant cell types
    if(cor_age$p.value < 0.05) {
      print(celltype)
      print(cor_age)
      sigcelltypes_age <- append(sigcelltypes_age, celltype)
    }
    if(cor_bmi$p.value < 0.05){
      print(celltype)
      print(cor_bmi)
      sigcelltypes_bmi <- append(sigcelltypes_bmi, celltype)
    }
  }
  
  # For the cell types whose proportions are significantly correlated with age
  # Plot the correlation between that cell type and age and save
  for(celltype in sigcelltypes_age){
    plot <- ggplot(data, aes(y= as.numeric(data[,celltype]), x=as.numeric(`Age at enrollment`))) +
      geom_point(size=2)+
      xlab("Age at enrollment")+
      ylab(paste0("Proportion ", celltype, " ", propstring))+ggtitle(celltype)+ theme(plot.title = element_text(size=22, hjust = 0.5))+
      theme(axis.title=element_text(size=14,face="bold"))+theme(text = element_text(size = 12))
    ggsave(paste0(celltype, ".pdf"), plot, device = "pdf", height = 7, width = 7)
  }
  
  # For the cell types whose proportions are significantly correlated with BMI
  # Plot the correlation between that cell type and BMI and save
  for(celltype in sigcelltypes_bmi){
    plot <- ggplot(data, aes(y= as.numeric(data[,celltype]), x=as.numeric(`BMI`))) +
      geom_point(size=2)+
      xlab("BMI")+
      ylab(paste0("Proportion ", celltype, " ", propstring))+ggtitle(celltype)+ theme(plot.title = element_text(size=22, hjust = 0.5))+
      theme(axis.title=element_text(size=14,face="bold"))+theme(text = element_text(size = 12))
    ggsave(paste0(celltype, ".pdf"), plot, device = "pdf", height = 7, width = 7)
  } 
}

# Run the correlation functions
cor_func(all_celltypes, all_data_45, "/Users/addie/Desktop/20230406 Correlations", "as a proportion of all CD45+")
cor_func(all_lymphocytes, all_data_lymph, "/Users/addie/desktop/flow_data2/lymphocyte cors", "as a proportion of all lymphocytes")
cor_func(all_myeloid, all_data_myeloid, "/Users/addie/desktop/flow_data2/myeloid cors", "as a proportion of all myeloid cells")
cor_func(all_progenitors, all_data_prog, "/Users/addie/desktop/flow_data2/progenitor cors", "as a proportion of all progenitors")
cor_func(all_tcells, all_data_tcells, "/Users/addie/desktop/flow_data2/tcell cors", "as a proportion of all T cells")

# Get a list of conditions to examine in the volcano plot
conditions <- colnames(cohort_data[c(1, 4:10, 12:13)])

# Volcano plot generating function
volcano_plot <- function(data, cell_list, conditions, dir){
  
  setwd(dir)
  
  
  # For each clinical variable (or condition), print the condition and make a dataframe
  # the df will be populated with each cell type, the Log2FC, and the p-value
  for(condition in conditions){
    p_list <- list()
    print(condition)
    comparison <- as.data.frame(matrix(NA, ncol = 3, nrow = length(cell_list)))
    colnames(comparison) <- c("Cell", "LogFC", "p-value")
    
    # For each cell type, subset the data to be one of the two binary options
    # These are typically "Yes" and "No", but may also be "Female" and "Male"
    for(celltype in 1:length(cell_list)){
      comparison1 <- subset(data, subset = data[,condition] == "Yes" | data[,condition] == "Female")
      comparison2 <- subset(data, subset = data[,condition] == "No" | data[,condition] == "Male")
      
      # Get the mean proportions for the subsetted data and calculate the Log2FC
      mean1 <- mean(as.numeric(comparison1[,cell_list[celltype]]))
      mean2 <- mean(as.numeric(comparison2[,cell_list[celltype]]))
      
      logFC <- log2(mean1/mean2)
      
      
      # Perform a t test on the proportions of each cell type to get a p-value
      ttest <- t.test(as.numeric(comparison1[,cell_list[celltype]]), as.numeric(comparison2[,cell_list[celltype]]))
      p.value <- ttest$p.value
      
      p_list <- append(p_list, p.value)
      
      
      # Add the cell types, Log2FC, and -log10p-value to a dataframe for each cell type
      comparison[celltype,1] <- cell_list[celltype]
      comparison[celltype,2] <- logFC
      comparison[celltype,3] <- -log(p.value)
      
    }
    
    # For each row in this final data frame, if the log10 p-value is above 1.301 (<0.05)
    # and the log2FC is >1 (meaning a 2-fold change difference), mark the cell type as
    # either increased compared to controls, down compared to controls, or not sig.
    adj.p <- p.adjust(p_list, method = "BH")
    adj.logp <- -log(adj.p)
    comparison$adjusted <- adj.logp
    comparison$p <- p_list
    
    for(row in 1:nrow(comparison)){
      if(comparison$`p-value`[row] > 1.301 & comparison$LogFC[row] > 0.3){
        comparison$color[row] <- "Increased"
      }
      else if(comparison$`p-value`[row]> 1.301 & comparison$LogFC[row] < -0.3){
        comparison$color[row] <- "Decreased"
      }
      else(comparison$color[row] <- "Not Significant")
    }
    
    # Plot the volcano plot
    p <- ggplot(data=comparison, aes(x=LogFC, y=`p-value`, label = Cell, colour = color)) + geom_point()+geom_text(hjust=1, vjust=1)+
      ggtitle(paste(condition))+ theme(plot.title = element_text(size=22, hjust = 0.5))+
      theme(axis.title=element_text(size=14,face="bold"))+theme(text = element_text(size = 16))+
      theme(aspect.ratio = 1) + xlim(max(abs(comparison$LogFC)+0.5*max(comparison$LogFC))*-1, max(abs(comparison$LogFC)+0.5*max(comparison$LogFC)))+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    ggsave(paste0(condition, "_volcano.pdf"), p, device = "pdf", width = 12, height = 10)
    
  }
}

# Run the volcano plot function
volcano_plot(all_data_45, all_celltypes, conditions, "/Users/addie/Desktop/20230428")
volcano_plot(all_data_lymph, all_lymphocytes, conditions, "/Users/addie/Desktop/20230405_NewVolcanoPlots/lymphocytes")
volcano_plot(all_data_myeloid, all_myeloid, conditions, "/Users/addie/Desktop/20230405_NewVolcanoPlots/myeloid")
volcano_plot(all_data_prog, all_progenitors, conditions, "/Users/addie/Desktop/20230405_NewVolcanoPlots/progenitor")
volcano_plot(all_data_tcells, all_tcells, conditions, "/Users/addie/Desktop/20230405_NewVolcanoPlots/t cells")


##################################################################################
# Age/obesity heat map

data_for_heatmap <- all_data_45[,colnames(all_data_45) %in% all_celltypes]

ages <- sort(unique(all_data_45$`Age at enrollment`))
bmis <- sort(unique(all_data_45$BMI))

data_for_heatmap$age_rank <- NA
data_for_heatmap$bmi_rank <- NA

for(row in 1:nrow(all_data_45)){
  for(age in ages){
    data_for_heatmap$bmi_rank[row] <-  which(bmis == all_data_45$BMI[row])
    data_for_heatmap$age_rank[row] <-  which(ages == all_data_45$`Age at enrollment`[row])
  }
}

library(reshape2)
data_for_heatmap2 <- melt(data_for_heatmap, id = c("age_rank", "bmi_rank"), variable.name = "Cell type")
data_for_heatmap3 <- melt(data_for_heatmap2, id = c("Cell type", "value"), variable.name = "Rank category")


colnames(data_for_heatmap3) <- c("Cell type", "proportion", "Rank category", "Rank")

breaks_to_remove <- seq(3, 63, by = 3)
breaks <- c(1:63)

breaks <- breaks[-breaks_to_remove]


heatmap <- ggplot(data = data_for_heatmap3, mapping = aes(x = `Cell type`,
                                                       y = `Rank`,
                                                       fill = `proportion`)) +
  geom_tile(color = "white") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_distiller(values = c(0, 0.05, seq(0.05, 1, 0.05))) +
  facet_grid(`Rank category`~., scales = "free_y")
  
  
heatmap


library(pheatmap)

p <- pheatmap(data_for_heatmap3, 
         annotation_row = `Rank Category`, 
         breaks = breaks)


