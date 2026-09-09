# Vaishali Kaushal and Peter van Galen, 260722
# Run linear regression to determine effect size

# Setup ------------------------------------------------------------------------

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisVaishali"))

# Clear environment variables
rm(list = ls())

# Read in the merged CSV file
df <- read.csv("../Merged_Flow_Data 2.csv")


# Linear regression ------------------------------------------------------------

# Make a new data frame of just the clinical covariates
clin_df <- df[,1:23]

# Change the column names by getting rid of extraneous text
colnames(clin_df) <- gsub("Primary.pre.operative.diagnosis...checkboxes..choice.", "", colnames(clin_df))

# Remove columns that will not be used in the regression
clin_df$Record.ID <- NULL # Okay to do because the flow data is in the same order as the clinical data in merged df
clin_df$Procedure.s....notes <- NULL
clin_df$Cancer.type..years.since.treatment..other <- NULL
clin_df$Ethnicity <- NULL
clin_df$Current.smoker <- NULL
clin_df$Endocarditis. <- NULL # All are "unchecked," nothing to model

# Identify the columns that need to be binarized to 1's and 0's
lapply(clin_df, unique)

# Change sex column to binary numbers
clin_df$Sex.assigned.at.birth[clin_df$Sex.assigned.at.birth == "Male"] <- 0
clin_df$Sex.assigned.at.birth[clin_df$Sex.assigned.at.birth == "Female"] <- 1

# Change checks to 1 and unchecks to 0
clin_df[clin_df == "Unchecked"] <- 0
clin_df[clin_df == "Checked"] <- 1

# Change yes/no to 1/0
clin_df[clin_df == "No"] <- 0
clin_df[clin_df == "Yes"] <- 1

# Change race to binary columns or remove it (see below)
#clin_df$White <- 0
#clin_df$White[clin_df$Race == "White"] <- 1

#clin_df$Asian <- 0
#clin_df$Asian[clin_df$Race == "Asian"] <- 1

#clin_df$Multi <- 0
#clin_df$Multi[clin_df$Race == "More than One Race"] <- 1

#clin_df$Black <- 0
#clin_df$Black[clin_df$Race == "Black or African American"] <- 1

# Remove the race column
clin_df$Race <- NULL

# Check data classes
lapply(clin_df, class)

# Change everything to numeric because some columns were character class
clin_df <- as.data.frame(sapply(clin_df, as.numeric))

# Make a new data frame of the cell types we want
prop_data <- data.frame(
  "Pro_B" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD4.CD56..CD20.CD123..CD14.CD16..CD11b.CD11c..CD34.CD38..Pro.B.count,
  "Pre_Pro_B" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD4.CD56..CD20.CD123..CD14.CD16..CD11b.CD11c..CD34.CD38..Pre.pro.B.count,
  "B.cells" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..B.Cells.count,
  "Eary.NK" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..Early.NK.count,
  "Mature.NK" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..Mature.NK.count,
  "Non_classical.monocyte" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..Non.Classical.Monocyte.count,
  "Classical.monocyte" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..Classical.Monocyte.count,
  "MDSC_like" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..CD14..HLA.DR..MDSC.like.count,
  "DCs" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..Dendritic.Cells.count,
  "pDCs" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..Dendritic.Cells.pDC.count,
  "cDCs" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..Dendritic.Cells.cDC.count,
  "CD16pos_cDC" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..Dendritic.Cells.cDC.CD16..cDC.count,
  "CD16neg_cDC" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..Dendritic.Cells.cDC.CD16..cDC.count.1,
  "ILC" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..CD19..CD20..CD14..HLA.DR..ILC.count,
  "CD8neg_NKT" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..NKT.CD8..count,
  "CD8pos_NKT" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..NKT.CD8..count.1,
  "T_cell" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.count,
  "CD4_T" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD4..T.cell.count,
  "Tregs" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD4..T.cell.Tregs.count,
  "Naive_CD4" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD4..T.cell.CD4..CD197..Naive.CD4..T.Cell.count,
  "CM_CD4" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD4..T.cell.CD4..CD197..Central.Memory.CD4..T.Cell.count,
  "Effector_CD4" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD4..T.cell.CD4..CD197..Effector.CD4..T.Cell.count,
  "CD4_TPex" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD4..T.cell.CD279..CD4..T.Cell.CD4..TPex.count,
  "CD8_T" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD8..T.Cell.count,
  "CD8_TPex" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD8..T.Cell.CD279..CD8..T.Cell.CD8..TPex.count,
  "CM_CD8" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD8..T.Cell.CD197..CD8..T.Cell.Central.Memory.CD8..T.Cell.count,
  "Naive_CD8" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD8..T.Cell.CD197..CD8..T.Cell.Naive.CD8..T.Cell.count,
  "Effector_CD8" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD8..T.Cell.CD197..CD8..T.Cell.Effector.CD8..T.Cell.count,
  "gd_T" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..gd.T.cell.count,
  "progenitors" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..count
)

# Initial check of cell type effects in multiple regression

# Initialize an empty list to populate with cell type regression results
results_list <- list()

for(i in 1:ncol(prop_data)){ # For each column index:
  
  flow_data <- data.frame("cell_type" = prop_data[,i]) # Get the normalized counts data for that column
  flow_data <- cbind(flow_data, clin_df) # Add the clinical data
  
  fit <- lm(cell_type ~ ., data = flow_data) # Run linear regression
  
  results_list[[i]] <- summary(fit) # add it to the list
 
}

names(results_list) <- colnames(prop_data) # Add the cell type labels to your results

# View results
print(names(results_list[25])) # Change to any number 1:29
results_list[[25]] # Change to same number 1:29

View(results_list[[25]]$coefficients)

for(n in 1:length(results_list)){
  print(names(results_list[n]))
  write.csv(results_list[[n]]$coefficients,
            file = file.path("03.1_Regression_Results", paste0(names(results_list[n]), ".csv")))
}


# Generate plots ---------------------------------------------------------------


# Make Effect size box plots from regression results

# Prerequisite packages
library(dplyr)
library(ggplot2)

# 1. READ IN DATA
# First, set your working directory to where the regression results are stored and get the filenames

#setwd("<YOUR DROPBOX PATH>/Sternum_BM/Sternum_BM_Flow/AnalysisVaishali/Linear Regression/Regression Results 1_14_2026/regression_results")
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(file.path(repo_root, "AnalysisVaishali/03.1_Regression_Results"))
csv_files <- list.files(".")

#merged_data <- read.csv("<YOUR DROPBOX PATH>/Sternum_BM/Sternum_BM_Flow/AnalysisVaishali/Merged_Flow_Data 2.csv")
merged_data <- read.csv(file.path(repo_root, "Merged_Flow_Data 2.csv"))
# Make a new data frame of just the clinical covariates
clin_df <- merged_data[,1:23]

# Remove columns that will not be used in the regression
clin_df$Record.ID <- NULL # Okay to do because the flow data is in the same order as the clinical data in merged df
clin_df$Procedure.s....notes <- NULL
clin_df$Cancer.type..years.since.treatment..other <- NULL
clin_df$Ethnicity <- NULL
clin_df$Current.smoker <- NULL


# Change the column names by getting rid of extraneous text
colnames(clin_df) <- gsub("Primary.pre.operative.diagnosis...checkboxes..choice.", "", colnames(clin_df))

# Identify the columns that need to be binarized to 1's and 0's
for(i in 1:ncol(clin_df)){ # For each column index,
  
  print(colnames(clin_df)[i]) # Print the column name
  
  print(unique(clin_df[,i])) # Print its unique values
}

# change sex column to binary numbers
clin_df$Sex.assigned.at.birth[clin_df$Sex.assigned.at.birth == "Male"] <- 0
clin_df$Sex.assigned.at.birth[clin_df$Sex.assigned.at.birth == "Female"] <- 1


# Change checks to 1 and unchecks to 0
clin_df[clin_df == "Unchecked"] <- 0
clin_df[clin_df == "Checked"] <- 1

# Change yes/no to 1/0
clin_df[clin_df == "No"] <- 0
clin_df[clin_df == "Yes"] <- 1

# Change race to binary columns NOTE: THESE ARE COMMENTED OUT BECAUSE THEY WERE NOT USED FOR RREGRESSION
#clin_df$White <- 0
#clin_df$White[clin_df$Race == "White"] <- 1

#clin_df$Asian <- 0
#clin_df$Asian[clin_df$Race == "Asian"] <- 1

#clin_df$Multi <- 0
#clin_df$Multi[clin_df$Race == "More than One Race"] <- 1

#clin_df$Black <- 0
#clin_df$Black[clin_df$Race == "Black or African American"] <- 1

# Get rid of the race column with character strings
clin_df$Race <- NULL


# Check numerics are actually numeric
for(i in 1:ncol(clin_df)){ # For each column
  print(colnames(clin_df)[i]) # Print the column name
  print(class(clin_df[,i])) # Print the data class of the values in that column
}

# Change everything to numeric because some columns were character class
clin_df <- as.data.frame(sapply(clin_df, as.numeric))

# Remove the endocarditis column NOTE: FOR SOME REASON THIS WASN'T WORKING UNLESS I MOVED IT HERE
clin_df$Endocarditis. <- NULL

# -----------------------------------------------------------------------------------------
# 3. PLOTIING AND SAVING
# Okay, now we've got the clinical data exactly how you input it into the regression script
# Now, get the standard deviations for each clinical variable
sds <- apply(clin_df, 2, sd, na.rm = TRUE)

# Then, format and make the variables their own column for merging
sds <- as.data.frame(sds)
sds$X <- rownames(sds)

celltypes <- list() #generate an empty list to make saving easier

# Then loop through each csv file and generate a plot:
for(csv in csv_files){
  print(csv)
  data <- read.csv(csv) # Read in the csv file
  
  celltype <- gsub(".csv", "", csv) # get a character string of the cell type name
  celltypes <- unlist(append(celltypes, celltype)) # Add the character string to your list for saving
  
  data <- right_join(data, sds) # Merge the standard deviations with the csv (right join removes intercept)
  data$Effect.Size <- data$Estimate * data$sds # Calculate Effect Size
  data$p.value <- ifelse(data$Pr...t.. < 0.05, "Significant", "Not Significant") # Make a significance annotation for plotting
  
  # plot
  plt <- ggplot(data, aes(x = reorder(X, Effect.Size), y = Effect.Size, fill = p.value)) +
    geom_col() +
    coord_flip() +
    labs(y = "Effect Size", x = "Clinical Covariate", title = celltype) +
    theme_bw() +
    theme(legend.title = element_blank())
  
  
  assign(celltype, plt) # Save plot to environment
}

# Save PDFs
save_pdfs <- function(string){
  plt <- get(string) # Access the plot from the environment
  ggsave(file.path(repo_root, "AnalysisVaishali/03.2_Regression_Plots", paste0(string, ".pdf")), plt, device = "pdf")}

lapply(celltypes, save_pdfs) # Use custom function to save all plots in one line!
