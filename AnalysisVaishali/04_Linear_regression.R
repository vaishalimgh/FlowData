 # Vaishali Kaushal and Peter van Galen, 260722
# Run linear regression to determine effect size

# Setup ------------------------------------------------------------------------

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisVaishali"))

# Clear environment variables
rm(list = ls())

# Read in the merged CSV file
df <- read.csv("../Merged_Flow_Data.csv")


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
# Define the cell type columns to model (same panel as the volcano plot script)
all_cells <- c("CD45",
               "HSPCs", "Pro_B", "Pre_Pro_B", "B.cells", "Early.NK", "Mature.NK",
               "Non_classical.monocyte", "Classical.monocyte", "MDSC_like", "pDCs", "cDCs",
               "ILC", "CD8neg_NKT", "CD8pos_NKT", "CD4_T", "Naive_CD4", "CM_CD4",
               "Effector_CD4", "PD1_CD4", "CD4_TPex", "Tregs", "CD8_T", "Naive_CD8",
               "CM_CD8", "Effector_CD8", "PD1_CD8", "CD8_TPex", "gd_T")

missing_cols <- setdiff(all_cells, colnames(df))
if(length(missing_cols) > 0){
  stop("These expected columns are missing from the CSV: ", paste(missing_cols, collapse = ", "))
}

cell_cols <- setdiff(all_cells, "CD45")
prop_data <- df[, cell_cols]

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
            file = file.path("04.1_Regression_Results", paste0(names(results_list[n]), ".csv")))
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
dir.create(file.path(repo_root, "AnalysisVaishali/04.1_Regression_Results"), showWarnings = FALSE, recursive = TRUE)
setwd(file.path(repo_root, "AnalysisVaishali/04.1_Regression_Results"))
csv_files <- list.files(".")

#merged_data <- read.csv("<YOUR DROPBOX PATH>/Sternum_BM/Sternum_BM_Flow/AnalysisVaishali/Merged_Flow_Data 2.csv")
merged_data <- read.csv(file.path(repo_root, "Merged_Flow_Data.csv"))
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
dir.create("04.1_Regression_Results", showWarnings = FALSE)
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
dir.create(file.path(repo_root, "AnalysisVaishali/04.2_Regression_Plots"), showWarnings = FALSE, recursive = TRUE)
save_pdfs <- function(string){
  plt <- get(string) # Access the plot from the environment
  ggsave(file.path(repo_root, "AnalysisVaishali/04.2_Regression_Plots", paste0(string, ".pdf")), plt, device = "pdf")}

lapply(celltypes, save_pdfs) # Use custom function to save all plots in one line!
