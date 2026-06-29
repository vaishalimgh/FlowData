# Read in the merged CSV file
df <- read.csv("/Users/ritikajain/Library/Group Containers/UBF8T346G9.Office/Outlook/Outlook 15 Profiles/Main Profile/Files/S0/2/Attachments/0/Merged_Flow_Data[8068].csv")

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

# Change race to binary columns
# IF YOU WANT TO REMOVE RACE FROM MODELING, COMMENT OUT lines 44-54
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
  "CD8_Tcell" = df$FlowCut.passed.Cells.Single.Cells.Live.Cells.CD45..CD3..CD34..TCRab..T.Cell.CD8..T.Cell.count,
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

setwd("/Users/ritikajain/Desktop/regression_results")

for(n in 1:length(results_list)){
  print(names(results_list[n]))
  write.csv(file = paste0(names(results_list[n]), ".csv"), results_list[[n]]$coefficients)
}
