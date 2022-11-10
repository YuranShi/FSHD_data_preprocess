library(readr)
library(dplyr)
library(tibble)

## Read annotations
tx2gene <-
  read.table(
    "human_annotation.csv",
    header = TRUE,
    sep = ",",
    colClasses = c("NULL", "character", "character") # read only gene symbol and gene_id
  )
# Since multiple transcript correspond to one gene, eliminate genes that are not unique
tx2gene <- distinct(tx2gene)
View(tx2gene) # look at the annotation file

GSE_list <- c("GSE140261", "GSE115650", "GSE56787")
for (GSE in GSE_list) {  # Remove the for loop if running a single dataset
  ## Read data into dataframe
  normalized_data <-
    read.csv(file = sprintf("data/%s/%s_raw_counts.csv", GSE,GSE),
             header = T)
  
  ## Disable the version number (after .) on Gene ID
  for (i in 1:length(normalized_data$X)) {
    normalized_data$X[i] <- substr(normalized_data$X[i], 1, 15)
  }
  View(normalized_data)
  
  ## Match the gene names and add new column to dataframe
  idx <- match(normalized_data[, 1], tx2gene[, 1])
  reorder_symbol <- tx2gene[, 2][idx]
  complete_df <- add_column(normalized_data, reorder_symbol, .before = 1)
  names(complete_df)[1] <- "Gene_Symbol" # Change column names
  names(complete_df)[2] <- "Ensembl_Gene_ID"
  
  ## If missing Gene_Symbol, replace with Ensembl_Gene_ID
  for (i in 1:nrow(complete_df)) {
    if (complete_df[i, 1] == "" || is.na(complete_df[i, 1])) {
      complete_df[i, 1] <- complete_df[i, 2]
    }
  }
  
  write.csv(complete_df, file = sprintf("data/%s/%s_raw_counts.csv", GSE, GSE), row.names = F)
}
