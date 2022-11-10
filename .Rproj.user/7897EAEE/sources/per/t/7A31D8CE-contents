library(readr)
library(dplyr)
library(stringr)
library(tibble)

## Read annotations
tx2gene <-
  read.table(
    "human_annotation.csv",
    header = TRUE,
    sep = ",",
    colClasses = c("NULL", "character", "character") # read only gene symbol and gene_id
  )
# Since multiple transcript correspond to one gene, eliminate genes that're not unique
tx2gene <- distinct(tx2gene)
View(tx2gene) # look at the annotation file

GSE_list <- c("GSE2820","GSE3307","GSE9397","GSE36398","GSE10760","GSE15090","GSE26852")
## Read data into dataframe
GSE <- "GSE36398"
normalized_data <-
  read.csv(file = sprintf("data/%s/%s_symbol.csv", GSE, GSE),
           header = T)
View(normalized_data)

## if multiple gene symbol, keep the first one
for (i in 1:nrow(normalized_data)) {
  if (grepl("///", normalized_data[i, "reordered_gene_symbol"], fixed = TRUE) == TRUE) {
    genes <-
      str_split(normalized_data[i, "reordered_gene_symbol"], " /// ")
    normalized_data[i, "reordered_gene_symbol"] <- genes[[1]][1]
  }
}

## Match the gene names and add new column to dataframe
idx <- match(normalized_data[, 2], tx2gene[, 2])
reoreder_ENST <- tx2gene[, 1][idx]
complete_df <- add_column(normalized_data, reoreder_ENST, .after = 2)
names(complete_df)[1] <- "Probe_ID" # Change column names
names(complete_df)[2] <- "Gene_Symbol"
names(complete_df)[3] <- "Ensembl_Gene_ID"

for (i in 1:nrow(complete_df)) {
  # If missing, replace Gene_Symbol and Ensembl_Gene_ID with probe_id
  if (complete_df[i, 2] == "" || is.na(complete_df[i, 2])) {
    complete_df[i, 2] <- complete_df[i, 1]
    complete_df[i, 3] <- complete_df[i, 1]
  }
  # Replace Ensembl_Gene_ID with probe_ID if no matching to Gene_symbol
  if (complete_df[i, 3] == "" || is.na(complete_df[i, 3])) {
    complete_df[i, 3] <- complete_df[i, 1]
  }
}
# View(complete_df)

write.csv(complete_df, file = sprintf("data/%s/%s.csv", GSE, GSE))
