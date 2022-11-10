### Filter samples unrelated to FSHD based on keywords in metadata

## these datasets contains samples unrelated to FSHD
GSE_list = c("GSE9397", "GSE3307", "GSE26852")
GSE = "GSE26852"

# keywords to match samples related to FSHD
keywords = c("control",
             "FSHD",
             "normal",
             "Fascioscapulohumeral muscular dystrophy") 

## Read Metadata
meta <- read.csv(file = sprintf("meta/%s_meta.csv", GSE))
# View(meta)

## Get related sample using regular expression
pattern <- paste(keywords,collapse="|") # get or statement in regex
idx <- grep(pattern, meta[, 2], ignore.case = T) # index for matching samples
sample_list <-meta[idx, 1]

## Rewrite the metadata file
write.csv(meta, file = sprintf("meta/%s_meta_old.csv", GSE), row.names=FALSE) # write the original file as old
new_meta <- meta[idx,]
write.csv(new_meta, file = sprintf("meta/%s_meta.csv", GSE), row.names=FALSE)

## Read and rewrite expression files
dataset <- read.csv(file = sprintf("data/%s/%s.csv", GSE, GSE), header = T, row.names=1)
dataset <- dataset[, names(dataset) %in% c("Probe_ID", "Gene_Symbol", "Ensembl_Gene_ID", sample_list)]
View(dataset)
write.csv(dataset, file = sprintf("data/%s/%s.csv", GSE, GSE))
