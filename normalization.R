# Reference: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

library(DESeq2)

GSE = "GSE56787"

## Read data and metadata into dataframe
data <-
  read.csv(
    file = sprintf("data/%s/%s_raw_counts.csv", GSE, GSE),
    header = T,
    row.names = 1
  )
meta <- read.csv(file = sprintf("meta/%s_meta.csv", GSE),
                 header = T,
                 row.names = 1)

### Check classes of the data we just brought in
class(data)
class(meta)
View(data)
View(meta)

### Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

# Match meta and data orders if needed
if (all(colnames(data) == rownames(meta)) == F) {
  idx <- match(colnames(data), rownames(meta))
  meta$index <- idx
  meta <- meta[meta$index, ]
}

## Create DESeq2Dataset object
dds <-
  DESeqDataSetFromMatrix(countData = data,
                         colData = meta,
                         design = ~ sampletype)
View(counts(dds))

## Generate the normalized counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(
  normalized_counts,
  file = sprintf("data/%s/%s_normalized_from_raw.csv", GSE, GSE),
  col.names = NA
)
