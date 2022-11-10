# Modified from GEO_Genlist_heatmap.Rmd

library(GEOquery)

GSE_list <- c("GSE3307", "GSE9397", "GSE10760", "GSE15090", "GSE26852", 
              "GSE36398", "GSE56787", "GSE115650", "GSE140261")
GPL_list <- c("GPL96", "GPL96", "GPL96", "GPL570", "GPL6947", 
              "GPL6244", "GPL16791", "GPL16791", "GPL16791")

## Download series data from GEO ##
GSE = "GSE140261"
data.series = getGEO(GEO = GSE,
                     AnnotGPL = FALSE,
                     getGPL = FALSE)

## Download platform data from GEO and get sample (phenotype) information ##
GPL = "GPL16791"
data.platform = getGEO(GPL)
data.index = match(GPL, sapply(data.series, annotation))
data.p = pData(data.series[[data.index]])

data.expr = exprs(data.series[[data.index]])
common = intersect(colnames(data.expr), rownames(data.p))
m1 = match(common, colnames(data.expr))
m2 = match(common, rownames(data.p))
data.expr = data.expr[, m1]
data.p = data.p[m2, ]

## Write data
write.csv(
  data.expr,
  file = sprintf("data/%s/%s_expression.csv", GSE, GSE),
  quote = F
)

## Write metadata
View(data.series)
# Get the accession
series_matrix_file = sprintf("%s_series_matrix.txt.gz", GSE)
accession <- data.series[[series_matrix_file]]@phenoData@data[["geo_accession"]]
# Fill the patient ID
patient_id <- data.series[["GSE140261_series_matrix.txt.gz"]]@phenoData@data[["patient_id:ch1"]]
# Fill the sample type
sampletype <- data.series[["GSE140261_series_matrix.txt.gz"]]@phenoData@data[["disease state:ch1"]]
# Fill age
# age <- data.series[["GSE36398_series_matrix.txt.gz"]]@phenoData@data[["age (y):ch1"]]
# create dataframe and write to csv
metadata_df = data.frame(accession, sampletype, patient_id, row.names = TRUE)
write.csv(metadata_df,
          file = sprintf("meta/%s_meta.csv", GSE))

## Match gene symbol
probe_id <- rownames(data.expr)
all(probe_id == data.platform@dataTable@table[["ID"]]) # add if statement
idx = match(probe_id, data.platform@dataTable@table[["ID"]])
# Different platforms have different title for gene symbol
View(data.platform)
# reordered_gene_symbol = data.platform@dataTable@table[["Gene Symbol"]][idx]
# The following lines match probe_id for GSE36398, comment out for other datasets
gene_symbol_list <- c()
for (i in 1:length(data.platform@dataTable@table[["gene_assignment"]])) {
  gene_symbol <- strsplit(data.platform@dataTable@table[["gene_assignment"]][i], " // ")[[1]][2]
  gene_symbol_list <- append(gene_symbol_list, gene_symbol)
}
reordered_gene_symbol = gene_symbol_list[idx]

## Covert to dataframe and write .csv
expression_df <- as.data.frame(data.expr)
expression_df <- cbind(reordered_gene_symbol, expression_df)
View(expression_df)
write.csv(
  expression_df,
  file = sprintf("data/%s/%s_symbol.csv", GSE, GSE),
  quote = F
)
