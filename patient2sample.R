library(GEOquery)
library(stringi)

GSE_list <- c("GSE15090", "GSE36398", "GSE56787", "GSE115650", "GSE140261")

## Download series data from GEO ##
GSE <- "GSE15090"
data.series = getGEO(GEO = GSE, AnnotGPL = FALSE, getGPL = FALSE)
View(data.series)

# TODO: Update assignment when change GSE
geo_accession <-
  data.series[[sprintf("%s_series_matrix.txt.gz", GSE)]]@phenoData@data[["geo_accession"]]
# TODO: Update assignment when change GSE
patient_id <-
  data.series[[sprintf("%s_series_matrix.txt.gz", GSE)]]@phenoData@data[["title"]]

# Create sample2patient dataframe containing GSM accession and patient_id
sample2patient_df <-
  as.data.frame(list(Sample = geo_accession, Patient = patient_id))
View(sample2patient_df)

# Get patient characteristics data
patient_df <- read.csv(file = sprintf("patient_characteristics/%s_patient.csv", GSE), row.names = 1)
View(patient_df)

if (GSE == "GSE15090") {
  title <- stri_sub(sample2patient_df[,2], -3, -1)
} else {
  title <- substr(sample2patient_df[,2], 1, 7)
}
rows <- rownames(patient_df)

all(rows %in% title)
idx <- match(rows, title)
reordered_sample <- sample2patient_df[,1][idx]
# Check and remove NA in reordered sample
if (NA %in% reordered_sample) {
  null_index <- which(sapply(reordered_sample, function(reordered_sample) NA %in% reordered_sample))
  reordered_sample <- reordered_sample[-null_index]
  patient_df <- patient_df[-null_index, ]
}
rownames(patient_df) <- reordered_sample

write.csv(
  patient_df,
  file = sprintf("patient_characteristics_GSM/%s_patient.csv", GSE)
)

