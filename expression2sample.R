library(GEOquery)

# GSE_list <- c("GSE3307", "GSE9397", "GSE10760", "GSE15090", "GSE26852", 
#               "GSE36398", "GSE56787", "GSE115650", "GSE140261")
# GPL_list <- c("GPL96", "GPL96", "GPL96", "GPL570", "GPL6947", 
#               "GPL6244", "GPL16791", "GPL16791", "GPL16791")

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
View(data.p)

# Get the expression data
expression_df <- read.csv(file = sprintf("data/%s/%s_raw_counts.csv", GSE, GSE), row.names = 1)
View(expression_df)

if (GSE == "GSE115650" || GSE == "GSE140261") {
  title <- substr(data.p[,1], 4, 7)
  # use when processing raw counts
  # columns <- colnames(df) 
  # Use when processing the final expression
  # columns <- colnames(expression_df)[3:length(colnames(expression_df))]
  columns <- substr(columns, 4, 8)
} else if (GSE == "GSE56787"){
  # title <- data.p[,1]
  # title <- sapply(strsplit(title,"_"), `[`, 1) # split based on "_" and take the first element
  # columns <- colnames(expression_df)[3:length(colnames(expression_df))]
  title <- data.p[,1]
  columns <- colnames(expression_df)
}
idx <- match(columns, title)
reordered_sample <- data.p[,2][idx]
expression_df <- setNames(df, reordered_sample)
write.csv(
  expression_df,
  file = sprintf("data_GSM/%s_raw_counts.csv", GSE),
  quote = F
)


