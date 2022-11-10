# Reference: https://hbctraining.github.io/DGE_workshop_salmon/lessons/genomic_annotation.html

# Load libraries
library(AnnotationHub)
library(ensembldb)
library(dplyr)

# Connect to AnnotationHub
ah <- AnnotationHub()

# Query AnnotationHub
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
human_ens <- human_ens[["AH98047"]]

# Create a transcript dataframe
txdb <- transcripts(human_ens, return.type = "data.frame") %>%
  dplyr::select(tx_id, gene_id)
txdb <- txdb[grep("ENST", txdb$tx_id),]

# Create a gene-level dataframe
genedb <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, symbol)

# Merge the two dataframes together
annotations <- inner_join(txdb, genedb)

# Write to annotations dataframe to csv
write.csv(annotations, "human_annotation.csv", row.names=FALSE)
