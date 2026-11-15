# Transcript-to-gene table for tximport (rnaseq_analysis_RSEM_unified_genes.Rmd) from the NCBI RefSeq annotation
# GCF_048594095.1-RS_2025_03 (NCBI Datasets genomic.gff, saved as input_data/09_RNASeq/paramormyrops_ncbi_0.4.gff).
# Rscript code/09_RNASeq/transcript_to_gene_mapping.R
library(GenomicFeatures)

root <- system("git rev-parse --show-toplevel", intern = TRUE)
txdb <- makeTxDbFromGFF(file.path(root, "input_data/09_RNASeq/paramormyrops_ncbi_0.4.gff"))
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, columns = "TXNAME", keytype = "GENEID")
tx2gene <- df[, 2:1]
head(tx2gene)
write.csv(tx2gene, file.path(root, "input_data/09_RNASeq/paramormyrops_ncbi_0.4_tx2gene.csv"))
