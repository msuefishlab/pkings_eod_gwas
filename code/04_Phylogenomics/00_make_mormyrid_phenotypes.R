#!/usr/bin/env Rscript
# Build the tip table read by phylogeny_and_character_reconstruction_mormyrids.Rmd from Supplementary Data S3.
#
# Inputs:
#   input_data/04_Phylogenomics/SupplementalDataS3_MormyridEODs.xlsx   Supplementary Data S3 of the paper
#   input_data/04_Phylogenomics/ssAll_iqtree_timtree.tree             Osteoglossiformes time tree (Peterson et al. 2022)
# Output:
#   input_data/04_Phylogenomics/mormyrid_specimens_with_phenotypes.csv
#
# Data S3 lists the tips of the Mormyroidea subtree that were assigned a phenotype. Every other tip of the subtree is
# written with phenotype "??" (not assessed) and its genus and species taken from the tip label; the Rmd drops these
# tips before the reconstruction. Tip Number is the tip's index in the subtree, as in Data S3.
# Usage: Rscript code/04_Phylogenomics/00_make_mormyrid_phenotypes.R [DataS3.xlsx] [out.csv] [tree]
suppressPackageStartupMessages({
  library(readxl)
  library(ape)
})

root <- system("git rev-parse --show-toplevel", intern = TRUE)
ind <- file.path(root, "input_data", "04_Phylogenomics")
args <- commandArgs(trailingOnly = TRUE)
s3_path <- if (length(args) >= 1) args[1] else file.path(ind, "SupplementalDataS3_MormyridEODs.xlsx")
out_path <- if (length(args) >= 2) args[2] else file.path(ind, "mormyrid_specimens_with_phenotypes.csv")
tree_path <- if (length(args) >= 3) args[3] else file.path(ind, "ssAll_iqtree_timtree.tree")
if (!file.exists(s3_path)) stop("Supplementary Data S3 not found at ", s3_path, "; download it from the paper")

s3 <- as.data.frame(read_excel(s3_path, sheet = 1))

# The same Mormyroidea subtree as in phylogeny_and_character_reconstruction_mormyrids.Rmd
osteotree <- read.tree(tree_path)$`UTREE1=`
tips <- extract.clade(osteotree, 220, root.edge = 0, collapse.singles = TRUE)$tip.label
stopifnot(all(s3[["Tip Name"]] %in% tips), all(s3[["Tip Number"]] == match(s3[["Tip Name"]], tips)))

na0 <- function(x) ifelse(is.na(x), "", as.character(x))
known <- data.frame(
  "Tip Number" = as.integer(s3[["Tip Number"]]),
  "Tip Name" = s3[["Tip Name"]],
  Genus = s3$Genus,
  Species = s3$Species,
  "Corrected Species Name" = s3[["Species Name"]],
  Phenotype = s3[["Assigned Phenotype"]],
  "Accession Number" = na0(s3[["Macaulay Library Accession Number"]]),
  Source = na0(s3$Source),
  check.names = FALSE, stringsAsFactors = FALSE)

other <- setdiff(tips, s3[["Tip Name"]])
tok <- strsplit(other, "_")
unknown <- data.frame(
  "Tip Number" = match(other, tips),
  "Tip Name" = other,
  Genus = vapply(tok, `[`, "", 1),
  Species = vapply(tok, `[`, "", 2),
  "Corrected Species Name" = vapply(tok, function(x) paste(x[1], x[2]), ""),
  Phenotype = "??",
  "Accession Number" = "",
  Source = "",
  check.names = FALSE, stringsAsFactors = FALSE)

out <- rbind(known, unknown)
out <- out[order(out[["Tip Number"]]), ]
stopifnot(nrow(out) == length(tips), !anyDuplicated(out[["Tip Name"]]))
write.csv(out, out_path, row.names = FALSE)
message("wrote ", out_path, " (", nrow(known), " tips from Data S3, ", nrow(unknown), " not assessed)")
