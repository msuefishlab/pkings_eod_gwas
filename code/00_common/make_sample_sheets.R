#!/usr/bin/env Rscript
# Build the sample sheets the pipeline reads from Supplementary Data S1 (the authoritative specimen table).
#
# Input (see README, "Sample metadata"):
#   input_data/01_Terra/data_model/SupplementalDataS1_Fish_Specimens.xlsx   Supplementary Data S1 of the paper
# Outputs:
#   input_data/01_Terra/data_model/fish_data_2024a.txt   all 305 genome-sequenced fish
#   input_data/01_Terra/data_model/fish_data_2023.txt    the same fish minus the three outgroup species added later
#                                                         (P. longicaudatus, P. ntotom, P. sp. SN2)
#
# Data S1 labels mixed-waveform fish "Mixed"; the pipeline code calls them "Wobble", so they are relabelled here.
# EOD_FILE holds the Data S2 plot name for each fish; the raw recordings read by 02_Phenotyping/EOD_Analysis.m are
# not part of the public data.
# Usage: Rscript code/00_common/make_sample_sheets.R [path/to/DataS1.xlsx]
suppressPackageStartupMessages(library(readxl))

root <- system("git rev-parse --show-toplevel", intern = TRUE)
dm <- file.path(root, "input_data", "01_Terra", "data_model")
args <- commandArgs(trailingOnly = TRUE)
s1_path <- if (length(args) >= 1) args[1] else file.path(dm, "SupplementalDataS1_Fish_Specimens.xlsx")
if (!file.exists(s1_path)) stop("Supplementary Data S1 not found at ", s1_path, "; see README, 'Sample metadata'")

s1 <- as.data.frame(read_excel(s1_path, sheet = "fish_samples"))
s1 <- s1[s1[["Sample Type"]] == "Genomic", ]
stopifnot(nrow(s1) == 305, !anyDuplicated(s1$PN))

# m/d/yy, no leading zeros, as in the original sheets (read by 11_Drift_vs_Selection/local_ancestry/01_panels.py)
d <- s1$Date
d <- if (inherits(d, c("POSIXct", "Date"))) as.Date(d) else as.Date(as.numeric(d), origin = "1899-12-30")
date_str <- ifelse(is.na(d), "", paste0(as.integer(format(d, "%m")), "/", as.integer(format(d, "%d")), "/",
                                        format(d, "%y")))

eod <- s1$EOD_FILE

na0 <- function(x) ifelse(is.na(x), "", as.character(x))
out <- data.frame(
  "entity:participant_id" = s1$PN,
  Date = date_str,
  EOD_FILE = na0(eod),
  Lat = na0(s1$Lat),
  Long = na0(s1$Long),
  NO = sub("^[A-Z]+-", "", na0(s1$SpecNo)),
  Notes = "",
  PN = s1$PN,
  POP = s1$POP,
  Phenotype = sub("^Mixed$", "Wobble", s1$Phenotype),
  SL = na0(s1$SL),
  Sex = na0(s1$Sex),
  Site = s1$Site,
  SiteName = na0(s1$SiteName),
  SpecNo = na0(s1$SpecNo),
  Species = s1$Species,
  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(names(out)[9] == "POP")   # 07_Phasing/08_split_vcf_by_pop.sb reads column 9

dir.create(dm, recursive = TRUE, showWarnings = FALSE)
write.table(out, file.path(dm, "fish_data_2024a.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

later <- c("Paramormyrops longicaudatus", "Paramormyrops ntotom", "Paramormyrops SN2")
o23 <- out[!out$Species %in% later, ]
# Keep the original column layout (POP in column 9); the workflow-file columns stay empty.
for (col in c("analysis_ready_bam", "analysis_ready_bam_index", "analysis_ready_bam_md5", "output_vcf",
              "output_vcf_index", "set_of_unaligned_bams", "stage1_gvcf", "stage1_gvcf_index", "unaligned_bam_list"))
  o23[[col]] <- ""
o23 <- o23[, c("entity:participant_id", "analysis_ready_bam", "analysis_ready_bam_index", "analysis_ready_bam_md5",
               "output_vcf", "output_vcf_index", "Phenotype", "PN", "POP", "NO", "SpecNo", "Species", "Site",
               "SiteName", "Lat", "Long", "Date", "SL", "Sex", "EOD_FILE", "Notes", "set_of_unaligned_bams",
               "stage1_gvcf", "stage1_gvcf_index", "unaligned_bam_list")]
stopifnot(names(o23)[9] == "POP")
write.table(o23, file.path(dm, "fish_data_2023.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote fish_data_2024a.txt (", nrow(out), " fish) and fish_data_2023.txt (", nrow(o23), " fish) in ", dm)
