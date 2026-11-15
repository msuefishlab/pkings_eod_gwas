#!/usr/bin/env Rscript
# ESM Table S1: the six core GWAS peaks of the main LRT run (PKINGS_ALL_WOB_EXCLUDED).
#
# Peaks:      PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed (core peaks: SNPs with -log10 p > 5 merged within 75 kb).
#             BED start/end are the 1-based positions of the outermost peak SNPs (gemma.peak.bed: min/max ps),
#             not 0-based BED starts; they are reported as is.
# SNPs:       PKINGS_ALL_WOB_EXCLUDED.lmm2.assoc.txt (GEMMA -lmm 2, p_lrt). SNP counts are taken within the core
#             bounds, after the pre-peak-calling filter af > 0.15 (af = ALT frequency) and without it.
# Thresholds: 100 permutation minima (min.txt, p_lrt = V10), same rule as gwas_report.Rmd:
#             95th percentile = significant, 67th = suggestive.
# Genes:      RefSeq GCF_048594095.1-RS_2025_03 (paramormyrops_ncbi_0.4.renamed.gff), all `gene` and
#             `pseudogene` features (every biotype). A gene is listed when it overlaps the core interval or lies
#             within 25 kb of it (GenomicRanges::distance <= 25000), the candidate-gene rule of the ESM RNA-seq
#             Methods. Symbols used in the paper come from input_data/08_Peak_Analysis/pretty_gene_names.csv, then
#             input_data/06_Association/..._GENES_IN_PEAKS_manual_annos_add.txt (hand-curated names), then the NCBI
#             symbol when it is not a LOC ID.
# Run after gwas_report.Rmd (PEAKS.bed) and peak_analysis_new.Rmd (..._GENES_IN_PEAKS.txt).
#
# Rscript code/06_Association/make_table_s1.R
# Writes to output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/:
#   Table_S1_core_peaks.tsv       one row per peak (bounds, index SNP and its context, SNP counts, exon overlap)
#   Table_S1_core_peak_genes.tsv  one row per gene and peak (ID, symbol, description, coordinates, relation)

suppressMessages({
  library(data.table)
  library(GenomicRanges)
})

root <- rprojroot::find_root(".git/index")
run <- "PKINGS_ALL_WOB_EXCLUDED"
run_dir <- file.path(root, "output_data", "06_Association", run)
peaks_file <- file.path(run_dir, paste0(run, "_PEAKS.bed"))
assoc_file <- file.path(run_dir, paste0(run, ".lmm2.assoc.txt"))
min_file <- file.path(run_dir, paste0(run, "_permution"), "min.txt")
old_genes_file <- file.path(run_dir, paste0(run, "_GENES_IN_PEAKS.txt"))
manual_file <- file.path(root, "input_data", "06_Association", paste0(run, "_GENES_IN_PEAKS_manual_annos_add.txt"))
pretty_file <- file.path(root, "input_data", "08_Peak_Analysis", "pretty_gene_names.csv")
gff_file <- file.path(root, "input_data", "00_Reference_Genome", "paramormyrops_ncbi_0.4.renamed.gff")
out_peaks <- file.path(run_dir, "Table_S1_core_peaks.tsv")
out_genes <- file.path(run_dir, "Table_S1_core_peak_genes.tsv")

MAF <- 0.15 # gwas_report.Rmd: filter(af > MAF_thresh), as run
MARGIN <- 25000 # ESM RNA-seq Methods: genes overlapping a peak or within 25 kb

# ---- thresholds: same rule as gwas_report.Rmd ----
p_perm <- sort(fread(min_file, header = FALSE)$V10)
n_perm <- length(p_perm)
significant <- -log10(p_perm[n_perm - floor(n_perm * 0.95)])
suggestive <- -log10(p_perm[n_perm - floor(n_perm * 0.67)])
cat(sprintf("%d permutations: significant -log10 p > %.2f, suggestive > %.2f\n", n_perm, significant, suggestive))

# ---- core peaks ----
peaks <- fread(peaks_file, header = FALSE, col.names = c("chr", "start", "end", "peak", "max_log_p", "index_pos"))
peaks[, chr := paste0("chr", chr)]
setorder(peaks, peak)

# ---- SNPs inside the core bounds ----
assoc <- fread(assoc_file, select = c("chr", "ps", "af", "p_lrt"))
assoc <- assoc[chr %in% peaks$chr]
assoc[, log_p := -log10(p_lrt)]
snps <- assoc[peaks, on = .(chr, ps >= start, ps <= end), nomatch = NULL, .(peak, chr, ps = x.ps, af, log_p)]
snps <- snps[!is.na(log_p)]
snps[, pass_af := af > MAF]
snps[, pass_sym_maf := pmin(af, 1 - af) > MAF]

counts <- snps[, .(
  n_snps_tested = .N,
  n_snps_gt5 = sum(pass_af & log_p > 5),
  n_sig = sum(pass_af & log_p > significant),
  n_sugg = sum(pass_af & log_p > suggestive),
  n_sig_no_af_filter = sum(log_p > significant),
  n_sugg_no_af_filter = sum(log_p > suggestive),
  n_sig_symmetric_maf = sum(pass_sym_maf & log_p > significant),
  n_sugg_symmetric_maf = sum(pass_sym_maf & log_p > suggestive),
  max_log_p_check = max(log_p[pass_af]),
  index_pos_check = ps[pass_af][which.max(log_p[pass_af])],
  max_log_p_no_af_filter = max(log_p),
  index_pos_no_af_filter = ps[which.max(log_p)]
), by = peak]
peaks <- counts[peaks, on = "peak"]
stopifnot(all.equal(peaks$max_log_p, peaks$max_log_p_check), all(peaks$index_pos == peaks$index_pos_check))
# BED bounds are SNP positions: the outermost af-filtered SNPs with -log10 p > 5 sit exactly on them
bounds <- snps[pass_af & log_p > 5, .(lo = min(ps), hi = max(ps)), by = peak][peaks, on = "peak"]
stopifnot(all(bounds$lo == bounds$start), all(bounds$hi == bounds$end))

cat("\nSNPs in core bounds removed by the af filter that pass the suggestive threshold:\n")
print(snps[!pass_af & log_p > suggestive, .(peak, chr, ps, af, log_p = round(log_p, 2))])

# ---- annotation: genes, transcripts, exons and CDS on the peak chromosomes ----
gff <- fread(gff_file, sep = "\t", header = FALSE, quote = "",
  col.names = c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attr"))
gff <- gff[chr %in% peaks$chr & type != "region" & type != "cDNA_match"]
attr_val <- function(attr, key) {
  v <- regmatches(attr, regexpr(paste0("(^|;)", key, "=[^;]*"), attr))
  out <- rep(NA_character_, length(attr))
  out[grepl(paste0("(^|;)", key, "="), attr)] <- sub(paste0("^;?", key, "="), "", v)
  out
}
gff[, `:=`(ID = attr_val(attr, "ID"), Parent = sub(",.*", "", attr_val(attr, "Parent")))]

genes <- gff[type %in% c("gene", "pseudogene")]
genes[, `:=`(
  gene_name = attr_val(attr, "Name"),
  gene_id = sub(".*GeneID:([0-9]+).*", "\\1", attr_val(attr, "Dbxref")),
  biotype = attr_val(attr, "gene_biotype"),
  description = vapply(attr_val(attr, "description"), URLdecode, "", USE.NAMES = FALSE)
)]
# exons/CDS -> transcript -> gene (pseudogene exons can hang off the gene directly)
tx2gene <- gff[!is.na(Parent) & Parent %in% genes$ID, .(ID, gene_ID = Parent)]
tx2gene <- unique(rbind(tx2gene, genes[, .(ID, gene_ID = ID)]))
feat <- gff[type %in% c("exon", "CDS")]
feat <- tx2gene[feat, on = c(ID = "Parent"), nomatch = NULL][, .(chr, type, start, end, gene_ID)]

genes_gr <- GRanges(genes$chr, IRanges(genes$start, genes$end), strand = genes$strand, ID = genes$ID)
exon_gr <- with(feat[type == "exon"], GRanges(chr, IRanges(start, end), gene_ID = gene_ID))
cds_gr <- with(feat[type == "CDS"], GRanges(chr, IRanges(start, end), gene_ID = gene_ID))

# ---- symbols used in the paper ----
pretty <- fread(pretty_file, encoding = "UTF-8")
setnames(pretty, 1, "Name") # header carries a BOM
manual <- fread(manual_file, fill = TRUE, sep = "\t")
manual <- manual[!is.na(manual_name) & manual_name != "", .(Name = gene, manual_name)]
old_ids <- fread(old_genes_file)$ID

genes <- pretty[genes, on = c(Name = "gene_name")]
setnames(genes, c("Name", "display_name"), c("gene_name", "symbol_pretty"))
genes <- manual[genes, on = c(Name = "gene_name")]
setnames(genes, "Name", "gene_name")
genes[, symbol_paper := fcoalesce(symbol_pretty, manual_name,
  fifelse(grepl("^LOC[0-9]+$", gene_name), NA_character_, gene_name))]
genes[, symbol_source := fcase(
  !is.na(symbol_pretty), "pretty_gene_names.csv",
  !is.na(manual_name), "manual_annos_add",
  !is.na(symbol_paper), "NCBI symbol",
  default = "none (LOC ID only)")]
genes[, symbol_conflict := fifelse(!is.na(symbol_pretty) & !is.na(manual_name) &
  tolower(symbol_pretty) != tolower(manual_name),
paste0("manual_annos_add says ", manual_name), "")]
genes[, in_current_table_S1 := ID %in% old_ids]
setkey(genes, ID)

label <- function(id) {
  g <- genes[J(id)]
  fifelse(is.na(g$symbol_paper) | g$symbol_paper == g$gene_name, g$gene_name,
    paste0(g$symbol_paper, " (", g$gene_name, ")"))
}

# position of x relative to gene g: "5'" when x lies on the TSS side of the gene, "3'" on the far side
end_side <- function(x_left_of_gene, strand) {
  fifelse(x_left_of_gene == (strand == "+"), "5' (upstream)", "3' (downstream)")
}

# ---- genes overlapping each core peak or within 25 kb ----
gene_rows <- rbindlist(lapply(seq_len(nrow(peaks)), function(i) {
  p <- peaks[i]
  pg <- GRanges(p$chr, IRanges(p$start, p$end))
  d <- distance(genes_gr, pg)
  keep <- which(!is.na(d) & d <= MARGIN)
  g <- genes[J(genes_gr$ID[keep])]
  g[, dist_bp := d[keep]]
  g[, overlaps_core := dist_bp == 0 & start <= p$end & end >= p$start]
  g[, overlap_bp := pmax(0L, pmin(end, p$end) - pmax(start, p$start) + 1L)]
  g[, overlap_detail := fcase(
    !overlaps_core, "",
    start >= p$start & end <= p$end, "gene inside core",
    start <= p$start & end >= p$end, "core inside gene",
    start < p$start, "gene overlaps core start",
    default = "gene overlaps core end")]
  g[, side := fcase(overlaps_core, "", end < p$start, "left", default = "right")]
  # for a margin gene: where the core peak lies relative to the gene (strand-aware)
  g[, peak_relative_to_gene := fifelse(overlaps_core, "", end_side(side == "right", strand))]
  g[, relation := fifelse(overlaps_core,
    sprintf("overlaps core (%s; %.1f kb overlap)", overlap_detail, overlap_bp / 1000),
    sprintf("25-kb margin only: %.1f kb %s of core; core is %s of gene", dist_bp / 1000, side,
      peak_relative_to_gene))]
  # exons of this gene inside the core interval, and the index SNP
  ex <- exon_gr[exon_gr$gene_ID %in% g$ID]
  ex_in <- unique(ex$gene_ID[overlapsAny(ex, pg)])
  g[, exons_in_core := ID %in% ex_in]
  g[, contains_index_snp := start <= p$index_pos & end >= p$index_pos]
  g[, `:=`(peak = p$peak, margin_only = !overlaps_core)]
  g
}))
setorder(gene_rows, peak, start)

# ---- index SNP context and exonic SNPs ----
snp_context <- function(chr, pos) {
  sg <- GRanges(chr, IRanges(pos, pos))
  in_gene <- genes_gr$ID[overlapsAny(genes_gr, sg)]
  if (length(in_gene)) {
    in_cds <- unique(cds_gr$gene_ID[overlapsAny(cds_gr, sg)])
    in_exon <- unique(exon_gr$gene_ID[overlapsAny(exon_gr, sg)])
    parts <- vapply(in_gene, function(id) {
      bt <- genes[J(id)]$biotype
      kind <- if (id %in% in_cds) "exonic (CDS)" else if (id %in% in_exon) {
        if (bt == "protein_coding") "exonic (UTR)" else paste0("exonic (", bt, ")")
      } else "intronic"
      paste0(kind, ": ", label(id))
    }, "")
    return(paste(parts, collapse = "; "))
  }
  "intergenic"
}

# nearest protein-coding gene on each side of a position, strand-aware
flanks <- function(chr_, pos) {
  pc <- genes[chr == chr_ & biotype == "protein_coding"]
  left <- pc[end < pos][which.max(end)]
  right <- pc[start > pos][which.min(start)]
  fmt <- function(g, snp_left_of_gene) {
    if (!nrow(g)) return(NA_character_)
    dist <- if (snp_left_of_gene) g$start - pos else pos - g$end
    sprintf("%.1f kb %s of %s (%s strand)", dist / 1000, end_side(snp_left_of_gene, g$strand), label(g$ID), g$strand)
  }
  paste(c(fmt(left, FALSE), fmt(right, TRUE)), collapse = "; ")
}

peaks[, index_snp_context := mapply(snp_context, chr, index_pos)]
peaks[, index_snp_flanking_genes := mapply(flanks, chr, index_pos)]

# SNPs above the suggestive threshold (af-filtered) that fall in exons or CDS
sugg <- snps[pass_af & log_p > suggestive]
sugg_gr <- GRanges(sugg$chr, IRanges(sugg$ps, sugg$ps))
sugg[, in_exon := overlapsAny(sugg_gr, exon_gr)]
sugg[, in_cds := overlapsAny(sugg_gr, cds_gr)]
sugg[, context := mapply(snp_context, chr, ps)]
exonic <- sugg[, .(n_sugg_snps_exonic = sum(in_exon), n_sugg_snps_cds = sum(in_cds),
  exonic_sugg_snps = paste(sprintf("%s:%d %s", chr[in_exon], ps[in_exon], context[in_exon]), collapse = " | ")),
by = peak]
peaks <- exonic[peaks, on = "peak"]
# where the suggestive SNPs of each peak fall (e.g. "intergenic 10; intronic: c6orf132 (LOC111844745) 3")
sugg_ctx <- sugg[, .N, by = .(peak, context)][order(peak, -N)][
  , .(suggestive_snp_contexts = paste(sprintf("%s %d", context, N), collapse = "; ")), by = peak]
peaks <- sugg_ctx[peaks, on = "peak"]
# suggestive SNPs inside each listed gene's span
gene_rows[, n_sugg_in_gene := mapply(function(pk, s, e) sugg[peak == pk & ps >= s & ps <= e, .N], peak, start, end)]
genes_exons <- gene_rows[exons_in_core == TRUE, .(genes_with_exons_in_core = paste(label(ID), collapse = ", ")), by = peak]
peaks <- genes_exons[peaks, on = "peak"]
peaks[is.na(genes_with_exons_in_core), genes_with_exons_in_core := "none"]

# ---- write ----
peak_out <- peaks[, .(
  Peak = peak, Chromosome = chr, Start = start, End = end,
  Length_kb = round((end - start + 1) / 1000, 1),
  Index_SNP = paste0(chr, ":", index_pos),
  Max_neg_log10_p = round(max_log_p, 2),
  nSNPs_significant = n_sig, nSNPs_suggestive = n_sugg,
  nSNPs_significant_no_af_filter = n_sig_no_af_filter, nSNPs_suggestive_no_af_filter = n_sugg_no_af_filter,
  nSNPs_gt5 = n_snps_gt5, nSNPs_tested = n_snps_tested,
  Index_SNP_context = index_snp_context, Index_SNP_flanking_genes = index_snp_flanking_genes,
  Suggestive_SNP_contexts = suggestive_snp_contexts,
  Genes_with_exons_in_core = genes_with_exons_in_core,
  nSNPs_suggestive_exonic = n_sugg_snps_exonic, nSNPs_suggestive_CDS = n_sugg_snps_cds,
  Exonic_suggestive_SNPs = exonic_sugg_snps,
  Thresholds = sprintf("significant > %.2f, suggestive > %.2f (-log10 p; %d permutations)", significant, suggestive, n_perm)
)]
gene_out <- gene_rows[, .(
  Peak = peak, Chromosome = chr, Gene_ID = gene_name, NCBI_GeneID = gene_id, Symbol = symbol_paper,
  Symbol_source = symbol_source, Symbol_conflict = symbol_conflict, Biotype = biotype, Description = description,
  Strand = strand, Gene_start = start, Gene_end = end, Relation = relation, Distance_kb = round(dist_bp / 1000, 1),
  Margin_only = margin_only, Exons_in_core = exons_in_core, Contains_index_SNP = contains_index_snp,
  nSNPs_suggestive_in_gene = n_sugg_in_gene,
  In_current_Table_S1 = in_current_table_S1
)]
fwrite(peak_out, out_peaks, sep = "\t")
fwrite(gene_out, out_genes, sep = "\t")
cat("\nWrote", out_peaks, "\nWrote", out_genes, "\n")

cat("\nGenes in the current Table S1 that fail the core-or-25-kb rule:\n")
print(genes[J(setdiff(old_ids, gene_rows$ID)), .(chr, gene_name, symbol_paper, start, end)])
cat("\nBiotypes listed:", paste(sort(unique(gene_rows$biotype)), collapse = ", "), "\n")
cat("Symmetric MAF > 0.15 counts (significant / suggestive):",
  paste(sprintf("P%d %d/%d", peaks$peak, peaks$n_sig_symmetric_maf, peaks$n_sugg_symmetric_maf), collapse = "; "), "\n")
