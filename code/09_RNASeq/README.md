# 09_RNASeq: electric organ RNA-seq expression of genes at the GWAS peaks

Produces: Fig. 4A (`rnaseq_analysis_RSEM_unified_genes.Rmd`) and the RNA-seq section of the ESM Methods.

## Inputs
- Reads: 39 electric organ libraries, one SRA run per library: 31 new (BioProject PRJNA1218444), 3 from
  PRJNA417677 and 5 from PRJNA573805. Run accessions are in Supplementary Data S1 (rows with Sample Type "EO RNA",
  matched to `input_data/09_RNASeq/rnaseq_metadata.txt` by PN).
- `input_data/09_RNASeq/rnaseq_metadata.txt` — one row per sequencing file: SpecNo, Population, R1, R2, Phenotype,
  EOD_FILE, PN, NewPhenotype (shipped). The scripts loop over its SpecNo values; the notebook uses one row per PN.
- `input_data/09_RNASeq/paramormyrops_ncbi_0.4_tx2gene.csv` — transcript-to-gene table (shipped), written by
  `transcript_to_gene_mapping.R` from the NCBI RefSeq annotation GCF_048594095.1-RS_2025_03 (`genomic.gff` from
  NCBI Datasets, saved as `input_data/09_RNASeq/paramormyrops_ncbi_0.4.gff`; not shipped).
- `input_data/09_RNASeq/paramormyrops_ncbi_0.4.renamed.gtf` (not shipped) — the same RefSeq annotation in GTF format,
  with the RefSeq chromosome accessions (NC_132797.1, ...) renamed to chr1–chr25 as for
  `input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.gff` (top-level README, "Reference genome").
- `input_data/00_Reference_Genome/${reference}` — reference genome with renamed chromosomes.
- `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed` (core peaks) and
  `..._GENES_IN_PEAKS.txt` (06_Association, 08_Peak_Analysis); `input_data/08_Peak_Analysis/pretty_gene_names.csv`.
- Container: `${rnaseq_tools_image}` (Trimmomatic 0.39, STAR, RSEM).

## Steps
Run from the repository root after `source config.env`. `META=input_data/09_RNASeq/rnaseq_metadata.txt`. Edit the
`#SBATCH` lines of the `.sb` scripts for your cluster.

1. Reads. From SRA, write each library's pair directly to
   `output_data/09_RNASeq/raw_reads/<SpecNo>_R1_combined.fastq.gz` and `_R2_combined.fastq.gz` (SpecNo as in `META`,
   e.g. `160` for APA_160), for example
   `fasterq-dump --split-files <SRR> -O tmp && gzip -c tmp/<SRR>_1.fastq > output_data/09_RNASeq/raw_reads/160_R1_combined.fastq.gz`
   (and `_2` for R2). With the original per-lane files in `input_data/09_RNASeq/reads/` (names in columns R1 and R2),
   `bash code/09_RNASeq/02_cat_by_specno.sh $META` concatenates them per specimen instead.
2. `bash code/09_RNASeq/03_run_trim.sh $META` — Trimmomatic (`ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5
   TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25`; `run_trim.sb`).
3. `bash code/09_RNASeq/04_run_star_index.sh` — `rsem-prepare-reference --star` on the renamed RefSeq GTF and the
   reference (`run_star_index.sb`).
4. `bash code/09_RNASeq/05_run_star_align.sh $META` — STAR alignment to the transcriptome (`--quantMode
   TranscriptomeSAM`; `run_star_rnaseq.sb`).
5. `bash code/09_RNASeq/06_run_rsem_quantify.sh $META` — `rsem-calculate-expression` (`run_rsem.sb`).
6. Render `code/09_RNASeq/rnaseq_analysis_RSEM_unified_genes.Rmd` (e.g.
   `Rscript -e 'rmarkdown::render("code/09_RNASeq/rnaseq_analysis_RSEM_unified_genes.Rmd")'`): tximport to gene
   level, DESeq2 VST, genes overlapping or within 25 kb of a core peak, heatmaps and boxplots (Fig. 4A).

`Rscript code/09_RNASeq/transcript_to_gene_mapping.R` rebuilds the shipped tx2gene table from the NCBI GFF; it gives
a byte-identical file.

## Outputs
- `output_data/09_RNASeq/trimmed_reads/`, `aligned_reads_transcriptome/`, `pk_ncbi_0.4/` (RSEM/STAR index).
- `output_data/09_RNASeq/quantification/<SpecNo>.isoforms.results` — RSEM estimates.
- `output_data/09_RNASeq/RNAseq_peak_genes.pdf` and the rendered notebook — Fig. 4A.
