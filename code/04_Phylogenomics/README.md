# 04_Phylogenomics: Paramormyrops species tree and ancestral-state reconstruction of P0

Produces: Fig. 2B (Paramormyrops chronogram and stochastic character maps), Fig. S2 (mormyrid reconstruction),
Fig. S3 (species tree with gene and site concordance factors).

## Inputs
- `output_data/03_QC/${invariant_and_snps_vcf}` — filtered biallelic SNPs plus invariant sites (03_VCF_QC_Filter).
- `input_data/00_Reference_Genome/${gff_file}` — Liftoff gene models (see `input_data/00_Reference_Genome/README.md`).
- `input_data/04_Phylogenomics/representitives.txt` — one representative individual per species or
  *P. kingsleyae* population.
- `input_data/04_Phylogenomics/ssAll_iqtree_timtree.tree` — time-calibrated Osteoglossiformes tree of Peterson et
  al. (2022, Syst. Biol. 71:1032), used for the mormyrid reconstruction and for the Paramormyrops calibration ages.
- `input_data/04_Phylogenomics/SupplementalDataS3_MormyridEODs.xlsx` — Supplementary Data S3, downloaded from the paper.
- `input_data/01_Terra/data_model/fish_data_2024a.txt` — sample sheet (`Rscript code/00_common/make_sample_sheets.R`).

## Steps
Run from the repository root after `source config.env`; edit the `#SBATCH` headers for your cluster.
1. `bash code/04_Phylogenomics/00_generate_regions_list.sh input_data/00_Reference_Genome/${gff_file} > input_data/00_Reference_Genome/pkings_gene_regions.txt`
   — one `gene@chr:start-end` line per annotated gene.
2. `bash code/04_Phylogenomics/01_SplitVCF_to_Fasta.sh` — splits the gene list into 200 chunks; each array task
   (`split_vcf_to_genes.sb`) extracts its genes from the VCF and writes one whole-gene alignment per gene for all
   individuals (`create_fasta_alignment.py`, IUPAC codes at heterozygous sites).
3. `bash code/04_Phylogenomics/02_Extract_Rep_Sequences.sh`, then `bash code/04_Phylogenomics/02b_RemoveEmptySeqs.sh`
   — keep the representatives and drop empty alignments.
4. `bash code/04_Phylogenomics/03_Run_IQTree_LocusTrees_together.sh` — IQ-TREE 2 model selection and a maximum
   likelihood tree for each gene.
5. `bash code/04_Phylogenomics/06_Run_ASTRAL.sh` — ASTRAL-IV species tree from the gene trees. `run_astral.sb`
   calls `astral4` from ASTER (https://github.com/chaoszhang/ASTER) at `~/ASTER/bin/astral4`; edit the path for
   your installation.
6. `bash code/04_Phylogenomics/07_Run_concord.sh` — gene and site concordance factors on the species tree.
7. `Rscript code/04_Phylogenomics/00_make_mormyrid_phenotypes.R` — writes
   `input_data/04_Phylogenomics/mormyrid_specimens_with_phenotypes.csv` (tip names, species names and phenotypes)
   from Data S3 and the tree.
8. `bash code/04_Phylogenomics/08_Render_PHYLO_report.sh` — renders
   `phylogeny_and_character_reconstruction_mormyrids.Rmd` (Mk model fits and 1,000 stochastic character maps on
   the mormyrid subtree) and `phylogeny_and_character_reconstruction_paramormyrops.Rmd` (penalized-likelihood
   chronogram with `chronos`, Mk fits and stochastic character maps).

## Outputs
- `output_data/04_Phylogenomics/gene_alignments/`, `extracted_complete_alignments/` — per-gene alignments
- `output_data/04_Phylogenomics/pkings_iqtree_genes.treefile`, `pkings_iqtree_genes_astral`,
  `pkings_iqtree_genes_astral_gcf_scf_tree.cf.*` — gene trees, species tree, concordance factors
- `output_data/04_Phylogenomics/*_mormyrid_phylogeny.report.html`, `*_paramormyrops_phylogeny.report.html` and the
  `.svg` figures written by the two Rmds
