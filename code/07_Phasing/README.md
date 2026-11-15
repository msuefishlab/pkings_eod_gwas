# 07_Phasing: read-backed and statistical phasing of the SNP call set

Produces: the phased SNP VCF (`${phased_vcf}`) described in the ESM (SHAPEIT4, no genetic map or reference panel),
used by the haplotype analyses (Fig. S4, Table S2), the sweep scans, local ancestry, the TFBS scan and the
Wald vs LRT comparison in 06.

## Inputs
- `output_data/03_QC/${snps_only_for_phasing}` — filtered biallelic SNPs with original scaffold names (03_VCF_QC_Filter).
- Analysis-ready BAMs from the 01_Terra workflow, one per sample, at
  `${scratch_dir}/output_data/07_Phasing/<SAMPLE>.MIC4273_HAP2_NEWREF_MSU618.aligned.duplicate_marked.sorted.bam`.
  Reads are in SRA (BioProject PRJNA1218444); rerun 01_Terra to make the BAMs. `00_copy_ar_bams_from_terra.sh` copies
  them from a Terra bucket and reads `input_data/07_Phasing/AR_BAMLIST_2024.txt`, which is not shipped: write one
  bucket URL per line from your own Terra sample table. Skip it if the BAMs are already in place.
- `input_data/00_Reference_Genome/${original_reference}` (+ `.fai`, `.dict`) and `${chr_name_map}`.
- `input_data/01_Terra/chrom_intervals.list` and `remaining_intervals.list` — the 25 chromosome-scale scaffolds and
  the remaining scaffolds (made as in `code/01_Terra/README.md`).
- Containers: `${gwas_tools_image}` (bcftools, samtools, WhatsHap), `${gatk_image}`; SHAPEIT 4.2.2 from the cluster
  module system (`module load SHAPEIT4/4.2.2-foss-2023a`; edit for your cluster).

## Steps
Run from the repository root after `source config.env`. Edit the `#SBATCH` lines of the `.sb` scripts for your cluster.

0. `bash code/07_Phasing/00_copy_ar_bams_from_terra.sh` — optional, copy BAMs from Terra (`copy_ar_bams.sb`).
1. `bash code/07_Phasing/01_create_region_lists.sh` — one region list per chromosome-scale scaffold plus one for all
   other scaffolds (`input_data/07_Phasing/region_lists/`, `pseudo_chrs.txt`).
2. `bash code/07_Phasing/02_create_sample_batches.sh` — sample batches of 25 (`input_data/07_Phasing/sample_names_bNN`).
3. `bash code/07_Phasing/03_whatshap_submit.sh input_data/07_Phasing/sample_names_bNN` — WhatsHap read-backed phasing
   per sample and region (`submit_whatshap.sb`); run once per batch.
4. `bash code/07_Phasing/04_squish_outputs.sh input_data/07_Phasing/sample_names_bNN` — merge each sample's
   unplaced-scaffold outputs into one VCF (`squish_outputs.sb`); once per batch.
5. `bash code/07_Phasing/05_run_shapeit.sh` — merge samples per region and phase with SHAPEIT4 (`--use-PS 0.0001`;
   `run_shapeit.sb`).
6. `bash code/07_Phasing/06_run_merge_phased.sh` — concatenate the phased regions (`merge_phased.sb`).
7. `bash code/07_Phasing/07_rename_vcfs.sh` — rename scaffolds to chr1–chr25 with `${chr_name_map}` (`rename_vcfs.sb`).
8. `bash code/07_Phasing/08_split_vcf_by_pop.sb` — optional: per-population phased VCFs (MAF ≥ 0.05 within each
   population), using the POP column of `input_data/01_Terra/data_model/fish_data_2024a.txt`
   (from `Rscript code/00_common/make_sample_sheets.R`). Needs GNU sed (`gsed`) and GNU parallel.

The region lists, sample batches and `shapeit_regions.txt` that steps 1, 2 and 5 write to `input_data/07_Phasing/`
are regenerated from the reference index and the VCF and are not shipped. Intermediate phased chunks are written under
`${scratch_dir}/output_data/07_Phasing/`.

## Outputs
- `output_data/07_Phasing/${phased_vcf}` (+ `.tbi`) — phased SNPs, chromosomes named chr1–chr25.
- `output_data/07_Phasing/<POP>.vcf.gz` — per-population phased VCFs (step 8).
