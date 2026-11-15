# 01_Terra: read alignment, per-sample GVCFs and joint genotyping on Terra (GATK)

Produces: the raw all-sites VCF `$vcf_file` read by 03_VCF_QC_Filter (Methods, Variant Calling).

## Inputs
- Whole-genome reads: NCBI SRA BioProject PRJNA1218444, one run per fish; run accessions are in the `SRA` column of
  Supplementary Data S1 (rows with `Sample Type` = Genomic).
- `input_data/01_Terra/data_model/fish_data_2024a.txt` — the 305 genome-sequenced fish, generated from Supplementary
  Data S1 by `Rscript code/00_common/make_sample_sheets.R`. Its first column is `entity:participant_id`, so it loads
  into Terra as the participant table without changes.
- `input_data/00_Reference_Genome/${original_reference}` with `.fai`, `.dict` and bwa index files — see
  `input_data/00_Reference_Genome/README.md`. Terra works on the original scaffold names; chromosomes are renamed
  after filtering (03_VCF_QC_Filter).
- GATK interval lists, regenerated from the reference (step 2). `chrom_intervals.list` and
  `remaining_intervals.list` (the 25 chromosome-scale scaffolds and the rest) are shipped; they list the scaffolds of
  the assembly used for the paper, including `Scaffold_269`, which the public assembly lacks.

## Steps
Set `bucket` in `config.env` to your Terra workspace bucket, then run from the repository root after
`source config.env`. The WDLs and their input/output configurations are in `wdl/`; load each into the workspace
and run it with the configuration named below.
1. Reads. Download each run (for example `fasterq-dump --split-files <SRR>`, then gzip), write the local FASTQ paths
   one per line to a file, and run `bash code/01_Terra/01_run_upload.sh <file_list>` — a SLURM array
   (`submit_copy.sb`) that copies them to `gs://${bucket}/reads/`.
2. Interval lists (GATK SplitIntervals, run in `$gatk_image`):
   ```bash
   cd input_data/01_Terra
   fai=../00_Reference_Genome/${original_reference}.fai
   cut -f1 $fai > intervals.list
   cut -f1,2 $fai | sort -nrk 2 | head -25 | cut -f1 > chrom_intervals.list
   cut -f1,2 $fai | sort -nrk 2 | tail -n +26 | cut -f1 > remaining_intervals.list
   gatk="singularity exec --bind $root:/project_root ${gatk_image} /gatk/gatk"
   R=/project_root/input_data/00_Reference_Genome/${original_reference}
   O=/project_root/input_data/01_Terra
   $gatk SplitIntervals -R $R -L $O/intervals.list -scatter 50 -O $O/interval-files \
       -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
   $gatk SplitIntervals -R $R -L $O/chrom_intervals.list --scatter-count 100 -O $O/chrom-interval-files \
       --subdivision-mode INTERVAL_SUBDIVISION
   $gatk SplitIntervals -R $R -L $O/remaining_intervals.list --scatter-count 100 -O $O/remaining-interval-files \
       --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
   ls -1 interval-files | awk -v b=$bucket '{print "gs://" b "/interval-files/" $0}' > scattered-intervals-file.txt
   { ls -1 chrom-interval-files | awk -v b=$bucket '{print "gs://" b "/chrom-interval-files/" $0}'
     ls -1 remaining-interval-files | awk -v b=$bucket '{print "gs://" b "/remaining-interval-files/" $0}'; } > unpadded_intervals_file.txt
   gsutil -m cp -n scattered-intervals-file.txt unpadded_intervals_file.txt gs://${bucket}/
   gsutil -m cp -r -n interval-files chrom-interval-files remaining-interval-files gs://${bucket}/
   cd $root
   ```
   `interval-files/` scatters HaplotypeCaller; `chrom-interval-files/` and `remaining-interval-files/` scatter joint
   genotyping and are also read by 03_VCF_QC_Filter.
3. Reference. Copy `${original_reference}` and its `.fai`, `.dict`, `.amb`, `.ann`, `.bwt`, `.pac` and `.sa` files
   to `gs://${bucket}/reference/`.
4. Workspace attributes: `ref_fasta`, `ref_fasta_index`, `ref_dict`, `ref_amb`, `ref_ann`, `ref_bwt`, `ref_pac`,
   `ref_sa` (the files of step 3); `scattered_calling_intervals_list` = `gs://${bucket}/scattered-intervals-file.txt`;
   `unpadded_intervals_list` = `gs://${bucket}/unpadded_intervals_file.txt`; `gatk_docker` =
   `broadinstitute/gatk:4.0.1.0` (the HaplotypeCaller image used for the paper).
5. Data tables. Upload `fish_data_2024a.txt` as the participant table. Make a sample table with one row per read
   pair and the columns `entity:sample_id`, `participant` and `SN` (both the fish's `PN`), `R1`, `R2`
   (`gs://${bucket}/reads/...`), `LB`, `Lane`, `barcode` and `PL` (`Illumina`); these fill the read group of each
   unaligned BAM. The paper's libraries were split by flowcell lane; for SRA runs use the run accession as `LB` and
   `1` as `Lane` unless you split the reads by lane yourself.
6. `wdl/01_FastqtoSam.1.wdl` with `01_FastqToSam.1.json`, on the sample table — FASTQ to unaligned BAM (Picard
   FastqToSam).
7. For each participant, write the `gs://` paths of its unaligned BAMs to `unaligned_bams_<PN>.txt`, upload these
   files to the bucket and set the participant attribute `unaligned_bam_txt` to them; then run
   `wdl/02_collect_unaligned_bams_by_participant.1.wdl` with `02_collect_unaligned_bams_by_participant.json`.
8. `wdl/03_preprocess_data.2.wdl` with `03_preprocess_data_inputs.json` and `03_preprocess_data_outputs.json`, on
   participants — bwa mem alignment, merge with the unaligned BAM, sort and Picard MarkDuplicates.
9. `wdl/04_haplotype_caller.2.wdl` with `04_haplotype_caller.2.json`, on participants — per-sample GVCFs; map the
   outputs to the participant attributes `gvcf` and `gvcf_index`. About a day for 300 samples.
10. `wdl/05_msu_efishlab_joint-discovery-gatk4.10.wdl` with `05_msu_efishlab_joint-discovery-gatk4.10_inputs.json`
    and `..._outputs.json`, on a participant set of all fish — GenomicsDBImport and GenotypeGVCFs with all sites
    kept, scattered over the unpadded intervals. 12–24 h for 300 samples. `wdl/05b_index_vcf.wdl` indexes the
    final VCF if needed.
11. Download the final VCF and its index to `output_data/01_Terra/${vcf_file}` (and `.tbi`).

## Outputs
- `output_data/01_Terra/${vcf_file}` (+ `.tbi`) — raw joint-genotyped VCF with variant and invariant sites
