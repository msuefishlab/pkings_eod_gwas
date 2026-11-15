# Reference genome and annotation

The *P. kingsleyae* assembly (PKINGS_0.4) is GenBank JBLIQM000000000 (BioProject PRJNA1218457); its RefSeq assembly
and annotation are GCF_048594095.1 (annotation GCF_048594095.1-RS_2025_03). The genome files are not distributed
here; download them from NCBI and build the files below in this directory. File names are the variables in
`config.env`.

Shipped here:
- `jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.chrnames` (`$chr_name_map`) — tab-separated,
  with a header: original scaffold name (`V1`), chromosome name (`chr_name`: `chr1`–`chr25` for the 25 chromosomes,
  `upsN` for unplaced scaffolds) and chromosome number (`chrum`).
- `refseq_accession_to_chr.csv` — RefSeq accession to chromosome name (`NC_132797.1,chr1`), no header; the form
  read by `code/03_VCF_QC_Filter/00_rename_gff.py`.

## Build
Run from this directory after `source ../../config.env`, with the NCBI files
`GCF_048594095.1_PKINGS_0.4_genomic.fna` and `GCF_048594095.1_PKINGS_0.4_genomic.gff` (unzipped) in it.

1. Reference with chromosome names (`$reference`), used by most modules:
   ```bash
   awk -F, 'NR == FNR {m[$1] = $2; next}
        /^>/ {split(substr($0, 2), h, " "); keep = (h[1] in m); if (keep) print ">" m[h[1]]; else print "skipped " h[1] > "/dev/stderr"; next}
        keep' refseq_accession_to_chr.csv GCF_048594095.1_PKINGS_0.4_genomic.fna > ${reference}
   ```
2. The same sequences under the original scaffold names (`$original_reference`), used for read alignment and joint
   genotyping (01_Terra), the raw-VCF steps of 03_VCF_QC_Filter and read-backed phasing (07_Phasing):
   ```bash
   awk -F'\t' 'NR == FNR {if (FNR > 1) m[$2] = $1; next} /^>/ {print ">" m[substr($1, 2)]; next} {print}' \
       ${chr_name_map} ${reference} > ${original_reference}
   ```
3. Index both: `samtools faidx`, GATK `CreateSequenceDictionary`, and for alignment `bwa index ${original_reference}`.
4. RefSeq annotation with chromosome names, and its position-sorted copy:
   ```bash
   python3 ../../code/03_VCF_QC_Filter/00_rename_gff.py GCF_048594095.1_PKINGS_0.4_genomic.gff \
       refseq_accession_to_chr.csv paramormyrops_ncbi_0.4.renamed.gff
   LC_ALL=C sort -s -t $'\t' -k1,1 -k4,4n paramormyrops_ncbi_0.4.renamed.gff > paramormyrops_ncbi_0.4.renamed.sorted.gff
   ```
   (`00_rename_gff.py` needs Python 3 with pandas.) From the RS_2025_03 GFF, these two commands reproduce the files
   used for the paper byte for byte.

The original assembly was renamed the other way round, from scaffold to chromosome names, with
`code/03_VCF_QC_Filter/00_rename_reference.sh` (seqkit, using columns 1–2 of `$chr_name_map`); the filtered VCFs are
renamed with the same map in `code/03_VCF_QC_Filter/08_rename_vcfs.sh`.

All 407 RefSeq sequences have the lengths of the corresponding sequences in the assembly used for the analyses.
The original assembly had one more sequence, `Scaffold_269__1_contigs__length_38835` (`ups269`, 38.8 kb), which
carries tandem copies of the mitochondrial genome and is not part of GCF_048594095.1; a rerun from the public
assembly has no calls on it.

## Liftoff annotations
Two further annotations are gene models of the earlier *P. kingsleyae* assembly (GCF_002872115.1, PKINGS_0.1)
transferred to this assembly with Liftoff (container: `images/LIFTOFF_Dockerfile`), then renamed to chromosome names
with `00_rename_gff.py` / `00_rename_gtf.py` (map: columns 1–2 of `$chr_name_map` as `old,new` CSV):
- `$gff_file` (`...fasta.renamed.sorted.genes_only.gff`, `gene` features only) — NCBI Release 100 gene models;
  `code/04_Phylogenomics/00_generate_regions_list.sh` turns it into `pkings_gene_regions.txt` for the gene
  alignments.
- `$annotation_gff` (`...ensembl.liftoff.renamed.gtf`) — Ensembl gene models, used for the STAR genome index in
  09_RNASeq.
