# images: container definitions

The pipeline runs its tools from Singularity images in this directory, named by the `*_image` variables in
`config.env`. Build each image with Docker from the definition below, then convert it to a `.sif` file here
(`.sif` files are not tracked).

| `.sif` (config variable) | Definition | Docker Hub image used for the paper | Used by |
|---|---|---|---|
| `gatk_4.3.0.0.sif` (`gatk_image`) | pulled | `broadinstitute/gatk:4.3.0.0` | 01, 03, 07 |
| `gwas_tools.sif` (`gwas_tools_image`) | `Dockerfile` (+ `dnaml`) | `jasongallant/gwas_tools:latest` | 03–08, 11, 13, 14 |
| `rnaseq_tools.sif` (`rnaseq_tools_image`) | `RNAseqDockerfile` (+ `docker_helper.sh`) | `jasongallant/rnaseq_tools:latest` | 09 |
| `liftoff.sif` (`liftoff_image`) | `LIFTOFF_Dockerfile` (Liftoff, RagTag, BLAST, NCBI datasets) | `jasongallant/liftoff:latest` | Liftoff annotations (`input_data/00_Reference_Genome/README.md`) |
| `dsuite_tools.sif` | `DsuiteDocker` | `jasongallant/dsuite_tools:latest` | 10 |
| `dadi_tools.sif` (`dadi_image`) | `DadiDocker` | `jasongallant/dadi_tools:latest` | 11 (demographic modeling) |
| `sweepfinder2.sif` (`sweepfinder2_image`) | `Dockerfile.sweepfinder2` | `jasongallant/sweepfinder2:latest` | 11 (sweep detection) |
| `sweed.sif` (`sweed_image`) | `Dockerfile.sweed` | none (build it locally) | 11 (folded CLR) |
| `rehh.sif` (`rehh_image`) | `Dockerfile.rehh` | `jasongallant/rehh:latest` | 11 (EHH) |
| `smcpp.sif` | pulled | `terhorst/smcpp:latest` | 12 |

`config.env` also defines `rmats_image` and `ggsashimi_image`; no script in this repository uses them.

## Build
Pulled images:

    cd images
    singularity pull docker://broadinstitute/gatk:4.3.0.0              # writes gatk_4.3.0.0.sif
    singularity pull smcpp.sif docker://terhorst/smcpp:latest

Images built from a definition (`<definition>` and `<name>` from the table; `<name>` is the `.sif` name without
the extension):

    cd images
    docker build --platform linux/amd64 -f <definition> -t <name> .
    # then either push to a registry and build the .sif from it on the cluster ...
    docker tag <name> <your_registry_user>/<name> && docker push <your_registry_user>/<name>
    singularity build <name>.sif docker://<your_registry_user>/<name>:latest
    # ... or, on a machine with both Docker and Singularity
    singularity build <name>.sif docker-daemon://<name>:latest

For example, `docker build --platform linux/amd64 -f Dockerfile -t gwas_tools .` and then
`singularity build gwas_tools.sif docker://<your_registry_user>/gwas_tools:latest`. To use the images the paper
was run with instead, build from the Docker Hub names in the table, e.g.
`singularity build gwas_tools.sif docker://jasongallant/gwas_tools:latest`.

## Versions
- `Dockerfile` (gwas_tools): Bioconductor 3.18 base; pixy 1.2.10.beta2, GEMMA at commit 8132837, IQ-TREE 2.3.1,
  seqkit 2.6.1, PLINK 1.9 (build 20231211), ADMIXTURE 1.3.0, qctool 2.2.0, SHAPEIT v2.r904; bcftools, vcftools and
  samtools from the Debian packages of the base image. Other tools (SNPhylo, vcf2phylip, gemma-wrapper, gff3sort,
  cflib) are cloned at their current default branch.
- `DadiDocker`: dadi 2.4.3 with numpy < 2. The image used for the paper holds dadi 2.4.3, numpy 1.26.4,
  scipy 1.13.1, demes 0.2.3 and nlopt 2.7.1.
- `RNAseqDockerfile`: STAR 2.7.11b, Trimmomatic 0.39.
- Not pinned, so a rebuild can differ from the images used for the paper: the `latest` Docker Hub tags above, the
  `rocker/tidyverse:latest` base of `Dockerfile.rehh`, the untagged `continuumio/miniconda3` base of
  `LIFTOFF_Dockerfile`, and the Dsuite and SweeD sources, which `DsuiteDocker` and `Dockerfile.sweed` clone from
  their default branches.
