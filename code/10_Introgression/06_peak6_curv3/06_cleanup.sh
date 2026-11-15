## 06_cleanup.sh
## Deletes the intermediate VCFs (_scratch: the per-chromosome GT-only SNP VCFs from 02_make_chr_vcfs.sh) and gzips
## the kept Dinvestigate window files and Dsuite logs (05_summarise.R reads .txt or .txt.gz). Run after 05_summarise.R.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
outdir=${root}/output_data/10_Introgression/peak6_curv3

du -sh ${outdir}/_scratch 2>/dev/null
rm -rf ${outdir}/_scratch
gzip -f ${outdir}/logs/*.log 2>/dev/null
gzip -f ${outdir}/dinvestigate_curv3/*_localFstats_*_50_25.txt ${outdir}/dinvestigate_published_check/*_localFstats_*_50_25.txt 2>/dev/null
du -sh ${outdir}
