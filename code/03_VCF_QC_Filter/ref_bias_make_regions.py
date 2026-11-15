#!/usr/bin/env python3
"""Build the region BED for the reference-bias check (the reference is from a BP fish; do TP reads map worse?).

Three non-overlapping region classes on the 25 chromosomes, written in raw-VCF scaffold names:
  peakN  - the six main GWAS peaks (PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed), whole
  hidiv  - 1000 random 5-kb pixy windows from the top 5% of mean pi across the core
           P. kingsleyae populations, with no annotated gene and no peak overlap
           ("highly variable non-coding")
  genome - 2-kb blocks every 50 kb, minus anything overlapping the two classes above

Usage: python3 ref_bias_make_regions.py <root> <out.bed>
"""
import sys, random
import pandas as pd

root, out = sys.argv[1], sys.argv[2]
ref = f"{root}/input_data/00_Reference_Genome"
CORE = ["APA", "BAM", "BAVA", "BEN", "BIK", "BIR", "DOG", "DOV"]
BLOCK, STEP, NHIDIV = 2000, 50000, 1000
random.seed(1)

fai = pd.read_csv(f"{ref}/jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.renamed.fa.fai",
                  sep="\t", header=None, usecols=[0, 1], names=["chr", "len"])
fai = fai[fai.chr.str.match(r"^chr\d+$")]
cmap = pd.read_csv(f"{ref}/jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.chrnames",  # $chr_name_map
                   sep="\t")
chr2scaf = dict(zip(cmap.chr_name, cmap.V1))

peaks = pd.read_csv(f"{root}/output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed",
                    sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "id"])
peaks = [("chr%d" % r.chr, r.start - 1, r.end, "peak%d" % r.id) for r in peaks.itertuples()]

genes = {}
for l in open(f"{ref}/paramormyrops_ncbi_0.4.renamed.sorted.gff"):
    f = l.split("\t")
    if len(f) > 4 and f[2] == "gene":
        genes.setdefault(f[0], []).append((int(f[3]) - 1, int(f[4])))

def overlaps(c, s, e, ivs):
    return any(s < ie and is_ < e for is_, ie in ivs.get(c, []))

peak_ivs = {}
for c, s, e, _ in peaks:
    peak_ivs.setdefault(c, []).append((s, e))

pi = pd.read_csv(f"{root}/output_data/05_PopGen/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites_pi.txt", sep="\t")
pi = pi[pi["pop"].isin(CORE) & pi.chromosome.isin(fai.chr) & (pi.no_sites >= 2500)]
w = pi.groupby(["chromosome", "window_pos_1", "window_pos_2"], as_index=False).avg_pi.mean()
w = w[w.avg_pi >= w.avg_pi.quantile(0.95)]
cand = [(r.chromosome, r.window_pos_1 - 1, r.window_pos_2) for r in w.itertuples()
        if not overlaps(r.chromosome, r.window_pos_1 - 1, r.window_pos_2, genes)
        and not overlaps(r.chromosome, r.window_pos_1 - 1, r.window_pos_2, peak_ivs)]
print(f"top-5% pi windows: {len(w)}; gene- and peak-free: {len(cand)}", file=sys.stderr)
hidiv = [(c, s, e, "hidiv") for c, s, e in random.sample(cand, min(NHIDIV, len(cand)))]

taken = dict(peak_ivs)
for c, s, e, _ in hidiv:
    taken.setdefault(c, []).append((s, e))
genome = [(r.chr, s, s + BLOCK, "genome") for r in fai.itertuples()
          for s in range(0, r.len - BLOCK, STEP) if not overlaps(r.chr, s, s + BLOCK, taken)]

rows = sorted(peaks + hidiv + genome, key=lambda x: (int(x[0][3:]), x[1]))
with open(out, "w") as o:
    for c, s, e, lab in rows:
        o.write(f"{chr2scaf[c]}\t{s}\t{e}\t{lab}\t{c}\n")
n = pd.Series([r[3] for r in rows]).value_counts()
print(n.to_string(), file=sys.stderr)
