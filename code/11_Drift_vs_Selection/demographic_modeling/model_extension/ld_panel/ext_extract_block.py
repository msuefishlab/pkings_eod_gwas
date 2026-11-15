#!/usr/bin/env python3
"""ext_extract_block.py -- genotypes of one mask block for the union sample list, joined to the site-table flags.

Genotypes: 23a_observed_ld_decay.genotypes() imported unchanged (bcftools view -r -S -v snps -m2 -M2 | query GT),
i.e. exactly the SNP set and genotype coding simulation/23a used. Column order = the order bcftools writes with -S
(read from the header and stored as `samples`).

Usage: ext_extract_block.py <code_sim_dir> <vcf> <samples.txt> <chunks.tsv> <sitetable_npz_dir> <out_dir> <chrom> <block> <f>
Writes <out_dir>/<chrom>_<block>.npz: pos, G, samples, callable, pol, in_table, an_b, ac_b, an_t, ac_t, f, L
"""
import importlib.util, os, subprocess, sys
import numpy as np

code, vcf, samples, chunks_fn, st_dir, out_dir, ch, b, f = sys.argv[1:10]
b, f = int(b), float(f)
spec = importlib.util.spec_from_file_location("a23", os.path.join(code, "23a_observed_ld_decay.py"))
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
s, e = b * 1_000_000, (b + 1) * 1_000_000
pos, G = A.genotypes(vcf, samples, ch, s, e)
hdr = subprocess.run(["bcftools", "view", "-h", "-S", samples, vcf], capture_output=True, text=True, check=True).stdout
names = hdr.rstrip("\n").split("\n")[-1].split("\t")[9:]
assert G.shape[1] in (0, len(names)), (G.shape, len(names))
# site-table rows of this block
ids = []
for line in open(chunks_fn):
    p = line.rstrip("\n").split("\t")
    if p[0] == "chunk_id" or p[1] != ch:
        continue
    if int(p[2]) <= e and int(p[3]) > s:
        ids.append(p[0])
tp, tc, tpol, tab = [], [], [], {k: [] for k in ("an_b", "ac_b", "an_t", "ac_t")}
for i in ids:
    z = np.load(os.path.join(st_dir, f"st.{i}.npz"))
    k = (z["chrom"] == ch) & (z["pos"] > s) & (z["pos"] <= e)
    tp.append(z["pos"][k]); tc.append(z["callable"][k]); tpol.append(z["pol"][k])
    for kk in tab:
        tab[kk].append(z[kk][k])
tp = np.concatenate(tp) if tp else np.zeros(0, np.int64)
tc = np.concatenate(tc) if tc else np.zeros(0, bool)
tpol = np.concatenate(tpol) if tpol else np.zeros(0, np.int8)
tab = {kk: (np.concatenate(v) if v else np.zeros(0, np.int16)) for kk, v in tab.items()}
assert np.all(np.diff(tp) > 0), "site table positions not unique/sorted"
j = np.searchsorted(tp, pos)
j = np.clip(j, 0, max(len(tp) - 1, 0))
hit = (len(tp) > 0) & (tp[j] == pos) if len(tp) else np.zeros(len(pos), bool)
cal = np.where(hit, tc[j], False) if len(tp) else np.zeros(len(pos), bool)
pol = np.where(hit, tpol[j], 0).astype(np.int8) if len(tp) else np.zeros(len(pos), np.int8)
extra = {kk: (np.where(hit, v[j], -1).astype(np.int16) if len(tp) else np.full(len(pos), -1, np.int16))
         for kk, v in tab.items()}
L = int(round(f * 1_000_000))
np.savez_compressed(os.path.join(out_dir, f"{ch}_{b}.tmp.npz"), pos=pos, G=G, samples=np.array(names),
                    callable=cal, pol=pol, in_table=hit, f=f, L=L, n_table_rows=len(tp),
                    n_table_callable=int(tc.sum()), **extra)
os.replace(os.path.join(out_dir, f"{ch}_{b}.tmp.npz"), os.path.join(out_dir, f"{ch}_{b}.npz"))
print(f"{ch}:{b} snps={len(pos)} in_table={int(hit.sum())} callable={int(cal.sum())} table_rows={len(tp)} "
      f"table_callable={int(tc.sum())}", flush=True)
