#!/usr/bin/env python3
"""06_fit_dip.py -- diploid three-state local-ancestry fits (sensitivity analysis).

Same fish, markers (primary panel) and grid as 04_fit_hap.py primary, with the 0/1/2-BAM-copy HMM:
  phasedgeno  genotypes from the phased VCF, phase ignored (tests switch errors);
  unphased    calls from the unphased filtered callset at the marker positions: PASS biallelic SNPs, REF/ALT
              matching the phased VCF, genotype missing if FMT/DP < 8 or GQ < 20 (tests imputation as well).
Per fish: ML (T, m) on the genome-wide markers with chromosome-bootstrap intervals, q_m on the genotypes
(dosage/2), call rate. Writes <out_dir>/dip_<source>/fits.tsv and surfaces.npy; the unphased genotypes go to
<out_dir>/unphased_calls/<chr>.tsv.gz.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/06_fit_dip.py <out_dir> <unphased_vcf> [workers]
"""
import gzip, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import la_lib as L

ROLES = ("PB", "PA", "BN")
PH = "output_data/11_Drift_vs_Selection/sweep_detection/APA.polarized.vcf.gz"


def extract(out, vcf, fish, M, workers):
    d = f"{out}/unphased_calls"; os.makedirs(d, exist_ok=True)
    open(f"{d}/samples.txt", "w").write("\n".join(fish) + "\n")

    def one(c):
        ch = f"chr{c}"
        s = M["chrom"] == c
        reg = f"{d}/{ch}.pos.tsv"
        open(reg, "w").write("".join(f"{ch}\t{p}\n" for p in M["pos"][s]))
        dst = f"{d}/{ch}.tsv.gz"
        if not os.path.exists(dst):
            cmd = (f"bcftools view -R {reg} -S {d}/samples.txt -f PASS -v snps -m2 -M2 -Ou {vcf} 2>/dev/null | "
                   f"bcftools filter -S . -e 'FMT/DP<8 | FMT/GQ<20' -Ou 2>/dev/null | "
                   f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' | gzip > {dst}.tmp && mv {dst}.tmp {dst}")
            subprocess.run(cmd, shell=True, check=True)
        ref = subprocess.run(f"bcftools query -R {reg} -f '%POS\\t%REF\\t%ALT\\n' {PH} 2>/dev/null",
                             shell=True, capture_output=True, text=True, check=True).stdout
        return c, ref

    with ThreadPoolExecutor(workers) as ex:
        refs = dict(ex.map(one, range(1, 26)))
    return refs


def load_unphased(out, M, refs, nf):
    code = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2, "0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2}
    G = np.full((M["pos"].size, nf), -1, np.int8)
    idx = {(int(c), int(p)): i for i, (c, p) in enumerate(zip(M["chrom"], M["pos"]))}
    n_found = n_match = 0
    for c in range(1, 26):
        alle = {}
        for line in refs[c].splitlines():
            p, r, a = line.split("\t"); alle[int(p)] = (r, a)
        with gzip.open(f"{out}/unphased_calls/chr{c}.tsv.gz", "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                p = int(f[1]); key = (c, p)
                if key not in idx:
                    continue
                n_found += 1
                if alle.get(p) != (f[2], f[3]):
                    continue
                n_match += 1
                G[idx[key]] = [code.get(g, -1) for g in f[4:]]
    return G, n_found, n_match


def main():
    out, vcf = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    info = {}
    for line in open(f"{out}/samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sample":
            info[f[0]] = dict(group=f[2], q=f[3], site=f[4], roles=set(f[5].split(",")))
    M = L.load_markers(out, "primary")
    hs = M["hap_sample"]
    fish = [s for s in dict.fromkeys(hs) if s in info and not (info[s]["roles"] & set(ROLES))]
    gw = M["gw"].astype(bool)
    Gph = np.stack([M["H"][:, np.flatnonzero(hs == s)].sum(1) for s in fish], 1).astype(np.int8)
    refs = extract(out, vcf, fish, M, min(workers, 8))
    Gun, n_found, n_match = load_unphased(out, M, refs, len(fish))
    print(f"unphased: {n_found} of {M['pos'].size} marker sites found as PASS biallelic SNPs; {n_match} with matching "
          f"REF/ALT", flush=True)
    called = Gun >= 0
    agree = (Gun == Gph) | ~called
    print(f"unphased vs phased-VCF genotype agreement where called: {(Gun[called] == Gph[called]).mean():.4f}", flush=True)
    for src, G in (("phasedgeno", Gph), ("unphased", Gun)):
        od = f"{out}/dip_{src}"; os.makedirs(od, exist_ok=True)
        S = L.grid_surfaces("dip", L.per_chrom_split(M["chrom"], M["pos"], G, M["pBe"], M["pTe"], gw), workers)
        np.save(f"{od}/surfaces.npy", S.astype(np.float32))
        fits = L.fit_from_surfaces(S)
        rows = []
        for i, s in enumerate(fish):
            f = fits[i]
            x = np.where(G[gw, i] >= 0, G[gw, i] / 2, -1.0)[:, None]
            cr = float((G[gw, i] >= 0).mean())
            conc = float((Gun[gw, i] == Gph[gw, i])[Gun[gw, i] >= 0].mean()) if src == "unphased" else 1.0
            rows.append("\t".join([s, info[s]["group"], info[s]["site"], info[s]["q"],
                                   *(f"{v:.1f}" for r in L.R_BRACKET for v in (f["T"] * 1.4 / r, f["T_lo"] * 1.4 / r,
                                                                                f["T_hi"] * 1.4 / r)),
                                   f"{f['m']:.3f}", f"{f['m_lo']:.3f}", f"{f['m_hi']:.3f}", str(int(f["T_edge"])),
                                   str(int(f["m_edge"])), f"{L.projection_q(x, M['pBe'][gw], M['pTe'][gw]):.4f}",
                                   f"{np.log10(f['T']):.4f}", f"{cr:.4f}", f"{conc:.4f}"]))
        hdr = ["sample", "group", "site", "q26c"] + [f"T{x}_r{r}" for r in L.R_BRACKET for x in ("", "_lo", "_hi")] + \
              ["m", "m_lo", "m_hi", "T_at_grid_edge", "m_at_grid_edge", "q_m", "logT", "call_rate",
               "concordance_with_phased"]
        open(f"{od}/fits.tsv", "w").write("\t".join(hdr) + "\n" + "\n".join(rows) + "\n")
        print(f"wrote {od}/fits.tsv", flush=True)


if __name__ == "__main__":
    main()
