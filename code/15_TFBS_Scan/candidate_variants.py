#!/usr/bin/env python3
"""
Candidate regulatory variants at the three TFBS peaks: SNPs and indels in LD with the index SNP, scanned on
real phased local haplotypes for predicted motif changes. This is variant prioritisation,
not an enrichment test: no p-values are computed across variants.

For each locus (sptbn4b, cntnap5, itih3):
  1. Candidates are biallelic SNPs and indels in the union of the TFBS window and the extended GWAS peak with
     r^2 >= --r2 to the index SNP and a genotype call rate >= --min_call. r^2 is the squared Pearson correlation of
     ALT dosages in BP + TP fish (wobble and other phenotypes excluded), pairwise-complete.
  2. The BP allele is the allele more frequent in BP fish.
  3. Local haplotypes span the variant +-F (F = 32; the widest JASPAR motif is 33 bp). Neighbouring SNPs come from
     the phased SNP VCF: the most common local haplotype among BP chromosomes that carry the BP allele (and likewise
     for TP). Indels and SNPs absent from the phased VCF cannot be conditioned, so they use the most common local
     haplotype in each group (hap_conditioned = False). Other indels in the context stay as reference (counted).
  4. Both haplotypes are scanned with FIMO (--text --no-pgc --thresh 1e-2, 0-order genome background). For every
     motif, p_bp and p_tp are the best p of hits overlapping the variant footprint (indels: allele + 1 bp, to
     include the junction); no hit means p = 1e-2 (the scan floor).
  5. A strong change is min(p_bp, p_tp) < 1e-4 and |log10 p_bp - log10 p_tp| >= 1 (>= 10-fold), as in motifbreakR's
     strict mode.

Outputs (output_data/15_TFBS_Scan/candidates/): candidates_{locus}.tsv (one row per variant), motif_changes_{locus}.tsv
(one row per strong variant x motif change), data_table_{locus}.tsv (Supplementary Data S4-S6), candidate_haps.fa,
summary.txt.

Usage:
  source config.env
  python3 code/15_TFBS_Scan/candidate_variants.py
"""

import argparse
import collections
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

F = 32
P_FLOOR = 1e-2
P_PRESENT = 1e-4
MIN_FOLD_LOG10 = 1.0
COMP = str.maketrans("ACGTN", "TGCAN")

# TFBS windows (1-based inclusive) as in README.md, extended peaks from input_data/10_Introgression/peaks.bed,
# and index SNPs (the peak maxima in peaks.bed).
LOCI = {
    "sptbn4b": dict(chr="chr17", win=(11642887, 11742887), peak=(11494062, 11781225), index=11692887, peak_no=5),
    "cntnap5": dict(chr="chr16", win=(10840000, 11040000), peak=(10728257, 11082159), index=10916883, peak_no=4),
    "itih3": dict(chr="chr8", win=(15736956, 15836956), peak=(15771736, 15805670), index=15786956, peak_no=2),
}

# Named groups (README.md, statistic 4b)
NAMED_GROUPS = {
    "BCL6": {"MA0463.3", "MA0731.1"},
    "CEBP": {"MA0102.5", "MA0466.4", "MA0836.3", "MA0837.3", "MA0838.1", "MA1636.2", "MA0019.2"},
}


def env(name, fallback=None):
    return os.environ.get(name, fallback)


def parse_args():
    root = Path(subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True).stdout.strip())
    out_base = root / "output_data/15_TFBS_Scan"
    p = argparse.ArgumentParser(description="Candidate regulatory variants at TFBS peaks (prioritisation).")
    p.add_argument("--vcf", default=str(root / "output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites."
                                         "FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz"),
                   help="Filtered biallelic SNP + indel VCF from 03_VCF_QC_Filter (genotypes for LD and allele frequencies)")
    p.add_argument("--phased", default=env("phased_VCF"), help="Phased SNP VCF (default: $phased_VCF)")
    p.add_argument("--assoc", default=str(root / "output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/"
                                           "PKINGS_ALL_WOB_EXCLUDED.lmm2.assoc.txt"), help="GEMMA LRT output")
    p.add_argument("--fish_data", default=str(root / "input_data/01_Terra/data_model/fish_data_2023.txt"))
    p.add_argument("--ref", default=env("TFBS_REF"), help="Reference FASTA (default: $TFBS_REF)")
    p.add_argument("--gff", default=env("TFBS_GFF"), help="Annotation GFF (default: $TFBS_GFF)")
    p.add_argument("--motifs", default=env("JASPAR_MOTIFS"), help="JASPAR MEME file (default: $JASPAR_MOTIFS)")
    p.add_argument("--families", default=str(root / "input_data/15_TFBS_Scan/jaspar2024_families.tsv"))
    p.add_argument("--bfile", default=str(out_base / "genome.bg"))
    p.add_argument("--fimo", default="fimo", help="FIMO executable (MEME Suite 5.5.8; default: fimo on the PATH)")
    p.add_argument("--r2", type=float, default=0.8)
    p.add_argument("--min_call", type=float, default=0.8, help="Minimum genotype call rate in BP+TP fish")
    p.add_argument("--out_dir", default=str(out_base / "candidates"))
    a = p.parse_args()
    for k in ["phased", "ref", "gff", "motifs"]:
        if not getattr(a, k):
            p.error(f"--{k} not set (source config.env)")
    return a


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"[ERROR] {' '.join(cmd[:4])}...: {r.stderr[-500:]}")
    return r.stdout


def fetch(ref, chrom, start, end):
    out = run(["samtools", "faidx", ref, f"{chrom}:{start}-{end}"])
    return "".join(out.split("\n")[1:]).upper()


def bp_tp_samples(fish_data, vcf):
    pheno = {}
    with open(fish_data) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            pheno[f[0]] = f[6]
    samples = run(["bcftools", "query", "-l", vcf]).split()
    keep = [s for s in samples if pheno.get(s) in ("BP", "TP")]
    return keep, np.array([pheno[s] for s in keep])


def query_genotypes(vcf, region, samples):
    """Rows of (pos, ref, alt, [GT strings]) for biallelic sites."""
    out = run(["bcftools", "query", "-r", region, "-s", ",".join(samples),
               "-f", "%POS\t%REF\t%ALT[\t%GT]\n", vcf], env={**os.environ, "HTS_LOG": "0"})
    rows = []
    for line in out.splitlines():
        f = line.split("\t")
        if "," in f[2] or "*" in f[2]:
            continue
        rows.append((int(f[0]), f[1].upper(), f[2].upper(), f[3:]))
    return rows


def dosage(gts):
    d = np.full(len(gts), np.nan)
    for i, g in enumerate(gts):
        if "." in g:
            continue
        a = g.replace("|", "/").split("/")
        d[i] = sum(x != "0" for x in a)
    return d


def r2_to(index_d, D):
    """Pairwise-complete squared Pearson correlation of each row of D with index_d."""
    ok = ~np.isnan(D) & ~np.isnan(index_d)[None, :]
    x = np.where(ok, D, 0.0)
    y = np.where(ok, index_d[None, :], 0.0)
    n = ok.sum(1)
    sx, sy = x.sum(1), y.sum(1)
    sxx, syy, sxy = (x * x).sum(1), (y * y).sum(1), (x * y).sum(1)
    cov = sxy - sx * sy / n
    vx, vy = sxx - sx ** 2 / n, syy - sy ** 2 / n
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((vx > 0) & (vy > 0), cov ** 2 / (vx * vy), np.nan)


def gemma_logp(assoc, loci):
    """-log10 p_lrt keyed by (chr, pos, ref, alt) for SNPs inside each locus search range."""
    conds = " || ".join(f'($1=="{d["chr"]}" && $3>={min(d["win"][0], d["peak"][0])} && $3<={max(d["win"][1], d["peak"][1])})'
                        for d in loci.values())
    out = run(["awk", "-F", "\t", f"NR>1 && ({conds}) {{print $1\"\\t\"$3\"\\t\"$5\"\\t\"$6\"\\t\"$10}}", assoc])
    res = {}
    for line in out.splitlines():
        c, pos, a1, a0, p = line.split("\t")
        res[(c, int(pos), a1, a0)] = -math.log10(float(p)) if float(p) > 0 else np.nan
    return res


def read_gff(gff, chroms):
    tss, cds, exon, tx = (collections.defaultdict(list) for _ in range(4))
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in chroms:
                continue
            s, e = int(f[3]), int(f[4])
            if f[2] == "CDS":
                cds[f[0]].append((s, e))
            elif f[2] == "exon":
                exon[f[0]].append((s, e))
            elif f[2] in ("mRNA", "transcript", "lnc_RNA"):
                attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                gene = attrs.get("gene", attrs.get("Name", "."))
                tx[f[0]].append((s, e))
                tss[f[0]].append((s if f[6] == "+" else e, gene, f[6]))
    return tss, cds, exon, tx


def annotate(chrom, pos, end, gffd):
    tss, cds, exon, tx = gffd

    def hit(iv):
        return any(s <= end and pos <= e for s, e in iv[chrom])
    cls = "CDS" if hit(cds) else "UTR/ncRNA exon" if hit(exon) else "intron" if hit(tx) else "intergenic"
    best = min(tss[chrom], key=lambda t: abs(t[0] - pos)) if tss[chrom] else (pos, ".", "+")
    # signed distance: negative = upstream of the TSS on the transcript's strand
    d = (pos - best[0]) if best[2] == "+" else (best[0] - pos)
    return cls, best[1], d


def local_haps(ph, lo, hi, cand_pos, cand_ref, cand_alt, groups, bp_code, is_snp):
    """Most common local phased haplotype per group (dict group -> (alleles by pos, freq, conditioned, n))."""
    pos = ph["pos"]
    sel = np.where((pos >= lo) & (pos <= hi))[0]
    cj = None
    if is_snp:
        m = np.where((pos == cand_pos) & (ph["ref"] == cand_ref) & (ph["alt"] == cand_alt))[0]
        cj = int(m[0]) if len(m) else None
    out = {}
    for grp in ("BP", "TP"):
        cols = ph["cols"][grp]
        H = ph["H"][np.ix_(sel, cols)] if len(sel) else np.zeros((0, len(cols)), int)
        ok = (H >= 0).all(0)
        conditioned = cj is not None
        if conditioned:
            want = bp_code if grp == "BP" else 1 - bp_code
            ok &= ph["H"][cj, cols] == want
        n = int(ok.sum())
        if n == 0:
            out[grp] = ({}, 0.0, conditioned, 0)
            continue
        cnt = collections.Counter(map(tuple, H[:, ok].T))
        hap, c = cnt.most_common(1)[0]
        alleles = {int(pos[j]): (ph["ref"][j], ph["alt"][j], v) for j, v in zip(sel, hap) if pos[j] != cand_pos}
        out[grp] = (alleles, c / n, conditioned, n)
    return out


def build_seq(refseq, seg_lo, neighbours, cand_pos, cand_ref, allele):
    s = list(refseq)
    for p, (r, a, v) in neighbours.items():
        i = p - seg_lo
        if s[i] != r:
            raise AssertionError(f"neighbour REF mismatch at {p}: {s[i]} vs {r}")
        if v == 1:
            s[i] = a
    i = cand_pos - seg_lo
    if "".join(s[i:i + len(cand_ref)]) != cand_ref:
        raise AssertionError(f"candidate REF mismatch at {cand_pos}")
    s[i:i + len(cand_ref)] = list(allele)
    return "".join(s), i + 1  # sequence, 1-based footprint start


def main():
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fam = pd.read_csv(args.families, sep="\t").set_index("motif_id")
    group_of = {m: g for g, ms in NAMED_GROUPS.items() for m in ms}
    samples, pheno = bp_tp_samples(args.fish_data, args.vcf)
    ph_samples, ph_pheno = bp_tp_samples(args.fish_data, args.phased)
    print(f"[INFO] BP+TP fish: {len(samples)} in the SNP+indel VCF ({(pheno == 'BP').sum()} BP, "
          f"{(pheno == 'TP').sum()} TP); {len(ph_samples)} in the phased VCF")
    logp = gemma_logp(args.assoc, LOCI)
    gffd = read_gff(args.gff, {d["chr"] for d in LOCI.values()})

    cands, seqs = [], []
    for locus, d in LOCI.items():
        chrom = d["chr"]
        lo, hi = min(d["win"][0], d["peak"][0]), max(d["win"][1], d["peak"][1])
        rows = query_genotypes(args.vcf, f"{chrom}:{lo}-{hi}", samples)
        D = np.vstack([dosage(r[3]) for r in rows])
        ix = [i for i, r in enumerate(rows) if r[0] == d["index"] and len(r[1]) == 1 and len(r[2]) == 1]
        if len(ix) != 1:
            sys.exit(f"[ERROR] index SNP {chrom}:{d['index']} not found once in {args.vcf}")
        r2 = r2_to(D[ix[0]], D)
        call_rate = (~np.isnan(D)).mean(1)
        # a call-rate floor stops r^2 being inflated on the few fish genotyped at poorly called (mostly indel) sites
        keep = np.where((r2 >= args.r2) & (call_rate >= args.min_call))[0]
        indel_pos = np.array([r[0] for r in rows if len(r[1]) != len(r[2])])

        prow = query_genotypes(args.phased, f"{chrom}:{lo - F - 1}-{hi + F + 1}", ph_samples)
        H = np.array([[(-1 if "." in g else int(x)) for g in r[3] for x in g.replace("/", "|").split("|")] for r in prow])
        hcols = {g: np.concatenate([[2 * i, 2 * i + 1] for i in np.where(ph_pheno == g)[0]]) for g in ("BP", "TP")}
        ph = dict(pos=np.array([r[0] for r in prow]), ref=np.array([r[1] for r in prow]),
                  alt=np.array([r[2] for r in prow]), H=H, cols=hcols)
        refseq = fetch(args.ref, chrom, lo - F - 5, hi + F + 50)
        ref_lo = lo - F - 5

        for i in keep:
            pos, ref, alt, _ = rows[i]
            dv = D[i]
            af_bp = np.nanmean(dv[pheno == "BP"]) / 2
            af_tp = np.nanmean(dv[pheno == "TP"]) / 2
            bp_code = 1 if af_bp > af_tp else 0
            is_snp = len(ref) == 1 and len(alt) == 1
            end = pos + len(ref) - 1
            seg_lo, seg_hi = pos - F, end + F
            haps = local_haps(ph, seg_lo, seg_hi, pos, ref, alt, None, bp_code, is_snp)
            seg = refseq[seg_lo - ref_lo: seg_hi - ref_lo + 1]
            bp_allele, tp_allele = (alt, ref) if bp_code else (ref, alt)
            cid = len(cands)
            hs = {}
            for grp, allele in (("BP", bp_allele), ("TP", tp_allele)):
                s, off = build_seq(seg, seg_lo, haps[grp][0], pos, ref, allele)
                flen = len(allele) + (0 if is_snp else 1)
                hs[grp] = (s, off, flen)
                seqs.append((f"C{cid:05d}_{grp}", s))
            ctx_bp = {p: v[2] for p, v in haps["BP"][0].items()}
            ctx_tp = {p: v[2] for p, v in haps["TP"][0].items()}
            n_ctx = sum(ctx_bp.get(p, 0) != ctx_tp.get(p, 0) for p in set(ctx_bp) | set(ctx_tp))
            ann, gene, dist = annotate(chrom, pos, end, gffd)
            cands.append(dict(
                cid=cid, locus=locus, peak=d["peak_no"], var_id=f"{chrom}_{pos}_{ref}_{alt}", chr=chrom, pos=pos,
                type="SNP" if is_snp else ("insertion" if len(alt) > len(ref) else "deletion"), ref=ref, alt=alt,
                bp_allele=bp_allele, tp_allele=tp_allele, r2_index=round(float(r2[i]), 4),
                call_rate=round(float(call_rate[i]), 3),
                is_index=(i == ix[0]), af_alt_bp=round(float(af_bp), 4), af_alt_tp=round(float(af_tp), 4),
                delta_af=round(abs(float(af_bp - af_tp)), 4),
                gemma_log_p=round(logp.get((chrom, pos, ref, alt), np.nan), 3) if is_snp else np.nan,
                in_tfbs_window=d["win"][0] <= pos <= d["win"][1], annotation=ann, nearest_gene=gene,
                dist_to_tss=dist, hap_conditioned=haps["BP"][2], bp_hap_freq=round(haps["BP"][1], 3),
                tp_hap_freq=round(haps["TP"][1], 3), n_context_snp_diffs=n_ctx,
                n_other_indels_in_context=int(((indel_pos >= seg_lo) & (indel_pos <= seg_hi) & (indel_pos != pos)).sum()),
                _haps=hs))
        print(f"[INFO] {locus}: {len(rows)} variants searched in {chrom}:{lo}-{hi}; {len(keep)} with r2 >= {args.r2} "
              f"({sum(1 for c in cands if c['locus'] == locus and c['type'] == 'SNP')} SNPs)")

    fa = out / "candidate_haps.fa"
    with open(fa, "w") as fh:
        for name, s in seqs:
            fh.write(f">{name}\n{s}\n")
    fimo = run([args.fimo, "--text", "--no-pgc", "--thresh", str(P_FLOOR), "--bfile", args.bfile, args.motifs, str(fa)])
    seqd = dict(seqs)
    best = {}
    for line in fimo.splitlines():
        if line.startswith(("motif_id", "#")) or not line.strip():
            continue
        f = line.split("\t")
        name, start, stop = f[2], int(f[3]), int(f[4])
        frag = seqd[name][start - 1:stop]
        if f[5] == "-":
            frag = frag.translate(COMP)[::-1]
        assert frag == f[9].upper(), f"matched_sequence mismatch: {line}"
        cid, grp = int(name[1:6]), name[7:]
        _, off, flen = cands[cid]["_haps"][grp]
        if stop < off or start > off + flen - 1:
            continue
        k = (cid, grp, f[0])
        best[k] = min(best.get(k, 1.0), float(f[7]))

    motifs_by_cand = collections.defaultdict(set)
    for (cid, _, m) in best:
        motifs_by_cand[cid].add(m)
    changes = []
    motif_names = {}
    for line in open(args.motifs):
        if line.startswith("MOTIF"):
            x = line.split()
            motif_names[x[1]] = x[2]
    for c in cands:
        # sanity: the allele sits at the footprint of each haplotype
        for grp, allele in (("BP", c["bp_allele"]), ("TP", c["tp_allele"])):
            s, off, _ = c["_haps"][grp]
            assert s[off - 1: off - 1 + len(allele)] == allele, f"allele not at footprint for {c['var_id']}"
        motifs = motifs_by_cand.get(c["cid"], set())
        for m in sorted(motifs):
            pb, pt = best.get((c["cid"], "BP", m), P_FLOOR), best.get((c["cid"], "TP", m), P_FLOOR)
            dl = math.log10(pt) - math.log10(pb)  # > 0: stronger on the BP allele
            if min(pb, pt) < P_PRESENT and abs(dl) >= MIN_FOLD_LOG10:
                changes.append(dict(locus=c["locus"], var_id=c["var_id"], pos=c["pos"], type=c["type"],
                                    r2_index=c["r2_index"], gemma_log_p=c["gemma_log_p"], annotation=c["annotation"],
                                    nearest_gene=c["nearest_gene"], motif_id=m, motif_name=motif_names[m],
                                    tf_family=fam.loc[m, "tf_family"], named_group=group_of.get(m, ""),
                                    p_bp=pb, p_tp=pt, delta_log10p=round(dl, 3),
                                    stronger_on="BP" if dl > 0 else "TP",
                                    present_on=("both" if max(pb, pt) < P_PRESENT else "BP only" if pb < P_PRESENT
                                                else "TP only")))
    ch = pd.DataFrame(changes)
    cdf = pd.DataFrame([{k: v for k, v in c.items() if not k.startswith("_")} for c in cands])
    if len(ch):
        agg = ch.groupby("var_id").agg(n_strong_changes=("motif_id", "size"),
                                      n_stronger_on_bp=("stronger_on", lambda x: int((x == "BP").sum())),
                                      max_abs_delta=("delta_log10p", lambda x: float(np.abs(x).max())),
                                      named_groups=("named_group", lambda x: ";".join(sorted({g for g in x if g}))))
        cdf = cdf.merge(agg, left_on="var_id", right_index=True, how="left")
    for col in ["n_strong_changes", "n_stronger_on_bp", "max_abs_delta"]:
        cdf[col] = cdf[col].fillna(0) if col in cdf else 0
    cdf[["n_strong_changes", "n_stronger_on_bp"]] = cdf[["n_strong_changes", "n_stronger_on_bp"]].astype(int)
    cdf = cdf.sort_values(["locus", "r2_index", "max_abs_delta"], ascending=[True, False, False])

    lines = ["TFBS candidate variants (prioritisation, no enrichment test)",
             f"LD: r2 >= {args.r2} with the index SNP in {len(samples)} BP+TP fish (call rate >= {args.min_call}); "
             f"FIMO p floor {P_FLOOR}; "
             f"strong change = best p < {P_PRESENT} and >= 10-fold difference", ""]
    for locus in LOCI:
        c = cdf[cdf.locus == locus].drop(columns=["cid"])
        x = ch[ch.locus == locus] if len(ch) else ch
        c.to_csv(out / f"candidates_{locus}.tsv", sep="\t", index=False)
        (x.sort_values(["r2_index", "delta_log10p"], ascending=[False, False]) if len(x) else x).to_csv(
            out / f"motif_changes_{locus}.tsv", sep="\t", index=False)
        # Long-format table for Supplementary Data S4-S6: one row per candidate x strong change; candidates without
        # a strong change keep one row with empty motif fields.
        mcols = ["motif_id", "motif_name", "tf_family", "named_group", "p_bp", "p_tp", "delta_log10p", "present_on"]
        base = c.drop(columns=["n_strong_changes", "n_stronger_on_bp", "max_abs_delta", "named_groups"], errors="ignore")
        long = base.merge(x[["var_id"] + mcols] if len(x) else pd.DataFrame(columns=["var_id"] + mcols),
                          on="var_id", how="left")
        long["_abs"] = long["delta_log10p"].abs()
        long = long.sort_values(["r2_index", "pos", "_abs"], ascending=[False, True, False]).drop(columns="_abs")
        long.to_csv(out / f"data_table_{locus}.tsv", sep="\t", index=False)
        n_chg = int((c.n_strong_changes > 0).sum())
        lines.append(f"{locus} (peak {LOCI[locus]['peak_no']}, index {LOCI[locus]['chr']}:{LOCI[locus]['index']}): "
                     f"{len(c)} candidates ({(c.type == 'SNP').sum()} SNPs, {(c.type != 'SNP').sum()} indels); "
                     f"{n_chg} with >= 1 strong motif change; {len(x)} strong changes "
                     f"({int((x.stronger_on == 'BP').sum()) if len(x) else 0} stronger on BP)")
        lines.append("   annotation: " + ", ".join(f"{k} {v}" for k, v in c.annotation.value_counts().items()))
        if len(x):
            ng = x[x.named_group != ""]
            lines.append("   named groups: " + ("; ".join(f"{r.motif_name} at {r.var_id} ({r.stronger_on})"
                                                         for r in ng.itertuples()) or "none"))
        top = c[c.n_strong_changes > 0].head(8)
        for r in top.itertuples():
            xs = x[x.var_id == r.var_id].sort_values("delta_log10p", key=np.abs, ascending=False).head(4)
            mot = ", ".join(f"{m.motif_name} {m.present_on} ({m.delta_log10p:+.1f})" for m in xs.itertuples())
            lines.append(f"   {r.var_id} {r.type} r2 {r.r2_index:.2f} logp {r.gemma_log_p} {r.annotation} "
                         f"{r.nearest_gene} TSS {r.dist_to_tss:+d}: {mot}")
        lines.append("")
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
