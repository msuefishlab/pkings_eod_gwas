#!/usr/bin/env python3
"""records_full.py -- model extension, Stage 3: per-record site table of one 1-Mb M-sum block.

Derived from the coordinate pilot's records_mask.py: SAME bcftools command chain, SAME A1w callable
rule, SAME missing-genotype rule (GT not a diagnosed diploid code after +setGT DP<1 -> missing), SAME first-record-per-
position de-duplication. Added per record x sample (needed for the strict-genotype and heterozygote-dropout sensitivity rows): FMT/DP, FMT/GQ (variant
records) and FMT/RGQ (invariant records, GATK's GQ analogue).

Per block writes <out_dir>/<chrom>_<block>.rec.npz:
  rpos     (n,) int64   record positions (sorted, unique)
  code     (n, 114) uint8  bits 0-2 DP bin (0: DP 0-3, 1: 4-7, 2: 8-11, 3: 12-19, 4: >= 20; DP '.' -> 0)
                           bit 3 missing (drift-test / simulation 23a rule)
                           bit 4 missing under the strict-genotype rule: missing OR DP < 8 OR Q < 20, Q = GQ at variant records,
                                 RGQ at invariant records; Q '.' counts as failing
  callable (n,) bool, is_snp (n,) bool (ALT != '.'), pol (n,) int8 (as records_mask.py)
  vrow     (m,) int64 row index of the variant records (ALT != '.'); gq (m, 114) uint8 GQ (255 = '.')
Checks (<chrom>_<block>.check.tsv): callable records = L_scan; all Stage 1 block-npz SNVs found with identical callable
flag and identical per-genotype missingness; if a pilot records npz exists for the block: rpos / miss / callable identical.
Usage: records_full.py <vcf> <groups.tsv> <group_sizes.tsv> <samples_union.txt> <deserts.tsv> <blocks_npz_dir>
                       <pilot_records_dir> <chrom> <block> <L_expected> <out_dir>
"""
import os, subprocess, sys, time
import numpy as np

vcf, groups_fn, gsizes_fn, union_fn, des_fn, bdir, pilot_dir, ch, b, L_exp, out_dir = sys.argv[1:12]
b, L_exp = int(b), int(L_exp)
s, e = b * 1_000_000, (b + 1) * 1_000_000
N_BP, N_TP, CAP, EXCHET_MIN = 60, 60, 60, 0.001
OUT_MIN_AN, OUT_MIN_SPECIES = 8, 2
NI = 114
t0 = time.time()

groups = [l.split("\t")[0] for l in open(gsizes_fn).read().splitlines()[1:] if l.strip()]
exchet = [g for g in groups if g.startswith("w") and g not in ("wTP1", "wTP2", "wBP2")]
outg = [g for g in groups if g.startswith("OUT_")]
union = [l.strip() for l in open(union_fn) if l.strip()]
assert len(union) == NI
des = set()
for line in open(des_fn):
    f = line.rstrip("\n").split("\t")
    if f[0] == "A1w":
        des.add((f[1], int(f[2])))

tmpd = os.environ.get("TMPDIR", "/tmp")
allsamp = os.path.join(tmpd, f"s3rec_allsamples_{ch}_{b}.txt")
with open(allsamp, "w") as fo:
    fo.write("\n".join(l.split("\t")[0] for l in open(groups_fn).read().splitlines() if l.strip()) + "\n")
a1w = os.path.join(tmpd, f"s3rec_a1w_{ch}_{b}.txt")
with open(a1w, "w") as fo:
    fo.write("\n".join(union) + "\n")
fmt = "%POS\\t%REF\\t%ALT\\t%INFO/AN_wBAM\\t%INFO/AN_wTP1"
for g in outg:
    fmt += f"\\t%INFO/AN_{g}\\t%INFO/AC_{g}"
for g in exchet:
    fmt += f"\\t%INFO/ExcHet_{g}"
fmt += "[\\t%GT][\\t%DP][\\t%GQ][\\t%RGQ]\\n"
e_q = s + int(os.environ["PROBE_BP"]) if os.environ.get("PROBE_BP") else e    # parser smoke test only
cmd = (f"set -o pipefail; bcftools view -r '{ch}:{s + 1}-{e_q}' -S '{allsamp}' -Ou '{vcf}' "
       f"| bcftools +setGT -Ou -- -t q -n . -i 'FMT/DP<1' "
       f"| bcftools +fill-tags -Ou -- -S '{groups_fn}' -t AN,AC,ExcHet "
       f"| bcftools view -I -S '{a1w}' -Ou "
       f"| bcftools query -f '{fmt}'")
k0 = 5 + 2 * len(outg) + len(exchet)
hcmd = f"set -o pipefail; bcftools view -h -S '{allsamp}' '{vcf}' | bcftools view -h -I -S '{a1w}' - | tail -1"
names = subprocess.run(["bash", "-c", hcmd], capture_output=True, text=True, check=True).stdout.rstrip("\n").split("\t")[9:]
assert sorted(names) == sorted(union), (len(names), names[:3])
col = np.array([names.index(x) for x in union])
CODE = {"0/0", "0|0", "0/1", "0|1", "1|0", "1/0", "1/1", "1|1"}

cap = 1_100_000
code = np.zeros((cap, NI), np.uint8)
gq_rows, vrow = [], []
cpos = np.zeros(cap, np.int64); csnp = np.zeros(cap, bool); cpol = np.zeros(cap, np.int8); ccal = np.zeros(cap, bool)
st = dict(n_records=0, n_outside_block=0, n_invariant=0, n_snv=0, n_skipped_non_snv=0, n_callable_invariant=0,
          n_callable_snv=0, n_q_missing_called=0)
p = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
to_int = lambda x: 0 if x in (".", "") else int(x)
BINS = np.array([4, 8, 12, 20])
n = 0
for line in p.stdout:
    f = line.rstrip("\n").split("\t")
    pos = int(f[0])
    if not (s < pos <= e):
        st["n_outside_block"] += 1
        continue
    st["n_records"] += 1
    ref, alt = f[1], f[2]
    Nb, Nt = to_int(f[3]), to_int(f[4])
    spp = []
    for i in range(len(outg)):
        an = to_int(f[5 + 2 * i])
        if an >= 2:
            spp.append((an, to_int(f[6 + 2 * i])))
    out_ok = len(spp) >= OUT_MIN_SPECIES and sum(a for a, _ in spp) >= OUT_MIN_AN
    aa = "."
    if out_ok:
        if all(ac == 0 for _, ac in spp):
            aa = "R"
        elif all(ac == an for an, ac in spp):
            aa = "A"
    ok_bin = (ch, (pos - 1) // 100_000) not in des
    if alt == ".":
        st["n_invariant"] += 1
        cal = out_ok and min(Nb, CAP) >= N_BP and min(Nt, CAP) >= N_TP and ok_bin
        snp = False
    elif len(ref) != 1 or len(alt) != 1:
        st["n_skipped_non_snv"] += 1
        cal = False
        snp = True
    else:
        st["n_snv"] += 1
        ex = f[5 + 2 * len(outg):k0]
        ok_exc = not any(x not in (".", "") and float(x) < EXCHET_MIN for x in ex)
        cal = ok_exc and aa != "." and Nb >= N_BP and Nt >= N_TP and ok_bin
        snp = True
    if cal:
        st["n_callable_snv" if snp else "n_callable_invariant"] += 1
    gt = f[k0:k0 + NI]
    dp = np.array([to_int(x) for x in f[k0 + NI:k0 + 2 * NI]], np.int32)[col]
    qraw = f[k0 + 2 * NI:k0 + 3 * NI] if snp else f[k0 + 3 * NI:k0 + 4 * NI]
    q = np.array([-1 if x in (".", "") else int(x) for x in qraw], np.int32)[col]
    miss = np.fromiter((g not in CODE for g in gt), bool, NI)[col]
    st["n_q_missing_called"] += int(np.sum((q < 0) & ~miss))
    d019 = miss | (dp < 8) | (q < 20)
    if n >= cap:
        cap2 = cap + 200_000
        code = np.vstack([code, np.zeros((cap2 - cap, NI), np.uint8)])
        cpos, csnp, cpol, ccal = (np.concatenate([a, np.zeros(cap2 - cap, a.dtype)]) for a in (cpos, csnp, cpol, ccal))
        cap = cap2
    code[n] = np.searchsorted(BINS, dp, side="right").astype(np.uint8) | (miss.astype(np.uint8) << 3) | \
        (d019.astype(np.uint8) << 4)
    cpos[n] = pos; csnp[n] = snp; cpol[n] = 1 if aa == "R" else (-1 if aa == "A" else 0); ccal[n] = cal
    if snp:
        g = np.array([255 if x in (".", "") else min(int(x), 254) for x in f[k0 + 2 * NI:k0 + 3 * NI]], np.uint8)[col]
        gq_rows.append(g); vrow.append(n)
    n += 1
rc = p.wait()
assert rc == 0, f"bcftools pipeline exit {rc}"
code, cpos, csnp, cpol, ccal = code[:n], cpos[:n], csnp[:n], cpol[:n], ccal[:n]
gq = np.vstack(gq_rows) if gq_rows else np.zeros((0, NI), np.uint8)
vrow = np.array(vrow, np.int64)
n_cal = int(ccal.sum())
order = np.argsort(cpos, kind="stable")
inv = np.empty_like(order); inv[order] = np.arange(len(order))
cpos, code, csnp, cpol, ccal = cpos[order], code[order], csnp[order], cpol[order], ccal[order]
_, first = np.unique(cpos, return_index=True)
keepmask = np.zeros(len(cpos), bool); keepmask[first] = True
newidx = np.cumsum(keepmask) - 1
rpos, codeu, snpu, polu, calu = cpos[first], code[first], csnp[first], cpol[first], ccal[first]
vr_sorted = inv[vrow]                                   # sorted-row index of each variant record
vk = keepmask[vr_sorted]
vrow_u = newidx[vr_sorted[vk]]; gq_u = gq[vk]
o2 = np.argsort(vrow_u, kind="stable"); vrow_u, gq_u = vrow_u[o2], gq_u[o2]
missu = (codeu >> 3) & 1

# cross-check against the Stage 1 block npz
z = np.load(os.path.join(bdir, f"{ch}_{b}.npz"))
zn = list(z["samples"]); zc = np.array([zn.index(x) for x in union])
zpos, zcal = z["pos"], z["callable"]
zG = z["G"][:, zc] if z["G"].size else np.zeros((0, NI), np.int8)
j = np.clip(np.searchsorted(rpos, zpos), 0, max(len(rpos) - 1, 0))
hit = (rpos[j] == zpos) if len(rpos) else np.zeros(len(zpos), bool)
chk = dict(chrom=ch, block=b, L_expected=L_exp, n_records_kept=n, n_record_positions=len(rpos),
           n_callable_records=n_cal, n_callable_unique=int(calu.sum()), L_match=int(int(calu.sum()) == L_exp),
           npz_snps=len(zpos), npz_snps_found=int(hit.sum()),
           callable_flag_match=int(np.sum(calu[j][hit] == zcal[hit])),
           geno_compared=int(hit.sum()) * NI, geno_miss_match=int(np.sum(missu[j][hit].astype(bool) == (zG[hit] < 0))),
           npz_snps_have_gq=int(np.isin(j[hit], vrow_u).sum()),
           miss_rate_records=float(missu.mean()) if len(missu) else float("nan"),
           d019_rate_records=float(((codeu >> 4) & 1).mean()) if len(codeu) else float("nan"),
           seconds=0.0, **st)
pf = os.path.join(pilot_dir, f"{ch}_{b}.records.npz")
if os.path.exists(pf):
    zp = np.load(pf)
    pm = np.unpackbits(zp["miss"], axis=1, count=NI).astype(bool)
    chk["pilot_compared"] = 1
    chk["pilot_rpos_identical"] = int(np.array_equal(zp["rpos"], rpos))
    chk["pilot_miss_identical"] = int(chk["pilot_rpos_identical"] and np.array_equal(pm, missu.astype(bool)))
    chk["pilot_callable_identical"] = int(chk["pilot_rpos_identical"] and np.array_equal(zp["callable"], calu))
else:
    chk["pilot_compared"] = 0; chk["pilot_rpos_identical"] = chk["pilot_miss_identical"] = chk["pilot_callable_identical"] = -1
chk["seconds"] = round(time.time() - t0, 1)
os.makedirs(out_dir, exist_ok=True)
tmp = os.path.join(out_dir, f"{ch}_{b}.rec.tmp.npz")
np.savez_compressed(tmp, rpos=rpos, code=codeu, callable=calu, is_snp=snpu, pol=polu, vrow=vrow_u, gq=gq_u,
                    samples=np.array(union), n_callable_records=int(calu.sum()))
os.replace(tmp, os.path.join(out_dir, f"{ch}_{b}.rec.npz"))
with open(os.path.join(out_dir, f"{ch}_{b}.check.tsv.tmp"), "w") as fo:
    fo.write("\t".join(chk) + "\n" + "\t".join(str(v) for v in chk.values()) + "\n")
os.replace(os.path.join(out_dir, f"{ch}_{b}.check.tsv.tmp"), os.path.join(out_dir, f"{ch}_{b}.check.tsv"))
for fn in (allsamp, a1w):
    os.remove(fn)
print(chk, flush=True)
