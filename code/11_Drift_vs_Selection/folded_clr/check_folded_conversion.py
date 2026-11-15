#!/usr/bin/env python3
"""A3 precondition: verify the folded SweepFinder2 conversion is genuinely polarization-free
and did not inherit a REF-as-ancestral bug.  Checks P1-P5 of METHODS.md on one chromosome.
Writes a TSV of the checks; exits 1 if any required check fails."""
import subprocess, sys, os

# Paths derive from the repository root, exported by common.sh (which sources config.env).
SD     = os.environ.get("SD")
SF2DIR = os.environ.get("SF2DIR")
if not SD or not SF2DIR:
    raise SystemExit("SD and SF2DIR are not set: source "
                     "code/11_Drift_vs_Selection/folded_clr/common.sh before running this check")
POP    = sys.argv[1]
CHROM  = sys.argv[2]
OUTTSV = sys.argv[3]

vcf  = f"{SD}/by_chrom/{POP}.polarized.{CHROM}.vcf.gz"
fold = f"{SF2DIR}/sfs_input_folded/{POP}.{CHROM}.sf2.txt"
unf  = f"{SF2DIR}/sfs_input_unfolded/{POP}.{CHROM}.sf2.txt"

def read_sf2(path):
    d = {}
    flags = set()
    with open(path) as fh:
        hdr = fh.readline()
        for line in fh:
            p = line.split()
            d[int(p[0])] = (int(p[1]), int(p[2]))
            flags.add(p[3])
    return d, flags

F, fflags = read_sf2(fold)
U, uflags = read_sf2(unf)

# truth from the VCF, recomputed independently of the converter
q = subprocess.run(["bcftools", "query", "-r", CHROM,
                    "-f", "%POS\t%REF\t%ALT\t%INFO/AA\t%INFO/AC\t%INFO/AN\n", vcf],
                   capture_output=True, text=True, check=True)
truth_minor, truth_derived_folded = {}, {}
n_aa_alt = n_aa_ref = n_aa_missing = 0
for line in q.stdout.splitlines():
    pos, ref, alt, aa, ac, an = line.split("\t")
    if aa in (".", "N", ""):
        n_aa_missing += 1
        continue
    if an in (".", "0", "") or ac in (".", ""):
        continue
    try:
        ac_i, an_i = int(ac), int(an)
    except ValueError:
        continue
    if aa == ref:
        d = ac_i; n_aa_ref += 1
    elif aa == alt:
        d = an_i - ac_i; n_aa_alt += 1
    else:
        continue
    p = int(pos)
    truth_minor[p] = min(ac_i, an_i - ac_i)
    truth_derived_folded[p] = min(d, an_i - d)

rows = []
def add(name, ok, detail, required=True):
    rows.append((name, "PASS" if ok else "FAIL", detail, "required" if required else "diagnostic"))
    return ok

# P1 -- folded count equals the minor-allele count computed straight from AC/AN (no AA involved)
bad1 = [p for p, (x, n) in F.items() if truth_minor.get(p) is not None and x != truth_minor[p]]
miss1 = [p for p in F if p not in truth_minor]
ok1 = add("P1_folded_equals_minor_AC_over_AN", len(bad1) == 0 and len(miss1) == 0,
          f"{len(F)} folded sites; {len(bad1)} mismatched; {len(miss1)} absent from VCF recomputation")

# P2 -- folded count equals min(d, AN-d) with d from the AA rule
bad2 = [p for p, (x, n) in F.items() if truth_derived_folded.get(p) is not None and x != truth_derived_folded[p]]
ok2 = add("P2_folded_equals_fold_of_derived", len(bad2) == 0,
          f"{len(bad2)} of {len(F)} folded sites disagree with min(d, AN-d)")

# P3 -- identical emitted position set folded vs unfolded
only_f = set(F) - set(U); only_u = set(U) - set(F)
ok3 = add("P3_identical_site_set_vs_unfolded", len(only_f) == 0 and len(only_u) == 0,
          f"folded {len(F)}, unfolded {len(U)}, folded-only {len(only_f)}, unfolded-only {len(only_u)}")

# P4 -- the AA call is not a no-op: some sites have AA == ALT, i.e. REF is the derived allele
ok4 = add("P4_AA_equals_ALT_sites_exist", n_aa_alt > 0,
          f"AA==REF {n_aa_ref}, AA==ALT {n_aa_alt}, AA missing {n_aa_missing}")

# P5 -- the folded flag column
ok5 = add("P5_folded_flag_column", fflags == {"1"} and uflags == {"0"},
          f"folded col4 values {sorted(fflags)}; unfolded col4 values {sorted(uflags)}")

# diagnostics
shared = set(F) & set(U)
changed = sum(1 for p in shared if F[p][0] != U[p][0])
add("D1_sites_where_folding_changed_the_count", True,
    f"{changed} of {len(shared)} shared sites ({100.0*changed/max(len(shared),1):.2f}%) have folded x != unfolded x",
    required=False)
an_mismatch = sum(1 for p in shared if F[p][1] != U[p][1])
add("D2_allele_number_identical", an_mismatch == 0,
    f"{an_mismatch} shared sites with a different AN", required=False)

with open(OUTTSV, "w") as o:
    o.write("check\tresult\tdetail\tkind\tpopulation\tchromosome\n")
    for n, res, det, kind in rows:
        o.write(f"{n}\t{res}\t{det}\t{kind}\t{POP}\t{CHROM}\n")

print(open(OUTTSV).read())
if not all([ok1, ok2, ok3, ok4, ok5]):
    print("PRECONDITION FAILED -- folded scan not run", file=sys.stderr)
    sys.exit(1)
print("PRECONDITION PASSED")
