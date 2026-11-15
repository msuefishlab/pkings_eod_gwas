"""Tests for matched_null.py (TFBS Stage B, README.md). Synthetic fixtures only.

Run any of:
  python3 code/15_TFBS_Scan/tests/test_tfbs_null.py
  python3 -m unittest discover -s code/15_TFBS_Scan/tests
  pytest code/15_TFBS_Scan/tests   (if pytest is installed)
"""
import contextlib
import io
import shutil
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import matched_null as mn  # noqa: E402
import synth_fixtures as sf  # noqa: E402

N_DRAWS = 2000
N_CALIB = 500


def run_fixture(paths, out_dir, units=None, polarized=True, extra=()):
    argv = ["--units", str(units or paths["units"]), "--hits", str(paths["hits"]), "--families",
            str(paths["families"]), "--out_dir", str(out_dir), "--n_draws", str(N_DRAWS), "--n_calib", str(N_CALIB),
            "--min_stratum", "20", "--no_plots"]
    if polarized and paths.get("polarized"):
        argv += ["--polarized", str(paths["polarized"])]
    args = mn.parse_args(argv + list(extra))
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return mn.run(args)


def lt_row(res, s, rule, scope, stat="f", null="matched"):
    lt = res["locus_tests"]
    r = lt[(lt.set == s) & (lt.rule == rule) & (lt.scope == scope) & (lt.stat == stat) & (lt.null == null)]
    assert len(r) == 1, (s, rule, scope, stat, null)
    return r.iloc[0]


class TmpDirCase(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="tfbs_null_test_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)


# ---------------------------------------------------------------------------------------------------------------------
# (1) BP/TP label swap at the peak: f -> 1 - f, family (and named-group) nets negate, p-values unchanged
# ---------------------------------------------------------------------------------------------------------------------
class TestLabelSwap(TmpDirCase):
    def test_swap_mirrors_everything(self):
        paths, *_ = sf.make_fixture(self.tmp / "fx", seed=3)
        a = run_fixture(paths, self.tmp / "a", polarized=False)
        swapped = self.tmp / "fx" / "units_swapped.tsv"
        sf.swap_bp_tp(paths["units"], swapped)
        b = run_fixture(paths, self.tmp / "b", units=swapped, polarized=False)
        for s in a["sets"]:
            for rule in mn.RULES:
                for sc in a["loci"] + ["pooled"]:
                    ra, rb = lt_row(a, s, rule, sc), lt_row(b, s, rule, sc)
                    self.assertAlmostEqual(rb.observed, 1 - ra.observed, places=12)
                    self.assertEqual(ra.peak_nBP, rb.peak_nTP)
                    self.assertEqual(ra.peak_nTP, rb.peak_nBP)
                    self.assertAlmostEqual(ra.p_two_sided, rb.p_two_sided, places=12)
                    fa, fb = a["null_f"][(s, rule, sc)][0], b["null_f"][(s, rule, sc)][0]
                    np.testing.assert_allclose(fb, 1 - fa, atol=1e-12)          # the matched null mirrors exactly
                    sa, sb = lt_row(a, s, rule, sc, "f", "sign_flip"), lt_row(b, s, rule, sc, "f", "sign_flip")
                    self.assertAlmostEqual(sa.p_two_sided, sb.p_two_sided, places=12)
                    ea, eb = lt_row(a, s, rule, sc, "events_per_snp"), lt_row(b, s, rule, sc, "events_per_snp")
                    self.assertAlmostEqual(ea.observed, eb.observed, places=12)
                    oa, na = a["fam_null"][(s, rule, sc)]
                    ob, nb = b["fam_null"][(s, rule, sc)]
                    np.testing.assert_array_equal(ob, -oa)                      # family nets negate
                    np.testing.assert_array_equal(nb, -na)
                    ga, gb = a["A"][(s, sc, "peak")][rule], b["A"][(s, sc, "peak")][rule]
                    np.testing.assert_array_equal(gb["grpBP"] - gb["grpTP"], -(ga["grpBP"] - ga["grpTP"]))
        fa = a["family_tests"].set_index(["set", "rule", "scope", "tf_family"])
        fb = b["family_tests"].set_index(["set", "rule", "scope", "tf_family"]).loc[fa.index]
        np.testing.assert_array_equal(fb["peak_net"].to_numpy(), -fa["peak_net"].to_numpy())
        np.testing.assert_allclose(fb["p_two_sided"].to_numpy(), fa["p_two_sided"].to_numpy())
        # balance: bp_is_ref mirrors (peak and null)
        ba = a["balance"].set_index(["set", "scope", "covariate"])
        bb = b["balance"].set_index(["set", "scope", "covariate"])
        for col in ("peak", "null_mean"):
            np.testing.assert_allclose(bb.xs("bp_is_ref", level=2)[col], 1 - ba.xs("bp_is_ref", level=2)[col])


# ---------------------------------------------------------------------------------------------------------------------
# (2) Exchangeable synthetic null: calibration p ~ uniform, SD guard passes; a degenerate null trips the guard
# ---------------------------------------------------------------------------------------------------------------------
class TestCalibration(TmpDirCase):
    def test_exchangeable_null_is_calibrated(self):
        paths, *_ = sf.make_fixture(self.tmp / "fx", seed=5)
        res = run_fixture(paths, self.tmp / "out")
        cal = res["calibration"]
        prim = cal[(cal.set == "extended") & (cal.rule == "standard")]
        ks = prim[prim.check == "ks_uniform_randomized"]
        self.assertEqual(len(ks), 4)
        self.assertTrue((ks.p_value > 0.01).all(), ks[["scope", "stat", "p_value"]])
        self.assertTrue((ks.stat < 0.08).all())
        raw = prim[prim.check == "ks_uniform_raw"]
        self.assertTrue((raw.p_value > 0.001).all(), raw[["scope", "stat", "p_value"]])
        guard = cal[cal.check == "sd_guard"]
        self.assertTrue(guard[guard.passed.notna()].passed.astype(bool).all())
        self.assertGreater(lt_row(res, "extended", "standard", "pooled").p_two_sided, 0.01)
        for f in ["summary.txt", "locus_tests.tsv", "family_tests.tsv", "named_group_tests.tsv", "balance.tsv",
                  "calibration.tsv", "fallbacks.tsv", "per_unit_calls_sptbn4b.tsv", "per_unit_calls_cntnap5.tsv",
                  "per_unit_calls_itih3.tsv"]:
            self.assertTrue((self.tmp / "out" / f).exists(), f)

    def test_degenerate_null_trips_sd_guard(self):
        # every haplotype pair gives exactly one BP-only and one TP-only call -> f = 0.5 in every draw, SD 0
        paths, df, motifs, _ = sf.make_fixture(self.tmp / "fx", seed=2, n_pool1=1500, n_chain=150)
        k = df["k"].to_numpy()
        peak = (df["source"] == "peak").to_numpy()
        u, m = [], []
        for i in range(len(df)):
            if peak[i]:
                b = int(df["bp_mask"].iloc[i])
                masks = [b, b ^ (2 ** k[i] - 1)]
            else:
                masks = range(2 ** k[i])
            for mm in masks:
                u.append(i)
                m.append(mm)
        u, m = np.array(u), np.array(m)
        mot = np.where(m & 1, 0, 1)     # motif 0 when bit 0 is set, motif 1 otherwise
        sf.write_hits(dict(unit_idx=u, mask=m, motif_idx=mot, mlog10p=np.full(len(u), 6.0)), motifs, paths["hits"])
        with self.assertRaises(mn.CalibrationError):
            run_fixture(paths, self.tmp / "out")
        self.assertTrue((self.tmp / "out" / "calibration.tsv").exists())
        self.assertIn("ABORTED", (self.tmp / "out" / "summary.txt").read_text())
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()), \
                warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with self.assertRaises(SystemExit) as cm:
                mn.main(["--units", str(paths["units"]), "--hits", str(paths["hits"]), "--families",
                         str(paths["families"]), "--out_dir", str(self.tmp / "out2"), "--n_draws", "600",
                         "--n_calib", "100", "--min_stratum", "20", "--no_plots"])
        self.assertEqual(cm.exception.code, 2)

    def test_sd_guard_formula(self):
        g = mn.sd_guard(np.array([0.5, 0.5, 0.5]), np.array([100, 100, 100]))
        self.assertFalse(g["passed"])
        self.assertAlmostEqual(g["threshold"], 0.8 * 0.5 / 10)
        rng = np.random.default_rng(0)
        n = rng.poisson(100, 5000)
        f = rng.binomial(n, 0.5) / n
        self.assertTrue(mn.sd_guard(f, n)["passed"])


# ---------------------------------------------------------------------------------------------------------------------
# (3) Collapsed-ID failure mode and other loader checks
# ---------------------------------------------------------------------------------------------------------------------
class TestLoaders(TmpDirCase):
    def setUp(self):
        super().setUp()
        self.paths, self.df, self.motifs, self.h = sf.make_fixture(self.tmp / "fx", seed=4, n_pool1=500, n_chain=50)

    def test_duplicate_unit_idx_raises(self):
        u = pd.read_csv(self.paths["units"], sep="\t", dtype=str, keep_default_na=False)
        u.loc[5, "unit_idx"] = u.loc[4, "unit_idx"]
        bad = self.tmp / "dup_units.tsv"
        u.to_csv(bad, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "duplicate unit_idx"):
            mn.load_units(bad)

    def test_duplicate_hit_triple_raises(self):
        h = {key: np.r_[v, v[:1]] for key, v in self.h.items()}
        sf.write_hits(h, self.motifs, self.paths["hits"])
        units = mn.load_units(self.paths["units"])
        hits = mn.load_hits(self.paths["hits"], units)
        with self.assertRaisesRegex(ValueError, "duplicated"):
            mn.hits_to_rows(hits, units, mn.build_rows(units))

    def test_peak_hit_on_foreign_mask_raises(self):
        units = mn.load_units(self.paths["units"])
        big = int(np.where(units["is_peak"] & (units["k"] >= 2))[0][0])
        b = int(units["bp_mask"][big])
        other = [x for x in range(2 ** int(units["k"][big])) if x not in (b, b ^ (2 ** int(units["k"][big]) - 1))][0]
        h = {key: v.copy() for key, v in self.h.items()}
        h["unit_idx"] = np.r_[h["unit_idx"], big]
        h["mask"] = np.r_[h["mask"], other]
        h["motif_idx"] = np.r_[h["motif_idx"], 0]
        h["mlog10p"] = np.r_[h["mlog10p"], 5.0]
        sf.write_hits(h, self.motifs, self.paths["hits"])
        hits = mn.load_hits(self.paths["hits"], units)
        with self.assertRaisesRegex(ValueError, "out of sync"):
            mn.hits_to_rows(hits, units, mn.build_rows(units))

    def test_polarized_disagreement_raises(self):
        pol = pd.read_csv(self.paths["polarized"], sep="\t", dtype=str, keep_default_na=False)
        pol.loc[0, "bp_is_ref"] = "False" if pol.loc[0, "bp_is_ref"] == "True" else "True"
        bad = self.tmp / "bad_pol.tsv"
        pol.to_csv(bad, sep="\t", index=False)
        units = mn.load_units(self.paths["units"])
        with self.assertRaisesRegex(ValueError, "bp_mask disagrees"):
            mn.check_polarized(units, mn.load_polarized(bad))

    def test_family_table(self):
        units = mn.load_units(self.paths["units"])
        hits = mn.load_hits(self.paths["hits"], units)
        fam = mn.load_families(self.paths["families"], hits["motif_ids"], hits["motif_names"])
        self.assertEqual(fam["primary"][-1], "Fam00")       # "Fam00;Fam01" -> first family
        self.assertEqual(fam["primary"][-2], "Unknown")     # empty family
        self.assertIn("Unknown", fam["names"])


# ---------------------------------------------------------------------------------------------------------------------
# (4) Planted signal: 20 BP-only calls in one family -> q < 0.05; (4b) planted BCL6B -> BCL6 named group
# ---------------------------------------------------------------------------------------------------------------------
class TestPlantedSignal(TmpDirCase):
    def test_planted_family(self):
        motifs = sf.motif_table()
        fam_motifs = motifs.loc[motifs.tf_family == "Fam03", "motif_id"].tolist()
        spec = dict(locus="cntnap5", n=20, plant_motif=dict(motif=fam_motifs[0], clear=fam_motifs))
        paths, *_ = sf.make_fixture(self.tmp / "fx", seed=6, plant_spec=spec)
        res = run_fixture(paths, self.tmp / "out")
        ft = res["family_tests"]
        r = ft[(ft.set == "extended") & (ft.rule == "standard") & (ft.scope == "cntnap5") & (ft.tf_family == "Fam03")]
        self.assertEqual(len(r), 1)
        self.assertGreaterEqual(int(r.peak_nBP.iloc[0]), 20)
        self.assertLess(float(r.q_bh.iloc[0]), 0.05)
        # unplanted loci stay unremarkable for this family
        other = ft[(ft.set == "extended") & (ft.rule == "standard") & (ft.scope == "sptbn4b") & (ft.tf_family == "Fam03")]
        self.assertGreater(float(other.q_bh.iloc[0]), 0.05)

    def test_planted_bcl6b_named_group(self):
        motifs = sf.motif_table()
        bcl6 = mn.NAMED_GROUPS["BCL6"]
        spec = dict(locus="itih3", n=20, plant_motif=dict(motif="MA0731.1", clear=bcl6))
        paths, *_ = sf.make_fixture(self.tmp / "fx", seed=7, plant_spec=spec)
        res = run_fixture(paths, self.tmp / "out")
        ng = pd.read_csv(self.tmp / "out" / "named_group_tests.tsv", sep="\t")
        prim = ng[(ng.set == "extended") & (ng.rule == "standard")]
        r = prim[(prim.scope == "itih3") & (prim.group == "BCL6") & (prim.stat == "net")].iloc[0]
        self.assertGreaterEqual(r.observed, 20)
        self.assertLess(r.p_bonferroni, 0.05)
        self.assertAlmostEqual(r.p_bonferroni, min(1.0, 2 * r.p_two_sided))
        mem = prim[(prim.scope == "itih3") & (prim.stat == "member_count:MA0731.1 BCL6B")].iloc[0]
        self.assertGreaterEqual(mem.peak_nBP, 20)
        self.assertTrue(np.isnan(mem.p_two_sided))         # member counts are description only
        r_sp = prim[(prim.scope == "sptbn4b") & (prim.group == "BCL6") & (prim.stat == "net")].iloc[0]
        self.assertGreater(r_sp.p_bonferroni, 0.05)
        # named groups stay out of the family BH: the family table has no BCL6/CEBP group rows
        self.assertFalse(res["family_tests"].tf_family.isin(["BCL6", "CEBP"]).any())
        self.assertIn("4b. Named groups", (self.tmp / "out" / "summary.txt").read_text())

    def test_named_groups_constant(self):
        self.assertEqual(mn.NAMED_GROUPS["BCL6"], ["MA0463.3", "MA0731.1"])
        self.assertEqual(sorted(mn.NAMED_GROUPS["CEBP"]),
                         sorted(["MA0102.5", "MA0466.4", "MA0836.3", "MA0837.3", "MA0838.1", "MA1636.2", "MA0019.2"]))


# ---------------------------------------------------------------------------------------------------------------------
# (5) Strand-collapse class table, (6) k >= 2 orientation rule
# ---------------------------------------------------------------------------------------------------------------------
class TestClasses(unittest.TestCase):
    def test_directed_class_table(self):
        expected = {("C", "A"): "C>A", ("G", "T"): "C>A", ("C", "G"): "C>G", ("G", "C"): "C>G",
                    ("C", "T"): "C>T", ("G", "A"): "C>T", ("T", "A"): "T>A", ("A", "T"): "T>A",
                    ("T", "C"): "T>C", ("A", "G"): "T>C", ("T", "G"): "T>G", ("A", "C"): "T>G"}
        for (bp, tp), cls in expected.items():
            self.assertEqual(mn.directed_class(bp, tp), cls, (bp, tp))
            self.assertEqual(mn.CLASSES[mn.DC_TABLE[mn.BASE_CODE[bp], mn.BASE_CODE[tp]]], cls)
        self.assertEqual(mn.directed_class("G", "A"), "C>T")
        self.assertEqual(mn.directed_class("A", "C"), "T>G")
        rev = {c: mn.CLASSES[mn.REV_CLASS[i]] for i, c in enumerate(mn.CLASSES)}
        self.assertEqual(rev, {"C>A": "T>G", "C>G": "C>G", "C>T": "T>C", "T>A": "T>A", "T>C": "C>T", "T>G": "C>A"})
        with self.assertRaises(ValueError):
            mn.directed_class("A", "A")

    def test_orientation_rule_k2plus(self):
        c = mn.BASE_CODE
        # peak unit, 3 SNPs: BP->TP = G>A (C>T), C>G (self-complementary), A>C (T>G); bits from bp_mask
        peak_bp = np.array([c["G"], c["C"], c["A"]])
        peak_tp = np.array([c["A"], c["G"], c["C"]])
        peak_bits = np.array([1, 0, 1])
        peak_dc = [mn.DC_TABLE[b, t] for b, t in zip(peak_bp, peak_tp)]
        rng = np.random.default_rng(0)
        ref = rng.integers(0, 4, (2000, 3))
        alt = (ref + rng.integers(1, 4, (2000, 3))) % 4
        masks = mn.orient_chain_masks(ref, alt, peak_bp, peak_tp, peak_bits)
        n_match = 0
        for i in range(len(ref)):
            for j in range(3):
                bit = (masks[i] >> j) & 1
                bp, tp = (alt[i, j], ref[i, j]) if bit else (ref[i, j], alt[i, j])
                pool_dc = mn.DC_TABLE[bp, tp]
                if mn.UND_CLASS[mn.DC_TABLE[ref[i, j], alt[i, j]]] == mn.UND_CLASS[peak_dc[j]]:
                    self.assertEqual(pool_dc, peak_dc[j])          # reproduces the peak directed class
                    if not mn.SELF_COMP[peak_dc[j]]:
                        n_match += 1
                    else:
                        self.assertEqual(bit, peak_bits[j])       # C>G / T>A: copy bp_is_ref
                else:
                    self.assertEqual(bit, peak_bits[j])           # no match: copy bp_is_ref
        self.assertGreater(n_match, 100)
        # explicit case: pool C>T SNP read REF->ALT as T>C must be flipped to reproduce the peak's C>T
        m = mn.orient_chain_masks(np.array([[c["T"], c["C"]]]), np.array([[c["C"], c["A"]]]),
                                  np.array([c["C"], c["C"]]), np.array([c["T"], c["A"]]), np.array([0, 1]))
        self.assertEqual(int(m[0]) & 1, 1)          # BP = ALT (C), TP = REF (T): C>T
        self.assertEqual((int(m[0]) >> 1) & 1, 0)   # C>A peak vs pool C>A as REF->ALT: BP = REF

    def test_k1_pool_orientations_are_two_strata(self):
        paths, *_ = sf.make_fixture(Path(tempfile.mkdtemp(prefix="tfbs_k1_")), seed=8, n_pool1=300, n_chain=30)
        units = mn.load_units(paths["units"])
        rows = mn.build_rows(units)
        st = mn.build_strata(units, rows)
        keys, uniq, start, cnt, mem = st["idx1"].tables[-1]         # (class, bp_is_ref)
        self.assertEqual(int(cnt.sum()), 2 * st["n_pool1"])         # each singleton once per orientation
        u = rows["row_unit"][mem]
        for i in range(len(uniq)):
            self.assertEqual(len(np.unique(u[start[i]:start[i] + cnt[i]])), cnt[i])
        shutil.rmtree(paths["dir"], ignore_errors=True)


# ---------------------------------------------------------------------------------------------------------------------
# (7) p-value formula and BH
# ---------------------------------------------------------------------------------------------------------------------
class TestPValues(unittest.TestCase):
    def test_two_sided_formula(self):
        null = np.arange(1, 11, dtype=float)                       # N = 10
        self.assertAlmostEqual(mn.two_sided_p(10, null), 2 * 2 / 11)   # upper tail (1+1)/11
        self.assertAlmostEqual(mn.two_sided_p(100, null), 2 * 1 / 11)  # nothing as extreme: (0+1)/11
        self.assertAlmostEqual(mn.two_sided_p(1, null), 2 * 2 / 11)
        self.assertEqual(mn.two_sided_p(5.5, null), 1.0)              # 2 * 6/11 > 1 -> capped at 1
        self.assertAlmostEqual(mn.two_sided_p(3, null), 2 * 4 / 11)
        p = mn.two_sided_p(np.array([10.0, 3.0]), np.c_[null, null])
        np.testing.assert_allclose(p, [4 / 11, 8 / 11])
        self.assertAlmostEqual(mn.two_sided_p(10, np.r_[null, np.nan]), 4 / 11)   # NaN draws leave N
        self.assertTrue(np.isnan(mn.two_sided_p(np.nan, null)))

    def test_bh_hand_example(self):
        q = mn.bh_fdr(np.array([0.01, 0.04, 0.03, 0.005]))
        np.testing.assert_allclose(q, [0.02, 0.04, 0.04, 0.02])
        q = mn.bh_fdr(np.array([0.01, 0.011, 0.5]))                 # needs the monotonicity step
        np.testing.assert_allclose(q, [0.0165, 0.0165, 0.5])
        q = mn.bh_fdr(np.array([0.2, np.nan, 0.9, 0.95]))
        np.testing.assert_allclose(q, [0.6, np.nan, 0.95, 0.95])
        self.assertTrue((mn.bh_fdr(np.array([0.9, 0.99])) <= 1).all())

    def test_calibration_pvalues_uniform(self):
        rng = np.random.default_rng(1)
        null = rng.binomial(40, 0.5, 4000) / 40.0                  # heavy ties
        p_raw, p_rand = mn.calibration_pvalues(null, 1000, np.random.default_rng(2))
        from scipy import stats
        self.assertGreater(stats.kstest(p_rand, "uniform").pvalue, 0.01)
        self.assertTrue((p_raw >= p_rand - 1e-12).mean() > 0.9)    # raw (inclusive ties) is conservative


# ---------------------------------------------------------------------------------------------------------------------
# (8) Fallback order af -> gc -> cpg; sampling without replacement
# ---------------------------------------------------------------------------------------------------------------------
class TestFallback(unittest.TestCase):
    @staticmethod
    def _index(groups):
        f = {key: [] for key in ["cls", "cpg", "gcq", "afq", "bpref"]}
        for (cls, cpg, gcq, afq, bpref), n in groups:
            for key, v in zip(f, (cls, cpg, gcq, afq, bpref)):
                f[key] += [v] * n
        n = len(f["cls"])
        return mn.StratumIndex(np.arange(n), {key: np.array(v) for key, v in f.items()}, mn.RADIX, mn.K1_LEVELS)

    def test_order_af_then_gc_then_cpg(self):
        q = dict(cls=2, cpg=0, gcq=1, afq=1, bpref=1)
        # full 10; same gc other af 30 (drop_af 40); other gc same af 100. Dropping gc first would give 110 at the
        # second level; the specified order gives drop_af 40 (< 50) and then drop_af_gc 140.
        idx = self._index([((2, 0, 1, 1, 1), 10), ((2, 0, 1, 3, 1), 30), ((2, 0, 4, 1, 1), 100),
                           ((2, 0, 1, 1, 0), 500), ((2, 1, 1, 1, 1), 500)])
        lev, mem, _, sizes, below = idx.lookup(q, 50)
        self.assertEqual(mn.K1_LEVELS[lev][0], "drop_af_gc")
        self.assertEqual(sizes, [10, 40, 140])
        self.assertEqual(len(mem), 140)
        self.assertFalse(below)
        # cpg is dropped last: drop_af_gc 45 (< 50) -> drop_af_gc_cpg 145; dropping cpg first would stop at 120
        q2 = dict(cls=0, cpg=1, gcq=2, afq=2, bpref=0)
        idx2 = self._index([((0, 1, 2, 2, 0), 20), ((0, 1, 2, 4, 0), 20), ((0, 1, 3, 1, 0), 5),
                            ((0, 0, 2, 2, 0), 100), ((0, 1, 2, 2, 1), 999)])
        lev, mem, _, sizes, _ = idx2.lookup(q2, 50)
        self.assertEqual(mn.K1_LEVELS[lev][0], "drop_af_gc_cpg")
        self.assertEqual(sizes, [20, 40, 45, 145])
        # bp_is_ref and class are never dropped: the 999 bp_is_ref = 1 members are never candidates
        self.assertTrue((mem < 145).all())
        # enough members at the full key -> no fallback
        lev, mem, _, sizes, _ = self._index([((1, 0, 0, 0, 1), 60)]).lookup(dict(cls=1, cpg=0, gcq=0, afq=0, bpref=1), 50)
        self.assertEqual((lev, sizes), (0, [60]))
        # nothing reaches the minimum -> coarsest level, flagged
        lev, mem, _, sizes, below = self._index([((1, 0, 0, 0, 1), 5), ((1, 1, 3, 3, 1), 5)]).lookup(
            dict(cls=1, cpg=0, gcq=0, afq=0, bpref=1), 50)
        self.assertEqual((lev, len(mem), below), (3, 10, True))

    def test_fallbacks_logged(self):
        tmp = Path(tempfile.mkdtemp(prefix="tfbs_fb_"))
        try:
            paths, *_ = sf.make_fixture(tmp / "fx", seed=9, n_pool1=3000, n_chain=300)
            res = run_fixture(paths, tmp / "out")
            fb = pd.read_csv(tmp / "out" / "fallbacks.tsv", sep="\t")
            n_peak = int(res["units"]["is_peak"].sum())
            self.assertEqual(len(fb), n_peak)                       # one row per peak unit
            self.assertTrue(set(fb.level) <= {x[0] for x in mn.K1_LEVELS + mn.K2_LEVELS})
            self.assertTrue((fb.loc[fb.level != "full", "fallback"]).all())
            capped = fb[fb.k > 4]
            self.assertTrue(len(capped) >= 2 and capped.k_capped.all() and (capped.pool_k == 4).all())
            k1 = fb[fb.stratum_type == "k1"]
            self.assertTrue(k1.stratum_key.str.contains("bp_is_ref=").all())
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_sample_without_replacement(self):
        rng = np.random.default_rng(3)
        idx, repl = mn.sample_without_replacement(rng, 30, 12, 4000)
        self.assertFalse(repl)
        self.assertTrue(all(len(set(r)) == 12 for r in idx))
        counts = np.bincount(idx.ravel(), minlength=30) / idx.size
        np.testing.assert_allclose(counts, 1 / 30, atol=0.006)      # uniform over members
        pos = np.bincount(idx[:, 0], minlength=30) / len(idx)       # every slot is uniform too
        np.testing.assert_allclose(pos, 1 / 30, atol=0.012)
        _, repl = mn.sample_without_replacement(rng, 3, 5, 10)
        self.assertTrue(repl)


# ---------------------------------------------------------------------------------------------------------------------
# Directional calls on a hand-made unit (standard vs strong-change rule; orientation from bp_mask, not REF)
# ---------------------------------------------------------------------------------------------------------------------
class TestDirectionalCalls(TmpDirCase):
    def test_hand_unit(self):
        rng = np.random.default_rng(0)
        df = sf.build_unit_records(rng, 60, 0, 1, {"sptbn4b": [1]}, {})
        df.loc[0, "bp_mask"] = 1                                   # BP allele is ALT
        motifs = sf.motif_table()
        sf.write_units(df, self.tmp / "units.tsv")
        # peak unit 0: motif 0 BP 5.0 / TP 3.5 -> BP-only, not strong; motif 1 TP 4.5 only -> TP-only, strong;
        # motif 2 on both -> none; motif 3 BP 3.5 only -> none (below 1e-4)
        h = dict(unit_idx=np.array([0, 0, 0, 0, 0, 0]), mask=np.array([1, 0, 0, 1, 0, 1]),
                 motif_idx=np.array([0, 0, 1, 2, 2, 3]), mlog10p=np.array([5.0, 3.5, 4.5, 5.0, 5.0, 3.5]))
        sf.write_hits(h, motifs, self.tmp / "hits.npz")
        units = mn.load_units(self.tmp / "units.tsv")
        hits = mn.load_hits(self.tmp / "hits.npz", units)
        rows = mn.build_rows(units)
        hap = mn.hits_to_rows(hits, units, rows)
        motif_fam = np.arange(len(hits["motif_ids"])) % 5
        T = mn.directional_tables(hap, hits["motif_idx"], hits["mlog10p"], rows, motif_fam, 5,
                                  np.full(len(hits["motif_ids"]), -1), 2, len(hits["motif_ids"]))
        r = rows["hap_off"][0]                                     # BP row (mask = bp_mask = 1)
        self.assertEqual(int(rows["row_mask"][r]), 1)
        self.assertEqual((T["standard"]["nBP"][r], T["standard"]["nTP"][r]), (1, 1))
        self.assertEqual((T["strong"]["nBP"][r], T["strong"]["nTP"][r]), (0, 1))
        # reading the same unit with BP = REF flips the calls
        self.assertEqual((T["standard"]["nBP"][r + 1], T["standard"]["nTP"][r + 1]), (1, 1))
        self.assertEqual((T["strong"]["nBP"][r + 1], T["strong"]["nTP"][r + 1]), (1, 0))
        B = T["standard"]["fam"]["B"]
        self.assertEqual(B[r].toarray().ravel().tolist(), [1, 0, 0, 0, 0])
        self.assertEqual(T["standard"]["fam"]["Bt"][r].toarray().ravel().tolist(), [0, 1, 0, 0, 0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
