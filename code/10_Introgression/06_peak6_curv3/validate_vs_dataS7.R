## validate_vs_dataS7.R
## Compares the rerun of the main analysis's Peak 6 Dtrios (01_validate_peak6.sh) with the Peak 6 rows of
## Supplementary Data S7 (download it to input_data/10_Introgression/SupplementalDataS7_peak_ABBA_BABA.tsv).
## Trio choice and signing are copied from code/10_Introgression/DsuiteAnalysis.Rmd (hyp_defs, orient_trios).
## D_genome is the main analysis's whole_genome_combined_tree.txt; if 03_run_dtrios.sh has also rerun the main sets
## genome-wide (genome_published_sets/whole_genome_combined_tree.txt), that value is reported alongside.
## q_rerun: Data S7 q is BH over all 48 peak tests, so the rerun's 8 Peak 6 p values replace Data S7's (rounded) Peak 6
## p values and the 48 are adjusted again.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr); library(tibble)
})

root <- system("git rev-parse --show-toplevel", intern = TRUE)
intro_dir <- file.path(root, "output_data", "10_Introgression")
out_dir   <- file.path(intro_dir, "peak6_curv3")

read_dsuite <- function(f, peak_label = NA_character_) {
  df <- read_tsv(f, show_col_types = FALSE, progress = FALSE)
  names(df)[1:10] <- c("P1","P2","P3","D","Z","p","f4ratio","BBAA","ABBA","BABA")
  df %>% mutate(Peak = peak_label)
}

hyp_defs <- tribble(
  ~Hypothesis,     ~focal, ~sister,  ~P3,      ~type,
  "BP1–BP2",       "BP1",  "TP1",    "BP2",    "BP–BP",
  "BP1–BP3",       "BP1",  "TP1",    "BP3",    "BP–BP",
  "BP2–BP3",       "BP2",  "TP2",    "BP3",    "BP–BP",
  "BP1–BP_ANC",    "BP1",  "TP_OUT", "BP_ANC", "BP–other biphasic",
  "BP2–BP_ANC",    "BP2",  "TP_OUT", "BP_ANC", "BP–other biphasic",
  "BP3–BP_ANC",    "BP3",  "TP_OUT", "BP_ANC", "BP–other biphasic",
  "TP1–BP_ANC",    "TP1",  "TP_OUT", "BP_ANC", "TP control",
  "TP2–BP_ANC",    "TP2",  "TP_OUT", "BP_ANC", "TP control"
)

orient_trios <- function(df) {
  hyp_defs %>%
    pmap_dfr(function(Hypothesis, focal, sister, P3, type) {
      hit <- df[df$P3 == P3 & ((df$P1 == focal & df$P2 == sister) | (df$P1 == sister & df$P2 == focal)), ]
      stopifnot(nrow(hit) == n_distinct(df$Peak))
      hit %>% mutate(Hypothesis = Hypothesis, type = type, focal = focal, sister = sister,
                     sign = ifelse(P2 == focal, 1, -1))
    }) %>%
    mutate(D_o = sign * D, Z_o = sign * Z, f4_o = sign * f4ratio,
           SE = ifelse(Z != 0, abs(D / Z), NA_real_),
           nSites = ABBA + BABA)
}

gw <- read_dsuite(file.path(intro_dir, "whole_genome_combined_tree.txt"), "genome") %>%
  orient_trios() %>% select(Hypothesis, D_gw = D_o)

rerun <- read_dsuite(file.path(out_dir, "validation_published_sets", "peak_6_tree.txt"), "Peak 6") %>%
  orient_trios() %>%
  left_join(gw, by = "Hypothesis") %>%
  mutate(Z_diff = (D_o - D_gw) / SE)

s7 <- read_tsv(file.path(root, "input_data/10_Introgression/SupplementalDataS7_peak_ABBA_BABA.tsv"),
               show_col_types = FALSE)

# q over the 48 peak tests, with the rerun's Peak 6 p values swapped in
p48 <- s7 %>% select(Hypothesis, Peak, p) %>%
  left_join(rerun %>% transmute(Hypothesis, Peak, p_rerun = p), by = c("Hypothesis", "Peak")) %>%
  mutate(p_use = ifelse(is.na(p_rerun), p, p_rerun), q_rerun = p.adjust(p_use, method = "BH")) %>%
  filter(Peak == "Peak 6") %>% select(Hypothesis, q_rerun)

s7_6 <- s7 %>% filter(Peak == "Peak 6")

val <- s7_6 %>%
  select(Hypothesis, P1, P2, P3, D_S7 = D, Z_S7 = Z, p_S7 = p, q_S7 = q, f4_S7 = f4_ratio,
         D_genome_S7 = D_genome, Z_vs_genome_S7 = Z_vs_genome, sites_S7 = informative_sites) %>%
  left_join(rerun %>% transmute(Hypothesis, P1_rerun = P1, P2_rerun = P2, D_rerun = round(D_o, 3),
                                Z_rerun = round(Z_o, 2), p_rerun = signif(p, 2), f4_rerun = round(f4_o, 3),
                                Z_vs_genome_rerun = round(Z_diff, 2), sites_rerun = round(nSites, 1)),
            by = "Hypothesis") %>%
  left_join(p48 %>% mutate(q_rerun = signif(q_rerun, 2)), by = "Hypothesis")

# Optional: genome-wide rerun of the main sets
gw_local_f <- file.path(out_dir, "genome_published_sets", "whole_genome_combined_tree.txt")
if (file.exists(gw_local_f)) {
  gw_local <- read_dsuite(gw_local_f, "genome") %>% orient_trios() %>%
    transmute(Hypothesis, D_genome_rerun = round(D_o, 3), Z_genome_rerun = round(Z_o, 2))
  gw_s7 <- read_dsuite(file.path(intro_dir, "whole_genome_combined_tree.txt"), "genome") %>% orient_trios() %>%
    transmute(Hypothesis, Z_genome_published = round(Z_o, 2))
  val <- val %>% left_join(gw_local, by = "Hypothesis") %>% left_join(gw_s7, by = "Hypothesis") %>%
    left_join(rerun %>% transmute(Hypothesis, D_o, SE), by = "Hypothesis") %>%
    mutate(Z_vs_genome_rerun_localgenome = round((D_o - D_genome_rerun) / SE, 2)) %>%
    select(-D_o, -SE)
}

val <- val %>%
  mutate(dD = round(D_rerun - D_S7, 3), dZ = round(Z_rerun - Z_S7, 2), dSites = round(sites_rerun - sites_S7, 1),
         same_orientation = (P1_rerun == P1 & P2_rerun == P2))

write_tsv(val, file.path(out_dir, "validation_vs_dataS7.tsv"))
print(as.data.frame(val %>% select(Hypothesis, D_S7, D_rerun, Z_S7, Z_rerun, sites_S7, sites_rerun, f4_S7, f4_rerun,
                                   Z_vs_genome_S7, Z_vs_genome_rerun, q_S7, q_rerun, same_orientation)))

# Raw check: every trio of the rerun against the main analysis's peaks_extended/peak_6_tree.txt (all 35 trios).
# The f4-ratio splits P3 at random (Dsuite seeds std::mt19937 from std::random_device), so it is not reproducible
# to the last digit; D, Z and the pattern counts are deterministic.
raw_pub <- read_dsuite(file.path(intro_dir, "peaks_extended", "peak_6_tree.txt"), "Peak 6")
raw_new  <- read_dsuite(file.path(out_dir, "validation_published_sets", "peak_6_tree.txt"), "Peak 6")
raw_cmp <- inner_join(raw_pub, raw_new, by = c("P1", "P2", "P3"), suffix = c("_published", "_rerun")) %>%
  transmute(P1, P2, P3, D_published, D_rerun, Z_published, Z_rerun, f4_published = f4ratio_published, f4_rerun = f4ratio_rerun,
            BBAA_published, BBAA_rerun, ABBA_published, ABBA_rerun, BABA_published, BABA_rerun,
            max_abs_diff_D_Z_counts = pmax(abs(D_published - D_rerun), abs(Z_published - Z_rerun), abs(BBAA_published - BBAA_rerun),
                                           abs(ABBA_published - ABBA_rerun), abs(BABA_published - BABA_rerun)),
            abs_diff_f4 = abs(f4_published - f4_rerun))
stopifnot(nrow(raw_cmp) == nrow(raw_pub), nrow(raw_cmp) == nrow(raw_new))
write_tsv(raw_cmp, file.path(out_dir, "validation_peak6_all_trios_vs_published.tsv"))
cat("All", nrow(raw_cmp), "trios: max |diff| in D, Z, BBAA, ABBA, BABA =", max(raw_cmp$max_abs_diff_D_Z_counts),
    "; max |diff| f4-ratio =", round(max(raw_cmp$abs_diff_f4), 4), "\n")

# Genome-wide baseline: the DtriosCombine of the rerun main sets (03_run_dtrios.sh) against the main analysis's
# whole_genome_combined_tree.txt, all 35 trios (written only once the genome-wide rerun exists)
if (file.exists(gw_local_f)) {
  g_old <- read_dsuite(file.path(intro_dir, "whole_genome_combined_tree.txt"), "genome")
  g_new <- read_dsuite(gw_local_f, "genome")
  g_cmp <- inner_join(g_old, g_new, by = c("P1", "P2", "P3"), suffix = c("_published", "_rerun")) %>%
    transmute(P1, P2, P3, D_published, D_rerun, Z_published, Z_rerun, BBAA_published, BBAA_rerun, ABBA_published, ABBA_rerun,
              BABA_published, BABA_rerun, f4_published = f4ratio_published, f4_rerun = f4ratio_rerun,
              rel_diff_counts = pmax(abs(BBAA_published / BBAA_rerun - 1), abs(ABBA_published / ABBA_rerun - 1), abs(BABA_published / BABA_rerun - 1)),
              abs_diff_D = abs(D_published - D_rerun), abs_diff_Z = abs(Z_published - Z_rerun))
  stopifnot(nrow(g_cmp) == nrow(g_old))
  write_tsv(g_cmp, file.path(out_dir, "validation_genome_published_sets_vs_published.tsv"))
  cat("Genome-wide, all", nrow(g_cmp), "trios: max |dD| =", signif(max(g_cmp$abs_diff_D), 3),
      "; max |dZ| =", signif(max(g_cmp$abs_diff_Z), 3), "; max relative diff in counts =", signif(max(g_cmp$rel_diff_counts), 3), "\n")
}
