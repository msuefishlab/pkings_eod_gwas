# 16_IHC: SPTBN4 immunofluorescence quantification

Produces: the quantification of SPTBN4 labelling in Apassa (TP) vs Bambomo (BP) electrocytes (Fig. 4B).

## Inputs
- Leica `.lif` files in `$sptbn4_ihc_dir` (set in `config.env`): `PKINGS_SPTBN_IHC.lif` (2025-10-14),
  `PKINGS_SPTBN_IHC_OCT23.lif` (2025-10-23) and `PKINGS_SPTBN_IHC_NOV11_REP_SECTS_CHECKED.lif` (2025-11-11). One APA (TP)
  and one BAM (BP) fish; APA and BAM slides were stained together in each session. The raw images are not distributed
  with this repository.

Channels: L5 cube = anti-SPTBN4 with Alexa 488, DAP = DAPI, Ana = DIC; all 8-bit.

## Steps
Run from the repository root after `source config.env` (python3 with numpy, scipy and Pillow).
1. `python3 code/16_IHC/05_parse_lif.py "$sptbn4_ihc_dir" output_data/16_IHC/sptbn4_acquisition_settings.csv` — reads the
   `.lif` files without readlif and tabulates per-field L5/DAPI exposure, lamp and gain, with a `comparable_group` key
   (same session, objective, image size and L5 settings).
2. `python3 code/16_IHC/06_quantify_sptbn4_ihc.py "$sptbn4_ihc_dir" output_data/16_IHC/sptbn4_acquisition_settings.csv output_data/16_IHC`
   — band signal (L5 in the DAPI-defined electrocyte band minus the matched no-primary control) per field, BAM/APA
   ratios per comparable group with bootstrap intervals, sensitivity columns and QC panels. `ihc_masks.py` holds the
   band and blob masks.

Notes: four comparable groups have no-primary controls for both populations (2025-10-14 40×, 2025-10-23 20×,
2025-11-11 40× and 20×); the 2025-10-14 100× fields have no control and are reported on raw intensity only. One fish
per population, so fields are technical replicates and there is no population test. The blob mask is not a vessel
mask: in APA it mostly catches stalks.

## Outputs
- `output_data/16_IHC/sptbn4_acquisition_settings.csv` — per-image settings and `comparable_group`
- `output_data/16_IHC/sptbn4_field_metrics.csv`, `sptbn4_group_summary.csv` — per-field metrics and per-group ratios
- `output_data/16_IHC/qc/sptbn4_*.png` — DAPI, L5 and mask overlays per field
