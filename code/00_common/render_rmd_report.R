#!/usr/bin/env Rscript


# Render an R Markdown report that declares the params data_path and MAF_thresh.
# Usage: Rscript code/00_common/render_rmd_report.R DATA_PATH OUTPUT_PREFIX REPORT.Rmd [MAF_THRESHOLD]
# Writes OUTPUT_PREFIX.report.html.

args <- commandArgs(trailingOnly = TRUE)

root <- rprojroot::find_root(".git/index")
require(yaml)
require(tools)


rmarkdown::render(
  input = file.path(args[3]),
  output_file = basename(paste0(args[2], ".report.html")),
  output_dir = dirname(args[2]),
  params = list(
    data_path = args[1],
    MAF_thresh = args[4]
  )
)
