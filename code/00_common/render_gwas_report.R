#!/usr/bin/env Rscript


# Render a GWAS-style R Markdown report (params data_path, MAF_thresh, test).
# Usage: Rscript code/00_common/render_gwas_report.R DATA_PATH OUTPUT_PREFIX REPORT.Rmd MAF_THRESHOLD [TEST]
# Writes OUTPUT_PREFIX.report.html. TEST is lrt (default, GEMMA -lmm 2; the published run) or wald (-lmm 1).

args <- commandArgs(trailingOnly = TRUE)
test <- if (length(args) >= 5) args[5] else "lrt"

root <- rprojroot::find_root(".git/index")
require(yaml)
require(tools)


rmarkdown::render(
    input = file.path(args[3]),
    output_file = basename(paste0(args[2], ".report.html")),
    output_dir = dirname(args[2]),
    params = list(
        data_path = args[1],
        MAF_thresh = as.numeric(args[4]),
        test = test
    )
)
