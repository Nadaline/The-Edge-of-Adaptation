#!/usr/bin/env Rscript
# run_pipeline.R
# Simple wrapper to run pipeline steps: preprocess, rda, smoke_test, all

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-s", "--step"), type = "character", default = "all",
              help = "Step to run: preprocess, rda, smoke_test, all", metavar = "step")
)
opt <- parse_args(OptionParser(option_list = option_list))
step <- opt$step

message(Sys.time(), " - Starting pipeline step: ", step)

if (step %in% c("preprocess", "all", "smoke_test")) {
  source("code/utils.R")
  source("code/01_preprocess_dartR.R")
  run_preprocess()
}

if (step %in% c("rda", "all")) {
  source("code/utils.R")
  source("code/02_rda.R")
  run_rda()
}

if (step == "smoke_test") {
  message("smoke_test completed")
}

message(Sys.time(), " - Pipeline step finished: ", step)
