# utils.R
# Helper functions for the reproducible pipeline

library(yaml)
library(withr)

load_config <- function(file = "config.yml") {
  if (!file.exists(file)) stop("config.yml not found. Please create it.")
  cfg <- yaml::read_yaml(file)
  return(cfg)
}

ensure_dirs <- function(cfg){
  dirs <- c(cfg$paths$processed_dir, cfg$paths$results_dir, cfg$paths$logs_dir, "data/sample")
  for (d in dirs){
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
}

safe_save_rds <- function(obj, path){
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, file = path)
  message("Saved: ", path)
}

log_msg <- function(...){
  message(Sys.time(), " - ", paste(...))
}

# simple function to read config and ensure dirs
init_pipeline <- function(){
  cfg <- load_config()
  ensure_dirs(cfg)
  invisible(cfg)
}

