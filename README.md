# The Edge of Adaptation

This repository documents a clear step-by-step workflow to run the main analyses for Ameiva ameiva:
- initial filtering (same steps used in the Lampropholis project),
- dartR-based preprocessing,
- RDA (partial RDA) and outlier detection,
- projection of adaptive indices and genomic offset.

The README below is a single file you can copy to `README.md`. It contains runnable R code snippets (filtering, preprocess, RDA, offset) so you can paste them into scripts or run interactively. Edit `config.yml` first to point to your files.

---

## Table of contents

1. Requirements
2. Repository layout
3. Configure paths and thresholds (config.yml example)
4. Step 1 — Filtering (STACKS)
5. Step 2 — Preprocessing wrapper (dartR)
6. Step 3 — RDA, variable selection and outlier detection
7. Step 4 — Genomic offset (projection)
8. Run commands (sequence)
9. Outputs to check
10. Notes & troubleshooting
11. License

---

## 1. Requirements

System
- Linux or macOS (some commands assume POSIX shell)
- R >= 4.0 (4.5.3 recommended if matching used packages)
- System libraries as needed (libxml2-dev, libssl-dev, libcurl4-openssl-dev, gdal-bin, libgdal-dev)

R packages (install as needed):
- dartR, vcfR, vegan, Boruta, raster, sf, qvalue, robust, ggplot2, yaml, optparse, withr

Example to install common packages:
```r
install.packages(c("vcfR","raster","sf","yaml","optparse","withr","ggplot2"))
# install dartR and other packages from CRAN or GitHub as appropriate
```

## 2. Repository layout (recommended)
config.yml — configuration of input paths and thresholds
code/01_preprocess_dartR.R — preprocessing script (dartR)
code/02_rda.R — RDA & outlier detection script
code/03_offset.R — genomic offset projection script
data/raw/ — raw files (VCF, metadata CSV, rasters, shapefiles) — do not commit large files
data/processed/ — processed R objects (RDS)
results/ — results outputs (RDA summary, outliers, offsets, figures)


## 3. Step 1 — Filtering
This exact sequence is applied prior to downstream analyses. Run interactively or save as code/00_filtering_example.R.

Locus call rate (loc)
Individual call rate (ind)
Remove monomorphic loci
Remove loci with all NA
MAF >= 0.05
LD pruning: keep 1 SNP per 1000 bp per contig (if mapping info available, For Ameiva we use a draft genome for Reference)
```R
# code/00_filtering_example.R
library(vcfR)
library(dartR)

vcf_path <- "data/raw/populations.snps.vcf"
meta_path <- "data/raw/Ameiva_Metadata.csv"  # optional

vcf <- read.vcfR(vcf_path)
gl <- vcfR2genlight(vcf)

# optional metadata
if (file.exists(meta_path)) gl <- gl.add.indmetrics(gl, meta_path)

# 1) locus call rate
gl <- gl.filter.callrate(gl, method = "loc", threshold = 0.7, verbose = 0)

# 2) individual call rate
gl <- gl.filter.callrate(gl, method = "ind", threshold = 0.7, verbose = 0)

# 3) remove monomorphic loci
gl <- gl.filter.monomorphs(gl)

# 4) remove loci where all genotypes are NA
gl <- gl.filter.allna(gl)

# 5) MAF filter (custom)
geno_mat <- as.matrix(gl)
maf_calc <- function(col) {
  p <- mean(col, na.rm = TRUE) / 2
  min(p, 1 - p)
}
mafs <- apply(geno_mat, 2, maf_calc)
keep_loci <- names(mafs)[mafs >= 0.05]
gl <- gl.keep.loc(gl, loc.list = keep_loci)

# 6) LD pruning (1 SNP per 1000 bp per contig)
if (!is.null(gl@position) && !is.null(gl@chromosome)) {
  loc_info <- data.frame(loc = locNames(gl), chr = as.character(gl@chromosome), pos = gl@position, stringsAsFactors = FALSE)
  loc_info <- loc_info[order(loc_info$chr, loc_info$pos), ]
  keep <- character(0)
  last_chr <- ""
  last_pos <- -Inf
  for (i in seq_len(nrow(loc_info))) {
    chr_i <- loc_info$chr[i]
    pos_i <- as.numeric(loc_info$pos[i])
    if (chr_i != last_chr || (pos_i - last_pos) >= 1000) {
      keep <- c(keep, loc_info$loc[i])
      last_chr <- chr_i
      last_pos <- pos_i
    }
  }
  gl <- gl.keep.loc(gl, loc.list = keep)
} else {
  message("No map info; skipping LD window pruning.")
}

# Save results
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(gl, "data/processed/genlight_filtered.rds")
pcoa <- gl.pcoa(gl)
saveRDS(pcoa, "data/processed/pcoa.rds")
AllFreqi <- as.matrix(gl)
saveRDS(AllFreqi, "data/processed/AllFreqi.rds")

if (!is.null(pop(gl))) {
  AllFreqp <- aggregate(AllFreqi, by = list(pop(gl)), FUN = function(x) mean(x, na.rm = TRUE)/2)
  row.names(AllFreqp) <- as.character(AllFreqp$Group.1)
  saveRDS(AllFreqp, "data/processed/AllFreqp.rds")
}
```
## 5. Step 2 — Preprocessing wrapper (dartR)
Place as code/01_preprocess_dartR.R. This loads config.yml, runs the filters above and writes processed RDS files.

```R
# code/01_preprocess_dartR.R
library(yaml)
library(dartR)
library(vcfR)
library(dplyr)

cfg <- yaml::read_yaml("config.yml")

run_preprocess <- function(cfg) {
  vcf_path <- cfg$data$raw_vcf
  meta_path <- cfg$data$metadata
  processed_dir <- cfg$paths$processed_dir
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

  processed_rds <- file.path(processed_dir, "genlight_filtered.rds")
  if (file.exists(processed_rds)) {
    message("Loading existing processed genlight")
    ref_gl <- readRDS(processed_rds)
  } else {
    message("Reading VCF")
    vcf <- read.vcfR(vcf_path)
    ref_gl <- vcfR2genlight(vcf)
    if (file.exists(meta_path)) ref_gl <- gl.add.indmetrics(ref_gl, meta_path)

    # filtering
    ref_gl <- gl.filter.callrate(ref_gl, method = "loc", threshold = cfg$filters$callrate_loc, verbose = 0)
    ref_gl <- gl.filter.callrate(ref_gl, method = "ind", threshold = cfg$filters$callrate_ind, verbose = 0)
    ref_gl <- gl.filter.monomorphs(ref_gl)
    ref_gl <- gl.filter.allna(ref_gl)

    geno_mat <- as.matrix(ref_gl)
    maf_calc <- function(col) { p <- mean(col, na.rm = TRUE)/2; min(p, 1-p) }
    mafs <- apply(geno_mat, 2, maf_calc)
    keep_loci <- names(mafs)[which(mafs >= cfg$filters$maf)]
    ref_gl <- gl.keep.loc(ref_gl, loc.list = keep_loci)

    if (!is.null(ref_gl@position) && !is.null(ref_gl@chromosome)) {
      loc_info <- data.frame(loc = locNames(ref_gl), chr = as.character(ref_gl@chromosome), pos = ref_gl@position, stringsAsFactors = FALSE)
      loc_info <- loc_info[order(loc_info$chr, loc_info$pos), ]
      keep <- character(0); last_chr <- ""; last_pos <- -Inf
      for (i in seq_len(nrow(loc_info))) {
        chr_i <- loc_info$chr[i]; pos_i <- as.numeric(loc_info$pos[i])
        if (chr_i != last_chr || (pos_i - last_pos) >= cfg$filters$ld_window_bp) {
          keep <- c(keep, loc_info$loc[i]); last_chr <- chr_i; last_pos <- pos_i
        }
      }
      ref_gl <- gl.keep.loc(ref_gl, loc.list = keep)
    } else {
      message("No map info; skip LD window.")
    }

    saveRDS(ref_gl, processed_rds)
  }

  # PCoA / tables
  pcoa <- gl.pcoa(ref_gl); saveRDS(pcoa, file.path(processed_dir,"pcoa.rds"))
  AllFreqi <- as.matrix(ref_gl); saveRDS(AllFreqi, file.path(processed_dir,"AllFreqi.rds"))
  if (!is.null(pop(ref_gl))) {
    AllFreqp <- aggregate(AllFreqi, by = list(pop(ref_gl)), FUN = function(x) mean(x, na.rm = TRUE)/2)
    row.names(AllFreqp) <- as.character(AllFreqp$Group.1)
    saveRDS(AllFreqp, file.path(processed_dir,"AllFreqp.rds"))
  }

  message("Preprocess finished")
  return(invisible(TRUE))
}

if (interactive() == FALSE) run_preprocess(cfg)
```

## 6. Step 3 — RDA, variable selection and outlier detection
Place as code/02_rda.R. This runs Boruta variable selection, builds a partial RDA (Condition on neutral structure if available), and runs the rdadapt-like outlier scan on RDA loadings. This script was build based on Capblanq RDA Swiss Knife original paper and also with Josue Azevedo mentoring.

```R
# code/02_rda.R
library(yaml)
library(vegan)
library(Boruta)
library(qvalue)
library(robust)
library(raster)

cfg <- yaml::read_yaml("config.yml")
processed_dir <- cfg$paths$processed_dir
results_dir <- cfg$paths$results_dir
dir.create(file.path(results_dir,"rda"), recursive = TRUE, showWarnings = FALSE)

AllFreq.i <- readRDS(file.path(processed_dir,"AllFreqi.rds"))
if (file.exists(file.path(processed_dir,"Variables.f.rds"))) {
  Variables.f <- readRDS(file.path(processed_dir,"Variables.f.rds"))
} else {
  Variables.f <- read.csv("data/sample/variables_sample.csv", row.names = 1)
}

# Boruta variable selection (example)
numeric_preds <- Variables.f[, sapply(Variables.f, is.numeric), drop = FALSE]
Y <- rowMeans(AllFreq.i, na.rm = TRUE)
df_boruta <- cbind(Y = Y, numeric_preds)
set.seed(7)
boruta_res <- tryCatch(Boruta::Boruta(Y ~ ., data = df_boruta, doTrace = 0, pValue = 0.05, maxRuns = 100),
                       error = function(e) NULL)
if (!is.null(boruta_res)) {
  boruta_final <- TentativeRoughFix(boruta_res)
  selected_vars <- getSelectedAttributes(boruta_final, withTentative = FALSE)
} else {
  selected_vars <- colnames(numeric_preds)[1:min(6,ncol(numeric_preds))]
}
message("Selected variables: ", paste(selected_vars, collapse = ", "))

# Build partial RDA formula
cond_vars <- NULL
if (all(c("PC1","PC2","PC3") %in% colnames(Variables.f))) cond_vars <- c("PC1","PC2","PC3")
if (is.null(cond_vars) && all(c("clust1","clust2","clust3") %in% colnames(Variables.f))) cond_vars <- c("clust1","clust2","clust3")

formula_str <- paste0("AllFreq.i ~ ", paste(selected_vars, collapse = " + "))
if (!is.null(cond_vars)) formula_str <- paste0(formula_str, " + Condition(", paste(cond_vars, collapse = " + "), ")")
RDAenv <- rda(as.formula(formula_str), Variables.f)
saveRDS(RDAenv, file.path(results_dir,"rda","RDA_env.rds"))

# Summarize
varex <- summary(RDAenv)$cont$importance
rda1_pct <- round(varex[2,1]*100,2); rda2_pct <- round(varex[2,2]*100,2)
model_p <- tryCatch(anova(RDAenv, permutations = cfg$rda$permutation)$"Pr(>F)"[1], error = function(e) NA)
yaml::write_yaml(list(rda1_pct=rda1_pct, rda2_pct=rda2_pct, model_p=model_p), file.path(results_dir,"rda","summary.yml"))

# rdadapt-like outlier detection
rdadapt <- function(rda_obj, K) {
  zscores <- rda_obj$CCA$v[,1:K]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- median(resmaha) / qchisq(0.5, df = K)
  pvals <- pchisq(resmaha / lambda, df = K, lower.tail = FALSE)
  qvals <- qvalue::qvalue(pvals)$qvalues
  data.frame(p.values = pvals, q.values = qvals)
}

rd_res <- rdadapt(RDAenv, cfg$rda$rda_axes)
thres <- 0.01 / nrow(rd_res)
outliers_idx <- which(rd_res$p.values < thres)
outlier_names <- if (length(outliers_idx)>0) colnames(AllFreq.i)[outliers_idx] else character(0)
saveRDS(outlier_names, file.path(results_dir,"rda","outliers_list.rds"))
message("Outliers:", length(outlier_names))
```
## 7. Step 4 — Genomic offset (projection)
Place as code/03_offset.R. This uses RDA loadings to compute per-pixel adaptive index for present and future and the weighted euclidean distance (global genomic offset).
```R
# code/03_offset.R
library(raster)
library(yaml)

cfg <- yaml::read_yaml("config.yml")
results_dir <- cfg$paths$results_dir
processed_dir <- cfg$paths$processed_dir

RDAenv <- readRDS(file.path(results_dir,"rda","RDA_env.rds"))
center_env <- if (file.exists(file.path(processed_dir,"center_env.rds"))) readRDS(file.path(processed_dir,"center_env.rds")) else NULL
scale_env  <- if (file.exists(file.path(processed_dir,"scale_env.rds")))  readRDS(file.path(processed_dir,"scale_env.rds"))  else NULL
if (is.null(center_env) || is.null(scale_env)) stop("center_env/scale_env required for projection. Save them during env preparation.")

# helper to project
adaptive_index <- function(RDA, K, env_stack, center_env, scale_env, mask = NULL) {
  vars <- row.names(RDA$CCA$biplot)
  stopifnot(all(vars %in% names(env_stack)))
  if (!is.null(mask)) env_stack <- mask(env_stack, mask)
  pts <- rasterToPoints(env_stack[[1]])[,1:2]
  vals <- as.data.frame(rasterToPoints(env_stack)[, -(1:2)])
  colnames(vals) <- names(env_stack)
  vals_std <- scale(vals, center = center_env[vars], scale = scale_env[vars])
  Proj_list <- list()
  for (i in 1:K) {
    load_i <- RDA$CCA$biplot[,i]
    scores_vec <- as.vector(as.matrix(vals_std[,vars]) %*% as.numeric(load_i))
    Proj_list[[i]] <- rasterFromXYZ(data.frame(pts, Z = scores_vec), crs = crs(env_stack))
  }
  Proj_list
}

genomic_offset <- function(RDA, K, env_pres, env_fut, center_env, scale_env, mask = NULL) {
  pres_list <- adaptive_index(RDA, K, env_pres, center_env, scale_env, mask)
  fut_list  <- adaptive_index(RDA, K, env_fut,  center_env, scale_env, mask)
  eig <- RDA$CCA$eig; weights <- eig[1:K] / sum(eig[1:K])
  pts <- rasterToPoints(pres_list[[1]])[,1:2]
  pres_mat <- do.call(cbind, lapply(pres_list, raster::values))
  fut_mat  <- do.call(cbind, lapply(fut_list,  raster::values))
  pres_w <- sweep(pres_mat, 2, weights, "*")
  fut_w  <- sweep(fut_mat,  2, weights, "*")
  diffs <- sqrt(rowSums((pres_w - fut_w)^2, na.rm = TRUE))
  rasterFromXYZ(data.frame(pts, Z = diffs), crs = crs(env_pres))
}

# load rasters (present and future)
ras_pres_files <- list.files(cfg$data$rasters_dir, pattern = cfg$climate$present_pattern, full.names = TRUE)
ras_fut_files  <- list.files(cfg$data$rasters_dir, pattern = cfg$climate$fut2100_pattern, full.names = TRUE)
ras_pres <- stack(ras_pres_files); ras_fut <- stack(ras_fut_files)

K <- cfg$rda$rda_axes
offset_ras <- genomic_offset(RDAenv, K, ras_pres, ras_fut, center_env, scale_env)
dir.create(file.path(results_dir,"offset"), recursive = TRUE, showWarnings = FALSE)
writeRaster(offset_ras, file.path(results_dir,"offset","genomic-offset-2100.tif"), overwrite = TRUE)
message("Saved: results/offset/genomic-offset-2100.tif")
```
## Notes & troubleshooting
If Boruta or other packages fail on large data, run variable selection on a subset of predictors or increase memory.
If rasters do not align or have NA areas, use the helper to mask rasters with species range shapefile before projection.
Save center_env and scale_env used to standardize predictors during model fitting — those must be used for future projections.
For large raster projections, aggregate or run on HPC.

