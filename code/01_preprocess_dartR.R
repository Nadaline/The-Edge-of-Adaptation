# 01_preprocess_dartR.R
# Preprocessing adapted from Ameiva_dartRverse.R
# This script defines run_preprocess() which performs filtering and saves processed objects.

run_preprocess <- function(){
  cfg <- init_pipeline()
  log_msg("Loading packages")
  suppressPackageStartupMessages({
    library(dartR)
    library(vcfR)
    library(dplyr)
  })

  # Paths
  vcf_path <- cfg$data$raw_vcf
  meta_path <- cfg$data$metadata
  processed_dir <- cfg$paths$processed_dir

  # If user provided processed RDS files, prefer those
  if (file.exists(file.path(processed_dir, "genlight_filtered.rds"))){
    log_msg("Processed genlight found. Loading existing RDS.")
    ref_gl <- readRDS(file.path(processed_dir, "genlight_filtered.rds"))
  } else {
    # Check for raw VCF
    if (!file.exists(vcf_path)){
      log_msg("VCF not found at ", vcf_path, ". For a smoke test, run data/sample/generate_sample.R first.")
      return(invisible(NULL))
    }
    log_msg("Reading VCF: ", vcf_path)
    vcf <- read.vcfR(vcf_path)
    ref_gl <- vcfR2genlight(vcf)

    # add metadata if available
    if (file.exists(meta_path)){
      log_msg("Adding metadata from: ", meta_path)
      ref_gl <- gl.add.indmetrics(ref_gl, meta_path)
    }

    # Filtering steps (Lampropholis style)
    log_msg("Filtering: locus call rate (loc)")
    ref_gl <- gl.filter.callrate(ref_gl, method = "loc", threshold = cfg$filters$callrate_loc, verbose = 0)
    log_msg("Filtering: individual call rate (ind)")
    ref_gl <- gl.filter.callrate(ref_gl, method = "ind", threshold = cfg$filters$callrate_ind, verbose = 0)
    log_msg("Removing monomorphic loci")
    ref_gl <- gl.filter.monomorphs(ref_gl)
    log_msg("Removing loci with all NA")
    ref_gl <- gl.filter.allna(ref_gl)

    # MAF filter
    log_msg("Applying MAF filter: ", cfg$filters$maf)
    # dartR does not have a direct gl.filter.maf in older versions; implement custom
    geno_mat <- as.matrix(ref_gl)
    maf_calc <- function(col){
      p <- mean(col, na.rm = TRUE)/2
      maf <- min(p, 1-p)
      return(maf)
    }
    mafs <- apply(geno_mat, 2, maf_calc)
    keep_loci <- names(mafs)[which(mafs >= cfg$filters$maf)]
    ref_gl <- gl.keep.loc(ref_gl, loc.list = keep_loci)

    # LD pruning: keep one SNP per 1000 bp per contig
    log_msg("LD pruning: keep one SNP per ", cfg$filters$ld_window_bp, " bp per contig")
    if (!is.null(ref_gl@position) & !is.null(ref_gl@chromosome)){
      loc_info <- data.frame(loc = locNames(ref_gl), chr = as.character(ref_gl@chromosome), pos = ref_gl@position, stringsAsFactors = FALSE)
      loc_info <- loc_info[order(loc_info$chr, loc_info$pos),]
      keep <- c()
      last_chr <- ""
      last_pos <- -Inf
      for (i in seq_len(nrow(loc_info))){
        chr_i <- loc_info$chr[i]
        pos_i <- as.numeric(loc_info$pos[i])
        if (chr_i != last_chr || (pos_i - last_pos) >= cfg$filters$ld_window_bp){
          keep <- c(keep, loc_info$loc[i])
          last_chr <- chr_i
          last_pos <- pos_i
        }
      }
      ref_gl <- gl.keep.loc(ref_gl, loc.list = keep)
    } else {
      log_msg("No chromosome/position information found; skipping LD window filter. You may supply mapped VCF." )
    }

    # Save processed genlight
    safe_save_rds(ref_gl, file.path(processed_dir, "genlight_filtered.rds"))
  }

  # Create basic summaries and PCoA
  log_msg("Computing PCoA and saving results")
  pcoa <- gl.pcoa(ref_gl)
  safe_save_rds(pcoa, file.path(processed_dir, "pcoa.rds"))

  # Export allele frequency tables used by RDA
  log_msg("Creating allele frequency tables (AllFreq.i and AllFreq.p)")
  geno_mat <- as.matrix(ref_gl)
  AllFreq.i <- geno_mat
  # make a simple population freq table if pop information exists
  if (!is.null(pop(ref_gl))){
    AllFreq.p <- aggregate(AllFreq.i, by = list(pop(ref_gl)), function(x) mean(x, na.rm = TRUE)/2)
    row.names(AllFreq.p) <- as.character(AllFreq.p$Group.1)
  } else {
    AllFreq.p <- NULL
  }
  safe_save_rds(AllFreq.i, file.path(processed_dir, "AllFreqi.rds"))
  if (!is.null(AllFreq.p)) safe_save_rds(AllFreq.p, file.path(processed_dir, "AllFreqp.rds"))

  log_msg("Preprocess finished")
  return(invisible(TRUE))
}
