source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

target_state_formal <- "Clec7a+ lipid-handling Microglia"
target_state_legacy <- "Clec7a+ lipid-associated Microglia"

candidate_tfs <- c("Mef2a", "Nr3c1", "Fli1", "Zfp148", "Bclaf1")

functional_genes <- c("Clec7a", "Trem2", "Cst7", "Cx3cr1", "Axl", "C1qa", "Apoe", "Lpl", "Abca1", "Apoc1",
  "Apoc4", "Npc2", "Soat1", "Tyrobp", "Csf1r", "Itgam", "Hexb", "Ctsd", "Lamp1", "Lamp2")

qc_mtThreshold <- 0.1
qc_minLSize <- 1000L
qc_minCells <- 25L

max_network_genes <- 5000L

nc_lambda <- 0
nc_nNet <- 10L
nc_nCells <- 500L
nc_nComp <- 3L
nc_scaleScores <- TRUE
nc_symmetric <- FALSE
nc_q <- 0.9

td_K <- 3L
td_maxIter <- 1000L
td_maxError <- 1e-05
td_nDecimal <- 3L

ma_nDim <- 2L

n_cores <- min(8L, max(1L, parallel::detectCores() - 1L))

seed_WT <- 9101L
seed_MA_base <- 9200L

fdr_cutoff <- 0.05
top_n_export <- 50L
top_n_plot <- 10L

png_dpi <- 600

source_rds <- paste0(project_path("results/20_AUCell_lipid_metabolism_scoring/"), "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

frozen_signature_file <- paste0(project_path("results/21_Clec7a_DEG_volcano_myeloid/"), "21_Clec7a_vs_all_myeloid_DEG_up_v01.csv")

step90_functional_file <- paste0(project_path("results/90_pySCENIC_top10_TF_network_signature_correlation/"),
  "06_90_Figure6E_20gene_functional_support_audit_v05.csv")

step90_qc_file <- paste0(project_path("results/90_pySCENIC_top10_TF_network_signature_correlation/"), "14_90_FINAL_QC_v05.csv")

outdir <- project_path("results/91_scTenifoldKnk_top5_virtual_KO")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

output_index_file <- file.path(outdir, "00_91_OUTPUT_INDEX_v04.csv")

environment_file <- file.path(outdir, "01_91_environment_v04.csv")

input_qc_file <- file.path(outdir, "02_91_input_QC_v04.csv")

gene_selection_file <- file.path(outdir, "03_91_network_gene_selection_v04.csv.gz")

wt_network_file <- file.path(outdir, "04_91_WT_network_v04.rds")

combined_dr_file <- file.path(outdir, "05_91_top5_KO_diffRegulation_all_v04.csv.gz")

significant_dr_file <- file.path(outdir, "06_91_top5_KO_diffRegulation_FDR05_v04.csv.gz")

top50_file <- file.path(outdir, "07_91_top5_KO_top50_perturbed_genes_v04.csv")

ko_summary_file <- file.path(outdir, "08_91_top5_KO_summary_v04.csv")

figure6f_pdf <- file.path(outdir, "09_91_Figure6F_top5_virtual_KO_diagnostic_v04.pdf")

figure6f_svg <- file.path(outdir, "09_91_Figure6F_top5_virtual_KO_diagnostic_v04.svg")

figure6f_png <- file.path(outdir, "09_91_Figure6F_top5_virtual_KO_diagnostic_v04.png")

bundle_file <- file.path(outdir, "10_91_scTenifoldKnk_top5_virtual_KO_bundle_v04.rds")

final_qc_file <- file.path(outdir, "11_91_FINAL_QC_v04.csv")

report_file <- file.path(outdir, "12_91_REPORT_v04.txt")

checkpoint_dir <- file.path(outdir, "checkpoints_v04")

if (!dir.exists(checkpoint_dir)) {
  dir.create(checkpoint_dir, recursive = TRUE)
}

cp_counts_network <- file.path(checkpoint_dir, "CP01_counts_network_v04.rds")

cp_gene_selection <- file.path(checkpoint_dir, "CP01_gene_selection_v04.csv.gz")

cp_WT_ensemble <- file.path(checkpoint_dir, "CP02_WT_ensemble_v04.rds")

cp_WT_network <- file.path(checkpoint_dir, "CP03_WT_network_v04.rds")

timing_file <- file.path(checkpoint_dir, "CP_timing_history_v04.csv")

run_log_file <- file.path(checkpoint_dir, "CP_run_log_v04.txt")

make_qc <- qc_row

fmt_duration <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) {
    return("unknown")
  }
  h <- floor(seconds/3600)
  m <- floor((seconds%%3600)/60)
  s <- round(seconds%%60)
  if (h > 0) {
    sprintf("%dh %02dm %02ds", h, m, s)
  } else if (m > 0) {
    sprintf("%dm %02ds", m, s)
  } else {
    sprintf("%ds", s)
  }
}

fmt_clock <- function(x) {
  format(x, "%Y-%m-%d %H:%M:%S")
}

log_msg <- function(...) {
  msg <- paste0(...)
  line <- paste0("[", fmt_clock(Sys.time()), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = run_log_file, append = TRUE)
  invisible(line)
}

read_timing <- function() {
  if (!file.exists(timing_file)) {
    return(data.frame(stage = character(0), seconds = numeric(0), completed_at = character(0), stringsAsFactors = FALSE))
  }
  read.csv(timing_file, stringsAsFactors = FALSE, check.names = FALSE)
}

record_timing <- function(stage, seconds) {
  hist <- read_timing()
  new_row <- data.frame(stage = stage, seconds = as.numeric(seconds), completed_at = fmt_clock(Sys.time()),
    stringsAsFactors = FALSE)
  hist <- rbind(hist[hist$stage != stage, , drop = FALSE], new_row)
  write.csv(hist, timing_file, row.names = FALSE)
  invisible(hist)
}

show_eta_from_seconds <- function(label, seconds_remaining) {
  if (!is.finite(seconds_remaining) || seconds_remaining < 0) {
    return(invisible(NULL))
  }
  eta_time <- Sys.time() + seconds_remaining
  log_msg(label, ": remaining ~", fmt_duration(seconds_remaining), "; estimated finish ", fmt_clock(eta_time))
  invisible(eta_time)
}

checkpoint_exists <- function(path) {
  file.exists(path) && is.finite(file.info(path)$size) && file.info(path)$size > 0
}

safe_log2_fc <- function(x) {
    out <- rep(NA_real_, length(x))
    ok <- is.finite(x) & x > 0
    out[ok] <- log2(x[ok])
    out
}

rank_perturbation <- function(dr, ko_gene) {

  z <- dr

  if (!"gene" %in% colnames(z)) {
    z$gene <- rownames(z)
  }

  z$KO_TF <- ko_gene

  z$log2_perturbation_FC <- safe_log2_fc(z$FC)

  z <- z[order(-z$distance, z$p.adj, z$gene), , drop = FALSE]

  z$perturbation_rank <- seq_len(nrow(z))

  z$is_KO_gene <- z$gene == ko_gene

  z
}

log_msg("Step 91 v04 START/RESUME. Workstation tuning: ", n_cores, " workers on ", parallel::detectCores(),
  " detected logical cores.")

required_packages <- c("Matrix", "SeuratObject", "scTenifoldKnk", "scTenifoldNet", "ggplot2")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

knk_version <- as.character(packageVersion("scTenifoldKnk"))

net_version <- as.character(packageVersion("scTenifoldNet"))

if (knk_version != "1.0.3") {
  stop("This script is for scTenifoldKnk 1.0.3; observed ", knk_version, ".")
}

if (net_version != "1.3") {
  stop("This script is for scTenifoldNet 1.3; observed ", net_version, ".")
}

scQC_v103 <- getFromNamespace("scQC", "scTenifoldKnk")

scqc_formals <- names(formals(scQC_v103))

required_scqc_args <- c("X", "mtThreshold", "minLSize")

if (!all(required_scqc_args %in% scqc_formals)) {
  stop("Unexpected scQC() signature in scTenifoldKnk 1.0.3. Observed formals: ", paste(scqc_formals,
    collapse = ", "))
}

env <- data.frame(item = c("R", "scTenifoldKnk", "scTenifoldNet", "scQC_access", "n_cores", "WT_seed",
  "network_gene_cap"), value = c(R.version.string, knk_version, net_version, "getFromNamespace(scQC, scTenifoldKnk)",
  n_cores, seed_WT, max_network_genes), stringsAsFactors = FALSE)

write.csv(env, environment_file, row.names = FALSE)

required_files <- c(source_rds, frozen_signature_file, step90_functional_file, step90_qc_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
}

step90_qc <- read.csv(step90_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

if (any(step90_qc$status == "FAIL")) {
  stop("Step 90 v05 contains blocking FAIL item(s).")
}

sig <- read.csv(frozen_signature_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!"gene" %in% colnames(sig)) {
  stop("Frozen target-state signature lacks a 'gene' column.")
}

frozen_signature <- unique(as.character(sig$gene))

if (length(frozen_signature) != 1320L) {
  stop("Expected frozen 1,320-gene signature; observed ", length(frozen_signature), ".")
}

func_audit <- read.csv(step90_functional_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(functional_genes %in% func_audit$gene)) {
  stop("The Step 90 v05 functional panel does not match this script.")
}

if (checkpoint_exists(cp_counts_network) && checkpoint_exists(cp_gene_selection)) {

  log_msg("Checkpoint CP01 found: loading prepared network matrix; skipping source-object reload/QC/gene-selection.")

  counts_network <- readRDS(cp_counts_network)

  gene_selection <- read.csv(gzfile(cp_gene_selection), stringsAsFactors = FALSE, check.names = FALSE)

  if (!all(candidate_tfs %in% rownames(counts_network))) {
    stop("CP01 is invalid: one or more candidate TFs are missing.")
  }

  if (!all(intersect(functional_genes, gene_selection$gene) %in% rownames(counts_network))) {
    stop("CP01 is invalid: retained functional genes do not match checkpoint metadata.")
  }

  target_cells <- colnames(counts_network)

  n_target_cells_frozen <- 17713L
  hvg_capped <- nrow(counts_network) >= max_network_genes

} else {

  t_cp01 <- Sys.time()

  obj <- readRDS(source_rds)

  if (!"RNA" %in% names(obj@assays)) {
    stop("RNA assay not found.")
  }

  if (!"myeloid_state" %in% colnames(obj@meta.data)) {
    stop("myeloid_state metadata column not found.")
  }

  counts_all <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")

  cell_names <- colnames(counts_all)

  state <- as.character(obj@meta.data[cell_names, "myeloid_state"])

  state[state == target_state_legacy] <- target_state_formal

  target_cells_frozen <- cell_names[state == target_state_formal]

  if (length(target_cells_frozen) != 17713L) {
    stop("target-state count mismatch: ", length(target_cells_frozen), " / 17713.")
  }

  n_target_cells_frozen <- length(target_cells_frozen)

  counts_target <- counts_all[, target_cells_frozen, drop = FALSE]

  rm(counts_all, obj)

  gc()

  counts_cellqc <- scQC_v103(counts_target, mtThreshold = qc_mtThreshold, minLSize = qc_minLSize)

  if (ncol(counts_cellqc) < nc_nCells) {
    stop("Fewer than ", nc_nCells, " cells remain after QC.")
  }

  gene_detected_cells <- Matrix::rowSums(counts_cellqc != 0)

  keep_gene <- gene_detected_cells >= qc_minCells

  counts_genefilter <- counts_cellqc[keep_gene, , drop = FALSE]

  detected_after <- gene_detected_cells[keep_gene]

  if (!all(candidate_tfs %in% rownames(counts_genefilter))) {
    stop("Candidate TF(s) lost after >=25-cell gene QC: ", paste(setdiff(candidate_tfs, rownames(counts_genefilter)),
      collapse = ", "))
  }

  libsize <- Matrix::colSums(counts_genefilter)

  if (any(!is.finite(libsize)) || any(libsize <= 0)) {
    stop("Invalid library size after scQC.")
  }

  scale_factor <- 10000/libsize

  lognorm <- counts_genefilter %*% Matrix::Diagonal(x = scale_factor)

  lognorm@x <- log1p(lognorm@x)

  mu <- Matrix::rowMeans(lognorm)

  mu2 <- Matrix::rowMeans(lognorm^2)

  gene_var <- pmax(mu2 - mu^2, 0)

  names(gene_var) <- rownames(counts_genefilter)

  required_keep_genes <- unique(c(candidate_tfs, functional_genes))

  required_keep_present <- intersect(required_keep_genes, rownames(counts_genefilter))

  if (nrow(counts_genefilter) <= max_network_genes) {

    selected_genes <- rownames(counts_genefilter)

    hvg_capped <- FALSE

  } else {

    ord <- order(-gene_var, names(gene_var))

    selected_genes <- names(gene_var)[ord[seq_len(max_network_genes)]]

    selected_genes <- unique(c(required_keep_present, selected_genes))

    if (length(selected_genes) > max_network_genes) {

      nonrequired <- setdiff(selected_genes, required_keep_present)

      nonrequired <- nonrequired[order(-gene_var[nonrequired], nonrequired)]

      n_nonrequired <- max_network_genes - length(required_keep_present)

      selected_genes <- c(required_keep_present, head(nonrequired, n_nonrequired))
    }
  }

  if (!all(candidate_tfs %in% selected_genes)) {
    stop("One or more candidate TFs are absent from the network gene set.")
  }

  counts_network <- counts_genefilter[selected_genes, , drop = FALSE]

  gene_selection <- data.frame(gene = rownames(counts_genefilter), detected_cells = as.numeric(detected_after[rownames(counts_genefilter)]),
    logCPM_variance = as.numeric(gene_var[rownames(counts_genefilter)]), selected_for_network = rownames(counts_genefilter) %in%
      selected_genes, is_candidate_TF = rownames(counts_genefilter) %in% candidate_tfs, is_final20_functional_gene = rownames(counts_genefilter) %in%
      functional_genes, is_frozen1320_signature = rownames(counts_genefilter) %in% frozen_signature,
    stringsAsFactors = FALSE)

  saveRDS(counts_network, cp_counts_network, compress = FALSE)

  write.csv(gene_selection, gzfile(cp_gene_selection), row.names = FALSE)

  write.csv(gene_selection, gzfile(gene_selection_file), row.names = FALSE)

  elapsed_cp01 <- as.numeric(difftime(Sys.time(), t_cp01, units = "secs"))

  record_timing("CP01_input_preparation", elapsed_cp01)

  log_msg("CP01 saved: ", nrow(counts_network), " genes x ", ncol(counts_network), " cells; elapsed ",
    fmt_duration(elapsed_cp01), ".")

  rm(lognorm, counts_target, counts_cellqc, counts_genefilter)

  gc()
}

if (!file.exists(gene_selection_file)) {
  write.csv(gene_selection, gzfile(gene_selection_file), row.names = FALSE)
}

strictDirection <- getFromNamespace("strictDirection", "scTenifoldKnk")

dRegulation_v103 <- getFromNamespace("dRegulation", "scTenifoldKnk")

wt_runtime_min <- NA_real_

if (checkpoint_exists(cp_WT_network)) {

  log_msg("Checkpoint CP03 found: loading final WT network; skipping makeNetworks() and tensor decomposition.")

  WT <- as.matrix(readRDS(cp_WT_network))

  if (!all(candidate_tfs %in% rownames(WT))) {
    stop("CP03 WT checkpoint is invalid: candidate TF(s) missing.")
  }

} else {

  if (checkpoint_exists(cp_WT_ensemble)) {

    log_msg("Checkpoint CP02 found: loading WT network ensemble; skipping makeNetworks().")

    WT_ensemble <- readRDS(cp_WT_ensemble)

  } else {

    t_make <- Sys.time()

    set.seed(seed_WT)

    WT_ensemble <- scTenifoldNet::makeNetworks(X = counts_network, q = nc_q, nNet = nc_nNet, nCells = nc_nCells,
      scaleScores = nc_scaleScores, symmetric = nc_symmetric, nComp = nc_nComp, nCores = n_cores)

    elapsed_make <- as.numeric(difftime(Sys.time(), t_make, units = "secs"))

    saveRDS(WT_ensemble, cp_WT_ensemble, compress = FALSE)

    record_timing("CP02_makeNetworks", elapsed_make)

    log_msg("CP02 saved; makeNetworks elapsed ", fmt_duration(elapsed_make), ".")
  }

  t_td <- Sys.time()

  set.seed(seed_WT)

  WT_td <- scTenifoldNet::tensorDecomposition(xList = WT_ensemble, K = td_K, maxError = td_maxError,
    maxIter = td_maxIter, nDecimal = td_nDecimal)

  WT <- WT_td$X

  WT <- strictDirection(WT, lambda = nc_lambda)

  WT <- as.matrix(WT)

  diag(WT) <- 0

  WT <- t(WT)

  elapsed_td <- as.numeric(difftime(Sys.time(), t_td, units = "secs"))

  saveRDS(Matrix::Matrix(WT, sparse = TRUE), cp_WT_network, compress = "gzip")

  saveRDS(Matrix::Matrix(WT, sparse = TRUE), wt_network_file, compress = "gzip")

  record_timing("CP03_tensor_WT", elapsed_td)

  log_msg("CP03 saved; tensor/WT elapsed ", fmt_duration(elapsed_td), ".")

  rm(WT_ensemble, WT_td)

  gc()
}

if (!file.exists(wt_network_file)) {
  saveRDS(Matrix::Matrix(WT, sparse = TRUE), wt_network_file, compress = "gzip")
}

dr_list <- vector("list", length(candidate_tfs))

names(dr_list) <- candidate_tfs

ko_summary_list <- vector("list", length(candidate_tfs))

names(ko_summary_list) <- candidate_tfs

completed_ko_seconds <- numeric(0)

for (i in seq_along(candidate_tfs)) {

  tf <- candidate_tfs[i]

  per_tf_all <- file.path(outdir, paste0("KO_", tf, "_diffRegulation_v04.csv.gz"))

  per_tf_sig <- file.path(outdir, paste0("KO_", tf, "_FDR05_v04.csv"))

  per_tf_summary <- file.path(checkpoint_dir, paste0("KO_", tf, "_summary_checkpoint_v04.csv"))

  if (checkpoint_exists(per_tf_all) && checkpoint_exists(per_tf_summary)) {

    log_msg("KO ", tf, ": checkpoint found -> SKIP and load.")

    DR <- read.csv(gzfile(per_tf_all), stringsAsFactors = FALSE, check.names = FALSE)

    dr_list[[tf]] <- DR

    ko_summary_list[[tf]] <- read.csv(per_tf_summary, stringsAsFactors = FALSE, check.names = FALSE)

    hist <- read_timing()
    stage_name <- paste0("KO_", tf)

    if (stage_name %in% hist$stage) {
      completed_ko_seconds <- c(completed_ko_seconds, hist$seconds[match(stage_name, hist$stage)])
    }

    next
  }

  if (length(completed_ko_seconds) > 0) {

    remaining_including_current <- sum(!vapply(candidate_tfs, function(z) {
      checkpoint_exists(file.path(outdir, paste0("KO_", z, "_diffRegulation_v04.csv.gz"))) && checkpoint_exists(file.path(checkpoint_dir,
        paste0("KO_", z, "_summary_checkpoint_v04.csv")))
    }, logical(1)))

    eta_sec <- mean(completed_ko_seconds) * remaining_including_current

    show_eta_from_seconds("KO-stage ETA based on completed KOs", eta_sec)

  } else {

  }

  outgoing_nonzero <- sum(WT[tf, , drop = FALSE] != 0)

  if (outgoing_nonzero == 0) {
    stop(tf, " has no outgoing edges in the WT network.")
  }

  KO <- WT

  KO[tf, ] <- 0

  set.seed(seed_MA_base + i)

  t_start <- Sys.time()

  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = ma_nDim, nCores = n_cores)

  DR <- dRegulation_v103(MA, tf)

  runtime_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  runtime_min <- runtime_sec/60

  DR <- rank_perturbation(DR, tf)

  DR$is_frozen1320_signature <- DR$gene %in% frozen_signature

  DR$is_final20_functional_gene <- DR$gene %in% functional_genes

  DR$is_candidate_TF <- DR$gene %in% candidate_tfs

  DR$FDR05 <- !is.na(DR$p.adj) & DR$p.adj < fdr_cutoff

  dr_list[[tf]] <- DR

  ko_row <- DR[DR$gene == tf, , drop = FALSE]

  nonko <- DR[DR$gene != tf, , drop = FALSE]

  sig_nonko <- nonko[nonko$FDR05, , drop = FALSE]

  top50_nonko <- head(nonko, top_n_export)

  summary_row <- data.frame(KO_TF = tf, WT_outgoing_edges = outgoing_nonzero, KO_gene_rank = if (nrow(ko_row) ==
    1) {
    ko_row$perturbation_rank
  } else {
    NA_integer_
  }, KO_gene_distance = if (nrow(ko_row) == 1) {
    ko_row$distance
  } else {
    NA_real_
  }, KO_gene_FC = if (nrow(ko_row) == 1) {
    ko_row$FC
  } else {
    NA_real_
  }, n_FDR05_nonKO = nrow(sig_nonko), n_top50_in_frozen1320_signature = sum(top50_nonko$is_frozen1320_signature),
    n_top50_in_final20_functional = sum(top50_nonko$is_final20_functional_gene), max_nonKO_distance = max(nonko$distance,
      na.rm = TRUE), median_top50_log2_perturbation_FC = stats::median(top50_nonko$log2_perturbation_FC,
      na.rm = TRUE), runtime_min = runtime_min, stringsAsFactors = FALSE)

  ko_summary_list[[tf]] <- summary_row

  write.csv(DR, gzfile(per_tf_all), row.names = FALSE)

  write.csv(sig_nonko, per_tf_sig, row.names = FALSE)

  write.csv(summary_row, per_tf_summary, row.names = FALSE)

  record_timing(paste0("KO_", tf), runtime_sec)

  completed_ko_seconds <- c(completed_ko_seconds, runtime_sec)

  log_msg("KO ", tf, " COMPLETE: outgoing=", outgoing_nonzero, "; FDR<0.05 non-KO=", nrow(sig_nonko),
    "; elapsed=", fmt_duration(runtime_sec), ".")

  remaining_after <- length(candidate_tfs) - i

  if (remaining_after > 0 && length(completed_ko_seconds) > 0) {
    show_eta_from_seconds("Remaining KO ETA", mean(completed_ko_seconds) * remaining_after)
  }

  rm(KO, MA, DR, nonko, sig_nonko, top50_nonko, ko_row)

  gc()
}

combined_dr <- do.call(rbind, dr_list)

rownames(combined_dr) <- NULL

write.csv(combined_dr, gzfile(combined_dr_file), row.names = FALSE)

significant_dr <- combined_dr[!combined_dr$is_KO_gene & combined_dr$FDR05, , drop = FALSE]

write.csv(significant_dr, gzfile(significant_dr_file), row.names = FALSE)

top50 <- do.call(rbind, lapply(candidate_tfs, function(tf) {
    z <- combined_dr[combined_dr$KO_TF == tf & !combined_dr$is_KO_gene, , drop = FALSE]
    head(z, top_n_export)
}))

rownames(top50) <- NULL

write.csv(top50, top50_file, row.names = FALSE)

ko_summary <- do.call(rbind, ko_summary_list)

rownames(ko_summary) <- NULL

write.csv(ko_summary, ko_summary_file, row.names = FALSE)

library(ggplot2)

plot_df <- do.call(rbind, lapply(candidate_tfs, function(tf) {
    z <- combined_dr[combined_dr$KO_TF == tf & !combined_dr$is_KO_gene, , drop = FALSE]
    head(z, top_n_plot)
}))

rownames(plot_df) <- NULL

plot_df$label <- paste0(plot_df$gene, " | ", plot_df$KO_TF)

plot_df$label <- factor(plot_df$label, levels = rev(unique(plot_df$label)))

p6f <- ggplot(plot_df, aes(x = log2_perturbation_FC, y = label)) + geom_col(width = 0.75) + facet_wrap(~KO_TF,
  scales = "free_y", ncol = 5) + labs(x = "log2 perturbation FC", y = NULL) + theme_classic(base_size = 9) +
  theme(strip.text = element_text(face = "bold"), axis.text.y = element_text(size = 6), axis.title.x = element_text(face = "bold"))

ggsave(figure6f_pdf, p6f, width = 13, height = 6, units = "in", device = cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(figure6f_svg, p6f, width = 13, height = 6, units = "in", device = svglite::svglite)
} else {
    svg(figure6f_svg, width = 13, height = 6)
    print(p6f)
    dev.off()
}

ggsave(figure6f_png, p6f, width = 13, height = 6, units = "in", dpi = png_dpi)

bundle <- list(candidate_TFs = candidate_tfs, functional_genes = functional_genes, parameters = list(scTenifoldKnk_version = knk_version,
  scTenifoldNet_version = net_version, qc_mtThreshold = qc_mtThreshold, qc_minLSize = qc_minLSize,
  qc_minCells_manual = qc_minCells, max_network_genes = max_network_genes, nc_lambda = nc_lambda, nc_nNet = nc_nNet,
  nc_nCells = nc_nCells, nc_nComp = nc_nComp, nc_scaleScores = nc_scaleScores, nc_symmetric = nc_symmetric,
  nc_q = nc_q, td_K = td_K, td_maxIter = td_maxIter, td_maxError = td_maxError, td_nDecimal = td_nDecimal,
  ma_nDim = ma_nDim, n_cores = n_cores, seed_WT = seed_WT), gene_selection = gene_selection, KO_summary = ko_summary,
  top50_perturbed = top50)

saveRDS(bundle, bundle_file, compress = "gzip")

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("scTenifoldKnk version", if (knk_version == "1.0.3") "PASS" else "FAIL",
  knk_version)

qc_list[[length(qc_list) + 1]] <- make_qc("scTenifoldNet version", if (net_version == "1.3") "PASS" else "FAIL",
  net_version)

qc_list[[length(qc_list) + 1]] <- make_qc("Internal scQC resolution", if (is.function(scQC_v103) && all(required_scqc_args %in%
  scqc_formals)) "PASS" else "FAIL", paste0("getFromNamespace; formals=", paste(scqc_formals, collapse = ",")))

qc_list[[length(qc_list) + 1]] <- make_qc("Step90 v05 QC", "PASS", "No blocking Step90 v05 FAIL items.")

qc_list[[length(qc_list) + 1]] <- make_qc("Target-state input", if (n_target_cells_frozen == 17713L) "PASS" else "FAIL",
  paste0(n_target_cells_frozen, " / 17713"))

qc_list[[length(qc_list) + 1]] <- make_qc("Cells after scQC", if (ncol(counts_network) >= nc_nCells) "PASS" else "FAIL",
  paste0(ncol(counts_network), " cells"))

qc_list[[length(qc_list) + 1]] <- make_qc("Manual >=25-cell gene filter", "PASS", "Applied explicitly to work around the v1.0.3 qc_minCells assignment bug.")

qc_list[[length(qc_list) + 1]] <- make_qc("Network gene count", if (nrow(counts_network) <= max_network_genes) "PASS" else "FAIL",
  paste0(nrow(counts_network), " / cap ", max_network_genes))

qc_list[[length(qc_list) + 1]] <- make_qc("All 5 KO TFs in network", if (all(candidate_tfs %in% rownames(WT))) "PASS" else "FAIL",
  paste(candidate_tfs, collapse = ", "))

qc_list[[length(qc_list) + 1]] <- make_qc("Five independent KOs completed", if (length(dr_list) == 5 &&
  all(vapply(dr_list, function(x) {
    is.data.frame(x) && nrow(x) > 0
  }, logical(1)))) "PASS" else "FAIL", paste(candidate_tfs, collapse = ", "))

qc_list[[length(qc_list) + 1]] <- make_qc("KO genes represented in DR output", if (all(ko_summary$KO_gene_rank >=
  1, na.rm = TRUE)) "PASS" else "FAIL", paste0("ranks=", paste(ko_summary$KO_gene_rank, collapse = ", ")))

qc_list[[length(qc_list) + 1]] <- make_qc("Candidate TFs have outgoing WT edges", if (all(ko_summary$WT_outgoing_edges >
  0)) "PASS" else "FAIL", paste0(paste(ko_summary$KO_TF, ko_summary$WT_outgoing_edges, sep = "=", collapse = "; ")))

qc_list[[length(qc_list) + 1]] <- make_qc("Figure6F diagnostic", if (all(file.exists(c(figure6f_pdf,
  figure6f_svg, figure6f_png)))) "PASS" else "FAIL", "PDF + SVG + 600dpi PNG")

qc_list[[length(qc_list) + 1]] <- make_qc("Directionality interpretation", "PASS", "scTenifoldKnk FC is treated as perturbation magnitude, not expression up/down.")

qc_list[[length(qc_list) + 1]] <- make_qc("Checkpoint CP01 network input", if (checkpoint_exists(cp_counts_network)) "PASS" else "FAIL",
  basename(cp_counts_network))

qc_list[[length(qc_list) + 1]] <- make_qc("Checkpoint CP02 WT ensemble", if (checkpoint_exists(cp_WT_ensemble)) "PASS" else "WARN",
  "If makeNetworks was skipped from an existing CP03, CP02 may be unnecessary; otherwise expected.")

qc_list[[length(qc_list) + 1]] <- make_qc("Checkpoint CP03 WT network", if (checkpoint_exists(cp_WT_network)) "PASS" else "FAIL",
  basename(cp_WT_network))

qc_list[[length(qc_list) + 1]] <- make_qc("ETA/timing log", if (file.exists(timing_file) && file.exists(run_log_file)) "PASS" else "WARN",
  "Stage timings and estimated finish times were recorded.")

final_qc <- do.call(rbind, qc_list)

rownames(final_qc) <- NULL

write.csv(final_qc, final_qc_file, row.names = FALSE)

n_pass <- sum(final_qc$status == "PASS")

n_warn <- sum(final_qc$status == "WARN")

n_fail <- sum(final_qc$status == "FAIL")

overall_status <- if (n_fail > 0) {
  "FAIL"
} else if (n_warn > 0) {
  "PASS_WITH_WARNINGS"
} else {
  "PASS"
}

summary_order <- ko_summary[order(-ko_summary$n_FDR05_nonKO, -ko_summary$n_top50_in_frozen1320_signature,
  -ko_summary$max_nonKO_distance), , drop = FALSE]

writeLines("Numerical results are provided in the accompanying output tables.", 
    report_file)

input_qc <- data.frame(item = c("Source myeloid RDS", "Frozen 1320-gene signature", "Step90 v05 functional panel",
  "Step90 v05 QC"), path = c(source_rds, frozen_signature_file, step90_functional_file, step90_qc_file),
  exists = file.exists(c(source_rds, frozen_signature_file, step90_functional_file, step90_qc_file)),
  stringsAsFactors = FALSE)

write.csv(input_qc, input_qc_file, row.names = FALSE)

output_files <- c(environment_file, input_qc_file, gene_selection_file, wt_network_file, combined_dr_file,
  significant_dr_file, top50_file, ko_summary_file, figure6f_pdf, figure6f_svg, figure6f_png, bundle_file,
  final_qc_file, report_file)

output_index <- data.frame(file = basename(output_files), exists = file.exists(output_files), size_MB = ifelse(file.exists(output_files),
  file.info(output_files)$size/1024^2, NA_real_), full_path = normalizePath(output_files, winslash = "/",
  mustWork = FALSE), stringsAsFactors = FALSE)

write.csv(output_index, output_index_file, row.names = FALSE)

{
  hist_tmp <- read_timing()
  wt_sec_tmp <- sum(hist_tmp$seconds[hist_tmp$stage %in% c("CP02_makeNetworks", "CP03_tensor_WT")],
    na.rm = TRUE)
  cat("WT runtime    = ", if (is.finite(wt_sec_tmp) && wt_sec_tmp > 0) {
    round(wt_sec_tmp/60, 2)
  } else {
    NA_real_
  }, " min\n", sep = "")
}
cat("WARN          = ", n_warn, "\n", sep = "")
cat("FAIL          = ", n_fail, "\n", sep = "")

print(ko_summary, row.names = FALSE)

if (n_fail > 0) {
  stop("Step 91 contains blocking FAIL item(s).")
}
