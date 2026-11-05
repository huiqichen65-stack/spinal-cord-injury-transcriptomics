source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, width = 160, scipen = 50)

ROOT <- PROJECT_ROOT
VERSION <- "v02"
STEP <- "135"
OUT_FOLDER <- "135_TF_candidate_comparison_resampling_v02"
OUT <- file.path(ROOT, "results", OUT_FOLDER)

RUN_NETWORK <- TRUE
NETWORK_WORKERS <- 6L
NETWORK_NNET <- 10L
NETWORK_NCELLS <- 500L
BOOTSTRAP_B <- 500L
BOOTSTRAP_SEED <- 135001L
MIN_FOCUS_CELLS <- 10L
MIN_COMPARATOR_CELLS <- 20L
MIN_LIPID_CORR_CELLS <- 30L
FDR <- 0.05
TOP_K_STABILITY <- 50L
OVERWRITE_ANALYTIC_CSV <- TRUE

RUN_SOURCE_SPECIFIC <- TRUE
IMPORT_V01_S01_IF_COMPLETE <- TRUE
REQUIRE_V01_S00_REFERENCE <- TRUE

TARGET <- "Clec7a+ lipid-handling Microglia"
TARGET_LEGACY <- "Clec7a+ lipid-associated Microglia"
ORIGINAL_FIVE <- c("Mef2a", "Nr3c1", "Fli1", "Zfp148", "Bclaf1")
EXPECTED_TEN <- c("Mef2a", "Nr3c1", "Pura", "Ikzf1", "Bach1", "Fli1", "Zfp148", "Ep300", "Bclaf1", "Chd1")

P <- list(meta = file.path(ROOT, "results/84_pySCENIC_input_preparation/09_84_pySCENIC_cell_metadata_v01.csv.gz"),
  regulon_auc = file.path(ROOT, "results/88_pySCENIC_regulon_AUCell_full_atlas/11_88_pySCENIC_AUCell_full_atlas_v02.csv.gz"),
  rss = file.path(ROOT, "results/89_pySCENIC_RSS_target_regulon_ranking/05_89_Clec7a_lipid_handling_MG_regulon_RSS_ranking_v01.csv"),
  evidence = file.path(ROOT, "results/90_pySCENIC_top10_TF_network_signature_correlation/10_90_top10_TF_evidence_summary_v05.csv"),
  targets = file.path(ROOT, "results/90_pySCENIC_top10_TF_network_signature_correlation/03_90_top10_regulon_targets_full_v03.csv"),
  lipid_auc = file.path(ROOT, "results/20_AUCell_lipid_metabolism_scoring/20_AUCell_cell_scores_v03.csv"),
  lipid_genes = file.path(ROOT, "results/20_AUCell_lipid_metabolism_scoring/20_AUCell_lipid_gene_set_manifest_v03.csv"),
  ko_summary = file.path(ROOT, "results/91_scTenifoldKnk_top5_virtual_KO/08_91_top5_KO_summary_v04.csv"),
  ko_lipid = file.path(ROOT, "results/92_scTenifoldKnk_Mef2a_functional_enrichment/20_92_top5_TF_lipid_perturbation_enrichment_v02.csv"),
  cp_counts = file.path(ROOT, "results/91_scTenifoldKnk_top5_virtual_KO/checkpoints_v04/CP01_counts_network_v04.rds"),
  cp_genes = file.path(ROOT, "results/91_scTenifoldKnk_top5_virtual_KO/checkpoints_v04/CP01_gene_selection_v04.csv.gz"),
  cp_wt = file.path(ROOT, "results/91_scTenifoldKnk_top5_virtual_KO/checkpoints_v04/CP03_WT_network_v04.rds"))

if (!dir.exists(ROOT)) stop("Project ROOT does not exist: ", ROOT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NETDIR <- file.path(OUT, "network_checkpoints_v02")
dir.create(NETDIR, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "135_TF_candidate_comparison_resampling_v02.log")

log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG, append = TRUE)
  invisible(msg)
}
stopf <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
read_tab <- function(path, select = NULL) {
  if (!file.exists(path))
    stopf("Missing required input: ", path)
  x <- NULL
  if (requireNamespace("data.table", quietly = TRUE)) {
    x <- tryCatch(data.table::fread(path, select = select, data.table = FALSE, showProgress = FALSE),
      error = function(e) {
        log_msg("fread fallback for ", basename(path), ": ", conditionMessage(e))
        NULL
      })
  }
  if (is.null(x)) {
    x <- if (grepl("\\.gz$", path, ignore.case = TRUE))
      utils::read.csv(gzfile(path), check.names = FALSE) else utils::read.csv(path, check.names = FALSE)
    if (!is.null(select))
      x <- x[, select, drop = FALSE]
  }
  as.data.frame(x, check.names = FALSE)
}
write_tab <- function(x, basename) {
  path <- file.path(OUT, basename)
  if (!OVERWRITE_ANALYTIC_CSV && file.exists(path))
    stopf("Refusing to overwrite: ", path)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}
require_cols <- function(x, cols, name) {
  missing <- setdiff(cols, names(x))
  if (length(missing))
    stopf(name, " missing columns: ", paste(missing, collapse = ", "))
}
clean_state <- function(x) {
  x <- as.character(x)
  x[x == TARGET_LEGACY] <- TARGET
  x
}
as_num <- function(x) suppressWarnings(as.numeric(x))
finite_med <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x))
    NA_real_ else median(x)
}
finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x))
    NA_real_ else mean(x)
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5L || length(unique(x[ok])) < 3L || length(unique(y[ok])) < 3L)
    return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}
partial_spearman <- function(x, y, z_depth, z_features) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z_depth) & is.finite(z_features)
  if (sum(ok) < 10L || sd(x[ok]) == 0 || sd(y[ok]) == 0)
    return(NA_real_)
  rx <- rank(x[ok], ties.method = "average")
  ry <- rank(y[ok], ties.method = "average")
  rz_depth <- rank(z_depth[ok], ties.method = "average")
  rz_features <- rank(z_features[ok], ties.method = "average")
  design <- cbind(1, rz_depth, rz_features)
  xr <- stats::lm.fit(design, rx)$residuals
  yr <- stats::lm.fit(design, ry)$residuals
  if (!is.finite(sd(xr)) || !is.finite(sd(yr)) || sd(xr) == 0 || sd(yr) == 0)
    return(NA_real_)
  stats::cor(xr, yr)
}
write_text <- function(lines, basename) {
  writeLines(as.character(lines), file.path(OUT, basename), useBytes = TRUE)
}
if (file.exists(LOG)) file.remove(LOG)

primary <- c("meta", "regulon_auc", "rss", "evidence", "targets", "lipid_auc", "lipid_genes", "ko_summary",
  "ko_lipid")
audit <- data.frame(input = names(P), original_path = unlist(P, use.names = FALSE), exists = file.exists(unlist(P,
  use.names = FALSE)), size_bytes = unname(file.info(unlist(P, use.names = FALSE))$size))
write_tab(audit, "00_135_INPUT_AUDIT_v02.csv")
if (any(!audit$exists[match(primary, audit$input)])) {
  stopf("Required original project inputs missing: ", paste(primary[!audit$exists[match(primary, audit$input)]],
    collapse = ", "))
}

evi <- read_tab(P$evidence)
require_cols(evi, c("TF", "regulon", "RSS_rank", "target_RSS", "RSS_margin_vs_best_other", "target_mean_AUC",
  "target_has_highest_mean_AUC", "regulon_size"), "Step 90")
evi$RSS_rank <- as.integer(evi$RSS_rank)
evi <- evi[order(evi$RSS_rank, evi$TF), , drop = FALSE]
if (anyDuplicated(evi$TF) || nrow(evi) != 10L || !setequal(evi$TF, EXPECTED_TEN)) {
  stopf("The top-ten roster changed; inspect Step 90 before running 135.")
}
TFS <- as.character(evi$TF)
roster <- data.frame(evi, selected_in_Step91 = evi$TF %in% ORIGINAL_FIVE, selection_rule = "Top 10 target-state RSS, common comparison; Step91 five descriptive only",
  stringsAsFactors = FALSE)
write_tab(roster, "01_135_top10_RSS_frozen_candidate_table_v02.csv")

reg_targets <- read_tab(P$targets)
require_cols(reg_targets, c("TF", "target", "regulon_weight", "target_in_frozen_UP_signature"), "Step 90 targets")
reg_targets <- unique(reg_targets[, c("TF", "target", "regulon_weight", "target_in_frozen_UP_signature")])
reg_targets <- reg_targets[reg_targets$TF %in% TFS & !is.na(reg_targets$target), , drop = FALSE]
if (!all(TFS %in% reg_targets$TF)) stopf("Missing a TF in Step90 target list")
lipid_manifest <- read_tab(P$lipid_genes)
require_cols(lipid_manifest, c("gene_set", "gene"), "Step 20 lipid manifest")
lipid_manifest <- unique(lipid_manifest[, c("gene_set", "gene")])
lipid_manifest <- lipid_manifest[!is.na(lipid_manifest$gene) & nzchar(lipid_manifest$gene), , drop = FALSE]
lipid_gene_union <- unique(lipid_manifest$gene)
reg_targets$in_lipid_manifest <- reg_targets$target %in% lipid_gene_union
reg_targets$lipid_categories <- vapply(reg_targets$target, function(g) {
  paste(sort(unique(lipid_manifest$gene_set[lipid_manifest$gene == g])), collapse = ";")
}, character(1))
write_tab(reg_targets, "02_135_top10_motif_pruned_targets_lipid_annotation_v02.csv")

target_stat <- do.call(rbind, lapply(TFS, function(tf) {
  a <- reg_targets[reg_targets$TF == tf, , drop = FALSE]
  unique_genes <- unique(a$target)
  data.frame(TF = tf, n_motif_pruned_targets = length(unique_genes), n_targets_in_frozen_signature = sum(a$target_in_frozen_UP_signature %in%
    c(TRUE, "TRUE"), na.rm = TRUE), n_lipid_manifest_targets = sum(a$in_lipid_manifest), fraction_lipid_targets = mean(a$in_lipid_manifest),
    median_regulon_weight = finite_med(as_num(a$regulon_weight)))
}))
target_categories <- do.call(rbind, lapply(TFS, function(tf) {
  genes <- unique(reg_targets$target[reg_targets$TF == tf])
  do.call(rbind, lapply(sort(unique(lipid_manifest$gene_set)), function(gs) {
    hits <- intersect(genes, unique(lipid_manifest$gene[lipid_manifest$gene_set == gs]))
    data.frame(TF = tf, lipid_gene_set = gs, n_set_genes = length(unique(lipid_manifest$gene[lipid_manifest$gene_set ==
      gs])), n_regulon_targets = length(genes), n_overlap = length(hits), overlap_genes = paste(sort(hits),
      collapse = ";"))
  }))
}))
write_tab(target_stat, "03_135_top10_target_set_properties_v02.csv")
write_tab(target_categories, "04_135_top10_lipid_gene_set_target_overlap_v02.csv")
overlap_mat <- do.call(rbind, lapply(TFS, function(a) do.call(rbind, lapply(TFS, function(b) {
  x <- unique(reg_targets$target[reg_targets$TF == a])
  y <- unique(reg_targets$target[reg_targets$TF == b])
  u <- union(x, y)
  data.frame(TF_A = a, TF_B = b, n_shared = length(intersect(x, y)), n_union = length(u), jaccard = if (length(u))
    length(intersect(x, y))/length(u) else NA_real_)
}))))
write_tab(overlap_mat, "05_135_top10_regulon_target_pairwise_Jaccard_v02.csv")

meta_cols <- c("cell", "source_dataset", "sample_id", "sample_id_source", "time", "treatment", "myeloid_state_frozen",
  "nCount_RNA", "nFeature_RNA")
meta <- read_tab(P$meta, select = meta_cols)
require_cols(meta, meta_cols, "Step84 metadata")
if (anyDuplicated(meta$cell)) stopf("Step84 duplicate cell barcodes")
meta$myeloid_state_frozen <- clean_state(meta$myeloid_state_frozen)
meta$source_dataset <- as.character(meta$source_dataset)
meta$sample_id_source <- as.character(meta$sample_id_source)
meta$sample_id <- as.character(meta$sample_id)
meta$time <- as.character(meta$time)
meta$unit <- paste(meta$source_dataset, meta$sample_id, sep = "::")
meta$unit_category <- ifelse(meta$source_dataset == "GSE162610" & meta$sample_id_source == "metadata:sample_id",
  "sample_label_independence_unconfirmed", ifelse(meta$sample_id_source == "derived_from_timepoint",
    "timepoint_pseudounit_NOT_replication", "library_or_pool_independence_unconfirmed"))
meta$depth <- as_num(meta$nCount_RNA)
meta$log_depth <- log1p(pmax(meta$depth, 0))
meta$log_features <- log1p(pmax(as_num(meta$nFeature_RNA), 0))
meta$is_focus <- meta$myeloid_state_frozen == TARGET
meta$is_microglia <- grepl("Microglia", meta$myeloid_state_frozen, ignore.case = TRUE)

auc_columns <- paste0(TFS, "(+)")
auc <- read_tab(P$regulon_auc, select = c("Cell", auc_columns))
require_cols(auc, c("Cell", auc_columns), "Step88 AUCell")
if (anyDuplicated(auc$Cell)) stopf("Step88 duplicate Cell identifiers")
match_auc <- match(meta$cell, auc$Cell)
if (anyNA(match_auc) || length(match_auc) != nrow(auc)) stopf("Step84 and Step88 cell IDs do not align one-to-one")
for (i in seq_along(TFS)) meta[[paste0("TF_", TFS[i])]] <- as_num(auc[[auc_columns[i]]][match_auc])
rm(auc)
gc(verbose = FALSE)

lip <- read_tab(P$lipid_auc)
require_cols(lip, "cell", "Step20 lipid AUCell")
lipid_cols <- grep("^AUC_", names(lip), value = TRUE)
lipid_programs <- setdiff(lipid_cols, "AUC_LipidCompositeZ")
if (length(lipid_programs) < 8L) stopf("Unexpected lipid program schema: ", paste(lipid_cols, collapse = ","))
if (anyDuplicated(lip$cell)) stopf("Step20 duplicate cell IDs")
match_lip <- match(meta$cell, lip$cell)
if (anyNA(match_lip)) stopf("Step20 lipid scores missing cells from Step84")
for (nm in lipid_programs) meta[[nm]] <- as_num(lip[[nm]][match_lip])
rm(lip)
gc(verbose = FALSE)

unit_list <- split(seq_len(nrow(meta)), meta$unit)
unit_meta <- do.call(rbind, lapply(unit_list, function(ix) {
  z <- meta[ix, , drop = FALSE]
  data.frame(source_dataset = z$source_dataset[1], sample_id = z$sample_id[1], unit = z$unit[1], unit_category = z$unit_category[1],
    time = z$time[1], treatment = z$treatment[1], n_cells = nrow(z), n_focus = sum(z$is_focus), n_other_microglia = sum(z$is_microglia &
      !z$is_focus), n_other_myeloid = sum(!z$is_focus), qualified_focus = sum(z$is_focus) >= MIN_FOCUS_CELLS,
    stringsAsFactors = FALSE)
}))
row.names(unit_meta) <- NULL
write_tab(unit_meta, "06_135_unit_provenance_and_cell_eligibility_v02.csv")
stage_counts <- as.data.frame(table(meta$source_dataset, meta$time, meta$unit_category), responseName = "n_cells")
names(stage_counts)[1:3] <- c("source_dataset", "time", "unit_category")
write_tab(stage_counts, "07_135_source_stage_unit_type_counts_v02.csv")

unit_diffs <- list()
j <- 0L
for (u in names(unit_list)) {
  ix <- unit_list[[u]]
  fo <- ix[meta$is_focus[ix]]
  allothers <- ix[!meta$is_focus[ix]]
  mg <- ix[meta$is_microglia[ix] & !meta$is_focus[ix]]
  for (type in c("other_microglia", "other_myeloid")) {
    comp <- if (type == "other_microglia")
      mg else allothers
    eligible <- length(fo) >= MIN_FOCUS_CELLS && length(comp) >= MIN_COMPARATOR_CELLS
    for (tf in TFS) {
      v <- meta[[paste0("TF_", tf)]]
      j <- j + 1L
      unit_diffs[[j]] <- data.frame(TF = tf, unit = u, source_dataset = meta$source_dataset[ix[1]],
        time = meta$time[ix[1]], unit_category = meta$unit_category[ix[1]], comparator = type,
        n_focus = length(fo), n_comparator = length(comp), eligible = eligible, mean_auc_focus = if (length(fo))
          finite_mean(v[fo]) else NA_real_, mean_auc_comparator = if (length(comp))
          finite_mean(v[comp]) else NA_real_, delta_mean_auc = if (eligible)
          finite_mean(v[fo]) - finite_mean(v[comp]) else NA_real_)
    }
  }
}
unit_diffs <- do.call(rbind, unit_diffs)
row.names(unit_diffs) <- NULL
write_tab(unit_diffs, "08_135_per_unit_top10_regulon_activity_differences_v02.csv")

unit_summary <- do.call(rbind, lapply(split(unit_diffs, interaction(unit_diffs$TF, unit_diffs$comparator,
  unit_diffs$unit_category, drop = TRUE)), function(a) {
  ok <- a[a$eligible & is.finite(a$delta_mean_auc), , drop = FALSE]
  data.frame(TF = a$TF[1], comparator = a$comparator[1], unit_category = a$unit_category[1], n_eligible_units = nrow(ok),
    n_source_datasets = length(unique(ok$source_dataset)), n_stages = length(unique(ok$time)), mean_unit_delta = finite_mean(ok$delta_mean_auc),
    median_unit_delta = finite_med(ok$delta_mean_auc), fraction_positive = if (nrow(ok))
      mean(ok$delta_mean_auc > 0) else NA_real_)
}))
row.names(unit_summary) <- NULL
write_tab(unit_summary, "09_135_per_unit_activity_summary_by_provenance_v02.csv")

valid_boot <- subset(unit_diffs, unit_category == "sample_label_independence_unconfirmed" & eligible &
  is.finite(delta_mean_auc))
set.seed(BOOTSTRAP_SEED)
boot_rows <- list()
bcount <- 0L
for (comp in c("other_microglia", "other_myeloid")) {
  for (tf in TFS) {
    a <- valid_boot[valid_boot$TF == tf & valid_boot$comparator == comp, , drop = FALSE]
    stages <- sort(unique(a$time))
    if (!length(stages))
      next

    estimate <- mean(vapply(stages, function(s) mean(a$delta_mean_auc[a$time == s]), numeric(1)))
    per_stage_n <- vapply(stages, function(s) sum(a$time == s), integer(1))

    has_repeat <- length(stages) >= 2L && any(per_stage_n >= 2L)
    bs <- rep(NA_real_, BOOTSTRAP_B)
    if (has_repeat)
      for (b in seq_len(BOOTSTRAP_B)) {
        st <- vapply(stages, function(s) {
          v <- a$delta_mean_auc[a$time == s]
          mean(sample(v, length(v), replace = TRUE))
        }, numeric(1))
        bs[b] <- mean(st)
      }
    bcount <- bcount + 1L
    boot_rows[[bcount]] <- data.frame(TF = tf, comparator = comp, n_eligible_labelled_units = nrow(a),
      n_stages = length(stages), stage_count_detail = paste(vapply(stages, function(s) paste0(s,
        ":", sum(a$time == s)), character(1)), collapse = ";"), stage_equal_weight_delta = estimate,
      n_stages_with_at_least_two_labels = sum(per_stage_n >= 2L), n_singleton_stages = sum(per_stage_n ==
        1L), bootstrapped_B = if (has_repeat)
        BOOTSTRAP_B else 0L, bootstrap_q025 = if (has_repeat)
        unname(quantile(bs, 0.025)) else NA_real_, bootstrap_q975 = if (has_repeat)
        unname(quantile(bs, 0.975)) else NA_real_, bootstrap_positive_fraction = if (has_repeat)
        mean(bs > 0) else NA_real_, inference = "DESCRIPTIVE: sample labels are not verified animals; singleton stages held fixed in bootstrap")
  }
}
boot_out <- if (length(boot_rows)) do.call(rbind, boot_rows) else data.frame()
write_tab(boot_out, "10_135_stage_stratified_labelled_unit_bootstrap_v02.csv")

assoc <- list()
j <- 0L
for (u in names(unit_list)) {
  ix <- unit_list[[u]]
  ix <- ix[meta$is_focus[ix]]
  if (!length(ix))
    next
  for (tf in TFS) {
    x <- meta[[paste0("TF_", tf)]][ix]
    for (pr in lipid_programs) {
      y <- meta[[pr]][ix]
      valid <- is.finite(x) & is.finite(y) & is.finite(meta$log_depth[ix]) & is.finite(meta$log_features[ix])
      eligible <- sum(valid) >= MIN_LIPID_CORR_CELLS
      j <- j + 1L
      assoc[[j]] <- data.frame(TF = tf, lipid_program = pr, unit = u, source_dataset = meta$source_dataset[ix[1]],
        time = meta$time[ix[1]], unit_category = meta$unit_category[ix[1]], n_focus_cells = length(ix),
        n_complete_cells = sum(valid), eligible = eligible, rho = if (eligible)
          safe_cor(x[valid], y[valid]) else NA_real_, partial_rho_depth_adjusted = if (eligible)
          partial_spearman(x[valid], y[valid], meta$log_depth[ix][valid], meta$log_features[ix][valid]) else NA_real_)
    }
  }
}
assoc <- if (length(assoc)) do.call(rbind, assoc) else data.frame()
write_tab(assoc, "11_135_within_unit_top10_regulon_lipid_correlations_v02.csv")

summary_assoc <- do.call(rbind, lapply(split(assoc, interaction(assoc$TF, assoc$lipid_program, assoc$unit_category,
  drop = TRUE)), function(a) {
  ok <- a[a$eligible & is.finite(a$partial_rho_depth_adjusted), , drop = FALSE]
  data.frame(TF = a$TF[1], lipid_program = a$lipid_program[1], unit_category = a$unit_category[1],
    n_units = nrow(ok), n_sources = length(unique(ok$source_dataset)), mean_raw_rho = finite_mean(ok$rho),
    median_raw_rho = finite_med(ok$rho), mean_depth_adjusted_rho = finite_mean(ok$partial_rho_depth_adjusted),
    median_depth_adjusted_rho = finite_med(ok$partial_rho_depth_adjusted), positive_unit_fraction = if (nrow(ok))
      mean(ok$partial_rho_depth_adjusted > 0) else NA_real_)
}))
row.names(summary_assoc) <- NULL
write_tab(summary_assoc, "12_135_within_unit_lipid_association_summary_v02.csv")
source_assoc <- do.call(rbind, lapply(split(assoc, interaction(assoc$TF, assoc$lipid_program, assoc$source_dataset,
  drop = TRUE)), function(a) {
  ok <- a[a$eligible & is.finite(a$partial_rho_depth_adjusted), , drop = FALSE]
  data.frame(TF = a$TF[1], lipid_program = a$lipid_program[1], source_dataset = a$source_dataset[1],
    n_units = nrow(ok), n_stages = length(unique(ok$time)), median_partial_rho = finite_med(ok$partial_rho_depth_adjusted),
    median_raw_rho = finite_med(ok$rho))
}))
row.names(source_assoc) <- NULL
write_tab(source_assoc, "13_135_source_resolved_regulon_lipid_association_v02.csv")

ko_old <- read_tab(P$ko_summary)
enrich_old <- read_tab(P$ko_lipid)
write_tab(ko_old, "14_135_frozen_five_KO_summary_reference_v02.csv")
write_tab(enrich_old, "15_135_frozen_five_lipid_enrichment_reference_v02.csv")

basic <- merge(evi, target_stat, by = "TF", all.x = TRUE, sort = FALSE)
bs_main <- boot_out[boot_out$comparator == "other_microglia", intersect(c("TF", "n_eligible_labelled_units",
  "stage_equal_weight_delta", "bootstrap_q025", "bootstrap_q975", "bootstrapped_B"), names(boot_out)),
  drop = FALSE]
basic <- merge(basic, bs_main, by = "TF", all.x = TRUE, sort = FALSE)
unit_lip <- summary_assoc[summary_assoc$unit_category == "sample_label_independence_unconfirmed", , drop = FALSE]
lip_roll <- do.call(rbind, lapply(TFS, function(tf) {
  a <- unit_lip[unit_lip$TF == tf, , drop = FALSE]
  data.frame(TF = tf, n_lipid_programs_reported = nrow(a), median_association_over_programs = finite_med(a$median_depth_adjusted_rho),
    n_programs_positive = if (nrow(a))
      sum(a$median_depth_adjusted_rho > 0, na.rm = TRUE) else NA_integer_)
}))
basic <- merge(basic, lip_roll, by = "TF", all.x = TRUE, sort = FALSE)
basic <- basic[order(basic$RSS_rank), , drop = FALSE]
basic$interpretation <- "Exploratory multi-criterion comparison, NOT a prespecified composite ranking"
write_tab(basic, "16_135_top10_unified_five_criterion_comparison_PRENETWORK_v02.csv")

NETDIR <- file.path(OUT, "network_checkpoints_v02")
dir.create(NETDIR, recursive = TRUE, showWarnings = FALSE)
V01_OUT <- file.path(ROOT, "results", "135_TF_candidate_comparison_resampling")
V01_NET <- file.path(V01_OUT, "network_checkpoints_v01")

scenarios <- data.frame(id = c("H00_v01_original5000", "H01_v01_extra5000_IF_COMPLETE", "B00_seed9101_3000",
  "B01_seed13521_3000", "B02_seed13522_3000", "C01_cells80_13523", "C02_cells80_13524", "C03_cells80_13525",
  "D01_noForced20", "E01_q085", "E02_K4", "F01_source_17702045", "F02_source_GSE196928"), type = c("historical",
  "historical_optional", rep("new", 11)), seed = c(9101L, 13501L, 9101L, 13521L, 13522L, 13523L, 13524L,
  13525L, 9101L, 9101L, 9101L, 13527L, 13528L), cell_fraction = c(1, 1, 1, 1, 1, 0.8, 0.8, 0.8, 1,
  1, 1, 1, 1), gene_cap = c(5000L, 5000L, rep(3000L, 11)), q = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9,
  0.9, 0.9, 0.85, 0.9, 0.9, 0.9), K = c(3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 3L, 3L), nNet = rep(NETWORK_NNET,
  13L), nCells = rep(NETWORK_NCELLS, 13L), force_functional20 = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
  TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE), source_constraint = c(rep("all", 11), "17702045", "GSE196928"),
  comparison_family = c("historical", "historical_optional", rep("seed", 3), rep("cell80", 3), "forced_panel",
    "q", "K", "source", "source"), stringsAsFactors = FALSE)
write_tab(scenarios, "17_135_network_sensitivity_scenario_plan_v02.csv")
status <- data.frame(id = scenarios$id, status = "NOT_STARTED", n_TFs_completed = 0L, info = "", timestamp = "",
  stringsAsFactors = FALSE)
status_file <- file.path(OUT, "18_135_network_scenario_status_v02.csv")
write_status <- function(id, label, msg = "", n = NA_integer_) {
  i <- match(id, status$id)
  status$status[i] <<- label
  if (!is.na(n))
    status$n_TFs_completed[i] <<- n
  status$info[i] <<- substr(as.character(msg), 1, 400)
  status$timestamp[i] <<- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  utils::write.csv(status, status_file, row.names = FALSE, na = "")
}
utils::write.csv(status, status_file, row.names = FALSE, na = "")

scenario_signature <- function(s, genes, cell_ids, version_net, version_knk) {
  c(paste0("id=", s$id), paste0("seed=", s$seed), paste0("cell_fraction=", s$cell_fraction), paste0("q=",
    s$q), paste0("K=", s$K), paste0("nNet=", s$nNet), paste0("nCells=", s$nCells), paste0("force_functional20=",
    s$force_functional20), paste0("source=", s$source_constraint), paste0("n_gene=", length(genes)),
    paste0("n_cell=", length(cell_ids)), paste0("gene_list=", paste(genes, collapse = ";")), paste0("cell_list=",
      paste(cell_ids, collapse = ";")), paste0("scTenifoldNet=", version_net), paste0("scTenifoldKnk=",
      version_knk))
}
write_atomic_gzip <- function(x, path) {
  tmp <- paste0(path, ".incomplete.gz")
  con <- gzfile(tmp, "wt")
  tryCatch(utils::write.csv(x, con, row.names = FALSE, na = ""), finally = close(con))
  if (file.exists(path))
    file.remove(path)
  if (!file.rename(tmp, path))
    stopf("Cannot finalize gzipped CSV: ", path)
}

import_historical <- function(s, subdir) {
  paths <- setNames(file.path(subdir, paste0("KO_", TFS, "_v01.csv.gz")), TFS)
  stat_path <- file.path(subdir, "KO_summary_top10_v01.csv")
  if (!all(file.exists(paths)) || !file.exists(stat_path))
    return(NULL)
  st <- read_tab(stat_path)
  if (!all(TFS %in% st$TF) || nrow(st) != 10L)
    return(NULL)
  drlist <- lapply(paths, read_tab)
  names(drlist) <- TFS
  good <- vapply(TFS, function(tf) {
    d <- drlist[[tf]]
    all(c("gene", "Z", "distance", "p.adj") %in% names(d)) && sum(!duplicated(d$gene)) == nrow(d) &&
      all(is.finite(as_num(d$Z))) && nrow(d) == 5000L
  }, logical(1))
  if (!all(good))
    return(NULL)
  gene_universe <- sort(as.character(drlist[[TFS[1]]]$gene))
  if (!all(vapply(drlist, function(d) identical(sort(as.character(d$gene)), gene_universe), logical(1))))
    return(NULL)
  list(drs = drlist, source = "READ_ONLY v01 5000-gene network checkpoint")
}

if (RUN_NETWORK) {
  needed <- c("cp_counts", "cp_genes", "cp_wt")
  missed <- needed[!file.exists(unlist(P[needed]))]
  if (length(missed))
    stopf("Network checkpoint files missing: ", paste(missed, collapse = ", "))
  for (pkg in c("Matrix", "scTenifoldNet", "scTenifoldKnk", "fgsea")) if (!requireNamespace(pkg, quietly = TRUE))
    stopf("Install required R package: ", pkg)
  if (as.character(utils::packageVersion("scTenifoldNet")) != "1.3" || as.character(utils::packageVersion("scTenifoldKnk")) !=
    "1.0.3")
    stopf("Expected versions scTenifoldNet 1.3, scTenifoldKnk 1.0.3")

  counts <- readRDS(P$cp_counts)
  gene_stats <- read_tab(P$cp_genes)
  require_cols(gene_stats, c("gene", "logCPM_variance"), "Step91 gene selection")
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)))
    stopf("Duplicate gene/cell ID in CP01")
  if (!all(TFS %in% rownames(counts)))
    stopf("Ten TFs not contained in CP01")
  if (length(rownames(counts)) != 5000L)
    stopf("Expected Step91 CP01 5000-gene universe")
  if (anyNA(match(rownames(counts), gene_stats$gene)))
    stopf("Missing logCPM variance annotation")
  old_wt <- readRDS(P$cp_wt)
  if (!setequal(rownames(old_wt), rownames(counts)))
    stopf("Step91 WT genes differ from CP01")
  rm(old_wt)
  gc(verbose = FALSE)
  ix_meta <- match(colnames(counts), meta$cell)
  if (anyNA(ix_meta))
    stopf("CP01 cells missing Step84 metadata")

  functional20 <- c("Clec7a", "Trem2", "Cst7", "Cx3cr1", "Axl", "C1qa", "Apoe", "Lpl", "Abca1", "Apoc1",
    "Apoc4", "Npc2", "Soat1", "Tyrobp", "Csf1r", "Itgam", "Hexb", "Ctsd", "Lamp1", "Lamp2")
  build_genes <- function(force_panel = TRUE, n = 3000L) {
    gs <- gene_stats[match(rownames(counts), gene_stats$gene), , drop = FALSE]
    ord <- order(-as_num(gs$logCPM_variance), as.character(gs$gene))
    ranked <- as.character(gs$gene[ord])
    mandatory <- unique(c(TFS, if (force_panel) intersect(functional20, rownames(counts))))
    if (length(mandatory) > n)
      stopf("Mandatory gene set is larger than cap")
    genes <- c(mandatory, head(setdiff(ranked, mandatory), n - length(mandatory)))
    if (length(genes) != n || anyDuplicated(genes))
      stopf("Gene selection failed")
    genes
  }
  dregfun <- getFromNamespace("dRegulation", "scTenifoldKnk")
  dirfun <- getFromNamespace("strictDirection", "scTenifoldKnk")
  version_net <- as.character(utils::packageVersion("scTenifoldNet"))
  version_knk <- as.character(utils::packageVersion("scTenifoldKnk"))
}

scenario_drs <- list()
run_info <- list()
for (jj in seq_len(nrow(scenarios))) {
  s <- scenarios[jj, , drop = FALSE]
  id <- s$id
  if (!RUN_NETWORK) {
    write_status(id, "SKIPPED_BY_USER", "RUN_NETWORK=FALSE")
    next
  }
  if (s$type != "new") {
    if (s$type == "historical_optional" && !IMPORT_V01_S01_IF_COMPLETE) {
      write_status(id, "OPTIONAL_NOT_REQUESTED", "v01 S01 import disabled")
      next
    }
    folder <- if (s$type == "historical")
      "S00_original_WT91" else "S01_seed13501_full"
    old <- import_historical(s, file.path(V01_NET, folder))
    if (is.null(old)) {
      if (s$type == "historical" && REQUIRE_V01_S00_REFERENCE)
        stopf("Step 135 v01 S00 reference is not complete at ", file.path(V01_NET, folder), ". Save/finish S00 before running v02.")
      write_status(id, "NOT_AVAILABLE", "Historical checkpoint is absent/incomplete; skipped")
      next
    }
    scenario_drs[[id]] <- old$drs
    run_info[[id]] <- data.frame(scenario = id, TF = TFS, n_genes = 5000L, n_source_cells = 12055L,
      seed = s$seed, comparison_family = s$comparison_family, network_provenance = old$source,
      stringsAsFactors = FALSE)
    write_status(id, "COMPLETE_IMPORTED", "Read-only historical 5000-gene KO data", length(TFS))
    next
  }
  if (s$comparison_family == "source" && !RUN_SOURCE_SPECIFIC) {
    write_status(id, "OPTIONAL_NOT_REQUESTED", "RUN_SOURCE_SPECIFIC=FALSE")
    next
  }
  use_cells <- colnames(counts)
  if (s$source_constraint != "all") {
    use_cells <- use_cells[meta$source_dataset[ix_meta] == s$source_constraint]
  }
  if (length(use_cells) < max(500L, NETWORK_NCELLS)) {
    write_status(id, "NOT_ESTIMABLE", paste("Only", length(use_cells), "source cells; require at least",
      NETWORK_NCELLS))
    next
  }
  set.seed(s$seed)
  if (s$cell_fraction < 1) {
    n_sub <- max(NETWORK_NCELLS, floor(length(use_cells) * s$cell_fraction))
    use_cells <- sample(use_cells, n_sub, replace = FALSE)
  }
  genes <- build_genes(force_panel = s$force_functional20, n = s$gene_cap)
  cp <- file.path(NETDIR, id)
  dir.create(cp, recursive = TRUE, showWarnings = FALSE)
  gene_file <- file.path(cp, "00_gene_list_v02.txt")
  cell_file <- file.path(cp, "00_cell_list_v02.txt")
  signature_file <- file.path(cp, "00_parameter_signature_v02.txt")
  WT_file <- file.path(cp, "CP03_WT_network_v02.rds")
  ensemble_file <- file.path(cp, "CP02_WT_ensemble_v02.rds")
  tf_paths <- setNames(file.path(cp, paste0("KO_", TFS, "_v02.csv.gz")), TFS)
  signature <- scenario_signature(s, genes, use_cells, version_net, version_knk)
  existing <- file.exists(signature_file)
  if (existing) {
    sig_old <- readLines(signature_file, warn = FALSE, encoding = "UTF-8")
    if (!identical(sig_old, signature)) {
      stopf("Cached scenario parameters mismatch in ", cp, ". DO NOT reuse; choose a fresh folder/name.")
    }
  } else {
    if (file.exists(WT_file) || file.exists(ensemble_file) || any(file.exists(tf_paths)))
      stopf("Orphaned network checkpoint without signature in ", cp, "; inspect manually")
    writeLines(signature, signature_file, useBytes = TRUE)
    writeLines(genes, gene_file, useBytes = TRUE)
    writeLines(use_cells, cell_file, useBytes = TRUE)
  }
  tryCatch({
    write_status(id, "RUNNING", paste(length(genes), "genes", length(use_cells), "cells"), sum(file.exists(tf_paths)))
    if (file.exists(WT_file)) {
      WT <- as.matrix(readRDS(WT_file))
      if (!identical(rownames(WT), genes) || !identical(colnames(WT), genes))
        stopf("WT gene order does not match scenario signature")
    } else {
      if (file.exists(ensemble_file)) {
        WT_ensemble <- readRDS(ensemble_file)
      } else {
        set.seed(s$seed)
        WT_ensemble <- scTenifoldNet::makeNetworks(X = counts[genes, use_cells, drop = FALSE],
          q = s$q, nNet = s$nNet, nCells = s$nCells, scaleScores = TRUE, symmetric = FALSE, nComp = 3L,
          nCores = min(NETWORK_WORKERS, max(1L, parallel::detectCores() - 1L)))
        tmp <- paste0(ensemble_file, ".incomplete")
        saveRDS(WT_ensemble, tmp, compress = FALSE)
        if (!file.rename(tmp, ensemble_file))
          stopf("Cannot save completed WT ensemble")
      }
      set.seed(s$seed)
      td <- scTenifoldNet::tensorDecomposition(xList = WT_ensemble, K = s$K, maxError = 1e-05,
        maxIter = 1000L, nDecimal = 3L)
      WT <- as.matrix(dirfun(td$X, lambda = 0))
      diag(WT) <- 0
      WT <- t(WT)
      if (!identical(rownames(WT), genes) || !identical(colnames(WT), genes)) {

        if (!setequal(rownames(WT), genes) || !setequal(colnames(WT), genes))
          stopf("Tensor output missing expected genes")
        WT <- WT[genes, genes, drop = FALSE]
      }
      tmp <- paste0(WT_file, ".incomplete")
      saveRDS(Matrix::Matrix(WT, sparse = TRUE), tmp, compress = "gzip")
      if (!file.rename(tmp, WT_file))
        stopf("Cannot save completed WT checkpoint")
      rm(td, WT_ensemble)
      gc(verbose = FALSE)
    }
    if (!setequal(rownames(WT), colnames(WT)))
      stopf("WT not square/gene-labelled")
    if (!all(TFS %in% rownames(WT)))
      stopf("Ten TFs missing from WT")
    drs <- list()
    for (i in seq_along(TFS)) {
      tf <- TFS[i]
      if (file.exists(tf_paths[tf]) && file.info(tf_paths[tf])$size > 0L) {
        d <- read_tab(tf_paths[tf])
      } else {
        if (sum(WT[tf, ] != 0) == 0L)
          stopf(tf, " has no WT outgoing edges")
        KO <- WT
        KO[tf, ] <- 0

        set.seed(9101L + i * 103L)
        MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2L, nCores = min(NETWORK_WORKERS,
          max(1L, parallel::detectCores() - 1L)))
        d <- as.data.frame(dregfun(MA, tf))
        if (!"gene" %in% names(d))
          d$gene <- rownames(d)
        require_cols(d, c("gene", "distance", "Z", "p.adj"), paste0(id, "/", tf))
        d$TF <- tf
        d$scenario <- id
        write_atomic_gzip(d, tf_paths[tf])
        rm(KO, MA)
        gc(verbose = FALSE)
      }
      require_cols(d, c("gene", "distance", "Z", "p.adj"), paste0(id, "/", tf))
      if (nrow(d) != length(genes) || !setequal(as.character(d$gene), genes) || anyDuplicated(d$gene))
        stopf(id, "/", tf, " does not cover expected network gene universe")
      drs[[tf]] <- d
      write_status(id, "RUNNING", paste("KO", i, "/", length(TFS), tf), i)
    }
    scenario_drs[[id]] <- drs
    run_info[[id]] <- data.frame(scenario = id, TF = TFS, n_genes = length(genes), n_source_cells = length(use_cells),
      seed = s$seed, comparison_family = s$comparison_family, network_provenance = "new v02 parameter-verified checkpoint",
      stringsAsFactors = FALSE)
    write_status(id, "COMPLETE", paste(length(genes), "genes", length(use_cells), "cells"), length(TFS))
    rm(WT)
    gc(verbose = FALSE)
  }, error = function(e) {
    write_status(id, "FAILED", conditionMessage(e), sum(file.exists(tf_paths)))
  })
}
if (RUN_NETWORK) {
  rm(counts)
  gc(verbose = FALSE)
}

ko_all <- list()
ranks_all <- list()
for (id in names(scenario_drs)) {
  for (tf in TFS) {
    d <- scenario_drs[[id]][[tf]]
    if (is.null(d))
      next
    d <- d[d$gene != tf & is.finite(as_num(d$Z)) & is.finite(as_num(d$distance)), , drop = FALSE]
    d$Z <- as_num(d$Z)
    d$distance <- as_num(d$distance)
    d$p.adj <- as_num(d$p.adj)
    d <- d[order(-d$Z, d$gene), , drop = FALSE]
    ranks_all[[paste(id, tf, sep = "::")]] <- d[, c("gene", "Z", "distance", "p.adj"), drop = FALSE]
    ko_all[[paste(id, tf, sep = "::")]] <- data.frame(scenario = id, TF = tf, n_evaluable_nonKO = nrow(d),
      n_nonKO_FDR05 = sum(d$p.adj < FDR, na.rm = TRUE), median_top50_Z = finite_med(head(d$Z, TOP_K_STABILITY)),
      n_lipid_genes_in_network = sum(d$gene %in% lipid_gene_union), stringsAsFactors = FALSE)
  }
}
ko_summary <- if (length(ko_all)) do.call(rbind, ko_all) else data.frame()
write_tab(ko_summary, "19_135_all_scenario_top10_virtual_KO_comparison_v02.csv")

reference_ids <- intersect(c("H00_v01_original5000", "B00_seed9101_3000"), names(scenario_drs))
pair_rows <- list()
n_pair <- 0L
for (tf in TFS) {
  for (ref_id in reference_ids) {
    x <- ranks_all[[paste(ref_id, tf, sep = "::")]]
    if (is.null(x))
      next
    for (id in setdiff(names(scenario_drs), ref_id)) {
      y <- ranks_all[[paste(id, tf, sep = "::")]]
      if (is.null(y))
        next
      shared <- intersect(x$gene, y$gene)
      if (length(shared) < 100L)
        next
      xx <- x$Z[match(shared, x$gene)]
      yy <- y$Z[match(shared, y$gene)]
      r_x <- rank(-xx, ties.method = "average")
      r_y <- rank(-yy, ties.method = "average")
      top_x <- shared[order(-xx, shared)][seq_len(min(TOP_K_STABILITY, length(shared)))]
      top_y <- shared[order(-yy, shared)][seq_len(min(TOP_K_STABILITY, length(shared)))]
      n_pair <- n_pair + 1L
      pair_rows[[n_pair]] <- data.frame(TF = tf, reference = ref_id, comparison = id, n_shared_nonKO = length(shared),
        spearman_common_Z = safe_cor(xx, yy), spearman_common_rank = safe_cor(r_x, r_y), topK = length(top_x),
        jaccard_common_universe_topK = length(intersect(top_x, top_y))/length(union(top_x, top_y)),
        topK_intersection = length(intersect(top_x, top_y)), stringsAsFactors = FALSE)
    }
  }
}
pair_df <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame()
write_tab(pair_df, "20_135_common_universe_network_ranking_stability_v02.csv")

fg_rows <- list()
h <- 0L
if (RUN_NETWORK && length(ranks_all)) {
  lipid_sets <- lapply(split(lipid_manifest$gene, lipid_manifest$gene_set), unique)
  for (id in names(scenario_drs)) {
    scenario_result <- list()
    for (i in seq_along(TFS)) {
      tf <- TFS[i]
      d <- ranks_all[[paste(id, tf, sep = "::")]]
      if (is.null(d))
        next
      z <- d$Z
      names(z) <- d$gene
      z <- sort(z, decreasing = TRUE)
      pways <- lapply(lipid_sets, function(gs) intersect(gs, names(z)))
      pways <- pways[vapply(pways, length, integer(1)) >= 5L & vapply(pways, length, integer(1)) <=
        500L]
      if (!length(pways))
        next
      set.seed(9202L + i)
      fg <- tryCatch(fgsea::fgseaMultilevel(pathways = pways, stats = z, minSize = 5L, maxSize = 500L,
        eps = 0), error = function(e) {
        log_msg("FGSEA failed ", id, "/", tf, ": ", conditionMessage(e))
        NULL
      })
      if (is.null(fg) || !nrow(fg))
        next
      fg <- as.data.frame(fg)
      fg$leading_edge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
      fg$leadingEdge <- NULL
      fg$scenario <- id
      fg$TF <- tf
      fg$rank_semantics <- "scTenifoldKnk distance-standardized Z, NOT expression direction"
      scenario_result[[tf]] <- fg
    }
    if (length(scenario_result)) {
      aa <- do.call(rbind, scenario_result)
      aa$FDR_global_scenario <- p.adjust(aa$pval, method = "BH")
      h <- h + 1L
      fg_rows[[h]] <- aa
    }
  }
}
fg_all <- if (length(fg_rows)) do.call(rbind, fg_rows) else data.frame()
write_tab(fg_all, "21_135_top10_lipid_enrichment_all_scenarios_v02.csv")

fg_summary <- list()
h <- 0L
if (nrow(fg_all)) {
  for (tf in TFS) {
    for (gs in sort(unique(fg_all$pathway))) {
      a <- fg_all[fg_all$TF == tf & fg_all$pathway == gs, , drop = FALSE]
      if (!nrow(a))
        next
      h <- h + 1L
      fg_summary[[h]] <- data.frame(TF = tf, lipid_program = gs, n_scenarios_evaluable = nrow(a),
        n_scenarios_FDR05 = sum(a$FDR_global_scenario < 0.05, na.rm = TRUE), fraction_FDR05 = mean(a$FDR_global_scenario <
          0.05, na.rm = TRUE), n_positive_NES = sum(a$NES > 0, na.rm = TRUE), n_negative_NES = sum(a$NES <
          0, na.rm = TRUE), median_NES = finite_med(a$NES), min_NES = if (all(is.na(a$NES)))
          NA_real_ else min(a$NES, na.rm = TRUE), max_NES = if (all(is.na(a$NES)))
          NA_real_ else max(a$NES, na.rm = TRUE), stringsAsFactors = FALSE)
    }
  }
}
fg_sum <- if (length(fg_summary)) do.call(rbind, fg_summary) else data.frame()
write_tab(fg_sum, "22_135_lipid_enrichment_scenario_stability_v02.csv")

gene_baseline_file <- file.path(NETDIR, "B00_seed9101_3000", "00_gene_list_v02.txt")
gene_panel_off_file <- file.path(NETDIR, "D01_noForced20", "00_gene_list_v02.txt")
gene_panel_check <- data.frame(pair = "B00_seed9101_3000 vs D01_noForced20", n_unique_only_baseline = NA_integer_,
  n_unique_only_noForced = NA_integer_, n_changed_total = NA_integer_, status = "NOT_EVALUABLE", stringsAsFactors = FALSE)
if (all(file.exists(c(gene_baseline_file, gene_panel_off_file)))) {
  bg <- readLines(gene_baseline_file, warn = FALSE)
  ng <- readLines(gene_panel_off_file, warn = FALSE)
  n_a <- length(setdiff(bg, ng))
  n_b <- length(setdiff(ng, bg))
  gene_panel_check$n_unique_only_baseline <- n_a
  gene_panel_check$n_unique_only_noForced <- n_b
  gene_panel_check$n_changed_total <- n_a + n_b
  gene_panel_check$status <- if (n_a + n_b > 0)
    "INFORMATIVE" else "UNINFORMATIVE_IDENTICAL_GENE_SET"
}
write_tab(gene_panel_check, "22a_135_forced20_gene_selection_change_audit_v02.csv")
if (gene_panel_check$status == "UNINFORMATIVE_IDENTICAL_GENE_SET") log_msg("WARNING: no-forced-panel and baseline selected identical genes; cannot claim panel sensitivity")

resample_ids <- c("B00_seed9101_3000", "B01_seed13521_3000", "B02_seed13522_3000", "C01_cells80_13523",
  "C02_cells80_13524", "C03_cells80_13525")
resampled <- intersect(resample_ids, names(scenario_drs))
stability <- do.call(rbind, lapply(TFS, function(tf) {
  df <- ko_summary[ko_summary$TF == tf & ko_summary$scenario %in% resampled, , drop = FALSE]
  pair <- pair_df[pair_df$TF == tf & pair_df$reference == "B00_seed9101_3000" & pair_df$comparison %in%
    setdiff(resampled, "B00_seed9101_3000"), , drop = FALSE]
  data.frame(TF = tf, n_resampling_scenarios_completed = nrow(df), n_seed_replicates = sum(df$scenario %in%
    resample_ids[1:3]), n_cell80_replicates = sum(df$scenario %in% resample_ids[4:6]), median_FDR05_nonKO = finite_med(df$n_nonKO_FDR05),
    median_common_universe_Z_spearman = finite_med(pair$spearman_common_Z), median_common_universe_top50_Jaccard = finite_med(pair$jaccard_common_universe_topK),
    stability_claim_ready = (sum(df$scenario %in% resample_ids[1:3]) == 3L && sum(df$scenario %in%
      resample_ids[4:6]) == 3L), stringsAsFactors = FALSE)
}))
write_tab(stability, "23_135_top10_network_resampling_stability_summary_v02.csv")

final <- merge(basic, stability, by = "TF", all.x = TRUE, sort = FALSE)
if (nrow(fg_sum)) {
  fg_tf <- do.call(rbind, lapply(TFS, function(tf) {
    a <- fg_sum[fg_sum$TF == tf, , drop = FALSE]
    data.frame(TF = tf, n_lipid_terms_evaluable = nrow(a), max_lipid_FDR05_scenario_frequency = if (nrow(a))
      max(a$fraction_FDR05, na.rm = TRUE) else NA_real_, stringsAsFactors = FALSE)
  }))
  final <- merge(final, fg_tf, by = "TF", all.x = TRUE, sort = FALSE)
}
final <- final[order(final$RSS_rank), , drop = FALSE]
final$network_stability_complete <- final$stability_claim_ready
write_tab(final, "24_135_top10_unified_five_criterion_comparison_FINAL_v02.csv")

status_counts <- table(status$status)
qc <- data.frame(check = c("all_top10_TFs_present", "RSS_1_to_10", "Step88_AUCell_cell_ID_match", "Step20_lipid_AUCell_cell_ID_match",
  "independent_animals_confirmed", "historical_5000_network_imported", "all_3_seed_networks_complete",
  "all_3_cell80_networks_complete", "fixed_internal_nNet_10", "gene_selection_panel_perturbation_complete",
  "q_parameter_only_complete", "K_parameter_only_complete", "source_specific_two_available", "de_novo_whole_transcriptome_gene_sensitivity",
  "motif_regulon_resampling", "lipid_enrichment_all_10TFs_completed"), status = c(if (setequal(final$TF,
  TFS)) "PASS" else "FAIL", if (identical(sort(as.integer(final$RSS_rank)), 1:10)) "PASS" else "FAIL",
  "PASS", "PASS", "NOT_ESTABLISHED", if ("H00_v01_original5000" %in% names(scenario_drs)) "PASS" else "MISSING",
  if (all(resample_ids[1:3] %in% names(scenario_drs))) "PASS" else "INCOMPLETE", if (all(resample_ids[4:6] %in%
    names(scenario_drs))) "PASS" else "INCOMPLETE", if (NETWORK_NNET == 10L) "PASS" else "CHANGED",
  if ("D01_noForced20" %in% names(scenario_drs) && gene_panel_check$status == "INFORMATIVE") "PASS" else "INCOMPLETE_OR_UNINFORMATIVE",
  if ("E01_q085" %in% names(scenario_drs)) "PASS" else "INCOMPLETE", if ("E02_K4" %in% names(scenario_drs)) "PASS" else "INCOMPLETE",
  if (all(c("F01_source_17702045", "F02_source_GSE196928") %in% names(scenario_drs))) "PASS" else "LIMITED_OR_NOT_ESTIMABLE",
  "NOT_DONE_FROZEN_5000_POOL", "NOT_DONE", if (nrow(fg_all) > 0 && all(ko_summary$scenario %in% fg_all$scenario) &&
    all(TFS %in% fg_all$TF)) "PASS" else "INCOMPLETE"), detail = c("Same top ten RSS TFs for every analysis",
  "Step90 order", "Cell IDs verified in v02 descriptive phase", "Cell IDs verified in v02 descriptive phase",
  "Sample labels may not represent independent animals", "Step135 v01 S00 full non-KO results are required",
  "3 independent 3000-gene WT random seeds", "3 distinct 80% cell subsets", "Network nNet is fixed at the Step91 value",
  "No forced 20 genes, 10 TFs always retained", "Paired with B00: only q changes", "Paired with B00: only tensor K changes",
  "Depends on post-QC cell number; source and injury stages confounded", "3000 selected only within Step91's 5000 genes",
  "scTenifoldNet network sensitivity does not rerun pySCENIC motif pruning", "Every completed scenario must have lipid FGSEA outputs for all ten TFs"),
  stringsAsFactors = FALSE)
write_tab(qc, "25_135_FINAL_QC_v02.csv")
write_tab(data.frame(item = c("R", "root", "new_output_dir", "network_workers", "nNet", "nCells", "n_replicate_seeds",
  "n_cell80_subsamples", "run_source_specific", "scTenifoldNet", "scTenifoldKnk", "n_computational_scenarios_completed",
  "timestamp"), value = c(R.version.string, ROOT, OUT, NETWORK_WORKERS, NETWORK_NNET, NETWORK_NCELLS,
  3, 3, as.character(RUN_SOURCE_SPECIFIC), if (requireNamespace("scTenifoldNet", quietly = TRUE)) as.character(utils::packageVersion("scTenifoldNet")) else "missing",
  if (requireNamespace("scTenifoldKnk", quietly = TRUE)) as.character(utils::packageVersion("scTenifoldKnk")) else "missing",
  length(scenario_drs), format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "26_135_environment_parameters_v02.csv")
write_text("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    "27_135_REPORT_v02.txt")
files_out <- list.files(OUT, recursive = TRUE, full.names = TRUE)
files_out <- files_out[!is.na(file.info(files_out)$isdir) & !file.info(files_out)$isdir]
rel <- substring(files_out, nchar(OUT) + 2L)
write_tab(data.frame(relative_path = gsub("\\\\", "/", rel), size_bytes = unname(file.info(files_out)$size)),
  "00_135_OUTPUT_INDEX_v02.csv")
