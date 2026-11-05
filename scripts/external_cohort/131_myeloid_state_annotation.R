options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
ROOT <- "D:/zz"
STEM <- "GSE281205"
D124 <- file.path(ROOT, "results", paste0("124_", STEM, "_lineage_review_myeloid_extraction"))
D128 <- file.path(ROOT, "results", paste0("128_", STEM, "_reviewed_lineage_finalization"))
D130 <- file.path(ROOT, "results", paste0("130_", STEM, "_myeloid_reclustering"))
OUT <- file.path(ROOT, "results", paste0("131_", STEM, "_myeloid_subtype_annotation_and_microglia_groups"))
F_DIR <- file.path(OUT, "figures")
R_DIR <- file.path(OUT, "RDS")
L_DIR <- file.path(OUT, "logs")

IN_128_IDS <- file.path(D128, "06_128_myeloid_include_cell_ids_v02.csv")
IN_130_MEM <- file.path(D130, "05_130_per_cell_myeloid_cluster_membership_v02.csv")
IN_130_MARKERS <- file.path(D130, "10_130_canonical_marker_evidence_v02.csv")
IN_130_RDS <- file.path(D130, "RDS", "130_GSE281205_myeloid_LogNorm_PCA_Harmony_clusters_v02.rds")
IN_FEATURES <- file.path(D124, "12_124_gene_feature_mapping_v01.csv")
OUT_MARKER_CACHE <- file.path(R_DIR, "131_per_cell_marker_flags_cache_v01.rds")
OUT_MANIFEST_RDS <- file.path(R_DIR, "131_analysis_cell_manifest_v01.rds")
LOG_FILE <- file.path(L_DIR, "00_131_run_log_v01.log")
EXPECTED_GSM <- paste0("GSM", 8614662:8614667)
EXPECTED_N <- 28716L
PRIMARY_CONTAM_LT <- 20
SENS_CONTAM_LT <- 50
MIN_MICROGLIA_CORE <- 2L
MIN_PB_CELLS <- 20L
MIN_REPS_PER_STAGE <- 3L
RESOLUTION_PURITY_FLAG <- 0.8
SEED <- 13101L
set.seed(SEED)

pkg <- c("SeuratObject", "Matrix")
missing_pkgs <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
for (p in c(OUT, F_DIR, R_DIR, L_DIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

log_it <- function(...) {
  z <- paste0(..., collapse = "")
  cat(z, "\n", sep = "")
  cat(z, "\n", sep = "", file = LOG_FILE, append = TRUE)
  invisible(z)
}
assert <- function(condition, message) {
  if (!isTRUE(condition))
    stop(message, call. = FALSE)
}
read_csv <- function(f) {
  assert(file.exists(f), paste("Required input missing:", f))
  utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}
write_csv <- function(d, filename) {
  assert(is.data.frame(d), paste("Expected data.frame:", filename))
  utils::write.csv(d, file.path(OUT, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  log_it("SAVED ", filename, " [", nrow(d), " rows]")
}
require_cols <- function(d, columns, label) {
  absent <- setdiff(columns, names(d))
  assert(!length(absent), paste(label, "missing required columns:", paste(absent, collapse = ", ")))
}
int_dpi <- function(x) {
  z <- suppressWarnings(as.integer(as.character(x)))
  assert(!anyNA(z) && all(z %in% c(14L, 28L)), "dpi must be 14 or 28")
  z
}
median_or_na <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x)))
    return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

partition_review <- function(d, field, label) {
  tab <- as.data.frame(table(r010 = factor(as.character(d$myeloid_harmony_r010), levels = as.character(0:4)),
    refined = as.character(d[[field]])), stringsAsFactors = FALSE)
  names(tab)[3] <- "n_cells"
  tab$analysis <- label
  denom_refined <- tapply(tab$n_cells, tab$refined, sum)
  tab$fraction_within_refined <- ifelse(as.numeric(denom_refined[tab$refined]) > 0, tab$n_cells/as.numeric(denom_refined[tab$refined]),
    NA_real_)
  domin <- tapply(tab$fraction_within_refined, tab$refined, max)
  tab$refined_cluster_parent_purity <- as.numeric(domin[tab$refined])
  tab[, c("analysis", "r010", "refined", "n_cells", "fraction_within_refined", "refined_cluster_parent_purity")]
}

IDENTITY_MARKERS <- c("Fcrls", "Hexb", "Cx3cr1")
STATE_MARKERS <- c("P2ry12", "Tmem119", "Clec7a", "Apoe", "Spp1", "Lpl", "Gpnmb", "Cst7", "Ccr2", "Ly6c2",
  "Ms4a7", "Flt3", "Cd209a", "Mki67", "Pclaf", "Birc5", "Cd74", "H2-Aa")
MARKERS <- unique(c(IDENTITY_MARKERS, STATE_MARKERS))

get_marker_flags <- function(membership, features) {
  desired_ids <- membership$cell_id
  rds_info <- file.info(IN_130_RDS)
  rds_size <- as.numeric(rds_info$size)
  rds_mtime <- as.numeric(rds_info$mtime)
  mapping <- do.call(rbind, lapply(MARKERS, function(symbol) {
    ids <- as.character(features$gene_id[as.character(features$gene_symbol) == symbol])
    data.frame(gene_symbol = symbol, gene_id = if (length(ids))
      ids[1] else NA_character_, n_mapping_hits = length(ids), stringsAsFactors = FALSE)
  }))
  assert(all(mapping$n_mapping_hits == 1L), paste("Core marker gene mapping missing or ambiguous:",
    paste(mapping$gene_symbol[mapping$n_mapping_hits != 1L], collapse = ", ")))
  if (file.exists(OUT_MARKER_CACHE)) {
    cached <- tryCatch(readRDS(OUT_MARKER_CACHE), error = function(e) NULL)
    if (is.list(cached) && identical(cached$version, "v01") && identical(cached$cell_id, desired_ids) &&
      identical(cached$marker_gene_ids, mapping$gene_id) && identical(cached$source_RDS_size, rds_size) &&
      identical(cached$source_RDS_mtime, rds_mtime) && is.data.frame(cached$flags)) {
      log_it("RESUME: validated existing per-cell marker flag cache")
      return(list(flags = cached$flags, mapping = mapping))
    }
  }
  log_it("Loading Step130 clustered Seurat object for cell-level marker evidence: ", IN_130_RDS)
  seu <- readRDS(IN_130_RDS)
  assert(inherits(seu, "Seurat"), "Step130 RDS is not a Seurat object")
  assert(ncol(seu) == EXPECTED_N, "Unexpected Seurat cell count in Step130")
  assert(setequal(colnames(seu), desired_ids), "Step130 RDS cell IDs do not match Step130 CSV")
  rds_meta <- seu[[]]
  mi <- match(desired_ids, rownames(rds_meta))
  assert(!anyNA(mi), "Step130 RDS metadata is missing some cell IDs")
  assert("myeloid_harmony_r010" %in% colnames(rds_meta) && identical(as.character(rds_meta$myeloid_harmony_r010[mi]),
    as.character(membership$myeloid_harmony_r010)), "Step130 RDS and CSV contain DIFFERENT r0.10 cluster assignments")
  rm(rds_meta, mi)
  assert("DecontX" %in% SeuratObject::Assays(seu), "Step130 Seurat RDS has no DecontX assay")
  assay <- seu[["DecontX"]]
  count_mat <- SeuratObject::LayerData(assay, layer = "counts")
  assert(all(mapping$gene_id %in% rownames(count_mat)), "Some mapped marker features are absent from the Step130 DecontX assay")
  ci <- match(desired_ids, colnames(count_mat))
  assert(!anyNA(ci), "DecontX assay count matrix cannot be aligned to Step130 CSV")
  flags <- data.frame(cell_id = desired_ids, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(mapping))) {
    symbol <- mapping$gene_symbol[i]
    x <- as.numeric(count_mat[mapping$gene_id[i], ci, drop = TRUE])
    flags[[paste0("marker_", symbol, "_ge1")]] <- is.finite(x) & x >= 1
  }
  rm(count_mat, assay, seu)
  gc(verbose = FALSE)
  saveRDS(list(version = "v01", cell_id = desired_ids, marker_gene_ids = mapping$gene_id, source_RDS_size = rds_size,
    source_RDS_mtime = rds_mtime, flags = flags), OUT_MARKER_CACHE, compress = TRUE)
  log_it("SAVED compact marker flag cache: ", OUT_MARKER_CACHE)
  list(flags = flags, mapping = mapping)
}

run <- function() {
  t0 <- Sys.time()
  log_it("=====================================================")
  log_it("STEP131 v01 | High-confidence lineage/state groups | ", format(t0, "%Y-%m-%d %H:%M:%S"))
  log_it("STRICT no ZIP; no signature access; no cell-level DE; no pseudobulk here")
  input_paths <- c(IN_128_IDS, IN_130_MEM, IN_130_MARKERS, IN_130_RDS, IN_FEATURES)
  assert(all(file.exists(input_paths)), paste("Missing files:", paste(input_paths[!file.exists(input_paths)],
    collapse = "; ")))
  write_csv(data.frame(input_role = c("Step128_myeloid_INCLUDE_IDs", "Step130_recluster_membership",
    "Step130_cluster_marker_evidence", "Step130_clustered_Seurat_DecontX", "Step124_gene_mapping"),
    input_path = input_paths, n_bytes = as.numeric(file.info(input_paths)$size), stringsAsFactors = FALSE),
    "01_131_input_provenance_v01.csv")

  inc <- read_csv(IN_128_IDS)
  d <- read_csv(IN_130_MEM)
  marker_summary <- read_csv(IN_130_MARKERS)
  features <- read_csv(IN_FEATURES)
  require_cols(inc, c("cell_id", "gsm", "dpi", "curated_major_lineage_128"), "Step128")
  require_cols(d, c("cell_id", "gsm", "dpi", "replicate", "parent_cluster_123", "curated_major_lineage_128",
    "myeloid_harmony_r010", "myeloid_harmony_r015", "myeloid_harmony_r020", "myeloid_pca_r010", "contamination_percent"),
    "Step130")
  require_cols(features, c("gene_id", "gene_symbol"), "Step124 feature mapping")
  require_cols(marker_summary, c("cluster", "gene_symbol", "pct_cells_corrected_counts_ge1"), "Step130 marker evidence")
  assert(nrow(d) == EXPECTED_N && nrow(inc) == EXPECTED_N, "Expected 28,716 candidate myeloid cells")
  assert(!anyNA(d$cell_id) && !anyDuplicated(d$cell_id) && !anyNA(inc$cell_id) && !anyDuplicated(inc$cell_id) &&
    setequal(d$cell_id, inc$cell_id), "Step128/130 cell IDs disagree or duplicate")
  ix <- match(d$cell_id, inc$cell_id)
  assert(!anyNA(ix) && identical(as.character(d$gsm), as.character(inc$gsm[ix])) && identical(int_dpi(d$dpi),
    int_dpi(inc$dpi[ix])) && identical(as.character(d$curated_major_lineage_128), as.character(inc$curated_major_lineage_128[ix])),
    "Step128 vs Step130 per-cell GSM/dpi/lineage mismatch")
  assert(setequal(d$gsm, EXPECTED_GSM), "GSM IDs differ from six validated independent samples")
  assert(all(is.finite(d$contamination_percent) & d$contamination_percent >= 0 & d$contamination_percent <=
    100), "contamination_percent must be within 0..100")
  assert(setequal(as.character(d$myeloid_harmony_r010), as.character(0:4)), "Expected r0.10 clusters 0,1,2,3,4; cannot silently relabel another run")
  actual_n <- table(as.character(d$myeloid_harmony_r010))
  expected_n <- c(`0` = 13236L, `1` = 13110L, `2` = 1926L, `3` = 317L, `4` = 127L)
  assert(identical(as.integer(actual_n[names(expected_n)]), as.integer(expected_n)), "Step130 r0.10 cell counts changed; review before using fixed label map")
  assert(!anyNA(d$myeloid_harmony_r015) && !anyNA(d$myeloid_harmony_r020), "Missing Harmony r0.15/r0.20 clusters")
  assert(setequal(unique(d$curated_major_lineage_128), c("Microglia", "Monocyte_macrophage", "APC_dendritic")),
    "Step128 lineage categories unexpectedly changed")
  d$dpi <- int_dpi(d$dpi)
  for (g in EXPECTED_GSM) {
    expected_dpi <- if (g %in% EXPECTED_GSM[1:3])
      14L else 28L
    assert(all(d$dpi[d$gsm == g] == expected_dpi), paste("GSM vs stage mismatch:", g))
  }
  d$myeloid_harmony_r010 <- as.character(d$myeloid_harmony_r010)
  d$myeloid_harmony_r015 <- as.character(d$myeloid_harmony_r015)
  d$myeloid_harmony_r020 <- as.character(d$myeloid_harmony_r020)
  d$myeloid_pca_r010 <- as.character(d$myeloid_pca_r010)
  d$stage_label <- paste0(d$dpi, "dpi")

  state_map <- data.frame(myeloid_harmony_r010 = as.character(0:4), provisional_transcriptomic_state = c("P2ry12_high_microglia_like",
    "Clec7a_high_lipid_associated_microglia_like", "Monocyte_macrophage_like", "APC_dendritic_like",
    "Proliferating_microglia_like"), defining_evidence = c("P2ry12/Tmem119/Fcrls/Hexb/Cx3cr1", "Clec7a/Spp1/Lpl/Gpnmb/Cst7 + Fcrls/Hexb/Cx3cr1",
    "Ms4a7/Lyz2/Mrc1 with macrophage/APC genes", "Cd74/H2-Aa/Flt3/Cd209a, heterogeneous quality",
    "Mki67/Pclaf/Birc5 with microglial identity markers"), reference_lineage_128 = c("Microglia",
    "Microglia", "Monocyte_macrophage", "APC_dendritic", "Microglia"), notes = c("Not Clec7a-negative: Clec7a is detected in many cells",
    "Working activation state, not developmental ontogeny or 1320 signature", "Do not merge all this cluster with resident microglia",
    "Includes discordant parent labels and potential cycling/high-contamination cells", "Separate state; broad sensitivity can include with other microglia"),
    stringsAsFactors = FALSE)
  state_map$n_cluster_cells <- as.integer(actual_n[state_map$myeloid_harmony_r010])
  write_csv(state_map, "02_131_myeloid_subcluster_working_annotations_v01.csv")

  refinements <- rbind(partition_review(d, "myeloid_harmony_r015", "Harmony_r0.15"), partition_review(d,
    "myeloid_harmony_r020", "Harmony_r0.20"), partition_review(d, "myeloid_pca_r010", "PCA_uncorrected_r0.10"))
  write_csv(refinements, "07_131_cross_resolution_cluster_contingencies_v01.csv")
  for (pair in list(c("myeloid_harmony_r015", "Harmony_r0.15", "r015"), c("myeloid_harmony_r020", "Harmony_r0.20",
    "r020"), c("myeloid_pca_r010", "PCA_uncorrected_r0.10", "pca"))) {
    x <- refinements[refinements$analysis == pair[2], , drop = FALSE]
    q <- tapply(x$refined_cluster_parent_purity, x$refined, max)
    d[[paste0(pair[3], "_refined_parent_purity")]] <- as.numeric(q[d[[pair[1]]]])
  }
  d$refined_partition_caution <- d$r015_refined_parent_purity < RESOLUTION_PURITY_FLAG | d$r020_refined_parent_purity <
    RESOLUTION_PURITY_FLAG | d$pca_refined_parent_purity < RESOLUTION_PURITY_FLAG

  mapping_result <- get_marker_flags(d, features)
  f <- mapping_result$flags
  assert(identical(f$cell_id, d$cell_id), "Marker flags not aligned to manifest")
  write_csv(mapping_result$mapping, "03_131_predefined_marker_feature_mapping_v01.csv")
  d$microglia_core_markers_ge1 <- as.integer(f$marker_Fcrls_ge1) + as.integer(f$marker_Hexb_ge1) +
    as.integer(f$marker_Cx3cr1_ge1)
  d$microglia_core_supported <- d$microglia_core_markers_ge1 >= MIN_MICROGLIA_CORE
  d$monocyte_pair_Ccr2_Ly6c2 <- f$marker_Ccr2_ge1 & f$marker_Ly6c2_ge1
  d$microglia_identity_conflict <- d$monocyte_pair_Ccr2_Ly6c2
  d$has_Clec7a_ge1 <- f$marker_Clec7a_ge1
  d$has_P2ry12_ge1 <- f$marker_P2ry12_ge1
  d$has_Tmem119_ge1 <- f$marker_Tmem119_ge1
  d$lipid_state_markers_ge1 <- as.integer(f$marker_Spp1_ge1) + as.integer(f$marker_Lpl_ge1) + as.integer(f$marker_Gpnmb_ge1) +
    as.integer(f$marker_Cst7_ge1)
  d$cycling_markers_ge1 <- as.integer(f$marker_Mki67_ge1) + as.integer(f$marker_Pclaf_ge1) + as.integer(f$marker_Birc5_ge1)
  d$qc_lt20 <- d$contamination_percent < PRIMARY_CONTAM_LT
  d$qc_lt50 <- d$contamination_percent < SENS_CONTAM_LT
  parent <- as.character(d$curated_major_lineage_128)
  cl <- d$myeloid_harmony_r010
  d$parent_cluster_concordant <- parent == state_map$reference_lineage_128[match(cl, state_map$myeloid_harmony_r010)]
  d$annotation_131 <- ifelse(cl == "0", state_map$provisional_transcriptomic_state[1], ifelse(cl ==
    "1", state_map$provisional_transcriptomic_state[2], ifelse(cl == "2", state_map$provisional_transcriptomic_state[3],
    ifelse(cl == "3", state_map$provisional_transcriptomic_state[4], state_map$provisional_transcriptomic_state[5]))))
  d$cross_lineage_discordance_note <- ifelse(d$parent_cluster_concordant, "CONCORDANT_PARENT_CLUSTER",
    "PARENT_R010_DISCORDANT_INSPECT")
  is_mg <- parent == "Microglia"
  is_mono <- parent == "Monocyte_macrophage"
  is_apc <- parent == "APC_dendritic"
  mg_confident <- is_mg & d$microglia_core_supported & !d$microglia_identity_conflict

  eligible_mg <- mg_confident & d$qc_lt20
  eligible_nonmg <- d$qc_lt20
  d$primary_group_131 <- "EXCLUDED_AMBIGUOUS_OR_QC"
  d$primary_group_131[eligible_mg & cl == "1"] <- "TARGET_CLEC7A_HIGH_MICROGLIA"
  d$primary_group_131[eligible_mg & cl == "0"] <- "REFERENCE_P2RY12_HIGH_MICROGLIA"
  d$primary_group_131[eligible_mg & cl == "4"] <- "PROLIFERATING_MICROGLIA"
  d$primary_group_131[eligible_nonmg & is_mono & cl == "2"] <- "OTHER_MONOCYTE_MACROPHAGE"
  d$primary_group_131[eligible_nonmg & is_apc & cl == "3"] <- "OTHER_APC_DC_LIKE"
  d$included_primary <- d$primary_group_131 != "EXCLUDED_AMBIGUOUS_OR_QC"
  d$primary_exclusion_reason <- ""
  ex <- !d$included_primary
  d$primary_exclusion_reason[ex & !d$qc_lt20] <- "DecontX_contamination_ge20pct"
  d$primary_exclusion_reason[ex & d$qc_lt20 & is_mg & !d$microglia_core_supported] <- "microglia_core_markers_lt2of3"
  d$primary_exclusion_reason[ex & d$qc_lt20 & is_mg & d$microglia_identity_conflict] <- "monocyte_Ccr2_Ly6c2_coexpression"
  d$primary_exclusion_reason[ex & d$qc_lt20 & !d$parent_cluster_concordant] <- "Step128_vs_Step130_lineage_discordance"
  d$primary_exclusion_reason[ex & d$qc_lt20 & d$primary_exclusion_reason == ""] <- "noncanonical_microglia_state_or_other_unresolved"
  assert(all(d$primary_group_131 %in% c("EXCLUDED_AMBIGUOUS_OR_QC", "TARGET_CLEC7A_HIGH_MICROGLIA",
    "REFERENCE_P2RY12_HIGH_MICROGLIA", "PROLIFERATING_MICROGLIA", "OTHER_MONOCYTE_MACROPHAGE", "OTHER_APC_DC_LIKE")),
    "Unexpected primary analysis group")

  target <- d$primary_group_131 == "TARGET_CLEC7A_HIGH_MICROGLIA"
  reference <- d$primary_group_131 == "REFERENCE_P2RY12_HIGH_MICROGLIA"
  cycling <- d$primary_group_131 == "PROLIFERATING_MICROGLIA"
  nonmg <- d$primary_group_131 %in% c("OTHER_MONOCYTE_MACROPHAGE", "OTHER_APC_DC_LIKE")
  assign_contrast <- function(t, r) ifelse(t, "TARGET", ifelse(r, "OTHER", "NOT_IN_CONTRAST"))
  d$contrast_A_target_vs_other_microglia_131 <- assign_contrast(target, reference)
  d$contrast_B_target_vs_all_other_highconf_myeloid_131 <- assign_contrast(target, reference | cycling |
    nonmg)
  d$contrast_C_target_vs_nonmicroglial_myeloid_131 <- assign_contrast(target, nonmg)

  sens_target <- is_mg & cl == "1" & d$qc_lt50
  sens_other_mg <- is_mg & cl %in% c("0", "4") & d$qc_lt50
  sens_nonmg <- ((is_mono & cl == "2") | (is_apc & cl == "3")) & d$qc_lt50
  d$sensitivity_A_broad_microglia_lt50 <- assign_contrast(sens_target, sens_other_mg)
  d$sensitivity_B_broad_all_myeloid_lt50 <- assign_contrast(sens_target, sens_other_mg | sens_nonmg)

  d$sensitivity_A_partition_consistent <- assign_contrast(target & !d$refined_partition_caution, reference &
    !d$refined_partition_caution)

  flag_names <- c("marker_Fcrls_ge1", "marker_Hexb_ge1", "marker_Cx3cr1_ge1", "marker_Clec7a_ge1",
    "marker_P2ry12_ge1", "marker_Tmem119_ge1", "marker_Spp1_ge1", "marker_Lpl_ge1", "marker_Gpnmb_ge1",
    "marker_Cst7_ge1", "marker_Ccr2_ge1", "marker_Ly6c2_ge1")
  for (nm in flag_names) d[[nm]] <- f[[nm]]
  assert(nrow(d) == EXPECTED_N && !anyDuplicated(d$cell_id) && sum(d$included_primary) <= EXPECTED_N,
    "Post-gate manifest integrity failure")
  write_csv(d, "04_131_per_cell_subtype_and_analysis_manifest_v01.csv")
  saveRDS(d, OUT_MANIFEST_RDS, compress = TRUE)
  log_it("SAVED compact RDS manifest: ", OUT_MANIFEST_RDS)

  clabs <- c("TARGET_CLEC7A_HIGH_MICROGLIA", "REFERENCE_P2RY12_HIGH_MICROGLIA", "PROLIFERATING_MICROGLIA",
    "OTHER_MONOCYTE_MACROPHAGE", "OTHER_APC_DC_LIKE", "EXCLUDED_AMBIGUOUS_OR_QC")
  per_sample <- as.data.frame(table(gsm = factor(d$gsm, levels = EXPECTED_GSM), group_131 = factor(d$primary_group_131,
    levels = clabs)), stringsAsFactors = FALSE)
  names(per_sample)[3] <- "n_cells"
  per_sample$dpi <- ifelse(per_sample$gsm %in% EXPECTED_GSM[1:3], 14L, 28L)
  gsm_totals <- table(factor(d$gsm, levels = EXPECTED_GSM))
  per_sample$total_candidate_myeloid_in_GSM <- as.integer(gsm_totals[as.character(per_sample$gsm)])
  per_sample$fraction_of_candidate_myeloid_in_GSM <- per_sample$n_cells/per_sample$total_candidate_myeloid_in_GSM
  write_csv(per_sample, "05_131_GSM_primary_group_composition_v01.csv")

  stage_summary <- do.call(rbind, lapply(sort(unique(per_sample$group_131)), function(gr) {
    do.call(rbind, lapply(c(14L, 28L), function(st) {
      sub <- per_sample[per_sample$group_131 == gr & per_sample$dpi == st, , drop = FALSE]
      data.frame(group_131 = gr, dpi = st, n_libraries = nrow(sub), n_cells_pooled = sum(sub$n_cells),
        median_n_per_library = median_or_na(sub$n_cells), mean_fraction_per_library = mean(sub$fraction_of_candidate_myeloid_in_GSM),
        range_fraction_min = min(sub$fraction_of_candidate_myeloid_in_GSM), range_fraction_max = max(sub$fraction_of_candidate_myeloid_in_GSM),
        note = "Descriptive specimen-level summary; no biological replicate inference", stringsAsFactors = FALSE)
    }))
  }))
  write_csv(stage_summary, "06_131_stage_group_composition_descriptive_v01.csv")

  cross <- as.data.frame(table(parent_128 = as.character(d$curated_major_lineage_128), cluster_130_r010 = d$myeloid_harmony_r010),
    stringsAsFactors = FALSE)
  names(cross)[3] <- "n_cells"
  write_csv(cross, "08_131_parent_lineage_vs_myeloid_state_v01.csv")

  exc <- as.data.frame(table(gsm = factor(d$gsm, levels = EXPECTED_GSM), primary_exclusion_reason = ifelse(d$included_primary,
    "INCLUDED", d$primary_exclusion_reason)), stringsAsFactors = FALSE)
  names(exc)[3] <- "n_cells"
  exc$dpi <- ifelse(exc$gsm %in% EXPECTED_GSM[1:3], 14L, 28L)
  write_csv(exc, "09_131_exclusion_and_QC_audit_by_GSM_v01.csv")

  csum <- do.call(rbind, lapply(as.character(0:4), function(k) {
    x <- d[d$myeloid_harmony_r010 == k, , drop = FALSE]
    data.frame(cluster = k, n_cells = nrow(x), n_parent_concordant = sum(x$parent_cluster_concordant),
      n_primary_included = sum(x$included_primary), pct_parent_concordant = 100 * mean(x$parent_cluster_concordant),
      pct_core_markers_ge2 = 100 * mean(x$microglia_core_supported), pct_Clec7a_detected_ge1 = 100 *
        mean(x$has_Clec7a_ge1), pct_P2ry12_detected_ge1 = 100 * mean(x$has_P2ry12_ge1), pct_Tmem119_detected_ge1 = 100 *
        mean(x$has_Tmem119_ge1), pct_DecontX_ge20 = 100 * mean(!x$qc_lt20), pct_DecontX_ge50 = 100 *
        mean(!x$qc_lt50), pct_partition_caution = 100 * mean(x$refined_partition_caution), n_14dpi = sum(x$dpi ==
        14L), n_28dpi = sum(x$dpi == 28L), stringsAsFactors = FALSE)
  }))
  write_csv(csum, "10_131_subcluster_identity_QC_stability_summary_v01.csv")

  contrast_fields <- c("contrast_A_target_vs_other_microglia_131", "contrast_B_target_vs_all_other_highconf_myeloid_131",
    "contrast_C_target_vs_nonmicroglial_myeloid_131", "sensitivity_A_broad_microglia_lt50", "sensitivity_B_broad_all_myeloid_lt50",
    "sensitivity_A_partition_consistent")
  all_contrast_rows <- list()
  feasibility_rows <- list()
  for (field in contrast_fields) {
    tab <- as.data.frame(table(gsm = factor(d$gsm, levels = EXPECTED_GSM), group = factor(d[[field]],
      levels = c("TARGET", "OTHER"))), stringsAsFactors = FALSE)
    names(tab)[3] <- "n_cells"
    tab$comparison <- field
    tab$dpi <- ifelse(tab$gsm %in% EXPECTED_GSM[1:3], 14L, 28L)
    tab$specimen_meets_min_cells <- tab$n_cells >= MIN_PB_CELLS
    all_contrast_rows[[field]] <- tab[, c("comparison", "gsm", "dpi", "group", "n_cells", "specimen_meets_min_cells")]
    for (stage in c(14L, 28L)) {
      for (grp in c("TARGET", "OTHER")) {
        v <- tab[tab$dpi == stage & tab$group == grp, , drop = FALSE]
        feasibility_rows[[length(feasibility_rows) + 1L]] <- data.frame(comparison = field, dpi = stage,
          group = grp, n_GSM_available = nrow(v), n_GSM_above_min_cells = sum(v$specimen_meets_min_cells),
          min_cells_across_GSM = min(v$n_cells), max_cells_across_GSM = max(v$n_cells), ready_for_pseudobulk = nrow(v) >=
          MIN_REPS_PER_STAGE && sum(v$specimen_meets_min_cells) >= MIN_REPS_PER_STAGE, notes = "Feasibility only. Step132 will aggregate original Step119 integer UMIs.",
          stringsAsFactors = FALSE)
      }
    }
  }
  contrast_counts <- do.call(rbind, all_contrast_rows)
  feas <- do.call(rbind, feasibility_rows)
  write_csv(contrast_counts, "11_131_GSM_contrast_group_counts_v01.csv")
  write_csv(feas, "12_131_stage_contrast_pseudobulk_feasibility_v01.csv")

  output_cols <- c("cell_id", "gsm", "dpi", "replicate", "parent_cluster_123", "curated_major_lineage_128",
    "myeloid_harmony_r010", "myeloid_harmony_r015", "myeloid_harmony_r020", "myeloid_pca_r010", "primary_group_131",
    "microglia_core_markers_ge1", "microglia_core_supported", "contamination_percent", "refined_partition_caution",
    contrast_fields)
  write_csv(d[, output_cols, drop = FALSE], "13_131_READY_FOR_STEP132_original_UMI_cell_manifest_v01.csv")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p1 <- ggplot2::ggplot(per_sample[per_sample$group_131 != "EXCLUDED_AMBIGUOUS_OR_QC", ], ggplot2::aes(x = gsm,
      y = fraction_of_candidate_myeloid_in_GSM, fill = group_131)) + ggplot2::geom_col(width = 0.82) +
      ggplot2::facet_grid(~dpi, scales = "free_x", space = "free_x") + ggplot2::labs(title = "Specimen-level composition of curated high-confidence myeloid groups",
      subtitle = "Fractions denominator: all 28,716 candidate-myeloid cells by GSM; excluded cells retain their denominator",
      x = "Independent GSM library", y = "Fraction of candidate-myeloid cells", fill = "Working group") +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35,
      hjust = 1))
    tryCatch(ggplot2::ggsave(file.path(F_DIR, "Fig_131_curated_groups_per_GSM_v01.png"), plot = p1,
      width = 13, height = 6, dpi = 150), error = function(e) log_it("Figure warning: ", conditionMessage(e)))
    p2 <- ggplot2::ggplot(csum, ggplot2::aes(x = cluster, y = pct_core_markers_ge2)) + ggplot2::geom_col() +
      ggplot2::labs(title = "Predefined microglial lineage marker support by myeloid cluster",
        subtitle = "At least 2 of Fcrls/Hexb/Cx3cr1 (estimated DecontX count >=1)", x = "Harmony r0.10 cluster",
        y = "% cells passing core-marker rule") + ggplot2::theme_bw(base_size = 11)
    tryCatch(ggplot2::ggsave(file.path(F_DIR, "Fig_131_microglia_marker_support_v01.png"), plot = p2,
      width = 8, height = 5, dpi = 150), error = function(e) log_it("Figure warning: ", conditionMessage(e)))
  } else {
    log_it("ggplot2 unavailable: CSV/RDS/TXT results are unaffected")
  }

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, "14_131_REPORT_v01.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(OUT, "15_131_sessionInfo_v01.txt"))
  output_files <- list.files(OUT, recursive = TRUE, full.names = TRUE)
  output_files <- output_files[file.exists(output_files) & !dir.exists(output_files)]
  write_csv(data.frame(relative_path = substring(output_files, nchar(OUT) + 2L), bytes = as.numeric(file.info(output_files)$size),
    stringsAsFactors = FALSE), "16_131_OUTPUT_INDEX_v01.csv")
  log_it("SUCCESS | Step131 v01 | ", sum(d$included_primary), " high-confidence myeloid cells classified; ",
    sum(d$primary_group_131 == "TARGET_CLEC7A_HIGH_MICROGLIA"), " target cells; 6 GSM checked")
  log_it("Elapsed seconds: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  invisible(list(manifest = d, feasibility = feas))
}

tryCatch(run(), error = function(e) {
  log_it("FAILED | Step131 v01 | ", conditionMessage(e))
  stop(e)
})
