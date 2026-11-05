options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
ROOT <- "D:/zz"
STEM <- "GSE281205"
N_EXPECTED <- 42361L
CLUSTERS <- as.character(0:23)
GSMS <- paste0("GSM", 8614662:8614667)

D123 <- file.path(ROOT, "results", paste0("123_", STEM, "_major_lineage_clustering"))
D124 <- file.path(ROOT, "results", paste0("124_", STEM, "_lineage_review_myeloid_extraction"))
D125 <- file.path(ROOT, "results", paste0("125_", STEM, "_major_lineage_annotation"))
D126 <- file.path(ROOT, "results", paste0("126_", STEM, "_targeted_lineage_review"))
D127 <- file.path(ROOT, "results", paste0("127_", STEM, "_final_lineage_adjudication_audit"))
OUT <- file.path(ROOT, "results", paste0("128_", STEM, "_reviewed_lineage_finalization"))
LOG_DIR <- file.path(OUT, "logs")
RDS_DIR <- file.path(OUT, "RDS")
for (d in c(OUT, LOG_DIR, RDS_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

IN <- c(cell_manifest = file.path(D123, "03_123_per_cell_cluster_membership_v03.csv"), old_myeloid_ids = file.path(D124,
  "06_124_myeloid_extended_cell_ids_v01.csv"), old_myeloid_rds = file.path(D124, "RDS", "124_GSE281205_extended_myeloid_DecontX_counts_v01.rds"),
  step125_evidence = file.path(D125, "02_125_24cluster_evidence_summary_v01.csv"), step126_dossier = file.path(D126,
    "07_126_cluster_review_dossier_v01.csv"), step127_recommendations = file.path(D127, "06_127_24cluster_annotation_RECOMMENDATIONS_NOT_CONFIRMED_v01.csv"))
LOG <- file.path(LOG_DIR, "00_128_run_log_v02.log")

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = LOG, append = TRUE, sep = "")
  invisible(msg)
}
assert <- function(ok, text) {
  if (!isTRUE(ok))
    stop(text, call. = FALSE)
}
read_tab <- function(path) {
  assert(file.exists(path), paste0("Missing required input: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}
write_tab <- function(x, name) {
  f <- file.path(OUT, name)
  utils::write.csv(x, f, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(f)
}
req <- function(d, cols, label) {
  missing <- setdiff(cols, names(d))
  assert(length(missing) == 0L, paste0(label, " missing columns: ", paste(missing, collapse = ", ")))
}
ord24 <- function(d, label) {
  req(d, "cluster", label)
  k <- as.character(d$cluster)
  assert(nrow(d) == 24L && !anyDuplicated(k) && setequal(k, CLUSTERS), paste0(label, " must contain clusters 0:23 once each"))
  d[match(CLUSTERS, k), , drop = FALSE]
}

expected_lineage <- c("Microglia", "Microglia", "Microglia", "Microglia", "Ependymal", "Endothelial",
  "Monocyte_macrophage", "T_cell", "Unresolved", "OPC", "Pericyte_vSMC", "Neuron", "NK_cytotoxic",
  "Vascular_ambiguous", "Fibroblast_meningeal", "Neuron", "Ependymal", "Astrocyte", "APC_dendritic",
  "Unresolved", "T_cell", "Neuron", "Microglia", "Oligodendrocyte")
expected_gate <- ifelse(CLUSTERS %in% c("0", "1", "2", "3", "6", "18", "22"), "INCLUDE", ifelse(CLUSTERS %in%
  c("8", "13", "19"), "HOLD", "EXCLUDE"))

run <- function() {

  exists <- file.exists(IN)
  write_tab(data.frame(role = names(IN), path = unname(IN), exists = exists, size_bytes = as.numeric(file.info(IN)$size)),
    "01_128_input_provenance_v02.csv")
  assert(all(exists), paste("Missing Step123/124/125/126/127 file(s):", paste(names(IN)[!exists], collapse = ", ")))

  cells <- read_tab(IN[["cell_manifest"]])
  old_ids <- read_tab(IN[["old_myeloid_ids"]])
  e125 <- ord24(read_tab(IN[["step125_evidence"]]), "Step125 evidence")
  e126 <- ord24(read_tab(IN[["step126_dossier"]]), "Step126 dossier")
  rec <- ord24(read_tab(IN[["step127_recommendations"]]), "Step127 proposals")
  req(cells, c("cell_id", "gsm", "dpi", "replicate", "cluster_harmony_r04", "cluster_harmony_r06",
    "cluster_harmony_r08", "cluster_pca_r04", "contamination_percent"), "Step123 cells")
  req(old_ids, "cell_id", "Step124 expanded myeloid ID list")
  req(e125, c("n_cells", "top12_enriched_genes", "top12_detected_canonical", "review_flags", "median_nCount_RNA",
    "median_nFeature_RNA", "pct_contamination_ge50"), "Step125 evidence")
  req(e126, c("review_priority", "focus_question"), "Step126 dossier")
  req(rec, c("n_cells", "proposed_lineage_REQUIRES_HUMAN_REVIEW", "proposed_myeloid_gate_REQUIRES_HUMAN_REVIEW",
    "proposal_caution"), "Step127 proposal")

  assert(nrow(cells) == N_EXPECTED && !anyNA(cells$cell_id) && !anyDuplicated(cells$cell_id), "Step123 manifest must have 42361 unique cell IDs")
  assert(!anyNA(cells$gsm) && setequal(as.character(cells$gsm), GSMS), "GSM list differs from the six expected libraries")
  assert(!anyNA(cells$dpi) && all(as.integer(cells$dpi) %in% c(14L, 28L)), "Unexpected stage values")
  for (g in GSMS) {
    v <- unique(as.integer(cells$dpi[as.character(cells$gsm) == g]))
    want <- if (g %in% GSMS[1:3])
      14L else 28L
    assert(length(v) == 1L && !is.na(v) && v == want, paste0("GSM / dpi mismatch: ", g))
  }
  cl <- as.character(cells$cluster_harmony_r04)
  counts <- as.integer(table(factor(cl, levels = CLUSTERS)))
  assert(all(counts == as.integer(e125$n_cells)) && all(counts == as.integer(rec$n_cells)), "Step123/125/127 cluster counts inconsistent; STOP")

  in_lineage <- trimws(as.character(rec$proposed_lineage_REQUIRES_HUMAN_REVIEW))
  in_gate <- toupper(trimws(as.character(rec$proposed_myeloid_gate_REQUIRES_HUMAN_REVIEW)))
  assert(identical(in_lineage, expected_lineage) && identical(in_gate, expected_gate), "Step127 table differs from endorsed 24-cluster map; STOP to prevent silent changes")

  curated <- data.frame(cluster = CLUSTERS, n_cells = counts, curated_major_lineage = expected_lineage,
    myeloid_gate = expected_gate, evidence_decision_status = ifelse(expected_gate == "HOLD", "UNRESOLVED_REVIEWED_AND_HELD",
      "EVIDENCE_ADOPTED_WORKING_LABEL"), annotation_basis = "Step125 marker evidence + Step126 targeted review + Step127 raw-UMI audit",
    important_caution = as.character(rec$proposal_caution), key_enriched_genes = as.character(e125$top12_enriched_genes),
    key_detected_markers = as.character(e125$top12_detected_canonical), step125_review_flags = as.character(e125$review_flags),
    step126_priority = as.character(e126$review_priority), step126_focus_question = as.character(e126$focus_question),
    median_nCount_RNA = as.numeric(e125$median_nCount_RNA), median_nFeature_RNA = as.numeric(e125$median_nFeature_RNA),
    pct_contamination_ge50 = as.numeric(e125$pct_contamination_ge50), stringsAsFactors = FALSE)
  assert(sum(curated$n_cells) == N_EXPECTED, "Cluster totals do not equal all cells")
  write_tab(curated, "02_128_curated_24cluster_annotation_v02.csv")

  idx <- match(cl, curated$cluster)
  assert(!anyNA(idx), "Some cell clusters have no curated lineage")
  annotated <- cells
  annotated$curated_major_lineage_128 <- curated$curated_major_lineage[idx]
  annotated$myeloid_gate_128 <- curated$myeloid_gate[idx]
  annotated$evidence_decision_status_128 <- curated$evidence_decision_status[idx]
  annotated$curation_version_128 <- "Step128_v02_evidence_adopted"
  assert(nrow(annotated) == N_EXPECTED && !anyNA(annotated$curated_major_lineage_128), "Curated annotation mapping failed")
  write_tab(annotated, "03_128_per_cell_curated_lineage_v02.csv")

  included <- annotated[annotated$myeloid_gate_128 == "INCLUDE", , drop = FALSE]
  held <- annotated[annotated$myeloid_gate_128 == "HOLD", , drop = FALSE]
  assert(nrow(included) == 28716L && nrow(held) == 2017L && sum(annotated$myeloid_gate_128 == "EXCLUDE") ==
    11628L, "Unexpected 28,716/2,017/11,628 gate totals")
  assert(nrow(old_ids) == nrow(included) && !anyNA(old_ids$cell_id) && !anyDuplicated(old_ids$cell_id) &&
    setequal(as.character(old_ids$cell_id), as.character(included$cell_id)), "Step124 candidate IDs do not match final Step128 include IDs; STOP")
  write_tab(included[, c("cell_id", "gsm", "dpi", "replicate", "cluster_harmony_r04", "curated_major_lineage_128",
    "contamination_percent"), drop = FALSE], "06_128_myeloid_include_cell_ids_v02.csv")
  write_tab(held[, c("cell_id", "gsm", "dpi", "replicate", "cluster_harmony_r04", "curated_major_lineage_128",
    "contamination_percent"), drop = FALSE], "07_128_unresolved_HOLD_cell_ids_v02.csv")

  gsm_lineage <- as.data.frame(table(gsm = factor(as.character(annotated$gsm), levels = GSMS), lineage = factor(annotated$curated_major_lineage_128,
    levels = unique(curated$curated_major_lineage)), gate = factor(annotated$myeloid_gate_128, levels = c("INCLUDE",
    "EXCLUDE", "HOLD"))), stringsAsFactors = FALSE)
  names(gsm_lineage)[4] <- "n_cells"
  gsm_lineage <- gsm_lineage[gsm_lineage$n_cells > 0L, , drop = FALSE]
  gsm_lineage$dpi <- ifelse(as.character(gsm_lineage$gsm) %in% GSMS[1:3], 14L, 28L)
  write_tab(gsm_lineage[, c("gsm", "dpi", "lineage", "gate", "n_cells")], "04_128_GSM_lineage_counts_v02.csv")
  stage_lineage <- aggregate(n_cells ~ dpi + lineage + gate, data = gsm_lineage, FUN = sum)
  write_tab(stage_lineage[, c("dpi", "lineage", "gate", "n_cells")], "05_128_stage_lineage_counts_v02.csv")
  per_gsm <- do.call(rbind, lapply(GSMS, function(g) {
    d <- annotated[as.character(annotated$gsm) == g, , drop = FALSE]
    data.frame(gsm = g, dpi = unique(as.integer(d$dpi)), n_all_singlets = nrow(d), n_included_myeloid = sum(d$myeloid_gate_128 ==
      "INCLUDE"), n_working_microglia = sum(d$curated_major_lineage_128 == "Microglia"), n_working_monocyte_macrophage = sum(d$curated_major_lineage_128 ==
      "Monocyte_macrophage"), n_working_APC_dendritic = sum(d$curated_major_lineage_128 == "APC_dendritic"),
      n_excluded = sum(d$myeloid_gate_128 == "EXCLUDE"), n_held_unresolved = sum(d$myeloid_gate_128 ==
        "HOLD"), note = if (g %in% GSMS[4:6])
        "Young-28dpi by title/BioSample; GEO treatment field conflict retained" else "14dpi metadata consistent", stringsAsFactors = FALSE)
  }))
  assert(all(per_gsm$n_included_myeloid > 0L) && sum(per_gsm$n_included_myeloid) == nrow(included),
    "Per-GSM include counts do not reconcile")
  write_tab(per_gsm, "08_128_myeloid_per_GSM_feasibility_v02.csv")

  saveRDS(annotated, file = file.path(RDS_DIR, "128_GSE281205_curated_cell_metadata_v02.rds"), compress = "gzip")
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, "09_128_status_and_method_notes_v02.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(OUT, "10_128_sessionInfo_v02.txt"), useBytes = TRUE)
  idx <- list.files(OUT, recursive = TRUE, full.names = TRUE)
  write_tab(data.frame(path = normalizePath(idx, winslash = "/", mustWork = FALSE), size_bytes = as.numeric(file.info(idx)$size)),
    "11_128_OUTPUT_INDEX_v02.csv")
  invisible(TRUE)
}

tryCatch(run(), error = function(e) {
  log_msg("ERROR | Step128 v02 | ", conditionMessage(e))
  stop(e)
})
