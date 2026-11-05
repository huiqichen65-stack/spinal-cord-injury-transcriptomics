source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
ROOT <- PROJECT_ROOT
VERSION <- "v01"
SAMPLE_IDS <- paste0("GSM", 8614662:8614667)
EXPECTED_CELLS <- 42361L
EXPECTED_CLUSTERS <- as.character(0:23)
FOCUS_CLUSTERS <- c("5", "10", "13")
CLUSTER_COLUMN <- "cluster_harmony_r04"

DIR119 <- file.path(ROOT, "results", "119_GSE281205_scDblFinder")
DIR123 <- file.path(ROOT, "results", "123_GSE281205_major_lineage_clustering")
DIR125 <- file.path(ROOT, "results", "125_GSE281205_major_lineage_annotation")
DIR126 <- file.path(ROOT, "results", "126_GSE281205_targeted_lineage_review")
OUT <- file.path(ROOT, "results", "127_GSE281205_final_lineage_adjudication_audit")
FIG <- file.path(OUT, "figures")
LOGS <- file.path(OUT, "logs")
for (d in c(OUT, FIG, LOGS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

P119 <- file.path(DIR119, "02_119_per_cell_doublet_calls_v02.csv")
P123 <- file.path(DIR123, "03_123_per_cell_cluster_membership_v03.csv")
PMAP <- file.path(DIR123, "07_123_marker_gene_mapping_v03.csv")
P125 <- file.path(DIR125, "04_125_manual_decisions_v01.csv")
PEV <- file.path(DIR125, "02_125_24cluster_evidence_summary_v01.csv")
P126 <- file.path(DIR126, "07_126_cluster_review_dossier_v01.csv")
P126_GSM <- file.path(DIR126, "04_126_pairwise_codetection_by_cluster_GSM_v01.csv")
LOG_FILE <- file.path(LOGS, "00_127_run_log_v01.log")

stopif <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}
needed <- c(P119, P123, PMAP, P125, PEV, P126, P126_GSM)
missing <- needed[!file.exists(needed)]
stopif(length(missing) == 0L, paste("Missing required direct input file(s):", paste(missing, collapse = " ; ")))
stopif(requireNamespace("Matrix", quietly = TRUE), "Package Matrix is required: install.packages('Matrix')")
stopif(requireNamespace("ggplot2", quietly = TRUE), "Package ggplot2 is required: install.packages('ggplot2')")

cat("", file = LOG_FILE)
msg <- function(...) {
  line <- paste0(..., collapse = "")
  cat(line, "\n", sep = "")
  cat(line, "\n", file = LOG_FILE, sep = "", append = TRUE)
  invisible(line)
}
read_csv <- function(f) {
    utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}
write_csv <- function(x, name) {
    stopif(is.data.frame(x), paste("Expected data.frame:", name))
    utils::write.csv(x, file.path(OUT, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
stat_median <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}
stat_quantile <- function(x, p) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}
percent <- function(x) 100 * mean(x, na.rm = TRUE)

SYMS <- c("Pecam1", "Cldn5", "Kdr", "Vwf", "Rgs5", "Pdgfrb",
          "Acta2", "Tagln", "P2ry12", "Fcrls", "Clec7a", "Ccr2")
MUST_HAVE <- c("Pecam1", "Cldn5", "Rgs5", "Pdgfrb")

PROPOSED_LINEAGES <- c(
  "Microglia",
  "Microglia",
  "Microglia",
  "Microglia",
  "Ependymal",
  "Endothelial",
  "Monocyte_macrophage",
  "T_cell",
  "Unresolved",
  "OPC",
  "Pericyte_vSMC",
  "Neuron",
  "NK_cytotoxic",
  "Vascular_ambiguous",
  "Fibroblast_meningeal",
  "Neuron",
  "Ependymal",
  "Astrocyte",
  "APC_dendritic",
  "Unresolved",
  "T_cell",
  "Neuron",
  "Microglia",
  "Oligodendrocyte"
)
PROPOSED_GATES <- ifelse(PROPOSED_LINEAGES %in% c("Microglia", "Monocyte_macrophage", "APC_dendritic"), 
    "INCLUDE", "EXCLUDE")
PROPOSED_GATES[as.integer(c("8", "13", "19")) + 1L] <- "HOLD"
CAUTIONS <- rep("Review canonical markers, enriched genes and cross-GSM consistency.", 24L)
CAUTIONS[1] <- "Activated Clec7a-associated microglial candidate; transcriptomic identity does NOT prove ontogeny or target 1320-gene state."
CAUTIONS[4] <- "Activated Clec7a-associated microglial candidate; check macrophage markers and independence from validation gene set."
CAUTIONS[9] <- "Very high estimated DecontX contamination / low UMI; HOLD; do not assign a biological lineage without additional evidence."
CAUTIONS[11] <- "Mural-cell marker evidence with a high-contamination subset; inspect QC and within-GSM reliability."
CAUTIONS[14] <- "Endothelial and mural marker co-detection; possible mixed phenotype or undetected heterotypic doublets. HOLD pending raw-UMI audit."
CAUTIONS[19] <- "APC/dendritic evidence (Flt3/Cd209a/Cd74/H2-Aa); do not confuse with neutrophil by auto score."
CAUTIONS[20] <- "Unresolved lineage / low depth and elevated contamination; HOLD."
CAUTIONS[23] <- "Proliferating microglial candidate; cell-cycle markers alone do not establish lineage."

run <- function() {

  mem <- read_csv(P123)
  db <- read_csv(P119)
  mapping <- read_csv(PMAP)
  decisions <- read_csv(P125)
  ev <- read_csv(PEV)
  dossier <- read_csv(P126)
  gsm_evidence <- read_csv(P126_GSM)
  stopif(all(c("cell_id", "gsm", "dpi", CLUSTER_COLUMN,
               "contamination_percent") %in% names(mem)),
         "Step123 per-cell columns do not match expected schema")
  stopif(all(c("cell_id", "gsm", "predicted_singlet", "scDblFinder_score",
               "nCount_RNA", "nFeature_RNA") %in% names(db)),
         "Step119 per-cell doublet calls schema changed")
  stopif(all(c("gene_symbol", "gene_id", "mapping_status") %in% names(mapping)),
         "Step123 marker mapping schema changed")
  stopif(all(c("cluster", "manual_lineage", "review_status", "myeloid_gate")
             %in% names(decisions)),
         "Step125 manual decisions schema changed")
  stopif(all(c("cluster", "n_cells", "draft_lineage") %in% names(ev)),
         "Step125 evidence schema changed")
  stopif(all(c("cluster", "n_cells") %in% names(dossier)),
         "Step126 dossier schema changed")
  stopif(nrow(mem) == EXPECTED_CELLS && !anyDuplicated(mem$cell_id),
         "Step123 manifest must contain 42,361 distinct cells")
  stopif(setequal(as.character(mem$gsm), SAMPLE_IDS), "Unexpected GSM IDs")
  stopif(setequal(as.character(mem[[CLUSTER_COLUMN]]), EXPECTED_CLUSTERS),
         "Expected clusters 0 through 23")
  stopif(nrow(decisions) == 24L && !anyDuplicated(decisions$cluster) &&
         setequal(as.character(decisions$cluster), EXPECTED_CLUSTERS),
         "Step125 decision table is not exactly 24 unique clusters")
  stopif(nrow(ev) == 24L && setequal(as.character(ev$cluster), EXPECTED_CLUSTERS),
         "Step125 evidence table clusters disagree")
  stopif(nrow(dossier) == 24L && setequal(as.character(dossier$cluster), EXPECTED_CLUSTERS),
         "Step126 dossier clusters disagree")
  stopif(!anyDuplicated(db$cell_id), "Duplicate Step119 per-cell IDs")
  stopif(all(as.character(mem$cell_id) %in% as.character(db$cell_id)),
         "Some Step123 singlet IDs are missing in Step119 doublet calls")
  stopif(all(as.character(mem$gsm) ==
        as.character(db$gsm[match(mem$cell_id, db$cell_id)])),
        "Step119 / Step123 GSM misalignment")
  stopif(all(as.logical(db$predicted_singlet[match(mem$cell_id, db$cell_id)])),
         "Step123 retained a barcode classified as doublet in Step119")
  for (g in SAMPLE_IDS) {
    expected_dpi <- if (g %in% SAMPLE_IDS[1:3]) 14L else 28L
    stopif(all(as.integer(mem$dpi[mem$gsm == g]) == expected_dpi),
           paste("Stage/GSM mismatch:", g))
  }
  stopif(all(as.integer(table(factor(as.character(mem[[CLUSTER_COLUMN]]),
           levels = EXPECTED_CLUSTERS))) ==
           as.integer(ev$n_cells[match(EXPECTED_CLUSTERS,
                                            as.character(ev$cluster))])),
         "Cell counts disagree with Step125 table")
  md5_before <- unname(tools::md5sum(P125))

  map_out <- do.call(rbind, lapply(SYMS, function(sym) {
    idx <- which(as.character(mapping$gene_symbol) == sym &
                 as.character(mapping$mapping_status) == "unique")
    ids <- unique(as.character(mapping$gene_id[idx]))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    data.frame(gene_symbol = sym,
               gene_id = if (length(ids) == 1L) ids else NA_character_,
               status = if (length(ids) == 1L) "unique" else
                 if (!length(ids)) "not_unique_or_missing" else "ambiguous",
               n_unique_mapping_ids = length(ids), stringsAsFactors = FALSE)
  }))
  write_csv(map_out, "01_127_raw_marker_mapping_audit_v01.csv")
  stopif(all(map_out$status[match(MUST_HAVE, map_out$gene_symbol)] == "unique"),
         "Cannot audit cluster13 without unique Pecam1/Cldn5/Rgs5/Pdgfrb mapping")
  good <- map_out[map_out$status == "unique", , drop = FALSE]

  study <- mem[as.character(mem[[CLUSTER_COLUMN]]) %in% FOCUS_CLUSTERS,
               c("cell_id", "gsm", "dpi", CLUSTER_COLUMN,
                 "contamination_percent"), drop = FALSE]
  names(study)[names(study) == CLUSTER_COLUMN] <- "cluster"
  study$cluster <- as.character(study$cluster)
  score_idx <- match(study$cell_id, db$cell_id)
  stopif(!anyNA(score_idx), "Focused cells missing Step119 scores")
  study$scDblFinder_score <- as.numeric(db$scDblFinder_score[score_idx])
  study$nCount_RNA <- as.numeric(db$nCount_RNA[score_idx])
  study$nFeature_RNA <- as.numeric(db$nFeature_RNA[score_idx])
  stopif(all(is.finite(study$scDblFinder_score)),
         "Non-finite doublet score in candidate singlets")
  stopif(all(is.finite(study$contamination_percent)),
         "Non-finite DecontX contamination percentage")
  stopif(sum(study$cluster == "13") == 549L,
         "Expected 549 cluster13 cells; check changed Step123 output")

  raw_cols <- paste0(good$gene_symbol, "_UMI")

  per_gsm <- vector("list", length(SAMPLE_IDS))
  names(per_gsm) <- SAMPLE_IDS
  file_audit <- vector("list", length(SAMPLE_IDS))
  for (g in SAMPLE_IDS) {
    f <- file.path(DIR119, "RDS",
                   paste0("119_", g, "_singlet_raw_counts_v02.rds"))
    stopif(file.exists(f), paste("Missing direct Step119 raw RDS:", f))
    item <- readRDS(f)
    stopif(is.list(item) && all(c("counts", "cell_metadata",
          "feature_metadata") %in% names(item)),
          paste("Unexpected Step119 RDS structure:", g))
    xx <- item$counts
    stopif(inherits(xx, "sparseMatrix"), paste("Not sparse counts:", g))
    stopif(nrow(xx) == nrow(item$feature_metadata) &&
             ncol(xx) == nrow(item$cell_metadata),
           paste("RDS matrix/metadata dimension mismatch:", g))
    stopif(identical(colnames(xx), as.character(item$cell_metadata$cell_id)) &&
             identical(rownames(xx), as.character(item$feature_metadata$gene_id)),
           paste("RDS gene/cell IDs not aligned with metadata:", g))
    stopif(all(is.finite(xx@x)) && all(xx@x >= 0) &&
             all(abs(xx@x - round(xx@x)) < 1e-8),
           paste("Input must be raw nonnegative integer UMI:", g))
    keep <- which(study$gsm == g)
    stopif(length(keep) > 0L, paste("No focus cells for", g))
    ids <- study$cell_id[keep]
    stopif(all(ids %in% colnames(xx)), paste("Some focus cells absent:", g))
    mapped <- good$gene_id %in% rownames(xx)
    stopif(all(mapped[good$gene_symbol %in% MUST_HAVE]),
           paste("Required marker absent from raw matrix:", g))
    z <- matrix(NA_real_, nrow = length(keep), ncol = length(SYMS),
                dimnames = list(ids, paste0(SYMS, "_UMI")))
    if (any(mapped)) {
      mx <- as.matrix(xx[good$gene_id[mapped], ids, drop = FALSE])

      z[, paste0(good$gene_symbol[mapped], "_UMI")] <- t(mx)
    }
    part <- study[keep, , drop = FALSE]
    stopif(identical(rownames(z), as.character(part$cell_id)),
           paste("Focus count cell order mismatch:", g))
    part <- cbind(part, as.data.frame(z, check.names = FALSE))
    per_gsm[[g]] <- part
    file_audit[[match(g, SAMPLE_IDS)]] <- data.frame(
      gsm = g, raw_RDS = f, all_singlets_in_RDS = ncol(xx),
      focus_cells_analyzed = nrow(part), features = nrow(xx),
      required_markers_found = all(good$gene_id[good$gene_symbol %in%
              MUST_HAVE] %in% rownames(xx)), stringsAsFactors = FALSE)
    rm(item, xx, z, part)
    invisible(gc(verbose = FALSE))
  }
  focused <- do.call(rbind, per_gsm)
  rownames(focused) <- NULL
  stopif(nrow(focused) == sum(as.character(mem[[CLUSTER_COLUMN]]) %in%
                                 FOCUS_CLUSTERS),
         "Some focused cells were lost during six-GSM raw extraction")
  stopif(!anyDuplicated(focused$cell_id), "Repeated focus cell IDs")
  write_csv(do.call(rbind, file_audit), "02_127_six_GSM_raw_RDS_audit_v01.csv")
  write_csv(focused, "03_127_focus_clusters_rawUMI_doublet_per_cell_v01.csv")

  pairs <- list(
    c("Pecam1", "Cldn5"), c("Rgs5", "Pdgfrb"),
    c("Pecam1", "Rgs5"), c("Cldn5", "Pdgfrb")
  )
  groups <- rbind(
    data.frame(cluster = FOCUS_CLUSTERS, gsm = "ALL", stringsAsFactors = FALSE),
    expand.grid(cluster = FOCUS_CLUSTERS, gsm = SAMPLE_IDS,
                stringsAsFactors = FALSE)
  )
  co <- list()
  zz <- 0L
  for (i in seq_len(nrow(groups))) {
    k <- groups$cluster[i]; g <- groups$gsm[i]
    a <- focused[focused$cluster == k &
                 (g == "ALL" | focused$gsm == g), , drop = FALSE]
    for (pair in pairs) {
      for (cutoff in c(1L, 2L, 5L)) {
        x <- a[[paste0(pair[1], "_UMI")]]
        y <- a[[paste0(pair[2], "_UMI")]]
        n_eval <- sum(!is.na(x) & !is.na(y))
        xx <- !is.na(x) & (x >= cutoff)
        yy <- !is.na(y) & (y >= cutoff)
        zz <- zz + 1L
        co[[zz]] <- data.frame(
          cluster = k, gsm = g, dpi = if (g == "ALL") "ALL" else
                  if (g %in% SAMPLE_IDS[1:3]) "14" else "28",
          gene1 = pair[1], gene2 = pair[2], raw_UMI_threshold = cutoff,
          n_cells = nrow(a), n_evaluable = n_eval,
          n_co_detected = sum(xx & yy),
          pct_co_detected = if (n_eval) 100 * sum(xx & yy) / n_eval else NA_real_,
          pct_gene1_detected = if (n_eval) 100 * sum(xx) / n_eval else NA_real_,
          pct_gene2_detected = if (n_eval) 100 * sum(yy) / n_eval else NA_real_,
          stringsAsFactors = FALSE)
      }
    }
  }
  co <- do.call(rbind, co)
  write_csv(co, "04_127_rawUMI_pairwise_codetection_by_GSM_v01.csv")

  quality_groups <- list()
  ii <- 0L
  for (k in FOCUS_CLUSTERS) {
    for (g in c("ALL", SAMPLE_IDS)) {
      a <- focused[focused$cluster == k &
                   (g == "ALL" | focused$gsm == g), , drop = FALSE]
      for (cont in c("all", "DecontX_contamination_below20")) {
        aa <- if (cont == "all") a else
                  a[a$contamination_percent < 20, , drop = FALSE]
        ii <- ii + 1L
        quality_groups[[ii]] <- data.frame(
          cluster = k, gsm = g, stratum = cont, n_cells = nrow(aa),
          median_scDblFinder_score = stat_median(aa$scDblFinder_score),
          q90_scDblFinder_score = stat_quantile(aa$scDblFinder_score, 0.90),
          median_raw_UMI = stat_median(aa$nCount_RNA),
          median_raw_detected_genes = stat_median(aa$nFeature_RNA),
          median_DecontX_contamination_percent =
             stat_median(aa$contamination_percent),
          pct_contamination_ge50 = if (nrow(aa))
            percent(aa$contamination_percent >= 50) else NA_real_,
          stringsAsFactors = FALSE)
      }
    }
  }
  qs <- do.call(rbind, quality_groups)
  write_csv(qs, "05_127_cluster13_vs_references_doublet_QC_v01.csv")

  idx <- match(EXPECTED_CLUSTERS, as.character(ev$cluster))
  de <- ev[idx, , drop = FALSE]
  di <- dossier[match(EXPECTED_CLUSTERS, as.character(dossier$cluster)),
                , drop = FALSE]
  old <- decisions[match(EXPECTED_CLUSTERS,
                         as.character(decisions$cluster)), , drop = FALSE]
  proposal <- data.frame(
    cluster = EXPECTED_CLUSTERS,
    n_cells = as.integer(de$n_cells),
    step125_draft_lineage = as.character(de$draft_lineage),
    proposed_lineage_REQUIRES_HUMAN_REVIEW = PROPOSED_LINEAGES,
    proposed_myeloid_gate_REQUIRES_HUMAN_REVIEW = PROPOSED_GATES,
    proposal_caution = CAUTIONS,
    step125_existing_manual_lineage = as.character(old$manual_lineage),
    step125_existing_review_status = as.character(old$review_status),
    step125_existing_myeloid_gate = as.character(old$myeloid_gate),
    step126_review_priority = as.character(di$review_priority),
    evidence_reference_125 = basename(PEV),
    evidence_reference_126 = basename(P126),
    requires_human_adjudication = TRUE,
    stringsAsFactors = FALSE)
  stopif(all(proposal$requires_human_adjudication),
         "Never silently confirm cluster annotations")
  write_csv(proposal, "06_127_24cluster_annotation_RECOMMENDATIONS_NOT_CONFIRMED_v01.csv")

  c13 <- co[co$cluster == "13" & co$gsm == "ALL" &
            co$raw_UMI_threshold == 1L, , drop = FALSE]
  summary13 <- merge(c13, qs[qs$cluster == "13" & qs$gsm == "ALL" &
                                 qs$stratum == "all",
                             c("cluster", "n_cells", "median_scDblFinder_score",
                               "q90_scDblFinder_score", "median_raw_UMI",
                               "median_raw_detected_genes"), drop = FALSE],
                     by = c("cluster", "n_cells"), all.x = TRUE)
  write_csv(summary13, "07_127_cluster13_raw_evidence_summary_v01.csv")

  plot_data <- co[co$gsm != "ALL" & co$raw_UMI_threshold == 1L &
                   co$gene1 %in% c("Pecam1", "Rgs5") &
                   co$gene2 %in% c("Cldn5", "Pdgfrb"), , drop = FALSE]
  plot_data$pair <- paste(plot_data$gene1, "+", plot_data$gene2)
  plot_data$gsm <- factor(plot_data$gsm, levels = SAMPLE_IDS)
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = gsm, y = pct_co_detected,
                                                 color = cluster, group = cluster)) +
    ggplot2::geom_point(size = 2.0) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(~ pair, ncol = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::labs(title = "Raw UMI marker co-detection by original GSM",
       subtitle = "Clusters 5 (endothelial), 10 (mural) and 13 (ambiguous); >=1 raw UMI",
       x = "GSM biological library", y = "% cells co-detecting marker pair",
       color = "Cluster") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ggplot2::ggsave(filename = file.path(FIG,
    "Fig_127_rawUMI_marker_codetection_by_GSM_v01.png"), plot = p1,
    width = 12, height = 8, dpi = 160)

  p2 <- ggplot2::ggplot(focused, ggplot2::aes(x = factor(cluster,
                            levels = FOCUS_CLUSTERS), y = scDblFinder_score)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey85") +
    ggplot2::geom_jitter(width = 0.18, alpha = 0.12, size = 0.4) +
    ggplot2::facet_wrap(~ gsm, scales = "free_y") +
    ggplot2::labs(title = "scDblFinder score among retained singlets",
       subtitle = "A low score does not exclude undetected heterotypic doublets",
       x = "Whole-atlas cluster", y = "scDblFinder score") +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(filename = file.path(FIG,
    "Fig_127_cluster13_vs_vascular_references_doublet_scores_v01.png"),
    plot = p2, width = 12, height = 7, dpi = 160)

  p3 <- ggplot2::ggplot(focused, ggplot2::aes(x = factor(cluster,
                            levels = FOCUS_CLUSTERS), y = nCount_RNA)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey85") +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~ gsm) +
    ggplot2::labs(title = "Original UMI depth of vascular reference/ambiguous clusters",
       x = "Whole-atlas cluster", y = "Original UMI count (log10 scale)") +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(filename = file.path(FIG,
    "Fig_127_original_UMI_depth_by_cluster_GSM_v01.png"),
    plot = p3, width = 12, height = 7, dpi = 160)

writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(OUT, "08_127_REPORT_v01.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()),
             file.path(OUT, "09_127_sessionInfo_v01.txt"), useBytes = TRUE)
  src_manifest <- data.frame(
    source = c(needed, file.path(DIR119, "RDS",
                 paste0("119_", SAMPLE_IDS, "_singlet_raw_counts_v02.rds"))),
    purpose = c("scDblFinder per-cell scores", "Step123 cluster and cell IDs",
       "canonical marker Ensembl mapping", "Step125 human decisions (READ ONLY)",
       "Step125 atlas annotation evidence", "Step126 cluster review dossier",
       "Step126 per-GSM marker co-detection", rep("original integer UMI singlets", 6)),
    stringsAsFactors = FALSE)
  stopif(nrow(src_manifest) == 13L, "Unexpected number of source files")
  src_manifest$exists <- file.exists(src_manifest$source)
  write_csv(src_manifest, "10_127_input_provenance_v01.csv")
  md5_after <- unname(tools::md5sum(P125))
  stopif(identical(md5_before, md5_after),
         "Step125 manual decisions unexpectedly changed during Step127")
  invisible(TRUE)
}

tryCatch(run(), error = function(e) {
  msg("ERROR | Step127 v01 | ", conditionMessage(e))
  stop(e)
})
