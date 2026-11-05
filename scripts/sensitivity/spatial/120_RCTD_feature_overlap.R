source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1)
ROOT <- PROJECT_ROOT
VERSION <- "v02"
RUN_ID <- "120_RCTD_lipid_feature_overlap"
RECONSTRUCT_RCTD_PREPROCESSING <- TRUE

SAMPLES <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1",
  "28dpi_rep2")

D65 <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")
D73 <- file.path(ROOT, "results", "73_GSE195783_Clec7a_lipid_spatial_coupling")
D64 <- file.path(ROOT, "results", "64_GSE195783_spatial_cluster_marker_region_annotation")
OUT <- file.path(ROOT, "results", RUN_ID)
TABLES <- file.path(OUT, "tables")
dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)

LIPID_FILE <- file.path(D73, "01_73_lipid_program_gene_matching_v03.csv")
REF_FILE <- file.path(D65, "00_65_hierarchical_scRNA_reference_v02.rds")
SPATIAL_FILE <- file.path(D64, "12_64_GSE195783_spatial_region_annotated_v03.rds")
PREVIOUS_SAMPLE_QC <- file.path(D65, "02_65_RCTD_sample_QC_v03.csv")
LOG_FILE <- file.path(OUT, paste0(RUN_ID, "_", VERSION, ".log"))
REPORT_FILE <- file.path(OUT, paste0("08_120_REPORT_", VERSION, ".txt"))

write_csv <- function(x, stem) {
  path <- file.path(TABLES, paste0(stem, "_", VERSION, ".csv"))
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  path
}

log_message <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = LOG_FILE, sep = "", append = TRUE)
  invisible(msg)
}

path_check <- function(path, purpose) {
  data.frame(purpose = purpose, absolute_path = normalizePath(path, winslash = "/", mustWork = FALSE),
    exists = file.exists(path), size_bytes = if (file.exists(path))
      unname(file.info(path)$size) else NA_real_, stringsAsFactors = FALSE)
}

clean_genes <- function(x) {
  x <- unique(trimws(as.character(x)))
  x[!is.na(x) & nzchar(x) & x != "NA"]
}

extract_feature_lists <- function(preprocessed) {
  iv <- NULL
  if (is.list(preprocessed) && !is.null(preprocessed$internal_vars)) {
    iv <- preprocessed$internal_vars
  }
  if (is.null(iv) && methods::isS4(preprocessed) && "internal_vars" %in% methods::slotNames(preprocessed)) {
    iv <- methods::slot(preprocessed, "internal_vars")
  }
  if (is.null(iv)) {
    iv <- tryCatch(spacexr::internal_vars(preprocessed), error = function(e) NULL)
  }
  candidates <- list(iv, preprocessed)
  pull <- function(key) {
    for (obj in candidates) {
      if (is.null(obj))
        next
      z <- tryCatch(obj[[key]], error = function(e) NULL)
      if (is.character(z) || is.factor(z)) {
        z <- clean_genes(z)
        if (length(z) > 0L)
          return(z)
      }
    }
    character()
  }
  list(regression = pull("gene_list_reg"), platform = pull("gene_list_bulk"))
}

get_spatial_counts <- function(obj) {
  assay_name <- "Spatial"
  if (!assay_name %in% names(obj@assays)) {
    stop("Missing Spatial assay")
  }
  assay_obj <- obj[[assay_name]]
  m <- NULL
  if (inherits(assay_obj, "Assay5")) {

    layer_names <- SeuratObject::Layers(assay_obj, search = NA)
    count_layers <- layer_names[grepl("^counts($|\\.)", layer_names)]
    if (length(count_layers) == 0L) {
      stop("Assay5 Spatial assay has no counts layer. Available: ", paste(layer_names, collapse = ", "))
    }
    if (length(count_layers) == 1L) {

      m <- SeuratObject::LayerData(assay_obj, layer = count_layers[[1L]])
    } else {
      joined_assay <- tryCatch(SeuratObject::JoinLayers(object = assay_obj, layers = "counts",
        new = "counts"), error = function(e) {
        stop("Could not join multiple Spatial counts layers using ", "SeuratObject::JoinLayers: ",
          conditionMessage(e))
      })
      m <- SeuratObject::LayerData(joined_assay, layer = "counts")
      rm(joined_assay)
    }
  } else {

    m <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay_name, slot = "counts"), error = function(e) {
      stop("Could not extract Spatial v3 counts: ", conditionMessage(e))
    })
  }
  if (is.null(m) || length(dim(m)) != 2L) {
    stop("Spatial raw counts were not retrieved as a 2D matrix")
  }
  if (!inherits(m, "dgCMatrix"))
    m <- methods::as(m, "dgCMatrix")
  if (is.null(rownames(m)) || is.null(colnames(m)) || anyDuplicated(rownames(m)) || anyDuplicated(colnames(m))) {
    stop("Spatial counts have missing or duplicated gene/spot identifiers")
  }
  m
}

cat("Step120 startup: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = LOG_FILE, sep = "")

files <- do.call(rbind, list(path_check(LIPID_FILE, "Step73 gene-matching manifest"), path_check(REF_FILE,
  "Step65 reference RDS"), path_check(SPATIAL_FILE, "Step64 spatial RDS"), path_check(PREVIOUS_SAMPLE_QC,
  "Step65 original sample QC (optional)")))
write_csv(files, "00_120_input_inventory")
if (!file.exists(LIPID_FILE)) stop("Required file missing: ", LIPID_FILE)

lipid <- utils::read.csv(LIPID_FILE, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
needed_cols <- c("module", "requested_gene", "matched_gene", "detected")
if (!all(needed_cols %in% colnames(lipid))) {
  stop("Step73 gene manifest is missing required columns: ", paste(setdiff(needed_cols, names(lipid)),
    collapse = ", "))
}
core <- lipid[as.character(lipid$module) == "Core_lipid_handling_union", , drop = FALSE]
if (nrow(core) == 0L) stop("Missing Core_lipid_handling_union in Step73")
core$detected <- tolower(trimws(as.character(core$detected))) %in% c("true", "t", "1")
core$requested_gene <- trimws(as.character(core$requested_gene))
core$matched_gene <- trimws(as.character(core$matched_gene))
core$matched_gene[is.na(core$matched_gene)] <- ""
if (anyDuplicated(core$requested_gene) || anyDuplicated(core$matched_gene[core$detected])) {
  stop("Unexpected duplicate genes in Step73 core program")
}
scoring_genes <- clean_genes(core$matched_gene[core$detected])
requested_genes <- clean_genes(core$requested_gene)
if ("Clec7a" %in% scoring_genes) {
  stop("Clec7a unexpectedly present in Step73 core score")
}
core$included_in_original_spatial_score <- core$detected & nzchar(core$matched_gene)
write_csv(core, "01_120_frozen_Step73_lipid_union_genes")
if (length(scoring_genes) != 37L) {
  log_message("WARNING: scored genes != 37; check Step73 source versions")
}

source_status <- data.frame(sample_label = SAMPLES, feature_source = "not_available", sample_status = "PENDING",
  n_shared_input_genes = NA_integer_, n_RCTD_regression_features = NA_integer_, n_RCTD_platform_features = NA_integer_,
  diagnostic = "", stringsAsFactors = FALSE)
reg_by_sample <- setNames(vector("list", length(SAMPLES)), SAMPLES)
bulk_by_sample <- setNames(vector("list", length(SAMPLES)), SAMPLES)

if (RECONSTRUCT_RCTD_PREPROCESSING) {
  required <- c("Seurat", "SeuratObject", "Matrix", "SummarizedExperiment", "SpatialExperiment", "S4Vectors",
    "spacexr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("Missing R packages: ", paste(missing, collapse = ", "))
  if (!all(file.exists(c(REF_FILE, SPATIAL_FILE)))) {
    stop("Step65 reference / Step64 spatial input RDS missing. ", "No substitute inputs are allowed.")
  }
  log_message("Package versions: Seurat=", as.character(utils::packageVersion("Seurat")), "; SeuratObject=",
    as.character(utils::packageVersion("SeuratObject")), "; spacexr=", as.character(utils::packageVersion("spacexr")))
  reference <- readRDS(REF_FILE)
  spatial <- readRDS(SPATIAL_FILE)
  if (!inherits(reference, "SummarizedExperiment"))
    stop("Unexpected reference class")
  if (!inherits(spatial, "Seurat"))
    stop("Unexpected spatial object class")
  md <- spatial[[]]
  if (!all(c("sample_label", "GSM", "timepoint", "array_x", "array_y", "spatial_region_v03") %in% colnames(md))) {
    stop("Step64 spatial metadata missing a field required by Step65 v02")
  }
  counts <- get_spatial_counts(spatial)
  genes_common <- intersect(rownames(reference), rownames(counts))
  if (length(genes_common) < 5000L)
    stop("<5000 shared input genes; unexpected source")
  reference <- reference[genes_common, , drop = FALSE]
  counts <- counts[genes_common, , drop = FALSE]
  if (!identical(colnames(counts), rownames(md))) {
    if (!setequal(colnames(counts), rownames(md))) {
      stop("Spatial counts/metadata spot identifiers do not match")
    }
    md <- md[colnames(counts), , drop = FALSE]
  }

  if (file.exists(PREVIOUS_SAMPLE_QC)) {
    previous_qc <- utils::read.csv(PREVIOUS_SAMPLE_QC, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    previous_qc <- NULL
  }

  for (i in seq_along(SAMPLES)) {
    sid <- SAMPLES[[i]]
    one <- tryCatch({
      cell_names <- rownames(md)[as.character(md$sample_label) == sid]
      if (length(cell_names) < 100L)
        stop("Fewer than 100 spots")
      cts <- counts[, cell_names, drop = FALSE]
      mds <- md[cell_names, , drop = FALSE]
      coords <- as.matrix(mds[, c("array_x", "array_y"), drop = FALSE])
      storage.mode(coords) <- "numeric"
      colnames(coords) <- c("x", "y")
      rownames(coords) <- cell_names
      coldata <- S4Vectors::DataFrame(GSM = as.character(mds$GSM), sample_label = as.character(mds$sample_label),
        timepoint = as.character(mds$timepoint), spatial_region_v03 = as.character(mds$spatial_region_v03),
        nUMI = as.numeric(Matrix::colSums(cts)), row.names = cell_names)
      spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = cts), colData = coldata,
        spatialCoords = coords)

      set.seed(1234)
      prepared <- spacexr::createRctd(spatial_experiment = spe, reference_experiment = reference,
        cell_type_col = "cell_type", require_int = TRUE, ref_UMI_min = 100, ref_n_cells_min = 25,
        ref_n_cells_max = 2000)
      genes <- extract_feature_lists(prepared)
      rm(prepared, spe, cts, coldata)
      gc(verbose = FALSE)
      if (!length(genes$regression) || !length(genes$platform)) {
        stop("createRctd succeeded but neither/both feature lists were not exposed. ", "This version may not retain gene_list_reg/gene_list_bulk.")
      }
      if (!all(genes$regression %in% genes_common) || !all(genes$platform %in% genes_common)) {
        stop("Selected features absent from common input genes: inspect ID mapping")
      }
      if (!is.null(previous_qc) && all(c("sample_label", "common_genes") %in% colnames(previous_qc))) {
        old <- previous_qc[previous_qc$sample_label == sid, "common_genes"]
        if (length(old) == 1L && is.finite(as.numeric(old)) && as.integer(old) != length(genes_common)) {
          stop("common-gene count disagrees with original Step65 v03 QC")
        }
      }
      genes
    }, error = function(e) e)

    if (inherits(one, "error")) {
      source_status$sample_status[i] <- "FAILED"
      source_status$feature_source[i] <- "reconstruction_failed"
      source_status$diagnostic[i] <- conditionMessage(one)
      log_message("  FAILED: ", conditionMessage(one))
    } else {
      reg_by_sample[[sid]] <- sort(one$regression)
      bulk_by_sample[[sid]] <- sort(one$platform)
      source_status$sample_status[i] <- "COMPLETE_RECONSTRUCTED"
      source_status$feature_source[i] <- "Step65_v02_frozen_inputs_createRctd_reconstruction"
      source_status$n_shared_input_genes[i] <- length(genes_common)
      source_status$n_RCTD_regression_features[i] <- length(one$regression)
      source_status$n_RCTD_platform_features[i] <- length(one$platform)
    }
  }
  rm(reference, spatial, counts)
  gc(verbose = FALSE)
} else {
  source_status$sample_status <- "NOT_RUN_NO_RDS"
  source_status$diagnostic <- "Preprocessing reconstruction disabled by user setting"
}

write_csv(source_status, "02_120_feature_source_status")

complete_samples <- SAMPLES[source_status$sample_status == "COMPLETE_RECONSTRUCTED"]
all_selected <- list()
overlap_rows <- list()
pergene_rows <- list()

for (sid in SAMPLES) {
  reg <- reg_by_sample[[sid]]
  bulk <- bulk_by_sample[[sid]]
  ok <- sid %in% complete_samples
  all_features <- union(reg, bulk)
  if (ok) {
    all_selected[[sid]] <- data.frame(sample_label = sid, gene = sort(all_features), RCTD_regression_fit = sort(all_features) %in%
      reg, RCTD_platform_normalization = sort(all_features) %in% bulk, source_type = "RECONSTRUCTED_createRctd_not_original_persisted",
      stringsAsFactors = FALSE)
  }
  overlap_rows[[sid]] <- data.frame(sample_label = sid, feature_source_status = if (ok)
    "RECONSTRUCTED" else "UNAVAILABLE", n_original_scoring_genes = length(scoring_genes), n_RCTD_regression = if (ok)
    length(reg) else NA_integer_, n_RCTD_platform = if (ok)
    length(bulk) else NA_integer_, n_RCTD_either = if (ok)
    length(all_features) else NA_integer_, n_scoring_overlap_regression = if (ok)
    length(intersect(reg, scoring_genes)) else NA_integer_, n_scoring_overlap_platform = if (ok)
    length(intersect(bulk, scoring_genes)) else NA_integer_, n_scoring_overlap_either = if (ok)
    length(intersect(all_features, scoring_genes)) else NA_integer_, n_scoring_nonoverlap_both = if (ok)
    length(setdiff(scoring_genes, all_features)) else NA_integer_, fraction_scoring_overlap_either = if (ok)
    length(intersect(all_features, scoring_genes))/length(scoring_genes) else NA_real_, stringsAsFactors = FALSE)
  pergene_rows[[sid]] <- data.frame(sample_label = sid, gene = scoring_genes, is_scored_in_original_Step73 = TRUE,
    in_RCTD_regression_fit = if (ok)
      scoring_genes %in% reg else NA, in_RCTD_platform_normalization = if (ok)
      scoring_genes %in% bulk else NA, in_either_RCTD_selected_set = if (ok)
      scoring_genes %in% all_features else NA, eligible_for_nonoverlap_in_this_area = if (ok)
      !(scoring_genes %in% all_features) else NA, feature_evidence = if (ok)
      "RECONSTRUCTED" else "UNAVAILABLE", stringsAsFactors = FALSE)
}
selected_df <- if (length(all_selected)) do.call(rbind, all_selected) else data.frame(sample_label = character(),
  gene = character(), RCTD_regression_fit = logical(), RCTD_platform_normalization = logical(), source_type = character())
write_csv(selected_df, "03_120_reconstructed_RCTD_selected_feature_lists")
overlap_df <- do.call(rbind, overlap_rows)
write_csv(overlap_df, "04_120_overlap_by_capture_area")
pergene_df <- do.call(rbind, pergene_rows)
write_csv(pergene_df, "05_120_overlap_by_lipid_gene_and_area")

all_complete <- length(complete_samples) == length(SAMPLES)

if (all_complete) {
  global_reg <- unique(unlist(reg_by_sample, use.names = FALSE))
  global_bulk <- unique(unlist(bulk_by_sample, use.names = FALSE))
  global_either <- union(global_reg, global_bulk)
  count_reg <- vapply(scoring_genes, function(g) sum(vapply(reg_by_sample, function(v) g %in% v, logical(1))),
    integer(1))
  count_bulk <- vapply(scoring_genes, function(g) sum(vapply(bulk_by_sample, function(v) g %in% v,
    logical(1))), integer(1))
  count_either <- vapply(scoring_genes, function(g) sum(vapply(SAMPLES, function(sid) g %in% union(reg_by_sample[[sid]],
    bulk_by_sample[[sid]]), logical(1))), integer(1))
  candidate <- data.frame(gene = scoring_genes, n_areas_RCTD_regression = count_reg, n_areas_RCTD_platform = count_bulk,
    n_areas_either = count_either, any_RCTD_selected_feature_overlap = scoring_genes %in% global_either,
    nonoverlap_ALL_8_areas_and_both_sets = !(scoring_genes %in% global_either), evidence = "RECONSTRUCTED_ALL_8_AREAS",
    stringsAsFactors = FALSE)
} else {
  candidate <- data.frame(gene = scoring_genes, n_areas_RCTD_regression = NA_integer_, n_areas_RCTD_platform = NA_integer_,
    n_areas_either = NA_integer_, any_RCTD_selected_feature_overlap = NA, nonoverlap_ALL_8_areas_and_both_sets = NA,
    evidence = "INCOMPLETE_AUDIT_DO_NOT_CLAIM_NONOVERLAP", stringsAsFactors = FALSE)
}
write_csv(candidate, "06_120_global_nonoverlap_candidate_genes")

modules <- setdiff(unique(as.character(lipid$module)), "Core_lipid_handling_union")
module_rows <- list()
for (module_name in modules) {
  part <- lipid[as.character(lipid$module) == module_name & tolower(trimws(as.character(lipid$detected))) %in%
    c("true", "t", "1"), , drop = FALSE]
  genes <- clean_genes(part$matched_gene)
  for (sid in SAMPLES) {
    complete <- sid %in% complete_samples
    sel <- union(reg_by_sample[[sid]], bulk_by_sample[[sid]])
    module_rows[[length(module_rows) + 1L]] <- data.frame(module = module_name, sample_label = sid,
      n_detected_genes = length(genes), n_RCTD_selected_overlap = if (complete)
        length(intersect(genes, sel)) else NA_integer_, n_nonoverlap_genes = if (complete)
        length(setdiff(genes, sel)) else NA_integer_, evidence = if (complete)
        "RECONSTRUCTED" else "UNAVAILABLE", stringsAsFactors = FALSE)
  }
}
module_df <- if (length(module_rows)) do.call(rbind, module_rows) else data.frame(module = character(),
  sample_label = character(), n_detected_genes = integer(), n_RCTD_selected_overlap = integer(), n_nonoverlap_genes = integer(),
  evidence = character())
write_csv(module_df, "07_120_overlap_by_lipid_module")

qc <- data.frame(check = c("Step73 core program present", "Clec7a excluded from scoring", "Step73 unique detected scored genes = 37",
  "All eight sample-specific RCTD feature lists reconstructed", "Globally nonoverlapping set can be claimed from complete reconstruction",
  "runRctd NOT called", "Original Step64/65/73 files NOT overwritten"), pass = c(length(requested_genes) >
  0L, !("Clec7a" %in% scoring_genes), length(scoring_genes) == 37L, all_complete, all_complete, TRUE,
  TRUE), stringsAsFactors = FALSE)
write_csv(qc, "09_120_FINAL_QC")

lines <- c("============================================================", "Step 120 | RCTD lipid feature overlap",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "============================================================",
  paste("Reference source:", REF_FILE), paste("Spatial source:", SPATIAL_FILE), paste("Step73 gene manifest:",
    LIPID_FILE), paste("RCTD reconstruction enabled:", RECONSTRUCT_RCTD_PREPROCESSING), "", "IMPORTANT: gene_list_reg and gene_list_bulk recovered here, if present,",
  "are reconstructed by createRctd on inputs, NOT an exact export", "from the original Step65 v02 preprocessing object.",
  "No RCTD deconvolution, spot scoring or correlation was repeated.", "", paste("Step73 requested unique core lipid genes:",
    length(requested_genes)), paste("Step73 detected/scored unique genes:", length(scoring_genes)),
  paste("Clec7a excluded:", !("Clec7a" %in% scoring_genes)), paste("Sample-specific feature lists recovered:",
    length(complete_samples), "/8"), paste("Complete all-area audit:", all_complete), "")
for (sid in SAMPLES) {
  z <- overlap_df[overlap_df$sample_label == sid, , drop = FALSE]
  lines <- c(lines, if (identical(z$feature_source_status, "RECONSTRUCTED")) {
    sprintf("%s: RCTD regression=%d, platform=%d, scoring overlap either=%d/%d, nonoverlap=%d", sid,
      z$n_RCTD_regression, z$n_RCTD_platform, z$n_scoring_overlap_either, z$n_original_scoring_genes,
      z$n_scoring_nonoverlap_both)
  } else {
    paste(sid, ": feature list unavailable; ", source_status$diagnostic[source_status$sample_label ==
      sid])
  })
}
lines <- c(lines, "")
if (all_complete) {
  nonoverlap <- candidate$gene[candidate$nonoverlap_ALL_8_areas_and_both_sets]
  lines <- c(lines, paste("Candidate lipid genes not selected by either RCTD gene list in any area:",
    length(nonoverlap)), paste("Gene symbols:", if (length(nonoverlap)) paste(sort(nonoverlap), collapse = ", ") else "[NONE]"),
    "These candidates have NOT been rescored or revalidated; this audit only", "establishes feature overlap under reconstructed preprocessing.")
} else {
  lines <- c(lines, "Cannot conclude that any lipid gene is globally nonoverlapping because", "one or more sample feature lists are missing. See status and log files.")
}
lines <- c(lines, "", "INTERPRETATION LIMITATIONS", "- The Step73 gene set was not an a priori independent held-out panel.",
  "- Nonoverlap with feature selection does not guarantee independence from", "  the original single-cell reference or shared tissue-expression structure.",
  "- RCTD selected regression genes and platform-normalization genes", "  are separate feature sets; report both and their union.",
  "- Per-capture feature selection may differ; global exclusion should use", "  the UNION of both feature lists across all 8 capture areas.",
  "- Because spacexr preprocessing is reconstructed, verify installed", "  spacexr version and settings match the original analysis before",
  "  describing these as the exact historical Step65 gene selections.", "- Next step, only if the candidate set is adequate: Step65 v06",
  "  weights + Step64 log-expression to rescore and retest coupling (Step121).", "============================================================")
writeLines(lines, REPORT_FILE, useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT, paste0("10_120_sessionInfo_", VERSION, ".txt")),
  useBytes = TRUE)
if (!all_complete) {
}
