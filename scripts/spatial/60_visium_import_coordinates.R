source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = 1200)

set.seed(1234)

ROOT <- "D:/zz"

DATA_DIR <- file.path(ROOT, "data", "GSE195783")

STEP59_DIR <- file.path(ROOT, "results", "59_GSE195783_spatial_input_audit")

OUT_DIR <- file.path(ROOT, "results", "60_GSE195783_Visium_import_coordinate_reconstruction")

EXTRACT_DIR <- file.path(OUT_DIR, "00_extracted_processed_data_v05")

OBJECT_DIR <- file.path(OUT_DIR, "objects_v05")

STPIPE_TMP <- file.path(OUT_DIR, "_stPipe_V1_reconstruction_v05")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)

RAW_TAR <- file.path(DATA_DIR, "GSE195783_RAW.tar")

MANIFEST_FILE <- file.path(STEP59_DIR, "05_59_GSE195783_sample_manifest_v02.csv")

VERSION <- "v05"

COORD_VERSION <- "V1"

REQUIRED_MATCH_RATE <- 1

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "yaml", "stPipe")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step 60 v05.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(yaml)
  library(stPipe)
})

if (!exists("Run_Loc_Match", where = asNamespace("stPipe"), inherits = FALSE)) {
  stop("stPipe is installed, but stPipe::Run_Loc_Match is unavailable.")
}

write_csv_utf8 <- write_csv_utf8

strip_10x_suffix <- function(x) {
    sub("-[0-9]+$", "", x)
}

safe_median <- finite_median

read_10x_gene_expression <- function(matrix_dir) {
    x <- Seurat::Read10X(data.dir = matrix_dir, gene.column = 2, unique.features = TRUE, strip.suffix = FALSE)
    if (is.list(x)) {
        if ("Gene Expression" %in% names(x)) {
            x <- x[["Gene Expression"]]
        }
        else if (length(x) == 1L) {
            x <- x[[1]]
        }
        else {
            stop("Read10X returned multiple feature types, but ", "'Gene Expression' was not found.\nAvailable: ", 
                paste(names(x), collapse = ", "))
        }
    }
    if (!inherits(x, "dgCMatrix")) {
        x <- as(x, "dgCMatrix")
    }
    x
}

run_v1_coordinate_match <- function(barcodes_original, umi_counts, scratch_dir) {

  if (dir.exists(scratch_dir)) {
    unlink(scratch_dir, recursive = TRUE, force = TRUE)
  }

  dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)

  barcodes_stripped <- strip_10x_suffix(barcodes_original)

  if (anyDuplicated(barcodes_stripped) > 0L) {
    stop("Duplicated barcodes appeared after removing 10x suffix.")
  }

  synthetic_names <- sprintf("SPOT_%05d", seq_along(barcodes_original))

  sample_index <- data.frame(barcode_sequence = barcodes_stripped, cell_name = synthetic_names, stringsAsFactors = FALSE)

  write.csv(sample_index, file = file.path(scratch_dir, "sample_index.csv"), row.names = FALSE, quote = FALSE)

  mini_gene_count <- data.frame(gene_id = "__TOTAL_UMI__", stringsAsFactors = FALSE, check.names = FALSE)

  for (i in seq_along(synthetic_names)) {

    mini_gene_count[[synthetic_names[i]]] <- as.numeric(umi_counts[i])
  }

  write.csv(mini_gene_count, file = file.path(scratch_dir, "gene_count.csv"), row.names = FALSE, quote = FALSE)

  config_list <- list(output_directory = scratch_dir, data_directory = scratch_dir, technology_version = "Visium_polyA",
    visium_coordination = COORD_VERSION)

  config_file <- file.path(scratch_dir, "config_V1.yml")

  yaml::write_yaml(config_list, config_file)

  ans <- stPipe::Run_Loc_Match(config = config_file, pixel = FALSE, show.config = FALSE)

  matched <- ans$matched_data

  if (is.null(matched)) {
    stop("stPipe::Run_Loc_Match returned no matched_data.")
  }

  required_cols <- c("barcode_sequence")

  missing_cols <- setdiff(required_cols, colnames(matched))

  if (length(missing_cols) > 0L) {
    stop("stPipe matched_data is missing: ", paste(missing_cols, collapse = ", "))
  }

  x_col <- intersect(c("X_coordinate", "x_coordinate", "X", "x"), colnames(matched))

  y_col <- intersect(c("Y_coordinate", "y_coordinate", "Y", "y"), colnames(matched))

  if (length(x_col) == 0L || length(y_col) == 0L) {

    stop("Could not identify coordinate columns in stPipe matched_data.\n", "Columns: ", paste(colnames(matched),
      collapse = ", "))
  }

  x_col <- x_col[1]
  y_col <- y_col[1]

  matched$barcode_sequence <- as.character(matched$barcode_sequence)

  matched$barcode_original <- barcodes_original[match(matched$barcode_sequence, barcodes_stripped)]

  matched <- matched[!is.na(matched$barcode_original), , drop = FALSE]

  matched <- matched[!duplicated(matched$barcode_original), , drop = FALSE]

  coord <- data.frame(barcode = matched$barcode_original, array_x = as.numeric(matched[[x_col]]), array_y = as.numeric(matched[[y_col]]),
    stringsAsFactors = FALSE)

  coord
}

if (!file.exists(RAW_TAR)) {
  stop("Cannot find RAW TAR:\n", RAW_TAR)
}

if (!file.exists(MANIFEST_FILE)) {
  stop("Cannot find Step 59 v02 manifest:\n", MANIFEST_FILE)
}

manifest <- read.csv(MANIFEST_FILE, stringsAsFactors = FALSE, check.names = FALSE)

required_manifest_cols <- c("GSM", "inferred_timepoint", "inferred_replicate", "zip_file")

missing_manifest_cols <- setdiff(required_manifest_cols, colnames(manifest))

if (length(missing_manifest_cols) > 0L) {
  stop("Step 59 manifest is missing column(s): ", paste(missing_manifest_cols, collapse = ", "))
}

if (nrow(manifest) != 8L) {
  stop("Expected 8 GSE195783 samples, found ", nrow(manifest), ".")
}

if (any(is.na(manifest$inferred_timepoint))) {
  stop("Step 59 manifest still contains unresolved timepoints.")
}

observed_tp <- sort(unique(manifest$inferred_timepoint))

if (!setequal(observed_tp, TIME_LEVELS)) {
  stop("Unexpected timepoint structure.\nObserved: ", paste(observed_tp, collapse = ", "))
}

tp_counts <- table(manifest$inferred_timepoint)

if (!all(tp_counts[TIME_LEVELS] == 2L)) {
  stop("Expected exactly 2 samples per timepoint.")
}

manifest$inferred_timepoint <- factor(manifest$inferred_timepoint, levels = TIME_LEVELS)

manifest <- manifest %>% arrange(inferred_timepoint, inferred_replicate, GSM) %>% group_by(inferred_timepoint) %>% 
    mutate(timepoint_replicate = row_number()) %>% ungroup()

manifest$sample_label <- paste0(manifest$inferred_timepoint, "_rep", manifest$timepoint_replicate)

manifest$coordinate_version <- COORD_VERSION

print(manifest[, c("GSM", "inferred_timepoint", "timepoint_replicate", "sample_label", "zip_file")])

if (dir.exists(EXTRACT_DIR)) {
  unlink(EXTRACT_DIR, recursive = TRUE, force = TRUE)
}

dir.create(EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)

TAR_DIR <- file.path(EXTRACT_DIR, "_tar")

dir.create(TAR_DIR, recursive = TRUE, showWarnings = FALSE)

utils::untar(RAW_TAR, exdir = TAR_DIR)

zip_paths <- list.files(TAR_DIR, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

if (length(zip_paths) != 8L) {
  stop("Expected 8 processed ZIP files in RAW.tar, found ", length(zip_paths), ".")
}

zip_lookup <- setNames(zip_paths, basename(zip_paths))

missing_zip <- setdiff(manifest$zip_file, basename(zip_paths))

if (length(missing_zip) > 0L) {
  stop("ZIP file(s) listed in Step 59 manifest were not found:\n", paste(missing_zip, collapse = "\n"))
}

sample_dirs <- character(nrow(manifest))

matrix_dirs <- character(nrow(manifest))

for (i in seq_len(nrow(manifest))) {
    gsm <- manifest$GSM[i]
    sample_dir <- file.path(EXTRACT_DIR, gsm)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(zipfile = zip_lookup[manifest$zip_file[i]], exdir = sample_dir)
    candidates <- list.dirs(sample_dir, recursive = TRUE, full.names = TRUE)
    candidates <- candidates[basename(candidates) == "filtered_feature_bc_matrix"]
    if (length(candidates) != 1L) {
        stop(gsm, ": expected exactly one filtered_feature_bc_matrix directory; found ", length(candidates), 
            ".")
    }
    sample_dirs[i] <- sample_dir
    matrix_dirs[i] <- candidates[1]
}

manifest$extracted_sample_dir <- sample_dirs
manifest$matrix_dir <- matrix_dirs

write_csv_utf8(manifest, file.path(OUT_DIR, paste0("00_60_GSE195783_final_sample_manifest_", VERSION,
  ".csv")))

if (dir.exists(STPIPE_TMP)) {
  unlink(STPIPE_TMP, recursive = TRUE, force = TRUE)
}

dir.create(STPIPE_TMP, recursive = TRUE, showWarnings = FALSE)

counts_list <- vector("list", nrow(manifest))

names(counts_list) <- manifest$GSM

coordinate_list <- vector("list", nrow(manifest))

names(coordinate_list) <- manifest$GSM

match_audit_list <- vector("list", nrow(manifest))

preqc_list <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
    gsm <- manifest$GSM[i]
    sample_label <- manifest$sample_label[i]
    counts <- read_10x_gene_expression(manifest$matrix_dir[i])
    if (nrow(counts) == 0L || ncol(counts) == 0L) {
        stop(gsm, ": empty count matrix.")
    }
    if (is.null(colnames(counts))) {
        stop(gsm, ": missing spot barcodes.")
    }
    if (anyDuplicated(colnames(counts)) > 0L) {
        stop(gsm, ": duplicated count-matrix barcodes.")
    }
    counts_list[[gsm]] <- counts
    ncount <- Matrix::colSums(counts)
    nfeature <- Matrix::colSums(counts > 0)
    mt_idx <- grepl("^mt-", rownames(counts), ignore.case = TRUE)
    if (any(mt_idx)) {
        percent_mt <- 100 * Matrix::colSums(counts[mt_idx, , drop = FALSE])/pmax(ncount, 1)
    }
    else {
        percent_mt <- rep(NA_real_, ncol(counts))
    }
    coords <- run_v1_coordinate_match(barcodes_original = colnames(counts), umi_counts = ncount, scratch_dir = file.path(STPIPE_TMP, 
        gsm))
    n_total <- ncol(counts)
    n_matched <- sum(colnames(counts) %in% coords$barcode)
    match_rate <- n_matched/n_total
    n_extra_coord <- sum(!coords$barcode %in% colnames(counts))
    match_audit_list[[i]] <- data.frame(GSM = gsm, timepoint = as.character(manifest$inferred_timepoint[i]), 
        replicate = manifest$timepoint_replicate[i], sample_label = sample_label, coordinate_version = COORD_VERSION, 
        n_spots_matrix = n_total, n_spots_matched = n_matched, match_rate = match_rate, n_extra_coordinate_barcodes = n_extra_coord, 
        pass_100pct_match = identical(n_matched, n_total) && n_extra_coord == 0L, stringsAsFactors = FALSE)
    coordinate_list[[gsm]] <- coords
    preqc_list[[i]] <- data.frame(GSM = gsm, timepoint = as.character(manifest$inferred_timepoint[i]), 
        replicate = manifest$timepoint_replicate[i], sample_label = sample_label, coordinate_version = COORD_VERSION, 
        n_genes_matrix = nrow(counts), n_spots_filtered = ncol(counts), median_UMI = safe_median(as.numeric(ncount)), 
        median_detected_genes = safe_median(as.numeric(nfeature)), median_percent_mt = safe_median(as.numeric(percent_mt)), 
        min_UMI = min(ncount), max_UMI = max(ncount), stringsAsFactors = FALSE)
}

match_audit <- bind_rows(match_audit_list)

preqc_summary <- bind_rows(preqc_list)

write_csv_utf8(match_audit, file.path(OUT_DIR, paste0("01_60_GSE195783_V1_coordinate_match_audit_", VERSION,
  ".csv")))

write_csv_utf8(preqc_summary, file.path(OUT_DIR, paste0("02_60_GSE195783_preQC_sample_summary_", VERSION,
  ".csv")))

if (!all(match_audit$pass_100pct_match)) {
    failed <- match_audit[!match_audit$pass_100pct_match, , drop = FALSE]
    stop("V1 coordinate reconstruction did NOT achieve exact barcode coverage ", "for all samples.\nFailed sample(s):\n", 
        paste(paste0(failed$GSM, " [", failed$n_spots_matched, "/", failed$n_spots_matrix, "]"), collapse = "\n"), 
        "\nNo final Seurat objects were created.")
}

if (!all(abs(match_audit$match_rate - REQUIRED_MATCH_RATE) < .Machine$double.eps^0.5)) {
  stop("At least one sample has match_rate != 1.0.")
}

if (dir.exists(OBJECT_DIR)) {
  unlink(OBJECT_DIR, recursive = TRUE, force = TRUE)
}

dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)

seurat_list <- vector("list", nrow(manifest))

names(seurat_list) <- manifest$GSM

all_spot_meta_list <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {

  gsm <- manifest$GSM[i]
  counts <- counts_list[[gsm]]
  coords <- coordinate_list[[gsm]]

  idx <- match(colnames(counts), coords$barcode)

  if (any(is.na(idx))) {
    stop(gsm, ": barcode alignment failed after V1 validation.")
  }

  coords <- coords[idx, , drop = FALSE]

  if (!identical(coords$barcode, colnames(counts))) {
    stop(gsm, ": coordinate order does not exactly match matrix columns.")
  }

  ncount <- Matrix::colSums(counts)

  nfeature <- Matrix::colSums(counts > 0)

  mt_idx <- grepl("^mt-", rownames(counts), ignore.case = TRUE)

  if (any(mt_idx)) {

    percent_mt <- 100 * Matrix::colSums(counts[mt_idx, , drop = FALSE])/pmax(ncount, 1)

  } else {

    percent_mt <- rep(NA_real_, ncol(counts))
  }

  spot_meta <- data.frame(row.names = colnames(counts), barcode = colnames(counts), GSM = gsm, timepoint = as.character(manifest$inferred_timepoint[i]),
    replicate = manifest$timepoint_replicate[i], GEO_replicate = manifest$inferred_replicate[i],
    sample_label = manifest$sample_label[i], coordinate_system = "Visium_V1_array", histology_registered = FALSE,
    array_x = coords$array_x, array_y = coords$array_y, nCount_Spatial = as.numeric(ncount), nFeature_Spatial = as.numeric(nfeature),
    percent_mt = as.numeric(percent_mt), stringsAsFactors = FALSE)

  seu <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial", project = paste0("GSE195783_",
    gsm), meta.data = spot_meta, min.cells = 0, min.features = 0)

  spatial_embedding <- as.matrix(spot_meta[, c("array_x", "array_y")])

  colnames(spatial_embedding) <- c("ARRAY_1", "ARRAY_2")

  rownames(spatial_embedding) <- rownames(spot_meta)

  seu[["spatial_array"]] <- SeuratObject::CreateDimReducObject(embeddings = spatial_embedding, key = "ARRAY_",
    assay = "Spatial")

  seu@misc$GSE195783 <- list(dataset = "GSE195783", coordinate_system = "Visium V1 array coordinates",
    histology_registered = FALSE, normalization_performed = FALSE, filtering_performed_in_step60 = FALSE,
    source_matrix = "GEO filtered_feature_bc_matrix", step = 60L, version = VERSION)

  object_file <- file.path(OBJECT_DIR, paste0(gsm, "_", manifest$sample_label[i], "_VisiumV1_array_",
    VERSION, ".rds"))

  saveRDS(seu, file = object_file, compress = TRUE)

  seurat_list[[gsm]] <- seu

  all_spot_meta_list[[i]] <- data.frame(barcode = spot_meta$barcode, GSM = spot_meta$GSM, timepoint = spot_meta$timepoint,
    replicate = spot_meta$replicate, GEO_replicate = spot_meta$GEO_replicate, sample_label = spot_meta$sample_label,
    coordinate_system = spot_meta$coordinate_system, histology_registered = spot_meta$histology_registered,
    array_x = spot_meta$array_x, array_y = spot_meta$array_y, nCount_Spatial = spot_meta$nCount_Spatial,
    nFeature_Spatial = spot_meta$nFeature_Spatial, percent_mt = spot_meta$percent_mt, stringsAsFactors = FALSE)
}

all_spot_meta <- bind_rows(all_spot_meta_list)

all_spot_meta$timepoint <- factor(all_spot_meta$timepoint, levels = TIME_LEVELS)

all_spot_meta$sample_label <- factor(all_spot_meta$sample_label, levels = manifest$sample_label)

spot_meta_file <- file.path(OUT_DIR, paste0("03_60_GSE195783_final_spot_metadata_", VERSION, ".csv.gz"))

con_spot <- gzfile(spot_meta_file, open = "wt")

write.csv(all_spot_meta, file = con_spot, row.names = FALSE, na = "")

close(con_spot)

seurat_list_file <- file.path(OUT_DIR, paste0("05_60_GSE195783_spatial_Seurat_list_", VERSION, ".rds"))

saveRDS(seurat_list, file = seurat_list_file, compress = TRUE)

p_coord <- ggplot(all_spot_meta, aes(x = array_x, y = array_y)) + geom_point(size = 0.48, alpha = 0.92) +
  coord_equal() + scale_y_reverse() + facet_wrap(~sample_label, ncol = 4) + labs(title = "GSE195783 reconstructed Visium V1 array coordinates",
  subtitle = paste0("Filtered GEO spots; 100% barcode matching; ", "array coordinates only, without H&E/pixel registration"),
  x = "Visium array X", y = "Visium array Y") + theme_classic(base_family = "Arial", base_size = 10) +
  theme(plot.title = element_text(family = "Arial", size = 13), plot.subtitle = element_text(family = "Arial",
    size = 9), strip.text = element_text(family = "Arial", size = 8.5), axis.title = element_text(family = "Arial",
    size = 9), axis.text = element_text(family = "Arial", size = 7))

coord_pdf <- file.path(OUT_DIR, paste0("04_60_GSE195783_final_V1_array_coordinate_overview_", VERSION,
  ".pdf"))

coord_png <- file.path(OUT_DIR, paste0("04_60_GSE195783_final_V1_array_coordinate_overview_", VERSION,
  ".png"))

ggsave(filename = coord_pdf, plot = p_coord, width = 11, height = 7.2, units = "in", device = grDevices::cairo_pdf,
  bg = "white")

ggsave(filename = coord_png, plot = p_coord, width = 11, height = 7.2, units = "in", dpi = 300, bg = "white")

final_checks <- data.frame(check = c("Expected 8 samples", "Exactly 2 samples at each 3/7/14/28 dpi",
  "All matrices non-empty", "All samples use fixed Visium V1 coordinates", "All samples have 100% barcode-coordinate matching",
  "No extra coordinate barcodes", "All final spot metadata rows have coordinates", "No H&E registration claimed",
  "No normalization performed in Step 60", "Final Seurat list contains 8 objects"), pass = c(nrow(manifest) ==
  8L, all(table(manifest$inferred_timepoint)[TIME_LEVELS] == 2L), all(preqc_summary$n_genes_matrix >
  0 & preqc_summary$n_spots_filtered > 0), all(match_audit$coordinate_version == COORD_VERSION), all(match_audit$pass_100pct_match),
  all(match_audit$n_extra_coordinate_barcodes == 0L), all(is.finite(all_spot_meta$array_x) & is.finite(all_spot_meta$array_y)),
  all(all_spot_meta$histology_registered == FALSE), TRUE, length(seurat_list) == 8L), stringsAsFactors = FALSE)

write_csv_utf8(final_checks, file.path(OUT_DIR, paste0("06_60_GSE195783_FINAL_QC_", VERSION, ".csv")))

if (!all(final_checks$pass)) {
  warning("One or more final Step 60 checks did not pass. ", "Inspect 06_60_GSE195783_FINAL_QC_v05.csv before proceeding.")
}

report_file <- file.path(OUT_DIR, paste0("07_60_GSE195783_FINAL_REPORT_", VERSION, ".txt"))

con <- file(report_file, open = "wt", encoding = "UTF-8")

writeLines(c("============================================================", "GSE195783 Step 60 FINAL",
  "Visium matrix import + V1 array-coordinate reconstruction", paste0("Version: ", VERSION), "============================================================",
  "", paste0("RAW TAR: ", RAW_TAR), paste0("Step59 manifest: ", MANIFEST_FILE), paste0("Output directory: ",
    OUT_DIR), "", "SAMPLE STRUCTURE", paste0("  Total samples: ", nrow(manifest)), paste0("  3 dpi: ",
    sum(as.character(manifest$inferred_timepoint) == "3dpi")), paste0("  7 dpi: ", sum(as.character(manifest$inferred_timepoint) ==
    "7dpi")), paste0("  14 dpi: ", sum(as.character(manifest$inferred_timepoint) == "14dpi")), paste0("  28 dpi: ",
    sum(as.character(manifest$inferred_timepoint) == "28dpi")), "", "COORDINATE RECONSTRUCTION",
  paste0("  Fixed coordinate system: ", COORD_VERSION), paste0("  Samples with 100% barcode matching: ",
    sum(match_audit$pass_100pct_match), "/", nrow(match_audit)), paste0("  Minimum match rate: ",
    sprintf("%.6f", min(match_audit$match_rate))), paste0("  Maximum extra coordinate barcodes: ",
    max(match_audit$n_extra_coordinate_barcodes)), "", "DATA STATUS", "  Input expression: GEO filtered_feature_bc_matrix",
  "  Normalization performed: NO", "  Spot filtering performed in Step 60: NO", "  Clustering performed: NO",
  "  RCTD/cell2location performed: NO", "  Clec7a/lipid scoring performed: NO", "", "SPATIAL INTERPRETATION LIMIT",
  paste0("  Coordinates are reconstructed Visium V1 ARRAY coordinates."), paste0("  They are NOT H&E-registered pixel coordinates."),
  paste0("  No scalefactor or histology image was fabricated."), "", "FINAL OBJECT", paste0("  ", seurat_list_file),
  "", "NEXT STEP", paste0("  Step 61 can perform formal spatial QC/filtering and normalization ", "using these Step 60 objects."),
  "", "FINAL QC", paste(paste0(ifelse(final_checks$pass, "[PASS] ", "[CHECK] "), final_checks$check),
    collapse = "\n"), "", "============================================================"), con = con)

close(con)

print(match_audit)
print(preqc_summary)
print(final_checks)
