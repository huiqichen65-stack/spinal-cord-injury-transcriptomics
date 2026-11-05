source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = 1200)

set.seed(1234)

ROOT <- PROJECT_ROOT

SPATIAL_FILE <- file.path(ROOT, "results", "64_GSE195783_spatial_cluster_marker_region_annotation", "12_64_GSE195783_spatial_region_annotated_v03.rds")

GLOBAL_REF_FILE <- file.path(ROOT, "results", "09_Figure1A_major_celltype_UMAP", "09_Figure1A_major_celltype_UMAP_object_v04.rds")

MYELOID_REF_FILE <- file.path(ROOT, "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

OUT_DIR <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")

CHECKPOINT_DIR <- file.path(OUT_DIR, "checkpoints_v02")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v02"

TARGET_STATE <- "Clec7a+ lipid-handling Microglia"

MAX_REF_CELLS_PER_TYPE <- 2000L
MIN_REF_CELLS_PER_TYPE <- 25L

RCTD_MODE <- "full"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

MYELOID_MAJOR_TYPES <- c("Microglia", "Macrophage/Monocyte", "APC/DC-like myeloid")

NON_MYELOID_TYPES <- c("Neutrophil", "T/NK cell", "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Ependymal", 
    "Endothelial", "Pericyte/VSMC", "Fibroblast/Meningeal")

MYELOID_STATE_RAW <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

MYELOID_STATE_DISPLAY <- c(`Homeostatic Microglia` = "Homeostatic Microglia", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling Microglia", 
    `Activated inflammatory Microglia` = "Activated inflammatory Microglia", `Ccl12+ acute Microglia` = "Ccl12+ acute Microglia", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory Microglia", `Proliferating Microglia` = "Proliferating Microglia", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ antigen-presenting Macrophage", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb+ tissue-repair Macrophage", 
    `Acute-response Macrophage` = "Acute-response Macrophage")

MYELOID_STATE_DISPLAY_LEVELS <- unname(MYELOID_STATE_DISPLAY[MYELOID_STATE_RAW])

TARGET_SANITY_GENES <- c("Lpl", "Apoe", "Trem2", "Abca1", "Cst7", "Apoc1")

detected_cores <- parallel::detectCores(logical = FALSE)

if (is.na(detected_cores) || detected_cores < 2L) {
    detected_cores <- 2L
}

N_CORES <- max(1L, min(4L, as.integer(detected_cores - 1L)))

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "tidyr", "patchwork", "SummarizedExperiment", 
    "SpatialExperiment", "S4Vectors", "spacexr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall all required packages before running Step 65.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(S4Vectors)
  library(spacexr)
})

write_csv_utf8 <- write_csv_utf8

save_pdf_png <- function(p, stem, width, height) {
    ggsave(filename = file.path(OUT_DIR, paste0(stem, "_", VERSION, ".pdf")), plot = p, width = width, 
        height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
    ggsave(filename = file.path(OUT_DIR, paste0(stem, "_", VERSION, ".png")), plot = p, width = width, 
        height = height, units = "in", dpi = 300, bg = "white")
}

get_counts <- function(obj, assay = "RNA") {
    if (!assay %in% names(obj@assays)) {
        stop("Assay '", assay, "' was not found.")
    }
    obj2 <- obj
    if (inherits(obj2[[assay]], "Assay5")) {
        obj2 <- tryCatch({
            JoinLayers(object = obj2, assay = assay)
        }, error = function(e) {
            stop("JoinLayers failed for assay '", assay, "': ", conditionMessage(e))
        })
    }
    ans <- tryCatch({
        SeuratObject::LayerData(object = obj2, assay = assay, layer = "counts")
    }, error = function(e) {
        NULL
    })
    if (is.null(ans)) {
        ans <- Seurat::GetAssayData(object = obj2, assay = assay, slot = "counts")
    }
    if (!inherits(ans, "dgCMatrix")) {
        ans <- as(ans, "dgCMatrix")
    }
    ans
}

get_data_layer <- function(obj, assay = "Spatial") {
    obj2 <- obj
    if (inherits(obj2[[assay]], "Assay5")) {
        obj2 <- JoinLayers(object = obj2, assay = assay)
    }
    ans <- tryCatch({
        SeuratObject::LayerData(object = obj2, assay = assay, layer = "data")
    }, error = function(e) {
        NULL
    })
    if (is.null(ans)) {
        ans <- Seurat::GetAssayData(object = obj2, assay = assay, slot = "data")
    }
    ans
}

assert_integer_counts <- function(mat, label) {
    if (length(mat@x) == 0L) {
        stop(label, ": count matrix has no non-zero entries.")
    }
    delta <- abs(mat@x - round(mat@x))
    if (any(delta > 1e-08)) {
        stop(label, ": count matrix is not integer-valued. ", "RCTD requires untransformed raw counts.")
    }
    invisible(TRUE)
}

balanced_downsample <- function(metadata, label_col, strata_col, max_per_label, seed = 1234L) {
    set.seed(seed)
    labs <- as.character(metadata[[label_col]])
    strata <- as.character(metadata[[strata_col]])
    cells <- rownames(metadata)
    keep <- character()
    for (lab in sort(unique(labs[!is.na(labs)]))) {
        idx_lab <- which(labs == lab)
        if (length(idx_lab) <= max_per_label) {
            keep <- c(keep, cells[idx_lab])
            next
        }
        strata_lab <- strata[idx_lab]
        valid_strata <- sort(unique(strata_lab[!is.na(strata_lab) & nzchar(strata_lab)]))
        if (length(valid_strata) <= 1L) {
            keep <- c(keep, sample(cells[idx_lab], max_per_label, replace = FALSE))
            next
        }
        quota <- floor(max_per_label/length(valid_strata))
        selected <- character()
        for (st in valid_strata) {
            idx_st <- idx_lab[strata_lab == st]
            n_take <- min(quota, length(idx_st))
            if (n_take > 0L) {
                selected <- c(selected, sample(cells[idx_st], n_take, replace = FALSE))
            }
        }
        remaining <- setdiff(cells[idx_lab], selected)
        n_more <- max_per_label - length(selected)
        if (n_more > 0L && length(remaining) > 0L) {
            n_more <- min(n_more, length(remaining))
            selected <- c(selected, sample(remaining, n_more, replace = FALSE))
        }
        keep <- c(keep, selected)
    }
    unique(keep)
}

safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10L) {
        return(NA_real_)
    }
    suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

sanitize_rctd_label <- function(x) {
    y <- gsub("[^A-Za-z0-9_]", "_", x)
    y <- gsub("_+", "_", y)
    y <- gsub("^_+|_+$", "", y)
    starts_bad <- grepl("^[0-9]", y)
    y[starts_bad] <- paste0("X_", y[starts_bad])
    y
}

required_files <- c(SPATIAL_FILE, GLOBAL_REF_FILE, MYELOID_REF_FILE)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
    stop("Required input file(s) missing:\n", paste(missing_files, collapse = "\n"))
}

spatial_obj <- readRDS(SPATIAL_FILE)

global_ref <- readRDS(GLOBAL_REF_FILE)

myeloid_ref <- readRDS(MYELOID_REF_FILE)

if (!inherits(spatial_obj, "Seurat")) {
    stop("Spatial input is not a Seurat object.")
}

if (!inherits(global_ref, "Seurat")) {
    stop("Global reference is not a Seurat object.")
}

if (!inherits(myeloid_ref, "Seurat")) {
    stop("Myeloid reference is not a Seurat object.")
}

required_spatial_meta <- c("GSM", "sample_label", "timepoint", "array_x", "array_y", "spatial_region_v03", 
    "histology_registered")

missing_spatial_meta <- setdiff(required_spatial_meta, colnames(spatial_obj[[]]))

if (length(missing_spatial_meta) > 0L) {
    stop("Spatial object missing metadata: ", paste(missing_spatial_meta, collapse = ", "))
}

if (any(spatial_obj$histology_registered)) {
    stop("H&E registration must remain FALSE.")
}

if (!"major_celltype" %in% colnames(global_ref[[]])) {
    stop("Global scRNA object lacks 'major_celltype'.")
}

if (!"source_dataset" %in% colnames(global_ref[[]])) {
    stop("Global scRNA object lacks 'source_dataset'.")
}

if (!"myeloid_state" %in% colnames(myeloid_ref[[]])) {
    stop("Myeloid scRNA object lacks 'myeloid_state'.")
}

if (!"source_dataset" %in% colnames(myeloid_ref[[]])) {
    stop("Myeloid scRNA object lacks 'source_dataset'.")
}

global_md <- global_ref[[]]

global_md$major_celltype <- as.character(global_md$major_celltype)

nonmy_cells_all <- rownames(global_md)[global_md$major_celltype %in% NON_MYELOID_TYPES]

if (length(nonmy_cells_all) == 0L) {
    stop("No non-myeloid cells found in global reference.")
}

global_nonmy_md <- global_md[nonmy_cells_all, , drop = FALSE]

global_nonmy_md$ref_label <- global_nonmy_md$major_celltype

nonmy_cells_keep <- balanced_downsample(metadata = global_nonmy_md, label_col = "ref_label", strata_col = "source_dataset", 
    max_per_label = MAX_REF_CELLS_PER_TYPE, seed = 1234L)

myeloid_md <- myeloid_ref[[]]

myeloid_md$myeloid_state <- as.character(myeloid_md$myeloid_state)

unexpected_states <- setdiff(unique(myeloid_md$myeloid_state[!is.na(myeloid_md$myeloid_state)]), MYELOID_STATE_RAW)

if (length(unexpected_states) > 0L) {
    stop("Unexpected myeloid_state label(s): ", paste(unexpected_states, collapse = " | "))
}

myeloid_cells_all <- rownames(myeloid_md)[myeloid_md$myeloid_state %in% MYELOID_STATE_RAW]

myeloid_md2 <- myeloid_md[myeloid_cells_all, , drop = FALSE]

myeloid_md2$ref_label <- unname(MYELOID_STATE_DISPLAY[myeloid_md2$myeloid_state])

myeloid_cells_keep <- balanced_downsample(metadata = myeloid_md2, label_col = "ref_label", strata_col = "source_dataset", 
    max_per_label = MAX_REF_CELLS_PER_TYPE, seed = 5678L)

global_nonmy_sub <- subset(global_ref, cells = nonmy_cells_keep)

myeloid_sub <- subset(myeloid_ref, cells = myeloid_cells_keep)

global_counts <- get_counts(global_nonmy_sub, assay = "RNA")

myeloid_counts <- get_counts(myeloid_sub, assay = "RNA")

assert_integer_counts(global_counts, "Global non-myeloid reference")

assert_integer_counts(myeloid_counts, "Myeloid-state reference")

common_ref_genes <- intersect(rownames(global_counts), rownames(myeloid_counts))

if (length(common_ref_genes) < 5000L) {
    stop("Too few common genes between global and myeloid references: ", length(common_ref_genes))
}

global_counts <- global_counts[common_ref_genes, , drop = FALSE]

myeloid_counts <- myeloid_counts[common_ref_genes, , drop = FALSE]

colnames(global_counts) <- paste0("GLOBAL__", colnames(global_counts))

colnames(myeloid_counts) <- paste0("MYELOID__", colnames(myeloid_counts))

ref_counts <- cbind(global_counts, myeloid_counts)

global_labels <- as.character(global_nonmy_sub$major_celltype)

myeloid_raw_labels <- as.character(myeloid_sub$myeloid_state)

myeloid_labels <- unname(MYELOID_STATE_DISPLAY[myeloid_raw_labels])

ref_labels_display <- c(global_labels, myeloid_labels)

ref_labels <- sanitize_rctd_label(ref_labels_display)

if (anyDuplicated(unique(data.frame(display = ref_labels_display, internal = ref_labels, stringsAsFactors = FALSE)$internal)) > 
    0L) {
    stop("Sanitized RCTD labels are not unique. Review the label map.")
}

label_map <- unique(data.frame(display_label = ref_labels_display, RCTD_internal_label = ref_labels, 
    stringsAsFactors = FALSE))

if (anyDuplicated(label_map$RCTD_internal_label) > 0L) {
    stop("Two display labels map to the same sanitized RCTD internal label.")
}

global_sources <- as.character(global_nonmy_sub$source_dataset)

myeloid_sources <- as.character(myeloid_sub$source_dataset)

ref_sources <- c(global_sources, myeloid_sources)

names(ref_labels) <- colnames(ref_counts)

names(ref_labels_display) <- colnames(ref_counts)

names(ref_sources) <- colnames(ref_counts)

if (any(is.na(ref_labels))) {
    stop("NA labels present in hierarchical reference.")
}

ref_type_counts <- table(ref_labels)

if (any(ref_type_counts < MIN_REF_CELLS_PER_TYPE)) {
    bad <- ref_type_counts[ref_type_counts < MIN_REF_CELLS_PER_TYPE]
    stop("Reference type(s) with <", MIN_REF_CELLS_PER_TYPE, " cells:\n", paste(paste0(names(bad), " = ", 
        as.integer(bad)), collapse = "\n"))
}

if (length(ref_type_counts) != length(NON_MYELOID_TYPES) + length(MYELOID_STATE_RAW)) {
    stop("Hierarchical reference does not contain expected 19 types.")
}

internal_to_display <- setNames(label_map$display_label, label_map$RCTD_internal_label)

reference_manifest <- data.frame(RCTD_internal_type = names(ref_type_counts), reference_type = unname(internal_to_display[names(ref_type_counts)]), 
    n_cells = as.integer(ref_type_counts), stringsAsFactors = FALSE)

reference_manifest$reference_class <- ifelse(reference_manifest$reference_type %in% MYELOID_STATE_DISPLAY_LEVELS, 
    "myeloid_state", "broad_non_myeloid")

source_manifest <- as.data.frame(table(RCTD_internal_type = ref_labels, source_dataset = ref_sources), 
    stringsAsFactors = FALSE)

source_manifest <- source_manifest[source_manifest$Freq > 0, , drop = FALSE]

source_manifest$reference_type <- unname(internal_to_display[source_manifest$RCTD_internal_type])

source_manifest <- source_manifest[, c("reference_type", "RCTD_internal_type", "source_dataset", "Freq")]

reference_manifest <- reference_manifest %>% arrange(reference_class, reference_type)

write_csv_utf8(reference_manifest, file.path(OUT_DIR, paste0("01_65_reference_manifest_", VERSION, ".csv")))

write_csv_utf8(label_map[order(label_map$display_label), , drop = FALSE], file.path(OUT_DIR, paste0("01a_65_RCTD_internal_label_map_", 
    VERSION, ".csv")))

write_csv_utf8(source_manifest, file.path(OUT_DIR, paste0("01b_65_reference_source_composition_", VERSION, 
    ".csv")))

reference_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = ref_counts), colData = S4Vectors::DataFrame(cell_type = factor(ref_labels), 
    source_dataset = ref_sources, row.names = colnames(ref_counts)))

REFERENCE_RDS <- file.path(OUT_DIR, paste0("00_65_hierarchical_scRNA_reference_", VERSION, ".rds"))

saveRDS(reference_se, REFERENCE_RDS, compress = TRUE)

rm(global_ref, myeloid_ref, global_nonmy_sub, myeloid_sub, global_counts, myeloid_counts, ref_counts)

gc()

spatial_counts <- get_counts(spatial_obj, assay = "Spatial")

assert_integer_counts(spatial_counts, "Spatial data")

common_genes <- intersect(rownames(reference_se), rownames(spatial_counts))

if (length(common_genes) < 5000L) {
    stop("Too few common reference/spatial genes: ", length(common_genes))
}

reference_se <- reference_se[common_genes, , drop = FALSE]

spatial_counts <- spatial_counts[common_genes, , drop = FALSE]

spatial_md <- spatial_obj[[]]

spatial_md$timepoint <- factor(as.character(spatial_md$timepoint), levels = TIME_LEVELS)

sample_order <- unique(as.character(spatial_md$sample_label))

expected_samples <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", 
    "28dpi_rep1", "28dpi_rep2")

if (setequal(sample_order, expected_samples)) {
    sample_order <- expected_samples
}

rctd_result_list <- list()
weights_list <- list()
sample_qc_list <- list()

for (
  i in seq_along(
    sample_order
  )
) {

  samp <- sample_order[
    i
  ]

  cells <- rownames(
    spatial_md
  )[
    as.character(
      spatial_md$sample_label
    ) ==
      samp
  ]

  if (length(
    cells
  ) < 100L) {
    stop(
      samp,
      ": fewer than 100 spatial spots."
    )
  }

  md_samp <- spatial_md[
    cells,
    ,
    drop = FALSE
  ]

  counts_samp <- spatial_counts[
    ,
    cells,
    drop = FALSE
  ]

  coords <- as.matrix(
    md_samp[
      ,
      c(
        "array_x",
        "array_y"
      ),
      drop = FALSE
    ]
  )

  colnames(
    coords
  ) <- c(
    "x",
    "y"
  )

  rownames(
    coords
  ) <- cells

  nUMI <- Matrix::colSums(
    counts_samp
  )

  coldata_samp <- S4Vectors::DataFrame(
    GSM =
      as.character(
        md_samp$GSM
      ),
    sample_label =
      as.character(
        md_samp$sample_label
      ),
    timepoint =
      as.character(
        md_samp$timepoint
      ),
    spatial_region_v03 =
      as.character(
        md_samp$spatial_region_v03
      ),
    nUMI =
      as.numeric(
        nUMI
      ),
    row.names =
      cells
  )

  spatial_spe <- SpatialExperiment::SpatialExperiment(
    assays = list(
      counts =
        counts_samp
    ),
    colData =
      coldata_samp,
    spatialCoords =
      coords
  )

  checkpoint_file <- file.path(
    CHECKPOINT_DIR,
    paste0(
      samp,
      "_RCTD_full_",
      VERSION,
      ".rds"
    )
  )

  if (file.exists(
    checkpoint_file
  )) {

    result_spe <- readRDS(
      checkpoint_file
    )

  } else {

    rctd_data <- spacexr::createRctd(
      spatial_experiment =
        spatial_spe,
      reference_experiment =
        reference_se,
      cell_type_col =
        "cell_type",
      require_int =
        TRUE,
      ref_UMI_min =
        100,
      ref_n_cells_min =
        MIN_REF_CELLS_PER_TYPE,
      ref_n_cells_max =
        MAX_REF_CELLS_PER_TYPE
    )

    result_spe <- spacexr::runRctd(
      rctd_data,
      rctd_mode =
        RCTD_MODE,
      max_cores =
        N_CORES
    )

    saveRDS(
      result_spe,
      checkpoint_file,
      compress = TRUE
    )

    rm(
      rctd_data
    )
  }

  if (!"weights" %in%
      SummarizedExperiment::assayNames(
        result_spe
      )) {
    stop(
      samp,
      ": RCTD result lacks assay 'weights'."
    )
  }

  weights <- SummarizedExperiment::assay(
    result_spe,
    "weights"
  )

  unmapped_weight_types <- setdiff(
    rownames(
      weights
    ),
    names(
      internal_to_display
    )
  )

  if (length(
    unmapped_weight_types
  ) > 0L) {
    stop(
      samp,
      ": RCTD returned unmapped internal cell-type label(s): ",
      paste(
        unmapped_weight_types,
        collapse = " | "
      )
    )
  }

  rownames(
    weights
  ) <- unname(
    internal_to_display[
      rownames(
        weights
      )
    ]
  )

  if (!TARGET_STATE %in%
      rownames(
        weights
      )) {
    stop(
      samp,
      ": target state missing from RCTD weights:\n",
      TARGET_STATE,
      "\nAvailable types:\n",
      paste(
        rownames(
          weights
        ),
        collapse = " | "
      )
    )
  }

  spot_common <- intersect(
    colnames(
      weights
    ),
    cells
  )

  if (length(
    spot_common
  ) == 0L) {
    stop(
      samp,
      ": no RCTD result spot names match input spots."
    )
  }

  weights <- weights[
    ,
    spot_common,
    drop = FALSE
  ]

  weight_sums <- Matrix::colSums(
    weights
  )

  qc_row <- data.frame(
    sample_label = samp,
    timepoint =
      as.character(
        unique(
          md_samp$timepoint
        )[
          1
        ]
      ),
    n_spots_input =
      length(
        cells
      ),
    n_spots_RCTD =
      ncol(
        weights
      ),
    fraction_spots_returned =
      ncol(
        weights
      ) /
      length(
        cells
      ),
    n_reference_types =
      nrow(
        weights
      ),
    mean_weight_sum =
      mean(
        weight_sums
      ),
    median_weight_sum =
      median(
        weight_sums
      ),
    min_weight_sum =
      min(
        weight_sums
      ),
    max_weight_sum =
      max(
        weight_sums
      ),
    common_genes =
      length(
        common_genes
      ),
    checkpoint_file =
      checkpoint_file,
    stringsAsFactors = FALSE
  )

  sample_qc_list[[samp]] <- qc_row

  weights_df <- as.data.frame(
    t(
      as.matrix(
        weights
      )
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  weights_df$barcode <- rownames(
    weights_df
  )

  md_out <- spatial_md[
    weights_df$barcode,
    ,
    drop = FALSE
  ]

  weights_df$GSM <-
    as.character(
      md_out$GSM
    )

  weights_df$sample_label <-
    as.character(
      md_out$sample_label
    )

  weights_df$timepoint <-
    as.character(
      md_out$timepoint
    )

  weights_df$spatial_region_v03 <-
    as.character(
      md_out$spatial_region_v03
    )

  weights_df$array_x <-
    md_out$array_x

  weights_df$array_y <-
    md_out$array_y

  weights_list[[samp]] <- weights_df

  rctd_result_list[[samp]] <- checkpoint_file

  rm(
    result_spe,
    weights,
    weights_df,
    spatial_spe,
    counts_samp
  )

  gc()
}

sample_qc <- bind_rows(sample_qc_list)

sample_qc$timepoint <- factor(sample_qc$timepoint, levels = TIME_LEVELS)

sample_qc <- sample_qc %>% arrange(timepoint, sample_label)

write_csv_utf8(sample_qc, file.path(OUT_DIR, paste0("02_65_RCTD_sample_QC_", VERSION, ".csv")))

spot_weights <- bind_rows(weights_list)

spot_weights$timepoint <- factor(spot_weights$timepoint, levels = TIME_LEVELS)

spot_weights$sample_label <- factor(spot_weights$sample_label, levels = sample_order)

SPOT_WEIGHT_FILE <- file.path(OUT_DIR, paste0("03_65_RCTD_spot_weights_", VERSION, ".csv.gz"))

con_weights <- gzfile(SPOT_WEIGHT_FILE, open = "wt")

write.csv(spot_weights, file = con_weights, row.names = FALSE, na = "")

close(con_weights)

target_sample_summary <- spot_weights %>% group_by(timepoint, sample_label) %>% summarise(n_spots = n(), 
    mean_target_weight = mean(.data[[TARGET_STATE]], na.rm = TRUE), median_target_weight = median(.data[[TARGET_STATE]], 
        na.rm = TRUE), q25_target_weight = quantile(.data[[TARGET_STATE]], 0.25, na.rm = TRUE), q75_target_weight = quantile(.data[[TARGET_STATE]], 
        0.75, na.rm = TRUE), fraction_weight_gt_0_10 = mean(.data[[TARGET_STATE]] > 0.1, na.rm = TRUE), 
    .groups = "drop")

target_time_summary <- target_sample_summary %>% group_by(timepoint) %>% summarise(n_replicates = n(), 
    mean_of_sample_means = mean(mean_target_weight), sd_of_sample_means = stats::sd(mean_target_weight), 
    min_sample_mean = min(mean_target_weight), max_sample_mean = max(mean_target_weight), .groups = "drop")

target_time_output <- bind_rows(target_sample_summary %>% mutate(summary_level = "sample") %>% rename(mean_weight = mean_target_weight, 
    median_weight = median_target_weight) %>% mutate(n_replicates = NA_integer_, sd_of_sample_means = NA_real_, 
    min_sample_mean = NA_real_, max_sample_mean = NA_real_) %>% select(summary_level, timepoint, sample_label, 
    n_spots, n_replicates, mean_weight, median_weight, q25_target_weight, q75_target_weight, fraction_weight_gt_0_10, 
    sd_of_sample_means, min_sample_mean, max_sample_mean), target_time_summary %>% mutate(summary_level = "timepoint_across_replicates", 
    sample_label = NA_character_, n_spots = NA_integer_, mean_weight = mean_of_sample_means, median_weight = NA_real_, 
    q25_target_weight = NA_real_, q75_target_weight = NA_real_, fraction_weight_gt_0_10 = NA_real_) %>% 
    select(summary_level, timepoint, sample_label, n_spots, n_replicates, mean_weight, median_weight, 
        q25_target_weight, q75_target_weight, fraction_weight_gt_0_10, sd_of_sample_means, min_sample_mean, 
        max_sample_mean))

write_csv_utf8(target_time_output, file.path(OUT_DIR, paste0("04_65_Clec7a_timepoint_summary_", VERSION, 
    ".csv")))

target_region_sample <- spot_weights %>% group_by(timepoint, sample_label, spatial_region_v03) %>% summarise(n_spots = n(), 
    mean_target_weight = mean(.data[[TARGET_STATE]], na.rm = TRUE), median_target_weight = median(.data[[TARGET_STATE]], 
        na.rm = TRUE), fraction_weight_gt_0_10 = mean(.data[[TARGET_STATE]] > 0.1, na.rm = TRUE), .groups = "drop")

target_region_time <- target_region_sample %>% group_by(timepoint, spatial_region_v03) %>% summarise(n_replicates = n_distinct(sample_label), 
    mean_of_sample_means = mean(mean_target_weight), sd_of_sample_means = ifelse(n() > 1L, stats::sd(mean_target_weight), 
        NA_real_), min_sample_mean = min(mean_target_weight), max_sample_mean = max(mean_target_weight), 
    .groups = "drop")

target_region_output <- bind_rows(target_region_sample %>% mutate(summary_level = "sample_region", n_replicates = NA_integer_, 
    sd_of_sample_means = NA_real_, min_sample_mean = NA_real_, max_sample_mean = NA_real_) %>% select(summary_level, 
    timepoint, sample_label, spatial_region_v03, n_spots, n_replicates, mean_target_weight, median_target_weight, 
    fraction_weight_gt_0_10, sd_of_sample_means, min_sample_mean, max_sample_mean), target_region_time %>% 
    mutate(summary_level = "timepoint_region_across_replicates", sample_label = NA_character_, n_spots = NA_integer_, 
        mean_target_weight = mean_of_sample_means, median_target_weight = NA_real_, fraction_weight_gt_0_10 = NA_real_) %>% 
    select(summary_level, timepoint, sample_label, spatial_region_v03, n_spots, n_replicates, mean_target_weight, 
        median_target_weight, fraction_weight_gt_0_10, sd_of_sample_means, min_sample_mean, max_sample_mean))

write_csv_utf8(target_region_output, file.path(OUT_DIR, paste0("05_65_Clec7a_region_summary_", VERSION, 
    ".csv")))

plot_df <- spot_weights

plot_df$sample_label <- factor(plot_df$sample_label, levels = sample_order)

p_target_map <- ggplot(plot_df, aes(x = array_x, y = array_y, color = .data[[TARGET_STATE]])) + geom_point(size = 0.55, 
    alpha = 0.95) + coord_equal() + scale_y_reverse() + scale_color_viridis_c(option = "magma", limits = c(0, 
    1), oob = scales::squish) + facet_wrap(~sample_label, ncol = 4) + labs(title = "RCTD: Clec7a+ lipid-handling Microglia", 
    subtitle = "Estimated spot weight; Visium V1 array coordinates, no H&E registration", x = "Visium array X", 
    y = "Visium array Y", color = "RCTD\nweight") + theme_classic(base_family = "Arial", base_size = 9) + 
    theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8.5), strip.text = element_text(size = 8), 
        legend.position = "right")

save_pdf_png(p = p_target_map, stem = "06_65_Clec7a_RCTD_spatial_maps", width = 11.7, height = 7.5)

region_plot_df <- target_region_sample

region_plot_df$timepoint <- factor(region_plot_df$timepoint, levels = TIME_LEVELS)

p_target_region <- ggplot(region_plot_df, aes(x = timepoint, y = mean_target_weight, group = sample_label)) + 
    geom_point(size = 2.2) + geom_line(linewidth = 0.45, alpha = 0.65) + facet_wrap(~spatial_region_v03, 
    scales = "free_y", ncol = 3) + labs(title = "Clec7a+ lipid-handling Microglia RCTD weight by spatial region", 
    subtitle = "Each point is one biological spatial sample; descriptive only", x = NULL, y = "Mean RCTD weight per sample-region") + 
    theme_classic(base_family = "Arial", base_size = 8.5) + theme(axis.text.x = element_text(angle = 35, 
    hjust = 1), strip.text = element_text(size = 7.5), plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8))

save_pdf_png(p = p_target_region, stem = "07_65_Clec7a_region_time_diagnostic", width = 11, height = 9)

myeloid_cols_present <- intersect(MYELOID_STATE_DISPLAY_LEVELS, colnames(spot_weights))

if (length(myeloid_cols_present) != length(MYELOID_STATE_DISPLAY_LEVELS)) {
    stop("Not all 9 myeloid states are present in RCTD weights.")
}

myeloid_long <- spot_weights %>% select(sample_label, timepoint, all_of(myeloid_cols_present)) %>% pivot_longer(cols = all_of(myeloid_cols_present), 
    names_to = "myeloid_state", values_to = "weight")

myeloid_sample_mean <- myeloid_long %>% group_by(sample_label, timepoint, myeloid_state) %>% summarise(mean_weight = mean(weight, 
    na.rm = TRUE), .groups = "drop")

myeloid_time_mean <- myeloid_sample_mean %>% group_by(timepoint, myeloid_state) %>% summarise(mean_weight = mean(mean_weight), 
    .groups = "drop")

myeloid_time_mean$timepoint <- factor(myeloid_time_mean$timepoint, levels = TIME_LEVELS)

myeloid_time_mean$myeloid_state <- factor(myeloid_time_mean$myeloid_state, levels = rev(MYELOID_STATE_DISPLAY_LEVELS))

p_heat <- ggplot(myeloid_time_mean, aes(x = timepoint, y = myeloid_state, fill = mean_weight)) + geom_tile(linewidth = 0.35) + 
    scale_fill_viridis_c(option = "magma") + labs(title = "RCTD-estimated mean weights of the 9 myeloid states", 
    subtitle = "Timepoint values are means of the two spatial-sample means", x = NULL, y = NULL, fill = "Mean\nweight") + 
    theme_classic(base_family = "Arial", base_size = 9) + theme(axis.text.y = element_text(size = 7.5), 
    plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8))

save_pdf_png(p = p_heat, stem = "08_65_all_myeloid_state_timepoint_heatmap", width = 7.8, height = 5.5)

spatial_data <- get_data_layer(spatial_obj, assay = "Spatial")

sanity_genes <- intersect(TARGET_SANITY_GENES, rownames(spatial_data))

sanity_cor <- data.frame()

if (length(sanity_genes) >= 3L) {
    sanity_score <- Matrix::colMeans(spatial_data[sanity_genes, , drop = FALSE])
    target_weight_named <- setNames(spot_weights[[TARGET_STATE]], spot_weights$barcode)
    sanity_rows <- list()
    for (samp in sample_order) {
        cells <- rownames(spatial_md)[as.character(spatial_md$sample_label) == samp]
        cells <- intersect(cells, names(target_weight_named))
        sanity_rows[[samp]] <- data.frame(sample_label = samp, timepoint = as.character(unique(spatial_md[cells, 
            "timepoint"])[1]), n_spots = length(cells), n_signature_genes = length(sanity_genes), signature_genes = paste(sanity_genes, 
            collapse = "; "), spearman_RCTD_vs_marker_score = safe_cor(target_weight_named[cells], sanity_score[cells]), 
            stringsAsFactors = FALSE)
    }
    sanity_cor <- bind_rows(sanity_rows)
    write_csv_utf8(sanity_cor, file.path(OUT_DIR, paste0("08b_65_target_marker_sanity_correlation_", 
        VERSION, ".csv")))
}

weight_type_cols <- setdiff(colnames(spot_weights), c("barcode", "GSM", "sample_label", "timepoint", 
    "spatial_region_v03", "array_x", "array_y"))

weight_column_map <- data.frame(RCTD_reference_type = weight_type_cols, Seurat_metadata_column = paste0("RCTD_", 
    make.names(weight_type_cols, unique = TRUE)), stringsAsFactors = FALSE)

write_csv_utf8(weight_column_map, file.path(OUT_DIR, paste0("03b_65_RCTD_weight_column_map_", VERSION, 
    ".csv")))

for (j in seq_len(nrow(weight_column_map))) {
    type_name <- weight_column_map$RCTD_reference_type[j]
    meta_name <- weight_column_map$Seurat_metadata_column[j]
    vals <- setNames(spot_weights[[type_name]], spot_weights$barcode)
    spatial_obj[[meta_name]] <- vals[Cells(spatial_obj)]
}

target_meta_col <- weight_column_map$Seurat_metadata_column[weight_column_map$RCTD_reference_type == 
    TARGET_STATE]

if (length(target_meta_col) != 1L) {
    stop("Could not uniquely identify target RCTD metadata column.")
}

spatial_obj@misc$GSE195783_step65 <- list(step = 65L, version = VERSION, method = "Bioconductor spacexr RCTD", 
    rctd_mode = RCTD_MODE, reference_design = "10 broad non-myeloid classes + 9 finalized myeloid states", 
    reference_max_cells_per_type = MAX_REF_CELLS_PER_TYPE, target_state = TARGET_STATE, target_metadata_column = target_meta_col, 
    per_sample_RCTD = TRUE, per_spot_inference = FALSE, histology_registered = FALSE)

ANNOTATED_OBJECT_FILE <- file.path(OUT_DIR, paste0("11_65_GSE195783_RCTD_annotated_spatial_", VERSION, 
    ".rds"))

saveRDS(spatial_obj, ANNOTATED_OBJECT_FILE, compress = TRUE)

all_returned <- all(sample_qc$fraction_spots_returned >= 0.95)

weight_sums_reasonable <- all(sample_qc$mean_weight_sum > 0.9 & sample_qc$mean_weight_sum < 1.1)

target_nonzero <- any(spot_weights[[TARGET_STATE]] > 0)

final_checks <- data.frame(check = c("Spatial input is Step64 v03 annotated object", "Hierarchical reference contains 19 reference types", 
    "Reference contains all 9 finalized myeloid states", "Reference contains 10 broad non-myeloid classes", 
    "All reference types contain at least 25 cells", "Reference and spatial data share at least 5000 genes", 
    "RCTD completed for all 8 spatial samples", "At least 95% of input spots were returned for every sample", 
    "Mean RCTD weight sum is approximately 1 for every sample", "Clec7a+ lipid-handling Microglia weight is present and non-zero", 
    "RCTD weights were added to all returned spatial spots", "No per-spot significance testing was performed", 
    "No H&E registration is claimed"), pass = c(grepl("64_GSE195783", SPATIAL_FILE), nrow(reference_manifest) == 
    19L, all(MYELOID_STATE_DISPLAY_LEVELS %in% reference_manifest$reference_type), all(NON_MYELOID_TYPES %in% 
    reference_manifest$reference_type), all(reference_manifest$n_cells >= MIN_REF_CELLS_PER_TYPE), length(common_genes) >= 
    5000L, nrow(sample_qc) == 8L, all_returned, weight_sums_reasonable, target_nonzero, all(!is.na(spatial_obj[[target_meta_col]][spot_weights$barcode, 
    1])), TRUE, all(spatial_obj$histology_registered == FALSE)), stringsAsFactors = FALSE)

write_csv_utf8(final_checks, file.path(OUT_DIR, paste0("09_65_GSE195783_RCTD_FINAL_QC_", VERSION, ".csv")))

reference_lines <- apply(reference_manifest, 1, function(z) {
    paste0("  ", z[["reference_type"]], ": ", z[["n_cells"]], " cells [", z[["reference_class"]], "]")
})

sample_lines <- apply(sample_qc, 1, function(z) {
    paste0("  ", z[["sample_label"]], ": ", z[["n_spots_RCTD"]], "/", z[["n_spots_input"]], " spots; mean weight sum=", 
        sprintf("%.4f", as.numeric(z[["mean_weight_sum"]])))
})

time_rows_report <- target_time_summary

target_time_lines <- apply(time_rows_report, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": mean of two sample means=", sprintf("%.4f", as.numeric(z[["mean_of_sample_means"]])), 
        " (range ", sprintf("%.4f", as.numeric(z[["min_sample_mean"]])), "-", sprintf("%.4f", as.numeric(z[["max_sample_mean"]])), 
        ")")
})

REPORT_FILE <- file.path(OUT_DIR, paste0("10_65_GSE195783_RCTD_REPORT_", VERSION, ".txt"))

report_con <- file(REPORT_FILE, open = "wt", encoding = "UTF-8")

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = report_con)

close(report_con)

print(target_time_summary)
print(final_checks)
