source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)
set.seed(111)

ROOT <- PROJECT_ROOT

STEP46_FILE <- file.path(ROOT, "results", "46_Bulk_limma_temporal_DEG_analysis", "30_46_Bulk_limma_temporal_DEG_analysis_v06_analysis_object.rds")

GLOBAL_REF_FILE <- file.path(ROOT, "results", "09_Figure1A_major_celltype_UMAP", "09_Figure1A_major_celltype_UMAP_object_v04.rds")

MYELOID_REF_FILE <- file.path(ROOT, "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

SIGNATURE_FILE <- file.path(ROOT, "results", "21_Clec7a_DEG_volcano_myeloid", "21_Clec7a_vs_all_myeloid_DEG_up_v01.csv")

STEP49_FILE <- file.path(ROOT, "results", "49_Bulk_WGCNA_network_modules", "09_49_Bulk_WGCNA_network_modules_v01_analysis_object.rds")

STEP44_BULK_CSV <- file.path(ROOT, "data", "17702045", "sci_bulkrna_expression.csv")

OUT_DIR <- file.path(ROOT, "results", "111_Bulk_myeloid_marker_score_adjustment")

TABLES_DIR <- file.path(OUT_DIR, "tables")
FIGURES_DIR <- file.path(OUT_DIR, "figures")
RDS_DIR <- file.path(OUT_DIR, "RDS")

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "111_Bulk_myeloid_marker_score_adjustment_v01"
LOG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".log"))

if (file.exists(LOG_FILE)) {
  file.remove(LOG_FILE)
}

write_csv <- function(x, filename, subdir = TABLES_DIR) {
  utils::write.csv(x, file = file.path(subdir, filename), row.names = FALSE, fileEncoding = "UTF-8",
    na = "")
}

EARLY_LABEL <- "Early"
LATE_LABEL <- "Established late"

GLOBAL_TYPES <- c("Neutrophil", "T/NK cell", "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Ependymal",
  "Endothelial", "Pericyte/VSMC", "Fibroblast/Meningeal")

NON_NEUTROPHIL_GLOBAL_TYPES <- setdiff(GLOBAL_TYPES, "Neutrophil")

MYELOID_STATE_RAW <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia",
  "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage",
  "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

MYELOID_STATE_DISPLAY <- c(`Homeostatic Microglia` = "Homeostatic Microglia", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling Microglia",
  `Activated inflammatory Microglia` = "Activated inflammatory Microglia", `Ccl12+ acute Microglia` = "Ccl12+ acute Microglia",
  `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory Microglia", `Proliferating Microglia` = "Proliferating Microglia",
  `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ antigen-presenting Macrophage", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb+ tissue-repair Macrophage",
  `Acute-response Macrophage` = "Acute-response Macrophage")

MYELOID_STATE_DISPLAY_LEVELS <- unname(MYELOID_STATE_DISPLAY[MYELOID_STATE_RAW])

CORE_GENES <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")

MAX_REF_CELLS_PER_TYPE <- 1000L
MIN_CELLS_PER_SOURCE_STATE <- 25L
MIN_CELLS_PER_TYPE <- 25L
MIN_STATE_LOGCPM <- 1
MIN_MYELOID_MEAN_LOGCPM <- 0.75
MIN_SPECIFICITY_LOG2 <- 0.75
MIN_STATE_PREVALENCE <- 7L
MAX_MARKERS_PER_SET <- 100L
MIN_MARKERS_FOR_SCORE <- 15L

ADJ_P_CUTOFF <- 0.05
ABS_LOGFC_CUTOFF <- 0.5

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "limma", "fgsea", "ggplot2")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs)) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall the missing package(s) before rerunning. ",
    "BayesPrism is NOT required by Step111.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(limma)
  library(fgsea)
  library(ggplot2)
})

safe_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- sub("^X(?=[0-9])", "", x, perl = TRUE)
  trimws(x)
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L)
    return(NA_real_)
  if (length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

safe_rbind <- function(x) {
    x <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, x)
    if (!length(x)) 
        return(data.frame())
    all_cols <- unique(unlist(lapply(x, names), use.names = FALSE))
    x <- lapply(x, function(z) {
        missing_cols <- setdiff(all_cols, names(z))
        for (nm in missing_cols) z[[nm]] <- NA
        z <- z[, all_cols, drop = FALSE]
        rownames(z) <- NULL
        z
    })
    out <- do.call(rbind, x)
    rownames(out) <- NULL
    out
}

assert_reference_counts <- function(mat, label) {
    if (is.null(dim(mat)) || !length(mat)) {
        stop(label, ": empty reference count matrix.")
    }
    values <- if (inherits(mat, "sparseMatrix")) {
        mat@x
    }
    else {
        as.numeric(mat)
    }
    if (!length(values) || any(!is.finite(values))) {
        stop(label, ": reference counts contain invalid or non-finite values.")
    }
    if (any(values < 0)) {
        stop(label, ": negative values detected in the scRNA reference counts.")
    }
    frac <- mean(abs(values - round(values)) < 1e-08)
    if (!is.finite(frac) || frac < 0.995) {
        stop(label, ": scRNA reference counts are not sufficiently integer-valued; ", "integer-like fraction=", 
            signif(frac, 5))
    }
    invisible(frac)
}

get_counts <- function(obj, assay = "RNA") {
    if (!assay %in% names(obj@assays)) {
        stop("Assay '", assay, "' was not found in the Seurat reference.")
    }
    obj2 <- obj
    if (inherits(obj2[[assay]], "Assay5")) {
        obj2 <- tryCatch(SeuratObject::JoinLayers(object = obj2, assay = assay), error = function(e) {
            stop("JoinLayers failed for reference object: ", conditionMessage(e))
        })
    }
    ans <- tryCatch(SeuratObject::LayerData(object = obj2, assay = assay, layer = "counts"), error = function(e) NULL)
    if (is.null(ans)) {
        ans <- tryCatch(Seurat::GetAssayData(object = obj2, assay = assay, slot = "counts"), error = function(e) {
            stop("Could not extract reference RNA counts: ", conditionMessage(e))
        })
    }
    if (!inherits(ans, "dgCMatrix")) {
        ans <- methods::as(ans, "dgCMatrix")
    }
    ans
}

collapse_sparse_gene_rows <- function(mat, label) {
    genes <- safe_gene_symbol(rownames(mat))
    keep <- !is.na(genes) & nzchar(genes)
    if (!all(keep)) {
        mat <- mat[keep, , drop = FALSE]
        genes <- genes[keep]
    }
    if (!length(genes)) {
        stop(label, ": no valid gene symbols remain.")
    }
    if (anyDuplicated(genes)) {
        lev <- unique(genes)
        map <- Matrix::sparseMatrix(i = match(genes, lev), j = seq_along(genes), x = 1, dims = c(length(lev), 
            length(genes)))
        mat <- map %*% mat
        rownames(mat) <- lev
        mat <- methods::as(mat, "dgCMatrix")
    }
    else {
        rownames(mat) <- genes
    }
    mat
}

balanced_downsample <- function(metadata, label_col, strata_col, max_per_label = 1000L, seed = 1234L) {
    set.seed(seed)
    labels <- as.character(metadata[[label_col]])
    strata <- as.character(metadata[[strata_col]])
    cells <- rownames(metadata)
    if (is.null(cells) || anyNA(labels) || anyNA(strata) || any(!nzchar(labels)) || any(!nzchar(strata))) {
        stop("Reference metadata contains missing cell IDs, labels or source_dataset.")
    }
    keep <- character()
    for (lab in sort(unique(labels))) {
        idx_lab <- which(labels == lab)
        if (length(idx_lab) <= max_per_label) {
            keep <- c(keep, cells[idx_lab])
            next
        }
        source_lab <- strata[idx_lab]
        source_levels <- sort(unique(source_lab))
        quota <- floor(max_per_label/length(source_levels))
        selected <- character()
        for (src in source_levels) {
            idx_src <- idx_lab[source_lab == src]
            n_take <- min(quota, length(idx_src))
            if (n_take > 0L) {
                selected <- c(selected, sample(cells[idx_src], n_take, replace = FALSE))
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

select_reference_cells <- function(metadata, label_col, expected_labels, reference_name, seed = 1234L) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    if (is.null(rownames(metadata))) {
        stop(reference_name, ": metadata must have cell IDs as row names.")
    }
    if (!all(c(label_col, "source_dataset") %in% colnames(metadata))) {
        stop(reference_name, ": metadata lacks ", label_col, " or source_dataset.")
    }
    metadata$profile_label <- as.character(metadata[[label_col]])
    metadata$source_dataset <- as.character(metadata$source_dataset)
    if (!all(expected_labels %in% unique(metadata$profile_label))) {
        stop(reference_name, ": missing expected cell labels: ", paste(setdiff(expected_labels, unique(metadata$profile_label)), 
            collapse = "; "))
    }
    metadata <- metadata[metadata$profile_label %in% expected_labels, , drop = FALSE]
    strata_key <- paste(metadata$profile_label, metadata$source_dataset, sep = "|||")
    strata_n <- table(strata_key)
    valid_strata <- names(strata_n)[strata_n >= MIN_CELLS_PER_SOURCE_STATE]
    keep_strata <- strata_key %in% valid_strata
    metadata <- metadata[keep_strata, , drop = FALSE]
    counts_by_type_source <- as.data.frame(table(profile_label = metadata$profile_label, source_dataset = metadata$source_dataset), 
        stringsAsFactors = FALSE)
    counts_by_type_source <- counts_by_type_source[counts_by_type_source$Freq > 0, , drop = FALSE]
    names(counts_by_type_source)[3] <- "n_cells_before_downsample"
    counts_by_type_source$reference <- reference_name
    type_counts <- table(metadata$profile_label)
    missing_after_filter <- setdiff(expected_labels, names(type_counts))
    if (length(missing_after_filter)) {
        stop(reference_name, ": labels lost after source-stratum QC: ", paste(missing_after_filter, collapse = "; "), 
            ". No incomplete reference will be forced.")
    }
    if (any(type_counts[expected_labels] < MIN_CELLS_PER_TYPE)) {
        stop(reference_name, ": at least one cell type has fewer than the minimum cells.")
    }
    selected_cells <- balanced_downsample(metadata = metadata, label_col = "profile_label", strata_col = "source_dataset", 
        max_per_label = MAX_REF_CELLS_PER_TYPE, seed = seed)
    selected_metadata <- metadata[selected_cells, , drop = FALSE]
    counts_by_type_retained <- as.data.frame(table(profile_label = selected_metadata$profile_label), 
        stringsAsFactors = FALSE)
    names(counts_by_type_retained)[2] <- "n_cells_retained"
    counts_by_type_retained$reference <- reference_name
    list(selected_cells = selected_cells, metadata = selected_metadata, cell_counts_by_type_source = counts_by_type_source, 
        cell_counts_by_type_retained = counts_by_type_retained)
}

make_reference_profiles <- function(counts, metadata, label_col, expected_labels, reference_name) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    if (is.null(rownames(metadata))) {
        stop(reference_name, ": metadata must have cell IDs as row names.")
    }
    if (!all(c(label_col, "source_dataset") %in% colnames(metadata))) {
        stop(reference_name, ": metadata lacks ", label_col, " or source_dataset.")
    }
    if (!all(colnames(counts) %in% rownames(metadata))) {
        stop(reference_name, ": count-matrix cells do not all match metadata.")
    }
    metadata <- metadata[colnames(counts), , drop = FALSE]
    metadata$profile_label <- as.character(metadata[[label_col]])
    if (!all(expected_labels %in% unique(metadata$profile_label))) {
        stop(reference_name, ": selected counts are missing expected labels.")
    }
    lib_size <- Matrix::colSums(counts)
    if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
        stop(reference_name, ": zero or invalid single-cell library size detected.")
    }
    norm <- counts %*% Matrix::Diagonal(x = 1e+06/as.numeric(lib_size))
    norm <- methods::as(norm, "dgCMatrix")
    norm@x <- log2(norm@x + 1)
    mean_logcpm <- lapply(expected_labels, function(lab) {
        idx <- which(metadata$profile_label == lab)
        if (length(idx) < MIN_CELLS_PER_TYPE) {
            stop(reference_name, ": insufficient retained cells for ", lab)
        }
        Matrix::rowMeans(norm[, idx, drop = FALSE])
    })
    mean_logcpm <- do.call(cbind, mean_logcpm)
    rownames(mean_logcpm) <- rownames(norm)
    colnames(mean_logcpm) <- expected_labels
    detection <- lapply(expected_labels, function(lab) {
        idx <- which(metadata$profile_label == lab)
        Matrix::rowMeans(counts[, idx, drop = FALSE] > 0)
    })
    detection <- do.call(cbind, detection)
    rownames(detection) <- rownames(counts)
    colnames(detection) <- expected_labels
    list(mean_logcpm = as.matrix(mean_logcpm), detection = as.matrix(detection))
}

derive_marker_set <- function(myeloid_profile, nonmyeloid_profile, signature_genes, marker_set_name, 
    min_state_prevalence = MIN_STATE_PREVALENCE) {
    common <- intersect(rownames(myeloid_profile), rownames(nonmyeloid_profile))
    if (length(common) < 500L) {
        stop(marker_set_name, ": fewer than 500 genes shared by reference profiles.")
    }
    my <- myeloid_profile[common, , drop = FALSE]
    nonmy <- nonmyeloid_profile[common, , drop = FALSE]
    my_mean <- rowMeans(my, na.rm = TRUE)
    nonmy_mean <- rowMeans(nonmy, na.rm = TRUE)
    specificity <- my_mean - nonmy_mean
    prevalence <- rowSums(my >= MIN_STATE_LOGCPM, na.rm = TRUE)
    in_signature <- common %in% signature_genes
    tab <- data.frame(marker_set = marker_set_name, gene = common, myeloid_mean_logCPM = as.numeric(my_mean), 
        nonmyeloid_mean_logCPM = as.numeric(nonmy_mean), specificity_log2 = as.numeric(specificity), 
        n_myeloid_states_above_expression_threshold = as.integer(prevalence), n_myeloid_states = ncol(my), 
        in_frozen_Clec7a_UP_signature = in_signature, stringsAsFactors = FALSE)
    pass <- (tab$myeloid_mean_logCPM >= MIN_MYELOID_MEAN_LOGCPM & tab$specificity_log2 >= MIN_SPECIFICITY_LOG2 & 
        tab$n_myeloid_states_above_expression_threshold >= min_state_prevalence & !tab$in_frozen_Clec7a_UP_signature)
    candidate <- tab[pass, , drop = FALSE]
    candidate <- candidate[order(-candidate$specificity_log2, -candidate$myeloid_mean_logCPM, candidate$gene), 
        , drop = FALSE]
    if (nrow(candidate) > MAX_MARKERS_PER_SET) {
        candidate <- candidate[seq_len(MAX_MARKERS_PER_SET), , drop = FALSE]
    }
    if (nrow(candidate) < MIN_MARKERS_FOR_SCORE) {
        stop(marker_set_name, ": only ", nrow(candidate), " marker genes passed the prespecified criteria; at least ", 
            MIN_MARKERS_FOR_SCORE, " are required. Review reference labels/QC rather than silently relaxing thresholds.")
    }
    tab$selected <- tab$gene %in% candidate$gene
    tab$selection_rank <- match(tab$gene, candidate$gene)
    list(ranking = tab, selected = candidate, genes = candidate$gene, profile_genes = common)
}

calculate_marker_score <- function(expr, genes, score_name) {
    genes_present <- intersect(genes, rownames(expr))
    if (length(genes_present) < MIN_MARKERS_FOR_SCORE) {
        stop(score_name, ": only ", length(genes_present), " selected marker genes overlap Step46 expression; at least ", 
            MIN_MARKERS_FOR_SCORE, " are required.")
    }
    x <- expr[genes_present, , drop = FALSE]
    gene_sd <- apply(x, 1, stats::sd, na.rm = TRUE)
    keep <- is.finite(gene_sd) & gene_sd > 0
    x <- x[keep, , drop = FALSE]
    if (nrow(x) < MIN_MARKERS_FOR_SCORE) {
        stop(score_name, ": fewer than ", MIN_MARKERS_FOR_SCORE, " variable marker genes remain in the Bulk matrix.")
    }
    gene_mean <- rowMeans(x, na.rm = TRUE)
    gene_sd <- apply(x, 1, stats::sd, na.rm = TRUE)
    z <- sweep(x, 1, gene_mean, FUN = "-")
    z <- sweep(z, 1, gene_sd, FUN = "/")
    score <- colMeans(z, na.rm = TRUE)
    list(score = as.numeric(score), sample_ids = colnames(expr), genes_used = rownames(x), n_selected = length(genes), 
        n_present = length(genes_present), n_variable = nrow(x), gene_z_matrix = z)
}

required_files <- c(STEP46_FILE, GLOBAL_REF_FILE, MYELOID_REF_FILE, SIGNATURE_FILE)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required input file(s) are missing:\n", paste(missing_files, collapse = "\n"))
}

step46 <- readRDS(STEP46_FILE)

required_step46 <- c("expression_original", "metadata", "primary")

if (!all(required_step46 %in% names(step46))) {
  stop("Step46 object is missing required components: ", paste(setdiff(required_step46, names(step46)),
    collapse = ", "))
}

expr <- as.matrix(step46$expression_original)
storage.mode(expr) <- "numeric"
meta <- as.data.frame(step46$metadata, stringsAsFactors = FALSE)

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Step46 expression_original must have gene and sample names.")
}
if (is.null(rownames(meta))) {
  stop("Step46 metadata must have sample IDs as row names.")
}
if (!all(colnames(expr) %in% rownames(meta))) {
  stop("Step46 expression columns do not all match metadata row names.")
}

meta <- meta[colnames(expr), , drop = FALSE]

if (!all(c("phase", "timepoint") %in% colnames(meta))) {
  stop("Step46 metadata must contain phase and timepoint.")
}

if (any(!is.finite(expr))) {
  stop("Step46 expression_original contains non-finite values.")
}

if (nrow(expr) != 12977L || ncol(expr) != 39L) {
  stop("Step46 matrix dimensions differ from the expected 12,977 genes x 39 samples: ", nrow(expr),
    " x ", ncol(expr), ". Stop and review the object rather than silently changing the analysis.")
}

gene_names <- safe_gene_symbol(rownames(expr))
if (anyNA(gene_names) || any(!nzchar(gene_names))) {
  stop("Blank/invalid gene symbols in Step46 expression_original.")
}
if (anyDuplicated(gene_names)) {
  dup <- unique(gene_names[duplicated(gene_names)])
  stop("Duplicate gene symbols in Step46 expression_original after symbol cleanup: ", paste(head(dup,
    20), collapse = "; "), ". No automatic aggregation is applied to the Bulk expression.")
}
rownames(expr) <- gene_names

primary_ids <- as.character(step46$primary$sample_ids)
if (length(primary_ids) < 10L || !all(primary_ids %in% colnames(expr))) {
  stop("Step46 primary sample IDs are missing or unexpectedly few.")
}
if (anyDuplicated(primary_ids)) {
  stop("Duplicate sample IDs in Step46 primary sample set.")
}

meta$sample_id <- rownames(meta)
meta$phase <- as.character(meta$phase)
meta$timepoint <- as.character(meta$timepoint)

if (!all(meta[primary_ids, "phase"] %in% c(EARLY_LABEL, LATE_LABEL))) {
  stop("Step46 primary sample IDs contain phase labels other than Early / Established late.")
}

primary_group <- factor(ifelse(meta[primary_ids, "phase"] == LATE_LABEL, "EstablishedLate", "Early"),
  levels = c("Early", "EstablishedLate"))

if (anyNA(primary_group) || length(unique(primary_group)) != 2L) {
  stop("Primary sample set must contain both Early and Established late groups.")
}

expr_primary <- expr[, primary_ids, drop = FALSE]

sig_raw <- utils::read.csv(SIGNATURE_FILE, check.names = FALSE, stringsAsFactors = FALSE)

required_sig_cols <- c("gene", "avg_log2FC", "p_val_adj")
if (!all(required_sig_cols %in% colnames(sig_raw))) {
  stop("Step21 signature file is missing gene, avg_log2FC or p_val_adj.")
}

sig_raw$gene <- safe_gene_symbol(sig_raw$gene)
sig_raw <- sig_raw[!is.na(sig_raw$gene) & nzchar(sig_raw$gene), , drop = FALSE]
sig_raw <- sig_raw[!duplicated(sig_raw$gene), , drop = FALSE]

if (nrow(sig_raw) != 1320L) {
  stop("Expected the frozen Step21 Clec7a-UP signature to contain 1,320 unique genes; found ", nrow(sig_raw),
    ". Do not silently change the signature definition.")
}
if (any(sig_raw$p_val_adj >= 0.05, na.rm = TRUE) || any(sig_raw$avg_log2FC < 0.25, na.rm = TRUE)) {
  stop("Step21 signature does not satisfy the frozen adjusted-P/logFC definition.")
}

signature_genes <- unique(sig_raw$gene)

input_audit <- data.frame(input = c(STEP46_FILE, GLOBAL_REF_FILE, MYELOID_REF_FILE, SIGNATURE_FILE, STEP49_FILE,
  STEP44_BULK_CSV), required = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), exists = file.exists(c(STEP46_FILE,
  GLOBAL_REF_FILE, MYELOID_REF_FILE, SIGNATURE_FILE, STEP49_FILE, STEP44_BULK_CSV)), used = c(TRUE,
  TRUE, TRUE, TRUE, FALSE, FALSE), role = c("continuous Step46 expression_original, metadata and primary sample IDs",
  "Whole-tissue reference for ten broad cell classes including Neutrophil", "Refined reference for nine microglia/macrophage states",
  "Frozen 1,320-gene Clec7a-UP signature; excluded from marker selection", "Optional WGCNA module eigengenes",
  "Not read as raw counts; original Bulk matrix is continuous transformed/normalized expression"),
  stringsAsFactors = FALSE)
write_csv(input_audit, "01_input_audit.csv")

bulk_value_audit <- data.frame(metric = c("n_genes", "n_samples", "minimum_expression", "maximum_expression",
  "mean_expression", "negative_value_count", "non_integer_like_fraction", "frozen_signature_n", "primary_sample_n",
  "primary_Early_n", "primary_EstablishedLate_n"), value = c(nrow(expr), ncol(expr), min(expr), max(expr),
  mean(expr), sum(expr < 0), mean(abs(as.numeric(expr) - round(as.numeric(expr))) >= 1e-08), length(signature_genes),
  length(primary_ids), sum(primary_group == "Early"), sum(primary_group == "EstablishedLate")), stringsAsFactors = FALSE)
write_csv(bulk_value_audit, "02_bulk_matrix_and_primary_sample_audit.csv")

global_ref <- readRDS(GLOBAL_REF_FILE)
myeloid_ref <- readRDS(MYELOID_REF_FILE)

if (!inherits(global_ref, "Seurat")) {
  stop("Step09 global reference is not a Seurat object.")
}
if (!inherits(myeloid_ref, "Seurat")) {
  stop("Step20 refined myeloid reference is not a Seurat object.")
}

global_md <- global_ref[[]]
myeloid_md <- myeloid_ref[[]]

if (!all(c("major_celltype", "source_dataset") %in% colnames(global_md))) {
  stop("Step09 global reference lacks major_celltype or source_dataset.")
}
if (!all(c("myeloid_state", "source_dataset") %in% colnames(myeloid_md))) {
  stop("Step20 refined myeloid reference lacks myeloid_state or source_dataset.")
}

global_md$major_celltype <- as.character(global_md$major_celltype)
global_md$source_dataset <- as.character(global_md$source_dataset)
myeloid_md$myeloid_state <- as.character(myeloid_md$myeloid_state)
myeloid_md$source_dataset <- as.character(myeloid_md$source_dataset)

observed_global <- sort(unique(global_md$major_celltype[global_md$major_celltype %in% GLOBAL_TYPES]))
missing_global <- setdiff(GLOBAL_TYPES, observed_global)
if (length(missing_global)) {
  stop("Step09 global reference is missing expected cell classes: ", paste(missing_global, collapse = "; "))
}

observed_myeloid <- sort(unique(myeloid_md$myeloid_state[myeloid_md$myeloid_state %in% MYELOID_STATE_RAW]))
missing_myeloid <- setdiff(MYELOID_STATE_RAW, observed_myeloid)
if (length(missing_myeloid)) {
  stop("Step20 reference is missing expected myeloid states: ", paste(missing_myeloid, collapse = "; "))
}

global_cells <- rownames(global_md)[global_md$major_celltype %in% GLOBAL_TYPES]
myeloid_cells <- rownames(myeloid_md)[myeloid_md$myeloid_state %in% MYELOID_STATE_RAW]

global_selection <- select_reference_cells(metadata = global_md[global_cells, , drop = FALSE], label_col = "major_celltype",
  expected_labels = GLOBAL_TYPES, reference_name = "Step09_global", seed = 11101L)

myeloid_md_for_selection <- myeloid_md[myeloid_cells, , drop = FALSE]
myeloid_md_for_selection$myeloid_state_display <- unname(MYELOID_STATE_DISPLAY[myeloid_md_for_selection$myeloid_state])

myeloid_selection <- select_reference_cells(metadata = myeloid_md_for_selection, label_col = "myeloid_state_display",
  expected_labels = MYELOID_STATE_DISPLAY_LEVELS, reference_name = "Step20_myeloid", seed = 11102L)

global_sub <- subset(global_ref, cells = global_selection$selected_cells)
myeloid_sub <- subset(myeloid_ref, cells = myeloid_selection$selected_cells)

global_counts <- get_counts(global_sub, assay = "RNA")
myeloid_counts <- get_counts(myeloid_sub, assay = "RNA")

global_counts <- collapse_sparse_gene_rows(global_counts, "Step09 global reference")
myeloid_counts <- collapse_sparse_gene_rows(myeloid_counts, "Step20 refined myeloid reference")

assert_reference_counts(global_counts, "Step09 global reference")
assert_reference_counts(myeloid_counts, "Step20 refined myeloid reference")

global_md_selected <- global_sub[[]]
myeloid_md_selected <- myeloid_sub[[]]
myeloid_md_selected$myeloid_state_display <- unname(MYELOID_STATE_DISPLAY[as.character(myeloid_md_selected$myeloid_state)])

reference_cell_counts <- safe_rbind(list(global_selection$cell_counts_by_type_source, myeloid_selection$cell_counts_by_type_source))
reference_retained_counts <- safe_rbind(list(global_selection$cell_counts_by_type_retained, myeloid_selection$cell_counts_by_type_retained))

write_csv(reference_cell_counts, "03_reference_cell_counts_by_type_source.csv")
write_csv(reference_retained_counts, "04_reference_cell_counts_retained.csv")

global_profiles <- make_reference_profiles(counts = global_counts, metadata = global_md_selected, label_col = "major_celltype",
  expected_labels = GLOBAL_TYPES, reference_name = "Step09_global")

myeloid_profiles <- make_reference_profiles(counts = myeloid_counts, metadata = myeloid_md_selected,
  label_col = "myeloid_state_display", expected_labels = MYELOID_STATE_DISPLAY_LEVELS, reference_name = "Step20_myeloid")

rm(global_ref, myeloid_ref, global_sub, myeloid_sub, global_md, myeloid_md, global_md_selected, myeloid_md_selected,
  global_counts, myeloid_counts, global_selection, myeloid_selection, myeloid_md_for_selection)
invisible(gc())

common_reference_genes <- intersect(rownames(global_profiles$mean_logcpm), rownames(myeloid_profiles$mean_logcpm))

if (length(common_reference_genes) < 500L) {
  stop("Fewer than 500 genes are shared between Step09 and Step20 reference profiles: ", length(common_reference_genes))
}

global_profile <- global_profiles$mean_logcpm[common_reference_genes, GLOBAL_TYPES, drop = FALSE]
myeloid_profile <- myeloid_profiles$mean_logcpm[common_reference_genes, MYELOID_STATE_DISPLAY_LEVELS,
  drop = FALSE]

myeloid9_marker <- derive_marker_set(myeloid_profile = myeloid_profile, nonmyeloid_profile = global_profile,
  signature_genes = signature_genes, marker_set_name = "Myeloid9", min_state_prevalence = MIN_STATE_PREVALENCE)

myeloid10_profile <- cbind(myeloid_profile, Neutrophil = global_profile[, "Neutrophil"])
nonmyeloid9_profile <- global_profile[, NON_NEUTROPHIL_GLOBAL_TYPES, drop = FALSE]

myeloid10_marker <- derive_marker_set(myeloid_profile = myeloid10_profile, nonmyeloid_profile = nonmyeloid9_profile,
  signature_genes = signature_genes, marker_set_name = "Myeloid10", min_state_prevalence = MIN_STATE_PREVALENCE)

marker_rankings <- safe_rbind(list(myeloid9_marker$ranking, myeloid10_marker$ranking))

marker_selected <- safe_rbind(list(myeloid9_marker$selected, myeloid10_marker$selected))

if (any(marker_selected$in_frozen_Clec7a_UP_signature)) {
  stop("A selected myeloid marker overlaps the frozen Clec7a-UP signature.")
}

write_csv(marker_rankings, "05_myeloid_marker_gene_rankings.csv")
write_csv(marker_selected, "06_selected_myeloid_marker_genes.csv")

marker_set_summary <- data.frame(marker_set = c("Myeloid9", "Myeloid10"), n_reference_genes_tested = c(nrow(myeloid9_marker$ranking),
  nrow(myeloid10_marker$ranking)), n_marker_genes_selected = c(length(myeloid9_marker$genes), length(myeloid10_marker$genes)),
  n_marker_genes_in_bulk = c(length(intersect(myeloid9_marker$genes, rownames(expr))), length(intersect(myeloid10_marker$genes,
    rownames(expr)))), n_marker_genes_overlapping_frozen_signature = c(length(intersect(myeloid9_marker$genes,
    signature_genes)), length(intersect(myeloid10_marker$genes, signature_genes))), min_state_prevalence = MIN_STATE_PREVALENCE,
  min_myeloid_mean_logCPM = MIN_MYELOID_MEAN_LOGCPM, min_specificity_log2 = MIN_SPECIFICITY_LOG2, max_markers = MAX_MARKERS_PER_SET,
  stringsAsFactors = FALSE)
write_csv(marker_set_summary, "07_marker_set_summary.csv")

score9 <- calculate_marker_score(expr = expr, genes = myeloid9_marker$genes, score_name = "Myeloid9_marker_score")

score10 <- calculate_marker_score(expr = expr, genes = myeloid10_marker$genes, score_name = "Myeloid10_marker_score")

if (!identical(score9$sample_ids, colnames(expr)) || !identical(score10$sample_ids, colnames(expr))) {
  stop("Marker score sample order does not match Step46 expression columns.")
}

score_all <- data.frame(sample_id = colnames(expr), timepoint = as.character(meta[colnames(expr), "timepoint"]),
  phase = as.character(meta[colnames(expr), "phase"]), Myeloid9_marker_score = score9$score, Myeloid10_marker_score = score10$score,
  stringsAsFactors = FALSE)

score_all$primary_group <- NA_character_
score_all$primary_group[score_all$sample_id %in% primary_ids] <- as.character(primary_group)

write_csv(score_all, "08_bulk_myeloid_marker_scores_all_samples.csv")

score_primary <- score_all[match(primary_ids, score_all$sample_id), , drop = FALSE]

if (!identical(score_primary$sample_id, primary_ids)) {
  stop("Primary marker score table does not match Step46 primary sample order.")
}

score_summary_rows <- list()
for (score_name in c("Myeloid9_marker_score", "Myeloid10_marker_score")) {
  for (grp in levels(primary_group)) {
    x <- score_primary[[score_name]][score_primary$primary_group == grp]
    score_summary_rows[[length(score_summary_rows) + 1L]] <- data.frame(score_name = score_name,
      primary_group = grp, n = length(x), mean = mean(x, na.rm = TRUE), SD = stats::sd(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE), IQR = stats::IQR(x, na.rm = TRUE), stringsAsFactors = FALSE)
  }
}
score_group_summary <- safe_rbind(score_summary_rows)
write_csv(score_group_summary, "09_primary_group_marker_score_summary.csv")

score_phase_test <- function(score, score_name) {
  d <- data.frame(score = as.numeric(score), primary_group = primary_group)
  fit <- stats::lm(score ~ primary_group, data = d)
  sm <- summary(fit)$coefficients
  term <- "primary_groupEstablishedLate"

  data.frame(score_name = score_name, n = nrow(d), mean_Early = mean(d$score[d$primary_group == "Early"],
    na.rm = TRUE), mean_EstablishedLate = mean(d$score[d$primary_group == "EstablishedLate"], na.rm = TRUE),
    difference_late_minus_early = mean(d$score[d$primary_group == "EstablishedLate"], na.rm = TRUE) -
      mean(d$score[d$primary_group == "Early"], na.rm = TRUE), estimate = if (term %in% rownames(sm))
      sm[term, "Estimate"] else NA_real_, SE = if (term %in% rownames(sm))
      sm[term, "Std. Error"] else NA_real_, t = if (term %in% rownames(sm))
      sm[term, "t value"] else NA_real_, P.Value = if (term %in% rownames(sm))
      sm[term, "Pr(>|t|)"] else NA_real_, stringsAsFactors = FALSE)
}

score_phase_tests <- safe_rbind(list(score_phase_test(score_primary$Myeloid9_marker_score, "Myeloid9_marker_score"),
  score_phase_test(score_primary$Myeloid10_marker_score, "Myeloid10_marker_score")))
score_phase_tests$BH_FDR <- p.adjust(score_phase_tests$P.Value, method = "BH")
write_csv(score_phase_tests, "10_primary_marker_score_phase_tests.csv")

design_unadjusted <- stats::model.matrix(~0 + primary_group, data = data.frame(primary_group = primary_group))
group_cols <- grep("^primary_group", colnames(design_unadjusted))
if (length(group_cols) != 2L) {
  stop("Could not identify the two phase columns in the unadjusted design.")
}
colnames(design_unadjusted)[group_cols] <- c("Early", "EstablishedLate")

if (qr(design_unadjusted)$rank != ncol(design_unadjusted)) {
  stop("Unadjusted phase-only design is rank deficient.")
}

array_weights <- limma::arrayWeights(expr_primary, design = design_unadjusted)

if (length(array_weights) != ncol(expr_primary) || any(!is.finite(array_weights)) || any(array_weights <=
  0)) {
  stop("limma::arrayWeights returned invalid sample weights.")
}

make_contrast <- function(design) {
  contrast <- matrix(0, nrow = ncol(design), ncol = 1, dimnames = list(colnames(design), "EstablishedLate_vs_Early"))
  contrast["EstablishedLate", 1] <- 1
  contrast["Early", 1] <- -1
  contrast
}

fit_limma_model <- function(analysis_id, model_label, marker_score = NULL, marker_score_name = NA_character_) {
    d <- data.frame(primary_group = primary_group, stringsAsFactors = FALSE)
    if (is.null(marker_score)) {
        design <- stats::model.matrix(~0 + primary_group, data = d)
    }
    else {
        marker_score <- as.numeric(marker_score)
        if (length(marker_score) != nrow(d) || any(!is.finite(marker_score)) || stats::sd(marker_score) == 
            0) {
            return(list(table = NULL, diagnostics = data.frame(analysis_id = analysis_id, model = model_label, 
                adjustment = marker_score_name, status = "NOT_ESTIMABLE_INVALID_OR_CONSTANT_SCORE", n_samples = nrow(d), 
                n_parameters = NA_integer_, design_rank = NA_integer_, residual_df = NA_integer_, score_phase_correlation = NA_real_, 
                VIF_marker_score = NA_real_, condition_number = NA_real_, reason = "Marker score is missing, non-finite or has no variation.", 
                stringsAsFactors = FALSE)))
        }
        d$marker_score_z <- as.numeric(scale(marker_score))
        design <- stats::model.matrix(~0 + primary_group + marker_score_z, data = d)
    }
    group_cols <- grep("^primary_group", colnames(design))
    if (length(group_cols) != 2L) {
        stop("Could not identify two primary phase coefficients in ", analysis_id)
    }
    colnames(design)[group_cols] <- c("Early", "EstablishedLate")
    design_rank <- qr(design)$rank
    residual_df <- nrow(design) - design_rank
    condition_number <- tryCatch(kappa(design, exact = TRUE), error = function(e) NA_real_)
    score_phase_correlation <- NA_real_
    vif_marker_score <- NA_real_
    if (!is.null(marker_score)) {
        score_phase_correlation <- safe_cor(as.numeric(marker_score), as.numeric(primary_group == "EstablishedLate"), 
            method = "pearson")
        aux <- stats::lm(d$marker_score_z ~ d$primary_group)
        r2 <- summary(aux)$r.squared
        vif_marker_score <- if (is.finite(r2) && r2 < 1) {
            1/(1 - r2)
        }
        else {
            Inf
        }
    }
    if (design_rank != ncol(design) || residual_df < 2L) {
        return(list(table = NULL, diagnostics = data.frame(analysis_id = analysis_id, model = model_label, 
            adjustment = marker_score_name, status = if (design_rank != ncol(design)) {
                "NOT_ESTIMABLE_RANK_DEFICIENT"
            } else {
                "NOT_ESTIMABLE_LOW_RESIDUAL_DF"
            }, n_samples = nrow(d), n_parameters = ncol(design), design_rank = design_rank, residual_df = residual_df, 
            score_phase_correlation = score_phase_correlation, VIF_marker_score = vif_marker_score, condition_number = condition_number, 
            reason = "The design was not safely estimable; no model was forced.", stringsAsFactors = FALSE)))
    }
    fit0 <- limma::lmFit(expr_primary, design = design, weights = array_weights)
    fit1 <- limma::contrasts.fit(fit0, contrasts = make_contrast(design))
    fit1 <- limma::eBayes(fit1, trend = TRUE, robust = TRUE)
    tt <- limma::topTable(fit1, coef = "EstablishedLate_vs_Early", number = Inf, sort.by = "P", adjust.method = "BH")
    tt$gene_id <- safe_gene_symbol(rownames(tt))
    tt$analysis_id <- analysis_id
    tt$model <- model_label
    tt$direction <- "NS"
    tt$direction[tt$adj.P.Val < ADJ_P_CUTOFF & tt$logFC >= ABS_LOGFC_CUTOFF] <- "Up"
    tt$direction[tt$adj.P.Val < ADJ_P_CUTOFF & tt$logFC <= -ABS_LOGFC_CUTOFF] <- "Down"
    tt <- tt[, c("analysis_id", "model", "gene_id", setdiff(colnames(tt), c("analysis_id", "model", "gene_id"))), 
        drop = FALSE]
    reason <- if (is.finite(vif_marker_score) && vif_marker_score > 5) {
        "Full-rank model fitted; marker score is strongly associated with phase (VIF > 5). Interpret cautiously."
    }
    else {
        "Full-rank limma model fitted; no batch term was used."
    }
    diagnostics <- data.frame(analysis_id = analysis_id, model = model_label, adjustment = marker_score_name, 
        status = "FITTED", n_samples = nrow(d), n_parameters = ncol(design), design_rank = design_rank, 
        residual_df = residual_df, score_phase_correlation = score_phase_correlation, VIF_marker_score = vif_marker_score, 
        condition_number = condition_number, reason = reason, stringsAsFactors = FALSE)
    list(table = tt, diagnostics = diagnostics, design = design, fit = fit1)
}

model_unadjusted <- fit_limma_model(analysis_id = "UNADJUSTED_PHASE_ONLY", model_label = "Primary phase-only model; no batch term",
  marker_score = NULL)

model_myeloid9 <- fit_limma_model(analysis_id = "ADJ_MYELOID9_MARKER_SCORE", model_label = "Primary model adjusted for Myeloid9 marker score",
  marker_score = score_primary$Myeloid9_marker_score, marker_score_name = "Myeloid9_marker_score")

model_myeloid10 <- fit_limma_model(analysis_id = "ADJ_MYELOID10_MARKER_SCORE", model_label = "Primary model adjusted for Myeloid10 marker score",
  marker_score = score_primary$Myeloid10_marker_score, marker_score_name = "Myeloid10_marker_score")

model_results <- list(UNADJUSTED_PHASE_ONLY = model_unadjusted, ADJ_MYELOID9_MARKER_SCORE = model_myeloid9,
  ADJ_MYELOID10_MARKER_SCORE = model_myeloid10)

model_diagnostics <- safe_rbind(lapply(model_results, `[[`, "diagnostics"))
write_csv(model_diagnostics, "11_limma_model_diagnostics.csv")

array_weight_table <- data.frame(sample_id = primary_ids, primary_group = as.character(primary_group),
  Myeloid9_marker_score = score_primary$Myeloid9_marker_score, Myeloid10_marker_score = score_primary$Myeloid10_marker_score,
  array_weight = as.numeric(array_weights), stringsAsFactors = FALSE)
write_csv(array_weight_table, "12_limma_array_weights.csv")

gene_tables <- lapply(model_results, function(z) z$table)
names(gene_tables) <- names(model_results)
gene_tables_nonnull <- Filter(function(z) !is.null(z), gene_tables)
gene_level_results <- safe_rbind(gene_tables_nonnull)

if (!nrow(gene_level_results)) {
  stop("No limma gene-level result tables were generated.")
}
write_csv(gene_level_results, "13_limma_primary_gene_level_results.csv")

tt_unadjusted <- model_unadjusted$table
if (is.null(tt_unadjusted)) {
  stop("The unadjusted primary phase-only limma model was not estimable.")
}

collapse_rank_table <- function(tt, label) {
    if (is.null(tt) || !"t" %in% colnames(tt)) {
        stop(label, " has no moderated t statistic.")
    }
    gene_ids <- if ("gene_id" %in% colnames(tt)) {
        safe_gene_symbol(tt$gene_id)
    }
    else {
        safe_gene_symbol(rownames(tt))
    }
    out <- data.frame(gene = gene_ids, t = as.numeric(tt$t), logFC = as.numeric(tt$logFC), P.Value = as.numeric(tt$P.Value), 
        adj.P.Val = as.numeric(tt$adj.P.Val), stringsAsFactors = FALSE)
    out <- out[!is.na(out$gene) & nzchar(out$gene) & is.finite(out$t), , drop = FALSE]
    out <- out[order(abs(out$t), decreasing = TRUE), , drop = FALSE]
    out <- out[!duplicated(out$gene), , drop = FALSE]
    out <- out[order(out$t, decreasing = TRUE), , drop = FALSE]
    rownames(out) <- NULL
    out
}

run_signature_gsea <- function(tt, analysis_id, label) {
    if (is.null(tt)) 
        return(NULL)
    rank_table <- collapse_rank_table(tt, label)
    ranks <- rank_table$t
    names(ranks) <- rank_table$gene
    ranks <- sort(ranks, decreasing = TRUE)
    overlap <- intersect(signature_genes, names(ranks))
    if (length(overlap) < 15L) {
        stop("Insufficient overlap between frozen Step21 signature and ranked Bulk genes for ", analysis_id, 
            ": ", length(overlap))
    }
    set.seed(111)
    fg <- fgsea::fgseaMultilevel(pathways = list(Clec7a_UP_scRNA = overlap), stats = ranks, minSize = 15, 
        maxSize = 5000, eps = 0, scoreType = "std")
    leading <- fg$leadingEdge[[1]]
    summary <- data.frame(analysis_id = analysis_id, analysis = label, NES = as.numeric(fg$NES[1]), ES = as.numeric(fg$ES[1]), 
        pval = as.numeric(fg$pval[1]), padj_within_model = as.numeric(fg$padj[1]), pathway_size = as.integer(fg$size[1]), 
        signature_n = length(signature_genes), signature_n_in_ranked_bulk = length(overlap), signature_overlap_fraction = length(overlap)/length(signature_genes), 
        leading_edge_n = length(leading), leading_edge_genes = paste(leading, collapse = ";"), stringsAsFactors = FALSE)
    leading_table <- rank_table[match(leading, rank_table$gene), , drop = FALSE]
    leading_table$analysis_id <- analysis_id
    leading_table$analysis <- label
    leading_table$leading_edge_order <- seq_len(nrow(leading_table))
    leading_table <- leading_table[, c("analysis_id", "analysis", "leading_edge_order", "gene", "logFC", 
        "t", "P.Value", "adj.P.Val"), drop = FALSE]
    list(summary = summary, leading_edge = leading_table, rank_table = rank_table)
}

gsea_results <- list()
gsea_results[["UNADJUSTED_PHASE_ONLY"]] <- run_signature_gsea(tt_unadjusted, "UNADJUSTED_PHASE_ONLY",
  "Primary phase-only model; no batch term")

if (!is.null(model_myeloid9$table)) {
  gsea_results[["ADJ_MYELOID9_MARKER_SCORE"]] <- run_signature_gsea(model_myeloid9$table, "ADJ_MYELOID9_MARKER_SCORE",
    "Primary model adjusted for Myeloid9 marker score")
}

if (!is.null(model_myeloid10$table)) {
  gsea_results[["ADJ_MYELOID10_MARKER_SCORE"]] <- run_signature_gsea(model_myeloid10$table, "ADJ_MYELOID10_MARKER_SCORE",
    "Primary model adjusted for Myeloid10 marker score")
}

gsea_summary <- safe_rbind(lapply(gsea_results, `[[`, "summary"))
if (nrow(gsea_summary)) {
    gsea_summary$BH_FDR_across_planned_models <- p.adjust(gsea_summary$pval, method = "BH")
    unadj_nes <- gsea_summary$NES[gsea_summary$analysis_id == "UNADJUSTED_PHASE_ONLY"]
    gsea_summary$NES_attenuation_percent_vs_unadjusted <- NA_real_
    if (length(unadj_nes) == 1L && is.finite(unadj_nes) && unadj_nes != 0) {
        gsea_summary$NES_attenuation_percent_vs_unadjusted <- 100 * (1 - abs(gsea_summary$NES/unadj_nes))
    }
}
write_csv(gsea_summary, "14_Clec7a_signature_GSEA_unadjusted_vs_adjusted.csv")

leading_edge_all <- safe_rbind(lapply(gsea_results, `[[`, "leading_edge"))
if (nrow(leading_edge_all)) {
  write_csv(leading_edge_all, "15_Clec7a_signature_GSEA_leading_edge_genes.csv")
}

compare_gene_effects <- function(base_tt, adj_tt, analysis_id) {
    if (is.null(adj_tt)) 
        return(data.frame())
    base_rank <- collapse_rank_table(base_tt, "Unadjusted")
    adj_rank <- collapse_rank_table(adj_tt, analysis_id)
    base_use <- base_rank[, c("gene", "logFC", "t", "P.Value", "adj.P.Val")]
    adj_use <- adj_rank[, c("gene", "logFC", "t", "P.Value", "adj.P.Val")]
    names(base_use)[-1] <- paste0(names(base_use)[-1], "_unadjusted")
    names(adj_use)[-1] <- paste0(names(adj_use)[-1], "_adjusted")
    m <- merge(base_use, adj_use, by = "gene", all = FALSE)
    sig_base <- m$adj.P.Val_unadjusted < ADJ_P_CUTOFF & abs(m$logFC_unadjusted) >= ABS_LOGFC_CUTOFF
    sig_adj <- m$adj.P.Val_adjusted < ADJ_P_CUTOFF & abs(m$logFC_adjusted) >= ABS_LOGFC_CUTOFF
    data.frame(analysis_id = analysis_id, n_common_genes = nrow(m), pearson_logFC = safe_cor(m$logFC_unadjusted, 
        m$logFC_adjusted, "pearson"), spearman_logFC = safe_cor(m$logFC_unadjusted, m$logFC_adjusted, 
        "spearman"), pearson_moderated_t = safe_cor(m$t_unadjusted, m$t_adjusted, "pearson"), spearman_moderated_t = safe_cor(m$t_unadjusted, 
        m$t_adjusted, "spearman"), same_direction_fraction = mean(sign(m$logFC_unadjusted) == sign(m$logFC_adjusted), 
        na.rm = TRUE), unadjusted_effect_sig_n = sum(sig_base, na.rm = TRUE), adjusted_effect_sig_n = sum(sig_adj, 
        na.rm = TRUE), shared_effect_sig_n = sum(sig_base & sig_adj, na.rm = TRUE), stringsAsFactors = FALSE)
}

gene_effect_comparison <- safe_rbind(list(compare_gene_effects(tt_unadjusted, model_myeloid9$table, "ADJ_MYELOID9_MARKER_SCORE"),
  compare_gene_effects(tt_unadjusted, model_myeloid10$table, "ADJ_MYELOID10_MARKER_SCORE")))
write_csv(gene_effect_comparison, "16_gene_effect_concordance_vs_unadjusted.csv")

count_one <- function(tt, analysis_id, label) {
    if (is.null(tt)) {
        return(data.frame(analysis_id = analysis_id, analysis = label, n_genes_tested = NA_integer_, 
            n_adjP_lt_0.05 = NA_integer_, n_up_adjP_lt_0.05_logFC_ge_0.5 = NA_integer_, n_down_adjP_lt_0.05_logFC_le_minus0.5 = NA_integer_, 
            status = "NOT_ESTIMABLE", stringsAsFactors = FALSE))
    }
    data.frame(analysis_id = analysis_id, analysis = label, n_genes_tested = nrow(tt), n_adjP_lt_0.05 = sum(tt$adj.P.Val < 
        ADJ_P_CUTOFF, na.rm = TRUE), n_up_adjP_lt_0.05_logFC_ge_0.5 = sum(tt$adj.P.Val < ADJ_P_CUTOFF & 
        tt$logFC >= ABS_LOGFC_CUTOFF, na.rm = TRUE), n_down_adjP_lt_0.05_logFC_le_minus0.5 = sum(tt$adj.P.Val < 
        ADJ_P_CUTOFF & tt$logFC <= -ABS_LOGFC_CUTOFF, na.rm = TRUE), status = "FITTED", stringsAsFactors = FALSE)
}

deg_count_summary <- safe_rbind(list(count_one(tt_unadjusted, "UNADJUSTED_PHASE_ONLY", "Primary phase-only model; no batch term"),
  count_one(model_myeloid9$table, "ADJ_MYELOID9_MARKER_SCORE", "Adjusted for Myeloid9 marker score"),
  count_one(model_myeloid10$table, "ADJ_MYELOID10_MARKER_SCORE", "Adjusted for Myeloid10 marker score")))
write_csv(deg_count_summary, "17_DEG_count_summary_unadjusted_vs_adjusted.csv")

core_gene_rows <- list()
for (model_id in names(model_results)) {
    tt <- model_results[[model_id]]$table
    if (is.null(tt)) 
        next
    tt2 <- as.data.frame(tt)
    tt2$gene_id <- safe_gene_symbol(tt2$gene_id)
    for (gene in CORE_GENES) {
        idx <- match(gene, tt2$gene_id)
        if (is.na(idx)) {
            core_gene_rows[[length(core_gene_rows) + 1L]] <- data.frame(analysis_id = model_id, gene = gene, 
                found = FALSE, logFC = NA_real_, t = NA_real_, P.Value = NA_real_, adj.P.Val = NA_real_, 
                stringsAsFactors = FALSE)
        }
        else {
            z <- tt2[idx, , drop = FALSE]
            core_gene_rows[[length(core_gene_rows) + 1L]] <- data.frame(analysis_id = model_id, gene = gene, 
                found = TRUE, logFC = as.numeric(z$logFC[1]), t = as.numeric(z$t[1]), P.Value = as.numeric(z$P.Value[1]), 
                adj.P.Val = as.numeric(z$adj.P.Val[1]), stringsAsFactors = FALSE)
        }
    }
}
core_gene_summary <- safe_rbind(core_gene_rows)
write_csv(core_gene_summary, "18_core_gene_effects_unadjusted_vs_adjusted.csv")

signature_in_expr <- intersect(signature_genes, rownames(expr))
sig_expr <- expr[signature_in_expr, , drop = FALSE]
sig_sd <- apply(sig_expr, 1, stats::sd, na.rm = TRUE)
sig_expr <- sig_expr[is.finite(sig_sd) & sig_sd > 0, , drop = FALSE]

signature_score_status <- "SKIPPED_LT15_VARIABLE_SIGNATURE_GENES"
signature_score_table <- data.frame()
signature_score_models <- data.frame()

if (nrow(sig_expr) >= MIN_MARKERS_FOR_SCORE) {
    sig_mean <- rowMeans(sig_expr, na.rm = TRUE)
    sig_sd <- apply(sig_expr, 1, stats::sd, na.rm = TRUE)
    sig_z <- sweep(sig_expr, 1, sig_mean, FUN = "-")
    sig_z <- sweep(sig_z, 1, sig_sd, FUN = "/")
    signature_score_all <- colMeans(sig_z, na.rm = TRUE)
    signature_score_table <- data.frame(sample_id = primary_ids, primary_group = factor(as.character(primary_group), 
        levels = c("Early", "EstablishedLate")), Clec7a_UP_signature_score = as.numeric(signature_score_all[primary_ids]), 
        Myeloid9_marker_score = score_primary$Myeloid9_marker_score, Myeloid10_marker_score = score_primary$Myeloid10_marker_score, 
        stringsAsFactors = FALSE)
    fit_signature_score <- function(score_covariate, model_id) {
        d <- signature_score_table
        if (is.null(score_covariate)) {
            fit <- stats::lm(Clec7a_UP_signature_score ~ primary_group, data = d)
        }
        else {
            d$score_z <- as.numeric(scale(score_covariate))
            fit <- stats::lm(Clec7a_UP_signature_score ~ primary_group + score_z, data = d)
        }
        sm <- summary(fit)$coefficients
        term <- "primary_groupEstablishedLate"
        data.frame(analysis_id = model_id, estimate_late_vs_early = if (term %in% rownames(sm)) 
            sm[term, "Estimate"]
        else NA_real_, SE = if (term %in% rownames(sm)) 
            sm[term, "Std. Error"]
        else NA_real_, t = if (term %in% rownames(sm)) 
            sm[term, "t value"]
        else NA_real_, P.Value = if (term %in% rownames(sm)) 
            sm[term, "Pr(>|t|)"]
        else NA_real_, stringsAsFactors = FALSE)
    }
    signature_score_models <- safe_rbind(list(fit_signature_score(NULL, "UNADJUSTED_PHASE_ONLY"), fit_signature_score(signature_score_table$Myeloid9_marker_score, 
        "ADJ_MYELOID9_MARKER_SCORE"), fit_signature_score(signature_score_table$Myeloid10_marker_score, 
        "ADJ_MYELOID10_MARKER_SCORE")))
    write_csv(signature_score_table, "19_primary_sample_Clec7a_signature_score.csv")
    write_csv(signature_score_models, "20_Clec7a_signature_score_phase_effects.csv")
    signature_score_status <- "COMPLETED"
}

module_adjusted_results <- data.frame()
module_analysis_status <- "SKIPPED_STEP49_OBJECT_NOT_FOUND"

if (file.exists(STEP49_FILE)) {
  step49 <- readRDS(STEP49_FILE)

  if (!is.null(step49$module_eigengenes)) {
    MEs <- as.matrix(step49$module_eigengenes)

    if (all(primary_ids %in% rownames(MEs))) {
      MEs_primary <- MEs[primary_ids, , drop = FALSE]
    } else if (all(primary_ids %in% colnames(MEs))) {
      MEs_primary <- t(MEs[, primary_ids, drop = FALSE])
    } else {
      MEs_primary <- NULL
    }

    if (!is.null(MEs_primary)) {
      module_rows <- lapply(seq_len(ncol(MEs_primary)), function(j) {
        d <- data.frame(primary_group = primary_group, ME_value = as.numeric(MEs_primary[, j]),
          Myeloid9_z = as.numeric(scale(score_primary$Myeloid9_marker_score)), stringsAsFactors = FALSE)

        fit <- stats::lm(ME_value ~ primary_group + Myeloid9_z, data = d)
        sm <- summary(fit)$coefficients
        term <- "primary_groupEstablishedLate"

        data.frame(module = colnames(MEs_primary)[j], estimate_late_vs_early = if (term %in%
          rownames(sm))
          sm[term, "Estimate"] else NA_real_, SE = if (term %in% rownames(sm))
          sm[term, "Std. Error"] else NA_real_, t = if (term %in% rownames(sm))
          sm[term, "t value"] else NA_real_, P.Value = if (term %in% rownames(sm))
          sm[term, "Pr(>|t|)"] else NA_real_, stringsAsFactors = FALSE)
      })

      module_adjusted_results <- safe_rbind(module_rows)
      module_adjusted_results$BH_FDR <- p.adjust(module_adjusted_results$P.Value, method = "BH")
      write_csv(module_adjusted_results, "21_WGCNA_module_eigengene_phase_effects_adjusted_for_Myeloid9_score.csv")
      module_analysis_status <- "COMPLETED"
    } else {
      module_analysis_status <- "SKIPPED_STEP49_EIGENGENE_SAMPLE_IDS_NOT_MATCHED"
    }
  } else {
    module_analysis_status <- "SKIPPED_STEP49_MODULE_EIGENGENES_NOT_FOUND"
  }
}

plot_df <- score_primary
plot_df$primary_group <- factor(plot_df$primary_group, levels = c("Early", "EstablishedLate"), labels = c("Early",
  "Established late"))

score_long <- rbind(data.frame(sample_id = plot_df$sample_id, primary_group = plot_df$primary_group,
  score_name = "Myeloid9 marker score", score = plot_df$Myeloid9_marker_score), data.frame(sample_id = plot_df$sample_id,
  primary_group = plot_df$primary_group, score_name = "Myeloid10 marker score", score = plot_df$Myeloid10_marker_score))

p_score <- ggplot(score_long, aes(x = primary_group, y = score)) + geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, size = 2) + facet_wrap(~score_name, scales = "free_y") + labs(title = "Reference-derived myeloid marker scores",
  x = NULL, y = "Mean gene-wise z-score") + theme_classic(base_size = 11)

ggsave(file.path(FIGURES_DIR, paste0(PREFIX, "_myeloid_marker_scores_by_phase.pdf")), p_score, width = 8,
  height = 5, device = grDevices::cairo_pdf)
ggsave(file.path(FIGURES_DIR, paste0(PREFIX, "_myeloid_marker_scores_by_phase.png")), p_score, width = 8,
  height = 5, dpi = 300)

if (nrow(gsea_summary)) {
    p_gsea <- ggplot(gsea_summary, aes(x = analysis_id, y = NES)) + geom_col(width = 0.65) + geom_hline(yintercept = 0, 
        linetype = 2) + labs(title = "Clec7a-UP signature enrichment before and after marker-score adjustment", 
        x = NULL, y = "Normalized enrichment score (NES)") + theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 25, 
        hjust = 1))
    ggsave(file.path(FIGURES_DIR, paste0(PREFIX, "_GSEA_NES_comparison.pdf")), p_gsea, width = 8, height = 5, 
        device = grDevices::cairo_pdf)
    ggsave(file.path(FIGURES_DIR, paste0(PREFIX, "_GSEA_NES_comparison.png")), p_gsea, width = 8, height = 5, 
        dpi = 300)
}

analysis_object <- list(step = 111L, version = "v01", purpose = paste("Test whether bulk Clec7a-UP signature enrichment is robust to",
  "adjustment for independent scRNA-reference-derived myeloid marker scores."), input_files = input_audit,
  bulk_value_audit = bulk_value_audit, primary_sample_ids = primary_ids, primary_group = data.frame(sample_id = primary_ids,
    primary_group = as.character(primary_group), stringsAsFactors = FALSE), reference_cell_counts = reference_cell_counts,
  reference_retained_counts = reference_retained_counts, reference_profiles = list(global_mean_logcpm = global_profile,
    myeloid_mean_logcpm = myeloid_profile, global_detection_fraction = global_profiles$detection,
    myeloid_detection_fraction = myeloid_profiles$detection), marker_gene_rankings = marker_rankings,
  selected_marker_genes = marker_selected, marker_set_summary = marker_set_summary, marker_scores_all_samples = score_all,
  marker_scores_primary = score_primary, marker_score_phase_tests = score_phase_tests, array_weights = array_weight_table,
  model_diagnostics = model_diagnostics, limma_models = model_results, gsea_results = gsea_results,
  gsea_summary = gsea_summary, gene_effect_comparison = gene_effect_comparison, deg_count_summary = deg_count_summary,
  core_gene_summary = core_gene_summary, signature_score_table = signature_score_table, signature_score_models = signature_score_models,
  signature_score_status = signature_score_status, module_adjusted_results = module_adjusted_results,
  module_analysis_status = module_analysis_status, settings = list(bulk_input = "Step46 expression_original; continuous transformed/normalized expression",
    bulk_transformations_applied = "None", raw_bulk_counts_used = FALSE, BayesPrism_used = FALSE,
    batch_term_used = FALSE, marker_selection = paste("Myeloid mean log2CPM >=", MIN_MYELOID_MEAN_LOGCPM,
      "; specificity >=", MIN_SPECIFICITY_LOG2, "; prevalence >=", MIN_STATE_PREVALENCE, "myeloid states; target signature excluded"),
    max_reference_cells_per_type = MAX_REF_CELLS_PER_TYPE, min_cells_per_source_state = MIN_CELLS_PER_SOURCE_STATE,
    marker_score = "Mean gene-wise z-score across all 39 Bulk samples", limma = "Common arrayWeights + lmFit + contrasts.fit + eBayes(trend=TRUE, robust=TRUE)",
    primary_contrast = "Established late (14+28 dpi) vs Early", GSEA_ranking = "limma moderated t statistic",
    significance = "BH adjusted P < 0.05", DEG_effect_threshold = "|logFC| >= 0.5"))

ANALYSIS_OBJECT_FILE <- file.path(RDS_DIR, paste0(PREFIX, "_analysis_object.rds"))
saveRDS(analysis_object, ANALYSIS_OBJECT_FILE, compress = "xz")

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, paste0(PREFIX, "_REPORT.txt")), useBytes = TRUE)

writeLines(capture.output(utils::sessionInfo()), file.path(OUT_DIR, paste0(PREFIX, "_sessionInfo.txt")),
  useBytes = TRUE)

expected_outputs <- c("tables/01_input_audit.csv", "tables/02_bulk_matrix_and_primary_sample_audit.csv",
  "tables/03_reference_cell_counts_by_type_source.csv", "tables/04_reference_cell_counts_retained.csv",
  "tables/05_myeloid_marker_gene_rankings.csv", "tables/06_selected_myeloid_marker_genes.csv", "tables/07_marker_set_summary.csv",
  "tables/08_bulk_myeloid_marker_scores_all_samples.csv", "tables/09_primary_group_marker_score_summary.csv",
  "tables/10_primary_marker_score_phase_tests.csv", "tables/11_limma_model_diagnostics.csv", "tables/12_limma_array_weights.csv",
  "tables/13_limma_primary_gene_level_results.csv", "tables/14_Clec7a_signature_GSEA_unadjusted_vs_adjusted.csv",
  "tables/15_Clec7a_signature_GSEA_leading_edge_genes.csv", "tables/16_gene_effect_concordance_vs_unadjusted.csv",
  "tables/17_DEG_count_summary_unadjusted_vs_adjusted.csv", "tables/18_core_gene_effects_unadjusted_vs_adjusted.csv",
  "tables/19_primary_sample_Clec7a_signature_score.csv", "tables/20_Clec7a_signature_score_phase_effects.csv",
  "tables/21_WGCNA_module_eigengene_phase_effects_adjusted_for_Myeloid9_score.csv", paste0("RDS/",
    PREFIX, "_analysis_object.rds"), paste0("figures/", PREFIX, "_myeloid_marker_scores_by_phase.pdf"),
  paste0("figures/", PREFIX, "_myeloid_marker_scores_by_phase.png"), paste0("figures/", PREFIX, "_GSEA_NES_comparison.pdf"),
  paste0("figures/", PREFIX, "_GSEA_NES_comparison.png"), paste0(PREFIX, "_REPORT.txt"), paste0(PREFIX,
    "_sessionInfo.txt"), paste0(PREFIX, ".log"))

output_index <- data.frame(file = expected_outputs, exists = file.exists(file.path(OUT_DIR, expected_outputs)),
  purpose = c("input paths and use/exclusion audit", "Bulk dimensions, expression scale and primary sample set",
    "Reference cell counts by cell class and source dataset before downsampling", "Reference cell counts retained after QC/downsampling",
    "Full gene-level myeloid specificity ranking for both marker sets", "Selected target-signature-independent myeloid marker genes",
    "Marker set size, overlap and prespecified selection thresholds", "Myeloid9/Myeloid10 marker scores for all 39 samples",
    "Primary-group descriptive marker-score statistics", "Primary phase tests for marker scores",
    "Limma design rank, residual df, VIF and condition diagnostics", "Common sample-level array weights used across models",
    "Gene-level limma results for all estimable models", "Frozen Clec7a-UP signature GSEA before/after marker-score adjustment",
    "GSEA leading-edge genes for each model", "Gene-level effect concordance versus phase-only model",
    "DEG counts under the FDR/effect thresholds", "Six core gene effects under each model", "Secondary sample-level Clec7a-UP signature score",
    "Phase effects on the Clec7a-UP signature score", "Optional Step49 module eigengene tests adjusted for Myeloid9 score",
    "Complete reproducible Step111 analysis object", "Diagnostic marker-score distribution by primary phase",
    "Diagnostic marker-score distribution by primary phase", "Diagnostic signature NES comparison",
    "Diagnostic signature NES comparison", "Analysis report and interpretation safeguards", "R session and package versions",
    "Execution log"), stringsAsFactors = FALSE)

utils::write.csv(output_index, file = file.path(OUT_DIR, "00_OUTPUT_INDEX.csv"), row.names = FALSE, fileEncoding = "UTF-8",
  na = "")
