source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

required_pkgs <- c("SeuratObject", "Matrix")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), "\nInstall them before running Step 85.")
}

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

input_rds <- file.path("D:/zz/results/20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

step84_dir <- "D:/zz/results/84_pySCENIC_input_preparation"

step84_state_counts <- file.path(step84_dir, "03_84_pySCENIC_myeloid_state_counts_v01.csv")

step84_final_qc <- file.path(step84_dir, "13_84_pySCENIC_FINAL_QC_v01.csv")

outdir <- "D:/zz/results/85_pySCENIC_GRNBoost2_input_preparation"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

version_tag <- "v01"

expected_n_cells <- 58533L
expected_n_states <- 9L

state_col <- "myeloid_state"

formal_target <- "Clec7a+ lipid-handling Microglia"
legacy_target <- "Clec7a+ lipid-associated Microglia"

state_cap <- 3000L
sampling_seed <- 8501L

gene_detection_fraction <- 0.01

f_output_index <- file.path(outdir, paste0("00_85_pySCENIC_OUTPUT_INDEX_", version_tag, ".csv"))

f_input_audit <- file.path(outdir, paste0("01_85_pySCENIC_input_audit_", version_tag, ".csv"))

f_sampling_plan <- file.path(outdir, paste0("02_85_pySCENIC_state_sampling_plan_", version_tag, ".csv"))

f_selected_cells <- file.path(outdir, paste0("03_85_pySCENIC_GRN_discovery_selected_cells_", version_tag,
  ".csv.gz"))

f_composition_audit <- file.path(outdir, paste0("04_85_pySCENIC_sampling_composition_audit_", version_tag,
  ".csv"))

f_gene_filter_summary <- file.path(outdir, paste0("05_85_pySCENIC_gene_filter_summary_", version_tag,
  ".csv"))

f_gene_filter_stats <- file.path(outdir, paste0("06_85_pySCENIC_gene_filter_stats_", version_tag, ".csv.gz"))

f_genes <- file.path(outdir, paste0("07_85_pySCENIC_GRN_discovery_genes_", version_tag, ".tsv.gz"))

f_cells <- file.path(outdir, paste0("08_85_pySCENIC_GRN_discovery_cells_", version_tag, ".tsv.gz"))

f_matrix <- file.path(outdir, paste0("09_85_pySCENIC_GRN_discovery_counts_", version_tag, ".mtx.gz"))

f_discovery_metadata <- file.path(outdir, paste0("10_85_pySCENIC_GRN_discovery_metadata_", version_tag,
  ".csv.gz"))

f_full_manifest <- file.path(outdir, paste0("11_85_pySCENIC_full_58533_cell_manifest_", version_tag,
  ".csv.gz"))

f_matrix_manifest <- file.path(outdir, paste0("12_85_pySCENIC_matrix_manifest_", version_tag, ".csv"))

f_final_qc <- file.path(outdir, paste0("13_85_pySCENIC_FINAL_QC_", version_tag, ".csv"))

f_report <- file.path(outdir, paste0("14_85_pySCENIC_REPORT_", version_tag, ".txt"))

f_session <- file.path(outdir, paste0("15_85_pySCENIC_sessionInfo_", version_tag, ".txt"))

make_qc <- qc_row

write_csv_gz <- function(x, file) {
    con <- gzfile(file, open = "wt")
    on.exit(close(con), add = TRUE)
    write.csv(x, file = con, row.names = FALSE, quote = TRUE, na = "")
}

write_lines_gz <- function(x, file) {
    con <- gzfile(file, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(text = x, con = con, useBytes = TRUE)
}

gzip_file_base <- function(input_file, output_file, remove_input = TRUE, chunk_size = 8L * 1024L * 1024L) {
    in_con <- file(input_file, open = "rb")
    out_con <- gzfile(output_file, open = "wb")
    on.exit({
        try(close(in_con), silent = TRUE)
        try(close(out_con), silent = TRUE)
    }, add = TRUE)
    repeat {
        buf <- readBin(in_con, what = "raw", n = chunk_size)
        if (length(buf) == 0) {
            break
        }
        writeBin(buf, out_con)
    }
    close(in_con)
    close(out_con)
    if (!file.exists(output_file)) {
        stop("Failed to create gzip output:\n", output_file)
    }
    if (isTRUE(remove_input)) {
        unlink(input_file)
    }
    invisible(output_file)
}

allocate_proportional <- function(group_sizes, total_take, group_names = names(group_sizes)) {

  group_sizes <- as.integer(group_sizes)

  if (is.null(group_names)) {
    group_names <- as.character(seq_along(group_sizes))
  }

  names(group_sizes) <- group_names

  if (total_take >= sum(group_sizes)) {
    return(group_sizes)
  }

  desired <- group_sizes/sum(group_sizes) * total_take

  allocation <- floor(desired)

  allocation <- pmin(allocation, group_sizes)

  remainder <- total_take - sum(allocation)

  if (remainder > 0) {

    fractional <- desired - floor(desired)

    ord <- order(-fractional, group_names)

    for (i in ord) {

      if (remainder <= 0) {
        break
      }

      if (allocation[i] < group_sizes[i]) {
        allocation[i] <- allocation[i] + 1L
        remainder <- remainder - 1L
      }
    }
  }

  if (remainder > 0) {

    while (remainder > 0) {

      eligible <- which(allocation < group_sizes)

      if (length(eligible) == 0) {
        stop("Internal allocation error: no eligible strata remain ", "but sampling remainder is still positive.")
      }

      for (i in eligible) {

        if (remainder <= 0) {
          break
        }

        allocation[i] <- allocation[i] + 1L
        remainder <- remainder - 1L
      }
    }
  }

  if (sum(allocation) != total_take) {
    stop("Internal allocation error: allocated ", sum(allocation), " cells but expected ", total_take,
      ".")
  }

  allocation
}

sample_one_state <- function(meta_state, cap, strata_cols) {
    n_state <- nrow(meta_state)
    if (n_state <= cap) {
        return(rownames(meta_state))
    }
    if (length(strata_cols) == 0) {
        return(sample(rownames(meta_state), size = cap, replace = FALSE))
    }
    strata_df <- meta_state[, strata_cols, drop = FALSE]
    strata_df[] <- lapply(strata_df, function(x) {
        x <- as.character(x)
        x[is.na(x) | x == ""] <- "<NA>"
        x
    })
    stratum_key <- do.call(paste, c(strata_df, sep = " || "))
    names(stratum_key) <- rownames(meta_state)
    tab <- table(stratum_key)
    allocations <- allocate_proportional(group_sizes = as.integer(tab), total_take = cap, group_names = names(tab))
    names(allocations) <- names(tab)
    selected <- character(0)
    for (stratum in names(allocations)) {
        n_take <- allocations[[stratum]]
        cells_here <- names(stratum_key[stratum_key == stratum])
        if (n_take > 0) {
            selected <- c(selected, sample(cells_here, size = n_take, replace = FALSE))
        }
    }
    if (length(selected) != cap) {
        stop("Sampling error for state '", unique(meta_state$myeloid_state_frozen), "': selected ", length(selected), 
            " cells instead of ", cap, ".")
    }
    selected
}

sample_sparse_values <- function(mat, max_values = 100000L, seed = 8502L) {
    if (!inherits(mat, "Matrix")) {
        stop("Expected a sparse Matrix object.")
    }
    if (!("x" %in% slotNames(mat))) {
        stop("Sparse Matrix object does not contain an x slot.")
    }
    vals <- mat@x
    if (length(vals) == 0) {
        return(numeric(0))
    }
    if (length(vals) > max_values) {
        set.seed(seed)
        vals <- sample(vals, size = max_values, replace = FALSE)
    }
    as.numeric(vals)
}

required_files <- c(input_rds, step84_state_counts, step84_final_qc)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    stop("Required input file(s) are missing:\n", paste(missing_files, collapse = "\n"))
}

step84_qc <- read.csv(step84_final_qc, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("check", "status", "detail") %in% colnames(step84_qc))) {
  stop("Unexpected schema in Step 84 FINAL_QC file.")
}

step84_fail <- step84_qc[step84_qc$status == "FAIL", , drop = FALSE]

if (nrow(step84_fail) > 0) {
    stop("Step 84 contains blocking FAIL item(s). ", "Step 85 will not proceed.")
}

cat("Step 84 QC gate: PASS\n", "  blocking FAIL =", nrow(step84_fail), "\n\n")

obj <- readRDS(input_rds)

if (!inherits(obj, "Seurat")) {
    stop("Input RDS is not a Seurat object. Class: ", paste(class(obj), collapse = ", "))
}

if (ncol(obj) != expected_n_cells) {
    stop("cell count mismatch: observed ", ncol(obj), ", expected ", expected_n_cells, ".")
}

if (!state_col %in% colnames(obj[[]])) {
    stop("Required metadata column not found: ", state_col)
}

if (!"RNA" %in% Assays(obj)) {
    stop("RNA assay is missing.")
}

meta <- obj[[]]

if (!identical(rownames(meta), colnames(obj))) {
    meta <- meta[colnames(obj), , drop = FALSE]
}

meta$cell <- rownames(meta)

meta$myeloid_state_original <- as.character(meta[[state_col]])

meta$myeloid_state_frozen <- meta$myeloid_state_original

meta$myeloid_state_frozen[meta$myeloid_state_frozen == legacy_target] <- formal_target

observed_state_counts_original <- as.data.frame(table(meta$myeloid_state_original, useNA = "ifany"),
  stringsAsFactors = FALSE)

colnames(observed_state_counts_original) <- c("myeloid_state_original", "n_cells")

observed_state_counts_original$n_cells <- as.integer(observed_state_counts_original$n_cells)

step84_counts <- read.csv(step84_state_counts, check.names = FALSE, stringsAsFactors = FALSE)

step84_counts$n_cells <- as.integer(step84_counts$n_cells)

cmp_a <- observed_state_counts_original[order(observed_state_counts_original$myeloid_state_original),
  , drop = FALSE]

cmp_b <- step84_counts[order(step84_counts$myeloid_state_original), c("myeloid_state_original", "n_cells"),
  drop = FALSE]

rownames(cmp_a) <- NULL
rownames(cmp_b) <- NULL

if (!identical(cmp_a, cmp_b)) {
    stop("Current object's state counts do not match ", "Step 84 state counts.")
}

n_states_frozen <- length(unique(meta$myeloid_state_frozen))

if (n_states_frozen != expected_n_states) {
    stop("After formal name mapping, observed ", n_states_frozen, " myeloid states; expected ", expected_n_states, 
        ".")
}

if (!formal_target %in% meta$myeloid_state_frozen) {
    stop("Formal target state is absent after name mapping:\n", formal_target)
}

strata_cols <- character(0)

if ("dataset" %in% colnames(meta)) {
  strata_cols <- c(strata_cols, "dataset")
}

if ("time" %in% colnames(meta)) {
  strata_cols <- c(strata_cols, "time")
}

if (length(strata_cols) == 0) {
} else {
}

state_levels <- unique(meta$myeloid_state_frozen)

state_counts <- data.frame(myeloid_state = state_levels, n_total = vapply(state_levels, function(z) {
  sum(meta$myeloid_state_frozen == z)
}, integer(1)), stringsAsFactors = FALSE)

state_counts <- state_counts[order(state_counts$n_total, decreasing = TRUE), , drop = FALSE]

state_counts$n_target <- pmin(state_counts$n_total, state_cap)

state_counts$strategy <- ifelse(state_counts$n_total <= state_cap, "keep_all", paste0("proportional_stratified_cap_",
  state_cap))

expected_discovery_n <- sum(state_counts$n_target)

print(state_counts)

set.seed(sampling_seed)

selected_set <- character(0)

for (state_i in state_counts$myeloid_state) {
    idx <- which(meta$myeloid_state_frozen == state_i)
    meta_state <- meta[idx, , drop = FALSE]
    selected_i <- sample_one_state(meta_state = meta_state, cap = state_cap, strata_cols = strata_cols)
    selected_set <- c(selected_set, selected_i)
}

if (anyDuplicated(selected_set) > 0) {
    stop("Duplicate cell IDs were generated during sampling.")
}

if (length(selected_set) != expected_discovery_n) {
    stop("Selected discovery cell number = ", length(selected_set), "; expected = ", expected_discovery_n, 
        ".")
}

selected_cells <- colnames(obj)[colnames(obj) %in% selected_set]

if (length(selected_cells) != expected_discovery_n) {
    stop("Selected cell order reconstruction failed.")
}

selected_flag <- meta$cell %in% selected_cells

meta$is_GRN_discovery <- selected_flag

meta$GRN_discovery_order <- NA_integer_

meta$GRN_discovery_order[match(selected_cells, meta$cell)] <- seq_along(selected_cells)

state_counts$n_selected <- vapply(state_counts$myeloid_state, function(z) {
  sum(meta$myeloid_state_frozen == z & meta$is_GRN_discovery)
}, integer(1))

state_counts$retention_fraction <- state_counts$n_selected/state_counts$n_total

state_counts$state_cap <- state_cap
state_counts$sampling_seed <- sampling_seed

state_counts <- state_counts[, c("myeloid_state", "n_total", "n_target", "n_selected", "retention_fraction",
  "strategy", "state_cap", "sampling_seed")]

write.csv(state_counts, f_sampling_plan, row.names = FALSE)

audit_cols <- c("myeloid_state_frozen", strata_cols)

composition_all <- aggregate(rep(1L, nrow(meta)), by = meta[, audit_cols, drop = FALSE], FUN = sum)

colnames(composition_all)[ncol(composition_all)] <- "n_all"

meta_selected <- meta[selected_cells, , drop = FALSE]

composition_selected <- aggregate(rep(1L, nrow(meta_selected)), by = meta_selected[, audit_cols, drop = FALSE],
  FUN = sum)

colnames(composition_selected)[ncol(composition_selected)] <- "n_selected"

composition_audit <- merge(composition_all, composition_selected, by = audit_cols, all.x = TRUE, sort = FALSE)

composition_audit$n_selected[is.na(composition_audit$n_selected)] <- 0L

state_total_lookup <- setNames(state_counts$n_total, state_counts$myeloid_state)

state_selected_lookup <- setNames(state_counts$n_selected, state_counts$myeloid_state)

composition_audit$prop_all_within_state <- composition_audit$n_all/state_total_lookup[composition_audit$myeloid_state_frozen]

composition_audit$prop_selected_within_state <- composition_audit$n_selected/state_selected_lookup[composition_audit$myeloid_state_frozen]

composition_audit$delta_selected_minus_all <- composition_audit$prop_selected_within_state - composition_audit$prop_all_within_state

composition_audit$abs_delta <- abs(composition_audit$delta_selected_minus_all)

write.csv(composition_audit, f_composition_audit, row.names = FALSE)

selected_export_cols <- unique(c("cell", "myeloid_state_original", "myeloid_state_frozen", strata_cols,
  "is_GRN_discovery", "GRN_discovery_order"))

selected_manifest <- meta[selected_cells, selected_export_cols, drop = FALSE]

selected_manifest <- selected_manifest[order(selected_manifest$GRN_discovery_order), , drop = FALSE]

write_csv_gz(selected_manifest, f_selected_cells)

full_manifest_cols <- unique(c("cell", "myeloid_state_original", "myeloid_state_frozen", strata_cols,
  "is_GRN_discovery", "GRN_discovery_order"))

full_manifest <- meta[, full_manifest_cols, drop = FALSE]

write_csv_gz(full_manifest, f_full_manifest)

RNA_assay <- obj[["RNA"]]

rna_layers <- tryCatch(Layers(RNA_assay), error = function(e) character(0))

count_layers <- grep("^counts($|\\.)", rna_layers, value = TRUE)

if (length(count_layers) != 1) {
    stop("Step 85 expected exactly one RNA count layer, ", "but detected ", length(count_layers), ": ", 
        paste(count_layers, collapse = "; "))
}

if (count_layers != "counts") {
    warning("The single RNA count layer is named '", count_layers, "', not exactly 'counts'. ", "Proceeding because Step 84 confirmed a single count layer.")
}

counts_all <- LayerData(RNA_assay, layer = count_layers)

if (!inherits(counts_all, "Matrix")) {
    counts_all <- as(counts_all, "dgCMatrix")
}

if (!identical(colnames(counts_all), colnames(obj))) {
    if (!setequal(colnames(counts_all), colnames(obj))) {
        stop("RNA count cells do not match Seurat object cells.")
    }
    counts_all <- counts_all[, colnames(obj), drop = FALSE]
}

counts_discovery <- counts_all[, selected_cells, drop = FALSE]

if (!identical(colnames(counts_discovery), selected_cells)) {
    stop("Discovery count matrix column order does not match ", "the selected-cell manifest.")
}

if (ncol(counts_discovery) != expected_discovery_n) {
    stop("Discovery count matrix has ", ncol(counts_discovery), " cells; expected ", expected_discovery_n, 
        ".")
}

rm(counts_all, RNA_assay, obj)

invisible(gc(FALSE))

n_discovery_cells <- ncol(counts_discovery)

min_detected_cells <- ceiling(gene_detection_fraction * n_discovery_cells)

min_counts_per_gene <- ceiling(3 * gene_detection_fraction * n_discovery_cells)

detected_cells <- Matrix::rowSums(counts_discovery != 0)

total_counts <- Matrix::rowSums(counts_discovery)

gene_filter_stats <- data.frame(gene = rownames(counts_discovery), total_counts = as.numeric(total_counts),
  detected_cells = as.integer(detected_cells), detected_fraction = as.numeric(detected_cells/n_discovery_cells),
  pass_min_detected_cells = detected_cells >= min_detected_cells, pass_min_total_counts = total_counts >=
    min_counts_per_gene, stringsAsFactors = FALSE)

gene_filter_stats$keep_for_GRNBoost2 <- gene_filter_stats$pass_min_detected_cells & gene_filter_stats$pass_min_total_counts

n_genes_before <- nrow(gene_filter_stats)

n_genes_after <- sum(gene_filter_stats$keep_for_GRNBoost2)

if (n_genes_after == 0) {
    stop("Gene filtering retained zero genes. ", "Do not proceed.")
}

write_csv_gz(gene_filter_stats, f_gene_filter_stats)

gene_filter_summary <- data.frame(metric = c("discovery_cells", "gene_detection_fraction", "min_detected_cells",
  "min_counts_per_gene", "genes_before_filter", "genes_after_filter", "genes_removed", "gene_retention_fraction"),
  value = c(n_discovery_cells, gene_detection_fraction, min_detected_cells, min_counts_per_gene, n_genes_before,
    n_genes_after, n_genes_before - n_genes_after, n_genes_after/n_genes_before), stringsAsFactors = FALSE)

write.csv(gene_filter_summary, f_gene_filter_summary, row.names = FALSE)

keep_genes <- gene_filter_stats$gene[gene_filter_stats$keep_for_GRNBoost2]

counts_filtered <- counts_discovery[keep_genes, , drop = FALSE]

if (!identical(rownames(counts_filtered), keep_genes)) {
    stop("Filtered matrix gene order mismatch.")
}

if (!identical(colnames(counts_filtered), selected_cells)) {
    stop("Filtered matrix cell order mismatch.")
}

sampled_values <- sample_sparse_values(counts_filtered, max_values = 100000L, seed = 8502L)

counts_nonnegative <- if (length(sampled_values) > 0) {
  all(is.finite(sampled_values) & sampled_values >= 0)
} else {
  TRUE
}

counts_integer_like <- if (length(sampled_values) > 0) {
  all(abs(sampled_values - round(sampled_values)) < 1e-08)
} else {
  TRUE
}

write_lines_gz(rownames(counts_filtered), f_genes)

write_lines_gz(colnames(counts_filtered), f_cells)

discovery_metadata <- selected_manifest[match(colnames(counts_filtered), selected_manifest$cell), , drop = FALSE]

if (!identical(discovery_metadata$cell, colnames(counts_filtered))) {
    stop("Discovery metadata order does not match matrix columns.")
}

write_csv_gz(discovery_metadata, f_discovery_metadata)

tmp_mtx <- file.path(outdir, paste0("09_85_pySCENIC_GRN_discovery_counts_", version_tag, ".mtx"))

if (file.exists(tmp_mtx)) {
  unlink(tmp_mtx)
}

if (file.exists(f_matrix)) {
  unlink(f_matrix)
}

Matrix::writeMM(counts_filtered, file = tmp_mtx)

if (!file.exists(tmp_mtx)) {
  stop("Matrix::writeMM did not create the expected temporary file.")
}

gzip_file_base(input_file = tmp_mtx, output_file = f_matrix, remove_input = TRUE)

nnzero_filtered <- Matrix::nnzero(counts_filtered)

matrix_denominator <- as.numeric(nrow(counts_filtered)) * as.numeric(ncol(counts_filtered))

matrix_sparsity <- 1 - nnzero_filtered/matrix_denominator

matrix_manifest <- data.frame(item = c("matrix_file", "matrix_format", "matrix_orientation", "row_dimension",
  "column_dimension", "n_genes", "n_cells", "nnzero", "sparsity", "value_type", "sampled_values_nonnegative",
  "sampled_values_integer_like", "gene_file", "cell_file", "metadata_file", "full_atlas_manifest",
  "state_cap", "sampling_seed", "sampling_strata", "gene_detection_fraction", "min_detected_cells",
  "min_counts_per_gene", "formal_target_name"), value = c(basename(f_matrix), "Matrix Market coordinate sparse matrix, gzip-compressed",
  "genes x cells", "genes", "cells", nrow(counts_filtered), ncol(counts_filtered), nnzero_filtered,
  matrix_sparsity, "raw RNA UMI counts", counts_nonnegative, counts_integer_like, basename(f_genes),
  basename(f_cells), basename(f_discovery_metadata), basename(f_full_manifest), state_cap, sampling_seed,
  ifelse(length(strata_cols) == 0, "state only", paste(c("myeloid_state", strata_cols), collapse = " x ")),
  gene_detection_fraction, min_detected_cells, min_counts_per_gene, formal_target), stringsAsFactors = FALSE)

write.csv(matrix_manifest, f_matrix_manifest, row.names = FALSE)

input_audit <- data.frame(item = c("frozen_RDS", "Step84_state_counts", "Step84_FINAL_QC", "expected_total_cells",
  "observed_total_cells", "expected_states", "observed_states_after_formal_mapping", "formal_target",
  "formal_target_total_cells", "RNA_count_layer", "state_cap", "expected_discovery_cells", "observed_discovery_cells",
  "sampling_seed", "sampling_strata", "gene_filter_rule"), value = c(normalizePath(input_rds, winslash = "/",
  mustWork = TRUE), normalizePath(step84_state_counts, winslash = "/", mustWork = TRUE), normalizePath(step84_final_qc,
  winslash = "/", mustWork = TRUE), expected_n_cells, nrow(meta), expected_n_states, n_states_frozen,
  formal_target, sum(meta$myeloid_state_frozen == formal_target), count_layers, state_cap, expected_discovery_n,
  n_discovery_cells, sampling_seed, ifelse(length(strata_cols) == 0, "none", paste(strata_cols, collapse = " x ")),
  paste0("detected in >= ", gene_detection_fraction * 100, "% discovery cells AND total raw counts >= 3 * ",
    gene_detection_fraction * 100, "% of discovery-cell number")), stringsAsFactors = FALSE)

write.csv(input_audit, f_input_audit, row.names = FALSE)

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("Step 84 blocking QC", if (nrow(step84_fail) == 0) "PASS" else "FAIL",
  paste0("Step 84 blocking FAIL items = ", nrow(step84_fail)))

qc_list[[length(qc_list) + 1]] <- make_qc("total cell number", if (nrow(meta) == expected_n_cells) "PASS" else "FAIL",
  paste0("Observed = ", nrow(meta), "; expected = ", expected_n_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("state number after formal name mapping", if (n_states_frozen ==
  expected_n_states) "PASS" else "FAIL", paste0("Observed = ", n_states_frozen, "; expected = ", expected_n_states))

qc_list[[length(qc_list) + 1]] <- make_qc("Formal Clec7a target present", if (formal_target %in% meta$myeloid_state_frozen) "PASS" else "FAIL",
  paste0("Target = ", formal_target, "; n = ", sum(meta$myeloid_state_frozen == formal_target)))

qc_list[[length(qc_list) + 1]] <- make_qc("Discovery cell number", if (n_discovery_cells == expected_discovery_n) "PASS" else "FAIL",
  paste0("Observed = ", n_discovery_cells, "; expected from cap plan = ", expected_discovery_n))

qc_list[[length(qc_list) + 1]] <- make_qc("No duplicate discovery cells", if (anyDuplicated(selected_cells) ==
  0) "PASS" else "FAIL", paste0("Duplicated selected cell IDs = ", anyDuplicated(selected_cells)))

qc_list[[length(qc_list) + 1]] <- make_qc("All 9 states retained in discovery matrix", if (length(unique(meta_selected$myeloid_state_frozen)) ==
  expected_n_states) "PASS" else "FAIL", paste0("Discovery states = ", length(unique(meta_selected$myeloid_state_frozen))))

sampling_plan_valid <- all(state_counts$n_selected == state_counts$n_target)

qc_list[[length(qc_list) + 1]] <- make_qc("Per-state sampling targets achieved", if (sampling_plan_valid) "PASS" else "FAIL",
  paste0("All observed selected counts equal cap-based targets = ", sampling_plan_valid))

qc_list[[length(qc_list) + 1]] <- make_qc("Filtered matrix cell order", if (identical(colnames(counts_filtered),
  selected_cells)) "PASS" else "FAIL", "Matrix columns must exactly match exported discovery-cell order.")

qc_list[[length(qc_list) + 1]] <- make_qc("Filtered genes retained", if (n_genes_after > 0) "PASS" else "FAIL",
  paste0("Genes before = ", n_genes_before, "; after = ", n_genes_after))

qc_list[[length(qc_list) + 1]] <- make_qc("Filtered counts non-negative", if (isTRUE(counts_nonnegative)) "PASS" else "FAIL",
  paste0("Sampled-value result = ", counts_nonnegative))

qc_list[[length(qc_list) + 1]] <- make_qc("Filtered counts integer-like", if (isTRUE(counts_integer_like)) "PASS" else "WARN",
  paste0("Sampled-value result = ", counts_integer_like))

max_comp_delta <- if (nrow(composition_audit) > 0) {
  max(composition_audit$abs_delta, na.rm = TRUE)
} else {
  NA_real_
}

qc_list[[length(qc_list) + 1]] <- make_qc("Sampling composition preservation audit", "PASS", paste0("Maximum absolute within-state stratum proportion difference = ",
  signif(max_comp_delta, 5), ". This is descriptive and not used to force equal time-point sizes."))

required_export_files <- c(f_genes, f_cells, f_matrix, f_discovery_metadata, f_full_manifest, f_matrix_manifest)

export_files_ok <- all(file.exists(required_export_files))

qc_list[[length(qc_list) + 1]] <- make_qc("Required GRNBoost2 input files exported", if (export_files_ok) "PASS" else "FAIL",
  paste0("All required export files exist = ", export_files_ok))

final_qc <- do.call(rbind, qc_list)

rownames(final_qc) <- NULL

write.csv(final_qc, f_final_qc, row.names = FALSE)

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

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    f_report)

writeLines(capture.output(sessionInfo()), f_session)

output_files <- c(f_input_audit, f_sampling_plan, f_selected_cells, f_composition_audit, f_gene_filter_summary,
  f_gene_filter_stats, f_genes, f_cells, f_matrix, f_discovery_metadata, f_full_manifest, f_matrix_manifest,
  f_final_qc, f_report, f_session)

output_index <- data.frame(file = basename(output_files), exists = file.exists(output_files), size_MB = ifelse(file.exists(output_files),
  file.info(output_files)$size/1024^2, NA_real_), full_path = normalizePath(output_files, winslash = "/",
  mustWork = FALSE), stringsAsFactors = FALSE)

write.csv(output_index, f_output_index, row.names = FALSE)

print(state_counts[, c("myeloid_state", "n_total", "n_selected", "retention_fraction")])

cat("  WARN =", n_warn, "\n")
cat("  FAIL =", n_fail, "\n")

if (n_fail > 0) {
    stop("\nStep 85 generated QC outputs but contains ", n_fail, " blocking FAIL item(s).\n", "Do NOT run GRNBoost2 yet.\n", 
        "Review:\n", f_final_qc)
}
