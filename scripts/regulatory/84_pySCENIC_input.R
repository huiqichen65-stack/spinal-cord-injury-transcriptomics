source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 160)

required_pkgs <- c("SeuratObject", "Matrix")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), "\nInstall them before running Step 84.")
}

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

input_rds <- file.path(project_path("results/20_AUCell_lipid_metabolism_scoring"), "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

outdir <- project_path("results/84_pySCENIC_input_preparation")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

version_tag <- "v01"

expected_n_cells <- 58533L
expected_n_states <- 9L

formal_target <- "Clec7a+ lipid-handling Microglia"

legacy_target <- "Clec7a+ lipid-associated Microglia"

f_object_summary <- file.path(outdir, paste0("01_84_pySCENIC_object_summary_", version_tag, ".csv"))

f_metadata_fields <- file.path(outdir, paste0("02_84_pySCENIC_metadata_fields_", version_tag, ".csv"))

f_state_counts <- file.path(outdir, paste0("03_84_pySCENIC_myeloid_state_counts_", version_tag, ".csv"))

f_state_name_audit <- file.path(outdir, paste0("04_84_pySCENIC_state_name_audit_", version_tag, ".csv"))

f_assay_audit <- file.path(outdir, paste0("05_84_pySCENIC_assay_layer_audit_", version_tag, ".csv"))

f_counts_qc <- file.path(outdir, paste0("06_84_pySCENIC_RNA_counts_QC_", version_tag, ".csv"))

f_gene_summary <- file.path(outdir, paste0("07_84_pySCENIC_gene_identifier_summary_", version_tag, ".csv"))

f_duplicate_genes <- file.path(outdir, paste0("08_84_pySCENIC_duplicate_gene_names_", version_tag, ".csv"))

f_metadata_export <- file.path(outdir, paste0("09_84_pySCENIC_cell_metadata_", version_tag, ".csv.gz"))

f_umap_export <- file.path(outdir, paste0("10_84_pySCENIC_frozen_UMAP_coordinates_", version_tag, ".csv.gz"))

f_reduction_audit <- file.path(outdir, paste0("11_84_pySCENIC_reduction_audit_", version_tag, ".csv"))

f_manifest <- file.path(outdir, paste0("12_84_pySCENIC_input_manifest_", version_tag, ".csv"))

f_final_qc <- file.path(outdir, paste0("13_84_pySCENIC_FINAL_QC_", version_tag, ".csv"))

f_report <- file.path(outdir, paste0("14_84_pySCENIC_REPORT_", version_tag, ".txt"))

f_session <- file.path(outdir, paste0("15_84_pySCENIC_sessionInfo_", version_tag, ".txt"))

f_output_index <- file.path(outdir, paste0("00_84_pySCENIC_OUTPUT_INDEX_", version_tag, ".csv"))

write_csv_gz <- function(x, file) {
    con <- gzfile(file, open = "wt")
    on.exit(close(con), add = TRUE)
    write.csv(x, file = con, row.names = FALSE, quote = TRUE, na = "")
}

safe_character <- function(x) {
  if (is.factor(x)) {
    as.character(x)
  } else {
    x
  }
}

make_qc <- qc_row

sample_matrix_values <- function(mat, max_values = 100000L) {
    if (inherits(mat, "Matrix")) {
        if ("x" %in% slotNames(mat)) {
            vals <- mat@x
            if (length(vals) == 0) {
                return(numeric(0))
            }
            if (length(vals) > max_values) {
                set.seed(8401)
                vals <- sample(vals, max_values)
            }
            return(as.numeric(vals))
        }
    }
    total_n <- length(mat)
    if (total_n == 0) {
        return(numeric(0))
    }
    set.seed(8401)
    idx <- if (total_n > max_values) {
        sample.int(total_n, max_values)
    }
    else {
        seq_len(total_n)
    }
    as.numeric(mat[idx])
}

if (!file.exists(input_rds)) {
  stop("input RDS does not exist:\n", input_rds)
}

obj <- readRDS(input_rds)

if (!inherits(obj, "Seurat")) {
  stop("Input object is not recognized as a Seurat object. Class: ", paste(class(obj), collapse = ", "))
}

n_cells <- ncol(obj)
n_features_default <- nrow(obj)

assay_names <- Assays(obj)
default_assay <- DefaultAssay(obj)

meta <- obj[[]]

if (nrow(meta) != n_cells) {
  stop("Metadata row number does not equal object cell number.")
}

if (is.null(rownames(meta))) {
  stop("Metadata does not contain cell IDs as row names.")
}

if (!identical(rownames(meta), colnames(obj))) {
    warning("Metadata row order differs from object cell order. ", "Export will explicitly align metadata to colnames(obj).")
    meta <- meta[colnames(obj), , drop = FALSE]
}

object_summary <- data.frame(item = c("object_class", "n_cells", "n_features_default_assay", "default_assay",
  "n_assays", "assays", "SeuratObject_version", "R_version"), value = c(paste(class(obj), collapse = "; "),
  as.character(n_cells), as.character(n_features_default), default_assay, as.character(length(assay_names)),
  paste(assay_names, collapse = "; "), as.character(packageVersion("SeuratObject")), R.version.string),
  stringsAsFactors = FALSE)

write.csv(object_summary, f_object_summary, row.names = FALSE)

metadata_fields <- data.frame(column = colnames(meta), class = vapply(meta, function(x) paste(class(x),
  collapse = ";"), character(1)), n_unique = vapply(meta, function(x) length(unique(x)), integer(1)),
  n_missing = vapply(meta, function(x) sum(is.na(x)), integer(1)), stringsAsFactors = FALSE)

write.csv(metadata_fields, f_metadata_fields, row.names = FALSE)

state_col <- NULL

if ("myeloid_state" %in% colnames(meta)) {
  state_col <- "myeloid_state"
}

if (is.null(state_col)) {
    state_candidates <- grep("myeloid.*state|state.*myeloid", colnames(meta), value = TRUE, ignore.case = TRUE)
    print(state_candidates)
    stop("Exact column 'myeloid_state' was not found.\n", "No state column will be guessed automatically.")
}

meta[[state_col]] <- safe_character(meta[[state_col]])

state_counts <- as.data.frame(table(meta[[state_col]], useNA = "ifany"), stringsAsFactors = FALSE)

colnames(state_counts) <- c("myeloid_state_original", "n_cells")

state_counts <- state_counts[order(state_counts$n_cells, decreasing = TRUE), , drop = FALSE]

write.csv(state_counts, f_state_counts, row.names = FALSE)

n_states <- length(unique(meta[[state_col]]))

print(state_counts)

original_states <- unique(meta[[state_col]])

state_name_audit <- data.frame(original_state = original_states, frozen_export_state = original_states,
  n_cells = vapply(original_states, function(z) sum(meta[[state_col]] == z, na.rm = TRUE), integer(1)),
  action = "unchanged", stringsAsFactors = FALSE)

legacy_idx <- state_name_audit$original_state == legacy_target

if (any(legacy_idx)) {
  state_name_audit$frozen_export_state[legacy_idx] <- formal_target
  state_name_audit$action[legacy_idx] <- "legacy label mapped to formal label in exported metadata only"
}

formal_idx <- state_name_audit$original_state == formal_target

if (any(formal_idx)) {
  state_name_audit$action[formal_idx] <- "already uses formal label"
}

write.csv(state_name_audit, f_state_name_audit, row.names = FALSE)

meta_export <- meta

meta_export$myeloid_state_frozen <- meta_export[[state_col]]

meta_export$myeloid_state_frozen[meta_export$myeloid_state_frozen == legacy_target] <- formal_target

target_present_after_mapping <- formal_target %in% unique(meta_export$myeloid_state_frozen)

target_n <- sum(meta_export$myeloid_state_frozen == formal_target, na.rm = TRUE)

time_col <- if ("time" %in% colnames(meta)) {
  "time"
} else {
  NA_character_
}

dataset_candidates <- grep("dataset|source|study|GSE", colnames(meta), value = TRUE, ignore.case = TRUE)

sample_candidates <- grep("sample|orig.ident|replicate|animal|mouse", colnames(meta), value = TRUE, ignore.case = TRUE)

print(time_col)

print(dataset_candidates)

print(sample_candidates)

assay_audit_list <- list()

for (assay_nm in assay_names) {
    assay_obj <- obj[[assay_nm]]
    assay_class <- paste(class(assay_obj), collapse = ";")
    layer_names <- tryCatch(Layers(assay_obj), error = function(e) character(0))
    if (length(layer_names) > 0) {
        tmp <- data.frame(assay = assay_nm, assay_class = assay_class, layer = layer_names, n_features = NA_integer_, 
            n_cells = NA_integer_, nnzero = NA_real_, sparsity = NA_real_, stringsAsFactors = FALSE)
        for (j in seq_along(layer_names)) {
            layer_nm <- layer_names[j]
            layer_mat <- tryCatch(LayerData(object = assay_obj, layer = layer_nm), error = function(e) NULL)
            if (!is.null(layer_mat)) {
                tmp$n_features[j] <- nrow(layer_mat)
                tmp$n_cells[j] <- ncol(layer_mat)
                if (inherits(layer_mat, "Matrix")) {
                  nz <- Matrix::nnzero(layer_mat)
                  tmp$nnzero[j] <- nz
                  denominator <- as.numeric(nrow(layer_mat)) * as.numeric(ncol(layer_mat))
                  tmp$sparsity[j] <- if (denominator > 0) {
                    1 - nz/denominator
                  }
                  else {
                    NA_real_
                  }
                }
            }
            rm(layer_mat)
            invisible(gc(FALSE))
        }
    }
    else {
        tmp <- data.frame(assay = assay_nm, assay_class = assay_class, layer = "legacy_assay_slots", 
            n_features = nrow(assay_obj), n_cells = ncol(assay_obj), nnzero = NA_real_, sparsity = NA_real_, 
            stringsAsFactors = FALSE)
    }
    assay_audit_list[[assay_nm]] <- tmp
}

assay_audit <- do.call(rbind, assay_audit_list)

rownames(assay_audit) <- NULL

write.csv(assay_audit, f_assay_audit, row.names = FALSE)

rna_exists <- "RNA" %in% assay_names

counts_qc <- data.frame(metric = character(0), value = character(0), stringsAsFactors = FALSE)

rna_counts_available <- FALSE
count_layer_type <- NA_character_
count_layer_names <- character(0)

counts_cell_coverage <- NA_integer_
counts_feature_count <- NA_integer_
counts_duplicate_cell_membership <- NA_integer_

values_nonnegative <- NA
values_integer_like <- NA

if (rna_exists) {
    RNA_assay <- obj[["RNA"]]
    RNA_layers <- tryCatch(Layers(RNA_assay), error = function(e) character(0))
    if (length(RNA_layers) > 0) {
        count_layer_names <- grep("^counts($|\\.)", RNA_layers, value = TRUE)
        if (length(count_layer_names) == 1) {
            count_layer_type <- "single_count_layer"
            counts_mat <- LayerData(RNA_assay, layer = count_layer_names)
            rna_counts_available <- TRUE
            counts_feature_count <- nrow(counts_mat)
            counts_cell_coverage <- ncol(counts_mat)
            counts_duplicate_cell_membership <- 0L
            vals <- sample_matrix_values(counts_mat)
            if (length(vals) > 0) {
                values_nonnegative <- all(is.finite(vals) & vals >= 0)
                values_integer_like <- all(abs(vals - round(vals)) < 1e-08)
            }
            single_cell_set_match <- setequal(colnames(counts_mat), colnames(obj))
            single_cell_order_match <- identical(colnames(counts_mat), colnames(obj))
            counts_qc <- rbind(counts_qc, data.frame(metric = c("RNA_count_layer_structure", "RNA_count_layer_names", 
                "RNA_counts_n_features", "RNA_counts_n_cells", "RNA_counts_cell_set_matches_object", 
                "RNA_counts_cell_order_matches_object", "sampled_counts_nonnegative", "sampled_counts_integer_like"), 
                value = c(count_layer_type, paste(count_layer_names, collapse = "; "), counts_feature_count, 
                  counts_cell_coverage, single_cell_set_match, single_cell_order_match, values_nonnegative, 
                  values_integer_like), stringsAsFactors = FALSE))
            rm(counts_mat, vals)
            invisible(gc(FALSE))
        }
        else if (length(count_layer_names) > 1) {
            count_layer_type <- "multiple_count_layers"
            all_count_cells <- character(0)
            all_count_genes <- character(0)
            sampled_nonnegative <- logical(0)
            sampled_integer <- logical(0)
            for (layer_nm in count_layer_names) {
                mat_i <- LayerData(RNA_assay, layer = layer_nm)
                all_count_cells <- c(all_count_cells, colnames(mat_i))
                all_count_genes <- union(all_count_genes, rownames(mat_i))
                vals_i <- sample_matrix_values(mat_i)
                if (length(vals_i) > 0) {
                  sampled_nonnegative <- c(sampled_nonnegative, all(is.finite(vals_i) & vals_i >= 0))
                  sampled_integer <- c(sampled_integer, all(abs(vals_i - round(vals_i)) < 1e-08))
                }
                rm(mat_i, vals_i)
                invisible(gc(FALSE))
            }
            unique_count_cells <- unique(all_count_cells)
            counts_cell_coverage <- length(unique_count_cells)
            counts_feature_count <- length(all_count_genes)
            counts_duplicate_cell_membership <- length(all_count_cells) - length(unique_count_cells)
            values_nonnegative <- if (length(sampled_nonnegative) > 0) {
                all(sampled_nonnegative)
            }
            else {
                NA
            }
            values_integer_like <- if (length(sampled_integer) > 0) {
                all(sampled_integer)
            }
            else {
                NA
            }
            collective_cell_match <- setequal(unique_count_cells, colnames(obj))
            rna_counts_available <- collective_cell_match && counts_duplicate_cell_membership == 0
            counts_qc <- rbind(counts_qc, data.frame(metric = c("RNA_count_layer_structure", "RNA_count_layer_names", 
                "number_of_count_layers", "collective_unique_count_cells", "object_cells", "duplicate_cell_membership_across_count_layers", 
                "collective_count_cell_set_matches_object", "union_count_features", "sampled_counts_nonnegative", 
                "sampled_counts_integer_like", "Step84_action"), value = c(count_layer_type, paste(count_layer_names, 
                collapse = "; "), length(count_layer_names), counts_cell_coverage, n_cells, counts_duplicate_cell_membership, 
                collective_cell_match, counts_feature_count, values_nonnegative, values_integer_like, 
                paste0("QC only. Count layers were NOT joined. ", "Export strategy will be decided in the next step.")), 
                stringsAsFactors = FALSE))
            rm(all_count_cells, all_count_genes, unique_count_cells)
            invisible(gc(FALSE))
        }
        else {
            count_layer_type <- "no_count_layer_detected"
            counts_qc <- rbind(counts_qc, data.frame(metric = "RNA_count_layer_structure", value = count_layer_type, 
                stringsAsFactors = FALSE))
        }
    }
    else {
        counts_mat <- tryCatch(GetAssayData(object = obj, assay = "RNA", slot = "counts"), error = function(e) NULL)
        if (!is.null(counts_mat) && nrow(counts_mat) > 0 && ncol(counts_mat) > 0) {
            count_layer_type <- "legacy_counts_slot"
            rna_counts_available <- TRUE
            counts_feature_count <- nrow(counts_mat)
            counts_cell_coverage <- ncol(counts_mat)
            counts_duplicate_cell_membership <- 0L
            vals <- sample_matrix_values(counts_mat)
            if (length(vals) > 0) {
                values_nonnegative <- all(is.finite(vals) & vals >= 0)
                values_integer_like <- all(abs(vals - round(vals)) < 1e-08)
            }
            counts_qc <- rbind(counts_qc, data.frame(metric = c("RNA_count_layer_structure", "RNA_counts_n_features", 
                "RNA_counts_n_cells", "RNA_counts_cell_set_matches_object", "RNA_counts_cell_order_matches_object", 
                "sampled_counts_nonnegative", "sampled_counts_integer_like"), value = c(count_layer_type, 
                counts_feature_count, counts_cell_coverage, setequal(colnames(counts_mat), colnames(obj)), 
                identical(colnames(counts_mat), colnames(obj)), values_nonnegative, values_integer_like), 
                stringsAsFactors = FALSE))
            rm(counts_mat, vals)
            invisible(gc(FALSE))
        }
    }
}

write.csv(counts_qc, f_counts_qc, row.names = FALSE)

if (rna_exists) {
    RNA_genes <- rownames(obj[["RNA"]])
} else {
    RNA_genes <- character(0)
}

n_gene_names <- length(RNA_genes)

duplicate_genes <- unique(RNA_genes[duplicated(RNA_genes)])

empty_gene_names <- sum(is.na(RNA_genes) | trimws(RNA_genes) == "")

ensembl_mouse_like <- grepl("^ENSMUSG", RNA_genes, ignore.case = FALSE)

ensembl_like <- grepl("^ENS[A-Z]*G", RNA_genes, ignore.case = FALSE)

gene_summary <- data.frame(metric = c("RNA_assay_n_features", "n_duplicate_gene_names", "n_empty_gene_names",
  "n_ENSMUSG_like", "fraction_ENSMUSG_like", "n_any_Ensembl_gene_like", "fraction_any_Ensembl_gene_like",
  "first_30_gene_names"), value = c(n_gene_names, length(duplicate_genes), empty_gene_names, sum(ensembl_mouse_like),
  ifelse(n_gene_names > 0, sum(ensembl_mouse_like)/n_gene_names, NA), sum(ensembl_like), ifelse(n_gene_names >
    0, sum(ensembl_like)/n_gene_names, NA), paste(head(RNA_genes, 30), collapse = "; ")), stringsAsFactors = FALSE)

write.csv(gene_summary, f_gene_summary, row.names = FALSE)

duplicate_gene_table <- data.frame(gene = duplicate_genes, stringsAsFactors = FALSE)

write.csv(duplicate_gene_table, f_duplicate_genes, row.names = FALSE)

meta_export$cell <- rownames(meta_export)

meta_export <- meta_export[, c("cell", setdiff(colnames(meta_export), "cell")), drop = FALSE]

write_csv_gz(meta_export, f_metadata_export)

reduction_names <- Reductions(obj)

reduction_audit <- data.frame(reduction = reduction_names, selected_for_export = FALSE, reason = "",
  stringsAsFactors = FALSE)

preferred_umap <- "umap_myeloid_harmony30"

umap_candidates <- grep("umap", reduction_names, value = TRUE, ignore.case = TRUE)

selected_umap <- NA_character_

if (preferred_umap %in% reduction_names) {
    selected_umap <- preferred_umap
    reduction_audit$selected_for_export[reduction_audit$reduction == selected_umap] <- TRUE
    reduction_audit$reason[reduction_audit$reduction == selected_umap] <- "exact reduction name found"
} else if (length(umap_candidates) == 1) {
    selected_umap <- umap_candidates
    reduction_audit$selected_for_export[reduction_audit$reduction == selected_umap] <- TRUE
    reduction_audit$reason[reduction_audit$reduction == selected_umap] <- paste0("only one UMAP-like reduction found; ", 
        "exported for audit but verify before Figure 6 plotting")
} else if (length(umap_candidates) > 1) {
    reduction_audit$reason[reduction_audit$reduction %in% umap_candidates] <- paste0("multiple UMAP-like reductions found; ", 
        "none selected automatically because exact name was absent")
}

write.csv(reduction_audit, f_reduction_audit, row.names = FALSE)

umap_exported <- FALSE
umap_n_cells <- NA_integer_

if (!is.na(selected_umap)) {
    emb <- Embeddings(obj, reduction = selected_umap)
    umap_df <- as.data.frame(emb)
    umap_df$cell <- rownames(umap_df)
    umap_df <- umap_df[, c("cell", setdiff(colnames(umap_df), "cell")), drop = FALSE]
    umap_n_cells <- nrow(umap_df)
    write_csv_gz(umap_df, f_umap_export)
    umap_exported <- TRUE
    rm(emb, umap_df)
    invisible(gc(FALSE))
}

input_info <- file.info(input_rds)

manifest <- data.frame(item = c("frozen_input_RDS", "input_RDS_size_GB", "input_RDS_last_modified", "state_column",
  "formal_target_name", "target_cells_after_export_mapping", "time_column", "dataset_column_candidates",
  "sample_column_candidates", "RNA_assay_present", "RNA_count_structure", "RNA_count_layers", "selected_frozen_UMAP"),
  value = c(normalizePath(input_rds, winslash = "/", mustWork = TRUE), sprintf("%.3f", input_info$size/1024^3),
    as.character(input_info$mtime), state_col, formal_target, target_n, ifelse(is.na(time_col), "NOT FOUND",
      time_col), ifelse(length(dataset_candidates) == 0, "none detected", paste(dataset_candidates,
      collapse = "; ")), ifelse(length(sample_candidates) == 0, "none detected", paste(sample_candidates,
      collapse = "; ")), rna_exists, count_layer_type, ifelse(length(count_layer_names) == 0, "none detected",
      paste(count_layer_names, collapse = "; ")), ifelse(is.na(selected_umap), "NOT AUTOMATICALLY SELECTED",
      selected_umap)), stringsAsFactors = FALSE)

write.csv(manifest, f_manifest, row.names = FALSE)

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("myeloid cell number", if (n_cells == expected_n_cells) "PASS" else "FAIL",
  paste0("Observed = ", n_cells, "; expected = ", expected_n_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("myeloid state number", if (n_states == expected_n_states) "PASS" else "FAIL",
  paste0("Observed = ", n_states, "; expected = ", expected_n_states))

qc_list[[length(qc_list) + 1]] <- make_qc("myeloid_state column", if (identical(state_col, "myeloid_state")) "PASS" else "FAIL",
  paste0("Detected column = ", state_col))

qc_list[[length(qc_list) + 1]] <- make_qc("Clec7a target state", if (target_present_after_mapping) "PASS" else "FAIL",
  paste0("Formal exported target = ", formal_target, "; n = ", target_n))

qc_list[[length(qc_list) + 1]] <- make_qc("RNA assay present", if (rna_exists) "PASS" else "FAIL", paste0("Available assays = ",
  paste(assay_names, collapse = ", ")))

qc_list[[length(qc_list) + 1]] <- make_qc("RNA count data accessible", if (rna_counts_available) "PASS" else "FAIL",
  paste0("Count structure = ", count_layer_type, "; collective cells = ", counts_cell_coverage))

qc_list[[length(qc_list) + 1]] <- make_qc("RNA count values non-negative", if (isTRUE(values_nonnegative)) {
  "PASS"
} else if (is.na(values_nonnegative)) {
  "WARN"
} else {
  "FAIL"
}, paste0("Sampled-value result = ", values_nonnegative))

qc_list[[length(qc_list) + 1]] <- make_qc("RNA count values integer-like", if (isTRUE(values_integer_like)) {
  "PASS"
} else if (is.na(values_integer_like)) {
  "WARN"
} else {
  "WARN"
}, paste0("Sampled-value result = ", values_integer_like, ". This is a diagnostic rather than automatic data rewriting."))

qc_list[[length(qc_list) + 1]] <- make_qc("Duplicate RNA gene names", if (length(duplicate_genes) ==
  0) "PASS" else "FAIL", paste0("Number of duplicated gene names = ", length(duplicate_genes)))

qc_list[[length(qc_list) + 1]] <- make_qc("Empty RNA gene names", if (empty_gene_names == 0) "PASS" else "FAIL",
  paste0("Number of empty gene names = ", empty_gene_names))

qc_list[[length(qc_list) + 1]] <- make_qc("UMAP identified", if (preferred_umap %in% reduction_names) {
  "PASS"
} else if (!is.na(selected_umap)) {
  "WARN"
} else {
  "WARN"
}, paste0("Preferred = ", preferred_umap, "; selected = ", ifelse(is.na(selected_umap), "none", selected_umap)))

qc_list[[length(qc_list) + 1]] <- make_qc("time metadata column", if (!is.na(time_col)) "PASS" else "WARN",
  paste0("Detected = ", ifelse(is.na(time_col), "none", time_col)))

final_qc <- do.call(rbind, qc_list)

rownames(final_qc) <- NULL

write.csv(final_qc, f_final_qc, row.names = FALSE)

n_fail <- sum(final_qc$status == "FAIL")
n_warn <- sum(final_qc$status == "WARN")
n_pass <- sum(final_qc$status == "PASS")

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

all_outputs <- c(f_object_summary, f_metadata_fields, f_state_counts, f_state_name_audit, f_assay_audit,
  f_counts_qc, f_gene_summary, f_duplicate_genes, f_metadata_export, f_reduction_audit, f_manifest,
  f_final_qc, f_report, f_session)

if (file.exists(f_umap_export)) {
  all_outputs <- c(all_outputs, f_umap_export)
}

output_index <- data.frame(file = basename(all_outputs), full_path = normalizePath(all_outputs, winslash = "/",
  mustWork = FALSE), exists = file.exists(all_outputs), size_MB = ifelse(file.exists(all_outputs),
  file.info(all_outputs)$size/1024^2, NA_real_), stringsAsFactors = FALSE)

write.csv(output_index, f_output_index, row.names = FALSE)

cat("  WARN =", n_warn, "\n")
cat("  FAIL =", n_fail, "\n")

if (n_fail > 0) {
    stop("\nStep 84 generated all QC files, but formal QC contains ", n_fail, " FAIL item(s).\n", "Do NOT start GRNBoost2 until these are reviewed.\n", 
        "See:\n", f_final_qc)
}

cat("\nStep 84 passed formal blocking QC.\n", "Review warnings before designing Step 85 GRNBoost2 input.\n",
  sep = "")
