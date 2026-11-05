source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
N_CORES <- 1
TARGET_PARTITION <- "1"
N_BINS <- 50

RUN_GRAPH_TEST <- TRUE

GENE_GROUPS <- list(Homeostatic = c("P2ry12", "Tmem119", "Sall1", "Cx3cr1"), Lipid_handling = c("Clec7a",
  "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1"), Inflammatory = c("Ccl12", "Tnf", "Il1b"), Proliferative = c("Mki67",
  "Ccnb2"))

KEY_GENES <- unique(unlist(GENE_GROUPS, use.names = FALSE))

STATE_LEVELS <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia",
  "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage",
  "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

set.seed(SEED)

required_pkgs <- c("monocle3", "SingleCellExperiment", "Matrix", "ggplot2", "dplyr", "tidyr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(monocle3)
  library(SingleCellExperiment)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "25_Monocle3_all_myeloid_trajectory", "25_Monocle3_all_myeloid_trajectory_3scRNA_v03.rds")

if (file.exists(expected_input)) {
    in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
} else {
    hits <- list.files(project_dir, pattern = "^25_Monocle3_all_myeloid_trajectory_3scRNA_v03\\.rds$", 
        recursive = TRUE, full.names = TRUE)
    if (length(hits) == 0) {
        stop("Cannot find Step 25 v03 Monocle3 RDS.")
    }
    in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "26_Monocle3_keygene_pseudotime_dynamics")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "26_Monocle3_keygene_pseudotime_dynamics_v01.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

cds <- readRDS(in_file)

if (!inherits(cds, "cell_data_set")) {
  stop("Input is not a Monocle3 cell_data_set.")
}

required_meta <- c("myeloid_state", "source_dataset", "timepoint_main")

missing_meta <- setdiff(required_meta, colnames(colData(cds)))

if (length(missing_meta) > 0) {
  stop("Missing metadata: ", paste(missing_meta, collapse = ", "))
}

pt_all <- monocle3::pseudotime(cds)

part_all <- monocle3::partitions(cds)

part_all <- as.character(part_all[colnames(cds)])

pt_all <- as.numeric(pt_all[colnames(cds)])

if (!TARGET_PARTITION %in% unique(part_all)) {
  stop("Target partition ", TARGET_PARTITION, " not found.")
}

p1_cells <- colnames(cds)[part_all == TARGET_PARTITION & is.finite(pt_all)]

cds_p1 <- cds[, p1_cells]

pt_p1 <- monocle3::pseudotime(cds_p1)

pt_p1 <- as.numeric(pt_p1[colnames(cds_p1)])

if (length(pt_p1) != ncol(cds_p1) || any(!is.finite(pt_p1))) {
  stop("Partition-1 pseudotime is missing/non-finite after subsetting.")
}

genes_present <- intersect(KEY_GENES, rownames(cds_p1))

genes_missing <- setdiff(KEY_GENES, genes_present)

if (length(genes_missing) > 0) {
  log_msg("Missing genes: ", paste(genes_missing, collapse = ", "))
}

if (length(genes_present) == 0) {
  stop("None of the pre-specified key genes are present.")
}

counts_key <- SingleCellExperiment::counts(cds_p1[genes_present, ])

detect_n <- Matrix::rowSums(counts_key > 0)

detect_summary <- data.frame(gene = genes_present, gene_group = vapply(genes_present, function(g) {
  nm <- names(GENE_GROUPS)[vapply(GENE_GROUPS, function(x) g %in% x, logical(1))]
  if (length(nm) == 0) {
    NA_character_
  } else {
    nm[1]
  }
}, FUN.VALUE = character(1)), n_detected_cells = as.integer(detect_n[genes_present]), pct_detected_cells = 100 *
  as.integer(detect_n[genes_present])/ncol(cds_p1), stringsAsFactors = FALSE)

write.csv(detect_summary, file.path(out_dir, "26_keygene_detection_summary_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

capture.output(print(detect_summary), file = log_file, append = TRUE)

gene_order <- KEY_GENES[KEY_GENES %in% genes_present]

cds_key <- cds_p1[gene_order, ]

colData(cds_key)$myeloid_state <- factor(as.character(colData(cds_key)$myeloid_state), levels = STATE_LEVELS)

p_curve <- monocle3::plot_genes_in_pseudotime(cds_key, color_cells_by = "myeloid_state", min_expr = 0.5,
  ncol = 3) + ggplot2::theme_classic(base_size = 9) + ggplot2::theme(legend.position = "right")

ggplot2::ggsave(filename = file.path(out_dir, "26_keygene_pseudotime_curves_v01.pdf"), plot = p_curve,
  width = 12, height = 11, units = "in", limitsize = FALSE)

counts_all_p1 <- SingleCellExperiment::counts(cds_p1)

lib_size <- Matrix::colSums(counts_all_p1)

if (any(!is.finite(lib_size) | lib_size <= 0)) {
  stop("Non-positive library sizes detected in partition 1.")
}

counts_key <- counts_all_p1[gene_order, , drop = FALSE]

scale_factor <- 10000/lib_size

norm_key <- counts_key %*% Matrix::Diagonal(x = scale_factor)

norm_key@x <- log1p(norm_key@x)

rank_pt <- rank(pt_p1, ties.method = "average")

bin_id <- pmin(N_BINS, ceiling(rank_pt/length(rank_pt) * N_BINS))

bin_id[bin_id < 1] <- 1

bin_meta <- data.frame(cell = colnames(cds_p1), pseudotime = pt_p1, bin = bin_id, stringsAsFactors = FALSE)

bin_pt <- bin_meta %>% dplyr::group_by(bin) %>% dplyr::summarise(n_cells = dplyr::n(), pseudotime_mean = mean(pseudotime), 
    pseudotime_median = median(pseudotime), .groups = "drop")

bin_expr_list <- vector("list", length = N_BINS)

for (b in seq_len(N_BINS)) {
    idx <- which(bin_id == b)
    if (length(idx) == 0) {
        next
    }
    mu <- Matrix::rowMeans(norm_key[, idx, drop = FALSE])
    bin_expr_list[[b]] <- data.frame(bin = b, gene = names(mu), mean_log1p_norm_expr = as.numeric(mu), 
        stringsAsFactors = FALSE)
}

bin_expr <- dplyr::bind_rows(bin_expr_list) %>% dplyr::left_join(bin_pt, by = "bin") %>% dplyr::mutate(gene = factor(gene, 
    levels = gene_order))

write.csv(bin_expr, file.path(out_dir, "26_keygene_pseudotime_binned_expression_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

peak_table <- bin_expr %>% dplyr::group_by(gene) %>% dplyr::slice_max(order_by = mean_log1p_norm_expr, 
    n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>% dplyr::transmute(gene = as.character(gene), peak_bin = bin, 
    peak_pseudotime = pseudotime_median, peak_mean_log1p_norm_expr = mean_log1p_norm_expr) %>% dplyr::mutate(gene_group = vapply(gene, 
    function(g) {
        nm <- names(GENE_GROUPS)[vapply(GENE_GROUPS, function(x) g %in% x, logical(1))]
        if (length(nm) == 0) {
            NA_character_
        }
        else {
            nm[1]
        }
    }, FUN.VALUE = character(1))) %>% dplyr::arrange(peak_pseudotime)

write.csv(peak_table, file.path(out_dir, "26_keygene_peak_pseudotime_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")

capture.output(print(peak_table), file = log_file, append = TRUE)

heat_df <- bin_expr %>% dplyr::group_by(gene) %>% dplyr::mutate(z = as.numeric(scale(mean_log1p_norm_expr))) %>% 
    dplyr::ungroup() %>% dplyr::mutate(gene = factor(gene, levels = rev(gene_order)))

heat_df$z[!is.finite(heat_df$z)] <- 0

p_heat <- ggplot2::ggplot(heat_df, ggplot2::aes(x = pseudotime_median, y = gene, fill = z)) + ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, limits = c(-2.5,
    2.5), oob = scales::squish) + ggplot2::labs(title = "Key-gene dynamics along the partition-1 Monocle3 trajectory",
  x = "Monocle3 pseudotime", y = NULL, fill = "Row z-score") + ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(axis.ticks.y = ggplot2::element_blank())

ggplot2::ggsave(filename = file.path(out_dir, "26_keygene_pseudotime_heatmap_v01.pdf"), plot = p_heat,
  width = 9, height = 6, units = "in")

if (RUN_GRAPH_TEST) {
    counts_p1 <- SingleCellExperiment::counts(cds_p1)
    n_detected <- Matrix::rowSums(counts_p1 > 0)
    min_cells_graph_test <- max(50, ceiling(0.005 * ncol(cds_p1)))
    keep_graph_genes <- names(n_detected)[n_detected >= min_cells_graph_test]
    cds_gt <- cds_p1[keep_graph_genes, ]
    graph_test_res <- monocle3::graph_test(cds_gt, neighbor_graph = "principal_graph", cores = N_CORES)
    graph_test_res <- as.data.frame(graph_test_res)
    graph_test_res$gene <- rownames(graph_test_res)
    if ("q_value" %in% colnames(graph_test_res)) {
        graph_test_res <- graph_test_res %>% dplyr::arrange(q_value, dplyr::desc(morans_I))
    }
    write.csv(graph_test_res, file.path(out_dir, "26_graph_test_partition1_all_v01.csv"), row.names = FALSE, 
        fileEncoding = "UTF-8")
    graph_sig <- graph_test_res %>% dplyr::filter(is.finite(q_value), q_value < 0.05, is.finite(morans_I), 
        morans_I > 0.1)
    write.csv(graph_sig, file.path(out_dir, "26_graph_test_partition1_significant_v01.csv"), row.names = FALSE, 
        fileEncoding = "UTF-8")
    graph_key <- graph_test_res %>% dplyr::filter(gene %in% gene_order) %>% dplyr::mutate(gene = factor(gene, 
        levels = gene_order)) %>% dplyr::arrange(gene)
    write.csv(graph_key, file.path(out_dir, "26_graph_test_keygenes_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    saveRDS(graph_test_res, file.path(out_dir, "26_Monocle3_partition1_graph_test_v01.rds"), compress = FALSE)
} else {
    log_msg("graph_test skipped because RUN_GRAPH_TEST = FALSE.")
}

saveRDS(cds_p1, file.path(out_dir, "26_Monocle3_partition1_pseudotime_3scRNA_v01.rds"), compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)
