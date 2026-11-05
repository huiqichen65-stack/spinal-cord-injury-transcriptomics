source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
N_CORES <- 1
NUM_DIM <- 30
MIN_HOME_CONTROL_FOR_ROOT <- 10
UMAP_REDUCTION <- "umap_myeloid_harmony30"

STATE_LEVELS <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia",
  "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage",
  "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

TIME_LEVELS <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")

TIME_RANK <- c(Control = 0, `4 h` = 4/24, `1 dpi` = 1, `3 dpi` = 3, `7 dpi` = 7, `14 dpi` = 14, `28 dpi` = 28,
  `38 dpi` = 38)

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "monocle3", "ggplot2", "dplyr", "igraph")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them before running Step 25.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(monocle3)
  library(ggplot2)
  library(dplyr)
  library(igraph)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "24_CytoTRACE2_myeloid", "24_CytoTRACE2_myeloid_3scRNA_v04.rds")

if (!file.exists(expected_input)) {
  hits <- list.files(project_dir, pattern = "^24_CytoTRACE2_myeloid_3scRNA_v04\\.rds$", recursive = TRUE,
    full.names = TRUE)
  if (length(hits) == 0) {
    stop("Cannot find Step 24 v04 RDS:\n", "24_CytoTRACE2_myeloid_3scRNA_v04.rds")
  }
  in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
} else {
  in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "25_Monocle3_all_myeloid_trajectory")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "25_Monocle3_all_myeloid_trajectory_v03.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

required_meta <- c("myeloid_state", "source_dataset", "timepoint_main")

missing_meta <- setdiff(required_meta, colnames(obj[[]]))

if (length(missing_meta) > 0) {
  stop("Missing metadata: ", paste(missing_meta, collapse = ", "))
}

if (!UMAP_REDUCTION %in% names(obj@reductions)) {
  stop("Missing Seurat reduction: ", UMAP_REDUCTION)
}

obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) obj)

DefaultAssay(obj) <- "RNA"

if (ncol(obj) != 58533) {
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

lib_size <- Matrix::colSums(counts)

zero_count_flag <- !is.finite(lib_size) | lib_size <= 0

n_zero <- sum(zero_count_flag)

if (n_zero > 0) {
    zero_cells <- colnames(counts)[zero_count_flag]
    zero_meta <- obj[[]][zero_cells, , drop = FALSE]
    zero_export <- data.frame(cell = zero_cells, source_dataset = if ("source_dataset" %in% colnames(zero_meta)) 
        as.character(zero_meta$source_dataset)
    else NA_character_, myeloid_state = if ("myeloid_state" %in% colnames(zero_meta)) 
        as.character(zero_meta$myeloid_state)
    else NA_character_, timepoint_main = if ("timepoint_main" %in% colnames(zero_meta)) 
        as.character(zero_meta$timepoint_main)
    else NA_character_, total_RNA_counts = as.numeric(lib_size[zero_cells]), stringsAsFactors = FALSE)
    write.csv(zero_export, file.path(out_dir, "25_Monocle3_zero_count_cells_removed_v03.csv"), row.names = FALSE, 
        fileEncoding = "UTF-8")
    keep_cells <- colnames(counts)[!zero_count_flag]
    obj <- subset(obj, cells = keep_cells)
    counts <- counts[, keep_cells, drop = FALSE]
} else {
}

lib_size_after <- Matrix::colSums(counts)

if (any(!is.finite(lib_size_after) | lib_size_after <= 0)) {
  stop("Zero/non-finite RNA-count cells remain after filtering.")
}

cell_metadata <- obj[[]]
cell_metadata <- cell_metadata[colnames(counts), , drop = FALSE]

gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts), stringsAsFactors = FALSE)

cds <- monocle3::new_cell_data_set(expression_data = counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata)

cds <- monocle3::estimate_size_factors(cds)

sf <- monocle3::size_factors(cds)

if (is.null(sf) || length(sf) != ncol(cds) || any(!is.finite(sf)) || any(sf <= 0)) {
  stop("Monocle3 size-factor calculation failed after zero-count filtering.")
}

rm(counts, lib_size, lib_size_after, sf)
gc(verbose = FALSE)

set.seed(SEED)

cds <- monocle3::preprocess_cds(cds, method = "PCA", num_dim = NUM_DIM, norm_method = "log")

set.seed(SEED)

cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", umap.fast_sgd = FALSE,
  cores = N_CORES)

seurat_umap <- Embeddings(obj, reduction = UMAP_REDUCTION)

if (!setequal(rownames(seurat_umap), colnames(cds))) {
  stop("Cell names differ between Seurat UMAP and Monocle3 CDS.")
}

seurat_umap <- seurat_umap[colnames(cds), , drop = FALSE]

colnames(seurat_umap) <- c("UMAP_1", "UMAP_2")

SingleCellExperiment::reducedDim(cds, "UMAP") <- seurat_umap

set.seed(SEED)

cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")

partition_vec <- cds@clusters[["UMAP"]]$partitions
partition_vec <- partition_vec[colnames(cds)]

partition_summary <- data.frame(partition = as.character(partition_vec), myeloid_state = as.character(colData(cds)$myeloid_state), 
    timepoint_main = as.character(colData(cds)$timepoint_main), stringsAsFactors = FALSE) %>% dplyr::count(partition, 
    myeloid_state, name = "n_cells") %>% dplyr::arrange(partition, dplyr::desc(n_cells))

write.csv(partition_summary, file.path(out_dir, "25_Monocle3_partition_state_summary_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

set.seed(SEED)

cds <- monocle3::learn_graph(cds, use_partition = TRUE, learn_graph_control = list(euclidean_distance_ratio = 0.8))

closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex

closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])

closest_idx <- as.integer(closest_vertex[, 1])

graph_obj <- monocle3::principal_graph(cds)[["UMAP"]]

graph_node_names <- igraph::V(graph_obj)$name

if (any(closest_idx < 1 | closest_idx > length(graph_node_names))) {
  stop("Principal-node index mapping is invalid.")
}

closest_node <- graph_node_names[closest_idx]

node_df <- data.frame(cell = colnames(cds), principal_node = closest_node, partition = as.character(partition_vec),
  myeloid_state = as.character(colData(cds)$myeloid_state), timepoint_main = as.character(colData(cds)$timepoint_main),
  stringsAsFactors = FALSE)

if ("CytoTRACE2_Score" %in% colnames(colData(cds))) {
  node_df$CytoTRACE2_Score <- as.numeric(colData(cds)$CytoTRACE2_Score)
} else {
  node_df$CytoTRACE2_Score <- NA_real_
}

node_df$time_rank <- unname(TIME_RANK[node_df$timepoint_main])

root_rows <- list()
root_nodes <- character(0)

for (part in sort(unique(node_df$partition))) {
    x <- node_df[node_df$partition == part, , drop = FALSE]
    home_ctrl <- x[x$myeloid_state == "Homeostatic Microglia" & x$timepoint_main == "Control", , drop = FALSE]
    if (nrow(home_ctrl) >= MIN_HOME_CONTROL_FOR_ROOT) {
        tab <- sort(table(home_ctrl$principal_node), decreasing = TRUE)
        root_node <- names(tab)[1]
        rule <- "Control + Homeostatic Microglia"
        n_rule_cells <- as.integer(tab[1])
    }
    else {
        valid_time <- x[is.finite(x$time_rank), , drop = FALSE]
        if (nrow(valid_time) == 0) {
            stop("Partition ", part, " has no recognized timepoint for fallback root selection.")
        }
        earliest_rank <- min(valid_time$time_rank, na.rm = TRUE)
        earliest_cells <- valid_time[valid_time$time_rank == earliest_rank, , drop = FALSE]
        tab <- sort(table(earliest_cells$principal_node), decreasing = TRUE)
        root_node <- names(tab)[1]
        earliest_label <- unique(earliest_cells$timepoint_main)[1]
        rule <- paste0("Earliest available time: ", earliest_label)
        n_rule_cells <- as.integer(tab[1])
    }
    root_nodes <- c(root_nodes, root_node)
    root_cells <- x[x$principal_node == root_node, , drop = FALSE]
    root_rows[[length(root_rows) + 1]] <- data.frame(partition = part, root_pr_node = root_node, selection_rule = rule, 
        n_rule_cells_on_root_node = n_rule_cells, n_all_cells_on_root_node = nrow(root_cells), pct_homeostatic_on_root_node = 100 * 
            mean(root_cells$myeloid_state == "Homeostatic Microglia"), pct_control_on_root_node = 100 * 
            mean(root_cells$timepoint_main == "Control"), mean_CytoTRACE2_on_root_node = if (any(is.finite(root_cells$CytoTRACE2_Score))) {
            mean(root_cells$CytoTRACE2_Score, na.rm = TRUE)
        }
        else {
            NA_real_
        }, stringsAsFactors = FALSE)
}

root_table <- dplyr::bind_rows(root_rows)

root_nodes <- unique(root_nodes)

write.csv(root_table, file.path(out_dir, "25_Monocle3_root_node_selection_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

capture.output(print(root_table), file = log_file, append = TRUE)

cds <- monocle3::order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_nodes)

pt <- monocle3::pseudotime(cds)

n_finite_pt <- sum(is.finite(pt))

if (n_finite_pt < 0.95 * length(pt)) {
  log_msg("WARNING: <95% cells have finite pseudotime. ", "Inspect partition/root diagnostics before interpretation.")
}

cell_pt <- data.frame(cell = colnames(cds), partition = as.character(partition_vec), principal_node = closest_node,
  myeloid_state = as.character(colData(cds)$myeloid_state), source_dataset = as.character(colData(cds)$source_dataset),
  timepoint_main = as.character(colData(cds)$timepoint_main), pseudotime = as.numeric(pt), stringsAsFactors = FALSE)

if ("CytoTRACE2_Score" %in% colnames(colData(cds))) {
  cell_pt$CytoTRACE2_Score <- as.numeric(colData(cds)$CytoTRACE2_Score)
}

write.csv(cell_pt, file.path(out_dir, "25_Monocle3_cell_pseudotime_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

state_pt <- cell_pt %>% dplyr::filter(is.finite(pseudotime)) %>% dplyr::group_by(myeloid_state) %>% dplyr::summarise(n_cells = dplyr::n(), 
    mean_pseudotime = mean(pseudotime), median_pseudotime = median(pseudotime), q25 = quantile(pseudotime, 
        0.25), q75 = quantile(pseudotime, 0.75), .groups = "drop")

write.csv(state_pt, file.path(out_dir, "25_Monocle3_state_pseudotime_summary_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

time_pt <- cell_pt %>% dplyr::filter(is.finite(pseudotime)) %>% dplyr::group_by(timepoint_main) %>% dplyr::summarise(n_cells = dplyr::n(), 
    mean_pseudotime = mean(pseudotime), median_pseudotime = median(pseudotime), q25 = quantile(pseudotime, 
        0.25), q75 = quantile(pseudotime, 0.75), .groups = "drop")

write.csv(time_pt, file.path(out_dir, "25_Monocle3_time_pseudotime_summary_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

node_summary <- node_df %>% dplyr::group_by(partition, principal_node) %>% dplyr::summarise(n_cells = dplyr::n(), 
    n_homeostatic = sum(myeloid_state == "Homeostatic Microglia"), n_control = sum(timepoint_main == 
        "Control"), n_homeostatic_control = sum(myeloid_state == "Homeostatic Microglia" & timepoint_main == 
        "Control"), mean_CytoTRACE2 = mean(CytoTRACE2_Score, na.rm = TRUE), .groups = "drop")

write.csv(node_summary, file.path(out_dir, "25_Monocle3_principal_node_summary_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

state_plot_df <- cell_pt %>% dplyr::filter(is.finite(pseudotime)) %>% dplyr::mutate(myeloid_state = factor(myeloid_state, 
    levels = STATE_LEVELS))

time_plot_df <- cell_pt %>% dplyr::filter(is.finite(pseudotime), timepoint_main %in% TIME_LEVELS) %>% 
    dplyr::mutate(timepoint_main = factor(timepoint_main, levels = TIME_LEVELS))

pdf(file.path(out_dir, "25_Monocle3_diagnostics_v03.pdf"), width = 11, height = 8.5, onefile = TRUE)

print(monocle3::plot_cells(cds, reduction_method = "UMAP", color_cells_by = "myeloid_state", show_trajectory_graph = FALSE,
  label_cell_groups = FALSE) + ggtitle("All retained myeloid states on imported integrated UMAP"))

print(monocle3::plot_cells(cds, reduction_method = "UMAP", color_cells_by = "partition", show_trajectory_graph = FALSE,
  label_cell_groups = TRUE) + ggtitle("Monocle3 partitions before trajectory interpretation"))

print(monocle3::plot_cells(cds, reduction_method = "UMAP", color_cells_by = "myeloid_state", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE, label_roots = TRUE) + ggtitle("Learned principal graph colored by myeloid state"))

print(monocle3::plot_cells(cds, reduction_method = "UMAP", color_cells_by = "timepoint_main", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE, label_roots = TRUE) + ggtitle("Learned principal graph colored by injury time"))

print(monocle3::plot_cells(cds, reduction_method = "UMAP", color_cells_by = "pseudotime", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE, label_roots = TRUE) + ggtitle("Monocle3 pseudotime"))

print(ggplot(state_plot_df, aes(x = myeloid_state, y = pseudotime)) + geom_violin(scale = "width", trim = TRUE,
  linewidth = 0.25) + geom_boxplot(width = 0.1, outlier.shape = NA, linewidth = 0.25) + labs(title = "Pseudotime distribution across myeloid states",
  x = NULL, y = "Pseudotime") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 35,
  hjust = 1)))

print(ggplot(time_plot_df, aes(x = timepoint_main, y = pseudotime)) + geom_violin(scale = "width", trim = TRUE,
  linewidth = 0.25) + geom_boxplot(width = 0.1, outlier.shape = NA, linewidth = 0.25) + labs(title = "Pseudotime distribution across injury time points",
  x = NULL, y = "Pseudotime") + theme_classic(base_size = 10))

if ("CytoTRACE2_Score" %in% colnames(cell_pt)) {
    ct_plot <- cell_pt %>% dplyr::filter(is.finite(pseudotime), is.finite(CytoTRACE2_Score))
    if (nrow(ct_plot) > 20000) {
        set.seed(SEED)
        ct_plot <- ct_plot[sample(seq_len(nrow(ct_plot)), 20000), , drop = FALSE]
    }
    print(ggplot(ct_plot, aes(x = pseudotime, y = CytoTRACE2_Score)) + geom_point(alpha = 0.15, size = 0.35) + 
        geom_smooth(method = "loess", se = FALSE, linewidth = 0.6) + labs(title = "CytoTRACE2 vs Monocle3 pseudotime (independent diagnostics)", 
        x = "Monocle3 pseudotime", y = "CytoTRACE2 absolute score") + theme_classic(base_size = 10))
}

dev.off()

saveRDS(cds, file.path(out_dir, "25_Monocle3_all_myeloid_trajectory_3scRNA_v03.rds"), compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)
