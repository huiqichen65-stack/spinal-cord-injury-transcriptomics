source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
N_BINS <- 50
K_CANDIDATES <- 4:6
MAX_CLUSTER_GENES <- 3000
MIN_MORAN_I <- 0.1
MAX_Q <- 0.05
MIN_DETECTED_FRAC <- 0.005

KEY_GENES <- c("P2ry12", "Tmem119", "Sall1", "Cx3cr1", "Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1",
  "Ccl12", "Tnf", "Il1b", "Mki67", "Ccnb2")

set.seed(SEED)

required_pkgs <- c("monocle3", "SingleCellExperiment", "Matrix", "ggplot2", "dplyr", "tidyr", "cluster")

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
  library(cluster)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

input_dir <- file.path(project_dir, "results", "26_Monocle3_keygene_pseudotime_dynamics")

cds_file <- file.path(input_dir, "26_Monocle3_partition1_pseudotime_3scRNA_v01.rds")

gt_file <- file.path(input_dir, "26_graph_test_partition1_all_v01.csv")

if (!file.exists(cds_file)) {
  stop("Cannot find: ", cds_file)
}

if (!file.exists(gt_file)) {
  stop("Cannot find: ", gt_file)
}

out_dir <- file.path(project_dir, "results", "27_Monocle3_trajectory_gene_modules")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "27_trajectory_gene_modules_v01.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

cds <- readRDS(cds_file)

if (!inherits(cds, "cell_data_set")) {
  stop("Input RDS is not a Monocle3 cell_data_set.")
}

gt <- read.csv(gt_file, stringsAsFactors = FALSE, check.names = FALSE)

required_gt <- c("gene", "q_value", "morans_I")

missing_gt <- setdiff(required_gt, colnames(gt))

if (length(missing_gt) > 0) {
  stop("graph_test table missing column(s): ", paste(missing_gt, collapse = ", "))
}

pt <- monocle3::pseudotime(cds)
pt <- as.numeric(pt[colnames(cds)])

if (length(pt) != ncol(cds) || any(!is.finite(pt))) {
  stop("Non-finite or missing pseudotime in partition-1 CDS.")
}

counts_all <- SingleCellExperiment::counts(cds)

detect_n <- Matrix::rowSums(counts_all > 0)

min_detected_cells <- max(50, ceiling(MIN_DETECTED_FRAC * ncol(cds)))

detect_df <- data.frame(gene = names(detect_n), n_detected_cells = as.integer(detect_n), stringsAsFactors = FALSE)

gt2 <- gt %>% dplyr::left_join(detect_df, by = "gene") %>% dplyr::filter(is.finite(q_value), q_value < 
    MAX_Q, is.finite(morans_I), morans_I > MIN_MORAN_I, !is.na(n_detected_cells), n_detected_cells >= 
    min_detected_cells) %>% dplyr::arrange(dplyr::desc(morans_I))

if (nrow(gt2) < 100) {
  stop("Too few trajectory-associated genes after filtering: ", nrow(gt2))
}

if (nrow(gt2) > MAX_CLUSTER_GENES) {
    gt_use <- gt2 %>% dplyr::slice_head(n = MAX_CLUSTER_GENES)
} else {
    gt_use <- gt2
}

genes_use <- intersect(gt_use$gene, rownames(cds))

write.csv(gt_use, file.path(out_dir, "27_trajectory_genes_used_for_modules_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

rank_pt <- rank(pt, ties.method = "average")

bin_id <- pmin(N_BINS, ceiling(rank_pt/length(rank_pt) * N_BINS))

bin_id[bin_id < 1] <- 1

bin_meta <- data.frame(cell = colnames(cds), pseudotime = pt, bin = bin_id, stringsAsFactors = FALSE)

bin_pt <- bin_meta %>% dplyr::group_by(bin) %>% dplyr::summarise(n_cells = dplyr::n(), pseudotime_mean = mean(pseudotime), 
    pseudotime_median = median(pseudotime), .groups = "drop")

lib_size <- Matrix::colSums(counts_all)

if (any(!is.finite(lib_size) | lib_size <= 0)) {
  stop("Invalid library size detected.")
}

counts_use <- counts_all[genes_use, , drop = FALSE]

norm_use <- counts_use %*% Matrix::Diagonal(x = 10000/lib_size)

norm_use@x <- log1p(norm_use@x)

bin_mat <- matrix(NA_real_, nrow = length(genes_use), ncol = N_BINS, dimnames = list(genes_use, paste0("Bin",
  seq_len(N_BINS))))

for (b in seq_len(N_BINS)) {
    idx <- which(bin_id == b)
    if (length(idx) == 0) {
        next
    }
    bin_mat[, b] <- Matrix::rowMeans(norm_use[, idx, drop = FALSE])
}

if (anyNA(bin_mat)) {
  stop("At least one pseudotime bin has no cells; cannot construct full temporal matrix.")
}

row_sd <- apply(bin_mat, 1, stats::sd)

keep_var <- is.finite(row_sd) & row_sd > 0

bin_mat <- bin_mat[keep_var, , drop = FALSE]

genes_use <- rownames(bin_mat)

z_mat <- t(scale(t(bin_mat)))

z_mat[!is.finite(z_mat)] <- 0

smooth_mat <- z_mat

if (ncol(z_mat) >= 3) {
    for (j in 2:(ncol(z_mat) - 1)) {
        smooth_mat[, j] <- rowMeans(cbind(z_mat[, j - 1], z_mat[, j], z_mat[, j + 1]))
    }
}

set.seed(SEED)

if (nrow(smooth_mat) > 1500) {
  sil_idx <- seq_len(1500)
} else {
  sil_idx <- seq_len(nrow(smooth_mat))
}

sil_mat <- smooth_mat[sil_idx, , drop = FALSE]

dist_sil <- stats::dist(sil_mat)

k_eval <- vector("list", length(K_CANDIDATES))

for (i in seq_along(K_CANDIDATES)) {
    k <- K_CANDIDATES[i]
    set.seed(SEED)
    km_tmp <- stats::kmeans(sil_mat, centers = k, nstart = 50, iter.max = 100)
    sil <- cluster::silhouette(km_tmp$cluster, dist_sil)
    k_eval[[i]] <- data.frame(k = k, mean_silhouette = mean(sil[, "sil_width"]), stringsAsFactors = FALSE)
}

k_eval_df <- dplyr::bind_rows(k_eval)

best_k <- k_eval_df$k[which.max(k_eval_df$mean_silhouette)]

write.csv(k_eval_df, file.path(out_dir, "27_module_number_assessment_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")

capture.output(print(k_eval_df), file = log_file, append = TRUE)

rm(dist_sil, sil_mat)
gc(verbose = FALSE)

set.seed(SEED)

km <- stats::kmeans(smooth_mat, centers = best_k, nstart = 100, iter.max = 200)

raw_cluster <- km$cluster

peak_bin_gene <- max.col(smooth_mat, ties.method = "first")

gene_peak_pt <- bin_pt$pseudotime_median[peak_bin_gene]

tmp_assign <- data.frame(gene = rownames(smooth_mat), raw_cluster = as.integer(raw_cluster), peak_bin = peak_bin_gene,
  peak_pseudotime = gene_peak_pt, stringsAsFactors = FALSE)

cluster_order <- tmp_assign %>% dplyr::group_by(raw_cluster) %>% dplyr::summarise(median_peak_pseudotime = median(peak_pseudotime, 
    na.rm = TRUE), .groups = "drop") %>% dplyr::arrange(median_peak_pseudotime)

cluster_order$module <- paste0("Module ", seq_len(nrow(cluster_order)))

module_map <- stats::setNames(cluster_order$module, cluster_order$raw_cluster)

assignment <- tmp_assign %>% dplyr::mutate(module = unname(module_map[as.character(raw_cluster)])) %>% 
    dplyr::left_join(gt_use %>% dplyr::select(gene, q_value, morans_I, n_detected_cells), by = "gene") %>% 
    dplyr::arrange(factor(module, levels = paste0("Module ", seq_len(best_k))), peak_pseudotime, dplyr::desc(morans_I))

write.csv(assignment, file.path(out_dir, "27_trajectory_gene_module_assignment_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

module_summary <- assignment %>% dplyr::group_by(module) %>% dplyr::summarise(n_genes = dplyr::n(), median_peak_pseudotime = median(peak_pseudotime, 
    na.rm = TRUE), q25_peak_pseudotime = quantile(peak_pseudotime, 0.25, na.rm = TRUE), q75_peak_pseudotime = quantile(peak_pseudotime, 
    0.75, na.rm = TRUE), median_morans_I = median(morans_I, na.rm = TRUE), .groups = "drop") %>% dplyr::arrange(median_peak_pseudotime)

write.csv(module_summary, file.path(out_dir, "27_trajectory_module_summary_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

capture.output(print(module_summary), file = log_file, append = TRUE)

key_membership <- data.frame(gene = KEY_GENES, stringsAsFactors = FALSE) %>% dplyr::left_join(assignment %>% 
    dplyr::select(gene, module, peak_bin, peak_pseudotime, morans_I, q_value), by = "gene")

write.csv(key_membership, file.path(out_dir, "27_trajectory_module_keygene_membership_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

gene_order_df <- assignment %>% dplyr::arrange(factor(module, levels = paste0("Module ", seq_len(best_k))), 
    peak_pseudotime, dplyr::desc(morans_I)) %>% dplyr::mutate(y_order = dplyr::row_number())

z_df <- as.data.frame(smooth_mat)

z_df$gene <- rownames(z_df)

z_long <- z_df %>% tidyr::pivot_longer(cols = -gene, names_to = "bin_name", values_to = "z") %>% dplyr::mutate(bin = as.integer(sub("^Bin", 
    "", bin_name))) %>% dplyr::left_join(assignment %>% dplyr::select(gene, module), by = "gene") %>% 
    dplyr::left_join(gene_order_df %>% dplyr::select(gene, y_order), by = "gene") %>% dplyr::left_join(bin_pt, 
    by = "bin")

z_long$module <- factor(z_long$module, levels = paste0("Module ", seq_len(best_k)))

p_heat <- ggplot2::ggplot(z_long, ggplot2::aes(x = pseudotime_median, y = y_order, fill = z)) + ggplot2::geom_raster() +
  ggplot2::facet_grid(rows = ggplot2::vars(module), scales = "free_y", space = "free_y", switch = "y") +
  ggplot2::scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, limits = c(-2.5,
    2.5), oob = scales::squish) + ggplot2::labs(title = "Trajectory-associated gene modules along Monocle3 pseudotime",
  x = "Monocle3 pseudotime", y = NULL, fill = "Row z-score") + ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), strip.placement = "outside",
    strip.background = ggplot2::element_blank())

ggplot2::ggsave(filename = file.path(out_dir, "27_trajectory_module_heatmap_v01.pdf"), plot = p_heat,
  width = 10, height = 11, units = "in", limitsize = FALSE)

profile_df <- z_long %>% dplyr::group_by(module, bin, pseudotime_median) %>% dplyr::summarise(mean_z = mean(z, 
    na.rm = TRUE), se_z = stats::sd(z, na.rm = TRUE)/sqrt(dplyr::n()), .groups = "drop")

write.csv(profile_df, file.path(out_dir, "27_trajectory_module_profiles_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

p_profile <- ggplot2::ggplot(profile_df, ggplot2::aes(x = pseudotime_median, y = mean_z, group = module)) +
  ggplot2::geom_line(linewidth = 0.8) + ggplot2::facet_wrap(~module, ncol = 2, scales = "free_y") +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) + ggplot2::labs(title = "Mean expression profile of each trajectory-gene module",
  x = "Monocle3 pseudotime", y = "Mean row z-score") + ggplot2::theme_classic(base_size = 10)

ggplot2::ggsave(filename = file.path(out_dir, "27_trajectory_module_profiles_v01.pdf"), plot = p_profile,
  width = 9, height = 8, units = "in")

if (requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Mm.eg.db", quietly = TRUE)) {

  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
  })

  background_symbols <- unique(rownames(cds))

  bg_map <- tryCatch(clusterProfiler::bitr(background_symbols, fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Mm.eg.db), error = function(e) NULL)

  if (!is.null(bg_map) && nrow(bg_map) > 0) {

    go_list <- list()

    for (mod in paste0("Module ", seq_len(best_k))) {

      mod_symbols <- assignment$gene[assignment$module == mod]

      mod_map <- tryCatch(clusterProfiler::bitr(mod_symbols, fromType = "SYMBOL", toType = "ENTREZID",
        OrgDb = org.Mm.eg.db), error = function(e) NULL)

      if (is.null(mod_map) || nrow(mod_map) < 10) {
        next
      }

      ego <- tryCatch(clusterProfiler::enrichGO(gene = unique(mod_map$ENTREZID), universe = unique(bg_map$ENTREZID),
        OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
        qvalueCutoff = 0.2, readable = TRUE), error = function(e) NULL)

      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {

        go_df <- as.data.frame(ego)

        go_df$module <- mod

        go_list[[length(go_list) + 1]] <- go_df
      }
    }

    if (length(go_list) > 0) {

      go_all <- dplyr::bind_rows(go_list) %>%
        dplyr::arrange(factor(module, levels = paste0("Module ", seq_len(best_k))), p.adjust)

      write.csv(go_all, file.path(out_dir, "27_trajectory_module_GO_BP_v01.csv"), row.names = FALSE,
        fileEncoding = "UTF-8")

      go_top <- go_all %>%
        dplyr::group_by(module) %>%
        dplyr::slice_min(order_by = p.adjust, n = 8, with_ties = FALSE) %>%
        dplyr::ungroup()

      write.csv(go_top, file.path(out_dir, "27_trajectory_module_GO_BP_top_v01.csv"), row.names = FALSE,
        fileEncoding = "UTF-8")

      p_go <- ggplot2::ggplot(go_top, ggplot2::aes(x = -log10(p.adjust), y = reorder(Description,
        -log10(p.adjust)))) + ggplot2::geom_point(ggplot2::aes(size = Count)) + ggplot2::facet_wrap(~module,
        scales = "free_y", ncol = 2) + ggplot2::labs(title = "GO Biological Process enrichment of trajectory modules",
        x = "-log10 adjusted P", y = NULL, size = "Gene count") + ggplot2::theme_classic(base_size = 9)

      ggplot2::ggsave(filename = file.path(out_dir, "27_trajectory_module_GO_BP_top_v01.pdf"),
        plot = p_go, width = 10, height = 9, units = "in")

    } else {

    }

  } else {

    log_msg("GO BP skipped: background SYMBOL->ENTREZ mapping failed.")
  }

} else {

  log_msg("GO BP skipped: clusterProfiler and/or org.Mm.eg.db not installed.")
}

saveRDS(list(genes = assignment, module_summary = module_summary, k_assessment = k_eval_df, bin_pseudotime = bin_pt,
  z_matrix = smooth_mat), file.path(out_dir, "27_trajectory_gene_modules_3scRNA_v01.rds"), compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)
