source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = 1200)

set.seed(1234)

ROOT <- "D:/zz"

STEP62_DIR <- file.path(ROOT, "results", "62_GSE195783_joint_embedding_batch_assessment_clustering")

INPUT_RDS <- file.path(STEP62_DIR, "13_62_GSE195783_joint_unintegrated_clustered_v01.rds")

OUT_DIR <- file.path(ROOT, "results", "64_GSE195783_spatial_cluster_marker_region_annotation")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v01"

CLUSTER_COL <- "prelim_cluster_r0.4"

MARKER_MIN_PCT <- 0.10
MARKER_LOGFC_THRESHOLD <- 0.25
MARKER_PADJ_THRESHOLD <- 0.05

N_TOP_MARKERS_FOR_CONSISTENCY <- 50L
MIN_CLUSTER_SPOTS_PER_SAMPLE <- 10L
MIN_REST_SPOTS_PER_SAMPLE <- 50L

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "tidyr", "patchwork")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step 64.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

PRESTO_AVAILABLE <- requireNamespace("presto", quietly = TRUE)

write_csv_utf8 <- write_csv_utf8

save_pdf_png <- function(p, stem, width, height) {
    ggsave(filename = file.path(OUT_DIR, paste0(stem, "_", VERSION, ".pdf")), plot = p, width = width, 
        height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
    ggsave(filename = file.path(OUT_DIR, paste0(stem, "_", VERSION, ".png")), plot = p, width = width, 
        height = height, units = "in", dpi = 300, bg = "white")
}

get_spatial_data <- function(obj) {
    ans <- tryCatch({
        SeuratObject::LayerData(object = obj, assay = "Spatial", layer = "data")
    }, error = function(e) {
        NULL
    })
    if (is.null(ans)) {
        ans <- Seurat::GetAssayData(object = obj, assay = "Spatial", slot = "data")
    }
    ans
}

collapse_top_genes <- function(df, n = 10L) {
    if (nrow(df) == 0L) {
        return("")
    }
    paste(head(df$gene, n), collapse = "; ")
}

programs <- list(Neuronal_gray_matter = c("Snap25", "Syt1", "Rbfox3", "Tubb3", "Map2", "Slc17a6", "Slc17a7", 
    "Gad1", "Gad2"), Oligodendrocyte_white_matter = c("Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp", "Olig1", 
    "Olig2", "Sox10"), Astroglial = c("Gfap", "Aqp4", "Aldoc", "Slc1a3", "Glul", "S100b", "Sox9", "Clu"), 
    OPC = c("Pdgfra", "Cspg4", "Olig1", "Olig2", "Sox10", "Gpr17"), Myeloid_immune = c("Aif1", "Tyrobp", 
        "Laptm5", "Ctss", "C1qa", "C1qb", "C1qc", "Lgals3", "Apoe"), Lipid_handling_myeloid = c("Clec7a", 
        "Lpl", "Apoe", "Trem2", "Gpnmb", "Abca1", "Lgals3", "Ctsd", "Ctss"), Fibrotic_ECM = c("Col1a1", 
        "Col1a2", "Col3a1", "Fn1", "Dcn", "Lum", "Pdgfra", "Ccn2"), Vascular_endothelial = c("Pecam1", 
        "Kdr", "Emcn", "Cdh5", "Vwf", "Klf2", "Egfl7"), Pericyte_vSMC = c("Rgs5", "Pdgfrb", "Cspg4", 
        "Des", "Acta2", "Tagln", "Mcam"), Injury_response = c("Atf3", "Jun", "Fos", "Hmox1", "Lcn2", 
        "Spp1", "Gpnmb", "Timp1", "Serpina3n"))

program_display <- c(Neuronal_gray_matter = "Neuronal / gray-matter-like", Oligodendrocyte_white_matter = "Oligodendrocyte / white-matter-like", 
    Astroglial = "Astroglial", OPC = "OPC-like", Myeloid_immune = "Myeloid / immune", Lipid_handling_myeloid = "Lipid-handling myeloid", 
    Fibrotic_ECM = "Fibrotic / ECM", Vascular_endothelial = "Vascular / endothelial", Pericyte_vSMC = "Pericyte / vSMC", 
    Injury_response = "Injury-response")

if (!file.exists(INPUT_RDS)) {
    stop("Step 62 object not found:\n", INPUT_RDS)
}

if (PRESTO_AVAILABLE) {
} else {
}

obj <- readRDS(INPUT_RDS)

if (!inherits(obj, "Seurat")) {
    stop("Step 62 input is not a Seurat object.")
}

if (!"Spatial" %in% names(obj@assays)) {
    stop("Spatial assay is missing.")
}

md <- obj[[]]

required_md <- c(CLUSTER_COL, "GSM", "sample_label", "timepoint", "array_x", "array_y", "histology_registered")

missing_md <- setdiff(required_md, colnames(md))

if (length(missing_md) > 0L) {
    stop("Missing required metadata column(s): ", paste(missing_md, collapse = ", "))
}

if (any(md$histology_registered)) {
    stop("H&E registration must remain FALSE.")
}

if (ncol(obj) != 15204L) {
    stop("Expected 15,204 Step 62 spots, found ", ncol(obj), ".")
}

DefaultAssay(obj) <- "Spatial"

if (inherits(obj[["Spatial"]], "Assay5")) {
    obj <- JoinLayers(object = obj, assay = "Spatial")
}

Idents(obj) <- CLUSTER_COL

cluster_levels <- sort(unique(as.character(obj[[CLUSTER_COL]][, 1])))

md <- obj[[]]

md$cluster <- factor(as.character(md[[CLUSTER_COL]]), levels = cluster_levels)

md$timepoint <- factor(as.character(md$timepoint), levels = TIME_LEVELS)

input_validation <- data.frame(n_spots = ncol(obj), n_features = nrow(obj[["Spatial"]]), n_clusters = length(cluster_levels), 
    n_samples = length(unique(md$sample_label)), n_timepoints = length(unique(md$timepoint)), presto_available = PRESTO_AVAILABLE, 
    stringsAsFactors = FALSE)

write_csv_utf8(input_validation, file.path(OUT_DIR, paste0("00_64_input_validation_", VERSION, ".csv")))

all_markers <- Seurat::FindAllMarkers(object = obj, assay = "Spatial", only.pos = TRUE, test.use = "wilcox", 
    min.pct = MARKER_MIN_PCT, logfc.threshold = MARKER_LOGFC_THRESHOLD, return.thresh = MARKER_PADJ_THRESHOLD, 
    verbose = TRUE)

if (nrow(all_markers) == 0L) {
    stop("FindAllMarkers returned zero markers.")
}

fc_col_candidates <- c("avg_log2FC", "avg_logFC")

fc_col <- intersect(fc_col_candidates, colnames(all_markers))

if (length(fc_col) != 1L) {
    stop("Could not uniquely identify avg_log2FC/avg_logFC column.")
}

all_markers$avg_log2FC_std <- all_markers[[fc_col]]

all_markers <- all_markers %>% arrange(cluster, desc(avg_log2FC_std), p_val_adj)

write_csv_utf8(all_markers, file.path(OUT_DIR, paste0("01_64_all_positive_cluster_markers_", VERSION, 
    ".csv")))

top_marker_candidates <- all_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC_std, n = N_TOP_MARKERS_FOR_CONSISTENCY, 
    with_ties = FALSE) %>% ungroup()

norm_mat <- get_spatial_data(obj)

if (nrow(norm_mat) == 0L || ncol(norm_mat) == 0L) {
    stop("Spatial normalized data matrix is unavailable.")
}

if (!identical(colnames(norm_mat), rownames(md))) {
    idx <- match(rownames(md), colnames(norm_mat))
    if (any(is.na(idx))) {
        stop("Normalized matrix / metadata alignment failed.")
    }
    norm_mat <- norm_mat[, idx, drop = FALSE]
}

consistency_rows <- list()
row_counter <- 1L

for (cl in cluster_levels) {
    candidate_genes <- top_marker_candidates$gene[as.character(top_marker_candidates$cluster) == cl]
    candidate_genes <- intersect(candidate_genes, rownames(norm_mat))
    if (length(candidate_genes) == 0L) {
        next
    }
    for (samp in unique(as.character(md$sample_label))) {
        cl_cells <- which(as.character(md$cluster) == cl & as.character(md$sample_label) == samp)
        rest_cells <- which(as.character(md$cluster) != cl & as.character(md$sample_label) == samp)
        if (length(cl_cells) < MIN_CLUSTER_SPOTS_PER_SAMPLE || length(rest_cells) < MIN_REST_SPOTS_PER_SAMPLE) {
            next
        }
        cl_mean <- Matrix::rowMeans(norm_mat[candidate_genes, cl_cells, drop = FALSE])
        rest_mean <- Matrix::rowMeans(norm_mat[candidate_genes, rest_cells, drop = FALSE])
        delta <- cl_mean - rest_mean
        consistency_rows[[row_counter]] <- data.frame(cluster = cl, sample_label = samp, gene = candidate_genes, 
            mean_in_cluster = as.numeric(cl_mean), mean_in_rest = as.numeric(rest_mean), delta_mean_log_normalized = as.numeric(delta), 
            n_cluster_spots = length(cl_cells), n_rest_spots = length(rest_cells), stringsAsFactors = FALSE)
        row_counter <- row_counter + 1L
    }
}

if (length(consistency_rows) == 0L) {
    stop("No sample-aware marker comparisons were evaluable.")
}

consistency_long <- bind_rows(consistency_rows)

consistency_summary <- consistency_long %>% group_by(cluster, gene) %>% summarise(n_samples_evaluable = n(), 
    median_delta = median(delta_mean_log_normalized), min_delta = min(delta_mean_log_normalized), max_delta = max(delta_mean_log_normalized), 
    fraction_samples_positive = mean(delta_mean_log_normalized > 0), .groups = "drop")

marker_summary <- top_marker_candidates %>% select(cluster, gene, p_val, p_val_adj, pct.1, pct.2, avg_log2FC_std) %>% 
    left_join(consistency_summary, by = c("cluster", "gene")) %>% mutate(robust_directional_support = !is.na(n_samples_evaluable) & 
    n_samples_evaluable >= 2 & fraction_samples_positive >= 0.75 & median_delta > 0) %>% arrange(cluster, 
    desc(robust_directional_support), desc(avg_log2FC_std))

write_csv_utf8(marker_summary, file.path(OUT_DIR, paste0("02_64_top_markers_sample_consistency_", VERSION, 
    ".csv")))

cluster_time_counts <- as.data.frame(table(cluster = md$cluster, timepoint = md$timepoint), stringsAsFactors = FALSE)

cluster_time_counts <- cluster_time_counts[cluster_time_counts$Freq > 0, , drop = FALSE]

cluster_time_summary <- cluster_time_counts %>% group_by(cluster) %>% mutate(fraction_within_cluster = Freq/sum(Freq)) %>% 
    ungroup()

write_csv_utf8(cluster_time_summary, file.path(OUT_DIR, paste0("03_64_cluster_timepoint_composition_", 
    VERSION, ".csv")))

all_program_genes <- unique(unlist(programs, use.names = FALSE))

present_program_genes <- intersect(all_program_genes, rownames(norm_mat))

if (length(present_program_genes) < 20L) {
    stop("Too few predefined program genes are present in Spatial expression matrix.")
}

cluster_mean <- sapply(cluster_levels, function(cl) {
    cells <- which(as.character(md$cluster) == cl)
    Matrix::rowMeans(norm_mat[present_program_genes, cells, drop = FALSE])
})

if (is.null(dim(cluster_mean))) {
    stop("Cluster-average expression matrix could not be constructed.")
}

rownames(cluster_mean) <- present_program_genes

colnames(cluster_mean) <- cluster_levels

gene_z <- t(apply(cluster_mean, 1, function(x) {
    sx <- stats::sd(x)
    if (!is.finite(sx) || sx == 0) {
        return(rep(0, length(x)))
    }
    as.numeric(scale(x))
}))

rownames(gene_z) <- rownames(cluster_mean)

colnames(gene_z) <- colnames(cluster_mean)

program_score_rows <- list()
pcounter <- 1L

for (prog in names(programs)) {
    genes <- intersect(programs[[prog]], rownames(gene_z))
    if (length(genes) == 0L) {
        next
    }
    score <- colMeans(gene_z[genes, , drop = FALSE])
    program_score_rows[[pcounter]] <- data.frame(cluster = names(score), program = prog, program_display = unname(program_display[prog]), 
        score_z = as.numeric(score), n_genes_present = length(genes), genes_present = paste(genes, collapse = "; "), 
        stringsAsFactors = FALSE)
    pcounter <- pcounter + 1L
}

program_scores <- bind_rows(program_score_rows)

write_csv_utf8(program_scores, file.path(OUT_DIR, paste0("04_64_cluster_program_scores_", VERSION, ".csv")))

annotation_rows <- list()

for (cl in cluster_levels) {
    ps <- program_scores[program_scores$cluster == cl, , drop = FALSE]
    ps <- ps[order(ps$score_z, decreasing = TRUE), , drop = FALSE]
    top1 <- ps[1, , drop = FALSE]
    top2 <- ps[2, , drop = FALSE]
    mk <- marker_summary[as.character(marker_summary$cluster) == cl, , drop = FALSE]
    mk_robust <- mk[mk$robust_directional_support %in% TRUE, , drop = FALSE]
    mk_robust <- mk_robust[order(mk_robust$avg_log2FC_std, decreasing = TRUE), , drop = FALSE]
    mk_all <- mk[order(mk$avg_log2FC_std, decreasing = TRUE), , drop = FALSE]
    margin <- top1$score_z - top2$score_z
    annotation_status <- if (top1$score_z >= 0.5 && margin >= 0.25) {
        "provisional_clearer"
    }
    else {
        "provisional_mixed_or_ambiguous"
    }
    time_sub <- cluster_time_summary[as.character(cluster_time_summary$cluster) == cl, , drop = FALSE]
    time_sub <- time_sub[order(time_sub$fraction_within_cluster, decreasing = TRUE), , drop = FALSE]
    annotation_rows[[cl]] <- data.frame(cluster = cl, provisional_program = top1$program_display, top_program_score = top1$score_z, 
        second_program = top2$program_display, second_program_score = top2$score_z, score_margin = margin, 
        annotation_status = annotation_status, robust_top_markers = collapse_top_genes(mk_robust, n = 10L), 
        top_markers_overall = collapse_top_genes(mk_all, n = 10L), dominant_timepoint = as.character(time_sub$timepoint[1]), 
        dominant_timepoint_fraction = time_sub$fraction_within_cluster[1], stringsAsFactors = FALSE)
}

annotation_table <- bind_rows(annotation_rows)

annotation_table$cluster <- factor(annotation_table$cluster, levels = cluster_levels)

annotation_table <- annotation_table %>% arrange(cluster)

write_csv_utf8(annotation_table, file.path(OUT_DIR, paste0("05_64_provisional_transcriptional_region_annotation_", 
    VERSION, ".csv")))

program_scores$cluster <- factor(program_scores$cluster, levels = rev(cluster_levels))

program_scores$program_display <- factor(program_scores$program_display, levels = unname(program_display))

p_program_heatmap <- ggplot(program_scores, aes(x = program_display, y = cluster, fill = score_z)) + 
    geom_tile(linewidth = 0.25) + scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", 
    midpoint = 0) + labs(title = "Transcriptional program scores across spatial clusters", subtitle = "Cluster-average gene-wise z scores; descriptive annotation only", 
    x = NULL, y = "Spatial cluster", fill = "Program\nscore (z)") + theme_classic(base_family = "Arial", 
    base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), plot.title = element_text(size = 12), 
    plot.subtitle = element_text(size = 8.5))

save_pdf_png(p = p_program_heatmap, stem = "06_64_transcriptional_program_heatmap", width = 9.2, height = 5.8)

dot_marker_genes <- marker_summary %>% group_by(cluster) %>% arrange(desc(robust_directional_support), 
    desc(avg_log2FC_std), .by_group = TRUE) %>% slice_head(n = 3) %>% ungroup() %>% pull(gene) %>% unique()

dot_marker_genes <- intersect(dot_marker_genes, rownames(obj))

if (length(dot_marker_genes) > 0L) {
    Idents(obj) <- CLUSTER_COL
    p_marker_dot <- DotPlot(object = obj, assay = "Spatial", features = dot_marker_genes, cols = c("grey90", 
        "#B40426"), dot.scale = 5) + RotatedAxis() + labs(title = "Representative spatial-cluster markers", 
        x = NULL, y = "Spatial cluster") + theme_classic(base_family = "Arial", base_size = 8) + theme(axis.text.x = element_text(size = 7), 
        plot.title = element_text(size = 11))
    save_pdf_png(p = p_marker_dot, stem = "07_64_representative_marker_dotplot", width = 12.5, height = 5.8)
}

map_df <- obj[[]]

map_df$cluster <- factor(as.character(map_df[[CLUSTER_COL]]), levels = cluster_levels)

sample_order <- unique(as.character(md$sample_label))

map_df$sample_label <- factor(as.character(map_df$sample_label), levels = sample_order)

p_spatial_clusters <- ggplot(map_df, aes(x = array_x, y = array_y, color = cluster)) + geom_point(size = 0.48, 
    alpha = 0.95) + coord_equal() + scale_y_reverse() + facet_wrap(~sample_label, ncol = 4) + labs(title = "Step 62 spatial clusters across GSE195783 samples", 
    subtitle = "Visium V1 array coordinates only; no H&E/pixel registration", x = "Visium array X", y = "Visium array Y", 
    color = "Cluster") + theme_classic(base_family = "Arial", base_size = 9) + theme(plot.title = element_text(size = 12), 
    plot.subtitle = element_text(size = 8.5), strip.text = element_text(size = 8), legend.position = "right")

save_pdf_png(p = p_spatial_clusters, stem = "08_64_spatial_cluster_maps", width = 11.5, height = 7.4)

n_markers_per_cluster <- table(all_markers$cluster)

n_robust_per_cluster <- marker_summary %>% group_by(cluster) %>% summarise(n_robust = sum(robust_directional_support, 
    na.rm = TRUE), .groups = "drop")

final_checks <- data.frame(check = c("Primary input is Step62 unintegrated object", "Input contains 15204 retained spots", 
    "Working clustering contains 11 clusters", "Positive markers were detected for every cluster", "Sample-aware marker consistency was evaluable", 
    "All 10 predefined transcriptional programs were scored", "Provisional annotation table contains 11 clusters", 
    "Spatial cluster map was generated from array coordinates", "No H&E registration is claimed", "No scRNA-to-spatial deconvolution was performed"), 
    pass = c(grepl("62_GSE195783", INPUT_RDS), ncol(obj) == 15204L, length(cluster_levels) == 11L, length(unique(all_markers$cluster)) == 
        length(cluster_levels), nrow(consistency_summary) > 0L, length(unique(program_scores$program)) == 
        length(programs), nrow(annotation_table) == 11L, file.exists(file.path(OUT_DIR, paste0("08_64_spatial_cluster_maps_", 
        VERSION, ".pdf"))), all(obj$histology_registered == FALSE), TRUE), stringsAsFactors = FALSE)

write_csv_utf8(final_checks, file.path(OUT_DIR, paste0("09_64_GSE195783_FINAL_QC_", VERSION, ".csv")))

annotation_lines <- apply(annotation_table, 1, function(z) {
    paste0("  Cluster ", z[["cluster"]], ": ", z[["provisional_program"]], " [", z[["annotation_status"]], 
        "]", "; top markers = ", z[["top_markers_overall"]])
})

robust_lines <- apply(n_robust_per_cluster, 1, function(z) {
    paste0("  Cluster ", z[["cluster"]], ": ", z[["n_robust"]], " robust-directional top markers")
})

report_file <- file.path(OUT_DIR, paste0("10_64_GSE195783_cluster_marker_region_annotation_REPORT_", 
    VERSION, ".txt"))

report_con <- file(report_file, open = "wt", encoding = "UTF-8")

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = report_con)

close(report_con)

print(annotation_table)
print(final_checks)
