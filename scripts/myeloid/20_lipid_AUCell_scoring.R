source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
N_CORES <- 1

UMAP_REDUCTION <- "umap_myeloid_harmony30"

STATE_LEVELS <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

MICROGLIA_STATES <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia")

MAIN_TIMES <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")

FOCUS_TIMES <- c("Control", "7 dpi", "14 dpi", "28 dpi")

STATE_COLORS <- c(`Homeostatic Microglia` = "#3B6FB6", `Clec7a+ lipid-associated Microglia` = "#E08B2C", 
    `Activated inflammatory Microglia` = "#C94C4C", `Ccl12+ acute Microglia` = "#D66BA0", `Tnf+ inflammatory Microglia` = "#8E63B6", 
    `Proliferating Microglia` = "#2F9E8F", `MHC-II+ antigen-presenting Macrophage` = "#4F8A3B", `Gpnmb+ tissue-repair Macrophage` = "#8A6D5A", 
    `Acute-response Macrophage` = "#D9A514")

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "AUCell", "ggplot2", "patchwork", "dplyr", "tidyr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
    stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), "\nInstall AUCell with:\n", 
        "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n", "BiocManager::install('AUCell')")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(AUCell)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "19_myeloid_state_robustness_temporal_dynamics", 
    "19_myeloid_state_robustness_temporal_dynamics_3scRNA_v02.rds")

if (file.exists(expected_input)) {
    in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
} else {
    hits <- list.files(project_dir, pattern = "^19_myeloid_state_robustness_temporal_dynamics_3scRNA_v02\\.rds$", 
        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        stop("Input file not found:\n", "19_myeloid_state_robustness_temporal_dynamics_3scRNA_v02.rds")
    }
    in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "20_AUCell_lipid_metabolism_scoring")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "20_AUCell_lipid_metabolism_scoring_v03.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

required_meta <- c("myeloid_state", "source_dataset", "timepoint_main")

missing_meta <- setdiff(required_meta, colnames(obj[[]]))

if (length(missing_meta) > 0) {
    stop("Missing required metadata field(s): ", paste(missing_meta, collapse = ", "))
}

if (!UMAP_REDUCTION %in% names(obj@reductions)) {
    stop("Missing UMAP reduction: ", UMAP_REDUCTION)
}

obj$myeloid_state <- factor(as.character(obj$myeloid_state), levels = STATE_LEVELS)

obj$timepoint_main <- as.character(obj$timepoint_main)

obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) {
    log_msg("JoinLayers skipped/failed safely: ", conditionMessage(e))
    obj
})

DefaultAssay(obj) <- "RNA"

rna_layers <- tryCatch(Layers(obj[["RNA"]]), error = function(e) character())

if (!"data" %in% rna_layers) {
    obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, 
        verbose = TRUE)
}

lipid_gene_sets <- list(FA_uptake_transport = c("Cd36", "Slc27a1", "Slc27a2", "Slc27a3", "Slc27a4", "Slc27a5", 
    "Fabp3", "Fabp4", "Fabp5", "Lpl", "Lrp1", "Vldlr"), FA_synthesis = c("Acly", "Acaca", "Fasn", "Scd1", 
    "Scd2", "Elovl5", "Elovl6", "Acss2", "Me1", "Srebf1"), FA_beta_oxidation = c("Cpt1a", "Cpt2", "Acadm", 
    "Acadl", "Acads", "Hadha", "Hadhb", "Echs1", "Acox1", "Acaa1a", "Acaa2", "Eci1"), Cholesterol_biosynthesis = c("Hmgcr", 
    "Hmgcs1", "Mvk", "Pmvk", "Mvd", "Fdps", "Fdft1", "Sqle", "Lss", "Cyp51", "Msmo1", "Sc5d", "Dhcr7", 
    "Dhcr24", "Nsdhl", "Srebf2"), Cholesterol_efflux_transport = c("Abca1", "Abcg1", "Apoe", "Pltp", 
    "Npc1", "Npc2", "Lipa", "Soat1", "Lrp1", "Osbpl1a", "Nr1h3", "Nr1h2"), Lipoprotein_processing = c("Apoe", 
    "Apoc1", "Apoc2", "Apoc4", "Lpl", "Lrp1", "Vldlr", "Ldlr", "Abca1", "Abcg1", "Pltp", "Lipa", "Soat1"), 
    Lipid_droplet_storage = c("Plin2", "Plin3", "Dgat1", "Dgat2", "Acsl1", "Acsl3", "Acsl4", "Soat1", 
        "Pnpla2", "Abhd5", "Fitm2", "G0s2"), Phospholipid_metabolism = c("Chka", "Chkb", "Pcyt1a", "Pcyt1b", 
        "Chpt1", "Cept1", "Pemt", "Pla2g4a", "Pla2g6", "Lpcat1", "Lpcat2", "Lpcat3", "Mboat7", "Ptdss1", 
        "Ptdss2"), Sphingolipid_metabolism = c("Sptlc1", "Sptlc2", "Sptssa", "Kdsr", "Cers2", "Cers5", 
        "Cers6", "Ugcg", "Sgms1", "Sgms2", "Smpd1", "Smpd2", "Asah1", "Sphk1", "Sphk2"), Arachidonic_eicosanoid = c("Pla2g4a", 
        "Alox5", "Alox5ap", "Lta4h", "Ltc4s", "Ptgs1", "Ptgs2", "Hpgds", "Ptgds", "Ptges", "Ptges2", 
        "Ptges3"), Lysosomal_lipid_catabolism = c("Lipa", "Npc1", "Npc2", "Gba1", "Galc", "Hexa", "Hexb", 
        "Smpd1", "Asah1", "Ctsb", "Ctsd", "Lamp1", "Lamp2", "Lgmn"))

expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

if (!inherits(expr_mat, "dgCMatrix")) {
    expr_mat <- as(expr_mat, "dgCMatrix")
}

nonzero_genes <- Matrix::rowSums(expr_mat > 0) > 0

expr_mat <- expr_mat[nonzero_genes, , drop = FALSE]

gene_manifest <- bind_rows(lapply(names(lipid_gene_sets), function(gs) {
    requested <- unique(lipid_gene_sets[[gs]])
    present <- requested[requested %in% rownames(expr_mat)]
    missing <- setdiff(requested, present)
    data.frame(gene_set = gs, gene = requested, present_in_object = requested %in% present, stringsAsFactors = FALSE)
}))

write.csv(gene_manifest, file.path(out_dir, "20_AUCell_lipid_gene_set_manifest_v03.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

lipid_gene_sets_present <- lapply(lipid_gene_sets, function(gs) {
    intersect(gs, rownames(expr_mat))
})

gene_set_sizes <- data.frame(gene_set = names(lipid_gene_sets_present), n_requested = vapply(lipid_gene_sets, 
    length, FUN.VALUE = integer(1)), n_present = vapply(lipid_gene_sets_present, length, FUN.VALUE = integer(1)), 
    stringsAsFactors = FALSE)

capture.output(print(gene_set_sizes), file = log_file, append = TRUE)

too_small <- gene_set_sizes$gene_set[gene_set_sizes$n_present < 5]

if (length(too_small) > 0) {
    stop("These AUCell modules have fewer than 5 detected genes: ", paste(too_small, collapse = ", "), 
        "\nInspect 20_AUCell_lipid_gene_set_manifest_v03.csv.")
}

set.seed(SEED)

ranking_cache <- file.path(out_dir, "20_AUCell_rankings_v01.rds")

reuse_ranking <- FALSE

if (file.exists(ranking_cache)) {
    cached_rankings <- tryCatch(readRDS(ranking_cache), error = function(e) NULL)
    if (!is.null(cached_rankings)) {
        cached_cells <- colnames(cached_rankings)
        if (length(cached_cells) == ncol(expr_mat) && setequal(cached_cells, colnames(expr_mat))) {
            cells_rankings <- cached_rankings
            if (!identical(colnames(cells_rankings), colnames(expr_mat))) {
                cells_rankings <- cells_rankings[, colnames(expr_mat), drop = FALSE]
            }
            reuse_ranking <- TRUE
        }
    }
}

if (!reuse_ranking) {
    cells_rankings <- AUCell_buildRankings(exprMat = expr_mat, plotStats = FALSE, splitByBlocks = TRUE, 
        verbose = TRUE)
    saveRDS(cells_rankings, ranking_cache, compress = FALSE)
}

rank_n_genes <- nrow(cells_rankings)

auc_max_rank <- ceiling(0.05 * rank_n_genes)

cells_auc <- AUCell_calcAUC(geneSets = lipid_gene_sets_present, rankings = cells_rankings, normAUC = TRUE, 
    aucMaxRank = auc_max_rank, verbose = TRUE)

auc_mat <- as.matrix(getAUC(cells_auc))

if (!identical(colnames(auc_mat), Cells(obj))) {
    auc_mat <- auc_mat[, Cells(obj), drop = FALSE]
}

score_names <- paste0("AUC_", rownames(auc_mat))

rownames(auc_mat) <- score_names

auc_df <- as.data.frame(t(auc_mat), check.names = FALSE)

auc_z <- scale(as.matrix(auc_df))

auc_z[!is.finite(auc_z)] <- 0

auc_df$AUC_LipidCompositeZ <- rowMeans(auc_z)

rm(auc_z)

gc(verbose = FALSE)

obj <- AddMetaData(obj, metadata = auc_df)

cell_score_export <- data.frame(cell = rownames(auc_df), myeloid_state = as.character(obj$myeloid_state), 
    source_dataset = as.character(obj$source_dataset), timepoint_main = as.character(obj$timepoint_main), 
    auc_df, check.names = FALSE, stringsAsFactors = FALSE)

write.csv(cell_score_export, file.path(out_dir, "20_AUCell_cell_scores_v03.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

score_cols <- c(score_names, "AUC_LipidCompositeZ")

score_long <- cell_score_export %>% select(cell, myeloid_state, source_dataset, timepoint_main, all_of(score_cols)) %>% 
    pivot_longer(cols = all_of(score_cols), names_to = "score", values_to = "value")

state_summary <- score_long %>% group_by(myeloid_state, score) %>% summarise(n_cells = n(), mean = mean(value, 
    na.rm = TRUE), median = median(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), q25 = quantile(value, 
    0.25, na.rm = TRUE), q75 = quantile(value, 0.75, na.rm = TRUE), .groups = "drop") %>% mutate(myeloid_state = factor(myeloid_state, 
    levels = STATE_LEVELS)) %>% arrange(myeloid_state, score)

write.csv(state_summary, file.path(out_dir, "20_AUCell_state_summary_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

time_summary_main <- score_long %>% filter(timepoint_main %in% MAIN_TIMES) %>% mutate(timepoint_main = factor(timepoint_main, 
    levels = MAIN_TIMES)) %>% group_by(timepoint_main, score) %>% summarise(n_cells = n(), mean = mean(value, 
    na.rm = TRUE), median = median(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), q25 = quantile(value, 
    0.25, na.rm = TRUE), q75 = quantile(value, 0.75, na.rm = TRUE), .groups = "drop")

write.csv(time_summary_main, file.path(out_dir, "20_AUCell_time_summary_main_v03.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

state_time_summary_main <- score_long %>% filter(timepoint_main %in% MAIN_TIMES) %>% mutate(myeloid_state = factor(myeloid_state, 
    levels = STATE_LEVELS), timepoint_main = factor(timepoint_main, levels = MAIN_TIMES)) %>% group_by(myeloid_state, 
    timepoint_main, score) %>% summarise(n_cells = n(), mean = mean(value, na.rm = TRUE), median = median(value, 
    na.rm = TRUE), sd = sd(value, na.rm = TRUE), q25 = quantile(value, 0.25, na.rm = TRUE), q75 = quantile(value, 
    0.75, na.rm = TRUE), .groups = "drop")

write.csv(state_time_summary_main, file.path(out_dir, "20_AUCell_state_time_summary_main_v03.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

clec7a_focus <- state_time_summary_main %>% filter(as.character(myeloid_state) == "Clec7a+ lipid-associated Microglia", 
    as.character(timepoint_main) %in% FOCUS_TIMES) %>% mutate(timepoint_main = factor(as.character(timepoint_main), 
    levels = FOCUS_TIMES)) %>% arrange(score, timepoint_main)

write.csv(clec7a_focus, file.path(out_dir, "20_AUCell_Clec7a_focus_time_summary_v03.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

zscore_group_matrix <- function(summary_df, group_col, score_keep = score_names) {
    tmp <- summary_df %>% filter(score %in% score_keep) %>% select(all_of(group_col), score, mean) %>% 
        pivot_wider(names_from = all_of(group_col), values_from = mean)
    rn <- tmp$score
    mat <- as.matrix(tmp[, setdiff(colnames(tmp), "score"), drop = FALSE])
    rownames(mat) <- rn
    z <- t(scale(t(mat)))
    z[!is.finite(z)] <- 0
    z <- pmax(pmin(z, 2.5), -2.5)
    as.data.frame(as.table(z), stringsAsFactors = FALSE) %>% setNames(c("score", "group", "z"))
}

state_heat_df <- zscore_group_matrix(state_summary, "myeloid_state")

state_heat_df$group <- factor(state_heat_df$group, levels = STATE_LEVELS)

time_heat_df <- zscore_group_matrix(time_summary_main, "timepoint_main")

time_heat_df$group <- factor(time_heat_df$group, levels = MAIN_TIMES)

clec7a_heat_df <- zscore_group_matrix(clec7a_focus, "timepoint_main")

clec7a_heat_df$group <- factor(clec7a_heat_df$group, levels = FOCUS_TIMES)

diagnostic_pdf <- file.path(out_dir, "20_AUCell_lipid_metabolism_diagnostics_v03.pdf")

pdf(diagnostic_pdf, width = 12, height = 8.5, onefile = TRUE)

p_umap <- DimPlot(obj, reduction = UMAP_REDUCTION, group.by = "myeloid_state", cols = STATE_COLORS, pt.size = 0.1, 
    raster = TRUE) + ggtitle("Myeloid states used for AUCell lipid-metabolism scoring") + theme_classic(base_size = 10)

print(p_umap)

p_state_heat <- ggplot(state_heat_df, aes(x = group, y = score, fill = z)) + geom_tile(color = "white", 
    linewidth = 0.25) + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D94F4F", midpoint = 0, 
    limits = c(-2.5, 2.5), name = "Row Z") + labs(title = "AUCell lipid-module activity across myeloid states", 
    x = NULL, y = NULL) + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 35, 
    hjust = 1), axis.ticks = element_blank())

print(p_state_heat)

p_time_heat <- ggplot(time_heat_df, aes(x = group, y = score, fill = z)) + geom_tile(color = "white", 
    linewidth = 0.25) + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D94F4F", midpoint = 0, 
    limits = c(-2.5, 2.5), name = "Row Z") + labs(title = "AUCell lipid-module activity across the main SCI timeline", 
    x = NULL, y = NULL) + theme_classic(base_size = 9) + theme(axis.ticks = element_blank())

print(p_time_heat)

p_clec_time <- ggplot(clec7a_heat_df, aes(x = group, y = score, fill = z)) + geom_tile(color = "white", 
    linewidth = 0.25) + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D94F4F", midpoint = 0, 
    limits = c(-2.5, 2.5), name = "Row Z") + labs(title = "Clec7a+ candidate state: lipid-module dynamics", 
    subtitle = "Control -> 7 -> 14 -> 28 dpi", x = NULL, y = NULL) + theme_classic(base_size = 9) + theme(axis.ticks = element_blank())

print(p_clec_time)

composite_plot_df <- cell_score_export %>% mutate(myeloid_state = factor(myeloid_state, levels = STATE_LEVELS))

set.seed(SEED)

if (
  nrow(
    composite_plot_df
  ) > 30000
) {

  set.seed(SEED)

  composite_plot_df <- composite_plot_df %>%
    mutate(
      .rand_downsample = runif(n())
    ) %>%
    group_by(
      myeloid_state
    ) %>%
    arrange(
      .rand_downsample,
      .by_group = TRUE
    ) %>%
    filter(
      row_number() <= 3500
    ) %>%
    ungroup() %>%
    select(
      -.rand_downsample
    )
}

p_composite_state <- ggplot(composite_plot_df, aes(x = myeloid_state, y = AUC_LipidCompositeZ, fill = myeloid_state)) + 
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) + geom_boxplot(width = 0.12, outlier.shape = NA, 
    linewidth = 0.25) + scale_fill_manual(values = STATE_COLORS, drop = FALSE) + labs(title = "Composite lipid-program activity across myeloid states", 
    x = NULL, y = "Composite AUCell score (mean module Z)") + theme_classic(base_size = 9) + theme(legend.position = "none", 
    axis.text.x = element_text(angle = 35, hjust = 1))

print(p_composite_state)

focus_plot_df <- cell_score_export %>% filter(myeloid_state == "Clec7a+ lipid-associated Microglia", 
    timepoint_main %in% FOCUS_TIMES) %>% mutate(timepoint_main = factor(timepoint_main, levels = FOCUS_TIMES))

set.seed(SEED)

if (nrow(focus_plot_df) > 16000) {
    set.seed(SEED)
    focus_plot_df <- focus_plot_df %>% mutate(.rand_downsample = runif(n())) %>% group_by(timepoint_main) %>% 
        arrange(.rand_downsample, .by_group = TRUE) %>% filter(row_number() <= 3500) %>% ungroup() %>% 
        select(-.rand_downsample)
}

p_focus_time <- ggplot(focus_plot_df, aes(x = timepoint_main, y = AUC_LipidCompositeZ)) + geom_violin(scale = "width", 
    trim = TRUE, linewidth = 0.2) + geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25) + 
    labs(title = "Clec7a+ candidate state: composite lipid activity over time", x = NULL, y = "Composite AUCell score (mean module Z)") + 
    theme_classic(base_size = 10)

print(p_focus_time)

feature_scores <- c("AUC_FA_uptake_transport", "AUC_Cholesterol_efflux_transport", "AUC_Lipid_droplet_storage", 
    "AUC_Lysosomal_lipid_catabolism", "AUC_LipidCompositeZ")

feature_scores <- feature_scores[feature_scores %in% colnames(obj[[]])]

p_feature <- FeaturePlot(obj, features = feature_scores, reduction = UMAP_REDUCTION, ncol = 3, order = TRUE, 
    pt.size = 0.1, raster = TRUE) & theme_classic(base_size = 8)

print(p_feature)

dev.off()

out_rds <- file.path(out_dir, "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

saveRDS(obj, out_rds, compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)
