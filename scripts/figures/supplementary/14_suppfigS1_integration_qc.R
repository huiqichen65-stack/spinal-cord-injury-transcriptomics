options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrastr)
})

FIG_WIDTH_IN <- 7
FIG_HEIGHT_IN <- 8.8

PANEL_LETTER_SIZE <- 10
PANEL_TITLE_SIZE <- 8.5
AXIS_TITLE_SIZE <- 7.5
AXIS_TEXT_SIZE <- 6.5
LEGEND_TEXT_SIZE <- 6.5

FONT_FAMILY <- "sans"
USE_SHOWTEXT <- FALSE

arial_regular <- "C:/Windows/Fonts/arial.ttf"
arial_bold <- "C:/Windows/Fonts/arialbd.ttf"

if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE) && file.exists(arial_regular)) {
    if (file.exists(arial_bold)) {
        sysfonts::font_add(family = "ArialLocal", regular = arial_regular, bold = arial_bold)
    }
    else {
        sysfonts::font_add(family = "ArialLocal", regular = arial_regular)
    }
    showtext::showtext_auto(enable = TRUE)
    FONT_FAMILY <- "ArialLocal"
    USE_SHOWTEXT <- TRUE
}

project_dir <- normalizePath(if (dir.exists("D:/zz")) "D:/zz" else getwd(), winslash = "/", mustWork = TRUE)

final_file <- file.path(project_dir, "results", "09_Figure1A_major_celltype_UMAP", "09_Figure1A_major_celltype_UMAP_object_v04.rds")

harmony06_file <- file.path(project_dir, "results", "06_Harmony_cluster_resolution", "06_Harmony_cluster_resolution_3scRNA_v01.rds")

decontx_files <- c(`17702045` = file.path(project_dir, "results", "03_DecontX", "03_DecontX_17702045_v01.rds"),
  GSE162610 = file.path(project_dir, "results", "03_DecontX", "03_DecontX_GSE162610_v01.rds"), GSE196928 = file.path(project_dir,
    "results", "03_DecontX", "03_DecontX_GSE196928_Control_v01.rds"))

out_dir <- file.path(project_dir, "results", "14_SuppFigureS1_QC_integration")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

doublet_candidates <- list.files(file.path(project_dir, "results"), pattern = "^04_doublet_3scRNA_summary_v03\\.csv$",
  recursive = TRUE, full.names = TRUE)

if (length(doublet_candidates) == 0) {
  stop("未找到 04_doublet_3scRNA_summary_v03.csv。\n", "请确认第04步doublet结果仍在 results 目录下。")
}

doublet_csv <- doublet_candidates[1]

pre_cache_rds <- file.path(out_dir, "14_SuppFigureS1_preHarmony_UMAP_coordinates_v06.rds")

pre_cache_csv <- file.path(out_dir, "14_SuppFigureS1_preHarmony_UMAP_coordinates_v06.csv")

required_files <- c(final_file, harmony06_file, decontx_files)

if (any(!file.exists(required_files))) {
  stop("以下必要文件不存在：\n", paste(required_files[!file.exists(required_files)], collapse = "\n"))
}

dataset_levels <- c("17702045", "GSE162610", "GSE196928")

dataset_cols <- c(`17702045` = "#4C78A8", GSE162610 = "#F58518", GSE196928 = "#54A24B")

standardize_dataset <- function(x) {
    x <- as.character(x)
    out <- rep(NA_character_, length(x))
    out[grepl("17702045", x, ignore.case = TRUE)] <- "17702045"
    out[grepl("GSE162610", x, ignore.case = TRUE)] <- "GSE162610"
    out[grepl("GSE196928", x, ignore.case = TRUE)] <- "GSE196928"
    out
}

theme_pub <- theme_classic(base_family = FONT_FAMILY, base_size = AXIS_TEXT_SIZE) + theme(axis.title = element_text(family = FONT_FAMILY,
  size = AXIS_TITLE_SIZE, color = "black"), axis.text = element_text(family = FONT_FAMILY, size = AXIS_TEXT_SIZE,
  color = "black"), axis.line = element_line(linewidth = 0.35, color = "black"), axis.ticks = element_line(linewidth = 0.3,
  color = "black"), legend.title = element_blank(), legend.text = element_text(family = FONT_FAMILY,
  size = LEGEND_TEXT_SIZE, color = "black"), plot.title = element_text(family = FONT_FAMILY, size = PANEL_TITLE_SIZE,
  face = "bold", color = "black", hjust = 0), plot.tag = element_text(family = FONT_FAMILY, size = PANEL_LETTER_SIZE,
  face = "bold", color = "black"), plot.tag.position = c(0, 1), plot.margin = margin(4, 4, 4, 4, unit = "pt"))

obj <- readRDS(final_file)

if (!"source_dataset" %in% colnames(obj[[]])) {
  stop("final atlas中缺少 source_dataset。")
}

if (!"umap_harmony30" %in% names(obj@reductions)) {
  stop("final atlas中缺少 umap_harmony30。")
}

src_std <- standardize_dataset(obj$source_dataset)

if (anyNA(src_std)) {
  stop("final atlas中存在无法识别的 source_dataset。")
}

final_cells <- Cells(obj)

source_map <- src_std
names(source_map) <- final_cells

qc_md <- data.frame(cell = final_cells, source_dataset = factor(src_std, levels = dataset_levels), nFeature_RNA = obj$nFeature_RNA,
  nCount_RNA = obj$nCount_RNA, stringsAsFactors = FALSE)

qc_summary <- qc_md %>% group_by(source_dataset) %>% summarise(n_cells = n(), median_nFeature_RNA = median(nFeature_RNA, 
    na.rm = TRUE), median_nCount_RNA = median(nCount_RNA, na.rm = TRUE), .groups = "drop")

write.csv(qc_summary, file.path(out_dir, "14_SuppFigureS1_QC_summary_v07.csv"), row.names = FALSE, fileEncoding = "UTF-8")

set.seed(1234)

qc_plot_df <- qc_md %>% group_by(source_dataset) %>% group_modify(~slice_sample(.x, n = min(5000, nrow(.x)), 
    replace = FALSE)) %>% ungroup()

post_mat <- Embeddings(obj, "umap_harmony30")[, 1:2, drop = FALSE]

post_df_full <- data.frame(cell = rownames(post_mat), UMAP_1 = post_mat[, 1], UMAP_2 = post_mat[, 2],
  source_dataset = factor(source_map[rownames(post_mat)], levels = dataset_levels), stringsAsFactors = FALSE)

rm(obj, post_mat)
gc(verbose = FALSE)

pA1 <- ggplot(qc_plot_df, aes(x = source_dataset, y = nFeature_RNA, fill = source_dataset)) + geom_violin(scale = "width",
  trim = TRUE, linewidth = 0.25) + geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.25,
  fill = "white") + scale_fill_manual(values = dataset_cols) + labs(tag = "A", title = "Quality-control metrics",
  x = NULL, y = "Detected genes") + theme_pub + theme(legend.position = "none", axis.text.x = element_text(angle = 28,
  hjust = 1, size = 6))

pA2 <- ggplot(qc_plot_df, aes(x = source_dataset, y = nCount_RNA, fill = source_dataset)) + geom_violin(scale = "width",
  trim = TRUE, linewidth = 0.25) + geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.25,
  fill = "white") + scale_fill_manual(values = dataset_cols) + labs(x = NULL, y = "UMI counts") + theme_pub +
  theme(legend.position = "none", axis.text.x = element_text(angle = 28, hjust = 1, size = 6))

pA <- pA1 + pA2 + plot_layout(widths = c(1, 1))

find_contam_col <- function(md) {
    nm <- colnames(md)
    hits <- nm[grepl("decontx.*contam|contam.*decontx|^contamination$|contamination", nm, ignore.case = TRUE)]
    if (length(hits) == 0) {
        return(NA_character_)
    }
    numeric_hits <- hits[vapply(md[, hits, drop = FALSE], is.numeric, logical(1))]
    if (length(numeric_hits) == 0) {
        return(NA_character_)
    }
    numeric_hits[1]
}

decontx_list <- list()

for (ds in names(decontx_files)) {
    x <- readRDS(decontx_files[[ds]])
    mdx <- x[[]]
    cc <- find_contam_col(mdx)
    if (is.na(cc)) {
        stop("未在 ", basename(decontx_files[[ds]]), " 中找到DecontX contamination字段。")
    }
    decontx_list[[ds]] <- data.frame(cell = Cells(x), source_dataset = ds, contamination = as.numeric(mdx[[cc]]), 
        stringsAsFactors = FALSE)
    rm(x, mdx)
    gc(verbose = FALSE)
}

decontx_df <- bind_rows(decontx_list) %>% filter(is.finite(contamination))

rm(decontx_list)
gc(verbose = FALSE)

if (quantile(decontx_df$contamination, 0.95, na.rm = TRUE) <= 1.5) {
  decontx_df$contamination_percent <- decontx_df$contamination * 100
} else {
  decontx_df$contamination_percent <- decontx_df$contamination
}

decontx_df$source_dataset <- factor(decontx_df$source_dataset, levels = dataset_levels)

decontx_summary <- decontx_df %>% group_by(source_dataset) %>% summarise(n_cells = n(), median_contamination_percent = median(contamination_percent, 
    na.rm = TRUE), p95_contamination_percent = quantile(contamination_percent, 0.95, na.rm = TRUE), .groups = "drop")

write.csv(decontx_summary, file.path(out_dir, "14_SuppFigureS1_DecontX_summary_v07.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

set.seed(1234)

decontx_plot_df <- decontx_df %>% group_by(source_dataset) %>% group_modify(~slice_sample(.x, n = min(6000, 
    nrow(.x)), replace = FALSE)) %>% ungroup()

pB <- ggplot(decontx_plot_df, aes(x = source_dataset, y = contamination_percent, fill = source_dataset)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) + geom_boxplot(width = 0.13, outlier.shape = NA,
  fill = "white", linewidth = 0.25) + scale_fill_manual(values = dataset_cols) + labs(tag = "B", title = "DecontX contamination",
  x = NULL, y = "Estimated contamination (%)") + theme_pub + theme(legend.position = "none", axis.text.x = element_text(angle = 28,
  hjust = 1, size = 6))

doublet_raw <- read.csv(doublet_csv, stringsAsFactors = FALSE, check.names = FALSE)

required_doublet_cols <- c("dataset", "n_cells_before", "n_doublets", "n_singlets")

if (!all(required_doublet_cols %in% colnames(doublet_raw))) {
  stop("04 doublet summary缺少必要字段：", paste(setdiff(required_doublet_cols, colnames(doublet_raw)),
    collapse = ", "))
}

doublet_summary <- doublet_raw %>% filter(dataset %in% c("17702045", "GSE162610", "GSE196928_Control")) %>% 
    group_by(dataset) %>% summarise(input_cells = sum(n_cells_before, na.rm = TRUE), doublets = sum(n_doublets, 
    na.rm = TRUE), singlets = sum(n_singlets, na.rm = TRUE), .groups = "drop") %>% mutate(source_dataset = ifelse(dataset == 
    "GSE196928_Control", "GSE196928", dataset), singlet_percent = singlets/input_cells * 100, doublet_percent = doublets/input_cells * 
    100) %>% select(source_dataset, input_cells, singlets, doublets, singlet_percent, doublet_percent)

if (nrow(doublet_summary) != 3 || any(!is.finite(doublet_summary$doublet_percent))) {
  stop("S1C doublet summary异常，请检查04_doublet_3scRNA_summary_v03.csv。")
}

doublet_summary$source_dataset <- factor(doublet_summary$source_dataset, levels = dataset_levels)

write.csv(doublet_summary, file.path(out_dir, "14_SuppFigureS1_doublet_summary_v07.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

doublet_plot_df <- doublet_summary %>% select(source_dataset, Singlet = singlet_percent, Doublet = doublet_percent) %>% 
    pivot_longer(cols = c("Singlet", "Doublet"), names_to = "classification", values_to = "percent")

doublet_plot_df$classification <- factor(doublet_plot_df$classification, levels = c("Singlet", "Doublet"))

pC <- ggplot(doublet_plot_df, aes(x = source_dataset, y = percent, fill = classification)) + geom_col(width = 0.68,
  color = "white", linewidth = 0.25) + scale_fill_manual(values = c(Singlet = "#BDBDBD", Doublet = "#D95F5F")) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x,
    "%"), expand = c(0, 0)) + labs(tag = "C", title = "Doublet removal", x = NULL, y = "Cell proportion") +
  theme_pub + theme(legend.position = "top", legend.key.height = grid::unit(0.1, "in"), legend.key.width = grid::unit(0.1,
  "in"), axis.text.x = element_text(angle = 28, hjust = 1, size = 6))

if (file.exists(pre_cache_rds)) {
    pre_df_full <- readRDS(pre_cache_rds)
} else {
    obj06 <- readRDS(harmony06_file)
    if (!"pca3000" %in% names(obj06@reductions)) {
        stop("06对象中没有 pca3000 reduction。")
    }
    if (ncol(Embeddings(obj06, "pca3000")) < 30) {
        stop("06对象的pca3000少于30 PCs。")
    }
    overlap_cells <- intersect(final_cells, Cells(obj06))
    if (length(overlap_cells) != length(final_cells)) {
        stop("06对象与final atlas细胞集合不一致，", "为避免before/after比较失真，脚本停止。\n", 
            "final cells = ", length(final_cells), "; overlap = ", length(overlap_cells))
    }
    set.seed(1234)
    obj06 <- RunUMAP(object = obj06, reduction = "pca3000", dims = 1:30, reduction.name = "umap_preharmony_s1", 
        reduction.key = "preUMAP_", seed.use = 1234, verbose = TRUE)
    pre_mat <- Embeddings(obj06, "umap_preharmony_s1")[, 1:2, drop = FALSE]
    pre_df_full <- data.frame(cell = rownames(pre_mat), UMAP_1 = pre_mat[, 1], UMAP_2 = pre_mat[, 2], 
        source_dataset = factor(source_map[rownames(pre_mat)], levels = dataset_levels), stringsAsFactors = FALSE)
    saveRDS(pre_df_full, pre_cache_rds, compress = FALSE)
    write.csv(pre_df_full, pre_cache_csv, row.names = FALSE, fileEncoding = "UTF-8")
    rm(obj06, pre_mat)
    gc(verbose = FALSE)
}

needed_pre_cols <- c("cell", "UMAP_1", "UMAP_2", "source_dataset")

if (!all(needed_pre_cols %in% colnames(pre_df_full))) {
  stop("pre-Harmony缓存结构异常。")
}

common_cells <- Reduce(intersect, list(final_cells, pre_df_full$cell, post_df_full$cell))

if (length(common_cells) != length(final_cells)) {
  stop("pre/post UMAP未覆盖全部final atlas细胞。")
}

set.seed(1234)

plot_cells <- if (length(common_cells) > 50000) {
  sample(common_cells, 50000)
} else {
  common_cells
}

pre_df <- pre_df_full[match(plot_cells, pre_df_full$cell), , drop = FALSE]

post_df <- post_df_full[match(plot_cells, post_df_full$cell), , drop = FALSE]

pre_df$source_dataset <- factor(as.character(pre_df$source_dataset), levels = dataset_levels)

post_df$source_dataset <- factor(as.character(post_df$source_dataset), levels = dataset_levels)

pD1 <- ggplot(pre_df, aes(x = UMAP_1, y = UMAP_2, color = source_dataset)) + ggrastr::geom_point_rast(size = 0.18,
  alpha = 0.72, raster.dpi = 600) + scale_color_manual(values = dataset_cols, drop = FALSE) + labs(tag = "D",
  title = "Before Harmony", x = "UMAP 1", y = "UMAP 2") + theme_pub + theme(legend.position = "none",
  axis.text = element_blank(), axis.ticks = element_blank())

pD2 <- ggplot(post_df, aes(x = UMAP_1, y = UMAP_2, color = source_dataset)) + ggrastr::geom_point_rast(size = 0.18,
  alpha = 0.72, raster.dpi = 600) + scale_color_manual(values = dataset_cols, drop = FALSE) + labs(title = "After Harmony",
  x = "UMAP 1", y = "UMAP 2") + theme_pub + theme(axis.text = element_blank(), axis.ticks = element_blank(),
  legend.position = "right")

pD <- pD1 + pD2 + plot_layout(guides = "collect", widths = c(1, 1.12))

panel_e_stage_file <- file.path(project_dir, "results", "104_Harmony_source_sensitivity_v03", "01_104_source_stage_treatment_coverage_v03.csv")
panel_e_myeloid_file <- file.path(project_dir, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
panel_e_unit_cache <- file.path(out_dir, "14_SuppFigureS1_panelE_reporting_unit_composition_v07.csv")
panel_e_stage_out <- file.path(out_dir, "14_SuppFigureS1_panelE_source_stage_composition_v07.csv")
panel_e_audit_out <- file.path(out_dir, "14_SuppFigureS1_panelE_replication_audit_v07.csv")

if (!file.exists(panel_e_stage_file)) {
  stop("Panel E 缺少 Step104 来源分层结果：\n", panel_e_stage_file, "\n该文件是先前已计算的结果，不会自动重新运行分析。")
}

main_stage_order <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")
focus_e <- "Clec7a+ lipid-associated Microglia"

stage_raw <- read.csv(panel_e_stage_file, stringsAsFactors = FALSE, check.names = FALSE)
required_stage_cols <- c("source_dataset", "timepoint_main", "n_myeloid_cells", "n_frozen_focus", "pct_frozen_focus")
missing_stage <- setdiff(required_stage_cols, colnames(stage_raw))
if (length(missing_stage)) {
  stop("Panel E Step104 CSV 缺少字段：", paste(missing_stage, collapse = ", "))
}
stage_e <- stage_raw[stage_raw$source_dataset %in% dataset_levels & stage_raw$timepoint_main %in% main_stage_order,
  required_stage_cols, drop = FALSE]
if (!nrow(stage_e) || anyDuplicated(stage_e[c("source_dataset", "timepoint_main")])) {
  stop("Panel E Step104 来源×阶段记录缺失或重复。")
}
stage_e$n_myeloid_cells <- as.integer(stage_e$n_myeloid_cells)
stage_e$n_frozen_focus <- as.integer(stage_e$n_frozen_focus)
stage_e$pct_frozen_focus <- as.numeric(stage_e$pct_frozen_focus)
if (anyNA(stage_e$n_myeloid_cells) || any(stage_e$n_myeloid_cells <= 0) || anyNA(stage_e$n_frozen_focus) ||
  any(stage_e$n_frozen_focus < 0) || any(stage_e$n_frozen_focus > stage_e$n_myeloid_cells) || any(abs(100 *
  stage_e$n_frozen_focus/stage_e$n_myeloid_cells - stage_e$pct_frozen_focus) > 0.01)) {
  stop("Panel E Step104 来源分层分子/分母/比例不一致。")
}

cache_is_current <- file.exists(panel_e_unit_cache) && (!file.exists(panel_e_myeloid_file) || file.info(panel_e_unit_cache)$mtime >=
  file.info(panel_e_myeloid_file)$mtime)

if (cache_is_current) {
    unit_e <- read.csv(panel_e_unit_cache, stringsAsFactors = FALSE, check.names = FALSE)
} else {
    if (!file.exists(panel_e_myeloid_file)) {
        stop("Panel E 既没有可用缓存，也未找到 Step18 注释对象：\n", panel_e_myeloid_file)
    }
    myeloid_obj_e <- readRDS(panel_e_myeloid_file)
    md_e <- myeloid_obj_e[[]]
    rm(myeloid_obj_e)
    gc(verbose = FALSE)
    needed_md_e <- c("source_dataset", "timepoint_main", "sample_id", "myeloid_state")
    missing_md_e <- setdiff(needed_md_e, colnames(md_e))
    if (length(missing_md_e)) {
        stop("Step18 对象缺少 Panel E 必需元数据：", paste(missing_md_e, collapse = ", "))
    }
    md_e <- md_e[, needed_md_e, drop = FALSE]
    for (nm in needed_md_e) md_e[[nm]] <- trimws(as.character(md_e[[nm]]))
    if (anyNA(md_e) || any(!nzchar(as.matrix(md_e))) || any(!md_e$source_dataset %in% dataset_levels)) {
        stop("Panel E 发现缺失或不可识别的来源、时间、文库或状态注释。")
    }
    md_e$is_focus <- md_e$myeloid_state == focus_e
    if (!any(md_e$is_focus)) 
        stop("Step18 中找不到冻结的 Clec7a 目标状态标签。")
    unit_e <- md_e %>% group_by(source_dataset, timepoint_main, sample_id) %>% summarise(n_myeloid_cells = dplyr::n(), 
        n_frozen_focus = sum(is_focus), .groups = "drop") %>% mutate(pct_frozen_focus = 100 * n_frozen_focus/n_myeloid_cells) %>% 
        as.data.frame()
    rm(md_e)
    gc(verbose = FALSE)
    write.csv(unit_e, panel_e_unit_cache, row.names = FALSE, fileEncoding = "UTF-8")
}

unit_req_e <- c("source_dataset", "timepoint_main", "sample_id", "n_myeloid_cells", "n_frozen_focus",
  "pct_frozen_focus")
if (length(setdiff(unit_req_e, names(unit_e))) || anyDuplicated(unit_e[c("source_dataset", "timepoint_main",
  "sample_id")])) {
  stop("Panel E unit 缓存字段缺失或报告单位重复。")
}
unit_e <- unit_e[unit_e$source_dataset %in% dataset_levels & unit_e$timepoint_main %in% main_stage_order,
  , drop = FALSE]
if (!nrow(unit_e)) stop("Panel E 无主图六阶段对应的报告单位数据。")
unit_e$n_myeloid_cells <- as.integer(unit_e$n_myeloid_cells)
unit_e$n_frozen_focus <- as.integer(unit_e$n_frozen_focus)
unit_e$pct_frozen_focus <- as.numeric(unit_e$pct_frozen_focus)

unit_totals_e <- unit_e %>% group_by(source_dataset, timepoint_main) %>% summarise(n_units = n(), n_sum_myeloid = sum(n_myeloid_cells), 
    n_sum_focus = sum(n_frozen_focus), .groups = "drop")
check_e <- merge(stage_e, unit_totals_e, by = c("source_dataset", "timepoint_main"), all = TRUE)
if (nrow(check_e) != nrow(stage_e) || anyNA(check_e$n_sum_myeloid) || anyNA(check_e$n_myeloid_cells) ||
  any(check_e$n_sum_myeloid != check_e$n_myeloid_cells) || any(check_e$n_sum_focus != check_e$n_frozen_focus)) {
  write.csv(check_e, file.path(out_dir, "14_SuppFigureS1_panelE_MISMATCH_AUDIT_v07.csv"), row.names = FALSE)
  stop("Panel E Step18 与 Step104 的来源×阶段计数不一致：", "请检查 Step18 注释/缓存的版本，不可直接出图。")
}

unit_e$reporting_unit_type <- ifelse(unit_e$source_dataset == "17702045", "Stage-defined stratum", ifelse(unit_e$source_dataset ==
  "GSE196928", "Pooled Vehicle library", "GSM library / possible pool"))
unit_e$sample_id_provenance <- ifelse(unit_e$source_dataset == "17702045", "Derived stage label; NOT independent replicate",
  ifelse(unit_e$source_dataset == "GSE196928", "One pooled tissue library per condition", "GSM/library ID; animal pooling varies"))

expected_units_e <- data.frame(source_dataset = c(rep("GSE162610", 4), rep("GSE196928", 3)), timepoint_main = c("Control",
  "1 dpi", "3 dpi", "7 dpi", "Control", "7 dpi", "28 dpi"), expected_n = c(3L, 3L, 2L, 2L, 1L, 1L,
  1L))
audit_e <- merge(check_e, expected_units_e, by = c("source_dataset", "timepoint_main"), all.x = TRUE)
if (any(!is.na(audit_e$expected_n) & audit_e$n_units != audit_e$expected_n)) {
  stop("Panel E 中 GSE162610 或 GSE196928 文库数与 Table S2 不一致。")
}
audit_e$unit_definition <- ifelse(audit_e$source_dataset == "17702045", "Stage-defined cells; biological n not recoverable",
  ifelse(audit_e$source_dataset == "GSE196928", "One pooled Vehicle library per condition; no animal-level variance",
    "Source GSM libraries; some pooled animals"))
write.csv(audit_e, panel_e_audit_out, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(stage_e, panel_e_stage_out, row.names = FALSE, fileEncoding = "UTF-8")

unit_e$source_dataset <- factor(unit_e$source_dataset, levels = dataset_levels)
stage_e$source_dataset <- factor(stage_e$source_dataset, levels = dataset_levels)
unit_e$timepoint_main <- factor(unit_e$timepoint_main, levels = main_stage_order)
stage_e$timepoint_main <- factor(stage_e$timepoint_main, levels = main_stage_order)
unit_e$reporting_unit_type <- factor(unit_e$reporting_unit_type, levels = c("Stage-defined stratum",
  "GSM library / possible pool", "Pooled Vehicle library"))

pE <- ggplot() + geom_line(data = stage_e, aes(x = timepoint_main, y = pct_frozen_focus, group = source_dataset),
  color = "#9BA3AD", linewidth = 0.43) + geom_point(data = stage_e, aes(x = timepoint_main, y = pct_frozen_focus),
  shape = 95, color = "#363636", size = 5.2) + geom_point(data = unit_e, aes(x = timepoint_main, y = pct_frozen_focus,
  color = source_dataset, shape = reporting_unit_type), position = position_jitter(width = 0.105, height = 0,
  seed = 1234), size = 1.85, alpha = 0.91, stroke = 0.24) + facet_grid(~source_dataset, drop = FALSE) +
  scale_color_manual(values = dataset_cols, guide = "none") + scale_shape_manual(values = c(`Stage-defined stratum` = 15,
  `GSM library / possible pool` = 16, `Pooled Vehicle library` = 17), drop = FALSE) + scale_y_continuous(limits = c(0,
  100), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.015,
  0.035))) + scale_x_discrete(drop = FALSE) + labs(tag = "E", title = "Source- and reporting-unit-stratified Clec7a+ composition",
  x = NULL, y = "Clec7a+ / myeloid cells (%)", shape = NULL) + theme_pub + theme(strip.background = element_rect(fill = "#F2F3F5",
  color = NA), strip.text = element_text(size = 6.4, face = "bold", color = "black"), panel.spacing.x = grid::unit(0.12,
  "in"), plot.subtitle = element_blank(), plot.caption = element_blank(), axis.text.x = element_text(angle = 48,
  hjust = 1, vjust = 1, size = 5.7), axis.title.y = element_text(size = 7), legend.position = "bottom",
  legend.direction = "horizontal", legend.text = element_text(size = 5.9), legend.key.width = grid::unit(0.1,
    "in"), legend.key.height = grid::unit(0.1, "in"), legend.margin = margin(0, 0, 0, 0), legend.spacing.y = grid::unit(0.01,
    "in")) + guides(shape = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2)))

layout_s1_v08 <- "AABB\nCCDD\nEEEE"
combined <- patchwork::wrap_plots(A = patchwork::wrap_elements(full = pA), B = pB, C = pC, D = patchwork::wrap_elements(full = pD),
  E = pE, design = layout_s1_v08) + patchwork::plot_layout(widths = c(1, 1, 1.08, 1.08), heights = c(1.02,
  1.02, 1.14)) & theme(text = element_text(family = FONT_FAMILY))

stopifnot(length(strsplit(layout_s1_v08, "\n", fixed = TRUE)[[1]]) == 3L, length(unique(strsplit(gsub("\n",
  "", layout_s1_v08), "", fixed = TRUE)[[1]])) == 5L)

pdf_file <- file.path(out_dir, "14_SuppFigureS1_QC_integration_3scRNA_v09.pdf")

svg_file <- file.path(out_dir, "14_SuppFigureS1_QC_integration_3scRNA_v09.svg")

png_file <- file.path(out_dir, "14_SuppFigureS1_QC_integration_3scRNA_v09.png")

ggsave(filename = pdf_file, plot = combined, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, units = "in",
  device = cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = svg_file, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN)
    print(combined)
    dev.off()
} else {
    ggsave(filename = svg_file, plot = combined, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, units = "in")
}

if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = png_file, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, units = "in", res = 600, 
        background = "white")
    print(combined)
    dev.off()
} else {
    ggsave(filename = png_file, plot = combined, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, units = "in", 
        dpi = 600, bg = "white")
}

if (USE_SHOWTEXT) {
  showtext::showtext_auto(enable = FALSE)
}
