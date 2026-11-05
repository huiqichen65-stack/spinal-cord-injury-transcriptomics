source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
UMAP_REDUCTION <- "umap_myeloid_harmony30"

MAIN_TIMES <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")

FULL_TIMES <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")

STATE_LEVELS <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia",
  "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage",
  "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

STATE_COLORS <- c(`Homeostatic Microglia` = "#3B6FB6", `Clec7a+ lipid-associated Microglia` = "#E08B2C",
  `Activated inflammatory Microglia` = "#C94C4C", `Ccl12+ acute Microglia` = "#D66BA0", `Tnf+ inflammatory Microglia` = "#8E63B6",
  `Proliferating Microglia` = "#2F9E8F", `MHC-II+ antigen-presenting Macrophage` = "#4F8A3B", `Gpnmb+ tissue-repair Macrophage` = "#8A6D5A",
  `Acute-response Macrophage` = "#D9A514")

ROGUE_MIN_CELLS <- 50
ROGUE_MAX_CELLS_PER_STATE_DATASET <- 800
ROGUE_SPAN <- 0.5

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "dplyr", "tidyr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), "\nPlease install them before rerunning.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

HAS_ROGUE <- requireNamespace("ROGUE", quietly = TRUE)

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")

if (file.exists(expected_input)) {
    in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
} else {
    hits <- list.files(project_dir, pattern = "^18_myeloid_annotation_cleanup_3scRNA_v01\\.rds$", recursive = TRUE, 
        full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        stop("Input file not found:\n", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
    }
    in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "19_myeloid_state_robustness_temporal_dynamics")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "19_myeloid_state_robustness_temporal_dynamics_v02.log")

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

obj$timepoint_main <- factor(as.character(obj$timepoint_main), levels = FULL_TIMES)

if (anyNA(obj$myeloid_state)) {
  stop("Unexpected myeloid-state label detected.")
}

capture.output(print(table(obj$myeloid_state)), file = log_file, append = TRUE)

make_temporal_table <- function(meta, time_levels) {
    x <- meta %>% filter(as.character(timepoint_main) %in% time_levels) %>% mutate(timepoint_main = factor(as.character(timepoint_main), 
        levels = time_levels), myeloid_state = factor(as.character(myeloid_state), levels = STATE_LEVELS)) %>% 
        count(timepoint_main, myeloid_state, name = "observed", .drop = FALSE) %>% tidyr::complete(timepoint_main = factor(time_levels, 
        levels = time_levels), myeloid_state = factor(STATE_LEVELS, levels = STATE_LEVELS), fill = list(observed = 0))
    time_totals <- x %>% group_by(timepoint_main) %>% summarise(time_total = sum(observed), .groups = "drop")
    state_totals <- x %>% group_by(myeloid_state) %>% summarise(state_total = sum(observed), .groups = "drop")
    grand_total <- sum(x$observed)
    x %>% left_join(time_totals, by = "timepoint_main") %>% left_join(state_totals, by = "myeloid_state") %>% 
        mutate(grand_total = grand_total, fraction_within_time = ifelse(time_total > 0, observed/time_total, 
            NA_real_), fraction_within_state = ifelse(state_total > 0, observed/state_total, NA_real_), 
            expected = ifelse(grand_total > 0, time_total * state_total/grand_total, NA_real_), Roe = ifelse(!is.na(expected) & 
                expected > 0, observed/expected, NA_real_), log2_Roe = ifelse(!is.na(Roe) & Roe > 0, 
                log2(Roe), NA_real_))
}

meta_all <- obj[[]] %>% mutate(myeloid_state = factor(as.character(myeloid_state), levels = STATE_LEVELS), 
    timepoint_main = factor(as.character(timepoint_main), levels = FULL_TIMES))

full_temporal <- make_temporal_table(meta_all, FULL_TIMES)

main_temporal <- make_temporal_table(meta_all, MAIN_TIMES)

composition_cols <- c("timepoint_main", "myeloid_state", "observed", "time_total", "state_total", "grand_total",
  "fraction_within_time", "fraction_within_state")

roe_cols <- c("timepoint_main", "myeloid_state", "observed", "expected", "Roe", "log2_Roe")

write.csv(full_temporal[, composition_cols], file.path(out_dir, "19_myeloid_temporal_composition_full_v02.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")

write.csv(main_temporal[, composition_cols], file.path(out_dir, "19_myeloid_temporal_composition_main_v02.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")

write.csv(full_temporal[, roe_cols], file.path(out_dir, "19_myeloid_Roe_full_v02.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

write.csv(main_temporal[, roe_cols], file.path(out_dir, "19_myeloid_Roe_main_v02.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

rogue_long <- data.frame(source_dataset = character(), myeloid_state = character(), ROGUE = numeric(),
  n_cells_used = integer(), stringsAsFactors = FALSE)

if (!HAS_ROGUE) {

  log_msg("ROGUE package is not installed. ", "ROGUE section skipped; temporal analyses still completed.")

  write.csv(data.frame(status = "SKIPPED: install ROGUE and rerun Step 19 v02.", stringsAsFactors = FALSE),
    file.path(out_dir, "19_myeloid_ROGUE_v02.csv"), row.names = FALSE)

} else {

  set.seed(SEED)

  rogue_meta <- obj[[]] %>%
    mutate(cell = rownames(obj[[]]), source_dataset = as.character(source_dataset), myeloid_state = as.character(myeloid_state),
      group_id = paste(source_dataset, myeloid_state, sep = "|||"))

  sampled_cells <- unlist(lapply(split(rogue_meta$cell, rogue_meta$group_id), function(cells) {
    if (length(cells) > ROGUE_MAX_CELLS_PER_STATE_DATASET) {
      sample(cells, ROGUE_MAX_CELLS_PER_STATE_DATASET)
    } else {
      cells
    }
  }), use.names = FALSE)

  rogue_meta_sub <- rogue_meta[match(sampled_cells, rogue_meta$cell), , drop = FALSE]

  counts_rogue <- GetAssayData(obj, assay = "RNA", layer = "counts")[, sampled_cells, drop = FALSE]

  gene_keep <- Matrix::rowSums(counts_rogue > 0) >= 10

  counts_rogue <- counts_rogue[gene_keep, , drop = FALSE]

  rogue_result <- tryCatch({
    ROGUE::rogue(expr = counts_rogue, labels = rogue_meta_sub$myeloid_state, samples = rogue_meta_sub$source_dataset,
      platform = "UMI", min.cell.n = ROGUE_MIN_CELLS, remove.outlier.n = 2, span = ROGUE_SPAN,
      r = 1, mt.method = "fdr")
  }, error = function(e) {
    log_msg("ROGUE failed safely: ", conditionMessage(e))
    NULL
  })

  if (is.null(rogue_result)) {

    write.csv(data.frame(status = "ROGUE calculation failed; inspect Step 19 v02 log.", stringsAsFactors = FALSE),
      file.path(out_dir, "19_myeloid_ROGUE_v02.csv"), row.names = FALSE)

  } else {

    rogue_df <- as.data.frame(rogue_result, check.names = FALSE)

    rogue_df$source_dataset <- rownames(rogue_df)

    rownames(rogue_df) <- NULL

    rogue_long <- rogue_df %>%
      pivot_longer(cols = -source_dataset, names_to = "myeloid_state", values_to = "ROGUE") %>%
      mutate(myeloid_state = as.character(myeloid_state))

    n_used <- rogue_meta_sub %>%
      count(source_dataset, myeloid_state, name = "n_cells_used")

    rogue_long <- rogue_long %>%
      left_join(n_used, by = c("source_dataset", "myeloid_state")) %>%
      arrange(match(myeloid_state, STATE_LEVELS), source_dataset)

    write.csv(rogue_long, file.path(out_dir, "19_myeloid_ROGUE_v02.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  }

  rm(counts_rogue, rogue_meta, rogue_meta_sub, sampled_cells)

  gc(verbose = FALSE)
}

diagnostic_pdf <- file.path(out_dir, "19_myeloid_state_robustness_temporal_diagnostics_v02.pdf")

pdf(diagnostic_pdf, width = 12, height = 8.5, onefile = TRUE)

p_umap <- DimPlot(obj, reduction = UMAP_REDUCTION, group.by = "myeloid_state", cols = STATE_COLORS, pt.size = 0.1,
  raster = TRUE) + ggtitle("Consolidated myeloid-state atlas") + theme_classic(base_size = 10)

print(p_umap)

p_comp_main <- ggplot(main_temporal, aes(x = timepoint_main, y = fraction_within_time, fill = myeloid_state)) +
  geom_col(width = 0.82) + scale_fill_manual(values = STATE_COLORS, drop = FALSE) + scale_y_continuous(labels = function(x) {
  paste0(round(x * 100), "%")
}, expand = c(0, 0)) + labs(title = "Temporal composition of myeloid states | main timeline", x = NULL,
  y = "Fraction of myeloid cells", fill = "Myeloid state") + theme_classic(base_size = 10)

print(p_comp_main)

p_comp_full <- ggplot(full_temporal, aes(x = timepoint_main, y = fraction_within_time, fill = myeloid_state)) +
  geom_col(width = 0.82) + scale_fill_manual(values = STATE_COLORS, drop = FALSE) + scale_y_continuous(labels = function(x) {
  paste0(round(x * 100), "%")
}, expand = c(0, 0)) + labs(title = "Temporal composition of myeloid states | full timeline", x = NULL,
  y = "Fraction of myeloid cells", fill = "Myeloid state") + theme_classic(base_size = 10)

print(p_comp_full)

roe_main_plot <- main_temporal %>% mutate(plot_log2_Roe = pmax(pmin(log2_Roe, 2), -2))

p_roe_main <- ggplot(roe_main_plot, aes(x = timepoint_main, y = myeloid_state, fill = plot_log2_Roe)) +
  geom_tile(color = "white", linewidth = 0.25) + geom_text(aes(label = ifelse(is.na(Roe), "", sprintf("%.2f",
  Roe))), size = 2.5) + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D94F4F", midpoint = 0,
  limits = c(-2, 2), name = "log2(R/oe)") + labs(title = "Temporal enrichment of myeloid states | main timeline",
  x = NULL, y = NULL) + theme_classic(base_size = 9) + theme(axis.ticks = element_blank())

print(p_roe_main)

roe_full_plot <- full_temporal %>% mutate(plot_log2_Roe = pmax(pmin(log2_Roe, 2), -2))

p_roe_full <- ggplot(roe_full_plot, aes(x = timepoint_main, y = myeloid_state, fill = plot_log2_Roe)) +
  geom_tile(color = "white", linewidth = 0.25) + geom_text(aes(label = ifelse(is.na(Roe), "", sprintf("%.2f",
  Roe))), size = 2.4) + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D94F4F", midpoint = 0,
  limits = c(-2, 2), name = "log2(R/oe)") + labs(title = "Temporal enrichment of myeloid states | full timeline",
  x = NULL, y = NULL) + theme_classic(base_size = 9) + theme(axis.ticks = element_blank())

print(p_roe_full)

if (nrow(rogue_long) > 0 && "ROGUE" %in% colnames(rogue_long)) {
    rogue_plot_df <- rogue_long %>% filter(!is.na(ROGUE)) %>% mutate(myeloid_state = factor(myeloid_state, 
        levels = STATE_LEVELS))
    p_rogue <- ggplot(rogue_plot_df, aes(x = myeloid_state, y = ROGUE)) + geom_point(aes(shape = source_dataset), 
        size = 2.8, position = position_jitter(width = 0.08, height = 0)) + coord_cartesian(ylim = c(0, 
        1)) + labs(title = "ROGUE purity of consolidated states by source dataset", x = NULL, y = "ROGUE", 
        shape = "Dataset") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 35, 
        hjust = 1))
    print(p_rogue)
}

dev.off()

out_rds <- file.path(out_dir, "19_myeloid_state_robustness_temporal_dynamics_3scRNA_v02.rds")

saveRDS(obj, out_rds, compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)
