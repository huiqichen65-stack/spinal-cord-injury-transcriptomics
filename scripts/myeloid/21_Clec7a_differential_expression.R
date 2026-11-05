source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
set.seed(SEED)

TARGET_STATE <- "Clec7a+ lipid-associated Microglia"

LOG2FC_CUTOFF <- 0.25
PADJ_CUTOFF <- 0.05
MIN_PCT <- 0.1

required_pkgs <- c("Seurat", "SeuratObject", "dplyr", "ggplot2", "ggrepel")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), "\nPlease install them before rerunning.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

if (file.exists(expected_input)) {
    in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
} else {
    hits <- list.files(project_dir, pattern = "^20_AUCell_lipid_metabolism_scoring_3scRNA_v03\\.rds$", 
        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        stop("Input file not found:\n", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")
    }
    in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "21_Clec7a_DEG_volcano_myeloid")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "21_Clec7a_DEG_volcano_myeloid_v01.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

if (!"myeloid_state" %in% colnames(obj[[]])) {
  stop("Metadata column 'myeloid_state' is missing.")
}

state_values <- unique(as.character(obj$myeloid_state))

if (!TARGET_STATE %in% state_values) {
  stop("Target state not found in object:\n", TARGET_STATE)
}

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

Idents(obj) <- obj$myeloid_state

n_target <- sum(Idents(obj) == TARGET_STATE)

n_rest <- ncol(obj) - n_target

capture.output(print(sort(table(obj$myeloid_state), decreasing = TRUE)), file = log_file, append = TRUE)

deg <- FindMarkers(object = obj, ident.1 = TARGET_STATE, ident.2 = NULL, assay = "RNA", slot = "data",
  test.use = "wilcox", logfc.threshold = 0, min.pct = MIN_PCT, only.pos = FALSE, return.thresh = 1,
  verbose = TRUE)

deg$gene <- rownames(deg)
rownames(deg) <- NULL

fc_col <- NULL

if ("avg_log2FC" %in% colnames(deg)) {
  fc_col <- "avg_log2FC"
} else if ("avg_logFC" %in% colnames(deg)) {
  fc_col <- "avg_logFC"
} else {
  stop("Could not find avg_log2FC or avg_logFC in FindMarkers output.")
}

if (fc_col == "avg_logFC") {
  deg$avg_log2FC <- deg$avg_logFC
}

preferred_cols <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")

other_cols <- setdiff(colnames(deg), preferred_cols)

deg <- deg[, c(preferred_cols[preferred_cols %in% colnames(deg)], other_cols), drop = FALSE]

deg <- deg %>% mutate(comparison = "Clec7a+ lipid-associated Microglia vs all other myeloid states", 
    significance = case_when(p_val_adj < PADJ_CUTOFF & avg_log2FC >= LOG2FC_CUTOFF ~ "Up", p_val_adj < 
        PADJ_CUTOFF & avg_log2FC <= -LOG2FC_CUTOFF ~ "Down", TRUE ~ "NS"), direction = case_when(avg_log2FC > 
        0 ~ "Higher_in_Clec7a", avg_log2FC < 0 ~ "Lower_in_Clec7a", TRUE ~ "No_change")) %>% arrange(p_val_adj, 
    desc(abs(avg_log2FC)))

deg_sig <- deg %>% filter(significance != "NS")

deg_up <- deg %>% filter(significance == "Up") %>% arrange(p_val_adj, desc(avg_log2FC))

deg_down <- deg %>% filter(significance == "Down") %>% arrange(p_val_adj, avg_log2FC)

write.csv(deg, file.path(out_dir, "21_Clec7a_vs_all_myeloid_DEG_all_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(deg_sig, file.path(out_dir, "21_Clec7a_vs_all_myeloid_DEG_significant_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

write.csv(deg_up, file.path(out_dir, "21_Clec7a_vs_all_myeloid_DEG_up_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(deg_down, file.path(out_dir, "21_Clec7a_vs_all_myeloid_DEG_down_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

min_nonzero_padj <- suppressWarnings(min(deg$p_val_adj[deg$p_val_adj > 0 & is.finite(deg$p_val_adj)],
  na.rm = TRUE))

if (!is.finite(min_nonzero_padj)) {
  min_nonzero_padj <- 1e-300
}

deg$padj_for_plot <- deg$p_val_adj

deg$padj_for_plot[is.na(deg$padj_for_plot)] <- 1

deg$padj_for_plot[deg$padj_for_plot <= 0] <- max(min_nonzero_padj/10, 1e-300)

deg$neglog10_padj <- -log10(deg$padj_for_plot)

top_up_labels <- deg_up %>% arrange(p_val_adj, desc(avg_log2FC)) %>% slice_head(n = 10) %>% pull(gene)

top_down_labels <- deg_down %>% arrange(p_val_adj, avg_log2FC) %>% slice_head(n = 10) %>% pull(gene)

priority_genes <- c("Clec7a", "Lpl", "Apoe", "Apoc4", "Trem2", "Abca1", "Abcg1", "Lrp1", "Cd36", "Npc1",
  "Npc2", "Soat1", "Plin2")

priority_sig <- deg %>% filter(gene %in% priority_genes, significance != "NS") %>% pull(gene)

label_genes <- unique(c(top_up_labels, top_down_labels, priority_sig))

volcano_df <- deg %>% mutate(label = ifelse(gene %in% label_genes, gene, NA_character_))

volcano_colors <- c(Down = "#3B6FB6", NS = "#BDBDBD", Up = "#D94F4F")

p_volcano <- ggplot(volcano_df, aes(x = avg_log2FC, y = neglog10_padj)) + geom_point(aes(color = significance),
  size = 1.15, alpha = 0.68) + geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed",
  linewidth = 0.45) + geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", linewidth = 0.45) +
  ggrepel::geom_text_repel(data = volcano_df %>%
    filter(!is.na(label)), aes(label = label), size = 3.1, max.overlaps = Inf, box.padding = 0.35,
    point.padding = 0.2, min.segment.length = 0, segment.size = 0.3, seed = SEED, show.legend = FALSE) +
  scale_color_manual(values = volcano_colors, breaks = c("Up", "Down", "NS"), drop = FALSE) + labs(title = "Clec7a+ lipid-associated Microglia vs all other myeloid states",
  subtitle = paste0("Wilcoxon | adjusted P < ", PADJ_CUTOFF, " | |log2FC| >= ", LOG2FC_CUTOFF, " | Up = ",
    nrow(deg_up), ", Down = ", nrow(deg_down)), x = "Average log2 fold change", y = expression(-log[10](adjusted ~
    P)), color = "DEG") + theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold",
  size = 12), plot.subtitle = element_text(size = 9), legend.title = element_text(face = "bold"))

pdf_file <- file.path(out_dir, "21_Clec7a_vs_all_myeloid_volcano_v01.pdf")

pdf(pdf_file, width = 8, height = 7)

print(p_volcano)

dev.off()

png_file <- file.path(out_dir, "21_Clec7a_vs_all_myeloid_volcano_v01.png")

ggsave(filename = png_file, plot = p_volcano, width = 8, height = 7, units = "in", dpi = 600, bg = "white")

if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file.path(out_dir, "21_Clec7a_vs_all_myeloid_volcano_v01.svg"), width = 8, height = 7)
    print(p_volcano)
    dev.off()
}

obj$Step21_DEG_target <- ifelse(as.character(obj$myeloid_state) == TARGET_STATE, "Clec7a_target", "Other_myeloid")

saveRDS(obj, file.path(out_dir, "21_Clec7a_DEG_volcano_myeloid_3scRNA_v01.rds"), compress = FALSE)

capture.output(print(deg_up %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>% slice_head(n = 20)), 
    file = log_file, append = TRUE)

capture.output(print(deg_down %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>% slice_head(n = 20)), 
    file = log_file, append = TRUE)

capture.output(sessionInfo(), file = log_file, append = TRUE)
