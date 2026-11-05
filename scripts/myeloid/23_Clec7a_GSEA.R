source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
set.seed(SEED)

TARGET_STATE <- "Clec7a+ lipid-associated Microglia"

PADJ_CUTOFF <- 0.05
MIN_GS_SIZE <- 10
MAX_GS_SIZE <- 500
TOP_POS <- 10
TOP_NEG <- 10
TOP_CURVES_EACH_DIRECTION <- 5

required_pkgs <- c("Seurat", "SeuratObject", "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "enrichplot",
  "ggplot2", "dplyr", "readr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), "\nInstall before rerunning.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

if (file.exists(expected_input)) {
    in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
} else {
    hits <- list.files(project_dir, pattern = "^20_AUCell_lipid_metabolism_scoring_3scRNA_v03\\.rds$", 
        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        stop("Cannot find input object:\n", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")
    }
    in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

kegg_dir <- file.path(project_dir, "reference", "KEGG")

find_local_kegg_file <- function(exact_name, pattern) {
    exact_path <- file.path(kegg_dir, exact_name)
    if (file.exists(exact_path)) {
        return(normalizePath(exact_path, winslash = "/", mustWork = TRUE))
    }
    hits <- list.files(kegg_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        return(NA_character_)
    }
    normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

gene_pathway_file <- find_local_kegg_file("mmu_gene_pathway.txt", "mmu.*gene.*pathway|gene.*pathway.*mmu|link.*pathway")

pathway_list_file <- find_local_kegg_file("mmu_pathway_list.txt", "mmu.*pathway.*list|pathway.*list.*mmu|list.*pathway")

if (is.na(gene_pathway_file)) {
  stop("Cannot find local KEGG gene-pathway file under:\n", kegg_dir)
}

if (is.na(pathway_list_file)) {
  stop("Cannot find local KEGG pathway-list file under:\n", kegg_dir)
}

out_dir <- file.path(project_dir, "results", "23_Clec7a_GSEA_KEGG")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "23_Clec7a_GSEA_KEGG_v02.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

if (!"myeloid_state" %in% colnames(obj[[]])) {
  stop("Metadata column 'myeloid_state' is missing.")
}

if (!TARGET_STATE %in% unique(as.character(obj$myeloid_state))) {
  stop("Target state not found:\n", TARGET_STATE)
}

obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) {
  log_msg("JoinLayers skipped/failed safely: ", conditionMessage(e))
  obj
})

DefaultAssay(obj) <- "RNA"

rna_layers <- tryCatch(Layers(obj[["RNA"]]), error = function(e) {
  character()
})

if (!"data" %in% rna_layers) {
    obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, 
        verbose = TRUE)
}

Idents(obj) <- obj$myeloid_state

n_target <- sum(Idents(obj) == TARGET_STATE)

n_other <- ncol(obj) - n_target

if (requireNamespace("presto", quietly = TRUE)) {
} else {
}

full_rank_deg <- FindMarkers(object = obj, ident.1 = TARGET_STATE, ident.2 = NULL, assay = "RNA", slot = "data",
  test.use = "wilcox", logfc.threshold = 0, min.pct = 0, min.cells.feature = 1, min.cells.group = 1,
  only.pos = FALSE, return.thresh = 1, densify = FALSE, verbose = TRUE)

full_rank_deg$gene <- rownames(full_rank_deg)

rownames(full_rank_deg) <- NULL

if ("avg_log2FC" %in% colnames(full_rank_deg)) {
    fc_col <- "avg_log2FC"
} else if ("avg_logFC" %in% colnames(full_rank_deg)) {
    full_rank_deg$avg_log2FC <- full_rank_deg$avg_logFC
    fc_col <- "avg_log2FC"
} else {
    stop("FindMarkers output does not contain avg_log2FC/avg_logFC.")
}

full_rank_deg <- full_rank_deg %>% dplyr::filter(!is.na(avg_log2FC), is.finite(avg_log2FC)) %>% dplyr::arrange(dplyr::desc(avg_log2FC))

write.csv(full_rank_deg, file.path(out_dir, "23_Clec7a_GSEA_fullrank_FindMarkers_v02.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

gene_map <- suppressMessages(clusterProfiler::bitr(unique(full_rank_deg$gene), fromType = "SYMBOL", toType = "ENTREZID",
  OrgDb = org.Mm.eg.db))

gene_map <- gene_map %>% dplyr::distinct(SYMBOL, ENTREZID)

rank_df <- full_rank_deg %>% dplyr::inner_join(gene_map, by = c(gene = "SYMBOL")) %>% dplyr::mutate(abs_fc = abs(avg_log2FC))

rank_df <- rank_df %>% dplyr::arrange(dplyr::desc(abs_fc)) %>% dplyr::distinct(ENTREZID, .keep_all = TRUE)

gene_list <- rank_df$avg_log2FC

names(gene_list) <- as.character(rank_df$ENTREZID)

gene_list <- sort(gene_list, decreasing = TRUE)

if (anyDuplicated(names(gene_list)) > 0) {
  stop("Duplicate ENTREZ IDs remain after deduplication.")
}

if (length(gene_list) < 1000) {
  stop("Too few genes after ID mapping: ", length(gene_list))
}

write.csv(rank_df %>% dplyr::arrange(dplyr::desc(avg_log2FC)), file.path(out_dir, "23_Clec7a_GSEA_ranked_genes_v02.csv"), 
    row.names = FALSE, fileEncoding = "UTF-8")

raw_link <- readr::read_tsv(gene_pathway_file, col_names = FALSE, show_col_types = FALSE, progress = FALSE,
  trim_ws = TRUE)

if (ncol(raw_link) < 2) {
  stop("KEGG gene-pathway file does not contain two tab-separated columns.")
}

raw_link <- raw_link[, 1:2, drop = FALSE]

colnames(raw_link) <- c("gene_raw", "pathway_raw")

first_gene_like <- mean(grepl("^mmu:", raw_link$gene_raw))

second_gene_like <- mean(grepl("^mmu:", raw_link$pathway_raw))

if (second_gene_like > first_gene_like) {
    tmp <- raw_link$gene_raw
    raw_link$gene_raw <- raw_link$pathway_raw
    raw_link$pathway_raw <- tmp
}

term2gene <- raw_link %>% dplyr::transmute(term = sub("^path:", "", pathway_raw), gene = sub("^mmu:", 
    "", gene_raw)) %>% dplyr::filter(grepl("^mmu[0-9]+$", term), grepl("^[0-9]+$", gene)) %>% dplyr::distinct()

raw_name <- readr::read_tsv(pathway_list_file, col_names = FALSE, show_col_types = FALSE, progress = FALSE,
  trim_ws = TRUE)

if (ncol(raw_name) < 2) {
  stop("KEGG pathway-list file does not contain two tab-separated columns.")
}

raw_name <- raw_name[, 1:2, drop = FALSE]

colnames(raw_name) <- c("pathway_raw", "name_raw")

term2name <- raw_name %>% dplyr::transmute(term = sub("^path:", "", pathway_raw), name = sub(" - Mus musculus \\(house mouse\\)$", 
    "", name_raw)) %>% dplyr::filter(grepl("^mmu[0-9]+$", term), !is.na(name), name != "") %>% dplyr::distinct()

set.seed(SEED)

gsea_obj <- clusterProfiler::GSEA(geneList = gene_list, exponent = 1, minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE,
  eps = 0, pvalueCutoff = 1, pAdjustMethod = "BH", verbose = TRUE, seed = TRUE, by = "fgsea", TERM2GENE = term2gene,
  TERM2NAME = term2name)

gsea_all <- as.data.frame(gsea_obj)

if (nrow(gsea_all) == 0) {
  stop("GSEA returned no pathways.")
}

entrez_to_symbol <- stats::setNames(rank_df$gene, as.character(rank_df$ENTREZID))

gsea_all$core_enrichment_ENTREZ <- gsea_all$core_enrichment

gsea_all$core_enrichment_SYMBOL <- vapply(strsplit(gsea_all$core_enrichment, "/", fixed = TRUE), function(ids) {
    syms <- unname(entrez_to_symbol[ids])
    syms <- syms[!is.na(syms)]
    paste(unique(syms), collapse = "/")
}, FUN.VALUE = character(1))

gsea_all$direction <- ifelse(gsea_all$NES > 0, "Clec7a_enriched", "Other_myeloid_enriched")

gsea_all <- gsea_all %>% dplyr::arrange(p.adjust, dplyr::desc(abs(NES)))

gsea_sig <- gsea_all %>% dplyr::filter(p.adjust < PADJ_CUTOFF)

gsea_pos <- gsea_sig %>% dplyr::filter(NES > 0) %>% dplyr::arrange(p.adjust, dplyr::desc(NES))

gsea_neg <- gsea_sig %>% dplyr::filter(NES < 0) %>% dplyr::arrange(p.adjust, NES)

lipid_pattern <- paste(c("lipid", "fatty acid", "cholesterol", "glycerolipid", "glycerophospholipid",
  "sphingolipid", "arachidonic acid", "linoleic acid", "alpha-linolenic acid", "ether lipid", "PPAR",
  "biosynthesis of unsaturated fatty acids"), collapse = "|")

gsea_lipid <- gsea_all %>% dplyr::filter(grepl(lipid_pattern, Description, ignore.case = TRUE)) %>% dplyr::arrange(p.adjust, 
    dplyr::desc(abs(NES)))

write.csv(gsea_all, file.path(out_dir, "23_Clec7a_GSEA_KEGG_all_v02.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(gsea_sig, file.path(out_dir, "23_Clec7a_GSEA_KEGG_significant_v02.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

write.csv(gsea_pos, file.path(out_dir, "23_Clec7a_GSEA_KEGG_positive_v02.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(gsea_neg, file.path(out_dir, "23_Clec7a_GSEA_KEGG_negative_v02.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(gsea_lipid, file.path(out_dir, "23_Clec7a_GSEA_KEGG_lipid_related_v02.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

top_pos <- gsea_pos %>% dplyr::slice_head(n = min(TOP_POS, nrow(.)))

top_neg <- gsea_neg %>% dplyr::slice_head(n = min(TOP_NEG, nrow(.)))

plot_df <- dplyr::bind_rows(top_pos, top_neg)

if (nrow(plot_df) > 0) {
    plot_df <- plot_df %>% dplyr::arrange(NES) %>% dplyr::mutate(Description = factor(Description, levels = Description), 
        neglog10_padj = -log10(pmax(p.adjust, 1e-300)))
    p_summary <- ggplot(plot_df, aes(x = NES, y = Description)) + geom_vline(xintercept = 0, linewidth = 0.4) + 
        geom_point(aes(size = setSize, color = neglog10_padj), alpha = 0.9) + scale_color_gradient(low = "#4B9CD3", 
        high = "#C73E1D") + labs(title = "KEGG GSEA: Clec7a+ target state vs all other myeloid states", 
        subtitle = "Full-gene ranking | NES > 0: Clec7a+ enriched; NES < 0: other myeloid enriched", 
        x = "Normalized enrichment score (NES)", y = NULL, size = "Pathway size", color = expression(-log[10](adjusted ~ 
            P))) + theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 11), 
        plot.subtitle = element_text(size = 9))
    ggsave(file.path(out_dir, "23_Clec7a_GSEA_KEGG_summary_dotplot_v02.pdf"), p_summary, width = 10, 
        height = 7.5)
    ggsave(file.path(out_dir, "23_Clec7a_GSEA_KEGG_summary_dotplot_v02.png"), p_summary, width = 10, 
        height = 7.5, dpi = 600, bg = "white")
} else {
}

curve_ids <- c(gsea_pos %>% dplyr::slice_head(n = min(TOP_CURVES_EACH_DIRECTION, nrow(.))) %>% dplyr::pull(ID), 
    gsea_neg %>% dplyr::slice_head(n = min(TOP_CURVES_EACH_DIRECTION, nrow(.))) %>% dplyr::pull(ID))

curve_ids <- unique(curve_ids)

curve_pdf <- file.path(out_dir, "23_Clec7a_GSEA_KEGG_top_curves_v02.pdf")

pdf(curve_pdf, width = 8, height = 6, onefile = TRUE)

n_curve_ok <- 0

if (length(curve_ids) > 0) {
    for (pathway_id in curve_ids) {
        pathway_row <- gsea_all[gsea_all$ID == pathway_id, , drop = FALSE]
        pathway_title <- if (nrow(pathway_row) > 0) {
            paste0(pathway_row$Description[1], " | NES=", sprintf("%.2f", pathway_row$NES[1]), " | FDR=", 
                format(pathway_row$p.adjust[1], digits = 3, scientific = TRUE))
        }
        else {
            pathway_id
        }
        gp <- tryCatch({
            enrichplot::gseaplot2(gsea_obj, geneSetID = pathway_id, title = pathway_title, base_size = 10, 
                pvalue_table = FALSE)
        }, error = function(e) {
            log_msg("gseaplot2 failed for ", pathway_id, ": ", conditionMessage(e))
            NULL
        })
        if (!is.null(gp)) {
            print(gp)
            n_curve_ok <- n_curve_ok + 1
        }
    }
}

if (n_curve_ok == 0) {
    plot.new()
    text(0.5, 0.5, "No significant GSEA enrichment curves available.")
}

dev.off()

if (nrow(gsea_pos) > 0) {
    capture.output(print(gsea_pos %>% dplyr::select(ID, Description, NES, p.adjust, setSize) %>% dplyr::slice_head(n = 20)), 
        file = log_file, append = TRUE)
}

if (nrow(gsea_neg) > 0) {
    capture.output(print(gsea_neg %>% dplyr::select(ID, Description, NES, p.adjust, setSize) %>% dplyr::slice_head(n = 20)), 
        file = log_file, append = TRUE)
}

if (nrow(gsea_lipid) > 0) {
    capture.output(print(gsea_lipid %>% dplyr::select(ID, Description, NES, p.adjust, setSize)), file = log_file, 
        append = TRUE)
}

capture.output(sessionInfo(), file = log_file, append = TRUE)
