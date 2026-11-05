source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

P_ADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0.25
SHOW_GO_BP <- 20
SHOW_KEGG <- 20

required_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "ggplot2", "dplyr", "readr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), "\nInstall with:\n",
    "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n", "BiocManager::install(c('clusterProfiler','org.Mm.eg.db','AnnotationDbi'))\n",
    "install.packages(c('ggplot2','dplyr','readr'))")
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

deg_expected <- file.path(project_dir, "results", "21_Clec7a_DEG_volcano_myeloid", "21_Clec7a_vs_all_myeloid_DEG_all_v01.csv")

if (file.exists(deg_expected)) {
  deg_file <- normalizePath(deg_expected, winslash = "/", mustWork = TRUE)
} else {
  hits <- list.files(project_dir, pattern = "^21_Clec7a_vs_all_myeloid_DEG_all_v01\\.csv$", recursive = TRUE,
    full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) {
    stop("Cannot find 21_Clec7a_vs_all_myeloid_DEG_all_v01.csv")
  }
  deg_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

kegg_dir <- file.path(project_dir, "reference", "KEGG")

gene_pathway_expected <- file.path(kegg_dir, "mmu_gene_pathway.txt")

pathway_list_expected <- file.path(kegg_dir, "mmu_pathway_list.txt")

find_local_kegg_file <- function(expected, pattern) {
    if (file.exists(expected)) {
        return(normalizePath(expected, winslash = "/", mustWork = TRUE))
    }
    hits <- list.files(kegg_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        return(NA_character_)
    }
    normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

gene_pathway_file <- find_local_kegg_file(gene_pathway_expected, "mmu.*gene.*pathway|gene.*pathway.*mmu|link.*pathway")

pathway_list_file <- find_local_kegg_file(pathway_list_expected, "mmu.*pathway.*list|pathway.*list.*mmu|list.*pathway")

if (is.na(gene_pathway_file)) {
  stop("Local KEGG gene-pathway file not found under:\n", kegg_dir, "\nExpected name: mmu_gene_pathway.txt")
}

if (is.na(pathway_list_file)) {
  stop("Local KEGG pathway-name file not found under:\n", kegg_dir, "\nExpected name: mmu_pathway_list.txt")
}

out_dir <- file.path(project_dir, "results", "22_Clec7a_GO_KEGG_enrichment")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "22_Clec7a_GO_KEGG_enrichment_v03.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

deg <- readr::read_csv(deg_file, show_col_types = FALSE)

required_cols <- c("gene", "avg_log2FC", "p_val_adj")

missing_cols <- setdiff(required_cols, colnames(deg))

if (length(missing_cols) > 0) {
  stop("Missing DEG columns: ", paste(missing_cols, collapse = ", "))
}

deg <- deg %>% dplyr::filter(!is.na(gene), gene != "") %>% dplyr::distinct(gene, .keep_all = TRUE)

up_deg <- deg %>% dplyr::filter(p_val_adj < P_ADJ_CUTOFF, avg_log2FC >= LOG2FC_CUTOFF)

write.csv(up_deg, file.path(out_dir, "22_Clec7a_input_upregulated_DEGs_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

all_symbols <- unique(deg$gene)

up_symbols <- unique(up_deg$gene)

all_map <- suppressMessages(clusterProfiler::bitr(all_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
    OrgDb = org.Mm.eg.db)) %>% dplyr::distinct(SYMBOL, ENTREZID)

up_map <- suppressMessages(clusterProfiler::bitr(up_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
    OrgDb = org.Mm.eg.db)) %>% dplyr::distinct(SYMBOL, ENTREZID)

universe_entrez <- unique(as.character(all_map$ENTREZID))

up_entrez <- unique(as.character(up_map$ENTREZID))

write.csv(all_map, file.path(out_dir, "22_Clec7a_gene_mapping_background_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

write.csv(up_map, file.path(out_dir, "22_Clec7a_gene_mapping_up_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

run_go <- function(ont) {
    clusterProfiler::enrichGO(gene = up_entrez, universe = universe_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", 
        ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
}

ego_bp <- run_go("BP")
ego_cc <- run_go("CC")
ego_mf <- run_go("MF")
ego_all <- run_go("ALL")

go_bp_df <- as.data.frame(ego_bp)

go_cc_df <- as.data.frame(ego_cc)

go_mf_df <- as.data.frame(ego_mf)

go_all_df <- as.data.frame(ego_all)

write.csv(go_bp_df, file.path(out_dir, "22_Clec7a_GO_BP_up_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(go_cc_df, file.path(out_dir, "22_Clec7a_GO_CC_up_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(go_mf_df, file.path(out_dir, "22_Clec7a_GO_MF_up_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(go_all_df, file.path(out_dir, "22_Clec7a_GO_ALL_up_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

raw_link <- readr::read_tsv(gene_pathway_file, col_names = FALSE, show_col_types = FALSE, progress = FALSE,
  trim_ws = TRUE)

if (ncol(raw_link) < 2) {
  stop("mmu_gene_pathway.txt does not contain at least 2 tab-separated columns.")
}

raw_link <- raw_link[, 1:2, drop = FALSE]

colnames(raw_link) <- c("gene_raw", "pathway_raw")

first_gene_like <- mean(grepl("^mmu:", raw_link$gene_raw))

second_gene_like <- mean(grepl("^mmu:", raw_link$pathway_raw))

if (second_gene_like > first_gene_like) {
    raw_link <- raw_link %>% dplyr::rename(tmp = gene_raw) %>% dplyr::rename(gene_raw = pathway_raw) %>% 
        dplyr::rename(pathway_raw = tmp)
}

term2gene <- raw_link %>% dplyr::transmute(term = sub("^path:", "", pathway_raw), gene = sub("^mmu:", 
    "", gene_raw)) %>% dplyr::filter(grepl("^mmu[0-9]+$", term), grepl("^[0-9]+$", gene)) %>% dplyr::distinct()

raw_name <- readr::read_tsv(pathway_list_file, col_names = FALSE, show_col_types = FALSE, progress = FALSE,
  trim_ws = TRUE)

if (ncol(raw_name) < 2) {
  stop("mmu_pathway_list.txt does not contain at least 2 tab-separated columns.")
}

raw_name <- raw_name[, 1:2, drop = FALSE]

colnames(raw_name) <- c("pathway_raw", "name_raw")

term2name <- raw_name %>% dplyr::transmute(term = sub("^path:", "", pathway_raw), name = sub(" - Mus musculus \\(house mouse\\)$", 
    "", name_raw)) %>% dplyr::filter(grepl("^mmu[0-9]+$", term), !is.na(name), name != "") %>% dplyr::distinct()

write.csv(term2gene, file.path(out_dir, "22_Clec7a_local_KEGG_TERM2GENE_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

write.csv(term2name, file.path(out_dir, "22_Clec7a_local_KEGG_TERM2NAME_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

if (nrow(term2gene) < 1000) {
  stop("Too few local KEGG gene-pathway links were parsed. ", "Inspect the downloaded file format.")
}

ekegg <- clusterProfiler::enricher(gene = up_entrez, universe = universe_entrez, TERM2GENE = term2gene,
  TERM2NAME = term2name, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
  maxGSSize = 500)

kegg_df <- as.data.frame(ekegg)

if (nrow(kegg_df) > 0) {
    entrez_to_symbol <- stats::setNames(all_map$SYMBOL, all_map$ENTREZID)
    kegg_df$geneID_ENTREZ <- kegg_df$geneID
    kegg_df$geneID_SYMBOL <- vapply(strsplit(kegg_df$geneID, "/", fixed = TRUE), function(ids) {
        syms <- unname(entrez_to_symbol[ids])
        syms <- syms[!is.na(syms)]
        paste(unique(syms), collapse = "/")
    }, FUN.VALUE = character(1))
}

write.csv(kegg_df, file.path(out_dir, "22_Clec7a_KEGG_up_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

ratio_to_numeric <- function(x) {
    vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
        if (length(z) != 2) {
            return(NA_real_)
        }
        as.numeric(z[1])/as.numeric(z[2])
    }, FUN.VALUE = numeric(1))
}

make_dotplot <- function(df, title, show_n = 20) {
    if (nrow(df) == 0) {
        return(NULL)
    }
    x <- df %>% dplyr::filter(!is.na(p.adjust)) %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(show_n, 
        nrow(.))) %>% dplyr::mutate(GeneRatio_numeric = ratio_to_numeric(GeneRatio), Description = factor(Description, 
        levels = rev(Description)), neglog10_padj = -log10(pmax(p.adjust, 1e-300)))
    ggplot(x, aes(x = GeneRatio_numeric, y = Description)) + geom_point(aes(size = Count, color = neglog10_padj), 
        alpha = 0.9) + scale_color_gradient(low = "#4B9CD3", high = "#C73E1D") + labs(title = title, 
        x = "GeneRatio", y = NULL, size = "Count", color = expression(-log[10](adjusted ~ P))) + theme_classic(base_size = 10) + 
        theme(plot.title = element_text(face = "bold", size = 11))
}

make_barplot <- function(df, title, show_n = 20) {
    if (nrow(df) == 0) {
        return(NULL)
    }
    x <- df %>% dplyr::filter(!is.na(p.adjust)) %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(show_n, 
        nrow(.))) %>% dplyr::mutate(Description = factor(Description, levels = rev(Description)), neglog10_padj = -log10(pmax(p.adjust, 
        1e-300)))
    ggplot(x, aes(x = neglog10_padj, y = Description)) + geom_col() + labs(title = title, x = expression(-log[10](adjusted ~ 
        P)), y = NULL) + theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold", 
        size = 11))
}

if (nrow(go_bp_df) > 0) {
    p_go_dot <- make_dotplot(go_bp_df, "GO Biological Process enrichment\nClec7a+ target-state upregulated genes", 
        SHOW_GO_BP)
    p_go_bar <- make_barplot(go_bp_df, "GO Biological Process enrichment\nClec7a+ target-state upregulated genes", 
        SHOW_GO_BP)
    ggsave(file.path(out_dir, "22_Clec7a_GO_BP_dotplot_v03.pdf"), p_go_dot, width = 9, height = 7)
    ggsave(file.path(out_dir, "22_Clec7a_GO_BP_barplot_v03.pdf"), p_go_bar, width = 9, height = 7)
    ggsave(file.path(out_dir, "22_Clec7a_GO_BP_dotplot_v03.png"), p_go_dot, width = 9, height = 7, dpi = 600, 
        bg = "white")
}

if (nrow(kegg_df) > 0) {
    p_kegg_dot <- make_dotplot(kegg_df, "KEGG enrichment\nClec7a+ target-state upregulated genes", SHOW_KEGG)
    p_kegg_bar <- make_barplot(kegg_df, "KEGG enrichment\nClec7a+ target-state upregulated genes", SHOW_KEGG)
    ggsave(file.path(out_dir, "22_Clec7a_KEGG_dotplot_v03.pdf"), p_kegg_dot, width = 9, height = 7)
    ggsave(file.path(out_dir, "22_Clec7a_KEGG_barplot_v03.pdf"), p_kegg_bar, width = 9, height = 7)
    ggsave(file.path(out_dir, "22_Clec7a_KEGG_dotplot_v03.png"), p_kegg_dot, width = 9, height = 7, dpi = 600, 
        bg = "white")
} else {
    log_msg("No significant KEGG terms; KEGG plots skipped.")
}

top_go_bp <- go_bp_df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(30, nrow(.)))

write.csv(top_go_bp, file.path(out_dir, "22_Clec7a_GO_BP_top30_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

if (nrow(kegg_df) > 0) {
    top_kegg <- kegg_df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(30, nrow(.)))
    write.csv(top_kegg, file.path(out_dir, "22_Clec7a_KEGG_top30_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")
} else {
    top_kegg <- data.frame()
}

capture.output(print(top_go_bp %>%
  dplyr::select(ID, Description, GeneRatio, p.adjust, Count)), file = log_file, append = TRUE)

if (nrow(top_kegg) > 0) {
    capture.output(print(top_kegg %>% dplyr::select(ID, Description, GeneRatio, p.adjust, Count)), file = log_file, 
        append = TRUE)
}

capture.output(sessionInfo(), file = log_file, append = TRUE)
