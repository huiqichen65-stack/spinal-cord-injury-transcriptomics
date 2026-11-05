source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

target_state_formal <- "Clec7a+ lipid-handling Microglia"
target_state_legacy <- "Clec7a+ lipid-associated Microglia"

expected_cells <- 58533L
expected_target_cells <- 17713L
expected_signature_size <- 1320L
expected_top_tf <- 10L
expected_key_gene_n <- 20L

cor_method <- "spearman"
normalization_scale_factor <- 10000

png_dpi <- 600
heatmap_width_in <- 12
heatmap_height_in <- 6.5

outdir <- "D:/zz/results/90_pySCENIC_top10_TF_network_signature_correlation"

rss_ranking_file <- paste0("D:/zz/results/89_pySCENIC_RSS_target_regulon_ranking/", "05_89_Clec7a_lipid_handling_MG_regulon_RSS_ranking_v01.csv")

deg_up_file <- paste0("D:/zz/results/21_Clec7a_DEG_volcano_myeloid/", "21_Clec7a_vs_all_myeloid_DEG_up_v01.csv")

lipid_manifest_file <- paste0("D:/zz/results/20_AUCell_lipid_metabolism_scoring/", "20_AUCell_lipid_gene_set_manifest_v03.csv")

go_bp_file <- paste0("D:/zz/results/22_Clec7a_GO_KEGG_enrichment/", "22_Clec7a_GO_BP_up_v03.csv")

kegg_file <- paste0("D:/zz/results/22_Clec7a_GO_KEGG_enrichment/", "22_Clec7a_KEGG_up_v03.csv")

source_rds <- paste0("D:/zz/results/20_AUCell_lipid_metabolism_scoring/", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

v04_evidence_file <- file.path(outdir, "10_90_top10_TF_evidence_summary_v04.csv")

v04_final_qc_file <- file.path(outdir, "14_90_FINAL_QC_v04.csv")

v04_corr_target_file <- file.path(outdir, "08_90_Figure6E_TF_signature_correlation_target_cells_v04.csv")

v04_zero_library_file <- file.path(outdir, "18_90_zero_library_cells_excluded_from_correlation_v04.csv")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

output_index_file <- file.path(outdir, "00_90_OUTPUT_INDEX_v05.csv")

input_qc_file <- file.path(outdir, "01_90_input_QC_v05.csv")

top10_tf_file <- file.path(outdir, "02_90_top10_RSS_TFs_v05.csv")

functional_support_file <- file.path(outdir, "06_90_Figure6E_20gene_functional_support_audit_v05.csv")

key_gene_file <- file.path(outdir, "07_90_Figure6E_key_functional_gene_panel_v05.csv")

corr_target_file <- file.path(outdir, "08_90_Figure6E_TF_functional_gene_correlation_target_cells_v05.csv")

corr_all_file <- file.path(outdir, "09_90_Figure6E_TF_functional_gene_correlation_all_myeloid_sensitivity_v05.csv")

evidence_summary_file <- file.path(outdir, "10_90_top10_TF_evidence_summary_v05.csv")

heatmap_pdf <- file.path(outdir, "12_90_Figure6E_TF_functional_gene_correlation_v05.pdf")

heatmap_svg <- file.path(outdir, "12_90_Figure6E_TF_functional_gene_correlation_v05.svg")

heatmap_png <- file.path(outdir, "12_90_Figure6E_TF_functional_gene_correlation_v05.png")

bundle_file <- file.path(outdir, "13_90_Figure6E_20gene_functional_panel_bundle_v05.rds")

final_qc_file <- file.path(outdir, "14_90_FINAL_QC_v05.csv")

report_file <- file.path(outdir, "15_90_REPORT_v05.txt")

zero_library_file <- file.path(outdir, "18_90_zero_library_cells_excluded_from_correlation_v05.csv")

comparison_file <- file.path(outdir, "19_90_Figure6E_v04_vs_v05_TF_correlation_comparison_v05.csv")

go_term_audit_file <- file.path(outdir, "20_90_Figure6E_relevant_GO_term_audit_v05.csv")

kegg_audit_file <- file.path(outdir, "21_90_Figure6E_KEGG_context_audit_v05.csv")

state_defining_genes <- c("Clec7a", "Trem2", "Cst7", "Cx3cr1", "Axl", "C1qa")

lipid_handling_genes <- c("Apoe", "Lpl", "Abca1", "Apoc1", "Apoc4", "Npc2", "Soat1")

lysosomal_activation_genes <- c("Tyrobp", "Csf1r", "Itgam", "Hexb", "Ctsd", "Lamp1", "Lamp2")

key_gene_panel <- data.frame(selection_rank = seq_len(expected_key_gene_n), gene = c(state_defining_genes,
  lipid_handling_genes, lysosomal_activation_genes), functional_category = c(rep("State-defining / injury-associated",
  length(state_defining_genes)), rep("Lipid handling", length(lipid_handling_genes)), rep("Lysosomal / endocytic / activation",
  length(lysosomal_activation_genes))), stringsAsFactors = FALSE)

if (nrow(key_gene_panel) != expected_key_gene_n || anyDuplicated(key_gene_panel$gene)) {
  stop("Predeclared functional panel must contain exactly 20 unique genes.")
}

key_genes <- key_gene_panel$gene

relevant_go_terms <- c("microglial cell activation", "innate immune response", "endocytosis", "receptor-mediated endocytosis",
  "glycerolipid metabolic process", "neutral lipid catabolic process")

make_qc <- qc_row

split_gene_ids <- function(x) {
    if (is.na(x) || !nzchar(x)) {
        return(character(0))
    }
    unique(trimws(unlist(strsplit(x, "/", fixed = TRUE))))
}

spearman_pair_table <- function(tf_expr, key_expr, universe_label) {
    if (ncol(tf_expr) != ncol(key_expr)) {
        stop("TF and key-gene expression matrices have different cell counts.")
    }
    n <- ncol(tf_expr)
    rho <- stats::cor(t(tf_expr), t(key_expr), method = cor_method, use = "pairwise.complete.obs")
    out <- expand.grid(TF = rownames(tf_expr), key_gene = rownames(key_expr), stringsAsFactors = FALSE)
    out$rho <- mapply(function(tf, gene) {
        rho[tf, gene]
    }, out$TF, out$key_gene)
    rr <- pmin(pmax(out$rho, -0.999999999), 0.999999999)
    tstat <- rr * sqrt((n - 2)/pmax(1 - rr^2, .Machine$double.eps))
    out$p_value <- 2 * stats::pt(-abs(tstat), df = n - 2)
    out$FDR <- stats::p.adjust(out$p_value, method = "BH")
    out$n_cells <- n
    out$cell_universe <- universe_label
    out
}

summarize_tf_correlation <- function(corr_table, tf_order) {
    ans <- do.call(rbind, lapply(tf_order, function(tf) {
        z <- corr_table[corr_table$TF == tf, , drop = FALSE]
        data.frame(TF = tf, mean_rho_key_functional_genes = mean(z$rho, na.rm = TRUE), median_rho_key_functional_genes = stats::median(z$rho, 
            na.rm = TRUE), max_positive_rho = max(z$rho, na.rm = TRUE), min_rho = min(z$rho, na.rm = TRUE), 
            n_positive_FDR05 = sum(z$rho > 0 & z$FDR < 0.05, na.rm = TRUE), n_negative_FDR05 = sum(z$rho < 
                0 & z$FDR < 0.05, na.rm = TRUE), n_rho_ge_0.10_FDR05 = sum(z$rho >= 0.1 & z$FDR < 0.05, 
                na.rm = TRUE), n_rho_ge_0.20_FDR05 = sum(z$rho >= 0.2 & z$FDR < 0.05, na.rm = TRUE), 
            stringsAsFactors = FALSE)
    }))
    ans$correlation_rank_v05 <- rank(-ans$mean_rho_key_functional_genes, ties.method = "min")
    ans
}

required_packages <- c("Matrix", "SeuratObject", "ggplot2")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

required_files <- c(rss_ranking_file, deg_up_file, lipid_manifest_file, go_bp_file, kegg_file, source_rds,
  v04_evidence_file, v04_final_qc_file, v04_corr_target_file, v04_zero_library_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

v04_qc <- read.csv(v04_final_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

if (any(v04_qc$status == "FAIL")) {
  stop("Step 90 v04 contains a blocking FAIL item.")
}

deg_up <- read.csv(deg_up_file, stringsAsFactors = FALSE, check.names = FALSE)

required_deg_cols <- c("gene", "p_val_adj", "avg_log2FC", "pct.1", "pct.2")

if (!all(required_deg_cols %in% colnames(deg_up))) {
  stop("Step 21 UP DEG file lacks required columns.")
}

if (nrow(deg_up) != expected_signature_size || anyDuplicated(deg_up$gene)) {
  stop("Frozen signature mismatch: observed ", nrow(deg_up), "; expected ", expected_signature_size,
    ".")
}

if (any(deg_up$p_val_adj >= 0.05, na.rm = TRUE) || any(deg_up$avg_log2FC < 0.25, na.rm = TRUE)) {
  stop("Frozen signature violates p_adj<0.05 / avg_log2FC>=0.25.")
}

if (!all(key_genes %in% deg_up$gene)) {
    missing_panel <- setdiff(key_genes, deg_up$gene)
    stop("Predeclared key functional gene(s) are not in frozen Step 21 UP signature: ", paste(missing_panel, 
        collapse = ", "))
}

rss <- read.csv(rss_ranking_file, stringsAsFactors = FALSE, check.names = FALSE)

required_rss_cols <- c("regulon", "TF", "target_RSS", "RSS_rank")

if (!all(required_rss_cols %in% colnames(rss))) {
  stop("Step 89 RSS ranking lacks required columns.")
}

rss <- rss[order(rss$RSS_rank), , drop = FALSE]

top10 <- head(rss, expected_top_tf)

if (nrow(top10) != expected_top_tf || anyDuplicated(top10$TF)) {
  stop("Could not resolve exactly 10 unique RSS TFs.")
}

if (any(key_genes %in% top10$TF)) {
  stop("A key functional gene overlaps the top 10 TF list.")
}

write.csv(top10, top10_tf_file, row.names = FALSE)

lipid <- read.csv(lipid_manifest_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("gene_set", "gene") %in% colnames(lipid))) {
  stop("Step 20 lipid manifest lacks gene_set/gene columns.")
}

if ("present_in_object" %in% colnames(lipid)) {
  lipid <- lipid[lipid$present_in_object %in% TRUE, , drop = FALSE]
}

go_bp <- read.csv(go_bp_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("Description", "geneID") %in% colnames(go_bp))) {
  stop("Step 22 GO-BP file lacks Description/geneID columns.")
}

kegg <- read.csv(kegg_file, stringsAsFactors = FALSE, check.names = FALSE)

go_selected <- go_bp[go_bp$Description %in% relevant_go_terms, , drop = FALSE]

missing_go_terms <- setdiff(relevant_go_terms, go_selected$Description)

if (length(missing_go_terms) > 0) {
  stop("Expected GO term(s) not found: ", paste(missing_go_terms, collapse = ", "))
}

write.csv(go_selected, go_term_audit_file, row.names = FALSE)

write.csv(kegg, kegg_audit_file, row.names = FALSE)

lipid_support <- lapply(key_genes, function(g) {
    gs <- unique(as.character(lipid$gene_set[lipid$gene == g]))
    gs <- gs[!is.na(gs) & nzchar(gs)]
    list(n = length(gs), sets = paste(sort(gs), collapse = "; "))
})

go_support <- lapply(key_genes, function(g) {
    terms <- character(0)
    for (i in seq_len(nrow(go_selected))) {
        genes_this <- split_gene_ids(as.character(go_selected$geneID[i]))
        if (g %in% genes_this) {
            terms <- c(terms, as.character(go_selected$Description[i]))
        }
    }
    terms <- unique(terms)
    list(n = length(terms), terms = paste(sort(terms), collapse = "; "))
})

deg_match <- deg_up[match(key_genes, deg_up$gene), , drop = FALSE]

functional_support <- key_gene_panel

functional_support$avg_log2FC <- deg_match$avg_log2FC
functional_support$pct.1 <- deg_match$pct.1
functional_support$pct.2 <- deg_match$pct.2
functional_support$delta_pct <- deg_match$pct.1 - deg_match$pct.2
functional_support$p_val_adj <- deg_match$p_val_adj

functional_support$n_lipid_gene_sets <- vapply(lipid_support, function(x) x$n, integer(1))

functional_support$lipid_gene_sets <- vapply(lipid_support, function(x) x$sets, character(1))

functional_support$n_relevant_GO_terms <- vapply(go_support, function(x) x$n, integer(1))

functional_support$relevant_GO_terms <- vapply(go_support, function(x) x$terms, character(1))

functional_support$has_frozen_functional_support <- (functional_support$n_lipid_gene_sets > 0 | functional_support$n_relevant_GO_terms >
  0)

if (!all(functional_support$has_frozen_functional_support)) {
    no_support <- functional_support$gene[!functional_support$has_frozen_functional_support]
    stop("Key functional gene(s) lack both Step 20 lipid-set support and selected Step 22 GO support: ", 
        paste(no_support, collapse = ", "))
}

write.csv(functional_support, functional_support_file, row.names = FALSE)

write.csv(functional_support, key_gene_file, row.names = FALSE)

obj <- readRDS(source_rds)

if (!"RNA" %in% names(obj@assays)) {
  stop("RNA assay not found.")
}

if (!"myeloid_state" %in% colnames(obj@meta.data)) {
  stop("myeloid_state metadata column not found.")
}

counts <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")

if (ncol(counts) != expected_cells) {
  stop("Unexpected full-atlas cell count: ", ncol(counts), "; expected ", expected_cells, ".")
}

all_cells <- colnames(counts)

state_source <- as.character(obj@meta.data[all_cells, "myeloid_state"])

state_formal <- state_source

state_formal[state_formal == target_state_legacy] <- target_state_formal

target_cells_total <- all_cells[state_formal == target_state_formal]

if (length(target_cells_total) != expected_target_cells) {
  stop("target-state cell count mismatch: observed ", length(target_cells_total), "; expected ", expected_target_cells,
    ".")
}

genes_needed <- unique(c(top10$TF, key_genes))

missing_expr_genes <- setdiff(genes_needed, rownames(counts))

if (length(missing_expr_genes) > 0) {
  stop("Selected TF/key functional gene(s) missing from RNA counts: ", paste(missing_expr_genes, collapse = ", "))
}

library_size <- Matrix::colSums(counts)

invalid_library <- (!is.finite(library_size) | library_size <= 0)

invalid_cells <- all_cells[invalid_library]

zero_library_audit <- data.frame(cell = invalid_cells, myeloid_state = state_formal[invalid_library],
  RNA_library_size = as.numeric(library_size[invalid_library]), excluded_from_Step90_correlation = TRUE,
  stringsAsFactors = FALSE)

write.csv(zero_library_audit, zero_library_file, row.names = FALSE)

valid_cells <- all_cells[!invalid_library]

target_cells <- intersect(target_cells_total, valid_cells)

n_invalid_library <- length(invalid_cells)
n_invalid_target <- length(setdiff(target_cells_total, target_cells))

if (n_invalid_library >= 0.01 * expected_cells || n_invalid_target >= 0.01 * expected_target_cells) {
  stop("Zero/invalid RNA-library exclusion reached >=1%; review the source object.")
}

v04_zero <- read.csv(v04_zero_library_file, stringsAsFactors = FALSE, check.names = FALSE)

same_zero_cells_as_v04 <- setequal(as.character(v04_zero$cell), invalid_cells)

small_counts <- as.matrix(counts[genes_needed, valid_cells, drop = FALSE])

norm_expr <- log1p(sweep(small_counts, 2, library_size[valid_cells], "/") * normalization_scale_factor)

rownames(norm_expr) <- genes_needed
colnames(norm_expr) <- valid_cells

rm(counts, small_counts, obj)

gc()

tf_expr_target <- norm_expr[top10$TF, target_cells, drop = FALSE]

key_expr_target <- norm_expr[key_genes, target_cells, drop = FALSE]

corr_target <- spearman_pair_table(tf_expr_target, key_expr_target, "Clec7a+ lipid-handling Microglia cells only")

write.csv(corr_target, corr_target_file, row.names = FALSE)

tf_expr_all <- norm_expr[top10$TF, valid_cells, drop = FALSE]

key_expr_all <- norm_expr[key_genes, valid_cells, drop = FALSE]

corr_all <- spearman_pair_table(tf_expr_all, key_expr_all, "All myeloid cells with positive RNA library size")

write.csv(corr_all, corr_all_file, row.names = FALSE)

corr_summary_v05 <- summarize_tf_correlation(corr_target, top10$TF)

evidence_v04 <- read.csv(v04_evidence_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(top10$TF %in% evidence_v04$TF)) {
  stop("v04 evidence table lacks one or more top 10 TFs.")
}

drop_corr <- grepl("rho|correlation_rank|n_positive|n_negative|n_rho", colnames(evidence_v04), ignore.case = TRUE)

evidence_base <- evidence_v04[, !drop_corr, drop = FALSE]

evidence_v05 <- merge(evidence_base, corr_summary_v05, by = "TF", all.x = TRUE, sort = FALSE)

evidence_v05 <- evidence_v05[match(top10$TF, evidence_v05$TF), , drop = FALSE]

write.csv(evidence_v05, evidence_summary_file, row.names = FALSE)

v04_corr <- read.csv(v04_corr_target_file, stringsAsFactors = FALSE, check.names = FALSE)

v04_summary_raw <- do.call(rbind, lapply(top10$TF, function(tf) {
    z <- v04_corr[v04_corr$TF == tf, , drop = FALSE]
    data.frame(TF = tf, mean_rho_v04 = mean(z$rho, na.rm = TRUE), median_rho_v04 = stats::median(z$rho, 
        na.rm = TRUE), max_positive_rho_v04 = max(z$rho, na.rm = TRUE), stringsAsFactors = FALSE)
}))

v04_summary_raw$correlation_rank_v04 <- rank(-v04_summary_raw$mean_rho_v04, ties.method = "min")

comparison <- merge(v04_summary_raw, corr_summary_v05[, c("TF", "mean_rho_key_functional_genes", "median_rho_key_functional_genes",
  "max_positive_rho", "n_rho_ge_0.10_FDR05", "n_rho_ge_0.20_FDR05", "correlation_rank_v05"), drop = FALSE],
  by = "TF", all = TRUE, sort = FALSE)

comparison <- comparison[match(top10$TF, comparison$TF), , drop = FALSE]

comparison$delta_mean_rho_v05_minus_v04 <- comparison$mean_rho_key_functional_genes - comparison$mean_rho_v04

comparison$rank_shift_v05_minus_v04 <- comparison$correlation_rank_v05 - comparison$correlation_rank_v04

write.csv(comparison, comparison_file, row.names = FALSE)

library(ggplot2)

corr_plot <- corr_target

corr_plot$TF <- factor(corr_plot$TF, levels = rev(top10$TF))

corr_plot$key_gene <- factor(corr_plot$key_gene, levels = key_genes)

rho_limit <- max(abs(corr_plot$rho), na.rm = TRUE)

if (!is.finite(rho_limit) || rho_limit <= 0) {
  rho_limit <- 1
}

p_heat <- ggplot(corr_plot, aes(x = key_gene, y = TF, fill = rho)) + geom_tile(linewidth = 0.25, color = "white") +
  geom_vline(xintercept = c(length(state_defining_genes) + 0.5, length(state_defining_genes) + length(lipid_handling_genes) +
    0.5), linewidth = 0.45, linetype = 2) + scale_fill_gradient2(low = "#3B6FB6", mid = "white",
  high = "#B84A4A", midpoint = 0, limits = c(-rho_limit, rho_limit), name = "Spearman\nrho") + labs(x = "Key Clec7a+ lipid-handling Microglia functional genes",
  y = "Top 10 RSS-specific TFs") + theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 45,
  hjust = 1, vjust = 1), axis.text.y = element_text(face = "italic"), axis.title = element_text(face = "bold"),
  legend.title = element_text(face = "bold"), panel.grid = element_blank())

ggsave(heatmap_pdf, p_heat, width = heatmap_width_in, height = heatmap_height_in, units = "in", device = cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(heatmap_svg, p_heat, width = heatmap_width_in, height = heatmap_height_in, units = "in", device = svglite::svglite)
} else {
    svg(heatmap_svg, width = heatmap_width_in, height = heatmap_height_in)
    print(p_heat)
    dev.off()
}

ggsave(heatmap_png, p_heat, width = heatmap_width_in, height = heatmap_height_in, units = "in", dpi = png_dpi)

bundle <- list(frozen_top10_TFs = top10, frozen_signature = deg_up, final_key_functional_genes = functional_support,
  relevant_GO_terms = go_selected, KEGG_context = kegg, correlation_target_cells = corr_target, correlation_all_myeloid = corr_all,
  evidence_summary_v05 = evidence_v05, comparison_v04_vs_v05 = comparison, parameters = list(target_state = target_state_formal,
    key_functional_gene_n = expected_key_gene_n, state_defining_n = length(state_defining_genes),
    lipid_handling_n = length(lipid_handling_genes), lysosomal_activation_n = length(lysosomal_activation_genes),
    correlation_method = cor_method, normalization_scale_factor = normalization_scale_factor, no_correlation_used_for_gene_selection = TRUE,
    no_automatic_top5_TF_selection = TRUE, Figure6D_frozen_from_v03 = TRUE))

saveRDS(bundle, bundle_file, compress = "gzip")

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("Step 90 v04 QC", "PASS", "No blocking v04 FAIL items.")

qc_list[[length(qc_list) + 1]] <- make_qc("Frozen Step 21 signature", if (nrow(deg_up) == expected_signature_size) "PASS" else "FAIL",
  paste0(nrow(deg_up), " / ", expected_signature_size))

qc_list[[length(qc_list) + 1]] <- make_qc("Final key functional gene count", if (nrow(functional_support) ==
  expected_key_gene_n && !anyDuplicated(functional_support$gene)) "PASS" else "FAIL", paste0(nrow(functional_support),
  " unique genes"))

qc_list[[length(qc_list) + 1]] <- make_qc("All 20 genes in frozen target-state signature", if (all(functional_support$gene %in%
  deg_up$gene)) "PASS" else "FAIL", paste0(sum(functional_support$gene %in% deg_up$gene), " / ", expected_key_gene_n))

qc_list[[length(qc_list) + 1]] <- make_qc("All 20 genes have functional support", if (all(functional_support$has_frozen_functional_support)) "PASS" else "FAIL",
  paste0(sum(functional_support$has_frozen_functional_support), " / ", expected_key_gene_n))

qc_list[[length(qc_list) + 1]] <- make_qc("Functional category balance", if (length(state_defining_genes) ==
  6 && length(lipid_handling_genes) == 7 && length(lysosomal_activation_genes) == 7) "PASS" else "FAIL",
  "6 state-defining + 7 lipid-handling + 7 lysosomal/endocytic/activation")

qc_list[[length(qc_list) + 1]] <- make_qc("No correlation used for gene selection", "PASS", "The 20-gene panel was predeclared before v05 correlation.")

qc_list[[length(qc_list) + 1]] <- make_qc("target-state cells", if (length(target_cells_total) == expected_target_cells) "PASS" else "FAIL",
  paste0(length(target_cells_total), " / ", expected_target_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("RNA zero-library exclusion", if (n_invalid_library == 0) {
  "PASS"
} else if (n_invalid_library < 0.01 * expected_cells && n_invalid_target < 0.01 * expected_target_cells) {
  "WARN"
} else {
  "FAIL"
}, paste0("atlas=", n_invalid_library, "/", expected_cells, "; target=", n_invalid_target, "/", expected_target_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("Zero-library continuity with v04", if (same_zero_cells_as_v04) "PASS" else "WARN",
  paste0("same set=", same_zero_cells_as_v04))

qc_list[[length(qc_list) + 1]] <- make_qc("Primary correlation matrix", if (nrow(corr_target) == expected_top_tf *
  expected_key_gene_n && all(is.finite(corr_target$rho))) "PASS" else "FAIL", paste0(nrow(corr_target),
  " / ", expected_top_tf * expected_key_gene_n, " TF-gene pairs"))

qc_list[[length(qc_list) + 1]] <- make_qc("All-myeloid sensitivity matrix", if (nrow(corr_all) == expected_top_tf *
  expected_key_gene_n && all(is.finite(corr_all$rho))) "PASS" else "FAIL", paste0(nrow(corr_all), " TF-gene pairs"))

qc_list[[length(qc_list) + 1]] <- make_qc("Figure 6E outputs", if (all(file.exists(c(heatmap_pdf, heatmap_svg,
  heatmap_png)))) "PASS" else "FAIL", "PDF + SVG + 600 dpi PNG")

qc_list[[length(qc_list) + 1]] <- make_qc("Figure 6D remains ", "PASS", "v03 Figure 6D/network was not rerun or altered.")

qc_list[[length(qc_list) + 1]] <- make_qc("No automatic top-5 TF declaration", "PASS", "Top 5 will be only after reviewing RSS + Figure 6D + final v05 Figure 6E.")

final_qc <- do.call(rbind, qc_list)

rownames(final_qc) <- NULL

write.csv(final_qc, final_qc_file, row.names = FALSE)

n_pass <- sum(final_qc$status == "PASS")

n_warn <- sum(final_qc$status == "WARN")

n_fail <- sum(final_qc$status == "FAIL")

overall_status <- if (n_fail > 0) {
  "FAIL"
} else if (n_warn > 0) {
  "PASS_WITH_WARNINGS"
} else {
  "PASS"
}

top_corr <- corr_summary_v05[order(corr_summary_v05$correlation_rank_v05), , drop = FALSE]

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_file)

input_qc <- data.frame(item = c("Step89 RSS ranking", "Step21 frozen UP signature", "Step20 lipid manifest",
  "Step22 GO BP", "Step22 KEGG", "Full final myeloid RDS", "Step90 v04 evidence", "Step90 v04 QC",
  "Step90 v04 target correlation", "Step90 v04 zero-library audit"), path = c(rss_ranking_file, deg_up_file,
  lipid_manifest_file, go_bp_file, kegg_file, source_rds, v04_evidence_file, v04_final_qc_file, v04_corr_target_file,
  v04_zero_library_file), exists = file.exists(c(rss_ranking_file, deg_up_file, lipid_manifest_file,
  go_bp_file, kegg_file, source_rds, v04_evidence_file, v04_final_qc_file, v04_corr_target_file, v04_zero_library_file)),
  stringsAsFactors = FALSE)

write.csv(input_qc, input_qc_file, row.names = FALSE)

output_files <- c(input_qc_file, top10_tf_file, functional_support_file, key_gene_file, corr_target_file,
  corr_all_file, evidence_summary_file, heatmap_pdf, heatmap_svg, heatmap_png, bundle_file, final_qc_file,
  report_file, zero_library_file, comparison_file, go_term_audit_file, kegg_audit_file)

output_index <- data.frame(file = basename(output_files), exists = file.exists(output_files), size_MB = ifelse(file.exists(output_files),
  file.info(output_files)$size/1024^2, NA_real_), full_path = normalizePath(output_files, winslash = "/",
  mustWork = FALSE), stringsAsFactors = FALSE)

write.csv(output_index, output_index_file, row.names = FALSE)

cat("WARN                   = ", n_warn, "\n", sep = "")
cat("FAIL                   = ", n_fail, "\n", sep = "")

print(functional_support[, c("selection_rank", "gene", "functional_category", "n_lipid_gene_sets", "n_relevant_GO_terms"),
  drop = FALSE], row.names = FALSE)

print(top_corr[, c("TF", "mean_rho_key_functional_genes", "median_rho_key_functional_genes", "max_positive_rho",
  "n_rho_ge_0.10_FDR05", "correlation_rank_v05"), drop = FALSE], row.names = FALSE)

if (n_fail > 0) {
  stop("Step 90 v05 contains blocking FAIL item(s). Do not proceed.")
}
