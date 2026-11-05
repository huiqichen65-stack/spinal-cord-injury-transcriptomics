source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

target_state_formal <- "Clec7a+ lipid-handling Microglia"
target_state_legacy <- "Clec7a+ lipid-associated Microglia"

top_n_tf <- 10L

key_gene_min <- 8L
key_gene_max <- 12L

network_targets_per_tf <- 5L

cor_method <- "spearman"
normalization_scale_factor <- 10000

png_dpi <- 600
network_width_in <- 10
network_height_in <- 8
heatmap_width_in <- 8
heatmap_height_in <- 6.5

overwrite_completed_output <- FALSE

step89_dir <- project_path("results/89_pySCENIC_RSS_target_regulon_ranking")

rss_ranking_file <- file.path(step89_dir, "05_89_Clec7a_lipid_handling_MG_regulon_RSS_ranking_v01.csv")

step89_qc_file <- file.path(step89_dir, "10_89_pySCENIC_FINAL_QC_v01.csv")

step89_report_file <- file.path(step89_dir, "11_89_pySCENIC_REPORT_v01.txt")

regulon_file <- paste0(project_path("results/87_pySCENIC_cisTarget_pruning/"), "06_87_pySCENIC_regulons_v01.csv")

step87_qc_file <- paste0(project_path("results/87_pySCENIC_cisTarget_pruning/"), "09_87_pySCENIC_FINAL_QC_v01.csv")

adjacency_file <- paste0(project_path("results/86_pySCENIC_GRNBoost2/"), "09_86_pySCENIC_GRNBoost2_adjacencies_v03.tsv.gz")

step21_dir <- project_path("results/21_Clec7a_DEG_volcano_myeloid")

deg_up_file <- file.path(step21_dir, "21_Clec7a_vs_all_myeloid_DEG_up_v01.csv")

deg_all_file <- file.path(step21_dir, "21_Clec7a_vs_all_myeloid_DEG_all_v01.csv")

signature_audit_file <- paste0(project_path("results/47_Bulk_scRNA_Clec7a_signature_GSEA/"), "01_47_Bulk_scRNA_Clec7a_signature_GSEA_v02_signature_audit.csv")

figure2_marker_audit_file <- paste0(project_path("results/29_Figure2B_myeloid_marker_dotplot/"), "29_Figure2B_myeloid_marker_dotplot_v06_selected_marker_audit.csv")

figure2_volcano_label_file <- paste0(project_path("results/32_Figure2E_Clec7a_DEG_volcano/"), "32_Figure2E_Clec7a_DEG_volcano_v08_labeled_genes.csv")

source_rds <- paste0(project_path("results/20_AUCell_lipid_metabolism_scoring/"), "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

python_exe <- PYSCENIC_PYTHON

outdir <- project_path("results/90_pySCENIC_top10_TF_network_signature_correlation")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

helper_dir <- file.path(outdir, "_helpers_v03")

if (!dir.exists(helper_dir)) {
  dir.create(helper_dir, recursive = TRUE)
}

output_index_file <- file.path(outdir, "00_90_OUTPUT_INDEX_v03.csv")

input_qc_file <- file.path(outdir, "01_90_input_QC_v03.csv")

top10_tf_file <- file.path(outdir, "02_90_top10_RSS_TFs_v03.csv")

regulon_targets_file <- file.path(outdir, "03_90_top10_regulon_targets_full_v03.csv")

network_summary_file <- file.path(outdir, "04_90_top10_regulon_network_summary_v03.csv")

network_edges_file <- file.path(outdir, "05_90_Figure6D_network_edges_v03.csv")

network_nodes_file <- file.path(outdir, "06_90_Figure6D_network_nodes_v03.csv")

key_gene_audit_file <- file.path(outdir, "07_90_Figure6E_key_signature_gene_selection_v03.csv")

corr_target_file <- file.path(outdir, "08_90_Figure6E_TF_signature_correlation_target_cells_v03.csv")

corr_all_file <- file.path(outdir, "09_90_Figure6E_TF_signature_correlation_all_myeloid_sensitivity_v03.csv")

evidence_summary_file <- file.path(outdir, "10_90_top10_TF_evidence_summary_v03.csv")

network_pdf <- file.path(outdir, "11_90_Figure6D_regulatory_network_v03.pdf")

network_svg <- file.path(outdir, "11_90_Figure6D_regulatory_network_v03.svg")

network_png <- file.path(outdir, "11_90_Figure6D_regulatory_network_v03.png")

heatmap_pdf <- file.path(outdir, "12_90_Figure6E_TF_signature_correlation_v03.pdf")

heatmap_svg <- file.path(outdir, "12_90_Figure6E_TF_signature_correlation_v03.svg")

heatmap_png <- file.path(outdir, "12_90_Figure6E_TF_signature_correlation_v03.png")

bundle_file <- file.path(outdir, "13_90_top10_TF_network_correlation_bundle_v03.rds")

final_qc_file <- file.path(outdir, "14_90_FINAL_QC_v03.csv")

report_file <- file.path(outdir, "15_90_REPORT_v03.txt")

python_log_file <- file.path(outdir, "16_90_regulon_target_extraction_v03.log")

grn_mapping_audit_file <- file.path(outdir, "17_90_GRNBoost2_mapping_audit_v03.csv")

zero_library_audit_file <- file.path(outdir, "18_90_zero_library_cells_excluded_from_correlation_v03.csv")

extract_targets_py <- file.path(helper_dir, "90_extract_top10_regulon_targets_v03.py")

top10_tf_txt <- file.path(helper_dir, "90_top10_TFs_v03.txt")

expected_signature_size <- 1320L
expected_regulons <- 188L
expected_cells <- 58533L
expected_target_cells <- 17713L

q <- shell_quote

make_qc <- qc_row

extract_gene_column <- function(x) {
    nms <- colnames(x)
    exact <- c("gene", "Gene", "genes", "Genes", "symbol", "Symbol", "gene_symbol", "GeneSymbol", "feature", 
        "marker")
    hit <- exact[exact %in% nms]
    if (length(hit) > 0) {
        return(hit[1])
    }
    gene_like <- grep("gene|symbol|marker|feature", nms, ignore.case = TRUE, value = TRUE)
    if (length(gene_like) > 0) {
        return(gene_like[1])
    }
    NA_character_
}

read_optional_key_genes <- function(path, source_label, target_patterns = c("Clec7a", "lipid-associated",
  "lipid-handling")) {

  if (!file.exists(path)) {
    return(data.frame(gene = character(0), source = character(0), stringsAsFactors = FALSE))
  }

  x <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)

  if (is.null(x) || nrow(x) == 0) {
    return(data.frame(gene = character(0), source = character(0), stringsAsFactors = FALSE))
  }

  gene_col <- extract_gene_column(x)

  if (is.na(gene_col)) {
    return(data.frame(gene = character(0), source = character(0), stringsAsFactors = FALSE))
  }

  group_cols <- grep("state|cluster|cell.?type|group|identity|population|subpopulation", colnames(x),
    ignore.case = TRUE, value = TRUE)

  if (length(group_cols) > 0) {

    keep <- rep(FALSE, nrow(x))

    target_found <- FALSE

    for (gc in group_cols) {

      vals <- as.character(x[[gc]])

      this_keep <- rep(FALSE, length(vals))

      for (pat in target_patterns) {
        this_keep <- this_keep | grepl(pat, vals, ignore.case = TRUE)
      }

      if (any(this_keep, na.rm = TRUE)) {
        target_found <- TRUE
        keep <- keep | this_keep
      }
    }

    if (target_found) {
      x <- x[keep, , drop = FALSE]
    }
  }

  genes <- unique(trimws(as.character(x[[gene_col]])))

  genes <- genes[!is.na(genes) & nzchar(genes)]

  data.frame(gene = genes, source = source_label, stringsAsFactors = FALSE)
}

scale01 <- function(x) {
    x <- as.numeric(x)
    if (length(x) == 0 || all(!is.finite(x))) {
        return(rep(0, length(x)))
    }
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) == 0) {
        return(rep(0.5, length(x)))
    }
    (x - rng[1])/diff(rng)
}

spearman_pair_table <- function(tf_expr, key_expr, cells_used, universe_label) {

  if (ncol(tf_expr) != ncol(key_expr)) {
    stop("TF/key expression matrices have different cell counts.")
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

stopifnot(file.copy("scripts/helpers/90_extract_regulon_targets.py", extract_targets_py, overwrite = TRUE))

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("R package 'Matrix' is required.")
}

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("R package 'SeuratObject' is required.")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("R package 'ggplot2' is required.")
}

required_files <- c(python_exe, rss_ranking_file, step89_qc_file, step89_report_file, regulon_file, step87_qc_file,
  adjacency_file, deg_up_file, deg_all_file, signature_audit_file, source_rds)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

step89_qc <- read.csv(step89_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

step87_qc <- read.csv(step87_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

if (any(step89_qc$status == "FAIL")) {
  stop("Step 89 contains blocking FAIL item(s).")
}

if (any(step87_qc$status == "FAIL")) {
  stop("Step 87 contains blocking FAIL item(s).")
}

deg_up <- read.csv(deg_up_file, stringsAsFactors = FALSE, check.names = FALSE)

required_deg_cols <- c("gene", "p_val_adj", "avg_log2FC", "pct.1", "pct.2")

if (!all(required_deg_cols %in% colnames(deg_up))) {
  stop("Step 21 UP DEG file lacks required columns.")
}

if (nrow(deg_up) != expected_signature_size || anyDuplicated(deg_up$gene)) {
  stop("Frozen signature size/uniqueness mismatch. Observed ", nrow(deg_up), "; expected ", expected_signature_size,
    ".")
}

if (any(deg_up$p_val_adj >= 0.05, na.rm = TRUE) || any(deg_up$avg_log2FC < 0.25, na.rm = TRUE)) {
  stop("Step 21 UP signature violates frozen thresholds.")
}

sig_audit <- read.csv(signature_audit_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("metric", "value") %in% colnames(sig_audit))) {
  stop("Step 47 signature audit schema is unexpected.")
}

audit_map <- setNames(as.character(sig_audit$value), sig_audit$metric)

audit_size <- suppressWarnings(as.integer(audit_map["Frozen Step21 Clec7a-UP signature size"]))

if (length(audit_size) != 1 || is.na(audit_size) || audit_size != expected_signature_size) {
  stop("Step 47 audit does not confirm the frozen 1,320-gene signature.")
}

rss <- read.csv(rss_ranking_file, stringsAsFactors = FALSE, check.names = FALSE)

required_rss_cols <- c("regulon", "TF", "target_RSS", "RSS_rank")

if (!all(required_rss_cols %in% colnames(rss))) {
  stop("Step 89 RSS ranking lacks required columns.")
}

if (nrow(rss) != expected_regulons) {
  stop("Expected ", expected_regulons, " regulons in Step 89 ranking; observed ", nrow(rss), ".")
}

rss <- rss[order(rss$RSS_rank), , drop = FALSE]

top10 <- head(rss, top_n_tf)

if (nrow(top10) != top_n_tf || anyDuplicated(top10$TF)) {
  stop("Could not resolve exactly 10 unique top-RSS TFs.")
}

write.csv(top10, top10_tf_file, row.names = FALSE)

writeLines(top10$TF, top10_tf_txt, useBytes = TRUE)

if (file.exists(regulon_targets_file)) {
  unlink(regulon_targets_file)
}

if (file.exists(python_log_file)) {
  unlink(python_log_file)
}

extract_status <- suppressWarnings(system2(command = python_exe, args = c(q(extract_targets_py), q(regulon_file),
  q(top10_tf_txt), q(regulon_targets_file)), stdout = python_log_file, stderr = python_log_file))

if (is.null(extract_status)) {
  extract_status <- 0L
}

if (extract_status != 0 || !file.exists(regulon_targets_file)) {
  stop("Failed to extract motif-pruned targets.\n", "Review: ", python_log_file)
}

reg_targets <- read.csv(regulon_targets_file, stringsAsFactors = FALSE, check.names = FALSE)

if (length(unique(reg_targets$TF)) != top_n_tf) {
  stop("Not all top 10 TF regulons were extracted.")
}

if (!"regulon_weight" %in% colnames(reg_targets) || any(!is.finite(reg_targets$regulon_weight))) {
  stop("Final Step 87 regulon weights were not extracted correctly.")
}

if (requireNamespace("data.table", quietly = TRUE)) {
    adj <- data.table::fread(adjacency_file, data.table = FALSE, showProgress = TRUE)
} else {
    adj <- read.delim(gzfile(adjacency_file), stringsAsFactors = FALSE, check.names = FALSE)
}

required_adj_cols <- c("TF", "target", "importance")

if (!all(required_adj_cols %in% colnames(adj))) {
  stop("adjacency lacks TF/target/importance columns.")
}

adj_top <- adj[adj$TF %in% top10$TF, required_adj_cols, drop = FALSE]

adj_top <- stats::aggregate(importance ~ TF + target, data = adj_top, FUN = max, na.rm = TRUE)

colnames(adj_top)[colnames(adj_top) == "importance"] <- "GRNBoost2_importance"

rm(adj)
gc()

reg_targets <- merge(reg_targets, adj_top, by = c("TF", "target"), all.x = TRUE, sort = FALSE)

reg_targets$GRNBoost2_mapped <- !is.na(reg_targets$GRNBoost2_importance)

grn_mapping_audit <- do.call(rbind, lapply(top10$TF, function(tf) {
    z <- reg_targets[reg_targets$TF == tf, , drop = FALSE]
    data.frame(TF = tf, n_final_regulon_targets = nrow(z), n_exactly_mapped_to_Step86 = sum(z$GRNBoost2_mapped), 
        n_not_exactly_mapped_to_Step86 = sum(!z$GRNBoost2_mapped), mapping_fraction = mean(z$GRNBoost2_mapped), 
        n_unmapped_self_entries = sum(!z$GRNBoost2_mapped & z$target == z$TF), stringsAsFactors = FALSE)
}))

write.csv(grn_mapping_audit, grn_mapping_audit_file, row.names = FALSE)

n_unmapped_grn <- sum(!reg_targets$GRNBoost2_mapped)

deg_all <- read.csv(deg_all_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("gene", "avg_log2FC", "p_val_adj") %in% colnames(deg_all))) {
  stop("Step 21 all-DEG table lacks required columns.")
}

fc_map <- setNames(deg_all$avg_log2FC, deg_all$gene)

padj_map <- setNames(deg_all$p_val_adj, deg_all$gene)

up_signature <- unique(as.character(deg_up$gene))

reg_targets$target_log2FC <- unname(fc_map[reg_targets$target])

reg_targets$target_p_val_adj <- unname(padj_map[reg_targets$target])

reg_targets$target_in_frozen_UP_signature <- reg_targets$target %in% up_signature

write.csv(reg_targets, regulon_targets_file, row.names = FALSE)

network_summary <- do.call(rbind, lapply(top10$TF, function(tf) {
    z <- reg_targets[reg_targets$TF == tf, , drop = FALSE]
    zsig <- z[z$target_in_frozen_UP_signature, , drop = FALSE]
    data.frame(TF = tf, regulon = unique(z$regulon)[1], regulon_size = unique(z$regulon_size)[1], n_signature_targets = nrow(zsig), 
        fraction_targets_in_signature = nrow(zsig)/nrow(z), median_regulon_weight = stats::median(z$regulon_weight, 
            na.rm = TRUE), max_regulon_weight = max(z$regulon_weight, na.rm = TRUE), Step86_mapping_fraction = mean(z$GRNBoost2_mapped), 
        median_GRNBoost2_importance_mapped = if (any(z$GRNBoost2_mapped)) {
            stats::median(z$GRNBoost2_importance[z$GRNBoost2_mapped], na.rm = TRUE)
        }
        else {
            NA_real_
        }, stringsAsFactors = FALSE)
}))

network_summary <- merge(top10[, c("TF", "regulon", "target_RSS", "RSS_rank"), drop = FALSE], network_summary,
  by = c("TF", "regulon"), all.x = TRUE, sort = FALSE)

network_summary <- network_summary[match(top10$TF, network_summary$TF), , drop = FALSE]

write.csv(network_summary, network_summary_file, row.names = FALSE)

network_edges <- do.call(rbind, lapply(top10$TF, function(tf) {

  z <- reg_targets[reg_targets$TF == tf & reg_targets$target_in_frozen_UP_signature, , drop = FALSE]

  if (nrow(z) == 0) {
    z <- reg_targets[reg_targets$TF == tf, , drop = FALSE]
  }

  z <- z[order(-z$regulon_weight), , drop = FALSE]

  head(z, network_targets_per_tf)
}))

network_edges <- unique(network_edges[, c("TF", "target", "regulon_weight", "GRNBoost2_importance", "GRNBoost2_mapped",
  "target_log2FC", "target_p_val_adj", "target_in_frozen_UP_signature"), drop = FALSE])

write.csv(network_edges, network_edges_file, row.names = FALSE)

tf_nodes <- data.frame(name = top10$TF, node_type = "TF", regulon_size = network_summary$regulon_size[match(top10$TF,
  network_summary$TF)], log2FC = unname(fc_map[top10$TF]), stringsAsFactors = FALSE)

target_names <- setdiff(unique(network_edges$target), top10$TF)

target_nodes <- data.frame(name = target_names, node_type = "Target", regulon_size = NA_real_, log2FC = unname(fc_map[target_names]),
  stringsAsFactors = FALSE)

network_nodes <- rbind(tf_nodes, target_nodes)

write.csv(network_nodes, network_nodes_file, row.names = FALSE)

marker_candidates <- read_optional_key_genes(figure2_marker_audit_file, "Figure2B selected marker")

volcano_candidates <- read_optional_key_genes(figure2_volcano_label_file, "Figure2E volcano label")

prior_candidates <- rbind(marker_candidates, volcano_candidates)

if (nrow(prior_candidates) > 0) {
    prior_candidates <- prior_candidates[prior_candidates$gene %in% up_signature & !prior_candidates$gene %in% 
        top10$TF, , drop = FALSE]
    prior_candidates <- prior_candidates[!duplicated(prior_candidates$gene), , drop = FALSE]
}

deg_up$delta_pct <- deg_up$pct.1 - deg_up$pct.2

deg_up$marker_strength_score <- deg_up$avg_log2FC * pmax(deg_up$delta_pct, 0) * sqrt(pmax(deg_up$pct.1,
  0))

deg_up$key_source <- NA_character_

if (nrow(prior_candidates) > 0) {
    source_map <- setNames(prior_candidates$source, prior_candidates$gene)
    deg_up$key_source[deg_up$gene %in% names(source_map)] <- unname(source_map[deg_up$gene[deg_up$gene %in% 
        names(source_map)]])
}

prior_pool <- deg_up[!is.na(deg_up$key_source) & !deg_up$gene %in% top10$TF, , drop = FALSE]

prior_pool <- prior_pool[order(-prior_pool$marker_strength_score, -prior_pool$avg_log2FC), , drop = FALSE]

if (nrow(prior_pool) >= key_gene_max) {
    selected_key <- head(prior_pool, key_gene_max)
} else {
    selected_key <- prior_pool
    n_need <- max(key_gene_min - nrow(selected_key), 0)
    if (n_need > 0) {
        supplement_pool <- deg_up[!deg_up$gene %in% c(selected_key$gene, top10$TF) & deg_up$pct.1 >= 
            0.2 & deg_up$delta_pct >= 0.1, , drop = FALSE]
        supplement_pool <- supplement_pool[order(-supplement_pool$marker_strength_score, -supplement_pool$avg_log2FC), 
            , drop = FALSE]
        supplement <- head(supplement_pool, n_need)
        supplement$key_source <- "Step21 automatic robust-marker supplement"
        selected_key <- rbind(selected_key, supplement)
    }
}

if (nrow(selected_key) < key_gene_min) {
  stop("Could not identify at least ", key_gene_min, " key signature genes.")
}

if (nrow(selected_key) > key_gene_max) {
  selected_key <- head(selected_key, key_gene_max)
}

selected_key$selection_rank <- seq_len(nrow(selected_key))

key_gene_audit <- selected_key[, c("selection_rank", "gene", "key_source", "avg_log2FC", "pct.1", "pct.2",
  "delta_pct", "p_val_adj", "marker_strength_score"), drop = FALSE]

write.csv(key_gene_audit, key_gene_audit_file, row.names = FALSE)

key_genes <- key_gene_audit$gene

obj <- readRDS(source_rds)

if (!"RNA" %in% names(obj@assays)) {
  stop("RNA assay not found.")
}

if (!"myeloid_state" %in% colnames(obj@meta.data)) {
  stop("myeloid_state metadata column not found.")
}

counts <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")

if (ncol(counts) != expected_cells) {
  stop("Unexpected full-atlas cell count: ", ncol(counts))
}

all_cells <- colnames(counts)

state_source <- as.character(obj@meta.data[all_cells, "myeloid_state"])

state_formal <- state_source

state_formal[state_formal == target_state_legacy] <- target_state_formal

target_cells <- all_cells[state_formal == target_state_formal]

if (length(target_cells) != expected_target_cells) {
  stop("Target-cell count mismatch: observed ", length(target_cells), ", expected ", expected_target_cells,
    ".")
}

genes_needed <- unique(c(top10$TF, key_genes))

missing_expr_genes <- setdiff(genes_needed, rownames(counts))

if (length(missing_expr_genes) > 0) {
  stop("Selected TF/key genes missing from RNA counts:\n", paste(missing_expr_genes, collapse = ", "))
}

library_size <- Matrix::colSums(counts)

invalid_library <- (!is.finite(library_size) | library_size <= 0)

n_invalid_library <- sum(invalid_library)

invalid_cells <- names(library_size)[invalid_library]

if (is.null(invalid_cells)) {
  invalid_cells <- all_cells[invalid_library]
}

invalid_state <- state_formal[match(invalid_cells, all_cells)]

zero_library_audit <- data.frame(cell = invalid_cells, myeloid_state = invalid_state, RNA_library_size = as.numeric(library_size[invalid_library]),
  excluded_from_Step90_correlation = TRUE, stringsAsFactors = FALSE)

write.csv(zero_library_audit, zero_library_audit_file, row.names = FALSE)

valid_cells <- all_cells[!invalid_library]

target_cells_total <- target_cells

target_cells <- intersect(target_cells_total, valid_cells)

n_invalid_target <- length(setdiff(target_cells_total, target_cells))

if (length(valid_cells) < 0.99 * expected_cells) {
  stop("More than 1% of atlas cells have invalid/zero RNA library size: ", n_invalid_library, " / ",
    expected_cells, ". Review the source object before correlation analysis.")
}

if (length(target_cells) < 0.99 * expected_target_cells) {
  stop("More than 1% of target-state cells have invalid/zero RNA library size: ", n_invalid_target,
    " / ", expected_target_cells, ". Review the source object before correlation analysis.")
}

small_counts <- as.matrix(counts[genes_needed, valid_cells, drop = FALSE])

norm_expr <- log1p(sweep(small_counts, 2, library_size[valid_cells], "/") * normalization_scale_factor)

rownames(norm_expr) <- genes_needed
colnames(norm_expr) <- valid_cells

rm(counts, small_counts, obj)

gc()

tf_expr_target <- norm_expr[top10$TF, target_cells, drop = FALSE]

key_expr_target <- norm_expr[key_genes, target_cells, drop = FALSE]

corr_target <- spearman_pair_table(tf_expr_target, key_expr_target, target_cells, "Clec7a+ lipid-handling Microglia cells only")

write.csv(corr_target, corr_target_file, row.names = FALSE)

tf_expr_all <- norm_expr[top10$TF, valid_cells, drop = FALSE]

key_expr_all <- norm_expr[key_genes, valid_cells, drop = FALSE]

corr_all <- spearman_pair_table(tf_expr_all, key_expr_all, valid_cells, "All myeloid cells with positive RNA library size")

write.csv(corr_all, corr_all_file, row.names = FALSE)

corr_summary <- do.call(rbind, lapply(top10$TF, function(tf) {
    z <- corr_target[corr_target$TF == tf, , drop = FALSE]
    zsig <- z[z$FDR < 0.05, , drop = FALSE]
    data.frame(TF = tf, mean_rho_key_signature = mean(z$rho, na.rm = TRUE), median_rho_key_signature = stats::median(z$rho, 
        na.rm = TRUE), max_positive_rho = max(z$rho, na.rm = TRUE), n_positive_FDR05 = sum(z$rho > 0 & 
        z$FDR < 0.05, na.rm = TRUE), n_rho_ge_0.20_FDR05 = sum(z$rho >= 0.2 & z$FDR < 0.05, na.rm = TRUE), 
        stringsAsFactors = FALSE)
}))

evidence <- merge(top10, network_summary[, c("TF", "regulon_size", "n_signature_targets", "fraction_targets_in_signature",
  "median_regulon_weight", "max_regulon_weight", "Step86_mapping_fraction", "median_GRNBoost2_importance_mapped"),
  drop = FALSE], by = "TF", all.x = TRUE, sort = FALSE)

evidence <- merge(evidence, corr_summary, by = "TF", all.x = TRUE, sort = FALSE)

evidence <- evidence[match(top10$TF, evidence$TF), , drop = FALSE]

evidence$regulon_size_rank <- rank(-evidence$regulon_size, ties.method = "min")

evidence$signature_target_fraction_rank <- rank(-evidence$fraction_targets_in_signature, ties.method = "min")

evidence$correlation_rank <- rank(-evidence$mean_rho_key_signature, ties.method = "min")

write.csv(evidence, evidence_summary_file, row.names = FALSE)

library(ggplot2)

corr_target$TF <- factor(corr_target$TF, levels = rev(top10$TF))

corr_target$key_gene <- factor(corr_target$key_gene, levels = key_genes)

corr_target$significance_label <- ifelse(corr_target$FDR < 0.001, "***", ifelse(corr_target$FDR < 0.01,
  "**", ifelse(corr_target$FDR < 0.05, "*", "")))

p_heat <- ggplot(corr_target, aes(x = key_gene, y = TF, fill = rho)) + geom_tile(linewidth = 0.25) +
  geom_text(aes(label = significance_label), size = 3) + scale_fill_gradient2(low = "#3B6FB6", mid = "white",
  high = "#B84A4A", midpoint = 0, limits = c(-max(abs(corr_target$rho), na.rm = TRUE), max(abs(corr_target$rho),
    na.rm = TRUE)), name = "Spearman\nrho") + labs(x = "Key Clec7a+ lipid-handling Microglia signature genes",
  y = "Top 10 RSS-specific TFs") + theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 45,
  hjust = 1, vjust = 1), axis.text.y = element_text(face = "italic"), axis.title = element_text(face = "bold"),
  legend.title = element_text(face = "bold"))

ggsave(heatmap_pdf, p_heat, width = heatmap_width_in, height = heatmap_height_in, units = "in", device = cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(heatmap_svg, p_heat, width = heatmap_width_in, height = heatmap_height_in, units = "in", device = svglite::svglite)
} else {
  svg(heatmap_svg, width = heatmap_width_in, height = heatmap_height_in)
  print(p_heat)
  dev.off()
}

ggsave(heatmap_png, p_heat, width = heatmap_width_in, height = heatmap_height_in, units = "in", dpi = png_dpi)

network_plot_status <- "PASS"
network_plot_detail <- "PDF + SVG + 600 dpi PNG"

if (!requireNamespace("igraph", quietly = TRUE)) {

  network_plot_status <- "WARN"
  network_plot_detail <- paste0("igraph not installed; network tables were exported, ", "but Figure 6D image was not drawn.")

} else {

  library(igraph)

  g <- graph_from_data_frame(d = data.frame(from = network_edges$TF, to = network_edges$target, regulon_weight = network_edges$regulon_weight,
    stringsAsFactors = FALSE), directed = TRUE, vertices = network_nodes)

  set.seed(9001)

  lay <- layout_with_fr(g, niter = 1500)

  v_names <- V(g)$name

  node_info <- network_nodes[match(v_names, network_nodes$name), , drop = FALSE]

  is_tf <- node_info$node_type == "TF"

  tf_size_scaled <- 14 + 16 * scale01(node_info$regulon_size)

  vertex_size <- ifelse(is_tf, tf_size_scaled, 7)

  target_fc <- node_info$log2FC

  fc_abs <- max(abs(target_fc[!is_tf]), na.rm = TRUE)

  if (!is.finite(fc_abs) || fc_abs == 0) {
    fc_abs <- 1
  }

  palette_fun <- grDevices::colorRampPalette(c("#3B6FB6", "white", "#B84A4A"))

  pal <- palette_fun(101)

  fc_scaled <- round(1 + 100 * (pmin(pmax(target_fc, -fc_abs), fc_abs) + fc_abs)/(2 * fc_abs))

  fc_scaled[!is.finite(fc_scaled)] <- 51

  vertex_color <- ifelse(is_tf, "#D9D9D9", pal[pmin(pmax(fc_scaled, 1), 101)])

  vertex_shape <- ifelse(is_tf, "square", "circle")

  vertex_label_size <- ifelse(is_tf, 0.85, 0.55)

  edge_width <- 0.5 + 2.5 * scale01(E(g)$regulon_weight)

  draw_network <- function() {

    plot(g, layout = lay, vertex.shape = vertex_shape, vertex.size = vertex_size, vertex.color = vertex_color,
      vertex.frame.color = "#666666", vertex.label = V(g)$name, vertex.label.cex = vertex_label_size,
      vertex.label.color = "black", edge.width = edge_width, edge.color = grDevices::adjustcolor("grey45",
        alpha.f = 0.55), edge.arrow.size = 0.22, edge.curved = 0.08, main = "Top 10 RSS-specific TF regulatory network")

    legend("topleft", legend = c("TF", "Downstream target"), pch = c(15, 16), pt.cex = c(1.5, 1),
      bty = "n")
  }

  grDevices::cairo_pdf(network_pdf, width = network_width_in, height = network_height_in)

  draw_network()

  dev.off()

  svg(network_svg, width = network_width_in, height = network_height_in)

  draw_network()

  dev.off()

  png(network_png, width = network_width_in, height = network_height_in, units = "in", res = png_dpi)

  draw_network()

  dev.off()
}

bundle <- list(top10_TFs = top10, frozen_signature = deg_up, key_signature_genes = key_gene_audit, regulon_targets = reg_targets,
  network_summary = network_summary, network_edges = network_edges, network_nodes = network_nodes,
  correlation_target_cells = corr_target, correlation_all_myeloid = corr_all, evidence_summary = evidence,
  parameters = list(target_state = target_state_formal, top_n_tf = top_n_tf, key_gene_min = key_gene_min,
    key_gene_max = key_gene_max, network_targets_per_tf = network_targets_per_tf, correlation_method = cor_method,
    normalization_scale_factor = normalization_scale_factor, primary_correlation_universe = "target-state cells with positive RNA library size",
    sensitivity_correlation_universe = "all myeloid cells with positive RNA library size", atlas_cells_total = expected_cells,
    atlas_cells_excluded_zero_library = n_invalid_library, target_cells_total = expected_target_cells,
    target_cells_excluded_zero_library = n_invalid_target))

saveRDS(bundle, bundle_file, compress = "gzip")

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("Step 89 QC", "PASS", "No blocking Step 89 FAIL items.")

qc_list[[length(qc_list) + 1]] <- make_qc("Step 87 QC", "PASS", "No blocking Step 87 FAIL items.")

qc_list[[length(qc_list) + 1]] <- make_qc("Frozen Step 21 signature", if (nrow(deg_up) == expected_signature_size) "PASS" else "FAIL",
  paste0(nrow(deg_up), " genes; p_adj<0.05; avg_log2FC>=0.25"))

qc_list[[length(qc_list) + 1]] <- make_qc("Step 47 signature audit", if (audit_size == expected_signature_size) "PASS" else "FAIL",
  paste0("Frozen signature size=", audit_size))

qc_list[[length(qc_list) + 1]] <- make_qc("Top 10 RSS TFs", if (nrow(top10) == top_n_tf && !anyDuplicated(top10$TF)) "PASS" else "FAIL",
  paste(top10$TF, collapse = ", "))

qc_list[[length(qc_list) + 1]] <- make_qc("Motif-pruned target extraction", if (extract_status == 0 &&
  length(unique(reg_targets$TF)) == top_n_tf) "PASS" else "FAIL", paste0(nrow(reg_targets), " TF-target rows"))

qc_list[[length(qc_list) + 1]] <- make_qc("Final Step 87 regulon weights", if (all(is.finite(reg_targets$regulon_weight))) "PASS" else "FAIL",
  paste0("nonfinite=", sum(!is.finite(reg_targets$regulon_weight))))

qc_list[[length(qc_list) + 1]] <- make_qc("Step 86 exact TF-target remapping audit", if (n_unmapped_grn ==
  0) "PASS" else "WARN", paste0(sum(reg_targets$GRNBoost2_mapped), " / ", nrow(reg_targets), " mapped; unmatched=",
  n_unmapped_grn, ". Final Step 87 regulon retained as authoritative."))

qc_list[[length(qc_list) + 1]] <- make_qc("Key signature genes", if (length(key_genes) >= key_gene_min &&
  length(key_genes) <= key_gene_max && all(key_genes %in% up_signature)) "PASS" else "FAIL", paste(key_genes,
  collapse = ", "))

qc_list[[length(qc_list) + 1]] <- make_qc("target-state cell count", if (length(target_cells_total) ==
  expected_target_cells) "PASS" else "FAIL", paste0(length(target_cells_total), " / ", expected_target_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("RNA zero-library exclusion", if (n_invalid_library == 0) {
  "PASS"
} else if (n_invalid_library < 0.01 * expected_cells && n_invalid_target < 0.01 * expected_target_cells) {
  "WARN"
} else {
  "FAIL"
}, paste0("atlas excluded=", n_invalid_library, "/", expected_cells, "; target-state excluded=", n_invalid_target,
  "/", expected_target_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("Target-state correlation cells", if (length(target_cells) >=
  0.99 * expected_target_cells) "PASS" else "FAIL", paste0(length(target_cells), " positive-library cells used"))

qc_list[[length(qc_list) + 1]] <- make_qc("Primary correlation matrix", if (nrow(corr_target) == top_n_tf *
  length(key_genes) && all(is.finite(corr_target$rho))) "PASS" else "FAIL", paste0(top_n_tf, " TFs x ",
  length(key_genes), " genes"))

qc_list[[length(qc_list) + 1]] <- make_qc("All-myeloid sensitivity correlation", if (nrow(corr_all) ==
  top_n_tf * length(key_genes) && all(is.finite(corr_all$rho))) "PASS" else "FAIL", paste0(top_n_tf,
  " TFs x ", length(key_genes), " genes"))

qc_list[[length(qc_list) + 1]] <- make_qc("Figure 6D network", network_plot_status, network_plot_detail)

qc_list[[length(qc_list) + 1]] <- make_qc("Figure 6E heatmap", if (all(file.exists(c(heatmap_pdf, heatmap_svg,
  heatmap_png)))) "PASS" else "FAIL", "PDF + SVG + 600 dpi PNG")

qc_list[[length(qc_list) + 1]] <- make_qc("Top-10 evidence summary", if (nrow(evidence) == top_n_tf) "PASS" else "FAIL",
  paste0(nrow(evidence), " TFs"))

qc_list[[length(qc_list) + 1]] <- make_qc("No automatic top-5 declaration", "PASS", "Step 90 exports evidence; top 5 will be only after result review.")

qc_list[[length(qc_list) + 1]] <- make_qc("Step 90 bundle", if (file.exists(bundle_file) && file.info(bundle_file)$size >
  0) "PASS" else "FAIL", paste0(round(file.info(bundle_file)$size/1024^2, 3), " MB"))

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

writeLines("Numerical results are provided in the accompanying output tables.", 
    report_file)

input_qc <- data.frame(item = c("Step89 ranking", "Step87 regulons", "Step86 adjacency", "Step21 frozen UP signature",
  "Step21 all DEG", "Step47 signature audit", "Figure2B marker audit optional", "Figure2E volcano labels optional",
  "Full myeloid RDS"), path = c(rss_ranking_file, regulon_file, adjacency_file, deg_up_file, deg_all_file,
  signature_audit_file, figure2_marker_audit_file, figure2_volcano_label_file, source_rds), exists = file.exists(c(rss_ranking_file,
  regulon_file, adjacency_file, deg_up_file, deg_all_file, signature_audit_file, figure2_marker_audit_file,
  figure2_volcano_label_file, source_rds)), stringsAsFactors = FALSE)

write.csv(input_qc, input_qc_file, row.names = FALSE)

output_files <- c(input_qc_file, top10_tf_file, regulon_targets_file, network_summary_file, network_edges_file,
  network_nodes_file, key_gene_audit_file, corr_target_file, corr_all_file, evidence_summary_file,
  network_pdf, network_svg, network_png, heatmap_pdf, heatmap_svg, heatmap_png, bundle_file, final_qc_file,
  report_file, python_log_file, grn_mapping_audit_file, zero_library_audit_file)

output_index <- data.frame(file = basename(output_files), exists = file.exists(output_files), size_MB = ifelse(file.exists(output_files),
  file.info(output_files)$size/1024^2, NA_real_), full_path = normalizePath(output_files, winslash = "/",
  mustWork = FALSE), stringsAsFactors = FALSE)

write.csv(output_index, output_index_file, row.names = FALSE)

cat("WARN             = ", n_warn, "\n", sep = "")
cat("FAIL             = ", n_fail, "\n", sep = "")

if (n_fail > 0) {
  stop("Step 90 contains blocking FAIL item(s). ", "Do not proceed to scTenifoldKnk.")
}
