rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(20261118)

ROOT <- PROJECT_ROOT
OUT_DIR <- file.path( ROOT, "results", "147_Reviewer1_RemainingMajor4_bulk_scale_revision" )

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "147_Reviewer1_RemainingMajor4_bulk_scale_revision_v01"

FDR_CUTOFF <- 0.05
MIN_GS_SIZE <- 10L
MAX_GS_SIZE <- 500L

STEP46_DIR <- file.path( ROOT, "results", "46_Bulk_limma_temporal_DEG_analysis" )

STEP48_DIR <- file.path( ROOT, "results", "48_Bulk_functional_enrichment" )

PRIMARY_FILE <- file.path( STEP46_DIR, "04_46_Bulk_limma_temporal_DEG_analysis_v06_PRIMARY_EstablishedLate_vs_Early.csv" )

STEP48_MAPPING_FILE <- file.path( STEP48_DIR, "09_48_Bulk_functional_enrichment_v04_symbol_ENTREZ_mapping.csv" )

SUPP_CANDIDATES <- c( file.path(ROOT, "supplementary", "Supplementary_Tables.xlsx"), file.path(ROOT, "Molecular Neurobiology v02", "Supplementary_Tables.xlsx"),
  file.path(ROOT, "Molecular Neurobiology v01", "Supplementary_Tables.xlsx") )

SUPP_XLSX <- SUPP_CANDIDATES[file.exists(SUPP_CANDIDATES)]
if (length(SUPP_XLSX) > 0L) {
  SUPP_XLSX <- SUPP_XLSX[[1]]
} else {
  SUPP_XLSX <- NA_character_
}

TIMEPOINT_FILES <- data.frame( contrast = c( "15min_vs_Uninjured", "1dpi_vs_Uninjured", "3dpi_vs_Uninjured", "7dpi_vs_Uninjured", "14dpi_vs_Uninjured", "28dpi_vs_Uninjured",
  "42dpi_vs_Uninjured" ), display_label = c( "15 min", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "42 dpi" ), file = file.path( STEP46_DIR, c(
  "21_46_Bulk_limma_temporal_DEG_analysis_v06_15min_vs_Uninjured.csv", "22_46_Bulk_limma_temporal_DEG_analysis_v06_1dpi_vs_Uninjured.csv",
  "23_46_Bulk_limma_temporal_DEG_analysis_v06_3dpi_vs_Uninjured.csv", "24_46_Bulk_limma_temporal_DEG_analysis_v06_7dpi_vs_Uninjured.csv",
  "25_46_Bulk_limma_temporal_DEG_analysis_v06_14dpi_vs_Uninjured.csv", "26_46_Bulk_limma_temporal_DEG_analysis_v06_28dpi_vs_Uninjured.csv",
  "27_46_Bulk_limma_temporal_DEG_analysis_v06_42dpi_vs_Uninjured.csv" ) ), stringsAsFactors = FALSE )

ADDITIONAL_CONTRAST_FILES <- data.frame( contrast = c( "PRIMARY_EstablishedLate_vs_Early", "Transition7dpi_vs_Early", "Extended42dpi_vs_Early", "Sensitivity_28dpi_vs_14dpi", "28dpi_vs_14dpi"
  ), display_label = c( "Established late vs Early", "7 dpi vs Early", "42 dpi vs Early", "28 dpi vs 14 dpi sensitivity", "28 dpi vs 14 dpi" ), file = file.path( STEP46_DIR, c(
  "04_46_Bulk_limma_temporal_DEG_analysis_v06_PRIMARY_EstablishedLate_vs_Early.csv", "05_46_Bulk_limma_temporal_DEG_analysis_v06_Transition7dpi_vs_Early.csv",
  "06_46_Bulk_limma_temporal_DEG_analysis_v06_Extended42dpi_vs_Early.csv", "07_46_Bulk_limma_temporal_DEG_analysis_v06_Sensitivity_28dpi_vs_14dpi.csv",
  "28_46_Bulk_limma_temporal_DEG_analysis_v06_28dpi_vs_14dpi.csv" ) ), stringsAsFactors = FALSE )

ALL_CONTRAST_FILES <- rbind( ADDITIONAL_CONTRAST_FILES, TIMEPOINT_FILES[, c("contrast", "display_label", "file")] )

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages( pkg, repos = "https://cloud.r-project.org", dependencies = TRUE )
    }
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package unavailable: ", pkg)
  }
}

for (p in c("ggplot2", "ggrepel", "openxlsx")) {
  install_if_missing(p, bioc = FALSE)
}

for (p in c("clusterProfiler", "org.Mm.eg.db", "AnnotationDbi")) {
  install_if_missing(p, bioc = TRUE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

msg <- function(...) {
  cat( sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "" )
}

num_name <- function(num, suffix) {
  paste0( sprintf("%02d", as.integer(num)), "_", PREFIX, "_", suffix )
}

out_path <- function(num, suffix) {
  file.path(OUT_DIR, num_name(num, suffix))
}

write_csv <- function(x, file) {
  utils::write.csv( x, file = file, row.names = FALSE, na = "", fileEncoding = "UTF-8" )
}

write_csv_gz <- function(x, file) {
  con <- gzfile(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, file = con, row.names = FALSE, na = "")
}

assert_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop("Missing required input [", label, "]: ", path)
  }
  invisible(TRUE)
}

assert_no_serialized_input <- function(paths) {
  ext <- tolower(tools::file_ext(paths))
  forbidden <- c("rds", "rdata", "rda", "rdx", "rdb", "feather", "loom", "zip")
  bad <- paths[ext %in% forbidden]
  if (length(bad) > 0L) {
    stop( "Forbidden serialized/archive input detected. Step147 must use CSV/XLSX only:\n", paste(bad, collapse = "\n") )
  }
}

safe_neglog10 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x) | x <= 0] <- .Machine$double.xmin
  -log10(x)
}

REQUIRED_STEP46_COLUMNS <- c( "gene_id", "gene_symbol_restored", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B" )

standardize_step46 <- function(path, contrast, display_label) {
  assert_file(path, contrast)

  tab <- utils::read.csv( path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM" )

  missing_cols <- setdiff(REQUIRED_STEP46_COLUMNS, names(tab))
  if (length(missing_cols) > 0L) {
    stop( "Step46 schema mismatch in ", path, "\nMissing: ", paste(missing_cols, collapse = ", ") )
  }

  tab$logFC <- suppressWarnings(as.numeric(tab$logFC))
  tab$adj.P.Val <- suppressWarnings(as.numeric(tab$adj.P.Val))
  tab$P.Value <- suppressWarnings(as.numeric(tab$P.Value))
  tab$t <- suppressWarnings(as.numeric(tab$t))
  tab$AveExpr <- suppressWarnings(as.numeric(tab$AveExpr))
  tab$B <- suppressWarnings(as.numeric(tab$B))

  drop_cols <- intersect( c("stat_sig", "effect_sig", "direction"), names(tab) )
  if (length(drop_cols) > 0L) {
    tab <- tab[, setdiff(names(tab), drop_cols), drop = FALSE]
  }

  names(tab)[names(tab) == "logFC"] <- "transformed_scale_difference"

  tab$FDR_lt_0.05 <- ( is.finite(tab$adj.P.Val) & tab$adj.P.Val < FDR_CUTOFF )

  tab$scale_direction <- "Not significant"
  tab$scale_direction[ tab$FDR_lt_0.05 & is.finite(tab$transformed_scale_difference) & tab$transformed_scale_difference > 0 ] <- "Higher"

  tab$scale_direction[ tab$FDR_lt_0.05 & is.finite(tab$transformed_scale_difference) & tab$transformed_scale_difference < 0 ] <- "Lower"

  tab$scale_direction <- factor( tab$scale_direction, levels = c("Higher", "Lower", "Not significant") )

  tab$contrast <- contrast
  tab$display_label <- display_label

  front <- c( "contrast", "display_label", "gene_id", "gene_symbol_restored", "transformed_scale_difference", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "FDR_lt_0.05", "scale_direction" )

  tab <- tab[, c(front, setdiff(names(tab), front)), drop = FALSE]

  tab
}

count_summary_one <- function(tab) {
  data.frame( contrast = as.character(tab$contrast[[1]]), display_label = as.character(tab$display_label[[1]]), tested_genes_n = nrow(tab),
    FDR_lt_0.05_higher_n = sum(tab$scale_direction == "Higher", na.rm = TRUE), FDR_lt_0.05_lower_n = sum(tab$scale_direction == "Lower", na.rm = TRUE),
    FDR_lt_0.05_total_n = sum(tab$FDR_lt_0.05, na.rm = TRUE), absolute_effect_threshold_used = FALSE, coefficient_interpretation = "difference on the supplied transformed expression scale",
    stringsAsFactors = FALSE )
}

input_paths <- c( ALL_CONTRAST_FILES$file, STEP48_MAPPING_FILE )

assert_no_serialized_input(input_paths)

for (p in input_paths) assert_file(p)

input_audit <- data.frame( role = c( rep("Step46 complete all-gene contrast CSV", nrow(ALL_CONTRAST_FILES)), "Step48 SYMBOL-to-ENTREZ mapping CSV",
  "Supplementary workbook for terminology audit" ), path = c( ALL_CONTRAST_FILES$file, STEP48_MAPPING_FILE, ifelse(is.na(SUPP_XLSX), "NOT FOUND", SUPP_XLSX) ), exists = c(
  file.exists(ALL_CONTRAST_FILES$file), file.exists(STEP48_MAPPING_FILE), !is.na(SUPP_XLSX) && file.exists(SUPP_XLSX) ), read_RDS = FALSE, read_expression_matrix = FALSE,
  rerun_limma = FALSE, stringsAsFactors = FALSE )

write_csv(input_audit, out_path(1, "INPUT_AUDIT.csv"))

msg("Loading Step46 complete CSV statistics; no model refitting...")

contrast_tables <- vector("list", nrow(ALL_CONTRAST_FILES))
names(contrast_tables) <- ALL_CONTRAST_FILES$contrast

for (i in seq_len(nrow(ALL_CONTRAST_FILES))) {
  contrast_tables[[i]] <- standardize_step46( path = ALL_CONTRAST_FILES$file[[i]], contrast = ALL_CONTRAST_FILES$contrast[[i]], display_label = ALL_CONTRAST_FILES$display_label[[i]] )
}

contrast_count_summary <- do.call( rbind, lapply(contrast_tables, count_summary_one) )
rownames(contrast_count_summary) <- NULL

write_csv( contrast_count_summary, out_path(2, "CONTRAST_FDR_COUNT_SUMMARY.csv") )

primary <- contrast_tables[["PRIMARY_EstablishedLate_vs_Early"]]

primary_sig <- primary[ primary$FDR_lt_0.05, , drop = FALSE ]

write_csv( primary, out_path(3, "PRIMARY_all_genes_transformed_scale.csv") )

write_csv( primary_sig, out_path(4, "PRIMARY_FDR_significant_genes.csv") )

timepoint_list <- lapply( TIMEPOINT_FILES$contrast, function(cn) contrast_tables[[cn]] )
timepoint_all <- do.call(rbind, timepoint_list)
rownames(timepoint_all) <- NULL

write_csv_gz( timepoint_all, out_path(5, "TIMEPOINT_all_genes_transformed_scale.csv.gz") )

timepoint_count_summary <- contrast_count_summary[ contrast_count_summary$contrast %in% TIMEPOINT_FILES$contrast, , drop = FALSE ]

timepoint_count_summary <- timepoint_count_summary[ match(TIMEPOINT_FILES$contrast, timepoint_count_summary$contrast), , drop = FALSE ]

write_csv( timepoint_count_summary, out_path(6, "TIMEPOINT_FDR_COUNT_SUMMARY.csv") )

msg("Preparing threshold-free (FDR-only) ORA foregrounds...")

mapping <- utils::read.csv( STEP48_MAPPING_FILE, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM" )

required_map_cols <- c("gene_symbol", "ENTREZID")
missing_map_cols <- setdiff(required_map_cols, names(mapping))
if (length(missing_map_cols) > 0L) {
  stop( "Step48 mapping schema mismatch. Missing: ", paste(missing_map_cols, collapse = ", ") )
}

mapping$gene_symbol <- as.character(mapping$gene_symbol)
mapping$ENTREZID <- as.character(mapping$ENTREZID)

mapping <- mapping[ !is.na(mapping$gene_symbol) & nzchar(mapping$gene_symbol) & !is.na(mapping$ENTREZID) & nzchar(mapping$ENTREZID), , drop = FALSE ]

primary_symbols <- as.character(primary$gene_symbol_restored)

higher_symbols <- as.character( primary$gene_symbol_restored[ primary$scale_direction == "Higher" ] )

lower_symbols <- as.character( primary$gene_symbol_restored[ primary$scale_direction == "Lower" ] )

background_entrez <- unique( mapping$ENTREZID[mapping$gene_symbol %in% primary_symbols] )
higher_entrez <- unique( mapping$ENTREZID[mapping$gene_symbol %in% higher_symbols] )
lower_entrez <- unique( mapping$ENTREZID[mapping$gene_symbol %in% lower_symbols] )

mapping_qc <- data.frame( metric = c( "Step46_primary_input", "Step48_mapping_input", "tested_gene_n", "FDR_lt_0.05_higher_gene_n", "FDR_lt_0.05_lower_gene_n", "background_unique_ENTREZ_n",
  "higher_unique_ENTREZ_n", "lower_unique_ENTREZ_n", "FDR_cutoff", "absolute_effect_size_cutoff", "coefficient_interpretation", "GSEA_status" ), value = c( PRIMARY_FILE, STEP48_MAPPING_FILE,
  nrow(primary), length(higher_symbols), length(lower_symbols), length(background_entrez), length(higher_entrez), length(lower_entrez), FDR_CUTOFF, "NONE",
  "difference on the supplied transformed expression scale", "Not rerun; existing Step48 GSEA used moderated t ranking and does not depend on |logFC| >= 0.5" ), stringsAsFactors = FALSE )

write_csv( mapping_qc, out_path(7, "ORA_MAPPING_QC.csv") )

entrez_symbol_map <- mapping[ !duplicated(mapping$ENTREZID), c("ENTREZID", "gene_symbol"), drop = FALSE ]
entrez_to_symbol <- setNames( entrez_symbol_map$gene_symbol, entrez_symbol_map$ENTREZID )

ids_to_symbols <- function(x, sep = "/") {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return(NA_character_)
  ids <- strsplit(x, split = sep, fixed = TRUE)[[1]]
  syms <- unname(entrez_to_symbol[ids])
  syms <- syms[!is.na(syms) & nzchar(syms)]
  if (length(syms) == 0L) return(NA_character_)
  paste(unique(syms), collapse = "/")
}

format_ora <- function(obj, analysis, direction) {
  if (is.null(obj)) {
    return(data.frame( analysis = character(0), direction = character(0), stringsAsFactors = FALSE ))
  }

  df <- as.data.frame(obj)

  if (nrow(df) == 0L) {
    df$analysis <- character(0)
    df$direction <- character(0)
    return(df)
  }

  if ("geneID" %in% names(df)) {
    df$geneID_SYMBOL <- vapply( df$geneID, ids_to_symbols, character(1) )
  }

  df$analysis <- analysis
  df$direction <- direction

  front <- c("analysis", "direction")
  df[, c(front, setdiff(names(df), front)), drop = FALSE]
}

msg("Running revised GO BP ORA: Higher foreground...")

ora_higher_go_obj <- clusterProfiler::enrichGO( gene = higher_entrez, universe = background_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
  pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE, readable = FALSE )

msg("Running revised GO BP ORA: Lower foreground...")

ora_lower_go_obj <- clusterProfiler::enrichGO( gene = lower_entrez, universe = background_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
  pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE, readable = FALSE )

ora_higher_go <- format_ora( ora_higher_go_obj, "ORA_GO_BP", "Higher" )

ora_lower_go <- format_ora( ora_lower_go_obj, "ORA_GO_BP", "Lower" )

write_csv( ora_higher_go, out_path(8, "ORA_HIGHER_GO_BP.csv") )

write_csv( ora_lower_go, out_path(9, "ORA_LOWER_GO_BP.csv") )

run_kegg_safe <- function(gene_ids, direction) {
  ans <- try( clusterProfiler::enrichKEGG( gene = gene_ids, universe = background_entrez, organism = "mmu", keyType = "ncbi-geneid", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE, use_internal_data = FALSE ), silent = TRUE )

  if (inherits(ans, "try-error")) {
    warning( "KEGG ORA failed for ", direction, ". GO BP results and scale revision remain valid. Error: ", as.character(ans) )
    return(NULL)
  }

  ans
}

msg("Running revised KEGG ORA: Higher foreground...")
ora_higher_kegg_obj <- run_kegg_safe(higher_entrez, "Higher")

msg("Running revised KEGG ORA: Lower foreground...")
ora_lower_kegg_obj <- run_kegg_safe(lower_entrez, "Lower")

ora_higher_kegg <- format_ora( ora_higher_kegg_obj, "ORA_KEGG", "Higher" )

ora_lower_kegg <- format_ora( ora_lower_kegg_obj, "ORA_KEGG", "Lower" )

write_csv( ora_higher_kegg, out_path(10, "ORA_HIGHER_KEGG.csv") )

write_csv( ora_lower_kegg, out_path(11, "ORA_LOWER_KEGG.csv") )

msg("Generating revised Figure 3B...")

fig3b <- timepoint_all

fig3b$time_index <- match( fig3b$contrast, TIMEPOINT_FILES$contrast )

set.seed(20261118)
fig3b$x_plot <- fig3b$time_index + stats::runif( nrow(fig3b), min = -0.32, max = 0.32 )

fig3b$scale_direction <- factor( as.character(fig3b$scale_direction), levels = c("Higher", "Lower", "Not significant") )

select_labels <- function(tab, direction, n = 2L) {
  x <- tab[ tab$scale_direction == direction & is.finite(tab$adj.P.Val) & !is.na(tab$gene_symbol_restored) & nzchar(tab$gene_symbol_restored), , drop = FALSE ]

  if (nrow(x) == 0L) return(x)

  x <- x[ order( x$adj.P.Val, -abs(x$transformed_scale_difference), x$gene_symbol_restored ), , drop = FALSE ]

  head(x, n)
}

label_parts <- list()
k <- 1L

for (cn in TIMEPOINT_FILES$contrast) {
  tmp <- fig3b[fig3b$contrast == cn, , drop = FALSE]

  for (dr in c("Higher", "Lower")) {
    z <- select_labels(tmp, dr, n = 2L)
    if (nrow(z) > 0L) {
      label_parts[[k]] <- z
      k <- k + 1L
    }
  }
}

fig3b_labels <- if (length(label_parts) > 0L) {
  do.call(rbind, label_parts)
} else {
  fig3b[0, , drop = FALSE]
}

write_csv_gz( fig3b, out_path(12, "Figure3B_plot_data_all_genes.csv.gz") )

write_csv( fig3b_labels, out_path(13, "Figure3B_labeled_genes.csv") )

POINT_COLORS <- c( "Higher" = "#D84A4A", "Lower" = "#3F7FBF", "Not significant" = "#C7C7C7" )

p3b <- ggplot( fig3b, aes( x = x_plot, y = transformed_scale_difference, color = scale_direction ) ) +
  geom_hline( yintercept = 0, linewidth = 0.45, color = "black" ) +
  geom_point( data = fig3b[fig3b$scale_direction == "Not significant", , drop = FALSE], size = 0.45, alpha = 0.13 ) +
  geom_point( data = fig3b[fig3b$scale_direction != "Not significant", , drop = FALSE], size = 0.65, alpha = 0.70 ) +
  ggrepel::geom_text_repel( data = fig3b_labels, aes(label = gene_symbol_restored), size = 3.1, color = "black", seed = 20261118, max.overlaps = Inf, min.segment.length = 0,
    segment.size = 0.3, box.padding = 0.25, point.padding = 0.15, show.legend = FALSE ) +
  scale_color_manual( values = POINT_COLORS, drop = FALSE ) +
  scale_x_continuous( breaks = seq_len(nrow(TIMEPOINT_FILES)), labels = TIMEPOINT_FILES$display_label, expand = expansion(mult = c(0.02, 0.04)) ) +
  labs( x = NULL, y = "Difference on supplied transformed expression scale", color = NULL ) +
  theme_classic(base_size = 13) +
  theme( axis.text.x = element_text(size = 11), axis.title.y = element_text(size = 12), legend.position = "right", legend.text = element_text(size = 11) )

ggsave( filename = out_path(20, "Figure3B_Bulk_multiVolcano_transformed_scale.pdf"), plot = p3b, width = 12.5, height = 7.2, units = "in", device = cairo_pdf )

ggsave( filename = out_path(20, "Figure3B_Bulk_multiVolcano_transformed_scale.png"), plot = p3b, width = 12.5, height = 7.2, units = "in", dpi = 600 )

msg("Generating revised Figure 3C...")

fig3c <- primary
fig3c$minus_log10_FDR <- safe_neglog10(fig3c$adj.P.Val)

CORE_LIPID_GENES <- c( "Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1" )

STEROL_GENES <- c( "Cyp51", "Sqle", "Hmgcr", "Dhcr7", "Dhcr24", "Sc5d" )

label_candidates <- data.frame( gene_symbol_restored = c(CORE_LIPID_GENES, STEROL_GENES), biological_group = c( rep("Lipid-handling core", length(CORE_LIPID_GENES)),
  rep("Sterol biosynthesis", length(STEROL_GENES)) ), stringsAsFactors = FALSE )

idx <- match( label_candidates$gene_symbol_restored, fig3c$gene_symbol_restored )

label_audit <- label_candidates
label_audit$found <- !is.na(idx)
label_audit$transformed_scale_difference <- NA_real_
label_audit$adj.P.Val <- NA_real_
label_audit$FDR_lt_0.05 <- FALSE
label_audit$scale_direction <- NA_character_

found <- which(label_audit$found)

if (length(found) > 0L) {
  label_audit$transformed_scale_difference[found] <-
    fig3c$transformed_scale_difference[idx[found]]

  label_audit$adj.P.Val[found] <-
    fig3c$adj.P.Val[idx[found]]

  label_audit$FDR_lt_0.05[found] <-
    fig3c$FDR_lt_0.05[idx[found]]

  label_audit$scale_direction[found] <-
    as.character(fig3c$scale_direction[idx[found]])
}

label_genes <- label_audit$gene_symbol_restored[ label_audit$found & label_audit$FDR_lt_0.05 ]

fig3c_labels <- fig3c[ fig3c$gene_symbol_restored %in% label_genes, , drop = FALSE ]

write_csv( fig3c, out_path(14, "Figure3C_plot_data.csv") )

write_csv( label_audit, out_path(15, "Figure3C_label_audit.csv") )

p3c <- ggplot( fig3c, aes( x = transformed_scale_difference, y = minus_log10_FDR, color = scale_direction ) ) +
  geom_point( data = fig3c[fig3c$scale_direction == "Not significant", , drop = FALSE], size = 1.1, alpha = 0.22 ) +
  geom_point( data = fig3c[fig3c$scale_direction != "Not significant", , drop = FALSE], size = 1.2, alpha = 0.80 ) +
  geom_vline( xintercept = 0, linetype = "solid", linewidth = 0.4, color = "grey45" ) +
  geom_hline( yintercept = -log10(FDR_CUTOFF), linetype = "dashed", linewidth = 0.5, color = "grey40" ) +
  ggrepel::geom_text_repel( data = fig3c_labels, aes(label = gene_symbol_restored), size = 4.0, color = "black", seed = 20261118, max.overlaps = Inf, min.segment.length = 0,
    segment.size = 0.4, box.padding = 0.35, point.padding = 0.25, show.legend = FALSE ) +
  scale_color_manual( values = POINT_COLORS, drop = FALSE ) +
  labs( x = "Difference on supplied transformed expression scale", y = expression(-log[10]("BH-adjusted P value")), color = NULL ) +
  theme_classic(base_size = 14) +
  theme( legend.position = "right", axis.title = element_text(size = 13), axis.text = element_text(size = 11) )

ggsave( filename = out_path(21, "Figure3C_Bulk_PRIMARY_volcano_transformed_scale.pdf"), plot = p3c, width = 8.8, height = 7.5, units = "in", device = cairo_pdf )

ggsave( filename = out_path(21, "Figure3C_Bulk_PRIMARY_volcano_transformed_scale.png"), plot = p3c, width = 8.8, height = 7.5, units = "in", dpi = 600 )

msg("Generating revised Figure 3E...")

MANUAL_LOWER <- c( "cholesterol biosynthetic process via lathosterol", "cholesterol biosynthetic process via desmosterol", "zymosterol biosynthetic process",
  "cholesterol biosynthetic process", "sterol biosynthetic process" )

MANUAL_HIGHER <- c( "response to lipoprotein particle", "cellular response to lipoprotein particle stimulus", "lipid localization", "lipid transport", "positive regulation of endocytosis",
  "lysosome organization", "receptor-mediated endocytosis", "phagocytosis", "microglial cell activation", "myeloid leukocyte activation" )

select_manual_terms <- function(tab, terms, direction) {
  if (!("Description" %in% names(tab))) {
    stop("ORA table is missing required Description column for Figure 3E.")
  }

  if (nrow(tab) == 0L) {
    x <- tab
    x$direction_display <- character(0)
    x$manual_order <- integer(0)
    x$FDR_lt_0.05 <- logical(0)
    x$minus_log10_pvalue <- numeric(0)
    x$minus_log10_p_adjust <- numeric(0)
    return(x)
  }

  x <- tab[ tab$Description %in% terms, , drop = FALSE ]

  if (nrow(x) == 0L) {
    return(x)
  }

  x$direction_display <- direction
  x$manual_order <- match(x$Description, terms)
  x$FDR_lt_0.05 <- ( is.finite(as.numeric(x$p.adjust)) & as.numeric(x$p.adjust) < FDR_CUTOFF )

  x <- x[x$FDR_lt_0.05, , drop = FALSE]

  if (nrow(x) == 0L) return(x)

  x$minus_log10_pvalue <- safe_neglog10(x$pvalue)
  x$minus_log10_p_adjust <- safe_neglog10(x$p.adjust)

  x
}

fig3e_lower <- select_manual_terms( ora_lower_go, MANUAL_LOWER, "Lower" )

fig3e_higher <- select_manual_terms( ora_higher_go, MANUAL_HIGHER, "Higher" )

fig3e_data <- rbind(fig3e_lower, fig3e_higher)

manual_audit_one <- function(tab, terms, direction) {
  if (!("Description" %in% names(tab))) {
    return(data.frame( direction = direction, Description = terms, present_in_complete_ORA = FALSE, FDR_lt_0.05 = FALSE, stringsAsFactors = FALSE ))
  }

  data.frame(
    direction = direction,
    Description = terms,
    present_in_complete_ORA = terms %in% tab$Description,
    FDR_lt_0.05 = vapply(
      terms,
      function(term) {
        i <- match(term, tab$Description)
        if (is.na(i)) return(FALSE)
        p <- suppressWarnings(as.numeric(tab$p.adjust[[i]]))
        is.finite(p) && p < FDR_CUTOFF
      },
      logical(1)
    ),
    stringsAsFactors = FALSE
  )
}

manual_pathway_audit <- rbind( manual_audit_one(ora_lower_go, MANUAL_LOWER, "Lower"), manual_audit_one(ora_higher_go, MANUAL_HIGHER, "Higher") )

write_csv( fig3e_data, out_path(16, "Figure3E_selected_pathways.csv") )

write_csv( manual_pathway_audit, out_path(17, "Figure3E_manual_pathway_audit.csv") )

if (nrow(fig3e_data) > 0L) {

  fig3e_data$term_plot <- paste( fig3e_data$direction_display, fig3e_data$Description, sep = "|||" )

  level_list <- character(0)

  for (dr in c("Lower", "Higher")) {
    z <- fig3e_data[ fig3e_data$direction_display == dr, , drop = FALSE ]

    if (nrow(z) > 0L) {

      z <- z[order(z$minus_log10_pvalue, decreasing = FALSE), , drop = FALSE]
      level_list <- c(level_list, z$term_plot)
    }
  }

  fig3e_data$term_plot <- factor( fig3e_data$term_plot, levels = unique(level_list) )

  fig3e_data$direction_display <- factor( fig3e_data$direction_display, levels = c("Lower", "Higher") )

  p3e <- ggplot( fig3e_data, aes( x = minus_log10_pvalue, y = term_plot, fill = minus_log10_p_adjust ) ) +
    geom_col(width = 0.78) +
    facet_grid( direction_display ~ ., scales = "free_y", space = "free_y" ) +
    scale_y_discrete( labels = function(x) sub("^.*\\|\\|\\|", "", x) ) +
    scale_fill_gradient( low = "#E7DAD3", high = "#D95F3B", name = expression(-log[10]("BH-adjusted P")) ) +
    labs( title = "Established late vs Early GO biological process enrichment", x = expression(-log[10]("P value")), y = NULL ) +
    theme_bw(base_size = 12) +
    theme( plot.title = element_text(face = "bold", size = 13), strip.text.y = element_text(face = "bold", size = 12), axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.position = "right" )

  ggsave( filename = out_path(22, "Figure3E_Bulk_GO_BP_barplot_FDR_only.pdf"), plot = p3e, width = 7.5, height = 8.0, units = "in", device = cairo_pdf )

  ggsave( filename = out_path(22, "Figure3E_Bulk_GO_BP_barplot_FDR_only.png"), plot = p3e, width = 7.5, height = 8.0, units = "in", dpi = 600 )
} else {
  warning( "None of the pre-specified Figure 3E GO terms remained FDR<0.05 ", "under the revised FDR-only foreground definition. ",
    "No Figure 3E image was generated; inspect files 08, 09 and 17." )
}

msg("Auditing Supplementary_Tables.xlsx for scale-dependent terminology...")

audit_hits <- list()
audit_k <- 1L

if (!is.na(SUPP_XLSX) && file.exists(SUPP_XLSX)) {
  sheet_names <- openxlsx::getSheetNames(SUPP_XLSX)

  for (sn in sheet_names) {
    dat <- try( openxlsx::read.xlsx( SUPP_XLSX, sheet = sn, colNames = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE, check.names = FALSE ), silent = TRUE )

    if (inherits(dat, "try-error")) {
      audit_hits[[audit_k]] <- data.frame( workbook = SUPP_XLSX, sheet = sn, cell = NA_character_, value = NA_character_, issue = "Sheet could not be read for terminology audit",
        stringsAsFactors = FALSE )
      audit_k <- audit_k + 1L
      next
    }

    if (nrow(dat) == 0L || ncol(dat) == 0L) next

    for (r in seq_len(nrow(dat))) {
      for (cc in seq_len(ncol(dat))) {
        val <- dat[[cc]][[r]]

        if (is.na(val)) next

        txt <- as.character(val)
        if (!nzchar(txt)) next

        hit <- grepl( "log\\s*2\\s*fc|logfc|fold[- ]?change|fold[- ]?changes|absolute\\s*logfc|\\|\\s*logfc\\s*\\|", txt, ignore.case = TRUE, perl = TRUE )

        if (hit) {
          audit_hits[[audit_k]] <- data.frame( workbook = SUPP_XLSX, sheet = sn, cell = paste0(openxlsx::int2col(cc), r), value = txt,
            issue = "Potential unsupported fold-change / logFC wording", stringsAsFactors = FALSE )
          audit_k <- audit_k + 1L
        }
      }
    }
  }
}

if (length(audit_hits) > 0L) {
  supplementary_term_audit <- do.call(rbind, audit_hits)
} else {
  supplementary_term_audit <- data.frame( workbook = ifelse(is.na(SUPP_XLSX), "NOT FOUND", SUPP_XLSX), sheet = NA_character_, cell = NA_character_, value = NA_character_, issue = ifelse(
    is.na(SUPP_XLSX), "Supplementary workbook not found; terminology audit not performed", "No logFC/log2FC/fold-change wording detected by text scan" ), stringsAsFactors = FALSE )
}

write_csv( supplementary_term_audit, out_path(18, "SupplementaryTables_scale_term_audit.csv") )

replacement_map <- data.frame( old_term_or_rule = c( "logFC (bulk Step46 coefficient)", "log2FC (bulk Figure 3 / bulk supplementary tables)",
  "fold change / fold-change (bulk coefficient interpretation)", "adj.P.Val < 0.05 AND |logFC| >= 0.5", "logFC >= +0.5", "logFC <= -0.5", "UP / DOWN when defined by +/-0.5 threshold" ),
  revised_term_or_rule = c( "transformed_scale_difference", "Difference on supplied transformed expression scale", "Do not interpret as fold change",
  "BH-adjusted P value < 0.05; no absolute coefficient threshold", "BH-adjusted P < 0.05 AND transformed_scale_difference > 0", "BH-adjusted P < 0.05 AND transformed_scale_difference < 0",
  "Higher / Lower on the supplied transformed expression scale" ), rationale = c( rep(
  "Upstream normalization/transformation of the deposited bulk matrix is undocumented; coefficient magnitude cannot be assumed to be log2 fold change.", 7 ) ), stringsAsFactors = FALSE )

write_csv( replacement_map, out_path(19, "SCALE_INTERPRETATION_REPLACEMENT_MAP.csv") )

replacement_xlsx <- out_path( 25, "Supplementary_replacement_data.xlsx" )

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Revision_summary")
openxlsx::writeData( wb, "Revision_summary", data.frame( item = c( "Coefficient interpretation", "Significance rule", "Absolute coefficient threshold", "limma refit", "Bulk matrix loaded",
  "GSEA rerun", "ORA rerun" ), revised_status = c( "Difference on supplied transformed expression scale", "BH-adjusted P < 0.05", "None", "No", "No",
  "No - existing GSEA uses moderated t ranking", "Yes - Higher/Lower foregrounds defined by FDR and coefficient sign only" ), stringsAsFactors = FALSE ) )

openxlsx::addWorksheet(wb, "FDR_counts")
openxlsx::writeData(wb, "FDR_counts", contrast_count_summary)

openxlsx::addWorksheet(wb, "Primary_all_genes")
openxlsx::writeData(wb, "Primary_all_genes", primary)

openxlsx::addWorksheet(wb, "Primary_FDR_significant")
openxlsx::writeData(wb, "Primary_FDR_significant", primary_sig)

openxlsx::addWorksheet(wb, "ORA_Higher_GO_BP")
openxlsx::writeData(wb, "ORA_Higher_GO_BP", ora_higher_go)

openxlsx::addWorksheet(wb, "ORA_Lower_GO_BP")
openxlsx::writeData(wb, "ORA_Lower_GO_BP", ora_lower_go)

openxlsx::addWorksheet(wb, "Figure3E_pathway_audit")
openxlsx::writeData(wb, "Figure3E_pathway_audit", manual_pathway_audit)

openxlsx::addWorksheet(wb, "Supplement_term_audit")
openxlsx::writeData(wb, "Supplement_term_audit", supplementary_term_audit)

openxlsx::addWorksheet(wb, "Replacement_map")
openxlsx::writeData(wb, "Replacement_map", replacement_map)

for (sn in openxlsx::sheets(wb)) {
  openxlsx::freezePane(wb, sn, firstRow = TRUE)
  openxlsx::setColWidths(wb, sn, cols = 1:50, widths = "auto")
}

openxlsx::saveWorkbook( wb, replacement_xlsx, overwrite = TRUE )

primary_counts <- count_summary_one(primary)

report_lines <- c( "Reviewer 1 Remaining Major Comment 4 - bulk scale revision", paste0("Script: ", PREFIX, ".R"), "", "SCOPE", "- Reused Step46 complete all-gene CSV statistics.",
  "- Did NOT read RDS/RData/Rda.", "- Did NOT read/load the bulk expression matrix.", "- Did NOT rerun limma or any upstream statistical model.", "- Did NOT depend on any ZIP archive.", "",
  "COEFFICIENT INTERPRETATION", "- Source Step46 column name: logFC", "- Revised name: transformed_scale_difference",
  "- Revised interpretation: difference on the supplied transformed expression scale", "- No fold-change or log2-fold-change interpretation is assigned.", "", "SIGNIFICANCE",
  paste0("- BH-adjusted P cutoff: ", FDR_CUTOFF), "- Absolute transformed-scale coefficient cutoff: NONE", "", "PRIMARY Established late vs Early",
  paste0("- Tested genes: ", primary_counts$tested_genes_n), paste0("- FDR<0.05, Higher: ", primary_counts$FDR_lt_0.05_higher_n),
  paste0("- FDR<0.05, Lower: ", primary_counts$FDR_lt_0.05_lower_n), paste0("- FDR<0.05, Total: ", primary_counts$FDR_lt_0.05_total_n), "", "DOWNSTREAM",
  "- Figure 3B regenerated using all genes; significance is FDR-only plus coefficient sign.", "- Figure 3C regenerated; x-axis is transformed-scale difference and +/-0.5 lines are removed.",
  "- Figure 3E GO-BP ORA rerun using FDR-only Higher/Lower foregrounds.", "- KEGG ORA rerun when KEGG access is available; failure is nonfatal and recorded.",
  "- Existing GSEA is not rerun because it was based on moderated t ranking, not |logFC| cutoff.", "", "SUPPLEMENTARY TABLE CHECK",
  paste0("- Workbook audited: ", ifelse(is.na(SUPP_XLSX), "NOT FOUND", SUPP_XLSX)), paste0("- Scale/fold-change terminology hits: ", nrow(supplementary_term_audit)),
  "- Original Supplementary_Tables.xlsx was not modified.", "- A standalone replacement-data workbook was generated in the Step147 output directory.", "", "INTERPRETIVE LIMIT",
  "- Direction (Higher/Lower) is interpretable on the supplied transformed scale.", "- The numerical coefficient magnitude must not be interpreted as log2 fold change or fold change." )

writeLines( report_lines, out_path(23, "REPORT.txt"), useBytes = TRUE )

capture.output( sessionInfo(), file = out_path(24, "sessionInfo.txt") )

all_outputs <- list.files( OUT_DIR, full.names = TRUE, recursive = FALSE, all.files = FALSE )

output_index <- data.frame( file = basename(all_outputs), path = normalizePath( all_outputs, winslash = "/", mustWork = FALSE ), size_bytes = ifelse( file.exists(all_outputs),
  file.info(all_outputs)$size, NA_real_ ), stringsAsFactors = FALSE )

output_index <- output_index[ order(output_index$file), , drop = FALSE ]

write_csv( output_index, out_path(0, "OUTPUT_INDEX.csv") )

msg("Step147 complete.")
msg("Output directory: ", OUT_DIR)
msg("No RDS was read.")
msg("No expression matrix was loaded.")
msg("No limma model was rerun.")
msg("No absolute coefficient threshold was used.")
