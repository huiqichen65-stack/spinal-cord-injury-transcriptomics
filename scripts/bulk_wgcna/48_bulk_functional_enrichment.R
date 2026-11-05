rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260930)

STEP48_VERSION <- "v04"
OUTPUT_DIR <- "D:/zz/results/48_Bulk_functional_enrichment"

STEP46_PRIMARY_FILE <- paste0("D:/zz/results/46_Bulk_limma_temporal_DEG_analysis/", "04_46_Bulk_limma_temporal_DEG_analysis_v06_PRIMARY_EstablishedLate_vs_Early.csv")

EXPECTED_TESTED <- 12977L
EXPECTED_UP <- 1801L
EXPECTED_DOWN <- 249L

FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5
MIN_GS_SIZE <- 10L
MAX_GS_SIZE <- 500L

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cran_pkgs <- c("BiocManager")
bioc_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "fgsea", "AnnotationDbi")

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package installation/loading failed: ", pkg)
  }
}

for (p in cran_pkgs) install_if_missing(p, bioc = FALSE)
for (p in bioc_pkgs) install_if_missing(p, bioc = TRUE)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(fgsea)
})

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

write_csv <- function(x, file) {
  utils::write.csv(x, file = file, row.names = FALSE, na = "")
}

read_table_auto <- function(path, nrows = -1L) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, nrows = nrows, fileEncoding = "UTF-8-BOM"))
  }
  if (ext %in% c("tsv", "txt")) {

    x <- try(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, nrows = nrows,
      fileEncoding = "UTF-8-BOM"), silent = TRUE)
    if (!inherits(x, "try-error") && ncol(x) > 1L)
      return(x)
    return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, nrows = nrows, fileEncoding = "UTF-8-BOM"))
  }
  stop("Unsupported table extension: ", path)
}

norm_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

resolve_col <- function(nms, candidates, required = TRUE, regex_candidates = NULL) {
  nn <- norm_name(nms)
  cc <- norm_name(candidates)

  idx <- match(cc, nn)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L)
    return(nms[idx[1]])

  if (!is.null(regex_candidates)) {
    for (pat in regex_candidates) {
      hit <- grep(pat, nn, perl = TRUE)
      if (length(hit) == 1L)
        return(nms[hit])
      if (length(hit) > 1L) {

        lens <- nchar(nn[hit])
        return(nms[hit[which.min(lens)]])
      }
    }
  }

  if (required) {
    stop("Required column not found. Candidates: ", paste(candidates, collapse = ", "), "\nAvailable columns: ",
      paste(nms, collapse = ", "))
  }
  NA_character_
}

standardize_primary <- function(df, source_label = "") {
  if (nrow(df) == 0L)
    stop("Empty input table: ", source_label)

  required_cols <- c("gene_id", "gene_symbol_restored", "logFC", "t", "adj.P.Val")

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop("Step46 v06 Primary schema mismatch.\n", "Missing required column(s): ", paste(missing_cols,
      collapse = ", "), "\nAvailable columns: ", paste(names(df), collapse = ", "), "\nSource: ",
      source_label)
  }

  out <- data.frame(gene_id_input = as.character(df$gene_id), gene_symbol = as.character(df$gene_symbol_restored),
    logFC = suppressWarnings(as.numeric(df$logFC)), t = suppressWarnings(as.numeric(df$t)), adj.P.Val = suppressWarnings(as.numeric(df$adj.P.Val)),
    stringsAsFactors = FALSE)

  if ("P.Value" %in% names(df)) {
    out$P.Value <- suppressWarnings(as.numeric(df$P.Value))
  } else {
    out$P.Value <- NA_real_
  }

  bad_gene_id <- is.na(out$gene_id_input) | !nzchar(out$gene_id_input)
  if (any(bad_gene_id)) {
    stop("Missing/blank Step46 gene_id found in ", sum(bad_gene_id), " rows. Source: ", source_label)
  }

  if (anyDuplicated(out$gene_id_input)) {
    dups <- unique(out$gene_id_input[duplicated(out$gene_id_input)])
    stop("Step46 gene_id is unexpectedly duplicated. First duplicated IDs: ", paste(head(dups, 10),
      collapse = ", "), "\nSource: ", source_label)
  }

  bad_symbol <- is.na(out$gene_symbol) | !nzchar(out$gene_symbol)
  if (any(bad_symbol)) {
    stop("Missing/blank Step46 gene_symbol_restored found in ", sum(bad_symbol), " rows. Source: ",
      source_label)
  }

  bad_num <- !is.finite(out$logFC) | !is.finite(out$t) | !is.finite(out$adj.P.Val)
  if (any(bad_num)) {
    stop("Non-finite logFC/t/adj.P.Val found in ", sum(bad_num), " rows. Source: ", source_label)
  }

  out
}

find_contrast_col <- function(df) {
  resolve_col(names(df), c("contrast", "comparison", "coef", "coefficient", "test", "Comparison"),
    required = FALSE)
}

filter_primary_contrast_if_needed <- function(df, path) {
  ccol <- find_contrast_col(df)
  if (is.na(ccol))
    return(df)

  vals <- unique(as.character(df[[ccol]]))
  keep_label <- grepl("established", vals, ignore.case = TRUE) & grepl("late", vals, ignore.case = TRUE) &
    grepl("early", vals, ignore.case = TRUE) & !grepl("sens|exclude|review|flag", vals, ignore.case = TRUE)

  if (sum(keep_label) != 1L) {

    keep_label2 <- grepl("primary", vals, ignore.case = TRUE) & !grepl("sens|exclude|review|flag",
      vals, ignore.case = TRUE)
    if (sum(keep_label2) == 1L) {
      chosen <- vals[keep_label2]
    } else {
      stop("A contrast column was found but a unique Primary contrast could not be identified.\n",
        "File: ", path, "\n", "Contrast labels:\n  - ", paste(vals, collapse = "\n  - "))
    }
  } else {
    chosen <- vals[keep_label]
  }

  df[df[[ccol]] == chosen, , drop = FALSE]
}

validate_frozen_primary <- function(tab, source_path, stop_on_fail = TRUE) {
    n_tested <- nrow(tab)
    n_up <- sum(tab$adj.P.Val < FDR_CUTOFF & tab$logFC >= LOGFC_CUTOFF)
    n_down <- sum(tab$adj.P.Val < FDR_CUTOFF & tab$logFC <= -LOGFC_CUTOFF)
    frozen_markers <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")
    idx <- match(tolower(frozen_markers), tolower(tab$gene_symbol))
    marker_ok <- all(!is.na(idx)) && all(tab$logFC[idx] > 0) && all(tab$adj.P.Val[idx] < FDR_CUTOFF)
    ok <- (n_tested == EXPECTED_TESTED && n_up == EXPECTED_UP && n_down == EXPECTED_DOWN && marker_ok)
    info <- data.frame(source_file = source_path, tested_n = n_tested, up_n = n_up, down_n = n_down, 
        six_marker_check = marker_ok, frozen_match = ok, stringsAsFactors = FALSE)
    if (!ok && stop_on_fail) {
        stop("Step46 Primary validation failed.\n", "Source: ", source_path, "\n", "Observed tested/up/down = ", 
            n_tested, "/", n_up, "/", n_down, "\n", "Expected tested/up/down = ", EXPECTED_TESTED, "/", 
            EXPECTED_UP, "/", EXPECTED_DOWN, "\n", "Six marker check = ", marker_ok, "\n", "The script will not proceed with an unverified contrast.")
    }
    info
}

load_step46_primary_direct <- function() {
  if (!file.exists(STEP46_PRIMARY_FILE)) {
    stop("Step46 Primary file not found:\n", STEP46_PRIMARY_FILE, "\nPlease check that D:/zz/results/46_Bulk_limma_temporal_DEG_analysis/ exists.")
  }

  raw <- read_table_auto(STEP46_PRIMARY_FILE)

  schema_file <- file.path(OUTPUT_DIR, "13_48_Bulk_functional_enrichment_v04_Step46_input_schema.txt")
  schema_lines <- c(paste0("Input file: ", STEP46_PRIMARY_FILE), paste0("Rows: ", nrow(raw)), paste0("Columns: ",
    ncol(raw)), "", "Column names:", paste0(seq_along(names(raw)), ". ", names(raw)))
  writeLines(schema_lines, schema_file, useBytes = TRUE)

  raw <- filter_primary_contrast_if_needed(raw, STEP46_PRIMARY_FILE)
  tab <- standardize_primary(raw, STEP46_PRIMARY_FILE)

  audit_row <- validate_frozen_primary(tab, STEP46_PRIMARY_FILE, stop_on_fail = FALSE)
  audit_row$input_rows_raw <- nrow(raw)
  audit_row$input_columns_raw <- ncol(raw)

  audit_file <- file.path(OUTPUT_DIR, "14_48_Bulk_functional_enrichment_v04_Step46_input_validation_audit.csv")
  write_csv(audit_row, audit_file)

  if (!isTRUE(audit_row$frozen_match[1])) {
    frozen_markers <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")
    marker_rows <- tab[tolower(tab$gene_symbol) %in% tolower(frozen_markers), , drop = FALSE]
    write_csv(marker_rows, file.path(OUTPUT_DIR, "15_48_Bulk_functional_enrichment_v04_Step46_core_marker_validation.csv"))

    stop("The exact Step46 Primary file and schema were read successfully, but validation failed.\n",
      "Observed tested/up/down = ", audit_row$tested_n, "/", audit_row$up_n, "/", audit_row$down_n,
      "\n", "Expected tested/up/down = 12977/1801/249\n", "Six-marker check = ", audit_row$six_marker_check,
      "\n", "See files 13, 14 and 15 in the Step48 output directory for diagnostics.")
  }

  list(path = normalizePath(STEP46_PRIMARY_FILE, winslash = "/", mustWork = TRUE), table = tab, audit = audit_row)
}

primary_in <- load_step46_primary_direct()
primary_path <- primary_in$path
primary <- primary_in$table
step46_candidate_audit <- primary_in$audit

validate_frozen_primary(primary, primary_path, stop_on_fail = TRUE)

input_keytype <- "SYMBOL"
primary$gene_id_clean <- primary$gene_symbol

map_raw <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(primary$gene_id_clean), keytype = input_keytype,
  columns = "ENTREZID")

if (!(input_keytype %in% names(map_raw))) {
  stop("Annotation mapping did not return the expected key column: ", input_keytype, "\nReturned columns: ",
    paste(names(map_raw), collapse = ", "))
}

names(map_raw)[names(map_raw) == input_keytype] <- "gene_id_clean"
map_raw <- unique(map_raw[, c("gene_id_clean", "ENTREZID"), drop = FALSE])

map_raw <- map_raw[!is.na(map_raw$gene_id_clean) & nzchar(map_raw$gene_id_clean), , drop = FALSE]

mapped_input_ids <- unique(map_raw$gene_id_clean[!is.na(map_raw$ENTREZID) & nzchar(map_raw$ENTREZID)])

background_n <- nrow(primary)
mapped_background_n <- sum(primary$gene_id_clean %in% mapped_input_ids)

up_idx <- primary$adj.P.Val < FDR_CUTOFF & primary$logFC >= LOGFC_CUTOFF
down_idx <- primary$adj.P.Val < FDR_CUTOFF & primary$logFC <= -LOGFC_CUTOFF

up_ids <- primary$gene_id_clean[up_idx]
down_ids <- primary$gene_id_clean[down_idx]

mapped_up_n <- sum(up_ids %in% mapped_input_ids)
mapped_down_n <- sum(down_ids %in% mapped_input_ids)

background_entrez <- unique(map_raw$ENTREZID[map_raw$gene_id_clean %in% primary$gene_id_clean & !is.na(map_raw$ENTREZID) &
  nzchar(map_raw$ENTREZID)])
up_entrez <- unique(map_raw$ENTREZID[map_raw$gene_id_clean %in% up_ids & !is.na(map_raw$ENTREZID) & nzchar(map_raw$ENTREZID)])
down_entrez <- unique(map_raw$ENTREZID[map_raw$gene_id_clean %in% down_ids & !is.na(map_raw$ENTREZID) &
  nzchar(map_raw$ENTREZID)])

mapping_qc <- data.frame(metric = c("input_file", "annotation_keytype", "background_n_tested_genes",
  "background_n_mapped_genes", "background_unique_ENTREZ_n", "UP_foreground_n", "UP_foreground_mapped_genes_n",
  "UP_foreground_unique_ENTREZ_n", "DOWN_foreground_n", "DOWN_foreground_mapped_genes_n", "DOWN_foreground_unique_ENTREZ_n",
  "GSEA_rank_source", "ORA_FDR_cutoff", "ORA_abs_logFC_cutoff"), value = c(primary_path, input_keytype,
  background_n, mapped_background_n, length(background_entrez), sum(up_idx), mapped_up_n, length(up_entrez),
  sum(down_idx), mapped_down_n, length(down_entrez), "limma moderated t statistic", FDR_CUTOFF, LOGFC_CUTOFF),
  stringsAsFactors = FALSE)

mapping_detail <- merge(unique(primary[, c("gene_id_input", "gene_symbol", "gene_id_clean")]), map_raw,
  by = "gene_id_clean", all.x = TRUE, sort = FALSE)
mapping_detail <- mapping_detail[match(primary$gene_id_clean, mapping_detail$gene_id_clean), , drop = FALSE]

ora_up_go_obj <- clusterProfiler::enrichGO(gene = up_entrez, universe = background_entrez, OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GS_SIZE,
  maxGSSize = MAX_GS_SIZE, readable = FALSE)

ora_down_go_obj <- clusterProfiler::enrichGO(gene = down_entrez, universe = background_entrez, OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GS_SIZE,
  maxGSSize = MAX_GS_SIZE, readable = FALSE)

run_kegg_retry <- function(fun, label, attempts = 3L) {
  last_err <- NULL
  for (i in seq_len(attempts)) {
    ans <- try(fun(), silent = TRUE)
    if (!inherits(ans, "try-error"))
      return(ans)
    last_err <- ans
    msg(label, " failed on attempt ", i, "/", attempts, "; retrying...")
    Sys.sleep(2)
  }
  stop(label, " failed after ", attempts, " attempts.\n", "KEGG analysis requires access to KEGG resources.\n",
    as.character(last_err))
}

ora_up_kegg_obj <- run_kegg_retry(function() clusterProfiler::enrichKEGG(gene = up_entrez, universe = background_entrez,
  organism = "mmu", keyType = "ncbi-geneid", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE, use_internal_data = FALSE), "ORA UP KEGG")

ora_down_kegg_obj <- run_kegg_retry(function() clusterProfiler::enrichKEGG(gene = down_entrez, universe = background_entrez,
  organism = "mmu", keyType = "ncbi-geneid", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE, use_internal_data = FALSE), "ORA DOWN KEGG")

rank_long <- merge(primary[, c("gene_id_input", "gene_symbol", "gene_id_clean", "t", "logFC", "adj.P.Val")],
  map_raw, by = "gene_id_clean", all.x = FALSE, all.y = FALSE, sort = FALSE)

rank_long <- rank_long[!is.na(rank_long$ENTREZID) & nzchar(rank_long$ENTREZID) & is.finite(rank_long$t),
  , drop = FALSE]

rank_long <- rank_long[order(rank_long$ENTREZID, -abs(rank_long$t), rank_long$gene_id_input), , drop = FALSE]
rank_unique <- rank_long[!duplicated(rank_long$ENTREZID), , drop = FALSE]

geneList <- rank_unique$t
names(geneList) <- rank_unique$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

if (length(geneList) < 1000L) {
  stop("Too few mapped genes for GSEA: ", length(geneList))
}

tie_fraction <- 1 - length(unique(geneList))/length(geneList)

gsea_go_obj <- clusterProfiler::gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
  ont = "BP", exponent = 1, minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE, eps = 0, pvalueCutoff = 1,
  pAdjustMethod = "BH", verbose = FALSE, seed = TRUE, by = "fgsea")

gsea_kegg_obj <- run_kegg_retry(function() clusterProfiler::gseKEGG(geneList = geneList, organism = "mmu",
  keyType = "ncbi-geneid", exponent = 1, minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE, eps = 0,
  pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE, seed = TRUE, by = "fgsea", use_internal_data = FALSE),
  "GSEA KEGG")

entrez_symbol_map <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(c(background_entrez, names(geneList))),
  keytype = "ENTREZID", columns = c("ENTREZID", "SYMBOL"))
entrez_symbol_map <- entrez_symbol_map[!is.na(entrez_symbol_map$ENTREZID) & !is.na(entrez_symbol_map$SYMBOL),
  , drop = FALSE]
entrez_symbol_map <- entrez_symbol_map[!duplicated(entrez_symbol_map$ENTREZID), ]
entrez_to_symbol <- setNames(entrez_symbol_map$SYMBOL, entrez_symbol_map$ENTREZID)

ids_to_symbols <- function(x, sep = "/") {
  if (is.na(x) || !nzchar(x))
    return(NA_character_)
  ids <- strsplit(x, split = sep, fixed = TRUE)[[1]]
  syms <- unname(entrez_to_symbol[ids])
  syms <- syms[!is.na(syms) & nzchar(syms)]
  if (length(syms) == 0L)
    return(NA_character_)
  paste(unique(syms), collapse = "/")
}

format_ora <- function(obj, analysis, direction) {
  df <- as.data.frame(obj)
  if (nrow(df) == 0L) {
    df$analysis <- character(0)
    df$direction <- character(0)
    return(df)
  }
  if ("geneID" %in% names(df)) {
    df$geneID_SYMBOL <- vapply(df$geneID, ids_to_symbols, character(1))
  }
  df$analysis <- analysis
  df$direction <- direction
  front <- c("analysis", "direction")
  df[, c(front, setdiff(names(df), front)), drop = FALSE]
}

format_gsea <- function(obj, analysis) {
  df <- as.data.frame(obj)
  if (nrow(df) == 0L) {
    df$analysis <- character(0)
    df$direction <- character(0)
    return(df)
  }
  df$direction <- ifelse(df$NES > 0, "positive", "negative")
  if ("core_enrichment" %in% names(df)) {
    df$core_enrichment_SYMBOL <- vapply(df$core_enrichment, ids_to_symbols, character(1))
  }
  df$analysis <- analysis
  front <- c("analysis", "direction")
  df[, c(front, setdiff(names(df), front)), drop = FALSE]
}

ora_up_go <- format_ora(ora_up_go_obj, "ORA_GO_BP", "UP")
ora_up_kegg <- format_ora(ora_up_kegg_obj, "ORA_KEGG", "UP")
ora_down_go <- format_ora(ora_down_go_obj, "ORA_GO_BP", "DOWN")
ora_down_kegg <- format_ora(ora_down_kegg_obj, "ORA_KEGG", "DOWN")

gsea_go <- format_gsea(gsea_go_obj, "GSEA_GO_BP")
gsea_kegg <- format_gsea(gsea_kegg_obj, "GSEA_KEGG")

keyword_groups <- list(lipid_handling_transport = paste(c("lipid transport", "cholesterol transport",
  "cholesterol efflux", "lipoprotein", "lipid catabolic", "neutral lipid", "fatty acid uptake", "lipid localization",
  "lipid homeostasis"), collapse = "|"), lysosomal_processing = paste(c("lysosom", "phagosome", "endocyt",
  "endosome", "autophag"), collapse = "|"), myeloid_activation = paste(c("microgl", "macrophage", "myeloid",
  "phagocyt"), collapse = "|"), sterol_biosynthesis = paste(c("cholesterol biosynth", "sterol biosynth",
  "mevalonate", "isoprenoid biosynth", "isoprenoid metabolic", "terpenoid backbone", "lipid biosynthetic process"),
  collapse = "|"))

collect_key_matches <- function(df, source_name) {
  if (nrow(df) == 0L || !"Description" %in% names(df))
    return(NULL)

  out <- list()
  for (grp in names(keyword_groups)) {
    hit <- grepl(keyword_groups[[grp]], df$Description, ignore.case = TRUE, perl = TRUE)
    if (!any(hit))
      next

    tmp <- df[hit, , drop = FALSE]
    tmp$key_group <- grp
    tmp$source_table <- source_name

    if ("p.adjust" %in% names(tmp)) {
      tmp$BH_adjusted_P <- tmp$p.adjust
    } else {
      tmp$BH_adjusted_P <- NA_real_
    }

    tmp$significant_BH_0.05 <- !is.na(tmp$BH_adjusted_P) & tmp$BH_adjusted_P < 0.05
    out[[length(out) + 1L]] <- tmp
  }

  if (length(out) == 0L)
    return(NULL)

  z <- do.call(rbind, lapply(out, function(x) {
    data.frame(source_table = x$source_table, key_group = x$key_group, analysis = if ("analysis" %in%
      names(x))
      x$analysis else NA_character_, direction = if ("direction" %in% names(x))
      x$direction else NA_character_, ID = if ("ID" %in% names(x))
      x$ID else NA_character_, Description = x$Description, NES = if ("NES" %in% names(x))
      x$NES else NA_real_, pvalue = if ("pvalue" %in% names(x))
      x$pvalue else NA_real_, p.adjust = if ("p.adjust" %in% names(x))
      x$p.adjust else NA_real_, setSize = if ("setSize" %in% names(x))
      x$setSize else NA_real_, Count = if ("Count" %in% names(x))
      x$Count else NA_real_, core_enrichment = if ("core_enrichment" %in% names(x))
      x$core_enrichment else NA_character_, core_enrichment_SYMBOL = if ("core_enrichment_SYMBOL" %in% names(x))
      x$core_enrichment_SYMBOL else NA_character_, geneID = if ("geneID" %in% names(x))
      x$geneID else NA_character_, geneID_SYMBOL = if ("geneID_SYMBOL" %in% names(x))
      x$geneID_SYMBOL else NA_character_, significant_BH_0.05 = x$significant_BH_0.05, stringsAsFactors = FALSE)
  }))

  z
}

key_list <- Filter(Negate(is.null), list(collect_key_matches(ora_up_go, "ORA_UP_GO_BP"), collect_key_matches(ora_up_kegg,
  "ORA_UP_KEGG"), collect_key_matches(ora_down_go, "ORA_DOWN_GO_BP"), collect_key_matches(ora_down_kegg,
  "ORA_DOWN_KEGG"), collect_key_matches(gsea_go, "GSEA_GO_BP"), collect_key_matches(gsea_kegg, "GSEA_KEGG")))

if (length(key_list) > 0L) {
  key_summary <- do.call(rbind, key_list)
  key_summary <- key_summary[order(key_summary$key_group, !key_summary$significant_BH_0.05, key_summary$p.adjust,
    -abs(ifelse(is.na(key_summary$NES), 0, key_summary$NES))), , drop = FALSE]
} else {
  key_summary <- data.frame(source_table = character(0), key_group = character(0), analysis = character(0),
    direction = character(0), ID = character(0), Description = character(0), NES = numeric(0), pvalue = numeric(0),
    p.adjust = numeric(0), setSize = numeric(0), Count = numeric(0), core_enrichment = character(0),
    core_enrichment_SYMBOL = character(0), geneID = character(0), geneID_SYMBOL = character(0), significant_BH_0.05 = logical(0),
    stringsAsFactors = FALSE)
}

files <- list(mapping_qc = "01_48_Bulk_functional_enrichment_v04_input_and_mapping_QC.csv", ora_up_go = "02_48_Bulk_functional_enrichment_v04_ORA_UP_GO_BP.csv",
  ora_up_kegg = "03_48_Bulk_functional_enrichment_v04_ORA_UP_KEGG.csv", ora_down_go = "04_48_Bulk_functional_enrichment_v04_ORA_DOWN_GO_BP.csv",
  ora_down_kegg = "05_48_Bulk_functional_enrichment_v04_ORA_DOWN_KEGG.csv", gsea_go = "06_48_Bulk_functional_enrichment_v04_GSEA_GO_BP.csv",
  gsea_kegg = "07_48_Bulk_functional_enrichment_v04_GSEA_KEGG.csv", key_summary = "08_48_Bulk_functional_enrichment_v04_key_lipid_lysosome_sterol_pathway_summary.csv",
  mapping_detail = "09_48_Bulk_functional_enrichment_v04_symbol_ENTREZ_mapping.csv", analysis_object = "10_48_Bulk_functional_enrichment_v04_analysis_object.rds",
  report = "11_48_Bulk_functional_enrichment_v04_REPORT.txt", session = "12_48_Bulk_functional_enrichment_v04_sessionInfo.txt",
  step46_schema = "13_48_Bulk_functional_enrichment_v04_Step46_input_schema.txt", step46_audit = "14_48_Bulk_functional_enrichment_v04_Step46_input_validation_audit.csv",
  step46_markers = "15_48_Bulk_functional_enrichment_v04_Step46_core_marker_validation.csv")

write_csv(mapping_qc, file.path(OUTPUT_DIR, files$mapping_qc))
write_csv(ora_up_go, file.path(OUTPUT_DIR, files$ora_up_go))
write_csv(ora_up_kegg, file.path(OUTPUT_DIR, files$ora_up_kegg))
write_csv(ora_down_go, file.path(OUTPUT_DIR, files$ora_down_go))
write_csv(ora_down_kegg, file.path(OUTPUT_DIR, files$ora_down_kegg))
write_csv(gsea_go, file.path(OUTPUT_DIR, files$gsea_go))
write_csv(gsea_kegg, file.path(OUTPUT_DIR, files$gsea_kegg))
write_csv(key_summary, file.path(OUTPUT_DIR, files$key_summary))
write_csv(mapping_detail, file.path(OUTPUT_DIR, files$mapping_detail))

if (!file.exists(file.path(OUTPUT_DIR, files$step46_schema))) {
  writeLines(c(paste0("Input file: ", primary_path), paste0("Rows after standardization: ", nrow(primary)),
    paste0("Standardized columns: ", paste(names(primary), collapse = ", ")), "Annotation gene field: gene_symbol_restored -> standardized gene_symbol",
    "Stable unique row field: gene_id"), file.path(OUTPUT_DIR, files$step46_schema), useBytes = TRUE)
}

frozen_markers_for_output <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")
marker_validation_output <- primary[tolower(primary$gene_symbol) %in% tolower(frozen_markers_for_output),
  , drop = FALSE]
write_csv(marker_validation_output, file.path(OUTPUT_DIR, files$step46_markers))

if (is.null(step46_candidate_audit)) {
  step46_candidate_audit <- validate_frozen_primary(primary, primary_path, stop_on_fail = FALSE)
}
write_csv(step46_candidate_audit, file.path(OUTPUT_DIR, files$step46_audit))

analysis_object <- list(step = "Step48", version = STEP48_VERSION, frozen_contrast = "Established late (14+28 dpi) vs Early (Uninjured+15min+1dpi+3dpi)",
  step46_primary_source = primary_path, step46_gene_identity = list(stable_unique_id = "gene_id", annotation_symbol = "gene_symbol_restored",
    annotation_keytype = "SYMBOL"), step46_input_validation_audit = step46_candidate_audit, thresholds = list(ORA_FDR = FDR_CUTOFF,
    ORA_abs_logFC = LOGFC_CUTOFF, minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE), input_primary = primary,
  mapping_qc = mapping_qc, mapping_detail = mapping_detail, gene_sets_input = list(background_entrez = background_entrez,
    up_entrez = up_entrez, down_entrez = down_entrez, GSEA_geneList_moderated_t = geneList, GSEA_tie_fraction = tie_fraction,
    GSEA_duplicate_entrez_rule = "retain row with largest absolute moderated t"), ORA = list(UP_GO_BP = ora_up_go_obj,
    UP_KEGG = ora_up_kegg_obj, DOWN_GO_BP = ora_down_go_obj, DOWN_KEGG = ora_down_kegg_obj), GSEA = list(GO_BP = gsea_go_obj,
    KEGG = gsea_kegg_obj), result_tables = list(ORA_UP_GO_BP = ora_up_go, ORA_UP_KEGG = ora_up_kegg,
    ORA_DOWN_GO_BP = ora_down_go, ORA_DOWN_KEGG = ora_down_kegg, GSEA_GO_BP = gsea_go, GSEA_KEGG = gsea_kegg,
    key_pathway_summary = key_summary), keyword_groups = keyword_groups)

saveRDS(analysis_object, file = file.path(OUTPUT_DIR, files$analysis_object), compress = "xz")

sig_n <- function(df) {
  if (!"p.adjust" %in% names(df))
    return(NA_integer_)
  sum(!is.na(df$p.adjust) & df$p.adjust < 0.05)
}

gsea_sig_dir <- function(df, positive = TRUE) {
  if (nrow(df) == 0L)
    return(0L)
  dir_ok <- if (positive)
    df$NES > 0 else df$NES < 0
  sum(!is.na(df$p.adjust) & df$p.adjust < 0.05 & dir_ok)
}

count_key_sig <- function(group, direction = NULL, analysis_regex = NULL) {
  if (nrow(key_summary) == 0L)
    return(0L)
  keep <- key_summary$key_group == group & key_summary$significant_BH_0.05
  if (!is.null(direction))
    keep <- keep & key_summary$direction == direction
  if (!is.null(analysis_regex))
    keep <- keep & grepl(analysis_regex, key_summary$analysis)
  sum(keep, na.rm = TRUE)
}

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = file.path(OUTPUT_DIR, files$report), useBytes = TRUE)

sink(file.path(OUTPUT_DIR, files$session))
print(sessionInfo())
sink()

output_index <- data.frame(file = c("00_48_Bulk_functional_enrichment_v04_OUTPUT_INDEX.csv", files$mapping_qc,
  files$ora_up_go, files$ora_up_kegg, files$ora_down_go, files$ora_down_kegg, files$gsea_go, files$gsea_kegg,
  files$key_summary, files$mapping_detail, files$analysis_object, files$report, files$session, files$step46_schema,
  files$step46_audit, files$step46_markers), description = c("Index of Step48 v04 outputs", "Input/foreground/background mapping QC and analysis settings",
  "Complete ORA results: UP genes, GO Biological Process", "Complete ORA results: UP genes, KEGG",
  "Complete ORA results: DOWN genes, GO Biological Process", "Complete ORA results: DOWN genes, KEGG",
  "Complete preranked GSEA results: GO Biological Process; moderated-t ranking", "Complete preranked GSEA results: KEGG; moderated-t ranking",
  "Pre-specified lipid/lysosome/sterol/myeloid keyword summary; complete tables remain authoritative",
  "Detailed tested-gene to Mus musculus SYMBOL/ENTREZ mapping", "Complete Step48 analysis object including enrichment objects",
  "Plain-text Step48 QC/result report", "R sessionInfo", "Actual schema/column names of the exact Step46 Primary input",
  "Strict -input validation audit", "six-core-marker validation rows"), stringsAsFactors = FALSE)

write_csv(output_index, file.path(OUTPUT_DIR, "00_48_Bulk_functional_enrichment_v04_OUTPUT_INDEX.csv"))
