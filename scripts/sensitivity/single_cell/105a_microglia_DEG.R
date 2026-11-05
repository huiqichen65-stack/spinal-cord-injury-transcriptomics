source("config.R", local = TRUE)

root_dir <- PROJECT_ROOT
input_file <- file.path(root_dir, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
output_dir <- file.path(root_dir, "results", "105_microglia_only_DEG")
version <- "v02"

main <- function() {
  if (!dir.exists(root_dir))
    stop("Project root not found: ", root_dir)
  if (!file.exists(input_file))
    stop("Input RDS not found: ", input_file)

  needed <- c("Seurat", "SeuratObject")
  missing <- needed[!vapply(needed, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing))
    stop("Missing R package(s): ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(SeuratObject))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(output_dir, paste0("06_105_REPORT_", version, ".txt"))
  writeLines(c("105 microglia-only DEG", paste("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Input:", input_file), paste("Output:", output_dir), "No ZIP creation; no Desktop output; input RDS is not modified.",
    "WARNING: cell-level Wilcoxon is exploratory; biological replication is not modeled."), con = report_path,
    useBytes = TRUE)
  log_line <- function(...) {
    line <- paste0(..., collapse = "")
    cat(line, "\n", sep = "")
    cat(line, "\n", sep = "", file = report_path, append = TRUE)
  }
  write_table <- function(x, filename) {
    f <- file.path(output_dir, filename)
    utils::write.csv(x, file = f, row.names = FALSE, na = "")
    invisible(f)
  }

  obj <- readRDS(input_file)
  if (!inherits(obj, "Seurat"))
    stop("Input RDS is not a Seurat object")
  if (!("RNA" %in% SeuratObject::Assays(obj)))
    stop("RNA assay not found")
  meta <- obj[[]]
  if (!("myeloid_state" %in% colnames(meta))) {
    stop("Required metadata column not found: myeloid_state")
  }

  state_values <- as.character(meta$myeloid_state)
  microglia_keep <- !is.na(state_values) & grepl(" Microglia$", state_values)
  if (!any(microglia_keep)) {
    stop("No myeloid_state labels ending in ' Microglia' were found")
  }
  cells_keep <- rownames(meta)[microglia_keep]
  obj <- subset(obj, cells = cells_keep)
  rm(meta, state_values, cells_keep)
  invisible(gc(verbose = FALSE))
  DefaultAssay(obj) <- "RNA"

  join_message <- tryCatch({
    obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
    "RNA layers joined (or already compatible)."
  }, error = function(e) paste("JoinLayers skipped:", conditionMessage(e)))
  DefaultAssay(obj) <- "RNA"
  layers <- tryCatch(SeuratObject::Layers(obj[["RNA"]]), error = function(e) character())
  if (!("data" %in% layers)) {
    log_line("RNA 'data' layer missing; applying LogNormalize to retained microglia.")
    obj <- Seurat::NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000,
      verbose = FALSE)
  }
  layers <- tryCatch(SeuratObject::Layers(obj[["RNA"]]), error = function(e) character())
  if (length(layers) && !("data" %in% layers)) {
    stop("RNA normalized 'data' layer unavailable; inspect Seurat RNA layers")
  }

  Idents(obj) <- "myeloid_state"
  target <- "Clec7a+ lipid-associated Microglia"
  homeostatic <- "Homeostatic Microglia"
  counts <- table(as.character(SeuratObject::Idents(obj)))
  if (!(target %in% names(counts)) || !(homeostatic %in% names(counts))) {
    stop("Target/homeostatic label not found. Available labels: ", paste(names(counts), collapse = "; "))
  }
  target_n <- as.integer(counts[[target]])
  homeostatic_n <- as.integer(counts[[homeostatic]])
  other_n <- as.integer(ncol(obj) - target_n)
  if (min(target_n, other_n, homeostatic_n) < 3L) {
    stop("Insufficient cells for contrasts (target / other / homeostatic): ", paste(c(target_n, other_n,
      homeostatic_n), collapse = " / "))
  }

  state_counts <- as.data.frame(counts, stringsAsFactors = FALSE)
  names(state_counts) <- c("myeloid_state", "n_cells")
  write_table(state_counts, "01_105_microglia_state_counts_v02.csv")

  m <- obj[[]]
  source_candidates <- c("source_dataset", "dataset", "orig.dataset", "orig.ident")
  available_sources <- intersect(source_candidates, colnames(m))
  source_field <- if (length(available_sources))
    available_sources[1] else NA_character_
  if (!is.na(source_field)) {
    source_counts <- as.data.frame(table(source = as.character(m[[source_field]]), myeloid_state = as.character(m$myeloid_state)),
      stringsAsFactors = FALSE)
    names(source_counts)[3] <- "n_cells"
    source_counts$source_metadata_field <- source_field
    write_table(source_counts, "02_105_microglia_state_by_source_v02.csv")
  } else {
  }
  rm(m)

  run_contrast <- function(label, ident2 = NULL) {
    tab <- Seurat::FindMarkers(obj, ident.1 = target, ident.2 = ident2, assay = "RNA", slot = "data",
      test.use = "wilcox", logfc.threshold = 0, min.pct = 0.1, only.pos = FALSE, return.thresh = 1,
      verbose = FALSE)
    tab$gene <- rownames(tab)
    if (!("avg_log2FC" %in% names(tab)) && "avg_logFC" %in% names(tab)) {
      tab$avg_log2FC <- tab$avg_logFC
    }
    if (!("avg_log2FC" %in% names(tab))) {
      stop("No log2 fold-change column for: ", label)
    }
    tab$comparison <- label
    tab$significance <- ifelse(!is.na(tab$p_val_adj) & tab$p_val_adj < 0.05 & tab$avg_log2FC >= 0.25,
      "Up", ifelse(!is.na(tab$p_val_adj) & tab$p_val_adj < 0.05 & tab$avg_log2FC <= -0.25, "Down",
        "NS"))
    tab <- tab[, c("gene", setdiff(names(tab), "gene")), drop = FALSE]
    rownames(tab) <- NULL
    tab
  }

  deg_other <- run_contrast("Clec7a+ vs all other Microglia", ident2 = NULL)
  write_table(deg_other, "03_105_Clec7a_vs_other_microglia_DEG_all_v02.csv")
  deg_homeo <- run_contrast("Clec7a+ vs Homeostatic Microglia", ident2 = homeostatic)
  write_table(deg_homeo, "04_105_Clec7a_vs_homeostatic_microglia_DEG_all_v02.csv")

  one_summary <- function(tab, comparison, group2_n) {
    data.frame(comparison = comparison, target_n_cells = target_n, reference_n_cells = group2_n,
      tested_genes = nrow(tab), significant_up = sum(tab$significance == "Up", na.rm = TRUE), significant_down = sum(tab$significance ==
        "Down", na.rm = TRUE), nonsignificant = sum(tab$significance == "NS", na.rm = TRUE),
      statistical_unit = "cell (exploratory Wilcoxon)", stringsAsFactors = FALSE)
  }
  summary_tab <- rbind(one_summary(deg_other, "Clec7a+ vs all other Microglia", other_n), one_summary(deg_homeo,
    "Clec7a+ vs Homeostatic Microglia", homeostatic_n))
  write_table(summary_tab, "05_105_DEG_summary_v02.csv")
  log_line("Completed. Warning: these cell-level p-values do not model biological replicates.")

  session_file <- file.path(output_dir, "07_105_sessionInfo_v02.txt")
  writeLines(capture.output(utils::sessionInfo()), session_file, useBytes = TRUE)

  files <- c("01_105_microglia_state_counts_v02.csv", if (!is.na(source_field)) "02_105_microglia_state_by_source_v02.csv",
    "03_105_Clec7a_vs_other_microglia_DEG_all_v02.csv", "04_105_Clec7a_vs_homeostatic_microglia_DEG_all_v02.csv",
    "05_105_DEG_summary_v02.csv", "06_105_REPORT_v02.txt", "07_105_sessionInfo_v02.txt")
  idx <- data.frame(filename = files, size_bytes = as.numeric(file.info(file.path(output_dir, files))$size),
    stringsAsFactors = FALSE)
  utils::write.csv(idx, file.path(output_dir, "00_105_OUTPUT_INDEX_v02.csv"), row.names = FALSE)

  invisible(output_dir)
}

tryCatch(main(), error = function(e) {
  message("105 analysis failed: ", conditionMessage(e))
  if (dir.exists(output_dir)) {
    cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = file.path(output_dir, "06_105_REPORT_v02.txt"),
      append = TRUE)
  }
  stop(e)
})
