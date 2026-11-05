root_dir <- "D:/zz"
input_file <- file.path(root_dir, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
output_dir <- file.path(root_dir, "results", "105_microglia_only_DEG")

main <- function() {
  if (!file.exists(input_file))
    stop("Missing input RDS: ", input_file)
  for (pkg in c("Seurat", "SeuratObject")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Missing required R package: ", pkg)
  }
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(SeuratObject))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  report_file <- file.path(output_dir, "25_105_7dpi_REPORT_v04.txt")
  writeLines(c("105 7 dpi within-source microglia comparisons", paste0("Run started: ", format(Sys.time(),
    "%Y-%m-%d %H:%M:%S %Z")), paste0("Input: ", input_file), paste0("Output: ", output_dir), "These results are cell-level exploratory tests, not biological-replicate-aware.",
    "Do not interpret sample_id as an animal ID without original-design verification."), report_file,
    useBytes = TRUE)
  logmsg <- function(...) {
    s <- paste0(..., collapse = "")
    cat(s, "\n", sep = "")
    cat(s, "\n", sep = "", file = report_file, append = TRUE)
  }
  saved <- character(0)
  save_table <- function(x, filename) {
    utils::write.csv(x, file.path(output_dir, filename), row.names = FALSE, na = "")
    saved <<- c(saved, filename)
    logmsg("Saved: ", filename, " (", nrow(x), " rows)")
  }
  clean_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(trimws(x))] <- "MISSING"
    x
  }
  target <- "Clec7a+ lipid-associated Microglia"
  homeo <- "Homeostatic Microglia"
  marker_panel <- c("Clec7a", "Cst7", "Apoe", "Trem2", "Lpl", "Spp1", "Gpnmb", "Apoc1", "Fabp5", "Lgals3",
    "P2ry12", "Tmem119", "Hexb", "Cx3cr1", "Tyrobp", "Lyz2", "Csf1r")

  logmsg("Reading object...")
  obj <- readRDS(input_file)
  if (!inherits(obj, "Seurat"))
    stop("Input is not a Seurat object")
  if (!"RNA" %in% SeuratObject::Assays(obj))
    stop("RNA assay unavailable")
  m <- obj[[]]
  required <- c("myeloid_state", "time", "source_dataset")
  missing <- setdiff(required, colnames(m))
  if (length(missing))
    stop("Missing metadata field(s): ", paste(missing, collapse = ", "))

  src <- clean_chr(m$source_dataset)
  state <- clean_chr(m$myeloid_state)
  time <- clean_chr(m$time)

  keep <- (time == "7d") & grepl(" Microglia$", state) & src %in% c("17702045", "GSE162610", "GSE196928")
  if ("condition" %in% colnames(m)) {
    keep <- keep & clean_chr(m$condition) == "SCI"
  }
  if ("treatment" %in% colnames(m)) {

    keep <- keep & clean_chr(m$treatment) == "Control"
  }
  if (!any(keep))
    stop("No eligible 7d SCI Control microglia found; check metadata")
  obj <- subset(obj, cells = rownames(m)[keep])
  rm(m, src, state, time, keep)
  gc(verbose = FALSE)
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = "RNA"), error = function(e) {
    logmsg("WARNING: JoinLayers failed: ", conditionMessage(e))
    obj
  })
  DefaultAssay(obj) <- "RNA"
  layers <- tryCatch(SeuratObject::Layers(obj[["RNA"]]), error = function(e) character(0))
  if (!"data" %in% layers) {
    logmsg("RNA 'data' layer not found; performing log normalization on the in-memory copy")
    obj <- Seurat::NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000,
      verbose = FALSE)
  }
  if (!"data" %in% SeuratObject::Layers(obj[["RNA"]]))
    stop("RNA normalized layer unavailable")
  md <- obj[[]]
  sources <- c("17702045", "GSE162610", "GSE196928")
  missing_sources <- setdiff(sources, unique(as.character(md$source_dataset)))
  if (length(missing_sources))
    logmsg("WARNING missing source(s): ", paste(missing_sources, collapse = ", "))
  logmsg("Eligible 7 dpi microglia: ", ncol(obj))

  tab <- as.data.frame(table(source_dataset = as.character(md$source_dataset), myeloid_state = as.character(md$myeloid_state)),
    stringsAsFactors = FALSE)
  names(tab)[3] <- "n_cells"
  tab <- tab[tab$n_cells > 0, , drop = FALSE]
  tab$timepoint <- "7d"
  tab$condition <- "SCI"
  tab$treatment <- "Control"
  tab <- tab[, c("source_dataset", "timepoint", "condition", "treatment", "myeloid_state", "n_cells"),
    drop = FALSE]
  save_table(tab, "17_105_7dpi_matching_audit_v04.csv")

  sample_field <- if ("sample_id" %in% names(md))
    "sample_id" else if ("orig.ident" %in% names(md))
    "orig.ident" else NA_character_
  sample_info <- if ("sample_id_source" %in% names(md))
    clean_chr(md$sample_id_source) else rep("NOT_AVAILABLE", nrow(md))
  sample_lbl <- if (!is.na(sample_field))
    clean_chr(md[[sample_field]]) else rep("NOT_AVAILABLE", nrow(md))
  sample_data <- data.frame(source_dataset = clean_chr(md$source_dataset), sample_label = sample_lbl,
    sample_id_provenance = sample_info, myeloid_state = clean_chr(md$myeloid_state), n_cells = 1L,
    stringsAsFactors = FALSE)
  sample_data <- aggregate(n_cells ~ source_dataset + sample_label + sample_id_provenance + myeloid_state,
    data = sample_data, FUN = sum)
  sample_data <- sample_data[order(sample_data$source_dataset, sample_data$sample_label, sample_data$myeloid_state),
    , drop = FALSE]
  sample_data$sample_metadata_field <- if (is.na(sample_field))
    "NOT_AVAILABLE" else sample_field
  save_table(sample_data, "18_105_7dpi_sample_state_counts_v04.csv")

  results <- list()
  summaries <- list()
  panels <- list()
  for (source in sources) {
    logmsg("SOURCE: ", source)
    cells <- rownames(md)[as.character(md$source_dataset) == source]
    if (!length(cells)) {
      summaries[[source]] <- data.frame(source_dataset = source, target_n = 0L, comparator_n = 0L,
        tested_genes = 0L, n_up = 0L, n_down = 0L, status = "MISSING_SOURCE", stringsAsFactors = FALSE)
      next
    }
    small <- subset(obj, cells = cells)
    Idents(small) <- "myeloid_state"
    n1 <- sum(as.character(SeuratObject::Idents(small)) == target)
    n2 <- ncol(small) - n1
    if (n1 < 25L || n2 < 25L) {
      logmsg("Skipped low count: target=", n1, ", comparators=", n2)
      summaries[[source]] <- data.frame(source_dataset = source, target_n = n1, comparator_n = n2,
        tested_genes = 0L, n_up = 0L, n_down = 0L, status = "SKIPPED_LOW_CELL_COUNT", stringsAsFactors = FALSE)
      next
    }
    run_deg <- function(features = NULL, min_pct = 0.1) {
      Seurat::FindMarkers(object = small, ident.1 = target, ident.2 = NULL, assay = "RNA", slot = "data",
        test.use = "wilcox", features = features, logfc.threshold = 0, min.pct = min_pct, only.pos = FALSE,
        return.thresh = 1, verbose = FALSE)
    }
    deg <- tryCatch(run_deg(), error = function(e) {
      logmsg("ERROR source ", source, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(deg)) {
      summaries[[source]] <- data.frame(source_dataset = source, target_n = n1, comparator_n = n2,
        tested_genes = 0L, n_up = 0L, n_down = 0L, status = "FAILED", stringsAsFactors = FALSE)
      next
    }
    deg$gene <- rownames(deg)
    if (!"avg_log2FC" %in% names(deg) && "avg_logFC" %in% names(deg)) {
      deg$avg_log2FC <- deg$avg_logFC
    }
    if (!"avg_log2FC" %in% names(deg))
      stop("Seurat FindMarkers did not return logFC")
    deg$significance <- ifelse(!is.na(deg$p_val_adj) & deg$p_val_adj < 0.05 & deg$avg_log2FC >= 0.25,
      "Up", ifelse(!is.na(deg$p_val_adj) & deg$p_val_adj < 0.05 & deg$avg_log2FC <= -0.25, "Down",
        "NS"))
    deg$source_dataset <- source
    deg$timepoint <- "7d"
    results[[source]] <- deg[, c("source_dataset", "timepoint", "gene", setdiff(names(deg), c("source_dataset",
      "timepoint", "gene"))), drop = FALSE]
    summaries[[source]] <- data.frame(source_dataset = source, target_n = n1, comparator_n = n2,
      tested_genes = nrow(deg), n_up = sum(deg$significance == "Up"), n_down = sum(deg$significance ==
        "Down"), status = "OK", stringsAsFactors = FALSE)

    available <- intersect(marker_panel, rownames(small[["RNA"]]))
    panel <- if (length(available))
      tryCatch(run_deg(features = available, min_pct = 0), error = function(e) {
        logmsg("WARNING marker panel failed for ", source, ": ", conditionMessage(e))
        NULL
      }) else NULL
    if (!is.null(panel)) {
      panel$gene <- rownames(panel)
      if (!"avg_log2FC" %in% names(panel) && "avg_logFC" %in% names(panel)) {
        panel$avg_log2FC <- panel$avg_logFC
      }
    }
    one <- lapply(marker_panel, function(g) {
      present <- g %in% rownames(small[["RNA"]])
      r <- if (!is.null(panel))
        panel[panel$gene == g, , drop = FALSE] else NULL
      if (is.null(r) || nrow(r) == 0L) {
        return(data.frame(source_dataset = source, timepoint = "7d", gene = g, feature_in_RNA = present,
          avg_log2FC = NA_real_, p_val_adj = NA_real_, pct_target = NA_real_, pct_comparator = NA_real_,
          result = if (!present) "FEATURE_ABSENT" else "NOT_RETURNED", stringsAsFactors = FALSE))
      }
      significance <- if (is.na(r$p_val_adj[1]))
        "NOT_TESTED" else if (r$p_val_adj[1] < 0.05 && r$avg_log2FC[1] >= 0.25)
        "Up" else if (r$p_val_adj[1] < 0.05 && r$avg_log2FC[1] <= -0.25)
        "Down" else "NS"
      data.frame(source_dataset = source, timepoint = "7d", gene = g, feature_in_RNA = present,
        avg_log2FC = r$avg_log2FC[1], p_val_adj = r$p_val_adj[1], pct_target = r[["pct.1"]][1],
        pct_comparator = r[["pct.2"]][1], result = significance, stringsAsFactors = FALSE)
    })
    panels[[source]] <- do.call(rbind, one)
    logmsg(source, " completed; target=", n1, " comparator=", n2, " tested=", nrow(deg), " Up=",
      sum(deg$significance == "Up"), " Down=", sum(deg$significance == "Down"))
    rm(small)
    gc(verbose = FALSE)
  }
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  save_table(summary, "19_105_7dpi_source_DEG_summary_v04.csv")
  if (length(results)) {
    all_deg <- do.call(rbind, results)
    rownames(all_deg) <- NULL
    save_table(all_deg, "20_105_7dpi_source_DEG_all_v04.csv")
  } else {
    stop("No dataset-specific DEG comparison succeeded")
  }
  if (length(panels)) {
    panel_data <- do.call(rbind, panels)
    rownames(panel_data) <- NULL
    save_table(panel_data, "21_105_7dpi_core_marker_direct_tests_v04.csv")
  }

  successful <- names(results)
  if (length(successful) == 3L) {
    tested_sets <- lapply(results, function(x) as.character(x$gene))
    common_tested <- Reduce(intersect, tested_sets)
    all_up <- Reduce(intersect, lapply(results, function(x) x$gene[x$significance == "Up"]))
    all_down <- Reduce(intersect, lapply(results, function(x) x$gene[x$significance == "Down"]))
    summary_consistency <- data.frame(metric = c("sources_with_successful_DEG", "genes_tested_in_all_three",
      "up_in_all_three", "down_in_all_three"), n = c(3L, length(common_tested), length(all_up),
      length(all_down)), interpretation = c("Not three biological replicates", "Intersection of feature tests",
      "Cell-level concordance, not replicate proof", "Cell-level concordance, not replicate proof"),
      stringsAsFactors = FALSE)
  } else {
    summary_consistency <- data.frame(metric = "sources_with_successful_DEG", n = length(successful),
      interpretation = "Full three-source concordance not assessable", stringsAsFactors = FALSE)
  }
  save_table(summary_consistency, "22_105_7dpi_signature_consistency_v04.csv")
  limitations <- c("7 dpi / SCI / Control filtering was applied within each source dataset.", "No 7 dpi Clec7a+ vs homeostatic comparison: too few homeostatic cells per source.",
    "Cell-level Wilcoxon P values do not represent animal-level biological replication.", "sample_id and sample_id_source must be validated against original paper/GEO design.",
    "Matching 7 dpi does not harmonize injury model, anatomical sampling, protocol or sex.", "High between-source differences in Clec7a+ cell fractions remain unresolved.",
    "Any cell label based on Clec7a expression requires attention to circular selection.", "Reference DAM signatures remain literature-defined; do not infer a novel subtype solely",
    "from common or distinct gene markers.")
  writeLines(limitations, file.path(output_dir, "23_105_7dpi_methods_limitations_v04.txt"))
  saved <- c(saved, "23_105_7dpi_methods_limitations_v04.txt")
  capture.output(utils::sessionInfo(), file = file.path(output_dir, "24_105_7dpi_sessionInfo_v04.txt"))
  saved <- c(saved, "24_105_7dpi_sessionInfo_v04.txt")
  logmsg("Completed. 7 dpi stage-matched descriptive sensitivity; no desktop output.")
  index <- data.frame(filename = c(saved, basename(report_file)), size_bytes = as.numeric(file.info(file.path(output_dir,
    c(saved, basename(report_file))))$size), stringsAsFactors = FALSE)
  utils::write.csv(index, file.path(output_dir, "00_105_OUTPUT_INDEX_v04.csv"), row.names = FALSE,
    na = "")
  invisible(output_dir)
}

main()
