source("config.R", local = TRUE)

root_dir <- PROJECT_ROOT
input_file <- file.path(root_dir, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
output_dir <- file.path(root_dir, "results", "105_microglia_only_DEG")
version <- "v03"

main <- function() {
  if (!file.exists(input_file))
    stop("Missing input RDS: ", input_file)
  for (pkg in c("Seurat", "SeuratObject")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Required package missing: ", pkg)
  }
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(SeuratObject))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  report_file <- file.path(output_dir, "16_105_source_audit_REPORT_v03.txt")
  writeLines(c("105 source-stratified microglia audit", paste("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Input:", input_file), paste("Output:", output_dir), "Exploratory cell-level Wilcoxon only; NOT replicate-aware.",
    "Do NOT interpret source as biological replicate.", "Review conditions/timepoint before claiming cross-source reproducibility."),
    report_file, useBytes = TRUE)
  logmsg <- function(...) {
    line <- paste0(..., collapse = "")
    cat(line, "\n", sep = "")
    cat(line, "\n", sep = "", file = report_file, append = TRUE)
  }
  saved_files <- character(0)
  out <- function(data, name) {
    fn <- file.path(output_dir, name)
    utils::write.csv(data, file = fn, row.names = FALSE, na = "")
    saved_files <<- c(saved_files, name)
    logmsg("Saved ", name, "; rows=", nrow(data))
    invisible(fn)
  }

  logmsg("Reading Seurat object...")
  obj <- readRDS(input_file)
  if (!inherits(obj, "Seurat"))
    stop("Input is not a Seurat object")
  if (!("RNA" %in% SeuratObject::Assays(obj)))
    stop("RNA assay missing")
  m <- obj[[]]
  if (!("myeloid_state" %in% names(m)))
    stop("myeloid_state column missing")
  source_options <- c("source_dataset", "dataset", "orig.dataset", "orig.ident")
  source_fields <- intersect(source_options, names(m))
  if (!length(source_fields))
    stop("No source dataset column present")
  source_field <- source_fields[1L]
  keep <- !is.na(m$myeloid_state) & grepl(" Microglia$", as.character(m$myeloid_state)) & !is.na(m[[source_field]])
  if (!any(keep))
    stop("No annotated microglial cells with source labels")
  obj <- subset(obj, cells = rownames(m)[keep])
  rm(m)
  invisible(gc(verbose = FALSE))
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = "RNA"), error = function(e) {
    logmsg("WARNING: JoinLayers failed: ", conditionMessage(e))
    obj
  })
  DefaultAssay(obj) <- "RNA"
  layers <- tryCatch(SeuratObject::Layers(obj[["RNA"]]), error = function(e) character())
  if (!("data" %in% layers)) {
    logmsg("Normalizing microglia (RNA data layer unavailable)...")
    obj <- Seurat::NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000,
      verbose = FALSE)
  }
  if (!("data" %in% SeuratObject::Layers(obj[["RNA"]]))) {
    stop("RNA data layer not available after joining/normalization")
  }
  m <- obj[[]]
  m$myeloid_state <- as.character(m$myeloid_state)
  m$source_used <- as.character(m[[source_field]])
  logmsg("Source field: ", source_field)
  logmsg("Microglial cells: ", ncol(obj))

  md_names <- setdiff(names(m), "source_used")
  field_info <- lapply(md_names, function(field) {
    x <- as.character(m[[field]])
    ex <- unique(x[!is.na(x) & nzchar(x)])
    kind <- if (grepl("sample|replic|animal|mouse|specimen|orig.ident|gsm", field, ignore.case = TRUE)) {
      "possible_sample_or_replicate"
    } else if (grepl("dpi|time|day|stage|post.?injury", field, ignore.case = TRUE)) {
      "possible_time"
    } else if (grepl("condition|treat|group|control|plx|deplet", field, ignore.case = TRUE)) {
      "possible_condition"
    } else "other"
    data.frame(field = field, candidate_role = kind, n_unique_nonmissing = length(ex), example_values = paste(head(ex,
      4), collapse = " | "), stringsAsFactors = FALSE)
  })
  inventory <- do.call(rbind, field_info)
  out(inventory, "08_105_metadata_field_inventory_v03.csv")

  choose_field <- function(exact, pattern = NULL) {
    choices <- exact[exact %in% names(m)]
    if (length(choices))
      return(choices[1L])
    if (!is.null(pattern)) {
      choices <- grep(pattern, names(m), ignore.case = TRUE, value = TRUE)
      if (length(choices))
        return(choices[1L])
    }
    NA_character_
  }
  sample_field <- choose_field(c("sample_id", "sample_ID", "sampleID", "sample", "biological_sample",
    "replicate_id", "mouse_id", "animal_id", "orig.ident"), "sample|replic|animal|mouse|specimen|gsm")
  time_field <- choose_field(c("timepoint", "time_point", "Timepoint", "dpi", "days_post_injury", "injury_time",
    "time"), "timepoint|dpi|days_post|injury.time")
  condition_field <- choose_field(c("treatment", "treatment_group", "condition", "experimental_condition",
    "group", "experimental_group"), "treat|condition|plx|deplet|group")
  get_value <- function(f, placeholder) {
    if (is.na(f))
      return(rep(placeholder, nrow(m)))
    x <- as.character(m[[f]])
    x[is.na(x) | !nzchar(x)] <- "MISSING_VALUE"
    x
  }
  audit <- data.frame(source_dataset = m$source_used, sample_label = get_value(sample_field, "NO_SAMPLE_FIELD"),
    injury_time = get_value(time_field, "NO_TIME_FIELD"), condition = get_value(condition_field,
      "NO_CONDITION_FIELD"), myeloid_state = m$myeloid_state, n_cells = 1L, stringsAsFactors = FALSE)
  audit <- stats::aggregate(n_cells ~ source_dataset + sample_label + injury_time + condition + myeloid_state,
    data = audit, FUN = sum)
  audit <- audit[order(audit$source_dataset, audit$sample_label, audit$injury_time, audit$myeloid_state),
    , drop = FALSE]
  audit$sample_metadata_field <- if (is.na(sample_field))
    "UNAVAILABLE" else sample_field
  audit$time_metadata_field <- if (is.na(time_field))
    "UNAVAILABLE" else time_field
  audit$condition_metadata_field <- if (is.na(condition_field))
    "UNAVAILABLE" else condition_field
  out(audit, "09_105_sample_time_condition_state_counts_v03.csv")
  logmsg("Sample field (unverified replicate identity): ", audit$sample_metadata_field[1L])
  logmsg("Time field: ", audit$time_metadata_field[1L])
  logmsg("Condition field: ", audit$condition_metadata_field[1L])

  source_state <- as.data.frame(table(source_dataset = m$source_used, myeloid_state = m$myeloid_state),
    stringsAsFactors = FALSE)
  names(source_state)[3L] <- "n_cells"
  source_state <- source_state[source_state$n_cells > 0L, , drop = FALSE]
  source_state$proportion_of_source_microglia <- source_state$n_cells/as.numeric(tapply(source_state$n_cells,
    source_state$source_dataset, sum)[as.character(source_state$source_dataset)])
  out(source_state, "10_105_source_state_proportions_v03.csv")

  target <- "Clec7a+ lipid-associated Microglia"
  homeo <- "Homeostatic Microglia"
  marker_panel <- c("Clec7a", "Cst7", "Apoe", "Trem2", "Lpl", "Spp1", "Gpnmb", "Apoc1", "Fabp5", "Lgals3",
    "P2ry12", "Tmem119", "Hexb", "Cx3cr1", "Tyrobp", "Lyz2", "Csf1r", "Selenop", "Rack1", "H2az1")
  out_degs <- list()
  out_summary <- list()
  marker_rows <- list()
  sources <- sort(unique(m$source_used))
  for (source in sources) {
    logmsg("--- SOURCE: ", source, " ---")
    cells <- rownames(m)[m$source_used == source]
    small <- subset(obj, cells = cells)
    Idents(small) <- "myeloid_state"
    t <- table(as.character(SeuratObject::Idents(small)))
    n_target <- if (target %in% names(t))
      as.integer(t[[target]]) else 0L
    n_homeo <- if (homeo %in% names(t))
      as.integer(t[[homeo]]) else 0L
    n_other <- ncol(small) - n_target
    for (contrast in c("other_microglia", "homeostatic_microglia")) {
      n_ref <- if (contrast == "homeostatic_microglia")
        n_homeo else n_other
      status <- "OK"
      message_txt <- ""
      deg <- NULL
      if (n_target < 25L || n_ref < 25L) {
        status <- "SKIPPED_LOW_CELL_COUNT"
        message_txt <- "<25 cells in target or comparator"
      } else {
        deg <- tryCatch({
          Seurat::FindMarkers(small, ident.1 = target, ident.2 = if (contrast == "homeostatic_microglia")
          homeo else NULL, assay = "RNA", slot = "data", test.use = "wilcox", logfc.threshold = 0,
          min.pct = 0.1, only.pos = FALSE, return.thresh = 1, verbose = FALSE)
        }, error = function(e) {
          status <<- "FAILED"
          message_txt <<- conditionMessage(e)
          NULL
        })
      }
      if (!is.null(deg)) {
        deg$gene <- rownames(deg)
        if (!("avg_log2FC" %in% names(deg)) && "avg_logFC" %in% names(deg)) {
          deg$avg_log2FC <- deg$avg_logFC
        }
        if (!("avg_log2FC" %in% names(deg)))
          stop("Missing logFC column in FindMarkers")
        deg$significance <- ifelse(deg$p_val_adj < 0.05 & deg$avg_log2FC >= 0.25, "Up", ifelse(deg$p_val_adj <
          0.05 & deg$avg_log2FC <= -0.25, "Down", "NS"))
        deg$source_dataset <- source
        deg$comparison <- contrast
        deg <- deg[, c("source_dataset", "comparison", "gene", setdiff(names(deg), c("source_dataset",
          "comparison", "gene"))), drop = FALSE]
        rownames(deg) <- NULL
        out_degs[[length(out_degs) + 1L]] <- deg
        for (gene in marker_panel) {
          g <- deg[deg$gene == gene, , drop = FALSE]
          marker_rows[[length(marker_rows) + 1L]] <- data.frame(source_dataset = source, comparison = contrast,
          gene = gene, avg_log2FC = if (nrow(g))
            g$avg_log2FC[1L] else NA_real_, p_val_adj = if (nrow(g))
            g$p_val_adj[1L] else NA_real_, pct_target = if (nrow(g))
            g[["pct.1"]][1L] else NA_real_, pct_comparator = if (nrow(g))
            g[["pct.2"]][1L] else NA_real_, significance = if (nrow(g))
            g$significance[1L] else "NOT_TESTED", stringsAsFactors = FALSE)
        }
      }
      out_summary[[length(out_summary) + 1L]] <- data.frame(source_dataset = source, comparison = contrast,
        target_n_cells = n_target, comparator_n_cells = n_ref, n_tested = if (is.null(deg))
          0L else nrow(deg), n_up = if (is.null(deg))
          0L else sum(deg$significance == "Up", na.rm = TRUE), n_down = if (is.null(deg))
          0L else sum(deg$significance == "Down", na.rm = TRUE), status = status, details = message_txt,
        statistical_unit = "CELL; not replicate-aware", stringsAsFactors = FALSE)
      logmsg(source, " | ", contrast, " | ", status, " | target=", n_target, ", comparator=", n_ref,
        ", tested=", if (is.null(deg))
          0L else nrow(deg))
    }
    rm(small)
    invisible(gc(verbose = FALSE))
  }
  out(do.call(rbind, out_summary), "11_105_source_DEG_summary_v03.csv")
  if (length(out_degs)) {
    out(do.call(rbind, out_degs), "12_105_source_DEG_all_v03.csv")
  } else logmsg("WARNING: no source-specific DEG contrast succeeded")
  if (length(marker_rows)) {
    out(do.call(rbind, marker_rows), "13_105_core_markers_by_source_v03.csv")
  }
  methods_file <- file.path(output_dir, "14_105_methods_limitations_v03.txt")
  writeLines(c("Within-source microglia-only comparisons; target Clec7a+ annotated microglia.", "Seurat FindMarkers, RNA normalized data, Wilcoxon, logfc.threshold=0, min.pct=0.1.",
    "Significant if p_val_adj<0.05 and |avg_log2FC|>=0.25.", "Source datasets are NOT biological replicates.",
    "Sample label is inferred from metadata FIELD NAME only; must validate experimentally.", "Do not interpret significant cell-level P values as animal-level evidence.",
    "Timepoint, injury model, sex and treatment may differ across source datasets.", "Avoid interpreting cross-source variability as absent biology without design matching.",
    "No raw-expression data or cell-level identifiers are exported.", "For publication: use replicate-aware pseudobulk within comparable designs if possible."),
    methods_file, useBytes = TRUE)
  saved_files <- c(saved_files, basename(methods_file))
  session_file <- file.path(output_dir, "15_105_sessionInfo_v03.txt")
  writeLines(capture.output(utils::sessionInfo()), session_file, useBytes = TRUE)
  saved_files <- c(saved_files, basename(session_file), basename(report_file))
  logmsg("Completed. No Desktop output and no ZIP created.")
  idx <- data.frame(filename = unique(saved_files), size_bytes = as.numeric(file.info(file.path(output_dir,
    unique(saved_files)))$size), stringsAsFactors = FALSE)
  utils::write.csv(idx, file.path(output_dir, "00_105_OUTPUT_INDEX_v03.csv"), row.names = FALSE)
  invisible(output_dir)
}

tryCatch(main(), error = function(e) {
  message("105 source-stratified audit failed: ", conditionMessage(e))
  if (dir.exists(output_dir)) {
    cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = file.path(output_dir, "16_105_source_audit_REPORT_v03.txt"),
      append = TRUE)
  }
  stop(e)
})
