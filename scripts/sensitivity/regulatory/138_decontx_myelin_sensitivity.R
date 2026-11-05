options(stringsAsFactors = FALSE, width = 180)
root <- "D:/zz"
version <- "v01"
label <- "138_DecontX_multilayer_cell_ID_audit"
out <- file.path(root, "results", label)
tabs <- file.path(out, "tables")
logs <- file.path(out, "logs")
for (d in c(out, tabs, logs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(logs, "00_138_run_log_v01.log")
writeLines(paste("STEP 138 START", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), log_path)
write_csv <- function(x, filename) {
  if (!is.data.frame(x))
    stop("Internal error: output is not a data.frame: ", filename)
  utils::write.csv(x, file.path(tabs, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
loc <- function(...) file.path(root, ...)
input <- c(atlas = loc("results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds"),
  myelin_leading = loc("results", "102_Zfp148_GO_GSEA_fullrank_v01", "07_102_Zfp148_primary_GO_leading_edge_FDR05_v01.csv"),
  raw_17702045 = loc("results", "02_QC", "02_QC_17702045_v01.rds"), corr_17702045 = loc("results",
    "03_DecontX", "03_DecontX_17702045_v01.rds"), raw_GSE162610 = loc("results", "02_QC", "02_QC_GSE162610_v01.rds"),
  corr_GSE162610 = loc("results", "03_DecontX", "03_DecontX_GSE162610_v01.rds"), raw_GSE196928_Control = loc("results",
    "02_QC", "02_QC_GSE196928_Control_v01.rds"), corr_GSE196928_Control = loc("results", "03_DecontX",
    "03_DecontX_GSE196928_Control_v01.rds"), raw_GSE196928_PLX5622 = loc("results", "02_QC", "02_QC_GSE196928_PLX5622_v01.rds"),
  corr_GSE196928_PLX5622 = loc("results", "03_DecontX", "03_DecontX_GSE196928_PLX5622_v01.rds"))
write_csv(data.frame(input_key = names(input), path = unname(input), exists = file.exists(input), bytes = as.numeric(file.info(input)$size)),
  "00_138_input_provenance_v01.csv")
if (!all(file.exists(input[c("atlas", "myelin_leading")]))) {
  stop("Missing required atlas or GO input: consult 00_138_input_provenance_v01.csv")
}
for (p in c("Matrix", "SeuratObject")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Missing installed package ", p, "; script installs nothing")
}

first_col <- function(m, possibilities) {
  found <- possibilities[possibilities %in% names(m)]
  if (length(found))
    found[1L] else NA_character_
}
read_col <- function(m, nm) {
  if (is.na(nm) || !(nm %in% names(m)))
    return(rep(NA_character_, nrow(m)))
  z <- trimws(as.character(m[[nm]]))
  z[z %in% c("", "NA", "NaN", "NULL", "Unknown", "unknown")] <- NA_character_
  z
}
numeric_col <- function(m, nm) {
  if (!(nm %in% names(m)))
    return(rep(NA_real_, nrow(m)))
  suppressWarnings(as.numeric(as.character(m[[nm]])))
}
meta_of <- function(object, cells) {
  if (inherits(object, "Seurat")) {
    m <- object@meta.data
  } else if (inherits(object, "SummarizedExperiment")) {
    m <- as.data.frame(SummarizedExperiment::colData(object))
  } else {
    m <- data.frame(row.names = cells)
  }
  if (is.null(rownames(m)) || anyDuplicated(rownames(m)) || !all(cells %in% rownames(m)))
    stop("Missing/duplicated or unalignable cell metadata IDs")
  m[cells, , drop = FALSE]
}
assay_counts <- function(object, corrected = FALSE) {

  if (inherits(object, "Seurat")) {
    anames <- names(object@assays)
    selected_assay <- if (corrected && "decontXcounts" %in% anames)
      "decontXcounts" else "RNA"
    if (!(selected_assay %in% anames))
      stop("Required assay missing: ", selected_assay, " ; assays=", paste(anames, collapse = ";"))
    a <- object[[selected_assay]]
    layers <- tryCatch(SeuratObject::Layers(a, search = NA), error = function(e) character())
    chosen <- layers[grepl("^counts($|\\.)", layers)]
    if (!length(chosen))
      stop("No count layers in ", selected_assay, "; layers=", paste(layers, collapse = ";"))
    parts <- list()
    audit <- list()
    for (ln in chosen) {
      x <- tryCatch(SeuratObject::LayerData(a, layer = ln, fast = FALSE), error = function(e) stop("Cannot read layer ",
        ln, ": ", conditionMessage(e)))
      if (is.null(rownames(x)) || is.null(colnames(x)) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
        stop("Layer ", ln, " lacks unique feature/cell IDs")
      if (inherits(x, "sparseMatrix") && any(!is.finite(x@x)))
        stop("Layer ", ln, " contains non-finite counts")
      if (inherits(x, "sparseMatrix")) {
        if (any(x@x < 0))
          stop("Negative counts in layer ", ln)
        totals <- as.numeric(Matrix::colSums(x))
      } else {
        if (any(x < 0, na.rm = TRUE))
          stop("Negative counts in layer ", ln)
        totals <- as.numeric(base::colSums(x))
      }
      present <- intersect(genes, rownames(x))
      sub <- if (length(present))
        x[present, , drop = FALSE] else NULL
      if (!is.null(sub) && !inherits(sub, "sparseMatrix"))
        sub <- Matrix::Matrix(sub, sparse = TRUE)
      parts[[ln]] <- list(cells = colnames(x), nfeatures = nrow(x), present = present, panel = sub,
        total = totals)
      audit[[ln]] <- data.frame(layer = ln, n_cells = ncol(x), n_features = nrow(x), target_genes_present = length(present),
        example_cell = if (ncol(x))
          colnames(x)[1L] else NA_character_, UMI_sum = sum(totals), stringsAsFactors = FALSE)
      rm(x, sub)
      gc(verbose = FALSE)
    }
    allcells <- unlist(lapply(parts, `[[`, "cells"), use.names = FALSE)
    if (!length(allcells))
      stop("No cells in count layers")
    if (anyDuplicated(allcells))
      stop("Overlapping cell IDs across count layers; refusing arbitrary selection/merge")
    mat <- Matrix::Matrix(0, nrow = length(genes), ncol = length(allcells), sparse = TRUE)
    rownames(mat) <- genes
    colnames(mat) <- allcells
    avail <- matrix(FALSE, nrow = length(genes), ncol = length(allcells), dimnames = list(genes,
      allcells))
    total <- numeric(length(allcells))
    names(total) <- allcells
    by_layer <- character(length(allcells))
    names(by_layer) <- allcells
    offset <- 0L
    for (ln in names(parts)) {
      item <- parts[[ln]]
      at <- seq.int(offset + 1L, offset + length(item$cells))
      if (length(item$present)) {
        gi <- match(item$present, genes)
        mat[gi, at] <- item$panel
        avail[gi, at] <- TRUE
      }
      total[at] <- item$total
      by_layer[at] <- ln
      offset <- offset + length(item$cells)
    }
    audit <- do.call(rbind, audit)
    rownames(audit) <- NULL
    audit$assay <- selected_assay
    return(list(mat = mat, available = avail, total = total, by_layer = by_layer, audit = audit,
      nfeatures_union = NA_integer_, layer = paste0(selected_assay, "/", paste(chosen, collapse = ";")),
      assay_inventory = paste(anames, collapse = ";"), layer_inventory = paste(layers, collapse = ";")))
  }
  if (inherits(object, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
      stop("SummarizedExperiment not available")
    av <- SummarizedExperiment::assayNames(object)
    pref <- if (corrected)
      c("decontXcounts", "corrected", "counts") else c("counts", "raw")
    hit <- pref[pref %in% av]
    if (!length(hit))
      stop("No validated SCE count assay; found ", paste(av, collapse = ";"))
    x <- SummarizedExperiment::assay(object, hit[1L])
    if (is.null(rownames(x)) || is.null(colnames(x)) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
      stop("SCE assay has missing/duplicated IDs")
    if (inherits(x, "sparseMatrix")) {
      if (any(x@x < 0) || any(!is.finite(x@x)))
        stop("Invalid SCE counts")
      total <- as.numeric(Matrix::colSums(x))
    } else {
      if (any(x < 0, na.rm = TRUE))
        stop("Invalid SCE counts")
      total <- as.numeric(base::colSums(x))
    }
    present <- intersect(genes, rownames(x))
    mat <- Matrix::Matrix(0, nrow = length(genes), ncol = ncol(x), sparse = TRUE, dimnames = list(genes,
      colnames(x)))
    avail <- matrix(FALSE, length(genes), ncol(x), dimnames = list(genes, colnames(x)))
    if (length(present)) {
      mat[match(present, genes), ] <- x[present, , drop = FALSE]
      avail[match(present, genes), ] <- TRUE
    }
    names(total) <- colnames(x)
    aud <- data.frame(layer = hit[1L], n_cells = ncol(x), n_features = nrow(x), target_genes_present = length(present),
      example_cell = colnames(x)[1L], UMI_sum = sum(total), assay = "SummarizedExperiment", stringsAsFactors = FALSE)
    return(list(mat = mat, available = avail, total = total, by_layer = setNames(rep(hit[1L], ncol(x)),
      colnames(x)), audit = aud, nfeatures_union = nrow(x), layer = hit[1L], assay_inventory = paste(av,
      collapse = ";"), layer_inventory = paste(av, collapse = ";")))
  }
  stop("Unsupported assay object: ", paste(class(object), collapse = ";"))
}

get_barcode <- function(ids, metadata_barcode) {

  b <- trimws(as.character(metadata_barcode))
  b[!grepl("^[ACGT]{16}(-[0-9]+)?$", b)] <- NA_character_
  bare <- !is.na(b)
  b[bare] <- sub("-[0-9]+$", "", b[bare])

  miss <- which(is.na(b))
  for (i in miss) {
    hits <- regmatches(ids[i], gregexpr("(?<![ACGT])[ACGT]{16}(?![ACGT])", ids[i], perl = TRUE))[[1L]]
    if (length(hits) == 1L)
      b[i] <- hits[1L]
  }
  b
}
get_ids <- function(m, ids, origin, known_source) {

  original_col <- first_col(m, c("cell_id_original", "original_cell_id", "cellID_original", "raw_cell_id"))
  barcode_col <- first_col(m, c("barcode16", "barcode", "Barcode", "cell_barcode"))
  sample_col <- first_col(m, c("sample_id", "sample_id_source", "library_id", "library", "GSM", "sample"))
  source_col <- first_col(m, c("dataset", "source_dataset", "dataset_id", "source"))
  treatment_col <- first_col(m, c("treatment", "condition", "Treatment", "group"))
  stage_col <- first_col(m, c("time", "timepoint", "stage", "injury_stage", "dpi"))
  state_col <- first_col(m, c("myeloid_state"))
  original <- read_col(m, original_col)
  bar_meta <- read_col(m, barcode_col)
  source <- if (origin == "atlas")
    read_col(m, source_col) else rep(known_source, length(ids))
  sample <- read_col(m, sample_col)
  treatment <- read_col(m, treatment_col)

  data.frame(cell_id = ids, original_id = original, barcode = get_barcode(ids, bar_meta), sample = sample,
    source = source, treatment = treatment, stage = read_col(m, stage_col), state = read_col(m, state_col),
    nCount_RNA = numeric_col(m, "nCount_RNA"), nCount_before_DecontX = numeric_col(m, "nCount_RNA_before_DecontX"),
    decontX_contamination = numeric_col(m, "decontX_contamination"), decontX_removed_UMI = numeric_col(m,
      "decontX_removed_UMI"), stringsAsFactors = FALSE)
}
meaningful <- function(v) !is.na(v) & nzchar(trimws(as.character(v)))
key_sample_barcode <- function(df) {
  ok <- meaningful(df$sample) & meaningful(df$barcode)
  z <- rep(NA_character_, nrow(df))
  z[ok] <- paste(df$sample[ok], df$barcode[ok], sep = "::")
  z
}

link_candidates <- function(a, b, context) {
  rules <- list(exact_cell_id = list(a$cell_id, b$cell_id), a_original_to_b_cell_id = list(a$original_id,
    b$cell_id), a_cell_id_to_b_original = list(a$cell_id, b$original_id), original_id_to_original_id = list(a$original_id,
    b$original_id), explicit_sample_plus_barcode = list(key_sample_barcode(a), key_sample_barcode(b)))
  items <- list()
  k <- 0L
  for (rule in names(rules)) {
    aa <- as.character(rules[[rule]][[1L]])
    bb <- as.character(rules[[rule]][[2L]])
    aok <- meaningful(aa)
    bok <- meaningful(bb)
    va <- which(aok & !(duplicated(aa) | duplicated(aa, fromLast = TRUE)))
    vb <- which(bok & !(duplicated(bb) | duplicated(bb, fromLast = TRUE)))
    if (!length(va) || !length(vb))
      next
    i <- match(aa[va], bb[vb])
    keep <- which(!is.na(i))
    if (!length(keep))
      next
    k <- k + 1L
    items[[k]] <- data.frame(left = va[keep], right = vb[i[keep]], rule = rule, priority = match(rule,
      names(rules)), context = context, stringsAsFactors = FALSE)
  }
  if (!length(items))
    return(data.frame(left = integer(), right = integer(), rule = character(), priority = integer(),
      context = character()))
  do.call(rbind, items)
}

treatment_class <- function(x) {
  y <- tolower(trimws(as.character(x)))
  ans <- rep(NA_character_, length(y))
  ans[grepl("^(vehicle|veh|control|ctrl|untreated)$", y)] <- "Control"
  ans[grepl("plx5622|plx-5622|^plx$", y)] <- "PLX5622"
  ans
}
select_verified <- function(a, b, candidates, expected_treatment = NA_character_) {
  empty <- candidates[FALSE, , drop = FALSE]
  empty$status <- character()
  empty$reason <- character()
  if (!nrow(candidates))
    return(list(good = empty, all = empty))
  candidates$status <- "candidate"
  candidates$reason <- ""

  lr <- unique(candidates[c("left", "right")])
  nl <- table(lr$left)
  nr <- table(lr$right)
  conflict <- as.character(candidates$left) %in% names(nl[nl > 1L]) | as.character(candidates$right) %in%
    names(nr[nr > 1L])
  candidates$status[conflict] <- "REJECT"
  candidates$reason[conflict] <- "conflicting mappings across ID rules"
  for (i in which(!conflict)) {
    ai <- candidates$left[i]
    bi <- candidates$right[i]
    a_sample <- a$sample[ai]
    b_sample <- b$sample[bi]
    if (meaningful(a_sample) && meaningful(b_sample) && a_sample != b_sample) {
      candidates$status[i] <- "REJECT"
      candidates$reason[i] <- "explicit sample/library labels differ"
      next
    }
    ac <- treatment_class(a$treatment[ai])
    bc <- treatment_class(b$treatment[bi])
    if (!is.na(ac) && !is.na(bc) && ac != bc) {
      candidates$status[i] <- "REJECT"
      candidates$reason[i] <- "explicit raw/corrected or atlas/source treatment mismatch"
      next
    }
    if (!is.na(expected_treatment)) {
      if (is.na(ac) || ac != expected_treatment) {
        candidates$status[i] <- "REJECT"
        candidates$reason[i] <- "atlas treatment unknown/mismatched"
        next
      }
      if (!is.na(bc) && bc != expected_treatment) {
        candidates$status[i] <- "REJECT"
        candidates$reason[i] <- "source-pair treatment contradicts input directory"
        next
      }
    }
    candidates$status[i] <- "ACCEPT"
    candidates$reason[i] <- "source, uniqueness, and available sample/treatment checks passed"
  }
  good <- candidates[candidates$status == "ACCEPT", , drop = FALSE]
  good <- good[order(good$priority), , drop = FALSE]
  good <- good[!duplicated(good[c("left", "right")]), , drop = FALSE]

  if (anyDuplicated(good$left) || anyDuplicated(good$right))
    stop("Internal ambiguity in purported verified crosswalk")
  list(good = good, all = candidates)
}

go <- utils::read.csv(input[["myelin_leading"]], check.names = FALSE)
stopifnot(all(c("GO_ID", "gene") %in% names(go)))
edge <- sort(unique(as.character(go$gene[as.character(go$GO_ID) == "GO:0043209"])))
edge <- edge[meaningful(edge)]
if (!length(edge)) stop("GO:0043209 leading edge is absent from Step 102 CSV")
canonical <- c("Mbp", "Plp1", "Mog", "Mobp", "Mag", "Cnp", "Mal", "Pmp22")
oligo <- c("Olig1", "Olig2", "Sox10", "Myrf")
phago <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Tyrobp", "Cd68", "Cst7", "Gpnmb")
genes <- sort(unique(c(edge, canonical, oligo, phago)))
write_csv(data.frame(gene = genes, myelin_GO_leading_edge = genes %in% edge, canonical_myelin = genes %in%
  canonical, oligodendrocyte_context = genes %in% oligo, phagocytosis_context = genes %in% phago),
  "01_138_target_gene_manifest_v01.csv")

inventory <- list()
patterns <- list()
pair_qc <- list()
link_qc <- list()
layer_audit <- list()
layer_qc <- list()
all_crosswalk <- list()
all_rejections <- list()
all_gene <- list()
summary_flags <- character()

load_stage <- function(file, tag, known_source, corrected) {
  obj <- readRDS(file)
  assay <- assay_counts(obj, corrected)
  m <- assay$mat
  if (is.null(colnames(m)) || anyDuplicated(colnames(m)))
    stop(tag, " lacks unique count-layer cell IDs")
  meta <- meta_of(obj, colnames(m))
  metadata_source_field <- first_col(meta, c("dataset", "source_dataset", "dataset_id", "source"))
  declared_sources <- unique(read_col(meta, metadata_source_field))
  declared_sources <- declared_sources[meaningful(declared_sources)]
  if (tag != "atlas" && length(declared_sources)) {
    known_accessions <- c("17702045", "GSE162610", "GSE196928")
    recognized <- intersect(declared_sources, known_accessions)
    if (length(recognized) && !all(recognized == known_source))
      stop("Stage file metadata conflicts with provenance: ", paste(recognized, collapse = ";"))
  }
  ids <- get_ids(meta, colnames(m), if (tag == "atlas")
    "atlas" else tag, known_source)
  ids$matrix_total_UMI <- as.numeric(assay$total[ids$cell_id])
  ids$count_layer <- as.character(assay$by_layer[ids$cell_id])
  info <- data.frame(stage = tag, file = file, n_cells = ncol(m), n_features = NA_integer_, n_gene_panel_features = sum(rowSums(assay$available) >
    0), n_ID_original_present = sum(meaningful(ids$original_id)), n_barcode_present = sum(meaningful(ids$barcode)),
    n_sample_present = sum(meaningful(ids$sample)), n_treatment_present = sum(meaningful(ids$treatment)),
    n_count_layers = nrow(assay$audit), count_layer = assay$layer, assay_inventory = assay$assay_inventory,
    layer_inventory = assay$layer_inventory, declared_source_metadata = paste(declared_sources, collapse = ";"),
    metadata_field_names = paste(names(meta), collapse = ";"), stringsAsFactors = FALSE)
  layer_info <- assay$audit
  layer_info$stage <- tag
  layer_info$source <- known_source
  pattern <- data.frame(stage = tag, label = c("cell_id", "original_id", "sample", "barcode", "treatment"),
    n_available = c(length(ids$cell_id), sum(meaningful(ids$original_id)), sum(meaningful(ids$sample)),
      sum(meaningful(ids$barcode)), sum(meaningful(ids$treatment))), n_distinct = c(length(unique(ids$cell_id)),
      length(unique(ids$original_id[meaningful(ids$original_id)])), length(unique(ids$sample[meaningful(ids$sample)])),
      length(unique(ids$barcode[meaningful(ids$barcode)])), length(unique(ids$treatment[meaningful(ids$treatment)]))),
    example = c(paste(head(ids$cell_id, 3), collapse = " | "), paste(head(ids$original_id[meaningful(ids$original_id)],
      3), collapse = " | "), paste(head(unique(ids$sample[meaningful(ids$sample)]), 3), collapse = " | "),
      paste(head(ids$barcode[meaningful(ids$barcode)], 3), collapse = " | "), paste(head(unique(ids$treatment[meaningful(ids$treatment)]),
        3), collapse = " | ")), stringsAsFactors = FALSE)
  rm(obj, meta)
  gc(verbose = FALSE)
  list(ids = ids, counts = m, available = assay$available, info = info, pattern = pattern, layers = layer_info)
}

atlas <- load_stage(input[["atlas"]], "atlas", NA_character_, FALSE)
if (!all(c("state", "source", "cell_id") %in% names(atlas$ids))) stop("Atlas missing identity fields")
if (!any(meaningful(atlas$ids$state))) stop("Step 20 lacks myeloid_state; do not infer focal state")
if (!all(meaningful(atlas$ids$source))) stop("Atlas dataset/source labels incomplete; cannot link safely")
atlas$ids$focus <- atlas$ids$state %in% c("Clec7a+ lipid-handling Microglia", "Clec7a+ lipid-associated Microglia")
if (!any(atlas$ids$focus)) stop("Atlas focal annotation not found; cannot infer from marker expression")
inventory[["atlas"]] <- atlas$info
layer_audit[["atlas"]] <- atlas$layers
patterns[["atlas"]] <- atlas$pattern
atlas$counts <- NULL
atlas$available <- NULL

as <- atlas$ids
as$group <- ifelse(as$focus, "Clec7a_focus", "Other_myeloid")
meta_groups <- split(seq_len(nrow(as)), paste(as$source, as$group, sep = " | "))
meta_summary <- do.call(rbind, lapply(meta_groups, function(ii) {
  v <- as$decontX_contamination[ii]
  b <- as$nCount_before_DecontX[ii]
  after <- as$nCount_RNA[ii]
  rem <- as$decontX_removed_UMI[ii]
  data.frame(source = as$source[ii[1L]], cell_group = as$group[ii[1L]], n_cells = length(ii), n_with_contamination_score = sum(is.finite(v)),
    median_contamination_score = if (any(is.finite(v)))
      stats::median(v, na.rm = TRUE) else NA_real_, n_with_before_after_UMI_metadata = sum(is.finite(b) & is.finite(after)), n_with_removed_UMI_metadata = sum(is.finite(rem)),
    n_negative_metadata_removed_UMI = sum(is.finite(rem) & rem < 0), n_before_less_than_after = sum(is.finite(b) &
      is.finite(after) & b < after), interpretation = "Inherited QC metadata; not a gene-specific ambient-RNA estimate",
    stringsAsFactors = FALSE)
}))
write_csv(meta_summary, "11_138_Step20_inherited_DecontX_metadata_context_v01.csv")

source_groups <- c("17702045", "GSE162610", "GSE196928_Control", "GSE196928_PLX5622")
for (src in source_groups) {
  raw_key <- paste0("raw_", src)
  corr_key <- paste0("corr_", src)
  if (!file.exists(input[[raw_key]]) || !file.exists(input[[corr_key]])) {
    pair_qc[[src]] <- data.frame(source = src, raw_cells = NA_integer_, corrected_cells = NA_integer_,
      candidate_rows = 0L, n_verified_raw_corr = 0L, n_exact_ID_pairs = 0L, n_original_alias_pairs = 0L,
      n_sample_barcode_pairs = 0L, n_rejected_candidate_rows = 0L, status = "SKIP_missing_RDS",
      stringsAsFactors = FALSE)
    link_qc[[src]] <- data.frame(source = src, atlas_cells_in_source = NA_integer_, atlas_focus_in_source = NA_integer_,
      n_bridge_pairs = 0L, candidate_rows = 0L, n_verified_atlas_links = 0L, n_verified_atlas_focus = 0L,
      n_rejected_candidate_rows = 0L, status = "SKIP_missing_RDS", stringsAsFactors = FALSE)
    next
  }
  parent_src <- sub("_(Control|PLX5622)$", "", src)
  intended_treatment <- if (grepl("_Control$", src))
    "Control" else if (grepl("_PLX5622$", src))
    "PLX5622" else NA_character_
  result <- tryCatch((function() {

    raw <- load_stage(input[[raw_key]], paste0(src, "_raw"), parent_src, FALSE)
    corr <- load_stage(input[[corr_key]], paste0(src, "_corrected"), parent_src, TRUE)
    inventory[[paste0(src, "_raw")]] <<- raw$info
    layer_audit[[paste0(src, "_raw")]] <<- raw$layers
    inventory[[paste0(src, "_corrected")]] <<- corr$info
    layer_audit[[paste0(src, "_corrected")]] <<- corr$layers
    patterns[[paste0(src, "_raw")]] <<- raw$pattern
    patterns[[paste0(src, "_corrected")]] <<- corr$pattern
    cands <- link_candidates(raw$ids, corr$ids, paste0(src, ":raw_to_corrected"))
    check <- select_verified(raw$ids, corr$ids, cands)
    verified <- check$good
    if (nrow(check$all)) {
      rejected <- check$all[check$all$status == "REJECT", , drop = FALSE]
      if (nrow(rejected)) {
        rejected$left_cell_id <- raw$ids$cell_id[rejected$left]
        rejected$right_cell_id <- corr$ids$cell_id[rejected$right]
        all_rejections[[paste0(src, "_raw_corr")]] <<- rejected
      }
    }
    pair_qc[[src]] <<- data.frame(source = src, raw_cells = nrow(raw$ids), corrected_cells = nrow(corr$ids),
      candidate_rows = nrow(cands), n_verified_raw_corr = nrow(verified), n_exact_ID_pairs = sum(verified$rule ==
        "exact_cell_id"), n_original_alias_pairs = sum(grepl("original", verified$rule)), n_sample_barcode_pairs = sum(verified$rule ==
        "explicit_sample_plus_barcode"), n_rejected_candidate_rows = sum(check$all$status ==
        "REJECT"), status = if (nrow(verified))
        "VERIFIED_1_TO_1" else "NO_VERIFIABLE_RAW_CORR_PAIRS", stringsAsFactors = FALSE)
    if (!nrow(verified)) {
      link_qc[[src]] <<- data.frame(source = src, atlas_cells_in_source = sum(atlas$ids$source ==
        parent_src), atlas_focus_in_source = sum(atlas$ids$source == parent_src & atlas$ids$focus),
        n_bridge_pairs = 0L, candidate_rows = 0L, n_verified_atlas_links = 0L, n_verified_atlas_focus = 0L,
        n_rejected_candidate_rows = 0L, status = "NO_RAW_CORR_BRIDGE", stringsAsFactors = FALSE)
      return(invisible(NULL))
    }

    bridge <- corr$ids[verified$right, , drop = FALSE]
    bridge$original_id <- raw$ids$cell_id[verified$left]

    bridge$sample <- ifelse(meaningful(bridge$sample), bridge$sample, raw$ids$sample[verified$left])
    bridge$barcode <- ifelse(meaningful(bridge$barcode), bridge$barcode, raw$ids$barcode[verified$left])

    bridge$treatment <- ifelse(meaningful(bridge$treatment), bridge$treatment, raw$ids$treatment[verified$left])
    aidx <- which(atlas$ids$source == parent_src)
    a <- atlas$ids[aidx, , drop = FALSE]

    atlas_candidates <- link_candidates(a, bridge, paste0(src, ":atlas_to_bridge"))
    match_result <- select_verified(a, bridge, atlas_candidates, intended_treatment)
    matched <- match_result$good
    if (nrow(match_result$all)) {
      rejected <- match_result$all[match_result$all$status == "REJECT", , drop = FALSE]
      if (nrow(rejected)) {
        rejected$left_cell_id <- a$cell_id[rejected$left]
        rejected$right_cell_id <- bridge$cell_id[rejected$right]
        all_rejections[[paste0(src, "_atlas_bridge")]] <<- rejected
      }
    }
    link_qc[[src]] <<- data.frame(source = src, atlas_cells_in_source = nrow(a), atlas_focus_in_source = sum(a$focus),
      n_bridge_pairs = nrow(bridge), candidate_rows = nrow(atlas_candidates), n_verified_atlas_links = nrow(matched),
      n_verified_atlas_focus = if (nrow(matched))
        sum(a$focus[matched$left]) else 0L, n_rejected_candidate_rows = sum(match_result$all$status == "REJECT"), status = if (nrow(matched))
        "VERIFIED_ATLAS_BRIDGE" else "NO_VERIFIABLE_ATLAS_LINKS", stringsAsFactors = FALSE)
    if (!nrow(matched))
      return(invisible(NULL))
    arow <- a[matched$left, , drop = FALSE]
    pi <- matched$right
    ri <- verified$left[pi]
    ci <- verified$right[pi]
    cw <- data.frame(source = src, atlas_cell_id = arow$cell_id, atlas_original_id = arow$original_id,
      atlas_sample = arow$sample, atlas_barcode16 = arow$barcode, atlas_treatment = arow$treatment,
      atlas_stage = arow$stage, frozen_state = arow$state, is_Clec7a_focus = arow$focus, raw_cell_id = raw$ids$cell_id[ri],
      corrected_cell_id = corr$ids$cell_id[ci], raw_sample = raw$ids$sample[ri], corrected_sample = corr$ids$sample[ci],
      atlas_RNA_total_UMI = arow$matrix_total_UMI, raw_RNA_total_UMI = raw$ids$matrix_total_UMI[ri],
      corrected_RNA_total_UMI = corr$ids$matrix_total_UMI[ci], atlas_metadata_pre_DecontX_UMI = arow$nCount_before_DecontX,
      atlas_metadata_removed_UMI = arow$decontX_removed_UMI, raw_matches_atlas_pre_metadata = ifelse(is.finite(arow$nCount_before_DecontX),
        raw$ids$matrix_total_UMI[ri] == arow$nCount_before_DecontX, NA), corrected_matches_atlas_counts = ifelse(is.finite(arow$matrix_total_UMI),
        corr$ids$matrix_total_UMI[ci] == arow$matrix_total_UMI, NA), raw_corrected_link_rule = verified$rule[pi],
      atlas_to_pair_link_rule = matched$rule, sample_id_verified_both_sides = meaningful(arow$sample) &
        meaningful(bridge$sample[pi]) & arow$sample == bridge$sample[pi], interpretation = "verified cell-level linkage; not a proof of transcript origin",
      stringsAsFactors = FALSE)

    observed_count_ref <- (!is.na(cw$raw_matches_atlas_pre_metadata)) + (!is.na(cw$corrected_matches_atlas_counts))
    bad_reference <- (!is.na(cw$raw_matches_atlas_pre_metadata) & !cw$raw_matches_atlas_pre_metadata) |
      (!is.na(cw$corrected_matches_atlas_counts) & !cw$corrected_matches_atlas_counts)
    cw$has_UMI_reference <- observed_count_ref >= 1L
    cw$UMI_reference_consistent <- cw$has_UMI_reference & !bad_reference
    cw$corrected_total_not_higher <- cw$corrected_RNA_total_UMI <= cw$raw_RNA_total_UMI
    cw$count_decreased <- cw$corrected_RNA_total_UMI < cw$raw_RNA_total_UMI
    cw$eligible_for_expression <- cw$UMI_reference_consistent & cw$corrected_total_not_higher
    cw$eligibility_reason <- ifelse(!cw$has_UMI_reference, "NO_INDEPENDENT_UMI_REFERENCE", ifelse(bad_reference,
      "ATLAS_UMI_MISMATCH", ifelse(!cw$corrected_total_not_higher, "CORRECTED_TOTAL_EXCEEDS_RAW",
        "VALID")))
    branch_correction_observed <- any(cw$count_decreased & cw$eligible_for_expression)
    if (!branch_correction_observed) {
      cw$eligible_for_expression[] <- FALSE
      cw$eligibility_reason[cw$eligibility_reason == "VALID"] <- "NO_CORRECTION_OBSERVED"
    }
    all_crosswalk[[src]] <<- cw
    layer_qc[[src]] <<- data.frame(source = src, n_atlas_linked = nrow(cw), n_with_UMI_reference = sum(cw$has_UMI_reference),
      n_with_reference_consistent = sum(cw$UMI_reference_consistent), n_corrected_higher_than_raw = sum(!cw$corrected_total_not_higher),
      n_cells_with_count_decrease = sum(cw$count_decreased), n_eligible_for_expression = sum(cw$eligible_for_expression),
      n_eligible_focus = sum(cw$eligible_for_expression & cw$is_Clec7a_focus), correction_observed = branch_correction_observed,
      status = if (any(cw$eligible_for_expression))
        "VALIDATED" else "NOT_VALIDATED", stringsAsFactors = FALSE)
    good_idx <- which(cw$eligible_for_expression)
    if (!length(good_idx)) {
      summary_flags <<- c(summary_flags, paste(src, "no UMI-validated cell pairs"))
      return(invisible(NULL))
    }
    arow <- arow[good_idx, , drop = FALSE]
    pi <- pi[good_idx]
    ri <- ri[good_idx]
    ci <- ci[good_idx]
    both_genes <- intersect(rownames(raw$counts), rownames(corr$counts))
    if (!length(both_genes)) {
      summary_flags <<- c(summary_flags, paste(src, "no shared target-gene labels"))
      return(invisible(NULL))
    }
    group <- paste(arow$stage, ifelse(arow$focus, "Clec7a_focus", "Other_myeloid"), sep = " || ")
    idx_group <- split(seq_along(ri), group)
    results <- lapply(idx_group, function(q) {
      do.call(rbind, lapply(both_genes, function(g) {

        valid_gene <- raw$available[g, ri[q]] & corr$available[g, ci[q]]
        eligible <- q[valid_gene]
        if (!length(eligible))
          return(NULL)
        rr <- as.numeric(raw$counts[g, ri[eligible]])
        cc <- as.numeric(corr$counts[g, ci[eligible]])
        sr <- sum(rr)
        sc <- sum(cc)
        data.frame(source = src, stage = arow$stage[q[1L]], cell_group = ifelse(arow$focus[q[1L]],
          "Clec7a_focus", "Other_myeloid"), gene = g, myelin_GO_leading_edge = g %in% edge, n_verified_paired_cells_in_group = length(q),
          n_gene_assayed_both_stages = length(eligible), n_cells_missing_gene_in_either_layer = length(q) -
          length(eligible), raw_detected = sum(rr > 0), corrected_detected = sum(cc > 0), raw_detection_fraction = mean(rr >
          0), corrected_detection_fraction = mean(cc > 0), raw_sum_UMI = sr, corrected_sum_UMI = sc,
          n_cells_with_count_changed = sum(rr != cc), n_cells_with_corrected_count_higher_than_raw = sum(cc >
          rr), relative_UMI_removed = if (sr > 0)
          (sr - sc)/sr else NA_real_, mean_UMI_change_corrected_minus_raw = mean(cc - rr), raw_layer_set = paste(unique(raw$ids$count_layer[ri[eligible]]),
          collapse = ";"), corrected_layer_set = paste(unique(corr$ids$count_layer[ci[eligible]]),
          collapse = ";"), interpretation = "Descriptive UMI sensitivity; does not prove cellular RNA origin",
          stringsAsFactors = FALSE)
      }))
    })
    all_gene[[src]] <<- do.call(rbind, results)
    invisible(NULL)
  })(), error = function(e) {
    msg <- conditionMessage(e)
    if (is.null(pair_qc[[src]]))
      pair_qc[[src]] <<- data.frame(source = src, raw_cells = NA_integer_, corrected_cells = NA_integer_,
        candidate_rows = 0L, n_verified_raw_corr = 0L, n_exact_ID_pairs = 0L, n_original_alias_pairs = 0L,
        n_sample_barcode_pairs = 0L, n_rejected_candidate_rows = 0L, status = paste("ERROR",
          msg), stringsAsFactors = FALSE)
    if (is.null(link_qc[[src]]))
      link_qc[[src]] <<- data.frame(source = src, atlas_cells_in_source = sum(atlas$ids$source ==
        parent_src), atlas_focus_in_source = sum(atlas$ids$source == parent_src & atlas$ids$focus),
        n_bridge_pairs = 0L, candidate_rows = 0L, n_verified_atlas_links = 0L, n_verified_atlas_focus = 0L,
        n_rejected_candidate_rows = 0L, status = paste("ERROR", msg), stringsAsFactors = FALSE)
    summary_flags <<- c(summary_flags, paste(src, msg))
    NULL
  })
  rm(result)
  gc(verbose = FALSE)
}

combine <- function(xs, template = NULL) {
  if (length(xs))
    return(do.call(rbind, xs))
  if (!is.null(template))
    return(template)
  data.frame()
}
write_csv(combine(inventory), "02_138_object_metadata_and_assay_inventory_v01.csv")
write_csv(combine(patterns), "04_138_cell_ID_pattern_diagnostics_v01.csv")
write_csv(combine(layer_audit, data.frame(stage = character(), source = character(), assay = character(),
  layer = character(), n_cells = integer(), n_features = integer(), target_genes_present = integer(),
  example_cell = character(), UMI_sum = numeric())), "03_138_split_counts_layer_inventory_v01.csv")
pr <- combine(pair_qc)
lq <- combine(link_qc)
write_csv(pr, "05_138_raw_vs_corrected_pairing_QC_v01.csv")
write_csv(lq, "06_138_annotated_atlas_linkage_QC_v01.csv")
crosswalk <- combine(all_crosswalk, data.frame(source = character(), atlas_cell_id = character(), raw_cell_id = character(),
  corrected_cell_id = character(), is_Clec7a_focus = logical()))
write_csv(crosswalk, "07_138_verified_cell_ID_crosswalk_v01.csv")
write_csv(combine(layer_qc, data.frame(source = character(), n_atlas_linked = integer(), n_with_UMI_reference = integer(),
  n_with_reference_consistent = integer(), n_corrected_higher_than_raw = integer(), n_cells_with_count_decrease = integer(),
  n_eligible_for_expression = integer(), n_eligible_focus = integer(), correction_observed = logical(),
  status = character())), "08_138_count_provenance_validation_QC_v01.csv")
rejects <- combine(all_rejections, data.frame(context = character(), left = integer(), right = integer(),
  rule = character(), priority = integer(), status = character(), reason = character(), left_cell_id = character(),
  right_cell_id = character()))
write_csv(rejects, "09_138_ambiguous_or_rejected_ID_candidates_v01.csv")
gene_data <- combine(all_gene, data.frame(source = character(), stage = character(), cell_group = character(),
  gene = character(), n_verified_paired_cells_in_group = integer(), n_gene_assayed_both_stages = integer(),
  raw_detected = integer(), corrected_detected = integer()))
write_csv(gene_data, "10_138_myelin_gene_raw_vs_DecontX_validated_pairs_v01.csv")
limit <- data.frame(assessment = c("Raw-to-corrected ID pairs", "Atlas-to-correction linkage", "Gene-specific DecontX sensitivity",
  "Inherited Step20 contamination metadata", "Phagocytosed material", "Endogenous myelin expression",
  "Remyelination / repair"), conclusion = c("Only unique exact/original/sample+barcode pairs accepted; all conflicts excluded",
  "Must match exact source and treatment (GSE196928); missing identification is not imputed", "Only if ID and independent UMI references both pass, correction is observed, and gene is present in both layers",
  "Provides per-cell QC context but does not assign gene-specific ambient origin", "Cannot be directly identified with droplet RNA counts",
  "Cannot be distinguished conclusively from engulfed/ambient RNA with these data", "Not shown by GO term enrichment, co-detection, or inferred network effects"),
  stringsAsFactors = FALSE)
write_csv(limit, "12_138_interpretation_and_limitations_v01.csv")
verified_total <- if (nrow(lq)) sum(lq$n_verified_atlas_links, na.rm = TRUE) else 0L
focus_total <- if (nrow(lq)) sum(lq$n_verified_atlas_focus, na.rm = TRUE) else 0L
qc <- data.frame(check = c("Step20 focal state present", "GO:0043209 leading genes loaded", "Atlas source labels complete",
  "Raw/corrected sources processed", "Verified raw/corrected cell links", "Verified atlas-to-source links",
  "Verified focal Clec7a+ cell links", "Gene-specific paired comparison available", "Cell origin determined from expression",
  "Data and source packaging"), status = c(if (any(atlas$ids$focus)) "PASS" else "FAIL", if (length(edge)) "PASS" else "FAIL",
  if (all(meaningful(atlas$ids$source))) "PASS" else "FAIL", paste(sum(!grepl("SKIP|ERROR", pr$status)),
    "of", nrow(pr)), as.character(sum(pr$n_verified_raw_corr, na.rm = TRUE)), as.character(verified_total),
  as.character(focus_total), if (nrow(gene_data)) "YES (descriptive only)" else "NO (do not infer)",
  "NO - not identifiable", "Direct D:/zz paths, no ZIP"), interpretation = c("Step20 annotation", "Step102 existing CSV, no new GO-GSEA",
  "No source assigned from expression profiles", "Inspect source-specific QC table", "May include non-focal and/or excluded cells",
  "No fuzzy/position matching", "ID match alone does not qualify for gene comparison; consult 08_138 provenance QC",
  "No change in UMI is not proof of genuine cellular transcription", "Requires orthogonal evidence",
  "No archive dependency"), stringsAsFactors = FALSE)
write_csv(qc, "13_138_FINAL_QC_v01.csv")
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(out, "14_138_REPORT_v01.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "15_138_sessionInfo_v01.txt"), useBytes = TRUE)
write_csv(data.frame(filename = basename(list.files(tabs, full.names = TRUE)), bytes = as.numeric(file.info(list.files(tabs,
  full.names = TRUE))$size), stringsAsFactors = FALSE), "00_138_OUTPUT_INDEX_v01.csv")
