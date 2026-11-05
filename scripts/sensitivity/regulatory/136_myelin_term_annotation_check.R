source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, width = 170)
project_root <- PROJECT_ROOT
step_id <- "136"
version <- "v01"
output_dir <- file.path(project_root, "results", "136_myelin_sheath_annotation_audit")
table_dir <- file.path(output_dir, "tables")
log_dir <- file.path(output_dir, "logs")
for (d in c(output_dir, table_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "00_136_run_log_v01.log")
writeLines(paste("STEP 136", Sys.time(), "PROJECT", project_root), log_file)
out_csv <- function(df, filename) {
  dest <- file.path(table_dir, filename)
  utils::write.csv(df, dest, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(dest)
}
path <- function(...) file.path(project_root, ...)
input_files <- c(GO_CC = path("results", "102_Zfp148_GO_GSEA_fullrank_v01", "03_102_Zfp148_primary_filtered_GO_CC_FDR05_v01.csv"),
  GO_leading = path("results", "102_Zfp148_GO_GSEA_fullrank_v01", "07_102_Zfp148_primary_GO_leading_edge_FDR05_v01.csv"),
  GO_members = path("results", "103_Figure6G_Zfp148_GO_GSEA_ridgeplot_v03", "02_103_Zfp148_full_GO_member_genes_v03.csv"),
  focus_atlas = path("results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds"),
  raw_17702045 = path("results", "02_QC", "02_QC_17702045_v01.rds"), corrected_17702045 = path("results",
    "03_DecontX", "03_DecontX_17702045_v01.rds"), raw_GSE162610 = path("results", "02_QC", "02_QC_GSE162610_v01.rds"),
  corrected_GSE162610 = path("results", "03_DecontX", "03_DecontX_GSE162610_v01.rds"), raw_GSE196928_Control = path("results",
    "02_QC", "02_QC_GSE196928_Control_v01.rds"), corrected_GSE196928_Control = path("results", "03_DecontX",
    "03_DecontX_GSE196928_Control_v01.rds"))
manifest <- data.frame(input_id = names(input_files), path = unname(input_files), exists = file.exists(input_files),
  stringsAsFactors = FALSE)
out_csv(manifest, "00_136_input_provenance_v01.csv")
required <- c("GO_CC", "GO_leading", "GO_members", "focus_atlas")
if (!all(manifest$exists[match(required, manifest$input_id)])) stop("Missing required input(s): ", paste(required[!manifest$exists[match(required,
  manifest$input_id)]], collapse = ", "), " | See 00_136_input_provenance_v01.csv")
for (pkg in c("Matrix", "SeuratObject")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Required installed R package missing: ", pkg, ". Install it separately; this script never installs packages.")
}

cc <- utils::read.csv(input_files[["GO_CC"]], check.names = FALSE)
lead <- utils::read.csv(input_files[["GO_leading"]], check.names = FALSE)
members <- utils::read.csv(input_files[["GO_members"]], check.names = FALSE)
need_cols <- function(df, cols, which) {
  missed <- setdiff(cols, names(df))
  if (length(missed))
    stop(which, " is missing columns: ", paste(missed, collapse = ", "))
}
need_cols(cc, c("ID", "Description", "NES", "p.adjust", "core_enrichment"), "GO_CC")
need_cols(lead, c("GO_ID", "GO_term", "gene", "ENTREZID", "perturbation_Z", "Ontology"), "GO_leading")
need_cols(members, c("GO_ID", "gene", "ENTREZID", "Z"), "GO_members")
go_id <- "GO:0043209"
my_cc <- cc[as.character(cc$ID) == go_id, , drop = FALSE]
my_le <- lead[as.character(lead$GO_ID) == go_id, , drop = FALSE]
my_members <- members[as.character(members$GO_ID) == go_id, , drop = FALSE]
if (nrow(my_cc) != 1L || nrow(my_le) == 0L || nrow(my_members) == 0L) stop("Myelin sheath GO:0043209 absent/inconsistent in Step 102/103 inputs.")
my_le$gene <- trimws(as.character(my_le$gene))
my_members$gene <- trimws(as.character(my_members$gene))
my_le_names <- sort(unique(my_le$gene[nzchar(my_le$gene)]))
my_member_names <- sort(unique(my_members$gene[nzchar(my_members$gene)]))
if (length(setdiff(my_le_names, my_member_names))) stop("GO myelin leading edge not a subset of Step 103 membership. Check upstream inputs.")
my_term <- data.frame(GO_ID = go_id, term = as.character(my_cc$Description[1]), ontology = "CC", NES = as.numeric(my_cc$NES[1]),
  BH_FDR = as.numeric(my_cc$p.adjust[1]), annotated_members_in_rank = length(my_member_names), n_leading_edge = length(my_le_names),
  network_Z_not_expression = TRUE, source_GO_CC = basename(input_files[["GO_CC"]]), source_GO_members = basename(input_files[["GO_members"]]),
  stringsAsFactors = FALSE)
out_csv(my_term, "01_136_myelin_GO_term_definition_v01.csv")

my_gene <- my_members[!duplicated(my_members$gene), intersect(c("gene", "ENTREZID", "Z", "GSEA_rank"),
  names(my_members)), drop = FALSE]
my_gene$is_GSEA_leading_edge <- my_gene$gene %in% my_le_names
my_gene$gene_role <- ifelse(my_gene$is_GSEA_leading_edge, "GO leading edge (network Z rank)", "Other ranked GO member")
my_gene$Z <- suppressWarnings(as.numeric(my_gene$Z))
my_gene <- my_gene[order(!my_gene$is_GSEA_leading_edge, -my_gene$Z), , drop = FALSE]
out_csv(my_gene, "02_136_myelin_full_members_and_leading_edge_v01.csv")

leading_sets <- split(as.character(lead$gene), as.character(lead$GO_ID))
leading_sets <- lapply(leading_sets, function(x) unique(x[!is.na(x) & nzchar(x)]))
annotation_overlap <- lapply(seq_len(nrow(cc)), function(i) {
  id <- as.character(cc$ID[i])
  genes <- leading_sets[[id]]
  if (is.null(genes))
    genes <- character()
  shared <- intersect(my_le_names, genes)
  union <- union(my_le_names, genes)
  data.frame(myelin_GO_ID = go_id, other_GO_ID = id, other_GO_term = as.character(cc$Description[i]),
    other_GO_NES = as.numeric(cc$NES[i]), other_GO_FDR = as.numeric(cc$p.adjust[i]), n_myelin_leading = length(my_le_names),
    n_other_leading = length(genes), n_shared_leading = length(shared), jaccard_leading = if (length(union))
      length(shared)/length(union) else NA_real_, shared_leading_genes = paste(sort(shared), collapse = ";"), overlap_basis = "leading edge / GO CC significant terms",
    stringsAsFactors = FALSE)
})
annotation_overlap <- do.call(rbind, annotation_overlap)
annotation_overlap <- annotation_overlap[order(-annotation_overlap$n_shared_leading, -annotation_overlap$jaccard_leading),
  , drop = FALSE]
out_csv(annotation_overlap, "03_136_myelin_shared_GO_annotation_v01.csv")

canonical_myelin <- c("Mbp", "Plp1", "Mog", "Mobp", "Mag", "Cnp", "Mal", "Pmp22")
oligolineage <- c("Olig1", "Olig2", "Sox10", "Myrf", "Gpr17", "Pdgfra")
phago_context <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Tyrobp", "Cd68", "Cst7", "Gpnmb", "Ctsd", "Lamp1",
  "Axl")
broad_annotation <- c("Gapdh", "Actb", "Eno1", "Prdx1", "Gfap")
markers <- sort(unique(c(my_member_names, canonical_myelin, oligolineage, phago_context, broad_annotation)))
markers_df <- data.frame(gene = markers, is_myelin_GO_member = markers %in% my_member_names, is_myelin_GO_leading_edge = markers %in%
  my_le_names, canonical_myelin_transcript_marker = markers %in% canonical_myelin, oligodendrocyte_lineage_context_marker = markers %in%
  oligolineage, phagocytosis_lipid_handling_context_marker = markers %in% phago_context, broad_or_non_specific_GO_context = markers %in%
  broad_annotation, stringsAsFactors = FALSE)
out_csv(markers_df, "04_136_myelin_marker_and_annotation_manifest_v01.csv")

get_matrix <- function(obj, kind = "counts") {
  if (inherits(obj, "Seurat")) {
    if (!"RNA" %in% names(obj@assays))
      stop("Seurat RNA assay absent")
    a <- obj[["RNA"]]
    x <- tryCatch(SeuratObject::LayerData(a, layer = "counts"), error = function(e) NULL)
    if (is.null(x))
      x <- tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts"), error = function(e) NULL)
    if (is.null(x))
      stop("No unambiguous RNA counts layer in Seurat object")
    return(list(matrix = x, layer = "RNA/counts"))
  }
  if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
      stop("SummarizedExperiment required to read this object")
    aa <- SummarizedExperiment::assayNames(obj)
    order <- if (kind == "corrected")
      c("decontXcounts", "corrected", "counts") else c("counts", "raw")
    chosen <- order[order %in% aa][1]
    if (is.na(chosen))
      stop("No supported count assay: ", paste(aa, collapse = ","))
    return(list(matrix = SummarizedExperiment::assay(obj, chosen), layer = chosen))
  }
  if (is.matrix(obj) || inherits(obj, "Matrix"))
    return(list(matrix = obj, layer = "direct counts matrix: verify provenance"))
  stop("Unsupported RDS object class: ", paste(class(obj), collapse = ", "))
}
get_metadata <- function(obj, cells) {
  if (inherits(obj, "Seurat"))
    m <- obj@meta.data else if (inherits(obj, "SummarizedExperiment"))
    m <- as.data.frame(SummarizedExperiment::colData(obj)) else m <- data.frame(row.names = cells)
  if (is.null(rownames(m)))
    stop("Metadata lacks cell rownames")
  n_matched <- sum(cells %in% rownames(m))
  if (n_matched != length(cells))
    stop("Metadata cannot be aligned to all expression cell IDs: ", n_matched, "/", length(cells),
      " exact matches")
  m[cells, , drop = FALSE]
}
select_col <- function(m, options, exclude = character()) {
  nm <- names(m)
  hits <- options[options %in% nm & !(options %in% exclude)]
  if (length(hits))
    hits[1] else NA_character_
}
lookup <- function(m, key, fallback = "Unknown") {
  if (is.na(key))
    return(rep(fallback, nrow(m)))
  x <- as.character(m[[key]])
  x[is.na(x) | !nzchar(trimws(x))] <- fallback
  x
}
row_sums <- function(m) if (inherits(m, "Matrix")) Matrix::rowSums(m) else base::rowSums(m)
col_sums <- function(m) if (inherits(m, "Matrix")) Matrix::colSums(m) else base::colSums(m)

summarize_gene_groups <- function(x, meta, group_cols, group_label, universe) {
  keys <- lapply(group_cols, function(n) meta[[n]])
  key <- do.call(paste, c(keys, sep = "|"))
  groups <- split(seq_len(ncol(x)), key)
  out <- lapply(groups, function(idx) {
    y <- x[, idx, drop = FALSE]
    n <- length(idx)
    sums <- as.numeric(row_sums(y))
    detected <- as.numeric(row_sums(y > 0))
    result <- data.frame(gene = rownames(y), n_cells = n, n_detected = detected, pct_detected = detected/n,
      mean_UMI_all = sums/n, mean_UMI_detected = ifelse(detected > 0, sums/pmax(1, detected), NA_real_),
      sum_UMI = sums, data_layer = universe, grouping = group_label, stringsAsFactors = FALSE)
    for (col in group_cols) result[[col]] <- meta[[col]][idx[1]]
    result
  })
  if (!length(out))
    return(data.frame())
  do.call(rbind, out)
}
status <- character()
record_status <- function(where, ok, details) {
  status <<- c(status, paste(where, if (ok) "PASS" else "LIMITED", details, sep = " | "))
}

atlas <- readRDS(input_files[["focus_atlas"]])
atlas_data <- get_matrix(atlas, "counts")
X <- atlas_data$matrix
if (is.null(colnames(X)) || is.null(rownames(X))) stop("Step 20 counts lack gene/cell names")
if (anyDuplicated(colnames(X)) || anyDuplicated(rownames(X))) stop("Step 20 contains duplicate cell or gene names; cannot unambiguously summarize")
meta_raw <- get_metadata(atlas, colnames(X))
if (!"myeloid_state" %in% names(meta_raw)) stop("'myeloid_state' column missing from Step 20; do not infer focal state from expression")
state <- lookup(meta_raw, "myeloid_state")
state[state == "Clec7a+ lipid-associated Microglia"] <- "Clec7a+ lipid-handling Microglia"
focus <- state == "Clec7a+ lipid-handling Microglia"
if (!any(focus)) stop("No Clec7a+ lipid-handling Microglia in Step 20 metadata")
source_col <- select_col(meta_raw, c("dataset", "source_dataset", "dataset_id", "source", "Dataset",
  "dataset_name", "study", "data_source"))
stage_col <- select_col(meta_raw, c("timepoint", "stage", "injury_stage", "dpi", "time_point", "timepoint_label",
  "time", "Timepoint"))
sample_col <- select_col(meta_raw, c("sample_id", "sample", "sampleID", "orig.ident", "library", "library_id",
  "GSM", "animal_id"))
met <- data.frame(cell_id = colnames(X), cell_group = ifelse(focus, "Clec7a_focus", "Other_myeloid"),
  source = lookup(meta_raw, source_col), stage = lookup(meta_raw, stage_col), sample_or_library = lookup(meta_raw,
    sample_col), stringsAsFactors = FALSE)
record_status("Step20", TRUE, paste("N cells", ncol(X), "focus", sum(focus), "assay", atlas_data$layer))
record_status("Source field", !is.na(source_col), ifelse(is.na(source_col), "not identified", source_col))
record_status("Stage field", !is.na(stage_col), ifelse(is.na(stage_col), "not identified", stage_col))
record_status("Sample field", !is.na(sample_col), ifelse(is.na(sample_col), "not identified", sample_col))
meta_audit <- data.frame(column_name = names(meta_raw), detected_as = ifelse(names(meta_raw) == "myeloid_state",
  "state", ifelse(names(meta_raw) %in% c(source_col, stage_col, sample_col), "possible source/stage/sample metadata",
    "other")), stringsAsFactors = FALSE)
out_csv(meta_audit, "05_136_atlas_metadata_column_audit_v01.csv")

available <- intersect(markers, rownames(X))
missing_markers <- setdiff(markers, rownames(X))
Xsmall <- X[available, , drop = FALSE]
rm(X, atlas, atlas_data, meta_raw)
gc(verbose = FALSE)
markers_df$present_in_Step20_RNA_counts <- markers_df$gene %in% available
out_csv(markers_df, "04_136_myelin_marker_and_annotation_manifest_v01.csv")
record_status("Myelin GO coverage", length(intersect(my_le_names, available)) > 0, paste(length(intersect(my_le_names,
  available)), "/", length(my_le_names), "leading-edge genes in Step20"))

a <- summarize_gene_groups(Xsmall, met, c("source", "cell_group"), "source x cell_group", "Step20 RNA counts")
b <- summarize_gene_groups(Xsmall, met, c("source", "stage", "cell_group"), "source x stage x cell_group",
  "Step20 RNA counts")
c <- summarize_gene_groups(Xsmall, met, c("source", "sample_or_library", "cell_group"), "source x sample_or_library x cell_group",
  "Step20 RNA counts")
a$is_myelin_leading_edge <- a$gene %in% my_le_names
b$is_myelin_leading_edge <- b$gene %in% my_le_names
c$is_myelin_leading_edge <- c$gene %in% my_le_names
out_csv(a, "06_136_gene_detection_by_source_and_state_v01.csv")
out_csv(b, "07_136_gene_detection_by_source_stage_state_v01.csv")
out_csv(c, "08_136_gene_detection_by_sample_or_library_v01.csv")

per_cell_hits <- function(g) {
  valid <- intersect(g, rownames(Xsmall))
  if (!length(valid))
    return(rep(NA_integer_, ncol(Xsmall)))
  as.integer(col_sums(Xsmall[valid, , drop = FALSE] > 0))
}
myelin_hits <- per_cell_hits(canonical_myelin)
phago_hits <- per_cell_hits(phago_context)
oligo_hits <- per_cell_hits(oligolineage)
marker_group <- split(seq_len(nrow(met)), paste(met$source, met$stage, met$cell_group, sep = "|"))
co_res <- do.call(rbind, lapply(marker_group, function(idx) {
  safe_mean <- function(x) if (all(is.na(x)))
    NA_real_ else mean(x, na.rm = TRUE)
  data.frame(source = met$source[idx[1]], stage = met$stage[idx[1]], cell_group = met$cell_group[idx[1]],
    n_cells = length(idx), frac_any_canonical_myelin = safe_mean(myelin_hits[idx] >= 1), frac_at_least_2_myelin_markers = safe_mean(myelin_hits[idx] >=
      2), frac_any_oligolineage_marker = safe_mean(oligo_hits[idx] >= 1), frac_at_least_2_phago_context_markers = safe_mean(phago_hits[idx] >=
      2), frac_myelin_and_phago_codetection = safe_mean(myelin_hits[idx] >= 1 & phago_hits[idx] >=
      2), cell_origin_determined = FALSE, stringsAsFactors = FALSE)
}))
out_csv(co_res, "09_136_myelin_phagocytosis_marker_codetection_context_v01.csv")

qc_pairs <- list()
pair_tables <- list()
for (src in c("17702045", "GSE162610", "GSE196928_Control")) {
  raw_key <- paste0("raw_", src)
  corr_key <- paste0("corrected_", src)
  if (!file.exists(input_files[[raw_key]]) || !file.exists(input_files[[corr_key]])) {
    qc_pairs[[src]] <- data.frame(source = src, status = "SKIP: missing RDS input", n_raw = NA_integer_,
      n_corrected = NA_integer_, n_exact_matched = 0L, n_matched_Step20 = 0L, raw_layer = "", corrected_layer = "",
      stringsAsFactors = FALSE)
    next
  }
  one <- tryCatch((function() {
    obj0 <- readRDS(input_files[[raw_key]])
    dat0 <- get_matrix(obj0, "counts")
    x0 <- dat0$matrix
    obj1 <- readRDS(input_files[[corr_key]])
    dat1 <- get_matrix(obj1, "corrected")
    x1 <- dat1$matrix
    if (anyDuplicated(colnames(x0)) || anyDuplicated(colnames(x1)))
      stop("duplicate cell IDs in raw/corrected input")
    exact <- intersect(colnames(x0), colnames(x1))
    audit <- data.frame(source = src, status = if (length(exact))
      "paired exact IDs" else "NO exact matched cells", n_raw = ncol(x0), n_corrected = ncol(x1), n_exact_matched = length(exact),
      n_matched_Step20 = sum(exact %in% met$cell_id), raw_layer = dat0$layer, corrected_layer = dat1$layer,
      stringsAsFactors = FALSE)
    if (!length(exact))
      return(list(audit = audit, data = NULL))

    same_genes <- intersect(intersect(rownames(x0), rownames(x1)), markers)
    if (!length(same_genes)) {
      audit$status <- "NO common exact gene symbols"
      return(list(audit = audit, data = NULL))
    }
    row0 <- x0[same_genes, exact, drop = FALSE]
    row1 <- x1[same_genes, exact, drop = FALSE]

    check0 <- as.numeric(row0[seq_len(min(nrow(row0), 4L)), seq_len(min(ncol(row0), 4L)), drop = FALSE])
    check1 <- as.numeric(row1[seq_len(min(nrow(row1), 4L)), seq_len(min(ncol(row1), 4L)), drop = FALSE])
    if (any(!is.finite(c(check0, check1))) || any(c(check0, check1) < 0))
      stop("invalid/negative values in count assays")

    idx_atlas <- match(exact, met$cell_id)
    known <- !is.na(idx_atlas)
    compatible <- rep(FALSE, length(exact))
    if (any(known)) {
      exp_src <- sub("_Control$", "", src)
      label <- met$source[idx_atlas[known]]
      compatible[known] <- grepl(exp_src, label, fixed = TRUE)
    }
    group <- rep("source_mixed_unlinked_to_focus", length(exact))
    good <- known & compatible
    group[good] <- met$cell_group[idx_atlas[good]]
    audit$n_matched_Step20 <- sum(good)
    audit$status <- if (sum(good) > 0)
      "matched pair; some Step20 focal labels" else "matched pair; Step20 source-specific focal labels unavailable"
    out <- lapply(split(seq_along(exact), group), function(idx) {
      a0 <- row0[, idx, drop = FALSE]
      a1 <- row1[, idx, drop = FALSE]
      sum_raw <- as.numeric(row_sums(a0))
      sum_corr <- as.numeric(row_sums(a1))
      n_raw <- as.numeric(row_sums(a0 > 0))
      n_corr <- as.numeric(row_sums(a1 > 0))
      data.frame(source = src, cell_group = group[idx[1]], gene = same_genes, n_exact_paired_cells = length(idx),
        n_detect_raw = n_raw, n_detect_corrected = n_corr, frac_detect_raw = n_raw/length(idx),
        frac_detect_corrected = n_corr/length(idx), detection_fraction_change_corrected_minus_raw = (n_corr -
          n_raw)/length(idx), sum_raw_UMI = sum_raw, sum_corrected_UMI = sum_corr, relative_UMI_removed = ifelse(sum_raw >
          0, (sum_raw - sum_corr)/sum_raw, NA_real_), is_myelin_leading_edge = same_genes %in%
          my_le_names, direct_ambient_origin_proven = FALSE, stringsAsFactors = FALSE)
    })
    list(audit = audit, data = do.call(rbind, out))
  })(), error = function(e) {
    list(audit = data.frame(source = src, status = paste("SKIP:", conditionMessage(e)), n_raw = NA_integer_,
      n_corrected = NA_integer_, n_exact_matched = 0L, n_matched_Step20 = 0L, raw_layer = "", corrected_layer = "",
      stringsAsFactors = FALSE), data = NULL)
  })
  qc_pairs[[src]] <- one$audit
  if (!is.null(one$data))
    pair_tables[[src]] <- one$data
  gc(verbose = FALSE)
}
qc_ambient <- do.call(rbind, qc_pairs)
out_csv(qc_ambient, "10_136_DecontX_paired_cell_ID_QC_v01.csv")
if (length(pair_tables)) ambient <- do.call(rbind, pair_tables) else ambient <- data.frame(source = character(),
  cell_group = character(), gene = character(), n_exact_paired_cells = integer(), n_detect_raw = integer(),
  n_detect_corrected = integer(), frac_detect_raw = numeric(), frac_detect_corrected = numeric(), detection_fraction_change_corrected_minus_raw = numeric(),
  sum_raw_UMI = numeric(), sum_corrected_UMI = numeric(), relative_UMI_removed = numeric(), is_myelin_leading_edge = logical(),
  direct_ambient_origin_proven = logical())
out_csv(ambient, "11_136_DecontX_raw_vs_corrected_exact_matched_v01.csv")
record_status("DecontX matching", any(qc_ambient$n_exact_matched > 0), paste("exact matched cells total",
  sum(qc_ambient$n_exact_matched), "linked to atlas focus", sum(qc_ambient$n_matched_Step20)))

assessment <- data.frame(hypothesis_note = c("GO shared annotation", "Possible phagocytosed material",
  "Ambient RNA", "Genuine cell-intrinsic expression", "De novo myelin synthesis or repair", "Zfp148 knockout-induced expression change"),
  existing_evidence = c("GO leading-edge overlaps, exact shared gene symbols", "Phagocytosis/lipid-handling context marker detection and co-detection only",
    "Existing DecontX raw/corrected counts in exact-matched cell subsets where available", "Detection in Clec7a+ atlas cells and across known sources/stages",
    "No protein localisation, myelin formation assay or functional repair assay", "Step 102 GO-GSEA ranks scTenifoldKnk network-distance Z"),
  interpretation_status = c("Can quantify annotation redundancy descriptively", "Cannot establish ingested myelin-derived transcripts",
    "Can assess contamination sensitivity; origin of remaining RNA is unproven", "Supports detected RNA, not proof of cell-autonomous transcription",
    "Not established", "Not established"), forbidden_inference = c("Distinct GO terms must represent independent pathways",
    "Marker co-detection proves myelin debris phagocytosis", "DecontX-removed counts prove ambient origin or retained counts prove native transcription",
    "Detected myelin mRNA proves microglial myelination", "Clec7a+ microglia synthesize/repair myelin",
    "Positive network Z implies expression upregulation or pathway activation"), stringsAsFactors = FALSE)
out_csv(assessment, "12_136_interpretation_and_evidence_limits_v01.csv")
qc <- data.frame(check = c("GO:0043209 CC row", "GO leading-edge genes", "GO full members", "GO leading subset full",
  "Step20 focus cells", "Step20 unique cell IDs", "Step20 myelin leading gene coverage", "Matched DecontX cell-pairs",
  "Pairs linked to atlas by source-qualified exact ID"), value = c(nrow(my_cc), length(my_le_names),
  length(my_member_names), length(setdiff(my_le_names, my_member_names)) == 0, sum(focus), nrow(met) ==
    length(unique(met$cell_id)), length(intersect(my_le_names, available)), sum(qc_ambient$n_exact_matched),
  sum(qc_ambient$n_matched_Step20)), interpretation = c("should be exactly 1", "Step102 core", "Step103 GO membership",
  "TRUE required", "Step20 state", "TRUE required", "some may be absent", "zero is allowed; check pairing QC",
  "zero is allowed; no focal DecontX conclusion"), stringsAsFactors = FALSE)
out_csv(qc, "13_136_FINAL_QC_v01.csv")
file_index <- list.files(table_dir, pattern = "\\.csv$", full.names = TRUE)
idx <- data.frame(file = basename(file_index), bytes = file.info(file_index)$size, stringsAsFactors = FALSE)
out_csv(idx, "00_136_OUTPUT_INDEX_v01.csv")
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(output_dir, "14_136_REPORT_v01.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "15_136_sessionInfo_v01.txt"), useBytes = TRUE)
