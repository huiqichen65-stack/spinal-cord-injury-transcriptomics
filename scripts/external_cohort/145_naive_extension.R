source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
options(future.globals.maxSize = 16 * 1024^3)

ROOT <- PROJECT_ROOT
DATA_DIR <- file.path(ROOT, "data", "GSE281205")
OUT <- file.path(ROOT, "results", "145_Reviewer1_RemainingMajor1_GSE281205_naive_extension")
FIG_DIR <- file.path(OUT, "figures")
RDS_DIR <- file.path(OUT, "RDS")
LOG_DIR <- file.path(OUT, "logs")

D123 <- file.path(ROOT, "results", "123_Reviewer1_Major4_GSE281205_major_lineage_clustering")
D128 <- file.path(ROOT, "results", "128_Reviewer1_Major4_GSE281205_reviewed_lineage_finalization")
D130 <- file.path(ROOT, "results", "130_Reviewer1_Major4_GSE281205_myeloid_reclustering")
D131 <- file.path(ROOT, "results", "131_Reviewer1_Major4_GSE281205_myeloid_subtype_annotation_and_microglia_groups")

IN_REF_MAJOR_RDS <- file.path(D123, "RDS", "123_GSE281205_Harmony_major_lineage_clustering_v03.rds")
IN_REF_LINEAGE <- file.path(D128, "03_128_per_cell_curated_lineage_v02.csv")
IN_REF_MYE_RDS <- file.path(D130, "RDS", "130_GSE281205_myeloid_LogNorm_PCA_Harmony_clusters_v02.rds")
IN_REF_131_MANIFEST <- file.path(D131, "04_131_per_cell_subtype_and_analysis_manifest_v01.csv")
IN_131_COMPOSITION <- file.path(D131, "05_131_GSM_primary_group_composition_v01.csv")

NAIVE_PLAN <- data.frame( gsm = paste0("GSM", 8614659:8614661), group = "Young-naive", stage = "Naive", replicate = 1:3, stringsAsFactors = FALSE )

OLD_GSM_14 <- paste0("GSM", 8614662:8614664)
OLD_GSM_28 <- paste0("GSM", 8614665:8614667)
ALL_EXPECTED_GSM <- c(NAIVE_PLAN$gsm, OLD_GSM_14, OLD_GSM_28)

MIN_FEATURES <- 300L
MAX_MT_PERCENT <- 10
DBR_PER_1K <- 0.008
PRIMARY_CONTAM_LT <- 20
MIN_MICROGLIA_CORE <- 2L
TRANSFER_SCORE_SENSITIVITY <- 0.50
N_DIMS_TRANSFER <- 30L
SEED_BASE_DECONTX <- 145100L
SEED_BASE_DBL <- 145200L
SEED_BASE_TRANSFER <- 145300L

MAJOR_MYELOID_LINEAGES <- c("Microglia", "Monocyte_macrophage", "APC_dendritic")
IDENTITY_MARKERS <- c("Fcrls", "Hexb", "Cx3cr1", "Ccr2", "Ly6c2", "Clec7a", "P2ry12", "Tmem119", "Spp1", "Lpl", "Gpnmb", "Cst7")

for (d in c(OUT, FIG_DIR, RDS_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

FINAL_SENTINEL <- file.path(OUT, "16_145_OUTPUT_INDEX_v01.csv")
if (file.exists(FINAL_SENTINEL)) {
  stop("A completed Step145 v01 output already exists: ", FINAL_SENTINEL,
    "\nMove/rename the existing Step145 output directory before rerunning; this script does not silently overwrite a completed run.", call. = FALSE)
}

LOG_FILE <- file.path(LOG_DIR, "00_145_run_log_v01.log")
cat("", file = LOG_FILE)
log_msg <- function(...) {
  z <- paste0(..., collapse = "")
  cat(z, "\n", sep = "")
  cat(z, "\n", file = LOG_FILE, append = TRUE, sep = "")
  invisible(z)
}
stopif <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}
read_csv <- function(f) {
  stopif(file.exists(f), paste("Missing required input:", f))
  utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}
write_csv <- function(x, filename) {
  stopif(is.data.frame(x), paste("Expected data.frame:", filename))
  p <- file.path(OUT, filename)
  utils::write.csv(x, p, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  log_msg("SAVED ", filename, " [", nrow(x), " rows]")
  invisible(p)
}
write_csv_gz <- function(x, filename) {
  stopif(is.data.frame(x), paste("Expected data.frame:", filename))
  p <- file.path(OUT, filename)
  con <- gzfile(p, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "")
  close(con)
  on.exit(NULL, add = FALSE)
  log_msg("SAVED ", filename, " [", nrow(x), " rows]")
  invisible(p)
}
need_cols <- function(d, cols, label) {
  miss <- setdiff(cols, names(d))
  stopif(!length(miss), paste(label, "missing columns:", paste(miss, collapse = ", ")))
}
median_safe <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x))) return(NA_real_)
  stats::median(x[is.finite(x)])
}
mean_safe <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x))) return(NA_real_)
  mean(x[is.finite(x)])
}

required_pkgs <- c("Matrix", "Seurat", "SeuratObject", "SingleCellExperiment", "SummarizedExperiment", "scDblFinder", "BiocParallel")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

if (requireNamespace("decontX", quietly = TRUE)) {
  DECONTX_FUNCTION <- decontX::decontX
  DECONTX_PACKAGE <- "decontX"
} else if (requireNamespace("celda", quietly = TRUE)) {
  DECONTX_FUNCTION <- celda::decontX
  DECONTX_PACKAGE <- "celda"
} else {
  stop("Missing DecontX implementation. Install 'decontX' or 'celda'.", call. = FALSE)
}

resolve_naive_triplet <- function(gsm, replicate) {
  prefix <- paste0(gsm, "_Young-naive_", replicate, "_")
  ans <- c( barcodes = file.path(DATA_DIR, paste0(prefix, "barcodes.tsv.gz")), features = file.path(DATA_DIR, paste0(prefix, "features.tsv.gz")),
    matrix = file.path(DATA_DIR, paste0(prefix, "matrix.mtx.gz")) )
  stopif(all(file.exists(ans)), paste("Missing naive 10x file(s) for", gsm, ":", paste(ans[!file.exists(ans)], collapse = "; ")))
  ans
}

read_tsv_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  utils::read.delim(con, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
}

read_10x_triplet <- function(gsm, replicate) {
  f <- resolve_naive_triplet(gsm, replicate)
  feat <- read_tsv_gz(f[["features"]])
  bc <- read_tsv_gz(f[["barcodes"]])
  stopif(ncol(feat) >= 2L, paste("features.tsv.gz has <2 columns:", gsm))
  stopif(ncol(bc) >= 1L, paste("barcodes.tsv.gz empty:", gsm))
  con <- gzfile(f[["matrix"]], open = "rb")
  on.exit(close(con), add = TRUE)
  mat <- Matrix::readMM(con)
  close(con)
  on.exit(NULL, add = FALSE)
  mat <- as(mat, "dgCMatrix")
  stopif(nrow(mat) == nrow(feat) && ncol(mat) == nrow(bc), paste("10x dimensions disagree for", gsm))
  gene_id <- as.character(feat[[1]])
  gene_symbol <- as.character(feat[[2]])
  stopif(!anyNA(gene_id) && !anyDuplicated(gene_id), paste("Gene IDs missing/duplicated in", gsm))
  raw_barcode <- as.character(bc[[1]])
  stopif(!anyNA(raw_barcode) && !anyDuplicated(raw_barcode), paste("Barcodes missing/duplicated in", gsm))
  rownames(mat) <- gene_id
  cell_id <- paste0(gsm, "_", raw_barcode)
  colnames(mat) <- cell_id
  list( counts = mat, feature_metadata = data.frame(gene_id = gene_id, gene_symbol = gene_symbol, stringsAsFactors = FALSE), raw_barcode = raw_barcode, files = f )
}

cell_qc <- function(obj, gsm, replicate) {
  x <- obj$counts
  fm <- obj$feature_metadata
  symbols <- as.character(fm$gene_symbol)
  mt <- grepl("^mt-", symbols, ignore.case = TRUE)
  nCount <- as.numeric(Matrix::colSums(x))
  nFeature <- as.numeric(Matrix::colSums(x > 0))
  mtCount <- if (any(mt)) as.numeric(Matrix::colSums(x[mt, , drop = FALSE])) else rep(0, ncol(x))
  pctMt <- ifelse(nCount > 0, 100 * mtCount / nCount, 0)
  pass <- nFeature >= MIN_FEATURES & pctMt < MAX_MT_PERCENT
  md <- data.frame( cell_id = colnames(x), gsm = gsm, stage = "Naive", replicate = replicate, nCount_RNA = nCount, nFeature_RNA = nFeature, percent.mt = pctMt,
    QC_pass = pass, stringsAsFactors = FALSE )
  stopif(sum(pass) > 100L, paste("Too few QC-passed cells in", gsm))
  list(counts = x[, pass, drop = FALSE], metadata = md[pass, , drop = FALSE], feature_metadata = fm, summary = data.frame( gsm = gsm, stage = "Naive", replicate = replicate,
    n_raw_barcodes = ncol(x), n_QC_retained = sum(pass), retention_percent = 100 * mean(pass), median_nFeature_after = median_safe(nFeature[pass]),
    median_nCount_after = median_safe(nCount[pass]), median_percent_mt_after = median_safe(pctMt[pass]), stringsAsFactors = FALSE ))
}

run_decontx <- function(qc, gsm, sample_index) {
  raw <- qc$counts
  expressed <- which(Matrix::rowSums(raw) > 0)
  stopif(length(expressed) >= 1000L, paste("Too few expressed genes for DecontX:", gsm))
  model_input <- raw[expressed, , drop = FALSE]
  seed <- as.integer(SEED_BASE_DECONTX + sample_index)
  model_log <- file.path(LOG_DIR, paste0("145_", gsm, "_DecontX_model_v01.log"))
  log_msg("DecontX ", gsm, " | cells=", ncol(raw), " | expressed genes=", length(expressed), " | seed=", seed)
  set.seed(seed)
  res <- DECONTX_FUNCTION(x = model_input, z = NULL, seed = seed, logfile = model_log, verbose = TRUE)
  stopif(is.list(res) && all(c("decontXcounts", "contamination") %in% names(res)), paste("Unexpected DecontX return structure:", gsm))
  fit <- as(res$decontXcounts, "dgCMatrix")
  if (is.null(rownames(fit))) rownames(fit) <- rownames(model_input)
  if (is.null(colnames(fit))) colnames(fit) <- colnames(model_input)
  stopif(identical(dimnames(fit), dimnames(model_input)), paste("DecontX changed dimnames:", gsm))
  fit@x[fit@x < 0] <- 0
  fit <- Matrix::drop0(fit)
  corrected <- as(Matrix::Matrix(0, nrow = nrow(raw), ncol = ncol(raw), sparse = TRUE), "dgCMatrix")
  dimnames(corrected) <- dimnames(raw)
  corrected[expressed, ] <- fit
  corrected <- Matrix::drop0(corrected)
  contamination <- pmax(0, pmin(1, as.numeric(res$contamination)))
  stopif(length(contamination) == ncol(raw) && all(is.finite(contamination)), paste("Invalid DecontX contamination:", gsm))
  md <- qc$metadata
  stopif(identical(md$cell_id, colnames(raw)), paste("QC metadata order mismatch:", gsm))
  md$decontX_contamination_fraction <- contamination
  md$decontX_contamination_percent <- 100 * contamination
  list(corrected = corrected, raw = raw, metadata = md, feature_metadata = qc$feature_metadata)
}

run_doublets <- function(dc, gsm, sample_index) {
  raw <- dc$raw
  usable <- which(Matrix::rowSums(raw) > 0)
  stopif(length(usable) >= 1000L, paste("Too few expressed features for scDblFinder:", gsm))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = raw[usable, , drop = FALSE]))
  seed <- as.integer(SEED_BASE_DBL + sample_index)
  set.seed(seed)
  modeled <- scDblFinder::scDblFinder( sce, clusters = NULL, samples = NULL, dbr = NULL, dbr.per1k = DBR_PER_1K, BPPARAM = BiocParallel::SerialParam(), verbose = TRUE )
  cd <- SummarizedExperiment::colData(modeled)
  stopif(all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(cd)), paste("scDblFinder columns missing:", gsm))
  cls <- as.character(cd$scDblFinder.class)
  score <- as.numeric(cd$scDblFinder.score)
  keep <- cls == "singlet"
  stopif(sum(keep) > 100L, paste("Too few predicted singlets:", gsm))
  md <- dc$metadata
  md$scDblFinder_score <- score
  md$scDblFinder_class <- cls
  md$predicted_singlet <- keep
  list( corrected = dc$corrected[, keep, drop = FALSE], raw = raw[, keep, drop = FALSE], metadata = md[keep, , drop = FALSE], feature_metadata = dc$feature_metadata, summary = data.frame(
    gsm = gsm, stage = "Naive", n_QC_input = ncol(raw), n_predicted_singlets = sum(keep), n_predicted_doublets = sum(!keep), doublet_percent = 100 * mean(!keep),
    median_contamination_percent = median_safe(md$decontX_contamination_percent), stringsAsFactors = FALSE ) )
}

load_frozen_references <- function() {
  required <- c(IN_REF_MAJOR_RDS, IN_REF_LINEAGE, IN_REF_MYE_RDS, IN_REF_131_MANIFEST, IN_131_COMPOSITION)
  stopif(all(file.exists(required)), paste("Missing frozen reference input(s):", paste(required[!file.exists(required)], collapse = "; ")))

  log_msg("Loading frozen Step123 reference: ", IN_REF_MAJOR_RDS)
  ref_major <- readRDS(IN_REF_MAJOR_RDS)
  stopif(inherits(ref_major, "Seurat"), "Step123 reference is not a Seurat object")
  stopif("DecontX" %in% SeuratObject::Assays(ref_major), "Step123 reference lacks DecontX assay")
  stopif("pca" %in% SeuratObject::Reductions(ref_major), "Step123 reference lacks PCA")

  lineage <- read_csv(IN_REF_LINEAGE)
  need_cols(lineage, c("cell_id", "curated_major_lineage_128"), "Step128 lineage")
  mi <- match(colnames(ref_major), lineage$cell_id)
  stopif(!anyNA(mi), "Step123 reference cells are not fully represented in Step128 lineage table")
  ref_major$curated_major_lineage_128 <- as.character(lineage$curated_major_lineage_128[mi])
  stopif(!anyNA(ref_major$curated_major_lineage_128), "Missing Step128 lineage labels after join")
  SeuratObject::DefaultAssay(ref_major) <- "DecontX"

  log_msg("Loading frozen Step130 myeloid reference: ", IN_REF_MYE_RDS)
  ref_mye <- readRDS(IN_REF_MYE_RDS)
  stopif(inherits(ref_mye, "Seurat"), "Step130 reference is not a Seurat object")
  stopif("DecontX" %in% SeuratObject::Assays(ref_mye), "Step130 reference lacks DecontX assay")
  stopif("pca" %in% SeuratObject::Reductions(ref_mye), "Step130 reference lacks PCA")
  m131 <- read_csv(IN_REF_131_MANIFEST)
  need_cols(m131, c("cell_id", "myeloid_harmony_r010", "primary_group_131"), "Step131 manifest")
  mi2 <- match(colnames(ref_mye), m131$cell_id)
  stopif(!anyNA(mi2), "Step130 reference cells are not fully represented in Step131 manifest")
  ref_mye$frozen_myeloid_r010_131 <- as.character(m131$myeloid_harmony_r010[mi2])
  stopif(all(ref_mye$frozen_myeloid_r010_131 %in% as.character(0:4)), "Unexpected frozen myeloid r0.10 labels")
  SeuratObject::DefaultAssay(ref_mye) <- "DecontX"

  list(ref_major = ref_major, ref_mye = ref_mye)
}

pca_loading_features <- function(seu, n_min = 1000L) {
  ld <- SeuratObject::Loadings(seu[["pca"]])
  f <- rownames(ld)
  stopif(length(f) >= n_min, "Reference PCA contains too few loading features")
  f
}

transfer_labels_one_sample <- function(sing, gsm, refs, sample_index) {
  corr <- sing$corrected
  md <- sing$metadata
  stopif(identical(colnames(corr), md$cell_id), paste("Singlet cell order mismatch:", gsm))

  rownames(md) <- md$cell_id
  q <- Seurat::CreateSeuratObject(counts = corr, assay = "DecontX", project = gsm, min.cells = 0, min.features = 0, meta.data = md)
  SeuratObject::DefaultAssay(q) <- "DecontX"
  q <- Seurat::NormalizeData(q, assay = "DecontX", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

  ref_major <- refs$ref_major
  major_features <- intersect(pca_loading_features(ref_major), rownames(q))
  stopif(length(major_features) >= 1000L, paste("Insufficient shared Step123 PCA features for", gsm, ":", length(major_features)))
  set.seed(SEED_BASE_TRANSFER + sample_index)
  log_msg("Major-lineage transfer ", gsm, " | shared PCA features=", length(major_features))
  anchors_major <- Seurat::FindTransferAnchors( reference = ref_major, query = q, normalization.method = "LogNormalize", reference.reduction = "pca", reduction = "pcaproject",
    features = major_features, dims = seq_len(N_DIMS_TRANSFER), verbose = FALSE )
  pred_major <- Seurat::TransferData( anchorset = anchors_major, refdata = as.character(ref_major$curated_major_lineage_128), dims = seq_len(N_DIMS_TRANSFER), verbose = FALSE )
  pred_major <- as.data.frame(pred_major, stringsAsFactors = FALSE)
  stopif(nrow(pred_major) == ncol(q), paste("Major transfer row count mismatch:", gsm))
  q$predicted_major_lineage_128 <- as.character(pred_major$predicted.id)
  q$major_lineage_prediction_score <- as.numeric(pred_major$prediction.score.max)

  qmd <- q[[]]
  qmd$cell_id <- rownames(qmd)
  is_candidate <- qmd$predicted_major_lineage_128 %in% MAJOR_MYELOID_LINEAGES
  stopif(sum(is_candidate) >= 50L, paste("Too few predicted candidate myeloid cells:", gsm))
  candidate_ids <- qmd$cell_id[is_candidate]
  q_mye <- subset(q, cells = candidate_ids)

  ref_mye <- refs$ref_mye
  mye_features <- intersect(pca_loading_features(ref_mye), rownames(q_mye))
  stopif(length(mye_features) >= 1000L, paste("Insufficient shared Step130 PCA features for", gsm, ":", length(mye_features)))
  set.seed(SEED_BASE_TRANSFER + 100L + sample_index)
  log_msg("Myeloid-state transfer ", gsm, " | candidate cells=", ncol(q_mye), " | shared PCA features=", length(mye_features))
  anchors_mye <- Seurat::FindTransferAnchors( reference = ref_mye, query = q_mye, normalization.method = "LogNormalize", reference.reduction = "pca", reduction = "pcaproject",
    features = mye_features, dims = seq_len(N_DIMS_TRANSFER), verbose = FALSE )
  pred_state <- Seurat::TransferData( anchorset = anchors_mye, refdata = as.character(ref_mye$frozen_myeloid_r010_131), dims = seq_len(N_DIMS_TRANSFER), verbose = FALSE )
  pred_state <- as.data.frame(pred_state, stringsAsFactors = FALSE)
  stopif(nrow(pred_state) == ncol(q_mye), paste("State transfer row count mismatch:", gsm))

  out <- q_mye[[]]
  out$cell_id <- rownames(out)
  out$transferred_myeloid_r010_131 <- as.character(pred_state$predicted.id)
  out$myeloid_state_prediction_score <- as.numeric(pred_state$prediction.score.max)
  stopif(all(out$transferred_myeloid_r010_131 %in% as.character(0:4)), paste("Unexpected transferred state label:", gsm))

  fm <- sing$feature_metadata
  marker_map <- do.call(rbind, lapply(IDENTITY_MARKERS, function(sym) {
    ids <- as.character(fm$gene_id[as.character(fm$gene_symbol) == sym])
    data.frame(gene_symbol = sym, gene_id = if (length(ids) == 1L) ids else NA_character_, n_mapping_hits = length(ids), stringsAsFactors = FALSE)
  }))
  stopif(all(marker_map$n_mapping_hits == 1L), paste("Naive feature mapping missing/ambiguous for", gsm, ":",
    paste(marker_map$gene_symbol[marker_map$n_mapping_hits != 1L], collapse = ", ")))
  ci <- match(out$cell_id, colnames(corr))
  stopif(!anyNA(ci), paste("Candidate cells cannot be aligned to DecontX matrix:", gsm))
  for (ii in seq_len(nrow(marker_map))) {
    sym <- marker_map$gene_symbol[ii]
    gid <- marker_map$gene_id[ii]
    v <- as.numeric(corr[gid, ci, drop = TRUE])
    out[[paste0("marker_", sym, "_ge1")]] <- is.finite(v) & v >= 1
  }
  out$microglia_core_markers_ge1 <- as.integer(out$marker_Fcrls_ge1) +
    as.integer(out$marker_Hexb_ge1) + as.integer(out$marker_Cx3cr1_ge1)
  out$microglia_core_supported <- out$microglia_core_markers_ge1 >= MIN_MICROGLIA_CORE
  out$microglia_identity_conflict <- out$marker_Ccr2_ge1 & out$marker_Ly6c2_ge1
  out$qc_lt20 <- as.numeric(out$decontX_contamination_percent) < PRIMARY_CONTAM_LT

  out$primary_target_145 <-
    out$predicted_major_lineage_128 == "Microglia" &
    out$transferred_myeloid_r010_131 == "1" &
    out$microglia_core_supported &
    !out$microglia_identity_conflict &
    out$qc_lt20

  out$high_confidence_transfer_sensitivity <-
    out$major_lineage_prediction_score >= TRANSFER_SCORE_SENSITIVITY &
    out$myeloid_state_prediction_score >= TRANSFER_SCORE_SENSITIVITY
  out$primary_target_highconf_transfer_sensitivity <-
    out$primary_target_145 & out$high_confidence_transfer_sensitivity

  state_name <- c( `0` = "P2ry12_high_microglia_like", `1` = "Clec7a_high_lipid_associated_microglia_like", `2` = "Monocyte_macrophage_like", `3` = "APC_dendritic_like",
    `4` = "Proliferating_microglia_like" )
  out$transferred_working_state_145 <- unname(state_name[out$transferred_myeloid_r010_131])

  lineage_tab <- as.data.frame(table( predicted_major_lineage_128 = factor(qmd$predicted_major_lineage_128) ), stringsAsFactors = FALSE)
  names(lineage_tab)[2] <- "n_cells"
  lineage_tab$gsm <- gsm
  lineage_tab <- lineage_tab[, c("gsm", "predicted_major_lineage_128", "n_cells")]
  lineage_tab$mean_prediction_score <- NA_real_
  for (k in seq_len(nrow(lineage_tab))) {
    lab <- as.character(lineage_tab$predicted_major_lineage_128[k])
    lineage_tab$mean_prediction_score[k] <- mean_safe(qmd$major_lineage_prediction_score[ qmd$predicted_major_lineage_128 == lab])
  }

  state_tab <- as.data.frame(table( transferred_myeloid_r010_131 = factor(out$transferred_myeloid_r010_131, levels = as.character(0:4)) ), stringsAsFactors = FALSE)
  names(state_tab)[2] <- "n_cells"
  state_tab$gsm <- gsm
  state_tab <- state_tab[, c("gsm", "transferred_myeloid_r010_131", "n_cells")]
  state_tab$working_state <- unname(state_name[as.character(state_tab$transferred_myeloid_r010_131)])
  state_tab$mean_prediction_score <- NA_real_
  for (k in seq_len(nrow(state_tab))) {
    lab <- as.character(state_tab$transferred_myeloid_r010_131[k])
    state_tab$mean_prediction_score[k] <- mean_safe(out$myeloid_state_prediction_score[ out$transferred_myeloid_r010_131 == lab])
  }

  rm(q, q_mye, anchors_major, anchors_mye, pred_major, pred_state)
  gc(verbose = FALSE)

  list(assignments = out, lineage_summary = lineage_tab, state_summary = state_tab, marker_map = marker_map)
}

exact_two_sided_permutation <- function(group_a, group_b) {
  a <- as.numeric(group_a); b <- as.numeric(group_b)
  stopif(length(a) >= 2L && length(b) >= 2L && all(is.finite(c(a, b))), "Invalid exact permutation input")
  z <- c(a, b)
  n_a <- length(a)
  obs <- mean(b) - mean(a)
  cmb <- utils::combn(seq_along(z), n_a)
  perm <- apply(cmb, 2L, function(idx_a) {
    idx_b <- setdiff(seq_along(z), idx_a)
    mean(z[idx_b]) - mean(z[idx_a])
  })
  p <- mean(abs(perm) >= abs(obs) - 1e-12)
  c(observed_difference = obs, exact_two_sided_p = p, n_permutations = length(perm))
}

welch_ci_b_minus_a <- function(a, b, conf = 0.95) {
  tt <- stats::t.test(b, a, alternative = "two.sided", var.equal = FALSE, conf.level = conf)
  c(lower = unname(tt$conf.int[1]), upper = unname(tt$conf.int[2]))
}

pairwise_effect_row <- function(specimen_tab, a, b) {
  xa <- specimen_tab$proportion[specimen_tab$stage == a]
  xb <- specimen_tab$proportion[specimen_tab$stage == b]
  stopif(length(xa) == 3L && length(xb) == 3L, paste("Expected n=3 per group for", a, "vs", b))
  ex <- exact_two_sided_permutation(xa, xb)
  ci <- welch_ci_b_minus_a(xa, xb)
  data.frame( group_A = a, group_B = b, n_A = length(xa), n_B = length(xb), mean_A_percent = 100 * mean(xa), mean_B_percent = 100 * mean(xb),
    absolute_mean_difference_B_minus_A_percentage_points = 100 * (mean(xb) - mean(xa)), welch_95CI_lower_percentage_points = 100 * ci[["lower"]],
    welch_95CI_upper_percentage_points = 100 * ci[["upper"]], exact_two_sided_permutation_P = unname(ex[["exact_two_sided_p"]]),
    exact_number_of_label_permutations = as.integer(ex[["n_permutations"]]),
    interpretation = "Specimen-level effect estimate; n=3 per group; CI is imprecise and exact P has coarse resolution", stringsAsFactors = FALSE )
}

exact_three_group_omnibus <- function(x_naive, x_14, x_28) {
  z <- c(x_naive, x_14, x_28)
  stopif(length(z) == 9L && all(is.finite(z)), "Omnibus exact test expects 9 finite values")
  ss_between <- function(groups) {
    gm <- mean(z)
    lev <- unique(groups)
    sum(vapply(lev, function(g) {
      zz <- z[groups == g]
      length(zz) * (mean(zz) - gm)^2
    }, numeric(1)))
  }
  obs_groups <- rep(c("Naive", "14dpi", "28dpi"), each = 3L)
  obs <- ss_between(obs_groups)
  idx <- seq_along(z)
  g1 <- utils::combn(idx, 3L)
  vals <- numeric(0)
  for (j in seq_len(ncol(g1))) {
    rem <- setdiff(idx, g1[, j])
    g2 <- utils::combn(rem, 3L)
    for (k in seq_len(ncol(g2))) {
      groups <- rep("28dpi", 9L)
      groups[g1[, j]] <- "Naive"
      groups[g2[, k]] <- "14dpi"
      vals <- c(vals, ss_between(groups))
    }
  }

  p <- mean(vals >= obs - 1e-12)
  data.frame( statistic = "between_group_sum_of_squares", observed_statistic = obs, exact_two_sided_or_omnibus_P = p, n_unique_group_allocations = length(vals),
    note = "Exact specimen-label permutation across Naive/14dpi/28dpi; descriptive omnibus test, not a trajectory model", stringsAsFactors = FALSE )
}

make_figures <- function(specimen_tab, effects) {
  stage_levels <- c("Naive", "14dpi", "28dpi")
  specimen_tab$stage <- factor(specimen_tab$stage, levels = stage_levels)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p1 <- ggplot2::ggplot(specimen_tab, ggplot2::aes(x = stage, y = 100 * proportion)) +
      ggplot2::geom_jitter(width = 0.06, height = 0, size = 3) +
      ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.45, fatten = 0, linewidth = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = gsm), position = ggplot2::position_jitter(width = 0.06, height = 0), vjust = -0.9, size = 2.7, check_overlap = TRUE) +
      ggplot2::labs( title = "GSE281205 specimen-level Clec7a-high candidate-state proportions",
        subtitle = "Young-naive added in Step145; 14/28-dpi values imported unchanged from frozen Step131", x = NULL, y = "Clec7a-high candidate state (% of candidate myeloid cells)"
        ) + ggplot2::theme_bw(base_size = 11)
    ggplot2::ggsave(file.path(FIG_DIR, "Fig_145_GSE281205_naive_14_28_specimen_proportions_v01.png"), p1, width = 9, height = 6, dpi = 180)

    ef <- effects
    ef$comparison <- paste0(ef$group_B, " - ", ef$group_A)
    ef$comparison <- factor(ef$comparison, levels = rev(ef$comparison))
    p2 <- ggplot2::ggplot(ef, ggplot2::aes( x = absolute_mean_difference_B_minus_A_percentage_points, y = comparison)) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_errorbarh(ggplot2::aes( xmin = welch_95CI_lower_percentage_points, xmax = welch_95CI_upper_percentage_points), height = 0.15) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs( title = "Specimen-level absolute effect sizes", subtitle = "Points: mean difference; bars: Welch 95% CI (n=3 per group)",
        x = "Difference in candidate-state proportion (percentage points)", y = NULL ) + ggplot2::theme_bw(base_size = 11)
    ggplot2::ggsave(file.path(FIG_DIR, "Fig_145_GSE281205_pairwise_effect_sizes_v01.png"), p2, width = 8.5, height = 5, dpi = 180)
  } else {
    png(file.path(FIG_DIR, "Fig_145_GSE281205_naive_14_28_specimen_proportions_v01.png"), width = 1500, height = 1000, res = 180)
    on.exit(dev.off(), add = TRUE)
    xx <- as.numeric(specimen_tab$stage)
    plot(xx, 100 * specimen_tab$proportion, xaxt = "n", pch = 19, xlab = "", ylab = "Clec7a-high candidate state (%)", main = "GSE281205 specimen-level proportions", xlim = c(0.5, 3.5))
    axis(1, at = 1:3, labels = stage_levels)
    for (i in 1:3) segments(i - .25, 100 * mean(specimen_tab$proportion[xx == i]), i + .25, 100 * mean(specimen_tab$proportion[xx == i]), lwd = 2)
    dev.off(); on.exit(NULL, add = FALSE)
  }
}

run <- function() {
  t0 <- Sys.time()
  log_msg("============================================================")
  log_msg("STEP145 v01 | Reviewer1 Remaining Major Comment 1")
  log_msg("GSE281205 Young-naive extension with frozen 14/28 reference")
  log_msg("Start: ", format(t0, "%Y-%m-%d %H:%M:%S %Z"))
  log_msg("Data: ", DATA_DIR)
  log_msg("Output: ", OUT)
  log_msg("NO ZIP DEPENDENCY")
  log_msg("Existing 14/28 Step131 results are imported unchanged")
  log_msg("============================================================")

  ref_paths <- c( Step123_major_reference_RDS = IN_REF_MAJOR_RDS, Step128_curated_lineage_CSV = IN_REF_LINEAGE, Step130_myeloid_reference_RDS = IN_REF_MYE_RDS,
    Step131_analysis_manifest_CSV = IN_REF_131_MANIFEST, Step131_specimen_composition_CSV = IN_131_COMPOSITION )
  raw_rows <- list()
  for (i in seq_len(nrow(NAIVE_PLAN))) {
    fs <- resolve_naive_triplet(NAIVE_PLAN$gsm[i], NAIVE_PLAN$replicate[i])
    for (nm in names(fs)) {
      raw_rows[[length(raw_rows) + 1L]] <- data.frame( role = paste0("naive_", nm), gsm = NAIVE_PLAN$gsm[i], path = fs[[nm]],
        exists = file.exists(fs[[nm]]), bytes = as.numeric(file.info(fs[[nm]])$size), stringsAsFactors = FALSE)
    }
  }
  raw_audit <- do.call(rbind, raw_rows)
  write_csv(raw_audit, "01_145_naive_raw_file_audit_v01.csv")
  prov <- data.frame(role = names(ref_paths), path = unname(ref_paths), exists = file.exists(unname(ref_paths)), bytes = as.numeric(file.info(unname(ref_paths))$size),
    stringsAsFactors = FALSE)
  write_csv(prov, "02_145_frozen_reference_provenance_v01.csv")
  stopif(all(raw_audit$exists) && all(prov$exists), "One or more required inputs are missing")

  refs <- load_frozen_references()

  all_assign <- list()
  qc_summary <- list()
  dbl_summary <- list()
  lin_summary <- list()
  state_summary <- list()
  marker_maps <- list()

  for (i in seq_len(nrow(NAIVE_PLAN))) {
    gsm <- NAIVE_PLAN$gsm[i]
    repn <- NAIVE_PLAN$replicate[i]
    log_msg("------------------------------------------------------------")
    log_msg("Processing naive specimen ", gsm, " replicate ", repn)
    raw <- read_10x_triplet(gsm, repn)
    q <- cell_qc(raw, gsm, repn)
    qc_summary[[gsm]] <- q$summary
    log_msg("QC retained ", q$summary$n_QC_retained, "/", q$summary$n_raw_barcodes, " | ", round(q$summary$retention_percent, 2), "%")
    rm(raw); gc(verbose = FALSE)

    dc <- run_decontx(q, gsm, i)
    rm(q); gc(verbose = FALSE)
    sing <- run_doublets(dc, gsm, i)
    dbl_summary[[gsm]] <- sing$summary
    log_msg("scDblFinder singlets ", sing$summary$n_predicted_singlets, " | doublets ", sing$summary$n_predicted_doublets)
    rm(dc); gc(verbose = FALSE)

    tr <- transfer_labels_one_sample(sing, gsm, refs, i)
    all_assign[[gsm]] <- tr$assignments
    lin_summary[[gsm]] <- tr$lineage_summary
    state_summary[[gsm]] <- tr$state_summary
    marker_maps[[gsm]] <- data.frame(gsm = gsm, tr$marker_map, stringsAsFactors = FALSE, check.names = FALSE)
    log_msg("Candidate myeloid=", nrow(tr$assignments), " | primary target=", sum(tr$assignments$primary_target_145),
      " | primary proportion=", round(100 * mean(tr$assignments$primary_target_145), 3), "%")
    rm(sing, tr); gc(verbose = FALSE)
  }

  qc_tab <- do.call(rbind, qc_summary)
  dbl_tab <- do.call(rbind, dbl_summary)
  lin_tab <- do.call(rbind, lin_summary)
  state_tab <- do.call(rbind, state_summary)
  marker_tab <- do.call(rbind, marker_maps)
  assign <- do.call(rbind, all_assign)
  rownames(assign) <- NULL
  stopif(!anyDuplicated(assign$cell_id), "Duplicate naive cell IDs after combination")

  write_csv(qc_tab, "03_145_naive_cell_QC_summary_v01.csv")
  write_csv(dbl_tab, "04_145_naive_DecontX_scDblFinder_summary_v01.csv")
  write_csv(lin_tab, "05_145_naive_major_lineage_transfer_summary_v01.csv")
  write_csv(state_tab, "06_145_naive_myeloid_state_transfer_summary_v01.csv")
  write_csv(marker_tab, "07_145_naive_marker_mapping_v01.csv")
  write_csv_gz(assign, "08_145_naive_per_cell_frozen_rule_assignments_v01.csv.gz")
  saveRDS(assign, file.path(RDS_DIR, "145_GSE281205_naive_assignment_manifest_v01.rds"), compress = "gzip")
  log_msg("SAVED compact RDS metadata manifest")

  naive_spec <- do.call(rbind, lapply(NAIVE_PLAN$gsm, function(g) {
    d <- assign[assign$gsm == g, , drop = FALSE]
    stopif(nrow(d) > 0L, paste("No candidate myeloid cells for", g))
    highconf_den <- d$major_lineage_prediction_score >= TRANSFER_SCORE_SENSITIVITY &
      d$myeloid_state_prediction_score >= TRANSFER_SCORE_SENSITIVITY
    data.frame( gsm = g, stage = "Naive", n_candidate_myeloid = nrow(d), n_primary_target = sum(d$primary_target_145), proportion = mean(d$primary_target_145),
      percent = 100 * mean(d$primary_target_145), n_candidate_highconf_transfer = sum(highconf_den), n_target_highconf_transfer = sum(d$primary_target_highconf_transfer_sensitivity),
      proportion_highconf_transfer_sensitivity = if (sum(highconf_den) > 0) sum(d$primary_target_highconf_transfer_sensitivity) / sum(highconf_den) else NA_real_,
      mean_major_lineage_prediction_score = mean_safe(d$major_lineage_prediction_score), mean_myeloid_state_prediction_score = mean_safe(d$myeloid_state_prediction_score),
      source = "Step145 naive frozen-reference transfer", stringsAsFactors = FALSE )
  }))

  old <- read_csv(IN_131_COMPOSITION)
  need_cols(old, c("gsm", "group_131", "n_cells", "dpi", "total_candidate_myeloid_in_GSM", "fraction_of_candidate_myeloid_in_GSM"), "Step131 composition")
  old_target <- old[old$group_131 == "TARGET_CLEC7A_HIGH_MICROGLIA", , drop = FALSE]
  stopif(nrow(old_target) == 6L && setequal(old_target$gsm, c(OLD_GSM_14, OLD_GSM_28)), "Frozen Step131 target composition does not contain the expected six GSMs")
  old_spec <- data.frame( gsm = as.character(old_target$gsm), stage = ifelse(as.integer(old_target$dpi) == 14L, "14dpi", "28dpi"),
    n_candidate_myeloid = as.integer(old_target$total_candidate_myeloid_in_GSM), n_primary_target = as.integer(old_target$n_cells),
    proportion = as.numeric(old_target$fraction_of_candidate_myeloid_in_GSM), percent = 100 * as.numeric(old_target$fraction_of_candidate_myeloid_in_GSM),
    n_candidate_highconf_transfer = NA_integer_, n_target_highconf_transfer = NA_integer_, proportion_highconf_transfer_sensitivity = NA_real_,
    mean_major_lineage_prediction_score = NA_real_, mean_myeloid_state_prediction_score = NA_real_, source = "Frozen Step131 original 14/28 result (not recomputed)", stringsAsFactors = FALSE
    )
  specimen_tab <- rbind(naive_spec, old_spec)
  specimen_tab$stage <- factor(specimen_tab$stage, levels = c("Naive", "14dpi", "28dpi"))
  specimen_tab <- specimen_tab[order(specimen_tab$stage, specimen_tab$gsm), , drop = FALSE]
  specimen_tab$stage <- as.character(specimen_tab$stage)
  stopif(nrow(specimen_tab) == 9L && setequal(specimen_tab$gsm, ALL_EXPECTED_GSM), "Combined specimen table must contain exactly nine expected GSMs")
  write_csv(specimen_tab, "09_145_specimen_Clec7a_high_proportions_v01.csv")

  stage_summary <- do.call(rbind, lapply(c("Naive", "14dpi", "28dpi"), function(st) {
    x <- specimen_tab$proportion[specimen_tab$stage == st]
    data.frame( stage = st, n_specimens = length(x), mean_fraction = mean(x), median_fraction = stats::median(x), min_fraction = min(x), max_fraction = max(x),
      mean_percent = 100 * mean(x), median_percent = 100 * stats::median(x), min_percent = 100 * min(x), max_percent = 100 * max(x),
      note = "Specimen-level descriptive summary; cell counts are not biological replicates", stringsAsFactors = FALSE )
  }))
  write_csv(stage_summary, "10_145_stage_descriptive_summary_v01.csv")

  effects <- rbind( pairwise_effect_row(specimen_tab, "Naive", "14dpi"), pairwise_effect_row(specimen_tab, "Naive", "28dpi"), pairwise_effect_row(specimen_tab, "14dpi", "28dpi") )
  effects$BH_FDR_across_3_pairwise_exact_tests <- stats::p.adjust( effects$exact_two_sided_permutation_P, method = "BH")
  write_csv(effects, "11_145_pairwise_effect_size_exact_permutation_v01.csv")

  omnibus <- exact_three_group_omnibus( specimen_tab$proportion[specimen_tab$stage == "Naive"], specimen_tab$proportion[specimen_tab$stage == "14dpi"],
    specimen_tab$proportion[specimen_tab$stage == "28dpi"] )
  write_csv(omnibus, "12_145_three_group_exact_permutation_v01.csv")

  sensitivity <- naive_spec[, c( "gsm", "n_candidate_myeloid", "n_primary_target", "proportion", "n_candidate_highconf_transfer", "n_target_highconf_transfer",
    "proportion_highconf_transfer_sensitivity", "mean_major_lineage_prediction_score", "mean_myeloid_state_prediction_score" )]
  sensitivity$primary_percent <- 100 * sensitivity$proportion
  sensitivity$highconf_transfer_sensitivity_percent <-
    100 * sensitivity$proportion_highconf_transfer_sensitivity
  sensitivity$absolute_change_percentage_points <-
    sensitivity$highconf_transfer_sensitivity_percent - sensitivity$primary_percent
  sensitivity$note <- paste0( "Sensitivity threshold only: lineage and state prediction.score.max >= ", TRANSFER_SCORE_SENSITIVITY, "; primary analysis does not use this threshold")
  write_csv(sensitivity, "13_145_naive_transfer_confidence_sensitivity_v01.csv")

  make_figures(specimen_tab, effects)

  report <- c( "STEP145 v01 | Reviewer 1 Remaining Major Comment 1", "GSE281205 Young-naive extension using frozen 14/28-dpi reference",
    paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "", "SCIENTIFIC DESIGN:",
    "1. Three Young-naive GSM specimens were added without modifying the six existing 14/28-dpi Step131 results.",
    "2. Naive QC uses nFeature_RNA >=300 and percent.mt <10%, with no upper feature/UMI cutoff.",
    paste0("3. DecontX was run independently per naive GSM using ", DECONTX_PACKAGE, "::decontX, z=NULL."),
    "4. scDblFinder was run per GSM on original integer UMI counts; DecontX counts were not used for doublet calling.",
    "5. Major lineage was transferred from frozen Step123/128; myeloid r0.10 state was transferred from frozen Step130/131.",
    "6. No naive-informed reclustering, marker tuning or target-state threshold optimization was performed.",
    paste0("7. Frozen primary target gate: Microglia + r0.10 state 1 + contamination <", PRIMARY_CONTAM_LT,
    "% + >=", MIN_MICROGLIA_CORE, "/3 Fcrls/Hexb/Cx3cr1 + no Ccr2/Ly6c2 co-detection."),
    "8. Clec7a expression itself is descriptive and is NOT required for primary target membership, matching Step131.",
    "9. Statistical unit is GSM specimen. Pairwise exact permutation tests use all 20 possible 3-vs-3 label allocations.",
    "10. Welch 95% CIs quantify uncertainty of absolute specimen-level mean differences; with n=3/group they are expected to be wide.",
    "11. The three-group exact permutation test is an omnibus specimen-label test, not a continuous trajectory model.",
    "12. This analysis can test enrichment relative to uninjured tissue, but still does NOT independently test acute early-injury emergence.", "", "SPECIMEN-LEVEL PROPORTIONS:",
    capture.output(print(specimen_tab[, c("gsm", "stage", "n_candidate_myeloid", "n_primary_target", "percent", "source")], row.names = FALSE, digits = 5)), "", "STAGE SUMMARY:",
    capture.output(print(stage_summary, row.names = FALSE, digits = 5)), "", "PAIRWISE EFFECTS:", capture.output(print(effects, row.names = FALSE, digits = 5)), "", "OMNIBUS EXACT TEST:",
    capture.output(print(omnibus, row.names = FALSE, digits = 5)), "", "INTERPRETATION GUARDRAIL:", "Do not describe this analysis as validation of a complete early-to-late trajectory.",
    "If naive proportions are lower than injured proportions, the appropriate claim is support for enrichment relative to uninjured tissue within GSE281205.",
    "The 14-to-28-dpi comparison remains the same frozen six-specimen comparison as Step131.", "SUCCESS | STEP145 v01" )
  writeLines(report, file.path(OUT, "14_145_REPORT_v01.txt"), useBytes = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(OUT, "15_145_sessionInfo_v01.txt"), useBytes = TRUE)

  output_files <- list.files(OUT, recursive = TRUE, full.names = TRUE)
  output_files <- output_files[file.exists(output_files) & !dir.exists(output_files)]
  rel <- substring(normalizePath(output_files, winslash = "/", mustWork = FALSE), nchar(normalizePath(OUT, winslash = "/", mustWork = FALSE)) + 2L)
  idx <- data.frame(relative_path = rel, bytes = as.numeric(file.info(output_files)$size), stringsAsFactors = FALSE)

  write_csv(idx, "16_145_OUTPUT_INDEX_v01.csv")

  log_msg("============================================================")
  log_msg("SUCCESS | STEP145 v01")
  log_msg("Naive target proportions: ", paste0(NAIVE_PLAN$gsm, "=", sprintf("%.3f%%", 100 * naive_spec$proportion), collapse = "; "))
  log_msg("Elapsed minutes: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
  log_msg("============================================================")
  invisible(list(specimen = specimen_tab, stage = stage_summary, effects = effects, omnibus = omnibus))
}

tryCatch(run(), error = function(e) {
  log_msg("FAILED | STEP145 v01 | ", conditionMessage(e))
  stop(e)
})
