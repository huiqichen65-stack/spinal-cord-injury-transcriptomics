options(stringsAsFactors = FALSE, scipen = 50, width = 170)

ROOT <- "D:/zz"
STEP <- "143"
VER <- "v01"
PREFIX <- "143_de_novo_gene_selection_network_sensitivity"
OUT <- file.path(ROOT, "results", PREFIX)
CP <- file.path(OUT, "network_checkpoints_v01", "A01_whole_transcriptome_mean_bin_dispersion_3000")
LOG <- file.path(OUT, paste0(PREFIX, "_v01.log"))

NETWORK_GENE_CAP <- 3000L
MIN_DETECTED_CELLS <- 25L
N_MEAN_BINS <- 20L
MIN_ACTUAL_CHANGED_GENES <- 100L
MIN_NEW_OUTSIDE_OLD_5000 <- 1L
MIN_SHARED_NONKO <- 500L
WORKERS <- 6L
WT_SEED <- 9101L
N_NET <- 10L
N_CELLS <- 500L
Q_NETWORK <- 0.9
K_TENSOR <- 3L
N_COMPONENTS <- 3L
MA_DIMENSIONS <- 2L
TOP_K <- 50L
FDR_CUTOFF <- 0.05
USE_FGSEA <- TRUE
RUN_NETWORK <- TRUE

TFS <- c("Mef2a", "Nr3c1", "Pura", "Ikzf1", "Bach1", "Fli1", "Zfp148", "Ep300", "Bclaf1", "Chd1")
FUNCTIONAL20 <- c("Clec7a", "Trem2", "Cst7", "Cx3cr1", "Axl", "C1qa", "Apoe", "Lpl", "Abca1", "Apoc1",
  "Apoc4", "Npc2", "Soat1", "Tyrobp", "Csf1r", "Itgam", "Hexb", "Ctsd", "Lamp1", "Lamp2")

P <- list(full_RNA = file.path(ROOT, "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds"),
  step91_CP01 = file.path(ROOT, "results", "91_scTenifoldKnk_top5_virtual_KO", "checkpoints_v04", "CP01_counts_network_v04.rds"),
  step91_genes = file.path(ROOT, "results", "91_scTenifoldKnk_top5_virtual_KO", "checkpoints_v04",
    "CP01_gene_selection_v04.csv.gz"), step90_top10 = file.path(ROOT, "results", "90_pySCENIC_top10_TF_network_signature_correlation",
    "10_90_top10_TF_evidence_summary_v05.csv"), step20_lipid = file.path(ROOT, "results", "20_AUCell_lipid_metabolism_scoring",
    "20_AUCell_lipid_gene_set_manifest_v03.csv"), B00_genes = file.path(ROOT, "results", "135_TF_candidate_comparison_resampling_v02",
    "network_checkpoints_v02", "B00_seed9101_3000", "00_gene_list_v02.txt"), B00_cells = file.path(ROOT,
    "results", "135_TF_candidate_comparison_resampling_v02", "network_checkpoints_v02", "B00_seed9101_3000",
    "00_cell_list_v02.txt"), B00_signature = file.path(ROOT, "results", "135_TF_candidate_comparison_resampling_v02",
    "network_checkpoints_v02", "B00_seed9101_3000", "00_parameter_signature_v02.txt"))
B00_CP <- dirname(P$B00_genes)

stopf <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
if (!dir.exists(ROOT)) stopf("Project root is missing: ", ROOT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(CP, recursive = TRUE, showWarnings = FALSE)
logline <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = LOG, append = TRUE)
  invisible(line)
}
readcsv <- function(path) {
  if (!file.exists(path))
    stopf("Required file missing: ", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE))
    return(utils::read.csv(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
writecsv <- function(x, name) {
  fn <- file.path(OUT, name)
  if (grepl("\\.gz$", fn)) {
    con <- gzfile(fn, "wt")
    on.exit(close(con), add = TRUE)
    utils::write.csv(x, con, row.names = FALSE, na = "")
  } else utils::write.csv(x, fn, row.names = FALSE, na = "")
  logline("Saved ", name, " (", nrow(x), " rows)")
  invisible(fn)
}
atomicRDS <- function(x, path, compress = TRUE) {
  if (file.exists(path))
    stopf("Refusing to replace existing checkpoint: ", path)
  tmp <- paste0(path, ".incomplete")
  if (file.exists(tmp))
    stopf("Incomplete checkpoint exists; inspect: ", tmp)
  saveRDS(x, tmp, compress = compress)
  if (!file.rename(tmp, path))
    stopf("Failed to commit checkpoint: ", path)
  invisible(path)
}
atomicKO <- function(x, path) {
  if (file.exists(path))
    stopf("Refusing to overwrite an existing KO file: ", path)
  tmp <- paste0(path, ".incomplete")
  if (file.exists(tmp))
    stopf("Incomplete KO checkpoint exists: ", tmp)
  con <- gzfile(tmp, "wt")
  tryCatch(utils::write.csv(x, con, row.names = FALSE, na = ""), finally = close(con))
  if (!file.rename(tmp, path))
    stopf("Could not commit KO file: ", path)
  invisible(path)
}
checkcols <- function(x, cols, who) {
  missing <- setdiff(cols, names(x))
  if (length(missing))
    stopf(who, ": missing columns: ", paste(missing, collapse = ", "))
}
finite_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20L || length(unique(x[ok])) < 3L || length(unique(y[ok])) < 3L)
    return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}
file_signature <- function(path) {
  f <- file.info(path)
  paste0(normalizePath(path, winslash = "/"), "::bytes=", f$size, "::mtime=", format(f$mtime, "%Y-%m-%dT%H:%M:%S"))
}
write_status <- function(state, details = "") {
  x <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), status = state, details = details,
    stringsAsFactors = FALSE)
  write.csv(x, file.path(OUT, "01_143_run_status_v01.csv"), row.names = FALSE)
  logline(state, if (nzchar(details))
    paste0(": ", details) else "")
}
if (!file.exists(LOG)) writeLines(character(0), LOG)
logline("Step 143 v01 start; no zip input; ROOT=", ROOT)

input_files <- c(unlist(P, use.names = TRUE), setNames(file.path(B00_CP, paste0("KO_", TFS, "_v02.csv.gz")),
  paste0("B00_KO_", TFS)))
input_audit <- data.frame(name = names(input_files), path = unname(input_files), exists = file.exists(input_files),
  size_bytes = as.numeric(file.info(input_files)$size), stringsAsFactors = FALSE)
writecsv(input_audit, "00_143_INPUT_AUDIT_v01.csv")
if (any(!input_audit$exists)) stopf("Missing inputs in D:/zz: ", paste(input_audit$name[!input_audit$exists],
  collapse = ", "))
if (!RUN_NETWORK) stopf("RUN_NETWORK=FALSE; this script's objective requires building an alternative WT network.")
for (pkg in c("Matrix", "SeuratObject", "scTenifoldNet", "scTenifoldKnk")) if (!requireNamespace(pkg,
  quietly = TRUE)) stopf("Required R package missing: ", pkg)
if (USE_FGSEA && !requireNamespace("fgsea", quietly = TRUE)) stopf("Package fgsea is missing")
if (as.character(utils::packageVersion("scTenifoldNet")) != "1.3" || as.character(utils::packageVersion("scTenifoldKnk")) !=
  "1.0.3") stopf("Use scTenifoldNet 1.3 and scTenifoldKnk 1.0.3 (Step 91/135).")
if (N_NET != 10L || N_CELLS != 500L || Q_NETWORK != 0.9 || K_TENSOR != 3L) stopf("paired-network parameters changed; do not compare as gene-selection-only sensitivity.")
if (!is.finite(WORKERS) || WORKERS < 1) stopf("WORKERS must be >=1")

bgenes <- readLines(P$B00_genes, warn = FALSE, encoding = "UTF-8")
bcells <- readLines(P$B00_cells, warn = FALSE, encoding = "UTF-8")
bsig <- readLines(P$B00_signature, warn = FALSE, encoding = "UTF-8")
if (length(bgenes) != NETWORK_GENE_CAP || length(bcells) < N_CELLS || anyDuplicated(bgenes) || anyDuplicated(bcells)) stopf("135 v02 B00 gene/cell list malformed")
if (!all(TFS %in% bgenes)) stopf("Ten TFs missing from the B00 network")
for (token in c("seed=9101", "nNet=10", "nCells=500", "q=0.9", "K=3")) {
  if (!any(grepl(token, bsig, fixed = TRUE)))
    stopf("B00 parameters need manual inspection. Signature lacks: ", token)
}

evi <- readcsv(P$step90_top10)
checkcols(evi, c("TF", "RSS_rank"), "Step 90 top-ten")
if (!setequal(as.character(evi$TF), TFS)) stopf("Step 90 ten-TF roster differs from the preregistered comparison")

write_status("INPUTS_CHECKED", "All 10 B00 KOs, B00 signature/cells, original Seurat and Step91 CP01 exist")

cp01 <- readRDS(P$step91_CP01)
if (!inherits(cp01, "Matrix") && !is.matrix(cp01)) stopf("Step91 CP01 is not a matrix")
if (!identical(colnames(cp01), bcells)) stopf("B00 cell order != Step91 CP01 cell order; paired analysis invalid")
if (nrow(cp01) != 5000L) stopf("Expected Step91 5000-gene matrix")
old5000 <- rownames(cp01)
if (!all(bgenes %in% old5000)) stopf("B00 genes not contained in old Step91 universe")

logline("Reading full original Step20 RNA Seurat object (not 5k matrix) ...")
obj <- readRDS(P$full_RNA)
if (!"RNA" %in% names(obj@assays)) stopf("Full Seurat object lacks RNA assay")
all_counts <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
if (!inherits(all_counts, "dgCMatrix")) all_counts <- as(all_counts, "dgCMatrix")
if (anyDuplicated(rownames(all_counts)) || anyDuplicated(colnames(all_counts))) stopf("Duplicate gene or cell labels in full RNA assay")
if (!all(bcells %in% colnames(all_counts))) stopf("Step91/B00 cell IDs missing from full original RNA counts")
all_counts <- all_counts[, bcells, drop = FALSE]
rm(obj)
gc(verbose = FALSE)
if (nrow(all_counts) <= 5000L) stopf("Full RNA matrix does not contain >5000 genes")
set.seed(14301L)
testgenes <- sample(old5000, min(30L, length(old5000)), replace = FALSE)
delta <- all_counts[testgenes, , drop = FALSE] - cp01[testgenes, , drop = FALSE]
if (Matrix::nnzero(delta) != 0L) stopf("The two matrices have different raw counts; cannot isolate gene selection")
rm(delta, cp01)
gc(verbose = FALSE)

prevalence <- as.numeric(Matrix::rowSums(all_counts != 0))
keep <- prevalence >= MIN_DETECTED_CELLS
if (sum(keep) < NETWORK_GENE_CAP) stopf("Too few eligible genes after >=25-cell detection")
X <- all_counts[keep, , drop = FALSE]
rm(all_counts)
gc(verbose = FALSE)
libsize <- as.numeric(Matrix::colSums(X))
if (any(!is.finite(libsize) | libsize <= 0)) stopf("Invalid cell library size")
L <- X %*% Matrix::Diagonal(x = 10000/libsize)
L@x <- log1p(L@x)
mu <- as.numeric(Matrix::rowMeans(L))
vv <- pmax(as.numeric(Matrix::rowMeans(L^2)) - mu^2, 0)
rm(L)
gc(verbose = FALSE)
if (any(!is.finite(mu) | !is.finite(vv))) stopf("Nonfinite means/variances in HVG ranking")
rank_order <- order(log1p(mu), rownames(X))

bin <- integer(length(mu))
bin[rank_order] <- pmin(N_MEAN_BINS, ceiling(seq_along(rank_order)/(length(rank_order)/N_MEAN_BINS)))
score <- rep(NA_real_, length(mu))
for (b in sort(unique(bin))) {
  ii <- which(bin == b)
  center <- stats::median(vv[ii])
  scale <- stats::mad(vv[ii], center = center, constant = 1.4826)
  if (!is.finite(scale) || scale < 1e-08)
    scale <- max(stats::IQR(vv[ii])/1.349, 1e-08)
  score[ii] <- (vv[ii] - center)/scale
}
if (any(!is.finite(score))) stopf("Nonfinite adjusted-dispersion score")
fullstats <- data.frame(gene = rownames(X), detected_cells = prevalence[keep], mean_logCPM = mu, logCPM_variance = vv,
  mean_bin = bin, mean_bin_dispersion_z = score, stringsAsFactors = FALSE)
ord <- order(-fullstats$mean_bin_dispersion_z, fullstats$gene)
ranked <- fullstats$gene[ord]
if (!all(TFS %in% fullstats$gene)) stopf("Candidate TF(s) fail gene detection QC")

newgenes <- unique(c(TFS, head(setdiff(ranked, TFS), NETWORK_GENE_CAP - length(TFS))))
if (length(newgenes) != NETWORK_GENE_CAP || anyDuplicated(newgenes)) stopf("Failed to select 3000 unique genes")
fullstats$selected_A01 <- fullstats$gene %in% newgenes
fullstats$selected_B00 <- fullstats$gene %in% bgenes
fullstats$in_original5000 <- fullstats$gene %in% old5000
fullstats$is_TF <- fullstats$gene %in% TFS
fullstats$is_functional20 <- fullstats$gene %in% FUNCTIONAL20
writecsv(fullstats, "02_143_whole_transcriptome_gene_selection_v01.csv.gz")
writeLines(newgenes, file.path(OUT, "03_143_selected_gene_list_v01.txt"), useBytes = TRUE)

only_B00 <- setdiff(bgenes, newgenes)
only_A01 <- setdiff(newgenes, bgenes)
new_outside_5000 <- setdiff(newgenes, old5000)
shared <- intersect(bgenes, newgenes)
selection_audit <- data.frame(item = c("full_RNA_genes_all", "RNA_genes_passing_25_cells", "frozen_Step91_genes",
  "B00_genes", "A01_genes", "shared_B00_A01", "A01_only", "B00_only", "A01_outside_old5000", "gene_set_jaccard",
  "A01_functional20_present", "B00_functional20_present", "forced_TFs_present_A01", "same_cell_ids",
  "seed_matched", "q_matched", "K_matched", "nNet_matched", "nCells_matched"), value = as.character(c(length(prevalence),
  nrow(X), length(old5000), length(bgenes), length(newgenes), length(shared), length(only_A01), length(only_B00),
  length(new_outside_5000), length(shared)/length(union(bgenes, newgenes)), sum(FUNCTIONAL20 %in% newgenes),
  sum(FUNCTIONAL20 %in% bgenes), sum(TFS %in% newgenes), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)), stringsAsFactors = FALSE)
writecsv(selection_audit, "04_143_gene_selection_overlap_audit_v01.csv")
membership <- data.frame(gene = sort(union(bgenes, newgenes)), in_B00 = sort(union(bgenes, newgenes)) %in%
  bgenes, in_A01 = sort(union(bgenes, newgenes)) %in% newgenes, outside_old5000 = !(sort(union(bgenes,
  newgenes)) %in% old5000), in_functional20 = sort(union(bgenes, newgenes)) %in% FUNCTIONAL20, stringsAsFactors = FALSE)
writecsv(membership, "05_143_selected_gene_membership_v01.csv")
if (length(only_A01) < MIN_ACTUAL_CHANGED_GENES || length(new_outside_5000) < MIN_NEW_OUTSIDE_OLD_5000) stopf("Gene-selection comparison NOT INFORMATIVE: A01-only=",
  length(only_A01), "; outside original 5000=", length(new_outside_5000), ". No network computed. Inspect gene ranking and predeclared QC.")
logline("Gene selection changed: ", length(only_A01), "/3000 new genes; ", length(new_outside_5000),
  " outside 5000; overlap=", length(shared))
write_status("GENE_SELECTION_PASSED", paste(length(only_A01), "A01-only;", length(new_outside_5000),
  "outside old5000"))

Xnet <- X[newgenes, , drop = FALSE]
if (!identical(rownames(Xnet), newgenes) || !identical(colnames(Xnet), bcells)) stopf("New network count matrix row/column order mismatch")
rm(X)
gc(verbose = FALSE)

gene_cp <- file.path(CP, "00_gene_list_v01.txt")
cell_cp <- file.path(CP, "00_cell_list_v01.txt")
signature_cp <- file.path(CP, "00_parameter_signature_v01.txt")
ensemble_cp <- file.path(CP, "CP02_WT_ensemble_v01.rds")
wt_cp <- file.path(CP, "CP03_WT_network_v01.rds")
ko_paths <- setNames(file.path(CP, paste0("KO_", TFS, "_v01.csv.gz")), TFS)
sig <- c("step=143", "version=v01", "selection=mean_binned_logCPM_dispersion_whole_transcriptome", "score_bins=20",
  "gene_detected_cells>=25", "forced_candidates=tenTF_only", "forced_20_functional_genes=FALSE", "count_matrix=Step20_RNA_counts",
  paste0("source_RDS=", file_signature(P$full_RNA)), paste0("source_CP01=", file_signature(P$step91_CP01)),
  paste0("seed=", WT_SEED), paste0("q=", format(Q_NETWORK, trim = TRUE)), paste0("K=", K_TENSOR), paste0("nNet=",
    N_NET), paste0("nCells=", N_CELLS), paste0("nComp=", N_COMPONENTS), "scaleScores=TRUE", "symmetric=FALSE",
  "td_maxIter=1000", "td_maxError=1e-5", "td_nDecimal=3", "manifold_d=2", "manifold_KO_seed=9101+index*103",
  paste0("scTenifoldNet=", as.character(utils::packageVersion("scTenifoldNet"))), paste0("scTenifoldKnk=",
    as.character(utils::packageVersion("scTenifoldKnk"))), paste0("genes=", paste(newgenes, collapse = ";")),
  paste0("cells=", paste(bcells, collapse = ";")))
if (file.exists(signature_cp)) {
  oldsig <- readLines(signature_cp, warn = FALSE)
  if (!identical(oldsig, sig))
    stopf("Checkpoint signature mismatch; DO NOT reuse ", CP)
  if (!file.exists(gene_cp) || !file.exists(cell_cp) || !identical(readLines(gene_cp, warn = FALSE),
    newgenes) || !identical(readLines(cell_cp, warn = FALSE), bcells))
    stopf("Existing gene/cell checkpoint files disagree with the frozen signature")
} else {
  if (any(file.exists(c(ensemble_cp, wt_cp, ko_paths))))
    stopf("Orphaned network checkpoints without signature in ", CP)
  writeLines(newgenes, gene_cp, useBytes = TRUE)
  writeLines(bcells, cell_cp, useBytes = TRUE)
  writeLines(sig, signature_cp, useBytes = TRUE)
}

status_detail <- function(n) paste0("scenario=A01, finished TF KOs=", n, "/10")
if (file.exists(wt_cp)) {
  WT <- as.matrix(readRDS(wt_cp))
  if (!identical(rownames(WT), newgenes) || !identical(colnames(WT), newgenes))
    stopf("Cached WT network has wrong gene labels/order")
  logline("Resuming from saved CP03 WT network")
} else {
  if (file.exists(ensemble_cp)) {
    ens <- readRDS(ensemble_cp)
    logline("Resuming from saved CP02 WT ensemble")
  } else {
    write_status("NETWORK_RUNNING", "A01 makeNetworks; one full 3000-gene network")
    set.seed(WT_SEED)
    ens <- scTenifoldNet::makeNetworks(X = Xnet, q = Q_NETWORK, nNet = N_NET, nCells = N_CELLS, scaleScores = TRUE,
      symmetric = FALSE, nComp = N_COMPONENTS, nCores = min(WORKERS, max(1L, parallel::detectCores() -
        1L)))
    atomicRDS(ens, ensemble_cp, compress = FALSE)
    logline("CP02 WT ensemble saved")
  }
  set.seed(WT_SEED)
  td <- scTenifoldNet::tensorDecomposition(xList = ens, K = K_TENSOR, maxError = 1e-05, maxIter = 1000L,
    nDecimal = 3L)
  direction_fn <- getFromNamespace("strictDirection", "scTenifoldKnk")
  WT <- as.matrix(direction_fn(td$X, lambda = 0))
  diag(WT) <- 0
  WT <- t(WT)
  if (!setequal(rownames(WT), newgenes) || !setequal(colnames(WT), newgenes))
    stopf("A01 tensor-derived WT gene list differs from preselected 3000")
  WT <- WT[newgenes, newgenes, drop = FALSE]
  atomicRDS(Matrix::Matrix(WT, sparse = TRUE), wt_cp, compress = "gzip")
  rm(ens, td)
  gc(verbose = FALSE)
  logline("CP03 WT network saved")
}
rm(Xnet)
gc(verbose = FALSE)
if (!identical(dim(WT), c(NETWORK_GENE_CAP, NETWORK_GENE_CAP))) stopf("Wrong WT network dimensions")
if (!all(TFS %in% rownames(WT))) stopf("A01 WT lost candidate TFs")

regfun <- getFromNamespace("dRegulation", "scTenifoldKnk")
newdr <- list()
for (i in seq_along(TFS)) {
  tf <- TFS[i]
  path <- ko_paths[[tf]]
  if (file.exists(path)) {
    d <- readcsv(path)
    logline("Reused completed A01 KO: ", tf)
  } else {
    if (sum(WT[tf, ] != 0) == 0L)
      stopf("No outgoing WT edges for TF ", tf)
    KO <- WT
    KO[tf, ] <- 0
    set.seed(9101L + i * 103L)
    ma <- scTenifoldNet::manifoldAlignment(WT, KO, d = MA_DIMENSIONS, nCores = min(WORKERS, max(1L,
      parallel::detectCores() - 1L)))
    d <- as.data.frame(regfun(ma, tf))
    if (!"gene" %in% names(d))
      d$gene <- rownames(d)
    d$TF <- tf
    d$scenario <- "A01_whole_transcriptome_mean_bin_dispersion_3000"
    atomicKO(d, path)
    rm(KO, ma)
    gc(verbose = FALSE)
  }
  checkcols(d, c("gene", "distance", "Z", "p.adj"), tf)
  if (nrow(d) != NETWORK_GENE_CAP || anyDuplicated(as.character(d$gene)) || !setequal(as.character(d$gene),
    newgenes))
    stopf("Invalid KO output for ", tf, "; must contain exactly 3000 unique genes")
  newdr[[tf]] <- d
  write_status("NETWORK_RUNNING", status_detail(i))
}
rm(WT)
gc(verbose = FALSE)
write_status("NETWORK_COMPLETE", "All ten TF virtual knockouts completed in A01")

base_dr <- lapply(TFS, function(tf) readcsv(file.path(B00_CP, paste0("KO_", tf, "_v02.csv.gz"))))
names(base_dr) <- TFS
row_results <- list()
summary_rows <- list()
comp_z <- list()
for (tf in TFS) {
  bx <- base_dr[[tf]]
  ax <- newdr[[tf]]
  checkcols(bx, c("gene", "distance", "Z", "p.adj"), paste0("B00/", tf))
  checkcols(ax, c("gene", "distance", "Z", "p.adj"), paste0("A01/", tf))
  bx$gene <- as.character(bx$gene)
  ax$gene <- as.character(ax$gene)
  bx <- bx[bx$gene != tf & is.finite(suppressWarnings(as.numeric(bx$Z))), , drop = FALSE]
  ax <- ax[ax$gene != tf & is.finite(suppressWarnings(as.numeric(ax$Z))), , drop = FALSE]
  g <- intersect(bx$gene, ax$gene)
  if (length(g) < MIN_SHARED_NONKO)
    stopf(tf, ": <", MIN_SHARED_NONKO, " shared evaluable nonKO genes")
  bz <- as.numeric(bx$Z[match(g, bx$gene)])
  az <- as.numeric(ax$Z[match(g, ax$gene)])
  br <- rank(-bz, ties.method = "average")
  ar <- rank(-az, ties.method = "average")
  btop <- g[order(-bz, g)][seq_len(min(TOP_K, length(g)))]
  atop <- g[order(-az, g)][seq_len(min(TOP_K, length(g)))]
  row_results[[tf]] <- data.frame(TF = tf, reference = "135v02_B00", comparison = "143v01_A01", n_shared_nonKO = length(g),
    spearman_common_Z = finite_cor(bz, az), spearman_common_rank = finite_cor(br, ar), topK = length(btop),
    topK_overlap = length(intersect(btop, atop)), topK_jaccard = length(intersect(btop, atop))/length(union(btop,
      atop)), n_baseline_nonKO_FDR05_all = sum(as.numeric(bx$p.adj) < FDR_CUTOFF, na.rm = TRUE),
    n_alternative_nonKO_FDR05_all = sum(as.numeric(ax$p.adj) < FDR_CUTOFF, na.rm = TRUE), stringsAsFactors = FALSE)
  comp_z[[tf]] <- list(B00 = setNames(bz, g), A01 = setNames(az, g))
  summary_rows[[tf]] <- data.frame(TF = tf, n_B00_genes = nrow(bx), n_A01_genes = nrow(ax), n_common_nonKO = length(g),
    n_B00_FDR05_common = sum(as.numeric(bx$p.adj[match(g, bx$gene)]) < FDR_CUTOFF, na.rm = TRUE),
    n_A01_FDR05_common = sum(as.numeric(ax$p.adj[match(g, ax$gene)]) < FDR_CUTOFF, na.rm = TRUE),
    stringsAsFactors = FALSE)
}
rank_summary <- do.call(rbind, row_results)
writecsv(do.call(rbind, summary_rows), "06_143_tenTF_KO_nonKO_summary_v01.csv")
writecsv(rank_summary, "07_143_common_universe_ranking_stability_v01.csv")

fg_out <- data.frame()
fg_compare <- data.frame()
if (USE_FGSEA) {
  lipid <- readcsv(P$step20_lipid)
  checkcols(lipid, c("gene_set", "gene"), "Step20 lipid manifest")
  lipid <- lipid[!is.na(lipid$gene_set) & !is.na(lipid$gene) & nzchar(lipid$gene_set) & nzchar(lipid$gene),
    , drop = FALSE]
  pways <- lapply(split(as.character(lipid$gene), as.character(lipid$gene_set)), unique)
  fglist <- list()
  idx <- 0L
  for (mode in c("B00", "A01")) for (i in seq_along(TFS)) {
    tf <- TFS[i]
    z <- comp_z[[tf]][[mode]]
    z <- sort(z[is.finite(z)], decreasing = TRUE)
    eligible <- lapply(pways, function(g) intersect(g, names(z)))
    eligible <- eligible[vapply(eligible, length, integer(1)) >= 5L & vapply(eligible, length, integer(1)) <=
      500L]
    if (!length(eligible))
      next
    set.seed(14300L + i)
    fg <- fgsea::fgseaMultilevel(pathways = eligible, stats = z, minSize = 5L, maxSize = 500L, eps = 0)
    if (!nrow(fg))
      next
    fg <- as.data.frame(fg)
    fg$leading_edge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    fg$leadingEdge <- NULL
    fg$TF <- tf
    fg$case <- mode
    fg$rank_interpretation <- "distance-standardized virtual KO network Z, NOT gene expression direction"
    idx <- idx + 1L
    fglist[[idx]] <- fg
  }
  if (length(fglist)) {
    fg_out <- do.call(rbind, fglist)
    fg_out$FDR_global_10TF <- NA_real_
    for (mode in c("B00", "A01")) {
      ii <- which(fg_out$case == mode & is.finite(fg_out$pval))
      fg_out$FDR_global_10TF[ii] <- p.adjust(fg_out$pval[ii], method = "BH")
    }
    writecsv(fg_out, "08_143_lipid_enrichment_common_universe_all_v01.csv")
    fba <- fg_out[fg_out$case == "B00", c("TF", "pathway", "size", "NES", "pval", "FDR_global_10TF",
      "leading_edge")]
    faa <- fg_out[fg_out$case == "A01", c("TF", "pathway", "size", "NES", "pval", "FDR_global_10TF",
      "leading_edge")]
    names(fba)[-(1:2)] <- paste0(names(fba)[-(1:2)], "_B00")
    names(faa)[-(1:2)] <- paste0(names(faa)[-(1:2)], "_A01")
    fg_compare <- merge(fba, faa, by = c("TF", "pathway"), all = TRUE, sort = TRUE)
    fg_compare$NES_same_sign <- with(fg_compare, ifelse(is.finite(NES_B00) & is.finite(NES_A01),
      sign(NES_B00) == sign(NES_A01), NA))
    fg_compare$FDR05_both <- with(fg_compare, is.finite(FDR_global_10TF_B00) & is.finite(FDR_global_10TF_A01) &
      FDR_global_10TF_B00 < FDR_CUTOFF & FDR_global_10TF_A01 < FDR_CUTOFF)
    writecsv(fg_compare, "09_143_lipid_enrichment_paired_comparison_v01.csv")
  } else {
    logline("WARNING: no evaluable common-universe lipid pathways")
  }
}

qc <- data.frame(item = c("all_inputs_present", "Step90_top10_roster_matches", "same_cells_as_B00", "same_seed_and_network_parameters",
  "input_raw_counts_match_Step91_CP01", "network_selected_3000_unique", "A01_has_10_TFs", "A01_escaping_frozen_5000",
  "A01_changed_at_least_100_genes", "all_10_KOs_complete", "all_TF_shared_universe_at_least_500", "gene_selection_does_not_force_functional20",
  "paired_lipid_enrichment_available", "whole_transcriptome_gene_selection_performed", "separate_functional_panel_effect_identified"),
  status = c("PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", if (length(new_outside_5000) >
    0) "PASS" else "FAIL", if (length(only_A01) >= MIN_ACTUAL_CHANGED_GENES) "PASS" else "FAIL",
    if (length(newdr) == length(TFS)) "PASS" else "FAIL", if (all(rank_summary$n_shared_nonKO >=
      MIN_SHARED_NONKO)) "PASS" else "FAIL", "PASS", if (nrow(fg_out)) "PASS" else "NOT_RUN_OR_NO_ELIGIBLE_SET",
    "PASS", "NOT_TESTED_JOINT_STRATEGY_CHANGE"), detail = c("All original project inputs audited",
    "Step90 top-ten TF list", "Exactly the same B00 cells as Step91 CP01", "q=.90,K=3,nNet=10,nCells=500,seed=9101",
    "30 random genes have identical raw counts for all cells", paste(length(newgenes), "selected from",
      nrow(fullstats), "eligible genes"), "Ten candidate TFs retained for KO", paste(length(new_outside_5000),
      "genes outside old 5000"), paste(length(only_A01), "new vs B00"), "Ten complete KO files with 3000 unique genes",
    "Compared on intersection excluding each knocked-out TF", "Only ten candidate TFs mandatorily retained",
    "FGSEA within paired common-gene universes; BH across 10TF x program within case", "Alternative mean-binned dispersion computed from full Step20 RNA counts",
    "Different ranking and panel forcing changed together; not an isolated forced-panel contrast"),
  stringsAsFactors = FALSE)
writecsv(qc, "10_143_FINAL_QC_v01.csv")
writecsv(data.frame(parameter = c("version_R", "scTenifoldNet", "scTenifoldKnk", "workers", "seed", "gene_cap",
  "gene_min_detection", "mean_bins", "network_nNet", "network_nCells", "network_q", "tensor_K", "ma_d",
  "fgsea_common_universe", "comparand"), value = c(R.version.string, as.character(utils::packageVersion("scTenifoldNet")),
  as.character(utils::packageVersion("scTenifoldKnk")), as.character(WORKERS), as.character(WT_SEED),
  as.character(NETWORK_GENE_CAP), as.character(MIN_DETECTED_CELLS), as.character(N_MEAN_BINS), as.character(N_NET),
  as.character(N_CELLS), as.character(Q_NETWORK), as.character(K_TENSOR), as.character(MA_DIMENSIONS),
  as.character(USE_FGSEA), "135 v02 B00_seed9101_3000"), stringsAsFactors = FALSE), "11_143_environment_parameters_v01.csv")
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, "12_143_REPORT_v01.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "13_143_sessionInfo_v01.txt"), useBytes = TRUE)
index_files <- list.files(OUT, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
index <- data.frame(relative_path = substring(index_files, nchar(OUT) + 2L), size_bytes = as.numeric(file.info(index_files)$size),
  stringsAsFactors = FALSE)
writecsv(index, "14_143_OUTPUT_INDEX_v01.csv")
write_status("COMPLETE", "143 v01 done; inspect FINAL_QC and REPORT")
logline("STEP 143 v01 COMPLETE")
