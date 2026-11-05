source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1)
set.seed(1234L)

PROJECT_DIR <- PROJECT_ROOT
INPUT <- file.path(PROJECT_DIR, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
OUTPUT <- file.path(PROJECT_DIR, "results", "104_state_stability_v08")
SOURCES <- c("17702045", "GSE162610", "GSE196928")
MAIN_STAGES <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")
FULL_STAGES <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")
FOCUS_LABEL <- "Clec7a+ lipid-associated Microglia"
MARKERS <- c("Lpl", "Apoe", "Trem2", "Cst7", "Axl", "Apoc4")
SEED <- 1234L
N_HVG <- 3000L
N_PCS_COMPUTE <- 50L
N_PCS_USE <- 30L
RESOLUTION <- 0.1
HARMONY_MAX_ITER <- 20L
MAX_ZERO_COUNT_FRAC <- 0.02
MAX_COUNT_MISMATCH_FRAC <- 0.01

PROXY_THRESHOLDS <- data.frame(rule = c("primary", "stringent"), min_cluster_cells = c(50L, 50L), min_mean_module_score = c(0.58,
  0.62), min_fraction_module_high = c(0.28, 0.35), stringsAsFactors = FALSE)
RUN_SINGLE_SOURCE <- TRUE
RUN_UMAP_FULL <- FALSE
EXPORT_CELL_ASSIGNMENTS <- TRUE
MAKE_ZIP <- TRUE
ZIP_LIMIT_MB <- 30

pkgs <- c("Seurat", "SeuratObject", "Matrix", "harmony", "ggplot2")
not_installed <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed)) stop("Missing R packages: ", paste(not_installed, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})
if (!file.exists(INPUT)) stop("Input RDS not found: ", INPUT)
dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUTPUT, "104_state_stability_v08.log")
if (file.exists(LOG)) unlink(LOG)
logit <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = LOG, append = TRUE)
}
write_csv <- function(x, stem) {
  path <- file.path(OUTPUT, paste0(stem, "_v08.csv"))
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  logit("Wrote ", basename(path), ": ", nrow(x), " rows")
  invisible(path)
}
bind_rows <- function(x) {
  x <- Filter(function(a) is.data.frame(a) && nrow(a) > 0, x)
  if (!length(x))
    return(data.frame())
  ans <- do.call(rbind, x)
  rownames(ans) <- NULL
  ans
}
ari <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- a[ok]
  b <- b[ok]
  n <- length(a)
  if (n < 2)
    return(NA_real_)
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1)/2
  observed <- sum(choose2(tab))
  r <- sum(choose2(rowSums(tab)))
  c <- sum(choose2(colSums(tab)))
  expected <- r * c/choose2(n)
  maximum <- (r + c)/2
  if (abs(maximum - expected) < 1e-10)
    return(NA_real_)
  (observed - expected)/(maximum - expected)
}
jaccard <- function(a, b) {
  union <- sum(a | b)
  if (union == 0)
    return(NA_real_)
  sum(a & b)/union
}
rate <- function(x) if (!length(x)) NA_real_ else 100 * mean(x)

logit("104 v08 started ", as.character(Sys.time()))
logit("R ", as.character(getRversion()), "; Seurat ", as.character(packageVersion("Seurat")), "; harmony ",
  as.character(packageVersion("harmony")))
logit("Input: ", INPUT)
logit("Output: ", OUTPUT)

obj0 <- readRDS(INPUT)
if (!inherits(obj0, "Seurat")) stop("Input is not a Seurat object")
if (!"RNA" %in% Assays(obj0)) stop("RNA assay missing")
DefaultAssay(obj0) <- "RNA"
md <- obj0[[]]
required <- c("source_dataset", "timepoint_main", "myeloid_state")
if (length(setdiff(required, colnames(md)))) {
  stop("Missing metadata: ", paste(setdiff(required, names(md)), collapse = ", "))
}
md$cell <- rownames(md)
md$source_dataset <- as.character(md$source_dataset)
md$timepoint_main <- as.character(md$timepoint_main)
md$myeloid_state <- as.character(md$myeloid_state)
if (anyNA(md$source_dataset) || any(!md$source_dataset %in% SOURCES)) {
  stop("Source labels missing or unexpected")
}
if (anyNA(md$timepoint_main) || any(!md$timepoint_main %in% FULL_STAGES)) {
  stop("Stage labels missing or unexpected")
}
md$frozen_focus <- md$myeloid_state == FOCUS_LABEL
if (!any(md$frozen_focus, na.rm = TRUE)) stop("Original focus state absent")
if ("treatment" %in% colnames(md)) {
  tx <- unique(as.character(md$treatment[md$source_dataset == "GSE196928"]))
  logit("GSE196928 treatment labels: ", paste(tx, collapse = ";"))
  if (any(grepl("PLX5622", tx, ignore.case = TRUE))) {
    stop("PLX5622 unexpectedly found in retained atlas")
  }
}
if (nrow(md) != 58533L) warning("Final cell number differs from prior 58,533: ", nrow(md))

obj0 <- tryCatch(JoinLayers(obj0, assay = "RNA"), error = function(e) stop("JoinLayers failed: ", conditionMessage(e)))
counts <- tryCatch(GetAssayData(obj0, assay = "RNA", layer = "counts"), error = function(e) stop("RNA counts inaccessible: ",
  conditionMessage(e)))
if (!setequal(colnames(counts), md$cell) || ncol(counts) != nrow(md)) {
  stop("RNA counts/cell identifiers mismatch")
}
counts <- counts[, md$cell, drop = FALSE]
if (!all(MARKERS %in% rownames(counts))) {
  stop("Common marker rows missing: ", paste(setdiff(MARKERS, rownames(counts)), collapse = ", "))
}
if ("x" %in% methods::slotNames(counts)) {
  nonfinite_entries <- sum(!is.finite(methods::slot(counts, "x")))
} else if (is.matrix(counts)) {
  nonfinite_entries <- sum(!is.finite(counts))
} else {
  stop("Unrecognized count matrix representation")
}
if (nonfinite_entries > 0) stop("Nonfinite RNA matrix entries detected")
md$joined_nCount <- as.numeric(Matrix::colSums(counts))
md$positive_count <- is.finite(md$joined_nCount) & md$joined_nCount > 0
if ("nCount_RNA" %in% names(md)) {
  oldcount <- suppressWarnings(as.numeric(as.character(md$nCount_RNA)))
  md$count_mismatch <- is.finite(oldcount) & abs(oldcount - md$joined_nCount) > pmax(1, abs(oldcount) *
    0.001)
  if (any(!md$positive_count & is.finite(oldcount) & oldcount > 0)) {
    stop("Zero-count cells had positive original nCount; check RNA layers")
  }
} else {
  md$count_mismatch <- FALSE
}
qclist <- lapply(SOURCES, function(src) {
  d <- md[md$source_dataset == src, , drop = FALSE]
  data.frame(source_dataset = src, n_original = nrow(d), n_positive_count = sum(d$positive_count),
    n_excluded_from_recluster = sum(!d$positive_count), excluded_fraction = mean(!d$positive_count),
    n_count_mismatch = sum(d$count_mismatch), stringsAsFactors = FALSE)
})
qc <- bind_rows(qclist)
write_csv(qc, "00_104_RNA_input_QC")
if (any(qc$excluded_fraction > MAX_ZERO_COUNT_FRAC) || mean(md$count_mismatch) > MAX_COUNT_MISMATCH_FRAC) {
  stop("Excessive count anomalies; see 00_104_RNA_input_QC_v08.csv")
}

valid <- which(md$positive_count)
md$module_score <- NA_real_
md$module_high <- NA
expr <- as.matrix(counts[MARKERS, md$cell[valid], drop = FALSE])
expr <- log1p(sweep(expr, 2, md$joined_nCount[valid]/10000, "/"))
marker_audit <- list()
rank_mat <- matrix(NA_real_, nrow = length(MARKERS), ncol = length(valid), dimnames = list(MARKERS, md$cell[valid]))
valid_sources <- md$source_dataset[valid]
for (src in SOURCES) {
  ii <- which(valid_sources == src)
  for (j in seq_along(MARKERS)) {
    g <- MARKERS[j]
    vals <- expr[j, ii]
    n_detected <- sum(vals > 0)
    marker_audit[[length(marker_audit) + 1L]] <- data.frame(source_dataset = src, marker = g, n_valid_cells = length(ii),
      n_detected = n_detected, pct_detected = 100 * n_detected/length(ii), mean_logCP10K = mean(vals),
      status = if (n_detected == 0)
        "ALL_ZERO_IN_SOURCE" else "PRESENT_NONZERO", stringsAsFactors = FALSE)
    rank_mat[j, ii] <- rank(vals, ties.method = "average")/length(ii)
  }
}
marker_audit <- bind_rows(marker_audit)
write_csv(marker_audit, "01_104_common_marker_source_audit")
if (any(marker_audit$n_detected == 0)) {
  stop("At least one predeclared shared marker is completely absent in a source. ", "See 01_104_common_marker_source_audit_v08.csv; do not silently change panel.")
}
md$module_score[valid] <- colMeans(rank_mat)
rm(rank_mat, expr)
gc(verbose = FALSE)
for (src in SOURCES) {
  ii <- which(md$source_dataset == src & md$positive_count)

  score_rank <- rank(md$module_score[ii], ties.method = "average")/length(ii)
  md$module_high[ii] <- score_rank >= 0.8
}
md$marker_nondetect_limit <- !md$positive_count
write_csv(PROXY_THRESHOLDS, "02_104_proxy_threshold_specification")

base_idx <- split(which(md$positive_count & md$timepoint_main %in% MAIN_STAGES), interaction(md$source_dataset[md$positive_count &
  md$timepoint_main %in% MAIN_STAGES], md$timepoint_main[md$positive_count & md$timepoint_main %in%
  MAIN_STAGES], drop = TRUE))
base_summary <- bind_rows(lapply(base_idx, function(ii) {
  data.frame(source_dataset = md$source_dataset[ii[1]], timepoint_main = md$timepoint_main[ii[1]],
    n_expression_eligible = length(ii), n_frozen_focus = sum(md$frozen_focus[ii]), pct_frozen_focus = rate(md$frozen_focus[ii]),
    n_module_high = sum(md$module_high[ii]), pct_module_high = rate(md$module_high[ii]), mean_module_score = mean(md$module_score[ii]),
    note = "Fixed cells/expression; not independently re-clustered", stringsAsFactors = FALSE)
}))
write_csv(base_summary, "03_104_source_stage_fixed_score")

rm(counts)
gc(verbose = FALSE)

BRANCHES <- list(All_sources = SOURCES, Drop_17702045 = setdiff(SOURCES, "17702045"), Drop_GSE162610 = setdiff(SOURCES,
  "GSE162610"), Drop_GSE196928 = setdiff(SOURCES, "GSE196928"))
if (RUN_SINGLE_SOURCE) {
  for (src in SOURCES) BRANCHES[[paste0("Only_", src)]] <- src
}
cluster_audits <- list()
stage_audits <- list()
cell_assignments <- list()
arm_comp <- list()
statuses <- list()

cluster_audit <- function(d, branch, arm) {
  split_ids <- split(seq_len(nrow(d)), d$cluster)
  out <- lapply(names(split_ids), function(k) {
    ii <- split_ids[[k]]
    x <- data.frame(branch = branch, arm = arm, cluster = k, n_cells = length(ii), n_sources = length(unique(d$source_dataset[ii])),
      mean_module_score = mean(d$module_score[ii]), median_module_score = median(d$module_score[ii]),
      fraction_module_high = mean(d$module_high[ii]), pct_frozen_focus = rate(d$frozen_focus[ii]),
      stringsAsFactors = FALSE)
    for (r in seq_len(nrow(PROXY_THRESHOLDS))) {
      pr <- PROXY_THRESHOLDS[r, ]
      x[[paste0("selected_", pr$rule)]] <- length(ii) >= pr$min_cluster_cells && x$mean_module_score >=
        pr$min_mean_module_score && x$fraction_module_high >= pr$min_fraction_module_high
    }
    x
  })
  bind_rows(out)
}

apply_cluster_proxy <- function(d, table) {
  for (rule in PROXY_THRESHOLDS$rule) {
    col <- paste0("selected_", rule)
    chosen <- table$cluster[table[[col]]]
    d[[paste0("proxy_", rule)]] <- d$cluster %in% chosen
  }
  d
}

source_stage_summary <- function(d, branch, arm) {
  ix <- split(seq_len(nrow(d)), interaction(d$source_dataset, d$timepoint_main, drop = TRUE))
  bind_rows(lapply(ix, function(ii) {
    a <- data.frame(branch = branch, arm = arm, source_dataset = d$source_dataset[ii[1]], timepoint_main = d$timepoint_main[ii[1]],
      n_cells = length(ii), n_frozen_focus = sum(d$frozen_focus[ii]), pct_frozen_focus = rate(d$frozen_focus[ii]),
      n_module_high = sum(d$module_high[ii]), pct_module_high = rate(d$module_high[ii]), mean_module_score = mean(d$module_score[ii]),
      stringsAsFactors = FALSE)
    for (rule in PROXY_THRESHOLDS$rule) {
      k <- paste0("proxy_", rule)
      a[[paste0("n_", k)]] <- sum(d[[k]][ii])
      a[[paste0("pct_", k)]] <- rate(d[[k]][ii])
    }
    a
  }))
}

run_branch <- function(name, sources) {
  t0 <- Sys.time()
  logit("Starting ", name, " [", paste(sources, collapse = ";"), "]")
  use <- md$cell[md$source_dataset %in% sources & md$positive_count]
  if (length(use) < 200)
    stop("Insufficient cells for ", name)
  o <- subset(obj0, cells = use)
  DefaultAssay(o) <- "RNA"
  o@reductions <- list()
  o@graphs <- list()
  o@neighbors <- list()
  o <- JoinLayers(o, assay = "RNA")
  cm <- GetAssayData(o, assay = "RNA", layer = "counts")
  if (ncol(cm) != ncol(o) || !setequal(colnames(cm), Cells(o))) {
    stop("Incomplete joined counts for ", name)
  }
  rm(cm)
  o <- NormalizeData(o, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = FALSE)
  o <- FindVariableFeatures(o, assay = "RNA", selection.method = "vst", nfeatures = min(N_HVG, nrow(o)),
    verbose = FALSE)
  hvg <- VariableFeatures(o)
  if (length(hvg) < 100)
    stop("HVG count <100 in ", name)
  o <- ScaleData(o, assay = "RNA", features = hvg, verbose = FALSE)
  n_pc <- min(N_PCS_COMPUTE, length(hvg) - 1L, ncol(o) - 1L)
  use_pc <- min(N_PCS_USE, n_pc)
  if (use_pc < 5)
    stop("Too few PCs in ", name)
  o <- RunPCA(o, assay = "RNA", features = hvg, npcs = n_pc, reduction.name = "pca_104", reduction.key = "r104PC_",
    seed.use = SEED, verbose = FALSE)
  d <- md[Cells(o), , drop = FALSE]
  if (!identical(rownames(d), Cells(o)))
    stop("Cell order mismatch")

  o <- FindNeighbors(o, reduction = "pca_104", dims = seq_len(use_pc), graph.name = c("r104_pre_nn",
    "r104_pre_snn"), verbose = FALSE)
  o <- FindClusters(o, graph.name = "r104_pre_snn", resolution = RESOLUTION, cluster.name = "cluster_104_pre",
    random.seed = SEED, verbose = FALSE)
  pre <- as.character(o$cluster_104_pre)
  d$cluster <- pre
  pca_table <- cluster_audit(d, name, "PCA_only")
  dpca <- apply_cluster_proxy(d, pca_table)
  pca_stage <- source_stage_summary(dpca, name, "PCA_only")

  harmony_allowed <- length(sources) > 1L
  if (harmony_allowed) {
    o <- harmony::RunHarmony(object = o, group.by.vars = "source_dataset", reduction.use = "pca_104",
      dims.use = seq_len(use_pc), reduction.save = "harmony_104", max_iter = HARMONY_MAX_ITER,
      plot_convergence = FALSE, verbose = FALSE)
    if (!"harmony_104" %in% Reductions(o))
      stop("Harmony reduction missing")
    o <- FindNeighbors(o, reduction = "harmony_104", dims = seq_len(use_pc), graph.name = c("r104_harmony_nn",
      "r104_harmony_snn"), verbose = FALSE)
    o <- FindClusters(o, graph.name = "r104_harmony_snn", resolution = RESOLUTION, cluster.name = "cluster_104_harmony",
      random.seed = SEED, verbose = FALSE)
    har <- as.character(o$cluster_104_harmony)
    d$cluster <- har
    har_table <- cluster_audit(d, name, "Harmony")
    dhar <- apply_cluster_proxy(d, har_table)
    har_stage <- source_stage_summary(dhar, name, "Harmony")
  } else {
    har <- rep(NA_character_, nrow(d))
    har_table <- NULL
    dhar <- NULL
    har_stage <- NULL
  }

  celltab <- data.frame(branch = name, cell = d$cell, source_dataset = d$source_dataset, timepoint_main = d$timepoint_main,
    frozen_focus = d$frozen_focus, module_score = d$module_score, module_high = d$module_high, PCA_cluster = pre,
    Harmony_cluster = har, stringsAsFactors = FALSE)
  for (rule in PROXY_THRESHOLDS$rule) {
    key <- paste0("proxy_", rule)
    celltab[[paste0("PCA_", key)]] <- dpca[[key]]
    celltab[[paste0("Harmony_", key)]] <- if (harmony_allowed)
      dhar[[key]] else NA
  }

  comparisons <- list()
  if (harmony_allowed) {
    a <- data.frame(branch = name, arm_comparison = "PCA_only_vs_Harmony", n_cells = nrow(d), cluster_ARI = ari(pre,
      har), stringsAsFactors = FALSE)
    for (rule in PROXY_THRESHOLDS$rule) {
      key <- paste0("proxy_", rule)
      u <- dpca[[key]]
      v <- dhar[[key]]
      a[[paste0("Jaccard_", rule)]] <- jaccard(u, v)
      a[[paste0("agreement_", rule)]] <- mean(u == v)
      a[[paste0("n_PCA_selected_", rule)]] <- sum(u)
      a[[paste0("n_Harmony_selected_", rule)]] <- sum(v)
    }
    comparisons[[1]] <- a
  }

  cluster_audits[[paste0(name, "_PCA")]] <<- pca_table
  stage_audits[[paste0(name, "_PCA")]] <<- pca_stage
  if (harmony_allowed) {
    cluster_audits[[paste0(name, "_Harmony")]] <<- har_table
    stage_audits[[paste0(name, "_Harmony")]] <<- har_stage
    arm_comp[[name]] <<- comparisons[[1]]
  }
  cell_assignments[[name]] <<- celltab
  statuses[[name]] <<- data.frame(branch = name, status = "success", sources = paste(sources, collapse = ";"),
    n_cells = nrow(d), n_HVG = length(hvg), n_PC_used = use_pc, elapsed_minutes = round(as.numeric(difftime(Sys.time(),
      t0, units = "mins")), 2), comment = if (harmony_allowed)
      "PCA and Harmony re-clustered" else "PCA-only source control", stringsAsFactors = FALSE)
  logit("Done ", name, ": ", nrow(d), " cells, ", length(hvg), " HVGs")
  rm(o, d, dpca, dhar, celltab)
  invisible(gc(verbose = FALSE))
}

for (name in names(BRANCHES)) {
  tryCatch(run_branch(name, BRANCHES[[name]]), error = function(e) {
    logit("FAILED ", name, ": ", conditionMessage(e))
    statuses[[name]] <<- data.frame(branch = name, status = "failed", sources = paste(BRANCHES[[name]],
      collapse = ";"), n_cells = sum(md$source_dataset %in% BRANCHES[[name]] & md$positive_count),
      n_HVG = NA_integer_, n_PC_used = NA_integer_, elapsed_minutes = NA_real_, comment = conditionMessage(e),
      stringsAsFactors = FALSE)
    invisible(gc(verbose = FALSE))
  })
}

allcells <- bind_rows(cell_assignments)
write_csv(bind_rows(cluster_audits), "04_104_cluster_common_marker_audit")
write_csv(bind_rows(stage_audits), "05_104_recluster_source_stage_stability")
write_csv(bind_rows(arm_comp), "06_104_pre_vs_Harmony_agreement")
write_csv(bind_rows(statuses), "07_104_branch_status")

loo_comparisons <- list()
if ("All_sources" %in% names(cell_assignments)) {
  reference <- cell_assignments[["All_sources"]]
  for (name in setdiff(names(BRANCHES), "All_sources")) {
    if (!name %in% names(cell_assignments))
      next
    candidate <- cell_assignments[[name]]
    m <- match(candidate$cell, reference$cell)
    if (anyNA(m))
      stop("Cell ID not found in full reference for ", name)
    for (arm in c("PCA", "Harmony")) {
      col_cluster <- paste0(arm, "_cluster")
      if (all(is.na(candidate[[col_cluster]])))
        next
      z <- data.frame(branch = name, arm = arm, matched_cells = nrow(candidate), cluster_ARI_vs_all = ari(reference[[col_cluster]][m],
        candidate[[col_cluster]]), stringsAsFactors = FALSE)
      for (rule in PROXY_THRESHOLDS$rule) {
        key <- paste0(arm, "_proxy_", rule)
        x <- reference[[key]][m]
        y <- candidate[[key]]
        z[[paste0("Jaccard_", rule)]] <- jaccard(x, y)
        z[[paste0("agreement_", rule)]] <- mean(x == y)
        z[[paste0("n_full_selected_", rule)]] <- sum(x)
        z[[paste0("n_refit_selected_", rule)]] <- sum(y)
      }
      loo_comparisons[[length(loo_comparisons) + 1L]] <- z
    }
  }
}
write_csv(bind_rows(loo_comparisons), "08_104_leave_one_out_same_cell_agreement")

limits <- data.frame(contrast = c("14 dpi vs 28 dpi independent time effect", "14 dpi after omitting Figshare 17702045",
  "28 dpi after omitting GSE196928", "Animal-to-animal reproducibility for Figshare 17702045", "GSE196928 7 dpi vs 28 dpi between independent biological pools",
  "PLX5622 effect in retained atlas", "Sham vs uninjured as equivalent Control baseline", "A Clec7a-positive cell call in all three source datasets",
  "Fully independent state discovery from raw unselected cells"), identifiability = c("NOT IDENTIFIABLE: time and source completely confounded for these stages",
  "NOT ESTIMABLE: no 14 dpi coverage", "NOT ESTIMABLE: no 28 dpi coverage", "NOT TESTABLE: specimen IDs unavailable",
  "NOT TESTABLE: one pool per stage", "NOT ESTIMABLE: treatment group excluded", "NOT EQUALLY COMPARABLE",
  "NOT TESTABLE: original gene feature absent in two sources", "NOT PROVIDED: Step 18 cohort used"),
  stringsAsFactors = FALSE)
write_csv(limits, "09_104_nonidentifiable_comparisons")

if (EXPORT_CELL_ASSIGNMENTS && nrow(allcells) > 0) {
  fn <- file.path(OUTPUT, "10_104_cell_cluster_assignments_v08.csv.gz")
  con <- gzfile(fn, open = "wt", compression = 6)
  tryCatch({
    utils::write.table(allcells, file = con, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE,
      na = "", qmethod = "double")
  }, finally = close(con))
  logit("Saved compressed per-cell assignments: ", basename(fn))
}

plotfile <- file.path(OUTPUT, "11_104_state_stability_figures_v08.pdf")
pdf(plotfile, width = 11, height = 8, onefile = TRUE)
tryCatch({
  b <- base_summary
  p <- ggplot(b, aes(x = factor(timepoint_main, levels = MAIN_STAGES), y = pct_module_high, group = source_dataset,
    color = source_dataset)) + geom_point(size = 2) + geom_line(linewidth = 0.7) + labs(title = "Shared-marker module: source-relative top quintile by injury stage",
    subtitle = "Fixed gene panel; descriptive only; pooled/unknown specimen structure", x = "Stage (categorical)",
    y = "% module-high within source") + theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 35,
    hjust = 1))
  print(p)
  st <- bind_rows(stage_audits)
  if (nrow(st)) {
    z <- st[st$timepoint_main %in% MAIN_STAGES & st$branch %in% c("All_sources", "Drop_17702045",
      "Drop_GSE162610", "Drop_GSE196928"), , drop = FALSE]
    if (nrow(z)) {
      pp <- ggplot(z, aes(x = factor(timepoint_main, levels = MAIN_STAGES), y = pct_proxy_primary,
        group = arm, color = arm)) + geom_point(size = 1.4) + geom_line(linewidth = 0.6) + facet_grid(source_dataset ~
        branch, scales = "free_x") + labs(title = "Re-clustered shared-marker proxy: PCA vs Harmony / source exclusion",
        subtitle = "Exploratory cluster definition; no biological-replicate inference", x = "Injury stage",
        y = "% cells in proxy-selected clusters") + theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45,
        hjust = 1))
      print(pp)
    }
  }
}, error = function(e) logit("Plot error: ", conditionMessage(e)), finally = dev.off())

st <- bind_rows(statuses)
for (i in seq_len(nrow(st))) {
}
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUTPUT, "12_104_REPORT_v08.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUTPUT, "13_104_sessionInfo_v08.txt"))

if (MAKE_ZIP && requireNamespace("zip", quietly = TRUE)) {
  fs <- list.files(OUTPUT, pattern = "\\.(csv|csv.gz|pdf|txt|log)$", full.names = TRUE)
  zipfile <- file.path(OUTPUT, "104_state_stability_results_v08.zip")

  if (file.exists(zipfile))
    unlink(zipfile)
  oldwd <- getwd()
  tryCatch({
    setwd(OUTPUT)
    zip::zipr(zipfile = zipfile, files = basename(fs), recurse = FALSE)
  }, error = function(e) logit("ZIP creation failed: ", conditionMessage(e)), finally = setwd(oldwd))
  if (file.exists(zipfile)) {
    if (file.info(zipfile)$size > ZIP_LIMIT_MB * 1024^2) {
      unlink(zipfile)
      logit("ZIP exceeded ", ZIP_LIMIT_MB, " MB and was removed. CSVs remain.")
    } else {
      logit("ZIP OK: ", zipfile, " (", round(file.info(zipfile)$size/1024^2, 2), " MB)")
    }
  }
} else if (MAKE_ZIP) {
  logit("Optional package 'zip' is not installed. Install via install.packages('zip')")
}
logit("104 v08 finished. Upload the v08 ZIP and check 07_104_branch_status_v08.csv.")
