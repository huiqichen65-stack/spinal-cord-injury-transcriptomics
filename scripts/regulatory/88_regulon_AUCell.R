source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

python_exe <- "C:/Users/Jiang Wu/anaconda3/envs/pyscenic_env/python.exe"

num_workers <- 6L
aucell_seed <- 8801L
auc_threshold <- 0.05

use_weights <- FALSE

reuse_existing_loom <- TRUE

delete_temp_mtx_after_loom <- TRUE

overwrite_completed_auc <- FALSE

source_rds <- paste0("D:/zz/results/20_AUCell_lipid_metabolism_scoring/", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

step84_qc_file <- paste0("D:/zz/results/84_pySCENIC_input_preparation/", "13_84_pySCENIC_FINAL_QC_v01.csv")

step85_dir <- "D:/zz/results/85_pySCENIC_GRNBoost2_input_preparation"

discovery_gene_file <- file.path(step85_dir, "07_85_pySCENIC_GRN_discovery_genes_v01.tsv.gz")

step85_qc_file <- file.path(step85_dir, "13_85_pySCENIC_FINAL_QC_v01.csv")

step87_dir <- "D:/zz/results/87_pySCENIC_cisTarget_pruning"

regulon_file <- file.path(step87_dir, "06_87_pySCENIC_regulons_v01.csv")

step87_qc_file <- file.path(step87_dir, "09_87_pySCENIC_FINAL_QC_v01.csv")

step87_report_file <- file.path(step87_dir, "10_87_pySCENIC_REPORT_v01.txt")

outdir <- "D:/zz/results/88_pySCENIC_regulon_AUCell_full_atlas"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

helper_dir <- file.path(outdir, "_helpers_v02")

if (!dir.exists(helper_dir)) {
  dir.create(helper_dir, recursive = TRUE)
}

input_summary_file <- file.path(outdir, "01_88_pySCENIC_input_summary_v02.csv")

environment_file <- file.path(outdir, "02_88_pySCENIC_environment_versions_v02.csv")

full_gene_file <- file.path(outdir, "03_88_pySCENIC_full_atlas_genes_v02.tsv.gz")

full_cell_file <- file.path(outdir, "04_88_pySCENIC_full_atlas_cells_v02.tsv.gz")

metadata_file <- file.path(outdir, "05_88_pySCENIC_full_atlas_metadata_v02.csv.gz")

full_loom_file <- file.path(outdir, "06_88_pySCENIC_full_atlas_9505genes_v02.loom")

loom_qc_file <- file.path(outdir, "07_88_pySCENIC_full_atlas_loom_QC_v02.csv")

regulon_preflight_file <- file.path(outdir, "08_88_pySCENIC_regulon_preflight_QC_v02.csv")

command_file <- file.path(outdir, "09_88_pySCENIC_AUCell_command_v02.txt")

log_file <- file.path(outdir, "10_88_pySCENIC_AUCell_v02.log")

auc_partial_file <- file.path(outdir, "11_88_pySCENIC_AUCell_full_atlas_v02.partial.csv")

auc_final_file <- file.path(outdir, "11_88_pySCENIC_AUCell_full_atlas_v02.csv.gz")

auc_qc_file <- file.path(outdir, "12_88_pySCENIC_AUCell_QC_v02.csv")

bundle_rds_file <- file.path(outdir, "13_88_pySCENIC_AUCell_full_atlas_bundle_v02.rds")

state_mean_file <- file.path(outdir, "14_88_pySCENIC_regulon_state_mean_AUC_v02.csv")

final_qc_file <- file.path(outdir, "15_88_pySCENIC_FINAL_QC_v02.csv")

report_file <- file.path(outdir, "16_88_pySCENIC_REPORT_v02.txt")

output_index_file <- file.path(outdir, "00_88_pySCENIC_OUTPUT_INDEX_v02.csv")

temp_mtx_file <- file.path(outdir, "_TEMP_88_full_atlas_9505genes_v02.mtx")

env_probe_py <- file.path(helper_dir, "88_environment_probe_v02.py")

loom_builder_py <- file.path(helper_dir, "88_build_full_atlas_loom_v02.py")

regulon_preflight_py <- file.path(helper_dir, "88_regulon_preflight_v02.py")

auc_audit_py <- file.path(helper_dir, "88_auc_audit_and_compress_v02.py")

expected_cells <- 58533L
expected_genes <- 9505L
expected_states <- 9L
expected_regulons <- 188L
expected_target_cells <- 17713L

rna_assay_name <- "RNA"
counts_layer_name <- "counts"
state_column <- "myeloid_state"
time_column <- "time"
umap_reduction <- "umap_myeloid_harmony30"

legacy_target_name <- "Clec7a+ lipid-associated Microglia"
formal_target_name <- "Clec7a+ lipid-handling Microglia"

q <- shell_quote

make_qc <- qc_row

read_key_value <- function(path) {
    x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("metric", "value") %in% colnames(x))) {
        stop("Unexpected key/value schema:\n", path)
    }
    setNames(as.character(x$value), x$metric)
}

to_bool <- function(x) {
  toupper(as.character(x)) == "TRUE"
}

write_one_col_gz <- function(x, path) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(as.character(x), con = con, sep = "\n", useBytes = TRUE)
}

read_discovery_genes <- function(
  path,
  expected_n
) {

  x1 <- tryCatch(
    read.delim(
      gzfile(path),
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )

  if (!is.null(x1) && ncol(x1) >= 1) {
    v1 <- trimws(
      as.character(
        x1[[1]]
      )
    )

    v1 <- v1[
      nzchar(v1)
    ]

    if (
      length(v1) == expected_n &&
      !anyDuplicated(v1)
    ) {
      return(v1)
    }
  }

  x0 <- read.delim(
    gzfile(path),
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  v0 <- trimws(
    as.character(
      x0[[1]]
    )
  )

  v0 <- v0[
    nzchar(v0)
  ]

  if (
    length(v0) != expected_n ||
    anyDuplicated(v0)
  ) {
    stop(
      "Could not resolve exactly ",
      expected_n,
      " unique discovery genes from:\n",
      path
    )
  }

  v0
}

read_auc_table <- function(path) {
    if (requireNamespace("data.table", quietly = TRUE)) {
        x <- data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = TRUE)
    }
    else {
        x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    if (ncol(x) < 2) {
        stop("AUCell CSV contains fewer than two columns.")
    }
    cell_ids <- as.character(x[[1]])
    auc <- as.matrix(x[, -1, drop = FALSE])
    storage.mode(auc) <- "double"
    rownames(auc) <- cell_ids
    auc
}

stopifnot(file.copy("scripts/helpers/88_environment_probe.py", env_probe_py, overwrite = TRUE))

stopifnot(file.copy("scripts/helpers/88_build_atlas_loom.py", loom_builder_py, overwrite = TRUE))

writeLines(c("import csv", "import sys", "from pyscenic.cli.utils import load_signatures", "", "reg_path = sys.argv[1]", 
    "out_path = sys.argv[2]", "", "sigs = list(load_signatures(reg_path))", "", "sizes = [len(list(s.genes)) for s in sigs]", 
    "names = [str(s.name) for s in sigs]", "", "if len(sigs) == 0:", "    raise RuntimeError('No regulons loaded.')", 
    "", "rows = [", "    ('n_regulons', len(sigs)),", "    ('n_unique_names', len(set(names))),", "    ('n_empty_regulons', sum(x == 0 for x in sizes)),", 
    "    ('min_regulon_size', min(sizes)),", "    ('max_regulon_size', max(sizes)),", "]", "", "with open(out_path, 'w', newline='', encoding='utf-8') as f:", 
    "    w = csv.writer(f)", "    w.writerow(['metric','value'])", "    w.writerows(rows)", "", "if any(x == 0 for x in sizes):", 
    "    raise RuntimeError('At least one empty regulon was loaded.')", "", "print(f'regulons={len(sigs)}')"), 
    regulon_preflight_py, useBytes = TRUE)

stopifnot(file.copy("scripts/helpers/88_auc_compress.py", auc_audit_py, overwrite = TRUE))

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("R package 'Matrix' is required.")
}

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("R package 'SeuratObject' is required.")
}

required_files <- c(python_exe, source_rds, step84_qc_file, discovery_gene_file, step85_qc_file, regulon_file, 
    step87_qc_file, step87_report_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

step84_qc <- read.csv(step84_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

step85_qc <- read.csv(step85_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

step87_qc <- read.csv(step87_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

if (any(step84_qc$status == "FAIL")) {
  stop("Step 84 contains blocking FAIL item(s).")
}

if (any(step85_qc$status == "FAIL")) {
  stop("Step 85 contains blocking FAIL item(s).")
}

if (any(step87_qc$status == "FAIL")) {
  stop("Step 87 contains blocking FAIL item(s).")
}

obj <- readRDS(source_rds)

if (!rna_assay_name %in% names(obj@assays)) {
  stop("RNA assay not found in source object.")
}

if (!state_column %in% colnames(obj@meta.data)) {
    stop("Metadata column not found: ", state_column)
}

if (!time_column %in% colnames(obj@meta.data)) {
    stop("Metadata column not found: ", time_column)
}

if (!umap_reduction %in% names(obj@reductions)) {
    stop("UMAP reduction not found: ", umap_reduction)
}

full_cells <- colnames(obj)

if (length(full_cells) != expected_cells || anyDuplicated(full_cells)) {
    stop("Full atlas cell count/uniqueness mismatch. Observed ", length(full_cells), ", expected ", expected_cells, 
        ".")
}

discovery_genes <- read_discovery_genes(discovery_gene_file, expected_n = expected_genes)

counts <- tryCatch(SeuratObject::LayerData(obj[[rna_assay_name]], layer = counts_layer_name), error = function(e) {
    stop("Could not extract RNA counts layer with LayerData():\n", conditionMessage(e))
})

if (is.null(rownames(counts)) || is.null(colnames(counts))) {
  stop("Counts matrix lacks gene/cell names.")
}

if (!identical(colnames(counts), full_cells)) {
    stop("RNA counts cell order does not match the Seurat object.")
}

missing_genes <- setdiff(discovery_genes, rownames(counts))

if (length(missing_genes) > 0) {
    stop("discovery genes are missing from full-atlas RNA counts.\n", "Number missing: ", length(missing_genes), 
        "\nFirst examples: ", paste(head(missing_genes, 20), collapse = ", "))
}

counts_sub <- counts[discovery_genes, full_cells, drop = FALSE]

if (nrow(counts_sub) != expected_genes || ncol(counts_sub) != expected_cells) {
    stop("Full-atlas SCENIC matrix dimension mismatch.")
}

if (inherits(counts_sub, "sparseMatrix")) {
  nz <- counts_sub@x
} else {
  nz <- as.numeric(counts_sub)
  nz <- nz[nz != 0]
}

if (length(nz) == 0) {
  stop("Full-atlas counts matrix contains no nonzero values.")
}

counts_nonnegative <- all(nz >= 0)

counts_integer_like <- all(abs(nz - round(nz)) < 1e-08)

if (!counts_nonnegative) {
  stop("RNA counts contain negative values.")
}

if (!counts_integer_like) {
  stop("RNA counts are not integer-like.")
}

md <- obj@meta.data[full_cells, , drop = FALSE]

md$cell <- rownames(md)

md$myeloid_state_source <- as.character(md[[state_column]])

md$myeloid_state_formal <- md$myeloid_state_source

md$myeloid_state_formal[md$myeloid_state_formal == legacy_target_name] <- formal_target_name

n_states <- length(unique(md$myeloid_state_formal))

n_target <- sum(md$myeloid_state_formal == formal_target_name, na.rm = TRUE)

if (n_states != expected_states) {
    stop("Expected ", expected_states, " myeloid states; observed ", n_states, ".")
}

if (n_target != expected_target_cells) {
    stop("Expected ", expected_target_cells, " target cells; observed ", n_target, ".")
}

umap <- obj@reductions[[umap_reduction]]@cell.embeddings

if (!setequal(rownames(umap), full_cells)) {
 stop("UMAP cell set does not match full atlas.")
}

umap <- umap[full_cells, , drop = FALSE]

write_one_col_gz(discovery_genes, full_gene_file)

write_one_col_gz(full_cells, full_cell_file)

metadata_export <- data.frame(cell = full_cells, myeloid_state_source = md$myeloid_state_source, myeloid_state = md$myeloid_state_formal, 
    time = as.character(md[[time_column]]), UMAP_1 = as.numeric(umap[, 1]), UMAP_2 = as.numeric(umap[, 
        2]), stringsAsFactors = FALSE)

metadata_con <- gzfile(metadata_file, "wt")

write.csv(metadata_export, metadata_con, row.names = FALSE)

close(metadata_con)

input_summary <- data.frame(metric = c("full_atlas_cells", "scenic_gene_universe", "myeloid_states", 
    "target_cells", "counts_nonnegative", "counts_integer_like", "source_RDS", "regulon_file"), value = c(expected_cells, 
    expected_genes, n_states, n_target, counts_nonnegative, counts_integer_like, source_rds, regulon_file), 
    stringsAsFactors = FALSE)

write.csv(input_summary, input_summary_file, row.names = FALSE)

if (file.exists(environment_file)) {
  unlink(environment_file)
}

probe <- suppressWarnings(system2(command = python_exe, args = c(q(env_probe_py), q(environment_file)), 
    stdout = TRUE, stderr = TRUE))

probe_status <- attr(probe, "status")
if (is.null(probe_status)) probe_status <- 0L

if (probe_status != 0 || !file.exists(environment_file)) {
    stop("Python environment probe failed.")
}

env_versions <- read.csv(environment_file, stringsAsFactors = FALSE, check.names = FALSE)

get_ver <- function(name) {
    x <- env_versions$version[env_versions$component == name]
    if (length(x) != 1) {
        stop("Could not uniquely resolve version: ", name)
    }
    x
}

python_version <- get_ver("python")
pyscenic_version <- get_ver("pyscenic")
loompy_version <- get_ver("loompy")
numpy_version <- get_ver("numpy")
pandas_version <- get_ver("pandas")

if (pyscenic_version != "0.12.1") {
    stop("pySCENIC version mismatch: observed ", pyscenic_version, ", required 0.12.1.")
}

if (file.exists(regulon_preflight_file)) {
  unlink(regulon_preflight_file)
}

reg_pf <- suppressWarnings(system2(command = python_exe, args = c(q(regulon_preflight_py), q(regulon_file), 
    q(regulon_preflight_file)), stdout = TRUE, stderr = TRUE))

reg_pf_status <- attr(reg_pf, "status")
if (is.null(reg_pf_status)) reg_pf_status <- 0L

if (reg_pf_status != 0 || !file.exists(regulon_preflight_file)) {
    stop("Step 87 regulon preflight failed.")
}

reg_pf_qc <- read_key_value(regulon_preflight_file)

n_regulons_input <- as.integer(reg_pf_qc[["n_regulons"]])

n_empty_regulons <- as.integer(reg_pf_qc[["n_empty_regulons"]])

if (n_regulons_input != expected_regulons || n_empty_regulons != 0) {
    stop("Regulon preflight mismatch: expected ", expected_regulons, " non-empty regulons; observed ", 
        n_regulons_input, ".")
}

reuse_loom_ok <- FALSE

if (isTRUE(reuse_existing_loom) && file.exists(full_loom_file) && file.exists(loom_qc_file)) {
    old_lq <- tryCatch(read_key_value(loom_qc_file), error = function(e) NULL)
    if (!is.null(old_lq)) {
        reuse_loom_ok <- (as.integer(old_lq[["n_genes"]]) == expected_genes && as.integer(old_lq[["n_cells"]]) == 
            expected_cells && to_bool(old_lq[["shape_match"]]) && to_bool(old_lq[["gene_order_match"]]) && 
            to_bool(old_lq[["cell_order_match"]]))
    }
}

if (reuse_loom_ok) {
} else {
    if (file.exists(temp_mtx_file)) {
        unlink(temp_mtx_file)
    }
    if (file.exists(full_loom_file)) {
        unlink(full_loom_file)
    }
    if (file.exists(loom_qc_file)) {
        unlink(loom_qc_file)
    }
    Matrix::writeMM(counts_sub, temp_mtx_file)
    if (!file.exists(temp_mtx_file) || file.info(temp_mtx_file)$size <= 0) {
        stop("Failed to create temporary Matrix Market file.")
    }
    loom_build_log <- file.path(outdir, "_TEMP_88_loom_build_v02.log")
    if (file.exists(loom_build_log)) {
        unlink(loom_build_log)
    }
    loom_status <- suppressWarnings(system2(command = python_exe, args = c(q(loom_builder_py), q(temp_mtx_file), 
        q(full_gene_file), q(full_cell_file), q(full_loom_file), q(loom_qc_file)), stdout = loom_build_log, 
        stderr = loom_build_log))
    if (is.null(loom_status)) {
        loom_status <- 0L
    }
    if (loom_status != 0 || !file.exists(full_loom_file) || !file.exists(loom_qc_file)) {
        stop("Full-atlas loom creation failed.\n", "Review:\n", loom_build_log)
    }
    lq <- read_key_value(loom_qc_file)
    if (as.integer(lq[["n_genes"]]) != expected_genes || as.integer(lq[["n_cells"]]) != expected_cells || 
        !to_bool(lq[["shape_match"]]) || !to_bool(lq[["gene_order_match"]]) || !to_bool(lq[["cell_order_match"]])) {
        stop("Full-atlas loom failed post-build QC.")
    }
    if (isTRUE(delete_temp_mtx_after_loom) && file.exists(temp_mtx_file)) {
        unlink(temp_mtx_file)
    }
}

rm(counts_sub, counts, nz)

gc()

if (file.exists(auc_final_file) && !overwrite_completed_auc) {
    stop("Completed Step 88 AUCell output already exists:\n", auc_final_file, "\nIt will not be overwritten automatically.")
}

if (file.exists(auc_partial_file)) {
    unlink(auc_partial_file)
}

if (file.exists(auc_final_file) && overwrite_completed_auc) {
    unlink(auc_final_file)
}

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

aucell_args <- c("-m", "pyscenic.cli.pyscenic", "aucell", q(full_loom_file), q(regulon_file), "--output", 
    q(auc_partial_file), "--num_workers", as.character(num_workers), "--seed", as.character(aucell_seed), 
    "--auc_threshold", as.character(auc_threshold), "--cell_id_attribute", "CellID", "--gene_attribute", 
    "Gene")

if (isTRUE(use_weights)) {
    aucell_args <- c(aucell_args, "--weights")
}

command_text <- paste(q(python_exe), paste(aucell_args, collapse = " "))

writeLines(c("============================================================", "Step 88 pySCENIC AUCell exact command", 
    "============================================================", "", command_text, "", paste0("pySCENIC=", 
        pyscenic_version), paste0("workers=", num_workers), paste0("seed=", aucell_seed), paste0("auc_threshold=", 
        auc_threshold), paste0("weights=", use_weights), paste0("cells=", expected_cells), paste0("gene_universe=", 
        expected_genes), paste0("regulons=", expected_regulons), "", "Expression loom:", full_loom_file, 
    "", "Regulons:", regulon_file, "", "Partial output:", auc_partial_file, "", "Final compressed output:", 
    auc_final_file, "============================================================"), command_file)

if (file.exists(log_file)) {
  unlink(log_file)
}

start_time <- Sys.time()

aucell_status <- suppressWarnings(system2(command = python_exe, args = aucell_args, stdout = log_file, 
    stderr = log_file))

if (is.null(aucell_status)) {
  aucell_status <- 0L
}

end_time <- Sys.time()

runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))

if (aucell_status != 0 || !file.exists(auc_partial_file) || file.info(auc_partial_file)$size <= 0) {
    stop(paste0("Step 88 AUCell failed or was interrupted.\n", "Exit status: ", aucell_status, "\n", 
        "Runtime: ", round(runtime_minutes, 2), " min\n", "Log: ", log_file, "\n", "No official final AUCell output was created."))
}

if (file.exists(auc_qc_file)) {
  unlink(auc_qc_file)
}

auc_audit <- suppressWarnings(system2(command = python_exe, args = c(q(auc_audit_py), q(auc_partial_file), 
    q(full_cell_file), q(auc_final_file), q(auc_qc_file)), stdout = TRUE, stderr = TRUE))

auc_audit_status <- attr(auc_audit, "status")

if (is.null(auc_audit_status)) {
  auc_audit_status <- 0L
}

if (auc_audit_status != 0 || !file.exists(auc_qc_file) || !file.exists(auc_final_file)) {
    stop("AUCell output failed post-run audit.")
}

aq <- read_key_value(auc_qc_file)

auc_n_cells <- as.integer(aq[["n_cells"]])

auc_n_regulons <- as.integer(aq[["n_regulons"]])

auc_nonfinite <- as.numeric(aq[["nonfinite_values"]])

auc_negative <- as.numeric(aq[["negative_values"]])

auc_zero_variance <- as.integer(aq[["zero_variance_regulons"]])

auc_all_zero <- as.integer(aq[["all_zero_regulons"]])

if (auc_n_cells != expected_cells || auc_n_regulons != expected_regulons || !to_bool(aq[["cell_order_match"]]) || 
    !to_bool(aq[["cell_set_match"]]) || auc_nonfinite != 0 || auc_negative != 0 || auc_all_zero != 0) {
    stop("AUCell post-run QC contains a blocking failure.")
}

auc_mat <- read_auc_table(auc_final_file)

if (nrow(auc_mat) != expected_cells || ncol(auc_mat) != expected_regulons) {
    stop("R-loaded AUCell matrix dimension mismatch.")
}

if (!identical(rownames(auc_mat), full_cells)) {
    stop("R-loaded AUCell cell order does not match full atlas.")
}

if (any(!is.finite(auc_mat)) || any(auc_mat < 0)) {
    stop("R-loaded AUCell matrix contains invalid values.")
}

state_vec <- md$myeloid_state_formal

state_levels <- sort(unique(state_vec))

state_mean_auc <- do.call(rbind, lapply(state_levels, function(st) {
    colMeans(auc_mat[state_vec == st, , drop = FALSE])
}))

rownames(state_mean_auc) <- state_levels

state_mean_export <- data.frame(myeloid_state = rownames(state_mean_auc), state_mean_auc, check.names = FALSE, 
    stringsAsFactors = FALSE)

write.csv(state_mean_export, state_mean_file, row.names = FALSE)

bundle <- list(auc = auc_mat, metadata = metadata_export, umap = umap, regulon_names = colnames(auc_mat), 
    state_mean_auc = state_mean_auc, parameters = list(pySCENIC = pyscenic_version, auc_threshold = auc_threshold, 
        seed = aucell_seed, num_workers = num_workers, use_weights = use_weights, full_cells = expected_cells, 
        gene_universe = expected_genes, regulons = expected_regulons, source_state_column = state_column, 
        formal_target_name = formal_target_name, source_target_name = legacy_target_name), source_files = list(source_rds = source_rds, 
        discovery_gene_file = discovery_gene_file, regulon_file = regulon_file, expression_loom = full_loom_file))

saveRDS(bundle, bundle_rds_file, compress = "gzip")

if (!file.exists(bundle_rds_file) || file.info(bundle_rds_file)$size <= 0) {
    stop("Failed to save Step 88 downstream RDS bundle.")
}

if (file.exists(auc_partial_file)) {
  unlink(auc_partial_file)
}

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("Step 84 QC", "PASS", "No blocking Step 84 FAIL items.")

qc_list[[length(qc_list) + 1]] <- make_qc("Step 85 QC", "PASS", "No blocking Step 85 FAIL items.")

qc_list[[length(qc_list) + 1]] <- make_qc("Step 87 QC", "PASS", "No blocking Step 87 FAIL items.")

qc_list[[length(qc_list) + 1]] <- make_qc("Full atlas cells", if (length(full_cells) == expected_cells) "PASS" else "FAIL", 
    paste0(length(full_cells), " / ", expected_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("Myeloid states", if (n_states == expected_states) "PASS" else "FAIL", 
    paste0(n_states, " / ", expected_states))

qc_list[[length(qc_list) + 1]] <- make_qc("Target-state cells", if (n_target == expected_target_cells) "PASS" else "FAIL", 
    paste0(n_target, " / ", expected_target_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("SCENIC gene universe", if (length(discovery_genes) == expected_genes) "PASS" else "FAIL", 
    paste0(length(discovery_genes), " / ", expected_genes))

qc_list[[length(qc_list) + 1]] <- make_qc("Counts non-negative", if (counts_nonnegative) "PASS" else "FAIL", 
    paste0("nonnegative=", counts_nonnegative))

qc_list[[length(qc_list) + 1]] <- make_qc("Counts integer-like", if (counts_integer_like) "PASS" else "FAIL", 
    paste0("integer_like=", counts_integer_like))

qc_list[[length(qc_list) + 1]] <- make_qc("pySCENIC exact version", if (pyscenic_version == "0.12.1") "PASS" else "FAIL", 
    pyscenic_version)

qc_list[[length(qc_list) + 1]] <- make_qc("Input regulons", if (n_regulons_input == expected_regulons) "PASS" else "FAIL", 
    paste0(n_regulons_input, " / ", expected_regulons))

qc_list[[length(qc_list) + 1]] <- make_qc("Full-atlas loom", "PASS", paste0(expected_genes, " genes x ", 
    expected_cells, " cells"))

qc_list[[length(qc_list) + 1]] <- make_qc("AUCell process", if (aucell_status == 0) "PASS" else "FAIL", 
    paste0("exit=", aucell_status))

qc_list[[length(qc_list) + 1]] <- make_qc("AUCell dimensions", if (auc_n_cells == expected_cells && auc_n_regulons == 
    expected_regulons) "PASS" else "FAIL", paste0(auc_n_cells, " cells x ", auc_n_regulons, " regulons"))

qc_list[[length(qc_list) + 1]] <- make_qc("AUCell cell order", if (to_bool(aq[["cell_order_match"]])) "PASS" else "FAIL", 
    paste0("match=", aq[["cell_order_match"]]))

qc_list[[length(qc_list) + 1]] <- make_qc("Finite AUC values", if (auc_nonfinite == 0) "PASS" else "FAIL", 
    paste0("nonfinite=", auc_nonfinite))

qc_list[[length(qc_list) + 1]] <- make_qc("Non-negative AUC values", if (auc_negative == 0) "PASS" else "FAIL", 
    paste0("negative=", auc_negative))

qc_list[[length(qc_list) + 1]] <- make_qc("All-zero regulons", if (auc_all_zero == 0) "PASS" else "FAIL", 
    paste0("all_zero=", auc_all_zero))

qc_list[[length(qc_list) + 1]] <- make_qc("Zero-variance regulons", if (auc_zero_variance == 0) "PASS" else "WARN", 
    paste0("zero_variance=", auc_zero_variance))

qc_list[[length(qc_list) + 1]] <- make_qc("Downstream RDS bundle", if (file.exists(bundle_rds_file) && 
    file.info(bundle_rds_file)$size > 0) "PASS" else "FAIL", paste0(round(file.info(bundle_rds_file)$size/1024^2, 
    2), " MB"))

qc_list[[length(qc_list) + 1]] <- make_qc("Runtime recorded", "PASS", paste0(round(runtime_minutes, 2), 
    " minutes"))

final_qc <- do.call(rbind, qc_list)

rownames(final_qc) <- NULL

write.csv(final_qc, final_qc_file, row.names = FALSE)

n_pass <- sum(final_qc$status == "PASS")

n_warn <- sum(final_qc$status == "WARN")

n_fail <- sum(final_qc$status == "FAIL")

overall_status <- if (n_fail > 0) {
  "FAIL"
} else if (n_warn > 0) {
  "PASS_WITH_WARNINGS"
} else {
  "PASS"
}

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_file)

output_files <- c(input_summary_file, environment_file, full_gene_file, full_cell_file, metadata_file, 
    full_loom_file, loom_qc_file, regulon_preflight_file, command_file, log_file, auc_final_file, auc_qc_file, 
    bundle_rds_file, state_mean_file, final_qc_file, report_file)

output_index <- data.frame(file = basename(output_files), exists = file.exists(output_files), size_MB = ifelse(file.exists(output_files), 
    file.info(output_files)$size/1024^2, NA_real_), full_path = normalizePath(output_files, winslash = "/", 
    mustWork = FALSE), stringsAsFactors = FALSE)

write.csv(output_index, output_index_file, row.names = FALSE)

cat("WARN           = ", n_warn, "\n", sep = "")
cat("FAIL           = ", n_fail, "\n", sep = "")

if (n_fail > 0) {
    stop("Step 88 contains blocking FAIL item(s). ", "Do not proceed to Step 89.")
}
