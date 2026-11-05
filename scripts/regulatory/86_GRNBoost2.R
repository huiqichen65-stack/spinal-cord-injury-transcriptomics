source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

python_exe <- PYSCENIC_PYTHON

num_workers <- 6L
grn_seed <- 8601L
use_sparse <- TRUE

delete_stale_partial <- TRUE
overwrite_completed_output <- FALSE

outdir <- project_path("results/86_pySCENIC_GRNBoost2")

loom_file <- file.path(outdir, "06_86_pySCENIC_GRN_discovery_input_v01.loom")

loom_qc_file <- file.path(outdir, "07_86_pySCENIC_loom_conversion_QC_v01.csv")

tf_file <- file.path(outdir, "05_86_pySCENIC_mouse_TFs_in_discovery_matrix_v01.txt")

step85_qc_file <- file.path(project_path("results/85_pySCENIC_GRNBoost2_input_preparation"), "13_85_pySCENIC_FINAL_QC_v01.csv")

version_tag <- "v03"

adj_partial <- file.path(outdir, paste0("09_86_pySCENIC_GRNBoost2_adjacencies_", version_tag, ".partial.tsv.gz"))

adj_final <- file.path(outdir, paste0("09_86_pySCENIC_GRNBoost2_adjacencies_", version_tag, ".tsv.gz"))

command_file <- file.path(outdir, paste0("10_86_pySCENIC_GRNBoost2_command_", version_tag, ".txt"))

log_file <- file.path(outdir, paste0("11_86_pySCENIC_GRNBoost2_", version_tag, ".log"))

adj_qc_file <- file.path(outdir, paste0("12_86_pySCENIC_adjacency_QC_", version_tag, ".csv"))

final_qc_file <- file.path(outdir, paste0("14_86_pySCENIC_FINAL_QC_", version_tag, ".csv"))

report_file <- file.path(outdir, paste0("15_86_pySCENIC_REPORT_", version_tag, ".txt"))

env_file <- file.path(outdir, paste0("02_86_pySCENIC_environment_versions_", version_tag, ".csv"))

helper_dir <- file.path(outdir, "_helpers_v03")

if (!dir.exists(helper_dir)) {
  dir.create(helper_dir, recursive = TRUE)
}

env_probe_py <- file.path(helper_dir, "86_environment_probe_v03.py")

adj_audit_py <- file.path(helper_dir, "86_adjacency_audit_v03.py")

q <- shell_quote

make_qc <- qc_row

read_qc_map <- function(path) {
    x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    setNames(as.character(x$value), x$metric)
}

to_bool <- function(x) {
  toupper(as.character(x)) == "TRUE"
}

stopifnot(file.copy("scripts/helpers/86_environment_probe.py", env_probe_py, overwrite = TRUE))

stopifnot(file.copy("scripts/helpers/86_adjacency_check.py", adj_audit_py, overwrite = TRUE))

required <- c(python_exe, loom_file, loom_qc_file, tf_file, step85_qc_file)

missing <- required[!file.exists(required)]

if (length(missing) > 0) {
  stop("Missing required file(s):\n", paste(missing, collapse = "\n"))
}

step85_qc <- read.csv(step85_qc_file, check.names = FALSE, stringsAsFactors = FALSE)

if (any(step85_qc$status == "FAIL")) {
  stop("Step 85 contains blocking FAIL item(s).")
}

loom_qc <- read_qc_map(loom_qc_file)

n_genes <- as.integer(loom_qc[["n_genes"]])
n_cells <- as.integer(loom_qc[["n_cells"]])

if (n_genes != 9505L || n_cells != 23285L || !to_bool(loom_qc[["shape_match"]]) || !to_bool(loom_qc[["gene_order_match"]]) ||
  !to_bool(loom_qc[["cell_order_match"]])) {
  stop("Existing v01 loom failed the compatibility check.")
}

tfs <- unique(trimws(readLines(tf_file, warn = FALSE, encoding = "UTF-8")))

tfs <- tfs[nzchar(tfs)]

if (file.exists(env_file)) {
  unlink(env_file)
}

probe <- suppressWarnings(system2(command = python_exe, args = c(q(env_probe_py), q(env_file)), stdout = TRUE,
  stderr = TRUE))

probe_status <- attr(probe, "status")
if (is.null(probe_status)) probe_status <- 0L

if (probe_status != 0 || !file.exists(env_file)) {
  stop("Python environment probe failed.")
}

env <- read.csv(env_file, check.names = FALSE, stringsAsFactors = FALSE)

ver <- function(name) {
  x <- env$version[env$component == name]
  if (length(x) != 1) {
    stop("Cannot uniquely resolve version: ", name)
  }
  x
}

python_version <- ver("python")
pyscenic_version <- ver("pyscenic")
arboreto_version <- ver("arboreto")
numpy_version <- ver("numpy")

if (pyscenic_version != "0.12.1") {
  stop("pySCENIC must be exactly 0.12.1; observed ", pyscenic_version)
}

if (file.exists(adj_final) && !overwrite_completed_output) {
  stop("Completed v02 adjacency already exists:\n", adj_final)
}

if (file.exists(adj_partial) && delete_stale_partial) {
    unlink(adj_partial)
}

if (file.exists(adj_partial)) {
  stop("A partial adjacency exists; refusing to overwrite it.")
}

thread_env <- c("OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "NUMEXPR_NUM_THREADS=1")

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

grn_args <- c("-m", "pyscenic.cli.arboreto_with_multiprocessing", q(loom_file), q(tf_file), "--method",
  "grnboost2", "--output", q(adj_partial), "--num_workers", as.character(num_workers), "--seed", as.character(grn_seed),
  "--cell_id_attribute", "CellID", "--gene_attribute", "Gene")

if (isTRUE(use_sparse)) {
  grn_args <- c(grn_args, "--sparse")
}

command_text <- paste(q(python_exe), paste(grn_args, collapse = " "))

writeLines(c("Step 86 v03 exact command", "", command_text, "", paste0("workers=", num_workers), paste0("seed=",
  grn_seed), paste0("sparse=", use_sparse), "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1",
  "NUMEXPR_NUM_THREADS=1"), command_file)

if (file.exists(log_file)) {
  unlink(log_file)
}

start_time <- Sys.time()

grn_status <- suppressWarnings(system2(command = python_exe, args = grn_args, stdout = log_file, stderr = log_file))

if (is.null(grn_status)) {
  grn_status <- 0L
}

end_time <- Sys.time()

runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))

if (grn_status != 0 || !file.exists(adj_partial) || file.info(adj_partial)$size <= 0) {
  stop(paste0("GRNBoost2 v02 failed or was interrupted.\n", "Exit status: ", grn_status, "\n", "Runtime: ",
    round(runtime_minutes, 2), " min\n", "Log: ", log_file, "\n", "No official final adjacency was created."))
}

if (file.exists(adj_qc_file)) {
  unlink(adj_qc_file)
}

audit <- suppressWarnings(system2(command = python_exe, args = c(q(adj_audit_py), q(adj_partial), q(adj_qc_file)),
  stdout = TRUE, stderr = TRUE))

audit_status <- attr(audit, "status")
if (is.null(audit_status)) audit_status <- 0L

if (audit_status != 0 || !file.exists(adj_qc_file)) {
  stop("Adjacency post-run QC failed. ", "Partial file was NOT promoted.")
}

aq <- read_qc_map(adj_qc_file)

n_edges <- as.numeric(aq[["n_edges"]])
n_unique_tfs <- as.numeric(aq[["n_unique_TFs"]])
n_unique_targets <- as.numeric(aq[["n_unique_targets"]])
n_bad <- as.numeric(aq[["nonfinite_importance"]])
n_negative <- as.numeric(aq[["negative_importance"]])

if (n_edges <= 0 || n_unique_tfs <= 0 || n_unique_targets <= 0 || n_bad > 0 || !to_bool(aq[["required_columns_present"]])) {
  stop("Adjacency contains blocking QC failure(s).")
}

if (file.exists(adj_final) && overwrite_completed_output) {
  unlink(adj_final)
}

if (!file.rename(adj_partial, adj_final)) {
  stop("Adjacency passed QC but could not be renamed to final output.")
}

final_qc <- do.call(rbind, list(make_qc("Step 85 blocking QC", "PASS", "No blocking Step 85 failures."),
  make_qc("loom", "PASS", paste0(n_genes, " genes x ", n_cells, " cells")), make_qc("pySCENIC version",
    if (pyscenic_version == "0.12.1") "PASS" else "FAIL", pyscenic_version), make_qc("Sparse loading",
    if (use_sparse) "PASS" else "WARN", paste0("--sparse=", use_sparse)), make_qc("GRNBoost2 process",
    if (grn_status == 0) "PASS" else "FAIL", paste0("exit=", grn_status)), make_qc("Adjacency edges",
    if (n_edges > 0) "PASS" else "FAIL", paste0("edges=", n_edges)), make_qc("Adjacency TFs", if (n_unique_tfs >
    0) "PASS" else "FAIL", paste0("TFs=", n_unique_tfs)), make_qc("Adjacency targets", if (n_unique_targets >
    0) "PASS" else "FAIL", paste0("targets=", n_unique_targets)), make_qc("Finite importance", if (n_bad ==
    0) "PASS" else "FAIL", paste0("nonfinite=", n_bad)), make_qc("Non-negative importance", if (n_negative ==
    0) "PASS" else "WARN", paste0("negative=", n_negative)), make_qc("Runtime", "PASS", paste0(round(runtime_minutes,
    2), " minutes"))))

write.csv(final_qc, final_qc_file, row.names = FALSE)

n_pass <- sum(final_qc$status == "PASS")
n_warn <- sum(final_qc$status == "WARN")
n_fail <- sum(final_qc$status == "FAIL")

overall <- if (n_fail > 0) {
  "FAIL"
} else if (n_warn > 0) {
  "PASS_WITH_WARNINGS"
} else {
  "PASS"
}

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_file)

cat("WARN    = ", n_warn, "\n", sep = "")
cat("FAIL    = ", n_fail, "\n", sep = "")

if (n_fail > 0) {
  stop("Step 86 v03 contains blocking FAIL item(s).")
}
