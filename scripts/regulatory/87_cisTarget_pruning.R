source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

python_exe <- PYSCENIC_PYTHON

num_workers <- 6L

rank_threshold <- 5000L
auc_threshold  <- 0.05
nes_threshold  <- 3.0
min_genes      <- 20L
chunk_size     <- 100L

mode_ctx <- "custom_multiprocessing"

mask_dropouts <- FALSE

all_modules <- FALSE

delete_stale_partial <- TRUE
overwrite_completed_output <- FALSE

step86_dir <- project_path("results/86_pySCENIC_GRNBoost2")

adjacency_file <- file.path(step86_dir, "09_86_pySCENIC_GRNBoost2_adjacencies_v03.tsv.gz")

loom_file <- file.path(step86_dir, "06_86_pySCENIC_GRN_discovery_input_v01.loom")

step86_qc_file <- file.path(step86_dir, "14_86_pySCENIC_FINAL_QC_v03.csv")

resource_dir <- reference_path("SCENIC/cisTarget_mm10_v10")

db_10kb <- file.path(resource_dir, "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

db_500bp <- file.path(resource_dir, "mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

motif_annotations <- file.path(resource_dir, "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")

outdir <- project_path("results/87_pySCENIC_cisTarget_pruning")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

helper_dir <- file.path(outdir, "_helpers_v01")

if (!dir.exists(helper_dir)) {
  dir.create(helper_dir, recursive = TRUE)
}

version_tag <- "v01"

resource_qc_file <- file.path(outdir, "01_87_pySCENIC_resource_QC_v01.csv")

environment_versions_file <- file.path(outdir, "02_87_pySCENIC_environment_versions_v01.csv")

preflight_qc_file <- file.path(outdir, "03_87_pySCENIC_cisTarget_preflight_QC_v01.csv")

command_file <- file.path(outdir, "04_87_pySCENIC_ctx_command_v01.txt")

log_file <- file.path(outdir, "05_87_pySCENIC_ctx_v01.log")

regulons_partial <- file.path(outdir, "06_87_pySCENIC_regulons_v01.partial.csv")

regulons_final <- file.path(outdir, "06_87_pySCENIC_regulons_v01.csv")

regulon_qc_file <- file.path(outdir, "07_87_pySCENIC_regulon_QC_v01.csv")

regulon_summary_file <- file.path(outdir, "08_87_pySCENIC_regulon_summary_v01.csv")

final_qc_file <- file.path(outdir, "09_87_pySCENIC_FINAL_QC_v01.csv")

report_file <- file.path(outdir, "10_87_pySCENIC_REPORT_v01.txt")

output_index_file <- file.path(outdir, "00_87_pySCENIC_OUTPUT_INDEX_v01.csv")

env_probe_py <- file.path(helper_dir, "87_environment_probe_v01.py")

preflight_py <- file.path(helper_dir, "87_cistarget_preflight_v01.py")

regulon_audit_py <- file.path(helper_dir, "87_regulon_audit_v01.py")

q <- shell_quote

make_qc <- qc_row

read_key_value <- function(path) {
    x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("metric", "value") %in% colnames(x))) {
        stop("Unexpected key/value file schema:\n", path)
    }
    setNames(as.character(x$value), x$metric)
}

stopifnot(file.copy("scripts/helpers/87_environment_probe.py", env_probe_py, overwrite = TRUE))

writeLines(c("import csv, os, sys", "from ctxcore.rnkdb import opendb", "", "db1 = sys.argv[1]", "db2 = sys.argv[2]", 
    "annot = sys.argv[3]", "out = sys.argv[4]", "", "rows = []", "", "for label, path in [('db_10kb', db1), ('db_500bp', db2)]:", 
    "    db = opendb(fname=path, name=label)", "    genes = list(db.genes)", "    rows.extend([", "        (f'{label}_opened', True),", 
    "        (f'{label}_n_genes', len(genes)),", "        (f'{label}_first_gene', genes[0] if genes else ''),", 
    "        (f'{label}_size_bytes', os.path.getsize(path)),", "    ])", "", "with open(annot, 'rt', encoding='utf-8', errors='replace') as f:", 
    "    n_annot = sum(1 for _ in f)", "", "rows.extend([", "    ('annotation_exists', os.path.exists(annot)),", 
    "    ('annotation_size_bytes', os.path.getsize(annot)),", "    ('annotation_lines', n_annot),", "])", 
    "", "with open(out, 'w', newline='', encoding='utf-8') as f:", "    w = csv.writer(f)", "    w.writerow(['metric','value'])", 
    "    w.writerows(rows)", "", "print('cisTarget resources opened successfully')"), preflight_py, useBytes = TRUE)

stopifnot(file.copy("scripts/helpers/87_regulon_check.py", regulon_audit_py, overwrite = TRUE))

required_files <- c(python_exe, adjacency_file, loom_file, step86_qc_file, db_10kb, db_500bp, motif_annotations)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
}

step86_qc <- read.csv(step86_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

if (any(step86_qc$status == "FAIL")) {
    stop("Step 86 QC contains FAIL item(s). ", "Do not proceed to Step 87.")
}

resource_qc <- data.frame(resource = c("mm10_10kb", "mm10_500bp", "motif_annotation"), file = c(basename(db_10kb), 
    basename(db_500bp), basename(motif_annotations)), size_bytes = c(file.info(db_10kb)$size, file.info(db_500bp)$size, 
    file.info(motif_annotations)$size), size_MB = c(file.info(db_10kb)$size, file.info(db_500bp)$size, 
    file.info(motif_annotations)$size)/1024^2, stringsAsFactors = FALSE)

resource_qc$status <- ifelse(resource_qc$size_bytes > 10^7, "PASS", "FAIL")

write.csv(resource_qc, resource_qc_file, row.names = FALSE)

if (any(resource_qc$status == "FAIL")) {
  stop("One or more cisTarget resource files are implausibly small.")
}

if (file.exists(environment_versions_file)) {
  unlink(environment_versions_file)
}

probe <- suppressWarnings(system2(command = python_exe, args = c(q(env_probe_py), q(environment_versions_file)), 
    stdout = TRUE, stderr = TRUE))

probe_status <- attr(probe, "status")
if (is.null(probe_status)) probe_status <- 0L

if (probe_status != 0 || !file.exists(environment_versions_file)) {
    stop("Python environment probe failed.")
}

env_versions <- read.csv(environment_versions_file, stringsAsFactors = FALSE, check.names = FALSE)

get_ver <- function(x) {
    z <- env_versions$version[env_versions$component == x]
    if (length(z) != 1) {
        stop("Could not uniquely resolve version for: ", x)
    }
    z
}

python_version <- get_ver("python")
pyscenic_version <- get_ver("pyscenic")
ctxcore_version <- get_ver("ctxcore")
pyarrow_version <- get_ver("pyarrow")

if (pyscenic_version != "0.12.1") {
    stop("pySCENIC version mismatch: observed ", pyscenic_version, ", required 0.12.1.")
}

if (file.exists(preflight_qc_file)) {
  unlink(preflight_qc_file)
}

preflight <- suppressWarnings(system2(command = python_exe, args = c(q(preflight_py), q(db_10kb), q(db_500bp), 
    q(motif_annotations), q(preflight_qc_file)), stdout = TRUE, stderr = TRUE))

preflight_status <- attr(preflight, "status")
if (is.null(preflight_status)) preflight_status <- 0L

if (preflight_status != 0 || !file.exists(preflight_qc_file)) {
    stop("cisTarget resource preflight failed. ", "Do NOT start ctx.")
}

preflight_qc <- read_key_value(preflight_qc_file)

if (as.numeric(preflight_qc[["db_10kb_n_genes"]]) <= 0 || as.numeric(preflight_qc[["db_500bp_n_genes"]]) <= 
    0 || as.numeric(preflight_qc[["annotation_lines"]]) <= 0) {
    stop("cisTarget preflight returned invalid resource content.")
}

if (file.exists(regulons_final) && !overwrite_completed_output) {
    stop("Completed Step 87 regulon output already exists:\n", regulons_final, "\nIt will not be overwritten automatically.")
}

if (file.exists(regulons_partial) && delete_stale_partial) {
    unlink(regulons_partial)
}

if (file.exists(regulons_partial)) {
  stop("A partial Step 87 output already exists.")
}

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

ctx_args <- c("-m", "pyscenic.cli.pyscenic", "ctx", q(adjacency_file), q(db_10kb), q(db_500bp), "--annotations_fname", 
    q(motif_annotations), "--expression_mtx_fname", q(loom_file), "--mode", mode_ctx, "--output", q(regulons_partial), 
    "--num_workers", as.character(num_workers), "--chunk_size", as.character(chunk_size), "--rank_threshold", 
    as.character(rank_threshold), "--auc_threshold", as.character(auc_threshold), "--nes_threshold", 
    as.character(nes_threshold), "--min_genes", as.character(min_genes), "--cell_id_attribute", "CellID", 
    "--gene_attribute", "Gene")

if (isTRUE(mask_dropouts)) {
    ctx_args <- c(ctx_args, "--mask_dropouts")
}

if (isTRUE(all_modules)) {
    ctx_args <- c(ctx_args, "--all_modules")
}

command_text <- paste(q(python_exe), paste(ctx_args, collapse = " "))

writeLines(c("============================================================", "Step 87 pySCENIC ctx exact command", 
    "============================================================", "", command_text, "", paste0("pySCENIC=", 
        pyscenic_version), paste0("ctxcore=", ctxcore_version), paste0("workers=", num_workers), paste0("mode=", 
        mode_ctx), paste0("chunk_size=", chunk_size), paste0("rank_threshold=", rank_threshold), paste0("auc_threshold=", 
        auc_threshold), paste0("nes_threshold=", nes_threshold), paste0("min_genes=", min_genes), paste0("mask_dropouts=", 
        mask_dropouts), paste0("all_modules=", all_modules), "", "Input adjacency:", adjacency_file, 
    "", "Expression loom:", loom_file, "", "Ranking DBs:", db_10kb, db_500bp, "", "Motif annotation:", 
    motif_annotations, "", "Temporary output:", regulons_partial, "", "Final output:", regulons_final, 
    "============================================================"), command_file)

if (file.exists(log_file)) {
  unlink(log_file)
}

start_time <- Sys.time()

ctx_status <- suppressWarnings(system2(command = python_exe, args = ctx_args, stdout = log_file, stderr = log_file))

if (is.null(ctx_status)) {
  ctx_status <- 0L
}

end_time <- Sys.time()

runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))

if (ctx_status != 0 || !file.exists(regulons_partial) || file.info(regulons_partial)$size <= 0) {
    stop(paste0("Step 87 ctx failed or was interrupted.\n", "Exit status: ", ctx_status, "\n", "Runtime: ", 
        round(runtime_minutes, 2), " min\n", "Log: ", log_file, "\n", "Official regulon output was NOT created."))
}

if (file.exists(regulon_qc_file)) {
  unlink(regulon_qc_file)
}

if (file.exists(regulon_summary_file)) {
  unlink(regulon_summary_file)
}

audit <- suppressWarnings(system2(command = python_exe, args = c(q(regulon_audit_py), q(regulons_partial), 
    q(regulon_qc_file), q(regulon_summary_file)), stdout = TRUE, stderr = TRUE))

audit_status <- attr(audit, "status")
if (is.null(audit_status)) audit_status <- 0L

if (audit_status != 0 || !file.exists(regulon_qc_file) || !file.exists(regulon_summary_file)) {
    stop("Step 87 regulon output failed post-run QC. ", "Partial file was NOT promoted.")
}

rq <- read_key_value(regulon_qc_file)

n_regulons <- as.numeric(rq[["n_regulons"]])

n_unique_regulon_names <- as.numeric(rq[["n_unique_regulon_names"]])

n_unique_targets <- as.numeric(rq[["n_unique_target_genes"]])

min_regulon_size <- as.numeric(rq[["min_regulon_size"]])

max_regulon_size <- as.numeric(rq[["max_regulon_size"]])

if (n_regulons <= 0 || n_unique_regulon_names <= 0 || n_unique_targets <= 0 || min_regulon_size <= 0) {
    stop("Step 87 regulon QC contains a blocking failure.")
}

if (file.exists(regulons_final) && overwrite_completed_output) {
    unlink(regulons_final)
}

rename_ok <- file.rename(regulons_partial, regulons_final)

if (!rename_ok) {
    stop("Regulon output passed QC but could not be renamed to final output.")
}

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("Step 86 QC", "PASS", "Step 86 has no blocking FAIL items.")

qc_list[[length(qc_list) + 1]] <- make_qc("pySCENIC exact version", if (pyscenic_version == "0.12.1") "PASS" else "FAIL", 
    pyscenic_version)

qc_list[[length(qc_list) + 1]] <- make_qc("10kb ranking DB", "PASS", paste0(round(file.info(db_10kb)$size/1024^2, 
    2), " MB"))

qc_list[[length(qc_list) + 1]] <- make_qc("500bp ranking DB", "PASS", paste0(round(file.info(db_500bp)$size/1024^2, 
    2), " MB"))

qc_list[[length(qc_list) + 1]] <- make_qc("Motif annotation", "PASS", paste0(round(file.info(motif_annotations)$size/1024^2, 
    2), " MB"))

qc_list[[length(qc_list) + 1]] <- make_qc("cisTarget DB preflight", "PASS", paste0("10kb genes=", preflight_qc[["db_10kb_n_genes"]], 
    "; 500bp genes=", preflight_qc[["db_500bp_n_genes"]]))

qc_list[[length(qc_list) + 1]] <- make_qc("ctx exit status", if (ctx_status == 0) "PASS" else "FAIL", 
    paste0("exit=", ctx_status))

qc_list[[length(qc_list) + 1]] <- make_qc("Regulons loaded", if (n_regulons > 0) "PASS" else "FAIL", 
    paste0("n_regulons=", n_regulons))

qc_list[[length(qc_list) + 1]] <- make_qc("Unique regulon names", if (n_unique_regulon_names == n_regulons) "PASS" else "WARN", 
    paste0(n_unique_regulon_names, " / ", n_regulons))

qc_list[[length(qc_list) + 1]] <- make_qc("Regulon target coverage", if (n_unique_targets > 0) "PASS" else "FAIL", 
    paste0("unique targets=", n_unique_targets))

qc_list[[length(qc_list) + 1]] <- make_qc("Regulon size range", if (min_regulon_size > 0) "PASS" else "FAIL", 
    paste0("min=", min_regulon_size, "; max=", max_regulon_size))

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

output_files <- c(resource_qc_file, environment_versions_file, preflight_qc_file, command_file, log_file, 
    regulons_final, regulon_qc_file, regulon_summary_file, final_qc_file, report_file)

output_index <- data.frame(file = basename(output_files), exists = file.exists(output_files), size_MB = ifelse(file.exists(output_files), 
    file.info(output_files)$size/1024^2, NA_real_), full_path = normalizePath(output_files, winslash = "/", 
    mustWork = FALSE), stringsAsFactors = FALSE)

write.csv(output_index, output_index_file, row.names = FALSE)

cat("WARN      = ", n_warn, "\n", sep = "")
cat("FAIL      = ", n_fail, "\n", sep = "")

if (n_fail > 0) {
    stop("Step 87 contains blocking FAIL item(s). ", "Do not proceed to Step 88.")
}
