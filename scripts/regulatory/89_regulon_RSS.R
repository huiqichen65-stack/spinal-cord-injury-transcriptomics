source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

python_exe <- PYSCENIC_PYTHON

target_state <- "Clec7a+ lipid-handling Microglia"

top_n_export <- 30L
top_n_plot <- 20L
top_n_label <- 10L

plot_width_in <- 7.2
plot_height_in <- 5.2
png_dpi <- 600

overwrite_completed_output <- FALSE

step88_dir <- project_path("results/88_pySCENIC_regulon_AUCell_full_atlas")

auc_file <- file.path(step88_dir, "11_88_pySCENIC_AUCell_full_atlas_v02.csv.gz")

metadata_file <- file.path(step88_dir, "05_88_pySCENIC_full_atlas_metadata_v02.csv.gz")

state_mean_auc_file <- file.path(step88_dir, "14_88_pySCENIC_regulon_state_mean_AUC_v02.csv")

bundle_file <- file.path(step88_dir, "13_88_pySCENIC_AUCell_full_atlas_bundle_v02.rds")

step88_qc_file <- file.path(step88_dir, "15_88_pySCENIC_FINAL_QC_v02.csv")

step88_report_file <- file.path(step88_dir, "16_88_pySCENIC_REPORT_v02.txt")

outdir <- project_path("results/89_pySCENIC_RSS_target_regulon_ranking")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

helper_dir <- file.path(outdir, "_helpers_v01")

if (!dir.exists(helper_dir)) {
  dir.create(helper_dir, recursive = TRUE)
}

output_index_file <- file.path(outdir, "00_89_pySCENIC_OUTPUT_INDEX_v01.csv")

input_summary_file <- file.path(outdir, "01_89_pySCENIC_input_summary_v01.csv")

environment_file <- file.path(outdir, "02_89_pySCENIC_environment_versions_v01.csv")

rss_matrix_file <- file.path(outdir, "03_89_pySCENIC_RSS_all_states_v01.csv")

rss_long_file <- file.path(outdir, "04_89_pySCENIC_RSS_all_states_long_v01.csv")

target_ranking_file <- file.path(outdir, "05_89_Clec7a_lipid_handling_MG_regulon_RSS_ranking_v01.csv")

target_top_file <- file.path(outdir, "06_89_Clec7a_lipid_handling_MG_top_regulons_v01.csv")

plot_data_file <- file.path(outdir, "07_89_Figure6A_RSS_plot_data_v01.csv")

plot_pdf <- file.path(outdir, "08_89_Figure6A_target_RSS_ranking_v01.pdf")

plot_svg <- file.path(outdir, "08_89_Figure6A_target_RSS_ranking_v01.svg")

plot_png <- file.path(outdir, "08_89_Figure6A_target_RSS_ranking_v01.png")

rss_bundle_file <- file.path(outdir, "09_89_pySCENIC_RSS_bundle_v01.rds")

final_qc_file <- file.path(outdir, "10_89_pySCENIC_FINAL_QC_v01.csv")

report_file <- file.path(outdir, "11_89_pySCENIC_REPORT_v01.txt")

python_log_file <- file.path(outdir, "12_89_pySCENIC_RSS_python_v01.log")

env_probe_py <- file.path(helper_dir, "89_environment_probe_v01.py")

rss_py <- file.path(helper_dir, "89_calculate_RSS_v01.py")

expected_cells <- 58533L
expected_regulons <- 188L
expected_states <- 9L
expected_target_cells <- 17713L

q <- shell_quote

make_qc <- qc_row

rank_desc <- function(x) {
  rank(-x, ties.method = "min", na.last = "keep")
}

clean_regulon_tf <- function(x) {
  y <- sub("\\(\\+\\)$", "", x)
  y <- sub("\\(-\\)$", "", y)
  y <- trimws(y)
  y
}

stopifnot(file.copy("scripts/helpers/89_environment_probe.py", env_probe_py, overwrite = TRUE))

stopifnot(file.copy("scripts/helpers/89_calculate_RSS.py", rss_py, overwrite = TRUE))

required_files <- c(python_exe, auc_file, metadata_file, state_mean_auc_file, bundle_file, step88_qc_file,
  step88_report_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

step88_qc <- read.csv(step88_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

if (any(step88_qc$status == "FAIL")) {
  stop("Step 88 QC contains blocking FAIL item(s).")
}

bundle <- readRDS(bundle_file)

required_bundle_fields <- c("auc", "metadata", "umap", "regulon_names", "state_mean_auc")

if (!all(required_bundle_fields %in% names(bundle))) {
  stop("Step 88 bundle is missing required field(s): ", paste(setdiff(required_bundle_fields, names(bundle)),
    collapse = ", "))
}

auc_mat <- bundle$auc
metadata <- bundle$metadata
state_mean_auc <- bundle$state_mean_auc

if (nrow(auc_mat) != expected_cells || ncol(auc_mat) != expected_regulons) {
  stop("Unexpected AUCell matrix dimensions: ", nrow(auc_mat), " x ", ncol(auc_mat))
}

if (!"cell" %in% colnames(metadata)) {
  stop("Step 88 metadata lacks 'cell' column.")
}

if (!"myeloid_state" %in% colnames(metadata)) {
  stop("Step 88 metadata lacks formal 'myeloid_state' column.")
}

if (!identical(rownames(auc_mat), as.character(metadata$cell))) {
  stop("Step 88 AUCell and metadata cell order mismatch.")
}

states <- unique(as.character(metadata$myeloid_state))

if (length(states) != expected_states) {
  stop("Expected ", expected_states, " states; observed ", length(states), ".")
}

if (!target_state %in% states) {
  stop("Target state not present:\n", target_state)
}

n_target <- sum(metadata$myeloid_state == target_state)

if (n_target != expected_target_cells) {
  stop("Expected ", expected_target_cells, " target cells; observed ", n_target, ".")
}

if (file.exists(environment_file)) {
  unlink(environment_file)
}

probe <- suppressWarnings(system2(command = python_exe, args = c(q(env_probe_py), q(environment_file)),
  stdout = TRUE, stderr = TRUE))

probe_status <- attr(probe, "status")

if (is.null(probe_status)) {
  probe_status <- 0L
}

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
pandas_version <- get_ver("pandas")
numpy_version <- get_ver("numpy")
scipy_version <- get_ver("scipy")

if (pyscenic_version != "0.12.1") {
  stop("pySCENIC version mismatch: observed ", pyscenic_version, ", required 0.12.1.")
}

completed_outputs <- c(rss_matrix_file, target_ranking_file, rss_bundle_file, final_qc_file, report_file)

if (any(file.exists(completed_outputs)) && !overwrite_completed_output) {
  stop("One or more completed Step 89 outputs already exist.\n", "Set overwrite_completed_output=TRUE only if a deliberate rerun is intended.")
}

input_summary <- data.frame(metric = c("cells", "regulons", "states", "target_state", "target_cells",
  "pySCENIC", "AUC_file", "metadata_file", "bundle_file"), value = c(nrow(auc_mat), ncol(auc_mat),
  length(states), target_state, n_target, pyscenic_version, auc_file, metadata_file, bundle_file),
  stringsAsFactors = FALSE)

write.csv(input_summary, input_summary_file, row.names = FALSE)

if (file.exists(rss_matrix_file)) {
  unlink(rss_matrix_file)
}

if (file.exists(python_log_file)) {
  unlink(python_log_file)
}

rss_start <- Sys.time()

rss_status <- suppressWarnings(system2(command = python_exe, args = c(q(rss_py), q(auc_file), q(metadata_file),
  q("myeloid_state"), q(rss_matrix_file)), stdout = python_log_file, stderr = python_log_file))

if (is.null(rss_status)) {
  rss_status <- 0L
}

rss_end <- Sys.time()

runtime_seconds <- as.numeric(difftime(rss_end, rss_start, units = "secs"))

if (rss_status != 0 || !file.exists(rss_matrix_file) || file.info(rss_matrix_file)$size <= 0) {
  stop(paste0("RSS calculation failed.\n", "Exit status: ", rss_status, "\n", "Log: ", python_log_file))
}

rss_df <- read.csv(rss_matrix_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

rss_mat <- as.matrix(rss_df)

storage.mode(rss_mat) <- "double"

if (nrow(rss_mat) != expected_states || ncol(rss_mat) != expected_regulons) {
  stop("Unexpected RSS dimensions: ", nrow(rss_mat), " x ", ncol(rss_mat))
}

if (!target_state %in% rownames(rss_mat)) {
  stop("Target state absent from RSS matrix.")
}

if (any(!is.finite(rss_mat)) || any(rss_mat < 0) || any(rss_mat > 1)) {
  stop("RSS matrix contains invalid values.")
}

if (!setequal(colnames(rss_mat), colnames(auc_mat))) {
  stop("RSS regulon set differs from Step 88 AUCell regulon set.")
}

rss_mat <- rss_mat[, colnames(auc_mat), drop = FALSE]

rss_long <- do.call(rbind, lapply(rownames(rss_mat), function(st) {
  data.frame(myeloid_state = st, regulon = colnames(rss_mat), TF = clean_regulon_tf(colnames(rss_mat)),
    RSS = as.numeric(rss_mat[st, , drop = TRUE]), stringsAsFactors = FALSE)
}))

rownames(rss_long) <- NULL

write.csv(rss_long, rss_long_file, row.names = FALSE)

target_rss <- as.numeric(rss_mat[target_state, , drop = TRUE])

names(target_rss) <- colnames(rss_mat)

other_states <- setdiff(rownames(rss_mat), target_state)

max_other_rss <- apply(rss_mat[other_states, , drop = FALSE], 2, max)

best_other_state_rss <- apply(rss_mat[other_states, , drop = FALSE], 2, function(x) {
  names(which.max(x))[1]
})

target_rss_rank <- rank_desc(target_rss)

rss_margin <- target_rss - max_other_rss

rss_is_highest <- target_rss >= max_other_rss

state_mean_auc <- as.matrix(state_mean_auc)

storage.mode(state_mean_auc) <- "double"

if (!target_state %in% rownames(state_mean_auc)) {
  stop("Target state absent from Step 88 state-mean AUC matrix.")
}

if (!setequal(colnames(state_mean_auc), colnames(rss_mat))) {
  stop("State-mean AUC regulon set differs from RSS regulon set.")
}

state_mean_auc <- state_mean_auc[, colnames(rss_mat), drop = FALSE]

target_mean_auc <- as.numeric(state_mean_auc[target_state, , drop = TRUE])

names(target_mean_auc) <- colnames(state_mean_auc)

other_mean_auc <- state_mean_auc[setdiff(rownames(state_mean_auc), target_state), , drop = FALSE]

max_other_mean_auc <- apply(other_mean_auc, 2, max)

best_other_state_auc <- apply(other_mean_auc, 2, function(x) {
  names(which.max(x))[1]
})

mean_auc_delta <- target_mean_auc - max_other_mean_auc

target_auc_rank <- vapply(seq_len(ncol(state_mean_auc)), function(j) {
  rank_desc(state_mean_auc[, j])[which(rownames(state_mean_auc) == target_state)]
}, numeric(1))

target_auc_is_highest <- target_auc_rank == 1

target_auc_z <- vapply(seq_len(ncol(state_mean_auc)), function(j) {
    x <- state_mean_auc[, j]
    s <- stats::sd(x)
    if (!is.finite(s) || s == 0) {
        return(0)
    }
    (x[rownames(state_mean_auc) == target_state] - mean(x))/s
}, numeric(1))

ranking <- data.frame(regulon = colnames(rss_mat), TF = clean_regulon_tf(colnames(rss_mat)), target_state = target_state,
  target_RSS = target_rss, RSS_rank = as.integer(target_rss_rank), max_other_state_RSS = as.numeric(max_other_rss[colnames(rss_mat)]),
  RSS_margin_vs_best_other = as.numeric(rss_margin[colnames(rss_mat)]), target_has_highest_RSS = as.logical(rss_is_highest[colnames(rss_mat)]),
  best_other_state_by_RSS = as.character(best_other_state_rss[colnames(rss_mat)]), target_mean_AUC = target_mean_auc,
  max_other_state_mean_AUC = as.numeric(max_other_mean_auc[colnames(rss_mat)]), mean_AUC_delta_vs_best_other = as.numeric(mean_auc_delta[colnames(rss_mat)]),
  target_mean_AUC_rank_among_states = as.integer(target_auc_rank), target_has_highest_mean_AUC = as.logical(target_auc_is_highest),
  target_mean_AUC_z_across_states = as.numeric(target_auc_z), best_other_state_by_mean_AUC = as.character(best_other_state_auc[colnames(rss_mat)]),
  stringsAsFactors = FALSE)

ranking <- ranking[order(ranking$RSS_rank, -ranking$RSS_margin_vs_best_other, -ranking$target_mean_AUC),
  , drop = FALSE]

rownames(ranking) <- NULL

write.csv(ranking, target_ranking_file, row.names = FALSE)

top_regulons <- head(ranking, top_n_export)

write.csv(top_regulons, target_top_file, row.names = FALSE)

plot_data <- ranking

plot_data$plot_rank <- seq_len(nrow(plot_data))

plot_data$top_n_plot <- plot_data$plot_rank <= top_n_plot
plot_data$label_top <- plot_data$plot_rank <= top_n_label

write.csv(plot_data, plot_data_file, row.names = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("R package 'ggplot2' is required for the Step 89 RSS plot.")
}

library(ggplot2)

p <- ggplot(plot_data, aes(x = plot_rank, y = target_RSS)) + geom_point(aes(shape = top_n_plot), size = 2.2,
  alpha = 0.85) + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") + labs(x = "Regulon rank",
  y = "Regulon specificity score (RSS)", title = target_state) + theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold",
  size = 11, hjust = 0.5), axis.title = element_text(face = "bold"), plot.margin = margin(8, 12, 8,
  8))

label_df <- plot_data[plot_data$label_top, , drop = FALSE]

if (nrow(label_df) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = label_df, aes(label = TF), size = 3.4, box.padding = 0.35, 
        point.padding = 0.25, max.overlaps = Inf, min.segment.length = 0)
} else if (nrow(label_df) > 0) {
    p <- p + geom_text(data = label_df, aes(label = TF), size = 3.1, vjust = -0.7, check_overlap = TRUE)
}

ggsave(filename = plot_pdf, plot = p, width = plot_width_in, height = plot_height_in, units = "in", device = cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(filename = plot_svg, plot = p, width = plot_width_in, height = plot_height_in, units = "in", 
        device = svglite::svglite)
} else {
    svg(filename = plot_svg, width = plot_width_in, height = plot_height_in)
    print(p)
    dev.off()
}

ggsave(filename = plot_png, plot = p, width = plot_width_in, height = plot_height_in, units = "in", dpi = png_dpi)

rss_bundle <- list(rss = rss_mat, rss_long = rss_long, target_ranking = ranking, target_top = top_regulons,
  state_mean_auc = state_mean_auc, target_state = target_state, parameters = list(pySCENIC = pyscenic_version,
    RSS_method = "pyscenic.rss.regulon_specificity_scores", expected_cells = expected_cells, expected_regulons = expected_regulons,
    expected_states = expected_states, target_cells = n_target), source_files = list(auc_file = auc_file,
    metadata_file = metadata_file, step88_bundle = bundle_file))

saveRDS(rss_bundle, rss_bundle_file, compress = "gzip")

n_target_highest_rss <- sum(ranking$target_has_highest_RSS)

n_target_highest_auc <- sum(ranking$target_has_highest_mean_AUC)

top1 <- ranking[1, , drop = FALSE]

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("Step 88 QC", "PASS", "No blocking Step 88 FAIL items.")

qc_list[[length(qc_list) + 1]] <- make_qc("pySCENIC exact version", if (pyscenic_version == "0.12.1") "PASS" else "FAIL",
  pyscenic_version)

qc_list[[length(qc_list) + 1]] <- make_qc("AUCell input dimensions", if (nrow(auc_mat) == expected_cells &&
  ncol(auc_mat) == expected_regulons) "PASS" else "FAIL", paste0(nrow(auc_mat), " cells x ", ncol(auc_mat),
  " regulons"))

qc_list[[length(qc_list) + 1]] <- make_qc("Myeloid states", if (length(states) == expected_states) "PASS" else "FAIL",
  paste0(length(states), " / ", expected_states))

qc_list[[length(qc_list) + 1]] <- make_qc("Target-state cells", if (n_target == expected_target_cells) "PASS" else "FAIL",
  paste0(n_target, " / ", expected_target_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("RSS process", if (rss_status == 0) "PASS" else "FAIL", paste0("exit=",
  rss_status))

qc_list[[length(qc_list) + 1]] <- make_qc("RSS dimensions", if (nrow(rss_mat) == expected_states && ncol(rss_mat) ==
  expected_regulons) "PASS" else "FAIL", paste0(nrow(rss_mat), " states x ", ncol(rss_mat), " regulons"))

qc_list[[length(qc_list) + 1]] <- make_qc("RSS finite", if (all(is.finite(rss_mat))) "PASS" else "FAIL",
  paste0("nonfinite=", sum(!is.finite(rss_mat))))

qc_list[[length(qc_list) + 1]] <- make_qc("RSS range", if (min(rss_mat) >= 0 && max(rss_mat) <= 1) "PASS" else "FAIL",
  paste0(signif(min(rss_mat), 6), " to ", signif(max(rss_mat), 6)))

qc_list[[length(qc_list) + 1]] <- make_qc("RSS regulon set", if (setequal(colnames(rss_mat), colnames(auc_mat))) "PASS" else "FAIL",
  paste0(ncol(rss_mat), " / ", ncol(auc_mat)))

qc_list[[length(qc_list) + 1]] <- make_qc("Target ranking complete", if (nrow(ranking) == expected_regulons) "PASS" else "FAIL",
  paste0(nrow(ranking), " regulons ranked"))

qc_list[[length(qc_list) + 1]] <- make_qc("Target-highest RSS regulons", "PASS", paste0(n_target_highest_rss,
  " / ", expected_regulons))

qc_list[[length(qc_list) + 1]] <- make_qc("Target-highest mean-AUC regulons", "PASS", paste0(n_target_highest_auc,
  " / ", expected_regulons))

qc_list[[length(qc_list) + 1]] <- make_qc("Top-ranked regulon recorded", "PASS", paste0(top1$regulon,
  "; RSS=", signif(top1$target_RSS, 6)))

qc_list[[length(qc_list) + 1]] <- make_qc("RSS bundle", if (file.exists(rss_bundle_file) && file.info(rss_bundle_file)$size >
  0) "PASS" else "FAIL", paste0(round(file.info(rss_bundle_file)$size/1024^2, 3), " MB"))

qc_list[[length(qc_list) + 1]] <- make_qc("Figure 6A draft outputs", if (all(file.exists(c(plot_pdf,
  plot_svg, plot_png)))) "PASS" else "FAIL", "PDF + SVG + 600 dpi PNG")

qc_list[[length(qc_list) + 1]] <- make_qc("Runtime recorded", "PASS", paste0(round(runtime_seconds, 2),
  " seconds"))

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

output_files <- c(input_summary_file, environment_file, rss_matrix_file, rss_long_file, target_ranking_file,
  target_top_file, plot_data_file, plot_pdf, plot_svg, plot_png, rss_bundle_file, final_qc_file, report_file,
  python_log_file)

output_index <- data.frame(file = basename(output_files), exists = file.exists(output_files), size_MB = ifelse(file.exists(output_files),
  file.info(output_files)$size/1024^2, NA_real_), full_path = normalizePath(output_files, winslash = "/",
  mustWork = FALSE), stringsAsFactors = FALSE)

write.csv(output_index, output_index_file, row.names = FALSE)

cat("WARN             = ", n_warn, "\n", sep = "")
cat("FAIL             = ", n_fail, "\n", sep = "")

if (n_fail > 0) {
  stop("Step 89 contains blocking FAIL item(s). ", "Do not proceed to Step 90.")
}
