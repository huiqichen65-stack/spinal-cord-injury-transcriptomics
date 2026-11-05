options(stringsAsFactors = FALSE, scipen = 50, width = 170)

ROOT <- "D:/zz"
TFS <- c("Mef2a", "Nr3c1", "Pura", "Ikzf1", "Bach1", "Fli1", "Zfp148", "Ep300", "Bclaf1", "Chd1")

B00_CP <- file.path(ROOT, "results", "135_TF_candidate_comparison_resampling_v02", "network_checkpoints_v02",
  "B00_seed9101_3000")

A01_CP <- file.path(ROOT, "results", "143_de_novo_gene_selection_network_sensitivity", "network_checkpoints_v01",
  "A01_whole_transcriptome_mean_bin_dispersion_3000")

LIPID_MANIFEST <- file.path(ROOT, "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_gene_set_manifest_v03.csv")

COMMON_UNIVERSE_REF <- file.path(ROOT, "results", "143_de_novo_gene_selection_network_sensitivity_v02",
  "08_143_v02_lipid_enrichment_common_universe_all.csv")

OUT <- file.path(ROOT, "results", "144_full_universe_lipid_enrichment_sensitivity")
CP_DIR <- file.path(OUT, "checkpoints_v01")
LOG <- file.path(OUT, "144_full_universe_lipid_enrichment_sensitivity_v01.log")

EXPECTED_GENES <- 3000L
MIN_SIZE <- 5L
MAX_SIZE <- 500L
FGSEA_EPS <- 0
FDR_CUTOFF <- 0.05
SEED_BASE <- 14400L

stopf <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(CP_DIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(LOG)) writeLines(character(0), LOG)

read_csv_any <- function(path) {
  if (!file.exists(path))
    stopf("Required file does not exist: ", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

check_cols <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) {
    stopf(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".incomplete")
  if (file.exists(tmp))
    unlink(tmp)
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (file.exists(path))
    unlink(path)
  if (!file.rename(tmp, path))
    stopf("Failed to write output: ", path)
  invisible(path)
}

write_out <- function(x, filename) {
  path <- file.path(OUT, filename)
  write_csv_atomic(x, path)
  invisible(path)
}

write_status <- function(status, detail = "") {
  x <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), status = status, detail = detail,
    stringsAsFactors = FALSE)
  write_out(x, "01_144_run_status_v01.csv")
}

md5_file <- function(path) {
  as.character(unname(tools::md5sum(path)))
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

ko_paths <- c(setNames(file.path(B00_CP, paste0("KO_", TFS, "_v02.csv.gz")), paste0("B00_", TFS)), setNames(file.path(A01_CP,
  paste0("KO_", TFS, "_v01.csv.gz")), paste0("A01_", TFS)))

input_paths <- c(lipid_manifest = LIPID_MANIFEST, B00_gene_list = file.path(B00_CP, "00_gene_list_v02.txt"),
  A01_gene_list = file.path(A01_CP, "00_gene_list_v01.txt"), ko_paths)

input_audit <- data.frame(item = names(input_paths), path = unname(input_paths), exists = file.exists(unname(input_paths)),
  size_bytes = as.numeric(file.info(unname(input_paths))$size), stringsAsFactors = FALSE)
input_audit <- rbind(input_audit, data.frame(item = "optional_143_v02_common_universe_reference", path = COMMON_UNIVERSE_REF,
  exists = file.exists(COMMON_UNIVERSE_REF), size_bytes = if (file.exists(COMMON_UNIVERSE_REF)) as.numeric(file.info(COMMON_UNIVERSE_REF)$size) else NA_real_,
  stringsAsFactors = FALSE))
write_out(input_audit, "00_144_INPUT_AUDIT_v01.csv")

required_audit <- input_audit[input_audit$item != "optional_143_v02_common_universe_reference", , drop = FALSE]
if (any(!required_audit$exists)) {
  missing_items <- required_audit$item[!required_audit$exists]
  stopf("Missing required inputs: ", paste(missing_items, collapse = ", "))
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
  stopf("Required package 'fgsea' is not installed. Install fgsea, then rerun.")
}

lipid <- read_csv_any(LIPID_MANIFEST)
check_cols(lipid, c("gene_set", "gene"), "Step 20 lipid manifest")

lipid <- lipid[!is.na(lipid$gene_set) & !is.na(lipid$gene) & nzchar(as.character(lipid$gene_set)) & nzchar(as.character(lipid$gene)),
  c("gene_set", "gene"), drop = FALSE]
lipid$gene_set <- as.character(lipid$gene_set)
lipid$gene <- as.character(lipid$gene)
lipid <- unique(lipid)

lipid_sets <- lapply(split(lipid$gene, lipid$gene_set), function(x) unique(as.character(x)))
if (length(lipid_sets) == 0L) stopf("No lipid gene sets in Step 20 manifest.")

B00_GENE_LIST <- readLines(file.path(B00_CP, "00_gene_list_v02.txt"), warn = FALSE, encoding = "UTF-8")
A01_GENE_LIST <- readLines(file.path(A01_CP, "00_gene_list_v01.txt"), warn = FALSE, encoding = "UTF-8")
if (length(B00_GENE_LIST) != EXPECTED_GENES || anyDuplicated(B00_GENE_LIST)) {
  stopf("B00 gene list is not exactly 3000 unique genes.")
}
if (length(A01_GENE_LIST) != EXPECTED_GENES || anyDuplicated(A01_GENE_LIST)) {
  stopf("A01 gene list is not exactly 3000 unique genes.")
}
EXPECTED_GENE_LISTS <- list(B00 = B00_GENE_LIST, A01 = A01_GENE_LIST)

manifest_summary <- do.call(rbind, lapply(names(lipid_sets), function(gs) {
  data.frame(pathway = gs, manifest_genes = length(lipid_sets[[gs]]), stringsAsFactors = FALSE)
}))
write_out(manifest_summary, "02_144_frozen_lipid_manifest_summary_v01.csv")

rank_data <- list()
universe_rows <- list()
coverage_rows <- list()
ko_audit_rows <- list()

for (mode in c("B00", "A01")) {
    for (i in seq_along(TFS)) {
        tf <- TFS[i]
        ko_path <- if (mode == "B00") {
            file.path(B00_CP, paste0("KO_", tf, "_v02.csv.gz"))
        }
        else {
            file.path(A01_CP, paste0("KO_", tf, "_v01.csv.gz"))
        }
        dat <- read_csv_any(ko_path)
        check_cols(dat, c("gene", "Z", "p.adj"), paste(mode, tf, "KO output"))
        dat$gene <- as.character(dat$gene)
        dat$Z <- safe_numeric(dat$Z)
        dat$p.adj <- safe_numeric(dat$p.adj)
        if (nrow(dat) != EXPECTED_GENES) {
            stopf(mode, "/", tf, " KO file has ", nrow(dat), " rows; expected ", EXPECTED_GENES, ".")
        }
        if (anyNA(dat$gene) || any(!nzchar(dat$gene)) || anyDuplicated(dat$gene)) {
            stopf(mode, "/", tf, " KO file has missing or duplicated gene names.")
        }
        if (!setequal(dat$gene, EXPECTED_GENE_LISTS[[mode]])) {
            stopf(mode, "/", tf, " KO gene membership does not match its 3000-gene list.")
        }
        if (!tf %in% dat$gene) {
            stopf(mode, "/", tf, " KO gene is absent from its KO output; ", "cannot verify the expected non-KO exclusion.")
        }
        nonko <- dat[dat$gene != tf & is.finite(dat$Z), , drop = FALSE]
        if (nrow(nonko) < EXPECTED_GENES - 1L) {
            stopf(mode, "/", tf, " has fewer than 2999 finite non-KO Z scores.")
        }
        z <- nonko$Z
        names(z) <- nonko$gene
        z <- sort(z, decreasing = TRUE)
        if (anyDuplicated(names(z))) {
            stopf(mode, "/", tf, " has duplicated gene names in ranked universe.")
        }
        key <- paste(mode, tf, sep = "::")
        rank_data[[key]] <- z
        universe_rows[[key]] <- data.frame(network_strategy = mode, TF = tf, ko_output_rows = nrow(dat), 
            ko_gene_excluded = tf, full_nonKO_ranked_genes = length(z), finite_Z_genes = sum(is.finite(z)), 
            min_Z = min(z), median_Z = stats::median(z), max_Z = max(z), stringsAsFactors = FALSE)
        ko_audit_rows[[key]] <- data.frame(network_strategy = mode, TF = tf, file = ko_path, file_md5 = md5_file(ko_path), 
            n_rows = nrow(dat), n_unique_genes = length(unique(dat$gene)), n_finite_Z = sum(is.finite(dat$Z)), 
            KO_gene_present = tf %in% dat$gene, n_nonKO_ranked = length(z), stringsAsFactors = FALSE)
        for (gs in names(lipid_sets)) {
            gs_genes <- lipid_sets[[gs]]
            overlap <- intersect(gs_genes, names(z))
            eligible <- length(overlap) >= MIN_SIZE && length(overlap) <= MAX_SIZE
            coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(network_strategy = mode, TF = tf, 
                pathway = gs, manifest_size = length(gs_genes), genes_in_full_nonKO_universe = length(overlap), 
                coverage_fraction = length(overlap)/length(gs_genes), eligible_for_fgsea = eligible, 
                eligibility_reason = if (eligible) 
                  "ELIGIBLE"
                else if (length(overlap) < MIN_SIZE) 
                  "BELOW_MIN_SIZE"
                else "ABOVE_MAX_SIZE", stringsAsFactors = FALSE)
        }
    }
}

universe_df <- do.call(rbind, universe_rows)
coverage_df <- do.call(rbind, coverage_rows)
ko_audit <- do.call(rbind, ko_audit_rows)

write_out(universe_df, "03_144_full_nonKO_gene_universe_summary_v01.csv")
write_out(coverage_df, "04_144_lipid_pathway_coverage_by_full_universe_v01.csv")
write_out(ko_audit, "05_144_KO_input_validation_v01.csv")
write_status("INPUTS_VALIDATED", paste0("Validated 20 KO files; ", length(lipid_sets), " lipid programs; full non-KO universes loaded."))

manifest_md5 <- md5_file(LIPID_MANIFEST)
fg_cp_dir <- file.path(CP_DIR, "fgsea_full_universe")
dir.create(fg_cp_dir, recursive = TRUE, showWarnings = FALSE)

empty_fg <- function() {
  data.frame(pathway = character(), pval = numeric(), padj = numeric(), log2err = numeric(), ES = numeric(),
    NES = numeric(), size = integer(), leading_edge = character(), TF = character(), network_strategy = character(),
    n_ranked_genes = integer(), analysis_signature = character(), rank_interpretation = character(),
    stringsAsFactors = FALSE)
}

fg_results <- list()
total_jobs <- length(TFS) * 2L
job_n <- 0L

for (mode in c("B00", "A01")) {
    for (i in seq_along(TFS)) {
        tf <- TFS[i]
        job_n <- job_n + 1L
        key <- paste(mode, tf, sep = "::")
        ko_path <- if (mode == "B00") {
            file.path(B00_CP, paste0("KO_", tf, "_v02.csv.gz"))
        }
        else {
            file.path(A01_CP, paste0("KO_", tf, "_v01.csv.gz"))
        }
        signature <- paste(mode, tf, md5_file(ko_path), manifest_md5, MIN_SIZE, MAX_SIZE, FGSEA_EPS, 
            SEED_BASE + i, length(rank_data[[key]]), sep = "|")
        ck_path <- file.path(fg_cp_dir, paste0("FGSEA_", mode, "_", tf, "_full3000_v01.csv"))
        sig_path <- paste0(ck_path, ".signature.txt")
        cached <- NULL
        if (file.exists(ck_path) && file.exists(sig_path)) {
            saved_signature <- tryCatch(readLines(sig_path, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
            if (length(saved_signature) == 1L && identical(saved_signature, signature)) {
                cached <- tryCatch(utils::read.csv(ck_path, check.names = FALSE, stringsAsFactors = FALSE), 
                  error = function(e) NULL)
                if (!is.null(cached) && nrow(cached) > 0L && (!"analysis_signature" %in% names(cached) || 
                  !all(cached$analysis_signature == signature))) {
                  cached <- NULL
                }
            }
        }
        if (!is.null(cached)) {
            fg_results[[key]] <- cached
            next
        }
        z <- rank_data[[key]]
        z <- sort(z[is.finite(z)], decreasing = TRUE)
        eligible_sets <- lapply(lipid_sets, function(gs) intersect(gs, names(z)))
        eligible_sets <- eligible_sets[vapply(eligible_sets, length, integer(1)) >= MIN_SIZE & vapply(eligible_sets, 
            length, integer(1)) <= MAX_SIZE]
        if (length(eligible_sets) == 0L) {
            out <- empty_fg()
        }
        else {
            set.seed(SEED_BASE + i)
            fg <- fgsea::fgseaMultilevel(pathways = eligible_sets, stats = z, minSize = MIN_SIZE, maxSize = MAX_SIZE, 
                eps = FGSEA_EPS)
            fg <- as.data.frame(fg)
            if (nrow(fg) > 0L) {
                fg$leading_edge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
                fg$leadingEdge <- NULL
                fg$TF <- tf
                fg$network_strategy <- mode
                fg$n_ranked_genes <- length(z)
                fg$analysis_signature <- signature
                fg$rank_interpretation <- "scTenifoldKnk distance-standardized network Z; NOT gene expression direction"
                out <- fg[, names(empty_fg()), drop = FALSE]
            }
            else {
                out <- empty_fg()
            }
        }
        write_csv_atomic(out, ck_path)
        writeLines(signature, sig_path, useBytes = TRUE)
        fg_results[[key]] <- out
        write_status("FGSEA_RUNNING", paste0("Completed ", job_n, "/", total_jobs, " TF x strategy jobs; last=", 
            mode, "/", tf))
    }
}

nonempty <- fg_results[vapply(fg_results, nrow, integer(1)) > 0L]
if (!length(nonempty)) {
  stopf("FGSEA produced no result rows. Review pathway coverage and log.")
}
fg_all <- do.call(rbind, nonempty)
fg_all$FDR_global_10TF <- NA_real_
for (mode in c("B00", "A01")) {
  idx <- which(fg_all$network_strategy == mode & is.finite(fg_all$pval))
  fg_all$FDR_global_10TF[idx] <- p.adjust(fg_all$pval[idx], method = "BH")
}
write_out(fg_all, "06_144_full_universe_lipid_enrichment_all_v01.csv")

b <- fg_all[fg_all$network_strategy == "B00", , drop = FALSE]
a <- fg_all[fg_all$network_strategy == "A01", , drop = FALSE]

select_cols <- c("TF", "pathway", "size", "NES", "pval", "FDR_global_10TF", "leading_edge", "n_ranked_genes")
b <- b[, select_cols, drop = FALSE]
a <- a[, select_cols, drop = FALSE]
names(b)[-(1:2)] <- paste0(names(b)[-(1:2)], "_B00")
names(a)[-(1:2)] <- paste0(names(a)[-(1:2)], "_A01")

paired <- merge(b, a, by = c("TF", "pathway"), all = TRUE, sort = TRUE)

paired$NES_same_sign <- with(paired, ifelse(is.finite(NES_B00) & is.finite(NES_A01), sign(NES_B00) ==
  sign(NES_A01), NA))
paired$NES_difference_A01_minus_B00 <- paired$NES_A01 - paired$NES_B00
paired$FDR05_B00 <- is.finite(paired$FDR_global_10TF_B00) & paired$FDR_global_10TF_B00 < FDR_CUTOFF
paired$FDR05_A01 <- is.finite(paired$FDR_global_10TF_A01) & paired$FDR_global_10TF_A01 < FDR_CUTOFF
paired$FDR05_both <- paired$FDR05_B00 & paired$FDR05_A01

split_leading <- function(x) {
  if (is.na(x) || !nzchar(x))
    character(0) else unique(strsplit(x, ";", fixed = TRUE)[[1]])
}
paired$leading_edge_intersection_n <- mapply(function(x, y) length(intersect(split_leading(x), split_leading(y))),
  paired$leading_edge_B00, paired$leading_edge_A01)
paired$leading_edge_union_n <- mapply(function(x, y) length(union(split_leading(x), split_leading(y))),
  paired$leading_edge_B00, paired$leading_edge_A01)
paired$leading_edge_jaccard <- ifelse(paired$leading_edge_union_n > 0, paired$leading_edge_intersection_n/paired$leading_edge_union_n,
  NA_real_)
paired$comparison_interpretation <- "Descriptive cross-strategy comparison; each FGSEA used its own full non-KO gene universe"
write_out(paired, "07_144_full_universe_lipid_enrichment_paired_comparison_v01.csv")

if (file.exists(COMMON_UNIVERSE_REF)) {
    common_ref <- read_csv_any(COMMON_UNIVERSE_REF)
    check_cols(common_ref, c("pathway", "TF", "case", "NES", "pval", "FDR_global_10TF"), "143 v02 common-universe FGSEA reference")
    common_ref$network_strategy <- as.character(common_ref$case)
    common_ref$common_universe_evaluable <- TRUE
    common_ref <- common_ref[, c("TF", "pathway", "network_strategy", "NES", "pval", "FDR_global_10TF", 
        "common_universe_evaluable"), drop = FALSE]
    names(common_ref)[4:6] <- c("NES_common842", "pval_common842", "FDR_common842")
    full_ref <- fg_all[, c("TF", "pathway", "network_strategy", "NES", "pval", "FDR_global_10TF"), drop = FALSE]
    names(full_ref)[4:6] <- c("NES_full_nonKO", "pval_full_nonKO", "FDR_full_nonKO")
    common_vs_full <- merge(full_ref, common_ref, by = c("TF", "pathway", "network_strategy"), all = TRUE, 
        sort = TRUE)
    common_vs_full$full_universe_evaluable <- is.finite(common_vs_full$pval_full_nonKO)
    common_vs_full$common_universe_evaluable[is.na(common_vs_full$common_universe_evaluable)] <- FALSE
    common_vs_full$NES_same_sign_common_vs_full <- with(common_vs_full, ifelse(is.finite(NES_common842) & 
        is.finite(NES_full_nonKO), sign(NES_common842) == sign(NES_full_nonKO), NA))
    common_vs_full$comparison_note <- "Common842 used shared non-KO intersection; full_nonKO uses each strategy's own full non-KO universe"
    write_out(common_vs_full, "08_144_common842_vs_full_universe_enrichment_comparison_v01.csv")
} else {
}

summary_rows <- list()
for (tf in TFS) {
  x <- paired[paired$TF == tf, , drop = FALSE]
  summary_rows[[tf]] <- data.frame(TF = tf, n_pathways_compared = nrow(x), n_pathways_significant_B00 = sum(x$FDR05_B00,
    na.rm = TRUE), n_pathways_significant_A01 = sum(x$FDR05_A01, na.rm = TRUE), n_pathways_significant_both = sum(x$FDR05_both,
    na.rm = TRUE), n_NES_same_sign = sum(x$NES_same_sign %in% TRUE, na.rm = TRUE), median_abs_NES_difference = if (any(is.finite(x$NES_difference_A01_minus_B00))) {
    median(abs(x$NES_difference_A01_minus_B00), na.rm = TRUE)
  } else NA_real_, median_leading_edge_jaccard = if (any(is.finite(x$leading_edge_jaccard))) {
    median(x$leading_edge_jaccard, na.rm = TRUE)
  } else NA_real_, stringsAsFactors = FALSE)
}
tf_summary <- do.call(rbind, summary_rows)
write_out(tf_summary, "09_144_TF_level_full_universe_lipid_summary_v01.csv")

coverage_counts <- aggregate(eligible_for_fgsea ~ network_strategy + pathway, data = coverage_df, FUN = sum)
n_jobs_done <- length(fg_results)
all_jobs_present <- n_jobs_done == total_jobs
qc <- data.frame(item = c("all_required_input_files_present", "fgsea_package_available", "frozen_step20_manifest_loaded",
  "all_11_lipid_programs_in_manifest", "all_20_KO_files_validated", "each_KO_gene_excluded_from_its_own_rank",
  "full_nonKO_universe_approximately_2999_per_test", "full_universe_pathway_coverage_generated", "all_20_TF_strategy_jobs_processed_or_checkpoint_reused",
  "full_universe_FGSEA_results_nonempty", "paired_B00_A01_comparison_generated", "network_reconstruction_performed",
  "virtual_KO_rerun", "expression_matrix_loaded", "RDS_read", "ZIP_input_required"), status = c("PASS",
  "PASS", "PASS", if (length(lipid_sets) == 11L) "PASS" else "CHECK", if (nrow(ko_audit) == 20L &&
    all(ko_audit$n_rows == EXPECTED_GENES)) "PASS" else "FAIL", if (all(universe_df$full_nonKO_ranked_genes <=
    EXPECTED_GENES - 1L)) "PASS" else "FAIL", if (all(universe_df$full_nonKO_ranked_genes >= EXPECTED_GENES -
    1L)) "PASS" else "CHECK", if (nrow(coverage_df) == 20L * length(lipid_sets)) "PASS" else "FAIL",
  if (all_jobs_present) "PASS" else "FAIL", if (nrow(fg_all) > 0L) "PASS" else "FAIL", if (nrow(paired) >
    0L) "PASS" else "FAIL", "NO", "NO", "NO", "NO", "NO"), detail = c("All B00/A01 KO CSV.GZ files, gene lists and Step20 manifest checked",
  as.character(utils::packageVersion("fgsea")), paste("Manifest MD5:", manifest_md5), paste("Observed lipid programs:",
    length(lipid_sets)), "Ten B00 and ten A01 KO files; required columns and unique genes validated",
  "Each TF itself removed from its KO rank before enrichment", paste("Observed range:", min(universe_df$full_nonKO_ranked_genes),
    "-", max(universe_df$full_nonKO_ranked_genes)), paste("Coverage rows:", nrow(coverage_df)), paste("Jobs:",
    n_jobs_done, "/", total_jobs), paste("FGSEA result rows:", nrow(fg_all)), paste("Paired TF-pathway rows:",
    nrow(paired)), "Existing saved networks only", "Existing saved KO CSV.GZ outputs only", "No expression matrix read",
  "No RDS read", "All inputs are original D:/zz project files"), stringsAsFactors = FALSE)
write_out(qc, "10_144_FINAL_QC_v01.csv")

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, "11_144_REPORT_v01.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "12_144_sessionInfo_v01.txt"), useBytes = TRUE)

out_files <- list.files(OUT, recursive = TRUE, full.names = TRUE)
out_files <- out_files[!dir.exists(out_files)]
out_index <- data.frame(relative_path = substring(out_files, nchar(OUT) + 2L), size_bytes = as.numeric(file.info(out_files)$size),
  stringsAsFactors = FALSE)
write_out(out_index, "13_144_OUTPUT_INDEX_v01.csv")

write_status("COMPLETE", paste0("Full-universe FGSEA completed for B00/A01; ", nrow(fg_all), " result rows."))
