source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, scipen = 50, width = 170)

ROOT <- PROJECT_ROOT
PREFIX <- "143_de_novo_gene_selection_network_sensitivity"
V01_OUT <- file.path(ROOT, "results", PREFIX)
OUT <- file.path(ROOT, "results", paste0(PREFIX, "_v02"))
POST_CP <- file.path(OUT, "postprocess_checkpoints_v02")
LOG <- file.path(OUT, paste0(PREFIX, "_v02.log"))

TFS <- c("Mef2a", "Nr3c1", "Pura", "Ikzf1", "Bach1", "Fli1", "Zfp148", "Ep300", "Bclaf1", "Chd1")
FDR_CUTOFF <- 0.05
FGSEA_EPS <- 0
MIN_SHARED_NONKO <- 500L
EXPECTED_GENES <- 3000L

A01_CP <- file.path(V01_OUT, "network_checkpoints_v01", "A01_whole_transcriptome_mean_bin_dispersion_3000")
B00_CP <- file.path(ROOT, "results", "135_TF_candidate_comparison_resampling_v02", "network_checkpoints_v02",
  "B00_seed9101_3000")

P <- list(selected_genes = file.path(V01_OUT, "03_143_selected_gene_list_v01.txt"), selection_audit = file.path(V01_OUT,
  "04_143_gene_selection_overlap_audit_v01.csv"), ko_summary = file.path(V01_OUT, "06_143_tenTF_KO_nonKO_summary_v01.csv"),
  ranking_stability = file.path(V01_OUT, "07_143_common_universe_ranking_stability_v01.csv"), lipid_manifest = file.path(ROOT,
    "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_gene_set_manifest_v03.csv"),
  b00_genes = file.path(B00_CP, "00_gene_list_v02.txt"))

stopf <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(POST_CP, recursive = TRUE, showWarnings = FALSE)

logline <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = LOG, append = TRUE)
  invisible(line)
}
readcsv <- function(path) {
  if (!file.exists(path))
    stopf("Required file missing: ", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, "rt")
    on.exit(close(con), add = TRUE)
    return(utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE))
  }
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
checkcols <- function(x, cols, who) {
  miss <- setdiff(cols, names(x))
  if (length(miss))
    stopf(who, " missing columns: ", paste(miss, collapse = ", "))
}
writecsv_atomic <- function(x, path) {
  tmp <- paste0(path, ".incomplete")
  if (file.exists(tmp))
    unlink(tmp)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(tmp, "wt")
    tryCatch(utils::write.csv(x, con, row.names = FALSE, na = ""), finally = close(con))
  } else {
    utils::write.csv(x, tmp, row.names = FALSE, na = "")
  }
  if (file.exists(path))
    unlink(path)
  if (!file.rename(tmp, path))
    stopf("Could not commit output: ", path)
}
writecsv <- function(x, filename) {
  path <- file.path(OUT, filename)
  writecsv_atomic(x, path)
  logline("Saved ", filename, " (", nrow(x), " rows)")
  invisible(path)
}
write_status <- function(state, details = "") {
  x <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), status = state, details = details,
    stringsAsFactors = FALSE)
  writecsv(x, "01_143_v02_postprocess_status.csv")
}
finite_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20L || length(unique(x[ok])) < 3L || length(unique(y[ok])) < 3L)
    return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}
safe_read_checkpoint <- function(path) {
  if (!file.exists(path))
    return(NULL)
  tryCatch(readcsv(path), error = function(e) {
    logline("Invalid checkpoint ignored: ", basename(path), "; ", conditionMessage(e))
    NULL
  })
}
flatten_fgsea <- function(fg, tf, mode) {
  fg <- as.data.frame(fg)
  if (!nrow(fg))
    return(data.frame())
  if ("leadingEdge" %in% names(fg)) {
    fg$leading_edge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    fg$leadingEdge <- NULL
  } else if (!"leading_edge" %in% names(fg)) {
    fg$leading_edge <- ""
  }
  fg$TF <- tf
  fg$case <- mode
  fg$rank_interpretation <- "scTenifoldKnk distance-standardized network Z; NOT gene expression direction"
  fg
}

if (!file.exists(LOG)) writeLines(character(0), LOG)
logline("Step 143 v02 continuation started; no network reconstruction; no RDS input.")

input_paths <- c(unlist(P, use.names = TRUE), setNames(file.path(A01_CP, paste0("KO_", TFS, "_v01.csv.gz")),
  paste0("A01_KO_", TFS)), setNames(file.path(B00_CP, paste0("KO_", TFS, "_v02.csv.gz")), paste0("B00_KO_",
  TFS)))
input_audit <- data.frame(item = names(input_paths), path = unname(input_paths), exists = file.exists(input_paths),
  size_bytes = as.numeric(file.info(input_paths)$size), stringsAsFactors = FALSE)
writecsv(input_audit, "00_143_v02_INPUT_AUDIT.csv")
if (any(!input_audit$exists)) {
  stopf("Missing required inputs: ", paste(input_audit$item[!input_audit$exists], collapse = ", "))
}
if (!requireNamespace("fgsea", quietly = TRUE)) stopf("Package fgsea is required for the lipid enrichment step.")

a01_genes <- readLines(P$selected_genes, warn = FALSE, encoding = "UTF-8")
b00_genes <- readLines(P$b00_genes, warn = FALSE, encoding = "UTF-8")
if (length(a01_genes) != EXPECTED_GENES || anyDuplicated(a01_genes)) stopf("A01 selected gene list is not exactly 3000 unique genes.")
if (length(b00_genes) != EXPECTED_GENES || anyDuplicated(b00_genes)) stopf("B00 selected gene list is not exactly 3000 unique genes.")

sel_audit <- readcsv(P$selection_audit)
checkcols(sel_audit, c("item", "value"), "v01 gene-selection audit")
if (length(intersect(a01_genes, b00_genes)) < MIN_SHARED_NONKO) stopf("A01 and B00 gene sets have unexpectedly small overlap.")

lipid <- readcsv(P$lipid_manifest)
checkcols(lipid, c("gene_set", "gene"), "Step20 lipid manifest")
lipid <- lipid[!is.na(lipid$gene_set) & !is.na(lipid$gene) & nzchar(as.character(lipid$gene_set)) & nzchar(as.character(lipid$gene)),
  , drop = FALSE]
lipid_sets <- lapply(split(as.character(lipid$gene), as.character(lipid$gene_set)), unique)
if (!length(lipid_sets)) stopf("No lipid gene sets found in Step20 manifest.")

ko_data <- list(B00 = list(), A01 = list())
for (tf in TFS) {
  p_b <- file.path(B00_CP, paste0("KO_", tf, "_v02.csv.gz"))
  p_a <- file.path(A01_CP, paste0("KO_", tf, "_v01.csv.gz"))
  b <- readcsv(p_b)
  a <- readcsv(p_a)
  checkcols(b, c("gene", "distance", "Z", "p.adj"), paste0("B00 ", tf))
  checkcols(a, c("gene", "distance", "Z", "p.adj"), paste0("A01 ", tf))
  if (nrow(b) != EXPECTED_GENES || nrow(a) != EXPECTED_GENES)
    stopf(tf, ": B00/A01 KO output does not contain exactly 3000 rows.")
  if (anyDuplicated(as.character(b$gene)) || anyDuplicated(as.character(a$gene)))
    stopf(tf, ": duplicate genes in KO output.")
  if (!setequal(as.character(b$gene), b00_genes) || !setequal(as.character(a$gene), a01_genes))
    stopf(tf, ": KO gene labels disagree with gene lists.")
  ko_data$B00[[tf]] <- b
  ko_data$A01[[tf]] <- a
  logline("Validated saved KO files: ", tf, " (B00 and A01)")
}
write_status("KO_FILES_VALIDATED", "All 20 saved KO files validated; no network or RDS recomputation.")

rank_vectors <- list()
rank_rows <- list()
for (tf in TFS) {
    b <- ko_data$B00[[tf]]
    a <- ko_data$A01[[tf]]
    b$gene <- as.character(b$gene)
    a$gene <- as.character(a$gene)
    b <- b[b$gene != tf & is.finite(suppressWarnings(as.numeric(b$Z))), , drop = FALSE]
    a <- a[a$gene != tf & is.finite(suppressWarnings(as.numeric(a$Z))), , drop = FALSE]
    common <- intersect(b$gene, a$gene)
    if (length(common) < MIN_SHARED_NONKO) 
        stopf(tf, ": fewer than ", MIN_SHARED_NONKO, " common evaluable non-KO genes.")
    bz <- as.numeric(b$Z[match(common, b$gene)])
    az <- as.numeric(a$Z[match(common, a$gene)])
    names(bz) <- common
    names(az) <- common
    bz <- sort(bz[is.finite(bz)], decreasing = TRUE)
    az <- sort(az[is.finite(az)], decreasing = TRUE)
    rank_vectors[[paste("B00", tf, sep = "::")]] <- bz
    rank_vectors[[paste("A01", tf, sep = "::")]] <- az
    btop <- names(head(bz, 50L))
    atop <- names(head(az, 50L))
    rank_rows[[tf]] <- data.frame(TF = tf, n_shared_nonKO = length(common), spearman_common_Z = finite_cor(bz[common], 
        az[common]), top50_overlap = length(intersect(btop, atop)), top50_jaccard = length(intersect(btop, 
        atop))/length(union(btop, atop)), stringsAsFactors = FALSE)
}
rank_recheck <- do.call(rbind, rank_rows)
saved_rank <- readcsv(P$ranking_stability)
checkcols(saved_rank, c("TF", "n_shared_nonKO", "spearman_common_Z", "topK_jaccard"), "v01 ranking stability")
if (!setequal(as.character(saved_rank$TF), TFS)) stopf("Saved v01 ranking stability table does not contain all ten TFs.")

for (tf in TFS) {
  old <- saved_rank[match(tf, saved_rank$TF), ]
  new <- rank_recheck[match(tf, rank_recheck$TF), ]
  if (is.finite(old$spearman_common_Z) && abs(as.numeric(old$spearman_common_Z) - as.numeric(new$spearman_common_Z)) >
    1e-08)
    stopf("Recreated common-universe Z correlation differs from v01 for ", tf)
}
writecsv(rank_recheck, "02_143_v02_common_universe_stability_recheck.csv")
write_status("COMMON_UNIVERSE_RANKS_RECREATED", "Common-gene Z vectors recreated from saved KO CSV.GZ files.")

fg_checkpoint_dir <- file.path(POST_CP, "fgsea_pair_checkpoints")
dir.create(fg_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
fg_results <- list()
total_pairs <- length(TFS) * 2L
pair_counter <- 0L

for (mode in c("B00", "A01")) {
    for (i in seq_along(TFS)) {
        tf <- TFS[i]
        pair_counter <- pair_counter + 1L
        ck <- file.path(fg_checkpoint_dir, paste0("FGSEA_", mode, "_", tf, "_v02.csv"))
        marker <- file.path(fg_checkpoint_dir, paste0("FGSEA_", mode, "_", tf, "_v02.status.txt"))
        cached <- safe_read_checkpoint(ck)
        if (!is.null(cached)) {
            fg_results[[paste(mode, tf, sep = "::")]] <- cached
            logline("FGSEA checkpoint reused [", pair_counter, "/", total_pairs, "]: ", mode, "/", tf)
            next
        }
        if (file.exists(marker) && grepl("NO_ELIGIBLE_PATHWAY", paste(readLines(marker, warn = FALSE), 
            collapse = " "))) {
            fg_results[[paste(mode, tf, sep = "::")]] <- data.frame()
            logline("FGSEA no-pathway marker reused [", pair_counter, "/", total_pairs, "]: ", mode, 
                "/", tf)
            next
        }
        z <- rank_vectors[[paste(mode, tf, sep = "::")]]
        z <- sort(z[is.finite(z)], decreasing = TRUE)
        eligible <- lapply(lipid_sets, function(gs) intersect(gs, names(z)))
        eligible <- eligible[vapply(eligible, length, integer(1)) >= 5L & vapply(eligible, length, integer(1)) <= 
            500L]
        logline("FGSEA START [", pair_counter, "/", total_pairs, "]: ", mode, "/", tf, "; pathways=", 
            length(eligible), "; ranked genes=", length(z), "; eps=", FGSEA_EPS)
        if (!length(eligible)) {
            writeLines("NO_ELIGIBLE_PATHWAY", marker, useBytes = TRUE)
            fg_results[[paste(mode, tf, sep = "::")]] <- data.frame()
            logline("FGSEA SKIPPED: no eligible lipid pathways for ", mode, "/", tf)
            next
        }
        set.seed(14300L + i)
        fg <- tryCatch(fgsea::fgseaMultilevel(pathways = eligible, stats = z, minSize = 5L, maxSize = 500L, 
            eps = FGSEA_EPS), error = function(e) {
            logline("FGSEA ERROR ", mode, "/", tf, ": ", conditionMessage(e))
            stop(e)
        })
        fg <- flatten_fgsea(fg, tf, mode)
        if (!nrow(fg)) {
            writeLines("NO_RESULT_ROWS", marker, useBytes = TRUE)
            fg_results[[paste(mode, tf, sep = "::")]] <- data.frame()
            logline("FGSEA returned zero rows: ", mode, "/", tf)
            next
        }
        writecsv_atomic(fg, ck)
        fg_results[[paste(mode, tf, sep = "::")]] <- fg
        logline("FGSEA DONE [", pair_counter, "/", total_pairs, "]: ", mode, "/", tf, "; pathways returned=", 
            nrow(fg))
        write_status("FGSEA_RUNNING", paste0("Completed ", pair_counter, "/", total_pairs, " TF x case pairs; last=", 
            mode, "/", tf))
    }
}

nonempty <- fg_results[vapply(fg_results, nrow, integer(1)) > 0L]
if (!length(nonempty)) stopf("No FGSEA result rows were produced or recovered.")
fg_all <- do.call(rbind, nonempty)
fg_all$FDR_global_10TF <- NA_real_
for (mode in c("B00", "A01")) {
  ii <- which(fg_all$case == mode & is.finite(fg_all$pval))
  fg_all$FDR_global_10TF[ii] <- p.adjust(fg_all$pval[ii], method = "BH")
}
writecsv(fg_all, "08_143_v02_lipid_enrichment_common_universe_all.csv")

key_cols <- c("TF", "pathway")
b <- fg_all[fg_all$case == "B00", , drop = FALSE]
a <- fg_all[fg_all$case == "A01", , drop = FALSE]
keep_b <- c(key_cols, "size", "NES", "pval", "FDR_global_10TF", "leading_edge")
keep_a <- keep_b
b <- b[, intersect(keep_b, names(b)), drop = FALSE]
a <- a[, intersect(keep_a, names(a)), drop = FALSE]
names(b)[-(1:2)] <- paste0(names(b)[-(1:2)], "_B00")
names(a)[-(1:2)] <- paste0(names(a)[-(1:2)], "_A01")
fg_compare <- merge(b, a, by = key_cols, all = TRUE, sort = TRUE)
if (all(c("NES_B00", "NES_A01") %in% names(fg_compare))) {
  fg_compare$NES_same_sign <- with(fg_compare, ifelse(is.finite(NES_B00) & is.finite(NES_A01), sign(NES_B00) ==
    sign(NES_A01), NA))
}
if (all(c("FDR_global_10TF_B00", "FDR_global_10TF_A01") %in% names(fg_compare))) {
  fg_compare$FDR05_both <- with(fg_compare, is.finite(FDR_global_10TF_B00) & is.finite(FDR_global_10TF_A01) &
    FDR_global_10TF_B00 < FDR_CUTOFF & FDR_global_10TF_A01 < FDR_CUTOFF)
}
writecsv(fg_compare, "09_143_v02_lipid_enrichment_paired_comparison.csv")

ko_complete <- all(vapply(TFS, function(tf) {
  file.exists(file.path(A01_CP, paste0("KO_", tf, "_v01.csv.gz"))) && file.exists(file.path(B00_CP,
    paste0("KO_", tf, "_v02.csv.gz")))
}, logical(1)))
qc <- data.frame(item = c("v01_gene_selection_outputs_present", "A01_gene_list_3000_unique", "B00_gene_list_3000_unique",
  "all_20_saved_KO_files_present", "all_20_saved_KO_files_validated", "ten_TF_common_universe_stability_rechecked",
  "all_TF_common_universe_at_least_500", "lipid_manifest_readable", "lipid_FGSEA_results_available",
  "paired_B00_A01_lipid_comparison_available", "network_reconstruction_repeated_in_v02", "expression_matrix_loaded_in_v02",
  "RDS_loaded_in_v02"), status = c("PASS", if (length(a01_genes) == EXPECTED_GENES) "PASS" else "FAIL",
  if (length(b00_genes) == EXPECTED_GENES) "PASS" else "FAIL", if (ko_complete) "PASS" else "FAIL",
  "PASS", "PASS", if (all(rank_recheck$n_shared_nonKO >= MIN_SHARED_NONKO)) "PASS" else "FAIL", "PASS",
  if (nrow(fg_all) > 0L) "PASS" else "FAIL", if (nrow(fg_compare) > 0L) "PASS" else "FAIL", "NO", "NO",
  "NO"), detail = c("v01 selected genes, audits, and ranking stability reused", paste(length(a01_genes),
  "A01 genes"), paste(length(b00_genes), "B00 genes"), "Ten A01 and ten B00 KO CSV.GZ files", "All KO files have 3000 unique genes and required columns",
  "Recomputed from saved KO files and checked against v01", paste("Minimum shared non-KO genes:", min(rank_recheck$n_shared_nonKO)),
  paste("Gene sets:", length(lipid_sets)), paste("FGSEA rows:", nrow(fg_all)), paste("Paired rows:",
    nrow(fg_compare)), "No network rebuild or virtual KO rerun", "No expression matrix loaded", "No RDS read"),
  stringsAsFactors = FALSE)
writecsv(qc, "10_143_v02_FINAL_QC.csv")

params <- data.frame(parameter = c("R", "source_version", "postprocess_version", "fgsea_method", "fgsea_eps",
  "minSize", "maxSize", "FDR_adjustment", "gene_universe", "network_reconstruction", "virtual_KO_recomputation",
  "expression_matrix_read", "RDS_read"), value = c(R.version.string, "143 v01 saved outputs", "143 v02 continuation",
  "fgseaMultilevel", as.character(FGSEA_EPS), "5", "500", "BH within each case across all TF x pathways",
  "paired common non-KO gene universe per TF", "NO", "NO", "NO", "NO"), stringsAsFactors = FALSE)
writecsv(params, "11_143_v02_environment_parameters.csv")

writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(OUT, "12_143_v02_REPORT.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "13_143_v02_sessionInfo.txt"), useBytes = TRUE)

files_out <- list.files(OUT, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
files_out <- files_out[!is.na(file.info(files_out)$isdir) & !file.info(files_out)$isdir]
index <- data.frame(relative_path = substring(files_out, nchar(OUT) + 2L), size_bytes = as.numeric(file.info(files_out)$size),
  stringsAsFactors = FALSE)
writecsv(index, "14_143_v02_OUTPUT_INDEX.csv")
write_status("COMPLETE", "Post-processing complete; A01 network and ten KOs reused from v01.")
logline("STEP 143 v02 CONTINUATION COMPLETE.")
