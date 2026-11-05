format_number <- function (x) 
{
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

require_bioc_package <- function (pkg) 
{
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Required Bioconductor package unavailable: ", pkg)
    }
}

require_cran_package <- function (pkg) 
{
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Required CRAN package unavailable: ", pkg)
    }
}

require_package <- function (pkg) 
{
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Required package unavailable: ", pkg)
    }
}

log_separator <- function (char = "=", n = 90) 
{
    line <- paste(rep(char, n), collapse = "")
    if (exists("log_file", inherits = TRUE)) {
        write_log(line)
    }
    else {
        cat(line, "\n")
    }
}

write_log <- function (...) 
{
    x <- paste0(...)
    cat(x, "\n")
    cat(x, "\n", file = log_file, append = TRUE)
}

qc_row <- function (check, status, detail) 
{
    data.frame(check = check, status = status, detail = detail, stringsAsFactors = FALSE)
}

timestamp_message <- function (...) 
{
    cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

shell_quote <- function (x) 
{
    if (.Platform$OS.type == "windows") {
        shQuote(x, type = "cmd")
    }
    else {
        shQuote(x)
    }
}

read_csv_auto <- function (path) 
{
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        read.csv(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE)
    }
    else {
        read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
}

read_gzip_csv <- function (filename) 
{
    con <- gzfile(filename, open = "rt")
    on.exit(close(con), add = TRUE)
    read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
}

clean_gene_symbol <- function (x) 
{
    ifelse(grepl("^X[0-9]", x), sub("^X", "", x), x)
}

finite_median <- function (x) 
{
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    median(x)
}

write_result_csv <- function (x, filename) 
{
    utils::write.csv(x, file.path(OUT_DIR, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

write_gzip_csv <- function (x, filename) 
{
    con <- gzfile(filename, open = "wt")
    on.exit(close(con), add = TRUE)
    write.csv(x, file = con, row.names = FALSE, na = "")
}

write_gzip_csv_utf8 <- function (x, filename) 
{
    con <- gzfile(filename, open = "wt", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    write.table(x, file = con, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE, na = "")
}

write_csv_utf8 <- function (x, filename) 
{
    write.csv(x, file = filename, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
