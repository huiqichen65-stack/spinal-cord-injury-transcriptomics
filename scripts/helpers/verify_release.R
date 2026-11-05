try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE)
scripts <- list.files("scripts", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
stopifnot(length(scripts) == 146L)
invisible(lapply(scripts, parse, keep.source = FALSE, encoding = "UTF-8"))
helpers <- new.env()
sys.source("scripts/helpers/common.R", helpers)
stopifnot(identical(environment(helpers$write_log), helpers),
          identical(environment(helpers$write_result_csv), helpers))
helpers$log_file <- tempfile()
helpers$write_log("helper environment check")
stopifnot(any(grepl("helper environment check", readLines(helpers$log_file))))
unlink(helpers$log_file)

repaired <- 'sensitivity/spatial/122_matched_control_spatial_sensitivity.R'
ex <- parse(file.path('scripts', repaired), keep.source = FALSE)
block <- Filter(function(e) is.call(e) && identical(e[[1]], as.name('tryCatch')), as.list(ex))[[1]][[2]]
env <- new.env(parent = baseenv())
env$IN <- list(frozen_scores='scores', rctd='rctd', lipid_manifest='genes',
              portions117='portions', rctd_features120='features')
env$SAMPLES <- paste0('area', 1:8); env$EXPECTED_SECTIONS <- 32L
env$TARGET <- 'target'
scores <- data.frame(GSM=paste0('GSM', 1:32), array_x=1:32, array_y=rep(1,32),
                     sample_label=rep(env$SAMPLES, each=4), timepoint='7dpi',
                     section_uid_internal=paste0('section',1:32),
                     score_Core_lipid_handling=seq_len(32)/32, Clec7a_RCTDnorm=seq_len(32)/64)
rctd <- scores; rctd$target <- scores$Clec7a_RCTDnorm
rctd$spatial_region_v03 <- 'region'; rctd$Microglia <- 0.2
fixtures <- list(scores=scores,rctd=rctd,
                 genes=data.frame(module='m',requested_gene='g',matched_gene='g',detected=TRUE),
                 portions=data.frame(section_uid_internal=scores$section_uid_internal,
                                     inferred_uninjured=seq_len(32)<=5),
                 features=data.frame(sample_label=env$SAMPLES,gene='g',
                                     RCTD_regression_fit=TRUE,RCTD_platform_normalization=FALSE))
env$read_csv <- function(p) fixtures[[p]]
env$step_log <- function(...) invisible(NULL)
env$must <- function(x,cols,label) stopifnot(all(cols %in% names(x)))
env$key <- function(df) paste(df$GSM,df$array_x,df$array_y,sep='||')
body <- as.list(block)[-1L]
start <- which(vapply(body, function(e) is.call(e) && identical(e[[1]], as.name('<-')) && identical(e[[2]], as.name('spot_scores')), logical(1)))
end <- which(vapply(body,function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
                     identical(e[[2]],as.name('obj')),logical(1)))[1]-1L
stopifnot(length(start)==1L, end>start)
for(i in seq.int(start,end)) eval(body[[i]],env)
stopifnot(identical(env$spot_scores$Clec7a_RCTDnorm,scores$Clec7a_RCTDnorm),
          identical(env$spot_scores$section_uid_internal,scores$section_uid_internal),
          sum(env$spot_scores$inferred_uninjured)==5L,
          all(env$spot_scores$total_myeloid_weight==0.2))
cat('PASS: script syntax, direct helper environments, and synthetic Step122 input alignment.\n')
