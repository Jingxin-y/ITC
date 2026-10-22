# Plot one saved scenario; reuse the original statistical/drawing functions.
if(.Platform$OS.type=="windows") Sys.setlocale("LC_CTYPE","English_United States.utf8")
args <- commandArgs(TRUE)
if(length(args)<3L || length(args)>4L)
  stop("Usage: Rscript plotting/plot_one_csv.R RR|RD INPUT.csv OUTPUT_DIR [N_EXPECTED=1000]")
scale <- match.arg(args[1],c("RR","RD"));input <- normalizePath(args[2],winslash="/",mustWork=TRUE)
out <- args[3];expected <- if(length(args)==4L)suppressWarnings(as.integer(args[4]))else 1000L
if(length(expected)!=1L || is.na(expected) || expected<1L)stop("N_EXPECTED must be a positive integer")
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
home <- dirname(normalizePath(script,winslash="/",mustWork=TRUE))
p <- new.env(parent=globalenv())
sys.source(file.path(home,scale,if(scale=="RR")"plot_cluster_results.R"else"plot_rd_cluster_results.R"),p)
before <- unname(tools::md5sum(input));d <- read.csv(input,stringsAsFactors=FALSE,check.names=FALSE)
required <- c("seed","method","ci_type","component","status","estimate","se","se_status")
required <- c(required,if(scale=="RR")c("seed_complete","low","up","population_log_RR_AC_in_BC")else
  c("lower","upper","interval_status","population_RD_AC_in_BC"))
if(!nrow(d)||!all(required%in%names(d)))stop("Unsupported saved-result columns; see plotting/README_zh.txt")
if(anyNA(d[c("seed","method","ci_type","component","status")]) ||
   anyDuplicated(d[c("seed","method","ci_type")]) || any(d$component!="marginal") ||
   any(d$status%in%c("pending","running")))stop("Missing, duplicate, non-marginal or unfinished records")
ids <- sort(unique(d$seed))
if(any(!is.finite(ids)|ids<1|ids!=floor(ids)|ids>expected))stop("Seed IDs must lie in 1:N_EXPECTED")
for(k in intersect(c("dgm","event","n_AC","n_BC","ess_ratio","ess_source","weight_mode",
                    "em_level","estimation_mode","model_weighting","ci_mode"),names(d)))
  if(anyNA(d[[k]])||length(unique(d[[k]]))!=1L)stop("Multiple or missing scenario settings in column: ",k)
if("dgm"%in%names(d)&&any(d$dgm!=paste0("webster2024_",tolower(scale),"_ess_em_v2")))
  stop("The input is not from the requested scale's v2 DGM")
keys <- paste(d$method,d$ci_type,sep="|")
expected_keys <- if(scale=="RR")c("brm|Wald","brm|Exact","brm-FC|Wald","brm-FC|Exact",
  "LB|Wald","LP|Wald","RLP|Wald","RLB|Wald","CMH|Wald","CMH_sandwich|Wald",
  "brmad|Wald","brmad|Exact","LB-FC|Wald","RLB-FC|Wald","LP-FC|Wald","RLP-FC|Wald",
  "LP-BFC|Wald","RLP-BFC|Wald") else c("brm|Wald","brm-BC|BC","brm-FC|Wald","brm-FC-BC|BC",
  "IB|Wald","RIB|Wald","RLPM|Wald","MN_eff|MN-ESS","RD_MAIC_sandwich|Wald")
if(!setequal(unique(keys),expected_keys)||any(table(d$seed)!=length(expected_keys)))
  stop("Each saved seed must retain the complete original method/CI plan, including failures")
if(scale=="RR") {
  if(!is.logical(d$seed_complete)||anyNA(d$seed_complete)||!all(d$seed_complete))stop("Incomplete saved seeds")
  d$record_complete <- d$seed_complete
  d <- d[!p$plot_excluded_method(d$method),,drop=FALSE]
  metrics <- p$plot_metrics(d,expected)
} else {
  if(anyNA(d$interval_status)||any(d$interval_status%in%c("pending","running","BC_not_run")))
    stop("RD intervals must have finished processing")
  d$method <- sub("-BC$","",d$method)
  point <- d$status%in%c("available","boundary_truncated")&is.finite(d$estimate)
  ci <- d$interval_status%in%c("available","boundary_truncated","BC_boundary_truncated")&
    is.finite(d$lower)&is.finite(d$upper)&d$lower<=d$upper
  d$status <- ifelse(point,ifelse(ci,ifelse(d$interval_status%in%c("boundary_truncated","BC_boundary_truncated"),
    "boundary_truncated","available"),"interval_failed"),"fit_failed")
  d$record_complete <- TRUE;d$low <- d$lower;d$up <- d$upper
  d$evaluation_truth_saved <- d$population_RD_AC_in_BC
  metrics <- p$rdplot_metrics(d,expected)
}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
title <- if(scale=="RR")"BC population marginal A/C log RR"else"BC population marginal A-C RD"
subtitle <- paste("Saved scenario:",basename(dirname(input)),"| completed seeds",length(ids),"/",expected)
if(scale=="RR")p$plot_save(metrics,out,title,subtitle,c("png","pdf"),d,se_max=4,bias_max=5)else
  p$rdplot_save(metrics,out,title,subtitle,c("png","pdf"),d,se_max=1,bias_max=1)
write.csv(metrics,file.path(out,"metrics.csv"),row.names=FALSE,na="")
stopifnot(identical(unname(tools::md5sum(input)),before))
cat("Plotted",scale,"from saved CSV;",length(ids),"completed seeds; no fitting or CI calculation.\n")
