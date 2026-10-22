# Portable dispatcher for the unchanged saved-result plotting helpers.
if(.Platform$OS.type=="windows") Sys.setlocale("LC_CTYPE","English_United States.utf8")
args <- commandArgs(TRUE)
if(length(args)<2L) stop("Usage: Rscript plotting/plot_runs.R RR|RD RESULT_DIR [--out=DIR] [options]")
scale <- match.arg(args[1],c("RR","RD"))
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
home <- dirname(normalizePath(script,winslash="/",mustWork=TRUE))
helper <- file.path(home,scale,if(scale=="RR")"plot_cluster_results.R"else"plot_rd_cluster_results.R")
p <- new.env(parent=globalenv());sys.source(helper,p)
if(scale=="RR") p$plot_main(args[-1]) else {
  opt <- p$rdplot_options(args[-1])
  if(!is.null(opt)) p$plot_rd_cluster_results(opt$input,opt$out,opt$batch,
    opt$formats,opt$bias_max,opt$se_max,opt$dry_run)
}
