# Recompute displayed metrics from the saved Nibi CSVs; no models or CI inversion are run.
args <- commandArgs(TRUE)
if(length(args)!=1L) stop('Usage: Rscript code/plotting/plot_nibi_results.R REPOSITORY_ROOT')
root <- normalizePath(args[1],winslash='/',mustWork=TRUE)
if(.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','English_United States.utf8')
rr <- new.env(parent=globalenv()); rd <- new.env(parent=globalenv())
sys.source(file.path(root,'code/plotting/RR/plot_cluster_results.R'),rr)
sys.source(file.path(root,'code/plotting/RD/plot_rd_cluster_results.R'),rd)
raw <- file.path(root,'original result'); out <- file.path(root,'results')
config <- read.csv(file.path(raw,'config/SCENARIOS.csv'),stringsAsFactors=FALSE)
full <- list(); index <- list(); checks <- list()
close_number <- function(x,y,tol=1e-8) {
 if(is.na(x)&&is.na(y)) return(TRUE)
 is.finite(x)&&is.finite(y)&&abs(x-y)<=tol*(1+abs(y))
}
audit <- read.csv(file.path(out,'validation/independent_RR_RD_metrics.csv'),stringsAsFactors=FALSE)
audit_key <- with(audit,paste(scale,analysis,scenario_id,method,ci_type,sep='|'))
stopifnot(!anyDuplicated(audit_key))
for(i in seq_len(nrow(config))) {
 cfg <- config[i,]; d <- read.csv(file.path(raw,cfg$csv_path),stringsAsFactors=FALSE,check.names=FALSE)
 stopifnot(length(unique(d$seed))==1000L,nrow(d)==if(cfg$scale=='RR')18000L else 9000L)
 if(cfg$scale=='RR') {
  stopifnot(all(d$seed_complete));d$record_complete <- d$seed_complete
  m <- rr$plot_metrics(d,1000L)
  m$source_method <- m$method
 } else {
  d$method <- sub('-BC$','',d$method)
  point <- d$status%in%c('available','boundary_truncated')&is.finite(d$estimate)
  ci <- d$interval_status%in%c('available','boundary_truncated','BC_boundary_truncated') &
   is.finite(d$lower)&is.finite(d$upper)&d$lower<=d$upper
  boundary <- d$interval_status%in%c('boundary_truncated','BC_boundary_truncated')
  d$status <- ifelse(point,ifelse(ci,ifelse(boundary,'boundary_truncated','available'),'interval_failed'),'fit_failed')
  d$record_complete <- TRUE;d$low <- d$lower;d$up <- d$upper
  d$evaluation_truth_saved <- d$population_RD_AC_in_BC
  m <- rd$rdplot_metrics(d,1000L)
  m$source_method <- ifelse(m$ci_type=='BC',paste0(m$method,'-BC'),m$method)
 }
 m$scale <- cfg$scale;m$analysis <- cfg$analysis;m$scenario_id <- cfg$scenario_id
 m$cell <- cfg$cell;m$event <- cfg$event;m$n <- cfg$n;m$ess_ratio <- cfg$ess_ratio
 m$input_md5 <- cfg$md5
 for(j in seq_len(nrow(m))) {
  key <- paste(cfg$scale,cfg$analysis,cfg$scenario_id,m$source_method[j],m$ci_type[j],sep='|')
  a <- audit[match(key,audit_key),];stopifnot(nrow(a)==1L,!is.na(a$scale))
  map <- c(n_point='finite_point',n_ci='valid_intervals',bias='bias',rmse='rmse',
   empirical_sd='empirical_sd',mean_se='mean_se',empirical_sd_for_accuracy='sd_on_SE_subset',
   accuracy='mean_SE_over_SD',n_se='finite_SE_and_point',coverage='coverage',mean_width='mean_width',
   returned_and_covered_fraction='returned_and_covered',n_boundary='boundary_count')
  for(field in names(map)) {
   if(!close_number(m[[field]][j],a[[map[[field]]]]))
    stop('Metric differs from independent audit: ',key,' ',field,': ',m[[field]][j],' vs ',a[[map[[field]]]])
  }
 }
 folder <- file.path(out,cfg$scale,cfg$analysis,cfg$cell);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 title <- if(cfg$scale=='RR') 'BC population marginal A/C log RR' else 'BC population marginal A-C RD'
 subtitle <- paste(cfg$analysis,'|',cfg$cell,'| MC1000')
 if(cfg$scale=='RR') {
  keep <- !rr$plot_excluded_method(m$method)
  rr$plot_save(m[keep,,drop=FALSE],folder,title,subtitle,c('png','pdf'),d,se_max=4,bias_max=5)
 } else rd$rdplot_save(m,folder,title,subtitle,c('png','pdf'),d,se_max=1,bias_max=1)
 write.csv(m,file.path(folder,'metrics.csv'),row.names=FALSE,na='')
 writeLines(c('Input: saved per-replicate CSV; all 1000 seeds and failed rows retained.',
  paste('Input MD5:',cfg$md5),'Coverage = covered / valid intervals with truth.',
  'Interval return rate and returned-and-covered fraction use planned MC1000.',
  'Accuracy uses mean valid model SE / empirical SD on the same valid-SE subset.',
  'RR graphs omit brmad, LP-BFC and RLP-BFC for continuity with the existing plotting style; full tables retain all methods.'),
  file.path(folder,'analysis_notes.txt'))
 full[[i]] <- m
 index[[i]] <- data.frame(scale=cfg$scale,analysis=cfg$analysis,scenario_id=cfg$scenario_id,
  cell=cfg$cell,event=cfg$event,n=cfg$n,ess_ratio=cfg$ess_ratio,
  csv_path=paste0('../original result/',cfg$csv_path),
  results_path=paste(cfg$scale,cfg$analysis,cfg$cell,sep='/'),nsim=1000,input_md5=cfg$md5)
 cat('Plotted',i,'/48:',cfg$scale,cfg$analysis,cfg$cell,'\n');flush.console()
}
all <- do.call(rbind,full);stopifnot(nrow(all)==648L)
write.csv(all,file.path(out,'metrics_all.csv'),row.names=FALSE,na='')
for(scale in c('RR','RD'))write.csv(all[all$scale==scale,,drop=FALSE],
 file.path(out,scale,'metrics_all_runs.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,index),file.path(out,'SCENARIOS.csv'),row.names=FALSE)

stc_file <- file.path(raw,'STC_RR/report/combined_results.csv.gz')
stc <- read.csv(gzfile(stc_file),stringsAsFactors=FALSE)
stopifnot(nrow(stc)==192000L,length(unique(stc$seed))==1000L,
 setequal(unique(stc$uncertainty),c('fixed_bc_summary','bc_mean_delta')))
reference <- read.csv(file.path(raw,'STC_RR/report/method_summary.csv'),stringsAsFactors=FALSE)
reference <- reference[reference$path=='STC',]
stc_all <- list();stc_index <- list();k <- 0L
for(scope in sort(unique(stc$uncertainty))) for(s in sort(unique(stc$scenario_id))) {
 d <- stc[stc$uncertainty==scope & stc$scenario_id==s,,drop=FALSE]
 stopifnot(nrow(d)==12000L,!anyDuplicated(d[c('seed','method','ci_type')]))
 d$component <- 'marginal';d$record_complete <- TRUE
 d$evaluation_truth_saved <- d$population_log_RR_AC_in_BC
 m <- rr$plot_metrics(d,1000L)
 m$scale <- 'STC_RR';m$analysis <- scope;m$scenario_id <- s;m$cell <- unique(d$cell)
 m$event <- unique(d$event);m$n <- unique(d$n_AC);m$ess_ratio <- unique(d$ess_ratio)
 for(j in seq_len(nrow(m))) {
  a <- reference[reference$scenario_id==s & reference$uncertainty==scope & reference$method==m$method[j],]
  stopifnot(nrow(a)==1L)
  map <- c(n_point='points_returned',n_ci='intervals_returned',bias='bias',rmse='rmse',
   empirical_sd='empirical_sd',coverage='coverage_valid',mean_width='mean_width',
   returned_ci_fraction='interval_return_rate',returned_and_covered_fraction='returned_and_covered_rate')
  for(field in names(map))if(!close_number(m[[field]][j],a[[map[[field]]]]))
   stop('STC metric differs from saved report: ',scope,' ',s,' ',m$method[j],' ',field)
 }
 folder <- file.path(out,'STC_RR',scope,unique(d$cell));dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 rr$plot_save(m,folder,'STC marginal A/C log RR in BC target',
  paste(scope,'|',unique(d$cell),'| MC1000'),c('png','pdf'),d,se_max=4,bias_max=5)
 write.csv(m,file.path(folder,'metrics.csv'),row.names=FALSE,na='')
 writeLines(c('Saved STC postprocessing; existing score_unit fits were reused.',
  paste('Uncertainty scope:',scope),'Wald intervals only; no bootstrap inversion was rerun.',
  'Coverage is conditional on returned valid intervals; see return-rate columns.'),file.path(folder,'analysis_notes.txt'))
 k <- k+1L;stc_all[[k]] <- m
 stc_index[[k]] <- data.frame(uncertainty=scope,scenario_id=s,cell=unique(d$cell),
  results_path=paste('STC_RR',scope,unique(d$cell),sep='/'),nsim=1000)
 cat('Plotted STC',k,'/16:',scope,s,'\n');flush.console()
}
write.csv(do.call(rbind,stc_all),file.path(out,'STC_RR/metrics_all_runs.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,stc_index),file.path(out,'STC_RR/SCENARIOS.csv'),row.names=FALSE)
writeLines(c('PASS: all 48 RR/RD cell metrics matched the independent 2026-10-26 CSV audit.',
 'PASS: all 16 STC scope/cell metrics matched the source method_summary.csv STC rows.',
 'No model fitting, simulation or interval inversion was performed.'),file.path(out,'validation/METRIC_CHECKS.txt'))
cat('PASS: 48 RR/RD cells and 16 STC scope/cells plotted and numerically checked.\n')
