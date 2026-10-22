#!/usr/bin/env Rscript
# Standalone cluster analysis; base R only. Does not source or run model code.
# Rscript plot_cluster_results.R RESULT_DIR [--batch=TEXT] [--out=DIR] [--format=png,pdf]

plot_usage <- function() cat(paste(
 'Usage: Rscript plot_cluster_results.R RESULT_DIR [options]',
 'RESULT_DIR may be one run, a group such as weight_all_em_t, or a parent of groups.',
 '--batch=TEXT       Select run directory names containing TEXT (literal match).',
 '--out=DIR          Output root; default RESULT_DIR/analysis_plots_v2.',
 '--format=png,pdf   png, pdf, or png,pdf (default).',
 '--bias-max=NUMBER  Absolute critical value for Bias (default 5); cap only if exceeded.',
 '--se-max=NUMBER    Upper display limit for the SE boxplots (default 4).',
 '--dry-run          List selected runs without reading results or writing files.',
 'Each run is analyzed separately; repeated seeds in different runs are never pooled.',
 sep='\n'), '\n')

plot_options <- function(args) {
 if (!length(args) || any(args %in% c('-h','--help'))) { plot_usage(); return(NULL) }
 if (startsWith(args[1], '--')) stop('The first argument must be a result directory.')
 opt <- list(input=args[1], batch='', out=NULL, formats=c('png','pdf'), dry=FALSE, se_max=4, bias_max=5)
 for (arg in args[-1]) {
  if (arg=='--dry-run') opt$dry <- TRUE
  else if (startsWith(arg,'--batch=')) opt$batch <- substring(arg,9)
  else if (startsWith(arg,'--out=')) opt$out <- substring(arg,7)
  else if (startsWith(arg,'--format=')) opt$formats <- strsplit(substring(arg,10),',',fixed=TRUE)[[1]]
  else if (startsWith(arg,'--bias-max=')) opt$bias_max <- suppressWarnings(as.numeric(substring(arg,12)))
  else if (startsWith(arg,'--se-max=')) opt$se_max <- suppressWarnings(as.numeric(substring(arg,10)))
  else stop('Unknown argument: ',arg)
 }
 if (!dir.exists(opt$input)) stop('Result directory does not exist: ',opt$input)
 opt$input <- normalizePath(opt$input,winslash='/',mustWork=TRUE)
 if (is.null(opt$out)) opt$out <- file.path(opt$input,'analysis_plots_v2')
 if (!nzchar(opt$out)) stop('--out must not be empty.')
 opt$out <- normalizePath(opt$out,winslash='/',mustWork=FALSE)
 if (!length(opt$formats) || any(!opt$formats %in% c('png','pdf'))) stop('Use --format=png,pdf, png, or pdf.')
 if(length(opt$se_max)!=1 || !is.finite(opt$se_max) || opt$se_max<=0) stop('--se-max must be a positive finite number.')
 if(length(opt$bias_max)!=1 || !is.finite(opt$bias_max) || opt$bias_max<=0) stop('--bias-max must be a positive finite number.')
 opt
}

plot_discover <- function(root,batch='') {
 if (dir.exists(file.path(root,'by_seed'))) runs <- root else {
  dirs <- list.dirs(root,recursive=TRUE,full.names=TRUE)
  runs <- dirs[dir.exists(file.path(dirs,'by_seed'))]
 }
 runs <- sort(unique(normalizePath(runs,winslash='/',mustWork=TRUE)))
 if (nzchar(batch)) runs <- runs[grepl(batch,basename(runs),fixed=TRUE)]
 if (!length(runs)) stop('No run directories containing by_seed/ matched the input and --batch filter.')
 runs
}

plot_spec <- function(run) {
 path <- file.path(run,'submission.txt')
 if (!file.exists(path)) return(character())
 lines <- readLines(path,warn=FALSE)
 lines <- lines[grepl('=',lines,fixed=TRUE)]
 keys <- sub('=.*$','',lines)
 if (anyDuplicated(keys)) stop('Duplicate submission.txt fields: ',run)
 setNames(sub('^[^=]*=','',lines),keys)
}

plot_read_run <- function(run) {
 spec <- plot_spec(run)
 files <- sort(list.files(file.path(run,'by_seed'),pattern='^seed_[0-9]+[.]rds$',full.names=TRUE))
 if (!length(files)) stop('No seed_*.rds files found in ',file.path(run,'by_seed'))
 if(any(file.info(files)$size==0) && file.exists(file.path(run,'combined_results.csv')))return(plot_read_collected(run))
 ids <- as.integer(sub('^seed_([0-9]+)[.]rds$','\\1',basename(files)))
 if (anyNA(ids) || anyDuplicated(ids)) stop('Invalid or duplicate seed filenames in ',run)
 known_expected <- all(c('nsim','first_seed') %in% names(spec))
 if (known_expected) {
  n <- suppressWarnings(as.integer(spec[['nsim']]))
  first <- suppressWarnings(as.integer(spec[['first_seed']]))
  if (is.na(n) || n<1 || is.na(first) || first<1) stop('Invalid seed range in submission.txt.')
  expected <- seq.int(first,length.out=n)
 } else expected <- sort(ids)
 audit <- data.frame(seed=sort(unique(c(expected,ids))),expected=sort(unique(c(expected,ids))) %in% expected,
                     state='missing',message='',stringsAsFactors=FALSE)
 audit$state[!audit$expected] <- 'outside_requested_range'
 rows <- list(); conditionals <- list(); hashes <- character()
 keep <- c('method','ci_type','component','status','estimate','se','low','up',
           'se_status','population_log_RR_AC_in_BC','seed','event','n_AC','n_BC','ess_ratio','ci_mode',
           'dgm','dgp_mode','ess_source','em_level','weight_mode','realized_ess','realized_ess_ratio',
           'theoretical_ess_ratio_calibration','model_weighting','weight_covariates',
           'estimation_mode','estimation_em_covariates','estimation_non_em_covariates')
 required <- c('method','ci_type','component','status','estimate','se','low','up','seed',
               'population_log_RR_AC_in_BC')
 selected <- which(ids %in% expected)
 for (j in seq_along(selected)) {
  i <- selected[j]; a <- match(ids[i],audit$seed)
  record <- tryCatch({
   z <- readRDS(files[i]); d <- z$summary
   if (!is.data.frame(d) || !nrow(d) || !all(required %in% names(d))) stop('Unsupported summary schema.')
   if (anyNA(d$seed) || any(d$seed!=ids[i])) stop('Saved seed does not match its filename.')
   if (anyNA(d[c('method','ci_type','component','status')])) stop('Missing method/CI/status metadata.')
   if (anyDuplicated(d[c('seed','method','ci_type')])) stop('Duplicate seed/method/CI rows.')
   for (key in setdiff(keep,names(d))) d[[key]] <- NA
   d <- d[keep]; d$record_complete <- isTRUE(z$complete)
   h <- z$source_md5['estimation/helpers/CI_marginal_target.R']
   list(d=d,conditional=plot_conditional_rows(z),complete=isTRUE(z$complete),hash=if(length(h))as.character(h) else NA_character_)
  },error=function(e)e)
  if (inherits(record,'error')) {
   audit$state[a] <- 'read_error'; audit$message[a] <- conditionMessage(record)
  } else {
   rows[[length(rows)+1L]] <- record$d
   if(!is.null(record$conditional))conditionals[[length(conditionals)+1L]]<-record$conditional
   audit$state[a] <- if(record$complete)'complete' else 'incomplete'
   hashes <- c(hashes,record$hash)
  }
  if (j%%100==0 || j==length(selected)) {
   cat('  read',j,'/',length(selected),'seed files\n'); flush.console()
  }
 }
 if (!length(rows)) stop('No readable result summaries in ',run)
 dat <- do.call(rbind,rows); rownames(dat) <- NULL
 # A run should contain one design; do not silently pool scenarios or CI budgets.
 for (key in c('event','n_AC','n_BC','ess_ratio','ci_mode','dgm','ess_source','em_level','weight_mode')) {
  x <- unique(dat[[key]][!is.na(dat[[key]])])
  if (length(x)>1) stop('Mixed ',key,' values inside one run: ',run)
 }
 list(rows=dat,conditional=if(length(conditionals))do.call(rbind,conditionals) else NULL,
      audit=audit,spec=spec,n_expected=length(expected),
      expected_source=if(known_expected)'submission.txt' else 'observed filenames only',
      exact_hashes=paste(sort(unique(hashes[!is.na(hashes)&nzchar(hashes)])),collapse=';'))
}

plot_target_truth <- function(data) {
 if("evaluation_truth_saved" %in% names(data))return(data$evaluation_truth_saved)
 # The comparison target is BC-population marginal A/C log RR for all methods,
 # including CMH. The old exporter left this named truth column blank in CMH rows.
 # Recover it from other methods in the same seed; never borrow from another seed.
 truth<-rep(NA_real_,nrow(data))
 eligible<-data$component=='marginal' | data$method %in% c('CMH','CMH_sandwich','CMH_an','CMH_sandwich_an')
 for(idx in split(seq_len(nrow(data)),data$seed)) {
  saved<-unique(data$population_log_RR_AC_in_BC[idx])
  saved<-saved[is.finite(saved)]
  if(length(saved)>1)stop('Conflicting saved marginal truths within seed ',data$seed[idx[1]],'.')
  if(length(saved)==1)truth[idx[eligible[idx]]]<-saved
 }
 truth
}

plot_metrics <- function(data,n_expected) {
 data$evaluation_truth<-plot_target_truth(data)
 allowed_point <- c('available','boundary_truncated','interval_failed','se_failed','error','timeout','not_requested')
 keys <- unique(data[c('method','ci_type','component')])
 metrics <- lapply(seq_len(nrow(keys)),function(i) {
  key <- keys[i,];d <- data[data$method==key$method & data$ci_type==key$ci_type &
                          data$component==key$component & data$record_complete,,drop=FALSE]
  point <- d$status %in% allowed_point & is.finite(d$estimate)
  se <- point & d$ci_type=='Wald' & is.finite(d$se) & d$se>0 &
        (is.na(d$se_status) | d$se_status=='available')
  ci <- point & d$status %in% c('available','boundary_truncated') &
        is.finite(d$low) & is.finite(d$up) & d$low<=d$up
  true <- d$evaluation_truth
  truth <- is.finite(true)
  evaluated <- point & truth; covered <- ci & truth & d$low<=true & true<=d$up
  ncover <- sum(ci & truth); ncovered <- sum(covered,na.rm=TRUE)
  mean0 <- function(v) if(length(v))mean(v) else NA_real_
  sd0 <- function(v) if(length(v)>1)sd(v) else NA_real_
  sdall <- sd0(d$estimate[point]); sdse <- sd0(d$estimate[se])
  mse <- mean0(d$se[se]);coverage <- if(ncover)ncovered/ncover else NA_real_
  errors <- d$estimate[evaluated]-true[evaluated]
  cbind(key,data.frame(n_expected=n_expected,n_completed_rows=nrow(d),n_point=sum(point),
   n_se=sum(se),n_ci=sum(ci),n_boundary=sum(ci & d$status=='boundary_truncated'),
   n_evaluated=sum(evaluated),n_coverage=ncover,n_covered=ncovered,
   n_not_applicable=sum(grepl('^not_applicable',d$status)),
   mean_estimate=mean0(d$estimate[point]),mean_truth=mean0(true[evaluated]),
   bias=mean0(errors),bias_mcse=if(length(errors)>1)sd(errors)/sqrt(length(errors)) else NA_real_,
   rmse=if(length(errors))sqrt(mean(errors^2)) else NA_real_,empirical_sd=sdall,
   mean_se=mse,median_se=if(any(se))median(d$se[se]) else NA_real_,
   max_se=if(any(se))max(d$se[se]) else NA_real_,empirical_sd_for_accuracy=sdse,
   accuracy=if(is.finite(sdse)&&sdse>0)mse/sdse else NA_real_,coverage=coverage,
   coverage_mcse=if(ncover)sqrt(coverage*(1-coverage)/ncover) else NA_real_,
   returned_ci_fraction=sum(ci)/n_expected,
   returned_and_covered_fraction=if(any(truth))ncovered/n_expected else NA_real_,
   mean_width=mean0(d$up[ci]-d$low[ci]),stringsAsFactors=FALSE))
 })
 ans <- do.call(rbind,metrics);rownames(ans)<-NULL
 order_methods <- c('brm','brm-FC','LB','LP','RLP','RLB',
                    'CMH','CMH_sandwich','brmad','LB-FC','RLB-FC','LP-FC','RLP-FC','brm_an','LB_an','LP_an','RLP_an','RLB_sandwich_an',
                    'CMH_an','CMH_sandwich_an','brmad_an')
 rank <- match(ans$method,order_methods);rank[is.na(rank)]<-length(order_methods)+1
 ans[order(rank,ans$method,match(ans$ci_type,c('Wald','Exact'))),]
}

plot_samples <- function(data,m) {
 data$evaluation_truth<-plot_target_truth(data)
 allowed_point <- c('available','boundary_truncated','interval_failed','se_failed','error','timeout','not_requested')
 samples <- lapply(seq_len(nrow(m)),function(i) {
  d <- data[data$method==m$method[i] & data$ci_type==m$ci_type[i] &
            data$component==m$component[i] & data$record_complete,,drop=FALSE]
  point <- d$status %in% allowed_point & is.finite(d$estimate)
  evaluated <- point & is.finite(d$evaluation_truth)
  se <- point & d$ci_type=='Wald' & is.finite(d$se) & d$se>0 &
        (is.na(d$se_status) | d$se_status=='available')
  list(bias=d$estimate[evaluated]-d$evaluation_truth[evaluated],se=d$se[se])
 })
 ans <- list(bias=lapply(samples,`[[`,'bias'),se=lapply(samples,`[[`,'se'))
 # A boxplot must contain the same replicates as its summary, never a box of means.
 stopifnot(all(lengths(ans$bias)==m$n_evaluated),all(lengths(ans$se)==m$n_se))
 ans
}
# Publication-style drawing helpers. Base R only; no fitting or simulation.
plot_excluded_method <- function(method) {
 method<-toupper(gsub('_','-',trimws(as.character(method)),fixed=TRUE))
 grepl('^BRMAD',method) | method %in% c('LP-BFC','RLP-BFC')
}

plot_axis <- function(values,metric,se_max=4,bias_max=5) {
 bounds<-switch(metric,bias=c(-bias_max,bias_max),se=c(0,se_max),accuracy=c(0,2),coverage=c(0,1))
 reference<-switch(metric,bias=0,se=NA_real_,accuracy=1,coverage=.95)
 finite<-values[is.finite(values)]
 clipped<-c(any(finite<bounds[1]),any(finite>bounds[2]))
 boxes<-metric %in% c('bias','se')
 # Activate the full critical-value range only when a boxplot observation exceeds it.
 capped<-boxes && any(clipped)
 shown<-if(capped || !boxes)pmax(pmin(finite,bounds[2]),bounds[1]) else finite
 domain<-c(shown,reference[is.finite(reference)])
 if(!length(finite) || capped)limits<-bounds else {
  limits<-range(domain)
  if(diff(limits)==0)limits<-limits+c(-1,1)*max(.05,abs(limits[1])*.05)
  else limits<-limits+c(-.06,.06)*diff(limits)
  if(!boxes)limits<-pmax(pmin(limits,bounds[2]),bounds[1])
  if(metric=='se')limits[1]<-max(0,limits[1])
 }
 ticks<-pretty(limits,n=5)
 if(!capped && length(finite)) {
  limits<-range(c(limits,ticks))
  if(!boxes)limits<-pmax(pmin(limits,bounds[2]),bounds[1])
  if(metric=='se')limits[1]<-max(0,limits[1])
 }
 ticks<-ticks[ticks>=limits[1] & ticks<=limits[2]]
 if(capped)ticks<-sort(unique(c(ticks,bounds)))
 if(length(ticks)<2)ticks<-limits
 list(limits=limits,xlim=limits+c(-.04,.04)*diff(limits),ticks=ticks,
  reference=reference,clipped=clipped,capped=capped,bounds=bounds)
}

plot_metric_alpha <- function(values,metric) {
 target<-switch(metric,coverage=.95,accuracy=1,stop('Unsupported intensity metric.'))
 scale<-if(metric=='coverage').05 else .5
 distance<-abs(values-target)
 alpha<-.18+.82*distance/(distance+scale)
 alpha[!is.finite(distance)]<-.18
 alpha
}

plot_metric_base <- function(metric,ci_type='Wald') {
 if(metric=='coverage')ifelse(ci_type=='Exact','#D26920','#2474A6') else '#444444'
}

plot_metric_colours <- function(values,metric,ci_type=rep('Wald',length(values))) {
 if(!length(values))return(character())
 base<-rep_len(plot_metric_base(metric,ci_type),length(values))
 alpha<-plot_metric_alpha(values,metric)
 vapply(seq_along(values),function(i)adjustcolor(base[i],alpha.f=alpha[i]),character(1))
}

plot_method_order <- c('CMH','CMH_sandwich','LB','LP','RLP','RLB',
 'LB-FC','LP-FC','RLP-FC','RLB-FC','brm','brm-FC')

plot_palette <- function(method) {
 colour <- c(CMH='#A8A8A8',CMH_sandwich='#7F8C8D',LB='#D7AB87',LP='#CEBF86',
  RLP='#ADC18B',RLB='#83B99B','LB-FC'='#84BCA7','LP-FC'='#79BDB4',
  'RLP-FC'='#73BCC8','RLB-FC'='#87AECE',brm='#B5A0D3','brm-FC'='#C19ECC')
 x<-unname(colour[method]);x[is.na(x)]<-'#AAB5C0';x
}

plot_visible <- function(m,metric) {
 m<-m[!plot_excluded_method(m$method),,drop=FALSE]
 if(metric %in% c('bias','se')) {
  # Point estimates are shared by Wald and Exact; show one distribution per method.
  m<-m[order(match(m$ci_type,c('Wald','Exact'))),,drop=FALSE]
  m<-m[!duplicated(m$method),,drop=FALSE]
 } else if(metric=='accuracy') m<-m[m$ci_type=='Wald',,drop=FALSE]
 rank<-match(m$method,plot_method_order);rank[is.na(rank)]<-99
 m[order(rank,m$method,match(m$ci_type,c('Wald','Exact'))),,drop=FALSE]
}

plot_draw <- function(m,metric,title,subtitle,samples,se_max=4,bias_max=5) {
 n<-nrow(m);x<-seq_len(n);boxes<-metric %in% c('bias','se')
 conditional<-any(m$component=='conditional')
 colours<-plot_palette(m$method)
 values<-if(boxes)unlist(samples[[metric]],use.names=FALSE) else m[[metric]]
 ax<-plot_axis(values,metric,se_max,bias_max);reference<-ax$reference
 yrange<-ax$limits;span<-diff(yrange)
 # Space above/below the capped data range is reserved for count markers and labels.
 ylim<-yrange+c(-.12,.15)*span
 ylabels<-c(bias='Estimation error in log(RR)',se='Estimated SE of log(RR)',
  accuracy=if(conditional)'Mean SE / Monte Carlo SD of error' else 'Mean SE / Monte Carlo SD',
  coverage='Coverage rate')
 titles<-c(bias='All-method bias distributions',se='All-method SE distributions',
  accuracy='All-method SE calibration',coverage='All-method coverage rate')
 labels<-m$method;labels[labels=='CMH_sandwich']<-'CMH-SW'
 if(metric=='coverage')labels<-paste(labels,m$ci_type,sep=' / ')
 par(mar=c(8.1,5.2,5.4,1.3),mgp=c(2.7,.75,0),las=1,xpd=FALSE,
  family='sans',cex=1,col.axis='#555555',col.lab='#202020')
 plot(NA,xlim=c(.4,n+.6),ylim=ylim,xlab='',ylab=ylabels[[metric]],
  xaxt='n',yaxt='n',bty='n',xaxs='i',yaxs='i')
 minor<-sort(unique(c(ax$ticks,head(ax$ticks,-1)+diff(ax$ticks)/2)))
 abline(h=minor,col='#F0F0F0',lwd=.8)
 abline(h=ax$ticks,col='#E4E4E4',lwd=1.05)
 abline(v=x,col='#ECECEC',lwd=1)
 axis(2,at=ax$ticks,labels=if(metric=='coverage')sprintf('%.0f%%',100*ax$ticks) else
  format(ax$ticks,trim=TRUE),tick=FALSE,cex.axis=.95)
 if(is.finite(reference))abline(h=reference,lty=2,lwd=1.35,
  col=if(metric=='bias')'#3D7EAD' else '#333333')
 mtext(titles[[metric]],3,line=3.45,adj=0,font=2,cex=1.65)
 mtext(title,3,line=2.05,adj=0,cex=1.03)
 mtext(subtitle,3,line=.7,adj=0,cex=.8,col='#555555')
 par(xpd=NA)
 text(x,par('usr')[3]-.026*diff(par('usr')[3:4]),labels,srt=43,adj=1,cex=.85,col='#555555')
 mtext('Method',1,line=4.65,cex=1.03)
 par(xpd=FALSE)
 counts<-if(metric=='bias')m$n_evaluated else if(metric=='se')m$n_se else
  if(metric=='accuracy')m$n_se else m$n_coverage
 count_line<-paste0('n=',counts)
 fmt<-function(v)ifelse(abs(v)<1e4,sprintf('%.3f',v),formatC(v,digits=3,format='g'))
 if(boxes) {
  keep<-which(lengths(samples[[metric]])>0)
  if(length(keep)) {
   b<-boxplot(samples[[metric]][keep],plot=FALSE,range=1.5)
   b$stats<-pmax(pmin(b$stats,ax$limits[2]),ax$limits[1])
   visible<-b$out>=ax$limits[1] & b$out<=ax$limits[2]
   b$out<-b$out[visible];b$group<-b$group[visible]
   bxp(b,at=x[keep],horizontal=FALSE,add=TRUE,axes=FALSE,frame.plot=FALSE,
    show.names=FALSE,boxwex=.66,boxfill=adjustcolor(colours[keep],alpha.f=.9),
    boxcol='#4D4D4D',medcol='#333333',medlwd=2.1,whiskcol='#555555',
    staplecol='#555555',staplewex=0,whisklty=1,outcol=adjustcolor('#444444',alpha.f=.42),
    outpch=16,outcex=.35)
  }
  means<-if(metric=='bias')m$bias else m$mean_se
  points(x,pmax(pmin(means,ax$limits[2]),ax$limits[1]),pch=21,bg='white',
   col='#333333',cex=.76,lwd=.9)
  for(side in c(-1,1)) {
   if(metric=='se' && side<0)next
   bound<-if(side<0)-bias_max else if(metric=='bias')bias_max else se_max
   clipped<-vapply(samples[[metric]],function(v)sum(side*v>side*bound),integer(1))
   ii<-which(clipped>0)
   if(length(ii)) {
    points(x[ii],rep(bound,length(ii)),pch=21,bg='#D94C42',col='#A7322A',cex=1.9)
    text(x[ii],rep(bound,length(ii)),labels=clipped[ii],col='white',cex=.62,font=2)
   }
  }
  absent<-lengths(samples[[metric]])==0
  note<-'White points: full-data means. Red circles: display-clipped counts; statistics use all valid values.'
 } else {
  ok<-is.finite(m[[metric]]);v<-m[[metric]]
  stem_colours<-plot_metric_colours(v,metric,m$ci_type)
  shown<-pmax(pmin(v,ax$limits[2]),ax$limits[1]);clipped<-ok & v>ax$limits[2]
  if(any(ok)) {
   segments(x[ok],reference,x[ok],shown[ok],lwd=5.3,col=stem_colours[ok])
   points(x[ok],shown[ok],pch=ifelse(m$ci_type[ok]=='Exact',17,16),col='#222222',cex=.92)
  }
  labs<-if(metric=='coverage')sprintf('%.1f%%',100*v) else fmt(v)
  labs[clipped]<-paste0(labs[clipped],' >')
  above<-shown>=reference
  if(any(ok))text(x[ok],shown[ok]+ifelse(above[ok],.036,-.04)*span,labs[ok],cex=.73,
   col=ifelse(clipped[ok],'#A7322A','#333333'))
  absent<-!ok
  note<-if(metric=='accuracy')paste('Target = 1.',
   if(conditional)'Matched replicates; denominator SD(error).' else 'Matched SE-valid replicates; Exact SE not provided.',
   'Darker stems = larger deviation.') else
   'Target = 95%. Wald: blue; Exact: orange. Darker stems = larger deviation. Valid returned intervals only.'
 }
 if(any(absent))text(x[absent],yrange[1]+.06*span,'NA',col='#8A8A8A',cex=.8)
 text(x,rep(yrange[1]-.077*span,n),count_line,cex=.64,col='#777777')
 mtext(note,1,line=6.1,cex=.73,adj=0)
 if(metric=='bias')mtext(paste0('Bias reference = 0; critical values +/-',bias_max,
  if(ax$capped)'; capped range activated.' else '; adaptive axis, no observations clipped.'),
  1,line=6.85,cex=.7,adj=0)
 if(metric=='se')mtext(paste0('SE critical value = ',se_max,
  if(ax$capped)'; capped range activated.' else '; adaptive axis, no observations clipped.'),
  1,line=6.85,cex=.7,adj=0)
}

plot_save <- function(metrics,folder,title,subtitle,formats,rows,se_max=4,bias_max=5) {
 items<-c('bias','se','accuracy','coverage')
 draw<-function(metric) {
  mm<-plot_visible(metrics,metric)
  plot_draw(mm,metric,title,subtitle,plot_samples(rows,mm),se_max,bias_max)
 }
 if('png'%in%formats)for(metric in items) {
  type<-if(capabilities('cairo'))'cairo' else if(.Platform$OS.type=='windows')'windows' else
   stop('Headless PNG needs Cairo support. Use --format=pdf.')
  png(file.path(folder,paste0(metric,'.png')),width=4200,height=2220,res=300,type=type)
  tryCatch(draw(metric),finally=dev.off())
 }
 if('pdf'%in%formats) {
  pdf(file.path(folder,'metrics.pdf'),width=14,height=7.4,onefile=TRUE)
  tryCatch(for(metric in items)draw(metric),finally=dev.off())
 }
}

plot_conditional_rows <- function(z) {
 cf<-z$coefficients;ref<-z$result$reference_truth
 if(!is.data.frame(cf) || !nrow(cf))return(NULL)
 cf<-cf[(cf$component=='alpha' & cf$term=='(Intercept)') |
         (cf$component=='treatment' & cf$term=='x'),,drop=FALSE]
 cf<-cf[!plot_excluded_method(cf$method),,drop=FALSE]
 if(!nrow(cf))return(NULL)
 keys<-function(d)paste(d$method,d$component,d$term,sep='\034')
 true<-rep(NA_real_,nrow(cf))
 if(is.data.frame(ref) && all(c('method','component','term','truth')%in%names(ref))) {
  if(anyDuplicated(keys(ref)))stop('Duplicate conditional reference-truth keys.')
  true<-ref$truth[match(keys(cf),keys(ref))]
 }
 interval<-cf$interval_status
 interval[is.na(interval)]<-'interval_failed'
 status<-ifelse(cf$status=='available',
  ifelse(interval=='available','available',ifelse(interval=='BC_boundary_truncated','boundary_truncated','interval_failed')),
  cf$status)
 data.frame(method=cf$method,ci_type=ifelse(cf$interval=='Wald','Wald','Exact'),
  component='conditional',status=status,estimate=cf$estimate,se=cf$se,
  low=cf$lower,up=cf$upper,se_status=cf$se_status,seed=cf$seed,
  record_complete=isTRUE(z$complete),evaluation_truth_saved=true,
  population_log_RR_AC_in_BC=NA_real_,stringsAsFactors=FALSE)
}

plot_read_collected <- function(run) {
 spec<-plot_spec(run)
 d<-read.csv(file.path(run,'combined_results.csv'),stringsAsFactors=FALSE)
 required<-c('method','ci_type','component','status','estimate','se','low','up','seed',
  'population_log_RR_AC_in_BC','seed_complete')
 if(!nrow(d) || !all(required %in% names(d)))stop('Unsupported collected CSV schema.')
 if(anyNA(d[c('seed','method','ci_type','component','status','seed_complete')]) ||
  anyDuplicated(d[c('seed','method','ci_type')]))stop('Missing/duplicate collected CSV keys.')
 if(!is.logical(d$seed_complete))stop('CSV seed_complete must be an explicit logical flag.')
 if(any(d$seed!=as.integer(d$seed)) || any(d$seed<1))stop('Invalid CSV seed IDs.')
 perseed<-split(d$seed_complete,d$seed)
 if(any(vapply(perseed,function(x)length(unique(x))!=1L,logical(1))))stop('Conflicting CSV seed_complete flags.')
 known<-all(c('nsim','first_seed') %in% names(spec))
 if(known) {
  n<-as.integer(spec[['nsim']]);first<-as.integer(spec[['first_seed']])
  if(is.na(n) || n<1 || is.na(first) || first<1)stop('Invalid submission seed range.')
  expected<-seq.int(first,length.out=n)
 } else expected<-sort(unique(d$seed))
 seeds<-sort(unique(c(expected,d$seed)))
 audit<-data.frame(seed=seeds,expected=seeds %in% expected,state='missing',
  message='Local RDS empty or unavailable; validated combined_results.csv used.',rds_state='empty',stringsAsFactors=FALSE)
 for(i in seq_along(seeds)) {
  dd<-d[d$seed==seeds[i],,drop=FALSE]
  if(!audit$expected[i])audit$state[i]<-'outside_requested_range'
  else if(nrow(dd))audit$state[i]<-if(dd$seed_complete[1])'complete_from_csv' else 'incomplete_from_csv'
 }
 d<-d[d$seed %in% expected,,drop=FALSE];d$record_complete<-d$seed_complete
 files<-list.files(file.path(run,'by_seed'),pattern='^seed_[0-9]+[.]rds$',full.names=TRUE)
 files<-files[file.info(files)$size>0]
 conditionals<-list()
 rowkey<-function(x)paste(x$seed,x$method,x$ci_type,sep='\034')
 for(f in files) {
  zz<-readRDS(f);ss<-zz$summary
  if(!is.data.frame(ss) || !all(setdiff(required,'seed_complete') %in% names(ss)))stop('Unsupported nonempty RDS in CSV recovery.')
  ids<-match(rowkey(ss),rowkey(d))
  if(anyNA(ids))stop('CSV recovery: nonempty RDS keys do not match collected CSV.')
  fields<-c('method','ci_type','component','status','estimate','se','low','up','population_log_RR_AC_in_BC')
  if(!isTRUE(all.equal(ss[fields],d[ids,fields],check.attributes=FALSE,tolerance=1e-12)))
   stop('CSV recovery: nonempty RDS values differ from collected CSV.')
  if(any(d$seed_complete[ids]!=isTRUE(zz$complete)))stop('CSV/RDS completeness flags differ.')
  aa<-match(unique(ss$seed),audit$seed);audit$rds_state[aa]<-'available'
  audit$message[aa]<-'Nonempty local RDS matched against collected CSV; conditional truths read from RDS.'
  cf<-plot_conditional_rows(zz)
  if(!is.null(cf))conditionals[[length(conditionals)+1L]]<-cf
 }
 keep<-c('event','n_AC','n_BC','ess_ratio','ci_mode','dgm','dgp_mode','ess_source','em_level',
  'weight_mode','realized_ess','realized_ess_ratio','se_status','model_weighting','weight_covariates',
  'estimation_mode','estimation_em_covariates','estimation_non_em_covariates')
 for(k in setdiff(keep,names(d)))d[[k]]<-NA
 for(k in setdiff(keep,c('realized_ess','realized_ess_ratio','se_status'))) {
  if(length(unique(d[[k]][!is.na(d[[k]])]))>1L)stop('Mixed CSV design field: ',k)
 }
 snap<-file.path(run,'code_snapshot','estimation','helpers','CI_marginal_target.R')
 h<-if(file.exists(snap))tolower(unname(tools::md5sum(snap))) else ''
 cat('  using validated collected CSV;',length(files),'nonempty RDS verified; missing coefficient truths not invented.\n')
 list(rows=d,conditional=if(length(conditionals))do.call(rbind,conditionals) else NULL,audit=audit,spec=spec,n_expected=length(expected),
  expected_source=if(known)'submission.txt' else 'collected CSV seed IDs',
  exact_hashes=h,input_source='combined_results.csv; empty local RDS; nonempty RDS cross-validated')
}

plot_metadata <- function(m,z,run,target='marginal') {
 m$run_name<-basename(run);m$run_path<-run
 m$evaluation_target<-if(target=='marginal')'BC population marginal A/C log RR' else
  'Conditional A/C log RR at h_E=0 in each replicate analysis units'
 m$truth_source<-if(target=='marginal')'saved population_log_RR_AC_in_BC, matched within seed' else
  'saved result$reference_truth, matched within seed by method/component/term'
 keys<-c('event','n_AC','n_BC','ess_ratio','ci_mode','dgm','dgp_mode',
  'ess_source','em_level','weight_mode','model_weighting','weight_covariates',
  'estimation_mode','estimation_em_covariates','estimation_non_em_covariates')
 for(k in keys) {
  v<-unique(z$rows[[k]][!is.na(z$rows[[k]])])
  if(length(v)>1)stop('Mixed ',k,' in one run; comparison metadata are ambiguous.')
  m[[k]]<-if(length(v))as.character(v[1]) else if(k %in% names(z$spec))z$spec[[k]] else NA_character_
 }
 m$dataset<-basename(dirname(run))
 m$exact_source_md5<-z$exact_hashes
 for(k in c('n_AC','n_BC','ess_ratio'))m[[k]]<-as.numeric(m[[k]])
 good<-z$rows$record_complete & !duplicated(z$rows$seed)
 x<-z$rows$realized_ess[good];r<-z$rows$realized_ess_ratio[good]
 m$mean_realized_ess<-if(any(is.finite(x)))mean(x[is.finite(x)]) else NA_real_
 m$mean_realized_ess_ratio<-if(any(is.finite(r)))mean(r[is.finite(r)]) else NA_real_
 m$point_return_fraction<-m$n_point/m$n_expected
 m
}

plot_main <- function(args=commandArgs(TRUE)) {
 opt<-plot_options(args);if(is.null(opt))return(invisible(NULL))
 runs<-plot_discover(opt$input,opt$batch)
 cat('Selected',length(runs),'run(s); brmad, LP-BFC and RLP-BFC are excluded.\n')
 if(opt$dry){cat(paste(runs,collapse='\n'),'\n');return(invisible(runs))}
 dir.create(opt$out,recursive=TRUE,showWarnings=FALSE)
 opt$out<-normalizePath(opt$out,winslash='/',mustWork=TRUE)
 all_metrics<-conditional_metrics<-index<-list()
 for(run in runs) {
  cat('\nAnalyzing:',run,'\n');flush.console()
  relative<-if(identical(run,opt$input))basename(run) else substring(run,nchar(opt$input)+2)
  folder<-file.path(opt$out,'by_run',relative)
  status<-tryCatch({
   z<-plot_read_run(run)
   z$rows<-z$rows[!plot_excluded_method(z$rows$method),,drop=FALSE]
   m<-plot_metadata(plot_metrics(z$rows,z$n_expected),z,run)
   dir.create(folder,recursive=TRUE,showWarnings=FALSE)
   write.csv(z$audit,file.path(folder,'seed_audit.csv'),row.names=FALSE,na='')
   d<-z$rows[z$rows$record_complete,,drop=FALSE]
   statuses<-if(nrow(d))aggregate(list(count=rep(1L,nrow(d))),d[c('method','ci_type','status')],sum) else
    data.frame(method=character(),ci_type=character(),status=character(),count=integer())
   write.csv(statuses,file.path(folder,'status_counts.csv'),row.names=FALSE)
   write.csv(m,file.path(folder,'metrics.csv'),row.names=FALSE,na='')
   ncomplete<-sum(z$audit$state %in% c('complete','complete_from_csv'))
   title<-sprintf('Marginal A/C | %s | n = %s | ESS ratio = %s',m$event[1],m$n_AC[1],m$ess_ratio[1])
   subtitle<-sprintf('Adjustment: %s | fit weights: %s | calibration: %s | complete %d / %d',
    m$estimation_mode[1],m$model_weighting[1],m$weight_mode[1],ncomplete,z$n_expected)
   plot_save(m,folder,title,subtitle,opt$formats,z$rows,opt$se_max,opt$bias_max)
   cnote<-'Conditional plots skipped: no saved matching conditional truth.'
   if(is.data.frame(z$conditional) && any(is.finite(z$conditional$evaluation_truth_saved))) {
    cm<-plot_metrics(z$conditional,z$n_expected)
    # A replicate-specific coefficient truth varies with score centering/scale.
    # Calibration must use the SD of estimation errors in the same SE-valid subset.
    for(i in seq_len(nrow(cm))) {
     dd<-z$conditional[z$conditional$method==cm$method[i] & z$conditional$ci_type==cm$ci_type[i] &
      z$conditional$record_complete,,drop=FALSE]
     valid<-dd$status %in% c('available','boundary_truncated','interval_failed','se_failed','error','timeout','not_requested') &
      is.finite(dd$estimate) & is.finite(dd$evaluation_truth_saved) & is.finite(dd$se) & dd$se>0 &
      dd$ci_type=='Wald' & (is.na(dd$se_status) | dd$se_status=='available')
     er<-dd$estimate[valid]-dd$evaluation_truth_saved[valid]
     ss<-if(length(er)>1)sd(er) else NA_real_
     cm$empirical_sd_for_accuracy[i]<-ss
     cm$accuracy[i]<-if(is.finite(ss) && ss>0)mean(dd$se[valid])/ss else NA_real_
    }
    cm<-plot_metadata(cm,z,run,'conditional')
    cf<-file.path(folder,'conditional');dir.create(cf,showWarnings=FALSE)
    write.csv(cm,file.path(cf,'metrics.csv'),row.names=FALSE,na='')
    truth_seeds<-length(unique(z$conditional$seed[is.finite(z$conditional$evaluation_truth_saved)]))
    csubtitle<-paste0(subtitle,' | saved truth seeds ',truth_seeds,' / ',z$n_expected)
    plot_save(cm,cf,sub('Marginal','Conditional at h_E = 0',title),csubtitle,opt$formats,z$conditional,opt$se_max,opt$bias_max)
    conditional_metrics[[length(conditional_metrics)+1L]]<-cm
    cnote<-'Conditional plots use the saved per-seed coefficient truth; they are a different target from marginal A/C.'
   }
   writeLines(c('Saved-result analysis only; no models fitted. brmad, LP-BFC and RLP-BFC (including underscore aliases) excluded.',
    paste('Input:',run),paste('Expected seed range source:',z$expected_source),
    paste('Data source:',if(is.null(z$input_source))'by_seed RDS' else z$input_source),
    'Only complete saved records enter metrics; expected denominators retain missing/incomplete/unreadable seeds.',
    'Marginal truth: BC population A/C log RR, matched within seed, including CMH.',
    'Point/SE/CI status filters and marginal summary calculations follow the original script.',
    'Point distributions display one box per method (Wald preferred); Exact duplicates are omitted from boxplots.',
    'Exact SE and SE calibration are not supplied and are omitted from SE/calibration figures.',
    'Bias and SE critical ranges activate only when finite observations exceed the thresholds; otherwise axes are adaptive.',
    'Accuracy stems share one grey hue. Coverage stems share a blue Wald hue or an orange Exact hue.',
    'Stem opacity increases monotonically with absolute deviation from 1 (accuracy) or 0.95 (coverage).',
    'White boxplot points are full-data means. Red count circles represent finite values beyond display caps only.',
    'No trimming or winsorization in statistics. Missing/invalid values are not counted as clipped observations.',
    'Coverage conditions on returned valid intervals. Also inspect returned_ci_fraction and returned_and_covered_fraction.',
    'PDF contains four pages; PNGs are 4200 x 2220 at 300 dpi.',cnote,
    'Runs remain separate. Do not pool batches, fitting weights, calibration sets, or different estimands.'),
    file.path(folder,'analysis_notes.txt'))
   all_metrics[[length(all_metrics)+1L]]<-m
   data.frame(run_name=basename(run),run_path=run,output_dir=folder,status='analyzed',
    n_expected=z$n_expected,n_complete=ncomplete,n_missing=sum(z$audit$state=='missing'),
    n_incomplete=sum(z$audit$state %in% c('incomplete','incomplete_from_csv')),n_read_error=sum(z$audit$state=='read_error'),
    conditional_truth_available=grepl('^Conditional plots use',cnote),message='',stringsAsFactors=FALSE)
  },error=function(e)data.frame(run_name=basename(run),run_path=run,output_dir=folder,status='error',
   n_expected=NA,n_complete=NA,n_missing=NA,n_incomplete=NA,n_read_error=NA,
   conditional_truth_available=FALSE,message=conditionMessage(e),stringsAsFactors=FALSE))
  index[[length(index)+1L]]<-status
  cat(' ',status$status,':',if(status$status=='error')status$message else folder,'\n');flush.console()
 }
 ix<-do.call(rbind,index)
 write.csv(ix,file.path(opt$out,'run_index.csv'),row.names=FALSE,na='')
 if(length(all_metrics))write.csv(do.call(rbind,all_metrics),file.path(opt$out,'metrics_all_runs.csv'),row.names=FALSE,na='')
 if(length(conditional_metrics))write.csv(do.call(rbind,conditional_metrics),file.path(opt$out,'conditional_metrics_all_runs.csv'),row.names=FALSE,na='')
 writeLines(c('Publication-style plots, base R. brmad, LP-BFC and RLP-BFC excluded.',
  'Each by_run folder: bias/se/accuracy/coverage PNG, four-page metrics.pdf, metrics.csv and seed audit.',
  'conditional/ exists only when saved matching coefficient truths are available.',
  'metrics_all_runs.csv concatenates runs, without pooling seeds.',
  'See run_index.csv for all selected runs and errors.'),file.path(opt$out,'README.txt'))
 if(any(ix$status=='error'))stop(sum(ix$status=='error'),' run(s) failed; see run_index.csv.')
 invisible(list(index=ix,metrics=if(length(all_metrics))do.call(rbind,all_metrics) else data.frame()))
}
if(sys.nframe()==0L)plot_main()
