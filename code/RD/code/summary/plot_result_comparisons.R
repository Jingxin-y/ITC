#!/usr/bin/env Rscript
# Matched-design comparison figures; base R only. Never pools simulation batches.
# Rscript plot_result_comparisons.R PLOTS_V2_DIR [--out=DIR] [--batch=TEXT]
# [--focus=CMH,brm,brm-FC,RLP,RLP-FC]
args<-commandArgs(TRUE)
if(!length(args) || any(args %in% c('--help','-h'))) {
 cat('Rscript plot_result_comparisons.R PLOTS_V2_DIR [--out=DIR] [--batch=TEXT]\n',
  '[--focus=CMH,brm,brm-FC,RLP,RLP-FC]\n')
 quit(status=0)
}
script<-sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))][1])
source(file.path(dirname(normalizePath(script,winslash='/')),'plot_cluster_results.R'))
root<-normalizePath(args[1],winslash='/',mustWork=TRUE)
out<-file.path(root,'comparisons');batch<-''
focus<-c('CMH','brm','brm-FC','RLP','RLP-FC')
for(a in args[-1]) {
 if(startsWith(a,'--out='))out<-substring(a,7)
 else if(startsWith(a,'--batch='))batch<-substring(a,9)
 else if(startsWith(a,'--focus='))focus<-strsplit(substring(a,9),',',fixed=TRUE)[[1]]
 else stop('Unknown argument: ',a)
}
focus<-focus[!plot_excluded_method(focus)]
if(!length(focus))stop('At least one focus method is required.')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
metrics<-read.csv(file.path(root,'metrics_all_runs.csv'),stringsAsFactors=FALSE)
metrics<-metrics[!plot_excluded_method(metrics$method) & metrics$n_expected>=100,,drop=FALSE]
if(nzchar(batch))metrics<-metrics[grepl(batch,metrics$run_name,fixed=TRUE),,drop=FALSE]
if(!nrow(metrics))stop('No non-smoke scenarios matched.')
design<-c('dataset','dgm','dgp_mode','ess_source','em_level','weight_mode','model_weighting',
 'weight_covariates','estimation_mode','estimation_em_covariates','estimation_non_em_covariates',
 'ci_mode','exact_source_md5','ci_type')
if(!all(design %in% names(metrics)))stop('Use metrics produced by the revised plot_cluster_results.R.')
key<-function(d,cols)do.call(paste,c(lapply(d[cols],function(x)ifelse(is.na(x),'<missing>',as.character(x))),list(sep='\034')))
groups<-split(metrics,key(metrics,design))
audit<-list();manifest<-list()
metric_labels<-c(bias='Mean bias',rmse='RMSE',accuracy='Mean SE / Monte Carlo SD',coverage='Coverage rate',
 returned_ci_fraction='CI return fraction',point_return_fraction='Point return fraction')
ref<-function(metric)switch(metric,bias=0,accuracy=1,coverage=.95,NA_real_)
col_method<-c(CMH='#555555',brm='#0077BB','brm-FC'='#009988',RLP='#EE7733','RLP-FC'='#AA4499')
method_colour<-function(mm) {
 v<-unname(col_method[mm]);bad<-is.na(v);v[bad]<-plot_palette(mm[bad]);v
}
series_lines<-function(x,values,metric,ci_type,col,lty=1) {
 shown<-if(metric=='accuracy')pmin(values,2) else values
 if(!metric %in% c('accuracy','coverage')) {
  lines(x,shown,col=col,lty=lty,lwd=1.5)
  return(rep_len(col,length(values)))
 }
 point_colours<-plot_metric_colours(values,metric,ci_type)
 n<-length(values)
 if(n>1) {
  left<-seq_len(n-1L);right<-left+1L
  ok<-is.finite(x[left]) & is.finite(x[right]) & is.finite(shown[left]) & is.finite(shown[right])
  left<-left[ok];right<-right[ok]
  alpha<-plot_metric_alpha(values,metric)
  bases<-rep_len(plot_metric_base(metric,ci_type),n)
  if(length(left)) {
   colours<-vapply(seq_along(left),function(i)adjustcolor(bases[left[i]],
    alpha.f=mean(alpha[c(left[i],right[i])])),character(1))
   segments(x[left],shown[left],x[right],shown[right],col=colours,lty=lty,lwd=1.5)
  }
 }
 point_colours
}
limits<-function(d,metric) {
 v<-d[[metric]];se<-if(metric=='bias')d$bias_mcse else if(metric=='coverage')d$coverage_mcse else rep(0,nrow(d))
 v<-c(v,v-1.96*se,v+1.96*se,ref(metric));v<-v[is.finite(v)]
 if(!length(v))return(c(0,1))
 if(metric=='accuracy' && any(v>2))return(c(0,2.15))
 lim<-range(v)
 if(metric %in% c('coverage','returned_ci_fraction','point_return_fraction'))return(c(0,1.05))
 if(diff(lim)==0)lim<-lim+c(-.05,.05)
 lim+c(-.08,.1)*diff(lim)
}
grid<-function(xlim,ylim,metric,xlab='',title='') {
 plot(NA,xlim=xlim,ylim=ylim,xaxt='n',yaxt='n',bty='n',xlab=xlab,
  ylab=metric_labels[[metric]],main=title,cex.main=1.02)
 yy<-pretty(ylim,n=5);yy<-yy[yy>=ylim[1] & yy<=ylim[2]]
 abline(h=yy,col='#E8E8E8')
 axis(2,at=yy,labels=if(metric %in% c('coverage','returned_ci_fraction','point_return_fraction'))
  sprintf('%.0f%%',yy*100) else format(yy,trim=TRUE),tick=FALSE,cex.axis=.8)
 if(is.finite(ref(metric)))abline(h=ref(metric),lty=2,col='#666666')
}
errorbar<-function(x,d,metric,col) {
 if(!metric %in% c('bias','coverage'))return(invisible(NULL))
 ss<-if(metric=='bias')d$bias_mcse else d$coverage_mcse
 ok<-is.finite(d[[metric]]) & is.finite(ss) & ss>0
 lo<-d[[metric]][ok]-1.96*ss[ok];hi<-d[[metric]][ok]+1.96*ss[ok]
 if(metric=='coverage'){lo<-pmax(0,lo);hi<-pmin(1,hi)}
 if(any(ok))arrows(x[ok],lo,x[ok],hi,angle=90,code=3,length=.035,
  col=if(length(col)>1)col[ok] else col,lwd=.85)
}
save_figure<-function(stem,width,height,draw) {
 png(paste0(stem,'.png'),width=width*300,height=height*300,res=300,
  type=if(capabilities('cairo'))'cairo' else 'windows')
 tryCatch(draw(),finally=dev.off())
 pdf(paste0(stem,'.pdf'),width=width,height=height)
 tryCatch(draw(),finally=dev.off())
}
event_dashboard<-function(d,title,levels=c('common','rare'),labels=c('Common','Rare'),
 heading='Common versus rare events',extra_note='') {
 methods<-plot_method_order[plot_method_order %in% d$method]
 methods<-c(methods,setdiff(unique(d$method),methods))
 par(mfrow=c(2,2),mar=c(6.1,4.2,2,1),oma=c(2.8,0,4.3,0),mgp=c(2.5,.6,0),family='sans',las=1)
 items<-if(all(d$ci_type=='Exact'))c('bias','rmse','coverage','returned_ci_fraction') else
  c('bias','rmse','accuracy','coverage')
 for(metric in items) {
  grid(c(.4,length(methods)+.6),limits(d,metric),metric)
  abline(v=seq_along(methods),col='#EFEFEF')
  for(j in 1:2) {
   ev<-levels[j];dd<-d[d$event==ev,,drop=FALSE]
   xx<-match(dd$method,methods)+c(-.12,.12)[j];cc<-c('#0077BB','#EE7733')[j]
   if(metric %in% c('accuracy','coverage'))cc<-plot_metric_colours(dd[[metric]],metric,dd$ci_type)
   errorbar(xx,dd,metric,cc)
   vv<-if(metric=='accuracy')pmin(dd[[metric]],2) else dd[[metric]]
   points(xx,vv,pch=c(16,17)[j],col=cc,cex=.8)
   if(metric=='accuracy') {
    clipped<-is.finite(dd[[metric]]) & dd[[metric]]>2
    if(any(clipped)) {
     clip_col<-if(length(cc)>1)cc[clipped] else cc
     points(xx[clipped],vv[clipped],pch=24,col=clip_col,bg=clip_col,cex=.9)
     text(xx[clipped],vv[clipped]+.085,formatC(dd[[metric]][clipped],digits=3,format='g'),col=clip_col,cex=.63)
    }
   }
  }
  par(xpd=NA);text(seq_along(methods),par('usr')[3]-.025*diff(par('usr')[3:4]),
   sub('CMH_sandwich','CMH-SW',methods,fixed=TRUE),srt=45,adj=1,cex=.66,col='#555555');par(xpd=FALSE)
 }
 mtext(heading,outer=TRUE,3,line=2.85,adj=.05,font=2,cex=1.35)
 mtext(title,outer=TRUE,3,line=1.5,adj=.05,cex=.85)
 par(mfrow=c(1,1),fig=c(0,1,0,1),mar=rep(0,4),oma=rep(0,4),new=TRUE);plot.new()
 legend(.5,.97,xjust=.5,yjust=1,legend=labels,col=c('#333333','#333333'),
  pch=c(16,17),horiz=TRUE,bty='n',cex=.85)
 text(.035,.017,paste('Bars: +/-1.96 MCSE. Calibration >2 capped. Accuracy: grey; coverage: blue Wald/orange Exact; darker = larger deviation.',extra_note),adj=0,cex=.6)
}
trend_figure<-function(d,metric,title) {
 focus<-intersect(focus,unique(d$method))
 events<-c('common','rare');events<-events[events %in% d$event]
 esses<-sort(unique(d$ess_ratio));ns<-sort(unique(d$n_AC))
 lim<-limits(d,metric)
 par(mfrow=c(length(esses),length(events)),mar=c(3.5,4.1,2.3,1),
  oma=c(4.5,0,6.5,0),mgp=c(2.5,.65,0),family='sans',las=1)
 for(ee in esses)for(ev in events) {
  dd<-d[d$ess_ratio==ee & d$event==ev,,drop=FALSE]
  grid(range(ns)+c(-.07,.07)*max(1,diff(range(ns))),lim,metric,
   xlab='Sample size n',title=paste(ev,'| ESS ratio =',ee))
  abline(v=ns,col='#EFEFEF');axis(1,at=ns,tick=FALSE,cex.axis=.85)
  for(j in seq_along(focus)) {
   mm<-dd[dd$method==focus[j],,drop=FALSE];mm<-mm[order(mm$n_AC),,drop=FALSE]
   if(!nrow(mm))next
   vv<-if(metric=='accuracy')pmin(mm[[metric]],2) else mm[[metric]]
   cc<-series_lines(mm$n_AC,mm[[metric]],metric,mm$ci_type,method_colour(focus[j]),lty=j)
   errorbar(mm$n_AC,mm,metric,cc)
   points(mm$n_AC,vv,col=cc,pch=c(16,17,15,18,8)[(j-1)%%5+1],cex=.85)
   if(metric=='accuracy') {
    clipped<-is.finite(mm[[metric]]) & mm[[metric]]>2
    if(any(clipped)) {
     clip_col<-cc[clipped]
     points(mm$n_AC[clipped],vv[clipped],pch=24,col=clip_col,bg=clip_col,cex=.9)
     text(mm$n_AC[clipped],vv[clipped]+.08,formatC(mm[[metric]][clipped],digits=3,format='g'),col=clip_col,cex=.65)
    }
   }
  }
 }
 mtext(paste(metric_labels[[metric]],'by sample size and ESS'),outer=TRUE,3,line=4.75,adj=.055,font=2,cex=1.3)
 mtext(title,outer=TRUE,3,line=3.35,adj=.055,cex=.8)
 par(mfrow=c(1,1),fig=c(0,1,0,1),mar=rep(0,4),oma=rep(0,4),new=TRUE);plot.new()
 legend_col<-if(metric %in% c('accuracy','coverage'))rep_len(plot_metric_base(metric,d$ci_type[1]),length(focus)) else method_colour(focus)
 legend(.5,.956,xjust=.5,yjust=1,legend=focus,col=legend_col,
  lty=seq_along(focus),pch=c(16,17,15,18,8)[(seq_along(focus)-1)%%5+1],horiz=TRUE,bty='n',cex=.82)
 text(.035,.017,if(metric %in% c('accuracy','coverage'))'Shared hue within each CI type; darker = larger deviation from target. Methods differ by line type and shape.' else
  'Panels share the same y-axis. Focus methods are shown; all methods remain in detail figures and CSV.',
  adj=0,cex=.7)
}
ess_figure<-function(d,metric,title) {
 focus<-intersect(focus,unique(d$method))
 events<-c('common','rare');events<-events[events %in% d$event]
 ns<-sort(unique(d$n_AC));lim<-limits(d,metric)
 par(mfrow=c(length(events),length(ns)),mar=c(3.5,4.1,2.3,1),
  oma=c(4.5,0,6.5,0),mgp=c(2.5,.65,0),family='sans',las=1)
 for(ev in events)for(nn in ns) {
  dd<-d[d$event==ev & d$n_AC==nn,,drop=FALSE]
  xr<-range(d$mean_realized_ess[d$n_AC==nn],na.rm=TRUE)
  xr<-xr+c(-.07,.07)*max(1,diff(xr))
  grid(xr,lim,metric,xlab='Mean realized ESS',title=paste(ev,'| n =',nn))
  axis(1,at=pretty(xr,n=4),tick=FALSE,cex.axis=.85)
  for(j in seq_along(focus)) {
   mm<-dd[dd$method==focus[j] & is.finite(dd$mean_realized_ess),,drop=FALSE]
   mm<-mm[order(mm$mean_realized_ess),,drop=FALSE]
   if(!nrow(mm))next
   cc<-series_lines(mm$mean_realized_ess,mm[[metric]],metric,mm$ci_type,method_colour(focus[j]),lty=j)
   errorbar(mm$mean_realized_ess,mm,metric,cc)
   points(mm$mean_realized_ess,mm[[metric]],col=cc,pch=c(16,17,15,18,8)[(j-1)%%5+1],cex=.85)
  }
 }
 mtext(paste(metric_labels[[metric]],'by realized ESS'),outer=TRUE,3,line=4.75,adj=.055,font=2,cex=1.3)
 mtext(title,outer=TRUE,3,line=3.35,adj=.055,cex=.8)
 par(mfrow=c(1,1),fig=c(0,1,0,1),mar=rep(0,4),oma=rep(0,4),new=TRUE);plot.new()
 legend_col<-if(metric=='coverage')rep_len(plot_metric_base(metric,d$ci_type[1]),length(focus)) else method_colour(focus)
 legend(.5,.956,xjust=.5,yjust=1,legend=focus,col=legend_col,
  lty=seq_along(focus),pch=c(16,17,15,18,8)[(seq_along(focus)-1)%%5+1],horiz=TRUE,bty='n',cex=.82)
 text(.035,.017,'Sample size is held fixed within each panel. ESS is the mean across complete saved replicates.',
  adj=0,cex=.7)
}
for(i in seq_along(groups)) {
 d<-groups[[i]];id<-sprintf('S%02d',i)
 manifest[[i]]<-cbind(comparison_id=id,d[1,design,drop=FALSE])
 scenarios<-c('event','n_AC','n_BC','ess_ratio','method','ci_type')
 duplicate<-anyDuplicated(key(d,scenarios))>0
 if(duplicate) {
  audit[[i]]<-data.frame(comparison_id=id,status='skipped_duplicate_batches',
   message='Repeated design/method/CI rows. Select one batch with --batch; no seeds or runs were pooled.')
  next
 }
 folder<-file.path(out,id);dir.create(folder,showWarnings=FALSE)
 write.csv(d,file.path(folder,'comparison_metrics.csv'),row.names=FALSE,na='')
 context<-sprintf('%s | %s | adjustment: %s | fit: %s | calibration: %s',
  d$dataset[1],d$ci_type[1],d$estimation_mode[1],d$model_weighting[1],d$weight_mode[1])
 evgroups<-split(d,key(d,c('n_AC','n_BC','ess_ratio')))
 n_event<-0L
 for(gg in evgroups)if(all(c('common','rare') %in% gg$event)) {
  stem<-sprintf('common_rare_n%s_ess%03d',gg$n_AC[1],round(gg$ess_ratio[1]*100))
  save_figure(file.path(folder,stem),14,10,function()event_dashboard(gg,
   paste(context,'| n =',gg$n_AC[1],'| ESS ratio =',gg$ess_ratio[1])))
  n_event<-n_event+1L
 }
 varied<-length(unique(d$n_AC))>1 && length(unique(d$ess_ratio))>1
 trend_items<-if(all(d$ci_type=='Exact'))c('bias','rmse','coverage','returned_ci_fraction') else
  c('bias','rmse','accuracy','coverage','returned_ci_fraction')
 if(varied)for(metric in trend_items)
  save_figure(file.path(folder,paste0(metric,'_n_ess')),13,max(6,2.5*length(unique(d$ess_ratio))+1.5),
   function()trend_figure(d,metric,context))
 if(varied && all(is.finite(d$mean_realized_ess)))for(metric in c('rmse','coverage'))
  save_figure(file.path(folder,paste0(metric,'_realized_ess')),14,7.5,function()ess_figure(d,metric,context))
 audit[[i]]<-data.frame(comparison_id=id,status='generated',
  message=paste(n_event,'common/rare dashboards; n/ESS grid:',varied))
}
write.csv(do.call(rbind,manifest),file.path(out,'comparison_groups.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,audit),file.path(out,'comparison_audit.csv'),row.names=FALSE)
cfpath<-file.path(root,'conditional_metrics_all_runs.csv')
if(file.exists(cfpath)) {
 conditional<-read.csv(cfpath,stringsAsFactors=FALSE)
 folder<-file.path(out,'marginal_conditional');dir.create(folder,showWarnings=FALSE)
 target_audit<-list()
 for(rn in unique(metrics$run_name)) {
  mg<-metrics[metrics$run_name==rn & metrics$ci_type=='Wald',,drop=FALSE]
  cd<-conditional[conditional$run_name==rn & conditional$ci_type=='Wald',,drop=FALSE]
  if(!nrow(cd) || max(cd$n_evaluated)<100)next
  methods<-intersect(mg$method,cd$method)
  if(!length(methods))next
  mg<-mg[mg$method %in% methods,,drop=FALSE];cd<-cd[cd$method %in% methods,,drop=FALSE]
  ev<-mg$event[1];mg$event<-'marginal';cd$event<-'conditional'
  dd<-rbind(mg,cd)
  name<-gsub('[^a-zA-Z0-9_.-]','_',rn)
  write.csv(dd,file.path(folder,paste0(name,'.csv')),row.names=FALSE,na='')
  save_figure(file.path(folder,name),14,10,function()event_dashboard(dd,
   sprintf('%s | n = %s | ESS ratio = %s | Wald | matching methods only',ev,mg$n_AC[1],mg$ess_ratio[1]),
   levels=c('marginal','conditional'),labels=c('Marginal A/C','Conditional A/C at h_E=0'),
   heading='Performance for two distinct estimands',
   extra_note='Each estimate uses its own saved truth; conditional calibration uses SD(error).'))
  target_audit[[length(target_audit)+1L]]<-data.frame(run_name=rn,n_matching_methods=length(methods))
 }
 if(length(target_audit))write.csv(do.call(rbind,target_audit),file.path(folder,'target_comparison_index.csv'),row.names=FALSE)
}
writeLines(c('Only matched saved designs are compared. brmad, LP-BFC, RLP-BFC and smoke runs (n_expected < 100) excluded.',
 'Accuracy uses one grey hue. Coverage uses one blue Wald hue or one orange Exact hue; opacity increases with distance from target.',
 'Runs are never pooled; duplicate scenario rows cause a group to be skipped until --batch selects one batch.',
 'Design matching fixes DGM, dataset, adjustment, calibration, fitting weights, CI and Exact source MD5.',
 'common_rare: paired event comparison, with n and ESS fixed.',
 'n_ess: sample size on x; ESS in rows; event type in columns; shared y-axis; focus methods only.',
 'realized_ess: mean observed ESS on x; n fixed in columns; common/rare in rows.',
 'Bias and coverage use +/-1.96 Monte Carlo SE, not confidence intervals from a single fitted model.',
 'Coverage conditions on valid returned intervals. CI-return curves are supplied as a separate diagnostic.',
 'Adjustment effects are not automatically contrasted across designs with different fitting weights.',
 'Conditional effects use saved per-seed truths. marginal_conditional compares performance for two distinct targets, Wald only.',
 'See comparison_groups.csv and comparison_audit.csv for all grouping decisions.'),file.path(out,'README.txt'))
cat('Comparisons saved to:',out,'\n')
