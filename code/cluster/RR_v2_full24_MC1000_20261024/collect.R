a<-commandArgs(TRUE);stopifnot(length(a)==2L)
home<-normalizePath(a[1],winslash='/',mustWork=TRUE);out<-a[2]
source(file.path(home,'common.R'),encoding='UTF-8');runtime<-load.runtime(home);tab<-tasks(home)
report<-file.path(out,'report');dir.create(report,showWarnings=FALSE)
progress<-metrics<-list();combined<-file.path(report,'combined_results.csv');tmp<-paste0(combined,'.tmp_',Sys.getpid());first<-TRUE
for(i in 1:24){
    r<-tab[i,,drop=FALSE];d<-case.dir(out,r)
    if(!file.exists(file.path(d,'configuration.rds'))){
        progress[[i]]<-data.frame(task_id=i,cell=r$cell,records=0L,processed=0L,incomplete=0L,missing=1000L);next
    }
    check.config(d,config(r,runtime));processed<-maic_collect(d)
    missing<-nrow(read.csv(file.path(d,'missing_seeds.csv')));incomplete<-nrow(read.csv(file.path(d,'incomplete_seeds.csv')))
    progress[[i]]<-data.frame(task_id=i,cell=r$cell,records=1000L-missing,processed=processed,incomplete=incomplete,missing=missing)
    p<-file.path(d,'combined_results.csv');if(!file.exists(p))next
    x<-read.csv(p,stringsAsFactors=FALSE);x$task_id<-i;x$scenario_id<-r$scenario_id;x$analysis<-r$analysis
    write.table(x,tmp,sep=',',row.names=FALSE,col.names=first,append=!first,na='');first<-FALSE
    truth<-webster_rr_scenario_design(r$event,r$ess_ratio,'diffuse','high')$marginal_truth['BC','log_RR_AC']
    for(key in unique(paste(x$method,x$ci_type,sep='|'))){
        z<-x[paste(x$method,x$ci_type,sep='|')==key,,drop=FALSE]
        point<-is.finite(z$estimate)
        valid<-is.finite(z$low)&is.finite(z$up)&z$low<=z$up&z$status %in% c('available','boundary_truncated')
        covered<-valid&z$low<=truth&truth<=z$up;n<-sum(valid);pcover<-if(n)sum(covered)/n else NA_real_
        metrics[[length(metrics)+1L]]<-data.frame(task_id=i,scenario_id=r$scenario_id,event=r$event,n=r$n,
            ess_ratio=r$ess_ratio,analysis=r$analysis,method=z$method[1],ci_type=z$ci_type[1],planned=1000L,
            records=nrow(z),finite_point=sum(point),valid_interval=n,covered=sum(covered),
            coverage_valid_pct=100*pcover,coverage_MCSE_pp=if(n)100*sqrt(pcover*(1-pcover)/n)else NA_real_,
            returned_and_covered_pct=100*sum(covered)/1000,
            mean_logRR_width=if(n)mean(z$up[valid]-z$low[valid])else NA_real_,
            bias_finite_point=if(any(point))mean(z$estimate[point]-truth)else NA_real_,
            boundary_intervals=sum(valid&z$status=='boundary_truncated'))
    }
}
if(!first){if(.Platform$OS.type=='windows'&&file.exists(combined))unlink(combined);must(file.rename(tmp,combined),'Cannot finalize combined results')}
progress<-do.call(rbind,progress);write.csv.atomic(progress,file.path(report,'progress.csv'))
if(length(metrics))write.csv.atomic(do.call(rbind,metrics),file.path(report,'all_methods_summary.csv'))
complete<-all(progress$processed==1000L)
writeLines(c(if(complete)'ALL 24 CELLS PROCESSED: 1000 records each.'else 'INCOMPLETE: see progress.csv; resume with submit.sh.',
    'Processed records can include method or interval failures; valid_interval and coverage denominators are reported separately.',
    'Boundary truncation is a valid method outcome. No failed bootstrap responses are dropped or redrawn.',
    'All cells use the v2 DGM, diffuse/high/all, final multistart/nlminb solver and CI budget 600/1000/1000.'),file.path(report,'SUMMARY.txt'))
print(progress,row.names=FALSE);if(!complete)quit(status=2L)
