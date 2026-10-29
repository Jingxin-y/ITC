a<-commandArgs(TRUE);stopifnot(length(a)==3L)
home<-normalizePath(a[1],winslash='/',mustWork=TRUE);out<-a[2];requested<-a[3]
source(file.path(home,'common.R'),encoding='UTF-8');runtime<-load.runtime(home);tab<-tasks(home)
ids<-if(requested=='all')1:8 else as.integer(strsplit(requested,',',fixed=TRUE)[[1]])
must(length(ids)>0&&!anyNA(ids)&&!anyDuplicated(ids)&&all(ids %in% 1:8),'Task IDs must be all or unique values 1..8')
dir.create(out,recursive=TRUE,showWarnings=FALSE);plan<-list()
for(id in ids){
    r<-tab[id,,drop=FALSE];res<-resources(home,r);d<-case.dir(out,r);cfg<-config(r,runtime)
    dir.create(file.path(d,'by_seed'),recursive=TRUE,showWarnings=FALSE)
    p<-file.path(d,'configuration.rds')
    if(file.exists(p))check.config(d,cfg)else{
        must(!length(list.files(file.path(d,'by_seed'),'\\.rds$')),'Records exist without pinned configuration')
        saveRDS(cfg,p,version=2)
    }
    values<-c(nsim=1000,first_seed=1,ci_mode='full',event=r$event,estimation_mode=r$estimation_mode,
        n=r$n,ess_ratio=r$ess_ratio,cores=res$cpus,dgp_mode='ess_em',ess_source='diffuse',em_level='high',
        weight_mode='all',model_weighting=r$model_weighting)
    writeLines(paste(names(values),values,sep='='),file.path(d,'submission.txt'))
    marker<-file.path(d,'CASE_PROCESSED.txt')
    done<-file.exists(marker)&&identical(readLines(marker,warn=FALSE),md5(p))&&
        length(list.files(file.path(d,'by_seed'),'^seed_[0-9]{6}\\.rds$'))==1000L
    plan[[length(plan)+1L]]<-data.frame(task_id=id,scenario_id=r$scenario_id,analysis=r$analysis,n=r$n,
        event=r$event,ess_ratio=r$ess_ratio,nsim=1000,cpus=res$cpus,memory=res$memory,walltime=res$walltime,
        already_processed=done,to_submit=!done,cell=r$cell)
}
plan<-do.call(rbind,plan);write.csv.atomic(plan,file.path(out,'submission_plan.csv'));print(plan,row.names=FALSE)
todo<-plan[plan$to_submit,,drop=FALSE];groups<-unique(todo[c('cpus','memory','walltime')]);rows<-list()
for(i in seq_len(nrow(groups))){
    g<-groups[i,,drop=FALSE];keep<-todo$cpus==g$cpus&todo$memory==g$memory&todo$walltime==g$walltime
    rows[[i]]<-data.frame(group_id=i,cpus=g$cpus,memory=g$memory,walltime=g$walltime,
        tasks=paste(todo$task_id[keep],collapse=','))
}
groups<-if(length(rows))do.call(rbind,rows)else data.frame(group_id=integer(),cpus=integer(),memory=character(),walltime=character(),tasks=character())
write.table(groups,file.path(out,'submission_groups.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
cat('v2 DGM and final CI verified. Configuration only; no simulation or bootstrap executed.\n')
