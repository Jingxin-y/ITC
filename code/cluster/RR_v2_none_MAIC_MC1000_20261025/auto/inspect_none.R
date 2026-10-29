inspect.none<-function(home,out,state,round,seeds=1:1000,max.rounds=10L){
    source(file.path(home,'common.R'),encoding='UTF-8',local=.GlobalEnv)
    runtime<-load.runtime(home);tab<-tasks(home);tab<-tab[tab$analysis=='none_MAIC',,drop=FALSE]
    must(identical(tab$task_id,1:8)&&all(tab$model_weighting=='maic'),'Unexpected none cases')
    expected<-runtime$sources[grepl('\\.[Rr]$',names(runtime$sources))]
    dir.create(state,recursive=TRUE,showWarnings=FALSE);rows<-details<-list()
    for(i in seq_len(nrow(tab))){
        r<-tab[i,,drop=FALSE];d<-case.dir(out,r);check.config(d,config(r,runtime))
        count<-c(complete=0L,missing=0L,partial=0L,invalid=0L)
        for(seed in seeds){
            p<-file.path(d,'by_seed',sprintf('seed_%06d.rds',seed));error<-''
            status<-if(!file.exists(p))'missing'else tryCatch({
                z<-readRDS(p);s<-z$settings
                must(identical(z$source_md5[names(expected)],expected),'Scientific source checksum differs')
                must(isTRUE(s$seed==seed)&&isTRUE(s$exact_seed==1000000L+seed)&&
                    identical(s$dgm,'webster2024_rr_ess_em_v2')&&identical(s$ci_mode,'full')&&
                    identical(s$event,r$event)&&isTRUE(s$n_AC==r$n)&&isTRUE(s$n_BC==r$n)&&
                    isTRUE(all.equal(s$ess_ratio,r$ess_ratio))&&identical(s$estimation_mode,'none')&&
                    identical(s$model_weighting,'maic')&&identical(s$weight_covariates,paste0('z',1:6))&&
                    identical(s$ci_control,runtime$controls),'Seed identity or analysis settings differ')
                must(is.data.frame(z$summary)&&nrow(z$summary)>0L&&'seed_complete'%in%names(z$summary),'Invalid seed summary')
                if(isTRUE(z$complete)){
                    must(all(z$summary$seed_complete)&&!any(z$summary$status=='pending'),'Conflicting complete flags')
                    'complete'
                }else 'partial'
            },error=function(e){error<<-conditionMessage(e);'invalid'})
            count[status]<-count[status]+1L
            if(status!='complete')details[[length(details)+1L]]<-data.frame(task_id=r$task_id,cell=r$cell,seed=seed,status=status,error=error)
        }
        rows[[i]]<-data.frame(task_id=r$task_id,cell=r$cell,planned=length(seeds),
            complete=count['complete'],missing=count['missing'],partial=count['partial'],invalid=count['invalid'],row.names=NULL)
        cat('NONE CHECK',r$task_id,'complete',count['complete'],'missing',count['missing'],'partial',count['partial'],'invalid',count['invalid'],'\n');flush.console()
    }
    summary<-do.call(rbind,rows);detail<-if(length(details))do.call(rbind,details)else
        data.frame(task_id=integer(),cell=character(),seed=integer(),status=character(),error=character())
    write.csv.atomic(summary,file.path(state,'latest_summary.csv'));write.csv.atomic(detail,file.path(state,'pending_seeds.csv'))
    write.csv.atomic(summary,file.path(state,sprintf('check_%03d_summary.csv',round)))
    previous.path<-file.path(state,'progress.rds');previous<-if(file.exists(previous.path))readRDS(previous.path)else NULL
    completed<-sum(summary$complete);stalled<-0L
    if(!is.null(previous)&&round>previous$round&&completed<=previous$complete)stalled<-previous$stalled+1L
    saveRDS(list(round=round,complete=completed,stalled=stalled),previous.path,version=2)
    action<-if(sum(summary$invalid)>0L)'STOP_INVALID_RECORDS'else if(all(summary$complete==length(seeds)))'COMPLETE'else
        if(stalled>=2L)'STOP_NO_PROGRESS'else if(round>=max.rounds)'STOP_MAX_ROUNDS'else 'RESUME'
    todo<-summary$task_id[summary$missing+summary$partial>0L]
    writeLines(paste(todo,collapse=','),file.path(state,'pending_tasks.txt'));writeLines(action,file.path(state,'action.txt'))
    cat('NONE ACTION:',action,'complete',completed,'/',length(seeds)*8L,'round',round,'\n')
    invisible(list(action=action,summary=summary,detail=detail))
}
if(sys.nframe()==0L){
    a<-commandArgs(TRUE);stopifnot(length(a)==5L)
    inspect.none(normalizePath(a[1],winslash='/',mustWork=TRUE),a[2],a[3],as.integer(a[4]),max.rounds=as.integer(a[5]))
}
