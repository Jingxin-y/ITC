must<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
md5<-function(p)unname(tools::md5sum(p))
load.runtime<-function(home){
    m<-read.csv(file.path(home,'RUNTIME_FILES.csv'),stringsAsFactors=FALSE)
    must(all(md5(file.path(home,'code',m$path))==m$md5),'Runtime source checksum mismatch')
    code<-file.path(home,'code');rd<-file.path(code,'RD')
    source(file.path(rd,'load_anchored_rd.R'),encoding='UTF-8',local=.GlobalEnv)
    source_anchored_rd(code,rd)
    f<-formals(exact_marginal_rd_target);controls<-rd_full_bc_control()
    must(identical(as.integer(c(f$nsim.coarse,f$nsim.confirm,f$nsim.refine)),c(600L,1000L,1000L))&&
        isTRUE(f$optim.multistart)&&isTRUE(f$optim.fallback)&&isTRUE(f$optim.boundary.rescue)&&
        f$optim.multistart.max==3L,'Final RD defaults missing')
    txt<-gsub('[[:space:]]','',paste(deparse(body(exact_marginal_rd_target)),collapse=''))
    must(length(gregexpr('$acceptability>=(1-level)-1e-12',txt,fixed=TRUE)[[1]])==4L,'Final cutoff missing')
    must(identical(unname(unlist(controls[c('nsim.coarse','nsim.confirm','nsim.refine')])),c(600L,1000L,1000L))&&
        isTRUE(controls$optim.multistart)&&isTRUE(controls$optim.fallback)&&isTRUE(controls$optim.boundary.rescue)&&
        controls$nsim.p.value==500L&&controls$time.limit.sec==10800,'Formal full-mode CI controls differ')
    for(event in c('common','rare')){
        d<-webster_rd_scenario_design(event,.15,'diffuse','high')
        target<-if(event=='common')c(.418596226792393,.642197586315061,.502197586315061)else
            c(.0488309131792698,.1089376565366341,.0589376565366341)
        must(d$dgm=='webster2024_rd_ess_em_v2'&&d$scenario$version==2L&&
            max(abs(d$marginal_truth['BC',c('p_C','p_A','p_B')]-target))<1e-11,'Wrong RD v2 DGM or risk targets')
    }
    list(version=readLines(file.path(home,'VERSION.txt'),warn=FALSE),sources=setNames(m$md5,m$path),
         rd_sources=rd_source_hashes(code,rd),controls=controls)
}
tasks<-function(home){
    x<-read.csv(file.path(home,'tasks.csv'),stringsAsFactors=FALSE)
    must(nrow(x)==24L&&identical(x$task_id,1:24)&&all(x$nsim==1000L)&&all(x$first_seed==1L)&&
        all(x$ess_source=='diffuse')&&all(x$em_level=='high')&&all(x$weight_mode=='all')&&
        all(x$n %in% c(100L,500L))&&all(x$ess_ratio %in% c(.15,.5))&&
        all(table(x$scenario_id)==3L)&&!anyDuplicated(x$cell),'Invalid 24-cell design')
    x
}
settings.for<-function(r,runtime){
    d<-webster_rd_scenario_design(r$event,r$ess_ratio,'diffuse','high')
    analysis<-list(run_exact=TRUE,model_weighting=r$model_weighting,weight_covariates=paste0('z',1:6),
        score_spec=NULL,score_standardization='bc_estimated',model_covariates=NULL,effect_modifiers=NULL,
        fit_control=anchored_fit_control(list()),bc_control=runtime$controls,exact_seed_offset=1000000L,
        estimation_mode=r$estimation_mode,brm_effect_model='linear_atanh')
    list(version=1L,effect_scale='RD',target='BC marginal A-C RD',methods=rd_method_plan(),
        probability_mapping=.anchored_rd_solver_env$probability_mapping,dgm=d$dgm,design=d,
        n_AC=r$n,n_BC=r$n,analysis=analysis,source_md5=runtime$rd_sources)
}
config<-function(r,runtime)list(version=runtime$version,task=r,source_md5=runtime$sources,
    ci_control=runtime$controls,data_seeds=1:1000,ci_seed_rule='1000000 + data seed',
    dgm='webster2024_rd_ess_em_v2',methods=rd_method_plan(),rd_settings=settings.for(r,runtime))
case.dir<-function(out,r)file.path(out,'cells',r$cell)
seed.path<-function(d,seed)file.path(d,'by_seed',sprintf('rd_seed_%08d.rds',seed))
check.record<-function(z,seed,expected){
    s<-z$settings;s$requested_seeds<-NULL
    # Named list order differs after removing requested_seeds in the canonical producer.
    must(inherits(z,'anchored_rd_simulation')&&identical(z$effect_scale,'RD')&&
        setequal(names(s),names(expected))&&identical(s[names(expected)],expected),
        'Seed settings, effect scale or source hashes differ')
    must(identical(z$requested_seeds,as.integer(seed)),'Seed requested identifier differs')
    if(isTRUE(z$complete))must(identical(z$seeds,as.integer(seed))&&nrow(z$estimates)==9L&&
        identical(z$estimates$method,rd_method_plan()$method)&&all(z$estimates$r==seed),'Complete RD record schema differs')
    invisible(TRUE)
}
resources<-function(home,r){
    tab<-read.csv(file.path(home,'resources.csv'),stringsAsFactors=FALSE)
    z<-tab[tab$analysis==r$analysis&tab$n==r$n,,drop=FALSE];must(nrow(z)==1L,'Resource setting missing')
    prefix<-switch(r$analysis,none='NONE',score_unit='SCORE',score_MAIC='MAIC')
    z$cpus<-as.integer(Sys.getenv(paste0('RDFULL_',prefix,'_CPUS'),as.character(z$cpus)))
    z$memory<-Sys.getenv(paste0('RDFULL_',prefix,'_MEM'),z$memory)
    z$walltime<-Sys.getenv(if(prefix=='NONE')'RDFULL_NONE_TIME'else paste0('RDFULL_',prefix,'_TIME_N',r$n),z$walltime)
    must(is.finite(z$cpus)&&z$cpus>=1L&&z$cpus<=192L&&grepl('^[1-9][0-9]*[MG]$',z$memory)&&
        grepl('^([0-9]+-)?[0-9]+:[0-9]{2}:[0-9]{2}$',z$walltime),'Invalid requested resource setting')
    z
}
check.config<-function(d,cfg){
    p<-file.path(d,'configuration.rds');must(file.exists(p)&&identical(readRDS(p),cfg),'Incompatible formal run configuration')
}
write.csv.atomic<-function(x,p){
    tmp<-paste0(p,'.tmp_',Sys.getpid());write.csv(x,tmp,row.names=FALSE,na='')
    if(.Platform$OS.type=='windows'&&file.exists(p))unlink(p)
    must(file.rename(tmp,p),paste('Cannot finalize',p))
}
collect.cell<-function(d,r,runtime){
    expected<-settings.for(r,runtime);records<-list();missing<-incomplete<-integer();issues<-list()
    for(seed in 1:1000){
        p<-seed.path(d,seed)
        if(!file.exists(p)){missing<-c(missing,seed);next}
        z<-tryCatch(readRDS(p),error=function(e)e)
        if(inherits(z,'error')){
            incomplete<-c(incomplete,seed);issues[[length(issues)+1L]]<-data.frame(seed=seed,message=conditionMessage(z));next
        }
        check.record(z,seed,expected)
        if(!isTRUE(z$complete)){incomplete<-c(incomplete,seed);next}
        records[[length(records)+1L]]<-z
    }
    write.csv.atomic(data.frame(seed=missing),file.path(d,'missing_seeds.csv'))
    write.csv.atomic(data.frame(seed=incomplete),file.path(d,'incomplete_seeds.csv'))
    write.csv.atomic(if(length(issues))do.call(rbind,issues)else data.frame(seed=integer(),message=character()),
                     file.path(d,'record_issues.csv'))
    if(length(records)){
        combined<-records[[1]]
        combined$estimates<-do.call(rbind,lapply(records,function(z)z$estimates))
        combined$coefficient_estimates<-do.call(rbind,lapply(records,function(z)z$coefficient_estimates))
        combined$reference_truth<-do.call(rbind,lapply(records,function(z)z$reference_truth))
        combined$diagnostics<-unlist(lapply(records,function(z)z$diagnostics),recursive=FALSE)
        combined$seeds<-combined$requested_seeds<-unlist(lapply(records,function(z)z$seeds))
        combined$settings$requested_seeds<-combined$seeds
        combined$parameter_schema<-unique(combined$estimates[c('method','component','parameter','term')])
        rd_write_simulation(combined,d)
    }
    data.frame(records=1000L-length(missing),processed=length(records),incomplete=length(incomplete),missing=length(missing))
}
