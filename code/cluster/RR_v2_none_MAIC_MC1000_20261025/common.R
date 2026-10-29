must<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
md5<-function(p)unname(tools::md5sum(p))
load.runtime<-function(home){
    m<-read.csv(file.path(home,'RUNTIME_FILES.csv'),stringsAsFactors=FALSE)
    must(all(md5(file.path(home,'code',m$path))==m$md5),'Runtime source checksum mismatch')
    code<-file.path(home,'code');source(file.path(code,'load_anchored.R'),encoding='UTF-8',local=.GlobalEnv)
    source_anchored(code);source(file.path(code,'estimation/cluster_all_methods.R'),encoding='UTF-8',local=.GlobalEnv)
    f<-formals(exact_marginal_target)
    must(identical(as.integer(c(f$nsim.coarse,f$nsim.confirm,f$nsim.refine)),c(600L,1000L,1000L))&&
        isTRUE(f$optim.multistart)&&isTRUE(f$optim.fallback)&&f$optim.multistart.max==3L,'Final RR defaults missing')
    body.text<-gsub('[[:space:]]','',paste(deparse(body(exact_marginal_target)),collapse=''))
    must(length(gregexpr('$acceptability>=(1-level)-1e-12',body.text,fixed=TRUE)[[1]])==4L,'Final cutoff missing')
    z<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('controls')),as.list(body(maic_run_seed)))
    must(length(z)==1L,'Full CI control definition missing')
    controls<-eval(z[[1]][[3]],list(ci_mode='full'))
    must(identical(unname(unlist(controls[c('nsim.coarse','nsim.confirm','nsim.refine')])),c(600L,1000L,1000L))&&
        isTRUE(controls$optim.multistart)&&isTRUE(controls$optim.fallback)&&controls$nsim.p.value==500L,
        'Formal full-mode CI controls differ')
    for(event in c('common','rare')){
        d<-webster_rr_scenario_design(event,.15,'diffuse','high')
        target<-if(event=='common')c(.418596226792393,.514669087904676,.60)else c(.0488309131792698,.10,.15)
        must(d$dgm=='webster2024_rr_ess_em_v2'&&d$scenario$version==2L&&
            max(abs(d$marginal_truth['BC',c('p_C','p_A','p_B')]-target))<1e-11,'Wrong RR DGM or risk targets')
    }
    list(version=readLines(file.path(home,'VERSION.txt'),warn=FALSE),sources=setNames(m$md5,m$path),controls=controls)
}
tasks<-function(home){
    x<-read.csv(file.path(home,'tasks.csv'),stringsAsFactors=FALSE)
    must(nrow(x)==8L&&identical(x$task_id,1:8)&&identical(x$scenario_id,1:8)&&
        all(x$nsim==1000L)&&all(x$first_seed==1L)&&all(x$analysis=='none_MAIC')&&
        all(x$estimation_mode=='none')&&all(x$model_weighting=='maic')&&
        all(x$ess_source=='diffuse')&&all(x$em_level=='high')&&all(x$weight_mode=='all')&&
        identical(x$event,rep(c('common','rare'),each=4L))&&
        identical(x$n,rep(c(100L,100L,500L,500L),2L))&&
        isTRUE(all.equal(x$ess_ratio,rep(c(.15,.5),4L)))&&!anyDuplicated(x$cell),
        'Invalid eight-case none+MAIC design')
    x
}
set.analysis<-function(r){
    Sys.setenv(MAIC_DGP='ess_em',MAIC_EVENT=r$event,MAIC_ESS_SOURCE='diffuse',MAIC_EM_LEVEL='high',
        MAIC_WEIGHT_MODE='all',MAIC_EST_MODE=r$estimation_mode,MAIC_MODEL_WEIGHTING=r$model_weighting,
        MAIC_EST_EM='z1,z2,z3',MAIC_EST_NON_EM='z4,z5,z6')
}
config<-function(r,runtime)list(version=runtime$version,task=r,source_md5=runtime$sources,
    ci_control=runtime$controls,data_seeds=1:1000,ci_seed_rule='1000000 + data seed',
    dgm='webster2024_rr_ess_em_v2',risk_target_version='BC_CAB_20261015',methods='original formal all-method plan')
case.dir<-function(out,r)file.path(out,'cells',r$cell)
resources<-function(home,r){
    z<-read.csv(file.path(home,'resources.csv'),stringsAsFactors=FALSE)
    z<-z[z$analysis==r$analysis&z$n==r$n,,drop=FALSE]
    must(nrow(z)==1L&&z$cpus==64L&&z$memory=='10G'&&z$walltime=='04:00:00','Expected 64 CPUs / 10G / 4h')
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
