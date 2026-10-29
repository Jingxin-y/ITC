a<-commandArgs(TRUE);stopifnot(length(a)==4L)
home<-normalizePath(a[1],winslash='/',mustWork=TRUE);out<-a[2];id<-as.integer(a[3]);cores<-as.integer(a[4])
source(file.path(home,'common.R'),encoding='UTF-8');runtime<-load.runtime(home);tab<-tasks(home)
must(id %in% 1:24&&cores>=1L&&cores<=192L,'Invalid task or CPU count')
r<-tab[id,,drop=FALSE];d<-case.dir(out,r);cfg<-config(r,runtime);check.config(d,cfg)
dir.create(file.path(d,'logs'),showWarnings=FALSE)
capture.output(sessionInfo(),file=file.path(d,'logs',paste0('sessionInfo_',Sys.getenv('SLURM_JOB_ID','local'),'.txt')))
cat('FORMAL START',r$cell,'1000 seeds;',cores,'workers; RD v2; final CI 600/1000/1000\n');flush.console()
worker<-function(seed){
    started<-proc.time()[['elapsed']]
    error<-tryCatch({
        p<-seed.path(d,seed)
        if(file.exists(p)){
            old<-readRDS(p);check.record(old,seed,cfg$rd_settings)
            if(!isTRUE(old$complete)){
                backup<-file.path(d,'partial_history',paste0(basename(p),'.',md5(p)))
                dir.create(dirname(backup),showWarnings=FALSE)
                if(!file.exists(backup))must(file.copy(p,backup),'Cannot preserve partial checkpoint')
            }
        }
        rd_run_seed(code_dir=file.path(home,'code'),rd_dir=file.path(home,'code/RD'),run_dir=d,seed=seed,
            n=r$n,ess_ratio=r$ess_ratio,event=r$event,ess_source='diffuse',em_level='high',weight_mode='all',
            model_weighting=r$model_weighting,ci_mode='full',bc_control=runtime$controls,
            estimation_mode=r$estimation_mode,brm_effect_model='linear_atanh')
        z<-readRDS(p);check.record(z,seed,cfg$rd_settings);must(isTRUE(z$complete),'Pipeline record incomplete');NULL
    },error=function(e)conditionMessage(e))
    cat('FORMAL SEED',r$cell,seed,if(is.null(error))'PROCESSED'else error,
        'elapsed',proc.time()[['elapsed']]-started,'\n');flush.console()
    list(seed=seed,error=error)
}
# Linux fork workers share loaded code; each seed keeps the canonical data/CI RNG rules.
answer<-parallel::mclapply(1:1000,worker,mc.cores=min(cores,1000L),mc.preschedule=FALSE,mc.set.seed=FALSE)
bad<-vapply(answer,function(x)!is.list(x)||!is.null(x$error),logical(1))
write.csv.atomic(data.frame(seed=which(bad),message=vapply(answer[bad],function(x)
    if(is.list(x))x$error else 'Worker terminated without a result',character(1))),file.path(d,'logs/seed_errors.csv'))
progress<-collect.cell(d,r,runtime)
if(any(bad)||progress$processed!=1000L){cat('INCOMPLETE:',r$cell,'processed',progress$processed,'of 1000\n');quit(status=2L)}
writeLines(md5(file.path(d,'configuration.rds')),file.path(d,'CASE_PROCESSED.txt'))
cat('FORMAL FINISHED:',r$cell,'1000 records processed. Method/interval failures remain in result rows.\n')
