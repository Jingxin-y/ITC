args <- commandArgs(TRUE)
root <- normalizePath(if(length(args))args[1] else '..',winslash='/',mustWork=TRUE)
cfg <- readRDS(file.path(root,'original result/STC_RR/configuration.rds'))
metadata <- readRDS(file.path(root,'original result/STC_RR/report/run_metadata.rds'))
stopifnot(identical(cfg$version,'STC_RR_reuse_20261026_v1'),
 identical(as.integer(cfg$scenarios),1:8),identical(as.integer(cfg$seeds),1:1000),
 cfg$planned_per_scenario==1000,cfg$level==.95)
hash <- tools::md5sum(file.path(root,'code/STC_RR',names(cfg$runtime_hash)))
stopifnot(identical(unname(hash),unname(cfg$runtime_hash)),identical(metadata$configuration,cfg))
cat('PASS: STC runtime MD5 matches saved configuration; 8 scenarios, seeds 1..1000, level .95.\n')
