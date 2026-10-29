args <- commandArgs(trailingOnly = TRUE)
own <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
rd_dir <- dirname(dirname(normalizePath(own, winslash = "/")))
code_dir <- if (length(args)) args[1] else dirname(rd_dir)
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source(file.path(rd_dir, "load_anchored_rd.R"), encoding = "UTF-8")
source_anchored_rd(code_dir, rd_dir)
check <- function(x, message) { if (!isTRUE(x)) stop(message); cat("PASS:", message, "\n") }
fails <- function(x) inherits(tryCatch(x, error = identity), "error")
out <- Sys.getenv("RD_TEST_OUTPUT", tempfile("rd_effect_models_"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
design <- webster_rd_scenario_design("common", .5, "diffuse", "high")
linear <- run.anchored.rd(design, 500, seed = 19)
saturated <- run.anchored.rd(design, 500, seed = 19, brm_effect_model = "saturated")
none <- run.anchored.rd(design, 500, seed = 19, estimation_mode = "none")
check(identical(linear$variables$brm_effect_model, "linear_atanh") &&
  all(vapply(linear$fits[c("brm", "brm-FC")], function(f)
    isTRUE(f$converged) && ncol(f$va) == 2L && ncol(f$vb) == 4L &&
    identical(colnames(f$va), c("(Intercept)", "h_E")), logical(1))),
  "default BRM/FC use two linear atanh-RD effects and four saturated OP parameters")
check(identical(saturated$variables$brm_effect_model, "saturated") &&
  ncol(saturated$fits$brm$va) == 4L && ncol(saturated$fits$brm$vb) == 4L &&
  identical(none$variables$brm_effect_model, "intercept_only") && ncol(none$fits$brm$va) == 1L,
  "explicit saturated effects and none intercept-only mode remain available")
raw <- c(paste0("z",1:6),"x","y")
other <- c("IB","RIB","RLPM","MN_eff","RD_MAIC_sandwich")
fields <- c("method","estimate","se","lower","upper","p.value","status","interval_status")
check(identical(linear$data$data.IPD$data[raw], saturated$data$data.IPD$data[raw]) &&
  identical(linear$data$data.AD, saturated$data$data.AD) &&
  identical(linear$weights, saturated$weights) &&
  identical(linear$estimates[linear$estimates$method %in% other,fields],
            saturated$estimates[saturated$estimates$method %in% other,fields]),
  "effect model changes preserve generated data, BC summaries, MAIC and all other methods")
for (method in c("brm","brm-FC")) {
  f <- linear$fits[[method]]; target <- marginal_rd_target(f$point.est,f$va,f$vb,linear$weights$weight,TRUE)
  independent <- weighted.mean(tanh(drop(f$va %*% f$point.est[1:2])),linear$weights$weight)
  numeric_gradient <- vapply(seq_along(f$point.est),function(j) {
    up <- down <- f$point.est; up[j] <- up[j]+1e-5; down[j] <- down[j]-1e-5
    (marginal_rd_target(up,f$va,f$vb,linear$weights$weight)-
      marginal_rd_target(down,f$va,f$vb,linear$weights$weight))/2e-5
  },numeric(1))
  check(abs(target$value-independent)<1e-12 && max(abs(target$gradient-numeric_gradient))<1e-8,
    paste(method,"marginal RD and full gradient match independent calculation"))
}
records <- list()
for (strength in c("low","high")) for (seed in c(11L,30L,33L,37L,58L)) {
  difficult <- webster_rd_scenario_design("rare",.15,"diffuse",strength)
  a <- run.anchored.rd(difficult,100,seed=seed)
  b <- run.anchored.rd(difficult,100,seed=seed,brm_effect_model="saturated")
  good <- a$estimates$method %in% c("brm","brm-FC")
  bad <- b$estimates$method %in% c("brm","brm-FC")
  check(all(a$estimates$status[good]=="available") && all(is.finite(a$estimates$estimate[good])) &&
    all(b$estimates$status[bad]=="fit_failed"),paste(strength,"seed",seed,"missing-arm cell recovered by linear effects"))
  records[[length(records)+1L]] <- data.frame(strength=strength,seed=seed,
    linear_BRM=a$estimates$estimate[a$estimates$method=="brm"],
    linear_FC=a$estimates$estimate[a$estimates$method=="brm-FC"])
}
write.csv(do.call(rbind,records),file.path(out,"missing_cell_comparison.csv"),row.names=FALSE)
failed <- run.anchored.rd(webster_rd_scenario_design("rare",.15,"diffuse","high"),100,seed=4)
check(all(failed$estimates$status=="weighting_failed"),"infeasible MAIC weights retain independent weighting-failed status")
check(fails(run.anchored.rd(design,500,seed=19,brm_effect_model="linear_rd")),
  "unsupported natural-linear RD option rejected")
path <- file.path(out,"checkpoint.rds")
sim <- simulate_anchored_rd(design,500,seeds=19L,checkpoint_file=path,code_dir=code_dir,rd_dir=rd_dir)
resumed <- simulate_anchored_rd(design,500,seeds=19L,checkpoint_file=path,code_dir=code_dir,rd_dir=rd_dir)
check(identical(sim$estimates,resumed$estimates) &&
  identical(sim$settings$analysis$brm_effect_model,"linear_atanh") &&
  all(sim$estimates$brm_effect_model=="linear_atanh") &&
  identical(sim$diagnostics[[1]]$brm_effect_model,"linear_atanh") &&
  fails(simulate_anchored_rd(design,500,seeds=19L,checkpoint_file=path,
    brm_effect_model="saturated",code_dir=code_dir,rd_dir=rd_dir)),
  "checkpoint records effect model, resumes and rejects saturated mixing")
parallel <- rd_run_scenario(code_dir,rd_dir,file.path(out,"parallel_linear"),n=100,
  ess_ratio=.15,seeds=c(11L,30L),cores=2L,event="rare",ci_mode="wald")
check(isTRUE(parallel$complete) && all(parallel$estimates$brm_effect_model=="linear_atanh") &&
  all(parallel$estimates$status[parallel$estimates$method %in% c("brm","brm-FC")]=="available"),
  "two-worker cluster forwards linear effects and collects recovered missing-cell point estimates")
one <- file.path(out,"mixed","by_seed"); dir.create(one,recursive=TRUE,showWarnings=FALSE)
saveRDS(sim,file.path(one,"rd_seed_00000019.rds"))
sat_sim <- simulate_anchored_rd(design,500,seeds=20L,brm_effect_model="saturated",code_dir=code_dir,rd_dir=rd_dir)
saveRDS(sat_sim,file.path(one,"rd_seed_00000020.rds"))
check(fails(rd_collect(dirname(one))),"collector rejects mixed linear and saturated models")
smoke <- rd_run_scenario(code_dir,rd_dir,file.path(out,"smoke_linear"),n=100,
  ess_ratio=.15,seeds=11L,cores=1L,event="rare",ci_mode="smoke")
bc <- smoke$estimates$ci_type=="BC"
check(all(is.finite(smoke$estimates$estimate[bc])) &&
  all(smoke$estimates$interval_status[bc] %in% c("available","BC_boundary_truncated")) &&
  all(vapply(smoke$diagnostics[[1]]$bc,function(z)
    identical(z$diagnostics$settings$optimizer.version,"rd_bc_multistart_nlminb_v1"),logical(1))),
  "real small-budget linear-effect BC retains existing optimizer rescue and natural RD intervals")
cat("RD EFFECT MODEL OUTPUT:",out,"\n")
