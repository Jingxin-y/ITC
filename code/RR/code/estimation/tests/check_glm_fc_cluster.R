## Verify cluster checkpoint/export and fitting-weight forwarding; the CI spy
## verifies the call contract, not numerical interval accuracy.
source("load_anchored.R",encoding="UTF-8")
source_anchored()
source("estimation/cluster_all_methods.R",encoding="UTF-8")
check_cluster_fc <- function() {
  keys <- c("MAIC_EVENT","MAIC_EST_MODE","MAIC_MODEL_WEIGHTING","MAIC_WEIGHT_MODE","MAIC_DGP")
  old <- Sys.getenv(keys,unset=NA_character_)
  real_exact <- exact_marginal_target
  on.exit({
    exact_marginal_target <<- real_exact
    for (key in keys) if (is.na(old[key])) Sys.unsetenv(key) else do.call(Sys.setenv,setNames(list(old[key]),key))
  })
  Sys.setenv(MAIC_EVENT="common",MAIC_EST_MODE="none",MAIC_WEIGHT_MODE="all",MAIC_DGP="legacy")
  calls <- list()
  exact_marginal_target <<- function(...) {
    a <- list(...); calls[[length(calls)+1L]] <<- a
    point <- marginal_logrr_target(a$pars,a$va,a$vb,a$standardization.weight)
    list(estimate=point,low=point-.5,up=point+.5,p=.5,status="available",
      interval.available=TRUE,p.status="available",diagnostics=list(call_spy=TRUE))
  }
  for (mode in c("unit","maic")) {
    Sys.setenv(MAIC_MODEL_WEIGHTING=mode); calls <- list()
    folder <- tempfile("glm_fc_cluster_"); dir.create(folder)
    invisible(capture.output(z <- maic_run_seed(normalizePath("."),folder,811L,
      ci_mode="smoke",n=200L,ess_ratio=.5)))
    expected_design <- webster_rr_design(event="common",ess_ratio=.5)
    stopifnot(z$complete,nrow(z$summary)==18L,z$settings$version==10L,
      z$settings$counterfactual_risk_check=="finite_positive_standardized_means",
      identical(z$settings$dgm,expected_design$dgm),
      identical(z$settings$dgm_parameters,expected_design$parameters),
      identical(z$result$simulation$design$risk_table,expected_design$risk_table),
      all(z$summary$dgm==expected_design$dgm),
      all(z$coefficients$dgm==expected_design$dgm),
      z$settings$model_weighting==mode,
      all(anchored_glm_fc_method_names() %in% z$summary$method),length(calls)==3L)
    for (i in seq_along(calls)) {
      actual <- calls[[i]]$weight
      expected <- if (mode=="unit") rep(1,200) else z$result$model_weights$weight
      stopifnot(isTRUE(all.equal(actual,expected,check.attributes=FALSE)),
        isTRUE(all.equal(calls[[i]]$standardization.weight,z$result$weights$weight,check.attributes=FALSE)))
    }
    saved <- readRDS(file.path(folder,"by_seed","seed_000811.rds"))
    stopifnot(saved$complete,identical(saved$summary,z$summary))
    checkpoint <- file.path(folder,"by_seed","seed_000811.rds")
    ## Version 8 filtered Bernoulli FC targets; version 9 also filtered LB/RLB.
    ## Neither policy can silently resume as an unfiltered version 10 run.
    for (old_version in c(8L, 9L)) {
      old_record <- saved
      old_record$settings$version <- old_version
      old_record$settings$counterfactual_risk_check <-
        if (old_version == 9L) "strict_bernoulli_before_standardization" else NULL
      saveRDS(old_record,checkpoint)
      rejected <- tryCatch(maic_run_seed(normalizePath("."),folder,811L,
        ci_mode="smoke",n=200L,ess_ratio=.5),error=identity)
      stopifnot(inherits(rejected,"error"),
        grepl("Existing result has different settings",conditionMessage(rejected),fixed=TRUE))
    }
    saveRDS(saved,checkpoint)
    calls_before_resume <- length(calls)
    invisible(capture.output(resumed <- maic_run_seed(normalizePath("."),folder,811L,
      ci_mode="smoke",n=200L,ess_ratio=.5)))
    stopifnot(identical(resumed$summary,saved$summary),
      identical(resumed$settings,saved$settings),length(calls)==calls_before_resume)
    cat("PASS:",mode,"Webster DGM, cluster weights, FC exports, Exact forwarding and checkpoint compatibility\n")
  }
  ## Check the scenario grid without substituting the data generator or fitting.
  Sys.setenv(MAIC_MODEL_WEIGHTING="unit")
  for(event in c("common","rare")) for(ess in c(.15,.5,.8)) {
    Sys.setenv(MAIC_EVENT=event)
    folder <- tempfile("webster_cluster_grid_"); dir.create(folder)
    invisible(capture.output(z <- maic_run_seed(normalizePath("."),folder,811L,
      ci_mode="smoke",n=200L,ess_ratio=ess)))
    expected <- webster_rr_design(event=event,ess_ratio=ess)
    stopifnot(z$complete,identical(z$result$simulation$design$risk_table,expected$risk_table),
      identical(z$settings$dgm_parameters,expected$parameters),
      all(z$summary$dgm==expected$dgm),
      all(z$summary$population_log_RR_AC_in_BC[z$summary$component=="marginal"]==
        unname(expected$marginal_truth["BC","log_RR_AC"])))
  }
  cat("PASS: common/rare and ESS .15/.5/.8 cluster scenarios use Webster risks and truth\n")
}
check_cluster_fc()
