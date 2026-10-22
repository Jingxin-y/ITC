## Check new export/configuration and deterministic data sharing using a CI spy.
## Real bootstrap CI execution is tested by run_cluster_ess_em_smoke.R.
source("load_anchored.R", encoding = "UTF-8")
source_anchored()
source("estimation/cluster_all_methods.R", encoding = "UTF-8")
check_new_cluster_settings <- function() {
  keys <- c("MAIC_DGP", "MAIC_EVENT", "MAIC_ESS_SOURCE", "MAIC_EM_LEVEL", "MAIC_WEIGHT_MODE",
            "MAIC_MODEL_WEIGHTING", "MAIC_EST_MODE", "MAIC_EST_EM", "MAIC_EST_NON_EM")
  old <- Sys.getenv(keys, unset = NA_character_)
  real <- exact_marginal_target
  calls <- list()
  on.exit({
    exact_marginal_target <<- real
    for (key in keys) if (is.na(old[key])) Sys.unsetenv(key) else
      do.call(Sys.setenv, setNames(list(old[key]), key))
  })
  exact_marginal_target <<- function(...) {
    a <- list(...); calls[[length(calls) + 1L]] <<- a
    point <- marginal_logrr_target(a$pars, a$va, a$vb, a$standardization.weight)
    list(estimate = point, low = point - .5, up = point + .5, p = .5,
         status = "available", p.status = "available", interval.available = TRUE,
         diagnostics = list(call_spy = TRUE))
  }
  Sys.setenv(MAIC_DGP = "ess_em", MAIC_EVENT = "common", MAIC_ESS_SOURCE = "em_driven",
             MAIC_EM_LEVEL = "high", MAIC_MODEL_WEIGHTING = "unit", MAIC_EST_MODE = "scores",
             MAIC_EST_EM = "z1,z2,z3", MAIC_EST_NON_EM = "z4,z5,z6")
  data <- list()
  for (mode in c("all", "em", "nonem")) {
    Sys.setenv(MAIC_WEIGHT_MODE = mode); calls <- list()
    folder <- tempfile("ess_em_cluster_"); dir.create(folder)
    invisible(capture.output(z <- maic_run_seed(normalizePath("."), folder, 1L,
      "smoke", n = 200L, ess_ratio = .25)))
    data[[mode]] <- z$result$data
    expected <- webster_rr_scenario_design("common", .25, "em_driven", "high")
    stopifnot(z$complete, z$settings$version == 10L,
      z$settings$counterfactual_risk_check == "finite_positive_standardized_means",
      z$settings$cmh_weighting == "maic",
      z$settings$cmh_sandwich_uncertainty == "fixed",
      identical(z$result$simulation$design, expected),
      all(z$summary$ess_source == "em_driven"), all(z$summary$em_level == "high"),
      all(z$summary$weight_mode == mode), all(z$result$model_weights$weight == 1),
      length(calls) == 2L)
    for (a in calls) stopifnot(all(a$weight == 1),
      isTRUE(all.equal(a$standardization.weight, z$result$weights$weight)))
    record <- readRDS(file.path(folder, "by_seed", "seed_000001.rds"))
    Sys.setenv(MAIC_EM_LEVEL = "low")
    rejected <- try(maic_run_seed(normalizePath("."), folder, 1L, "smoke", n = 200L, ess_ratio = .25), silent = TRUE)
    stopifnot(inherits(rejected, "try-error"), grepl("different settings", as.character(rejected), fixed = TRUE))
    Sys.setenv(MAIC_EM_LEVEL = "high")
    before <- length(calls)
    invisible(capture.output(resumed <- maic_run_seed(normalizePath("."), folder, 1L, "smoke", n = 200L, ess_ratio = .25)))
    stopifnot(identical(resumed$summary, record$summary), length(calls) == before)
  }
  for (mode in c("em", "nonem")) for (key in c("data.IPD", "data.AD.or", "data.AD"))
    stopifnot(identical(data$all[[key]], data[[mode]][[key]]))
  cat("PASS: all/em/nonem share generated data; unit/MAIC CI weights forwarded; metadata and checkpoint isolation\n")
}
check_new_cluster_settings()
