## Run from anchored_v1: Rscript estimation/tests/check_webster_pipeline.R
source("load_anchored.R",encoding="UTF-8")
source_anchored()
close_to <- function(a,b,tol=1e-9) stopifnot(isTRUE(all.equal(a,b,tolerance=tol)))
d <- webster_rr_design()
scores <- make_brm_score_spec(c("z1","z2","z3"),c("z4","z5","z6"))
set.seed(811)
sim <- generate_anchored_data(1000,1200,d)
resolved <- brm_scores_from_bc(scores,sim$data.AD)
variables <- list(score_spec=resolved,score_standardization="bc_estimated",
                  effect_modifiers="h_E",bc_target="marginal")
truth <- brm_score_coefficient_truth(d,variables)
expected <- c(log(d$parameters$rr_AC)+unname(resolved$center["h_E"])*log(d$parameters$em_rr_AC),
              unname(resolved$scale["h_E"])*log(d$parameters$em_rr_AC))
brm <- subset(truth,method=="brm")
close_to(brm$truth[match(c("(Intercept)","h_E"),brm$term)],expected)
wrong <- make_brm_score_spec(c("z1","z4"),c("z2","z3","z5","z6"))
variables$score_spec <- brm_scores_from_bc(wrong,sim$data.AD)
stopifnot(is.null(brm_score_coefficient_truth(d,variables)))
cat("PASS: per-replicate truth uses new DGM and rejects misspecified EM score\n")

batch <- suppressWarnings(simulate_anchored(d,n_AC=1000,n_BC=1200,seeds=c(811,812),
  score_spec=scores,weight_covariates=paste0("z",1:6),run_exact=FALSE))
stopifnot(batch$complete,identical(batch$design,d),
          all(c(811,812) %in% batch$reference_truth$r),
          identical(batch$variables$score_spec$em_covariates,paste0("z",1:3)),
          identical(batch$variables$score_spec$non_em_covariates,paste0("z",4:6)))
summary <- summarize_anchored_simulation(batch)
lb <- subset(batch$estimates,method=="LP" & component=="treatment")
stopifnot(nrow(lb)==4L,all(lb$status=="available"),all(is.finite(lb$estimate)))
metrics <- subset(summary$metrics,method=="LP" & component=="treatment")
stopifnot(all(metrics$n_evaluated==2L),all(is.finite(metrics$bias)))
first_truth <- subset(batch$reference_truth,r==811,select=-r)
variables$score_spec <- resolved
expected_truth <- brm_score_coefficient_truth(d,variables)
for (name in names(expected_truth)) close_to(first_truth[[name]],expected_truth[[name]])
cat("PASS: new DGM runs through two-replicate estimation and truth-aligned summaries\n")
print(summary$method_status)

## Default configuration must preserve all three EMs and all three non-EMs.
config <- new.env(parent=globalenv())
sys.source("estimation/simulation_config.R",envir=config)
stopifnot(identical(config$design$dgm,"webster2024_rr_3em3nonem"),
          identical(config$scores$em_covariates,paste0("z",1:3)),
          identical(config$scores$non_em_covariates,paste0("z",4:6)))
cat("PASS: default configuration uses requested three-plus-three log-risk design\n")
