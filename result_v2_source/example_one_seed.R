# Interface example only: one seed, no bootstrap BC, no files written.
args <- commandArgs(trailingOnly=TRUE)
scale <- match.arg(if(length(args)) args[1L] else "RR",c("RR","RD"))
if(.Platform$OS.type=="windows") Sys.setlocale("LC_CTYPE","English_United States.utf8")
code <- normalizePath(file.path(scale,"code"),winslash="/",mustWork=TRUE)
if(scale=="RR") {
  source(file.path(code,"load_anchored.R"),encoding="UTF-8")
  source_anchored(code)
  design <- webster_rr_scenario_design("rare",.5,"diffuse","high")
  result <- run.anchored(design=design,n_AC=100,n_BC=100,seed=1,
    score_spec=make_brm_score_spec(),weight_covariates=paste0("z",1:6),
    model_weighting="unit",run_exact=FALSE)
  print(result$marginal_estimates)
} else {
  source(file.path(code,"RD/load_anchored_rd.R"),encoding="UTF-8")
  source_anchored_rd(code,file.path(code,"RD"))
  design <- webster_rd_scenario_design("rare",.5,"diffuse","high")
  result <- run.anchored.rd(design=design,n_AC=100,n_BC=100,seed=1,
    weight_covariates=paste0("z",1:6),model_weighting="unit",estimation_mode="scores",
    brm_effect_model="linear_atanh",run_exact=FALSE)
  print(result$estimates)
}
cat("One-seed interface example completed; BC disabled.\n")
