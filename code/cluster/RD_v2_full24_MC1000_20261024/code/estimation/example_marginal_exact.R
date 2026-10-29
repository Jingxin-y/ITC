## From anchored_v1:
## source('estimation/example_marginal_exact.R', encoding='UTF-8')
## ans <- example_marginal_exact(method='brm')
## ans$ci[c('estimate','low','up','p','status')]
## ans$ci$diagnostics$counts
## No simulation runs until example_marginal_exact() is explicitly called.

example_marginal_exact <- function(code_dir='.', method=c('brm','brm-FC'),
    n=200L, seed=1L, weight_covariates=paste0('z',1:3), ci_control=list()) {
  method <- match.arg(method)
  source(file.path(code_dir,'load_anchored.R'),encoding='UTF-8')
  source_anchored(code_dir)
  design <- joint_rr_gop_design(event='common',ess_ratio=.5)
  result <- run.anchored(design,n_AC=n,n_BC=n,seed=seed,run_exact=FALSE,
    model_weighting='unit',
    weight_covariates=weight_covariates,model_covariates=c('h_E','h_P'),
    effect_modifiers='h_E',non_effect_modifiers='h_P',op_model='marginal',
    score_spec=make_brm_score_spec(em_covariates=paste0('z',1:3),
                                  non_em_covariates=paste0('z',4:6)),
    return_full=TRUE,return_data=TRUE)
  fit <- result$fits[[method]]
  if(result$weights$status!='available' || is.null(fit) || !isTRUE(fit$converged))
    stop('Weight calibration or the requested outcome fit did not converge.')
  dat <- result$data$data.IPD$data
  ## The likelihood uses the outcome model's unit weights. MAIC weights enter
  ## both the marginal risk average and the intercept constraint, including in
  ## every bootstrap refit. The shared adapter also protects these weights from
  ## being overwritten by ci_control.
  ci <- marginal_exact_from_fit(fit,dat$y,dat$x,
    weight=result$model_weights$weight,
    standardization.weight=result$weights$weight,
    run_exact=TRUE,seed=seed+1000000L,bc_control=ci_control,
    label=paste(method,'marginal BC'))
  list(ci=ci,fit=fit,weights=result$weights,model_weights=result$model_weights,
       population_truth=design$marginal_truth['BC','log_RR_AC'])
}
