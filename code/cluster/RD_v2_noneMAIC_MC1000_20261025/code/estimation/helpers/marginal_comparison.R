## Main-pipeline adapter for marginal log RR bootstrap intervals.
## The fitting weights and target standardization weights have separate roles.
validate_marginal_bc_control <- function(control) {
  allowed <- c("max.step", "optim.maxit", "optim.reltol", "score.tol", "level",
    "psi.se", "psi.range", "search.mult", "thres.dicho",
    "bisection.max.step", "grid.mult", "nsim.coarse", "nsim.confirm",
    "nsim.refine", "nsim.p.value", "alpha.bound", "beta.bound", "root.tol",
    "root.maxit", "constraint.tol", "time.limit.sec", "compute.p.value",
    "compute.interval", "psi.grid", "optim.retries", "optim.retry.maxit",
    "optim.multistart", "optim.multistart.max", "optim.fallback")
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("bc_control must be a named list without duplicate names.")
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) stop("Unsupported marginal bc_control: ", paste(unknown, collapse=", "),
    ". Use exact_marginal_target tuning parameters; coefficient-only controls require bc_target='coefficient'.")
  invisible(control)
}

make_na_marginal_exact <- function(status = "not_requested", estimate = NA_real_) {
  list(estimate=estimate, low=NA_real_, up=NA_real_, p=NA_real_,
    effect.estimate=exp(estimate), effect.low=NA_real_, effect.up=NA_real_,
    target="marginal log RR", status=status, interval.available=FALSE,
    p.status=status, diagnostics=list())
}

marginal_exact_from_fit <- function(fit, y, x, weight, standardization.weight,
                                    run_exact, seed=NULL, bc_control=list(), label="marginal BC") {
  validate_marginal_bc_control(bc_control)
  if (is.null(fit$va) || is.null(fit$vb) || is.null(fit$point.est))
    return(make_na_marginal_exact("fit_failed"))
  if (!isTRUE(fit$converged)) return(make_na_marginal_exact("fit_not_converged"))
  target <- tryCatch(marginal_logrr_target(fit$point.est, fit$va, fit$vb,
                                         standardization.weight, TRUE), error=identity)
  if (inherits(target,"error")) {
    out <- make_na_marginal_exact("target_failed")
    out$message <- conditionMessage(target)
    return(out)
  }
  out <- make_na_marginal_exact(estimate=target$value)
  if (!run_exact) return(out)
  args <- list(param="RR", y=y, x=x, va=fit$va, vb=fit$vb, pars=fit$point.est,
    weight=weight, standardization.weight=standardization.weight, seed=seed)
  ## Covariance sets only the search width. Do not report it as a bootstrap SE.
  if (is.null(bc_control$psi.se) && is.null(bc_control$psi.range)) {
    V <- fit$cov
    valid <- is.matrix(V) && identical(dim(V), c(length(fit$point.est),length(fit$point.est))) &&
      all(is.finite(V)) && max(abs(V-t(V))) <= 1e-8*max(1,max(abs(V)))
    variance <- if (valid) drop(crossprod(target$gradient,V%*%target$gradient)) else NA_real_
    if (is.finite(variance) && variance>0) args$vcov.pars <- V else {
      ## Without EM the target is exactly alpha0, including the scalar Bayes correction.
      se <- if (ncol(fit$va)==1L) fit$se.est[1] else NA_real_
      if (length(se)!=1L || !is.finite(se) || se<=0) {
        warning(label, " has no usable marginal SE; using search width SE=1.",call.=FALSE)
        se <- 1
      }
      args$psi.se <- se
    }
  }
  args[names(bc_control)] <- bc_control
  ci <- tryCatch(do.call(exact_marginal_target,args), error=identity)
  if (inherits(ci,"error")) {
    out$status <- out$p.status <- if (inherits(ci,"marginal_exact_timeout")) "timeout" else "error"
    out$message <- conditionMessage(ci)
    out$diagnostics <- list(error=out$message, error.class=class(ci))
    warning(label, ": ", out$message, call.=FALSE)
    return(out)
  }
  ci
}

comparison_marginal_rows <- function(method, fit, bc, weights, run_exact,
                                     status="available", retain_nonconverged=FALSE) {
  out <- data.frame(method=method, component="marginal", parameter="log_RR",
    term="marginal log RR", estimate=NA_real_, se=NA_real_, lower=NA_real_, upper=NA_real_,
    p.value=NA_real_, interval=if(run_exact) "BC" else "BC_not_run", status=status,
    se_status="not_provided", interval_status=if(run_exact) "BC_failed" else "BC_not_run",
    p_status=if(run_exact) "BC_failed" else "BC_not_run", stringsAsFactors=FALSE)
  if (status=="available") {
    if (is.null(fit$va) || is.null(fit$vb) || is.null(fit$point.est)) out$status <- "fit_failed" else {
      if (!isTRUE(fit$converged)) out$status <- "fit_not_converged"
      if (isTRUE(fit$converged) || retain_nonconverged) {
        out$estimate <- tryCatch(marginal_logrr_target(fit$point.est,fit$va,fit$vb,weights),
                                 error=function(e) NA_real_)
        if (!is.finite(out$estimate) && out$status=="available") out$status <- "estimate_nonfinite"
      }
    }
  }
  if (out$status!="available") {
    out$se_status <- out$interval_status <- out$p_status <- out$status
    return(out)
  }
  if (!run_exact || is.null(bc)) return(out)
  bounds <- c(bc$low,bc$up)
  if (isTRUE(bc$interval.available) && bc$status %in% c("available","boundary_truncated") &&
      length(bounds)==2L && all(is.finite(bounds)) && bounds[1]<=bounds[2]) {
    out$lower <- bounds[1]; out$upper <- bounds[2]
    out$interval_status <- if(bc$status=="boundary_truncated") "BC_boundary_truncated" else "available"
  } else if (!is.null(bc$status)) out$interval_status <- paste0("BC_",bc$status)
  if (identical(bc$p.status,"available") && length(bc$p)==1L && is.finite(bc$p) && bc$p>=0 && bc$p<=1) {
    out$p.value <- bc$p; out$p_status <- "available"
  } else if (!is.null(bc$p.status)) out$p_status <- paste0("BC_",bc$p.status)
  out
}
