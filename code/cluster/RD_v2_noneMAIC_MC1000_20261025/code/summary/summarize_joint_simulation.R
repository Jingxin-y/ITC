## Point/SE/CI/P summaries use explicit statuses. Truth can vary by replicate r.
summarize_anchored_simulation <- function(result, truth = NULL) {
  if (!inherits(result,"anchored_simulation")) stop("result must be returned by simulate_anchored().")
  dat <- result$estimates
  if (is.null(dat) || !nrow(dat)) return(list(metrics=data.frame(),method_status=data.frame(),parameter_status=data.frame()))
  keys <- c("method","component","term")
  key <- function(x) do.call(paste,c(x[keys],list(sep="\034")))
  repkey <- function(x) paste(x$r,key(x),sep="\034")
  schema <- unique(rbind(dat[keys],result$parameter_schema[keys]))
  schema <- schema[complete.cases(schema),,drop=FALSE]
  if (is.null(truth)) truth <- result$reference_truth
  by_rep <- !is.null(truth) && "r" %in% names(truth)
  if (!is.null(truth)) {
    if (!is.data.frame(truth) || !all(c(keys,"truth") %in% names(truth)) ||
        !is.numeric(truth$truth) || any(!is.finite(truth$truth)) || anyNA(truth[keys]))
      stop("truth must have method/component/term keys and finite numeric truth values.")
    if (by_rep && (!is.numeric(truth$r) || anyNA(truth$r) || any(!truth$r %in% result$seeds)))
      stop("Per-replicate truth must use completed simulation seeds in column r.")
    if (anyDuplicated(if (by_rep) repkey(truth) else key(truth))) stop("truth contains duplicate keys.")
    if (any(!key(truth) %in% key(schema))) stop("truth contains an unknown method/component/term.")
    if (!by_rep && !is.null(result$reference_truth) &&
        any(key(truth) %in% key(result$reference_truth)))
      stop("BC-estimated score units vary by replicate. Use truth=NULL for aligned reference_truth or supply truth with column r; do not reuse fixed DGM coefficients.")
  }
  ## A finite BC endpoint in an older file without status diagnostics is unverified.
  for (field in c("se_status","interval_status","p_status")) {
    if (!field %in% names(dat)) {
      dat[[field]] <- dat$status
      if (field != "se_status") dat[[field]][dat$interval %in% c("BC","BC_not_run")] <- "BC_unverified"
    }
  }
  effects <- dat[dat$component %in% c("alpha","treatment","common_effect","marginal"),,drop=FALSE]
  targets <- schema[schema$component %in% c("alpha","treatment","common_effect","marginal"),keys,drop=FALSE]
  finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  finite_sd <- function(x) if (sum(is.finite(x))>1L) sd(x[is.finite(x)]) else NA_real_
  metrics <- lapply(seq_len(nrow(targets)),function(i) {
    target <- targets[i,,drop=FALSE]
    d <- effects[key(effects)==key(target),,drop=FALSE]
    true <- if (is.null(truth)) rep(NA_real_,nrow(d)) else
      truth$truth[match(if (by_rep) repkey(d) else key(d),if (by_rep) repkey(truth) else key(truth))]
    valid_est <- d$status %in% "available" & is.finite(d$estimate)
    valid_se <- valid_est & d$se_status %in% "available" & is.finite(d$se) & d$se>0
    valid_ci <- valid_est & d$interval_status %in% c("available","BC_boundary_truncated") &
      is.finite(d$lower) & is.finite(d$upper) & d$lower<=d$upper
    valid_p <- valid_est & d$p_status %in% "available" & is.finite(d$p.value) & d$p.value>=0 & d$p.value<=1
    evaluate <- valid_est & is.finite(true)
    coverage_rows <- valid_ci & is.finite(true)
    errors <- d$estimate[evaluate]-true[evaluate]
    sd_est <- finite_sd(d$estimate[valid_est])
    sd_error <- finite_sd(errors)
    has_truth <- any(is.finite(true))
    mc_sd <- if (has_truth) sd_error else sd_est
    avg_se <- finite_mean(d$se[valid_se])
    aligned_se <- finite_mean(d$se[valid_se & if (has_truth) is.finite(true) else TRUE])
    n_not_applicable <- sum(grepl("^not_applicable",d$status))
    n_applicable <- length(result$seeds)-n_not_applicable
    truth_values <- true[is.finite(true)]
    varies <- length(truth_values)>1L && diff(range(truth_values))>1e-12
    cbind(target,data.frame(
      n_total=length(result$seeds),n_estimate=sum(valid_est),n_se=sum(valid_se),
      n_interval=sum(valid_ci),n_boundary_truncated=sum(valid_ci & d$interval_status=="BC_boundary_truncated"),
      n_p=sum(valid_p),n_evaluated=sum(evaluate),n_coverage=sum(coverage_rows),
      n_not_applicable=n_not_applicable,n_estimate_failed=n_applicable-sum(valid_est),
      estimate_failure_rate=if (n_applicable>0) 1-sum(valid_est)/n_applicable else NA_real_,
      truth=if (length(truth_values) && !varies) truth_values[1] else NA_real_,
      mean_truth=finite_mean(truth_values),truth_varies=varies,
      mean_estimate=finite_mean(d$estimate[valid_est]),mean_se=avg_se,
      empirical_sd=mc_sd,empirical_sd_estimate=sd_est,empirical_sd_error=sd_error,
      se_ratio=if (is.finite(mc_sd) && mc_sd>0) aligned_se/mc_sd else NA_real_,
      bias=finite_mean(errors),
      coverage=if (any(coverage_rows)) mean(d$lower[coverage_rows]<=true[coverage_rows] &
        true[coverage_rows]<=d$upper[coverage_rows]) else NA_real_,
      rejection_rate_0=if (any(valid_p)) mean(d$p.value[valid_p]<=.05) else NA_real_))
  })
  status <- do.call(rbind,lapply(split(dat,interaction(dat$r,dat$method,drop=TRUE)),function(d) {
    effect <- d$component %in% c("alpha","treatment","common_effect","marginal")
    s <- unique(d$status[if (any(effect)) effect else rep(TRUE,nrow(d))])
    label <- if (length(s)==1L) s else if ("available" %in% s) "partially_available" else "fit_failed"
    data.frame(r=d$r[1],method=d$method[1],status=label)
  }))
  method_status <- aggregate(list(n=rep(1L,nrow(status))),status[c("method","status")],sum)
  parameter_status <- lapply(c("status","se_status","interval_status","p_status"),function(field) {
    d <- unique(dat[c("r",keys,field)])
    d <- d[complete.cases(d),,drop=FALSE]
    if (!nrow(d)) return(NULL)
    names(d)[names(d)==field] <- "status"
    cbind(output=field,aggregate(list(n=rep(1L,nrow(d))),d[c(keys,"status")],sum))
  })
  list(metrics=if (length(metrics)) do.call(rbind,metrics) else data.frame(),
       method_status=method_status,parameter_status=do.call(rbind,parameter_status))
}
