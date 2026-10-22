## 仅定义估计阶段的分数。数据生成器 rr_gop_scores() 及 DGM 真值保持独立。
make_brm_score_spec <- function(
    em_covariates = c("z1", "z2", "z3"),
    non_em_covariates = setdiff(paste0("z", 1:6), em_covariates),
    reverse_covariates = intersect(c("z1", "z3", "z5"), c(em_covariates, non_em_covariates)),
    coefficients = NULL,
    reference_prob = c(z1=.7, z2=.3, z3=.7, z4=.3, z5=.7, z6=.3),
    data.ad = NULL) {
  check_names <- function(x, label) {
    if (!is.character(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
      stop(label, " must be a vector of distinct column names.", call.=FALSE)
    if (any(x %in% c("x", "y", "S_E", "S_P", "h_E", "h_P")))
      stop(label, " must name raw baseline covariates, not treatment, outcome or score columns.", call.=FALSE)
  }
  check_names(em_covariates, "em_covariates")
  check_names(non_em_covariates, "non_em_covariates")
  check_names(reverse_covariates, "reverse_covariates")
  if (length(intersect(em_covariates, non_em_covariates)))
    stop("EM and non-EM score components must be disjoint.", call.=FALSE)
  cols <- c(em_covariates, non_em_covariates)
  if (any(!reverse_covariates %in% cols))
    stop("reverse_covariates must be selected score components.", call.=FALSE)
  if (length(cols)) {
    if (!is.numeric(reference_prob) || is.null(names(reference_prob)) ||
        anyDuplicated(names(reference_prob)) || any(!cols %in% names(reference_prob)))
      stop("reference_prob must supply named probabilities for all selected components.", call.=FALSE)
    reference_prob <- reference_prob[cols]
    if (any(!is.finite(reference_prob)) || any(reference_prob < 0 | reference_prob > 1))
      stop("reference_prob values must be finite and in [0,1].", call.=FALSE)
  } else reference_prob <- setNames(numeric(0), character(0))
  if (is.null(coefficients)) coefficients <- setNames(rep(1, length(cols)), cols)
  if (!is.numeric(coefficients) || (length(cols) &&
      (is.null(names(coefficients)) || anyDuplicated(names(coefficients)) ||
       any(!cols %in% names(coefficients)))))
    stop("coefficients must supply named numeric coefficients for all selected components.", call.=FALSE)
  coefficients <- coefficients[cols]
  if (any(!is.finite(coefficients)) || any(coefficients == 0))
    stop("Score coefficients must be finite and nonzero; omit unused components.", call.=FALSE)
  groups <- list(h_E=em_covariates, h_P=non_em_covariates)
  center <- scale <- setNames(numeric(2), names(groups))
  for (name in names(groups)) {
    selected <- groups[[name]]
    p <- reference_prob[selected]
    a <- coefficients[selected]
    p[selected %in% reverse_covariates] <- 1 - p[selected %in% reverse_covariates]
    center[name] <- sum(a * p)
    scale[name] <- if (length(selected)) sqrt(sum(a^2 * p * (1-p))) else 1
    if (scale[name] <= 0) stop("Selected score has zero reference variance: ", name, call.=FALSE)
  }
  spec <- structure(list(em_covariates=em_covariates, non_em_covariates=non_em_covariates,
    reverse_covariates=reverse_covariates, coefficients=coefficients,
    reference_prob=reference_prob, center=center, scale=scale), class="brm_score_spec")
  if (!is.null(data.ad)) spec <- brm_scores_from_bc(spec,data.ad)
  spec
}

brm_score_transform <- function(score_spec, score) {
  if (!inherits(score_spec, "brm_score_spec"))
    stop("score_spec must be created by make_brm_score_spec().", call.=FALSE)
  if (!score %in% c("S_E", "S_P", "h_E", "h_P")) stop("Unknown score column.")
  key <- if (score %in% c("S_E", "h_E")) "h_E" else "h_P"
  cols <- if (key == "h_E") score_spec$em_covariates else score_spec$non_em_covariates
  a <- score_spec$coefficients[cols]
  reverse <- cols %in% score_spec$reverse_covariates
  offset <- sum(a[reverse])
  a[reverse] <- -a[reverse]
  if (startsWith(score, "h_")) {
    offset <- (offset - score_spec$center[key]) / score_spec$scale[key]
    a <- a / score_spec$scale[key]
  }
  list(coefficients=a, offset=unname(offset))
}

## Resolve analysis units using reported BC summaries only. DGM scores stay untouched.
brm_scores_from_bc <- function(score_spec, data.ad) {
  for (group in c("E", "P")) {
    key <- paste0("h_", group)
    tr <- brm_score_transform(score_spec, paste0("S_", group))
    a <- tr$coefficients
    if (!length(a)) {
      score_spec$center[key] <- 0; score_spec$scale[key] <- 1
      next
    }
    mu <- data.ad$anchored.mean.x
    V <- data.ad$anchored.cov.x
    cols <- names(a)
    if (is.null(names(mu)) || any(!cols %in% names(mu)) || any(!is.finite(mu[cols])))
      stop("BC score standardization requires reported means for all selected raw covariates.")
    if (!is.matrix(V) || any(!cols %in% rownames(V)) || any(!cols %in% colnames(V)))
      stop("BC score standardization requires the reported covariance matrix; marginal variances alone are insufficient for a sum score.")
    V <- V[cols,cols,drop=FALSE]
    if (any(!is.finite(V)) || max(abs(V-t(V))) > 1e-8)
      stop("BC score covariance must be finite and symmetric.")
    variance <- as.numeric(crossprod(a,V %*% a))
    if (!is.finite(variance) || variance <= 0)
      stop("Selected score has nonpositive reported BC variance: ", key)
    score_spec$center[key] <- tr$offset + sum(a * mu[cols])
    score_spec$scale[key] <- sqrt(variance)
  }
  score_spec$standardization <- "bc_estimated"
  score_spec$reference_n <- data.ad$anchored.n
  score_spec
}

## Evaluation only, after fitting: true coefficients in this replicate's analysis units.
## A misspecified EM score has no exact linear RR coefficient truth; do not invent one.
brm_score_coefficient_truth <- function(design, variables) {
  spec <- variables$score_spec
  if (is.null(spec) || !identical(variables$score_standardization,"bc_estimated") ||
      !identical(spec$standardization,"bc_estimated") ||
      !identical(variables$effect_modifiers,"h_E")) return(NULL)
  z <- expand.grid(setNames(rep(list(0:1),6),paste0("z",1:6)))
  ## Read the selected DGM's risk table only for post-fit truth evaluation.
  ## Direct log-risk and RR-GOP designs have different parameterizations.
  scores <- rr_gop_scores(z)
  index <- 1L + scores$S_E + 4L * scores$S_P
  true_rr <- design$risk_table$log_RR_AC[index]
  h <- apply_brm_scores(z,spec)$h_E
  X <- cbind(1,h)
  fit <- lm.fit(X,true_rr)
  if (fit$rank != 2L || any(!is.finite(fit$coefficients)) || max(abs(fit$residuals)) > 1e-10)
    return(NULL)
  brm <- expand.grid(method=c("brm","brm-BC","brm-FC","brm-FC-BC"),
                    component="alpha",term=c("(Intercept)","h_E"),stringsAsFactors=FALSE)
  glm <- expand.grid(method=c("LB","LP","RLP","RLB",
                              anchored_glm_fc_method_names()),
                    component="treatment",term=c("x","x:h_E"),stringsAsFactors=FALSE)
  truth <- rbind(brm,glm)
  if (identical(variables$bc_target,"marginal"))
    truth <- truth[!truth$method %in% c("brm-BC","brm-FC-BC"),,drop=FALSE]
  truth$truth <- unname(fit$coefficients[ifelse(truth$term %in% c("h_E","x:h_E"),2L,1L)])
  truth
}

## 接受 IPD/AD 个体数据对象或 data.frame；返回估计用副本，不修改传入对象或其 truth。
apply_brm_scores <- function(data, score_spec) {
  brm_score_transform(score_spec, "h_E") # validate constructor
  is_frame <- is.data.frame(data)
  dat <- if (is_frame) data else data$data
  if (!is.data.frame(dat)) stop("data must be a data.frame or a list containing data.")
  cols <- c(score_spec$em_covariates, score_spec$non_em_covariates)
  if (any(!cols %in% names(dat)))
    stop("Unknown score components: ", paste(setdiff(cols,names(dat)),collapse=", "), call.=FALSE)
  if (any(!vapply(dat[cols], function(z)
      (is.numeric(z) || is.logical(z)) && !anyNA(z) && all(z %in% c(0,1)), logical(1))))
    stop("Custom score components must be observed binary 0/1 covariates.", call.=FALSE)
  for (group in c("E", "P")) {
    selected <- if (group == "E") score_spec$em_covariates else score_spec$non_em_covariates
    s <- rep(0, nrow(dat))
    for (name in selected) {
      z <- dat[[name]]
      if (name %in% score_spec$reverse_covariates) z <- 1-z
      s <- s + score_spec$coefficients[name] * z
    }
    ## 先求离散总分再标准化，避免等价取值因浮点消去误差被分成不同 OP 格。
    key <- paste0("h_",group)
    dat[[paste0("S_",group)]] <- s
    dat[[key]] <- (s - score_spec$center[key]) / score_spec$scale[key]
  }
  if (is_frame) {
    attr(dat,"brm_score_spec") <- score_spec
    return(dat)
  }
  data$data <- dat
  data$score_spec <- score_spec
  data
}

prepare_brm_score_model <- function(data, score_spec, effect_modifiers,
                                    non_effect_modifiers, model_covariates) {
  if (is.null(score_spec)) score_spec <- data$score_spec
  if (is.null(score_spec)) score_spec <- attr(data$data,"brm_score_spec")
  if (!is.null(score_spec)) {
    data <- apply_brm_scores(data, score_spec)
    em <- if (length(score_spec$em_covariates)) "h_E" else character(0)
    non_em <- if (length(score_spec$non_em_covariates)) "h_P" else character(0)
    if (is.null(effect_modifiers)) effect_modifiers <- em
    if (is.null(non_effect_modifiers)) non_effect_modifiers <- non_em
    if (is.null(model_covariates)) model_covariates <- c(em, non_em)
  }
  list(data=data, score_spec=score_spec, effect_modifiers=effect_modifiers,
       non_effect_modifiers=non_effect_modifiers, model_covariates=model_covariates)
}

## 自定义分数用于 weighting 时，用原始 AD 均值/协方差作仿射变换，避免沿用旧分数汇总。
brm_score_summary <- function(data.ad, score_spec, score_columns=c("h_E","h_P")) {
  if (!length(score_columns)) return(data.ad)
  transforms <- lapply(score_columns, function(name) brm_score_transform(score_spec,name))
  mu <- data.ad$anchored.mean.x
  if (is.null(names(mu)) || anyDuplicated(names(mu)))
    stop("Custom score weighting requires named raw AD means.", call.=FALSE)
  raw <- setdiff(names(mu),c("S_E","S_P","h_E","h_P"))
  needed <- unique(unlist(lapply(transforms,function(t)names(t$coefficients))))
  if (any(!needed %in% raw))
    stop("Custom score weighting requires raw AD summaries for: ",
         paste(setdiff(needed,raw),collapse=", "),call.=FALSE)
  V <- data.ad$anchored.cov.x
  if (!is.matrix(V) || is.null(rownames(V)) || is.null(colnames(V)) ||
      any(!raw %in% rownames(V)) || any(!raw %in% colnames(V)))
    stop("Custom score weighting requires the named raw AD covariance matrix.", call.=FALSE)
  B <- matrix(0,length(raw)+length(score_columns),length(raw),
              dimnames=list(c(raw,score_columns),raw))
  if (length(raw)) B[seq_along(raw),] <- diag(length(raw))
  offset <- rep(0,nrow(B))
  for (i in seq_along(transforms)) {
    a <- transforms[[i]]$coefficients
    B[length(raw)+i,names(a)] <- a
    offset[length(raw)+i] <- transforms[[i]]$offset
  }
  data.ad$anchored.mean.x <- setNames(as.vector(B %*% mu[raw]) + offset,rownames(B))
  data.ad$anchored.cov.x <- B %*% V[raw,raw,drop=FALSE] %*% t(B)
  data.ad$anchored.var.x <- diag(data.ad$anchored.cov.x)
  data.ad
}
