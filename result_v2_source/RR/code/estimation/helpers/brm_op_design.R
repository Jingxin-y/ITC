## BRM 的 reduced outcome model：在当前加权人群中 marginalize 非 EM，
## 保留 effect_modifiers 的联合取值，每个观察到的 EM 格自由估计一个 log OP。
## 直接拟合 Y | treatment, EM；不平均条件 OP，不使用 DGM 的真实风险。
brm_op_design <- function(data, effect_modifiers = NULL,
                          non_effect_modifiers = NULL, model_covariates = NULL,
                          op_model = c("marginal", "linear")) {
  op_model <- match.arg(op_model)
  validate_names <- function(cols, label) {
    if (is.null(cols)) return(character(0))
    if (!is.character(cols) || anyNA(cols) || any(!nzchar(cols)) || anyDuplicated(cols))
      stop(label, " must contain distinct column names.", call. = FALSE)
    if (any(!cols %in% names(data)))
      stop(label, " contains unknown columns: ", paste(setdiff(cols, names(data)), collapse = ", "),
           call. = FALSE)
    if (any(cols %in% c("x", "y")))
      stop(label, " must contain baseline covariates, not x or y.", call. = FALSE)
    cols
  }
  em <- validate_names(effect_modifiers, "effect_modifiers")
  model_covariates <- validate_names(
    if (is.null(model_covariates)) covariate_names(data) else model_covariates,
    "model_covariates")
  non_em <- validate_names(
    if (is.null(non_effect_modifiers)) setdiff(model_covariates, em) else non_effect_modifiers,
    "non_effect_modifiers")
  if (length(intersect(em, non_em)))
    stop("effect_modifiers and non_effect_modifiers must be disjoint.", call. = FALSE)
  numeric_columns <- function(cols) {
    if (any(!vapply(data[cols], function(z)
      (is.numeric(z) || is.logical(z)) && all(is.finite(z)), logical(1))))
      stop("BRM design columns must be finite numeric columns; encode categorical EM first.",
           call. = FALSE)
  }
  numeric_columns(em)
  n <- nrow(data)
  if (!n) stop("BRM needs at least one observation.", call. = FALSE)
  va <- cbind("(Intercept)" = rep(1, n), as.matrix(data[, em, drop = FALSE]))
  storage.mode(va) <- "double"
  cells <- NULL
  cell_id <- NULL
  if (op_model == "linear") {
    numeric_columns(model_covariates)
    vb <- cbind("(Intercept)" = rep(1, n), covariate_matrix(data, model_covariates))
  } else {
    ## 整数编码防止浮点值格式化后碰撞；只保留实际出现的联合格，避免空 dummy。
    if (length(em)) {
      codes <- lapply(data[em], function(z) match(z, sort(unique(z))))
      key <- do.call(paste, c(codes, list(sep = ":")))
      rows <- which(!duplicated(key))
      rows <- rows[do.call(order, unname(data[rows, em, drop = FALSE]))]
      cell_id <- match(key, key[rows])
      cells <- data[rows, em, drop = FALSE]
      rownames(cells) <- NULL
    } else {
      cell_id <- rep.int(1L, n)
      cells <- data.frame(row.names = 1L)
    }
    k <- nrow(cells)
    if (k + ncol(va) > n)
      stop("Saturated OP has more parameters than observations; use discrete EM or a prespecified grouped score.",
           call. = FALSE)
    vb <- matrix(0, n, k, dimnames = list(NULL,
      c("(Intercept)", if (k > 1L) paste0("OP_cell_", 2:k) else character(0))))
    vb[, 1L] <- 1
    non_reference <- which(cell_id > 1L)
    if (length(non_reference)) vb[cbind(non_reference, cell_id[non_reference])] <- 1
    cells <- cbind(data.frame(cell = seq_len(k)), cells,
                   n = tabulate(cell_id, nbins = k))
  }
  storage.mode(vb) <- "double"
  list(va = va, vb = vb, cell_id = cell_id,
    spec = list(mode = op_model, effect_modifiers = em, non_effect_modifiers = non_em,
      model_covariates = model_covariates,
      conditioning_covariates = if (op_model == "marginal") em else union(em, model_covariates),
      cells = cells, n_nuisance = ncol(vb)))
}
