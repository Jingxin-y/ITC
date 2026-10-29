## RD-only compatibility bridge to the unchanged shared solvers.
## Functions are rebound in a private environment; global RR/legacy RD entry points stay intact.
rd_getProbRD <- function(atanhrd, logop = NA_real_) {
  if (is.matrix(atanhrd) && ncol(atanhrd) == 2L) {
    logop <- atanhrd[, 2L]; atanhrd <- atanhrd[, 1L]
  } else if (length(logop) == 1L && is.na(logop) && length(atanhrd) == 2L) {
    logop <- atanhrd[2L]; atanhrd <- atanhrd[1L]
  }
  p <- rd_probabilities(atanhrd, logop)
  cbind(p0 = p$p0, p1 = p$p1)
}
rd_getProbScalarRD <- function(atanhrd, logop = NA_real_) {
  if (length(atanhrd) == 2L) {
    logop <- atanhrd[2L]; atanhrd <- atanhrd[1L]
  }
  if (length(atanhrd) != 1L || length(logop) != 1L)
    stop("rd_getProbScalarRD requires one RD/OP pair.")
  as.vector(rd_getProbRD(atanhrd, logop))
}

initialize_rd_solver_adapter <- function() {
  shared <- environment(get.estimate)
  private <- new.env(parent = shared)
  names <- c("get.estimate", "MLEst", "max.likelihood", "var.mle.rd", "var.mle.rd.joint",
    "max.likelihood.firth.rd", "max.likelihood.firth", "firth_rd_derivatives")
  for (key in names) {
    f <- get(key, envir = shared, inherits = TRUE)
    if (!is.function(f)) stop("Shared RD solver is not a function: ", key)
    environment(f) <- private
    assign(key, f, envir = private)
  }
  for (key in c("rd_probabilities", "rd_getProbRD", "rd_getProbScalarRD")) {
    f <- get(key, envir = environment(initialize_rd_solver_adapter), inherits = TRUE)
    environment(f) <- private
    assign(key, f, envir = private)
  }
  private$getProbRD <- private$rd_getProbRD
  private$getProbScalarRD <- private$rd_getProbScalarRD
  private$probability_mapping <- "stable_rd_op_v1"
  private
}
.anchored_rd_solver_env <- initialize_rd_solver_adapter()

rd_get_estimate <- function(..., param = "RD") {
  if (!identical(param, "RD")) stop("rd_get_estimate supports RD only.")
  fit <- do.call(.anchored_rd_solver_env$get.estimate, c(list(...), list(param = param)))
  fit$probability_mapping <- .anchored_rd_solver_env$probability_mapping
  fit
}
rd_solver_derivatives <- function(...) .anchored_rd_solver_env$firth_rd_derivatives(...)
