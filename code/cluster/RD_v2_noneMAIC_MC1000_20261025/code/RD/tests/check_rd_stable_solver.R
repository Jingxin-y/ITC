## Verify stable probabilities AND their use by actual shared-solver adapters.
args <- commandArgs(trailingOnly = TRUE)
own <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
rd_dir <- dirname(dirname(normalizePath(own, winslash = "/")))
code_dir <- if (length(args)) args[1] else dirname(rd_dir)
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source(file.path(code_dir, "load_anchored.R"), encoding = "UTF-8"); source_anchored(code_dir)
keys <- c("getProbRR", "getProbScalarRR", "getProbRD", "getProbScalarRD", "get.estimate",
  "max.likelihood", "var.mle.rd", "firth_rd_derivatives", "max.likelihood.firth.rd")
fingerprint <- function() lapply(setNames(keys, keys), function(key) {
  f <- get(key, envir = .GlobalEnv)
  list(body = deparse(body(f)), formals = deparse(formals(f)), environment = environment(f))
})
original <- fingerprint()
rr_ref <- getProbRR(c(-3, 0, .2, 5), c(-15, -1, 0, 20))
source(file.path(rd_dir, "load_anchored_rd.R"), encoding = "UTF-8")
source_anchored_rd(code_dir, rd_dir)
check <- function(x, message) { if (!isTRUE(x)) stop(message); cat("PASS:", message, "\n") }
fails <- function(x) inherits(tryCatch(x, error = identity), "error")
check(identical(original, fingerprint()) && identical(rr_ref, getProbRR(c(-3, 0, .2, 5), c(-15, -1, 0, 20))),
  "loading RD stable adapter leaves all global RR/shared functions and RR outputs unchanged")

grid <- expand.grid(theta = c(-20, -3, -.8, -1e-14, 0, 1e-14, .8, 3, 20),
  phi = c(-1e308, -1500, -1000, -100, -40, -13, -12.001, -12, -11.999,
    -1, -1e-12, 0, 1e-12, 1, 11.999, 12, 12.001, 13, 40, 100, 1000, 1e308))
p <- rd_getProbRD(grid$theta, grid$phi)
check(all(is.finite(p)) && all(p >= 0 & p <= 1) &&
  max(abs(p[, 2] - p[, 1] - tanh(grid$theta))) < 1e-12,
  "RD probabilities are finite/legal and preserve RD across zero OP, old thresholds and extreme predictors")
interior <- apply(p, 1, function(q) min(q) > 1e-10 && max(q) < 1 - 1e-10)
logop <- rowSums(log(p[interior, , drop = FALSE]) - log1p(-p[interior, , drop = FALSE]))
check(max(abs(logop - grid$phi[interior])) < 1e-5, "interior mapped risks preserve log odds product")
for (theta in c(-.8, 0, .8)) for (phi in c(-8, -1, 0, 1, 8)) {
  d <- tanh(theta)
  oracle <- uniroot(function(q) log(q) + log(q + d) - log1p(-q) - log1p(-q - d) - phi,
    c(max(0, -d), min(1, 1 - d)), tol = 1e-12)$root
  stopifnot(max(abs(rd_getProbScalarRD(theta, phi) - c(oracle, oracle + d))) < 2e-11)
}
check(TRUE, "probability mapping agrees with an independent risk-space root oracle")
check(identical(rd_getProbScalarRD(c(.4, -13)), rd_getProbScalarRD(.4, -13)) &&
  identical(rd_getProbRD(matrix(c(.4, -13), 1)), rd_getProbRD(.4, -13)) &&
  identical(rd_getProbRD(c(.4, -13)), rd_getProbRD(.4, -13)) &&
  identical(rd_getProbRD(0, c(-1, 0, 1)), rd_getProbRD(rep(0, 3), c(-1, 0, 1))) &&
  identical(rd_getProbRD(c(-.4, 0, .4), 0), rd_getProbRD(c(-.4, 0, .4), rep(0, 3))) &&
  identical(dim(rd_getProbRD(numeric(), numeric())), c(0L, 2L)),
  "scalar, vector, packed-pair, matrix, broadcasting and empty interfaces work")
check(fails(rd_getProbRD(NA_real_, 0)) && fails(rd_getProbRD(0, Inf)) &&
  fails(rd_getProbRD(1i, 0)) && fails(rd_getProbRD(1:3, 1:2)) &&
  fails(rd_getProbScalarRD(1:3, 0)), "invalid probability inputs are rejected")
check(all(rd_getProbScalarRD(0, -1000) > 0) &&
  max(abs(log(rd_getProbScalarRD(0, -1000)) + 500)) < 1e-12 &&
  rd_getProbScalarRD(1e-200, -1200)[1] > 0 &&
  getProbScalarRD(.4, -13)[1] == 0 && rd_getProbScalarRD(.4, -13)[1] > 0 &&
  getProbScalarRD(.4, 13)[2] == 1 && rd_getProbScalarRD(.4, 13)[2] < 1,
  "stable adapter removes old +/-12 approximations and prevents zero-RD/tiny-OP 0/0")

v <- matrix(1, 1, 1)
for (theta in c(-.4, 0, .4)) for (phi in c(-13, -1e-12, 0, 1e-12, 13)) for (x in 0:1) {
  pars <- c(theta, phi); step <- 1e-4
  f <- function(a) rd_getProbScalarRD(a[1], a[2])[x + 1L]
  e1 <- c(step, 0); e2 <- c(0, step)
  first <- c((f(pars + e1) - f(pars - e1)) / (2 * step),
    (f(pars + e2) - f(pars - e2)) / (2 * step))
  ## Independent implicit differentiation of RD and log-OP constraints.
  ## Difference these first derivatives to avoid cancellation of nearly equal risks in a Hessian.
  implicit_gradient <- function(a) {
    pp <- rd_getProbScalarRD(a[1], a[2]); vv <- pp * (1 - pp)
    J <- rbind(c(-1, 1), 1 / vv)
    c(solve(J, c(1 - tanh(a[1])^2, 0))[x + 1L], solve(J, c(0, 1))[x + 1L])
  }
  he1 <- c(1e-6, 0); he2 <- c(0, 1e-6)
  H <- cbind((implicit_gradient(pars + he1) - implicit_gradient(pars - he1)) / 2e-6,
    (implicit_gradient(pars + he2) - implicit_gradient(pars - he2)) / 2e-6)
  second <- c(H[1, 1], mean(c(H[1, 2], H[2, 1])), H[2, 2])
  z <- rd_solver_derivatives(theta, phi, x, v, v)
  stopifnot(isTRUE(z$interior), max(abs(first - c(z$d.a, z$d.b))) < 1e-7,
    max(abs(second - c(z$d.aa, z$d.ab, z$d.bb)) / (1 + abs(second))) < 1e-5)
}
check(TRUE, "FC first/second derivatives agree with independent differences of the stable map")

local({
  poisoned <- c("getProbRD", "getProbScalarRD", "firth_rd_derivatives")
  saved <- mget(poisoned, envir = .GlobalEnv)
  on.exit(for (key in poisoned) assign(key, saved[[key]], envir = .GlobalEnv))
  for (key in poisoned) assign(key, function(...) stop("Legacy global RD function reached"), envir = .GlobalEnv)
  design <- webster_rd_scenario_design("common", .5, "diffuse", "high")
  for (weighting in c("unit", "maic")) {
    fit <- run.anchored.rd(design, 500, seed = 19, model_weighting = weighting)
    stopifnot(all(fit$estimates$status == "available"),
      all(vapply(fit$fits[c("brm", "brm-FC")], function(z)
        identical(z$probability_mapping, "stable_rd_op_v1"), logical(1))))
  }
})
check(identical(original, fingerprint()), "actual BRM/FC point, fixed/joint variance and marginal calculations bypass legacy globals")
cat("PASS: RD STABLE SOLVER CHECKS COMPLETE\n")
out <- Sys.getenv("RD_TEST_OUTPUT")
if (nzchar(out)) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("PASS: RD stable solver checks complete", "PASS: legacy globals and RR outputs unchanged",
    "PASS: extreme probabilities, input interfaces, independent oracle and FC derivatives",
    "PASS: actual BRM/FC point and variance bypass legacy globals"),
    file.path(out, "check_rd_stable_solver.txt"))
}
