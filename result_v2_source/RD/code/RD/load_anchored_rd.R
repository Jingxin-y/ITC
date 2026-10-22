## All additions live in RD; existing RR and shared modules are read only.
.anchored_rd_module_dir <- local({
  paths <- vapply(sys.frames(), function(frame)
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1], character(1))
  path <- tail(paths[nzchar(paths)], 1)
  if (!length(path)) stop("Source load_anchored_rd.R from a saved file.")
  dirname(normalizePath(path, winslash = "/", mustWork = TRUE))
})
source_anchored_rd <- function(code_dir = dirname(.anchored_rd_module_dir),
                               rd_dir = .anchored_rd_module_dir) {
  code_dir <- normalizePath(code_dir, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(code_dir, "load_anchored.R")))
    stop("code_dir must be the original anchored_v1 code directory.")
  source(file.path(code_dir, "load_anchored.R"), encoding = "UTF-8")
  source_anchored(code_dir)
  for (f in c("data/generate_webster_rd_scenarios.R", "estimation/rd_methods.R", "estimation/rd_solver_adapter.R",
              "estimation/CI_marginal_rd_target.R", "estimation/rd_comparison.R",
              "estimation/rd_simulation.R", "estimation/rd_cluster.R"))
    source(file.path(rd_dir, f), encoding = "UTF-8")
  invisible(list(code_dir = code_dir, rd_dir = rd_dir))
}
