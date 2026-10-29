#!/usr/bin/env Rscript
# Read saved RDS fits; write only to a separate STC output tree.
stc_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("^--file=", "", arg[1]), winslash = "/", mustWork = TRUE)))
  frames <- sys.frames()
  paths <- vapply(frames, function(x) if (is.null(x$ofile)) "" else x$ofile, character(1))
  if (any(nzchar(paths))) return(dirname(normalizePath(tail(paths[nzchar(paths)], 1), winslash = "/", mustWork = TRUE)))
  getwd()
}

stc_canonical_path <- function(path) {
  path <- path.expand(path)
  if (file.exists(path) || dir.exists(path))
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  parent <- dirname(path)
  if (identical(parent, path)) stop("Cannot resolve path: ", path)
  file.path(stc_canonical_path(parent), basename(path))
}
stc_path_contains <- function(parent, child) {
  parent <- sub("/+$", "", gsub("\\\\", "/", parent))
  child <- sub("/+$", "", gsub("\\\\", "/", child))
  if (.Platform$OS.type == "windows") { parent <- tolower(parent); child <- tolower(child) }
  identical(parent, child) || startsWith(child, paste0(parent, "/"))
}
stc_check_output <- function(input_root, output_root) {
  input <- stc_canonical_path(input_root); output <- stc_canonical_path(output_root)
  if (stc_path_contains(input, output) || stc_path_contains(output, input))
    stop("Input and output trees must be separate, with neither containing the other.", call. = FALSE)
  list(input = input, output = output)
}
stc_atomic <- function(object, path, output_root, csv = FALSE) {
  # Resolve existing symlink/junction ancestors before every write.
  path <- stc_canonical_path(path)
  if (!stc_path_contains(output_root, path) || identical(output_root, path))
    stop("Refusing write outside the validated STC output tree: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".stc_write_", tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  if (csv) utils::write.csv(object, tmp, row.names = FALSE, na = "") else saveRDS(object, tmp)
  if (.Platform$OS.type == "windows" && file.exists(path)) {
    if (unlink(path) != 0L) stop("Cannot replace owned STC output: ", path)
  }
  if (!file.rename(tmp, path)) stop("Cannot finalize STC output: ", path)
}
stc_parse_ids <- function(x, max_id) {
  bits <- strsplit(x, ",", fixed = TRUE)[[1]]
  ids <- unlist(lapply(bits, function(s) {
    if (!grepl("^[0-9]+(:[0-9]+)?$", s)) stop("Invalid integer selection: ", x)
    a <- as.integer(strsplit(s, ":", fixed = TRUE)[[1]])
    if (anyNA(a) || any(a < 1 | a > max_id) || (length(a) == 2 && a[1] > a[2]))
      stop("Selection out of range: ", x)
    if (length(a) == 1L) a else seq.int(a[1], a[2])
  }))
  sort(unique(ids))
}
stc_parse_args <- function(args) {
  opts <- list(input_root = "/scratch/yanjin41/RR_v2_full24_MC1000_results_20261024",
    output_root = "/scratch/yanjin41/STC_RR_reuse_results_20261026",
    scenarios = 1:8, seeds = 1:1000, level = .95)
  if ("--help" %in% args) {
    cat("Rscript run_stc_reuse.R --input-root DIR --output-root NEW_DIR\n",
        "  [--scenarios 1:8] [--seeds 1:1000] [--level 0.95]\n",
        "Reads only score_unit saved RDS fits. Never fits or simulates.\n", sep = "")
    return(NULL)
  }
  if (length(args) %% 2L) stop("Arguments must be --name value pairs; see --help.")
  seen <- character(0)
  if (length(args)) for (i in seq.int(1L, length(args), by = 2L)) {
    key <- args[i]; value <- args[i + 1L]
    if (!key %in% c("--input-root", "--output-root", "--scenarios", "--seeds", "--level") || key %in% seen)
      stop("Unknown or repeated option: ", key)
    seen <- c(seen, key)
    name <- gsub("-", "_", sub("^--", "", key), fixed = TRUE)
    opts[[name]] <- switch(name,
      scenarios = stc_parse_ids(value, 8L), seeds = stc_parse_ids(value, 1000L),
      level = suppressWarnings(as.numeric(value)), value)
  }
  if (length(opts$level) != 1L || !is.finite(opts$level) || opts$level <= 0 || opts$level >= 1)
    stop("--level must be in (0,1).")
  opts
}

stc_process_file <- function(path, expected, level, runtime_hash, checkpoint, output_root) {
  start <- proc.time()[["elapsed"]]
  before <- if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_
  if (file.exists(checkpoint)) {
    saved <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
    if (is.list(saved) && identical(saved$runtime_hash, runtime_hash) &&
        identical(saved$source_md5, before) && isTRUE(saved$audit$source_unchanged) &&
        identical(saved$level, level)) {
      cache_after <- if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_
      if (!is.na(before) && identical(before, cache_after)) {
        saved$audit$source_md5_before <- before
        saved$audit$source_md5_after <- cache_after
        saved$audit$cached <- TRUE
        saved$audit$elapsed_seconds <- proc.time()[["elapsed"]] - start
        return(saved)
      }
    }
  }
  audit <- data.frame(scenario_id = expected$scenario_id, cell = expected$cell,
    seed = expected$seed, source_file = path, source_md5_before = before,
    source_md5_after = NA_character_, source_unchanged = NA, record_complete = NA,
    input_status = "pending", issue = "", source_settings_version = NA_integer_,
    source_code_files = NA_integer_, models_with_points = 0L,
    intervals_returned = 0L, cached = FALSE, elapsed_seconds = NA_real_,
    version = stc_version, stringsAsFactors = FALSE)
  record <- NULL
  answer <- tryCatch({
    stc_require(file.exists(path), "source_file_missing", "Original seed RDS is missing.")
    stc_require(is.character(before) && !is.na(before), "source_hash_failed", "Cannot checksum source RDS.")
    record <- tryCatch(readRDS(path), error = function(e) stc_stop("source_read_failed", conditionMessage(e)))
    audit$record_complete <- isTRUE(record$complete)
    v <- record$settings$version
    if (length(v) == 1L) audit$source_settings_version <- as.integer(v)
    audit$source_code_files <- length(record$source_md5)
    out <- stc_from_record(record, expected = expected, level = level, source_file = path)
    audit$input_status <- if (isTRUE(record$complete)) "available" else "source_incomplete"
    out
  }, error = identity)
  after <- if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_
  audit$source_md5_after <- after
  audit$source_unchanged <- !is.na(before) && identical(before, after)
  if (!is.na(before) && !isTRUE(audit$source_unchanged))
    answer <- structure(list(message = "Source RDS changed while being read; results discarded.",
      status = "source_changed"), class = c("stc_error", "error", "condition"))
  if (inherits(answer, "error")) {
    audit$input_status <- stc_error_status(answer)
    audit$issue <- conditionMessage(answer)
    rows <- stc_empty_rows(expected, expected$seed, audit$input_status, audit$issue, level)
    rows$source_file <- path
    # An incompatible/unreadable input is not used even as the MAIC baseline.
    answer <- list(summary = rows, paired = stc_pair_rows(rows), diagnostics = data.frame())
  }
  answer$summary$source_md5 <- before; answer$paired$source_md5 <- before
  audit$models_with_points <- length(unique(answer$summary$method[answer$summary$point_status == "available"]))
  audit$intervals_returned <- sum(answer$summary$status == "available")
  audit$elapsed_seconds <- proc.time()[["elapsed"]] - start
  saved <- c(answer, list(audit = audit, source_md5 = before, runtime_hash = runtime_hash,
                          level = level, version = stc_version))
  stc_atomic(saved, checkpoint, output_root)
  saved
}

stc_run <- function(opts, code_dir) {
  if (!dir.exists(opts$input_root)) stop("Input root does not exist: ", opts$input_root)
  paths <- stc_check_output(opts$input_root, opts$output_root)
  input <- paths$input; output <- paths$output
  runtime <- c("stc_reuse.R", "run_stc_reuse.R")
  runtime_hash <- setNames(unname(tools::md5sum(file.path(code_dir, runtime))), runtime)
  if (anyNA(runtime_hash)) stop("Cannot checksum STC runtime files.")
  cfg <- list(version = stc_version, input_root = input, scenarios = opts$scenarios,
    seeds = opts$seeds, level = opts$level, runtime_hash = runtime_hash,
    planned_per_scenario = length(opts$seeds))
  if (dir.exists(output)) {
    config_path <- file.path(output, "configuration.rds")
    if (!file.exists(config_path) || !identical(readRDS(config_path), cfg))
      stop("Existing output has no matching STC configuration; choose a new output directory.")
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  stc_check_output(input, output)
  lock <- file.path(output, ".stc_reuse_lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Output is locked by another run: ", lock)
  # Only this newly created, empty lock directory is removed.
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  stc_atomic(cfg, file.path(output, "configuration.rds"), output)
  scenarios <- stc_scenarios(); scenarios <- scenarios[scenarios$scenario_id %in% opts$scenarios, , drop = FALSE]
  combined <- paired <- audits <- list(); count <- 0L
  for (i in seq_len(nrow(scenarios))) {
    r <- scenarios[i, , drop = FALSE]
    for (seed in opts$seeds) {
      expected <- r; expected$seed <- seed
      source <- file.path(input, "cells", r$cell, "by_seed", sprintf("seed_%06d.rds", seed))
      source <- stc_canonical_path(source)
      if (!stc_path_contains(input, source)) stop("Source symlink escapes the input root: ", source)
      checkpoint <- file.path(output, "cells", r$cell, "by_seed", sprintf("stc_seed_%06d.rds", seed))
      z <- stc_process_file(source, expected, opts$level, runtime_hash, checkpoint, output)
      count <- count + 1L
      combined[[count]] <- z$summary; paired[[count]] <- z$paired; audits[[count]] <- z$audit
      if (seed == opts$seeds[1] || seed == tail(opts$seeds, 1) || count %% 25L == 0L) {
        cat(sprintf("STC REUSE %s seed=%d input=%s points=%d/12 intervals=%d/24 cached=%s\n",
          r$cell, seed, z$audit$input_status, z$audit$models_with_points, z$audit$intervals_returned, z$audit$cached))
        flush.console()
      }
    }
  }
  rows <- do.call(rbind, combined); pair <- do.call(rbind, paired); audit <- do.call(rbind, audits)
  summaries <- stc_summarize(pair, length(opts$seeds))
  report <- file.path(output, "report")
  tables <- c(list(combined_results = rows, paired_results = pair, input_audit = audit), summaries)
  issues <- rows[rows$point_status != "available" | rows$se_status != "available",
    c("scenario_id", "cell", "seed", "method", "uncertainty", "status", "message", "source_file")]
  tables$missing_or_failed <- issues
  tables$status_counts <- as.data.frame(table(scenario_id = rows$scenario_id, method = rows$method,
    uncertainty = rows$uncertainty, status = rows$status), stringsAsFactors = FALSE)
  tables$status_counts <- tables$status_counts[tables$status_counts$Freq > 0, ]
  for (name in names(tables)) stc_atomic(tables[[name]], file.path(report, paste0(name, ".csv")), output, csv = TRUE)
  stc_atomic(list(configuration = cfg, source_files = count,
    available_sources = sum(audit$input_status == "available"), rows = nrow(rows),
    returned_points = sum(rows$point_status == "available"), returned_intervals = sum(rows$status == "available"),
    source_hashes_unchanged = all(audit$source_unchanged[!is.na(audit$source_md5_before)]),
    session_info = sessionInfo()), file.path(report, "run_metadata.rds"), output)
  cat(sprintf("FINISHED: %d source slots, %d STC rows. Report: %s\n", count, nrow(rows), report))
  bad <- any(!audit$input_status %in% "available")
  if (bad) cat("INPUT ISSUES: see input_audit.csv and missing_or_failed.csv; no simulation or refit was attempted.\n")
  invisible(list(exit_status = if (bad) 2L else 0L, audit = audit, summary = summaries$method_summary))
}

stc_main <- function() {
  code_dir <- stc_script_dir()
  source(file.path(code_dir, "stc_reuse.R"), encoding = "UTF-8", local = .GlobalEnv)
  opts <- stc_parse_args(commandArgs(TRUE))
  if (is.null(opts)) return(invisible(NULL))
  result <- stc_run(opts, code_dir)
  if (result$exit_status != 0L) quit(status = result$exit_status)
  invisible(result)
}
if (sys.nframe() == 0L) stc_main()
