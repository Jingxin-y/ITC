args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript RD/summary/summarize_rd.R INPUT_RDS [OUTPUT_DIR] [SHARED_CODE_DIR]")
own <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
rd_dir <- dirname(dirname(normalizePath(own, winslash = "/")))
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source(file.path(rd_dir, "load_anchored_rd.R"), encoding = "UTF-8")
source_anchored_rd(if (length(args) >= 3) args[3] else dirname(rd_dir), rd_dir)
input <- if (dir.exists(args[1])) file.path(args[1], "anchored_rd.rds") else args[1]
rd_write_simulation(readRDS(input), if (length(args) >= 2) args[2] else dirname(input))
