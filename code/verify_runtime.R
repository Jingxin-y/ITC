args <- commandArgs(trailingOnly=TRUE)
scale <- match.arg(if(length(args)) args[1L] else "RR", c("RR","RD"))
if(.Platform$OS.type=="windows") Sys.setlocale("LC_CTYPE","English_United States.utf8")
root <- normalizePath(".", winslash="/", mustWork=TRUE)
code <- file.path(root,scale,"code")
files <- list.files(code, "\\.[Rr]$", recursive=TRUE, full.names=TRUE)
for(p in files) parse(p, encoding="UTF-8")
if(scale=="RR") {
  source(file.path(code,"load_anchored.R"), encoding="UTF-8")
  source_anchored(code)
  for(event in c("common","rare")) {
    d <- webster_rr_scenario_design(event,.5,"diffuse","high")
    stopifnot(identical(d$dgm,"webster2024_rr_ess_em_v2"),d$scenario$version==2L)
    expected <- if(event=="common") c(.418596226792393,.514669087904676,.60) else
      c(.0488309131792698,.10,.15)
    stopifnot(max(abs(d$marginal_truth["BC",c("p_C","p_A","p_B")]-expected))<1e-10)
  }
  stopifnot(exists(".firth_rr_rescue"),formals(exact_marginal_target)$bisection.max.step==15L)
} else {
  source(file.path(code,"RD/load_anchored_rd.R"), encoding="UTF-8")
  source_anchored_rd(code,file.path(code,"RD"))
  for(event in c("common","rare")) {
    d <- webster_rd_scenario_design(event,.5,"diffuse","high")
    stopifnot(identical(d$dgm,"webster2024_rd_ess_em_v2"),
      abs(d$marginal_truth["BC","p_A"]-d$marginal_truth["BC","p_B"]-
        if(event=="common") .14 else .05)<1e-10)
  }
  stopifnot(formals(exact_marginal_rd_target)$bisection.max.step==15L,
            formals(exact_marginal_rd_target)$thres.dicho==.001)
}
stopifnot(requireNamespace("MASS",quietly=TRUE))
cat("PASS",scale,":",length(files),"R source files parsed; loaders, v2 DGM and published solver defaults verified.\n")
cat("nleqslv available:",requireNamespace("nleqslv",quietly=TRUE),
    "(required for RR FC rescue; that branch is not exercised by this loading check).\n")
sessionInfo()
