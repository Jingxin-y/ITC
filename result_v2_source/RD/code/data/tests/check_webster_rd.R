## Run from anchored_v1: Rscript data/tests/check_webster_rd.R
source("load_anchored.R",encoding="UTF-8")
source_anchored_data()
close_to <- function(a,b,tol=1e-11)
  stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
bad <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
z <- as.matrix(expand.grid(rep(list(0:1),6)))
E <- (1-z[,1])+z[,2]+(1-z[,3])
P <- z[,4]+(1-z[,5])+z[,6]
for (event in c("common","rare")) {
  d <- webster_rd_design(event=event)
  par <- d$parameters
  tab <- d$risk_table[1+E+4*P,]
  old <- joint_rr_gop_design(event=event,ess_ratio=.5)
  close_to(d$marginal_truth[,names(old$marginal_truth[1,])],old$marginal_truth)
  stopifnot(d$calibration$matches_reference_targets,identical(d$effect_scale,"RD"),
            identical(d$covariate_roles$effect_modifiers,paste0("z",1:3)),
            identical(d$covariate_roles$non_effect_modifiers,paste0("z",4:6)),
            !"coefficient_truth" %in% names(d))
  for (arm in c("A","B")) {
    rd <- tab[[paste0("p_",arm)]]-tab$p_C
    close_to(rd,par[[paste0("rd_",arm,"C")]]+par[[paste0("em_rd_",arm,"C")]]*E)
    for (j in 1:6) {
      flip <- z; flip[,j] <- 1-flip[,j]
      EF <- (1-flip[,1])+flip[,2]+(1-flip[,3])
      PF <- flip[,4]+(1-flip[,5])+flip[,6]
      tf <- d$risk_table[1+EF+4*PF,]
      rd_flip <- tf[[paste0("p_",arm)]]-tf$p_C
      stopifnot(all(abs(tf$p_C-tab$p_C)>1e-6))
      if (j<=3) stopifnot(all(abs(rd_flip-rd)>1e-6)) else close_to(rd_flip,rd)
    }
  }
  for (pop in c("AC","BC")) {
    prob <- if (pop=="AC") d$q else d$e
    mass <- apply(sweep(z,2,prob,function(x,p)p^x*(1-p)^(1-x)),1,prod)
    risks <- colSums(tab[c("p_A","p_B","p_C")]*mass)
    close_to(risks,d$marginal_truth[pop,names(risks)])
    mean_E <- sum(mass*E)
    close_to(d$marginal_truth[pop,c("RD_AC","RD_BC")],
      c(par$rd_AC+par$em_rd_AC*mean_E,par$rd_BC+par$em_rd_BC*mean_E))
  }
  ## Change SP while holding the complete SE distribution fixed, allowing
  ## SE/SP dependence: marginal RD is unchanged without independence.
  base <- rep(1/64,64)
  correlated <- base*(1+.3*P+.1*E*P)
  correlated <- correlated / ave(correlated,E,FUN=sum) * ave(base,E,FUN=sum)
  close_to(tapply(base,E,sum),tapply(correlated,E,sum))
  for (arm in c("A","B")) {
    rd <- tab[[paste0("p_",arm)]]-tab$p_C
    close_to(sum(base*rd),sum(correlated*rd))
  }
  stopifnot(abs(sum(base*tab$p_C)-sum(correlated*tab$p_C))>1e-5)
  ## ESS affects populations; it must not silently change the risk formula.
  for (ess in c(1,.8,.25,.1)) {
    other <- webster_rd_design(event=event,ess_ratio=ess)
    stopifnot(identical(d$risk_table,other$risk_table))
    close_to(other$ess_ratio,ess)
  }
  set.seed(401)
  sim <- generate_anchored_data(4000,5000,d)
  for (name in c("data.IPD","data.AD.or")) {
    trial <- sim[[name]]
    tt <- d$risk_table[trial$truth$risk_table_row,]
    close_to(trial$truth$p0,tt$p_C)
    close_to(trial$truth$p1,tt[[if(name=="data.IPD")"p_A" else "p_B"]])
    p <- trial$truth$observed_risk
    stopifnot(abs(sum(trial$data$y-p))<6*sqrt(sum(p*(1-p))))
  }
  cat("PASS:",event,"RD calibration, three EMs / three non-EMs, correlated marginal invariance, sampling\n")
}
default <- webster_rd_design()
alternative <- webster_rd_design(rd_BC=.01,em_rd_BC=.05)
close_to(default$risk_table$p_C,alternative$risk_table$p_C)
close_to(default$risk_table$p_A,alternative$risk_table$p_A)
stopifnot(!alternative$calibration$matches_reference_targets,
          identical(alternative$calibration$overridden_parameters,c("rd_BC","em_rd_BC")))
for (args in list(list(baseline_risk=-.1),list(rd_AC=2),list(beta_P=1),
                 list(em_rd_BC=NA_real_),list(beta_E=c(.1,.2)),
                 list(ess_ratio=0),list(q=rep(1,6)),list(event="bad")))
  bad(do.call(webster_rd_design,args))
cat("PASS: explicit overrides preserve other arms; invalid probabilities rejected without clipping\n")
