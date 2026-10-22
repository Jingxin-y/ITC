## Run from anchored_v1: Rscript data/tests/check_webster_rr.R
source("load_anchored.R", encoding="UTF-8")
source_anchored_data()
close_to <- function(a,b,tol=1e-11)
  stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
bad <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))

d <- webster_rr_design()
par <- d$parameters
z <- as.matrix(expand.grid(rep(list(0:1),6)))
E <- (1-z[,1])+z[,2]+(1-z[,3])
P <- z[,4]+(1-z[,5])+z[,6]
index <- 1+E+4*P
tab <- d$risk_table[index,]
close_to(tab$p_C,par$baseline_risk*par$baseline_rr_E^E*par$baseline_rr_P^P)
close_to(tab$p_A,tab$p_C*par$rr_AC*par$em_rr_AC^E)
close_to(tab$p_B,tab$p_C*par$rr_BC*par$em_rr_BC^E)
stopifnot(!"log_GOP" %in% names(tab),
          identical(d$covariate_roles$effect_modifiers,paste0("z",1:3)),
          identical(d$covariate_roles$non_effect_modifiers,paste0("z",4:6)))
close_to(tab$log_RR_AC,log(par$rr_AC)+log(par$em_rr_AC)*E)
## Flip each raw variable separately: exactly three alter the conditional RR.
for (j in 1:6) {
  flipped <- z; flipped[,j] <- 1-flipped[,j]
  Ef <- (1-flipped[,1])+flipped[,2]+(1-flipped[,3])
  Pf <- flipped[,4]+(1-flipped[,5])+flipped[,6]
  tf <- d$risk_table[1+Ef+4*Pf,]
  stopifnot(all(abs(tf$p_C-tab$p_C)>1e-6))
  if (j<=3) stopifnot(all(abs(tf$log_RR_AC-tab$log_RR_AC)>1e-6)) else
    close_to(tf$log_RR_AC,tab$log_RR_AC)
}
cat("PASS: direct log-risk formulas and exactly three EMs / three prognostic non-EMs\n")

## Independent 64-cell integration AND analytic Bernoulli product moments.
for (population in c("AC","BC")) {
  prob <- if (population=="AC") d$q else d$e
  mass <- apply(sweep(z,2,prob,function(x,p) p^x*(1-p)^(1-x)),1,prod)
  risk <- colSums(tab[c("p_A","p_B","p_C")]*mass)
  close_to(risk,d$marginal_truth[population,names(risk)])
  u <- c(1-prob[1],prob[2],1-prob[3])
  v <- c(prob[4],1-prob[5],prob[6])
  control <- par$baseline_risk*prod(1-u+u*par$baseline_rr_E)*prod(1-v+v*par$baseline_rr_P)
  active <- par$baseline_risk*par$rr_AC*prod(1-u+u*par$baseline_rr_E*par$em_rr_AC)*
    prod(1-v+v*par$baseline_rr_P)
  active_B <- par$baseline_risk*par$rr_BC*prod(1-u+u*par$baseline_rr_E*par$em_rr_BC)*
    prod(1-v+v*par$baseline_rr_P)
  close_to(c(active,active_B,control),risk)
  close_to(d$marginal_truth[population,"log_RR_AC"],log(active/control))
}
stopifnot(all(d$marginal_truth[,"log_RR_AB"]>0))
alt <- webster_rr_design(rr_BC=.8,em_rr_BC=.9)
close_to(alt$risk_table$p_C,d$risk_table$p_C)
close_to(alt$risk_table$p_A,d$risk_table$p_A)
stopifnot(all(abs(alt$marginal_truth[,"log_RR_AB"])>.01))
close_to(alt$risk_table$p_B/alt$risk_table$p_C,.8*.9^alt$risk_table$S_E)
rare <- webster_rr_design(event="Rare")
for (current in list(d,rare)) {
  old <- joint_rr_gop_design(event=current$event,ess_ratio=.5)
  close_to(current$marginal_truth,old$marginal_truth)
  stopifnot(current$calibration$matches_reference_targets,
            all(as.matrix(current$risk_table[c("p_A","p_B","p_C")])>0),
            all(as.matrix(current$risk_table[c("p_A","p_B","p_C")])<1))
  ## Change only independent prognostic-variable distribution: RR stays fixed.
  changed <- current$q; changed[4:6] <- c(.1,.2,.9)
  other <- rr_gop_marginal_truth(current$risk_table,changed)
  close_to(other[c("log_RR_AC","log_RR_BC","log_RR_AB")],
           current$marginal_truth["AC",c("log_RR_AC","log_RR_BC","log_RR_AB")])
  stopifnot(abs(other["p_C"]-current$marginal_truth["AC","p_C"])>1e-4)
}
## The previously used .67 recipe remains available through explicit inputs.
paper <- webster_rr_design(baseline_risk=.9,rr_AC=.67,rr_BC=.67,
  em_rr_AC=.67,em_rr_BC=.67,baseline_rr_E=.67,baseline_rr_P=.67)
close_to(paper$risk_table$p_A,.9*.67^(1+2*paper$risk_table$S_E+paper$risk_table$S_P))
stopifnot(!paper$calibration$matches_reference_targets)
cat("PASS: both old populations' common/rare risks recovered; non-EM marginal RR invariant\n")
cat("PASS: exact marginal risks, analytic product moments, B/C alternatives and rare scenario\n")

## Preserve the original population/overlap construction.
old <- joint_rr_gop_design(event="common",ess_ratio=.5)
close_to(d$e,old$e); close_to(d$q,old$q)
for (ess in c(1,.8,.5,.25,.1)) {
  current <- webster_rr_design(ess_ratio=ess)
  close_to(current$ess_ratio,ess)
  stopifnot(identical(current$risk_table,d$risk_table))
}
q <- c(z6=.2,z5=.8,z4=.4,z3=.7,z2=.5,z1=.6)
close_to(webster_rr_design(q=q)$q,q[paste0("z",1:6)])
for (args in list(list(rr_AC=0),list(em_rr_BC=NA_real_),list(baseline_risk=1.1),
                 list(baseline_risk=.9,rr_AC=2),list(ess_ratio=0),
                 list(direction=rep(0,6)),list(q=rep(1,6)),list(event="bad")))
  bad(do.call(webster_rr_design,args))
cat("PASS: original ESS/covariate design and rejection of invalid probabilities\n")

set.seed(591)
sim <- generate_anchored_data(40000,40000,d)
for (name in c("data.IPD","data.AD.or")) {
  trial <- sim[[name]]
  dat <- trial$data
  close_to(trial$count[c("Na0","Na1")],c(20000,20000))
  stopifnot(max(abs(colMeans(dat[paste0("z",1:6)])-trial$design$covariate_prob))<.015)
  E <- (1-dat$z1)+dat$z2+(1-dat$z3)
  P <- dat$z4+(1-dat$z5)+dat$z6
  pc <- par$baseline_risk*par$baseline_rr_E^E*par$baseline_rr_P^P
  pa <- if (name=="data.IPD") pc*par$rr_AC*par$em_rr_AC^E else
    pc*par$rr_BC*par$em_rr_BC^E
  close_to(trial$truth$p0,pc); close_to(trial$truth$p1,pa)
  close_to(trial$truth$observed_risk,ifelse(dat$x==1,pa,pc))
  for (arm in 0:1) {
    keep <- dat$x==arm
    p <- trial$truth$observed_risk[keep]
    stopifnot(abs(sum(dat$y[keep]-p))<6*sqrt(sum(p*(1-p))))
  }
}
close_to(sim$data.AD$anchored.mean.x,colMeans(sim$data.AD.or$data[paste0("z",1:6)]))
set.seed(25); a <- generate_anchored_data(100,150,d)
set.seed(25); b <- generate_anchored_data(100,150,d)
stopifnot(identical(a,b))
cat("PASS: sampling, balanced randomization, observed outcome risks, AD summaries and reproducibility\n")
print(d$marginal_truth)
